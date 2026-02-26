#include "VzhuhSolver.h"

#include <algorithm>
#include <queue>

VzhuhSolver::VzhuhSolver(const std::vector<std::tuple<int, int, int> > &edge_list_,
                         const SolverParameters &params_) : primal_objective(INT64_MAX), dual_objective(INT64_MIN),
                                                            params(params_),
                                                            num_vertices_elementary(InitNumVertices(edge_list_)),
                                                            current_round(0) {
    // TODO check validity of edge_list_

    nodes.reserve(2 * num_vertices_elementary);
    for (int i = 0; i < num_vertices_elementary; ++i) {
        nodes.emplace_back(std::make_unique<Node>(i));
    }

    // reserve for adj list
    std::vector<int> degrees(num_vertices_elementary, 0);
    for (auto [from, to, weight] : edge_list_) {
        ++degrees[from];
        ++degrees[to];
    }

    adj_list = std::vector<std::vector<int> >(num_vertices_elementary, std::vector<int>());
    for (int i = 0; i < num_vertices_elementary; ++i) {
        adj_list[i].reserve(degrees[i]);
    }

    edges.reserve(edge_list_.size());
    for (int i = 0; i < static_cast<int>(edge_list_.size()); ++i) {
        int head = std::get<0>(edge_list_[i]);
        int tail = std::get<1>(edge_list_[i]);
        edges.emplace_back(std::make_unique<Edge>(head, tail, std::get<2>(edge_list_[i]), i));
        adj_list[head].emplace_back(i);
        adj_list[tail].emplace_back(i);
    }

    nodes_label_cnt = 1;
    num_trees_alive = INT32_MAX;

    aux_counter1 = 0;
    aux_counter2 = 0;
    aux_counter3 = 0;
    aux_counter4 = 0;
}

void VzhuhSolver::FindMinPerfectMatching() {
    GreedyInit();
    InitializeTrees();

    if (params.verbose) {
        PrintGraph();
    }

    while (!alive_trees.empty()) {
        if (params.print_statistics) {
            std::cout << "round " << current_round << std::endl;
            std::cout << "trees left: " << alive_trees.size() << std::endl;
            std::cout << "sample tree sizes: " << trees[alive_trees[0]]->tree_nodes.size() << " " << trees[alive_trees[
                    alive_trees.size() / 2]]->tree_nodes.size() << " " << trees[alive_trees.back()]->tree_nodes.size()
                <<
                std::endl;

            int min_children = INT32_MAX;
            for (std::unique_ptr<Node> &node : nodes) {
                if (!node->blossom_children.empty() && node->blossom_children.size() < min_children) {
                    min_children = node->blossom_children.size();
                }
            }
            std::cout << "min num children: " << min_children << std::endl;
        }

        ++current_round;
        bool primal_action = MakePrimalUpdates();
        bool dual_action = MakeDualUpdates();

        if (!primal_action && !dual_action) {
            throw std::runtime_error("No progress made in either of primal and dual updates");
        }
    }

    // compute the objectives and recover the matching
    ComputeDualObjectiveQuadrupled();
    if (dual_objective % 4 != 0) {
        throw std::runtime_error("Dual objective not integer");
    }
    dual_objective /= 4;

    if (params.verbose) {
        std::cout << "Dual objective:\t\t" << dual_objective << std::endl;
    }

    if (params.compute_dual_certificate) {
        ComputeDualCertificate();
    }

    DestroyBlossoms();
    ComputeMatching();
    ComputePrimalObjective();

    // std::cout << aux_counter1 << " " << aux_counter2 << " " << aux_counter3 << " " << aux_counter4 << std::endl;

    if (params.verbose) {
        std::cout << "Primal objective:\t" << primal_objective << std::endl;
    }
}

const std::vector<std::pair<int, int> > &VzhuhSolver::Matching() const {
    return matching;
}

const std::vector<std::tuple<int, int, int> > &VzhuhSolver::DualCertificate() const {
    if (!params.compute_dual_certificate) {
        throw std::runtime_error("Parameter \"compute_dual_certificate\" is set to false");
    }
    return dual_certificate;
}

VzhuhSolver::Edge::Edge(int head_, int tail_, int weight_, int index_) : heap(nullptr), heap_child(nullptr),
                                                                         heap_next(nullptr),
                                                                         heap_prev(nullptr), index(index_),
                                                                         head(head_),
                                                                         tail(tail_),
                                                                         elementary_head(head_),
                                                                         elementary_tail(tail_) {
    queue_index = -1;
    weight = weight_;
    slack_quadrupled_amortized_ = 4 * weight_;
    matched = false;
    maybe_has_zero_slack = false;
    must_be_updated = false;
    maybe_was_loop = true;
    slack_diff = 0;
}

int VzhuhSolver::Edge::Key() const {
    return slack_quadrupled_amortized_;
}

VzhuhSolver::Node::Node(int index_) : heap(nullptr), heap_child(nullptr), heap_next(nullptr), heap_prev(nullptr),
                                      blossom_parent(-1),
                                      old_blossom_parent(-1),
                                      matched_edge(-1),
                                      minus_parent(-1), receptacle_(index_), tree(-1), old_tree(-1),
                                      tree_var_at_birth(0), round_0slack_neighbors_updated(-1), is_in_record(false),
                                      index(index_) {
    is_alive = true;
    dual_var_quadrupled_amortized_ = 0;
    plus = false;
    old_plus = false;
    queue_index = -1;
    label = 0;
}

int VzhuhSolver::Node::Key() const {
    return dual_var_quadrupled_amortized_;
}

VzhuhSolver::Tree::Tree(int root_,
                        int minus_blossoms_,
                        int plus_empty_edges_,
                        int plus_plus_internal_edges_) : root(root_), minus_blossoms(minus_blossoms_),
                                                         plus_empty_edges(plus_empty_edges_),
                                                         plus_plus_internal_edges(plus_plus_internal_edges_) {
    is_alive = true;
    dual_var_quadrupled = 0;
    alive_index = 0;
    tree_nodes = {root};
}

void VzhuhSolver::PrintGraph() const {
    std::cout << "Adjacency list (to, weight, matched):" << std::endl;

    for (int i = 0; i < num_vertices_elementary; ++i) {
        std::cout << i << ": ";
        for (int edge : adj_list[i]) {
            std::cout << "(" << OtherElementaryEnd(edge, i) << " " << edges[edge]->weight << " "
                << edges[edge]->matched << ") ";
        }
        std::cout << std::endl;
    }

    std::cout << "Node structure: " << std::endl;
    for (int i(0); i < nodes.size(); ++i) {
        if (nodes[i]->is_alive) {
            PrintNode(i);
        }
    }

    std::cout << "Tree structure: " << std::endl;
    for (int tree : alive_trees) {
        if (!trees[tree]->is_alive) {
            continue;
        }
        std::cout << "root: " << trees[tree]->root << " var: " << trees[tree]->dual_var_quadrupled / 4.
            << std::endl;
    }
}

void VzhuhSolver::PrintNode(int node) const {
    std::cout << node << ": y_v = " << DualVariableQuadrupled(node) / 4.;
    if (nodes[node]->blossom_parent >= 0) {
        std::cout << " blossom_parent: " << nodes[node]->blossom_parent;
    }
    std::cout << std::endl;
}

void VzhuhSolver::GreedyInit() {
    // first, make all the slacks non-negative
    for (int i = 0; i < num_vertices_elementary; ++i) {
        if (adj_list[i].empty()) {
            std::cout << "vertex " << i << ": no neighbors" << std::endl;
            throw std::runtime_error("Found an isolated vertex => no perfect matching exists");
        }

        int min_weight = edges[adj_list[i].front()]->weight;
        for (int edge : adj_list[i]) {
            if (edges[edge]->weight < min_weight) {
                min_weight = edges[edge]->weight;
            }
        }

        nodes[i]->dual_var_quadrupled_amortized_ += 2 * min_weight;
        for (int edge : adj_list[i]) {
            edges[edge]->slack_quadrupled_amortized_ -= 2 * min_weight;
        }
    }

    for (int i = 0; i < num_vertices_elementary; ++i) {
        int node_idx(i);

        if (nodes[node_idx]->matched_edge >= 0) {
            continue;
        }

        int smallest_slack_edge = adj_list[i].front();
        for (const int edge : adj_list[i]) {
            if (edges[edge]->slack_quadrupled_amortized_ < edges[smallest_slack_edge]->slack_quadrupled_amortized_) {
                smallest_slack_edge = edge;
            }
        }

        int diff = edges[smallest_slack_edge]->slack_quadrupled_amortized_;
        nodes[node_idx]->dual_var_quadrupled_amortized_ += diff;
        for (int edge : adj_list[i]) {
            edges[edge]->slack_quadrupled_amortized_ -= diff;
        }
        if (nodes[OtherElementaryEnd(smallest_slack_edge, i)]->matched_edge < 0) {
            // if the other vertex is also unmatched, match the edge
            edges[smallest_slack_edge]->matched = true;
            nodes[node_idx]->matched_edge = smallest_slack_edge;
            nodes[OtherElementaryEnd(smallest_slack_edge, i)]->matched_edge =
                smallest_slack_edge;
        }
    }
}

void VzhuhSolver::InitializeTrees() {
    // is called after GreedyInit
    // also initializes queues and actionable_edges

    // TODO no need to populate queues in tree initialization

    std::vector<int> roots;
    for (int root_index(0); root_index < num_vertices_elementary; ++root_index) {
        if (nodes[root_index]->matched_edge >= 0) {
            continue;
        }
        roots.push_back(root_index);
    }

    num_trees_alive = static_cast<int>(roots.size());

    node_heaps.reserve(roots.size());
    edge_heaps.reserve(2 * roots.size());
    trees.reserve(roots.size());
    alive_trees.reserve(roots.size());
    for (int i = 0; i < static_cast<int>(roots.size()); ++i) {
        trees.emplace_back(std::make_unique<Tree>(roots[i], i, 2 * i, 2 * i + 1));
        alive_trees.emplace_back(i);

        node_heaps.emplace_back(std::make_unique<PairingHeap<Node> >());
        edge_heaps.emplace_back(std::make_unique<PairingHeap<Edge> >());
        edge_heaps.emplace_back(std::make_unique<PairingHeap<Edge> >());

        nodes[roots[i]]->tree = i;
        nodes[roots[i]]->plus = true;
        nodes[roots[i]]->old_tree = i;
        nodes[roots[i]]->old_plus = true;
    }

    // initialize queues and actionable_edges
    for (int root_index : roots) {
        for (int edge_to_neighbor : adj_list[root_index]) {
            int neighbor = OtherEnd(edge_to_neighbor, root_index);

            if (nodes[neighbor]->tree < 0) {
                AddEdgeToThisQueue(edge_to_neighbor, trees[nodes[root_index]->tree]->plus_empty_edges);
            } else {
                // neighbor is another root
                AddPQPlusPlus(nodes[root_index]->tree, nodes[neighbor]->tree, edge_to_neighbor);
            }

            if (edges[edge_to_neighbor]->slack_quadrupled_amortized_ == 0) {
                // safe because no dual updates has been made yet
                actionable_edges.push(edge_to_neighbor);
            }
        }
    }
}

void VzhuhSolver::ComputeMatching() {
    // must be called after the blossoms are dissolved
    if (!matching.empty()) {
        throw std::runtime_error("ComputeMatching: matching is already non-empty");
    }

    matching.reserve(num_vertices_elementary / 2);

    for (int i(0); i < static_cast<int>(edges.size()); ++i) {
        if (edges[i]->matched) {
            matching.emplace_back(edges[i]->elementary_head, edges[i]->elementary_tail);
        }
    }
}

void VzhuhSolver::ComputePrimalObjective() {
    primal_objective = 0;
    for (int i(0); i < static_cast<int>(edges.size()); ++i) {
        if (edges[i]->matched) {
            primal_objective += edges[i]->weight;
        }
    }
}

void VzhuhSolver::DestroyBlossoms() {
    RestoreFinalEdgeEnds();

    for (int blossom(nodes.size() - 1); !IsElementary(blossom); --blossom) {
        if (!nodes[blossom]->is_alive) {
            continue;
        }
        if (nodes[blossom]->blossom_parent >= 0) {
            throw std::runtime_error("In DestroyBlossoms: not going top down");
        }

        if (params.verbose) {
            std::cout << "destroying " << blossom << std::endl;
            PrintGraph();
        }

        int receptacle = FindFinalReceptacle(blossom);
        for (int child : nodes[blossom]->blossom_children) {
            nodes[child]->blossom_parent = -1;
        }
        UpdateMatching(blossom, receptacle);
        nodes[receptacle]->matched_edge = nodes[blossom]->matched_edge;
    }

    if (params.verbose) {
        std::cout << "destroyed all" << std::endl;
        PrintGraph();
    }
}

void VzhuhSolver::RestoreFinalEdgeEnds() {
    std::vector<int> depths(nodes.size(), 0);
    for (int node(nodes.size() - 1); node >= 0; --node) {
        if (!nodes[node]->is_alive) {
            continue;
        }
        if (nodes[node]->blossom_parent >= 0) {
            depths[node] = depths[nodes[node]->blossom_parent] + 1;
        }
    }

    std::vector<int> edges_to_restore;
    edges_to_restore.reserve(nodes.size() * 2);
    for (int node(0); node < nodes.size(); ++node) {
        if (nodes[node]->matched_edge >= 0) {
            edges_to_restore.push_back(nodes[node]->matched_edge);
        }
        if (nodes[node]->minus_parent >= 0) {
            edges_to_restore.push_back(nodes[node]->minus_parent);
        }
    }

    for (int edge : edges_to_restore) {
        int head = edges[edge]->elementary_head;
        int tail = edges[edge]->elementary_tail;

        while (depths[head] > 0 || depths[tail] > 0) {
            if (depths[head] > depths[tail]) {
                head = nodes[head]->blossom_parent;
            } else if (depths[tail] > depths[head]) {
                tail = nodes[tail]->blossom_parent;
            } else if (nodes[head]->blossom_parent == nodes[tail]->blossom_parent) {
                break;
            } else {
                head = nodes[head]->blossom_parent;
                tail = nodes[tail]->blossom_parent;
            }
        }

        edges[edge]->head = head;
        edges[edge]->tail = tail;
    }
}

int VzhuhSolver::FindFinalReceptacle(int blossom) const {
    int new_receptacle(-1);
    int head = edges[nodes[blossom]->matched_edge]->elementary_head;
    while (nodes[head]->blossom_parent >= 0) {
        if (nodes[head]->blossom_parent == blossom) {
            new_receptacle = head;
            break;
        }
        head = nodes[head]->blossom_parent;
    }
    if (new_receptacle < 0) {
        int tail = edges[nodes[blossom]->matched_edge]->elementary_tail;
        while (nodes[tail]->blossom_parent >= 0) {
            if (nodes[tail]->blossom_parent == blossom) {
                new_receptacle = tail;
                break;
            }
            tail = nodes[tail]->blossom_parent;
        }
    }
    if (new_receptacle < 0) {
        throw std::runtime_error("DestroyBlossoms: no new receptacle found");
    }

    return new_receptacle;
}

void VzhuhSolver::ComputeDualCertificate() {
    // can be called only before we destroy the blossoms in the end

    if (!dual_certificate.empty()) {
        throw std::runtime_error("Dual certificate is already non-empty");
    }

    std::vector<int> node_index_alive(nodes.size(), -1);
    dual_certificate.reserve(nodes.size());
    int index = 0;
    for (int i(0); i < static_cast<int>(nodes.size()); ++i) {
        if (nodes[i]->is_alive) {
            node_index_alive[i] = index;
            dual_certificate.emplace_back(index, DualVariableQuadrupled(i), -1);
            for (int child : nodes[i]->blossom_children) {
                std::get<2>(dual_certificate[node_index_alive[child]]) = index;
            }
            ++index;
        }
    }
}

void VzhuhSolver::ComputeDualObjectiveQuadrupled() {
    dual_objective = 0;
    for (int i(0); i < static_cast<int>(nodes.size()); ++i) {
        if (nodes[i]->is_alive) {
            dual_objective += DualVariableQuadrupled(i);
        }
    }
}

bool VzhuhSolver::MakePrimalUpdates() {
    bool action_taken = false;
    PrimalUpdateRecord record;
    record.changed_sign.reserve(num_vertices_elementary);

    std::vector<int> variables;
    std::vector<int> slacks;
    if (params.debug) {
        variables = NodeVariables();
        slacks = EdgeSlacks();
    }

    // expand phase
    for (int tree : alive_trees) {
        int blossom = PopExpandableBlossom(tree);
        while (blossom >= 0) {
            Expand(blossom, &record);
            blossom = PopExpandableBlossom(tree);
        }
    }

    // TODO make an early return if augmented the last pair of trees
    // unweighted solver like phase
    while (!actionable_edges.empty()) {
        int edge = actionable_edges.front();
        actionable_edges.pop();
        MakePrimalUpdate(edge, &record);
    }
    while (!actionable_nodes.empty()) {
        int node = actionable_nodes.front();
        actionable_nodes.pop();
        MakePrimalUpdateForNode(node, &record);
    }

    if (!record.changed_sign.empty()) {
        action_taken = true;
    }

    if (params.debug) {
        ValidateEvenOddPaths();
    }

    // update queues and amortized variables phase
    UpdateQueues(record);
    if (params.debug) {
        ValidateQueues();
    }

    if (params.debug) {
        std::vector<int> new_variables = NodeVariables();
        std::vector<int> new_slacks = EdgeSlacks();
        for (int i = 0; i < static_cast<int>(variables.size()); ++i) {
            if (!nodes[i]->is_alive) {
                continue;
            }
            if (variables[i] != new_variables[i]) {
                std::cout << i << " " << variables[i] << " " << new_variables[i] << std::endl;
                throw std::runtime_error("variables do not match");
            }
        }
        for (int i = 0; i < static_cast<int>(slacks.size()); ++i) {
            if (slacks[i] != new_slacks[i]) {
                int edge_idx = i;
                std::cout << "edge " << edges[edge_idx]->elementary_head << " " << edges[edge_idx]->elementary_tail <<
                    std::endl;
                std::cout << "edge index: " << i << std::endl;
                std::cout << "head: " << Head(edge_idx) << std::endl;
                std::cout << "tail: " << Tail(edge_idx) << std::endl;
                std::cout << "slacks before/after: " << slacks[i] << " " << new_slacks[i] << std::endl;
                throw std::runtime_error("slacks do not match");
            }
        }
    }

    // shrink phase
    if (params.verbose) {
        std::cout << "shrinking phase" << std::endl;
    }
    std::vector<std::vector<int> > future_blossoms = OrganizeBlossomChildren(record);
    for (std::vector<int> &children : future_blossoms) {
        // if (children.size() > 1) {
        Shrink(children);
        // }
    }

    return action_taken;
}

void VzhuhSolver::MakePrimalUpdate(int edge, PrimalUpdateRecord *record) {
    int parent = Head(edge);
    int child = Tail(edge);

    if (parent == child) {
        return;
    }

    if (nodes[parent]->tree < 0) {
        std::swap(parent, child);
    }

    if (nodes[parent]->tree < 0) {
        // edge not adjacent to any trees
        return;
    }

    if (nodes[parent]->tree == nodes[child]->tree) {
        if (nodes[parent]->plus && nodes[child]->plus) {
            MakeCherryBlossom(edge, record);
            return;
        }
    }

    if (nodes[child]->tree < 0) {
        if (nodes[parent]->plus) {
            Grow(parent, edge, record);
            return;
        }
    }

    if (nodes[parent]->tree != nodes[child]->tree) {
        if (nodes[parent]->plus && nodes[child]->plus) {
            Augment(edge, record);
        }
    }
}

void VzhuhSolver::MakePrimalUpdateForNode(int node, PrimalUpdateRecord *record) {
    if (!nodes[node]->plus || !nodes[node]->is_alive) {
        return;
    }

    for (int edge : NonLoopNeighbors(node)) {
        if (edges[edge]->maybe_has_zero_slack) {
            if (OldSlackQuadrupled(edge) == 0) {
                MakePrimalUpdate(edge, record);
            }
        }
    }
}

void VzhuhSolver::Expand(int blossom, PrimalUpdateRecord *record) {
    ++aux_counter2;

    if (params.verbose) {
        std::cout << "EXPAND " << blossom << std::endl;
    }
    if (nodes[blossom]->old_blossom_parent >= 0) {
        throw std::runtime_error("Expand: two-level expanding");
    }

    nodes[blossom]->is_alive = false;

    RestoreEdgeEndsBeforeExpand(blossom);
    // now OtherEnd is safe

    int new_receptacle = edges[nodes[blossom]->matched_edge]->head;
    if (nodes[new_receptacle]->old_blossom_parent != blossom) {
        new_receptacle = edges[nodes[blossom]->matched_edge]->tail;
    }
    int old_receptacle = Receptacle(new_receptacle);
    int elder_child = edges[nodes[blossom]->minus_parent]->head;
    if (nodes[elder_child]->old_blossom_parent != blossom) {
        elder_child = edges[nodes[blossom]->minus_parent]->tail;
    }

    UpdateInternalStructure(blossom, old_receptacle, new_receptacle, elder_child, record);

    for (int child : nodes[blossom]->blossom_children) {
        AddNodeToRecord(child, record);
        actionable_nodes.push(child);
        if (nodes[child]->tree >= 0) {
            trees[nodes[child]->tree]->tree_nodes.push_back(child);
        }
    }
}

void VzhuhSolver::RestoreEdgeEndsBeforeExpand(int blossom) {
    // makes head and tail of the adjacent edges to point to the children of the blossom,
    // make children top blossoms

    // TODO make better

    for (int child : nodes[blossom]->blossom_children) {
        std::vector<int> elementary_descendants = ElementaryBlossomDescendants(child);
        for (int descendant : elementary_descendants) {
            for (int edge : adj_list[descendant]) {
                if (edges[edge]->elementary_head == descendant) {
                    edges[edge]->head = child;
                } else {
                    edges[edge]->tail = child;
                }
            }
        }
    }

    for (int child : nodes[blossom]->blossom_children) {
        nodes[child]->blossom_parent = -1;
    }
}

void VzhuhSolver::RotateReceptacle(int blossom, int new_receptacle) {
    // changes matched edges and minus_parents

    if (params.verbose) {
        std::cout << "ROTATE RECEPTACLE " << blossom << " " << new_receptacle << std::endl;
    }

    int old_receptacle = Receptacle(new_receptacle);
    if (old_receptacle == new_receptacle) {
        nodes[new_receptacle]->matched_edge = nodes[blossom]->matched_edge;
        return;
    }

    std::vector<int> even_path = EvenPathToReceptacle(new_receptacle);
    std::vector<int> odd_path = OddPathToReceptacle(new_receptacle);

    if (even_path.back() == odd_path.back()) {
        int i = 0;
        int common_node = old_receptacle;
        while (even_path[even_path.size() - 1 - i] == odd_path[odd_path.size() - 1 - i]) {
            common_node = OtherEnd(even_path[even_path.size() - 1 - i], common_node);
            ++i;
        }

        RotateReceptacle(blossom, common_node);
        RotateReceptacle(blossom, new_receptacle);
        return;
    }

    int cur_node = new_receptacle;
    bool match = false;
    for (int edge : even_path) {
        if (match) {
            MakeEdgeMatched(edge);
            cur_node = OtherEnd(edge, cur_node);
        } else {
            MakeEdgeUnmatched(edge);
            nodes[cur_node]->minus_parent = edge;
            cur_node = OtherEnd(edge, cur_node);
            nodes[cur_node]->minus_parent = edge;
        }
        match = !match;
    }

    cur_node = new_receptacle;
    match = false;
    for (int edge : odd_path) {
        if (match) {
            cur_node = OtherEnd(edge, cur_node);
        } else {
            nodes[cur_node]->minus_parent = edge;
            cur_node = OtherEnd(edge, cur_node);
            nodes[cur_node]->minus_parent = edge;
        }
        match = !match;
    }

    nodes[new_receptacle]->minus_parent = -1;
    nodes[new_receptacle]->matched_edge = nodes[blossom]->matched_edge;
    nodes[new_receptacle]->receptacle_ = new_receptacle;
    nodes[old_receptacle]->receptacle_ = new_receptacle;
}

void VzhuhSolver::UpdateMatching(int blossom, int new_receptacle) {
    std::vector<int> path = EvenPathToReceptacle(new_receptacle);
    bool match = false;
    for (int edge : path) {
        if (match) {
            MakeEdgeMatched(edge);
        } else {
            MakeEdgeUnmatched(edge);
        }
        match = !match;
    }
    nodes[new_receptacle]->matched_edge = nodes[blossom]->matched_edge;
}

void VzhuhSolver::UpdateInternalStructure(int blossom,
                                          int old_receptacle,
                                          int new_receptacle,
                                          int elder_child,
                                          PrimalUpdateRecord *record) {
    // after expand, the remaining part of the blossom is a path in the tree

    RotateReceptacle(blossom, new_receptacle);
    std::vector<int> path = EvenPathToReceptacle(elder_child);

    // the path that stays in the tree
    int cur_node = elder_child;
    nodes[elder_child]->minus_parent = nodes[blossom]->minus_parent;
    nodes[elder_child]->plus = false;
    nodes[elder_child]->tree = nodes[blossom]->tree;
    for (int edge : path) {
        cur_node = OtherEnd(edge, cur_node);

        if (edges[edge]->matched) {
            nodes[cur_node]->minus_parent = -1;
            nodes[cur_node]->plus = true;
        } else {
            nodes[cur_node]->minus_parent = edge;
            nodes[cur_node]->plus = false;
        }
        nodes[cur_node]->tree = nodes[blossom]->tree;
    }

    // mark the nodes that stay + insert them into the record
    ++nodes_label_cnt;
    cur_node = elder_child;
    for (int edge : path) {
        nodes[cur_node]->label = nodes_label_cnt;
        AddNodeToRecord(cur_node, record);
        cur_node = OtherEnd(edge, cur_node);
    }
    nodes[new_receptacle]->label = nodes_label_cnt;
    AddNodeToRecord(new_receptacle, record);

    // clear the part that goes to waste
    for (int child : nodes[blossom]->blossom_children) {
        if (nodes[child]->label != nodes_label_cnt) {
            nodes[child]->plus = false;
            nodes[child]->minus_parent = -1;
            nodes[child]->tree = -1;
        }
    }

    // reset receptacles
    for (int child : nodes[blossom]->blossom_children) {
        nodes[child]->receptacle_ = child;
    }
}

std::vector<int> VzhuhSolver::EvenPathToReceptacle(int node) {
    int receptacle = Receptacle(node);
    std::vector<int> path;

    while (node != receptacle) {
        // std::cout << edges[nodes[node]->matched_edge]->head << " " << edges[nodes[node]->matched_edge]->tail << std::endl;
        int parent = OtherEnd(nodes[node]->matched_edge, node);
        int grandparent = OtherEnd(nodes[parent]->minus_parent, parent);

        if (params.verbose) {
            std::cout << "EvenPathToReceptacle, parent/grandparent: " << parent << " " << grandparent << std::endl;
        }

        path.push_back(nodes[node]->matched_edge);
        path.push_back(nodes[parent]->minus_parent);

        node = grandparent;

        if (path.size() > num_vertices_elementary) {
            std::cout << "receptacle: " << receptacle << ", node: " << node << std::endl;

            for (int pathedge : path) {
                std::cout << "(" << Head(pathedge) << ", " << Tail(pathedge) << ") ";
            }
            std::cout << std::endl;
            for (int pathedge : path) {
                std::cout << "(" << edges[pathedge]->elementary_head << ", " << edges[pathedge]->elementary_tail <<
                    ") ";
            }
            std::cout << std::endl;

            throw std::runtime_error("EvenPathToReceptacle: infinite loop");
        }
    }

    return path;
}

std::vector<int> VzhuhSolver::OddPathToReceptacle(int node) {
    int receptacle = Receptacle(node);
    if (receptacle == node) {
        throw std::runtime_error("OddPathToReceptacle: node is receptacle");
    }

    std::vector<int> path = {nodes[node]->minus_parent};
    node = OtherEnd(nodes[node]->minus_parent, node);
    while (node != receptacle) {
        int parent = OtherEnd(nodes[node]->matched_edge, node);
        int grandparent = OtherEnd(nodes[parent]->minus_parent, parent);

        path.push_back(nodes[node]->matched_edge);
        path.push_back(nodes[parent]->minus_parent);

        node = grandparent;
    }

    return path;
}

void VzhuhSolver::ExpandChildBeforeGrow(int blossom, PrimalUpdateRecord *record) {
    ++aux_counter2;

    if (nodes[blossom]->old_blossom_parent >= 0) {
        // avoid two-level expansion
        return;
    }

    nodes[blossom]->is_alive = false;

    RestoreEdgeEndsBeforeExpand(blossom);
    // now OtherEnd, Head, Tail is safe

    int new_receptacle = edges[nodes[blossom]->matched_edge]->head;
    if (nodes[new_receptacle]->old_blossom_parent != blossom) {
        new_receptacle = edges[nodes[blossom]->matched_edge]->tail;
    }
    UpdateMatching(blossom, new_receptacle);

    for (int child : nodes[blossom]->blossom_children) {
        nodes[child]->receptacle_ = child;
        nodes[child]->plus = false;
        nodes[child]->minus_parent = -1;
        nodes[child]->tree = -1;
    }

    for (int child : nodes[blossom]->blossom_children) {
        AddNodeToRecord(child, record);
        if (nodes[child]->tree >= 0) {
            actionable_nodes.push(child);
        }
    }
}

void VzhuhSolver::Grow(int parent, int edge, PrimalUpdateRecord *record) {
    ++aux_counter1;

    int child = OtherEnd(edge, parent);
    int tree = nodes[parent]->tree;

    if (!nodes[parent]->plus) {
        throw std::runtime_error("In Grow: parent is not a plus");
    }
    if (nodes[child]->tree >= 0) {
        throw std::runtime_error("In Grow: child vertex is not free");
    }
    if (nodes[child]->matched_edge < 0) {
        std::cout << "parent " << parent << ", child " << child << std::endl;
        throw std::runtime_error("Grow: child has no matched_edge");
    }

    int grandchild = OtherEnd(nodes[child]->matched_edge, child);
    if (params.verbose) {
        std::cout << "GROW " << parent << " " << child << " " << grandchild << std::endl;
    }

    if (!IsElementary(child)) {
        if (DualVariableQuadrupled(child,
                                   nodes[child]->old_tree,
                                   nodes[child]->old_plus,
                                   nodes[child]->old_blossom_parent) == 0) {
            if (params.verbose) {
                std::cout << "EXPAND CHILD, new child/grandchild: ";
            }
            ExpandChildBeforeGrow(child, record);
            child = OtherEnd(edge, parent);
            grandchild = OtherEnd(nodes[child]->matched_edge, child);
            if (params.verbose) {
                std::cout << child << " " << grandchild << std::endl;
            }
        }
    }

    nodes[child]->minus_parent = edge;
    nodes[child]->tree = tree;
    nodes[child]->plus = false;

    nodes[grandchild]->tree = tree;
    nodes[grandchild]->plus = true;

    AddNodeToRecord(child, record);
    AddNodeToRecord(grandchild, record);

    trees[tree]->tree_nodes.push_back(child);
    trees[tree]->tree_nodes.push_back(grandchild);

    actionable_nodes.push(grandchild);
}

void VzhuhSolver::MakeCherryBlossom(int edge_plus_plus, PrimalUpdateRecord *record) {
    int head = Head(edge_plus_plus);
    int tail = Tail(edge_plus_plus);

    if (Receptacle(head) == Receptacle(tail)) {
        return;
    }

    if (params.verbose) {
        std::cout << "MAKE CHERRY BLOSSOM " << head << " " << tail << std::endl;
    }

    auto [first_bound, second_bound] = CherryPathBounds(head, tail);

    UpdateCherryPath(head, first_bound, record);
    UpdateCherryPath(tail, second_bound, record);

    if (head != first_bound) {
        nodes[head]->minus_parent = edge_plus_plus;
    }
    if (tail != second_bound) {
        nodes[tail]->minus_parent = edge_plus_plus;
    }
}

std::pair<int, int> VzhuhSolver::CherryPathBounds(int first_vertex,
                                                  int second_vertex) {
    int lca = PlusPlusLCA(first_vertex, second_vertex);
    int lca_receptacle = Receptacle(lca);

    while (first_vertex != lca) {
        if (Receptacle(first_vertex) == lca_receptacle) {
            break;
        }
        first_vertex = OtherEnd(nodes[first_vertex]->matched_edge, first_vertex);
        first_vertex = OtherEnd(nodes[first_vertex]->minus_parent, first_vertex);
    }
    while (second_vertex != lca) {
        if (Receptacle(second_vertex) == lca_receptacle) {
            break;
        }
        second_vertex = OtherEnd(nodes[second_vertex]->matched_edge, second_vertex);
        second_vertex = OtherEnd(nodes[second_vertex]->minus_parent, second_vertex);
    }

    return {first_vertex, second_vertex};
}

void VzhuhSolver::UpdateCherryPath(int lower_node, int upper_node, PrimalUpdateRecord *record) {
    int receptacle = Receptacle(upper_node);
    nodes[Receptacle(upper_node)]->receptacle_ = receptacle; // TODO does this line do anything?

    AddNodeToRecord(lower_node, record);
    while (lower_node != upper_node) {
        int parent = OtherEnd(nodes[lower_node]->matched_edge, lower_node);
        int grandparent = OtherEnd(nodes[parent]->minus_parent, parent);

        AddNodeToRecord(parent, record);
        AddNodeToRecord(grandparent, record);

        nodes[Receptacle(lower_node)]->receptacle_ = receptacle;
        nodes[Receptacle(parent)]->receptacle_ = receptacle;

        if (!nodes[parent]->plus) {
            nodes[parent]->plus = true;
            actionable_nodes.push(parent);
        }

        if (grandparent != upper_node) {
            nodes[grandparent]->minus_parent = nodes[parent]->minus_parent;
        }

        lower_node = grandparent;
    }
}

void VzhuhSolver::Augment(int edge_plus_plus, PrimalUpdateRecord *record) {
    num_trees_alive -= 2;

    int head = Head(edge_plus_plus);
    int tail = Tail(edge_plus_plus);

    if (params.verbose) {
        std::cout << "AUGMENT " << head << " " << tail << ", elementary ends: " << edges[edge_plus_plus]->
            elementary_head
            << " " << edges[edge_plus_plus]->elementary_tail <<
            ", is in zero slack set: " << edges[edge_plus_plus]->maybe_has_zero_slack << std::endl;
    }

    std::vector<int> first_path = PathToRoot(head);
    std::vector<int> second_path = PathToRoot(tail);

    std::vector<int> path;
    path.reserve(first_path.size() + second_path.size() + 1);
    for (int i = static_cast<int>(first_path.size()) - 1; i >= 0; --i) {
        path.push_back(first_path[i]);
    }
    path.push_back(edge_plus_plus);
    for (const auto &edge : second_path) {
        path.push_back(edge);
    }

    int first_tree = nodes[head]->tree;
    int second_tree = nodes[tail]->tree;
    AugmentPath(path);
    ClearTree(first_tree, record);
    ClearTree(second_tree, record);
}

std::vector<int> VzhuhSolver::PathToRoot(int node_plus) {
    int root = TopBlossom(trees[nodes[node_plus]->tree]->root);

    std::vector<int> path;
    while (node_plus != root) {
        path.push_back(nodes[node_plus]->matched_edge);
        node_plus = OtherEnd(nodes[node_plus]->matched_edge, node_plus);
        path.push_back(nodes[node_plus]->minus_parent);
        node_plus = OtherEnd(nodes[node_plus]->minus_parent, node_plus);
    }

    return path;
}

void VzhuhSolver::AugmentPath(const std::vector<int> &path) {
    bool match = true;
    for (int edge : path) {
        if (match) {
            MakeEdgeMatched(edge);
        } else {
            MakeEdgeUnmatched(edge);
        }
        match = !match;
    }
}

void VzhuhSolver::ClearTree(int tree, PrimalUpdateRecord *record) {
    trees[tree]->is_alive = false;

    // TODO have a special routine for the last two trees
    // TODO amortize

    for (int node : trees[tree]->tree_nodes) {
        if (nodes[node]->blossom_parent >= 0 || nodes[node]->tree != tree) {
            continue;
        }

        // TODO make better

        AddNodeToRecord(node, record);
        nodes[node]->tree = -1;
        nodes[node]->minus_parent = -1;
        nodes[node]->receptacle_ = node;
        nodes[node]->plus = false;
    }
}

void VzhuhSolver::UpdateQueues(const PrimalUpdateRecord &record) {
    // updates amortized slams and variables,
    // edge_heaps, node_heaps
    // old_plus, old_tree, old_blossom_parent, is_in_record for nodes

    std::vector<int> edges_to_update;
    int num_edges_to_update = 0;
    for (int node : record.changed_sign) {
        if (!nodes[node]->is_alive) {
            continue;
        }
        if (nodes[node]->old_blossom_parent < 0 &&
            nodes[node]->old_plus == nodes[node]->plus &&
            nodes[node]->old_tree == nodes[node]->tree) {
            continue;
        }

        num_edges_to_update += NonLoopNeighbors(node).size();
    }
    edges_to_update.reserve(num_edges_to_update);

    for (int node : record.changed_sign) {
        if (!nodes[node]->is_alive) {
            continue;
        }
        if (nodes[node]->blossom_parent >= 0) {
            throw std::runtime_error("UpdateQueues: non-top node in record");
        }
        if (nodes[node]->old_blossom_parent < 0 &&
            nodes[node]->old_plus == nodes[node]->plus &&
            nodes[node]->old_tree == nodes[node]->tree) {
            continue;
        }

        // update dual_var_quadrupled_amortized_
        int old_variable = DualVariableQuadrupled(node,
                                                  nodes[node]->old_tree,
                                                  nodes[node]->old_plus,
                                                  nodes[node]->old_blossom_parent);
        int new_variable = DualVariableQuadrupled(node);
        nodes[node]->dual_var_quadrupled_amortized_ += (old_variable - new_variable);
        if (params.verbose) {
            std::cout << "updating node " << node << " old var, would-be-new var: " << old_variable << " " <<
                new_variable << std::endl;
        }

        // remember the edges to update
        int diff = 0;
        if (nodes[node]->tree >= 0) {
            if (nodes[node]->plus) {
                diff += trees[nodes[node]->tree]->dual_var_quadrupled;
            } else {
                diff -= trees[nodes[node]->tree]->dual_var_quadrupled;
            }
        }
        int old_top_node = node;
        if (nodes[node]->old_blossom_parent >= 0) {
            old_top_node = nodes[node]->old_blossom_parent;
        }
        if (nodes[old_top_node]->old_tree >= 0) {
            if (nodes[old_top_node]->old_plus) {
                diff -= trees[nodes[old_top_node]->old_tree]->dual_var_quadrupled;
            } else {
                diff += trees[nodes[old_top_node]->old_tree]->dual_var_quadrupled;
            }
        }
        for (int edge : NonLoopNeighbors(node)) {
            edges[edge]->slack_diff += diff;
            if (!edges[edge]->must_be_updated) {
                edges_to_update.push_back(edge);
                edges[edge]->must_be_updated = true;
            }
            if (old_top_node == node) {
                edges[edge]->maybe_was_loop = false;
            }
        }

        // update minus_blossoms queues
        if (!IsElementary(node)) {
            if (nodes[node]->old_plus != nodes[node]->plus || nodes[node]->old_tree != nodes[node]->tree) {
                // the status might have changed
                if (nodes[node]->old_tree >= 0 && !nodes[node]->old_plus) {
                    // were a minus blossom
                    RemoveNodeFromQueue(node);
                }
                if (nodes[node]->tree >= 0 && !nodes[node]->plus) {
                    // became a minus blossom
                    AddNodeToQueue(node, trees[nodes[node]->tree]->minus_blossoms);
                }
            }
        }
    }

    // update the edges
    for (int edge : edges_to_update) {
        UpdateEdgeSlack(edge);
    }

    // update old_plus, old_tree, old_blossom_parent
    for (int node : record.changed_sign) {
        nodes[node]->old_blossom_parent = nodes[node]->blossom_parent;
        nodes[node]->old_plus = nodes[node]->plus;
        nodes[node]->old_tree = nodes[node]->tree;
        nodes[node]->is_in_record = false;
    }
}

void VzhuhSolver::UpdateEdgeSlack(int edge) {
    RemoveEdgeFromQueue(edge);

    if (edges[edge]->maybe_was_loop) {
        int old_slack = OldSlackQuadrupled(edge);
        int new_slack = SlackQuadrupled(edge);
        edges[edge]->slack_quadrupled_amortized_ += (old_slack - new_slack);
    } else {
        edges[edge]->slack_quadrupled_amortized_ += edges[edge]->slack_diff;
    }

    edges[edge]->maybe_was_loop = true;
    edges[edge]->slack_diff = 0;
    edges[edge]->must_be_updated = false;

    if (Receptacle(Head(edge)) != Receptacle(Tail(edge))) {
        AddEdgeToQueue(edge);
    }
}

std::vector<std::vector<int> > VzhuhSolver::OrganizeBlossomChildren(const PrimalUpdateRecord &record) {
    ++nodes_label_cnt;
    int label_zero = nodes_label_cnt;

    for (int node : record.changed_sign) {
        if (nodes[node]->tree < 0) {
            continue;
        }
        if (nodes[node]->label >= label_zero) {
            // already seen this node
            continue;
        }

        int receptacle = Receptacle(node);
        if (nodes[receptacle]->label < label_zero) {
            nodes[receptacle]->label = nodes_label_cnt;
            ++nodes_label_cnt;
        }
        if (nodes[node]->label < label_zero) {
            nodes[node]->label = nodes[receptacle]->label;
        }
    }

    std::vector<int> sizes(nodes_label_cnt - label_zero, 0);
    std::vector<int> label_diff_to_index(nodes_label_cnt - label_zero, -1);

    for (int node : record.changed_sign) {
        if (nodes[node]->label < label_zero) {
            continue;
        }
        ++sizes[nodes[node]->label - label_zero];
    }

    int index = 0;
    for (int i = 0; i < static_cast<int>(sizes.size()); ++i) {
        if (sizes[i] > 1) {
            label_diff_to_index[i] = index;
            ++index;
        }
    }

    std::vector<std::vector<int> > result(index);
    for (int i = 0; i < static_cast<int>(sizes.size()); ++i) {
        if (label_diff_to_index[i] >= 0) {
            result[label_diff_to_index[i]].reserve(sizes[i]);
        }
    }

    for (int node : record.changed_sign) {
        if (nodes[node]->label < label_zero) {
            continue;
        }
        if (label_diff_to_index[nodes[node]->label - label_zero] >= 0) {
            result[label_diff_to_index[nodes[node]->label - label_zero]].push_back(node);
        }
    }

    return result;
}

void VzhuhSolver::Shrink(std::vector<int> &children) {
    ++aux_counter3;

    if (params.verbose) {
        std::cout << "SHRINK ";
        for (int node : children) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }

    if (children.size() % 2 == 0) {
        throw std::runtime_error("In Node: blossom has to have an odd number of vertices");
    }

    int new_index = static_cast<int>(nodes.size());
    int receptacle = Receptacle(children.front());

    nodes.emplace_back(std::make_unique<Node>(new_index));
    nodes.back()->plus = true;
    nodes.back()->old_plus = true;
    nodes.back()->blossom_children = std::move(children);
    nodes.back()->matched_edge = nodes[receptacle]->matched_edge;
    nodes.back()->tree = nodes[receptacle]->tree;
    nodes.back()->old_tree = nodes.back()->tree;
    nodes.back()->dual_var_quadrupled_amortized_ = -trees[nodes.back()->tree]->dual_var_quadrupled;
    nodes.back()->tree_var_at_birth = trees[nodes.back()->tree]->dual_var_quadrupled;

    // update blossom_parent of the children, set labels to empty, update variables
    for (int child : nodes.back()->blossom_children) {
        if (nodes[child]->blossom_parent >= 0) {
            throw std::runtime_error("In Node: some child node already has a parent");
        }
        nodes[child]->blossom_parent = new_index;
        nodes[child]->old_blossom_parent = new_index;
        nodes[child]->dual_var_quadrupled_amortized_ += trees[nodes.back()->tree]->dual_var_quadrupled;
    }

    // add the new blossom to the list of vertices in the tree
    trees[nodes.back()->tree]->tree_nodes.push_back(new_index);
}

void VzhuhSolver::ValidateQueues() {
    // checks that the queues are in a correct state, throws if not

    // every queue must hold only the edges that belong to this queue

    // for (int tree : alive_trees) {
    //     if (!trees[tree]->is_alive) {
    //         continue;
    //     }
    //
    //     // check plus empty
    //     edge_heaps[trees[tree]->plus_empty_edges]->ValidateHeap("plus empty");
    //     for (auto heap_node : edge_heaps[trees[tree]->plus_empty_edges]->heap_) {
    //         int edge = heap_node.value;
    //
    //         int first = Head(edge); // plus
    //         int second = Tail(edge); // empty
    //         if (nodes[first]->tree != tree) {
    //             std::swap(first, second);
    //         }
    //         if (nodes[first]->tree != tree || !nodes[first]->plus) {
    //             throw std::runtime_error("Incorrect plus empty queue");
    //         }
    //         if (nodes[second]->tree >= 0) {
    //             throw std::runtime_error("Incorrect plus empty queue");
    //         }
    //     }
    //
    //     // check (plus, plus) internal
    //     edge_heaps[trees[tree]->plus_plus_internal_edges]->ValidateHeap("+ + int");
    //     for (auto heap_node : edge_heaps[trees[tree]->plus_plus_internal_edges]->heap_) {
    //         int edge = heap_node.value;
    //
    //         int first = Head(edge);
    //         int second = Tail(edge);
    //
    //         if (nodes[first]->tree != tree || !nodes[first]->plus) {
    //             throw std::runtime_error("Incorrect plus plus internal queue");
    //         }
    //         if (nodes[second]->tree != tree || !nodes[second]->plus) {
    //             throw std::runtime_error("Incorrect plus plus internal queue");
    //         }
    //     }
    //
    //     // check minus blossoms
    //     node_heaps[trees[tree]->minus_blossoms]->ValidateHeap("minus blossoms");
    //     for (auto heap_node : node_heaps[trees[tree]->minus_blossoms]->heap_) {
    //         int node = heap_node.value;
    //
    //         if (IsElementary(node) || (nodes[node]->tree != tree) || nodes[node]->plus) {
    //             throw std::runtime_error("Incorrect minus blossom queue");
    //         }
    //     }
    //
    //     // check (plus, plus) external
    //     for (auto [other_tree, queue_index] : trees[tree]->pq_plus_plus) {
    //         if (!trees[other_tree]->is_alive) {
    //             continue;
    //         }
    //         for (auto heap_node : edge_heaps[queue_index]->heap_) {
    //             int edge = heap_node.value;
    //
    //             int first = Head(edge); // in this tree
    //             int second = Tail(edge); // in the other tree
    //             if (nodes[first]->tree != tree) {
    //                 std::swap(first, second);
    //             }
    //             if (nodes[first]->tree != tree || !nodes[first]->plus) {
    //                 throw std::runtime_error("Incorrect plus plus queue");
    //             }
    //             if (nodes[second]->tree != other_tree || !nodes[second]->plus) {
    //                 throw std::runtime_error("Incorrect plus plus queue");
    //             }
    //         }
    //         edge_heaps[queue_index]->ValidateHeap("+ + external");
    //     }
    //
    //     // check plus minus external
    //     for (auto [other_tree, queue_index] : trees[tree]->pq_plus_minus) {
    //         if (!trees[other_tree]->is_alive) {
    //             continue;
    //         }
    //         edge_heaps[queue_index]->ValidateHeap("+ - ext");
    //         for (auto heap_node : edge_heaps[queue_index]->heap_) {
    //             int edge = heap_node.value;
    //
    //             int first = Head(edge); // in this tree
    //             int second = Tail(edge); // in the other tree
    //             if (nodes[first]->tree != tree) {
    //                 std::swap(first, second);
    //             }
    //             if (nodes[first]->tree != tree || !nodes[first]->plus) {
    //                 throw std::runtime_error("Incorrect plus minus queue");
    //             }
    //             if (nodes[second]->tree != other_tree || nodes[second]->plus) {
    //                 throw std::runtime_error("Incorrect plus minus queue");
    //             }
    //         }
    //     }
    //
    //     // check minus plus external
    //     for (auto [other_tree, queue_index] : trees[tree].pq_minus_plus) {
    //         if (!trees[other_tree]->is_alive) {
    //             continue;
    //         }
    //         edge_heaps[queue_index]->ValidateHeap("- + ext");
    //         for (auto heap_node : edge_heaps[queue_index]->heap_) {
    //             int edge = heap_node.value;
    //
    //             int first = Head(edge); // in this tree
    //             int second = Tail(edge); // in the other tree
    //             if (nodes[first]->tree != tree) {
    //                 std::swap(first, second);
    //             }
    //             if (nodes[first]->tree != tree || nodes[first]->plus) {
    //                 throw std::runtime_error("Incorrect minus plus queue");
    //             }
    //             if (nodes[second]->tree != other_tree || !nodes[second]->plus) {
    //                 throw std::runtime_error("Incorrect minus plus queue");
    //             }
    //         }
    //     }
    // }
}

std::vector<int> VzhuhSolver::NodeVariables() const {
    std::vector<int> result;
    for (int i(0); i < static_cast<int>(nodes.size()); ++i) {
        if (nodes[i]->is_alive) {
            result.push_back(DualVariableQuadrupled(i));
        } else {
            result.push_back(0);
        }
    }
    return result;
}

std::vector<int> VzhuhSolver::EdgeSlacks() {
    std::vector<int> depths(nodes.size(), 0);
    for (int node(nodes.size() - 1); node >= 0; --node) {
        if (!nodes[node]->is_alive) {
            continue;
        }
        if (nodes[node]->blossom_parent >= 0) {
            depths[node] = depths[nodes[node]->blossom_parent] + 1;
        }
    }

    std::vector<int> result;
    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
        if (Head(i) != Tail(i)) {
            result.push_back(SlackQuadrupled(i));
        } else {
            int edge(i);

            int head = edges[edge]->elementary_head;
            int tail = edges[edge]->elementary_tail;

            while (head != tail) {
                if (depths[head] > depths[tail]) {
                    head = nodes[head]->blossom_parent;
                } else {
                    tail = nodes[tail]->blossom_parent;
                }
            }

            result.push_back(edges[edge]->slack_quadrupled_amortized_ - 2 * nodes[head]->tree_var_at_birth);
        }
    }

    return result;
}

void VzhuhSolver::ValidateEvenOddPaths() {
    for (int node(0); node < static_cast<int>(nodes.size()); ++node) {
        if (!nodes[node]->is_alive || nodes[node]->blossom_parent >= 0) {
            continue;
        }

        int receptacle = Receptacle(node);

        auto even_path = EvenPathToReceptacle(node);
        if (even_path.size() % 2 != 0) {
            throw std::runtime_error("ValidateEvenOddPaths: even path is not even");
        }
        int cur_node = node;
        for (int edge : even_path) {
            cur_node = OtherEnd(edge, cur_node);
        }
        if (cur_node != receptacle) {
            throw std::runtime_error("ValidateEvenOddPaths: incorrect even path");
        }

        if (node != receptacle) {
            auto odd_path = OddPathToReceptacle(node);
            if (odd_path.size() % 2 != 1) {
                throw std::runtime_error("ValidateEvenOddPaths: odd path is not odd");
            }
            cur_node = node;
            for (int edge : odd_path) {
                cur_node = OtherEnd(edge, cur_node);
            }
            if (cur_node != receptacle) {
                throw std::runtime_error("ValidateEvenOddPaths: incorrect odd path");
            }
        }
    }
}

bool VzhuhSolver::MakeDualUpdates() {
    UpdateAliveTreesList();

    if (num_trees_alive == 0) {
        return true;
    }

    if (params.verbose) {
        std::cout << "DUAL UPDATE" << std::endl;
    }

    std::vector<int> deltas = VariableDeltas();
    for (int i = 0; i < static_cast<int>(alive_trees.size()); ++i) {
        trees[alive_trees[i]]->dual_var_quadrupled += deltas[i];
    }

    InitNextRoundActionable();

    if (params.verbose) {
        std::cout << "after dual update:" << std::endl;
        PrintGraph();
    }

    return !std::all_of(deltas.begin(), deltas.end(), [](auto x) { return x == 0; });
}

void VzhuhSolver::UpdateAliveTreesList() {
    for (int i = 0; i < static_cast<int>(alive_trees.size()); ++i) {
        if (!trees[alive_trees[i]]->is_alive) {
            alive_trees[i] = alive_trees.back();
            alive_trees.pop_back();
            --i;
        }
    }
}

std::vector<int> VzhuhSolver::VariableDeltas() {
    DualConstraints dual_constraints = GetDualConstraints();

    // for (auto b : dual_constraints.upper_bound) {
    //     if (b < 0) {
    //         throw std::runtime_error("b < 0");
    //     }
    // }
    // for (auto c : dual_constraints->plus_plus_constraints) {
    //     for (auto [ind, sl] : c) {
    //         if (sl < 0) {
    //             throw std::runtime_error("constraint plus plus < 0");
    //         }
    //     }
    // }
    // for (auto c : dual_constraints->plus_minus_constraints) {
    //     for (auto [ind, sl] : c) {
    //         if (sl < 0) {
    //             throw std::runtime_error("constraint plus minus < 0");
    //         }
    //     }
    // }

    std::vector<std::vector<int> > connected_components = ConnectedComponentsTreeTree(dual_constraints);
    if (params.verbose) {
        std::cout << "num of connected_components : " << connected_components.size() << std::endl;
    }
    std::vector<int> deltas(alive_trees.size(), 0);

    std::vector<int> component_index(alive_trees.size(), 0);
    for (int i = 0; i < static_cast<int>(connected_components.size()); ++i) {
        for (int v : connected_components[i]) {
            component_index[v] = i;
        }
    }

    for (int cc_index = 0; cc_index < static_cast<int>(connected_components.size()); ++cc_index) {
        int delta = INT32_MAX;

        // find delta for this connected component
        for (int v : connected_components[cc_index]) {
            if (dual_constraints.upper_bound[v] < delta) {
                delta = dual_constraints.upper_bound[v];
            }

            for (auto [w, slack] : dual_constraints.plus_plus_constraints[v]) {
                if (component_index[w] == cc_index) {
                    if (slack % 2 != 0) {
                        throw std::runtime_error("In VariableDeltas: slack is not divisible by 2");
                    }
                    if (slack / 2 < delta) {
                        delta = slack / 2;
                    }
                } else {
                    if (slack - deltas[w] < delta) {
                        delta = slack - deltas[w];
                    }
                }
            }

            for (auto [w, slack] : dual_constraints.plus_minus_constraints[v]) {
                if (component_index[w] == cc_index) {
                    continue;
                }
                if (slack + deltas[w] < delta) {
                    delta = slack + deltas[w];
                }
            }
        }

        // apply delta to this connected component
        for (int v : connected_components[cc_index]) {
            deltas[v] = delta;
        }
    }

    for (int d : deltas) {
        if (d < 0) {
            throw std::runtime_error("In VariableDeltas: found negative delta");
        }
    }

    return deltas;
}

std::vector<std::vector<int> > VzhuhSolver::ConnectedComponentsTreeTree(const DualConstraints &dual_constraints) const {
    // returns a vector of index_of_connected_component

    std::vector<bool> visited(alive_trees.size(), false);
    std::vector<std::vector<int> > components;

    // Build undirected adjacency list for weak connectivity
    std::vector<std::vector<int> > adj_list_tree_tree(alive_trees.size(), std::vector<int>());
    for (int u = 0; u < static_cast<int>(alive_trees.size()); ++u) {
        for (auto &[v, slack] : dual_constraints.plus_minus_constraints[u]) {
            if (slack == 0) {
                adj_list_tree_tree[u].emplace_back(v);
                adj_list_tree_tree[v].emplace_back(u);
            }
        }
    }

    for (int start = 0; start < static_cast<int>(alive_trees.size()); ++start) {
        if (!visited[start]) {
            std::vector<int> comp;
            std::queue<int> q;

            visited[start] = true;
            q.push(start);

            while (!q.empty()) {
                int u = q.front();
                q.pop();
                comp.push_back(u);

                for (int v : adj_list_tree_tree[u]) {
                    if (!visited[v]) {
                        visited[v] = true;
                        q.push(v);
                    }
                }
            }

            components.push_back(std::move(comp));
        }
    }

    return components;
}

VzhuhSolver::DualConstraints VzhuhSolver::GetDualConstraints() {
    // get upper_bound
    std::vector<int> upper_bound(alive_trees.size(), INT32_MAX);
    for (int i = 0; i < static_cast<int>(alive_trees.size()); ++i) {
        int tree = alive_trees[i];

        int plus_empty = PlusEmptySlack(tree);
        int plus_plus_internal = PlusPlusInternalSlack(tree);
        int blossom_var = MinMinusBlossomVariable(tree);

        // if (plus_empty < 0) {
        //     throw std::runtime_error("In VariableDeltas: plus_empty is negative");
        // }
        // if (plus_plus_internal < 0) {
        //     throw std::runtime_error("In VariableDeltas: plus_plus_internal is not positive");
        // }
        // if (blossom_var < 0) {
        //     throw std::runtime_error("In VariableDeltas: blossom_var is not positive");
        // }
        // if (plus_plus_internal < INT32_MAX && plus_plus_internal % 2 != 0) {
        //     throw std::runtime_error("PlusPlusInternalSlack is not divisible by 2");
        // }

        if (plus_empty < upper_bound[i]) {
            upper_bound[i] = plus_empty;
        }
        if (plus_plus_internal / 2 < upper_bound[i]) {
            upper_bound[i] = plus_plus_internal / 2;
        }
        if (blossom_var < upper_bound[i]) {
            upper_bound[i] = blossom_var;
        }
    }

    // update alive indices
    for (int i = 0; i < static_cast<int>(alive_trees.size()); ++i) {
        trees[alive_trees[i]]->alive_index = i;
    }

    // get plus_plus_constraints and plus_minus_constraints
    std::vector plus_plus_constraints(alive_trees.size(), std::vector<std::pair<int, int> >());
    std::vector plus_minus_constraints(alive_trees.size(), std::vector<std::pair<int, int> >());
    for (int i = 0; i < static_cast<int>(alive_trees.size()); ++i) {
        int tree = alive_trees[i];

        std::vector<std::pair<int, int> > plus_plus_slacks = PlusPlusExternalSlacks(tree);
        std::vector<std::pair<int, int> > plus_minus_slacks = PlusMinusExternalSlacks(tree);
        for (auto [other_tree, slack] : plus_plus_slacks) {
            plus_plus_constraints[i].emplace_back(trees[other_tree]->alive_index, slack);
        }
        for (auto [other_tree, slack] : plus_minus_slacks) {
            plus_minus_constraints[i].emplace_back(trees[other_tree]->alive_index, slack);
        }
    }

    return std::move(DualConstraints{
        std::move(upper_bound), std::move(plus_plus_constraints), std::move(plus_minus_constraints)
    });
}

void VzhuhSolver::InitNextRoundActionable() {
    for (int tree : alive_trees) {
        if (params.verbose) {
            std::cout << "updating zero slack set next to tree " << tree << std::endl;
        }

        for (auto [other_tree, queue_index] : trees[tree]->pq_plus_plus) {
            if (tree < other_tree) {
                // to avoid double work
                AddZeroSlackEdgesFromQueue(queue_index, true);
            }
        }
        for (auto [_, queue_index] : trees[tree]->pq_plus_minus) {
            AddZeroSlackEdgesFromQueue(queue_index, false);
        }
        AddZeroSlackEdgesFromQueue(trees[tree]->plus_empty_edges, true);
        CleanLoopsFromQueueTop(tree);
        AddZeroSlackEdgesFromQueue(trees[tree]->plus_plus_internal_edges, true);
    }
}

void VzhuhSolver::AddZeroSlackEdgesFromQueue(int queue_index, bool add_to_actionable) {
    Edge *top = edge_heaps[queue_index]->GetMin();
    if (!top) {
        return;
    }
    if (SlackQuadrupled(top->index) > 0) {
        return;
    }

    edge_heaps[queue_index]->RemoveMin();
    edge_heaps[queue_index]->Insert(top);

    int min_key = top->slack_quadrupled_amortized_;

    std::vector<Edge *> stack;
    stack.push_back(edge_heaps[queue_index]->GetMin());

    while (!stack.empty()) {
        Edge *edge = stack.back();
        stack.pop_back();

        if (!edges[edge->index]->maybe_has_zero_slack) {
            edges[edge->index]->maybe_has_zero_slack = true;
        }
        if (add_to_actionable) {
            actionable_edges.push(edge->index);
        }

        for (Edge *child = edge->heap_child; child; child = child->heap_next) {
            if (child->Key() == min_key) {
                stack.push_back(child);
            }
        }
    }
}

void VzhuhSolver::CleanLoopsFromQueueTop(int tree) {
    int queue_idx = trees[tree]->plus_plus_internal_edges;

    while (edge_heaps[queue_idx]->GetMin()) {
        Edge *top = edge_heaps[queue_idx]->GetMin();
        if (Head(top->index) == Tail(top->index)) {
            RemoveEdgeFromQueue(top->index);
        } else {
            break;
        }
    }
}

bool VzhuhSolver::IsElementary(int node) const {
    return node < num_vertices_elementary;
}

int VzhuhSolver::TopBlossom(int node) const {
    while (nodes[node]->blossom_parent >= 0) {
        node = nodes[node]->blossom_parent;
    }
    return node;
}

int VzhuhSolver::Receptacle(int node) {
    int grandparent = node;
    while (nodes[grandparent]->receptacle_ != grandparent) {
        grandparent = nodes[grandparent]->receptacle_;
    }

    int cur_node = node;
    while (nodes[cur_node]->receptacle_ != grandparent) {
        int next_node = nodes[cur_node]->receptacle_;
        nodes[cur_node]->receptacle_ = grandparent;
        cur_node = next_node;
    }

    return grandparent;
}

int VzhuhSolver::DualVariableQuadrupled(int node) const {
    return DualVariableQuadrupled(node,
                                  nodes[node]->tree,
                                  nodes[node]->plus,
                                  nodes[node]->blossom_parent);
}

int VzhuhSolver::DualVariableQuadrupled(int node, int tree, bool plus, int blossom_parent) const {
    // if (!nodes[node]->is_alive) {
    //     throw std::runtime_error("DualVariableQuadrupled: node is not alive");
    // }

    if (tree < 0 || blossom_parent >= 0) {
        return nodes[node]->dual_var_quadrupled_amortized_;
    }
    if (plus) {
        return nodes[node]->dual_var_quadrupled_amortized_ + trees[tree]->dual_var_quadrupled;
    }
    return nodes[node]->dual_var_quadrupled_amortized_ - trees[tree]->dual_var_quadrupled;
}

std::vector<int> &VzhuhSolver::NonLoopNeighbors(int node) {
    // TODO avoid calling Head and Tail for every edge on the list

    if (IsElementary(node)) {
        return adj_list[node];
    }

    if (!nodes[node]->neighbors.empty()) {
        return nodes[node]->neighbors;
    }

    // otherwise, compute neighbors

    // mark the vertices, collect the lists
    std::queue<int> queue;
    std::vector<std::vector<int> *> lists;
    int total_length = 0;
    queue.push(node);
    ++nodes_label_cnt;
    nodes[node]->label = nodes_label_cnt;
    while (!queue.empty()) {
        int cur = queue.front();
        queue.pop();

        for (int child : nodes[cur]->blossom_children) {
            nodes[child]->label = nodes_label_cnt;

            if (IsElementary(child)) {
                lists.push_back(&adj_list[child]);
                total_length += adj_list[child].size();
            } else if (!nodes[child]->neighbors.empty()) {
                lists.push_back(&nodes[child]->neighbors);
                total_length += nodes[child]->neighbors.size();
            } else {
                queue.push(child);
            }
        }
    }

    // add non-loops to neighbors
    nodes[node]->neighbors.reserve(total_length);
    for (auto list_ptr : lists) {
        for (int edge : *list_ptr) {
            while (nodes[edges[edge]->head]->blossom_parent >= 0 && nodes[edges[edge]->head]->label !=
                nodes_label_cnt) {
                edges[edge]->head = nodes[edges[edge]->head]->blossom_parent;
            }
            while (nodes[edges[edge]->tail]->blossom_parent >= 0 && nodes[edges[edge]->tail]->label !=
                nodes_label_cnt) {
                edges[edge]->tail = nodes[edges[edge]->tail]->blossom_parent;
            }
            if (nodes[edges[edge]->head]->label != nodes_label_cnt || nodes[edges[edge]->tail]->label !=
                nodes_label_cnt) {
                nodes[node]->neighbors.push_back(edge);
            }
        }
    }

    return nodes[node]->neighbors;
}

std::vector<int> &VzhuhSolver::NonLoopZeroSlackNeighbors(int node) {
    ++aux_counter1;

    // if (IsElementary(node)) {
    //     return zero_slack_adj_list[node ];
    // }

    if (nodes[node]->round_0slack_neighbors_updated == current_round) {
        return nodes[node]->zero_slack_neighbors;
    }

    ++aux_counter2;

    for (int edge : NonLoopNeighbors(node)) {
        if (edges[edge]->maybe_has_zero_slack) {
            nodes[node]->zero_slack_neighbors.push_back(edge);
        }
    }

    return nodes[node]->zero_slack_neighbors;
}

std::vector<int> VzhuhSolver::ElementaryBlossomDescendants(int node) const {
    if (IsElementary(node)) {
        return {node};
    }

    std::vector<int> elementary_descendants;
    std::queue<int> queue;
    queue.push(node);
    while (!queue.empty()) {
        int cur = queue.front();
        queue.pop();
        for (int child : nodes[cur]->blossom_children) {
            if (IsElementary(child)) {
                elementary_descendants.push_back(child);
            } else {
                queue.push(child);
            }
        }
    }
    return elementary_descendants;
}

int VzhuhSolver::SlackQuadrupled(int edge) {
    int head = Head(edge);
    int tail = Tail(edge);

    // if (head == tail) {
    //     throw std::runtime_error("SlackQuadrupled called for a loop");
    // }

    int slack = edges[edge]->slack_quadrupled_amortized_;

    if (nodes[head]->tree >= 0) {
        if (nodes[head]->plus) {
            slack -= trees[nodes[head]->tree]->dual_var_quadrupled;
        } else {
            slack += trees[nodes[head]->tree]->dual_var_quadrupled;
        }
    }
    if (nodes[tail]->tree >= 0) {
        if (nodes[tail]->plus) {
            slack -= trees[nodes[tail]->tree]->dual_var_quadrupled;
        } else {
            slack += trees[nodes[tail]->tree]->dual_var_quadrupled;
        }
    }

    return slack;
}

int VzhuhSolver::OldSlackQuadrupled(int edge) {
    int head = Head(edge);
    int tail = Tail(edge);

    // if (head == tail) {
    //     throw std::runtime_error("OldSlackQuadrupled called for a loop");
    // }

    if (nodes[head]->old_blossom_parent >= 0) {
        head = nodes[head]->old_blossom_parent;
    }
    if (nodes[tail]->old_blossom_parent >= 0) {
        tail = nodes[tail]->old_blossom_parent;
    }

    if (head == tail) {
        // used to be a loop, but the node got expanded
        return edges[edge]->slack_quadrupled_amortized_ - 2 * nodes[head]->tree_var_at_birth;
    }

    int slack = edges[edge]->slack_quadrupled_amortized_;

    if (nodes[head]->old_tree >= 0) {
        if (nodes[head]->old_plus) {
            slack -= trees[nodes[head]->old_tree]->dual_var_quadrupled;
        } else {
            slack += trees[nodes[head]->old_tree]->dual_var_quadrupled;
        }
    }
    if (nodes[tail]->old_tree >= 0) {
        if (nodes[tail]->old_plus) {
            slack -= trees[nodes[tail]->old_tree]->dual_var_quadrupled;
        } else {
            slack += trees[nodes[tail]->old_tree]->dual_var_quadrupled;
        }
    }

    return slack;
}

int VzhuhSolver::OtherEnd(int edge, int node) {
    // if (Head(edge) == Tail(edge)) {
    //     std::cout << Head(edge) << " " << Tail(edge) << std::endl;
    //     throw std::runtime_error("In OtherEnd: querying for a loop");
    // }

    if (node == Tail(edge)) {
        return Head(edge);
    }
    if (node == Head(edge)) {
        return Tail(edge);
    }
    throw std::runtime_error("In OtherEnd");
}

int VzhuhSolver::OtherElementaryEnd(int edge, int node) const {
    if (node == edges[edge]->elementary_head) {
        return edges[edge]->elementary_tail;
    }
    return edges[edge]->elementary_head;
}

int VzhuhSolver::Head(int edge) {
    while (nodes[edges[edge]->head]->blossom_parent >= 0) {
        edges[edge]->head = nodes[edges[edge]->head]->blossom_parent;
    }
    return edges[edge]->head;
}

int VzhuhSolver::Tail(int edge) {
    while (nodes[edges[edge]->tail]->blossom_parent >= 0) {
        edges[edge]->tail = nodes[edges[edge]->tail]->blossom_parent;
    }
    return edges[edge]->tail;
}

int VzhuhSolver::PlusPlusLCA(int first_vertex, int second_vertex) {
    ++nodes_label_cnt;
    first_vertex = Receptacle(first_vertex);
    second_vertex = Receptacle(second_vertex);

    nodes[first_vertex]->label = nodes_label_cnt;
    nodes[second_vertex]->label = nodes_label_cnt;

    while (first_vertex != second_vertex) {
        if (nodes[first_vertex]->matched_edge >= 0) {
            first_vertex = OtherEnd(nodes[first_vertex]->matched_edge, first_vertex);
            first_vertex = OtherEnd(nodes[first_vertex]->minus_parent, first_vertex);
            first_vertex = Receptacle(first_vertex);
            if (nodes[first_vertex]->label == nodes_label_cnt) {
                return first_vertex;
            }
            nodes[first_vertex]->label = nodes_label_cnt;
        }

        if (nodes[second_vertex]->matched_edge >= 0) {
            second_vertex = OtherEnd(nodes[second_vertex]->matched_edge, second_vertex);
            second_vertex = OtherEnd(nodes[second_vertex]->minus_parent, second_vertex);
            second_vertex = Receptacle(second_vertex);
            if (nodes[second_vertex]->label == nodes_label_cnt) {
                return second_vertex;
            }
            nodes[second_vertex]->label = nodes_label_cnt;
        }
    }

    return first_vertex;
}

void VzhuhSolver::MakeEdgeMatched(int edge) {
    if (Head(edge) == Tail(edge)) {
        throw std::runtime_error("In MakeEdgeMatched: trying to match a loop");
    }
    nodes[Head(edge)]->matched_edge = edge;
    nodes[Tail(edge)]->matched_edge = edge;
    edges[edge]->matched = true;
}

void VzhuhSolver::MakeEdgeUnmatched(int edge) {
    edges[edge]->matched = false;
}

void VzhuhSolver::AddEdgeToQueue(int edge) {
    // adds the edge to the right queue

    int head = Head(edge);
    int tail = Tail(edge);

    // same tree or (empty, empty)
    if (nodes[head]->tree == nodes[tail]->tree) {
        if (nodes[head]->tree < 0) {
            // (empty, empty)
            RemoveEdgeFromQueue(edge);
            return;
        }
        if (nodes[head]->plus && nodes[tail]->plus) {
            // (+, +)
            AddEdgeToThisQueue(edge, trees[nodes[head]->tree]->plus_plus_internal_edges);
            return;
        }

        // (+, -), (-, +), (-, -)
        RemoveEdgeFromQueue(edge);
        return;
    }

    // different trees, make head always in some tree
    if (nodes[head]->tree < 0) {
        std::swap(head, tail);
    }

    // (+, empty)
    if (nodes[head]->plus && nodes[tail]->tree < 0) {
        AddEdgeToThisQueue(edge, trees[nodes[head]->tree]->plus_empty_edges);
        return;
    }

    // (+, +)
    if (nodes[head]->plus && nodes[tail]->plus) {
        AddPQPlusPlus(nodes[head]->tree, nodes[tail]->tree, edge);
        return;
    }

    // (+, -)
    if (nodes[head]->plus && !nodes[tail]->plus) {
        AddPQPlusMinus(nodes[head]->tree, nodes[tail]->tree, edge);
        return;
    }

    // (-, +)
    if (!nodes[head]->plus && nodes[tail]->plus) {
        AddPQPlusMinus(nodes[tail]->tree, nodes[head]->tree, edge);
        return;
    }
}

void VzhuhSolver::AddEdgeToThisQueue(int edge, int queue_index) {
    ++aux_counter4;

    if (edges[edge]->queue_index == queue_index) {
        return;
    }

    if (edges[edge]->queue_index >= 0) {
        edge_heaps[edges[edge]->queue_index]->Remove(edges[edge].get());
    }
    edge_heaps[queue_index]->Insert(edges[edge].get());
    edges[edge]->queue_index = queue_index;
}

void VzhuhSolver::RemoveEdgeFromQueue(int edge) {
    if (edges[edge]->queue_index >= 0) {
        edge_heaps[edges[edge]->queue_index]->Remove(edges[edge].get());
        edges[edge]->queue_index = -1;
    }
}

void VzhuhSolver::AddNodeToQueue(int node, int queue_index) {
    ++aux_counter4;

    if (nodes[node]->queue_index == queue_index) {
        return;
    }

    if (nodes[node]->queue_index >= 0) {
        node_heaps[nodes[node]->queue_index]->Remove(nodes[node].get());
    }
    node_heaps[queue_index]->Insert(nodes[node].get());
    nodes[node]->queue_index = queue_index;
}

void VzhuhSolver::RemoveNodeFromQueue(int node) {
    if (nodes[node]->queue_index >= 0) {
        node_heaps[nodes[node]->queue_index]->Remove(nodes[node].get());
        nodes[node]->queue_index = -1;
    }
}

void VzhuhSolver::AddPQPlusPlus(int first, int second, int edge) {
    int queue_index = TreeTreeQueueIndex(second, &trees[first]->pq_plus_plus);
    if (queue_index > 0) {
        AddEdgeToThisQueue(edge, queue_index);
    } else {
        edge_heaps.emplace_back(std::make_unique<PairingHeap<Edge> >());
        queue_index = static_cast<int>(edge_heaps.size()) - 1;
        AddEdgeToThisQueue(edge, queue_index);

        trees[first]->pq_plus_plus.emplace_back(second, queue_index);
        trees[second]->pq_plus_plus.emplace_back(first, queue_index);
    }
}

void VzhuhSolver::AddPQPlusMinus(int tree_plus, int tree_minus, int edge) {
    int queue_index = TreeTreeQueueIndex(tree_minus, &trees[tree_plus]->pq_plus_minus);
    if (queue_index > 0) {
        AddEdgeToThisQueue(edge, queue_index);
    } else {
        edge_heaps.emplace_back(std::make_unique<PairingHeap<Edge> >());
        queue_index = static_cast<int>(edge_heaps.size()) - 1;
        AddEdgeToThisQueue(edge, queue_index);

        trees[tree_plus]->pq_plus_minus.emplace_back(tree_minus, queue_index);
        trees[tree_minus]->pq_minus_plus.emplace_back(tree_plus, queue_index);
    }
}

int VzhuhSolver::TreeTreeQueueIndex(int other_tree,
                                    std::vector<std::pair<int, int> > *tree_neighbors) const {
    // may delete dead trees from the tree_neighbors
    for (int i = 0; i < static_cast<int>(tree_neighbors->size()); ++i) {
        if ((*tree_neighbors)[i].first == other_tree) {
            return (*tree_neighbors)[i].second;
        }
        if (!trees[(*tree_neighbors)[i].first]->is_alive) {
            (*tree_neighbors)[i] = tree_neighbors->back();
            tree_neighbors->pop_back();
            --i;
        }
    }
    return -1;
}

int VzhuhSolver::MinPlusEmptyEdge(int queue_index) const {
    if (edge_heaps[queue_index]->GetMin()) {
        Edge *edge = edge_heaps[queue_index]->GetMin();
        return edge->index;
    }
    return -1;
}

int VzhuhSolver::MinPlusPlusInternalEdge(int queue_index) {
    while (edge_heaps[queue_index]->GetMin()) {
        Edge *edge = edge_heaps[queue_index]->GetMin();

        int head = Head(edge->index);
        int tail = Tail(edge->index);

        // if (!nodes[head]->plus || !nodes[tail]->plus || (nodes[head]->tree != nodes[tail]->tree)) {
        //     throw std::runtime_error("MinPlusPlusInternalEdge: incorrect top");
        // }

        if (head != tail) {
            return edge->index;
        }
        edge_heaps[queue_index]->RemoveMin();
        edges[edge->index]->queue_index = -1;
    }
    return -1;
}

int VzhuhSolver::MinPlusPlusExternalEdge(int queue_index) const {
    if (edge_heaps[queue_index]->GetMin()) {
        Edge *edge = edge_heaps[queue_index]->GetMin();
        return edge->index;

        // int head = Head(edge);
        // int tail = Tail(edge);
        // if (nodes[head]->plus && nodes[tail]->plus && (nodes[head]->tree != nodes[tail]->tree)) {
        //     return edge;
        // }
        // throw std::runtime_error("MinPlusPlusExternalEdge: incorrect top");
    }
    return -1;
}

int VzhuhSolver::MinPlusMinusExternalEdge(int queue_index) const {
    if (edge_heaps[queue_index]->GetMin()) {
        Edge *edge = edge_heaps[queue_index]->GetMin();
        return edge->index;

        // int head = Head(edge); // plus
        // int tail = Tail(edge); // minus
        // if (!nodes[head]->plus) {
        //     std::swap(head, tail);
        // }
        // if (nodes[head]->plus && !nodes[tail]->plus && (nodes[head]->tree != nodes[tail]->tree) &&
        //     nodes[tail]->tree) {
        //     return edge;
        // }
        // throw std::runtime_error("MinPlusMinusExternalEdge: incorrect top");
    }
    return -1;
}

int VzhuhSolver::MinMinusBlossom(int queue_index) const {
    if (node_heaps[queue_index]->GetMin()) {
        Node *node = node_heaps[queue_index]->GetMin();
        // if (!nodes[node]->is_alive || !nodes[node]->tree || nodes[node]->plus || nodes[node].
        //     blossom_parent) {
        //     throw std::runtime_error("In MinMinusBlossom: incorrect top");
        // }
        return node->index;
    }
    return -1;
}

int VzhuhSolver::PopExpandableBlossom(int tree) const {
    int node = MinMinusBlossom(trees[tree]->minus_blossoms);
    if (node >= 0) {
        if (DualVariableQuadrupled(node) == 0) {
            node_heaps[trees[tree]->minus_blossoms]->RemoveMin();
            return node;
        }
    }
    return -1;
}

int VzhuhSolver::PlusEmptySlack(int tree) {
    int edge = MinPlusEmptyEdge(trees[tree]->plus_empty_edges);
    if (edge >= 0) {
        int slack = SlackQuadrupled(edge);
        // if (slack < 0) {
        //     int head = Head(edge);
        //     int tail = Tail(edge);
        //
        //     std::cout << "head: " << head << " plus: " << nodes[head]->plus << " tree: " << nodes[head]->tree <<
        //         std::endl;
        //     std::cout << "tail: " << tail << " plus: " << nodes[tail]->plus << " tree: " << nodes[tail]->tree <<
        //         std::endl;
        //
        //     throw std::runtime_error("PlusEmptySlack: negative slack");
        // }
        return slack;
    }
    return INT32_MAX;
}

int VzhuhSolver::PlusPlusInternalSlack(int tree) {
    int edge = MinPlusPlusInternalEdge(trees[tree]->plus_plus_internal_edges);
    if (edge >= 0) {
        return SlackQuadrupled(edge);
    }
    return INT32_MAX;
}

int VzhuhSolver::MinMinusBlossomVariable(int tree) const {
    int node = MinMinusBlossom(trees[tree]->minus_blossoms);
    if (node >= 0) {
        return DualVariableQuadrupled(node);
    }
    return INT32_MAX;
}

std::vector<std::pair<int, int> > VzhuhSolver::PlusPlusExternalSlacks(int tree) {
    std::vector<std::pair<int, int> > result;
    result.reserve(trees[tree]->pq_plus_plus.size());

    for (int i = 0; i < static_cast<int>(trees[tree]->pq_plus_plus.size()); ++i) {
        auto [tree_neighbor, queue_index] = trees[tree]->pq_plus_plus[i];
        if (trees[tree_neighbor]->is_alive) {
            int edge = MinPlusPlusExternalEdge(queue_index);
            if (edge >= 0) {
                result.emplace_back(tree_neighbor, SlackQuadrupled(edge_heaps[queue_index]->GetMin()->index));
            }
        } else {
            trees[tree]->pq_plus_plus[i] = trees[tree]->pq_plus_plus.back();
            trees[tree]->pq_plus_plus.pop_back();
            --i;
        }
    }
    return result;
}

std::vector<std::pair<int, int> > VzhuhSolver::PlusMinusExternalSlacks(int tree) {
    // TODO avoid code duplication

    std::vector<std::pair<int, int> > result;
    result.reserve(trees[tree]->pq_plus_minus.size());

    for (int i = 0; i < static_cast<int>(trees[tree]->pq_plus_minus.size()); ++i) {
        auto [tree_neighbor, queue_index] = trees[tree]->pq_plus_minus[i];
        if (trees[tree_neighbor]->is_alive) {
            int edge = MinPlusMinusExternalEdge(queue_index);
            if (edge >= 0) {
                result.emplace_back(tree_neighbor, SlackQuadrupled(edge_heaps[queue_index]->GetMin()->index));
            }
        } else {
            trees[tree]->pq_plus_minus[i] = trees[tree]->pq_plus_minus.back();
            trees[tree]->pq_plus_minus.pop_back();
            --i;
        }
    }
    return result;
}

void VzhuhSolver::AddNodeToRecord(int node, PrimalUpdateRecord *record) {
    if (nodes[node]->is_in_record) {
        return;
    }
    record->changed_sign.push_back(node);
    nodes[node]->is_in_record = true;
}

int VzhuhSolver::InitNumVertices(const std::vector<std::tuple<int, int, int> > &edge_list_) {
    int n = 0;
    for (const std::tuple<int, int, int> &edge : edge_list_) {
        if (std::get<0>(edge) > n) {
            n = std::get<0>(edge);
        }
        if (std::get<1>(edge) > n) {
            n = std::get<1>(edge);
        }
    }
    ++n;

    return n;
}
