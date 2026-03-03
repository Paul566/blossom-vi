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
    blossom_parents.reserve(2 * num_vertices_elementary);
    for (int i = 0; i < num_vertices_elementary; ++i) {
        nodes.emplace_back(Node(i));
        blossom_parents.push_back(-1);
    }

    // reserve for adj list
    std::vector<int> degrees(num_vertices_elementary, 0);
    for (auto [from, to, weight] : edge_list_) {
        ++degrees[from];
        ++degrees[to];
    }
    for (int i = 0; i < num_vertices_elementary; ++i) {
        nodes[i].neighbors.reserve(degrees[i]);
    }
    
    edges.reserve(edge_list_.size());
    edge_weights.reserve(edge_list_.size());
    for (int i = 0; i < static_cast<int>(edge_list_.size()); ++i) {
        int head = std::get<0>(edge_list_[i]);
        int tail = std::get<1>(edge_list_[i]);

        edges.emplace_back(Edge(head, tail, std::get<2>(edge_list_[i])));
        edge_weights.push_back(std::get<2>(edge_list_[i]));

        nodes[head].neighbors.emplace_back(i * 2);
        nodes[tail].neighbors.emplace_back(i * 2 + 1);
    }

    primal_update_record.reserve(num_vertices_elementary);

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
            std::cout << "sample tree sizes: " << trees[alive_trees[0]].tree_nodes.size() << " " << trees[alive_trees[
                    alive_trees.size() / 2]].tree_nodes.size() << " " << trees[alive_trees.back()].tree_nodes.size()
                <<
                std::endl;

            int min_children = INT32_MAX;
            for (Node &node : nodes) {
                if (!node.blossom_children.empty() && node.blossom_children.size() < min_children) {
                    min_children = node.blossom_children.size();
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

VzhuhSolver::Edge::Edge(int head_, int tail_, int weight_) : queue_index(-1), heap_child(-1),
                                                             heap_next(-1),
                                                             heap_prev(-1),
                                                             head(head_),
                                                             tail(tail_),
                                                             elementary_head(head_),
                                                             elementary_tail(tail_), last_round_updated(-1) {
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

VzhuhSolver::Node::Node(int index_) : queue_index(-1), heap_child(-1), heap_next(-1), heap_prev(-1),
                                      old_blossom_parent(-1),
                                      matched_edge(-1), minus_parent(-1), receptacle_(index_), tree(-1),
                                      old_tree(-1), tree_var_at_birth(0),
                                      slack_diff(0),
                                      is_in_record(false) {
    is_alive = true;
    dual_var_quadrupled_amortized_ = 0;
    plus = false;
    old_plus = false;
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
    std::cout << "Adjacency list (to, matched):" << std::endl;

    for (int i = 0; i < num_vertices_elementary; ++i) {
        std::cout << i << ": ";
        for (ArcIndex arc : nodes[i].neighbors) {
            std::cout << "(" << OtherElementaryEnd(arc) << " "
                << edges[arc.index / 2].matched << ") ";
        }
        std::cout << std::endl;
    }

    std::cout << "Node structure: " << std::endl;
    for (int i(0); i < nodes.size(); ++i) {
        if (nodes[i].is_alive) {
            PrintNode(i);
        }
    }

    std::cout << "Tree structure: " << std::endl;
    for (int tree : alive_trees) {
        if (!trees[tree].is_alive) {
            continue;
        }
        std::cout << "root: " << trees[tree].root << " var: " << trees[tree].dual_var_quadrupled / 4.
            << std::endl;
    }
}

void VzhuhSolver::PrintNode(int node) const {
    std::cout << node << ": y_v = " << DualVariableQuadrupled(node) / 4.;
    if (blossom_parents[node] >= 0) {
        std::cout << " blossom_parent: " << blossom_parents[node];
    }
    std::cout << std::endl;
}

void VzhuhSolver::GreedyInit() {
    // first, make all the slacks non-negative
    for (int i = 0; i < num_vertices_elementary; ++i) {
        if (nodes[i].neighbors.empty()) {
            std::cout << "vertex " << i << ": no neighbors" << std::endl;
            throw std::runtime_error("Found an isolated vertex => no perfect matching exists");
        }

        int min_weight = edge_weights[nodes[i].neighbors.front().index / 2];
        for (ArcIndex arc : nodes[i].neighbors) {
            if (edge_weights[arc.index / 2] < min_weight) {
                min_weight = edge_weights[arc.index / 2];
            }
        }

        nodes[i].dual_var_quadrupled_amortized_ += 2 * min_weight;
        for (ArcIndex arc : nodes[i].neighbors) {
            edges[arc.index / 2].slack_quadrupled_amortized_ -= 2 * min_weight;
        }
    }

    for (int i = 0; i < num_vertices_elementary; ++i) {
        if (nodes[i].matched_edge.index >= 0) {
            continue;
        }

        ArcIndex smallest_slack_arc = nodes[i].neighbors.front();
        for (ArcIndex arc : nodes[i].neighbors) {
            if (edges[arc.index / 2].slack_quadrupled_amortized_ < edges[smallest_slack_arc.index / 2].slack_quadrupled_amortized_) {
                smallest_slack_arc = arc;
            }
        }

        int diff = edges[smallest_slack_arc.index / 2].slack_quadrupled_amortized_;
        nodes[i].dual_var_quadrupled_amortized_ += diff;
        for (ArcIndex arc : nodes[i].neighbors) {
            edges[arc.index / 2].slack_quadrupled_amortized_ -= diff;
        }
        if (nodes[OtherElementaryEnd(smallest_slack_arc)].matched_edge.index < 0) {
            // if the other vertex is also unmatched, match the edge
            edges[smallest_slack_arc.index / 2].matched = true;
            nodes[i].matched_edge = smallest_slack_arc;
            nodes[OtherElementaryEnd(smallest_slack_arc)].matched_edge = ReverseArc(smallest_slack_arc);
        }
    }
}

void VzhuhSolver::InitializeTrees() {
    // is called after GreedyInit
    // also initializes queues and actionable_edges

    // TODO no need to populate queues in tree initialization

    std::vector<int> roots;
    for (int root_index(0); root_index < num_vertices_elementary; ++root_index) {
        if (nodes[root_index].matched_edge.index >= 0) {
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
        trees.emplace_back(Tree(roots[i], i, 2 * i, 2 * i + 1));
        alive_trees.emplace_back(i);

        node_heaps.emplace_back(NodeHeap());
        edge_heaps.emplace_back(EdgeHeap());
        edge_heaps.emplace_back(EdgeHeap());

        nodes[roots[i]].tree = i;
        nodes[roots[i]].plus = true;
        nodes[roots[i]].old_tree = i;
        nodes[roots[i]].old_plus = true;
    }

    // initialize queues and actionable_edges
    for (int root_index : roots) {
        for (ArcIndex arc : nodes[root_index].neighbors) {
            int edge_to_neighbor = arc.index / 2;
            int neighbor = OtherEnd(arc);

            if (nodes[neighbor].tree < 0) {
                RemoveEdgeFromQueue(edge_to_neighbor);
                AddEdgeToThisQueue(edge_to_neighbor, trees[nodes[root_index].tree].plus_empty_edges);
            } else {
                // neighbor is another root
                int queue_index = TreeTreeQueueIndex(nodes[neighbor].tree, &trees[nodes[root_index].tree].pq_plus_plus);
                if (queue_index >= 0) {
                    RemoveEdgeFromQueue(edge_to_neighbor);
                    AddEdgeToThisQueue(edge_to_neighbor, queue_index);
                } else {
                    edge_heaps.emplace_back(EdgeHeap());
                    queue_index = static_cast<int>(edge_heaps.size()) - 1;
                    RemoveEdgeFromQueue(edge_to_neighbor);
                    AddEdgeToThisQueue(edge_to_neighbor, queue_index);

                    trees[nodes[root_index].tree].pq_plus_plus.emplace_back(nodes[neighbor].tree, queue_index);
                    trees[nodes[neighbor].tree].pq_plus_plus.emplace_back(nodes[root_index].tree, queue_index);
                }
            }

            if (edges[edge_to_neighbor].slack_quadrupled_amortized_ == 0) {
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
        if (edges[i].matched) {
            matching.emplace_back(edges[i].elementary_head, edges[i].elementary_tail);
        }
    }
}

void VzhuhSolver::ComputePrimalObjective() {
    primal_objective = 0;
    for (int i(0); i < static_cast<int>(edges.size()); ++i) {
        if (edges[i].matched) {
            primal_objective += edge_weights[i];
        }
    }
}

void VzhuhSolver::DestroyBlossoms() {
    RestoreFinalEdgeEnds();

    for (int blossom(nodes.size() - 1); !IsElementary(blossom); --blossom) {
        if (!nodes[blossom].is_alive) {
            continue;
        }
        if (blossom_parents[blossom] >= 0) {
            throw std::runtime_error("In DestroyBlossoms: not going top down");
        }

        if (params.verbose) {
            std::cout << "destroying " << blossom << std::endl;
            PrintGraph();
        }

        int receptacle = FindFinalReceptacle(blossom);
        for (int child : nodes[blossom].blossom_children) {
            blossom_parents[child] = -1;
        }
        UpdateMatching(blossom, receptacle);
        nodes[receptacle].matched_edge = nodes[blossom].matched_edge;
        nodes[blossom].is_alive = false;
    }

    if (params.verbose) {
        std::cout << "destroyed all" << std::endl;
        PrintGraph();
    }
}

void VzhuhSolver::RestoreFinalEdgeEnds() {
    std::vector<int> depths(nodes.size(), 0);
    for (int node(nodes.size() - 1); node >= 0; --node) {
        if (!nodes[node].is_alive) {
            continue;
        }
        if (blossom_parents[node] >= 0) {
            depths[node] = depths[blossom_parents[node]] + 1;
        }
    }

    std::vector<int> edges_to_restore;
    edges_to_restore.reserve(nodes.size() * 2);
    for (int node(0); node < nodes.size(); ++node) {
        if (nodes[node].matched_edge.index >= 0) {
            edges_to_restore.push_back(nodes[node].matched_edge.index / 2);
        }
        if (nodes[node].minus_parent.index >= 0) {
            edges_to_restore.push_back(nodes[node].minus_parent.index / 2);
        }
    }

    for (int edge : edges_to_restore) {
        int head = edges[edge].elementary_head;
        int tail = edges[edge].elementary_tail;

        while (depths[head] > 0 || depths[tail] > 0) {
            if (depths[head] > depths[tail]) {
                head = blossom_parents[head];
            } else if (depths[tail] > depths[head]) {
                tail = blossom_parents[tail];
            } else if (blossom_parents[head] == blossom_parents[tail]) {
                break;
            } else {
                head = blossom_parents[head];
                tail = blossom_parents[tail];
            }
        }

        edges[edge].head = head;
        edges[edge].tail = tail;
    }
}

int VzhuhSolver::FindFinalReceptacle(int blossom) const {
    // TODO double work, use arcs?
    int new_receptacle(-1);
    int head = edges[nodes[blossom].matched_edge.index / 2].elementary_head;
    while (blossom_parents[head] >= 0) {
        if (blossom_parents[head] == blossom) {
            new_receptacle = head;
            break;
        }
        head = blossom_parents[head];
    }
    if (new_receptacle < 0) {
        int tail = edges[nodes[blossom].matched_edge.index / 2].elementary_tail;
        while (blossom_parents[tail] >= 0) {
            if (blossom_parents[tail] == blossom) {
                new_receptacle = tail;
                break;
            }
            tail = blossom_parents[tail];
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
        if (nodes[i].is_alive) {
            node_index_alive[i] = index;
            dual_certificate.emplace_back(index, DualVariableQuadrupled(i), -1);
            for (int child : nodes[i].blossom_children) {
                std::get<2>(dual_certificate[node_index_alive[child]]) = index;
            }
            ++index;
        }
    }
}

void VzhuhSolver::ComputeDualObjectiveQuadrupled() {
    dual_objective = 0;
    for (int i(0); i < static_cast<int>(nodes.size()); ++i) {
        if (nodes[i].is_alive) {
            dual_objective += DualVariableQuadrupled(i);
        }
    }
}

bool VzhuhSolver::MakePrimalUpdates() {
    bool action_taken = false;
    primal_update_record.clear();

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
            Expand(blossom);
            blossom = PopExpandableBlossom(tree);
        }
    }

    // TODO make an early return if augmented the last pair of trees
    // unweighted solver like phase
    while (!actionable_edges.empty()) {
        int edge = actionable_edges.front();
        actionable_edges.pop();
        MakePrimalUpdate(edge);
    }
    while (!actionable_nodes.empty()) {
        int node = actionable_nodes.front();
        actionable_nodes.pop();
        MakePrimalUpdateForNode(node);
    }

    if (!primal_update_record.empty()) {
        action_taken = true;
    }

    if (params.debug) {
        ValidateEvenOddPaths();
        ValidateArcs();
    }

    // update queues and amortized variables phase
    // UpdateQueues();
    UpdateQueuesRecordTraversal();
    if (params.debug) {
        ValidateQueues();
    }

    if (params.debug) {
        std::vector<int> new_variables = NodeVariables();
        std::vector<int> new_slacks = EdgeSlacks();
        for (int i = 0; i < static_cast<int>(variables.size()); ++i) {
            if (!nodes[i].is_alive) {
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
                std::cout << "edge " << edges[edge_idx].elementary_head << " " << edges[edge_idx].elementary_tail <<
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
    std::vector<std::vector<int> > future_blossoms = OrganizeBlossomChildren();
    for (std::vector<int> &children : future_blossoms) {
        // if (children.size() > 1) {
        Shrink(children);
        // }
    }

    return action_taken;
}

void VzhuhSolver::MakePrimalUpdate(int edge) {
    // TODO make it for arcs

    ArcIndex arc(edge * 2);
    int parent = Head(edge);
    int child = Tail(edge);

    if (parent == child) {
        return;
    }

    if (nodes[parent].tree < 0) {
        std::swap(parent, child);
        arc = ReverseArc(arc);
    }

    if (nodes[parent].tree < 0) {
        // edge not adjacent to any trees
        return;
    }

    if (nodes[parent].tree == nodes[child].tree) {
        if (nodes[parent].plus && nodes[child].plus) {
            MakeCherryBlossom(edge);
            return;
        }
    }

    if (nodes[child].tree < 0) {
        if (nodes[parent].plus) {
            Grow(parent, arc);
            return;
        }
    }

    if (nodes[parent].tree != nodes[child].tree) {
        if (nodes[parent].plus && nodes[child].plus) {
            Augment(edge);
        }
    }
}

void VzhuhSolver::MakePrimalUpdateForNode(int node) {
    if (!nodes[node].plus || !nodes[node].is_alive) {
        return;
    }

    int old_top = node;
    if (nodes[old_top].old_blossom_parent >= 0) {
        old_top = nodes[old_top].old_blossom_parent;
    }
    int delta_slack = 0;
    if (nodes[old_top].old_tree >= 0) {
        if (nodes[old_top].old_plus) {
            delta_slack = -trees[nodes[old_top].old_tree].dual_var_quadrupled;
        } else {
            delta_slack = trees[nodes[old_top].old_tree].dual_var_quadrupled;
        }
    }
    for (ArcIndex arc : NonLoopNeighbors(node)) {
        int edge = arc.index / 2;
        if (edges[edge].maybe_has_zero_slack) {
            // compute the old slack to compare it to 0
            int slack = edges[edge].slack_quadrupled_amortized_;

            int neighbor = OtherEnd(arc);
            if (nodes[neighbor].old_blossom_parent >= 0) {
                neighbor = nodes[neighbor].old_blossom_parent;
            }

            if (old_top == neighbor) {
                // used to be a loop, but the node got expanded
                slack -= -2 * nodes[old_top].tree_var_at_birth;
            } else {
                slack += delta_slack;
                if (nodes[neighbor].old_tree >= 0) {
                    if (nodes[neighbor].old_plus) {
                        slack -= trees[nodes[neighbor].old_tree].dual_var_quadrupled;
                    } else {
                        slack += trees[nodes[neighbor].old_tree].dual_var_quadrupled;
                    }
                }
            }

            if (slack == 0) {
                // TODO use the information about head/tail that we know here
                MakePrimalUpdate(edge);
            } else {
                edges[edge].maybe_has_zero_slack = false;
            }
        }
    }
}

void VzhuhSolver::Expand(int blossom) {
    if (params.verbose) {
        std::cout << "EXPAND " << blossom << std::endl;
    }
    if (nodes[blossom].old_blossom_parent >= 0) {
        throw std::runtime_error("Expand: two-level expanding");
    }

    nodes[blossom].is_alive = false;

    RestoreEdgeEndsBeforeExpand(blossom);
    // now OtherEnd is safe

    // TODO use arc
    int new_receptacle = edges[nodes[blossom].matched_edge.index / 2].head;
    if (nodes[new_receptacle].old_blossom_parent != blossom) {
        new_receptacle = edges[nodes[blossom].matched_edge.index / 2].tail;
    }
    int old_receptacle = Receptacle(new_receptacle);
    int elder_child = edges[nodes[blossom].minus_parent.index / 2].head;
    if (nodes[elder_child].old_blossom_parent != blossom) {
        elder_child = edges[nodes[blossom].minus_parent.index / 2].tail;
    }

    UpdateInternalStructure(blossom, old_receptacle, new_receptacle, elder_child);

    for (int child : nodes[blossom].blossom_children) {
        AddNodeToRecord(child);
        actionable_nodes.push(child);
        if (nodes[child].tree >= 0) {
            trees[nodes[child].tree].tree_nodes.push_back(child);
        }
    }
}

void VzhuhSolver::RestoreEdgeEndsBeforeExpand(int blossom) {
    // makes head and tail of the adjacent edges to point to the children of the blossom,
    // make children top blossoms

    // TODO make better

    for (int child : nodes[blossom].blossom_children) {
        std::vector<int> elementary_descendants = ElementaryBlossomDescendants(child);
        for (int descendant : elementary_descendants) {
            for (ArcIndex arc : nodes[descendant].neighbors) {
                // TODO use arcs
                int edge = arc.index / 2;
                if (edges[edge].elementary_head == descendant) {
                    edges[edge].head = child;
                } else {
                    edges[edge].tail = child;
                }
            }
        }
    }

    for (int child : nodes[blossom].blossom_children) {
        blossom_parents[child] = -1;
    }
}

void VzhuhSolver::RotateReceptacle(int blossom, int new_receptacle) {
    // changes matched edges and minus_parents

    if (params.verbose) {
        std::cout << "ROTATE RECEPTACLE " << blossom << " " << new_receptacle << std::endl;
    }

    int old_receptacle = Receptacle(new_receptacle);
    if (old_receptacle == new_receptacle) {
        nodes[new_receptacle].matched_edge = nodes[blossom].matched_edge;
        return;
    }

    std::vector<ArcIndex> even_path = EvenPathToReceptacle(new_receptacle);
    std::vector<ArcIndex> odd_path = OddPathToReceptacle(new_receptacle);

    if (even_path.back().index == odd_path.back().index) {
        int i = 0;
        int common_node = old_receptacle;
        while (even_path[even_path.size() - 1 - i].index == odd_path[odd_path.size() - 1 - i].index) {
            common_node = ThisEnd(even_path[even_path.size() - 1 - i]);
            ++i;
        }

        RotateReceptacle(blossom, common_node);
        RotateReceptacle(blossom, new_receptacle);
        return;
    }

    int cur_node = new_receptacle;
    bool match = false;
    for (ArcIndex arc : even_path) {
        int edge = arc.index / 2;
        if (match) {
            MakeEdgeMatched(edge);
            cur_node = OtherEnd(arc);
        } else {
            MakeEdgeUnmatched(edge);
            nodes[cur_node].minus_parent = arc;
            cur_node = OtherEnd(arc);
            nodes[cur_node].minus_parent = ReverseArc(arc);
        }
        match = !match;
    }

    cur_node = new_receptacle;
    match = false;
    for (ArcIndex arc : odd_path) {
        if (match) {
            cur_node = OtherEnd(arc);
        } else {
            nodes[cur_node].minus_parent = arc;
            cur_node = OtherEnd(arc);
            nodes[cur_node].minus_parent = ReverseArc(arc);
        }
        match = !match;
    }

    nodes[new_receptacle].minus_parent = ArcIndex(-1);
    nodes[new_receptacle].matched_edge = nodes[blossom].matched_edge;
    nodes[new_receptacle].receptacle_ = new_receptacle;
    nodes[old_receptacle].receptacle_ = new_receptacle;
}

void VzhuhSolver::UpdateMatching(int blossom, int new_receptacle) {
    std::vector<ArcIndex> path = EvenPathToReceptacle(new_receptacle);
    bool match = false;
    for (ArcIndex arc : path) {
        int edge = arc.index / 2;
        if (match) {
            MakeEdgeMatched(edge);
        } else {
            MakeEdgeUnmatched(edge);
        }
        match = !match;
    }
    nodes[new_receptacle].matched_edge = nodes[blossom].matched_edge;
}

void VzhuhSolver::UpdateInternalStructure(int blossom,
                                          int old_receptacle,
                                          int new_receptacle,
                                          int elder_child) {
    // after expand, the remaining part of the blossom is a path in the tree

    RotateReceptacle(blossom, new_receptacle);
    std::vector<ArcIndex> path = EvenPathToReceptacle(elder_child);

    // the path that stays in the tree
    int cur_node = elder_child;
    nodes[elder_child].minus_parent = nodes[blossom].minus_parent;
    nodes[elder_child].plus = false;
    nodes[elder_child].tree = nodes[blossom].tree;
    for (ArcIndex arc : path) {
        cur_node = OtherEnd(arc);

        if (edges[arc.index / 2].matched) {
            nodes[cur_node].minus_parent = ArcIndex(-1);
            nodes[cur_node].plus = true;
        } else {
            nodes[cur_node].minus_parent = ReverseArc(arc);
            nodes[cur_node].plus = false;
        }
        nodes[cur_node].tree = nodes[blossom].tree;
    }

    // mark the nodes that stay + insert them into the record
    ++nodes_label_cnt;
    cur_node = elder_child;
    for (ArcIndex arc : path) {
        nodes[cur_node].label = nodes_label_cnt;
        AddNodeToRecord(cur_node);
        cur_node = OtherEnd(arc);
    }
    nodes[new_receptacle].label = nodes_label_cnt;
    AddNodeToRecord(new_receptacle);

    // clear the part that goes to waste
    for (int child : nodes[blossom].blossom_children) {
        if (nodes[child].label != nodes_label_cnt) {
            nodes[child].plus = false;
            nodes[child].minus_parent = ArcIndex(-1);
            nodes[child].tree = -1;
        }
    }

    // reset receptacles
    for (int child : nodes[blossom].blossom_children) {
        nodes[child].receptacle_ = child;
    }
}

std::vector<VzhuhSolver::ArcIndex> VzhuhSolver::EvenPathToReceptacle(int node) {
    int receptacle = Receptacle(node);
    std::vector<ArcIndex> path;

    while (node != receptacle) {
        int parent = OtherEnd(nodes[node].matched_edge);
        int grandparent = OtherEnd(nodes[parent].minus_parent);

        if (params.verbose) {
            std::cout << "EvenPathToReceptacle, parent/grandparent: " << parent << " " << grandparent << std::endl;
        }

        path.push_back(nodes[node].matched_edge);
        path.push_back(nodes[parent].minus_parent);

        node = grandparent;

        if (path.size() > num_vertices_elementary) {
            std::cout << "receptacle: " << receptacle << ", node: " << node << std::endl;
            throw std::runtime_error("EvenPathToReceptacle: infinite loop");
        }
    }

    return path;
}

std::vector<VzhuhSolver::ArcIndex> VzhuhSolver::OddPathToReceptacle(int node) {
    int receptacle = Receptacle(node);
    if (receptacle == node) {
        throw std::runtime_error("OddPathToReceptacle: node is receptacle");
    }

    std::vector<ArcIndex> path = {nodes[node].minus_parent};
    node = OtherEnd(nodes[node].minus_parent);
    while (node != receptacle) {
        int parent = OtherEnd(nodes[node].matched_edge);
        int grandparent = OtherEnd(nodes[parent].minus_parent);

        path.push_back(nodes[node].matched_edge);
        path.push_back(nodes[parent].minus_parent);

        node = grandparent;
    }

    return path;
}

void VzhuhSolver::ExpandChildBeforeGrow(int blossom) {
    if (nodes[blossom].old_blossom_parent >= 0) {
        // avoid two-level expansion
        return;
    }

    nodes[blossom].is_alive = false;

    RestoreEdgeEndsBeforeExpand(blossom);
    // now OtherEnd, Head, Tail is safe

    // TODO use arcs
    int new_receptacle = edges[nodes[blossom].matched_edge.index / 2].head;
    if (nodes[new_receptacle].old_blossom_parent != blossom) {
        new_receptacle = edges[nodes[blossom].matched_edge.index / 2].tail;
    }
    UpdateMatching(blossom, new_receptacle);

    for (int child : nodes[blossom].blossom_children) {
        nodes[child].receptacle_ = child;
        nodes[child].plus = false;
        nodes[child].minus_parent = ArcIndex(-1);
        nodes[child].tree = -1;
    }

    for (int child : nodes[blossom].blossom_children) {
        AddNodeToRecord(child);
        if (nodes[child].tree >= 0) {
            actionable_nodes.push(child);
        }
    }
}

void VzhuhSolver::Grow(int parent, ArcIndex arc) {
    int child = OtherEnd(arc);
    int tree = nodes[parent].tree;

    if (!nodes[parent].plus) {
        throw std::runtime_error("In Grow: parent is not a plus");
    }
    if (nodes[child].tree >= 0) {
        throw std::runtime_error("In Grow: child vertex is not free");
    }
    if (nodes[child].matched_edge.index < 0) {
        std::cout << "parent " << parent << ", child " << child << std::endl;
        throw std::runtime_error("Grow: child has no matched_edge");
    }

    int grandchild = OtherEnd(nodes[child].matched_edge);
    if (params.verbose) {
        std::cout << "GROW " << parent << " " << child << " " << grandchild << std::endl;
    }

    if (!IsElementary(child)) {
        if (DualVariableQuadrupled(child,
                                   nodes[child].old_tree,
                                   nodes[child].old_plus,
                                   nodes[child].old_blossom_parent) == 0) {
            if (params.verbose) {
                std::cout << "EXPAND CHILD, new child/grandchild: ";
            }
            ExpandChildBeforeGrow(child);
            child = OtherEnd(arc);
            grandchild = OtherEnd(nodes[child].matched_edge);
            if (params.verbose) {
                std::cout << child << " " << grandchild << std::endl;
            }
        }
    }

    nodes[child].minus_parent = ReverseArc(arc);
    nodes[child].tree = tree;
    nodes[child].plus = false;

    nodes[grandchild].tree = tree;
    nodes[grandchild].plus = true;

    AddNodeToRecord(child);
    AddNodeToRecord(grandchild);

    trees[tree].tree_nodes.push_back(child);
    trees[tree].tree_nodes.push_back(grandchild);

    actionable_nodes.push(grandchild);
}

void VzhuhSolver::MakeCherryBlossom(int edge_plus_plus) {
    int head = Head(edge_plus_plus);
    int tail = Tail(edge_plus_plus);

    if (Receptacle(head) == Receptacle(tail)) {
        return;
    }

    if (params.verbose) {
        std::cout << "MAKE CHERRY BLOSSOM " << head << " " << tail << std::endl;
    }

    auto [first_bound, second_bound] = CherryPathBounds(head, tail);

    UpdateCherryPath(head, first_bound);
    UpdateCherryPath(tail, second_bound);

    if (head != first_bound) {
        nodes[head].minus_parent = ArcIndex(edge_plus_plus * 2);
    }
    if (tail != second_bound) {
        nodes[tail].minus_parent = ArcIndex(edge_plus_plus * 2 + 1);
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
        first_vertex = OtherEnd(nodes[first_vertex].matched_edge);
        first_vertex = OtherEnd(nodes[first_vertex].minus_parent);
    }
    while (second_vertex != lca) {
        if (Receptacle(second_vertex) == lca_receptacle) {
            break;
        }
        second_vertex = OtherEnd(nodes[second_vertex].matched_edge);
        second_vertex = OtherEnd(nodes[second_vertex].minus_parent);
    }

    return {first_vertex, second_vertex};
}

void VzhuhSolver::UpdateCherryPath(int lower_node, int upper_node) {
    int receptacle = Receptacle(upper_node);
    nodes[Receptacle(upper_node)].receptacle_ = receptacle; // TODO does this line do anything?

    AddNodeToRecord(lower_node);
    while (lower_node != upper_node) {
        int parent = OtherEnd(nodes[lower_node].matched_edge);
        int grandparent = OtherEnd(nodes[parent].minus_parent);

        AddNodeToRecord(parent);
        AddNodeToRecord(grandparent);

        nodes[Receptacle(lower_node)].receptacle_ = receptacle;
        nodes[Receptacle(parent)].receptacle_ = receptacle;

        if (!nodes[parent].plus) {
            nodes[parent].plus = true;
            actionable_nodes.push(parent);
        }

        if (grandparent != upper_node) {
            nodes[grandparent].minus_parent = ReverseArc(nodes[parent].minus_parent);
        }

        lower_node = grandparent;
    }
}

void VzhuhSolver::Augment(int edge_plus_plus) {
    num_trees_alive -= 2;

    int head = Head(edge_plus_plus);
    int tail = Tail(edge_plus_plus);

    if (params.verbose) {
        std::cout << "AUGMENT " << head << " " << tail << ", elementary ends: " << edges[edge_plus_plus].
            elementary_head
            << " " << edges[edge_plus_plus].elementary_tail <<
            ", is in zero slack set: " << edges[edge_plus_plus].maybe_has_zero_slack << std::endl;
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

    int first_tree = nodes[head].tree;
    int second_tree = nodes[tail].tree;
    AugmentPath(path);
    ClearTree(first_tree);
    ClearTree(second_tree);
}

std::vector<int> VzhuhSolver::PathToRoot(int node_plus) {
    int root = TopBlossom(trees[nodes[node_plus].tree].root);

    std::vector<int> path;
    while (node_plus != root) {
        path.push_back(nodes[node_plus].matched_edge.index / 2);
        node_plus = OtherEnd(nodes[node_plus].matched_edge);
        path.push_back(nodes[node_plus].minus_parent.index / 2);
        node_plus = OtherEnd(nodes[node_plus].minus_parent);
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

void VzhuhSolver::ClearTree(int tree) {
    trees[tree].is_alive = false;

    // TODO have a special routine for the last two trees
    // TODO amortize

    for (int node : trees[tree].tree_nodes) {
        if (blossom_parents[node] >= 0 || nodes[node].tree != tree) {
            continue;
        }

        // TODO make better

        AddNodeToRecord(node);
        nodes[node].tree = -1;
        nodes[node].minus_parent = ArcIndex(-1);
        nodes[node].receptacle_ = node;
        nodes[node].plus = false;
    }
}

void VzhuhSolver::UpdateQueuesRecordTraversal() {
    for (int node : primal_update_record) {
        int old_tree = nodes[node].old_tree;
        int old_blossom_parent = nodes[node].old_blossom_parent;
        bool old_plus = nodes[node].old_plus;
        bool plus = nodes[node].plus;
        int tree = nodes[node].tree;

        if (!nodes[node].is_alive || (old_blossom_parent < 0 &&
            old_plus == plus &&
            old_tree == tree)) {
            nodes[node].is_in_record = false;
            continue;
        }

        // update dual_var_quadrupled_amortized_
        if (old_tree != tree || old_plus != plus || old_blossom_parent >= 0) {
            if (old_tree >= 0 && old_blossom_parent < 0) {
                if (old_plus) {
                    nodes[node].dual_var_quadrupled_amortized_ += trees[old_tree].dual_var_quadrupled;
                } else {
                    nodes[node].dual_var_quadrupled_amortized_ -= trees[old_tree].dual_var_quadrupled;
                }
            }
            if (tree >= 0) {
                if (plus) {
                    nodes[node].dual_var_quadrupled_amortized_ -= trees[tree].dual_var_quadrupled;
                } else {
                    nodes[node].dual_var_quadrupled_amortized_ += trees[tree].dual_var_quadrupled;
                }
            }
        }

        // update minus_blossoms queues
        if (!IsElementary(node)) {
            if (old_plus != plus || old_tree != tree) {
                if (old_tree >= 0 && !old_plus) {
                    // was a minus blossom
                    RemoveNodeFromQueue(node);
                }
                if (tree >= 0 && !plus) {
                    // became a minus blossom
                    AddNodeToQueue(node, trees[tree].minus_blossoms);
                }
            }
        }

        nodes[node].slack_diff = 0;
        if (tree >= 0) {
            if (plus) {
                nodes[node].slack_diff += trees[tree].dual_var_quadrupled;
            } else {
                nodes[node].slack_diff -= trees[tree].dual_var_quadrupled;
            }
        }
        int old_top_node = node;
        if (old_blossom_parent >= 0) {
            old_top_node = old_blossom_parent;
        }
        if (nodes[old_top_node].old_tree >= 0) {
            if (nodes[old_top_node].old_plus) {
                nodes[node].slack_diff -= trees[nodes[old_top_node].old_tree].dual_var_quadrupled;
            } else {
                nodes[node].slack_diff += trees[nodes[old_top_node].old_tree].dual_var_quadrupled;
            }
        }
    }

    for (int node : primal_update_record) {
        if (!nodes[node].is_alive || (nodes[node].old_blossom_parent < 0 &&
            nodes[node].old_plus == nodes[node].plus &&
            nodes[node].old_tree == nodes[node].tree)) {
            nodes[node].is_in_record = false;
            continue;
        }

        for (ArcIndex arc : NonLoopNeighbors(node)) {
            // TODO better logic here
            // compute the queue index
            int queue_index = -1;

            int other_end = OtherEnd(arc);

            if (nodes[node].tree == nodes[other_end].tree) {
                if (nodes[node].tree >= 0) {
                    if (nodes[node].plus && nodes[other_end].plus) {
                        queue_index = trees[nodes[node].tree].plus_plus_internal_edges;
                    }
                }
            } else {
                // different trees, make head always in some tree
                int head = node;
                int tail = other_end;
                if (nodes[head].tree < 0) {
                    std::swap(head, tail);
                }

                if (nodes[head].plus) {
                    if (nodes[tail].tree < 0) { // (+, 0)
                        queue_index = trees[nodes[head].tree].plus_empty_edges;
                    } else {
                        if (nodes[tail].plus) { // (+, +)
                            queue_index = TreeTreeQueueIndex(nodes[tail].tree, &trees[nodes[head].tree].pq_plus_plus);
                            if (queue_index < 0) {
                                edge_heaps.emplace_back(EdgeHeap());
                                queue_index = static_cast<int>(edge_heaps.size()) - 1;
                                trees[nodes[head].tree].pq_plus_plus.emplace_back(nodes[tail].tree, queue_index);
                                trees[nodes[tail].tree].pq_plus_plus.emplace_back(nodes[head].tree, queue_index);
                            }
                        } else {    // (+, -)
                            queue_index = TreeTreeQueueIndex(nodes[tail].tree, &trees[nodes[head].tree].pq_plus_minus);
                            if (queue_index < 0) {
                                edge_heaps.emplace_back(EdgeHeap());
                                queue_index = static_cast<int>(edge_heaps.size()) - 1;
                                trees[nodes[tail].tree].pq_minus_plus.emplace_back(nodes[head].tree, queue_index);
                                trees[nodes[head].tree].pq_plus_minus.emplace_back(nodes[tail].tree, queue_index);
                            }
                        }
                    }
                } else {
                    if (nodes[tail].plus) { // (-, +)
                        queue_index = TreeTreeQueueIndex(nodes[head].tree, &trees[nodes[tail].tree].pq_plus_minus);
                        if (queue_index < 0) {
                            edge_heaps.emplace_back(EdgeHeap());
                            queue_index = static_cast<int>(edge_heaps.size()) - 1;
                            trees[nodes[tail].tree].pq_plus_minus.emplace_back(nodes[head].tree, queue_index);
                            trees[nodes[head].tree].pq_minus_plus.emplace_back(nodes[tail].tree, queue_index);
                        }
                    }
                }
            }

            if (edges[arc.index / 2].last_round_updated < current_round) {
                UpdateEdgeInfo(arc.index / 2, node, other_end, queue_index);
            }
        }
    }

    for (int node : primal_update_record) {
        nodes[node].old_blossom_parent = blossom_parents[node];
        nodes[node].old_plus = nodes[node].plus;
        nodes[node].old_tree = nodes[node].tree;
        nodes[node].is_in_record = false;
        nodes[node].slack_diff = 0;
    }
}

void VzhuhSolver::UpdateEdgeInfo(int edge, int endpoint, int other_endpoint, int queue_index) {
    edges[edge].last_round_updated = current_round;

    if (nodes[endpoint].old_blossom_parent >= 0 && nodes[endpoint].old_blossom_parent == nodes[other_endpoint].
        old_blossom_parent) {
        // edge was a loop
        edges[edge].slack_quadrupled_amortized_ -= 2 * nodes[nodes[endpoint].old_blossom_parent].tree_var_at_birth;

        if (nodes[endpoint].tree >= 0) {
            if (nodes[endpoint].plus) {
                edges[edge].slack_quadrupled_amortized_ += trees[nodes[endpoint].tree].dual_var_quadrupled;
            } else {
                edges[edge].slack_quadrupled_amortized_ -= trees[nodes[endpoint].tree].dual_var_quadrupled;
            }
        }

        if (nodes[other_endpoint].tree >= 0) {
            if (nodes[other_endpoint].plus) {
                edges[edge].slack_quadrupled_amortized_ += trees[nodes[other_endpoint].tree].dual_var_quadrupled;
            } else {
                edges[edge].slack_quadrupled_amortized_ -= trees[nodes[other_endpoint].tree].dual_var_quadrupled;
            }
        }
    } else {
        edges[edge].slack_quadrupled_amortized_ += nodes[endpoint].slack_diff;
        edges[edge].slack_quadrupled_amortized_ += nodes[other_endpoint].slack_diff;
    }

    RemoveEdgeFromQueue(edge);
    if (Receptacle(endpoint) != Receptacle(other_endpoint) && queue_index >= 0) {
        AddEdgeToThisQueue(edge, queue_index);
    }
}

std::vector<std::vector<int> > VzhuhSolver::OrganizeBlossomChildren() {
    ++nodes_label_cnt;
    int label_zero = nodes_label_cnt;

    for (int node : primal_update_record) {
        if (nodes[node].tree < 0) {
            continue;
        }
        if (nodes[node].label >= label_zero) {
            // already seen this node
            continue;
        }

        int receptacle = Receptacle(node);
        if (nodes[receptacle].label < label_zero) {
            nodes[receptacle].label = nodes_label_cnt;
            ++nodes_label_cnt;
        }
        if (nodes[node].label < label_zero) {
            nodes[node].label = nodes[receptacle].label;
        }
    }

    std::vector<int> sizes(nodes_label_cnt - label_zero, 0);
    std::vector<int> label_diff_to_index(nodes_label_cnt - label_zero, -1);

    for (int node : primal_update_record) {
        if (nodes[node].label < label_zero) {
            continue;
        }
        ++sizes[nodes[node].label - label_zero];
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

    for (int node : primal_update_record) {
        if (nodes[node].label < label_zero) {
            continue;
        }
        if (label_diff_to_index[nodes[node].label - label_zero] >= 0) {
            result[label_diff_to_index[nodes[node].label - label_zero]].push_back(node);
        }
    }

    return result;
}

void VzhuhSolver::Shrink(std::vector<int> &children) {
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

    nodes.emplace_back(Node(new_index));
    blossom_parents.push_back(-1);
    nodes.back().plus = true;
    nodes.back().old_plus = true;
    nodes.back().blossom_children = std::move(children);
    nodes.back().matched_edge = nodes[receptacle].matched_edge;
    nodes.back().tree = nodes[receptacle].tree;
    nodes.back().old_tree = nodes.back().tree;
    nodes.back().dual_var_quadrupled_amortized_ = -trees[nodes.back().tree].dual_var_quadrupled;
    nodes.back().tree_var_at_birth = trees[nodes.back().tree].dual_var_quadrupled;

    // update blossom_parent of the children, set labels to empty, update variables
    for (int child : nodes.back().blossom_children) {
        if (blossom_parents[child] >= 0) {
            throw std::runtime_error("In Node: some child node already has a parent");
        }
        blossom_parents[child] = new_index;
        nodes[child].old_blossom_parent = new_index;
        nodes[child].dual_var_quadrupled_amortized_ += trees[nodes.back().tree].dual_var_quadrupled;
    }

    // add the new blossom to the list of vertices in the tree
    trees[nodes.back().tree].tree_nodes.push_back(new_index);
}

void VzhuhSolver::ValidateQueues() {
    // checks that the queues are in a correct state, throws if not

    // every queue must hold only the edges that belong to this queue

    // for (int tree : alive_trees) {
    //     if (!trees[tree].is_alive) {
    //         continue;
    //     }
    //
    //     // check plus empty
    //     edge_heaps[trees[tree].plus_empty_edges].ValidateHeap("plus empty");
    //     for (auto heap_node : edge_heaps[trees[tree].plus_empty_edges].heap_) {
    //         int edge = heap_node.value;
    //
    //         int first = Head(edge); // plus
    //         int second = Tail(edge); // empty
    //         if (nodes[first].tree != tree) {
    //             std::swap(first, second);
    //         }
    //         if (nodes[first].tree != tree || !nodes[first].plus) {
    //             throw std::runtime_error("Incorrect plus empty queue");
    //         }
    //         if (nodes[second].tree >= 0) {
    //             throw std::runtime_error("Incorrect plus empty queue");
    //         }
    //     }
    //
    //     // check (plus, plus) internal
    //     edge_heaps[trees[tree].plus_plus_internal_edges].ValidateHeap("+ + int");
    //     for (auto heap_node : edge_heaps[trees[tree].plus_plus_internal_edges].heap_) {
    //         int edge = heap_node.value;
    //
    //         int first = Head(edge);
    //         int second = Tail(edge);
    //
    //         if (nodes[first].tree != tree || !nodes[first].plus) {
    //             throw std::runtime_error("Incorrect plus plus internal queue");
    //         }
    //         if (nodes[second].tree != tree || !nodes[second].plus) {
    //             throw std::runtime_error("Incorrect plus plus internal queue");
    //         }
    //     }
    //
    //     // check minus blossoms
    //     node_heaps[trees[tree].minus_blossoms].ValidateHeap("minus blossoms");
    //     for (auto heap_node : node_heaps[trees[tree].minus_blossoms].heap_) {
    //         int node = heap_node.value;
    //
    //         if (IsElementary(node) || (nodes[node].tree != tree) || nodes[node].plus) {
    //             throw std::runtime_error("Incorrect minus blossom queue");
    //         }
    //     }
    //
    //     // check (plus, plus) external
    //     for (auto [other_tree, queue_index] : trees[tree].pq_plus_plus) {
    //         if (!trees[other_tree].is_alive) {
    //             continue;
    //         }
    //         for (auto heap_node : edge_heaps[queue_index].heap_) {
    //             int edge = heap_node.value;
    //
    //             int first = Head(edge); // in this tree
    //             int second = Tail(edge); // in the other tree
    //             if (nodes[first].tree != tree) {
    //                 std::swap(first, second);
    //             }
    //             if (nodes[first].tree != tree || !nodes[first].plus) {
    //                 throw std::runtime_error("Incorrect plus plus queue");
    //             }
    //             if (nodes[second].tree != other_tree || !nodes[second].plus) {
    //                 throw std::runtime_error("Incorrect plus plus queue");
    //             }
    //         }
    //         edge_heaps[queue_index].ValidateHeap("+ + external");
    //     }
    //
    //     // check plus minus external
    //     for (auto [other_tree, queue_index] : trees[tree].pq_plus_minus) {
    //         if (!trees[other_tree].is_alive) {
    //             continue;
    //         }
    //         edge_heaps[queue_index].ValidateHeap("+ - ext");
    //         for (auto heap_node : edge_heaps[queue_index].heap_) {
    //             int edge = heap_node.value;
    //
    //             int first = Head(edge); // in this tree
    //             int second = Tail(edge); // in the other tree
    //             if (nodes[first].tree != tree) {
    //                 std::swap(first, second);
    //             }
    //             if (nodes[first].tree != tree || !nodes[first].plus) {
    //                 throw std::runtime_error("Incorrect plus minus queue");
    //             }
    //             if (nodes[second].tree != other_tree || nodes[second].plus) {
    //                 throw std::runtime_error("Incorrect plus minus queue");
    //             }
    //         }
    //     }
    //
    //     // check minus plus external
    //     for (auto [other_tree, queue_index] : trees[tree].pq_minus_plus) {
    //         if (!trees[other_tree].is_alive) {
    //             continue;
    //         }
    //         edge_heaps[queue_index].ValidateHeap("- + ext");
    //         for (auto heap_node : edge_heaps[queue_index].heap_) {
    //             int edge = heap_node.value;
    //
    //             int first = Head(edge); // in this tree
    //             int second = Tail(edge); // in the other tree
    //             if (nodes[first].tree != tree) {
    //                 std::swap(first, second);
    //             }
    //             if (nodes[first].tree != tree || nodes[first].plus) {
    //                 throw std::runtime_error("Incorrect minus plus queue");
    //             }
    //             if (nodes[second].tree != other_tree || !nodes[second].plus) {
    //                 throw std::runtime_error("Incorrect minus plus queue");
    //             }
    //         }
    //     }
    // }
}

std::vector<int> VzhuhSolver::NodeVariables() const {
    std::vector<int> result;
    for (int i(0); i < static_cast<int>(nodes.size()); ++i) {
        if (nodes[i].is_alive) {
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
        if (!nodes[node].is_alive) {
            continue;
        }
        if (blossom_parents[node] >= 0) {
            depths[node] = depths[blossom_parents[node]] + 1;
        }
    }

    std::vector<int> result;
    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
        if (Head(i) != Tail(i)) {
            result.push_back(SlackQuadrupled(i));
        } else {
            int edge(i);

            int head = edges[edge].elementary_head;
            int tail = edges[edge].elementary_tail;

            while (head != tail) {
                if (depths[head] > depths[tail]) {
                    head = blossom_parents[head];
                } else {
                    tail = blossom_parents[tail];
                }
            }

            result.push_back(edges[edge].slack_quadrupled_amortized_ - 2 * nodes[head].tree_var_at_birth);
        }
    }

    return result;
}

void VzhuhSolver::ValidateEvenOddPaths() {
    for (int node(0); node < static_cast<int>(nodes.size()); ++node) {
        if (!nodes[node].is_alive || blossom_parents[node] >= 0) {
            continue;
        }

        int receptacle = Receptacle(node);

        auto even_path = EvenPathToReceptacle(node);
        if (even_path.size() % 2 != 0) {
            throw std::runtime_error("ValidateEvenOddPaths: even path is not even");
        }
        int cur_node = node;
        for (ArcIndex arc : even_path) {
            cur_node = OtherEnd(arc);
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
            for (ArcIndex arc : odd_path) {
                cur_node = OtherEnd(arc);
            }
            if (cur_node != receptacle) {
                throw std::runtime_error("ValidateEvenOddPaths: incorrect odd path");
            }
        }
    }
}

void VzhuhSolver::ValidateArcs() {
    for (int node = 0; node < static_cast<int>(nodes.size()); ++node) {
        if (blossom_parents[node] >= 0 || !nodes[node].is_alive) {
            continue;
        }

        if (nodes[node].minus_parent.index >= 0) {
            if (node != ThisEnd(nodes[node].minus_parent)) {
                throw std::runtime_error("ValidateArcs: incorrect minus parent");
            }
        }
        if (nodes[node].matched_edge.index >= 0) {
            if (node != ThisEnd(nodes[node].matched_edge)) {
                throw std::runtime_error("ValidateArcs: incorrect matched edge");
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
        trees[alive_trees[i]].dual_var_quadrupled += deltas[i];
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
        if (!trees[alive_trees[i]].is_alive) {
            alive_trees[i] = alive_trees.back();
            alive_trees.pop_back();
            --i;
        }
    }
}

std::vector<int> VzhuhSolver::VariableDeltas() {
    DualConstraints dual_constraints = GetDualConstraints();

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
    DualConstraints result(std::vector<int>(alive_trees.size(), INT32_MAX),
                           std::vector(alive_trees.size(), std::vector<std::pair<int, int> >()),
                           std::vector(alive_trees.size(), std::vector<std::pair<int, int> >()));

    // update alive indices
    for (int i = 0; i < static_cast<int>(alive_trees.size()); ++i) {
        trees[alive_trees[i]].alive_index = i;
    }

    for (int i = 0; i < static_cast<int>(alive_trees.size()); ++i) {
        int tree = alive_trees[i];

        // get upper_bound
        int plus_empty = PlusEmptySlack(tree);
        int plus_plus_internal = PlusPlusInternalSlack(tree);
        int blossom_var = MinMinusBlossomVariable(tree);

        if (plus_empty < result.upper_bound[i]) {
            result.upper_bound[i] = plus_empty;
        }
        if (plus_plus_internal / 2 < result.upper_bound[i]) {
            result.upper_bound[i] = plus_plus_internal / 2;
        }
        if (blossom_var < result.upper_bound[i]) {
            result.upper_bound[i] = blossom_var;
        }

        // get plus_plus_constraints and plus_minus_constraints

        result.plus_plus_constraints[i] = PlusPlusExternalSlacks(tree);
        result.plus_minus_constraints[i] = PlusMinusExternalSlacks(tree);
    }

    return result;
}

void VzhuhSolver::InitNextRoundActionable() {
    for (int tree : alive_trees) {
        if (params.verbose) {
            std::cout << "updating zero slack set next to tree " << tree << std::endl;
        }

        int tree_var = trees[tree].dual_var_quadrupled;
        for (auto [other_tree, queue_index] : trees[tree].pq_plus_plus) {
            if (tree < other_tree) {
                int top_edge = GetMinEdgeHeap(queue_index);
                if (top_edge >= 0) {
                    if (edges[top_edge].slack_quadrupled_amortized_ - tree_var - trees[other_tree].dual_var_quadrupled
                        == 0) {
                        AddZeroSlackEdgesFromQueue(queue_index, true);
                    }
                }
            }
        }
        for (auto [other_tree, queue_index] : trees[tree].pq_plus_minus) {
            int top_edge = GetMinEdgeHeap(queue_index);
            if (top_edge >= 0) {
                if (edges[top_edge].slack_quadrupled_amortized_ - tree_var + trees[other_tree].dual_var_quadrupled ==
                    0) {
                    AddZeroSlackEdgesFromQueue(queue_index, false);
                }
            }
        }

        int queue_index = trees[tree].plus_empty_edges;
        int top_edge = GetMinEdgeHeap(queue_index);
        if (top_edge >= 0) {
            if (edges[top_edge].slack_quadrupled_amortized_ - tree_var == 0) {
                AddZeroSlackEdgesFromQueue(queue_index, true);
            }
        }

        CleanLoopsFromQueueTop(tree);
        queue_index = trees[tree].plus_plus_internal_edges;
        top_edge = GetMinEdgeHeap(queue_index);
        if (top_edge >= 0) {
            if (edges[top_edge].slack_quadrupled_amortized_ - 2 * tree_var == 0) {
                AddZeroSlackEdgesFromQueue(queue_index, true);
            }
        }
    }
}

void VzhuhSolver::AddZeroSlackEdgesFromQueue(int queue_index, bool add_to_actionable) {
    int top_edge = GetMinEdgeHeap(queue_index);
    RemoveMinEdgeHeap(queue_index);
    InsertEdgeHeap(top_edge, queue_index);

    int min_key = edges[top_edge].slack_quadrupled_amortized_;

    std::stack<int> stack;
    stack.push(GetMinEdgeHeap(queue_index));

    while (!stack.empty()) {
        int edge = stack.top();
        stack.pop();

        if (!edges[edge].maybe_has_zero_slack) {
            edges[edge].maybe_has_zero_slack = true;
        }
        if (add_to_actionable) {
            actionable_edges.push(edge);
        }

        for (int child = edges[edge].heap_child; child >= 0; child = edges[child].heap_next) {
            if (edges[child].slack_quadrupled_amortized_ == min_key) {
                stack.push(child);
            }
        }
    }
}

void VzhuhSolver::CleanLoopsFromQueueTop(int tree) {
    int queue_idx = trees[tree].plus_plus_internal_edges;

    int top = GetMinEdgeHeap(queue_idx);
    while (top >= 0) {
        if (Head(top) == Tail(top)) {
            RemoveEdgeFromQueue(top);
        } else {
            break;
        }

        top = GetMinEdgeHeap(queue_idx);
    }
}

bool VzhuhSolver::IsElementary(int node) const {
    return node < num_vertices_elementary;
}

int VzhuhSolver::TopBlossom(int node) const {
    while (blossom_parents[node] >= 0) {
        node = blossom_parents[node];
    }
    return node;
}

int VzhuhSolver::Receptacle(int node) {
    int grandparent = node;
    while (nodes[grandparent].receptacle_ != grandparent) {
        grandparent = nodes[grandparent].receptacle_;
    }

    int cur_node = node;
    while (nodes[cur_node].receptacle_ != grandparent) {
        int next_node = nodes[cur_node].receptacle_;
        nodes[cur_node].receptacle_ = grandparent;
        cur_node = next_node;
    }

    return grandparent;
}

int VzhuhSolver::DualVariableQuadrupled(int node) const {
    return DualVariableQuadrupled(node,
                                  nodes[node].tree,
                                  nodes[node].plus,
                                  blossom_parents[node]);
}

int VzhuhSolver::DualVariableQuadrupled(int node, int tree, bool plus, int blossom_parent) const {
    // if (!nodes[node].is_alive) {
    //     throw std::runtime_error("DualVariableQuadrupled: node is not alive");
    // }

    if (tree < 0 || blossom_parent >= 0) {
        return nodes[node].dual_var_quadrupled_amortized_;
    }
    if (plus) {
        return nodes[node].dual_var_quadrupled_amortized_ + trees[tree].dual_var_quadrupled;
    }
    return nodes[node].dual_var_quadrupled_amortized_ - trees[tree].dual_var_quadrupled;
}

std::vector<VzhuhSolver::ArcIndex> &VzhuhSolver::NonLoopNeighbors(int node) {
    if (!nodes[node].neighbors.empty()) {
        return nodes[node].neighbors;
    }

    // otherwise, compute neighbors

    // mark the vertices, collect the lists
    std::queue<int> queue;
    std::vector<std::vector<ArcIndex> *> lists;
    int total_length = 0;
    queue.push(node);
    ++nodes_label_cnt;
    nodes[node].label = nodes_label_cnt;
    while (!queue.empty()) {
        int cur = queue.front();
        queue.pop();

        for (int child : nodes[cur].blossom_children) {
            nodes[child].label = nodes_label_cnt;

            if (!nodes[child].neighbors.empty()) {
                lists.push_back(&nodes[child].neighbors);
                total_length += nodes[child].neighbors.size();
            } else {
                queue.push(child);
            }
        }
    }

    // add non-loops to neighbors
    nodes[node].neighbors.reserve(total_length);
    for (auto list_ptr : lists) {
        for (ArcIndex arc : *list_ptr) {
            // TODO only update the right end
            int edge = arc.index / 2;
            while (blossom_parents[edges[edge].head] >= 0 && nodes[edges[edge].head].label !=
                nodes_label_cnt) {
                edges[edge].head = blossom_parents[edges[edge].head];
            }
            while (blossom_parents[edges[edge].tail] >= 0 && nodes[edges[edge].tail].label !=
                nodes_label_cnt) {
                edges[edge].tail = blossom_parents[edges[edge].tail];
            }
            if (nodes[edges[edge].head].label != nodes_label_cnt) {
                nodes[node].neighbors.push_back(arc);
                edges[edge].tail = node;
            } else if (nodes[edges[edge].tail].label != nodes_label_cnt) {
                nodes[node].neighbors.push_back(arc);
                edges[edge].head = node;
            }
        }
    }

    return nodes[node].neighbors;
}

// std::vector<int> &VzhuhSolver::NonLoopZeroSlackNeighbors(int node) {
//     if (nodes[node].round_0slack_neighbors_updated == current_round) {
//         return nodes[node].zero_slack_neighbors;
//     }
//
//     for (int edge : NonLoopNeighbors(node)) {
//         if (edges[edge].maybe_has_zero_slack) {
//             nodes[node].zero_slack_neighbors.push_back(edge);
//         }
//     }
//
//     return nodes[node].zero_slack_neighbors;
// }

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
        for (int child : nodes[cur].blossom_children) {
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

    int slack = edges[edge].slack_quadrupled_amortized_;

    if (nodes[head].tree >= 0) {
        if (nodes[head].plus) {
            slack -= trees[nodes[head].tree].dual_var_quadrupled;
        } else {
            slack += trees[nodes[head].tree].dual_var_quadrupled;
        }
    }
    if (nodes[tail].tree >= 0) {
        if (nodes[tail].plus) {
            slack -= trees[nodes[tail].tree].dual_var_quadrupled;
        } else {
            slack += trees[nodes[tail].tree].dual_var_quadrupled;
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

    if (nodes[head].old_blossom_parent >= 0) {
        head = nodes[head].old_blossom_parent;
    }
    if (nodes[tail].old_blossom_parent >= 0) {
        tail = nodes[tail].old_blossom_parent;
    }

    if (head == tail) {
        // used to be a loop, but the node got expanded
        return edges[edge].slack_quadrupled_amortized_ - 2 * nodes[head].tree_var_at_birth;
    }

    int slack = edges[edge].slack_quadrupled_amortized_;

    if (nodes[head].old_tree >= 0) {
        if (nodes[head].old_plus) {
            slack -= trees[nodes[head].old_tree].dual_var_quadrupled;
        } else {
            slack += trees[nodes[head].old_tree].dual_var_quadrupled;
        }
    }
    if (nodes[tail].old_tree >= 0) {
        if (nodes[tail].old_plus) {
            slack -= trees[nodes[tail].old_tree].dual_var_quadrupled;
        } else {
            slack += trees[nodes[tail].old_tree].dual_var_quadrupled;
        }
    }

    return slack;
}

int VzhuhSolver::OtherEnd(ArcIndex arc) {
    if (arc.index % 2 == 1) {
        return Head(arc.index / 2);
    }
    return Tail(arc.index / 2);
}

int VzhuhSolver::ThisEnd(ArcIndex arc) {
    if (arc.index % 2 == 0) {
        return Head(arc.index / 2);
    }
    return Tail(arc.index / 2);
}

VzhuhSolver::ArcIndex VzhuhSolver::ReverseArc(ArcIndex arc) {
    if (arc.index % 2 == 0) {
        return ArcIndex(arc.index + 1);
    }
    return ArcIndex(arc.index - 1);
}

int VzhuhSolver::OtherElementaryEnd(ArcIndex arc) const {
    if (arc.index % 2 == 1) {
        return edges[arc.index / 2].elementary_head;
    }
    return edges[arc.index / 2].elementary_tail;
}

int VzhuhSolver::Head(int edge) {
    int head = edges[edge].head;

    if (blossom_parents[head] >= 0) {
        while (blossom_parents[head] >= 0) {
            head = blossom_parents[head];
        }
        edges[edge].head = head;
    }

    return head;
}

int VzhuhSolver::Tail(int edge) {
    int tail = edges[edge].tail;

    if (blossom_parents[tail] >= 0) {
        while (blossom_parents[tail] >= 0) {
            tail = blossom_parents[tail];
        }
        edges[edge].tail = tail;
    }

    return tail;
}

int VzhuhSolver::PlusPlusLCA(int first_vertex, int second_vertex) {
    ++nodes_label_cnt;
    first_vertex = Receptacle(first_vertex);
    second_vertex = Receptacle(second_vertex);

    nodes[first_vertex].label = nodes_label_cnt;
    nodes[second_vertex].label = nodes_label_cnt;

    while (first_vertex != second_vertex) {
        if (nodes[first_vertex].matched_edge.index >= 0) {
            first_vertex = OtherEnd(nodes[first_vertex].matched_edge);
            first_vertex = OtherEnd(nodes[first_vertex].minus_parent);
            first_vertex = Receptacle(first_vertex);
            if (nodes[first_vertex].label == nodes_label_cnt) {
                return first_vertex;
            }
            nodes[first_vertex].label = nodes_label_cnt;
        }

        if (nodes[second_vertex].matched_edge.index >= 0) {
            second_vertex = OtherEnd(nodes[second_vertex].matched_edge);
            second_vertex = OtherEnd(nodes[second_vertex].minus_parent);
            second_vertex = Receptacle(second_vertex);
            if (nodes[second_vertex].label == nodes_label_cnt) {
                return second_vertex;
            }
            nodes[second_vertex].label = nodes_label_cnt;
        }
    }

    return first_vertex;
}

void VzhuhSolver::MakeEdgeMatched(int edge) {
    if (Head(edge) == Tail(edge)) {
        throw std::runtime_error("In MakeEdgeMatched: trying to match a loop");
    }
    nodes[Head(edge)].matched_edge = ArcIndex(edge * 2);
    nodes[Tail(edge)].matched_edge = ArcIndex(edge * 2 + 1);
    edges[edge].matched = true;
}

void VzhuhSolver::MakeEdgeUnmatched(int edge) {
    edges[edge].matched = false;
}

void VzhuhSolver::AddEdgeToThisQueue(int edge, int queue_index) {
    InsertEdgeHeap(edge, queue_index);
    edges[edge].queue_index = queue_index;
}

void VzhuhSolver::RemoveEdgeFromQueue(int edge) {
    if (edges[edge].queue_index >= 0) {
        RemoveEdgeHeap(edge);
        edges[edge].queue_index = -1;
    }
}

void VzhuhSolver::AddNodeToQueue(int node, int queue_index) {
    if (nodes[node].queue_index == queue_index) {
        return;
    }

    if (nodes[node].queue_index >= 0) {
        RemoveNodeHeap(node);
    }
    InsertNodeHeap(node, queue_index);
    nodes[node].queue_index = queue_index;
}

void VzhuhSolver::RemoveNodeFromQueue(int node) {
    if (nodes[node].queue_index >= 0) {
        RemoveNodeHeap(node);
        nodes[node].queue_index = -1;
    }
}

int VzhuhSolver::TreeTreeQueueIndex(int other_tree,
                                    std::vector<std::pair<int, int> > *tree_neighbors) const {
    // may delete dead trees from the tree_neighbors
    for (int i = 0; i < static_cast<int>(tree_neighbors->size()); ++i) {
        if ((*tree_neighbors)[i].first == other_tree) {
            return (*tree_neighbors)[i].second;
        }
        if (!trees[(*tree_neighbors)[i].first].is_alive) {
            (*tree_neighbors)[i] = tree_neighbors->back();
            tree_neighbors->pop_back();
            --i;
        }
    }
    return -1;
}

int VzhuhSolver::MinPlusPlusInternalEdge(int queue_index) {
    while (GetMinEdgeHeap(queue_index) >= 0) {
        int edge = GetMinEdgeHeap(queue_index);

        int head = Head(edge);
        int tail = Tail(edge);

        // TODO maybe add IsLoop(int edge) function
        if (head != tail) {
            return edge;
        }
        RemoveMinEdgeHeap(queue_index);
        edges[edge].queue_index = -1;
    }
    return -1;
}

int VzhuhSolver::PopExpandableBlossom(int tree) {
    int node = GetMinNodeHeap(trees[tree].minus_blossoms);
    if (node >= 0) {
        if (nodes[node].dual_var_quadrupled_amortized_ - trees[tree].dual_var_quadrupled == 0) {
            RemoveMinNodeHeap(trees[tree].minus_blossoms);
            return node;
        }
    }
    return -1;
}

int VzhuhSolver::PlusEmptySlack(int tree) {
    int edge = GetMinEdgeHeap(trees[tree].plus_empty_edges);
    if (edge >= 0) {
        return edges[edge].slack_quadrupled_amortized_ - trees[tree].dual_var_quadrupled;
    }
    return INT32_MAX;
}

int VzhuhSolver::PlusPlusInternalSlack(int tree) {
    int edge = MinPlusPlusInternalEdge(trees[tree].plus_plus_internal_edges);
    if (edge >= 0) {
        return edges[edge].slack_quadrupled_amortized_ - 2 * trees[tree].dual_var_quadrupled;
    }
    return INT32_MAX;
}

int VzhuhSolver::MinMinusBlossomVariable(int tree) const {
    int node = GetMinNodeHeap(trees[tree].minus_blossoms);
    if (node >= 0) {
        return nodes[node].dual_var_quadrupled_amortized_ - trees[tree].dual_var_quadrupled;
    }
    return INT32_MAX;
}

std::vector<std::pair<int, int> > VzhuhSolver::PlusPlusExternalSlacks(int tree) {
    std::vector<std::pair<int, int> > result;
    result.reserve(trees[tree].pq_plus_plus.size());

    int tree_var = trees[tree].dual_var_quadrupled;
    for (int i = 0; i < static_cast<int>(trees[tree].pq_plus_plus.size()); ++i) {
        auto [tree_neighbor, queue_index] = trees[tree].pq_plus_plus[i];
        if (trees[tree_neighbor].is_alive) {
            int edge = GetMinEdgeHeap(queue_index);
            if (edge >= 0) {
                result.emplace_back(trees[tree_neighbor].alive_index,
                                    edges[edge].slack_quadrupled_amortized_ - tree_var - trees[tree_neighbor].
                                    dual_var_quadrupled);
            }
        } else {
            trees[tree].pq_plus_plus[i] = trees[tree].pq_plus_plus.back();
            trees[tree].pq_plus_plus.pop_back();
            --i;
        }
    }
    return result;
}

std::vector<std::pair<int, int> > VzhuhSolver::PlusMinusExternalSlacks(int tree) {
    // TODO avoid code duplication?

    std::vector<std::pair<int, int> > result;
    result.reserve(trees[tree].pq_plus_minus.size());

    int tree_var = trees[tree].dual_var_quadrupled;
    for (int i = 0; i < static_cast<int>(trees[tree].pq_plus_minus.size()); ++i) {
        auto [tree_neighbor, queue_index] = trees[tree].pq_plus_minus[i];
        if (trees[tree_neighbor].is_alive) {
            int edge = GetMinEdgeHeap(queue_index);
            if (edge >= 0) {
                result.emplace_back(trees[tree_neighbor].alive_index,
                                    edges[edge].slack_quadrupled_amortized_ - tree_var + trees[tree_neighbor].
                                    dual_var_quadrupled);
            }
        } else {
            trees[tree].pq_plus_minus[i] = trees[tree].pq_plus_minus.back();
            trees[tree].pq_plus_minus.pop_back();
            --i;
        }
    }
    return result;
}

void VzhuhSolver::AddNodeToRecord(int node) {
    if (nodes[node].is_in_record) {
        return;
    }
    primal_update_record.push_back(node);
    nodes[node].is_in_record = true;
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

int VzhuhSolver::GetMinEdgeHeap(int heap_index) const {
    return edge_heaps[heap_index].root;
}

void VzhuhSolver::InsertEdgeHeap(int edge, int heap_index) {
    edges[edge].heap_child = -1;
    edges[edge].heap_next = -1;
    edges[edge].heap_prev = -1;

    edge_heaps[heap_index].root = MeldEdgeHeap(edge_heaps[heap_index].root, edge);
}

void VzhuhSolver::RemoveMinEdgeHeap(int heap_index) {
    if (edge_heaps[heap_index].root < 0) {
        return;
    }

    int old_root = edge_heaps[heap_index].root;
    int children = edges[edge_heaps[heap_index].root].heap_child;

    if (children >= 0) {
        edges[children].heap_prev = -1; // new root candidates
    }

    edge_heaps[heap_index].root = TwoPassMergeEdgeHeap(children);
    edges[old_root].heap_child = -1;
}

void VzhuhSolver::RemoveEdgeHeap(int edge) {
    int heap_index = edges[edge].queue_index;
    if (edge == edge_heaps[heap_index].root) {
        RemoveMinEdgeHeap(heap_index);
        return;
    }

    CutEdgeHeap(edge);

    int subtree = TwoPassMergeEdgeHeap(edges[edge].heap_child);
    if (subtree >= 0) {
        edges[subtree].heap_prev = -1;
    }

    edge_heaps[heap_index].root = MeldEdgeHeap(edge_heaps[heap_index].root, subtree);

    edges[edge].heap_child = -1;
}

int VzhuhSolver::MeldEdgeHeap(int edge_a, int edge_b) {
    if (edge_a < 0) {
        return edge_b;
    }
    if (edge_b < 0) {
        return edge_a;
    }

    if (edges[edge_b].slack_quadrupled_amortized_ < edges[edge_a].slack_quadrupled_amortized_) {
        std::swap(edge_a, edge_b);
    }

    edges[edge_b].heap_prev = edge_a;
    edges[edge_b].heap_next = edges[edge_a].heap_child;
    if (edges[edge_a].heap_child >= 0) {
        edges[edges[edge_a].heap_child].heap_prev = edge_b;
    }

    edges[edge_a].heap_child = edge_b;

    return edge_a;
}

int VzhuhSolver::TwoPassMergeEdgeHeap(int edge_first) {
    if (edge_first < 0 || edges[edge_first].heap_next < 0) {
        return edge_first;
    }

    int a = edge_first;
    int b = edges[a].heap_next;
    int rest = edges[b].heap_next;

    edges[a].heap_next = -1;
    edges[b].heap_next = -1;

    int merged = MeldEdgeHeap(a, b);
    // TODO avoid recursion
    int remaining = TwoPassMergeEdgeHeap(rest);

    return MeldEdgeHeap(merged, remaining);
}

void VzhuhSolver::CutEdgeHeap(int edge) {
    if (edges[edge].heap_prev < 0) {
        return; // already root
    }

    if (edges[edges[edge].heap_prev].heap_child == edge) {
        // node is leftmost child
        int parent = edges[edge].heap_prev;
        edges[parent].heap_child = edges[edge].heap_next;
        if (edges[edge].heap_next >= 0) {
            edges[edges[edge].heap_next].heap_prev = parent;
        }
    } else {
        // node has left sibling
        int left_sibling = edges[edge].heap_prev;
        edges[left_sibling].heap_next = edges[edge].heap_next;
        if (edges[edge].heap_next >= 0) {
            edges[edges[edge].heap_next].heap_prev = left_sibling;
        }
    }

    edges[edge].heap_prev = -1;
    edges[edge].heap_next = -1;
}

int VzhuhSolver::GetMinNodeHeap(int heap_index) const {
    return node_heaps[heap_index].root;
}

void VzhuhSolver::InsertNodeHeap(int node, int heap_index) {
    nodes[node].heap_child = -1;
    nodes[node].heap_next = -1;
    nodes[node].heap_prev = -1;

    node_heaps[heap_index].root = MeldNodeHeap(node_heaps[heap_index].root, node);
}

void VzhuhSolver::RemoveMinNodeHeap(int heap_index) {
    if (node_heaps[heap_index].root < 0) {
        return;
    }

    int old_root = node_heaps[heap_index].root;
    int children = nodes[node_heaps[heap_index].root].heap_child;

    if (children >= 0) {
        nodes[children].heap_prev = -1; // new root candidates
    }

    node_heaps[heap_index].root = TwoPassMergeNodeHeap(children);
    nodes[old_root].heap_child = -1;
}

void VzhuhSolver::RemoveNodeHeap(int node) {
    int heap_index = nodes[node].queue_index;
    if (node == node_heaps[heap_index].root) {
        RemoveMinNodeHeap(heap_index);
        return;
    }

    CutNodeHeap(node);

    int subtree = TwoPassMergeNodeHeap(nodes[node].heap_child);
    if (subtree >= 0) {
        nodes[subtree].heap_prev = -1;
    }

    node_heaps[heap_index].root = MeldNodeHeap(node_heaps[heap_index].root, subtree);

    nodes[node].heap_child = -1;
}

int VzhuhSolver::MeldNodeHeap(int node_a, int node_b) {
    if (node_a < 0) {
        return node_b;
    }
    if (node_b < 0) {
        return node_a;
    }

    if (nodes[node_b].dual_var_quadrupled_amortized_ < nodes[node_a].dual_var_quadrupled_amortized_) {
        std::swap(node_a, node_b);
    }

    nodes[node_b].heap_prev = node_a;
    nodes[node_b].heap_next = nodes[node_a].heap_child;
    if (nodes[node_a].heap_child >= 0) {
        nodes[nodes[node_a].heap_child].heap_prev = node_b;
    }

    nodes[node_a].heap_child = node_b;

    return node_a;
}

int VzhuhSolver::TwoPassMergeNodeHeap(int node_first) {
    if (node_first < 0 || nodes[node_first].heap_next < 0) {
        return node_first;
    }

    int a = node_first;
    int b = nodes[a].heap_next;
    int rest = nodes[b].heap_next;

    nodes[a].heap_next = -1;
    nodes[b].heap_next = -1;

    int merged = MeldNodeHeap(a, b);
    // TODO avoid recursion
    int remaining = TwoPassMergeNodeHeap(rest);

    return MeldNodeHeap(merged, remaining);
}

void VzhuhSolver::CutNodeHeap(int node) {
    if (nodes[node].heap_prev < 0) {
        return; // already root
    }

    if (nodes[nodes[node].heap_prev].heap_child == node) {
        // node is leftmost child
        int parent = nodes[node].heap_prev;
        nodes[parent].heap_child = nodes[node].heap_next;
        if (nodes[node].heap_next >= 0) {
            nodes[nodes[node].heap_next].heap_prev = parent;
        }
    } else {
        // node has left sibling
        int left_sibling = nodes[node].heap_prev;
        nodes[left_sibling].heap_next = nodes[node].heap_next;
        if (nodes[node].heap_next >= 0) {
            nodes[nodes[node].heap_next].heap_prev = left_sibling;
        }
    }

    nodes[node].heap_prev = -1;
    nodes[node].heap_next = -1;
}
