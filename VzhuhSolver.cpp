#include "VzhuhSolver.h"

#include <queue>
#include <ranges>
#include <unordered_set>

VzhuhSolver::VzhuhSolver(const std::vector<std::tuple<int, int, int> > &edge_list_,
                         const SolverParameters &params_) : primal_objective(INT64_MAX), dual_objective(INT64_MIN),
                                                            params(params_),
                                                            num_vertices_elementary(InitNumVertices(edge_list_)) {
    nodes.reserve(2 * num_vertices_elementary);
    for (int i = 0; i < num_vertices_elementary; ++i) {
        nodes.emplace_back(i);
    }

    adj_list = std::vector<std::vector<EdgeIndex> >(num_vertices_elementary, std::vector<EdgeIndex>());
    zero_slack_adj_list = std::vector<std::vector<EdgeIndex> >(num_vertices_elementary, std::vector<EdgeIndex>());
    edges.reserve(edge_list_.size());
    for (int i = 0; i < static_cast<int>(edge_list_.size()); ++i) {
        int head = std::get<0>(edge_list_[i]);
        int tail = std::get<1>(edge_list_[i]);
        edges.emplace_back(head, tail, std::get<2>(edge_list_[i]));
        adj_list[head].emplace_back(i);
        adj_list[tail].emplace_back(i);
    }

    nodes_label_cnt = 1;
    num_trees_alive = INT32_MAX;
}

void VzhuhSolver::FindMinPerfectMatching() {
    GreedyInit();
    InitializeTrees();

    if (params.verbose) {
        PrintGraph();
    }

    int num_rounds = 0;
    while (!alive_trees.empty()) {
        ++num_rounds;
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

bool VzhuhSolver::NodeComparator::operator()(const NodeIndex &a, const NodeIndex &b) const {
    return solver->nodes[a.index].dual_var_quadrupled_amortized_ < solver->nodes[b.index].
        dual_var_quadrupled_amortized_;
}

bool VzhuhSolver::EdgeComparator::operator()(const EdgeIndex &a, const EdgeIndex &b) const {
    return solver->edges[a.index].slam_quadrupled_amortized_ < solver->edges[b.index].slam_quadrupled_amortized_;
}

VzhuhSolver::Edge::Edge(int head_, int tail_, int weight_) : head(head_), tail(tail_), elementary_head(head_),
                                                             elementary_tail(tail_) {
    handle = nullptr;
    queue_index = -1;
    weight = weight_;
    slam_quadrupled_amortized_ = 4 * weight_;
    matched = false;
    is_in_zero_slack_set = false;
}

VzhuhSolver::Node::Node(int index) : matched_edge(-1), blossom_parent(-1), minus_parent(-1),
                                     receptacle_(index), tree(-1), old_tree(-1), old_blossom_parent(-1) {
    is_alive = true;
    dual_var_quadrupled_amortized_ = 0;
    plus = false;
    old_plus = false;
    queue_index = -1;
    handle = nullptr;
    label = 0;
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

void VzhuhSolver::PrintGraph() {
    std::cout << "Adjacency list (to, weight, slam, matched):" << std::endl;

    for (int i = 0; i < num_vertices_elementary; ++i) {
        std::cout << i << ": ";
        for (EdgeIndex edge : adj_list[i]) {
            std::cout << "(" << OtherElementaryEnd(edge, NodeIndex(i)).index << " " << edges[edge.index].weight << " "
                << SlamQuadrupled(edge) / 4. << " " << edges[edge.index].matched << ") ";
        }
        std::cout << std::endl;
    }

    std::cout << "Zero slack adjacency list (to, weight, slam, matched):" << std::endl;

    for (int i = 0; i < num_vertices_elementary; ++i) {
        std::cout << i << ": ";
        for (EdgeIndex edge : zero_slack_adj_list[i]) {
            std::cout << "(" << OtherElementaryEnd(edge, NodeIndex(i)).index << " " << edges[edge.index].weight << " "
                << SlamQuadrupled(edge) / 4. << " " << edges[edge.index].matched << ") ";
        }
        std::cout << std::endl;
    }

    std::cout << "Node structure: " << std::endl;
    for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
        if (nodes[i].is_alive) {
            PrintNode(NodeIndex(i));
        }
    }

    std::cout << "Tree structure: " << std::endl;
    for (TreeIndex tree : alive_trees) {
        if (!trees[tree.index].is_alive) {
            continue;
        }
        std::cout << "root: " << trees[tree.index].root.index << " var: " << trees[tree.index].dual_var_quadrupled / 4.
            << std::endl;
    }
}

void VzhuhSolver::PrintNode(NodeIndex node) {
    std::cout << node.index << ": y_v = " << DualVariableQuadrupled(node) / 4.;
    if (nodes[node.index].blossom_parent) {
        std::cout << " blossom_parent: " << nodes[node.index].blossom_parent.index;
    }
    std::cout << std::endl;
}

void VzhuhSolver::GreedyInit() {
    // first, make all the slacks non-negative
    for (int node_index = 0; node_index < num_vertices_elementary; ++node_index) {
        if (adj_list[node_index].empty()) {
            std::cout << "vertex " << node_index << ": no neighbors" << std::endl;
            throw std::runtime_error("Found an isolated vertex => no perfect matching exists");
        }

        int min_weight = edges[adj_list[node_index].front().index].weight;
        for (EdgeIndex edge : adj_list[node_index]) {
            if (edges[edge.index].weight < min_weight) {
                min_weight = edges[edge.index].weight;
            }
        }

        nodes[node_index].dual_var_quadrupled_amortized_ += 2 * min_weight;
        for (EdgeIndex edge : adj_list[node_index]) {
            edges[edge.index].slam_quadrupled_amortized_ -= 2 * min_weight;
        }
    }

    for (int node_index = 0; node_index < num_vertices_elementary; ++node_index) {
        if (nodes[node_index].matched_edge) {
            continue;
        }

        EdgeIndex smallest_slack_edge = adj_list[node_index].front();
        for (const EdgeIndex edge : adj_list[node_index]) {
            if (edges[edge.index].slam_quadrupled_amortized_ < edges[smallest_slack_edge.index].
                slam_quadrupled_amortized_) {
                smallest_slack_edge = edge;
            }
        }

        int diff = edges[smallest_slack_edge.index].slam_quadrupled_amortized_;
        nodes[node_index].dual_var_quadrupled_amortized_ += diff;
        for (EdgeIndex edge : adj_list[node_index]) {
            edges[edge.index].slam_quadrupled_amortized_ -= diff;
        }
        if (!nodes[OtherElementaryEnd(smallest_slack_edge, NodeIndex(node_index)).index].matched_edge) {
            // if the other vertex is also unmatched, match the edge
            edges[smallest_slack_edge.index].matched = true;
            nodes[node_index].matched_edge = smallest_slack_edge;
            nodes[OtherElementaryEnd(smallest_slack_edge, NodeIndex(node_index)).index].matched_edge =
                smallest_slack_edge;
        }
    }

    // update zero_slack_adj_list
    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
        EdgeIndex edge(i);
        if (edges[edge.index].slam_quadrupled_amortized_ == 0) {
            zero_slack_adj_list[edges[edge.index].head.index].push_back(edge);
            zero_slack_adj_list[edges[edge.index].tail.index].push_back(edge);
            edges[edge.index].is_in_zero_slack_set = true;
        }
    }
}

void VzhuhSolver::InitializeTrees() {
    // is called after GreedyInit
    // also initializes queues and actionable_edges

    std::vector<int> roots;
    for (int root_index = 0; root_index < num_vertices_elementary; ++root_index) {
        if (nodes[root_index].matched_edge) {
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
        trees.emplace_back(roots[i], i, 2 * i, 2 * i + 1);
        alive_trees.emplace_back(i);

        node_heaps.emplace_back(std::make_unique<NodeHeap>(params.heap_arity, NodeComparator{this}));
        edge_heaps.emplace_back(std::make_unique<EdgeHeap>(params.heap_arity, EdgeComparator{this}));
        edge_heaps.emplace_back(std::make_unique<EdgeHeap>(params.heap_arity, EdgeComparator{this}));

        nodes[roots[i]].tree = TreeIndex(i);
        nodes[roots[i]].plus = true;
        nodes[roots[i]].old_tree = TreeIndex(i);
        nodes[roots[i]].old_plus = true;
    }

    // initialize queues and actionable_edges
    for (int root_index : roots) {
        for (EdgeIndex edge_to_neighbor : adj_list[root_index]) {
            NodeIndex neighbor = OtherEnd(edge_to_neighbor, NodeIndex(root_index));

            if (!nodes[neighbor.index].tree) {
                AddEdgeToThisQueue(edge_to_neighbor, trees[nodes[root_index].tree.index].plus_empty_edges);
            } else {
                // neighbor is another root
                AddPQPlusPlus(nodes[root_index].tree, nodes[neighbor.index].tree, edge_to_neighbor);
            }

            if (edges[edge_to_neighbor.index].slam_quadrupled_amortized_ == 0) {
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

    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
        if (edges[i].matched) {
            matching.emplace_back(edges[i].elementary_head.index, edges[i].elementary_tail.index);
        }
    }
}

void VzhuhSolver::ComputePrimalObjective() {
    primal_objective = 0;
    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
        if (edges[i].matched) {
            primal_objective += edges[i].weight;
        }
    }
}

void VzhuhSolver::DestroyBlossoms() {
    for (int i = static_cast<int>(nodes.size()) - 1; i >= num_vertices_elementary; --i) {
        if (!nodes[i].is_alive) {
            continue;
        }
        if (nodes[i].blossom_parent) {
            throw std::runtime_error("In DestroyBlossoms: not going top down");
        }

        if (params.verbose) {
            std::cout << "destroying " << i << std::endl;
            PrintGraph();
        }

        NodeIndex blossom = NodeIndex(i);
        nodes[blossom.index].is_alive = false;

        // find the new receptacle
        NodeIndex new_receptacle(-1);
        NodeIndex head = edges[nodes[blossom.index].matched_edge.index].elementary_head;
        while (nodes[head.index].blossom_parent) {
            if (nodes[head.index].blossom_parent == blossom) {
                new_receptacle = head;
                break;
            }
            head = nodes[head.index].blossom_parent;
        }
        if (!new_receptacle) {
            NodeIndex tail = edges[nodes[blossom.index].matched_edge.index].elementary_tail;
            while (nodes[tail.index].blossom_parent) {
                if (nodes[tail.index].blossom_parent == blossom) {
                    new_receptacle = tail;
                    break;
                }
                tail = nodes[tail.index].blossom_parent;
            }
        }
        if (!new_receptacle) {
            throw std::runtime_error("DestroyBlossoms: no new receptacle found");
        }

        RestoreEdgeEndsBeforeExpand(blossom);

        UpdateMatching(blossom, new_receptacle);
        nodes[new_receptacle.index].matched_edge = nodes[blossom.index].matched_edge;
    }

    if (params.verbose) {
        std::cout << "destroyed all"<< std::endl;
        PrintGraph();
    }
}

void VzhuhSolver::ComputeDualCertificate() {
    // can be called only before we destroy the blossoms in the end

    if (!dual_certificate.empty()) {
        throw std::runtime_error("Dual certificate is already non-empty");
    }

    std::vector<int> node_index_alive(nodes.size(), -1);
    dual_certificate.reserve(nodes.size());
    int index = 0;
    for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
        if (nodes[i].is_alive) {
            node_index_alive[i] = index;
            dual_certificate.emplace_back(index, DualVariableQuadrupled(NodeIndex(i)), -1);
            for (NodeIndex child : nodes[i].blossom_children) {
                std::get<2>(dual_certificate[node_index_alive[child.index]]) = index;
            }
            ++index;
        }
    }
}

void VzhuhSolver::ComputeDualObjectiveQuadrupled() {
    dual_objective = 0;
    for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
        if (nodes[i].is_alive) {
            dual_objective += DualVariableQuadrupled(NodeIndex(i));
        }
    }
}

bool VzhuhSolver::MakePrimalUpdates() {
    bool action_taken = false;
    PrimalUpdateRecord record;

    std::vector<int> variables;
    std::vector<int> slams;
    if (params.debug) {
        variables = NodeVariables();
        slams = EdgeSlams();
    }

    // expand phase
    for (TreeIndex tree : alive_trees) {
        NodeIndex blossom = PopExpandableBlossom(tree);
        while (blossom) {
            Expand(blossom, &record);
            blossom = PopExpandableBlossom(tree);
        }
    }

    // unweighted solver like phase
    while (!actionable_edges.empty()) {
        EdgeIndex edge = actionable_edges.front();
        actionable_edges.pop();
        MakePrimalUpdate(edge, &record);
        // if (num_trees_alive == 0) {
        //     return;
        // }
    }

    if (!record.changed_sign.empty()) {
        action_taken = true;
    }

    if (params.debug) {
        ValidateEvenOddPaths();
    }

    // shrink phase
    if (params.verbose) {
        std::cout << "shrinking phase" << std::endl;
    }
    DeleteDuplicates(&record.changed_sign);
    std::vector<std::vector<NodeIndex> > future_blossoms = OrganizeBlossomChildren(record);
    for (std::vector<NodeIndex> &children : future_blossoms) {
        Shrink(children);
    }

    // update queues and amortized variables phase
    UpdateQueues(record);
    if (params.debug) {
        ValidateQueues();
    }

    if (params.debug) {
        std::vector<int> new_variables = NodeVariables();
        std::vector<int> new_slams = EdgeSlams();
        for (int i = 0; i < static_cast<int>(variables.size()); ++i) {
            if (!nodes[i].is_alive) {
                continue;
            }
            if (variables[i] != new_variables[i]) {
                std::cout << i << std::endl;
                throw std::runtime_error("variables do not match");
            }
        }
        for (int i = 0; i < static_cast<int>(slams.size()); ++i) {
            if (slams[i] != new_slams[i]) {
                std::cout << "edge " << edges[i].elementary_head.index << " " << edges[i].elementary_tail.index <<
                    std::endl;
                throw std::runtime_error("slams do not match");
            }
        }
    }

    return action_taken;
}

void VzhuhSolver::MakePrimalUpdate(EdgeIndex edge, PrimalUpdateRecord *record) {
    NodeIndex parent = Head(edge);
    NodeIndex child = Tail(edge);

    if (parent == child || OldSlamQuadrupled(edge) != 0) {
        return;
    }

    if (!nodes[parent.index].tree) {
        std::swap(parent, child);
    }

    if (!nodes[parent.index].tree) {
        // edge not adjacent to any trees
        return;
    }

    if (nodes[parent.index].tree == nodes[child.index].tree) {
        if (nodes[parent.index].plus && nodes[child.index].plus) {
            MakeCherryBlossom(edge, record);
            return;
        }
    }

    if (!nodes[child.index].tree) {
        if (nodes[parent.index].plus) {
            Grow(parent, edge, record);
            return;
        }
    }

    if (nodes[parent.index].tree != nodes[child.index].tree) {
        if (nodes[parent.index].plus && nodes[child.index].plus) {
            Augment(edge, record);
        }
    }
}

void VzhuhSolver::Expand(NodeIndex blossom, PrimalUpdateRecord *record) {
    if (params.verbose) {
        std::cout << "EXPAND " << blossom.index << std::endl;
    }

    nodes[blossom.index].is_alive = false;

    RestoreEdgeEndsBeforeExpand(blossom);
    // now OtherEnd is safe

    NodeIndex new_receptacle = edges[nodes[blossom.index].matched_edge.index].head;
    if (nodes[new_receptacle.index].old_blossom_parent != blossom) {
        new_receptacle = edges[nodes[blossom.index].matched_edge.index].tail;
    }
    NodeIndex old_receptacle = Receptacle(new_receptacle);
    NodeIndex elder_child = edges[nodes[blossom.index].minus_parent.index].head;
    if (nodes[elder_child.index].old_blossom_parent != blossom) {
        elder_child = edges[nodes[blossom.index].minus_parent.index].tail;
    }

    UpdateInternalStructure(blossom, old_receptacle, new_receptacle, elder_child, record);

    for (NodeIndex child : nodes[blossom.index].blossom_children) {
        record->changed_sign.push_back(child);
        AddNeighborsToActionable(child);
        if (nodes[child.index].tree) {
            trees[nodes[child.index].tree.index].tree_nodes.push_back(child);
        }
    }
}

void VzhuhSolver::RestoreEdgeEndsBeforeExpand(NodeIndex blossom) {
    // makes head and tail of the adjacent edges to point to the children of the blossom,
    // make children top blossoms

    for (NodeIndex child : nodes[blossom.index].blossom_children) {
        std::vector<NodeIndex> elementary_descendants = ElementaryBlossomDescendants(child);
        for (NodeIndex descendant : elementary_descendants) {
            for (EdgeIndex edge : adj_list[descendant.index]) {
                if (edges[edge.index].elementary_head == descendant) {
                    edges[edge.index].head = child;
                } else {
                    edges[edge.index].tail = child;
                }
            }
        }
    }

    for (NodeIndex child : nodes[blossom.index].blossom_children) {
        nodes[child.index].blossom_parent = NodeIndex(-1);
    }
}

void VzhuhSolver::RotateReceptacle(NodeIndex blossom, NodeIndex new_receptacle) {
    // changes matched edges and minus_parents

    if (params.verbose) {
        std::cout << "ROTATE RECEPTACLE " << blossom.index << " " << new_receptacle.index << std::endl;
    }

    NodeIndex old_receptacle = Receptacle(new_receptacle);
    if (old_receptacle == new_receptacle) {
        nodes[new_receptacle.index].matched_edge = nodes[blossom.index].matched_edge;
        return;
    }

    std::vector<EdgeIndex> even_path = EvenPathToReceptacle(new_receptacle);
    std::vector<EdgeIndex> odd_path = OddPathToReceptacle(new_receptacle);

    if (even_path.back() == odd_path.back()) {
        int i = 0;
        NodeIndex common_node = old_receptacle;
        while (even_path[even_path.size() - 1 - i] == odd_path[odd_path.size() - 1 - i]) {
            common_node = OtherEnd(even_path[even_path.size() - 1 - i], common_node);
            ++i;
        }

        RotateReceptacle(blossom, common_node);
        RotateReceptacle(blossom, new_receptacle);
        return;
    }

    NodeIndex cur_node = new_receptacle;
    bool match = false;
    for (EdgeIndex edge : even_path) {
        if (match) {
            MakeEdgeMatched(edge);
            cur_node = OtherEnd(edge, cur_node);
        } else {
            MakeEdgeUnmatched(edge);
            nodes[cur_node.index].minus_parent = edge;
            cur_node = OtherEnd(edge, cur_node);
            nodes[cur_node.index].minus_parent = edge;
        }
        match = !match;
    }

    cur_node = new_receptacle;
    match = false;
    for (EdgeIndex edge : odd_path) {
        if (match) {
            cur_node = OtherEnd(edge, cur_node);
        } else {
            nodes[cur_node.index].minus_parent = edge;
            cur_node = OtherEnd(edge, cur_node);
            nodes[cur_node.index].minus_parent = edge;
        }
        match = !match;
    }

    nodes[new_receptacle.index].minus_parent = EdgeIndex(-1);
    nodes[new_receptacle.index].matched_edge = nodes[blossom.index].matched_edge;
    nodes[new_receptacle.index].receptacle_ = new_receptacle;
    nodes[old_receptacle.index].receptacle_ = new_receptacle;
}

void VzhuhSolver::UpdateMatching(NodeIndex blossom, NodeIndex new_receptacle) {
    NodeIndex old_receptacle = Receptacle(new_receptacle);
    NodeIndex cur_node = new_receptacle;

    std::vector<EdgeIndex> path = EvenPathToReceptacle(new_receptacle);
    bool match = false;
    for (EdgeIndex edge : path) {
        if (match) {
            MakeEdgeMatched(edge);
        } else {
            MakeEdgeUnmatched(edge);
        }
        match = !match;
    }
    nodes[new_receptacle.index].matched_edge = nodes[blossom.index].matched_edge;
}

void VzhuhSolver::UpdateInternalStructure(NodeIndex blossom,
                                          NodeIndex old_receptacle,
                                          NodeIndex new_receptacle,
                                          NodeIndex elder_child,
                                          PrimalUpdateRecord *record) {
    // after expand, the remaining part of the blossom is a path in the tree

    RotateReceptacle(blossom, new_receptacle);
    std::vector<EdgeIndex> path = EvenPathToReceptacle(elder_child);

    // the path that stays in the tree
    NodeIndex cur_node = elder_child;
    nodes[elder_child.index].minus_parent = nodes[blossom.index].minus_parent;
    nodes[elder_child.index].plus = false;
    nodes[elder_child.index].tree = nodes[blossom.index].tree;
    for (EdgeIndex edge : path) {
        cur_node = OtherEnd(edge, cur_node);

        if (edges[edge.index].matched) {
            nodes[cur_node.index].minus_parent = EdgeIndex(-1);
            nodes[cur_node.index].plus = true;
        } else {
            nodes[cur_node.index].minus_parent = edge;
            nodes[cur_node.index].plus = false;
        }
        nodes[cur_node.index].tree = nodes[blossom.index].tree;
    }

    // mark the nodes that stay + insert them into the record
    ++nodes_label_cnt;
    cur_node = elder_child;
    for (EdgeIndex edge : path) {
        nodes[cur_node.index].label = nodes_label_cnt;
        record->changed_sign.push_back(cur_node);
        cur_node = OtherEnd(edge, cur_node);
    }
    nodes[new_receptacle.index].label = nodes_label_cnt;
    record->changed_sign.push_back(new_receptacle);

    // clear the part that goes to waste
    for (NodeIndex child : nodes[blossom.index].blossom_children) {
        if (nodes[child.index].label != nodes_label_cnt) {
            nodes[child.index].plus = false;
            nodes[child.index].minus_parent = EdgeIndex(-1);
            nodes[child.index].tree = TreeIndex(-1);
        }
    }

    // reset receptacles
    for (NodeIndex child : nodes[blossom.index].blossom_children) {
        nodes[child.index].receptacle_ = child;
    }
}

std::vector<VzhuhSolver::EdgeIndex> VzhuhSolver::EvenPathToReceptacle(NodeIndex node) {
    NodeIndex receptacle = Receptacle(node);
    std::vector<EdgeIndex> path;

    while (node != receptacle) {
        NodeIndex parent = OtherEnd(nodes[node.index].matched_edge, node);
        NodeIndex grandparent = OtherEnd(nodes[parent.index].minus_parent, parent);

        path.push_back(nodes[node.index].matched_edge);
        path.push_back(nodes[parent.index].minus_parent);

        node = grandparent;

        if (path.size() > num_vertices_elementary) {
            std::cout << "receptacle: " << receptacle.index << ", node: " << node.index << std::endl;
            throw std::runtime_error("EvenPathToReceptacle: infinite loop");
        }
    }

    return path;
}

std::vector<VzhuhSolver::EdgeIndex> VzhuhSolver::OddPathToReceptacle(NodeIndex node) {
    NodeIndex receptacle = Receptacle(node);
    if (receptacle == node) {
        throw std::runtime_error("OddPathToReceptacle: node is receptacle");
    }

    std::vector<EdgeIndex> path = {nodes[node.index].minus_parent};
    node = OtherEnd(nodes[node.index].minus_parent, node);
    while (node != receptacle) {
        NodeIndex parent = OtherEnd(nodes[node.index].matched_edge, node);
        NodeIndex grandparent = OtherEnd(nodes[parent.index].minus_parent, parent);

        path.push_back(nodes[node.index].matched_edge);
        path.push_back(nodes[parent.index].minus_parent);

        node = grandparent;
    }

    return path;
}

void VzhuhSolver::ExpandChildBeforeGrow(NodeIndex blossom, PrimalUpdateRecord *record) {
    nodes[blossom.index].is_alive = false;

    RestoreEdgeEndsBeforeExpand(blossom);
    // now OtherEnd, Head, Tail is safe

    NodeIndex new_receptacle = edges[nodes[blossom.index].matched_edge.index].head;
    if (nodes[new_receptacle.index].old_blossom_parent != blossom) {
        new_receptacle = edges[nodes[blossom.index].matched_edge.index].tail;
    }
    UpdateMatching(blossom, new_receptacle);

    for (NodeIndex child : nodes[blossom.index].blossom_children) {
        nodes[child.index].receptacle_ = child;
        nodes[child.index].plus = false;
        nodes[child.index].minus_parent = EdgeIndex(-1);
        nodes[child.index].tree = TreeIndex(-1);
    }

    for (NodeIndex child : nodes[blossom.index].blossom_children) {
        record->changed_sign.push_back(child);
        if (nodes[child.index].tree) {
            AddNeighborsToActionable(child);
        }
    }
}

void VzhuhSolver::Grow(NodeIndex parent, EdgeIndex edge, PrimalUpdateRecord *record) {
    NodeIndex child = OtherEnd(edge, parent);
    TreeIndex tree = nodes[parent.index].tree;

    if (!nodes[parent.index].plus) {
        throw std::runtime_error("In Grow: parent is not a plus");
    }
    if (nodes[child.index].tree) {
        throw std::runtime_error("In Grow: child vertex is not free");
    }
    if (!nodes[child.index].matched_edge) {
        std::cout << "parent " << parent.index << ", child " << child.index << std::endl;
        throw std::runtime_error("Grow: child has no matched_edge");
    }

    NodeIndex grandchild = OtherEnd(nodes[child.index].matched_edge, child);
    if (params.verbose) {
        std::cout << "GROW " << parent.index << " " << child.index << " " << grandchild.index << std::endl;
    }

    if (!IsElementary(child)) {
        if (DualVariableQuadrupled(child, nodes[child.index].old_tree, nodes[child.index].old_plus, nodes[child.index].old_blossom_parent) == 0) {
            if (params.verbose) {
                std::cout << "EXPAND CHILD, new child/grandchild: ";
            }
            ExpandChildBeforeGrow(child, record);
            child = OtherEnd(edge, parent);
            grandchild = OtherEnd(nodes[child.index].matched_edge, child);
            if (params.verbose) {
                std::cout << child.index << " " << grandchild.index << std::endl;
            }
        }
    }

    nodes[child.index].minus_parent = edge;
    nodes[child.index].tree = tree;
    nodes[child.index].plus = false;

    nodes[grandchild.index].tree = tree;
    nodes[grandchild.index].plus = true;

    record->changed_sign.push_back(child);
    record->changed_sign.push_back(grandchild);

    trees[tree.index].tree_nodes.push_back(child);
    trees[tree.index].tree_nodes.push_back(grandchild);

    AddNeighborsToActionable(grandchild);
}

void VzhuhSolver::MakeCherryBlossom(EdgeIndex edge_plus_plus, PrimalUpdateRecord *record) {
    NodeIndex head = Head(edge_plus_plus);
    NodeIndex tail = Tail(edge_plus_plus);

    if (Receptacle(head) == Receptacle(tail)) {
        return;
    }

    if (params.verbose) {
        std::cout << "MAKE CHERRY BLOSSOM " << head.index << " " << tail.index << std::endl;
    }

    auto [first_bound, second_bound] = CherryPathBounds(head, tail);

    UpdateCherryPath(head, first_bound, record);
    UpdateCherryPath(tail, second_bound, record);

    if (head != first_bound) {
        nodes[head.index].minus_parent = edge_plus_plus;
    }
    if (tail != second_bound) {
        nodes[tail.index].minus_parent = edge_plus_plus;
    }
}

std::pair<VzhuhSolver::NodeIndex, VzhuhSolver::NodeIndex> VzhuhSolver::CherryPathBounds(NodeIndex first_vertex,
    NodeIndex second_vertex) {
    NodeIndex lca = PlusPlusLCA(first_vertex, second_vertex);
    NodeIndex lca_receptacle = Receptacle(lca);

    while (first_vertex != lca) {
        if (Receptacle(first_vertex) == lca_receptacle) {
            break;
        }
        first_vertex = OtherEnd(nodes[first_vertex.index].matched_edge, first_vertex);
        first_vertex = OtherEnd(nodes[first_vertex.index].minus_parent, first_vertex);
    }
    while (second_vertex != lca) {
        if (Receptacle(second_vertex) == lca_receptacle) {
            break;
        }
        second_vertex = OtherEnd(nodes[second_vertex.index].matched_edge, second_vertex);
        second_vertex = OtherEnd(nodes[second_vertex.index].minus_parent, second_vertex);
    }

    return {first_vertex, second_vertex};
}

void VzhuhSolver::UpdateCherryPath(NodeIndex lower_node, NodeIndex upper_node, PrimalUpdateRecord *record) {
    NodeIndex receptacle = Receptacle(upper_node);
    nodes[Receptacle(upper_node).index].receptacle_ = receptacle; // TODO does this line do anything?

    record->changed_sign.push_back(lower_node);
    while (lower_node != upper_node) {
        NodeIndex parent = OtherEnd(nodes[lower_node.index].matched_edge, lower_node);
        NodeIndex grandparent = OtherEnd(nodes[parent.index].minus_parent, parent);

        record->changed_sign.push_back(parent);
        record->changed_sign.push_back(grandparent);

        nodes[Receptacle(lower_node).index].receptacle_ = receptacle;
        nodes[Receptacle(parent).index].receptacle_ = receptacle;

        if (!nodes[parent.index].plus) {
            nodes[parent.index].plus = true;
            AddNeighborsToActionable(parent);
        }

        if (grandparent != upper_node) {
            nodes[grandparent.index].minus_parent = nodes[parent.index].minus_parent;
        }

        lower_node = grandparent;
    }
}

void VzhuhSolver::Augment(EdgeIndex edge_plus_plus, PrimalUpdateRecord *record) {
    num_trees_alive -= 2;

    NodeIndex head = Head(edge_plus_plus);
    NodeIndex tail = Tail(edge_plus_plus);

    if (params.verbose) {
        std::cout << "AUGMENT " << head.index << " " << tail.index << ", elementary ends: " << edges[edge_plus_plus.
                index].elementary_head.index << " " << edges[edge_plus_plus.index].elementary_tail.index <<
            ", is in zero slack set: " << edges[edge_plus_plus.index].is_in_zero_slack_set << std::endl;
    }

    std::vector<EdgeIndex> first_path = PathToRoot(head);
    std::vector<EdgeIndex> second_path = PathToRoot(tail);

    std::vector<EdgeIndex> path;
    path.reserve(first_path.size() + second_path.size() + 1);
    for (int i = static_cast<int>(first_path.size()) - 1; i >= 0; --i) {
        path.push_back(first_path[i]);
    }
    path.push_back(edge_plus_plus);
    for (const auto &edge : second_path) {
        path.push_back(edge);
    }

    TreeIndex first_tree = nodes[head.index].tree;
    TreeIndex second_tree = nodes[tail.index].tree;
    AugmentPath(path);
    ClearTree(first_tree, record);
    ClearTree(second_tree, record);
}

std::vector<VzhuhSolver::EdgeIndex> VzhuhSolver::PathToRoot(NodeIndex node_plus) {
    NodeIndex root = TopBlossom(trees[nodes[node_plus.index].tree.index].root);

    std::vector<EdgeIndex> path;
    while (node_plus != root) {
        path.push_back(nodes[node_plus.index].matched_edge);
        node_plus = OtherEnd(nodes[node_plus.index].matched_edge, node_plus);
        path.push_back(nodes[node_plus.index].minus_parent);
        node_plus = OtherEnd(nodes[node_plus.index].minus_parent, node_plus);
    }

    return path;
}

void VzhuhSolver::AugmentPath(const std::vector<EdgeIndex> &path) {
    bool match = true;
    for (EdgeIndex edge : path) {
        if (match) {
            MakeEdgeMatched(edge);
        } else {
            MakeEdgeUnmatched(edge);
        }
        match = !match;
    }
}

void VzhuhSolver::ClearTree(TreeIndex tree, PrimalUpdateRecord *record) {
    trees[tree.index].is_alive = false;

    // TODO don't clean up if last pair of trees, and amortize this

    for (NodeIndex node : trees[tree.index].tree_nodes) {
        if (nodes[node.index].blossom_parent || nodes[node.index].tree != tree) {
            continue;
        }

        record->changed_sign.push_back(node);
        nodes[node.index].tree = TreeIndex(-1);
        nodes[node.index].minus_parent = EdgeIndex(-1);
        nodes[node.index].receptacle_ = node;
        nodes[node.index].plus = false;

        // add edges to actionable
        // TODO make better
        AddNeighborsToActionable(node);
    }
}

void VzhuhSolver::UpdateQueues(const PrimalUpdateRecord &record) {
    // updates amortized slams and variables,
    // edge_heaps, node_heaps
    // old_plus, old_tree, old_blossom_parent for nodes

    for (NodeIndex node : record.changed_sign) {
        if (!nodes[node.index].is_alive) {
            continue;
        }
        if (params.verbose) {
            std::cout << "updating node " << node.index << std::endl;
        }

        if (nodes[node.index].old_blossom_parent == nodes[node.index].blossom_parent &&
            nodes[node.index].old_plus == nodes[node.index].plus &&
            nodes[node.index].old_tree == nodes[node.index].tree) {
            continue;
        }

        // update dual_var_quadrupled_amortized_
        int old_variable = DualVariableQuadrupled(node,
                                                  nodes[node.index].old_tree,
                                                  nodes[node.index].old_plus,
                                                  nodes[node.index].old_blossom_parent);
        int new_variable = DualVariableQuadrupled(node);
        nodes[node.index].dual_var_quadrupled_amortized_ += (old_variable - new_variable);

        // update amortized slams and queues
        NodeIndex old_top_node = node;
        if (nodes[node.index].old_blossom_parent) {
            old_top_node = nodes[node.index].old_blossom_parent;
        }
        NodeIndex top_node = node;
        if (nodes[node.index].blossom_parent) {
            top_node = nodes[node.index].blossom_parent;
        }
        int diff = 0;
        if (nodes[top_node.index].tree) {
            if (nodes[top_node.index].plus) {
                diff += trees[nodes[top_node.index].tree.index].dual_var_quadrupled;
            } else {
                diff -= trees[nodes[top_node.index].tree.index].dual_var_quadrupled;
            }
        }
        if (nodes[old_top_node.index].old_tree) {
            if (nodes[old_top_node.index].old_plus) {
                diff -= trees[nodes[old_top_node.index].old_tree.index].dual_var_quadrupled;
            } else {
                diff += trees[nodes[old_top_node.index].old_tree.index].dual_var_quadrupled;
            }
        }
        if (diff != 0) {
            for (EdgeIndex edge : NeighborsWLoops(node)) {
                RemoveEdgeFromQueue(edge);
                edges[edge.index].slam_quadrupled_amortized_ += diff;
                AddEdgeToQueue(edge);
            }
        }
        if (nodes[node.index].old_plus != nodes[node.index].plus || nodes[node.index].old_tree != nodes[node.index].
            tree) {
            for (EdgeIndex edge : NonLoopNeighbors(node)) {
                AddEdgeToQueue(edge);
            }
        }

        // update minus_blossoms queues
        if (!IsElementary(node)) {
            if (nodes[node.index].old_plus != nodes[node.index].plus || nodes[node.index].old_tree != nodes[node.index].
                tree) {
                // the status might have changed
                if (nodes[node.index].old_tree && !nodes[node.index].old_plus) {
                    // were a minus blossom
                    RemoveNodeFromQueue(node);
                }
                if (nodes[node.index].tree && !nodes[node.index].plus) {
                    // became a minus blossom
                    AddNodeToQueue(node, trees[nodes[node.index].tree.index].minus_blossoms);
                }
            }
        }

        // update old_plus, old_tree, old_blossom_parent
        nodes[node.index].old_blossom_parent = nodes[node.index].blossom_parent;
        nodes[node.index].old_plus = nodes[node.index].plus;
        nodes[node.index].old_tree = nodes[node.index].tree;
    }
}

void VzhuhSolver::AddNeighborsToActionable(NodeIndex node) {
    std::vector<EdgeIndex> grandchild_neighbors = NonLoopZeroSlackNeighbors(node);
    for (EdgeIndex neighbor_edge : grandchild_neighbors) {
        actionable_edges.push(neighbor_edge);
    }
}

std::vector<std::vector<VzhuhSolver::NodeIndex> > VzhuhSolver::
OrganizeBlossomChildren(const PrimalUpdateRecord &record) {
    std::vector<std::vector<NodeIndex> > result;

    ++nodes_label_cnt;
    int label_zero = nodes_label_cnt;

    for (NodeIndex node : record.changed_sign) {
        if (!nodes[node.index].tree) {
            continue;
        }
        if (nodes[node.index].label >= label_zero) {
            // already seen this node
            continue;
        }

        NodeIndex receptacle = Receptacle(node);
        if (nodes[receptacle.index].label < label_zero) {
            nodes[receptacle.index].label = nodes_label_cnt;
            ++nodes_label_cnt;
            result.push_back({receptacle});
        }
        if (nodes[node.index].label < label_zero) {
            result[nodes[receptacle.index].label - label_zero].push_back(node);
            nodes[node.index].label = nodes[receptacle.index].label;
        }
    }

    // remove singleton lists
    for (int i = 0; i < static_cast<int>(result.size()); ++i) {
        if (result[i].size() <= 1) {
            result[i] = std::move(result.back());
            result.pop_back();
            --i;
        }
    }

    return result;
}

void VzhuhSolver::Shrink(std::vector<NodeIndex> &children) {
    if (params.verbose) {
        std::cout << "SHRINK ";
        for (NodeIndex node : children) {
            std::cout << node.index << " ";
        }
        std::cout << std::endl;
    }

    if (children.size() % 2 == 0) {
        throw std::runtime_error("In Node: blossom has to have an odd number of vertices");
    }

    int new_index = static_cast<int>(nodes.size());
    NodeIndex receptacle = Receptacle(children.front());

    nodes.emplace_back(new_index);
    nodes.back().plus = true;
    nodes.back().old_plus = true;
    nodes.back().blossom_children = std::move(children);
    nodes.back().matched_edge = nodes[receptacle.index].matched_edge;
    nodes.back().tree = nodes[receptacle.index].tree;
    nodes.back().old_tree = nodes.back().tree;
    nodes.back().dual_var_quadrupled_amortized_ = -trees[nodes.back().tree.index].dual_var_quadrupled;

    // update blossom_parent of the children, set labels to empty
    for (NodeIndex child : nodes.back().blossom_children) {
        if (nodes[child.index].blossom_parent) {
            throw std::runtime_error("In Node: some child node already has a parent");
        }
        nodes[child.index].blossom_parent = NodeIndex(new_index);
    }

    // add the new blossom to the list of vertices in the tree
    trees[nodes.back().tree.index].tree_nodes.push_back(NodeIndex(new_index));
}

void VzhuhSolver::ValidateQueues() {
    // checks that the queues are in a correct state, throws if not

    // every queue must hold only the edges that belong to this queue

    for (TreeIndex tree : alive_trees) {
        if (!trees[tree.index].is_alive) {
            continue;
        }

        // check plus empty
        edge_heaps[trees[tree.index].plus_empty_edges]->ValidateHeap("plus empty");
        for (auto heap_node : edge_heaps[trees[tree.index].plus_empty_edges]->heap_) {
            EdgeIndex edge = heap_node.value;

            NodeIndex first = Head(edge); // plus
            NodeIndex second = Tail(edge); // empty
            if (nodes[first.index].tree != tree) {
                std::swap(first, second);
            }
            if (nodes[first.index].tree != tree || !nodes[first.index].plus) {
                throw std::runtime_error("Incorrect plus empty queue");
            }
            if (nodes[second.index].tree) {
                throw std::runtime_error("Incorrect plus empty queue");
            }
        }

        // check plus plus internal
        edge_heaps[trees[tree.index].plus_plus_internal_edges]->ValidateHeap("+ + int");
        for (auto heap_node : edge_heaps[trees[tree.index].plus_plus_internal_edges]->heap_) {
            EdgeIndex edge = heap_node.value;

            NodeIndex first = Head(edge);
            NodeIndex second = Tail(edge);

            if (nodes[first.index].tree != tree || !nodes[first.index].plus) {
                throw std::runtime_error("Incorrect plus plus internal queue");
            }
            if (nodes[second.index].tree != tree || !nodes[second.index].plus) {
                throw std::runtime_error("Incorrect plus plus internal queue");
            }
        }

        // check minus blossoms
        node_heaps[trees[tree.index].minus_blossoms]->ValidateHeap("minus blossoms");
        for (auto heap_node : node_heaps[trees[tree.index].minus_blossoms]->heap_) {
            NodeIndex node = heap_node.value;

            if (IsElementary(node) || (nodes[node.index].tree != tree) || nodes[node.index].plus) {
                throw std::runtime_error("Incorrect minus blossom queue");
            }
        }

        // check plus plus external
        for (auto [other_tree, queue_index] : trees[tree.index].pq_plus_plus) {
            if (!trees[other_tree.index].is_alive) {
                continue;
            }
            for (auto heap_node : edge_heaps[queue_index]->heap_) {
                EdgeIndex edge = heap_node.value;

                NodeIndex first = Head(edge); // in this tree
                NodeIndex second = Tail(edge); // in the other tree
                if (nodes[first.index].tree != tree) {
                    std::swap(first, second);
                }
                if (nodes[first.index].tree != tree || !nodes[first.index].plus) {
                    throw std::runtime_error("Incorrect plus plus queue");
                }
                if (nodes[second.index].tree != other_tree || !nodes[second.index].plus) {
                    throw std::runtime_error("Incorrect plus plus queue");
                }
            }
            edge_heaps[queue_index]->ValidateHeap("+ + external");
        }

        // check plus minus external
        for (auto [other_tree, queue_index] : trees[tree.index].pq_plus_minus) {
            if (!trees[other_tree.index].is_alive) {
                continue;
            }
            edge_heaps[queue_index]->ValidateHeap("+ - ext");
            for (auto heap_node : edge_heaps[queue_index]->heap_) {
                EdgeIndex edge = heap_node.value;

                NodeIndex first = Head(edge); // in this tree
                NodeIndex second = Tail(edge); // in the other tree
                if (nodes[first.index].tree != tree) {
                    std::swap(first, second);
                }
                if (nodes[first.index].tree != tree || !nodes[first.index].plus) {
                    throw std::runtime_error("Incorrect plus minus queue");
                }
                if (nodes[second.index].tree != other_tree || nodes[second.index].plus) {
                    throw std::runtime_error("Incorrect plus minus queue");
                }
            }
        }

        // check minus plus external
        for (auto [other_tree, queue_index] : trees[tree.index].pq_minus_plus) {
            if (!trees[other_tree.index].is_alive) {
                continue;
            }
            edge_heaps[queue_index]->ValidateHeap("- + ext");
            for (auto heap_node : edge_heaps[queue_index]->heap_) {
                EdgeIndex edge = heap_node.value;

                NodeIndex first = Head(edge); // in this tree
                NodeIndex second = Tail(edge); // in the other tree
                if (nodes[first.index].tree != tree) {
                    std::swap(first, second);
                }
                if (nodes[first.index].tree != tree || nodes[first.index].plus) {
                    throw std::runtime_error("Incorrect minus plus queue");
                }
                if (nodes[second.index].tree != other_tree || !nodes[second.index].plus) {
                    throw std::runtime_error("Incorrect minus plus queue");
                }
            }
        }
    }
}

std::vector<int> VzhuhSolver::NodeVariables() const {
    std::vector<int> result;
    for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
        if (nodes[i].is_alive) {
            result.push_back(DualVariableQuadrupled(NodeIndex(i)));
        } else {
            result.push_back(0);
        }
    }
    return result;
}

std::vector<int> VzhuhSolver::EdgeSlams() {
    std::vector<int> result;
    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
        result.push_back(SlamQuadrupled(EdgeIndex(i)));
    }
    return result;
}

void VzhuhSolver::ValidateZeroSlackAdjList() {
    // check that every non-loop zero slack edge is in the list
    // check that there are no duplicates

    std::vector<std::unordered_set<int>> indices (zero_slack_adj_list.size());
    for (int i = 0; i < static_cast<int>(zero_slack_adj_list.size()); ++i) {
        for (EdgeIndex edge : zero_slack_adj_list[i]) {
            indices[i].insert(edge.index);
        }
        if (indices[i].size() != zero_slack_adj_list[i].size()) {
            std::cout << "vertex " << i << std::endl;
            throw std::runtime_error("zero_slack_adj_list contains duplicate edges");
        }
    }

    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
        EdgeIndex edge(i);
        if (Head(edge) == Tail(edge) || SlamQuadrupled(edge) > 0) {
            continue;
        }
        if (!indices[edges[i].elementary_head.index].contains(edge.index) || !indices[edges[i].elementary_tail.index].contains(edge.index)) {
            std::cout << "edge " << i << ", head/tail: " << edges[i].elementary_head.index << " " << edges[i].elementary_tail.index << std::endl;
            throw std::runtime_error("zero_slack_adj_list does not contain the edge");
        }
    }
}

void VzhuhSolver::ValidateEvenOddPaths() {
    for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
        NodeIndex node(i);
        if (!nodes[i].is_alive || nodes[i].blossom_parent) {
            continue;
        }

        NodeIndex receptacle = Receptacle(node);

        auto even_path = EvenPathToReceptacle(node);
        if (even_path.size() % 2 != 0) {
            throw std::runtime_error("ValidateEvenOddPaths: even path is not even");
        }
        NodeIndex cur_node = node;
        for (EdgeIndex edge : even_path) {
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
            for (EdgeIndex edge : odd_path) {
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
        trees[alive_trees[i].index].dual_var_quadrupled += deltas[i];
    }

    UpdateZeroSlackSetAndActionable();

    if (params.verbose) {
        std::cout << "after dual update:" << std::endl;
        PrintGraph();
    }
    if (params.debug) {
        ValidateZeroSlackAdjList();
    }

    return !std::all_of(deltas.begin(), deltas.end(), [](auto x){ return x == 0; });
}

void VzhuhSolver::UpdateAliveTreesList() {
    for (int i = 0; i < static_cast<int>(alive_trees.size()); ++i) {
        if (!trees[alive_trees[i].index].is_alive) {
            alive_trees[i] = alive_trees.back();
            alive_trees.pop_back();
            --i;
        }
    }
}

std::vector<int> VzhuhSolver::VariableDeltas() {
    DualConstraints dual_constraints = GetDualConstraints();

    for (auto b : dual_constraints.upper_bound) {
        if (b < 0) {
            throw std::runtime_error("b < 0");
        }
    }
    for (auto c : dual_constraints.plus_plus_constraints) {
        for (auto [ind, sl] : c) {
            if (sl < 0) {
                throw std::runtime_error("constraint plus plus < 0");
            }
        }
    }
    for (auto c : dual_constraints.plus_minus_constraints) {
        for (auto [ind, sl] : c) {
            if (sl < 0) {
                throw std::runtime_error("constraint plus minus < 0");
            }
        }
    }

    std::vector<std::vector<int> > connected_components = ConnectedComponentsTreeTree(dual_constraints);
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
        TreeIndex tree = alive_trees[i];

        int plus_empty = PlusEmptySlack(tree);
        int plus_plus_internal = PlusPlusInternalSlack(tree);
        int blossom_var = MinMinusBlossomVariable(tree);

        if (plus_empty < 0) {
            throw std::runtime_error("In VariableDeltas: plus_empty is negative");
        }
        if (plus_plus_internal < 0) {
            throw std::runtime_error("In VariableDeltas: plus_plus_internal is not positive");
        }
        if (blossom_var < 0) {
            throw std::runtime_error("In VariableDeltas: blossom_var is not positive");
        }
        if (plus_plus_internal < INT32_MAX && plus_plus_internal % 2 != 0) {
            throw std::runtime_error("PlusPlusInternalSlack is not divisible by 2");
        }

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
        trees[alive_trees[i].index].alive_index = i;
    }

    // get plus_plus_constraints and plus_minus_constraints
    std::vector plus_plus_constraints(alive_trees.size(), std::vector<std::pair<int, int> >());
    std::vector plus_minus_constraints(alive_trees.size(), std::vector<std::pair<int, int> >());
    for (int i = 0; i < static_cast<int>(alive_trees.size()); ++i) {
        TreeIndex tree = alive_trees[i];

        std::vector<std::pair<TreeIndex, int> > plus_plus_slacks = PlusPlusExternalSlacks(tree);
        std::vector<std::pair<TreeIndex, int> > plus_minus_slacks = PlusMinusExternalSlacks(tree);
        for (auto [other_tree, slack] : plus_plus_slacks) {
            plus_plus_constraints[i].emplace_back(trees[other_tree.index].alive_index, slack);
        }
        for (auto [other_tree, slack] : plus_minus_slacks) {
            plus_minus_constraints[i].emplace_back(trees[other_tree.index].alive_index, slack);
        }
    }

    return std::move(DualConstraints{
        std::move(upper_bound), std::move(plus_plus_constraints), std::move(plus_minus_constraints)
    });
}

void VzhuhSolver::UpdateZeroSlackSetAndActionable() {
    /*for (int i = 0; i < static_cast<int>(zero_slack_adj_list.size()); ++i) {
        for (int j = 0; j < static_cast<int>(zero_slack_adj_list[i].size()); ++j) {
            if (Head(zero_slack_adj_list[i][j]) == Tail(zero_slack_adj_list[i][j])) {
                continue;
            }
            if (SlamQuadrupled(zero_slack_adj_list[i][j]) != 0) {
                edges[zero_slack_adj_list[i][j].index].is_in_zero_slack_set = false;
                zero_slack_adj_list[i][j] = zero_slack_adj_list[i].back();
                zero_slack_adj_list[i].pop_back();
                --j;
            }
        }
    }

    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
        EdgeIndex edge(i);
        if (SlamQuadrupled(edge) == 0 && Head(edge) != Tail(edge) && !edges[edge.index].is_in_zero_slack_set) {
            zero_slack_adj_list[edges[edge.index].elementary_head.index].push_back(edge);
            zero_slack_adj_list[edges[edge.index].elementary_tail.index].push_back(edge);
            edges[edge.index].is_in_zero_slack_set = true;
            actionable_edges.push(edge);
        }
    }*/

    for (TreeIndex tree : alive_trees) {
        if (params.verbose) {
            std::cout << "updating zero slack set next to tree " << tree.index << std::endl;
        }

        AddZeroSlackEdgesFromQueue(trees[tree.index].plus_empty_edges, true);
        CleanLoopsFromQueueTop(tree);
        AddZeroSlackEdgesFromQueue(trees[tree.index].plus_plus_internal_edges, true);
        for (auto [other_tree, queue_index] : trees[tree.index].pq_plus_plus) {
            if (other_tree.index > tree.index) {
                // to avoid double work
                AddZeroSlackEdgesFromQueue(queue_index, true);
            }
        }
        for (auto [_, queue_index] : trees[tree.index].pq_plus_minus) {
            AddZeroSlackEdgesFromQueue(queue_index, false);
        }
    }
}

void VzhuhSolver::AddZeroSlackEdgesFromQueue(int queue_index, bool add_to_actionable) {
    if (edge_heaps[queue_index]->Empty()) {
        return;
    }

    EdgeIndex top = edge_heaps[queue_index]->Top();
    if (SlamQuadrupled(top) > 0) {
        return;
    }

    if (SlamQuadrupled(top) < 0) {
        throw std::runtime_error("AddZeroSlackEdgesFromQueue: top with negative slam");
    }

    std::vector<EdgeIndex> zero_slack_edges = edge_heaps[queue_index]->ElementsEqualToTop();
    for (EdgeIndex edge : zero_slack_edges) {
        if (!edges[edge.index].is_in_zero_slack_set) {
            edges[edge.index].is_in_zero_slack_set = true;

            zero_slack_adj_list[edges[edge.index].elementary_head.index].push_back(edge);
            zero_slack_adj_list[edges[edge.index].elementary_tail.index].push_back(edge);
        }

        if (add_to_actionable) {
            actionable_edges.push(edge);
        }
    }
}

void VzhuhSolver::CleanLoopsFromQueueTop(TreeIndex tree) {
    int queue_idx = trees[tree.index].plus_plus_internal_edges;

    while (!edge_heaps[queue_idx]->Empty()) {
        EdgeIndex top = edge_heaps[queue_idx]->Top();
        if (Head(top) == Tail(top)) {
            RemoveEdgeFromQueue(top);
        } else {
            break;
        }
    }
}

bool VzhuhSolver::IsElementary(NodeIndex node) const {
    return node.index < num_vertices_elementary;
}

VzhuhSolver::NodeIndex VzhuhSolver::TopBlossom(NodeIndex node) const {
    while (nodes[node.index].blossom_parent) {
        node = nodes[node.index].blossom_parent;
    }
    return node;
}

VzhuhSolver::NodeIndex VzhuhSolver::Receptacle(NodeIndex node) {
    NodeIndex grandparent = node;
    while (nodes[grandparent.index].receptacle_ != grandparent) {
        grandparent = nodes[grandparent.index].receptacle_;
    }

    NodeIndex cur_node = node;
    while (nodes[cur_node.index].receptacle_ != grandparent) {
        NodeIndex next_node = nodes[cur_node.index].receptacle_;
        nodes[cur_node.index].receptacle_ = grandparent;
        cur_node = next_node;
    }

    return grandparent;
}

int VzhuhSolver::DualVariableQuadrupled(NodeIndex node) const {
    return DualVariableQuadrupled(node,
                                  nodes[node.index].tree,
                                  nodes[node.index].plus,
                                  nodes[node.index].blossom_parent);
}

int VzhuhSolver::DualVariableQuadrupled(NodeIndex node, TreeIndex tree, bool plus, NodeIndex blossom_parent) const {
    if (!nodes[node.index].is_alive) {
        throw std::runtime_error("DualVariableQuadrupled: node is not alive");
    }

    if (!tree || blossom_parent) {
        return nodes[node.index].dual_var_quadrupled_amortized_;
    }
    if (plus) {
        return nodes[node.index].dual_var_quadrupled_amortized_ + trees[tree.index].dual_var_quadrupled;
    }
    return nodes[node.index].dual_var_quadrupled_amortized_ - trees[tree.index].dual_var_quadrupled;
}

std::vector<VzhuhSolver::EdgeIndex> VzhuhSolver::NonLoopNeighbors(NodeIndex node) {
    // TODO make it a field and compute lazily?

    if (IsElementary(node)) {
        return adj_list[node.index];
    }

    std::vector<NodeIndex> elementary_descendants = ElementaryBlossomDescendants(node);

    // label the elementary descendants
    ++nodes_label_cnt;
    for (NodeIndex descendant : elementary_descendants) {
        nodes[descendant.index].label = nodes_label_cnt;
    }

    std::vector<EdgeIndex> neighbors;
    for (NodeIndex descendant : elementary_descendants) {
        for (EdgeIndex edge : adj_list[descendant.index]) {
            if (nodes[OtherElementaryEnd(edge, descendant).index].label != nodes_label_cnt) {
                neighbors.push_back(edge);
            }
        }
    }

    return neighbors;
}

std::vector<VzhuhSolver::EdgeIndex> VzhuhSolver::NeighborsWLoops(NodeIndex node) {
    if (IsElementary(node)) {
        return adj_list[node.index];
    }

    std::vector<NodeIndex> elementary_descendants = ElementaryBlossomDescendants(node);

    std::vector<EdgeIndex> neighbors;
    for (NodeIndex descendant : elementary_descendants) {
        for (EdgeIndex edge : adj_list[descendant.index]) {
            neighbors.push_back(edge);
        }
    }

    return neighbors;
}

std::vector<VzhuhSolver::EdgeIndex> VzhuhSolver::NonLoopZeroSlackNeighbors(NodeIndex node) {
    if (IsElementary(node)) {
        return zero_slack_adj_list[node.index];
    }

    std::vector<NodeIndex> elementary_descendants = ElementaryBlossomDescendants(node);

    // label the elementary descendants
    ++nodes_label_cnt;
    for (NodeIndex descendant : elementary_descendants) {
        nodes[descendant.index].label = nodes_label_cnt;
    }

    std::vector<EdgeIndex> neighbors;
    for (NodeIndex descendant : elementary_descendants) {
        for (EdgeIndex edge : zero_slack_adj_list[descendant.index]) {
            if (nodes[OtherElementaryEnd(edge, descendant).index].label != nodes_label_cnt) {
                neighbors.push_back(edge);
            }
        }
    }

    return neighbors;
}

std::vector<VzhuhSolver::NodeIndex> VzhuhSolver::ElementaryBlossomDescendants(NodeIndex node) const {
    if (IsElementary(node)) {
        return {node};
    }

    std::vector<NodeIndex> elementary_descendants;
    std::queue<NodeIndex> queue;
    queue.push(node);
    while (!queue.empty()) {
        NodeIndex cur = queue.front();
        queue.pop();
        for (NodeIndex child : nodes[cur.index].blossom_children) {
            if (IsElementary(child)) {
                elementary_descendants.push_back(child);
            } else {
                queue.push(child);
            }
        }
    }
    return elementary_descendants;
}

int VzhuhSolver::SlamQuadrupled(EdgeIndex edge) {
    NodeIndex head = Head(edge);
    NodeIndex tail = Tail(edge);

    int slam = edges[edge.index].slam_quadrupled_amortized_;

    if (nodes[head.index].tree) {
        if (nodes[head.index].plus) {
            slam -= trees[nodes[head.index].tree.index].dual_var_quadrupled;
        } else {
            slam += trees[nodes[head.index].tree.index].dual_var_quadrupled;
        }
    }
    if (nodes[tail.index].tree) {
        if (nodes[tail.index].plus) {
            slam -= trees[nodes[tail.index].tree.index].dual_var_quadrupled;
        } else {
            slam += trees[nodes[tail.index].tree.index].dual_var_quadrupled;
        }
    }

    return slam;
}

int VzhuhSolver::OldSlamQuadrupled(EdgeIndex edge) {
    NodeIndex head = Head(edge);
    if (nodes[head.index].old_blossom_parent) {
        head = nodes[head.index].old_blossom_parent;
    }
    NodeIndex tail = Tail(edge);
    if (nodes[tail.index].old_blossom_parent) {
        tail = nodes[tail.index].old_blossom_parent;
    }

    int slam = edges[edge.index].slam_quadrupled_amortized_;

    if (nodes[head.index].old_tree) {
        if (nodes[head.index].old_plus) {
            slam -= trees[nodes[head.index].old_tree.index].dual_var_quadrupled;
        } else {
            slam += trees[nodes[head.index].old_tree.index].dual_var_quadrupled;
        }
    }
    if (nodes[tail.index].old_tree) {
        if (nodes[tail.index].old_plus) {
            slam -= trees[nodes[tail.index].old_tree.index].dual_var_quadrupled;
        } else {
            slam += trees[nodes[tail.index].old_tree.index].dual_var_quadrupled;
        }
    }

    return slam;
}

VzhuhSolver::NodeIndex VzhuhSolver::OtherEnd(EdgeIndex edge, NodeIndex node) {
    if (node == Tail(edge)) {
        return Head(edge);
    }
    if (node == Head(edge)) {
        return Tail(edge);
    }
    if (Head(edge) != Tail(edge)) {
        throw std::runtime_error("In OtherEnd");
    }
    return node;
}

VzhuhSolver::NodeIndex VzhuhSolver::OtherElementaryEnd(EdgeIndex edge, NodeIndex node) const {
    if (node == edges[edge.index].elementary_head) {
        return edges[edge.index].elementary_tail;
    }
    return edges[edge.index].elementary_head;
}

VzhuhSolver::NodeIndex VzhuhSolver::Head(EdgeIndex edge) {
    while (nodes[edges[edge.index].head.index].blossom_parent) {
        // throw std::runtime_error("Head is not amortized now");
        edges[edge.index].head = nodes[edges[edge.index].head.index].blossom_parent;
    }
    return edges[edge.index].head;
}

VzhuhSolver::NodeIndex VzhuhSolver::Tail(EdgeIndex edge) {
    while (nodes[edges[edge.index].tail.index].blossom_parent) {
        // throw std::runtime_error("tail is not amortized now");
        edges[edge.index].tail = nodes[edges[edge.index].tail.index].blossom_parent;
    }
    return edges[edge.index].tail;
}

VzhuhSolver::NodeIndex VzhuhSolver::PlusPlusLCA(NodeIndex first_vertex, NodeIndex second_vertex) {
    ++nodes_label_cnt;
    first_vertex = Receptacle(first_vertex);
    second_vertex = Receptacle(second_vertex);

    nodes[first_vertex.index].label = nodes_label_cnt;
    nodes[second_vertex.index].label = nodes_label_cnt;

    while (first_vertex != second_vertex) {
        if (nodes[first_vertex.index].matched_edge) {
            first_vertex = OtherEnd(nodes[first_vertex.index].matched_edge, first_vertex);
            first_vertex = OtherEnd(nodes[first_vertex.index].minus_parent, first_vertex);
            first_vertex = Receptacle(first_vertex);
            if (nodes[first_vertex.index].label == nodes_label_cnt) {
                return first_vertex;
            }
            nodes[first_vertex.index].label = nodes_label_cnt;
        }

        if (nodes[second_vertex.index].matched_edge) {
            second_vertex = OtherEnd(nodes[second_vertex.index].matched_edge, second_vertex);
            second_vertex = OtherEnd(nodes[second_vertex.index].minus_parent, second_vertex);
            second_vertex = Receptacle(second_vertex);
            if (nodes[second_vertex.index].label == nodes_label_cnt) {
                return second_vertex;
            }
            nodes[second_vertex.index].label = nodes_label_cnt;
        }
    }

    return first_vertex;
}

void VzhuhSolver::MakeEdgeMatched(EdgeIndex edge) {
    if (Head(edge) == Tail(edge)) {
        throw std::runtime_error("In MakeEdgeMatched: trying to match a loop");
    }
    nodes[Head(edge).index].matched_edge = edge;
    nodes[Tail(edge).index].matched_edge = edge;
    edges[edge.index].matched = true;
}

void VzhuhSolver::MakeEdgeUnmatched(EdgeIndex edge) {
    edges[edge.index].matched = false;
}

void VzhuhSolver::AddEdgeToQueue(EdgeIndex edge) {
    // adds the edge to the right queue

    NodeIndex head = Head(edge);
    NodeIndex tail = Tail(edge);

    // same tree or (empty, empty)
    if (nodes[head.index].tree == nodes[tail.index].tree) {
        if (!nodes[head.index].tree) {
            // (empty, empty)
            RemoveEdgeFromQueue(edge);
            return;
        }
        if (nodes[head.index].plus && nodes[tail.index].plus) {
            // (+, +)
            AddEdgeToThisQueue(edge, trees[nodes[head.index].tree.index].plus_plus_internal_edges);
            return;
        }

        // (+, -), (-, +), (-, -)
        RemoveEdgeFromQueue(edge);
        return;
    }

    // different trees, make head always in some tree
    if (!nodes[head.index].tree) {
        std::swap(head, tail);
    }

    // (+, empty)
    if (nodes[head.index].plus && !nodes[tail.index].tree) {
        AddEdgeToThisQueue(edge, trees[nodes[head.index].tree.index].plus_empty_edges);
        return;
    }

    // (+, +)
    if (nodes[head.index].plus && nodes[tail.index].plus) {
        AddPQPlusPlus(nodes[head.index].tree, nodes[tail.index].tree, edge);
        return;
    }

    // (+, -)
    if (nodes[head.index].plus && !nodes[tail.index].plus) {
        AddPQPlusMinus(nodes[head.index].tree, nodes[tail.index].tree, edge);
        return;
    }

    // (-, +)
    if (!nodes[head.index].plus && nodes[tail.index].plus) {
        AddPQPlusMinus(nodes[tail.index].tree, nodes[head.index].tree, edge);
        return;
    }
}

void VzhuhSolver::AddEdgeToThisQueue(EdgeIndex edge, int queue_index) {
    if (edges[edge.index].queue_index == queue_index) {
        return;
    }

    if (edges[edge.index].queue_index >= 0) {
        edge_heaps[edges[edge.index].queue_index]->Erase(edges[edge.index].handle);
    }
    edges[edge.index].handle = edge_heaps[queue_index]->Push(edge);
    edges[edge.index].queue_index = queue_index;
}

void VzhuhSolver::RemoveEdgeFromQueue(EdgeIndex edge) {
    if (edges[edge.index].queue_index >= 0) {
        edge_heaps[edges[edge.index].queue_index]->Erase(edges[edge.index].handle);
        edges[edge.index].queue_index = -1;
    }
}

void VzhuhSolver::AddNodeToQueue(NodeIndex node, int queue_index) {
    if (nodes[node.index].queue_index == queue_index) {
        return;
    }

    if (nodes[node.index].queue_index >= 0) {
        node_heaps[nodes[node.index].queue_index]->Erase(nodes[node.index].handle);
    }
    nodes[node.index].handle = node_heaps[queue_index]->Push(node);
    nodes[node.index].queue_index = queue_index;
}

void VzhuhSolver::RemoveNodeFromQueue(NodeIndex node) {
    if (nodes[node.index].queue_index >= 0) {
        node_heaps[nodes[node.index].queue_index]->Erase(nodes[node.index].handle);
        nodes[node.index].queue_index = -1;
    }
}

void VzhuhSolver::AddPQPlusPlus(TreeIndex first, TreeIndex second, EdgeIndex edge) {
    int queue_index = TreeTreeQueueIndex(second, &trees[first.index].pq_plus_plus);
    if (queue_index > 0) {
        AddEdgeToThisQueue(edge, queue_index);
    } else {
        edge_heaps.emplace_back(std::make_unique<EdgeHeap>(params.heap_arity, EdgeComparator{this}));
        queue_index = static_cast<int>(edge_heaps.size()) - 1;
        AddEdgeToThisQueue(edge, queue_index);

        trees[first.index].pq_plus_plus.emplace_back(second, queue_index);
        trees[second.index].pq_plus_plus.emplace_back(first, queue_index);
    }
}

void VzhuhSolver::AddPQPlusMinus(TreeIndex tree_plus, TreeIndex tree_minus, EdgeIndex edge) {
    int queue_index = TreeTreeQueueIndex(tree_minus, &trees[tree_plus.index].pq_plus_minus);
    if (queue_index > 0) {
        AddEdgeToThisQueue(edge, queue_index);
    } else {
        edge_heaps.emplace_back(std::make_unique<EdgeHeap>(params.heap_arity, EdgeComparator{this}));
        queue_index = static_cast<int>(edge_heaps.size()) - 1;
        AddEdgeToThisQueue(edge, queue_index);

        trees[tree_plus.index].pq_plus_minus.emplace_back(tree_minus, queue_index);
        trees[tree_minus.index].pq_minus_plus.emplace_back(tree_plus, queue_index);
    }
}

int VzhuhSolver::TreeTreeQueueIndex(TreeIndex other_tree,
                                    std::vector<std::pair<TreeIndex, int> > *tree_neighbors) {
    // may delete dead trees from the tree_neighbors
    for (int i = 0; i < static_cast<int>(tree_neighbors->size()); ++i) {
        if ((*tree_neighbors)[i].first == other_tree) {
            return (*tree_neighbors)[i].second;
        }
        if (!trees[(*tree_neighbors)[i].first.index].is_alive) {
            (*tree_neighbors)[i] = tree_neighbors->back();
            tree_neighbors->pop_back();
            --i;
        }
    }
    return -1;
}

VzhuhSolver::EdgeIndex VzhuhSolver::MinPlusEmptyEdge(int queue_index) {
    if (!edge_heaps[queue_index]->Empty()) {
        EdgeIndex edge = edge_heaps[queue_index]->Top();

        NodeIndex head = Head(edge); // plus
        NodeIndex tail = Tail(edge); // empty
        if (!nodes[head.index].tree) {
            std::swap(head, tail);
        }
        if (head == tail) {
            throw std::runtime_error("MinPlusEmptyEdge is a loop");
        }
        if (nodes[head.index].plus && !nodes[tail.index].tree) {
            return edge;
        }
        throw std::runtime_error("MinPlusEmptyEdge: incorrect top");
    }
    return EdgeIndex(-1);
}

VzhuhSolver::EdgeIndex VzhuhSolver::MinPlusPlusInternalEdge(int queue_index) {
    while (!edge_heaps[queue_index]->Empty()) {
        EdgeIndex edge = edge_heaps[queue_index]->Top();

        NodeIndex head = Head(edge);
        NodeIndex tail = Tail(edge);

        if (!nodes[head.index].plus || !nodes[tail.index].plus || (nodes[head.index].tree != nodes[tail.index].tree)) {
            throw std::runtime_error("MinPlusPlusInternalEdge: incorrect top");
        }

        if (head != tail) {
            return edge;
        }
        edge_heaps[queue_index]->Pop();
        edges[edge.index].queue_index = -1;
    }
    return EdgeIndex(-1);
}

VzhuhSolver::EdgeIndex VzhuhSolver::MinPlusPlusExternalEdge(int queue_index) {
    if (!edge_heaps[queue_index]->Empty()) {
        EdgeIndex edge = edge_heaps[queue_index]->Top();

        NodeIndex head = Head(edge);
        NodeIndex tail = Tail(edge);
        if (nodes[head.index].plus && nodes[tail.index].plus && (nodes[head.index].tree != nodes[tail.index].tree)) {
            return edge;
        }
        throw std::runtime_error("MinPlusPlusExternalEdge: incorrect top");
    }
    return EdgeIndex(-1);
}

VzhuhSolver::EdgeIndex VzhuhSolver::MinPlusMinusExternalEdge(int queue_index) {
    if (!edge_heaps[queue_index]->Empty()) {
        EdgeIndex edge = edge_heaps[queue_index]->Top();

        NodeIndex head = Head(edge); // plus
        NodeIndex tail = Tail(edge); // minus
        if (!nodes[head.index].plus) {
            std::swap(head, tail);
        }
        if (nodes[head.index].plus && !nodes[tail.index].plus && (nodes[head.index].tree != nodes[tail.index].tree) &&
            nodes[tail.index].tree) {
            return edge;
        }
        throw std::runtime_error("MinPlusMinusExternalEdge: incorrect top");
    }
    return EdgeIndex(-1);
}

VzhuhSolver::NodeIndex VzhuhSolver::MinMinusBlossom(int queue_index) const {
    if (!node_heaps[queue_index]->Empty()) {
        NodeIndex node = node_heaps[queue_index]->Top();
        if (!nodes[node.index].is_alive || !nodes[node.index].tree || nodes[node.index].plus || nodes[node.index].
            blossom_parent) {
            throw std::runtime_error("In MinMinusBlossom: incorrect top");
        }
        return node;
    }
    return NodeIndex(-1);
}

VzhuhSolver::NodeIndex VzhuhSolver::PopExpandableBlossom(TreeIndex tree) {
    NodeIndex node = MinMinusBlossom(trees[tree.index].minus_blossoms);
    if (node) {
        if (DualVariableQuadrupled(node) == 0) {
            node_heaps[trees[tree.index].minus_blossoms]->Pop();
            return node;
        }
    }
    return NodeIndex(-1);
}

int VzhuhSolver::PlusEmptySlack(TreeIndex tree) {
    EdgeIndex edge = MinPlusEmptyEdge(trees[tree.index].plus_empty_edges);
    if (edge) {
        return SlamQuadrupled(edge);
    }
    return INT32_MAX;
}

int VzhuhSolver::PlusPlusInternalSlack(TreeIndex tree) {
    EdgeIndex edge = MinPlusPlusInternalEdge(trees[tree.index].plus_plus_internal_edges);
    if (edge) {
        return SlamQuadrupled(edge);
    }
    return INT32_MAX;
}

int VzhuhSolver::MinMinusBlossomVariable(TreeIndex tree) {
    NodeIndex node = MinMinusBlossom(trees[tree.index].minus_blossoms);
    if (node) {
        return DualVariableQuadrupled(node);
    }
    return INT32_MAX;
}

std::vector<std::pair<VzhuhSolver::TreeIndex, int> > VzhuhSolver::PlusPlusExternalSlacks(TreeIndex tree) {
    std::vector<std::pair<TreeIndex, int> > result;
    result.reserve(trees[tree.index].pq_plus_plus.size());

    if (params.verbose) {
    }

    for (int i = 0; i < static_cast<int>(trees[tree.index].pq_plus_plus.size()); ++i) {
        auto [tree_neighbor, queue_index] = trees[tree.index].pq_plus_plus[i];
        if (trees[tree_neighbor.index].is_alive) {
            EdgeIndex edge = MinPlusPlusExternalEdge(queue_index);
            if (edge) {
                result.emplace_back(tree_neighbor, SlamQuadrupled(edge_heaps[queue_index]->Top()));
            }
        } else {
            trees[tree.index].pq_plus_plus[i] = trees[tree.index].pq_plus_plus.back();
            trees[tree.index].pq_plus_plus.pop_back();
            --i;
        }
    }
    return result;
}

std::vector<std::pair<VzhuhSolver::TreeIndex, int> > VzhuhSolver::PlusMinusExternalSlacks(TreeIndex tree) {
    // TODO avoid code duplication

    std::vector<std::pair<TreeIndex, int> > result;
    result.reserve(trees[tree.index].pq_plus_minus.size());

    for (int i = 0; i < static_cast<int>(trees[tree.index].pq_plus_minus.size()); ++i) {
        auto [tree_neighbor, queue_index] = trees[tree.index].pq_plus_minus[i];
        if (trees[tree_neighbor.index].is_alive) {
            EdgeIndex edge = MinPlusMinusExternalEdge(queue_index);
            if (edge) {
                result.emplace_back(tree_neighbor, SlamQuadrupled(edge_heaps[queue_index]->Top()));
            }
        } else {
            trees[tree.index].pq_plus_minus[i] = trees[tree.index].pq_plus_minus.back();
            trees[tree.index].pq_plus_minus.pop_back();
            --i;
        }
    }
    return result;
}

void VzhuhSolver::DeleteDuplicates(std::vector<NodeIndex> *node_list) {
    ++nodes_label_cnt;
    for (int i = 0; i < static_cast<int>(node_list->size()); ++i) {
        if (nodes[(*node_list)[i].index].label != nodes_label_cnt) {
            nodes[(*node_list)[i].index].label = nodes_label_cnt;
        } else {
            (*node_list)[i] = node_list->back();
            node_list->pop_back();
            --i;
        }
    }
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
