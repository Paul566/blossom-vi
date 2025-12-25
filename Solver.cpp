#include "Solver.h"

#include <iostream>
#include <queue>
#include <unordered_set>

bool Solver::NodeComparator::operator()(const NodeIndex &a, const NodeIndex &b) const {
    int a_var = solver->DualVariableQuadrupled(a);
    int b_var = solver->DualVariableQuadrupled(b);
    if (a_var != b_var) {
        return a_var < b_var;
    }
    return a.index < b.index;
}

bool Solver::EdgeComparator::operator()(const EdgeIndex &a, const EdgeIndex &b) const {
    int a_var = solver->SlackQuadrupled(a);
    int b_var = solver->SlackQuadrupled(b);
    if (a_var != b_var) {
        return a_var < b_var;
    }
    return a.index < b.index;
}

bool Solver::NodeValidator::operator()(const NodeIndex &node) const {
    return true;
}

bool Solver::EdgeValidator::operator()(const EdgeIndex &edge) const {
    return true;
}

Solver::Solver(const std::vector<std::tuple<int, int, int> > &edge_list_,
               const SolverParameters &params_) : primal_objective(INT64_MAX), dual_objective(INT64_MIN),
                                                  num_vertices_elementary(InitNumVertices(edge_list_)),
                                                  num_trees_alive(0),
                                                  first_alive_tree(-1), params(params_) {
    // TODO try permuting nodes and edges
    // TODO reserve more space for nodes

    // initialize nodes
    nodes.is_alive = std::vector<bool>(num_vertices_elementary, true);
    nodes.neighbors = std::vector<std::list<EdgeIndex> >(num_vertices_elementary, std::list<EdgeIndex>());
    nodes.children_neighbors_boundaries = std::vector<std::vector<std::list<EdgeIndex>::iterator> >(
        num_vertices_elementary,
        std::vector<std::list<EdgeIndex>::iterator>());
    nodes.dual_var_quadrupled_amortized = std::vector<int>(num_vertices_elementary, 0);
    nodes.matched_edge = std::vector<EdgeIndex>(num_vertices_elementary, EdgeIndex(-1));
    nodes.blossom_parent = std::vector<NodeIndex>(num_vertices_elementary, NodeIndex(-1));
    nodes.blossom_children = std::vector<std::vector<NodeIndex> >(num_vertices_elementary, std::vector<NodeIndex>());
    nodes.blossom_brother_clockwise = std::vector<EdgeIndex>(num_vertices_elementary, EdgeIndex(-1));
    nodes.blossom_brother_anticlockwise = std::vector<EdgeIndex>(num_vertices_elementary, EdgeIndex(-1));
    nodes.plus = std::vector<bool>(num_vertices_elementary, false);
    nodes.tree_parent = std::vector<EdgeIndex>(num_vertices_elementary, EdgeIndex(-1));
    nodes.tree_children = std::vector<std::vector<EdgeIndex> >(num_vertices_elementary, std::vector<EdgeIndex>());
    nodes.tree = std::vector<TreeIndex>(num_vertices_elementary, TreeIndex(-1));
    nodes.queue_index = std::vector<int>(num_vertices_elementary, -1);
    nodes.handle = std::vector<NodeHeap::Handle *>(num_vertices_elementary);

    // initialize edges and fill the adjacency list
    edges.weight.reserve(edge_list_.size());
    edges.matched = std::vector<bool>(edge_list_.size(), false);
    edges.slack_quadrupled_amortized.reserve(edge_list_.size());
    edges.head.reserve(edge_list_.size());
    edges.tail.reserve(edge_list_.size());
    edges.queue_index = std::vector<int>(edge_list_.size(), -1);
    edges.handle = std::vector<EdgeHeap::Handle *>(edge_list_.size());
    for (int i = 0; i < static_cast<int>(edge_list_.size()); ++i) {
        edges.weight.push_back(std::get<2>(edge_list_[i]));
        edges.slack_quadrupled_amortized.push_back(4 * std::get<2>(edge_list_[i]));
        edges.head.push_back(NodeIndex(std::get<0>(edge_list_[i])));
        edges.tail.push_back(NodeIndex(std::get<1>(edge_list_[i])));

        nodes.neighbors[std::get<0>(edge_list_[i])].push_back(EdgeIndex(i));
        nodes.neighbors[std::get<1>(edge_list_[i])].push_back(EdgeIndex(i));
    }
}

void Solver::FindMinPerfectMatching() {
    GreedyInit();
    InitializeTrees();

    if (params.verbose) {
        PrintGraph();
        PrintTrees();
    }

    int num_rounds = 0;
    while (num_trees_alive > 0) {
        ++num_rounds;
        MakePrimalUpdates();
        // ValidateQueues();
        MakeDualUpdates();
        // ValidatePositiveSlacks();

        if (params.print_statistics) {
            std::cout << "round " << num_rounds << std::endl;
            PrintStatistics();
        }
    }

    // compute the objectives and recover the matching
    int64_t quadrupled_dual = DualObjectiveQuadrupled();
    if (quadrupled_dual % 4 != 0) {
        throw std::runtime_error("Dual objective not integer");
    }
    dual_objective = quadrupled_dual / 4;
    if (params.compute_dual_certificate) {
        ComputeDualCertificate();
    }
    if (params.verbose) {
        std::cout << "Dual objective:\t\t" << dual_objective << std::endl;
    }

    if (params.print_statistics) {
        std::cout << "number of rounds: " << num_rounds << std::endl;
        PrintStatistics();
    }

    DestroyBlossoms();
    ComputeMatching();

    primal_objective = PrimalObjective();

    if (params.verbose) {
        std::cout << "Primal objective:\t" << primal_objective << std::endl;
    }
}

void Solver::PrintGraph() const {
    std::cout << "Adjacency list (to, weight, slack, matched):" << std::endl;
    for (int i = 0; i < static_cast<int>(nodes.plus.size()); ++i) {
        if (nodes.is_alive[i]) {
            PrintNode(NodeIndex(i));
        }
    }
}

void Solver::PrintNode(NodeIndex node) const {
    std::cout << node.index << ": y_v = " << DualVariableQuadrupled(node) / 4.;
    if (!IsTopBlossom(node)) {
        std::cout << " blossom_parent: " << nodes.blossom_parent[node.index].index << std::endl;
        return;
    }
    for (const EdgeIndex edge : nodes.neighbors[node.index]) {
        std::cout << " (" << OtherEnd(edge, node).index << ", " << edges.weight[edge.index] << ", " <<
            SlackQuadrupled(edge) / 4. << ", " << edges.matched[edge.index] << ")";
    }
    std::cout << std::endl;
}

void Solver::PrintTrees() const {
    std::cout << "Trees:" << std::endl;

    for (int tree_index = 0; tree_index < static_cast<int>(trees.is_alive.size()); ++tree_index) {
        if (!trees.is_alive[tree_index]) {
            continue;
        }

        std::cout << "Root: " << trees.root[tree_index].index;
        std::cout << ", dual variable: " << trees.dual_var_quadrupled[tree_index] / 4. << std::endl;

        struct Frame {
            const NodeIndex node;
            std::string prefix; // prefix to print before connector / node value
            bool hasConnector; // whether to print "├── " or "└── " before the node value
            bool isLast; // if hasConnector==true, whether this node is the last sibling
        };

        std::vector<Frame> stack;
        stack.push_back({TopBlossom(trees.root[tree_index]), "", false, false});

        while (!stack.empty()) {
            Frame f = stack.back();
            stack.pop_back();

            if (f.hasConnector) {
                std::cout << f.prefix << (f.isLast ? "└── " : "├── ");
            } else {
                std::cout << f.prefix; // usually empty for root
            }
            std::cout << f.node.index << '\n';

            const std::vector<EdgeIndex> &ch = nodes.tree_children[f.node.index];
            for (int i = static_cast<int>(ch.size()) - 1; i >= 0; --i) {
                bool childIsLast = (i == static_cast<int>(ch.size()) - 1);

                std::string childPrefix = f.prefix;
                if (f.hasConnector) {
                    childPrefix += (f.isLast ? "    " : "│   ");
                }

                stack.push_back({OtherEnd(ch[i], f.node), childPrefix, true, childIsLast});
            }
        }
    }
}

void Solver::PrintStatistics() {
    int alive_trees = 0;
    for (int i = 0; i < static_cast<int>(trees.is_alive.size()); ++i) {
        if (trees.is_alive[i]) {
            ++alive_trees;
        }
    }
    std::cout << "number of trees: " << alive_trees << std::endl;

    int max_tree_deg = 0;
    for (int i = 0; i < static_cast<int>(trees.is_alive.size()); ++i) {
        if (trees.is_alive[i]) {
            if (static_cast<int>(trees.pq_minus_plus[i].size()) > max_tree_deg) {
                max_tree_deg = static_cast<int>(trees.pq_minus_plus[i].size());
            }
            if (static_cast<int>(trees.pq_plus_minus[i].size()) > max_tree_deg) {
                max_tree_deg = static_cast<int>(trees.pq_plus_minus[i].size());
            }
            if (static_cast<int>(trees.pq_plus_plus[i].size()) > max_tree_deg) {
                max_tree_deg = static_cast<int>(trees.pq_plus_plus[i].size());
            }
        }
    }
    std::cout << "max degree of a tree: " << max_tree_deg << std::endl;

    int num_blossoms = 0;
    for (int i = num_vertices_elementary; i < static_cast<int>(nodes.is_alive.size()); ++i) {
        if (nodes.is_alive[i]) {
            ++num_blossoms;
        }
    }
    std::cout << "number of blossoms: " << num_blossoms << std::endl;

    int max_deg = 0;
    for (int i = num_vertices_elementary; i < static_cast<int>(nodes.is_alive.size()); ++i) {
        if (nodes.is_alive[i]) {
            if (max_deg < static_cast<int>(nodes.neighbors[i].size())) {
                max_deg = static_cast<int>(nodes.neighbors[i].size());
            }
        }
    }
    std::cout << "max degree of a blossom: " << max_deg << std::endl;

    int max_depth = 0;
    for (int i = 0; i < num_vertices_elementary; ++i) {
        int depth = NodeDepth(NodeIndex(i));
        if (max_depth < depth) {
            max_depth = depth;
        }
    }
    std::cout << "max depth of a blossom: " << max_depth << std::endl;
}

const std::vector<std::pair<int, int> > &Solver::Matching() const {
    return matching;
}

const std::vector<std::tuple<int, int, int> > &Solver::DualCertificate() const {
    return dual_certificate;
}

void Solver::GreedyInit() {
    // first, make all the slacks non-negative
    for (int node_index = 0; node_index < num_vertices_elementary; ++node_index) {
        if (nodes.neighbors[node_index].empty()) {
            std::cout << "vertex " << node_index << ": no neighbors" << std::endl;
            throw std::runtime_error("Found an isolated vertex => no perfect matching exists");
        }

        int min_weight = edges.weight[nodes.neighbors[node_index].front().index];
        for (const EdgeIndex edge : nodes.neighbors[node_index]) {
            if (edges.weight[edge.index] < min_weight) {
                min_weight = edges.weight[edge.index];
            }
        }

        nodes.dual_var_quadrupled_amortized[node_index] += min_weight * 2;
    }

    for (int node_index = 0; node_index < num_vertices_elementary; ++node_index) {
        if (nodes.matched_edge[node_index]) {
            continue;
        }

        EdgeIndex smallest_slack_edge = nodes.neighbors[node_index].front();
        for (const EdgeIndex edge : nodes.neighbors[node_index]) {
            if (SlackQuadrupled(edge) < SlackQuadrupled(smallest_slack_edge)) {
                smallest_slack_edge = edge;
            }
        }

        nodes.dual_var_quadrupled_amortized[node_index] += SlackQuadrupled(smallest_slack_edge);
        if (!nodes.matched_edge[OtherEnd(smallest_slack_edge, NodeIndex(node_index)).index]) {
            // if the other vertex is also unmatched, match the edge
            MakeEdgeMatched(smallest_slack_edge);
        }
    }
}

void Solver::InitializeTrees() {
    // is called after GreedyInit

    for (int root_index = 0; root_index < num_vertices_elementary; ++root_index) {
        if (nodes.matched_edge[root_index]) {
            continue;
        }
        nodes.tree[root_index] = TreeIndex(static_cast<int>(trees.root.size()));
        nodes.plus[root_index] = true;
        trees.root.push_back(NodeIndex(root_index));
    }

    num_trees_alive = static_cast<int>(trees.root.size());
    first_alive_tree = TreeIndex(0);
    trees.is_alive = std::vector<bool>(num_trees_alive, true);
    trees.dual_var_quadrupled = std::vector<int>(num_trees_alive, 0);

    trees.next_alive_tree.reserve(num_trees_alive);
    trees.alive_index.reserve(num_trees_alive);
    for (int i = 0; i < num_trees_alive - 1; ++i) {
        trees.next_alive_tree.push_back(TreeIndex(i + 1));
        trees.alive_index.push_back(i);
    }
    trees.next_alive_tree.push_back(TreeIndex(0));
    trees.alive_index.push_back(num_trees_alive - 1);

    trees.minus_blossoms.reserve(num_trees_alive);
    trees.plus_empty_edges.reserve(num_trees_alive);
    trees.plus_plus_internal_edges.reserve(num_trees_alive);
    trees.pq_plus_plus = std::vector<std::vector<std::pair<TreeIndex, int> > >(
        num_trees_alive,
        std::vector<std::pair<TreeIndex, int> >());
    trees.pq_plus_minus = std::vector<std::vector<std::pair<TreeIndex, int> > >(
        num_trees_alive,
        std::vector<std::pair<TreeIndex, int> >());
    trees.pq_minus_plus = std::vector<std::vector<std::pair<TreeIndex, int> > >(
        num_trees_alive,
        std::vector<std::pair<TreeIndex, int> >());
    queues.node_heaps.reserve(num_trees_alive);
    queues.edge_heaps.reserve(2 * num_trees_alive);
    for (int i = 0; i < num_trees_alive; ++i) {
        trees.minus_blossoms.push_back(i);
        trees.plus_empty_edges.push_back(2 * i);
        trees.plus_plus_internal_edges.push_back(2 * i + 1);

        queues.node_heaps.emplace_back(std::make_unique<NodeHeap>(params.heap_arity, NodeComparator{this}, NodeValidator{this}));
        queues.edge_heaps.emplace_back(std::make_unique<EdgeHeap>(params.heap_arity, EdgeComparator{this}, EdgeValidator{this}));
        queues.edge_heaps.emplace_back(std::make_unique<EdgeHeap>(params.heap_arity, EdgeComparator{this}, EdgeValidator{this}));
    }

    InitializeQueues();
}

void Solver::InitializeQueues() {
    for (NodeIndex root : trees.root) {
        TreeIndex tree = nodes.tree[root.index];

        for (EdgeIndex edge_to_neighbor : nodes.neighbors[root.index]) {
            NodeIndex neighbor = OtherEnd(edge_to_neighbor, root);

            if (!nodes.tree[neighbor.index]) {
                AddEdgeToQueue(edge_to_neighbor, trees.plus_empty_edges[tree.index]);
            } else {
                // neighbor is another root
                AddPQPlusPlus(tree, nodes.tree[neighbor.index], edge_to_neighbor);
            }
        }
    }
}

void Solver::MakeDualUpdates() {
    // makes the dual updates for the whole collection of the trees

    if (num_trees_alive == 0) {
        return;
    }

    if (params.verbose) {
        std::cout << "DUAL UPDATE" << std::endl;
    }

    std::vector<int> deltas = VariableDeltas();
    TreeIndex tree = FirstAliveTree();
    int i = 0;
    do {
        trees.dual_var_quadrupled[tree.index] += deltas[i];
        tree = NextAliveTree(tree);
        ++i;
    } while (tree != FirstAliveTree());

    if (params.verbose) {
        std::cout << "after dual update:" << std::endl;
        PrintGraph();
    }
}

std::vector<int> Solver::VariableDeltas() {
    DualConstraints dual_constraints = GetDualConstraints();
    std::vector<std::vector<int> > connected_components = ConnectedComponentsTreeTree(dual_constraints);
    std::vector<int> deltas(num_trees_alive, 0);

    std::vector<int> component_index(num_trees_alive, 0);
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

    return deltas;
}
std::vector<std::vector<int> > Solver::ConnectedComponentsTreeTree(const DualConstraints &dual_constraints) const {
    // returns a vector of index_of_connected_component

    std::vector<bool> visited(num_trees_alive, false);
    std::vector<std::vector<int> > components;

    // Build undirected adjacency list for weak connectivity
    std::vector<std::vector<int> > adj_list(num_trees_alive, std::vector<int>());
    for (int u = 0; u < num_trees_alive; ++u) {
        for (auto &[v, slack] : dual_constraints.plus_minus_constraints[u]) {
            if (slack == 0) {
                adj_list[u].emplace_back(v);
                adj_list[v].emplace_back(u);
            }
        }
    }

    for (int start = 0; start < num_trees_alive; ++start) {
        if (!visited[start]) {
            std::vector<int> comp;
            std::queue<int> q;

            visited[start] = true;
            q.push(start);

            while (!q.empty()) {
                int u = q.front();
                q.pop();
                comp.push_back(u);

                for (int v : adj_list[u]) {
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

Solver::DualConstraints Solver::GetDualConstraints() {
    UpdateAliveIndices();

    // get upper_bound
    std::vector<int> upper_bound(num_trees_alive, INT32_MAX);
    int i = 0;
    TreeIndex tree = FirstAliveTree();
    do {
        int plus_empty = PlusEmptySlack(tree);
        int plus_plus_internal = PlusPlusInternalSlack(tree);
        int blossom_var = MinMinusBlossomVariable(tree);

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

        ++i;
        tree = NextAliveTree(tree);
    } while (tree != FirstAliveTree());

    // get plus_plus_constraints and plus_minus_constraints
    std::vector plus_plus_constraints(num_trees_alive, std::vector<std::pair<int, int> >());
    std::vector plus_minus_constraints(num_trees_alive, std::vector<std::pair<int, int> >());
    i = 0;
    tree = FirstAliveTree();
    do {
        std::vector<std::pair<TreeIndex, int> > plus_plus_slacks = PlusPlusExternalSlacks(tree);
        std::vector<std::pair<TreeIndex, int> > plus_minus_slacks = PlusMinusExternalSlacks(tree);
        for (auto [other_tree, slack] : plus_plus_slacks) {
            plus_plus_constraints[i].emplace_back(trees.alive_index[other_tree.index], slack);
        }
        for (auto [other_tree, slack] : plus_minus_slacks) {
            plus_minus_constraints[i].emplace_back(trees.alive_index[other_tree.index], slack);
        }
        ++i;
        tree = NextAliveTree(tree);
    } while (tree != FirstAliveTree());

    return std::move(DualConstraints{
        std::move(upper_bound), std::move(plus_plus_constraints), std::move(plus_minus_constraints)
    });
}

void Solver::MakePrimalUpdates() {
    // makes primal updates for the whole collection of the trees

    // TODO check if the sizes of the trees become imbalanced

    if (params.verbose) {
        std::cout << "PRIMAL UPDATES" << std::endl;
        PrintGraph();
        PrintTrees();
    }

    UpdateAliveIndices();
    std::vector<TreeIndex> next_active_tree;
    std::vector<bool> tree_is_active(num_trees_alive, true);
    next_active_tree.reserve(num_trees_alive);
    TreeIndex tree = FirstAliveTree();
    do {
        next_active_tree.push_back(NextAliveTree(tree));
        tree = NextAliveTree(tree);
    } while (tree != FirstAliveTree());

    tree = FirstAliveTree();
    while (tree) {
        bool success = true;
        TreeIndex augmented_tree = MakePrimalUpdate(tree, &success);

        // if (params.verbose) {
        //     std::cout << "VALIDATION" << std::endl;
        //     ValidateQueues();
        // }

        if (augmented_tree) {
            trees.is_alive[tree.index] = false;
            trees.is_alive[augmented_tree.index] = false;
            num_trees_alive -= 2;

            tree_is_active[trees.alive_index[tree.index]] = false;
            tree_is_active[trees.alive_index[augmented_tree.index]] = false;
        }

        if (!success) {
            tree_is_active[trees.alive_index[tree.index]] = false;
        }

        TreeIndex next_tree = next_active_tree[trees.alive_index[tree.index]];
        while (!tree_is_active[trees.alive_index[next_tree.index]]) {
            next_tree = next_active_tree[trees.alive_index[next_tree.index]];
            if (next_tree == tree) {
                next_tree = TreeIndex(-1);
                break;
            }
        }
        next_active_tree[trees.alive_index[tree.index]] = next_tree;
        tree = next_tree;
    }

    if (params.verbose) {
        std::cout << "after updates:" << std::endl;
        PrintGraph();
        PrintTrees();
    }
}

Solver::TreeIndex Solver::MakePrimalUpdate(TreeIndex tree, bool *success) {
    if (EdgeIndex augmentable_edge = AugmentableEdge(tree)) {
        TreeIndex other = Augment(tree, augmentable_edge);
        return other;
    }

    if (EdgeIndex growable_edge = GrowableEdge(tree)) {
        Grow(tree, growable_edge);
        return TreeIndex(-1);
    }

    if (EdgeIndex shrinkable_edge = ShrinkableEdge(tree)) {
        Shrink(tree, shrinkable_edge);
        return TreeIndex(-1);
    }

    if (NodeIndex expandable_blossom = ExpandableBlossom(tree)) {
        Expand(tree, expandable_blossom);
        return TreeIndex(-1);
    }

    *success = false;
    return TreeIndex(-1);
}

void Solver::Grow(TreeIndex tree, EdgeIndex edge) {
    NodeIndex parent = edges.head[edge.index];
    NodeIndex child = edges.tail[edge.index];
    if (nodes.tree[child.index] == tree) {
        std::swap(parent, child);
    }

    if (nodes.tree[parent.index] != tree) {
        throw std::runtime_error("In Tree::Grow: parent vertex is not in this tree");
    }
    if (!nodes.plus[parent.index]) {
        throw std::runtime_error("In Tree::Grow: parent is not a plus");
    }
    if (nodes.tree[child.index]) {
        throw std::runtime_error("In Tree::Grow: child vertex is not free");
    }

    NodeIndex grandchild = OtherEnd(nodes.matched_edge[child.index], child);
    if (params.verbose) {
        std::cout << "GROW " << parent.index << " " << child.index << " " << grandchild.index << std::endl;
    }

    nodes.tree_children[parent.index].push_back(edge);

    nodes.tree_parent[child.index] = edge;
    nodes.tree[child.index] = tree;
    nodes.plus[child.index] = false;
    nodes.dual_var_quadrupled_amortized[child.index] += trees.dual_var_quadrupled[tree.index];
    nodes.tree_children[child.index] = {nodes.matched_edge[child.index]};

    nodes.tree_parent[grandchild.index] = nodes.matched_edge[child.index];
    nodes.tree[grandchild.index] = tree;
    nodes.plus[grandchild.index] = true;
    nodes.dual_var_quadrupled_amortized[grandchild.index] -= trees.dual_var_quadrupled[tree.index];

    UpdateQueuesAfterGrow(child, grandchild);
}

void Solver::Shrink(TreeIndex tree, EdgeIndex edge_plus_plus) {
    NodeIndex first = edges.head[edge_plus_plus.index];
    NodeIndex second = edges.tail[edge_plus_plus.index];

    if (params.verbose) {
        std::cout << "SHRINK " << first.index << " " << second.index << std::endl;
    }

    NodeIndex lca = LCA(first, second);

    // TODO maybe make LCA also compute the blossom edges
    std::vector<EdgeIndex> blossom_edges;
    while (first != lca) {
        blossom_edges.push_back(nodes.tree_parent[first.index]);
        first = OtherEnd(nodes.tree_parent[first.index], first);
    }
    std::reverse(blossom_edges.begin(), blossom_edges.end());
    blossom_edges.push_back(edge_plus_plus);
    while (second != lca) {
        blossom_edges.push_back(nodes.tree_parent[second.index]);
        second = OtherEnd(nodes.tree_parent[second.index], second);
    }

    NodeIndex new_node(static_cast<int>(nodes.plus.size()));
    MakeBlossom(blossom_edges, lca);
}

void Solver::MakeBlossom(std::vector<EdgeIndex> blossom_edges, NodeIndex lca) {
    // blossom_edges need to be consecutive and start and end at the receptacle

    if (blossom_edges.size() % 2 == 0) {
        throw std::runtime_error("In Node: blossom has to have an odd number of vertices");
    }

    NodeIndex new_blossom(static_cast<int>(nodes.plus.size()));

    nodes.is_alive.emplace_back(true);
    nodes.blossom_parent.emplace_back(-1);
    nodes.blossom_brother_clockwise.emplace_back(-1);
    nodes.blossom_brother_anticlockwise.emplace_back(-1);
    nodes.plus.emplace_back(true); // after Shrink, we must be a plus
    nodes.queue_index.emplace_back(-1);
    nodes.handle.emplace_back();

    // update blossom_children
    nodes.blossom_children.emplace_back();
    nodes.blossom_children.back().reserve(blossom_edges.size());
    NodeIndex cur_vertex = SharedNode(blossom_edges.front(), blossom_edges.back());
    TreeIndex this_tree = nodes.tree[cur_vertex.index];
    for (EdgeIndex edge : blossom_edges) {
        cur_vertex = OtherEnd(edge, cur_vertex);
        nodes.blossom_children.back().push_back(cur_vertex);
    }

    // update the queues
    std::vector<NodeIndex> minus_children;
    for (NodeIndex child : nodes.blossom_children.back()) {
        if (!nodes.plus[child.index]) {
            minus_children.push_back(child);
        }
    }
    UpdateQueuesBeforeShrink(minus_children);

    // initialize blossom brothers of the children
    for (EdgeIndex edge : blossom_edges) {
        nodes.blossom_brother_clockwise[cur_vertex.index] = edge;
        cur_vertex = OtherEnd(edge, cur_vertex);
        nodes.blossom_brother_anticlockwise[cur_vertex.index] = edge;
    }

    // update blossom_parent and dual_var_quadrupled_amortized of the children
    for (NodeIndex vertex : nodes.blossom_children.back()) {
        if (nodes.blossom_parent[vertex.index]) {
            throw std::runtime_error("In Node: some child node already has a parent");
        }
        nodes.blossom_parent[vertex.index] = new_blossom;
        if (nodes.plus[vertex.index]) {
            nodes.dual_var_quadrupled_amortized[vertex.index] += trees.dual_var_quadrupled[this_tree.index];
        } else {
            nodes.dual_var_quadrupled_amortized[vertex.index] -= trees.dual_var_quadrupled[this_tree.index];
        }
    }

    std::unordered_set<int> blossom_set;
    // TODO try vector + linear search, check if blossoms usually have few children
    blossom_set.reserve(nodes.blossom_children.back().size());
    for (NodeIndex vertex : nodes.blossom_children.back()) {
        blossom_set.insert(vertex.index);
    }

    // find tree_children
    nodes.tree_children.emplace_back();
    for (NodeIndex child : nodes.blossom_children.back()) {
        for (EdgeIndex child_edge : nodes.tree_children[child.index]) {
            if (!blossom_set.contains(OtherEnd(child_edge, child).index)) {
                nodes.tree_children.back().push_back(child_edge);
            }
        }
    }

    // update edges TODO amortize this
    for (NodeIndex child : nodes.blossom_children.back()) {
        for (EdgeIndex edge : nodes.neighbors[child.index]) {
            if (!blossom_set.contains(OtherEnd(edge, child).index)) {
                UpdateEdgeAfterShrink(edge, child);
            }
        }
    }

    // make neighbors, children_neighbors_boundaries
    nodes.neighbors.emplace_back();
    nodes.children_neighbors_boundaries.emplace_back();
    for (NodeIndex child : nodes.blossom_children.back()) {
        auto segment_start = nodes.neighbors[child.index].begin();
        nodes.neighbors.back().splice(nodes.neighbors.back().end(), nodes.neighbors[child.index]);
        nodes.children_neighbors_boundaries.back().push_back(segment_start);
    }

    nodes.matched_edge.emplace_back(nodes.matched_edge[lca.index]);
    nodes.tree_parent.emplace_back(nodes.matched_edge.back());
    nodes.tree.emplace_back(this_tree);
    nodes.dual_var_quadrupled_amortized.emplace_back(-trees.dual_var_quadrupled[this_tree.index]);
    // the true dual variable must be zero

    RemoveLoopsFromQueues(new_blossom);
}

void Solver::Expand(TreeIndex tree, NodeIndex blossom) {
    if (params.verbose) {
        std::cout << "EXPAND " << blossom.index << std::endl;
    }

    if (IsElementary(blossom)) {
        throw std::runtime_error("In Tree::Expand: supervertex must be not elementary");
    }
    if (nodes.blossom_parent[blossom.index]) {
        throw std::runtime_error("In Tree::Expand: supervertex must be a top blossom");
    }
    if (nodes.plus[blossom.index]) {
        throw std::runtime_error("In Tree::Expand: supervertex is not a minus");
    }
    if (DualVariableQuadrupled(blossom) != 0) {
        throw std::runtime_error("In Tree::Expand: the supervertex has to have zero dual variable");
    }

    std::vector<NodeIndex> children = nodes.blossom_children[blossom.index];
    Dissolve(blossom);
    UpdateQueuesAfterExpand(tree, blossom, children);
}

Solver::TreeIndex Solver::Augment(TreeIndex tree, EdgeIndex edge) {
    // TODO amortize tree dissolves

    NodeIndex parent = edges.head[edge.index];
    NodeIndex child = edges.tail[edge.index];
    if (nodes.tree[parent.index] != tree) {
        std::swap(parent, child);
    }

    if (params.verbose) {
        std::cout << "AUGMENT " << parent.index << " " << child.index << std::endl;
    }

    if (nodes.tree[parent.index] != tree) {
        throw std::runtime_error("In Augment: parent vertex is not in this tree");
    }
    if (!nodes.tree[child.index]) {
        throw std::runtime_error("In Augment: child is not in a tree");
    }

    MakeEdgeMatched(edge);

    AugmentFromNode(tree, parent);
    TreeIndex other_tree = nodes.tree[child.index];
    AugmentFromNode(other_tree, child);

    return other_tree;
}

void Solver::AugmentFromNode(TreeIndex tree, NodeIndex node) {
    bool match = false;
    const std::vector<EdgeIndex> path = PathToRoot(tree, node);
    for (EdgeIndex edge_from_path : path) {
        if (match) {
            MakeEdgeMatched(edge_from_path);
        } else {
            MakeEdgeUnmatched(edge_from_path);
        }
        match = !match;
    }

    DissolveTree(tree);
}

void Solver::DissolveTree(TreeIndex tree) {
    // clears the tree structure from the nodes
    // updates priority queues of the adjacent trees

    trees.is_alive[tree.index] = false;

    std::vector<NodeIndex> all_vertices;
    std::queue<NodeIndex> queue;
    queue.push(TopBlossom(trees.root[tree.index]));
    while (!queue.empty()) {
        NodeIndex cur_vertex = queue.front();
        queue.pop();
        all_vertices.push_back(cur_vertex);
        for (EdgeIndex child_edge : nodes.tree_children[cur_vertex.index]) {
            queue.push(OtherEnd(child_edge, cur_vertex));
        }
    }

    // update plus_empty_edges of other trees
    for (NodeIndex vertex : all_vertices) {
        for (EdgeIndex neighbor_edge : nodes.neighbors[vertex.index]) {
            NodeIndex neighbor = OtherEnd(neighbor_edge, vertex);
            if ((nodes.tree[neighbor.index] != tree) && nodes.plus[neighbor.index]) {
                AddEdgeToQueue(neighbor_edge, trees.plus_empty_edges[nodes.tree[neighbor.index].index]);
            }
        }
    }

    for (NodeIndex vertex : all_vertices) {
        ClearNodeDuringTreeDissolve(tree, vertex);
    }
}
void Solver::ClearNodeDuringTreeDissolve(TreeIndex tree, NodeIndex node) {
    if (nodes.plus[node.index]) {
        nodes.dual_var_quadrupled_amortized[node.index] += trees.dual_var_quadrupled[tree.index];
    } else {
        nodes.dual_var_quadrupled_amortized[node.index] -= trees.dual_var_quadrupled[tree.index];
    }
    nodes.tree[node.index] = TreeIndex(-1);
    nodes.tree_children[node.index].clear();
    nodes.tree_parent[node.index] = EdgeIndex(-1);
    nodes.plus[node.index] = false;
}

void Solver::Dissolve(NodeIndex node) {
    // TODO amortize

    if (nodes.blossom_children[node.index].empty()) {
        throw std::runtime_error("In Dissolve: can't dissolve an elementary vertex");
    }
    if (!nodes.matched_edge[node.index]) {
        throw std::runtime_error("In Dissolve: there is no matched_edge -- this blossom can't contain the root");
    }
    if (nodes.blossom_parent[node.index]) {
        throw std::runtime_error("In Dissolve: this must be a top blossom");
    }

    nodes.is_alive[node.index] = false;
    RotateReceptacle(node, DeeperNode(node, nodes.matched_edge[node.index]));

    if (nodes.tree[node.index]) {
        UpdateNodeInternalTreeStructure(node);
        for (NodeIndex child : nodes.blossom_children[node.index]) {
            if (nodes.plus[child.index]) {
                nodes.dual_var_quadrupled_amortized[child.index] -= trees.dual_var_quadrupled[nodes.tree[node.index].
                    index];
            }
            if (nodes.tree[child.index] && !nodes.plus[child.index]) {
                nodes.dual_var_quadrupled_amortized[child.index] += trees.dual_var_quadrupled[nodes.tree[node.index].
                    index];
            }
        }
    } else {
        ClearNodeInternalTreeStructure(node);
    }

    for (NodeIndex child : nodes.blossom_children[node.index]) {
        if (!child) {
            throw std::runtime_error("Invalid child");
        }
        nodes.blossom_brother_clockwise[child.index] = EdgeIndex(-1);
        nodes.blossom_brother_anticlockwise[child.index] = EdgeIndex(-1);
        nodes.blossom_parent[child.index] = NodeIndex(-1);
    }

    // restore neighbors of the children
    for (int i = 0; i < static_cast<int>(nodes.blossom_children[node.index].size()); ++i) {
        auto it_left = nodes.children_neighbors_boundaries[node.index][i];
        auto it_right = nodes.neighbors[node.index].end();
        if (i + 1 < static_cast<int>(nodes.blossom_children[node.index].size())) {
            it_right = nodes.children_neighbors_boundaries[node.index][i + 1];
        }

        nodes.neighbors[nodes.blossom_children[node.index][i].index].splice(
            nodes.neighbors[nodes.blossom_children[node.index][i].index].end(),
            nodes.neighbors[node.index],
            it_left,
            it_right);
    }

    // update the edges
    for (NodeIndex child : nodes.blossom_children[node.index]) {
        for (EdgeIndex edge : nodes.neighbors[child.index]) {
            if (node == edges.tail[edge.index]) {
                edges.tail[edge.index] = child;
                edges.slack_quadrupled_amortized[edge.index] += DualVariableQuadrupled(edges.tail[edge.index]);
            } else if (node == edges.head[edge.index]) {
                edges.head[edge.index] = child;
                edges.slack_quadrupled_amortized[edge.index] += DualVariableQuadrupled(edges.head[edge.index]);
            }
        }
    }
}

void Solver::UpdateEdgeAfterShrink(EdgeIndex edge, NodeIndex node) {
    edges.slack_quadrupled_amortized[edge.index] -= DualVariableQuadrupled(node);
    if (node == edges.tail[edge.index]) {
        edges.tail[edge.index] = nodes.blossom_parent[node.index];
        return;
    }
    edges.head[edge.index] = nodes.blossom_parent[node.index];
}

void Solver::UpdateNodeInternalTreeStructure(NodeIndex node) {
    // updates tree_children, tree_parents, tree_root for the vertices inside this blossom before the Expand
    // updates the plus and minus markers
    if (!nodes.matched_edge[node.index]) {
        throw std::runtime_error("In UpdateInternalTreeStructure: there is no matched_edge");
    }
    if (!nodes.tree_parent[node.index]) {
        throw std::runtime_error("In UpdateInternalTreeStructure: there is no tree_parent");
    }
    NodeIndex elder_child = DeeperNode(node, nodes.tree_parent[node.index]);
    NodeIndex receptacle = DeeperNode(node, nodes.matched_edge[node.index]);

    std::function next_edge = [this](const NodeIndex current) -> EdgeIndex {
        return nodes.blossom_brother_clockwise[current.index];
    };
    if (edges.matched[nodes.blossom_brother_anticlockwise[elder_child.index].index]) {
        next_edge = [this](const NodeIndex current) -> EdgeIndex {
            return nodes.blossom_brother_anticlockwise[current.index];
        };
    }

    // the part that stays in the tree
    NodeIndex cur_vertex = elder_child;
    EdgeIndex prev_edge = nodes.tree_parent[node.index];
    bool is_plus = false;
    while (cur_vertex != receptacle) {
        nodes.tree_parent[cur_vertex.index] = prev_edge;
        nodes.tree[cur_vertex.index] = nodes.tree[node.index];
        nodes.tree_children[cur_vertex.index] = {next_edge(cur_vertex)};
        nodes.plus[cur_vertex.index] = is_plus;

        prev_edge = next_edge(cur_vertex);
        cur_vertex = OtherEnd(prev_edge, cur_vertex);
        is_plus = !is_plus;
    }
    nodes.tree_parent[cur_vertex.index] = prev_edge;
    nodes.tree[cur_vertex.index] = nodes.tree[node.index];
    nodes.tree_children[cur_vertex.index] = {nodes.matched_edge[node.index]};
    nodes.plus[cur_vertex.index] = false;

    // the part that goes to waste
    cur_vertex = OtherEnd(next_edge(cur_vertex), cur_vertex);
    while (cur_vertex != elder_child) {
        nodes.tree_parent[cur_vertex.index] = EdgeIndex(-1);
        nodes.tree[cur_vertex.index] = TreeIndex(-1);
        nodes.tree_children[cur_vertex.index].clear();
        nodes.plus[cur_vertex.index] = false;

        cur_vertex = OtherEnd(next_edge(cur_vertex), cur_vertex);
    }
}

void Solver::ClearNodeInternalTreeStructure(NodeIndex node) {
    for (NodeIndex child : nodes.blossom_children[node.index]) {
        nodes.tree_parent[child.index] = EdgeIndex(-1);
        nodes.tree[child.index] = TreeIndex(-1);
        nodes.tree_children[child.index].clear();
        nodes.plus[child.index] = false;
    }
}

Solver::EdgeIndex Solver::GrowableEdge(TreeIndex tree) {
    EdgeIndex edge = MinPlusEmptyEdge(trees.plus_empty_edges[tree.index]);
    if (edge) {
        if (SlackQuadrupled(edge) == 0) {
            return edge;
        }
    }
    return EdgeIndex(-1);
}

Solver::EdgeIndex Solver::AugmentableEdge(TreeIndex tree) {
    for (auto [other_tree, queue_index] : trees.pq_plus_plus[tree.index]) {
        EdgeIndex edge = MinPlusPlusExternalEdge(queue_index);
        if (edge) {
            if (SlackQuadrupled(edge) == 0) {
                return edge;
            }
        }
    }
    return EdgeIndex(-1);
}

Solver::EdgeIndex Solver::ShrinkableEdge(TreeIndex tree) {
    EdgeIndex edge = MinPlusPlusInternalEdge(trees.plus_plus_internal_edges[tree.index]);
    if (edge) {
        if (SlackQuadrupled(edge) == 0) {
            return edge;
        }
    }
    return EdgeIndex(-1);
}

Solver::NodeIndex Solver::ExpandableBlossom(TreeIndex tree) {
    NodeIndex node = MinMinusBlossom(trees.minus_blossoms[tree.index]);
    if (node) {
        if (DualVariableQuadrupled(node) == 0) {
            return node;
        }
    }
    return NodeIndex(-1);
}

Solver::EdgeIndex Solver::MinPlusEmptyEdge(int queue_index) {
    while (!queues.edge_heaps[queue_index]->Empty()) {
        EdgeIndex edge = queues.edge_heaps[queue_index]->Top();

        NodeIndex head = edges.head[edge.index]; // plus
        NodeIndex tail = edges.tail[edge.index]; // empty
        if (!nodes.tree[head.index]) {
            std::swap(head, tail);
        }
        if (nodes.plus[head.index] && !nodes.tree[tail.index]) {
            return edge;
        }
        queues.edge_heaps[queue_index]->Pop();
        edges.queue_index[edge.index] = -1;
    }
    return EdgeIndex(-1);
}

Solver::EdgeIndex Solver::MinPlusPlusInternalEdge(int queue_index) {
    while (!queues.edge_heaps[queue_index]->Empty()) {
        EdgeIndex edge = queues.edge_heaps[queue_index]->Top();

        NodeIndex head = edges.head[edge.index];
        NodeIndex tail = edges.tail[edge.index];
        if (nodes.plus[head.index] && nodes.plus[tail.index] && (nodes.tree[head.index] == nodes.tree[tail.index]) && (head != tail) && IsTopBlossom(head) && IsTopBlossom(tail)) {
            return edge;
        }
        queues.edge_heaps[queue_index]->Pop();
        edges.queue_index[edge.index] = -1;
    }
    return EdgeIndex(-1);
}

Solver::EdgeIndex Solver::MinPlusPlusExternalEdge(int queue_index) {
    while (!queues.edge_heaps[queue_index]->Empty()) {
        EdgeIndex edge = queues.edge_heaps[queue_index]->Top();

        NodeIndex head = edges.head[edge.index];
        NodeIndex tail = edges.tail[edge.index];
        if (nodes.plus[head.index] && nodes.plus[tail.index] && (nodes.tree[head.index] != nodes.tree[tail.index])) {
            return edge;
        }
        queues.edge_heaps[queue_index]->Pop();
        edges.queue_index[edge.index] = -1;
    }
    return EdgeIndex(-1);
}

Solver::EdgeIndex Solver::MinPlusMinusExternalEdge(int queue_index) {
    while (!queues.edge_heaps[queue_index]->Empty()) {
        EdgeIndex edge = queues.edge_heaps[queue_index]->Top();

        NodeIndex head = edges.head[edge.index]; // plus
        NodeIndex tail = edges.tail[edge.index]; // minus
        if (!nodes.plus[head.index]) {
            std::swap(head, tail);
        }
        if (nodes.plus[head.index] && !nodes.plus[tail.index] && (nodes.tree[head.index] != nodes.tree[tail.index]) &&
            nodes.tree[tail.index]) {
            return edge;
        }
        queues.edge_heaps[queue_index]->Pop();
        edges.queue_index[edge.index] = -1;
    }
    return EdgeIndex(-1);
}

Solver::NodeIndex Solver::MinMinusBlossom(int queue_index) {
    while (!queues.node_heaps[queue_index]->Empty()) {
        NodeIndex node = queues.node_heaps[queue_index]->Top();
        if (nodes.is_alive[node.index] && nodes.tree[node.index] && !nodes.plus[node.index] && !nodes.blossom_parent[node.index]) {
            return node;
        }
        queues.node_heaps[queue_index]->Pop();
        nodes.queue_index[node.index] = -1;
    }
    return NodeIndex(-1);
}

void Solver::RotateReceptacle(NodeIndex blossom, NodeIndex child) {
    // update the structure of matched edges inside the blossom and at the new receptacle

    if (nodes.blossom_children[blossom.index].empty()) {
        throw std::runtime_error("In RotateReceptacle: the vertex must be a blossom");
    }
    if (nodes.blossom_parent[child.index] != blossom) {
        throw std::runtime_error("In RotateReceptacle: the new receptacle is not a blossom child of this node");
    }

    NodeIndex cur_node = OtherEnd(nodes.blossom_brother_clockwise[child.index], child);
    if (!cur_node) {
        throw std::runtime_error("Invalid cur_node");
    }

    bool match = false;
    while (cur_node != child) {
        if (match) {
            MakeEdgeMatched(nodes.blossom_brother_anticlockwise[cur_node.index]);
        } else {
            MakeEdgeUnmatched(nodes.blossom_brother_anticlockwise[cur_node.index]);
        }
        match = !match;
        cur_node = OtherEnd(nodes.blossom_brother_clockwise[cur_node.index], cur_node);
    }
    MakeEdgeUnmatched(nodes.blossom_brother_anticlockwise[child.index]);
    nodes.matched_edge[child.index] = nodes.matched_edge[blossom.index]; // never a nullptr
}

Solver::NodeIndex Solver::DeeperNode(NodeIndex blossom, EdgeIndex edge) const {
    // returns a node that is a blossom child of blossom and is adjacent to edge
    // TODO do something smarter maybe

    for (int i = 0; i + 1 < static_cast<int>(nodes.children_neighbors_boundaries[blossom.index].size()); ++i) {
        auto it_left = nodes.children_neighbors_boundaries[blossom.index][i];
        auto it_right = nodes.children_neighbors_boundaries[blossom.index][i + 1];

        for (auto it = it_left; it != it_right; ++it) {
            if (*it == edge) {
                return nodes.blossom_children[blossom.index][i];
            }
        }
    }
    return nodes.blossom_children[blossom.index].back();
}

void Solver::UpdateQueuesAfterGrow(NodeIndex child, NodeIndex grandchild) {
    TreeIndex tree = nodes.tree[child.index];

    if (!IsElementary(child)) {
        AddNodeToQueue(child, trees.minus_blossoms[tree.index]);
    }

    for (EdgeIndex edge_to_neighbor : nodes.neighbors[grandchild.index]) {
        NodeIndex grandchild_neighbor = OtherEnd(edge_to_neighbor, grandchild);
        if (grandchild_neighbor == grandchild) {
            continue;
        }

        if (!nodes.tree[grandchild_neighbor.index]) {
            AddEdgeToQueue(edge_to_neighbor, trees.plus_empty_edges[tree.index]);
            continue;
        }
        if (nodes.tree[grandchild_neighbor.index] == tree && nodes.plus[grandchild_neighbor.index]) {
            AddEdgeToQueue(edge_to_neighbor, trees.plus_plus_internal_edges[tree.index]);
            continue;
        }

        if (nodes.tree[grandchild_neighbor.index] != tree) {
            if (nodes.plus[grandchild_neighbor.index]) {
                AddPQPlusPlus(tree, nodes.tree[grandchild_neighbor.index], edge_to_neighbor);
            } else {
                AddPQPlusMinus(tree, nodes.tree[grandchild_neighbor.index], edge_to_neighbor);
            }
        }
    }

    for (EdgeIndex edge_to_neighbor : nodes.neighbors[child.index]) {
        NodeIndex child_neighbor = OtherEnd(edge_to_neighbor, child);
        if (child_neighbor == child) {
            continue;
        }

        if ((nodes.tree[child_neighbor.index] == tree) && nodes.plus[child_neighbor.index]) {
            // the edge is no longer plus empty
            RemoveEdgeFromQueue(edge_to_neighbor);
            continue;
        }

        if (nodes.plus[child_neighbor.index] && (nodes.tree[child_neighbor.index] != tree)) {
            // child_neighbor is in some other tree
            AddPQPlusMinus(nodes.tree[child_neighbor.index], tree, edge_to_neighbor);
        }
    }
}

void Solver::UpdateQueuesBeforeShrink(const std::vector<NodeIndex> &minus_children) {
    if (minus_children.empty()) {
        return;
    }
    TreeIndex tree = nodes.tree[minus_children.front().index];

    for (NodeIndex child : minus_children) {
        for (EdgeIndex edge : nodes.neighbors[child.index]) {
            NodeIndex neighbor = OtherEnd(edge, child);

            if (!nodes.tree[neighbor.index]) {
                AddEdgeToQueue(edge, trees.plus_empty_edges[tree.index]);
                continue;
            }
            if ((nodes.tree[neighbor.index] == tree) && nodes.plus[neighbor.index]) {
                AddEdgeToQueue(edge, trees.plus_plus_internal_edges[tree.index]);
                continue;
            }

            if (nodes.tree[neighbor.index] != tree) {
                if (nodes.plus[neighbor.index]) {
                    AddPQPlusPlus(tree, nodes.tree[neighbor.index], edge);
                } else {
                    AddPQPlusMinus(tree, nodes.tree[neighbor.index], edge);
                }
            }
        }
    }
}

void Solver::RemoveLoopsFromQueues(NodeIndex blossom) {
    for (EdgeIndex edge : nodes.neighbors[blossom.index]) {
        if (!IsTopBlossom(edges.head[edge.index])) {
            RemoveEdgeFromQueue(edge);
        }
    }
}

void Solver::UpdateQueuesAfterExpand(TreeIndex tree, NodeIndex blossom, const std::vector<NodeIndex> &children) {
    RemoveNodeFromQueue(blossom);

    for (NodeIndex child : children) {
        if (nodes.tree[child.index] && !nodes.plus[child.index] && !IsElementary(child)) {
            AddNodeToQueue(child, trees.minus_blossoms[tree.index]);
        }

        if (nodes.plus[child.index]) {
            for (EdgeIndex edge : nodes.neighbors[child.index]) {
                NodeIndex other_end = OtherEnd(edge, child);

                if (!nodes.tree[other_end.index]) {
                    AddEdgeToQueue(edge, trees.plus_empty_edges[tree.index]);
                    continue;
                }
                if (nodes.tree[other_end.index] == tree && nodes.plus[other_end.index]) {
                    AddEdgeToQueue(edge, trees.plus_plus_internal_edges[tree.index]);
                    continue;
                }
                if (nodes.tree[other_end.index] && (nodes.tree[other_end.index] != tree)) {
                    // other_end is in another tree
                    if (nodes.plus[other_end.index]) {
                        AddPQPlusPlus(tree, nodes.tree[other_end.index], edge);
                    } else {
                        AddPQPlusMinus(tree, nodes.tree[other_end.index], edge);
                    }
                }
            }
        }

        if (!nodes.tree[child.index]) {
            for (EdgeIndex edge : nodes.neighbors[child.index]) {
                RemoveEdgeFromQueue(edge);
                NodeIndex other_end = OtherEnd(edge, child);
                if (nodes.plus[other_end.index]) {
                    AddEdgeToQueue(edge, trees.plus_empty_edges[nodes.tree[other_end.index].index]);
                }
            }
        }
    }
}

int Solver::TreeTreeQueueIndex(TreeIndex other_tree, const std::vector<std::pair<TreeIndex, int> > &tree_neighbors) {
    // TODO consider traversing from back to front

    for (auto [tree, queue_index] : tree_neighbors) {
        if (tree == other_tree) {
            return queue_index;
        }
    }
    return -1;
}

void Solver::AddPQPlusPlus(TreeIndex first, TreeIndex second, EdgeIndex edge) {
    int queue_index = TreeTreeQueueIndex(second, trees.pq_plus_plus[first.index]);
    if (queue_index > 0) {
        AddEdgeToQueue(edge, queue_index);
    } else {
        queues.edge_heaps.emplace_back(std::make_unique<EdgeHeap>(params.heap_arity, EdgeComparator{this}, EdgeValidator{this}));
        queue_index = static_cast<int>(queues.edge_heaps.size()) - 1;
        AddEdgeToQueue(edge, queue_index);

        trees.pq_plus_plus[first.index].emplace_back(second, queue_index);
        trees.pq_plus_plus[second.index].emplace_back(first, queue_index);
    }
}

void Solver::AddPQPlusMinus(TreeIndex tree_plus, TreeIndex tree_minus, EdgeIndex edge) {
    int queue_index = TreeTreeQueueIndex(tree_minus, trees.pq_plus_minus[tree_plus.index]);
    if (queue_index > 0) {
        AddEdgeToQueue(edge, queue_index);
    } else {
        queues.edge_heaps.emplace_back(std::make_unique<EdgeHeap>(params.heap_arity, EdgeComparator{this}, EdgeValidator{this}));
        queue_index = static_cast<int>(queues.edge_heaps.size()) - 1;
        AddEdgeToQueue(edge, queue_index);

        trees.pq_plus_minus[tree_plus.index].emplace_back(tree_minus, queue_index);
        trees.pq_minus_plus[tree_minus.index].emplace_back(tree_plus, queue_index);
    }
}

/*void Solver::ValidateQueues() {
    // checks that the queues are in a correct state, throws if not

    // every queue must hold only the edges that belong to this queue

    for (int tree_index = 0; tree_index < static_cast<int>(trees.is_alive.size()); ++tree_index) {
        if (!trees.is_alive[tree_index]) {
            continue;
        }

        // check plus empty
        for (EdgeIndex edge : *queues.edge_heaps[trees.plus_empty_edges[tree_index]]) {
            NodeIndex first = edges.head[edge.index];    // plus
            NodeIndex second = edges.tail[edge.index];   // empty
            if (nodes.tree[first.index].index != tree_index) {
                std::swap(first, second);
            }
            if (nodes.tree[first.index].index != tree_index || !nodes.plus[first.index]) {
                throw std::runtime_error("Incorrect plus empty queue");
            }
            if (nodes.tree[second.index]) {
                throw std::runtime_error("Incorrect plus empty queue");
            }
        }

        // check plus plus internal
        for (EdgeIndex edge : *queues.edge_heaps[trees.plus_plus_internal_edges[tree_index]]) {
            NodeIndex first = edges.head[edge.index];
            NodeIndex second = edges.tail[edge.index];

            if (nodes.tree[first.index].index != tree_index || !nodes.plus[first.index]) {
                throw std::runtime_error("Incorrect plus plus internal queue");
            }
            if (nodes.tree[second.index].index != tree_index || !nodes.plus[second.index]) {
                throw std::runtime_error("Incorrect plus plus internal queue");
            }
        }

        // check minus blossoms
        for (NodeIndex node : *queues.node_heaps[trees.minus_blossoms[tree_index]]) {
            if (IsElementary(node) || (nodes.tree[node.index].index != tree_index) || nodes.plus[node.index]) {
                throw std::runtime_error("Incorrect minus blossom queue");
            }
        }

        // check plus plus external
        for (auto [other_tree, queue_index] : trees.pq_plus_plus[tree_index]) {
            if (!trees.is_alive[other_tree.index]) {
                continue;
            }
            for (EdgeIndex edge : *queues.edge_heaps[queue_index]) {
                NodeIndex first = edges.head[edge.index];    // in this tree
                NodeIndex second = edges.tail[edge.index];   // in the other tree
                if (nodes.tree[first.index].index != tree_index) {
                    std::swap(first, second);
                }
                if (nodes.tree[first.index].index != tree_index || !nodes.plus[first.index]) {
                    throw std::runtime_error("Incorrect plus plus queue");
                }
                if (nodes.tree[second.index] != other_tree || !nodes.plus[second.index]) {
                    throw std::runtime_error("Incorrect plus plus queue");
                }
            }
        }

        // check plus minus external
        for (auto [other_tree, queue_index] : trees.pq_plus_minus[tree_index]) {
            if (!trees.is_alive[other_tree.index]) {
                continue;
            }
            for (EdgeIndex edge : *queues.edge_heaps[queue_index]) {
                NodeIndex first = edges.head[edge.index];    // in this tree
                NodeIndex second = edges.tail[edge.index];   // in the other tree
                if (nodes.tree[first.index].index != tree_index) {
                    std::swap(first, second);
                }
                if (nodes.tree[first.index].index != tree_index || !nodes.plus[first.index]) {
                    throw std::runtime_error("Incorrect plus minus queue");
                }
                if (nodes.tree[second.index] != other_tree || nodes.plus[second.index]) {
                    throw std::runtime_error("Incorrect plus minus queue");
                }
            }
        }

        // check minus plus external
        for (auto [other_tree, queue_index] : trees.pq_minus_plus[tree_index]) {
            if (!trees.is_alive[other_tree.index]) {
                continue;
            }
            for (EdgeIndex edge : *queues.edge_heaps[queue_index]) {
                NodeIndex first = edges.head[edge.index];    // in this tree
                NodeIndex second = edges.tail[edge.index];   // in the other tree
                if (nodes.tree[first.index].index != tree_index) {
                    std::swap(first, second);
                }
                if (nodes.tree[first.index].index != tree_index || nodes.plus[first.index]) {
                    throw std::runtime_error("Incorrect minus plus queue");
                }
                if (nodes.tree[second.index] != other_tree || !nodes.plus[second.index]) {
                    throw std::runtime_error("Incorrect minus plus queue");
                }
            }
        }
    }
}*/

void Solver::ValidatePositiveSlacks() {
    for (int edge_index = 0; edge_index < static_cast<int>(edges.head.size()); ++edge_index) {
        if (SlackQuadrupled(EdgeIndex(edge_index)) < 0) {
            std::cout << "edge " << edges.head[edge_index].index << " " << edges.tail[edge_index].index << std::endl;
            std::cout << "blossom parents: " << nodes.blossom_parent[edges.head[edge_index].index].index << " " << nodes.blossom_parent[edges.tail[edge_index].index].index << std::endl;
            std::cout << "in trees: " << nodes.tree[edges.head[edge_index].index].index << " " << nodes.tree[edges.tail[edge_index].index].index << std::endl;
            std::cout << "pluses: " << nodes.plus[edges.head[edge_index].index] << " " << nodes.plus[edges.tail[edge_index].index] << std::endl;

            throw std::runtime_error("Negative slack edge found");
        }
    }
}

void Solver::ComputeMatching() {
    // must be called after the blossoms are dissolved
    if (!matching.empty()) {
        throw std::runtime_error("ComputeMatching: matching is already non-empty");
    }

    matching.reserve(num_vertices_elementary / 2);

    for (int i = 0; i < static_cast<int>(edges.matched.size()); ++i) {
        if (edges.matched[i]) {
            matching.emplace_back(edges.head[i].index, edges.tail[i].index);
        }
    }
}

void Solver::ComputeDualCertificate() {
    // can be called only before we destroy the blossoms in the end

    if (!dual_certificate.empty()) {
        throw std::runtime_error("Dual certificate is already non-empty");
    }

    std::vector<int> node_index_alive(nodes.is_alive.size(), -1);
    dual_certificate.reserve(nodes.is_alive.size());
    int index = 0;
    for (int i = 0; i < static_cast<int>(nodes.is_alive.size()); ++i) {
        if (nodes.is_alive[i]) {
            node_index_alive[i] = index;
            dual_certificate.emplace_back(index, DualVariableQuadrupled(NodeIndex(i)), -1);
            for (NodeIndex child : nodes.blossom_children[i]) {
                std::get<2>(dual_certificate[node_index_alive[child.index]]) = index;
            }
            ++index;
        }
    }
}

void Solver::DestroyBlossoms() {
    for (int i = static_cast<int>(nodes.is_alive.size()) - 1; i >= num_vertices_elementary; --i) {
        if (nodes.is_alive[i]) {
            if (!IsTopBlossom(NodeIndex(i))) {
                throw std::runtime_error("In SolverWeighted::RotateMatchingInsideBlossoms: not going top down");
            }

            Dissolve(NodeIndex(i));
        }
    }
}

int64_t Solver::DualObjectiveQuadrupled() const {
    int64_t objective = 0;

    for (int i = 0; i < static_cast<int>(nodes.is_alive.size()); ++i) {
        if (nodes.is_alive[i]) {
            objective += DualVariableQuadrupled(NodeIndex(i));
        }
    }

    return objective;
}

int64_t Solver::PrimalObjective() const {
    int objective = 0;

    for (int i = 0; i < static_cast<int>(edges.matched.size()); ++i) {
        if (edges.matched[i]) {
            objective += edges.weight[i];
        }
    }

    return objective;
}

bool Solver::IsElementary(const NodeIndex node) const {
    return node.index < num_vertices_elementary;
}

bool Solver::IsTopBlossom(const NodeIndex node) const {
    return nodes.blossom_parent[node.index].index < 0;
}

Solver::NodeIndex Solver::TopBlossom(NodeIndex node) const {
    NodeIndex cur_vertex = node;
    while (nodes.blossom_parent[cur_vertex.index]) {
        cur_vertex = nodes.blossom_parent[cur_vertex.index];
    }
    return cur_vertex;
}

int Solver::NodeDepth(NodeIndex node) const {
    int depth = 0;
    NodeIndex cur_vertex = node;
    while (nodes.blossom_parent[cur_vertex.index]) {
        cur_vertex = nodes.blossom_parent[cur_vertex.index];
        ++depth;
    }
    return depth;
}

int Solver::DualVariableQuadrupled(NodeIndex node) const {
    if (!nodes.is_alive[node.index]) {
        throw std::runtime_error("DualVariableQuadrupled: node is not alive");
    }

    if (!nodes.tree[node.index]) {
        return nodes.dual_var_quadrupled_amortized[node.index];
    }
    if (nodes.blossom_parent[node.index]) {
        return nodes.dual_var_quadrupled_amortized[node.index];
    }
    if (nodes.plus[node.index]) {
        return nodes.dual_var_quadrupled_amortized[node.index] + trees.dual_var_quadrupled[nodes.tree[node.index].
            index];
    }
    return nodes.dual_var_quadrupled_amortized[node.index] - trees.dual_var_quadrupled[nodes.tree[node.index].index];
}

Solver::NodeIndex Solver::LCA(NodeIndex first, NodeIndex second) const {
    // TODO make better

    if (nodes.tree[first.index] != nodes.tree[second.index]) {
        throw std::runtime_error("In LCA: first and second are in different trees");
    }
    if (!nodes.tree[first.index]) {
        throw std::runtime_error("In LCA: node is not in any tree");
    }

    NodeIndex root_blossom = TopBlossom(trees.root[nodes.tree[first.index].index]);

    std::unordered_set<int> visited;
    visited.insert(first.index);
    while (first != root_blossom) {
        first = OtherEnd(nodes.tree_parent[first.index], first);
        visited.insert(first.index);
    }

    while (second != root_blossom) {
        if (visited.contains(second.index)) {
            return second;
        }
        second = OtherEnd(nodes.tree_parent[second.index], second);
    }

    return root_blossom;
}

std::vector<Solver::EdgeIndex> Solver::PathToRoot(TreeIndex tree, NodeIndex node) {
    NodeIndex cur_vertex = node;
    std::vector<EdgeIndex> path;
    while (nodes.tree_parent[cur_vertex.index]) {
        path.push_back(nodes.tree_parent[cur_vertex.index]);
        cur_vertex = OtherEnd(nodes.tree_parent[cur_vertex.index], cur_vertex);
    }

    return path;
}

int Solver::SlackQuadrupled(EdgeIndex edge) const {
    return edges.slack_quadrupled_amortized[edge.index] - DualVariableQuadrupled(edges.head[edge.index]) -
        DualVariableQuadrupled(edges.tail[edge.index]);
}

Solver::NodeIndex Solver::OtherEnd(EdgeIndex edge, NodeIndex node) const {
    if (node == edges.tail[edge.index]) {
        return edges.head[edge.index];
    }
    if (node == edges.head[edge.index]) {
        return edges.tail[edge.index];
    }
    // the edge is a loop
    return node;
}

Solver::NodeIndex Solver::SharedNode(EdgeIndex first_edge, EdgeIndex second_edge) {
    // returns the max blossom adjacent to both edges

    NodeIndex first_head = edges.head[first_edge.index];
    NodeIndex first_tail = edges.tail[first_edge.index];
    NodeIndex second_head = edges.head[second_edge.index];
    NodeIndex second_tail = edges.tail[second_edge.index];

    if ((first_head == second_head) || (first_head == second_tail)) {
        return first_head;
    }
    if ((first_tail == second_head) || (first_tail == second_tail)) {
        return first_tail;
    }
    throw std::runtime_error("In Node::SharedNode: the edges don't share an endpoint");
}

void Solver::MakeEdgeMatched(EdgeIndex edge) {
    nodes.matched_edge[edges.head[edge.index].index] = edge;
    nodes.matched_edge[edges.tail[edge.index].index] = edge;
    edges.matched[edge.index] = true;
}

void Solver::MakeEdgeUnmatched(EdgeIndex edge) {
    edges.matched[edge.index] = false;
}

void Solver::AddEdgeToQueue(EdgeIndex edge, int queue_index) {
    // removes the edge from the current queue (if exists), adds to the new one

    // if (edges.queue_index[edge.index] == queue_index) {
    //     return;
    // }

    if (edges.queue_index[edge.index] >= 0) {
        queues.edge_heaps[edges.queue_index[edge.index]]->Erase(edges.handle[edge.index]);
    }
    edges.handle[edge.index] = queues.edge_heaps[queue_index]->Push(edge);
    edges.queue_index[edge.index] = queue_index;
}

void Solver::RemoveEdgeFromQueue(EdgeIndex edge) {
    // removes the edge from the current queue (if exists)
    if (edges.queue_index[edge.index] >= 0) {
        queues.edge_heaps[edges.queue_index[edge.index]]->Erase(edges.handle[edge.index]);
        edges.queue_index[edge.index] = -1;
    }
}

void Solver::AddNodeToQueue(NodeIndex node, int queue_index) {
    // if (nodes.queue_index[node.index] == queue_index) {
    //     return;
    // }

    if (nodes.queue_index[node.index] >= 0) {
        queues.node_heaps[nodes.queue_index[node.index]]->Erase(nodes.handle[node.index]);
    }
    nodes.handle[node.index] = queues.node_heaps[queue_index]->Push(node);
    nodes.queue_index[node.index] = queue_index;
}

void Solver::RemoveNodeFromQueue(NodeIndex node) {
    if (nodes.queue_index[node.index] >= 0) {
        queues.node_heaps[nodes.queue_index[node.index]]->Erase(nodes.handle[node.index]);
        nodes.queue_index[node.index] = -1;
    }
}

Solver::TreeIndex Solver::NextAliveTree(TreeIndex current) {
    TreeIndex ans = trees.next_alive_tree[current.index];
    if (trees.is_alive[ans.index]) {
        return ans;
    }

    while (!trees.is_alive[ans.index]) {
        ans = trees.next_alive_tree[ans.index];
    }
    trees.next_alive_tree[current.index] = ans;
    return ans;
}

Solver::TreeIndex Solver::FirstAliveTree() {
    TreeIndex ans = first_alive_tree;
    if (trees.is_alive[ans.index]) {
        return ans;
    }

    while (!trees.is_alive[ans.index]) {
        ans = trees.next_alive_tree[ans.index];
    }
    first_alive_tree = ans;
    return ans;
}

void Solver::UpdateAliveIndices() {
    if (num_trees_alive == 0) {
        return;
    }

    TreeIndex tree = FirstAliveTree();
    int index = 0;
    do {
        trees.alive_index[tree.index] = index;
        ++index;
        tree = NextAliveTree(tree);
    } while (tree != FirstAliveTree());
}

int Solver::PlusEmptySlack(TreeIndex tree) {
    EdgeIndex edge = MinPlusEmptyEdge(trees.plus_empty_edges[tree.index]);
    if (edge) {
        return SlackQuadrupled(edge);
    }
    return INT32_MAX;
}

int Solver::PlusPlusInternalSlack(TreeIndex tree) {
    EdgeIndex edge = MinPlusPlusInternalEdge(trees.plus_plus_internal_edges[tree.index]);
    if (edge) {
        return SlackQuadrupled(edge);
    }
    return INT32_MAX;
}

int Solver::MinMinusBlossomVariable(TreeIndex tree) {
    NodeIndex node = MinMinusBlossom(trees.minus_blossoms[tree.index]);
    if (node) {
        return DualVariableQuadrupled(node);
    }
    return INT32_MAX;
}

std::vector<std::pair<Solver::TreeIndex, int> > Solver::PlusPlusExternalSlacks(TreeIndex tree) {
    std::vector<std::pair<TreeIndex, int> > result;
    result.reserve(trees.pq_plus_plus[tree.index].size());
    for (auto [tree_neighbor, queue_index] : trees.pq_plus_plus[tree.index]) {
        if (trees.is_alive[tree_neighbor.index]) {
            EdgeIndex edge = MinPlusPlusExternalEdge(queue_index);
            if (edge) {
                result.emplace_back(tree_neighbor, SlackQuadrupled(queues.edge_heaps[queue_index]->Top()));
            }
        }
    }
    return result;
}

std::vector<std::pair<Solver::TreeIndex, int> > Solver::PlusMinusExternalSlacks(TreeIndex tree) {
    std::vector<std::pair<TreeIndex, int> > result;
    result.reserve(trees.pq_plus_minus[tree.index].size());
    for (auto [tree_neighbor, queue_index] : trees.pq_plus_minus[tree.index]) {
        if (trees.is_alive[tree_neighbor.index]) {
            EdgeIndex edge = MinPlusMinusExternalEdge(queue_index);
            if (edge) {
                result.emplace_back(tree_neighbor, SlackQuadrupled(queues.edge_heaps[queue_index]->Top()));
            }
        }
    }
    return result;
}

int Solver::InitNumVertices(const std::vector<std::tuple<int, int, int> > &edge_list_) {
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