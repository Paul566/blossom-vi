#include "MWPMSolver.h"

#include <algorithm>
#include <queue>

MWPMSolver::MWPMSolver(const std::vector<std::tuple<int, int, int> > &edge_list_,
                         const SolverParameters &params_) : primal_objective(INT64_MAX), dual_objective(INT64_MIN),
                                                            params(params_),
                                                            num_vertices_elementary(InitNumVertices(edge_list_)),
                                                            current_round(0), actionable_edges_head(0),
                                                            actionable_nodes_head(0) {
    // TODO check validity of edge_list_

    nodes.reserve(2 * num_vertices_elementary);
    blossom_parents.reserve(2 * num_vertices_elementary);
    blossom_ancestors.reserve(2 * num_vertices_elementary);
    adj_list.reserve(2 * num_vertices_elementary);
    blossom_structures.reserve(2 * num_vertices_elementary);
    node_heap_infos.reserve(2 * num_vertices_elementary);
    for (int i = 0; i < num_vertices_elementary; ++i) {
        nodes.emplace_back(Node(i));
        blossom_parents.emplace_back(-1);
        blossom_ancestors.emplace_back(-1);
        adj_list.emplace_back();
        blossom_structures.emplace_back();
        node_heap_infos.emplace_back(-1, -1, -1, -1, 0);
    }

    // reserve for adj list
    std::vector<int> degrees(num_vertices_elementary, 0);
    for (auto [from, to, weight] : edge_list_) {
        ++degrees[from];
        ++degrees[to];
    }
    for (int i = 0; i < num_vertices_elementary; ++i) {
        adj_list[i].reserve(degrees[i]);
    }

    // rearrange the edges
    std::vector<std::vector<std::pair<int, int> > > incident_edges(num_vertices_elementary);
    for (int i = 0; i < num_vertices_elementary; ++i) {
        incident_edges[i].reserve(degrees[i]);
    }
    for (int i = 0; i < static_cast<int>(edge_list_.size()); ++i) {
        int head = std::get<0>(edge_list_[i]);
        int tail = std::get<1>(edge_list_[i]);
        incident_edges[head].emplace_back(tail, std::get<2>(edge_list_[i]));
    }

    edges.reserve(edge_list_.size());
    edge_weights.reserve(edge_list_.size());
    elementary_heads.reserve(edge_list_.size());
    elementary_tails.reserve(edge_list_.size());
    matched = std::vector<uint8_t>(edge_list_.size(), false);
    // maybe_has_zero_slack = std::vector<uint8_t>(edge_list_.size(), true);
    int i = 0;
    for (int from = 0; from < num_vertices_elementary; ++from) {
        for (auto [to, weight] : incident_edges[from]) {
            edges.emplace_back(Edge(from, to, weight));
            edge_weights.push_back(weight);
            elementary_heads.push_back(from);
            elementary_tails.push_back(to);

            adj_list[from].emplace_back(i * 2);
            adj_list[to].emplace_back(i * 2 + 1);

            ++i;
        }
    }

    primal_update_record.reserve(num_vertices_elementary);
    actionable_edges.reserve(edge_list_.size());
    actionable_nodes.reserve(num_vertices_elementary);

    nodes_label_cnt = 1;
    num_trees_alive = INT32_MAX;

    aux_counter1 = 0;
    aux_counter2 = 0;
    aux_counter3 = 0;
    aux_counter4 = 0;

    even_path_tmp.reserve(num_vertices_elementary);
    odd_path_tmp.reserve(num_vertices_elementary);
    path_to_root.reserve(num_vertices_elementary);
    traversal_nodes_tmp.reserve(num_vertices_elementary);
    traversal_lists_tmp.reserve(num_vertices_elementary);
    edges_stack_tmp.reserve(64);
    organize_sizes_tmp.reserve(num_vertices_elementary);
    organize_label_to_index_tmp.reserve(num_vertices_elementary);
    future_blossoms_tmp.reserve(num_vertices_elementary >> 1);
}

void MWPMSolver::FindMinPerfectMatching() {
    Init();
    InitializeTrees();

    if (params.verbose) {
        PrintGraph();
    }

    while (!alive_trees.empty()) {
        if (params.print_statistics) {
            std::cout << "round " << current_round << std::endl;
            std::cout << "trees left: " << alive_trees.size() << std::endl;

            // int avgpp = 0;
            // int avgpm = 0;
            // for (int tree : alive_trees) {
            //     avgpp += tree_heap_infos[tree].pq_plus_plus.size();
            //     avgpm += tree_heap_infos[tree].pq_plus_minus.size();
            // }
            // std::cout << avgpp * 1. / num_trees_alive << " " << avgpm * 1. / num_trees_alive << std::endl;
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

    if (params.print_statistics) {
        std::cout << current_round << " " << primal_objective << std::endl;
        std::cout << aux_counter1 << " " << aux_counter2 << " " << aux_counter3 << " " << aux_counter4 << std::endl;
    }

    if (params.verbose) {
        std::cout << "Primal objective:\t" << primal_objective << std::endl;
    }
}

const std::vector<std::pair<int, int> > &MWPMSolver::Matching() const {
    return matching;
}

const std::vector<std::tuple<int, int, int> > &MWPMSolver::DualCertificate() const {
    if (!params.compute_dual_certificate) {
        throw std::runtime_error("Parameter \"compute_dual_certificate\" is set to false");
    }
    return dual_certificate;
}

MWPMSolver::Edge::Edge(int head_, int tail_, int weight_) : head(head_), tail(tail_),
                                                             last_round_updated(-1),
                                                             queue_index(-1),
                                                             heap_child(-1), heap_next(-1), heap_prev(-1) {
    slack_quadrupled_amortized_ = 4 * weight_;
}

MWPMSolver::Node::Node(int index_) : matched_edge(-1), minus_parent(-1), receptacle_(index_), tree(-1),
                                      old_tree(-1), tree_var_at_birth(0),
                                      slack_diff(0),
                                      is_in_record(false) {
    is_alive = true;
    plus = false;
    old_plus = false;
    label = 0;
}

MWPMSolver::Tree::Tree(int root_) {
    is_alive = true;
    dual_var_quadrupled = 0;
    alive_index = 0;
}

void MWPMSolver::PrintGraph() const {
    std::cout << "Adjacency list (to, matched, weight, amortized slack):" << std::endl;

    for (int i = 0; i < num_vertices_elementary; ++i) {
        std::cout << i << ": ";
        for (ArcIndex arc : adj_list[i]) {
            std::cout << "(" << OtherElementaryEnd(arc) << " " << static_cast<bool>(matched[arc.index >> 1]) << " " <<
                edge_weights[arc.index >> 1] << " " << edges[arc.index >> 1].slack_quadrupled_amortized_ << ") ";
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
        std::cout << "root: " << roots[tree] << " var: " << trees[tree].dual_var_quadrupled / 4. << std::endl;
    }
}

void MWPMSolver::PrintNode(int node) const {
    std::cout << node << ": y_v = " << DualVariableQuadrupled(node) / 4.;
    if (blossom_parents[node] >= 0) {
        std::cout << " blossom_parent: " << blossom_parents[node];
    }
    std::cout << std::endl;
}

void MWPMSolver::Init() {
    InitMakeSlacksNonnegative();

    InitGreedyIncreaseVars();

    // InitFindLengthThreeAugmentations();
    FractionalMatchingInit();
}

void MWPMSolver::InitMakeSlacksNonnegative() {
    for (int i = 0; i < num_vertices_elementary; ++i) {
        if (adj_list[i].empty()) {
            std::cout << "vertex " << i << ": no neighbors" << std::endl;
            throw std::runtime_error("Found an isolated vertex => no perfect matching exists");
        }

        int min_weight = edge_weights[adj_list[i].front().index >> 1];
        for (ArcIndex arc : adj_list[i]) {
            if (edge_weights[arc.index >> 1] < min_weight) {
                min_weight = edge_weights[arc.index >> 1];
            }
        }

        node_heap_infos[i].dual_var_quadrupled_amortized_ += 2 * min_weight;
        for (ArcIndex arc : adj_list[i]) {
            edges[arc.index >> 1].slack_quadrupled_amortized_ -= 2 * min_weight;
        }
    }
}

void MWPMSolver::InitGreedyIncreaseVars() {
    for (int i = 0; i < num_vertices_elementary; ++i) {
        if (nodes[i].matched_edge.index >= 0) {
            continue;
        }

        int smallest_slack = INT32_MAX;
        for (ArcIndex arc : adj_list[i]) {
            if (edges[arc.index >> 1].slack_quadrupled_amortized_ < smallest_slack) {
                smallest_slack = edges[arc.index >> 1].slack_quadrupled_amortized_;
            }
        }

        if (smallest_slack < 0) {
            throw std::runtime_error("InitGreedyIncreaseVars: smallest_slack < 0");
        }
        node_heap_infos[i].dual_var_quadrupled_amortized_ += smallest_slack;
        for (ArcIndex arc : adj_list[i]) {
            edges[arc.index >> 1].slack_quadrupled_amortized_ -= smallest_slack;

            // if (edges[arc.index >> 1].slack_quadrupled_amortized_ == 0) {
            //     maybe_has_zero_slack[arc.index >> 1] = true;
            // } else {
            //     maybe_has_zero_slack[arc.index >> 1] = false;
            // }
        }

        for (ArcIndex arc : adj_list[i]) {
            if (edges[arc.index >> 1].slack_quadrupled_amortized_ == 0) {
                if (nodes[OtherElementaryEnd(arc)].matched_edge.index < 0) {
                    // if the other vertex is also unmatched, match the edge
                    matched[arc.index >> 1] = true;
                    nodes[i].matched_edge = arc;
                    nodes[OtherElementaryEnd(arc)].matched_edge = ReverseArc(arc);
                    break;
                }
            }
        }
    }
}

void MWPMSolver::InitFindLengthThreeAugmentations() {
    for (int root = 0; root < num_vertices_elementary; ++root) {
        if (nodes[root].matched_edge.index >= 0) {
            continue;
        }

        int child = -1;
        int edge_to_child = -1;
        for (ArcIndex arc : adj_list[root]) {
            int slack = edges[arc.index >> 1].slack_quadrupled_amortized_;
            if (slack == 0) {
                edge_to_child = arc.index >> 1;
                child = OtherElementaryEnd(arc);
                break;
            }
        }

        if (child < 0) {
            throw std::runtime_error("InitFindLengthThreeAugmentations: no child");
            // TODO this is possible?
        }
        int grandchild = OtherElementaryEnd(nodes[child].matched_edge);
        for (ArcIndex arc : adj_list[grandchild]) {
            int edge = arc.index >> 1;
            if (!matched[edge]) {
                int slack = edges[edge].slack_quadrupled_amortized_;
                if (slack == 0) {
                    int final_end = OtherElementaryEnd(arc);
                    if (final_end == root || nodes[final_end].matched_edge.index >= 0) {
                        continue;
                    }

                    // augment
                    MakeEdgeUnmatched(nodes[child].matched_edge.index >> 1);
                    MakeEdgeMatched(edge_to_child);
                    MakeEdgeMatched(edge);
                    break;
                }
            }
        }
    }
}

void MWPMSolver::FractionalMatchingInit() {
    // TODO use a different slack amortization

    if (params.verbose) {
        PrintGraph();
    }

    half_int_clockwise = std::vector<ArcIndex>(num_vertices_elementary, {-1});
    int tree_cnt = 0;
    init_tree_nodes.resize(params.init_max_tree_size);

    for (int root = 0; root < num_vertices_elementary; ++root) {
        if (nodes[root].matched_edge.index >= 0 || half_int_clockwise[root].index >= 0) {
            continue;
        }

        ++tree_cnt; // marking the nodes in the tree
        nodes[root].tree = tree_cnt;
        nodes[root].plus = true;
        init_tree_nodes.clear();
        init_tree_nodes.emplace_back(root);

        int tree_var = 0;
        init_plus_empty = {};
        init_min_plus_plus_edge = -1;
        init_min_plus_plus_slack_amortized = INT32_MAX; // must be correct

        // might contain (+, +) or (+, -) arcs
        for (ArcIndex arc : adj_list[root]) {
            init_plus_empty.push({edges[arc.index >> 1].slack_quadrupled_amortized_, arc});
        }

        while (init_tree_nodes.size() < params.init_max_tree_size) {
            if (init_plus_empty.empty()) {
                break;
            }
            auto [slack_amortized, arc] = init_plus_empty.top();
            init_plus_empty.pop();

            if (slack_amortized < init_min_plus_plus_slack_amortized / 2) {
                int neighbor = OtherElementaryEnd(arc);
                if (nodes[neighbor].tree == tree_cnt) {
                    // in this tree
                    if (nodes[neighbor].plus) {
                        // ++
                        throw std::runtime_error("FractionalMatchingInit: ++ edge not discovered");
                    } else {
                        // +-
                        continue;
                    }
                } else {
                    // not in this tree
                    if (nodes[neighbor].matched_edge.index >= 0) {
                        tree_var = slack_amortized;
                        if (params.verbose) {
                            std::cout << "tree var: " << tree_var << std::endl;
                        }
                        FracInitGrow(arc, tree_cnt, tree_var);
                    } else {
                        if (half_int_clockwise[neighbor].index >= 0) {
                            tree_var = slack_amortized;
                            if (params.verbose) {
                                std::cout << "tree var: " << tree_var << std::endl;
                            }
                            FracInitAugmentCycle(arc, root);
                            break;
                        } else {
                            tree_var = slack_amortized;
                            if (params.verbose) {
                                std::cout << "tree var: " << tree_var << std::endl;
                            }
                            FracInitAugment(ThisElementaryEnd(arc), root);
                            MakeEdgeMatched(arc.index >> 1);
                            break;
                        }
                    }
                }
            } else {
                tree_var = init_min_plus_plus_slack_amortized / 2;
                if (params.verbose) {
                    std::cout << "tree var: " << tree_var << std::endl;
                }
                FracInitMakeCycle(init_min_plus_plus_edge, root);
                break;
            }
        }

        FracInitCleanTree(tree_var);
    }

    FracInitCleanup();
    if (params.verbose) {
        PrintGraph();
    }
}

void MWPMSolver::FracInitGrow(ArcIndex arc, int tree_cnt, int tree_var) {
    int child = OtherElementaryEnd(arc);
    int grandchild = OtherElementaryEnd(nodes[child].matched_edge);

    if (params.verbose) {
        std::cout << "INIT GROW " << ThisElementaryEnd(arc) << " " << child << " " << grandchild << std::endl;
    }

    nodes[child].minus_parent = ReverseArc(arc);

    nodes[child].plus = false;
    nodes[grandchild].plus = true;

    nodes[child].tree = tree_cnt;
    nodes[grandchild].tree = tree_cnt;

    init_tree_nodes.push_back(child);
    init_tree_nodes.push_back(grandchild);

    node_heap_infos[child].dual_var_quadrupled_amortized_ += tree_var;
    node_heap_infos[grandchild].dual_var_quadrupled_amortized_ -= tree_var;

    for (ArcIndex grand_arc : adj_list[grandchild]) {
        edges[grand_arc.index >> 1].slack_quadrupled_amortized_ += tree_var;

        int neighbor = OtherElementaryEnd(grand_arc);
        if (nodes[neighbor].tree != tree_cnt) {
            init_plus_empty.push({edges[grand_arc.index >> 1].slack_quadrupled_amortized_, grand_arc});
        } else {
            if (nodes[neighbor].plus && edges[grand_arc.index >> 1].slack_quadrupled_amortized_ <
                init_min_plus_plus_slack_amortized) {
                init_min_plus_plus_slack_amortized = edges[grand_arc.index >> 1].slack_quadrupled_amortized_;
                init_min_plus_plus_edge = grand_arc.index >> 1;
            }
        }
    }
    for (ArcIndex child_arc : adj_list[child]) {
        edges[child_arc.index >> 1].slack_quadrupled_amortized_ -= tree_var;
    }
}

void MWPMSolver::FracInitAugment(int node_plus, int this_root) {
    if (params.verbose) {
        std::cout << "INIT AUGMENT " << node_plus << std::endl;
    }

    path_to_root.clear();
    while (node_plus != this_root) {
        path_to_root.push_back(nodes[node_plus].matched_edge.index >> 1);
        node_plus = OtherEnd(nodes[node_plus].matched_edge);
        path_to_root.push_back(nodes[node_plus].minus_parent.index >> 1);
        node_plus = OtherEnd(nodes[node_plus].minus_parent);
    }

    for (int i = 0; i < path_to_root.size(); ++i) {
        if (i % 2 == 0) {
            MakeEdgeUnmatched(path_to_root[i]);
        } else {
            MakeEdgeMatched(path_to_root[i]);
        }
    }
}

void MWPMSolver::FracInitAugmentCycle(ArcIndex arc, int this_root) {
    if (params.verbose) {
        std::cout << "INIT AUGMENT CYCLE " << ThisElementaryEnd(arc) << " " << OtherElementaryEnd(arc) << std::endl;
    }

    FracInitAugment(ThisElementaryEnd(arc), this_root);
    MakeEdgeMatched(arc.index >> 1);

    int receptacle = OtherElementaryEnd(arc);

    MakeEdgeUnmatched(half_int_clockwise[receptacle].index >> 1);
    int cur_node = OtherElementaryEnd(half_int_clockwise[receptacle]);
    half_int_clockwise[receptacle] = ArcIndex(-1);
    bool match = true;

    while (cur_node != receptacle) {
        if (match) {
            MakeEdgeMatched(half_int_clockwise[cur_node].index >> 1);
        } else {
            MakeEdgeUnmatched(half_int_clockwise[cur_node].index >> 1);
        }
        match = !match;
        int next_node = OtherElementaryEnd(half_int_clockwise[cur_node]);
        half_int_clockwise[cur_node] = ArcIndex(-1);
        cur_node = next_node;
    }
}

void MWPMSolver::FracInitMakeCycle(int edge, int this_root) {
    ++nodes_label_cnt;
    int first_vertex = elementary_heads[edge];
    int second_vertex = elementary_tails[edge];

    if (params.verbose) {
        std::cout << "INIT MAKE CYCLE " << first_vertex << " " << second_vertex << std::endl;
    }

    nodes[first_vertex].label = nodes_label_cnt;
    nodes[second_vertex].label = nodes_label_cnt;

    int lca = this_root;
    while (first_vertex != second_vertex) {
        if (nodes[first_vertex].matched_edge.index >= 0) {
            first_vertex = OtherElementaryEnd(nodes[first_vertex].matched_edge);
            first_vertex = OtherElementaryEnd(nodes[first_vertex].minus_parent);
            if (nodes[first_vertex].label == nodes_label_cnt) {
                lca = first_vertex;
                break;
            }
            nodes[first_vertex].label = nodes_label_cnt;
        }

        if (nodes[second_vertex].matched_edge.index >= 0) {
            second_vertex = OtherElementaryEnd(nodes[second_vertex].matched_edge);
            second_vertex = OtherElementaryEnd(nodes[second_vertex].minus_parent);
            if (nodes[second_vertex].label == nodes_label_cnt) {
                lca = second_vertex;
                break;
            }
            nodes[second_vertex].label = nodes_label_cnt;
        }
    }

    FracInitAugment(lca, this_root);

    first_vertex = elementary_heads[edge];
    second_vertex = elementary_tails[edge];
    while (second_vertex != lca) {
        half_int_clockwise[second_vertex] = nodes[second_vertex].matched_edge;
        MakeEdgeUnmatched(nodes[second_vertex].matched_edge.index >> 1);
        nodes[second_vertex].matched_edge = ArcIndex(-1);
        second_vertex = OtherElementaryEnd(half_int_clockwise[second_vertex]);
        nodes[second_vertex].matched_edge = ArcIndex(-1);

        half_int_clockwise[second_vertex] = nodes[second_vertex].minus_parent;
        second_vertex = OtherElementaryEnd(nodes[second_vertex].minus_parent);
    }
    half_int_clockwise[first_vertex] = ArcIndex(edge * 2);
    while (first_vertex != lca) {
        ArcIndex upwards_arc = nodes[first_vertex].matched_edge;
        nodes[first_vertex].matched_edge = ArcIndex(-1);
        first_vertex = OtherElementaryEnd(upwards_arc);
        nodes[first_vertex].matched_edge = ArcIndex(-1);
        MakeEdgeUnmatched(upwards_arc.index >> 1);
        half_int_clockwise[first_vertex] = ReverseArc(upwards_arc);

        upwards_arc = nodes[first_vertex].minus_parent;
        first_vertex = OtherElementaryEnd(nodes[first_vertex].minus_parent);
        half_int_clockwise[first_vertex] = ReverseArc(upwards_arc);
    }
    nodes[lca].matched_edge = ArcIndex(-1);
}

void MWPMSolver::FracInitCleanTree(int tree_var) {
    if (params.verbose) {
        std::cout << "INIT CLEAN TREE" << std::endl;
    }
    // change the slacks and vars by tree_var

    for (int node : init_tree_nodes) {
        int diff = tree_var;
        if (!nodes[node].plus) {
            diff = -diff;
        }

        node_heap_infos[node].dual_var_quadrupled_amortized_ += diff;
        for (ArcIndex arc : adj_list[node]) {
            edges[arc.index >> 1].slack_quadrupled_amortized_ -= diff;
        }
    }

    if (params.debug) {
        for (int edge = 0; edge < edges.size(); ++edge) {
            if (edges[edge].slack_quadrupled_amortized_ < 0) {
                throw std::runtime_error("FracInitCleanTree: got negative slack");
            }
        }
    }
}

void MWPMSolver::FracInitCleanup() {
    for (int node = 0; node < num_vertices_elementary; ++node) {
        nodes[node].tree = -1;
        nodes[node].minus_parent = ArcIndex(-1);

        if (half_int_clockwise[node].index >= 0) {
            int cur_node = OtherElementaryEnd(half_int_clockwise[node]);
            half_int_clockwise[node] = ArcIndex(-1);
            bool match = true;
            while (cur_node != node) {
                if (match) {
                    MakeEdgeMatched(half_int_clockwise[cur_node].index >> 1);
                }
                match = !match;
                int next_node = OtherElementaryEnd(half_int_clockwise[cur_node]);
                half_int_clockwise[cur_node] = ArcIndex(-1);
                cur_node = next_node;
            }
        }
    }
}

void MWPMSolver::InitializeTrees() {
    // is called after GreedyInit
    // also initializes actionable_edges and primal update record

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
    tree_heap_infos.reserve(roots.size());
    alive_trees.reserve(roots.size());
    tree_nodes.reserve(roots.size());
    for (int i = 0; i < static_cast<int>(roots.size()); ++i) {
        trees.emplace_back(Tree(roots[i]));
        tree_nodes.emplace_back(std::deque<int>{roots[i]});
        alive_trees.emplace_back(i);
        tree_heap_infos.emplace_back(i, 2 * i, 2 * i + 1);

        node_heaps.emplace_back(NodeHeap());
        edge_heaps.emplace_back(EdgeHeap());
        edge_heaps.emplace_back(EdgeHeap());
        edge_heap_alive.push_back(1);
        edge_heap_alive.push_back(1);

        nodes[roots[i]].tree = i;
        nodes[roots[i]].plus = true;
    }

    // initialize actionable_edges and record
    for (int root_index : roots) {
        for (ArcIndex arc : adj_list[root_index]) {
            int edge_to_neighbor = arc.index >> 1;

            AddNodeToRecord(root_index);

            if (edges[edge_to_neighbor].slack_quadrupled_amortized_ == 0) {
                // safe because no dual updates has been made yet
                actionable_edges.push_back(edge_to_neighbor);
            }
        }
    }
}

void MWPMSolver::ComputeMatching() {
    // must be called after the blossoms are dissolved
    if (!matching.empty()) {
        throw std::runtime_error("ComputeMatching: matching is already non-empty");
    }

    matching.reserve(num_vertices_elementary >> 1);

    for (int i(0); i < static_cast<int>(edges.size()); ++i) {
        if (matched[i]) {
            matching.emplace_back(elementary_heads[i], elementary_tails[i]);
        }
    }
}

void MWPMSolver::ComputePrimalObjective() {
    primal_objective = 0;
    for (int i(0); i < static_cast<int>(edges.size()); ++i) {
        if (matched[i]) {
            primal_objective += edge_weights[i];
        }
    }
}

void MWPMSolver::DestroyBlossoms() {
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
        for (int child : blossom_structures[blossom].blossom_children) {
            blossom_parents[child] = -1;
            blossom_ancestors[child] = -1;
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

void MWPMSolver::RestoreFinalEdgeEnds() {
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
            edges_to_restore.push_back(nodes[node].matched_edge.index >> 1);
        }
        if (nodes[node].minus_parent.index >= 0) {
            edges_to_restore.push_back(nodes[node].minus_parent.index >> 1);
        }
    }

    for (int edge : edges_to_restore) {
        int head = elementary_heads[edge];
        int tail = elementary_tails[edge];

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

int MWPMSolver::FindFinalReceptacle(int blossom) const {
    int new_receptacle = ThisElementaryEnd(nodes[blossom].matched_edge);
    while (blossom_parents[new_receptacle] != blossom) {
        new_receptacle = blossom_parents[new_receptacle];
    }

    if (new_receptacle < 0) {
        throw std::runtime_error("DestroyBlossoms: no new receptacle found");
    }

    return new_receptacle;
}

void MWPMSolver::ComputeDualCertificate() {
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
            for (int child : blossom_structures[i].blossom_children) {
                std::get<2>(dual_certificate[node_index_alive[child]]) = index;
            }
            ++index;
        }
    }
}

void MWPMSolver::ComputeDualObjectiveQuadrupled() {
    dual_objective = 0;
    for (int i(0); i < static_cast<int>(nodes.size()); ++i) {
        if (nodes[i].is_alive) {
            dual_objective += DualVariableQuadrupled(i);
        }
    }
}

bool MWPMSolver::MakePrimalUpdates() {
    bool action_taken = false;

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
    while (actionable_edges_head < actionable_edges.size()) {
        int edge = actionable_edges[actionable_edges_head++];
        MakePrimalUpdate(edge);
    }
    actionable_edges.clear();
    actionable_edges_head = 0;

    while (actionable_nodes_head < actionable_nodes.size()) {
        int node = actionable_nodes[actionable_nodes_head++];
        MakePrimalUpdateForNode(node);
    }
    actionable_nodes.clear();
    actionable_nodes_head = 0;

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
                std::cout << "edge " << elementary_heads[i] << " " << elementary_tails[i] <<
                    std::endl;
                std::cout << "edge index: " << i << std::endl;
                std::cout << "head: " << Head(i) << std::endl;
                std::cout << "tail: " << Tail(i) << std::endl;
                std::cout << "slacks before/after: " << slacks[i] << " " << new_slacks[i] << std::endl;
                throw std::runtime_error("slacks do not match");
            }
        }
    }

    // shrink phase
    if (params.verbose) {
        std::cout << "shrinking phase" << std::endl;
    }
    OrganizeBlossomChildren();
    for (std::vector<int> &children : future_blossoms_tmp) {
        // if (children.size() > 1) {
        Shrink(children);
        // }
    }
    future_blossoms_tmp.clear();

    primal_update_record.clear();
    return action_taken;
}

void MWPMSolver::MakePrimalUpdate(int edge) {
    int parent = Head(edge);
    int child = Tail(edge);
    if (parent == child) {
        return;
    }

    ArcIndex arc(edge * 2);
    if (nodes[parent].tree < 0) {
        std::swap(parent, child);
        arc = ReverseArc(arc);
    }

    if (nodes[parent].tree < 0) {
        // edge not adjacent to any trees
        return;
    }

    if (nodes[parent].plus) {
        if (nodes[child].tree < 0) {
            Grow(parent, arc);
            return;
        }

        if (nodes[child].plus) {
            if (nodes[parent].tree == nodes[child].tree) {
                MakeCherryBlossom(edge);
            } else {
                Augment(edge);
            }
        }
    }
}

void MWPMSolver::MakePrimalUpdateForNode(int node) {
    if (!nodes[node].plus || !nodes[node].is_alive) {
        return;
    }

    int delta_slack = 0;
    if (nodes[node].old_tree >= 0) {
        if (nodes[node].old_plus) {
            delta_slack = -trees[nodes[node].old_tree].dual_var_quadrupled;
        } else {
            delta_slack = trees[nodes[node].old_tree].dual_var_quadrupled;
        }
    }

    UpdateNonLoopNeighbors(node);
    for (ArcIndex arc : adj_list[node]) {
        int edge = arc.index >> 1;
        // if (maybe_has_zero_slack[edge]) {
        if (true) {
            // compute the old slack to compare it to 0
            int slack = edges[edge].slack_quadrupled_amortized_;
            slack += delta_slack;

            int neighbor = OtherEnd(arc);
            if (nodes[neighbor].old_tree >= 0) {
                if (nodes[neighbor].old_plus) {
                    slack -= trees[nodes[neighbor].old_tree].dual_var_quadrupled;
                } else {
                    slack += trees[nodes[neighbor].old_tree].dual_var_quadrupled;
                }
            }

            if (slack == 0) {
                // make a primal update
                if (node == neighbor) {
                    throw std::runtime_error("In MakePrimalUpdateForNode: loop");
                }
                if (nodes[neighbor].tree < 0) {
                    Grow(node, arc);
                } else if (nodes[neighbor].plus) {
                    if (nodes[node].tree == nodes[neighbor].tree) {
                        MakeCherryBlossom(edge);
                    } else {
                        Augment(edge);
                        return;
                    }
                }
            } else {
                // maybe_has_zero_slack[edge] = false;
            }
        }
    }
}

void MWPMSolver::Expand(int blossom) {
    ++aux_counter2;

    if (params.verbose) {
        std::cout << "EXPAND " << blossom << std::endl;
    }

    nodes[blossom].is_alive = false;

    // TODO make better
    RestoreEdgeEndsBeforeExpand(blossom);
    ChangeLoopSlacksBeforeExpand(blossom);
    // now ThisEnd/OtherEnd is safe

    int new_receptacle = ThisEnd(nodes[blossom].matched_edge);
    int old_receptacle = Receptacle(new_receptacle);
    int elder_child = ThisEnd(nodes[blossom].minus_parent);

    UpdateInternalStructure(blossom, old_receptacle, new_receptacle, elder_child);

    for (int child : blossom_structures[blossom].blossom_children) {
        AddNodeToRecord(child);
        if (nodes[child].tree >= 0) {
            tree_nodes[nodes[child].tree].push_back(child);
            if (nodes[child].plus) {
                actionable_nodes.push_back(child);
            }
        }
    }
}

void MWPMSolver::RestoreEdgeEndsBeforeExpand(int blossom) {
    // makes head and tail of the adjacent edges to point to the children of the blossom,
    // make children top blossoms

    // TODO make better?

    for (int child : blossom_structures[blossom].blossom_children) {
        if (IsElementary(child)) {
            for (ArcIndex arc : adj_list[child]) {
                if (arc.index % 2 == 1) {
                    edges[arc.index >> 1].tail = child;
                } else {
                    edges[arc.index >> 1].head = child;
                }
            }
            continue;
        }

        traversal_nodes_tmp.clear();
        traversal_nodes_tmp.push_back(child);
        for (std::size_t queue_index = 0; queue_index < traversal_nodes_tmp.size(); ++queue_index) {
            int cur = traversal_nodes_tmp[queue_index];
            for (int potential_leaf : blossom_structures[cur].blossom_children) {
                blossom_ancestors[potential_leaf] = child;
                if (IsElementary(potential_leaf)) {
                    for (ArcIndex arc : adj_list[potential_leaf]) {
                        if (arc.index % 2 == 1) {
                            edges[arc.index >> 1].tail = child;
                        } else {
                            edges[arc.index >> 1].head = child;
                        }
                    }
                } else {
                    traversal_nodes_tmp.push_back(potential_leaf);
                }
            }
        }
    }

    for (int child : blossom_structures[blossom].blossom_children) {
        blossom_parents[child] = -1;
        blossom_ancestors[child] = -1;
    }
}

void MWPMSolver::ChangeLoopSlacksBeforeExpand(int blossom) {
    // TODO make better?

    int old_tree_of_children = nodes[blossom_structures[blossom].blossom_children[0]].old_tree;
    int old_tree_var_of_children = trees[old_tree_of_children].dual_var_quadrupled;

    int old_tree_blossom = nodes[blossom].old_tree;
    int old_tree_var_blossom = 0;
    if (old_tree_blossom >= 0) {
        old_tree_var_blossom = trees[old_tree_blossom].dual_var_quadrupled;
    }

    for (int child : blossom_structures[blossom].blossom_children) {
        UpdateNonLoopNeighbors(child);
    }

    ++nodes_label_cnt;
    for (int child : blossom_structures[blossom].blossom_children) {
        nodes[child].label = nodes_label_cnt;
        node_heap_infos[child].dual_var_quadrupled_amortized_ -= old_tree_var_of_children;
    }

    for (int child : blossom_structures[blossom].blossom_children) {
        for (ArcIndex arc : adj_list[child]) {
            int edge = arc.index >> 1;

            if (nodes[OtherEnd(arc)].label == nodes_label_cnt) {
                edges[edge].slack_quadrupled_amortized_ -= nodes[blossom].tree_var_at_birth;
                edges[edge].slack_quadrupled_amortized_ += old_tree_var_of_children;
            } else {
                if (old_tree_blossom >= 0) {
                    if (nodes[blossom].old_plus) {
                        edges[edge].slack_quadrupled_amortized_ -= old_tree_var_blossom;
                    } else {
                        edges[edge].slack_quadrupled_amortized_ += old_tree_var_blossom;
                    }
                }
                edges[edge].slack_quadrupled_amortized_ += old_tree_var_of_children;
            }

            if (!nodes[child].old_plus || OtherEnd(arc) == child) {
                throw std::runtime_error("In ChangeLoopSlacksBeforeExpand");
            }
        }
    }
}

void MWPMSolver::RotateReceptacle(int blossom, int new_receptacle) {
    // changes matched edges and minus_parents

    if (params.verbose) {
        std::cout << "ROTATE RECEPTACLE " << blossom << " " << new_receptacle << std::endl;
    }

    int old_receptacle = Receptacle(new_receptacle);
    if (old_receptacle == new_receptacle) {
        nodes[new_receptacle].matched_edge = nodes[blossom].matched_edge;
        return;
    }

    EvenPathToReceptacle(new_receptacle);
    OddPathToReceptacle(new_receptacle);

    if (even_path_tmp.back().index == odd_path_tmp.back().index) {
        int i = 0;
        int common_node = old_receptacle;
        while (even_path_tmp[even_path_tmp.size() - 1 - i].index == odd_path_tmp[odd_path_tmp.size() - 1 - i].index) {
            common_node = ThisEnd(even_path_tmp[even_path_tmp.size() - 1 - i]);
            ++i;
        }

        RotateReceptacle(blossom, common_node);
        RotateReceptacle(blossom, new_receptacle);
        return;
    }

    int cur_node = new_receptacle;
    bool match = false;
    for (ArcIndex arc : even_path_tmp) {
        int edge = arc.index >> 1;
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
    for (ArcIndex arc : odd_path_tmp) {
        if (match) {
            cur_node = OtherEnd(arc);
        } else {
            if (!matched[arc.index >> 1]) {
                nodes[cur_node].minus_parent = arc;
                cur_node = OtherEnd(arc);
                nodes[cur_node].minus_parent = ReverseArc(arc);
            }
        }
        match = !match;
    }

    nodes[new_receptacle].minus_parent = ArcIndex(-1);
    nodes[new_receptacle].matched_edge = nodes[blossom].matched_edge;
    nodes[new_receptacle].receptacle_ = new_receptacle;
    nodes[old_receptacle].receptacle_ = new_receptacle;

    if (params.debug) {
        for (int child : blossom_structures[blossom].blossom_children) {
            if (nodes[child].minus_parent.index >= 0) {
                if (matched[nodes[child].minus_parent.index >> 1]) {
                    std::cout << blossom << " " << child << " " << nodes[child].minus_parent.index << std::endl;
                    throw std::runtime_error("In RotateReceptacle");
                }
            }
        }
    }
}

void MWPMSolver::UpdateMatching(int blossom, int new_receptacle) {
    EvenPathToReceptacle(new_receptacle);
    bool match = false;
    for (ArcIndex arc : even_path_tmp) {
        int edge = arc.index >> 1;
        if (match) {
            MakeEdgeMatched(edge);
        } else {
            MakeEdgeUnmatched(edge);
        }
        match = !match;
    }
    nodes[new_receptacle].matched_edge = nodes[blossom].matched_edge;
}

void MWPMSolver::UpdateInternalStructure(int blossom,
                                          int old_receptacle,
                                          int new_receptacle,
                                          int elder_child) {
    // after expand, the remaining part of the blossom is a path in the tree

    RotateReceptacle(blossom, new_receptacle);
    EvenPathToReceptacle(elder_child);

    // the path that stays in the tree
    int cur_node = elder_child;
    nodes[elder_child].minus_parent = nodes[blossom].minus_parent;
    nodes[elder_child].plus = false;
    nodes[elder_child].tree = nodes[blossom].tree;
    for (ArcIndex arc : even_path_tmp) {
        cur_node = OtherEnd(arc);

        if (matched[arc.index >> 1]) {
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
    for (ArcIndex arc : even_path_tmp) {
        nodes[cur_node].label = nodes_label_cnt;
        AddNodeToRecord(cur_node);
        cur_node = OtherEnd(arc);
    }
    nodes[new_receptacle].label = nodes_label_cnt;
    AddNodeToRecord(new_receptacle);

    // clear the part that goes to waste
    for (int child : blossom_structures[blossom].blossom_children) {
        if (nodes[child].label != nodes_label_cnt) {
            nodes[child].plus = false;
            nodes[child].minus_parent = ArcIndex(-1);
            nodes[child].tree = -1;
        }
    }

    // reset receptacles
    for (int child : blossom_structures[blossom].blossom_children) {
        nodes[child].receptacle_ = child;
    }
}

void MWPMSolver::EvenPathToReceptacle(int node) {
    int receptacle = Receptacle(node);
    even_path_tmp.clear();

    while (node != receptacle) {
        int parent = OtherEnd(nodes[node].matched_edge);
        int grandparent = OtherEnd(nodes[parent].minus_parent);

        if (params.verbose) {
            std::cout << "EvenPathToReceptacle, parent/grandparent: " << parent << " " << grandparent << std::endl;
        }

        even_path_tmp.push_back(nodes[node].matched_edge);
        even_path_tmp.push_back(nodes[parent].minus_parent);

        node = grandparent;

        // if (even_path_tmp.size() > 3 * num_vertices_elementary) {
        //     std::cout << "receptacle: " << receptacle << ", node: " << node << std::endl;
        //     std::cout << nodes[receptacle].matched_edge.index << " " << nodes[receptacle].minus_parent.index << " " <<
        //         nodes[receptacle].tree << " " << nodes[receptacle].old_tree << std::endl;
        //     std::cout << nodes[node].is_alive << " " << nodes[receptacle].is_alive << " " << blossom_parents[node] <<
        //         " " << blossom_parents[receptacle] << std::endl;
        //     std::cout << nodes[node].matched_edge.index << " " << nodes[parent].minus_parent.index << " " << bool(matched[
        //             nodes[node].matched_edge.index >> 1]) << " " << bool(matched[nodes[parent].minus_parent.index >> 1]) <<
        //         std::endl;
        //     throw std::runtime_error("EvenPathToReceptacle: infinite loop");
        // }
    }
}

void MWPMSolver::OddPathToReceptacle(int node) {
    int receptacle = Receptacle(node);
    if (receptacle == node) {
        throw std::runtime_error("OddPathToReceptacle: node is receptacle");
    }

    odd_path_tmp.clear();
    odd_path_tmp.push_back(nodes[node].minus_parent);
    node = OtherEnd(nodes[node].minus_parent);
    while (node != receptacle) {
        int parent = OtherEnd(nodes[node].matched_edge);
        int grandparent = OtherEnd(nodes[parent].minus_parent);

        odd_path_tmp.push_back(nodes[node].matched_edge);
        odd_path_tmp.push_back(nodes[parent].minus_parent);

        node = grandparent;
    }
}

void MWPMSolver::ExpandChildBeforeGrow(int blossom) {
    ++aux_counter2;

    nodes[blossom].is_alive = false;

    RestoreEdgeEndsBeforeExpand(blossom);
    ChangeLoopSlacksBeforeExpand(blossom);
    // now OtherEnd, Head, Tail is safe

    int new_receptacle = ThisEnd(nodes[blossom].matched_edge);
    UpdateMatching(blossom, new_receptacle);

    for (int child : blossom_structures[blossom].blossom_children) {
        nodes[child].receptacle_ = child;
        nodes[child].plus = false;
        nodes[child].minus_parent = ArcIndex(-1);
        nodes[child].tree = -1;
        AddNodeToRecord(child);
    }
}

void MWPMSolver::Grow(int parent, ArcIndex arc) {
    ++aux_counter1;

    int child = OtherEnd(arc);
    int tree = nodes[parent].tree;

    // if (!nodes[parent].plus) {
    //     throw std::runtime_error("In Grow: parent is not a plus");
    // }
    // if (nodes[child].tree >= 0) {
    //     throw std::runtime_error("In Grow: child vertex is not free");
    // }
    // if (nodes[child].matched_edge.index < 0) {
    //     std::cout << "parent " << parent << ", child " << child << std::endl;
    //     throw std::runtime_error("Grow: child has no matched_edge");
    // }

    int grandchild = OtherEnd(nodes[child].matched_edge);
    if (params.verbose) {
        std::cout << "GROW " << parent << " " << child << " " << grandchild << std::endl;
    }

    if (!IsElementary(child)) {
        if (DualVariableQuadrupled(child,
                                   nodes[child].old_tree,
                                   nodes[child].old_plus,
                                   -1) == 0) {
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

    tree_nodes[tree].push_back(child);
    tree_nodes[tree].push_back(grandchild);

    actionable_nodes.push_back(grandchild);
}

void MWPMSolver::MakeCherryBlossom(int edge_plus_plus) {
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

std::pair<int, int> MWPMSolver::CherryPathBounds(int first_vertex,
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

void MWPMSolver::UpdateCherryPath(int lower_node, int upper_node) {
    int receptacle = Receptacle(upper_node);

    // TODO only add node to record if it changes sign, have an alternative record for the future shrinks

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
            actionable_nodes.push_back(parent);
        }

        if (grandparent != upper_node) {
            nodes[grandparent].minus_parent = ReverseArc(nodes[parent].minus_parent);
        }

        lower_node = grandparent;
    }
}

void MWPMSolver::Augment(int edge_plus_plus) {
    num_trees_alive -= 2;

    int head = Head(edge_plus_plus);
    int tail = Tail(edge_plus_plus);

    if (params.verbose) {
        std::cout << "AUGMENT " << head << " " << tail << ", elementary ends: " << elementary_heads[edge_plus_plus]
            << " " << elementary_tails[edge_plus_plus] << std::endl;
    }

    PathToRoot(head);
    AugmentPathToRoot();
    PathToRoot(tail);
    AugmentPathToRoot();
    MakeEdgeMatched(edge_plus_plus);

    int first_tree = nodes[head].tree;
    int second_tree = nodes[tail].tree;
    ClearTree(first_tree);
    ClearTree(second_tree);
}

void MWPMSolver::PathToRoot(int node_plus) {
    int root = TopBlossom(roots[nodes[node_plus].tree]);

    path_to_root.clear();
    while (node_plus != root) {
        path_to_root.push_back(nodes[node_plus].matched_edge.index >> 1);
        if (nodes[node_plus].matched_edge.index < 0) {
            throw std::runtime_error("");
        }
        node_plus = OtherEnd(nodes[node_plus].matched_edge);
        path_to_root.push_back(nodes[node_plus].minus_parent.index >> 1);
        // if (nodes[node_plus].minus_parent.index < 0) {
        //     throw std::runtime_error("");
        // }
        node_plus = OtherEnd(nodes[node_plus].minus_parent);

        // if (path_to_root.size() > 2 * num_vertices_elementary) {
        //     throw std::runtime_error("PathToRoot: loop");
        // }
    }
}

void MWPMSolver::AugmentPathToRoot() {
    bool match = false;
    for (int edge : path_to_root) {
        if (match) {
            MakeEdgeMatched(edge);
        } else {
            MakeEdgeUnmatched(edge);
        }
        match = !match;
    }
}

void MWPMSolver::ClearTree(int tree) {
    trees[tree].is_alive = false;

    // TODO have a special routine for the last two trees
    // TODO amortize

    for (int node : tree_nodes[tree]) {
        if (blossom_parents[node] >= 0 || nodes[node].tree != tree || !nodes[node].is_alive) {
            continue;
        }

        // TODO make better

        AddNodeToRecord(node);
        nodes[node].tree = -1;
        nodes[node].minus_parent = ArcIndex(-1);
        nodes[node].receptacle_ = node;
        nodes[node].plus = false;
    }

    edge_heap_alive[tree_heap_infos[tree].plus_empty_edges] = 0;
    edge_heap_alive[tree_heap_infos[tree].plus_plus_internal_edges] = 0;
    for (auto [other_tree, heap] : tree_heap_infos[tree].pq_plus_plus) {
        edge_heap_alive[heap] = 0;
    }
    for (auto [other_tree, heap] : tree_heap_infos[tree].pq_plus_minus) {
        edge_heap_alive[heap] = 0;
    }
}

void MWPMSolver::UpdateQueuesRecordTraversal() {
    UpdateQueuesFirstPass();
    UpdateQueuesSecondPass();
    UpdateQueuesThirdPass();
}

void MWPMSolver::UpdateQueuesFirstPass() {
    for (int node : primal_update_record) {
        int old_tree = nodes[node].old_tree;
        bool old_plus = nodes[node].old_plus;
        bool plus = nodes[node].plus;
        int tree = nodes[node].tree;

        if (!nodes[node].is_alive || (old_plus == plus && old_tree == tree)) {
            continue;
        }

        UpdateNonLoopNeighbors(node);

        // update dual_var_quadrupled_amortized_
        if (old_tree != tree || old_plus != plus) {
            if (old_tree >= 0) {
                if (old_plus) {
                    node_heap_infos[node].dual_var_quadrupled_amortized_ += trees[old_tree].dual_var_quadrupled;
                } else {
                    node_heap_infos[node].dual_var_quadrupled_amortized_ -= trees[old_tree].dual_var_quadrupled;
                }
            }
            if (tree >= 0) {
                if (plus) {
                    node_heap_infos[node].dual_var_quadrupled_amortized_ -= trees[tree].dual_var_quadrupled;
                } else {
                    node_heap_infos[node].dual_var_quadrupled_amortized_ += trees[tree].dual_var_quadrupled;
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
                    AddNodeToQueue(node, tree_heap_infos[tree].minus_blossoms);
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
        if (nodes[node].old_tree >= 0) {
            if (nodes[node].old_plus) {
                nodes[node].slack_diff -= trees[nodes[node].old_tree].dual_var_quadrupled;
            } else {
                nodes[node].slack_diff += trees[nodes[node].old_tree].dual_var_quadrupled;
            }
        }
    }
}

void MWPMSolver::UpdateQueuesSecondPass() {
    for (int node : primal_update_record) {
        if (!nodes[node].is_alive || (nodes[node].old_plus == nodes[node].plus &&
            nodes[node].old_tree == nodes[node].tree)) {
            continue;
        }

        // TODO maybe first update all slacks, then all queues
        if (nodes[node].tree >= 0) {
            if (nodes[node].plus) {
                HandleIncidentPlus(node);
            } else {
                HandleIncidentMinus(node);
            }
        } else {
            HandleIncidentEmpty(node);
        }
    }
}

void MWPMSolver::HandleIncidentEmpty(int node) {
    for (ArcIndex arc : adj_list[node]) {
        int edge = arc.index >> 1;
        if (edges[edge].last_round_updated < current_round) {
            edges[edge].last_round_updated = current_round;

            int other_end = OtherEnd(arc);

            edges[edge].slack_quadrupled_amortized_ += nodes[node].slack_diff;
            edges[edge].slack_quadrupled_amortized_ += nodes[other_end].slack_diff;

            RemoveEdgeFromQueue(edge);
            if (nodes[other_end].tree >= 0 && nodes[other_end].plus) {
                int queue_index = tree_heap_infos[nodes[other_end].tree].plus_empty_edges;
                if (queue_index >= 0) {
                    AddEdgeToThisQueue(edge, queue_index);
                }
            }
        }
    }
}

void MWPMSolver::HandleIncidentPlus(int node) {
    int tree = nodes[node].tree;
    int receptacle_node = Receptacle(node);

    for (ArcIndex arc : adj_list[node]) {
        int edge = arc.index >> 1;
        if (edges[arc.index >> 1].last_round_updated < current_round) {
            edges[arc.index >> 1].last_round_updated = current_round;
            int queue_index = -1;
            int other_end = OtherEnd(arc);
            int other_tree = nodes[other_end].tree;

            edges[edge].slack_quadrupled_amortized_ += nodes[node].slack_diff;
            edges[edge].slack_quadrupled_amortized_ += nodes[other_end].slack_diff;

            RemoveEdgeFromQueue(edge);
            if (tree == other_tree) {
                if (nodes[other_end].plus) {
                    queue_index = tree_heap_infos[tree].plus_plus_internal_edges;
                    int receptacle_other = Receptacle(other_end);
                    if (receptacle_node != receptacle_other) {
                        AddEdgeToThisQueue(edge, queue_index);
                    }
                }
            } else {
                if (other_tree < 0) {
                    // (+, 0)
                    queue_index = tree_heap_infos[tree].plus_empty_edges;
                } else {
                    if (nodes[other_end].plus) {
                        // (+, +)
                        queue_index = TreeTreeQueueIndex(other_tree, &tree_heap_infos[tree].pq_plus_plus);
                        if (queue_index < 0) {
                            edge_heaps.emplace_back(EdgeHeap());
                            edge_heap_alive.push_back(1);
                            queue_index = static_cast<int>(edge_heaps.size()) - 1;
                            tree_heap_infos[tree].pq_plus_plus.emplace_back(other_tree, queue_index);
                            tree_heap_infos[other_tree].pq_plus_plus.emplace_back(tree, queue_index);
                        }
                    } else {
                        // (+, -)
                        queue_index = TreeTreeQueueIndex(other_tree, &tree_heap_infos[tree].pq_plus_minus);
                        if (queue_index < 0) {
                            edge_heaps.emplace_back(EdgeHeap());
                            edge_heap_alive.push_back(1);
                            queue_index = static_cast<int>(edge_heaps.size()) - 1;
                            tree_heap_infos[tree].pq_plus_minus.emplace_back(other_tree, queue_index);
                        }
                    }
                }

                AddEdgeToThisQueue(edge, queue_index);
            }
        }
    }
}

void MWPMSolver::HandleIncidentMinus(int node) {
    int tree = nodes[node].tree;

    for (ArcIndex arc : adj_list[node]) {
        int edge = arc.index >> 1;
        if (edges[arc.index >> 1].last_round_updated < current_round) {
            edges[arc.index >> 1].last_round_updated = current_round;
            int queue_index = -1;
            int other_end = OtherEnd(arc);
            int other_tree = nodes[other_end].tree;

            edges[edge].slack_quadrupled_amortized_ += nodes[node].slack_diff;
            edges[edge].slack_quadrupled_amortized_ += nodes[other_end].slack_diff;

            RemoveEdgeFromQueue(edge);
            if (other_tree >= 0) {
                if (tree != other_tree) {
                    if (nodes[other_end].plus) {
                        // (-, +)
                        queue_index = TreeTreeQueueIndex(tree, &tree_heap_infos[other_tree].pq_plus_minus);
                        if (queue_index < 0) {
                            edge_heaps.emplace_back(EdgeHeap());
                            edge_heap_alive.push_back(1);
                            queue_index = static_cast<int>(edge_heaps.size()) - 1;
                            tree_heap_infos[other_tree].pq_plus_minus.emplace_back(tree, queue_index);
                        }
                        AddEdgeToThisQueue(edge, queue_index);
                    }
                }
            }
        }
    }
}

void MWPMSolver::UpdateQueuesThirdPass() {
    for (int node : primal_update_record) {
        nodes[node].old_plus = nodes[node].plus;
        nodes[node].old_tree = nodes[node].tree;
        nodes[node].is_in_record = false;
        nodes[node].slack_diff = 0;
    }
}

void MWPMSolver::OrganizeBlossomChildren() {
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

    organize_sizes_tmp.assign(nodes_label_cnt - label_zero, 0);
    organize_label_to_index_tmp.assign(nodes_label_cnt - label_zero, -1);

    for (int node : primal_update_record) {
        if (nodes[node].label < label_zero) {
            continue;
        }
        ++organize_sizes_tmp[nodes[node].label - label_zero];
    }

    int index = 0;
    for (int i = 0; i < static_cast<int>(organize_sizes_tmp.size()); ++i) {
        if (organize_sizes_tmp[i] > 1) {
            organize_label_to_index_tmp[i] = index;
            ++index;
        }
    }

    future_blossoms_tmp.resize(index);
    for (std::vector<int> &children : future_blossoms_tmp) {
        children.clear();
    }
    for (int i = 0; i < static_cast<int>(organize_sizes_tmp.size()); ++i) {
        if (organize_label_to_index_tmp[i] >= 0) {
            future_blossoms_tmp[organize_label_to_index_tmp[i]].reserve(organize_sizes_tmp[i]);
        }
    }

    for (int node : primal_update_record) {
        if (nodes[node].label < label_zero) {
            continue;
        }
        int group_index = organize_label_to_index_tmp[nodes[node].label - label_zero];
        if (group_index >= 0) {
            future_blossoms_tmp[group_index].push_back(node);
        }
    }
}

void MWPMSolver::Shrink(std::vector<int> &children) {
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

    // if (params.debug) {
    //     for (int child : children) {
    //         if (child == 1005615) {
    //             std::cout << "higher blossom " << new_index << std::endl;
    //         }
    //     }
    // }

    nodes.emplace_back(Node(new_index));
    blossom_parents.push_back(-1);
    blossom_ancestors.push_back(-1);
    adj_list.emplace_back();
    blossom_structures.emplace_back();
    nodes.back().plus = true;
    nodes.back().old_plus = true;
    blossom_structures.back().blossom_children = std::move(children);
    nodes.back().matched_edge = nodes[receptacle].matched_edge;
    nodes.back().tree = nodes[receptacle].tree;
    nodes.back().old_tree = nodes.back().tree;
    nodes.back().tree_var_at_birth = trees[nodes.back().tree].dual_var_quadrupled;
    node_heap_infos.emplace_back(-1, -1, -1, -1, -trees[nodes.back().tree].dual_var_quadrupled);

    // update blossom_parent of the children, set labels to empty, update variables
    for (int child : blossom_structures.back().blossom_children) {
        if (blossom_parents[child] >= 0) {
            throw std::runtime_error("In Node: some child node already has a parent");
        }
        blossom_parents[child] = new_index;
        blossom_ancestors[child] = new_index;
        node_heap_infos[child].dual_var_quadrupled_amortized_ += trees[nodes.back().tree].dual_var_quadrupled;
    }

    // add the new blossom to the list of vertices in the tree
    tree_nodes[nodes.back().tree].push_back(new_index);
}

std::vector<int> MWPMSolver::NodeVariables() const {
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

std::vector<int> MWPMSolver::EdgeSlacks() {
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
    for (int edge = 0; edge < static_cast<int>(edges.size()); ++edge) {
        if (Head(edge) != Tail(edge)) {
            result.push_back(SlackQuadrupled(edge));
        } else {
            int head = elementary_heads[edge];
            int tail = elementary_tails[edge];

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

void MWPMSolver::ValidateEvenOddPaths() {
    for (int node(0); node < static_cast<int>(nodes.size()); ++node) {
        if (!nodes[node].is_alive || blossom_parents[node] >= 0) {
            continue;
        }

        int receptacle = Receptacle(node);

        EvenPathToReceptacle(node);
        if (even_path_tmp.size() % 2 != 0) {
            throw std::runtime_error("ValidateEvenOddPaths: even path is not even");
        }
        int cur_node = node;
        for (ArcIndex arc : even_path_tmp) {
            cur_node = OtherEnd(arc);
        }
        if (cur_node != receptacle) {
            throw std::runtime_error("ValidateEvenOddPaths: incorrect even path");
        }

        if (node != receptacle) {
            OddPathToReceptacle(node);
            if (odd_path_tmp.size() % 2 != 1) {
                throw std::runtime_error("ValidateEvenOddPaths: odd path is not odd");
            }
            cur_node = node;
            for (ArcIndex arc : odd_path_tmp) {
                cur_node = OtherEnd(arc);
            }
            if (cur_node != receptacle) {
                throw std::runtime_error("ValidateEvenOddPaths: incorrect odd path");
            }
        }
    }
}

void MWPMSolver::ValidateArcs() {
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

bool MWPMSolver::MakeDualUpdates() {
    UpdateAliveTreesList();

    if (num_trees_alive == 0) {
        return true;
    }

    if (params.verbose) {
        std::cout << "DUAL UPDATE" << std::endl;
    }

    std::vector<DualConstraintsNode> dual_constraints = GetDualConstraints();
    DualUpdater dual_updater(std::move(dual_constraints));
    dual_updater.FindDeltas();
    const std::vector<int> &deltas = dual_updater.Deltas();
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

void MWPMSolver::UpdateAliveTreesList() {
    for (int i = 0; i < static_cast<int>(alive_trees.size()); ++i) {
        if (!trees[alive_trees[i]].is_alive) {
            alive_trees[i] = alive_trees.back();
            alive_trees.pop_back();
            --i;
        }
    }
}

std::vector<DualConstraintsNode> MWPMSolver::GetDualConstraints() {
    std::vector<DualConstraintsNode> result(alive_trees.size());

    // update alive indices
    for (int i = 0; i < static_cast<int>(alive_trees.size()); ++i) {
        trees[alive_trees[i]].alive_index = i;
    }

    for (int i = 0; i < static_cast<int>(alive_trees.size()); ++i) {
        int tree = alive_trees[i];
        int tree_var = trees[tree].dual_var_quadrupled;

        // get upper_bound
        int edge = GetMinEdgeHeap(tree_heap_infos[tree].plus_empty_edges);
        if (edge >= 0) {
            result[i].upper_bound = edges[edge].slack_quadrupled_amortized_ - tree_var;
        }

        edge = MinPlusPlusInternalEdge(tree_heap_infos[tree].plus_plus_internal_edges);
        if (edge >= 0) {
            int plus_plus_internal_halved = edges[edge].slack_quadrupled_amortized_ / 2 - tree_var;
            if (plus_plus_internal_halved < result[i].upper_bound) {
                result[i].upper_bound = plus_plus_internal_halved;
            }
        }

        int node = GetMinNodeHeap(tree_heap_infos[tree].minus_blossoms);
        if (node >= 0) {
            int blossom_var = node_heap_infos[node].dual_var_quadrupled_amortized_ - tree_var;
            if (blossom_var < result[i].upper_bound) {
                result[i].upper_bound = blossom_var;
            }
        }

        // get plus_plus_constraints and plus_minus_constraints

        result[i].plus_plus_neighbors.reserve(tree_heap_infos[tree].pq_plus_plus.size());
        result[i].plus_plus_constraints.reserve(tree_heap_infos[tree].pq_plus_plus.size());
        for (int j = 0; j < static_cast<int>(tree_heap_infos[tree].pq_plus_plus.size()); ++j) {
            auto [tree_neighbor, queue_index] = tree_heap_infos[tree].pq_plus_plus[j];
            if (trees[tree_neighbor].is_alive) {
                edge = GetMinEdgeHeap(queue_index);
                if (edge >= 0) {
                    result[i].plus_plus_neighbors.emplace_back(trees[tree_neighbor].alive_index);
                    result[i].plus_plus_constraints.emplace_back(edges[edge].slack_quadrupled_amortized_ - tree_var -
                        trees[tree_neighbor].dual_var_quadrupled);
                }
            } else {
                tree_heap_infos[tree].pq_plus_plus[j] = tree_heap_infos[tree].pq_plus_plus.back();
                tree_heap_infos[tree].pq_plus_plus.pop_back();
                --j;
            }
        }

        result[i].plus_minus_neighbors.reserve(tree_heap_infos[tree].pq_plus_minus.size());
        result[i].plus_minus_constraints.reserve(tree_heap_infos[tree].pq_plus_minus.size());
        for (int j = 0; j < static_cast<int>(tree_heap_infos[tree].pq_plus_minus.size()); ++j) {
            auto [tree_neighbor, queue_index] = tree_heap_infos[tree].pq_plus_minus[j];
            if (trees[tree_neighbor].is_alive) {
                edge = GetMinEdgeHeap(queue_index);
                if (edge >= 0) {
                    result[i].plus_minus_neighbors.emplace_back(trees[tree_neighbor].alive_index);
                    result[i].plus_minus_constraints.emplace_back(edges[edge].slack_quadrupled_amortized_ - tree_var +
                        trees[tree_neighbor].dual_var_quadrupled);
                }
            } else {
                tree_heap_infos[tree].pq_plus_minus[j] = tree_heap_infos[tree].pq_plus_minus.back();
                tree_heap_infos[tree].pq_plus_minus.pop_back();
                --j;
            }
        }
    }

    return result;
}

void MWPMSolver::InitNextRoundActionable() {
    for (int tree : alive_trees) {
        int tree_var = trees[tree].dual_var_quadrupled;
        for (auto [other_tree, queue_index] : tree_heap_infos[tree].pq_plus_plus) {
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
    }

    for (int tree : alive_trees) {
        int tree_var = trees[tree].dual_var_quadrupled;
        int queue_index = tree_heap_infos[tree].plus_empty_edges;
        int top_edge = GetMinEdgeHeap(queue_index);
        if (top_edge >= 0) {
            if (edges[top_edge].slack_quadrupled_amortized_ - tree_var == 0) {
                AddZeroSlackEdgesFromQueue(queue_index, true);
            }
        }
    }

    for (int tree : alive_trees) {
        int tree_var = trees[tree].dual_var_quadrupled;
        CleanLoopsFromQueueTop(tree);
        int queue_index = tree_heap_infos[tree].plus_plus_internal_edges;
        int top_edge = GetMinEdgeHeap(queue_index);
        if (top_edge >= 0) {
            if (edges[top_edge].slack_quadrupled_amortized_ - 2 * tree_var == 0) {
                AddZeroSlackEdgesFromQueue(queue_index, true);
            }
        }
    }
}

void MWPMSolver::AddZeroSlackEdgesFromQueue(int queue_index, bool add_to_actionable) {
    int top_edge = GetMinEdgeHeap(queue_index);

    int min_key = edges[top_edge].slack_quadrupled_amortized_;

    edges_stack_tmp.clear();
    edges_stack_tmp.push_back(top_edge);

    while (!edges_stack_tmp.empty()) {
        int edge = edges_stack_tmp.back();
        edges_stack_tmp.pop_back();

        if (add_to_actionable) {
            actionable_edges.push_back(edge);
        }

        for (int child = edges[edge].heap_child; child >= 0; child = edges[child].heap_next) {
            if (edges[child].slack_quadrupled_amortized_ == min_key) {
                edges_stack_tmp.push_back(child);
            }
        }
    }
}

void MWPMSolver::CleanLoopsFromQueueTop(int tree) {
    int queue_idx = tree_heap_infos[tree].plus_plus_internal_edges;

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

int MWPMSolver::DualVariableQuadrupled(int node) const {
    return DualVariableQuadrupled(node,
                                  nodes[node].tree,
                                  nodes[node].plus,
                                  blossom_parents[node]);
}

int MWPMSolver::DualVariableQuadrupled(int node, int tree, bool plus, int blossom_parent) const {
    if (tree < 0 || blossom_parent >= 0) {
        return node_heap_infos[node].dual_var_quadrupled_amortized_;
    }
    if (plus) {
        return node_heap_infos[node].dual_var_quadrupled_amortized_ + trees[tree].dual_var_quadrupled;
    }
    return node_heap_infos[node].dual_var_quadrupled_amortized_ - trees[tree].dual_var_quadrupled;
}

void MWPMSolver::UpdateNonLoopNeighbors(int node) {
    if (!adj_list[node].empty()) {
        return;
    }

    traversal_nodes_tmp.clear();
    traversal_lists_tmp.clear();
    int total_length = 0;
    traversal_nodes_tmp.push_back(node);
    ++nodes_label_cnt;
    nodes[node].label = nodes_label_cnt;
    for (std::size_t queue_index = 0; queue_index < traversal_nodes_tmp.size(); ++queue_index) {
        int cur = traversal_nodes_tmp[queue_index];
        for (int child : blossom_structures[cur].blossom_children) {
            nodes[child].label = nodes_label_cnt;

            if (!adj_list[child].empty()) {
                traversal_lists_tmp.push_back(child);
                total_length += adj_list[child].size();
            } else {
                traversal_nodes_tmp.push_back(child);
            }
        }
    }

    adj_list[node].reserve(total_length);
    for (int descendant : traversal_lists_tmp) {
        for (ArcIndex arc : adj_list[descendant]) {
            int edge = arc.index >> 1;
            if (arc.index % 2 == 1) {
                while (blossom_ancestors[edges[edge].head] >= 0 && nodes[edges[edge].head].label != nodes_label_cnt) {
                    edges[edge].head = blossom_ancestors[edges[edge].head];
                }
                if (nodes[edges[edge].head].label != nodes_label_cnt) {
                    adj_list[node].push_back(arc);
                    edges[edge].tail = node;
                }
            } else {
                while (blossom_ancestors[edges[edge].tail] >= 0 && nodes[edges[edge].tail].label != nodes_label_cnt) {
                    edges[edge].tail = blossom_ancestors[edges[edge].tail];
                }
                if (nodes[edges[edge].tail].label != nodes_label_cnt) {
                    adj_list[node].push_back(arc);
                    edges[edge].head = node;
                }
            }
        }
    }
}

int MWPMSolver::SlackQuadrupled(int edge) {
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

int MWPMSolver::PlusPlusLCA(int first_vertex, int second_vertex) {
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

void MWPMSolver::MakeEdgeMatched(int edge) {
    if (Head(edge) == Tail(edge)) {
        throw std::runtime_error("In MakeEdgeMatched: trying to match a loop");
    }
    nodes[Head(edge)].matched_edge = ArcIndex(edge * 2);
    nodes[Tail(edge)].matched_edge = ArcIndex(edge * 2 + 1);
    matched[edge] = true;
}

void MWPMSolver::MakeEdgeUnmatched(int edge) {
    matched[edge] = false;
}

void MWPMSolver::AddEdgeToThisQueue(int edge, int queue_index) {
    InsertEdgeHeap(edge, queue_index);
    edges[edge].queue_index = queue_index;
}

void MWPMSolver::RemoveEdgeFromQueue(int edge) {
    if (edges[edge].queue_index >= 0) {
        RemoveEdgeHeap(edge);
        edges[edge].queue_index = -1;
    }
}

void MWPMSolver::AddNodeToQueue(int node, int queue_index) {
    if (node_heap_infos[node].queue_index == queue_index) {
        return;
    }

    if (node_heap_infos[node].queue_index >= 0) {
        RemoveNodeHeap(node);
    }
    InsertNodeHeap(node, queue_index);
    node_heap_infos[node].queue_index = queue_index;
}

void MWPMSolver::RemoveNodeFromQueue(int node) {
    if (node_heap_infos[node].queue_index >= 0) {
        RemoveNodeHeap(node);
        node_heap_infos[node].queue_index = -1;
    }
}

int MWPMSolver::TreeTreeQueueIndex(int other_tree,
                                    boost::container::small_vector<std::pair<int, int>, 4> *tree_neighbors) const {
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

int MWPMSolver::MinPlusPlusInternalEdge(int queue_index) {
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

int MWPMSolver::PopExpandableBlossom(int tree) {
    int node = GetMinNodeHeap(tree_heap_infos[tree].minus_blossoms);
    if (node >= 0) {
        if (node_heap_infos[node].dual_var_quadrupled_amortized_ - trees[tree].dual_var_quadrupled == 0) {
            RemoveMinNodeHeap(tree_heap_infos[tree].minus_blossoms);
            return node;
        }
    }
    return -1;
}

void MWPMSolver::AddNodeToRecord(int node) {
    if (nodes[node].is_in_record) {
        return;
    }
    primal_update_record.push_back(node);
    nodes[node].is_in_record = true;
}

int MWPMSolver::InitNumVertices(const std::vector<std::tuple<int, int, int> > &edge_list_) {
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

int MWPMSolver::GetMinEdgeHeap(int heap_index) const {
    return edge_heaps[heap_index].root;
}

void MWPMSolver::InsertEdgeHeap(int edge, int heap_index) {
    // ++aux_counter4;

    edges[edge].heap_child = -1;
    edges[edge].heap_next = -1;
    edges[edge].heap_prev = -1;

    edge_heaps[heap_index].root = MeldEdgeHeap(edge_heaps[heap_index].root, edge);
}

void MWPMSolver::RemoveMinEdgeHeap(int heap_index) {
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

void MWPMSolver::RemoveEdgeHeap(int edge) {
    int heap_index = edges[edge].queue_index;
    if (edge_heap_alive[heap_index]) {
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
}

int MWPMSolver::MeldEdgeHeap(int edge_a, int edge_b) {
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

int MWPMSolver::TwoPassMergeEdgeHeap(int edge_first) {
    if (edge_first < 0) {
        return -1;
    }

    int stack = -1;

    // First pass: pair siblings and push meld results onto stack
    while (edge_first >= 0) {
        int a = edge_first;
        int b = edges[edge_first].heap_next;

        edge_first = (b >= 0) ? edges[b].heap_next : -1;
        edges[a].heap_next = -1;

        int merged;
        if (b >= 0) {
            edges[b].heap_next = -1;
            merged = MeldEdgeHeap(a, b);
        } else {
            merged = a;
        }

        // push onto stack (reverse order)
        edges[merged].heap_next = stack;
        stack = merged;
    }

    // Second pass: meld stack from left to right
    int result = stack;
    stack = edges[stack].heap_next;
    edges[result].heap_next = -1;

    while (stack >= 0) {
        int next = edges[stack].heap_next;
        edges[stack].heap_next = -1;
        result = MeldEdgeHeap(stack, result);
        stack = next;
    }

    return result;
}

void MWPMSolver::CutEdgeHeap(int edge) {
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

int MWPMSolver::GetMinNodeHeap(int heap_index) const {
    return node_heaps[heap_index].root;
}

void MWPMSolver::InsertNodeHeap(int node, int heap_index) {
    node_heap_infos[node].heap_child = -1;
    node_heap_infos[node].heap_next = -1;
    node_heap_infos[node].heap_prev = -1;

    node_heaps[heap_index].root = MeldNodeHeap(node_heaps[heap_index].root, node);
}

void MWPMSolver::RemoveMinNodeHeap(int heap_index) {
    if (node_heaps[heap_index].root < 0) {
        return;
    }

    int old_root = node_heaps[heap_index].root;
    int children = node_heap_infos[node_heaps[heap_index].root].heap_child;

    if (children >= 0) {
        node_heap_infos[children].heap_prev = -1; // new root candidates
    }

    node_heaps[heap_index].root = TwoPassMergeNodeHeap(children);
    node_heap_infos[old_root].heap_child = -1;
}

void MWPMSolver::RemoveNodeHeap(int node) {
    int heap_index = node_heap_infos[node].queue_index;
    if (node == node_heaps[heap_index].root) {
        RemoveMinNodeHeap(heap_index);
        return;
    }

    CutNodeHeap(node);

    int subtree = TwoPassMergeNodeHeap(node_heap_infos[node].heap_child);
    if (subtree >= 0) {
        node_heap_infos[subtree].heap_prev = -1;
    }

    node_heaps[heap_index].root = MeldNodeHeap(node_heaps[heap_index].root, subtree);

    node_heap_infos[node].heap_child = -1;
}

int MWPMSolver::MeldNodeHeap(int node_a, int node_b) {
    if (node_a < 0) {
        return node_b;
    }
    if (node_b < 0) {
        return node_a;
    }

    if (node_heap_infos[node_b].dual_var_quadrupled_amortized_ < node_heap_infos[node_a].
        dual_var_quadrupled_amortized_) {
        std::swap(node_a, node_b);
    }

    node_heap_infos[node_b].heap_prev = node_a;
    node_heap_infos[node_b].heap_next = node_heap_infos[node_a].heap_child;
    if (node_heap_infos[node_a].heap_child >= 0) {
        node_heap_infos[node_heap_infos[node_a].heap_child].heap_prev = node_b;
    }

    node_heap_infos[node_a].heap_child = node_b;

    return node_a;
}

int MWPMSolver::TwoPassMergeNodeHeap(int node_first) {
    if (node_first < 0) {
        return -1;
    }

    int stack = -1;

    while (node_first >= 0) {
        int a = node_first;
        int b = node_heap_infos[node_first].heap_next;

        node_first = (b >= 0) ? node_heap_infos[b].heap_next : -1;
        node_heap_infos[a].heap_next = -1;

        int merged;
        if (b >= 0) {
            node_heap_infos[b].heap_next = -1;
            merged = MeldNodeHeap(a, b);
        } else {
            merged = a;
        }

        node_heap_infos[merged].heap_next = stack;
        stack = merged;
    }

    int result = stack;
    stack = node_heap_infos[stack].heap_next;
    node_heap_infos[result].heap_next = -1;

    while (stack >= 0) {
        int next = node_heap_infos[stack].heap_next;
        node_heap_infos[stack].heap_next = -1;
        result = MeldNodeHeap(stack, result);
        stack = next;
    }

    return result;
}

void MWPMSolver::CutNodeHeap(int node) {
    if (node_heap_infos[node].heap_prev < 0) {
        return; // already root
    }

    if (node_heap_infos[node_heap_infos[node].heap_prev].heap_child == node) {
        // node is leftmost child
        int parent = node_heap_infos[node].heap_prev;
        node_heap_infos[parent].heap_child = node_heap_infos[node].heap_next;
        if (node_heap_infos[node].heap_next >= 0) {
            node_heap_infos[node_heap_infos[node].heap_next].heap_prev = parent;
        }
    } else {
        // node has left sibling
        int left_sibling = node_heap_infos[node].heap_prev;
        node_heap_infos[left_sibling].heap_next = node_heap_infos[node].heap_next;
        if (node_heap_infos[node].heap_next >= 0) {
            node_heap_infos[node_heap_infos[node].heap_next].heap_prev = left_sibling;
        }
    }

    node_heap_infos[node].heap_prev = -1;
    node_heap_infos[node].heap_next = -1;
}
