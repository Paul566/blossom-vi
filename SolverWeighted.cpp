#include "SolverWeighted.h"
#include "Tree.h"

#include <iostream>
#include <cstdint>
#include <queue>
#include <unordered_set>

SolverWeighted::SolverWeighted(const std::vector<std::tuple<int, int, int> > &edge_list_,
                               const SolverParams &params_) : primal_objective(INT64_MAX),
                                                              dual_objective(0), params{params_} {
    // TODO try permuting nodes and edges

    // measure n
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

    elementary_iters.reserve(n);
    for (int i = 0; i < n; ++i) {
        elementary_nodes_list.emplace_back(i);
        elementary_iters.push_back(std::prev(elementary_nodes_list.end()));
    }

    for (const std::tuple<int, int, int> &edge : edge_list_) {
        edges.emplace_back(*elementary_iters[std::get<0>(edge)],
                           *elementary_iters[std::get<1>(edge)],
                           std::get<2>(edge));
        elementary_iters[std::get<0>(edge)]->neighbors.push_back(&edges.back());
        elementary_iters[std::get<1>(edge)]->neighbors.push_back(&edges.back());
    }
}

void SolverWeighted::FindMinPerfectMatching() {
    GreedyInit();

    // initialize trees
    for (Node &root : elementary_nodes_list) {
        if (root.IsMatched()) {
            continue;
        }
        trees.emplace_back(&root, &blossoms, &iter_to_blossom, static_cast<int>(elementary_nodes_list.size()), params.verbose);
        iter_to_tree[&trees.back()] = std::prev(trees.end());
    }

    while (!trees.empty()) {
        MakePrimalUpdates();
        MakeDualUpdates();
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

    DestroyBlossoms();

    primal_objective = PrimalObjective();

    if (params.verbose) {
        std::cout << "Primal objective:\t" << primal_objective << std::endl;
    }
}

void SolverWeighted::PrintGraph() const {
    std::cout << "Adjacency list (to, weight, slack, matched):" << std::endl;
    for (const Node &vertex : elementary_nodes_list) {
        vertex.PrintNode();
    }
    for (const Node &vertex : blossoms) {
        vertex.PrintNode();
    }
}

void SolverWeighted::GreedyInit() {
    // first, make all the slacks non-negative
    for (Node &vertex : elementary_nodes_list) {
        vertex.MakeSlackNonnegativeInInit();
    }

    for (Node &vertex : elementary_nodes_list) {
        vertex.InitVarGreedily();
    }
}

int SolverWeighted::OptimalSingleDelta() {
    // returns the dual variable increment in the single delta approach

    int delta = INT32_MAX;

    for (Tree &tree : trees) {
        int plus_empty = tree.PlusEmptySlack();
        int plus_plus_external = tree.PlusPlusExternalSlack();
        int plus_plus_internal = tree.PlusPlusInternalSlack();
        int blossom_var = tree.MinMinusBlossomVariable();

        if (plus_plus_internal < INT32_MAX && plus_plus_internal % 2 != 0) {
            throw std::runtime_error("PlusPlusInternalSlack is not divisible by 2");
        }
        if (plus_plus_external < INT32_MAX && plus_plus_external % 2 != 0) {
            throw std::runtime_error("PlusPlusExternalSlack is not divisible by 2");
        }

        if (plus_empty < delta) {
            delta = plus_empty;
        }
        if (plus_plus_external / 2 < delta) {
            delta = plus_plus_external / 2;
        }
        if (plus_plus_internal / 2 < delta) {
            delta = plus_plus_internal / 2;
        }
        if (blossom_var < delta) {
            delta = blossom_var;
        }
    }

    return delta;
}

std::vector<int> SolverWeighted::VariableDeltas() {
    DualConstraints dual_constraints = GetDualConstraints();
    std::vector<std::vector<int>> connected_components = ConnectedComponentsTreeTree(dual_constraints);
    std::vector<int> deltas(trees.size(), 0);

    std::vector<int> component_index(trees.size(), 0);
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

std::vector<std::vector<int>> SolverWeighted::ConnectedComponentsTreeTree(const DualConstraints & dual_constraints) {
    // returns a vector of index_of_connected_component

    int n = static_cast<int>(trees.size());
    std::vector<bool> visited(trees.size(), false);
    std::vector<std::vector<int>> components;

    // Build undirected adjacency list for weak connectivity
    std::vector<std::vector<int>> adj_list(n, std::vector<int>());
    for (int u = 0; u < n; ++u) {
        for (auto &[v, slack] : dual_constraints.plus_minus_constraints[u]) {
            if (slack == 0) {
                adj_list[u].emplace_back(v);
                adj_list[v].emplace_back(u);
            }
        }
    }

    for (int start = 0; start < n; ++start) {
        if (!visited[start]) {
            std::vector<int> comp;
            std::queue<int> q;

            visited[start] = true;
            q.push(start);

            while (!q.empty()) {
                int u = q.front(); q.pop();
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

DualConstraints SolverWeighted::GetDualConstraints() {
    std::unordered_map<Tree*, int> tree_indices;    // TODO avoid this unordered map maybe
    int i = 0;
    for (Tree &tree : trees) {
        tree_indices[&tree] = i;
        ++i;
    }

    // get upper_bound
    std::vector<int> upper_bound(trees.size(), INT32_MAX);
    i = 0;
    for (Tree & tree : trees) {
        int plus_empty = tree.PlusEmptySlack();
        int plus_plus_internal = tree.PlusPlusInternalSlack();
        int blossom_var = tree.MinMinusBlossomVariable();

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
    }

    // get plus_plus_constraints and plus_minus_constraints
    std::vector plus_plus_constraints(trees.size(), std::vector<std::pair<int, int> >());
    std::vector plus_minus_constraints(trees.size(), std::vector<std::pair<int, int> >());
    i = 0;
    for (Tree &tree : trees) {
        std::vector<std::pair<Tree *, int>> plus_plus_slacks = tree.PlusPlusExternalSlacks();
        std::vector<std::pair<Tree *, int>> plus_minus_slacks = tree.PlusMinusExternalSlacks();
        for (auto [other_tree, slack] : plus_plus_slacks) {
            plus_plus_constraints[i].emplace_back(tree_indices[other_tree], slack);
        }
        for (auto [other_tree, slack] : plus_minus_slacks) {
            plus_minus_constraints[i].emplace_back(tree_indices[other_tree], slack);
        }
        ++i;
    }

    return {upper_bound, plus_plus_constraints, plus_minus_constraints};
    // TODO avoid copying here
}

std::vector<std::pair<int, int> > SolverWeighted::Matching() const {
    std::vector<std::pair<int, int> > matching;
    matching.reserve(elementary_nodes_list.size() / 2);

    for (const EdgeWeighted &edge : edges) {
        if (edge.matched) {
            auto [head, tail] = edge.ElementaryEndpoints();
            matching.emplace_back(head.index, tail.index);
        }
    }

    return matching;
}

const std::vector<std::tuple<int, int, int> > &SolverWeighted::DualCertificate() const {
    if (!params.compute_dual_certificate) {
        throw std::runtime_error(
            "In SolverWeighted::DualCertificate: params.compute_dual_certificate is false, so the dual certificate was not computed");
    }
    return dual_certificate;
}

int64_t SolverWeighted::DualObjectiveQuadrupled() const {
    int64_t objective = 0;

    for (const Node &vertex : elementary_nodes_list) {
        objective += vertex.DualVariableQuadrupled();
    }
    for (const Node &vertex : blossoms) {
        objective += vertex.DualVariableQuadrupled();
    }

    return objective;
}

int64_t SolverWeighted::PrimalObjective() const {
    int objective = 0;

    for (const EdgeWeighted &edge : edges) {
        if (edge.matched) {
            objective += edge.weight;
        }
    }

    return objective;
}

void SolverWeighted::ComputeDualCertificate() {
    // can be called only before we destroy the blossoms in the end

    if (!dual_certificate.empty()) {
        throw std::runtime_error("Dual certificate is already non-empty");
    }

    std::unordered_map<Node *, int> vtx_to_index;

    dual_certificate.reserve(elementary_nodes_list.size() + blossoms.size());
    int index = 0;
    for (Node &vertex : elementary_nodes_list) {
        dual_certificate.emplace_back(index, vertex.DualVariableQuadrupled(), -1);
        vtx_to_index[&vertex] = index;
        ++index;
    }
    for (Node &blossom : blossoms) {
        dual_certificate.emplace_back(index, blossom.DualVariableQuadrupled(), -1);
        vtx_to_index[&blossom] = index;
        for (Node *child : blossom.BlossomChildren()) {
            std::get<2>(dual_certificate[vtx_to_index[child]]) = index;
        }
        ++index;
    }
}

void SolverWeighted::MakeDualUpdates() {
    // makes the dual updates for the whole collection of the trees

    if (!params.multiple_trees) {
        int delta = OptimalSingleDelta();
        for (Tree &tree : trees) {
            tree.dual_var_quadrupled += delta;
        }

        if (params.verbose) {
            std::cout << "single delta dual progress: " << delta << " number of trees: " << trees.size() << std::endl;
        }
    } else {
        std::vector<int> deltas = VariableDeltas();
        int i = 0;
        for (Tree &tree : trees) {
            tree.dual_var_quadrupled += deltas[i];
            ++i;
        }
    }

    if (params.verbose) {
        std::cout << "after dual update:" << std::endl;
        PrintGraph();
    }
}

void SolverWeighted::MakePrimalUpdates() {
    // makes primal updates for the whole collection of the trees

    // TODO check if the sizes of the trees become imbalanced

    if (params.verbose) {
        std::cout << "making primal updates" << std::endl;
        PrintGraph();
        std::cout << "trees before:" << std::endl;
        for (Tree &tree : trees) {
            tree.PrintTree();
            std::cout << std::endl;
        }
    }

    std::list<Tree*> active_trees;
    std::unordered_map<Tree *, std::list<Tree*>::iterator> active_tree_pointer;
    for (Tree &tree : trees) {
        active_trees.push_back(&tree);
        active_tree_pointer[&tree] = std::prev(active_trees.end());
    }

    auto it = active_trees.begin();
    while (!active_trees.empty()) {
        auto next_it = std::next(it);

        bool success = true;
        Tree *augmented_tree = (*it)->MakePrimalUpdate(&success);
        if (augmented_tree) {
            if (next_it != active_trees.end()) {
                if (*next_it == augmented_tree) {
                    next_it = std::next(next_it);
                }
            }

            trees.erase(iter_to_tree[*it]);
            active_trees.erase(it);
            trees.erase(iter_to_tree[augmented_tree]);
            if (active_tree_pointer.contains(augmented_tree)) {
                active_trees.erase(active_tree_pointer[augmented_tree]);
            }
        }

        if (!success) {
            active_tree_pointer.erase(*it);
            active_trees.erase(it);
        }

        it = next_it;
        if (it == active_trees.end()) {
            it = active_trees.begin();
        }
    }

    if (params.verbose) {
        std::cout << "state after:" << std::endl;
        PrintGraph();
        std::cout << "trees after:" << std::endl;
        for (Tree &tree : trees) {
            tree.PrintTree();
            std::cout << std::endl;
        }
    }
}

void SolverWeighted::DestroyBlossoms() {
    if (Matching().size() * 2 != elementary_nodes_list.size()) {
        throw std::runtime_error("In SolverWeighted: must only destroy blossoms after a perfect matching is found");
    }

    // TODO maybe assert that we are really going top down
    for (auto it = blossoms.rbegin(); it != blossoms.rend(); ++it) {
        if (!it->IsTopBlossom()) {
            throw std::runtime_error("In SolverWeighted::RotateMatchingInsideBlossoms: not going top down");
        }

        it->Dissolve();
    }
}
