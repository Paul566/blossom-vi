#include "SolverWeighted.h"
#include "Tree.h"

#include <iostream>
#include <stdint.h>
#include <unordered_set>

SolverWeighted::SolverWeighted(const std::vector<std::tuple<int, int, int> > &edge_list_,
                               const SolverParams &params_) : primal_objective(INT64_MAX),
                                                              dual_objective(0), params{params_} {
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

    // PrintElementaryAdjList();

    for (Node &root : elementary_nodes_list) {
        if (root.matched_edge) {
            continue;
        }

        // trees.emplace_back(&root, &blossoms);
        // Tree & tree = trees.back();
        Tree tree(&root, &blossoms, &iter_to_self);
        bool is_augmented = false;
        while (!is_augmented) {
            while (tree.MakePrimalUpdate(&is_augmented)) {
            }
            if (!is_augmented) {
                tree.MakeDualUpdate();
            }
        }
    }

    int64_t quadrupled_dual = DualObjectiveQuadrupled();
    if (quadrupled_dual % 4 != 0) {
        throw std::runtime_error("Dual objective not integer");
    }
    dual_objective = quadrupled_dual / 4;

    if (params.verbose) {
        std::cout << "Dual objective:\t\t" << dual_objective << std::endl;
    }

    DestroyBlossoms();

    primal_objective = PrimalObjective();

    if (params.verbose) {
        std::cout << "Primal objective:\t" << primal_objective << std::endl;
    }
}

void SolverWeighted::PrintElementaryAdjList() const {
    std::cout << "Adjacency list (to, weight, slack, matched):" << std::endl;
    for (const Node &vertex : elementary_nodes_list) {
        std::cout << vertex.index << ": y_v = " << vertex.DualVariableQuadrupled() / 4. << ", ";
        for (const EdgeWeighted *edge : vertex.neighbors) {
            std::cout << "(" << edge->OtherElementaryEnd(vertex).index << ", " << edge->weight << ", " << edge->
                slack_quadrupled / 4. << ", " << edge->matched << ") ";
        }
        std::cout << std::endl;
    }
}

void SolverWeighted::GreedyInit() {
    // first, make all the slacks non-negative
    for (Node &vertex : elementary_nodes_list) {
        if (vertex.neighbors.empty()) {
            continue;
        }

        int min_weight = vertex.neighbors.front()->weight;
        for (const EdgeWeighted *edge : vertex.neighbors) {
            if (edge->weight < min_weight) {
                min_weight = edge->weight;
            }
        }

        vertex.IncreaseDualVariableQuadrupled(min_weight * 2);
    }

    for (Node &vertex : elementary_nodes_list) {
        if ((vertex.matched_edge) || (vertex.neighbors.empty())) {
            continue;
        }

        EdgeWeighted *smallest_slack_edge = vertex.neighbors.front();
        for (EdgeWeighted *edge : vertex.neighbors) {
            if (edge->slack_quadrupled < smallest_slack_edge->slack_quadrupled) {
                smallest_slack_edge = edge;
            }
        }

        const int smallest_slack_quadrupled = smallest_slack_edge->slack_quadrupled;
        vertex.IncreaseDualVariableQuadrupled(smallest_slack_quadrupled);

        if (!smallest_slack_edge->OtherEnd(vertex).matched_edge) {
            // if the other vertex is also unmatched, match the edge
            smallest_slack_edge->matched = true;
            vertex.matched_edge = smallest_slack_edge;
            smallest_slack_edge->OtherEnd(vertex).matched_edge = smallest_slack_edge;
        }
    }
}

std::vector<int> SolverWeighted::RootIndices() const {
    std::vector<int> roots;

    for (const Node &node : elementary_nodes_list) {
        if (!node.matched_edge) {
            roots.push_back(node.index);
        }
    }

    return roots;
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

const std::vector<std::tuple<int, int, std::vector<int>>>& SolverWeighted::DualCertificate() const {
    if (!params.compute_dual_certificate) {
        throw std::runtime_error("In SolverWeighted::DualCertificate: params.compute_dual_certificate is false");
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

void SolverWeighted::DestroyBlossoms() {
    if (Matching().size() * 2 != elementary_nodes_list.size()) {
        throw std::runtime_error("In SolverWeighted: must only destroy blossoms after a perfect matching is found");
    }

    // TODO maybe write guarantees that we are really going top down
    for (auto it = blossoms.rbegin(); it != blossoms.rend(); ++it) {
        if (it->blossom_parent) {
            throw std::runtime_error("In SolverWeighted::RotateMatchingInsideBlossoms: not going top down");
        }

        it->Dissolve();
    }
}
