#include "SolverWeighted.h"
#include "Tree.h"

#include <iostream>


SolverWeighted::SolverWeighted(const std::vector<std::tuple<int, int, int>> &edge_list_) {
    int n = 0;
    for (const auto &edge : edge_list_) {
        if (std::get<0>(edge) > n) {
            n = std::get<0>(edge);
        }
        if (std::get<1>(edge) > n) {
            n = std::get<1>(edge);
        }
    }
    ++n;

    elementary_nodes = std::vector<std::shared_ptr<Node>>();
    elementary_nodes.reserve(n);
    for (int index = 0; index < n; ++index) {
        Node node = Node(index);
        auto node_ptr = std::make_shared<Node>(node);
        elementary_nodes.push_back(node_ptr);
    }

    for (const auto &edge : edge_list_) {
        EdgeWeighted edge_weighted = EdgeWeighted(elementary_nodes[std::get<0>(edge)], elementary_nodes[std::get<1>(edge)], std::get<2>(edge));
        auto edge_ptr = std::make_shared<EdgeWeighted>(edge_weighted);
        elementary_nodes[std::get<0>(edge)]->neighbors.push_back(edge_ptr);
        elementary_nodes[std::get<1>(edge)]->neighbors.push_back(edge_ptr);
    }
}

void SolverWeighted::FindMinPerfectMatching() {
    GreedyInit();

    auto roots = Roots();
    Tree tree = Tree(roots.front());

    // TODO don't forget to dissolve blossoms in the end
}

void SolverWeighted::PrintElementaryAdjList() {
    std::cout << "Adjacency list (to, weight, slack, matched):" << std::endl;
    for (const auto& vertex : elementary_nodes) {
        std::cout << vertex->index << ": y_v = " << vertex->DualVariableQuadrupled() / 4. << ", ";
        for (const auto& edge : vertex->neighbors) {
            std::cout << "(" << edge->OtherElementary(vertex)->index << ", " << edge->weight << ", " << edge->slack_quadrupled / 4. << ", " << edge->Matched() << ") ";
        }
        std::cout << std::endl;
    }
}

void SolverWeighted::GreedyInit() {
    // first, make all the slacks non-negative
    for (const auto& vertex : elementary_nodes) {
        if (vertex->neighbors.empty()) {
            continue;
        }

        int min_weight = vertex->neighbors.front()->weight;
        for (const auto& edge : vertex->neighbors) {
            if (edge->weight < min_weight) {
                min_weight = edge->weight;
            }
        }

        vertex->IncreaseDualVariableQuadrupled(min_weight * 2);
    }

    for (const auto& vertex : elementary_nodes) {
        if ((vertex->matched_edge) || (vertex->neighbors.empty())) {
            continue;
        }

        std::shared_ptr<EdgeWeighted> smallest_slack_edge = vertex->neighbors.front();
        for (const auto& edge : vertex->neighbors) {
            if (edge->slack_quadrupled < smallest_slack_edge->slack_quadrupled) {
                smallest_slack_edge = edge;
            }
        }

        int smallest_slack_quadrupled = smallest_slack_edge->slack_quadrupled;
        vertex->IncreaseDualVariableQuadrupled(smallest_slack_quadrupled);

        if (!smallest_slack_edge->OtherElementary(vertex)->matched_edge) {
            // if the other vertex is also unmatched, match the edge
            smallest_slack_edge->MakeMatched();
        }
    }
}

std::vector<std::shared_ptr<Node>> SolverWeighted::Roots() {
    std::vector<std::shared_ptr<Node>> roots;

    for (const auto& vertex : elementary_nodes) {
        if (vertex->parent_blossom) {
            throw std::runtime_error("Roots() must be called if there are no supernodes");
        }

        if (!vertex->matched_edge) {
            roots.push_back(vertex);
        }
    }

    return roots;
}

std::vector<std::pair<int, int>> SolverWeighted::Matching() {
    std::vector<std::pair<int, int>> matching;
    matching.reserve(elementary_nodes.size() / 2);

    for (const auto& vertex : elementary_nodes) {
        if (!vertex->matched_edge) {
            throw std::runtime_error("In Matching(): Matching is not perfect");
        }

        auto to = vertex->matched_edge->OtherElementary(vertex);
        if (to > vertex) {
            matching.emplace_back(vertex->index, to->index);
        }
    }

    return matching;
}
