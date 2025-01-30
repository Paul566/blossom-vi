#include "SolverWeighted.h"

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

void SolverWeighted::PrintElementaryAdjList() {
    std::cout << "Adjacency list (to, weight, slack, matched):" << std::endl;
    for (auto vertex : elementary_nodes) {
        std::cout << vertex->index << ": y_v = " << vertex->dual_variable_quadrupled / 4. << ", ";
        for (const auto& edge : vertex->neighbors) {
            std::cout << "(" << edge->OtherNode(vertex)->index << ", " << edge->weight << ", " << edge->slack_quadrupled / 4. << ", " << edge->matched << ") ";
        }
        std::cout << std::endl;
    }
}

std::shared_ptr<Node> SolverWeighted::TopBlossom(std::shared_ptr<Node> vertex) {
    if (vertex->parent_blossom == nullptr) {
        return vertex;
    }
    return TopBlossom(vertex->parent_blossom);
}

void SolverWeighted::GreedyInit() {
    // first, make all the slacks non-negative
    for (auto vertex : elementary_nodes) {
        if (vertex->neighbors.empty()) {
            continue;
        }

        int min_weight = vertex->neighbors.front()->weight;
        for (auto edge : vertex->neighbors) {
            if (edge->weight < min_weight) {
                min_weight = edge->weight;
            }
        }

        vertex->dual_variable_quadrupled = min_weight * 2;
        for (auto edge : vertex->neighbors) {
            edge->slack_quadrupled -= vertex->dual_variable_quadrupled;
        }
    }

    for (auto vertex : elementary_nodes) {
        if ((vertex->matched_edge) || (vertex->neighbors.empty())) {
            continue;
        }

        std::shared_ptr<EdgeWeighted> smallest_slack_edge = vertex->neighbors.front();
        for (auto edge : vertex->neighbors) {
            if (edge->slack_quadrupled < smallest_slack_edge->slack_quadrupled) {
                smallest_slack_edge = edge;
            }
        }

        int smallest_slack_quadrupled = smallest_slack_edge->slack_quadrupled;
        vertex->dual_variable_quadrupled += smallest_slack_quadrupled;
        for (auto edge : vertex->neighbors) {
            edge->slack_quadrupled -= smallest_slack_quadrupled;
        }

        if (!smallest_slack_edge->OtherNode(vertex)->matched_edge) {
            smallest_slack_edge->matched = true;
            vertex->matched_edge = smallest_slack_edge;
            smallest_slack_edge->OtherNode(vertex)->matched_edge = smallest_slack_edge;
        }
    }
}

std::vector<std::shared_ptr<Node>> SolverWeighted::Roots() {
    std::vector<std::shared_ptr<Node>> roots;

    for (auto vertex : elementary_nodes) {
        if (vertex->parent_blossom) {
            throw std::runtime_error("Roots() must be called if there are no supernodes");
        }

        if (!vertex->matched_edge) {
            roots.push_back(vertex);
        }
    }

    return roots;
}
