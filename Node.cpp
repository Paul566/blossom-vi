#include "Node.h"

#include <unordered_set>

#include "EdgeWeighted.h"


Node::Node(const int index_) : index(index_), children_blossom(std::vector<std::shared_ptr<Node>>()) {
    dual_variable_quadrupled = 0;
    parent_blossom = nullptr;
    matched_edge = nullptr;
    tree_parent = nullptr;
    tree_root = nullptr;
    tree_children = std::vector<std::shared_ptr<EdgeWeighted>>();
    neighbors = std::vector<std::shared_ptr<EdgeWeighted>>();
    plus = false;
    minus = false;
}

Node::Node(std::vector<std::shared_ptr<Node>> blossom_vertices) : index(-1), children_blossom(blossom_vertices) {
    if (blossom_vertices.size() % 2 == 0) {
        throw std::runtime_error("In Node: blossom has to have an odd number of vertices");
    }

    dual_variable_quadrupled = 0;
    parent_blossom = nullptr;
    plus = true;    // after Shrink, we must be a plus
    minus = false;

    // update the parent of the children
    for (const auto& vertex : blossom_vertices) {
        if (vertex->parent_blossom) {
            throw std::runtime_error("In Node: some vertex is not a top blossom");
        }
        vertex->parent_blossom = std::shared_ptr<Node>(this);
    }

    // find tree_root
    tree_root = blossom_vertices.front()->tree_root;
    for (const auto& vertex : blossom_vertices) {
        if (vertex->tree_root != tree_root) {
            throw std::runtime_error("In Node: not all the vertices have the same tree_root");
        }
    }

    std::unordered_set<std::shared_ptr<Node>> blossom_set;
    blossom_set.reserve(blossom_vertices.size());
    for (const auto& vertex : blossom_vertices) {
        blossom_set.insert(vertex);
    }

    // find neighbors
    for (const auto& vertex : blossom_vertices) {
        for (const auto& edge : vertex->neighbors) {
            if (!blossom_set.contains(edge->OtherBlossom(vertex))) {
                neighbors.push_back(edge);
            }
        }
    }

    // find matched_edge and tree_parent
    for (const auto& vertex : blossom_vertices) {
        if ((matched_edge) && (!blossom_set.contains(vertex->matched_edge->OtherBlossom(vertex)))) {
            throw std::runtime_error("In Node: there must be only one outgoing matched edge");
        }

        if (!blossom_set.contains(vertex->matched_edge->OtherBlossom(vertex))) {
            matched_edge = vertex->matched_edge;
            tree_parent = vertex->matched_edge;
        }
    }

    // find tree_children
    for (const auto& vertex : blossom_vertices) {
        for (auto& child_edge : vertex->tree_children) {
            if (!blossom_set.contains(child_edge->OtherBlossom(vertex))) {
                tree_children.push_back(child_edge);
            }
        }
    }

}

std::shared_ptr<Node> Node::TopBlossom() {
    if (parent_blossom == nullptr) {
        return std::shared_ptr<Node>(this);
    }
    return parent_blossom->TopBlossom();
}

int Node::DualVariableQuadrupled() const {
    return dual_variable_quadrupled;
}

void Node::IncreaseDualVariableQuadrupled(const int increment) {
    if (parent_blossom) {
        throw std::runtime_error("In IncreaseDualVariableQuadrupled: the vertex is not a top blossom");
    }

    dual_variable_quadrupled += increment;

    for (const auto& edge : neighbors) {
        edge->slack_quadrupled -= increment;
    }
}
