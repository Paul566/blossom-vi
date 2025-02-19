#include <queue>
#include <unordered_set>
#include "Tree.h"

#include <algorithm>

Tree::Tree(const std::shared_ptr<Node> &root_) {
    if (root_->index == -1) {
        throw std::runtime_error("In Tree: root must be elementary");
    }

    root = root_;
}

void Tree::Grow(const std::shared_ptr<EdgeWeighted> &edge) const {
    auto [parent, child] = edge->VerticesTopBlossoms();
    if (child->tree_root == root) {
        std::swap(parent, child);
    }

    if (parent->tree_root != root) {
        throw std::runtime_error("In Tree::Grow: parent vertex is not in this tree");
    }
    if (!parent->plus) {
        throw std::runtime_error("In Tree::Grow: parent is not a plus");
    }
    if (child->tree_root) {
        throw std::runtime_error("In Tree::Grow: child vertex is not free");
    }

    parent->tree_children.push_back(edge);

    child->tree_parent = edge;
    child->tree_root = root;
    child->minus = true;

    auto grandchild = child->matched_edge->OtherBlossom(child);
    grandchild->tree_parent = child->matched_edge;
    grandchild->tree_root = root;
    grandchild->plus = true;
}

void Tree::Shrink(const std::shared_ptr<EdgeWeighted> &edge_plus_plus) const {
    auto lca = LCA(edge_plus_plus);

    std::vector<std::shared_ptr<EdgeWeighted>> blossom_edges;
    auto [first, second] = edge_plus_plus->VerticesTopBlossoms();

    while (first != lca) {
        blossom_edges.push_back(first->tree_parent);
        first = first->tree_parent->OtherBlossom(first);
    }
    std::reverse(blossom_edges.begin(), blossom_edges.end());

    blossom_edges.push_back(edge_plus_plus);

    while (second != lca) {
        blossom_edges.push_back(second->tree_parent);
        second = second->tree_parent->OtherBlossom(second);
    }

    Node blossom(blossom_edges);
}

void Tree::Expand(const std::shared_ptr<Node> &supervertex) const {
    if (supervertex->tree_root != root) {
        throw std::runtime_error("In Tree::Expand: supervertex is not in this tree");
    }
    if (supervertex->index != -1) {
        throw std::runtime_error("In Tree::Expand: supervertex must be not elementary");
    }
    if (supervertex->parent_blossom) {
        throw std::runtime_error("In Tree::Expand: supervertex must be a top blossom");
    }
    if ((!supervertex->minus) || (supervertex->plus)) {
        throw std::runtime_error("In Tree::Expand: supervertex is not a minus");
    }
    if (supervertex->DualVariableQuadrupled() != 0) {
        throw std::runtime_error("In Tree::Expand: the supervertex has to have zero dual variable");
    }

    supervertex->Dissolve();
}

void Tree::Augment(const std::shared_ptr<EdgeWeighted> &edge) {
    // TODO change it when moving to multiple trees
    auto [parent, child] = edge->VerticesTopBlossoms();
    if (child->tree_root) {
        std::swap(parent, child);
    }
    if (child->tree_root) {
        throw std::runtime_error("In Tree::Augment: child must be not in a tree");
    }
    if (child->matched_edge) {
        throw std::runtime_error("In Tree::Augment: child must be unmatched");
    }

    auto path = PathToRoot(parent);
    edge->MakeMatched();
    bool match = false;
    for (const auto& edge_path : path) {
        if (match) {
            edge_path->MakeMatched();
        } else {
            edge_path->MakeUnmatched();
        }
        match = !match;
    }

    DissolveTree();
}

void Tree::DissolveTree() const {
    std::vector<std::shared_ptr<Node>> all_vertices;

    std::queue<std::shared_ptr<Node>> queue;
    queue.push(root->TopBlossom());
    while (!queue.empty()) {
        auto cur_vertex = queue.front();
        queue.pop();
        all_vertices.push_back(cur_vertex);
        for (const auto& child : cur_vertex->tree_children) {
            queue.push(child->OtherBlossom(cur_vertex));
        }
    }

    for (const auto& vertex : all_vertices) {
        vertex->tree_root = nullptr;
        vertex->tree_children.clear();
        vertex->tree_parent = nullptr;
        vertex->plus = false;
        vertex->minus = false;
    }
}

std::vector<std::shared_ptr<EdgeWeighted>> Tree::PathToRoot(std::shared_ptr<Node> vertex) {
    std::vector<std::shared_ptr<EdgeWeighted>> path;
    while (vertex->tree_parent) {
        path.push_back(vertex->tree_parent);
        vertex = vertex->tree_parent->OtherBlossom(vertex);
    }

    return path;
}

std::shared_ptr<Node> Tree::LCA(const std::shared_ptr<EdgeWeighted> &edge_plus_plus) const {
    // TODO can be made better
    auto [first, second] = edge_plus_plus->VerticesTopBlossoms();

    auto root_blossom = root->TopBlossom();

    std::unordered_set<std::shared_ptr<Node> > visited;
    visited.insert(first);
    while (first != root_blossom) {
        first = first->tree_parent->OtherBlossom(first);
        visited.insert(first);
    }

    while (second != root_blossom) {
        if (visited.contains(second)) {
            return second;
        }
        second = second->tree_parent->OtherBlossom(second);
    }

    return root_blossom;
}

std::shared_ptr<EdgeWeighted> Tree::MinSlackEdgeFromPlus() const {
    int min_slack = INT_MAX;
    std::shared_ptr<EdgeWeighted> min_slack_edge = nullptr;

    std::queue<std::shared_ptr<Node> > queue;
    queue.push(root);
    while (!queue.empty()) {
        std::shared_ptr<Node> node = queue.front();
        queue.pop();

        if (node->plus) {
            for (const auto &edge : node->neighbors) {
                if (edge->slack_quadrupled < min_slack) {
                    min_slack = edge->slack_quadrupled;
                    min_slack_edge = edge;
                }
            }
        }

        for (const auto &child : node->tree_children) {
            queue.push(child->OtherBlossom(node));
        }
    }

    return min_slack_edge;
}
