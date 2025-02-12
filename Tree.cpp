#include <queue>
#include <unordered_set>
#include "Tree.h"


Tree::Tree(const std::shared_ptr<Node> &root_) {
    if (root_->index == -1) {
        throw std::runtime_error("In Tree: root must be elementary");
    }

    root = root_;
}

void Tree::Grow(const std::shared_ptr<EdgeWeighted>& edge) {
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

void Tree::Shrink(const std::shared_ptr<EdgeWeighted>& edge_plus_plus) {
    auto lca = LCA(edge_plus_plus);

    std::vector<std::shared_ptr<Node>> blossom_vertices;
    auto [first, second] = edge_plus_plus->VerticesTopBlossoms();
    while (first != lca) {
        blossom_vertices.push_back(first);
        first = first->tree_parent->OtherBlossom(first);
    }
    while (second != lca) {
        blossom_vertices.push_back(second);
        second = second->tree_parent->OtherBlossom(second);
    }
    blossom_vertices.push_back(lca);

    Node blossom(blossom_vertices);
}

void Tree::Expand(const std::shared_ptr<Node> &supervertex) {
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

    for (const auto& blossom_child : supervertex->children_blossom) {
        blossom_child->parent_blossom = nullptr;
    }

}

std::shared_ptr<Node> Tree::LCA(const std::shared_ptr<EdgeWeighted>& edge_plus_plus) const {
    // TODO can be made better
    auto [first, second] = edge_plus_plus->VerticesTopBlossoms();

    auto root_blossom = root->TopBlossom();

    std::unordered_set<std::shared_ptr<Node>> visited;
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

    std::queue<std::shared_ptr<Node>> queue;
    queue.push(root);
    while (!queue.empty()) {
        std::shared_ptr<Node> node = queue.front();
        queue.pop();

        if (node->plus) {
            for (const auto& edge : node->neighbors) {
                if (edge->slack_quadrupled < min_slack) {
                    min_slack = edge->slack_quadrupled;
                    min_slack_edge = edge;
                }
            }
        }

        for (const auto& child : node->tree_children) {
            queue.push(child->OtherBlossom(node));
        }
    }

    return min_slack_edge;
}
