#include <queue>
#include <unordered_set>
#include "Tree.h"

#include <algorithm>
#include <iostream>

Tree::Tree(const std::shared_ptr<Node> &root_) {
    if (root_->index == -1) {
        throw std::runtime_error("In Tree: root must be elementary");
    }

    root = root_;
    root->plus = true;
    root->tree_root = root;
}

void Tree::Grow(const std::shared_ptr<EdgeWeighted> &edge) const {
    auto [parent, child] = edge->VerticesTopBlossoms();
    if (child->tree_root == root) {
        std::swap(parent, child);
    }

    std::cout << "grow " << parent->index << " " << child->index << std::endl;

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

    const auto grandchild = child->matched_edge->OtherBlossom(child);
    grandchild->tree_parent = child->matched_edge;
    child->tree_children = {child->matched_edge};
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
    std::cout << "augment" << std::endl;

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

bool Tree::MakePrimalUpdate() {
    // tries to make a primal update, returns true if the update was made

    // TODO when you will have multiple trees, consider the case of a (+,-) edge for edge_min_slack

    auto edge_min_slack = MinSlackEdgeFromPlus();
    if (edge_min_slack->slack_quadrupled == 0) {
        auto [vertex_in_tree, vertex_other] = edge_min_slack->VerticesTopBlossoms();
        if (vertex_in_tree->tree_root != root) {
            std::swap(vertex_in_tree, vertex_other);
        }

        if (vertex_other->tree_root == root) {
            Shrink(edge_min_slack);
            return true;
        }

        if ((!vertex_other->tree_root) && (vertex_other->matched_edge)) {
            Grow(edge_min_slack);
            return true;
        }

        if (!vertex_other->matched_edge) {
            // TODO change when multiple trees
            Augment(edge_min_slack);
            return true;
        }
    }

    auto expandable_blossom = ExpandableBlossom();

    if (!expandable_blossom) {
        return false;
    }

    Expand(expandable_blossom);
    return true;
}

void Tree::MakeDualUpdate() {
    // TODO this will change with multiple trees

    auto edge_min_slack = MinSlackEdgeFromPlus();
    auto blossom_min_var = MinYMinusBlossom();

    int increment_quadrupled = INT32_MAX;
    if (blossom_min_var) {
        increment_quadrupled = blossom_min_var->DualVariableQuadrupled();
    }

    auto [first, second] = edge_min_slack->VerticesTopBlossoms();
    if ((first->tree_root == root) && (second->tree_root == root) && (first->plus) && (second->plus)) {
        if (edge_min_slack->slack_quadrupled % 2 != 0) {
            throw std::runtime_error("In MakeDualUpdate: a ++ edge with odd quadrupled slack");
        }

        if (increment_quadrupled > edge_min_slack->slack_quadrupled / 2) {
            increment_quadrupled = edge_min_slack->slack_quadrupled / 2;
        }
    } else {
        if (increment_quadrupled > edge_min_slack->slack_quadrupled) {
            increment_quadrupled = edge_min_slack->slack_quadrupled;
        }
    }

    ChangeDualVariables(increment_quadrupled);
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
    // returns a (+, smth) edge with minimal slack that is not an edge in the tree
    // TODO refactor into min slack edges of different types

    int min_slack = INT_MAX;
    std::shared_ptr<EdgeWeighted> min_slack_edge = nullptr;

    std::queue<std::shared_ptr<Node> > queue;
    queue.push(root->TopBlossom());
    while (!queue.empty()) {
        std::shared_ptr<Node> node = queue.front();
        queue.pop();

        if (node->plus) {
            for (const auto &edge : node->neighbors) {
                if ((edge->OtherBlossom(node)->minus) && (edge->OtherBlossom(node)->tree_root == root)) {
                    // this is a (+, -) edge inside the tree
                    continue;
                }

                if (edge->slack_quadrupled < min_slack) {
                    min_slack = edge->slack_quadrupled;
                    min_slack_edge = edge;

                    if (min_slack == 0) {
                        return min_slack_edge;
                    }
                }
            }
        }

        for (const auto &child : node->tree_children) {
            queue.push(child->OtherBlossom(node));
        }
    }

    if (!min_slack_edge) {
        throw std::runtime_error("In MinSlackEdgeFromPlus: edge from plus not found");
    }

    return min_slack_edge;
}

std::shared_ptr<Node> Tree::ExpandableBlossom() {
    // returns a blossom that is a minus and has zero dual variable
    // if there is no such blossom, returns nullptr

    std::queue<std::shared_ptr<Node> > queue;
    queue.push(root->TopBlossom());
    while (!queue.empty()) {
        std::shared_ptr<Node> node = queue.front();
        queue.pop();

        if ((node->minus) && (node->index == -1) && (node->DualVariableQuadrupled() == 0)) {
            return node;
        }

        for (const auto &child : node->tree_children) {
            queue.push(child->OtherBlossom(node));
        }
    }

    return nullptr;
}

std::shared_ptr<Node> Tree::MinYMinusBlossom() {
    // returns a blossom that is a minus and has minimal dual variable
    // if there are no minus blossoms, returns nullptr

    std::shared_ptr<Node> answer = nullptr;
    int min_variable = INT32_MAX;

    std::queue<std::shared_ptr<Node> > queue;
    queue.push(root->TopBlossom());
    while (!queue.empty()) {
        std::shared_ptr<Node> node = queue.front();
        queue.pop();

        if ((node->minus) && (node->index == -1) && (node->DualVariableQuadrupled() < min_variable)) {
            min_variable = node->DualVariableQuadrupled();
            answer = node;
        }

        for (const auto &child : node->tree_children) {
            queue.push(child->OtherBlossom(node));
        }
    }

    return answer;
}

void Tree::ChangeDualVariables(const int increment) const {
    std::queue<std::shared_ptr<Node> > queue;
    queue.push(root->TopBlossom());
    while (!queue.empty()) {
        std::shared_ptr<Node> node = queue.front();
        queue.pop();

        if (node->plus) {
            node->IncreaseDualVariableQuadrupled(increment);
        } else {
            node->IncreaseDualVariableQuadrupled(-increment);
        }

        for (const auto &child : node->tree_children) {
            queue.push(child->OtherBlossom(node));
        }
    }
}
