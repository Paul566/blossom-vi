#include <queue>
#include <unordered_set>
#include "Tree.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <list>

Tree::Tree(Node *root_,
           std::list<Node> *blossom_storage_,
           std::unordered_map<Node *, std::list<Node>::iterator> *iter_to_self_) : root(root_),
    blossom_storage(blossom_storage_), iter_to_self(iter_to_self_) {
    if (root_->index == -1) {
        throw std::runtime_error("In Tree: root must be elementary");
    }

    root->plus = true;
    root->tree_root = root;
}

void Tree::Grow(EdgeWeighted &edge) const {
    auto endpoints = edge.Endpoints();
    Node *parent = &endpoints.first;
    Node *child = &endpoints.second;
    if (child->tree_root == root) {
        std::swap(parent, child);
    }

    // std::cout << "grow " << parent->index << " " << child->index << " " << child->matched_edge->OtherEnd(*child).index
    // << std::endl;

    if (parent->index == 5 && child->index == -1) {
        // std::cout << "a;lkjfds" << std::endl;
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

    parent->tree_children.push_back(&edge);

    child->tree_parent = &edge;
    child->tree_root = root;
    child->minus = true;

    Node &grandchild = child->matched_edge->OtherEnd(*child);
    grandchild.tree_parent = child->matched_edge;
    child->tree_children = {child->matched_edge};
    grandchild.tree_root = root;
    grandchild.plus = true;
}

void Tree::Shrink(EdgeWeighted &edge_plus_plus) const {
    auto endpoints = edge_plus_plus.Endpoints();
    const Node *first = &endpoints.first;
    const Node *second = &endpoints.second;

    // std::cout << "shrink " << first->index << " " << second->index << std::endl;

    const Node &lca = LCA(edge_plus_plus);

    std::vector<EdgeWeighted *> blossom_edges;
    while (first != &lca) {
        blossom_edges.push_back(first->tree_parent);
        first = &first->tree_parent->OtherEnd(*first);
    }
    std::reverse(blossom_edges.begin(), blossom_edges.end());

    blossom_edges.push_back(&edge_plus_plus);

    while (second != &lca) {
        blossom_edges.push_back(second->tree_parent);
        second = &second->tree_parent->OtherEnd(*second);
    }

    const int max_blossom_index = blossom_storage->back().index;
    blossom_storage->emplace_back(blossom_edges, max_blossom_index + 1);
    (*iter_to_self)[&blossom_storage->back()] = std::prev(blossom_storage->end());
}

void Tree::Expand(Node &blossom) const {
    // std::cout << "expand" << std::endl;

    if (blossom.tree_root != root) {
        throw std::runtime_error("In Tree::Expand: supervertex is not in this tree");
    }
    if (blossom.index != -1) {
        throw std::runtime_error("In Tree::Expand: supervertex must be not elementary");
    }
    if (blossom.blossom_parent) {
        throw std::runtime_error("In Tree::Expand: supervertex must be a top blossom");
    }
    if ((!blossom.minus) || (blossom.plus)) {
        throw std::runtime_error("In Tree::Expand: supervertex is not a minus");
    }
    if (blossom.DualVariableQuadrupled() != 0) {
        throw std::runtime_error("In Tree::Expand: the supervertex has to have zero dual variable");
    }

    blossom.Dissolve();
    std::list<Node>::iterator it = (*iter_to_self)[&blossom];
    iter_to_self->erase(&blossom);
    blossom_storage->erase(it);
}

void Tree::Augment(EdgeWeighted &edge) {
    // TODO change it when moving to multiple trees

    auto endpoints = edge.Endpoints();
    Node *parent = &endpoints.first;
    Node *child = &endpoints.second;
    if (child->tree_root) {
        std::swap(parent, child);
    }

    // std::cout << "augment " << parent->index << " " << child->index << std::endl;

    if (parent->index == -1 && child->index == 4) {
        // std::cout << "da;lskjf" << std::endl;
    }

    if (child->tree_root) {
        throw std::runtime_error("In Tree::Augment: child must be not in a tree");
    }
    if (child->matched_edge) {
        throw std::runtime_error("In Tree::Augment: child must be unmatched");
    }

    Node::MakeMatched(edge);

    bool match = false;
    const std::vector<EdgeWeighted *> path = PathToRoot(*parent);
    for (EdgeWeighted *edge_from_path : path) {
        if (match) {
            Node::MakeMatched(*edge_from_path);
        } else {
            Node::MakeUnmatched(*edge_from_path);
        }
        match = !match;
    }

    DissolveTree();
}

bool Tree::MakePrimalUpdate(bool *is_augmented) {
    // tries to make a primal update, returns true if a non-augmentation update was made

    // TODO when you will have multiple trees, consider the case of a (+,-) edge for edge_min_slack

    EdgeWeighted &edge_min_slack = MinSlackEdgeFromPlus();
    if (edge_min_slack.slack_quadrupled == 0) {
        auto endpoints = edge_min_slack.Endpoints();
        Node *vertex_in_tree = &endpoints.first;
        Node *vertex_other = &endpoints.second;
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
            *is_augmented = true;
            Augment(edge_min_slack);
            return false;
        }
    }

    Node *expandable_blossom = ExpandableBlossom();

    if (expandable_blossom == nullptr) {
        return false;
    }

    Expand(*expandable_blossom);
    return true;
}

void Tree::MakeDualUpdate() {
    // TODO this will change with multiple trees

    EdgeWeighted &edge_min_slack = MinSlackEdgeFromPlus();
    const Node *blossom_min_var = MinYMinusBlossom();

    int increment_quadrupled = INT32_MAX;
    if (blossom_min_var) {
        increment_quadrupled = blossom_min_var->DualVariableQuadrupled();
    }

    auto [first, second] = edge_min_slack.Endpoints();
    if ((first.tree_root == root) && (second.tree_root == root) && (first.plus) && (second.plus)) {
        if (edge_min_slack.slack_quadrupled % 2 != 0) {
            throw std::runtime_error("In MakeDualUpdate: a ++ edge with odd quadrupled slack");
        }

        if (increment_quadrupled > edge_min_slack.slack_quadrupled / 2) {
            increment_quadrupled = edge_min_slack.slack_quadrupled / 2;
        }
    } else {
        if (increment_quadrupled > edge_min_slack.slack_quadrupled) {
            increment_quadrupled = edge_min_slack.slack_quadrupled;
        }
    }

    ChangeDualVariables(increment_quadrupled);
}

void Tree::DissolveTree() {
    // clears the tree structure from the nodes

    std::vector<Node *> all_vertices;

    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        auto cur_vertex = queue.front();
        queue.pop();
        all_vertices.push_back(cur_vertex);
        for (const EdgeWeighted *child : cur_vertex->tree_children) {
            queue.push(&child->OtherEnd(*cur_vertex));
        }
    }

    for (Node *vertex : all_vertices) {
        vertex->tree_root = nullptr;
        vertex->tree_children.clear();
        vertex->tree_parent = nullptr;
        vertex->plus = false;
        vertex->minus = false;
    }
}

std::vector<EdgeWeighted *> Tree::PathToRoot(const Node &vertex) {
    // // std::cout << "path to root" << std::endl;

    if (vertex.blossom_parent) {
        throw std::runtime_error("In Tree::PathToRoot: vertex is not a top blossom");
    }

    const Node *cur_vertex = &vertex;
    std::vector<EdgeWeighted *> path;
    while (cur_vertex->tree_parent) {
        auto [head, tail] = cur_vertex->tree_parent->Endpoints();
        if (head.blossom_parent || tail.blossom_parent) {
            if (head.blossom_parent != tail.blossom_parent) {
                throw std::runtime_error("In PathToRoot: found an incorrect edge");
            }
        }

        path.push_back(cur_vertex->tree_parent);
        cur_vertex = &cur_vertex->tree_parent->OtherEnd(*cur_vertex);
    }

    return path;
}

Node &Tree::LCA(EdgeWeighted &edge_plus_plus) const {
    // TODO can be made better
    auto endpoints = edge_plus_plus.Endpoints();
    Node *first = &endpoints.first;
    Node *second = &endpoints.second;

    Node &root_blossom = root->TopBlossom();

    std::unordered_set<Node *> visited;
    visited.insert(first);
    while (first != &root_blossom) {
        first = &first->tree_parent->OtherEnd(*first);
        visited.insert(first);
    }

    while (second != &root_blossom) {
        if (visited.contains(second)) {
            return *second;
        }
        second = &second->tree_parent->OtherEnd(*second);
    }

    return root_blossom;
}

EdgeWeighted &Tree::MinSlackEdgeFromPlus() const {
    // returns a (+, smth) edge with minimal slack that is not an edge in the tree
    // TODO refactor into min slack edges of different types

    int min_slack = INT32_MAX;
    EdgeWeighted *min_slack_edge = nullptr;

    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        const Node *node = queue.front();
        queue.pop();

        if (node->plus) {
            for (EdgeWeighted *edge : node->neighbors) {
                if ((edge->OtherEnd(*node).minus) && (edge->OtherEnd(*node).tree_root == root)) {
                    // this is a (+, -) edge inside the tree
                    continue;
                }

                if (edge->slack_quadrupled < min_slack) {
                    min_slack = edge->slack_quadrupled;
                    min_slack_edge = edge;

                    if (min_slack == 0) {
                        return *min_slack_edge;
                    }
                }
            }
        }

        for (const EdgeWeighted *child : node->tree_children) {
            queue.push(&child->OtherEnd(*node));
        }
    }

    if (!min_slack_edge) {
        throw std::runtime_error("No perfect matching exists (In MinSlackEdgeFromPlus: edge from plus not found)");
    }

    return *min_slack_edge;
}

Node *Tree::ExpandableBlossom() {
    // returns a blossom that is a minus and has zero dual variable
    // if there is no such blossom, returns nullptr

    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        Node *node = queue.front();
        queue.pop();

        if ((node->minus) && (node->index == -1) && (node->DualVariableQuadrupled() == 0)) {
            return node;
        }

        for (EdgeWeighted *child_edge : node->tree_children) {
            queue.push(&child_edge->OtherEnd(*node));
        }
    }

    return nullptr;
}

Node *Tree::MinYMinusBlossom() const {
    // returns a blossom that is a minus and has minimal dual variable
    // if there are no minus blossoms, returns nullptr

    Node *answer = nullptr;
    int min_variable = INT32_MAX;

    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        Node *node = queue.front();
        queue.pop();

        if ((node->minus) && (node->index == -1) && (node->DualVariableQuadrupled() < min_variable)) {
            min_variable = node->DualVariableQuadrupled();
            answer = node;
        }

        for (const EdgeWeighted *child : node->tree_children) {
            queue.push(&child->OtherEnd(*node));
        }
    }

    return answer;
}

void Tree::ChangeDualVariables(const int increment) const {
    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        Node *node = queue.front();
        queue.pop();

        if (node->plus) {
            node->IncreaseDualVariableQuadrupled(increment);
        } else {
            node->IncreaseDualVariableQuadrupled(-increment);
        }

        for (EdgeWeighted *child : node->tree_children) {
            queue.push(&child->OtherEnd(*node));
        }
    }
}
