#include <queue>
#include <unordered_set>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <list>

#include "Tree.h"
#include "EdgeWeighted.h"
#include "Node.h"

Tree::Tree(Node *root_,
           std::list<Node> *blossom_storage_,
           std::unordered_map<Node *, std::list<Node>::iterator> *iter_to_self_,
           int num_elementary_nodes_) : root(root_),
                                        blossom_storage(blossom_storage_), iter_to_self(iter_to_self_),
                                        num_elementary_nodes(num_elementary_nodes_), dual_var_quadrupled(0) {
    if (!root->blossom_children.empty()) {
        throw std::runtime_error("In Tree: root must be elementary");
    }

    root->plus = true;
    root->tree = this;
}

void Tree::Grow(EdgeWeighted &edge) {
    auto endpoints = edge.Endpoints();
    Node *parent = &endpoints.first;
    Node *child = &endpoints.second;
    if (child->tree == this) {
        std::swap(parent, child);
    }

    //std::cout << "grow " << parent->index << " " << child->index << " " << child->matched_edge->OtherEnd(*child).index << std::endl;

    if (parent->tree != this) {
        throw std::runtime_error("In Tree::Grow: parent vertex is not in this tree");
    }
    if (!parent->plus) {
        throw std::runtime_error("In Tree::Grow: parent is not a plus");
    }
    if (child->tree) {
        throw std::runtime_error("In Tree::Grow: child vertex is not free");
    }

    parent->tree_children.push_back(&edge);

    child->tree_parent = &edge;
    child->tree = this;
    child->minus = true;
    child->dual_var_quadrupled_amortized += dual_var_quadrupled;

    Node &grandchild = child->matched_edge->OtherEnd(*child);
    grandchild.tree_parent = child->matched_edge;
    child->tree_children = {child->matched_edge};
    grandchild.tree = this;
    grandchild.plus = true;
    grandchild.dual_var_quadrupled_amortized -= dual_var_quadrupled;
}

void Tree::Shrink(EdgeWeighted &edge_plus_plus) const {
    auto endpoints = edge_plus_plus.Endpoints();
    const Node *first = &endpoints.first;
    const Node *second = &endpoints.second;

    //std::cout << "shrink " << first->index << " " << second->index << std::endl;

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

    int next_blossom_index = num_elementary_nodes;
    if (!blossom_storage->empty()) {
        next_blossom_index = blossom_storage->back().index + 1;
    }
    blossom_storage->emplace_back(blossom_edges, next_blossom_index);
    (*iter_to_self)[&blossom_storage->back()] = std::prev(blossom_storage->end());
}

void Tree::Expand(Node &blossom) const {
    //std::cout << "expand" << std::endl;

    if (blossom.tree != this) {
        throw std::runtime_error("In Tree::Expand: supervertex is not in this tree");
    }
    if (blossom.blossom_children.empty()) {
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

Tree *Tree::Augment(EdgeWeighted &edge) {
    // TODO amortize tree dissolves

    //std::cout << "augmenting edge with slack " << edge.SlackQuadrupled() << std::endl;

    auto endpoints = edge.Endpoints();
    Node *parent = &endpoints.first;
    Node *child = &endpoints.second;
    if (parent->tree != this) {
        std::swap(parent, child);
    }

    if (parent->tree != this) {
        throw std::runtime_error("In Augment: parent vertex is not in this tree");
    }

    Node::MakeMatched(edge);
    // std::cout << "make matched " << parent->index << " " << child->index << std::endl;

    AugmentFromNode(*parent);
    Tree *other_tree = child->tree;
    other_tree->AugmentFromNode(*child);

    return other_tree;
}

void Tree::PrintTree() {
    std::cout << "Tree structure:" << std::endl;
    std::cout << "Root: " << root->index;
    std::cout << ", dual variable: " << 1. * dual_var_quadrupled / 4 << std::endl;

    struct Frame {
        const Node *node;
        std::string prefix; // prefix to print before connector / node value
        bool hasConnector; // whether to print "├── " or "└── " before the node value
        bool isLast; // if hasConnector==true, whether this node is the last sibling
    };

    std::vector<Frame> stack;
    stack.push_back({&root->TopBlossom(), "", false, false});

    while (!stack.empty()) {
        Frame f = stack.back();
        stack.pop_back();

        if (f.hasConnector) {
            std::cout << f.prefix << (f.isLast ? "└── " : "├── ");
        } else {
            std::cout << f.prefix; // usually empty for root
        }
        std::cout << f.node->index << '\n';

        const std::vector<EdgeWeighted *> &ch = f.node->tree_children;
        for (int i = static_cast<int>(ch.size()) - 1; i >= 0; --i) {
            bool childIsLast = (i == static_cast<int>(ch.size()) - 1);

            std::string childPrefix = f.prefix;
            if (f.hasConnector) {
                childPrefix += (f.isLast ? "    " : "│   ");
            }

            stack.push_back({&ch[i]->OtherEnd(*f.node), childPrefix, true, childIsLast});
        }
    }
}

Tree *Tree::MakePrimalUpdates() {
    // makes primal updates while we can, or until an augmentation was made
    // returns a pointer to the other tree in case of augmentation

    if (root->TopBlossom().tree == nullptr) {
        throw std::runtime_error("In Tree::MakePrimalUpdates: we are in the deleted tree");
    }

    int operation_cnt = 0;
    bool success = true;
    while (success) {
        // TODO experiment with the order of the operations

        ++operation_cnt;

        EdgeWeighted *augmentable_edge = AugmentableEdge();
        if (augmentable_edge) {
            Tree *other = Augment(*augmentable_edge);
            // std::cout << operation_cnt << " primal updates for this tree" << std::endl;
            return other;
        }

        EdgeWeighted *growable_edge = GrowableEdge();
        if (growable_edge) {
            Grow(*growable_edge);
            continue;
        }

        EdgeWeighted *shrinkable_edge = ShrinkableEdge();
        if (shrinkable_edge) {
            Shrink(*shrinkable_edge);
            continue;
        }

        Node *expandable_blossom = ExpandableBlossom();
        if (expandable_blossom) {
            Expand(*expandable_blossom);
            continue;
        }

        --operation_cnt;

        success = false;
    }

    // std::cout << operation_cnt << " primal updates for this tree" << std::endl;

    return nullptr;
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
        vertex->ClearDuringTreeDissolve();
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

Node &Tree::LCA(const EdgeWeighted &edge_plus_plus) const {
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

void Tree::AugmentFromNode(Node &vertex) {
    // std::cout << "augment from node " << vertex.index << std::endl;

    if (vertex.tree != this) {
        throw std::runtime_error("In Tree::AugmentFromNode: vertex is not in this tree");
    }

    bool match = false;
    const std::vector<EdgeWeighted *> path = PathToRoot(vertex);
    for (EdgeWeighted *edge_from_path : path) {
        if (match) {
            Node::MakeMatched(*edge_from_path);
            // std::cout << "make matched " << edge_from_path->ElementaryEndpoints().first.index << " " << edge_from_path->
            //     ElementaryEndpoints().second.index << std::endl;
        } else {
            Node::MakeUnmatched(*edge_from_path);
            // std::cout << "make unmatched " << edge_from_path->ElementaryEndpoints().first.index << " " << edge_from_path
            //     ->ElementaryEndpoints().second.index << std::endl;
        }
        match = !match;
    }

    DissolveTree();
}

EdgeWeighted *Tree::GrowableEdge() const {
    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        const Node *node = queue.front();
        queue.pop();

        if (node->plus) {
            for (EdgeWeighted *edge : node->neighbors) {
                if ((edge->SlackQuadrupled() == 0) && (edge->OtherEnd(*node).tree == nullptr)) {
                    return edge;
                }
            }
        }

        for (const EdgeWeighted *child : node->tree_children) {
            queue.push(&child->OtherEnd(*node));
        }
    }

    return nullptr;
}

EdgeWeighted *Tree::AugmentableEdge() const {
    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        const Node *node = queue.front();
        queue.pop();

        if (node->plus) {
            for (EdgeWeighted *edge : node->neighbors) {
                if (edge->SlackQuadrupled() == 0) {
                    Node &other_end = edge->OtherEnd(*node);
                    if ((other_end.tree != this) && (other_end.plus) && (other_end.tree != nullptr)) {
                        return edge;
                    }
                }
            }
        }

        for (const EdgeWeighted *child : node->tree_children) {
            queue.push(&child->OtherEnd(*node));
        }
    }

    return nullptr;
}

EdgeWeighted *Tree::ShrinkableEdge() const {
    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        const Node *node = queue.front();
        queue.pop();

        if (node->plus) {
            for (EdgeWeighted *edge : node->neighbors) {
                if (edge->SlackQuadrupled() == 0) {
                    Node &other_end = edge->OtherEnd(*node);
                    if ((other_end.tree == this) && (other_end.plus)) {
                        return edge;
                    }
                }
            }
        }

        for (const EdgeWeighted *child : node->tree_children) {
            queue.push(&child->OtherEnd(*node));
        }
    }

    return nullptr;
}

Node *Tree::ExpandableBlossom() {
    // returns a blossom that is a minus and has zero dual variable
    // if there is no such blossom, returns nullptr

    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        Node *node = queue.front();
        queue.pop();

        if ((node->minus) && (!node->blossom_children.empty()) && (node->DualVariableQuadrupled() == 0)) {
            return node;
        }

        for (EdgeWeighted *child_edge : node->tree_children) {
            queue.push(&child_edge->OtherEnd(*node));
        }
    }

    return nullptr;
}

int Tree::PlusEmptySlack() const {
    // returns the minimal quadrupled slack of a (+, empty) edge

    if (root->TopBlossom().tree == nullptr) {
        throw std::runtime_error("In Tree::MakePrimalUpdates: we are in the deleted tree");
    }

    int answer = INT32_MAX;

    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        const Node *node = queue.front();
        queue.pop();

        if (node->plus) {
            for (EdgeWeighted *edge : node->neighbors) {
                if ((edge->SlackQuadrupled() < answer) && (edge->OtherEnd(*node).tree == nullptr)) {
                    answer = edge->SlackQuadrupled();
                }
            }
        }

        for (const EdgeWeighted *child : node->tree_children) {
            queue.push(&child->OtherEnd(*node));
        }
    }

    return answer;
}

int Tree::PlusPlusExternalSlack() const {
    // returns the minimal quadrupled slack of a (+, +) edge, where the other vertex is in another tree
    int answer = INT32_MAX;

    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        const Node *node = queue.front();
        queue.pop();

        if (node->plus) {
            for (EdgeWeighted *edge : node->neighbors) {
                if (edge->SlackQuadrupled() < answer) {
                    Node &other_end = edge->OtherEnd(*node);
                    if ((other_end.tree != this) && (other_end.plus) && (other_end.tree != nullptr)) {
                        answer = edge->SlackQuadrupled();
                    }
                }
            }
        }

        for (const EdgeWeighted *child : node->tree_children) {
            queue.push(&child->OtherEnd(*node));
        }
    }

    return answer;
}

int Tree::PlusPlusInternalSlack() const {
    // returns the minimal quadrupled slack of a (+, +) edge, where the other vertex is in the same tree
    int answer = INT32_MAX;

    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        const Node *node = queue.front();
        queue.pop();

        if (node->plus) {
            for (EdgeWeighted *edge : node->neighbors) {
                if (edge->SlackQuadrupled() < answer) {
                    Node &other_end = edge->OtherEnd(*node);
                    if ((other_end.tree == this) && (other_end.plus)) {
                        answer = edge->SlackQuadrupled();
                    }
                }
            }
        }

        for (const EdgeWeighted *child : node->tree_children) {
            queue.push(&child->OtherEnd(*node));
        }
    }

    return answer;
}

int Tree::PlusMinusExternalSlack() const {
    // returns the minimal quadrupled slack of a (+, -) edge, where the other vertex is in another tree
    int answer = INT32_MAX;

    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        const Node *node = queue.front();
        queue.pop();

        if (node->plus) {
            for (EdgeWeighted *edge : node->neighbors) {
                if (edge->SlackQuadrupled() < answer) {
                    Node &other_end = edge->OtherEnd(*node);
                    if ((other_end.tree != this) && (other_end.minus)) {
                        answer = edge->SlackQuadrupled();
                    }
                }
            }
        }

        for (const EdgeWeighted *child : node->tree_children) {
            queue.push(&child->OtherEnd(*node));
        }
    }

    return answer;
}

int Tree::MinMinusBlossomVariable() const {
    // returns the minimum quadrupled dual variable of a (-) non-elementary blossom

    int answer = INT32_MAX;

    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        Node *node = queue.front();
        queue.pop();

        if ((node->minus) && (!node->blossom_children.empty()) && (node->DualVariableQuadrupled() < answer)) {
            answer = node->DualVariableQuadrupled();
        }

        for (const EdgeWeighted *child : node->tree_children) {
            queue.push(&child->OtherEnd(*node));
        }
    }

    return answer;
}
