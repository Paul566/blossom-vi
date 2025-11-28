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
                                        dual_var_quadrupled(0), num_elementary_nodes(num_elementary_nodes_),
                                        blossom_storage(blossom_storage_), iter_to_self(iter_to_self_) {
    if (!root->IsElementary()) {
        throw std::runtime_error("In Tree: root must be elementary");
    }

    root->MakeRootOfTree(*this);
}

void Tree::Grow(EdgeWeighted &edge) {
    auto endpoints = edge.Endpoints();
    Node *parent = &endpoints.first;
    Node *child = &endpoints.second;
    if (child->IsInThisTree(*this)) {
        std::swap(parent, child);
    }

    //std::cout << "grow " << parent->index << " " << child->index << " " << child->matched_edge->OtherEnd(*child).index << std::endl;

    if (!parent->IsInThisTree(*this)) {
        throw std::runtime_error("In Tree::Grow: parent vertex is not in this tree");
    }
    if (!parent->Plus()) {
        throw std::runtime_error("In Tree::Grow: parent is not a plus");
    }
    if (child->IsInSomeTree()) {
        throw std::runtime_error("In Tree::Grow: child vertex is not free");
    }

    child->MakeATreeChild(edge);

    if (!child->IsElementary()) {
        minus_blossoms.insert(child);
    }
}

void Tree::Shrink(EdgeWeighted &edge_plus_plus) {
    auto endpoints = edge_plus_plus.Endpoints();
    Node *first = &endpoints.first;
    Node *second = &endpoints.second;

    //std::cout << "shrink " << first->index << " " << second->index << std::endl;

    Node &lca = Node::LCA(*first, *second);

    // also need to update the minus_blossoms
    std::vector<EdgeWeighted *> blossom_edges;
    while (first != &lca) {
        blossom_edges.push_back(first->TreeParentEdge());
        first = &first->TreeParentEdge()->OtherEnd(*first);
        if (first->Minus() && !first->IsElementary()) {
            minus_blossoms.erase(first);
        }
    }
    std::reverse(blossom_edges.begin(), blossom_edges.end());
    blossom_edges.push_back(&edge_plus_plus);
    while (second != &lca) {
        blossom_edges.push_back(second->TreeParentEdge());
        second = &second->TreeParentEdge()->OtherEnd(*second);
        if (second->Minus() && !second->IsElementary()) {
            minus_blossoms.erase(second);
        }
    }

    int next_blossom_index = num_elementary_nodes;
    if (!blossom_storage->empty()) {
        next_blossom_index = blossom_storage->back().index + 1;
    }
    blossom_storage->emplace_back(blossom_edges, next_blossom_index);
    (*iter_to_self)[&blossom_storage->back()] = std::prev(blossom_storage->end());
}

void Tree::Expand(Node &blossom) {
    //std::cout << "Expand " << blossom.index << std::endl;

    if (!blossom.IsInThisTree(*this)) {
        throw std::runtime_error("In Tree::Expand: supervertex is not in this tree");
    }
    if (blossom.IsElementary()) {
        throw std::runtime_error("In Tree::Expand: supervertex must be not elementary");
    }
    if (!blossom.IsTopBlossom()) {
        throw std::runtime_error("In Tree::Expand: supervertex must be a top blossom");
    }
    if ((!blossom.Minus()) || (blossom.Plus())) {
        throw std::runtime_error("In Tree::Expand: supervertex is not a minus");
    }
    if (blossom.DualVariableQuadrupled() != 0) {
        throw std::runtime_error("In Tree::Expand: the supervertex has to have zero dual variable");
    }

    std::vector<Node *> children = blossom.BlossomChildren();

    blossom.Dissolve();

    // update minus_blossoms
    for (Node * child : children) {
        if (child->Minus() && !child->IsElementary()) {
            minus_blossoms.insert(child);
        }
    }
    minus_blossoms.erase(&blossom);

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
    if (!parent->IsInThisTree(*this)) {
        std::swap(parent, child);
    }

    if (!parent->IsInThisTree(*this)) {
        throw std::runtime_error("In Augment: parent vertex is not in this tree");
    }
    if (!child->IsInSomeTree()) {
        throw std::runtime_error("In Augment: child is not in a tree");
    }

    Node::MakeMatched(edge);
    // std::cout << "make matched " << parent->index << " " << child->index << std::endl;

    AugmentFromNode(*parent);
    Tree *other_tree = child->TreeOf();
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

        const std::vector<EdgeWeighted *> &ch = f.node->TreeChildren();
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

    if (!root->TopBlossom().IsInSomeTree()) {
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
        for (const EdgeWeighted *child : cur_vertex->TreeChildren()) {
            queue.push(&child->OtherEnd(*cur_vertex));
        }
    }

    for (Node *vertex : all_vertices) {
        vertex->ClearDuringTreeDissolve();
    }
}

void Tree::AugmentFromNode(Node &vertex) {
    // std::cout << "augment from node " << vertex.index << std::endl;

    if (!vertex.IsInThisTree(*this)) {
        throw std::runtime_error("In Tree::AugmentFromNode: vertex is not in this tree");
    }

    bool match = false;
    const std::vector<EdgeWeighted *> path = vertex.PathToRoot();
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

EdgeWeighted *Tree::GrowableEdge() const {
    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        const Node *node = queue.front();
        queue.pop();

        if (node->Plus()) {
            for (EdgeWeighted *edge : node->neighbors) {
                if ((edge->SlackQuadrupled() == 0) && (!edge->OtherEnd(*node).IsInSomeTree())) {
                    return edge;
                }
            }
        }

        for (const EdgeWeighted *child : node->TreeChildren()) {
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

        if (node->Plus()) {
            for (EdgeWeighted *edge : node->neighbors) {
                if (edge->SlackQuadrupled() == 0) {
                    Node &other_end = edge->OtherEnd(*node);
                    if ((!other_end.IsInThisTree(*this)) && (other_end.Plus()) && (other_end.IsInSomeTree())) {
                        return edge;
                    }
                }
            }
        }

        for (const EdgeWeighted *child : node->TreeChildren()) {
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

        if (node->Plus()) {
            for (EdgeWeighted *edge : node->neighbors) {
                if (edge->SlackQuadrupled() == 0) {
                    Node &other_end = edge->OtherEnd(*node);
                    if ((other_end.IsInThisTree(*this)) && (other_end.Plus())) {
                        return edge;
                    }
                }
            }
        }

        for (const EdgeWeighted *child : node->TreeChildren()) {
            queue.push(&child->OtherEnd(*node));
        }
    }

    return nullptr;
}

Node *Tree::ExpandableBlossom() const {
    // returns a blossom that is a minus and has zero dual variable
    // if there is no such blossom, returns nullptr

    if (minus_blossoms.empty()) {
        return nullptr;
    }
    if ((*minus_blossoms.begin())->DualVariableQuadrupled() == 0) {
        return *minus_blossoms.begin();
    }
    return nullptr;
}

int Tree::PlusEmptySlack() const {
    // returns the minimal quadrupled slack of a (+, empty) edge

    if (!root->TopBlossom().IsInSomeTree()) {
        throw std::runtime_error("In Tree::MakePrimalUpdates: we are in the deleted tree");
    }

    int answer = INT32_MAX;

    std::queue<Node *> queue;
    queue.push(&root->TopBlossom());
    while (!queue.empty()) {
        const Node *node = queue.front();
        queue.pop();

        if (node->Plus()) {
            for (EdgeWeighted *edge : node->neighbors) {
                if ((edge->SlackQuadrupled() < answer) && (!edge->OtherEnd(*node).IsInSomeTree())) {
                    answer = edge->SlackQuadrupled();
                }
            }
        }

        for (const EdgeWeighted *child : node->TreeChildren()) {
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

        if (node->Plus()) {
            for (EdgeWeighted *edge : node->neighbors) {
                if (edge->SlackQuadrupled() < answer) {
                    Node &other_end = edge->OtherEnd(*node);
                    if ((!other_end.IsInThisTree(*this)) && (other_end.Plus())) {
                        answer = edge->SlackQuadrupled();
                    }
                }
            }
        }

        for (const EdgeWeighted *child : node->TreeChildren()) {
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

        if (node->Plus()) {
            for (EdgeWeighted *edge : node->neighbors) {
                if (edge->SlackQuadrupled() < answer) {
                    Node &other_end = edge->OtherEnd(*node);
                    if ((other_end.IsInThisTree(*this)) && (other_end.Plus())) {
                        answer = edge->SlackQuadrupled();
                    }
                }
            }
        }

        for (const EdgeWeighted *child : node->TreeChildren()) {
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

        if (node->Plus()) {
            for (EdgeWeighted *edge : node->neighbors) {
                if (edge->SlackQuadrupled() < answer) {
                    Node &other_end = edge->OtherEnd(*node);
                    if ((!other_end.IsInThisTree(*this)) && (other_end.Minus())) {
                        answer = edge->SlackQuadrupled();
                    }
                }
            }
        }

        for (const EdgeWeighted *child : node->TreeChildren()) {
            queue.push(&child->OtherEnd(*node));
        }
    }

    return answer;
}

int Tree::MinMinusBlossomVariable() const {
    // returns the minimum quadrupled dual variable of a (-) non-elementary blossom

    if (minus_blossoms.empty()) {
        return INT32_MAX;
    }
    return (*minus_blossoms.begin())->DualVariableQuadrupled();
}
