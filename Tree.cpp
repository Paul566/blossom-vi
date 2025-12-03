#include <queue>
#include <unordered_set>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <list>

#include "Tree.h"
#include "EdgeWeighted.h"
#include "Node.h"

bool NodeComparator::operator()(const Node *a, const Node *b) const {
    int a_var = a->DualVariableQuadrupled();
    int b_var = b->DualVariableQuadrupled();
    if (a_var != b_var) {
        return a_var < b_var;
    }
    return a < b;
}

bool EdgeComparator::operator()(const EdgeWeighted *a, const EdgeWeighted *b) const {
    int a_key = a->SlackQuadrupled();
    int b_key = b->SlackQuadrupled();
    if (a_key != b_key) {
        return a_key < b_key;
    }
    return a < b;
}

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
    UpdateQueuesAfterInit();
}

void Tree::Grow(EdgeWeighted &edge) {
    auto endpoints = edge.Endpoints();
    Node *parent = &endpoints.first;
    Node *child = &endpoints.second;
    if (child->IsInThisTree(*this)) {
        std::swap(parent, child);
    }

    if (!parent->IsInThisTree(*this)) {
        throw std::runtime_error("In Tree::Grow: parent vertex is not in this tree");
    }
    if (!parent->Plus()) {
        throw std::runtime_error("In Tree::Grow: parent is not a plus");
    }
    if (child->IsInSomeTree()) {
        throw std::runtime_error("In Tree::Grow: child vertex is not free");
    }

    // std::cout << "GROW " << parent->index << " " << child->index << " " << child->MatchedNeighbor()->index << std::endl;

    child->MakeATreeChild(edge);
    UpdateQueuesAfterGrow(*child);
}

void Tree::Shrink(EdgeWeighted &edge_plus_plus) {
    auto endpoints = edge_plus_plus.Endpoints();
    Node *first = &endpoints.first;
    Node *second = &endpoints.second;

    // std::cout << "SHRINK " << first->index << " " << second->index << std::endl;

    Node &lca = Node::LCA(*first, *second);

    std::vector<EdgeWeighted *> blossom_edges;
    while (first != &lca) {
        blossom_edges.push_back(first->TreeParentEdge());
        first = &first->TreeParentEdge()->OtherEnd(*first);
    }
    std::reverse(blossom_edges.begin(), blossom_edges.end());
    blossom_edges.push_back(&edge_plus_plus);
    while (second != &lca) {
        blossom_edges.push_back(second->TreeParentEdge());
        second = &second->TreeParentEdge()->OtherEnd(*second);
    }

    int next_blossom_index = num_elementary_nodes;
    if (!blossom_storage->empty()) {
        next_blossom_index = blossom_storage->back().index + 1;
    }
    blossom_storage->emplace_back(blossom_edges, next_blossom_index);
    Node &new_blossom = blossom_storage->back();
    (*iter_to_self)[&new_blossom] = std::prev(blossom_storage->end());

    UpdateQueuesAfterShrink(new_blossom);
}

void Tree::Expand(Node &blossom) {
    // std::cout << "EXPAND " << blossom.index << std::endl;

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
    minus_blossoms.erase(&blossom);
    UpdateQueuesAfterExpand(children);

    std::list<Node>::iterator it = (*iter_to_self)[&blossom];
    iter_to_self->erase(&blossom);
    blossom_storage->erase(it);
}

Tree *Tree::Augment(EdgeWeighted &edge) {
    // TODO amortize tree dissolves

    auto endpoints = edge.Endpoints();
    Node *parent = &endpoints.first;
    Node *child = &endpoints.second;
    if (!parent->IsInThisTree(*this)) {
        std::swap(parent, child);
    }

    // std::cout << "AUGMENT " << parent->index << " " << child->index << std::endl;

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

void Tree::PrintTree() const {
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

void Tree::DissolveTree() const {
    // clears the tree structure from the nodes
    // updates priority queues of the adjacent trees

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

    // update plus_empty_edges of other trees
    for (Node *vertex : all_vertices) {
        for (EdgeWeighted *neighbor_edge : vertex->neighbors) {
            Node &neighbor = neighbor_edge->OtherEnd(*vertex);
            if (neighbor.Plus()) {
                neighbor.TreeOf()->plus_empty_edges.insert(neighbor_edge);
                if (vertex->Plus()) {
                    neighbor.TreeOf()->plus_plus_external_edges.erase(neighbor_edge);
                }
            }
        }
    }

    for (Node *vertex : all_vertices) {
        vertex->ClearDuringTreeDissolve();
    }
}

void Tree::UpdateQueuesAfterGrow(Node &child) {
    if (!child.IsElementary()) {
        minus_blossoms.insert(&child);
    }

    // neighboring edges of child are no longer (+, empty)
    for (EdgeWeighted *neighbor_edge : child.neighbors) {
        Node &neighbor = neighbor_edge->OtherEnd(child);
        if (neighbor.Plus()) {
            neighbor.TreeOf()->plus_empty_edges.erase(neighbor_edge);
        }
    }

    // neighboring edges of grandchild can be (+, empty)
    Node &grandchild = child.TreeChildren().front()->OtherEnd(child);
    for (EdgeWeighted *edge_to_neighbor : grandchild.neighbors) {
        Node &grandchild_neighbor = edge_to_neighbor->OtherEnd(grandchild);
        if (!grandchild_neighbor.IsInSomeTree()) {
            plus_empty_edges.insert(edge_to_neighbor);
        }
    }

    // neighboring edges of grandchild are no longer (+, empty)
    for (EdgeWeighted *neighbor_edge : grandchild.neighbors) {
        Node &neighbor = neighbor_edge->OtherEnd(grandchild);
        if (neighbor.Plus()) {
            neighbor.TreeOf()->plus_empty_edges.erase(neighbor_edge);
        }
    }

    // maybe new (+, +) internal edges adjacent to grandchild
    for (EdgeWeighted *edge_to_neighbor : grandchild.neighbors) {
        Node &grandchild_neighbor = edge_to_neighbor->OtherEnd(grandchild);
        if (grandchild_neighbor.IsInThisTree(*this) && grandchild_neighbor.Plus()) {
            plus_plus_internal_edges.insert(edge_to_neighbor);
        }
    }

    // maybe new (+, +) external edges adjacent to grandchild
    for (EdgeWeighted *edge_to_neighbor : grandchild.neighbors) {
        Node &grandchild_neighbor = edge_to_neighbor->OtherEnd(grandchild);
        if (!grandchild_neighbor.IsInThisTree(*this) && grandchild_neighbor.Plus()) {
            plus_plus_external_edges.insert(edge_to_neighbor);
            grandchild_neighbor.TreeOf()->plus_plus_external_edges.insert(edge_to_neighbor);
        }
    }
}

void Tree::UpdateQueuesAfterShrink(const Node &blossom) {
    // update minus_blossoms
    for (Node *child : blossom.BlossomChildren()) {
        if (child->Minus()) {
            minus_blossoms.erase(child);
        }
    }

    // update plus_empty edges
    for (EdgeWeighted *edge : blossom.neighbors) {
        if (!edge->DeeperNode(blossom).Plus()) {
            if (!edge->OtherEnd(blossom).IsInSomeTree()) {
                plus_empty_edges.insert(edge);
            }
        }
    }

    // maybe new (+, +) internal edges
    for (EdgeWeighted *edge : blossom.neighbors) {
        Node &neighbor = edge->OtherEnd(blossom);
        if (neighbor.IsInThisTree(*this) && neighbor.Plus() && !edge->DeeperNode(blossom).Plus()) {
            plus_plus_internal_edges.insert(edge);
        }
    }

    // erase the (+, +) edges that were inside the blossom
    for (Node *child : blossom.BlossomChildren()) {
        if (child->Plus()) {
            for (EdgeWeighted *neighbor_edge : child->neighbors) {
                if (neighbor_edge->IsInsideBlossom()) {
                    if (neighbor_edge->OtherEnd(*child).Plus()) {
                        plus_plus_internal_edges.erase(neighbor_edge);
                    }
                }
            }
        }
    }

    // maybe new (+, +) external edges
    for (EdgeWeighted *edge : blossom.neighbors) {
        Node &neighbor = edge->OtherEnd(blossom);
        if (!neighbor.IsInThisTree(*this) && neighbor.Plus() && !edge->DeeperNode(blossom).Plus()) {
            plus_plus_external_edges.insert(edge);
            neighbor.TreeOf()->plus_plus_external_edges.insert(edge);
        }
    }
}

void Tree::UpdateQueuesAfterExpand(const std::vector<Node *> &children) {
    // update minus_blossoms
    for (Node *child : children) {
        if (child->Minus() && !child->IsElementary()) {
            minus_blossoms.insert(child);
        }
    }

    // update plus_empty for the part that stays in the tree
    for (Node *child : children) {
        if (child->Plus()) {
            for (EdgeWeighted *edge : child->neighbors) {
                if (!edge->OtherEnd(*child).IsInSomeTree()) {
                    plus_empty_edges.insert(edge);
                }
            }
        }
    }

    // update plus_empty_edges of the trees adjacent to the part that goes to waste
    for (Node *child : children) {
        if (!child->IsInSomeTree()) {
            for (EdgeWeighted *edge : child->neighbors) {
                Node &other_end = edge->OtherEnd(*child);
                if (other_end.Plus()) {
                    other_end.TreeOf()->plus_empty_edges.insert(edge);
                }
            }
        }
    }

    // update (+, +) internal
    for (Node *child : children) {
        if (child->Plus()) {
            for (EdgeWeighted *edge : child->neighbors) {
                Node &other_end = edge->OtherEnd(*child);
                if (other_end.IsInThisTree(*this) && other_end.Plus()) {
                    plus_plus_internal_edges.insert(edge);
                }
            }
        }
    }

    // update (+, +) external
    for (Node *child : children) {
        if (child->Plus()) {
            for (EdgeWeighted *edge : child->neighbors) {
                Node &other_end = edge->OtherEnd(*child);
                if (!other_end.IsInThisTree(*this) && other_end.Plus()) {
                    plus_plus_external_edges.insert(edge);
                    other_end.TreeOf()->plus_plus_external_edges.insert(edge);
                }
            }
        }
    }
}

void Tree::UpdateQueuesAfterInit() {
    for (EdgeWeighted *edge : root->neighbors) {
        if (!edge->OtherEnd(*root).IsInSomeTree()) {
            plus_empty_edges.insert(edge);
        }
    }

    for (EdgeWeighted *edge : root->neighbors) {
        Node &neighbor = edge->OtherEnd(*root);
        if (neighbor.Plus()) {
            neighbor.TreeOf()->plus_empty_edges.erase(edge);
            plus_plus_external_edges.insert(edge);
            neighbor.TreeOf()->plus_plus_external_edges.insert(edge);
        }
    }
}

void Tree::ValidatePlusEmpty() const {
    // validates that all the edges in plus_empty_edges are really (+, empty)
    for (EdgeWeighted *edge : plus_empty_edges) {
        auto endpoints = edge->Endpoints();
        Node *first = &endpoints.first;
        Node *second = &endpoints.second;
        if (second->IsInThisTree(*this)) {
            std::swap(first, second);
        }

        if (!first->IsInThisTree(*this) || !first->Plus() || !first->IsTopBlossom()) {
            throw std::runtime_error("Invalid plus_empty_edges");
        }
        if (second->IsInSomeTree()) {
            throw std::runtime_error("Invalid plus_empty_edges");
        }
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

EdgeWeighted *Tree::GrowableEdge() {
    if (plus_empty_edges.empty()) {
        return nullptr;
    }
    if ((*plus_empty_edges.begin())->SlackQuadrupled() == 0) {
        return *plus_empty_edges.begin();
    }
    return nullptr;
}

EdgeWeighted *Tree::AugmentableEdge() const {
    if (plus_plus_external_edges.empty()) {
        return nullptr;
    }
    if ((*plus_plus_external_edges.begin())->SlackQuadrupled() == 0) {
        return *plus_plus_external_edges.begin();
    }
    return nullptr;
}

EdgeWeighted *Tree::ShrinkableEdge() const {
    if (plus_plus_internal_edges.empty()) {
        return nullptr;
    }
    if ((*plus_plus_internal_edges.begin())->SlackQuadrupled() == 0) {
        return *plus_plus_internal_edges.begin();
    }
    return nullptr;
}

Node *Tree::ExpandableBlossom() {
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

int Tree::PlusEmptySlack() {
    // returns the minimal quadrupled slack of a (+, empty) edge

    if (plus_empty_edges.empty()) {
        return INT32_MAX;
    }
    return (*plus_empty_edges.begin())->SlackQuadrupled();
}

int Tree::PlusPlusExternalSlack() const {
    // returns the minimal quadrupled slack of a (+, +) edge, where the other vertex is in another tree

    if (plus_plus_external_edges.empty()) {
        return INT32_MAX;
    }
    return (*plus_plus_external_edges.begin())->SlackQuadrupled();
}

int Tree::PlusPlusInternalSlack() const {
    // returns the minimal quadrupled slack of a (+, +) edge, where the other vertex is in the same tree

    if (plus_plus_internal_edges.empty()) {
        return INT32_MAX;
    }
    return (*plus_plus_internal_edges.begin())->SlackQuadrupled();
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

int Tree::MinMinusBlossomVariable() {
    // returns the minimum quadrupled dual variable of a (-) non-elementary blossom

    if (minus_blossoms.empty()) {
        return INT32_MAX;
    }
    return (*minus_blossoms.begin())->DualVariableQuadrupled();
}
