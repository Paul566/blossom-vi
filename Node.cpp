#include "Node.h"

#include <functional>
#include <unordered_set>

#include "EdgeWeighted.h"
#include "Tree.h"

Node::Node(const int index_) : index(index_) {
    dual_var_quadrupled_amortized = 0;
    blossom_parent = nullptr;
    matched_edge = nullptr;
    tree_parent = nullptr;
    tree = nullptr;
    plus = false;
    blossom_brother_clockwise = nullptr;
    blossom_brother_anticlockwise = nullptr;
}

Node::Node(const std::vector<EdgeWeighted *> &blossom_edges, int index_) : index(index_) {
    // blossom_edges need to be consecutive and start and end at the receptacle

    if (blossom_edges.size() % 2 == 0) {
        throw std::runtime_error("In Node: blossom has to have an odd number of vertices");
    }

    blossom_parent = nullptr;
    plus = true; // after Shrink, we must be a plus

    // update blossom_children
    blossom_children.reserve(blossom_edges.size());
    Node *cur_vertex = &SharedNode(*blossom_edges.front(), *blossom_edges.back());
    for (const EdgeWeighted *edge : blossom_edges) {
        cur_vertex = &edge->OtherEnd(*cur_vertex);
        blossom_children.push_back(cur_vertex);
    }

    // initialize blossom brothers of the children
    for (EdgeWeighted *edge : blossom_edges) {
        cur_vertex->blossom_brother_clockwise = edge;
        cur_vertex = &edge->OtherEnd(*cur_vertex);
        cur_vertex->blossom_brother_anticlockwise = edge;
    }

    // update blossom_parent and dual_var_quadrupled_amortized of the children
    for (Node *vertex : blossom_children) {
        if (vertex->blossom_parent) {
            throw std::runtime_error("In Node: some child node already has a parent");
        }
        vertex->blossom_parent = this;
        if (vertex->plus) {
            vertex->dual_var_quadrupled_amortized += vertex->tree->dual_var_quadrupled;
        } else {
            vertex->dual_var_quadrupled_amortized -= vertex->tree->dual_var_quadrupled;
        }
    }

    std::unordered_set<Node *> blossom_set;
    blossom_set.reserve(blossom_children.size());
    for (Node *vertex : blossom_children) {
        blossom_set.insert(vertex);
    }

    // find tree_children
    for (const auto &vertex : blossom_children) {
        for (auto &child_edge : vertex->tree_children) {
            if (!blossom_set.contains(&child_edge->OtherEnd(*vertex))) {
                tree_children.push_back(child_edge);
            }
        }
    }

    // find neighbors + update edges
    for (const Node *vertex : blossom_children) {
        for (EdgeWeighted *edge : vertex->neighbors) {
            if (!blossom_set.contains(&edge->OtherEnd(*vertex))) {
                neighbors.push_back(edge);
                edge->UpdateAfterShrink(*vertex);
            }
        }
    }

    // find matched_edge and tree_parent
    matched_edge = nullptr;
    tree_parent = nullptr;
    int matched_edge_cnt_debug = 0;
    for (EdgeWeighted *edge : neighbors) {
        if (edge->matched) {
            matched_edge = edge;
            tree_parent = edge;
            ++matched_edge_cnt_debug;
        }
    }
    if (matched_edge_cnt_debug > 1) {
        throw std::runtime_error("In Node: found several adjacent matched edges");
    }

    // find tree_root
    tree = blossom_children.front()->tree;
    for (const auto &vertex : blossom_children) {
        if (vertex->tree != tree) {
            throw std::runtime_error("In Node: not all the vertices have the same tree_root");
        }
    }

    dual_var_quadrupled_amortized = -tree->dual_var_quadrupled;
    // the true dual variable must be zero
}

Node &Node::TopBlossom() {
    Node *cur_vertex = this;
    while (cur_vertex->blossom_parent) {
        cur_vertex = cur_vertex->blossom_parent;
    }
    return *cur_vertex;
}

bool Node::IsInSomeTree() const {
    return tree != nullptr;
}

bool Node::IsInThisTree(const Tree &tree_) const {
    return tree == &tree_;
}

Tree *Node::TreeOf() const {
    return tree;
}

bool Node::Plus() const {
    if (plus && !IsInSomeTree()) {
        throw std::runtime_error("In Node::Plus: incorrect state discovered");
    }
    return plus;
}

bool Node::Minus() const {
    return !plus && (tree != nullptr);
}

bool Node::IsMatched() const {
    return matched_edge != nullptr;
}

Node * Node::MatchedNeighbor() const {
    if (matched_edge) {
        return &matched_edge->OtherEnd(*this);
    }
    return nullptr;
}

EdgeWeighted *Node::TreeParentEdge() const {
    return tree_parent;
}

const std::vector<EdgeWeighted *> &Node::TreeChildren() const {
    return tree_children;
}

Node *Node::BlossomParent() const {
    return blossom_parent;
}

const std::vector<Node *> &Node::BlossomChildren() const {
    return blossom_children;
}

std::vector<EdgeWeighted *> Node::PathToRoot() const {
    if (!IsTopBlossom()) {
        throw std::runtime_error("In Node::PathToRoot: vertex is not a top blossom");
    }
    if (!IsInSomeTree()) {
        throw std::runtime_error("In Node::PathToRoot: vertex is not in a tree");
    }

    const Node *cur_vertex = this;
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

Node &Node::LCA(Node &first, Node &second) {
    // TODO can be made better

    if (first.tree != second.tree) {
        throw std::runtime_error("In LCA: first and second are in different trees");
    }
    if (!first.IsInSomeTree()) {
        throw std::runtime_error("In LCA: node is not in any tree");
    }

    Node &root_blossom = first.tree->root->TopBlossom();
    Node *first_ptr = &first;
    Node *second_ptr = &second;

    std::unordered_set<Node *> visited;
    visited.insert(first_ptr);
    while (first_ptr != &root_blossom) {
        first_ptr = &first_ptr->tree_parent->OtherEnd(*first_ptr);
        visited.insert(first_ptr);
    }

    while (second_ptr != &root_blossom) {
        if (visited.contains(second_ptr)) {
            return *second_ptr;
        }
        second_ptr = &second_ptr->tree_parent->OtherEnd(*second_ptr);
    }

    return root_blossom;
}

void Node::Dissolve() {
    if (blossom_children.empty()) {
        throw std::runtime_error("In Dissolve: can't dissolve an elementary vertex");
    }
    if (!matched_edge) {
        throw std::runtime_error("In Dissolve: there is no matched_edge -- this blossom can't contain the root");
    }
    if (blossom_parent) {
        throw std::runtime_error("In Dissolve: this must be a top blossom");
    }

    RotateReceptacle(&matched_edge->DeeperNode(*this));

    if (tree) {
        UpdateInternalTreeStructure();
        for (Node *child : blossom_children) {
            if (child->plus) {
                child->dual_var_quadrupled_amortized -= tree->dual_var_quadrupled;
            }
            if (child->Minus()) {
                child->dual_var_quadrupled_amortized += tree->dual_var_quadrupled;
            }
        }
    } else {
        ClearInternalTreeStructure();
    }

    for (Node *child : blossom_children) {
        child->blossom_brother_clockwise = nullptr;
        child->blossom_brother_anticlockwise = nullptr;
        child->blossom_parent = nullptr;
    }

    for (EdgeWeighted *edge : neighbors) {
        edge->UpdateAfterDissolve(*this);
    }
}

void Node::RotateReceptacle(Node *new_receptacle) const {
    // update the structure of matched edges inside the blossom and at the new receptacle

    if (blossom_children.empty()) {
        throw std::runtime_error("In RotateReceptacle: the vertex must be a blossom");
    }
    if (new_receptacle->blossom_parent != this) {
        throw std::runtime_error("In RotateReceptacle: the new receptacle is not a blossom child of this node");
    }

    Node *cur_vertex = &new_receptacle->blossom_brother_clockwise->OtherEnd(*new_receptacle);
    bool match = false;
    while (cur_vertex != new_receptacle) {
        if (match) {
            MakeMatched(*cur_vertex->blossom_brother_anticlockwise);
        } else {
            MakeUnmatched(*cur_vertex->blossom_brother_anticlockwise);
        }
        match = !match;
        cur_vertex = &cur_vertex->blossom_brother_clockwise->OtherEnd(*cur_vertex);
    }
    MakeUnmatched(*new_receptacle->blossom_brother_anticlockwise);

    new_receptacle->matched_edge = matched_edge; // never a nullptr
}

void Node::PrintNode() const {
    if (IsElementary()) {
        std::cout << index << ": y_v = " << DualVariableQuadrupled() / 4.;
        if (!IsTopBlossom()) {
            std::cout << " blossom_parent: " << blossom_parent->index << std::endl;
            return;
        }
        for (const EdgeWeighted *edge : neighbors) {
            std::cout << " (" << edge->OtherEnd(*this).index << ", " << edge->weight << ", " << edge->
                SlackQuadrupled() / 4. << ", " << edge->matched << ")";
        }
        std::cout << std::endl;
        return;
    }

    std::cout << index << ": y_v = " << DualVariableQuadrupled() / 4. << ", blossom_children: ";
    for (const Node *child : blossom_children) {
        std::cout << child->index << " ";
    }

    if (!IsTopBlossom()) {
        std::cout << "blossom_parent: " << blossom_parent->index << std::endl;
        return;
    }

    for (const EdgeWeighted *edge : neighbors) {
        std::cout << "(" << edge->OtherEnd(*this).index << ", " << edge->weight << ", " << edge->
            SlackQuadrupled() / 4. << ", " << edge->matched << ") ";
    }
    std::cout << std::endl;
}

int Node::DualVariableQuadrupled() const {
    if (tree == nullptr) {
        return dual_var_quadrupled_amortized;
    }
    if (blossom_parent) {
        return dual_var_quadrupled_amortized;
    }
    if (plus) {
        return dual_var_quadrupled_amortized + tree->dual_var_quadrupled;
    }
    return dual_var_quadrupled_amortized - tree->dual_var_quadrupled;
}

int Node::DualVariableQuadrupledAmortized() const {
    return dual_var_quadrupled_amortized;
}

bool Node::IsElementary() const {
    return blossom_children.empty();
}

bool Node::IsTopBlossom() const {
    return blossom_parent == nullptr;
}

void Node::ClearDuringTreeDissolve() {
    if (blossom_parent) {
        throw std::runtime_error("ClearDuringTreeDissolve can only be called for top blossoms");
    }
    if (plus) {
        dual_var_quadrupled_amortized += tree->dual_var_quadrupled;
    } else {
        dual_var_quadrupled_amortized -= tree->dual_var_quadrupled;
    }
    tree = nullptr;
    tree_children.clear();
    tree_parent = nullptr;
    plus = false;
}

void Node::MakeRootOfTree(Tree &tree_) {
    plus = true;
    tree = &tree_;
}

void Node::MakeATreeChild(EdgeWeighted &edge_to_parent) {
    Node &parent = edge_to_parent.OtherEnd(*this);
    parent.tree_children.push_back(&edge_to_parent);

    tree_parent = &edge_to_parent;
    tree = parent.tree;
    plus = false;
    dual_var_quadrupled_amortized += tree->dual_var_quadrupled;
    tree_children = {matched_edge};

    Node &grandchild = matched_edge->OtherEnd(*this);
    grandchild.tree_parent = matched_edge;
    grandchild.tree = tree;
    grandchild.plus = true;
    grandchild.dual_var_quadrupled_amortized -= tree->dual_var_quadrupled;
}

void Node::MakeSlackNonnegativeInInit() {
    // is called once for each vertex in GreedyInit

    if (neighbors.empty()) {
        std::cout << "vertex " << index << ": no neighbors" << std::endl;
        throw std::runtime_error("Found an isolated vertex => no perfect matching exists");
    }

    int min_weight = neighbors.front()->weight;
    for (const EdgeWeighted *edge : neighbors) {
        if (edge->weight < min_weight) {
            min_weight = edge->weight;
        }
    }

    dual_var_quadrupled_amortized += min_weight * 2;
}

void Node::InitVarGreedily() {
    // is called once for each vertex in GreedyInit

    if (neighbors.empty()) {
        throw std::runtime_error("Found an isolated vertex => no perfect matching exists");
    }

    if (IsMatched()) {
        return;
    }

    EdgeWeighted *smallest_slack_edge = neighbors.front();
    for (EdgeWeighted *edge : neighbors) {
        if (edge->SlackQuadrupled() < smallest_slack_edge->SlackQuadrupled()) {
            smallest_slack_edge = edge;
        }
    }

    const int smallest_slack_quadrupled = smallest_slack_edge->SlackQuadrupled();
    dual_var_quadrupled_amortized += smallest_slack_quadrupled;

    if (!smallest_slack_edge->OtherEnd(*this).IsMatched()) {
        // if the other vertex is also unmatched, match the edge
        MakeMatched(*smallest_slack_edge);
    }
}

Node &Node::SharedNode(const EdgeWeighted &first_edge, const EdgeWeighted &second_edge) {
    // returns the max blossom adjacent to both edges

    auto [first_head, first_tail] = first_edge.Endpoints();
    auto [second_head, second_tail] = second_edge.Endpoints();
    if ((&first_head == &second_head) || (&first_head == &second_tail)) {
        return first_head;
    }
    if ((&first_tail == &second_head) || (&first_tail == &second_tail)) {
        return first_tail;
    }
    throw std::runtime_error("In Node::SharedNode: the edges don't share an endpoint");
}

void Node::UpdateInternalTreeStructure() {
    // updates tree_children, tree_parents, tree_root for the vertices inside this blossom before the Expand
    // updates the plus and minus markers
    if (!matched_edge) {
        throw std::runtime_error("In UpdateInternalTreeStructure: there is no matched_edge");
    }
    if (!tree_parent) {
        throw std::runtime_error("In UpdateInternalTreeStructure: there is no tree_parent");
    }
    Node &elder_child = tree_parent->DeeperNode(*this);
    const Node &receptacle = matched_edge->DeeperNode(*this);

    std::function next_edge = [](const Node *current) -> EdgeWeighted * {
        return current->blossom_brother_clockwise;
    };
    if (elder_child.blossom_brother_anticlockwise->matched) {
        next_edge = [](const Node *current) -> EdgeWeighted * {
            return current->blossom_brother_anticlockwise;
        };
    }

    // the part that stays in the tree
    Node *cur_vertex = &elder_child;
    EdgeWeighted *prev_edge = tree_parent;
    bool is_plus = false;
    while (cur_vertex != &receptacle) {
        cur_vertex->tree_parent = prev_edge;
        cur_vertex->tree = tree;
        cur_vertex->tree_children = {next_edge(cur_vertex)};
        cur_vertex->plus = is_plus;

        prev_edge = next_edge(cur_vertex);
        cur_vertex = &prev_edge->OtherEnd(*cur_vertex);
        is_plus = !is_plus;
    }
    cur_vertex->tree_parent = prev_edge;
    cur_vertex->tree = tree;
    cur_vertex->tree_children = {matched_edge};
    cur_vertex->plus = false;

    // the part that goes to waste
    cur_vertex = &next_edge(cur_vertex)->OtherEnd(*cur_vertex);
    while (cur_vertex != &elder_child) {
        cur_vertex->tree_parent = nullptr;
        cur_vertex->tree = nullptr;
        cur_vertex->tree_children.clear();
        cur_vertex->plus = false;

        cur_vertex = &next_edge(cur_vertex)->OtherEnd(*cur_vertex);
    }
}

void Node::ClearInternalTreeStructure() const {
    for (Node *child_blossom : blossom_children) {
        child_blossom->tree_parent = nullptr;
        child_blossom->tree = nullptr;
        child_blossom->tree_children.clear();
        child_blossom->plus = false;
    }
}

void Node::MakeMatched(EdgeWeighted &edge) {
    auto [head, tail] = edge.Endpoints();
    head.matched_edge = &edge;
    tail.matched_edge = &edge;
    edge.matched = true;
}

void Node::MakeUnmatched(EdgeWeighted &edge) {
    edge.matched = false;
}
