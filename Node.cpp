#include "Node.h"

#include <functional>
#include <unordered_set>

#include "EdgeWeighted.h"

Node::Node(const int index_) : index(index_) {
    dual_variable_quadrupled = 0;
    blossom_parent = nullptr;
    matched_edge = nullptr;
    tree_parent = nullptr;
    tree_root = nullptr;
    plus = false;
    minus = false;
    blossom_brother_clockwise = nullptr;
    blossom_brother_anticlockwise = nullptr;
}

Node::Node(const std::vector<EdgeWeighted *> &blossom_edges, int index_) : index(index_) {
    // blossom_edges need to be consecutive and start and end at the receptacle

    if (blossom_edges.size() % 2 == 0) {
        throw std::runtime_error("In Node: blossom has to have an odd number of vertices");
    }

    dual_variable_quadrupled = 0;
    blossom_parent = nullptr;
    plus = true; // after Shrink, we must be a plus
    minus = false;

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

    // update the blossom_parent of the children
    for (Node *vertex : blossom_children) {
        if (vertex->blossom_parent) {
            throw std::runtime_error("In Node: some child node already has a parent");
        }
        vertex->blossom_parent = this;
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
    tree_root = blossom_children.front()->tree_root;
    for (const auto &vertex : blossom_children) {
        if (vertex->tree_root != tree_root) {
            throw std::runtime_error("In Node: not all the vertices have the same tree_root");
        }
    }
}

Node &Node::TopBlossom() {
    Node *cur_vertex = this;
    while (cur_vertex->blossom_parent) {
        cur_vertex = cur_vertex->blossom_parent;
    }
    return *cur_vertex;
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

    if (tree_root) {
        UpdateInternalTreeStructure();
    } else {
        ClearInternalTreeStructure();
    }

    for (Node *child : blossom_children) {
        child->blossom_brother_clockwise = nullptr;
        child->blossom_brother_anticlockwise = nullptr;
        child->blossom_parent = nullptr;
    }

    for (EdgeWeighted * edge : neighbors) {
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

    new_receptacle->matched_edge = matched_edge;    // never a nullptr
}

int Node::DualVariableQuadrupled() const {
    return dual_variable_quadrupled;
}

void Node::IncreaseDualVariableQuadrupled(const int increment) {
    if (blossom_parent) {
        throw std::runtime_error("In IncreaseDualVariableQuadrupled: the vertex is not a top blossom");
    }

    dual_variable_quadrupled += increment;

    for (EdgeWeighted *edge : neighbors) {
        edge->slack_quadrupled -= increment;
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

    std::function next_edge = [] (const Node * current) -> EdgeWeighted * {
        return current->blossom_brother_clockwise;
    };
    if (elder_child.blossom_brother_anticlockwise->matched) {
        next_edge = [] (const Node * current) -> EdgeWeighted * {
            return current->blossom_brother_anticlockwise;
        };
    }

    // the part that stays in the tree
    Node *cur_vertex = &elder_child;
    EdgeWeighted *prev_edge = tree_parent;
    bool is_plus = false;
    while (cur_vertex != &receptacle) {
        cur_vertex->tree_parent = prev_edge;
        cur_vertex->tree_root = tree_root;
        cur_vertex->tree_children = {next_edge(cur_vertex)};
        cur_vertex->plus = is_plus;
        cur_vertex->minus = !is_plus;

        prev_edge = next_edge(cur_vertex);
        cur_vertex = &prev_edge->OtherEnd(*cur_vertex);
        is_plus = !is_plus;
    }
    cur_vertex->tree_parent = prev_edge;
    cur_vertex->tree_root = tree_root;
    cur_vertex->tree_children = {matched_edge};
    cur_vertex->plus = false;
    cur_vertex->minus = true;

    // the part that goes to waste
    cur_vertex = &next_edge(cur_vertex)->OtherEnd(*cur_vertex);
    while (cur_vertex != &elder_child) {
        cur_vertex->tree_parent = nullptr;
        cur_vertex->tree_root = nullptr;
        cur_vertex->tree_children.clear();
        cur_vertex->plus = false;
        cur_vertex->minus = false;

        cur_vertex = &next_edge(cur_vertex)->OtherEnd(*cur_vertex);
    }

}

void Node::ClearInternalTreeStructure() const {
    for (Node * child_blossom : blossom_children) {
        child_blossom->tree_parent = nullptr;
        child_blossom->tree_root = nullptr;
        child_blossom->tree_children.clear();
        child_blossom->plus = false;
        child_blossom->minus = false;
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
