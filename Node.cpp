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
    blossom_brother_clockwise = nullptr;
    blossom_brother_anticlockwise = nullptr;
}

Node::Node(std::vector<std::shared_ptr<EdgeWeighted>> blossom_edges) : index(-1) {
    // blossom_edges need to be consecutive and start and end at the receptacle

    if (blossom_edges.size() % 2 == 0) {
        throw std::runtime_error("In Node: blossom has to have an odd number of vertices");
    }

    dual_variable_quadrupled = 0;
    parent_blossom = nullptr;
    plus = true;    // after Shrink, we must be a plus
    minus = false;

    // find minus_parents for the children
    std::vector<std::shared_ptr<Node>> blossom_vertices;
    blossom_vertices.reserve(blossom_edges.size());
    const auto receptacle = SharedMaxDistinctBlossom(blossom_edges.front(), blossom_edges.back());
    auto cur_vertex = receptacle;
    if (cur_vertex == nullptr) {
        throw std::runtime_error("In Node: first and last edges don't share a top node");
    }
    for (const auto& edge : blossom_edges) {
        cur_vertex = edge->OtherBlossom(cur_vertex);
        blossom_vertices.push_back(cur_vertex);
    }

    children_blossom = blossom_vertices;

    // initialize blossom brothers of the children
    cur_vertex = receptacle;
    for (const auto& edge : blossom_edges) {
        cur_vertex->blossom_brother_clockwise = edge;
        cur_vertex = edge->OtherBlossom(cur_vertex);
        cur_vertex->blossom_brother_anticlockwise = edge;
    }

    // update the parent of the children
    for (const auto& vertex : blossom_vertices) {
        if (vertex->parent_blossom) {
            throw std::runtime_error("In Node: some vertex is not a top blossom");
        }
        vertex->parent_blossom = shared_from_this();
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
    if (!parent_blossom) {
        return shared_from_this();
    }
    return parent_blossom->TopBlossom();
}

std::shared_ptr<Node> Node::NextToTopBlossom() {
    if (!parent_blossom) {
        throw std::runtime_error("In NextToTopBlossom: this must not be a top blossom");
    }

    auto cur_vertex = shared_from_this();
    while (cur_vertex->parent_blossom->parent_blossom) {
        cur_vertex = cur_vertex->parent_blossom;
    }
    return cur_vertex;
}

void Node::Dissolve() {
    if (index != -1) {
        throw std::runtime_error("In Dissolve: can't dissolve an elementary vertex");
    }
    if (!matched_edge) {
        throw std::runtime_error("In Dissolve: there is no matched_edge");
    }
    if (parent_blossom) {
        throw std::runtime_error("In Dissolve: this must be a top blossom");
    }

    MakeChildrenTopBlossoms();

    RotateReceptacle(ReceptacleChild());

    if (!tree_root) { // if we are not in any tree
        for (const auto& child : children_blossom) {
            child->tree_parent = nullptr;
            child->tree_root = nullptr;
            child->tree_children.clear();
            child->plus = false;
            child->minus = false;
        }
        MakeChildrenBrotherless();
        return;
    }

    ClearChildrenTreeChildren();

    UpdateInternalTreeStructure();

    UpdateChildrenExternalTreeChildren();

    UpdateElderChildTreeParent();

    MakeChildrenBrotherless();
}

void Node::RotateReceptacle(const std::shared_ptr<Node> &new_receptacle) const {
    if (index != -1) {
        throw std::runtime_error("In RotateReceptacle: the vertex must be a blossom");
    }

    auto cur_vertex = new_receptacle->blossom_brother_clockwise->OtherMaxDistinctBlossom(new_receptacle);
    bool match = false;
    while (cur_vertex != new_receptacle) {
        if (match) {
            cur_vertex->blossom_brother_anticlockwise->MakeMatched();
        } else {
            cur_vertex->blossom_brother_anticlockwise->MakeUnmatched();
        }
        match = !match;
    }
    new_receptacle->blossom_brother_anticlockwise->MakeUnmatched();
}

bool Node::Ancestor(const std::shared_ptr<Node> &potential_ancestor) {
    // returns true if potential_ancestor is an ancestor of this

    if (potential_ancestor == shared_from_this()) {
        return true;
    }
    if (!parent_blossom) {
        return false;
    }
    return parent_blossom->Ancestor(potential_ancestor);
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

        if (edge->slack_quadrupled < 0) {
            throw std::runtime_error("In IncreaseDualVariableQuadrupled: got a negative slack edge");
        }
    }
}

std::shared_ptr<Node> Node::SharedMaxDistinctBlossom(const std::shared_ptr<EdgeWeighted> &first_edge,
    const std::shared_ptr<EdgeWeighted> &second_edge) {
    // returns the max distinct blossom adjacent to both edges
    // if the edges don't share a max distinct blossom, returns nullptr

    auto [first_head, first_tail] = first_edge->VerticesMaxDistinctBlossoms();
    auto [second_head, second_tail] = second_edge->VerticesMaxDistinctBlossoms();
    if ((first_head == second_head) || (first_head == second_tail)) {
        return first_head;
    }
    if ((first_tail == second_head) || (first_tail == second_tail)) {
        return first_tail;
    }
    return nullptr;
}

std::shared_ptr<Node> Node::ReceptacleChild() {
    // returns the blossom child vertex that is in the outgoing matched edge

    if (parent_blossom) {
        throw std::runtime_error("In ReceptacleChild: must be called only for a top blossom");
    }
    if (!matched_edge) {
        throw std::runtime_error("In ReceptacleChild: must be called only for a matched supervertex");
    }

    auto [vertex_receptacle, vertex_out] = matched_edge->VerticesElementary();
    if (vertex_out->Ancestor(shared_from_this())) {
        std::swap(vertex_receptacle, vertex_out);
    }

    return vertex_receptacle->NextToTopBlossom();
}

std::shared_ptr<Node> Node::TreeElderChild() {
    // returns the blossom child vertex that is next to tree_parent

    if (parent_blossom) {
        throw std::runtime_error("In TreeElderChild: must be called only for a top blossom");
    }
    if (!tree_parent) {
        throw std::runtime_error("In TreeElderChild: no tree_parent");
    }

    auto [vertex_elder, vertex_out] = matched_edge->VerticesElementary();
    if (vertex_out->Ancestor(shared_from_this())) {
        std::swap(vertex_elder, vertex_out);
    }

    return vertex_elder->NextToTopBlossom();
}

void Node::MakeChildrenTopBlossoms() {
    // rewrites the parent_blossom of the blossom children

    for (const auto& child : children_blossom) {
        child->parent_blossom = nullptr;
    }
}

void Node::MakeChildrenBrotherless() {
    for (const auto& child : children_blossom) {
        child->blossom_brother_clockwise = nullptr;
        child->blossom_brother_anticlockwise = nullptr;
    }
}

void Node::ClearChildrenTreeChildren() {
    for (const auto& subvertex : children_blossom) {
        subvertex->tree_children.clear();
    }
}

void Node::UpdateChildrenExternalTreeChildren() {
    // this blossom had some tree_children, assign them to subblossoms before the Expand of this blossom
    // must be called after UpdateInternalTreeStructure()
    // also update the corresponding tree_parents

    std::unordered_set<std::shared_ptr<EdgeWeighted>> tree_children_set;
    for (const auto& child : tree_children) {
        tree_children_set.insert(child);
    }

    for (const auto& subvertex : children_blossom) {
        if (!subvertex->tree_root) { // we are in the part of the blossom that goes to waste
            // maybe we are never here, since the minus blossom only has one tree_child
            continue;
        }

        for (const auto& edge : subvertex->neighbors) {
            if (tree_children_set.contains(edge)) {
                subvertex->tree_children.push_back(edge);
                edge->OtherBlossom(subvertex)->tree_parent = edge;
            }
        }
    }
}

void Node::UpdateElderChildTreeParent() {
    // updates the tree_parent of the elder child before the Expand
    // also updates tree_children of this tree_parent

    auto elder_child = TreeElderChild();
    auto parent_edge = tree_parent;

    if (!parent_edge) {
        throw std::runtime_error("In UpdateElderChildTreeParent: no tree_parent");
    }

    elder_child->tree_parent = parent_edge;
    auto parent_vertex = parent_edge->OtherBlossom(elder_child);
    for (auto & tree_child_candidate : parent_vertex->tree_children) {
        if (tree_child_candidate->OtherBlossom(parent_vertex) == shared_from_this()) {
            tree_child_candidate = parent_edge;
            break;
        }
    }
}

void Node::UpdateInternalTreeStructure() {
    // updates tree_children, tree_parents, tree_root for the vertices inside this blossom before the Expand
    // updates the plus and minus markers

    auto elder_child = TreeElderChild();
    auto receptacle = ReceptacleChild();

    auto cur_vertex = receptacle;
    int clockwise_dist = 0;
    while (cur_vertex != elder_child) {
        cur_vertex = cur_vertex->blossom_brother_clockwise->OtherMaxDistinctBlossom(cur_vertex);
        ++clockwise_dist;
    }

    cur_vertex = receptacle->blossom_brother_clockwise->OtherMaxDistinctBlossom(receptacle);
    bool is_plus = true;
    while (cur_vertex != elder_child) {
        if (clockwise_dist % 2 == 0) {
            cur_vertex->plus = is_plus;
            cur_vertex->minus = !is_plus;
            cur_vertex->tree_root = tree_root;
            cur_vertex->tree_parent = cur_vertex->blossom_brother_clockwise;
        } else {
            cur_vertex->plus = false;
            cur_vertex->minus = false;
            cur_vertex->tree_root = nullptr;
            cur_vertex->tree_parent = nullptr;
        }

        cur_vertex = cur_vertex->blossom_brother_clockwise->OtherMaxDistinctBlossom(cur_vertex);
        is_plus = !is_plus;

        if (clockwise_dist % 2 == 0) {
            cur_vertex->tree_children.push_back(cur_vertex->blossom_brother_anticlockwise);
        } else {
            cur_vertex->tree_children.clear();
        }
    }

    cur_vertex = receptacle;
    is_plus = true;
    while (cur_vertex != elder_child) {
        if (clockwise_dist % 2 != 0) {
            cur_vertex->plus = is_plus;
            cur_vertex->minus = !is_plus;
            cur_vertex->tree_root = tree_root;
            cur_vertex->tree_parent = cur_vertex->blossom_brother_anticlockwise;
        } else {
            cur_vertex->plus = false;
            cur_vertex->minus = false;
            cur_vertex->tree_root = nullptr;
            cur_vertex->tree_parent = nullptr;
        }

        cur_vertex = cur_vertex->blossom_brother_anticlockwise->OtherMaxDistinctBlossom(cur_vertex);
        is_plus = !is_plus;

        if (clockwise_dist % 2 != 0) {
            cur_vertex->tree_children.push_back(cur_vertex->blossom_brother_clockwise);
        } else {
            cur_vertex->tree_children.clear();
        }
    }

}
