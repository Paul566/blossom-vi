#include "EdgeWeighted.h"

#include "Tree.h"

EdgeWeighted::EdgeWeighted(Node &head_, Node &tail_, const int weight_) : weight(weight_), matched(false),
                                                                          slack_quadrupled_amortized(4 * weight_) {
    head = &head_;
    tail = &tail_;
}

int EdgeWeighted::SlackQuadrupled() const {
    return slack_quadrupled_amortized - head->DualVariableQuadrupled() -
        tail->DualVariableQuadrupled();
}

std::pair<Node &, Node &> EdgeWeighted::Endpoints() const {
    return {*head, *tail};
}

Node & EdgeWeighted::OtherEnd(const Node &vertex) const {
    if (&vertex == tail) {
        return *head;
    }
    if (&vertex == head) {
        return *tail;
    }
    throw std::runtime_error("EdgeWeighted::OtherEnd: vertex is not adjacent to this edge");
}

Node & EdgeWeighted::DeeperNode(const Node &vertex) {
    // returns a node that is a blossom child of vertex and is adjacent to this edge
    if (&vertex == tail) {
        if (tail->IsElementary()) {
            throw std::runtime_error("EdgeWeighted::DeeperNode: vertex is elementary");
        }
        return *tail->edge_to_deeper_node[this];
    }
    if (&vertex == head) {
        if (head->IsElementary()) {
            throw std::runtime_error("EdgeWeighted::DeeperNode: vertex is elementary");
        }
        return *head->edge_to_deeper_node[this];
    }
    throw std::runtime_error("EdgeWeighted::DeeperNode: vertex is not on top of either stack");
}

bool EdgeWeighted::IsInsideBlossom() const {
    return !head->IsTopBlossom();
}

void EdgeWeighted::UpdateAfterShrink(const Node &vertex) {
    if (vertex.IsTopBlossom()) {
        throw std::runtime_error("EdgeWeighted::UpdateAfterShrink: vertex has no blossom parent");
    }
    if (!vertex.IsInSomeTree()) {
        throw std::runtime_error("EdgeWeighted::UpdateAfterShrink: vertex is not in a tree");
    }

    slack_quadrupled_amortized -= vertex.DualVariableQuadrupled();
    if (&vertex == tail) {
        tail = vertex.BlossomParent();
        return;
    }
    if (&vertex == head) {
        head = vertex.BlossomParent();
        return;
    }

    std::cout << "EdgeWeighted::UpdateAfterShrink: edge: " << head->index << " " << tail->index << std::endl;
    std::cout << "vertex: " << vertex.index << std::endl;
    throw std::runtime_error("EdgeWeighted::UpdateAfterShrink: vertex is not on top of either stack");
}

void EdgeWeighted::UpdateAfterDissolve(const Node &vertex) {
    if (&vertex == tail) {
        tail = tail->edge_to_deeper_node[this];
        slack_quadrupled_amortized += tail->DualVariableQuadrupled();
        return;
    }
    if (&vertex == head) {
        head = head->edge_to_deeper_node[this];
        slack_quadrupled_amortized += head->DualVariableQuadrupled();
        return;
    }
    throw std::runtime_error("EdgeWeighted::UpdateAfterDissolve: vertex is not on top of either stack");
}
