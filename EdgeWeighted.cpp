#include "EdgeWeighted.h"

EdgeWeighted::EdgeWeighted(Node &head_, Node &tail_, const int weight_) : weight(weight_), matched(false),
                                                                          slack_quadrupled_amortized(4 * weight_) {
    head_stack.push_back(&head_);
    tail_stack.push_back(&tail_);
}

int EdgeWeighted::SlackQuadrupled() const {
    return slack_quadrupled_amortized - head_stack.back()->DualVariableQuadrupled() -
        tail_stack.back()->DualVariableQuadrupled();
}

std::pair<Node &, Node &> EdgeWeighted::Endpoints() const {
    return {*head_stack.back(), *tail_stack.back()};
}

std::pair<Node &, Node &> EdgeWeighted::ElementaryEndpoints() const {
    return {*head_stack.front(), *tail_stack.front()};
}

Node & EdgeWeighted::OtherEnd(const Node &vertex) const {
    if (head_stack.empty() || tail_stack.empty()) {
        throw std::runtime_error("In OtherEnd: invalid edge");
    }
    if (&vertex == tail_stack.back()) {
        return *head_stack.back();
    }
    if (&vertex == head_stack.back()) {
        return *tail_stack.back();
    }
    throw std::runtime_error("EdgeWeighted::OtherEnd: vertex is not on top of either stack");
}

Node & EdgeWeighted::OtherElementaryEnd(const Node &vertex) const {
    if (vertex.index == -1) {
        throw std::runtime_error("In OtherElementaryEnd: vertex must be elementary");
    }
    if (&vertex == tail_stack.front()) {
        return *head_stack.front();
    }
    if (&vertex == head_stack.front()) {
        return *tail_stack.front();
    }
    throw std::runtime_error("EdgeWeighted::OtherEnd: vertex is not in the front of either stack");
}

Node & EdgeWeighted::DeeperNode(const Node &vertex) const {
    // returns a node that is a blossom child of vertex and is adjacent to this edge
    if (&vertex == tail_stack.back()) {
        if (tail_stack.size() < 2) {
            throw std::runtime_error("EdgeWeighted::DeeperNode: vertex is at the bottom of the stack");
        }
        return *tail_stack[tail_stack.size() - 2];
    }
    if (&vertex == head_stack.back()) {
        if (head_stack.size() < 2) {
            throw std::runtime_error("EdgeWeighted::DeeperNode: vertex is at the bottom of the stack");
        }
        return *head_stack[head_stack.size() - 2];
    }
    throw std::runtime_error("EdgeWeighted::DeeperNode: vertex is not on top of either stack");
}

bool EdgeWeighted::IsInsideBlossom() const {
    return !head_stack.back()->IsTopBlossom();
}

void EdgeWeighted::UpdateAfterShrink(const Node &vertex) {
    if (vertex.IsTopBlossom()) {
        throw std::runtime_error("EdgeWeighted::UpdateAfterShrink: vertex has no blossom parent");
    }
    if (!vertex.IsInSomeTree()) {
        throw std::runtime_error("EdgeWeighted::UpdateAfterShrink: vertex is not in a tree");
    }

    slack_quadrupled_amortized -= vertex.DualVariableQuadrupled();
    if (&vertex == tail_stack.back()) {
        tail_stack.push_back(vertex.BlossomParent());
        return;
    }
    if (&vertex == head_stack.back()) {
        head_stack.push_back(vertex.BlossomParent());
        return;
    }

    const auto [head, tail] = ElementaryEndpoints();
    std::cout << "EdgeWeighted::UpdateAfterShrink: edge: " << head.index << " " << tail.index << std::endl;
    std::cout << "vertex: " << vertex.index << std::endl;
    throw std::runtime_error("EdgeWeighted::UpdateAfterShrink: vertex is not on top of either stack");
}

void EdgeWeighted::UpdateAfterDissolve(const Node &vertex) {
    if (&vertex == tail_stack.back()) {
        tail_stack.pop_back();
        slack_quadrupled_amortized += tail_stack.back()->DualVariableQuadrupled();
        return;
    }
    if (&vertex == head_stack.back()) {
        head_stack.pop_back();
        slack_quadrupled_amortized += head_stack.back()->DualVariableQuadrupled();
        return;
    }
    throw std::runtime_error("EdgeWeighted::UpdateAfterDissolve: vertex is not on top of either stack");
}
