#ifndef EDGEWEIGHTED_H
#define EDGEWEIGHTED_H

#include "Node.h"

class EdgeWeighted {
public:
    bool matched;
    int weight;
    int slack_quadrupled;

    EdgeWeighted(const std::shared_ptr<Node> head_, const std::shared_ptr<Node> tail_,
                 const int weight_) : head(head_), tail(tail_), weight(weight_), slack_quadrupled(4 * weight_) {
        matched = false;
    }

    std::shared_ptr<Node> OtherElementary(const std::shared_ptr<Node> &vertex) const {
        if (vertex->index == -1) {
            throw std::runtime_error("In EdgeWeighted::OtherElementary: vertex must be elementary");
        }

        if (vertex == head) {
            return tail;
        }
        if (vertex == tail) {
            return head;
        }

        throw std::runtime_error("In EdgeWeighted::OtherElementary: vertex is neither head nor tail");
    }

    std::shared_ptr<Node> OtherBlossom(const std::shared_ptr<Node> &vertex) const {
        auto vertex_top = vertex->TopBlossom();
        auto head_top = head->TopBlossom();
        auto tail_top = tail->TopBlossom();

        if (vertex_top == head_top) {
            return tail_top;
        }
        if (vertex_top == tail_top) {
            return head_top;
        }

        throw std::runtime_error("In EdgeWeighted::OtherBlossom: vertex.TopBlossom is not in the edge");
    }

    std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>> VerticesElementary() const {
        return {head, tail};
    }

    std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>> VerticesTopBlossoms() const {
        return {head->TopBlossom(), tail->TopBlossom()};
    }

private:
    const std::shared_ptr<Node> head, tail;
};

#endif //EDGEWEIGHTED_H
