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

    std::shared_ptr<Node> OtherNode(const std::shared_ptr<Node> vertex) const {
        if (vertex == head) {
            return tail;
        }
        if (vertex == tail) {
            return head;
        }
        throw std::runtime_error("In OtherNode weighted: vertex is not in the edge");
    }

    std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>> Vertices() const {
        return {head, tail};
    }

private:
    std::shared_ptr<Node> head, tail;
};

#endif //EDGEWEIGHTED_H
