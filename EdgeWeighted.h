#ifndef EDGEWEIGHTED_H
#define EDGEWEIGHTED_H

#include <iostream>
#include <stdexcept>
#include <vector>

#include "Node.h"

class EdgeWeighted {
    public:
        const int weight;
        bool matched;

        EdgeWeighted(Node &head_, Node &tail_, int weight_);

        EdgeWeighted(const EdgeWeighted &other) = delete;
        EdgeWeighted(EdgeWeighted &&other) = delete;
        EdgeWeighted &operator=(const EdgeWeighted &other) = delete;
        EdgeWeighted &operator=(EdgeWeighted &&other) = delete;

        int SlackQuadrupled() const;

        std::pair<Node &, Node &> Endpoints() const;

        Node &OtherEnd(const Node &vertex) const;
        Node &DeeperNode(const Node &vertex);
        bool IsInsideBlossom() const;

        void UpdateAfterShrink(const Node &vertex);
        void UpdateAfterDissolve(const Node &vertex);

    private:
        int slack_quadrupled_amortized;
        Node * head;
        Node * tail;
        // head, tail: top blossoms
};

#endif //EDGEWEIGHTED_H
