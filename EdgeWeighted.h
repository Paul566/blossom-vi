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

        EdgeWeighted(Node &head_, Node &tail_, const int weight_);

        EdgeWeighted(const EdgeWeighted &other) = delete;
        EdgeWeighted(EdgeWeighted &&other) = delete;
        EdgeWeighted &operator=(const EdgeWeighted &other) = delete;
        EdgeWeighted &operator=(EdgeWeighted &&other) = delete;

        int SlackQuadrupled() const;

        std::pair<Node &, Node &> Endpoints() const;

        std::pair<Node &, Node &> ElementaryEndpoints() const;

        Node &OtherEnd(const Node &vertex) const;

        Node &OtherElementaryEnd(const Node &vertex) const;

        Node &DeeperNode(const Node &vertex) const;

        bool IsInsideBlossom() const;

        void UpdateAfterShrink(const Node &vertex);

        void UpdateAfterDissolve(const Node &vertex);

    private:
        int slack_quadrupled_amortized;
        std::vector<Node *> head_stack;
        std::vector<Node *> tail_stack;
        // stack tops: current top blossoms or children of the blossom that contains this edge
        // in other words, maximal blossoms that contain the elementary endpoints but don't contain the entire edge
        // stack bottoms: elementary nodes
};

#endif //EDGEWEIGHTED_H
