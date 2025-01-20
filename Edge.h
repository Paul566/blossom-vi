#ifndef BLOSSOM_VI_EDGE_H
#define BLOSSOM_VI_EDGE_H

class Edge {
public:
    bool matched;

    Edge(const int head_, const int tail_) : head(head_), tail(tail_) {
        matched = false;
    }

    int OtherNode(const int vertex) const {
        if (vertex == head) {
            return tail;
        }
        return head;
    }

    std::pair<int, int> Vertices() const {
        return {head, tail};
    }

private:
    const int head, tail;
};

#endif
