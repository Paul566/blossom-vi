class Edge {
public:
    bool matched;

    Edge(int head_, int tail_) : head(head_), tail(tail_) {
        matched = false;
    }

    int OtherNode(int vertex) const {
        if (vertex == head) {
            return tail;
        }
        if (vertex == tail) {
            return head;
        }
        throw std::runtime_error("in OtherNode: the vertex is not in the edge");
        return -1;
    }

    std::pair<int, int> Vertices() {
        return {head, tail};
    }

private:
    const int head, tail;
};