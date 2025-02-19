#ifndef EDGEWEIGHTED_H
#define EDGEWEIGHTED_H

#include <unordered_set>

#include "Node.h"

class EdgeWeighted {
public:
    int weight;
    int slack_quadrupled;

    EdgeWeighted(const std::shared_ptr<Node> head_, const std::shared_ptr<Node> tail_,
                 const int weight_) : head(head_), tail(tail_), weight(weight_), slack_quadrupled(4 * weight_) {
        matched = false;
    }

    bool Matched() const {
        return matched;
    }

    void MakeUnmatched() {
        matched = false;
    }

    void MakeMatched() {
        matched = true;

        // now update the matched_edge of the vertices
        const auto lca = LCABlossom();
        auto first = head;
        auto second = tail;
        while (first->parent_blossom != lca) {
            first->matched_edge = std::shared_ptr<EdgeWeighted>(this);
            first = first->parent_blossom;
        }
        while (second->parent_blossom != lca) {
            second->matched_edge = std::shared_ptr<EdgeWeighted>(this);
            second = second->parent_blossom;
        }
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

        if (head_top == tail_top) {
            throw std::runtime_error("In OtherBlossom: top blossoms are the same, you must call OtherMaxDistinctBlossom instead");
        }

        if (vertex_top == head_top) {
            return tail_top;
        }
        if (vertex_top == tail_top) {
            return head_top;
        }

        throw std::runtime_error("In EdgeWeighted::OtherBlossom: vertex.TopBlossom is not in the edge");
    }

    std::shared_ptr<Node> OtherMaxDistinctBlossom(const std::shared_ptr<Node> &vertex) const {
        auto [head_max, tail_max] = VerticesMaxDistinctBlossoms();

        if (vertex->Ancestor(head_max)) {
            return tail_max;
        }
        if (vertex->Ancestor(tail_max)) {
            return head_max;
        }

        throw std::runtime_error("In EdgeWeighted::OtherBlossom: vertex is not an ancestor of either of VerticesMaxDistinctBlossoms");
    }

    std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>> VerticesElementary() const {
        return {head, tail};
    }

    std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>> VerticesTopBlossoms() const {
        return {head->TopBlossom(), tail->TopBlossom()};
    }

    std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>> VerticesMaxDistinctBlossoms() const {
        // returns blossoms u, v, containing the edge such that
        // either both u, v are top blossoms, or
        // u and v share a parent blossom

        auto lca = LCABlossom();

        if (lca == nullptr) { // top blossoms of head and tail are distinct
            return VerticesTopBlossoms();
        }

        auto first = head;
        auto second = tail;
        while (first->parent_blossom != lca) {
            first = first->parent_blossom;
        }
        while (second->parent_blossom != lca) {
            second = second->parent_blossom;
        }

        if (first == second) {
            throw std::runtime_error("In VerticesMaxDistinctBlossoms: first == second");
        }
        if (first->parent_blossom != second->parent_blossom) {
            throw std::runtime_error("In VerticesMaxDistinctBlossoms: the parent is not common");
        }

        return {first, second};
    }

private:
    const std::shared_ptr<Node> head, tail;

    bool matched;

    std::shared_ptr<Node> LCABlossom() const {
        // returns nullptr if head and tail are in distinct top blossoms

        auto first = head;
        auto second = tail;

        std::unordered_set<std::shared_ptr<Node>> visited;
        visited.insert(first);
        while (first->parent_blossom) {
            first = first->parent_blossom;
            visited.insert(first);
        }

        if (visited.contains(second)) {
            throw std::runtime_error("In LCABlossom: second is an ancestor of first");
        }

        while (second->parent_blossom) {
            second = second->parent_blossom;
            if (visited.contains(second)) {
                return second;
            }
        }
        return nullptr;
    }
};

#endif //EDGEWEIGHTED_H
