#ifndef BLOSSOM_VI_NODE_H
#define BLOSSOM_VI_NODE_H

#include <vector>
#include <memory>

class EdgeWeighted;

class Node {
public:
    const int index;  // -1 if it is a supernode (blossom)
    bool plus;
    bool minus;
    std::shared_ptr<EdgeWeighted> matched_edge;
    std::shared_ptr<Node> parent_blossom;
    const std::vector<std::shared_ptr<Node>> children_blossom;
    std::vector<std::shared_ptr<EdgeWeighted>> neighbors;
    std::shared_ptr<EdgeWeighted> tree_parent;
    std::vector<std::shared_ptr<EdgeWeighted>> tree_children;
    std::shared_ptr<Node> tree_root;    // tree_root is always an elementary vertex (maybe inside a blossom)

    explicit Node(int index_);

    explicit Node(std::vector<std::shared_ptr<Node>> blossom_vertices);

    std::shared_ptr<Node> TopBlossom();

    int DualVariableQuadrupled() const;

    void IncreaseDualVariableQuadrupled(int increment);

private:
    int dual_variable_quadrupled;

};

#endif
