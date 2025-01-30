#ifndef BLOSSOM_VI_NODE_H
#define BLOSSOM_VI_NODE_H

#include <vector>
#include <memory>

class EdgeWeighted;

class Node {
public:
    int dual_variable_quadrupled;
    const int index;  // -1 if it is a supernode (blossom)
    std::shared_ptr<EdgeWeighted> matched_edge;
    std::shared_ptr<Node> parent_blossom;
    std::vector<std::shared_ptr<Node>> children_blossom;
    std::vector<std::shared_ptr<EdgeWeighted>> neighbors;

    Node(int index_);
};

#endif
