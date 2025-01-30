#ifndef BLOSSOM_VI_NODE_H
#define BLOSSOM_VI_NODE_H

#include <vector>
#include <memory>

class EdgeWeighted;

class Node {
public:
    int dual_variable;
    std::shared_ptr<Node> parent_blossom;
    std::vector<std::shared_ptr<Node>> children_blossom;
    std::vector<std::shared_ptr<EdgeWeighted>> neighbors;

    Node();
};

#endif
