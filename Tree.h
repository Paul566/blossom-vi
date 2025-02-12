#ifndef TREE_H
#define TREE_H

#include "EdgeWeighted.h"
#include "Node.h"


class Tree {
public:
    std::shared_ptr<Node> root; // must be an elementary vertex

    explicit Tree(const std::shared_ptr<Node> &root_);

    void Grow(const std::shared_ptr<EdgeWeighted>& edge);

    void Shrink(const std::shared_ptr<EdgeWeighted>& edge_plus_plus);

    void Expand(const std::shared_ptr<Node> &supervertex);

    std::shared_ptr<EdgeWeighted> MinSlackEdgeFromPlus() const;

private:
    std::shared_ptr<Node> LCA(const std::shared_ptr<EdgeWeighted>& edge_plus_plus) const;
};

#endif //TREE_H
