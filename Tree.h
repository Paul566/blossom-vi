#ifndef TREE_H
#define TREE_H

#include "EdgeWeighted.h"
#include "Node.h"

class Tree {
    public:
        std::shared_ptr<Node> root; // must be an elementary vertex

        explicit Tree(const std::shared_ptr<Node> &root_);

        void Grow(const std::shared_ptr<EdgeWeighted> &edge) const;

        void Shrink(const std::shared_ptr<EdgeWeighted> &edge_plus_plus) const;

        void Expand(const std::shared_ptr<Node> &supervertex) const;

        void Augment(const std::shared_ptr<EdgeWeighted> &edge);

        bool MakePrimalUpdate();

        void MakeDualUpdate();

    private:
        std::shared_ptr<EdgeWeighted> MinSlackEdgeFromPlus() const;

        std::shared_ptr<Node> ExpandableBlossom();

        std::shared_ptr<Node> MinYMinusBlossom();

        void ChangeDualVariables(int increment) const;

        std::shared_ptr<Node> LCA(const std::shared_ptr<EdgeWeighted> &edge_plus_plus) const;

        void DissolveTree() const;

        static std::vector<std::shared_ptr<EdgeWeighted>> PathToRoot(std::shared_ptr<Node> vertex) ;
};

#endif //TREE_H
