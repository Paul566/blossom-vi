#ifndef TREE_H
#define TREE_H

#include <list>
#include <unordered_map>

#include "EdgeWeighted.h"
#include "Node.h"

class Tree {
    public:
        Node *root; // must be an elementary vertex

        Tree(Node *root_, std::list<Node> *blossom_storage_, std::unordered_map<Node *, std::list<Node>::iterator> * iter_to_self_);

        Tree(const Tree &other) = delete;
        Tree(Tree &&other) = delete;
        Tree &operator=(const Tree &other) = delete;
        Tree &operator=(Tree &&other) = delete;

        void Grow(EdgeWeighted &edge);

        void Shrink(EdgeWeighted &edge_plus_plus) const;

        void Expand(Node &blossom) const;

        void Augment(EdgeWeighted &edge);

        bool MakePrimalUpdate(bool * is_augmented);

        void MakeDualUpdate();

    private:
        std::list<Node> *blossom_storage;
        std::unordered_map<Node *, std::list<Node>::iterator> * iter_to_self;

        // void AugmentFromNode(Node &vertex);

        EdgeWeighted &MinSlackEdgeFromPlus() const;

        Node *ExpandableBlossom();

        Node *MinYMinusBlossom() const;

        void ChangeDualVariables(int increment) const;

        Node &LCA(const EdgeWeighted &edge_plus_plus) const;

        void DissolveTree();

        static std::vector<EdgeWeighted *> PathToRoot(const Node &vertex);
};

#endif //TREE_H
