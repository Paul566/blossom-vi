#ifndef TREE_H
#define TREE_H

#include <list>
#include <unordered_map>

// #include "EdgeWeighted.h"
// #include "Node.h"
class EdgeWeighted;
class Node;

class Tree {
    public:
        Node *root; // must be an elementary vertex
        // TODO maybe the top blossom of the root makes more sense here
        int dual_var_quadrupled;

        Tree(Node *root_,
             std::list<Node> *blossom_storage_,
             std::unordered_map<Node *, std::list<Node>::iterator> *iter_to_self_,
             int num_elementary_nodes_);

        Tree(const Tree &other) = delete;
        Tree(Tree &&other) = delete;
        Tree &operator=(const Tree &other) = delete;
        Tree &operator=(Tree &&other) = delete;

        void PrintTree();

        Tree *MakePrimalUpdates();

        int PlusEmptySlack() const;

        int PlusPlusExternalSlack() const;

        int PlusPlusInternalSlack() const;

        int PlusMinusExternalSlack() const;

        int MinMinusBlossomVariable() const;

    private:
        const int num_elementary_nodes; // TODO get rid of this field
        std::list<Node> *blossom_storage;
        std::unordered_map<Node *, std::list<Node>::iterator> *iter_to_self;

        void Grow(EdgeWeighted &edge);
        void Shrink(EdgeWeighted &edge_plus_plus) const;    // TODO maybe remove/tie parallel edges after shrinking
        void Expand(Node &blossom) const;
        Tree *Augment(EdgeWeighted &edge);

        void AugmentFromNode(Node &vertex);

        EdgeWeighted *GrowableEdge() const;
        EdgeWeighted *AugmentableEdge() const;
        EdgeWeighted *ShrinkableEdge() const;
        Node *ExpandableBlossom();

        Node &LCA(const EdgeWeighted &edge_plus_plus) const;

        void DissolveTree();

        static std::vector<EdgeWeighted *> PathToRoot(const Node &vertex);
};

#endif //TREE_H
