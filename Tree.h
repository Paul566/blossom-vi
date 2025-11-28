#ifndef TREE_H
#define TREE_H

#include <list>
#include <unordered_map>
#include <set>

#include "Node.h"

class EdgeWeighted;
class Node;

struct NodeComparator {
    bool operator()(const Node * a, const Node * b) const {
        if (a->DualVariableQuadrupled() != b->DualVariableQuadrupled()) {
            return a->DualVariableQuadrupled() < b->DualVariableQuadrupled();
        }
        return a < b;
    }
};

class Tree {
    public:
        Node *root; // must be an elementary vertex
        // TODO maybe the top blossom of the root makes more sense here
        int dual_var_quadrupled;
        std::set<Node *, NodeComparator> minus_blossoms;
        // TODO update wrt minus_blossoms

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
        void Shrink(EdgeWeighted &edge_plus_plus);
        void Expand(Node &blossom);
        Tree *Augment(EdgeWeighted &edge);

        void AugmentFromNode(Node &vertex);

        EdgeWeighted *GrowableEdge() const;
        EdgeWeighted *AugmentableEdge() const;
        EdgeWeighted *ShrinkableEdge() const;
        Node *ExpandableBlossom() const;

        void DissolveTree();
};

#endif //TREE_H
