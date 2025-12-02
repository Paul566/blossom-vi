#ifndef TREE_H
#define TREE_H

#include <list>
#include <unordered_map>
#include <set>


class EdgeWeighted;
class Node;

struct NodeComparator {
    bool operator()(const Node * a, const Node * b) const;
};
struct EdgeComparator {
    bool operator()(const EdgeWeighted * a, const EdgeWeighted * b) const;
};

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

        int PlusEmptySlack();
        int PlusPlusExternalSlack() const;
        int PlusPlusInternalSlack() const;
        int PlusMinusExternalSlack() const;
        int MinMinusBlossomVariable();

        void ValidatePlusEmpty() const;   // debugging purposes

    private:
        const int num_elementary_nodes; // TODO get rid of this field
        std::list<Node> *blossom_storage;
        std::unordered_map<Node *, std::list<Node>::iterator> *iter_to_self;

        std::set<Node *, NodeComparator> minus_blossoms;
        std::set<EdgeWeighted *, EdgeComparator> plus_empty_edges;
        std::set<EdgeWeighted *, EdgeComparator> plus_plus_internal_edges;
        std::set<EdgeWeighted *, EdgeComparator> plus_plus_external_edges;

        void Grow(EdgeWeighted &edge);
        void Shrink(EdgeWeighted &edge_plus_plus);
        void Expand(Node &blossom);
        Tree *Augment(EdgeWeighted &edge);

        void AugmentFromNode(Node &vertex);

        EdgeWeighted *GrowableEdge();
        EdgeWeighted *AugmentableEdge() const;
        EdgeWeighted *ShrinkableEdge() const;
        auto ExpandableBlossom() -> Node *;

        void DissolveTree();
        void UpdateQueuesAfterGrow(Node & child);
        void UpdateQueuesAfterShrink(const Node & blossom);
        void UpdateQueuesAfterExpand(const std::vector<Node *> & children);
        void UpdateQueuesAfterInit();
};

#endif //TREE_H
