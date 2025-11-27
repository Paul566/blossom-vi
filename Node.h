#ifndef BLOSSOM_VI_NODE_H
#define BLOSSOM_VI_NODE_H

#include <list>
#include <vector>

class EdgeWeighted;
class Tree;

class Node {
    public:
        const int index; // >=n if it is a supernode (blossom)
        std::vector<EdgeWeighted *> neighbors;

        explicit Node(int index_);
        explicit Node(const std::vector<EdgeWeighted *> &blossom_edges, int index_);

        Node(const Node &other) = delete;
        Node(Node &&other) = delete;
        Node &operator=(const Node &other) = delete;
        Node &operator=(Node &&other) = delete;

        void PrintNode() const;

        int DualVariableQuadrupled() const;
        bool IsElementary() const;
        bool IsTopBlossom() const;
        Node &TopBlossom(); // TODO make const
        bool IsInSomeTree() const;
        bool IsInThisTree(const Tree & tree_) const;
        Tree * TreeOf() const;
        bool Plus() const;
        bool Minus() const;
        bool IsMatched() const;
        EdgeWeighted * TreeParentEdge() const;
        const std::vector<EdgeWeighted *> & TreeChildren() const;
        Node * BlossomParent() const;
        const std::vector<Node *> & BlossomChildren() const;

        std::vector<EdgeWeighted *> PathToRoot() const;
        static Node & LCA(Node & first, Node & second);

        static void MakeMatched(EdgeWeighted & edge);
        static void MakeUnmatched(EdgeWeighted & edge);

        void Dissolve();
        void RotateReceptacle(Node *new_receptacle) const;
        void ClearDuringTreeDissolve();

        void MakeRootOfTree(Tree & tree_);
        void MakeATreeChild(EdgeWeighted &edge_to_parent);

        void MakeSlackNonnegativeInInit();
        void InitVarGreedily();

    private:
        int dual_var_quadrupled_amortized;
        EdgeWeighted *matched_edge;

        // blossom related fields
        Node *blossom_parent;
        std::vector<Node *> blossom_children;
        EdgeWeighted *blossom_brother_clockwise;
        EdgeWeighted *blossom_brother_anticlockwise;

        // tree related fields
        bool plus;
        bool minus; // TODO get rid of this, can be inferred from plus and tree==nullptr
        EdgeWeighted *tree_parent;
        std::vector<EdgeWeighted *> tree_children;
        Tree *tree;

        static Node &SharedNode(const EdgeWeighted &first_edge, const EdgeWeighted &second_edge);

        void UpdateInternalTreeStructure();
        void ClearInternalTreeStructure() const;
};

#endif
