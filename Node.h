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
        EdgeWeighted *matched_edge;

        // blossom related fields
        Node *blossom_parent;
        std::vector<Node *> blossom_children;
        EdgeWeighted *blossom_brother_clockwise;
        EdgeWeighted *blossom_brother_anticlockwise;

        // tree related fields
        bool plus;
        bool minus;
        EdgeWeighted *tree_parent;
        std::vector<EdgeWeighted *> tree_children;
        Tree *tree;

        explicit Node(int index_);

        explicit Node(const std::vector<EdgeWeighted *> &blossom_edges, int index_);

        Node(const Node &other) = delete;
        Node(Node &&other) = delete;
        Node &operator=(const Node &other) = delete;
        Node &operator=(Node &&other) = delete;

        Node &TopBlossom();

        void Dissolve();

        void RotateReceptacle(Node *new_receptacle) const;

        int DualVariableQuadrupled() const;

        void IncreaseDualVariableQuadrupled(int increment);

        static void MakeMatched(EdgeWeighted & edge);

        static void MakeUnmatched(EdgeWeighted & edge);

    private:
        int dual_variable_quadrupled;

        Node &SharedNode(const EdgeWeighted &first_edge, const EdgeWeighted &second_edge);

        void UpdateInternalTreeStructure();

        void ClearInternalTreeStructure() const;
};

#endif
