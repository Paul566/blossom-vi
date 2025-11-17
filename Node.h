#ifndef BLOSSOM_VI_NODE_H
#define BLOSSOM_VI_NODE_H

#include <vector>
#include <memory>

class EdgeWeighted;

class Node : public std::enable_shared_from_this<Node> {
    public:
        const int index; // -1 if it is a supernode (blossom)
        std::vector<std::shared_ptr<EdgeWeighted> > neighbors;
        std::shared_ptr<EdgeWeighted> matched_edge;

        // blossom related fields
        std::shared_ptr<Node> parent_blossom;
        std::vector<std::shared_ptr<Node> > children_blossom;
        std::shared_ptr<EdgeWeighted> blossom_brother_clockwise;
        std::shared_ptr<EdgeWeighted> blossom_brother_anticlockwise;

        // tree related fields
        bool plus;
        bool minus;
        std::shared_ptr<EdgeWeighted> tree_parent;
        std::vector<std::shared_ptr<EdgeWeighted> > tree_children;
        std::shared_ptr<Node> tree_root; // tree_root is always an elementary vertex (maybe inside a blossom)

        explicit Node(int index_);

        explicit Node(std::vector<std::shared_ptr<EdgeWeighted> > blossom_edges);

        std::shared_ptr<Node> TopBlossom();

        std::shared_ptr<Node> NextToTopBlossom();

        void Dissolve();

        bool Ancestor(const std::shared_ptr<Node> &potential_ancestor);

        int DualVariableQuadrupled() const;

        void IncreaseDualVariableQuadrupled(int increment);

    private:
        int dual_variable_quadrupled;

        static std::shared_ptr<Node> SharedMaxDistinctBlossom(const std::shared_ptr<EdgeWeighted> &first_edge,
                                                       const std::shared_ptr<EdgeWeighted> &second_edge) ;

        void RotateReceptacle(const std::shared_ptr<Node> &new_receptacle) const;

        std::shared_ptr<Node> ReceptacleChild();

        std::shared_ptr<Node> TreeElderChild();

        void MakeChildrenTopBlossoms();

        void MakeChildrenBrotherless();

        void ClearChildrenTreeChildren();

        void UpdateChildrenExternalTreeChildren();

        void UpdateElderChildTreeParent();

        void UpdateInternalTreeStructure();
};

#endif
