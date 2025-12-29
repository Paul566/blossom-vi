#ifndef BLOSSOM_VI_SOLVER_H
#define BLOSSOM_VI_SOLVER_H

#include <list>
#include <memory>
#include <vector>
#include "Heap.h"

struct SolverParameters {
    bool compute_dual_certificate = false;
    bool verbose = false;
    bool print_statistics = true;
    int heap_arity = 2;
};

class Solver {
    public:
        int64_t primal_objective;
        int64_t dual_objective;
        const SolverParameters params;

        explicit Solver(const std::vector<std::tuple<int, int, int> > &edge_list_,
                        const SolverParameters &params_ = {});

        void FindMinPerfectMatching();

        const std::vector<std::pair<int, int> > &Matching() const;
        const std::vector<std::tuple<int, int, int> > &DualCertificate() const;
        // a vector of (quadrupled dual variable, index of the blossom parent or -1)

    private:
        // to not shoot myself in the foot:
        struct NodeIndex {
            int index;
            explicit NodeIndex(int index_) : index(index_) {
            }
            explicit operator bool() const {
                return index >= 0;
            }
            bool operator==(const NodeIndex &other) const {
                return index == other.index;
            }
        };
        struct EdgeIndex {
            int index;
            explicit EdgeIndex(int index_) : index(index_) {
            }
            explicit operator bool() const {
                return index >= 0;
            }
            bool operator==(const EdgeIndex &other) const {
                return index == other.index;
            }
        };
        struct TreeIndex {
            int index;
            explicit TreeIndex(int index_) : index(index_) {
            }
            explicit operator bool() const {
                return index >= 0;
            }
            bool operator==(const TreeIndex &other) const {
                return index == other.index;
            }
        };

        struct NodeComparator {
            const Solver *solver;
            bool operator()(const NodeIndex &a, const NodeIndex &b) const;
        };
        struct EdgeComparator {
            const Solver *solver;
            bool operator()(const EdgeIndex &a, const EdgeIndex &b) const;
        };
        using EdgeHeap = Heap<
            EdgeIndex,
            EdgeComparator
        >;
        using NodeHeap = Heap<
            NodeIndex,
            NodeComparator
        >;
        struct Queues {
            // TODO maybe optimize the heaps a bit
            std::vector<std::unique_ptr<NodeHeap> > node_heaps;
            std::vector<std::unique_ptr<EdgeHeap> > edge_heaps;
        };

        struct Edges {
            std::vector<int> weight;
            std::vector<bool> matched;
            std::vector<int> slack_quadrupled_amortized;
            std::vector<NodeIndex> head;
            std::vector<NodeIndex> tail;
            std::vector<int> queue_index;
            std::vector<EdgeHeap::Handle *> handle;
        };

        struct Nodes {
            // TODO consider not initializing fields irrelevant for elementary vertices for elementary vertices

            std::vector<bool> is_alive;
            std::vector<std::list<EdgeIndex> > neighbors;
            // concatenated neighbors of blossom_children, might contain loops
            // TODO try caching the neighbors into std::vector if we already went through them once
            std::vector<std::vector<std::list<EdgeIndex>::iterator> > children_neighbors_boundaries;
            // a vector of breakpoints to assign neighbors to children in Dissolve
            std::vector<int> dual_var_quadrupled_amortized;
            std::vector<EdgeIndex> matched_edge;

            // blossom related fields
            std::vector<NodeIndex> blossom_parent; // Node(-1) if no blossom_parent
            std::vector<std::vector<NodeIndex> > blossom_children;
            std::vector<EdgeIndex> blossom_edge_clockwise;
            std::vector<EdgeIndex> blossom_edge_anticlockwise;
            std::vector<NodeIndex> blossom_sibling_clockwise;
            std::vector<NodeIndex> blossom_sibling_anticlockwise;

            // tree related fields
            std::vector<bool> plus;
            std::vector<EdgeIndex> tree_parent;
            std::vector<std::vector<EdgeIndex> > tree_children;
            std::vector<TreeIndex> tree;

            // queue related fields
            std::vector<int> queue_index;
            std::vector<NodeHeap::Handle *> handle;
        };

        struct Trees {
            std::vector<bool> is_alive;
            std::vector<NodeIndex> root; // top blossom
            std::vector<int> dual_var_quadrupled;
            std::vector<TreeIndex> next_alive_tree;
            std::vector<int> alive_index;
            // TODO use std::vector of size num_trees_alive that holds the indices of alive trees instead of next_alive_tree

            // indices of the heaps
            std::vector<int> minus_blossoms;
            std::vector<int> plus_empty_edges;
            std::vector<int> plus_plus_internal_edges;
            std::vector<std::vector<std::pair<TreeIndex, int> > > pq_plus_plus;
            std::vector<std::vector<std::pair<TreeIndex, int> > > pq_plus_minus;
            std::vector<std::vector<std::pair<TreeIndex, int> > > pq_minus_plus;
            // TODO delete the queues that lead to dead trees
        };

        struct DualConstraints {
            std::vector<int> upper_bound;
            // upper bound on delta_T quadrupled
            std::vector<std::vector<std::pair<int, int> > > plus_plus_constraints;
            // (other_tree_index_alive, slack_plus_plus quadrupled)
            std::vector<std::vector<std::pair<int, int> > > plus_minus_constraints;
            // (other_tree_index_alive, slack_plus_minus quadrupled)
        };

        Nodes nodes;
        Edges edges;
        Trees trees;
        Queues queues;

        const int num_vertices_elementary; // the original number of vertices (i.e. not counting blossoms)
        int num_trees_alive;
        TreeIndex first_alive_tree;

        std::vector<std::pair<int, int> > matching;
        std::vector<std::tuple<int, int, int> > dual_certificate;
        // a vector of (quadrupled dual variable, index of the blossom parent or -1)
        // dual_certificate is empty unless params.compute_dual_certificate is true

        void PrintGraph() const;
        void PrintNode(NodeIndex node) const;
        void PrintTrees() const;
        void PrintStatistics();

        void GreedyInit();
        void InitializeTrees();
        void InitializeQueues();

        void ComputeMatching();
        void ComputeDualCertificate();
        void DestroyBlossoms();
        int64_t DualObjectiveQuadrupled() const;
        int64_t PrimalObjective() const;

        void MakeDualUpdates();
        std::vector<int> VariableDeltas();
        std::vector<std::vector<int> > ConnectedComponentsTreeTree(const DualConstraints &dual_constraints) const;
        DualConstraints GetDualConstraints();

        void MakePrimalUpdates();
        TreeIndex MakePrimalUpdate(TreeIndex tree, bool *success);
        void Grow(TreeIndex tree, EdgeIndex edge);
        void Shrink(TreeIndex tree, EdgeIndex edge_plus_plus);
        void MakeBlossom(std::vector<EdgeIndex> blossom_edges, NodeIndex lca);
        void Expand(TreeIndex tree, NodeIndex blossom);
        TreeIndex Augment(TreeIndex tree, EdgeIndex edge);
        void AugmentFromNode(TreeIndex tree, NodeIndex node);
        void DissolveTree(TreeIndex tree);
        void ClearNodeDuringTreeDissolve(TreeIndex tree, NodeIndex node);
        void Dissolve(NodeIndex blossom);

        void UpdateEdgeAfterShrink(EdgeIndex edge);

        void UpdateNodeInternalTreeStructure(NodeIndex node, NodeIndex receptacle, NodeIndex elder_child);
        void ClearNodeInternalTreeStructure(NodeIndex node);
        void RotateReceptacle(NodeIndex blossom, NodeIndex child);
        NodeIndex DeeperNode(NodeIndex blossom, EdgeIndex edge) const;

        void UpdateQueuesAfterGrow(NodeIndex child, NodeIndex grandchild);
        void UpdateQueuesBeforeShrink(const std::vector<NodeIndex> &minus_children);
        void RemoveLoopsFromQueues(NodeIndex blossom); // TODO get rid of this
        void UpdateQueuesAfterExpand(TreeIndex tree, NodeIndex blossom, const std::vector<NodeIndex> &children);
        int TreeTreeQueueIndex(TreeIndex other_tree, const std::vector<std::pair<TreeIndex, int> > &tree_neighbors);
        void AddPQPlusPlus(TreeIndex first, TreeIndex second, EdgeIndex edge);
        void AddPQPlusMinus(TreeIndex tree_plus, TreeIndex tree_minus, EdgeIndex edge);
        void ValidatePositiveSlacks(); // debugging purposes
        void ValidatePositiveVars(); // debugging purposes
        void ValidateQueues();  // debugging purposes


        EdgeIndex MinPlusEmptyEdge(int queue_index);
        EdgeIndex MinPlusPlusInternalEdge(int queue_index);
        EdgeIndex MinPlusPlusExternalEdge(int queue_index);
        EdgeIndex MinPlusMinusExternalEdge(int queue_index);
        NodeIndex MinMinusBlossom(int queue_index);
        EdgeIndex GrowableEdge(TreeIndex tree);
        EdgeIndex AugmentableEdge(TreeIndex tree);
        EdgeIndex ShrinkableEdge(TreeIndex tree);
        NodeIndex ExpandableBlossom(TreeIndex tree);
        int PlusEmptySlack(TreeIndex tree);
        int PlusPlusInternalSlack(TreeIndex tree);
        int MinMinusBlossomVariable(TreeIndex tree);
        std::vector<std::pair<TreeIndex, int> > PlusPlusExternalSlacks(TreeIndex tree);
        std::vector<std::pair<TreeIndex, int> > PlusMinusExternalSlacks(TreeIndex tree);

        bool IsElementary(NodeIndex node) const;
        bool IsTopBlossom(NodeIndex node) const;
        NodeIndex TopBlossom(NodeIndex node) const;
        int NodeDepth(NodeIndex node) const;
        int DualVariableQuadrupled(NodeIndex node) const;

        NodeIndex LCA(NodeIndex first, NodeIndex second) const;
        std::vector<EdgeIndex> PathToRoot(TreeIndex tree, NodeIndex node);

        int SlackQuadrupled(EdgeIndex edge) const;
        NodeIndex OtherEnd(EdgeIndex edge, NodeIndex node) const;
        NodeIndex SharedNode(EdgeIndex first_edge, EdgeIndex second_edge);

        void MakeEdgeMatched(EdgeIndex edge);
        void MakeEdgeUnmatched(EdgeIndex edge);
        void AddEdgeToQueue(EdgeIndex edge, int queue_index);
        void RemoveEdgeFromQueue(EdgeIndex edge);
        void AddNodeToQueue(NodeIndex node, int queue_index);
        void RemoveNodeFromQueue(NodeIndex node);

        TreeIndex NextAliveTree(TreeIndex current);
        TreeIndex FirstAliveTree();
        void UpdateAliveIndices();

        static int InitNumVertices(const std::vector<std::tuple<int, int, int> > &edge_list_);
};

#endif //BLOSSOM_VI_SOLVER_H
