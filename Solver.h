#ifndef BLOSSOM_VI_SOLVER_H
#define BLOSSOM_VI_SOLVER_H

#include <list>
#include <unordered_map>
#include <vector>
#include <boost/heap/d_ary_heap.hpp>

struct SolverParameters {
    bool compute_dual_certificate = false;
    bool verbose = false;
    bool print_statistics = true;
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
        struct Node {
            int index;
            explicit Node(int index_) : index(index_) {
            }
            explicit operator bool() const {
                return index >= 0;
            }
            bool operator==(const Node &other) const {
                return index == other.index;
            }
        };
        struct Edge {
            int index;
            explicit Edge(int index_) : index(index_) {
            }
            explicit operator bool() const {
                return index >= 0;
            }
            bool operator==(const Edge &other) const {
                return index == other.index;
            }
        };
        struct Tree {
            int index;
            explicit Tree(int index_) : index(index_) {
            }
            explicit operator bool() const {
                return index >= 0;
            }
            bool operator==(const Tree &other) const {
                return index == other.index;
            }
        };

        struct NodeComparator {
            const Solver *solver;
            bool operator()(const Node &a, const Node &b) const;
        };
        struct EdgeComparator {
            const Solver *solver;
            bool operator()(const Edge &a, const Edge &b) const;
        };
        using EdgeHeap = boost::heap::d_ary_heap<
            Edge,
            boost::heap::arity<2>,
            boost::heap::compare<EdgeComparator>,
            boost::heap::mutable_<true>
        >;
        using NodeHeap = boost::heap::d_ary_heap<
            Node,
            boost::heap::arity<2>,
            boost::heap::compare<NodeComparator>,
            boost::heap::mutable_<true>
        >;
        struct Queues {
            std::vector<std::unique_ptr<NodeHeap>> node_heaps;
            std::vector<std::unique_ptr<EdgeHeap>> edge_heaps;
        };

        struct Edges {
            std::vector<int> weight;
            std::vector<bool> matched;
            std::vector<int> slack_quadrupled_amortized;
            std::vector<Node> head;
            std::vector<Node> tail;
            std::vector<int> queue_index;
            std::vector<EdgeHeap::handle_type> handle;
        };

        struct Nodes {
            // TODO consider not initializing fields irrelevant for elementary vertices for elementary vertices

            std::vector<bool> is_alive;
            std::vector<std::list<Edge> > neighbors; // concatenated neighbors of blossom_children, might contain loops
            // TODO try caching the neighbors into std::vector if we already went through them once
            std::vector<std::vector<std::list<Edge>::iterator> > children_neighbors_boundaries;
            // a vector of breakpoints to assign neighbors to children in Dissolve
            std::vector<int> dual_var_quadrupled_amortized;
            std::vector<Edge> matched_edge;

            // blossom related fields
            std::vector<Node> blossom_parent; // Node(-1) if no blossom_parent
            std::vector<std::vector<Node> > blossom_children;
            std::vector<Edge> blossom_brother_clockwise;
            std::vector<Edge> blossom_brother_anticlockwise;

            // tree related fields
            std::vector<bool> plus;
            std::vector<Edge> tree_parent;
            std::vector<std::vector<Edge> > tree_children;
            std::vector<Tree> tree;

            // queue related fields
            std::vector<int> queue_index;
            std::vector<NodeHeap::handle_type> handle;
        };

        struct Trees {
            std::vector<bool> is_alive;
            std::vector<Node> root; // top blossom
            std::vector<int> dual_var_quadrupled;
            std::vector<Tree> next_alive_tree;
            std::vector<int> alive_index;
            // TODO use std::vector of size num_trees_alive that holds the indices of alive trees instead of next_alive_tree

            // indices of the heaps
            std::vector<int> minus_blossoms;
            std::vector<int> plus_empty_edges;
            std::vector<int> plus_plus_internal_edges;
            std::vector<std::vector<std::pair<Tree, int> > > pq_plus_plus;
            std::vector<std::vector<std::pair<Tree, int> > > pq_plus_minus;
            std::vector<std::vector<std::pair<Tree, int> > > pq_minus_plus;
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
        Tree first_alive_tree;

        std::vector<std::pair<int, int> > matching;
        std::vector<std::tuple<int, int, int> > dual_certificate;
        // a vector of (quadrupled dual variable, index of the blossom parent or -1)
        // dual_certificate is empty unless params.compute_dual_certificate is true

        void PrintGraph() const;
        void PrintNode(Node node) const;
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
        Tree MakePrimalUpdate(Tree tree, bool *success);
        void Grow(Tree tree, Edge edge);
        void Shrink(Tree tree, Edge edge_plus_plus);
        void MakeBlossom(std::vector<Edge> blossom_edges, Node lca);
        void Expand(Tree tree, Node blossom);
        Tree Augment(Tree tree, Edge edge);
        void AugmentFromNode(Tree tree, Node node);
        void DissolveTree(Tree tree);
        void ClearNodeDuringTreeDissolve(Tree tree, Node node);
        void Dissolve(Node node);

        void UpdateEdgeAfterShrink(Edge edge, Node node);

        void UpdateNodeInternalTreeStructure(Node node);
        void ClearNodeInternalTreeStructure(Node node);
        void RotateReceptacle(Node blossom, Node child);
        Node DeeperNode(Node blossom, Edge edge) const;

        void UpdateQueuesAfterGrow(Node child, Node grandchild);
        void UpdateQueuesBeforeShrink(const std::vector<Node> &minus_children);
        void RemoveLoopsFromQueues(Node blossom);   // TODO get rid of this
        void UpdateQueuesAfterExpand(Tree tree, Node blossom, const std::vector<Node> &children);
        int TreeTreeQueueIndex(Tree other_tree, const std::vector<std::pair<Tree, int> > &tree_neighbors);
        void AddPQPlusPlus(Tree first, Tree second, Edge edge);
        void AddPQPlusMinus(Tree tree_plus, Tree tree_minus, Edge edge);
        void ValidateQueues();          // debugging purposes
        void ValidatePositiveSlacks();  // debugging purposes

        Edge MinPlusEmptyEdge(int queue_index);
        Edge MinPlusPlusInternalEdge(int queue_index);
        Edge MinPlusPlusExternalEdge(int queue_index);
        Edge MinPlusMinusExternalEdge(int queue_index);
        Node MinMinusBlossom(int queue_index);
        Edge GrowableEdge(Tree tree);
        Edge AugmentableEdge(Tree tree);
        Edge ShrinkableEdge(Tree tree);
        Node ExpandableBlossom(Tree tree);
        int PlusEmptySlack(Tree tree);
        int PlusPlusInternalSlack(Tree tree);
        int MinMinusBlossomVariable(Tree tree);
        std::vector<std::pair<Tree, int> > PlusPlusExternalSlacks(Tree tree);
        std::vector<std::pair<Tree, int> > PlusMinusExternalSlacks(Tree tree);

        bool IsElementary(Node node) const;
        bool IsTopBlossom(Node node) const;
        Node TopBlossom(Node node) const;
        int NodeDepth(Node node) const;
        int DualVariableQuadrupled(Node node) const;

        Node LCA(Node first, Node second) const;
        std::vector<Edge> PathToRoot(Tree tree, Node node);

        int SlackQuadrupled(Edge edge) const;
        Node OtherEnd(Edge edge, Node node) const;
        Node SharedNode(Edge first_edge, Edge second_edge);

        void MakeEdgeMatched(Edge edge);
        void MakeEdgeUnmatched(Edge edge);
        void AddEdgeToQueue(Edge edge, int queue_index);
        void RemoveEdgeFromQueue(Edge edge);
        void AddNodeToQueue(Node node, int queue_index);
        void RemoveNodeFromQueue(Node node);

        Tree NextAliveTree(Tree current);
        Tree FirstAliveTree();
        void UpdateAliveIndices();

        static int InitNumVertices(const std::vector<std::tuple<int, int, int> > &edge_list_);
};

#endif //BLOSSOM_VI_SOLVER_H
