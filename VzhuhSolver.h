#ifndef BLOSSOM_VI_VZHUHSOLVER_H
#define BLOSSOM_VI_VZHUHSOLVER_H

#include <memory>
#include <queue>
#include <deque>
#include <iostream>
#include <stack>

struct SolverParameters {
    bool compute_dual_certificate = false;
    bool verbose = false;
    bool print_statistics = true;
    bool debug = false;
    int heap_arity = 4;
};

class VzhuhSolver {
    public:
        int64_t primal_objective;
        int64_t dual_objective;
        const SolverParameters params;

        explicit VzhuhSolver(const std::vector<std::tuple<int, int, int> > &edge_list_,
                             const SolverParameters &params_ = {});

        void FindMinPerfectMatching();

        const std::vector<std::pair<int, int> > &Matching() const;
        const std::vector<std::tuple<int, int, int> > &DualCertificate() const;
        // a vector of (quadrupled dual variable, index of the blossom parent or -1)

    private:
        struct Edge {
            int queue_index;
            int heap_child;
            int heap_next;
            int heap_prev;

            // const int index;
            int weight;
            int slack_quadrupled_amortized_;
            int slack_diff;
            int head;
            int tail;
            int elementary_head;
            int elementary_tail;
            int last_round_updated;
            bool matched;
            bool maybe_has_zero_slack;
            bool must_be_updated;
            bool maybe_was_loop;
            Edge(int head_, int tail_, int weight_);
            int Key() const;
        };

        struct ArcIndex {
            int index = -1;
        };

        struct Node {
            int queue_index;
            int heap_child;
            int heap_next;
            int heap_prev;

            std::vector<int> blossom_children;
            std::vector<ArcIndex> neighbors; // TODO make sure we don't use too much memory

            int blossom_parent;
            int old_blossom_parent;
            ArcIndex matched_edge;
            ArcIndex minus_parent;
            int receptacle_; // by default, a node is its own receptacle
            int tree;
            int old_tree;

            int dual_var_quadrupled_amortized_;
            int tree_var_at_birth;
            int label;
            int slack_diff;

            bool is_alive;
            bool plus;
            bool old_plus;
            bool is_in_record;

            explicit Node(int index_);
            int Key() const;
        };

        struct Tree {
            std::vector<std::pair<int, int> > pq_plus_plus;
            std::vector<std::pair<int, int> > pq_plus_minus;
            std::vector<std::pair<int, int> > pq_minus_plus;

            std::deque<int> tree_nodes;

            const int root; // elementary node
            int dual_var_quadrupled;
            int alive_index;

            // indices of the heaps
            int minus_blossoms;
            int plus_empty_edges;
            int plus_plus_internal_edges;

            bool is_alive;

            Tree(int root_, int minus_blossoms_, int plus_empty_edges_, int plus_plus_internal_edges_);
        };

        struct EdgeHeap {
            int root = -1;
        };
        struct NodeHeap {
            int root = -1;
        };

        struct DualConstraints {
            std::vector<int> upper_bound;
            // upper bound on delta_T quadrupled
            std::vector<std::vector<std::pair<int, int> > > plus_plus_constraints;
            // (other_tree_index_alive, slack_plus_plus quadrupled)
            std::vector<std::vector<std::pair<int, int> > > plus_minus_constraints;
            // (other_tree_index_alive, slack_plus_minus quadrupled)
        };

        const int num_vertices_elementary; // the original number of vertices (i.e. not counting blossoms)
        int num_trees_alive;
        int current_round;
        int aux_counter1;
        int aux_counter2;
        int aux_counter3;
        int aux_counter4;
        int nodes_label_cnt; // TODO make int64_t

        std::vector<Node> nodes;
        std::vector<Edge> edges;
        std::vector<Tree> trees;
        std::vector<int> alive_trees;
        std::vector<EdgeHeap> edge_heaps;
        std::vector<NodeHeap> node_heaps;

        std::vector<int> primal_update_record;

        std::queue<int> actionable_edges;
        std::queue<int> actionable_nodes;

        std::vector<std::pair<int, int> > matching;
        std::vector<std::tuple<int, int, int> > dual_certificate;
        // a vector of (quadrupled dual variable, index of the blossom parent or -1)
        // dual_certificate is empty unless params.compute_dual_certificate is true

        void PrintGraph() const;
        void PrintNode(int node) const;

        void GreedyInit();
        void InitializeTrees();

        void ComputeMatching();
        void ComputePrimalObjective();
        void UpdateMatching(int blossom, int new_receptacle);
        void ComputeDualCertificate();
        void ComputeDualObjectiveQuadrupled();

        void DestroyBlossoms();
        void RestoreFinalEdgeEnds();
        int FindFinalReceptacle(int blossom) const;

        bool MakePrimalUpdates();
        // first phase: expand
        // second phase: grow, make cherry blossoms, augment
        // third phase: shrink the cherry blossoms
        // fourth phase: update the queues
        // TODO make record a field, don't make a new vector every time
        void MakePrimalUpdate(int edge);
        void MakePrimalUpdateForNode(int node);

        void Expand(int blossom);
        void RestoreEdgeEndsBeforeExpand(int blossom);
        void RotateReceptacle(int blossom, int new_receptacle);
        void UpdateInternalStructure(int blossom,
                                     int old_receptacle,
                                     int new_receptacle,
                                     int elder_child);
        std::vector<ArcIndex> EvenPathToReceptacle(int node);
        std::vector<ArcIndex> OddPathToReceptacle(int node);
        void ExpandChildBeforeGrow(int blossom);

        void Grow(int parent, ArcIndex arc);

        void MakeCherryBlossom(int edge_plus_plus);
        std::pair<int, int> CherryPathBounds(int first_vertex, int second_vertex);
        void UpdateCherryPath(int lower_node, int upper_node);

        void Augment(int edge_plus_plus);
        std::vector<int> PathToRoot(int node_plus);
        void AugmentPath(const std::vector<int> &path);
        void ClearTree(int tree);

        // updates amortized slam and variables, edge_heaps and node_heaps, old_tree, old_plus, old_blossom_parent
        void UpdateQueuesRecordTraversal();
        void UpdateEdgeInfo(int edge, int endpoint, int other_endpoint, int queue_index);
        std::vector<std::vector<int> > OrganizeBlossomChildren();
        void Shrink(std::vector<int> &children);

        // for debug:
        void ValidateQueues();
        auto NodeVariables() const -> std::vector<int>;
        std::vector<int> EdgeSlacks();
        void ValidateEvenOddPaths();
        void ValidateArcs();

        bool MakeDualUpdates();
        void UpdateAliveTreesList();
        std::vector<int> VariableDeltas();
        std::vector<std::vector<int> > ConnectedComponentsTreeTree(const DualConstraints &dual_constraints) const;
        DualConstraints GetDualConstraints();
        void InitNextRoundActionable();
        void AddZeroSlackEdgesFromQueue(int queue_index, bool add_to_actionable);
        void CleanLoopsFromQueueTop(int tree); // makes top of plus_plus_internal_edges a non-loop

        bool IsElementary(int node) const;
        int TopBlossom(int node) const;
        int Receptacle(int node);
        int DualVariableQuadrupled(int node) const;
        int DualVariableQuadrupled(int node, int tree, bool plus, int blossom_parent) const;
        std::vector<ArcIndex> &NonLoopNeighbors(int node);
        std::vector<int> ElementaryBlossomDescendants(int node) const;

        int SlackQuadrupled(int edge);
        int OldSlackQuadrupled(int edge);
        int OtherEnd(ArcIndex arc);
        int ThisEnd(ArcIndex arc);
        ArcIndex ReverseArc(ArcIndex arc);
        int OtherElementaryEnd(ArcIndex arc) const;
        int Head(int edge);
        int Tail(int edge);
        int PlusPlusLCA(int first_vertex, int second_vertex);

        void MakeEdgeMatched(int edge);
        void MakeEdgeUnmatched(int edge);

        void AddEdgeToThisQueue(int edge, int queue_index);
        void RemoveEdgeFromQueue(int edge);
        void AddNodeToQueue(int node, int queue_index);
        void RemoveNodeFromQueue(int node);
        int TreeTreeQueueIndex(int other_tree, std::vector<std::pair<int, int> > *tree_neighbors) const;

        int MinPlusPlusInternalEdge(int queue_index);
        int PopExpandableBlossom(int tree);
        int PlusEmptySlack(int tree);
        int PlusPlusInternalSlack(int tree);
        int MinMinusBlossomVariable(int tree) const;
        std::vector<std::pair<int, int> > PlusPlusExternalSlacks(int tree);
        std::vector<std::pair<int, int> > PlusMinusExternalSlacks(int tree);

        void AddNodeToRecord(int node);

        static int InitNumVertices(const std::vector<std::tuple<int, int, int> > &edge_list_);

        int GetMinEdgeHeap(int heap_index) const;
        void InsertEdgeHeap(int edge, int heap_index);
        void RemoveMinEdgeHeap(int heap_index);
        void RemoveEdgeHeap(int edge);
        int MeldEdgeHeap(int edge_a, int edge_b);
        int TwoPassMergeEdgeHeap(int edge_first);
        void CutEdgeHeap(int edge);
        // TODO avoid duplication
        int GetMinNodeHeap(int heap_index) const;
        void InsertNodeHeap(int node, int heap_index);
        void RemoveMinNodeHeap(int heap_index);
        void RemoveNodeHeap(int node);
        int MeldNodeHeap(int node_a, int node_b);
        int TwoPassMergeNodeHeap(int node_first);
        void CutNodeHeap(int node);
};

#endif //BLOSSOM_VI_VZHUHSOLVER_H
