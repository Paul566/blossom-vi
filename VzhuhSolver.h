#ifndef BLOSSOM_VI_VZHUHSOLVER_H
#define BLOSSOM_VI_VZHUHSOLVER_H

#include <memory>
#include <queue>
#include <deque>
#include <stack>

#include "Heap.h"

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
            Heap<int>::Handle *handle;
            int queue_index;
            int weight;
            int slack_quadrupled_amortized_;
            int slack_diff;
            int head;
            int tail;
            int elementary_head;
            int elementary_tail;
            bool matched;
            bool maybe_has_zero_slack;
            bool must_be_updated;
            bool maybe_was_loop;
            Edge(int head_, int tail_, int weight_);
        };

        struct Node {
            Heap<int>::Handle *handle;
            std::vector<int> blossom_children;
            std::vector<int> neighbors;   // TODO make sure we don't use too much memory
            std::vector<int> zero_slack_neighbors;

            int blossom_parent;
            int old_blossom_parent;
            int matched_edge;
            int minus_parent;
            int receptacle_; // by default, a node is its own receptacle
            int tree;
            int old_tree;

            int queue_index;
            int dual_var_quadrupled_amortized_;
            int tree_var_at_birth;
            int label;
            int round_0slack_neighbors_updated;

            bool is_alive;
            bool plus;
            bool old_plus;
            bool is_in_record;

            explicit Node(int index);
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
            const int minus_blossoms;
            const int plus_empty_edges;
            const int plus_plus_internal_edges;

            bool is_alive;

            Tree(int root_, int minus_blossoms_, int plus_empty_edges_, int plus_plus_internal_edges_);
        };

        struct DualConstraints {
            std::vector<int> upper_bound;
            // upper bound on delta_T quadrupled
            std::vector<std::vector<std::pair<int, int> > > plus_plus_constraints;
            // (other_tree_index_alive, slack_plus_plus quadrupled)
            std::vector<std::vector<std::pair<int, int> > > plus_minus_constraints;
            // (other_tree_index_alive, slack_plus_minus quadrupled)
        };

        struct PrimalUpdateRecord {
            std::vector<int> changed_sign;
        };

        const int num_vertices_elementary; // the original number of vertices (i.e. not counting blossoms)
        int num_trees_alive;
        int current_round;
        int aux_counter1;
        int aux_counter2;
        int aux_counter3;
        int aux_counter4;

        std::vector<Node> nodes;
        std::vector<Edge> edges;
        std::vector<Tree> trees;
        std::vector<int> alive_trees;
        std::vector<std::unique_ptr<Heap<int>> > node_heaps;
        std::vector<std::unique_ptr<Heap<int>> > edge_heaps;

        std::vector<std::vector<int> > adj_list;
        // std::vector<std::vector<int> > zero_slack_adj_list;
        std::queue<int> actionable_edges;
        std::queue<int> actionable_nodes;

        int nodes_label_cnt; // TODO make int64_t

        std::vector<std::pair<int, int> > matching;
        std::vector<std::tuple<int, int, int> > dual_certificate;
        // a vector of (quadrupled dual variable, index of the blossom parent or -1)
        // dual_certificate is empty unless params.compute_dual_certificate is true

        void PrintGraph() const ;
        void PrintNode(int node) const ;

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
        void MakePrimalUpdate(int edge, PrimalUpdateRecord *record);
        void MakePrimalUpdateForNode(int node, PrimalUpdateRecord *record);

        void Expand(int blossom, PrimalUpdateRecord *record);
        void RestoreEdgeEndsBeforeExpand(int blossom);
        void RotateReceptacle(int blossom, int new_receptacle);
        void UpdateInternalStructure(int blossom,
                                     int old_receptacle,
                                     int new_receptacle,
                                     int elder_child,
                                     PrimalUpdateRecord *record);
        std::vector<int> EvenPathToReceptacle(int node);
        std::vector<int> OddPathToReceptacle(int node);
        void ExpandChildBeforeGrow(int blossom, PrimalUpdateRecord *record);

        void Grow(int parent, int edge, PrimalUpdateRecord *record);

        void MakeCherryBlossom(int edge_plus_plus, PrimalUpdateRecord *record);
        std::pair<int, int> CherryPathBounds(int first_vertex, int second_vertex);
        void UpdateCherryPath(int lower_node, int upper_node, PrimalUpdateRecord *record);

        void Augment(int edge_plus_plus, PrimalUpdateRecord *record);
        std::vector<int> PathToRoot(int node_plus);
        void AugmentPath(const std::vector<int> &path);
        void ClearTree(int tree, PrimalUpdateRecord *record);

        void UpdateQueues(const PrimalUpdateRecord &record);
        // updates amortized slam and variables, edge_heaps and node_heaps, old_tree, old_plus, old_blossom_parent
        void UpdateEdgeSlack(int edge);
        std::vector<std::vector<int> > OrganizeBlossomChildren(const PrimalUpdateRecord &record);
        void Shrink(std::vector<int> &children);

        // for debug:
        void ValidateQueues();
        auto NodeVariables() const -> std::vector<int>;
        std::vector<int> EdgeSlacks();
        void ValidateEvenOddPaths();

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
        std::vector<int>& NonLoopNeighbors(int node);
        std::vector<int>& NonLoopZeroSlackNeighbors(int node);
        std::vector<int> ElementaryBlossomDescendants(int node) const;

        int SlackQuadrupled(int edge);
        int OldSlackQuadrupled(int edge);
        int OtherEnd(int edge, int node);
        int OtherElementaryEnd(int edge, int node) const;
        int Head(int edge);
        int Tail(int edge);
        int PlusPlusLCA(int first_vertex, int second_vertex);

        void MakeEdgeMatched(int edge);
        void MakeEdgeUnmatched(int edge);

        void AddEdgeToQueue(int edge);
        void AddEdgeToThisQueue(int edge, int queue_index);
        void RemoveEdgeFromQueue(int edge);
        void AddNodeToQueue(int node, int queue_index);
        void RemoveNodeFromQueue(int node);
        void AddPQPlusPlus(int first, int second, int edge);
        void AddPQPlusMinus(int tree_plus, int tree_minus, int edge);
        int TreeTreeQueueIndex(int other_tree, std::vector<std::pair<int, int> > *tree_neighbors) const;

        int MinPlusEmptyEdge(int queue_index) const;
        int MinPlusPlusInternalEdge(int queue_index);
        int MinPlusPlusExternalEdge(int queue_index) const;
        int MinPlusMinusExternalEdge(int queue_index) const;
        int MinMinusBlossom(int queue_index) const;
        int PopExpandableBlossom(int tree) const;
        int PlusEmptySlack(int tree);
        int PlusPlusInternalSlack(int tree);
        int MinMinusBlossomVariable(int tree) const;
        std::vector<std::pair<int, int> > PlusPlusExternalSlacks(int tree);
        std::vector<std::pair<int, int> > PlusMinusExternalSlacks(int tree);

        void AddNodeToRecord(int node, PrimalUpdateRecord *record);

        static int InitNumVertices(const std::vector<std::tuple<int, int, int> > &edge_list_);
};

#endif //BLOSSOM_VI_VZHUHSOLVER_H
