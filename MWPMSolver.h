#ifndef BLOSSOM_VI_MWPMSOLVER_H
#define BLOSSOM_VI_MWPMSOLVER_H

#include <memory>
#include <queue>
#include <deque>
#include <iostream>
#include <boost/container/small_vector.hpp>

#include "DualUpdater.h"

struct SolverParameters {
    bool compute_dual_certificate = false;
    bool verbose = false;
    bool print_statistics = false;
    bool debug = false;
    int init_max_tree_size = 100;
};

class MWPMSolver {
    public:
        int64_t primal_objective;
        int64_t dual_objective;
        const SolverParameters params;

        explicit MWPMSolver(const std::vector<std::tuple<int, int, int> > &edge_list_,
                             const SolverParameters &params_ = {});

        void FindMinPerfectMatching();

        const std::vector<std::pair<int, int> > &Matching() const;
        const std::vector<std::tuple<int, int, int> > &DualCertificate() const;
        // a vector of (quadrupled dual variable, index of the blossom parent or -1)

    private:
        struct Edge {
            int head;
            int tail;
            int slack_quadrupled_amortized_;
            int last_round_updated;

            int queue_index;
            int heap_child;
            int heap_next;
            int heap_prev;

            Edge(int head_, int tail_, int weight_);
        };

        struct ArcIndex {
            int index = -1;

            bool operator>=(const ArcIndex & other) const {
                return index >= other.index;
            }
        };

        struct Node {
            int old_blossom_parent;
            ArcIndex matched_edge;
            ArcIndex minus_parent;
            int receptacle_; // by default, a node is its own receptacle
            // TODO make a labeled UnionFind
            int tree;
            int old_tree;

            int tree_var_at_birth;
            int label;
            int slack_diff;

            bool is_alive;
            bool plus;
            bool old_plus;
            bool is_in_record;

            explicit Node(int index_);
        };

        struct NodeHeapInfo {
            int queue_index;
            int heap_child;
            int heap_next;
            int heap_prev;

            int dual_var_quadrupled_amortized_;
        };

        struct NodeBlossomStructure {
            std::vector<int> blossom_children;
        };

        struct Tree {
            int dual_var_quadrupled;
            int alive_index;
            bool is_alive;

            Tree(int root_);
        };

        struct TreeHeapInfo {
            // indices of the heaps
            int minus_blossoms;
            int plus_empty_edges;
            int plus_plus_internal_edges;
            boost::container::small_vector<std::pair<int, int>, 4> pq_plus_plus;
            boost::container::small_vector<std::pair<int, int>, 4> pq_plus_minus;
        };

        struct EdgeHeap {
            int root = -1;
        };
        struct NodeHeap {
            int root = -1;
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
        std::vector<NodeHeapInfo> node_heap_infos;
        std::vector<int> blossom_parents;
        std::vector<int> blossom_ancestors;
        std::vector<boost::container::small_vector<ArcIndex, 8> > adj_list;
        // TODO store in another way
        std::vector<NodeBlossomStructure> blossom_structures;

        std::vector<Edge> edges;
        std::vector<int> edge_weights;
        std::vector<int> elementary_heads;
        std::vector<int> elementary_tails;
        std::vector<uint8_t> matched;
        std::vector<uint8_t> maybe_has_zero_slack;

        std::vector<Tree> trees;
        std::vector<TreeHeapInfo> tree_heap_infos;
        std::vector<int> alive_trees;
        std::vector<std::deque<int> > tree_nodes;
        std::vector<int> roots; // elementary node

        std::vector<EdgeHeap> edge_heaps;
        std::vector<int8_t> edge_heap_alive;
        std::vector<NodeHeap> node_heaps;

        std::vector<int> primal_update_record;
        std::vector<int> actionable_edges;
        int actionable_edges_head;
        std::vector<int> actionable_nodes;
        int actionable_nodes_head;

        std::vector<ArcIndex> even_path_tmp;
        std::vector<ArcIndex> odd_path_tmp;
        std::vector<int> path_to_root;
        std::vector<int> traversal_nodes_tmp;
        std::vector<int> traversal_lists_tmp;
        std::vector<int> edges_stack_tmp;
        std::vector<int> organize_sizes_tmp;
        std::vector<int> organize_label_to_index_tmp;
        std::vector<std::vector<int> > future_blossoms_tmp;

        std::vector<std::pair<int, int> > matching;
        std::vector<std::tuple<int, int, int> > dual_certificate;
        // a vector of (quadrupled dual variable, index of the blossom parent or -1)
        // dual_certificate is empty unless params.compute_dual_certificate is true

        struct InitComparator {
            bool operator()(const std::pair<int, ArcIndex> & lhs, const std::pair<int, ArcIndex> & rhs) const {
                return lhs.first > rhs.first;
            }
        };
        std::priority_queue<std::pair<int, ArcIndex>, std::vector<std::pair<int, ArcIndex>>, InitComparator> init_plus_empty;
        std::vector<ArcIndex> half_int_clockwise;
        std::vector<int> init_tree_nodes;
        int init_min_plus_plus_edge = -1;
        int init_min_plus_plus_slack_amortized = INT32_MAX;

        void PrintGraph() const;
        void PrintNode(int node) const;

        void Init();
        void InitMakeSlacksNonnegative();
        void InitGreedyIncreaseVars();
        void FractionalMatchingInit();
        void FracInitGrow(ArcIndex arc, int tree_cnt, int tree_var);
        void FracInitAugment(int node_plus, int this_root);
        void FracInitAugmentCycle(ArcIndex arc, int this_root);
        void FracInitMakeCycle(int edge, int this_root);
        void FracInitCleanTree(int tree_var);
        void FracInitCleanup(); // clear minus parents, clear .tree, round the half integral edges


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
        void MakePrimalUpdate(int edge);
        void MakePrimalUpdateForNode(int node);

        void Expand(int blossom);
        void RestoreEdgeEndsBeforeExpand(int blossom);
        void RotateReceptacle(int blossom, int new_receptacle);
        void UpdateInternalStructure(int blossom,
                                     int old_receptacle,
                                     int new_receptacle,
                                     int elder_child);
        void EvenPathToReceptacle(int node);
        void OddPathToReceptacle(int node);
        void ExpandChildBeforeGrow(int blossom);

        void Grow(int parent, ArcIndex arc);

        void MakeCherryBlossom(int edge_plus_plus);
        std::pair<int, int> CherryPathBounds(int first_vertex, int second_vertex);
        void UpdateCherryPath(int lower_node, int upper_node);

        void Augment(int edge_plus_plus);
        void PathToRoot(int node_plus);
        void AugmentPathToRoot();
        void ClearTree(int tree);
        void ClearTreeForLastPair(int tree);

        // updates amortized slack and variables, edge_heaps and node_heaps, old_tree, old_plus,
        // old_blossom_parent
        void UpdateQueuesRecordTraversal();
        void UpdateQueuesFirstPass();
        void UpdateQueuesSecondPass();
        void HandleIncidentEmpty(int node);
        void HandleIncidentPlus(int node);
        void HandleIncidentMinus(int node);
        void UpdateQueuesThirdPass();

        void OrganizeBlossomChildren();
        void Shrink(std::vector<int> &children);

        // for debug:
        auto NodeVariables() const -> std::vector<int>;
        std::vector<int> EdgeSlacks();
        void ValidateEvenOddPaths();
        void ValidateArcs();

        bool MakeDualUpdates();
        void UpdateAliveTreesList();
        std::vector<DualConstraintsNode> GetDualConstraints();
        void InitNextRoundActionable();
        void AddZeroSlackEdgesFromQueue(int queue_index, bool add_to_actionable);
        void CleanLoopsFromQueueTop(int tree); // makes top of plus_plus_internal_edges a non-loop

        bool IsElementary(int node) const {
            return node < num_vertices_elementary;
        }
        int TopBlossom(int node) const {
            while (blossom_parents[node] >= 0) {
                node = blossom_parents[node];
            }
            return node;
        }
        int Receptacle(int node) {
            int grandparent = node;
            while (nodes[grandparent].receptacle_ != grandparent) {
                grandparent = nodes[grandparent].receptacle_;
            }

            int cur_node = node;
            while (nodes[cur_node].receptacle_ != grandparent) {
                int next_node = nodes[cur_node].receptacle_;
                nodes[cur_node].receptacle_ = grandparent;
                cur_node = next_node;
            }

            return grandparent;
        }
        int DualVariableQuadrupled(int node) const;
        int DualVariableQuadrupled(int node, int tree, bool plus, int blossom_parent) const;
        void UpdateNonLoopNeighbors(int node);

        int SlackQuadrupled(int edge);
        int PlusPlusLCA(int first_vertex, int second_vertex);

        void MakeEdgeMatched(int edge);
        void MakeEdgeUnmatched(int edge);

        void AddEdgeToThisQueue(int edge, int queue_index);
        void RemoveEdgeFromQueue(int edge);
        void AddNodeToQueue(int node, int queue_index);
        void RemoveNodeFromQueue(int node);
        int TreeTreeQueueIndex(int other_tree,
                               boost::container::small_vector<std::pair<int, int>, 4> *tree_neighbors) const;

        int MinPlusPlusInternalEdge(int queue_index);
        int PopExpandableBlossom(int tree);

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

        int OtherEnd(ArcIndex arc) {
            if (arc.index % 2 == 1) {
                return Head(arc.index >> 1);
            }
            return Tail(arc.index >> 1);
        }
        int ThisEnd(ArcIndex arc) {
            if (arc.index % 2 == 0) {
                return Head(arc.index >> 1);
            }
            return Tail(arc.index >> 1);
        }
        ArcIndex ReverseArc(ArcIndex arc) {
            if (arc.index % 2 == 0) {
                return ArcIndex(arc.index + 1);
            }
            return ArcIndex(arc.index - 1);
        }
        int OtherElementaryEnd(ArcIndex arc) const {
            if (arc.index % 2 == 1) {
                return elementary_heads[arc.index >> 1];
            }
            return elementary_tails[arc.index >> 1];
        }
        int ThisElementaryEnd(ArcIndex arc) const {
            if (arc.index % 2 == 0) {
                return elementary_heads[arc.index >> 1];
            }
            return elementary_tails[arc.index >> 1];
        }
        int Head(int edge) {
            int head = edges[edge].head;

            if (blossom_ancestors[head] < 0) {
                return head;
            }

            int top_blossom = head;
            while (blossom_ancestors[top_blossom] >= 0) {
                top_blossom = blossom_ancestors[top_blossom];
            }
            int next_head = blossom_ancestors[head];
            while (next_head != top_blossom) {
                blossom_ancestors[head] = top_blossom;
                head = next_head;
                next_head = blossom_ancestors[next_head];
            }
            edges[edge].head = top_blossom;

            return top_blossom;
        }
        int Tail(int edge) {
            // TODO code duplication
            int tail = edges[edge].tail;

            if (blossom_ancestors[tail] < 0) {
                return tail;
            }

            int top_blossom = tail;
            while (blossom_ancestors[top_blossom] >= 0) {
                top_blossom = blossom_ancestors[top_blossom];
            }
            int next_tail = blossom_ancestors[tail];
            while (next_tail != top_blossom) {
                blossom_ancestors[tail] = top_blossom;
                tail = next_tail;
                next_tail = blossom_ancestors[next_tail];
            }
            edges[edge].tail = top_blossom;

            return top_blossom;
        }
};

#endif //BLOSSOM_VI_MWPMSOLVER_H
