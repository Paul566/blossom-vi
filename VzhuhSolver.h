#ifndef BLOSSOM_VI_VZHUHSOLVER_H
#define BLOSSOM_VI_VZHUHSOLVER_H

#include <cstdint>
#include <memory>
#include <queue>

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
        struct Index {
            int index;
            explicit Index(int index_) : index(index_) {
            }
            explicit operator bool() const {
                return index >= 0;
            }
            bool operator==(const Index &other) const {
                return index == other.index;
            }
            Index &operator++() {
                ++index;
                return *this;
            }
            Index &operator--() {
                --index;
                return *this;
            }
            friend bool operator<(Index lhs, std::size_t rhs) {
                return lhs.index < rhs;
            }
            friend bool operator<(Index lhs, Index rhs) {
                return lhs.index < rhs.index;
            }
            friend std::ostream &operator<<(std::ostream &os, const Index &idx) {
                return os << idx.index;
            }
        };
        struct NodeIndex : public Index {
            explicit NodeIndex(int index_) : Index(index_) {
            }
        };
        struct EdgeIndex : public Index {
            explicit EdgeIndex(int index_) : Index(index_) {
            }
        };
        struct TreeIndex : public Index {
            explicit TreeIndex(int index_) : Index(index_) {
            }
        };

        struct NodeComparator {
            const VzhuhSolver *solver;
            bool operator()(const NodeIndex &a, const NodeIndex &b) const;
        };
        struct EdgeComparator {
            const VzhuhSolver *solver;
            bool operator()(const EdgeIndex &a, const EdgeIndex &b) const;
        };
        using EdgeHeap = Heap<EdgeIndex, EdgeComparator>;
        using NodeHeap = Heap<NodeIndex, NodeComparator>;

        struct Edge {
            EdgeHeap::Handle *handle;
            int queue_index;
            int weight;
            int slam_quadrupled_amortized_;
            NodeIndex head;
            NodeIndex tail;
            NodeIndex elementary_head;
            NodeIndex elementary_tail;
            bool matched;
            bool is_in_zero_slack_set;
            Edge(int head_, int tail_, int weight_);
        };

        struct Node {
            bool is_alive;
            int dual_var_quadrupled_amortized_;
            EdgeIndex matched_edge;

            // blossom related fields
            NodeIndex blossom_parent;
            std::vector<NodeIndex> blossom_children;

            // tree related fields
            bool plus;
            EdgeIndex minus_parent;
            NodeIndex receptacle_; // by default, a node is its own receptacle
            TreeIndex tree;
            bool old_plus;
            TreeIndex old_tree;
            NodeIndex old_blossom_parent;

            // queue related fields
            int queue_index;
            NodeHeap::Handle *handle;

            int label;

            explicit Node(int index);
        };

        struct Tree {
            bool is_alive;
            const NodeIndex root; // elementary node
            int dual_var_quadrupled;
            int alive_index;
            std::vector<NodeIndex> tree_nodes;

            // indices of the heaps
            const int minus_blossoms;
            const int plus_empty_edges;
            const int plus_plus_internal_edges;
            std::vector<std::pair<TreeIndex, int> > pq_plus_plus;
            std::vector<std::pair<TreeIndex, int> > pq_plus_minus;
            std::vector<std::pair<TreeIndex, int> > pq_minus_plus;

            Tree(int root_, int minus_blossoms_, int plus_empty_edges_, int plus_plus_internal_edges_);
        };

        template <class Obj, class Idx>
        class Vector {
            public:
                Obj &operator[](const Idx idx) {
                    return obj_vector[idx.index];
                }
                const Obj &operator[](const Idx idx) const {
                    return obj_vector[idx.index];
                }
                Obj &Back() {
                    return obj_vector.back();
                }
                std::size_t Size() const {
                    return obj_vector.size();
                }
                void PushBack(const Obj &obj) {
                    obj_vector.push_back(obj);
                }
                template <typename... Args>
                void EmplaceBack(Args&&... args) {
                    obj_vector.emplace_back(std::forward<Args>(args)...);
                }
                void Reserve(std::size_t size) {
                    obj_vector.reserve(size);
                }

            private:
                std::vector<Obj> obj_vector;
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
            std::vector<NodeIndex> changed_sign;
        };

        const int num_vertices_elementary; // the original number of vertices (i.e. not counting blossoms)
        int num_trees_alive;

        Vector<Node, NodeIndex> nodes;
        Vector<Edge, EdgeIndex> edges;
        Vector<Tree, TreeIndex> trees;
        std::vector<TreeIndex> alive_trees;
        std::vector<std::unique_ptr<NodeHeap> > node_heaps;
        std::vector<std::unique_ptr<EdgeHeap> > edge_heaps;

        std::vector<std::vector<EdgeIndex> > adj_list;
        std::vector<std::vector<EdgeIndex> > zero_slack_adj_list;
        // TODO zero_slack_adj_list might also decrease
        std::queue<EdgeIndex> actionable_edges;

        int nodes_label_cnt; // TODO make int64_t

        std::vector<std::pair<int, int> > matching;
        std::vector<std::tuple<int, int, int> > dual_certificate;
        // a vector of (quadrupled dual variable, index of the blossom parent or -1)
        // dual_certificate is empty unless params.compute_dual_certificate is true

        void PrintGraph();
        void PrintNode(NodeIndex node);

        void GreedyInit();
        void InitializeTrees();

        void ComputeMatching();
        void ComputePrimalObjective();
        void DestroyBlossoms();
        void UpdateMatching(NodeIndex blossom, NodeIndex new_receptacle);
        void ComputeDualCertificate();
        void ComputeDualObjectiveQuadrupled();

        bool MakePrimalUpdates();
        // first phase: expand
        // second phase: grow, make cherry blossoms, augment
        // third phase: shrink the cherry blossoms
        // fourth phase: update the queues
        void MakePrimalUpdate(EdgeIndex edge, PrimalUpdateRecord *record);

        void Expand(NodeIndex blossom, PrimalUpdateRecord *record);
        void RestoreEdgeEndsBeforeExpand(NodeIndex blossom);
        void RotateReceptacle(NodeIndex blossom, NodeIndex new_receptacle);
        void UpdateInternalStructure(NodeIndex blossom,
                                     NodeIndex old_receptacle,
                                     NodeIndex new_receptacle,
                                     NodeIndex elder_child,
                                     PrimalUpdateRecord *record);
        std::vector<EdgeIndex> EvenPathToReceptacle(NodeIndex node);
        std::vector<EdgeIndex> OddPathToReceptacle(NodeIndex node);
        void ExpandChildBeforeGrow(NodeIndex blossom, PrimalUpdateRecord *record);

        void Grow(NodeIndex parent, EdgeIndex edge, PrimalUpdateRecord *record);
        void AddNeighborsToActionable(NodeIndex node);

        void MakeCherryBlossom(EdgeIndex edge_plus_plus, PrimalUpdateRecord *record);
        std::pair<NodeIndex, NodeIndex> CherryPathBounds(NodeIndex first_vertex, NodeIndex second_vertex);
        void UpdateCherryPath(NodeIndex lower_node, NodeIndex upper_node, PrimalUpdateRecord *record);

        void Augment(EdgeIndex edge_plus_plus, PrimalUpdateRecord *record);
        std::vector<EdgeIndex> PathToRoot(NodeIndex node_plus);
        void AugmentPath(const std::vector<EdgeIndex> &path);
        void ClearTree(TreeIndex tree, PrimalUpdateRecord *record);

        void UpdateQueues(const PrimalUpdateRecord &record);
        // updates amortized slam and variables, edge_heaps and node_heaps, old_tree, old_plus, old_blossom_parent
        std::vector<std::vector<NodeIndex> > OrganizeBlossomChildren(const PrimalUpdateRecord &record);
        void Shrink(std::vector<NodeIndex> &children);

        // for debug:
        void ValidateQueues();
        auto NodeVariables() const -> std::vector<int>;
        std::vector<int> EdgeSlams();
        void ValidateZeroSlackAdjList();
        void ValidateEvenOddPaths();

        bool MakeDualUpdates();
        void UpdateAliveTreesList();
        std::vector<int> VariableDeltas();
        std::vector<std::vector<int> > ConnectedComponentsTreeTree(const DualConstraints &dual_constraints) const;
        DualConstraints GetDualConstraints();
        void UpdateZeroSlackSetAndActionable();
        void AddZeroSlackEdgesFromQueue(int queue_index, bool add_to_actionable);
        void CleanLoopsFromQueueTop(TreeIndex tree); // makes top of plus_plus_internal_edges a non-loop

        bool IsElementary(NodeIndex node) const;
        NodeIndex TopBlossom(NodeIndex node) const;
        NodeIndex Receptacle(NodeIndex node);
        int DualVariableQuadrupled(NodeIndex node) const;
        int DualVariableQuadrupled(NodeIndex node, TreeIndex tree, bool plus, NodeIndex blossom_parent) const;
        std::vector<EdgeIndex> NonLoopNeighbors(NodeIndex node);
        std::vector<EdgeIndex> NeighborsWLoops(NodeIndex node);
        std::vector<EdgeIndex> NonLoopZeroSlackNeighbors(NodeIndex node);
        std::vector<NodeIndex> ElementaryBlossomDescendants(NodeIndex node) const;

        int SlamQuadrupled(EdgeIndex edge);
        int OldSlamQuadrupled(EdgeIndex edge);
        NodeIndex OtherEnd(EdgeIndex edge, NodeIndex node);
        NodeIndex OtherElementaryEnd(EdgeIndex edge, NodeIndex node) const;
        NodeIndex Head(EdgeIndex edge);
        NodeIndex Tail(EdgeIndex edge);
        NodeIndex PlusPlusLCA(NodeIndex first_vertex, NodeIndex second_vertex);

        void MakeEdgeMatched(EdgeIndex edge);
        void MakeEdgeUnmatched(EdgeIndex edge);

        void AddEdgeToQueue(EdgeIndex edge);
        void AddEdgeToThisQueue(EdgeIndex edge, int queue_index);
        void RemoveEdgeFromQueue(EdgeIndex edge);
        void AddNodeToQueue(NodeIndex node, int queue_index);
        void RemoveNodeFromQueue(NodeIndex node);
        void AddPQPlusPlus(TreeIndex first, TreeIndex second, EdgeIndex edge);
        void AddPQPlusMinus(TreeIndex tree_plus, TreeIndex tree_minus, EdgeIndex edge);
        int TreeTreeQueueIndex(TreeIndex other_tree, std::vector<std::pair<TreeIndex, int> > *tree_neighbors);

        EdgeIndex MinPlusEmptyEdge(int queue_index);
        EdgeIndex MinPlusPlusInternalEdge(int queue_index);
        EdgeIndex MinPlusPlusExternalEdge(int queue_index);
        EdgeIndex MinPlusMinusExternalEdge(int queue_index);
        NodeIndex MinMinusBlossom(int queue_index) const;
        NodeIndex PopExpandableBlossom(TreeIndex tree);
        int PlusEmptySlack(TreeIndex tree);
        int PlusPlusInternalSlack(TreeIndex tree);
        int MinMinusBlossomVariable(TreeIndex tree);
        std::vector<std::pair<TreeIndex, int> > PlusPlusExternalSlacks(TreeIndex tree);
        std::vector<std::pair<TreeIndex, int> > PlusMinusExternalSlacks(TreeIndex tree);

        void DeleteDuplicates(std::vector<NodeIndex> *node_list);

        static int InitNumVertices(const std::vector<std::tuple<int, int, int> > &edge_list_);
};

#endif //BLOSSOM_VI_VZHUHSOLVER_H
