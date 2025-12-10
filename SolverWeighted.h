#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <list>
#include <unordered_map>

#include "EdgeWeighted.h"
#include "Tree.h"

struct SolverParams {
    bool compute_dual_certificate = false;
    bool verbose = false;
    bool multiple_trees = true;
};

struct DualConstraints {
    std::vector<int> upper_bound;
    // upper bound on delta_T quadrupled
    std::vector<std::vector<std::pair<int, int> > > plus_plus_constraints;
    // (other_tree_index, slack_plus_plus quadrupled)
    std::vector<std::vector<std::pair<int, int> > > plus_minus_constraints;
    // (other_tree_index, slack_plus_minus quadrupled)
};

class SolverWeighted {
    public:
        int64_t primal_objective;
        int64_t dual_objective;
        const SolverParams params;

        explicit SolverWeighted(const std::vector<std::tuple<int, int, int> > &edge_list_,
                                const SolverParams &params_ = {});

        void FindMinPerfectMatching();

        void PrintGraph() const;

        std::vector<std::pair<int, int> > Matching() const;
        const std::vector<std::tuple<int, int, int>>& DualCertificate() const;
        // (index, quadrupled dual variable, index of the blossom parent or -1)

    private:
        std::list<Node> elementary_nodes_list;
        std::vector<std::list<Node>::iterator> elementary_iters;
        std::list<EdgeWeighted> edges;
        std::list<Node> blossoms;
        std::unordered_map<Node *, std::list<Node>::iterator> iter_to_blossom;
        std::list<Tree> trees;
        std::unordered_map<Tree *, std::list<Tree>::iterator> iter_to_tree;
        // TODO do lazy deletions from the list
        // TODO store the iterators in some other way
        // TODO consider std::vector<std::unique_ptr<Obj>> instead of std::list<Obj>

        std::vector<std::tuple<int, int, int>> dual_certificate;
        // (index, quadrupled dual variable, index of the blossom parent or -1)
        // dual_certificate is empty unless params.compute_dual_certificate is true

        void GreedyInit();

        void MakeDualUpdates();
        void MakePrimalUpdates();

        int OptimalSingleDelta();
        std::vector<int> VariableDeltas();
        std::vector<std::vector<int>> ConnectedComponentsTreeTree(const DualConstraints & dual_constraints);
        DualConstraints GetDualConstraints();

        void DestroyBlossoms();

        int64_t DualObjectiveQuadrupled() const;
        int64_t PrimalObjective() const;

        void ComputeDualCertificate();
};

#endif //SOLVER_H
