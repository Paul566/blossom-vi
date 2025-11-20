#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <list>
#include <unordered_map>

#include "EdgeWeighted.h"
#include "Tree.h"

class SolverWeighted {
    public:
        explicit SolverWeighted(const std::vector<std::tuple<int, int, int> > &edge_list_);

        void FindMinPerfectMatching();

        void PrintElementaryAdjList() const;

        void GreedyInit();

        std::vector<int> RootIndices() const;

        std::vector<std::pair<int, int> > Matching() const;

        int DualObjectiveQuadrupled() const;

        int PrimalObjectiveQuadrupled() const;

    private:
        std::list<Node> elementary_nodes_list; // never reallocated after initialization
        std::vector<std::list<Node>::iterator> elementary_iters;
        std::list<EdgeWeighted> edges;    // never reallocated after initialization
        std::list<Node> blossoms;
        std::list<Tree> trees;
        std::unordered_map<Node *, std::list<Node>::iterator> iter_to_self;
        // TODO store the iterators in some other way
        // TODO maybe make all nodes live inside a single std::list

        void DestroyBlossoms();
};

#endif //SOLVER_H
