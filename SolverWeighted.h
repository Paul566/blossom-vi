#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <memory>
#include "EdgeWeighted.h"

class SolverWeighted {
public:
    explicit SolverWeighted(const std::vector<std::tuple<int, int, int>> &edge_list_);

    void FindMinPerfectMatching();

    void PrintElementaryAdjList();

    void GreedyInit();

    std::vector<std::shared_ptr<Node>> Roots();

    std::vector<std::pair<int, int>> Matching();

private:
    std::vector<std::shared_ptr<Node>> elementary_nodes;

    int DualObjectiveQuadrupled();

    int PrimalObjectiveQuadrupled() const;

    void DissolveBlossoms();

    static void DestroyBlossomTopDown(const std::shared_ptr<Node> &blossom);
};

#endif //SOLVER_H
