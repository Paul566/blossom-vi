#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <memory>
#include "EdgeWeighted.h"

class SolverWeighted {
public:
    explicit SolverWeighted(const std::vector<std::tuple<int, int, int>> &edge_list_);

    void PrintElementaryAdjList();

    void GreedyInit();

    std::vector<std::shared_ptr<Node>> Roots();

private:
    std::vector<std::shared_ptr<Node>> elementary_nodes;

    std::shared_ptr<Node> TopBlossom(std::shared_ptr<Node> vertex);
};

#endif //SOLVER_H
