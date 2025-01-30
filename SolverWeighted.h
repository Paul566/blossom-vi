#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <memory>
#include "EdgeWeighted.h"

class SolverWeighted {
public:
    explicit SolverWeighted(const std::vector<std::tuple<int, int, int>> &edge_list_);

private:
    std::vector<std::shared_ptr<Node>> elementary_nodes;
};



#endif //SOLVER_H
