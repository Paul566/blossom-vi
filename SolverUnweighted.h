#include <vector>
#include <memory>
#include <iostream>
#include <queue>
#include "Edge.h"


class SolverUnweighted {
public:
    std::vector<std::vector<std::shared_ptr<Edge>>> adj_list;

    SolverUnweighted(const std::vector<std::vector<int>> &adj_list_);

    void Print();

    void Solve();

private:
    const int n;  // the number of vertices
    std::vector<std::shared_ptr<Edge>> matched_edge;

    void GreedyInit();

    void Augment(std::vector<std::shared_ptr<Edge>> path, int start);

    // first_vertex and second_vertex are pluses
    // walks up the tree until the lca and makes everyone plusminus, 
    // adding the new pluses to the queue and updating parent pointers
    void UpdateLabels(std::shared_ptr<Edge> edge_plus_plus, int root, 
        std::queue<int> & queue, 
        std::vector<std::shared_ptr<Edge>> & plus_parents,
        std::vector<std::shared_ptr<Edge>> & minus_parents,
        std::vector<bool> & plus,
        std::vector<bool> & minus,
        std::vector<int> & depth_plus,
        std::vector<int> & depth_minus);
};