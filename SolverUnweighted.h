#ifndef BLOSSOM_VI_SOLVERUNWEIGHTED_H
#define BLOSSOM_VI_SOLVERUNWEIGHTED_H


#include <vector>
#include <memory>
#include <queue>
#include "Edge.h"
#include "LabeledDisjointSets.h"


class SolverUnweighted {
public:
    std::vector<std::vector<std::shared_ptr<Edge>>> adj_list;

    explicit SolverUnweighted(const std::vector<std::vector<int>> &adj_list_);

    void PrintMatching();

    void PrintAdjList();

    void PrintTreeData();

    void Solve();

private:
    const int n;  // the number of vertices
    std::vector<std::shared_ptr<Edge>> matched_edge;
    std::vector<std::shared_ptr<Edge>> minus_parents;
    std::vector<bool> plus;
    std::vector<bool> minus;
    std::queue<int> growable_vertices;
    LabeledDisjointSets cherry_blossoms; // labels are receptacles

    void GreedyInit();

    void Augment(std::vector<std::shared_ptr<Edge>> path, int start);

    // first_vertex and second_vertex are pluses
    // walks up the tree until the lca and makes everyone plusminus, 
    // adding the new pluses to the queue and updating parent pointers
    void MakeCherryBlossom(std::shared_ptr<Edge> edge_plus_plus, int root);

    int PlusPlusLCA(int first_vertex, int second_vertex) const;

    std::pair<int, int> PathUpperBounds(int first_vertex, int second_vertex);

    void UpdatePath(int lower_vertex, int upper_vertex);

    void ClearTrees();
};

#endif
