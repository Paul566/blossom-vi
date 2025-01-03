#ifndef BLOSSOM_VI_TESTER_H
#define BLOSSOM_VI_TESTER_H


#include <vector>
#include <iostream>
#include <chrono>
#include <unordered_map>
#include <queue>
#include "SolverUnweighted.h"


class Tester {
public:
    double runtime;

    Tester(const std::vector<std::vector<int>> &adj_list_);

    bool Validate();

private:
    std::vector<std::vector<std::tuple<int, bool, int>>> adj_list_out; // (destination, is_matched, edge_index)
    std::vector<std::vector<std::tuple<int, bool, int>>> adj_list_in;

    bool IsAMatching();

    bool AugmentingPathsExist();

    int MatchingSize();

    std::vector<int> UnmatchedVertices();

    // returns the adjacency list of the auxiliary graph and the number of terminal vertices
    std::pair<std::vector<std::vector<int>>, int> AuxGraph();

    int NumberOfEdges();
};

#endif
