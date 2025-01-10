#ifndef BLOSSOM_VI_TESTER_H
#define BLOSSOM_VI_TESTER_H

#include <vector>


class Tester {
public:
    double runtime;

    explicit Tester(const std::vector<std::vector<int>> &adj_list_,
                                  int greedy_init_type_, bool delete_edges_in_cherries_, bool verbose_ = false);

    bool Validate();

private:
    std::vector<std::vector<int>> adj_list; // (destination, is_matched)
    std::vector<int> matched_to;
    std::vector<int> blossom_index;
    bool verbose;
    int greedy_init_type;
    bool delete_edges_in_cherries;

    int LCA(int first_vertex, int second_vertex, const std::vector<int> &parents) const;

    bool AugmentingPathsExist();

    void MakeBlossom(const std::vector<int> &blossom, std::vector<int> &parents);

    int TopBlossom(int vertex) const;
};

#endif
