#ifndef BLOSSOM_VI_TESTERWEIGHTED_H
#define BLOSSOM_VI_TESTERWEIGHTED_H

#include <random>
#include <vector>

#include "SolverWeighted.h"

class TesterWeighted {
    public:
        TesterWeighted(bool verify_output_, int seed = 239);

        void RunRandomCliques(int num_vertices, int weight_min, int weight_max, int num_iter = 1, bool verbose = true);

    private:
        bool verify_output;
        std::mt19937 generator;

        std::vector<std::tuple<int, int, int> > RandomClique(int max_num_vertices,
                                                                     int weight_min,
                                                                     int weight_max);

        std::vector<std::tuple<int, int, int> > RandomCliqueFixedSize(int num_vertices,
                                                                     int weight_min,
                                                                     int weight_max);
};

#endif //BLOSSOM_VI_TESTERWEIGHTED_H
