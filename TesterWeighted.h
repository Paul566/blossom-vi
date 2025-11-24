#ifndef BLOSSOM_VI_TESTERWEIGHTED_H
#define BLOSSOM_VI_TESTERWEIGHTED_H

#include <random>
#include <vector>

#include "SolverWeighted.h"

struct PairHash {
    std::size_t operator()(const std::pair<int, int> &p) const noexcept {
        uint64_t x = (static_cast<uint64_t>(static_cast<uint32_t>(p.first)) << 32) | static_cast<uint32_t>(p.second);
        return std::hash<uint64_t>{}(x);
    }
};

class TesterWeighted {
    public:
        bool verify_output;

        explicit TesterWeighted(bool verify_output_, int seed = 239);

        void RunRandomCliques(int num_vertices, int weight_min, int weight_max, int num_iter = 1, bool verbose = true);

    private:
        std::mt19937 generator;

        std::vector<std::tuple<int, int, int> > RandomClique(int max_num_vertices,
                                                             int weight_min,
                                                             int weight_max);

        std::vector<std::tuple<int, int, int> > RandomCliqueFixedSize(int num_vertices,
                                                                      int weight_min,
                                                                      int weight_max);

        void Verify(const std::vector<std::tuple<int, int, int> > &edge_list, const SolverWeighted &solver);

        static bool IsPerfectMatching(const std::vector<std::vector<std::pair<int, int> > > &adj_list,
                                      const std::vector<std::pair<int, int> > &matching);

        bool IsCorrectDualSolution(const std::vector<std::vector<std::pair<int, int> > > &adj_list,
                                   const std::vector<std::tuple<int, int, int> > &dual_solution);

        static int64_t PrimalObjective(const std::vector<std::vector<std::pair<int, int> > > &adj_list,
                                       const std::vector<std::pair<int, int> > &matching);

        static int64_t QuadrupledDualObjective(const std::vector<std::vector<std::pair<int, int> > > &adj_list,
                                               const std::vector<std::tuple<int, int, int> > &dual_solution);
        // dual solution comes in the form (index, quadrupled dual variable, index of the blossom parent or -1)

        static std::vector<std::vector<std::pair<int, int> > > AdjList(
            const std::vector<std::tuple<int, int, int> > &edge_list);

        int64_t SlackQuadrupled(int vtx1,
                                int vtx2,
                                int weight,
                                const std::vector<std::tuple<int, int, int> > &dual_solution);
};

#endif //BLOSSOM_VI_TESTERWEIGHTED_H
