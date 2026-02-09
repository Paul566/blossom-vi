#ifndef BLOSSOM_VI_TESTERWEIGHTED_H
#define BLOSSOM_VI_TESTERWEIGHTED_H

#include <functional>
#include <random>
#include <vector>

// #include "SolverWeighted.h"
// #include "Solver.h"
#include "VzhuhSolver.h"

using EdgeListType = std::vector<std::tuple<int, int, int> >;
using AdjListType = std::vector<std::vector<std::pair<int, int> > >;

struct PairHash {
    std::size_t operator()(const std::pair<int, int> &p) const noexcept {
        uint64_t x = (static_cast<uint64_t>(static_cast<uint32_t>(p.first)) << 32) | static_cast<uint32_t>(p.second);
        return std::hash<uint64_t>{}(x);
    }
};

class GraphGenerator {
    public:
        virtual ~GraphGenerator() = default;
        virtual EdgeListType Generate(std::mt19937 * generator) const = 0;
};
class CliqueGenerator : public GraphGenerator {
    // generates graphs that are unions of a random perfect matching and a random graph
    public:
        CliqueGenerator(int num_vertices_, int weight_min_, int weight_max_);
        EdgeListType Generate(std::mt19937 * generator) const override;
        int num_vertices;
        int weight_min;
        int weight_max;
};
class MatchingPlusGraphGenerator : public GraphGenerator {
    // generates graphs that are unions of a random perfect matching and a random graph
    public:
        MatchingPlusGraphGenerator(int num_vertices_, int num_edges_, int weight_min_, int weight_max_);
        EdgeListType Generate(std::mt19937 * generator) const override;
        int num_vertices;
        int num_edges;
        int weight_min;
        int weight_max;
};

class TesterWeighted {
    public:
        bool verify_output;

        explicit TesterWeighted(bool verify_output_, int seed = 239);

        void RunInstances(const GraphGenerator &graph_generator,
                          int num_iter = 1,
                          bool verbose = true);

        static void MeasureBenchmark(const std::string &path, int num_iter, double max_time_per_instance = 60.);
        static void MeasureInstance(const std::string &filename, int num_iter, double max_time_per_instance = 60.,
                                     bool with_debug = false);

    private:
        std::mt19937 generator;

        void Verify(const EdgeListType &edge_list, const VzhuhSolver &solver);
        static bool IsPerfectMatching(const AdjListType &adj_list,
                                      const std::vector<std::pair<int, int> > &matching);
        bool IsCorrectDualSolution(const AdjListType &adj_list,
                                   const std::vector<std::tuple<int, int, int> > &dual_solution);
        static int64_t PrimalObjective(const AdjListType &adj_list,
                                       const std::vector<std::pair<int, int> > &matching);
        static int64_t QuadrupledDualObjective(const AdjListType &adj_list,
                                               const std::vector<std::tuple<int, int, int> > &dual_solution);
        // dual solution comes in the form (index, quadrupled dual variable, index of the blossom parent or -1)

        static AdjListType AdjList(
            const EdgeListType &edge_list);
        int64_t SlackQuadrupled(int vtx1,
                                int vtx2,
                                int weight,
                                const std::vector<std::tuple<int, int, int> > &dual_solution);

        static std::vector<std::string> AllFiles(const std::string &directory_path);
        static EdgeListType ReadWeightedEdgeList(const std::string &filename);
};

#endif //BLOSSOM_VI_TESTERWEIGHTED_H
