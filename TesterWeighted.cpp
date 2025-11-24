#include "TesterWeighted.h"

#include <chrono>
#include <iomanip>

TesterWeighted::TesterWeighted(bool verify_output_, int seed) : generator(seed), verify_output(verify_output_) {
}

void TesterWeighted::RunRandomCliques(int num_vertices, int weight_min, int weight_max, int num_iter, bool verbose) {
    std::cout << "Testing random weighted cliques" << std::endl;

    std::cout << std::setprecision(3);

    std::vector<double> runtimes;
    runtimes.reserve(num_iter);
    std::vector<double> init_times;
    init_times.reserve(num_iter);

    for (int i = 0; i < num_iter; ++i) {
        std::vector<std::tuple<int, int, int> > edge_list = RandomCliqueFixedSize(num_vertices, weight_min, weight_max);

        auto start = std::chrono::high_resolution_clock::now();
        SolverWeighted solver = SolverWeighted(edge_list, {.verbose = verbose});
        solver.FindMinPerfectMatching();
        auto stop = std::chrono::high_resolution_clock::now();

        double runtime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
            stop - start).count()) / 1'000'000;

        std::cout << "iter " << i << std::endl;

        if (solver.dual_objective != solver.primal_objective) {
            solver.PrintElementaryAdjList();
            throw std::runtime_error("");
        }

        runtimes.push_back(runtime);
    }

    std::cout << "All tests passed!" << std::endl;
    std::cout << "runtimes:\t\t\t\t\t";
    for (int i = 0; i < num_iter; ++i) {
        std::cout << runtimes[i] << "\t";
    }
    std::cout << "\naverage runtime: " << std::accumulate(runtimes.begin(), runtimes.end(), 0.) / num_iter << std::endl;
}

std::vector<std::tuple<int, int, int> > TesterWeighted::RandomClique(int max_num_vertices,
                                                                     int weight_min,
                                                                     int weight_max) {
    // returns the edge list
    std::uniform_int_distribution<> dist_n(1, max_num_vertices / 2);
    int num_vertices = dist_n(generator);
    num_vertices *= 2;

    return RandomCliqueFixedSize(num_vertices, weight_min, weight_max);
}

std::vector<std::tuple<int, int, int> > TesterWeighted::RandomCliqueFixedSize(int num_vertices,
                                                                              int weight_min,
                                                                              int weight_max) {
    std::vector<std::tuple<int, int, int> > edge_list;
    edge_list.reserve(num_vertices * (num_vertices - 1) / 2);

    std::uniform_int_distribution<> dist_weight(weight_min, weight_max);

    for (int vertex = 0; vertex < num_vertices; ++vertex) {
        for (int to = vertex + 1; to < num_vertices; ++to) {
            int weight = dist_weight(generator);
            edge_list.emplace_back(vertex, to, weight);
        }
    }

    return edge_list;
}
