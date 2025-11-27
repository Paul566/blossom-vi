#include "TesterWeighted.h"

#include <chrono>
#include <iomanip>
#include <unordered_set>

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
        std::cout << "------------------------------------------------------------\niter " << i << std::endl;
        if (i == 95) {
            std::cout << "a;ldsfj" << std::endl;
        }

        std::vector<std::tuple<int, int, int> > edge_list = RandomCliqueFixedSize(num_vertices, weight_min, weight_max);

        auto start = std::chrono::high_resolution_clock::now();
        SolverWeighted solver = SolverWeighted(edge_list,{.compute_dual_certificate = verify_output, .verbose = verbose});
        solver.FindMinPerfectMatching();
        auto stop = std::chrono::high_resolution_clock::now();
        double runtime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
            stop - start).count()) / 1'000'000;
        runtimes.push_back(runtime);

        if (verify_output) {
            Verify(edge_list, solver);
        }
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

void TesterWeighted::Verify(const std::vector<std::tuple<int, int, int> > &edge_list, const SolverWeighted &solver) {
    std::vector<std::pair<int, int> > matching = solver.Matching();
    std::vector<std::tuple<int, int, int> > dual_solution = solver.DualCertificate();
    std::vector<std::vector<std::pair<int, int> > > adj_list = AdjList(edge_list);

    if (!IsPerfectMatching(adj_list, matching)) {
        solver.PrintGraph();
        throw std::runtime_error("Not a perfect matching");
    }

    if (!IsCorrectDualSolution(adj_list, dual_solution)) {
        solver.PrintGraph();
        throw std::runtime_error("Not a correct dual solution");
    }

    if (4 * PrimalObjective(adj_list, matching) != QuadrupledDualObjective(adj_list, dual_solution)) {
        solver.PrintGraph();
        throw std::runtime_error("The dual objective is not equal the primal objective");
    }
}

bool TesterWeighted::IsPerfectMatching(const std::vector<std::vector<std::pair<int, int> > > &adj_list,
                                       const std::vector<std::pair<int, int> > &matching) {
    // check if the matched edges exist
    std::unordered_set<std::pair<int, int>, PairHash> edge_set;
    for (int i = 0; i < static_cast<int>(adj_list.size()); ++i) {
        for (auto [to, weight] : adj_list[i]) {
            edge_set.emplace(i, to);
        }
    }
    for (auto pair : matching) {
        if (!edge_set.contains(pair)) {
            return false;
        }
    }

    // check if the matching is correct
    std::vector<bool> matched_vtx(adj_list.size(), false);
    for (auto [from, to] : matching) {
        if (matched_vtx[from]) {
            return false;
        }
        matched_vtx[from] = true;
        if (matched_vtx[to]) {
            return false;
        }
        matched_vtx[to] = true;
    }

    // check if the matching is perfect
    if (2 * matching.size() != adj_list.size()) {
        return false;
    }

    return true;
}

bool TesterWeighted::IsCorrectDualSolution(const std::vector<std::vector<std::pair<int, int> > > &adj_list,
                                           const std::vector<std::tuple<int, int, int> > &dual_solution) {
    // check that the dual variables of non-trivial blossoms are non-negative
    for (int i = static_cast<int>(adj_list.size()); i < static_cast<int>(dual_solution.size()); ++i) {
        if (std::get<1>(dual_solution[i]) < 0) {
            return false;
        }
    }

    // check that the edge slacks are non-negative
    for (int i = 0; i < static_cast<int>(adj_list.size()); ++i) {
        for (auto [to, weight] : adj_list[i]) {
            if (SlackQuadrupled(i, to, weight, dual_solution) < 0) {
                return false;
            }
        }
    }

    return true;
}

int64_t TesterWeighted::PrimalObjective(const std::vector<std::vector<std::pair<int, int> > > &adj_list,
                                        const std::vector<std::pair<int, int> > &matching) {
    std::unordered_map<std::pair<int, int>, int, PairHash> edge_weights;
    for (int i = 0; i < static_cast<int>(adj_list.size()); ++i) {
        for (auto [to, weight] : adj_list[i]) {
            edge_weights[{i, to}] = weight;
        }
    }

    int64_t answer = 0;
    for (auto [from, to] : matching) {
        answer += edge_weights[{from, to}];
    }
    return answer;
}

int64_t TesterWeighted::QuadrupledDualObjective(const std::vector<std::vector<std::pair<int, int> > > &adj_list,
                                                const std::vector<std::tuple<int, int, int> > &dual_solution) {
    int64_t answer = 0;
    for (auto [vertex, var, parent] : dual_solution) {
        answer += var;
    }
    return answer;
}

std::vector<std::vector<std::pair<int, int> > > TesterWeighted::AdjList(
    const std::vector<std::tuple<int, int, int> > &edge_list) {
    int n = 0;
    for (auto [from, to, weight] : edge_list) {
        if (from > n) {
            n = from;
        }
        if (to > n) {
            n = to;
        }
    }
    ++n;

    std::vector<std::vector<std::pair<int, int> > > adj_list(n, std::vector<std::pair<int, int> >());
    for (auto [from, to, weight] : edge_list) {
        adj_list[from].emplace_back(to, weight);
        adj_list[to].emplace_back(from, weight);
    }
    return adj_list;
}

int64_t TesterWeighted::SlackQuadrupled(int vtx1,
                                        int vtx2,
                                        int weight,
                                        const std::vector<std::tuple<int, int, int> > &dual_solution) {
    int64_t answer = 4 * weight;

    std::unordered_set<int> visited;
    int first = vtx1;
    while (first != -1) {
        answer -= std::get<1>(dual_solution[first]);
        visited.insert(first);
        first = std::get<2>(dual_solution[first]);
    }

    int second = vtx2;
    while (second != -1) {
        if (visited.contains(second)) {
            answer += std::get<1>(dual_solution[second]);
        } else {
            answer -= std::get<1>(dual_solution[second]);
        }
        second = std::get<2>(dual_solution[second]);
    }

    return answer;
}
