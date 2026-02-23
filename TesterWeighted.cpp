#include "TesterWeighted.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_set>

CliqueGenerator::CliqueGenerator(int num_vertices_, int weight_min_, int weight_max_) : num_vertices(num_vertices_),
    weight_min(weight_min_), weight_max(weight_max_) {
    if (num_vertices % 2 != 0) {
        throw std::invalid_argument("In CliqueGenerator: num_vertices must be even");
    }
}

EdgeListType CliqueGenerator::Generate(std::mt19937 *generator) const {
    std::vector<std::tuple<int, int, int> > edge_list;
    edge_list.reserve(num_vertices * (num_vertices - 1) / 2);

    std::uniform_int_distribution<> dist_weight(weight_min, weight_max);

    for (int vertex = 0; vertex < num_vertices; ++vertex) {
        for (int to = vertex + 1; to < num_vertices; ++to) {
            int weight = dist_weight(*generator);
            edge_list.emplace_back(vertex, to, weight);
        }
    }

    return edge_list;
}

MatchingPlusGraphGenerator::MatchingPlusGraphGenerator(int num_vertices_,
                                                       int num_edges_,
                                                       int weight_min_,
                                                       int weight_max_) : num_vertices(num_vertices_),
                                                                          num_edges(num_edges_),
                                                                          weight_min(weight_min_),
                                                                          weight_max(weight_max_) {
    if (num_vertices % 2 != 0) {
        throw std::invalid_argument("In MatchingPlusGraphGenerator: num_vertices must be even");
    }

    int min_edges_needed = num_vertices / 2; // size of the perfect matching
    if (num_edges < min_edges_needed) {
        throw std::invalid_argument("In MatchingPlusGraphGenerator: num_edges must be >= num_vertices/2");
    }
    if (num_edges > num_vertices * (num_vertices - 1) / 2) {
        throw std::invalid_argument("In MatchingPlusGraphGenerator: too many edges");
    }
}

EdgeListType MatchingPlusGraphGenerator::Generate(std::mt19937 *generator) const {
    std::vector<std::tuple<int, int, int> > edges;
    edges.reserve(num_edges);

    std::uniform_int_distribution<> weight_dist(weight_min, weight_max);

    // random perfect matching
    std::vector<int> vertices(num_vertices);
    for (int i = 0; i < num_vertices; i++) {
        vertices[i] = i;
    }

    std::shuffle(vertices.begin(), vertices.end(), *generator);

    std::unordered_set<std::pair<int, int>, PairHash> edge_set;
    edge_set.reserve(num_edges * 2);

    for (int i = 0; i < num_vertices; i += 2) {
        int u = vertices[i];
        int v = vertices[i + 1];
        if (u > v) {
            std::swap(u, v);
        }

        int w = weight_dist(*generator);
        edges.emplace_back(u, v, w);
        edge_set.insert({u, v});
    }

    // add random edges until total reaches num_edges
    std::uniform_int_distribution<> vertex_dist(0, num_vertices - 1);

    while (static_cast<int>(edges.size()) < num_edges) {
        int u = vertex_dist(*generator);
        int v = vertex_dist(*generator);
        if (u == v) continue;
        if (u > v) {
            std::swap(u, v);
        }

        std::pair<int, int> e = {u, v};
        if (edge_set.contains(e)) {
            continue;
        }

        int w = weight_dist(*generator);
        edges.emplace_back(u, v, w);
        edge_set.insert(e);
    }

    return edges;
}

TesterWeighted::TesterWeighted(bool verify_output_, int seed) : generator(seed), verify_output(verify_output_) {
}

void TesterWeighted::RunInstances(const GraphGenerator &graph_generator,
                                  int num_iter,
                                  bool verbose) {
    std::cout << std::setprecision(3);

    std::vector<double> runtimes;
    runtimes.reserve(num_iter);
    std::vector<double> init_times;
    init_times.reserve(num_iter);

    for (int i = 0; i < num_iter; ++i) {
        std::cout << "------------------------------------------------------------\niter " << i << std::endl;

        std::vector<std::tuple<int, int, int> > edge_list = graph_generator.Generate(&generator);

        auto start = std::chrono::high_resolution_clock::now();
        VzhuhSolver solver = VzhuhSolver(edge_list,
                                         {
                                             .compute_dual_certificate = verify_output,
                                             .verbose = (verbose || (i == 10318)), .print_statistics = false,
                                             .debug = verify_output
                                         });
        solver.FindMinPerfectMatching();
        auto stop = std::chrono::high_resolution_clock::now();
        double runtime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
            stop - start).count()) / 1'000'000;
        runtimes.push_back(runtime);
        if (verbose) {
            std::cout << "runtime: " << runtime << std::endl;
        }

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

void TesterWeighted::MeasureBenchmark(const std::string &path, int num_iter, double max_time_per_instance) {
    auto files = AllFiles(path);
    std::cout << std::setprecision(3);

    for (const auto &file : files) {
        MeasureInstance(file, num_iter, max_time_per_instance);
    }
}

void TesterWeighted::MeasureInstance(const std::string &filename,
                                     int num_iter,
                                     double max_time_per_instance,
                                     bool with_debug) {
    std::cout << filename << std::endl;
    EdgeListType edge_list = ReadWeightedEdgeList(filename);

    std::vector<double> runtimes;
    std::vector<double> init_times;

    std::cout << "runtimes:\t\t\t\t\t" << std::flush;
    double total_time = 0.;
    int real_iters = 0;
    for (int i = 0; i < num_iter; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        VzhuhSolver solver = VzhuhSolver(edge_list,
                                         {
                                             .compute_dual_certificate = false, .verbose = false,
                                             .print_statistics = false, .debug = with_debug
                                         });
        solver.FindMinPerfectMatching();
        auto stop = std::chrono::high_resolution_clock::now();
        double runtime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
            stop - start).count()) / 1'000'000;
        runtimes.push_back(runtime);
        total_time += runtime;
        ++real_iters;

        std::cout << runtime << "\t" << std::flush;
        if (total_time > max_time_per_instance) {
            break;
        }
    }

    std::cout << "\naverage runtime: " << std::accumulate(runtimes.begin(), runtimes.end(), 0.) / real_iters <<
        std::endl;
    std::cout << std::endl;
}

void TesterWeighted::Verify(const std::vector<std::tuple<int, int, int> > &edge_list, const VzhuhSolver &solver) {
    const std::vector<std::pair<int, int> > &matching = solver.Matching();
    const std::vector<std::tuple<int, int, int> > &dual_solution = solver.DualCertificate();
    std::vector<std::vector<std::pair<int, int> > > adj_list = AdjList(edge_list);

    if (!IsPerfectMatching(adj_list, matching)) {
        throw std::runtime_error("Not a perfect matching");
    }

    if (!IsCorrectDualSolution(adj_list, dual_solution)) {
        throw std::runtime_error("Not a correct dual solution");
    }

    if (4 * PrimalObjective(adj_list, matching) != QuadrupledDualObjective(adj_list, dual_solution)) {
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

std::vector<std::string> TesterWeighted::AllFiles(const std::string &directory_path) {
    std::vector<std::string> files;

    try {
        for (const auto &entry : std::filesystem::directory_iterator(directory_path)) {
            if (std::filesystem::is_regular_file(entry.status())) {
                files.push_back(entry.path().string());
            }
        }
    } catch (const std::filesystem::filesystem_error &e) {
        std::cerr << "Cannot access directory: " << e.what() << std::endl;
    }

    std::sort(files.begin(), files.end());

    return files;
}

EdgeListType TesterWeighted::ReadWeightedEdgeList(const std::string &filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    std::vector<std::tuple<int, int, int> > edge_list;
    int u, v, weight;

    while (infile >> u >> v >> weight) {
        edge_list.emplace_back(u, v, weight);
    }

    return edge_list;
}
