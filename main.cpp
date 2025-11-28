#include <algorithm>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <random>
#include <set>

#include "unweighted/SolverUnweighted.h"
#include "unweighted/Tester.h"
#include "SolverWeighted.h"
#include "TesterWeighted.h"

void PrintVector(const std::vector<int> &numbers) {
    for (const int number : numbers) {
        std::cout << number << " ";
    }
    std::cout << std::endl;
}

std::vector<std::vector<int> > ReadUnweightedEdgeList(const std::string &path) {
    std::vector<std::vector<int> > adj_list;
    std::fstream input_file(path);

    if (input_file.is_open()) {
        std::string line;
        std::getline(input_file, line);
        std::istringstream stream(line);

        int num_vertices, num_edges;
        stream >> num_vertices >> num_edges;

        adj_list.reserve(num_vertices);
        for (int i = 0; i < num_vertices; ++i) {
            adj_list.emplace_back();
        }

        for (int i = 0; i < num_edges; ++i) {
            std::getline(input_file, line);
            stream = std::istringstream(line);

            int head, tail;
            stream >> head >> tail;
            adj_list[head].push_back(tail);
            adj_list[tail].push_back(head);
        }
    } else {
        throw std::runtime_error("input file " + path + " not found");
    }
    input_file.close();

    return adj_list;
}

std::vector<std::vector<int> > RandomUnweightedGraph(const int num_vertices,
                                                     const int num_edges,
                                                     std::mt19937 &generator) {
    std::vector<std::vector<int> > adj_list(num_vertices, std::vector<int>());

    if (static_cast<long>(num_edges) > static_cast<long>(num_vertices) * static_cast<long>(num_vertices - 1) / 2) {
        throw std::runtime_error("In RandomGraph: too many edges");
    }

    std::set<std::pair<int, int> > edges;

    std::uniform_int_distribution<> dist(0, num_vertices - 1);
    while (edges.size() < static_cast<size_t>(num_edges)) {
        int first_vertex = dist(generator); // Random vertex u
        int second_vertex = dist(generator);

        if (first_vertex != second_vertex) {
            auto edge = std::minmax(first_vertex, second_vertex);
            if (!edges.contains(edge)) {
                edges.insert(edge);
                adj_list[edge.first].push_back(edge.second);
                adj_list[edge.second].push_back(edge.first);
            }
        }
    }

    return adj_list;
}

std::vector<std::tuple<int, int, int> > RandomWeightedGraph(const int num_vertices,
                                                            const int num_edges,
                                                            std::mt19937 &generator,
                                                            const int weight_min,
                                                            const int weight_max) {
    // returns the edge list
    std::vector<std::tuple<int, int, int> > edge_list;

    std::uniform_int_distribution<> dist_weight(weight_min, weight_max);

    auto adj_list = RandomUnweightedGraph(num_vertices, num_edges, generator);
    for (int vertex = 0; vertex < num_vertices; ++vertex) {
        for (int to : adj_list[vertex]) {
            if (to > vertex) {
                int weight = dist_weight(generator);
                edge_list.emplace_back(vertex, to, weight);
            }
        }
    }

    return edge_list;
}

std::vector<std::tuple<int, int, int> > RandomWeightedClique(const int num_vertices,
                                                            std::mt19937 &generator,
                                                            const int weight_min,
                                                            const int weight_max) {
    // returns the edge list
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

void RunSavedUnweightedTests(int init_type, bool delete_edges_in_cherries) {
    const std::string prefix = std::filesystem::current_path().string() + "/../tests/random-graphs/";

    for (int n = 3; n <= 15; ++n) {
        for (int m = n; m <= n * (n - 1) / 2; ++m) {
            std::string filename_graph = std::to_string(n) + "-" + std::to_string(m) + ".txt";

            std::cout << filename_graph << ":" << std::endl;

            auto adj_list = ReadUnweightedEdgeList(prefix + filename_graph);

            bool verbose = false;
            if ((n == 9) && (m == 14)) {
                verbose = true;
            }

            Tester tester = Tester(adj_list, init_type, delete_edges_in_cherries, verbose);
            std::cout << "runtime: " << tester.runtime << " seconds" << std::endl;
            if (!tester.Validate()) {
                throw std::runtime_error("test failed");
            }
            std::cout << "\n";
        }
    }

    std::cout << "All tests passed!\n";
}

void RunRandomUnweightedTests(int max_vertices,
                              int num_tests,
                              std::mt19937 &generator,
                              int init_type,
                              bool delete_edges_in_cherries) {
    std::uniform_int_distribution<> dist_vertices(1, max_vertices);

    for (int i = 0; i < num_tests; ++i) {
        std::cout << "test " << i << "\n";

        int num_vertices = dist_vertices(generator);
        std::uniform_int_distribution<> dist_edges(0, num_vertices * (num_vertices - 1) / 2);
        int num_edges = dist_edges(generator);

        auto adj_list = RandomUnweightedGraph(num_vertices, num_edges, generator);

        bool verbose = false;
        if (i == 14080) {
            verbose = true;
        }
        Tester tester = Tester(adj_list, init_type, delete_edges_in_cherries, verbose);

        if (!tester.Validate()) {
            throw std::runtime_error("test failed");
        }
    }

    std::cout << "All tests passed!\n";
}

void MeasureUnweightedTime(const std::string &path, int init_type, bool delete_edges_in_cherries, int num_iter = 1) {
    std::cout << std::setprecision(3);

    std::vector<double> runtimes;
    runtimes.reserve(num_iter);
    std::vector<double> init_times;
    init_times.reserve(num_iter);

    for (int i = 0; i < num_iter; ++i) {
        auto adj_list = ReadUnweightedEdgeList(path);
        Tester tester = Tester(adj_list, init_type, delete_edges_in_cherries);
        runtimes.push_back(tester.runtime);
        init_times.push_back(tester.init_time);
    }

    std::cout << "runtimes:\t\t\t\t\t";
    for (int i = 0; i < num_iter; ++i) {
        std::cout << runtimes[i] << "\t";
    }
    std::cout << "\nruntime/init_time ratios:\t";
    for (int i = 0; i < num_iter; ++i) {
        std::cout << (runtimes[i] / init_times[i]) << "\t";
    }
    std::cout << "\naverage runtime: " << std::accumulate(runtimes.begin(), runtimes.end(), 0.) / num_iter << std::endl;
}

std::vector<std::string> AllFiles(const std::string &directory_path) {
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

void MeasureAllUnweighted(const std::string &directory_path, int init_type, bool delete_edges_in_cherries) {
    auto files = AllFiles(directory_path);
    for (const auto &file : files) {
        std::cout << file << std::endl;
        MeasureUnweightedTime(file, init_type, delete_edges_in_cherries);
        std::cout << std::endl;
    }
}

int main() {
    // std::mt19937 gen(239);
    TesterWeighted tester(true, 239);
    // tester.RunRandomCliques(8, 0, 5, 10000, true);
    tester.RunRandomCliques(100, -1000, 1000, 100, false);

    return 0;
}
