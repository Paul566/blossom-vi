#include <algorithm>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <random>
#include <set>

#include "SolverUnweighted.h"
#include "Tester.h"

void PrintVector(const std::vector<int> &numbers) {
    for (const int number : numbers) {
        std::cout << number << " ";
    }
    std::cout << std::endl;
}

std::vector<std::vector<int>> ReadEdgeList(const std::string &path) {
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

std::vector<std::vector<int>> RandomGraph(const int num_vertices, const int num_edges, std::mt19937 &generator) {
    std::vector<std::vector<int>> adj_list(num_vertices, std::vector<int>());

    if (static_cast<long>(num_edges) > static_cast<long>(num_vertices) * static_cast<long>(num_vertices - 1) / 2) {
        throw std::runtime_error("In RandomGraph: too many edges");
    }

    std::set<std::pair<int, int>> edges;

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

void RunSavedTests(int init_type, bool delete_edges_in_cherries) {
    const std::string prefix = std::filesystem::current_path().string() + "/../tests/random-graphs/";

    for (int n = 3; n <= 15; ++n) {
        for (int m = n; m <= n * (n - 1) / 2; ++m) {
            std::string filename_graph = std::to_string(n) + "-" + std::to_string(m) + ".txt";

            std::cout << filename_graph << ":" << std::endl;

            auto adj_list = ReadEdgeList(prefix + filename_graph);

            // SolverUnweighted solver = SolverUnweighted(adj_list);
            // solver.Solve();
            // solver.PrintMatching();
            // solver.PrintAdjList();

            Tester tester = Tester(adj_list, init_type, delete_edges_in_cherries);
            std::cout << "runtime: " << tester.runtime << " seconds" << std::endl;
            if (!tester.Validate()) {
                throw std::runtime_error("test failed");
            }
            std::cout << "\n";
        }
    }

    std::cout << "All tests passed!\n";
}

void RunRandomTests(int max_vertices, int num_tests, std::mt19937 &generator, int init_type, bool delete_edges_in_cherries) {
    std::uniform_int_distribution<> dist_vertices(1, max_vertices);

    for (int i = 0; i < num_tests; ++i) {
        std::cout << "test " << i << "\n";

        int num_vertices = dist_vertices(generator);
        std::uniform_int_distribution<> dist_edges(0, num_vertices * (num_vertices - 1) / 2);
        int num_edges = dist_edges(generator);

        auto adj_list = RandomGraph(num_vertices, num_edges, generator);

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

void MeasureTime(const std::string& path, int init_type, bool delete_edges_in_cherries) {
    auto adj_list = ReadEdgeList(path);
    Tester tester = Tester(adj_list, init_type, delete_edges_in_cherries);
    std::cout << "runtime: " << tester.runtime << " seconds" << std::endl;
    std::cout << "runtime/init_time ratio: " << tester.runtime / tester.init_time << std::endl;
}

std::vector<std::string> AllFiles(const std::string& directory_path) {
    std::vector<std::string> files;

    try {
        for (const auto& entry : std::filesystem::directory_iterator(directory_path)) {
            if (std::filesystem::is_regular_file(entry.status())) {
                files.push_back(entry.path().string());
            }
        }
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Cannot access directory: " << e.what() << std::endl;
    }

    std::sort(files.begin(), files.end());

    return files;
}

void MeasureAll(const std::string& directory_path, int init_type, bool delete_edges_in_cherries) {
    auto files = AllFiles(directory_path);
    for (const auto& file : files) {
        std::cout << file << std::endl;
        MeasureTime(file, init_type, delete_edges_in_cherries);
        std::cout << std::endl;
    }
}

int main() {
    std::mt19937 gen(239);
    int init_type = 1;
    bool delete_edges_in_cherries = false;

    MeasureAll("../tests", init_type, delete_edges_in_cherries);
    // MeasureTime("../tests/spatial2d-100000-391934.txt", init_type, delete_edges_in_cherries);
    // MeasureTime("../tests/random-100000-300000.txt", init_type, delete_edges_in_cherries);

    // RunSavedTests(init_type, delete_edges_in_cherries);
    // RunRandomTests(50, 100000, gen, init_type, delete_edges_in_cherries);

    // std::string prefix = std::filesystem::current_path().string() + "/../tests/random-graphs/";
    // int n = 5;
    // int m = 8;
    // std::string filename_graph = std::to_string(n) + "-" + std::to_string(m) + ".txt";
    // auto adj_list = ReadEdgeList(prefix + filename_graph);
    // SolverUnweighted solver = SolverUnweighted(adj_list, init_type, delete_edges_in_cherries, true);
    // solver.PrintAdjList();
    // solver.Solve();
    // solver.PrintMatching();

    return 0;
}
