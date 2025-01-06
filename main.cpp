#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>
#include "SolverUnweighted.h"
#include "Tester.h"


void PrintVector(const std::vector<int>& numbers) {
    for (const int number : numbers) {
        std::cout << number << " ";
    }
    std::cout << std::endl;
}

std::vector<std::vector<int>> ReadEdgeList(const std::string &path) {
    std::vector<std::vector<int>> adj_list;
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

void RunRandomTests() {
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

            Tester tester = Tester(adj_list);
            std::cout << "runtime: " << tester.runtime << " seconds" << std::endl;
            if (!tester.Validate()) {
                throw std::runtime_error("test failed");
            }
            std::cout << "\n";
        }
    }

    std::cout << "All tests passed!\n";
}

int main() {
    RunRandomTests();

    // std::string prefix = std::filesystem::current_path().string() + "/../tests/random-graphs/";
    // int n = 5;
    // int m = 8;
    // std::string filename_graph = std::to_string(n) + "-" + std::to_string(m) + ".txt";
    // auto adj_list = ReadEdgeList(prefix + filename_graph);
    // Tester tester = Tester(adj_list);
    // tester.Validate();


    // std::vector<std::vector<int>> adj_list;
    // adj_list.push_back(std::vector<int>({1, 2, 3}));
    // adj_list.push_back(std::vector<int>({0, 2, 3}));
    // adj_list.push_back(std::vector<int>({0, 1}));
    // adj_list.push_back(std::vector<int>({0, 1}));

    // SolverUnweighted solver = SolverUnweighted(adj_list);
    // solver.Solve();
    // solver.Print();

    return 0;
}