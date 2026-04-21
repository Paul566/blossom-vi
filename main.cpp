#include <cstdint>
#include <cstdlib>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "MWPMSolver.h"

namespace {

struct GraphInput {
    int num_vertices = 0;
    int num_edges = 0;
    std::vector<std::tuple<int, int, int>> edges;
};

struct CommandLineOptions {
    std::filesystem::path graph_path;
    std::filesystem::path matching_output_path;
    bool print_matching = false;
    bool print_runtime = false;
    bool save_matching = false;
};

void PrintUsage(const char* program_name) {
    std::cerr << "Usage: " << program_name
              << " [--print-matching] [--print-runtime] <graph-file> [matching-output-file]\n"
              << "\n"
              << "Input format:\n"
              << "  n m\n"
              << "  endpoint1 endpoint2 weight\n"
              << "  ...\n"
              << "\n"
              << "By default stdout contains only the optimal weight.\n"
              << "Use --print-matching to also print the matched edges to stdout.\n"
              << "Use --print-runtime to print the time spent in FindMinPerfectMatching().\n"
              << "Use matching-output-file to save the matched edges to a file.\n";
}

CommandLineOptions ParseCommandLine(int argc, char** argv) {
    CommandLineOptions options;
    std::vector<std::string> positional_args;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--print-matching") {
            options.print_matching = true;
        } else if (arg == "--print-runtime") {
            options.print_runtime = true;
        } else if (arg == "--help" || arg == "-h") {
            PrintUsage(argv[0]);
            std::exit(0);
        } else if (!arg.empty() && arg[0] == '-') {
            throw std::invalid_argument("Unknown option: " + arg);
        } else {
            positional_args.push_back(arg);
        }
    }

    if (positional_args.empty() || positional_args.size() > 2) {
        throw std::invalid_argument("Expected <graph-file> and optional [matching-output-file]");
    }

    options.graph_path = positional_args[0];
    if (positional_args.size() == 2) {
        options.matching_output_path = positional_args[1];
        options.save_matching = true;
    }
    return options;
}

GraphInput ReadGraph(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open input file: " + path.string());
    }

    GraphInput graph;
    if (!(input >> graph.num_vertices >> graph.num_edges)) {
        throw std::runtime_error("Failed to read graph header: " + path.string());
    }
    if (graph.num_vertices < 0 || graph.num_edges < 0) {
        throw std::runtime_error("Graph header must contain non-negative integers");
    }
    if (graph.num_vertices % 2 != 0) {
        throw std::runtime_error("The number of vertices must be even");
    }

    graph.edges.reserve(graph.num_edges);
    std::vector<int> degrees(graph.num_vertices, 0);
    for (int i = 0; i < graph.num_edges; ++i) {
        int from = 0;
        int to = 0;
        int weight = 0;
        if (!(input >> from >> to >> weight)) {
            throw std::runtime_error("Failed to read edge " + std::to_string(i));
        }
        if (from < 0 || from >= graph.num_vertices || to < 0 || to >= graph.num_vertices) {
            throw std::runtime_error("Edge endpoint out of range at edge " + std::to_string(i));
        }
        if (from == to) {
            throw std::runtime_error("Self-loops are not supported");
        }

        graph.edges.emplace_back(from, to, weight);
        ++degrees[from];
        ++degrees[to];
    }

    int extra = 0;
    if (input >> extra) {
        throw std::runtime_error("Input contains more data than declared by m");
    }

    for (int vertex = 0; vertex < graph.num_vertices; ++vertex) {
        if (degrees[vertex] == 0) {
            throw std::runtime_error("Vertex " + std::to_string(vertex) + " is isolated");
        }
    }

    return graph;
}

void PrintMatching(const MWPMSolver& solver) {
    for (const auto& [from, to] : solver.Matching()) {
        std::cout << "match " << from << ' ' << to << '\n';
    }
}

void SaveMatching(const MWPMSolver& solver, const std::filesystem::path& output_path) {
    std::ofstream output(output_path);
    if (!output.is_open()) {
        throw std::runtime_error("Failed to open matching output file: " + output_path.string());
    }

    for (const auto& [from, to] : solver.Matching()) {
        output << from << ' ' << to << '\n';
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const CommandLineOptions options = ParseCommandLine(argc, argv);
        const GraphInput graph = ReadGraph(options.graph_path);
        MWPMSolver solver(graph.edges, {
            .compute_dual_certificate = false,
            .verbose = false,
            .print_statistics = false,
            .debug = false,
        });
        const auto start_time = std::chrono::steady_clock::now();
        solver.FindMinPerfectMatching();
        const auto stop_time = std::chrono::steady_clock::now();
        const std::chrono::duration<double> runtime = stop_time - start_time;

        std::cout << "Optimal weight: " << solver.primal_objective << '\n';
        if (options.print_runtime) {
            std::cout << "Runtime: " << runtime.count() << " seconds" << '\n';
        }
        if (options.print_matching) {
            PrintMatching(solver);
        }
        if (options.save_matching) {
            SaveMatching(solver, options.matching_output_path);
        }
    } catch (const std::exception& error) {
        PrintUsage(argv[0]);
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
