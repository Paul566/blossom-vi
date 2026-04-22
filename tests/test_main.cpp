#include <iostream>
#include <stdexcept>
#include <string>

#include "TesterWeighted.h"

namespace {

int ParseIntArg(const char* value, const char* arg_name) {
    try {
        return std::stoi(value);
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string("Invalid integer for ") + arg_name + ": " + value);
    }
}

void PrintUsage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " [num_vertices] [num_edges] [iterations] [seed]\n"
              << "\n"
              << "Runs randomized correctness tests with primal/dual verification.\n";
}

}  // namespace

int main(int argc, char** argv) {
    // TODO add testing with no fractional matching initialization
    // TODO add another family of tests
    // TODO add tests for throwing if there is no perfect matching
    try {
        if (argc > 5 || (argc > 1 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h"))) {
            PrintUsage(argv[0]);
            return argc > 5 ? 1 : 0;
        }

        const int num_vertices = argc > 1 ? ParseIntArg(argv[1], "num_vertices") : 100;
        const int num_edges = argc > 2 ? ParseIntArg(argv[2], "num_edges") : 500;
        const int iterations = argc > 3 ? ParseIntArg(argv[3], "iterations") : 100;
        const int seed = argc > 4 ? ParseIntArg(argv[4], "seed") : 239;

        TesterWeighted tester(true, seed);
        tester.RunInstances(MatchingPlusGraphGenerator(num_vertices, num_edges, -1000, 1000), iterations, false);
    } catch (const std::exception& error) {
        PrintUsage(argv[0]);
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
