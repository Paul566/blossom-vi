#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "TesterWeighted.h"

namespace {

std::filesystem::path ResolveBenchmarkPath(const std::string& user_path) {
    const std::filesystem::path direct_path(user_path);
    if (std::filesystem::exists(direct_path)) {
        return direct_path;
    }

    const std::filesystem::path parent_path = std::filesystem::path("..") / user_path;
    if (std::filesystem::exists(parent_path)) {
        return parent_path;
    }

    throw std::runtime_error("Benchmark path does not exist: " + user_path);
}

int ParseIntArg(const char* value, const char* arg_name) {
    try {
        return std::stoi(value);
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string("Invalid integer for ") + arg_name + ": " + value);
    }
}

double ParseDoubleArg(const char* value, const char* arg_name) {
    try {
        return std::stod(value);
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string("Invalid number for ") + arg_name + ": " + value);
    }
}

} // namespace

int main(int argc, char** argv) {
    TesterWeighted tester(true, 239);
    // tester.RunInstances(MatchingPlusGraphGenerator(6, 10, 0, 10), 100000, false);
    // tester.RunInstances(MatchingPlusGraphGenerator(100, 500, -1000, 1000), 1000, false);

    const std::string benchmark_path_arg = argc > 1 ? argv[1] : "tests-weighted";
    const int iterations = argc > 2 ? ParseIntArg(argv[2], "iterations") : 1;
    const double max_time_per_instance = argc > 3 ? ParseDoubleArg(argv[3], "max_time_per_instance") : 20.;

    tester.MeasureBenchmark(ResolveBenchmarkPath(benchmark_path_arg).string(), iterations, max_time_per_instance);

    // tester.MeasureInstance("../tests-weighted/delaunay-100000-299968", 1, 20, false);
    // tester.MeasureInstance("../tests-weighted/delaunay-1000000-2999962", 1, 20, false);
    // tester.MeasureInstance("../tests-weighted/delaunay-1000000-2999965", 1, 20, false, true);
    // tester.MeasureInstance("../tests-weighted/dan59296-177299", 10, 20, false);
    // tester.MeasureInstance("../tests-weighted/sra104815-314222", 10, 20, false);
    // tester.MeasureInstance("../tests-weighted/ara238025-713594", 1, 20, false);
    // tester.MeasureInstance("../tests-weighted/lrb744710-2233725", 1, 20, false);
    // tester.MeasureInstance("../tests-weighted/lra498378-1494967", 1, 20, false);
    // tester.MeasureInstance("../tests-weighted/random-10000-100000", 1, 20, false);
    // tester.MeasureInstance("../tests-weighted/random-100000-1000000", 1, 20, false);

    return 0;
}
