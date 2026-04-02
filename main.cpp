#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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

struct GraphSize {
    int n = 0;
    int m = 0;
};

GraphSize ReadGraphSize(const std::filesystem::path& filename) {
    std::ifstream input(filename);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open benchmark file: " + filename.string());
    }

    GraphSize size;
    if (!(input >> size.n >> size.m)) {
        throw std::runtime_error("Failed to read benchmark header: " + filename.string());
    }

    return size;
}

void ExportRuntimes(const std::string& benchmark_path_arg, int iterations, double max_time_per_instance) {
    const std::filesystem::path benchmark_dir = ResolveBenchmarkPath(benchmark_path_arg);
    const std::filesystem::path output_dir =
        std::filesystem::path("runtimes") / "blossom-vi" / benchmark_dir.filename();
    const std::filesystem::path output_file = output_dir / "runtimes.csv";

    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::directory_iterator(benchmark_dir)) {
        if (entry.is_regular_file()) {
            files.push_back(entry.path());
        }
    }
    std::sort(files.begin(), files.end());

    std::filesystem::create_directories(output_dir);
    std::ofstream output(output_file);
    if (!output.is_open()) {
        throw std::runtime_error("Failed to open output file: " + output_file.string());
    }

    output << "instance,n,m,iterations,runtime_seconds,sigma_seconds\n";
    for (const auto& file : files) {
        const GraphSize size = ReadGraphSize(file);
        const TesterWeighted::MeasurementResult result =
            TesterWeighted::MeasureInstance(file.string(), iterations, max_time_per_instance);
        output << file.filename().string() << ','
               << size.n << ','
               << size.m << ','
               << iterations << ','
               << std::setprecision(17) << result.mean_runtime << ','
               << result.sigma_runtime << '\n';
    }

    std::cout << "saved runtimes to " << output_file << std::endl;
}

} // namespace

int main(int argc, char** argv) {
    TesterWeighted tester(true, 239);
    // tester.RunInstances(MatchingPlusGraphGenerator(6, 10, 0, 10), 100000, false);
    tester.RunInstances(MatchingPlusGraphGenerator(100, 500, -1000, 1000), 1000, false);

    const std::string benchmark_path_arg = argc > 1 ? argv[1] : "tests-weighted";

    if (benchmark_path_arg == "--export-runtimes") {
        if (argc < 3) {
            throw std::invalid_argument("Expected benchmark directory after --export-runtimes");
        }

        const std::string export_path_arg = argv[2];
        const int export_iterations = argc > 3 ? ParseIntArg(argv[3], "iterations") : 1;
        const double export_max_time_per_instance =
            argc > 4 ? ParseDoubleArg(argv[4], "max_time_per_instance") : 20.;
        ExportRuntimes(export_path_arg, export_iterations, export_max_time_per_instance);
        return 0;
    }

    const int iterations = argc > 2 ? ParseIntArg(argv[2], "iterations") : 1;
    const double max_time_per_instance = argc > 3 ? ParseDoubleArg(argv[3], "max_time_per_instance") : 20.;

    tester.MeasureBenchmark(ResolveBenchmarkPath(benchmark_path_arg).string(), iterations, max_time_per_instance);

    // tester.MeasureInstance("../tests-weighted/maxcut-big-weights/sphere-maxcut-300000-2699982", 1, 20, false);
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
