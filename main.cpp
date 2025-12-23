#include <iostream>

#include "TesterWeighted.h"

int main() {
    // std::mt19937 gen(239);
    TesterWeighted tester(true, 239);
    // tester.RunInstances(CliqueGenerator(4, 0, 10), 10000, true);
	// tester.RunInstances(CliqueGenerator(100, -1000, 1000), 1000, false);
    // tester.RunInstances(MatchingPlusGraphGenerator(10000, 60000, -1000, 1000), 10, false);
	// tester.MeasureBenchmark("../tests-weighted", 1, 30);
	// tester.MeasureInstance("../tests-weighted/delaunay-100000-299968", 1);
    tester.MeasureInstance("../tests-weighted/sra104815-314222", 1);

    return 0;
}
