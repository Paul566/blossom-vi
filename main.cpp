#include <iostream>

#include "TesterWeighted.h"

int main() {
    TesterWeighted tester(true, 239);
	// tester.RunInstances(MatchingPlusGraphGenerator(8, 16, 0, 10), 100000, false);
	// tester.RunInstances(MatchingPlusGraphGenerator(100, 500, -1000, 1000), 1000, false);
	tester.MeasureBenchmark("../tests-weighted", 1, 30);
	// tester.MeasureInstance("../tests-weighted/delaunay-100000-299968", 1, 60, false);
    // tester.MeasureInstance("../tests-weighted/sra104815-314222", 1, 60, false);
	// tester.MeasureInstance("../tests-weighted/ara238025-713594", 1, 60, false);
	// tester.MeasureInstance("../tests-weighted/lrb744710-2233725", 1, 60, false);

    return 0;
}
