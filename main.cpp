#include <iostream>

#include "TesterWeighted.h"

int main() {
    TesterWeighted tester(true, 239);
	// tester.RunInstances(MatchingPlusGraphGenerator(8, 16, 0, 10), 100000, false);
	tester.RunInstances(MatchingPlusGraphGenerator(100, 500, -1000, 1000), 1000, false);

	tester.MeasureBenchmark("../tests-weighted", 10, 20);

	// tester.MeasureInstance("../tests-weighted/delaunay-100000-299968", 10, 20, false);
	// tester.MeasureInstance("../tests-weighted/dan59296-177299", 10, 20, false);
    // tester.MeasureInstance("../tests-weighted/sra104815-314222", 10, 20, false);
	// tester.MeasureInstance("../tests-weighted/ara238025-713594", 10, 20, false);
	// tester.MeasureInstance("../tests-weighted/lrb744710-2233725", 10, 20, false);

    return 0;
}
