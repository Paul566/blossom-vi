#include <iostream>

#include "TesterWeighted.h"
// #include "Heap.h"

int main() {
    TesterWeighted tester(true, 239);
	// tester.RunInstances(MatchingPlusGraphGenerator(6, 10, 0, 10), 10000, true);
	tester.RunInstances(MatchingPlusGraphGenerator(100, 500, -1000, 1000), 1000, false);
	// tester.MeasureBenchmark("../tests-weighted", 1, 30);
	// tester.MeasureInstance("../tests-weighted/delaunay-100000-299968", 1);
    tester.MeasureInstance("../tests-weighted/sra104815-314222", 1);

    return 0;
}
