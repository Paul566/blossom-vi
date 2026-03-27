#include <format>
#include <iomanip>
#include <iostream>

#include "TesterWeighted.h"

int main() {
    TesterWeighted tester(true, 239);
    // tester.RunInstances(MatchingPlusGraphGenerator(6, 10, 0, 10), 100000, false);
    tester.RunInstances(MatchingPlusGraphGenerator(100, 500, -1000, 1000), 1000, false);

    int iterations = 1;
    tester.MeasureBenchmark("../tests-weighted", iterations, 20);

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
