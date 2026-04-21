# Blossom VI

Blossom VI is a C++ implementation of a **minimum weight perfect matching** solver. Given an undirected weighted graph, the solver finds a perfect matching of minimum total edge weight. If the graph does not contain a perfect matching, an exception is thrown.

## Build

Build:

```bash
cmake -S . -B cmake-build-release -DCMAKE_BUILD_TYPE=Release
cmake --build cmake-build-release -j2
```

The executable is written to `cmake-build-release/blossom_vi`.

Requirements:

- CMake 3.30 or newer
- A C++20 compiler
- Boost headers

## Solve an instance

Run the solver on a graph file:

```bash
./cmake-build-release/blossom_vi examples/example.txt
```

Output:

```text
Optimal weight: 7
```

To save the matching to a file `matching.txt`:

```bash
./cmake-build-release/blossom_vi examples/example.txt matching.txt
```

The matching file contains one matched edge per line:

```text
0 2
1 3
```

By default, stdout contains only the optimal matching weight. To also print the matching:

```bash
./cmake-build-release/blossom_vi --print-matching examples/example.txt
```

To also print the time spent in the main solver function call:

```bash
./cmake-build-release/blossom_vi --print-runtime examples/example.txt
```

## Input format

Graph instances use the following format.

The first line contains:

```text
n m
```

where `n` is the number of vertices and `m` is the number of edges.

The next `m` lines contain one undirected weighted edge per line:

```text
endpoint1 endpoint2 weight
```

Vertices are numbered from `0` to `n - 1`. Weights are integers and may be negative. The absolute values of the edge weights must not exceed $2^{29}-1$, and the difference between the maximal and the minimal weight also must not exceed this number, otherwise there might be overflows.

Example:

```text
4 5
0 1 10
0 2 3
1 3 4
2 3 7
1 2 8
```

## Use from C++

The public solver class is `MWPMSolver`. Construct it from a vector of `(endpoint1, endpoint2, weight)` tuples, run `FindMinPerfectMatching()`, then read the matching:

```cpp
#include <tuple>
#include <utility>
#include <vector>

#include "MWPMSolver.h"

int main() {
    std::vector<std::tuple<int, int, int>> edges = {
        {0, 1, 10},
        {0, 2, 3},
        {1, 3, 4},
        {2, 3, 7},
        {1, 2, 8},
    };

    MWPMSolver solver(edges);
    solver.FindMinPerfectMatching();

    const std::vector<std::pair<int, int>>& matching = solver.Matching();
    const int64_t optimal_weight = solver.primal_objective;
}
```

## Validation and testing

There is a correctness test that runs randomized instances. The correctness of the output is validated through verification of the primal feasibility, the dual feasibility, and the complementary slackness conditions.

```bash
./cmake-build-release/blossom_vi_test
```

You can also choose the random graph size, number of iterations, and seed:

```bash
./cmake-build-release/blossom_vi_test 100 500 1000 67
```

The arguments are `num_vertices`, `num_edges`, `iterations`, and `seed`.
`num_vertices` must be even, and `num_edges` must be at least
`num_vertices / 2`.
