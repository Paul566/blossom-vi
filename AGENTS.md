# Repository Guidelines

This is a high-performance code for solving the minimum weight perfect matching problem. An (incomplete) high-level description of the algorithm is contained in `paper-unfinished.tex`.

## Project Structure & Module Organization
Core solver code lives at the repository root. `MWPMSolver.*` implements the matching algorithm, `DualUpdater.*` handles dual updates, `TesterWeighted.*` provides benchmark and verification helpers, and `main.cpp` is the local entrypoint for running benchmarks. Benchmark inputs are stored in `tests-weighted/`. Historical or reference material such as `blossom5-v2.05.src/`, `MetagraphMatchingBenchmark/`, and `paper-unfinished.tex` should be treated as supporting assets unless a task explicitly targets them.

## Solver remarks

The containers `nodes` and `trees` may contain lazily deleted nodes and trees that were destroyed earlier.

The number of remaining trees is stored in the variable `num_trees_alive`, which gets decreased upon augmentations.

`alive_trees` stores indices of trees that have not been destroyed. `alive_trees` gets updated after every primal phase.

## Build, Test, and Development Commands
Build locally with CMake:

```bash
cmake -S . -B cmake-build-release -DCMAKE_BUILD_TYPE=Release
cmake --build cmake-build-release -j2
```

Run the benchmark binary from the repo root:

```bash
./cmake-build-release/blossom_vi tests-weighted 1 20
```

Arguments are `path`, `iterations`, and `max_time_per_instance`. Use a one-file directory in `/tmp` for quick A/B timing when testing optimizations.

## Coding Style & Naming Conventions
Use C++20 and Google styleguide. You may add short, direct comments only where the code is not self-explanatory.

## Testing Guidelines
Validation is done through property based testing by including/uncommenting `tester.RunInstances(MatchingPlusGraphGenerator(100, 500, -1000, 1000), 1000, false);` in main. It validates the feasibility of the primal and dual solutions as well as the complementary slackness conditions, thus validating the optimality of the output. For solver changes, rebuild and run the random tests.

For testing the performance, use the benchmark in `tests-weighted`. Compare before/after timings on the same input when claiming a speedup.

## Commit & Pull Request Guidelines
Recent commits use short, imperative subjects such as `avoid extra container allocations` and `remove old files`. Keep commit titles concise and action-oriented. Pull requests should explain the algorithmic or performance motivation, list the benchmark commands used, and include before/after timings for any optimization-related change.
