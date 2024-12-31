#include <iostream>
#include <vector>
#include "Node.h"
#include "SolverUnweighted.h"

void PrintVector(std::vector<int> numbers) {
    for (int number : numbers) {
        std::cout << number << " ";
    }
    std::cout << std::endl;
}

int main() {
    std::vector<std::vector<int>> adj_list;
    adj_list.push_back(std::vector<int>({1, 2, 3}));
    adj_list.push_back(std::vector<int>({0, 2, 3}));
    adj_list.push_back(std::vector<int>({0, 1}));
    adj_list.push_back(std::vector<int>({0, 1}));

    SolverUnweighted solver = SolverUnweighted(adj_list);
    solver.Solve();
    solver.Print();

    return 0;
}