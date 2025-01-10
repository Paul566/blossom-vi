#include "Tester.h"

#include <chrono>
#include "SolverUnweighted.h"

Tester::Tester(const std::vector<std::vector<int> > &adj_list_,
               int greedy_init_type_,
               bool verbose_) : verbose(verbose_), greedy_init_type(greedy_init_type_),
                                adj_list(adj_list_) {

    SolverUnweighted solver = SolverUnweighted(adj_list_, verbose_, greedy_init_type_);

    const auto start = std::chrono::high_resolution_clock::now();
    solver.Solve();
    const auto stop = std::chrono::high_resolution_clock::now();
    runtime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
        stop - start).count()) / 1'000'000;

    matched_to = std::vector<int>(adj_list_.size(), -1);
    blossom_index.reserve(adj_list_.size());
    for (int i = 0; i < static_cast<int>(adj_list_.size()); ++i) {
        blossom_index.push_back(i);
    }

    for (int vertex = 0; vertex < static_cast<int>(solver.adj_list.size()); ++vertex) {
        for (const auto &edge : solver.adj_list[vertex]) {
            if (edge->matched) {
                if (matched_to[vertex] != -1) {
                    throw std::runtime_error("Not a matching");
                }
                matched_to[vertex] = edge->OtherNode(vertex);
            }
        }
    }
}

bool Tester::Validate() {
    if (AugmentingPathsExist()) {
        throw std::runtime_error("there is an augmenting path");
        return false;
    }
    return true;
}

int Tester::LCA(int first_vertex, int second_vertex, const std::vector<int> &parents) const {
    std::vector<bool> visited_first(adj_list.size(), false);
    visited_first[first_vertex] = true;

    while (parents[first_vertex] != -1) {
        first_vertex = parents[first_vertex];
        if (first_vertex == second_vertex) {
            return first_vertex;
        }
        visited_first[first_vertex] = true;
    }

    while (parents[second_vertex] != -1) {
        second_vertex = parents[second_vertex];
        if (visited_first[second_vertex]) {
            return second_vertex;
        }
    }

    return first_vertex;
}

bool Tester::AugmentingPathsExist() {
    for (int root = 0; root < static_cast<int>(adj_list.size()); ++root) {
        if (matched_to[root] != -1) {
            continue;
        }

        std::queue<int> queue;
        queue.push(root);

        std::vector<bool> visited(adj_list.size(), false);
        std::vector<bool> plus(adj_list.size(), false);
        std::vector<int> parents(adj_list.size(), -1);
        visited[root] = true;
        plus[root] = true;

        while (!queue.empty()) {
            int cur_vertex = queue.front(); // is a plus
            cur_vertex = TopBlossom(cur_vertex);
            queue.pop();
            for (int next_vertex : adj_list[cur_vertex]) {
                next_vertex = TopBlossom(next_vertex);
                if (cur_vertex == next_vertex) {
                    continue;
                }
                if (!visited[next_vertex]) {
                    if (matched_to[next_vertex] == -1) {
                        return true;
                    }
                    visited[next_vertex] = true;
                    parents[next_vertex] = cur_vertex;
                    int next_plus = matched_to[next_vertex];
                    visited[next_plus] = true;
                    parents[next_plus] = next_vertex;
                    plus[next_plus] = true;
                    queue.push(next_plus);
                } else {
                    if (!plus[next_vertex]) {
                        continue;
                    }

                    std::vector<int> future_blossom;
                    int lca = LCA(cur_vertex, next_vertex, parents);
                    int cur_end = cur_vertex;
                    while (cur_end != lca) {
                        future_blossom.push_back(cur_end);
                        cur_end = parents[cur_end];
                    }
                    while (next_vertex != lca) {
                        future_blossom.push_back(next_vertex);
                        next_vertex = parents[next_vertex];
                    }
                    future_blossom.push_back(lca);

                    MakeBlossom(future_blossom, parents);

                    int new_blossom = static_cast<int>(adj_list.size()) - 1;
                    visited.push_back(true);
                    plus.push_back(true);
                    queue.push(new_blossom);

                    cur_vertex = new_blossom;
                }
            }
        }
    }

    return false;
}

void Tester::MakeBlossom(const std::vector<int> &blossom, std::vector<int> &parents) {
    // the receptacle is the last element of blossom

    std::vector<bool> in_blossom(adj_list.size(), false);
    for (int vertex : blossom) {
        in_blossom[vertex] = true;
    }

    int back_index = static_cast<int>(adj_list.size());
    adj_list.emplace_back();

    matched_to.emplace_back(matched_to[blossom.back()]);
    if (matched_to[blossom.back()] != -1) {
        matched_to[matched_to[blossom.back()]] = back_index;
    }

    for (int vertex = 0; vertex < static_cast<int>(adj_list.size()) - 1; ++vertex) {
        for (int i = 0; i < static_cast<int>(adj_list[vertex].size()); ++i) {
            if ((!in_blossom[vertex]) && (in_blossom[adj_list[vertex][i]])) {
                adj_list[vertex][i] = back_index;
                adj_list.back().emplace_back(vertex);
            }
        }
    }

    for (const int vertex : blossom) {
        adj_list[vertex].clear();
    }

    for (int vertex = 0; vertex < static_cast<int>(adj_list.size()) - 1; ++vertex) {
        if (parents[vertex] != -1) {
            if (in_blossom[parents[vertex]]) {
                parents[vertex] = back_index;
            }
        }
    }
    parents.emplace_back(parents[blossom.back()]);

    for (const int vertex : blossom) {
        blossom_index[vertex] = back_index;
    }
    blossom_index.emplace_back(back_index);
}

int Tester::TopBlossom(int vertex) const {
    while (blossom_index[vertex] != vertex) {
        vertex = blossom_index[vertex];
    }
    return vertex;
}
