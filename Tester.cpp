#include "Tester.h"


Tester::Tester(const std::vector<std::vector<int>> &adj_list_) {
    SolverUnweighted solver = SolverUnweighted(adj_list_);

    auto start = std::chrono::high_resolution_clock::now();
    solver.Solve();
    auto stop = std::chrono::high_resolution_clock::now();
    runtime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
            stop - start).count()) / 1'000'000;

    adj_list_out = std::vector<std::vector<std::tuple<int, bool, int>>>(adj_list_.size(), 
                                            std::vector<std::tuple<int, bool, int>>());
    adj_list_in = std::vector<std::vector<std::tuple<int, bool, int>>>(adj_list_.size(), 
                                            std::vector<std::tuple<int, bool, int>>());
    int edge_counter = 0;
    for (int vertex = 0; vertex < static_cast<int>(solver.adj_list.size()); ++vertex) {
        for (auto edge : solver.adj_list[vertex]) {
            adj_list_out[vertex].emplace_back(edge->OtherNode(vertex), edge->matched, edge_counter);
            adj_list_in[edge->OtherNode(vertex)].emplace_back(vertex, edge->matched, edge_counter);
            std::cout << vertex << " " << edge->OtherNode(vertex) << " " << edge->matched << " " << edge_counter << std::endl;
            edge_counter++;
        }
    }
}

bool Tester::Validate() {
    if (!IsAMatching()) {
        throw std::runtime_error("not a matching");
        return false;
    }
    if (AugmentingPathsExist()) {
        throw std::runtime_error("there is an augmenting path");
        return false;
    }
    return true;
}

bool Tester::IsAMatching() {
    for (auto neighbors : adj_list_out) {
        int matched_neighbors = 0;
        for (auto edge : neighbors) {
            if (std::get<1>(edge)) {
                matched_neighbors++;
            }
        }
        if (matched_neighbors > 1) {
            return false;
        }
    }
    return true;
}

bool Tester::AugmentingPathsExist() {
    std::vector<std::vector<int>> aux_adj_list;
    int num_terminals;
    std::tie(aux_adj_list, num_terminals) = AuxGraph();

    for (int terminal = 0; terminal < num_terminals; ++terminal) {
        std::queue<int> queue;
        queue.push(terminal);

        std::vector<bool> visited(aux_adj_list.size(), false);
        visited[terminal] = true;

        while(!queue.empty()) {
            int cur_vertex = queue.front();
            queue.pop();
            for (int next_vertex : aux_adj_list[cur_vertex]) {
                if (visited[next_vertex]) {
                    continue;
                }
                if (next_vertex < num_terminals) { // found another terminal
                    return true;
                }
                visited[next_vertex] = true;
                queue.push(next_vertex);
            }
        }
    }
    
    return false;
}

int Tester::MatchingSize() {
    int matching_size = 0;

    for (auto neighbors : adj_list_out) {
        for (auto edge : neighbors) {
            if (std::get<1>(edge)) {
                matching_size++;
            }
        }
    }

    return matching_size / 2;
}

std::vector<int> Tester::UnmatchedVertices() {
    std::vector<int> roots;

    for (int vertex = 0; vertex < static_cast<int>(adj_list_out.size()); ++vertex) {
        bool matched = false;
        for (auto edge : adj_list_out[vertex]) {
            if (std::get<1>(edge)) {
                matched = true;
                break;
            }
        }
        if (!matched) {
            roots.push_back(vertex);
        }
    }

    return roots;
}

std::pair<std::vector<std::vector<int>>, int> Tester::AuxGraph() {
    auto roots = UnmatchedVertices();
    std::vector<bool> is_root(adj_list_out.size(), false);
    for (int root : roots) {
        is_root[root] = true;
    }
    int num_terminals = static_cast<int>(roots.size());

    std::vector<std::vector<int>> aux_adj_list(num_terminals + NumberOfEdges(), std::vector<int>());

    int roots_counter = 0;
    for (int vertex = 0; vertex < static_cast<int>(adj_list_out.size()); ++vertex) {
        int matched_edge_out = -1;
        for (auto edge : adj_list_out[vertex]) {
            if (std::get<1>(edge)) {
                matched_edge_out = std::get<2>(edge);
                break;
            }
        }

        if (matched_edge_out == -1) { // vertex is unmatched, create a terminal
            for (auto edge : adj_list_out[vertex]) {
                aux_adj_list[roots_counter].push_back(std::get<2>(edge) + num_terminals);
            }
            for (auto edge : adj_list_in[vertex]) {
                aux_adj_list[std::get<2>(edge) + num_terminals].push_back(roots_counter);
            }
            roots_counter++;
        } else { // vertex is matched
            int matched_edge_in = -1;
            for (auto edge : adj_list_in[vertex]) {
                if (std::get<1>(edge)) {
                    matched_edge_in = std::get<2>(edge);
                    break;
                }
            }
            if (matched_edge_in == -1) {
                throw std::runtime_error("in AuxGraph: matched_edge_out exists but matched_edge_in does not");
            }

            for (auto edge : adj_list_out[vertex]) {
                if (!std::get<1>(edge)) {
                    aux_adj_list[num_terminals + matched_edge_in].push_back(std::get<2>(edge) + num_terminals);
                }
            }
            for (auto edge : adj_list_in[vertex]) {
                if (!std::get<1>(edge)) {
                    aux_adj_list[num_terminals + matched_edge_out].push_back(std::get<2>(edge) + num_terminals);
                }
            }
        }
    }

    std::cout << "num_terminals: " << num_terminals << std::endl;
    std::cout << "aux adj list:" << std::endl;
    for (int i = 0; i < static_cast<int>(aux_adj_list.size()); ++i) {
        std::cout << i << ": ";
        for (int j : aux_adj_list[i]) {
            std::cout << j << " ";
        }
        std::cout << std::endl;
    }

    return {aux_adj_list, num_terminals};
}

int Tester::NumberOfEdges() {
    int num_edges = 0;

    for (auto neighbors : adj_list_out) {
        num_edges += static_cast<int>(neighbors.size());
    }

    return num_edges;
}
