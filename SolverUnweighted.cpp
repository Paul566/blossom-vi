#include "SolverUnweighted.h"


SolverUnweighted::SolverUnweighted(const std::vector<std::vector<int>> &adj_list_) : 
                                                n(static_cast<int>(adj_list_.size())){
    adj_list = std::vector<std::vector<std::shared_ptr<Edge>>>(n);
    matched_edge = std::vector<std::shared_ptr<Edge>>(n, nullptr);

    for (int i = 0; i < n; ++i) {
        adj_list[i].reserve(n);
        for (int to: adj_list_[i]) {
            if (to > i) {
                Edge new_edge(i, to);
                std::shared_ptr<Edge> edge_ptr = std::make_shared<Edge>(new_edge);
                adj_list[i].push_back(edge_ptr);
                adj_list[to].push_back(edge_ptr);
            }
        }
    }
}

void SolverUnweighted::Print() {
    std::cout << "Matching:" << std::endl;
    for (int vertex = 0; vertex < n; ++vertex) {
        for (auto edge : adj_list[vertex]) {
            if (edge->matched) {
                if (vertex < edge->OtherNode(vertex)) {
                    std::cout << vertex << " " << edge->OtherNode(vertex) << std::endl;
                }
            }
        }
    }

    std::cout << "\nAdjacency list:" << std::endl;
    for (int vertex = 0; vertex < n; ++vertex) {
        std::cout << vertex << ": ";
        for (auto edge : adj_list[vertex]) {
            std::cout << edge->OtherNode(vertex) << " ";
        }
        std::cout << std::endl;
    }
}

void SolverUnweighted::Solve() {
    GreedyInit();

    for (int root = 0; root < n; ++root) {
        if (matched_edge[root]) {
            continue;
        }

        std::queue<int> queue;
        queue.push(root);

        std::vector<std::shared_ptr<Edge>> plus_parents(n, nullptr);
        std::vector<std::shared_ptr<Edge>> minus_parents(n, nullptr);
        std::vector<bool> plus(n, false);
        std::vector<bool> minus(n, false);
        std::vector<int> depth_plus(n, -1);
        std::vector<int> depth_minus(n, -1);
        plus[root] = true;
        depth_plus[root] = 0;
        bool augmented = false;

        while((!queue.empty()) && (!augmented)) {
            int cur_vertex = queue.front();
            queue.pop();

            for (auto edge : adj_list[cur_vertex]) {
                int to = edge->OtherNode(cur_vertex);

                if (to == root) {
                    continue;
                }

                if ((!plus[to]) && (!minus[to])) { // to is not in the tree yet
                    if (matched_edge[to]) { // if to is matched, add it
                        minus[to] = true;
                        plus[matched_edge[to]->OtherNode(to)] = true;
                        minus_parents[to] = edge;
                        plus_parents[matched_edge[to]->OtherNode(to)] = matched_edge[to];
                        depth_minus[to] = depth_plus[cur_vertex] + 1;
                        depth_plus[matched_edge[to]->OtherNode(to)] = depth_plus[cur_vertex] + 2;
                        queue.push(matched_edge[to]->OtherNode(to));
                    } else { // if to is unmatched, augment
                        std::vector<std::shared_ptr<Edge>> path;
                        path.push_back(edge);
                        int path_vertex = cur_vertex;
                        bool now_plus = true;
                        while (path_vertex != root) {
                            if (now_plus) {
                                path.push_back(plus_parents[path_vertex]);
                                path_vertex = plus_parents[path_vertex]->OtherNode(path_vertex);
                                now_plus = false;
                            } else {
                                path.push_back(minus_parents[path_vertex]);
                                path_vertex = minus_parents[path_vertex]->OtherNode(path_vertex);
                                now_plus = true;
                            }
                        }

                        Augment(path, to);
                        augmented = true;
                        break;
                    }
                } else { // to is in the tree
                    if (!minus[to]) { // we will create a cherry blossom
                        UpdateLabels(edge, root, queue, plus_parents, minus_parents, 
                                    plus, minus, depth_plus, depth_minus);
                    }
                }
            }
        }
    }
}

void SolverUnweighted::GreedyInit() {
    for (int vertex = 0; vertex < n; ++vertex) {
        if (!matched_edge[vertex]) {
            for (auto edge : adj_list[vertex]) {
                int to = edge->OtherNode(vertex);
                if (!matched_edge[to]) {
                    edge->matched = true;
                    matched_edge[vertex] = edge;
                    matched_edge[to] = edge;
                    break;
                }
            }
        }
    }
}

void SolverUnweighted::Augment(std::vector<std::shared_ptr<Edge>> path, int start) {
    if (path.size() % 2 == 0) {
        throw std::runtime_error("in Augment: the length of the path must be odd");
    }
    if (matched_edge[start]) {
        throw std::runtime_error("in Augment: start must be unmatched");
    }

    for (int i = 0; i < static_cast<int>(path.size()); ++i) {
        if ((i % 2 == 0) && (path[i]->matched)) {
            throw std::runtime_error("in Augment: not an alternating path");
        }
        if ((i % 2 == 1) && (!path[i]->matched)) {
            throw std::runtime_error("in Augment: not an alternating path");
        }
    }

    int cur_vertex = start;
    for (int i = 0; i < static_cast<int>(path.size()); ++i) {
        int next_vertex = path[i]->OtherNode(cur_vertex);

        if (i == static_cast<int>(path.size()) - 1) {
            if (matched_edge[next_vertex]) {
                throw std::runtime_error("in Augment: end must be unmatched");
            }
        }

        if (i % 2 == 0) {
            matched_edge[cur_vertex] = path[i];
            matched_edge[next_vertex] = path[i];
        }
        path[i]->matched = !path[i]->matched;
        cur_vertex = next_vertex;
    }
}

void SolverUnweighted::UpdateLabels(std::shared_ptr<Edge> edge_plus_plus, int root, 
    std::queue<int> &queue, 
    std::vector<std::shared_ptr<Edge>> &plus_parents, 
    std::vector<std::shared_ptr<Edge>> &minus_parents, 
    std::vector<bool> &plus, std::vector<bool> &minus,
    std::vector<int> &depth_plus,
    std::vector<int> &depth_minus) {

    int first_vertex, second_vertex = -1;
    std::tie(first_vertex, second_vertex) = edge_plus_plus->Vertices();

    minus_parents[first_vertex] = edge_plus_plus;
    minus_parents[second_vertex] = edge_plus_plus;
    depth_minus[first_vertex] = depth_plus[second_vertex] + 1;
    depth_minus[second_vertex] = depth_plus[first_vertex] + 1;

    while(first_vertex != second_vertex) {
        // std::cout << "while in UpdLabels:" << first_vertex << " " << second_vertex << std::endl;

        int cur_vertex = first_vertex;
        if (depth_plus[first_vertex] < depth_plus[second_vertex]) {
            cur_vertex = second_vertex;
        } 

        minus[cur_vertex] = true;

        int parent = plus_parents[first_vertex]->OtherNode(first_vertex);
        if (!plus[parent]) {
            plus[parent] = true;
            queue.push(parent);
            plus_parents[parent] = plus_parents[first_vertex];
            depth_plus[parent] = depth_minus[first_vertex] + 1;
        }

        int grandparent = minus_parents[parent]->OtherNode(parent);
        if (!minus[grandparent]) {
            minus[grandparent] = true;
            minus_parents[grandparent] = minus_parents[parent];
            depth_minus[grandparent] = depth_plus[parent] + 1;
        }

        if (cur_vertex == first_vertex) {
            first_vertex = grandparent;
        } else {
            second_vertex = grandparent;
        }
    }
}
