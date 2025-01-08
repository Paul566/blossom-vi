#include "SolverUnweighted.h"

#include <iostream>
#include <unordered_set>

SolverUnweighted::SolverUnweighted(const std::vector<std::vector<int> > &adj_list_) : n(static_cast<int>(adj_list_.
    size())), cherry_blossoms(LabeledDisjointSets(n)) {
    adj_list = std::vector<std::vector<std::shared_ptr<Edge> > >(n);
    matched_edge = std::vector<std::shared_ptr<Edge> >(n, nullptr);

    for (int i = 0; i < n; ++i) {
        adj_list[i].reserve(n);
        for (int to : adj_list_[i]) {
            if (to > i) {
                Edge new_edge(i, to);
                std::shared_ptr<Edge> edge_ptr = std::make_shared<Edge>(new_edge);
                adj_list[i].push_back(edge_ptr);
                adj_list[to].push_back(edge_ptr);
            }
        }
    }

    minus_parents = std::vector<std::shared_ptr<Edge> >(n, nullptr);
    plus = std::vector<bool>(n, false);
    minus = std::vector<bool>(n, false);
    growable_vertices = std::queue<int>();
    children = std::vector<std::vector<int> >(n, std::vector<int>());
}

void SolverUnweighted::PrintMatching() {
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
}

void SolverUnweighted::PrintAdjList() {
    std::cout << "Adjacency list:" << std::endl;
    for (int vertex = 0; vertex < n; ++vertex) {
        std::cout << vertex << ": ";
        for (auto edge : adj_list[vertex]) {
            std::cout << edge->OtherNode(vertex) << " ";
        }
        std::cout << std::endl;
    }
}

void SolverUnweighted::PrintTreeData() {
    std::cout << "Tree data:\n";
    for (int v = 0; v < n; ++v) {
        std::cout << v << ": ";
        if (plus[v]) {
            if (matched_edge[v]) {
                std::cout << "(+ parent: " << matched_edge[v]->OtherNode(v) << ") ";
            } else {
                std::cout << "(is root) ";
            }
        }
        if (minus[v]) {
            std::cout << "(- parent: " << minus_parents[v]->OtherNode(v) << ") ";
        }
        std::cout << "(receptacle: " << cherry_blossoms.Label(v) << ")";
        std::cout << std::endl;
    }
}

void SolverUnweighted::Solve() {
    GreedyInit();
    // PrintAdjList();
    // std::cout << "after greedy init:\n";
    // PrintMatching();
    // std::cout << std::endl;

    for (int root = 0; root < n; ++root) {
        if (matched_edge[root]) {
            continue;
        }

        growable_vertices.push(root);
        plus[root] = true;
        bool augmented = false;

        while ((!growable_vertices.empty()) && (!augmented)) {
            int cur_vertex = growable_vertices.front();
            growable_vertices.pop();
            if ((!plus[cur_vertex]) && (!minus[cur_vertex])) {
                // the vertex is in the queue but its tree was deleted earlier
                continue;;
            }
            // std::cout << "root: " << root << " " << cur_vertex << std::endl;
            // for (int v = 0; v < n; ++v) {
            //     std::cout << v << " " << plus[v] << " " << minus[v] << std::endl;
            // }

            for (auto edge : adj_list[cur_vertex]) {
                int to = edge->OtherNode(cur_vertex);

                if (to == root) {
                    continue;
                }

                if ((!plus[to]) && (!minus[to])) {
                    // to is not in the tree yet
                    if (matched_edge[to]) {
                        // if to is matched, add it
                        minus[to] = true;
                        plus[matched_edge[to]->OtherNode(to)] = true;
                        minus_parents[to] = edge;
                        children[cur_vertex].push_back(to);
                        children[to].push_back(matched_edge[to]->OtherNode(to));
                        growable_vertices.push(matched_edge[to]->OtherNode(to));
                    } else {
                        // if to is unmatched, augment
                        std::vector<std::shared_ptr<Edge> > path;
                        path.push_back(edge);
                        int path_vertex = cur_vertex;
                        bool now_plus = true;
                        while (path_vertex != root) {
                            if (now_plus) {
                                path.push_back(matched_edge[path_vertex]);
                                path_vertex = matched_edge[path_vertex]->OtherNode(path_vertex);
                                now_plus = false;
                            } else {
                                path.push_back(minus_parents[path_vertex]);
                                path_vertex = minus_parents[path_vertex]->OtherNode(path_vertex);
                                now_plus = true;
                            }
                            // std::cout << path_vertex << " ";
                        }
                        // std::cout << std::endl;

                        Augment(path, to);
                        augmented = true;
                        ClearTree(root);
                        break;
                    }
                } else {
                    // to is in the tree
                    if (plus[to]) {
                        // we will create a cherry blossom
                        MakeCherryBlossom(edge, root);
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

void SolverUnweighted::Augment(std::vector<std::shared_ptr<Edge> > path, int start) {
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

void SolverUnweighted::MakeCherryBlossom(std::shared_ptr<Edge> edge_plus_plus, int root) {
    int first_vertex, second_vertex = -1;
    std::tie(first_vertex, second_vertex) = edge_plus_plus->Vertices();
    if (first_vertex == second_vertex) {
        throw std::runtime_error("in UpdateLabels: first_vertex == second_vertex");
    }
    if (cherry_blossoms.Representative(first_vertex) == cherry_blossoms.Representative(second_vertex)) {
        return;
    }

    // std::cout << "making cherry blossom " << first_vertex << " " << second_vertex << std::endl;
    // PrintTreeData();

    int first_bound, second_bound;
    std::tie(first_bound, second_bound) = PathUpperBounds(first_vertex, second_vertex);
    // std::cout << "bounds: " << first_bound << " " << second_bound << std::endl;

    UpdatePath(first_vertex, first_bound);
    UpdatePath(second_vertex, second_bound);

    if (first_vertex != first_bound) {
        // first_vertex is not a root
        minus[first_vertex] = true;
        minus_parents[first_vertex] = edge_plus_plus;
    }
    if (second_vertex != second_bound) {
        // second_vertex is not a root
        minus[second_vertex] = true;
        minus_parents[second_vertex] = edge_plus_plus;
    }
}

int SolverUnweighted::PlusPlusLCA(int first_vertex, int second_vertex) const {
    std::unordered_set<int> visited_first;
    visited_first.insert(first_vertex);

    while (matched_edge[first_vertex]) {
        first_vertex = matched_edge[first_vertex]->OtherNode(first_vertex);
        first_vertex = minus_parents[first_vertex]->OtherNode(first_vertex);
        if (first_vertex == second_vertex) {
            return first_vertex;
        }
        visited_first.insert(first_vertex);
    }

    while (matched_edge[second_vertex]) {
        second_vertex = matched_edge[second_vertex]->OtherNode(second_vertex);
        second_vertex = minus_parents[second_vertex]->OtherNode(second_vertex);
        if (visited_first.contains(second_vertex)) {
            return second_vertex;
        }
    }

    return first_vertex;
}

std::pair<int, int> SolverUnweighted::PathUpperBounds(int first_vertex, int second_vertex) {
    int lca = PlusPlusLCA(first_vertex, second_vertex);
    int lca_receptacle = cherry_blossoms.Label(lca);

    while (first_vertex != lca) {
        if (cherry_blossoms.Label(first_vertex) == lca_receptacle) {
            break;
        }
        first_vertex = matched_edge[first_vertex]->OtherNode(first_vertex);
        first_vertex = minus_parents[first_vertex]->OtherNode(first_vertex);
    }
    while (second_vertex != lca) {
        if (cherry_blossoms.Label(second_vertex) == lca_receptacle) {
            break;
        }
        second_vertex = matched_edge[second_vertex]->OtherNode(second_vertex);
        second_vertex = minus_parents[second_vertex]->OtherNode(second_vertex);
    }

    return {first_vertex, second_vertex};
}

void SolverUnweighted::UpdatePath(int lower_vertex, int upper_vertex) {
    int receptacle = cherry_blossoms.Label(upper_vertex);

    while (lower_vertex != upper_vertex) {
        // std::cout << lower_vertex << std::endl;

        int parent = matched_edge[lower_vertex]->OtherNode(lower_vertex);
        int grandparent = minus_parents[parent]->OtherNode(parent);

        cherry_blossoms.Unite(lower_vertex, upper_vertex, receptacle);
        cherry_blossoms.Unite(parent, upper_vertex, receptacle);

        if (!plus[parent]) {
            plus[parent] = true;
            growable_vertices.push(parent);
        }

        if (grandparent != upper_vertex) {
            minus[grandparent] = true;
            minus_parents[grandparent] = minus_parents[parent];
        }

        lower_vertex = grandparent;
    }
}

void SolverUnweighted::ClearTree(int root) {
    std::vector<int> tree_vertices;
    std::queue<int> queue;
    queue.push(root);
    while (!queue.empty()) {
        int cur_vertex = queue.front();
        queue.pop();
        tree_vertices.push_back(cur_vertex);
        for (int child : children[cur_vertex]) {
            queue.push(child);
        }
    }

    for (int vertex : tree_vertices) {
        minus_parents[vertex] = nullptr;
        plus[vertex] = false;
        minus[vertex] = false;
        children[vertex].clear();
        cherry_blossoms.Detach(vertex);
    }
    // TODO some of the vertices may still be in the growable_vertices queue
    // growable_vertices = std::queue<int>();
}
