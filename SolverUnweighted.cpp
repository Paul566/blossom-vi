#include "SolverUnweighted.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <utility>

SolverUnweighted::SolverUnweighted(const std::vector<std::vector<int> > &adj_list_,
                                   int greedy_init_type_,
                                   bool delete_edges_in_cherries_,
                                   bool verbose_) : n(static_cast<int>(adj_list_.size())),
                                                    cherry_blossoms(LabeledDisjointSets(n)), verbose(verbose_),
                                                    greedy_init_type(greedy_init_type_),
                                                    delete_edges_in_cherries(delete_edges_in_cherries_) {
    adj_list = std::vector<std::vector<std::shared_ptr<Edge> > >(n);
    matched_edge = std::vector<std::shared_ptr<Edge> >(n, nullptr);

    for (int i = 0; i < n; ++i) {
        adj_list[i].reserve(adj_list_[i].size());
    }

    for (int i = 0; i < n; ++i) {
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
    growable_vertices = std::queue<int>();
    _root_of_vertex = std::vector<int>(n, -1);
    tree_plus_neighbors = std::vector<std::vector<int> >(n, std::vector<int>());
    lca_markers = std::vector<int>(n, -1);
    lca_count = 0;
}

std::vector<std::pair<int, int> > SolverUnweighted::Matching() {
    std::vector<std::pair<int, int> > result;
    result.reserve(n / 2);
    for (int vertex = 0; vertex < n; ++vertex) {
        if (matched_edge[vertex]) {
            int other_vertex = matched_edge[vertex]->OtherNode(vertex);
            if (vertex < other_vertex) {
                result.emplace_back(vertex, other_vertex);
            }
        }
    }
    return result;
}

void SolverUnweighted::PrintMatching() {
    std::cout << "Matching:" << std::endl;

    auto matching = Matching();
    for (auto edge : matching) {
        std::cout << edge.first << " " << edge.second << std::endl;
    }
}

void SolverUnweighted::PrintAdjList() {
    std::cout << "Adjacency list:" << std::endl;
    for (int vertex = 0; vertex < n; ++vertex) {
        std::cout << vertex << ": ";
        for (const auto& edge : adj_list[vertex]) {
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
        if (minus_parents[v]) {
            std::cout << "(- parent: " << minus_parents[v]->OtherNode(v) << ") ";
        }
        std::cout << "(receptacle: " << cherry_blossoms.Label(v) << ") ";
        std::cout << "(root: " << RootOfVertex(v) << ") ";
        std::cout << "(is a plus: " << plus[v] << ") ";
        std::cout << std::endl;
    }

    std::cout << "queue: ";
    auto queue_copy = growable_vertices;
    while (!queue_copy.empty())
    {
        std::cout << queue_copy.front() << " ";
        queue_copy.pop();
    }
    std::cout << std::endl;
}

void SolverUnweighted::Solve() {
    GreedyInit();

    if (verbose) {
        PrintAdjList();
        PrintMatching();
    }

    for (int root = 0; root < n; ++root) {
        if (!matched_edge[root]) {
            growable_vertices.push(root);
            plus[root] = true;;
            _root_of_vertex[root] = root;
        }
    }

    while (!growable_vertices.empty()) {
        int cur_vertex = growable_vertices.front();
        if (verbose) {
            PrintAdjList();
            PrintMatching();
            PrintTreeData();
            std::cout << "cur_vertex: " << cur_vertex << std::endl;
        }
        growable_vertices.pop();
        HandleVertex(cur_vertex);
    }
}

void SolverUnweighted::GreedyInit() {
    if (greedy_init_type == 0) {
        for (int vertex = 0; vertex < n; ++vertex) {
            if (!matched_edge[vertex]) {
                for (const auto &edge : adj_list[vertex]) {
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
        return;
    }

    if (greedy_init_type == 1) {
        std::vector<std::pair<int, int> > vertices;
        vertices.reserve(n);
        for (int i = 0; i < n; ++i) {
            vertices.emplace_back(adj_list[i].size(), i);
        }
        std::sort(vertices.begin(), vertices.end());

        for (int i = 0; i < n; ++i) {
            int vertex = vertices[i].second;
            if (!matched_edge[vertex]) {
                for (const auto &edge : adj_list[vertex]) {
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
        return;
    }

    if (greedy_init_type == 2) {
        std::vector<std::pair<int, std::shared_ptr<Edge> > > edges;
        for (int i = 0; i < n; ++i) {
            for (const auto &edge : adj_list[i]) {
                if (i < edge->OtherNode(i)) {
                    edges.emplace_back(adj_list[i].size() + adj_list[edge->OtherNode(i)].size(), edge);
                }
            }
        }
        std::sort(edges.begin(), edges.end());
        for (auto &edge : edges) {
            auto [head, tail] = edge.second->Vertices();
            if ((!matched_edge[head]) && (!matched_edge[tail])) {
                edge.second->matched = true;
                matched_edge[head] = edge.second;
                matched_edge[tail] = edge.second;
            }
        }
        return;
    }

    if (greedy_init_type == 3) {
        std::vector<std::pair<int, std::shared_ptr<Edge> > > edges;
        for (int i = 0; i < n; ++i) {
            for (const auto &edge : adj_list[i]) {
                if (i < edge->OtherNode(i)) {
                    edges.emplace_back(std::max(adj_list[i].size(), adj_list[edge->OtherNode(i)].size()), edge);
                }
            }
        }
        std::sort(edges.begin(), edges.end());
        for (auto &edge : edges) {
            auto [head, tail] = edge.second->Vertices();
            if ((!matched_edge[head]) && (!matched_edge[tail])) {
                edge.second->matched = true;
                matched_edge[head] = edge.second;
                matched_edge[tail] = edge.second;
            }
        }
    }
}

bool SolverUnweighted::HandleVertex(const int cur_vertex) {
    if (!plus[cur_vertex]) {
        return false;
    }
    int cur_root = RootOfVertex(cur_vertex);
    if (cur_root == -1) {
        // the vertex is in the queue but its tree was deleted earlier
        return false;
    }

    for (int i = 0; i < static_cast<int>(adj_list[cur_vertex].size()); ++i) {
        // auto edge = adj_list[cur_vertex][i];
        int to = adj_list[cur_vertex][i]->OtherNode(cur_vertex);
        int to_root = RootOfVertex(to);

        if ((to_root == -1) && (matched_edge[to])) {
            // to is not in the tree yet
            AddVertex(to, cur_vertex, adj_list[cur_vertex][i]);
        } else {
            // to is in the tree
            if (plus[to]) {
                if (cur_root == to_root) {
                    if (!cherry_blossoms.SameSet(cur_vertex, to)) {
                        MakeCherryBlossom(adj_list[cur_vertex][i]);
                    } else {
                        if ((delete_edges_in_cherries) && (matched_edge[cur_vertex] != adj_list[cur_vertex][i]) && (minus_parents[
                            cur_vertex] != adj_list[cur_vertex][i]) && (minus_parents[to] != adj_list[cur_vertex][i])) {
                            // if the edge is not anyone's parent, delete it
                            adj_list[cur_vertex][i] = adj_list[cur_vertex].back();
                            adj_list[cur_vertex].pop_back();
                            --i;
                        }
                    }
                } else {
                    Augment(adj_list[cur_vertex][i], cur_vertex, to);
                    return true;
                }
            } else {
                tree_plus_neighbors[to_root].push_back(cur_vertex);
            }
        }
    }

    return false;
}

void SolverUnweighted::Augment(const std::shared_ptr<Edge> &edge_plus_plus, int cur_vertex, int to) {
    auto [first_vertex, second_vertex] = edge_plus_plus->Vertices();
    auto first_path = PathToRoot(first_vertex);
    auto second_path = PathToRoot(second_vertex);

    std::vector<std::shared_ptr<Edge> > path;
    path.reserve(first_path.size() + second_path.size() + 1);
    for (int i = static_cast<int>(first_path.size()) - 1; i >= 0; --i) {
        path.push_back(first_path[i]);
    }
    path.push_back(edge_plus_plus);
    for (const auto &edge : second_path) {
        path.push_back(edge);
    }

    int first_root = UnmatchedVertex(path.front());
    int second_root = UnmatchedVertex(path.back());

    AugmentPath(path);
    ClearTree(first_root, second_root);
    ClearTree(second_root, first_root);
}

std::vector<std::shared_ptr<Edge> > SolverUnweighted::PathToRoot(int vertex_plus) {
    if (!matched_edge[vertex_plus]) {
        // already in a root
        return {};
    }

    std::vector<std::shared_ptr<Edge> > path;
    while (matched_edge[vertex_plus]) {
        path.push_back(matched_edge[vertex_plus]);
        vertex_plus = matched_edge[vertex_plus]->OtherNode(vertex_plus);
        path.push_back(minus_parents[vertex_plus]);
        vertex_plus = minus_parents[vertex_plus]->OtherNode(vertex_plus);
    }

    return path;
}

void SolverUnweighted::AugmentPath(std::vector<std::shared_ptr<Edge> > path) {
    int cur_vertex = UnmatchedVertex(path.front());

    for (int i = 0; i < static_cast<int>(path.size()); ++i) {
        int next_vertex = path[i]->OtherNode(cur_vertex);

        if (i % 2 == 0) {
            matched_edge[cur_vertex] = path[i];
            matched_edge[next_vertex] = path[i];
        }
        path[i]->matched = !path[i]->matched;
        cur_vertex = next_vertex;
    }
}

void SolverUnweighted::MakeCherryBlossom(const std::shared_ptr<Edge> &edge_plus_plus) {
    int first_vertex, second_vertex = -1;
    std::tie(first_vertex, second_vertex) = edge_plus_plus->Vertices();

    if (cherry_blossoms.SameSet(first_vertex, second_vertex)) {
        return;
    }

    if (verbose) {
        PrintTreeData();
        std::cout << "making cherry blossom " << first_vertex << " " << second_vertex << std::endl;
    }

    int first_bound, second_bound;
    std::tie(first_bound, second_bound) = PathUpperBounds(first_vertex, second_vertex);
    if (verbose) {
        std::cout << "upper bounds: " << first_bound << " " << second_bound << std::endl;
    }

    UpdatePath(first_vertex, first_bound);
    UpdatePath(second_vertex, second_bound);

    if (first_vertex != first_bound) {
        // first_vertex is not a root
        minus_parents[first_vertex] = edge_plus_plus;
    }
    if (second_vertex != second_bound) {
        // second_vertex is not a root
        minus_parents[second_vertex] = edge_plus_plus;
    }
}

void SolverUnweighted::AddVertex(int new_minus, int parent, const std::shared_ptr<Edge> &edge) {
    int next_plus = matched_edge[new_minus]->OtherNode(new_minus);

    plus[new_minus] = false;

    plus[next_plus] = true;

    minus_parents[new_minus] = edge;
    minus_parents[next_plus] = nullptr;

    cherry_blossoms.Detach(new_minus);
    cherry_blossoms.Detach(next_plus);

    _root_of_vertex[new_minus] = RootOfVertex(parent);
    _root_of_vertex[next_plus] = RootOfVertex(parent);

    growable_vertices.push(next_plus);
}

int SolverUnweighted::PlusPlusLCA(int first_vertex, int second_vertex) {
    if (verbose) {
        PrintTreeData();
    }

    ++lca_count;
    first_vertex = cherry_blossoms.Label(first_vertex);
    second_vertex = cherry_blossoms.Label(second_vertex);

    lca_markers[first_vertex] = lca_count;
    lca_markers[second_vertex] = lca_count;

    while (first_vertex != second_vertex) {
        if (matched_edge[first_vertex]) {
            first_vertex = matched_edge[first_vertex]->OtherNode(first_vertex);
            first_vertex = minus_parents[first_vertex]->OtherNode(first_vertex);
            first_vertex = cherry_blossoms.Label(first_vertex);
            if (lca_markers[first_vertex] == lca_count) {
                return first_vertex;
            }
            lca_markers[first_vertex] = lca_count;
        }

        if (matched_edge[second_vertex]) {
            if (verbose) {
                PrintTreeData();
                std::cout << "in lca: " << first_vertex << " " << second_vertex << std::endl;
            }
            second_vertex = matched_edge[second_vertex]->OtherNode(second_vertex);
            second_vertex = minus_parents[second_vertex]->OtherNode(second_vertex);
            if (verbose) {
                std::cout << first_vertex << " " << second_vertex << std::endl;
            }
            second_vertex = cherry_blossoms.Label(second_vertex);
            if (lca_markers[second_vertex] == lca_count) {
                return second_vertex;
            }
            lca_markers[second_vertex] = lca_count;
        }
    }

    return first_vertex;
}

std::pair<int, int> SolverUnweighted::PathUpperBounds(int first_vertex, int second_vertex) {
    int lca = PlusPlusLCA(first_vertex, second_vertex);
    int lca_receptacle = cherry_blossoms.Label(lca);

    if (verbose) {
        std::cout << "lca: " << lca << std::endl;
    }

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
        if (verbose) {
            std::cout << lower_vertex << " " << upper_vertex << std::endl;
        }

        int parent = matched_edge[lower_vertex]->OtherNode(lower_vertex);
        int grandparent = minus_parents[parent]->OtherNode(parent);

        cherry_blossoms.Unite(lower_vertex, upper_vertex, receptacle);
        cherry_blossoms.Unite(parent, upper_vertex, receptacle);

        if (!plus[parent]) {
            plus[parent] = true;
            growable_vertices.push(parent);
        }

        if (grandparent != upper_vertex) {
            // minus[grandparent] = true;
            minus_parents[grandparent] = minus_parents[parent];
        }

        lower_vertex = grandparent;
    }
}

void SolverUnweighted::ClearTree(int root, int other_root) {
    _root_of_vertex[root] = -1;

    for (int vertex : tree_plus_neighbors[root]) {
        if (RootOfVertex(vertex) != other_root) {
            growable_vertices.push(vertex);
        }
    }
}

int SolverUnweighted::UnmatchedVertex(const std::shared_ptr<Edge> &edge) const {
    auto [first_vertex, second_vertex] = edge->Vertices();
    if (!matched_edge[first_vertex]) {
        return first_vertex;
    }
    return second_vertex;
}

int SolverUnweighted::RootOfVertex(const int vertex) {
    int root = _root_of_vertex[vertex];

    if (root == -1) {
        return -1;
    }
    if (matched_edge[root]) {
        _root_of_vertex[vertex] = -1;
        return -1;
    }

    return root;
}
