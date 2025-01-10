#include "SolverUnweighted.h"

#include <algorithm>
#include <iostream>
#include <unordered_set>

SolverUnweighted::SolverUnweighted(const std::vector<std::vector<int> > &adj_list_,
                                   int greedy_init_type_,
                                   bool verbose_) : n(static_cast<int>(adj_list_.size())),
                                                            cherry_blossoms(LabeledDisjointSets(n)), verbose(verbose_),
                                                            greedy_init_type(greedy_init_type_) {
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
    root_of_vertex = std::vector<int>(n, -1);
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
        std::cout << "(receptacle: " << cherry_blossoms.Label(v) << ") ";
        std::cout << "(root: " << root_of_vertex[v] << ") ";
        std::cout << "(is a plus: " << plus[v] << ") ";
        std::cout << std::endl;
    }
}

void SolverUnweighted::Solve() {
    GreedyInit();
    if (verbose) {
        PrintMatching();
    }

    for (int root = 0; root < n; ++root) {
        if (!matched_edge[root]) {
            growable_vertices.push(root);
            plus[root] = true;;
            root_of_vertex[root] = root;
        }
    }

    while (!growable_vertices.empty()) {
        int cur_vertex = growable_vertices.front();
        if (verbose) {
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
        std::vector<std::pair<int, int>> vertices;
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
        std::vector<std::pair<int, std::shared_ptr<Edge>>> edges;
        for (int i = 0; i < n; ++i) {
            for (const auto& edge : adj_list[i]) {
                if (i < edge->OtherNode(i)) {
                    edges.emplace_back(adj_list[i].size() + adj_list[edge->OtherNode(i)].size(), edge);
                }
            }
        }
        std::sort(edges.begin(), edges.end());
        for (auto & edge : edges) {
            auto [head, tail] = edge.second->Vertices();
            if ((!matched_edge[head]) && (matched_edge[tail])) {
                edge.second->matched = true;
                matched_edge[head] = edge.second;
                matched_edge[tail] = edge.second;
            }
        }
        return;
    }

    if (greedy_init_type == 3) {
        std::vector<std::pair<int, std::shared_ptr<Edge>>> edges;
        for (int i = 0; i < n; ++i) {
            for (const auto& edge : adj_list[i]) {
                if (i < edge->OtherNode(i)) {
                    edges.emplace_back(std::max(adj_list[i].size(), adj_list[edge->OtherNode(i)].size()), edge);
                }
            }
        }
        std::sort(edges.begin(), edges.end());
        for (auto & edge : edges) {
            auto [head, tail] = edge.second->Vertices();
            if ((!matched_edge[head]) && (matched_edge[tail])) {
                edge.second->matched = true;
                matched_edge[head] = edge.second;
                matched_edge[tail] = edge.second;
            }
        }
    }
}

bool SolverUnweighted::HandleVertex(const int cur_vertex) {
    if (!plus[cur_vertex]) {
        // the vertex is in the queue but its tree was deleted earlier
        return false;
    }

    for (auto edge : adj_list[cur_vertex]) {
        int to = edge->OtherNode(cur_vertex);

        if ((!plus[to]) && (!minus[to])) {
            // to is not in the tree yet
            int next_plus = matched_edge[to]->OtherNode(to);
            minus[to] = true;
            plus[next_plus] = true;
            minus_parents[to] = edge;
            children[cur_vertex].push_back(to);
            children[to].push_back(next_plus);
            root_of_vertex[to] = root_of_vertex[cur_vertex];
            root_of_vertex[next_plus] = root_of_vertex[cur_vertex];
            growable_vertices.push(next_plus);
        } else {
            // to is in the tree
            if (plus[to]) {
                if (root_of_vertex[cur_vertex] == root_of_vertex[to]) {
                    MakeCherryBlossom(edge);
                } else {
                    Augment(edge, cur_vertex, to);
                    return true;
                }
            }
        }
    }

    return false;
}

void SolverUnweighted::Augment(std::shared_ptr<Edge> edge_plus_plus, int cur_vertex, int to) {
    auto [first_vertex, second_vertex] = edge_plus_plus->Vertices();
    auto first_path = PathToRoot(first_vertex);
    auto second_path = PathToRoot(second_vertex);

    std::vector<std::shared_ptr<Edge> > path;
    path.reserve(first_path.size() + second_path.size() + 1);
    for (int i = static_cast<int>(first_path.size()) - 1; i >= 0; --i) {
        path.push_back(first_path[i]);
    }
    path.push_back(edge_plus_plus);
    for (auto edge : second_path) {
        path.push_back(edge);
    }

    int first_root = UnmatchedVertex(path.front());
    int second_root = UnmatchedVertex(path.back());

    AugmentPath(path);
    ClearTree(first_root);
    ClearTree(second_root);
}

std::vector<std::shared_ptr<Edge> > SolverUnweighted::PathToRoot(int vertex_plus) {
    if (!matched_edge[vertex_plus]) {
        // already in a root
        return std::vector<std::shared_ptr<Edge> >();
    }

    if (!plus[vertex_plus]) {
        throw std::runtime_error("In PathToRoot: vertex is not a plus");
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
    if (path.size() % 2 == 0) {
        throw std::runtime_error("in Augment: the length of the path must be odd");
    }

    for (int i = 0; i < static_cast<int>(path.size()); ++i) {
        if ((i % 2 == 0) && (path[i]->matched)) {
            throw std::runtime_error("in Augment: not an alternating path");
        }
        if ((i % 2 == 1) && (!path[i]->matched)) {
            throw std::runtime_error("in Augment: not an alternating path");
        }
    }

    int cur_vertex = UnmatchedVertex(path.front());

    if (matched_edge[cur_vertex]) {
        throw std::runtime_error("in Augment: start must be unmatched");
    }

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

void SolverUnweighted::MakeCherryBlossom(std::shared_ptr<Edge> edge_plus_plus) {
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
        root_of_vertex[vertex] = -1;

        // have to make adjacent plus vertices from other trees growable again
        for (auto edge : adj_list[vertex]) {
            int to = edge->OtherNode(vertex);
            if ((plus[to]) && (root_of_vertex[to] != root)) {
                growable_vertices.push(to);
            }
        }
    }
}

int SolverUnweighted::UnmatchedVertex(const std::shared_ptr<Edge> &edge) {
    auto [first_vertex, second_vertex] = edge->Vertices();
    if (!matched_edge[first_vertex]) {
        return first_vertex;
    }
    if (matched_edge[second_vertex]) {
        throw std::runtime_error("In UnmatchedVertex: no unmatched vertex");
    }
    return second_vertex;
}
