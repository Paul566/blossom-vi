#include "DualUpdater.h"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <queue>

DualUpdater::DualUpdater(std::vector<DualConstraintsNode> &&constraints,
                         const Parameters &params_) : params(params_), num_nodes(constraints.size()) {
    this->constraints = std::move(constraints);
    deltas = std::vector<int>(num_nodes, 0);

    // ValidateConstraintNonNegativity();
}

void DualUpdater::FindDeltas() {
    if (params.update_type == Parameters::UpdateType::ConnectedComponents) {
        for (int i = 0; i < params.repetitions; ++i) {
            FindDeltasCC();
            // std::cout << std::accumulate(deltas.begin(), deltas.end(), 0) << "\t";
        }
        // std::cout << std::endl;
        return;
    }

    if (params.update_type == Parameters::UpdateType::ShortestPaths) {
        FindDeltasCC();
        std::cout << std::accumulate(deltas.begin(), deltas.end(), 0) << "\t";
        for (int u = 0; u < num_nodes; ++u) {
            deltas[u] = 0;
        }

        FindDeltasShortestPaths();
        ValidateDeltasFeasibility();
        std::cout << std::accumulate(deltas.begin(), deltas.end(), 0) << std::endl;
    }

    if (params.update_type == Parameters::UpdateType::LP) {
        if (num_nodes > 30) {
            FindDeltasCC();
        } else {
            FindDeltasMinCostFlow();
            ValidateDeltasFeasibility();
            // std::cout << std::accumulate(deltas.begin(), deltas.end(), 0);
            FindDeltasCC();
            // std::cout << " " << std::accumulate(deltas.begin(), deltas.end(), 0) << std::endl;
            ValidateDeltasFeasibility();
        }
    }
}

const std::vector<int> &DualUpdater::Deltas() {
    return deltas;
}

void DualUpdater::ValidateConstraintNonNegativity() {
    for (const DualConstraintsNode &node : constraints) {
        if (node.upper_bound < 0) {
            throw std::runtime_error("ValidateConstraintNonNegativity: upper_bound");
        }
        for (int c : node.plus_plus_constraints) {
            if (c < 0) {
                throw std::runtime_error("ValidateConstraintNonNegativity: plus_plus_constraints");
            }
        }
        for (int c : node.plus_minus_constraints) {
            if (c < 0) {
                throw std::runtime_error("ValidateConstraintNonNegativity: plus_minus_constraints");
            }
        }
    }
}

void DualUpdater::ValidateDeltasFeasibility() {
    for (int u = 0; u < num_nodes; ++u) {
        if (deltas[u] < 0) {
            throw std::runtime_error("ValidateDeltasFeasibility: negative delta");
        }

        if (deltas[u] > constraints[u].upper_bound) {
            throw std::runtime_error("ValidateDeltasFeasibility: upper_bound");
        }

        for (int i = 0; i < constraints[u].plus_plus_neighbors.size(); ++i) {
            int v = constraints[u].plus_plus_neighbors[i];
            if (deltas[u] + deltas[v] > constraints[u].plus_plus_constraints[i]) {
                throw std::runtime_error("ValidateDeltasFeasibility: plus plus");
            }
        }

        for (int i = 0; i < constraints[u].plus_minus_neighbors.size(); ++i) {
            int v = constraints[u].plus_minus_neighbors[i];
            if (deltas[u] - deltas[v] > constraints[u].plus_minus_constraints[i]) {
                throw std::runtime_error("ValidateDeltasFeasibility: plus minus");
            }
        }
    }
}

void DualUpdater::FindDeltasCC() {
    std::vector<std::vector<int> > connected_components = ConnectedComponents();

    std::vector<int> delta_increase(num_nodes, 0);

    std::vector<int> component_index(num_nodes, 0);
    for (int i = 0; i < static_cast<int>(connected_components.size()); ++i) {
        for (int v : connected_components[i]) {
            component_index[v] = i;
        }
    }

    for (int cc_index = 0; cc_index < static_cast<int>(connected_components.size()); ++cc_index) {
        int delta = INT32_MAX;

        // find delta for this connected component
        for (int v : connected_components[cc_index]) {
            if (constraints[v].upper_bound - deltas[v] < delta) {
                delta = constraints[v].upper_bound - deltas[v];
            }

            for (int i = 0; i < constraints[v].plus_plus_constraints.size(); ++i) {
                int w = constraints[v].plus_plus_neighbors[i];
                int slack = constraints[v].plus_plus_constraints[i] - deltas[v] - deltas[w];

                if (component_index[w] == cc_index) {
                    if (slack % 2 != 0) {
                        throw std::runtime_error("In VariableDeltas: slack is not divisible by 2");
                    }
                    if (slack / 2 < delta) {
                        delta = slack / 2;
                    }
                } else {
                    if (slack - delta_increase[w] < delta) {
                        delta = slack - delta_increase[w];
                    }
                }
            }

            for (int i = 0; i < constraints[v].plus_minus_constraints.size(); ++i) {
                int w = constraints[v].plus_minus_neighbors[i];
                int slack = constraints[v].plus_minus_constraints[i] - deltas[v] + deltas[w];
                if (component_index[w] == cc_index) {
                    continue;
                }
                if (slack + delta_increase[w] < delta) {
                    delta = slack + delta_increase[w];
                }
            }
        }

        // apply delta to this connected component
        for (int v : connected_components[cc_index]) {
            delta_increase[v] = delta;
        }
    }

    for (int v = 0; v < num_nodes; ++v) {
        deltas[v] += delta_increase[v];
        if (delta_increase[v] < 0) {
            throw std::runtime_error("Negative delta increase");
        }
    }

    for (int d : deltas) {
        if (d < 0) {
            throw std::runtime_error("In VariableDeltas: found negative delta");
        }
    }
}

std::vector<std::vector<int> > DualUpdater::ConnectedComponents() {
    std::vector<bool> visited(num_nodes, false);
    std::vector<std::vector<int> > components;

    // Build undirected adjacency list for weak connectivity
    std::vector<std::vector<int> > adj_list_tree_tree(num_nodes, std::vector<int>());
    for (int u = 0; u < static_cast<int>(num_nodes); ++u) {
        for (int i = 0; i < constraints[u].plus_minus_constraints.size(); ++i) {
            int v = constraints[u].plus_minus_neighbors[i];
            int slack = constraints[u].plus_minus_constraints[i] - deltas[u] + deltas[v];
            if (slack == 0) {
                adj_list_tree_tree[u].emplace_back(v);
                adj_list_tree_tree[v].emplace_back(u);
            }
        }
    }

    for (int start = 0; start < static_cast<int>(num_nodes); ++start) {
        if (!visited[start]) {
            std::vector<int> comp;
            std::queue<int> q;

            visited[start] = true;
            q.push(start);

            while (!q.empty()) {
                int u = q.front();
                q.pop();
                comp.push_back(u);

                for (int v : adj_list_tree_tree[u]) {
                    if (!visited[v]) {
                        visited[v] = true;
                        q.push(v);
                    }
                }
            }

            components.push_back(std::move(comp));
        }
    }

    return components;
}

void DualUpdater::FindDeltasShortestPaths() {
    RefineUpperBounds();
    SetDeltasToDists();

    // make plus-plus constraints respected
    for (int u = 0; u < num_nodes; ++u) {
        for (int i = 0; i < constraints[u].plus_plus_neighbors.size(); ++i) {
            int v = constraints[u].plus_plus_neighbors[i];
            int violation = deltas[u] + deltas[v] - constraints[u].plus_plus_constraints[i];
            if (violation > 0) {
                int half_violation = violation / 2;
                deltas[u] -= half_violation;
                deltas[v] -= (violation - half_violation);

                if (deltas[u] < 0 || deltas[v] < 0) {
                    throw std::runtime_error("In FindDeltasShortestPaths: negative deltas");
                }
            }
        }
    }

    // make the plus minus constraints respected once again
    std::vector<std::vector<std::pair<int, int>>> reversed_pm(num_nodes);
    for (int u = 0; u < num_nodes; ++u) {
        for (int i = 0; i < constraints[u].plus_minus_neighbors.size(); ++i) {
            reversed_pm[constraints[u].plus_minus_neighbors[i]].emplace_back(u, constraints[u].plus_minus_constraints[i]);
        }
    }
    std::priority_queue<
        std::pair<int, int>,
        std::vector<std::pair<int, int> >,
        std::greater<>
    > pq;
    for (int v = 0; v < num_nodes; ++v) {
        pq.push({deltas[v], v});
    }

    while (!pq.empty()) {
        auto [d,u] = pq.top();
        pq.pop();

        if (d > deltas[u]) {
            continue;
        }

        for (auto [v, w] : reversed_pm[u]) {
            if (deltas[v] > d + w) {
                deltas[v] = d + w;
                pq.push({deltas[v], v});
            }
        }
    }

    if (std::all_of(deltas.begin(), deltas.end(), [](auto x) { return x == 0; })) {
        FindDeltasCC();
    }
}

void DualUpdater::SetDeltasToDists() {
    std::priority_queue<
        std::pair<int, int>,
        std::vector<std::pair<int, int> >,
        std::greater<>
    > pq;

    std::vector<std::vector<std::pair<int, int>>> reversed_pm(num_nodes);
    for (int u = 0; u < num_nodes; ++u) {
        for (int i = 0; i < constraints[u].plus_minus_neighbors.size(); ++i) {
            reversed_pm[constraints[u].plus_minus_neighbors[i]].emplace_back(u, constraints[u].plus_minus_constraints[i]);
        }
    }


    for (int v = 0; v < num_nodes; ++v) {
        deltas[v] = constraints[v].upper_bound;
        pq.push({deltas[v], v});
    }

    while (!pq.empty()) {
        auto [d,u] = pq.top();
        pq.pop();

        if (d > deltas[u]) {
            continue;
        }

        for (auto [v, w] : reversed_pm[u]) {
            if (deltas[v] > d + w) {
                deltas[v] = d + w;
                pq.push({deltas[v], v});
            }
        }
    }
}

void DualUpdater::RefineUpperBounds() {
    for (int u = 0; u < num_nodes; ++u) {
        for (int i = 0; i < constraints[u].plus_plus_neighbors.size(); ++i) {
            int c = constraints[u].plus_plus_constraints[i];
            if (c < constraints[u].upper_bound) {
                constraints[u].upper_bound = c;
            }
        }
    }
}

void DualUpdater::FindDeltasMinCostFlow() {
    int source = 2 * num_nodes;
    int sink = 2 * num_nodes + 1;
    MinCostFlow mcf(2 * num_nodes + 2);
    // (delta pluses, delta minuses, source, sink)

    const int INF = 1000000000; // TODO int32max?

    for (int i = 0; i < num_nodes; ++i) {
        mcf.AddEdge(source, num_nodes + i, 1, 0);
        mcf.AddEdge(i, sink, 1, 0);

        if (constraints[i].upper_bound < INT32_MAX / 2) {
            mcf.AddEdge(num_nodes + i, i, INF, 2 * constraints[i].upper_bound);
        } else {
            mcf.AddEdge(num_nodes + i, i, INF, INT32_MAX);
        }
        mcf.AddEdge(i, num_nodes + i, INF, 0);
    }

    for (int u = 0; u < num_nodes; ++u) {
        for (int i = 0; i < constraints[u].plus_plus_neighbors.size(); ++i) {
            int c = constraints[u].plus_plus_constraints[i];
            int v = constraints[u].plus_plus_neighbors[i];
            if (u < v) {
                mcf.AddEdge(num_nodes + v, u, INF, c);
                mcf.AddEdge(num_nodes + u, v, INF, c);
            }
        }

        for (int i = 0; i < constraints[u].plus_minus_neighbors.size(); ++i) {
            int c = constraints[u].plus_minus_constraints[i];
            int v = constraints[u].plus_minus_neighbors[i];

            mcf.AddEdge(v, u, INF, c);
            mcf.AddEdge(num_nodes + u, num_nodes + v, INF, c);
        }
    }

    auto [flow, cost] = mcf.MinCostMaxFlow(source, sink);
    std::vector<long long> pot = mcf.GetPotentials();

    int source_pot = pot[source];
    // std::cout << cost / 2. << std::endl;
    // std::cout << pot[source] << " " << pot[sink] << " " << source_pot << std::endl;
    for (int u = 0; u < num_nodes; ++u) {
        deltas[u] = (pot[u] - pot[num_nodes + u]) / 2;
        // std::cout << pot[u] << " " << pot[num_nodes + u] << " " << constraints[u].upper_bound << std::endl;
    }
}
