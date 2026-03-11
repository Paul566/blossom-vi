#include "DualUpdater.h"

#include <iostream>
#include <numeric>
#include <queue>

DualUpdater::DualUpdater(std::vector<DualConstraintsNode> &&constraints,
                         const Parameters &params_) : params(params_), num_nodes(constraints.size()) {
    
    this->constraints = std::move(constraints);
    deltas = std::vector<int> (num_nodes, 0);
}

void DualUpdater::FindDeltas() {
    if (params.update_type == Parameters::UpdateType::ConnectedComponents) {
        for (int i = 0; i < params.repetitions; ++i) {
            FindDeltasCC();
            // std::cout << std::accumulate(deltas.begin(), deltas.end(), 0) << "\t";
        }
        // std::cout << std::endl;
    }
}

const std::vector<int> & DualUpdater::Deltas() {
    return deltas;
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
