#ifndef BLOSSOM_VI_DUALUPDATER_H
#define BLOSSOM_VI_DUALUPDATER_H
#include <iostream>
#include <queue>
#include <vector>
#include <boost/container/small_vector.hpp>

struct Parameters {
    enum class UpdateType {
        ConnectedComponents,
        StronglyConnectedComponents,
        ShortestPaths,
        LP
    };

    UpdateType update_type = UpdateType::LP;
    int repetitions = 1;
    int LP_threshold = 100; // use CC if more than this many variables
};

struct DualConstraintsNode {
    int upper_bound = INT32_MAX;
    boost::container::small_vector<int, 4> plus_plus_neighbors;
    boost::container::small_vector<int, 4> plus_plus_constraints;
    boost::container::small_vector<int, 4> plus_minus_neighbors;
    boost::container::small_vector<int, 4> plus_minus_constraints;
};

struct Edge {
    int to;
    int rev;
    int cap;
    int cost;
};

class MinCostFlow {
    public:
        explicit MinCostFlow(int n)
            : n_(n),
              graph_(n),
              dist_(n),
              pot_(n, 0),
              prev_node_(n),
              prev_edge_(n) {
        }

        void AddEdge(int from, int to, int cap, int cost) {
            graph_[from].push_back({to, static_cast<int>(graph_[to].size()), cap, cost});
            graph_[to].push_back({from, static_cast<int>(graph_[from].size()) - 1, 0, -cost});
        }

        void MinCostMaxFlow(int source, int sink) {
            const int64_t INF = INT64_MAX;

            while (true) {
                std::fill(dist_.begin(), dist_.end(), INF);
                dist_[source] = 0;

                std::priority_queue<
                    std::pair<int64_t, int>,
                    std::vector<std::pair<int64_t, int> >,
                    std::greater<std::pair<int64_t, int> >
                > pq;

                pq.push({0, source});

                while (!pq.empty()) {
                    auto current = pq.top();
                    pq.pop();

                    int64_t d = current.first;
                    int v = current.second;

                    if (d != dist_[v]) {
                        continue;
                    }

                    for (int i = 0; i < static_cast<int>(graph_[v].size()); ++i) {
                        Edge &e = graph_[v][i];
                        if (e.cap == 0) {
                            continue;
                        }

                        int64_t nd = d + e.cost + pot_[v] - pot_[e.to];

                        if (nd < dist_[e.to]) {
                            dist_[e.to] = nd;
                            prev_node_[e.to] = v;
                            prev_edge_[e.to] = i;
                            pq.push({nd, e.to});
                        }
                    }
                }

                if (dist_[sink] == INF) {
                    break;
                }

                for (int i = 0; i < n_; ++i) {
                    if (dist_[i] < INF) {
                        pot_[i] += dist_[i];
                    }
                }

                for (int v = sink; v != source; v = prev_node_[v]) {
                    Edge &e = graph_[prev_node_[v]][prev_edge_[v]];
                    e.cap -= 1; // I know that max capacity is 1 on the path
                    graph_[v][e.rev].cap += 1;
                }
            }
        }

        const std::vector<int64_t> &GetPotentials() const {
            return pot_;
        }

    private:
        int n_;
        std::vector<std::vector<Edge> > graph_;
        std::vector<int64_t> dist_;
        std::vector<int64_t> pot_;
        std::vector<int> prev_node_;
        std::vector<int> prev_edge_;
};

class DualUpdater {
    public:
        const Parameters params;

        explicit DualUpdater(std::vector<DualConstraintsNode> &&constraints, const Parameters &params_ = {});
        void FindDeltas();
        const std::vector<int> &Deltas();

    private:
        int num_nodes;
        std::vector<DualConstraintsNode> constraints;
        std::vector<int> deltas;

        void ValidateConstraintNonNegativity();
        void ValidateDeltasFeasibility();

        void FindDeltasCC();
        std::vector<std::vector<int> > ConnectedComponents();

        void FindDeltasShortestPaths();
        void SetDeltasToDists();
        void RefineUpperBounds();

        void FindDeltasMinCostFlow();
};

#endif //BLOSSOM_VI_DUALUPDATER_H
