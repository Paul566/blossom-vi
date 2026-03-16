#ifndef BLOSSOM_VI_DUALUPDATER_H
#define BLOSSOM_VI_DUALUPDATER_H
#include <iostream>
#include <limits>
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

    UpdateType update_type = UpdateType::ConnectedComponents;
    int repetitions = 1;
};

struct DualConstraintsNode {
    int upper_bound = INT32_MAX;
    boost::container::small_vector<int, 4>  plus_plus_neighbors;
    boost::container::small_vector<int, 4>  plus_plus_constraints;
    boost::container::small_vector<int, 4>  plus_minus_neighbors;
    boost::container::small_vector<int, 4>  plus_minus_constraints;
};

struct Edge {
    int to;
    int rev;
    int cap;
    long long cost;
};

class MinCostFlow {
public:
    explicit MinCostFlow(int n)
        : n_(n),
          graph_(n),
          dist_(n),
          pot_(n, 0),
          prev_node_(n),
          prev_edge_(n) {}

    void AddEdge(int from, int to, int cap, long long cost) {
        graph_[from].push_back({to, static_cast<int>(graph_[to].size()), cap, cost});
        graph_[to].push_back({from, static_cast<int>(graph_[from].size()) - 1, 0, -cost});
    }

    std::pair<int, long long> MinCostMaxFlow(int source, int sink) {
        int flow = 0;
        long long cost = 0;

        const long long INF = std::numeric_limits<long long>::max();

        int num_paths = 0;
        while (true) {
            std::fill(dist_.begin(), dist_.end(), INF);
            dist_[source] = 0;

            std::priority_queue<
                std::pair<long long,int>,
                std::vector<std::pair<long long,int>>,
                std::greater<std::pair<long long,int>>
            > pq;

            pq.push({0, source});

            while (!pq.empty()) {
                auto current = pq.top();
                pq.pop();

                long long d = current.first;
                int v = current.second;

                if (d != dist_[v]) continue;

                for (int i = 0; i < static_cast<int>(graph_[v].size()); i++) {
                    Edge &e = graph_[v][i];
                    if (e.cap == 0) continue;

                    long long nd = d + e.cost + pot_[v] - pot_[e.to];
                    if (nd < 0) {
                        throw std::runtime_error("mincost flow loop");
                    }

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
            ++num_paths;

            for (int i = 0; i < n_; i++) {
                if (dist_[i] < INF) {
                    pot_[i] += dist_[i];
                }
            }

            int aug = std::numeric_limits<int>::max();

            for (int v = sink; v != source; v = prev_node_[v]) {
                Edge &e = graph_[prev_node_[v]][prev_edge_[v]];
                aug = std::min(aug, e.cap);
            }

            for (int v = sink; v != source; v = prev_node_[v]) {
                Edge &e = graph_[prev_node_[v]][prev_edge_[v]];
                e.cap -= aug;
                graph_[v][e.rev].cap += aug;
                cost += static_cast<long long>(aug) * e.cost;
            }

            flow += aug;
            if (flow > n_) {
                throw std::runtime_error("too much flow");
            }
        }

        // if (num_paths != (n_ - 2) / 2) {
        //     std::cout << "num_paths: " << num_paths << " " << n_ << std::endl;
        //     throw std::runtime_error("Flow is not maximum");
        // }

        return {flow, cost};
    }

    const std::vector<long long>& GetPotentials() const {
        return pot_;
    }

private:
    int n_;
    std::vector<std::vector<Edge>> graph_;
    std::vector<long long> dist_;
    std::vector<long long> pot_;
    std::vector<int> prev_node_;
    std::vector<int> prev_edge_;
};

class DualUpdater {
    public:
        const Parameters params;

        explicit DualUpdater(std::vector<DualConstraintsNode>&& constraints, const Parameters &params_ = {});
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
