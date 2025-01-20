#ifndef BLOSSOM_VI_SOLVERUNWEIGHTED_H
#define BLOSSOM_VI_SOLVERUNWEIGHTED_H

#include <vector>
#include <memory>
#include <queue>
#include "Edge.h"
#include "LabeledDisjointSets.h"

class SolverUnweighted {
    public:
        std::vector<std::vector<std::shared_ptr<Edge> > > adj_list;

        explicit SolverUnweighted(const std::vector<std::vector<int> > &adj_list_,
                                  int greedy_init_type_, bool delete_edges_in_cherries_, bool verbose_ = false);

        std::vector<std::pair<int, int>> Matching();

        void PrintMatching();

        void PrintAdjList();

        void PrintTreeData();

        void PrintTreeStats();

        void Solve();

    private:
        const int n; // the number of vertices
        std::vector<std::shared_ptr<Edge> > matched_edge;
        std::vector<std::shared_ptr<Edge> > minus_parents;
        std::vector<std::vector<int> > children;
        std::vector<bool> plus;
        std::vector<bool> minus;
        std::queue<int> growable_vertices;
        LabeledDisjointSets cherry_blossoms; // labels are receptacles
        std::vector<int> _root_of_vertex;
        std::vector<int> lca_markers;
        int lca_count;
        std::vector<std::vector<int>> tree_plus_neighbors; // roots remember plus vertices adjacent to their trees
        bool verbose;
        int greedy_init_type;
        // 0: usual greedy init
        // 1: go through vertices from the lowest degree to the highest degree
        // 2: go through edges (u, v) sorted by deg(u)+deg(v)
        // 3: go through edges (u, v) sorted by max(deg(u), deg(v))
        bool delete_edges_in_cherries;

        void GreedyInit();

        bool HandleVertex(int cur_vertex);
        // returns true if performed an augmentation

        void Augment(const std::shared_ptr<Edge>& edge_plus_plus, int cur_vertex, int to);

        std::vector<std::shared_ptr<Edge> > PathToRoot(int vertex_plus);
        // a vector of edges that lead to the root, starting from a plus vertex

        void AugmentPath(std::vector<std::shared_ptr<Edge> > path);

        void MakeCherryBlossom(const std::shared_ptr<Edge>& edge_plus_plus);
        // first_vertex and second_vertex are pluses
        // walks up the tree until the lca and makes everyone plusminus,
        // adding the new pluses to the queue and updating parent pointers

        void AddVertex(int new_minus, int parent, std::shared_ptr<Edge> &edge);

        int PlusPlusLCA(int first_vertex, int second_vertex);

        std::pair<int, int> PathUpperBounds(int first_vertex, int second_vertex);

        void UpdatePath(int lower_vertex, int upper_vertex);

        void ClearTree(int root, int other_root);

        int UnmatchedVertex(const std::shared_ptr<Edge> &edge);

        int RootOfVertex(const int vertex);
};

#endif
