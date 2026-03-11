#ifndef BLOSSOM_VI_DUALUPDATER_H
#define BLOSSOM_VI_DUALUPDATER_H
#include <vector>
#include <boost/container/small_vector.hpp>

struct Parameters {
    enum class UpdateType {
        ConnectedComponents,
        StronglyConnectedComponents
    };

    UpdateType update_type = UpdateType::ConnectedComponents;
    int repetitions = 3;
};

struct DualConstraintsNode {
    int upper_bound = INT32_MAX;
    boost::container::small_vector<int, 4>  plus_plus_neighbors;
    boost::container::small_vector<int, 4>  plus_plus_constraints;
    boost::container::small_vector<int, 4>  plus_minus_neighbors;
    boost::container::small_vector<int, 4>  plus_minus_constraints;
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

        void FindDeltasCC();
        std::vector<std::vector<int> > ConnectedComponents();
};

#endif //BLOSSOM_VI_DUALUPDATER_H
