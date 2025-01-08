#ifndef BLOSSOM_VI_LABELEDDISJOINTSETS_H
#define BLOSSOM_VI_LABELEDDISJOINTSETS_H


#include <vector>

class LabeledDisjointSets {
public:
    explicit LabeledDisjointSets(int size);

    int Representative(int element);

    int Label(int element);

    void Unite(int first, int second, int new_label);

    void Detach(int element);

private:
    std::vector<int> ranks;
    std::vector<int> parents;
    std::vector<int> labels;

    void UniteRepresentatives(int first, int second, int new_label);
};


#endif
