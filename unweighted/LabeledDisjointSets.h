#ifndef BLOSSOM_VI_LABELEDDISJOINTSETS_H
#define BLOSSOM_VI_LABELEDDISJOINTSETS_H


#include <vector>

class LabeledDisjointSets {
public:
    explicit LabeledDisjointSets(int size);

    int Label(int element);

    bool SameSet(int first, int second);

    void Unite(int first, int second, int new_label);

    void Detach(int element);

private:
    std::vector<int> ranks;
    std::vector<int> parents;
    std::vector<int> labels;

    int Representative(int element);

    void UniteRepresentatives(int first, int second, int new_label);
};


#endif
