#include <stdexcept>
#include "LabeledDisjointSets.h"

LabeledDisjointSets::LabeledDisjointSets(int size) {
    ranks = std::vector<int>(size, 0);
    parents.reserve(size);
    labels.reserve(size);
    for (int i = 0; i < size; ++i) {
        parents.push_back(i);
        labels.push_back(i);
    }
}

int LabeledDisjointSets::Representative(int element) {
    if ((element < 0) || (element >= static_cast<int>(ranks.size()))) {
        throw std::runtime_error("element not in disjoint sets");
    }

    int grandparent = element;
    while (parents[grandparent] != grandparent) {
        grandparent = parents[grandparent];
    }

    int current_element = element;
    while (parents[current_element] != grandparent) {
        int next_element = parents[current_element];
        parents[current_element] = grandparent;
        current_element = next_element;
    }

    return grandparent;
}

int LabeledDisjointSets::Label(int element) {
    return labels[Representative(element)];
}

bool LabeledDisjointSets::SameSet(int first, int second) {
    return Representative(first) == Representative(second);
}

void LabeledDisjointSets::UniteRepresentatives(int first, int second, int new_label) {
    if (first == second) {
        labels[first] = new_label;
        return;
    }

    if (ranks[first] < ranks[second]) {
        parents[first] = second;
        labels[second] = new_label;
        return;
    }

    if (ranks[first] > ranks[second]) {
        parents[second] = first;
        labels[first] = new_label;
        return;
    }

    parents[first] = second;
    labels[second] = new_label;
    ++ranks[second];
}

void LabeledDisjointSets::Unite(int first, int second, int new_label) {
    UniteRepresentatives(Representative(first), Representative(second), new_label);
}

void LabeledDisjointSets::Detach(const int element) {
    // TODO this should not be here later
    parents[element] = element;
    labels[element] = element;
    ranks[element] = 0;
}
