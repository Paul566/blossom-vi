#include "Node.h"
#include "EdgeWeighted.h"


Node::Node(int index_) : index(index_) {
    dual_variable_quadrupled = 0;
    parent_blossom = nullptr;
    matched_edge = nullptr;
    children_blossom = std::vector<std::shared_ptr<Node>>();
}
