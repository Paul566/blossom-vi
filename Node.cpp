#include "Node.h"
#include "EdgeWeighted.h"


Node::Node() {
    dual_variable = 0;
    parent_blossom = nullptr;
    children_blossom = std::vector<std::shared_ptr<Node>>();
}