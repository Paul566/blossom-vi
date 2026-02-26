#ifndef BLOSSOM_VI_PAIRINGHEAP_H
#define BLOSSOM_VI_PAIRINGHEAP_H

#include <utility>

template<class Obj>
class PairingHeap {
    public:
        Obj *GetMin() const {
            return root;
        }

        void Insert(Obj *node) {
            node->heap_child = nullptr;
            node->heap_next = nullptr;
            node->heap_prev = nullptr;

            root = Meld(root, node);
        }

        void RemoveMin() {
            if (!root)
                return;

            Obj *old_root = root;
            Obj *children = root->heap_child;

            if (children) {
                children->heap_prev = nullptr; // new root candidates
            }

            root = TwoPassMerge(children);

            old_root->heap_child = nullptr;
        }

        void Remove(Obj *node) {
            if (node == root) {
                RemoveMin();
                return;
            }

            Cut(node);

            Obj *subtree = TwoPassMerge(node->heap_child);
            if (subtree) {
                subtree->heap_prev = nullptr;
            }

            root = Meld(root, subtree);

            node->heap_child = nullptr;
        }

        std::vector<Obj *> ElementsEqualToTop() {
            Obj *cur_root = root;
            RemoveMin();
            Insert(cur_root);

            std::vector<Obj *> result;

            if (!root) {
                return result;
            }

            int min_key = root->Key();

            std::vector<Obj *> stack;
            stack.push_back(root);

            while (!stack.empty()) {
                Obj *node = stack.back();
                stack.pop_back();

                result.push_back(node);

                for (Obj *child = node->heap_child; child; child = child->heap_next) {
                    if (child->Key() == min_key) {
                        stack.push_back(child);
                    }
                }
            }

            return result;
        }

    private:
        Obj *root = nullptr;

        static Obj *Meld(Obj *a, Obj *b) {
            if (!a) return b;
            if (!b) return a;

            if (b->Key() < a->Key()) {
                std::swap(a, b);
            }

            b->heap_prev = a;
            b->heap_next = a->heap_child;
            if (a->heap_child) {
                a->heap_child->heap_prev = b;
            }

            a->heap_child = b;

            return a;
        }

        static Obj *TwoPassMerge(Obj *first) {
            if (!first || !first->heap_next) {
                return first;
            }

            Obj *a = first;
            Obj *b = a->heap_next;
            Obj *rest = b->heap_next;

            a->heap_next = nullptr;
            b->heap_next = nullptr;

            Obj *merged = Meld(a, b);
            // TODO avoid recursion
            Obj *remaining = TwoPassMerge(rest);

            return Meld(merged, remaining);
        }

        void Cut(Obj *node) {
            if (!node->heap_prev) {
                return; // already root
            }

            if (node->heap_prev->heap_child == node) {
                // node is leftmost child
                Obj *parent = node->heap_prev;
                parent->heap_child = node->heap_next;
                if (node->heap_next) {
                    node->heap_next->heap_prev = parent;
                }
            } else {
                // node has left sibling
                Obj *leftSibling = node->heap_prev;
                leftSibling->heap_next = node->heap_next;
                if (node->heap_next) {
                    node->heap_next->heap_prev = leftSibling;
                }
            }

            node->heap_prev = nullptr;
            node->heap_next = nullptr;
        }
};

#endif //BLOSSOM_VI_PAIRINGHEAP_H
