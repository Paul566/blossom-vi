#ifndef BLOSSOM_VI_HEAP_H
#define BLOSSOM_VI_HEAP_H
#include <functional>
#include <iostream>
#include <stdexcept>

template<typename T>
class Heap {
    public:
        struct Handle {
            int index = 0;
        };

        explicit Heap(int d) : d_(d) {
            if (d_ < 2) {
                throw std::invalid_argument("d-arity must be >= 2");
            }
        }

        bool Empty() {
            return heap_.empty();
        }

        Handle *Push(T value, int key=0) { // TODO remove the default value for key
            auto handle = std::make_unique<Handle>();
            handle->index = heap_.size();

            Handle *handle_ptr = handle.get();
            // TODO consider using raw pointers, reallocation of a vector of unique_ptr is a bit more expensive
            handles_.emplace_back(std::move(handle));

            heap_.emplace_back(value, key, handle_ptr);
            SiftUp(heap_.size() - 1);

            return handle_ptr;
        }

        const T &Top() {
            if (heap_.empty()) {
                throw std::runtime_error("heap is empty");
            }
            return heap_[0].value;
        }

        void Pop() {
            if (!heap_.empty()) {
                EraseAt(0);
            }
        }

        std::vector<T> ElementsEqualToTop() {
            if (Empty()) {
                return {};
            }

            std::vector<int> indices;
            indices.push_back(0);
            int left = 0;
            int top_key = heap_[0].key;
            while (left < indices.size()) {
                int index = indices[left];
                ++left;

                for (int k = 0; k < d_; ++k) {
                    int c = Child(index, k);
                    if (c < heap_.size()) {
                        if (heap_[c].key == top_key) {
                            indices.push_back(c);
                        }
                    }
                }
            }

            std::vector<T> result;
            result.reserve(indices.size());
            for (int i : indices) {
                result.push_back(heap_[i].value);
            }
            return result;
        }

        bool Erase(Handle *h) {
            if (!h) {
                return false;
            }
            EraseAt(h->index);
            return true;
        }

        // debugging purposes
        void ValidateHeap(const std::string &msg = "") {
            for (std::size_t i = 1; i < heap_.size(); ++i) {
                if (heap_[i].key < heap_[Parent(i)].key ) {
                    std::cout << msg << std::endl;
                    std::cout << i << " " << Parent(i) << std::endl;
                    throw std::runtime_error("Incorrect heap");
                }
            }
        }

    // private:
        struct Node {
            T value;
            int key;
            Handle *handle;
        };

        int d_;
        std::vector<Node> heap_;
        std::vector<std::unique_ptr<Handle> > handles_;

        int Parent(int i) const {
            return (i - 1) / d_;
        }

        int Child(int i, int k) const {
            return d_ * i + k + 1;
        }

        void SwapNodes(int a, int b) {
            std::swap(heap_[a], heap_[b]);
            heap_[a].handle->index = a;
            heap_[b].handle->index = b;
        }

        void SiftUp(int i) {
            while (i > 0) {
                int p = Parent(i);
                if (heap_[i].key >= heap_[p].key) {
                    return;
                }
                SwapNodes(i, p);
                i = p;
            }
        }

        void SiftDown(int i) {
            while (true) {
                int best = i;

                for (int k = 0; k < d_; ++k) {
                    int c = Child(i, k);
                    if (c < heap_.size()) {
                        if (heap_[c].key < heap_[best].key) {
                            best = c;
                        }
                    }
                }

                if (best == i) {
                    return;
                }

                SwapNodes(i, best);
                i = best;
            }
        }

        void EraseAt(int i) {
            int last = heap_.size() - 1;
            if (i != last) {
                heap_[i] = std::move(heap_[last]);
                heap_[i].handle->index = i;
            }

            heap_.pop_back();

            if (i < heap_.size()) {
                SiftDown(i);
                SiftUp(i);
            }
        }

};
#endif //BLOSSOM_VI_HEAP_H
