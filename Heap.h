#ifndef BLOSSOM_VI_HEAP_H
#define BLOSSOM_VI_HEAP_H
#include <functional>
#include <iostream>
#include <stdexcept>
#include <queue>

template<typename T, typename Compare = std::less<T>>
class Heap {
    public:
        struct Handle {
            int index = 0;
            bool alive = false;
        };

        explicit Heap(int d, Compare comp = Compare{}) : d_(d), comp_(comp) {
            if (d_ < 2) {
                throw std::invalid_argument("d-arity must be >= 2");
            }
        }

        bool Empty() {
            return heap_.empty();
        }

        Handle *Push(T value) {
            auto handle = std::make_unique<Handle>();
            handle->index = heap_.size();
            handle->alive = true;

            Handle *handle_ptr = handle.get();
            handles_.push_back(std::move(handle));

            heap_.push_back({std::move(value), handle_ptr});
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

            // TODO make better
            std::vector<int> indices;
            std::queue<int> to_process;
            to_process.push(0);
            T top = heap_[0].value;
            while (!to_process.empty()) {
                int index = to_process.front();
                to_process.pop();
                indices.push_back(index);

                for (int k = 0; k < d_; ++k) {
                    int c = Child(index, k);
                    if (c < heap_.size()) {
                        if (!comp_(heap_[c].value, top) && !comp_(top, heap_[c].value)) {
                            // TODO make the equality test work in another way, maybe have int Key instead of Comparator
                            to_process.push(c);
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
            if (!h || !h->alive) {
                return false;
            }
            EraseAt(h->index);
            return true;
        }

        // debugging purposes
        void ValidateHeap(const std::string &msg = "") {
            for (std::size_t i = 1; i < heap_.size(); ++i) {
                if (comp_(heap_[i].value, heap_[Parent(i)].value )) {
                    std::cout << msg << std::endl;
                    std::cout << i << " " << Parent(i) << std::endl;
                    throw std::runtime_error("Incorrect heap");
                }
            }
        }

    // private:
        struct Node {
            T value;
            Handle *handle;
        };

        int d_;
        Compare comp_;
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
                if (!comp_(heap_[i].value, heap_[p].value)) {
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
                        if (comp_(heap_[c].value, heap_[best].value)) {
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
            Node &victim = heap_[i];
            victim.handle->alive = false;

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
