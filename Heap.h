#ifndef BLOSSOM_VI_HEAP_H
#define BLOSSOM_VI_HEAP_H
#include <functional>
#include <stdexcept>

template<
    typename T,
    typename Compare = std::less<T>,
    typename Validity = std::function<bool(const T &)> >
class Heap {
    public:
        struct Handle {
            int index = 0;
            bool alive = false;
        };

    private:
        struct Node {
            T value;
            Handle *handle;
        };

        std::vector<Node> heap_;
        std::vector<std::unique_ptr<Handle> > handles_;

        int d_;
        Compare comp_;
        Validity valid_;

    public:
        explicit Heap(
            int d,
            Compare comp = Compare{},
            Validity valid = [](const T &) { return true; }
        ) : d_(d), comp_(comp), valid_(std::move(valid)) {
            if (d_ < 2)
                throw std::invalid_argument("d-arity must be >= 2");
        }

        bool Empty() {
            CleanupInvalids();
            return heap_.empty();
        }

        int Size() {
            CleanupInvalids();
            return heap_.size();
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
            CleanupInvalids();
            if (heap_.empty()) {
                throw std::runtime_error("heap is empty");
            }
            return heap_[0].value;
        }

        void Pop() {
            CleanupInvalids();
            if (!heap_.empty()) {
                EraseAt(0);
            }
        }

        bool Erase(Handle *h) {
            if (!h || !h->alive) {
                return false;
            }
            EraseAt(h->index);
            return true;
        }

    private:
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
                if (!valid_(heap_[i].value)) {
                    EraseAt(i);
                    return;
                }

                int p = Parent(i);
                if (!valid_(heap_[p].value)) {
                    EraseAt(p);
                    continue; // recheck current i
                }

                if (!comp_(heap_[i].value, heap_[p].value)) {
                    break;
                }

                SwapNodes(i, p);
                i = p;
            }
        }

        void SiftDown(int i) {
            while (true) {
                if (!valid_(heap_[i].value)) {
                    EraseAt(i);
                    return;
                }

                int best = i;

                for (int k = 0; k < d_; ++k) {
                    int c = Child(i, k);
                    if (c < heap_.size()) {
                        if (!valid_(heap_[c].value)) {
                            EraseAt(c);
                            if (c < heap_.size()) {
                                k--; // recheck swapped element
                            }
                            continue;
                        }
                        if (comp_(heap_[c].value, heap_[best].value)) {
                            best = c;
                        }
                    }
                }

                if (best == i) {
                    break;
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

        void CleanupInvalids() {
            while (!heap_.empty() && !valid_(heap_[0].value)) {
                EraseAt(0);
            }
        }
};
#endif //BLOSSOM_VI_HEAP_H
