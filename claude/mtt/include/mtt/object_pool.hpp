#pragma once
#include <stack>
#include <memory>
#include <mutex>
#include <atomic>
#include <cstddef>
#include <new>

namespace mtt {

// ─────────────────────────────────────────────────────────────────────────────
// Fixed-capacity object pool — zero heap allocation at runtime.
//
// Objects are allocated once via ::operator new[] into a raw byte buffer
// and constructed in-place with placement-new. This sidesteps the
// std::vector<T>::resize() requirement that T be move-constructible,
// allowing T to contain std::mutex or other immovable members.
//
// T must be: default-constructible, and expose a reset() method.
// ─────────────────────────────────────────────────────────────────────────────
template<typename T, std::size_t Capacity = 10'000>
class ObjectPool {
public:
    ObjectPool() {
        // Allocate raw storage (no construction yet)
        storage_ = static_cast<T*>(
            ::operator new[](Capacity * sizeof(T), std::align_val_t{alignof(T)}));

        // Construct all T objects in-place and seed the free stack
        for (std::size_t i = 0; i < Capacity; ++i) {
            new (&storage_[i]) T{};
            free_.push(&storage_[i]);
        }
    }

    ~ObjectPool() {
        // Destroy all objects in reverse order, then free raw storage
        for (std::size_t i = Capacity; i-- > 0; )
            storage_[i].~T();
        ::operator delete[](storage_,
                            std::align_val_t{alignof(T)});
    }

    // Non-copyable, non-movable
    ObjectPool(const ObjectPool&)            = delete;
    ObjectPool& operator=(const ObjectPool&) = delete;

    // Acquire an idle object. Returns nullptr if the pool is exhausted.
    T* acquire() {
        std::lock_guard<std::mutex> lk(mx_);
        if (free_.empty()) return nullptr;
        T* obj = free_.top();
        free_.pop();
        ++active_count_;
        return obj;
    }

    // Return an object to the pool. Calls obj->reset() to clear state.
    void release(T* obj) {
        if (!obj) return;
        obj->reset();
        std::lock_guard<std::mutex> lk(mx_);
        free_.push(obj);
        --active_count_;
    }

    std::size_t active()   const noexcept { return active_count_.load(std::memory_order_relaxed); }
    std::size_t capacity() const noexcept { return Capacity; }
    std::size_t available() const {
        std::lock_guard<std::mutex> lk(mx_);
        return free_.size();
    }

private:
    T*                               storage_ = nullptr;
    mutable std::mutex               mx_;
    std::stack<T*>                   free_;
    std::atomic<std::size_t>         active_count_{0};
};

}  // namespace mtt
