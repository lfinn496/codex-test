#pragma once
#include <queue>
#include <mutex>
#include <condition_variable>
#include <optional>
#include <chrono>

namespace mtt {

// ─────────────────────────────────────────────────────────────────────────────
// Thread-safe blocking queue (bounded capacity).
// Used between the I/O reader thread and the tracker thread(s).
// ─────────────────────────────────────────────────────────────────────────────
template<typename T>
class BlockingQueue {
public:
    explicit BlockingQueue(std::size_t max_size = 100'000)
        : max_size_(max_size) {}

    // Push an item; blocks if full.
    void push(T item) {
        std::unique_lock lk(mx_);
        not_full_.wait(lk, [&]{ return q_.size() < max_size_ || stopped_; });
        if (stopped_) return;
        q_.push(std::move(item));
        not_empty_.notify_one();
    }

    // Push without blocking (drops if full); returns false if dropped.
    bool try_push(T item) {
        std::lock_guard lk(mx_);
        if (q_.size() >= max_size_) return false;
        q_.push(std::move(item));
        not_empty_.notify_one();
        return true;
    }

    // Pop with timeout; returns nullopt on timeout or stop.
    std::optional<T> pop(std::chrono::milliseconds timeout = std::chrono::milliseconds(100)) {
        std::unique_lock lk(mx_);
        if (!not_empty_.wait_for(lk, timeout, [&]{ return !q_.empty() || stopped_; }))
            return std::nullopt;
        if (q_.empty()) return std::nullopt;
        T item = std::move(q_.front());
        q_.pop();
        not_full_.notify_one();
        return item;
    }

    // Signal shutdown.
    void stop() {
        std::lock_guard lk(mx_);
        stopped_ = true;
        not_empty_.notify_all();
        not_full_.notify_all();
    }

    bool is_stopped() const {
        std::lock_guard lk(mx_);
        return stopped_;
    }

    std::size_t size() const {
        std::lock_guard lk(mx_);
        return q_.size();
    }

private:
    mutable std::mutex      mx_;
    std::condition_variable not_empty_, not_full_;
    std::queue<T>           q_;
    std::size_t             max_size_;
    bool                    stopped_ = false;
};

}  // namespace mtt
