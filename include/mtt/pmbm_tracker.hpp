#pragma once

#include "mtt/filter.hpp"

#include <memory>

namespace mtt {

// PMBM-inspired interface. The current implementation is a lightweight scaffold
// intended for extension with MBM global hypothesis management.
class PmbmTracker final : public MultiTargetTracker {
  public:
    explicit PmbmTracker(TrackerConfig config);
    ~PmbmTracker() override;

    void predict_to(double timestamp_s) override;
    void update(const UpdateBatch& batch) override;
    std::vector<Track> tracks() const override;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace mtt
