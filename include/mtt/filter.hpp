#pragma once

#include "mtt/types.hpp"

#include <vector>

namespace mtt {

class MultiTargetTracker {
  public:
    virtual ~MultiTargetTracker() = default;
    virtual void predict_to(double timestamp_s) = 0;
    virtual void update(const UpdateBatch& batch) = 0;
    virtual std::vector<Track> tracks() const = 0;
};

} // namespace mtt
