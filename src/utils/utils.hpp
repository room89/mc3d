#pragma once

#include "pcg32.hpp"

#include <cstdint>
#include <random>
#include <memory>
#include <vector>

namespace utils {

// Streams belong to a simulation/cell, never to a worker thread.
// A stream has one active consumer; concurrent operations must use distinct streams.
class RandomStream {
 public:
  inline static constexpr const char* kAlgorithm = "pcg32-xsh-rr-v1";
  explicit RandomStream(uint32_t seed = std::random_device{}(), uint64_t id = 0)
      : engine_(seed, id) {}
  void Seed(uint32_t seed, uint64_t id = 0) {
    engine_.seed(seed, id);
  }
  Pcg32& Engine() { return engine_; }
  double Uniform() { return std::uniform_real_distribution<double>(0, 1)(engine_); }
 private:
  Pcg32 engine_;
};

class RandomStreams {
 public:
  explicit RandomStreams(uint32_t seed) : seed_(seed) {}
  // Bind streams before starting worker tasks; returned streams have stable addresses.
  std::shared_ptr<RandomStream> At(std::size_t index) {
    if (index > Pcg32::kMaxStreamId)
      throw std::invalid_argument("PCG32 stream id exceeds 2^63 - 1");
    if (index >= streams_.max_size())
      throw std::length_error("Random stream index exceeds container capacity");
    if (streams_.size() <= index) streams_.resize(index + 1);
    if (!streams_[index]) streams_[index] = std::make_shared<RandomStream>(seed_, index);
    return streams_[index];
  }
 private:
  uint32_t seed_;
  std::vector<std::shared_ptr<RandomStream>> streams_;
};

// Moments of a drifting Maxwellian crossing a plane, per unit number density.
double IncomingFlux(double inward_drift, double temperature);
// Density proportional to w * exp(-(w-inward_drift)^2/(2*T)), w > 0.
double IncomingNormalSpeed(double inward_drift, double temperature, RandomStream& random);

}  // namespace utils
