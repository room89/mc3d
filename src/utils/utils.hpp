#pragma once

#include <cstdint>
#include <random>
#include <memory>
#include <vector>

namespace utils {

// Streams belong to a simulation/cell, never to a worker thread.
class RandomStream {
 public:
  explicit RandomStream(uint32_t seed = std::random_device{}(), uint64_t id = 0) { Seed(seed, id); }
  void Seed(uint32_t seed, uint64_t id = 0) {
    std::seed_seq sequence{seed, static_cast<uint32_t>(id), static_cast<uint32_t>(id >> 32)};
    engine_.seed(sequence);
  }
  std::mt19937& Engine() { return engine_; }
  double Uniform() { return std::uniform_real_distribution<double>(0, 1)(engine_); }
 private:
  std::mt19937 engine_;
};

class RandomStreams {
 public:
  explicit RandomStreams(uint32_t seed) : seed_(seed) {}
  // Bind streams before starting worker tasks; returned streams have stable addresses.
  std::shared_ptr<RandomStream> At(std::size_t index) {
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
