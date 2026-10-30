/*
 * PCG Random Number Generation for C.
 *
 * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * For additional information about the PCG random number generation scheme,
 * including its license and other licensing options, visit
 *
 *       http://www.pcg-random.org
 */

// Adapted from imneme/pcg-c-basic, revision
// bc39cd76ac3d541e618606bcc6e1e5ba5e5e6aa3, pcg_basic.c.
// Changes: C++ URBG wrapper, std::rotr, checked stream IDs; no global RNG.
// License text: third_party/pcg/LICENSE.txt.
#pragma once

#include <bit>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace utils {
class Pcg32 {
 public:
  using result_type = std::uint32_t;
  static constexpr std::uint64_t kMaxStreamId = (std::uint64_t{1} << 63) - 1;

  explicit Pcg32(std::uint64_t initial_state = 0, std::uint64_t stream_id = 0) {
    seed(initial_state, stream_id);
  }

  void seed(std::uint64_t initial_state, std::uint64_t stream_id = 0) {
    // Validate before mutating either word. PCG uses only 63 stream-id bits.
    if (stream_id > kMaxStreamId)
      throw std::invalid_argument("PCG32 stream id exceeds 2^63 - 1");
    state_ = 0;
    increment_ = (stream_id << 1u) | 1u;
    (*this)();
    state_ += initial_state;
    (*this)();
  }

  result_type operator()() noexcept {
    const auto old_state = state_;
    state_ = old_state * 6364136223846793005ULL + increment_;
    const auto shifted = static_cast<result_type>(((old_state >> 18u) ^ old_state) >> 27u);
    const auto rotation = static_cast<int>(old_state >> 59u);
    return std::rotr(shifted, rotation);
  }

  static constexpr result_type min() noexcept { return 0; }
  static constexpr result_type max() noexcept { return std::numeric_limits<result_type>::max(); }

 private:
  std::uint64_t state_;
  std::uint64_t increment_;
};
static_assert(sizeof(Pcg32) == 16, "PCG32 must retain its compact two-word state");
}  // namespace utils
