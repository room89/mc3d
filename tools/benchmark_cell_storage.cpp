// Compile against either revision's mc3d_core; no production instrumentation.
// c++ -std=c++23 -O3 -I src tools/benchmark_cell_storage.cpp
// build/libmc3d_core.a -pthread -o build/cell_storage_benchmark
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <new>
#include <type_traits>

#include "cell_cluster.h"

namespace {
thread_local bool counting = false;
thread_local std::size_t allocations = 0;
}  // namespace
void* operator new(std::size_t size) {
  if (counting) ++allocations;
  if (void* p = std::malloc(size ? size : 1)) return p;
  throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

template <class Cell>
std::size_t SeededAllocations() {
  alignas(Cell) unsigned char storage[sizeof(Cell)];
  allocations = 0;
  counting = true;
  Cell* cell;
  if constexpr (std::is_constructible_v<Cell, uint32_t, uint64_t>)
    cell = new (storage) Cell(123, 0);
  else
    cell = new (storage) Cell(std::make_shared<utils::RandomStream>(123, 0));
  counting = false;
  const auto result = allocations;
  cell->~Cell();
  return result;
}
int main(int argc, char** argv) {
  const unsigned n = argc > 1 ? std::stoul(argv[1]) : 40;
  std::cout << "cell_bytes=" << sizeof(mc3d::Cell)
            << " seeded_cell_allocations=" << SeededAllocations<mc3d::Cell>()
            << '\n';
  mc3d::CellCluster cluster(1, 123);
  cluster.SetApex({0, 0, 0});
  cluster.SetSize(1, 1, 1);
  auto body = std::make_unique<mc3d::Geometry>();
  const auto start = std::chrono::steady_clock::now();
  cluster.Initialize(n, n, n, 0, .04, .6, std::move(body), 0, .05, 1.2);
  const auto elapsed = std::chrono::duration<double, std::milli>(
                           std::chrono::steady_clock::now() - start)
                           .count();
  std::cout << "initialization_ms=" << elapsed << '\n';
}
