// File: initial_conditions.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 19.04.11.
// Description: Program for calculation of rarefaid flows.

#pragma once

#include <deque>
#include <fstream>
#include <typeinfo>
#include "boundary.h"
#include "cell_cluster.h"

namespace mc3d {
class initial_conditions {
 private:
  std::deque<boundary> boundary;
  // friend cell_cluster;
 public:
  initial_conditions(void);
  ~initial_conditions(void);
  void create_cluster_file();
};
}  // namespace mc3d