// File: warning.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 12.04.11.
// Description: Program for calculation of rarefied flows.

#pragma once

#include "warning.h"

namespace mc3d {
namespace unusual_situations {

warning::warning() {
  this->log_message = 0;
  this->warning_message = 0;
  this->code = warning::UNKNOWN_WARNING;
}

warning::warning(int code) { this->code = code; }

warning::~warning() {
  delete[] warning_message;
  delete[] log_message;
}
}  // namespace unusual_situations
}  // namespace mc3d