// File: warning.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 12.04.11.
// Description: Program for calculation of rarefied flows.

#include "warning.h"

namespace mc3d {
namespace unusual_situations {

Warning::Warning() {
  this->log_message = 0;
  this->warning_message = 0;
  this->code = Warning::UNKNOWN_WARNING;
}

Warning::Warning(int code) { this->code = code; }

Warning::~Warning() {
  delete[] warning_message;
  delete[] log_message;
}
}  // namespace unusual_situations
}  // namespace mc3d