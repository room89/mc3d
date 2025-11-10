// File: exit_code.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.4.0
// Last modified: 11.10.10.
// Description: Program for calculation of rarefied flows.

// В данном файле обявляются коды возврата ошибок

#pragma once

namespace mc3d {
namespace unusual_situations {
namespace exit_code {
const int MEMORY_NOT_ALLOCATED = -10;
const int CONSTRUCTOR_ERROR = -1;
const int BOUNDARY_CONSTRUCTOR_ERROR = -2;
const int FATAL_ERROR = -100;
}  // namespace exit_code
}  // namespace unusual_situations
}  // namespace mc3d