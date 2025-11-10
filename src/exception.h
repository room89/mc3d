#pragma once

#include <iostream>
#include <string>

#include "exit_code.h"

namespace mc3d {
namespace unusual_situations {
class Exception {
 protected:
  int code;
  std::string error_message;

 public:
  Exception();
  Exception(int code);
  Exception(int code, std::string error_message);
  Exception(int code, void* error_object_ptr);
  virtual void PrintErrorMessage();
  ~Exception();
};
}  // namespace unusual_situations
}  // namespace mc3d
