#pragma once

#include <iostream>
#include <string>
#include "exit_code.h"

namespace mc3d {
namespace unusual_situations {
class exception {
 protected:
  int code;
  std::string error_message;

 public:
  exception(void);
  exception(int code);
  exception(int code, std::string error_message);
  exception(int code, void* error_object_ptr);
  virtual void print_error_message();
  ~exception(void);
};
}  // namespace unusual_situations
}  // namespace mc3d
