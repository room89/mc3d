#include <typeinfo>

#include "exception.h"

namespace mc3d {
namespace unusual_situations {
exception::exception(void) {}

exception::exception(int code) { this->code = code; }

exception::exception(int code, std::string error_message) {
  this->code = code;
  this->error_message = error_message;
}

exception::exception(int code, void* error_object_ptr) {
  std::cout << "Exception is throwing. Error code: " << code
            << "Class: " << typeid(error_object_ptr).name() << std::endl;
}

void exception::print_error_message() {
  std::clog << error_message << std::endl;
}

exception::~exception(void) {}
}  // namespace unusual_situations
}  // namespace mc3d
