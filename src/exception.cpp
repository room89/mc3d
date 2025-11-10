#include "exception.h"

#include <typeinfo>

namespace mc3d {
namespace unusual_situations {
Exception::Exception() {}

Exception::Exception(int code) { this->code = code; }

Exception::Exception(int code, std::string error_message) {
  this->code = code;
  this->error_message = error_message;
}

Exception::Exception(int code, void* error_object_ptr) {
  std::cout << "Exception is throwing. Error code: " << code
            << "Class: " << typeid(error_object_ptr).name() << std::endl;
}

void Exception::PrintErrorMessage() { std::clog << error_message << std::endl; }

Exception::~Exception() {}
}  // namespace unusual_situations
}  // namespace mc3d
