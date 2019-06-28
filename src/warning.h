// File: warning.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 12.04.11.
// Description: Program for calculation of rarefied flows.

//В данном файле обявляются класс предепреждений

#pragma once

namespace mc3d {
namespace unusual_situations {
class warning {
 private:
  unsigned int level;
  int code;
  char* warning_message;
  char* log_message;

 public:
  warning();  //конструктор неизвестного предупреждения(UNKNOWN_WARNING)
  warning(int code);  //конструктор предупреждения с заданым кодом
  warning(int code, char* warning_message);  //конструктор предупреждения с
                                             //заданым кодом и сообщением
  warning(int code, char* warning_message,
          char* log_message);  //конструктор предупреждения с заданым кодом и
                               //сообщением и записью в log-файл(пока нет
                               //поддержки log-файлов)
  ~warning();
  static const int CRITICAL_WARNING = 1;
  static const int DANGER_WARNING = 2;
  static const int NON_DANGER_WARNING = 4;
  static const int INFORMATION_WARNING = 8;
  static const int UNKNOWN_WARNING = 16;
};
}  // namespace unusual_situations
}  // namespace mc3d