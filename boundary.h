// File: boundary.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.3.2
// Last modified: 21.04.09.
// Description: Program for calculation of freemolecular flows.

#pragma once

#include "particle_array.h"
#include "point.h"

class boundary {
 private:
 public:
  point nrml;  // ����������� ��������� �������
  point pstn;  // ��������� ��������� �������

  virtual int bondary_condition(
      particle_array *a) = 0;  // ����������� ����� ���������� ��������� �������
                               // � ������� ������ a
  void set_position(point a);  // ��������� ����� ��������� �������
  void set_normal(point a);  // ��������� ����������� ��������� �������
};
