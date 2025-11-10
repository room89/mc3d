// File: particle_array.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.3.1
// Last modified: 15.02.09.
// Description: Program for calculation of freemolecular flows.

#pragma once
#include <iostream>

#include "p_elem.h"
#include "point.h"

class particle_array {
 private:
  p_elem *data;
  p_elem *end;
  p_elem *start;
  unsigned int N;

 public:
  particle_array(void);
  ~particle_array(void);
  bool gen_array(int a);
  bool add(particle add);
  bool add(p_elem *add);
  bool add_particle_array(particle_array *add);
  bool del();
  p_elem *remove(void);
  bool go_back();
  bool go_forward();
  bool go_to_start();
  int num_of_elem();
  p_elem *get_end();
  p_elem *get_start();
  p_elem *get_data();
  unsigned int get_N();
  unsigned int calc_N();
  void genrate_random(double t, double u, double v, double w, point apex,
                      point d);

 private:
  bool clear();
  void set_end(p_elem *new_end);
  void set_start(p_elem *new_start);
  void set_data(p_elem *new_data);
};
