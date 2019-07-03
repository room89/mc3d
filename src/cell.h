// File: cell.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.5.1
// Last modified: 18.05.10.
// Description: Program for calculation of rarefaid flows.

#pragma once

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iostream>
//#include "free_boundary.h"
#include "exception.h"
#include "geometry.h"
#include "inner_boundary.h"
#include "particle.h"
#include "point.h"
//#include <mpi.h>
//#include "free_boundary.h
//#include "Header.h"

using namespace std;
using namespace mc3d;

namespace mc3d {
class cell {
 private:
  deque<particle> particles;  //частицы находящиеся внутри ячейки
  deque<particle> particle_buffer;  //буффер частиц вылетевших из ячейки
  point apex;         //опорная точка
  double lx, ly, lz;  //размеры
  double T;           //температура внутри ячейки
  // double u;
  // //средняя скорость частиц в ячейке, координата X double v;
  // //средняя скорость частиц в ячейке, координата Y double w;
  // //средняя скорость частиц в ячейке, координата Z
  point velocity;
  point mass_center;
  double volume;  //объем ячейки(пока не используеться)
  size_t n;       //количество частиц в ячейке
  // unsigned int total_n;
  // //
  double Kn_l;       //локальый кнутсен
  double Kn;         //
  double L;          //характерный размер
  double dt;         //
  double time;       //
  double calc_time;  //время вычисления в ячейке за один временной шаг dt
  double start_time;
  bool body_mark;   //содержание тела внутри ячейки
  unsigned int np;  //
  inner_boundary body_boundary;
  deque<cell*> neighbors;  //
  int* thread_mark;        //
 public:
  cell(void);                        //
  ~cell(void);                       //
  virtual void particle_move(void);  //перемещение частиц
  virtual void collisions(void);  //соударения между частицами
  void neighbor_sort(void);  //
  unsigned int N(void);  //возвращает количство частиц в ячейке
  double t();                                      //
  void set_size(double lx, double ly, double lz);  //установка размера ячейки
  void set_size(point dl);
  void set_apex(point a);  //установка "опорной" точки
  void set_L(double L);  //установка характерного размера ячейки(одинаковый для
                         //всех 3х измерений)
  void set_t(double t);  //
  void set_param(double S, double alpha, double T);  //
  void set_vel(point velocity);
  void set_Kn(double Kn);  //
  void set_inner_boundary(inner_boundary bound);
  void set_body_mark(bool mark);
  virtual bool initialazition(
      size_t N,
      const std::unique_ptr<geometry>&
          bbody);  //инициализация ячейки (N-количество частиц)
  virtual double generate_random(size_t N);                           //
  virtual double generate_random(unsigned int N, double T, point V);  //
  virtual double generate_free_random(
      unsigned int N, double T, point V,
      point nrml);  //генерация частиц для свободной границы
  virtual double generate_giper_free_random(
      unsigned int N, point V,
      double T);  //генерация частиц для свободной границы
  void sort(void);  //перемешение вылетевштх изячейки частиц в буфер
  deque<particle>* get_buffer(void);  //возвращает буфер ячейки
  void add_particle(
      deque<particle>* a);  //забирает частицы находящиеся внутри ячейки
  // void add_particle_nb(deque<particle> *a);
  // //
  void add_particle(particle* a, int n);                 //
  void fix_particle(unsigned int N, double T, point V);  //
  double calc_Kn();                                      //													//
  point get_apex() const;                                //
  point get_center() const;                              //
  point get_mass_center() const;
  point get_size() const;  //
  double get_Kn();         //
  double get_dt();         //
  double get_t();          //
  double get_u();          //
  double get_v();          //
  double get_w();          //
  double get_volume() const;
  point get_particle_mass_center();
  point get_velocity();
  double get_energy();  //
  double get_L();
  bool get_body_mark();
  cell* get_ptr();   //
  point calc_vel();  //
  double calc_T();   //
  double calc_volume();
  void set_dt(double dt);                                 //
  double calc_dt();                                       //
  void add_neighbor(cell* neighbor);                      //
  friend ostream& operator<<(ostream& o, const cell& c);  //
  bool write_file(ofstream* file);                        //
  bool write_file();                                      //
  void calc();                                            //
  void attach_thread_mark(int& ptr);                      //
  deque<cell*> fragmentation(const std::unique_ptr<geometry>& body);
  void clean_inner_particle(const geometry& body);

  bool _dbg_test_particle();

  friend bool operator>(const cell& a, const cell& b);
  friend bool operator<(const cell& a, const cell& b);
};
}  // namespace mc3d
