// File: cell_cluster.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.5.1
// Last modified: 16.12.10.
// Description: Program for calculation of rarefaid flows.

#pragma once

//#include "mpi.h"
//#include "mpicxx.h"
//#include "boost\mpi.hpp"
#include "cell.h"
#include "inner_cell.h"
#include "point.h"
//#include "Header.h"
#include "exit_code.h"
#include "free_boundary.h"
#include "geometry.h"
#include "giper_free_boundary.h"
#include "mirror_boundary.h"
#include "pereodic_boundary.h"
//#include "typeinfo.h"
#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>

const int NUM_CPU = 4;

/*
*** 	 cell_cluster - основной класс, в нем происходит инициализация, обмен
*частицами между ячейками, выполнение граничных условий, обмен между
*кластерами(если их больше одного)
*/
namespace mc3d {
class cell_cluster {
 private:
  int numproc;  //количество процессов(не путать с нитями), специально для MPI
  int proc_id;  //номер процесса
  unsigned int i_j_k;  //переменная используется при инициализации кластера|
                       //todo: удалить
  unsigned int np;     //уровень статистики
  unsigned int ncx, ncy, ncz;  //начальное количество ячеек по x, y, z
  unsigned int N;  //начальное количество частиц в кластере(в ячейках кластера)
  unsigned int step;
  double Lx, Ly, Lz;  //начальные размеры кластера
  point apex;  //координаты опорной точки кластера.
  mutable int*
      thread_mark;  //указатель на массив флагов, определяющих состояние нити
  double Kn;        //Кнудсен
  double Cu;        //число Куранта
  double t;  //время с начала расчёта(физическое)
  double t_end;  //время окончания рассчета
  double dt;  //шаг по времени(должен быть равен во всех кластерах)
  double begin_time;  //время начала расчета(для вычисления времени работы)
  double begin_iter_time;  //время начала выполнения шага по времени
  deque<double> times;  //дэк с временами выполнения шага по времени
  deque<double>
      send_times;  //дэк с временами синхронизации данных между кластерами
  deque<double> bound_times;  //дэк с временами выполнения граничных условий
  deque<double> sort_times;  //дэк с временами сортировки частиц по ячейкам
  deque<double>
      calc_times;  //дэк с временами расчета соударений между частицами
  deque<cell*> cells;  //двусвязный список ячеек кластера.
  deque<cell*>::iterator cell_iter;  //итератор для обхода ячеек
  deque<particle> partile_buffer;  //двусвязный список буфера кластера
  inline void boundary_condition(void);  //выполнение граничных условий
  deque<boundary*> boundary_cond_outer;  //дэк с внешними граничными условиями
  deque<boundary*> boundary_cond_inner;  //дэк с внутреними граничными условиями
  inline bool
  send_data();  //посылка данных на все класстеры за исключением кластера из
                //которого посылают. не используется, возможно стоит удалить
  inline bool recv_data();  //прием данных с класстеров. не используеться,
                            //возможно стоит удалить
  inline void sync_data();  //синхронизация данных между кластерами. замена для
                            //методов:  send_data(), recv_data()
  inline bool
  send_dt();  //передача шага по времени для синхронизации его на всех кластерах
  inline bool recv_dt();  //прием шага по времени, и выбор наименьшего
  inline void sync_dt();  //синхронизация шага по времени. тоже самое что и
                          //методы: send_dt() и recv_dt(), только в одном методе
  ofstream log;  //переменная для вывода логов, пока не реализовано
  geometry* body;
  double data_dt;
  double data_t;
  double density;
  double volume;

 public:
  double* send_buffer_u;  //буфер для передачи частиц кластеру, скорость по x
  double* send_buffer_v;  //буфер для передачи частиц кластеру, скорость по y
  double* send_buffer_w;  //буфер для передачи частиц кластеру, скорость по z
  double* send_buffer_x;  //буфер для передачи частиц кластеру, координата x
  double* send_buffer_y;  //буфер для передачи частиц кластеру, координата y
  double* send_buffer_z;  //буфер для передачи частиц кластеру, координата z
  double* recv_buffer_u;  //буфер для приёма частиц кластером, скорость по x
  double* recv_buffer_v;  //буфер для приёма частиц кластером, скорость по y
  double* recv_buffer_w;  //буфер для приёма частиц кластером, скорость по z
  double* recv_buffer_x;  //буфер для приёма частиц кластером, координата x
  double* recv_buffer_y;  //буфер для приёма частиц кластером, координата y
  double* recv_buffer_z;  //буфер для приёма частиц кластером, координата w
  inline void calc_dt();  //вычисление шага по времени для данного кластера
  inline void set_dt_in_cells(
      double dt);  //установка шага по времени во все ячейки
  inline bool time_step(
      void);  //выполнение шага по времени(соударения, перемещения, сортировка и
              //синхронизация частиц)
  cell_cluster(void);   //конструктор по умолчанию
  ~cell_cluster(void);  //деструктор
  void set_apex(point apex);  //задание опорной точки кластера
  void set_size(double Lx, double Ly, double Lz);  //задание размеров кластера
  void set_end_time(double t_end);  //задание времени окончания рассчета
  void set_data_save_dtime(double data_dt);
  bool initialazition(unsigned int ncx, unsigned int ncy, unsigned int ncz,
                      unsigned int np, double Kn, double Cu, geometry* body,
                      int numproc = 1, int my_id = 0);
  //инициализация кластера(параметры:	ncx, ncy, ncz	- количество ячеек по
  //координатам, np - уровень статистик, Kn - Кнудсен
  bool initialazition(unsigned int ncx, unsigned int ncy, unsigned int ncz,
                      double density, double Kn, double Cu, geometry* body,
                      bool fragmentation, double S, double alpha, double T,
                      int numproc = 1, int my_id = 0);
  bool initialazition(const char* init_file, double Kn, double Cu, double L,
                      geometry* body, int numproc = 1, int my_id = 0);
  // bool read_cell_file(const char *init_file);
  bool write_cell_file(const char* init_file);
  bool write_file();  //запись данных в файл(плотность и энэргия), данные
                      //собираються со всех кластеров
  bool write_speed_file();  //запись данных в файл(скорость), данные собираються
                            //со всех кластеров
  bool write_file(const char* file_name);  //запись данных в файл(плотность и
                                           //энэргия) с заданным именем, данные
                                           //собираються со всех кластеров
  bool write_speed_file(const char* file_name);
  bool write_times();  //запись данных о врем ени выполнения каждой итерации
  void set_boundary_condition(
      boundary** a, int i);  //установка граничных условий a - указатель на
                             //массив, i - количество элементов в массиве.
  void set_boundary_condition(
      boundary* a);  //установка граничных условий a - указатель на гр. условие
  void set_geometry(geometry* body);
  void computation(void);  //вычисления
  void cells_fragmentation();
  void cells_test();
  void clean_inner_particle();
};
};  // namespace mc3d
