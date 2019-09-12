#pragma once

#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>

#include <cell.h>
#include <exit_code.h>
#include <free_boundary.h>
#include <geometry.h>
#include <giper_free_boundary.h>
#include <inner_cell.h>
#include <mirror_boundary.h>
#include <pereodic_boundary.h>
#include <point.h>
#include <utils/threadpool.hpp>

const int NUM_CPU = 8;

namespace mc3d {

class CellCluster {
 public:
  inline void calc_dt();  //вычисление шага по времени для данного кластера
  inline void set_dt_in_cells(
      double dt);  //установка шага по времени во все ячейки
  inline bool time_step();
  CellCluster(void);   //конструктор по умолчанию
  ~CellCluster(void);  //деструктор
  void set_apex(point apex);  //задание опорной точки кластера
  void set_size(double Lx, double Ly, double Lz);  //задание размеров кластера
  void set_end_time(double t_end);  //задание времени окончания рассчета
  void set_data_save_dtime(double data_dt);
  //инициализация кластера(параметры:	ncx, ncy, ncz	- количество ячеек по
  //координатам, np - уровень статистик, Kn - Кнудсен
  bool initialazition(unsigned int ncx, unsigned int ncy, unsigned int ncz,
                      double density, double Kn, double Cu,
                      std::unique_ptr<geometry>&& body_, double S, double alpha,
                      double T);
  bool write_cell_file(const string& init_file);
  bool WriteFile();  //запись данных в файл(плотность и энэргия), данные
                     //собираються со всех кластеров
  bool write_speed_file();  //запись данных в файл(скорость), данные собираються
                            //со всех кластеров
  bool WriteFile(const string& file_name);  //запись данных в файл(плотность и
                                            //энэргия) с заданным именем, данные
                                            //собираються со всех кластеров
  bool write_speed_file(const char* file_name);
  bool write_times();  //запись данных о врем ени выполнения каждой итерации
  void set_boundary_condition(
      boundary** a, int i);  //установка граничных условий a - указатель на
                             //массив, i - количество элементов в массиве.
  void set_boundary_condition(
      boundary* a);  //установка граничных условий a - указатель на гр. условие
  void computation(void);  //вычисления
  void cells_fragmentation();
  void cells_test();
  void clean_inner_particle();

 private:
  int numproc;  //количество процессов(не путать с нитями), специально для MPI
  int proc_id;  //номер процесса
  unsigned int i_j_k;  //переменная используется при инициализации кластера|
                       // todo: удалить
  unsigned int np;     //уровень статистики
  unsigned int ncx, ncy, ncz;  //начальное количество ячеек по x, y, z
  unsigned int N;  //начальное количество частиц в кластере(в ячейках кластера)
  unsigned int step;
  double Lx, Ly, Lz;  //начальные размеры кластера
  point apex;  //координаты опорной точки кластера.
  double Kn;   //Кнудсен
  double Cu;   //число Куранта
  double t;    //время с начала расчёта(физическое)
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
  deque<cell> cells;  //двусвязный список ячеек кластера.
  deque<cell>::iterator cell_iter;  //итератор для обхода ячеек
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
  inline void sync_dt();  //синхронизация шага по времени. тоже самое что и
                          //методы: send_dt() и recv_dt(), только в одном методе
  ofstream log;  //переменная для вывода логов, пока не реализовано
  std::unique_ptr<geometry> body_;
  double data_dt;
  double data_t;
  double density;
  double volume;

  utils::ThreadPool thread_pool_;
};
};  // namespace mc3d
