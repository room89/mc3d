// File: cell_cluster.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.5.1
// Last modified: 16.12.10.
// Description: Program for calculation of rarefaid flows.

#pragma once

#include <thread>

#include "cell_cluster.h"

using namespace std;
using namespace mc3d;

namespace mc3d {
/*double V(double x, double y, double R)
{
  return 0.3 * (x * x + y * y) * std::exp(0.5 * (1 - (x * x + y * y) / (R * R)))
/ R;
}*/

cell_cluster::cell_cluster(void) {
  Kn = 0;
  np = 0;
  ncx = 0;
  ncy = 0;
  ncz = 0;
  N = 0;
  Lx = 0;
  Ly = 0;
  Lz = 0;
  t = 0;
  dt = 10000;
  step = 0;

  thread_mark = new int[NUM_CPU];

  for (int i = 0; i < NUM_CPU; i++) {
    thread_mark[i] =
        0;  //устанавилваем для каждой нити состояние готовности рассчета
  }
  /*char *log_name = "proc00.log";
  char a = '0' + proc_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  log.open(log_name);*/
}

cell_cluster::~cell_cluster(void) {
  cell_iter = cells.begin();
  while (cell_iter != cells.end()) {
    delete *cell_iter;
    cell_iter++;
  }
  cells.clear();
  // log.close();
  // MPI_File_close(&mpi_file);
}

bool cell_cluster::initialazition(unsigned int ncx, unsigned int ncy,
                                  unsigned int ncz, double density, double Kn,
                                  double Cu, geometry* body, bool fragmentation,
                                  double S, double alpha, double T, int numproc,
                                  int proc_id) {
  this->numproc = numproc;
  this->proc_id = proc_id;
  this->body = body;

  srand((unsigned int)time(0));

  // char  _fname[] = "mpi_data.dat";

  int mpi_init = 1;

  //    MPI_Initialized(&mpi_init);

  if (!mpi_init)  //проверка инициализации MPI
  {
    int arg1 = 0;
    char** arg2 = NULL;
    //      MPI_Init(&arg1, &arg2);                     //инициализируем MPI
    //      MPI_Comm_size(MPI_COMM_WORLD, &numproc);    //устанавливаем кол-во
    //      процессов MPI_Comm_rank(MPI_COMM_WORLD, &proc_id); //устанавливаем
    //      номер процесса
  }

  //    begin_time = MPI_Wtime();	                    //устанавливаем
  //    время начала работы

  begin_iter_time = begin_time;  //устанавливаем время начала нулевой
                                 //итерации(нулевая итерация - инициализация)

  // MPI_File_open(MPI_COMM_WORLD, mpi_fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
  // MPI_INFO_NULL, &mpi_file);

  this->ncx = ncx;
  this->ncy = ncy;
  this->ncz = ncz;

  this->density = density;

  this->Kn = Kn;
  this->Cu = Cu;

  unsigned int N = 0;

  double dx = Lx / double(ncx);  //устанавливаем размер ячейки
  double dy = Ly / double(ncy);  //устанавливаем размер ячейки
  double dz = Lz / double(ncz);  //устанавливаем размер ячейки

  this->t = 0;  //устанавливаем "физическое" время на ноль
  this->dt = 100000;

  double fragm_dist = std::min(dx, min(dy, dz)) * 3.1;

  cell_iter = cells.begin();
  for (int i = 0; i < int(ncx);
       i++)  //цыклы генерации ячеек	по 3м координатам
  {
    for (int j = 0; j < int(ncy); j++) {
      for (int k = 0; k < int(ncz); k++) {
        if (fragmentation) {
          point a, dc(dx / 2, dy / 2, dz / 2);
          a.x = apex.x + i * dx;
          a.y = apex.y + j * dy;
          a.z = apex.z + k * dz;

          if (body != NULL) {
            // double di = body->dist_to_point(a + dc);
            if (body->dist_to_point(a + dc) < fragm_dist) {
              // point a(0, 0, 0);
              point dl(dx / 2, dy / 2, dz / 2);

              unsigned int N_ =
                  static_cast<unsigned int>(density * dl.volume());

              cell* temp_cells[8];

              for (int i = 0; i < 8; i++) {
                temp_cells[i] = new cell;
                temp_cells[i]->set_param(S, alpha, T);
              }

              for (int i = 0; i < 8; i++) {
                temp_cells[i]->set_size(dl);
              }

              point cell_apex = a;
              temp_cells[0]->set_apex(cell_apex);

              point shift(dx / 2, 0, 0);

              cell_apex = a + shift;
              temp_cells[1]->set_apex(cell_apex);

              shift.set(0, dy / 2, 0);
              cell_apex = a + shift;
              temp_cells[2]->set_apex(cell_apex);

              shift.set(0, 0, dz / 2);
              cell_apex = a + shift;
              temp_cells[3]->set_apex(cell_apex);

              shift.set(dx / 2, dy / 2, 0);
              cell_apex = a + shift;
              temp_cells[4]->set_apex(cell_apex);

              shift.set(dx / 2, 0, dz / 2);
              cell_apex = a + shift;
              temp_cells[5]->set_apex(cell_apex);

              shift.set(0, dy / 2, dz / 2);
              cell_apex = a + shift;
              temp_cells[6]->set_apex(cell_apex);

              shift.set(dx / 2, dy / 2, dz / 2);
              cell_apex = a + shift;
              temp_cells[7]->set_apex(cell_apex);

              for (int i = 0; i < 8; i++) {
                temp_cells[i]->initialazition(N_, body);
                temp_cells[i]->set_L(
                    Lx);  //задание характерного размера для ячейки
                temp_cells[i]->set_Kn(Kn);  //задание кнудсена для ячейки
                // temp_cells[i]->set_body_mark(true);

                cell_iter = cells.insert(cell_iter, temp_cells[i]);
              }

              N += N_ * 8;
            } else {
              cell* temp_cell;

              temp_cell = new cell;

              unsigned int i_j_k = k * ncx * ncy + j * ncx + i;

              unsigned int N_;  //количество частиц в добавляймой ячейке
              point a;
              a.x = apex.x + i * dx;
              a.y = apex.y + j * dy;
              a.z = apex.z + k * dz;

              N_ = static_cast<unsigned int>(density * dx * dy * dz);
              temp_cell->set_param(S, alpha, T);

              temp_cell->set_apex(a);  //задание опорной точки
              temp_cell->set_size(dx, dy, dz);  //задание размеров

              temp_cell->initialazition(
                  N_,
                  body);  //инициализация ячейки, с заданным количеством частиц

              temp_cell->set_L(Lx);  //задание характерного размера для ячейки
              temp_cell->set_Kn(Kn);  //задание кнудсена для ячейки

              cell_iter = cells.insert(
                  cell_iter, temp_cell);  //копирование и вставка подготовленной
                                          //ячейки в дэк класстера

              N += N_;

              double dtt = (*cell_iter)->get_dt();
              dt = min(
                  this->dt,
                  dtt);  //вычисляем минимальный шаг по времени для всех ячеек
            }
          }
        } else {
          cell* temp_cell;

          temp_cell = new cell;

          unsigned int i_j_k = k * ncx * ncy + j * ncx + i;

          unsigned int N_;  //количество частиц в добавляймой ячейке
          point a;
          a.x = apex.x + i * dx;
          a.y = apex.y + j * dy;
          a.z = apex.z + k * dz;

          N_ = static_cast<unsigned int>(density * dx * dy * dz);
          temp_cell->set_param(S, alpha, T);

          temp_cell->set_apex(a);  //задание опорной точки
          temp_cell->set_size(dx, dy, dz);  //задание размеров

          temp_cell->initialazition(
              N_, body);  //инициализация ячейки, с заданным количеством частиц

          temp_cell->set_L(Lx);  //задание характерного размера для ячейки
          temp_cell->set_Kn(Kn);  //задание кнудсена для ячейуи

          cell_iter = cells.insert(
              cell_iter, temp_cell);  //копирование и вставка подготовленной
                                      //ячейки в дэк класстера

          N += N_;

          double dtt = (*cell_iter)->get_dt();
          dt = min(this->dt,
                   dtt);  //вычисляем минимальный шаг по времени для всех ячеек
        }
      }
    }
  }

  cells_test();

  // cells_fragmentation();

  // cells_test();

  cout << endl << "dobavlyaem \"sosedey\"" << endl;
  cout << endl << "|--------------------------------------|" << endl;
  int i = 0, j = 1;
  //поиск "соседей"(соседних ячеек) для каждой ячейки:
  cell_iter = cells.begin();
  while (cell_iter != cells.end()) {
    i++;
    if (int(i * 41 / cells.size()) > j) {
      cout << "*" << flush;
      j++;
    }

    deque<cell*>::iterator cell_iter2 = cells.begin();

    point cell_c1 = (*cell_iter)->get_center();

    while (cell_iter2 != cells.end()) {
      point cell_c2 = (*cell_iter2)->get_center();

      double LL = std::sqrt((cell_c1 - cell_c2) * (cell_c1 - cell_c2));

      if ((LL <= 1.8 * dx) && (*cell_iter != *cell_iter2))
        (*cell_iter)->add_neighbor(*cell_iter2);  //добавление "соседей"

      cell_iter2++;
    }

    cell_iter++;
  }

  cout << endl;

  this->N = N;

  this->density = N;

  cell_iter = cells.begin();
  while (cell_iter != cells.end()) {
    (*cell_iter)->clean_inner_particle(body);
    cell_iter++;
  }

  /*cell_iter = cells.begin();
  while(cell_iter != cells.end())
  {
    if(*cell_iter == NULL)
    {
      cell_iter = cells.erase(cell_iter);
    }
    else cell_iter++;
  }*/

  this->sync_dt();
  cout << "Cluster inicializirovan" << endl;

  return true;
}

bool cell_cluster::initialazition(unsigned int ncx, unsigned int ncy,
                                  unsigned int ncz, unsigned int np, double Kn,
                                  double Cu, geometry* body, int numproc,
                                  int proc_id) {
  this->numproc = numproc;
  this->proc_id = proc_id;
  this->body = body;

  srand((unsigned int)time(0));

  // char mpi_fname[] = "mpi_data.dat";

  int mpi_init = 1;

  //    MPI_Initialized(&mpi_init);

  if (!mpi_init)  //проверка инициализации MPI
  {
    int arg1 = 0;
    char** arg2 = NULL;
    //      MPI_Init(&arg1, &arg2);                     //инициализируем MPI
    //      MPI_Comm_size(MPI_COMM_WORLD, &numproc);    //устанавливаем кол-во
    //      процессов MPI_Comm_rank(MPI_COMM_WORLD, &proc_id); //устанавливаем
    //      номер процесса
  }

  //    begin_time = MPI_Wtime();	                    //устанавливаем
  //    время начала работы

  begin_iter_time = begin_time;  //устанавливаем время начала нулевой
                                 //итерации(нулевая итерация - инициализация)

  // MPI_File_open(MPI_COMM_WORLD, mpi_fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
  // MPI_INFO_NULL, &mpi_file);

  this->ncx = ncx;
  this->ncy = ncy;
  this->ncz = ncz;

  this->np = np;

  this->Kn = Kn;
  this->Cu = Cu;

  unsigned int N = 0;

  double dx = Lx / double(ncx);  //устанавливаем размер ячейки
  double dy = Ly / double(ncy);  //устанавливаем размер ячейки
  double dz = Lz / double(ncz);  //устанавливаем размер ячейки

  this->t = 0;  //устанавливаем "физическое" время на ноль
  this->dt = 100000;

  cell* temp_cell;

  cell_iter = cells.begin();
  for (int i = 0; i < int(ncx);
       i++)  //цыклы генерации ячеек	по 3м координатам
  {
    for (int j = 0; j < int(ncy); j++) {
      for (int k = 0; k < int(ncz); k++) {
        temp_cell = new cell;

        unsigned int i_j_k = k * ncx * ncy + j * ncx + i;

        unsigned int N_;  //количество частиц в добавляймой ячейке
        point a;
        a.x = apex.x + i * dx;
        a.y = apex.y + j * dy;
        a.z = apex.z + k * dz;
        // N_ = unsigned int(np * (1 + 10 * F(a.x + 0.5 * dx, a.y + 0.5 * dy,
        // a.z + 0.5 * dz))); N_ = unsigned int(np * (F2(a.x + 0.5 * dx, a.y +
        // 0.5 * dy, .05)));

        N_ = np;
        temp_cell->set_param(10, 0, 1);

        temp_cell->set_apex(a);  //задание опорной точки
        temp_cell->set_size(dx, dy, dz);  //задание размеров

        temp_cell->initialazition(
            N_, body);  //инициализация ячейки, с заданным количеством частиц

        temp_cell->set_L(Lx);  //задание характерного размера для ячейки
        temp_cell->set_Kn(Kn);  //задание кнудсена для ячейуи

        // cout << "New cell vel: " << temp_cell->calc_vel() <<endl;

        /*N_ = (unsigned int)(np * (1. + 0.1 * sin(2. * Pi * (a.x + 0.5 * dx +
        0.5)) * sin(2. * Pi * (a.y + 0.5 * dy + 0.5)))); temp_cell->set_vel(0.5
        * sin(2. * Pi * (a.x + 0.5 * dx + 0.5)) * cos(2. * Pi * (a.y + 0.5 * dy
        + 0.5)), -0.4 * cos(2. * Pi * (a.x + 0.5 * dx + 0.5)) * sin(2. * Pi *
        (a.y + 0.5 * dy + 0.5)), 0); temp_cell->set_t(1. + 0.1 * cos(2. * Pi *
        (a.x + 0.5*dx + 0.5)) * cos(2. * Pi * (a.y + 0.5*dy + 0.5)));*/

        /*N_ = (unsigned int)(np * (1. + 0.1 * sin(2. * Pi * (a.x + 0.5*dx +
        0.5)) * sin(2. * Pi * (a.y + 0.5*dy + 0.5)))); temp_cell->set_vel(0.2 -
        (a.y + 0.5 * dy) * mc3d::V(a.x + 0.5 * dx, a.y + 0.5*dy, 0.1) / ((a.x +
        0.5*dx) * (a.x + 0.5*dx) + (a.y + 0.5*dy) * (a.y + 0.5*dy)), 0.35 + (a.x
        + 0.5*dx) * mc3d::V(a.x + 0.5*dx, a.y + 0.5*dy, 0.1) / ((a.x + 0.5*dx) *
        (a.x + 0.5*dx) + (a.y + 0.5*dy) * (a.y + 0.5*dy)), 0);
        temp_cell->set_t(1 - 0.2 * mc3d::V(a.x + 0.5*dx, a.y + 0.5*dy, 0.1) *
        mc3d::V(a.x + 0.5*dx, a.y + 0.5*dy, 0.1));*/

        cell_iter = cells.insert(
            cell_iter, temp_cell);  //копирование и вставка подготовленной
                                    //ячейки в дэк класстера

        N += N_;

        double dtt = (*cell_iter)->get_dt();
        dt = min(this->dt,
                 dtt);  //вычисляем минимальный шаг по времени для всех ячеек
      }
    }
  }

  this->set_dt_in_cells(dt);  //устанавливаем шаг по времени

  deque<cell*> cells_swap;

  cells_test();

  // cells_fragmentation();

  // cells_test();

  //поиск "соседей"(соседних ячеек) для каждой ячейки:
  cell_iter = cells.begin();
  while (cell_iter != cells.end()) {
    deque<cell*>::iterator cell_iter2 = cells.begin();

    point cell_c1 = (*cell_iter)->get_center();

    while (cell_iter2 != cells.end()) {
      point cell_c2 = (*cell_iter2)->get_center();

      double LL = std::sqrt((cell_c1 - cell_c2) * (cell_c1 - cell_c2));

      if ((LL <= 1.8 * dx) && (*cell_iter != *cell_iter2))
        (*cell_iter)->add_neighbor(*cell_iter2);  //добавление "соседей"

      cell_iter2++;
    }

    cell_iter++;
  }

  this->N = N;

  this->density = N;

  cell_iter = cells.begin();
  while (cell_iter != cells.end()) {
    if (*cell_iter == NULL) {
      cell_iter = cells.erase(cell_iter);
    } else
      cell_iter++;
  }

  this->sync_dt();

  return 0;
}

void cell_cluster::set_apex(point apex) { this->apex = apex; }

void cell_cluster::set_size(double Lx, double Ly, double Lz) {
  this->Lx = Lx;
  this->Ly = Ly;
  this->Lz = Lz;
}

bool cell_cluster::time_step() {
  //    times.push_back(MPI_Wtime() - begin_iter_time);

  //    begin_iter_time = MPI_Wtime();

  sync_dt();

  if (!proc_id) cout << "dt = " << dt << endl;

  cell_iter = cells.begin();

  // point cell_vel = point(0, 0, 0);

  // cout << "first cell vel: " << cell_vel << endl;

  deque<particle>* cell_buffer;

  //    double start_section_time = MPI_Wtime();
  int cell_number = 0;
  while (cell_iter != cells.end()) {
    // cell_vel += (*cell_iter)->calc_vel();
    // for(int i = 0; i < NUM_CPU; i++)
    //{
    //	if(cell_iter == cells.end()) break;
    //	if(thread_mark[i] == 0)															//проверка
    //не занятости нити(точнее проверка возможности запустить нить(максимальное
    //количество нитей - константа NUM_CPU))
    //	{
    //		thread_mark[i] = 2;
    ////устанавливаем семафор нити в состояние вычисления
    //		(*cell_iter)->attach_thread_mark(thread_mark[i]);
    ////привязываем флаг(семафор)	к ячейке. 		my_thread[i] =
    ///new
    // boost::thread(&cell::calc, (*cell_iter)->get_ptr());
    // cell_iter++;
    //	}
    //}
    // for(int i = 0; i < NUM_CPU; i++)
    //{
    //	if(cell_iter == cells.end()) break;
    //	if(thread_mark[i] == 1)															//проверка
    //на окончания вычисления нити(1 - нить закончила вычисления и готова к
    //завершению работы)
    //	{
    //		thread_mark[i] = 0;
    //		my_thread[i]->join();
    //		delete my_thread[i];
    //	}
    //}
    std::vector<std::thread> threads(NUM_CPU);
    for (size_t i = 0; i < NUM_CPU; i++) {
      if (cell_iter == cells.end()) break;
      threads[i] = std::thread(&cell::calc, (*cell_iter)->get_ptr());
      //        threads.add_thread(my_thread[i]);
      //(*cell_iter)->_dbg_test_particle();
      cell_iter++;
    }
    for (auto& thread : threads) thread.join();
    // delete [] my_thread;
  }

  // cout << "cells vel: " << cell_vel / cells.size() << endl;

  //    calc_times.push_back(MPI_Wtime() - start_section_time);

  //    start_section_time = MPI_Wtime();

  //сортировка частиц по ячейкам и буфферу кластера:
  cell_iter = cells.begin();
  while (cell_iter != cells.end()) {
    //(*cell_iter)->_dbg_test_particle();
    (*cell_iter)->calc_vel();
    (*cell_iter)->sort();
    (*cell_iter)->neighbor_sort();
    cell_buffer = (*cell_iter)->get_buffer();
    partile_buffer.insert(partile_buffer.end(), cell_buffer->begin(),
                          cell_buffer->end());
    cell_buffer->clear();
    cell_iter++;
  }

  cell_iter = cells.begin();
  while (cell_iter != cells.end()) {
    (*cell_iter)->add_particle(&partile_buffer);
    cell_iter++;
  }
  //окончание сортировки.

  //    sort_times.push_back(MPI_Wtime() - start_section_time);

  //    start_section_time = MPI_Wtime();

  boundary_condition();

  //сортировака частиц оставшихся в буффере кластера после граничных условий
  cell_iter = cells.begin();
  while (cell_iter != cells.end()) {
    (*cell_iter)->add_particle(&partile_buffer);
    cell_iter++;
  }
  //окончание сортировки.

  //    double section_time = MPI_Wtime() - start_section_time;
  //    bound_times.push_back(section_time);

  //    start_section_time = MPI_Wtime();
  // sync_data();
  //    send_times.push_back(MPI_Wtime() - start_section_time);

  t += dt;

  if (!proc_id) cout << "t = " << t << endl;
  //    double work_time = MPI_Wtime() - begin_time;
  //    if(!proc_id) cout << "Work time: " << int(work_time / 60) << "(" <<
  //    int(100 * t / t_end) << "%)\tTime to end: " << int((work_time) * (t_end
  //    / t - 1)) / 60 << "( " << int((t_end - t) * section_time / dt) / 60 <<
  //    ") min "
  //      << int((work_time) * (t_end  / t - 1)) % 60 << " sec." << endl;

  partile_buffer
      .clear();  //очистка буффера кластера, т.к. в нем остались только
                 //"лишние"(напр. дублированые в других кластерах) частицы

  if (data_dt > 0) {
    if (t > data_t) {
      std::string file_name = "data t=";
      std::ostringstream ost;
      ost << t << ".dat";
      file_name += ost.str();
      data_t += data_dt;
      write_file(file_name.c_str());
      /*file_name.clear();
      file_name = "velocity t=";
      file_name += ost.str();
      write_speed_file(file_name.c_str());*/
    }
  }

  //сортировка ячеек по времени рассчета
  if (!step % 3) {
    // std::rotate(cells.begin(), cells.begin() + 1 + cells.size() / (NUM_CPU *
    // 2), cells.end());

    std::vector<std::thread> threads;

    for (size_t i = 0; i < NUM_CPU; i++) {
      int hj = cells.size() / NUM_CPU;
      auto sort_iter = cells.begin() + i * cells.size() / NUM_CPU;
      threads[i] = std::thread(&sort<deque<cell*>::iterator>, sort_iter,
                               sort_iter + cells.size() / NUM_CPU);
    }
    for (auto& thread : threads) thread.join();
  }

  //окончание сортировки ячеек по времени рассчета

  if (t >= t_end)
    return true;
  else
    return false;
}

void cell_cluster::boundary_condition() {
  deque<boundary*>::iterator a = boundary_cond_outer.begin();
  while (a != boundary_cond_outer.end()) {
    (*a)->bondary_condition(&partile_buffer, dt);
    a++;
  }
}

bool cell_cluster::send_data() {
  if (numproc > 0) {
    int n = partile_buffer.size();

    //буффер используются т.к. mpi не "умеет" посылать deque, соответственно
    //создаем буффер для обмена
    send_buffer_u = new double[n];
    send_buffer_v = new double[n];
    send_buffer_w = new double[n];
    send_buffer_x = new double[n];
    send_buffer_y = new double[n];
    send_buffer_z = new double[n];

    int j = 0;
    deque<particle>::iterator particle_iter = partile_buffer.begin();
    while (particle_iter != partile_buffer.end()) {
      send_buffer_u[j] = particle_iter->velocity.x;
      send_buffer_v[j] = particle_iter->velocity.y;
      send_buffer_w[j] = particle_iter->velocity.z;
      send_buffer_x[j] = particle_iter->velocity.x;
      send_buffer_y[j] = particle_iter->velocity.y;
      send_buffer_z[j] = particle_iter->velocity.z;

      j++;

      particle_iter++;
    }

    for (int i = 0; i < numproc; i++) {
      if (i == proc_id) continue;
      //        MPI_Send(send_buffer_u, n, MPI_DOUBLE, i, 100, MPI_COMM_WORLD);
      //        MPI_Send(send_buffer_v, n, MPI_DOUBLE, i, 101, MPI_COMM_WORLD);
      //        MPI_Send(send_buffer_w, n, MPI_DOUBLE, i, 102, MPI_COMM_WORLD);
      //        MPI_Send(send_buffer_x, n, MPI_DOUBLE, i, 103, MPI_COMM_WORLD);
      //        MPI_Send(send_buffer_y, n, MPI_DOUBLE, i, 104, MPI_COMM_WORLD);
      //        MPI_Send(send_buffer_z, n, MPI_DOUBLE, i, 105, MPI_COMM_WORLD);
    }

    partile_buffer.clear();
    delete[] send_buffer_u;
    delete[] send_buffer_v;
    delete[] send_buffer_w;
    delete[] send_buffer_x;
    delete[] send_buffer_y;
    delete[] send_buffer_z;
  }
  return true;
}

bool cell_cluster::recv_data() {
  if (numproc > 1) {
    //      MPI_Status *status = NULL;
    bool* step1 = NULL;
    bool* step2 = NULL;
    int* recv_elem = NULL;

    try {
      bool recv = true;  // true - прием сообщений, false - сообщения приняты.
      step1 = new bool[numproc];  //флаг готовности приема сообщения
      step2 =
          new bool[numproc];  //флаг приема сообщения(true - сообщение принято)
      int k;
      recv_elem = new int[numproc];  //размер принимаемого сообщения

      for (int i = 0; i < numproc; i++) {
        step1[i] = step2[i] = false;
      }

      step1[proc_id] = true;
      step2[proc_id] = true;

      while (recv) {
        for (int i = 0; i < numproc; i++) {
          if (!step1[i]) {
            //              MPI_Iprobe(i, 100, MPI_COMM_WORLD, &k, &status[i]);
            //              //k == true - сообщение готово к приему(тип int)
            step1[i] = (bool)k;
          }
          if (step1[i] && !step2[i]) {
            //              MPI_Get_count(&status[i], MPI_DOUBLE,
            //              &recv_elem[i]);

            recv_buffer_u = new double[recv_elem[i]];
            recv_buffer_v = new double[recv_elem[i]];
            recv_buffer_w = new double[recv_elem[i]];
            recv_buffer_x = new double[recv_elem[i]];
            recv_buffer_y = new double[recv_elem[i]];
            recv_buffer_z = new double[recv_elem[i]];

            //              MPI_Recv(recv_buffer_u, recv_elem[i], MPI_DOUBLE, i,
            //              100, MPI_COMM_WORLD, &status[i]);
            //              MPI_Recv(recv_buffer_v, recv_elem[i], MPI_DOUBLE, i,
            //              101, MPI_COMM_WORLD, &status[i]);
            //              MPI_Recv(recv_buffer_w, recv_elem[i], MPI_DOUBLE, i,
            //              102, MPI_COMM_WORLD, &status[i]);
            //              MPI_Recv(recv_buffer_x, recv_elem[i], MPI_DOUBLE, i,
            //              103, MPI_COMM_WORLD, &status[i]);
            //              MPI_Recv(recv_buffer_y, recv_elem[i], MPI_DOUBLE, i,
            //              104, MPI_COMM_WORLD, &status[i]);
            //              MPI_Recv(recv_buffer_z, recv_elem[i], MPI_DOUBLE, i,
            //              105, MPI_COMM_WORLD, &status[i]);

            particle* temp = new particle[recv_elem[i]];

            for (int j = 0; j < recv_elem[i]; j++) {
              temp[j].velocity.x = recv_buffer_u[j];
              temp[j].velocity.y = recv_buffer_v[j];
              temp[j].velocity.z = recv_buffer_w[j];
              temp[j].velocity.x = recv_buffer_x[j];
              temp[j].velocity.y = recv_buffer_y[j];
              temp[j].velocity.z = recv_buffer_z[j];
            }
            cell_iter = cells.begin();
            while (cell_iter != cells.end()) {
              (*cell_iter)->add_particle(temp, recv_elem[i]);
              cell_iter++;
            }

            delete[] temp;
            delete[] recv_buffer_u;
            delete[] recv_buffer_v;
            delete[] recv_buffer_w;
            delete[] recv_buffer_x;
            delete[] recv_buffer_y;
            delete[] recv_buffer_z;

            step2[i] = true;
          }
        }

        for (int i = 0; i < numproc; i++) {
          if (i != proc_id) recv = recv && (!step2[i]);
        }
      }

      return true;
    } catch (bad_alloc& ba) {
      cout << ba.what() << endl;
      //        if(status != NULL) delete [] status;
      if (step1 != NULL) delete[] step1;
      if (step2 != NULL) delete[] step2;
      if (recv_elem != NULL) delete[] recv_elem;
      exit(unusual_situations::exit_code::MEMORY_NOT_ALLOCATED);
    }
  }
  return true;
}

bool cell_cluster::write_file(const char* file_name) {
  // if(proc_id != 0) return false;
  std::ofstream file1(file_name);
  cell_iter = cells.begin();
  point ap;
  unsigned int N;
  while (cell_iter != cells.end()) {
    ap = (*cell_iter)->get_center();
    N = (*cell_iter)->N();
    if (double(N) / ((*cell_iter)->get_volume() * density) < 0) {
      // cout<< "N = " << N << "\tvolume = " << (*cell_iter)->get_volume() <<
      // "\tdensity = " << density << endl;
      bool a = body->is_inner_point(ap);
      cout << a << endl;
      a = body->is_inner_point((*cell_iter)->get_apex());
      cout << a << endl;
    }
    file1 << ap << N << "\t"
          << double(N) / ((*cell_iter)->get_volume() * density) << "\t"
          << (*cell_iter)->get_t() << "\t" << (*cell_iter)->get_velocity()
          << endl;
    cell_iter++;
  }
  file1.close();
  return true;
}

bool cell_cluster::write_speed_file(const char* file_name) {
  // if(proc_id != 0) return false;
  std::ofstream file1(file_name);
  cell_iter = cells.begin();
  point ap;
  unsigned int N;
  while (cell_iter != cells.end()) {
    ap = (*cell_iter)->get_center();
    N = (*cell_iter)->N();
    if (double(N) / ((*cell_iter)->get_volume() * density) < 0) {
      // cout<< "N = " << N << "\tvolume = " << (*cell_iter)->get_volume() <<
      // "\tdensity = " << density << endl;
      bool a = body->is_inner_point(ap);
      cout << a << endl;
      a = body->is_inner_point((*cell_iter)->get_apex());
      cout << a << endl;
    }
    file1 << ap << (*cell_iter)->get_velocity() << endl;
    cell_iter++;
  }
  file1.close();
  return true;
}

bool cell_cluster::write_file() {
  if (proc_id == 0) {
    // density = 10000;
    std::ofstream file1("data.dat");
    cell_iter = cells.begin();
    point ap;
    unsigned int N;
    double av_den = 0;
    point cell_size;
    double volume = 1;
    while (cell_iter != cells.end()) {
      ap = (*cell_iter)->get_center();
      // ap = (*cell_iter)->get_mass_center();
      /*volume = (*cell_iter)->get_volume();
      cell_size = (*cell_iter)->get_size();
      if(volume <= 0) volume = cell_size.x * cell_size.y * cell_size.z;
      else volume = (*cell_iter)->get_volume() / cell_size.x * cell_size.y *
      cell_size.z;*/
      // ap = (*cell_iter)->get_particle_mass_center();
      N = (*cell_iter)->N();
      av_den += double(N) / np;
      file1 << ap << "\t" << double(N) / ((*cell_iter)->get_volume() * density)
            << "\t" << (*cell_iter)->get_t() << endl;
      cell_iter++;
    }

    /*for(int i = 1; i < numproc; i++)
    {
      MPI_Status status;
      int k;
      MPI_Probe(i, 10, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &k);
      double *data_x = new double[k];
      double *data_y = new double[k];
      double *data_z = new double[k];
      double *data_density = new double[k];
      double *data_T = new double[k];
      MPI_Recv(data_x, k, MPI_DOUBLE, i, 10, MPI_COMM_WORLD, &status);
      MPI_Recv(data_y, k, MPI_DOUBLE, i, 11, MPI_COMM_WORLD, &status);
      MPI_Recv(data_z, k, MPI_DOUBLE, i, 12, MPI_COMM_WORLD, &status);
      MPI_Recv(data_density, k, MPI_DOUBLE, i, 13, MPI_COMM_WORLD, &status);
      MPI_Recv(data_T, k, MPI_DOUBLE, i, 14, MPI_COMM_WORLD, &status);

      for(int j = 0; j < k; j++)
      {
        file1 << data_x[j] << "\t" << data_y[j] << "\t" << data_z[j] << "\t" <<
    data_density[j] << "\t" << data_T[j] << endl;;
      }
      delete [] data_x;
      delete [] data_y;
      delete [] data_z;
      delete [] data_density;
      delete [] data_T;
    }*/
    file1.close();
  } else {
    int n = cells.size();
    double* data_x = new double[n];
    double* data_y = new double[n];
    double* data_z = new double[n];
    double* data_density = new double[n];
    double* data_T = new double[n];

    point p;
    int i = 0;

    cell_iter = cells.begin();
    while (cell_iter != cells.end()) {
      p = (*cell_iter)->get_center();
      data_x[i] = p.x;
      data_y[i] = p.y;
      data_z[i] = p.z;
      data_density[i] = double((*cell_iter)->N()) / np;
      data_T[i] = (*cell_iter)->get_t();

      i++;
      cell_iter++;
    }
    //      MPI_Send(data_x, n, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
    //      MPI_Send(data_y, n, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD);
    //      MPI_Send(data_z, n, MPI_DOUBLE, 0, 12, MPI_COMM_WORLD);
    //      MPI_Send(data_density, n, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
    //      MPI_Send(data_T, n, MPI_DOUBLE, 0, 14, MPI_COMM_WORLD);
    delete[] data_x;
    delete[] data_y;
    delete[] data_z;
    delete[] data_density;
    delete[] data_T;
  }
  return true;
}

bool cell_cluster::write_speed_file() {
  if (proc_id == 0) {
    std::ofstream file1("speed.dat");
    cell_iter = cells.begin();
    point ap;
    unsigned int N;
    double av_den = 0;
    while (cell_iter != cells.end()) {
      ap = (*cell_iter)->get_center();
      N = (*cell_iter)->N();
      av_den += double(N) / np;
      file1 << ap << "\t" << (*cell_iter)->get_velocity() << endl;
      cell_iter++;
    }

    for (int i = 1; i < numproc; i++) {
      //        MPI_Status status;
      int k;
      //        MPI_Probe(i, 15, MPI_COMM_WORLD, &status);
      //        MPI_Get_count(&status, MPI_DOUBLE, &k);
      double* data_x = new double[k];
      double* data_y = new double[k];
      double* data_z = new double[k];
      double* data_u = new double[k];
      double* data_v = new double[k];
      double* data_w = new double[k];
      //        MPI_Recv(data_x, k, MPI_DOUBLE, i, 15, MPI_COMM_WORLD, &status);
      //        MPI_Recv(data_y, k, MPI_DOUBLE, i, 16, MPI_COMM_WORLD, &status);
      //        MPI_Recv(data_z, k, MPI_DOUBLE, i, 17, MPI_COMM_WORLD, &status);
      //        MPI_Recv(data_u, k, MPI_DOUBLE, i, 18, MPI_COMM_WORLD, &status);
      //        MPI_Recv(data_v, k, MPI_DOUBLE, i, 19, MPI_COMM_WORLD, &status);
      //        MPI_Recv(data_w, k, MPI_DOUBLE, i, 20, MPI_COMM_WORLD, &status);

      for (int j = 0; j < k; j++) {
        file1 << data_x[j] << "\t" << data_y[j] << "\t" << data_z[j] << "\t"
              << data_u[j] << "\t" << data_v[j] << "\t" << data_w[j] << endl;
      }
      delete[] data_x;
      delete[] data_y;
      delete[] data_z;
      delete[] data_u;
      delete[] data_v;
      delete[] data_w;
    }
    file1.close();
  } else {
    int n = cells.size();
    double* data_x = new double[n];
    double* data_y = new double[n];
    double* data_z = new double[n];
    double* data_u = new double[n];
    double* data_v = new double[n];
    double* data_w = new double[n];

    point p;
    int i = 0;

    cell_iter = cells.begin();
    while (cell_iter != cells.end()) {
      p = (*cell_iter)->get_center();
      data_x[i] = p.x;
      data_y[i] = p.y;
      data_z[i] = p.z;
      data_u[i] = (*cell_iter)->get_u();
      data_v[i] = (*cell_iter)->get_v();
      data_w[i] = (*cell_iter)->get_w();

      i++;
      cell_iter++;
    }
    //      MPI_Send(data_x, n, MPI_DOUBLE, 0, 15, MPI_COMM_WORLD);
    //      MPI_Send(data_y, n, MPI_DOUBLE, 0, 16, MPI_COMM_WORLD);
    //      MPI_Send(data_z, n, MPI_DOUBLE, 0, 17, MPI_COMM_WORLD);
    //      MPI_Send(data_u, n, MPI_DOUBLE, 0, 18, MPI_COMM_WORLD);
    //      MPI_Send(data_v, n, MPI_DOUBLE, 0, 19, MPI_COMM_WORLD);
    //      MPI_Send(data_w, n, MPI_DOUBLE, 0, 20, MPI_COMM_WORLD);
    delete[] data_x;
    delete[] data_y;
    delete[] data_z;
    delete[] data_u;
    delete[] data_v;
    delete[] data_w;
  }
  return true;
}

bool cell_cluster::write_times() {
  if (proc_id == 0) {
    std::ofstream times_file("time.dat");

    deque<double>::iterator time_iter = times.begin();

    int i = 0;

    while (time_iter != times.end()) {
      times_file << i << "\t" << times[i] << "\t" << sort_times[i] << "\t"
                 << send_times[i] << "\t" << bound_times[i] << "\t"
                 << calc_times[i] << endl;

      i++;
      time_iter++;
    }

    times_file.close();
  }
  return 0;
}

void cell_cluster::set_boundary_condition(boundary** a, int n) {
  for (int i = 0; i < n; i++) {
    boundary_cond_outer.push_back(a[i]);
  }

  //добавляем ссылки на ячейки в свободные границы
  deque<boundary*>::iterator bc = boundary_cond_outer.begin();
  while (bc != boundary_cond_outer.end()) {
    if (typeid(**bc) == typeid(free_boundary)) {
      free_boundary* b = dynamic_cast<free_boundary*>(*bc);
      b->add_cell(&cells);
      // b->set_np(np);
    } else if (typeid(**bc) == typeid(giper_free_boundary)) {
      free_boundary* b = dynamic_cast<free_boundary*>(*bc);
      b->add_cell(&cells);
      // b->set_np(np);
    }

    bc++;
  }
}

void cell_cluster::set_boundary_condition(boundary* a) {
  if (typeid(a) == typeid(free_boundary)) {
    free_boundary* b = dynamic_cast<free_boundary*>(a);
    b->add_cell(&cells);
  }

  boundary_cond_outer.push_back(a);
}

void cell_cluster::computation(void) {
  // unsigned int a = 0;
  while (t < t_end) {
    if (!proc_id) cout << "Time step:" << ++step << endl;
    this->time_step();
    /*if(!(a%5))
    {
      std::sort(cells.begin(), cells.end());
    }*/
  }
}

void cell_cluster::set_end_time(double t_end) {
  this->data_t = 0;
  this->data_dt = t_end / 100 + 0.000001;
  this->t_end = t_end;
}
void cell_cluster::calc_dt() {
  dt = 1000000.;

  cell_iter = cells.begin();
  // cout << "T = " << (*cell_iter)->calc_T() << endl;
  while (cell_iter != cells.end()) {
    double dtt = (*cell_iter)->calc_dt();
    if (dtt < dt) dt = dtt;
    cell_iter++;
  }
  dt *= Cu;
  if (t + dt > t_end) dt = t_end - t + 0.000000001;
}

void cell_cluster::set_dt_in_cells(double dt) {
  cell_iter = cells.begin();
  while (cell_iter != cells.end()) {
    (*cell_iter)->set_dt(dt);
    cell_iter++;
  }
}

bool cell_cluster::send_dt() {
  for (int i = 0; i < numproc; i++) {
    int snd = false;
    //      if(i != proc_id) snd = MPI_Send(&dt, 1, MPI_DOUBLE, i, 8,
    //      MPI_COMM_WORLD);
  }
  return true;
}

bool cell_cluster::recv_dt() {
  if (numproc > 1) {
    for (int i = 0; i < numproc; i++) {
      if (i != proc_id) {
        double dtt;
        //          MPI_Status status;
        int rcv = false;
        //          if(i != proc_id) MPI_Recv(&dtt, 1, MPI_DOUBLE, i, 8,
        //          MPI_COMM_WORLD, &status);
        if (dtt < dt) dt = dtt;
      }
    }
  }
  return true;
}

void cell_cluster::sync_dt() {
  calc_dt();
  send_dt();
  recv_dt();
  set_dt_in_cells(dt);
}

void cell_cluster::sync_data() {
  if (numproc > 0) {
    int n = partile_buffer.size();
    int* recv_el = new int[numproc];
    for (int i = 0; i < numproc; i++) {
      if (i == proc_id) continue;
      //        MPI_Send(&n, 1, MPI_INT, i, 40, MPI_COMM_WORLD);
    }
    for (int i = 0; i < numproc; i++) {
      if (i == proc_id) continue;
      //        MPI_Status status;
      //        MPI_Recv(&recv_el[i], 1, MPI_INT, i, 40, MPI_COMM_WORLD,
      //        &status);
    }

    send_buffer_u = new double[n];
    send_buffer_v = new double[n];
    send_buffer_w = new double[n];
    send_buffer_x = new double[n];
    send_buffer_y = new double[n];
    send_buffer_z = new double[n];

    int k = 0;

    deque<particle>::iterator particle_iter = partile_buffer.begin();
    while (particle_iter != partile_buffer.end()) {
      send_buffer_u[k] = particle_iter->get_u();
      send_buffer_v[k] = particle_iter->get_v();
      send_buffer_w[k] = particle_iter->get_w();
      send_buffer_x[k] = particle_iter->get_x();
      send_buffer_y[k] = particle_iter->get_y();
      send_buffer_z[k] = particle_iter->get_z();

      k++;

      particle_iter++;
    }

    partile_buffer.clear();

    for (int i = 0; i < numproc; i++) {
      if (i == proc_id) continue;
      //        MPI_Status status;
      recv_buffer_u = new double[recv_el[i]];
      recv_buffer_v = new double[recv_el[i]];
      recv_buffer_w = new double[recv_el[i]];
      recv_buffer_x = new double[recv_el[i]];
      recv_buffer_y = new double[recv_el[i]];
      recv_buffer_z = new double[recv_el[i]];
      //        MPI_Sendrecv(send_buffer_u, n, MPI_DOUBLE, i, 100,
      //        recv_buffer_u, recv_el[i], MPI_DOUBLE, i, 100, MPI_COMM_WORLD,
      //        &status); MPI_Sendrecv(send_buffer_v, n, MPI_DOUBLE, i, 101,
      //        recv_buffer_v, recv_el[i], MPI_DOUBLE, i, 101, MPI_COMM_WORLD,
      //        &status); MPI_Sendrecv(send_buffer_w, n, MPI_DOUBLE, i, 102,
      //        recv_buffer_w, recv_el[i], MPI_DOUBLE, i, 102, MPI_COMM_WORLD,
      //        &status); MPI_Sendrecv(send_buffer_x, n, MPI_DOUBLE, i, 103,
      //        recv_buffer_x, recv_el[i], MPI_DOUBLE, i, 103, MPI_COMM_WORLD,
      //        &status); MPI_Sendrecv(send_buffer_y, n, MPI_DOUBLE, i, 104,
      //        recv_buffer_y, recv_el[i], MPI_DOUBLE, i, 104, MPI_COMM_WORLD,
      //        &status); MPI_Sendrecv(send_buffer_z, n, MPI_DOUBLE, i, 105,
      //        recv_buffer_z, recv_el[i], MPI_DOUBLE, i, 105, MPI_COMM_WORLD,
      //        &status);

      for (int j = 0; j < recv_el[i]; j++) {
        particle a(recv_buffer_x[j], recv_buffer_y[j], recv_buffer_z[j],
                   recv_buffer_u[j], recv_buffer_v[j], recv_buffer_w[j]);
        /*a.u = recv_buffer_u[j];
        a.v = recv_buffer_v[j];
        a.w = recv_buffer_w[j];
        a.x = recv_buffer_x[j];
        a.y = recv_buffer_y[j];
        a.z = recv_buffer_z[j];*/
        partile_buffer.push_back(a);
      }

      cell_iter = cells.begin();

      while (cell_iter != cells.end()) {
        (*cell_iter)->add_particle(&partile_buffer);
        cell_iter++;
      }

      delete[] recv_buffer_u;
      delete[] recv_buffer_v;
      delete[] recv_buffer_w;
      delete[] recv_buffer_x;
      delete[] recv_buffer_y;
      delete[] recv_buffer_z;
    }

    partile_buffer.clear();

    delete[] send_buffer_u;
    delete[] send_buffer_v;
    delete[] send_buffer_w;
    delete[] send_buffer_x;
    delete[] send_buffer_y;
    delete[] send_buffer_z;

    delete[] recv_el;

    /*cout << "Sending data..." << flush;
    send_data();
    cout << "is over." << endl << "Reciving data...";
    recv_data();
    cout << "is over." << endl;*/
  }
}

void cell_cluster::set_geometry(geometry* body) { this->body = body; }

void cell_cluster::cells_fragmentation() {
  deque<cell*> temp_cells;

  cell_iter = cells.begin();
  while (cell_iter != cells.end()) {
    if ((*cell_iter)->get_body_mark()) {
      deque<cell*> new_cells = ((*cell_iter)->fragmentation(body));
      temp_cells.insert(temp_cells.begin(), new_cells.begin(), new_cells.end());
      delete (*cell_iter);
      *cell_iter = NULL;

      // cell_iter = cells.erase(cell_iter);
      // continue;
    }
    cell_iter++;
  }

  cells.clear();
  cells.insert(cells.end(), temp_cells.begin(), temp_cells.end());

  temp_cells.clear();
}

void cell_cluster::cells_test() {
  deque<cell*> cells_swap;

  cout << "Number of cells: " << cells.size() << endl
       << "|--------------------------------------|" << endl;
  int i = 0, j = 1;

  cell_iter = cells.begin();
  while (cell_iter != cells.end()) {
    i++;
    if (int(i * 41 / cells.size()) > j) {
      cout << "*" << flush;
      j++;
    }
    /*if(cell_iter == cells.end())
    {
      break;
    }*/
    /*if(*cell_iter == NULL)
    {
      cell_iter = cells.erase(cell_iter);
      continue;
    }
    else if((*cell_iter)->calc_volume() <= 0)
    {
      if(body != NULL)
      {
        if(body->is_inner_point((*cell_iter)->get_center()))
        {
          //delete *cell_iter;
          //*cell_iter = NULL;
          //cell_iter = cells.erase(cell_iter);
          //cell_iter = cells.begin();
          //continue;
        }
      }
    }*/

    (*cell_iter)->calc_volume();

    cell_iter++;
  }
}

bool cell_cluster::initialazition(const char* init_file, double Kn, double Cu,
                                  double L, geometry* body, int numproc,
                                  int my_id) {
  //    begin_time = MPI_Wtime();	                    //устанавливаем
  //    время начала работы
  begin_iter_time = begin_time;  //устанавливаем время начала нулевой
                                 //итерации(нулевая итерация - инициализация)

  this->numproc = numproc;
  this->proc_id = my_id;

  this->np = np;

  this->Kn = Kn;
  this->Cu = Cu;

  this->t = 0;  //устанавливаем "физическое" время на ноль
  this->dt = 100000;

  this->volume = 0;
  this->N = 0;

  this->set_geometry(body);

  double vol = 0;

  std::cout << endl
            << "...reading cell cluster file " << init_file << "... " << endl;
  string word;
  vector<string> v;
  std::ifstream file(init_file);
  if (!file.is_open()) {
    std::cout << endl << "Can't open file!!!" << endl;
    throw unusual_situations::exception(
        mc3d::unusual_situations::exit_code::BOUNDARY_CONSTRUCTOR_ERROR, this);
  }

  while (file >> word) {
    v.push_back(word);
  }

  if (v[0] == "cluster" && v[v.size() - 1] == "endcluster") {
    std::vector<string>::size_type i = 1;
    while (i < v.size()) {
      if (v[i] == "comment") {
        while (v[i] != "endcomment") {
          i++;
        }
      }
      if (v[i] == "cellgroup") {
        point apex_;
        point size;
        unsigned int ncx, ncy, ncz;
        point velocity;
        unsigned int N_;
        double temperature;
        double S;
        double alpha;
        double density_ = -1;
        i++;
        while (v[i] != "endcellgroup") {
          if (v[i] == "apex") {
            apex_.x = std::atof(v[i + 1].c_str());
            apex_.y = std::atof(v[i + 2].c_str());
            apex_.z = std::atof(v[i + 3].c_str());
            i += 4;
          }

          if (v[i] == "size") {
            size.x = std::atof(v[i + 1].c_str());
            size.y = std::atof(v[i + 2].c_str());
            size.z = std::atof(v[i + 3].c_str());
            i += 4;
          }

          if (v[i] == "nc") {
            ncx = (unsigned int)std::atoi(v[i + 1].c_str());
            ncy = (unsigned int)std::atoi(v[i + 2].c_str());
            ncz = (unsigned int)std::atoi(v[i + 3].c_str());
            i += 4;
          }

          if (v[i] == "speed_ratio") {
            S = std::atof(v[++i].c_str());
            i++;
          }

          if (v[i] == "alpha") {
            alpha = std::atof(v[++i].c_str());
            i++;
          }

          if (v[i] == "density") {
            density_ = (unsigned int)std::atoi(v[++i].c_str());
            i++;
          }

          if (v[i] == "particles") {
            N_ = (unsigned int)std::atoi(v[++i].c_str());
            i++;
          }

          if (v[i] == "temperature") {
            temperature = std::atof(v[++i].c_str());
            i++;
          }
        }

        double dx = size.x / double(ncx);  //устанавливаем размер ячейки
        double dy = size.y / double(ncy);  //устанавливаем размер ячейки
        double dz = size.z / double(ncz);  //устанавливаем размер ячейки

        if (density_ > 0) {
          N_ = density_ * dx * dy * dz;
        }

        this->ncx = ncx;
        this->ncy = ncy;
        this->ncz = ncz;

        for (int i = 0; i < int(ncx);
             i++)  //цыклы генерации ячеек	по 3м координатам
        {
          for (int j = 0; j < int(ncy); j++) {
            for (int k = 0; k < int(ncz); k++) {
              point a;
              a.x = apex.x + i * dx;
              a.y = apex.y + j * dy;
              a.z = apex.z + k * dz;

              cell* new_cell = new cell;
              new_cell->set_param(S, alpha, 1);

              new_cell->set_apex(a);  //задание опорной точки
              new_cell->set_size(dx, dy, dz);  //задание размеров

              new_cell->initialazition(
                  N_,
                  body);  //инициализация ячейки, с заданным количеством частиц

              new_cell->set_L(L);  //задание характерного размера для ячейки
              new_cell->set_Kn(Kn);  //задание кнудсена для ячейуи
              cells.push_back(new_cell);  //копирование и вставка подготовленной
                                          //ячейки в дэк класстера

              N += N_;

              vol += new_cell->get_size().volume();

              double dtt = new_cell->get_dt();
              dt = min(this->dt, dtt);
            }
          }
        }
      }

      if (v[i] == "cell") {
        point apex_;
        point size;
        point velocity;
        unsigned int N_;
        double temperature;
        double S;
        double alpha;
        double L;
        i++;
        while (v[i] != "endcell") {
          if (v[i] == "apex") {
            apex_.x = std::atof(v[i + 1].c_str());
            apex_.y = std::atof(v[i + 2].c_str());
            apex_.z = std::atof(v[i + 3].c_str());
            i += 4;
          }

          if (v[i] == "size") {
            size.x = std::atof(v[i + 1].c_str());
            size.y = std::atof(v[i + 2].c_str());
            size.z = std::atof(v[i + 3].c_str());
            i += 4;
          }

          if (v[i] == "S") {
            S = std::atof(v[i + 1].c_str());
            i++;
          }

          if (v[i] == "alpha") {
            alpha = std::atof(v[++i].c_str());
            i++;
          }

          if (v[i] == "particles") {
            N_ = (unsigned int)std::atoi(v[++i].c_str());
            i++;
          }

          if (v[i] == "temperature") {
            temperature = std::atof(v[++i].c_str());
            i++;
          }

          if (v[i] == "L") {
            L = std::atof(v[++i].c_str());
            i++;
          }
        }

        // S = sqrt((velocity.x

        cell* new_cell = new cell;

        new_cell->set_param(S, alpha, 1);

        new_cell->set_apex(apex_);  //задание опорной точки
        new_cell->set_size(size);   //задание размеров

        new_cell->initialazition(
            N_, body);  //инициализация ячейки, с заданным количеством частиц

        new_cell->set_L(L);  //задание характерного размера для ячейки
        new_cell->set_Kn(Kn);  //задание кнудсена для ячейуи
        cells.push_back(new_cell);  //копирование и вставка подготовленной
                                    //ячейки в дэк класстера

        N += N_;

        double dtt = new_cell->get_dt();
        dt = min(this->dt, dtt);
      }

      i++;
    }
  }

  this->density = double(N) / vol;

  this->cells_test();

  // this->clean_inner_particle();

  return true;
}

/*bool cell_cluster::read_cell_file(const char *init_file)
{
}*/

bool cell_cluster::write_cell_file(const char* init_file) {
  std::ofstream file(init_file);
  if (!file.is_open()) {
    std::cout << endl << "Can't open file!!!" << endl;
    return false;
  }

  file << "cluster" << std::endl;

  cell_iter = cells.begin();
  while (cell_iter != cells.end()) {
    point V = (*cell_iter)->get_velocity();
    double S = sqrt((V.x * V.x + V.y * V.y) / 2 * (*cell_iter)->get_t());
    double alpha = 3.141592654 / 2 - std::atan(V.x / V.y);
    if (V.x > 100000 * V.y) {
      V.x > 0 ? alpha = 0 : alpha = 3.141592654;
    }

    file << "\tcell" << endl;

    file << "\t\tapex\t" << (*cell_iter)->get_apex() << endl;

    file << "\t\tsize\t" << (*cell_iter)->get_size() << endl;

    file << "\t\tS\t" << S << endl;

    file << "\t\talpha\t" << alpha << endl;

    file << "\t\tparticles\t" << (*cell_iter)->N() << endl;

    file << "\t\ttemperature\t" << (*cell_iter)->get_t() << endl;

    file << "\t\tL\t" << (*cell_iter)->get_L() << endl;

    file << "\tendcell" << endl;

    cell_iter++;
  }

  return true;
}

void cell_cluster::set_data_save_dtime(double data_dt) {
  this->data_dt = data_dt;
}

void cell_cluster::clean_inner_particle() {
  cell_iter = cells.begin();
  while (cell_iter != cells.end()) {
    (*cell_iter)->clean_inner_particle(body);
    cell_iter++;
  }
}
}  // namespace mc3d
