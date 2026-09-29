#pragma once

#include <cell.h>
#include <exit_code.h>
#include <free_boundary.h>
#include <geometry.h>
#include <giper_free_boundary.h>
#include <inner_cell.h>
#include <mirror_boundary.h>
#include <pereodic_boundary.h>
#include <point.h>

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <limits>
#include <memory>
#include <utils/threadpool.hpp>
#include <vector>

namespace mc3d {

class CellClusterTestAccess;

class CellCluster {
 public:
  inline static constexpr std::size_t kDefaultThreadPoolSize = 4;
  friend class CellClusterTestAccess;
  inline void CalculateDt();  // вычисление шага по времени для данного кластера
  inline void SetDtInCells(
      double dt);  // установка шага по времени во все ячейки
  bool TimeStep();
  explicit CellCluster(std::size_t thread_pool_size =
                           kDefaultThreadPoolSize);  // конструктор по умолчанию
  ~CellCluster();                                    // деструктор
  void SetApex(Point apex);  // задание опорной точки кластера
  void SetSize(double Lx, double Ly, double Lz);  // задание размеров кластера
  void SetEndTime(double t_end);  // задание времени окончания рассчета
  void SetDataSaveDtime(double data_dt);
  // инициализация кластера(параметры:	ncx, ncy, ncz	- количество ячеек по
  // координатам, np - уровень статистик, Kn - Кнудсен
  bool Initialize(unsigned int ncx, unsigned int ncy, unsigned int ncz,
                  double density, double Kn, double Cu,
                  std::unique_ptr<Geometry>&& body_, double S, double alpha,
                  double T, double wall_temperature = 1.0);
  bool WriteCellFile(const string& init_file);
  void SetBinaryOutput(bool enabled);
  bool WriteFile();       // запись данных в файл(плотность и энэргия), данные
                          // собираються со всех кластеров
  bool WriteSpeedFile();  // запись данных в файл(скорость), данные
                          // собираются со всех кластеров
  bool WriteFile(
      const string& file_name);  // запись данных в файл(плотность и
                                 // энэргия) с заданным именем, данные
                                 // собираются со всех кластеров
  bool WriteSpeedFile(const char* file_name);
  void SetSnapshotInterval(double interval);
  bool WriteTimes();  // запись данных о времени выполнения каждой итерации
  void SetBoundaryCondition(
      std::vector<std::unique_ptr<Boundary>>&& boundaries);
  void SetBoundaryCondition(std::unique_ptr<Boundary> boundary);
  void Compute();  // вычисления
  void FragmentCells();
  void TestCells();
  void CleanInnerParticles();

 private:
  int numproc;  // количество процессов(не путать с нитями), специально для MPI
  int proc_id;  // номер процесса
  unsigned int i_j_k;  // переменная используется при инициализации кластера|
                       //  todo: удалить
  unsigned int np;     // уровень статистики
  unsigned int ncx, ncy, ncz;  // начальное количество ячеек по x, y, z
  unsigned int N;  // начальное количество частиц в кластере(в ячейках кластера)
  unsigned int step;
  double Lx, Ly, Lz;  // начальные размеры кластера
  Point apex;         // координаты опорной точки кластера.
  double Kn;          // Кнудсен
  double Cu;          // число Куранта
  double t;           // время с начала расчёта(физическое)
  // No end-time limit until SetEndTime is called (including during Initialize).
  double t_end = std::numeric_limits<double>::infinity();
  double dt;          // шаг по времени(должен быть равен во всех кластерах)
  double begin_time;  // время начала расчета(для вычисления времени работы)
  double begin_iter_time;  // время начала выполнения шага по времени
  deque<double> times;     // дэк с временами выполнения шага по времени
  deque<double>
      send_times;  // дэк с временами синхронизации данных между кластерами
  deque<double> bound_times;  // дэк с временами выполнения граничных условий
  deque<double> sort_times;   // дэк с временами сортировки частиц по ячейкам
  deque<double>
      calc_times;     // дэк с временами расчета соударений между частицами
  deque<Cell> cells;  // двусвязный список ячеек кластера.
  deque<Cell>::iterator cell_iter;       // итератор для обхода ячеек
  std::vector<Particle> partile_buffer;  // буфер кластера
  inline void BoundaryCondition();       // выполнение граничных условий
  std::vector<std::unique_ptr<Boundary>>
      boundary_cond_outer;  // внешние граничные условия
  inline bool
  SendData();  // посылка данных на все класстеры за исключением кластера из
               // которого посылают. не используется, возможно стоит удалить
  inline bool RecvData();  // прием данных с класстеров. не используеться,
                           // возможно стоит удалить
  inline void SyncData();  // синхронизация данных между кластерами. замена для
                           // методов:  send_data(), recv_data()
  inline void SyncDt();    // синхронизация шага по времени. тоже самое что и
                         // методы: send_dt() и recv_dt(), только в одном методе
  ofstream log;  // переменная для вывода логов, пока не реализовано
  void DistributeParticles(std::vector<Particle>& buffer);
  bool FindCellIndex(const Point& position, size_t& cell_index) const;
  std::unique_ptr<Geometry> body_;
  double data_dt;
  double data_t;
  double density;
  double volume;
  bool binary_output_ = false;
  double snapshot_interval_ = 0.0;
  double next_snapshot_time_ = 0.0;

  utils::ThreadPool thread_pool_;
  std::vector<Cell*> cell_lookup_;
  double cell_dx_ = 0.0;
  double cell_dy_ = 0.0;
  double cell_dz_ = 0.0;
  bool UpdateSnapshotSchedule();
  bool WriteTextSnapshot(std::ofstream& file);
  bool WriteBinarySnapshot(std::ofstream& file);
};
};  // namespace mc3d
