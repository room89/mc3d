// File: cell.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.5.1
// Last modified: 18.05.10.
// Description: Program for calculation of rarefaid flows.

#pragma once

#include <cstdlib>
#include <deque>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <vector>
// #include "free_boundary.h"
#include "geometry.h"
#include "inner_boundary.h"
#include "particle.h"
#include "point.h"
// #include <mpi.h>
// #include "free_boundary.h
// #include "Header.h"

using namespace std;
using namespace mc3d;

namespace mc3d {
class Cell {
 private:
  std::vector<Particle> particles;
  std::vector<Particle>
      particle_buffer;  // Particle that departing from the cell
  Point apex;
  double lx, ly, lz;  // size of cell
  double T;           // temperature in cell
  Point velocity;
  Point mass_center;
  double volume_;
  double Kn_l;  // локальый кнутсен
  double Kn;
  double L;   // характерный размер
  double dt;  //
  double time;
  double calc_time;
  double start_time;
  bool body_mark;  // is body inside cell
  unsigned int np;
  InnerBoundary body_boundary;
  std::vector<std::reference_wrapper<Cell>> neighbors;
  std::optional<std::reference_wrapper<int>> thread_mark;
  std::shared_ptr<std::mutex> particles_mutex_;

 public:
  Cell();
  virtual ~Cell();

  virtual void MoveParticles();  // перемещение частиц
  virtual void Collisions();     // соударения между частицами
  void SortNeighbors();
  unsigned int GetParticleCount() const;  // возвращает количество частиц
  std::size_t CountInnerParticles(const Geometry& body) const;
  double GetTemperatureRaw();
  void SetSize(double lx, double ly, double lz);  // установка размера ячейки
  void SetSize(Point dl);
  void SetApex(Point a);  // установка "опорной" точки
  void SetCharacteristicLength(double length);
  void SetTemperature(double t);
  void SetParameters(double S, double alpha, double T);
  void SetVelocity(Point velocity);
  void SetKn(double Kn);
  void SetInnerBoundary(InnerBoundary bound);
  void SetBodyMark(bool mark);
  virtual bool Initialize(
      size_t N,
      const std::unique_ptr<Geometry>&
          body);  // инициализация ячейки (N-количество частиц)
  virtual double GenerateRandom(size_t N);
  virtual double GenerateRandom(unsigned int N, double T, Point V);
  virtual double GenerateFreeRandom(unsigned int N, double T, Point V,
                                    Point normal);
  virtual double GenerateHyperFreeRandom(unsigned int N, Point V, double T);
  void Sort();
  std::vector<Particle>& GetBuffer();
  void AddParticle(std::vector<Particle>& particles);
  bool TryAcceptParticle(Particle& particle);
  void FixParticle(unsigned int N, double T, Point V);
  double CalculateKn();
  Point GetApex() const;
  Point GetCenter() const;
  Point GetMassCenter() const;
  Point GetSize() const;
  double GetKn();
  double GetDt();
  double PeekDt() const;
  double GetTemperature();
  double GetU();
  double GetV();
  double GetW();
  double GetVolume() const;
  Point GetParticleMassCenter();
  Point GetVelocity();
  double GetEnergy();
  double GetCharacteristicLength();
  bool GetBodyMark();
  Point CalculateVelocity() const;
  double CalculateTemperature();
  double CalculateVolume();
  void SetDt(double dt);
  double CalculateDt();
  void AddNeighbor(Cell& neighbor);
  friend ostream& operator<<(ostream& o, const Cell& c);
  bool WriteFile(ofstream* file);
  bool WriteFile();
  void Calculate();
  void AttachThreadMark(int& ptr);
  std::deque<std::unique_ptr<Cell>> Fragment(
      const std::unique_ptr<Geometry>& body);
  void CleanInnerParticles(const Geometry& body);

  bool DebugTestParticle();

  friend bool operator>(const Cell& a, const Cell& b);
  friend bool operator<(const Cell& a, const Cell& b);
};
}  // namespace mc3d
