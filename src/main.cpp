// File: main.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.5.1
// Last modified: 18.05.10.
// Description: Program for calculation of rarefaid flows.

#include "cell_cluster.h"
#include "point.h"

using namespace std;

const double Lx = 1, Ly = 1, Lz = 1;

int main(int argc, char* argv[]) {
  auto cone = std::make_unique<Geometry>();
  cone->CreateWedge(-0.2, 1.1, 0.3, 3.1415 / 12);
  //  cone->CreateCube();
  cone->FixPolygons();
  //  cone->ReverseNormals();
  cone->WriteGeometryFile("wadge.stl", "wadge");
  // cone->Scale(-1);
  cone->Fragment(0.3);

  double Kn = .05, Cu = 0.9, T = 1;
  unsigned int ncx = 10, ncy = 10, ncz = 3;
  unsigned int np = 50;
  double S = 10;
  Point a(-0.5 * Lx, -0.5 * Ly, -0.5 * Lz);

  mc3d::CellCluster cluster;

  cluster.SetApex(a);

  cluster.SetSize(Lx, Ly, Lz);

  cluster.Initialize(ncx, ncy, ncz, np * ncx * ncy * ncz, Kn, Cu,
                     std::move(cone), S, 0, 1);

  cluster.WriteFile();

  Boundary* boundary_cond[6];

  a.x = -Lx / 2;
  a.y = 0;
  a.z = 0;
  Point b(-1, 0, 0);
  Point c(Ly, 0, 0);

  // boundary_cond[0] = new MirrorBoundary(a, b);
  // boundary_cond[1] = new MirrorBoundary(a * -1, b * -1);

  // boundary_cond[0] = new PeriodicBoundary(a, b, c);
  // boundary_cond[1] = new PeriodicBoundary(a * -1, b * -1, c * -1);

  //  boundary_cond[0] = new FreeBoundary(a, b, np, S, T);
  //  boundary_cond[1] = new FreeBoundary(a * -1, b * -1, np, S, T);

  boundary_cond[0] = new HyperFreeBoundary(a, b, np, S, T);
  boundary_cond[1] = new HyperFreeBoundary(a * -1, b * -1, np, S, T);

  a.Set(0, -Ly * 0.5, 0);
  b.Set(0, -1, 0);
  c.Set(0, Ly, 0);

  //  boundary_cond[2] = new MirrorBoundary(a, b);
  //  boundary_cond[3] = new MirrorBoundary(a * -1, b * -1);

  // boundary_cond[2] = new PeriodicBoundary(a, b, c);
  // boundary_cond[3] = new PeriodicBoundary(a * -1, b * -1, c * -1);

  //  boundary_cond[2] = new FreeBoundary(a, b, np, 0, T, 0);
  //  boundary_cond[3] = new FreeBoundary(a * -1, b * -1, np, 0, T, 0);

  boundary_cond[2] = new HyperFreeBoundary(a, b, np, S, T);
  boundary_cond[3] = new HyperFreeBoundary(a * -1, b * -1, np, S, T);

  a.Set(0, 0, -Lz * 0.5);
  b.Set(0, 0, -1);
  c.Set(0, 0, Lz);

  // boundary_cond[4] = new MirrorBoundary(a, b);
  // boundary_cond[5] = new MirrorBoundary(a * -1, b * -1);

  //  boundary_cond[4] = new PeriodicBoundary(a, b, c);
  //  boundary_cond[5] = new PeriodicBoundary(a * -1, b * -1, c * -1);

  //  boundary_cond[4] = new FreeBoundary(a, b, np, 0, T, 0);
  //  boundary_cond[5] = new FreeBoundary(a * -1, b * -1, np, 0, T, 0);

  boundary_cond[4] = new HyperFreeBoundary(a, b, np, S, T);
  boundary_cond[5] = new HyperFreeBoundary(a * -1, b * -1, np, S, T);

  cluster.SetBoundaryCondition(boundary_cond, 6);
  // cluster.SetBoundaryCondition(boundary_cond[5]);

  cluster.WriteSpeedFile();

  cluster.SetEndTime(T);

  cluster.WriteFile("data_NU.dat");

  cluster.Compute();

  // while(!cluster.TimeStep());
  std::cout << "Computation is over." << endl;
  cluster.WriteCellFile("end_cell.net");
  cluster.WriteFile();
  cluster.WriteTimes();
  cluster.WriteSpeedFile();

  for (int i = 0; i < 6; i++) delete boundary_cond[i];

  // delete cone;

  //	work_time = MPI_Wtime() - work_time;

  //	MPI_Barrier(MPI_COMM_WORLD);

  //	MPI_Finalize();

  //	std::cin >> N;

  return 0;
}
