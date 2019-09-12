// File: main.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.5.1
// Last modified: 18.05.10.
// Description: Program for calculation of rarefaid flows.

#include "cell.h"
#include "cell_cluster.h"
#include "point.h"

using namespace std;

const double Lx = 1, Ly = 1, Lz = 1;

int main(int argc, char* argv[]) {
  auto cone = std::make_unique<geometry>();
  cone->create_wedge(-0.2, 1.1, 0.3, 3.1415 / 12);
  //  cone->create_cube();
  cone->fix_poligons();
  //  cone->rev_nrml();
  cone->write_geometry_file("wadge.stl", "wadge");
  // cone->scaling(-1);
  cone->fragmentation(0.3);

  double Kn = .05, Cu = 0.9, T = 1;
  unsigned int ncx = 10, ncy = 10, ncz = 3;
  unsigned int np = 50;
  double S = 10;
  point a(-0.5 * Lx, -0.5 * Ly, -0.5 * Lz);

  mc3d::CellCluster cluster;

  cluster.set_apex(a);

  cluster.set_size(Lx, Ly, Lz);

  cluster.initialazition(ncx, ncy, ncz, np * ncx * ncy * ncz, Kn, Cu,
                         std::move(cone), S, 0, 1);

  cluster.WriteFile();

  boundary* boundary_cond[6];

  a.x = -Lx / 2;
  a.y = 0;
  a.z = 0;
  point b(-1, 0, 0);
  point c(Ly, 0, 0);

  // boundary_cond[0] = new mirror_boundary(a, b);
  // boundary_cond[1] = new mirror_boundary(a * -1, b * -1);

  // boundary_cond[0] = new pereodic_boundary(a, b, c);
  // boundary_cond[1] = new pereodic_boundary(a * -1, b * -1, c * -1);

  //  boundary_cond[0] = new FreeBoundary(a, b, np, S, T);
  //  boundary_cond[1] = new FreeBoundary(a * -1, b * -1, np, S, T);

  boundary_cond[0] = new giper_free_boundary(a, b, np, S, T);
  boundary_cond[1] = new giper_free_boundary(a * -1, b * -1, np, S, T);

  a.set(0, -Ly * 0.5, 0);
  b.set(0, -1, 0);
  c.set(0, Ly, 0);

  //  boundary_cond[2] = new mirror_boundary(a, b);
  //  boundary_cond[3] = new mirror_boundary(a * -1, b * -1);

  // boundary_cond[2] = new pereodic_boundary(a, b, c);
  // boundary_cond[3] = new pereodic_boundary(a * -1, b * -1, c * -1);

  //  boundary_cond[2] = new FreeBoundary(a, b, np, 0, T, 0);
  //  boundary_cond[3] = new FreeBoundary(a * -1, b * -1, np, 0, T, 0);

  boundary_cond[2] = new giper_free_boundary(a, b, np, S, T);
  boundary_cond[3] = new giper_free_boundary(a * -1, b * -1, np, S, T);

  a.set(0, 0, -Lz * 0.5);
  b.set(0, 0, -1);
  c.set(0, 0, Lz);

  // boundary_cond[4] = new mirror_boundary(a, b);
  // boundary_cond[5] = new mirror_boundary(a * -1, b * -1);

  //  boundary_cond[4] = new pereodic_boundary(a, b, c);
  //  boundary_cond[5] = new pereodic_boundary(a * -1, b * -1, c * -1);

  //  boundary_cond[4] = new FreeBoundary(a, b, np, 0, T, 0);
  //  boundary_cond[5] = new FreeBoundary(a * -1, b * -1, np, 0, T, 0);

  boundary_cond[4] = new giper_free_boundary(a, b, np, S, T);
  boundary_cond[5] = new giper_free_boundary(a * -1, b * -1, np, S, T);

  cluster.set_boundary_condition(boundary_cond, 6);
  // cluster.set_boundary_condition(boundary_cond[5]);

  cluster.write_speed_file();

  cluster.set_end_time(T);

  cluster.WriteFile("data_NU.dat");

  cluster.computation();

  // while(!cluster.time_step());
  std::cout << "Computation is over." << endl;
  cluster.write_cell_file("end_cell.net");
  cluster.WriteFile();
  cluster.write_times();
  cluster.write_speed_file();

  for (int i = 0; i < 6; i++) delete boundary_cond[i];

  // delete cone;

  //	work_time = MPI_Wtime() - work_time;

  //	MPI_Barrier(MPI_COMM_WORLD);

  //	MPI_Finalize();

  //	std::cin >> N;

  return 0;
}
