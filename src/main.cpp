//File: main.h 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.9
//Last modified: 18.05.10.
//Description: Program for calculation of freemolecular flows.

#pragma once
#include <mpi.h>
#include <iostream>
#include "timer.h"
#include "cell.h"
#include "point.h"
#include "cell_claster.h"
#include "pereodic_boundary.h"


//#include "header.h"
using namespace std;

const double Lx = 1,  Ly = 1, Lz = 1;
typedef unsigned int uint;

int main(int argc, char * argv[])
{
	int *argc1 = 0;
	char ***argv1 = NULL;
	MPI_Init(argc1, argv1);

	double Kn, Cu, T;
	unsigned int N = 0;
	unsigned int ncx = 0, ncy = 0, ncz = 0;
	//double dx, dy, dz;
	unsigned int np = 100;
	int numproc, proc_id;

	MPI_Comm_size(MPI_COMM_WORLD, &numproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);


	Kn = .05;
	cout<<"Input Knudsen number: "<< Kn <<endl;

	Cu = .9;
	cout<<"Input Courant number, Cu<=1: "<< Cu <<endl;

	T = .2;
	cout<<"Input time: "<< T;
	
	ncx = 25;
	cout<<"Input dimension, x: "<< ncx <<endl;
	ncy = 25;
	cout<<"Input dimension, y: "<< ncy <<endl;
	ncz = 1;
	cout<<"Input dimension, z: "<< ncz <<endl;
	
	np = 400;
	cout<<"Input statistics level: "<< np <<endl;

	cout<<"Numproc: "<<numproc<<"   "<<"proc_id: "<<proc_id<<endl;
	

	point a;
	a.x = -Lx * 0.5 + Lx * proc_id / numproc;
	a.y = -Ly * 0.5;
	a.z = -Lz * 0.5;

	cell_claster claster;

	claster.set_apex(a);
	claster.set_size(Lx / numproc, Ly, Lz);

	claster.initialazition(ncx / numproc, ncy, ncz, np, Kn, Cu, numproc, proc_id);
	
	boundary *boundary_cond[6];

	a.x = -Lx / 2;
	a.y = 0;
	a.z = 0;
	point b(-1,0,0);
	point c(Ly,0,0);

	boundary_cond[0] = new pereodic_boundary(a, b, c);
	boundary_cond[1] = new pereodic_boundary(a * -1, b * -1, c * -1);

	a.set(0, -Ly * 0.5, 0);
	b.set(0, -1, 0);
	c.set(0, Ly, 0);

	boundary_cond[2] = new pereodic_boundary(a, b, c);
	boundary_cond[3] = new pereodic_boundary(a * -1, b * -1, c * -1);

	a.set(0, 0, -Lz * 0.5);
	b.set(0, 0, -1);
	c.set(0, 0, Lz);

	boundary_cond[4] = new pereodic_boundary(a, b, c);
	boundary_cond[5] = new pereodic_boundary(a * -1, b * -1, c * -1);

	claster.set_boundary_condition(boundary_cond, 6);

	claster.set_end_time(T);

	cout<<" t1 = "<<global_timer.calc_av()<<endl;

	claster.computation();
	claster.write_file();

	for(int i = 0; i < 6; i++)
		delete boundary_cond[i];

	MPI_Finalize();
	//cin>>np;

	return 0;
}