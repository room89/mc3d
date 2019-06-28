//File: cell_claster.h 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.9
//Last modified: 18.05.10.
//Description: Program for calculation of freemolecular flows.

#pragma once
#include "mpi.h"
#include "cell.h"
#include "point.h"
#include "boundary_ref.h"

class cell_claster
{
private:
	int numproc;
	int proc_id;
	deque<cell> cells;
	deque<cell>::iterator cell_iter;
	deque<particle> partile_buffer;
	unsigned int i_j_k;
	unsigned int np;
	unsigned int ncx, ncy, ncz;
	unsigned int N;
	double Lx, Ly, Lz;
	point apex;
	double Kn;
	double Cu;
	double t;
	double t_end;
	double dt;
	void collisions(void);										//соудрения
	void boundary_condition(void);
	deque<boundary_ref> boundary_cond;
	bool send_data();
	bool recv_data();
	bool send_dt();
	bool recv_dt();
public:
	void calc_dt();
	void set_dt_in_cells();
	void time_step(void);
	cell_claster(void);
	~cell_claster(void);
	void set_apex(point apex);
	void set_size(double Lx, double Ly, double Lz);
	void set_end_time(double t_end);
	double min_Kn(void);
	bool initialazition(unsigned int ncx, unsigned int ncy, unsigned int ncz, unsigned int np, double Kn, double Cu, int numproc, int my_id);
	bool write_file();
	void set_boundary_condition(boundary **a, int i);
	void set_boundary_condition(boundary *a);
	void computation(void);
};
