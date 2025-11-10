//File: cell_claster.h 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.1
//Last modified: 16.01.09.
//Description: Program for calculation of freemolecular flows.

#pragma once
#include "cell.h"
#include "point.h"
#include "boundary.h"

class cell_claster
{
private:
	cell *cells;
	cell *calc_cell;
	particle_array partile_buffer;
	unsigned int i_j_k;
	unsigned int np;
	unsigned int ncx, ncy, ncz;
	unsigned int N;
	double Lx, Ly, Lz;
	point apex;
	double Kn;
	double t;
	double t_end;
	double dt;
	void collisions(void);										//соудрения
	cell* first_cell(void);
	void boundary_condition(void);
	boundary *boundary_cond[6];
	bool go_to_next_cell(void);
	bool send_data(void);
public:
	void time_step(void);
	cell_claster(void);
	~cell_claster(void);
	void set_apex(point apex);
	void set_size(double Lx, double Ly, double Lz);
	void set_end_time(double t_end);
	double min_Kn(void);
	bool initialazition(unsigned int ncx, unsigned int ncy, unsigned int ncz, unsigned int np);
	bool write_file();
	void set_boundary_condition(boundary *a, int i);
	void computation(void);
};
