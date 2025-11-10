//File: cell_claster.cpp
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.1
//Last modified: 15.05.09.
//Description: Program for calculation of freemolecular flows.

#include "cell_claster.h"
#include "Header.h"
#include <iostream>
#include <cstdlib>
#include <fstream>

#pragma once

using namespace std;

//double F(double x, double y, double z)
//{
//	if(x <= 0.5) return 1.2;
//	return 1;
//}

cell_claster::cell_claster(void)
{
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
}

cell_claster::~cell_claster(void)
{
	delete [] cells;
}

bool cell_claster::initialazition(unsigned int ncx, unsigned int ncy, unsigned int ncz, unsigned int np)
{
	this->ncx = ncx;
	this->ncy = ncy;
	this->ncz = ncz;

	this->np = np;
	
	unsigned int N = 0;

	double dx = Lx / double(ncx);
	double dy = Ly / double(ncy);
	double dz = Lz / double(ncz);

	cells = new cell[ncx * ncy * ncz];

	this->t = 0;
	this->dt = 100000;

	for(unsigned int i = 0; i < ncx; i++)
	{
		for(unsigned int j = 0; j < ncy; j++)
		{
			for(unsigned int k = 0; k < ncz; k++)
			{
				unsigned int i_j_k = k * ncx * ncy + j * ncx + i;

				unsigned int N_;
				point a;
				a.x = apex.x + i * dx;
				a.y = apex.y + j * dy;
				a.z = apex.z + k * dz;
				N_ = unsigned int(np * F(a.x + 0.5 * dx, a.y + 0.5 * dy, a.z + 0.5 * dz));

				cells[i_j_k].set_apex(a);
				cells[i_j_k].set_size(dx, dy, dz);

				if(!cells[i_j_k].initialazition(N_))
				{
					cout<<"Error, bad cell initialazition."<<endl;
					exit(1);
				}

				cells[i_j_k].generate_random();

				N += N_;
				
				double dtt = cells[i_j_k].get_dt();
				this->dt = min(this->dt, dtt); 
			}
		}
	}

	this->first_cell();

	while(true)
	{
		calc_cell->set_dt(this->dt);
		if(!go_to_next_cell()) break;
	}

	this->N = N;

	return 0;
}

//double cell_claster::min_Kn()
//{
//	double minKn = 1000;
//	
//	for(unsigned int i = 0; i < ncz; i++)
//	{
//		for(unsigned int j = 0; j < ncy; j++)
//		{
//			for(unsigned int k = 0; k < ncx; k++)
//			{
//				unsigned int i_j_k = i * ncz * ncy + j * ncz + k;
//
//				double Kn_l = Kn * double(np) / double(cells[i_j_k].N());
//				if(minKn > Kn_l) minKn = Kn_l;
//			}
//		}
//	}
//	return minKn;
//}

void cell_claster::set_apex(point apex)
{
	this->apex = apex;
}

void cell_claster::set_size(double Lx, double Ly, double Lz)
{
	this->Lx = Lx;
	this->Ly = Ly;
	this->Lz = Lz;
}

void cell_claster::time_step()
{
	calc_cell = first_cell();
	unsigned int a = 0;
	while(true)
	{
		//cout<<"Cell number:"<<a++<<endl;
		calc_cell->time_step();
		//cout<<"Particle moveing is over"<<endl;
		partile_buffer.add_particle_array(calc_cell->get_buffer());
		collisions();
		//cout<<"Collisions are over"<<endl;
		boundary_condition();
		//cout<<"Boundary condition is done. Go to the next cell."<<endl;
		if(!go_to_next_cell()) break;
	}
	calc_cell = first_cell();
	int j = 0;
	while(true)
	{
		partile_buffer.calc_N();
		calc_cell->add_particle(&partile_buffer);
		if(!go_to_next_cell()) break;
	}
	send_data();
	t += dt;
}

void cell_claster::collisions()
{
	t = 0;
	dt = 0;
}

void cell_claster::boundary_condition()
{
	for(int i = 0; i < 6; i++)
	{
		boundary_cond[i]->bondary_condition(&partile_buffer);
	}
}

bool cell_claster::go_to_next_cell()
{
	if(++i_j_k >= ncx * ncy *ncz)
	{
		i_j_k--;
		return false;
	}
	calc_cell = &cells[i_j_k];
	return true;
}

cell* cell_claster::first_cell()
{
	i_j_k = 0;
	calc_cell = &cells[i_j_k];
	return calc_cell;
}

bool cell_claster::send_data()
{
	return true;
}

bool cell_claster::write_file()
{
	std::ofstream file("data.dat");
	this->first_cell();
	while(go_to_next_cell())
	{
		file<<calc_cell<<endl;
	}
	return true;
}

void cell_claster::set_boundary_condition(boundary *a, int i)
{
	if((i >= 0) && (i < 6))	boundary_cond[i] = a;
}

void cell_claster::computation(void)
{
	unsigned int a = 0;
	while(t < t_end)
	{
		this->time_step();
		std::cout<<"Time step:"<<a++<<endl;
	}
	//this->write_file();
}

void cell_claster::set_end_time(double t_end)
{
	this->t_end = t_end;
}