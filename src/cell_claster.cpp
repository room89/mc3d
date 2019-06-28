//File: cell_claster.cpp
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.9
//Last modified: 18.05.10.
//Description: Program for calculation of freemolecular flows.

#pragma once

#include "cell_claster.h"
#include "Header.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "point.h"
#include <complex>


const double Pi = 3.14159265358979;


using namespace std;

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
}

bool cell_claster::initialazition(unsigned int ncx, unsigned int ncy, unsigned int ncz, unsigned int np, double Kn, double Cu, int numproc, int proc_id)
{	
	this->numproc = numproc;
	this->proc_id = proc_id;

	this->ncx = ncx;
	this->ncy = ncy;
	this->ncz = ncz;

	this->np = np;
	
	this->Kn = Kn;
	this->Cu = Cu;

	unsigned int N = 0;

	double dx = Lx / double(ncx);
	double dy = Ly / double(ncy);
	double dz = Lz / double(ncz);

	this->t = 0;
	this->dt = 100000;

	cell_iter = cells.begin();
	cell temp;

	for(int i = 0; i < int(ncx); i++)
	{
		for(int j = 0; j < int(ncy); j++)
		{
			for(int k = 0; k < int(ncz); k++)
			{
				unsigned int i_j_k = k * ncx * ncy + j * ncx + i;

				unsigned int N_;
				point a(0, 0, 0);
				a.x = apex.x + i * dx;
				a.y = apex.y + j * dy;
				a.z = apex.z + k * dz;
				//N_ = unsigned int(np * F(a.x + 0.5 * dx, a.y + 0.5 * dy, a.z + 0.5 * dz));

				N_ = unsigned int(np * (1. + 0.1 * sin(2. * Pi * (a.x + 0.5*dx + 0.5)) * sin(2. * Pi * (a.y + 0.5*dy + 0.5))));
				

				temp.set_apex(a);
				temp.set_size(dx, dy, dz);
				temp.set_vel(0.5 * sin(2. * Pi * (a.x + 0.5*dx + 0.5)) * cos(2. * Pi * (a.y + 0.5*dy + 0.5)),
					-0.4 * cos(2. * Pi * (a.x + 0.5*dx + 0.5)) * sin(2. * Pi * (a.y + 0.5*dy + 0.5)), 0);
				temp.set_t(1. + 0.1 * cos(2. * Pi * (a.x + 0.5*dx + 0.5)) * cos(2. * Pi * (a.y + 0.5*dy + 0.5)));

				cell_iter = cells.insert(cell_iter, temp);

				cell_iter->initialazition(N_);

				N += N_;
				
				double dtt = cell_iter->get_dt();
				dt = min(this->dt, dtt); 
			}
		}
	}

	cell_iter = cells.begin();

	while(cell_iter != cells.end())
	{
		cell_iter->set_dt(this->dt);
		cell_iter++;
	}

	

	this->N = N;

	//this->write_file("data.dat");

	return 0;
}

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
	calc_dt();
	send_dt();
	recv_dt();
	set_dt_in_cells();
	
	cout<<"dt(cluster) = "<<dt<<endl;

	cell_iter = cells.begin();

	
	unsigned int a = 0;

	cout<<"cell[0].size() = "<<cell_iter->N()<<endl;

	deque<particle> *cell_buffer;

	while(cell_iter != cells.end())
	{
		//cout<<"Cell number:"<<a++<<endl;
		cell_iter->time_step();
		cell_buffer = cell_iter->get_buffer();
		partile_buffer.insert(partile_buffer.end(), cell_buffer->begin(), cell_buffer->end());
		cell_buffer->clear();
		collisions();
		//cout<<"Collisions are over"<<endl;
		cell_iter++;
	}

	boundary_condition();
	
	cout<<"Boundary condition is done."<<endl;
	
	cell_iter = cells.begin();

	//int i = 0;

	while(cell_iter != cells.end())
	{
		cell_iter->add_particle(&partile_buffer);
		cell_iter++;
		//cout<<"cell["<<i<<"]"<<" partile_buffer.size() = "<<partile_buffer.size()<<endl;
	}

	cout<<"Sending data...";
	send_data();
	cout<<"is over."<<endl<<"Reciving data...";
	recv_data();
	cout<<"is over."<<endl;




	t += dt;

	cout<<"t = "<<t<<"\tdt = "<<dt<<endl;
}

void cell_claster::collisions()
{
}

void cell_claster::boundary_condition()
{
	deque<boundary_ref>::iterator a = boundary_cond.begin();
	while(a != boundary_cond.end())
	{
		a->ref->bondary_condition(&partile_buffer);
		a++;
	}
}


bool cell_claster::send_data()
{
	int n = partile_buffer.size();
	particle *temp = new particle[n];
	int i = 0;
	deque<particle>::iterator particle_iter = partile_buffer.begin();
	while(particle_iter != partile_buffer.end())
	{
		temp[i].u = particle_iter->u;
		temp[i].v = particle_iter->v;
		temp[i].w = particle_iter->w;
		temp[i].x = particle_iter->x;
		temp[i].y = particle_iter->y;
		temp[i].z = particle_iter->z;

		i++;
		particle_iter++;
	}

	for(int i = 0; i < numproc; i++)
	{
		if(i != proc_id) MPI_Send(temp, n * 6, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
	}

	partile_buffer.clear();
	delete [] temp;

	return true;
}

bool cell_claster::recv_data()
{
	if(numproc > 1)
	{
		MPI_Status *status = new MPI_Status[numproc];
		bool recv = false;
		bool *step1 = new bool[numproc];
		bool *step2 = new bool[numproc];
		int k;
		int *recv_elem = new int[numproc];

		for(int i = 0; i < numproc; i++)
		{
				step1 = step2 = false;
		}

		step2[proc_id] = true;
		
		while(!recv)
		{
			for(int i = 0; i < numproc; i++)
			{
				if((i != proc_id) && !step1[i])
				{
					MPI_Iprobe(i, 1, MPI_COMM_WORLD, &k, &status[i]);
					step1[i] = (bool)k;
				}
				if((i != proc_id) && step1[i] && !step2[i])
				{
					MPI_Get_count(&status[i], MPI_DOUBLE, &recv_elem[i]);
					particle *temp = new particle[recv_elem[i]];
					MPI_Recv(temp, recv_elem[i], MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status[i]);
					
					while(cell_iter != cells.end())
					{
						cell_iter->add_particle(temp, recv_elem[i]);
						cell_iter++;
					}

					step2[i] = true;

					delete [] temp;
				}
			}
			for(int i = 0; i < numproc; i++)
			{
				recv = recv&&step2[i];
			}

		}
		
		return true;
	}
	return true;
}

bool cell_claster::write_file()
{
	std::ofstream file1("data1.dat");
	cell_iter = cells.begin();
	point ap;
	unsigned int N;
	while(cell_iter != cells.end())
	{
		ap = cell_iter->get_center();
		N = cell_iter->N();
		file1<<ap.x<<"\t"<<ap.y<<"\t"<<ap.z<<"\t"<<double(N) / np<<"\t"<<cell_iter->get_energy()<<endl;
		cell_iter++;
	}
	file1.close();
	return true;
}

void cell_claster::set_boundary_condition(boundary **a, int n)
{
	boundary_ref b;
	for(int i = 0; i < n; i++)
	{
		b.ref = a[i];
		boundary_cond.push_back(b);
		cout<<i<<" bound cond added"<<endl;
	}
}

void cell_claster::set_boundary_condition(boundary *a)
{
	boundary_ref b;
	b.ref = a;
	boundary_cond.push_back(b);
}

void cell_claster::computation(void)
{
	unsigned int a = 0;
	while(t < t_end)
	{

		if(a == 0) cout<<"dt = "<<t<<endl;
		this->time_step();
		std::cout<<"Time step:"<<a++<<endl;
	}
	//this->write_file();
}

void cell_claster::set_end_time(double t_end)
{
	this->t_end = t_end;
}
void cell_claster::calc_dt()
{
	dt = 1000000.;

	double dx = Lx / double(ncx), dy = Ly / double(ncy), dz = Lz / double(ncz);

	int i = 0;
	cell_iter = cells.begin();
	while(cell_iter != cells.end())
	{
		double dtt = cell_iter->calc_dt();
		if(dtt < dt) dt = dtt;
		cell_iter++;
	}
	dt *= Cu;
}

void cell_claster::set_dt_in_cells()
{
	cell_iter = cells.begin();
	while(cell_iter != cells.end())
	{
		cell_iter->set_dt(dt);
		cell_iter++;
	}
}

bool cell_claster::send_dt()
{
	for(int i = 0; i < numproc; i++)
	{
		MPI_Request req;
		int snd = false;
		if(i != proc_id) snd = MPI_Send(&dt, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
		cout<<"snd = "<<snd<<endl;
	}
	cout<<"dt was send."<<endl;
	return true;
}

bool cell_claster::recv_dt()
{
	if(numproc > 1)
	{
		for(int i = 0; i < numproc; i++)
		{
			double dtt;
			MPI_Status status;
			int rcv = false;
			if(i != proc_id) rcv = MPI_Recv(&dtt, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);
			cout<<"dtt = "<<dtt<<"\t rcv = "<<rcv<<endl;
			if(dtt < dt) dt = dtt;
		}
	}
	return true;
}