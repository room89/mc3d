//File: main.h 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.1
//Last modified: 17.04.08.
//Description: Program for calculation of freemolecular flows.

#pragma once
#include "particle_array.h"
#include <iostream>
#include "cell.h"
#include "point.h"
#include "cell_claster.h"
#include "pereodic_boundary.h"
//#include "header.h"
using namespace std;

const double Lx = 1, Ly = 1, Lz = 1;
typedef unsigned int uint;

int main()
{
	cell_claster claster;

	double Kn, Cu, T;
	unsigned int N = 0;
	unsigned int ncx = 0, ncy = 0, ncz = 0;
	//double dx, dy, dz;
	unsigned int np = 100;

	bool mark = true;
	while(mark)
	{
		cout<<"Input Knudsen number: ";
		cin>>Kn;
		if(Kn <= 0)
		{
			cout<<"Knudsen number must be > 0, input other Knudsen number."<<endl;
		}
		else mark = false;
	}

	mark = true;
	while(mark)
	{
		cout<<"Input Courant number, Cu<=1: ";
		cin>>Cu;
		if(Cu > 1) cout<<"Courant number must be <= 1, input other Courant number."<<endl;
		else if(Cu <= 0) cout<<"Courant number must be > 0, input other Courant number."<<endl;
		else mark = false;
	}

	cout<<"Input time: ";
	cin>>T;
	if(T <= 0)
	{
		cout<<"T <= 0, breaking the calculation."<<endl;
		return 0;
	}

	mark = true;
	while(mark)
	{
		cout<<"Input dimension, x: ";
		cin>>ncx;
		if(ncx < 10) cout<<"ncx is too small, input other dimension."<<endl;
		else if(ncx < unsigned int(Lx / Kn)) cout<<"number of cells is not sufficient, min ncx = "
			<<unsigned int(Lx / Kn)<<endl;
		else mark = false;
	}

	mark = true;
	while(mark)
	{
		cout<<"Input dimension, y: ";
		cin>>ncy;
		if(ncy < 10) cout<<"ncy is too small, input other dimension."<<endl;
		else if(ncy < unsigned int(Ly / Kn)) cout<<"number of cells is not sufficient, min ncx = "
			<<unsigned int(Lx / Kn)<<endl;
		else mark = false;
	}

	mark = true;
	while(mark)
	{
		cout<<"Input dimension, z: ";
		cin>>ncz;
		if(ncz < 10) cout<<"ncz is too small, input other dimension."<<endl;
		else if(ncz < unsigned int(Lz / Kn)) cout<<"number of cells is not sufficient, min ncx = "
			<<unsigned int(Lx / Kn)<<endl;
		else mark = false;
	}

	mark = true;
	while(mark)
	{
		cout<<"Input statistics level: ";
		cin>>np;
		if(np < 3) cout<<"np is too small, input other statistics level"<<endl;
		else mark = false;
	}

	point a;
	a.x = -Lx * 0.5;
	a.y = -Ly * 0.5;
	a.z = -Lz * 0.5;

	claster.set_apex(a);
	claster.set_size(Lx, Ly, Lz);

	claster.initialazition(ncx, ncy, ncz, np);
	
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

	for(int i = 0; i < 6; i++) claster.set_boundary_condition(boundary_cond[i], i);

	claster.set_end_time(T);

	claster.computation();
	claster.write_file();
	
	bool rjr;
	cin>>rjr;

	//delete [] cells;

	return 0;
}