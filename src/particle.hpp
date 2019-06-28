//File: particle.h 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.2.1
//Last modified: 24.12.08.
//Description: Program for calculation of freemolecular flows.

#pragma once

struct particle
{
	double x, y, z;
	double u, v, w;
	//double r;
	//friend void collision(particle a, particle b);
	particle(void);
	~particle(void);
	particle(double x, double y, double z, double u, double v, double w, double r);
	particle(double x, double y, double z, double u, double v, double w);
	void print();
	void move(double dt);
};

bool collision(particle *a, particle *b, double &g_max_ch, double &frequency_t_ch, double factor);
