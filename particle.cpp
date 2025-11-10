//File: particle.cpp 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.2.1
//Last modified: 24.12.08.
//Description: Program for calculation of freemolecular flows.

#include "particle.h"
#include <iostream>
#include <complex>

particle::particle(void)
{
	this->r = 0.001;
	this->u = this->v = this->w = this->x = this->y = this->z = 0;
}

particle::particle(double x, double y, double z, double u, double v, double w, double r)
{
	this->r = r;
	this->u = u;
	this->v = v;
	this->w = w;
	this->x = x;
	this->y = y;
	this->z = z;
}

particle::particle(double x, double y, double z, double u, double v, double w)
{
	this->r = 0.001;
	this->u = u;
	this->v = v;
	this->w = w;
	this->x = x;
	this->y = y;
	this->z = z;
}

particle::~particle(void)
{
}

void particle::print()
{
	std::cout<<" x = "<<x<<" y = "<<y<<" z = "<<z<<" u = "<<u<<" v = "<<v<<" w = "<<w<<std::endl;
}

void particle::move(double dt)
{
	x += dt * u;
	y += dt * v;
	z += dt * w;
}
bool collision(particle *a, particle *b, double &g_max, double &frequency_t, double factor)
{
	double rmt = 1 / RAND_MAX;
	double Pi = 3.1415923565;

	double u1 = a->u, v1 = a->v, w1 = a->w;
	double u2 = b->u, v2 = b->v, w2 = b->w;
	double gx = u2 - u1, gy = v2 - v1, gz = w2 - w1;
	double g = sqrt(gx * gx + gy * gy + gz * gz);

	bool rtrn = true;

	if(g > g_max)
	{
		g_max = g;
		frequency_t = factor / g_max;
		rtrn = false;
	}

	double r1 = rand() * rmt;
	double r2 = rand() * rmt;

	double g1z = g * cos(Pi * r1), g1y = g * sin(Pi * r1) * sin(2 * Pi * r2),
		g1x = g * sin(Pi * r1) * cos(2 * Pi * r2);

	double rr = rand() * rmt;
	
	if(g/g_max>rr)
	{
		double r1=rand()*rmt;
		double r2=rand()*rmt;

		double g1z=g*cos(Pi*r1), g1y=g*sin(Pi*r1)*sin(2.*Pi*r2), g1x=g*sin(Pi*r1)*cos(2.*Pi*r2);

		a->u = 0.5 * (u1 + u2) - 0.5 * g1x;
		a->v = 0.5 * (v1 + v2) - 0.5 * g1y;
		a->w = 0.5 * (w1 + w2) - 0.5 * g1z;

		b->u = 0.5 * (u1 + u2) + 0.5 * g1x;
		b->v = 0.5 * (v1 + v2) + 0.5 * g1y;
		b->w = 0.5 * (w1 + w2) + 0.5 * g1z;
	}
	

	return rtrn;
}