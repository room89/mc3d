//File: var.cpp 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.3
//Last modified: 9.07.09.
//Description: Program for calculation of freemolecular flows.

#include "var.h"

var::var(void)
{
}

var::var(double u, double v, double w, double t, unsigned int n)
{
	this->set(u, v, w, t, n);
}

var::~var(void)
{
}

void var::set(double u, double v, double w, double t, unsigned int n)
{
	this->n = n;
	this->u = u;
	this->v = v;
	this->w = w;
	this->t = t;
}

ostream & operator<<(ostream &o, const var &c)
{
	o<<double(c.n)/c.t<<" "<<double(c.n);
	return o;
}