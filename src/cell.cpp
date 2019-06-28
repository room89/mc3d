//File: cell.cpp 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.9
//Last modified: 18.05.10.
//Description: Program for calculation of freemolecular flows.

#include "cell.h"
#include <complex>
#include "particle.h"
#include "var.h"
#include <iostream>
#include <fstream>
#include <algorithm>

const double Pi = 3.1415926535;

cell::cell(void)
{
	np = 1;
	body = false;
}

cell::~cell(void)
{
}

void cell::set_size(double lx, double ly, double lz)
{
	if((lx <= 0) || (ly <= 0) || (lz <= 0)) return;
	this->lx = lx;
	this->ly = ly;
	this->lz = lz;
}

void cell::set_apex(point a)
{
	apex = a;
}

bool cell::initialazition(unsigned int N)
{
	this->n = N;
	this->generate_random(n);
	
	dt = 100000;
	double c = sqrt(2 * T);
	double dtt = lx / (fabs(u) + c) < ly / (fabs(v) + c) ? lx / (fabs(u) + c) : ly / (fabs(v) + c);
	dtt = dtt < lz / (fabs(w) + c) ? dtt : lz / (fabs(w) + c);
	this->dt = dtt;
	return true;
}

unsigned int cell::N()
{
	return n = particles.size();
}

double cell::generate_random(unsigned int N)
{
	point d(lx, ly, lz);
	
	double rmt = 1 / double(RAND_MAX);
	double ti = 0;
	deque<particle>::iterator data;

	unsigned int nn = N;
	double nt = 1 / double(N);
	
	data = particles.begin();

	unsigned int nn2 = (nn - nn%2) / 2;

	particle p1, p2;

	for(unsigned int i = 0; i < nn2; i++)
	{
		double rn1 = double(rand()) * rmt;
		double rn2 = double(rand()) * rmt;

		if(rn1 <= 0) rn1 = 0.00001;

		double slg = sqrt(2 * fabs(log(rn1)));

		p1.u = slg * cos(2 * Pi * rn2);
		p1.v = slg * sin(2 * Pi * rn2);


		p2.u = -slg * cos(2 * Pi * rn2);
		p2.v = -slg * sin(2 * Pi * rn2);

		double rn3 = rand() * rmt;
		double rn4 = rand() * rmt;

		if(rn3 <= 0) rn3 = 0.00001;

		slg = sqrt(2 * fabs(log(rn3)));

		p1.w = slg * cos(2 * Pi * rn4);
		p2.w = -slg * cos(2 * Pi * rn4);

		particles.push_back(p1);
		particles.push_back(p2);

		ti += 2 * (p1.u * p1.u + p1.v * p1.v + p1.w + p1.w);
	}

	data = particles.begin();

	if(N%2 == 1)
	{
		p1.u = p1.v = p1.w = 0;
		particles.push_back(p1);
	}

	ti = ti * nt / 3.;

	double sf = sqrt(1 / ti);

	data = particles.begin();

	while(data != particles.end())
	{
		data->u *= sf;
		data->v *= sf;
		data->w *= sf;
		data++;
	}

	ti = 0;
	data = particles.begin();

	while(data != particles.end())
	{
		ti += (data->u) * (data->u) + (data->v) * (data->v) + (data->w) * (data->w);
		data++;
	}

	data = particles.begin();

	while(data != particles.end())
	{
		data->u = u + sqrt(T) * data->u;
		data->v = v + sqrt(T) * data->v;
		data->w = w + sqrt(T) * data->w;

		double rnx = rand() * rmt;
		double rny = rand() * rmt;
		double rnz = rand() * rmt;

		data->x = apex.x + d.x * rnx;
		data->y = apex.y + d.y * rny;
		data->z = apex.z + d.z * rnz;

		data++;
	}

	random_shuffle(particles.begin(), particles.end());

	return 0;
}

void cell::time_step()
{
	this->collisions();
	this->particle_move();
	this->sort();
}

void cell::collisions()
{
	calc_Kn();
	double rmt = 1 / RAND_MAX;
	double g_max = 2 * sqrt(T);
	double lambda_t = 1 / (L * Kn_l);
	double factor = 1 / (0.5 * double(n) * lambda_t);
	double frequency_t = factor / g_max;

	double t = 0;
	double tau_mean = 0;

	deque<particle>::iterator a, b;

	if(particles.size() <= 1) return;

	a = b = particles.begin();
	b++;

	while(t <= dt)
	{
		double r = rand() * rmt;

		if(r <= 0) r = 0.00001;

		double tau = -frequency_t * log(r);
		t += tau;
		if(t > dt) break;

		double rmt = 1 / RAND_MAX;
		double Pi = 3.1415923565;

		double u1 = a->u, v1 = a->v, w1 = a->w;
		double u2 = b->u, v2 = b->v, w2 = b->w;
		double gx = u2 - u1, gy = v2 - v1, gz = w2 - w1;
		double g = sqrt(gx * gx + gy * gy + gz * gz);

		if(g_max < g)
		{
			g_max = g;
			t -=tau;
			frequency_t = factor / g_max;
		}

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
		if(b + 1 != particles.end())
		{
			if(b + 2 != particles.end())
			{
				a += 2;
				b += 2;
			}
			a = b = particles.begin();
			b++;
			random_shuffle(particles.begin(), particles.end());
		}
		else
		{
			a = b = particles.begin();
			b++;
			random_shuffle(particles.begin(), particles.end());
		}
	}
}


void cell::particle_move()
{
	deque<particle>::iterator a;
	a = particles.begin();
	u = v = w = 0;
	while(a != particles.end())
	{
		a->x += a->u * dt;
		a->y += a->v * dt;
		a->z += a->w * dt;

		a++;
	}

	this->calc_T();
}

void cell::sort(void)
{
	if(!body)
	{
		particle_buffer.insert(particle_buffer.end(), particles.begin(), particles.end());
		particles.clear();
		return;
	}
		
	deque<particle>::iterator data = particles.begin();

	while(data != particles.end())
	{
		if((data->x < apex.x) || (data->x > apex.x + lx) || (data->y < apex.y) || (data->y > apex.y + ly)
			|| (data->z < apex.z) || (data->z > apex.z + lz))
		{
			particle_buffer.insert(particle_buffer.end(), data, data + 1);
			data = particles.erase(data);
		}
		else data++;
	}
}

void cell::set_L(double L)
{
	this->L = L;
}

void cell::add_particle(deque<particle> *a)
{
	deque<particle>::iterator temp;
	temp = a->begin();
	while(temp != a->end())
	{
		if((temp->x > apex.x) && (temp->x < apex.x + lx) && (temp->y > apex.y) && (temp->y < apex.y + ly)
			&& (temp->z > apex.z) && (temp->z < apex.z + lz))
		{
			particles.insert(particles.end(), temp, temp + 1);
			temp = a->erase(temp);
		}
		else temp++;
	}
}

void cell::add_particle(particle *a, int n)
{

	for(int i = 0; i < n; i++)
	{
		if((a[i].x > apex.x) && (a[i].x < apex.x + lx) && (a[i].y > apex.y) && (a[i].y < apex.y + ly)
			&& (a[i].z > apex.z) && (a[i].z < apex.z + lz))
		{
			particles.push_back(a[i]);
		}
	}
}


deque<particle>* cell::get_buffer()
{
	return &particle_buffer;
}

var cell::get_var()
{
	var a(u, v, w, T, n);
	return a;
}

point cell::get_apex()
{
	return apex;
}

point cell::get_center()
{
	point a(lx / 2, ly / 2, lz / 2);
	return a = a + apex;
}

double cell::get_dt()
{
	return dt;
}

void cell::set_dt(double dt)
{
	this->dt = dt;
}

double cell::t()
{
	return T;
}


void cell::set_vel(double u, double v, double w)
{
	this->u = u;
	this->v = v;
	this->w = w;
}

void cell::set_t(double t)
{
	this->T = t;
}

double cell::get_Kn()
{
	return Kn;
}

double cell::calc_Kn()
{
	Kn = Kn * np / particles.size();
	return Kn;
}

double cell::get_t()
{
	calc_T();
	return T;
}

double cell::get_u()
{
	deque<particle>::iterator data = particles.begin();
	double av_u = 0;
	while(data != particles.end())
	{
		av_u += data->u;
		data++;
	}
	av_u /= particles.size();
	u = av_u;
	return u;
}

double cell::get_v()
{
	deque<particle>::iterator data = particles.begin();
	double av_v = 0;
	while(data != particles.end())
	{
		av_v += data->v;
		data++;
	}
	av_v /= particles.size();
	v = av_v;
	return av_v;
}

double cell::get_w()
{
	deque<particle>::iterator data = particles.begin();
	double av_w = 0;
	while(data != particles.end())
	{
		av_w += data->w;
		data++;
	}
	av_w /= particles.size();
	w = av_w;
	return av_w;
}


double cell::get_energy()
{
	double E = 0;
	int N = particles.size();
	deque<particle>::iterator data = particles.begin();
	while(data != particles.end())
	{
		E += data->u * data->u + data->w * data->w + data->v * data->v;
		data++;
	}
	E /= N;
	return E;
}

void cell::calc_vel()
{
	deque<particle>::iterator data = particles.begin();
	double av_u = 0, av_v = 0, av_w = 0;
	while(data != particles.end())
	{
		av_u += data->u;
		av_v += data->v;
		av_w += data->w;
		data++;
	}
	av_u /= particles.size();
	av_v /= particles.size();
	av_w /= particles.size();
	u = av_u;
	v = av_v;
	w = av_w;
}

void cell::calc_T()
{
	calc_vel();
	T = (get_energy() - u * u - v * v - w * w) / 3;
}

double cell::calc_dt()
{
	double c = sqrt(2 * T);
	double dtt = min(lx / (fabs(u) + c), ly / fabs(v) + c);
	dtt = min(dtt, fabs(lz / (w)) + c);
	if(dtt < dt) dt = dtt;
	return dt;
}