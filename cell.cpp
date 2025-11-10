//File: cell.cpp 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.1
//Last modified: 15.05.09.
//Description: Program for calculation of freemolecular flows.

#include "cell.h"
#include <complex>
#include "particle.h"
#include "var.h"

cell::cell(void)
{
	np = 1;
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
	bool a = this->particles.gen_array(N);
	if(a)
	{
		this->u = 0;
		this->v = 0;
		this->w = 0;
		this->T = 1;
		this->n = N;
		this->generate_random();
	}
	dt = 100000;
	double c = sqrt(2 * T);
	double dtt = lx / (fabs(u) + c) > ly / (fabs(v) + c) ? lx / (fabs(u) + c) : ly / (fabs(v) + c);
	dtt = dtt > lz / (fabs(w) + c) ? dtt : lz / (fabs(w) + c);
	this->dt = dtt;
	return a;
}

inline unsigned int cell::N()
{
	return n;
}

double cell::generate_random()
{
	point d(lx, ly, lz);
	particles.genrate_random(T, u, v, w, apex, d);
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
	if(particles.get_N() < 2) return;
	double rmt = 1 / RAND_MAX;
	double g_max = 2 * sqrt(T);
	double lambda_t = 1 / (L * Kn_l);
	double factor = 1 / (0.5 * double(n) * lambda_t);
	double frequency_t = factor / g_max;

	double t = 0;
	double tau_mean = 0;

	p_elem *a, *b;

	this->particles.go_to_start();
	a = this->particles.get_data();
	this->particles.go_forward();
	b = this->particles.get_data();

	while(t < dt)
	{
		bool tau_bool;
		double r = rand() * rmt;

		if(r <= 0) r = 0.00001;

		if(!particles.go_to_start()) particles.go_to_start(); 
		a = this->particles.get_data();
		if(!particles.go_forward()) particles.go_to_start();
		b = this->particles.get_data();

		double tau = -frequency_t * log(r);
		tau_bool = collision(&a->p, &b->p, g_max, frequency_t, factor);

		if(tau_bool) t += tau;
	}
}

void cell::particle_move()
{
	p_elem *a;
	particles.go_to_start();
	u = v = w =0;
	for(unsigned int i = 0; i < n; i++)
	{
		a = particles.get_data();
		a->move_particle(dt);
		u += a->p.u;
		v += a->p.v;
		w += a->p.w;
		T = u * u + v * v + w * w;
		if(!particles.go_forward()) break;
	}
}

void cell::sort(void)
{
	if(body)
	{
		particle_buffer.add_particle_array(&particles);
		return;
	}
	particles.go_to_start();
	while(true)
	{
		p_elem *a;
		a = particles.get_data();
		if((a->p.x < apex.x) || (a->p.x > apex.x + lx) || (a->p.y < apex.y) || (a->p.y > apex.y + ly)
			|| (a->p.z < apex.z) || (a->p.z > apex.z + lz))
		{
			a = particles.remove();
			particle_buffer.add(a);
		}
		if(!particles.go_forward()) return;
	}
}

void cell::set_L(double L)
{
	this->L = L;
}

void cell::add_particle(particle_array* a)
{
	p_elem* temp;
	a->go_to_start();
	int i = 0;
	while(true)
	{
		//cout<<"work with particle "<<i++<<endl;
		temp = a->get_data();
		if((temp->p.x > apex.x) || (temp->p.x < apex.x + lx) || (temp->p.y > apex.y) || (temp->p.y < apex.y + ly)
			|| (temp->p.z > apex.z) || (temp->p.z < apex.z + lz))
		{
			this->particles.add(temp);
		}
		a->calc_N();
		if(!a->go_forward())
		{
			break;
		}
		return;
	}
}

particle_array* cell::get_buffer()
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

ostream & operator<<(ostream &o, const cell &c)
{
	return o<<c.apex.x<<"\t"<<c.apex.y<<"\t"<<c.apex.z<<"\t"<<double(c.n)*c.T<<"\t"<<double(c.n);
}