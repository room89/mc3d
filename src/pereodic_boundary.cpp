//File: pereodic_boundary.cpp 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.1
//Last modified: 21.12.09.
//Description: Program for calculation of freemolecular flows.

#include "pereodic_boundary.h"

pereodic_boundary::pereodic_boundary(void)
{
}

pereodic_boundary::~pereodic_boundary(void)
{
}

int pereodic_boundary::bondary_condition(deque<particle> *particles)
{
	deque<particle>::iterator data = particles->begin();
	while(data != particles->end())
	{
		while( (data->x - pstn.x) * nrml.x + (data->y - pstn.y) * nrml.y + (data->z - pstn.z) * nrml.z > 0)
		{
			data->x += mixing.x;
			data->y += mixing.y;
			data->z += mixing.z;
		}
		data++;
	}
	return 0;
}

pereodic_boundary::pereodic_boundary(point pstn, point nrml, point mixing)
{
	this->nrml = nrml;
	this->pstn = pstn;
	this->mixing = mixing;
}

void pereodic_boundary::set_mixing(point b)
{
	mixing = b;
}