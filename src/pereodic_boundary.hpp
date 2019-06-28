//File: pereodic_boundary.h 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.1
//Last modified: 15.01.10.
//Description: Program for calculation of freemolecular flows.


#pragma once
#include "boundary.h"

class pereodic_boundary :
	public boundary
{
private:
	point mixing;
public:
	pereodic_boundary(point pstn, point nrml, point mixing);
	void set_mixing(point b);
	int bondary_condition(deque<particle> *cluster_particle);
	pereodic_boundary(void);
	~pereodic_boundary(void);
};
