//File: p_elem.h 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.2.1
//Last modified: 24.12.08.
//Description: Program for calculation of freemolecular flows.

#pragma once
#include "particle.h"


struct p_elem
{
public:
	particle p;
	p_elem *prev;
	p_elem *next;
	void move_particle(double dt);
};
