//File: p_elem.cpp 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.2.1
//Last modified: 24.12.08.
//Description: Program for calculation of freemolecular flows.

#include "p_elem.h"

void p_elem::move_particle(double dt)
{
	p.move(dt);
}