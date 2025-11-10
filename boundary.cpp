//File: boundary.cpp
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.2.1
//Last modified: 16.01.09.
//Description: Program for calculation of freemolecular flows.


#include "boundary.h"

void boundary::set_normal(point a)
{
	nrml = a;
}

void boundary::set_position(point a)
{
	pstn = a;
}