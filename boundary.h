//File: boundary.h
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.2
//Last modified: 21.04.09.
//Description: Program for calculation of freemolecular flows.

#pragma once

#include "point.h"
#include "particle_array.h"

virtual class boundary
{
private:
public:
	point nrml;												//направление граничных условий
	point pstn;												//положение граничных условий

	virtual int bondary_condition(particle_array *a) = 0;	//абстракнтый метод выполнения граничных условий с масивом частиц a
	void set_position(point a);								//установка места граничных условий
	void set_normal(point a);								//установка направления граничных условий
};
