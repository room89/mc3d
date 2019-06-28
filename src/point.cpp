//File: point.cpp
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.3
//Last modified: 9.07.09.
//Description: Program for calculation of freemolecular flows.

#include "point.h"
#include <iostream>
#include <fstream>

point::point(void)
{
}

point::~point(void)
{
}

point::point(double x, double y, double z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

void point::print()
{
	std::cout<<"("<<this->x<<","<<this->y<<","<<this->z<<")"<<std::endl;
}

void point::set(double x, double y, double z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}


point operator+(point a, point b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}

point operator*(point a, point b)
{
	a.x *= b.x;
	a.y *= b.y;
	a.z *= b.z;
	return a;
}

point operator*(point a, double b)
{
	point c(b, b, b);
	return a * c;
}

point operator-(point a, point b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	return a;
}

std::ostream & operator<<(std::ostream &o, const point &c)
{
	o<<c.x<<" "<<c.y<<" "<<c.z;
	return o;
}