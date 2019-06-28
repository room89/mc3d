//File: point.h 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.3
//Last modified: 9.07.09.
//Description: Program for calculation of freemolecular flows.

#pragma once

#include <iostream>
using namespace std;

struct point
{
	point(void);
	point(double x, double y, double z);
	~point(void);
	double x, y, z;
	void print();
	void set(double x, double y, double z);
	friend point operator+(point a, point b);
	friend point operator*(point a, point b);
	friend point operator*(point a, double b);
	friend point operator-(point a, point b);
	friend ostream & operator<<(ostream &o, const point &c);
};
