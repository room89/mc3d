//File: timer.h 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.9
//Last modified: 18.05.10.
//Description: Program for calculation of freemolecular flows.

#pragma once
#include <mpi.h>
#include <deque>
#include "time.h"

using namespace std;

class timer
{
private:
	deque<time> time_deque;
public:
	timer(void);
	~timer(void);
	void checkpoint();
	double calc_av();
};
