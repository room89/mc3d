//File: cell.h 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.1
//Last modified: 15.05.09.
//Description: Program for calculation of freemolecular flows.

#pragma once

#include "particle_array.h"
#include "point.h"
#include "var.h"
#include <iostream>

using namespace std;

class cell
{
private:
	point apex;
	double lx, ly, lz;
	double T;
	double u;
	double v;
	double w;
	double volume;
	unsigned int n;
	double Kn_l;										//локальый кнутсен
	double L;											//характерный размер
	double dt;
	double time;										
	double calc_time;									//врем€ вычислени€ в €чейке за один временной шаг dt
	bool body;											//содержание тела внутри €чейки
	unsigned int np;
	particle_array particles;
	particle_array particle_buffer;
	void particle_move(void);							//перемещение частиц
	void collisions(void);								//соударени€ между частицами
public:
	cell(void);
	~cell(void);
	unsigned int N(void);								//возвращает количство частиц в €чейке
	void set_size(double lx, double ly, double lz);		//установка размера €чейки
	void set_apex(point a);								//установка "опорной" точки
	void set_L(double L);								//установка характерного размера €чейки(одинаковый дл€ всех 3х измерений)
	bool initialazition(unsigned int N);				//инициализаци€ €чейки (N-количество частиц)
	double generate_random(void);
	void time_step(void);								//шаг по времени
	void sort(void);									//перемешение вылетевштх из€чейки частиц в буфер
	particle_array* get_buffer(void);					//возвращает буфер €чейки
	void add_particle(particle_array *a);				//забирает частицы наход€щиес€ внутри €чейки
	var get_var();
	point get_apex();
	point get_center();
	double get_dt();
	void set_dt(double dt);
	friend ostream & operator<<(ostream &o, const cell &c);
};
