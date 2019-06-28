//File: cell.h 
//Program: MC3D
//Author: Khokhlov "AAsad" Ivan
//Version: 0.3.9
//Last modified: 18.05.10.
//Description: Program for calculation of freemolecular flows.

#pragma once

#include "particle.h"
#include "point.h"
#include "var.h"
#include <iostream>
#include <deque>
#include <algorithm>

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
	unsigned int total_n;
	double Kn_l;										//локальый кнутсен
	double Kn;
	double L;											//характерный размер
	double dt;
	double time;										
	double calc_time;									//врем€ вычислени€ в €чейке за один временной шаг dt
	bool body;											//содержание тела внутри €чейки
	unsigned int np;
	deque<particle> particles;
	deque<particle> particle_buffer;
	void particle_move(void);							//перемещение частиц
	void collisions(void);								//соударени€ между частицами
public:
	cell(void);
	~cell(void);
	unsigned int N(void);								//возвращает количство частиц в €чейке
	double t();
	void set_size(double lx, double ly, double lz);		//установка размера €чейки
	void set_apex(point a);								//установка "опорной" точки
	void set_L(double L);								//установка характерного размера €чейки(одинаковый дл€ всех 3х измерений)
	void set_t(double t);
	void set_vel(double u, double v, double w);
	bool initialazition(unsigned int N);				//инициализаци€ €чейки (N-количество частиц)
	double generate_random(unsigned int N);
	void time_step(void);								//шаг по времени
	void sort(void);									//перемешение вылетевштх из€чейки частиц в буфер
	deque<particle>* get_buffer(void);					//возвращает буфер €чейки
	void add_particle(deque<particle> *a);				//забирает частицы наход€щиес€ внутри €чейки
	void add_particle(particle *a, int n);
	double calc_Kn();
	var get_var();
	point get_apex();
	point get_center();
	double get_Kn();
	double get_dt();
	double get_t();
	double get_u();
	double get_v();
	double get_w();
	double get_energy();
	void calc_vel();
	void calc_T();
	void set_dt(double dt);
	double calc_dt();
	friend ostream & operator<<(ostream &o, const cell &c);
	bool write_file(ofstream *file);
	bool write_file();
};
