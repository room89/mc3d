// File: particle_array.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.3.1
// Last modified: 16.01.09.
// Description: Program for calculation of freemolecular flows.

#include "particle_array.h"

#include <complex>
#include <cstdlib>

const double Pi = 3.1415926535;

particle_array::particle_array(void) {
  data = NULL;
  end = NULL;
  start = NULL;
  N = 0;
}

bool particle_array::gen_array(int n) {
  particle a;
  for (int i = 0; i < n; i++) {
    if (!this->add(a)) break;
  }
  if (this->N == n) {
    return true;
  } else
    return false;
}

particle_array::~particle_array(void) {
  if (this->data == NULL) return;
  p_elem *del;
  this->data = this->start;
  while (this->data->next) {
    del = this->data;
    this->data = this->data->next;
    delete del;
  }
  delete this->data;
  end = NULL;
  start = NULL;
  data = NULL;
  N = 0;
}

bool particle_array::add(particle add_p) {
  if (N == 0) {
    data = new p_elem;
    data->p = add_p;
    data->next = 0;
    data->prev = 0;
    end = data;
    start = data;
    N++;
    return true;
  }
  p_elem *temp_a;
  p_elem *temp_b;
  temp_a = data;
  temp_b = data->next;
  data->next = new p_elem;
  data = data->next;
  data->p = add_p;
  data->prev = temp_a;
  data->next = temp_b;
  N++;
  if (data->next == 0) this->end = data;
  if (data->prev == 0) this->start = data;
  return true;
}

bool particle_array::add(p_elem *add) {
  /*if((add->next != NULL) && (add->prev != NULL))
  {
          add->prev->next = add->next;
          add->next->prev = add->prev;
  }
  else if(add->next == NULL)
  {
          add->prev->next = NULL;
  }*/
  /*if((add->next != NULL) || (add->prev != NULL))
  {
          if(add->prev != NULL) add->prev->next = add->next;
          add->next->prev = add->prev;
  }
  if(data == NULL)
  {
          data = add;
          end = add;
          start = add;
          add->next = add->prev = NULL;
  }
  else
  {
          end->next = add;
          add->prev = end;
          end = add;
  }
  return true;*/
  if ((add->next != NULL) && (add->prev != NULL)) {
    add->next->prev = add->prev;
    add->prev->next = add->next;
  } else if (add->next == NULL) {
    add->prev->next = NULL;
  } else if (add->prev == NULL) {
    add->next->prev == NULL;
  }
  if (data == NULL) {
    data = add;
    start = add;
    end = add;
    N++;
    return true;
  }
  p_elem *temp = data->next;
  data->next = add;
  data->next->next = temp;
  data->next->prev = data;
  temp->prev = data->next;
  N++;
  return true;
}

bool particle_array::go_back() {
  if (data->prev != 0) {
    data = data->prev;
    return true;
  } else
    return false;
}

bool particle_array::go_forward() {
  if (data->next != NULL) {
    data = data->next;
    return true;
  } else
    return false;
}

bool particle_array::del() {
  if (N == 0) return false;
  p_elem *del;
  del = data;
  if (data->next)
    data->next->prev = data->prev;
  else
    this->end = data->prev;
  if (data->prev) {
    data->prev->next = data->next;
    data = data->prev;
  } else if (this->data->next)
    data = data->next;
  else
    return false;
  delete del;
  N--;
  return true;
}

int particle_array::num_of_elem() { return this->N; }

p_elem *particle_array::get_data() { return data; }

p_elem *particle_array::get_end() { return end; }

p_elem *particle_array::get_start() { return start; }

bool particle_array::add_particle_array(
    particle_array *a)  // возможна утечка памяти...?
{
  if (this->N == 0) {
    this->end = a->end;
    start = a->start;
    data = a->start;
    N = a->N;
    a->clear();
    return true;
  }
  end->next = a->data;
  a->data->prev = end;
  N += a->num_of_elem();
  a->clear();
  return true;
}

void particle_array::genrate_random(double t, double u, double v, double w,
                                    point apex, point d) {
  double rmt = 1 / double(RAND_MAX);
  double ti = 0;
  unsigned int nn = N;
  double nt = 1 / double(N);
  p_elem *data2 = this->end;

  if ((N % 2) != 0) {
    nn -= 1;
    this->end->p.u = 0;
    this->end->p.v = 0;
    this->end->p.w = 0;
    data2 = this->end->prev;
  }

  unsigned int nn2 = nn / 2;

  this->go_to_start();

  for (unsigned int i = 0; i < nn2; i++) {
    double rn1 = double(rand()) * rmt;
    double rn2 = double(rand()) * rmt;

    if (rn1 <= 0) rn1 = 0.00001;

    double slg = sqrt(2 * fabs(log(rn1)));

    this->data->p.u = slg * cos(2 * Pi * rn2);
    this->data->p.v = slg * sin(2 * Pi * rn2);

    data2->p.u = -slg * cos(2 * Pi * rn2);
    data2->p.v = -slg * sin(2 * Pi * rn2);

    double rn3 = rand() * rmt;
    double rn4 = rand() * rmt;

    if (rn3 <= 0) rn3 = 0.00001;

    slg = sqrt(2 * fabs(log(rn3)));
    this->data->p.w = slg * cos(2 * Pi * rn4);
    data2->p.w = -slg * cos(2 * Pi * rn4);

    ti += 2 * (this->data->p.u * this->data->p.u +
               this->data->p.v * this->data->p.v + this->data->p.w +
               this->data->p.w);

    this->data = this->data->next;
    data2 = data2->prev;
  }

  ti = ti * nt / 3.;

  double sf = sqrt(1 / ti);

  data = start;

  for (unsigned int i = 0; i < N; i++) {
    this->data->p.v *= sf;
    this->data->p.u *= sf;
    this->data->p.w *= sf;
    if (i != N - 1) this->data = this->data->next;
  }

  ti = 0;

  this->data = this->start;

  for (unsigned int i = 0; i < N; i++) {
    ti += (this->data->p.u) * (this->data->p.u) +
          (this->data->p.v) * (this->data->p.v) +
          (this->data->p.w) * (this->data->p.w);
    this->data = this->data->next;
  }

  this->data = this->start;

  for (unsigned int i = 0; i < N; i++) {
    data->p.u = sqrt(t) * data->p.u + data->p.u;
    data->p.v = sqrt(t) * data->p.v + data->p.v;
    data->p.w = sqrt(t) * data->p.w + data->p.w;

    double rnx = rand() * rmt;
    double rny = rand() * rmt;
    double rnz = rand() * rmt;

    data->p.x = apex.x + d.x * rnx;
    data->p.y = apex.y + d.y * rny;
    data->p.z = apex.z + d.z * rnz;
  }
}  // не доделано

bool particle_array::clear() {
  end = NULL;
  start = NULL;
  data = NULL;
  N = 0;
  return true;
}

bool particle_array::go_to_start() {
  if (this->start) {
    this->data = this->start;
    return true;
  } else
    return false;
}

p_elem *particle_array::remove(void) {
  if (data) {
    p_elem *rem, *temp1, *temp2;
    rem = data;
    temp1 = data->next;
    temp2 = data->prev;
    temp1->prev = temp2;
    temp2->next = temp1;
    return rem;
  } else
    return NULL;
}

unsigned int particle_array::get_N() { return N; }

unsigned int particle_array::calc_N() {
  p_elem *temp = this->data;
  this->go_to_start();
  unsigned int i = 0;
  while (this->go_forward()) i++;
  this->N = i;
  this->data = temp;
  return i;
}
