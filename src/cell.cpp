// File: cell.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.5.1
// Last modified: 18.05.10.
// Description: Program for calculation of rarefaid flows.

#include "cell.h"

namespace mc3d {
namespace {
const double Pi = 3.14159265358979;
}

cell::cell(void) {
  this->calc_time = 0;
  //		this->start_time = MPI_Wtime();
  np = 1;
  body_mark = false;
}

cell::~cell(void) {
  particles.clear();
  particle_buffer.clear();
  neighbors.clear();
}

void cell::set_size(double lx, double ly, double lz) {
  if ((lx <= 0) || (ly <= 0) || (lz <= 0)) return;
  this->lx = lx;
  this->ly = ly;
  this->lz = lz;
  L = lx > ly ? lx : ly;
  L = L > lz ? L : lz;
}

void cell::set_size(point dl) {
  if ((dl.x <= 0) || (dl.y <= 0) || (dl.z <= 0)) return;
  lx = dl.x;
  ly = dl.y;
  lz = dl.z;
  L = lx > ly ? lx : ly;
  L = L > lz ? L : lz;
}

void cell::set_apex(point a) { apex = a; }

bool cell::initialazition(size_t N, const std::unique_ptr<geometry>& bbody) {
  this->n = N;
  this->generate_random(n);

  dt = 100000;
  calc_T();
  double c = sqrt(2 * T);
  double dtt = lx / (fabs(velocity.x) + c) < ly / (fabs(velocity.y) + c)
                   ? lx / (fabs(velocity.x) + c)
                   : ly / (fabs(velocity.y) + c);
  dtt = dtt < lz / (fabs(velocity.z) + c) ? dtt : lz / (fabs(velocity.z) + c);
  this->dt = dtt;
  random_shuffle(particles.begin(), particles.end());

  if (bbody) {
    body_mark = body_boundary.add_poligon(bbody.get(), get_center(), L);
    body_boundary.set_geometry(bbody.get());
  } else
    body_mark = false;

  return true;
}

unsigned int cell::N() { return n = particles.size(); }

double cell::generate_random(size_t N) {
  point d(lx, ly, lz);

  double rmt = 1 / double(RAND_MAX);
  double ti = 0;
  deque<particle>::iterator data;

  size_t nn = N;
  double nt = 1 / double(N);

  data = particles.begin();

  size_t nn2 = (nn - nn % 2) / 2;

  particle p1, p2;

  for (unsigned int i = 0; i < nn2; i++) {
    double rn1 = double(std::rand()) * rmt;
    double rn2 = double(std::rand()) * rmt;

    if (rn1 <= 0) rn1 = 0.00000001;

    double slg = sqrt(2 * fabs(log(rn1)));

    /*p1.u = slg * cos(2 * Pi * rn2);
    p1.v = slg * sin(2 * Pi * rn2);*/

    p1.set_u(slg * cos(2 * Pi * rn2));
    p1.set_v(slg * sin(2 * Pi * rn2));

    /*p2.u = -slg * cos(2 * Pi * rn2);
    p2.v = -slg * sin(2 * Pi * rn2);*/

    p2.set_u(-slg * cos(2 * Pi * rn2));
    p2.set_v(-slg * sin(2 * Pi * rn2));

    double rn3 = std::rand() * rmt;
    double rn4 = std::rand() * rmt;

    if (rn3 <= 0) rn3 = 0.00000001;

    slg = sqrt(2 * fabs(log(rn3)));

    /*p1.w = slg * cos(2 * Pi * rn4);
    p2.w = -slg * cos(2 * Pi * rn4);*/

    p1.set_w(slg * cos(2 * Pi * rn4));
    p2.set_w(-slg * cos(2 * Pi * rn4));

    particles.push_front(p1);
    particles.push_front(p2);

    /*ti += 2 * (p1.u * p1.u + p1.v * p1.v + p1.w * p1.w);*/

    ti += 2 * (p1.get_velocity() * p1.get_velocity());
  }

  data = particles.begin();

  if (N % 2 == 1) {
    // p1.u = p1.v = p1.w = 0;
    p1.set_velocity(point(0, 0, 0));
    particles.push_front(p1);
  }

  ti = ti * nt / 3.;

  double sf = sqrt(1 / ti);

  data = particles.begin();

  for (unsigned int i = 0; i < nn; i++) {
    if (data == particles.end()) break;
    /*data->u *= sf;
    data->v *= sf;
    data->w *= sf;*/
    data->velocity *= sf;
    data++;
  }

  data = particles.begin();

  for (unsigned int i = 0; i < nn; i++) {
    if (data == particles.end()) break;
    data->velocity = data->velocity * sqrt(T) + velocity;
    // data->u = u + sqrt(T) * data->u;
    // data->v = v + sqrt(T) * data->v;
    // data->w = w + sqrt(T) * data->w;

    double rnx = std::rand() * rmt;
    double rny = std::rand() * rmt;
    double rnz = std::rand() * rmt;

    point pos(d.x * rnx, d.y * rny, d.z * rnz);

    data->position = apex + pos;

    /*data->x = apex.x + d.x * rnx;
    data->y = apex.y + d.y * rny;
    data->z = apex.z + d.z * rnz;*/

    data++;
  }

  calc_T();

  return 0;
}

double cell::generate_random(unsigned int N, double T, point V) {
  point d(lx, ly, lz);

  double rmt = 1 / double(RAND_MAX);
  double ti = 0;
  deque<particle>::iterator data;

  unsigned int nn = N;
  double nt = 1 / double(N);

  data = particles.begin();

  unsigned int nn2 = (nn - nn % 2) / 2;

  particle p1, p2;

  for (unsigned int i = 0; i < nn2; i++) {
    double rn1 = double(std::rand()) * rmt;
    double rn2 = double(std::rand()) * rmt;

    if (rn1 <= 0) rn1 = 0.00001;

    double slg = sqrt(2 * fabs(log(rn1)));

    /*p1.u = slg * cos(2 * Pi * rn2);
    p1.v = slg * sin(2 * Pi * rn2);*/

    p1.set_u(slg * cos(2 * Pi * rn2));
    p1.set_v(slg * sin(2 * Pi * rn2));

    /*p2.u = -slg * cos(2 * Pi * rn2);
    p2.v = -slg * sin(2 * Pi * rn2);*/

    p2.set_u(-slg * cos(2 * Pi * rn2));
    p2.set_v(-slg * sin(2 * Pi * rn2));

    double rn3 = std::rand() * rmt;
    double rn4 = std::rand() * rmt;

    if (rn3 <= 0) rn3 = 0.00001;

    slg = sqrt(2 * fabs(log(rn3)));

    /*p1.w = slg * cos(2 * Pi * rn4);
    p2.w = -slg * cos(2 * Pi * rn4);*/

    p1.set_w(slg * cos(2 * Pi * rn4));
    p2.set_w(-slg * cos(2 * Pi * rn4));

    particles.push_front(p1);
    particles.push_front(p2);

    // ti += 2 * (p1.u * p1.u + p1.v * p1.v + p1.w * p1.w);

    ti += 2 * (p1.get_velocity() * p1.get_velocity());
  }

  data = particles.begin();

  if (N % 2 == 1) {
    // p1.u = p1.v = p1.w = 0;
    p1.set_velocity(point(0, 0, 0));
    particles.push_front(p1);
  }

  ti = ti * nt / 3.;

  double sf = sqrt(1 / ti);

  data = particles.begin();

  for (unsigned int i = 0; i < nn; i++) {
    if (data == particles.end()) break;
    /*data->u *= sf;
    data->v *= sf;
    data->w *= sf;*/

    data->velocity *= sf;
    data++;
  }

  data = particles.begin();

  for (unsigned int i = 0; i < nn; i++) {
    if (data == particles.end()) break;
    /*data->u = V.x + sqrt(T) * data->u;
    data->v = V.y + sqrt(T) * data->v;
    data->w = V.z + sqrt(T) * data->w;*/

    data->velocity *= sqrt(T);
    data->velocity += V;

    double rnx = std::rand() * rmt;
    double rny = std::rand() * rmt;
    double rnz = std::rand() * rmt;

    /*data->x = apex.x + d.x * rnx;
    data->y = apex.y + d.y * rny;
    data->z = apex.z + d.z * rnz;*/

    data->position = apex + point(d.x * rnx, d.y * rny, d.z * rnz);

    data++;
  }

  random_shuffle(particles.begin(), particles.end());

  return 0;
}

double cell::generate_free_random(unsigned int N, double T, point V,
                                  point nrml) {
  point d(lx, ly, lz);

  double rmt = 1 / double(RAND_MAX);
  double ti = 0;
  deque<particle>::iterator data;
  deque<particle> added_particles;

  double nt = 1 / double(N);

  data = added_particles.begin();
  particle p1;
  point dvel(0, 0, 0);

  unsigned int i = 0;

  while (i < N) {
    double rn1 = double(rand()) * rmt;
    double rn2 = double(rand()) * rmt;

    if (rn1 <= 0) rn1 = 0.00001;

    double slg = sqrt(2 * T * fabs(log(rn1)));

    p1.set_u(slg * cos(2 * Pi * rn2));
    p1.set_v(slg * sin(2 * Pi * rn2));

    double rn3 = rand() * rmt;
    double rn4 = rand() * rmt;

    if (rn3 <= 0) rn3 = 0.00001;

    slg = sqrt(2 * T * fabs(log(rn3)));

    p1.set_w(slg * cos(2 * Pi * rn4));

    if ((p1.get_velocity() + V) * nrml < 0) {
      added_particles.push_front(p1);
      ti += (p1.velocity * p1.velocity);
      dvel += p1.velocity;
      i++;
    }
  }

  dvel *= nt;
  ti = (ti * nt - dvel * dvel) / 3;
  /*
                  double sf = sqrt(T / ti);
  //		sf = 1;

                  data = added_particles.begin();

                  for(i = 0; i < N; i++)
                  {
                          if(data == added_particles.end()) break;
                          data->velocity = (data->velocity - dvel) * sf + V;
                          if(data->velocity * nrml < 0) data->velocity = 2 * V -
  data->velocity; // data++;
                  }
  */
  data = added_particles.begin();

  for (i = 0; i < N; i++) {
    if (data == added_particles.end()) break;

    data->velocity += V;  //

    double rnx = rand() * rmt;
    double rny = rand() * rmt;
    double rnz = rand() * rmt;

    data->position = apex + point(d.x * rnx, d.y * rny, d.z * rnz);

    data++;
  }

  //              random_shuffle(added_particles.begin(),
  //              added_particles.end());

  particles.insert(particles.end(), added_particles.begin(),
                   added_particles.end());

  return 0;
}

/*double cell::generate_giper_free_random(unsigned int N, point V, double T)
{
         deque<particle> added_particles;
         double rmt = 1 / double(RAND_MAX);
         point position;
         point d(lx, ly, lz);

         for(unsigned int i = 0; i < N; i++)
         {
                double rnx = rand() * rmt;
                double rny = rand() * rmt;
                double rnz = rand() * rmt;

                position = apex + point(d.x * rnx, d.y * rny, d.z * rnz);

                particle new_particle;

                new_particle.position = position;
                point noise(rand() * rmt - .5, rand() * rmt - .5, rand() * rmt -
.5); new_particle.velocity = V + sqrt(2 * T) * noise;

                added_particles.push_back(new_particle);
         }

         particles.insert(particles.end(), added_particles.begin(),
added_particles.end());

         return 0;
}*/

double cell::generate_giper_free_random(unsigned int N, point V, double T) {
  deque<particle> added_particles;
  double rmt = 1 / double(RAND_MAX);
  point position;
  point d(lx, ly, lz);
  double A = sqrt(3 * T);

  for (unsigned int i = 0; i < N; i++) {
    double rnx = rand() * rmt;
    double rny = rand() * rmt;
    double rnz = rand() * rmt;

    position = apex + point(d.x * rnx, d.y * rny, d.z * rnz);

    particle new_particle;

    new_particle.position = position;
    point noise(A * (2 * rand() * rmt - 1), A * (2 * rand() * rmt - 1),
                A * (2 * rand() * rmt - 1));
    new_particle.velocity = V + noise;

    added_particles.push_back(new_particle);
  }

  particles.insert(particles.end(), added_particles.begin(),
                   added_particles.end());

  return 0;
}

void cell::collisions() {
  if (particles.size() <= 5) return;

  // cout << calc_vel() << endl;

  calc_Kn();

  double rmt = 1 / double(RAND_MAX);
  double g_max = 2 * sqrt(calc_T());
  double factor = 2 * sqrt(2.) * L * Kn_l / particles.size();
  double frequency_t = factor / g_max;

  double t = 0;
  double tau_mean = 0;

  deque<particle>::iterator particle_1 = particles.begin();
  deque<particle>::iterator particle_2 = particles.begin();

  particle_2++;

  while (t <= dt) {
    double r = rand() * rmt;

    if (r <= 0) r = 3.e-5;

    double tau = -frequency_t * log(r);
    t += tau;
    if (t > dt) break;

    point v1 = particle_1->get_velocity();
    point v2 = particle_2->get_velocity();

    point velocity_sum = v1 + v2;

    point g = v2 - v1;
    double gmod = g.mod();

    if (g_max < gmod) {
      g_max = gmod;
      t -= tau;
      frequency_t = factor / g_max;
      continue;
    }

    double rr = rand() * rmt;

    if (gmod / g_max > rr)  // g == gmod?
    {
      double r1 = rand() * rmt;
      double r2 = rand() * rmt;

      double Fi = 2. * Pi * r1;
      double Eta = 2. * acos(r2);

      double gx = g.x;
      double gy = g.y;
      double gz = g.z;
      double gxz = sqrt(gx * gx + gz * gz);

      point g1;

      //                              point g1(gmod * sin(Pi * r1) * cos(2. * Pi
      //                              * r2),
      //                                       gmod * sin(Pi * r1) * sin(2. * Pi
      //                                       * r2), gmod * cos(Pi * r1));

      if (gxz > 1.e-6) {
        g1.x = gx * cos(Eta) -
               sin(Eta) * (gmod * gz * cos(Fi) - gx * gy * sin(Fi)) / gxz;
        g1.y = gy * cos(Eta) - gxz * sin(Fi) * sin(Eta);
        g1.z = gz * cos(Eta) +
               sin(Eta) * (gmod * gx * cos(Fi) + gy * gz * sin(Fi)) / gxz;
      } else {
        g1.x = -sin(Eta) * gmod * (cos(Fi) - sin(Fi)) / std::sqrt(2.);
        g1.y = gmod * cos(Eta);
        g1.z = sin(Eta) * gmod * (cos(Fi) + sin(Fi)) / std::sqrt(2.);
      }

      v1 = (velocity_sum - g1) * 0.5;
      v2 = (velocity_sum + g1) * 0.5;

      particle_1->velocity = v1;
      particle_2->velocity = v2;
    }

    particle_2++;

    if (particle_2 == particles.end()) {
      std::random_shuffle(particles.begin(), particles.end());

      particle_1 = particles.begin();
      particle_2 = particles.begin();
      particle_2++;
    } else {
      particle_1 = particle_2;
      particle_2++;
      if (particle_2 == particles.end()) {
        std::random_shuffle(particles.begin(), particles.end());

        particle_1 = particles.begin();
        particle_2 = particles.begin();
        particle_2++;
      }
    }
  }
}

// void cell::collisions()
//{
//	if(particles.size() <= 10) return;

//	//cout << calc_vel() << endl;

//	calc_Kn();

//	double rmt = 1 / double(RAND_MAX);
//	double g_max = 2 * sqrt(calc_T());
//	double factor = 2 * sqrt(2.) * L * Kn_l / particles.size();
//	double frequency_t = factor / g_max;

//	double t = 0;
//	double tau_mean = 0;

//	deque<particle>::iterator particle_1 = particles.begin();
//	deque<particle>::iterator particle_2 = particles.begin();

//	particle_2++;

//	while(t <= dt)
//	{
//		double r = rand() * rmt;

//		if(r <= 0) r = 0.00001;

//		double tau = -frequency_t * log(r);
//		t += tau;
//		if(t > dt) break;

//		point velocity_sum = particle_1->get_velocity() +
// particle_2->get_velocity();

//		point v1 = particle_1->get_velocity();
//		point v2 = particle_2->get_velocity();

//		velocity_sum = v1 + v2;

//		double g = (v2 - v1).mod();

//		if(g_max < g)
//		{
//			g_max = g;
//			t -=tau;
//			frequency_t = factor / g_max;
//			continue;
//		}

//		double rr = std::rand() * rmt;
//
//		if(g / g_max > rr)
//		{
//			double r1 = double(std::rand()) * rmt;
//			double r2 = double(std::rand()) * rmt;

//			point g1(	g * sin(Pi * r1) * cos(2. * Pi * r2),
//						g * sin(Pi * r1) * sin(2. * Pi *
// r2), 						g * cos(Pi * r1));

//			v1 = (velocity_sum - g1) * 0.5;
//			v2 = (velocity_sum + g1) * 0.5;

//			particle_1->velocity = v1;
//			particle_2->velocity = v2;
//		}
//
//		particle_2++;
//		if(particle_2 == particles.end())
//		{
//			std::random_shuffle(particles.begin(), particles.end());

//			particle_1 = particles.begin();
//			particle_2 = particles.begin();
//			particle_2++;
//		}
//		else
//		{
//			particle_1 = particle_2;
//			particle_2++;
//			if(particle_2 == particles.end())
//			{
//				std::random_shuffle(particles.begin(),
// particles.end());

//				particle_1 = particles.begin();
//				particle_2 = particles.begin();
//				particle_2++;
//			}
//		}
//	}
//}

// void cell::collisions()
//{
//	point temp_vel = calc_vel();
//	calc_Kn();
//	double rmt = 1 / double(RAND_MAX);
//	double g_max = 2 * sqrt(T);
//	double factor = 2 * sqrt(2.) * L * Kn_l / particles.size();
//	double frequency_t = factor / g_max;

//	double t = 0;
//	double tau_mean = 0;

//	deque<particle>::iterator a, b;
//	const deque<particle>::iterator end = particles.end();

//	if(particles.size() <= 1) return;

//	a = particles.begin();
//	b = particles.begin();
//	b++;

//	while(t <= dt)
//	{
//		try
//		{
//			if(a == particles.end())
//			{
//				std::random_shuffle(particles.begin(),
// particles.end()); 				a = particles.begin();
// b = particles.begin(); 				b++;
//			}
//
//			if(b == particles.end())
//			{
//				std::random_shuffle(particles.begin(),
// particles.end()); 				a = particles.begin();
// b = particles.begin(); 				b++;
//			}
//
//			double r = std::rand() * rmt;

//			if(r <= 0) r = 0.00001;

//			double tau = -frequency_t * log(r);
//			t += tau;
//			if(t > dt) break;

//			point vel_1;
//			point vel_2;
//
//			vel_1 = a->get_velocity();
//			vel_2 = b->get_velocity();
//
//			double g = (vel_2 - vel_1).mod();

//			if(g_max < g)
//			{
//				g_max = g;
//				t -=tau;
//				frequency_t = factor / g_max;
//			}

//			double rr = std::rand() * rmt;
//
//			if(g / g_max > rr)
//			{
//				double r1 = double(std::rand()) * rmt;
//				double r2 = double(std::rand()) * rmt;

//				point g1(	g * sin(Pi * r1) * cos(2. * Pi *
// r2), 							g * sin(Pi * r1) * sin(2. * Pi
// * r2), g * cos(Pi * r1));

//				double g1x = g * sin(Pi * r1) * cos(2. * Pi *
// r2); 				double g1y = g * sin(Pi * r1) * sin(2. * Pi
// * r2); double g1z = g * cos(Pi * r1);

//				//cout << "g1 = " << point(g1x, g1y, g1z) <<
// endl;

//				/*a->velocity.x = 0.5 * (vel_1.x + vel_2.x -
// g1.x); 				a->velocity.y = 0.5 * (vel_1.y + vel_2.y
// - g1.y); a->velocity.z = 0.5 * (vel_1.z + vel_2.z - g1.z);

//				b->velocity.x = 0.5 * (vel_1.x + vel_2.x +
// g1.x); 				b->velocity.y = 0.5 * (vel_1.y + vel_2.y
// + g1.y); b->velocity.z = 0.5 * (vel_1.z + vel_2.z + g1.z);*/
//				//cout << (vel_1 + vel_2 - g1) * .5 << endl <<
//(vel_1
//+ vel_2 + g1) * .5 << endl;

//				point v1 = (vel_1 + vel_2 - g1) * .5, v2 =
//(vel_1
//+ vel_2 + g1) * .5; 				a->velocity = v1;
// b->velocity = v2;

//				//vel_1 /= 2;
//				//vel_2 /= 2;

//				//cout << a->velocity << endl << b->velocity <<
// endl;

//				vel_1 /= 2;

//				/*a->velocity.x *= .5;
//				a->velocity.y *= .5;
//				a->velocity.z *= .5;
//				b->velocity.x *= .5;
//				b->velocity.y *= .5;
//				b->velocity.z *= .5;*/
//			}

//			a++;
//			if(a == particles.end())
//			{
//				std::random_shuffle(particles.begin(),
// particles.end());

//				a = particles.begin();
//				b = particles.begin();
//				b++;
//			}
//			else
//			{
//				b = a;
//				b++;
//			}
//			if(b == particles.end())
//			{
//				std::random_shuffle(particles.begin(),
// particles.end());

//				a = particles.begin();
//				b = particles.begin();
//				b++;
//			}
//		}
//		catch(...)
//		{
//			std::cerr << "Exception in cell::collisions()!" << endl;
//			exit(mc3d::unusual_situations::FATAL_ERROR);
//		}
//	}
//	temp_vel -= calc_vel();
//	//if(temp_vel.mod() > 0.0001)cout << temp_vel << endl;
//}

void cell::particle_move() {
  deque<particle>::iterator a;
  a = particles.begin();
  velocity = point(0, 0, 0);
  while (a != particles.end()) {
    /*a->x += a->u * dt;
    a->y += a->v * dt;
    a->z += a->w * dt;*/

    a->position += a->velocity * dt;

    a++;
  }
}

void cell::sort(void) {
  /*if(body)
  {
          particle_buffer.insert(particle_buffer.end(), particles.begin(),
  particles.end()); particles.clear(); return;
  }*/

  deque<particle>::iterator data = particles.begin();

  while (data != particles.end()) {
    if ((data->position.x < apex.x) || (data->position.x > apex.x + lx) ||
        (data->position.y < apex.y) || (data->position.y > apex.y + ly) ||
        (data->position.z < apex.z) || (data->position.z > apex.z + lz)) {
      particle_buffer.insert(particle_buffer.end(), data, data + 1);
      data = particles.erase(data);
    } else
      data++;
  }
}

void cell::set_L(double L) { this->L = L; }

void cell::add_particle(deque<particle>* a) {
  deque<particle>::iterator temp;
  temp = a->begin();
  while (temp != a->end()) {
    if ((temp->position.x > apex.x) && (temp->position.x < apex.x + lx) &&
        (temp->position.y > apex.y) && (temp->position.y < apex.y + ly) &&
        (temp->position.z > apex.z) && (temp->position.z < apex.z + lz)) {
      particles.insert(particles.end(), temp, temp + 1);
      temp = a->erase(temp);
    } else
      temp++;
  }
}

void cell::add_particle(particle* a, int n) {
  for (int i = 0; i < n; i++) {
    if ((a[i].position.x > apex.x) && (a[i].position.x < apex.x + lx) &&
        (a[i].position.y > apex.y) && (a[i].position.y < apex.y + ly) &&
        (a[i].position.z > apex.z) && (a[i].position.z < apex.z + lz)) {
      particles.push_back(a[i]);
    }
  }
}

deque<particle>* cell::get_buffer() { return &particle_buffer; }

/*var cell::get_var()
{
        var a(velocity.x, velocity.y, velocity.z, T, n);
        return a;
}*/

point cell::get_apex() const { return apex; }

point cell::get_center() const { return apex + this->get_size() / 2; }

point cell::get_mass_center() const { return mass_center; }

double cell::get_dt() {
  this->calc_dt();
  return dt;
}

void cell::set_dt(double dt) { this->dt = dt; }

double cell::t() { return T; }

void cell::set_param(double S, double alpha, double T) {
  this->T = T;
  velocity =
      point(S * sqrt(2 * T) * cos(alpha), S * sqrt(2 * T) * sin(alpha), 0);
}

void cell::set_vel(point velocity) { this->velocity = velocity; }

void cell::set_t(double t) { this->T = t; }

double cell::get_Kn() { return Kn; }

double cell::calc_Kn() {
  Kn_l = Kn * np / particles.size();
  return Kn;
}

double cell::get_t() {
  calc_T();
  return T;
}

double cell::get_u() {
  deque<particle>::iterator data = particles.begin();
  // double av_vel = 0;
  velocity = point(0, 0, 0);
  while (data != particles.end()) {
    velocity += data->velocity;
    data++;
  }
  velocity /= double(particles.size());
  return velocity.x;
}

double cell::get_v() {
  deque<particle>::iterator data = particles.begin();
  while (data != particles.end()) {
    velocity += data->velocity;
    data++;
  }
  velocity /= double(particles.size());
  return velocity.y;
}

double cell::get_w() {
  deque<particle>::iterator data = particles.begin();
  while (data != particles.end()) {
    velocity += data->velocity;
    data++;
  }
  velocity /= double(particles.size());
  return velocity.z;
}

double cell::get_energy() {
  double E = 0;
  deque<particle>::iterator data = particles.begin();
  while (data != particles.end()) {
    /*E += data->u * data->u + data->w * data->w + data->v * data->v;*/
    E += data->velocity * data->velocity;
    data++;
  }
  this->n = particles.size();
  E /= 2 * n;
  return E;
}

point cell::calc_vel() {
  deque<particle>::iterator data = particles.begin();

  velocity = point(0, 0, 0);
  while (data != particles.end()) {
    velocity += data->velocity;
    data++;
  }

  n = particles.size();

  velocity /= double(particles.size());

  return velocity;
}

double cell::calc_T() {
  deque<particle>::iterator data = particles.begin();
  double E = 0;
  point av_vel(0, 0, 0);
  while (data != particles.end()) {
    av_vel += data->velocity;
    E += data->velocity * data->velocity;
    data++;
  }
  this->n = particles.size();
  av_vel /= double(n);
  E /= double(n);
  velocity = av_vel;
  T = (E - velocity * velocity) / 3;
  return T;
}

double cell::calc_dt() {
  calc_T();
  dt = 1000000;
  double c = sqrt(2 * T);
  double dtt = min(lx / (fabs(velocity.x) + c), ly / fabs(velocity.y) + c);
  dtt = min(dtt, lz / (fabs(velocity.z) + c));
  if (dtt < dt) dt = dtt;
  // if(dt < 0.00001) cin>>c;
  return dt;
}

void cell::set_Kn(double Kn) { this->Kn = Kn; }

void cell::add_neighbor(cell* neighbor) { this->neighbors.push_back(neighbor); }
cell* cell::get_ptr() { return this; }

void cell::neighbor_sort() {
  deque<cell*>::iterator cell_iter = neighbors.begin();

  while (cell_iter != neighbors.end()) {
    cell* a = *cell_iter;
    a->add_particle(&particle_buffer);
    cell_iter++;
  }
}

void cell::calc() {
  //*thread_mark = 2;
  //		double start_time = MPI_Wtime();
  collisions();
  if (body_mark)
    body_boundary.bondary_condition(
        &particles, dt);  //при выполнении граничных условий на теле выполняется
                          //и передвижение частиц!!!
  else
    particle_move();
  //		calc_time += MPI_Wtime() - start_time;
  //*thread_mark = 1;
}

void cell::attach_thread_mark(int& ptr) { thread_mark = &ptr; }

point cell::get_size() const { return point(abs(lx), abs(ly), abs(lz)); }

double cell::get_L() { return L; }

deque<cell*> cell::fragmentation(const std::unique_ptr<geometry>& body) {
  deque<cell*> new_cells(8);

  for (size_t i = 0; i < 8; i++) {
    new_cells[i] = new cell;
  }

  point a(0, 0, 0);
  point dl(lx / 2, ly / 2, lz / 2);

  for (size_t i = 0; i < 8; i++) {
    new_cells[i]->set_size(dl);
    // new_cells[i]->set_inner_boundary(body_boundary);
  }

  a = apex;
  new_cells[0]->set_apex(a);
  // new_cells[0]->add_particle(&particles);

  point shift(lx / 2, 0, 0);

  a = apex + shift;
  new_cells[1]->set_apex(a);
  // new_cells[1]->add_particle(&particles);

  shift.set(0, ly / 2, 0);
  a = apex + shift;
  new_cells[2]->set_apex(a);
  // new_cells[2]->add_particle(&particles);

  shift.set(0, 0, lz / 2);
  a = apex + shift;
  new_cells[3]->set_apex(a);
  // new_cells[3]->add_particle(&particles);

  shift.set(lx / 2, ly / 2, 0);
  a = apex + shift;
  new_cells[4]->set_apex(a);
  // new_cells[4]->add_particle(&particles);

  shift.set(lx / 2, 0, lz / 2);
  a = apex + shift;
  new_cells[5]->set_apex(a);
  // new_cells[5]->add_particle(&particles);

  shift.set(0, ly / 2, lz / 2);
  a = apex + shift;
  new_cells[6]->set_apex(a);
  // new_cells[6]->add_particle(&particles);

  shift.set(lx / 2, ly / 2, lz / 2);
  a = apex + shift;
  new_cells[7]->set_apex(a);
  // new_cells[7]->add_particle(&particles);

  for (size_t i = 0; i < 8; i++) {
    new_cells[i]->initialazition(n / 8, body);
    // new_cells[i]->set_inner_boundary(body_boundary);
  }

  return new_cells;
}

void cell::set_inner_boundary(inner_boundary bound) { body_boundary = bound; }

double cell::calc_volume() {
  volume = body_boundary.calc_cell_volume(apex, get_size(), &mass_center);
  if (volume < 0) {
    return volume = get_size().volume();
  }
  return volume;
}

double cell::get_volume() const { return volume; }

bool cell::get_body_mark() { return body_mark; }

point cell::get_velocity() {
  deque<particle>::iterator data = particles.begin();
  velocity = point(0, 0, 0);
  while (data != particles.end()) {
    velocity += data->velocity;
    data++;
  }
  this->n = particles.size();
  velocity /= double(particles.size());
  return velocity;
}

point cell::get_particle_mass_center() {
  deque<particle>::iterator data = particles.begin();
  point particle_mass_center(0, 0, 0);
  while (data != particles.end()) {
    particle_mass_center += data->position;
    data++;
  }
  this->n = particles.size();
  particle_mass_center /= double(particles.size());
  return particle_mass_center;
}

void cell::clean_inner_particle(const geometry& body) {
  if (particles.size() < 100) {
    point ad = point(0, 0, 0);
  }
  this->particle_buffer.swap(particles);

  deque<particle>::iterator data = particle_buffer.begin();
  while (data != particle_buffer.end()) {
    if (!body.is_inner_point(data->position)) {
      particles.push_back(*data);
    }
    data++;
  }

  particle_buffer.clear();
}

bool operator>(const cell& a, const cell& b) {
  return a.calc_time > b.calc_time;
}

bool operator<(const cell& a, const cell& b) {
  return b.calc_time > a.calc_time;
}

void cell::set_body_mark(bool mark) { this->body_mark = mark; }

bool cell::_dbg_test_particle() {
  deque<particle>::iterator data = particles.begin();
  int i = 0;

  while (data != particles.end()) {
    if (this->body_boundary.get_geometry_ptr()->is_inner_point(
            data->get_position())) {
      cout << data->get_position() << "\t" << data->get_velocity() << endl;
      this->body_boundary.get_geometry_ptr()->is_inner_point(
          data->get_position());
      i++;
    }
    data++;
  }
  if (i > 0)
    return true;
  else
    return false;
}
}  // namespace mc3d
