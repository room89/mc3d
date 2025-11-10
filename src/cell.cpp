#include "cell.h"

#include <algorithm>
#include <random>
#include <utils/logger.hpp>
#include <utils/utils.hpp>

namespace mc3d {
namespace {
const double Pi = 3.14159265358979;
}

Cell::Cell() {
  this->calc_time = 0;
  //		this->start_time = MPI_Wtime();
  np = 1;
  body_mark = false;
}

Cell::~Cell() {
  particles.clear();
  particle_buffer.clear();
  neighbors.clear();
}

void Cell::SetSize(double lx, double ly, double lz) {
  if ((lx <= 0) || (ly <= 0) || (lz <= 0)) return;
  this->lx = lx;
  this->ly = ly;
  this->lz = lz;
  L = lx > ly ? lx : ly;
  L = L > lz ? L : lz;
}

void Cell::SetSize(Point dl) {
  if ((dl.x <= 0) || (dl.y <= 0) || (dl.z <= 0)) return;
  lx = dl.x;
  ly = dl.y;
  lz = dl.z;
  L = lx > ly ? lx : ly;
  L = L > lz ? L : lz;
}

void Cell::SetApex(Point a) { apex = a; }

bool Cell::Initialize(size_t N, const std::unique_ptr<Geometry>& bbody) {
  this->GenerateRandom(N);

  dt = 100000;
  CalculateTemperature();
  double c = sqrt(2 * T);
  double dtt = lx / (fabs(velocity.x) + c) < ly / (fabs(velocity.y) + c)
                   ? lx / (fabs(velocity.x) + c)
                   : ly / (fabs(velocity.y) + c);
  dtt = dtt < lz / (fabs(velocity.z) + c) ? dtt : lz / (fabs(velocity.z) + c);
  this->dt = dtt;
  static std::mt19937 rng(std::random_device{}());
  std::shuffle(particles.begin(), particles.end(), rng);

  if (bbody) {
    body_mark = body_boundary.AddPolygon(bbody.get(), GetCenter(), L);
    body_boundary.SetGeometry(bbody.get());
    //    if(body_boundary.Empty())
    //      particles.clear();
  } else
    body_mark = false;

  return true;
}

unsigned int Cell::GetParticleCount() const { return particles.size(); }

double Cell::GenerateRandom(size_t N) {
  Point d(lx, ly, lz);

  const double rmt = 1.0 / static_cast<double>(RAND_MAX);
  double ti = 0;
  deque<Particle>::iterator data;

  size_t nn = N;
  double nt = 1 / double(N);

  data = particles.begin();

  size_t nn2 = (nn - nn % 2) / 2;

  Particle p1, p2;

  for (unsigned int i = 0; i < nn2; i++) {
    double rn1 = double(std::rand()) * rmt;
    double rn2 = double(std::rand()) * rmt;

    if (rn1 <= 0) rn1 = 0.00000001;

    double slg = sqrt(2 * fabs(log(rn1)));

    /*p1.u = slg * cos(2 * Pi * rn2);
    p1.v = slg * sin(2 * Pi * rn2);*/

    p1.SetU(slg * cos(2 * Pi * rn2));
    p1.SetV(slg * sin(2 * Pi * rn2));

    /*p2.u = -slg * cos(2 * Pi * rn2);
    p2.v = -slg * sin(2 * Pi * rn2);*/

    p2.SetU(-slg * cos(2 * Pi * rn2));
    p2.SetV(-slg * sin(2 * Pi * rn2));

    double rn3 = std::rand() * rmt;
    double rn4 = std::rand() * rmt;

    if (rn3 <= 0) rn3 = 0.00000001;

    slg = sqrt(2 * fabs(log(rn3)));

    /*p1.w = slg * cos(2 * Pi * rn4);
    p2.w = -slg * cos(2 * Pi * rn4);*/

    p1.SetW(slg * cos(2 * Pi * rn4));
    p2.SetW(-slg * cos(2 * Pi * rn4));

    particles.push_front(p1);
    particles.push_front(p2);

    /*ti += 2 * (p1.u * p1.u + p1.v * p1.v + p1.w * p1.w);*/

    ti += 2 * (p1.GetVelocity() * p1.GetVelocity());
  }

  data = particles.begin();

  if (N % 2 == 1) {
    // p1.u = p1.v = p1.w = 0;
    p1.SetVelocity(Point(0, 0, 0));
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

    Point pos(d.x * rnx, d.y * rny, d.z * rnz);

    data->position = apex + pos;

    /*data->x = apex.x + d.x * rnx;
    data->y = apex.y + d.y * rny;
    data->z = apex.z + d.z * rnz;*/

    data++;
  }

  CalculateTemperature();

  return 0;
}

double Cell::GenerateRandom(unsigned int N, double T, Point V) {
  Point d(lx, ly, lz);

  const double rmt = 1.0 / static_cast<double>(RAND_MAX);
  double ti = 0;
  deque<Particle>::iterator data;

  unsigned int nn = N;
  double nt = 1 / double(N);

  data = particles.begin();

  unsigned int nn2 = (nn - nn % 2) / 2;

  Particle p1, p2;

  for (unsigned int i = 0; i < nn2; i++) {
    double rn1 = double(std::rand()) * rmt;
    double rn2 = double(std::rand()) * rmt;

    if (rn1 <= 0) rn1 = 0.00001;

    double slg = sqrt(2 * fabs(log(rn1)));

    /*p1.u = slg * cos(2 * Pi * rn2);
    p1.v = slg * sin(2 * Pi * rn2);*/

    p1.SetU(slg * cos(2 * Pi * rn2));
    p1.SetV(slg * sin(2 * Pi * rn2));

    /*p2.u = -slg * cos(2 * Pi * rn2);
    p2.v = -slg * sin(2 * Pi * rn2);*/

    p2.SetU(-slg * cos(2 * Pi * rn2));
    p2.SetV(-slg * sin(2 * Pi * rn2));

    double rn3 = std::rand() * rmt;
    double rn4 = std::rand() * rmt;

    if (rn3 <= 0) rn3 = 0.00001;

    slg = sqrt(2 * fabs(log(rn3)));

    /*p1.w = slg * cos(2 * Pi * rn4);
    p2.w = -slg * cos(2 * Pi * rn4);*/

    p1.SetW(slg * cos(2 * Pi * rn4));
    p2.SetW(-slg * cos(2 * Pi * rn4));

    particles.push_front(p1);
    particles.push_front(p2);

    // ti += 2 * (p1.u * p1.u + p1.v * p1.v + p1.w * p1.w);

    ti += 2 * (p1.GetVelocity() * p1.GetVelocity());
  }

  data = particles.begin();

  if (N % 2 == 1) {
    // p1.u = p1.v = p1.w = 0;
    p1.SetVelocity(Point(0, 0, 0));
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

    data->position = apex + Point(d.x * rnx, d.y * rny, d.z * rnz);

    data++;
  }

  static std::mt19937 rng(std::random_device{}());
  std::shuffle(particles.begin(), particles.end(), rng);

  return 0;
}

double Cell::GenerateFreeRandom(unsigned int N, double T, Point V, Point nrml) {
  Point d(lx, ly, lz);

  const double rmt = 1.0 / static_cast<double>(RAND_MAX);
  double ti = 0;
  deque<Particle>::iterator data;
  deque<Particle> added_particles;

  double nt = 1 / double(N);

  data = added_particles.begin();
  Particle p1;
  Point dvel(0, 0, 0);

  unsigned int i = 0;

  while (i < N) {
    double rn1 = rand() * rmt;
    double rn2 = rand() * rmt;

    if (rn1 <= 0) rn1 = 0.00001;

    double slg = sqrt(2 * T * fabs(log(rn1)));

    p1.SetU(slg * cos(2 * Pi * rn2));
    p1.SetV(slg * sin(2 * Pi * rn2));

    double rn3 = rand() * rmt;
    double rn4 = rand() * rmt;

    if (rn3 <= 0) rn3 = 0.00001;

    slg = sqrt(2 * T * fabs(log(rn3)));

    p1.SetW(slg * cos(2 * Pi * rn4));

    if ((p1.GetVelocity() + V) * nrml < 0) {
      added_particles.push_front(p1);
      ti += (p1.velocity * p1.velocity);
      dvel += p1.velocity;
      i++;
    }
  }

  dvel *= nt;
  ti = (ti * nt - dvel * dvel) / 3;

  double sf = sqrt(T / ti);
  //		sf = 1;

  data = added_particles.begin();

  for (i = 0; i < N; i++) {
    if (data == added_particles.end()) break;
    data->velocity = (data->velocity - dvel) * sf + V;
    if (data->velocity * nrml < 0)
      data->velocity = 2 * V - data->velocity;  // data++;
  }

  data = added_particles.begin();

  for (i = 0; i < N; i++) {
    if (data == added_particles.end()) break;

    data->velocity += V;  //

    double rnx = rand() * rmt;
    double rny = rand() * rmt;
    double rnz = rand() * rmt;

    data->position = apex + Point(d.x * rnx, d.y * rny, d.z * rnz);

    data++;
  }

  //              random_shuffle(added_particles.begin(),
  //              added_particles.end());

  particles.insert(particles.end(), added_particles.begin(),
                   added_particles.end());

  return 0;
}

/*double Cell::GenerateHyperFreeRandom(unsigned int N, Point V, double T)
{
         deque<Particle> added_particles;
        const double rmt = 1.0 / static_cast<double>(RAND_MAX);
         Point position;
         Point d(lx, ly, lz);

         for(unsigned int i = 0; i < N; i++)
         {
                double rnx = rand() * rmt;
                double rny = rand() * rmt;
                double rnz = rand() * rmt;

                position = apex + Point(d.x * rnx, d.y * rny, d.z * rnz);

                Particle new_particle;

                new_particle.position = position;
                Point noise(rand() * rmt - .5, rand() * rmt - .5, rand() * rmt -
.5); new_particle.velocity = V + sqrt(2 * T) * noise;

                added_particles.push_back(new_particle);
         }

         particles.insert(particles.end(), added_particles.begin(),
added_particles.end());

         return 0;
}*/

double Cell::GenerateHyperFreeRandom(unsigned int N, Point V, double T) {
  deque<Particle> added_particles;
  const double rmt = 1.0 / static_cast<double>(RAND_MAX);
  Point position;
  Point d(lx, ly, lz);
  double A = sqrt(3 * T);

  for (unsigned int i = 0; i < N; i++) {
    double rnx = utils::RandomDouble(0, 1);
    double rny = utils::RandomDouble(0, 1);
    double rnz = utils::RandomDouble(0, 1);

    position = apex + Point(d.x * rnx, d.y * rny, d.z * rnz);

    Particle new_particle;

    new_particle.position = position;
    Point noise(A * utils::RandomDouble(-1, 1), A * utils::RandomDouble(-1, 1),
                A * utils::RandomDouble(-1, 1));
    new_particle.velocity = V + noise;

    added_particles.push_back(new_particle);
  }

  Point av_velocity;
  for (const auto& Particle : added_particles) {
    av_velocity += Particle.velocity;
  }
  av_velocity /= particles.size();

  particles.insert(particles.end(), added_particles.begin(),
                   added_particles.end());

  return 0;
}

void Cell::Collisions() {
  if (particles.size() <= 5) return;

  // cout << CalculateVelocity() << endl;

  CalculateKn();

  const double rmt = 1.0 / static_cast<double>(RAND_MAX);
  double g_max = 2 * sqrt(CalculateTemperature());
  double factor = 2 * sqrt(2.) * L * Kn_l / particles.size();
  double frequency_t = factor / g_max;

  double t = 0;

  static std::mt19937 rng(std::random_device{}());

  deque<Particle>::iterator particle_1 = particles.begin();
  deque<Particle>::iterator particle_2 = particles.begin();

  particle_2++;

  while (t <= dt) {
    double r = rand() * rmt;

    if (r <= 0) r = 3.e-5;

    double tau = -frequency_t * log(r);
    t += tau;
    if (t > dt) break;

    Point v1 = particle_1->GetVelocity();
    Point v2 = particle_2->GetVelocity();

    Point velocity_sum = v1 + v2;

    Point g = v2 - v1;
    double gmod = g.Mod();

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

      Point g1;

      //                              Point g1(gmod * sin(Pi * r1) * cos(2. * Pi
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
      std::shuffle(particles.begin(), particles.end(), rng);

      particle_1 = particles.begin();
      particle_2 = particles.begin();
      particle_2++;
    } else {
      particle_1 = particle_2;
      particle_2++;
      if (particle_2 == particles.end()) {
        std::shuffle(particles.begin(), particles.end(), rng);

        particle_1 = particles.begin();
        particle_2 = particles.begin();
        particle_2++;
      }
    }
  }
}

// void Cell::Collisions()
//{
//	if(particles.size() <= 10) return;

//	//cout << CalculateVelocity() << endl;

//	CalculateKn();

//	double rmt = 1 / double(RAND_MAX);
//	double g_max = 2 * sqrt(CalculateTemperature());
//	double factor = 2 * sqrt(2.) * L * Kn_l / particles.size();
//	double frequency_t = factor / g_max;

//	double t = 0;
//	double tau_mean = 0;

//	deque<Particle>::iterator particle_1 = particles.begin();
//	deque<Particle>::iterator particle_2 = particles.begin();

//	particle_2++;

//	while(t <= dt)
//	{
//		double r = rand() * rmt;

//		if(r <= 0) r = 0.00001;

//		double tau = -frequency_t * log(r);
//		t += tau;
//		if(t > dt) break;

//		Point velocity_sum = particle_1->GetVelocity() +
// particle_2->GetVelocity();

//		Point v1 = particle_1->GetVelocity();
//		Point v2 = particle_2->GetVelocity();

//		velocity_sum = v1 + v2;

//		double g = (v2 - v1).Mod();

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

//			Point g1(	g * sin(Pi * r1) * cos(2. * Pi * r2),
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

// void Cell::Collisions()
//{
//	Point temp_vel = CalculateVelocity();
//	CalculateKn();
//	double rmt = 1 / double(RAND_MAX);
//	double g_max = 2 * sqrt(T);
//	double factor = 2 * sqrt(2.) * L * Kn_l / particles.size();
//	double frequency_t = factor / g_max;

//	double t = 0;
//	double tau_mean = 0;

//	deque<Particle>::iterator a, b;
//	const deque<Particle>::iterator end = particles.end();

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

//			Point vel_1;
//			Point vel_2;
//
//			vel_1 = a->GetVelocity();
//			vel_2 = b->GetVelocity();
//
//			double g = (vel_2 - vel_1).Mod();

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

//				Point g1(	g * sin(Pi * r1) * cos(2. * Pi *
// r2), 							g * sin(Pi * r1)
// * sin(2.
// * Pi
// * r2), g * cos(Pi * r1));

//				double g1x = g * sin(Pi * r1) * cos(2. * Pi *
// r2); 				double g1y = g * sin(Pi * r1) * sin(2. *
// Pi
// * r2); double g1z = g * cos(Pi * r1);

//				//cout << "g1 = " << Point(g1x, g1y, g1z) <<
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

//				Point v1 = (vel_1 + vel_2 - g1) * .5, v2 =
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
//			std::cerr << "Exception in Cell::Collisions()!" << endl;
//			exit(mc3d::unusual_situations::FATAL_ERROR);
//		}
//	}
//	temp_vel -= CalculateVelocity();
//	//if(temp_vel.Mod() > 0.0001)cout << temp_vel << endl;
//}

void Cell::MoveParticles() {
  deque<Particle>::iterator a;
  a = particles.begin();
  velocity = Point(0, 0, 0);
  while (a != particles.end()) {
    /*a->x += a->u * dt;
    a->y += a->v * dt;
    a->z += a->w * dt;*/

    a->position += a->velocity * dt;

    a++;
  }
}

void Cell::Sort() {
  /*if(body)
  {
          particle_buffer.insert(particle_buffer.end(), particles.begin(),
  particles.end()); particles.clear(); return;
  }*/

  deque<Particle>::iterator data = particles.begin();

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

void Cell::SetCharacteristicLength(double L) { this->L = L; }

void Cell::AddParticle(deque<Particle>* a) {
  deque<Particle>::iterator temp;
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

void Cell::AddParticle(Particle* a, int n) {
  for (int i = 0; i < n; i++) {
    if ((a[i].position.x > apex.x) && (a[i].position.x < apex.x + lx) &&
        (a[i].position.y > apex.y) && (a[i].position.y < apex.y + ly) &&
        (a[i].position.z > apex.z) && (a[i].position.z < apex.z + lz)) {
      particles.push_back(a[i]);
    }
  }
}

deque<Particle>* Cell::GetBuffer() { return &particle_buffer; }

/*var Cell::get_var()
{
        var a(velocity.x, velocity.y, velocity.z, T, n);
        return a;
}*/

Point Cell::GetApex() const { return apex; }

Point Cell::GetCenter() const { return apex + this->GetSize() / 2; }

Point Cell::GetMassCenter() const { return mass_center; }

double Cell::GetDt() {
  this->CalculateDt();
  return dt;
}

void Cell::SetDt(double dt) { this->dt = dt; }

double Cell::GetTemperatureRaw() { return T; }

void Cell::SetParameters(double S, double alpha, double T) {
  this->T = T;
  velocity =
      Point(S * sqrt(2 * T) * cos(alpha), S * sqrt(2 * T) * sin(alpha), 0);
}

void Cell::SetVelocity(Point velocity) { this->velocity = velocity; }

void Cell::SetTemperature(double t) { this->T = t; }

double Cell::GetKn() { return Kn; }

double Cell::CalculateKn() {
  Kn_l = Kn * np / particles.size();
  return Kn;
}

double Cell::GetTemperature() {
  CalculateTemperature();
  return T;
}

double Cell::GetU() {
  deque<Particle>::iterator data = particles.begin();
  // double av_vel = 0;
  velocity = Point(0, 0, 0);
  while (data != particles.end()) {
    velocity += data->velocity;
    data++;
  }
  velocity /= double(particles.size());
  return velocity.x;
}

double Cell::GetV() {
  deque<Particle>::iterator data = particles.begin();
  while (data != particles.end()) {
    velocity += data->velocity;
    data++;
  }
  velocity /= double(particles.size());
  return velocity.y;
}

double Cell::GetW() {
  deque<Particle>::iterator data = particles.begin();
  while (data != particles.end()) {
    velocity += data->velocity;
    data++;
  }
  velocity /= double(particles.size());
  return velocity.z;
}

double Cell::GetEnergy() {
  double E = 0;
  deque<Particle>::iterator data = particles.begin();
  while (data != particles.end()) {
    /*E += data->u * data->u + data->w * data->w + data->v * data->v;*/
    E += data->velocity * data->velocity;
    data++;
  }
  E /= 2 * particles.size();
  return E;
}

Point Cell::CalculateVelocity() const {
  auto velocity = Point(0, 0, 0);
  for (const auto& Particle : particles) {
    velocity += Particle.GetVelocity();
  }

  velocity /= double(particles.size());

  return velocity;
}

double Cell::CalculateTemperature() {
  deque<Particle>::iterator data = particles.begin();
  double E = 0;
  Point av_vel(0, 0, 0);
  while (data != particles.end()) {
    av_vel += data->velocity;
    E += data->velocity * data->velocity;
    data++;
  }
  av_vel /= double(particles.size());
  E /= double(particles.size());
  velocity = av_vel;
  T = (E - velocity * velocity) / 3;
  return T;
}

double Cell::CalculateDt() {
  CalculateTemperature();
  dt = 1000000;
  double c = sqrt(2 * T);
  double dtt = min(lx / (fabs(velocity.x) + c), ly / fabs(velocity.y) + c);
  dtt = min(dtt, lz / (fabs(velocity.z) + c));
  if (dtt < dt) dt = dtt;
  // if(dt < 0.00001) cin>>c;
  return dt;
}

void Cell::SetKn(double Kn) { this->Kn = Kn; }

void Cell::AddNeighbor(Cell* neighbor) { this->neighbors.push_back(neighbor); }
Cell* Cell::GetPtr() { return this; }

void Cell::SortNeighbors() {
  deque<Cell*>::iterator cell_iter = neighbors.begin();

  while (cell_iter != neighbors.end()) {
    Cell* a = *cell_iter;
    a->AddParticle(&particle_buffer);
    cell_iter++;
  }
}

void Cell::Calculate() {
  Collisions();
  if (body_mark)
    body_boundary.BoundaryCondition(&particles, dt);
  else
    MoveParticles();
}

void Cell::AttachThreadMark(int& ptr) { thread_mark = &ptr; }

Point Cell::GetSize() const { return Point(abs(lx), abs(ly), abs(lz)); }

double Cell::GetCharacteristicLength() { return L; }

deque<Cell*> Cell::Fragment(const std::unique_ptr<Geometry>& body) {
  deque<Cell*> new_cells(8);

  for (size_t i = 0; i < 8; i++) {
    new_cells[i] = new Cell;
  }

  Point a(0, 0, 0);
  Point dl(lx / 2, ly / 2, lz / 2);

  for (size_t i = 0; i < 8; i++) {
    new_cells[i]->SetSize(dl);
    // new_cells[i]->SetInnerBoundary(body_boundary);
  }

  a = apex;
  new_cells[0]->SetApex(a);
  // new_cells[0]->AddParticle(&particles);

  Point shift(lx / 2, 0, 0);

  a = apex + shift;
  new_cells[1]->SetApex(a);
  // new_cells[1]->AddParticle(&particles);

  shift.Set(0, ly / 2, 0);
  a = apex + shift;
  new_cells[2]->SetApex(a);
  // new_cells[2]->AddParticle(&particles);

  shift.Set(0, 0, lz / 2);
  a = apex + shift;
  new_cells[3]->SetApex(a);
  // new_cells[3]->AddParticle(&particles);

  shift.Set(lx / 2, ly / 2, 0);
  a = apex + shift;
  new_cells[4]->SetApex(a);
  // new_cells[4]->AddParticle(&particles);

  shift.Set(lx / 2, 0, lz / 2);
  a = apex + shift;
  new_cells[5]->SetApex(a);
  // new_cells[5]->AddParticle(&particles);

  shift.Set(0, ly / 2, lz / 2);
  a = apex + shift;
  new_cells[6]->SetApex(a);
  // new_cells[6]->AddParticle(&particles);

  shift.Set(lx / 2, ly / 2, lz / 2);
  a = apex + shift;
  new_cells[7]->SetApex(a);
  // new_cells[7]->AddParticle(&particles);

  for (size_t i = 0; i < 8; i++) {
    new_cells[i]->Initialize(particles.size() / 8, body);
    // new_cells[i]->SetInnerBoundary(body_boundary);
  }

  return new_cells;
}

void Cell::SetInnerBoundary(InnerBoundary bound) { body_boundary = bound; }

double Cell::CalculateVolume() {
  volume_ = body_boundary.CalcCellVolume(apex, GetSize(), &mass_center);
  if (volume_ < 0) {
    return volume_ = GetSize().Volume();
  }
  return volume_;
}

double Cell::GetVolume() const { return volume_; }

bool Cell::GetBodyMark() { return body_mark; }

Point Cell::GetVelocity() {
  deque<Particle>::iterator data = particles.begin();
  velocity = Point(0, 0, 0);
  while (data != particles.end()) {
    velocity += data->velocity;
    data++;
  }
  velocity /= double(particles.size());
  return velocity;
}

Point Cell::GetParticleMassCenter() {
  deque<Particle>::iterator data = particles.begin();
  Point particle_mass_center(0, 0, 0);
  while (data != particles.end()) {
    particle_mass_center += data->position;
    data++;
  }
  particle_mass_center /= double(particles.size());
  return particle_mass_center;
}

void Cell::CleanInnerParticles(const Geometry& body) {
  if (particles.size() < 100) {
    Point ad = Point(0, 0, 0);
  }
  this->particle_buffer.swap(particles);

  deque<Particle>::iterator data = particle_buffer.begin();
  while (data != particle_buffer.end()) {
    if (!body.IsInnerPoint(data->position)) {
      particles.push_back(*data);
    }
    data++;
  }

  particle_buffer.clear();
}

bool operator>(const Cell& a, const Cell& b) {
  return a.calc_time > b.calc_time;
}

bool operator<(const Cell& a, const Cell& b) {
  return b.calc_time > a.calc_time;
}

void Cell::SetBodyMark(bool mark) { this->body_mark = mark; }

bool Cell::DebugTestParticle() {
  deque<Particle>::iterator data = particles.begin();
  int i = 0;

  while (data != particles.end()) {
    if (this->body_boundary.GetGeometryPtr()->IsInnerPoint(
            data->GetPosition())) {
      LOG_DEBUG() << data->GetPosition() << "\t" << data->GetVelocity();
      this->body_boundary.GetGeometryPtr()->IsInnerPoint(data->GetPosition());
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
