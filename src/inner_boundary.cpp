#include "inner_boundary.h"

namespace mc3d {
const double Pi = 3.1415926535;

InnerBoundary::InnerBoundary() {
  Tw = 1;
  eps = 0.00001;
}

InnerBoundary::~InnerBoundary() {}

bool InnerBoundary::AddPolygon(Geometry* body, Point cell_center, double L) {
  /*Tw = 1;
  eps = 0.00001;

  this->cell_center = cell_center;

  deque<Polygon*>::iterator poligon_iter = bbody->poligons.begin();

  while(poligon_iter != bbody->poligons.end())
  {
          if(((*poligon_iter)->GetP1() - cell_center) *
  (*poligon_iter)->GetNormal() < L)
  //проверка растояния до плоскости полигона
          {
                  if( (((*poligon_iter)->GetGmt() - cell_center).Mod() < 2 * L)
  ||
  //проверка растояния до ГМТ полигона(полигон - треугольник)
                          (((*poligon_iter)->GetP1()  - cell_center).Mod() < 2
  * L) ||
  //проверка растояния до первого угла полигона
                          (((*poligon_iter)->GetP2()  - cell_center).Mod() < 2
  * L) ||
  //проверка растояния до второго угла полигона
                          (((*poligon_iter)->GetP3()  - cell_center).Mod() < 2
  * L))
  //проверка растояния до третьего угла полигона
                  {
                                  poligon_ptrs.push_back((*poligon_iter)->GetPtr());
  //добавление указателя на полигон
                  }
          }

          poligon_iter++;
  }

  if(poligon_ptrs.size() > 0) return true;
  else return false;*/
  this->cell_center = cell_center;

  deque<Polygon*>::iterator poligon_iter = body->poligons.begin();

  while (poligon_iter != body->poligons.end()) {
    if ((((*poligon_iter)->GetGmt() - cell_center).Mod() <
         L) ||  // проверка растояния до ГМТ полигона(полигон - треугольник)
        (((*poligon_iter)->GetP1() - cell_center).Mod() <
         L) ||  // проверка растояния до первого угла полигона
        (((*poligon_iter)->GetP2() - cell_center).Mod() <
         L) ||  // проверка растояния до второго угла полигона
        (((*poligon_iter)->GetP3() - cell_center).Mod() <
         L))  // проверка растояния до третьего угла полигона
    {
      poligon_ptrs.push_back(
          (*poligon_iter)->GetPtr());  // добавление указателя на полигон
    } else {
      double dist =
          (cell_center - (*poligon_iter)->GetGmt()) * (*poligon_iter)->normal;
      if ((dist < L) && (dist > -L)) {
        Point point_on_poligon = cell_center + (*poligon_iter)->normal * dist;

        (*poligon_iter)->Fix();

        Point a = (*poligon_iter)->GetP2() - (*poligon_iter)->GetP1();
        Point b = (*poligon_iter)->GetP3() - (*poligon_iter)->GetP2();
        Point c = (*poligon_iter)->GetP1() - (*poligon_iter)->GetP3();

        Point d1 = (*poligon_iter)->GetP1() - point_on_poligon;
        Point d2 = (*poligon_iter)->GetP2() - point_on_poligon;
        Point d3 = (*poligon_iter)->GetP3() - point_on_poligon;

        if ((*poligon_iter)->GetNormal() * d1.Cross(a) > 0. &&
            (*poligon_iter)->GetNormal() * d2.Cross(b) > 0. &&
            (*poligon_iter)->GetNormal() * d3.Cross(c) > 0.) {
          poligon_ptrs.push_back(
              (*poligon_iter)->GetPtr());  // добавление указателя на полигон
        }
      }
    }

    poligon_iter++;
  }

  if (poligon_ptrs.size() >= 1)
    return true;
  else
    return false;
}

int InnerBoundary::BoundaryCondition(deque<Particle>* cluster_particle,
                                     double dt) {
  // if(poligon_ptrs.size() <= 0) return 1;

  deque<Particle>::iterator particle_iter = cluster_particle->begin();

  while (particle_iter != cluster_particle->end()) {
    deque<Polygon*>::iterator poligon_iterator = this->poligon_ptrs.begin();
    deque<Polygon*>::iterator collision_poligon_iterator =
        this->poligon_ptrs.end();

    double dtt = dt;
    double particle_dt = dt;

    bool collision_mark = false;

    while (particle_dt > 0) {
      while (poligon_iterator !=
             poligon_ptrs
                 .end())  // ищем полигон с которым будет соударятся
                          // частица(полигон к которому частица прелетит первой)
      {
        if ((*poligon_iterator)->GetNormal() * particle_iter->GetVelocity() <
            0)  // проверка направления скорости частицы на возможность
                // соударения
        {
          Point poligon_gmt = (*poligon_iterator)->GetGmt();

          double A =
              particle_iter->GetVelocity() * (*poligon_iterator)->GetNormal();

          double tc =
              ((*poligon_iterator)->GetP1() - particle_iter->GetPosition()) *
              (*poligon_iterator)->GetNormal() / A;

          if (tc <= 0) {
            poligon_iterator++;
            continue;
          } else if (tc > dtt) {
            poligon_iterator++;
            continue;
          }
          // else if(tc > min_dt) continue;
          else {
            Point collision_pstn = particle_iter->GetPosition() +
                                   particle_iter->GetVelocity() * tc;

            Point a =
                (*poligon_iterator)->GetP2() - (*poligon_iterator)->GetP1();
            Point b =
                (*poligon_iterator)->GetP3() - (*poligon_iterator)->GetP2();
            Point c =
                (*poligon_iterator)->GetP1() - (*poligon_iterator)->GetP3();

            Point d1 = (*poligon_iterator)->GetP1() - collision_pstn;
            Point d2 = (*poligon_iterator)->GetP2() - collision_pstn;
            Point d3 = (*poligon_iterator)->GetP3() - collision_pstn;

            if ((*poligon_iterator)->GetNormal() * d1.Cross(a) > 0. &&
                (*poligon_iterator)->GetNormal() * d2.Cross(b) > 0. &&
                (*poligon_iterator)->GetNormal() * d3.Cross(c) > 0.) {
              dtt = tc;
              collision_poligon_iterator = poligon_iterator;
              collision_mark = true;
            }
          }
        }

        poligon_iterator++;
      }  // while(poligon_iterator != poligon_ptrs.end())

      if (collision_mark) {
        Point normal = (*collision_poligon_iterator)->GetNormal();
        // particle_iter->position += particle_iter->velocity * dtt + normal *
        // 0.000001;
        particle_iter->position += particle_iter->velocity * dtt;

        (*collision_poligon_iterator)->force += particle_iter->velocity;

        const double rmt = 1.0 / static_cast<double>(RAND_MAX);

        double r1 = std::rand() * rmt;
        double r2 = std::rand() * rmt;
        double r3 = std::rand() * rmt;

        if (r1 == 0.) r1 = 0.00001;
        if (r2 == 0.) r2 = 0.00001;

        Point vel = particle_iter->velocity;

        double sp = vel * normal;
        double r = sqrt(2. * Tw * fabs(log(r1)));

        Point vn;
        vn = normal * r;

        Point vni;
        vni = normal * sp;

        Point vt;
        vt = vel - vni;
        double lvt = vt.Mod();
        if (lvt < eps) lvt = eps;
        vt /= lvt;

        r = sqrt(2. * Tw * fabs(log(r2)));
        double teta = 2. * Pi * r3;
        double vt1m = r * cos(teta);
        double vt2m = r * sin(teta);

        Point vt1, vt2;

        vt1 *= vt1m;

        vt2 = normal.Cross(vt);

        vt2 *= vt2m;

        particle_iter->velocity = vt1 + vt2 + vn;

        (*collision_poligon_iterator)->force -= particle_iter->velocity;

        /*
                                                deque<Particle>::iterator
           particle_iter2 = cluster_particle->begin(); while(particle_iter2 !=
           cluster_particle->end())
                                                {
                                                        if(particle_iter2 !=
           particle_iter)
                                                        {
                                                                particle_iter2->position
           += particle_iter2->velocity * dtt;
                                                        }

                                                        particle_iter2++;
                                                }*/

        particle_dt -= dtt;
      } else {
        particle_iter->position += particle_iter->velocity * dtt;
        particle_dt = 0;
      }

      collision_mark = false;
    }  // while(dt > 0)

    particle_iter++;
  }

  return 0;
}

double InnerBoundary::CalcCellVolume(
    Point cell_apex, Point size,
    Point* mass_center_out)  // вычисление отсеченного обёма ячейки, аргументы:
                             // входные параметры:cell_apex - опорный угол
                             // ячейки, size - размеры ячейки, выходные
                             // параметры: mass_center - центр мас ячейки.
{
  double volume = 0;
  double tot_vol = size.x * size.y * size.z;
  double dtx = size.x / 1;
  double dty = size.y / 1;
  double dtz = size.z / 1;

  Point mass_center(0, 0, 0);

  if (poligon_ptrs.size() == 0) {
    // Заглушка!!!
    return tot_vol;
  }

  size_t N = 100000;
  size_t outer_N = N;

  std::vector<Point> rand_points;
  rand_points.reserve(N);

  for (size_t i = 0; i < N; i++) {
    rand_points.emplace_back(cell_apex + size.RandomPoint());
  }

  auto poligon_col_iterator_x = poligon_ptrs.end();
  auto poligon_col_iterator_y = poligon_ptrs.end();
  auto poligon_col_iterator_z = poligon_ptrs.end();

  Point collision_pstn_x;
  Point collision_pstn_y;
  Point collision_pstn_z;

  for (size_t i = 0; i < N; i++) {
    double dttx = size.x / 1;
    double dtty = size.y / 1;
    double dttz = size.z / 1;

    auto poligon_iterator = poligon_ptrs.begin();

    while (poligon_iterator !=
           poligon_ptrs.end())  // ищем полигон с которым будет соударятся
                                // виртуальная частица
    {
      Point poligon_gmt = (*poligon_iterator)->GetGmt();

      Point a = (*poligon_iterator)->GetP2() - (*poligon_iterator)->GetP1();
      Point b = (*poligon_iterator)->GetP3() - (*poligon_iterator)->GetP2();
      Point c = (*poligon_iterator)->GetP1() - (*poligon_iterator)->GetP3();

      double Ax = Point(1, 0, 0) * (*poligon_iterator)->GetNormal();
      double tcx = ((*poligon_iterator)->GetP1() - rand_points[i]) *
                   (*poligon_iterator)->GetNormal() / Ax;

      if (tcx < dttx && tcx > -dttx) {
        Point collision_pstn = rand_points[i] + Point(1, 0, 0) * tcx;

        Point d1 = (*poligon_iterator)->GetP1() - collision_pstn;
        Point d2 = (*poligon_iterator)->GetP2() - collision_pstn;
        Point d3 = (*poligon_iterator)->GetP3() - collision_pstn;

        if ((*poligon_iterator)->GetNormal() * d1.Cross(a) > 0. &&
            (*poligon_iterator)->GetNormal() * d2.Cross(b) > 0. &&
            (*poligon_iterator)->GetNormal() * d3.Cross(c) > 0.) {
          poligon_col_iterator_x = poligon_iterator;
          collision_pstn_x = collision_pstn;
          dttx = tcx;
        }
      }

      double Ay = Point(0, 1, 0) * (*poligon_iterator)->GetNormal();
      double tcy = ((*poligon_iterator)->GetP1() - rand_points[i]) *
                   (*poligon_iterator)->GetNormal() / Ay;

      if (tcy < dtty && tcy > -dtty) {
        Point collision_pstn = rand_points[i] + Point(1, 0, 0) * tcy;

        Point d1 = (*poligon_iterator)->GetP1() - collision_pstn;
        Point d2 = (*poligon_iterator)->GetP2() - collision_pstn;
        Point d3 = (*poligon_iterator)->GetP3() - collision_pstn;

        if ((*poligon_iterator)->GetNormal() * d1.Cross(a) > 0. &&
            (*poligon_iterator)->GetNormal() * d2.Cross(b) > 0. &&
            (*poligon_iterator)->GetNormal() * d3.Cross(c) > 0.) {
          poligon_col_iterator_y = poligon_iterator;
          collision_pstn_y = collision_pstn;
          dtty = tcy;
        }
      }

      double Az = Point(0, 0, 1) * (*poligon_iterator)->GetNormal();
      double tcz = ((*poligon_iterator)->GetP1() - rand_points[i]) *
                   (*poligon_iterator)->GetNormal() / Az;

      if (tcz < dttz && tcz > -dttz) {
        Point collision_pstn = rand_points[i] + Point(1, 0, 0) * tcz;

        Point d1 = (*poligon_iterator)->GetP1() - collision_pstn;
        Point d2 = (*poligon_iterator)->GetP2() - collision_pstn;
        Point d3 = (*poligon_iterator)->GetP3() - collision_pstn;

        if ((*poligon_iterator)->GetNormal() * d1.Cross(a) > 0. &&
            (*poligon_iterator)->GetNormal() * d2.Cross(b) > 0. &&
            (*poligon_iterator)->GetNormal() * d3.Cross(c) > 0.) {
          poligon_col_iterator_z = poligon_iterator;
          collision_pstn_z = collision_pstn;
          dttz = tcz;
        }
      }
      poligon_iterator++;
    }  // while(poligon_iterator != poligon_ptrs.end())

    if (poligon_col_iterator_x != poligon_ptrs.end()) {
      if ((*poligon_col_iterator_x)->GetNormal() *
              (collision_pstn_x - rand_points[i]) <=
          0) {
        outer_N--;
        break;
      }
    }

    if (poligon_col_iterator_y != poligon_ptrs.end()) {
      if ((*poligon_col_iterator_y)->GetNormal() *
              (collision_pstn_y - rand_points[i]) <=
          0) {
        outer_N--;
        break;
      }
    }

    if (poligon_col_iterator_z != poligon_ptrs.end()) {
      if ((*poligon_col_iterator_z)->GetNormal() *
              (collision_pstn_z - rand_points[i]) <=
          0) {
        outer_N--;
        break;
      }
    }

    mass_center += rand_points[i];
  }

  volume = tot_vol * (double(outer_N) / double(N));
  if (mass_center_out != NULL) {
    mass_center /= outer_N;
    *mass_center_out = mass_center;
  }

  if (outer_N == N) {
    return -1;
  }

  return volume;
}

void InnerBoundary::SetGeometry(Geometry* bbody) { this->body = bbody; }

Geometry* InnerBoundary::GetGeometryPtr() { return body; }

bool InnerBoundary::Empty() const { return poligon_ptrs.empty(); }
}  // namespace mc3d
