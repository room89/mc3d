#include "inner_boundary.h"

namespace mc3d {
const double Pi = 3.1415926535;

inner_boundary::inner_boundary(void) {
  Tw = 1;
  eps = 0.00001;
}

inner_boundary::~inner_boundary(void) {}

bool inner_boundary::add_poligon(geometry* body, point cell_center, double L) {
  /*Tw = 1;
  eps = 0.00001;

  this->cell_center = cell_center;
  
  deque<poligon*>::iterator poligon_iter = bbody->poligons.begin();

  while(poligon_iter != bbody->poligons.end())
  {
          if(((*poligon_iter)->get_p1() - cell_center) *
  (*poligon_iter)->get_normal() < L)
  //проверка растояния до плоскости полигона
          {
                  if( (((*poligon_iter)->get_gmt() - cell_center).mod() < 2 * L)
  ||
  //проверка растояния до ГМТ полигона(полигон - треугольник)
                          (((*poligon_iter)->get_p1()  - cell_center).mod() < 2
  * L) ||
  //проверка растояния до первого угла полигона
                          (((*poligon_iter)->get_p2()  - cell_center).mod() < 2
  * L) ||
  //проверка растояния до второго угла полигона
                          (((*poligon_iter)->get_p3()  - cell_center).mod() < 2
  * L))
  //проверка растояния до третьего угла полигона
                  {
                                  poligon_ptrs.push_back((*poligon_iter)->get_ptr());
  //добавление указателя на полигон
                  }
          }

          poligon_iter++;
  }

  if(poligon_ptrs.size() > 0) return true;
  else return false;*/
  this->cell_center = cell_center;

  deque<poligon*>::iterator poligon_iter = body->poligons.begin();

  while (poligon_iter != body->poligons.end()) {
    if ((((*poligon_iter)->get_gmt() - cell_center).mod() <
         L) ||  //проверка растояния до ГМТ полигона(полигон - треугольник)
        (((*poligon_iter)->get_p1() - cell_center).mod() <
         L) ||  //проверка растояния до первого угла полигона
        (((*poligon_iter)->get_p2() - cell_center).mod() <
         L) ||  //проверка растояния до второго угла полигона
        (((*poligon_iter)->get_p3() - cell_center).mod() <
         L))  //проверка растояния до третьего угла полигона
    {
      poligon_ptrs.push_back(
          (*poligon_iter)->get_ptr());  //добавление указателя на полигон
    } else {
      double dist =
          (cell_center - (*poligon_iter)->get_gmt()) * (*poligon_iter)->normal;
      if ((dist < L) && (dist > -L)) {
        point point_on_poligon = cell_center + (*poligon_iter)->normal * dist;

        (*poligon_iter)->fix();

        point a = (*poligon_iter)->get_p2() - (*poligon_iter)->get_p1();
        point b = (*poligon_iter)->get_p3() - (*poligon_iter)->get_p2();
        point c = (*poligon_iter)->get_p1() - (*poligon_iter)->get_p3();

        point d1 = (*poligon_iter)->get_p1() - point_on_poligon;
        point d2 = (*poligon_iter)->get_p2() - point_on_poligon;
        point d3 = (*poligon_iter)->get_p3() - point_on_poligon;

        if ((*poligon_iter)->get_normal() * d1.vec_mult(a) > 0. &&
            (*poligon_iter)->get_normal() * d2.vec_mult(b) > 0. &&
            (*poligon_iter)->get_normal() * d3.vec_mult(c) > 0.) {
          poligon_ptrs.push_back(
              (*poligon_iter)->get_ptr());  //добавление указателя на полигон
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

int inner_boundary::bondary_condition(deque<particle>* cluster_particle,
                                      double dt) {
  // if(poligon_ptrs.size() <= 0) return 1;

  deque<particle>::iterator particle_iter = cluster_particle->begin();

  while (particle_iter != cluster_particle->end()) {
    deque<poligon*>::iterator poligon_iterator = this->poligon_ptrs.begin();
    deque<poligon*>::iterator collision_poligon_iterator =
        this->poligon_ptrs.end();

    double dtt = dt;
    double particle_dt = dt;

    bool collision_mark = false;

    while (particle_dt > 0) {
      while (poligon_iterator !=
             poligon_ptrs
                 .end())  //ищем полигон с которым будет соударятся
                          //частица(полигон к которому частица прелетит первой)
      {
        if ((*poligon_iterator)->get_normal() * particle_iter->get_velocity() <
            0)  //проверка направления скорости частицы на возможность
                //соударения
        {
          point poligon_gmt = (*poligon_iterator)->get_gmt();

          double A =
              particle_iter->get_velocity() * (*poligon_iterator)->get_normal();

          double tc =
              ((*poligon_iterator)->get_p1() - particle_iter->get_position()) *
              (*poligon_iterator)->get_normal() / A;

          if (tc <= 0) {
            poligon_iterator++;
            continue;
          } else if (tc > dtt) {
            poligon_iterator++;
            continue;
          }
          // else if(tc > min_dt) continue;
          else {
            point collision_pstn = particle_iter->get_position() +
                                   particle_iter->get_velocity() * tc;

            point a =
                (*poligon_iterator)->get_p2() - (*poligon_iterator)->get_p1();
            point b =
                (*poligon_iterator)->get_p3() - (*poligon_iterator)->get_p2();
            point c =
                (*poligon_iterator)->get_p1() - (*poligon_iterator)->get_p3();

            point d1 = (*poligon_iterator)->get_p1() - collision_pstn;
            point d2 = (*poligon_iterator)->get_p2() - collision_pstn;
            point d3 = (*poligon_iterator)->get_p3() - collision_pstn;

            if ((*poligon_iterator)->get_normal() * d1.vec_mult(a) > 0. &&
                (*poligon_iterator)->get_normal() * d2.vec_mult(b) > 0. &&
                (*poligon_iterator)->get_normal() * d3.vec_mult(c) > 0.) {
              dtt = tc;
              collision_poligon_iterator = poligon_iterator;
              collision_mark = true;
            }
          }
        }

        poligon_iterator++;
      }  // while(poligon_iterator != poligon_ptrs.end())

      if (collision_mark) {
        point normal = (*collision_poligon_iterator)->get_normal();
        // particle_iter->position += particle_iter->velocity * dtt + normal *
        // 0.000001;
        particle_iter->position += particle_iter->velocity * dtt;

        (*collision_poligon_iterator)->force += particle_iter->velocity;

        double rmt = 1 / double(RAND_MAX);

        double r1 = std::rand() * rmt;
        double r2 = std::rand() * rmt;
        double r3 = std::rand() * rmt;

        if (r1 == 0.) r1 = 0.00001;
        if (r2 == 0.) r2 = 0.00001;

        point vel = particle_iter->velocity;

        double sp = vel * normal;
        double r = sqrt(2. * Tw * fabs(log(r1)));

        point vn;
        vn = normal * r;

        point vni;
        vni = normal * sp;

        point vt;
        vt = vel - vni;
        double lvt = vt.mod();
        if (lvt < eps) lvt = eps;
        vt /= lvt;

        r = sqrt(2. * Tw * fabs(log(r2)));
        double teta = 2. * Pi * r3;
        double vt1m = r * cos(teta);
        double vt2m = r * sin(teta);

        point vt1, vt2;

        vt1 *= vt1m;

        vt2 = normal.vec_mult(vt);

        vt2 *= vt2m;

        particle_iter->velocity = vt1 + vt2 + vn;

        (*collision_poligon_iterator)->force -= particle_iter->velocity;

        /*
                                                deque<particle>::iterator
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

double inner_boundary::calc_cell_volume(
    point cell_apex, point size,
    point* mass_center_out)  //вычисление отсеченного обёма ячейки, аргументы:
                             //входные параметры:cell_apex - опорный угол
                             //ячейки, size - размеры ячейки, выходные
                             //параметры: mass_center - центр мас ячейки.
{
  double volume = 0;
  double tot_vol = size.x * size.y * size.z;
  double dtx = size.x / 1;
  double dty = size.y / 1;
  double dtz = size.z / 1;

  point mass_center(0, 0, 0);

  if (poligon_ptrs.size() == 0) {
    //Заглушка!!!
    return tot_vol;
  }

  int N = 100000;
  int outer_N = N;

  point* rand_point = new point[N];

  // int *outer_point_index = new int[N];

  for (int i = 0; i < N; i++) {
    rand_point[i] = cell_apex + size.rand_point();
  }

  auto poligon_col_iterator_x = poligon_ptrs.end();
  auto poligon_col_iterator_y = poligon_ptrs.end();
  auto poligon_col_iterator_z = poligon_ptrs.end();

  point collision_pstn_x;
  point collision_pstn_y;
  point collision_pstn_z;

  for (int i = 0; i < N; i++) {
    double dttx = size.x / 1;
    double dtty = size.y / 1;
    double dttz = size.z / 1;

    auto poligon_iterator = poligon_ptrs.begin();

    while (poligon_iterator !=
           poligon_ptrs.end())  //ищем полигон с которым будет соударятся
                                //виртуальная частица
    {
      point poligon_gmt = (*poligon_iterator)->get_gmt();

      point a = (*poligon_iterator)->get_p2() - (*poligon_iterator)->get_p1();
      point b = (*poligon_iterator)->get_p3() - (*poligon_iterator)->get_p2();
      point c = (*poligon_iterator)->get_p1() - (*poligon_iterator)->get_p3();

      double Ax = point(1, 0, 0) * (*poligon_iterator)->get_normal();
      double tcx = ((*poligon_iterator)->get_p1() - rand_point[i]) *
                   (*poligon_iterator)->get_normal() / Ax;

      if (tcx < dttx && tcx > -dttx) {
        point collision_pstn = rand_point[i] + point(1, 0, 0) * tcx;

        point d1 = (*poligon_iterator)->get_p1() - collision_pstn;
        point d2 = (*poligon_iterator)->get_p2() - collision_pstn;
        point d3 = (*poligon_iterator)->get_p3() - collision_pstn;

        if ((*poligon_iterator)->get_normal() * d1.vec_mult(a) > 0. &&
            (*poligon_iterator)->get_normal() * d2.vec_mult(b) > 0. &&
            (*poligon_iterator)->get_normal() * d3.vec_mult(c) > 0.) {
          poligon_col_iterator_x = poligon_iterator;
          collision_pstn_x = collision_pstn;
          dttx = tcx;
        }
      }

      double Ay = point(0, 1, 0) * (*poligon_iterator)->get_normal();
      double tcy = ((*poligon_iterator)->get_p1() - rand_point[i]) *
                   (*poligon_iterator)->get_normal() / Ay;

      if (tcy < dtty && tcy > -dtty) {
        point collision_pstn = rand_point[i] + point(1, 0, 0) * tcy;

        point d1 = (*poligon_iterator)->get_p1() - collision_pstn;
        point d2 = (*poligon_iterator)->get_p2() - collision_pstn;
        point d3 = (*poligon_iterator)->get_p3() - collision_pstn;

        if ((*poligon_iterator)->get_normal() * d1.vec_mult(a) > 0. &&
            (*poligon_iterator)->get_normal() * d2.vec_mult(b) > 0. &&
            (*poligon_iterator)->get_normal() * d3.vec_mult(c) > 0.) {
          poligon_col_iterator_y = poligon_iterator;
          collision_pstn_y = collision_pstn;
          dtty = tcy;
        }
      }

      double Az = point(0, 0, 1) * (*poligon_iterator)->get_normal();
      double tcz = ((*poligon_iterator)->get_p1() - rand_point[i]) *
                   (*poligon_iterator)->get_normal() / Az;

      if (tcz < dttz && tcz > -dttz) {
        point collision_pstn = rand_point[i] + point(1, 0, 0) * tcz;

        point d1 = (*poligon_iterator)->get_p1() - collision_pstn;
        point d2 = (*poligon_iterator)->get_p2() - collision_pstn;
        point d3 = (*poligon_iterator)->get_p3() - collision_pstn;

        if ((*poligon_iterator)->get_normal() * d1.vec_mult(a) > 0. &&
            (*poligon_iterator)->get_normal() * d2.vec_mult(b) > 0. &&
            (*poligon_iterator)->get_normal() * d3.vec_mult(c) > 0.) {
          poligon_col_iterator_z = poligon_iterator;
          collision_pstn_z = collision_pstn;
          dttz = tcz;
        }
      }
      poligon_iterator++;
    }  // while(poligon_iterator != poligon_ptrs.end())

    if (poligon_col_iterator_x != poligon_ptrs.end()) {
      if ((*poligon_col_iterator_x)->get_normal() *
              (collision_pstn_x - rand_point[i]) <=
          0) {
        outer_N--;
        break;
      }
    }

    if (poligon_col_iterator_y != poligon_ptrs.end()) {
      if ((*poligon_col_iterator_y)->get_normal() *
              (collision_pstn_y - rand_point[i]) <=
          0) {
        outer_N--;
        break;
      }
    }

    if (poligon_col_iterator_z != poligon_ptrs.end()) {
      if ((*poligon_col_iterator_z)->get_normal() *
              (collision_pstn_z - rand_point[i]) <=
          0) {
        outer_N--;
        break;
      }
    }

    mass_center += rand_point[i];
  }

  volume = tot_vol * (double(outer_N) / double(N));
  if (mass_center_out != NULL) {
    mass_center /= outer_N;
    *mass_center_out = mass_center;
  }

  delete[] rand_point;

  if (outer_N == N) {
    return -1;
  }

  return volume;
}

void inner_boundary::set_geometry(geometry* bbody) { this->body = bbody; }

geometry* inner_boundary::get_geometry_ptr() { return body; }
}  // namespace mc3d