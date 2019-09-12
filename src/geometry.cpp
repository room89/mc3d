#include "geometry.h"

namespace mc3d {
geometry::geometry(void) {}

geometry::geometry(const char* file_name) {
  std::cout << endl
            << "...reading geometry from file " << file_name << "... " << endl;
  // unsigned int n = 0;
  string word;
  vector<string> v;
  std::ifstream file(file_name);
  if (!file.is_open()) {
    std::cout << endl << "Can't open file!!!" << endl;
    throw unusual_situations::exception(
        mc3d::unusual_situations::exit_code::BOUNDARY_CONSTRUCTOR_ERROR, this);
  }
  while (file >> word) {
    v.push_back(word);
  }
  if (v[0] == "solid" && v[v.size() - 1] == "endsolid") {
    unsigned int i = 1;
    unsigned int n = 0;
    if (v[1] != "facet") {
      std::cout << endl << "Solid: " << v[1] << endl;
      i++;
    }
    while (i < v.size()) {
      if (v[i] == "facet" && v[i + 1] == "normal" && v[i + 6] == "loop" &&
          v[i + 19] == "endloop" && v[i + 20] == "endfacet") {
        if (v[i] == "endsolid") break;
        mc3d::poligon* temp = new poligon;
        mc3d::point p1, p2, p3, normal;

        p1.x = std::atof(v[i + 8].c_str());
        p1.y = std::atof(v[i + 9].c_str());
        p1.z = std::atof(v[i + 10].c_str());

        p2.x = std::atof(v[i + 12].c_str());
        p2.y = std::atof(v[i + 13].c_str());
        p2.z = std::atof(v[i + 14].c_str());

        p3.x = std::atof(v[i + 16].c_str());
        p3.y = std::atof(v[i + 17].c_str());
        p3.z = std::atof(v[i + 18].c_str());

        normal.x = std::atof(v[i + 2].c_str());
        normal.y = std::atof(v[i + 3].c_str());
        normal.z = std::atof(v[i + 4].c_str());

        normal.nrmlz();

        temp->set(p1, p2, p3, normal);

        poligons.push_back(temp);

        i += 21;
      } else {
        i++;
      }
    }
  }

  file.close();
}

geometry::~geometry(void) {
  deque<poligon*>::iterator iter = poligons.begin();
  while (iter != poligons.end()) {
    delete (*iter);
    iter++;
  }
  poligons.clear();
}

//заготовка
// int geometry::bondary_condition(deque<particle> *cluster_particle, double dt)
//{
//	deque<particle>::iterator particle_iter = cluster_particle->begin();
//
//	while(particle_iter != cluster_particle->end())
//	{
//		deque<poligon>::iterator poligon_iterator = poligons.begin();
//
//		double dtt = dt;
//
//		while(dtt > 0)
//		{
//			while(poligon_iterator != poligons.end())
//			{
//				if((*poligon_iterator)->get_normal() *
// particle_iter->get_velocity() < 0)	//проверка направления скорости частицы
//на возможность соударения
//				{
//					point poligon_gmt =
//(*poligon_iterator)->get_gmt();
//
//					double A = particle_iter->get_velocity()
//*
//(*poligon_iterator)->get_normal();
//
//					double tc =
//((*poligon_iterator)->get_p1()
//- particle_iter->get_position()) * (*poligon_iterator)->get_normal() / A;
//
//					if(tc <= 0) continue;
//					if(tc > dt) continue;
//				}
//
//				poligon_iterator++;
//			}
//		}
//
//		particle_iter++;
//	}
//
//	return 0;
//}

void geometry::write_geometry_file(const char* file_name,
                                   const char* body_name) {
  std::ofstream file(file_name);

  if (!file.is_open()) return;

  file << "solid " << body_name << std::endl;

  std::deque<poligon*>::iterator i = poligons.begin();

  while (i != poligons.end()) {
    file << "facet normal " << std::setw(15) << (*i)->get_normal().x
         << std::setw(15) << (*i)->get_normal().y << std::setw(15)
         << (*i)->get_normal().z << std::endl;
    file << "\touter loop" << std::endl;
    file << "\t\tvertex" << std::setw(15) << (*i)->get_p1().x << std::setw(15)
         << (*i)->get_p1().y << std::setw(15) << (*i)->get_p1().z << std::endl;
    file << "\t\tvertex" << std::setw(15) << (*i)->get_p2().x << std::setw(15)
         << (*i)->get_p2().y << std::setw(15) << (*i)->get_p2().z << std::endl;
    file << "\t\tvertex" << std::setw(15) << (*i)->get_p3().x << std::setw(15)
         << (*i)->get_p3().y << std::setw(15) << (*i)->get_p3().z << std::endl;
    file << "\tendloop" << std::endl << "endfacet" << std::endl;

    i++;
  }

  file << "endsolid" << body_name << std::endl;

  file.close();
}

point geometry::mass_center() {
  point mc(0, 0, 0);

  double Stot = 0;

  std::deque<poligon*>::iterator i = poligons.begin();

  while (i != poligons.end()) {
    point a = (*i)->get_p2() - (*i)->get_p1();
    point b = (*i)->get_p3() - (*i)->get_p2();

    double S = 0.5 * a.vec_mult(b).mod();

    Stot += S;
    mc += S * ((*i)->get_p1() + (*i)->get_p2() + (*i)->get_p3()) / 3;

    i++;
  }

  if (Stot > 0) {
    mc /= Stot;
  }

  return mc;
}

void geometry::move(point a) {
  std::deque<poligon*>::iterator i = poligons.begin();

  while (i != poligons.end()) {
    (*i)->move(a);
    i++;
  }
}

void geometry::scaling(double e) {
  std::deque<poligon*>::iterator i = poligons.begin();

  while (i != poligons.end()) {
    (*i)->scaling(e);

    i++;
  }
}

bool geometry::is_inner_point(point test_point) const {
  size_t n = 0;
  bool res = false;

  auto poligon_iter = poligons.begin();

  while (poligon_iter != poligons.end()) {
    /*if(dist == (*poligon_iter)->dist_to_point(test_point))
    {
            res = res;
    }*/
    /*if(dist > (*poligon_iter)->dist_to_point(test_point))
    {
            dist = (*poligon_iter)->dist_to_point(test_point);
            poligon_col = poligon_iter;
    }*/

    point gmt = (*poligon_iter)->get_gmt();
    point asd = test_point - (*poligon_iter)->get_gmt();
    if ((*poligon_iter)->get_normal() *
            (test_point - (*poligon_iter)->get_gmt()) <
        0) {
      n++;
    }

    poligon_iter++;
  }

  if (n == poligons.size()) {
    return true;
  }

  return res;
}

/*bool geometry::is_inner_point(point test_point)
{
        deque<poligon*>::iterator poligon_col_iterator_x = poligons.end();
        deque<poligon*>::iterator poligon_col_iterator_y = poligons.end();
        deque<poligon*>::iterator poligon_col_iterator_z = poligons.end();

        point collision_pstn_x(0, 0, 0);
        point collision_pstn_y(0, 0, 0);
        point collision_pstn_z(0, 0, 0);

        deque<poligon*>::iterator poligon_iterator = poligons.begin();

        double dttx = 100000;
        double dtty = 100000;
        double dttz = 100000;

        while(poligon_iterator != poligons.end())
//ищем полигон с которым будет соударятся виртуальная частица
        {
                point poligon_gmt = (*poligon_iterator)->get_gmt();

                point a = (*poligon_iterator)->get_p2() -
(*poligon_iterator)->get_p1(); point b = (*poligon_iterator)->get_p3() -
(*poligon_iterator)->get_p2(); point c = (*poligon_iterator)->get_p1() -
(*poligon_iterator)->get_p3();

                double Ax = point(1, 0, 0) * (*poligon_iterator)->get_normal();
                double tcx = 0;
                if(Ax != 0) tcx	 = ((*poligon_iterator)->get_p1() - test_point)
* (*poligon_iterator)->get_normal() / Ax;

                if(tcx < dttx && tcx > -dttx)
                {
                        point collision_pstn = test_point + point(1, 0, 0) *
tcx;
                         
                        point d1 = (*poligon_iterator)->get_p1() -
collision_pstn; point d2 = (*poligon_iterator)->get_p2() - collision_pstn; point
d3 = (*poligon_iterator)->get_p3() - collision_pstn;

                        if( (*poligon_iterator)->get_normal() * d1.vec_mult(a) >
0. &&
                                (*poligon_iterator)->get_normal() *
d2.vec_mult(b) > 0. &&
                                (*poligon_iterator)->get_normal() *
d3.vec_mult(c) > 0.)
                        {
                                poligon_col_iterator_x = poligon_iterator;
                                collision_pstn_x = 	collision_pstn;
                                dttx = tcx;
                        }
                }
                        
                        
                double Ay = point(0, 1, 0) * (*poligon_iterator)->get_normal();
                double tcy = ((*poligon_iterator)->get_p1() - test_point) *
(*poligon_iterator)->get_normal() / Ay;

                if(tcy < dtty && tcy > -dtty)
                {
                        point collision_pstn = test_point + point(1, 0, 0) *
tcy;
                         
                        point d1 = (*poligon_iterator)->get_p1() -
collision_pstn; point d2 = (*poligon_iterator)->get_p2() - collision_pstn; point
d3 = (*poligon_iterator)->get_p3() - collision_pstn;

                        if( (*poligon_iterator)->get_normal() * d1.vec_mult(a) >
0. &&
                                (*poligon_iterator)->get_normal() *
d2.vec_mult(b) > 0. &&
                                (*poligon_iterator)->get_normal() *
d3.vec_mult(c) > 0.)
                        {
                                poligon_col_iterator_y = poligon_iterator;
                                collision_pstn_y = 	collision_pstn;
                                dtty = tcy;
                        }
                }

                double Az = point(0, 0, 1) * (*poligon_iterator)->get_normal();
                double tcz = ((*poligon_iterator)->get_p1() - test_point) *
(*poligon_iterator)->get_normal() / Az;

                if(tcz < dttz && tcz > -dttz)
                {
                        point collision_pstn = test_point + point(1, 0, 0) *
tcz;
                         
                        point d1 = (*poligon_iterator)->get_p1() -
collision_pstn; point d2 = (*poligon_iterator)->get_p2() - collision_pstn; point
d3 = (*poligon_iterator)->get_p3() - collision_pstn;

                        if( (*poligon_iterator)->get_normal() * d1.vec_mult(a) >
0. &&
                                (*poligon_iterator)->get_normal() *
d2.vec_mult(b) > 0. &&
                                (*poligon_iterator)->get_normal() *
d3.vec_mult(c) > 0.)
                        {
                                poligon_col_iterator_z = poligon_iterator;
                                collision_pstn_z = 	collision_pstn;
                                dttz = tcz;
                        }
                }
                poligon_iterator++;
        }//while(poligon_iterator != poligon_ptrs.end())

        if(poligon_col_iterator_x != poligons.end())
        {
                if((*poligon_col_iterator_x)->get_normal() * (collision_pstn_x -
test_point) < 0)
                {
                        return true;
                }
        }

        if(poligon_col_iterator_y != poligons.end())
        {
                if((*poligon_col_iterator_y)->get_normal() * (collision_pstn_y -
test_point) < 0)
                {
                        return true;
                }
        }

        if(poligon_col_iterator_z != poligons.end())
        {
                if((*poligon_col_iterator_z)->get_normal() * (collision_pstn_z -
test_point) < 0)
                {
                        return true;
                }
        }
        if((poligon_col_iterator_x == poligons.end()) && (poligon_col_iterator_y
== poligons.end()) && (poligon_col_iterator_z == poligons.end()))
        {
                return true;
        }
        return false;
}*/

/*void geometry::create_cone(point apex1, point apex2, double H)
{
        poligons.clear();
        point vertices[6];

        vertices[0] = apex1;
        vertices[1] = apex2;
        vertices[2] = point(apex1.x, apex1.y, apex2.z);
        vertices[3] = point(apex2.x, apex2.y, apex1.z);
        vertices[4] = apex2 - point(0, 0, H);
        vertices[5] = vertices[3] - point(0, 0, H);

        poligon *temp_poligon = new poligon;
        temp_poligon->set(vertices[0], vertices[1], vertices[2], (vertices[1] -
vertices[0]).vec_mult(vertices[2] - vertices[1]).nrmlz());
        poligons.push_back(temp_poligon);

        temp_poligon = new poligon;
        temp_poligon->set(vertices[0], vertices[1], vertices[3], (vertices[3] -
vertices[1]).vec_mult(vertices[1] - vertices[0]).nrmlz());
        poligons.push_back(temp_poligon);

        temp_poligon = new poligon;
        temp_poligon->set(vertices[2], vertices[1], vertices[4], (vertices[1] -
vertices[2]).vec_mult(vertices[4] - vertices[2]).nrmlz());
        poligons.push_back(temp_poligon);
        
        temp_poligon = new poligon;
        temp_poligon->set(vertices[0], vertices[3], vertices[5], (vertices[5] -
vertices[0]).vec_mult(vertices[3] - vertices[0]).nrmlz());
        poligons.push_back(temp_poligon);

        temp_poligon = new poligon;
        temp_poligon->set(vertices[0], vertices[2], vertices[4], (vertices[0] -
vertices[2]).vec_mult(vertices[4] - vertices[2]).nrmlz());
        poligons.push_back(temp_poligon);

        temp_poligon = new poligon;
        temp_poligon->set(vertices[0], vertices[4], vertices[5], (vertices[5] -
vertices[0]).vec_mult(vertices[4] - vertices[0]).nrmlz());
        poligons.push_back(temp_poligon);
}*/
void geometry::create_wedge(double x, double width, double length,
                            double alpha) {
  point vertices[6];

  vertices[0].set(x, 0, -(width / 2));
  vertices[1].set(x + length, length * sin(alpha), -width / 2);
  vertices[2].set(x + length, -length * sin(alpha), -width / 2);

  vertices[3].set(x, 0, width / 2);
  vertices[4].set(x + length, length * sin(alpha), width / 2);
  vertices[5].set(x + length, -length * sin(alpha), width / 2);

  poligon* temp_poligon = new poligon;
  temp_poligon->set(vertices[0], vertices[1], vertices[2],
                    -1 * (vertices[0] - vertices[1])
                             .vec_mult(vertices[2] - vertices[1])
                             .nrmlz());
  temp_poligon->flux = 0;
  temp_poligon->force = point(0, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->set(vertices[3], vertices[4], vertices[5],
                    -1 * (vertices[4] - vertices[3])
                             .vec_mult(vertices[5] - vertices[4])
                             .nrmlz());
  temp_poligon->flux = 0;
  temp_poligon->force = point(0, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->set(vertices[0], vertices[1], vertices[3],
                    -1 * (vertices[1] - vertices[0])
                             .vec_mult(vertices[3] - vertices[0])
                             .nrmlz());
  temp_poligon->flux = 0;
  temp_poligon->force = point(0, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->set(vertices[1], vertices[4], vertices[3],
                    -1 * (vertices[4] - vertices[1])
                             .vec_mult(vertices[3] - vertices[1])
                             .nrmlz());
  temp_poligon->flux = 0;
  temp_poligon->force = point(0, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->set(vertices[0], vertices[2], vertices[3],
                    -1 * (vertices[3] - vertices[0])
                             .vec_mult(vertices[2] - vertices[0])
                             .nrmlz());
  temp_poligon->flux = 0;
  temp_poligon->force = point(0, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->set(vertices[2], vertices[3], vertices[5],
                    -1 * (vertices[3] - vertices[2])
                             .vec_mult(vertices[5] - vertices[2])
                             .nrmlz());
  temp_poligon->flux = 0;
  temp_poligon->force = point(0, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->set(vertices[1], vertices[2], vertices[5],
                    -1 * (vertices[2] - vertices[1])
                             .vec_mult(vertices[4] - vertices[1])
                             .nrmlz());
  temp_poligon->flux = 0;
  temp_poligon->force = point(0, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->set(vertices[1], vertices[5], vertices[4],
                    -1 * (vertices[5] - vertices[1])
                             .vec_mult(vertices[4] - vertices[1])
                             .nrmlz());
  temp_poligon->flux = 0;
  temp_poligon->force = point(0, 0, 0);
  poligons.push_back(temp_poligon);
}

void geometry::fragmentation(double Lmax) {
  double max_Lmax = 10000.;
  while (max_Lmax > Lmax) {
    max_Lmax = 0;
    std::deque<poligon*>::iterator poligon_iter = poligons.begin();

    while (poligon_iter != poligons.end()) {
      if ((*poligon_iter)->get_Lmax() > Lmax) {
        max_Lmax > (*poligon_iter)->get_Lmax()
            ? max_Lmax = max_Lmax
            : max_Lmax = (*poligon_iter)->get_Lmax();
        std::pair<poligon*, poligon*> new_poligon = (*poligon_iter)->divide();
        delete *poligon_iter;
        poligon_iter = poligons.erase(poligon_iter);
        poligons.push_back(new_poligon.first);
        poligons.push_back(new_poligon.second);
        poligon_iter = poligons.begin();
        continue;
      }

      poligon_iter++;
    }
  }
}

std::pair<point, point> geometry::size() {
  point p1(0, 0, 0), p2(0, 0, 0);

  std::deque<poligon*>::iterator poligon_iter = poligons.begin();

  while (poligon_iter != poligons.end()) {
    p1.x = std::max(p1.x, (*poligon_iter)->p1.x);
    p1.y = std::max(p1.y, (*poligon_iter)->p1.y);
    p1.z = std::max(p1.z, (*poligon_iter)->p1.z);

    p1.x = std::max(p1.x, (*poligon_iter)->p2.x);
    p1.y = std::max(p1.y, (*poligon_iter)->p2.y);
    p1.z = std::max(p1.z, (*poligon_iter)->p2.z);

    p1.x = std::max(p1.x, (*poligon_iter)->p3.x);
    p1.y = std::max(p1.y, (*poligon_iter)->p3.y);
    p1.z = std::max(p1.z, (*poligon_iter)->p3.z);

    p2.x = std::min(p2.x, (*poligon_iter)->p1.x);
    p2.y = std::min(p2.y, (*poligon_iter)->p1.y);
    p2.z = std::min(p2.z, (*poligon_iter)->p1.z);

    p2.x = std::min(p2.x, (*poligon_iter)->p2.x);
    p2.y = std::min(p2.y, (*poligon_iter)->p2.y);
    p2.z = std::min(p2.z, (*poligon_iter)->p2.z);

    p2.x = std::min(p2.x, (*poligon_iter)->p3.x);
    p2.y = std::min(p2.y, (*poligon_iter)->p3.y);
    p2.z = std::min(p2.z, (*poligon_iter)->p3.z);

    poligon_iter++;
  }

  return pair<point, point>(p2, p1);
}

void geometry::create_pyramid(double x, double width, double length, double H) {
}

double geometry::dist_to_point(point p) {
  double dist = 100000.;

  std::deque<poligon*>::iterator poligon_iter = poligons.begin();

  while (poligon_iter != poligons.end()) {
    dist = std::min(dist, (*poligon_iter)->dist_to_point(p));
    poligon_iter++;
  }

  return dist;
}

void geometry::create_cube(double x, double width, double length, double H) {
  point p[8];
  for (int i = 0; i < 4; i++) {
    p[i].x = x;
    p[i + 4].x = x + length;
  }

  p[0].y = p[3].y = p[4].y = p[7].y = -H / 2;
  p[1].y = p[2].y = p[5].y = p[6].y = H / 2;

  p[0].z = p[1].z = p[4].z = p[5].z = width / 2;
  p[2].z = p[3].z = p[6].z = p[7].z = -width / 2;

  poligon* temp_poligon;

  temp_poligon = new poligon;
  temp_poligon->p1 = p[0];
  temp_poligon->p2 = p[1];
  temp_poligon->p3 = p[2];
  temp_poligon->normal = point(-1, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->p1 = p[0];
  temp_poligon->p2 = p[2];
  temp_poligon->p3 = p[3];
  temp_poligon->normal = point(-1, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->p1 = p[4];
  temp_poligon->p2 = p[5];
  temp_poligon->p3 = p[6];
  temp_poligon->normal = point(1, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->p1 = p[4];
  temp_poligon->p2 = p[6];
  temp_poligon->p3 = p[7];
  temp_poligon->normal = point(1, 0, 0);
  poligons.push_back(temp_poligon);

  /////////////////////////////////////////

  temp_poligon = new poligon;
  temp_poligon->p1 = p[1];
  temp_poligon->p2 = p[2];
  temp_poligon->p3 = p[6];
  temp_poligon->normal = point(0, 1, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->p1 = p[1];
  temp_poligon->p2 = p[6];
  temp_poligon->p3 = p[5];
  temp_poligon->normal = point(0, 1, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->p1 = p[0];
  temp_poligon->p2 = p[3];
  temp_poligon->p3 = p[7];
  temp_poligon->normal = point(0, -1, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->p1 = p[0];
  temp_poligon->p2 = p[7];
  temp_poligon->p3 = p[4];
  temp_poligon->normal = point(0, -1, 0);
  poligons.push_back(temp_poligon);

  /////////////////////////////////////////

  temp_poligon = new poligon;
  temp_poligon->p1 = p[0];
  temp_poligon->p2 = p[1];
  temp_poligon->p3 = p[5];
  temp_poligon->normal = point(0, 0, -1);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->p1 = p[0];
  temp_poligon->p2 = p[5];
  temp_poligon->p3 = p[4];
  temp_poligon->normal = point(0, 0, -1);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->p1 = p[3];
  temp_poligon->p2 = p[2];
  temp_poligon->p3 = p[6];
  temp_poligon->normal = point(0, 0, 1);
  poligons.push_back(temp_poligon);

  temp_poligon = new poligon;
  temp_poligon->p1 = p[3];
  temp_poligon->p2 = p[6];
  temp_poligon->p3 = p[7];
  temp_poligon->normal = point(0, 0, 1);
  poligons.push_back(temp_poligon);
}

int geometry::fix_poligons() {
  int res = 0;
  deque<poligon*>::iterator poligon_iter = poligons.begin();

  while (poligon_iter != poligons.end()) {
    if ((*poligon_iter)->fix()) res++;
    poligon_iter++;
  }

  return res;
}

void geometry::rev_nrml() {
  deque<poligon*>::iterator poligon_iter = poligons.begin();

  while (poligon_iter != poligons.end()) {
    (*poligon_iter)->normal = (*poligon_iter)->normal * -1.;
    poligon_iter++;
  }
}
}  // namespace mc3d
