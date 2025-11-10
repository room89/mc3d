#include "geometry.h"

namespace mc3d {
Geometry::Geometry() {}

Geometry::Geometry(const char* file_name) {
  std::cout << endl
            << "...reading Geometry from file " << file_name << "... " << endl;
  // unsigned int n = 0;
  string word;
  vector<string> v;
  std::ifstream file(file_name);
  if (!file.is_open()) {
    std::cout << endl << "Can't open file!!!" << endl;
    throw unusual_situations::Exception(
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
        mc3d::Polygon* temp = new Polygon;
        mc3d::Point p1, p2, p3, normal;

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

        normal.Normalize();

        temp->Set(p1, p2, p3, normal);

        poligons.push_back(temp);

        i += 21;
      } else {
        i++;
      }
    }
  }

  file.close();
}

Geometry::~Geometry() {
  deque<Polygon*>::iterator iter = poligons.begin();
  while (iter != poligons.end()) {
    delete (*iter);
    iter++;
  }
  poligons.clear();
}

// заготовка
//  int Geometry::bondary_condition(deque<Particle> *cluster_particle, double
//  dt)
//{
//	deque<Particle>::iterator particle_iter = cluster_particle->begin();
//
//	while(particle_iter != cluster_particle->end())
//	{
//		deque<Polygon>::iterator poligon_iterator = poligons.begin();
//
//		double dtt = dt;
//
//		while(dtt > 0)
//		{
//			while(poligon_iterator != poligons.end())
//			{
//				if((*poligon_iterator)->GetNormal() *
//  particle_iter->GetVelocity() < 0)	//проверка направления скорости частицы
// на возможность соударения
//				{
//					Point poligon_gmt =
//(*poligon_iterator)->GetGmt();
//
//					double A = particle_iter->GetVelocity()
//*
//(*poligon_iterator)->GetNormal();
//
//					double tc =
//((*poligon_iterator)->GetP1()
//- particle_iter->GetPosition()) * (*poligon_iterator)->GetNormal() / A;
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
// }

void Geometry::WriteGeometryFile(const char* file_name, const char* body_name) {
  std::ofstream file(file_name);

  if (!file.is_open()) return;

  file << "solid " << body_name << std::endl;

  std::deque<Polygon*>::iterator i = poligons.begin();

  while (i != poligons.end()) {
    file << "facet normal " << std::setw(15) << (*i)->GetNormal().x
         << std::setw(15) << (*i)->GetNormal().y << std::setw(15)
         << (*i)->GetNormal().z << std::endl;
    file << "\touter loop" << std::endl;
    file << "\t\tvertex" << std::setw(15) << (*i)->GetP1().x << std::setw(15)
         << (*i)->GetP1().y << std::setw(15) << (*i)->GetP1().z << std::endl;
    file << "\t\tvertex" << std::setw(15) << (*i)->GetP2().x << std::setw(15)
         << (*i)->GetP2().y << std::setw(15) << (*i)->GetP2().z << std::endl;
    file << "\t\tvertex" << std::setw(15) << (*i)->GetP3().x << std::setw(15)
         << (*i)->GetP3().y << std::setw(15) << (*i)->GetP3().z << std::endl;
    file << "\tendloop" << std::endl << "endfacet" << std::endl;

    i++;
  }

  file << "endsolid" << body_name << std::endl;

  file.close();
}

Point Geometry::MassCenter() {
  Point mc(0, 0, 0);

  double Stot = 0;

  std::deque<Polygon*>::iterator i = poligons.begin();

  while (i != poligons.end()) {
    Point a = (*i)->GetP2() - (*i)->GetP1();
    Point b = (*i)->GetP3() - (*i)->GetP2();

    double S = 0.5 * a.Cross(b).Mod();

    Stot += S;
    mc += S * ((*i)->GetP1() + (*i)->GetP2() + (*i)->GetP3()) / 3;

    i++;
  }

  if (Stot > 0) {
    mc /= Stot;
  }

  return mc;
}

void Geometry::Move(Point a) {
  std::deque<Polygon*>::iterator i = poligons.begin();

  while (i != poligons.end()) {
    (*i)->Move(a);
    i++;
  }
}

void Geometry::Scale(double e) {
  std::deque<Polygon*>::iterator i = poligons.begin();

  while (i != poligons.end()) {
    (*i)->Scale(e);

    i++;
  }
}

bool Geometry::IsInnerPoint(Point test_point) const {
  size_t n = 0;
  bool res = false;

  auto poligon_iter = poligons.begin();

  while (poligon_iter != poligons.end()) {
    /*if(dist == (*poligon_iter)->DistanceToPoint(test_point))
    {
            res = res;
    }*/
    /*if(dist > (*poligon_iter)->DistanceToPoint(test_point))
    {
            dist = (*poligon_iter)->DistanceToPoint(test_point);
            poligon_col = poligon_iter;
    }*/

    Point gmt = (*poligon_iter)->GetGmt();
    Point asd = test_point - (*poligon_iter)->GetGmt();
    if ((*poligon_iter)->GetNormal() *
            (test_point - (*poligon_iter)->GetGmt()) <
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

/*bool Geometry::IsInnerPoint(Point test_point)
{
        deque<Polygon*>::iterator poligon_col_iterator_x = poligons.end();
        deque<Polygon*>::iterator poligon_col_iterator_y = poligons.end();
        deque<Polygon*>::iterator poligon_col_iterator_z = poligons.end();

        Point collision_pstn_x(0, 0, 0);
        Point collision_pstn_y(0, 0, 0);
        Point collision_pstn_z(0, 0, 0);

        deque<Polygon*>::iterator poligon_iterator = poligons.begin();

        double dttx = 100000;
        double dtty = 100000;
        double dttz = 100000;

        while(poligon_iterator != poligons.end())
//ищем полигон с которым будет соударятся виртуальная частица
        {
                Point poligon_gmt = (*poligon_iterator)->GetGmt();

                Point a = (*poligon_iterator)->GetP2() -
(*poligon_iterator)->GetP1(); Point b = (*poligon_iterator)->GetP3() -
(*poligon_iterator)->GetP2(); Point c = (*poligon_iterator)->GetP1() -
(*poligon_iterator)->GetP3();

                double Ax = Point(1, 0, 0) * (*poligon_iterator)->GetNormal();
                double tcx = 0;
                if(Ax != 0) tcx	 = ((*poligon_iterator)->GetP1() - test_point)
* (*poligon_iterator)->GetNormal() / Ax;

                if(tcx < dttx && tcx > -dttx)
                {
                        Point collision_pstn = test_point + Point(1, 0, 0) *
tcx;

                        Point d1 = (*poligon_iterator)->GetP1() -
collision_pstn; Point d2 = (*poligon_iterator)->GetP2() - collision_pstn; Point
d3 = (*poligon_iterator)->GetP3() - collision_pstn;

                        if( (*poligon_iterator)->GetNormal() * d1.Cross(a) >
0. &&
                                (*poligon_iterator)->GetNormal() *
d2.Cross(b) > 0. &&
                                (*poligon_iterator)->GetNormal() *
d3.Cross(c) > 0.)
                        {
                                poligon_col_iterator_x = poligon_iterator;
                                collision_pstn_x = 	collision_pstn;
                                dttx = tcx;
                        }
                }


                double Ay = Point(0, 1, 0) * (*poligon_iterator)->GetNormal();
                double tcy = ((*poligon_iterator)->GetP1() - test_point) *
(*poligon_iterator)->GetNormal() / Ay;

                if(tcy < dtty && tcy > -dtty)
                {
                        Point collision_pstn = test_point + Point(1, 0, 0) *
tcy;

                        Point d1 = (*poligon_iterator)->GetP1() -
collision_pstn; Point d2 = (*poligon_iterator)->GetP2() - collision_pstn; Point
d3 = (*poligon_iterator)->GetP3() - collision_pstn;

                        if( (*poligon_iterator)->GetNormal() * d1.Cross(a) >
0. &&
                                (*poligon_iterator)->GetNormal() *
d2.Cross(b) > 0. &&
                                (*poligon_iterator)->GetNormal() *
d3.Cross(c) > 0.)
                        {
                                poligon_col_iterator_y = poligon_iterator;
                                collision_pstn_y = 	collision_pstn;
                                dtty = tcy;
                        }
                }

                double Az = Point(0, 0, 1) * (*poligon_iterator)->GetNormal();
                double tcz = ((*poligon_iterator)->GetP1() - test_point) *
(*poligon_iterator)->GetNormal() / Az;

                if(tcz < dttz && tcz > -dttz)
                {
                        Point collision_pstn = test_point + Point(1, 0, 0) *
tcz;

                        Point d1 = (*poligon_iterator)->GetP1() -
collision_pstn; Point d2 = (*poligon_iterator)->GetP2() - collision_pstn; Point
d3 = (*poligon_iterator)->GetP3() - collision_pstn;

                        if( (*poligon_iterator)->GetNormal() * d1.Cross(a) >
0. &&
                                (*poligon_iterator)->GetNormal() *
d2.Cross(b) > 0. &&
                                (*poligon_iterator)->GetNormal() *
d3.Cross(c) > 0.)
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
                if((*poligon_col_iterator_x)->GetNormal() * (collision_pstn_x -
test_point) < 0)
                {
                        return true;
                }
        }

        if(poligon_col_iterator_y != poligons.end())
        {
                if((*poligon_col_iterator_y)->GetNormal() * (collision_pstn_y -
test_point) < 0)
                {
                        return true;
                }
        }

        if(poligon_col_iterator_z != poligons.end())
        {
                if((*poligon_col_iterator_z)->GetNormal() * (collision_pstn_z -
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

/*void Geometry::create_cone(Point apex1, Point apex2, double H)
{
        poligons.clear();
        Point vertices[6];

        vertices[0] = apex1;
        vertices[1] = apex2;
        vertices[2] = Point(apex1.x, apex1.y, apex2.z);
        vertices[3] = Point(apex2.x, apex2.y, apex1.z);
        vertices[4] = apex2 - Point(0, 0, H);
        vertices[5] = vertices[3] - Point(0, 0, H);

        Polygon *temp_poligon = new Polygon;
        temp_poligon->Set(vertices[0], vertices[1], vertices[2], (vertices[1] -
vertices[0]).Cross(vertices[2] - vertices[1]).Normalize());
        poligons.push_back(temp_poligon);

        temp_poligon = new Polygon;
        temp_poligon->Set(vertices[0], vertices[1], vertices[3], (vertices[3] -
vertices[1]).Cross(vertices[1] - vertices[0]).Normalize());
        poligons.push_back(temp_poligon);

        temp_poligon = new Polygon;
        temp_poligon->Set(vertices[2], vertices[1], vertices[4], (vertices[1] -
vertices[2]).Cross(vertices[4] - vertices[2]).Normalize());
        poligons.push_back(temp_poligon);

        temp_poligon = new Polygon;
        temp_poligon->Set(vertices[0], vertices[3], vertices[5], (vertices[5] -
vertices[0]).Cross(vertices[3] - vertices[0]).Normalize());
        poligons.push_back(temp_poligon);

        temp_poligon = new Polygon;
        temp_poligon->Set(vertices[0], vertices[2], vertices[4], (vertices[0] -
vertices[2]).Cross(vertices[4] - vertices[2]).Normalize());
        poligons.push_back(temp_poligon);

        temp_poligon = new Polygon;
        temp_poligon->Set(vertices[0], vertices[4], vertices[5], (vertices[5] -
vertices[0]).Cross(vertices[4] - vertices[0]).Normalize());
        poligons.push_back(temp_poligon);
}*/
void Geometry::CreateWedge(double x, double width, double length,
                           double alpha) {
  Point vertices[6];

  vertices[0].Set(x, 0, -(width / 2));
  vertices[1].Set(x + length, length * sin(alpha), -width / 2);
  vertices[2].Set(x + length, -length * sin(alpha), -width / 2);

  vertices[3].Set(x, 0, width / 2);
  vertices[4].Set(x + length, length * sin(alpha), width / 2);
  vertices[5].Set(x + length, -length * sin(alpha), width / 2);

  Polygon* temp_poligon = new Polygon;
  temp_poligon->Set(vertices[0], vertices[1], vertices[2],
                    -1 * (vertices[0] - vertices[1])
                             .Cross(vertices[2] - vertices[1])
                             .Normalize());
  temp_poligon->flux = 0;
  temp_poligon->force = Point(0, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->Set(vertices[3], vertices[4], vertices[5],
                    -1 * (vertices[4] - vertices[3])
                             .Cross(vertices[5] - vertices[4])
                             .Normalize());
  temp_poligon->flux = 0;
  temp_poligon->force = Point(0, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->Set(vertices[0], vertices[1], vertices[3],
                    -1 * (vertices[1] - vertices[0])
                             .Cross(vertices[3] - vertices[0])
                             .Normalize());
  temp_poligon->flux = 0;
  temp_poligon->force = Point(0, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->Set(vertices[1], vertices[4], vertices[3],
                    -1 * (vertices[4] - vertices[1])
                             .Cross(vertices[3] - vertices[1])
                             .Normalize());
  temp_poligon->flux = 0;
  temp_poligon->force = Point(0, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->Set(vertices[0], vertices[2], vertices[3],
                    -1 * (vertices[3] - vertices[0])
                             .Cross(vertices[2] - vertices[0])
                             .Normalize());
  temp_poligon->flux = 0;
  temp_poligon->force = Point(0, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->Set(vertices[2], vertices[3], vertices[5],
                    -1 * (vertices[3] - vertices[2])
                             .Cross(vertices[5] - vertices[2])
                             .Normalize());
  temp_poligon->flux = 0;
  temp_poligon->force = Point(0, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->Set(vertices[1], vertices[2], vertices[5],
                    -1 * (vertices[2] - vertices[1])
                             .Cross(vertices[4] - vertices[1])
                             .Normalize());
  temp_poligon->flux = 0;
  temp_poligon->force = Point(0, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->Set(vertices[1], vertices[5], vertices[4],
                    -1 * (vertices[5] - vertices[1])
                             .Cross(vertices[4] - vertices[1])
                             .Normalize());
  temp_poligon->flux = 0;
  temp_poligon->force = Point(0, 0, 0);
  poligons.push_back(temp_poligon);
}

void Geometry::Fragment(double Lmax) {
  double max_Lmax = 10000.;
  while (max_Lmax > Lmax) {
    max_Lmax = 0;
    std::deque<Polygon*>::iterator poligon_iter = poligons.begin();

    while (poligon_iter != poligons.end()) {
      if ((*poligon_iter)->GetLmax() > Lmax) {
        max_Lmax > (*poligon_iter)->GetLmax()
            ? max_Lmax = max_Lmax
            : max_Lmax = (*poligon_iter)->GetLmax();
        std::pair<Polygon*, Polygon*> new_poligon = (*poligon_iter)->Divide();
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

std::pair<Point, Point> Geometry::Size() {
  Point p1(0, 0, 0), p2(0, 0, 0);

  std::deque<Polygon*>::iterator poligon_iter = poligons.begin();

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

  return pair<Point, Point>(p2, p1);
}

void Geometry::CreatePyramid(double x, double width, double length, double H) {}

double Geometry::DistanceToPoint(Point p) {
  double dist = 100000.;

  std::deque<Polygon*>::iterator poligon_iter = poligons.begin();

  while (poligon_iter != poligons.end()) {
    dist = std::min(dist, (*poligon_iter)->DistanceToPoint(p));
    poligon_iter++;
  }

  return dist;
}

void Geometry::CreateCube(double x, double width, double length, double H) {
  Point p[8];
  for (int i = 0; i < 4; i++) {
    p[i].x = x;
    p[i + 4].x = x + length;
  }

  p[0].y = p[3].y = p[4].y = p[7].y = -H / 2;
  p[1].y = p[2].y = p[5].y = p[6].y = H / 2;

  p[0].z = p[1].z = p[4].z = p[5].z = width / 2;
  p[2].z = p[3].z = p[6].z = p[7].z = -width / 2;

  Polygon* temp_poligon;

  temp_poligon = new Polygon;
  temp_poligon->p1 = p[0];
  temp_poligon->p2 = p[1];
  temp_poligon->p3 = p[2];
  temp_poligon->normal = Point(-1, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->p1 = p[0];
  temp_poligon->p2 = p[2];
  temp_poligon->p3 = p[3];
  temp_poligon->normal = Point(-1, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->p1 = p[4];
  temp_poligon->p2 = p[5];
  temp_poligon->p3 = p[6];
  temp_poligon->normal = Point(1, 0, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->p1 = p[4];
  temp_poligon->p2 = p[6];
  temp_poligon->p3 = p[7];
  temp_poligon->normal = Point(1, 0, 0);
  poligons.push_back(temp_poligon);

  /////////////////////////////////////////

  temp_poligon = new Polygon;
  temp_poligon->p1 = p[1];
  temp_poligon->p2 = p[2];
  temp_poligon->p3 = p[6];
  temp_poligon->normal = Point(0, 1, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->p1 = p[1];
  temp_poligon->p2 = p[6];
  temp_poligon->p3 = p[5];
  temp_poligon->normal = Point(0, 1, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->p1 = p[0];
  temp_poligon->p2 = p[3];
  temp_poligon->p3 = p[7];
  temp_poligon->normal = Point(0, -1, 0);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->p1 = p[0];
  temp_poligon->p2 = p[7];
  temp_poligon->p3 = p[4];
  temp_poligon->normal = Point(0, -1, 0);
  poligons.push_back(temp_poligon);

  /////////////////////////////////////////

  temp_poligon = new Polygon;
  temp_poligon->p1 = p[0];
  temp_poligon->p2 = p[1];
  temp_poligon->p3 = p[5];
  temp_poligon->normal = Point(0, 0, -1);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->p1 = p[0];
  temp_poligon->p2 = p[5];
  temp_poligon->p3 = p[4];
  temp_poligon->normal = Point(0, 0, -1);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->p1 = p[3];
  temp_poligon->p2 = p[2];
  temp_poligon->p3 = p[6];
  temp_poligon->normal = Point(0, 0, 1);
  poligons.push_back(temp_poligon);

  temp_poligon = new Polygon;
  temp_poligon->p1 = p[3];
  temp_poligon->p2 = p[6];
  temp_poligon->p3 = p[7];
  temp_poligon->normal = Point(0, 0, 1);
  poligons.push_back(temp_poligon);
}

int Geometry::FixPolygons() {
  int res = 0;
  deque<Polygon*>::iterator poligon_iter = poligons.begin();

  while (poligon_iter != poligons.end()) {
    if ((*poligon_iter)->Fix()) res++;
    poligon_iter++;
  }

  return res;
}

void Geometry::ReverseNormals() {
  deque<Polygon*>::iterator poligon_iter = poligons.begin();

  while (poligon_iter != poligons.end()) {
    (*poligon_iter)->normal = (*poligon_iter)->normal * -1.;
    poligon_iter++;
  }
}
}  // namespace mc3d
