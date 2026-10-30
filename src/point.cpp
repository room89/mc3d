// File: point.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 12.10.10.
// Description: Program for calculation of freemolecular flows.

#include "point.h"

#include <cmath>
#include <iostream>
#include <utils/utils.hpp>

namespace mc3d {
void Point::Print() {
  std::cout << "(" << x << "," << y << "," << this->z << ")" << std::endl;
}

std::ostream& operator<<(std::ostream& o, const Point c) {
  o << c.x << ", " << c.y << ", " << c.z;
  return o;
}

/*Point& Point::operator=(const Point &a)
{
        x = a.x;
        y = a.y;
        z = a.z;
        return *this;
}*/

Point Point::RandomPoint(
    Point size, utils::RandomStream& random)  // возврощает рандомную точку с коорд.: x от 0 до size.x, y от
                 // 0 до size.y, z от 0 до size.z
{
  double rnx = random.Uniform();
  double rny = random.Uniform();
  double rnz = random.Uniform();

  return Point(size.x * rnx, size.y * rny, size.z * rnz);
}

Point Point::RandomPoint(utils::RandomStream& random)  // возврощает рандомную точку с коорд.: x от 0 до
                            //  this->x, y от 0 до this->y, z от 0 до this->z
{
  double rnx = random.Uniform();
  double rny = random.Uniform();
  double rnz = random.Uniform();

  return Point(x * rnx, y * rny, z * rnz);
}

}  // namespace mc3d
