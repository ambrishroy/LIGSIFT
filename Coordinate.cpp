#include "Coordinate.h"

Coordinate::Coordinate(void):
  x(0.0),
  y(0.0),
  z(0.0)
{
}
/*************************************/         
Coordinate::Coordinate(double x, double y, double z):
  x(x),
  y(y),
  z(z)
{
};
/*************************************/
std::ostream&
operator<< (std::ostream& os, const Coordinate& A){
  os << "(" << A.x << "," << A.y << "," << A.z << ")";
  return os;
};
