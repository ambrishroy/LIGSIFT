#ifndef COORDINATE_H
#define COORDINATE_H

#include <iostream>

class Coordinate{
 public:
  double		x;
  double		y;
  double		z;
  Coordinate(void);
  Coordinate(double, double, double);
};


/*************************************/
std::ostream& operator<< (std::ostream&, const Coordinate&);

#endif // COORDINATE_H
