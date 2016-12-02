#ifndef COMMON_H
#define COMMON_H
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string.h>
#include <stdio.h>
using namespace std;

#ifndef TRUE
#define TRUE            1
#endif

#ifndef FALSE
#define FALSE           0
#endif

#ifndef epsilon
#define epsilon       1e-7
#endif

#ifndef PI
#define PI		3.14159265
#endif

#ifndef pi
#define pi		3.14159265
#endif

#ifndef ROUND
#define ROUND(x)      ((int) ((x) + 0.5))
#endif

#ifndef MAX
#define MAX(x, y)     ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x, y)     ((x) > (y) ? (y) : (x))
#endif

#ifndef GCI
#define GCI		2.828427125
#endif

#ifndef GCI2
#define GCI2		7.999999999
#endif

#ifndef GCI3
#define GCI3		2.42433 
#endif


#ifndef SPR
#define SPR		(4/3)*3.14159265
#endif

#ifndef OCON
#define OCON		3.14159265/12
#endif

#ifndef ALPHAC
#define ALPHAC		2.414824925
#endif

#ifndef NULL
#define NULL		0
#endif

#ifndef TAU
#define TAU		1e-12
#endif

#ifndef INF
#define INF		HUGE_VAL
#endif

#ifndef TRUE
#define TRUE		1
#endif

#ifndef FALSE
#define FALSE		0
#endif

#ifndef kB
#define kB		1.987206504191549E-003
#endif 

#ifndef kBT
#define kBT		298.15*1.987206504191549E-003 
#endif

#ifndef RAD_TO_DEG
#define RAD_TO_DEG	57.2957
#endif

#ifndef MY_4PI
#define MY_4PI 		12.56637061435917295384
#endif

#ifndef KCAL_TO_KJ
#define KCAL_TO_KJ     4.1868
#endif


#ifndef EC
#define EC               14.3996
#endif

double rand(const double, const double);
double distsq(double x[3], double y[3]);
double dist(double x[3], double y[3]);
double transform(double t[3], double r[3][3], double *x, unsigned i);
void   PVal(const double val, const double v1, const double v2, const int SC, const int Opt, double *SIM, double *PVal);
string cleanS(string S);

void rotateX(double**, double**, const double, const double, unsigned);
void rotateY(double**, double**, const double, const double, unsigned);
void rotateZ(double**, double**, const double, const double, unsigned);
void RotationMatrix(double M[3][3], const double, const double, double, double);
void RandRotationMatrix(double R[3][3]);
void AlignVector(double x[3], double y[3], double R[3][3]);
void crossProduct (double x[3], double y[3], double z[3]);
double dotProduct (double x[3], double y[3]);
double vectorAngle(double x[3], double y[3]);


string GetFilename(string  &  s);
string GetDirectory(string &  s);
string removeExtension(string &  s);
string toLower(string s);



template <class A> void CreateArray(A *** array, int Narray1, int Narray2)
{
  *array=new A* [Narray1];
  for(int i=0; i<Narray1; i++) *(*array+i)=new A [Narray2];
};

template <class A> void DeleteArray(A *** array, int Narray)
{
  for(int i=0; i<Narray; i++)
    if(*(*array+i)) delete [] *(*array+i);
  if(Narray) delete [] (*array);
  (*array)=NULL;
};

typedef std::pair<string, double> MyHash;
//Usage:std::vector< MyHash > v;
////std::sort(v.begin(), v.end(), SortKey());
struct SortKey {
 bool operator() (const MyHash& a, const MyHash& b) const {
 return a.first < b.first;
 };
};
//Usage:std::vector< MyHash > v;
//std::sort(v.begin(), v.end(), SortValue());
struct SortValue{
  bool operator() (const MyHash& a, const MyHash& b) const{
  return a.second < b.second;
  };
};
struct RevSortValue{
  bool operator() (const MyHash& a, const MyHash& b) const{
  return a.second > b.second;
  };
};




#endif // COMMON_H

