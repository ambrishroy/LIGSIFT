#ifndef RMSD_H
#define RMSD_H

#include <math.h>
//using namespace std;


bool Kabsch(double **x, double **y, int n, int mode, double *rms, double t[3], double u[3][3]);
 


#endif  // RMSD_H
