//
//  LAP.h
//  pickPocket
//
//  Created by Ambrish Roy on 11/11/15.
//  Copyright (c) 2015 Ambrish Roy. All rights reserved.
//

#ifndef pickPocket_LAP_h
#define pickPocket_LAP_h
#include <stdio.h>

/*************** CONSTANTS  *******************/
#define BIG 100000


/*************** TYPES      *******************/
typedef int row;
typedef int col;
typedef double cost;


/*************** FUNCTIONS  *******************/
extern double    lap(int dim, double **assigncost, int *rowsol, int *colsol, double *u, double *v);

extern void checklap(int dim, double **assigncost, int *rowsol, int *colsol, double *u, double *v);

#endif
