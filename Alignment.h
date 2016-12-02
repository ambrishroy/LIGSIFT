#ifndef ALIGNMENT_H
#define ALIGNMENT_H
#include <map>
#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include "Pharmacophores.h"
#include "Coordinate.h"
#include "Common.h"
#include "RMSD.h"
#include "LAP.h"
#include "Eigen.h"


void Overlap(vector<PharmacophorePoint>&,vector<PharmacophorePoint>&, vector<double>& Scores, vector<double>& Vol, std::map<int, int>& BestAln, int Opt);

// Cost Matrix for finding optimal alignment
void FillMatrix(double **M, vector<PharmacophorePoint>&, vector<PharmacophorePoint>&, int Opt);
// Overlap Optimization
void OverlapScore(vector<PharmacophorePoint>&, vector<PharmacophorePoint>&, double *Shape, double *Chem, double *Charge, double **M, double *U, double *V, int *Row, int *Col, double **VEC1, double **VEC2,std::map<int, int>& Aln, int Opt);
// Overlap Score
double Score(double ShapeS, double ChemS, double V1shape, double V2shape, double V1chem, double V2chem, int Opt);
// Find directionality of pharmacophoric points
double norm(Coordinate & p);
// Find directionality cosine similarity of pharmacophoric points
double cosine(Coordinate& p1, Coordinate& p2);
// Rotate and Traslate
void transformation(vector<PharmacophorePoint>&, double u[3][3],double t[3]);

#endif  // ALIGNMENT_H
