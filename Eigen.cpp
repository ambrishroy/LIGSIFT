#include "Eigen.h"
using namespace std;
/**************************************************************
     Fast Calculation of eigenvalue and eignevectors using
            Moment of inertia tensor matrix
---------------------------------------------------------------
     Cov[3][3]-MI Tensor matrix 
     EigVal[3]-eigenvalue 
     EigVec[3][3]-eigenvector
     E(1)=<x^2>; E(1)=<y^2>; E(1)=<z^2> in the rotated system
***************************************************************/
bool Eigen(double A[3][3], double EigVec[3][3], double EigVal[3]){
  
  double p1=0, p2=0, p3=0, p4=0, p5=0, p6=0, p7=0, p8=0, p9=0, p10=0, p11=0;
  double ap5=0, x=0, y=0, fnorm1=0;

  p1=-1;
  p2 =A[0][0] + A[1][1] + A[2][2];
  p3=-A[0][0]*A[1][1] - A[0][0]*A[2][2] - A[1][1]*A[2][2] \
    + A[1][2]*A[1][2] + A[0][2]*A[0][2] + A[0][1]*A[0][1];
  p4= A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] \
    + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0]	  \
    - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2];
  p5 = (-(1.0/3) * pow((p2/p1),2) + (p3/p1))/3;
  ap5= sqrt(-p5);
  p6 = ( (2.0/27) * pow((p2/p1),3) - (1.0/3)*(p2*p3/pow(p1,2)) + (p4/p1) )/2;
  p11= p2/(3*p1);
    
  if((abs(-p6/sqrt(-pow(p5,3)))) > 1){   
    p7 = 1;
  }
  else{
    p7=acos(-p6/sqrt(-pow(p5,3)));
  }
  p8 =  2*ap5*cos(p7/3.0);
  p9 = -2*ap5*cos((p7+PI)/3.0);
  p10= -2*ap5*cos((p7-PI)/3.0);
  EigVal[0]= p8-p11;               //eigenvalue
  EigVal[1]= p9-p11;               //eigenvalue
  EigVal[2]= p10-p11;              //eigenvalue

  /***Normalized eigenvectors ****/
  for(int i=0; i < 3; ++i){
    fnorm1=A[1][0]*A[0][1]-(A[0][0]-EigVal[i])*(A[1][1]-EigVal[i]);
      x=((A[1][1]-EigVal[i])*A[0][2]-A[0][1]*A[1][2])/fnorm1;
      y=((A[0][0]-EigVal[i])*A[1][2]-A[1][0]*A[0][2])/fnorm1;
      EigVec[i][2]=1/sqrt( x*x + y*y + 1);
      EigVec[i][0]=x*EigVec[i][2];
      EigVec[i][1]=y*EigVec[i][2];
  }
  return true;  
}

