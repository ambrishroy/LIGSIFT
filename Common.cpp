#include "Common.h"


double distsq(double x[3], double y[3]){
  double d1=x[0]-y[0];
  double d2=x[1]-y[1];
  double d3=x[2]-y[2];
  
  return (d1*d1 + d2*d2 + d3*d3);
}
/*************************************/
double dist(double x[3], double y[3]){
  double d1=x[0]-y[0];
  double d2=x[1]-y[1];
  double d3=x[2]-y[2];	
  return sqrt(d1*d1 + d2*d2 + d3*d3);
}
/*************************************/
double transform(double t[3], double r[3][3], double *x, unsigned i){
  double coor= t[i] + (r[i][0]*x[0]) + (r[i][1]*x[1]) + (r[i][2]*x[2]);
  return coor;
}
/*************************************/
void rotateX(double** oc, double** nc, const double c, const double s, unsigned N){
  for (unsigned int i = 0; i < N; ++i){
    nc[i][0] =  oc[i][0];
    nc[i][1] =  oc[i][1] * c + oc[i][2] * s;
    nc[i][2] = -oc[i][1] * s + oc[i][2] * c;
  }
}
/*************************************/
void rotateY(double** oc, double** nc, const double c, const double s, unsigned N){
  for (unsigned int i = 0; i < N; ++i){
    nc[i][0] =  oc[i][0] * c + oc[i][2] * s;
    nc[i][1] =  oc[i][1];
    nc[i][2] = -oc[i][0] * s + oc[i][2] * c;
  }
}
/*************************************/
void rotateZ(double** oc, double** nc, const double c, const double s, unsigned N){
  for (unsigned int i = 0; i < N; ++i){
    nc[i][0] =  oc[i][0] * c + oc[i][1] * s;
    nc[i][1] = -oc[i][0] * s + oc[i][1] * c;
    nc[i][2] =  oc[i][2];
  }
}
/*************************************/ 
void RotationMatrix(double rotationMatrix[3][3], double angle, double u, double v, double w){
  //(u, v, w) is a vector that specifies the axis about which the object is to be rotated.
  double L = (u*u + v * v + w * w);
  angle = angle * PI / 180.0; //converting to radian value
  double u2 = u * u;
  double v2 = v * v;
  double w2 = w * w; 
   
  rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(angle)) / L;
  rotationMatrix[0][1] = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
  rotationMatrix[0][2] = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
  
  rotationMatrix[1][0] = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
  rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(angle)) / L;
  rotationMatrix[1][2] = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
  
  rotationMatrix[2][0] = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
  rotationMatrix[2][1] = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
  rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(angle)) / L;
} 
/*************************************/ 
void RandRotationMatrix(double rotationMatrix[3][3]){
  double x1= rand(0,1);
  double x2= rand(0,1);
  double x3= rand(0,1);
  double rotMatrix[3][3]={0};
  double vv[3][3]={0};
  double H[3][3] ={0};
  rotMatrix[0][0] = cos(2.0*PI*x1);
  rotMatrix[0][1] = sin(2.0*PI*x1);
  rotMatrix[0][2] = 0.0;     
  rotMatrix[1][0] = -1.0 * sin(2.0*PI*x1);
  rotMatrix[1][1] = cos(2.0*PI*x1);
  rotMatrix[1][2] = 0.0;     
  rotMatrix[2][0] = 0.0;
  rotMatrix[2][1] = 0.0;
  rotMatrix[2][2] = 1.0;
  double v[3]     = { sqrt(x3)*cos(2.0*pi*x2), sqrt(x3)*sin(2.0*pi*x2), sqrt(1.0 - x3)};
  
  vv[0][0] = v[0] * v[0];
  vv[0][1] = v[0] * v[1];
  vv[0][2] = v[0] * v[2];  
  vv[1][0] = v[1] * v[0];
  vv[1][1] = v[1] * v[1];
  vv[1][2] = v[1] * v[2];  
  vv[2][0] = v[2] * v[0];
  vv[2][1] = v[2] * v[1];
  vv[2][2] = v[2] * v[2];  
  
  H[0][0] = 1.0 - 2*vv[0][0];
  H[0][1] = 0.0 - 2*vv[0][1];
  H[0][2] = 0.0 - 2*vv[0][2];
  
  H[1][0] = 0.0 - 2*vv[1][0];
  H[1][1] = 1.0 - 2*vv[1][1];
  H[1][2] = 0.0 - 2*vv[1][2];
  
  H[2][0] = 0.0 - 2*vv[2][0];
  H[2][1] = 0.0 - 2*vv[2][1];
  H[2][2] = 1.0 - 2*vv[2][2];
  
  for(int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j ){
      rotationMatrix[i][j] = 0.0;
      for(int k = 0; k < 3; ++k){
	rotationMatrix[i][j] = (rotationMatrix[i][j] + H[i][k] * rotMatrix[k][j]);
      }
      rotationMatrix[i][j] *= -1.0;
    }	
  }
}

/*************************************/
void crossProduct(double v1[3], double v2[3], double tv3[3]){
    tv3[0] =   v1[1]*v2[2] - v1[2]*v2[1];
    tv3[1] = - v1[0]*v2[2] + v1[2]*v2[0];
    tv3[2] =   v1[0]*v2[1] - v1[1]*v2[0];
}
/*************************************/
double dotProduct(double v1[3], double v2[3]){
  return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

/*************************************/
double vectorAngle (double v1[3], double v2[3]){
    double dp;
    double l1=sqrt((v1[0]*v1[0]) + (v1[1]*v1[1]) + (v1[2]*v1[2]));
    double l2=sqrt((v2[0]*v2[0]) + (v2[1]*v2[1]) + (v2[2]*v2[2]));
    dp = dotProduct(v1,v2)/(l1 * l2);
    if (dp < -0.999999)
      dp = -0.9999999;

    if (dp > 0.9999999)
      dp = 0.9999999;

    return((RAD_TO_DEG * acos(dp)));
}
/*************************************/
void AlignVector(double CordVec1[3], double CordVec2[3], double rotationMatrix[3][3]){
   double CordVec3[3]={0}, CordVec4[3]={0}, RM[3][3]={0};

   crossProduct(CordVec1,CordVec2,CordVec3);   
   double VecLength=sqrt((CordVec3[0]*CordVec3[0]) + (CordVec3[1]*CordVec3[1]) + (CordVec3[2]*CordVec3[2]));
   //Normalize Vec3
   if(VecLength>0){
    CordVec3[0] /=VecLength; CordVec3[1] /=VecLength; CordVec3[2] /=VecLength;
    double angle=vectorAngle(CordVec1,CordVec2);
    //cout << angle << " " << CordVec3[0] << " " << CordVec3[1] << " " << CordVec3[2] << endl;
    RotationMatrix(RM, -angle, CordVec3[0], CordVec3[1], CordVec3[2]);
    rotationMatrix[0][0]=RM[0][0];
    rotationMatrix[0][1]=RM[0][1]; 
    rotationMatrix[0][2]=RM[0][2]; 
    rotationMatrix[1][0]=RM[1][0];  
    rotationMatrix[1][1]=RM[1][1];
    rotationMatrix[1][2]=RM[1][2];
    rotationMatrix[2][0]=RM[2][0];  
    rotationMatrix[2][1]=RM[2][1];
    rotationMatrix[2][2]=RM[2][2];
   }
   else{
    rotationMatrix[0][0]= 1;
    rotationMatrix[0][1]= 0;
    rotationMatrix[0][2]= 0;
    rotationMatrix[1][0]= 0;
    rotationMatrix[1][1]= 1;
    rotationMatrix[1][2]= 0;
    rotationMatrix[2][0]= 0;
    rotationMatrix[2][1]= 0;
    rotationMatrix[2][2]= 1;
   }
}
 
/*************************************/
std::string cleanS(string S ){
  string C;
  for(int j=0; j<S.length(); j++){
    if(S[j]==' '){
      continue;
    }
    else if(S[j]=='\n'){
      continue;
    }
    else if(S[j]=='\t'){
      continue;
    }
    else{
      C += S[j];
    }
  }
  return C;
}

/*************************************/
std::string GetFilename(string& s){
  char sep = '/';
  size_t i = s.rfind(sep, s.length( ));
  if (i != string::npos) {
    return(s.substr(i+1, s.length( ) - i));
  }  
  return("");  
}
/*************************************/
std::string GetDirectory(string& path ){
  string dir=path.substr( 0, path.find_last_of( '\\' ) +1 );
  return dir;
}
/*************************************/
std::string removeExtension(string& filename) {
  size_t lastdot = filename.find_last_of(".");
  if (lastdot == std::string::npos) return filename;
  return filename.substr(0, lastdot); 
}
/*************************************/
string toLower(string strr){	
  char str[100];
  string ret;
  strcpy(str,strr.c_str());
  int differ = 'A'-'a';
  char ch;
  int ii = strlen(str);
  for (int i=0; i <ii;i++){
    strncpy(&ch,str+i,1);
    if (ch>='A' && ch<='Z'){
      ch = ch-differ;
      memcpy(str+i,&ch,1);
    }
  }
  ret = str;
  return ret;
}
/*************************************/
// function to return a random number between 'start' to 'end'
double rand(const double start, const double end){
  return (end-start) * (double(rand()) / RAND_MAX) + start;
}
/*************************************/
void
PVal(const double val, const double v1, const double v2, const int SC, const int Optim, double *SIM, double *pval){
  // Optim 0 : Shape based ; 1: Shape + Chemistry based ; 2: Chemistry
  // SC    0 : Shape Similarity ; 1: Chemical Similarity
  double vs   = v1 + v2;
  double vd   = abs(v1-v2);
  double con = 0.0;      
  double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
  double k1= 0.0, k2= 0.0, k3= 0.0, k4= 0.0, k5=0.0; 
  double k6= 0.0, k7= 0.0, k8= 0.0, k9=0.0, k10=0.0;
  double SIMp =0.0; double S0 =0.0; double mue=0; double sig=0.0;
  
  if((Optim == 0) && (SC==0)){
    con = 0.49;
    a = 1.452;        b = 0.278;       c = 0.279;       d=-0.669;
    k1= 2.063959e-01;k2=-2.434431e-05; k3=-2.545694e-05; k4= 4.222393e-02; k5= 1.225259e-05;	
    k6= 0.35817104;  k7= 0.0463052;    k8=0.0462624;     k9=-0.1294002;
  }
  else if((Optim == 0) && (SC==1)){
    con = 0.35;
    a = 1.141;      b= 0.259;       c= 0.259;       d=-0.600;
    k1=-4.835488e-01;k2=-1.136479e-04; k3=-1.136376e-04; k4= 1.350711e-01; k5= 1.547031e-05;	
    k6= 0.38259663;  k7= 0.05740400;   k8=0.05742364;    k9=-0.15130673;
  }
  else if((Optim == 2) && (SC==0)){	
    con = 0.43;
    a = 1.251; b = 0.230; c= 0.230; d=-0.559;
    k1= 1.039311e-01;k2=-3.641794e-05; k3=-3.633077e-05; k4= 5.060486e-02; k5= 1.336170e-05;	
    k6= 0.28001779;  k7= 0.02844948;   k8=0.02846297;    k9=-0.08681664;
  }
  else if((Optim == 2) && (SC==1)){
    con = 0.4;	
    a = 1.325; b = 0.2574; c= 0.2574; d=-0.6193; 
    k1= 1.101308e-01;k2=-3.469808e-05; k3=-3.469274e-05; k4= 4.809996e-02; k5= 1.382582e-05;	
    k6= 0.28753880;  k7= 0.02963015;   k8=0.02967537;    k9=-0.08995957;
  }
  else if((Optim == 1) && (SC==0)){
    con = 0.55;
    a=  1.621 ;        b= 0.273;         c= 0.273;         d=-0.681; 
    k1= 3.645590e-01;k2=-2.174805e-05; k3=-2.175689e-05; k4= 2.795199e-02; k5= 1.009904e-05;	
    k6= 0.20025365;  k7= 0.03226759;   k8=0.03228038;    k9=-0.08302792;
  }
  else if((Optim == 1) && (SC==1)){
    con = 0.38;
    a= 1.077; b= 0.202; c= 0.202; d=-0.487;
    k1=-1.407326e-01;k2=-6.373657e-05; k3=-6.373825e-05; k4= 8.241650e-02; k5= 1.091803e-05;	
    k6= 0.31108818;  k7= 0.03466335;   k8=0.03466188;    k9=-0.10165057;
  }
  SIMp = a + (b*log(v1)) + (c*log(v2)) + (d*log(vs));
  S0   = (SIMp-con)/(con-1); if(S0 < 0.001){S0=0.001;}
  *SIM  = (val+S0)/(1+S0);
  mue  = k1 + (k2*v1) + (k3*v2) + (k4*log(vs)) + (k5*vd);
  sig  = k6 + (k7*log(v1)) + (k8*log(v2)) + (k9*log(vs));
  double zs  =((*SIM-mue)/sig);
  *pval=1 - exp(-exp(-zs));
}


