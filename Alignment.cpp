#include "Alignment.h"
using namespace std;



void 
Overlap(vector<PharmacophorePoint>& P1, vector<PharmacophorePoint>& P2, vector<double>& Overlap_Scores, vector<double>& Volumes, std::map<int,int>& BestAlignment,int Opt){
  // Initialization
  double W1= 0.7, W2=(1-W1); // Weights for calculating Tversky scores
  double STanimoto = 0, CTanimoto = 0, ChTanimoto = 0, TVshape = 0, TVchem=0, INTF=0;
  int npharm1=P1.size();
  int npharm2=P2.size();
  int npharm_max =2*MAX(npharm1, npharm2);
  if(npharm_max < 3)npharm_max=3;

  vector<PharmacophorePoint> QPharmacophores;
  vector<PharmacophorePoint> TPharmacophores;
  vector<PharmacophorePoint> TPharmacophores2;
  std::map<int,int> Alignment;
  int NEXCL1=0, NEXCL2=0;  
  
  
  Coordinate P1Center, P2Center;
  double V1shape(0.0), V2shape(0.0), V1chem(0.0), V2chem(0.0), maxoverlap(0.00);
  double V1charge(0.0), V2charge(0.0), INTERFACE(0.0);
  double CordVec1[3], CordVec2[3], CordVec3[3], CordVec4[3];  
  double RMSD; double RM[3][3]={0};
  
  
  double **MAT, **VEC1, **VEC2, **DIST1, **DIST2;
  double *U, *V;
  int    *Row,  *Col;
  
  CreateArray(&MAT, npharm_max, npharm_max);
  CreateArray(&DIST1, npharm_max, npharm_max);
  CreateArray(&DIST2, npharm_max, npharm_max);
  CreateArray(&VEC1, npharm_max, 3); CreateArray(&VEC2, npharm_max, 3);
  U   = new double[npharm_max]; V   = new double[npharm_max];
  Row = new int   [npharm_max]; Col = new int [npharm_max];
  QPharmacophores= P1; 

  
  P1Center.x =0; P1Center.y=0; P1Center.z=0;  
  for(int i=0; i <npharm1; ++i){       
    if(P1[i].includedBefore){continue;}
    double D=0;   
    double alpha1 = P1[i].alpha;    
    double v1     = GCI * pow((PI/alpha1), 1.5);      
    P1[i].vol=v1;
    
//    if(P1[i].func < 8){ V1chem     += v1; V1charge    += v1;}
    if(P1[i].func < 8){V1charge    += v1;}
    V1shape     += v1;  
     V1chem     += v1;   
    P1Center.x  += v1 * P1[i].point.x;
    P1Center.y  += v1 * P1[i].point.y;
    P1Center.z  += v1 * P1[i].point.z;    
  } 
  P1Center.x /= V1shape;
  P1Center.y /= V1shape;
  P1Center.z /= V1shape;
  

  double TM1[3][3]={0.00}, EVec1[3][3]={0}, EVal1[3]={0};
  for(int i=0; i < npharm1; i++){    
    P1[i].point.x -= P1Center.x;
    P1[i].point.y -= P1Center.y;
    P1[i].point.z -= P1Center.z;    
    if(P1[i].hasNormal){
      P1[i].normal.x -= P1Center.x;
      P1[i].normal.y -= P1Center.y;
      P1[i].normal.z -= P1Center.z;
    }
    //cout << i << " " << P1[i].func << " " << P1[i].vdw << endl; 
    if(P1[i].includedBefore){continue;}    
    double v1    = P1[i].vol;
    TM1[0][0] += v1 * P1[i].point.x * P1[i].point.x;
    TM1[0][1] += v1 * P1[i].point.x * P1[i].point.y;
    TM1[0][2] += v1 * P1[i].point.x * P1[i].point.z;
    TM1[1][1] += v1 * P1[i].point.y * P1[i].point.y;
    TM1[1][2] += v1 * P1[i].point.y * P1[i].point.z;
    TM1[2][2] += v1 * P1[i].point.z * P1[i].point.z;    
  }  
  TM1[1][0] = TM1[0][1];
  TM1[2][0] = TM1[0][2];
  TM1[2][1] = TM1[1][2];
  for(int k=0; k<3; k++)
    for(int l=0; l< 3; l++)
      TM1[k][l] /= V1shape;
  
  QPharmacophores = P1;
  Eigen(TM1, EVec1, EVal1);

  
  P2Center.x =0; P2Center.y=0; P2Center.z=0; 
  for(int j=0; j <npharm2; ++j){
    if(P2[j].includedBefore){continue;}
    double D=0;
    double alpha1 = P2[j].alpha;
    double v1     = GCI * pow((PI/alpha1), 1.5);     
    INTF +=P2[j].interactions * v1;
    P2[j].vol=v1;
//    if(P2[j].func < 8){ V2chem     += v1;V2charge    += v1;}   
    if(P2[j].func < 8){V2charge    += v1;}
    V2chem      += v1;
    V2shape     += v1;   
    P2Center.x  += v1 * P2[j].point.x;
    P2Center.y  += v1 * P2[j].point.y;
    P2Center.z  += v1 * P2[j].point.z;  
  }   
  P2Center.x /= V2shape;
  P2Center.y /= V2shape;
  P2Center.z /= V2shape;
 
  double TM2[3][3]={0.00}, EVec2[3][3]={0}, EVal2[3]={0.00};
  for(int j=0; j < npharm2; j++){   
    P2[j].point.x -= P2Center.x;
    P2[j].point.y -= P2Center.y;
    P2[j].point.z -= P2Center.z;
    if(P2[j].hasNormal){
      P2[j].normal.x -= P2Center.x;
      P2[j].normal.y -= P2Center.y;
      P2[j].normal.z -= P2Center.z;
    } 
    //cout << j << " " << P2[j].func << " " << P2[j].vdw << endl;
    if(P2[j].includedBefore){continue;}
    double v2    = P2[j].vol;   
    TM2[0][0] += v2 * P2[j].point.x * P2[j].point.x;
    TM2[0][1] += v2 * P2[j].point.x * P2[j].point.y;
    TM2[0][2] += v2 * P2[j].point.x * P2[j].point.z;
    TM2[1][1] += v2 * P2[j].point.y * P2[j].point.y;
    TM2[1][2] += v2 * P2[j].point.y * P2[j].point.z;
    TM2[2][2] += v2 * P2[j].point.z * P2[j].point.z;    
  }  
  TM2[1][0] = TM2[0][1];
  TM2[2][0] = TM2[0][2];
  TM2[2][1] = TM2[1][2];
  for(int k=0; k<3; k++)
    for(int l=0; l< 3; l++)
      TM2[k][l] /= V2shape;
  
  Eigen(TM2, EVec2, EVal2);

  double Rg1=0.00; double Rg2=0.0;
  for(int i=0; i<npharm1;++i){
    CordVec1[0]=P1[i].point.x;
    CordVec1[1]=P1[i].point.y;
    CordVec1[2]=P1[i].point.z;
    CordVec2[0]=0;
    CordVec2[1]=0;
    CordVec2[2]=0;
    double D1=distsq(CordVec1, CordVec2);
    Rg1 += D1;    
    for(int ii=i+1; ii<npharm1;++ii){
      CordVec2[0]=P1[ii].point.x;
      CordVec2[1]=P1[ii].point.y;
      CordVec2[2]=P1[ii].point.z;
      double D= distsq(CordVec1,CordVec2);    
      DIST1[i][ii]=D;
      DIST1[ii][i]=D;
    }
  }
  
  for(int i=0; i<npharm2;++i){
    CordVec1[0]=P2[i].point.x;
    CordVec1[1]=P2[i].point.y;
    CordVec1[2]=P2[i].point.z;    
    CordVec2[0]=0;
    CordVec2[1]=0;
    CordVec2[2]=0;
    double D1=distsq(CordVec1, CordVec2);
    Rg2 += D1;    
    for(int ii=i+1; ii<npharm2;++ii){
      CordVec2[0]=P2[ii].point.x;
      CordVec2[1]=P2[ii].point.y;
      CordVec2[2]=P2[ii].point.z;
      double D= distsq(CordVec1, CordVec2);      
      DIST2[i][ii]=D;
      DIST2[ii][i]=D;
    }
  }  
  Rg1 /=npharm1; Rg2 /= npharm2;
   
  double T[3]={0}; double R[3][3]={0}; int bestPAxis=0;;  
  double ShapeScore=0.0; double ChemScore=0; double ChargeScore=0.00; double INTFscore=0;
  double Delta=3.0; double maxR=0; int rotationStep=1;
 
  for(int pa1=0;pa1<3;++pa1){         
    CordVec1[0]=EVec1[pa1][0];
    CordVec1[1]=EVec1[pa1][1];
    CordVec1[2]=EVec1[pa1][2];        
    for(int pa2=0;pa2<3;++pa2){          
      CordVec2[0]=EVec2[pa2][0];
      CordVec2[1]=EVec2[pa2][1];
      CordVec2[2]=EVec2[pa2][2];           
      AlignVector(CordVec1, CordVec2, R);              
      TPharmacophores2.clear();
      TPharmacophores2=P2; 
      transformation(TPharmacophores2,R,T);
      OverlapScore(QPharmacophores, TPharmacophores2, &ShapeScore, &ChemScore, &ChargeScore, MAT, U, V, Row, Col, VEC1, VEC2, Alignment,Opt);    
      double S = Score(ShapeScore,ChemScore,V1shape,V2shape,V1chem,V2chem,Opt);
      if(S > maxoverlap){
	
	bestPAxis= pa2;
	TPharmacophores=TPharmacophores2;   
	BestAlignment.clear();
	for (std::map<int,int>::iterator it=Alignment.begin(); it!=Alignment.end(); ++it){
	  BestAlignment.insert(std::make_pair(it->first,it->second));	 
	}
	maxoverlap= S;
	STanimoto = ShapeScore/(V1shape + V2shape - ShapeScore);
	CTanimoto = ChemScore/(V1chem + V2chem - ChemScore);
	TVshape   = (2*ShapeScore)/((W1*V1shape) + (W2*V2shape) + ShapeScore);
	TVchem    = (2*ChemScore)/((W1*V1chem)  + (W2*V2chem) + ChemScore);
	ChTanimoto= ChargeScore/(V1charge + V2charge - ChargeScore);	
	if(INTF >0){
	  INTERFACE = INTFscore/INTF;
	}
      }     
    }
    for(int pa2=0;pa2<3;++pa2){          
      CordVec2[0]=-EVec2[pa2][0];
      CordVec2[1]=-EVec2[pa2][1];
      CordVec2[2]=-EVec2[pa2][2];           
      AlignVector(CordVec1, CordVec2, R);              
      TPharmacophores2.clear();
      TPharmacophores2=P2; 
      transformation(TPharmacophores2,R,T);
      OverlapScore(QPharmacophores, TPharmacophores2, &ShapeScore, &ChemScore, &ChargeScore, MAT, U, V, Row, Col, VEC1, VEC2, Alignment,Opt);    
      double S = Score(ShapeScore,ChemScore,V1shape,V2shape,V1chem,V2chem,Opt);
      if(S > maxoverlap){	
	bestPAxis= pa2;
	TPharmacophores=TPharmacophores2;   
	BestAlignment.clear();
	for (std::map<int,int>::iterator it=Alignment.begin(); it!=Alignment.end(); ++it){
	  BestAlignment.insert(std::make_pair(it->first,it->second));	 
	}
	maxoverlap= S;
	STanimoto = ShapeScore/(V1shape + V2shape - ShapeScore);
	CTanimoto = ChemScore/(V1chem + V2chem - ChemScore);
	TVshape   = (2*ShapeScore)/((W1*V1shape) + (W2*V2shape) + ShapeScore);
	TVchem    = (2*ChemScore)/((W1*V1chem)  + (W2*V2chem) + ChemScore);
	ChTanimoto= ChargeScore/(V1charge + V2charge - ChargeScore);	
	if(INTF >0){
	  INTERFACE = INTFscore/INTF;
	}
      }     
    }
  }
  P2=TPharmacophores;
  maxR=maxoverlap; 
  T[0]=0; T[1]=0; T[2]=0;
  for(int alphaR=0; alphaR<360; alphaR +=rotationStep){  
    TPharmacophores.clear();
    TPharmacophores=P2;     
    RotationMatrix(RM, alphaR, EVec1[bestPAxis][0], EVec1[bestPAxis][1], EVec1[bestPAxis][2]);      
    transformation(TPharmacophores,RM,T);
    for(int randomStep=0; randomStep<10; ++randomStep){	  
      TPharmacophores2.clear();
      TPharmacophores2=TPharmacophores;    
      double Tx= Delta*rand(-1,1);
      double Ty= Delta*rand(-1,1);
      double Tz= Delta*rand(-1,1);
      for(int k=0; k < npharm2; ++k){	  
	TPharmacophores2[k].point.x  += Tx;
	TPharmacophores2[k].point.y  += Ty;
	TPharmacophores2[k].point.z  += Tz;
	TPharmacophores2[k].normal.x += Tx;
	TPharmacophores2[k].normal.y += Ty;
	TPharmacophores2[k].normal.z += Tz;
      }       
            
      OverlapScore(QPharmacophores, TPharmacophores2, &ShapeScore, &ChemScore, &ChargeScore, MAT, U, V, Row, Col, VEC1, VEC2, Alignment,Opt);
      double S = Score(ShapeScore,ChemScore,V1shape,V2shape,V1chem,V2chem,Opt);
     
      double dE= (maxoverlap-S)*100;
      double r =rand(0,1); 
      if((r < exp(-dE/15)) || (S > maxR)){
	maxR=S;
	TPharmacophores=TPharmacophores2;
	if(S >maxoverlap){	    
	  maxoverlap= S;
	  BestAlignment.clear();
	  for (std::map<int,int>::iterator it=Alignment.begin(); it!=Alignment.end(); ++it){
	    BestAlignment.insert(std::make_pair(it->first,it->second));	 
	  }
	  STanimoto = ShapeScore/(V1shape + V2shape - ShapeScore);
	  CTanimoto = ChemScore/(V1chem + V2chem - ChemScore);
	  TVshape   = (2*ShapeScore)/((W1*V1shape) + (W2*V2shape) + ShapeScore);
	  TVchem    = (2*ChemScore)/((W1*V1chem)  + (W2*V2chem) + ChemScore);  
	  ChTanimoto= ChargeScore/(V1charge + V2charge - ChargeScore);
	  if(INTF >0){
	    INTERFACE = INTFscore/INTF;
	  }
	}
      } 
      
    }
  }  
  for(int i=0; i<npharm1;++i){
    int ii=i+1;
    int iii=i+2;
    if((ii >=npharm1)||(iii >=npharm1)){continue;}
    if(((P1[i].func < 2) && (P1[ii].func < 2)) &&   ((P1[i].normal.x == P1[ii].normal.x ) || (P1[i].normal.y == P1[ii].normal.y))){continue;}
    if(((P1[ii].func < 2) && (P1[iii].func < 2)) && ((P1[ii].normal.x== P1[iii].normal.x) || (P1[ii].normal.y== P1[iii].normal.y))){continue;}
    if(((P1[i].func < 2) && (P1[iii].func < 2)) &&  ((P1[i].normal.x == P1[iii].normal.x) || (P1[i].normal.y == P1[iii].normal.y))){continue;}
    VEC1[0][0]=P1[i].point.x;   VEC1[0][1]=P1[i].point.y;   VEC1[0][2]=P1[i].point.z;
    VEC1[1][0]=P1[ii].point.x;  VEC1[1][1]=P1[ii].point.y;  VEC1[1][2]=P1[ii].point.z;
    VEC1[2][0]=P1[iii].point.x; VEC1[2][1]=P1[iii].point.y; VEC1[2][2]=P1[iii].point.z;

    for(int j=0; j<npharm2; ++j){
      int jj =j+1;
      int jjj=j+2;
      if((jj >=npharm2)||(jjj >=npharm2)){continue;}
      if(((P2[j].func < 2) && (P2[jj].func < 2)) && (P2[j].normal.x == P2[jj].normal.x)){continue;}
      if(((P2[jj].func < 2) && (P2[jjj].func < 2)) && (P2[jj].normal.x == P2[jjj].normal.x)){continue;}
      if(((P2[j].func < 2) && (P2[jjj].func < 2)) && (P2[j].normal.x == P2[jjj].normal.x)){continue;}

      VEC2[0][0]=P2[j].point.x;   VEC2[0][1]=P2[j].point.y;   VEC2[0][2]=P2[j].point.z;
      VEC2[1][0]=P2[jj].point.x;  VEC2[1][1]=P2[jj].point.y;  VEC2[1][2]=P2[jj].point.z;
      VEC2[2][0]=P2[jjj].point.x; VEC2[2][1]=P2[jjj].point.y; VEC2[2][2]=P2[jjj].point.z;

      Kabsch(VEC2, VEC1, 3, 1, &RMSD, T, R);
      if(RMSD >25){continue;}
	  
      TPharmacophores2.clear();
      TPharmacophores2=P2;
      transformation(TPharmacophores2,R,T); 
      OverlapScore(QPharmacophores, TPharmacophores2, &ShapeScore, &ChemScore, &ChargeScore, MAT, U, V, Row, Col, VEC1, VEC2, Alignment,Opt);
      double S = Score(ShapeScore,ChemScore,V1shape,V2shape,V1chem,V2chem,Opt); 
      if(S >maxoverlap){
	TPharmacophores=TPharmacophores2;
	BestAlignment.clear();
	for (std::map<int,int>::iterator it=Alignment.begin(); it!=Alignment.end(); ++it){
	  BestAlignment.insert(std::make_pair(it->first,it->second));	 
	}
	maxoverlap= S;
	STanimoto = ShapeScore/(V1shape + V2shape - ShapeScore);
	CTanimoto = ChemScore/(V1chem + V2chem - ChemScore);
	TVshape   = (2*ShapeScore)/((W1*V1shape) + (W2*V2shape) + ShapeScore);
	TVchem    = (2*ChemScore)/((W1*V1chem)  + (W2*V2chem) + ChemScore);	
	ChTanimoto= ChargeScore/(V1charge + V2charge - ChargeScore);	
	if(INTF >0){
	  INTERFACE = INTFscore/INTF;
	}
      }	 
    }
  }
  
  P2=TPharmacophores;
  maxR=maxoverlap; 
  for(int randomStep=0; randomStep<500; ++randomStep){	         
    TPharmacophores2.clear();
    TPharmacophores2=TPharmacophores;     
    double Tx=Delta*rand(-1,1);
    double Ty=Delta*rand(-1,1);
    double Tz=Delta*rand(-1,1);
    for(int k=0; k < npharm2; ++k){	  
      TPharmacophores2[k].point.x  += Tx;
      TPharmacophores2[k].point.y  += Ty;
      TPharmacophores2[k].point.z  += Tz;
      TPharmacophores2[k].normal.x += Tx;
      TPharmacophores2[k].normal.y += Ty;
      TPharmacophores2[k].normal.z += Tz;
    }          
    OverlapScore(QPharmacophores, TPharmacophores2, &ShapeScore, &ChemScore, &ChargeScore, MAT, U, V, Row, Col, VEC1, VEC2, Alignment,Opt);
    double S = Score(ShapeScore,ChemScore,V1shape,V2shape,V1chem,V2chem,Opt); 
    double dE= (maxoverlap-S)*100;
    double r=rand(0,1);
    
    if((r < exp(-dE/15)) || (S > maxR)){
      maxR=S;      
      TPharmacophores=TPharmacophores2;
      if(S >maxoverlap){
	P2=TPharmacophores2;	
	BestAlignment.clear();
	for (std::map<int,int>::iterator it=Alignment.begin(); it!=Alignment.end(); ++it){
	  BestAlignment.insert(std::make_pair(it->first,it->second));	 
	}
	maxoverlap= S;
	STanimoto = ShapeScore/(V1shape + V2shape - ShapeScore);
	CTanimoto = ChemScore/(V1chem + V2chem - ChemScore);
	TVshape   = (2*ShapeScore)/((W1*V1shape) + (W2*V2shape) + ShapeScore);
	TVchem    = (2*ChemScore)/((W1*V1chem)  + (W2*V2chem) + ChemScore);	
	ChTanimoto= ChargeScore/(V1charge + V2charge - ChargeScore);
	if(INTF >0){
	  INTERFACE = INTFscore/INTF;
	}
      }
    }   
  }
  
  Overlap_Scores[0]=STanimoto;
  Overlap_Scores[1]=CTanimoto;
  Overlap_Scores[2]=TVshape;
  Overlap_Scores[3]=TVchem;
  Overlap_Scores[4]=ChTanimoto;
  Overlap_Scores[5]=INTERFACE;
  Volumes[0]=V1shape;
  Volumes[1]=V2shape;
  //	clear memory;
  DeleteArray(&VEC1, npharm_max);
  DeleteArray(&VEC2, npharm_max);
  delete[] U;
  delete[] V;
  delete[] Row;
  delete[] Col;
  DeleteArray(&MAT, npharm_max);
  DeleteArray(&DIST1, npharm_max);
  DeleteArray(&DIST2, npharm_max);
}

void OverlapScore(vector<PharmacophorePoint>& P1, vector<PharmacophorePoint>& P2, double *ShapeScore, double *ChemScore, double *ChargeScore, double **MAT, double *U, double *V, int *Row, int *Col, double **Vec1, double **Vec2, map<int,int>& Aln,int Opt){  
  
  int npharm1=P1.size();
  int npharm2=P2.size();
  int npharm_max =MAX(npharm1, npharm2); 
  double CordVec1[3]={0}, CordVec2[3]={0}; 
  double RMSD=0; double T[3]={0}; double R[3][3]={0};
  *ShapeScore=0; // Initialize the score to 0
  *ChemScore =0; // Initialize the score to 0
  *ChargeScore=0;// Initialize the score to 0

  double Shape = 0.00, Chem = 0.00, Charge=0.00, INTF=0.00;
  Aln.clear();
  FillMatrix(MAT, P1, P2, Opt);
  lap(npharm_max, MAT, Row, Col, U, V); 
  
  int k=0;
  for(int i=0; i < npharm_max; ++i){
    int j       = Row[i];    
    if(i >= npharm1){continue;}
    if(j >= npharm2){continue;}       
    if(P1[i].includedBefore){continue;}
    if(P2[j].includedBefore){continue;}
    
    CordVec1[0]=P1[i].point.x;
    CordVec1[1]=P1[i].point.y;
    CordVec1[2]=P1[i].point.z;    
    CordVec2[0]=P2[j].point.x;
    CordVec2[1]=P2[j].point.y;
    CordVec2[2]=P2[j].point.z;
    double D = distsq(CordVec1,CordVec2);   
    double D2= P1[i].vdw+P2[j].vdw;
    if(D < 6.5){
      Aln.insert(std::make_pair(i,j));
      Vec1[k][0] = P1[i].point.x;
      Vec1[k][1] = P1[i].point.y;
      Vec1[k][2] = P1[i].point.z;   
      Vec2[k][0] = P2[j].point.x;
      Vec2[k][1] = P2[j].point.y;
      Vec2[k][2] = P2[j].point.z;  
      k++;     
      
      if((P1[i].hasNormal) && (P2[j].hasNormal)){
	Vec1[k][0] = P1[i].normal.x;
        Vec1[k][1] = P1[i].normal.y;
        Vec1[k][2] = P1[i].normal.z;
        Vec2[k][0] = P2[j].normal.x;
        Vec2[k][1] = P2[j].normal.y;
        Vec2[k][2] = P2[j].normal.z;
        k++;	
	
      }       
    }
    
  }
  
  if(k >2){  
    Kabsch(Vec2, Vec1,k, 1, &RMSD, T, R);   
    for(int j=0; j < npharm2; ++j){
      CordVec2[0]=P2[j].point.x;
      CordVec2[1]=P2[j].point.y;
      CordVec2[2]=P2[j].point.z;
      P2[j].point.x= (R[0][0]*CordVec2[0]) + (R[0][1]*CordVec2[1]) + (R[0][2]*CordVec2[2]) + T[0];
      P2[j].point.y= (R[1][0]*CordVec2[0]) + (R[1][1]*CordVec2[1]) + (R[1][2]*CordVec2[2]) + T[1];
      P2[j].point.z= (R[2][0]*CordVec2[0]) + (R[2][1]*CordVec2[1]) + (R[2][2]*CordVec2[2]) + T[2];       
      CordVec2[0]=P2[j].normal.x;
      CordVec2[1]=P2[j].normal.y;
      CordVec2[2]=P2[j].normal.z;
      P2[j].normal.x= (R[0][0]*CordVec2[0]) + (R[0][1]*CordVec2[1]) + (R[0][2]*CordVec2[2]) + T[0];
      P2[j].normal.y= (R[1][0]*CordVec2[0]) + (R[1][1]*CordVec2[1]) + (R[1][2]*CordVec2[2]) + T[1];
      P2[j].normal.z= (R[2][0]*CordVec2[0]) + (R[2][1]*CordVec2[1]) + (R[2][2]*CordVec2[2]) + T[2];
    }
  }
  for(int i=0; i < npharm_max; ++i){
    int j       = Row[i];   
    if(i >= npharm1){continue;}
    if(j >= npharm2){continue;}
    if(P1[i].includedBefore){continue;}
    if(P2[j].includedBefore){continue;}  
    
    CordVec1[0] = P1[i].point.x;
    CordVec1[1] = P1[i].point.y;
    CordVec1[2] = P1[i].point.z;
    double alphaX1 = P1[i].alpha;     
    CordVec2[0] = P2[j].point.x;
    CordVec2[1] = P2[j].point.y;
    CordVec2[2] = P2[j].point.z;	    
    double alphaX2 = P2[j].alpha;
 
    double D = distsq(CordVec1, CordVec2);
    double D2= P1[i].vdw+P2[j].vdw;
    
    
    double OShape  = 1.00, OChem = 0.00, OCharge=1.00;     
    double GaussO  = (GCI2 * pow(PI/(alphaX1 + alphaX2), 1.5)) * exp(-(alphaX1 * alphaX2 * D)/(alphaX1 + alphaX2));
    if(((P1[i].func == 0) && (P2[j].func == 0)) || ((P1[i].func == 1) && (P2[j].func == 1)) ||
       ((P1[i].func == 2) && ((P2[j].func == 2) || (P2[j].func == 3) || (P2[j].func == 4))) ||
       ((P2[j].func == 2) && ((P1[i].func == 2) || (P1[i].func == 3) || (P1[i].func == 4))) ||
       ((P1[i].func == P2[j].func) && (P1[i].func < 8)) ||
       ((P1[i].alpha == P2[j].alpha) && ((P1[i].func <= 8) || (P2[j].func <= 8)))){
         OChem       = 1.00;
    }  
   
    if(P1[i].hasNormal && P2[j].hasNormal){      
      double NormalD =cosine(P1[i].normal, P2[j].normal);  
      if(P1[i].func >1)
	OChem  *= NormalD;
    }    
    
    OShape   *=  GaussO;
    OChem    *=  GaussO; 
    OCharge  *=  exp(-1*abs(P1[i].epotential - P2[j].epotential))*OChem;
    Chem   += OChem;
    Shape  += OShape;
    Charge += OCharge;
  }
  *ChemScore  =Chem;  
  *ShapeScore =Shape;
  *ChargeScore=Charge;
}


void FillMatrix(double **M, vector<PharmacophorePoint>& P1 ,vector<PharmacophorePoint>& P2, int Opt){
  // Initialize the scoring matrix
  double def=100.00;
  int npharm1=P1.size();
  int npharm2=P2.size();
  int npharm=MAX(npharm1,npharm2);
  double Vec1[3]={0}, Vec2[3]={0};
  for(int i=0;i<npharm;++i){
    for(int j=0;j <npharm; ++j){
      M[i][j]=def; // initializing all the cells to default value
    }
  }

  for(int i=0;i < npharm1;++i){
    if(P1[i].includedBefore){continue;}
    Vec1[0]=P1[i].point.x;
    Vec1[1]=P1[i].point.y;
    Vec1[2]=P1[i].point.z;    
    double alphaX1 = P1[i].alpha;   
    for(int j=0;j < npharm2;++j){
      if(P2[j].includedBefore){continue;}
      Vec2[0]=P2[j].point.x;
      Vec2[1]=P2[j].point.y;
      Vec2[2]=P2[j].point.z;
      double alphaX2 = P2[j].alpha;      
      double D = distsq(Vec1, Vec2);
      double D2= P1[i].vdw+P2[j].vdw;
      
      double OShape  = 1.00, OChem = 0.00, OCharge=0.00;
      double GaussO  = (GCI2 * pow(PI/(alphaX1 + alphaX2), 1.5)) * exp(-(alphaX1 * alphaX2 * D)/(alphaX1 + alphaX2));            
      if(((P1[i].func  == 0) && (P2[j].func == 0)) || ((P1[i].func == 1) && (P2[j].func == 1)) ||
         ((P1[i].func  == 2) && ((P2[j].func == 2) || (P2[j].func == 3) || (P2[j].func == 4))) ||
         ((P2[j].func  == 2) && ((P1[i].func == 2) || (P1[i].func == 3) || (P1[i].func == 4))) ||
	 ((P1[i].func  == P2[j].func) &&  (P1[i].func < 8)) ||
	 ((P1[i].alpha == P2[j].alpha) && ((P1[i].func <= 8) || (P2[j].func <= 8)))){
	  OChem       = 1.00;
      } 
      if(P1[i].hasNormal && P2[j].hasNormal){
        double NormalD =cosine(P1[i].normal, P2[j].normal);
        if(P1[i].func >1)
          OChem  *= NormalD;
      }

      double ESP_A=P1[i].epotential;
      double ESP_B=P2[j].epotential;
      OShape  *=  GaussO;
      OChem   *=  GaussO;        
      OCharge *=  exp(-1*abs(ESP_A - ESP_B))*GaussO;    
      if(Opt == 0){
	M[i][j] -=  OShape;
      }
      if(Opt == 1){
	M[i][j] -= (OShape + OChem);
      }
      if(Opt == 2){
	M[i][j] -=  OChem;
      }
    }
  }  
}
double
Score(double ShapeS, double ChemS, double V1shape, double V2shape, double V1chem, double V2chem, int Opt){
  double S=0;
  if(Opt == 0){
    S=(ShapeS/(V1shape + V2shape - ShapeS));
  }
  if(Opt == 1){
    S=(ShapeS/(V1shape + V2shape - ShapeS)) + (ChemS/(V1chem + V2chem - ChemS)); 
  }
  if(Opt == 2){
    S=(ChemS/(V1chem + V2chem - ChemS)); 
  }  
  return S;
}

double
cosine(Coordinate& p1, Coordinate& p2){
  double c(p1.x*p2.x + p1.y*p2.y + p1.z*p2.z);
  c /= norm(p1);
  c /= norm(p2);
  return c;
}

double
norm(Coordinate& p){
  return sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
}

double transformation(vector<PharmacophorePoint>& P, double R[3][3], double T[3]){    
  int npharm=P.size();
  double CordVec1[3]={0}, CordVec2[3]={0};
  for(int k=0; k < npharm; ++k){
    CordVec1[0]=P[k].point.x;
    CordVec1[1]=P[k].point.y;
    CordVec1[2]=P[k].point.z;
    P[k].point.x= (R[0][0]*CordVec1[0]) + (R[0][1]*CordVec1[1]) + (R[0][2]*CordVec1[2]) + T[0];
    P[k].point.y= (R[1][0]*CordVec1[0]) + (R[1][1]*CordVec1[1]) + (R[1][2]*CordVec1[2]) + T[1]; 
    P[k].point.z= (R[2][0]*CordVec1[0]) + (R[2][1]*CordVec1[1]) + (R[2][2]*CordVec1[2]) + T[2];	  
    CordVec2[0]=P[k].normal.x;
    CordVec2[1]=P[k].normal.y;
    CordVec2[2]=P[k].normal.z;
    P[k].normal.x = (R[0][0]*CordVec2[0]) + (R[0][1]*CordVec2[1]) + (R[0][2]*CordVec2[2]) + T[0];
    P[k].normal.y = (R[1][0]*CordVec2[0]) + (R[1][1]*CordVec2[1]) + (R[1][2]*CordVec2[2]) + T[1];
    P[k].normal.z = (R[2][0]*CordVec2[0]) + (R[2][1]*CordVec2[1]) + (R[2][2]*CordVec2[2]) + T[2];	  
  }
}
