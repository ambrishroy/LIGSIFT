/*######################################################################
  #   .____    .___  ________  _________._________________________     #
  #   |    |   |   |/  _____/ /   _____/|   \_   _____/\__    ___/     #
  #   |    |   |   /   \  ___ \_____  \ |   ||    __)    |    |        #
  #   |    |___|   \    \_\  \/        \|   ||     \     |    |        #
  #   |_______ \___|\______  /_______  /|___|\___  /     |____|        #
  #           \/           \/        \/          \/                    #
  #                                                                    #
  ######################################################################
*/
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include "Pharmacophores.h"
#include "Smallmolecule.h"
#include "Alignment.h"
#include "RMSD.h"
#include "Common.h"


using namespace std;
using namespace OpenBabel;

int main(int argc, char *argv[]){
  string version ("LIGSIFT (v1.3)");

  cout << " ********************************************************************************************* " << endl
       << " * "<<version <<": An open-source tool for ligand structural alignment and virtual screening * " << endl
       << " * Ambrish Roy & Jeffrey Skolnick                                                            * " << endl
       << " * Please mail your comments to Ambrish Roy (ambrish.roy@gmail.com)                          * " << endl
       << " ********************************************************************************************* " << endl;
  
  std::string qFilename, dFilename, oFilename, sFilename, OptimInp;
  
  bool vOpt= false;
  bool qOpt= false; bool dbOpt= false; bool oOpt= false;
  int Optim=2; int sOpt=0;
  if(argc < 7){
    cout << endl
         << "Usage:" << endl
         << "LIGSIFT -q    <query molecules, MOL2>" << endl
         << "        -db   <database molecules,MOL2>" << endl 
         << "        -o    <Output score file>"<< endl
         << "Additional options:"<< endl
         << "        -s    <Superposed database molecules on query,MOL2>" << endl
         << "        -v    For verbose output" << endl
         << "        -opt  0(Overlap optimized based on Shape similarity)" << endl
         << "              1(Overlap optimized based on Shape + Chemical similarity)" << endl
         << "              2(Overlap optimized based on Chemical similarity:Default method)" << endl
         << "For example:" << endl << endl
         << "./LIGSIFT -q ZINCXXX.mol2   -db DB.mol2   -o alignment_scores.txt   -s superposed.mol2" << endl << endl;
    return 1;
  }
  
  for(int i = 0; i < argc; i++ ){
    if(!strcmp(argv[i],"-q")  && i < argc ) { qFilename      = string(argv[i+1]); qOpt = true;}
    if(!strcmp(argv[i],"-Q")  && i < argc ) { qFilename      = string(argv[i+1]); qOpt = true;}
    if(!strcmp(argv[i],"-db") && i < argc ) { dFilename      = string(argv[i+1]); dbOpt= true;}
    if(!strcmp(argv[i],"-DB") && i < argc ) { dFilename      = string(argv[i+1]); dbOpt= true;}
    if(!strcmp(argv[i],"-o")  && i < argc ) { oFilename      = string(argv[i+1]); oOpt = true;}
    if(!strcmp(argv[i],"-O")  && i < argc ) { oFilename      = string(argv[i+1]); oOpt = true;}
    if(!strcmp(argv[i],"-s")  && i < argc ) { sFilename      = string(argv[i+1]); sOpt = 1;}
    if(!strcmp(argv[i],"-opt")&& i < argc ) { OptimInp       = string(argv[i+1]); Optim= atoi(OptimInp.c_str());}
    if(!strcmp(argv[i],"-v")  && i < argc ) { vOpt           = true;}
    if(!strcmp(argv[i],"-V")  && i < argc ) { vOpt           = true;}
  }
  
  if (!qOpt){
    cout << "Please provide query molecules (preferably in MOL2 format)" << endl;
    exit(EXIT_FAILURE);
  }
  if(!dbOpt){
    cout << "Please provide database molecules (preferably in MOL2 format)" << endl;
    exit(EXIT_FAILURE);
  }
  if(!oOpt){
    cout << "Please provide a filename for recording alignment scores." << endl;
    exit(EXIT_FAILURE);
  }
  if(Optim >2){Optim=2;}

  
  ifstream qFile(qFilename.c_str());
  ifstream dFile(dFilename.c_str());
  if( !qFile.good() ){
    cout << "ERROR: " << endl
	 << "Cannot open query file " << qFilename << endl;
    exit(EXIT_FAILURE);
  }	   
  if( !dFile.good() ){
    cout << "ERROR: " << endl
	 << "Cannot open database file " << dFilename << endl;
    exit(EXIT_FAILURE);
  }
  
  OBMol OBmol1, OBmol2, Dmol;
  SmallMolecule SM1, SM2;
  Pharmacophore QPharm, TPharm;
  OBFormat     *inFormat1, *inFormat2,*outFormat;
  OBConversion  OBconv1, OBconv2;
  vector<std::string> Qnames, Tnames;
  vector<PharmacophorePoint> PH;
  string ff1  = "MMFF94";
  string ff2  = "GAFF";
  
  
  string line;
  int Nmol1=0, Nmol2=0;
  OBForceField* pFF1    = OBForceField::FindForceField(ff1);
  OBForceField* pFF2    = OBForceField::FindForceField(ff1); 
  pFF1->SetLogFile(&cerr);
  pFF1->SetLogLevel(OBFF_LOGLVL_NONE);
  pFF2->SetLogFile(&cerr);
  pFF2->SetLogLevel(OBFF_LOGLVL_NONE);

  
  vector<PharmacophorePoint> QPharmacophores;
  vector<PharmacophorePoint> TPharmacophores;
  vector<OBMol> DBmol;
  ifstream ifs1(qFilename.c_str());
  ofstream ofs1;  
  inFormat1 = OBconv1.FormatFromExt(qFilename.c_str());
  outFormat = OBconv2.FindFormat("mol2");  
  OBconv1.SetInAndOutFormats(inFormat1,outFormat);
  inFormat2 = OBconv2.FormatFromExt(dFilename.c_str());  
  OBconv2.SetInAndOutFormats(inFormat2,outFormat);
  //##############Start with query Ligand file###############//
  bool notatend = OBconv1.ReadFile(&OBmol1,qFilename);
  while (notatend){    
    (bool) OBmol1.AddPolarHydrogens();
    string MolID = OBmol1.GetTitle();
    if(!pFF1->Setup(OBmol1)){     
      pFF1 = OBForceField::FindForceField(ff2);
    }
    if(!pFF1->Setup(OBmol1)){
      cout << "ERROR: Force-field set up failed for " << MolID << endl;
    }
    else{
      pFF1->GetPartialCharges(OBmol1);        
      SM1.getMoleculeData(OBmol1);   
      SM1.InitPharmacophore(OBmol1,PH);       
      QPharm.push_back(PH);
      Qnames.push_back(MolID);
      SM1.clear_molecule(); PH.clear();
      ++Nmol1;
    }
    OBmol1.Clear();
    notatend = OBconv1.Read(&OBmol1);   
  }
   
  //##############Now read database molecules###############//
  bool notatend2 = OBconv2.ReadFile(&OBmol2,dFilename);
  while (notatend2){    
    (bool) OBmol2.AddPolarHydrogens();
    string MolID = OBmol2.GetTitle();
    if(!pFF2->Setup(OBmol2)){
      pFF2 = OBForceField::FindForceField(ff2);
    }
    if(!pFF2->Setup(OBmol2)){
      cout << "ERROR: Force-field set up failed for " << MolID << endl;
    }
    else{
      pFF2->GetPartialCharges(OBmol2);
      SM2.getMoleculeData(OBmol2);
      SM2.InitPharmacophore(OBmol2,PH);     
      TPharm.push_back(PH);
      Tnames.push_back(MolID);
      
      SM2.clear_molecule(); PH.clear();   
      if(sOpt){
	DBmol.push_back(OBmol2);	
      }
      ++Nmol2; 
    }
    OBmol2.Clear();
    notatend2 = OBconv2.Read(&OBmol2);   
  }
  
  if(vOpt){
    cout << "No. of molecules read from query file : "<< QPharm.size() << endl;
    cout << "No. of molecules read from database   : "<< TPharm.size() << endl;  
    if(Optim==2){
      cout << "Using -opt 2 (default) Chemical similarity for finding the best overlap." << endl;
    }
    if(Optim==0){
      cout << "Using -opt 0 Shape similarity for finding the best overlap." << endl;
    }
    if(Optim==1){
      cout << "Using -opt 1 combination of Shape & Chemical similarity (1:1) for finding the best overlap." << endl;
    }
  }
 
  

  
  FILE * oFILE;
  int ic=0;
  std::map<int,int> Alignment;
  double RMSD; 
  double T[3]={0}; double R[3][3]={0}; double CordVec[3];
  std::vector<double> ScoreArr;
  std::vector<double> Volumes;
  ScoreArr.resize(6);
  Volumes.resize(2);
  oFILE= fopen(oFilename.c_str(),"w");
  if(vOpt){
     printf("%-15s\t%-15s\t%-13s  %-13s   %-10s %-8s  %-13s %-13s %-13s %-13s\n","Database_name","Query_name","ShapeTanimoto","ChemTanimoto","ShapeSim","ChemSim","ShapeSimPval","ChemSimPval","TverskyShape","TverskyChem");    
  }
  fprintf (oFILE,"%-15s\t%-15s\t%-13s  %-13s   %-10s %-8s  %-13s %-13s %-13s %-13s\n","Database_name","Query_name","ShapeTanimoto","ChemTanimoto","ShapeSim","ChemSim","ShapeSimPval","ChemSimPval","TverskyShape","TverskyChem");
  
  double **VECQ, **VECT;
  for(unsigned int i=0;i<QPharm.size();++i){
    QPharmacophores=QPharm[i];
    for(unsigned int j=0;j<TPharm.size();++j){
      TPharmacophores=TPharm[j];	
      ScoreArr[0]=0;ScoreArr[1]=0;ScoreArr[2]=0; ScoreArr[3]=0;ScoreArr[4]=0;ScoreArr[5]=0;
      Volumes[0] =0; Volumes[1]=0;      	
      Overlap(QPharmacophores,TPharmacophores,ScoreArr, Volumes, Alignment, Optim);
      double v1  = Volumes[0];
      double v2  = Volumes[1];
      if((v1 ==0) || (v2==0))
	continue;
      double TCshape=0.0, TCchem=0.0, sTCshape=0.0, sTCchem=0.0, pval_shape=1.00, pval_chem=1.0;
      TCshape=ScoreArr[0];
      TCchem =ScoreArr[1];
      PVal(TCshape, v1, v2, 0, Optim, &sTCshape, &pval_shape);
      PVal(TCchem,  v1, v2, 1, Optim, &sTCchem,  &pval_chem);
      if(vOpt){
	printf("%-15s\t%-15s\t   %-13.3f  %-13.3f %-10.3f %-8.3f %e  %e    %-13.3f %-13.3f\n",Tnames[j].c_str(),Qnames[i].c_str(),TCshape, TCchem, sTCshape, sTCchem, pval_shape, pval_chem, ScoreArr[2], ScoreArr[3]); 
      }
      fprintf (oFILE,"%-15s\t%-15s\t   %-13.3f  %-13.3f %-10.3f %-8.3f %e  %e    %-13.3f %-13.3f\n",Tnames[j].c_str(),Qnames[i].c_str(),TCshape, TCchem, sTCshape, sTCchem, pval_shape, pval_chem, ScoreArr[2], ScoreArr[3]);
     
      if(sOpt==1){	
	QPharmacophores.clear();TPharmacophores.clear();
	int pno=0;
	QPharmacophores=QPharm[i];
	TPharmacophores=TPharm[j];	
	int L1= QPharm[i].size(); int L2=TPharm[j].size();	
	if((L1>0) && (L2>0)){	  
	  CreateArray(&VECQ,L1, 3); CreateArray(&VECT,L2, 3);	  
	  for (std::map<int,int>::iterator it=Alignment.begin(); it!=Alignment.end(); ++it){	  
	    int iQ=it->first;
	    int iT=it->second;
	    VECQ[pno][0]=QPharmacophores[iQ].point.x; VECQ[pno][1]=QPharmacophores[iQ].point.y; VECQ[pno][2]=QPharmacophores[iQ].point.z; 
	    VECT[pno][0]=TPharmacophores[iT].point.x; VECT[pno][1]=TPharmacophores[iT].point.y; VECT[pno][2]=TPharmacophores[iT].point.z;
	    ++pno;	  
	  }
	  Kabsch(VECT, VECQ, 3, 1, &RMSD, T, R);
	  if(vOpt){
	      printf("Transformation matrix to rotate Database molecule on query -----\n");
              printf("i\t%18s %15s %15s %15s\n", "T[i]", "R[i][0]", "R[i][1]", "R[i][2]");
              for (unsigned int k = 0; k < 3; k++){
                 printf("%d\t%18.10f %15.10f %15.10f %15.10f\n", k, T[k], R[k][0], R[k][1], R[k][2]);
              }
          }
	  Dmol=DBmol[j];
	  FOR_ATOMS_OF_MOL(atom,Dmol){   
	    CordVec[0]=atom->x();
	    CordVec[1]=atom->y();
	    CordVec[2]=atom->z();    
	    double nX=transform(T,R,CordVec,0);
	    double nY=transform(T,R,CordVec,1);    
	    double nZ=transform(T,R,CordVec,2);	  	   
	    atom->SetVector(nX,nY,nZ);	  
	  }      
	  ofs1.open(sFilename.c_str(),ios::app);
	  OBconv2.Write(&Dmol,&ofs1);
	  ofs1.close();	  
	  DeleteArray(&VECQ,L1);
	  DeleteArray(&VECT,L2);	  
	}
      }
      Alignment.clear();
      ++ic;        
    }
  }
  fclose(oFILE);

  return 0;
}


