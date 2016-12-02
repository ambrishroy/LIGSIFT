#include "Smallmolecule.h"
using namespace std;
using namespace OpenBabel;

	
void
SmallMolecule::getMoleculeData(OpenBabel::OBMol& OBmol){	
  int natoms=0, nbonds=0, atomno=0, nrings;
  natoms=OBmol.NumAtoms();
  nbonds=OBmol.NumBonds(); 
  std::vector<OpenBabel::OBRing*> rings = OBmol.GetSSSR();
  nrings  =  rings.size();
  allocate_arrays(natoms, nbonds, nrings);
  //##Center molecule around its COG/
  //OBmol.Center();
  
  unsigned int n, a(0);
  //#### Radii & Rings in molecule/
  FOR_ATOMS_OF_MOL(atom, OBmol){
    a = atom->GetIdx();
    n = (unsigned int) atom->GetAtomicNum();
    switch (n){
    case 1:  // H
      radii[a] = 1.20;
      break;
    case 3:  // Li
      radii[a] = 1.82;
      break;
    case 5:  // B
      radii[a] = 2.00;
      break;
    case 6:  // C
      radii[a] = 1.70;
      break;
    case 7:  // N
      radii[a] = 1.55;
      break;
    case 8:  // O
      radii[a] = 1.52;
      break;
    case 9:  // F
      radii[a] = 1.47;
      break;
    case 11: // Na
      radii[a] = 2.27;
      break;
    case 12: // Mg
      radii[a] = 1.73;
      break;
    case 14: // Si
      radii[a] = 2.10;
      break;
    case 15: // P
      radii[a] = 1.80;
      break;
    case 16: // S
      radii[a] = 1.80;
      break;
    case 17: // Cl
      radii[a] = 1.75;
      break;
    case 19: // K
      radii[a] = 2.75;
      break;
    case 20: // Ca
      radii[a] = 2.00;
      break;
    case 26: // Fe
      radii[a] = 1.10;
      break;
    case 29: // Cu
      radii[a] = 1.40;
      break;
    case 30: // Zn
      radii[a] = 1.39;
      break;
    case 35: // Br
      radii[a] = 1.85;
      break;
    case 53: // I
      radii[a] = 1.98;
      break;
    default:
      radii[a]= 1.50;
    }
    //#### Counter/
    if(atom->IsInRing()){
      flag_ring[atom->GetIdx()]=true;
      if(atom->IsAromatic()){
	flag_ring_aromatic[atom->GetIdx()]=true;
      }
      else{
	flag_ring_nonaromatic[atom->GetIdx()]=true;
      }
    }
  }	
  for (unsigned int i(0); i < rings.size(); ++i){
    int rsize = rings[i]->Size();
    double var = 180/rsize;
    double radii = 1.4/(2*sin(var*PI/180));
    ringradii[i] = radii;		
  }
  
  OBSmartsPattern  OH, HBA1, HBA2, HBD1, HBD2;
  OBSmartsPattern  Pho1, Pho2, Pho3, Pho4, Pho5, Pho6, Pho7, Pho8, Pho9, Pho10, Pho11, Pho12, Pho13;
  /*##SMARTS pattern to search in molecule*/   
  OH.Init("[OX2H]");// Hydroxyl can act as both H-bond donor and acceptor
  HBA1.Init("[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]");
  HBA2.Init("[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]");
  HBD1.Init("[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]");    
  HBD2.Init("[#8!H0&!$([OH][C,S,P]=O)]");
  //##branched terminals as one point/
  Pho1.Init("[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]");
  Pho2.Init("[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]");
  Pho3.Init("*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]");
  //##simple rings only; need to combine points to get good results for 3d structures/    
  Pho4.Init("[C&r3]1~[C&r3]~[C&r3]1");
  Pho5.Init("[C&r4]1~[C&r4]~[C&r4]~[C&r4]1");
  Pho6.Init("[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1");
  Pho7.Init("[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1");
  Pho8.Init("[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1");
  Pho9.Init("[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1");
  //##aliphatic chains/
  Pho10.Init("[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]");
  Pho11.Init("[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]");
  Pho12.Init("[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]");
  //##sulfur (apparently)/
  Pho13.Init("[$([S]~[#6])&!$(S~[!#6])]");

 
  
  vector< vector <int> > maplist;
  vector< vector <int> >::iterator matches;
  
  if(OH.Match(OBmol)){
    maplist    = OH.GetUMapList();
    for(unsigned int m(0); m < maplist.size(); ++m){
      for(unsigned int a(0);a<maplist[m].size();++a){
	OBAtom* atom = OBmol.GetAtom(maplist[m][a]);
	atomno=atom->GetIdx();				
	flag_hydroxyl[atomno]=true;
      }
    }		
  }

  std::vector<OpenBabel::OBAtom*>::iterator ai;
  for(OpenBabel::OBAtom* a = OBmol.BeginAtom(ai); a; a = OBmol.NextAtom(ai)){
    atomno=a->GetIdx();
    if (a->GetAtomicNum() == 7 || a->GetAtomicNum() == 8){
      if (a->GetFormalCharge() >= 0 && ((a->GetImplicitValence() - a->GetHvyValence()) !=0)){
            flag_donor[atomno]=true;
      }
      if(a->GetFormalCharge() <= 0){
	flag_acceptor[atomno]=true;
      }
    }
    std::string S(a->GetType());
    if((a->GetAtomicNum() == 0) && (strncmp (S.c_str(),"N",1)==0) && (a->GetValence()!=3)){
       flag_donor[atomno]=true;
    }
  }    
  if(HBA1.Match(OBmol)){
    maplist    = HBA1.GetUMapList();
    for(unsigned int m(0); m < maplist.size(); ++m){
      for(unsigned int a(0);a<maplist[m].size();++a){
	OBAtom* atom = OBmol.GetAtom(maplist[m][a]);
	atomno=atom->GetIdx();
	if( flag_hydroxyl[atomno]){
	  continue;
	}
	else{
	  flag_acceptor[atomno]=true;
	}
      }
    }
  }
  if(HBA2.Match(OBmol)){
    maplist    = HBA2.GetUMapList();
    for(unsigned int m(0); m < maplist.size(); ++m){
      for(unsigned int a(0);a<maplist[m].size();++a){
	OBAtom* atom = OBmol.GetAtom(maplist[m][a]);
	atomno=atom->GetIdx();
	if( flag_hydroxyl[atomno]){
	  continue;
	}
	else{	  
	  flag_acceptor[atomno]=true;
	}
      }
    }
  }  
  if(HBD1.Match(OBmol)){
    maplist    = HBD1.GetUMapList();
    for(unsigned int m(0); m < maplist.size(); ++m){
      for(unsigned int a(0);a<maplist[m].size();++a){
	OBAtom* atom = OBmol.GetAtom(maplist[m][a]);
	atomno=atom->GetIdx();
	if( flag_hydroxyl[atomno]){
	  continue;
	}
	else{
	  flag_donor[atomno]=true;
	}
      }
    }
  }  
  if(HBD2.Match(OBmol)){
    maplist    = HBD2.GetUMapList();
    for(unsigned int m(0); m < maplist.size(); ++m){
      for(unsigned int a(0);a<maplist[m].size();++a){
	OBAtom* atom = OBmol.GetAtom(maplist[m][a]);
	atomno=atom->GetIdx();
	if( flag_hydroxyl[atomno]){
	  continue;
	}
	else{
	  flag_donor[atomno]=true;
	}	
      }
    }
  }
  if(Pho1.Match(OBmol)){		
    maplist = Pho1.GetUMapList();
    std::vector<int> groups;
    groups.clear();
    for(unsigned int m(0); m < maplist.size(); ++m){
      for(unsigned int a(0);a<maplist[m].size();++a){
	OBAtom* atom = OBmol.GetAtom(maplist[m][a]);
	atomno=atom->GetIdx();
	groups.push_back(atomno);
      }
    }
    hydrophobic_groups.push_back(groups);
  }
  if(Pho2.Match(OBmol)){          
    maplist = Pho2.GetUMapList(); 
    std::vector<int> groups; 
    groups.clear(); 
    for(unsigned int m(0); m < maplist.size(); ++m){ 
      for(unsigned int a(0);a<maplist[m].size();++a){ 
        OBAtom* atom = OBmol.GetAtom(maplist[m][a]); 
        atomno=atom->GetIdx(); 
        groups.push_back(atomno); 
      } 
    } 
    hydrophobic_groups.push_back(groups);
  }
  if(Pho3.Match(OBmol)){           
    maplist = Pho3.GetUMapList();  
    std::vector<int> groups;  
    groups.clear();  
    for(unsigned int m(0); m < maplist.size(); ++m){  
      for(unsigned int a(0);a<maplist[m].size();++a){  
        OBAtom* atom = OBmol.GetAtom(maplist[m][a]);  
        atomno=atom->GetIdx();  
        groups.push_back(atomno);  
      }  
    }  
    hydrophobic_groups.push_back(groups); 
  }
  if(Pho4.Match(OBmol)){           
    maplist = Pho4.GetUMapList();  
    std::vector<int> groups;  
    groups.clear();  
    for(unsigned int m(0); m < maplist.size(); ++m){  
      for(unsigned int a(0);a<maplist[m].size();++a){  
        OBAtom* atom = OBmol.GetAtom(maplist[m][a]);  
        atomno=atom->GetIdx();  
        groups.push_back(atomno);  
      }  
    }  
    hydrophobic_groups.push_back(groups); 
  } 
  if(Pho5.Match(OBmol)){           
    maplist = Pho5.GetUMapList();  
    std::vector<int> groups;  
    groups.clear();  
    for(unsigned int m(0); m < maplist.size(); ++m){  
      for(unsigned int a(0);a<maplist[m].size();++a){  
        OBAtom* atom = OBmol.GetAtom(maplist[m][a]);  
        atomno=atom->GetIdx();  
        groups.push_back(atomno);  
      }  
    }  
    hydrophobic_groups.push_back(groups); 
  } 
  if(Pho6.Match(OBmol)){           
    maplist = Pho6.GetUMapList();  
    std::vector<int> groups;  
    groups.clear();  
    for(unsigned int m(0); m < maplist.size(); ++m){  
      for(unsigned int a(0);a<maplist[m].size();++a){  
        OBAtom* atom = OBmol.GetAtom(maplist[m][a]);  
        atomno=atom->GetIdx();  
        groups.push_back(atomno);  
      }  
    }  
    hydrophobic_groups.push_back(groups); 
  } 
  if(Pho7.Match(OBmol)){           
    maplist = Pho7.GetUMapList();  
    std::vector<int> groups;  
    groups.clear();  
    for(unsigned int m(0); m < maplist.size(); ++m){  
      for(unsigned int a(0);a<maplist[m].size();++a){  
        OBAtom* atom = OBmol.GetAtom(maplist[m][a]);  
        atomno=atom->GetIdx();  
        groups.push_back(atomno);  
      }  
    }  
    hydrophobic_groups.push_back(groups); 
  } 
  if(Pho8.Match(OBmol)){           
    maplist = Pho8.GetUMapList();  
    std::vector<int> groups;  
    groups.clear();  
    for(unsigned int m(0); m < maplist.size(); ++m){  
      for(unsigned int a(0);a<maplist[m].size();++a){  
        OBAtom* atom = OBmol.GetAtom(maplist[m][a]);  
        atomno=atom->GetIdx();  
        groups.push_back(atomno);  
      }  
    }  
    hydrophobic_groups.push_back(groups); 
  } 
  if(Pho9.Match(OBmol)){           
    maplist = Pho9.GetUMapList();  
    std::vector<int> groups;  
    groups.clear();  
    for(unsigned int m(0); m < maplist.size(); ++m){  
      for(unsigned int a(0);a<maplist[m].size();++a){  
        OBAtom* atom = OBmol.GetAtom(maplist[m][a]);  
        atomno=atom->GetIdx();  
        groups.push_back(atomno);  
      }  
    }  
    hydrophobic_groups.push_back(groups); 
  } 
  if(Pho10.Match(OBmol)){           
    maplist = Pho10.GetUMapList();  
    std::vector<int> groups;  
    groups.clear();  
    for(unsigned int m(0); m < maplist.size(); ++m){  
      for(unsigned int a(0);a<maplist[m].size();++a){  
        OBAtom* atom = OBmol.GetAtom(maplist[m][a]);  
        atomno=atom->GetIdx();  
        groups.push_back(atomno);  
      }  
    }  
    hydrophobic_groups.push_back(groups); 
  } 
  if(Pho11.Match(OBmol)){           
    maplist = Pho11.GetUMapList();  
    std::vector<int> groups;  
    groups.clear();  
    for(unsigned int m(0); m < maplist.size(); ++m){  
      for(unsigned int a(0);a<maplist[m].size();++a){  
        OBAtom* atom = OBmol.GetAtom(maplist[m][a]);  
        atomno=atom->GetIdx();  
        groups.push_back(atomno);  
      }  
    }  
    hydrophobic_groups.push_back(groups); 
  } 
  if(Pho12.Match(OBmol)){            
    maplist = Pho12.GetUMapList();   
    std::vector<int> groups;   
    groups.clear();   
    for(unsigned int m(0); m < maplist.size(); ++m){   
      for(unsigned int a(0);a<maplist[m].size();++a){   
        OBAtom* atom = OBmol.GetAtom(maplist[m][a]);   
        atomno=atom->GetIdx();   
        groups.push_back(atomno);   
      }   
    }   
    hydrophobic_groups.push_back(groups);  
  }  
  if(Pho13.Match(OBmol)){             
    maplist = Pho13.GetUMapList();    
    std::vector<int> groups;    
    groups.clear();    
    for(unsigned int m(0); m < maplist.size(); ++m){    
      for(unsigned int a(0);a<maplist[m].size();++a){    
        OBAtom* atom = OBmol.GetAtom(maplist[m][a]);    
        atomno=atom->GetIdx();    
        groups.push_back(atomno);    
      }    
    }    
    hydrophobic_groups.push_back(groups);   
  }
}
/*************************************/
void
SmallMolecule::InitPharmacophore(OpenBabel::OBMol& OBmol, vector<PharmacophorePoint>& Phar){
  Phar.clear();

  int atomno=0;
  std::vector<OpenBabel::OBAtom*>::iterator ai;
  for (OpenBabel::OBAtom* a = OBmol.BeginAtom(ai); a; a = OBmol.NextAtom(ai)){
    //    OBPairData *chg = (OBPairData*) a->GetData("FFPartialCharge");
    int charge = a->GetFormalCharge();
    atomno     = a->GetIdx();   
    double vd  = radii[atomno];
    double pc  = a->GetPartialCharge();
    /*## 1. Anion*/
    if (charge < 0){
      PharmacophorePoint p;
      p.func = NEGC;
      p.point.x = a->x();
      p.point.y = a->y();
      p.point.z = a->z();
      p.alpha = GCI3/pow(radii[atomno],2);
      p.hasNormal = false;
      p.pcharge    = pc;
      p.vdw        = vd;
      if(flag_in_pharmacophore[atomno]){p.includedBefore = true;}
      Phar.push_back(p);
      flag_in_pharmacophore[atomno]=true;
      //cout << atomno << " NEGC  " << p.func << endl;
    }
    /*## 2. Cation*/
    else if (charge > 0){
      PharmacophorePoint p;
      p.func = POSC;
      p.point.x = a->x();
      p.point.y = a->y();
      p.point.z = a->z();
      p.alpha = GCI3/pow(radii[atomno],2);
      p.hasNormal = false;
      p.pcharge    = pc;
      p.vdw        = vd;
      if(flag_in_pharmacophore[atomno]){p.includedBefore = true;}
      Phar.push_back(p);
      flag_in_pharmacophore[atomno]=true;
      //cout << atomno << " POSC  " << p.func << endl;
    }
    /*## 3. OH Group*/
    if(flag_hydroxyl[atomno]){
      PharmacophorePoint p;
      p.func = OHG;
      p.point.x = a->x();
      p.point.y = a->y();
      p.point.z = a->z();
      p.hasNormal = true;
      p.pcharge    = pc;
      p.vdw        = vd;
      if(flag_in_pharmacophore[atomno]){p.includedBefore = true;}
      p.alpha = GCI3/pow(radii[atomno],2);
      p.normal = CalcNormal(a);
      Phar.push_back(p);
      flag_in_pharmacophore[atomno]=true;
      //cout << atomno << " OHG  " << p.func << endl;
    }
    /*## 4. HBond Donor */ 	
    else if((a->IsHbondDonor()) || (flag_donor[atomno])){
      PharmacophorePoint p;
      p.func = HDON;
      p.point.x = a->x();
      p.point.y = a->y();
      p.point.z = a->z();
      p.hasNormal = true;
      if(flag_in_pharmacophore[atomno]){p.includedBefore = true;}
      p.alpha = GCI3/pow(radii[atomno],2);
      p.normal = CalcNormal(a);
      p.pcharge    = pc;
      p.vdw        = vd;
      Phar.push_back(p);
      flag_in_pharmacophore[atomno]=true;
      //cout << atomno << " HDON  " << p.func << endl;
    }    
    /*## 5. HBond Acceptor */
    else if((a->IsHbondAcceptor()) || (flag_acceptor[atomno])){
      PharmacophorePoint p;
      p.func = HACC;
      p.point.x = a->x();
      p.point.y = a->y();
      p.point.z = a->z();
      p.hasNormal = true;
      if(flag_in_pharmacophore[atomno]){p.includedBefore = true;}
      p.alpha  = GCI3/pow(radii[atomno],2);
      p.normal = CalcNormal(a);
      p.pcharge = pc;
      p.vdw     = vd;
      Phar.push_back(p);
      flag_in_pharmacophore[atomno]=true;   
      //cout << atomno << " HACC  " << p.func << endl;
    }
  }
  /*## 6. Rings atoms */
  OpenBabel::vector3 center;
  OpenBabel::vector3 norm1;
  OpenBabel::vector3 norm2;
  std::vector<OpenBabel::OBRing*> rings = OBmol.GetSSSR();
  vector<OBRing*>::iterator i;
  vector<int>::iterator j; 
  for (unsigned int i(0); i < rings.size(); ++i){
   
    (bool) rings[i]->findCenterAndNormal(center, norm1, norm2);
    for(j = rings[i]->_path.begin(); j != rings[i]->_path.end(); ++j){      
      int atomno     =*j;
      OBAtom *atom = OBmol.GetAtom(atomno);
      PharmacophorePoint p;
      p.func = NAROM;     
      double pc  = atom->GetPartialCharge();
      double vd  = radii[atomno];
      
      if(rings[i]->IsAromatic()){
	p.func = AROM;
      }       
      p.point.x = atom->x();
      p.point.y = atom->y();
      p.point.z = atom->z();
      p.hasNormal = true;
      p.pcharge   = pc; 
      p.vdw       = vd;
      if(flag_in_pharmacophore[atomno]){continue;}       
      p.normal.x = center.x();
      p.normal.y = center.y();
      p.normal.z = center.z();
      p.alpha    = GCI3/pow(radii[atomno],2);
      p.pcharge    = pc;

      Phar.push_back(p);
      flag_in_pharmacophore[atomno]=true;
    }
  }
  /*## 7. Hydrophobic groups atoms*/
  for(unsigned int h1=0;h1< hydrophobic_groups.size();++h1){
    PharmacophorePoint p;
    p.func = HPHOB;            
    p.hasNormal = false;                                
    vector<int> group_elements=hydrophobic_groups[h1];
    for(unsigned int h2=0;h2 < group_elements.size();++h2){			
      OBAtom *atom = OBmol.GetAtom(group_elements[h2]);      
      if(flag_in_pharmacophore[group_elements[h2]]){continue;}      
      int atomno     =group_elements[h2]; 
      double pc  = atom->GetPartialCharge();
      double vd  = radii[atomno];
      PharmacophorePoint p;
      p.func = HPHOB;
      p.hasNormal = false;      
      p.point.x   = atom->x();
      p.point.y   = atom->y();
      p.point.z   = atom->z();
      p.alpha     = GCI3/pow(radii[atomno],2);
      p.pcharge   = pc;
      p.vdw       = vd;
      p.hasNormal = false;
      Phar.push_back(p);
      flag_in_pharmacophore[group_elements[h2]]=true;      
    }
  }
  /*## 8. Mark remaining non-hydrogen atoms as excluded volume*/
  for (OpenBabel::OBAtom* a = OBmol.BeginAtom(ai); a; a = OBmol.NextAtom(ai)){
    atomno = a->GetIdx();
    if((!flag_in_pharmacophore[atomno]) &&(a->GetAtomicNum() > 1)){
      double pc  = a->GetPartialCharge();
      double vd  = radii[atomno];      
      PharmacophorePoint p;
      p.func = EXCL;
      p.point.x = a->x();
      p.point.y = a->y();
      p.point.z = a->z();
      p.alpha  = GCI3/pow(radii[atomno],2);
      p.pcharge = pc;
      p.vdw     = vd;
      p.hasNormal = false;
      Phar.push_back(p);
      flag_in_pharmacophore[atomno]=true; 
    }    
  }
  /*## 9. Calculate electrostatic potential for each atom*/
  for(unsigned int i=0; i <Phar.size(); ++i){        
    for(unsigned int j=0; j <Phar.size(); ++j){             
      if(Phar[j].includedBefore){continue;}
      double Xd=Phar[i].point.x-Phar[j].point.x;
      double Yd=Phar[i].point.y-Phar[j].point.y;
      double Zd=Phar[i].point.z-Phar[j].point.z;
      double D =sqrt((Xd*Xd) + (Yd*Yd) + (Zd*Zd));
      //if(D<=Phar[i].vdw)D=Phar[i].vdw;
      if(D<Phar[i].vdw)continue;
      Phar[i].epotential += (Phar[j].pcharge)/D;            
    }
    Phar[i].epotential *=EC;       
  }  
}
/*************************************/
Coordinate
SmallMolecule::CalcNormal(OpenBabel::OBAtom* a){
  int nbrBonds(0);
  Coordinate normal;
  
  std::vector<OpenBabel::OBBond*>::iterator bi;
  for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi)){
    OpenBabel::OBAtom* aa = b->GetNbrAtom(a);
    if (aa->GetAtomicNum() == 1){
      continue;
    }
    ++nbrBonds;
    normal.x += (aa->x() - a->x());
    normal.y += (aa->y() - a->y());
    normal.z += (aa->z() - a->z());
  }
  double length(sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z));
  if(length >0){
      normal.x /= length;
      normal.y /= length;
      normal.z /= length;

      normal.x = -normal.x;
      normal.y = -normal.y;
      normal.z = -normal.z;

      normal.x += a->x();
      normal.y += a->y();
      normal.z += a->z();
  }
  else{
      normal.x = a->x();
      normal.y = a->y();
      normal.z = a->z();
  }
  return normal;
}
/*************************************/
void
SmallMolecule::allocate_arrays(int natoms, int nbonds, int nrings){
  arrays_allocated = false;	
  clear_molecule();

  /* Keep the same mapping as ligand atom no. i.e Starting from 1*/
  flag_hydroxyl    = new bool[natoms+1];
  flag_cation      = new bool[natoms+1];
  flag_anion       = new bool[natoms+1];
  flag_acceptor    = new bool[natoms+1];
  flag_donor       = new bool[natoms+1];    
  flag_ring        = new bool[natoms+1];
  flag_ring_aromatic    =new bool[natoms+1];
  flag_ring_nonaromatic =new bool[natoms+1];
  flag_in_pharmacophore =new bool[natoms+1];
  radii                 =new double[natoms+1];
  ringradii		=new double[nrings+1]; 
  for (int i=1;i<=natoms;i++){
    flag_hydroxyl[i]  = false;
    flag_cation[i]    = false;
    flag_anion[i]     = false;
    flag_acceptor[i]  = false;
    flag_donor[i]     = false;
    flag_ring[i]      = false;
    flag_ring_aromatic[i]    = false;
    flag_ring_nonaromatic[i] = false;
    flag_in_pharmacophore[i] = false;
    radii[i]          = 1.50;
  }
  for (int i=1;i<=nrings;i++){
    ringradii[i]=1.4;
  }
  arrays_allocated = true;
}
/*************************************/
void
SmallMolecule::clear_molecule(){
  if(arrays_allocated){
    delete[]flag_hydroxyl;    
    delete[]flag_cation;
    delete[]flag_anion;
    delete[]flag_acceptor;
    delete[]flag_donor;
    delete[]flag_ring;
    delete[]flag_ring_aromatic;
    delete[]flag_ring_nonaromatic;    
    delete[]flag_in_pharmacophore;
    delete[]radii;	
    delete[]ringradii;
    hydrophobic_groups.clear();
  }
  arrays_allocated = false;
}
/*************************************/
