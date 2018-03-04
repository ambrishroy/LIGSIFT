#ifndef smallmolecule_H
#define smallmolecule_H 
#include <iosfwd>
#include <string>
#include <vector>
#include <set>
#include <list>
#include <math.h>
#include <stdio.h>
#include <openbabel/atom.h>
#include <openbabel/mol.h>
#include <openbabel/ring.h>
#include <openbabel/data.h>
#include <openbabel/math/vector3.h>
#include "Coordinate.h"
#include "Pharmacophores.h"
#include "Common.h"

class SmallMolecule{
 public:
  void           getMoleculeData(OpenBabel::OBMol& OBmol);
  void           InitPharmacophore(OpenBabel::OBMol&, vector<PharmacophorePoint>&);
  void           clear_molecule();
  SmallMolecule(): arrays_allocated(false){ } //initialize them in the constructor
//  Pharmacophore	 pharmacophore;
 private:
  
  void           allocate_arrays(int , int, int);
  bool           arrays_allocated;
  bool           *flag_hydroxyl;
  bool           *flag_cation;
  bool           *flag_anion;
  bool           *flag_acceptor;
  bool           *flag_donor;
  bool           *flag_ring_aromatic;
  bool           *flag_ring_nonaromatic;
  bool           *flag_ring;
  bool           *flag_in_pharmacophore;
  double         *radii;
  double         *ringradii;
  std::vector< std::vector<int> > hydrophobic_groups;			
 
  Coordinate		CalcNormal(OpenBabel::OBAtom*);  
};



#endif  // smallmolecule_H
