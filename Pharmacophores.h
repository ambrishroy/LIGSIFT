#ifndef pharmacophores_H
#define pharmacophores_H 
#include <iosfwd>
#include <string>
#include <vector>
#include "Coordinate.h"

using namespace std;

enum FuncGroup{
  AROM,  ///< Aromatic atoms,
  NAROM, ///< Non-aromatic atoms,
  OHG,   ///< Hydroxyl, 
  HDON,  ///< Hydrogen donor,
  HACC,  ///< Hydrogen acceptor,
  POSC,  ///< Positive charge,
  NEGC,  ///< Negative charge,
  HPHOB, ///< Hydrophobic,
  EXCL,  ///< Excluded volume 
  UNDEF, ///< Undefined value (used for initialisation)
};
/*************************************/
class  PharmacophorePoint{
 public:		
  /*#### Member variables*/
  Coordinate     point;   // coordinates of the pharmacophoric point
  FuncGroup      func;    // type of functional group
  Coordinate     normal;  // coordinates of the directionality of the pharmacophore
  double         alpha;   // spread of the gaussian is inverse proportional to radius squared
  double         vol;
  double         pcharge;
  double	 vdw;
  double         interactions;
  double         epotential;
  bool           hasNormal;  // does the pharmacophoric point has directionality  
  string         name;
  bool           includedBefore;
  /*#### Member functions */
  PharmacophorePoint(void); //Default Constructor (intialize to NULL)
  PharmacophorePoint(const PharmacophorePoint&);
  PharmacophorePoint(const PharmacophorePoint*);   
};
/*************************************/
typedef std::vector< vector<PharmacophorePoint> > Pharmacophore;

#endif  // pharmacophores_H
