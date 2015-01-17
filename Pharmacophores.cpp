#include "Pharmacophores.h"

using namespace std;

PharmacophorePoint::PharmacophorePoint(){
  name   = "NA";
  func   = UNDEF;
  alpha  = 1.0;
  vol    = 0.0;
  pcharge= 0.00;
  vdw    = 1.50;
  interactions= 0.00;
  epotential=0.0;
  normal.x = 0.0;
  normal.y = 0.0;
  normal.z = 0.0;	
  point.x = 0.0;
  point.y = 0.0;
  point.z = 0.0;  
  hasNormal      = false;
  includedBefore = false;
}
/*************************************/
PharmacophorePoint::PharmacophorePoint(const PharmacophorePoint & p){
  name   = p.name;
  point  = p.point;
  func   = p.func;
  alpha  = p.alpha;
  vol    = p.vol;
  pcharge      = p.pcharge;
  vdw          = p.vdw;
  interactions = p.interactions;
  epotential   = p.epotential;
  normal = p.normal;
  hasNormal = p.hasNormal;
  includedBefore = p.includedBefore;
}
/*************************************/
PharmacophorePoint::PharmacophorePoint(const PharmacophorePoint * p){
  name  = p->name;
  point = p->point;
  func  = p->func;
  vol   = p->vol;
  alpha = p->alpha;
  pcharge     = p->pcharge;
  vdw         = p->vdw;
  interactions= p->interactions;
  epotential  = p->epotential;
  normal = p->normal;
  hasNormal = p->hasNormal;
  includedBefore = p->includedBefore;
}
/*************************************/
