#ifndef STRUCTURE_TYPEDEF_H
#define STRUCTURE_TYPEDEF_H

#include "main_typedef.h"
#include "Sst.h"
#include "Interface.h"
#include "Boundary.h"
#include "SstCarac_InterCarac.h"
#include "VolumicForces.h"
#include "ThermalLoad.h"


typedef LMT::Vec<Sst> Substructures;
typedef LMT::Vec<Interface> VecInterfaces;
typedef LMT::Vec< VecPointedValues<Sst> > PointedSubstructures;
typedef LMT::Vec< VecPointedValues<Interface> > PointedInterfaces;
typedef LMT::Vec<Boundary> Boundaries;
typedef LMT::Vec<SstCarac> Materials;
typedef LMT::Vec<InterCarac> Links;

#endif //STRUCTURE_TYPEDEF_H