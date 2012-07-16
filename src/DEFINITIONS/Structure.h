#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "structure_typedef.h"

struct Structure{
    //attributs==============================================================================================
    /// Attributs geometriques
    Substructures S;            /// vecteur des sous structures
    VecInterfaces Inter;        /// vecteur des interfaces
    
    /// Attributs pour MPI
    PointedInterfaces SubS;     /// MPI : vecteur de pointeur vers les Sst
    PointedInterfaces Stot;     /// MPI : vecteur de pointeur vers les Sst
    PointedInterfaces SubI;     /// MPI : vecteur de pointeur vers les Interfaces
    
    /// attributs de comportement
    Materials materials;        /// vecteur des propriete materiaux
    Links links;                /// vecteur des proprietes des liaisons
    
    /// Attributs de chargement
    Boundaries CL;              /// vecteur des conditions limites
    VolumicForces Fvol;         /// efforts volumiques sur la structure
};


#endif // STRUCTURE_H