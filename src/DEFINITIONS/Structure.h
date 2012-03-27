#ifndef Structure_H
#define Structure_H
//librairie Hugo
#include "../../LMT/include/containers/mat.h"
#include "../../LMT/include/containers/vecpointedvalues.h"

#include "Process.h"
#include "Sst.h"
#include "Interface.h"
#include "SstCarac_InterCarac.h"
#include "MacroProblem.h"
#include "Boundary.h"

#include "../UTILS/SCtime.h"
#include "../UTILS/Sc2String.h"

//biblioteque venant de SC_create_2
#include <Metil/BasicVec.h>
#include <Metil/StructCompactor.h>
#include <Metil/CudaMetil.h>
#include <Metil/Hdf.h>

#include "../GEOMETRY/GeometryUser.h"
#include "../COMPUTE/DataUser.h"
#include "../COMPUTE/FieldStructureUser.h"


struct Structure{
    //attributs==============================================================================================
    //attribut de mise en données
    Process process;
    TicTac tic1;
    
    DataUser data_user;
    GeometryUser geometry_user;
    FieldStructureUser field_structure_user;
    
    //attribut géométriques
    Vec<Sst> S;                                 // vecteur des sous structures
    Vec<Interface> Inter;                       // vecteur des interfaces
    
    Vec< VecPointedValues<Sst> > SubS;          // MPI : vecteur de pointeur vers les Sst
    Vec< VecPointedValues<Sst> > Stot;          // MPI : vecteur de pointeur vers les Sst
    Vec< VecPointedValues<Interface> > SubI;    // MPI : vecteur de pointeur vers les Interfaces
    
    // attributs de comportement
    Vec<Boundary> CL;                           // vecteur des conditions limites
    Vec<SstCarac> matprops;                     // vecteur des propriété matériaux
    
    // opérateurs de la structure
    MacroProblem Global;                        // propriété macro de la structure
    



    //méthodes================================================================================================
    //lecture des données
    void initialisation_MPI(int argc,char **argv);
    void finalisation_MPI();
    void lecture_fichiers(Sc2String id_model, Sc2String id_calcul);
    void chargement_donnees();
    
    void multiscale();
    void multiscale_calculation();
    void synchronize_and_print_duration(TicTac& tic, const char* msg);
    
};



#endif // Structure_H