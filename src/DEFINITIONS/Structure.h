#ifndef Structure_H
#define Structure_H
//librairie Hugo
#include "containers/mat.h"
#include "containers/vecpointedvalues.h"

#include "definition_PARAM.h"
#include "definition_SST_time.h"
#include "definition_INTER_time.h"
#include "definition_materials_property.h"
#include "definition_GLOB.h"
#include "definition_CL.h"

#include "allocation_memoire_PARAM.h"
#include "SCtime.h"

//biblioteque venant de SC_create_2
#include <Metil/BasicVec.h>
#include <Metil/StructCompactor.h>
#include <Metil/CudaMetil.h>
#include <Metil/Hdf.h>

#include "GeometryUser.h"
#include "DataUser.h"
#include "FieldStructureUser.h"


struct Structure{
    //attributs==============================================================================================
    //attribut de mise en données
    Param process;
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
    Glob Global;                                // propriété macro de la structure
    



    //méthodes================================================================================================
    //lecture des données
    void initialisation_param_MPI(int argc,char **argv);
    void finalisation_MPI();
    void initialisation_SC_create_2(std::string id_model, std::string id_calcul);
    
    void multiscale();
    void multiscale_calculation();
    
    
};



#endif // Structure_H