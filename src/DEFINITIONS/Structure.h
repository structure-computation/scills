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
    ///attribut de mise en donnees
    Process process;
    TicTac tic1;
    
    DataUser data_user;
    GeometryUser geometry_user;
    FieldStructureUser field_structure_user;
    
    ///attribut geometriques
    Vec<Sst> S;                                 /// vecteur des sous structures
    Vec<Interface> Inter;                       /// vecteur des interfaces
    
    Vec< VecPointedValues<Sst> > SubS;          /// MPI : vecteur de pointeur vers les Sst
    Vec< VecPointedValues<Sst> > Stot;          /// MPI : vecteur de pointeur vers les Sst
    Vec< VecPointedValues<Interface> > SubI;    /// MPI : vecteur de pointeur vers les Interfaces
    
    /// attributs de comportement
    Vec<Boundary> CL;                           /// vecteur des conditions limites
    Vec<SstCarac> matprops;                     /// vecteur des propriete materiaux
    
    /// operateurs de la structure
    MacroProblem Global;                        /// propriete macro de la structure
    



    //methodes================================================================================================
    int algorithme_incremental(int argc,char **argv);
    void initialisation_MPI(int argc,char **argv);
    void finalisation_MPI();
    void lecture_fichiers(Sc2String id_model, Sc2String id_calcul);
    void chargement_donnees();
    
    void boucle_multi_resolution();
    void boucle_steps_client();
    void synchronize_and_print_duration(TicTac& tic, const char* msg);
    void display(const char* msg);
};



#endif // Structure_H