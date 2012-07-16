#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "main_typedef.h"
#include "../PARAMETERS/Parameters.h"

/** \defgroup Conditions_limites
Classe relative aux conditions limites
*/

//******************************************
// Condition aux limites
//*******************************************
/** \ingroup Conditions_limites
\brief Classe definissant les conditions aux limites
*/
struct Boundary {
    Boundary();
    ~Boundary();
    void read_data_user(int index,DataUser &data_user);
    void affiche();
    
    static void prepareParameters();
    static void updateParameters();
    static ParameterGroup CL_parameters;
    static MainParameter x;
    static MainParameter y;
#if DIM == 3
    static MainParameter z;
#endif
    
    Id id;                                  ///  id de la condition limite, la même que dans data_user
    Sc2String comp;                         /// type de condition aux limites
    Vec<UserParameter,DIM> fcts_spatiales;    /// Vecteur des parametres pour evaluer la CL
};

#endif //BOUNDARY_H

struct stat;
