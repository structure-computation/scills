#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "../UTILS/Sc2String.h"
#include "../../LMT/include/containers/vec.h"
//#include "../../LMT/include/codegen/codegen.h"
//#include "../../LMT/include/containers/basicops.h"
//#include "Process.h"
//#include "TimeParameters.h"

/** \defgroup Conditions_limites
Classe relative aux conditions limites
*/

using namespace LMT;
//******************************************
// Condition aux limites
//*******************************************
/** \ingroup Conditions_limites
\brief Classe definissant les conditions aux limites
*/
struct Boundary {
    Boundary() {}
    int id;                             ///<  id de la condition limite, la même que dans data_user
    Sc2String comp;                        ///< type de condition aux limites
    Vec<Vec<TYPEREEL,DIM>, 2> box, box1;     ///< boite incluant la CL dans laquelle sont cherchees les interfaces, box1 la boite de la 2eme interface pour les périodiques
    Vec<Vec<Sc2String,DIM> > fcts_spatiales; ///< fonctions chargement spatial : 2 ou 3 fcts analytiques (x, y, z) separees par des ; (premier vecteur égal au nombre de step)
    Vec<Sc2String > fcts_temporelles;      ///< fonctions chargement temporel : 1 fct de "t" par step pour chaque coordonnée      
    Vec<Vec<TYPEREEL,2> > intervalles_temps;   ///< intervalles des differentes fonctions temporelles
    Vec<TYPEREEL,DIM> ft;                      ///< valeur de la fonction temporelle pour le piquet de temps considere pour chaque coordonnee
    Vec<unsigned> sst_num; ///< numero des sous-structures où chercher la CL
};

#endif //BOUNDARY_H
