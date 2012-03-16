#ifndef CL_H
#define CL_H

#include "codegen/codegen.h"
#include "containers/basicops.h"
#include "definition_PARAM.h"
#include "definition_PARAM_TEMPS.h"

/** \defgroup Conditions_limites
Classe relative aux conditions limites
*/

using namespace Codegen;
using namespace std;
using namespace LMT;
//******************************************
// Condition aux limites
//*******************************************
/** \ingroup Conditions_limites
\brief Classe parametrable definissant les conditions aux limites
*/
struct Boundary {
    Boundary() {
      ////modif DAVID 02-09-2007
      //ft=1;
      ////fin modif
    }
    static const unsigned dim = DIM;   ///< dimension utilisable pour toutes les instances
    typedef Vec<TYPEREEL,dim> Pvec;            ///< type des points
    int id;                             ///<  id de la condition limite, la même que dans data_user
    string comp;                        ///< type de condition aux limites
    Vec<Pvec, 2> box, box1;             ///< boite incluant la CL dans laquelle sont cherchees les interfaces, box1 la boite de la 2eme interface pour les périodiques
    Vec<Vec<string,dim> > fcts_spatiales; ///< fonctions chargement spatial : 2 ou 3 fcts analytiques (x, y, z) separees par des ; (premier vecteur égal au nombre de step)
    Vec<string > fcts_temporelles;      ///< fonctions chargement temporel : 1 fct de "t" par step pour chaque coordonnée      
    Vec<Vec<TYPEREEL,2> > intervalles_temps;   ///< intervalles des differentes fonctions temporelles
    //// modif DAVID 02-09-2007
    Vec<TYPEREEL,dim> ft;                      ///< valeur de la fonction temporelle pour le piquet de temps considere pour chaque coordonnee
    //// fin modif
    Vec<unsigned> sst_num; ///< numero des sous-structures où chercher la CL
};

#endif //CL_H
