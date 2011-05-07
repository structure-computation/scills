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
template<unsigned dim_, class TT_>
struct Boundary {
    Boundary() {
      ////modif DAVID 02-09-2007
      //ft=1;
      ////fin modif
    }
    static const unsigned dim = dim_;///< dimension utilisable pour toutes les instances
    typedef TT_ T; ///< type des flottants de la classe
    typedef Vec<T,dim> Pvec;///< type des points
    string comp; ///< type de condition aux limites
    Vec<Pvec, 2> box, box1; ///< boite incluant la CL dans laquelle sont cherchees les interfaces, box1 la boite de la 2eme interface pour les périodiques
    Vec<string,dim> fcts_spatiales; ///< fonctions chargement spatial : 2 ou 3 fcts analytiques (x, y, z) separees par des ;
    Vec<Vec<string> > fcts_temporelles; ///< fonctions chargement temporel : 1 fct de "t"
    Vec<Vec<T,2> > intervalles_temps; ///< intervalles des differentes fonctions temporelles
    //// modif DAVID 02-09-2007
    Vec<T> ft; ///< valeur de la fonction temporelle pour le piquet de temps considere pour chaque coordonnee
    //T ft; ///< valeur de la fonction temporelle pour le piquet de temps considere
    //// fin modif
   Vec<unsigned> sst_num; ///< numero des sous-structures où chercher la CL
};

#endif //CL_H
