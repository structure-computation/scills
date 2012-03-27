#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/SstCarac_InterCarac.h"

using namespace LMT;

/** \ingroup Operateurs_sst_Qstat
\brief Assignation des valeurs temporelles aux sous-structures à partir du champ Time de process

Cet operateur() permet d'assigner la valeur dt (pas de temps) et theta (paramètre de la theta methode) au maillage des sous-structures pour que ceux ci soit accessible sans avoir à passer Param en argument.
 */
struct assign_val_temporelle {
    void operator() (Sst &S, TimeParameters &temps) const{
      S.matprop.dt = temps.dt;
      S.matprop.stat = 1;
      //S.mesh.theta = temps.theta;
   }
};


/** \ingroup Operateurs_sst_Qstat
\brief Assignation d'un booleen par sst pour specifier si on effectue un calcul en statique ou quasistatique

 */
struct assign_bool_stat {
   void operator() (Sst &S) const{
      S.mesh.stat = 0;
      
      //S.mesh.theta = temps.theta;
   }
};


