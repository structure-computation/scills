#ifndef ITERATION_DECLARATIONS_H
#define ITERATION_DECLARATIONS_H


#include "../../LMT/include/containers/vec.h"
#include "../../LMT/include/containers/vecpointedvalues.h"

#include "../COMPUTE/FieldStructureUser.h"
#include "../COMPUTE/DataUser.h"

#include "../DEFINITIONS/Process.h"
#include "../DEFINITIONS/Sst.h"
#include "../DEFINITIONS/Interface.h"
#include "../DEFINITIONS/Boundary.h"
#include "../DEFINITIONS/MacroProblem.h"

using namespace LMT;
using namespace Metil;


/** \defgroup Strategie_iterative Strategie Iterative
\brief Strategie iterative                                                 *

Dans cette partie, on cherche a resoudre le probleme de statique ou quasistatique. Deux methodes sont envisagees : 
- on effectue une resolution latin a savoir une boucle iterative dans laquelle on cherche la solution de maniere incrementale sur tous les pas de temps : multiscale_iterate_latin(). 
- on effectue la boucle en temps et pour chaque piquet de temps on cherche une solution de maniere iterative : multiscale_iterate_incr()

*/

/** \ingroup   LATIN
\brief Procedure principale iterative pour la resolution LATIN             *

Cette procedure est constituee des fonctions suivantes :
- allocate_quantities() : On alloue dans un premier temps la memoire pour chaque quantites de Sst, Interface et Global pour chaque piquet de temps:
- assign_CL_values_space_time() : on assigne les valeurs des quantites chapeau sur le bord a partir des conditions aux limites pour chaque piquet de temps.
- iterate_latin() : on applique la strategie iterative latin : \ref LATIN*
*/
void multiscale_iterate_latin(Vec<Sst>                          &S,
                              Vec<VecPointedValues<Sst> >       &SubS, 
                              Vec<Interface>                    &Inter, 
                              Vec<VecPointedValues<Interface> > &SubI, 
                              Process                           &process, 
                              MacroProblem                      &Global, 
                              Vec<Boundary>                     &CL, 
                              DataUser                          &data_user);

/** \ingroup   Incrementale
\brief Procedure principale iterative pour la resolution du probleme de maniere incre*mentale

Cette procedure est constituee des etapes suivantes :
- allocate_quantities() : On alloue dans un premier temps la memoire pour chaque quantites de Sst, Interface et Global (pas de temps 0 et 1):
- pour chaque piquet de temps :
- assign_CL_values_space_time() : on assigne les valeurs des quantites chapeau sur le bord a partir des conditions aux limites pour le piquet de temps concerne
- iterate_incr() : on effectue une boucle iterative \ref Incrementale
- assign_quantities_current_to_old() : on met a jour les valeurs 0 a partir des valeurs 1
*/
void multiscale_iterate_incr(Vec<Sst>                          &S,
                             Vec<VecPointedValues<Sst> >       &SubS, 
                             Vec<Interface>                    &Inter, 
                             Vec<VecPointedValues<Interface> > &SubI, 
                             Process                           &process, 
                             MacroProblem                      &Global, 
                             Vec<Boundary>                     &CL, 
                             DataUser                          &data_user, 
                             GeometryUser                      &geometry_user, 
                             FieldStructureUser                &field_structure_user);

#include "manipulate_quantities.h"

#endif //ITERATION_DECLARATIONS_H