#include "../DEFINITIONS/Sst.h"
#include "../DEFINITIONS/Interface.h"

#include "../../LMT/include/containers/vec.h"
#include "../../LMT/include/containers/vecpointedvalues.h"

using namespace LMT;


//allocation de t_post uniquement pour les sorties de resultat
void allocate_quantities_post(Vec<VecPointedValues<Sst > > &S, Vec<VecPointedValues<Interface > > &Inter,Param &process);


//allocation des differentes quantites
/** \ingroup   Strategie_iterative
\brief Allocation de mémoire pour les sous-structures, interfaces et vecteurs macro. 

Pour chaque pas de temps les différents vecteurs sont initialisés à 0. On n'alloue pas la mémoire pour les différents vecteurs en chaque pas de temps des ssts car ils ne sont stockés qu'à convergence.
*/
void allocate_quantities(Vec<VecPointedValues<Sst > > &S, Vec<VecPointedValues<Interface > > &Inter,Param &process,Glob &Global);

