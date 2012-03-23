#include "../DEFINITIONS/Sst.h"
#include "../DEFINITIONS/Interface.h"

#include "../../LMT/include/containers/vec.h"
#include "../../LMT/include/containers/vecpointedvalues.h"

using namespace LMT;


//allocation de t_post uniquement pour les sorties de resultat
void allocate_quantities_post(Vec<VecPointedValues<Sst > > &S, Vec<VecPointedValues<Interface > > &Inter,Param &process);


//allocation des differentes quantites
/** \ingroup   Strategie_iterative
\brief Allocation de m�moire pour les sous-structures, interfaces et vecteurs macro. 

Pour chaque pas de temps les diff�rents vecteurs sont initialis�s � 0. On n'alloue pas la m�moire pour les diff�rents vecteurs en chaque pas de temps des ssts car ils ne sont stock�s qu'� convergence.
*/
void allocate_quantities(Vec<VecPointedValues<Sst > > &S, Vec<VecPointedValues<Interface > > &Inter,Param &process,Glob &Global);

