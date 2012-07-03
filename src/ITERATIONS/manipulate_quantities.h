#ifndef MANIPULATE_QUANTITIES_H
#define MANIPULATE_QUANTITIES_H

#include "../COMPUTE/DataUser.h"

#include "../DEFINITIONS/Process.h"
#include "../DEFINITIONS/Structure.h"
#include "../DEFINITIONS/MacroProblem.h"

using namespace LMT;


//allocation des differentes quantites
/** \ingroup   Strategie_iterative
 * \brief Allocation de mémoire pour les sous-structures, interfaces et vecteurs macro. 
 * 
 * Pour chaque pas de temps les différents vecteurs sont initialisés à 0. On n'alloue pas la mémoire pour les différents vecteurs en chaque pas de temps des ssts car ils ne sont stockés qu'à convergence.
 */
void allocate_quantities_Sst_Inter(PointedSubstructures &SubS, PointedInterfaces &SubI,Process &process);

///allocation de t_post uniquement pour les sorties de resultat
void allocate_quantities_post(PointedSubstructures &SubS, PointedInterfaces &SubI,Process &process);


void assign_quantities_current_to_old(PointedSubstructures &SubS, PointedInterfaces &SubI, Process &process);

void assign_t_post(PointedSubstructures &SubS, PointedInterfaces &SubI, Process &process);

/** \ingroup   Latin
 * \brief Pour la reprise d'un calcul, on recopie les donnees relues dans les quantités old
 * 
 */
void recopie_old_from_new(VecInterfaces &SubI,Process &process);

/** \ingroup   Latin
 * \brief Pour la reprise d'un calcul, on recopie les donnees relues dans les quantités old
 * 
 */
void recopie_old_from_new_post(VecInterfaces &SubI,Process &process);

void rebuild_state(Sst &S,Sst::Time &t, DataUser &data_user);

#endif //MANIPULATE_QUANTITIES_H