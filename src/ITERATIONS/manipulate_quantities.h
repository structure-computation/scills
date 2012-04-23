#ifndef MANIPULATE_QUANTITIES_H
#define MANIPULATE_QUANTITIES_H

#include "../COMPUTE/DataUser.h"

#include "../DEFINITIONS/Process.h"
#include "../DEFINITIONS/Sst.h"
#include "../DEFINITIONS/Interface.h"
#include "../DEFINITIONS/MacroProblem.h"

#include "../../LMT/include/containers/vec.h"
#include "../../LMT/include/containers/vecpointedvalues.h"

using namespace LMT;


//allocation des differentes quantites
/** \ingroup   Strategie_iterative
 * \brief Allocation de m�moire pour les sous-structures, interfaces et vecteurs macro. 
 * 
 * Pour chaque pas de temps les diff�rents vecteurs sont initialis�s � 0. On n'alloue pas la m�moire pour les diff�rents vecteurs en chaque pas de temps des ssts car ils ne sont stock�s qu'� convergence.
 */
void allocate_quantities(Vec<VecPointedValues<Sst > > &S, Vec<VecPointedValues<Interface > > &Inter,Process &process,MacroProblem &Global);

//allocation de t_post uniquement pour les sorties de resultat
void allocate_quantities_post(Vec<VecPointedValues<Sst > > &S, Vec<VecPointedValues<Interface > > &Inter,Process &process);


void assign_quantities_current_to_old(Vec<VecPointedValues<Sst> > &S, Vec<VecPointedValues<Interface> > &Inter, Process &process);

void assign_t_post(Vec<VecPointedValues<Sst> > &S, Vec<VecPointedValues<Interface> > &Inter, Process &process);

/** \ingroup   Latin
 * \brief Pour la reprise d'un calcul, on recopie les donnees relues dans les quantit�s old
 * 
 */
void recopie_old_from_new(Vec<Interface> &Inter,Process &process);

/** \ingroup   Latin
 * \brief Pour la reprise d'un calcul, on recopie les donnees relues dans les quantit�s old
 * 
 */
void recopie_old_from_new_post(Vec<Interface> &Inter,Process &process);

void assign_dep_cont_slave(Sst &S,Sst::Time &t, DataUser &data_user);

#endif //MANIPULATE_QUANTITIES_H