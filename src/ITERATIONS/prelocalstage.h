#ifndef PRELOCALSTAGE_H
#define PRELOCALSTAGE_H


#include "../DEFINITIONS/Process.h"

struct Interface;//#include "../DEFINITIONS/Interface.h"

#include "../DEFINITIONS/Boundary.h"

#include "../COMPUTE/DataUser.h"

#include "../../LMT/include/containers/vec.h"
#include "../../LMT/include/containers/vecpointedvalues.h"

using namespace LMT;

//************************************************************************************************************
//Procedures faisant intervenir le temps : modulation des fonctions spatiales par une fonction scalaire temporelle
//************************************************************************************************************

/** Initialise les grandeurs chapeau des interfaces au valeurs initiales donnees
 * Initialise Wchap dans le cas d'un deplacement impose
 */
void initialise_CL_values_space_time(Vec<VecPointedValues<Interface> > &Inter, Vec<Boundary> &CL, Process &process, DataUser &data_user );

/** \ingroup   LATIN
 \brief Assignation pour tous les pas de temps des conditions limites aux quantites chapeau des interfaces de bord
 
 Pour les interfaces a effort impose, on utilise assign_CL_eff().
 Pour les interfaces a deplacement impose, on utilise assign_CL_depl().
 Pour les interfaces de type symetrie, on utilise assign_CL_sym()
 Pour les interfaces de type deplacement normal donne, on utilise assign_CL_depl_normal()
 Pour les interfaces interieures de type contact avec jeu, on utilise assign_jeu()
 Sinon une erreur est envoyee.
 */
void assign_CL_values_space_time_latin(Vec<VecPointedValues<Interface> > &Inter, Vec<Boundary> &CL, Process &process, DataUser &data_user);


/** \ingroup   Incrementale
 \brief Assignation pour chaque pas de temps des conditions limites aux quantites chapeau des interfaces de bord
 
 Pour les interfaces a effort impose, on utilise assign_CL_eff().
 Pour les interfaces a deplacement impose, on utilise assign_CL_depl().
 Pour les interfaces de type symetrie, on utilise assign_CL_sym()
 Pour les interfaces de type deplacement normal donne, on utilise assign_CL_depl_normal()
 Pour les interfaces interieures de type contact avec jeu, on utilise assign_jeu() uniquement pour le piquet de temps initial
 Sinon une erreur est envoyee.
 */
void assign_CL_values_space_time_incr(Vec<VecPointedValues<Interface> > &Inter, Vec<Boundary> &CL, Process &process, DataUser &data_user );


#endif //PRELOCALSTAGE_H
