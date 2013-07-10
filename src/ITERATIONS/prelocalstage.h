#ifndef PRELOCALSTAGE_H
#define PRELOCALSTAGE_H

#include "../DEFINITIONS/structure_typedef.h"
#include "Process.h"
#include "../COMPUTE/DataUser.h"


//************************************************************************************************************
//Procedures faisant intervenir le temps : modulation des fonctions spatiales par une fonction scalaire temporelle
//************************************************************************************************************

/** Initialise les grandeurs chapeau des interfaces au valeurs initiales donnees
 * Initialise Wchap dans le cas d'un deplacement impose
 */
void initialise_CL_values(PointedInterfaces &Inter, Boundaries &CL);


/** \ingroup   Incrementale
 \brief Assignation pour chaque pas de temps des conditions limites aux quantites chapeau des interfaces de bord
 
 Pour les interfaces a effort impose, on utilise assign_CL_eff().
 Pour les interfaces a deplacement impose, on utilise assign_CL_depl().
 Pour les interfaces de type symetrie, on utilise assign_CL_sym()
 Pour les interfaces de type deplacement normal donne, on utilise assign_CL_depl_normal()
 Pour les interfaces interieures de type contact avec jeu, on utilise assign_jeu() uniquement pour le piquet de temps initial
 Sinon une erreur est envoyee.
 */
void update_CL_values(PointedInterfaces &Inter, Boundaries &CL, Process &process, DataUser &data_user );


#endif //PRELOCALSTAGE_H
