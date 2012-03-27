#ifndef PRELOCALSTAGE_H
#define PRELOCALSTAGE_H


#include "../DEFINITIONS/Process.h"

struct Interface;//#include "../DEFINITIONS/Interface.h"

#include "../DEFINITIONS/Boundary.h"

#include "../COMPUTE/DataUser.h"

#include "../../LMT/include/containers/vec.h"
#include "../../LMT/include/containers/vecpointedvalues.h"
#include "../../LMT/include/codegen/codegen.h"
#include "../../LMT/include/containers/basicops.h"

using namespace LMT;
using namespace Codegen;

//************************************************************************************************************
//Procedures faisant intervenir le temps : modulation des fonctions spatiales par une fonction scalaire temporelle
//************************************************************************************************************

void calc_CL_time(Process &process,Vec<Boundary> &CL, DataUser &data_user );


///assignation de valeurs variables en fonction des variables x, y et z et du temps t
void assign_CL_spatial_temporel(Vec<TYPEREEL> &V, Vec<Vec<TYPEREEL,DIM> > &nodeeq, Boundary &CL, int i_step, DataUser &data_user);


/** \ingroup   LATIN
 \brief Assignation pour tous les pas de temps des conditions limites aux quantités chapeau des interfaces de bord
 
 Pour les interfaces à effort imposé, on utilise assign_CL_eff().
 Pour les interfaces à déplacement imposé, on utilise assign_CL_depl().
 Pour les interfaces de type symétrie, on utilise assign_CL_sym()
 Pour les interfaces de type déplacement normal donné, on utilise assign_CL_depl_normal()
 Pour les interfaces intérieures de type contact avec jeu, on utilise assign_jeu()
 Sinon une erreur est envoyée.
 */
void assign_CL_values_space_time_latin(Vec<VecPointedValues<Interface> > &Inter, Vec<Boundary> &CL, Process &process, DataUser &data_user);


/** \ingroup   Incrementale
 \brief Assignation pour chaque pas de temps des conditions limites aux quantités chapeau des interfaces de bord
 
 Pour les interfaces à effort imposé, on utilise assign_CL_eff().
 Pour les interfaces à déplacement imposé, on utilise assign_CL_depl().
 Pour les interfaces de type symétrie, on utilise assign_CL_sym()
 Pour les interfaces de type déplacement normal donné, on utilise assign_CL_depl_normal()
 Pour les interfaces intérieures de type contact avec jeu, on utilise assign_jeu() uniquement pour le piquet de temps initial
 Sinon une erreur est envoyée.
 */
void assign_CL_values_space_time_incr(Vec<VecPointedValues<Interface> > &Inter, Vec<Boundary> &CL, Process &process, DataUser &data_user );


#endif //PRELOCALSTAGE_H
