#ifndef CREATE_OP_INTER_H
#define CREATE_OP_INTER_H

//librairies Hugo
#include "../../../LMT/include/containers/mat.h"
#include "../../../LMT/include/containers/vec.h"
#include "../../../LMT/include/containers/vecpointedvalues.h"

using namespace LMT;

//fichiers de definition des variables 
#include "Process.h"
#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/Interface.h"

/** \defgroup Operateurs_inter  Creation des operateurs par interface
\ingroup Operateurs
 
 Cette procédure permet de construire les opérateurs relatifs aux interfaces, à savoir la matrice de masse M et de sous-intégration N le projecteur macro PM et base macro eM, les normales.
 
 Elle fait appel aux fonctions définies soit dans op_inter.h.
 Le déroulement des opérations est le suivant :
 - calcul de la matrice de masse et de sous-intégration : CalcMN
 - recherche des correspondances entre les ddls de part et d'autre de l'interface (non nécessaire pour des maillages compatibles créés automatiquement) : Corresp_ddlinter
 - détermination des caractéristiques géométriques de l'interface : CalcMeasure_G
 - En multiéchelle seulement, 
     - création de la base principale d'inertie : CalcBPI
     - création des projecteurs macro et micro : CreateProjMacro
 - Création des normales  : CalcNormales
*/

/** \ingroup Operateurs_inter
 * \brief Procédure principale permettant la création des opérateurs d'interface
 * 
 * Calcul des matrices de masse et de sousintegration pour chaque interface
 */
void create_op_INTER(Vec<VecPointedValues<Sst> > &S,Vec<Interface> &Inter,Vec<VecPointedValues<Interface> > &SubI,Process &process);

#endif // CREATE_OP_INTER_H