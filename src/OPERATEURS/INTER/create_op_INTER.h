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
 
 Cette proc�dure permet de construire les op�rateurs relatifs aux interfaces, � savoir la matrice de masse M et de sous-int�gration N le projecteur macro PM et base macro eM, les normales.
 
 Elle fait appel aux fonctions d�finies soit dans op_inter.h.
 Le d�roulement des op�rations est le suivant :
 - calcul de la matrice de masse et de sous-int�gration : CalcMN
 - recherche des correspondances entre les ddls de part et d'autre de l'interface (non n�cessaire pour des maillages compatibles cr��s automatiquement) : Corresp_ddlinter
 - d�termination des caract�ristiques g�om�triques de l'interface : CalcMeasure_G
 - En multi�chelle seulement, 
     - cr�ation de la base principale d'inertie : CalcBPI
     - cr�ation des projecteurs macro et micro : CreateProjMacro
 - Cr�ation des normales  : CalcNormales
*/

/** \ingroup Operateurs_inter
 * \brief Proc�dure principale permettant la cr�ation des op�rateurs d'interface
 * 
 * Calcul des matrices de masse et de sousintegration pour chaque interface
 */
void create_op_INTER(Vec<VecPointedValues<Sst> > &S,Vec<Interface> &Inter,Vec<VecPointedValues<Interface> > &SubI,Process &process);

#endif // CREATE_OP_INTER_H