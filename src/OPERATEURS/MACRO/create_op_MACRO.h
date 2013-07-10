//librairies Hugo
#include "../../../LMT/include/containers/vec.h"
#include "../../../LMT/include/containers/vecpointedvalues.h"

//fichiers de definition des variables
#include "Process.h"
#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/Interface.h"
#include "../../DEFINITIONS/MacroProblem.h"


using namespace LMT;


/** \defgroup Operateurs_macro Creation du probl�me macro global
\ingroup Operateurs
 
 Cette proc�dure permet de construire les op�rateurs relatifs au probl�me macro global. Celle ci n'est accessible qu'en multi�chelle.
 
 Elle fait appel aux fonctions d�finies soit dans operateur_macro.cpp.
 Le d�roulement des op�rations est le suivant :
 - rep�rage des ddls des interfaces dans le probl�me macro : Repere_ddl_Inter()
 - Assemblage du probl�me macro : Assem_prob_macro()
 - prise en compte des conditions aux limites : macro_CL()
 - p�nalisation puis factorisation de la matrice macro : penalisation(), get_factorization()
 
*/

/** \ingroup Operateurs_macro
\brief Proc�dure principale permettant la cr�ation des quantit�s utiles au probl�me macro
 */
//************************************
// Procedure macro globale
//************************************
void create_op_MACRO(Vec<VecPointedValues<Sst> > &S, Vec<Interface> &Inter, Process &process,  MacroProblem &Global);
