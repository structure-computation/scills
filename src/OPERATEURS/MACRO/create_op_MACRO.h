//librairies Hugo
#include "../../../LMT/include/containers/vec.h"
#include "../../../LMT/include/containers/vecpointedvalues.h"

//fichiers de definition des variables
#include "Process.h"
#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/Interface.h"
#include "../../DEFINITIONS/MacroProblem.h"


using namespace LMT;


/** \defgroup Operateurs_macro Creation du problème macro global
\ingroup Operateurs
 
 Cette procédure permet de construire les opérateurs relatifs au problème macro global. Celle ci n'est accessible qu'en multiéchelle.
 
 Elle fait appel aux fonctions définies soit dans operateur_macro.cpp.
 Le déroulement des opérations est le suivant :
 - repérage des ddls des interfaces dans le problème macro : Repere_ddl_Inter()
 - Assemblage du problème macro : Assem_prob_macro()
 - prise en compte des conditions aux limites : macro_CL()
 - pénalisation puis factorisation de la matrice macro : penalisation(), get_factorization()
 
*/

/** \ingroup Operateurs_macro
\brief Procédure principale permettant la création des quantités utiles au problème macro
 */
//************************************
// Procedure macro globale
//************************************
void create_op_MACRO(Vec<VecPointedValues<Sst> > &S, Vec<Interface> &Inter, Process &process,  MacroProblem &Global);
