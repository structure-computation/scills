//librairies Hugo
#include "../../../LMT/include/containers/mat.h"
#include "../../../LMT/include/containers/vec.h"
#include "../../../LMT/include/containers/vecpointedvalues.h"

//#include "mesh/mesh.h"
//#include "mesh/problem.h"
//#include "mesh/displayparaview.h"
#include <iostream>

//fichiers de definition des variables
#include "../../DEFINITIONS/Process.h"
#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/MultiScaleParameters.h"
#include "../../DEFINITIONS/MacroProblem.h"
#include "../../DEFINITIONS/Boundary.h"


// fonctions speciales math et autre
#include "../../UTILITAIRES/utilitaires.h"

#include "op_macro.h"

using namespace LMT;


/** \defgroup Operateurs_macro Creation du problème macro global
\ingroup Operateurs
 
 Cette procédure permet de construire les opérateurs relatifs au problème macro global. Celle ci n'est accessible qu'en multiéchelle.
 
 Elle fait appel aux fonctions définies soit dans op_macro.h.
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
void create_op_MACRO(Vec<VecPointedValues<Sst> > &S, Vec<Interface> &Inter, Process &process,  MacroProblem &Global) {
    ///Reperage des ddls macro dans le probleme macro pour tous le monde
    process.print("\t Reperage ddl macro");
    Repere_ddl_Inter(S,Inter,process);
    
    /// Puis seul le master (processeur 0) se charge du pb macro
    if(process.parallelisation->is_master_cpu()){
        
        /// Creation de la matrice de raideur macroscopique
        Mat<TYPEREEL, Sym<>, SparseLine<> > bigK;
        process.print("\t Assemblage probleme macro");
        Assem_prob_macro(S,Inter,process,bigK);
        
        /// Blocage des ddls imposes du probleme macro
        process.print("\t Blocage du probleme macro");
        macro_CL(Inter,process,Global.repddlMbloq);
        
        /// Penalisation de la matrice macro
        process.print("\t Penalisation du probleme macro");
        penalisation(bigK,Global.repddlMbloq,Global.coefpenalisation);

        process.print_data("Taille du probleme macro : ",bigK.nb_rows());
        
        /// Factorisation de la matrice macro
        process.print("\t Factorisation matrice macro");
        Global.l.get_factorization( bigK, true, true );
    }
};
