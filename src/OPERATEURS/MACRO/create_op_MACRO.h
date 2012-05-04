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
    if(process.rank == 0) std::cout <<  "\t Reperage ddl macro" << endl;
    Repere_ddl_Inter(S,Inter,process);
    
    /// Puis seul le master (processeur 0) se charge du pb macro
    if(process.rank == 0){
        /// Creation de la matrice de raideur macroscopique
        Mat<TYPEREEL, Sym<>, SparseLine<> > bigK;
        std::cout <<  "\t Assemblage probleme macro" << endl;
        Assem_prob_macro(S,Inter,process,bigK);
        
        /// Blocage des ddls imposes du probleme macro
        std::cout <<  "\t Blocage du probleme macro" << endl;
        macro_CL(Inter,process,Global.repddlMbloq);
        
        /// Penalisation de la matrice macro
        std::cout <<  "\t Penalisation du probleme macro" << endl;
        penalisation(bigK,Global.repddlMbloq,Global.coefpenalisation);

        std::cout << "Taille du probleme macro : " << bigK.nb_rows() << endl;
        
        /// Factorisation de la matrice macro
        std::cout <<  "\t Factorisation matrice macro" << endl;
        Global.l.get_factorization( bigK, true, true );
    }
};
