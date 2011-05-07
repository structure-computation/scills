//librairies Hugo
#include "containers/mat.h"

//#include "mesh/mesh.h"
//#include "mesh/problem.h"
//#include "mesh/displayparaview.h"
#include <fstream>
#include <map>

//#include "../entete.h"
//fichiers de definition des variables
#include "definition_PARAM.h"
#include "definition_PARAM_MULTI.h"
//#include "definition_PARAM_LATIN.h"
#include "definition_GLOB.h"
#include "definition_CL.h"


// fonctions speciales math et autre
#include "utilitaires.h"

#include "op_macro.h"

using namespace LMT;
using namespace std;


/** \defgroup Operateurs_macro Creation du problème macro global
\ingroup Operateurs
 
 Cette procédure permet de construire les opérateurs relatifs au problème macro global. Celle ci n'est accessible qu'en multiéchelle.
 
 Elle fait appel aux fonctions définies soit dans op_macro.h.
 Le déroulement des opérations est le suivant :
 - repérage des ddls des interfaces dans le problème macro : Repere_ddl_Inter()
 - Assemblage du problème macro : Assem_prob_macro()
 - prise en compte des conditions aux limites : macro_CL()
 - pénalisation puis factorisation de la matrice macro : penalisation() get_factorization()
 
*/

/** \ingroup Operateurs_macro
\brief Procédure principale permettant la création des quantités utiles au problème macro
 */
//************************************
// Procedure macro globale
//************************************
template<class TV1, class TV2, class GLOB >
void create_op_MACRO(TV1 &S, TV2 &Inter, Param &process,  GLOB &Global) {

    //reperage des ddls macro dans le probleme macro
    if (process.rank == 0)
        cout <<  "\t Reperage ddl macro" << endl;
    Repere_ddl_Inter(S,Inter,process);

    Mat<TYPEREEL, Sym<>, SparseLine<> > bigK;
    if (process.rank == 0)
        cout <<  "\t Assemblage probleme macro" << endl;
    if (process.rank == 0)
        bigK=Assem_prob_macro(S,Inter,process);
    if (process.rank == 0)
        cout <<  "\t Blocage du probleme macro" << endl;
    if (process.rank == 0)
        Global.repddlMbloq=macro_CL(Inter,process);
    //penalisation de la matrice macro
    if (process.rank == 0)
        cout <<  "\t Penalisation du probleme macro" << endl;
    if (process.rank == 0)
        penalisation(bigK,Global.repddlMbloq,Global.coefpenalisation);

    if (process.rank == 0)
        cout << "Taille du probleme macro : " << bigK.nb_rows() << endl;
    //factorisation de la matrice macro
    if (process.rank == 0)
        cout <<  "\t Factorisation matrice macro" << endl;
    //Global.bigL=inv(bigK);
//    Ecrire la matrice macro dans un fichier
/*    if (process.rank == 0) {
       for(unsigned i=0;i<bigK.nb_rows();++i) {
           if (bigK(i,i)==0){
               bigK(i,i)=Global.coefpenalisation;
               cout << "\t\tModification du terme " << i << endl;
           }
       }
    }
   string name="Kmacro";
   ofstream f(name.c_str());
   for(unsigned i=0;i<bigK.data.size();++i) {
       for(unsigned j=0;j<bigK.data[i].indices.size();j++){
           f << (i+1) << " " << (bigK.data[i].indices[j]+1) << " " << bigK.data[i].data[j] << endl ;
   }
   }*/
//    f << bigK << endl ;


     if (process.rank == 0)
       Global.l.get_factorization( bigK, true, true );


};
