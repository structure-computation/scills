//librairies Hugo
#include "containers/mat.h"

//fichiers de definition des variables 
#include "Param.h"
#include "definition_PARAM_MULTI.h"

// fonctions speciales math et autre
#include "algebre.h"
#include "utilitaires.h"

#include "containers/vec_mt.h"
#include "op_inter.h"

using namespace LMT;
using namespace std;

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
\brief Procédure principale permettant la création des opérateurs d'interface
 */
// calcul des matrices de masse et de sousintegration pour chaque interface
template<class TV1, class TV2, class TV3> void create_op_INTER(TV1 &S,TV2 &Inter,TV3 &SubI,Param &process)
{

  // calcul des matrices de masse et de sousintegration pour chaque interface
   if (process.rank == 0) std::cout << "\t Matrice de masse et sous-integration " << std::endl;
   apply_mt(SubI,process.nb_threads,CalcMN(),S);
  
   if (process.rank == 0) std::cout << "\t Correspondance ddl de chaque cote de l'interface" << std::endl;       
   apply_mt(SubI,process.nb_threads,Corresp_ddlinter());
   
   if (process.rank == 0) std::cout << "\t Centre de gravite des interfaces" << std::endl;
  
   if (process.multiscale->multiechelle==1){  
      if (process.rank == 0) std::cout << "\t Calcul des BPI" << std::endl;
     // calcul Base Principale d'Inertie
      apply_mt(SubI,process.nb_threads,CalcBPI());
     // Application du nombre de fct de base macro par interface 
      apply_mt(SubI,process.nb_threads,Apply_nb_macro(),process);
     // creation des projecteurs macro et micro
      if (process.rank == 0) std::cout << "\t Calcul des projecteurs macro " << std::endl;
      apply_mt(SubI,process.nb_threads,CreateProjMacro(),process);
   }
   if (process.rank == 0) std::cout << std::endl;
};

