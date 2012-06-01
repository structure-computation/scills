//librairies Hugo
#include "../../../LMT/include/containers/mat.h"
#include "../../../LMT/include/containers/vec.h"
#include "../../../LMT/include/containers/vecpointedvalues.h"

//fichiers de definition des variables 
#include "../../DEFINITIONS/Process.h"
#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/Interface.h"
#include "../../DEFINITIONS/MultiScaleParameters.h"

// fonctions speciales math et autre
#include "../../UTILITAIRES/algebre.h"
#include "../../UTILITAIRES/utilitaires.h"

#include "op_inter.h"

using namespace LMT;

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
\brief Proc�dure principale permettant la cr�ation des op�rateurs d'interface
 */
// calcul des matrices de masse et de sousintegration pour chaque interface
void create_op_INTER(Vec<VecPointedValues<Sst> > &S,Vec<Interface> &Inter,Vec<VecPointedValues<Interface> > &SubI,Process &process)
{
    // calcul des matrices de masse et de sousintegration pour chaque interface
    if (process.parallelisation->is_master_cpu()) std::cout << "\t Matrice de masse et sous-integration " << std::endl;
    apply_mt(SubI,process.parallelisation->nb_threads,CalcMN(),S);

    if (process.parallelisation->is_master_cpu()) std::cout << "\t Correspondance ddl de chaque cote de l'interface" << std::endl;       
    apply_mt(SubI,process.parallelisation->nb_threads,Corresp_ddlinter());

    if (process.parallelisation->is_master_cpu()) std::cout << "\t Centre de gravite des interfaces" << std::endl;

    if (process.multiscale->multiechelle==1){  
        if (process.parallelisation->is_master_cpu()) std::cout << "\t Calcul des BPI" << std::endl;
        // calcul Base Principale d'Inertie
        apply_mt(SubI,process.parallelisation->nb_threads,CalcBPI());
        // Application du nombre de fct de base macro par interface 
        apply_mt(SubI,process.parallelisation->nb_threads,Apply_nb_macro(),process);
        // creation des projecteurs macro et micro
        if (process.parallelisation->is_master_cpu()) std::cout << "\t Calcul des projecteurs macro " << std::endl;
        apply_mt(SubI,process.parallelisation->nb_threads,CreateProjMacro(),process);
    }
};

