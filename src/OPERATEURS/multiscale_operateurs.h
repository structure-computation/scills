#ifndef MULTISCALE_OPERATEURS_H
#define MULTISCALE_OPERATEURS_H

#include "../COMPUTE/DataUser.h"
#include "../GEOMETRY/GeometryUser.h"

#include "../DEFINITIONS/Process.h"
#include "../DEFINITIONS/structure_typedef.h"

using namespace Metil;
using namespace LMT;


/** \defgroup Operateurs Creation des operateurs
 \brief  Creation des operateurs pour le problème sous-structure         *
 
L'ordre d'utilisation des fonctions est le suivant :
- operateurs_multiechelle() : fonction principale valable en 2D ou 3D
    - create_op_INTER() : creation des operateurs d'interface \ref Operateurs_inter
        - CalcMN : calcul des matrices de masse M et sous-integration N
        - Corresp_ddlinter : correspondance des ddls de chaque cote de l'interface
        - Apply_nb_macro : application du nombre de fonctions macro par interface
        - CalcBPI : calcul de la base principale
        - CreateProjMacro : creation des projecteurs macro 
    - create_op_SST() : creation des operateurs par sous-structure \ref Operateurs_sst 
        - assign_val_temporelle ou assign_bool_stat : assignation des paramètres temporels aux sst
        - Calc_SST_rigidite_K0 : calcul de la matrice de raideur EF
        - Calc_SST_Correspddl : calcul des correspondances entre les ddls de bord et ddls de la sst
        - optimise_direction : optimisation des directions de recherche
        - Calc_SST_rigidite_K0_k : ajout des directions de recherche
        - repere_ind_interface_LE : Reperage des ddl macro de bord des Ssts dans l'operateur homogeneise
        - Calc_SST_LE : Creation de l'operateur homogeneise
    - create_op_MACRO() : creation des operateurs macro globaux \ref Operateurs_macro
        - Repere_ddl_Inter() : reperage des ddls des interfaces dans le problème macro 
        - Assem_prob_macro() : Assemblage du problème macro 
        - macro_CL() : prise en compte des conditions aux limites 
        - penalisation() get_factorization() : penalisation puis factorisation de la matrice macro 
 
 */
void multiscale_operateurs(PointedSubstructures &Stot,
                           PointedSubstructures &SubS, 
                           VecInterfaces        &Inter, 
                           PointedInterfaces    &SubI, 
                           Process              &process, 
                           MacroProblem         &Global, 
                           DataUser             &data_user);

#endif //MULTISCALE_OPERATEURS_H