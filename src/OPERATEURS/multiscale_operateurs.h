#ifndef MULTISCALE_OPERATEURS_H
#define MULTISCALE_OPERATEURS_H

#include "../COMPUTE/DataUser.h"
#include "../GEOMETRY/GeometryUser.h"

#include "../DEFINITIONS/Param.h"
#include "../DEFINITIONS/Glob.h"
#include "../DEFINITIONS/Sst.h"
#include "../DEFINITIONS/Interface.h"

#include "../../LMT/include/containers/vec.h"
#include "../../LMT/include/containers/vecpointedvalues.h"

using namespace Metil;
using namespace LMT;


/** \defgroup Operateurs Création des opérateurs
 \brief  Creation des opérateurs pour le problème sous-structuré         *
 
L'ordre d'utilisation des fonctions est le suivant :
- operateurs_multiechelle() : fonction principale valable en 2D ou 3D
    - create_op_INTER() : création des opérateurs d'interface \ref Operateurs_inter
        - CalcMN : calcul des matrices de masse M et sous-intégration N
        - Corresp_ddlinter : correspondance des ddls de chaque coté de l'interface
        - Apply_nb_macro : application du nombre de fonctions macro par interface
        - CalcBPI : calcul de la base principale
        - CreateProjMacro : création des projecteurs macro 
    - create_op_SST() : création des opérateurs par sous-structure \ref Operateurs_sst 
        - assign_val_temporelle ou assign_bool_stat : assignation des paramètres temporels aux sst
        - Calc_SST_rigidite_K0 : calcul de la matrice de raideur EF
        - Calc_SST_Correspddl : calcul des correspondances entre les ddls de bord et ddls de la sst
        - optimise_direction : optimisation des directions de recherche
        - Calc_SST_rigidite_K0_k : ajout des directions de recherche
        - repere_ind_interface_LE : Repérage des ddl macro de bord des Ssts dans l'opérateur homogénéisé
        - Calc_SST_LE : Création de l'opérateur homogénéisé
    - create_op_MACRO() : création des opérateurs macro globaux \ref Operateurs_macro
        - Repere_ddl_Inter() : repérage des ddls des interfaces dans le problème macro 
        - Assem_prob_macro() : Assemblage du problème macro 
        - macro_CL() : prise en compte des conditions aux limites 
        - penalisation() get_factorization() : pénalisation puis factorisation de la matrice macro 
 
 */
void multiscale_operateurs(Vec<VecPointedValues<Sst> >       &Stot,
                           Vec<VecPointedValues<Sst> >       &SubS, 
                           Vec<Interface>                    &Inter, 
                           Vec<VecPointedValues<Interface> > &SubI, 
                           Param                             &process, 
                           Glob                              &Global, 
                           DataUser                          &data_user);

#endif //MULTISCALE_OPERATEURS_H