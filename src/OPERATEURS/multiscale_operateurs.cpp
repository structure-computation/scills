#include "multiscale_operateurs.h"

//librairies Hugo
#include "../../LMT/include/containers/mat.h"
#include "../../LMT/include/containers/vecpointedvalues.h"
#include "../../LMT/include/containers/allocator.h"
#include "../UTILS/SCtime.h"
#ifndef INFO_TIME
#define INFO_TIME
#endif 

// fonction utilisees pour la creation des operateurs
#include "INTER/create_op_INTER.h"
#include "MACRO/create_op_MACRO.h"
#include "SST/create_op_SST.h"
#include "../DEFINITIONS/main_typedef.h"
#include "../../LMT/include/containers/matcholamd.h"

//pour l'affichage : inclure ce .h
// #include "affichage.h"

#include "mpi.h"
#include "../DEFINITIONS/MultiResolutionData.h"



using namespace LMT;
extern Crout crout;


/** \defgroup Operateurs Création des opérateurs
\brief  Creation des opérateurs pour le problème sous-structuré
 
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

//***************************************************************************************************
//  Procedure de creation des operateurs utilises dans la strategie multiechelle
//***************************************************************************************************
void multiscale_operateurs(PointedSubstructures &Stot,
                           PointedSubstructures &SubS, 
                           VecInterfaces        &Inter, 
                           PointedInterfaces    &SubI, 
                           Process              &process, 
                           MacroProblem         &Global, 
                           DataUser             &data_user) {

    TicToc2 tic;
    tic.start();
#ifdef INFO_TIME
    process.parallelisation->synchronisation();
    TicTac tic1;
    if (process.parallelisation->is_master_cpu()) tic1.start();
#endif

#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.parallelisation->rank)+" : Verifie memoire avant create_op_INTER : ").c_str(),1);
#endif

    if(process.multiresolution->type=="Off" or process.multiresolution->m==0){  
        if (process.parallelisation->is_master_cpu())
            std::cout << "Calcul des Quantites d'interfaces" << endl;
        create_op_INTER(Stot,Inter,SubI,process);
        crout << process.parallelisation->rank << " : Inter.size :" <<Inter.size() << "  :  Op inter : " ;
        tic.stop();
        process.parallelisation->synchronisation();
        tic.start();
#ifdef PRINT_ALLOC
        disp_alloc((to_string(process.parallelisation->rank)+" : Verifie memoire avant create_op_SST : ").c_str(),1);
#endif
#ifdef INFO_TIME
        process.parallelisation->synchronisation();
        if (process.parallelisation->is_master_cpu()){
            std::cout << "Operateurs d interface : " ;
            tic1.stop();
            std::cout << std::endl;
            tic1.start();
        }
#endif
    }

    if (process.parallelisation->is_master_cpu()) std::cout << "Calcul des Quantites par SST" << endl;
    create_op_SST(Stot,Inter,SubS,SubI,process, data_user);
#ifdef INFO_TIME
    process.parallelisation->synchronisation();
    if (process.parallelisation->is_master_cpu()){
        std::cout << "Creation OP SST : " ;
        tic1.stop();
        std::cout << std::endl;
        tic1.start();
    }
#endif

    crout << process.parallelisation->rank << " : Op SST : "  ;
    tic.stop();
    tic.start();

#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.parallelisation->rank)+" : Verifie memoire avant create_op_MACRO : ").c_str(),1);
#endif

    if (process.multiscale->multiechelle==1) { //cas multiechelle
        if (process.parallelisation->is_master_cpu()) std::cout << "Creation du probleme macro" << endl;
        if (process.parallelisation->is_local_cpu()){
            create_op_MACRO(Stot,Inter,process,Global);
        } else {
            create_op_MACRO(SubS,Inter,process,Global);//juste pour faire repddl pour savoir où on balance le macro dans bigF...
        }
#ifdef INFO_TIME
        process.parallelisation->synchronisation();
        if (process.parallelisation->is_master_cpu()){
            std::cout << "Creation OP MACRO : " ;
            tic1.stop();
            std::cout << std::endl;
            tic1.start();
        }
#endif
    }

    crout << process.parallelisation->rank << " : Op Macro : "  ;
    tic.stop();
    process.parallelisation->synchronisation();
    for(int i_sst = 0; i_sst < SubS.size(); i_sst++){
        //std::cout << "Matrice de raideur :" << std::endl;
        //display(*(SubS[i_sst].K));
        std
        std::cout << std::endl << "Matrice homogénéisée : (" << SubS[i_sst].LE.nb_rows() << "," << SubS[i_sst].LE.nb_cols() << ")" << std::endl;
        display(std::cout, SubS[i_sst].LE);
        std::cout << std::endl;
    }
    process.print("");
}
