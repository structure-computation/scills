//librairies Hugo
#include "containers/mat.h"
#include "containers/vecpointedvalues.h"

#ifndef INFO_TIME
#define INFO_TIME
#endif 
//fichiers de definition des variables
#include "definition_PARAM.h"
#include "definition_PARAM_MULTI.h"
#include "definition_PARAM_TEMPS.h"
#include "definition_PARAM_COMP_INTER.h"
#include "definition_GLOB.h"
#include "definition_SST_time.h"
#include "definition_INTER_time.h"

// fonction utilisees pour la creation des operateurs
#include "create_op_INTER.h"
#include "create_op_MACRO.h"
#include "create_op_SST.h"

//pour l'affichage : inclure ce .h
// #include "affichage.h"

#include "mpi.h"


using namespace LMT;
using namespace std;
extern Crout crout;

const double Apply_nb_macro::eps;

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
template <class TV1, class TV2, class TV3, class TV4>
void multiscale_operateurs(TV1 &S, TV1 &SubS,TV2 &Inter, TV3 &SubI,Param &process,  TV4 &Global) {

    TicToc2 tic;
    tic.start();
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    TicToc tic1;
    if (process.rank==0) tic1.start();
#endif

#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Verifie memoire avant create_op_INTER : ").c_str(),1);
#endif

    if (process.rank == 0)
        std::cout << "Calcul des Quantites d'interfaces" << endl;
    create_op_INTER(S,Inter,SubI,process);
    crout << process.rank << " : Inter.size :" <<Inter.size() << "  :  Op inter : " ;
    tic.stop();
    if (process.size>1)
        MPI_Barrier(MPI_COMM_WORLD);
    tic.start();
#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Verifie memoire avant create_op_SST : ").c_str(),1);
#endif
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Operateurs d interface : " ;
    if (process.rank==0) tic1.stop();;
    if (process.rank==0) tic1.start();
#endif

    if (process.rank == 0)
        std::cout << "Calcul des Quantites par SST" << endl;
    create_op_SST(S,Inter,SubS,SubI,process);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Creation OP SST : " ;
    if (process.rank==0) tic1.stop();;
    if (process.rank==0) tic1.start();
#endif

    crout << process.rank << " : Op SST : "  ;
    tic.stop();
    tic.start();

#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Verifie memoire avant create_op_MACRO : ").c_str(),1);
#endif

    if (process.multiscale->multiechelle==1) { //cas multiechelle
        if (process.rank == 0)
            std::cout << "Creation du probleme macro" << endl;
        if (process.rank == 0 and process.size>1)
            create_op_MACRO(S,Inter,process,Global);
        else
            create_op_MACRO(SubS,Inter,process,Global);//juste pour faire repddl pour savoir où on balance le macro dans bigF...
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Creation OP MACRO : " ;
    if (process.rank==0) tic1.stop();;
    if (process.rank==0) tic1.start();
#endif
    }

    crout << process.rank << " : Op Macro : "  ;
    tic.stop();
    if (process.size>1)
        MPI_Barrier(MPI_COMM_WORLD);


    if (process.rank == 0)
        std::cout << endl;
};



void fake_multiscale_operateurs() {
    Param process;

    Vec<Interface<DIM,TYPEREEL> > Inter3;
    Vec<VecPointedValues<Sst<DIM,TYPEREEL> > > SubS3,S3;
    Vec<VecPointedValues<Interface<DIM,TYPEREEL> > > SubI3;
    Glob<DIM,TYPEREEL> Global3;

    multiscale_operateurs(S3, SubS3,Inter3, SubI3, process, Global3);


}
