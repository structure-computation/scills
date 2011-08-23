//librairies Hugo
#include "containers/mat.h"
#include "containers/evaluate_nb_cycles.h"
#include "containers/vec_mt.h"

//fichiers de definition des variables
#include "definition_PARAM.h"
#include "definition_PARAM_LATIN.h"
#include "definition_PARAM_MULTI.h"
#include "definition_PARAM_MPI.h"

// fonctions speciales
#include "utilitaires.h"

// fonctions utilisees dans le programme
#include "definition_PARAM_TEMPS.h"
#include "optimisation_direction.h"
#include "op_SST.h"

#include "mpi.h"

using namespace LMT;
using namespace std;

/** \defgroup Operateurs_sst Creation des operateurs par sous-structure
\ingroup Operateurs
 
 Cette procédure permet de construire les opérateurs relatifs aux sous-structures, à savoir la matrice de rigidité assemblée avec les directions de recherche et l'opérateur homogénéisé.
 
 Elle fait appel aux fonctions définies soit dans op_sst.h et dans op_sst_time.h.
 Le déroulement des opérations est le suivant :
 - création de la matrice de rigidité EF : Calc_SST_rigidite_K0
 - recherche des ddls de bord et repérage dans la matrice EF : Calc_SST_Correspddl
 - Optimisation des directions de recherche : optimise_direction
 - Ajout des directions de recherche à la matrice de rigidité : Calc_SST_rigidite_K0_k
 - Repérage des ddl macro de bord des Ssts dans l'opérateur homogénéisé : repere_ind_interface_LE
 - Création de l'opérateur homogénéisé : Calc_SST_LE
 
*/
//*******************************************************
// Creation des operateurs par SST : PROCEDURE PRINCIPALE
//******************************************************
/** \ingroup  Operateurs_sst
\brief Procédure principale permettant la création des opérateurs par sous-structure.
 */
template<class TV1, class TV2,class TV3, class TV4>
void create_op_SST(TV1 &S, TV2 &Inter,TV3 &SubS,TV4 &SubI,Param &process, DataUser &data_user) {

/*#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Verifie memoire avant Calc_SST_Correspddl : ").c_str(),1);
#endif

    if (process.rank == 0)
        std::cout << "\t Reperage des ddls de bord" << endl;
    apply_mt(SubS,process.nb_threads,Calc_SST_Correspddl(), process);
*/
    
#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Verifie memoire avant optimise_direction : ").c_str(),1);
#endif

    if (process.rank == 0)
        std::cout << "\t Optimisation des directions de recherche" << endl;
    apply_mt(SubI,process.nb_threads,optimise_direction(),S,*process.latin);


#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Verifie memoire avant Calc_SST_rigidite_K0_k : ").c_str(),1);
#endif

    if (process.rank == 0)
        std::cout << "\t Rigidite totale par SST" << endl;
    if (process.size == 1 or process.rank>0)
        apply_mt(SubS,process.nb_threads,Calc_SST_rigidite_K0_k(), Inter,process, data_user);

    if (process.sousint == 1) apply_mt(S,process.nb_threads,efface_mesh_edge(),Inter);
    if (process.size > 1 )
        MPI_Barrier(MPI_COMM_WORLD);//a virer

#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Verifie memoire avant LE : ").c_str(),1);
#endif

    if (process.multiscale->multiechelle ==1) {//cas multiechelle
    if (process.size>1) {
        for( unsigned i=0;i<S.size() ;i++ ) {
            for( unsigned k1=0;k1 <process.multi_mpi->repartition_sst.size()  ;k1++ )
                if (find(process.multi_mpi->repartition_sst[k1],LMT::_1==(int)i)) {
                for( unsigned jj=0;jj<S[i].edge.size() ;jj++ ){
                    if (S[i].edge[jj].datanum == 0){
                        MPI_Bcast(&Inter[S[i].edge[jj].internum].nb_macro_espace,1,MPI_INT, k1 ,MPI_COMM_WORLD);
                    }
                }

                break;
                }
        }
    }
    if (process.rank == 0)
            std::cout << "\t Operateurs homogeneises" << endl;
        if (process.size > 1 and process.rank==0 ) {
            apply_mt(S,process.nb_threads,repere_ind_interface_LE(), Inter, process);
        } else {
            apply_mt(SubS,process.nb_threads,repere_ind_interface_LE(), Inter, process);
            apply_mt(SubS,process.nb_threads,Calc_SST_LE(), Inter, process);
        }
        if (process.size > 1 )
            MPI_Barrier(MPI_COMM_WORLD);//a virer

        if (process.size>1) {
            for( unsigned i=0;i<S.size() ;i++ ) {
                for( unsigned k1=0;k1 <process.multi_mpi->repartition_sst.size()  ;k1++ )
                        if (find(process.multi_mpi->repartition_sst[k1],LMT::_1==(int)i)) {
                        if (process.rank == (int)k1 ) MPI_Send(&S[i].nb_macro,1,MPI_INT, 0 ,201,MPI_COMM_WORLD);
                        if (process.rank == 0 ) {MPI_Status status;MPI_Recv(&S[i].nb_macro,1,MPI_INT, k1 ,201,MPI_COMM_WORLD,&status);}
                        if (process.rank == 0 ) S[i].LE.resize(S[i].nb_macro);
                        if (process.rank == (int)k1 ) MPI_Send(S[i].LE.data.ptr(),S[i].LE.data.size(),MPI_DOUBLE, 0 ,201,MPI_COMM_WORLD);
                        if (process.rank == 0 ) {MPI_Status status;MPI_Recv(S[i].LE.data.ptr(),S[i].LE.data.size(),MPI_DOUBLE, k1 ,201,MPI_COMM_WORLD,&status);}
                        if (process.rank == (int)k1 ) S[i].LE.free();
                        break;
                    }
            }
        }
    }
};
