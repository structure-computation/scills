//librairies Hugo
#include "../../../LMT/include/containers/mat.h"
#include "../../../LMT/include/containers/evaluate_nb_cycles.h"
#include "../../../LMT/include/containers/vec_mt.h"

//fichiers de definition des variables
#include "../../DEFINITIONS/Process.h"
#include "../../DEFINITIONS/LatinData.h"
#include "../../DEFINITIONS/MultiScaleData.h"
#include "../../DEFINITIONS/ParallelisationData.h"
#include "../../DEFINITIONS/TimeData.h"

// fonctions speciales
#include "../../UTILITAIRES/utilitaires.h"

// fonctions utilisees dans le programme
#include "optimisation_direction.h"
#include "op_SST.h"

#include "mpi.h"

using namespace LMT;

/** \defgroup Operateurs_sst Creation des operateurs par sous-structure
\ingroup Operateurs
 
 Cette proc�dure permet de construire les op�rateurs relatifs aux sous-structures, � savoir la matrice de rigidit� assembl�e avec les directions de recherche et l'op�rateur homog�n�is�.
 
 Elle fait appel aux fonctions d�finies soit dans op_sst.h et dans op_sst_time.h.
 Le d�roulement des op�rations est le suivant :
 - cr�ation de la matrice de rigidit� EF : Calc_SST_rigidite_K0
 - recherche des ddls de bord et rep�rage dans la matrice EF : Calc_SST_Correspddl
 - Optimisation des directions de recherche : optimise_direction
 - Ajout des directions de recherche � la matrice de rigidit� : Calc_SST_rigidite_K0_k
 - Rep�rage des ddl macro de bord des Ssts dans l'op�rateur homog�n�is� : repere_ind_interface_LE
 - Cr�ation de l'op�rateur homog�n�is� : Calc_SST_LE
 
*/
//*******************************************************
// Creation des operateurs par SST : PROCEDURE PRINCIPALE
//******************************************************
/** \ingroup  Operateurs_sst
\brief Proc�dure principale permettant la cr�ation des op�rateurs par sous-structure.
 */
void create_op_SST(PointedSubstructures &S, VecInterfaces &Inter,PointedSubstructures &SubS,PointedInterfaces &SubI,Process &process, DataUser &data_user) {

/*#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.parallelisation->rank)+" : Verifie memoire avant Calc_SST_Correspddl : ").c_str(),1);
#endif

    if (process.parallelisation->is_master_cpu())
        std::cout << "\t Reperage des ddls de bord" << endl;
    apply_mt(SubS,process.parallelisation->nb_threads,Calc_SST_Correspddl(), process);
*/
    
#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.parallelisation->rank)+" : Verifie memoire avant optimise_direction : ").c_str(),1);
#endif

    if (process.parallelisation->is_master_cpu())
        std::cout << "\t Optimisation des directions de recherche" << endl;
    apply_mt(SubI,process.parallelisation->nb_threads,optimise_direction(),S,*process.latin);


#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.parallelisation->rank)+" : Verifie memoire avant Calc_SST_rigidite_K0_k : ").c_str(),1);
#endif

    if (process.parallelisation->is_master_cpu())
        std::cout << "\t Rigidite totale par SST" << endl;
    if (process.parallelisation->is_local_cpu())
        apply_mt(SubS,process.parallelisation->nb_threads,Calc_SST_rigidite_K0_k(), Inter,process, data_user);

    if (process.sousint == 1) apply_mt(S,process.parallelisation->nb_threads,efface_mesh_edge(),Inter);
    process.parallelisation->synchronisation();//a virer

#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.parallelisation->rank)+" : Verifie memoire avant LE : ").c_str(),1);
#endif

    if (process.multiscale->multiechelle ==1) {///cas multiechelle
        if (process.parallelisation->is_multi_cpu()) {
            for( unsigned i=0;i<S.size() ;i++ ) {
                for( unsigned k1=0;k1 <process.parallelisation->repartition_sst.size()  ;k1++ ){
                    if (find(process.parallelisation->repartition_sst[k1],LMT::_1==(int)i)) {
                        for( unsigned jj=0;jj<S[i].edge.size() ;jj++ ){
                            if (S[i].edge[jj].datanum == 0){
                                MPI_Bcast(&Inter[S[i].edge[jj].internum].nb_macro_espace,1,MPI_INT, k1 ,MPI_COMM_WORLD);
                            }
                        }
                        break;
                    }
                }
            }
        }
        process.print("\t Operateurs homogeneises");
        if (process.parallelisation->is_multi_cpu() and process.parallelisation->is_master_cpu()) {
            apply_mt(S,process.parallelisation->nb_threads,repere_ind_interface_LE(), Inter, process);
        } else {
            apply_mt(SubS,process.parallelisation->nb_threads,repere_ind_interface_LE(), Inter, process);
            apply_mt(SubS,process.parallelisation->nb_threads,Calc_SST_LE(), Inter, process);
        }
        process.parallelisation->synchronisation();//a virer

        if (process.parallelisation->is_multi_cpu()) {
            for( unsigned i=0;i<S.size() ;i++ ) {
                for( unsigned k1=0;k1 <process.parallelisation->repartition_sst.size()  ;k1++ ){
                    if (find(process.parallelisation->repartition_sst[k1],LMT::_1==(int)i)) {
                        if (process.parallelisation->rank == (int)k1) MPI_Send(&S[i].nb_macro,1,MPI_INT, 0 ,201,MPI_COMM_WORLD);
                        if (process.parallelisation->is_master_cpu()) {MPI_Status status;MPI_Recv(&S[i].nb_macro,1,MPI_INT, k1 ,201,MPI_COMM_WORLD,&status);}
                        if (process.parallelisation->is_master_cpu()) S[i].LE.resize(S[i].nb_macro);
                        if (process.parallelisation->rank == (int)k1) MPI_Send(S[i].LE.data.ptr(),S[i].LE.data.size(),MPI_DOUBLE, 0 ,201,MPI_COMM_WORLD);
                        if (process.parallelisation->is_master_cpu()) {MPI_Status status;MPI_Recv(S[i].LE.data.ptr(),S[i].LE.data.size(),MPI_DOUBLE, k1 ,201,MPI_COMM_WORLD,&status);}
                        if (process.parallelisation->rank == (int)k1) S[i].LE.free();
                        break;
                    }
                }
            }
        }
    }
};
