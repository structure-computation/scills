

#include "DEFINITIONS/Structure.h"
#include "Process.h"
#include "DEFINITIONS/Sst.h"
#include "DEFINITIONS/Interface.h"
#include "DEFINITIONS/Boundary.h"
#include "COMPUTE/DataUser.h"

#include "../LMT/include/containers/vec.h"
using namespace LMT;

#include "OPERATEURS/multiscale_operateurs.h"
#include "ITERATIONS/prelocalstage.h"
#include "ITERATIONS/LINEAR/linearstage.h"
#include "ITERATIONS/LOCAL/localstage.h"
#include "ITERATIONS/ERROR/calculate_error.h"


int Structure::algorithme_incremental(int argc, char** argv){/// Boucle temporelle
    Sc2String id_model = argv[ 1 ];
    Sc2String id_calcul = argv[ 2 ];

    initialisation_MPI(argc, argv);
    lecture_fichiers(id_model, id_calcul);
    
    if(data_user.options.mode=="visu_CL"){
        if(process.parallelisation->is_master_cpu()) std::cout << "Mode de visualisation des bords" << std::endl;
        geometry_user.visualize_group_edges_within_geometry(data_user);
    }
    else{
        geometry_user.split_group_edges_within_geometry(data_user);
        chargement_donnees();
        
        for(process.temps->init(); process.temps->not_ended(); process.temps->next()){
            /// Mis a jour des parametres dependant du temps
            //Parametre::update_time_parameters(process.temps->current_time);
            
            /// Mis a jour des operateurs (si necessaire)
            if(process.compute_operators or process.temps->step_changed()){
                multiscale_operateurs(Stot,SubS,Inter,SubI,process,Global, data_user);
                process.compute_operators = false;
            }
            
            /// Calcul des conditions aux limites
            if(process.parallelisation->is_local_cpu()){
                if(process.temps->pt_cur == 0){
                    initialise_CL_values_space_time(SubI, CL, process, data_user);
                }
                assign_CL_values_space_time_incr(SubI, CL, process, data_user);
            }
            
            /// Boucle iterative
            unsigned i_iter = 0;
            bool iteration_ended = false;
            while(not iteration_ended){
                /// Etape lineaire
                etape_lineaire(S,Inter,process,Global,data_user);
                
                /// Relaxation
                if (process.latin->iter == 0){
                    process.latin->mu=1.0;
                }else{
                    process.latin->mu=process.latin->facteur_relaxation;
                }
                apply_mt(S,process.parallelisation->nb_threads,relaxation_quantites(),Inter,process);
                
                /// Echange des grandeurs macro calculees
                if (process.parallelisation->is_multi_cpu()){
                    SendRecvInter(process.parallelisation->intertoexchangebypro,Inter,process);
                }
                
                /// Etape locale
                if (process.parallelisation->is_local_cpu()){
                    etape_locale(SubI,S,process);
                }
                
                /// Calcul de l'erreur
                calcul_erreur_incr(S,Inter,process,Global);
                
                /// Test de convergence de l'iteration
                ...
            }
        }
        
    }
    finalisation_MPI();
    if(process.parallelisation->is_master_cpu()) std::cout << "End of SC_multi_" << DIM << ".exe " << id_model << " " << id_calcul << std::endl;
    
    
}