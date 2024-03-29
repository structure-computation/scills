#ifdef qfmlqsdnkfqmdlsfnkqmelsfnkdmljcxvn // pour que rien ne soit compile
//librairies Hugo

#include "ITERATIONS_declarations.h"

//fichiers de definition des variables
#include "../DEFINITIONS/LatinData.h"
#include "../DEFINITIONS/TimeData.h"
#include "../DEFINITIONS/SavingData.h"

// fonctions speciales math
#include "../UTILITAIRES/algebre.h"
#include "../UTILITAIRES/utilitaires.h"

// fonction utilisees pour la phase d'initialisation ou affectation des valeurs des CL puis phase iterative
#include "manipulate_quantities.h"
#include "iterate.h"
#include "prelocalstage.h"
#include "ERROR/calculate_error.h"

#include "modification_sst_inter_behaviour.h"

#include "../POSTTRAITEMENTS/affichage.h"

//fonctions pour les envoies des donnees par mpi
#include "../MPI/crout.h"
#include "../POSTTRAITEMENTS/save_read_data.h"
#include "../POSTTRAITEMENTS/save_hdf_data.h"

#include "../../LMT/include/containers/mat.h"

//using namespace LMT;
//using namespace Metil;


void multiscale_iterate_incr(ScSstVec           &S,
                             ScSstRef           &SubS, 
                             ScInterVec         &Inter, 
                             ScInterRef         &SubI, 
                             Process            &process, 
                             MacroProblem       &Global, 
                             ScCLVec            &CL, 
                             DataUser           &data_user, 
                             GeometryUser       &geometry_user, 
                             FieldStructureUser &field_structure_user) {
    process.temps->pt=1;        /// On reecrit toujours dans les memes Sst::Time et Interface::Edge::Time en incremental
    /// Presence d'interface Breakable ?
    int nb_breakable=0;
    if (process.parallelisation->is_master_cpu())
        for(unsigned q=0; q <Inter.size();q++)
            if (Inter[q].comp =="Breakable")
                nb_breakable++;
            if (process.parallelisation->is_multi_cpu())
        MPI_Bcast(&nb_breakable,1, MPI_INT, 0, MPI_COMM_WORLD);
    process.nb_breakable = nb_breakable ;
    
    /// Assignation des grandeurs a l'instant initial
    if(process.temps->step_cur == 0){
        if(process.parallelisation->is_master_cpu()) std::cout << "Assignation des grandeurs a l'instant initial" << std::endl;
        initialise_CL_values(SubI, CL);
    }
    
    /* REPRISE CALCUL NE SEMBLE PLUS UTILISE
    if(process.reprise_calcul==1) {
        read_data_sst(S, process);
        read_data_inter(Inter, process);
        calcul_erreur_latin(SubS, Inter, process, Global);
    }

    if (process.reprise_calcul>0){
        calcul_erreur_latin(SubS, Inter, process, Global);
        if (process.parallelisation->is_master_cpu())
            std::cout << "Erreur initiale avant reprise : " << process.latin->error[process.latin->iter] << endl;
        for(unsigned i=0;i<Inter.size();i++){
            for(unsigned j=0;j<Inter[i].side.size();j++) {
                Inter[i].side[j].t[1]=Inter[i].side[j].t_post[1];
                Inter[i].side[j].t[0]=Inter[i].side[j].t_post[0];
                if (process.reprise_calcul==1)
                    recopie_old_from_new_post(Inter,process);
            }
        }
        calcul_erreur_incr(SubS, Inter, process, Global);
        if (process.parallelisation->is_master_cpu())
            std::cout << "Erreur initiale avant reprise : " << process.latin->error[process.latin->iter] << endl;
    }
    //*/
    process.parallelisation->synchronisation();
    
    /// Reactualisation du 2nd membre de micro1 pour le calcul de l'erreur
    //apply_mt(SubS,process.parallelisation->nb_threads,Calcul_2nd_membre_micro1_sst(),process, data_user);
    
    /// Boucle sur les pas de temps
    unsigned i_step=process.temps->step_cur;
    for(unsigned i_pt = 0 ; i_pt < process.temps->time_step[i_step].nb_time_step; i_pt++){
        /// Actualisation des indicateurs temporels
        process.temps->time_step[i_step].pt_cur=i_pt;
        process.temps->pt_cur+=1;
        process.temps->t_cur=process.temps->time_step[i_step].t_ini+(i_pt+1)*process.temps->time_step[i_step].dt;
        process.print_data("----Piquet de temps courant ",process.temps->t_cur);

        /// Assignation des Conditions aux limites pour chaque intervalle de temps et chaque interface de bord
        process.print("\n - Comportements des interfaces :");
        if (process.parallelisation->is_local_cpu()) update_CL_values(SubI, CL, process, data_user);
        process.print("\n - Conditions aux limites :");
        
        /* REPRISE CALCUL NE SEMBLE PLUS UTILISE
        if(process.reprise_calcul>0){
            calcul_erreur_incr(SubS, Inter, process, Global);
            if (process.parallelisation->is_master_cpu())
                std::cout << "Erreur initiale apres assignation : " << process.latin->error[process.latin->iter] << endl;
        }
        //*/
        /// Calcul sur le pas de temps
        if (nb_breakable>0) {
            int nb_change = 0;
            int sous_iter = 1;
            while(nb_change != 0 or sous_iter == 1) {
                if (process.parallelisation->is_local_cpu()){
                    for(unsigned q=0; q < SubI.size();q++){
                        if (SubI[q].comp == "Breakable")
                            SubI[q].convergence = -1; 
                    }
                }
                if (process.parallelisation->is_master_cpu()) std::cout << "          Sous iteration interface cassable : " << sous_iter << std::endl;
                iterate_incr(process,SubS,Inter,SubI,Global,data_user);
                if (process.parallelisation->is_local_cpu()){
                    for(unsigned q=0; q < SubI.size();q++){
                        if (SubI[q].comp == "Breakable")
                            nb_change += SubI[q].convergence ; 
                    }
                }
            }
        } else {
            iterate_incr(process,SubS,Inter,SubI,Global,data_user);
        }
        ///assignation ptcur au ptold
        process.print(" - Reactualisation des valeurs pour le pas de temps suivant");
        assign_quantities_current_to_old(SubS,SubI,process);
        
        /// Sauvegarde des resultats
        if(process.save_data==1){
            process.print(" - Sauvegarde des resultats au format HDF"); 
            if (process.parallelisation->is_local_cpu()) {
                // write_hdf_fields_SST_INTER(SubS, Inter, process , data_user);  BUG
                convert_fields_to_field_structure_user(SubS, Inter, process , data_user, field_structure_user, geometry_user);
                Sc2String rank; rank << process.parallelisation->rank;
                Sc2String file_output_hdf5 = process.affichage->name_hdf + "_" + rank + ".h5";
                field_structure_user.write_hdf5_in_parallel(file_output_hdf5, geometry_user, process.affichage->name_fields, process.temps->pt_cur, process.temps->t_cur, process.parallelisation->rank);
            }
        }
        
        /// Modification du comportement des entites
        //modification_sst_inter_behaviour(S,Inter,process.temps);  A TESTER
        
        process.print_data("----Fin piquet de temps ",process.temps->t_cur);
    }

    ///Affichage des energies
    if (process.affichage->trac_ener_imp == 1) {
        process.affichage->param_ener[0]=1; process.affichage->param_ener[1]=0;
        affichage_energie(SubS,Inter,process,data_user);
        process.affichage->param_ener[0]=1; process.affichage->param_ener[1]=1;
        affichage_energie(SubS,Inter,process,data_user);
    }
    if (process.affichage->trac_ener_diss == 1) {
        process.affichage->param_ener[0]=0; process.affichage->param_ener[1]=0;
        affichage_energie(SubS,Inter,process,data_user);
        process.affichage->param_ener[0]=0; process.affichage->param_ener[1]=1;
        affichage_energie(SubS,Inter,process,data_user);
    }

};
#endif
