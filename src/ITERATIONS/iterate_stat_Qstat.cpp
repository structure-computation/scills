//librairies Hugo

#include "ITERATIONS_declarations.h"

//fichiers de definition des variables
#include "../DEFINITIONS/LatinParameters.h"
#include "../DEFINITIONS/TimeParameters.h"
#include "../DEFINITIONS/SaveParameters.h"

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


void multiscale_iterate_latin(Vec<Sst>                          &S,
                              Vec<VecPointedValues<Sst> >       &SubS, 
                              Vec<Interface>                    &Inter, 
                              Vec<VecPointedValues<Interface> > &SubI, 
                              Process                           &process, 
                              MacroProblem                      &Global, 
                              Vec<Boundary>                     &CL, 
                              DataUser                          &data_user) {
    if(process.latin->alloc_quantites==1) {
        ///1ere phase : allocations et initialisation des quantites
        if (process.rank == 0) std::cout << " Allocations des quantites d'interfaces et SST" << endl;
        //if(process.reprise_calcul!=2)  //REPRISE CALCUL NE SEMBLE PLUS UTILISE
            allocate_quantities(SubS,SubI,process,Global);//on alloue si on reprend un calcul a partir d un resultat en memoire
#ifdef PRINT_ALLOC
        disp_alloc(to_string(process.rank)+" : Memoire apres allocations : ",1);
#endif
        ///Assignation des Conditions aux limites pour chaque intervalle de temps et chaque interface de bord
        /* REPRISE CALCUL NE SEMBLE PLUS UTILISE
        if(process.reprise_calcul==1) {
            read_data_sst(S, process);
            read_data_inter(Inter, process);
            calcul_erreur_latin(SubS, Inter, process, Global);
            if (process.rank == 0)
                std::cout << "Erreur initiale apres relecture : " << process.latin->error[process.latin->iter] << endl;
        }
        //*/
        if (process.size == 1 or process.rank>0)
            assign_CL_values_space_time_latin(SubI, CL, process, data_user);
        /* REPRISE CALCUL NE SEMBLE PLUS UTILISE
        if(process.reprise_calcul>0){
            calcul_erreur_latin(SubS, Inter, process, Global);
            recopie_old_from_new(Inter,process);
            if (process.rank == 0)
                std::cout << "Erreur initiale apres assignation : " << process.latin->error[process.latin->iter] << endl;
        }
        //*/
    }
    if (process.rank == 0) std::cout << " Processus iteratif " << endl;
    iterate_latin(process,SubS,Inter,SubI,Global, data_user);
};


void multiscale_iterate_incr(Vec<Sst>                          &S,
                             Vec<VecPointedValues<Sst> >       &SubS, 
                             Vec<Interface>                    &Inter, 
                             Vec<VecPointedValues<Interface> > &SubI, 
                             Process                           &process, 
                             MacroProblem                      &Global, 
                             Vec<Boundary>                     &CL, 
                             DataUser                          &data_user, 
                             GeometryUser                      &geometry_user, 
                             FieldStructureUser                &field_structure_user) {

    process.temps->pt=1;        /// On reecrit toujours dans les memes Sst::Time et Interface::Edge::Time
    /// Présence d'interface Breakable ?
    int nb_breakable=0;
    if (process.rank == 0)
        for(unsigned q=0; q <Inter.size();q++)
            if (Inter[q].comp =="Breakable")
                nb_breakable++;
    if (process.size>1)
        MPI_Bcast(&nb_breakable,1, MPI_INT, 0, MPI_COMM_WORLD);
    process.nb_breakable = nb_breakable ;
  
    /// Allocations et initialisation des quantites
    if (process.rank == 0) std::cout << "Allocations des vecteurs de stockage des resultats" << endl;
    //if (process.reprise_calcul!=2) //REPRISE CALCUL NE SEMBLE PLUS UTILISE
        allocate_quantities(SubS,SubI,process,Global);
#ifdef PRINT_ALLOC
    disp_alloc(to_string(process.rank)+" : Memoire apres allocations : ",1);
#endif
    /* REPRISE CALCUL NE SEMBLE PLUS UTILISE
    if(process.reprise_calcul==1) {
        read_data_sst(S, process);
        read_data_inter(Inter, process);
        calcul_erreur_latin(SubS, Inter, process, Global);
    }

    if (process.reprise_calcul>0){
        calcul_erreur_latin(SubS, Inter, process, Global);
        if (process.rank == 0)
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
        if (process.rank == 0)
            std::cout << "Erreur initiale avant reprise : " << process.latin->error[process.latin->iter] << endl;
    }
    //*/
    if (process.size>1)
        MPI_Barrier(MPI_COMM_WORLD);

    /// Reactualisation du 2nd membre de micro1 pour le calcul de l'erreur
    apply_mt(SubS,process.nb_threads,Calcul_2nd_membre_micro1_sst(),process, data_user);
    
    /// Initialisation des grandeurs initiales
    std::cout << "Initialisation des grandeurs" << std::endl;
    initialise_CL_values_space_time(SubI, CL, process, data_user);

    /// Boucle sur les pas de temps
    unsigned i_step=process.temps->step_cur;
    for(unsigned i_pt = 0 ; i_pt < process.temps->time_step[i_step].nb_time_step; i_pt++){
        /// Actualisation des indicateurs temporels
        process.temps->time_step[i_step].pt_cur=i_pt;
        process.temps->pt_cur+=1;
        process.temps->current_time=process.temps->time_step[i_step].t_ini+(i_pt+1)*process.temps->time_step[i_step].dt;
        if (process.rank == 0) std::cout << std::endl << "----Piquet de temps courant " << process.temps->current_time << std::endl;

        /// Assignation des Conditions aux limites pour chaque intervalle de temps et chaque interface de bord
        if (process.rank == 0) std::cout << " - Comportements des interfaces :" << std::endl;
        if (process.size == 1 or process.rank>0)
            assign_CL_values_space_time_incr(SubI, CL, process, data_user);
        if (process.rank == 0) std::cout << std::endl << " - Conditions aux limites :" << std::endl;
        for(int ic=0;ic<CL.size();ic++){
            if (process.rank == 0) std::cout << "id : " << CL[ic].id << " : " << CL[ic].ft << std::endl;
                //std::cout <<"fspace " << CL[ic].fcts_spatiales[i_step]<< endl;
        }
        if (process.rank == 0) std::cout << std::endl;
        
        /* REPRISE CALCUL NE SEMBLE PLUS UTILISE
        if(process.reprise_calcul>0){
            calcul_erreur_incr(SubS, Inter, process, Global);
            if (process.rank == 0)
                std::cout << "Erreur initiale apres assignation : " << process.latin->error[process.latin->iter] << endl;
        }
        //*/
        
        /// Calcul sur le pas de temps
        if (nb_breakable>0) {
            int nb_change = 0;
            int sous_iter = 1;
            while(nb_change != 0 or sous_iter == 1) {
                if (process.size == 0 or process.rank > 0){
                    for(unsigned q=0; q < SubI.size();q++){
                        if (SubI[q].comp == "Breakable")
                            SubI[q].param_comp->convergence = -1; 
                    }
                }
                if (process.rank == 0) std::cout << "          Sous iteration interface cassable : " << sous_iter << std::endl;
                iterate_incr(process,SubS,Inter,SubI,Global,data_user);
                if (process.size == 0 or process.rank > 0){
                    for(unsigned q=0; q < SubI.size();q++){
                        if (SubI[q].comp == "Breakable")
                            nb_change += SubI[q].param_comp->convergence ; 
                    }
                }
            }
        } else {
            iterate_incr(process,SubS,Inter,SubI,Global,data_user);
        }
        ///assignation ptcur au ptold
        if (process.rank == 0) std::cout << " - Reactualisation des valeurs pour le pas de temps suivant" << endl;
        assign_quantities_current_to_old(SubS,SubI,process);
        
        /// Sauvegarde des resultats
        if(process.save_data==1){
            if (process.rank == 0) std::cout << " - Sauvegarde des resultats au format HDF" << std::endl; 
            if (process.size == 1 or process.rank>0) {
                // write_hdf_fields_SST_INTER(SubS, Inter, process , data_user);  BUG
                convert_fields_to_field_structure_user(SubS, Inter, process , data_user, field_structure_user, geometry_user);
                Sc2String rank; rank << process.rank;
                Sc2String file_output_hdf5 = process.affichage->name_hdf + "_" + rank + ".h5";
                field_structure_user.write_hdf5_in_parallel(file_output_hdf5, geometry_user, process.affichage->name_fields, process.temps->pt_cur, process.temps->current_time, process.rank);
            }
        }
        
        /// Modification du comportement des entites
        //modification_sst_inter_behaviour(S,Inter,process.temps);  A TESTER
        
        if (process.rank == 0) std::cout << "----Fin piquet de temps " << process.temps->current_time << std::endl;
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
// #endif
