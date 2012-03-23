//librairies Hugo
#include "../../LMT/include/containers/mat.h"
#include "../../LMT/include/containers/vecpointedvalues.h"

//fichiers de definition des variables
#include "../DEFINITIONS/definition_PARAM_COMP_INTER.h"
#include "../DEFINITIONS/Param.h"
#include "../DEFINITIONS/LATIN.h"
#include "../DEFINITIONS/MULTI.h"
#include "../DEFINITIONS/TEMPS.h"
#include "../DEFINITIONS/Glob.h"
#include "../DEFINITIONS/Sst.h"
#include "../DEFINITIONS/Interface.h"
#include "../DEFINITIONS/Boundary.h"
#include "../DEFINITIONS/AFFICHAGE.h"

// fonctions speciales math
#include "../UTILITAIRES/algebre.h"
#include "../UTILITAIRES/utilitaires.h"

// fonction utilisees pour la phase d'initialisation ou affectation des valeurs des CL puis phase iterative
#include "allocate.h"
#include "iterate.h"
#include "prelocalstage.h"
#include "assign_quantities_current_to_old.h"
#include "ERROR/calculate_error.h"

#include "modification_sst_inter_behaviour.h"

#include "../POSTTRAITEMENTS/affichage.h"

//fonctions pour les envoies des donn�es par mpi
#include "../MPI/crout.h"
#include "../POSTTRAITEMENTS/save_read_data.h"
#include "../POSTTRAITEMENTS/save_hdf_data.h"

#include "../COMPUTE/FieldStructureUser.h"

using namespace LMT;
using namespace Metil;


void multiscale_iterate_latin(Vec<Sst>                          &S,
                              Vec<VecPointedValues<Sst> >       &SubS, 
                              Vec<Interface>                    &Inter, 
                              Vec<VecPointedValues<Interface> > &SubI, 
                              Param                             &process, 
                              Glob                              &Global, 
                              Vec<Boundary>                     &CL, 
                              DataUser                          &data_user) {
    if(process.latin->alloc_quantites==1) {
        //1ere phase : allocations et initialisation des quantites
        if (process.rank == 0)
            std::cout << " Allocations des quantites d'interfaces et SST" << endl;
        if(process.reprise_calcul!=2) 
            allocate_quantities(SubS,SubI,process,Global);//on alloue si on reprend un calcul a partir d un resultat en memoire
#ifdef PRINT_ALLOC
        disp_alloc(to_string(process.rank)+" : Memoire apres allocations : ",1);
#endif
        //Assignation des Conditions aux limites pour chaque intervalle de temps et chaque interface de bord
        if(process.reprise_calcul==1) {
            read_data_sst(S, process);
            read_data_inter(Inter, process);
            calcul_erreur_latin(SubS, Inter, process, Global);
            if (process.rank == 0)
                std::cout << "Erreur initiale apres relecture : " << process.latin->error[process.latin->iter] << endl;
        }
        if (process.size == 1 or process.rank>0)
            assign_CL_values_space_time_latin(SubI, CL, process, data_user);
        if(process.reprise_calcul>0)
            calcul_erreur_latin(SubS, Inter, process, Global);
        if(process.reprise_calcul>0)
            recopie_old_from_new(Inter,process);
        if(process.reprise_calcul>0)
            if (process.rank == 0)
                std::cout << "Erreur initiale apres assignation : " << process.latin->error[process.latin->iter] << endl;
    }
    if (process.rank == 0)
        std::cout << " Processus iteratif " << endl;
    iterate_latin(process,SubS,Inter,SubI,Global, data_user);
};


void multiscale_iterate_incr(Vec<Sst>                          &S,
                             Vec<VecPointedValues<Sst> >       &SubS, 
                             Vec<Interface>                    &Inter, 
                             Vec<VecPointedValues<Interface> > &SubI, 
                             Param                             &process, 
                             Glob                              &Global, 
                             Vec<Boundary>                     &CL, 
                             DataUser                          &data_user, 
                             GeometryUser                      &geometry_user, 
                             FieldStructureUser                &field_structure_user) {

    //Présence d'interface Breakable ?
    int nb_breakable=0;
    if (process.rank == 0)
        for(unsigned q=0; q <Inter.size();q++)
            if (Inter[q].comp =="Breakable")
                nb_breakable++;
    if (process.size>1)
        MPI_Bcast(&nb_breakable,1, MPI_INT, 0, MPI_COMM_WORLD);
    process.nb_breakable = nb_breakable ;
  
    //1ere phase : allocations et initialisation des quantites
    if (process.rank == 0)
        std::cout << " Allocations des quantites d'interfaces et SST" << endl;
    if(process.reprise_calcul!=2) allocate_quantities(SubS,SubI,process,Global);
#ifdef PRINT_ALLOC
    disp_alloc(to_string(process.rank)+" : Memoire apres allocations : ",1);
#endif
    if(process.reprise_calcul==1) {
        read_data_sst(S, process);
        read_data_inter(Inter, process);
        calcul_erreur_latin(SubS, Inter, process, Global);
    }
    if(process.reprise_calcul>0)
      calcul_erreur_latin(SubS, Inter, process, Global);
    if(process.reprise_calcul>0)
      if (process.rank == 0)
        std::cout << "Erreur initiale avant reprise : " << process.latin->error[process.latin->iter] << endl;


    if (process.reprise_calcul>0)
        for(unsigned i=0;i<Inter.size();i++)
            for(unsigned j=0;j<Inter[i].side.size();j++) {
                Inter[i].side[j].t[1]=Inter[i].side[j].t_post[1];
                Inter[i].side[j].t[0]=Inter[i].side[j].t_post[0];
                if (process.reprise_calcul==1)
                    recopie_old_from_new_post(Inter,process);
            }

    if(process.reprise_calcul>0)
        calcul_erreur_incr(SubS, Inter, process, Global);
    if(process.reprise_calcul>0)
        if (process.rank == 0)
            std::cout << "Erreur initiale avant reprise : " << process.latin->error[process.latin->iter] << endl;
    if (process.size>1)
        MPI_Barrier(MPI_COMM_WORLD);

    // A mettre ici pour la thermique seulement
    apply_mt(SubS,process.nb_threads,calcul_secmemb_micro_sst(),process, data_user);

    if (process.rank == 0)
        std::cout<<"Processus iteratif incremental" << std::endl;


    unsigned i_step=process.temps->step_cur;
    for(unsigned i_pt = 0 ; i_pt < process.temps->time_step[i_step].nb_time_step; i_pt++){
        process.temps->time_step[i_step].pt_cur=i_pt;
        process.temps->pt_cur+=1;
        process.temps->current_time=process.temps->time_step[i_step].t_ini+(i_pt+1)*process.temps->time_step[i_step].dt;
        
        if (process.rank == 0)
            std::cout << "----Piquet de temps courant " << process.temps->current_time << std::endl;

//         if (process.size >1 ) MPI_Barrier(MPI_COMM_WORLD);
//         std::cout << process.rank<< " : barrier2 : " << endl;

        //Assignation des Conditions aux limites pour chaque intervalle de temps et chaque interface de bord
        if (process.size == 1 or process.rank>0)
            assign_CL_values_space_time_incr(SubI, CL, process, data_user);

        for(int ic=0;ic<CL.size();ic++){
            if (process.rank == 0) 
                std::cout << "ft " << CL[ic].ft << std::endl;
/*            std::cout <<"fspace " << CL[ic].fcts_spatiales[i_step]<< endl;*/
        }
        
/*        if (process.size >1 ) MPI_Barrier(MPI_COMM_WORLD);
        std::cout << process.rank<< " : barrier3 : " << endl;*/
        
        if(process.reprise_calcul>0)
            calcul_erreur_incr(SubS, Inter, process, Global);
        if(process.reprise_calcul>0)
            if (process.rank == 0)
                std::cout << "Erreur initiale apres assignation : " << process.latin->error[process.latin->iter] << endl;
        // A mettre ici dans le cas d'un comportement dependant du temps pour les sst
        //apply_mt(S,process.nb_threads,calcul_secmemb_micro_sst(),process);
            
        std::cout << "          iterate_incr : " << endl;
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
                if (process.rank == 0) std::cout << "\t\tSous iteration interface cassable : " << sous_iter << std::endl;
                    iterate_incr(process,SubS,Inter,SubI,Global);
                if (process.size == 0 or process.rank > 0){
                    for(unsigned q=0; q < SubI.size();q++){
                        if (SubI[q].comp == "Breakable")
                            nb_change += SubI[q].param_comp->convergence ; 
                    }
                }
            }
        } else {
            iterate_incr(process,SubS,Inter,SubI,Global);
        }
        //assignation ptcur au ptold
        std::cout << "          assign_quantities_current_to_old : " << endl;
        assign_quantities_current_to_old(SubS,SubI,process);
        
        if(process.save_data==1) 
            if (process.size == 1 or process.rank>0) {
                // write_hdf_fields_SST_INTER(SubS, Inter, process , data_user);
                convert_fields_to_field_structure_user(SubS, Inter, process , data_user, field_structure_user, geometry_user);
                Sc2String rank; rank << process.rank;
                Sc2String file_output_hdf5 = process.affichage->name_hdf + "_" + rank + ".h5";
                //file_output_hdf5 << process.affichage->name_hdf <<"_"<< process.rank<<".h5";
                field_structure_user.write_hdf5_in_parallel(file_output_hdf5, geometry_user, process.affichage->name_fields, process.temps->pt_cur, process.temps->current_time, process.rank);
            }
        
        //modification de certaines interfaces ou sst (exemple endommagement)
        //modification_sst_inter_behaviour(S,Inter,param_incr);
        if (process.rank == 0)
            std::cout << "---- fin piquet de temps " << process.temps->current_time << std::endl;
    }
//     calcul_erreur_latin(SubS, Inter, process, Global);
//     if (process.rank == 0)
//         std::cout << "Erreur : " << process.latin->error[process.latin->iter] << endl;

    //Affichage des energies
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
