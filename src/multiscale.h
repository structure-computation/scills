#include <string>

#include "read_data_process.h"
#include "read_data_structure.h"
#include "read_CL.h"

//fichiers de definition des variables
#include "definition_PARAM_COMP_INTER.h"
#include "definition_PARAM.h"
#include "definition_PARAM_LATIN.h"
#include "definition_PARAM_MULTI.h"
#include "definition_PARAM_TEMPS.h"
#include "definition_GLOB.h"
#include "definition_SST_time.h"
#include "definition_INTER_time.h"
#include "definition_CL.h"
#include "allocate.h"

//pour le post traitement
#include "containers/gnuplot.h"
#include "containers/vecpointedvalues.h"
//pour l'affichage : inclure ce .h
#include "affichage.h"
//pour l'interactivite
#include "interactivite.h"
// fcts pont entre les fichiers cpp
#include "definition_fcts.h"

#include "save_read_data.h"
#include "mpi_sendrecvall.h"

//fcts MPI
#include "calculate_error.h"
#include "prelocalstage.h"

#ifndef INFO_TIME
#define INFO_TIME
#endif 

#include "containers/evaluate_nb_cycles.h"

//biblioteque venant de SC_create_2
#include <Metil/BasicVec.h>
#include <Metil/StructCompactor.h>
#include <Metil/CudaMetil.h>
#include <Metil/Hdf.h>

#include "GeometryUser.h"
#include "DataUser.h"

using namespace Metil;
using namespace LMT;
using namespace std;



/** \ingroup
\brief Fonction principale pour un calcul sous-structuré
*/
template<class TV1,class TV2,class TV5,class GLOB>

void multiscale(DataUser &data_user, GeometryUser &geometry_user, TV1 &S, TV2 &Inter, Param &process,  TV5 &CL, GLOB &Global) {

    /// lecture des donnees de calcul
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    TicToc tic1;
    if (process.rank==0) tic1.start();
#endif
    // read_data_process(process,n);
    read_data_process(process, data_user);
    // donnees associees a la geometrie, maillage...
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Lecture data_process : " << std::endl;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) tic1.start();
#endif
    //read_data_structure(process,n);
    read_data_structure(process, data_user);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Lecture data_structure : " << std::endl;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) tic1.start();
#endif

    /// lecture des conditions aux limites
    //std::cout << "Lecture des conditions aux limites" << std::endl;
    //read_CL(n,CL,process);
    read_CL(data_user,CL,process);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Lecture des CL : " << std::endl;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) tic1.start();
#endif



    if (process.rank == 0) std::cout <<  std::endl;
    if (process.rank == 0) std::cout << "******************************************" << std::endl;
    if (process.rank == 0) std::cout << " Construction du maillage micro et macro" << std::endl;
    if (process.rank == 0) std::cout << "******************************************" << std::endl;
    /// construction de la geometrie et des maillages
    Vec<VecPointedValues<typename TV1::template SubType<0>::T> > Stot,SubS;
    Vec<VecPointedValues<typename TV2::template SubType<0>::T> > SubI;
    
    multiscale_geometry_mesh( data_user, geometry_user, S, Inter, process, CL, Stot, SubS, SubI );
    //multiscale_geometry_mesh(n,S,Inter,process,CL,Stot,SubS,SubI);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Construction maillages : " ;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) tic1.start();
#endif
    
//     //verification
//     for(int i_inter=0; i_inter<Inter.size(); i_inter++){
//         Inter[i_inter].affiche();
//     }
    
    /// assignation des proprietes materiaux aux maillages des sst

    //assignation_materials_property_SST(n, S, Inter,process);//pas de SubS, pour les directions de recherches, besoin de connaitre les E de SST pas sur le pro
    assignation_materials_property_SST(data_user, S, Inter,process);//pas de SubS, pour les directions de recherches, besoin de connaitre les E de SST pas sur le pro
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Assignation materiaux SST : " ;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) tic1.start();
#endif

    /// modification du comportement des interfaces interieures si besoin : contact...
    assignation_materials_property_INTER(data_user,Inter,S,process);//pas de SubI on verifie juste que les interfaces a modifier sont utiles pour le pro
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Assignation materiaux INTER : " ;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) tic1.start();
#endif

//     //verification
//     for(int i_inter=0; i_inter<Inter.size(); i_inter++){
//         Inter[i_inter].affiche();
//     }

    /// affichage du maillage si necessaire
    affichage_maillage(SubS,SubI,S,process);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Affichage maillage : " ;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) tic1.start();
#endif

    process.temps->dt=0;
    if (process.rank == 0)  std::cout << "Nombre de pas de temps total " << process.temps->nbpastemps << std::endl;
    process.temps->pt_cur=0;
    
    if (process.rank == 0)
        std::cout << " Allocations des quantites d'interfaces et SST" << std::endl;
    if(process.reprise_calcul!=2) allocate_quantities_post(SubS,SubI,process);
#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Memoire apres allocations : ").c_str(),1);
#endif    
    
    for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
        if (process.rank == 0)
            std::cout<<"Step : " << i_step << std::endl;
        process.temps->step_cur=i_step;
        if(process.temps->dt!=process.temps->time_step[i_step].dt){
            process.temps->dt=process.temps->time_step[i_step].dt;
            
        

            if (process.rank == 0) std::cout <<  std::endl;
            if (process.rank == 0) std::cout << "******************************************" << std::endl;
            if (process.rank == 0) std::cout << " Construction des operateurs" << std::endl;
            if (process.rank == 0) std::cout << "******************************************" << std::endl;
    
#ifdef PRINT_ALLOC
            disp_alloc((to_string(process.rank)+" : Verifie memoire avant construction : ").c_str(),1);
#endif
            if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
            /// construction des operateurs
            multiscale_operateurs(Stot,SubS,Inter,SubI,process,Global);
        }

#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Construction operateurs : " ;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) tic1.start();
#endif

            /// destruction des quantites surperflues en MPI
#ifdef PRINT_ALLOC
            disp_alloc((to_string(process.rank)+" : Verifie memoire avant free : ").c_str(),1);
#endif
            //memory_free(S,Inter,process);
#ifdef PRINT_ALLOC
            disp_alloc((to_string(process.rank)+" : Verifie memoire apres free : ").c_str(),1);
#endif
        
    if(process.read_data==1) {
      if (process.rank == 0) std::cout << "Allocations des quantites d'interfaces et SST" << std::endl;
        allocate_quantities(SubS,SubI,process,Global);
        if (process.rank == 0) std::cout << "Lecture des resultats dans les fichiers save_sst et save_inter" << std::endl;
        read_data_sst(S, process);
        read_data_inter(Inter, process);
        if (process.size > 1 ) SendRecvInterAll(process.multi_mpi->intertoexchangebypro,Inter,process);
        calcul_erreur_latin(SubS, Inter, process, Global);
        if (process.rank == 0) std::cout << "Erreur : " << process.latin->error[process.latin->iter] << std::endl;
    } else {

        if (process.rank == 0) std::cout <<  std::endl;
        if (process.rank == 0) std::cout << "**************************" << std::endl;
        if (process.rank == 0) std::cout << " DEBUT DU CALCUL ITERATIF "<< std::endl;
        if (process.rank == 0) std::cout << "**************************" << std::endl;

        /// calcul itératif
        if (process.nom_calcul=="latin") {
            if (process.rank == 0) std::cout << "Calcul latin" << std::endl;
            multiscale_iterate_latin(S,SubS,Inter, SubI,process, Global,CL);
        } else if(process.nom_calcul=="incr") {
            if (process.rank == 0) std::cout << "Calcul incremental" << std::endl;

            multiscale_iterate_incr(S,SubS, Inter,SubI,process, Global,CL);
        } else {
            if (process.rank == 0) std::cout << "nom de calcul non defini : choix entre latin ou incremental" << std::endl;
            assert(0);
        }

        if(process.save_data==1) {
          if (process.rank == 0) std::cout << "Sauvegarde des résultats dans les fichiers save_sst et save_inter" << std::endl;
            Vec<string> fields_to_save("Fchap","F","Wpchap","Wp","Wchap","W");
            save_data_inter(Inter,SubS, process, fields_to_save);
            Vec<string> fields_to_save2("q");
            save_data_sst(SubS, process, fields_to_save2);
        }
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Duree de la partie iterative : " ;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) tic1.start();
#endif
    }
    
    }

    //post traitements
    if(process.rank==0){
    /// post traitements
    // affichage sous gnuplot de la courbe d'erreur
      if (process.affichage->display_error==1) {
        GnuPlot gp;
        gp.plot(log10(process.latin->error[range(process.latin->iter)]));
        gp.wait();
      }
    }
    // affichage sous paraview du resultat
    affichage_resultats(SubS,process);
}

