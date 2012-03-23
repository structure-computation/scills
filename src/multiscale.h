#include "../UTILS/Sc2String.h"


//fichiers de definition des variables
#include "DEFINITIONS/definition_PARAM_COMP_INTER.h"
#include "DEFINITIONS/Param.h"
#include "DEFINITIONS/LATIN.h"
#include "DEFINITIONS/MULTI.h"
#include "DEFINITIONS/TEMPS.h"
#include "DEFINITIONS/Glob.h"
#include "DEFINITIONS/Sst.h"
#include "DEFINITIONS/Interface.h"
#include "DEFINITIONS/Boundary.h"



//pour le post traitement
#include "../LMT/include/containers/gnuplot.h"
#include "../LMT/include/containers/vecpointedvalues.h"
//pour l'affichage : inclure ce .h
#include "POSTTRAITEMENTS/affichage.h"
//pour l'interactivite
#include "POSTTRAITEMENTS/interactivite.h"
// fcts pont entre les fichiers cpp
#include "DEFINITIONS/read_CL.h"
#include "ITERATIONS/ITERATIONS_declarations.h"
#include "ITERATIONS/ERROR/calculate_error.h"
#include "MATERIAUX/MATERIAUX_declarations.h"
#include "OPERATEURS/multiscale_operateurs.h"
#include "POSTTRAITEMENTS/save_read_data.h"

//fcts MPI
#include "MPI/mpi_lmt_functions.h"
#include "MPI/mpi_sendrecvall.h"
#include "MPI/assignation_mpi.h"


#ifndef INFO_TIME
#define INFO_TIME
#endif 

// #include "containers/evaluate_nb_cycles.h"
#include "SCtime.h"

//biblioteque venant de SC_create_2
#include <Metil/BasicVec.h>
#include <Metil/StructCompactor.h>
#include <Metil/CudaMetil.h>
#include <Metil/Hdf.h>

#include "GEOMETRY/GeometryUser.h"
#include "COMPUTE/DataUser.h"
#include "COMPUTE/FieldStructureUser.h"

using namespace Metil;
using namespace LMT;

#include "write_xdmf_geometry_fields.h"
#include "save_hdf_data.h"

/** \ingroup
\brief Fonction principale pour un calcul sous-structuré. Cette routine est appelée plusieurs fois dans le cas d'une multirésolution
*/
template<class Matprops, class TV1,class TV2,class TV5,class GLOB>
void multiscale_calculation(DataUser &data_user, GeometryUser &geometry_user, Matprops& matprops, FieldStructureUser &field_structure_user, TV1 &S, TV2 &Inter, Param &process,  TV5 &CL, GLOB &Global, Vec<VecPointedValues<typename TV1::template SubType<0>::T> > &SubS,  Vec<VecPointedValues<typename TV1::template SubType<0>::T> > &Stot,  Vec<VecPointedValues<typename TV2::template SubType<0>::T> > &SubI) {
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    TicTac tic1;
    if (process.rank==0) {tic1.init();tic1.start();}
#endif
    
    //reevaluation des materiaux si utilisation de la multiresolution
    bool calculate_operator=0;
    if(data_user.options.Multiresolution_on==1 and data_user.options.Multiresolution_material_link_CL_CLvolume[0]==1){
        if (process.rank==0) std::cout << "Reevaluation (multiresolution) des materiaux SST : " ;
        assignation_materials_property_SST(data_user,matprops,S,process,field_structure_user);
        assignation_materials_property_INTER(data_user,Inter,process,field_structure_user);
        calculate_operator=1;
    }
    
    //process.temps->dt=0;
    if (process.rank == 0)  std::cout << "Nombre de pas de temps total " << process.temps->nbpastemps << std::endl;
    process.temps->pt_cur=0;
    std::cout << " PAS DE TEMPS : " << process.temps->dt << " " << process.temps->time_step[0].dt<< std::endl;
    for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
        if (process.rank == 0)
            std::cout<<"Step : " << i_step << std::endl;
        process.temps->step_cur=i_step;
        if(process.temps->dt!=process.temps->time_step[i_step].dt or calculate_operator==1){
            process.temps->dt=process.temps->time_step[i_step].dt;
            calculate_operator=0;
            
            if (process.rank == 0) std::cout <<  std::endl;
            if (process.rank == 0) std::cout << "******************************************" << std::endl;
            if (process.rank == 0) std::cout << " Construction des operateurs" << std::endl;
            if (process.rank == 0) std::cout << "******************************************" << std::endl;
    
#ifdef PRINT_ALLOC
            disp_alloc((to_string(process.rank)+" : Verifie memoire avant construction : ").c_str(),1);
#endif
            if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
            /// construction des operateurs
            multiscale_operateurs(Stot,SubS,Inter,SubI,process,Global, data_user);
#ifdef INFO_TIME
            if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
            if (process.rank==0) std::cout << "Construction operateurs : " ;
            if (process.rank==0) tic1.stop();
            if (process.rank==0) std::cout << std::endl;
            if (process.rank==0) tic1.start();
#endif            
        }

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

        /// calcul iteratif
        if (process.nom_calcul=="latin") {
            if (process.rank == 0) std::cout << "Calcul latin" << std::endl;
            multiscale_iterate_latin(S,SubS,Inter, SubI,process, Global,CL, data_user);
        } else if(process.nom_calcul=="incr") {
            if (process.rank == 0) std::cout << "Calcul incremental" << std::endl;
            multiscale_iterate_incr(S,SubS, Inter,SubI,process, Global,CL, data_user, geometry_user, field_structure_user);
        } else {
            if (process.rank == 0) std::cout << "nom de calcul non defini : choix entre latin ou incremental" << std::endl;
            assert(0);
        }


        if(process.save_data==1) {
//           if (process.rank == 0) std::cout << "Sauvegarde des resultats dans les fichiers save_sst et save_inter" << std::endl;
//             Vec<Sc2String> fields_to_save("Fchap","F","Wpchap","Wp","Wchap","W");
//             save_data_inter(Inter,SubS, process, fields_to_save);
//             Vec<Sc2String> fields_to_save2("q");
//             save_data_sst(SubS, process, fields_to_save2);
        }
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Duree de la partie iterative : " ;
    if (process.rank==0) tic1.stop(); 
    if (process.rank==0) std::cout << std::endl;
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
    
    //sortie xdmf à partir du fichier hdf5 créé au fur et à mesure du calcul
    if(process.rank==0 and process.save_data==1){
        //write_xdmf_file_compute(process, data_user);   
        
    }  
    
}



/** \ingroup
\brief Fonction principale pour un calcul sous-structurï¿½
*/
void multiscale(DataUser       &data_user,
                GeometryUser   &geometry_user,
                Vec<SstCarac>  &matprops,
                Vec<Sst>       &S,
                Vec<Interface> &Inter,
                Param          &process,
                Vec<Boundary>  &CL,
                Glob           &Global) {

    /// lecture des donnees de calcul
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    TicTac tic1;
    if (process.rank==0) {tic1.init();tic1.start();}
#endif
    // read_data_process(process,n);
    process.read_data_user(data_user);
    // donnees associees a la geometrie, maillage...
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Lecture de DataUser : " ; 
    if (process.rank==0) tic1.stop();
    if (process.rank==0) std::cout << std::endl;
    if (process.rank==0) tic1.start();
#endif

    /// lecture des conditions aux limites
    //std::cout << "Lecture des conditions aux limites" << std::endl;
    //read_CL(n,CL,process);
    read_CL(data_user,CL,process);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Lecture des CL : " ;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) std::cout << std::endl;
    if (process.rank==0) tic1.start();
#endif



    if (process.rank == 0) std::cout <<  std::endl;
    if (process.rank == 0) std::cout << "******************************************" << std::endl;
    if (process.rank == 0) std::cout << " Construction du maillage micro et macro" << std::endl;
    if (process.rank == 0) std::cout << "******************************************" << std::endl;
    /// construction de la geometrie et des maillages
    Vec<VecPointedValues<Sst> > Stot,SubS;
    Vec<VecPointedValues<Interface> > SubI;
    
    multiscale_geometry_mesh( data_user, geometry_user, S, Inter, process, CL, Stot, SubS, SubI );
    
    //multiscale_geometry_mesh(n,S,Inter,process,CL,Stot,SubS,SubI);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Construction maillages : " ;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) std::cout << std::endl;
    if (process.rank==0) tic1.start();
#endif
    
//     //verification
//     for(int i_inter=0; i_inter<Inter.size(); i_inter++){
//         Inter[i_inter].affiche();
//     }
    
    
    process.save_data=1;
    
    
    //recherche des donnees utilisant les parametres de multiresolution
    data_user.find_Multiresolution_parameters();

    FieldStructureUser field_structure_user(geometry_user);

    //assignation des proprietes materiaux aux maillages des sst
    assignation_materials_property_SST(data_user, matprops, S, process, field_structure_user);//pas de SubS, pour les directions de recherches, besoin de connaitre les E de SST pas sur le pro
    #ifdef INFO_TIME
        if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
        if (process.rank==0) std::cout << "Assignation materiaux SST : " ;
        if (process.rank==0) tic1.stop();
        if (process.rank==0) std::cout << std::endl;
        if (process.rank==0) tic1.start();
    #endif

    /// modification du comportement des interfaces interieures si besoin : contact...
    assignation_materials_property_INTER(data_user,Inter,process, field_structure_user);//pas de SubI on verifie juste que les interfaces a modifier sont utiles pour le pro
    #ifdef INFO_TIME
        if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
        if (process.rank==0) std::cout << "Assignation materiaux INTER : " ;
        if (process.rank==0) tic1.stop();
        if (process.rank==0) std::cout << std::endl;
        if (process.rank==0) tic1.start();
    #endif
            
//     //verification
//     for(int i_inter=0; i_inter<Inter.size(); i_inter++){
//         Inter[i_inter].affiche();
//     }

    /// sauvegarde du maillage pour la visualisation des resultats
    if (process.rank == 0) std::cout << "Sauvegarde de la geometrie du maillage de peau au format hdf pour la visualisation des resultats" << std::endl;
    process.affichage->name_hdf << data_user.name_directory  << "/calcul_" << data_user.id_calcul << "/results/";
    int temp=system(("mkdir -p "+process.affichage->name_hdf).c_str());//Il faut créer le répertoire results
    process.affichage->name_hdf << "geometry_fields";   
    process.affichage->name_geometry = "/Level_0/Geometry";
    process.affichage->name_fields = "/Level_0/Fields";
    if (process.size == 1 or process.rank>0) {
        if(process.save_data==1){
            // write_hdf_geometry_SST_INTER(SubS,Inter,process, geometry_user);
            Sc2String file_output_hdf ; file_output_hdf << process.affichage->name_hdf <<"_"<< process.rank<<".h5";
            geometry_user.write_hdf5_in_parallel(file_output_hdf,process.rank);
        }
    }
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    
    //ecriture du fichier de sortie xdmf
    if(process.rank==0 and process.save_data==1){
        std::cout << "Sortie xdmf" << std::endl;
        //write_xdmf_file_geometry(process, data_user);   
    }
        
    /// affichage du maillage si necessaire
    //process.affichage->affich_mesh=1;
    process.affichage->type_affichage= "all";
    affichage_maillage(SubS,SubI,S,process, data_user);
/*    process.affichage->type_affichage= "Sbord";*/
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Affichage maillage : " ;
    if (process.rank==0) std::cout << std::endl;
    if (process.rank==0) tic1.start();
#endif
    
    std::cout << std::endl;
    std::cout << "data_user.options.mode = " << data_user.options.mode << std::endl;
    if(data_user.options.mode == "test"){
       PRINT("fin de la mise en donnée, mode test, on ne va pas plus loin !");
       return;
    }
    

    
    if (process.rank == 0)
        std::cout << " Allocations des quantites d'interfaces et SST" << std::endl;
    if(process.reprise_calcul!=2) allocate_quantities_post(SubS,SubI,process);
#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Memoire apres allocations : ").c_str(),1);
#endif    
    
    
    //boucle multiresolution (avant assignation des proprietes materiaux et liaisons, les CL sont assignees plus tard)
    if(data_user.options.Multiresolution_on==1){
        if (process.rank==0) std::cout << "Calcul parametrique : " << data_user.options.Multiresolution_type <<std::endl;
        if(data_user.options.Multiresolution_type=="fatigue"){
            data_user.options.Multiresolution_nb_calcul=data_user.Multiresolution_parameters[0].nb_values;
            std::cout << "Nombre de calculs a faire : " << data_user.options.Multiresolution_nb_calcul<< std::endl;
            for(int i_res=0 ;i_res<  data_user.options.Multiresolution_nb_calcul ;i_res++){
                data_user.options.Multiresolution_current_resolution=i_res;
                for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
                    data_user.Multiresolution_parameters[i_par].current_value=data_user.options.Multiresolution_nb_cycle*i_res+data_user.Multiresolution_parameters[i_par].min_value;
                    PRINT(data_user.Multiresolution_parameters[i_par].current_value);
                }
                multiscale_calculation(data_user, geometry_user, matprops, field_structure_user, S, Inter, process,  CL, Global, SubS, Stot, SubI);   
                affichage_resultats(SubS,process, data_user);
                affichage_resultats_inter(SubI, S ,process, data_user); //sortie paraview pour les interfaces
            }
        }
        else {
            std::cerr <<"TODO : plan d'experience" << std::endl;
            assert(0);
        }
    }
    else{
        process.temps->dt=0;
        multiscale_calculation(data_user, geometry_user, matprops, field_structure_user, S, Inter, process,  CL, Global, SubS, Stot, SubI);
        affichage_resultats(SubS,process, data_user); //sortie paraview pour les sst (volume et peau)
        affichage_resultats_inter(SubI, S ,process, data_user); //sortie paraview pour les interfaces
    }
    if (process.size>1)
        MPI_Barrier(MPI_COMM_WORLD);
    
    memory_free(S,Inter,process);

}

