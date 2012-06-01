#include "Structure.h"
#include "../MPI/assignation_mpi.h"

//fichiers de definition des variables
#include "definition_PARAM_COMP_INTER.h"
#include "LatinParameters.h"
#include "MultiScaleParameters.h"
#include "TimeParameters.h"

//pour l'affichage et le post-traitement
#include "../../LMT/include/containers/gnuplot.h"
#include "../POSTTRAITEMENTS/affichage.h"
#include "../POSTTRAITEMENTS/write_xdmf_geometry_fields.h"
#include "../POSTTRAITEMENTS/save_hdf_data.h"

// fcts pont entre les fichiers cpp
#include "../ITERATIONS/ITERATIONS_declarations.h"
#include "../ITERATIONS/manipulate_quantities.h"
#include "../MAILLAGE/multiscale_geometry_mesh.h"
#include "../MATERIAUX/MATERIAUX_declarations.h"
#include "../OPERATEURS/multiscale_operateurs.h"
#include "read_CL.h"

#include "../POSTTRAITEMENTS/save_read_data.h"
#include "../MPI/mpi_sendrecvall.h"
#include "../ITERATIONS/ERROR/calculate_error.h"
#include "../ITERATIONS/prelocalstage.h"

#ifndef INFO_TIME
#define INFO_TIME
#endif 

using namespace Metil;
using namespace LMT;

void Structure::synchronize_and_print_duration(TicTac& tic, const char* msg){
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << msg << " : " ; 
    if (process.rank==0) tic.stop();
    if (process.rank==0) std::cout << std::endl;
    if (process.rank==0) tic.start();
#endif
}

void Structure::initialisation_MPI(int argc,char **argv){
    process.affichage->name_data= argv[2];
    
    /// on fait ce qu il y a faire quand on est en MPI
    if (argc == 4 and strcmp(argv[3], "mpi")==0 ) {
        definition_mpi_param(process,argc,argv);
        if (process.rank == 0) std::cout << "Calcul MPI" << std::endl;
    } else {
        process.rank=0;
        process.size=1;
    }
    crout.open(process.rank);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) {tic1.init();tic1.start();}
#endif
}


void Structure::finalisation_MPI(){
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Duree complete du programme : " ;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) std::cout << std::endl;
#endif
    
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.size > 1 ) MPI_Finalize();
}

void Structure::lecture_fichiers(Sc2String id_model, Sc2String id_calcul){
    if (process.rank==0)std::cout << std::endl;
    if (process.rank==0)std::cout << "*********************************************" << std::endl;
    if (process.rank==0)std::cout << "* Lecture des données venant de SC-create_2 *" << std::endl;
    if (process.rank==0)std::cout << "*********************************************" << std::endl;
        
    if (process.rank==0)std::cout << std::endl << "Lecture du fichier de calcul-----------------------------------------------------" << std::endl;
    data_user.initialisation(id_model, id_calcul);
    //data_user.read_json_calcul();
    data_user.read_json_calcul_v2();
    
    if (process.rank==0)std::cout << std::endl << "Lecture du fichier du modele-----------------------------------------------------" << std::endl;
    geometry_user.initialisation(id_model, id_calcul);
    if (process.rank==0)std::cout << data_user.options.mode << std::endl;
    geometry_user.read_hdf5(false,true,data_user.options.mode);                       // true si on lit les info micro, true si on lit toutes les infos
    S.resize(geometry_user.nb_group_elements);
}

void Structure::chargement_donnees(){
    if (process.rank==0)std::cout << std::endl;
    if (process.rank == 0 ) std::cout << "****************************" << std::endl;
    if (process.rank == 0 ) std::cout << " Lecture des donnees " << std::endl;
    if (process.rank == 0 ) std::cout << "****************************" << std::endl;
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) {tic1.init();tic1.start();}
#endif
    /// lecture des donnees de calcul
    if (process.rank==0)std::cout << std::endl << "Lecture du DataUser -------------------------------------------------------------" << std::endl;
    process.read_data_user(data_user);
    synchronize_and_print_duration(tic1,"Lecture du DataUser");
    
    /// lecture des conditions aux limites
    if (process.rank==0)std::cout << std::endl << "Lecture des Conditions aux Limites ----------------------------------------------" << std::endl;
    read_CL(data_user,CL,process);
    synchronize_and_print_duration(tic1,"Lecture de Conditions aux Limites");
    
    
    if (process.rank == 0) std::cout <<  std::endl;
    if (process.rank == 0) std::cout << "******************************************" << std::endl;
    if (process.rank == 0) std::cout << " Construction du maillage micro et macro" << std::endl;
    if (process.rank == 0) std::cout << "******************************************" << std::endl;
    /// construction de la geometrie et des maillages
    multiscale_geometry_mesh( data_user, geometry_user, S, Inter, process, CL, Stot, SubS, SubI );
    synchronize_and_print_duration(tic1,"Construction de la geometrie et des maillages");
    
    process.save_data=1;
    
    ///recherche des donnees utilisant les parametres de multiresolution
    data_user.find_Multiresolution_parameters();
    
    field_structure_user.load_geometry_user(geometry_user);
    
    ///assignation des proprietes materiaux aux maillages des sst
    assignation_materials_property_SST(data_user, matprops, S, process, field_structure_user);//pas de SubS, pour les directions de recherches, besoin de connaitre les E de SST pas sur le pro
    synchronize_and_print_duration(tic1,"Assignation des proprietes materiaux aux maillages des sst");
    
    /// modification du comportement des interfaces interieures si besoin : contact...
    assignation_materials_property_INTER(data_user,Inter,process, field_structure_user);//pas de SubI on verifie juste que les interfaces a modifier sont utiles pour le pro
    synchronize_and_print_duration(tic1,"Modification du comportement des interfaces interieures si besoin : contact...");
}

/** \ingroup
\brief Fonction principale pour un calcul sous-structuré. Cette routine est appelée plusieurs fois dans le cas d'une multirésolution
*/
void Structure::multiscale_calculation() {
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    TicTac tic2;
    if (process.rank==0) {tic2.init();tic2.start();}
#endif
    
    bool calculate_operator=false;  /// indique si les operateurs doivent etre recalcules
    /// Reevaluation des proprietes materiaux en MultiResolution
    if(data_user.options.Multiresolution_on==1 and data_user.options.Multiresolution_material_link_CL_CLvolume[0]==1){
        if (process.rank==0) std::cout << "Reevaluation (multiresolution) des materiaux SST : ";
        assignation_materials_property_SST(data_user,matprops,S,process,field_structure_user);
        assignation_materials_property_INTER(data_user,Inter,process,field_structure_user);
        calculate_operator=true;
    }
    
    if (process.rank == 0)  std::cout << "Nombre de pas de temps total " << process.temps->nbpastemps << std::endl;
    process.temps->pt_cur=0;
    std::cout << " PAS DE TEMPS : " << process.temps->dt << " " << process.temps->time_step[0].dt<< std::endl;
    /// Boucle sur les steps temporels
    for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
        if (process.rank == 0) std::cout<<"Step : " << i_step << std::endl;
        process.temps->step_cur=i_step;
        
        /// Calcul des operateurs
        if(process.temps->dt!=process.temps->time_step[i_step].dt or calculate_operator==true){
            process.temps->dt=process.temps->time_step[i_step].dt;
            calculate_operator=false;
            
            if (process.rank == 0) std::cout <<  std::endl;
            if (process.rank == 0) std::cout << "******************************************" << std::endl;
            if (process.rank == 0) std::cout << " Calcul des operateurs" << std::endl;
            if (process.rank == 0) std::cout << "******************************************" << std::endl;
            
#ifdef PRINT_ALLOC
            disp_alloc((to_string(process.rank)+" : Verifie memoire avant construction : ").c_str(),1);
#endif
            if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);

            multiscale_operateurs(Stot,SubS,Inter,SubI,process,Global, data_user);
            synchronize_and_print_duration(tic2,"Construction operateurs");
        }
        
        /// destruction des quantites surperflues en MPI
#ifdef PRINT_ALLOC
        disp_alloc((to_string(process.rank)+" : Verifie memoire avant free : ").c_str(),1);
#endif
        //memory_free(S,Inter,process);
#ifdef PRINT_ALLOC
        disp_alloc((to_string(process.rank)+" : Verifie memoire apres free : ").c_str(),1);
#endif
        /*    STAND BY
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
        //*/
            if (process.rank == 0) std::cout <<  std::endl;
            if (process.rank == 0) std::cout << "**************************" << std::endl;
            if (process.rank == 0) std::cout << " DEBUT DU CALCUL ITERATIF "<< std::endl;
            if (process.rank == 0) std::cout << "**************************" << std::endl;
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
                //if (process.rank == 0) std::cout << "Sauvegarde des resultats dans les fichiers save_sst et save_inter" << std::endl;
                //Vec<Sc2String> fields_to_save("Fchap","F","Wpchap","Wp","Wchap","W");
                //save_data_inter(Inter,SubS, process, fields_to_save);
                //Vec<Sc2String> fields_to_save2("q");
                //save_data_sst(SubS, process, fields_to_save2);
            }
            synchronize_and_print_duration(tic2,"Duree de la partie iterative");
        //}
    }
    
    /// Affichage en temps reel de la convergence
    if(process.rank==0){
        if (process.affichage->display_error==1) {
            GnuPlot gp;
            gp.plot(log10(process.latin->error[range(process.latin->iter)]));
            gp.wait();
        }
    }
    
    ///sortie xdmf à partir du fichier hdf5 créé au fur et à mesure du calcul
    if(process.rank==0 and process.save_data==1){
        //write_xdmf_file_compute(process, data_user);
    }
    
    affichage_resultats(SubS,process, data_user);           ///sortie paraview pour les sst (volumes et peaux)
    affichage_resultats_inter(SubI, S ,process, data_user); ///sortie paraview pour les interfaces
}



/** \ingroup
\brief Fonction principale pour un calcul sous-structure.
*/
void Structure::multiscale() {
    
    /// sauvegarde du maillage pour la visualisation des resultats
    if (process.rank == 0) std::cout << "Sauvegarde de la geometrie du maillage de peau au format hdf pour la visualisation des resultats" << std::endl;
    process.affichage->name_hdf << data_user.name_directory  << "/calcul_" << data_user.id_calcul << "/results/";
    int temp=system(("mkdir -p "+process.affichage->name_hdf).c_str());//Il faut créer le répertoire results
    process.affichage->name_hdf << "geometry_fields";   
    process.affichage->name_geometry = "/Level_0/Geometry";
    process.affichage->name_fields = "/Level_0/Fields";
    if (process.size == 1 or process.rank>0) {
        if(process.save_data==1){
            //write_hdf_geometry_SST_INTER(SubS,Inter,process, geometry_user);
            Sc2String file_output_hdf ; file_output_hdf << process.affichage->name_hdf <<"_"<< process.rank<<".h5";
            geometry_user.write_hdf5_in_parallel(file_output_hdf,process.rank);
        }
    }
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    
    /// ecriture du fichier de sortie xdmf
    if(process.rank==0 and process.save_data==1){
        std::cout << "Sortie xdmf" << std::endl;
        //write_xdmf_file_geometry(process, data_user);
    }
    
    /// affichage du maillage si necessaire
    //process.affichage->affich_mesh=1;
    process.affichage->type_affichage= "all";
    affichage_maillage(SubS,SubI,S,process, data_user);
    synchronize_and_print_duration(tic1,"Affichage du maillage");
    
    /// Verification du mode de calcul ("test" ou "normal")
    std::cout << std::endl;
    std::cout << "data_user.options.mode = " << data_user.options.mode << std::endl;
    if(data_user.options.mode == "test"){
        PRINT("fin de la mise en donnée, mode test, on ne va pas plus loin !");
        return;
    }
    
    if(process.rank == 0) std::cout << " Allocations des quantites d'interfaces et SST" << std::endl;
    
    if(process.reprise_calcul!=2)
        allocate_quantities_post(SubS,SubI,process);
    
#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Memoire apres allocations : ").c_str(),1);
#endif    
    
    
    /// Lancement du calcul
    if(data_user.options.Multiresolution_on==1){
        /// Calcul avec MultiResolution
        if (process.rank==0) std::cout << "Calcul parametrique : " << data_user.options.Multiresolution_type <<std::endl;
        if(data_user.options.Multiresolution_type=="fatigue"){
            data_user.options.Multiresolution_nb_calcul=data_user.Multiresolution_parameters[0].nb_values;
            std::cout << "Nombre de calculs a faire : " << data_user.options.Multiresolution_nb_calcul<< std::endl;
            /// Boucle sur la MultiResolution
            for(int i_res=0 ;i_res<  data_user.options.Multiresolution_nb_calcul ;i_res++){
                /// Evaluation des parametres
                data_user.options.Multiresolution_current_resolution=i_res;
                for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
                    data_user.Multiresolution_parameters[i_par].current_value=data_user.options.Multiresolution_nb_cycle*i_res+data_user.Multiresolution_parameters[i_par].min_value;
                    PRINT(data_user.Multiresolution_parameters[i_par].current_value);
                }
                multiscale_calculation();
            }
        }
        else {
            std::cerr <<"TODO : plan d'experience" << std::endl;
            assert(0);
        }
    }
    else{
        process.temps->dt=0;
        multiscale_calculation();
    }
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    
    memory_free(S,Inter,process);
}