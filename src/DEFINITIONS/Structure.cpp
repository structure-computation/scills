#include "Structure.h"/*
#include "../MPI/assignation_mpi.h"

//fichiers de definition des variables
#include "definition_PARAM_COMP_INTER.h"
#include "LatinData.h"
#include "MultiScaleData.h"
#include "MultiResolutionData.h"
#include "TimeData.h"

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

void iterate_incr(Process &process,ScSstRef &SubS,ScInterVec &Inter,ScInterRef &SubI,MacroProblem &Global,DataUser &data_user);

#include "../POSTTRAITEMENTS/save_read_data.h"
#include "../MPI/mpi_sendrecvall.h"
#include "../ITERATIONS/ERROR/calculate_error.h"
#include "../ITERATIONS/prelocalstage.h"

#ifndef INFO_TIME
#define INFO_TIME
#endif 

using namespace Metil;
using namespace LMT;

void Structure::initialisation_MPI(int argc,char **argv){
    process.affichage->name_data= argv[2];
    
    /// on fait ce qu il y a faire quand on est en MPI
    if (argc == 4 and strcmp(argv[3], "mpi")==0 ) {
        definition_mpi_param(process,argc,argv);
        if (process.parallelisation->is_master_cpu()) std::cout << "Calcul MPI" << std::endl;
    } else {
        process.parallelisation->rank=0;
        process.parallelisation->size=1;
    }
    crout.open(process.parallelisation->rank);
#ifdef INFO_TIME
    process.parallelisation->synchronisation();
    if (process.parallelisation->is_master_cpu()) {tic1.init();tic1.start();}
#endif
}



void Structure::lecture_fichiers(Sc2String id_model, Sc2String id_calcul){
    process.print_title(0,Sc2String("LANCEMENT DU CALCUL ")+id_model+" "+id_calcul);
    process.print_title(1,"Lecture des donnees venant de SC-create_2");
    /// Lecture du JSON
    process.print_title(2,"Lecture du fichier de calcul");
    data_user.initialisation(id_model, id_calcul);
    //data_user.read_json_calcul();
    data_user.read_json_calcul_v2();
    
    /// Lecture du HDF5
    process.print_title(2,"Lecture du fichier du modele");
    geometry_user.initialisation(id_model, id_calcul);
    geometry_user.read_hdf5(false,true,data_user.options.mode);                       // true si on lit les info micro, true si on lit toutes les infos
}



void Structure::mise_en_donnees(){
    /// Lecture du DataUser
    process.print_title(1,"Lecture du DataUser");
#ifdef INFO_TIME
    process.parallelisation->synchronisation();
    if (process.parallelisation->is_master_cpu()) {tic1.init();tic1.start();}
#endif
    process.read_data_user(data_user);
#ifdef INFO_TIME
    process.print_duration(tic1);
#endif
    
    /// construction de la geometrie et des maillages
    process.print_title(1,"Construction du maillage micro et macro");
    multiscale_geometry_mesh( data_user, geometry_user, S, Inter, process, CL, Stot, SubS, SubI );
#ifdef INFO_TIME
    process.print_duration(tic1);
#endif
    
    /// Remplissage du FieldStructureUser
    /// Chargement des maillages
    field_structure_user.load_geometry_user(geometry_user);      
    /// Assignation des proprietes aux group_elements (pour GPU)
    field_structure_user.assign_material_id_to_group_elements(data_user);
    field_structure_user.assign_material_properties_to_group_elements(data_user,mat_prop_temp);
    ///idem pour les group_interfaces (pour GPU)
    field_structure_user.assign_link_id_to_group_interfaces(data_user);
    field_structure_user.assign_link_properties_to_group_interfaces(data_user,link_prop_temp);
    
    /// sauvegarde du maillage pour la visualisation des resultats
    process.print("Sauvegarde de la geometrie du maillage de peau au format hdf pour la visualisation des resultats");
    process.affichage->name_hdf << data_user.name_directory  << "/calcul_" << data_user.id_calcul << "/results/";
    int temp=system(("mkdir -p "+process.affichage->name_hdf).c_str());//Il faut creer le repertoire results
    process.affichage->name_hdf << "geometry_fields";   
    process.affichage->name_geometry = "/Level_0/Geometry";
    process.affichage->name_fields = "/Level_0/Fields";
    if (process.parallelisation->is_local_cpu()) {
        if(process.save_data==1){
            //write_hdf_geometry_SST_INTER(SubS,Inter,process, geometry_user);
            Sc2String file_output_hdf ; file_output_hdf << process.affichage->name_hdf <<"_"<< process.parallelisation->rank<<".h5";
            geometry_user.write_hdf5_in_parallel(file_output_hdf,process.parallelisation->rank);
        }
    }
    process.parallelisation->synchronisation();
    
    /// ecriture du fichier de sortie xdmf
    if(process.parallelisation->is_master_cpu() and process.save_data==1){
        process.print("Sortie XDMF");
        //write_xdmf_file_geometry(process, data_user);
    }
    
    /// affichage du maillage si necessaire
    //process.affichage->affich_mesh=1;
    process.affichage->type_affichage= "all";
    affichage_maillage(SubS,SubI,S,process, data_user);
    #ifdef INFO_TIME
    process.print_duration(tic1);
    #endif
    
    /// Verification du mode de calcul ("test" ou "normal")
    if(data_user.options.mode == "test"){
        process.print("FIN DU MODE TEST.");
        return;
    }
    
    /// Allocations et initialisation des quantites
    process.print("Allocations des vecteurs de stockage des resultats");
    allocate_quantities(SubS,SubI,process,Global);
    #ifdef PRINT_ALLOC
    disp_alloc(to_string(process.parallelisation->rank)+" : Memoire apres allocations : ",1);
    #endif
}



void Structure::boucle_multi_resolution() {
    /// Lancement du calcul
    MultiResolutionData* multiresolution = process.multiresolution;
    process.print_data("Calcul parametrique : ",multiresolution->type);
    process.print_data("Nombre de calculs : ",multiresolution->nb_calculs);
    for(multiresolution->init();multiresolution->has_next();multiresolution->next()){
        multiresolution->update_parameters();                                           /// Mise a jour des parametres de multi-resolution
        SstCarac::sst_materials_parameters.update_all(multiresolution->m.value);        /// Mise a jour des parametres materiaux des sst
        InterCarac::inter_materials_parameters.update_all(multiresolution->m.value);    /// Mise a jour des parametres materiaux des interfaces
        boucle_temporelle();
    }
    process.parallelisation->synchronisation();
    
    memory_free(S,Inter,process);
}



void Structure::boucle_temporelle(){
    #ifdef INFO_TIME
    process.parallelisation->synchronisation();
    TicTac tic2;
    if (process.parallelisation->is_master_cpu()) {tic2.init();tic2.start();}
    #endif
    
    /// Boucle sur les steps temporels
    TimeData* temps = process.temps;
    process.print_title(1,"DEBUT DU CALCUL ITERATIF ");
    process.print_data("Nombre de pas de temps total : ",temps->nbpastemps);
    for(temps->init();temps->has_next();temps->next()){
        if(temps->step_changed()){
            process.print_data("********************Step : ",temps->step_cur);
        }
        process.print_data("----Piquet de temps courant ",temps->t_cur);
        temps->update_parameters();
        
        /// Calcul des operateurs
        process.print_title(1,"Mise a jour des operateurs");
        for(int i_sst = 0; i_sst < S.size(); i_sst++){
            if(S[i_sst].update_operator){
                process.print_data("Solide ",S[i_sst].id);
                #ifdef PRINT_ALLOC
                disp_alloc((to_string(process.parallelisation->rank)+" : Verifie memoire avant construction : ").c_str(),1);
                #endif
                multiscale_operateurs(Stot,SubS,Inter,SubI,process,Global, data_user);
                #ifdef PRINT_ALLOC
                disp_alloc((to_string(process.parallelisation->rank)+" : Verifie memoire apres construction : ").c_str(),1);
                #endif
                break;  /// multiscale_operateurs a remis a jour les operateurs de tout le monde
            }
        }
        #ifdef INFO_TIME
        process.print_duration(tic2);
        #endif
        
        if(process.nom_calcul=="incr") {
            process.print("Calcul incremental");
            /// Presence d'interface Breakable ?
            int nb_breakable=0;
            if (process.parallelisation->is_master_cpu()){
                for(unsigned q=0; q <Inter.size();q++){
                    if (Inter[q].comp =="Breakable"){
                        nb_breakable++;
                    }
                    if (process.parallelisation->is_multi_cpu()){
                        MPI_Bcast(&nb_breakable,1, MPI_INT, 0, MPI_COMM_WORLD);
                    }
                    process.nb_breakable = nb_breakable ;
                }
            }
            
            /// Mise a jour des conditions aux limites
            if(process.temps->pt_cur == 0 and process.parallelisation->is_local_cpu()){
                process.print(" - Initialisation des Conditions aux limites :");
                initialise_CL_values(SubI, CL);
            }
            process.print("\n - Mise a jour des Conditions aux limites :");
            if (process.parallelisation->is_local_cpu()) update_CL_values(SubI, CL, process, data_user);
        
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
            if(process.save_data){
                process.print(" - Sauvegarde des resultats au format HDF"); 
                if (process.parallelisation->is_local_cpu()) {
                    // write_hdf_fields_SST_INTER(SubS, Inter, process , data_user);  BUG
                    convert_fields_to_field_structure_user(SubS, Inter, process , data_user, field_structure_user, geometry_user);
                    Sc2String rank; rank << process.parallelisation->rank;
                    Sc2String file_output_hdf5 = process.affichage->name_hdf + "_" + rank + ".h5";
                    field_structure_user.write_hdf5_in_parallel(file_output_hdf5, geometry_user, process.affichage->name_fields, temps->pt_cur, temps->t_cur, process.parallelisation->rank);
                }
            }
            
            /// Modification du comportement des entites
            //modification_sst_inter_behaviour(S,Inter,temps);  A TESTER
            
            process.print_data("----Fin piquet de temps ",temps->t_cur);
        
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
        } else {
            std::cerr << "Nom de calcul non defini : incremental uniquement" << std::endl;
            assert(0);
        }
    }

    
    #ifdef INFO_TIME
    process.print_duration(tic2);
    #endif
    
    ///sortie xdmf à partir du fichier hdf5 cree au fur et à mesure du calcul
    if(process.parallelisation->is_master_cpu() and process.save_data==1){
        //write_xdmf_file_compute(process, data_user);
    }
    
    affichage_resultats(SubS,process, data_user);           ///sortie paraview pour les sst (volumes et peaux)
    affichage_resultats_inter(SubI, S ,process, data_user); ///sortie paraview pour les interfaces
}



void Structure::finalisation_MPI(){
    process.print("Programme complet : ",true);
#ifdef INFO_TIME
    process.print_duration(tic1);
#endif
    if (process.parallelisation->is_multi_cpu()) MPI_Finalize();
}
*/
