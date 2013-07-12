#include "Process.h"



#include <Metil/BasicVec.h>
#include <Metil/StructCompactor.h>
#include <Metil/CudaMetil.h>
#include <Metil/Hdf.h>

#include "MPI/assignation_mpi.h"
#include "MAILLAGE/multiscale_geometry_mesh.h"
#include "MATERIAUX/assignation_material_properties_Sst.h"
#include "OPERATEURS/multiscale_operateurs.h"
#include "ITERATIONS/prelocalstage.h"
#include "ITERATIONS/iterate.h"
#include "POSTTRAITEMENTS/affichage.h"
#include "POSTTRAITEMENTS/calculs_resultantes.h"
#include "POSTTRAITEMENTS/save_hdf_data.h"


Process::Process()
{
    // initialisation des valeurs
    sousint=1;
    type_sousint="p";
    nom_calcul="latin";
    recopie_t_post=0;
    save_data=1;
    allocate();
}



Process::~Process(){
    free();
}

void Process::allocate(){
    affichage               = new SavingData;
    structure               = new GeneralData;
    latin                   = new LatinData;
    multiscale              = new MultiScaleData();
    temps                   = new TimeData;
    //properties              = new PROPERTY;
    parallelisation         = new ParallelisationData;
    multiresolution         = new MultiResolutionData;
    data_user               = new DataUser;
    geometry_user           = new GeometryUser;
    field_structure_user    = new FieldStructureUser;
    Global                  = new MacroProblem;
    S                       = new Substructures;
    Inter                   = new VecInterfaces;
    SubS                    = new PointedSubstructures;
    Stot                    = new PointedSubstructures;
    SubI                    = new PointedInterfaces;
    CL                      = new Boundaries;
    Fvol                    = new VolumicForces;
    Tload                   = new ThermalLoad;
    sst_materials           = new Materials;
    inter_materials         = new Links;
    #ifdef PRINT_ALLOC
    total_allocated[ typeid(SavingData).name() ]          += sizeof(SavingData);
    total_allocated[ typeid(GeneralData).name() ]         += sizeof(GeneralData);
    total_allocated[ typeid(LatinData).name() ]           += sizeof(LatinData);
    total_allocated[ typeid(MultiScaleData).name() ]      += sizeof(MultiScaleData);
    total_allocated[ typeid(TimeData).name() ]            += sizeof(TimeData);
    //total_allocated[ typeid(PROPERTY).name() ]            += sizeof(PROPERTY);
    total_allocated[ typeid(ParallelisationData).name() ] += sizeof(ParallelisationData);
    #endif
}

void Process::free(){
    if (affichage               != NULL) delete affichage;
    if (structure               != NULL) delete structure;
    if (latin                   != NULL) delete latin;
    if (multiscale              != NULL) delete multiscale;
    if (temps                   != NULL) delete temps;
    //if (properties              != NULL) delete properties;
    if (parallelisation         != NULL) delete parallelisation;
    if (data_user               != NULL) delete data_user;
    if (geometry_user           != NULL) delete geometry_user;
    if (field_structure_user    != NULL) delete field_structure_user;
    if (Global                  != NULL) delete Global;
    if (S                       != NULL) delete S;
    if (Inter                   != NULL) delete Inter;
    if (SubS                    != NULL) delete SubS;
    if (Stot                    != NULL) delete Stot;
    if (SubI                    != NULL) delete SubI;
    if (CL                      != NULL) delete CL;
    if (Fvol                    != NULL) delete Fvol;
    if (Tload                   != NULL) delete Tload;
    if (sst_materials           != NULL) delete sst_materials;
    if (inter_materials         != NULL) delete inter_materials;
}

void Process::initialisation_MPI_for_scwal(int argc,char **argv){
    affichage->name_data= "test";

    /// on fait ce qu il y a faire quand on est en MPI
    if (argc == 3 and strcmp(argv[2], "mpi")==0 ) {
        definition_mpi_param(*this,argc,argv);
        if (parallelisation->is_master_cpu()) std::cout << "Calcul MPI" << std::endl;
    } else {
        parallelisation->rank=0;
        parallelisation->size=1;
    }
    crout.open(parallelisation->rank);
    #ifdef INFO_TIME
    parallelisation->synchronisation();
    if (parallelisation->is_master_cpu()) {tic1.init();tic1.start();}
    #endif
}

void Process::initialisation_MPI(int argc,char **argv){
    affichage->name_data= argv[2];
    
    /// on fait ce qu il y a faire quand on est en MPI
    if (argc == 4 and strcmp(argv[3], "mpi")==0 ) {
        definition_mpi_param(*this,argc,argv);
        if (parallelisation->is_master_cpu()) std::cout << "Calcul MPI" << std::endl;
    } else {
        parallelisation->rank=0;
        parallelisation->size=1;
    }
    crout.open(parallelisation->rank);
    #ifdef INFO_TIME
    parallelisation->synchronisation();
    if (parallelisation->is_master_cpu()) {tic1.init();tic1.start();}
    #endif
}

void Process::test_MPI(){
    // ceci est un test ----------
    //sleep(2);
    parallelisation->synchronisation();
    if(parallelisation->size > 1){ 
      std::cout << "rank = " << parallelisation->rank << " ; size = " << parallelisation->size << endl; 
    }
    //sleep(2);
    parallelisation->synchronisation();
    Vec<double> test;
    test.resize(10,0.0);
    if(parallelisation->is_master_cpu()) test.set(1.0);
    
    PRINT(parallelisation->rank);
    PRINT(test);
    std::cout << endl;
    
    if(parallelisation->is_master_cpu()) sleep(5);
    //if (process.parallelisation->size > 1) MPI_Barrier( MPI_COMM_WORLD );
    if (parallelisation->size > 1) MPI_Bcast(test.ptr(),test.size() , MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    PRINT(parallelisation->rank);
    PRINT(test);
    std::cout << endl;
    
    //sleep(2);
    parallelisation->synchronisation();
    if(parallelisation->size > 1){ 
      std::cout << "rank = " << parallelisation->rank << " ; size = " << parallelisation->size << endl;
    }
    //sleep(5);
    parallelisation->synchronisation();
    // ceci est un test ----------
}


void Process::lecture_fichiers(Sc2String id_model, Sc2String id_calcul){
    print_title(0,Sc2String("LANCEMENT DU CALCUL ")+id_model+" "+id_calcul);
    print_title(1,"Lecture des donnees venant de SC-create_2");
    /// Lecture du JSON
    print_title(2,"Lecture du fichier de calcul");
    data_user->initialisation(id_model, id_calcul);
    data_user->read_json_calcul();
    
    /// Lecture du HDF5
    print_title(2,"Lecture du fichier du modele");
    geometry_user->initialisation(id_model, id_calcul);
    geometry_user->read_hdf5(false,true,data_user->options.mode);                       // true si on lit les info micro, true si on lit toutes les infos
}


void Process::test_load_data(){
    /// Lecture du DataUser
    print_title(1,"Lecture du DataUser");
    #ifdef INFO_TIME
    parallelisation->synchronisation();
    if (parallelisation->is_master_cpu()) {tic1.init();tic1.start();}
    #endif
    read_data_user();
    #ifdef INFO_TIME
    print_duration(tic1);
    #endif
    
    /// construction de la geometrie et des maillages + repartition MPI
    print_title(1,"Construction du maillage micro et macro");
    multiscale_geometry_mesh( *data_user, *geometry_user, *S, *Inter, *this, *CL, *Stot, *SubS, *SubI );
    #ifdef INFO_TIME
    print_duration(tic1);
    #endif
    
    /// Remplissage du FieldStructureUser  A REVOIR
    /// Chargement des maillages
    //field_structure_user->load_geometry_user(geometry_user);      
    /// Assignation des proprietes aux group_elements (pour GPU)
    //field_structure_user->assign_material_properties_to_group_elements(data_user,mat_prop_temp);
    /// Idem pour les group_interfaces (pour GPU)
    //field_structure_user->assign_link_properties_to_group_interfaces(data_user,link_prop_temp);
    
    /// Sauvegarde du maillage pour la visualisation des resultats
//     print("Sauvegarde de la geometrie du maillage de peau au format hdf pour la visualisation des resultats");
//     affichage->name_hdf = data_user->result_path;
//     int temp=system(("mkdir -p "+affichage->name_hdf).c_str());//Il faut creer le repertoire results
//     affichage->name_hdf << "geometry_fields";
//     affichage->name_geometry = "/Level_0/Geometry";
//     affichage->name_fields = "/Level_0/Fields";
//     if (parallelisation->is_local_cpu()) {
//         write_hdf_geometry_SST_INTER(*SubS,*Inter,*this, *geometry_user);
//         Sc2String file_output_hdf ; file_output_hdf << affichage->name_hdf <<"_"<< parallelisation->rank<<".h5";
//         geometry_user->write_hdf5_in_parallel(file_output_hdf,parallelisation->rank);
//     }
    parallelisation->synchronisation();
    
    /// ecriture du fichier de sortie xdmf
    if(parallelisation->is_master_cpu() and save_data==1){
        print("Sortie XDMF");
        //write_xdmf_file_geometry(*this, data_user);
    }
    
    /// affichage du maillage si necessaire
    //affichage->affich_mesh=1;
    affichage->type_affichage= "all";
    affichage->affich_mesh= 1;
    affichage_maillage(*SubS,*SubI,*S,*this, *data_user);
    #ifdef INFO_TIME
    print_duration(tic1);
    #endif
}


void Process::preparation_calcul(){
    /// Lecture du DataUser
    print_title(1,"Lecture du DataUser");
    #ifdef INFO_TIME
    parallelisation->synchronisation();
    if (parallelisation->is_master_cpu()) {tic1.init();tic1.start();}
    #endif
    read_data_user();
    #ifdef INFO_TIME
    print_duration(tic1);
    #endif
    
    /// construction de la geometrie et des maillages + repartition MPI
    print_title(1,"Construction du maillage micro et macro");
    multiscale_geometry_mesh( *data_user, *geometry_user, *S, *Inter, *this, *CL, *Stot, *SubS, *SubI );
    #ifdef INFO_TIME
    print_duration(tic1);
    #endif
    
    /// Remplissage du FieldStructureUser  A REVOIR
    /// Chargement des maillages
    //field_structure_user->load_geometry_user(geometry_user);      
    /// Assignation des proprietes aux group_elements (pour GPU)
    //field_structure_user->assign_material_properties_to_group_elements(data_user,mat_prop_temp);
    /// Idem pour les group_interfaces (pour GPU)
    //field_structure_user->assign_link_properties_to_group_interfaces(data_user,link_prop_temp);
    
    /// Sauvegarde du maillage pour la visualisation des resultats
//     print("Sauvegarde de la geometrie du maillage de peau au format hdf pour la visualisation des resultats");
//     affichage->name_hdf = data_user->result_path;
//     int temp=system(("mkdir -p "+affichage->name_hdf).c_str());//Il faut creer le repertoire results
//     affichage->name_hdf << "geometry_fields";
//     affichage->name_geometry = "/Level_0/Geometry";
//     affichage->name_fields = "/Level_0/Fields";
//     if (parallelisation->is_local_cpu()) {
//         if(save_data==1){
//             write_hdf_geometry_SST_INTER(*SubS,*Inter,*this, *geometry_user);
//             Sc2String file_output_hdf;
//             file_output_hdf << affichage->name_hdf <<"_"<< parallelisation->rank<<".h5";
//             geometry_user->write_hdf5_in_parallel(file_output_hdf,parallelisation->rank);
//         }
//     }
    parallelisation->synchronisation();
    
    /// ecriture du fichier de sortie xdmf
//     if(parallelisation->is_master_cpu() and save_data==1){
//         print("Sortie XDMF");
//         write_xdmf_file_geometry(*this, data_user);
//     }
    
    /// Creation des liens vers les materiaux et les formulations
    print("assignation des comportements matériaux");
    apply(*S,assignation_material_to_SST(),*sst_materials,plasticite,endommagement);
    
    print("assignation des comportements liaisons");
    for(unsigned i = 0; i < Inter->size(); i++){
        PRINT((*Inter)[i].id_link);
        if((*Inter)[i].id_link >= 0){
            int index_link = (*data_user).find_links_index((*Inter)[i].id_link);
            (*Inter)[i].matprop = &((*inter_materials)[index_link]);
            PRINT((*Inter)[i].matprop->comp);
            PRINT((*Inter)[i].matprop->type_num);
            /// Assignation du type de comportement
            switch((*Inter)[i].matprop->type_num){
                case 0:
                    (*Inter)[i].comp = Interface::comp_parfait;
                    break;
                case 1:
                    (*Inter)[i].comp = Interface::comp_elastique;
                    break;
                case 2:
                    (*Inter)[i].comp = Interface::comp_contact_parfait;
                    break;
                case 3:
                    (*Inter)[i].comp = Interface::comp_cassable_parfait;
                    break;
                case 4:
                    (*Inter)[i].comp = Interface::comp_cassable_elastique;
                    break;
                case 5:
                    (*Inter)[i].comp = Interface::comp_cohesive;
                    break;
                default:
                    break;
            }
        }
    }
    
    for(unsigned i = 0; i < Inter->size(); i++){
        (*Inter)[i].affiche();
    }
    
    /// affichage du maillage si necessaire
    //affichage->affich_mesh=1;
    print("affichage du maillage");
    affichage->type_affichage= "all";
    affichage_maillage(*SubS,*SubI,*S,*this, *data_user);
    ///creation des fichiers pvd
    create_pvd_geometry(*SubS,*S,*Inter,*this);
    
    #ifdef INFO_TIME
    print_duration(tic1);
    #endif
    
    
    /// Verification du mode de calcul ("test" ou "normal")
    if(data_user->options.mode == "test"){
        print("FIN DU MODE TEST.");
        assert(0);
    }
    /*
    /// Allocations et initialisation des quantites
    print("Allocations des vecteurs de stockage des resultats");
    allocate_quantities_Sst_Inter(*SubS,*SubI,*this);
    allocate_quantities_post(*SubS,*SubI,*this);
    #ifdef PRINT_ALLOC
    disp_alloc(to_string(parallelisation->rank)+" : Memoire apres allocations : ",1);
    #endif
    print("Repartition des solides");
    if(parallelisation->is_master_cpu()){
        for(int i = 0; i < parallelisation->repartition_sst.size(); i++){
            std::cout << "\t" << i << " :";
            unsigned nb_nodes = 0;
            for(int j = 0; j < parallelisation->repartition_sst[i].size(); j++){
                std::cout << "\t" << parallelisation->repartition_sst[i][j];
                nb_nodes += (*S)[parallelisation->repartition_sst[i][j]].mesh.elem_list_size;
            }
            std::cout << "\t" << nb_nodes << std::endl;
        }
    }
    //*/
}


// template<class Updater_, class MP_>
// void Process::boucle_multi_resolution(Updater_ updater, MP_ mp) {
void Process::boucle_multi_resolution() {
    /// Lancement des calculs parametriques
//     updater.add_message( mp, updater.ET_Info, "DEBUT DES CALCULS PARAMETRIQUES" );
    print_title(1,"DEBUT DES CALCULS PARAMETRIQUES ");
    print_data("Calcul parametrique : ",multiresolution->type);
    print_data("Nombre de calculs : ",multiresolution->nb_calculs);
    for(multiresolution->init();multiresolution->has_next();multiresolution->next()){
        print_data("************************************************************ Calcul : ",multiresolution->calcul_cur);
        multiresolution->updateParameters();    /// Mise a jour des parametres de multi-resolution
        SstCarac::updateParameters();           /// Mise a jour des parametres materiaux des sst
        boucle_temporelle();
        print_data("******************************************************** Fin Calcul : ",multiresolution->calcul_cur);
        parallelisation->synchronisation();
    }
    
    
    memory_free(*S,*Inter,*CL,*sst_materials,*inter_materials,*this);
    PRINT("fin de la multiresolution");
}



void Process::boucle_temporelle(){
    #ifdef INFO_TIME
    parallelisation->synchronisation();
    TicTac tic2;
    if (parallelisation->is_master_cpu()) {tic2.init();tic2.start();}
    #endif
    
    /// Allocations et initialisation des quantites
    print("Allocations des vecteurs de stockage des resultats");
    allocate_quantities_Sst_Inter(*SubS,*SubI,*this);
    allocate_quantities_post(*SubS,*SubI,*this);
    #ifdef PRINT_ALLOC
    disp_alloc(to_string(parallelisation->rank)+" : Memoire apres allocations : ",1);
    #endif
    print("Repartition des solides");
    if(parallelisation->is_master_cpu()){
        for(int i = 0; i < parallelisation->repartition_sst.size(); i++){
            std::cout << "\t" << i << " :";
            unsigned nb_nodes = 0;
            for(int j = 0; j < parallelisation->repartition_sst[i].size(); j++){
                std::cout << "\t" << parallelisation->repartition_sst[i][j];
                nb_nodes += (*S)[parallelisation->repartition_sst[i][j]].mesh.elem_list_size;
            }
            std::cout << "\t" << nb_nodes << std::endl;
        }
    }
    
    /// Boucle sur les steps temporels
    print_title(1,"DEBUT DU CALCUL ITERATIF ");
    print_data("Nombre de pas de temps total : ",temps->nbpastemps);
    /// Calcul des operateurs  A DEPLACER VERS LE DEBUT DE LA BOUCLE ITERATIVE
    print_title(2,"Mise a jour des operateurs");
    for(int i_sst = 0; i_sst < S->size(); i_sst++){
        if((*S)[i_sst].update_operator){
            #ifdef PRINT_ALLOC
            disp_alloc((to_string(parallelisation->rank)+" : Verifie memoire avant construction : ").c_str(),1);
            #endif
            for(int i = 0; i < (*sst_materials).size(); i++){
                //(*sst_materials)[i].affiche();
            }
            multiscale_operateurs(*Stot,*SubS,*Inter,*SubI,*this,*Global, *data_user);
            Global->allocations(multiscale->sizeM);
            #ifdef PRINT_ALLOC
            disp_alloc((to_string(parallelisation->rank)+" : Verifie memoire apres construction : ").c_str(),1);
            #endif
            break;  /// multiscale_operateurs a remis a jour les operateurs de tout le monde
        }
    }
    #ifdef INFO_TIME
    print_duration(tic2);
    #endif
    
    
    
    for(temps->init();temps->has_next();temps->next()){
        if(temps->step_changed()){
            print_data("****************************** Step : ",temps->step_cur);
        }
        print_data("*************** Time : ",temps->t_cur);
        print_title(2,"Mise a jour des parametres");
        //updater.add_message( mp, updater.ET_Info, "pas de temps i" );
        
        temps->updateParameters();              /// Mise a jour des parametres temporels utilisateur
        Boundary::updateParameters();           /// Mise a jour des CL (PENSER A ENLEVER PLUS BAS LORSQUE PRET)
        InterCarac::updateParameters();         /// Mise a jour des parametres materiaux des interfaces
        
        
//         /// Calcul des operateurs  A DEPLACER VERS LE DEBUT DE LA BOUCLE ITERATIVE
//         print_title(2,"Mise a jour des operateurs");
//         for(int i_sst = 0; i_sst < S->size(); i_sst++){
//             if((*S)[i_sst].update_operator){
//                 #ifdef PRINT_ALLOC
//                 disp_alloc((to_string(parallelisation->rank)+" : Verifie memoire avant construction : ").c_str(),1);
//                 #endif
//                 for(int i = 0; i < (*sst_materials).size(); i++){
//                     //(*sst_materials)[i].affiche();
//                 }
//                 multiscale_operateurs(*Stot,*SubS,*Inter,*SubI,*this,*Global, *data_user);
//                 Global->allocations(multiscale->sizeM);
//                 #ifdef PRINT_ALLOC
//                 disp_alloc((to_string(parallelisation->rank)+" : Verifie memoire apres construction : ").c_str(),1);
//                 #endif
//                 break;  /// multiscale_operateurs a remis a jour les operateurs de tout le monde
//             }
//         }
//         #ifdef INFO_TIME
//         print_duration(tic2);
//         #endif
        print_title(2,"Mise a jour des Conditions d'interface");
        
        parallelisation->synchronisation();
        if(nom_calcul=="incr") {
            /// Presence d'interface Breakable ?
            int nb_breakable=0;
            if (parallelisation->is_master_cpu()){
                for(unsigned q=0; q <Inter->size();q++){
                    if ((*Inter)[q].comp =="Breakable"){
                        nb_breakable++;
                    }
                }
            }
            if (parallelisation->is_multi_cpu()){
                MPI_Bcast(&nb_breakable,1, MPI_INT, 0, MPI_COMM_WORLD);
            }
            nb_breakable = nb_breakable ;
            /// Mise a jour des conditions aux limites
            if(temps->pt_cur == 1 and parallelisation->is_local_cpu()){
                print_title(2,"    Initialisation des Conditions aux limites :");
                //for(int i = 0; i < SubI->size(); i++){
                //    (*SubI)[i].init();
                //}
                initialise_CL_values(*SubI, *CL);
            }
            parallelisation->synchronisation();
            if (parallelisation->is_local_cpu()){
                print_title(2,"    Mise a jour des Caracteristiques des interfaces :");
                for(int i = 0; i < SubI->size(); i++){
                    (*SubI)[i].init(temps->pt);
                }
                print_title(2,"    Mise a jour des Conditions aux limites :");
                update_CL_values(*SubI, *CL, *this, *data_user);
            }
            
            /// Calcul sur le pas de temps
            if (nb_breakable>0) {
                int nb_change = 0;
                int sous_iter = 1;
                while(nb_change != 0 or sous_iter == 1) {
                    if (parallelisation->is_local_cpu()){
                        for(unsigned q=0; q < SubI->size();q++){
                            if ((*SubI)[q].comp == "Breakable") {
                                (*SubI)[q].convergence = -1;
                            }
                        }
                    }
                    print_data("          Sous iteration interface cassable : ",sous_iter);
                    iterate_incr(*this,*SubS,*Inter,*SubI,*Global,*data_user);
                    if (parallelisation->is_local_cpu()) {
                        for(unsigned q=0; q < SubI->size();q++) {
                            if ((*SubI)[q].comp == "Breakable") {
                                nb_change += (*SubI)[q].convergence ;
                            }
                        }
                    }
                }
            }
            else {
                iterate_incr(*this,*SubS,*Inter,*SubI,*Global,*data_user);
            }
            
            ///assignation ptcur au ptold
            parallelisation->synchronisation();
            print_title(2,"Reactualisation des valeurs pour le pas de temps suivant");
            assign_quantities_current_to_old(*SubS,*SubI,*this);
            
            parallelisation->synchronisation();
	    //if (multiresolution->calcul_cur==0%5){
            affichage_resultats(*SubS,*this, *data_user);
            affichage_resultats_inter(*SubI, *S ,*this, *data_user);
	    //}
            /// Sauvegarde des resultats
            if(save_data){
                print_title(2,"Sauvegarde des resultats au format HDF"); 
                if (parallelisation->is_local_cpu()) {/* A REVOIR APRES MODIFICATION DATAUSER (pour field_structure_user) + MODIFIER format HDF5 pour la multi-resolution
                    write_hdf_fields_SST_INTER(*SubS, *Inter, *this , *data_user);//  BUG
                    convert_fields_to_field_structure_user(*SubS, *Inter, *this , *data_user, *field_structure_user, *geometry_user);
                    Sc2String rank; rank << parallelisation->rank;
                    Sc2String file_output_hdf5 = affichage->name_hdf + "_" + rank + ".h5";
                    field_structure_user->write_hdf5_in_parallel(file_output_hdf5, *geometry_user, affichage->name_fields, temps->pt_cur, temps->t_cur, parallelisation->rank);
                //*/
                }
            }
            
            /// Modification du comportement des entites
            //modification_sst_inter_behaviour(S,Inter,temps);  A TESTER
            
            print_data("*************** End time : ",temps->t_cur);
            
            
            ///Affichage des energies
            if (affichage->trac_ener_imp == 1) {
                affichage->param_ener[0]=1;
                affichage->param_ener[1]=0;
                affichage_energie(*SubS,*Inter,*this,*data_user);
                affichage->param_ener[0]=1;
                affichage->param_ener[1]=1;
                affichage_energie(*SubS,*Inter,*this,*data_user);
            }
            if (affichage->trac_ener_diss == 1) {
                affichage->param_ener[0]=0;
                affichage->param_ener[1]=0;
                affichage_energie(*SubS,*Inter,*this,*data_user);
                affichage->param_ener[0]=0;
                affichage->param_ener[1]=1;
                affichage_energie(*SubS,*Inter,*this,*data_user);
            }
        } else {
            std::cerr << "Nom de calcul non defini : incremental uniquement" << std::endl;
            assert(0);
        }
        /*
        std::cout << "******************************************************************************" << std::endl;
        const int side = affichage->side;
        const int pt = temps->pt_cur;
        for(int i = 0; i < (*Inter).size(); i++){
            std::cout << "Interface " << (*Inter)[i].id;
            std::cout << std::endl << "\tFchap : ";
            for(int j = 0; j < (*Inter)[i].side[side].t_post[pt].Fchap.size(); j++){
                std::cout << "\t" << (*Inter)[i].side[side].t_post[pt].Fchap[j] << std::endl;
            }
            std::cout << std::endl << "\tWpchap : ";
            for(int j = 0; j < (*Inter)[i].side[side].t_post[pt].Wpchap.size(); j++){
                std::cout << "\t" << (*Inter)[i].side[side].t_post[pt].Wpchap[j] << std::endl;
            }
        }
        std::cout << "******************************************************************************" << std::endl;*/
    }
    ///Affichage des résultantes sur les interfaces
    calcul_resultante(*SubS,*S,*Inter,*this);
    //if (multiresolution->calcul_cur==0%5){
    ///creation des fichiers pvd
    create_pvd_results(*SubS,*S,*Inter,*this);
    //}
    
    #ifdef INFO_TIME
    print_duration(tic2);
    #endif
    
    ///sortie xdmf à partir du fichier hdf5 cree au fur et à mesure du calcul
    if(parallelisation->is_master_cpu() and save_data==1){
        //write_xdmf_file_compute(*this, data_user);
    }
    
    //affichage_resultats(*SubS,*this, *data_user);            ///sortie paraview pour les sst (volumes et peaux)   /// TMP, test sauvegarde a la fin de chaque pas de temps
    //affichage_resultats_inter(*SubI, *S ,*this, *data_user); ///sortie paraview pour les interfaces               /// TMP, test sauvegarde a la fin de chaque pas de temps
}



void Process::finalisation_MPI(){
    print("Programme complet : ",true);
    #ifdef INFO_TIME
    print_duration(tic1);
    #endif
    if (parallelisation->is_multi_cpu()) MPI_Finalize();
}


void Process::read_data_user() {
    sousint = false;
    type_sousint = "h";
    rbm.bloq = false;
    rbm.mvts_bloques.resize(3);
    rbm.mvts_bloques[0]= "Ty";
    rbm.mvts_bloques[1]= "Tx";
    rbm.mvts_bloques[2]= "Rz";
    save_data = true;
    read_data = false;
    reprise_calcul = 0;
    //properties->deltaT = 0;
    
    if (parallelisation->is_master_cpu())
        if(data_user->time_steps.time_scheme == "statique"){
            std::cout << "Type de cacul : STATIQUE" << std::endl;
        }else if(data_user->time_steps.time_scheme == "quasistatique"){
            std::cout << "Type de cacul : QUASISTATIQUE" << std::endl;
        }
    nom_calcul = "incr";
    
    multiresolution->read_data_user(*data_user);
    structure->read_data_user(*data_user);
    multiscale->read_data_user(*data_user);
    latin->read_data_user(*data_user);
    affichage->read_data_user(*data_user);
    temps->read_data_user(*data_user);
    sst_materials->resize(data_user->materials_vec.size());
    for(unsigned i = 0; i < sst_materials->size(); i++){
        (*sst_materials)[i].read_data_user(i,*data_user);
    }
    S->resize(data_user->pieces_vec.size());
    for(unsigned i = 0; i < S->size(); i++){
        (*S)[i].read_data_user(i,*data_user,*geometry_user);
    }
    inter_materials->resize(data_user->links_vec.size());
    for(unsigned i = 0; i < inter_materials->size(); i++){
        (*inter_materials)[i].read_data_user(i,*data_user);
    }
    Inter->resize(data_user->interfaces_vec.size());
    for(unsigned i = 0; i < Inter->size(); i++){
        (*Inter)[i].read_data_user(i,*data_user,*geometry_user);
    }
    CL->resize(data_user->boundary_conditions_vec.size());
    for(unsigned i = 0; i < CL->size(); i++){
        (*CL)[i].read_data_user(i,*data_user);
    }
    Fvol->read_data_user(*data_user);
    Tload->read_data_user(*data_user);
    
    /// Creation des liens de parente entre groupes de parametres
    SstCarac::sst_materials_parameters.addParent(&(multiresolution->parameters));
    InterCarac::inter_materials_parameters.addParent(&(Boundary::CL_parameters));
    temps->parameters.addParent(&(multiresolution->parameters));
    Boundary::CL_parameters.addParent(&(temps->parameters));
    Fvol->parameters.addParent(&(temps->parameters));
    Tload->parameters.addParent(&(temps->parameters));
    
    /// Preparation des parametres
    multiresolution->prepareParameters();
    temps->prepareParameters();
    SstCarac::prepareParameters();
    InterCarac::prepareParameters();
    Boundary::prepareParameters();
    Fvol->prepareParameters();
    Tload->prepareParameters();
    
    /// Debuggage
    //multiresolution->affiche();
    //temps->affiche();
    //for(int i = 0; i < inter_materials->size(); i++){
    //    (*inter_materials)[i].affiche();
    //}
    for(int i = 0; i < CL->size(); i++){
        (*CL)[i].affiche();
    }
    //Fvol->affiche();
    //Tload->affiche();
};

void Process::print_title(int level,Sc2String title){
    if(not parallelisation->is_master_cpu())
        return;
    if(level == 0){
        /// Titre principal
        const int n = 78 - title.size(); /// Nombre d'espaces pour completer la ligne
        std::cout << std::endl << std::endl;
        std::cout << "********************************************************************************" << std::endl;
        std::cout << "*                                                                              *" << std::endl;
        std::cout << "*";
        const int n2 = n/2;
        for(int i = 0; i < n2; i++)
            std::cout << " ";
        std::cout << title;
        const int n3 = n - n2;
        for(int i = 0; i < n3; i++)
            std::cout << " ";
        std::cout << "*" << std::endl;
        std::cout << "*                                                                              *" << std::endl;
        std::cout << "********************************************************************************" << std::endl;
        std::cout << std::endl << std::endl;
    } else if(level == 1) {
        /// Grande partie
        const int n = title.size()+8;
        std::cout << std::endl << std::endl;
        for(int i = 0; i < n; i++)
            std::cout << "*";
        std::cout << std::endl << "*   " << title << "   *" << std::endl;
        for(int i = 0; i < n; i++)
            std::cout << "*";
        std::cout << std::endl;
    } else {
        /// Petite partie
        const int n = 79 - title.size();
        std::cout << std::endl << title << " ";
        for(int i = 0; i < n; i++)
            std::cout << "-";
        std::cout << std::endl;
    }
}

/// void print(Sc2String statement,bool no_endl = false) affiche 'statement' avec ('no_endl'=false) ou sans ('no_endl'=true) retour a la ligne
void Process::print(Sc2String statement,bool no_endl){
    if(parallelisation->is_master_cpu()){
        std::cout << statement;
        if(not no_endl) std::cout << std::endl;
    }
}

/// void print_duration(TicTac& tic) synchronise les CPU et affiche la duree du timer 'tic'
void Process::print_duration(TicTac& tic){
    #ifdef INFO_TIME
    parallelisation->synchronisation();
    if (parallelisation->is_master_cpu()){
        std::cout << "duree : ";
        tic.stop();
        std::cout << " secondes" << std::endl;
        tic.start();
    }
    #endif
}