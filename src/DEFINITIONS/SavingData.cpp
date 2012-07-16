#include "SavingData.h"

// Pour les fonctions de debuggage
#include "../UTILS/utils_2.h"


SavingData::SavingData(){    
    affich_resultat=0;
    affich_mesh=0;
    display_error=0;
    type_affichage="Sinterieur";
    save="save";
    repertoire_save="./tmp";
    interactivite=0;
    pt=1;
    side=0;
    num_inter_select=Vec<int>(-1);
    fichiers_paraview_sst_crees=0;
    fichiers_paraview_inter_crees=0;
    param_ener.set(0);   
    trac_ener_imp = 0;
    trac_ener_diss = 0;
}

void SavingData::read_data_user(DataUser &data_user){
    interactivite= 0;
    affich_resultat= 1;
    if(data_user.options.mode == "test"){
        type_affichage= "Inter";
        affich_mesh= 1;
    }else{
        type_affichage= "Sbord";
        affich_mesh= 0;
    }
    display_error= 0; 
    save= "save";

    display_fields_sst_bulk.resize(15);
    display_fields_sst_bulk[0]= "dep";
    display_fields_sst_bulk[1]= "qtrans";
    display_fields_sst_bulk[2]= "sigma";
    display_fields_sst_bulk[3]= "epsilon";
    display_fields_sst_bulk[4]= "ener";
    display_fields_sst_bulk[5]= "sigma_von_mises";
    display_fields_sst_bulk[6]= "numsst";
    display_fields_sst_bulk[7]= "f_vol_e";
    display_fields_sst_bulk[8]= "num_proc";
    display_fields_sst_bulk[9]= "plast_cumulee";
    display_fields_sst_bulk[10]= "plast_ecrouissage";
    display_fields_sst_bulk[11]= "epsilon_p";
    display_fields_sst_bulk[12]= "d1";
    display_fields_sst_bulk[13]= "d2";
    display_fields_sst_bulk[14]= "df";

    display_fields_sst_skin.resize(8);
    display_fields_sst_skin[0]= "dep";
    display_fields_sst_skin[1]= "qtrans";
    display_fields_sst_skin[2]= "sigma_skin";
    display_fields_sst_skin[3]= "epsilon_skin";
    display_fields_sst_skin[5]= "sigma_mises_skin";
    display_fields_sst_skin[6]= "numsst_skin";
    display_fields_sst_skin[7]= "num_proc_skin";

    if(type_affichage== "Sinterieur"){
        display_fields=display_fields_sst_bulk;
    } else if(type_affichage== "Sbord"){
        display_fields=display_fields_sst_skin;
    }
            
    repertoire_save= data_user.calcul_path + "/";
    name_data= "result";
    command_file= "No";
    //trac_ener_imp  = data_user.options.trac_ener_imp;     A REVOIR : N'EXISTE PLUS DANS LE DATAUSER
    //trac_ener_diss = data_user.options.trac_ener_diss;    A REVOIR : IDEM
}


void SavingData::display_all_data(){
    std::cout << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "******************* Debug SavingData : ******************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    debug("bool affich_resultat    ",affich_resultat);
    debug("bool affich_mesh        ",affich_mesh);
    debug("Sc2String type_affichage",type_affichage);
    debug("bool display_error      ",display_error);
    debug("Vec<Sc2String> display_fields         ",display_fields);
    debug("Vec<Sc2String> display_fields_sst_bulk",display_fields_sst_bulk);
    debug("Vec<Sc2String> display_fields_sst_skin",display_fields_sst_skin);
    debug("Vec<Sc2String> display_fields_inter   ",display_fields_inter);
    debug("Sc2String save           ",save);
    debug("Sc2String repertoire_save",repertoire_save);
    debug("Sc2String name_data      ",name_data);
    debug("Sc2String command_file   ",command_file);
    debug("bool interactivite       ",interactivite);
    debug("bool affich_inter_data   ",affich_inter_data);
    debug("Vec<int> num_inter_select",num_inter_select);
    debug("unsigned side            ",side);
    debug("unsigned pt              ",pt);
    debug("bool affich_depl_pt      ",affich_depl_pt);
    debug("Vec<TYPEREEL> coor_point   ",coor_point);
    debug("bool fichiers_paraview_sst_crees  ",fichiers_paraview_sst_crees);
    debug("bool fichiers_paraview_inter_crees",fichiers_paraview_inter_crees);
    debug("Vec<int,3> param_ener",param_ener);
    debug("bool trac_ener_imp   ",trac_ener_imp);
    debug("bool trac_ener_diss  ",trac_ener_diss);
    debug("Sc2String name_hdf          ",name_hdf);
    debug("Sc2String name_geometry     ",name_geometry);
    debug("Sc2String name_fields       ",name_fields);
    debug("Sc2String name_xdmf_geometry",name_xdmf_geometry);
    debug("Sc2String name_xdmf_fields  ",name_xdmf_fields);
    std::cout << "*************************************************************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << std::endl;
}
