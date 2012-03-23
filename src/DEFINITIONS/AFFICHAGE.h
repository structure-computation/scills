#ifndef PARAM_AFFICHAGE_H
#define PARAM_AFFICHAGE_H

#include <Metil/StructCompactor.h>
#include <Metil/CudaMetil.h>
#include <Metil/Hdf.h>

using namespace Metil;
using namespace LMT;


#include "../COMPUTE/DataUser.h"
#include "../UTILS/Sc2String.h"

// Pour les fonctions de debuggage
#include "../../LMT/include/io/ioexception.h"
#include "../UTILS/utils_2.h"


/**\ingroup Parametres
\brief Parametres d'affichage accessibles dans la structure Param.
*/
struct AFFICHAGE{
    // parametres pour l'affichage
    bool affich_resultat; ///< affichage du resultat de calcul
    bool affich_mesh;   ///< affichage du maillage apres construction
    Sc2String type_affichage;   ///< type d'affichage "Sinterieur", "Sbord", "Interface"
    bool display_error; ///< affichage graphique de l'erreur au cours des iterations
    Vec<Sc2String> display_fields; ///< champs a afficher pour le resultat du calcul (ex : "dep sigma epsilon")
    Vec<Sc2String> display_fields_sst_bulk; ///< liste de tous les champs a afficher pour les sst volumiques
    Vec<Sc2String> display_fields_sst_skin; ///< liste de tous les champs a afficher pour les sst peau
    Vec<Sc2String> display_fields_inter; ///< champs a afficher pour le resultat du calcul des interfaces, par defaut c est tout donc on extrait l ensemble des donnees au element et on remplit ce vecteur (ex : "dep qtrans")
    Sc2String save;         ///< sauvegarde ou affichage des resultats (save ou display)
    Sc2String repertoire_save;   ///< repertoire de sauvegarde (ex : ./tmp/data/)
    Sc2String name_data;   ///< nom generique des fichiers de sauvegarde
    Sc2String command_file;///< nom du fichier de commande pour l interactivite
    bool interactivite; ///< Booléen permettant de post traiter les resultats de maniere interactive sans relancer un calcul. Taper help pour connaitre les mots clés
    bool affich_inter_data; ///< affichage des efforts et deplacements d'interface
    Vec<int> num_inter_select; ///< numero des interfaces sélectionnées pour l'affichage des efforts et déplacement (-1 pour toutes les interfaces)
    unsigned side; ///< cote selectionne pour affichage des efforts et deplacement
    unsigned pt; ///< pas de temps a afficher
    bool affich_depl_pt; ///< affichage du deplacement d'un point sous gnuplot
    Vec<double> coor_point; ///< coordonnees du point a tracer
    bool fichiers_paraview_sst_crees;///< booleen indiquant si les fichiers paraview des sst ont ete crees
    bool fichiers_paraview_inter_crees;///< booleen indiquant si les fichiers paraview des interfaces ont ete crees
    Vec<int,3> param_ener; ///< parametres pour l'affichage des energies dissipees ou imposees
    bool trac_ener_imp;
    bool trac_ener_diss;

    Sc2String name_hdf; ///< nom du fichier hdf5 permettant la sauvegarde des donnees
    Sc2String name_geometry; ///< nom du dataset contenant les elements de geometry dans le fichier de sauvegarde
    Sc2String name_fields; ///< nom du dataset contenant les elements de geometry dans le fichier de sauvegarde
    Sc2String name_xdmf_geometry;   ///< nom du fichier xdmf de sortie pour visualisation des resultats sous paraview, geometrie uniquement
    Sc2String name_xdmf_fields;   ///< nom du fichier xdmf de sortie pour visualisation des resultats sous paraview, geometrie + champs solutions
    
    
    //constructeur : initialisation
    AFFICHAGE(){    
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
    
    void read_data_user(DataUser &data_user){
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
        
        display_fields_sst_bulk.resize(10);
        display_fields_sst_bulk[0]= "dep";
        display_fields_sst_bulk[1]= "qtrans";
        display_fields_sst_bulk[2]= "sigma";
        display_fields_sst_bulk[3]= "epsilon";
        display_fields_sst_bulk[4]= "ener";
        display_fields_sst_bulk[5]= "sigma_mises";
        display_fields_sst_bulk[6]= "numsst";
        display_fields_sst_bulk[7]= "f_vol_e";
        display_fields_sst_bulk[8]= "num_proc";
        
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
        }
        else if(type_affichage== "Sbord"){
            display_fields=display_fields_sst_skin;
        }
        
        repertoire_save= data_user.calcul_path + "/";
        name_data= "result";
        command_file= "No";
        trac_ener_imp  = data_user.options.trac_ener_imp;
        trac_ener_diss = data_user.options.trac_ener_diss;
    }
    
    void debug_all(){
        std::cout << std::endl << std::endl;
        std::cout << "***********************************************************" << std::endl;
        std::cout << "******************** Debug AFFICHAGE : ********************" << std::endl;
        std::cout << "***********************************************************" << std::endl;
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
        debug("Vec<double> coor_point   ",coor_point);
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
        std::cout << "***********************************************************" << std::endl;
        std::cout << "***********************************************************" << std::endl;
        std::cout << std::endl << std::endl;
    }
};

#endif // PARAM_AFFICHAGE_H

