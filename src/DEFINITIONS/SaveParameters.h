#ifndef PARAM_SC2AFFICHAGE_H
#define PARAM_SC2AFFICHAGE_H

#include <Metil/StructCompactor.h>
#include <Metil/CudaMetil.h>
#include <Metil/Hdf.h>

#include "../COMPUTE/DataUser.h"
#include "../UTILS/Sc2String.h"

#include "../../LMT/include/containers/vec.h"

using namespace Metil;
using namespace LMT;


/**\ingroup Parametres
\brief Parametres d'affichage accessibles dans la structure Param.
*/
struct SaveParameters{
    //attributs==============================================================================================
    bool affich_resultat;   ///< affichage du resultat de calcul
    bool affich_mesh;       ///< affichage du maillage apres construction
    Sc2String type_affichage;   ///< type d'affichage "Sinterieur", "Sbord", "Interface"
    bool display_error; ///< affichage graphique de l'erreur au cours des iterations
    Vec<Sc2String> display_fields;          ///< champs a afficher pour le resultat du calcul (ex : "dep sigma epsilon")
    Vec<Sc2String> display_fields_sst_bulk; ///< liste de tous les champs a afficher pour les sst volumiques
    Vec<Sc2String> display_fields_sst_skin; ///< liste de tous les champs a afficher pour les sst peau
    Vec<Sc2String> display_fields_inter;    ///< champs a afficher pour le resultat du calcul des interfaces, par defaut c est tout donc on extrait l ensemble des donnees au element et on remplit ce vecteur (ex : "dep qtrans")
    Sc2String save;             ///< sauvegarde ou affichage des resultats (save ou display)
    Sc2String repertoire_save;  ///< repertoire de sauvegarde (ex : ./tmp/data/)
    Sc2String name_data;        ///< nom generique des fichiers de sauvegarde
    Sc2String command_file;     ///< nom du fichier de commande pour l interactivite
    bool interactivite;         ///< Booléen permettant de post traiter les resultats de maniere interactive sans relancer un calcul. Taper help pour connaitre les mots clés
    bool affich_inter_data;     ///< affichage des efforts et deplacements d'interface
    Vec<int> num_inter_select;  ///< numero des interfaces selectionnees pour l'affichage des efforts et déplacement (-1 pour toutes les interfaces)
    unsigned side;          ///< cote selectionne pour affichage des efforts et deplacement
    unsigned pt;            ///< pas de temps a afficher
    bool affich_depl_pt;    ///< affichage du deplacement d'un point sous gnuplot
    Vec<TYPEREEL> coor_point; ///< coordonnees du point a tracer
    bool fichiers_paraview_sst_crees;   ///< booleen indiquant si les fichiers paraview des sst ont ete crees
    bool fichiers_paraview_inter_crees; ///< booleen indiquant si les fichiers paraview des interfaces ont ete crees
    Vec<int,3> param_ener; ///< parametres pour l'affichage des energies dissipees ou imposees
    bool trac_ener_imp;
    bool trac_ener_diss;

    Sc2String name_hdf;             ///< nom du fichier hdf5 permettant la sauvegarde des donnees
    Sc2String name_geometry;        ///< nom du dataset contenant les elements de geometry dans le fichier de sauvegarde
    Sc2String name_fields;          ///< nom du dataset contenant les elements de geometry dans le fichier de sauvegarde
    Sc2String name_xdmf_geometry;   ///< nom du fichier xdmf de sortie pour visualisation des resultats sous paraview, geometrie uniquement
    Sc2String name_xdmf_fields;     ///< nom du fichier xdmf de sortie pour visualisation des resultats sous paraview, geometrie + champs solutions
    
    
    //methodes===============================================================================================
    SaveParameters();     ///< Constructeur
    void read_data_user(DataUser &data_user);   ///< Chargement des donnees depuis le DataUser
    void display_all_data();   ///< Affiche toutes les informations stockees dans cet objet
};

#endif // PARAM_SC2AFFICHAGE_H

