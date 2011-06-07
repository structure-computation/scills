#ifndef PARAM_AFFICHAGE_H
#define PARAM_AFFICHAGE_H
using namespace LMT;
using namespace std;

/**\ingroup Parametres
\brief Parametres d'affichage accessibles dans la structure Param.
*/
struct AFFICHAGE{
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
    }
  // parametres pour l'affichage
  bool affich_resultat; ///< affichage du resultat de calcul
  bool affich_mesh;   ///< affichage du maillage apres construction
  string type_affichage;   ///< type d'affichage "Sinterieur", "Sbord", "Interface"
  bool display_error; ///< affichage graphique de l'erreur au cours des iterations
  Vec<string> display_fields; ///< champs a afficher pour le resultat du calcul (ex : "dep sigma epsilon")
  Vec<string> display_fields_inter; ///< champs a afficher pour le resultat du calcul des interfaces, par defaut c est tout donc on extrait l ensemble des donnees au element et on remplit ce vecteur (ex : "dep qtrans")
  string save;         ///< sauvegarde ou affichage des resultats (save ou display)
  string repertoire_save;   ///< repertoire de sauvegarde (ex : ./tmp/data/)
  string name_data;   ///< nom generique des fichiers de sauvegarde
  string command_file;///< nom du fichier de commande pour l interactivite
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
};

#endif // PARAM_AFFICHAGE_H

