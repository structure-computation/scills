#ifndef PARAM_STRUCTURE_H
#define PARAM_STRUCTURE_H


#include "../UTILS/Sc2String.h"
using namespace LMT;


#include "../COMPUTE/DataUser.h"


/** \ingroup Parametres
\brief Paramètres lies à la structure étudiée : nom des répertoires contenant les maillages et autres ...
*/
struct STRUCTURE{
   double scale;  /// parametre d'homothetie des differentes quantites (pas utilisé) 
   Sc2String nom_fichier_qualification_materiaux;  /// nom du fichier contenant le numero du matériau associé à chaque sst
   Sc2String repertoire_des_maillages; /// répertoire contenant les maillages
   unsigned nb_maillages; /// nombre de maillages lus (nbre de ssts)
   Sc2String nom_des_maillages; /// nom générique pour les différents maillages (ex : S pour S1, S2 ...)
   Sc2String extension; /// extension des fichiers à lire (.avs par défaut)
   //unsigned nb_maillages_CL; /// nombre de maillages de bord de forme compliquee (peu utilise)
   //Sc2String nom_des_maillages_CL; /// nom generique pour les maillages de bord(ex : L pour L1, L2 ...)
   
   unsigned jeu_physique; /// flag indiquant si l'on doit lire des maillages pour des jeux physiques
   Sc2String nom_maillages_jeu; /// nom generique pour les maillages d'interface avec jeu physique
   Vec< Vec<int> > inter_jeu;/// vecteur contenant une paire de numéros de sous-structures voisines des interfaces avec jeu physique

   double volumetot; /// mesure totale de la structure etudiee
   
   
   ///constructeur
   STRUCTURE(){
       scale=1;
       extension=".avs";
   }
   
   void read_data_user(DataUser &data_user) {
       repertoire_des_maillages = data_user.model_path + "/MESH/";
       nom_fichier_qualification_materiaux = "";
       nom_des_maillages = "";
       nb_maillages = 0;
       extension = "";
       jeu_physique = 0;   
   }
};

#endif //PARAM_STRUCTURE_H
