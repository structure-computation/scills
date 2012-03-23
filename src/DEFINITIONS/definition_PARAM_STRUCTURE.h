#ifndef PARAM_STRUCTURE_H
#define PARAM_STRUCTURE_H


#include <string>
using namespace LMT;


#include "DataUser.h"


/** \ingroup Parametres
\brief Param�tres lies � la structure �tudi�e : nom des r�pertoires contenant les maillages et autres ...
*/
struct STRUCTURE{
   //constructeur
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
   
   double scale;  /// parametre d'homothetie des differentes quantites (pas utilis�) 
   std::string nom_fichier_qualification_materiaux;  /// nom du fichier contenant le numero du mat�riau associ� � chaque sst
   std::string repertoire_des_maillages; /// r�pertoire contenant les maillages
   unsigned nb_maillages; /// nombre de maillages lus (nbre de ssts)
   std::string nom_des_maillages; /// nom g�n�rique pour les diff�rents maillages (ex : S pour S1, S2 ...)
   std::string extension; /// extension des fichiers � lire (.avs par d�faut)
   //unsigned nb_maillages_CL; /// nombre de maillages de bord de forme compliquee (peu utilise)
   //string nom_des_maillages_CL; /// nom generique pour les maillages de bord(ex : L pour L1, L2 ...)
   
   unsigned jeu_physique; /// flag indiquant si l'on doit lire des maillages pour des jeux physiques
   std::string nom_maillages_jeu; /// nom generique pour les maillages d'interface avec jeu physique
   Vec< Vec<int> > inter_jeu;/// vecteur contenant une paire de num�ros de sous-structures voisines des interfaces avec jeu physique

   double volumetot; /// mesure totale de la structure etudiee
};

#endif //PARAM_STRUCTURE_H
