#ifndef PARAM_STRUCTURE_H
#define PARAM_STRUCTURE_H
using namespace LMT;
using namespace std;

/** \ingroup Parametres
\brief Param�tres lies � la structure �tudi�e : nom des r�pertoires contenant les maillages et autres ...
*/
struct STRUCTURE{
   //constructeur
   STRUCTURE(){
      scale=1;
      extension=".avs";
   }
   
   double scale;  /// parametre d'homothetie des differentes quantites (pas utilis�) 
   string nom_fichier_qualification_materiaux;  /// nom du fichier contenant le numero du mat�riau associ� � chaque sst
   string repertoire_des_maillages; /// r�pertoire contenant les maillages
   unsigned nb_maillages; /// nombre de maillages lus (nbre de ssts)
   string nom_des_maillages; /// nom g�n�rique pour les diff�rents maillages (ex : S pour S1, S2 ...)
   string extension; /// extension des fichiers � lire (.avs par d�faut)
   //unsigned nb_maillages_CL; /// nombre de maillages de bord de forme compliquee (peu utilise)
   //string nom_des_maillages_CL; /// nom generique pour les maillages de bord(ex : L pour L1, L2 ...)
   
   unsigned jeu_physique; /// flag indiquant si l'on doit lire des maillages pour des jeux physiques
   string nom_maillages_jeu; /// nom generique pour les maillages d'interface avec jeu physique
   Vec< Vec<int> > inter_jeu;/// vecteur contenant une paire de num�ros de sous-structures voisines des interfaces avec jeu physique

   double volumetot; /// mesure totale de la structure etudiee
};

#endif //PARAM_STRUCTURE_H
