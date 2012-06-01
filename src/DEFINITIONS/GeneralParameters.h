#ifndef GENERALPARAMETERS_H
#define GENERALPARAMETERS_H


#include "../UTILS/Sc2String.h"
#include "../COMPUTE/DataUser.h"
#include "../../LMT/include/containers/vec.h"
using namespace LMT;



/** \ingroup Parametres
\brief Param�tres lies � la structure etudiee : nom des repertoires contenant les maillages et autres ...
*/
struct GeneralParameters{
    //attributs==============================================================================================
    //TYPEREEL scale;       /// parametre d'homothetie des differentes quantites (pas utilise) 
    Sc2String nom_fichier_qualification_materiaux;  /// nom du fichier contenant le numero du materiau associe � chaque sst
    Sc2String repertoire_des_maillages;             /// repertoire contenant les maillages
    unsigned nb_maillages;          /// nombre de maillages lus (nbre de ssts)
    Sc2String nom_des_maillages;    /// nom generique pour les differents maillages (ex : S pour S1, S2 ...)
    Sc2String extension;            /// extension des fichiers � lire (.avs par defaut)
    
    unsigned jeu_physique; /// flag indiquant si l'on doit lire des maillages pour des jeux physiques
    Sc2String nom_maillages_jeu; /// nom generique pour les maillages d'interface avec jeu physique
    Vec< Vec<int> > inter_jeu;/// vecteur contenant une paire de numeros de sous-structures voisines des interfaces avec jeu physique

    TYPEREEL volumetot; /// mesure totale de la structure etudiee
   
   
    //methodes===============================================================================================
    GeneralParameters();
    void read_data_user(DataUser &data_user);
    void display_all_data();
};

#endif //GENERALPARAMETERS_H
