#ifndef GENERALPARAMETERS_H
#define GENERALPARAMETERS_H


#include "../UTILS/Sc2String.h"
#include "../COMPUTE/DataUser.h"
#include "../../LMT/include/containers/vec.h"
using namespace LMT;



/** \ingroup Parametres
\brief Paramètres lies à la structure étudiée : nom des répertoires contenant les maillages et autres ...
*/
struct GeneralParameters{
    //attributs==============================================================================================
    double scale;       /// parametre d'homothetie des differentes quantites (pas utilisé) 
    Sc2String nom_fichier_qualification_materiaux;  /// nom du fichier contenant le numero du matériau associé à chaque sst
    Sc2String repertoire_des_maillages;             /// répertoire contenant les maillages
    unsigned nb_maillages;          /// nombre de maillages lus (nbre de ssts)
    Sc2String nom_des_maillages;    /// nom générique pour les différents maillages (ex : S pour S1, S2 ...)
    Sc2String extension;            /// extension des fichiers à lire (.avs par défaut)
    
    unsigned jeu_physique; /// flag indiquant si l'on doit lire des maillages pour des jeux physiques
    Sc2String nom_maillages_jeu; /// nom generique pour les maillages d'interface avec jeu physique
    Vec< Vec<int> > inter_jeu;/// vecteur contenant une paire de numéros de sous-structures voisines des interfaces avec jeu physique

    double volumetot; /// mesure totale de la structure etudiee
   
   
    //methodes===============================================================================================
    GeneralParameters();
    void read_data_user(DataUser &data_user);
    void display_all_data();
};

#endif //GENERALPARAMETERS_H
