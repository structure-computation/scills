#ifndef PROCESS_H
#define PROCESS_H


#include "GeneralParameters.h"
#include "LatinParameters.h"
#include "MultiScaleParameters.h"
#include "SaveParameters.h"
#include "TimeParameters.h"
#include "PROPERTY.h"
#include "MPIParameters.h"

#include "../COMPUTE/DataUser.h"
#include "../UTILS/Sc2String.h"

#include "../../LMT/include/containers/vec.h"
using namespace LMT;

/** \defgroup Parametres
\brief Paramètres du programme multiéchelle 
*/

//****************************
// parametres de calcul
//****************************
/** \ingroup Parametres
\brief classe Parametres de calcul contient toutes les informations sur l'affichage, les parametres de calcul, ...
                   
Cette classe regroupe toutes les donnees d'un calcul multi-echelle de la methode LaTIn.
*/
struct Process
{
    //attributs==============================================================================================
    /// Principales structures de donnees
    GeneralParameters *structure;       ///< parametres donnant la geometrie et maillages de la structure
    LatinParameters *latin;             ///< parametres lies a la strategie latin
    MultiScaleParameters *multiscale;   ///< parametres multiechelle
    SaveParameters *affichage;          ///< parametres d'affichage et de sauvegarde
    TimeParameters *temps;              ///< parametres temporels
    PROPERTY *properties;               ///< parametres pour les proprietes materielles
    MPIParameters *multi_mpi;           ///< parametres MPI
    
    /// Donnees pour le parallelisme
    int size;                           ///< nombre de processeurs alloues au calcul
    int rank;                           ///< numero du processeur gerant ce programme
    int nb_threads;                     ///< nombre de threads par machine
    
    /// donnees globales pour le calcul 
    bool sousint;           ///< sousintegration oui = 1 non=0
    Sc2String type_sousint; ///< type de sous-intégration h ou p
    Sc2String nom_calcul;   ///< nom du calcul effectué : latin ou incremental (valable pour le quasistatique ou statique)
    bool recopie_t_post;    ///< booleen pour recopier les quantites obtenues par un calcul incremental aux quantites pour un calcul latin
    bool save_data;         ///< booléen permettant de sauvegarder dans des fichiers textes les données aux interfaces et aux sst
    bool read_data;         ///< booléen permettant de lire dans des fichiers textes les données aux interfaces et aux sst
    bool plasticite;        ///< indique si au moins un solide a une loi de comportement plastique
    int reprise_calcul;     ///< booleen permettant de dire si on initialise le calcul avec les donnees sauvegarder dans save_sst et save_inter
    int nb_breakable;       ///< nombre d'interfaces cassables
    
    ///structure pour bloquer les modes de corps rigides manuellement
    struct RBM
    {
        bool bloq;                      ///< blocage ou non des modes de corps rigide
        Vec<Sc2String> mvts_bloques;    ///< liste des modes de corps rigides bloques (ex : "Tx Ty Rz")
    };
    RBM rbm;
    

    //methodes===============================================================================================
    Process();  ///< Constructeur
    ~Process(); ///< Destructeur
    void read_data_user(DataUser &data_user);   ///< Charge les donnees du DataUser dans les structures de donnees ci-dessous
    void allocation_memoire();      ///< Alloue la memoire des differentes structures de donnees (appele par le constructeur)
    void desallocation_memoire();   ///< Libere la memoire des differentes structures de donnees (appele par le destructeur)
};


#endif // PROCESS_H
