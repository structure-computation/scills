#ifndef PROCESS_H
#define PROCESS_H

// #if DIM == 2
//     #include "Scills2DUpdater.h"
// #else
//     #include "Scills3DUpdater.h"
// #endif

#include "../COMPUTE/DataUser.h"
#include "../COMPUTE/FieldStructureUser.h"
#include "../GEOMETRY/GeometryUser.h"

#include "structure_typedef.h"
#include "SavingData.h"
#include "GeneralData.h"
#include "MultiScaleData.h"
#include "LatinData.h"
#include "TimeData.h"
#include "PROPERTY.h"
#include "ParallelisationData.h"
#include "MultiResolutionData.h"
#include "MacroProblem.h"

#include "../UTILS/SCtime.h"


/** \defgroup Parametres
\brief Param?tres du programme multi?chelle 
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
    /// Structures de donnees associees a la methode LaTIn
    GeneralData         *structure;                     /// parametres donnant la geometrie et maillages de la structure
    LatinData           *latin;                         /// parametres lies a la strategie latin
    MultiScaleData      *multiscale;                    /// parametres multiechelle
    SavingData          *affichage;                     /// parametres d'affichage et de sauvegarde
    TimeData            *temps;                         /// parametres temporels
    //PROPERTY *properties;                   /// parametres pour les proprietes materielles
    ParallelisationData *parallelisation;               /// parametres MPI
    MultiResolutionData *multiresolution;               /// parametres de multi-resolution
    DataUser            *data_user;                     /// structure de stockage des informations du fichier JSON
    GeometryUser        *geometry_user;                 /// structure de stockage des informations du fichier HDF5
    FieldStructureUser  *field_structure_user;
    
    MacroProblem *Global;        /// propriete macro de la structure
    
    /// Donnees associees a la structure etudiee
    Substructures *S;            /// vecteur des sous structures (contient toutes les Sst du pb, mais seules celles dans SubS seront traitees par l'instance du programme)
    VecInterfaces *Inter;        /// vecteur des interfaces (idem que S mais pour les interfaces)
    PointedSubstructures *SubS;  /// MPI : vecteur de pointeur vers les Sst a traiter dans cette instance du programme
    PointedSubstructures *Stot;  /// MPI : vecteur de pointeur vers toutes les Sst (pour avoir le meme type que SubS)
    PointedInterfaces *SubI;     /// MPI : vecteur de pointeur vers les Interfaces a traiter dans cette instance du programme
    Boundaries *CL;              /// vecteur des conditions limites
    VolumicForces *Fvol;         /// efforts volumiques sur la structure
    ThermalLoad *Tload;          /// chargement thermique sur la structure
    Materials *sst_materials;    /// vecteur des propriete materiaux des Sst
    Links *inter_materials;      /// vecteur des propriete materiaux des Interfaces
    
    /// Operateurs de la structure
    
    /// donnees globales pour le calcul 
    bool sousint;           /// sousintegration oui = 1 non=0
    Sc2String type_sousint; /// type de sous-int?gration h ou p
    Sc2String nom_calcul;   /// nom du calcul effectu? : latin ou incremental (valable pour le quasistatique ou statique)
    bool recopie_t_post;    /// booleen pour recopier les quantites obtenues par un calcul incremental aux quantites pour un calcul latin
    bool save_data;         /// bool?en permettant de sauvegarder dans des fichiers textes les donn?es aux interfaces et aux sst
    bool read_data;         /// bool?en permettant de lire dans des fichiers textes les donn?es aux interfaces et aux sst
    int reprise_calcul;     /// booleen permettant de dire si on initialise le calcul avec les donnees sauvegarder dans save_sst et save_inter
    bool compute_operators; /// indique s'il faut (re)calculer les operateurs
    bool plasticite;        /// indique si au moins un solide a une loi de comportement plastique
    bool endommagement;     /// indique si au moins un solide a une loi de comportement endommageable
    int nb_breakable;       /// nombre d'interfaces cassables
    
    ///structure pour bloquer les modes de corps rigides manuellement
    struct RBM
    {
        bool bloq;                      /// blocage ou non des modes de corps rigide
        Vec<Sc2String> mvts_bloques;    /// liste des modes de corps rigides bloques (ex : "Tx Ty Rz")
    };
    RBM rbm;
    

    //methodes===============================================================================================
    /// Constructeur
    Process();
    
    /// Destructeur
    ~Process();
    /// Alloue la memoire des differentes structures de donnees (appele par le constructeur)
    void allocate();
    /// Libere la memoire des differentes structures de donnees (appele par le destructeur)
    void free();
    
    /// Initialisation de MPI
    void initialisation_MPI(int argc,char **argv);
    void initialisation_MPI_for_scwal(int argc,char **argv);
    /// test de synchronisation de MPI
    void test_MPI();
    /// Lecture des fichier JSON et HDF5
    void lecture_fichiers(Sc2String id_model, Sc2String id_calcul);
    /// test des structures de donnees
    void test_load_data();
    /// Preparation des structures de donnees
    void preparation_calcul();
    /// Gestion de la boucle de multi-resolution
    //template<class Updater_, class MP_> void boucle_multi_resolution(Updater_ updater, MP_ mp);
    void boucle_multi_resolution();
    template<class Updater_, class MP_> void boucle_multi_resolution(MP_ mp, Updater_ &updater);
    /// Gestion de la boucle temporelle
    void boucle_temporelle();
    template<class Updater_, class MP_> void boucle_temporelle(MP_ mp, Updater_ &updater);
    /// Terminaison de MPI
    void finalisation_MPI();
    
    /// Charge les donnees du DataUser dans les structures de donnees
    void read_data_user();
    
    /// Affiche le titre 'title', de niveau 'level'
    void print_title(int level,Sc2String title);
    /// Affiche 'statement' avec ('no_endl'=false) ou sans ('no_endl'=true) retour a la ligne
    void print(Sc2String statement,bool no_endl = false);
    /// Synchronise les CPU et affiche la duree du timer 'tic'
    void print_duration(TicTac& tic);
    /// Affiche la valeur 'value' precedee du message 'msg'
    template<typename Data> void print_data(const char* msg,const Data& value){
        if(parallelisation->is_master_cpu()) std::cout << msg << value << std::endl;
    }
    void system();
};


#endif // PROCESS_H
