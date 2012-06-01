#ifndef PROCESS_H
#define PROCESS_H


#include "GeneralParameters.h"
#include "LatinParameters.h"
#include "MultiScaleParameters.h"
#include "SaveParameters.h"
#include "TimeParameters.h"
#include "PROPERTY.h"
#include "ParallelisationParameters.h"

#include "../COMPUTE/DataUser.h"
#include "../UTILS/Sc2String.h"
#include "../UTILS/SCtime.h"

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
    GeneralParameters *structure;                   /// parametres donnant la geometrie et maillages de la structure
    LatinParameters *latin;                         /// parametres lies a la strategie latin
    MultiScaleParameters *multiscale;               /// parametres multiechelle
    SaveParameters *affichage;                      /// parametres d'affichage et de sauvegarde
    TimeParameters *temps;                          /// parametres temporels
    PROPERTY *properties;                           /// parametres pour les proprietes materielles
    ParallelisationParameters *parallelisation;     /// parametres MPI
    
    /// donnees globales pour le calcul 
    bool sousint;           /// sousintegration oui = 1 non=0
    Sc2String type_sousint; /// type de sous-intégration h ou p
    Sc2String nom_calcul;   /// nom du calcul effectué : latin ou incremental (valable pour le quasistatique ou statique)
    bool recopie_t_post;    /// booleen pour recopier les quantites obtenues par un calcul incremental aux quantites pour un calcul latin
    bool save_data;         /// booléen permettant de sauvegarder dans des fichiers textes les données aux interfaces et aux sst
    bool read_data;         /// booléen permettant de lire dans des fichiers textes les données aux interfaces et aux sst
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
    Process();  /// Constructeur
    ~Process(); /// Destructeur
    void read_data_user(DataUser &data_user);   /// Charge les donnees du DataUser dans les structures de donnees ci-dessous
    void allocation_memoire();      /// Alloue la memoire des differentes structures de donnees (appele par le constructeur)
    void desallocation_memoire();   /// Libere la memoire des differentes structures de donnees (appele par le destructeur)
    
    /// void print_title(int level,Sc2String title) affiche le titre 'title', de niveau 'level'
    void print_title(int level,Sc2String title){
        if(not parallelisation->is_master_cpu())
            return;
        if(level == 0){
            /// Titre principal
            const int n = 78 - title.size(); /// Nombre d'espaces pour completer la ligne
            std::cout << std::endl << std::endl;
            std::cout << "********************************************************************************" << std::endl;
            std::cout << "*                                                                              *" << std::endl;
            std::cout << "*";
            const int n2 = n/2;
            for(int i = 0; i < n2; i++)
                std::cout << " ";
            std::cout << title;
            const int n3 = n - n2;
            for(int i = 0; i < n3; i++)
                std::cout << " ";
            std::cout << "*" << std::endl;
            std::cout << "*                                                                              *" << std::endl;
            std::cout << "********************************************************************************" << std::endl;
            std::cout << std::endl << std::endl;
        } else if(level == 1) {
            /// Grande partie
            const int n = title.size()+8;
            std::cout << std::endl << std::endl;
            for(int i = 0; i < n; i++)
                std::cout << "*";
            std::cout << std::endl << "*   " << title << "   *" << std::endl;
            for(int i = 0; i < n; i++)
                std::cout << "*";
            std::cout << std::endl;
        } else {
            /// Petite partie
            const int n = 79 - title.size();
            std::cout << std::endl << title << " ";
            for(int i = 0; i < n; i++)
                std::cout << "-";
            std::cout << std::endl;
        }
    }
    
    /// void print(Sc2String statement,bool no_endl = false) affiche 'statement' avec ('no_endl'=false) ou sans ('no_endl'=true) retour a la ligne
    void print(Sc2String statement,bool no_endl = false){
        if(parallelisation->is_master_cpu()){
            std::cout << statement;
            if(not no_endl) std::cout << std::endl;
        }
    }
    
    /// void print_duration(TicTac& tic) synchronise les CPU et affiche la duree du timer 'tic'
    void print_duration(TicTac& tic){
        parallelisation->synchronisation();
        if (parallelisation->is_master_cpu()){
            std::cout << "duree : ";
            tic.stop();
            std::cout << " secondes" << std::endl;
            tic.start();
        }
    }
    
    /// template< typename Data > void print_data(const char* msg,const Data& value) affiche la valeur 'value' precedee du message 'msg'
    template<typename Data> void print_data(const char* msg,const Data& value){
        if(parallelisation->is_master_cpu()) std::cout << msg << value << std::endl;
    }
};


#endif // PROCESS_H
