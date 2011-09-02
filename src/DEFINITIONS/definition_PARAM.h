#ifndef PARAM_H
#define PARAM_H
using namespace LMT;
using namespace std;

/** \defgroup Parametres
\brief Paramètres du programme multiéchelle 
*/

struct AFFICHAGE; // parametres d'affichage
struct STRUCTURE; // parametres de maillages et geometries
struct LATIN; // parametres associes a la resolution LATIN
struct MULTI; // parametres multiechelles
struct TEMPS; // parametres temporels
struct PROPERTY; // parametres de proprietes materielles globales
struct MULTI_MPI; // parametres de MPI

//****************************
// parametres de calcul
//****************************
/** \ingroup Parametres
\brief classe Parametres de calcul contient toutes les informations sur l'affichage, les parametres de calcul, ...
                   
Pour les classes telles que structure, latin, multiscale, affichage, temps et properties qui sont definies en tant que pointeurs, les donnees sont accessibles par la commande :
Param process;
process.structure->nb_maillages par exemple
*/
struct Param
{
  Param()
  {
    // initialisation des valeurs
    nb_threads=1;
    sousint=1;
    type_sousint="p";
    dim=2;
    nom_calcul="latin";
    recopie_t_post=0;
    save_data=0;
  }

  /// parametres donnant la geometrie et maillages de la structure
  STRUCTURE *structure;
  /// parametres lies a la strategie latin
  LATIN *latin;
  /// parametres multiechelle
  MULTI *multiscale;
  /// parametres affichage et sauvegarde
  AFFICHAGE *affichage;

  //donnees globales pour le calcul 
  unsigned dim; ///< Dimension pour le calcul multiéchelle permettant de sélectionner sans recompiler la dimension du problème
  bool sousint;///< sousintegration oui = 1 non=0
  string type_sousint; ///< type de sous-intégration h ou p
  ///structure pour bloquer les modes de corps rigides manuellement
  struct RBM
  {
    bool bloq;      ///< blocage ou non des modes de corps rigide
    Vec<string> mvts_bloques; ///< liste des modes de corps rigides bloques (ex : "Tx Ty Rz")
  };
  RBM rbm;
  int nb_threads; ///< nombre de threads par machine pour parallelisme
   
  ///parametres temporels
  TEMPS *temps;
    
  /// parametres pour les proprietes materielles
  PROPERTY *properties;

  ///parametres MPI
  MULTI_MPI *multi_mpi;
  int rank,size;
  
  string nom_calcul; ///< nom du calcul effectué : latin ou incremental (valable pour le quasistatique ou statique)
  bool recopie_t_post; ///< booleen pour recopier les quantites obtenues par un calcul incremental aux quantites pour un calcul latin
  bool save_data; ///< booléen permettant de sauvegarder dans des fichiers textes les données aux interfaces et aux sst
  bool read_data; ///< booléen permettant de lire dans des fichiers textes les données aux interfaces et aux sst
  int reprise_calcul; ///< booleen permettant de dire si on initialise le calcul avec les donnees sauvegarder dans save_sst et save_inter
  
  int nb_breakable; ///nombre d'interfaces cassables
  
};


#endif // PARAM_H
