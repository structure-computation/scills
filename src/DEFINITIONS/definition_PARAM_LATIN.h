#ifndef PARAM_LATIN_H
#define PARAM_LATIN_H
using namespace LMT;
using namespace std;

/** \ingroup Parametres
\brief Paramètres associés à la stratégie LATIN
*/
struct LATIN{
  ///Constructeur : valeurs par défaut pour les paramètres
  LATIN(){
    iter=0;
    nbitermax=100;
    mu=0.8;
    ktype="scalaire_auto_CL";
    kfact=200e3;
    copydirection=1;
    critere_erreur=1e-4;
    facteur_relaxation=mu;
    type_error="ddr";
    list_error=0;
    save_depl_SST=1;
    alloc_quantites=1;
  }
  int iter;         ///< numero de l'itération courante
  int nbitermax;      ///< nombre d'itérations max
  double mu;         ///< facteur de relaxation courant  (1 à la première itération puis égal au facteur de relaxation)
  double critere_erreur; ///< valeur du critère d'erreur à atteindre
  double critere_erreur_diss; ///< valeur du critère d'erreur en dissipation à atteindre (si = 0 on ne la calcul pas)
  Vec<double> error; ///< erreurs au cours des iterations
  string type_error; ///< type d'erreur : choix possible entre "ddr" : erreur basée sur les directions de recherche, "energie" : erreur basée sur un calcul d'énergie par interface 
  bool list_error; ///< listage de l'erreur a l'écran au cours des itérations
  double facteur_relaxation; ///< facteur de relaxation entré par l'utilisateur
  bool save_depl_SST; ///< booléen pour sauvegarder le vecteur de déplacement à chaque pas de temps
  bool alloc_quantites; ///< booléen indiquant s'il est necessaire d'allouer les quantites et de realiser une etape prelocale (defaut 1)
  
  /// parametres pour les directions de recherche
  double kfact;    ///< facteur multiplicatif de la direction de recherche 
  unsigned copydirection; ///< copie de la direction de recherche de part et d'autre d'une interface
  /** 
  \brief Type de direction de recherche (utilise kfact)
  
   Plusieurs choix sont possibles selon le nom donné à ce paramètre : 
   - "scalaire_donne" : un paramètre scalaire identique pour toutes les interfaces donné par kfact
   - "scalaire_donne_CL" : un paramètre scalaire identique pour toutes les interfaces avec modification pour les interfaces de CL ou de type contact
   - "scalaire_auto" : un paramètre scalaire identique pour toutes les interfaces obtenu en prenant le module d'young de la sous-structure voisine divisé par une longueur caractéristique de l'interface considérée. Ce facteur est multiplié par kfact
   - "scalaire_auto_CL" : idem en modificant les valeurs pour les interfaces à CL ou de type contact
   
   Dans chaque cas, on calcule une direction normale et tangentielle (identiques pour la plupart des cas) : \f$ k_n \f$ et \f$ k_t \f$
   Les modifications pour les CL ou les interfaces de type contact sont les suivantes :
    - pour les interfaces de type effort, on multiplie les scalaires par 1/1000
    - pour les interfaces de type déplacement  on multiplie par 1000 les scalaires \f$ k_t \f$ et \f$ k_n \f$
    - pour les interfaces de type symétrie ou déplacement normal donné, on multiplie par 1000 \f$ k_n \f$ et par 1/1000 \f$ k_t \f$
    - pour les interfaces de type Contact, on multiplie par la coefficient de frottement \f$ k_t = f * k_n \f$   
  */
  string ktype;
};

#endif //PARAM_LATIN_H

