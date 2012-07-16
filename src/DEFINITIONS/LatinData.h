#ifndef LATINDATA_H
#define LATINDATA_H

#include "main_typedef.h"
#include "../COMPUTE/DataUser.h"

/** \ingroup Parametres
\brief Param�tres associes � la strategie LATIN
*/
struct LatinData{
    //attributs==============================================================================================
    int iter;               ///< numero de l'iteration courante
    int nbitermax;          ///< nombre d'iterations max
    
    Scalar facteur_relaxation;    ///< facteur de relaxation entre par l'utilisateur
    Scalar mu;                    ///< facteur de relaxation courant  (1 � la premi�re iteration puis facteur_relaxation)
    
    Scalar critere_erreur;              ///< valeur du crit�re d'erreur � atteindre
    Scalar critere_erreur_diss;         ///< valeur du crit�re d'erreur en dissipation � atteindre (si = 0 on ne la calcule pas)
    Scalar critere_erreur_auto_stop;    ///< valeur pour laquelle l'algorithme s'arrete si la variation d'erreur est inferieure 5 fois de suite � celle ci (si = 0 on ne fait pas d'arret automatique))
    
    bool list_error;        ///< listage de l'erreur a l'ecran au cours des iterations
    Sc2String type_error;   ///< type d'erreur : choix possible entre "ddr" : erreur basee sur les directions de recherche, "energie" : erreur basee sur un calcul d'energie par interface 
    Vector error;           ///< erreurs au cours des iterations
    
    bool save_depl_SST;             ///< booleen pour sauvegarder le vecteur de deplacement � chaque pas de temps
    //bool alloc_quantites;           ///< booleen indiquant s'il est necessaire d'allouer les quantites et de realiser une etape prelocale (defaut 1) PLUS UTILISE POUR LE MOMENT

    /** Parametres pour les directions de recherche
    \brief Type de direction de recherche (utilise kfact)

    Plusieurs choix sont possibles selon le nom donne � ce param�tre : 
    - "scalaire_donne" : un param�tre scalaire identique pour toutes les interfaces donne par kfact
    - "scalaire_donne_CL" : un param�tre scalaire identique pour toutes les interfaces avec modification pour les interfaces de CL ou de type contact
    - "scalaire_auto" : un param�tre scalaire identique pour toutes les interfaces obtenu en prenant le module d'young de la sous-structure voisine divise par une longueur caracteristique de l'interface consideree. Ce facteur est multiplie par kfact
    - "scalaire_auto_CL" : idem en modificant les valeurs pour les interfaces � CL ou de type contact

    Dans chaque cas, on calcule une direction normale et tangentielle (identiques pour la plupart des cas) : \f$ k_n \f$ et \f$ k_t \f$
    Les modifications pour les CL ou les interfaces de type contact sont les suivantes :
    - pour les interfaces de type effort, on multiplie les scalaires par 1/1000
    - pour les interfaces de type deplacement  on multiplie par 1000 les scalaires \f$ k_t \f$ et \f$ k_n \f$
    - pour les interfaces de type symetrie ou deplacement normal donne, on multiplie par 1000 \f$ k_n \f$ et par 1/1000 \f$ k_t \f$
    - pour les interfaces de type Contact, on multiplie par la coefficient de frottement \f$ k_t = f * k_n \f$   
    */
    Sc2String ktype;        ///< type de direction de recherche
    Scalar kfact;           ///< facteur multiplicatif de la direction de recherche 
    unsigned copydirection; ///< copie de la direction de recherche de part et d'autre d'une interface
    
    //methodes===============================================================================================
    LatinData(); ///< Constructeur : valeurs par defaut pour les parametres
    void read_data_user(const Metil::DataUser &data_user); ///< Charge les donnees depuis le DataUser
};

#endif //LATINDATA_H

