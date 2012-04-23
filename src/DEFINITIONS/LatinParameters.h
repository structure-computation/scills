#ifndef PARAM_SC2LATIN_H
#define PARAM_SC2LATIN_H


#include "../UTILS/Sc2String.h"
#include "../COMPUTE/DataUser.h"
#include "../../LMT/include/containers/vec.h"
using namespace LMT;

/** \ingroup Parametres
\brief Param�tres associ�s � la strat�gie LATIN
*/
struct LatinParameters{
    //attributs==============================================================================================
    int iter;               ///< numero de l'it�ration courante
    int nbitermax;          ///< nombre d'it�rations max
    double mu;              ///< facteur de relaxation courant  (1 � la premi�re it�ration puis �gal au facteur de relaxation)
    double critere_erreur;              ///< valeur du crit�re d'erreur � atteindre
    double critere_erreur_diss;         ///< valeur du crit�re d'erreur en dissipation � atteindre (si = 0 on ne la calcule pas)
    double critere_erreur_auto_stop;    ///< valeur pour laquelle l'algorithme s'arrete si la variation d'erreur est inf�rieure 5 fois de suite � celle ci (si = 0 on ne fait pas d'arret automatique))
    Vec<double> error;  ///< erreurs au cours des iterations
    Sc2String type_error;  ///< type d'erreur : choix possible entre "ddr" : erreur bas�e sur les directions de recherche, "energie" : erreur bas�e sur un calcul d'�nergie par interface 
    bool list_error;    ///< listage de l'erreur a l'�cran au cours des it�rations
    double facteur_relaxation;  ///< facteur de relaxation entr� par l'utilisateur
    bool save_depl_SST;         ///< bool�en pour sauvegarder le vecteur de d�placement � chaque pas de temps
    bool alloc_quantites;       ///< bool�en indiquant s'il est necessaire d'allouer les quantites et de realiser une etape prelocale (defaut 1)

    /// parametres pour les directions de recherche
    double kfact;    ///< facteur multiplicatif de la direction de recherche 
    unsigned copydirection; ///< copie de la direction de recherche de part et d'autre d'une interface
    /** 
    \brief Type de direction de recherche (utilise kfact)

    Plusieurs choix sont possibles selon le nom donn� � ce param�tre : 
    - "scalaire_donne" : un param�tre scalaire identique pour toutes les interfaces donn� par kfact
    - "scalaire_donne_CL" : un param�tre scalaire identique pour toutes les interfaces avec modification pour les interfaces de CL ou de type contact
    - "scalaire_auto" : un param�tre scalaire identique pour toutes les interfaces obtenu en prenant le module d'young de la sous-structure voisine divis� par une longueur caract�ristique de l'interface consid�r�e. Ce facteur est multipli� par kfact
    - "scalaire_auto_CL" : idem en modificant les valeurs pour les interfaces � CL ou de type contact

    Dans chaque cas, on calcule une direction normale et tangentielle (identiques pour la plupart des cas) : \f$ k_n \f$ et \f$ k_t \f$
    Les modifications pour les CL ou les interfaces de type contact sont les suivantes :
    - pour les interfaces de type effort, on multiplie les scalaires par 1/1000
    - pour les interfaces de type d�placement  on multiplie par 1000 les scalaires \f$ k_t \f$ et \f$ k_n \f$
    - pour les interfaces de type sym�trie ou d�placement normal donn�, on multiplie par 1000 \f$ k_n \f$ et par 1/1000 \f$ k_t \f$
    - pour les interfaces de type Contact, on multiplie par la coefficient de frottement \f$ k_t = f * k_n \f$   
    */
    Sc2String ktype;
    
    //methodes===============================================================================================
    LatinParameters(); ///< Constructeur : valeurs par d�faut pour les parametres
    void read_data_user(DataUser &data_user); ///< Charge les donnees depuis le DataUser
};

#endif //SC2PARAM_LATIN_H
