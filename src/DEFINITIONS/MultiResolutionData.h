#ifndef MULTIRESOLUTIONDATA_H
#define MULTIRESOLUTIONDATA_H

#include "main_typedef.h"
#include "../PARAMETERS/Parameters.h"
#include "../COMPUTE/DataUser.h"
#include <Metil/BasicVec.h>

struct MultiResolutionData{
    
    struct ParameterData
    {
        enum Type{Function,List};
        Sc2String name;                 /// Le nom du parametre
        Type type;                      /// Type de source (liste ou fonction)
        BasicVec<Sc2String> values;     /// La/Les valeur(s) du parametres
        UserParameter* user_parameter;  /// Pointeur sur le UserParameter associe
    };
    
    Vec<ParameterData> parameters_data; /// Les donnees "brutes" des parametres
    Sc2String type;                     /// Type de calcul ("off","sequential" ...)
    MainParameter m;                    /// Parametre de controle representant le numero du calcul parametrique
    ParameterGroup parameters;          /// Groupe des parametres de multi-resolution
    Vec<Id> id_materials;               /// Liste des Id des materiaux utilisant au moins un des parametres
    Vec<Id> id_links;                   /// Liste des Id des interfaces utilisant au moins un des parametres
    Vec<Id> id_CL;                      /// Liste des Id des CL utilisant au moins un des parametres
    Vec<Id> id_CLvolume;                /// Liste des Id des CL volumiques utilisant au moins un des parametres
    unsigned nb_calculs;                /// Nombre total de calculs a realiser
    unsigned calcul_cur;                /// Indice du calcul courant
    
    MultiResolutionData();
    void free();
    void read_data_user(const DataUser &data_user);   /// Lecture du DataUser
    void prepareParameters();
    void updateParameters();
    
    void init();        /// Initialiser les variables
    void next();        /// Passer au calcul suivant
    bool has_next();    /// Indique s'il reste des calculs a realiser
    
    void affiche();
};

#endif  //MULTIRESOLUTIONDATA_H
