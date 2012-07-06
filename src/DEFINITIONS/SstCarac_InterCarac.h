#ifndef MATERIALS_H
#define MATERIALS_H

// Ajout class PARAM_DAMAGE_SST contenant les grandeurs materiaux associees a l'endommagement
#include "definition_PARAM_COMP_INTER.h"
#include "main_typedef.h"
#include "../PARAMETERS/Parameters.h"


//*******************************************
// Caracteristiques materielles des SST
//*******************************************
/** \ingroup Materiaux
\brief Classe parametrable definissant les proprietes materielles par sous-structure
*/
struct SstCarac
{
    static ParameterGroup sst_materials_parameters;
    
    /// Attributs communs
    int id;                         /// identite du materiaux dans data_user
    int type_num;                   /// numero d'identite du comportement materiaux : 0=isotrope elastique
    Sc2String type;                 /// type de formulation : isotrope, orthotrope, orthotrope endommageable
    Sc2String comp;                 /// type de comportement : elastique, endommageable, plastique, mesomodele...
    bool resolution;                /// type de resolution contrainte_plane (1) ou deformation_plane (0) : utilise en 2d
    Point v1,v2;        /// direction pour les materiaux orthotropes
    UserParameter density;          /// densite du materiaux
    Point f_vol;        /// Champs de force volumique constant
    LMT::Vec<Sc2String,DIM> f_vol_e;     /// Champs de force volumique par element
    Scalar dt;                    /// pas de temps lu uniquement pour la quasistatique (obtenu a partir de process.temps->dt)
    
    /// Comportement elastique
    /// Isotrope
    UserParameter elastic_modulus;
    UserParameter poisson_ratio;
    /// Anisotrope
    UserParameter elastic_modulus_1;
    UserParameter elastic_modulus_2;
    UserParameter elastic_modulus_3;
    UserParameter poisson_ratio_12;
    UserParameter poisson_ratio_13;
    UserParameter poisson_ratio_23;
    UserParameter shear_modulus_12;
    UserParameter shear_modulus_13;
    UserParameter shear_modulus_23;
    
    /// Comportement thermique
    Scalar deltaT;    /// Variable pour le stockage de la variation de temperature
    /// Isotrope
    UserParameter alpha;
    /// Anisotrope
    UserParameter alpha_1;
    UserParameter alpha_2;
    UserParameter alpha_3;
    
    /// Compoprtement plastique
    Sc2String type_plast;   /// Options decrivant l'evolution du domaine elastique (parfaite, isotrope, cinematique) et la contrainte equivalente (von_mises)
    /// Coefficients de la loi d'ecrouissage : R(p) = R0 + k_p * p ^ m_p
    UserParameter plast_ecrouissage_init;    /// Limite d'elasticite initiale
    UserParameter plast_ecrouissage_mult;    /// Coefficient multiplicateur de la loi d'ecrouissage
    UserParameter plast_ecrouissage_expo;    /// Exposant de la loi d'ecrouissage
    UserParameter plast_cinematique_coef;    /// Coefficient multiplicateur du modele cinematique : dX = C*depsilon_p
    VoigtMatrix coeff_seq;     /// Coefficients multiplicateurs pour la fonction seuil (matrice de la forme bilineaire associee a la norme energetique)
                                                                /// COEFF_SEQ N'EST PAS UTILISE POUR L'INSTANT
    
    /// Comportement endommageable de type mesomodele de composite
    Sc2String type_endo;    /// Options decrivant l'endommagement A SPECIFIER!!!
    UserParameter Yo;
    UserParameter Yc;
    UserParameter Ycf;      /// Efforts limites de l'endommagement de la rupture
    Scalar dmax;            /// Maximum possible de l'endommagement avant rupture (pour eviter les divisions par 0)
    Scalar couplage;        /// Coefficient pour le calcul de la contrainte equivalente
    Scalar b_c;             /// Couplage entre micro-fissuration et decohesion fibres/matrice
    bool effet_retard;      /// Indique si un effet retard doit etre applique a l'endommagement
    Scalar a,tau_c;         /// Coefficients de l'effet_retard de l'endommagement
    
    /// Comportement visqueux
    UserParameter viscosite;     ///< Viscosite du materiau
    
    
    SstCarac();
    void read_data_user(int index,DataUser &data_user);
    static void prepareParameters();
    static void updateParameters();
    void affiche();
};

//**************************************************************************
// Caracteristiques des interfaces pour application de parametres 
//**************************************************************************
/** \ingroup Materiaux
\brief Classe parametrable contenant les parametres materiau associes a des interfaces donnees

Deux possibilites sont offertes pour selectionner les interfaces particulieres :
- Utilisation d'une boite definie par ses deux points extremes
- Utilisation des numeros des deux sous-structures adjacentes

*/

struct InterCarac
{
    static ParameterGroup inter_materials_parameters;
    
    unsigned id;                    /// numero d'identification de cette caracteristique
    unsigned type_num;                 /// numero d'identification du comportement d'interface
    Sc2String name;                    /// nom de la caracteristique d'interface (informatif)
    Sc2String type;                    /// type d'interfaces : contact_box, contact_sst, contact_jeu_box, contact_jeu_sst, contact_jeu_physique, jeu_impose_sst, jeu_impose_box, cohesive, discrete, contact_ep
    Sc2String comp;                    /// comportement des interfaces incluses
    LMT::Vec<unsigned,2> num_sst;   /// numero des sous-structures adjacentes
    bool degradable;                /// indique si l'interface est degradable
    
    UserParameter f_jeu ;               /// parametre donnant le jeu en fonction des variables d'espace par une fonction analytique
    UserParameter f_coeffrottement;     /// parametre representant le coefficient de frottement
    //Scalar coeffrottement;        /// coefficient de frottement sur l'interface
    UserParameter f_raideur;            /// parametre donnant la raideur d'une interface elastique
    UserParameter Gcrit;                /// valeur de taux de restitution critique pour les interfaces cassables
    
    UserParameter kn;       /// raideur normale (en traction si knc definie)
    UserParameter kt;       /// raideur tangentielle
    UserParameter knc;      /// raideur normale en compression
    
    UserParameter gamma;    /// coefficients pour le calcul des efforts associes a l'endommagement
    UserParameter alpha;    /// Y(t) = sup[tau < t]{ ( Y3(tau)^alpha + (gamma*Y1(tau))^alpha + (gamma*Y2(tau))^alpha )^(1/alpha) } si la normale est le vecteur indice 3
    
    UserParameter n;        /// coefficients pour le calcul de l'endommagement:
    UserParameter Yc;       /// d = min( (n/(n+1))*abs(Y-Yo)/(Yc-Yo) , 1 )
    UserParameter Yo;       /// Yo et Yc sont les seuils limite de l'endommagement
    
    unsigned nbpastempsimpos;   /// Voir comportement "jeu impose"
    
    InterCarac();
    void allocate();
    void free();
    void read_data_user(int index,DataUser &data_user);
    static void prepareParameters();
    static void updateParameters();
    void affiche();
};

/*
struct InterCarac
{
    typedef Vec<Scalar,DIM> Pvec; ///< type des points
    
    ///selection des interfaces pour application des parametres materiau : 2 possibilites : boite ou numero des sst voisines
    Vec<Pvec,2> box;            /// boite dans laquelle les interfaces possedant cette propriete doivent se trouver
    Vec<unsigned,2> num_sst;    /// numero des sous-structures adjacentes
    Scalar coeffrottement;    /// coefficient de frottement 
    String jeu ;             /// fonction donnant le jeu en fonction des variables d'espace par une fonction analytique
    String f_coeffrottement; /// fonction analytique donnant le coeficient de frottement
    String f_R;              /// fonction analytique donnant la raideur d'une interface elastique
    
    unsigned nbpastempsimpos;   /// Voir comportement "jeu impose"
    
    unsigned id;       /// numero d'identification de cette caracteristique
    String name;    /// nom de la caracteristique d'interface (informatif)
    Scalar Gcrit;    /// valeur de taux de restitution critique pour les interfaces discretes
    String type;    /// type d'interfaces : contact_box, contact_sst, contact_jeu_box, contact_jeu_sst, contact_jeu_physique, jeu_impose_sst, jeu_impose_box, cohesive, discrete, contact_ep
    String comp;    /// comportement des interfaces incluses
    
    ///parametres d'endommagement pour les interfaces cohesives
    PARAM_DAMAGE_INTER param_damage;
    
    
    void affiche(){
        std::cout << "----------------------------------------------------------------------"<< std::endl;
        std::cout << "id = " << id << std::endl;
        std::cout << "coeffrottement = " << coeffrottement << std::endl;
        std::cout << "jeu = " << jeu << std::endl;
        std::cout << "Gcrit = " << Gcrit << std::endl;
        std::cout << "type = " << type << std::endl;
        std::cout << "comp = " << comp << std::endl;
        std::cout << "----------------------------------------------------------------------"<< std::endl;
    }
    
};
//*/
#endif //MATERIALS_H

