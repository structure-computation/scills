#ifndef MATERIALS_H
#define MATERIALS_H

// Ajout class PARAM_DAMAGE_SST contenant les grandeurs materiaux associées a l'endommagement
#include "definition_PARAM_COMP_INTER.h"
#include <boost/concept_check.hpp>
#include "../UTILS/Sc2String.h"
using namespace LMT;


//*******************************************
// Caracteristiques materielles des SST
//*******************************************
/** \ingroup Materiaux
\brief Classe parametrable definissant les proprietes materielles par sous-structure
*/
struct SstCarac
{
    SstCarac(){
        dt=1.0;
        resolution=1;
        /// Initialisation des grandeurs a -1 pour detecter les erreurs de chargement
        density = -1;
        elastic_modulus = -1;
        poisson_ratio = -1;
        elastic_modulus_1 = -1;
        elastic_modulus_2 = -1;
        elastic_modulus_3 = -1;
        poisson_ratio_12 = -1;
        poisson_ratio_13 = -1;
        poisson_ratio_23 = -1;
        shear_modulus_12 = -1;
        shear_modulus_13 = -1;
        shear_modulus_23 = -1;
        alpha = -1;
        alpha_1 = -1;
        alpha_2 = -1;
        alpha_3 = -1;
        deltaT = -1;
        k_p = -1;
        m_p = -1;
        R0 = -1;
        Yo = -1;
        Yc = -1;
        Ycf = -1;
        dmax = -1;
        b_c = -1;
        a = -1;
        tau_c = -1;
        viscosite = -1;
        
        ///Initialisation des coefficients de la contrainte equivalente a Von Mises standart
        for(int i = 0; i < DIM*(DIM+1)/2; i++)
            coeff_seq[i].set(0.0);
        if(DIM == 2){
            ;
        } else if(DIM == 3){
            coeff_seq[0][0] = 1;
            coeff_seq[1][1] = 1;
            coeff_seq[2][2] = 1;
            coeff_seq[0][1] = -1;
            coeff_seq[0][2] = -1;
            coeff_seq[1][2] = -1;
            coeff_seq[3][3] = 3;
            coeff_seq[4][4] = 3;
            coeff_seq[5][5] = 3;
        }
    }
    
    /// Attributs communs
    int id;                         ///< identite du materiaux dans data_user
    int type_num;                   ///< numero d'identité du comportement materiaux : 0=isotrope elastique
    Sc2String type;                 ///< type de formulation : isotrope, orthotrope, orthotrope endommageable
    Sc2String comp;                 ///< type de comportement : elastique, endommageable, plastique, mesomodele...
    bool resolution;                ///< type de resolution contrainte_plane (1) ou deformation_plane (0) : utilise en 2d
    Vec<TYPEREEL,DIM> v1,v2;        ///< direction pour les materiaux orthotropes
    TYPEREEL density;               ///< densite du materiaux
    Vec<TYPEREEL,DIM> f_vol;        ///< Champs de force volumique constant
    Vec<Sc2String,DIM> f_vol_e;     ///< Champs de force volumique par element
    TYPEREEL dt;                    ///< pas de temps lu uniquement pour la quasistatique (obtenu a partir de process.temps->dt)
    
    /// Comportement elastique
    TYPEREEL elastic_modulus,poisson_ratio; /// Isotrope
    TYPEREEL elastic_modulus_1,elastic_modulus_2,elastic_modulus_3,poisson_ratio_12,poisson_ratio_13,poisson_ratio_23,shear_modulus_12,shear_modulus_13,shear_modulus_23; /// Anisotrope
    
    /// Comportement thermique
    TYPEREEL alpha;                     ///< Isotrope
    TYPEREEL alpha_1,alpha_2,alpha_3;   ///< Anisotrope
    TYPEREEL deltaT;                    ///< Variable pour le stockage de la variation de temperature
    
    /// Compoprtement plastique
    TYPEREEL R0;        ///< Contrainte limite d'elasticite
    TYPEREEL k_p,m_p;   ///< Coefficients de la loi d'ecrouissage
    TYPEREEL coefvm_composite;                                    ///< Coefficient de couplage pour la contrainte equivalente (cf mesomodele de composite)
    Vec<Vec<TYPEREEL,DIM*(DIM+1)/2>,DIM*(DIM+1)/2> coeff_seq;     ///< Coefficients multiplicateurs pour la fonction seuil (matrice de la forme bilineaire associee a la norme energetique)
    
    /// Comportement endommageable de type mesomodele de composite
    TYPEREEL Yo,Yc,Ycf;     ///< Efforts limites de l'endommagement de la rupture
    TYPEREEL dmax;          ///< Maximum possible de l'endommagement avant rupture (pour eviter les divisions par 0)
    TYPEREEL b_c;           ///< Couplage entre micro-fissuration et decohesion fibres/matrice
    bool effet_retard;      ///< Indique si un effet retard doit etre applique a l'endommagement
    TYPEREEL a,tau_c;       ///< Coefficients de l'effet_retard de l'endommagement
    
    /// Comportement visqueux
    TYPEREEL viscosite;     ///< Viscosite du materiau
    
    void affiche(){
        std::cout << std::endl << std::endl;
        std::cout << "******************SSTCARAC_DEBUG*********************" << std::endl;
        std::cout << "density : " << density << std::endl;
        if (type == "isotrope"){
            std::cout << "elastic_modulus : " << elastic_modulus << std::endl;
            std::cout << "poisson_ratio : " << poisson_ratio << std::endl;
            std::cout << "alpha : " << alpha << std::endl;
        } else if (type == "orthotrope"){
            std::cout << "elastic_modulus_1 : " << elastic_modulus_1 << std::endl;
            std::cout << "elastic_modulus_2 : " << elastic_modulus_2 << std::endl;
            std::cout << "elastic_modulus_3 : " << elastic_modulus_3 << std::endl;
            std::cout << "poisson_ratio_12 : " << poisson_ratio_12 << std::endl;
            std::cout << "poisson_ratio_13 : " << poisson_ratio_13 << std::endl;
            std::cout << "poisson_ratio_23 : " << poisson_ratio_23 << std::endl;
            std::cout << "shear_modulus_12 : " << shear_modulus_12 << std::endl;
            std::cout << "shear_modulus_13 : " << shear_modulus_13 << std::endl;
            std::cout << "shear_modulus_23 : " << shear_modulus_23 << std::endl;
            std::cout << "alpha_1 : " << alpha_1 << std::endl;
            std::cout << "alpha_2 : " << alpha_2 << std::endl;
            std::cout << "alpha_3 : " << alpha_3 << std::endl;
        }
        std::cout << "deltaT  : " << deltaT << std::endl;
        if (comp == "plastique" or comp == "mesomodele"){
            std::cout << "k_p      : " << k_p << std::endl;
            std::cout << "m_p      : " << m_p << std::endl;
            std::cout << "R0       : " << R0 << std::endl;
            std::cout << "couplage : " << coefvm_composite << std::endl;
        }
        if (comp == "mesomodele"){
            std::cout << "Yo           : " << Yo << std::endl;
            std::cout << "Yc           : " << Yc << std::endl;
            std::cout << "Ycf          : " << Ycf << std::endl;
            std::cout << "dmax         : " << dmax << std::endl;
            std::cout << "b_c          : " << b_c << std::endl;
            std::cout << "effet_retard : " << effet_retard << std::endl;
            std::cout << "a            : " << a << std::endl;
            std::cout << "tau_c        : " << tau_c << std::endl;
        }
    }
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
    typedef Vec<TYPEREEL,DIM> Pvec; ///< type des points
    
    //selection des interfaces pour application des parametres materiau : 2 possibilites : boite ou numero des sst voisines
    Vec<Pvec,2> box; ///< boite dans laquelle les interfaces possedant cette propriete doivent se trouver
    Vec<unsigned,2> num_sst; ///< numero des sous-structures adjacentes
    TYPEREEL coeffrottement; ///< coefficient de frottement 
    Sc2String jeu ; ///< fonction donnant le jeu en fonction des variables d'espace par une fonction analytique
    Sc2String f_coeffrottement; ///< fonction analytique donnant le coeficient de frottement
    Sc2String f_R; ///< fonction analytique donnant la raideur d'une interface elastique
    
    unsigned nbpastempsimpos;

    unsigned id;    ///< numero d'identification de cette caracteristique
    Sc2String name;    ///< nom de la caracteristique d'interface (informatif)
    TYPEREEL Gcrit;     ///< valeur de taux de restitution critique pour les interfaces discretes
    Sc2String type;    ///< type d'interfaces : contact_box, contact_sst, contact_jeu_box, contact_jeu_sst, contact_jeu_physique, jeu_impose_sst, jeu_impose_box, cohesive, discrete, contact_ep
    Sc2String comp;    ///< comportement des interfaces incluses
    
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

#endif //MATERIALS_H

