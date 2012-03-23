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
    SstCarac(){dt=1.;resolution=1;}
    
    int id;             ///< identite du materiaux dans data_user
    int type_num;       ///< numero d'identité du comportement materiaux : 0=isotrope elastique
    Sc2String type;        ///< type de formulation : isotrope, orthotrope, orthotrope endommageable, mesomodele
    Sc2String comp;        ///< type de comportement : elastique, endommageable, plastique...
    bool resolution;    ///< type de resolution contrainte_plane (1) ou deformation_plane (0) : utilise en 2d
    Vec<TYPEREEL,DIM> v1,v2;    ///< direction pour les materiaux orthotropes
    TYPEREEL density;      ///< densite du materiaux
    TYPEREEL elastic_modulus,poisson_ratio, elastic_modulus_1,elastic_modulus_2,elastic_modulus_3,poisson_ratio_12,poisson_ratio_13,poisson_ratio_23,shear_modulus_12,shear_modulus_13,shear_modulus_23;
    TYPEREEL alpha, alpha_1,alpha_2,alpha_3, deltaT;    ///< Coefficients de la thermique
    TYPEREEL k_p,m_p,R0,coefvm_composite;          ///< Coefficients de la plasticite
    TYPEREEL Yo,Yc,Ycf,dmax,b_c,a,tau_c;           ///< Coefficients de l'endommagement
    bool effet_retard;
    Vec<TYPEREEL,DIM> f_vol;       ///< champs de force volumique constant
    Vec<Sc2String,DIM> f_vol_e;    ///< champs de force volumique par element
    TYPEREEL dt;           ///< pas de temps lu uniquement pour la quasistatique (obtenu a partir de process.temps->dt)
    TYPEREEL viscosite;
    
    void affiche(){
        std::cout << std::endl << std::endl;
        std::cout << "******************SSTCARAC_DEBUG*********************" << std::endl;
        std::cout << "density : " << density << std::endl;
        if (type == "isotrope"){
            std::cout << "elastic_modulus : " << elastic_modulus << std::endl;
            std::cout << "poisson_ratio : " << poisson_ratio << std::endl;
            std::cout << "alpha : " << alpha << std::endl;
        }
        if (type == "orthotrope" or type == "mesomodele"){
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
        if (type == "mesomodele"){
            std::cout << "k_p      : " << k_p << std::endl;
            std::cout << "m_p      : " << m_p << std::endl;
            std::cout << "R0       : " << R0 << std::endl;
            std::cout << "couplage : " << coefvm_composite << std::endl;
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

