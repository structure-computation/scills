#ifndef MATERIALS_H
#define MATERIALS_H

// Ajout class PARAM_DAMAGE_SST contenant les grandeurs materiaux associées a l'endommagement
#include "definition_PARAM_COMP_INTER.h"
#include <boost/concept_check.hpp>

using namespace LMT;
using namespace std;


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
    string type;        ///< type de formulation : isotrope, orthotrope, orthotrope endommageable, mesomodele
    string comp;        ///< type de comportement : elastique, endommageable, plastique...
    bool resolution;    ///< type de resolution contrainte_plane (1) ou deformation_plane (0) : utilise en 2d
    Vec<TYPEREEL,DIM> v1,v2;    ///< direction pour les materiaux orthotropes
    TYPEREEL density;      ///< densite du materiaux
    TYPEREEL elastic_modulus,poisson_ratio, elastic_modulus_1,elastic_modulus_2,elastic_modulus_3,poisson_ratio_12,poisson_ratio_13,poisson_ratio_23,shear_modulus_12,shear_modulus_13,shear_modulus_23;
    TYPEREEL alpha, alpha_1,alpha_2,alpha_3, deltaT;
    TYPEREEL k_p,m_p,R0,coefvm_composite;          ///< Coefficients de la plasticite
    TYPEREEL Yo,Yc,Ycf,dmax,b_c,a,tau_c;           ///< Coefficients de l'endommagement
    bool effet_retard;
    Vec<TYPEREEL,DIM> f_vol;       ///< champs de force volumique constant
    Vec<string,DIM> f_vol_e;    ///< champs de force volumique par element
    TYPEREEL dt;           ///< pas de temps lu uniquement pour la quasistatique (obtenu a partir de process.temps->dt)
    TYPEREEL viscosite;
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

template<unsigned dim_, class Reel_> struct InterCarac
{
    static const unsigned dim=dim_; ///< dimension 2 ou 3
    typedef Vec<Reel_,dim_> Pvec; ///< type des points
    
    //selection des interfaces pour application des parametres materiau : 2 possibilites : boite ou numero des sst voisines
    Vec<Pvec,2> box; ///< boite dans laquelle les interfaces possedant cette propriete doivent se trouver
    Vec<unsigned,2> num_sst; ///< numero des sous-structures adjacentes
    Reel_ coeffrottement; ///< coefficient de frottement 
    string jeu ; ///< fonction donnant le jeu en fonction des variables d'espace par une fonction analytique
    string f_coeffrottement; ///< fonction analytique donnant le coeficient de frottement
    string f_R; ///< fonction analytique donnant la raideur d'une interface elastique
    
    unsigned nbpastempsimpos;

    unsigned id;    ///< numero d'identification de cette caracteristique
    string name;    ///< nom de la caracteristique d'interface (informatif)
    Reel_ Gcrit;     ///< valeur de taux de restitution critique pour les interfaces discretes
    string type;    ///< type d'interfaces : contact_box, contact_sst, contact_jeu_box, contact_jeu_sst, contact_jeu_physique, jeu_impose_sst, jeu_impose_box, cohesive, discrete, contact_ep
    string comp;    ///< comportement des interfaces incluses
    
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

