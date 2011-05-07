#ifndef MATERIALS_H
#define MATERIALS_H

// Ajout class PARAM_DAMAGE_SST contenant les grandeurs materiaux associées a l'endommagement
#include "definition_PARAM_DAMAGE_SST.h"

using namespace LMT;
using namespace std;


//*******************************************
// Caracteristiques materielles des SST
//*******************************************
/** \ingroup Materiaux
\brief Classe parametrable definissant les proprietes materielles par sous-structure
*/
template<unsigned dim_, class T_> struct SstCarac
{
  SstCarac(){dt=1.;resolution=1;}
  typedef  T_ T; ///< type de flottant
  static const unsigned dim=dim_; ///< dimension 2 ou 3
  string type;     ///< type de formulation : isotrope, orthotrope, orthotrope endommageable
  bool resolution; ///< type de resolution contrainte_plane (1) ou deformation_plane (0) : utilise en 2d
  Vec<double> coef; ///< coefficients materiau : isotrope E, v, orthotrope : E1, E2, E3, v12, v13, v23, G12, G13, G23
  Vec<Vec<double,3>,2 > direction; ///< direction pour les materiaux orthotropes
  Vec<double> coefth; ///< coefficients thermiques : isotrope : alpha, orthotrope : alpha_1, alpha_2, alpha_3
  double caract;   ///< caracteristiques (epaisseur par exemple), non pris en compte
  double dt; ///< pas de temps lu uniquement pour la quasistatique (obtenu a partir de process.temps->dt)
  Vec<double,dim> f_vol;///force volumique constante par sst
  //Vec<double,3> epshorsplan; // deformation horsplan pour les def planes generalisees
  ///parametres pour l'endommagement
  PARAM_DAMAGE_SST<dim,T> param_damage;
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

template<unsigned dim_, class T_> struct InterCarac
{
  typedef  T_ T; ///< type de flottant
  static const unsigned dim=dim_; ///< dimension 2 ou 3
  typedef Vec<T_,dim_> Pvec; ///< type des points
  
  //selection des interfaces pour application des parametres materiau : 2 possibilites : boite ou numero des sst voisines
  Vec<Pvec,2> box; ///< boite dans laquelle les interfaces possedant cette propriete doivent se trouver
  Vec<unsigned,2> num_sst; ///< numero des sous-structures adjacentes
  T coeffrottement; ///< coefficient de frottement 
  string jeu ; ///< fonction donnant le jeu en fonction des variables d'espace par une fonction analytique
  unsigned nbpastempsimpos;

  unsigned id;  ///< numero d'identification de cette caracteristique
  string name;  ///< nom de la caracteristique d'interface (informatif)
  T Gcrit;  ///< valeur de taux de restitution critique pour les interfaces discretes
  string type; ///< type d'interfaces : contact_box, contact_sst, contact_jeu_box, contact_jeu_sst, contact_jeu_physique, jeu_impose_sst, jeu_impose_box, cohesive, discrete
  string comp; ///< comportement des interfaces incluses
  
  ///parametres d'endommagement pour les interfaces cohesives
  PARAM_DAMAGE_INTER param_damage;

};

#endif //MATERIALS_H

