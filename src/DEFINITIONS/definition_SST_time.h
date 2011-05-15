#ifndef SST_H
#define SST_H

//definition automatique des formulations ancetres et filles pour chaque formulation et element donnes dans SConstruct
#include "problem_pb_elast/problem.h"
// definition du maillage sur les bords de la sous-structure
#include "meshcaracinter.h"
#include "meshmulti.h"

#include "util/solveLDL.h"
using namespace LMT;
using namespace std;

template<unsigned dim,class T> struct PARAM_DAMAGE_SST;

//**********************************
// Sous structure
//*********************************
/** \ingroup Sous_structures
\brief Classe Sous-structure parametrable

Paramètres de la classe :.
dim_ dimension du probleme 2 ou 3d.
TT_ type de flottant.
*/
template<unsigned dim_, class TT_> struct Sst
{

  Sst() : pb(*mesh.m,true) {} ///< constructeur de la formulation pour la sous-structure
  
  void free(){///destructeur de la SST - on ne peut pas libérer les entiers et flottants
    vois.free();
    //box.free();
    //pb.free();
    //cout << mesh.node_list.size() << endl;
    //cout << mesh.elem_list.size() << endl;

/*    for( unsigned i=0;i<edge.size() ;i++ ){
        delete edge[i].mesh;
#ifdef PRINT_ALLOC
    total_allocated[ typeid(typename Sst::TMESHedge).name() ] -= sizeof(typename Sst::TMESHedge);
#endif
    }*/
    edge.free();
    //cout << LE.data.size() << endl;
    LE.free();
    //delete f;
    //delete K;//K est libéré apres factorisation de toute facon et non calculer où y a pas besoin
    fvol.free();
    //l.free();//non defini
    
    //t.free();
    //t_post.free();
    
    //delete param_damage;
  }

  static const unsigned dim = dim_; ///< variable dim accessible de l'exterieur, constante, obtenu a partir du param dim_ de Sst<2,double>
  // types connus de l'exterieur
  typedef  TT_ T; ///< Type des flottants
  typedef  Meshmulti<Mesh_carac_pb_elast<T,dim> > TMESH; ///< type de maillage pour la sous-structure
  typedef  Mesh<Meshcaracinter<dim,dim-1> > TMESHedge; ///< type de maillage pour le bord des sous-structures
  typedef  Mat<T, Sym<>, SparseCholMod > TMATS; ///< type de matrice sparse (pour le solveur CholMod)
  //typedef  Formulation<TM,elasticity_isotropy> TF;
  typedef FormulationAncestor<T> TF; ///< formulation generique choisie par le type de materiau et le type de resolution etudiee
  typedef Vec<TT_,dim_> Pvec; ///< Type des points
  
 // donnees geometriques
  Pvec G; ///< centre de gravite
  T measure; ///< mesure de la sst
  int num; ///< numero de la sous-structure
  int num_proc; ///< numero du processeur qui traite la sous-structure
  Vec<int> vois; ///< sous-structures voisines (-1 pour le bord)
  Vec<Pvec,2 > box; ///< boite incluant le maillage de la sst
  
 // comportement de la SST
  unsigned typmat; ///< type de materiau par sst
  Problem_pb_elast<T,dim> pb; ///< problème étudié (permet de définir toutes les formulations disponibles)
  TF *f; ///< pointeur vers une formulation

  TMESH mesh; ///< maillage de la sst : contient les champs indiques dans Mesh_carac_pb_elast

///Structure contenant les informations sur chaque bord de la sous-structure
  struct Edge
  {
      Edge(){
          mesh=NULL;
/*          mesh.sub_mesh(Number<1>()).elem_list.change_hash_size( mesh,1);
          mesh.sub_mesh(Number<2>()).elem_list.change_hash_size( mesh,1);
          mesh.sub_mesh(Number<0>()).elem_list.change_hash_size( mesh,1);
          mesh.elem_list.change_hash_size( mesh,1);*/
      }

     int datanum; ///< cote correspondant de l'interface (0 ou 1)
     int internum; ///< numero de l'interface
     Vec<unsigned> repddledge; ///< reperage des ddls de bords dans les ddls du maillage
     Vec<unsigned> repLE; ///< reperage des ddls macro dans l'operateur homogeneise     
     TMESHedge *mesh; ///< maillage de bord
     Vec<Pvec,2> box; ///< boite englobant le cote
     Pvec G; ///< centre de gravite du cote
  };

  Vec<Edge> edge; ///< Bords de la sous-structure

 // operateurs
  unsigned nb_macro; ///< nbre de fcts de base macro par sst (somme des fonctions macro des interfaces autour de la sst)
  unsigned nb_macro_espace; ///< nbre de fcts de base macro spatiale par sst
  Mat <T , Gen<>, Dense<> > LE; ///< matrice de comportement homogeneise
  TMATS *K; ///< matrice de rigidite + directions de recherche (utile pour le solver Cholmod)

  Vec<T> fvol ; ///< second membre provenant des efforts ou quantites chapeaux volumiques par sous-structure
  LDL_solver l; ///< solver sparse direct plus rapide que skyline

  /// Structure temporelle contenant differents vecteurs
  struct Time{
  Vec<T> q; ///< vecteurs solution q
  void allocations(unsigned nbnode,Param &process); ///< defini dans allocate.h
  };
  
  Vec<Time> t;///< Vecteur contenant pour chaque piquet de temps les solutions
  Vec<Time> t_post; ///< Vecteurs piquet de temps pour sauvegarder les donnees pour chaque piquet de temps (en incremental, non utile en latin)
  
  Vec<T> Fadd;///vecteur effort macro sur le bord d'une SST Fadd
  //ajout mesomodele
  PARAM_DAMAGE_SST<dim,T> *param_damage; ///< paramètres d'endommagement pour le mésomodèle
  
};

#endif // SST_H

