#ifndef SST_H
#define SST_H

#include "Process.h"
//definition automatique des formulations ancetres et filles pour chaque formulation et element donnes dans SConstruct
#include "../../build/problem_pb_elast/problem.h"
// definition du maillage sur les bords de la sous-structure
#include "../MAILLAGE/meshcaracinter.h"
#include "meshmulti.h"
#include "SstCarac_InterCarac.h"

#include "../../LMT/include/util/solveLDL.h"
#include "../../LMT/include/containers/vec.h"
using namespace LMT;


//**********************************
// Sous structure
//*********************************
/** \ingroup Sous_structures
\brief Classe Sous-structure parametrable
*/
struct Sst
{
    typedef Meshmulti<Mesh_carac_pb_elast<TYPEREEL,DIM> > TMESH;    ///< type de maillage pour la sous-structure
    typedef Mesh<Meshcaracinter<DIM,DIM-1> > TMESHedge;             ///< type de maillage pour le bord des sous-structures
    typedef Mat<TYPEREEL, Sym<>, SparseCholMod > TMATS;             ///< type de matrice sparse (pour le solveur CholMod)
    typedef FormulationAncestor<TYPEREEL> TF;                       ///< formulation generique choisie par le type de materiau et le type de resolution etudiee
    typedef Vec<TYPEREEL,DIM> Pvec;                                 ///< Type des points
    
    // donnees geometriques
    Pvec G;             ///< centre de gravite
    TYPEREEL measure;   ///< mesure de la sst
    int num;            ///< numero de la sous-structure
    int id;             ///< id de la sous-structure
    int num_proc;       ///< numero du processeur qui traite la sous-structure
    Vec<int> vois;      ///< sous-structures voisines (-1 pour le bord)
    Vec<Pvec,2 > box;   ///< boite incluant le maillage de la sst
    
    // comportement de la SST
    unsigned id_material;               ///< id du materiau dans le fichier json
    unsigned typmat;                    ///< type de materiau par sst
    bool plastique;                     ///< indique si le materiau doit etre plastique
    SstCarac matprop;                   ///< pointeur vers les caracteristiques du materiau
    Problem_pb_elast<TYPEREEL,DIM> pb;  ///< problème étudié (permet de définir toutes les formulations disponibles)
    TF *f;                              ///< pointeur vers une formulation

    TMESH mesh;     ///< maillage de la sst : contient les champs indiques dans Mesh_carac_pb_elast

    ///Structure contenant les informations sur chaque bord de la sous-structure
    struct Edge
    {
        Edge(): mesh(NULL) {}

        int datanum;                ///< cote correspondant de l'interface (0 ou 1)
        int internum;               ///< numero de l'interface
        Vec<unsigned> repddledge;   ///< reperage des ddls de bords dans les ddls du maillage
        Vec<unsigned> repLE;        ///< reperage des ddls macro dans l'operateur homogeneise     
        TMESHedge *mesh;            ///< maillage de bord
        Vec<Pvec,2> box;            ///< boite englobant le cote
        Pvec G;                     ///< centre de gravite du cote
    };

    Vec<Edge> edge;     ///< Bords de la sous-structure

    // operateurs
    unsigned nb_macro;                      ///< nbre de fcts de base macro par sst (somme des fonctions macro des interfaces autour de la sst)
    unsigned nb_macro_espace;               ///< nbre de fcts de base macro spatiale par sst
    Mat <TYPEREEL , Gen<>, Dense<> > LE;    ///< matrice de comportement homogeneise
    TMATS *K;                               ///< matrice de rigidite + directions de recherche (utile pour le solver Cholmod)

    Vec<TYPEREEL> fvol ;   ///< second membre provenant des efforts ou quantites chapeaux volumiques par sous-structure
    LDL_solver l;          ///< solver sparse direct plus rapide que skyline

    /// Structure temporelle contenant differents vecteurs
    struct Time{
        Vec<TYPEREEL> q;                              ///< vecteur des deplacements aux noeuds
        Vec<TYPEREEL> p;                              ///< vecteur de la plasticite cumulee aux points de Gauss
        Vec<TYPEREEL> R_p;                            ///< vecteur de l'ecrouissage aux points de Gauss
        Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > epsilon_p;  ///< vecteur des deformations plastiques aux points de Gauss
        Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > sigma;      ///< vecteur des contraintes aux points de Gauss
        Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > X_p;        ///< vecteur des centres du domaine elastique de chaque elements
        //Vec<TYPEREEL,2> erreur_plasticite;  ///< vecteur des erreurs en plasticite cumulee calculees aux points de Gauss PRISE EN COMPTE ERREUR_PLASTICITE NON IMPLEMENTEE
        
        void allocations(unsigned nbddl,unsigned nbelem,const Process &process, Sst &S);   ///< Alloue les vecteurs de Time
        void display_all_data();
    };
    
    Vec<Time> t;        ///< Vecteur contenant pour chaque piquet de temps les solutions
    Vec<Time> t_post;   ///< Vecteurs piquet de temps pour sauvegarder les donnees pour chaque piquet de temps (en incremental, non utile en latin)
    
    Vec<TYPEREEL> Fadd;    ///< vecteur effort macro sur le bord d'une SST Fadd

#warning "on suppose qu'en 2D, la connectivite vaut 2, en 3D elle vaut 3  pour les elements de peau"
#if DIM==2
    static const   int nb_nodes_by_element=2;
#else
    static const   int nb_nodes_by_element=3; 
#endif  

    BasicVec<BasicVec<TYPEREEL>,DIM> nodes;                         ///< coordonnées des noeuds de peau d'une sst pour la sortie hdf
    BasicVec<BasicVec<TYPEREEL>,DIM > dep_nodes_skin;               ///< deplacements aux noeuds de peau d'une sst
    BasicVec<BasicVec<TYPEREEL>,DIM > dep_nodes;                    ///< deplacements aux noeuds de peau d'une sst
    BasicVec<BasicVec<int> > mesh_connectivities_skin;              ///< connectivites du maillage de peau d'une sst pour la sortie hdf (tient compte de la numérotation globale des noeuds)
    BasicVec<BasicVec<int> > mesh_connectivities;                   ///< connectivites du maillage d'une sst pour la sortie hdf (tient compte de la numérotation globale des noeuds)
    BasicVec<BasicVec< TYPEREEL > , DIM*(DIM+1)/2 > sigma_skin;     ///< vecteur contrainte du maillage de peau d'une sst pour la sortie hdf 
    BasicVec<BasicVec< TYPEREEL > , DIM*(DIM+1)/2 > epsilon_skin;   ///< vecteur deformation du maillage de peau d'une sst pour la sortie hdf 
    BasicVec< TYPEREEL > sigma_mises_skin;                          ///< vecteur contrainte de von mises du maillage de peau d'une sst pour la sortie hdf 
    BasicVec<BasicVec< TYPEREEL > , DIM*(DIM+1)/2 > sigma;          ///< vecteur contrainte du maillage de peau d'une sst pour la sortie hdf 
    BasicVec<BasicVec< TYPEREEL > , DIM*(DIM+1)/2 > epsilon;        ///< vecteur deformation du maillage de peau d'une sst pour la sortie hdf 
    BasicVec< TYPEREEL > sigma_mises;   ///< vecteur contrainte de von mises du maillage de peau d'une sst pour la sortie hdf 
    BasicVec< int > num_processor;      ///< vecteur numéro du processeur d'une sst pour la sortie hdf 
    BasicVec< int > num_group;          ///< vecteur numéro des sst pour la sortie hdf 
    BasicVec< int > material;           ///< vecteur numéro des materiaux des sst pour la sortie hdf 
    Vec<int,4> nb_elements_with_type;   ///< utilisé pour connaitre le nombre d'elements d'un type donné pour une SST : Triangle, Quadrilateral, Tetrahedron, Hexahedron
    int nb_nodes_by_element_sst, nb_nodes_by_element_sst_skin, pattern_id;
    Sc2String type_elements_sst, type_elements_sst_skin;
  
    Sst();     ///< constructeur de la formulation pour la sous-structure
    ~Sst();    ///< destructeur de la SST
    
    void free();  ///< Libere tous les pointeurs de la Sst
    
    /// Trouver une sst à partir de son id----------------------------------------------
    static Sst* find_sst(Vec<Sst> &S,int id_);
    
    /// Trouver l'index d'une sst à partir de son id -----------------------------------
    static int find_index_sst(Vec<Sst> &S, int id_);
    
    /// Assigne les proprietes materiau au maillage associe a la Sst
    void assign_material_on_element(const DataUser &data_user);
};

#endif // SST_H

