#ifndef SST_H
#define SST_H

#include "Param.h"
//definition automatique des formulations ancetres et filles pour chaque formulation et element donnes dans SConstruct
#include "../../build/problem_pb_elast/problem.h"
// definition du maillage sur les bords de la sous-structure
#include "../MAILLAGE/meshcaracinter.h"
#include "meshmulti.h"
#include "SstCarac_InterCarac.h"

#include "../../LMT/include/util/solveLDL.h"
using namespace LMT;


//**********************************
// Sous structure
//*********************************
/** \ingroup Sous_structures
\brief Classe Sous-structure parametrable

Paramètres de la classe :.
dim_ dimension du probleme 2 ou 3d.
Reel_ type de flottant.
*/
struct Sst
{
    typedef  Meshmulti<Mesh_carac_pb_elast<TYPEREEL,DIM> > TMESH;   ///< type de maillage pour la sous-structure
    typedef  Mesh<Meshcaracinter<DIM,DIM-1> > TMESHedge;            ///< type de maillage pour le bord des sous-structures
    typedef  Mat<TYPEREEL, Sym<>, SparseCholMod > TMATS;            ///< type de matrice sparse (pour le solveur CholMod)
    typedef FormulationAncestor<TYPEREEL> TF;                       ///< formulation generique choisie par le type de materiau et le type de resolution etudiee
    typedef Vec<TYPEREEL,DIM> Pvec;                                 ///< Type des points
    
    // donnees geometriques
    Pvec G;             ///< centre de gravite
    TYPEREEL measure;      ///< mesure de la sst
    int num;            ///< numero de la sous-structure
    int id;             ///< id de la sous-structure
    int num_proc;       ///< numero du processeur qui traite la sous-structure
    Vec<int> vois;      ///< sous-structures voisines (-1 pour le bord)
    Vec<Pvec,2 > box;   ///< boite incluant le maillage de la sst
    
    // comportement de la SST
    unsigned id_material;               ///< id du materiau dans le fichier json
    unsigned typmat;                    ///< type de materiau par sst
    SstCarac matprop;       ///< pointeur vers les caracteristiques du materiau
    Problem_pb_elast<TYPEREEL,DIM> pb;     ///< problème étudié (permet de définir toutes les formulations disponibles)
    TF *f;                              ///< pointeur vers une formulation

    TMESH mesh;     ///< maillage de la sst : contient les champs indiques dans Mesh_carac_pb_elast

    ///Structure contenant les informations sur chaque bord de la sous-structure
    struct Edge
    {
        Edge(){
            mesh=NULL;
        }

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
    unsigned nb_macro;                  ///< nbre de fcts de base macro par sst (somme des fonctions macro des interfaces autour de la sst)
    unsigned nb_macro_espace;           ///< nbre de fcts de base macro spatiale par sst
    Mat <TYPEREEL , Gen<>, Dense<> > LE;   ///< matrice de comportement homogeneise
    TMATS *K;                           ///< matrice de rigidite + directions de recherche (utile pour le solver Cholmod)

    Vec<TYPEREEL> fvol ;   ///< second membre provenant des efforts ou quantites chapeaux volumiques par sous-structure
    LDL_solver l;       ///< solver sparse direct plus rapide que skyline

    /// Structure temporelle contenant differents vecteurs
    struct Time{
        Vec<TYPEREEL> q;   ///< vecteurs solution q
        void allocations(unsigned nbnode,Param &process){
            if (process.rank>0 or process.size==1) q.resize(nbnode);
            if (process.rank>0 or process.size==1) q.set(0.0);
            
        }
    };
    
    Vec<Time> t;        ///< Vecteur contenant pour chaque piquet de temps les solutions
    Vec<Time> t_post;   ///< Vecteurs piquet de temps pour sauvegarder les donnees pour chaque piquet de temps (en incremental, non utile en latin)
    
    Vec<TYPEREEL> Fadd;    ///vecteur effort macro sur le bord d'une SST Fadd

#warning "on suppose qu'en 2D, la connectivite vaut 2, en 3D elle vaut 3  pour les elements de peau"
#if DIM==2
    static const   int nb_nodes_by_element=2;
#else
    static const   int nb_nodes_by_element=3; 
#endif  
    BasicVec<BasicVec<TYPEREEL>,DIM> nodes;                           ///< coordonnées des noeuds de peau d'une sst pour la sortie hdf
    BasicVec<BasicVec<TYPEREEL>,DIM > dep_nodes_skin;                 ///< deplacements aux noeuds de peau d'une sst
    BasicVec<BasicVec<TYPEREEL>,DIM > dep_nodes;                      ///< deplacements aux noeuds de peau d'une sst
    BasicVec<BasicVec<int> > mesh_connectivities_skin;              ///< connectivites du maillage de peau d'une sst pour la sortie hdf (tient compte de la numérotation globale des noeuds)
    BasicVec<BasicVec<int> > mesh_connectivities;                   ///< connectivites du maillage d'une sst pour la sortie hdf (tient compte de la numérotation globale des noeuds)
    BasicVec<BasicVec< TYPEREEL > , DIM*(DIM+1)/2 > sigma_skin;      ///< vecteur contrainte du maillage de peau d'une sst pour la sortie hdf 
    BasicVec<BasicVec< TYPEREEL > , DIM*(DIM+1)/2 > epsilon_skin;    ///< vecteur deformation du maillage de peau d'une sst pour la sortie hdf 
    BasicVec< TYPEREEL > sigma_mises_skin;                             ///< vecteur contrainte de von mises du maillage de peau d'une sst pour la sortie hdf 
    BasicVec<BasicVec< TYPEREEL > , DIM*(DIM+1)/2 > sigma;           ///< vecteur contrainte du maillage de peau d'une sst pour la sortie hdf 
    BasicVec<BasicVec< TYPEREEL > , DIM*(DIM+1)/2 > epsilon;         ///< vecteur deformation du maillage de peau d'une sst pour la sortie hdf 
    BasicVec< TYPEREEL > sigma_mises;      ///< vecteur contrainte de von mises du maillage de peau d'une sst pour la sortie hdf 
    BasicVec< int > num_processor;      ///< vecteur numéro du processeur d'une sst pour la sortie hdf 
    BasicVec< int > num_group;          ///< vecteur numéro des sst pour la sortie hdf 
    BasicVec< int > material;           ///< vecteur numéro des materiaux des sst pour la sortie hdf 
    Vec<int,4> nb_elements_with_type;   ///< utilisé pour connaitre le nombre d'elements d'un type donné pour une SST : Triangle, Quadrilateral, Tetrahedron, Hexahedron
    int nb_nodes_by_element_sst, nb_nodes_by_element_sst_skin, pattern_id;
    Sc2String type_elements_sst, type_elements_sst_skin;
  
    Sst() : pb(*mesh.m,true) {}     ///< constructeur de la formulation pour la sous-structure
    
    ~Sst() {free();}    ///destructeur de la SST
    
    void free(){
        vois.free();
        edge.free();
        LE.free();
        //K est libéré apres factorisation de toute facon et non calculer où y a pas besoin
        fvol.free();
        t.free();
        t_post.free();
    }
    
    /// Trouver une sst à partir de son id----------------------------------------------
    static Sst* find_sst(Vec<Sst> &S,int id_) {
        for (int i_group=0; i_group<S.size(); i_group++) {
            if (S[i_group].id == id_) {
                return &S[i_group];
                break;
            }
        }
    }
    
    /// Trouver l'index d'une sst à partir de son id -----------------------------------
    static int find_index_sst(Vec<Sst> &S, int id_) {
        for (int i_group=0; i_group<S.size(); i_group++) {
            if (S[i_group].id == id_) {
                return i_group;
                break;
            }
        }
    }
    
    void assign_material_on_element(DataUser &data_user){
        //formulation isotrope 
        if (mesh.type_formulation=="isotrope") {
            mesh->elastic_modulus = matprop.elastic_modulus ; 
            mesh->poisson_ratio   = matprop.poisson_ratio   ; 
            mesh->deltaT          = matprop.deltaT          ; 
            mesh->resolution      = matprop.resolution      ; 
            mesh->alpha           = matprop.alpha           ; 
            mesh->f_vol           = matprop.f_vol           ; 
            mesh->density         = matprop.density         ; 
            mesh.load_f_vol_e(matprop.f_vol_e,data_user);
        }
        //formulation orthotrope 
        if (mesh.type_formulation=="orthotrope") {
            mesh->elastic_modulus_1 = matprop.elastic_modulus_1;
            mesh->elastic_modulus_2 = matprop.elastic_modulus_2;
            mesh->elastic_modulus_3 = matprop.elastic_modulus_3;
            mesh->poisson_ratio_12  = matprop.poisson_ratio_12 ;
            mesh->poisson_ratio_13  = matprop.poisson_ratio_13 ;
            mesh->poisson_ratio_23  = matprop.poisson_ratio_23 ;
            mesh->shear_modulus_12  = matprop.shear_modulus_12 ;
            mesh->shear_modulus_13  = matprop.shear_modulus_13 ;
            mesh->shear_modulus_23  = matprop.shear_modulus_23 ;
            mesh->v1                = matprop.v1               ;
            mesh->v2                = matprop.v2               ;
            mesh->deltaT            = matprop.deltaT           ;
            mesh->resolution        = matprop.resolution       ;
            mesh->alpha_1           = matprop.alpha_1          ;
            mesh->alpha_2           = matprop.alpha_2          ;
            mesh->alpha_3           = matprop.alpha_3          ;
            mesh->f_vol             = matprop.f_vol            ;
            mesh->density           = matprop.density          ;
            mesh.load_f_vol_e(matprop.f_vol_e,data_user);
        }
        //formulation mesomodele 
        if (mesh.type_formulation=="mesomodele") {
            mesh->elastic_modulus_1 = matprop.elastic_modulus_1;
            mesh->elastic_modulus_2 = matprop.elastic_modulus_2;
            mesh->elastic_modulus_3 = matprop.elastic_modulus_3;
            mesh->poisson_ratio_12  = matprop.poisson_ratio_12 ;
            mesh->poisson_ratio_13  = matprop.poisson_ratio_13 ;
            mesh->poisson_ratio_23  = matprop.poisson_ratio_23 ;
            mesh->shear_modulus_12  = matprop.shear_modulus_12 ;
            mesh->shear_modulus_13  = matprop.shear_modulus_13 ;
            mesh->shear_modulus_23  = matprop.shear_modulus_23 ;
            mesh->v1                = matprop.v1               ;
            mesh->v2                = matprop.v2               ;
            mesh->deltaT            = matprop.deltaT           ;
            mesh->resolution        = matprop.resolution       ;
            mesh->alpha_1           = matprop.alpha_1          ;
            mesh->alpha_2           = matprop.alpha_2          ;
            mesh->alpha_3           = matprop.alpha_3          ;
            mesh->f_vol             = matprop.f_vol            ;
            mesh->density           = matprop.density          ;
            
            mesh->k_p      = matprop.k_p;
            mesh->m_p      = matprop.m_p;
            mesh->R0       = matprop.R0;
            mesh->couplage = matprop.coefvm_composite;
            
            mesh->Yo           = matprop.Yo;
            mesh->Yc           = matprop.Yc;
            mesh->Ycf          = matprop.Ycf;
            mesh->dmax         = matprop.dmax;
            mesh->b_c          = matprop.b_c;
            mesh->effet_retard = matprop.effet_retard;
            mesh->a            = matprop.a;
            mesh->tau_c        = matprop.tau_c;
            
            mesh.load_f_vol_e(matprop.f_vol_e,data_user);
        }
    }
};

#endif // SST_H

