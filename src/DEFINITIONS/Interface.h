#ifndef INTERFACE_H
#define INTERFACE_H

#include "main_typedef.h"
#include "SstCarac_InterCarac.h"
#include "../COMPUTE/DataUser.h"
#include "../GEOMETRY/GeometryUser.h"
/*
namespace Metil {
template<class n >
struct N;
}
*/
//***********************************
// classe Interface parametrable
//***********************************
/** \ingroup Interfaces
\brief Classe definissant les interfaces

Classe de stockage des variables associees a une interface
*/
struct Interface
{
    static Sc2String type_ext;  /// Nom pour le type interface exterieure (CL)
    static Sc2String type_int;  /// Nom pour le type interface interieure (liaison)
    
    static Sc2String comp_parfait;              /// Nom pour le comportement parfait
    static Sc2String comp_elastique;            /// Nom pour le comportement elastique
    static Sc2String comp_contact_parfait;      /// Nom pour le comportement contact de type parfait
    static Sc2String comp_contact_elastique;    /// Nom pour le comportement contact de type elastique
    static Sc2String comp_cohesive;             /// Nom pour le comportement cohesif
    static Sc2String comp_cassable_parfait;     /// Nom pour le comportement cassable
    static Sc2String comp_cassable_elastique;   /// Nom pour le comportement cassable
    
    static Sc2String comp_deplacement;          /// Nom pour le comportement deplacement impose
    static Sc2String comp_deplacement_normal;   /// Nom pour le comportement deplacement normal
    static Sc2String comp_deplacement_nul;      /// Nom pour le comportement deplacement nul
    static Sc2String comp_vitesse;              /// Nom pour le comportement vitesse imposee
    static Sc2String comp_vitesse_normale;      /// Nom pour le comportement vitesse normale
    static Sc2String comp_vitesse_nulle;        /// Nom pour le comportement vitesse nulle
    static Sc2String comp_effort;               /// Nom pour le comportement effort impose
    static Sc2String comp_effort_normal;        /// Nom pour le comportement effort normal (pression)
    static Sc2String comp_symetrie;             /// Nom pour le comportement symetrie
    static Sc2String comp_periodique;           /// Nom pour le comportement periodique
    
    /// Geometrie
    int id;     /// id du group d'interface de geometry_user
    int num;    /// numero de l'interface 
    
    Point G;                /// centre de gravite
    Vec<Point, DIM> BPI;    /// Base principale d'inertie de l'interface
    Point Moments_inertie;  /// Moments d'inertie pour savoir si l'interface est plane ou non
    Vec<int> vois;          /// sous-structures voisines + cote correspondant (ex: sst 0 cote 1 sst 3 cote 2 - > 0 1 3 2)
    Scalar measure;         /// mesure de l'interface (aire ou longueur)

    /// comportement de l'interface
    Sc2String type;         /// type d'interface (exterieure : Ext ou interieure : Int )
    Sc2String comp;         /// comportement d'interface (effort : effort, deplacement : depl, symetrie : sym, deplacement normal : depl_normal parfaite : Parfait, contact : Contact, contact epais : Contact_ep
    int refCL;              /// numero de la condition aux limites concernee pour les interfaces exterieures (index dans la liste)
    int edge_id;            /// identite du group edge equivalent dans le json. Ce numero peut être commun a plusieurs interfaces de bord "ext"
    int id_link;            /// identite du behaviour_links dans le data_user issu du json
    int id_bc;              /// identite du behaviour_bc dans le data_user issu du json
    InterCarac* matprop;    /// Caracteristiques materiaux de l'interface
    
    
    /// donnees liees au pb macro
    Vec<int> repddl;        /// reperage des ddls dans le probleme macro
    int nb_macro_total;     /// nombre de fct macro total
    int nb_macro_espace;    /// nombre de fct de base macro en espace



    /// classe decrivant les grandeurs ou champs d'un cote d'une interface
    struct Side
    {
        Side();
        /// maillage et caracteristiques
        InterfaceMesh *mesh;            /// maillage d'interface par cote
        LMT::Vec<Point> nodeeq;     /// noeuds equivalents sur l'interface
        Vector neq;    /// normales equivalentes (par element) (nx1,ny1,nx2,ny2...)
        Vec<int,2> vois;        /// sous-structure voisine et cote correspondant

        /// operateurs
        /// M masse, N operateur de sous_integration 
        /// eM matrice permettant de renvoyer les valeurs pour chaque noeuds à partir des quantites macro, MeM la meme avec multiplication par la matrice de masse :permet d'extraire les composantes macro d'une distribution donnee,
        /// kloc direction de recherche, cloc direction de recherche inverse (sous forme de matrice globale sur l'interface: pas utiliser pour l'instant)
        SparseMatrix M;        /// Operateur de masse
        SparseMatrix N;        /// Operateur de sous-integration
        SparseMatrix Nt;       /// 
        SparseMatrix hglo;     
        SparseMatrix kglo;     
        DenseMatrix eM;       /// Projecteur sur la base macroscopique
        DenseMatrix MeM;      /// Matrice eM premultipliee par la matrice de masse
        /// fonction projecteur macro . La définition d'un projecteur macro et micro sous forme d'une fonction permet de ne pas construire et stocker de nouvel opérateur et prend autant de temps que l'utilisation d'une matrice.
        Vector PM(Vector &f);
        /// fonction projecteur micro
        Vector Pm(Vector &f);
        /// fonction projecteur normal
        Vector Pn(Vector &f);
        /// fonction projecteur tangentiel
        Vector Pt(Vector &f);

        /// scalaire pour direction de recherche normale tangentielle, ou direction unique
        Scalar kn,kt,hn,ht,h,k;

        Vec<unsigned> ddlcorresp; /// correspondance des ddls des deux cotes de l'interface (pas encore bien teste)

        /// Structure temporelle contenant les vecteurs nodaux
        struct Time{
            Vector F;       /// efforts sur la face (etape lineaire)
            Vector W;       /// deplacements sur la face (etape lineaire)
            Vector Wp;      /// vitesses sur la face (etape lineaire)
            Vector Fchap;   /// efforts sur la face (etape locale)
            Vector Wchap;   /// deplacements sur la face (etape locale)
            Vector Wpchap;  /// vitesses sur la face (etape locale)
            Vector WtildeM; /// pseudo multiplicateur de Lagrange (c.f. pb macro)
            Vector oldF;    /// efforts a l'iteration (en latin) ou au pas de temps (en incremental) precedent
            Vector oldW;    /// deplacements a l'iteration (en latin) ou au pas de temps (en incremental) precedent
            Vector oldWp;   /// vitesses a l'iteration (en latin) ou au pas de temps (en incremental) precedent
            Vector d;       /// endommagement des elements de l'interface
            
            void allocations(unsigned sizenodeeq,bool endommageable);
        };
        Vec<Time> t;        /// Vecteurs piquet de temps
        Vec<Time> t_post;   /// Vecteurs piquet de temps pour sauvegarder les donnees pour chaque piquet de temps (en incremental, non utile en latin)
        
        BasicVec<BasicVec<int> > mesh_connectivities; /// connectivites du maillage de peau d'une sst pour la sortie hdf (tient compte de la numérotation globale des noeuds)
        BasicVec<int> nature; /// type d'interface, c.f. int get_type_elem()
        BasicVec<int> number; /// numéro de l'interface
        BasicVec< BasicVec<Scalar>, DIM > F, Fchap, W, Wchap, Wp, Wpchap; /// champs de déplacement et contrainte
    };
    Vec<Side> side; ///< cotes de l'interface

    Scalar coeffrottement;      /// coefficient de frottement global
    Vector coeffrottement_vec;  /// vecteur des valeurs du coefficient de frottement
    Vector Ep_imposee;          /// valeurs du jeu ou de l'epaisseur imposee par l'utilisateur sur l'interface
    Vector old_Ep_imposee;      /// valeurs de Ep_impose au pas de temps precedent
    Vector Ep_elastique;        /// jeu ou epaisseur cree par l'elasticite de l'interface
    Vector old_Ep_elastique;    /// valeurs de Ep_elastique au pas de temps precedent
    Vector precharge;           /// précharge imposee par l'utilisateur sur l'interface
    Vector old_precharge;       /// précharge imposee par l'utilisateur sur l'interface
    Vec<unsigned> comportement;     /// indique pour chaque element s'il y a modification du comportement, si on casse ou non
    int convergence;            ///< =-1 si le calcul du pas de temps ne converge pas, >=0 sinon
                                ///< =0 après l'etape locale si aucun comportement d'element est mis à jour, >0 sinon
    
    /// Structure permettant de definir l'operateur de direction de recherche local a partir de direction normale et tangentielle
    /// S'utilise ensuite comme une matrice
    struct LocalOperator{
        /// Rigidites normale (kn) et tangentielle (kt)
        Scalar kn,kt;
        /// Vecteur normal
        Point n;
        /// Surcharge pour la multiplication d'un vecteur W
        template<class TV> Point operator*(TV &W){return kn*n*dot(n,W)+kt*(W-n*dot(n,W));}
        /// Surcharge pour la multiplication d'un vecteur W constant
        template<class TV> Point operator*(const TV &W){return kn*n*dot(n,W)+kt*(W-n*dot(n,W));}
    };
    
    /// Structure pour simplifier l'acces aux donnees en un noeud de l'interface
    struct NodalState {
        Interface &interface;       /// Pointeur sur l'interface associee
        InterCarac *matprop;        /// Pointeur sur les caracteristiques associees
        unsigned id;                /// Numero du noeud
        unsigned i_time;            /// Numero du noeud
        Scalar dt;                  /// Numero du noeud
        
        /// Vecteurs globaux renumerotees (_F1 = side[0].F[list1] par exemple)
        Vec<unsigned> _comportement;/// Cohesion entre les cotes?
        Vector _coeffrottement;     /// Coefficient de frottement
        Vector _n1;                 /// Normale du bord 1 vers le bord 2
        Vector _F1;                 /// Efforts sur le bord 1 (etape lineaire)
        Vector _F2;                 /// Efforts sur le bord 2 (etape lineaire)
        Vector _Wp1;                /// Vitessse de deplacement du bord 1 (etape lineaire)
        Vector _Wp2;                /// Vitessse de deplacement du bord 2 (etape lineaire)
        Vector _old_Wchap1;         /// Ancienne valeur du deplacement du bord 1 (etape locale au pas de temps precedent)
        Vector _old_Wchap2;         /// Ancienne valeur du deplacement du bord 2 (etape locale au pas de temps precedent)
        Vector _Fchap1;             /// Efforts sur le bord 1 (etape locale)
        Vector _Fchap2;             /// Efforts sur le bord 2 (etape locale)
        Vector _Wpchap1;            /// Vitessse de deplacement du bord 1 (etape locale)
        Vector _Wpchap2;            /// Vitessse de deplacement du bord 2 (etape locale)
        Vector _Wchap1;             /// Deplacement du bord 1 (etape locale)
        Vector _Wchap2;             /// Deplacement du bord 2 (etape locale)
        Vector _Precharge;          /// Precharge dans l'interface
        Vector _old_Precharge;      /// Ancienne valeur de la precharge
        Vector _Ep_imposee;         /// Epaisseur imposee par l'utilisateur
        Vector _old_Ep_imposee;     /// Ancienne valeur de l'epaisseur imposee par l'utilisateur
        Vector _Ep_elastique;       /// Modification de l'epaisseur de l'interface due a l'elasticite
        Vector _old_Ep_elastique;   /// Ancienne valeur de la modification de l'epaisseur due a l'elasticite
        Vector _d;                  /// Endommagement
        Vector _old_d;              /// Ancienne valeur de l'endommagement (pas de temps precedent)
        
        /// Valeurs locals des grandeurs de l'interface
        unsigned comportement;      /// Cohesion entre les cotes?
        Scalar coeffrottement;      /// Coefficient de frottement
        Point n1;                   /// Normale du bord 1 vers le bord 2
        Point F1;                   /// Efforts sur le bord 1 (etape lineaire)
        Point F2;                   /// Efforts sur le bord 2 (etape lineaire)
        Point Wp1;                  /// Vitessse de deplacement du bord 1 (etape lineaire)
        Point Wp2;                  /// Vitessse de deplacement du bord 2 (etape lineaire)
        Point Fchap1;               /// Efforts sur le bord 1 (etape locale)
        Point Fchap2;               /// Efforts sur le bord 2 (etape locale)
        Point Wpchap1;              /// Vitessse de deplacement du bord 1 (etape locale)
        Point Wpchap2;              /// Vitessse de deplacement du bord 2 (etape locale)
        Point Wchap1;               /// Deplacement du bord 1 (etape locale)
        Point Wchap2;               /// Deplacement du bord 2 (etape locale)
        Point Ep_imposee;           /// Epaisseur imposee par l'utilisateur
        Point old_Ep_imposee;       /// Ancienne valeur de l'epaisseur imposee par l'utilisateur
        Point Ep_elastique;         /// Modification de l'epaisseur de l'interface due a l'elasticite
        Point old_Ep_elastique;     /// Ancienne valeur de la modification de l'epaisseur due a l'elasticite
        Point Precharge;            /// Precharge dans l'interface
        Point old_Precharge;        /// Ancienne valeur de la precharge
        Scalar d;                   /// Endommagement
        Point old_Wchap1;           /// Ancienne valeur du deplacement du bord 1 (etape locale au pas de temps precedent)
        Point old_Wchap2;           /// Ancienne valeur du deplacement du bord 2 (etape locale au pas de temps precedent)
        Scalar old_d;               /// Ancienne valeur de l'endommagement (pas de temps precedent)
        
        LocalOperator k1;           /// Direction (en raideur) de recherche sur le bord 1
        LocalOperator k2;           /// Direction (en raideur) de recherche sur le bord 2
        LocalOperator h1;           /// Direction (en souplesse) de recherche sur le bord 1
        LocalOperator h2;           /// Direction (en souplesse) de recherche sur le bord 2
        LocalOperator K;            /// Operateur d'elasticite
        
        NodalState(Interface &I,unsigned pt,Scalar dt_);    /// Construit un NodalState associee a l'interface I, au pas de temps pt (dt_ permet de calculer les derivees)
        void set_node(unsigned i_node);                     /// Positionne le NodalState au noeud indice i_node
        void store_results();                               /// Conserve les resultats dans un stockage intermediaire
        void save_results();                                /// Sauvegarde les resultats sur l'interface
        
        void comportement_parfait();
        void comportement_elastique();
        void comportement_cassable();
        void comportement_cohesif();
        void comportement_contact_parfait();
        void comportement_contact_elastique();
        
        void check_ddr();
        void check_comportement_parfait();
        void check_comportement_elastique();
        void check_comportement_cassable();
        void check_comportement_cohesif();
        void check_comportement_contact_parfait();
        void check_comportement_contact_elastique();
    };
    
    #if DIM==2
    static const int nb_nodes_by_element=2;
    #else
    static const int nb_nodes_by_element=3;
    #endif  
    BasicVec<BasicVec<Scalar>,DIM> nodes; ///< coordonnées des noeuds de peau d'une sst pour la sortie hdf
    BasicVec<BasicVec<int> > mesh_connectivities; ///< connectivites du maillage de peau d'une sst pour la sortie hdf (tient compte de la numérotation globale des noeuds)
    BasicVec<int> nature; ///< type d'interface : 0 : deplacement imposé, 1 : effort imposé, 2 : symetrie, 3 : depl normal imposé, 4 : parfait, 5 : contact
    BasicVec<int> number; ///< numéro de l'interface
    BasicVec< BasicVec<Scalar>, DIM > F, Fchap, W, Wchap, Wp, Wpchap; ///< champs de déplacement et contrainte
    
    //*******************************************************************************************
    // methodes de la class
    //*******************************************************************************************
    void allocate(unsigned nbpastemps);
    void free();
    void init(unsigned pt);
    void affiche();
    void read_data_user(int index,const DataUser &data_user, const GeometryUser &geometry_user);
    
    /// Retourne le type d'element de l'interface sous forme d'un entier
    int get_type_elem() const;

    Interface();
    ~Interface();
};


#endif // INTERFACE_H

class stat;
