#ifndef INTERFACE_H
#define INTERFACE_H

#include "main_typedef.h"
#include "SstCarac_InterCarac.h"
#include "../COMPUTE/DataUser.h"
#include "../GEOMETRY/GeometryUser.h"

//***********************************
// classe Interface parametrable
//***********************************
/** \ingroup Interfaces
\brief Classe definissant les interfaces

Classe de stockage des variables associees a une interface
*/
struct Interface
{
    /// Geometrie
    int id;     /// id du group d'interface de geometry_user
    int num;    /// numero de l'interface 
    
    Point G;                             /// centre de gravite
    Vec<Point, DIM> BPI;                 /// Base principale d'inertie de l'interface
    Point Moments_inertie;  /// Moments d'inertie pour savoir si l'interface est plane ou non
    Vec<int> vois;                      /// sous-structures voisines + cote correspondant (ex: sst 0 cote 1 sst 3 cote 2 - > 0 1 3 2)
    Scalar measure;                   /// mesure de l'interface (aire ou longueur)

    /// comportement de l'interface
    Sc2String type; /// type d'interface (exterieure : Ext ou interieure : Int )
    Sc2String comp; /// comportement d'interface (effort : effort, deplacement : depl, symetrie : sym, deplacement normal : depl_normal parfaite : Parfait, contact : Contact, contact epais : Contact_ep
    int refCL;      /// numero de la condition aux limites concernee pour les interfaces exterieures (index dans la liste)
    int edge_id;    /// identite du group edge equivalent dans le json. Ce numero peut être commun a plusieurs interfaces de bord "ext"
    int id_link;    /// identite du behaviour_links dans le data_user issu du json
    int id_bc;      /// identite du behaviour_bc dans le data_user issu du json
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
        TYPEREEL kn,kt,hn,ht,h,k;

        Vec<unsigned> ddlcorresp; /// correspondance des ddls des deux cotes de l'interface (pas encore bien teste)

        /// Structure temporelle contenant les vecteurs nodaux
        struct Time{
            Vector F;        /// efforts sur la face
            Vector W;        /// deplacements sur la face
            Vector Wp;       /// vitesses sur la face
            Vector Fchap;    /// 
            Vector Wchap;    ///
            Vector Wpchap;   ///
            Vector WtildeM;  /// pseudo multiplicateur de
            Vector oldF;     /// efforts a l'iteration (en incremental) ou au pas de temps (en latin) precedent (pour la relaxation)
            Vector oldW;     /// deplacements a l'iteration (en incremental) ou au pas de temps (en latin) precedent (pour la relaxation)
            Vector oldWp;    /// vitesses a l'iteration (en incremental) ou au pas de temps (en latin) precedent (pour la relaxation)
            void allocations(unsigned sizenodeeq);
        };
        Vec<Time> t;        /// Vecteurs piquet de temps
        Vec<Time> t_post;   /// Vecteurs piquet de temps pour sauvegarder les donnees pour chaque piquet de temps (en incremental, non utile en latin)
        
        BasicVec<BasicVec<int> > mesh_connectivities; /// connectivites du maillage de peau d'une sst pour la sortie hdf (tient compte de la numérotation globale des noeuds)
        BasicVec<int> nature; /// type d'interface : 0 : deplacement imposé, 1 : effort imposé, 2 : symetrie, 3 : depl normal imposé, 4 : parfait, 5 : contact
        BasicVec<int> number; /// numéro de l'interface
        BasicVec< BasicVec<TYPEREEL>, DIM > F, Fchap, W, Wchap, Wp, Wpchap; /// champs de déplacement et contrainte
        
    };
    Vec<Side> side; ///< cotes de l'interface

    Vector jeu ;             /// vecteur des valeurs de jeu sur l'interface
    Vector coeffrottement;   /// vecteur des valeurs de coeficient de frottement sur l'interface
    Vector raideur;          /// vecteur des valeurs de raideur sur l'interface
    Vec<bool> comportement;         /// indique s'il y a modification du comportement d'un element
    int convergence;                ///< =-1 si le calcul du pas de temps ne converge pas, >=0 sinon
                                    ///< =0 après l'etape locale si aucun comportement d'element est mis à jour, >0 sinon
    
    struct Time{
        Vector d;    /// endommagement des elements de l'interface
        void allocate(const Interface& interface){
            if(interface.matprop->degradable){
                d.resize(interface.side[0].nodeeq.size());
            }
        }
        void free(){d.free();}
    };
    Vec<Time> t;
                                    
    
    #if DIM==2
    static const   int nb_nodes_by_element=2;
    #else
    static const   int nb_nodes_by_element=3; 
    #endif  
    BasicVec<BasicVec<TYPEREEL>,DIM> nodes; ///< coordonnées des noeuds de peau d'une sst pour la sortie hdf
    BasicVec<BasicVec<int> > mesh_connectivities; ///< connectivites du maillage de peau d'une sst pour la sortie hdf (tient compte de la numérotation globale des noeuds)
    BasicVec<int> nature; ///< type d'interface : 0 : deplacement imposé, 1 : effort imposé, 2 : symetrie, 3 : depl normal imposé, 4 : parfait, 5 : contact
    BasicVec<int> number; ///< numéro de l'interface
    BasicVec< BasicVec<TYPEREEL>, DIM > F, Fchap, W, Wchap, Wp, Wpchap; ///< champs de déplacement et contrainte
    
    //*******************************************************************************************
    // methodes de la class
    //*******************************************************************************************
    void free();
    void affiche();
    void read_data_user(int index,const DataUser &data_user, const GeometryUser &geometry_user);

    Interface();
    ~Interface();
};


#endif // INTERFACE_H
