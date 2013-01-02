#ifndef SST_H
#define SST_H

#include "main_typedef.h"
#include "SstCarac_InterCarac.h"
/// definition du solveur LDL
#include "../../LMT/include/containers/matcholamd.h"
#include "../../LMT/include/util/solveLDL.h"

#include "../COMPUTE/DataUser.h"
#include "../GEOMETRY/GeometryUser.h"

//**********************************
// Sous structure
//**********************************
/** \ingroup Sous_structures
\brief Classe Sous-structure parametrable
*/
struct Sst
{
    typedef LMT::Mat<Scalar, LMT::Sym<>, LMT::SparseCholMod > CholModMatrix;    /// type de matrice sparse (pour le solveur CholMod)
    
    /// Donnees geometriques
    Point G;            /// centre de gravite
    Scalar measure;     /// mesure de la sst
    int num;            /// numero de la sous-structure
    int id;             /// id de la sous-structure
    int num_proc;       /// numero du processeur qui traite la sous-structure
    LMT::Vec<int> vois; /// sous-structures voisines (-1 pour le bord)
    
    
    /// Comportement de la SST
    unsigned id_material;   /// id du materiau dans le fichier json
    unsigned typmat;        /// index du SstCarac correspondant au materiau assigne
    bool plastique;         /// indique si le materiau doit etre plastique
    bool endommageable;     /// indique si le materiau doit etre endommageable
    bool update_operator;   /// indique si les operateurs doivent etre reactualises
    SstCarac *matprop;      /// pointeur vers les caracteristiques du materiau
    Problem pb;             /// problème étudié (permet de définir toutes les formulations disponibles)
    Formulation *f;         /// pointeur vers une formulation
    
    SstMesh mesh;           /// maillage de la sst : contient les champs indiques dans Mesh_carac_pb_elast

    
    /// Structure contenant les informations sur chaque bord de la sous-structure
    struct Edge
    {
        Edge(): mesh(NULL) {}
        void affiche() const{
            std::cout << "**************************** Edge ***************************" << std::endl;
            std::cout << "datanum    : " << datanum << std::endl;
            std::cout << "internum   : " << internum << std::endl;
            std::cout << "repddledge : " << repddledge.size() << std::endl;
            std::cout << "repLE      : " << repLE.size() << std::endl;
            std::cout << "G          : " << G << std::endl;
//             std::cout << "nb ddl     : " << mesh->list_group_node.size() << std::endl;
//             std::cout << "nb elem    : " << mesh->list_group_elem.size() << std::endl;
            std::cout << "*************************************************************" << std::endl;
        }

        int datanum;                    /// cote correspondant de l'interface (0 ou 1)
        int internum;                   /// numero de l'interface
        LMT::Vec<unsigned> repddledge;  /// reperage des ddls de bords dans les ddls du maillage
        LMT::Vec<unsigned> repLE;       /// reperage des ddls macro dans l'operateur homogeneise     
        EdgeMesh *mesh;                 /// maillage de bord
        Point G;                        /// centre de gravite du cote
    };
    LMT::Vec<Edge> edge;     /// Bords de la sous-structure
    
    
    /// Operateurs
    unsigned nb_macro;          /// nbre de fcts de base macro par sst (somme des fonctions macro des interfaces autour de la sst)
    unsigned nb_macro_espace;   /// nbre de fcts de base macro spatiale par sst
    DenseMatrix LE;             /// matrice de comportement homogeneise
    CholModMatrix *K;           /// matrice de rigidite + directions de recherche (utile pour le solver Cholmod)
    LDL_solver l;               /// solver sparse direct plus rapide que skyline

    
    /// Structure temporelle contenant les differents vecteurs de stockage des resultats
    struct Time{
        Vector q;                               /// vecteur des deplacements aux noeuds
        Vector p;                               /// vecteur de la plasticite cumulee aux points de Gauss
        Vector R_p;                             /// vecteur de l'ecrouissage aux points de Gauss
        LMT::Vec<VoigtVector> epsilon_p;        /// vecteur des deformations plastiques aux points de Gauss
        LMT::Vec<VoigtVector> X_p;              /// vecteur des centres du domaine elastique de chaque elements
        //Vec<TYPEREEL,2> erreur_plasticite;      /// vecteur des erreurs en plasticite cumulee calculees aux points de Gauss PRISE EN COMPTE ERREUR_PLASTICITE NON IMPLEMENTEE
        Vector d1;                              /// vecteur des endommagements (modele standart ou micro-fissuration du mesomodele) aux points de Gauss
        Vector d2;                              /// vecteur des endommagements (decohesion fibres/matrice du mesomodele) aux points de Gauss
        Vector df;                              /// vecteur des endommagements (rupture fragile des fibres du mesomodele) aux points de Gauss
        LMT::Vec<Vec<Scalar,3> > Yd;            /// vecteur des forces thermodynamyques d'endommagement (type mesomodele) aux points de Gauss
        //Vec<TYPEREEL,2> erreur_endommagement;   /// vecteur des erreurs en endommagement calculees aux points de Gauss PRISE EN COMPTE ERREUR_ENDOMMAGEMENT NON IMPLEMENTEE
        
        Time();
        ~Time();
        void allocate(Sst &S);   /// Alloue les vecteurs de Time
        void free();
        void affiche() const;
    };
    Vec<Time> t;    /// Vecteur contenant pour chaque piquet de temps les solutions
    Vec<Time> t_post;
    
    Vector fvol ;   /// second membre provenant des efforts ou quantites chapeaux volumiques par sous-structure
    Vector Fadd;    /// vecteur effort macro sur le bord d'une SST Fadd

//#warning "on suppose qu'en 2D, la connectivite vaut 2, en 3D elle vaut 3  pour les elements de peau"
    static const int nb_nodes_by_element = DIM; 

    BasicVec<BasicVec<Scalar>,DIM> nodes;                         /// coordonnées des noeuds de peau d'une sst pour la sortie hdf
    BasicVec<BasicVec<Scalar>,DIM > dep_nodes_skin;               /// deplacements aux noeuds de peau d'une sst
    BasicVec<BasicVec<Scalar>,DIM > dep_nodes;                    /// deplacements aux noeuds de peau d'une sst
    BasicVec<BasicVec<int> > mesh_connectivities_skin;            /// connectivites du maillage de peau d'une sst pour la sortie hdf (tient compte de la numérotation globale des noeuds)
    BasicVec<BasicVec<int> > mesh_connectivities;                 /// connectivites du maillage d'une sst pour la sortie hdf (tient compte de la numérotation globale des noeuds)
    BasicVec<BasicVec< Scalar > , DIM*(DIM+1)/2 > sigma_skin;     /// vecteur contrainte du maillage de peau d'une sst pour la sortie hdf 
    BasicVec<BasicVec< Scalar > , DIM*(DIM+1)/2 > epsilon_skin;   /// vecteur deformation du maillage de peau d'une sst pour la sortie hdf 
    BasicVec< Scalar > sigma_mises_skin;                          /// vecteur contrainte de von mises du maillage de peau d'une sst pour la sortie hdf 
    BasicVec<BasicVec< Scalar > , DIM*(DIM+1)/2 > sigma;          /// vecteur contrainte du maillage de peau d'une sst pour la sortie hdf 
    BasicVec<BasicVec< Scalar > , DIM*(DIM+1)/2 > epsilon;        /// vecteur deformation du maillage de peau d'une sst pour la sortie hdf 
    BasicVec< Scalar > sigma_mises;     /// vecteur contrainte de von mises du maillage de peau d'une sst pour la sortie hdf 
    BasicVec< int > num_processor;      /// vecteur numéro du processeur d'une sst pour la sortie hdf 
    BasicVec< int > num_group;          /// vecteur numéro des sst pour la sortie hdf 
    BasicVec< int > material;           /// vecteur numéro des materiaux des sst pour la sortie hdf 
    Vec<int,4> nb_elements_with_type;   /// utilisé pour connaitre le nombre d'elements d'un type donné pour une SST : Triangle, Quadrilateral, Tetrahedron, Hexahedron
    int nb_nodes_by_element_sst, nb_nodes_by_element_sst_skin, pattern_id;
    Sc2String type_elements_sst, type_elements_sst_skin;
  
    Sst();     /// constructeur de la formulation pour la sous-structure
    ~Sst();    /// destructeur de la Sst
    
    //void resize_results(unsigned nb_t,bool allocate_t,bool recopie_t_post);   /// Alloue les vecteurs de resultats de la Sst
    void free();                                                        /// Libere la Sst
    
    /// Lis les donnes du DataUser (A REVOIR)
    void read_data_user(int index,const DataUser &data_user,GeometryUser &geometry_user);
    
    /// Trouver une sst à partir de son id----------------------------------------------
    static Sst* find_sst(LMT::Vec<Sst> &S,int id_);
    
    /// Trouver l'index d'une sst à partir de son id -----------------------------------
    static int find_index_sst(LMT::Vec<Sst> &S, int id_);
    
    /// Assigne les proprietes materiau au maillage associe a la Sst
    void apply_behavior();
    
    /// Intersecte le maillage de la sst et de ses bords, etablis les tables de correspondance (avec verification si elles sont completes)
    void calc_SST_Correspddl(); // A REMETTRE A CA PLACE
    
    /// Fonction d'affichage pour le debug
    void affiche() const;
};


#endif // SST_H

