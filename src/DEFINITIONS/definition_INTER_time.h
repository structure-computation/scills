#ifndef INTER_H
#define INTER_H

#include <string>
using namespace LMT;

//definition du maillage des interfaces
#include "meshcaracinter.h"

//struct PARAM_COMP_INTER;
#include "definition_PARAM_COMP_INTER.h"

//***********************************
// classe Interface parametrable
//***********************************
/** \ingroup Interfaces
\brief Classe Interface parametrable 

dim_ dimension du probleme (accessible tout le temps sans test)
TT_ type des donnees pour les flottants de l'interface
*/
struct Interface
{
    typedef Mat <TYPEREEL , Gen<>, SparseLine<> > TMATS; ///< type de matrice sparseline
    typedef Mat <TYPEREEL , Gen<>, Dense<> > TMATF; ///< type de matrice pleine
    typedef Mesh<Meshcaracinter<DIM,DIM-1> > TMESH; ///< type de maillage des interfaces
    typedef Vec<TYPEREEL,DIM> Pvec; ///< type des points
    typedef PARAM_COMP_INTER PARAM_COMP; ///< type des parametres supplementaires pour les comportements
   
    //Geometrie
    int id;     ///< id du group d'interface de geometry_user
    int num;     ///< numero de l'interface 
    
    Pvec G;  ///< centre de gravite
    Vec<Pvec, DIM> BPI; ///< Base principale d'inertie de l'interface
    Vec<TYPEREEL,DIM> Moments_inertie; ///< Moments d'inertie pour savoir si l'interface est plane ou non
    Vec<int> vois; ///< sous-structures voisines + cote correspondant (ex: sst 0 cote 1 sst 3 cote 2 - > 0 1 3 2)
    Vec<Pvec,2> box; ///< boite contenant l'interface
    TYPEREEL measure; ///< mesure de l'interface (aire ou longueur)

    //comportement de l'interface
    std::string type; ///< type d'interface (exterieure : Ext ou interieure : Int )
    std::string comp; ///< comportement d'interface (effort : effort, deplacement : depl, symetrie : sym, deplacement normal : depl_normal parfaite : Parfait, contact : Contact, vontact epais : Contact_ep
    int refCL;   ///< numero de la condition aux limites concernee pour les interfaces exterieures (index dans la liste)
    int edge_id;   ///< identite du group edge equivalent dans le json. Ce numero peut être commun a plusieurs interfaces de bord "ext"
    int id_link;   ///< identite du behaviour_links dans le data_user issu du json
    int id_bc;     ///< identite du behaviour_bc dans le data_user issu du json
    
    
    // donnees liees au macro
    Vec<int> repddl; ///< reperage des ddls dans le probleme macro
    int nb_macro_total; ///< nombre de fct macro total
    int nb_macro_espace; ///< nombre de fct de base macro en espace



    /// classe decrivant les grandeurs ou champs d'un cote d'une interface
    struct Side
    {
        Side(){
            t.resize(1);
            mesh=NULL;
        }
        // maillage et caracteristiques
        TMESH *mesh; ///< maillage d'interface par cote
        Vec< Pvec > nodeeq; ///< noeuds equivalents sur l'interface
        Vec< TYPEREEL > neq;   ///< normales equivalentes (par element) (nx1,ny1,nx2,ny2...)
        Vec<int,2> vois;      ///< sous-structure voisine et cote correspondant

        /// operateurs
        /// M masse, N operateur de sous_integration 
        /// eM matrice permettant de renvoyer les valeurs pour chaque noeuds à partir des quantites macro, MeM la meme avec multiplication par la matrice de masse :permet d'extraire les composantes macro d'une distribution donnee,
        /// kloc direction de recherche, cloc direction de recherche inverse (sous forme de matrice globale sur l'interface: pas utiliser pour l'instant)
        TMATS M,N,Nt,hglo,kglo;
        TMATF eM,MeM;
        /// fonction projecteur macro . La définition d'un projecteur macro et micro sous forme d'une fonction permet de ne pas construire et stocker de nouvel opérateur et prend autant de temps que l'utilisation d'une matrice.
        Vec<TYPEREEL,-1,void> PM(Vec<TYPEREEL> &f) {
            Vec<TYPEREEL,-1,void> fM; fM.resize(f.size()); fM.set(0.);
            for(unsigned i=0;i<eM.nb_cols();i++) fM+=dot(eM.col(i),f)*eM.col(i);
            return fM;
        }
        /// fonction projecteur micro
        Vec<TYPEREEL,-1,void> Pm(Vec<TYPEREEL> &f) {
            Vec<TYPEREEL,-1,void> fm; fm.resize(f.size()); fm.set(0.);
            for(unsigned i=0;i<eM.nb_cols();i++) fm+=dot(eM.col(i),f)*eM.col(i);
            return f-fm;
        }
        /// fonction projecteur normal
        Vec<TYPEREEL,-1,void> Pn(Vec<TYPEREEL> &f) {
            Vec<TYPEREEL,-1,void> fn; fn.resize(f.size()); fn.set(0.);
            for(unsigned i=0;i<nodeeq.size();++i) { Vec<TYPEREEL,DIM> rep=range(i*DIM,(i+1)*DIM); fn[rep]=dot(neq[rep],f[rep])*neq[rep]; }
            return fn;
        }
        /// fonction projecteur tangentiel
        Vec<TYPEREEL,-1,void> Pt(Vec<TYPEREEL> &f) {
            Vec<TYPEREEL,-1,void> ft; ft.resize(f.size()); ft.set(0.);
            for(unsigned i=0;i<nodeeq.size();++i) { Vec<TYPEREEL,DIM> rep=range(i*DIM,(i+1)*DIM); ft[rep]=dot(neq[rep],f[rep])*neq[rep]; }
            return f-ft;
        }

        /// scalaire pour direction de recherche normale tangentielle, ou direction unique
        TYPEREEL kn,kt,hn,ht,h,k;

        Vec<unsigned> ddlcorresp; ///< correspondance des ddls des deux cotes de l'interface (pas encore bien teste)

        ///< Structure temporelle contenant les vecteurs nodaux
        struct Time{
            Vec<TYPEREEL> F,Wp,W,Fchap,Wpchap,Wchap,WtildeM,oldF,oldWp,oldW;
            void allocations(unsigned sizenodeeq,Param &process){
                if (process.rank>0 or process.size==1) F.resize(sizenodeeq);
                if (process.rank>0 or process.size==1) F.set(0.0);
                if (process.rank>0 or process.size==1) Wp.resize(sizenodeeq);
                if (process.rank>0 or process.size==1) Wp.set(0.0);
                if (process.rank>0 or process.size==1) W.resize(sizenodeeq);
                if (process.rank>0 or process.size==1) W.set(0.0);
                if (process.rank>0 or process.size==1) oldW.resize(sizenodeeq);
                if (process.rank>0 or process.size==1) oldW.set(0.0);
                if (process.rank>0 or process.size==1) Fchap.resize(sizenodeeq);
                if (process.rank>0 or process.size==1) Fchap.set(0.0);
                if (process.rank>0 or process.size==1) Wchap.resize(sizenodeeq);
                if (process.rank>0 or process.size==1) Wchap.set(0.0);
                if (process.rank>0 or process.size==1) Wpchap.resize(sizenodeeq);
                if (process.rank>0 or process.size==1) Wpchap.set(0.0);
                if (process.rank>0 or process.size==1) WtildeM.resize(sizenodeeq);
                if (process.rank>0 or process.size==1) WtildeM.set(0.0);
                if (process.rank>0 or process.size==1) oldF.resize(sizenodeeq);
                if (process.rank>0 or process.size==1) oldF.set(0.0);
                if (process.rank>0 or process.size==1) oldWp.resize(sizenodeeq);
                if (process.rank>0 or process.size==1) oldWp.set(0.0);
            }
        };
        Vec<Time> t; ///< Vecteurs piquet de temps
        Vec<Time> t_post; ///< Vecteurs piquet de temps pour sauvegarder les donnees pour chaque piquet de temps (en incremental, non utile en latin)
        
        BasicVec<BasicVec<int> > mesh_connectivities; ///< connectivites du maillage de peau d'une sst pour la sortie hdf (tient compte de la numérotation globale des noeuds)
        BasicVec<int> nature; ///< type d'interface : 0 : deplacement imposé, 1 : effort imposé, 2 : symetrie, 3 : depl normal imposé, 4 : parfait, 5 : contact
        BasicVec<int> number; ///< numéro de l'interface
        BasicVec< BasicVec<TYPEREEL>, DIM > F, Fchap, W, Wchap, Wp, Wpchap; ///< champs de déplacement et contrainte
        
    };
    Vec<Side> side; ///< cotes de l'interface

    PARAM_COMP_INTER *param_comp; ///<  supplementaires pour le comportement des interfaces (jeu coefficient de frottement Gcrit...)
   
   
    //*******************************************************************************************
    // methodes de la class
    //*******************************************************************************************
   
    void free(){///suppression des interfaces
        vois.free();
        repddl.free();

    #ifdef PRINT_ALLOC
        if (side[0].mesh != NULL) total_allocated[ typeid(typename Interface::TMESH).name() ] -= sizeof(typename Interface::TMESH);
    #endif

        //les maillages sont les memes des 2 cotés...
        if (side.size() != 0){
                if (side[0].mesh != NULL) delete side[0].mesh;
            side[0].mesh=NULL;
                side.free();
                if (param_comp != NULL ) param_comp->free();
            if (param_comp != NULL ) delete param_comp;
        }
    }
    
    #include <boost/concept_check.hpp>
    void affiche(){
        std::cout << "------------------- interface ------------------------" << std::endl;
        PRINT(id);
        std::cout << type << std::endl;
        std::cout << comp << std::endl;
        
        PRINT(refCL);
        PRINT(edge_id);
        PRINT(id_bc);
        PRINT(id_link);
        
        std::cout << "------------------- fin interface ------------------------" << std::endl;
    }

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

    Interface() {param_comp = NULL;}
    ~Interface() {Interface::free();}

};


#endif // INTER_H
