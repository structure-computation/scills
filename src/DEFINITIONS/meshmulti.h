//
// C++ Interface: meshmulti
//
// Description: 
//
//
// Author: Alain CAIGNOT <caignot@lmt.ens-cachan.fr>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MESHMULTI_H
#define MESHMULTI_H
#include "create_new_elem_p.h"
#include <ext/hash_map>

#include "mesh/read_avs.h"
#include "mesh/read_geof.h"

#include "GeometryUser.h"
#include "DataUser.h"
#include "codegen/codegen.h"
#include <boost/concept_check.hpp>

using namespace std;
using namespace LMT;
using namespace Codegen;
using namespace __gnu_cxx;

///application de flag indiquant le type de materiau et le numero de la sst pour chaque element (utile pour l'affichage)
struct apply_mat_elem {
  template<class TE>
      void operator() ( TE &e, const int &typmat, const int &numsst,const int &num_proc) const {
    e.typmat=typmat;
    e.numsst=numsst;
    e.num_proc=num_proc;
      }
};

/** \ingroup  maillages 
\brief Structure pour le hashage du maillage cf. p_surdiscretise_SST
*/
struct MyHash {
    template<class TV>
    unsigned operator()(const TV &vec1) const {
        return vec1.num_hash;
    }
};

/** \ingroup  maillages 
\brief Structure pour reperer les noeuds dans la hash_map
*/
struct NodesEq {
    template<class TV>
    bool operator()(const TV &n1, const TV &n2) const {
        double eps=1e-6;
        return (length(n1.pos-n2.pos)<=eps);
    }
};
/** \ingroup  maillages 
\brief Structure pour le hashage du maillage
*/
template<class T, unsigned dim>
struct Noeud_Hash {
    Noeud_Hash() {
        num_hash=0;
        pos.set(0);
    }
    Vec<T, dim> pos;
    unsigned num_hash;
};

/** \ingroup  maillages 
\brief Ajout de noeuds sur les aretes d'un maillage
*/
struct add_nodes {
    template<class TE,class TM, class TH>
    void operator() (TE &e, TM &m, TH &hm) const {
        typedef typename TM::Pvec Pvec;
        //ajout du noeud au maillage
        Pvec G = center(e);
        m.add_node(G);
        //creation de la structure pour hashage
        Noeud_Hash<typename TE::T,TM::dim> newnoeud;
        newnoeud.pos=G;

        //reperage des numeros des noeuds extremes de l'element
        unsigned num_hash=0;
        for(unsigned i=0;i<e.nb_nodes;i++)
            num_hash+=e.node(i)->number_in_original_mesh();
        //assignation de la valeur de hashage (somme des numeros des noeuds extremes)
        newnoeud.num_hash=num_hash;
        //ajout à la table de hashage de l'entree newnoeud correspondant au numero dans le maillage du dernier noeud entre.
        hm[newnoeud]=m.node_list.size()-1;
    }
};

/** \ingroup  Sous_structures 
\brief Modification de l'element
*/
template<unsigned number>
struct ModifTypeElem {
    template<class TE, class TM, class TNH>
    void operator()(TE &e,TM &m, hash_map<TNH, unsigned, MyHash, NodesEq> &hm, TM &m2) const {
        unsigned nb_children = NbChildrenElement<typename TE::NE,number>::res;

        //boucle sur les elements enfants
        Vec<unsigned> rep_nodes;
        rep_nodes.resize(nb_children+e.nb_nodes,0);
        for(unsigned i=0;i<e.nb_nodes;i++)
            rep_nodes[i]=e.node(i)->number_in_original_mesh();

        for(unsigned i=0;i<nb_children;i++) {
            typename TM::EA *ea = m.get_children_of(e,LMT::Number<number>())[i]; //element ancestor i;
            typename TM::Pvec G = center(*ea);

            Noeud_Hash<typename TE::T,TM::dim> newnoeud;
            newnoeud.pos=G;
            unsigned num_hash=0;
            for(unsigned j=0;j<ea->nb_nodes_virtual();j++)
                num_hash+=ea->node_virtual(j)->number_in_original_mesh();

            //assignation de la valeur de hashage (somme des numeros des noeuds extremes)
            newnoeud.num_hash=num_hash;
            rep_nodes[i+e.nb_nodes]=hm[newnoeud];
        }

        create_new_elem_p(e,m2,rep_nodes);
    }
};

/** \ingroup  Sous_structures 
\brief lecture du maillage à paritr de geometry_user issue de SC_create_2
*/
template<class TM>
void read_mesh_sst_geometry_user(TM &mesh, GeometryUser &geometry_user, int id_sst) throw(std::runtime_error) {
    //TM mesh;
    typedef typename TM::Tpos T;
    typedef typename TM::Pvec Pvec;
    typedef typename TM::TNode TNode;
    typedef typename TM::EA EA;
  
    // obtaining nbnode, nbelem
    unsigned nbnode = geometry_user.find_group_elements(id_sst)->map_mesh_nodes.size();
    unsigned nbelem = geometry_user.find_group_elements(id_sst)->nb_elements;
  
    
    //ajout des noeuds au maillage
    map<int,TNode *> map_num_node;
    Vec<TYPE,DIM> vec;
    for(int i_node=0; i_node<nbnode; i_node++){
        for(unsigned d=0; d<DIM; d++){
            vec[d] = geometry_user.find_group_elements(id_sst)->local_nodes[d][i_node];
        }
        map_num_node[i_node] = mesh.add_node(vec);
    }
    
    //ajout des elements
    switch (geometry_user.find_group_elements(id_sst)->pattern_base_id){
        //for Triangle
        case 0 :{
//             PRINT("Triangle");
            int nb_node_elem = 3;
            Vec<TNode *> vn;
            vn.resize(nb_node_elem);
            for(int i_elem=0; i_elem<nbelem; i_elem++) {
                for(int i_node=0; i_node<nb_node_elem; i_node++) {
                    vn[i_node] = map_num_node[geometry_user.find_group_elements(id_sst)->local_connectivities[i_node][i_elem]];
                }
                typename TM::EA *ne = reinterpret_cast<typename TM::EA *>(mesh.add_element(Triangle(),DefaultBehavior(),&vn[0]));
                ne->group = id_sst;
            }
            break;
        }
        //for Triangle_6
        case 1 :{
//             PRINT("Triangle_6");
            int nb_node_elem = 6;
            Vec<TNode *> vn;
            vn.resize(nb_node_elem);
            for(int i_elem=0; i_elem<nbelem; i_elem++) {
                for(int i_node=0; i_node<nb_node_elem; i_node++) {
                    vn[i_node] = map_num_node[geometry_user.find_group_elements(id_sst)->local_connectivities[i_node][i_elem]];
                }
                typename TM::EA *ne = reinterpret_cast<typename TM::EA *>(mesh.add_element(Triangle_6(),DefaultBehavior(),&vn[0]));
                ne->group = id_sst;
            }
            break;
        }
        //for Tetra
        case 2 :{
//             PRINT("Tetra");
            int nb_node_elem = 4;
            Vec<TNode *> vn;
            vn.resize(nb_node_elem);
            for(int i_elem=0; i_elem<nbelem; i_elem++) {
                for(int i_node=0; i_node<nb_node_elem; i_node++) {
                    vn[i_node] = map_num_node[geometry_user.find_group_elements(id_sst)->local_connectivities[i_node][i_elem]];
                }
                typename TM::EA *ne = reinterpret_cast<typename TM::EA *>(mesh.add_element(Tetra(),DefaultBehavior(),&vn[0]));
                ne->group = id_sst;
            }
            break;
        }
        //for Tetra_10
        case 3 :{
//             PRINT("Tetra_10");
            int nb_node_elem = 10;
            Vec<TNode *> vn;
            vn.resize(nb_node_elem);
            for(int i_elem=0; i_elem<nbelem; i_elem++) {
                for(int i_node=0; i_node<nb_node_elem; i_node++) {
                    vn[i_node] = map_num_node[geometry_user.find_group_elements(id_sst)->local_connectivities[i_node][i_elem]];
                }
                typename TM::EA *ne = reinterpret_cast<typename TM::EA *>(mesh.add_element(Tetra_10(),DefaultBehavior(),&vn[0]));
                ne->group = id_sst;
            }
            break;
        }
        default :{
            std::cerr << "type de pattern non implemente" << std::endl; assert(0);                    
        }
        mesh.sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( mesh, mesh.elem_list.size() / 2 + 1 );
    }
}

/** \ingroup  maillages 
\brief assignation de la force volumique par elements
*/
struct assigne_f_vol_e {
    template<class TE,class TM>
    void operator() (TE &e, TM &m, Vec<std::string> force_volumique) const {
        typedef typename TM::Pvec Pvec;
        //ajout du noeud au maillage
	
        Pvec G = center(e);
	
        
        std::vector<Ex> symbols;
        if (DIM==2) {
            symbols.push_back("x");
            symbols.push_back("y");
        }
        else if (DIM==3) {
            symbols.push_back("x");
            symbols.push_back("y");
            symbols.push_back("z");
        }
        
        Vec<Ex> expr;
        expr.resize(DIM);
        for(unsigned d2=0;d2<DIM;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
            expr[d2] = read_ex(force_volumique[d2],symbols);
        }
        
        Vec<double,DIM> data;
        Ex::MapExNum var;
        for(unsigned d2=0;d2<DIM;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
            var[symbols[d2]]= G[d2];
        }
        
        for(unsigned d2=0;d2<DIM;++d2){//boucle sur les inconnues possibles (dimension des vecteurs)
            //e.f_vol_e[d2] = measure(e) * m.density * (double)expr[d2].subs_numerical(var);
	    e.f_vol_e[d2] = m.density * (double)expr[d2].subs_numerical(var);
        }
        
//         if(e.number<10){
// 	    PRINT("-----------element------------" );
// 	    PRINT(e.number);
// 	    for(unsigned d2=0;d2<DIM;++d2) {
// 		std::cout << "        " << G[d2] << std::endl;
// 	    }
// 	    for(unsigned d2=0;d2<DIM;++d2) {
// 		std::cout << " 		force_volumique = " << force_volumique[d2] << std::endl;
// 		std::cout << " 		e.f_vol_e[d2] = " << e.f_vol_e[d2] << std::endl;
// 	    }
// 	}
    }
};


/** \ingroup  Sous_structures 
\relates Mesh
\brief Classe de maillage des sous-structures

Cette classe permet d'avoir quelques caractéristiques utiles pour le calcul LATIN.
Elle possède un maillage m de type pointeur, on définit l'opérateur -> qui renvoie le maillage directement sans avoir a faire mesh.m-> ...

On lui attache les méthodes suivantes :
- load() 
- unload()
- sousint()
- sousint_p()

*/

template<class Carac>
struct Meshmulti {
    typedef LMT::Mesh<Carac> TM;//type de maillage utilsé pour les sous-structures
    TM *m;//pointeur du maillage
    string name;//nom du maillage a lire
    string sousintegration;//type de sousintegration
    bool flag;//flag permettant de savoir si le maillage est chargé ou non
    int typmat,numsst,num_proc;// caractéristiques associés a la sous-structure que l'on remet automatiquement à la relecteur d'un maillages
    unsigned node_list_size,elem_list_size;// caractéristiques utilisés sans besoin de travailler sur le maillage, on peut ne pas charger le maillage
    Vec<double,Carac::dim> f_vol;//champs de force volumique
    Vec<string,Carac::dim> f_vol_e;//champs de force volumique par element
    double elastic_modulus,poisson_ratio,density,deltaT,resolution,alpha,elastic_modulus_1,elastic_modulus_2,elastic_modulus_3,poisson_ratio_12,poisson_ratio_13,poisson_ratio_23,shear_modulus_12,shear_modulus_13,shear_modulus_23,alpha_1,alpha_2,alpha_3,viscosite;
    Vec<double,Carac::dim> v1,v2;
    string type_formulation;
    
    //ajout pour les données venant de SC_create_2
    int id_sst;
    GeometryUser *geometry_user;  //pointeur vers geometry_user
    //
    Meshmulti() {
        flag=0;
        typmat=0;
        numsst=0;
        num_proc=0;
    }
    //
    void unload() {///dechargement du maillage
        if (flag !=0) {
            delete m;
#ifdef PRINT_ALLOC
            total_allocated[ typeid(TM).name() ] -= sizeof(TM);
#endif
            flag = 0;
        }
    }
    //
    TM *operator->() { load(); return m; }///operateur -> qui permet de retourner directement le mesh mesh.m->pouet remplacer par mesh->pouet
    //
    void load() {///chargement du maillage AVS seulement pour le moment
        if (name != "" and flag == 0) {
#ifdef PRINT_ALLOC
            total_allocated[ typeid(TM).name() ] += sizeof(TM);
#endif
            m = new TM;
            read_mesh_sst_geometry_user(*m, *geometry_user, id_sst);
            flag=1;
            //affiche();
            if (typmat!=0 or numsst!=0 or num_proc!=0) apply(m->elem_list,apply_mat_elem(),typmat,numsst,num_proc);
            if (sousintegration == "p") sousint();
        }
    }
    void load(GeometryUser &geometry_user_,int id_sst_) {///chargement du maillage à partir de geometry_user
        if (name != "" and flag == 0) {
#ifdef PRINT_ALLOC
            total_allocated[ typeid(TM).name() ] += sizeof(TM);
#endif
            //ajout pour SC_create_2
            id_sst = id_sst_;
            geometry_user = &geometry_user_;
            
            node_list_size = geometry_user->find_group_elements(id_sst)->map_mesh_nodes.size();
            elem_list_size = geometry_user->find_group_elements(id_sst)->nb_elements;
        }
    }
    void load_f_vol_e() {///application du chargement à chaque noeud
        std::cout << "fvole " << f_vol_e << endl;
        apply(m->elem_list,assigne_f_vol_e(),*m,f_vol_e);
    }
    //    
    void sousint() {///sousintegration de type p
        sousintegration = "p";
        m->update_elem_children();
        m->update_elem_parents();
        //1ere etape : ajout des noeuds et creation de la table de hashage
        typedef Noeud_Hash<typename Meshmulti<Carac>::TM::Tpos, Carac::dim> TNH;
        hash_map<TNH, unsigned, MyHash, NodesEq> hm;
        apply(m->sub_mesh(LMT::Number<Carac::dim-1>()).elem_list,add_nodes(),*m,hm);
        //2eme etape : boucle sur les noeuds et creation de nouveaux éléments
#ifdef PRINT_ALLOC
        total_allocated[ typeid(TM).name() ] += sizeof(TM);
#endif
        TM *m2 = new TM;
        for(unsigned i=0;i<m->node_list.size();i++)
            m2->add_node(m->node_list[i].pos);
        apply(m->elem_list,ModifTypeElem<Carac::dim-1>(),*m,hm,*m2);
        unload();
        m = m2;
        flag=1;
        m->sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( *m, m->elem_list.size() /2 +1);
        m->sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( *m, m->elem_list.size() /2 +1);
        if (typmat!=0 or numsst!=0 or num_proc!=0) apply(m->elem_list,apply_mat_elem(),typmat,numsst,num_proc);
        node_list_size=m->node_list.size();
        elem_list_size=m->elem_list.size();
    }
    
    void affiche(){/// affiche les infos du maillage pour vérification
        PRINT("---------------affich maillage-------------------");
        PRINT(typmat);
        PRINT(numsst);
        PRINT(num_proc);
        PRINT(m->elem_list.size());
        PRINT(m->node_list.size());
    }
    
    
};




#endif //MESHMULTI_H

