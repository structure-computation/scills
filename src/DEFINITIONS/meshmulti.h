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

using namespace std;
using namespace LMT;
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
            typename TM::EA *ea = m.get_children_of(e,Number<number>())[i]; //element ancestor i;
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
    double elastic_modulus,poisson_ratio,deltaT,resolution,alpha,elastic_modulus_1,elastic_modulus_2,elastic_modulus_3,poisson_ratio_12,poisson_ratio_13,poisson_ratio_23,shear_modulus_12,shear_modulus_13,shear_modulus_23,v1,v2,alpha_1,alpha_2,alpha_3,viscosite;
    string type_formulation;
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
            if ((const int &)name.find(".avs")!=-1) read_avs(*m,name.c_str());
            if ((const int &)name.find(".geof")!=-1) read_geof(*m,name.c_str());
            flag=1;
            if (typmat!=0 or numsst!=0 or num_proc!=0) apply(m->elem_list,apply_mat_elem(),typmat,numsst,num_proc);
            node_list_size=m->node_list.size();
            elem_list_size=m->elem_list.size();
            if (sousintegration == "p") sousint();
        }
    }
    //    
    void sousint() {///sousintegration de type p
        sousintegration = "p";
        m->update_elem_children();
        m->update_elem_parents();
        //1ere etape : ajout des noeuds et creation de la table de hashage
        typedef Noeud_Hash<typename Meshmulti<Carac>::TM::Tpos, Carac::dim> TNH;
        hash_map<TNH, unsigned, MyHash, NodesEq> hm;
        apply(m->sub_mesh(Number<Carac::dim-1>()).elem_list,add_nodes(),*m,hm);
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
        m->sub_mesh(Number<1>()).elem_list.change_hash_size( *m, m->elem_list.size() /2 +1);
        m->sub_mesh(Number<2>()).elem_list.change_hash_size( *m, m->elem_list.size() /2 +1);
        if (typmat!=0 or numsst!=0 or num_proc!=0) apply(m->elem_list,apply_mat_elem(),typmat,numsst,num_proc);
        node_list_size=m->node_list.size();
        elem_list_size=m->elem_list.size();
    }
};




#endif //MESHMULTI_H

