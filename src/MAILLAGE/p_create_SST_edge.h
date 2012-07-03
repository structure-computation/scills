#include "create_new_elem_p.h"
#include "../DEFINITIONS/structure_typedef.h"

/** \ingroup  Sous_structures 
\brief Modification de l'element
*/
// struct ModifTypeElemInter{
// template<class TE, class TM, class TNH> void operator()(TE &e, TM &m, hash_map<TNH, unsigned, MyHash, NodesEq> &hm, TM &m2) const{
//    unsigned nb_children = NbChildrenElement<typename TE::NE,dim-2>::res;
// 
//    //boucle sur les elements enfants
//    Vec<unsigned> rep_nodes;
//    rep_nodes.resize(nb_children+e.nb_nodes,0);
//    for(unsigned i=0;i<e.nb_nodes;i++)
//       rep_nodes[i]=e.node(i)->number_in_original_mesh();
//    
//    for(unsigned i=0;i<nb_children;i++){
//      typename TM::EA *ea = m.get_children_of(e,LMT::Number<dim-2>())[i]; //element ancestor i;
//      typename TM::Pvec G = center(*ea);
// 
//      Noeud_Hash<typename TE::T,TM::dim> newnoeud;
//      newnoeud.pos=G;
//      unsigned num_hash=0;
//      for(unsigned j=0;j<ea->nb_nodes_virtual();j++)
//         num_hash+=ea->node_virtual(j)->number_in_original_mesh();
//         
//      //assignation de la valeur de hashage (somme des numeros des noeuds extremes)
//      newnoeud.num_hash=num_hash;
//      rep_nodes[i+e.nb_nodes]=hm[newnoeud];    
//    }
// 
//    create_new_elem_p(e,m2,rep_nodes);
// }
// };

struct ModifTypeElemInter{
template<class TE, class TM, class TNH> void operator()(TE &e, TM &m, hash_map<TNH, unsigned, MyHash, NodesEq> &hm, TM &m2) const{
   //boucle sur les elements 
   Vec<unsigned> rep_nodes;
   rep_nodes.resize(e.nb_nodes+1,0);
   for(unsigned i=0;i<e.nb_nodes;i++)
      rep_nodes[i]=e.node(i)->number_in_original_mesh();
   
   typename TM::Pvec G = center(e);
   Noeud_Hash newnoeud;
   newnoeud.pos=G;
   unsigned num_hash=0;
   for(unsigned j=0;j<e.nb_nodes;j++)
        num_hash+=e.node(j)->number_in_original_mesh(); 
   newnoeud.num_hash=num_hash;
   rep_nodes[e.nb_nodes]=hm[newnoeud];         

   create_new_elem_p(e,m2,rep_nodes);
}
};

//pour les elements triangle ou quad (3D)
template<class TM> void create_new_mesh_p_3D(TM &m, TM &m2){
   m2.append(m);
   p_surdiscretise<TM::dim-2>(m2);
}

//pour les elements bar (2D)
template<class TM> void create_new_mesh_p_2D(TM &m, TM &m2){
    m.sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( m, m.elem_list.size() /2 +1);
    m.sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( m, m.elem_list.size() /2 +1);
    m.update_elem_children();
   m.update_elem_parents();
   //1ere etape : ajout des noeuds et creation de la table de hashage
   hash_map<Noeud_Hash, unsigned, MyHash, NodesEq> hm;
   apply(m.elem_list,add_nodes(),m,hm);

   //2eme etape : boucle sur les noeuds et creation de nouveaux éléments
   for(unsigned i=0;i<m.node_list.size();i++)
      m2.add_node(m.node_list[i].pos);
   
   apply(m.elem_list,ModifTypeElemInter(),m,hm,m2);

}




/// assignation du numero de l'element du bord de la sst correspondant a un element du maillage d'interface ( avec verification si le cdg est le meme)
struct copy_num_elem_verif{
  template<class TE, class TM> void operator()(TE &e,unsigned &incr, TM &m) const{
    //verification si l'element correspondant au numero donne par incr du maillage m a la meme cdg que l'element e
    TYPEREEL eps=1e-8;
    if( length(center(e)-center(*m.elem_list[incr]))<=eps ) {
      e.elem.resize(1);
      e.elem[0]=incr;
      incr+=1;
    }
    else{
      std::cout << "Les elements du maillage de bord ne sont pas numerotes de la meme facon que ceux de l'interface : TODO " << endl;
      assert(0);
    }
  }
};

/// copie des noeuds des interfaces sur les bords des SST + assignation du numero de l'element peau de bord de la sous-structure correspondant a chaque element de l'interface
struct p_decoup_INTER{
  void operator()(Interface &inter, Substructures &S) const{
    for(unsigned j=0;j<inter.side.size();++j){
      int ii=inter.side[j].vois[0];
      int jj=inter.side[j].vois[1];

      if (j == 0 or inter.comp=="Contact_jeu_physique" or inter.comp=="periodique") {
#ifdef PRINT_ALLOC
      //total_allocated[ typeid(typename TV1::template SubType<0>::T::TMESHedge).name() ] += sizeof(typename TV1::template SubType<0>::T::TMESHedge);
#endif
      S[ii].edge[jj].mesh=new EdgeMesh;
      S[ii].edge[jj].mesh->sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( *S[ii].edge[jj].mesh,1);
      S[ii].edge[jj].mesh->sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( *S[ii].edge[jj].mesh,1);
      S[ii].edge[jj].mesh->sub_mesh(LMT::Number<0>()).elem_list.change_hash_size( *S[ii].edge[jj].mesh,1);
      S[ii].edge[jj].mesh->elem_list.change_hash_size( *S[ii].edge[jj].mesh,1);
      
      // creation du maillage de bord de la sous-structure correspondant
      if (DIM==3)
      create_new_mesh_p_3D(*inter.side[j].mesh,*S[ii].edge[jj].mesh);
      else if(DIM==2)
      create_new_mesh_p_2D(*inter.side[j].mesh,*S[ii].edge[jj].mesh);
      else{ std::cout << "pas de multiechelle en 1D" << endl;
         assert(0);}
      
         unsigned incr=0;
         apply(inter.side[j].mesh->elem_list,copy_num_elem_verif(),incr,*S[ii].edge[jj].mesh);
         
         //creation des noeuds equivalents de l'interface
         apply(inter.side[j].mesh->elem_list,createnodeeq(),inter.side[j].nodeeq);

      } else {
          //on copie juste le pointeur vers l autre edge
          S[ii].edge[jj].mesh=S[inter.side[0].vois[0]].edge[inter.side[0].vois[1]].mesh;
          inter.side[1].nodeeq=inter.side[0].nodeeq;
      }
    }
  }
};

