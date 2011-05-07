using namespace LMT;
using namespace std;


//******************************************************************************
// Procedures de redecoupage des elements et construction des noeuds equivalents
//******************************************************************************

// sur-discretisation
template<class TM>
struct DivideOnSkin {
    template<class TE> bool operator()(TE &e) const {
        bool is_on_skin = false;
        typename TM::TNode *nn[ NbChildrenElement<typename TE::NE,1>::res ];
        for(unsigned i=0;i<NbChildrenElement<typename TE::NE,1>::res;++i) {
            if ( m->sub_mesh(Number<1>()).get_parents_of_EA( m->get_children_of(e,Number<1>())[i] ).size()==1 ) {
                is_on_skin = true;
                
                nn[i] = m->add_node();
                
                typename TM::EA *child_elem = m->get_children_of(e,Number<1>())[i];
                
                Vec<std::pair<typename TM::Tpos,const typename TM::TNode *> > pond_list;
                for(unsigned j=0;j<child_elem->nb_nodes_virtual();++j)
                    pond_list.push_back( std::pair<typename TM::Tpos,const typename TM::TNode *>( 1.0/child_elem->nb_nodes_virtual(), child_elem->node_virtual(j) ) );
                DM::get_ponderation( pond_list.ptr(), pond_list.size(), *nn[i] );
            }
            else
                nn[i] = 0;
            
        }
        if ( is_on_skin )
            divide_element_using_elem_children( e, *m, nn );
        //return false;
        return is_on_skin;
    }
    TM *m;
};

template<class TM,class DD>
struct MarkPos {
    template<class TE> void operator()(TE &e) const {
        m->elem_list.get_data( *d, e ) = center(e)[0];
    }
    DD *d;
    TM *m;
};

/**\ingroup Sous_structures
\relates surdiscretise_SST
\brief Sur-discretisation du maillage des Ssts.
*/
template<class TM> void surdiscretise(TM &m) {
    m.update_elem_parents(Number<1>());
    DivideOnSkin<TM> dos; dos.m = &m;
    m.remove_elements_if( dos );
}

/**\ingroup Sous_structures
\brief Sur-discretisation du maillage des Ssts.

Tous les éléments de bord de la sous-structure sont découpés de manière automatique en fonction de leur configuration (1, 2 ou plusieurs faces à découper). On appelle ici la fonction surdiscretise().

Une autre méthode en cours de développement est de raffiner les faces de bord, ceci permet de garder la même forme d'élément.
*/
struct surdiscretise_SST {
   template<class SST> void operator()(SST &S) const {
      surdiscretise(*S.mesh.m);
   }
};

// noeuds equivalents sur l'interface
struct createelem
{
  // Bar
  template<class TN,class TNG,class TD,unsigned NET,class TM> void operator()(Element<Bar,TN,TNG,TD,NET> &e, TM &m) const
  {
    // add center node
    unsigned nbnodes=m.node_list.size();
    m.add_node( center(e) );
    Vec<unsigned > rep;
    for(unsigned j=0;j<e.nb_nodes;++j)
      rep.push_back(e.node(j)->number_in_original_mesh() );
    e.elem.resize(2);
    unsigned nsupp=m.node_list[nbnodes].number_in_original_mesh();
    m.add_element(Bar() ,DefaultBehavior(),&m.node_list[rep[0]],&m.node_list[nsupp]);
    e.elem[0]=m.elem_list.size()-1;
    m.add_element(Bar() ,DefaultBehavior(), &m.node_list[nsupp],&m.node_list[rep[1]]);
    e.elem[1]=m.elem_list.size()-1;
  }
  // Triangle
  template<class TN,class TNG,class TD,unsigned NET,class TM> void operator()(Element<Triangle,TN,TNG,TD,NET> &e, TM &m) const
  {
    // add center node
    unsigned nbnodes=m.node_list.size();
    m.add_node( center(e) );
    Vec<unsigned > rep;
    for(unsigned j=0;j<e.nb_nodes;++j)
      rep.push_back(e.node(j)->number_in_original_mesh() );
    e.elem.resize(3);
    unsigned nsupp=m.node_list[nbnodes].number_in_original_mesh();
    m.add_element(Triangle() ,DefaultBehavior(),&m.node_list[rep[0]],&m.node_list[rep[1]],&m.node_list[nsupp]);
    e.elem[0]=m.elem_list.size()-1;
    m.add_element(Triangle() ,DefaultBehavior(), &m.node_list[rep[1]],&m.node_list[rep[2]],&m.node_list[nsupp]);
    e.elem[1]=m.elem_list.size()-1;
    m.add_element(Triangle() ,DefaultBehavior(), &m.node_list[rep[2]],&m.node_list[rep[0]],&m.node_list[nsupp]);
    e.elem[2]=m.elem_list.size()-1;
  }
  // Par defaut
  template<class TE,class TM> void operator()(TE &e, TM &m) const
  {
   cout << "Surdiscretisation de type h non implementee pour les triangle_6 quad quad_8 et bar_3" << endl;
   assert(0); 
  }
};

struct addnode_SSTedge
{
  template<class TN,class TM> void operator()(TN &n,TM &m) const{m.add_node(n.pos);}
};

struct createnodeeq
{
  template<class TE> void operator()(TE &e, Vec< Vec<typename TE::T,TE::TNode::dim> > &nodeeq) const {  nodeeq.push_back(center(e)); }
};


/** \ingroup  Maillage_geometrie 
\brief Création des maillages des bords des sous-structures et création des noeuds équivalents
*/
struct decoup_INTER
{
  template<class INTER,class TV1> void operator()(INTER &inter, TV1 &S) const{
    for(unsigned j=0;j<inter.side.size();++j){
      int ii=inter.side[j].vois[0];
      int jj=inter.side[j].vois[1];
      // recopie des noeuds du maillage
      apply(inter.side[j].mesh->node_list,addnode_SSTedge(),*S[ii].edge[jj].mesh);
      apply(inter.side[j].mesh->elem_list,createelem(),*S[ii].edge[jj].mesh);
      //creation des noeuds equivalents de l'interface
      apply(inter.side[j].mesh->elem_list,createnodeeq(),inter.side[j].nodeeq);

    }
  }
};


struct copyelem{
  template<class TE> void operator()(TE &e,unsigned &incr) const{
    e.elem.resize(1);
    e.elem[0]=incr;
    incr+=1;
    }
};

// copie des noeuds des interfaces sur les bords des SST
/** \ingroup  Maillage_geometrie 
\brief Création des maillages des bords des sous-structures par copie des maillages d'interface et création des noeuds équivalents 
*/
struct copy_INTER{
  template<class INTER,class TV1> void operator()(INTER &inter, TV1 &S) const{
    for(unsigned j=0;j<inter.side.size();++j){
      int ii=inter.side[j].vois[0];
      int jj=inter.side[j].vois[1];

      // recopie des noeuds du maillage sur les bord des SSTs
//       S[ii].edge[jj].mesh.append(inter.side[j].mesh);
//       unsigned incr=0;
//       apply(inter.side[j].mesh.elem_list,copyelem(),incr);
      
      
/*#ifdef PRINT_ALLOC
      total_allocated[ typeid(typename TV1::template SubType<0>::T::TMESHedge).name() ] += sizeof(typename TV1::template SubType<0>::T::TMESHedge);
#endif
      S[ii].edge[jj].mesh=new typename TV1::template SubType<0>::T::TMESHedge;
      S[ii].edge[jj].mesh->append(*inter.side[j].mesh);
      S[ii].edge[jj].mesh->sub_mesh(Number<1>()).elem_list.change_hash_size( *S[ii].edge[jj].mesh,1);
      S[ii].edge[jj].mesh->sub_mesh(Number<2>()).elem_list.change_hash_size( *S[ii].edge[jj].mesh,1);
      S[ii].edge[jj].mesh->sub_mesh(Number<0>()).elem_list.change_hash_size( *S[ii].edge[jj].mesh,1);
      S[ii].edge[jj].mesh->elem_list.change_hash_size( *S[ii].edge[jj].mesh,1);*/
      
      //on copie via le pointeur
      S[ii].edge[jj].mesh=inter.side[j].mesh;
      
      // recopie des noeuds du maillage sur les bord des SSTs
      unsigned incr=0;
      if (j==0 or inter.comp=="Contact_jeu_physique" or inter.comp=="periodique") apply(inter.side[j].mesh->elem_list,copyelem(),incr);

      //creation des noeuds equivalents de l'interface
      if (j==0 or inter.comp=="Contact_jeu_physique" or inter.comp=="periodique") apply(inter.side[j].mesh->elem_list,createnodeeq(),inter.side[j].nodeeq);
      else inter.side[1].nodeeq=inter.side[0].nodeeq;

    }
  }
};

//  TODO : renumerotation des noeuds du maillage : a faire apres surdiscretisation
struct renumerotation_noeuds {
  template<class SST> void operator()(SST &S) const { xyz_ordering(S.mesh);}
};
