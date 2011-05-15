#include <ext/hash_map>
using namespace __gnu_cxx;
using namespace std;


/** \ingroup  maillages 
\brief Structure pour le hashage du maillage cf. p_surdiscretise_SST
*/
// struct MyHash{
// template<class TV> unsigned operator()(const TV &vec1) const{
//    return vec1.num_hash;
// }
// };

/** \ingroup  maillages 
\brief Structure pour reperer les noeuds dans la hash_map
*/
// struct NodesEq{
//    template<class TV> bool operator()(const TV &n1, const TV &n2) const{
//       double eps=1e-6;
//       return (length(n1.pos-n2.pos)<=eps);
//    }
// };



/** \ingroup  maillages 
\brief Structure pour le hashage du maillage
*/
// template<class T, unsigned dim> struct Noeud_Hash{
//    Noeud_Hash(){num_hash=0;pos.set(0);}
//    Vec<T, dim> pos;
//    unsigned num_hash;
// };
// 
// template<class TE,class TM, class TH> void add_new_nodes(TE &e, TM &m, TH &hm){
//          typedef typename TM::Pvec Pvec;
//          //ajout du noeud au maillage
//          Pvec G = center(e);
//          m.add_node(G);
//          //creation de la structure pour hashage
//          Noeud_Hash<typename TE::T,TM::dim> newnoeud;
//          newnoeud.pos=G;
//          
//          //reperage des numeros des noeuds extremes de l'element
//          unsigned num_hash=0;
//          for(unsigned i=0;i<e.nb_nodes;i++)
//             num_hash+=e.node(i)->number_in_original_mesh();
//          //assignation de la valeur de hashage (somme des numeros des noeuds extremes)
//          newnoeud.num_hash=num_hash;
//          //ajout � la table de hashage de l'entree newnoeud correspondant au numero dans le maillage du dernier noeud entre.
//          hm[newnoeud]=m.node_list.size()-1;
//          //cout << "num " << m.node_list.size()-1 << " pos "<< G << endl;
// }

/** \ingroup  maillages 
\brief Ajout de noeuds sur les aretes d'un maillage
*/
// struct add_nodes {
//    template<class TE,class TM, class TH> void operator() (TE &e, TM &m, TH &hm) const{
//       typedef typename TM::Pvec Pvec;
//       //ajout du noeud au maillage
//       Pvec G = center(e);
//       m.add_node(G);
//       //creation de la structure pour hashage
//       Noeud_Hash<typename TE::T,TM::dim> newnoeud;
//       newnoeud.pos=G;
//       
//       //reperage des numeros des noeuds extremes de l'element
//       unsigned num_hash=0;
//       for(unsigned i=0;i<e.nb_nodes;i++)
//          num_hash+=e.node(i)->number_in_original_mesh();
//       //assignation de la valeur de hashage (somme des numeros des noeuds extremes)
//       newnoeud.num_hash=num_hash;
//       //ajout � la table de hashage de l'entree newnoeud correspondant au numero dans le maillage du dernier noeud entre.
//       hm[newnoeud]=m.node_list.size()-1;
//    }
// };

// #include "create_new_elem_p.h" 

/** \ingroup  Sous_structures 
\brief Modification de l'element
*/
// template<unsigned number>
// struct ModifTypeElem{
// template<class TE, class TM, class TNH> void operator()(TE &e,TM &m, hash_map<TNH, unsigned, MyHash, NodesEq> &hm, TM &m2) const{
//    unsigned nb_children = NbChildrenElement<typename TE::NE,number>::res;
// 
//    //boucle sur les elements enfants
//    Vec<unsigned> rep_nodes;
//    rep_nodes.resize(nb_children+e.nb_nodes,0);
//    for(unsigned i=0;i<e.nb_nodes;i++)
//       rep_nodes[i]=e.node(i)->number_in_original_mesh();
//    
//    for(unsigned i=0;i<nb_children;i++){
//      typename TM::EA *ea = m.get_children_of(e,LMT::Number<number>())[i]; //element ancestor i;
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


/** \ingroup  Sous_structures 
\relates p_surdiscretise_SST
\brief Proc�dure principale pour la p_surdiscr�tisation.

La proc�dure est la suivante :
- Dans un premier temps, on boucle sur le sous maillage (1 en 2D) (2 en 3D) et on cr�e les noeuds suppl�mentaires milieux des �l�ments Bar.
- Les nouveaux noeuds sonts stock�s sous la forme d'un vecteur de Noeud_Hash, structure particuli�re contenant la position des noeuds et une valeur de hashage (obtenue en sommant le num�ro des noeuds extr�mit�s de la Bar consid�r�e pour cr�er le noeud. Ceci permet ensuite de g�n�rer une Hash_map contenant comme cl� la position des noeud et renvoyant le num�ro du noeud dans le maillage. On aura donc au fur et � mesure de la cr�ation des noeuds, ajouter chaque noeud au maillage et renvoyer son num�ro dans la hash_map.
- On boucle ensuite sur tous les �l�ments du maillage et pour les enfants de niveau 1 en 2D ou 2 en 3D, on cr�e les noeuds milieux. En recherchant la position du noeud cr�� dans la hash_map, on sort la position du noeud correspondant dans le maillage. Cette technique permet d'obtenir une proc�dure tr�s rapide ne n�cessitant que peu de tests.
- Ayant ainsi rep�r� les num�ros de tous les nouveaux noeuds obtenus pour un �l�ment donn�, on peut modifier cet �l�ment. On l'�limine et on le remplace par un �l�ment de degr� plus �lev� s'appuyant sur les noeuds dont les num�ros ont �t� relev�s.
*/
template<unsigned number,class TM> void p_surdiscretise(TM &m){
   m.sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( m, m.elem_list.size() /2 +1);
   m.sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( m, m.elem_list.size() /2 +1);
   m.update_elem_children();
   m.update_elem_parents();   
   //1ere etape : ajout des noeuds et creation de la table de hashage
   typedef Noeud_Hash<typename TM::Tpos, TM::dim> TNH;
   hash_map<TNH, unsigned, MyHash, NodesEq> hm;
   apply(m.sub_mesh(LMT::Number<number>()).elem_list,add_nodes(),m,hm);
   
   //2eme etape : boucle sur les noeuds et creation de nouveaux �l�ments
   TM m2;
   for(unsigned i=0;i<m.node_list.size();i++)
      m2.add_node(m.node_list[i].pos);
      
   apply(m.elem_list,ModifTypeElem<number>(),m,hm,m2);
   m.free();
   m=m2;
}

/** \ingroup  Sous_structures 
\brief Sur discr�tisation des sous-structures : m�thode p

On modifie le type d'�l�ment du maillage de la sous-structure pour obtenir un degr� plus �lev�. Les �l�ments P1 deviennent des �l�ments P2. Ainsi sur chaque ar�te des �l�ments 2D ou 3D on vient ajouter un noeud suppl�mentaire et modifier le type d'�l�ment.

Il est donc n�cessaire lors de la compilation d'utiliser � la fois l'�l�ment 2D ou 3D g�n�rique correspondant aux maillages r�els lus (ex : Triangle) et l'�l�ment de degr� plus �lev� pour effectuer une sur-discretisation (ex : Triangle_6)

La proc�dure de modification des �l�ments est d�crite dans la fonction p_surdiscretise()
*/
struct p_surdiscretise_SST {
   template<class SST> void operator()(SST &S) const {
      p_surdiscretise<SST::dim-1>(S.mesh);
      apply(S.mesh.elem_list,apply_mat_elem(),S.typmat,S.num,S.num_proc);
   }
};


/** \ingroup  Interfaces 
\brief Extraction des �l�ments de bord correspondant aux interfaces.

En bouclant sur les �l�ments de l'interface, on rep�re l'�l�ment du maillage de peau de la sous-structure correspondante ayant le m�me centre de gravit� que l'�l�ment donn�. On cr�e alors le maillage de bord � partir des �l�ments s�lectionn�s. Le num�ro de l'�l�mentancestor cr�� est stock� dans le champ elem de l'interface.
*/
