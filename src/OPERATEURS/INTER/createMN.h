using namespace LMT;
using namespace std;

#include "generate_Ne.h"

/** \ingroup  Operateurs_inter 
\brief Calcul de la matrice �l�mentaire Ne de sous integration et assemblage

Cette proc�dure est �crite pour chaque type d'�l�ment possible de l'interface (� savoir Bar, Triangle ou Quad). Pour l'�l�ment consid�r�, on s�lectionne le(s) �l�ment(s) du maillage de bord de la sst dont le(s) num�ro(s) dans la liste (elem_list) est donn� par le champ e.elem. 

Apr�s s�lection du rep�rage dans la matrice N globale de la matrice �l�mentaire, on appelle la proc�dure add_elem_Ne() g�n�r�e en � partir du script python et sp�cifique � chaque type d'�l�ment du maillage de bord de la sst. Celle ci renvoit la matrice Ne.

On proc�de enfin � l'assemblage de cette matrice �l�mentaire.
*/
struct add_elem_N{
template<class TE, class TM, class TMAT> void operator()(TE &e,TM &m,TMAT &N, unsigned &incr) const {
      // recherche des elements de m (maillage de bord) pointes par e.elem
      for(unsigned i=0;i<e.elem.size();++i){

         typename TM::EA *e2 =m.elem_list[e.elem[i]]; //ElementAncestor

         //creation du reperage dans la matrice N des inconnues nodales du cote de la SST considere (colonnes)
         Vec<unsigned > repNc;
         repNc.resize((e2->nb_nodes_virtual())*TM::dim);
         for(unsigned i=0;i<e2->nb_nodes_virtual();++i){
            unsigned rep = e2->node_virtual(i)->number_in_original_mesh();
            repNc[range(TM::dim*i,TM::dim*(i+1))]=range(TM::dim*rep,TM::dim*(rep+1));
         }

         // creation du reperage dans la matrice N des inconnues aux noeuds equivalents de l'interface consideree (lignes)
         Vec<unsigned> repNr;
         repNr.resize(TM::dim);
         repNr=range(incr*TM::dim,(incr+1)*TM::dim);
  
         // assignation dans N[repNr,repNc] de Ne;
         Mat<typename TMAT::T,Gen<>, SparseLine<> > Ne;
         Ne.resize(repNr.size(),repNc.size());
         TM::TElemList::apply_on_down_cast(e2,add_elem_Ne(),Ne);

         N(repNr,repNc)+=Ne;
     }
     incr+=1;
   }
};


/** \ingroup  Operateurs_inter
\brief Calcul de la matrice de sous-int�gration pour l'interface.

Cette proc�dure cr�e la matrice N (1 de chaque cote de l'interface) pour chaque discretisation possible et chaque type d'�l�ment de la sous-structure. La discr�tisation de l'interface est n�cessairement P0. La matrice N est donc de la taille dim*(nbre de noeuds sur le cote de la sst) x dim*(nombre d'�l�ments de l'interface). On appelle ici pour chaque �l�ment d'interface la proc�dure add_elem_N.
*/
template<class TM1, class TM2, class TMAT>
void calcN(TM1 &minter, TM2 &medge, TMAT &N){
    N.resize(minter.elem_list.size()*TM1::dim,medge.node_list.size()*TM1::dim);
    unsigned incr=0;
    apply(minter.elem_list,add_elem_N(),medge,N,incr); 
};

/** \ingroup  Operateurs_inter
\relates calcM
\brief Creation de la matrice de masse par element.
*/
struct add_elem_M{
   template<class TE, class TMAT> void operator()(TE &e, TMAT &M, Vec<unsigned> &indice,unsigned dim) const{
   for(unsigned l=0;l<indice.size();++l)
      M(indice[l],indice[l])=measure(e);
   indice+=dim;
   }
};

/** \ingroup  Operateurs_inter
\brief Calcul de la matrice de masse pour l'interface.

Cette proc�dure n'est adapt�e que pour le cas d'une discr�tisation P0 pour l'interface. La matrice de masse est donc diagonale et est constitu�e de la mesure des �l�ments sur chaque terme de la diagonale. On appelle pour chaque �l�ment d'interface la proc�dure add_elem_M.
*/
template<class TM1, class TMAT>
void calcM(TM1 &m, TMAT &M){
    M.resize(m.elem_list.size()*TM1::dim,m.elem_list.size()*TM1::dim);
    Vec<unsigned> indice;
    indice=range(0,(int)TM1::dim);
    unsigned dim=TM1::dim;
    apply(m.elem_list,add_elem_M(),M,indice,dim);
};


template<class TMAT>
void calcMinvN(TMAT &M, TMAT &N){
    // marche pour une matrice diagonale M ( P0 par element )
    TMAT Minv;
    Minv.resize(M.nb_rows(),M.nb_cols());
    Minv.diag()=1./M.diag();
    N=Minv*N;
//     for(unsigned i=0;i<M.nb_rows();++i){
//       for(unsigned j=0;j<N.nb_cols();++j){
//          N(i,j)=N(i,j)/M(i,i);
//       }
//     }
    //cout << N << endl;
};
