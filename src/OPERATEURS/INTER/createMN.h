#include "../../../LMT/include/mesh/mesh.h"

using namespace LMT;

#include "../../DEFINITIONS/Interface.h"
#include "../../DEFINITIONS/Sst.h"

#include "generate_Ne.h"

/** \ingroup  Operateurs_inter 
 * \brief Calcul de la matrice elementaire Ne de sous integration et assemblage
 * 
 * Cette procedure est ecrite pour chaque type d'element possible de l'interface (à savoir Bar, Triangle ou Quad). Pour l'element considere, on selectionne le(s) element(s) du maillage de bord de la sst dont le(s) numero(s) dans la liste (elem_list) est donne par le champ e.elem. 
 * 
 * Après selection du reperage dans la matrice N globale de la matrice elementaire, on appelle la procedure add_elem_Ne() generee en à partir du script python et specifique à chaque type d'element du maillage de bord de la sst. Celle ci renvoit la matrice Ne.
 * 
 * On procède enfin à l'assemblage de cette matrice elementaire.
 */
struct add_elem_N{
    template<typename Te> void operator()(Te &e,EdgeMesh &m, SparseMatrix &N, unsigned &incr) const {
        // recherche des elements de m (maillage de bord) pointes par e.elem
        for(unsigned i=0;i<e.elem.size();++i){
            
            EdgeMesh::EA *e2 =m.elem_list[e.elem[i]]; //ElementAncestor
            
            //creation du reperage dans la matrice N des inconnues nodales du cote de la SST considere (colonnes)
            Vec<unsigned > repNc;
            repNc.resize((e2->nb_nodes_virtual())*DIM);
            for(unsigned i=0;i<e2->nb_nodes_virtual();++i){
                unsigned rep = e2->node_virtual(i)->number_in_original_mesh();
                repNc[range(DIM*i,DIM*(i+1))]=range(DIM*rep,DIM*(rep+1));
            }
            
            // creation du reperage dans la matrice N des inconnues aux noeuds equivalents de l'interface consideree (lignes)
            Vec<unsigned> repNr;
            repNr.resize(DIM);
            repNr=range(incr*DIM,(incr+1)*DIM);
            
            // assignation dans N[repNr,repNc] de Ne;
            SparseMatrix Ne;
            Ne.resize(repNr.size(),repNc.size());
            EdgeMesh::TElemList::apply_on_down_cast(e2,add_elem_Ne(),Ne);
            
            N(repNr,repNc)+=Ne;
        }
        incr+=1;
    }
};


/** \ingroup  Operateurs_inter
 * \brief Calcul de la matrice de sous-integration pour l'interface.
 * 
 * Cette procedure cree la matrice N (1 de chaque cote de l'interface) pour chaque discretisation possible et chaque type d'element de la sous-structure. La discretisation de l'interface est necessairement P0. La matrice N est donc de la taille dim*(nbre de noeuds sur le cote de la sst) x dim*(nombre d'elements de l'interface). On appelle ici pour chaque element d'interface la procedure add_elem_N.
 */
void calcN(InterfaceMesh &minter, EdgeMesh &medge, SparseMatrix &N){
    N.resize(minter.elem_list.size()*DIM,medge.node_list.size()*DIM);
    unsigned incr=0;
    apply(minter.elem_list,add_elem_N(),medge,N,incr); 
}

/** \ingroup  Operateurs_inter
 * \relates calcM
 * \brief Creation de la matrice de masse par element.
 */
struct add_elem_M{
    template<typename Te> void operator()(Te &e, SparseMatrix &M, Vec<unsigned> &indice,unsigned dim) const{
        for(unsigned l=0;l<indice.size();++l)
            M(indice[l],indice[l])=measure(e);
        indice+=dim;
    }
};

/** \ingroup  Operateurs_inter
 * \brief Calcul de la matrice de masse pour l'interface.
 * 
 * Cette procedure n'est adaptee que pour le cas d'une discretisation P0 pour l'interface. La matrice de masse est donc diagonale et est constituee de la mesure des elements sur chaque terme de la diagonale. On appelle pour chaque element d'interface la procedure add_elem_M.
 */
void calcM(InterfaceMesh &m, SparseMatrix &M){
    M.resize(m.elem_list.size()*DIM,m.elem_list.size()*DIM);
    Vec<unsigned> indice;
    indice=range(0,(int) DIM);
    apply(m.elem_list,add_elem_M(),M,indice,DIM);
}


void calcMinvN(SparseMatrix &M, SparseMatrix &N){
    // marche pour une matrice diagonale M ( P0 par element )
    SparseMatrix Minv;
    Minv.resize(M.nb_rows(),M.nb_cols());
    Minv.diag()=1./M.diag();
    N=Minv*N;
}
