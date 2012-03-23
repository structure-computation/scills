#include "mesh/calculate_measure.h"
#include "util/symrcm.h"
#ifndef UTILITAIRES_H
#define UTILITAIRES_H

using namespace LMT;

///creation d'une boite contenant tout le maillage donne
template<class TM>
Vec<typename TM::Pvec,2> create_box_mesh(TM &m) {
    typedef typename TM::Pvec Pvec;
    Pvec mini;
    Pvec maxi;
    Vec<Pvec,2> box;
    //creating box
    get_min_max(m.node_list,ExtractDM<pos_DM>(),mini,maxi);
    box[0]=mini;
    box[1]=maxi;
    return box;
}

template<class TV, class T,class STO, class STR> void tens(TV &n1, TV &n2, Mat<T, STO , STR>  &M){
   M.resize(n1.size(), n2.size());
   for(unsigned i=0;i<n1.size();i++)
      for(unsigned j=0;j<n2.size();j++)
         M(i,j)=n1[i]*n2[j];
}

template<class TV>
bool ComparVec(TV &V1, TV &V2) {
    unsigned nbtrue=0;
    bool res=0;
    for(unsigned i=0;i<V1.size();++i) {
        if (V1[i]<=V2[i])
            nbtrue+=1;
    }
    if (nbtrue==V1.size())
        res=1;
    std::cout << nbtrue << std::endl;
    return res;
}

template<class TV> bool inCL(Vec<TV,2> &box1, Vec<TV,2> &box2, const double &eps=1e-6) {
    TV P1=box1[0], P2=box1[1];
    TV Q1=box2[0], Q2=box2[1];

    bool testin=0;
    // verification si P1 et P2 sont dans la boite box2
    //verification si P1 est dans box2
    if ( bool_vec(P1<=Q2+eps and P1>=Q1-eps) ) {
        //verification si P2 est dans box2
        if ( bool_vec(P2<=Q2+eps and P2>=Q1-eps) )
            testin=1;
    }
    return testin;
};

/**
 *  Valable en dimension 2
 * @param box1 : vecteur contenant les coordonnees des noeuds extremite de la boite 1
 * @param box2 : vecteur contenant les coordonnees des noeuds extremite de la boite 2
 * @return : booleen =1 si la boite 1 est completement incluse dans la boite 2 (a epsilon prÃ¨s = 1e-12)
 */
template<class TT>
bool inCL(Vec<TT,4> &box1, Vec<TT,4> &box2) {
    TT eps=1e-6;
    Vec<TT,2> P1( box1[range(0,2)] );
    Vec<TT,2> P2=box1[range(2,4)];
    Vec<TT,2> Q1=box2[range(0,2)];
    Vec<TT,2> Q2=box2[range(2,4)];
    bool testin=0;

    // verification si P1 et P2 sont dans la boite box2
    //verification si P1 est dans box2
    if ( (Q1[0]-eps<= P1[0]) && (P1[0]<=Q2[0]+eps )  && (Q1[1]-eps<= P1[1]) && (P1[1]<=Q2[1]+eps ) ) {
        //verification si P2 est dans box2
        if ( (Q1[0]-eps<= P2[0]) && (P2[0]<=Q2[0]+eps )  && (Q1[1]-eps<= P2[1]) && (P2[1]<=Q2[1]+eps ) )
            testin=1;
    }
    return testin;
};

template<class TT>
bool inCL(Vec<TT,6> &box1, Vec<TT,6> &box2) {
    TT eps=1e-6;
    Vec<TT,3> P1=box1[range(0,3)];
    Vec<TT,3> P2=box1[range(3,6)];
    Vec<TT,3> Q1=box2[range(0,3)];
    Vec<TT,3> Q2=box2[range(3,6)];

    bool testin=0;
    // verification si P1 et P2 sont dans la boite box2
    //verification si P1 est dans box2
    if ( (Q1[0]-eps<= P1[0]) && (P1[0]<=Q2[0]+eps )  && (Q1[1]-eps<= P1[1]) && (P1[1]<=Q2[1]+eps ) && (Q1[2]-eps<= P1[2]) && (P1[2]<=Q2[2]+eps )) {
        //verification si P2 est dans box2
        if ( (Q1[0]-eps<= P2[0]) && (P2[0]<=Q2[0]+eps )  && (Q1[1]-eps<= P2[1]) && (P2[1]<=Q2[1]+eps )  && (Q1[2]-eps<= P2[2]) && (P2[2]<=Q2[2]+eps ))
            testin=1;
    }
    return testin;
};






template<class TV, class TM>
void assign_q_mesh(TV &q, TM &mesh) {
    unsigned dim=TM::dim;
    for(unsigned i=0;i<mesh.node_list.size();++i) {
        mesh.node_list[i].dep=q[range(i*dim,(i+1)*dim)];
    }
};

//penalisation de matrice
// a partir des ddls bloques a zero
template<class T,class Str,class Sto,class OP, class TV>
void penalisation(Mat<T,Str,Sto,OP>&K, TV &v, T &coef ) {
    coef=1e5*max(K.diag());
    for(unsigned i=0;i<v.size();++i)
        K(v[i],v[i])+= coef;
}

// a partir des ddls bloques, de la matrice A de relation entre ces ddls
template<class T,class Sto,class OP, class TV, class TMAT2>
void penalisation(Mat<T,Gen<>,Sto,OP> &K, TV &v, TMAT2 &A, T &coef ) {
    coef=1e5*max(A.diag());
    K(v,v)+=coef*trans(A)*A;
}

template<class T,class Sto,class OP, class TV, class TMAT2>
void penalisation(Mat<T,Sym<>,Sto,OP> &K, TV &v, TMAT2 &A, T &coef ) {
    coef=1e5*max(A.diag());
    K[v]+=coef*trans(A)*A;
}

//penalisation de vecteur
// a partir des ddls bloques a zero
template<class TV1, class TV2>
void penalisation(TV1 &v, TV2 &rep, typename TV1::T &coef) {
    for(unsigned i=0;i<v.size();++i)
        v(rep[i])+= coef;

}

template<class TT,int static_size_>
Vec<TT,static_size_,void> ProjT(Vec<TT,static_size_,void> &v1, Vec<TT,static_size_,void> &n) {
    Vec<TT,static_size_,void>  res;
    res=v1-dot(v1,n)*n;
    return res;
}


template<class TT,int static_size_>
TT norm_2(Vec<TT,static_size_,void> &v) {
    TT res;
    res=std::sqrt(dot(v,v));
    return res;
}



// ecriture d'une matrice dans un fichier texte
template<class TMAT>
void write_matrix(TMAT &M, const Sc2String &name) {
    std::ofstream f(name.c_str());
    for(unsigned i=0;i<M.nb_rows();++i) {
        for(unsigned j=0;j<M.nb_cols();++j) {
            typename TMAT::T a;
            a=M(i,j);
            f<< a<< " ";
        }
        f << "\n";
    }
}

struct add_elem_node {
    template<class TE, class TMAT>
    void operator()(TE &e, TMAT &K) const {
        Vec<unsigned> rep;
        rep.resize(e.nb_nodes);
        for(unsigned i=0;i<e.nb_nodes;++i) {
            rep[i]=e.node(i)->number_in_original_mesh();
        }
        K(rep,rep)+=1;
    }

};

// reorganisation des noeuds d'un maillage en fonction d'une renumerotation Cuthin Mackee
template<class TM>
void reordering_nodes(TM &m) {

    Mat<bool, Gen<>, SparseLine<> > K;
    K.resize(m.node_list.size(),m.node_list.size());

    // creation de la matrice d'indices
    apply(m.elem_list,add_elem_node(),K);

    // renumerotation de la matrice
    Vec<Vec<unsigned> > indices;
    for(unsigned i=0;i<K.nb_rows();++i)
        indices.push_back( K.data[i].indices );
    #ifndef INCOMPLETE
         Vec<unsigned> res = symrcm( indices );
    #else  
         Vec<unsigned> res = symamd( indices );
    #endif
   
    Vec<typename TM::TNode *> perm;
    perm.resize( res.size() );
    Vec<typename TM::TNode> nn;
    nn.resize( res.size() );
    for(unsigned i=0;i<res.size();++i) {
        perm[ res[i] ] = &m.node_list[ i ];
        nn[ i ] = m.node_list[ res[i] ];
    }
    m.change_node_ptr_in_elements( perm.ptr() );

    for(unsigned i=0;i<res.size();++i) {
        m.node_list[i] = nn[i];
        m.node_list[i].number = i;
    }

    m.signal_connectivity_change();

    /*Mat<bool, Gen<>, SparseLine<> > K2;
     K2.resize(m.node_list.size(),m.node_list.size());
     
     // creation de la matrice d'indices
     apply(m.elem_list,add_elem_node(),K2);
     
      */
}

struct reordering_nodes_SST {
    template<class SST>
    void operator()(SST &S) const {
        //reordering_nodes(S.mesh);
        xyz_ordering(S.mesh);
        //node_reorder( S.mesh );
    }
};


/** \ingroup  Maillages 
\brief Recherche si un pt est a l'interieur d'une boite définie par ses deux extrémités et par une base dont le dernier vecteur est la normale. Une procédure est écrite en 2d et une autre en 3d

Pour le cas 2d, on cherche à déterminer si le pt M est dans le segment [AB]. Pour cela, on calcul les distances MA et MB puis on vérifie si \f$ \frac{\vec{MA}.\vec{MB}}{MA MB} =-1 \f$ à epsilon près.

Pour le cas 3d, connaissant la base dont le dernier vecteur est la normale à la surface et les deux autres sont parallèles aux cotés, on cherche à savoir si le pt M est dans le rectangle défini par les deux points extrémités de la boite. On calcule donc :
\f$ (\vec{MA}.\vec{n})x (\vec{MB}.\vec{n}) <=0 \f$
\f$ (\vec{MA}.\vec{t1})x (\vec{MB}.\vec{t1}) <=0 \f$
\f$ (\vec{MA}.\vec{t2})x (\vec{MB}.\vec{t2}) <=0 \f$
Attention il faut nécessairement que t1 et t2 soient parallèles aux cotés du rectangle (à construire lorsque l'on crée la base). 
*/
template<class T> bool pt_in_box(Vec<T,2> &pt, Vec<Vec<T,2>,2> &box, Vec<Vec<T,2>,2> &base, const double eps=1e-6){
   T la = length(pt-box[0]);
   T lb = length(pt-box[1]);
   bool flag=0;
   if(la<=eps or lb<=eps or LMT::abs(dot(pt-box[0],pt-box[1])+la*lb)<=eps)
         flag=1;
   return flag;
}
template<class T> bool pt_in_box(Vec<T,3> &pt, Vec<Vec<T,3>,2> &box, Vec<Vec<T,3>,3> &base, const double eps=1e-6){
   bool flag=0;
   if( (dot(pt-box[0],base[2])*dot(pt-box[1],base[2])<=eps) and (dot(pt-box[0],base[0])*dot(pt-box[1],base[0])<=eps) and (dot(pt-box[0],base[1])*dot(pt-box[1],base[1])<=eps) )
         flag=1;     
   return flag;
}

/** \ingroup  Maillages 
\brief Recherche si un pt est a l'interieur d'une boite définie par ses deux extrémités (rectangle en 2d et parallélépipède en 3d)

Pour le cas 2d ou 3d, on cherche à savoir si le pt M est dans le rectangle ou (parallélépipède) défini par les deux points extrémités de la boite A et B. On calcule donc :
\f$ (\vec{MA}.\vec{x})x (\vec{MB}.\vec{x}) <=0 \f$
\f$ (\vec{MA}.\vec{y})x (\vec{MB}.\vec{y}) <=0 \f$
\f$ (\vec{MA}.\vec{z})x (\vec{MB}.\vec{z}) <=0 \f$
Attention il faut nécessairement que t1 et t2 soient parallèles aux cotés du rectangle (à construire lorsque l'on crée la base). 
*/
template<class T> bool pt_in_box(Vec<T,2> &pt, Vec<Vec<T,2>,2> &box, const double eps=1e-6){
   bool flag=0;
   Vec<T,2> x(1.,0.),y(0.,1.);
   if( (dot(pt-box[0],x)*dot(pt-box[1],x)<=eps) and (dot(pt-box[0],y)*dot(pt-box[1],y)<=eps) )
         flag=1;
   return flag;
}
template<class T> bool pt_in_box(Vec<T,3> &pt, Vec<Vec<T,3>,2> &box,  const double eps=1e-6){
   bool flag=0;
   Vec<T,3> x(1.,0.,0.),y(0.,1.,0.),z(0.,0.,1.);   
   if( (dot(pt-box[0],x)*dot(pt-box[1],x)<=eps) and (dot(pt-box[0],y)*dot(pt-box[1],y)<=eps) and (dot(pt-box[0],z)*dot(pt-box[1],z)<=eps) )
         flag=1;     
   return flag;
}
template<class T> bool pt_in_box(const Vec<T,2> &pt, Vec<Vec<T,2>,2> &box, const double eps=1e-6){
   bool flag=0;
   Vec<T,2> x(1.,0.),y(0.,1.);
   if( (dot(pt-box[0],x)*dot(pt-box[1],x)<=eps) and (dot(pt-box[0],y)*dot(pt-box[1],y)<=eps) )
         flag=1;
   return flag;
}
template<class T> bool pt_in_box(const Vec<T,3> &pt, Vec<Vec<T,3>,2> &box,  const double eps=1e-6){
   bool flag=0;
   Vec<T,3> x(1.,0.,0.),y(0.,1.,0.),z(0.,0.,1.);   
   if( (dot(pt-box[0],x)*dot(pt-box[1],x)<=eps) and (dot(pt-box[0],y)*dot(pt-box[1],y)<=eps) and (dot(pt-box[0],z)*dot(pt-box[1],z)<=eps) )
         flag=1;     
   return flag;
}


/** \brief Intersection entre deux boites definies par leurs noeuds extremite. Cette fonction renvoit 1 si les deux boites s'intersectent ou se touchent (à epsilon près), 0 sinon

Une fonction est écrite pour des boites 2d et une autre pour des boites 3d. Normalement aucune restriction n'est nécessaire dans l'ordre des points extrémités des boites. Les boites doivent simplement pouvoir être représentées par 2 points et donc alignées sur les axes x, y, z. Cette fonction ne nécessite en 3d que 3 comparaisons et 2 tests ainsi que l'utilisation de produits scalaires qui sont normalement optimisés.
**/
template<class T> bool intersection_box(Vec<Vec<T,2>,2> &box1, Vec<Vec<T,2>,2> &box2, const double &eps=1e-6){
   typedef Vec<T,2> TV;
   //boite 1
   TV P1 = (box1[0]+box1[1])/2.;
   TV E1 = LMT::abs(box1[1]-box1[0])/2.;
   //boite 2
   TV P2 = (box2[0]+box2[1])/2.;
   TV E2 = LMT::abs(box2[1]-box2[0])/2.;

   TV V = P2-P1;
   TV x(1,0);
   TV y(0,1);
   
   return ( (LMT::abs(dot(V,x))<=dot(E1,x)+dot(E2,x)+eps) and (LMT::abs(dot(V,y))<=dot(E1,y)+dot(E2,y)+eps) ) ;
}
template<class T> bool intersection_box(Vec<Vec<T,3>,2> &box1, Vec<Vec<T,3>,2> &box2, const double &eps=1e-6){
   typedef Vec<T,3> TV;
   //boite 1
   TV P1 = (box1[0]+box1[1])/2.;
   TV E1 = LMT::abs(box1[1]-box1[0])/2.;
   //boite 2
   TV P2 = (box2[0]+box2[1])/2.;
   TV E2 = LMT::abs(box2[1]-box2[0])/2.;

   TV V = P2-P1;
   TV x(1,0,0);
   TV y(0,1,0);
   TV z(0,0,1);
   
   return ( (LMT::abs(dot(V,x))<=dot(E1,x)+dot(E2,x)+eps) and (LMT::abs(dot(V,y))<=dot(E1,y)+dot(E2,y)+eps) and (LMT::abs(dot(V,z))<=dot(E1,z)+dot(E2,z)+eps) );
}


/** \ingroup  Maillages 
\brief Translation des noeuds d'un maillage

On passe le maillage et un vecteur de translation pour modifier la position des noeuds du maillage
*/
template<class TM,class TV> void translate_mesh(TM &m, const TV &trans){
   for(unsigned i=0;i<m.node_list.size();i++)
      m.node_list[i].pos+=trans;
}
template<class TM,class TV> void translate_mesh(TM &m, TV &trans){
   for(unsigned i=0;i<m.node_list.size();i++)
      m.node_list[i].pos+=trans;
}

#endif //UTILITAIRES_H
