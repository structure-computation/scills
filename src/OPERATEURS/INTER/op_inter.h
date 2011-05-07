
using namespace LMT;
using namespace std;

//***********************************************************
// Calcul des matrices de masse M et sous-integration N
//***********************************************************
#include "createMN.h"
#include "mesh/remove_doubles.h"
#include "containers/evaluate_nb_cycles.h"


//********************************************
// Calcul matrice de masse et sous-integration
//*******************************************
/** \ingroup Operateurs_inter
\brief Construction de la matrice de masse et de sous-int�gration (Interface::Side::M et Interface::Side::N)
 
On construit la matrice de masse Interface::Side::M par la fonction calcM() puis la matrice de sous-int�gration par calcN(). 
- La matrice de masse permet l'int�gration de F*W sur une interface et est d�finie de fa�on �l�mentaire par : \f$ M_{ij} = \int_{\Gamma} \psi_i \psi_j dS \f$ o� \f$ \psi_i \f$ est la fonction de forme associ�e aux �l�ments d'interface. Par d�faut, ces fonction de forme sont contantes par �l�ment (P0).
 
- La matrice de sous-int�gration N permet la projection de W (interface) � U restriction du d�placement sur le bord d'une sous-structure. Par d�faut, l'interpolation de la sous-structure est de degr� P1. La matrice �l�mentaire \f$ N_e \f$ est donc d�finie par : \f$ N_{e_{ij}}=\int_{\Gamma} \psi_i \phi_j \f$ o� \f$ \phi_j \f$ est la fonction de forme associ�e au noeud i de l'�l�ment de bord de la sous-structure. 
De plus pour �viter des oscillations parasites (notamment observables sur des interfaces de contact ou entre des plis orthotropes...), on assemble les matrices �l�mentaires \f$ N_e \f$ pour un ensemble d'�l�ments correspondant � un seul �l�ment d'interface.
 
Une autre mani�re serait d'avoir un degr� P2 sur le bord de la sous-structure. La sous-int�gration permettrait alors directement d'�viter les oscillations. (en cours de d�veloppement). 
 
Enfin la v�ritable matrice N est telle que W = N U, on multiplie donc N par \f$ M^{-1} \f$. On calcule en prime N transpos� pour acc�lerer les calculs par la suite.
 
*/
struct CalcMN {
    template<class TV1, class INTER >
    void operator()(INTER &Inter, TV1 &S) const {
        for(unsigned j=0;j<Inter.side.size();++j) {
            int ii=Inter.side[j].vois[0];
            int jj=Inter.side[j].vois[1];
            calcM(*Inter.side[j].mesh,Inter.side[j].M);
            calcN(*Inter.side[j].mesh,*S[ii].edge[jj].mesh , Inter.side[j].N);
            calcMinvN(Inter.side[j].M,Inter.side[j].N);
            Inter.side[j].Nt=trans(Inter.side[j].N);
            //testMN(Inter.side[j].N,S[ii].edge[jj].mesh,Inter.side[j].nodeeq,Inter);
        }
    }
};

template<class TMAT, class TM, class TVe, class INTER>
void testMN(TMAT &N,TM &mesh,TVe &nodeeq, INTER &inter) {
    unsigned dim = INTER::dim;
    typedef Vec<typename INTER::T> TV;

    TV xbord;
    xbord.resize(N.nb_cols(),0.0);
    for(unsigned i=0;i<mesh.node_list.size();++i) {
        xbord[range(i*dim,(i+1)*dim)]=mesh.node_list[i].pos;
        cout << mesh.node_list[i].pos << endl;
    }
    cout << "\n" << endl;
    cout << "xbord" << xbord <<endl;

    TV xinter;
    xinter.resize(N.nb_rows(),0.0);
    for(unsigned i=0;i<nodeeq.size();++i) {
        xinter[range(i*dim,(i+1)*dim)]=nodeeq[i];
    }
    //
    cout << "xinter " << xinter << endl;
    cout <<  "Nxbord-xinter" << N*xbord-xinter << endl;


};

//********************************************
// Correspondance des ddls d'interface de chaque cote
//*******************************************
/** \ingroup Operateurs_inter
\brief  Creation des listes de correspondances des ddls d'interface - compatibilite. (Interface::Side::ddlcorresp)
 
Ne sert pas r�ellement m�me si les correspondances sont utilis�es dans l'�tape locale, car la mani�re de construire les interfaces �vite un ordre diff�rent pour les noeuds des interfaces (cf. \ref Maillage_geometrie  ).
 
Cette fonction est par contre n�cessairement utile pour des maillages d'interface disjoints.  
*/
struct  Corresp_ddlinter {
    template<class INTER >
    void operator()(INTER &Inter) const {
        //cout << Inter.vois << endl;
        if (Inter.side.size()==2) {
            if (Inter.side[0].nodeeq.size() != Inter.side[1].nodeeq.size()) {
              cout << "attention maillages incompatibles - non implemente (nbnodeeq) - inter : " << Inter.num << " voisin " << Inter.vois << endl;
                assert(0);
            }

            Inter.side[0].ddlcorresp.resize(Inter.side[0].mesh->elem_list.size()*INTER::dim);
            Inter.side[1].ddlcorresp.resize(Inter.side[1].mesh->elem_list.size()*INTER::dim);
            Inter.side[0].ddlcorresp=range(0,int(Inter.side[1].nodeeq.size()*INTER::dim));
            Inter.side[1].ddlcorresp.set(0);

            int test=0;
            Vec<unsigned> candidats=range(Inter.side[0].mesh->elem_list.size());
            for (unsigned i=0;i<Inter.side[0].mesh->elem_list.size();i++) {
                typename INTER::TMESH::EA *e1=Inter.side[0].mesh->elem_list[i];
                typename INTER::Pvec g1 = center(*e1);
                typename INTER::Pvec n1 = Inter.side[0].neq[range(i*INTER::dim,(i+1)*INTER::dim)];

                for( unsigned j=0;j<candidats.size() ;j++ ) {
                    typename INTER::TMESH::EA *e2=Inter.side[1].mesh->elem_list[candidats[j]];
                    typename INTER::Pvec g2 = center(*e2);
                    if ( (length(vect_prod(g2-g1,n1)) <= 1E-3 )and (dot(g2-g1,n1)  >= -1E-6)) {
                        test=1;
                        //on rentre le jeu physique :)
                        if (Inter.comp=="Contact_jeu_physique" )
                            Inter.param_comp->jeu[range(INTER::dim*i,INTER::dim*(i+1))]=g2-g1;
                        //on ecrit la correspondance entre les elements du cote 0 et cote 1
                        Inter.side[1].ddlcorresp[range(i*INTER::dim,(i+1)*INTER::dim)]=range(candidats[j]*INTER::dim,(candidats[j]+1)*INTER::dim);
                        candidats.erase_elem_nb(j);
                        break;
                    } else if (Inter.comp=="periodique" and (length(vect_prod(g2-g1,n1)) <= 1E-6 ) and (dot(g2-g1,n1)  <= 1E-6)) {
                        Inter.side[1].ddlcorresp[range(i*INTER::dim,(i+1)*INTER::dim)]=range(candidats[j]*INTER::dim,(candidats[j]+1)*INTER::dim);
                        candidats.erase_elem_nb(j);
                        break;
                    }
//                    cout << "Num elem correspondant " << i << " " << candidats[j] << " " << length(g2-g1)<< " "  << dot(g2-g1,n1) << " " << length(vect_prod(g2-g1,n1)) << " " << g1 << " " << g2 << " " << n1 << endl;
                }
            }
            //if (test==0) cout << "On n'a jamais atteind le critere qui va bien" << endl;
            if (Inter.comp=="Contact_jeu_physique" and norm_2(Inter.side[0].ddlcorresp-Inter.side[1].ddlcorresp)>1e-6) cout << "Numero : " << Inter.num << "\n" << Inter.side[0].ddlcorresp << "\n" << Inter.side[1].ddlcorresp << endl;
            if (Inter.comp=="periodique" and norm_2(Inter.side[0].ddlcorresp-Inter.side[1].ddlcorresp)>1e-6) cout << "Numero : " << Inter.num << "\n" << Inter.side[0].ddlcorresp << "\n" << Inter.side[1].ddlcorresp << endl;
            if (candidats.size() != 0) {
              cout << "attention maillages incompatibles - non implemente - jeu inter : " << Inter.num << " voisin " << Inter.vois << endl;
                assert(0);
            }

            //cout << "Correspondance terminee" << endl;
        }
    }
};

template <class TV1>
unsigned nb_moment_inertie_inf_eps(TV1 V,double eps){
    unsigned zob=0;
    for(unsigned i=0;i<V.size();i++)
        if (abs(V[i])< eps)
            zob++;
    return zob;
}
//*************************************************
// Application du nombre de fonctions de base par interface
//*************************************************
/** \ingroup Operateurs_inter
\brief Application du nombre de fonctions de base macro pour une interface 
 
Selon le champ  process.multiscale->type_base_macro diff�rentes bases macro sont utilisees. Ce champ peut prendre la valeur 1 (resultante) et le nombre de fcts de base est alors de 2 en 2d et 3 en 3d. Lorsque le champ vaut 2 (resultantes et moments), on a 3 fcts de base en 2d et 6 en 3d. Quand le champ vaut 4 (resultante moment et partie lineaire), plusieurs cas sont consid�r�s selon la forme de l'interface. En 2d, si l'interface est droite on a au total 4 fcts de base, si celle ci est de forme quelconque on utilise 6 fcts. Pour le 3d, pour des interfaces planes on utilise 9 fcts de base et 12 pour des interfaces quelconques.
*/
struct Apply_nb_macro {
    static const double eps=1e-6;
    //cas 3d
    template<class T>
    void operator()(Interface<3,T> &Inter, Param &process) const {

        if(process.multiscale->type_base_macro==1) {
            //resultantes :
            Inter.nb_macro_espace=3;
        } else if(process.multiscale->type_base_macro==2) {
            //resultantes + moments:
            Inter.nb_macro_espace=6;
        } else if(process.multiscale->type_base_macro==3) {
           //resultantes + moments + extensions:
            if (Inter.side[0].mesh->elem_list.size()==0 or nb_moment_inertie_inf_eps(Inter.Moments_inertie,eps)>1) {
                Inter.nb_macro_espace=0;
            } else if(min(abs(Inter.Moments_inertie))<=eps) {
                Inter.nb_macro_espace=9;
            } else {
                Inter.nb_macro_espace=12;
            }
        } else {
            cout << "Nbre de fct non implemente : a modifier" << endl;
            assert(0);
        }
    }
    template<class T>
    void operator()(Interface<2,T> &Inter, Param &process) const {
        if(process.multiscale->type_base_macro==1) {
            //resultantes :
            Inter.nb_macro_espace=2;
        } else if(process.multiscale->type_base_macro==2) {
            //resultantes + moment :
            Inter.nb_macro_espace=3;
        } else if(process.multiscale->type_base_macro==3) {
            //resultantes + moment + extension :
            if (Inter.side[0].mesh->elem_list.size()==0) {
                Inter.nb_macro_espace=0;
            } else if(min(abs(Inter.Moments_inertie))<=eps) {
                Inter.nb_macro_espace=4;
            } else {
                Inter.nb_macro_espace=6;
            }
        } else {
            cout << "Nbre de fct non implemente : a modifier" << endl;
            assert(0);
        }
    }


    /*    template<class INTER >
         void operator()(INTER &Inter, Param &process) const {
            if(process.multiscale->nbmacro_identique==1)
                Inter.nb_macro_espace = process.multiscale->nbmacro;
            else {
                Vec<unsigned> ind;
                ind = find_with_index(process.multiscale->inter_correspondantes==(unsigned)Inter.num);
                if(ind.size()>0) {
                    Inter.nb_macro_espace=process.multiscale->nb_fcts_macro_par_inter[ind[0]];
                } else {
                    Inter.nb_macro_espace=process.multiscale->nbmacro;
                }
            }
        }*/
};

//*************************************************
// calcul BPI (base unique pour l'interface)
//*************************************************
/** \ingroup Operateurs_inter
\brief Cr�ation de la base principale d'inertie de l'interface Interface::BPI
 
Connaissant le centre de gravit� et la position des noeuds �quivalents de l'interface, on d�termine le vecteur position GM, puis on construit l'op�rateur 2x2 ou 3x3 : \f$ \int_{Gamma} \vec{GM} \wedge \vec{GM} dS = GM(i)*M*GM(j) \f$. La recherche des directions principale de cette matrice (class�es par ordre d�croissant des valeurs propres) permet d'obtenir la base principale d'inertie orthonormalis�e ainsi que les moments principaux (utilis�s pour savoir si l'interface est lin�ique ou plane ou non).
*/
struct CalcBPI {
    template<class INTER >
    void operator()(INTER &inter) const {
        typedef Vec<typename INTER::T,INTER::dim> TV;
        typedef typename INTER::T T;
        Vec<TV > nodeeq;
        nodeeq=inter.side[0].nodeeq;
        //calcul de la distance GM (entre G et un noeud de l'interface)
        Mat<T, Gen<>, Dense<> > GM;
        GM.resize(nodeeq.size(),INTER::dim);
        for(unsigned m=0 ;m<GM.nb_rows() ;++m)
            GM.row(m)=nodeeq[m]-inter.G;
        //calcul matrice de masse
        Mat<T, Gen<>, SparseLine<> > M1;
        M1.resize(nodeeq.size(),nodeeq.size());
        M1=inter.side[0].M(range(0,(int)INTER::dim,(int)inter.side[0].M.nb_rows()),range(0,(int)INTER::dim,(int)inter.side[0].M.nb_rows()));
        //calcul matrice d'inertie secondaire int(X'X);
        Mat<T, Gen<INTER::dim>, Dense<> > I1;
        for(unsigned i=0;i<INTER::dim;++i)
            for(unsigned j=0;j<INTER::dim;++j)
                I1(i,j)=dot(GM.col(i),M1*GM.col(j));
        //calcul valeurs et vecteurs propres
        Vec<Vec<T> > V;
        Vec<T> D;
        eig_jacobi(I1,V,D);

        // classement par ordre decroissant des valeurs propres
        Vec<unsigned,INTER::dim> ind;
        D*=-1.0;
        ind=sort_with_index(D);
        inter.BPI=V[ind];
        orthonormalisation_schmidt(inter.BPI);
        inter.Moments_inertie=-1.0*D;
    }
};

//pour affichage
template<unsigned dim>
struct assign_depl {
    template<class TE, class TV>
    void operator()(TE &e, unsigned &compt, TV &W) const {
        Vec<unsigned> rep=range(compt*dim,(compt+1)*dim);
        e.W=W[rep];
        compt+=1;
    }
};


//*************************************************
// creation des projecteurs macro et micro
//*************************************************
/** \ingroup Operateurs_inter
\brief Cr�ation des projecteurs macro et micro (Interface::Side::PM, Interface::Side::Pm, Interface::Side::eM)
 
On d�termine � nouveau le vecteur position GM et en utilisant la base principale d'inertie, on extrait les diff�rentes composantes macro normalis�es. La matrice Interface::Side::eM permet de d�terminer les valeurs de la fonction sur chaque noeud �quivalent � partir des grandeurs macro. La matrice Interface::Side::MeM permet d'extraire les composantes macro d'une distribution donn�e. Enfin l'op�rateur Interface::Side::PM permet d'extraire les valeurs macro aux noeuds �quivalentes d'une distribution donn�e.
 
Par d�faut, 4 fonctions de base macro sont construites en 2D (r�sultantes, moment extension) et 9 en 3D. Proc�dure � modifier pour prendre en compte d'autre fonctions de base macro (+ Repere_ddl_Inter() pour le r�p�rage des ddls macro)
*/
struct CreateProjMacro {
    static const double eps=1e-6;
    // dimension 2
    template<class TT >
    void operator()(Interface<2,TT> &inter, Param &process) const {

        for(unsigned q=0;q<inter.side.size();++q) {
            //if(q==0) {
                unsigned n=inter.side[q].nodeeq.size();
                //calcul de la distance GM (entre G et un noeud de l'interface)
                Vec<Vec<TT,2> > GM;
                GM.resize(n);
                for(unsigned m=0 ;m<n ;++m)
                    GM[m]=inter.side[q].nodeeq[m]-inter.G;

                //fcts de base
                Vec<Vec<TT> > T;
                T.resize(inter.nb_macro_espace);
                for(unsigned i=0;i<(unsigned)inter.nb_macro_espace;++i)
                    T[i].resize(n*2);

                if(inter.nb_macro_espace==2)
                    for(unsigned i=0;i<n;++i) {
                        // resultante selon N1
                        T[0][range(i*2,(i+1)*2)]=inter.BPI[0];
                        // resultante selon N2
                        T[1][range(i*2,(i+1)*2)]=inter.BPI[1];
                    }
                else if(inter.nb_macro_espace==3)
                    for(unsigned i=0;i<n;++i) {
                        // resultante selon N1
                        T[0][range(i*2,(i+1)*2)]=inter.BPI[0];
                        // resultante selon N2
                        T[1][range(i*2,(i+1)*2)]=inter.BPI[1];
                        // Moment selon N3
                        T[2][range(i*2,(i+1)*2)]=dot(GM[i],inter.BPI[0])*inter.BPI[1]-dot(GM[i],inter.BPI[1])*inter.BPI[0];
                    }
                else if (inter.nb_macro_espace==4)
                    for(unsigned i=0;i<n;++i) {
                        // resultante selon N1
                        T[0][range(i*2,(i+1)*2)]=inter.BPI[0];
                        // resultante selon N2
                        T[1][range(i*2,(i+1)*2)]=inter.BPI[1];
                        // Moment selon N3
                        T[2][range(i*2,(i+1)*2)]=dot(GM[i],inter.BPI[0])*inter.BPI[1]-dot(GM[i],inter.BPI[1])*inter.BPI[0];
                        // Extension selon N1
                        T[3][range(i*2,(i+1)*2)]=dot(inter.BPI[0],GM[i])*inter.BPI[0];
                    }
                else if (inter.nb_macro_espace==6) {
                    for(unsigned i=0;i<n;++i) {
                        // resultante selon N1
                        T[0][range(i*2,(i+1)*2)]=inter.BPI[0];
                        // resultante selon N2
                        T[1][range(i*2,(i+1)*2)]=inter.BPI[1];
                        // Moment selon N3
                        T[2][range(i*2,(i+1)*2)]=dot(GM[i],inter.BPI[0])*inter.BPI[1]-dot(GM[i],inter.BPI[1])*inter.BPI[0];
                        // Extension selon N1
                        T[3][range(i*2,(i+1)*2)]=dot(inter.BPI[0],GM[i])*inter.BPI[0];
                        // Extension selon N2
                        T[4][range(i*2,(i+1)*2)]=dot(inter.BPI[1],GM[i])*inter.BPI[1];
                        // cisaillement selon N1,N2
                        T[5][range(i*2,(i+1)*2)]=(dot(inter.BPI[0],GM[i])*inter.BPI[1]+dot(inter.BPI[1],GM[i])*inter.BPI[0])/2.;
                    }
                }

                //orthonormalisation des vecteurs de base
                for(unsigned i=0;i<(unsigned)inter.nb_macro_espace;i++) {
                    Vec<TT> Tsupp;
                    Tsupp=T[i];
                    for(unsigned j=0;j<i;j++)
                        Tsupp-=dot(T[j],inter.side[q].M*T[i])*T[j];
                    T[i] = Tsupp/sqrt(dot( Tsupp,inter.side[q].M*Tsupp ) );
                }

                // stockage matrice de base macro
                inter.side[q].eM.resize(2*n,inter.nb_macro_espace);
                for(unsigned i=0;i<(unsigned)inter.nb_macro_espace;++i)
                    inter.side[q].eM.col(i)=T[i];

                inter.side[q].MeM=inter.side[q].M*inter.side[q].eM;
/*            } else {
                inter.side[q].eM=inter.side[0].eM;
                inter.side[q].MeM=inter.side[0].MeM;
            }*/
        }
    }

    // dimension 3
    template<class TT >
    void operator()(Interface<3,TT> &inter, Param &process) const {
        for(unsigned q=0;q<inter.side.size();++q) {
            //if(q==0) {
                unsigned n=inter.side[q].nodeeq.size();

                //calcul de la distance GM (entre G et un noeud de l'interface)
                Vec<Vec<TT,3> > GM;
                GM.resize(n);
                for(unsigned m=0 ;m<n ;++m)
                    GM[m]=inter.side[q].nodeeq[m]-inter.G;

                //fcts de base
                Vec<Vec<TT> > T;
                T.resize(inter.nb_macro_espace);
                for(unsigned i=0;i<(unsigned)inter.nb_macro_espace;++i)
                    T[i].resize(n*3);

                if(inter.nb_macro_espace==3)
                    for(unsigned i=0;i<n;++i) {
                        // resultante selon N1
                        T[0][range(i*3,(i+1)*3)]=inter.BPI[0];
                        // resultante selon N2
                        T[1][range(i*3,(i+1)*3)]=inter.BPI[1];
                        // resultante selon N3
                        T[2][range(i*3,(i+1)*3)]=inter.BPI[2];
                    }
                else if(inter.nb_macro_espace==6)
                    for(unsigned i=0;i<n;++i) {
                        // resultante selon N1
                        T[0][range(i*3,(i+1)*3)]=inter.BPI[0];
                        // resultante selon N2
                        T[1][range(i*3,(i+1)*3)]=inter.BPI[1];
                        // resultante selon N3
                        T[2][range(i*3,(i+1)*3)]=inter.BPI[2];
                        // Moment selon N1
                        T[3][range(i*3,(i+1)*3)]=vect_prod(inter.BPI[0],GM[i]);
                        // Moment selon N2
                        T[4][range(i*3,(i+1)*3)]=vect_prod(inter.BPI[1],GM[i]);
                        // Moment selon N3
                        T[5][range(i*3,(i+1)*3)]=vect_prod(inter.BPI[2],GM[i]);
                    }
                else if(inter.nb_macro_espace==9)
                    for(unsigned i=0;i<n;++i) {
                        // resultante selon N1
                        T[0][range(i*3,(i+1)*3)]=inter.BPI[0];
                        // resultante selon N2
                        T[1][range(i*3,(i+1)*3)]=inter.BPI[1];
                        // resultante selon N3
                        T[2][range(i*3,(i+1)*3)]=inter.BPI[2];
                        // Moment selon N1
                        T[3][range(i*3,(i+1)*3)]=vect_prod(inter.BPI[0],GM[i]);
                        // Moment selon N2
                        T[4][range(i*3,(i+1)*3)]=vect_prod(inter.BPI[1],GM[i]);
                        // Moment selon N3
                        T[5][range(i*3,(i+1)*3)]=vect_prod(inter.BPI[2],GM[i]);
                        // Extension selon N1
                        T[6][range(i*3,(i+1)*3)]=dot(inter.BPI[0],GM[i])*inter.BPI[0];
                        // Extension selon N2
                        T[7][range(i*3,(i+1)*3)]=dot(inter.BPI[1],GM[i])*inter.BPI[1];
                        // Cisaillement selon N1, N2
                        T[8][range(i*3,(i+1)*3)]=(inter.BPI[0]*dot(inter.BPI[1],GM[i])+inter.BPI[1]*dot(inter.BPI[0],GM[i]))/2.;
                    }
                else if(inter.nb_macro_espace==12)
                    for(unsigned i=0;i<n;++i) {
                        // resultante selon N1
                        T[0][range(i*3,(i+1)*3)]=inter.BPI[0];
                        // resultante selon N2
                        T[1][range(i*3,(i+1)*3)]=inter.BPI[1];
                        // resultante selon N3
                        T[2][range(i*3,(i+1)*3)]=inter.BPI[2];
                        // Moment selon N1
                        T[3][range(i*3,(i+1)*3)]=vect_prod(inter.BPI[0],GM[i]);
                        // Moment selon N2
                        T[4][range(i*3,(i+1)*3)]=vect_prod(inter.BPI[1],GM[i]);
                        // Moment selon N3
                        T[5][range(i*3,(i+1)*3)]=vect_prod(inter.BPI[2],GM[i]);
                        // Extension selon N1
                        T[6][range(i*3,(i+1)*3)]=dot(inter.BPI[0],GM[i])*inter.BPI[0];
                        // Extension selon N2
                        T[7][range(i*3,(i+1)*3)]=dot(inter.BPI[1],GM[i])*inter.BPI[1];
                        // Extension selon N3
                        T[8][range(i*3,(i+1)*3)]=dot(inter.BPI[2],GM[i])*inter.BPI[2];
                        //cisaillement selon N1, N2
                        T[9][range(i*3,(i+1)*3)]=(inter.BPI[0]*dot(inter.BPI[1],GM[i])+inter.BPI[1]*dot(inter.BPI[0],GM[i]))/2.;
                        //cisaillement selon N1, N3
                        T[10][range(i*3,(i+1)*3)]=(inter.BPI[0]*dot(inter.BPI[2],GM[i])+inter.BPI[2]*dot(inter.BPI[0],GM[i]))/2.;
                        //cisaillement selon N2, N3
                        T[11][range(i*3,(i+1)*3)]=(inter.BPI[1]*dot(inter.BPI[2],GM[i])+inter.BPI[2]*dot(inter.BPI[1],GM[i]))/2.;
                    }

                // orthonormalisation des fonctions supplementaires
                for(unsigned i=0;i<(unsigned)inter.nb_macro_espace;i++) {
                    Vec<TT> Tsupp;
                    Tsupp=T[i];
                    for(unsigned j=0;j<i;j++)
                        Tsupp-=dot(T[j],inter.side[q].M*T[i])*T[j];
                    T[i] = Tsupp/sqrt(dot( Tsupp,inter.side[q].M*Tsupp ) );
                }

                // stockage matrice de base macro
                inter.side[q].eM.resize(3*n,inter.nb_macro_espace);
                for(unsigned i=0;i<(unsigned)inter.nb_macro_espace;++i)
                    inter.side[q].eM.col(i)=T[i];
                inter.side[q].MeM=inter.side[q].M*inter.side[q].eM;

//                         DisplayParaview dp;
//                         for(unsigned i=0;i<(unsigned)inter.nb_macro_espace;i++){
//                         unsigned numelem=0;
//                         apply(inter.side[0].mesh.elem_list,assign_depl<3>(),numelem,inter.side[0].eM.col(i));
//                         dp.add_mesh(inter.side[0].mesh);
//                         }
//                         dp.exec();

/*            }
            else {
                inter.side[q].eM=inter.side[0].eM;
                inter.side[q].MeM=inter.side[0].MeM;
            }*/
        }
    }

};



