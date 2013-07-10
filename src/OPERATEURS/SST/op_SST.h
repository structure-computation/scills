#ifndef OP_SST_H
#define OP_SST_H

#include "../../COMPUTE/DataUser.h"

#include "Process.h"
#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/Interface.h"

#include "../../MAILLAGE/correspondance_ddl_sst.h"

#include "../../../LMT/include/containers/vec.h"
#include "../../../LMT/include/containers/mat.h"
#include "../../../LMT/include/util/solveLDL.h"
using namespace LMT;


//********************************************
// calcul matrice de rigidite par SST
//********************************************
/** \ingroup Operateurs_sst
 * \brief Fonction permettant de determiner la matrice de raideur d'une SST. 
 * 
 * Pour modifier le comportement creer des fichiers formulation_... .py (ex orthotropy)
 */
/* A VOIR
struct Calc_SST_rigidite_K0 {
    template<class SST>
    void operator()(SST &S) const {
        S.mesh.load();
        S.f->set_mesh(S.mesh.m); 
        S.f->want_amd=false;
        S.f->allocate_matrices();
        S.f->assemble(true,true);
    }
};
 */



struct efface_mesh_edge{
    void operator()(Sst &S, VecInterfaces &Inter) const {
        for( unsigned i=0;i<S.edge.size() ;i++ ){
            if (S.edge[i].datanum == 0 or Inter[S.edge[i].internum].comp=="Contact_jeu_physique" or Inter[S.edge[i].internum].comp=="periodique") {
#ifdef PRINT_ALLOC
                if (S.edge[i].mesh != NULL) total_allocated[ typeid(EdgeMesh).name() ] -= sizeof(EdgeMesh);
#endif
               delete S.edge[i].mesh;
            }
            S.edge[i].mesh=NULL;
        }
    }
};


//********************************************
// calcul des rigidite globale des SST _ K0 +k
//********************************************
/**\ingroup Operateurs_sst
 *\brief Fonction permettant d'obtenir la matrice de raideur totale (incluant les directions de recherche) sur une SST
 
 Ayant obtenu la matrice élément fini classique correspondant à la formulation choisie, on ajoute la direction de recherche : 
 \f$ \int_{\Omega_E} Tr(\mathbf{K}\epsilon(u) \epsilon(u*) ) d\Omega + \int_{\partial \Omega_E} k W W* dS \f$
 ou en quasistatique :
 \f$ \int_{\Omega_E} Tr(\mathbf{K}\epsilon(\dot{u}) \epsilon(u*) ) d\Omega + \int_{\partial \Omega_E} k \dot{W} W* dS \f$
 
 De manière discrète, on obtient donc pour les termes sur le bord : 
 \f$ K(repddlbord,repddlbord) + \frac{N^t M k N}{\Delta t} \f$
 car \f$ \int_{\partial \Omega_E} k W W* \f$ correspond à \f$ W*^t M k W \f$ et \f$ W = N U[repddlbord] \f$
 */

struct Calc_SST_rigidite_K0_k {
    void operator()(Sst &S, VecInterfaces &Inter,Process &process, DataUser &data_user) const {
        /// Creation de l'operateur de rigidite de la formulation
        S.f->free_matrices();
        S.calc_SST_Correspddl();
        S.apply_behavior();
        process.Fvol->apply_on_sst(S);
        process.Tload->apply_on_sst(S);
        S.f->set_mesh(S.mesh.m);
        S.f->want_amd=false;
        S.f->allocate_matrices();
        S.f->assemble(true,true);
        /// Recuperation de la matrice de raideur
#if LDL
        SymetricMatrix *Kl;
        S.f->get_mat(Kl);
        SymetricMatrix &K = *Kl;
        //std::cout << "Kldl : " << K.nb_rows() << "," << K.nb_cols() << std::endl;
#else
        S.f->get_mat( S.K );
        Sst::CholModMatrix &K = *S.K;
        //std::cout << "Kchol : " << K.nb_rows() << "," << K.nb_cols() << std::endl;
#endif
        
#ifdef PRINT_ALLOC
        disp_alloc((to_string(process.parallelisation->rank)+" : Verifie memoire apres get_mat : ").c_str(),1);
#endif
        /// DEBUG : Recuperation de K
        //ofstream f(("K_sst"+to_string(S.num)+"_"+to_string(data_user.options.Multiresolution_current_resolution)).c_str());
        //f << K.nb_cols() << endl;
        //for(unsigned i=0;i<K.data.size();++i) {
        //    for(unsigned j=0;j<K.data[i].indices.size();j++){
        //        f << (i+1) << " " << (K.data[i].indices[j]+1) << " " << K.data[i].data[j] << "\n" ;
        //    }
        //}

        /// Ajout des directions de recherche sur chaque cote
        for(unsigned j=0;j<S.edge.size();++j) {
            unsigned q=S.edge[j].internum;
            unsigned data=S.edge[j].datanum;
            K[S.edge[j].repddledge] += Inter[q].side[data].Nt*(Inter[q].side[data].M*(Inter[q].side[data].kglo*Inter[q].side[data].N))/process.temps->dt; // a optimiser
        }
        
#ifdef PRINT_ALLOC
        disp_alloc((to_string(process.parallelisation->rank)+" : Verifie memoire apres ajout ddr : ").c_str(),1);
#endif
        /// DEBUG : Recuperation de K + k0
        //ofstream f2(("Kk0_"+to_string(S.num)).c_str());
        //f2 << K.nb_cols() << endl;
        //for(unsigned i=0;i<K.data.size();++i) {
        //    for(unsigned j=0;j<K.data[i].indices.size();j++){
        //        f2 << (i+1) << " " << (K.data[i].indices[j]+1) << " " << K.data[i].data[j] << "\n" ;
        //    }
        //}
        
        /// Factorisation de la matrice
#if LDL
        S.l.get_factorization( K, true, true );
#else
        K.get_factorization();
#endif

#ifdef PRINT_ALLOC
        disp_alloc((to_string(process.parallelisation->rank)+" : Verifie memoire apres factorisation : ").c_str(),1);
#endif
        /// Liberation du maillage LMTpp
        S.mesh.unload();
    }
};

//***************************************
// reperage des inconnues macro dans LE
//***************************************
/**\ingroup Operateurs_sst
 \brief Reperage des inconnues macro dans l'operateur homogeneise pour une SST
 
On détermine le vecteur d'indice permettant de repérer le coté de la sous-structure (ou l'interface) dans l'opérateur homogénéisé selon le nombre de composante macro des interfaces entourant la sous-structure.
 */
struct repere_ind_interface_LE {
    void operator()(Sst &S, Vec<Interface> &Inter, Process &process) const {
        unsigned nbmacroS=0;
        for(unsigned j=0;j<S.edge.size();++j) {
            unsigned q=S.edge[j].internum;
            unsigned nbmacro=Inter[q].nb_macro_espace;
            S.edge[j].repLE=nbmacroS+range(nbmacro);
            nbmacroS+=nbmacro;
        }
        S.nb_macro=nbmacroS;
    }
};


//**************************************************************
// blocage ddl en 3d
// repddl : reperage des ddls dans les ddls globaux du maillage
// ddlb : relation entre les ddls , equations supplementaires
//**************************************************************
/* PAS UTILISE
template<class TM,class TD,class TOP>
void ddlbloq(TM &m,  Mat<typename TM::Tpos, Gen<6,9>, TD,TOP> &ddlb, Vec<unsigned,9> &repddl) {
    typedef typename TM::Tpos T;
    T eps=1e-15;

    // reperage de 3 noeuds non alignes
    Vec<unsigned> rep(0,1,2);
    unsigned ind=2;
    Vec<T,3> V1,V2;
    V1=m.node_list[1].pos-m.node_list[0].pos;
    V2=m.node_list[ind].pos-m.node_list[0].pos;
    while( find(abs(vect_prod(V1,V2)),LMT::_1>=eps)==0 ) {
        ind+=1;
        V2=m.node_list[ind].pos-m.node_list[0].pos;
    }
    rep[2]=ind;
    for(unsigned i=0;i<3;++i) {
        repddl[range(i*3,(i+1)*3)]=range(rep[i]*3,(rep[i]+1)*3);
    }

    // calcul des normales
    Vec<T,3> normale1, normale2;
    normale1=vect_prod(V1,V2);
    normale2=vect_prod(normale1,V1);
    normale2=normale2/sqrt(dot(normale2,normale2));
    normale1=normale1/sqrt(dot(normale1,normale1));

    // initialisation ddlb
    ddlb.set(0.0);
    // blocage noeud1 tous ddls
    for(unsigned i=0;i<3;++i)
        ddlb(i,i)=1.0;
    // blocage noeud2 : ddl.n et ddl.n2
    ddlb.row(3)[range(3,6)]=normale1;
    ddlb.row(4)[range(3,6)]=normale2;
    //blocage noeud3 : ddl.n
    ddlb.row(5)[range(6,9)]=normale1;
};
//*/

//***********************
// blocage ddl en 2d
//***********************
/* PAS UTILISE
template<class TM,class TD,class TOP>
void ddlbloq(TM &m, Mat<typename TM::Tpos,Gen<3,4>,TD,TOP> &ddlb,Vec<unsigned,4> &repddl) {
    typedef typename TM::Tpos T;

    // reperage de 2 noeuds non alignes
    Vec<unsigned> rep(0,1);
    Vec<T,2> V1;
    V1=m.node_list[1].pos-m.node_list[0].pos;
    repddl=range(4);
    // calcul des normales
    Vec<T,2> normale1;
    normale1[0]=-1.0*V1[1];
    normale1[1]=V1[0];
    normale1=normale1/norm_2(normale1);

    // initialisation ddlb
    ddlb.set(0.0);
    // blocage noeud1 tous ddls
    for(unsigned i=0;i<2;++i)
        ddlb(i,i)=1.0;
    // blocage noeud2 : ddl.n
    ddlb.row(2)[range(2,4)]=normale1;
};
//*/

//*************************************************
// calcul operateur homogeneise en temps
//*************************************************
/** \ingroup Operateurs_sst
\brief Fonction permettant de construire l'opérateur homogénéisé LE pour une etude quasistatique
 
Pour chaque sous-structure, on boucle sur le nombre de coté et on détermine l'interface correspondante. On applique alors sucessivement une fonction macro non nulle sur l'interface selectionnee et nulle sur les autres interfaces.
 
On assemble ensuite le second membre correspondant à ce déplacement imposé sur chaque coté :
\f$ \int_{\Omega_E} Tr(\sigma \epsilon(u*)) d\Omega + \int_{\partial \Omega_E} k\dot{W} W* dS= \int_{\partial \Omega_E} k\tilde{W}^M W* dS  \f$
 
On résout donc ce problème pour en extraire les déplacements sur le bord de la sous-structure (par interface). Par la direction de recherche :
\f$ F =  k(\tilde{W}^M - k\dot{W} \f$, on en déduit l'effort sur chaque interface et on en prend la partie macro. On entre enfin chaque composante macro pour chaque interface dans la colonne de l'opérateur homogénéisé correspondant au multiplicateur non nul imposé. 
 
*/
struct Calc_SST_LE {
    void operator()(Sst &S, VecInterfaces &Inter, Process &process) const {
        /// Initialisation de LE : operateur homogeneise
        DenseMatrix LE;
        LE.resize(S.nb_macro);
        /// Creation des colonnes de LE en appliquant successivement des chargements macro en deplacement (multiplicateur Wtilde) sur chaque bord du solide
        unsigned repg= 0, repgj= 0;
        
        for(unsigned jj=0;jj<S.edge.size();++jj) {

            unsigned qi=S.edge[jj].internum;    /// Interface selectionnee
            unsigned nbmacro=Inter[qi].nb_macro_espace; /// Dimension de l'espace macro
            
            for(unsigned k=0;k<nbmacro;++k) {
                repg=repgj+k;
                /// Calcul du second membre Qd associe a une deplacement macro d'un cote donnee et assemblage du second membre
                Vector Qd,Wd;
                Qd.resize(DIM*S.mesh.node_list_size);
                Qd.set(0.0);
                for(unsigned j=0;j<S.edge.size();++j) {
                    unsigned q=S.edge[j].internum;
                    unsigned data=S.edge[j].datanum;
                    Wd.resize(Inter[q].side[data].M.nb_cols());
                    if(j==jj) Wd=Inter[q].side[data].eM.col(k);
                    else Wd.set(0.0);
                    /*
                    std::cout << "********** CREATION LE **********" << std::endl;
                    PRINT(jj);
                    PRINT(k << std::endl;
                    PRINT(j << std::endl;
                    PRINT(q << std::endl;
                    PRINT(data << std::endl;
                    std::cout << "kglo : " << Inter[q].side[data].kglo.nb_rows() << "," << Inter[q].side[data].kglo.nb_cols() << std::endl;
                    std::cout << "Nt : " << Inter[q].side[data].Nt.nb_rows() << "," << Inter[q].side[data].Nt.nb_cols() << std::endl;
                    display(std::cout,Inter[q].side[data].Nt,0);
                    std::cout << "M  : " << Inter[q].side[data].M.nb_rows()  << "," << Inter[q].side[data].M.nb_cols()  << std::endl;
                    std::cout << "eM : " << Inter[q].side[data].eM.nb_rows() << "," << Inter[q].side[data].eM.nb_cols() << std::endl;
                    display(std::cout,Inter[q].side[data].M,0);
                    std::cout << "Wd  : " << Wd.size() << std::endl;
                    for(int toto = 0; toto < Wd.size(); toto++)
                        std::cout << "    " << toto << " : " << Wd[toto];
                    */
                    Qd[S.edge[j].repddledge] += Inter[q].side[data].Nt * Inter[q].side[data].M * Inter[q].side[data].kglo * Wd;
                }
                
                /// Resolution
                Vector Sq;
#if LDL
                Sq=Qd;
                S.l.solve(Sq);
#else
                Sq = S.K->solve(Qd);
#endif

                /// Construction des colonnes de LE
                Vector colonne;

                colonne.resize(LE.nb_rows());
                for(unsigned j=0;j<S.edge.size();++j) {
                    unsigned q=S.edge[j].internum;
                    unsigned data=S.edge[j].datanum;
                    Wd.resize(Inter[q].side[data].M.nb_cols());
                    if(j==jj) Wd=Inter[q].side[data].eM.col(k);
                    else Wd.set(0);
                    Inter[q].side[data].t[0].F=Inter[q].side[data].kglo*(Wd-Inter[q].side[data].N * Sq[S.edge[j].repddledge]/process.temps->dt);
                    colonne[S.edge[j].repLE]=trans(Inter[q].side[data].MeM)*Inter[q].side[data].t[0].F;
                }
                LE.col(repg)=colonne;
            }
            repgj += nbmacro;
        }
        S.LE=LE;
    }
};

#endif //OP_SST_H
