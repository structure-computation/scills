using namespace LMT;
using namespace std;

#include "correspondance_ddl_sst.h"

//********************************************
// calcul matrice de rigidite par SST
//********************************************
/*\ingroup Operateurs_sst
 \brief Fonction permettant de determiner la matrice de raideur d'une SST. 
  
 Pour modifier le comportement creer des fichiers formulation_... .py (ex orthotropy)
 */
/*struct Calc_SST_rigidite_K0 {
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
    template <class SST, class TV2>
    void operator()(SST &S, TV2 &Inter) const {
        for( unsigned i=0;i<S.edge.size() ;i++ ){
            if (S.edge[i].datanum == 0 or Inter[S.edge[i].internum].comp=="Contact_jeu_physique" or Inter[S.edge[i].internum].comp=="periodique") {
#ifdef PRINT_ALLOC
                if (S.edge[i].mesh != NULL) total_allocated[ typeid(typename SST::TMESHedge).name() ] -= sizeof(typename SST::TMESHedge);
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
#include "util/solveLDL.h"
#include "assignation_materiaux_sst.h"

struct Calc_SST_rigidite_K0_k {
    template<class SST,class TV2>
    void operator()(SST &S, TV2 &Inter,Param &process, DataUser &data_user) const {
        S.f->free_matrices();
        //reperage des ddl de bords (chargement du maillage + non effacement)
        Calc_SST_Correspddl(S,process);
//         S.mesh.load();
        assign_material_on_element(S, data_user);
        S.f->set_mesh(S.mesh.m); 
        S.f->want_amd=false;
        S.f->allocate_matrices();
        S.f->assemble(true,true);
#if LDL
        Mat<TYPEREEL, Sym<>,SparseLine<> > *Kl;
        S.f->get_mat(Kl);
        Mat<TYPEREEL, Sym<>,SparseLine<> > &K = *Kl;
#else
        S.f->get_mat( S.K );
        typename SST::TMATS &K = *S.K;
#endif

#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Verifie memoire apres get_mat : ").c_str(),1);
#endif

//    ofstream f(("K_sst"+to_string(S.num)+"_"+to_string(data_user.options.Multiresolution_current_resolution)).c_str());
//    f << K.nb_cols() << endl;
//    for(unsigned i=0;i<K.data.size();++i) {
//        for(unsigned j=0;j<K.data[i].indices.size();j++){
//            f << (i+1) << " " << (K.data[i].indices[j]+1) << " " << K.data[i].data[j] << "\n" ;
//        }
//    }

        // ajout des directions de recherche sur chaque cote
        for(unsigned j=0;j<S.edge.size();++j) {
            unsigned q=S.edge[j].internum;
            unsigned data=S.edge[j].datanum;            
            K[S.edge[j].repddledge] += Inter[q].side[data].Nt*(Inter[q].side[data].M*(Inter[q].side[data].kglo*Inter[q].side[data].N))/process.temps->dt; // a optimiser
        }
#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Verifie memoire apres ajout ddr : ").c_str(),1);
#endif

//    ofstream f2(("Kk0_"+to_string(S.num)).c_str());
//    f2 << K.nb_cols() << endl;
//    for(unsigned i=0;i<K.data.size();++i) {
//    for(unsigned j=0;j<K.data[i].indices.size();j++){
//      f2 << (i+1) << " " << (K.data[i].indices[j]+1) << " " << K.data[i].data[j] << "\n" ;
//    }
//    }

//    ofstream f(("K_sst"+to_string(S.num)+"_"+to_string(data_user.options.Multiresolution_current_resolution)).c_str());
//    f << K.nb_cols() << endl;
//    for(unsigned i=0;i<K.data.size();++i) {
//        for(unsigned j=0;j<K.data[i].indices.size();j++){
//            f << (i+1) << " " << (K.data[i].indices[j]+1) << " " << K.data[i].data[j] << "\n" ;
//        }
//    }
#if LDL
        S.l.get_factorization( K, true, true );
#else

        K.get_factorization();
#endif

#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Verifie memoire apres factorisation : ").c_str(),1);
#endif
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
    template<class SST,class TV2>
    void operator()(SST &S, TV2 &Inter, Param &process) const {
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

//***********************
// blocage ddl en 2d
//***********************
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
    template<class SST,class TV2>
    void operator()(SST &S, TV2 &Inter, Param &process) const {
        typedef Mat<TYPEREEL, Gen<>, Dense<> > TMAT;

        unsigned nbincmacro=S.nb_macro;
        // initialisation de LE : operateur homogeneise
        TMAT LE;
        LE.resize(nbincmacro);
        // creation des colonnes de LE en appliquant successivement des chargements macro en deplacement (multiplicateur) sur chaque cote
        unsigned repg= 0, repgj= 0;
        
        for(unsigned jj=0;jj<S.edge.size();++jj) {

            unsigned qi=S.edge[jj].internum; //interface selectionnee
            unsigned nbmacro=Inter[qi].nb_macro_espace;
            
            for(unsigned k=0;k<nbmacro;++k) {
                repg=repgj+k;
                //calcul du second membre Qd associe a une deplacement macro d'un cote donnee et assemblage du second membre
                Vec<TYPEREEL> droitm,Wd;
                droitm.resize(DIM*S.mesh.node_list_size);
                droitm.set(0.0);
                for(unsigned j=0;j<S.edge.size();++j) {
                    unsigned q=S.edge[j].internum;
                    unsigned data=S.edge[j].datanum;
                    Wd.resize(Inter[q].side[data].M.nb_cols());
                    if(j==jj) Wd=Inter[q].side[data].eM.col(k);
                    else Wd.set(0);
                    
                    droitm[S.edge[j].repddledge] += Inter[q].side[data].Nt * Inter[q].side[data].M * Inter[q].side[data].kglo * Wd;
                }
                // resolution
                Vec<TYPEREEL> Sq;

#if LDL
                Sq=droitm;
                S.l.solve( Sq );
#else
                Sq = S.K->solve( droitm );
#endif

                // construction des colonnes de LE
                Vec<TYPEREEL> colonneg;

                colonneg.resize(LE.nb_rows());
                for(unsigned j=0;j<S.edge.size();++j) {
                    unsigned q=S.edge[j].internum;
                    unsigned data=S.edge[j].datanum;
                    Wd.resize(Inter[q].side[data].M.nb_cols());
                    if(j==jj) Wd=Inter[q].side[data].eM.col(k);
                    else Wd.set(0);
                    Inter[q].side[data].t[0].F=Inter[q].side[data].kglo*(Wd-Inter[q].side[data].N * Sq[S.edge[j].repddledge]/process.temps->dt);
                    colonneg[S.edge[j].repLE]=trans(Inter[q].side[data].MeM)*Inter[q].side[data].t[0].F;
                }
                LE.col(repg)=colonneg;
            }
            repgj += nbmacro;
        }
        S.LE=LE;

    }
};

