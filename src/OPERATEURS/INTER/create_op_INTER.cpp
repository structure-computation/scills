#include "create_op_INTER.h"
#include "../../DEFINITIONS/MultiScaleData.h"

#include "../../../LMT/include/containers/vec.h"
#include "../../../LMT/include/containers/mat.h"
#include "../../../LMT/include/mesh/mesh.h"
#include "../../../LMT/include/mesh/remove_doubles.h"
#include "../../../LMT/include/containers/evaluate_nb_cycles.h"

using namespace LMT;

//***********************************************************
// Calcul des matrices de masse M et sous-integration N
//***********************************************************
#include "../../UTILITAIRES/algebre.h"
#include "createMN.h"


/** \ingroup Operateurs_inter
\brief Construction de la matrice de masse et de sous-intégration (Interface::Side::M et Interface::Side::N)
 
On construit la matrice de masse Interface::Side::M par la fonction calcM() puis la matrice de sous-intégration par calcN(). 
- La matrice de masse permet l'intégration de F*W sur une interface et est définie de façon élémentaire par : \f$ M_{ij} = \int_{\Gamma} \psi_i \psi_j dS \f$ où \f$ \psi_i \f$ est la fonction de forme associée aux éléments d'interface. Par défaut, ces fonction de forme sont contantes par élément (P0).
 
- La matrice de sous-intégration N permet la projection de W (interface) à U restriction du déplacement sur le bord d'une sous-structure. Par défaut, l'interpolation de la sous-structure est de degré P1. La matrice élémentaire \f$ N_e \f$ est donc définie par : \f$ N_{e_{ij}}=\int_{\Gamma} \psi_i \phi_j \f$ où \f$ \phi_j \f$ est la fonction de forme associée au noeud i de l'élément de bord de la sous-structure. 
De plus pour éviter des oscillations parasites (notamment observables sur des interfaces de contact ou entre des plis orthotropes...), on assemble les matrices élémentaires \f$ N_e \f$ pour un ensemble d'éléments correspondant à un seul élément d'interface.
 
Une autre manière serait d'avoir un degré P2 sur le bord de la sous-structure. La sous-intégration permettrait alors directement d'éviter les oscillations. (en cours de développement). 
 
Enfin la véritable matrice N est telle que W = N U, on multiplie donc N par \f$ M^{-1} \f$. On calcule en prime N transposé pour accélerer les calculs par la suite.
 
*/
struct CalcMN {
    void operator()(Interface &Inter, Vec<VecPointedValues<Sst> > &S) const {
        for(unsigned j=0;j<Inter.side.size();++j) {
            int i_sst = Inter.side[j].vois[0];
            int i_edge = Inter.side[j].vois[1];
            calcM(*Inter.side[j].mesh,Inter.side[j].M);
            calcN(*Inter.side[j].mesh,*S[i_sst].edge[i_edge].mesh , Inter.side[j].N);
            calcMinvN(Inter.side[j].M,Inter.side[j].N);
            Inter.side[j].Nt=trans(Inter.side[j].N);
            //testMN(Inter.side[j].N,S[i_sst].edge[i_edge].mesh,Inter.side[j].nodeeq,Inter);
        }
    }
};

void testMN(SparseMatrix &N,EdgeMesh &mesh,Vec<Point> &nodeeq, Interface &inter) {

    Vector xbord;
    xbord.resize(N.nb_cols(),0.0);
    for(unsigned i = 0; i < mesh.node_list.size(); ++i) {
        xbord[range(i*DIM,(i+1)*DIM)] = mesh.node_list[i].pos;
        std::cout << mesh.node_list[i].pos << std::endl;
    }
    std::cout << "\nb_nodes" << std::endl;
    std::cout << "xbord" << xbord <<std::endl;

    Vector xinter;
    xinter.resize(N.nb_rows(),0.0);
    for(unsigned i = 0; i < nodeeq.size(); ++i) {
        xinter[range(i*DIM,(i+1)*DIM)] = nodeeq[i];
    }
    //
    std::cout << "xinter " << xinter << std::endl;
    std::cout <<  "Nxbord-xinter" << N*xbord-xinter << std::endl;


};

//********************************************
// Correspondance des ddls d'interface de chaque cote
//*******************************************
/** \ingroup Operateurs_inter
\brief  Creation des listes de correspondances des ddls d'interface - compatibilite. (Interface::Side::ddlcorresp)
 
Ne sert pas réellement même si les correspondances sont utilisées dans l'étape locale, car la manière de construire les interfaces évite un ordre différent pour les noeuds des interfaces (cf. \ref Maillage_geometrie  ).
 
Cette fonction est par contre nécessairement utile pour des maillages d'interface disjoints.  
*/
struct  Corresp_ddlinter {
    void operator()(Interface &Inter) const {
        //std::cout << Inter.vois << std::endl;
        if (Inter.side.size() == 2) {
            if (Inter.side[0].nodeeq.size() != Inter.side[1].nodeeq.size()) {
                std::cout << "attention maillages incompatibles - non implemente (nbnodeeq) - inter : " << Inter.num << " voisin " << Inter.vois << std::endl;
                assert(0);
            }

            Inter.side[0].ddlcorresp.resize(Inter.side[0].mesh->elem_list.size()*DIM);
            Inter.side[1].ddlcorresp.resize(Inter.side[1].mesh->elem_list.size()*DIM);
            Inter.side[0].ddlcorresp=range(0,int(Inter.side[1].nodeeq.size()*DIM));
            Inter.side[1].ddlcorresp.set(0);

            int test=0;
            Vec<unsigned> candidats=range(Inter.side[0].mesh->elem_list.size());
            for (unsigned i=0;i<Inter.side[0].mesh->elem_list.size();i++) {
                InterfaceMesh::EA *e1=Inter.side[0].mesh->elem_list[i];
                Point g1 = center(*e1);
                Point n1 = Inter.side[0].neq[range(i*DIM,(i+1)*DIM)];

                for( unsigned j=0;j<candidats.size() ;j++ ) {
                    InterfaceMesh::EA *e2=Inter.side[1].mesh->elem_list[candidats[j]];
                    Point g2 = center(*e2);
                    if ( (length(vect_prod(g2-g1,n1)) <= 1E-3 )and (dot(g2-g1,n1)  >= -1E-6)) {
                        test=1;
                        /// on rentre le jeu physique :)
                        if (Inter.comp == "Contact_jeu_physique" ){
                            Inter.jeu[range(DIM*i,DIM*(i+1))]=g2-g1;
                        }
                        /// on ecrit la correspondance entre les elements du cote 0 et cote 1
                        Inter.side[1].ddlcorresp[range(i*DIM,(i+1)*DIM)]=range(candidats[j]*DIM,(candidats[j]+1)*DIM);
                        candidats.erase_elem_nb(j);
                        break;
                    } else if (Inter.comp == "periodique" and (length(vect_prod(g2-g1,n1)) <= 1E-6 ) and (dot(g2-g1,n1)  <= 1E-6)) {
                        Inter.side[1].ddlcorresp[range(i*DIM,(i+1)*DIM)]=range(candidats[j]*DIM,(candidats[j]+1)*DIM);
                        candidats.erase_elem_nb(j);
                        break;
                    }
                    //std::cout << "Num elem correspondant " << i << " " << candidats[j] << " " << length(g2-g1)<< " "  << dot(g2-g1,n1) << " " << length(vect_prod(g2-g1,n1)) << " " << g1 << " " << g2 << " " << n1 << std::endl;
                }
            }
            //if (test==0) std::cout << "On nb_nodes'a jamais atteind le critere qui va bien" << std::endl;
            if (Inter.comp=="Contact_jeu_physique" and norm_2(Inter.side[0].ddlcorresp-Inter.side[1].ddlcorresp)>1e-6){
                std::cout << "Numero : " << Inter.num << "\nb_nodes" << Inter.side[0].ddlcorresp << "\nb_nodes" << Inter.side[1].ddlcorresp << std::endl;
            }
            if (Inter.comp=="periodique" and norm_2(Inter.side[0].ddlcorresp-Inter.side[1].ddlcorresp)>1e-6){
                std::cout << "Numero : " << Inter.num << "\nb_nodes" << Inter.side[0].ddlcorresp << "\nb_nodes" << Inter.side[1].ddlcorresp << std::endl;
            }
            if (candidats.size() != 0) {
                std::cout << "attention maillages incompatibles - non implemente - jeu inter : " << Inter.num << " voisin " << Inter.vois << std::endl;
                assert(0);
            }
            //std::cout << "Correspondance terminee" << std::endl;
        }
    }
};


unsigned nb_moment_inertie_inf_eps(Point moments,Scalar eps){
    unsigned zob=0;
    for(unsigned i=0;i<moments.size();i++)
        if (std::abs(moments[i])< eps)
            zob++;
    return zob;
}
//*************************************************
// Application du nombre de fonctions de base par interface
//*************************************************
/** \ingroup Operateurs_inter
\brief Application du nombre de fonctions de base macro pour une interface 
 
Selon le champ  process.multiscale->type_base_macro différentes bases macro sont utilisees. Ce champ peut prendre la valeur 1 (resultante) et le nombre de fcts de base est alors de 2 en 2d et 3 en 3d. Lorsque le champ vaut 2 (resultantes et moments), on a 3 fcts de base en 2d et 6 en 3d. Quand le champ vaut 4 (resultante moment et partie lineaire), plusieurs cas sont considérés selon la forme de l'interface. En 2d, si l'interface est droite on a au total 4 fcts de base, si celle ci est de forme quelconque on utilise 6 fcts. Pour le 3d, pour des interfaces planes on utilise 9 fcts de base et 12 pour des interfaces quelconques.
*/
struct Apply_nb_macro {
    static const Scalar eps=1e-6;
    //cas 3d
    void operator()(Interface &Inter, Process &process) const {
#if DIM == 2
        if(process.multiscale->type_base_macro == 1) {
            //resultantes :
            Inter.nb_macro_espace = 2;
        } else if(process.multiscale->type_base_macro == 2) {
            //resultantes + moment :
            Inter.nb_macro_espace = 3;
        } else if(process.multiscale->type_base_macro == 3) {
            //resultantes + moment + extension :
            if (Inter.side[0].mesh->elem_list.size() == 0) {
                Inter.nb_macro_espace = 0;
                std::cerr << "WARNING !! Aucun element dans l'interface " << Inter.id << std::endl;
            } else if(min(abs(Inter.Moments_inertie))<=eps) {
                Inter.nb_macro_espace = 4;
            } else {
                Inter.nb_macro_espace = 6;
            }
        } else {
            std::cout << "Nbre de fct non implemente : a modifier" << std::endl;
            assert(0);
        }
#elif DIM == 3
        if(process.multiscale->type_base_macro == 1) {
            //resultantes :
            Inter.nb_macro_espace = 3;
        } else if(process.multiscale->type_base_macro == 2) {
            //resultantes + moments:
            Inter.nb_macro_espace = 6;
        } else if(process.multiscale->type_base_macro == 3) {
            //resultantes + moments + extensions:
            unsigned nb_null = 0;
            for(unsigned i=0;i<Inter.Moments_inertie.size();i++){
                if (std::abs(Inter.Moments_inertie[i])< eps){
                    nb_null++;
                }
            }
            if (Inter.side[0].mesh->elem_list.size() == 0){
                Inter.nb_macro_espace = 0;
                std::cerr << "WARNING !! Aucun element dans l'interface " << Inter.id << std::endl;
            } /*else if(nb_null>1) {
                Inter.nb_macro_espace=0;
                std::cerr << "WARNING !! Base macro nulle dans l'interface " << Inter.id << std::endl;
            }*/ else if(min(abs(Inter.Moments_inertie))<=eps) {
                Inter.nb_macro_espace = 9;
            } else {
                Inter.nb_macro_espace = 12;
            }
        } else {
            std::cout << "Nbre de fct non implemente : a modifier" << std::endl;
            assert(0);
        }
#endif
    }
};

//*************************************************
// calcul BPI (base unique pour l'interface)
//*************************************************
/** \ingroup Operateurs_inter
\brief Création de la base principale d'inertie de l'interface Interface::BPI
 
Connaissant le centre de gravité et la position des noeuds équivalents de l'interface, on détermine le vecteur position GM, puis on construit l'opérateur 2x2 ou 3x3 : \f$ \int_{Gamma} \vec{GM} \wedge \vec{GM} dS = GM(i)*M*GM(j) \f$. La recherche des directions principale de cette matrice (classées par ordre décroissant des valeurs propres) permet d'obtenir la base principale d'inertie orthonormalisée ainsi que les moments principaux (utilisés pour savoir si l'interface est linéique ou plane ou non).
*/
struct CalcBPI {
    void operator()(Interface &inter) const {
        const Vec<Point > &nodeeq = inter.side[0].nodeeq;
        const unsigned nb_nodes = nodeeq.size();
        
        /// Calcul de la distance GM (entre G et un noeud de l'interface)
        DenseMatrix GM;
        GM.resize(nb_nodes,DIM);
        for(unsigned m = 0; m < nb_nodes; ++m){
            GM.row(m) = nodeeq[m] - inter.G;
        }
        
        /// Calcul matrice de masse
        SparseMatrix M1;
        M1.resize(nb_nodes,nb_nodes);
        M1=inter.side[0].M(range(0,(int)DIM,(int)inter.side[0].M.nb_rows()),range(0,(int)DIM,(int)inter.side[0].M.nb_rows()));
        
        /// Calcul matrice d'inertie secondaire int(X'X);
        Mat<Scalar, Gen<DIM>, Dense<> > I1;
        for(unsigned i = 0; i < DIM; ++i){
            for(unsigned j = 0; j < DIM; ++j){
                I1(i,j) = dot(GM.col(i),M1*GM.col(j));
            }
        }
        
        /// Calcul valeurs et vecteurs propres
        Vec<Vector > V;
        Vector D;
        eig_jacobi(I1,V,D);

        /// Classement par ordre decroissant des valeurs propres
        Vec<unsigned,DIM> ind;
        D *= -1.0;
        ind = sort_with_index(D);
        inter.BPI = V[ind];
        orthonormalisation_schmidt(inter.BPI);
        inter.Moments_inertie = -1.0*D;
    }
};


//*************************************************
// creation des projecteurs macro et micro
//*************************************************
/** \ingroup Operateurs_inter
\brief Création des projecteurs macro et micro (Interface::Side::PM, Interface::Side::Pm, Interface::Side::eM)
 
On détermine à nouveau le vecteur position GM et en utilisant la base principale d'inertie, on extrait les différentes composantes macro normalisées. La matrice Interface::Side::eM permet de déterminer les valeurs de la fonction sur chaque noeud équivalent à partir des grandeurs macro. La matrice Interface::Side::MeM permet d'extraire les composantes macro d'une distribution donnée. Enfin l'opérateur Interface::Side::PM permet d'extraire les valeurs macro aux noeuds équivalentes d'une distribution donnée.
 
Par défaut, 4 fonctions de base macro sont construites en 2D (résultantes, moment extension) et 9 en 3D. Procédure à modifier pour prendre en compte d'autre fonctions de base macro (+ Repere_ddl_Inter() pour le répérage des ddls macro)
*/
struct CreateProjMacro {
    static const Scalar eps=1e-6;
    void operator()(Interface &inter, Process &process) const {
        for(unsigned q = 0; q < inter.side.size(); ++q) {
            /// Nombre de noeuds sur l'interface
            unsigned nb_nodes = inter.side[q].nodeeq.size();
            /// Calcul de la distance GM (entre G et un noeud de l'interface)
            Vec<Point> GM;
            GM.resize(nb_nodes);
            for(unsigned m = 0; m < nb_nodes; ++m){
                GM[m]=inter.side[q].nodeeq[m]-inter.G;
            }
            /// Initialisation des vecteurs de base
            Vec<Vector> T;
            T.resize(inter.nb_macro_espace);
            for(unsigned i=0;i<(unsigned)inter.nb_macro_espace;++i){
                T[i].resize(nb_nodes*DIM);
            }
            /// Generation des vecteurs
#if DIM == 2
            for(unsigned i = 0; i < nb_nodes; ++i) {
                /// Attention : ordre decroissant et pas de break aux cas 12, 9 et 6 pour que les cas suivants (jusqu'au 3) soit egalement executes
                switch(inter.nb_macro_espace){
                    case 6:
                        // cisaillement selon N1,N2
                        T[5][range(i*2,(i+1)*2)] = (dot(inter.BPI[0],GM[i])*inter.BPI[1]+dot(inter.BPI[1],GM[i])*inter.BPI[0])/2.;
                        // Extension selon N2
                        T[4][range(i*2,(i+1)*2)] = dot(inter.BPI[1],GM[i])*inter.BPI[1];
                    case 4:
                        // Extension selon N1
                        T[3][range(i*2,(i+1)*2)] = dot(inter.BPI[0],GM[i])*inter.BPI[0];
                    case 3:
                        // Moment selon N3
                        T[2][range(i*2,(i+1)*2)] = dot(GM[i],inter.BPI[0])*inter.BPI[1]-dot(GM[i],inter.BPI[1])*inter.BPI[0];
                    case 2:
                        // resultante selon N2
                        T[1][range(i*2,(i+1)*2)] = inter.BPI[1];
                        // resultante selon N1
                        T[0][range(i*2,(i+1)*2)] = inter.BPI[0];
                        break;
                    default:
                        assert(0);  /// Dimension de la base macro incorrecte
                }
            }
#elif DIM == 3
            for(unsigned i = 0; i < nb_nodes; ++i) {
                /// Attention : ordre decroissant et pas de break aux cas 12, 9 et 6 pour que les cas suivants (jusqu'au 3) soit egalement executes
                switch(inter.nb_macro_espace){
                    case 12:
                        //cisaillement selon N2, N3
                        T[11][range(i*3,(i+1)*3)] = (inter.BPI[1]*dot(inter.BPI[2],GM[i])+inter.BPI[2]*dot(inter.BPI[1],GM[i]))/2.;
                        //cisaillement selon N1, N3
                        T[10][range(i*3,(i+1)*3)] = (inter.BPI[0]*dot(inter.BPI[2],GM[i])+inter.BPI[2]*dot(inter.BPI[0],GM[i]))/2.;
                        //cisaillement selon N1, N2
                        T[9][range(i*3,(i+1)*3)] = (inter.BPI[0]*dot(inter.BPI[1],GM[i])+inter.BPI[1]*dot(inter.BPI[0],GM[i]))/2.;
                        // Extension selon N3
                        T[8][range(i*3,(i+1)*3)] = dot(inter.BPI[2],GM[i])*inter.BPI[2];
                    case 9:
                        /// Pour eviter les conflits avec le cas 12
                        if(inter.nb_macro_espace == 9){
                            // Cisaillement selon N1, N2
                            T[8][range(i*3,(i+1)*3)] = (inter.BPI[0]*dot(inter.BPI[1],GM[i])+inter.BPI[1]*dot(inter.BPI[0],GM[i]))/2.;
                        }
                        // Extension selon N2
                        T[7][range(i*3,(i+1)*3)] = dot(inter.BPI[1],GM[i])*inter.BPI[1];
                        // Extension selon N1
                        T[6][range(i*3,(i+1)*3)] = dot(inter.BPI[0],GM[i])*inter.BPI[0];
                    case 6:
                        // Moment selon N3
                        T[5][range(i*3,(i+1)*3)] = vect_prod(inter.BPI[2],GM[i]);
                        // Moment selon N2
                        T[4][range(i*3,(i+1)*3)] = vect_prod(inter.BPI[1],GM[i]);
                        // Moment selon N1
                        T[3][range(i*3,(i+1)*3)] = vect_prod(inter.BPI[0],GM[i]);
                    case 3:
                        // resultante selon N3
                        T[2][range(i*3,(i+1)*3)] = inter.BPI[2];
                        // resultante selon N2
                        T[1][range(i*3,(i+1)*3)] = inter.BPI[1];
                        // resultante selon N1
                        T[0][range(i*3,(i+1)*3)] = inter.BPI[0];
                        break;
                    default:
                        std::cerr << "ERREUR : Interface (id " << inter.id << ") nb_macro_espace = " << inter.nb_macro_espace << std::endl; 
                        assert(0);  /// Dimension de la base macro incorrecte
                }
            }
#endif
            /// Orthonormalisation des vecteurs
            for(unsigned i = 0; i < (unsigned)inter.nb_macro_espace; i++){
                Vector Tsupp;
                Tsupp=T[i];
                for(unsigned j = 0; j < i; j++){
                    Tsupp -= dot(T[j],inter.side[q].M*T[i])*T[j];
                }
                T[i] = Tsupp / std::sqrt(dot( Tsupp,inter.side[q].M*Tsupp));
            }
            
            /// stockage matrice de la base macro
            inter.side[q].eM.resize(DIM*nb_nodes,inter.nb_macro_espace);
            for(unsigned i = 0; i < (unsigned)inter.nb_macro_espace; ++i){
                inter.side[q].eM.col(i) = T[i];
            }
            inter.side[q].MeM=inter.side[q].M*inter.side[q].eM;
        }
    }
};

void create_op_INTER(Vec<VecPointedValues<Sst> > &S,Vec<Interface> &Inter,Vec<VecPointedValues<Interface> > &SubI,Process &process)
{
    /// Calcul des matrices de masse et de sous-integration pour chaque interface
    process.print("\t Matrice de masse et sous-integration ");
    apply_mt(SubI,process.parallelisation->nb_threads,CalcMN(),S);
    
    process.print("\t Correspondance des ddl de chaque cote des interfaces");
    apply_mt(SubI,process.parallelisation->nb_threads,Corresp_ddlinter());
    
    //process.print("\t Centre de gravite des interfaces"); // A REVOIR : inutile ou destruction involontaire
    
    if (process.multiscale->multiechelle == 1){  
        /// Calcul de la Base Principale d'Inertie
        process.print("\t Calcul des BPI");
        apply_mt(SubI,process.parallelisation->nb_threads,CalcBPI());
        /// Application du nombre de fct de base macro par interface 
        apply_mt(SubI,process.parallelisation->nb_threads,Apply_nb_macro(),process);
        /// Creation des projecteurs macro et micro
        process.print("\t Calcul des projecteurs macro ");
        apply_mt(SubI,process.parallelisation->nb_threads,CreateProjMacro(),process);
    }
}
