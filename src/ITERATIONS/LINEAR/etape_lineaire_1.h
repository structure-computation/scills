#include "../../DEFINITIONS/Process.h"
#include "../../DEFINITIONS/Interface.h"
#include "../../DEFINITIONS/TimeData.h"
#include "../../DEFINITIONS/LatinData.h"

#include "../../../LMT/include/containers/vec.h"
#include "../../../LMT/include/containers/evaluate_nb_cycles.h"

//Fonctions utilisees a l'etape micro 1

/** \ingroup etape_lineaire
\relates semilinstage1
\brief Resolution des problemes micro 1 
 
 Cette etape permet de determiner les efforts et vitesses pour une sous-structure a partir des quantites connues pour un pas de temps donne.
 On effectue les operations suivantes :
 - Assemblage du second membre:
   - obtention du terme de second membre a  partir de la formulation ( thermique et quantites chapeau par sst)
   - assemblage des quantites chapeau des interfaces entourant la sous-structure : \f$ Q_d = \hat{F} + k*\dot{\hat{W}} + k \frac{W_n}{\Delta t} \f$ en quasistatique ou \f$ Q_d = \hat{F} + k*\hat{W} \f$ en statique
 - Resolution du problème pour obtenir le deplacement \f$ q_1 \f$: \f$ \mathbf{K}\dot{u} + k \frac{W}{\Delta t}  =  Q_d  \f$ en quasistatique ou \f$ \mathbf{K} u + k W =  Q_d \f$ en statique
 - Reconstruction du deplacement \f$ W_1 \f$  sur les interfaces a partir \f$ q_1 \f$
 - Reconstruction de l'effort \f$ F_1 \f$ a partir de la direction de recherche : \f$ F_1 = Q_d - k*\frac{W_1}{\Delta t} \f$ en quasistatique ou en statique
 
 Pour la partie multiechelle, on construit enfin le terme a ajouter au second membre macro \f$ \int_{\partial \Omega_E} \hat{F}_d \tilde{W}^* \f$ , provenant des differents cotes des interfaces. Ici, \f$ \hat{F}_d = F_1 \f$ et on effectue simplement l'integrale et la multiplication par la base macro pour extraire les composantes macro a assembler.
 */
struct semilinstage1 {
    void operator()(Sst &S,VecInterfaces &Inter,Process &process) const {
        unsigned pt=process.temps->pt;
        
        /// assemblage du second membre : droitm
        Vector droitm,Qd ;
        droitm.resize(DIM*S.mesh.node_list_size,0.);
        
        ///second membre prenant en compte le comportement thermique et les quantites chapeau:
        droitm = S.fvol;
        
        ///assemblage des quantites chapeau des interfaces
        for (unsigned j=0;j<S.edge.size();j++) {
            unsigned data=S.edge[j].datanum;
            unsigned q=S.edge[j].internum;
            Qd = Inter[q].side[data].t[pt].Fchap + Inter[q].side[data].kglo * (Inter[q].side[data].t[pt].Wpchap + Inter[q].side[data].t[pt-1].W/process.temps->dt);
//             if(Inter[q].id==8){
//                 PRINT("  ");
//                 PRINT("avant lineaire 1  ");
//                 PRINT(Inter[q].id);
//                 PRINT(data);
//                 PRINT(Inter[q].side[data].t[pt-1].W[LMT::range(0,DIM*1)]);
//                 PRINT(Inter[q].side[data].t[pt].Wpchap[LMT::range(0,DIM*1)]);
//                 PRINT(Inter[q].side[data].t[pt].Fchap[LMT::range(0,DIM*1)]);
//                 PRINT(Qd[LMT::range(0,DIM*1)]);
//                 PRINT("  ");
//             }
            droitm[S.edge[j].repddledge] +=  Inter[q].side[data].Nt * Inter[q].side[data].M * Qd ;
        };
        
        ///calcul du deplacement q1
        Vector Sq;
#if LDL
        Sq=droitm;
        S.l.solve(Sq);
#else
        Sq=S.K->solve(droitm);
#endif
        
        ///sauvegarde des resultats par sst
        if(process.latin->save_depl_SST==1)
            S.t[pt].q=Sq;
        
        /// reconstruction du champ de vitesse sur les cotes de la SST
        for (unsigned j=0; j<S.edge.size(); j++) {
            unsigned q= S.edge[j].internum ;
            unsigned data = S.edge[j].datanum ;
            Inter[q].side[data].t[pt].W=Inter[q].side[data].N*Sq[S.edge[j].repddledge];
            Qd = Inter[q].side[data].t[pt].Fchap + Inter[q].side[data].kglo * (Inter[q].side[data].t[pt].Wpchap + Inter[q].side[data].t[pt-1].W/process.temps->dt);
            Inter[q].side[data].t[pt].F = Qd - Inter[q].side[data].kglo * Inter[q].side[data].t[pt].W/process.temps->dt;
            //Inter[q].side[data].t[pt].Wp1=(Inter[q].side[data].t[pt].W1-Inter[q].side[data].t[pt-1].W)/process.temps->dt;
            
//             if(Inter[q].id==8){
//                 PRINT("  ");
//                 PRINT("après lineaire 1  ");
//                 PRINT(Inter[q].id);
//                 PRINT(data);
//                 PRINT(Inter[q].side[data].t[pt].W[LMT::range(0,DIM*1)]);
//                 PRINT(Inter[q].side[data].t[pt].F[LMT::range(0,DIM*1)]);
//                 PRINT("  ");
//             }
        }
        
        /// etape non utile en monoechelle
        if (process.multiscale->multiechelle ==1) {
            /// construction du vecteur d'efforts macro donnes pour le second membre du probleme macro par SST
            S.Fadd.resize(S.nb_macro);
            for (unsigned j=0; j <S.edge.size(); j++) {
                unsigned data=S.edge[j].datanum ;
                unsigned q=S.edge[j].internum ;
                S.Fadd[S.edge[j].repLE]=trans(Inter[q].side[data].MeM) * Inter[q].side[data].t[pt].F;
            }
        }
        //std::cout << "etape lineaire, pb micro 1, solide " << S.id << " :" << std::endl;
        //PRINT(S.t[pt].q);
    }
};

