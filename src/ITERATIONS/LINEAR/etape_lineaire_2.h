#include "../../../LMT/include/containers/matcholamd.h"
#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/Interface.h"

//Fonctions utilisees a l'etape micro 2
/** \ingroup etape_lineaire
\relates semilinstage2
 \brief R�solution des problemes micro 2 
 
 Cette �tape permet de d�terminer les efforts et vitesses (ou d�placements) pour une sous-structure � partir des quantit�s \f$ \tilde{W} \f$ connues pour un pas de temps donn�.
 On effectue les op�rations suivantes :
 - Assemblage du second membre:
   - assemblage des quantit�s \f$ \tilde{W} \f$  des interfaces entourant la sous-structure : \f$ Q_d = k*\dot{\tilde{W}} + \f$
 - R�solution du probl�me pour obtenir la vitesse (ou d�placement) \f$ q_2 \f$: \f$ \mathbf{K} u + k\frac{W}{\Delta t}  =  Q_d  \f$ en quasistatique et \f$ \mathbf{K} u + k W =  Q_d \f$ en statique
 - Reconstruction du d�placement \f$ W_2 \f$ sur les interfaces � partir de \f$ q_2 \f$
 - Reconstruction de l'effort \f$ F_2 \f$ � partir de la direction de recherche : \f$ F_2 = Q_d - k*\frac{W_2}{\Delta t} \f$ en quasistatique ou \f$ F_2 = Q_d - k*W_2 \f$ en statique
 */
struct semilinstage2 {
    void operator()(Sst &S,Vec<Interface> &Inter,Process &process) const {
        unsigned pt=process.temps->pt;
        Vec<TYPEREEL> droitm,Qd,W2 ;
        droitm.resize(DIM*S.mesh.node_list_size);
        
        // mise a  zero du second membre
        droitm.set(0.0);
        
        //assemblage des quantites chapeau des interfaces
        for (unsigned j=0;j<S.edge.size();j++) {
            unsigned data=S.edge[j].datanum;
            unsigned q=S.edge[j].internum;
            Qd=Inter[q].side[data].kglo * Inter[q].side[data].t[pt].WtildeM;
            droitm[S.edge[j].repddledge] +=  Inter[q].side[data].Nt * Inter[q].side[data].M * Qd ;
        }
        
        // resolution du probleme : choix lors de la compilation
        Vec<TYPEREEL> Sq;
#if LDL
        Sq=droitm;
        S.l.solve(Sq);
#else
        Sq=S.K->solve(droitm);
#endif
        //sauvegarde des resultats par sst
        if(process.latin->save_depl_SST==1)
            S.t[pt].q+=Sq;
        
        // reconstruction du champ de deplacement sur les cotes de la SST et des Efforts
        for (unsigned j=0; j<S.edge.size(); j++) {
            unsigned q= S.edge[j].internum ;
            unsigned data = S.edge[j].datanum ;
            W2=Inter[q].side[data].N*Sq[S.edge[j].repddledge];
            Inter[q].side[data].t[pt].W+=W2;
            Qd=Inter[q].side[data].kglo * Inter[q].side[data].t[pt].WtildeM;
            Inter[q].side[data].t[pt].F += Qd - Inter[q].side[data].kglo * W2/process.temps->dt;
            //Inter[q].side[data].t[pt].Wp2=(Inter[q].side[data].t[pt].W2-Inter[q].side[data].t[pt-1].W2)/temps.dt;
        }
    }
};

