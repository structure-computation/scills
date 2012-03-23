#include "../../../LMT/include/containers/evaluate_nb_cycles.h"

//Fonctions utilisees a l'etape micro 1

/** \ingroup etape_lineaire
\relates semilinstage1
 \brief Résolution des problemes micro 1 
 
 Cette étape permet de déterminer les efforts et vitesses pour une sous-structure à partir des quantités connues pour un pas de temps donné.
 On effectue les opérations suivantes :
 - Assemblage du second membre:
   - obtention du terme de second membre à  partir de la formulation ( thermique et quantites chapeau par sst)
   - assemblage des quantités chapeau des interfaces entourant la sous-structure : \f$ Q_d = \hat{F} + k*\dot{\hat{W}} + k \frac{W_n}{\Delta t} \f$ en quasistatique ou \f$ Q_d = \hat{F} + k*\hat{W} \f$ en statique
 - Résolution du problème pour obtenir le déplacement \f$ q_1 \f$: \f$ \mathbf{K}\dot{u} + k \frac{W}{\Delta t}  =  Q_d  \f$ en quasistatique ou \f$ \mathbf{K} u + k W =  Q_d \f$ en statique
 - Reconstruction du déplacement \f$ W_1 \f$  sur les interfaces à partir \f$ q_1 \f$
 - Reconstruction de l'effort \f$ F_1 \f$ à partir de la direction de recherche : \f$ F_1 = Q_d - k*\frac{W_1}{\Delta t} \f$ en quasistatique ou en statique
 
 Pour la partie multiéchelle, on construit enfin le terme à ajouter au second membre macro \f$ \int_{\partial \Omega_E} \hat{F}_d \tilde{W}^* \f$ , provenant des différents cotés des interfaces. Ici, \f$ \hat{F}_d = F_1 \f$ et on effectue simplement l'intégrale et la multiplication par la base macro pour extraire les composantes macro à assembler.
 */
struct semilinstage1 {
    void operator()(Sst &S,Vec<Interface> &Inter,Param &process) const {
        unsigned pt=process.temps->pt;
        // assemblage du second membre : droitm
        Vec<TYPEREEL> droitm,Qd ;
        droitm.resize(DIM*S.mesh.node_list_size,0.);
        
        //second membre prenant en compte le comportement thermique et les quantites chapeau:
        droitm = S.fvol;
        
        //assemblage des quantites chapeau des interfaces
        for (unsigned j=0;j<S.edge.size();j++) {
            unsigned data=S.edge[j].datanum;
            unsigned q=S.edge[j].internum;
            Qd = Inter[q].side[data].t[pt].Fchap + Inter[q].side[data].kglo * (Inter[q].side[data].t[pt].Wpchap + Inter[q].side[data].t[pt-1].W/process.temps->dt);
            droitm[S.edge[j].repddledge] +=  Inter[q].side[data].Nt * Inter[q].side[data].M * Qd ;
        };
        //calcul du deplacement q1
        Vec<TYPEREEL> Sq;
#if LDL
        Sq=droitm;
        S.l.solve(Sq);
#else
        Sq=S.K->solve(droitm);
#endif
        //sauvegarde des resultats par sst
        if(process.latin->save_depl_SST==1)
            S.t[pt].q=Sq;
        
        // reconstruction du champ de vitesse sur les cotes de la SST
        for (unsigned j=0; j<S.edge.size(); j++) {
            unsigned q= S.edge[j].internum ;
            unsigned data = S.edge[j].datanum ;
            Inter[q].side[data].t[pt].W=Inter[q].side[data].N*Sq[S.edge[j].repddledge];
            Qd = Inter[q].side[data].t[pt].Fchap + Inter[q].side[data].kglo * (Inter[q].side[data].t[pt].Wpchap + Inter[q].side[data].t[pt-1].W/process.temps->dt);
            Inter[q].side[data].t[pt].F = Qd - Inter[q].side[data].kglo * Inter[q].side[data].t[pt].W/process.temps->dt;
            //Inter[q].side[data].t[pt].Wp1=(Inter[q].side[data].t[pt].W1-Inter[q].side[data].t[pt-1].W)/process.temps->dt;
        }
        
        // etape non utile en monoechelle
        if (process.multiscale->multiechelle ==1) {
            // construction du vecteur d'efforts macro donnes pour le second membre du probleme macro par SST
            S.Fadd.resize(S.nb_macro);
            for (unsigned j=0; j <S.edge.size(); j++) {
                unsigned data=S.edge[j].datanum ;
                unsigned q=S.edge[j].internum ;
                S.Fadd[S.edge[j].repLE]=trans(Inter[q].side[data].MeM) * Inter[q].side[data].t[pt].F;
            }
        }
    }
};

