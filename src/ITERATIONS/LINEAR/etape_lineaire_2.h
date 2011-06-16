#include "containers/matcholamd.h"

//Fonctions utilisees a l'etape micro 2
/** \ingroup etape_lineaire
\relates semilinstage2
 \brief Résolution des problemes micro 2 
 
 Cette étape permet de déterminer les efforts et vitesses (ou déplacements) pour une sous-structure à partir des quantités \f$ \tilde{W} \f$ connues pour un pas de temps donné.
 On effectue les opérations suivantes :
 - Assemblage du second membre:
   - assemblage des quantités \f$ \tilde{W} \f$  des interfaces entourant la sous-structure : \f$ Q_d = k*\dot{\tilde{W}} + \f$
 - Résolution du problème pour obtenir la vitesse (ou déplacement) \f$ q_2 \f$: \f$ \mathbf{K} u + k\frac{W}{\Delta t}  =  Q_d  \f$ en quasistatique et \f$ \mathbf{K} u + k W =  Q_d \f$ en statique
 - Reconstruction du déplacement \f$ W_2 \f$ sur les interfaces à partir de \f$ q_2 \f$
 - Reconstruction de l'effort \f$ F_2 \f$ à partir de la direction de recherche : \f$ F_2 = Q_d - k*\frac{W_2}{\Delta t} \f$ en quasistatique ou \f$ F_2 = Q_d - k*W_2 \f$ en statique
 */
template<class SST,class TV2> void semilin2slave(SST &S,TV2 &Inter, Param &process) {

   unsigned pt=process.temps->pt;

   // assemblage du second membre : droitm
   Vec<typename SST::T> droitm,Qd,W2 ;
   droitm.resize(SST::dim*S.mesh.node_list_size);
   // mise a  zero du second membre
   droitm.set(0.0);

   //assemblage des quantites chapeau des interfaces
   for (unsigned j=0;j<S.edge.size();j++) {
      unsigned data=S.edge[j].datanum;
      unsigned q=S.edge[j].internum;
      Qd=Inter[q].side[data].kglo * Inter[q].side[data].t[pt].WtildeM;
      droitm[S.edge[j].repddledge] +=  Inter[q].side[data].Nt * Inter[q].side[data].M * Qd ;
   }

// TicToc2 tic1;
// tic1.start();
   // resolution du probleme : choix lors de la compilation
   Vec<typename SST::T> Sq;
#if LDL
   Sq=droitm;
   S.l.solve(Sq);
#else
   Sq=S.K->solve(droitm);
#endif
// cout << process.rank<<" : resolution systeme :"; tic1.stop();
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
   };


}

/** \ingroup etape_lineaire
 \brief Structure pour la parallelisation des etapes micro 2 :  Fait appel à la fonction  semilin2slave()
 */
struct semilinstage2 {
   template<class SST, class TV2> void operator()(SST &S,TV2 &Inter,Param &process) const {
      semilin2slave(S,Inter,process);
   }
};


