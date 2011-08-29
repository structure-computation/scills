using namespace LMT;
using namespace std;

#include "containers/vec_mt.h"
#include "etape_lineaire_1.h"
#include "etape_lineaire_2.h"
#include "etape_macro.h"

//fcts MPI
#include "mpi_transactions.h"
#include "containers/evaluate_nb_cycles.h"
extern Crout crout;

/** \defgroup etape_lineaire Etape linéaire
\ingroup LATIN
\brief Description des phases de l'étape linéaire
 
 Cette étape de la stratégie LATIN est constituée des phases suivantes :
 - pour chaque intervalle de temps on fait :
   - étape linéaire 1 : \ref semilinstage1
   - étape macro (seulement en multiéchelle)
      - construction du second membre macro :  \ref macroassemble
      - résolution du problème macro: Glob::resolmacro()
      - extraction des valeurs du multiplicateur macro et assignation aux interfaces : \ref interextrmacro
      - résolution du problème micro 2 : \ref semilinstage2
   - reconstruction des quantités à partir des deux étapes :\ref  reconstruction_quantites
 - une fois la boucle sur le temps terminée on effectue une relaxation des quantités :\ref  relaxation_quantites
 
 Rq : Pour la statique un seul intervalle de temps est réalisé. (Le pas de temps n'a alors pas d'importance)
 
*/

/** \ingroup etape_lineaire
\relates reconstruction_quantites
\brief Etape Lineaire : Reconstruction du déplacement et des quantités d'interface
 
En sommant le champ solution \f$ q_1 \f$ de l'étape linéaire 1 et le champ solution \f$ q_2 \f$ de l'étape linéaire 2, on reconstruit le déplacement \f$ q \f$ puis par dérivation on peut obtenir la vitesse pour l'intervalle de temps considéré.
De la même manière on reconstruit les efforts, déplacements et vitesses d'interface. Pour le monoéchelle, on n'utilise que le champ q1
 
Il est possible ici d'extraire les efforts macro ou micro, en utilisant le projecteur adéquat (cf : Interface).
*/
template<class SST, class TV2>
void reconstruction(SST &S,TV2 &Inter,Param &process) {
    unsigned  pt=process.temps->pt;
    
    //reconstruction du deplacement si necessaire
//     if(process.latin->save_depl_SST==1) {
//         if(process.multiscale->multiechelle==1)
//             S.t[pt].q=S.t[pt].q1+S.t[pt].q2;
//         else
//             S.t[pt].q=S.t[pt].q1;
//     }

    if(process.multiscale->multiechelle==1) {
        for(unsigned j=0;j<S.edge.size();++j) {
            unsigned data=S.edge[j].datanum;
            unsigned q=S.edge[j].internum;
            //Inter[q].side[data].t[pt].F=Inter[q].side[data].t[pt].F1+Inter[q].side[data].t[pt].F2;
            //Inter[q].side[data].t[pt].Wp=Inter[q].side[data].t[pt].Wp1 + Inter[q].side[data].t[pt].Wp2;
            //Inter[q].side[data].t[pt].W=Inter[q].side[data].t[pt].W1 + Inter[q].side[data].t[pt].W2;
            Inter[q].side[data].t[pt].Wp=(Inter[q].side[data].t[pt].W-Inter[q].side[data].t[pt-1].W)/process.temps->dt;
            //extraction des quantites micro et macro si besoin
            //Vec<typename INTER::T> FM;
            //FM = Inter[q].side[data].PM(Inter[q].side[data].t[pt].F);
            //Inter[q].side[data].Fm=Inter[q].side[data].Pm*Inter[q].side[data].F;
        }
    } else {
        for(unsigned j=0;j<S.edge.size();++j) {
            unsigned data=S.edge[j].datanum;
            unsigned q=S.edge[j].internum;
            //Inter[q].side[data].t[pt].F=Inter[q].side[data].t[pt].F1;
            //Inter[q].side[data].t[pt].W=Inter[q].side[data].t[pt].W1;
            Inter[q].side[data].t[pt].Wp=(Inter[q].side[data].t[pt].W-Inter[q].side[data].t[pt-1].W)/process.temps->dt;
        }
    }
}
/** \ingroup etape_lineaire
    \relates derivation
    \brief Etape Lineaire : Reconstruction de la vitesse
 
 */
template<class SST, class TV2>
void derivation(SST &S,TV2 &Inter,Param &process) {
  unsigned  pt=process.temps->pt;

  for(unsigned j=0;j<S.edge.size();++j) {
    unsigned data=S.edge[j].datanum;
    unsigned q=S.edge[j].internum;
    Inter[q].side[data].t[pt].Wp=(Inter[q].side[data].t[pt].W-Inter[q].side[data].t[pt-1].W)/process.temps->dt;
  }
}

/** \ingroup etape_lineaire
\relates relaxation_quantites
\brief Etape Lineaire : Relaxation
 
Ayant déterminé pour chaque intervalle de temps la solution q pour les Sst, la vitesse \f$ \dot{W} ( ou déplacement \f$ W \f$ ) \f$ et effort \f$ F \f$ sur chaque interface, on procède à la relaxation, permettant d'obtenir une meilleure convergence.
 
On effectue donc l'opération suivante, avec mu=0.8 par défaut (1 à la première itération) :
\code q = (1-mu)*oldq + mu*q \endcode
où oldq est le déplacement solution à l'itération précédente. 
On effectue de même avec les quantités d'interfaces, puis on met à jour les anciennes quantités.
*/
template<class SST, class TV2>
void relaxation (SST &S,TV2 &Inter,Param &process) {
    typename SST::T mu = process.latin->mu;
    unsigned pt = process.temps->pt;
//     if(process.latin->save_depl_SST==1) {
//         S.t[pt].q=(1-mu)*S.t[pt].oldq + mu*S.t[pt].q;
//         S.t[pt].oldq=S.t[pt].q;
//     }

    for(unsigned j=0;j<S.edge.size();++j) {
        unsigned q=S.edge[j].internum;
        unsigned data=S.edge[j].datanum;
        Inter[q].side[data].t[pt].F=(1-mu)*Inter[q].side[data].t[pt].oldF + mu*Inter[q].side[data].t[pt].F;
        Inter[q].side[data].t[pt].Wp=(1-mu)*Inter[q].side[data].t[pt].oldWp + mu*Inter[q].side[data].t[pt].Wp;
        Inter[q].side[data].t[pt].W=(1-mu)*Inter[q].side[data].t[pt].oldW + mu*Inter[q].side[data].t[pt].W;

        Inter[q].side[data].t[pt].oldF=Inter[q].side[data].t[pt].F;
        Inter[q].side[data].t[pt].oldWp=Inter[q].side[data].t[pt].Wp;
        Inter[q].side[data].t[pt].oldW=Inter[q].side[data].t[pt].W;

    }
}


/** \ingroup etape_lineaire
\relates integration_quantites
\brief Etape Lineaire : Integration
 
Connaissant la vitesse pour chaque pas de temps, il est possible de déterminer le déplacement à chaque piquet de temps en prenant en compte le déplacement initial (t=0) au départ.
*/
template<class INTER>
void integration (INTER &Inter,Param &process) {
    for(unsigned j=0;j<Inter.side.size();j++) {
        for(unsigned pt=1;pt<=process.temps->nbpastemps;pt++)
            Inter.side[j].t[pt].W = Inter.side[j].t[pt-1].W + process.temps->dt *  Inter.side[j].t[pt].Wp;
    }

}

/** \ingroup etape_lineaire
\relates derivation_quantites
\brief Etape Lineaire : Derivation
*/
template<class INTER>
void derivation (INTER &Inter,Param &process) {
    for(unsigned j=0;j<Inter.side.size();j++) {
        for(unsigned pt=1;pt<=process.temps->nbpastemps;pt++)
            Inter.side[j].t[pt].Wp = (Inter.side[j].t[pt].W - Inter.side[j].t[pt-1].W)/process.temps->dt ;
    }

}

/**\ingroup etape_lineaire
 \brief Etape Lineaire : Reconstruction du déplacement et des quantités d'interface : cf reconstruction()
 */
struct reconstruction_quantites {
    template<class SST, class TV2>
    void operator()(SST &S,TV2 &Inter,Param &process) const {
        reconstruction(S,Inter,process);
    }
};

/**\ingroup etape_lineaire
 \brief Etape Lineaire : Relaxation : cf relaxation()
 */
struct relaxation_quantites {
    template<class SST, class TV2>
    void operator()(SST &S,TV2 &Inter,Param &process) const {
        relaxation (S,Inter,process);
    }
};



/**\ingroup etape_lineaire
 \brief Etape Lineaire : Integration du deplacement sur les interfaces : cf integration_depl()
 */
struct integration_quantites {
    template< class INTER>
    void operator()(INTER &Inter,Param &process) const {
        integration(Inter,process);
    }
};

/**\ingroup etape_lineaire
 \brief Etape Lineaire : Dérivation du deplacement sur les interfaces : cf derivation()
 */
struct derivation_quantites {
    template< class INTER>
    void operator()(INTER &Inter,Param &process) const {
        derivation(Inter,process);
    }
};

/**\ingroup etape_lineaire
 \brief Etape Lineaire : Dérivation du deplacement sur les interfaces : cf derivation()
 */
struct derivation_quantites_sst {
template< class TV1,class TV2>
void operator()(TV1 &S,TV2 &Inter,Param &process) const {
  derivation(S,Inter,process);
 }
};

/** \ingroup  etape_lineaire
\brief Calcul du second membre micro provenant des sous-structures (mais pas des cotes). Ex : thermique, endommagement ...
 
*/
struct calcul_secmemb_micro_sst {
    template<class SST>
    void operator()(SST &S, Param &process, DataUser &data_user) const {
        //second membre prenant en compte le comportement thermique et la condition au pas de temps precedent ou les quantites chapeau:
        S.mesh.load();
        S.mesh->density=S.mesh.density;
        S.mesh.load_f_vol_e(data_user);
        S.mesh->f_vol=S.mesh.f_vol;
        S.f->set_mesh(S.mesh.m);
        S.f->assemble(false,true);
        S.fvol = S.f->get_sollicitation();
        S.mesh.unload();
    }
};

/** \ingroup etape_lineaire
\brief Programme principal pour l'étape Linéaire
 */
template<class TV1, class TV2, class GLOB>
void etape_lineaire(TV1 &S, TV2 &Inter,Param &process,GLOB &Global) {
    unsigned nb_threads=process.nb_threads;
     TicToc2 tic1,tic2;
     if (process.temps->pt==2) {tic1.start();tic2.start();}

    if (process.size == 1 or process.rank > 0) apply_mt(S,nb_threads,semilinstage1(),Inter,process);

   if (process.temps->pt==2) {crout << process.rank<<" : lineaire1 :"; tic1.stop();tic1.start();}
//         if (process.size > 1) MPI_Barrier( MPI_COMM_WORLD );
//         cout << "Barriere2" << endl;
    if (process.multiscale->multiechelle ==1) {
//         cout << "Assemblage du second membre macro" << endl;
        Global.bigF.set(0.0);
        if (process.size == 1 or process.rank > 0) apply_mt(S,nb_threads,macroassemble(),Inter,*process.temps,Global);
   if (process.temps->pt==2) {crout << process.rank<<" : macroassemble :"; tic1.stop();}
//        cout << "Attention : " << process.rank << " " << Global.bigF << endl;
         if (process.size > 1) MPI_Barrier( MPI_COMM_WORLD );
//         cout << "Barriere3" << endl;
        //Deploiement de bigF sur le master
         if (process.temps->pt==2) {tic1.start();}
         if (process.size > 1) SendbigF(process,Global.bigF);
   if (process.temps->pt==2) {crout << process.rank<<" : send bidF :"; tic1.stop();tic1.start();}
//         cout << "Attention : " << process.rank << " " << norm_2(Global.bigF) << endl;
//          if (process.size > 1) MPI_Barrier( MPI_COMM_WORLD );
//         cout << "Barriere4" << endl;
//         cout << "Resolution du probleme macro" << endl;
        if (process.rank==0) Global.resolmacro();

        if (process.temps->pt==2) {crout << process.rank<<" : resolmacro :"; tic1.stop();}
//         if (process.size > 1) MPI_Barrier( MPI_COMM_WORLD );
//         cout << "Barriere5" << endl;
        //Deploiement de bigW sur toutes les machines
//         if (process.size > 1) MPI_Barrier( MPI_COMM_WORLD );
        if (process.temps->pt==2) {tic1.start();}

//   if (process.temps->pt==2) {cout << process.rank<<" : synchro :"; tic1.stop();tic1.start();}
//          cout << process.rank << "  " << Global.bigW.size() << endl;
        if (process.size > 1) MPI_Bcast(Global.bigW.ptr(),Global.bigW.size() , MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (process.temps->pt==2) {crout << process.rank<<" : bcast bigW :"; tic1.stop();tic1.start();}
//         if (process.rank==0) cout << "DEBUG : bigW apres etape macro : " << norm_2(Global.bigW) << endl;
//         if (process.size > 1) MPI_Barrier( MPI_COMM_WORLD );
//         cout << "Barriere6" << endl;
        //double erreur=norm_2(Global.bigW)/Global.bigW.size();
        //Global.max_erreur=max(Global.max_erreur,erreur);
         if(process.multiscale->opti_multi==1 and (norm_2(Global.bigW)/Global.bigW.size()<=process.multiscale->erreur_macro) and process.latin->iter != 0)
//        if(process.multiscale->opti_multi==1 and erreur/Global.max_erreur<=process.multiscale->erreur_macro and process.latin->iter != 0)
          process.multiscale->multiechelle=0;
        //cout << norm_2(Global.bigW)/Global.bigW.size() << endl;
        //cout << "Extraction des deplacements macro - multiplicateur" << endl;
        if (process.size == 1 or process.rank > 0) apply_mt(S,nb_threads,interextrmacro(),Inter,*process.temps,Global);
   if (process.temps->pt==2) {crout << process.rank<<" : interextrmacro :"; tic1.stop();tic1.start();}
        //cout << "Etape micro 2" << endl;
//         if (process.size > 1) MPI_Barrier( MPI_COMM_WORLD );
//         cout << "Barriere7" << endl;
        if (process.size == 1 or process.rank > 0) apply_mt(S,nb_threads,semilinstage2(),Inter,process);
   if (process.temps->pt==2) {crout << process.rank<<" : lineaire2 :"; tic1.stop();tic1.start();}
    }
//         if (process.size > 1) MPI_Barrier( MPI_COMM_WORLD );
//         cout << "Barriere8" << endl;
//    if (process.size == 1 or process.rank > 0) apply_mt(S,nb_threads,reconstruction_quantites(),Inter,process);
    if (process.size == 1 or process.rank > 0) apply_mt(S,nb_threads,derivation_quantites_sst(),Inter,process);

   if (process.temps->pt==2) {crout << process.rank<<" : reconstruction :"; tic1.stop();crout << process.rank<<" : cout d une etape lineaire par pas de temps : "; tic2.stop();}
};

