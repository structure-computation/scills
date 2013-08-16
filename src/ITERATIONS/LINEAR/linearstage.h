#include "Process.h"
#include "../../DEFINITIONS/structure_typedef.h"
#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/Interface.h"
#include "../../DEFINITIONS/LatinData.h"
#include "../../DEFINITIONS/MacroProblem.h"
#include "../../DEFINITIONS/MultiScaleData.h"
#include "../../COMPUTE/DataUser.h"

#include "../LOCAL/plasticite.h"
#include "../LOCAL/mesomodele.h"
#include "etape_lineaire_1.h"
#include "etape_lineaire_2.h"
#include "etape_macro.h"

//fcts MPI
#include "../../MPI/mpi_transactions.h"
using namespace LMT;
extern Crout crout;

/** \defgroup etape_lineaire Etape linéaire
\ingroup LATIN
\brief Description des phases de l'étape linéaire
 
 Cette étape de la stratégie LATIN est constituée des phases suivantes :
 - pour chaque intervalle de temps on fait :
   - étape linéaire 1 : \ref semilinstage1
   - étape macro (seulement en multiéchelle)
      - construction du second membre macro :  \ref macroassemble
      - résolution du problème macro: MacroProblem::resolmacro()
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
struct reconstruction_quantites {
    void operator()(Sst &S,VecInterfaces &Inter,Process &process) const {
        unsigned  pt=process.temps->pt;
        
        if(process.multiscale->multiechelle==1) {
            for(unsigned j=0;j<S.edge.size();++j) {
                unsigned data=S.edge[j].datanum;
                unsigned q=S.edge[j].internum;
                Inter[q].side[data].t[pt].Wp=(Inter[q].side[data].t[pt].W-Inter[q].side[data].t[pt-1].W)/process.temps->dt;
            }
        } else {
            for(unsigned j=0;j<S.edge.size();++j) {
                unsigned data=S.edge[j].datanum;
                unsigned q=S.edge[j].internum;
                Inter[q].side[data].t[pt].Wp=(Inter[q].side[data].t[pt].W-Inter[q].side[data].t[pt-1].W)/process.temps->dt;
            }
        }
    }
};


/** \ingroup etape_lineaire
    \relates derivation
    \brief Etape Lineaire : Reconstruction de la vitesse
 */
struct derivation_quantites_sst {
    void operator()(Sst &S,VecInterfaces &Inter,Process &process) const{
        unsigned  pt=process.temps->pt;
        
        for(unsigned j=0;j<S.edge.size();++j) {
            unsigned data=S.edge[j].datanum;
            unsigned q=S.edge[j].internum;
            Inter[q].side[data].t[pt].Wp=(Inter[q].side[data].t[pt].W-Inter[q].side[data].t[pt-1].W)/process.temps->dt;
        }
    }
};


/** \ingroup etape_lineaire
\relates relaxation_quantites
\brief Etape Lineaire : Relaxation
 
Ayant déterminé pour chaque intervalle de temps la solution q pour les Sst, la vitesse \f$ \dot{W} ( ou déplacement \f$ W \f$ ) \f$ et effort \f$ F \f$ sur chaque interface, on procède à la relaxation, permettant d'obtenir une meilleure convergence.
 
On effectue donc l'opération suivante, avec mu=0.8 par défaut (1 à la première itération) :
\code q = (1-mu)*oldq + mu*q \endcode
où oldq est le déplacement solution à l'itération précédente. 
On effectue de même avec les quantités d'interfaces, puis on met à jour les anciennes quantités.
*/
struct relaxation_quantites {
    void operator()(Sst &S,VecInterfaces &Inter,Process &process) const{
        TYPEREEL mu = process.latin->mu;
        unsigned pt = process.temps->pt;
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
};


/** \ingroup etape_lineaire
\relates integration_quantites
\brief Etape Lineaire : Integration
 
Connaissant la vitesse pour chaque pas de temps, il est possible de déterminer le déplacement à chaque piquet de temps en prenant en compte le déplacement initial (t=0) au départ.
*/
struct integration_quantites {
    void operator()(Interface &Inter,Process &process) const{
        for(unsigned j=0;j<Inter.side.size();j++) {
            for(unsigned pt=1;pt<=process.temps->nbpastemps;pt++)
                Inter.side[j].t[pt].W = Inter.side[j].t[pt-1].W + process.temps->dt *  Inter.side[j].t[pt].Wp;
        }
    }
};


/** \ingroup etape_lineaire
\relates derivation_quantites
\brief Etape Lineaire : Derivation
*/
struct derivation_quantites {
    void operator()(Interface &Inter,Process &process) const{
        for(unsigned j=0;j<Inter.side.size();j++) {
            for(unsigned pt=1;pt<=process.temps->nbpastemps;pt++)
                Inter.side[j].t[pt].Wp = (Inter.side[j].t[pt].W - Inter.side[j].t[pt-1].W)/process.temps->dt ;
        }
    }
};


/** \ingroup  etape_lineaire
 * \brief Calcul du second membre du probleme micro1 provenant des sous-structures (mais pas des cotes)
 *
 * Calcul le second membre du probleme micro 1 associe a la sous-structure S et le stocke sous S.fvol.
 * Prend en compte les efforts volumiques definis par l'utilisateur, ainsi que ceux definis dans la formulation. Ex : thermique, endommagement ...
 */
struct Calcul_2nd_membre_micro1_sst {
    void operator()(Sst &S, const Process &process, const DataUser &data_user) const{
        S.mesh.load();
        S.apply_behavior();
        process.Fvol->apply_on_sst(S);
        process.Tload->apply_on_sst(S);
        if(S.f == S.pb.formulation_plasticity_isotropy_stat_Qstat){
            upload_epsilon_p(S,process.temps->pt);
        }
        S.f->set_mesh(S.mesh.m);
        S.f->assemble(false,true);
        S.fvol = S.f->get_sollicitation();
        //S.mesh.unload();
    }
};

/** \ingroup etape_lineaire
\brief Programme principal pour l'étape Linéaire
 */
void etape_lineaire(PointedSubstructures &S, VecInterfaces &Inter,Process &process,MacroProblem &Global,DataUser &data_user){
    
    bool get_durations = (process.temps->pt==2);    /// Faut-il chronometrer les operations?
    unsigned nb_threads=process.parallelisation->nb_threads;
    TicToc2 tic1,tic2;
    if (get_durations) {tic1.start();tic2.start();}
    
    /// Reactualisation des operateurs (si necessaire)
    if(process.endommagement)
        apply_mt(S,nb_threads,reactualisation_rigidite(),Inter,process,data_user);
    
    /// Reactualisation du 2nd membre du probleme micro1
    if(process.plasticite)  // Uniquement si le 2nd membre de micro1 depend de l'etat du materiau
        apply_mt(S,nb_threads,Calcul_2nd_membre_micro1_sst(),process,data_user);
    
    /// Resolution du probleme micro1
    if (process.parallelisation->is_local_cpu()) apply_mt(S,nb_threads,semilinstage1(),Inter,process);
    if (get_durations) {crout << process.parallelisation->rank<<" : lineaire1 :"; tic1.stop();tic1.start();}
    
    
    process.parallelisation->synchronisation();
    if (process.multiscale->multiechelle ==1) {
        /// Assemblage du probleme macro
        Global.bigF.set(0.0);
        if (process.parallelisation->is_local_cpu()) apply_mt(S,nb_threads,macroassemble(),Inter,process,Global);
        if (get_durations) {crout << process.parallelisation->rank<<" : macroassemble :"; tic1.stop();}
        
        /// Deploiement de bigF sur le master
        if (get_durations) {tic1.start();}
        if (process.parallelisation->is_multi_cpu()) SendbigF(process,Global.bigF);
        if (get_durations) {crout << process.parallelisation->rank<<" : send bidF :"; tic1.stop();tic1.start();}
        
        /// Resolution du probleme macro
        if (process.parallelisation->is_master_cpu()) {
          Global.resolmacro();
        }
        //if (process.parallelisation->rank != 0)
//         {
//             PRINT(process.parallelisation->size);
//             PRINT(process.parallelisation->rank);
//             PRINT(Global.bigW.size());
//             PRINT(Global.bigW.ptr());
//             std::cout << std::endl;
//         }
        if (get_durations) {crout << process.parallelisation->rank<<" : resolmacro :"; tic1.stop();}
        
        //process.parallelisation->synchronisation();
        /// Deploiement de bigW sur toutes les machines
        if (get_durations) {tic1.start();}
        //if (process.parallelisation->is_multi_cpu()) MPI_Bcast(Global.bigW.ptr(),Global.bigW.size() , MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (process.parallelisation->size>1){
            MPI_Bcast(Global.bigW.ptr(),Global.bigW.size() , MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        if (get_durations) {crout << process.parallelisation->rank<<" : bcast bigW :"; tic1.stop();tic1.start();}
        if(process.multiscale->opti_multi==1 and (norm_2(Global.bigW)/Global.bigW.size()<=process.multiscale->erreur_macro) and process.latin->iter != 0)
            process.multiscale->multiechelle=0;
        //if (process.parallelisation->rank != 0)
//         {
//             PRINT(process.parallelisation->size);
//             PRINT(process.parallelisation->rank);
//             PRINT(Global.bigW.size());
//             PRINT(Global.bigW);
//             std::cout << std::endl;
//         }
        
        /// Assemblage du multiplicateur WtildeM
        if (process.parallelisation->is_local_cpu()) apply_mt(S,nb_threads,interextrmacro(),Inter,process,Global);
        if (get_durations) {crout << process.parallelisation->rank<<" : interextrmacro :"; tic1.stop();tic1.start();}
        
        /// Resolution du probleme micro2
        if (process.parallelisation->is_local_cpu()) apply_mt(S,nb_threads,semilinstage2(),Inter,process);
        if (get_durations) {crout << process.parallelisation->rank<<" : lineaire2 :"; tic1.stop();tic1.start();}
    }
    
    /// Reconstruction de la vitesse aux interfaces
    if (process.parallelisation->is_local_cpu()) apply_mt(S,nb_threads,derivation_quantites_sst(),Inter,process);
    
    if (get_durations){
        crout << process.parallelisation->rank<<" : reconstruction :"; 
        tic1.stop();
        crout << process.parallelisation->rank<<" : cout d une etape lineaire par pas de temps : ";
        tic2.stop();
    }
}
