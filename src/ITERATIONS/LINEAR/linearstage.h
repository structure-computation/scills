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

/** \defgroup etape_lineaire Etape lin�aire
\ingroup LATIN
\brief Description des phases de l'�tape lin�aire
 
 Cette �tape de la strat�gie LATIN est constitu�e des phases suivantes :
 - pour chaque intervalle de temps on fait :
   - �tape lin�aire 1 : \ref semilinstage1
   - �tape macro (seulement en multi�chelle)
      - construction du second membre macro :  \ref macroassemble
      - r�solution du probl�me macro: MacroProblem::resolmacro()
      - extraction des valeurs du multiplicateur macro et assignation aux interfaces : \ref interextrmacro
      - r�solution du probl�me micro 2 : \ref semilinstage2
   - reconstruction des quantit�s � partir des deux �tapes :\ref  reconstruction_quantites
 - une fois la boucle sur le temps termin�e on effectue une relaxation des quantit�s :\ref  relaxation_quantites
 
 Rq : Pour la statique un seul intervalle de temps est r�alis�. (Le pas de temps n'a alors pas d'importance)
 
*/

/** \ingroup etape_lineaire
\relates reconstruction_quantites
\brief Etape Lineaire : Reconstruction du d�placement et des quantit�s d'interface
 
En sommant le champ solution \f$ q_1 \f$ de l'�tape lin�aire 1 et le champ solution \f$ q_2 \f$ de l'�tape lin�aire 2, on reconstruit le d�placement \f$ q \f$ puis par d�rivation on peut obtenir la vitesse pour l'intervalle de temps consid�r�.
De la m�me mani�re on reconstruit les efforts, d�placements et vitesses d'interface. Pour le mono�chelle, on n'utilise que le champ q1
 
Il est possible ici d'extraire les efforts macro ou micro, en utilisant le projecteur ad�quat (cf : Interface).
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
 
Ayant d�termin� pour chaque intervalle de temps la solution q pour les Sst, la vitesse \f$ \dot{W} ( ou d�placement \f$ W \f$ ) \f$ et effort \f$ F \f$ sur chaque interface, on proc�de � la relaxation, permettant d'obtenir une meilleure convergence.
 
On effectue donc l'op�ration suivante, avec mu=0.8 par d�faut (1 � la premi�re it�ration) :
\code q = (1-mu)*oldq + mu*q \endcode
o� oldq est le d�placement solution � l'it�ration pr�c�dente. 
On effectue de m�me avec les quantit�s d'interfaces, puis on met � jour les anciennes quantit�s.
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
 
Connaissant la vitesse pour chaque pas de temps, il est possible de d�terminer le d�placement � chaque piquet de temps en prenant en compte le d�placement initial (t=0) au d�part.
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
\brief Programme principal pour l'�tape Lin�aire
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
