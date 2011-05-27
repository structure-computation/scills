//librairies Hugo
#include "containers/mat.h"
#include "containers/vecpointedvalues.h"

//fichiers de definition des variables
#include "definition_PARAM_COMP_INTER.h"
#include "definition_PARAM.h"
#include "definition_PARAM_LATIN.h"
#include "definition_PARAM_MULTI.h"
#include "definition_PARAM_TEMPS.h"
#include "definition_GLOB.h"
#include "definition_SST_time.h"
#include "definition_INTER_time.h"
#include "definition_CL.h"
#include "definition_PARAM_AFFICHAGE.h"

// fonctions speciales math
#include "algebre.h"
#include "utilitaires.h"

// fonction utilisees pour la phase d'initialisation ou affectation des valeurs des CL puis phase iterative
#include "allocate.h"
#include "save_read_data.h"
#include "iterate.h"
#include "prelocalstage.h"
#include "assign_quantities_current_to_old.h"
#include "calculate_error.h"

#include "modification_sst_inter_behaviour.h"

//fonctions pour les envoies des données par mpi
#include "crout.h"


using namespace LMT;
using namespace std;

/**\defgroup Strategie_iterative Stratégie Itérative
\brief Stratégie itérative
 
Dans cette partie, on cherche à résoudre le problème de statique ou quasistatique. Deux méthodes sont envisagées : 
- on effectue une résolution latin à savoir une boucle itérative dans laquelle on cherche la solution de manière incrémentale sur tous les pas de temps : multiscale_iterate_latin(). 
- on effectue la boucle en temps et pour chaque piquet de temps on cherche une solution de manière itérative : multiscale_iterate_incr()
 
*/

/** \ingroup   LATIN
\brief Procédure principale itérative pour la résolution LATIN
 
Cette procédure est constituée des fonctions suivantes :
   - allocate_quantities() : On alloue dans un premier temps la mémoire pour chaque quantités de Sst, Interface et Global pour chaque piquet de temps:
   - assign_CL_values_space_time() : on assigne les valeurs des quantités chapeau sur le bord à partir des conditions aux limites pour chaque piquet de temps.
   - iterate_latin() : on applique la stratégie itérative latin : \ref LATIN*
*/
template <class TV1,class TV2, class TV3, class TV4, class TV5, class TV6>
void multiscale_iterate_latin(TV1 &S,TV2 &SubS, TV3 &Inter, TV4 &SubI, Param &process, TV5 &Global, TV6 &CL) {
  if(process.latin->alloc_quantites==1) {
        //1ere phase : allocations et initialisation des quantites
    if (process.rank == 0)
      std::cout << " Allocations des quantites d'interfaces et SST" << endl;
    if(process.reprise_calcul!=2) allocate_quantities(SubS,SubI,process,Global);//on alloue si on reprend un calcul a partir d un resultat en memoire
#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Memoire apres allocations : ").c_str(),1);
#endif
        //Assignation des Conditions aux limites pour chaque intervalle de temps et chaque interface de bord
    if(process.reprise_calcul==1) {
      read_data_sst(S, process);
      read_data_inter(Inter, process);
      calcul_erreur_latin(SubS, Inter, process, Global);
      if (process.rank == 0)
        std::cout << "Erreur initiale apres relecture : " << process.latin->error[process.latin->iter] << endl;
    }
    if (process.size == 1 or process.rank>0)
      assign_CL_values_space_time_latin(SubI, CL, process);
    if(process.reprise_calcul>0)
      calcul_erreur_latin(SubS, Inter, process, Global);
    if(process.reprise_calcul>0)
      recopie_old_from_new(Inter,process);
    if(process.reprise_calcul>0)
      if (process.rank == 0)
        std::cout << "Erreur initiale apres assignation : " << process.latin->error[process.latin->iter] << endl;
  }
  if (process.rank == 0)
    std::cout << " Processus iteratif " << endl;
  iterate_latin(process,SubS,Inter,SubI,Global);
};


/** \ingroup   Incrementale
\brief Procédure principale itérative pour la résolution du problème de manière incrémentale
 
Cette procédure est constituée des étapes suivantes :
- allocate_quantities() : On alloue dans un premier temps la mémoire pour chaque quantités de Sst, Interface et Global (pas de temps 0 et 1):
- pour chaque piquet de temps :
   - assign_CL_values_space_time() : on assigne les valeurs des quantités chapeau sur le bord à partir des conditions aux limites pour le piquet de temps concerné
   - iterate_incr() : on effectue une boucle itérative \ref Incrementale
   - assign_quantities_current_to_old() : on met à jour les valeurs 0 à partir des valeurs 1
*/
template <class TV1,class TV2, class TV3, class TV4, class TV5, class TV6>
void multiscale_iterate_incr(TV1 &S,TV2 &SubS, TV3 &Inter, TV4 &SubI, Param &process, TV5 &Global, TV6 &CL) {
    //1ere phase : allocations et initialisation des quantites
    if (process.rank == 0)
        std::cout << " Allocations des quantites d'interfaces et SST" << endl;
    if(process.reprise_calcul!=2) allocate_quantities(SubS,SubI,process,Global);
#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.rank)+" : Memoire apres allocations : ").c_str(),1);
#endif
    if(process.reprise_calcul==1) {
        read_data_sst(S, process);
        read_data_inter(Inter, process);
        calcul_erreur_latin(SubS, Inter, process, Global);
    }
    if(process.reprise_calcul>0)
      calcul_erreur_latin(SubS, Inter, process, Global);
    if(process.reprise_calcul>0)
      if (process.rank == 0)
        std::cout << "Erreur initiale avant reprise : " << process.latin->error[process.latin->iter] << endl;


    if (process.reprise_calcul>0)
        for(unsigned i=0;i<Inter.size();i++)
            for(unsigned j=0;j<Inter[i].side.size();j++) {
                Inter[i].side[j].t[1]=Inter[i].side[j].t_post[1];
                Inter[i].side[j].t[0]=Inter[i].side[j].t_post[0];
                if (process.reprise_calcul==1) recopie_old_from_new_post(Inter,process);
            }

    if(process.reprise_calcul>0)
      calcul_erreur_incr(SubS, Inter, process, Global);
    if(process.reprise_calcul>0)
      if (process.rank == 0)
        std::cout << "Erreur initiale avant reprise : " << process.latin->error[process.latin->iter] << endl;
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);

    // A mettre ici pour la thermique seulement
    apply_mt(SubS,process.nb_threads,calcul_secmemb_micro_sst(),process);

    if (process.rank == 0)
        std::cout<<"Processus iteratif incremental" << std::endl;


    unsigned i_step=process.temps->step_cur;
    for(unsigned i_pt = 0 ; i_pt < process.temps->time_step[i_step].nb_time_step; i_pt++){
        process.temps->time_step[i_step].pt_cur=i_pt;
        process.temps->pt_cur+=1;
        if (process.rank == 0)
            std::cout << "----Piquet de temps " << process.temps->time_step[i_step].t_ini+(i_pt+1)*process.temps->time_step[i_step].dt << std::endl;

//         if (process.size >1 ) MPI_Barrier(MPI_COMM_WORLD);
//         std::cout << process.rank<< " : barrier2 : " << endl;

        //Assignation des Conditions aux limites pour chaque intervalle de temps et chaque interface de bord
        if (process.size == 1 or process.rank>0)
            assign_CL_values_space_time_incr(SubI, CL, process);

        for(int ic=0;ic<CL.size();ic++){
            std::cout << "ft " << CL[ic].ft << std::endl;
/*            std::cout <<"fspace " << CL[ic].fcts_spatiales[i_step]<< endl;*/
        }
        
/*        if (process.size >1 ) MPI_Barrier(MPI_COMM_WORLD);
        std::cout << process.rank<< " : barrier3 : " << endl;*/
        
        if(process.reprise_calcul>0)
          calcul_erreur_incr(SubS, Inter, process, Global);
        if(process.reprise_calcul>0)
          if (process.rank == 0)
            std::cout << "Erreur initiale apres assignation : " << process.latin->error[process.latin->iter] << endl;
        // A mettre ici dans le cas d'un comportement dependant du temps pour les sst
        //apply_mt(S,process.nb_threads,calcul_secmemb_micro_sst(),process);

        iterate_incr(process,SubS,Inter,SubI,Global);
        //assignation ptcur au ptold
        assign_quantities_current_to_old(SubS,SubI,process);
        //modification de certaines interfaces ou sst (exemple endommagement)
        //modification_sst_inter_behaviour(S,Inter,param_incr);
        if (process.rank == 0)
            std::cout << "---- fin piquet de temps " << process.temps->time_step[i_step].t_ini+(i_pt+1)*process.temps->time_step[i_step].dt << std::endl;
    }
    calcul_erreur_latin(SubS, Inter, process, Global);
    if (process.rank == 0)
        std::cout << "Erreur : " << process.latin->error[process.latin->iter] << endl;

};
// #endif

void fake_iterate() {
    XmlNode n;
    Param process;

    Vec<Interface<DIM,TYPEREEL> > Inter3;
    Vec<Sst<DIM,TYPEREEL> > S3;
    Vec<VecPointedValues<Sst<DIM,TYPEREEL> > > SubS3;
    Vec<VecPointedValues<Interface<DIM,TYPEREEL> > > SubI3;
    Glob<DIM,TYPEREEL> Global3;
    Vec<Boundary<DIM,TYPEREEL> > CL3;
    
    multiscale_iterate_latin(S3,SubS3, Inter3, SubI3, process, Global3,CL3);
    multiscale_iterate_incr(S3,SubS3, Inter3, SubI3, process, Global3,CL3);

}
