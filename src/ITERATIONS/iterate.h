// fonction utilisees pour la procedure iterative
using namespace LMT;


#include "linearstage.h"
#include "localstage.h"
#include "calculate_error.h"

//fcts MPI
#include "crout.h"
#include "definition_PARAM_MPI.h"
#include "mpi_transactions.h"
#include "containers/evaluate_nb_cycles.h"

#include "create_file_pvd.h"
#include "containers/gnuplot.h"

extern Crout crout;

/** \defgroup LATIN Stratégie LATIN
\ingroup Strategie_iterative
\brief Procédure principale de la stratégie latin en quasistatique.
  
On effectue la boucle latin tant que le critère d'erreur n'est pas atteint (LATIN::critere_erreur) ou pour un nombre maximal d'itérations (LATIN::nbitermax). Celle ci est décrite dans les modules suivants :
- \ref etape_lineaire pour chaque piquet de temps successivement
- \ref relaxation_quantities pour chaque piquet de temps successivement
- \ref etape_locale pour chaque piquet de temps successivement
- \ref calcul_erreur_latin
 
Pendant la boucle latin, on assigne une valeur différente au coefficient LATIN::mu pour la relaxation. A la première itération on lui assigne 1.0 sinon on lui affecte la valeur de LATIN::facteur_relaxation (donné par l'utilisateur). 
Le paramètre LATIN::list_error permet de lister l'erreur latin au cours des itérations.
*/
template<class TV1, class TV2,class TV3, class GLOBAL>
void iterate_latin(Param &process, TV1 &S, TV2 &Inter,TV3 &SubI, GLOBAL &Global, DataUser &data_user) {
    //phase iterative
    bool flag_convergence=0;
    bool save_depl_SST=process.latin->save_depl_SST;
    process.latin->save_depl_SST=0;
    // A mettre ici pour la thermique seulement
    /*    if (process.size > 1 )
            MPI_Barrier(MPI_COMM_WORLD);*/
    apply_mt(S,process.nb_threads,calcul_secmemb_micro_sst(),process, data_user);
    bool multiechelle=process.multiscale->multiechelle;
    TicToc2 tic1;
    for(unsigned i=0;i<(unsigned)process.latin->nbitermax+1;++i) {
        if (process.rank ==0)
            tic1.start();
        process.latin->iter=i;

        TicToc2 tic;
        tic.start();

        for(unsigned pt=1;pt<process.temps->nbpastemps+1;pt++) {
            process.temps->pt=pt;
            process.temps->pt_cur=pt;
            // A mettre ici dans le cas d'un comportement dependant du temps pour les sst
            //apply_mt(S,process.nb_threads,calcul_secmemb_micro_sst(),process);
            etape_lineaire(S,Inter,process,Global);
        }
        crout << process.rank<< " : etape lineaire : " ;
        tic.stop();
        tic.start();

        //adaptation du facteur de relaxation
        if (process.latin->iter ==0 and process.reprise_calcul==0)
            process.latin->mu=1.0;
        else
            process.latin->mu=process.latin->facteur_relaxation;
        for(unsigned pt=1;pt<=process.temps->nbpastemps;pt++) {
            process.temps->pt=pt;
            process.temps->pt_cur=pt;
            apply_mt(S,process.nb_threads,relaxation_quantites(),Inter,process);
        }
        crout << process.rank<< " : relaxation : " ;
        tic.stop();
        //std::cout << "\t Etape locale " << endl;
        //         if (process.size > 1 )
        if (process.size > 1 )
            MPI_Barrier(MPI_COMM_WORLD);
        tic.start();
        //         if (process.rank>0)
        if (process.size > 1 )
            SendRecvInter(process.multi_mpi->intertoexchangebypro,Inter,process);
        crout << process.rank<< " : envoie des vecteurs d interface : ";
        tic.stop();
        tic.start();

        if (process.size==1 or process.rank>0)
            for(unsigned pt=1;pt<=process.temps->nbpastemps;pt++) {
                process.temps->pt=pt;
                process.temps->pt_cur=pt;
                etape_locale(SubI,S,process);
            }


        crout << process.rank<< " : etape locale : ";
        tic.stop();

        if (process.size > 1 )
            MPI_Barrier(MPI_COMM_WORLD);
        tic.start();
        //std::cout << "\t Calcul d'erreur " << endl;
        //                 std::cout << process.rank<< " : temps mort synchro : " ; tic.stop();tic.start();
        calcul_erreur_latin(S,Inter,process,Global);

        crout << process.rank<< " : calcul erreur : ";
        tic.stop();
        tic.start();

        if (process.latin->list_error==1)
            if (process.rank == 0)
                std::cout << "Erreur - iteration - " << i << " : " << process.latin->error[i] << endl;

        if(flag_convergence==1) {
            if (process.rank == 0)
                std::cout << "**** Sortie Critere atteint en "<< i <<" iterations ****" <<  process.latin->error[i] <<    endl;
            break;
        }
        if(i==5 and min(process.latin->error)>=0.9999) {
            if (process.rank == 0)
                std::cout << "**** Sortie car chargement nul impose ******" << endl;
        }

        // arret du processus si critere atteint (1 iteration supplementaire pour sauvegarder les deplacements)
        if (process.latin->error[i]<=process.latin->critere_erreur) {
            if (process.latin->critere_erreur_diss != 0) {
                Vec<double> frac;
                frac.resize(2+Inter.size());
                frac.set(0);
                Vec<double> fracfin;
                fracfin.resize(2+Inter.size());
                fracfin.set(0);
                apply(S,calcerror_dissi(),Inter,frac,process);
                if (process.size>1)
                  MPI_Allreduce(frac.ptr(),fracfin.ptr(),frac.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                else
                    fracfin=frac;
                double den=fracfin[1];//quantite n
                double num=fracfin[0];//quantite chapeau

                if (process.rank == 0){
                    std::cout << "\tErreur dissipation : " << std::abs(num-den)/std::abs(num) << endl;
//                     std::cout << "\tParticipation des interfaces : " << std::abs(fracfin[range(2,(int)Inter.size()+2)])/sum(std::abs(fracfin[range(2,(int)Inter.size()+2)]))*100 << endl;
                }
                if ( std::abs(den)>1e-9 and std::abs(num-den)/std::abs(num) <= process.latin->critere_erreur_diss) {
                    flag_convergence=1;
                    process.latin->save_depl_SST=save_depl_SST;
                }
            } else {
                flag_convergence=1;
                process.latin->save_depl_SST=save_depl_SST;
            }
        }
        // arret du processus nb d iteration max-1 atteint (1 iteration supplementaire pour sauvegarder les deplacements)
        if (process.latin->iter==process.latin->nbitermax-1) {
            flag_convergence=1;
            process.latin->save_depl_SST=save_depl_SST;
        }

        crout << process.rank<< " : post erreur : ";
        tic.stop();
        if (process.rank ==0) {
            crout << "\tCout d'une iteration LATIN : ";
            tic1.stop();
        }
    }
    if (process.rank ==0)
        if (process.latin->error[ process.latin->iter]>=process.latin->critere_erreur)
            std::cout << "**** Sortie Nb iter max (pour info err=  ****" << process.latin->error[process.latin->iter] << ")" << endl;
    process.multiscale->multiechelle=multiechelle;
    if (process.rank ==0 and (find(process.affichage->display_fields,LMT::_1==std::string("contact"))))
        for (unsigned pt=1;pt<=process.temps->nbpastemps;pt++)
            create_file_pvd(process,"contact_",pt);
}

/** \defgroup Incrementale Stratégie Incrémentale
\ingroup Strategie_iterative
\brief Procédure principale de la stratégie incrémentale en quasistatique.
 
 
Pour chaque piquet de temps, on effectue une résolution itérative. On assigne une valeur différente au coefficient LATIN::mu pour la relaxation. A la première itération on lui assigne 1.0 sinon on lui affecte la valeur de LATIN::facteur_relaxation (donné par l'utilisateur).
 
On effectue la boucle latin tant que le critère d'erreur n'est pas atteint (LATIN::critere_erreur) ou pour un nombre maximal d'itérations (LATIN::nbitermax). Celle ci est décrite dans les modules suivants :
- \ref etape_lineaire pour le piquet de temps 1
- \ref relaxation_quantites pour le piquet de temps 1
- \ref etape_locale pour le piquet de temps 1
- \ref calcul_erreur_incr
 
Le paramètre LATIN::list_error permet de lister l'erreur latin au cours des itérations.
*/
template<class TV1, class TV2,class TV3, class GLOBAL>
void iterate_incr(Param &process, TV1 &S, TV2 &Inter,TV3 &SubI, GLOBAL &Global) {

    //phase iterative
    process.temps->pt=1;
    bool flag_convergence=0;
    bool multiechelle=process.multiscale->multiechelle;
    bool save_depl_SST=process.latin->save_depl_SST;
    process.latin->save_depl_SST=0;
    unsigned d_err = 0;
    TicToc2 tic1;
    for(unsigned i=0;i<(unsigned)process.latin->nbitermax+1;++i) {
      if (process.rank ==0)
        tic1.start();
      
        process.latin->iter=i;
        
        TicToc2 tic;
        tic.start();
        
//         if (process.size >1 ) MPI_Barrier(MPI_COMM_WORLD);
//         std::cout << process.rank<< " : barrier1 : " << endl;

        //std::cout << "\t Etape lineaire " << endl;
        etape_lineaire(S,Inter,process,Global);
        
        
/*        if (process.size >1 ) MPI_Barrier(MPI_COMM_WORLD);
        crout << process.rank<< " : etape lineaire : " ;*/
        tic.stop();
        tic.start();

        //relaxation
//        if (process.latin->iter ==0  and process.reprise_calcul==0)
        if (process.latin->iter ==0 )
            process.latin->mu=1.0;
        else
            process.latin->mu=process.latin->facteur_relaxation;
        apply_mt(S,process.nb_threads,relaxation_quantites(),Inter,process);

        crout << process.rank<< " : relaxation : " ;
        tic.stop();
        tic.start();

        if (process.size > 1 )
            SendRecvInter(process.multi_mpi->intertoexchangebypro,Inter,process);
        crout << process.rank<< " : envoie des vecteurs d interface : ";
        tic.stop();
        tic.start();

        //Si on a converge, on le dit aux interfaces cassables pour qu'elles mettent à jour leur comportement
        if (flag_convergence == 1 && process.nb_breakable > 0) {
            if (process.size == 0 or process.rank > 0){
                for(unsigned q=0; q < SubI.size();q++){
                    if (SubI[q].comp == "Breakable")
                        SubI[q].param_comp->convergence = 0; 
                }
            }
        }
        //std::cout << "\t Etape locale " << endl;
#include "../../../LMTpp/include/util/unit_test.h"
        if (process.size==1 or process.rank>0)
            etape_locale(SubI,S,process);
        //etape locale sur les sous-structures
        //etape_locale(S);
        crout << process.rank<< " : etape locale : ";
        tic.stop();
        tic.start();
//         std::cout << "Avant lineaire" << endl;

        //std::cout << "\t Calcul d'erreur " << endl;
        calcul_erreur_incr(S,Inter,process,Global);

        crout << process.rank<< " : calcul erreur : ";
        tic.stop();
        tic.start();

        if (process.latin->list_error==1)
          if (process.rank == 0 )
                std::cout << "Erreur - iteration - " << i << " : " << process.latin->error[i] << endl;
        if(flag_convergence==1) {
            if (process.rank == 0 )
                std::cout << "**** Sortie Critere atteint en "<< i <<" iterations ****" <<  process.latin->error[i] <<    endl;
            break;
        }
        
        //attention le paramètre est à déterminer (1e-4) pour etre optimal
        if ( i > 0 && (process.latin->error[i-1]-process.latin->error[i] < process.latin->critere_erreur_auto_stop) )
	  d_err++;
	else
	  d_err=0;
        
        // arret du processus si critere atteint (1 iteration supplementaire pour sauvegarder les deplacements
        if (process.latin->error[i]<=process.latin->critere_erreur) {
            if (process.latin->critere_erreur_diss != 0) {
                assign_t_post(SubI, process);
                Vec<double> frac;
                frac.resize(2+Inter.size());
                frac.set(0);
                Vec<double> fracfin;
                fracfin.resize(2+Inter.size());
                fracfin.set(0);
                apply(S,calcerror_dissi_post(),Inter,frac,process);
                if (process.size>1)
                  MPI_Allreduce(frac.ptr(),fracfin.ptr(),frac.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                else
                    fracfin=frac;
                double den=fracfin[1];
                double num=fracfin[0];

                if (process.rank == 0) {
                    std::cout << "\tErreur dissipation : " << std::abs(num-den)/std::abs(num) << endl;
//                     std::cout << "\tParticipation des interfaces : " << std::abs(fracfin[range(2,(int)Inter.size()+2)])/sum(std::abs(fracfin[range(2,(int)Inter.size()+2)]))*100 << endl;
                }
                //if (process.rank == 0) std::cout << "\tnum " << num << " -- den " << den << endl;
                if ( std::abs(den)>1e-9 and std::abs(num-den)/std::abs(num) <= process.latin->critere_erreur_diss) {
                    flag_convergence=1;
                    process.latin->save_depl_SST=save_depl_SST;
                }
            } else {
                flag_convergence=1;
                process.latin->save_depl_SST=save_depl_SST;
            }
        // sinon si 5 itérations de suite le delta d'erreur n'est pas assez grand on s'arrete
        } else if ( d_err == 5 ) {
	      if (process.rank == 0) std::cout << "Arret par manque de convergence" << std::endl; 
                flag_convergence=1;
                process.latin->save_depl_SST=save_depl_SST;
	}
        // arret du processus nb d iteration max-1 atteint (1 iteration supplementaire pour sauvegarder les deplacements)
        if (process.latin->iter==process.latin->nbitermax-1) {
            flag_convergence=1;
            process.latin->save_depl_SST=save_depl_SST;
        }
        crout << process.rank<< " : post erreur : ";
        tic.stop();
        if (process.rank ==0) {
          crout << "\tCout d'une iteration LATIN : ";
          tic1.stop();
        }
    }
    if (process.rank == 0)
        if (process.latin->error[ process.latin->iter]>=process.latin->critere_erreur)
            std::cout << "**** Sortie Nb iter max (pour info err=  ****" << process.latin->error[process.latin->iter] << ")" << endl;
    process.multiscale->multiechelle=multiechelle;
    
//     if (process.rank == 0){
//         std::cout << process.latin->error[range(process.latin->iter+1)] << std::endl;
// 	GnuPlot gp;
// 	gp.plot(log10(process.latin->error[range(process.latin->iter+1)]));
// 	gp.wait();
//     }

}

