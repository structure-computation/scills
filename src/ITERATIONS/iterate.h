// fonction utilisees pour la procedure iterative
#include "../all_declarations.h"
using namespace LMT;


#include "LINEAR/linearstage.h"
#include "LOCAL/localstage.h"
#include "ERROR/calculate_error.h"

//fcts MPI
#include "../MPI/crout.h"
#include "../DEFINITIONS/ParallelisationParameters.h"
#include "../MPI/mpi_transactions.h"
#include "../../LMT/include/containers/evaluate_nb_cycles.h"

#include "../POSTTRAITEMENTS/create_file_pvd.h"
#include "../../LMT/include/containers/gnuplot.h"

extern Crout crout;

/** \defgroup Incrementale Strat�gie Incr�mentale
\ingroup Strategie_iterative
\brief Proc�dure principale de la strat�gie incr�mentale en quasistatique.
 
 
 Pour chaque piquet de temps, on effectue une r�solution it�rative. On assigne une valeur diff�rente au coefficient LatinParameters::mu pour la relaxation. A la premi�re it�ration on lui assigne 1.0 sinon on lui affecte la valeur de LatinParameters::facteur_relaxation (donn� par l'utilisateur).
 
 On effectue la boucle latin tant que le crit�re d'erreur n'est pas atteint (LatinParameters::critere_erreur) ou pour un nombre maximal d'it�rations (LatinParameters::nbitermax). Celle ci est d�crite dans les modules suivants :
- \ref etape_lineaire pour le piquet de temps 1
- \ref relaxation_quantites pour le piquet de temps 1
- \ref etape_locale pour le piquet de temps 1
- \ref calcul_erreur_incr
 
 Le param�tre LatinParameters::list_error permet de lister l'erreur latin au cours des it�rations.
*/
void iterate_incr(Process &process, Vec<VecPointedValues<Sst> > &S, Vec<Interface> &Inter,Vec<VecPointedValues<Interface> > &SubI, MacroProblem &Global,DataUser &data_user) {

    //phase iterative
    bool flag_convergence=0;
    bool multiechelle=process.multiscale->multiechelle;
    process.latin->save_depl_SST=true;
    bool save_depl_SST=process.latin->save_depl_SST;
    unsigned d_err = 0;
    TicToc2 tic1;
    if(not process.plasticite) // On n'initialise les seconds membres micro qu'une seul fois si ils ne dependent pas de l'etat de la Sst
        apply_mt(S,process.parallelisation->nb_threads,Calcul_2nd_membre_micro1_sst(),process, data_user);
    for(unsigned i_iter=0;i_iter<(unsigned)process.latin->nbitermax+1;++i_iter) {
        if (process.parallelisation->is_master_cpu())
            tic1.start();
      
        process.latin->iter=i_iter;
        TicToc2 tic;
        tic.start();
        
        /// Etape lineaire
        etape_lineaire(S,Inter,process,Global,data_user);
        crout << process.parallelisation->rank<< " : etape lineaire : " ;
        tic.stop();
        tic.start();
        
        /// Relaxation
        if (process.latin->iter ==0 )
            process.latin->mu=1.0;
        else
            process.latin->mu=process.latin->facteur_relaxation;
        apply_mt(S,process.parallelisation->nb_threads,relaxation_quantites(),Inter,process);
        crout << process.parallelisation->rank<< " : relaxation : " ;
        tic.stop();
        tic.start();
        
        /// Echange des grandeurs macro calculees
        if (process.parallelisation->is_multi_cpu())
            SendRecvInter(process.parallelisation->intertoexchangebypro,Inter,process);
        crout << process.parallelisation->rank<< " : envoie des vecteurs d interface : ";
        tic.stop();
        tic.start();
        
        ///Si on a converge, on le dit aux interfaces cassables pour qu'elles mettent � jour leur comportement
        if (flag_convergence == 1 and process.nb_breakable > 0) {
            if (process.parallelisation->is_local_cpu()){
                for(unsigned q=0; q < SubI.size();q++){
                    if (SubI[q].comp == "Breakable")
                        SubI[q].param_comp->convergence = 0; 
                }
            }
        }
        
        /// Etape locale
        if (process.parallelisation->is_local_cpu())
            etape_locale(SubI,S,process);
        crout << process.parallelisation->rank<< " : etape locale : ";
        tic.stop();
        tic.start();
        
        /// Calcul de l'erreur
        calcul_erreur_incr(S,Inter,process,Global);
        crout << process.parallelisation->rank<< " : calcul erreur : ";
        tic.stop();
        tic.start();
        if (process.latin->list_error==1)
          if (process.parallelisation->is_master_cpu())
              std::cout << "Erreur - iteration - " << i_iter << " : " << process.latin->error[i_iter] << endl;
        if(flag_convergence==1) {
            if (process.parallelisation->is_master_cpu())
                std::cout << "**** Sortie Critere atteint en "<< i_iter <<" iterations ****" <<  process.latin->error[i_iter] <<    endl;
            break;
        }
        
        ///attention le param�tre est � d�terminer (1e-4) pour etre optimal
        if ( i_iter > 0 and (process.latin->error[i_iter-1]-process.latin->error[i_iter] < process.latin->critere_erreur_auto_stop) )
            d_err++;
        else
            d_err=0;
        
        /// arret du processus si critere atteint (1 iteration supplementaire pour sauvegarder les deplacements
        if (process.latin->error[i_iter]<=process.latin->critere_erreur) {
            if (process.latin->critere_erreur_diss != 0) {
                assign_t_post(S,SubI, process);
                Vec<TYPEREEL> frac;
                frac.resize(2+Inter.size());
                frac.set(0);
                Vec<TYPEREEL> fracfin;
                fracfin.resize(2+Inter.size());
                fracfin.set(0);
                apply(S,calcerror_dissi_post(),Inter,frac,process);
                if (process.parallelisation->is_multi_cpu())
                    MPI_Allreduce(frac.ptr(),fracfin.ptr(),frac.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                else
                    fracfin=frac;
                TYPEREEL den=fracfin[1];
                TYPEREEL num=fracfin[0];

                if (process.parallelisation->is_master_cpu()) std::cout << "\tErreur dissipation : " << std::abs(num-den)/std::abs(num) << endl;
                
                if ( std::abs(den)>1e-9 and std::abs(num-den)/std::abs(num) <= process.latin->critere_erreur_diss) {
                    flag_convergence=1;
                    process.latin->save_depl_SST=save_depl_SST;
                }
            } else {
                flag_convergence=1;
                process.latin->save_depl_SST=save_depl_SST;
            }
        /// sinon si sur plusieurs it�rations de suite le delta d'erreur n'est pas assez grand ou negatif on s'arrete
        } else if ( d_err == 10 ) {
            if (process.parallelisation->is_master_cpu()) std::cout << "Arret par manque de convergence" << std::endl; 
                flag_convergence=1;
                process.latin->save_depl_SST=save_depl_SST;
        }
        /// arret du processus nb d iteration max-1 atteint (1 iteration supplementaire pour sauvegarder les deplacements)
        if (process.latin->iter==process.latin->nbitermax-1) {
            flag_convergence=1;
            process.latin->save_depl_SST=save_depl_SST;
        }
        crout << process.parallelisation->rank<< " : post erreur : ";
        tic.stop();
        if (process.parallelisation->is_master_cpu()) {
          crout << "\tCout d'une iteration LATIN : ";
          tic1.stop();
        }
    }
    if (process.parallelisation->is_master_cpu())
        if (process.latin->error[ process.latin->iter]>=process.latin->critere_erreur)
            std::cout << "**** Sortie Nb iter max (pour info err=  ****" << process.latin->error[process.latin->iter] << ")" << endl;
    process.multiscale->multiechelle=multiechelle;
    
//     if (process.parallelisation->is_master_cpu()){
//         std::cout << process.latin->error[range(process.latin->iter+1)] << std::endl;
// 	GnuPlot gp;
// 	gp.plot(log10(process.latin->error[range(process.latin->iter+1)]));
// 	gp.wait();
//     }

}

