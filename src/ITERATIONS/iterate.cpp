#include "iterate.h"

#include "LINEAR/linearstage.h"
#include "LOCAL/localstage.h"
#include "ERROR/calculate_error.h"

#include "../DEFINITIONS/MultiScaleData.h"

//fcts MPI
#include "../MPI/crout.h"
#include "../DEFINITIONS/ParallelisationData.h"
#include "../MPI/mpi_transactions.h"
#include "../../LMT/include/containers/evaluate_nb_cycles.h"

#include "../POSTTRAITEMENTS/create_file_pvd.h"
#include "../../LMT/include/containers/gnuplot.h"

using namespace LMT;
extern Crout crout;


void iterate_incr(Process &process, PointedSubstructures &S, VecInterfaces &Inter,PointedInterfaces &SubI, MacroProblem &Global,DataUser &data_user) {
    
    //phase iterative
    process.print_title(2,"    Début de la phase itérative :");
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
        if (process.parallelisation->is_master_cpu()) 
          std::cout << "    étape linéaire----------" << std::endl; 
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
        process.parallelisation->synchronisation();
        if (process.parallelisation->is_multi_cpu())
            SendRecvInter(process.parallelisation->intertoexchangebypro,Inter,process);
        crout << process.parallelisation->rank<< " : envoie des vecteurs d interface : ";
        tic.stop();
        tic.start();
        
        ///Si on a converge, on le dit aux interfaces cassables pour qu'elles mettent à jour leur comportement
        if (flag_convergence == 1 and process.nb_breakable > 0) {
            if (process.parallelisation->is_local_cpu()){
                for(unsigned q=0; q < SubI.size();q++){
                    if (SubI[q].comp == "Breakable")
                        SubI[q].convergence = 0; 
                }
            }
        }
        
        /// Etape locale
        if (process.parallelisation->is_master_cpu())
          std::cout << "    étape locale----------" << std::endl; 
        process.parallelisation->synchronisation();
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
            
            ///attention le paramètre est à déterminer (1e-4) pour etre optimal
            if ( i_iter > 0 and (process.latin->error[i_iter-1]-process.latin->error[i_iter] < process.latin->critere_erreur_auto_stop) )
                d_err++;
            else
                d_err=0;
            
            /// arret du processus si critere atteint (1 iteration supplementaire pour sauvegarder les deplacements
            if (process.latin->error[i_iter]<=process.latin->critere_erreur) {
                if (process.latin->critere_erreur_diss != 0) {
                    assign_t_post(S,SubI, process);
                    Vector frac;
                    frac.resize(2+Inter.size());
                    frac.set(0);
                    Vector fracfin;
                    fracfin.resize(2+Inter.size());
                    fracfin.set(0);
                    apply(S,calcerror_dissi_post(),Inter,frac,process);
                    if (process.parallelisation->is_multi_cpu())
                        MPI_Allreduce(frac.ptr(),fracfin.ptr(),frac.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                    else
                        fracfin=frac;
                    Scalar den=fracfin[1];
                    Scalar num=fracfin[0];
                    
                    process.print_data("\tErreur dissipation : ",std::abs(num-den)/std::abs(num));
                    
                    if ( std::abs(den)>1e-9 and std::abs(num-den)/std::abs(num) <= process.latin->critere_erreur_diss) {
                        flag_convergence=1;
                        process.latin->save_depl_SST=save_depl_SST;
                    }
                } else {
                    flag_convergence=1;
                    process.latin->save_depl_SST=save_depl_SST;
                }
            } else if ( d_err == 100 ) {
                /// sinon si sur plusieurs itérations de suite le delta d'erreur n'est pas assez grand ou negatif on s'arrete
                process.print("Arret par manque de convergence"); 
                flag_convergence=1;
                process.latin->save_depl_SST=save_depl_SST;
            }
            if (process.latin->iter==process.latin->nbitermax-1) {
                /// arret du processus nb d iteration max-1 atteint (1 iteration supplementaire pour sauvegarder les deplacements)
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
        //  GnuPlot gp;
        //  gp.plot(log10(process.latin->error[range(process.latin->iter+1)]));
        //  gp.wait();
        //     }
        
}