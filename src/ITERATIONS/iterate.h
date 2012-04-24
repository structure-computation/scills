// fonction utilisees pour la procedure iterative
using namespace LMT;


#include "LINEAR/linearstage.h"
#include "LOCAL/localstage.h"
#include "ERROR/calculate_error.h"

//fcts MPI
#include "../MPI/crout.h"
#include "../DEFINITIONS/MPIParameters.h"
#include "../MPI/mpi_transactions.h"
#include "../../LMT/include/containers/evaluate_nb_cycles.h"

#include "../POSTTRAITEMENTS/create_file_pvd.h"
#include "../../LMT/include/containers/gnuplot.h"

extern Crout crout;

/** \defgroup LATIN Stratégie LATIN
\ingroup Strategie_iterative
\brief Procédure principale de la stratégie latin en quasistatique.
  
  On effectue la boucle latin tant que le critère d'erreur n'est pas atteint (LatinParameters::critere_erreur) ou pour un nombre maximal d'itérations (LatinParameters::nbitermax). Celle ci est décrite dans les modules suivants :
- \ref etape_lineaire pour chaque piquet de temps successivement
- \ref relaxation_quantities pour chaque piquet de temps successivement
- \ref etape_locale pour chaque piquet de temps successivement
- \ref calcul_erreur_latin
 
 Pendant la boucle latin, on assigne une valeur différente au coefficient LatinParameters::mu pour la relaxation. A la première itération on lui assigne 1.0 sinon on lui affecte la valeur de LatinParameters::facteur_relaxation (donné par l'utilisateur). 
Le paramètre LATIN::list_error permet de lister l'erreur latin au cours des itérations.
*/
void iterate_latin(Process &process, Vec<VecPointedValues<Sst> > &S, Vec<Interface> &Inter,Vec<VecPointedValues<Interface > > &SubI, MacroProblem &Global, DataUser &data_user) {
    //phase iterative
    bool flag_convergence=0;
    bool save_depl_SST=process.latin->save_depl_SST;
    process.latin->save_depl_SST=0;
    // A mettre ici pour la thermique seulement
    /*if (process.size > 1 )
        M PI_Barrier(MPI_COMM_WORLD);*/
    if(not process.plasticite) // On n'initialise qu'une seul fois les efforts volumiques si ils sont 
        apply_mt(S,process.nb_threads,Calcul_2nd_membre_micro1_sst(),process, data_user);
    bool multiechelle=process.multiscale->multiechelle;
    TicToc2 tic1;
    for(unsigned i_iter=0;i_iter<(unsigned)process.latin->nbitermax+1;++i_iter) {
        if (process.rank ==0)
            tic1.start();
        process.latin->iter=i_iter;

        TicToc2 tic;
        tic.start();

        /// Etape lineaire
        for(unsigned i_pt=1;i_pt<process.temps->nbpastemps+1;i_pt++) {
            process.temps->pt=i_pt;
            process.temps->pt_cur=i_pt;
            etape_lineaire(S,Inter,process,Global,data_user);
        }
        crout << process.rank<< " : etape lineaire : " ;
        tic.stop();
        tic.start();

        /// Relaxation
        if (process.latin->iter ==0 and process.reprise_calcul==0)
            process.latin->mu=1.0;
        else
            process.latin->mu=process.latin->facteur_relaxation;
        for(unsigned i_pt=1;i_pt<=process.temps->nbpastemps;i_pt++) {
            process.temps->pt=i_pt;
            process.temps->pt_cur=i_pt;
            apply_mt(S,process.nb_threads,relaxation_quantites(),Inter,process);
        }
        crout << process.rank<< " : relaxation : " ;
        tic.stop();
        
        /// Echange des grandeurs macro calculees
        if (process.size > 1 )
            MPI_Barrier(MPI_COMM_WORLD);
        tic.start();
        if (process.size > 1 )
            SendRecvInter(process.multi_mpi->intertoexchangebypro,Inter,process);
        crout << process.rank<< " : envoie des vecteurs d interface : ";
        tic.stop();
        tic.start();

        /// Etape locale
        if (process.size==1 or process.rank>0)
            for(unsigned i_pt=1;i_pt<=process.temps->nbpastemps;i_pt++) {
                process.temps->pt=i_pt;
                process.temps->pt_cur=i_pt;
                etape_locale(SubI,S,process);
            }
        crout << process.rank<< " : etape locale : ";
        tic.stop();

        /// Calcul d'erreur
        if (process.size > 1 )
            MPI_Barrier(MPI_COMM_WORLD);
        tic.start();
        calcul_erreur_latin(S,Inter,process,Global);
        crout << process.rank<< " : calcul erreur : ";
        tic.stop();
        tic.start();

        if (process.latin->list_error==1)
            if (process.rank == 0)
                std::cout << "Erreur - iteration - " << i_iter << " : " << process.latin->error[i_iter] << endl;

        /// Test de convergence
        if(flag_convergence==1) {
            if (process.rank == 0)
                std::cout << "**** Sortie Critere atteint en "<< i_iter <<" iterations ****" <<  process.latin->error[i_iter] <<    endl;
            break;
        }
        if(i_iter==5 and min(process.latin->error)>=0.9999) {
            if (process.rank == 0)
                std::cout << "**** Sortie car chargement nul impose ******" << endl;
        }

        /// arret du processus si critere atteint (1 iteration supplementaire pour sauvegarder les deplacements)
        if (process.latin->error[i_iter]<=process.latin->critere_erreur) {
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
    if (process.rank ==0 and (find(process.affichage->display_fields,LMT::_1==Sc2String("contact"))))
        for (unsigned pt=1;pt<=process.temps->nbpastemps;pt++)
            create_file_pvd(process,"contact_",pt);
}

/** \defgroup Incrementale Stratégie Incrémentale
\ingroup Strategie_iterative
\brief Procédure principale de la stratégie incrémentale en quasistatique.
 
 
 Pour chaque piquet de temps, on effectue une résolution itérative. On assigne une valeur différente au coefficient LatinParameters::mu pour la relaxation. A la première itération on lui assigne 1.0 sinon on lui affecte la valeur de LatinParameters::facteur_relaxation (donné par l'utilisateur).
 
 On effectue la boucle latin tant que le critère d'erreur n'est pas atteint (LatinParameters::critere_erreur) ou pour un nombre maximal d'itérations (LatinParameters::nbitermax). Celle ci est décrite dans les modules suivants :
- \ref etape_lineaire pour le piquet de temps 1
- \ref relaxation_quantites pour le piquet de temps 1
- \ref etape_locale pour le piquet de temps 1
- \ref calcul_erreur_incr
 
 Le paramètre LatinParameters::list_error permet de lister l'erreur latin au cours des itérations.
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
        apply_mt(S,process.nb_threads,Calcul_2nd_membre_micro1_sst(),process, data_user);
    for(unsigned i_iter=0;i_iter<(unsigned)process.latin->nbitermax+1;++i_iter) {
        if (process.rank ==0)
            tic1.start();
      
        process.latin->iter=i_iter;
        TicToc2 tic;
        tic.start();
        
        /// Etape lineaire
        etape_lineaire(S,Inter,process,Global,data_user);
        crout << process.rank<< " : etape lineaire : " ;
        tic.stop();
        tic.start();
        
        /// Relaxation
        if (process.latin->iter ==0 )
            process.latin->mu=1.0;
        else
            process.latin->mu=process.latin->facteur_relaxation;
        apply_mt(S,process.nb_threads,relaxation_quantites(),Inter,process);
        crout << process.rank<< " : relaxation : " ;
        tic.stop();
        tic.start();
        
        /// Echange des grandeurs macro calculees
        if (process.size > 1 )
            SendRecvInter(process.multi_mpi->intertoexchangebypro,Inter,process);
        crout << process.rank<< " : envoie des vecteurs d interface : ";
        tic.stop();
        tic.start();
        
        ///Si on a converge, on le dit aux interfaces cassables pour qu'elles mettent à jour leur comportement
        if (flag_convergence == 1 and process.nb_breakable > 0) {
            if (process.size == 0 or process.rank > 0){
                for(unsigned q=0; q < SubI.size();q++){
                    if (SubI[q].comp == "Breakable")
                        SubI[q].param_comp->convergence = 0; 
                }
            }
        }
        
        /// Etape locale
        if (process.size==1 or process.rank>0)
            etape_locale(SubI,S,process);
        crout << process.rank<< " : etape locale : ";
        tic.stop();
        tic.start();
        
        /// Calcul de l'erreur
        calcul_erreur_incr(S,Inter,process,Global);
        crout << process.rank<< " : calcul erreur : ";
        tic.stop();
        tic.start();
        if (process.latin->list_error==1)
          if (process.rank == 0 )
              std::cout << "Erreur - iteration - " << i_iter << " : " << process.latin->error[i_iter] << endl;
        if(flag_convergence==1) {
            if (process.rank == 0 )
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

                if (process.rank == 0) std::cout << "\tErreur dissipation : " << std::abs(num-den)/std::abs(num) << endl;
                
                if ( std::abs(den)>1e-9 and std::abs(num-den)/std::abs(num) <= process.latin->critere_erreur_diss) {
                    flag_convergence=1;
                    process.latin->save_depl_SST=save_depl_SST;
                }
            } else {
                flag_convergence=1;
                process.latin->save_depl_SST=save_depl_SST;
            }
        /// sinon si sur plusieurs itérations de suite le delta d'erreur n'est pas assez grand ou negatif on s'arrete
        } else if ( d_err == 10 ) {
            if (process.rank == 0) std::cout << "Arret par manque de convergence" << std::endl; 
                flag_convergence=1;
                process.latin->save_depl_SST=save_depl_SST;
        }
        /// arret du processus nb d iteration max-1 atteint (1 iteration supplementaire pour sauvegarder les deplacements)
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

