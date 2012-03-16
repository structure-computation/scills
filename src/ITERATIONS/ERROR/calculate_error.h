#ifndef CALCUL_ERREUR
#define CALCUL_ERREUR
using namespace LMT;
#include "types_erreurs.h"
#include "assign_quantities_current_to_old.h"

/** \defgroup calcul_erreur Calcul de l'erreur
\ingroup LATIN
\brief Description des phases du calcul de l'erreur
 
Pour chaque pas de temps, on détermine l'erreur permettant de savoir si la stratégie converge. Plusieurs types d'erreurs ou indicateurs sont disponibles :
- erreur en énergie : (LATIN::type_error == "energie") : calcerror_ener
- erreur en direction de recherche : (LATIN::type_error == "ddr") : calcerror_ddr
Pour l'erreur choisie on détermine le numérateur et dénominateur que l'on somme sur tous les pas de temps puis on calcule l'erreur pour l'itération latin donnée : \f$ erreur = \frac{\std::sqrt{num}}{\std::sqrt{den}} \f$
*/

/** \ingroup   calcul_erreur
\brief Fonction principale pour le calcul de l'erreur
*/
template<class TV1, class TV2,class GLOB>
void calcul_erreur_latin(TV1 &S, TV2 &Inter,Param &process, GLOB &Global) {
    TYPEREEL num=0.0,den=0.0,eps=1e-9;
    Vec<TYPEREEL> frac;frac.resize(2+Inter.size());frac.set(0.);
    Vec<TYPEREEL> fracfin;fracfin.resize(2+Inter.size());fracfin.set(0.);
    if (process.size == 1 or process.rank >0) {
        if(process.nom_calcul=="latin") {
            for(unsigned imic=1;imic<=process.temps->nbpastemps;imic++) {
                process.temps->pt=imic;
                process.temps->pt_cur=imic;
                if(process.latin->type_error=="energie") {
                    //erreur en energie sur les interfaces
                    apply(S,calcerror_ener(),Inter,frac,process);
                } else if(process.latin->type_error=="ddr") {
                    //indicateur d'erreur sur les directions de recherche
                    apply(S,calcerror_ddr(),Inter,frac,process);
                } else if(process.latin->type_error=="dissipation") {
                    break;
                } else if(process.latin->type_error=="residu_depl") {
                    //indicateur d'erreur sur la dissipation
                    apply(S,calcerror_residu_depl(),Inter,frac,process);
                    if (process.size == 1 or process.rank == 1) frac[1]=1.0;
                } else {
                    std::cout << process.latin->type_error << " : Type d'erreur non reconnu : par defaut erreur sur les ddr" << std::endl;
                    apply(S,calcerror_ddr(),Inter,frac,process);
                }
            }
            if (process.latin->type_error=="dissipation") {
                apply(S,calcerror_dissi(),Inter,frac,process);
            }
        } else if (process.nom_calcul=="incr") {
//            for(unsigned imic=1;imic<=process.temps->nbpastemps;imic++) {//MOD A LA CON POUR FAIRE L ERREUR LATIN UNIQUEMENT SUR LES CYCLES ET PAS LE PRECHARGEMENT
            for(unsigned imic=101;imic<=process.temps->nbpastemps;imic++) {
                process.temps->pt=imic;
                process.temps->pt_cur=imic;
                if(process.latin->type_error=="energie") {
                    //erreur en energie sur les interfaces
                    apply(S,calcerror_ener(),Inter,frac,process);
                } else if(process.latin->type_error=="ddr") {
                    //indicateur d'erreur sur les directions de recherche
                    apply(S,calcerror_ddr_post(),Inter,frac,process);
                } else if(process.latin->type_error=="dissipation") {
                    //indicateur d'erreur sur la dissipation
                    apply(S,calcerror_ddr_post(),Inter,frac,process);
                } else if(process.latin->type_error=="residu_depl") {
                    //indicateur d'erreur sur la dissipation
                    apply(S,calcerror_residu_depl_post(),Inter,frac,process);
                    if (process.size == 1 or process.rank == 1) frac[1]=1.0;
                } else {
                    std::cout << process.latin->type_error << " : Type d'erreur non reconnu : par defaut erreur sur les ddr" << std::endl;
                    apply(S,calcerror_ddr_post(),Inter,frac,process);
                }
            }
        }
    }

    if (process.size>1) MPI_Allreduce(frac.ptr(),fracfin.ptr(),frac.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    else fracfin=frac;
        den=fracfin[1];
        num=fracfin[0]; 
        if (std::abs(den)>=eps)
          if (process.latin->type_error!="dissipation") {
            process.latin->error[process.latin->iter] =std::sqrt(num)/std::sqrt(den);
//             if (process.rank == 0) std::cout << "\tParticipation des interfaces : " << fracfin[range(2,(int)Inter.size()+2)] << std::endl;
          }
          else process.latin->error[process.latin->iter] =std::abs(num-den)/std::abs(num);
        else
            process.latin->error[process.latin->iter] =1;
};

/** \defgroup calcul_erreur_incr Calcul de l'erreur
\ingroup Incrementale
\brief Description des phases du calcul de l'erreur
 
Pour chaque pas de temps, on détermine l'erreur permettant de savoir si la stratégie converge. Plusieurs types d'erreurs ou indicateurs sont disponibles :
- erreur en énergie : (LATIN::type_error == "energie") : calcerror_ener
- erreur en direction de recherche : (LATIN::type_error == "ddr") : calcerror_ddr
Pour l'erreur choisie on détermine le numérateur et dénominateur que l'on somme sur tous les pas de temps puis on calcule l'erreur pour l'itération latin donnée : \f$ erreur = \frac{\std::sqrt{num}}{\std::sqrt{den}} \f$
*/

/** \ingroup   calcul_erreur_incr
\brief Fonction principale pour le calcul de l'erreur
*/
template<class TV1, class TV2,class GLOB>
void calcul_erreur_incr(TV1 &S, TV2 &Inter,Param &process, GLOB &Global) {
    TYPEREEL num=0.0,den=0.0,eps=1e-9;
    Vec<TYPEREEL> frac;frac.resize(2+Inter.size());frac.set(0.);
    Vec<TYPEREEL> fracfin;fracfin.resize(2+Inter.size());fracfin.set(0.);
    process.temps->pt=1;
    if(process.latin->type_error=="energie") {
        //erreur en energie sur les interfaces
        apply(S,calcerror_ener(),Inter,frac,process);
    } else if(process.latin->type_error=="ddr") {
        //indicateur d'erreur sur les directions de recherche
        apply(S,calcerror_ddr(),Inter,frac,process);
    } else if(process.latin->type_error=="dissipation") {
        //indicateur d'erreur sur la dissipation
        apply(S,calcerror_dissi_post(),Inter,frac,process);
    } else if(process.latin->type_error=="residu_depl") {
        //indicateur d'erreur sur la dissipation
        apply(S,calcerror_residu_depl(),Inter,frac,process);
        if (process.size == 1 or process.rank == 1) frac[1]=1.0;
    } else {
        std::cout << process.latin->type_error << " : Type d'erreur non reconnu : par defaut erreur sur les ddr" << std::endl;
        apply(S,calcerror_ddr(),Inter,frac,process);
    }
    

    if (process.size>1) MPI_Allreduce(frac.ptr(),fracfin.ptr(),frac.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    else fracfin=frac;
    den=fracfin[1];
    num=fracfin[0]; 
    if (std::abs(den)>=eps)
      if (process.latin->type_error!="dissipation") {
      process.latin->error[process.latin->iter] =std::sqrt(num)/std::sqrt(den);
//       if (process.rank == 0) std::cout << "\tParticipation des interfaces : " << fracfin[range(2,(int)Inter.size()+2)]/num*100 << std::endl;
      }
      else process.latin->error[process.latin->iter] =std::abs(num-den)/std::abs(num);
    else
      process.latin->error[process.latin->iter] =1;


};
#endif //CALCUL_ERREUR
