#include "../DEFINITIONS/TimeData.h"
#include "../DEFINITIONS/LatinData.h"
#include "manipulate_quantities.h"
#include "LOCAL/plasticite.h"

using namespace LMT;


//allocation des differentes quantites
/** \ingroup   Strategie_iterative
 * \brief Allocation de memoire pour les sous-structures, interfaces et vecteurs macro. 
 * 
 * Pour chaque pas de temps les differents vecteurs sont initialises à 0. On n'alloue pas la memoire pour les differents vecteurs en chaque pas de temps des ssts car ils ne sont stockes qu'à convergence.
 */
void allocate_quantities_Sst_Inter(PointedSubstructures &SubS, PointedInterfaces &SubI,Process &process){
    /// Recuperation du nombre de pas de temps
    unsigned nbpastemps;
    if(process.nom_calcul=="incr"){
        nbpastemps=1;
    }else if (process.nom_calcul=="latin"){
        nbpastemps=process.temps->nbpastemps;
    }else{
        std::cerr << "Nom de calcul invalide : " << process.nom_calcul << endl;
        assert(0);
    }
    
    /// Allocation des quantites des Sst
    if (process.latin->save_depl_SST==1){
        for(unsigned i=0;i<SubS.size();++i){
            SubS[i].t.resize(nbpastemps+1);
            if(process.recopie_t_post==1){
                SubS[i].t=SubS[i].t_post;
            }
            if(process.parallelisation->is_local_cpu()){
                for(unsigned pt=0;pt<(nbpastemps+1);pt++){  ///L'initialisation commence a 0 et non plus a 1, pour les conditions
                    SubS[i].t[pt].allocate(SubS[i]);
                }
            }
        }
    }
    /* Affichage de debug
    for(int i = 0; i < SubS.size(); i++){
        SubS[i].affiche();
    }//*/
    
    /// Allocation des quantites des Interface
    if (process.parallelisation->is_local_cpu()){
        for(unsigned i=0;i<SubI.size();++i){
            for(unsigned j=0;j<SubI[i].side.size();++j){
                SubI[i].side[j].t.resize(nbpastemps+1);
                if(process.recopie_t_post==1){
                    SubI[i].side[j].t=SubI[i].side[j].t_post;
                }else{
                    for(unsigned pt=0;pt<(nbpastemps+1);pt++){
                        if (process.parallelisation->is_local_cpu()) {
                            SubI[i].side[j].t[pt].allocations(SubI[i].side[j].nodeeq.size()*DIM);
                        }
                    }
                }
            }
            ///allocation des parametres de comportements dependant du temps (exemple endommagement)
            SubI[i].t.resize(nbpastemps+1);
            for(unsigned j=0;j<nbpastemps+1;j++){
                SubI[i].t[j].allocate(SubI[i]);
            }
        }
    }
    
    if(process.recopie_t_post!=1){
        process.latin->error.resize(process.latin->nbitermax+1);
        process.latin->error.set(0.0);
    }
}


//allocation de t_post uniquement pour les sorties de resultat
void allocate_quantities_post(PointedSubstructures &SubS, PointedInterfaces &SubI,Process &process){
    unsigned nbpastemps;
    if(process.nom_calcul=="incr"){
        nbpastemps=1;
    }else if (process.nom_calcul=="latin"){
        nbpastemps=process.temps->nbpastemps;
    }else{
        std::cout << "Nom de clacul non implemente dans allocate.h" << std::endl;
        assert(0);
    }
    
    if (process.latin->save_depl_SST==1){
        for(unsigned i=0;i<SubS.size();++i){
            SubS[i].t.resize(nbpastemps+1);
            if(process.recopie_t_post==1){
                SubS[i].t=SubS[i].t_post;
            }else if(process.parallelisation->is_local_cpu()){
                SubS[i].t_post.resize(process.temps->nbpastemps+1);//utile pour le post traitement en incremental
                for(unsigned pt=1;pt<(nbpastemps+1);pt++){
                    SubS[i].t[pt].allocate(SubS[i]);
                }
                if (process.nom_calcul=="incr"){
                    for(unsigned pt=1;pt<(process.temps->nbpastemps+1);pt++){
                        SubS[i].t_post[pt].allocate(SubS[i]);
                    }
                }
            }
                    
        }
        if (process.parallelisation->is_local_cpu()){
            for(unsigned i=0;i<SubI.size();++i){
                for(unsigned j=0;j<SubI[i].side.size();++j){
                    SubI[i].side[j].t.resize(nbpastemps+1);
                    if(process.recopie_t_post==1)
                        SubI[i].side[j].t=SubI[i].side[j].t_post;
                    else{
                        SubI[i].side[j].t_post.resize(process.temps->nbpastemps+1); //utile pour le post traitement en incremental
                        if(process.nom_calcul=="incr" and process.parallelisation->is_local_cpu()){
                            for(unsigned pt=0;pt<(nbpastemps+1);pt++){
                                SubI[i].side[j].t[pt].allocations(SubI[i].side[j].nodeeq.size()*DIM);
                            }
                            for(unsigned pt=0;pt<(process.temps->nbpastemps+1);pt++){
                                SubI[i].side[j].t_post[pt].allocations(SubI[i].side[j].nodeeq.size()*DIM);                        
                            }
                        }
                    }
                }
                /// allocation des parametres de comportements dependant du temps (exemple endommagement)
                SubI[i].t.resize(nbpastemps+1);
                for(unsigned j=0;j<nbpastemps+1;j++){
                    SubI[i].t[j].allocate(SubI[i]);
                }
            }
        }
    }
            
}


void assign_quantities_current_to_old(PointedSubstructures &SubS, PointedInterfaces &SubI, Process &process){
    for(unsigned i=0;i<SubS.size();i++){
        SubS[i].t_post[process.temps->pt_cur]=SubS[i].t[1];
        SubS[i].t[0]=SubS[i].t[1];
    }
    for(unsigned i=0;i<SubI.size();i++){
        for(unsigned j=0;j<SubI[i].side.size();j++){
            SubI[i].side[j].t[0].F=SubI[i].side[j].t[1].F;
            SubI[i].side[j].t[0].Wp=SubI[i].side[j].t[1].Wp;
            SubI[i].side[j].t[0].Fchap=SubI[i].side[j].t[1].Fchap;
            SubI[i].side[j].t[0].Wpchap=SubI[i].side[j].t[1].Wpchap;
            SubI[i].side[j].t[0].W=SubI[i].side[j].t[1].W;
            SubI[i].side[j].t[0].Wchap=SubI[i].side[j].t[1].Wchap;
            SubI[i].side[j].t_post[process.temps->pt_cur]=SubI[i].side[j].t[1];
        }
    }
}


void assign_t_post(PointedSubstructures &SubS, PointedInterfaces &SubI, Process &process){
    
    for(unsigned i=0;i<SubS.size();i++){
        SubS[i].t_post[process.temps->pt_cur]=SubS[i].t[1];
        SubS[i].t[0]=SubS[i].t[1];
    }
    for(unsigned i=0;i<SubI.size();i++){
        for(unsigned j=0;j<SubI[i].side.size();j++)
            SubI[i].side[j].t_post[process.temps->pt_cur]=SubI[i].side[j].t[1];
    }
}


void recopie_old_from_new(PointedInterfaces &SubI,Process &process) {
    unsigned nbpastemps=0;
    if (process.nom_calcul=="incr") nbpastemps=1;
    if (process.nom_calcul=="latin") nbpastemps=process.temps->nbpastemps;
    for(unsigned pt=1;pt<=nbpastemps;pt++)
        for( unsigned q=0;q<SubI.size() ;q++ )
            for( unsigned data=0;data<SubI[q].side.size() ;data++ ) {
                SubI[q].side[data].t[pt].oldF=SubI[q].side[data].t[pt].F;
                SubI[q].side[data].t[pt].oldWp=SubI[q].side[data].t[pt].Wp;
                SubI[q].side[data].t[pt].oldW=SubI[q].side[data].t[pt].W;
            }
}


void recopie_old_from_new_post(PointedInterfaces &SubI,Process &process) {
    for( unsigned q=0;q<SubI.size() ;q++ )
        for( unsigned data=0;data<SubI[q].side.size() ;data++ ) {
            SubI[q].side[data].t[1].oldF=SubI[q].side[data].t_post[1].F;
            SubI[q].side[data].t[1].oldWp=SubI[q].side[data].t_post[1].Wp;
            SubI[q].side[data].t[1].oldW=SubI[q].side[data].t_post[1].W;
        }
}


/** Assigne les resultats stockes dans le Sst::Time dans la formulation et force le calcul sur les autres grandeurs
 * Utilisee pour le stockage des resultats dans les fichiers
 */
void rebuild_state(Sst &S,Sst::Time &t, Process &process){
    S.mesh.load();
    S.apply_behavior();
    process.Fvol->apply_on_sst(S);
    process.Tload->apply_on_sst(S);
    S.f->set_mesh(S.mesh.m);
    upload_q(S,t);
    if(S.f == S.pb.formulation_plasticity_isotropy_stat_Qstat){
        upload_epsilon_p(S,t);  /// pour obtenir sigma
        upload_p(S,t);          /// pour obtenir R_p
    }
    if(S.f == S.pb.formulation_elasticity_damageable_isotropy_stat_Qstat or S.f == S.pb.formulation_mesomodele){
        upload_d1(S,t);
    }
    if(S.f == S.pb.formulation_mesomodele){
        upload_d2(S,t);
        upload_df(S,t);
    }
    S.f->update_variables();
    S.f->call_after_solve();
}

