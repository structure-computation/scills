#include "../DEFINITIONS/TimeParameters.h"
#include "../DEFINITIONS/LatinParameters.h"
#include "manipulate_quantities.h"

using namespace LMT;


//allocation des differentes quantites
/** \ingroup   Strategie_iterative
 * \brief Allocation de mémoire pour les sous-structures, interfaces et vecteurs macro. 
 * 
 * Pour chaque pas de temps les différents vecteurs sont initialisés à 0. On n'alloue pas la mémoire pour les différents vecteurs en chaque pas de temps des ssts car ils ne sont stockés qu'à convergence.
 */
void allocate_quantities(Vec<VecPointedValues<Sst > > &S, Vec<VecPointedValues<Interface > > &Inter,Process &process,MacroProblem &Global){
    /// Recuperation du nombre de pas de temps
    unsigned nbpastemps;
    if(process.nom_calcul=="incr")
        nbpastemps=1;
    else if (process.nom_calcul=="latin")
        nbpastemps=process.temps->nbpastemps;
    else{std::cout << "Nom de clacul non implemente dans allocate.h" << endl;assert(0);}
    
    /// Allocation des quantites des Sst
    if (process.latin->save_depl_SST==1){
        for(unsigned i=0;i<S.size();++i){
            S[i].t.resize(nbpastemps+1);
            if(process.recopie_t_post==1)
                S[i].t=S[i].t_post;
            for(unsigned pt=0;pt<(nbpastemps+1);pt++){  ///L'initialisation commence a 0 et non plus a 1, pour les initialisations
                unsigned nbddl = S[i].mesh.node_list_size*DIM;
                unsigned nbelem = S[i].mesh.elem_list_size;
                S[i].t[pt].allocations(nbddl,nbelem,process,S[i].plastique);
            }
        }
    }
    
    /// Allocation des quantites des Interface
    if (process.rank>0 or process.size==1){
        for(unsigned i=0;i<Inter.size();++i){
            for(unsigned j=0;j<Inter[i].side.size();++j){
                Inter[i].side[j].t.resize(nbpastemps+1);
                if(process.recopie_t_post==1)
                    Inter[i].side[j].t=Inter[i].side[j].t_post;
                else{
                    for(unsigned pt=0;pt<(nbpastemps+1);pt++)
                        Inter[i].side[j].t[pt].allocations(Inter[i].side[j].nodeeq.size()*DIM,process);
                }
            }
            ///allocation des parametres de comportements dependant du temps (exemple endommagement)
            Inter[i].param_comp->t.resize(nbpastemps+1);
            for(unsigned j=0;j<Inter[i].param_comp->t.size();j++)
                Inter[i].param_comp->t[j].allocate(Inter[i].side[0].nodeeq.size());
        }
    }
    
    if(process.recopie_t_post!=1){   
        Global.allocations(process);
        process.latin->error.resize(process.latin->nbitermax+1);
        process.latin->error.set(0.0);
    }
    
}


//allocation de t_post uniquement pour les sorties de resultat
void allocate_quantities_post(Vec<VecPointedValues<Sst > > &S, Vec<VecPointedValues<Interface > > &Inter,Process &process){
    unsigned nbpastemps;
    if(process.nom_calcul=="incr")
        nbpastemps=1;
    else if (process.nom_calcul=="latin")
        nbpastemps=process.temps->nbpastemps;
    else{std::cout << "Nom de clacul non implemente dans allocate.h" << std::endl;assert(0);}
    
    if (process.latin->save_depl_SST==1)
        for(unsigned i=0;i<S.size();++i){
            S[i].t.resize(nbpastemps+1);
            if(process.recopie_t_post==1)
                S[i].t=S[i].t_post;
            else
                S[i].t_post.resize(process.temps->nbpastemps+1);//utile pour le post traitement en incremental
                for(unsigned pt=1;pt<(nbpastemps+1);pt++){
                    unsigned nbddl = S[i].mesh.node_list_size*DIM;
                    unsigned nbelem = S[i].mesh.elem_list_size;
                    S[i].t[pt].allocations(nbddl,nbelem,process,S[i].plastique);
                }
                if (process.nom_calcul=="incr")
                    for(unsigned pt=1;pt<(process.temps->nbpastemps+1);pt++){
                        unsigned nbddl = S[i].mesh.node_list_size*DIM;
                        unsigned nbelem = S[i].mesh.elem_list_size;
                        S[i].t_post[pt].allocations(nbddl,nbelem,process,S[i].plastique);
                    }
                    
        }
        if (process.rank>0 or process.size==1)
            for(unsigned i=0;i<Inter.size();++i){
                for(unsigned j=0;j<Inter[i].side.size();++j){
                    Inter[i].side[j].t.resize(nbpastemps+1);
                    if(process.recopie_t_post==1)
                        Inter[i].side[j].t=Inter[i].side[j].t_post;
                    else{
                        Inter[i].side[j].t_post.resize(process.temps->nbpastemps+1); //utile pour le post traitement en incremental
                        for(unsigned pt=0;pt<(nbpastemps+1);pt++)
                            Inter[i].side[j].t[pt].allocations(Inter[i].side[j].nodeeq.size()*DIM,process);
                        for(unsigned pt=0;pt<(process.temps->nbpastemps+1);pt++)
                            if(process.nom_calcul=="incr")
                                Inter[i].side[j].t_post[pt].allocations(Inter[i].side[j].nodeeq.size()*DIM,process);                        
                            
                    }
                }
                //allocation des parametres de comportements dependant du temps (exemple endommagement)
                Inter[i].param_comp->t.resize(nbpastemps+1);
                for(unsigned j=0;j<Inter[i].param_comp->t.size();j++)
                    Inter[i].param_comp->t[j].allocate(Inter[i].side[0].nodeeq.size());
            }
            
}


void assign_quantities_current_to_old(Vec<VecPointedValues<Sst> > &S, Vec<VecPointedValues<Interface> > &Inter, Process &process){
    
    for(unsigned i=0;i<S.size();i++){
        S[i].t_post[process.temps->pt_cur]=S[i].t[1];
    }
    for(unsigned i=0;i<Inter.size();i++){
        for(unsigned j=0;j<Inter[i].side.size();j++){
            Inter[i].side[j].t[0].F=Inter[i].side[j].t[1].F;
            Inter[i].side[j].t[0].Wp=Inter[i].side[j].t[1].Wp;
            Inter[i].side[j].t[0].Fchap=Inter[i].side[j].t[1].Fchap;
            Inter[i].side[j].t[0].Wpchap=Inter[i].side[j].t[1].Wpchap;
            Inter[i].side[j].t[0].W=Inter[i].side[j].t[1].W;
            Inter[i].side[j].t[0].Wchap=Inter[i].side[j].t[1].Wchap;
            Inter[i].side[j].t_post[process.temps->pt_cur]=Inter[i].side[j].t[1];
        }
    }
}


void assign_t_post(Vec<VecPointedValues<Interface> > &Inter, Process &process){
    
    for(unsigned i=0;i<Inter.size();i++){
        for(unsigned j=0;j<Inter[i].side.size();j++)
            Inter[i].side[j].t_post[process.temps->pt_cur]=Inter[i].side[j].t[1];
    }
}


void recopie_old_from_new(Vec<Interface> &Inter,Process &process) {
    unsigned nbpastemps=0;
    if (process.nom_calcul=="incr") nbpastemps=1;
    if (process.nom_calcul=="latin") nbpastemps=process.temps->nbpastemps;
    for(unsigned pt=1;pt<=nbpastemps;pt++)
        for( unsigned q=0;q<Inter.size() ;q++ )
            for( unsigned data=0;data<Inter[q].side.size() ;data++ ) {
                Inter[q].side[data].t[pt].oldF=Inter[q].side[data].t[pt].F;
                Inter[q].side[data].t[pt].oldWp=Inter[q].side[data].t[pt].Wp;
                Inter[q].side[data].t[pt].oldW=Inter[q].side[data].t[pt].W;
            }
}


void recopie_old_from_new_post(Vec<Interface> &Inter,Process &process) {
    for( unsigned q=0;q<Inter.size() ;q++ )
        for( unsigned data=0;data<Inter[q].side.size() ;data++ ) {
            Inter[q].side[data].t[1].oldF=Inter[q].side[data].t_post[1].F;
            Inter[q].side[data].t[1].oldWp=Inter[q].side[data].t_post[1].Wp;
            Inter[q].side[data].t[1].oldW=Inter[q].side[data].t_post[1].W;
        }
}


void assign_dep_cont_slave(Sst &S,Vec<TYPEREEL> &q, DataUser &data_user){
    S.mesh.load();
    S.assign_material_on_element(data_user);
    S.f->set_mesh(S.mesh.m);
    S.f->get_result()=q;
    S.f->update_variables();
    S.f->call_after_solve();
}

