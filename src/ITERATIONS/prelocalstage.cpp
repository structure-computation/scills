#include "prelocalstage.h"
#include "../COMPUTE/DataUser.h"
#include "../DEFINITIONS/Process.h"
#include "../DEFINITIONS/Sst.h"
#include "../DEFINITIONS/Interface.h"
#include "../DEFINITIONS/Boundary.h"
#include "../UTILITAIRES/utilitaires.h"
#include "../../LMT/include/codegen/codegen.h"
#include "../../LMT/include/containers/basicops.h"

using namespace Codegen;


/// Calcul des valeurs numeriques des fonctions temporelles des CL au pas de temps courant
void calc_CL_time(Process &process, ScCLVec &CL, DataUser &data_user ) {
    /// Creation des symbols
    Ex t = symbol("t");
    std::vector<Ex> symbols;
    symbols.push_back(t);
    if(data_user.options.Multiresolution_on==1){
        ///ajout des variables de multiresolution aux symboles
        for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
            Sc2String var="V"; var<<i_par; ///nom de la variable de multiresolution
            symbols.push_back(var.c_str());
        }
    }
    
    unsigned i_step=process.temps->step_cur;
    unsigned tpas=process.temps->time_step[i_step].pt_cur;
    TYPEREEL ti=process.temps->current_time;
    /// Recuperation de l'instant precedent (pour le calcul des derivees)
    TYPEREEL told = ti - process.temps->time_step[i_step].dt;
    /// Evaluation des CL
    for(unsigned i_cl=0;i_cl<CL.size();++i_cl) {
        /// Recuperation des expression symboliques
        Ex expr,expr2;
        expr = read_ex(CL[i_cl].fcts_temporelles[i_step],symbols);
        expr2 = read_ex(CL[i_cl].fcts_temporelles[i_step],symbols);
        /// Construction de la table des (symbol, valeur)
        Ex::MapExNum var_temp,var_temp2;
        var_temp[symbols[0]]=ti;
        var_temp2[symbols[0]]=told;
        if(data_user.options.Multiresolution_on==1){
            /// Evaluation des variables de multiresolution
            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
                var_temp[symbols[1+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
                var_temp2[symbols[1+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
            }
        }
        /// Evaluation des fonctions temporelles
        for( unsigned dir=0;dir<DIM ;dir++ ){
            CL[i_cl].ft[dir]=(TYPEREEL)expr.subs_numerical(var_temp);
        }
        
        /// Cas particuliers ou il faut calculer une derivee
        if(CL[i_cl].comp == "depl" or CL[i_cl].comp == "depl_normal"){
            for( unsigned dir=0;dir<DIM ;dir++ ){
                CL[i_cl].ft[dir]-=(TYPEREEL)expr2.subs_numerical(var_temp2);
                CL[i_cl].ft[dir]/= ti-told;
            }
        }
    }
}


/// Calcul des valeurs numeriques des fonctions spatiales de la CL et assignation du produit par les fonctions temporelles dans le vecteur V associe au chargement defini (Fchap ou Wpchap)
void assign_CL_spatial_temporel(Vec<TYPEREEL> &V, Vec<ScPoint > &nodeeq, Boundary &CL, int i_step, DataUser &data_user) {
    std::vector<Ex> symbols;
    
#if DIM==2
    symbols.push_back("x");
    symbols.push_back("y");
#elif DIM==3
    symbols.push_back("x");
    symbols.push_back("y");
    symbols.push_back("z");
#endif
    if(data_user.options.Multiresolution_on==1){
        ///ajout des variables de multiresolution aux symboles
        for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
            Sc2String var="V"; var<<i_par; ///nom de la variable de multiresolution
            symbols.push_back(var.c_str());
        }
    }
    
    for(unsigned i=0;i<nodeeq.size();++i) {
        Vec<TYPEREEL, DIM> data;
        
        for(unsigned dir=0;dir<DIM;++dir) { /// boucle sur la dimension des vecteurs
            Ex expr;
            expr = read_ex(CL.fcts_spatiales[i_step][dir],symbols);
            Ex::MapExNum var;
            for(unsigned i_dir=0;i_dir<DIM;++i_dir)///boucle sur les inconnues possibles (dimension des vecteurs)
                var[symbols[i_dir]]= nodeeq[i][i_dir];
            if(data_user.options.Multiresolution_on==1){
                ///evaluation des variables de multiresolution
                for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                    var[symbols[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
            }
            data[dir] = (TYPEREEL)expr.subs_numerical(var);
            data[dir] *= CL.ft[dir];
        }
        V[range(i*DIM,(i+1)*DIM)]=data;
    }
}


/// Calcul des valeurs numeriques des fonctions spatiales de la CL en deplacement normal et assignation du produit par les fonctions temporelles dans le vecteur V associe (Wpchap)
void assign_CL_spatial_temporel_normale(Vec<TYPEREEL> &V, Vec<ScPoint > &nodeeq, Vec<TYPEREEL> &neqs, Boundary &CL, int i_step, DataUser &data_user) {
    std::vector<Ex> symbols;
    
#if DIM==2
    symbols.push_back("x");
    symbols.push_back("y");
#elif DIM==3
    symbols.push_back("x");
    symbols.push_back("y");
    symbols.push_back("z");
#endif
    if(data_user.options.Multiresolution_on==1){
        ///ajout des variables de multiresolution aux symboles
        for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
            Sc2String var="V"; var<<i_par; //nom de la variable de multiresolution
            symbols.push_back(var.c_str());
        }
    }
    
    for(unsigned i=0;i<nodeeq.size();++i) {
        TYPEREEL data;
        Vec<TYPEREEL,DIM> neq = neqs[range(i*DIM,(i+1)*DIM)];
        ///une seule expression dans le cas d'un deplacement normal
        Ex expr;
        expr = read_ex(CL.fcts_spatiales[i_step][0],symbols);
        Ex::MapExNum var;
        for(unsigned i_dir=0;i_dir<DIM;++i_dir)//boucle sur les inconnues possibles (dimension des vecteurs)
            var[symbols[i_dir]]= nodeeq[i][i_dir];
        if(data_user.options.Multiresolution_on==1){
            ///evaluation des variables de multiresolution
            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                var[symbols[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
        }
        data = (TYPEREEL)expr.subs_numerical(var);
        Vec<TYPEREEL,DIM> temp=V[range(i*DIM,(i+1)*DIM)];
        V[range(i*DIM,(i+1)*DIM)]=ProjT(temp,neq)+CL.ft[0]*data*neq;
        /*if(i == 0){
            std::cout << "Maj CL normale :" << std::endl;
            std::cout << "data : " << data << std::endl;
            std::cout << "Vecteur : " << V[range(0,3)] << std::endl;
        }//*/
    }
}


void initialise_CL_values_space_time(ScInterRef &Inter, ScCLVec &CL, Process &process, DataUser &data_user ){
    /// Donnees utiles
    for(unsigned i_inter = 0; i_inter < Inter.size(); i_inter++){
        /// Teste si la CL doit etre initialisee
        if(Inter[i_inter].type != "Ext" or (Inter[i_inter].comp != "depl" and Inter[i_inter].comp != "depl_normal"))
            continue;   /// Si inutile, passer a l'interface suivante
        const unsigned i_dir_max = (Inter[i_inter].comp == "depl_normal")? 1 : DIM;
        /// Creation des symbols
        Ex t = symbol("t");
        std::vector<Ex> symbols;
        symbols.push_back(t);
        if(data_user.options.Multiresolution_on==1){
            ///ajout des variables de multiresolution aux symboles
            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
                Sc2String var="V"; var<<i_par; ///nom de la variable de multiresolution
                symbols.push_back(var.c_str());
            }
        }
        
        /// Evaluation des CL
        for(unsigned i_cl=0;i_cl<CL.size();++i_cl) {
            /// Recuperation des expression symboliques
            Ex expr;
            expr = read_ex(CL[i_cl].fcts_temporelles[0],symbols);
            
            /// Construction de la table des (symbol, valeur)
            Ex::MapExNum var_temp;
            var_temp[symbols[0]]=0.0;
            if(data_user.options.Multiresolution_on==1){
                /// Evaluation des variables de multiresolution
                for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                    var_temp[symbols[1+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
            }
            /// Evaluation des fonctions temporelles
            std::cout << Inter[i_inter].comp << " : " << std::endl;
            for(unsigned i_dir=0;i_dir<i_dir_max ;i_dir++ ){
                CL[i_cl].ft[i_dir]=(TYPEREEL)expr.subs_numerical(var_temp);
                std::cout << CL[i_cl].id << ":    " << CL[i_cl].fcts_temporelles[0] << "*(" << CL[i_cl].fcts_spatiales[0][i_dir] << ")" << std::endl;
            }
        }
        
        std::vector<Ex> symbols2;
#if DIM==2
        symbols2.push_back("x");
        symbols2.push_back("y");
#elif DIM==3
        symbols2.push_back("x");
        symbols2.push_back("y");
        symbols2.push_back("z");
#endif
        if(data_user.options.Multiresolution_on==1){
            ///ajout des variables de multiresolution aux symboles
            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
                Sc2String var="V"; var<<i_par; //nom de la variable de multiresolution
                symbols2.push_back(var.c_str());
            }
        }
        
        for(unsigned i=0;i<Inter[i_inter].side[0].nodeeq.size();++i) {
            Vec<TYPEREEL,DIM> data;
            Vec<TYPEREEL,DIM> neq;
            if(Inter[i_inter].comp == "depl_normal"){
                neq = Inter[i_inter].side[0].neq[range(i*DIM,(i+1)*DIM)];
            }
            
            Vec<Ex,DIM> expr;
            for(unsigned i_dir = 0; i_dir < i_dir_max; i_dir++){
                expr[i_dir] = read_ex(CL[Inter[i_inter].refCL].fcts_spatiales[0][i_dir],symbols2);
            }
            
            Ex::MapExNum var;
            for(unsigned i_dir = 0; i_dir < DIM; ++i_dir)//boucle sur les inconnues possibles (dimension des vecteurs)
                var[symbols2[i_dir]]= Inter[i_inter].side[0].nodeeq[i][i_dir];
            if(data_user.options.Multiresolution_on==1){
                ///evaluation des variables de multiresolution
                for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                    var[symbols2[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
            }
            for(unsigned i_dir = 0; i_dir < i_dir_max; i_dir++){
                data[i_dir] = (TYPEREEL)expr[i_dir].subs_numerical(var);
            }
            if(Inter[i_inter].comp == "depl_normal"){
                Vec<TYPEREEL,DIM> temp=Inter[i_inter].side[0].t[0].W[range(i*DIM,(i+1)*DIM)];
                Inter[i_inter].side[0].t[0].W[range(i*DIM,(i+1)*DIM)]=ProjT(temp,neq)+CL[Inter[i_inter].refCL].ft[0]*data[0]*neq;
            }else if(Inter[i_inter].comp == "depl"){
                Inter[i_inter].side[0].t[0].W[range(i*DIM,(i+1)*DIM)]=data;
            }
        }
    }
}


/// Mise a jour des CL et assignation des nouvelles valeurs a Fchap et Wpchap (pour une iteration LATIN)
void assign_CL_values_space_time_latin(ScInterRef &Inter, ScCLVec &CL, Process &process, DataUser &data_user) {
    for(unsigned i_inter=0;i_inter<Inter.size();++i_inter) {
        if (Inter[i_inter].type=="Ext" and Inter[i_inter].comp != "periodique") {
            for(unsigned pt=1;pt<=process.temps->nbpastemps;pt++) {
                ///creation de la valeur de la fonction temporelle au pas de temps donne : modification de CL.ft
                process.temps->pt_cur=pt;
                calc_CL_time(process,CL,data_user);
                
                if (Inter[i_inter].comp=="effort") {
                    assign_CL_spatial_temporel(Inter[i_inter].side[0].t[pt].Fchap,Inter[i_inter].side[0].nodeeq,CL[Inter[i_inter].refCL],process.temps->step_cur,data_user);
                    //cout << "Eff " << Inter[i_inter].side[0].Fchap << endl;
                } else if (Inter[i_inter].comp=="depl" or Inter[i_inter].comp=="depl_nul" or Inter[i_inter].comp=="vit" or Inter[i_inter].comp=="vit_nulle") {
                    assign_CL_spatial_temporel(Inter[i_inter].side[0].t[pt].Wpchap,Inter[i_inter].side[0].nodeeq,CL[Inter[i_inter].refCL],process.temps->step_cur,data_user);
                    //cout << "vit " <<Inter[i_inter].side[0].Wchap << endl;
                } else if (Inter[i_inter].comp=="depl_normal" or Inter[i_inter].comp=="vit_normale") {
                    assign_CL_spatial_temporel_normale(Inter[i_inter].side[0].t[pt].Wpchap,Inter[i_inter].side[0].nodeeq,Inter[i_inter].side[0].neq,CL[Inter[i_inter].refCL],process.temps->step_cur,data_user);//il faut updater juste la nouvelle partie normal mais on laisse la partie tangentielle
                    if(process.reprise_calcul==0)
                        Inter[i_inter].side[0].t[pt].Fchap.set(0.0);
                } else if (Inter[i_inter].comp=="sym") {
                    if(process.reprise_calcul==0)
                        assign_CL_spatial_temporel(Inter[i_inter].side[0].t[pt].Wpchap,Inter[i_inter].side[0].nodeeq,CL[Inter[i_inter].refCL],process.temps->step_cur,data_user);//si on reprend le calcul les conditions sont toujours ok sinon on initialise à 0
                        if(process.reprise_calcul==0)
                            Inter[i_inter].side[0].t[pt].Fchap.set(0.0);
                } else {
                    std::cout << "Erreur d'interface ext - prelocalstage " << std::endl;
                    assert(0);
                }
            }
        } else if(Inter[i_inter].comp=="Contact_jeu" or Inter[i_inter].comp=="Contact_jeu_physique") {
            //le jeu est reparti en moyenne sur chacun des deplacements des cotes 1 et 2
            Inter[i_inter].side[1].t[0].Wchap=Inter[i_inter].param_comp->jeu[Inter[i_inter].side[0].ddlcorresp]/2.;
            Inter[i_inter].side[0].t[0].Wchap=-1.*Inter[i_inter].param_comp->jeu/2.;
        }
    }
}


/// Mise a jour des CL et assignation des nouvelles valeurs a Fchap et Wpchap (pour une iteration incrementale)
//*
void assign_CL_values_space_time_incr(ScInterRef &Inter, ScCLVec &CL, Process &process, DataUser &data_user ) {
    std::cout << "Mise a jour des Conditions aux Limites sur le processeur " << process.parallelisation->rank << " : " << std::endl;
    for(unsigned i_inter=0;i_inter<Inter.size();++i_inter) {
        std::cout << "    id : " << Inter[i_inter].id  << "    type : " << Inter[i_inter].type << "    comportement : " << Inter[i_inter].comp << "    valeur : ";
        if (Inter[i_inter].type=="Ext" and Inter[i_inter].comp != "periodique") {
            calc_CL_time(process,CL,data_user);
            if (Inter[i_inter].comp=="effort") {
                assign_CL_spatial_temporel(Inter[i_inter].side[0].t[1].Fchap,Inter[i_inter].side[0].nodeeq,CL[Inter[i_inter].refCL],process.temps->step_cur,data_user);
                std::cout << Inter[i_inter].side[0].t[1].Fchap[0] << " ; ";
                std::cout << Inter[i_inter].side[0].t[1].Fchap[1] << " ; ";
                std::cout << Inter[i_inter].side[0].t[1].Fchap[2] << " ; ";
            } else if (Inter[i_inter].comp=="effort_normal") {
                assign_CL_spatial_temporel_normale(Inter[i_inter].side[0].t[1].Fchap,Inter[i_inter].side[0].nodeeq,Inter[i_inter].side[0].neq,CL[Inter[i_inter].refCL],process.temps->step_cur,data_user);
                std::cout << Inter[i_inter].side[0].t[1].Fchap[0] << " ; ";
                std::cout << Inter[i_inter].side[0].t[1].Fchap[1] << " ; ";
                std::cout << Inter[i_inter].side[0].t[1].Fchap[2] << " ; ";
            } else if (Inter[i_inter].comp=="depl" or Inter[i_inter].comp=="depl_nul" or Inter[i_inter].comp=="vit" or Inter[i_inter].comp=="vit_nulle") {
                assign_CL_spatial_temporel(Inter[i_inter].side[0].t[1].Wpchap,Inter[i_inter].side[0].nodeeq,CL[Inter[i_inter].refCL],process.temps->step_cur,data_user);
                std::cout << Inter[i_inter].side[0].t[1].Wpchap[0] << " ; ";
                std::cout << Inter[i_inter].side[0].t[1].Wpchap[1] << " ; ";
                std::cout << Inter[i_inter].side[0].t[1].Wpchap[2] << " ; ";
            } else if (Inter[i_inter].comp=="sym") {
                if(process.temps->pt_cur==1) {
                    if(process.reprise_calcul==0){
                        assign_CL_spatial_temporel(Inter[i_inter].side[0].t[1].Wpchap,Inter[i_inter].side[0].nodeeq,CL[Inter[i_inter].refCL],process.temps->step_cur,data_user);
                    }
                    if(process.reprise_calcul==0)
                        Inter[i_inter].side[0].t[1].Fchap.set(0.0);
                }
                std::cout << Inter[i_inter].side[0].t[1].Wpchap[0] << " ; ";
                std::cout << Inter[i_inter].side[0].t[1].Wpchap[1] << " ; ";
                std::cout << Inter[i_inter].side[0].t[1].Wpchap[2] << " ; ";
            } else if (Inter[i_inter].comp=="depl_normal" or Inter[i_inter].comp=="vit_normale") {
                if(process.temps->pt_cur==1) {
                    assign_CL_spatial_temporel_normale(Inter[i_inter].side[0].t[1].Wpchap,Inter[i_inter].side[0].nodeeq,Inter[i_inter].side[0].neq,CL[Inter[i_inter].refCL],process.temps->step_cur,data_user);//le Wpchap evolue au cours des iterations donc si on reprend on initialise avec le resultat du calcul precedent donc on fait rien...
                    if(process.reprise_calcul==0)
                        Inter[i_inter].side[0].t[1].Fchap.set(0.0);
                } else {
                    Vec<TYPEREEL> Wpchapnormal;
                    Wpchapnormal.resize(Inter[i_inter].side[0].t[1].Wpchap.size());
                    assign_CL_spatial_temporel_normale(Wpchapnormal,Inter[i_inter].side[0].nodeeq,Inter[i_inter].side[0].neq,CL[Inter[i_inter].refCL],process.temps->step_cur,data_user);
                    Inter[i_inter].side[0].t[1].Wpchap = Inter[i_inter].side[0].Pt(Inter[i_inter].side[0].t[1].Wpchap)+Wpchapnormal;
                }
                std::cout << Inter[i_inter].side[0].t[1].Wpchap[0] << " ; ";
                std::cout << Inter[i_inter].side[0].t[1].Wpchap[1] << " ; ";
                std::cout << Inter[i_inter].side[0].t[1].Wpchap[2] << " ; ";
            } else {
                std::cout << "Erreur d'interface ext - prelocalstage " << std::endl;
                assert(0);
            }
        } else if(Inter[i_inter].comp=="Contact_jeu" or Inter[i_inter].comp=="Contact_jeu_physique") {
            ///le jeu est reparti en moyenne sur chacun des deplacements des cotes 1 et 2 et impose uniquement pour le premier pas de temps
            if(process.temps->pt_cur==1) {
                Inter[i_inter].side[1].t[0].Wchap[Inter[i_inter].side[1].ddlcorresp]=Inter[i_inter].param_comp->jeu/2.;
                Inter[i_inter].side[0].t[0].Wchap=-1.*Inter[i_inter].param_comp->jeu/2.;
                //if (Inter[i_inter].num == 15 ) std::cout << "Jeu : " << Inter[i_inter].side[0].t[0].Wchap << endl;
                //if (Inter[i_inter].num == 15 ) std::cout << "Jeu : " << Inter[i_inter].side[1].t[0].Wchap << endl;
            }
        }
        std::cout << std::endl;
    }
}//*/
/*
/// Mise a jour des CL et assignation des nouvelles valeurs a Fchap et Wpchap (pour une iteration incrementale)
void assign_CL_values_space_time_incr(Vec<VecPointedValues<Interface> > &Inter, Vec<Boundary> &CL, Process &process, DataUser &data_user ) {
    /// Construit la liste des symboles (si necessaire, c.f. Boundary.cpp)
    Boundary::buildBoundarySymbols(data_user);
    /// Associe les valeurs aux grandeurs dependant du temps
    Boundary::updateTimeValues(process,data_user);
    /// Evalue la CL de chaque interface
    for(unsigned i_inter=0;i_inter<Inter.size();++i_inter) {
        std::cout << "id : " << Inter[i_inter].id << "    comportement : " << Inter[i_inter].comp << std::endl;
        if (Inter[i_inter].type=="Ext" and Inter[i_inter].comp != "periodique") {
            if (Inter[i_inter].comp=="effort") {
                CL[Inter[i_inter].refCL].evaluate(process.temps->step_cur,Inter[i_inter].side[0].t[1].Fchap,Inter[i_inter].side[0].nodeeq);
            } else if (Inter[i_inter].comp=="effort_normal") {
                CL[Inter[i_inter].refCL].evaluate(process.temps->step_cur,Inter[i_inter].side[0].t[1].Fchap,Inter[i_inter].side[0].nodeeq,Inter[i_inter].side[0].neq);
            } else if (Inter[i_inter].comp=="vit" or Inter[i_inter].comp=="vit_nulle") {
                CL[Inter[i_inter].refCL].evaluate(process.temps->step_cur,Inter[i_inter].side[0].t[1].Wpchap,Inter[i_inter].side[0].nodeeq);
            } else if (Inter[i_inter].comp=="sym") {
                if(process.temps->pt_cur==1) {
                    if(process.reprise_calcul==0)
                        CL[Inter[i_inter].refCL].evaluate(process.temps->step_cur,Inter[i_inter].side[0].t[1].Wpchap,Inter[i_inter].side[0].nodeeq);
                    if(process.reprise_calcul==0)
                        Inter[i_inter].side[0].t[1].Fchap.set(0.0);
                }
            } else if (Inter[i_inter].comp=="vit_normale") {
                if(process.temps->pt_cur==1) {
                    CL[Inter[i_inter].refCL].evaluate(process.temps->step_cur,Inter[i_inter].side[0].t[1].Wpchap,Inter[i_inter].side[0].nodeeq,Inter[i_inter].side[0].neq);//le Wpchap evolue au cours des iterations donc si on reprend on initialise avec le resultat du calcul precedent donc on fait rien...
                    if(process.reprise_calcul==0)
                        Inter[i_inter].side[0].t[1].Fchap.set(0.0);
                } else {
                    Vec<TYPEREEL> Wpchapnormal;
                    Wpchapnormal.resize(Inter[i_inter].side[0].t[1].Wpchap.size());
                    CL[Inter[i_inter].refCL].evaluate(process.temps->step_cur,Wpchapnormal,Inter[i_inter].side[0].nodeeq,Inter[i_inter].side[0].neq);
                    Inter[i_inter].side[0].t[1].Wpchap = Inter[i_inter].side[0].Pt(Inter[i_inter].side[0].t[1].Wpchap)+Wpchapnormal;
                }
            } else {
                std::cout << "Erreur d'interface ext - prelocalstage " << std::endl;
                assert(0);
            }
        } else if(Inter[i_inter].comp=="Contact_jeu" or Inter[i_inter].comp=="Contact_jeu_physique") {
            ///le jeu est reparti en moyenne sur chacun des deplacements des cotes 1 et 2 et impose uniquement pour le premier pas de temps
            if(process.temps->pt_cur==1) {
                Inter[i_inter].side[1].t[0].Wchap[Inter[i_inter].side[1].ddlcorresp]=Inter[i_inter].param_comp->jeu/2.;
                Inter[i_inter].side[0].t[0].Wchap=-1.*Inter[i_inter].param_comp->jeu/2.;
                //if (Inter[i_inter].num == 15 ) std::cout << "Jeu : " << Inter[i_inter].side[0].t[0].Wchap << endl;
                //if (Inter[i_inter].num == 15 ) std::cout << "Jeu : " << Inter[i_inter].side[1].t[0].Wchap << endl;
            }
        }
    }
}//*/
