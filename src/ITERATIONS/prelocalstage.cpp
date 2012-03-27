#include "prelocalstage.h"
#include "../DEFINITIONS/Interface.h"
#include "../UTILITAIRES/utilitaires.h"

void calc_CL_time(Process &process,Vec<Boundary> &CL, DataUser &data_user ) {
    Ex t = symbol("t");
    unsigned i_step=process.temps->step_cur;
    unsigned tpas=process.temps->time_step[i_step].pt_cur;    
    std::vector<Ex> symbols;
    symbols.push_back(t);
    if(data_user.options.Multiresolution_on==1){
        //ajout des variables de multiresolution aux symboles
        for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
            Sc2String var="V"; var<<i_par; //nom de la variable de multiresolution
            symbols.push_back(var.c_str());
        }
    }
    
    double ti=process.temps->time_step[i_step].t_ini+(tpas+1)*process.temps->time_step[i_step].dt;
    Sc2String fcttemps;
    Ex res;
    for(unsigned j=0;j<CL.size();++j) {
        fcttemps = CL[j].fcts_temporelles[i_step];
        CL[j].ft.resize(DIM);
        Ex expr;
        expr = read_ex(fcttemps.c_str(),symbols);
        Ex::MapExNum var_temp;
        var_temp[symbols[0]]=ti;
        if(data_user.options.Multiresolution_on==1){
            //evaluation des variables de multiresolution
            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                var_temp[symbols[1+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
        }
        
        for( unsigned d1=0;d1<DIM ;d1++ ){
            //CL[j].ft[d1]=(T)expr.subs_numerical(t,ti);
            CL[j].ft[d1]=(TYPEREEL)expr.subs_numerical(var_temp);
        }
    }
}


void assign_CL_spatial_temporel(Vec<TYPEREEL> &V, Vec<Vec<TYPEREEL,DIM> > &nodeeq, Boundary &CL, int i_step, DataUser &data_user) {
    
    std::vector<Ex> symbols;
    
    if (DIM==2) {
        symbols.push_back("x");
        symbols.push_back("y");
    } else if (DIM==3) {
        symbols.push_back("x");
        symbols.push_back("y");
        symbols.push_back("z");
    }
    if(data_user.options.Multiresolution_on==1){
        //ajout des variables de multiresolution aux symboles
        for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
            Sc2String var="V"; var<<i_par; //nom de la variable de multiresolution
            symbols.push_back(var.c_str());
        }
    }
    
    for(unsigned i=0;i<nodeeq.size();++i) {
        Vec<TYPEREEL, DIM> data;
        
        for(unsigned d1=0;d1<DIM;++d1) { // boucle sur la dimension des vecteurs
            Ex expr;
            expr = read_ex(CL.fcts_spatiales[i_step][d1].c_str(),symbols);
            Ex::MapExNum var;
            for(unsigned d2=0;d2<DIM;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                var[symbols[d2]]= nodeeq[i][d2];
            if(data_user.options.Multiresolution_on==1){
                //evaluation des variables de multiresolution
                for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                    var[symbols[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
            }
            data[d1] = (TYPEREEL)expr.subs_numerical(var);
            data[d1] *= CL.ft[d1];
        }
        V[range(i*DIM,(i+1)*DIM)]=data;
    }
}


void assign_CL_spatial_temporel_normale(Vec<TYPEREEL> &V, Vec<Vec<TYPEREEL,DIM> > &nodeeq, Vec<TYPEREEL> &neqs, Boundary &CL, int i_step, DataUser &data_user) {
    
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
        //ajout des variables de multiresolution aux symboles
        for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
            Sc2String var="V"; var<<i_par; //nom de la variable de multiresolution
            symbols.push_back(var.c_str());
        }
    }
    
    for(unsigned i=0;i<nodeeq.size();++i) {
        TYPEREEL data;
        Vec<TYPEREEL,DIM> neq = neqs[range(i*DIM,(i+1)*DIM)];
        //une seule expression dans le cas d'un deplacement normal
        Ex expr;
        expr = read_ex(CL.fcts_spatiales[i_step][0].c_str(),symbols);
        Ex::MapExNum var;
        for(unsigned d2=0;d2<DIM;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
            var[symbols[d2]]= nodeeq[i][d2];
        if(data_user.options.Multiresolution_on==1){
            //evaluation des variables de multiresolution
            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                var[symbols[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
        }
        data = (TYPEREEL)expr.subs_numerical(var);
        Vec<TYPEREEL,DIM> temp=V[range(i*DIM,(i+1)*DIM)];
        V[range(i*DIM,(i+1)*DIM)]=ProjT(temp,neq)+CL.ft[0]*data*neq;
    }
}


void assign_CL_values_space_time_latin(Vec<VecPointedValues<Interface> > &Inter, Vec<Boundary> &CL, Process &process, DataUser &data_user) {
    for(unsigned q=0;q<Inter.size();++q) {
        if (Inter[q].type=="Ext" and Inter[q].comp != "periodique") {
            for(unsigned pt=1;pt<=process.temps->nbpastemps;pt++) {
                //creation de la valeur de la fonction temporelle au pas de temps donne : modification de CL.ft
                process.temps->pt_cur=pt;
                calc_CL_time(process,CL,data_user);
                
                if (Inter[q].comp=="effort") {
                    assign_CL_spatial_temporel(Inter[q].side[0].t[pt].Fchap,Inter[q].side[0].nodeeq,CL[Inter[q].refCL],process.temps->step_cur,data_user);
                    //cout << "Eff " << Inter[q].side[0].Fchap << endl;
                } else if (Inter[q].comp=="depl") {
                    assign_CL_spatial_temporel(Inter[q].side[0].t[pt].Wpchap,Inter[q].side[0].nodeeq,CL[Inter[q].refCL],process.temps->step_cur,data_user);
                    //cout << "depl " <<Inter[q].side[0].Wchap << endl;
                } else if (Inter[q].comp=="sym") {
                    if(process.reprise_calcul==0)
                        assign_CL_spatial_temporel(Inter[q].side[0].t[pt].Wpchap,Inter[q].side[0].nodeeq,CL[Inter[q].refCL],process.temps->step_cur,data_user);//si on reprend le calcul les conditions sont toujours ok sinon on initialise à 0
                        if(process.reprise_calcul==0)
                            Inter[q].side[0].t[pt].Fchap.set(0.0);
                } else if (Inter[q].comp=="depl_normal") {
                    assign_CL_spatial_temporel_normale(Inter[q].side[0].t[pt].Wpchap,Inter[q].side[0].nodeeq,Inter[q].side[0].neq,CL[Inter[q].refCL],process.temps->step_cur,data_user);//il faut updater juste la nouvelle partie normal mais on laisse la partie tangentielle
                    if(process.reprise_calcul==0)
                        Inter[q].side[0].t[pt].Fchap.set(0.0);
                } else {
                    std::cout << "Erreur d'interface ext - prelocalstage " << std::endl;
                    assert(0);
                }
            }
        } else if(Inter[q].comp=="Contact_jeu" or Inter[q].comp=="Contact_jeu_physique") {
            //le jeu est reparti en moyenne sur chacun des deplacements des cotes 1 et 2
            Inter[q].side[1].t[0].Wchap=Inter[q].param_comp->jeu[Inter[q].side[0].ddlcorresp]/2.;
            Inter[q].side[0].t[0].Wchap=-1.*Inter[q].param_comp->jeu/2.;
        }
    }
}


void assign_CL_values_space_time_incr(Vec<VecPointedValues<Interface> > &Inter, Vec<Boundary> &CL, Process &process, DataUser &data_user ) {
    for(unsigned q=0;q<Inter.size();++q) {
        std::cout << "Inter[q].id = " << Inter[q].id << std::endl;
        std::cout << "Inter[q].comp = " << Inter[q].comp << std::endl;
        if (Inter[q].type=="Ext" and Inter[q].comp != "periodique") {
            calc_CL_time(process,CL,data_user);
            if (Inter[q].comp=="effort") {
                assign_CL_spatial_temporel(Inter[q].side[0].t[1].Fchap,Inter[q].side[0].nodeeq,CL[Inter[q].refCL],process.temps->step_cur,data_user);
            } else if (Inter[q].comp=="effort_normal") {
                assign_CL_spatial_temporel_normale(Inter[q].side[0].t[1].Fchap,Inter[q].side[0].nodeeq,Inter[q].side[0].neq,CL[Inter[q].refCL],process.temps->step_cur,data_user);
            } else if (Inter[q].comp=="depl" or Inter[q].comp=="depl_nul") {
                assign_CL_spatial_temporel(Inter[q].side[0].t[1].Wpchap,Inter[q].side[0].nodeeq,CL[Inter[q].refCL],process.temps->step_cur,data_user);
            } else if (Inter[q].comp=="sym") {
                if(process.temps->pt_cur==1) {
                    if(process.reprise_calcul==0)
                        assign_CL_spatial_temporel(Inter[q].side[0].t[1].Wpchap,Inter[q].side[0].nodeeq,CL[Inter[q].refCL],process.temps->step_cur,data_user);
                    if(process.reprise_calcul==0)
                        Inter[q].side[0].t[1].Fchap.set(0.0);
                }
            } else if (Inter[q].comp=="depl_normal") {
                if(process.temps->pt_cur==1) {
                    assign_CL_spatial_temporel_normale(Inter[q].side[0].t[1].Wpchap,Inter[q].side[0].nodeeq,Inter[q].side[0].neq,CL[Inter[q].refCL],process.temps->step_cur,data_user);//le Wpchap evolue au cours des iterations donc si on reprend on initialise avec le resultat du calcul precedent donc on fait rien...
                    if(process.reprise_calcul==0)
                        Inter[q].side[0].t[1].Fchap.set(0.0);
                } else {
                    Vec<TYPEREEL> Wpchapnormal;
                    Wpchapnormal.resize(Inter[q].side[0].t[1].Wpchap.size());
                    assign_CL_spatial_temporel_normale(Wpchapnormal,Inter[q].side[0].nodeeq,Inter[q].side[0].neq,CL[Inter[q].refCL],process.temps->step_cur,data_user);
                    Inter[q].side[0].t[1].Wpchap = Inter[q].side[0].Pt(Inter[q].side[0].t[1].Wpchap)+Wpchapnormal;
                }
            } else {
                std::cout << "Erreur d'interface ext - prelocalstage " << std::endl;
                assert(0);
            }
        } else if(Inter[q].comp=="Contact_jeu" or Inter[q].comp=="Contact_jeu_physique") {
            //le jeu est reparti en moyenne sur chacun des deplacements des cotes 1 et 2 et impose uniquement pour le premier pas de temps
            if(process.temps->pt_cur==1) {
                Inter[q].side[1].t[0].Wchap[Inter[q].side[1].ddlcorresp]=Inter[q].param_comp->jeu/2.;
                Inter[q].side[0].t[0].Wchap=-1.*Inter[q].param_comp->jeu/2.;
                //                 if (Inter[q].num == 15 ) std::cout << "Jeu : " << Inter[q].side[0].t[0].Wchap << endl;
                //                 if (Inter[q].num == 15 ) std::cout << "Jeu : " << Inter[q].side[1].t[0].Wchap << endl;
            }
        }
    }
}
