#include "codegen/codegen.h"
#include "containers/basicops.h"
using namespace LMT;
using namespace std;

using namespace Codegen;

//assignation de valeurs variables en fonction des variables x, y et z
template<class TV, class TVN, class BOUNDARY>
void assign_CL_spatial(TV &V, TVN &nodeeq, BOUNDARY &CL) {

    typedef typename BOUNDARY::T T;
    std::vector<Ex> symbols;
    if (BOUNDARY::dim==2) {
        symbols.push_back("x");
        symbols.push_back("y");
    } else if (BOUNDARY::dim==3) {
        symbols.push_back("x");
        symbols.push_back("y");
        symbols.push_back("z");
    }

    for(unsigned i=0;i<nodeeq.size();++i) {
        Vec<T, BOUNDARY::dim> data;
        for(unsigned d1=0;d1<BOUNDARY::dim;++d1) { // boucle sur la dimension des vecteurs
            Ex expr;
            expr = read_ex(CL.fcts_spatiales[d1].c_str(),symbols);
            Ex::MapExNum var;
            for(unsigned d2=0;d2<BOUNDARY::dim;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                var[symbols[d2]]= nodeeq[i][d2];
            data[d1] = (T)expr.subs_numerical(var);
        }
        V[range(i*BOUNDARY::dim,(i+1)*BOUNDARY::dim)]=data;
    }
}

//************************************************************************************************************
//Procedures faisant intervenir le temps : modulation des fonctions spatiales par une fonction scalaire temporelle
//************************************************************************************************************

template<class BOUNDARY>
void calc_CL_time(Param &process,Vec<BOUNDARY> &CL) {
    typedef typename BOUNDARY::T T;
    Ex t = symbol("t");
    unsigned i_step=process.temps->step_cur;
    unsigned tpas=process.temps->time_step[i_step].pt_cur;    
    std::vector<Ex> symbols;
    symbols.push_back(t);
    double ti=process.temps->time_step[i_step].t_ini+(tpas+1)*process.temps->time_step[i_step].dt;
    string fcttemps;
    Ex res;
/*    std::cout << i_step << " " << tpas << " " << ti << endl;*/
    for(unsigned j=0;j<CL.size();++j) {
        fcttemps = CL[j].fcts_temporelles[i_step];
/*        std::cout << fcttemps << endl;*/
//         for(unsigned i=0;i<CL[j].fcts_temporelles.size();++i) {
//             if(ti>=CL[j].intervalles_temps[i][0] && ti<=CL[j].intervalles_temps[i][1]) {
//                 fcttemps = CL[j].fcts_temporelles[i];
//                 break;
//             }
//         }
        ////modif DAVID 02-09-2007
        CL[j].ft.resize(BOUNDARY::dim);
        Ex expr;
        expr = read_ex(fcttemps.c_str(),symbols);
        for( unsigned d1=0;d1<BOUNDARY::dim ;d1++ ){
                CL[j].ft[d1]=(T)expr.subs_numerical(t,ti);
        }    
  /*      CL[j].ft.resize(CL[j].fcts_temporelles[i_step].size());
        for(unsigned d1=0;d1<CL[j].fcts_temporelles[i_step].size();++d1) { // boucle sur la dimension des vecteurs
            Ex expr;
            expr = read_ex(fcttemps[d1].c_str(),symbols);
            CL[j].ft[d1]= (T)expr.subs_numerical(t,ti);
        }
        if (CL[j].ft.size() != BOUNDARY::dim) {
            CL[j].ft.resize(BOUNDARY::dim) ;
            for( unsigned d1=1;d1<BOUNDARY::dim ;d1++ ){
                CL[j].ft[d1]=CL[j].ft[0];
            }
        }*/
        //res = read_ex(fcttemps.c_str(),symbols);
         //CL[j].ft=res.subs_numerical(t,ti);
         // if (CL[j].comp == "effort") std::cout << CL[j].ft << endl;
         ////fin modif DAVID

    }
}

//assignation de valeurs variables en fonction des variables x, y et z et du temps t
template<class TV, class TVN, class BOUNDARY>
void assign_CL_spatial_temporel(TV &V, TVN &nodeeq, BOUNDARY &CL, int i_step) {

    typedef typename BOUNDARY::T T;
    std::vector<Ex> symbols;

    if (BOUNDARY::dim==2) {
        symbols.push_back("x");
        symbols.push_back("y");
    } else if (BOUNDARY::dim==3) {
        symbols.push_back("x");
        symbols.push_back("y");
        symbols.push_back("z");
    }

    for(unsigned i=0;i<nodeeq.size();++i) {
        Vec<T, BOUNDARY::dim> data;

        for(unsigned d1=0;d1<BOUNDARY::dim;++d1) { // boucle sur la dimension des vecteurs
            Ex expr;
            expr = read_ex(CL.fcts_spatiales[i_step][d1].c_str(),symbols);
            Ex::MapExNum var;
            for(unsigned d2=0;d2<BOUNDARY::dim;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                var[symbols[d2]]= nodeeq[i][d2];
            data[d1] = (T)expr.subs_numerical(var);
            ////modif DAVID 02-09-2007
            data[d1] *= CL.ft[d1];
            //data[d1] *= CL.ft;
            //// fin modif
        }
        V[range(i*BOUNDARY::dim,(i+1)*BOUNDARY::dim)]=data;
    }
}

template<class TV, class TVN, class TV2, class BOUNDARY>
void assign_CL_spatial_temporel_normale(TV &V, TVN &nodeeq, TV2 &neqs, BOUNDARY &CL, int i_step) {

    typedef typename BOUNDARY::T T;
    typedef Vec<T,BOUNDARY::dim> Pvec;
    std::vector<Ex> symbols;
    Pvec temp;

    if (BOUNDARY::dim==2) {
        symbols.push_back("x");
        symbols.push_back("y");
    } else if (BOUNDARY::dim==3) {
        symbols.push_back("x");
        symbols.push_back("y");
        symbols.push_back("z");
    }

    for(unsigned i=0;i<nodeeq.size();++i) {
        T data;
        Pvec neq = neqs[range(i*BOUNDARY::dim,(i+1)*BOUNDARY::dim)];
        //une seule expression dans le cas d'un deplacement normal
        Ex expr;
        expr = read_ex(CL.fcts_spatiales[i_step][0].c_str(),symbols);
        Ex::MapExNum var;
        for(unsigned d2=0;d2<BOUNDARY::dim;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
            var[symbols[d2]]= nodeeq[i][d2];
        data = (T)expr.subs_numerical(var);
        temp=V[range(i*BOUNDARY::dim,(i+1)*BOUNDARY::dim)];
         ////modif DAVID 02-09-2007
         V[range(i*BOUNDARY::dim,(i+1)*BOUNDARY::dim)]=ProjT(temp,neq)+CL.ft[0]*data*neq;

         //V[range(i*BOUNDARY::dim,(i+1)*BOUNDARY::dim)]=ProjT(temp,neq)+CL.ft*data*neq;
         //// fin modif
    }

}


/** \ingroup   LATIN
\brief Assignation pour tous les pas de temps des conditions limites aux quantit�s chapeau des interfaces de bord
 
Pour les interfaces � effort impos�, on utilise assign_CL_eff().
Pour les interfaces � d�placement impos�, on utilise assign_CL_depl().
Pour les interfaces de type sym�trie, on utilise assign_CL_sym()
Pour les interfaces de type d�placement normal donn�, on utilise assign_CL_depl_normal()
Pour les interfaces int�rieures de type contact avec jeu, on utilise assign_jeu()
Sinon une erreur est envoy�e.
 
*/
template<class TV2,class TV5>
void assign_CL_values_space_time_latin(TV2 &Inter, TV5 &CL, Param &process) {
    for(unsigned q=0;q<Inter.size();++q) {
        if (Inter[q].type=="Ext" and Inter[q].comp != "periodique") {
            for(unsigned pt=1;pt<=process.temps->nbpastemps;pt++) {
                //creation de la valeur de la fonction temporelle au pas de temps donne : modification de CL.ft
                process.temps->pt_cur=pt;
                calc_CL_time(process,CL);

                if (Inter[q].comp=="effort") {
                    assign_CL_spatial_temporel(Inter[q].side[0].t[pt].Fchap,Inter[q].side[0].nodeeq,CL[Inter[q].refCL],process.temps->step_cur);
                    //cout << "Eff " << Inter[q].side[0].Fchap << endl;
                } else if (Inter[q].comp=="depl") {
                    assign_CL_spatial_temporel(Inter[q].side[0].t[pt].Wpchap,Inter[q].side[0].nodeeq,CL[Inter[q].refCL],process.temps->step_cur);
                    //cout << "depl " <<Inter[q].side[0].Wchap << endl;
                } else if (Inter[q].comp=="sym") {
                    if(process.reprise_calcul==0)
                        assign_CL_spatial_temporel(Inter[q].side[0].t[pt].Wpchap,Inter[q].side[0].nodeeq,CL[Inter[q].refCL],process.temps->step_cur);//si on reprend le calcul les conditions sont toujours ok sinon on initialise � 0
                    if(process.reprise_calcul==0)
                        Inter[q].side[0].t[pt].Fchap.set(0.0);
                } else if (Inter[q].comp=="depl_normal") {
                    assign_CL_spatial_temporel_normale(Inter[q].side[0].t[pt].Wpchap,Inter[q].side[0].nodeeq,Inter[q].side[0].neq,CL[Inter[q].refCL],process.temps->step_cur);//il faut updater juste la nouvelle partie normal mais on laisse la partie tangentielle
                    if(process.reprise_calcul==0)
                        Inter[q].side[0].t[pt].Fchap.set(0.0);
                } else {
                    std::cout << "Erreur d'interface ext - prelocalstage " <<endl;
                    assert(0);
                }
            }
        } else if(Inter[q].comp=="Contact_jeu" or Inter[q].comp=="Contact_jeu_physique") {
            //le jeu est reparti en moyenne sur chacun des deplacements des cotes 1 et 2
                Inter[q].side[1].t[0].Wchap=Inter[q].param_comp->jeu[Inter[q].side[0].ddlcorresp]/2.;
                Inter[q].side[0].t[0].Wchap=-1.*Inter[q].param_comp->jeu/2.;
        }
    }
};

/** \ingroup   Incrementale
\brief Assignation pour chaque pas de temps des conditions limites aux quantit�s chapeau des interfaces de bord
 
Pour les interfaces � effort impos�, on utilise assign_CL_eff().
Pour les interfaces � d�placement impos�, on utilise assign_CL_depl().
Pour les interfaces de type sym�trie, on utilise assign_CL_sym()
Pour les interfaces de type d�placement normal donn�, on utilise assign_CL_depl_normal()
Pour les interfaces int�rieures de type contact avec jeu, on utilise assign_jeu() uniquement pour le piquet de temps initial
Sinon une erreur est envoy�e.
 
*/
template<class TV2,class TV5>
void assign_CL_values_space_time_incr(TV2 &Inter, TV5 &CL, Param &process) {
    for(unsigned q=0;q<Inter.size();++q) {
        std::cout << "Inter[q].id = " << Inter[q].id << std::endl;
        std::cout << "Inter[q].comp = " << Inter[q].comp << std::endl;
        if (Inter[q].type=="Ext" and Inter[q].comp != "periodique") {
            calc_CL_time(process,CL);
            std::cout << "calc_CL_time ok " << std::endl;
            if (Inter[q].comp=="effort") {
                assign_CL_spatial_temporel(Inter[q].side[0].t[1].Fchap,Inter[q].side[0].nodeeq,CL[Inter[q].refCL],process.temps->step_cur);
            } else if (Inter[q].comp=="effort_normal") {
                assign_CL_spatial_temporel_normale(Inter[q].side[0].t[1].Fchap,Inter[q].side[0].nodeeq,Inter[q].side[0].neq,CL[Inter[q].refCL],process.temps->step_cur);
            } else if (Inter[q].comp=="depl" or Inter[q].comp=="depl_nul") {
                assign_CL_spatial_temporel(Inter[q].side[0].t[1].Wpchap,Inter[q].side[0].nodeeq,CL[Inter[q].refCL],process.temps->step_cur);
            } else if (Inter[q].comp=="sym") {
                if(process.temps->pt_cur==1) {
                    if(process.reprise_calcul==0)
                        assign_CL_spatial_temporel(Inter[q].side[0].t[1].Wpchap,Inter[q].side[0].nodeeq,CL[Inter[q].refCL],process.temps->step_cur);
                    if(process.reprise_calcul==0)
                        Inter[q].side[0].t[1].Fchap.set(0.0);
                }
            } else if (Inter[q].comp=="depl_normal") {
                if(process.temps->pt_cur==1) {
                    assign_CL_spatial_temporel_normale(Inter[q].side[0].t[1].Wpchap,Inter[q].side[0].nodeeq,Inter[q].side[0].neq,CL[Inter[q].refCL],process.temps->step_cur);//le Wpchap evolue au cours des iterations donc si on reprend on initialise avec le resultat du calcul precedent donc on fait rien...
                    if(process.reprise_calcul==0)
                        Inter[q].side[0].t[1].Fchap.set(0.0);
                } else {
                    Vec<typename TV2::template SubType<0>::T::T> Wpchapnormal;
                    Wpchapnormal.resize(Inter[q].side[0].t[1].Wpchap.size());
                    assign_CL_spatial_temporel_normale(Wpchapnormal,Inter[q].side[0].nodeeq,Inter[q].side[0].neq,CL[Inter[q].refCL],process.temps->step_cur);
                    Inter[q].side[0].t[1].Wpchap = Inter[q].side[0].Pt(Inter[q].side[0].t[1].Wpchap)+Wpchapnormal;
                }
            } else {
                std::cout << "Erreur d'interface ext - prelocalstage " <<endl;
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

};


/** \ingroup   Latin
\brief Pour la reprise d'un calcul, on recopie les donnees relues dans les quantit�s old
 
 */
template<class TV2>
void recopie_old_from_new(TV2 &Inter,Param &process) {
  unsigned nbpastemps=0;
  if (process.nom_calcul=="incr") nbpastemps=1;
  if (process.nom_calcul=="latin") nbpastemps=process.temps->nbpastemps;
  for(unsigned pt=1;pt<=nbpastemps;pt++) {
        for( unsigned q=0;q<Inter.size() ;q++ ) {
            for( unsigned data=0;data<Inter[q].side.size() ;data++ ) {
                Inter[q].side[data].t[pt].oldF=Inter[q].side[data].t[pt].F;
                Inter[q].side[data].t[pt].oldWp=Inter[q].side[data].t[pt].Wp;
                Inter[q].side[data].t[pt].oldW=Inter[q].side[data].t[pt].W;
            }
        }
    }
}

/** \ingroup   Latin
\brief Pour la reprise d'un calcul, on recopie les donnees relues dans les quantit�s old
 
 */
template<class TV2>
void recopie_old_from_new_post(TV2 &Inter,Param &process) {
        for( unsigned q=0;q<Inter.size() ;q++ ) {
            for( unsigned data=0;data<Inter[q].side.size() ;data++ ) {
              Inter[q].side[data].t[1].oldF=Inter[q].side[data].t_post[1].F;
              Inter[q].side[data].t[1].oldWp=Inter[q].side[data].t_post[1].Wp;
              Inter[q].side[data].t[1].oldW=Inter[q].side[data].t_post[1].W;
            }
        }
}

