#include "containers/vecpointedvalues.h"

/** \ingroup  Post_Traitements
\brief Extraction de l'evolution du déplacement d'un point donné par ses coordonnées 
*/
template<class TT>
void extraction_depl_pt(Vec<VecPointedValues<Sst<3,TT> > > &S, Param &process) {
    //recherche du point dans une sous-structure
    for(unsigned i=0;i<S.size();i++) {
        for(unsigned j=0;j<S[i].mesh->node_list.size();j++) {
            if(length(S[i].mesh->node_list[j].pos-process.affichage->coor_point)<=0.001) {
                Vec<unsigned> repddl=range(j*3,(j+1)*3);
                Vec<double > evol_depl_pt_x,evol_depl_pt_y,evol_depl_pt_z;

                evol_depl_pt_x.resize(process.temps->nbpastemps+1,0.);
                evol_depl_pt_y.resize(process.temps->nbpastemps+1,0.);
                evol_depl_pt_z.resize(process.temps->nbpastemps+1,0.);
                //S[i].t.resize(process.temps->nbpastemps+1);
                if(process.nom_calcul=="incr") {
                    for(unsigned imic=1;imic<=process.temps->nbpastemps;imic++) {
                        evol_depl_pt_x[imic]=(S[i].t_post[imic].q[repddl[0]]);
                        evol_depl_pt_y[imic]=(S[i].t_post[imic].q[repddl[1]]);
                        evol_depl_pt_z[imic]=(S[i].t_post[imic].q[repddl[2]]);
                    }
                } else if(process.nom_calcul=="latin") {
                    for(unsigned imic=1;imic<=process.temps->nbpastemps;imic++) {
                        evol_depl_pt_x[imic]=(S[i].t[imic].q[repddl[0]]);
                        evol_depl_pt_y[imic]=(S[i].t[imic].q[repddl[1]]);
                        evol_depl_pt_z[imic]=(S[i].t[imic].q[repddl[2]]);
                    }
                }
                std::cout << "evol x" << evol_depl_pt_x << endl;
                std::cout << "evol y" << evol_depl_pt_y << endl;
                std::cout << "evol z" << evol_depl_pt_z << endl;
                if (process.affichage->command_file==""){
                  GnuPlot gp1;
                  gp1.plot(evol_depl_pt_x);
                  gp1.wait();
                  GnuPlot gp2;
                  gp2.plot(evol_depl_pt_y);
                  gp2.wait();
                  GnuPlot gp3;
                  gp3.plot(evol_depl_pt_z);
                  gp3.wait();
                  std::cout  << endl;
                }
            }
        }
    }
}

/** \ingroup  Post_Traitements
\brief Extraction de l'evolution des champs sur une interface a partir de son numéro et de l'element souhaite
*/
/*template<class TT>
void affichage_var_inter(Vec<VecPointedValues<Sst<2,TT> > > &S,Vec<Interface<2,TT>, Param &process) {*/
template<class TV1,class TV2>
void affichage_var_inter(TV1 &S,TV2 &Inter, Param &process) {
    typedef typename TV2::template SubType<0>::T INTER;
    //recherche du point dans une sous-structure
    for(unsigned i=0;i<S.size();i++) {
        for(unsigned j=0;j<S[i].edge.size();j++) {
            if (S[i].edge[j].internum==process.affichage->param_ener[0]) {
                unsigned internum=S[i].edge[j].internum;
                unsigned datanum=S[i].edge[j].datanum;
                unsigned elemnum=process.affichage->param_ener[1];
                
                if (datanum==0){
                    Vec<double> evolFchap,evolWchap,evolWpchap;
                    Vec<Vec<double> > evolFchap1,evolWchap1,evolWpchap1;
                evolFchap.resize(process.temps->nbpastemps+1,0.);
                evolWchap.resize(process.temps->nbpastemps+1,0.);
                evolWpchap.resize(process.temps->nbpastemps+1,0.);
                evolFchap1.resize(process.temps->nbpastemps+1);
                evolWchap1.resize(process.temps->nbpastemps+1);
                evolWpchap1.resize(process.temps->nbpastemps+1);
                Vec<unsigned> rep=range(elemnum*INTER::dim,(elemnum+1)*INTER::dim);
                if(process.nom_calcul=="incr") {
                    for( unsigned imic=0;imic<process.temps->nbpastemps+1 ;imic++ ){
                        Vec<double> tmp;
                        tmp=Inter[internum].side[datanum].Pt(Inter[internum].side[datanum].t_post[imic].Fchap);
                        evolFchap[imic]=norm_2(tmp[rep]);
                        tmp=Inter[internum].side[datanum].Pt(Inter[internum].side[datanum].t_post[imic].Wchap)-Inter[internum].side[1-datanum].Pt(Inter[internum].side[1-datanum].t_post[imic].Wchap);
                        evolWchap[imic]=norm_2(tmp[rep]);
                        tmp=Inter[internum].side[datanum].Pt(Inter[internum].side[datanum].t_post[imic].Wpchap)-Inter[internum].side[1-datanum].Pt(Inter[internum].side[1-datanum].t_post[imic].Wpchap);
                        evolWpchap[imic]=norm_2(tmp[rep]);
                    }
                }else if(process.nom_calcul=="latin") {
                    for( unsigned imic=0;imic<process.temps->nbpastemps+1 ;imic++ ){
                        Vec<double> tmp;
                        tmp=Inter[internum].side[datanum].t[imic].Fchap;
                        evolFchap1[imic].resize(INTER::dim,0.0);
                        evolFchap1[imic]=tmp[rep];
                        tmp=Inter[internum].side[datanum].Pt(tmp);
                        evolFchap[imic]=norm_2(tmp[rep]);
                        tmp=Inter[internum].side[datanum].t[imic].Wchap-Inter[internum].side[1-datanum].t[imic].Wchap;
                        evolWchap1[imic].resize(INTER::dim,0.0);
                        evolWchap1[imic]=tmp[rep];
                        tmp=Inter[internum].side[datanum].Pt(Inter[internum].side[datanum].t[imic].Wchap)-Inter[internum].side[1-datanum].Pt(Inter[internum].side[1-datanum].t[imic].Wchap);
                        evolWchap[imic]=norm_2(tmp[rep]);
                        tmp=Inter[internum].side[datanum].t[imic].Wpchap-Inter[internum].side[1-datanum].t[imic].Wpchap;
                        evolWpchap1[imic].resize(INTER::dim,0.0);
                        evolWpchap1[imic]=tmp[rep];
                        tmp=Inter[internum].side[datanum].Pt(Inter[internum].side[datanum].t[imic].Wpchap)-Inter[internum].side[1-datanum].Pt(Inter[internum].side[1-datanum].t[imic].Wpchap);
                        evolWpchap[imic]=norm_2(tmp[rep]);
                    }
                }


                if (process.affichage->param_ener[2] == 0 or process.affichage->param_ener[2] == 3) std::cout << "Cote : " << datanum << " , surface=" << Inter[internum].side[datanum].mesh->elem_list[elemnum]->measure_virtual() << " : " << evolFchap << endl;
                if (process.affichage->param_ener[2] == 1 or process.affichage->param_ener[2] == 3) std::cout << "Cote : " << datanum << " : " << evolWchap << endl;
                if (process.affichage->param_ener[2] == 2 ) std::cout << "Cote : " << datanum << " : " << evolWpchap << endl;

                if (process.affichage->command_file==""){
                  GnuPlot gp1;
                  if (process.affichage->param_ener[2] == 0) gp1.plot(evolFchap);
                  if (process.affichage->param_ener[2] == 1) gp1.plot(evolWchap);
                  if (process.affichage->param_ener[2] == 2) gp1.plot(evolWpchap);
                  if (process.affichage->param_ener[2] == 3) gp1.plot(evolWchap,evolFchap);
                  gp1.wait();
                  std::cout  << endl;
                }
                if (process.affichage->param_ener[2] == 3){
                    std::ofstream os("interfaces_valeur.m",std::ofstream::app);
                    os << "%Interface numero : "<<internum<<", element : " <<elemnum<<", surface : " << Inter[internum].side[datanum].mesh->elem_list[elemnum]->measure_virtual() <<endl;
                    os << "Fchap=["<<evolFchap1<<"];"<<endl;
                    os << "Wchap=["<<evolWchap1<<"];"<<endl;
                    os << "Wpchap=["<<evolWpchap1<<"];"<<endl;
                    os << "FchapTnorme=["<<evolFchap<<"];"<<endl;
                    os << "WchapTnorme=["<<evolWchap<<"];"<<endl;
                    os << "WpchapTnorme=["<<evolWpchap<<"];"<<endl;
                    os << "figure;plot(WchapTnorme,FchapTnorme);"<<endl;
                    os << "figure;plot(WpchapTnorme,FchapTnorme);"<<endl;
                    os.close();
                }
                }
            }
        }
    }
}



/** \ingroup  Post_Traitements
\brief Extraction de l'evolution du déplacement d'un point donné par ses coordonnées 
*/
template<class TT>
void extraction_depl_pt(Vec<VecPointedValues<Sst<2,TT> > > &S, Param &process) {
    //recherche du point dans une sous-structure
    for(unsigned i=0;i<S.size();i++) {
        for(unsigned j=0;j<S[i].mesh->node_list.size();j++) {
            if(length(S[i].mesh->node_list[j].pos-process.affichage->coor_point)<=0.001) {
                Vec<unsigned> repddl=range(j*2,(j+1)*2);
                Vec<double > evol_depl_pt_x,evol_depl_pt_y;

                evol_depl_pt_x.resize(process.temps->nbpastemps+1,0.);
                evol_depl_pt_y.resize(process.temps->nbpastemps+1,0.);
                //S[i].t.resize(process.temps->nbpastemps+1);
                if(process.nom_calcul=="incr") {
                    for(unsigned imic=1;imic<=process.temps->nbpastemps;imic++) {
                        evol_depl_pt_x[imic]=(S[i].t_post[imic].q[repddl[0]]);
                        evol_depl_pt_y[imic]=(S[i].t_post[imic].q[repddl[1]]);
                    }
                } else if(process.nom_calcul=="latin") {
                    for(unsigned imic=1;imic<=process.temps->nbpastemps;imic++) {
                        evol_depl_pt_x[imic]=(S[i].t[imic].q[repddl[0]]);
                        evol_depl_pt_y[imic]=(S[i].t[imic].q[repddl[1]]);
                    }
                }
                std::cout << "evol x" << evol_depl_pt_x << endl;
                std::cout << "evol y" << evol_depl_pt_y << endl;
                if (process.affichage->command_file==""){
                  GnuPlot gp1;
                  gp1.plot(evol_depl_pt_x);
                  gp1.wait();
                  GnuPlot gp2;
                  gp2.plot(evol_depl_pt_y);
                  gp2.wait();
                  std::cout  << endl;
                }
            }
        }
    }
}
