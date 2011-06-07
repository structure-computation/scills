//
// C++ Implementation: calculs_energies
//
// Description:
//
//
// Author: Alain CAIGNOT <caignot@lmt.ens-cachan.fr>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
//librairie Hugo
#include "containers/mat.h"

#include "mesh/mesh.h"
#include <fstream>
#include <map>

// fichiers de definition des variables
#include "definition_PARAM_MICRO_INTER.h"
#include "definition_PARAM.h"
#include "definition_PARAM_AFFICHAGE.h"
#include "definition_SST_time.h"
#include "definition_INTER_time.h"


/**
Fonction permettant le calcul de l'énergie dissipée à partir de la structure Inter et du Temps avec les quantités chapeaux.
*/
template <class TV1,class TV2>
void calcul_ener_dissi_chap(TV1 &S, TV2 &Inter,Vec<double> &dissipation,Param &process) {
    dissipation.set(0.);
    Vec<double> dissi_inter;
    dissi_inter.resize(dissipation.size());
    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned i=S[j].edge[e].internum;
            unsigned data=S[j].edge[e].datanum;
            if ((Inter[i].comp=="Contact" or Inter[i].comp=="Contact_jeu" or Inter[i].comp=="Contact_jeu_physique" or Inter[i].comp=="Contact_ep") and data==0) {
                dissi_inter.set(0.);
                if(process.nom_calcul=="incr") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        dissi_inter[j+1]=dissi_inter[j]+process.temps->dt*(
                                             dot((Inter[i].side[0].Pt(Inter[i].side[0].t_post[j+1].Fchap)+Inter[i].side[0].Pt(Inter[i].side[0].t_post[j].Fchap))/2.0,Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t_post[j+1].Wpchap))+
                                             dot((Inter[i].side[1].Pt(Inter[i].side[1].t_post[j+1].Fchap)+Inter[i].side[1].Pt(Inter[i].side[1].t_post[j].Fchap))/2.0,Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t_post[j+1].Wpchap)));
                    }
                } else if(process.nom_calcul=="latin") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        dissi_inter[j+1]=dissi_inter[j]+process.temps->dt*(
                                             dot((Inter[i].side[0].Pt(Inter[i].side[0].t[j+1].Fchap)+Inter[i].side[0].Pt(Inter[i].side[0].t[j].Fchap))/2.0,Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t[j+1].Wpchap))+
                                             dot((Inter[i].side[1].Pt(Inter[i].side[1].t[j+1].Fchap)+Inter[i].side[1].Pt(Inter[i].side[1].t[j].Fchap))/2.0,Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t[j+1].Wpchap)));
                    }
                } else {
                    std::cout << "Nom de calcul nom pris en compte" << endl;
                    assert(0);
                }
                dissipation+=dissi_inter;
                //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter <<  endl;
                
                if (norm_2(Inter[i].side[0].ddlcorresp-Inter[i].side[1].ddlcorresp)>1e-6) {
                    std::cout << "ATTENTION Interface : " << Inter[i].num << " maillage numerote differement : " << endl; 
                    std::cout << Inter[i].side[0].ddlcorresp << endl; 
                    std::cout  << Inter[i].side[1].ddlcorresp<< endl;
                }

            }
        }
    }
}
template <class TV1,class TV2>
        void calcul_ener_dissi_chap_petrus(TV1 &S, TV2 &Inter,Vec<double> &dissipation,Param &process) {
            dissipation.set(0.);
            Vec<double> dissi_inter;
            dissi_inter.resize(dissipation.size());
            for(unsigned j=0;j<S.size();j++) {
                for(unsigned e=0;e<S[j].edge.size();++e) {
                    unsigned i=S[j].edge[e].internum;
                    unsigned data=S[j].edge[e].datanum;
                    if ((Inter[i].comp=="Contact" or Inter[i].comp=="Contact_jeu" or Inter[i].comp=="Contact_jeu_physique" or Inter[i].comp=="Contact_ep") and data==0) {
                        dissi_inter.set(0.);
                        if(process.nom_calcul=="incr") {
                            for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                                dissi_inter[j+1]=dissi_inter[j]+process.temps->dt/2*(
                                        dot(Inter[i].side[1].t_post[j+1].Fchap,Inter[i].side[1].M*(Inter[i].side[1].t_post[j+1].Wpchap-Inter[i].side[0].t_post[j+1].Wpchap))+
                                        dot(Inter[i].side[1].t_post[j].Fchap  ,Inter[i].side[1].M*(Inter[i].side[1].t_post[j].Wpchap  -Inter[i].side[0].t_post[j].Wpchap)));
                            }
                        } else if(process.nom_calcul=="latin") {
                            for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                                dissi_inter[j+1]=dissi_inter[j]+process.temps->dt/2*(
                                        dot(Inter[i].side[1].t[j+1].Fchap,Inter[i].side[1].M*(Inter[i].side[1].t[j+1].Wpchap-Inter[i].side[0].t[j+1].Wpchap))+
                                        dot(Inter[i].side[1].t[j].Fchap  ,Inter[i].side[1].M*(Inter[i].side[1].t[j].Wpchap  -Inter[i].side[0].t[j].Wpchap)));
                            }
                        } else {
                            std::cout << "Nom de calcul nom pris en compte" << endl;
                            assert(0);
                        }
                        dissipation+=dissi_inter;
                        //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter <<  endl;
                    }
                }
            }
        }
        
template <class TV1,class TV2>
void calcul_ener_dissi_chap2(TV1 &S, TV2 &Inter,Vec<double> &dissipation,Param &process) {
    typedef typename TV2::template SubType<0>::T INTER;
    dissipation.set(0.);
    Vec<double> dissi_inter;
    dissi_inter.resize(dissipation.size());
    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned internum=S[j].edge[e].internum;
            unsigned datanum=S[j].edge[e].datanum;
            if ((Inter[internum].comp=="Contact" or Inter[internum].comp=="Contact_jeu" or Inter[internum].comp=="Contact_jeu_physique" or Inter[internum].comp=="Contact_ep") and datanum==0) {
                dissi_inter.set(0.);
                if(process.nom_calcul=="incr") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        dissi_inter[j+1]=dissi_inter[j]+process.temps->dt/2*(
                                dot(Inter[internum].side[1].t_post[j+1].Fchap,Inter[internum].side[1].M*(Inter[internum].side[1].t_post[j+1].Wpchap-Inter[internum].side[0].t_post[j+1].Wpchap))+
                                dot(Inter[internum].side[1].t_post[j].Fchap  ,Inter[internum].side[1].M*(Inter[internum].side[1].t_post[j].Wpchap  -Inter[internum].side[0].t_post[j].Wpchap)));
                    }
                } else if(process.nom_calcul=="latin") {
                    Vec<double> Fchap1,Fchap2,Wpchap1,Wpchap2;
                    Fchap1.resize(Inter[internum].side[datanum].t[0].Fchap.size());
                    Fchap2.resize(Inter[internum].side[datanum].t[0].Fchap.size());
                    Wpchap1.resize(Inter[internum].side[datanum].t[0].Fchap.size());
                    Wpchap2.resize(Inter[internum].side[datanum].t[0].Fchap.size());
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        double dissi_tmp=0;
                        Fchap1=Inter[internum].side[datanum].Pt(Inter[internum].side[datanum].t[j+1].Fchap);
                        Fchap2=Inter[internum].side[datanum].Pt(Inter[internum].side[datanum].t[j].Fchap);
                        Wpchap1=Inter[internum].side[datanum].Pt(Inter[internum].side[datanum].t[j+1].Wpchap);
                        Wpchap2=Inter[internum].side[1-datanum].Pt(Inter[internum].side[1-datanum].t[j+1].Wpchap);
                        for( unsigned elem=0;elem<Wpchap2.size()/INTER::dim ;elem++ ){
                            Vec<unsigned> rep=range(elem*INTER::dim,(elem+1)*INTER::dim);
                            dissi_tmp+=(norm_2(Fchap1[rep])+norm_2(Fchap2[rep]))/2.0*Inter[internum].side[datanum].mesh->elem_list[elem]->measure_virtual()*norm_2(Wpchap1[rep]-Wpchap2[rep]);

                        }
                        dissi_inter[j+1]=dissi_inter[j]+process.temps->dt*(dissi_tmp);
                    }
                } else {
                    std::cout << "Nom de calcul nom pris en compte" << endl;
                    assert(0);
                }
                dissipation+=dissi_inter;
                //std::cout << "Contribution interface " << Inter[internum].num << " entre les pieces " << Inter[internum].vois << " : " << dissi_inter <<  endl;
            }
        }
    }
}

/**
Fonction permettant le calcul de l'énergie dissipée à partir de la structure Inter et du Temps avec les quantités n.
*/
template <class TV1, class TI>
void calcul_ener_dissi_lin(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
    dissipation.set(0.);
    Vec<double> dissi_inter;
    dissi_inter.resize(dissipation.size());
    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned i=S[j].edge[e].internum;
            unsigned data=S[j].edge[e].datanum;
            if ((Inter[i].comp=="Contact" or Inter[i].comp=="Contact_jeu" or Inter[i].comp=="Contact_jeu_physique" or Inter[i].comp=="Contact_ep") and data==0) {
                dissi_inter.set(0.);
                if(process.nom_calcul=="incr") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        dissi_inter[j+1]=dissi_inter[j]+process.temps->dt*(
                                             dot((Inter[i].side[0].Pt(Inter[i].side[0].t_post[j+1].F)+Inter[i].side[0].Pt(Inter[i].side[0].t_post[j].F))/2.0,Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t_post[j+1].Wp))+
                                             dot((Inter[i].side[1].Pt(Inter[i].side[1].t_post[j+1].F)+Inter[i].side[1].Pt(Inter[i].side[1].t_post[j].F))/2.0,Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t_post[j+1].Wp)));
                    }
                } else if(process.nom_calcul=="latin") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        dissi_inter[j+1]=dissi_inter[j]+process.temps->dt*(
                                             dot((Inter[i].side[0].Pt(Inter[i].side[0].t[j+1].F)+Inter[i].side[0].Pt(Inter[i].side[0].t[j].F))/2.0,Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t[j+1].Wp))+
                                             dot((Inter[i].side[1].Pt(Inter[i].side[1].t[j+1].F)+Inter[i].side[1].Pt(Inter[i].side[1].t[j].F))/2.0,Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t[j+1].Wp)));
                    }
                } else {
                    std::cout << "Nom de calcul nom pris en compte" << endl;
                    assert(0);
                }
                dissipation+=dissi_inter;
                //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter << endl;
            }
        }
    }
}
template <class TV1, class TI>
        void calcul_ener_dissi_lin2(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
            dissipation.set(0.);
            Vec<double> dissi_inter;
            dissi_inter.resize(dissipation.size());
            for(unsigned j=0;j<S.size();j++) {
                for(unsigned e=0;e<S[j].edge.size();++e) {
                    unsigned i=S[j].edge[e].internum;
                    unsigned data=S[j].edge[e].datanum;
                    if ((Inter[i].comp=="Contact" or Inter[i].comp=="Contact_jeu" or Inter[i].comp=="Contact_jeu_physique" or Inter[i].comp=="Contact_ep") and data==0) {
                        dissi_inter.set(0.);
                        if(process.nom_calcul=="incr") {
                            for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                                dissi_inter[j+1]=dissi_inter[j]+process.temps->dt/2*(
                                        dot(Inter[i].side[0].t_post[j+1].F,Inter[i].side[0].M*Inter[i].side[0].t_post[j+1].Wp)+
                                        dot(Inter[i].side[1].t_post[j+1].F,Inter[i].side[1].M*Inter[i].side[1].t_post[j+1].Wp)+
                                        dot(Inter[i].side[0].t_post[j].F,Inter[i].side[0].M*Inter[i].side[0].t_post[j].Wp)+
                                        dot(Inter[i].side[1].t_post[j].F,Inter[i].side[1].M*Inter[i].side[1].t_post[j].Wp));
                            }
                        } else if(process.nom_calcul=="latin") {
                            for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                                dissi_inter[j+1]=dissi_inter[j]+process.temps->dt/2*(
                                        dot(Inter[i].side[0].t[j+1].F,Inter[i].side[0].M*Inter[i].side[0].t[j+1].Wp)+
                                        dot(Inter[i].side[1].t[j+1].F,Inter[i].side[1].M*Inter[i].side[1].t[j+1].Wp)+
                                        dot(Inter[i].side[0].t[j].F,Inter[i].side[0].M*Inter[i].side[0].t[j].Wp)+
                                        dot(Inter[i].side[1].t[j].F,Inter[i].side[1].M*Inter[i].side[1].t[j].Wp));
                            }
                        } else {
                            std::cout << "Nom de calcul nom pris en compte" << endl;
                            assert(0);
                        }
                        dissipation+=dissi_inter;
                        //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter << endl;
                    }
                }
            }
        }

/**
Fonction permettant le calcul de l'énergie imposée à partir de la structure Inter et du Temps avec les quantités chapeaux.
*/
template <class TV1, class TI>
void calcul_ener_imp_chap(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
    dissipation.set(0.);
    Vec<double> dissi_inter;
    dissi_inter.resize(dissipation.size());
    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned i=S[j].edge[e].internum;
            unsigned data=S[j].edge[e].datanum;
            if ((Inter[i].comp=="depl" or Inter[i].comp=="depl_normal" or Inter[i].comp=="effort"  or Inter[i].comp=="effort_normal"/*or Inter[i].comp=="Jeu_impose"*/) and data == 0) {
                dissi_inter.set(0.);
                if(process.nom_calcul=="incr") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        dissi_inter[j+1]=dissi_inter[j]+process.temps->dt*(
                                             dot((Inter[i].side[0].t_post[j+1].Fchap+Inter[i].side[0].t_post[j].Fchap)/2.0,Inter[i].side[0].M*Inter[i].side[0].t_post[j+1].Wpchap));
                    }
                } else if(process.nom_calcul=="latin") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        dissi_inter[j+1]=dissi_inter[j]+process.temps->dt*(
                                             dot((Inter[i].side[0].t[j+1].Fchap+Inter[i].side[0].t[j].Fchap)/2.0,Inter[i].side[0].M*Inter[i].side[0].t[j+1].Wpchap));
                    }
                } else {
                    std::cout << "Nom de calcul nom pris en compte" << endl;
                    assert(0);
                }
                dissipation+=dissi_inter;
                //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter << endl;
            }
        }
    }
}
// template <class TV1, class TI>
//         void calcul_ener_imp_chap(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
//             dissipation.set(0.);
//             Vec<double> dissi_inter;
//             dissi_inter.resize(dissipation.size());
//             for(unsigned j=0;j<S.size();j++) {
//                 for(unsigned e=0;e<S[j].edge.size();++e) {
//                     unsigned i=S[j].edge[e].internum;
//                     unsigned data=S[j].edge[e].datanum;
//                     if ((Inter[i].comp=="depl" or Inter[i].comp=="depl_normal" or Inter[i].comp=="effort" /*or Inter[i].comp=="Jeu_impose"*/) and data == 0) {
//                         dissi_inter.set(0.);
//                         if(process.nom_calcul=="incr") {
//                             for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
//                                 dissi_inter[j+1]=dissi_inter[j]+process.temps->dt/2*(
//                                         dot(Inter[i].side[0].t_post[j+1].Fchap,Inter[i].side[0].M*Inter[i].side[0].t_post[j+1].Wpchap)+
//                                         dot(Inter[i].side[0].t_post[j].Fchap,Inter[i].side[0].M*Inter[i].side[0].t_post[j].Wpchap));
//                             }
//                         } else if(process.nom_calcul=="latin") {
//                             for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
//                                 dissi_inter[j+1]=dissi_inter[j]+process.temps->dt/2*(
//                                         dot(Inter[i].side[0].t[j+1].Fchap,Inter[i].side[0].M*Inter[i].side[0].t[j+1].Wpchap)+
//                                         dot(Inter[i].side[0].t[j].Fchap,Inter[i].side[0].M*Inter[i].side[0].t[j].Wpchap));
//                             }
//                         } else {
//                             std::cout << "Nom de calcul nom pris en compte" << endl;
//                             assert(0);
//                         }
//                         dissipation+=dissi_inter;
//                         std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter << endl;
//                     }
//                 }
//             }
//         }

/**
Fonction permettant le calcul de l'énergie imposée à partir de la structure Inter et du Temps avec les quantités chapeaux.
*/
template <class TV1, class TI>
void calcul_ener_imp_lin(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
    dissipation.set(0.);
    Vec<double> dissi_inter;
    dissi_inter.resize(dissipation.size());
    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned i=S[j].edge[e].internum;
            unsigned data=S[j].edge[e].datanum;
            if ((Inter[i].comp=="depl" or Inter[i].comp=="depl_normal" or Inter[i].comp=="effort"  or Inter[i].comp=="effort_normale"/*or Inter[i].comp=="Jeu_impose"*/) and data == 0) {
                dissi_inter.set(0.);
                if(process.nom_calcul=="incr") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        dissi_inter[j+1]=dissi_inter[j]+process.temps->dt*(
                                             dot((Inter[i].side[0].t_post[j+1].F+Inter[i].side[0].t_post[j].F)/2.0,Inter[i].side[0].M*Inter[i].side[0].t_post[j+1].Wp));
                    }
                } else if(process.nom_calcul=="latin") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        dissi_inter[j+1]=dissi_inter[j]+process.temps->dt*(
                                             dot((Inter[i].side[0].t[j+1].F+Inter[i].side[0].t[j].F)/2.0,Inter[i].side[0].M*Inter[i].side[0].t[j+1].Wp));
                    }
                } else {
                    std::cout << "Nom de calcul nom pris en compte" << endl;
                    assert(0);
                }

                dissipation+=dissi_inter;
                //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter <<  endl;
            }
        }
    }
}
/*template <class TV1, class TI>
void calcul_ener_imp_lin(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
    dissipation.set(0.);
    Vec<double> dissi_inter;
    dissi_inter.resize(dissipation.size());
    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned i=S[j].edge[e].internum;
            unsigned data=S[j].edge[e].datanum;
            if ((Inter[i].comp=="depl" or Inter[i].comp=="depl_normal" or Inter[i].comp=="effort" or Inter[i].comp=="Jeu_impose") and data == 0) {
                dissi_inter.set(0.);
                if(process.nom_calcul=="incr") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        dissi_inter[j+1]=dissi_inter[j]+process.temps->dt/2*(
                                dot(Inter[i].side[0].t_post[j+1].F,Inter[i].side[0].M*Inter[i].side[0].t_post[j+1].Wp)+
                                dot(Inter[i].side[0].t_post[j].F,Inter[i].side[0].M*Inter[i].side[0].t_post[j].Wp));
                    }
                } else if(process.nom_calcul=="latin") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        dissi_inter[j+1]=dissi_inter[j]+process.temps->dt/2*(
                                dot(Inter[i].side[0].t[j+1].F,Inter[i].side[0].M*Inter[i].side[0].t[j+1].Wp)+
                                dot(Inter[i].side[0].t[j].F,Inter[i].side[0].M*Inter[i].side[0].t[j].Wp));
                    }
                } else {
                    std::cout << "Nom de calcul nom pris en compte" << endl;
                    assert(0);
                }

                dissipation+=dissi_inter;
                //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter <<  endl;
            }
        }
    }
}*/
        
//ajout de l'energie par element
struct add_ener_elem{
    template<class TE> void operator()(TE &e, typename TE::T &ener) const {
        ener+=e.ener;
    }
};

template <class TV1, class TI>
void calcul_energie_elastique(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
    dissipation.set(0.);
    Vec<double> dissi_inter;
    dissi_inter.resize(dissipation.size());
    if (process.latin->save_depl_SST!=1){
       if(process.rank == 0) std::cout << "ATTENTION il faut mettre save_depl_SST a 1 pour utiliser cette commande" << endl;
    } else {
     for(unsigned i=0;i<S.size();i++) {
        dissi_inter.set(0.);
        for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
            if(process.nom_calcul=="incr")
                assign_dep_cont_slave(S[i],S[i].t_post[j+1].q);
            else if(process.nom_calcul=="latin")
                assign_dep_cont_slave(S[i],S[i].t[j+1].q);

            apply(S[i].mesh->elem_list,add_ener_elem(),dissi_inter[j+1]);
            S[i].mesh.unload();
        }
        dissipation+=dissi_inter;
        //std::cout << "Contribution SST " << S[i].num << " : " << dissi_inter <<  endl;
    }
    }
}


/**
Fonction permettant le calcul de l'intégrale des efforts tangents au carré à partir de la structure Inter et du Temps avec les quantités chapeaux.
*/
template <class TV1, class TI>
void calcul_Ft2_chap(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
    dissipation.set(0.);
    Vec<double> dissi_inter,Ut,Un,tmp0,tmp1,Ft0,Ftj;
    dissi_inter.resize(dissipation.size());

    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned i=S[j].edge[e].internum;
            unsigned data=S[j].edge[e].datanum;
            Un.resize(Inter[i].side[0].nodeeq.size());
            Un.set(0.);
            Ut.resize(Inter[i].side[0].nodeeq.size());
            Ut.set(0.);
            tmp0.resize(Inter[i].side[0].nodeeq.size());
            tmp0.set(0.);
            tmp1.resize(Inter[i].side[0].nodeeq.size());
            tmp1.set(0.);
            Ft0.resize(Inter[i].side[0].nodeeq.size());
            Ft0.set(0.);
            Ftj.resize(Inter[i].side[0].nodeeq.size());
            Ftj.set(0.);
            if ((Inter[i].comp=="Contact" or Inter[i].comp=="Contact_jeu" or Inter[i].comp=="Contact_jeu_physique" or Inter[i].comp=="Contact_ep") and data == 0) {
                dissi_inter.set(0.);
                Vec<double,TI::template SubType<0>::T::dim> a0;
                a0.set(0.);
                std::cout << Inter[i].measure << endl;
                if(process.nom_calcul=="incr") {
                    tmp0=Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[1].Wchap);
                    tmp1=Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t_post[1].Wchap);
                    Un=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                    tmp0=Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t_post[1].Wchap);
                    tmp1=Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t_post[1].Wchap);
                    Ut=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                    Ft0=Inter[i].side[0].t_post[2].Fchap;
                    for( unsigned k=0;k< Un.size()/TI::template SubType<0>::T::dim; k++) {
                        if (norm_2(Un[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]) >
                                1e-6 or norm_2(Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)])
                                < 1e-6) {
                            Ft0[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                            =a0;
                        }
                    }
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].Wchap);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t_post[j+1].Wchap);
                        Un=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t_post[j+1].Wchap);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t_post[j+1].Wchap);
                        Ut=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        Ftj=Inter[i].side[0].t_post[j+1].Fchap;
                        for( unsigned k=0;k< Un.size()/TI::template SubType<0>::T::dim; k++) {
                            if (norm_2(Un[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]) >
                                    1e-6 or norm_2(Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)])
                                    < 1e-6) {
                                Ftj[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                                =a0;
                            }
                        }
                        dissi_inter[j+1]=dot(Inter[i].side[0].Pt(Ftj)-Inter[i].side[0].Pt(Ft0),Inter[i].side[0].M*(Inter[i].side[0].Pt(Ftj)-Inter[i].side[0].Pt(Ft0)));
                    }
                } else if(process.nom_calcul=="latin") {
                    tmp0=Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[1].Wchap);
                    tmp1=Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t[1].Wchap);
                    Un=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                    tmp0=Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t[1].Wchap);
                    tmp1=Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t[1].Wchap);
                    Ut=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                    Ft0=Inter[i].side[0].t[2].Fchap;
                    for( unsigned k=0;k< Un.size()/TI::template SubType<0>::T::dim; k++) {
                        if (norm_2(Un[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]) >
                                1e-6 or norm_2(Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)])
                                < 1e-6) {
                            Ft0[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                            =a0;
                        }
                    }
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].Wchap);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t[j+1].Wchap);
                        Un=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t[j+1].Wchap);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t[j+1].Wchap);
                        Ut=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        Ftj=Inter[i].side[0].t[j+1].Fchap;
                        for( unsigned k=0;k< Un.size()/TI::template SubType<0>::T::dim; k++) {
                            if (norm_2(Un[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]) >
                                    1e-6 or norm_2(Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)])
                                    < 1e-6) {
                                Ftj[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                                =a0;
                            }
                        }
                        dissi_inter[j+1]=dot(Inter[i].side[0].Pt(Ftj)-Inter[i].side[0].Pt(Ft0),Inter[i].side[0].M*(Inter[i].side[0].Pt(Ftj)-Inter[i].side[0].Pt(Ft0)));
                    }
                } else {
                    std::cout << "Nom de calcul nom pris en compte" << endl;
                    assert(0);
                }
                dissipation+=dissi_inter;
                //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter <<  endl;
            }
        }
    }
}

/**
Fonction permettant le calcul de l'intégrale des efforts tangents au carré à partir de la structure Inter et du Temps avec les quantités n.
*/
template <class TV1, class TI>
void calcul_Ft2_lin(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
    dissipation.set(0.);
    Vec<double> dissi_inter,Ut,Un,tmp0,tmp1,Ft0,Ftj;
    dissi_inter.resize(dissipation.size());

    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned i=S[j].edge[e].internum;
            unsigned data=S[j].edge[e].datanum;
            Un.resize(Inter[i].side[0].nodeeq.size());
            Un.set(0.);
            Ut.resize(Inter[i].side[0].nodeeq.size());
            Ut.set(0.);
            tmp0.resize(Inter[i].side[0].nodeeq.size());
            tmp0.set(0.);
            tmp1.resize(Inter[i].side[0].nodeeq.size());
            tmp1.set(0.);
            Ft0.resize(Inter[i].side[0].nodeeq.size());
            Ft0.set(0.);
            Ftj.resize(Inter[i].side[0].nodeeq.size());
            Ftj.set(0.);
            if ((Inter[i].comp=="Contact" or Inter[i].comp=="Contact_jeu" or Inter[i].comp=="Contact_jeu_physique" or Inter[i].comp=="Contact_ep") and data == 0) {
                dissi_inter.set(0.);
                Vec<double,TI::template SubType<0>::T::dim> a0;
                a0.set(0.);
                std::cout << Inter[i].measure << endl;
                if(process.nom_calcul=="incr") {
                    tmp0=Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[1].W);
                    tmp1=Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t_post[1].W);
                    Un=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                    tmp0=Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t_post[1].W);
                    tmp1=Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t_post[1].W);
                    Ut=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                    Ft0=Inter[i].side[0].t_post[2].F;
                    for( unsigned k=0;k< Un.size()/TI::template SubType<0>::T::dim; k++) {
                        if (norm_2(Un[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]) >
                                1e-6 or norm_2(Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)])
                                < 1e-6) {
                            Ft0[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                            =a0;
                        }
                    }
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].W);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t_post[j+1].W);
                        Un=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t_post[j+1].W);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t_post[j+1].W);
                        Ut=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        Ftj=Inter[i].side[0].t_post[j+1].F;
                        for( unsigned k=0;k< Un.size()/TI::template SubType<0>::T::dim; k++) {
                            if (norm_2(Un[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]) >
                                    1e-6 or norm_2(Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)])
                                    < 1e-6) {
                                Ftj[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                                =a0;
                            }
                        }
                        dissi_inter[j+1]=dot(Inter[i].side[0].Pt(Ftj)-Inter[i].side[0].Pt(Ft0),Inter[i].side[0].M*(Inter[i].side[0].Pt(Ftj)-Inter[i].side[0].Pt(Ft0)));
                    }
                } else if(process.nom_calcul=="latin") {
                    tmp0=Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[1].W);
                    tmp1=Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t[1].W);
                    Un=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                    tmp0=Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t[1].W);
                    tmp1=Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t[1].W);
                    Ut=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                    Ft0=Inter[i].side[0].t[2].F;
                    for( unsigned k=0;k< Un.size()/TI::template SubType<0>::T::dim; k++) {
                        if (norm_2(Un[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]) >
                                1e-6 or norm_2(Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)])
                                < 1e-6) {
                            Ft0[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                            =a0;
                        }
                    }
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].W);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t[j+1].W);
                        Un=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t[j+1].W);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t[j+1].W);
                        Ut=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        Ftj=Inter[i].side[0].t[j+1].F;
                        for( unsigned k=0;k< Un.size()/TI::template SubType<0>::T::dim; k++) {
                            if (norm_2(Un[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]) >
                                    1e-6 or norm_2(Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)])
                                    < 1e-6) {
                                Ftj[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                                =a0;
                            }
                        }
                        dissi_inter[j+1]=dot(Inter[i].side[0].Pt(Ftj)-Inter[i].side[0].Pt(Ft0),Inter[i].side[0].M*(Inter[i].side[0].Pt(Ftj)-Inter[i].side[0].Pt(Ft0)));
                    }

                } else {
                    std::cout << "Nom de calcul nom pris en compte" << endl;
                    assert(0);
                }
                dissipation+=dissi_inter;
                //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter << endl;
            }
        }
    }
}



/**
Fonction permettant le calcul de la moyenne des efforts normaux à partir de la structure Inter et du Temps avec les quantités chapeaux.
*/
template <class TV1, class TI>
void calcul_Fn_chap(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
    dissipation.set(0.);
    Vec<double> dissi_inter;
    dissi_inter.resize(dissipation.size());
    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned i=S[j].edge[e].internum;
            unsigned data=S[j].edge[e].datanum;
            if ((Inter[i].comp=="Contact" or Inter[i].comp=="Contact_jeu" or Inter[i].comp=="Contact_jeu_physique" or Inter[i].comp=="Contact_ep") and data == 0) {
                dissi_inter.set(0.);
                Vec<double> vecx,vecy,vecz;
                vecx.resize(Inter[i].side[0].nodeeq.size());
                vecx.set(0.);
                vecy.resize(Inter[i].side[0].nodeeq.size());
                vecy.set(0.);
                vecz.resize(Inter[i].side[0].nodeeq.size());
                vecz.set(0.);
                for( unsigned ii=0;ii<Inter[i].side[0].nodeeq.size()/TI::template SubType<0>::T::dim ;ii++ ) {
                    vecx[3*ii]=1.;
                    vecy[3*ii+1]=1.;
                    vecz[3*ii+2]=1.;
                }
                std::cout << Inter[i].measure << endl;
                if(process.nom_calcul=="incr") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        double a=dot(vecx,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].Fchap));
                        double b=dot(vecy,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].Fchap));
                        double c=dot(vecz,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].Fchap));
                        dissi_inter[j+1]=std::sqrt((pow(a,2)+pow(b,2)+pow(c,2)))/Inter[i].measure;
                    }
                } else if(process.nom_calcul=="latin") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        double a=dot(vecx,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].Fchap));
                        double b=dot(vecy,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].Fchap));
                        double c=dot(vecz,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].Fchap));
                        dissi_inter[j+1]=std::sqrt((pow(a,2)+pow(b,2)+pow(c,2)))/Inter[i].measure;
                    }
                } else {
                    std::cout << "Nom de calcul nom pris en compte" << endl;
                    assert(0);
                }
                dissipation+=dissi_inter;
                //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter <<  endl;
            }
        }
    }
}

/**
Fonction permettant le calcul de l'intégrale des efforts normaux à partir de la structure Inter et du Temps avec les quantités n.
*/
template <class TV1, class TI>
void calcul_Fn_lin(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
    dissipation.set(0.);
    Vec<double> dissi_inter;
    dissi_inter.resize(dissipation.size());
    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned i=S[j].edge[e].internum;
            unsigned data=S[j].edge[e].datanum;
            if ((Inter[i].comp=="Contact" or Inter[i].comp=="Contact_jeu" or Inter[i].comp=="Contact_jeu_physique" or Inter[i].comp=="Contact_ep") and data == 0) {
                dissi_inter.set(0.);
                Vec<double> vecx,vecy,vecz;
                vecx.resize(Inter[i].side[0].nodeeq.size());
                vecx.set(0.);
                vecy.resize(Inter[i].side[0].nodeeq.size());
                vecy.set(0.);
                vecz.resize(Inter[i].side[0].nodeeq.size());
                vecz.set(0.);
                for( unsigned ii=0;ii<Inter[i].side[0].nodeeq.size()/TI::template SubType<0>::T::dim ;ii++ ) {
                    vecx[3*ii]=1.;
                    vecy[3*ii+1]=1.;
                    vecz[3*ii+2]=1.;
                }
                std::cout << Inter[i].measure << endl;
                if(process.nom_calcul=="incr") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        double a=dot(vecx,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].F));
                        double b=dot(vecy,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].F));
                        double c=dot(vecz,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].F));
                        dissi_inter[j+1]=std::sqrt((pow(a,2)+pow(b,2)+pow(c,2)))/Inter[i].measure;
                    }
                } else if(process.nom_calcul=="latin") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        double a=dot(vecx,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].F));
                        double b=dot(vecy,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].F));
                        double c=dot(vecz,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].F));
                        dissi_inter[j+1]=std::sqrt((pow(a,2)+pow(b,2)+pow(c,2)))/Inter[i].measure;
                    }
                } else {
                    std::cout << "Nom de calcul nom pris en compte" << endl;
                    assert(0);
                }
                dissipation+=dissi_inter;
                //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter << endl;
            }
        }
    }
}







/**
Fonction permettant le calcul de la moyenne du saut de déplacement normal à partir de la structure Inter et du Temps avec les quantités chapeaux.
*/
template <class TV1, class TI>
void calcul_Un_chap(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
    dissipation.set(0.);
    Vec<double> dissi_inter;
    dissi_inter.resize(dissipation.size());
    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned i=S[j].edge[e].internum;
            unsigned data=S[j].edge[e].datanum;
            if ((Inter[i].comp=="Contact" or Inter[i].comp=="Contact_jeu" or Inter[i].comp=="Contact_jeu_physique" or Inter[i].comp=="Contact_ep") and data == 0) {
                dissi_inter.set(0.);
                Vec<double> vecx,vecy,vecz;
                vecx.resize(Inter[i].side[0].nodeeq.size());
                vecx.set(0.);
                vecy.resize(Inter[i].side[0].nodeeq.size());
                vecy.set(0.);
                vecz.resize(Inter[i].side[0].nodeeq.size());
                vecz.set(0.);
                for( unsigned ii=0;ii<Inter[i].side[0].nodeeq.size()/TI::template SubType<0>::T::dim ;ii++ ) {
                    vecx[3*ii]=1.;
                    vecy[3*ii+1]=1.;
                    vecz[3*ii+2]=1.;
                }
                std::cout << Inter[i].measure << endl;
                if(process.nom_calcul=="incr") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        double a=dot(vecx,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].Wchap)-
                                     Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t_post[j+1].Wchap));
                        double b=dot(vecy,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].Wchap)-
                                     Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t_post[j+1].Wchap));
                        double c=dot(vecz,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].Wchap)-
                                     Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t_post[j+1].Wchap));
                        dissi_inter[j+1]=std::sqrt((pow(a,2)+pow(b,2)+pow(c,2)))/Inter[i].measure;
                    }
                } else if(process.nom_calcul=="latin") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        double a=dot(vecx,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].Wchap)-
                                     Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t[j+1].Wchap));
                        double b=dot(vecy,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].Wchap)-
                                     Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t[j+1].Wchap));
                        double c=dot(vecz,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].Wchap)-
                                     Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t[j+1].Wchap));
                        dissi_inter[j+1]=std::sqrt((pow(a,2)+pow(b,2)+pow(c,2)))/Inter[i].measure;
                    }
                } else {
                    std::cout << "Nom de calcul nom pris en compte" << endl;
                    assert(0);
                }
                dissipation+=dissi_inter;
                //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter <<  endl;
            }
        }
    }
}

/**
Fonction permettant le calcul de la moyenne du saut de déplacement normal à partir de la structure Inter et du Temps avec les quantités n*/
template <class TV1, class TI>
void calcul_Un_lin(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
    dissipation.set(0.);
    Vec<double> dissi_inter;
    dissi_inter.resize(dissipation.size());
    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned i=S[j].edge[e].internum;
            unsigned data=S[j].edge[e].datanum;
            if ((Inter[i].comp=="Contact" or Inter[i].comp=="Contact_jeu" or Inter[i].comp=="Contact_jeu_physique" or Inter[i].comp=="Contact_ep") and data == 0) {
                dissi_inter.set(0.);
                Vec<double> vecx,vecy,vecz;
                vecx.resize(Inter[i].side[0].nodeeq.size());
                vecx.set(0.);
                vecy.resize(Inter[i].side[0].nodeeq.size());
                vecy.set(0.);
                vecz.resize(Inter[i].side[0].nodeeq.size());
                vecz.set(0.);
                for( unsigned ii=0;ii<Inter[i].side[0].nodeeq.size()/TI::template SubType<0>::T::dim ;ii++ ) {
                    vecx[3*ii]=1.;
                    vecy[3*ii+1]=1.;
                    vecz[3*ii+2]=1.;
                }
                std::cout << Inter[i].measure << endl;
                if(process.nom_calcul=="incr") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        double a=dot(vecx,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].W)-
                                     Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t_post[j+1].W));
                        double b=dot(vecy,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].W)-
                                     Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t_post[j+1].W));
                        double c=dot(vecz,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].W)-
                                     Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t_post[j+1].W));
                        dissi_inter[j+1]=std::sqrt((pow(a,2)+pow(b,2)+pow(c,2)))/Inter[i].measure;
                    }
                } else if(process.nom_calcul=="latin") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        double a=dot(vecx,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].W)-
                                     Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t[j+1].W));
                        double b=dot(vecy,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].W)-
                                     Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t[j+1].W));
                        double c=dot(vecz,Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].W)-
                                     Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t[j+1].W));
                        dissi_inter[j+1]=std::sqrt((pow(a,2)+pow(b,2)+pow(c,2)))/Inter[i].measure;
                    }
                } else {
                    std::cout << "Nom de calcul nom pris en compte" << endl;
                    assert(0);
                }
                dissipation+=dissi_inter;
                //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter << endl;
            }
        }
    }
}





/**
Fonction permettant le calcul de la moyenne du saut de déplacement tangent à partir de la structure Inter et du Temps avec les quantités chapeaux.
*/
template <class TV1, class TI>
void calcul_Ut_chap(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
    dissipation.set(0.);
    Vec<double> dissi_inter,surf;
    dissi_inter.resize(dissipation.size());
    surf.resize(dissipation.size());
    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned i=S[j].edge[e].internum;
            unsigned data=S[j].edge[e].datanum;
            if ((Inter[i].comp=="Contact" or Inter[i].comp=="Contact_jeu" or Inter[i].comp=="Contact_jeu_physique" or Inter[i].comp=="Contact_ep") and data == 0) {
                dissi_inter.set(0.);
                surf.set(0.);
                Vec<double> vecx,vecy,vecz,Un,Ut,tmp1,tmp0;
                vecx.resize(Inter[i].side[0].nodeeq.size());
                vecx.set(0.);
                vecy.resize(Inter[i].side[0].nodeeq.size());
                vecy.set(0.);
                vecz.resize(Inter[i].side[0].nodeeq.size());
                vecz.set(0.);
                Un.resize(Inter[i].side[0].nodeeq.size());
                Un.set(0.);
                Ut.resize(Inter[i].side[0].nodeeq.size());
                Ut.set(0.);
                tmp0.resize(Inter[i].side[0].nodeeq.size());
                tmp0.set(0.);
                tmp1.resize(Inter[i].side[0].nodeeq.size());
                tmp1.set(0.);
                for( unsigned ii=0;ii<Inter[i].side[0].nodeeq.size()/TI::template SubType<0>::T::dim ;ii++ ) {
                    vecx[3*ii]=1.;
                    vecy[3*ii+1]=1.;
                    vecz[3*ii+2]=1.;
                }
                std::cout << Inter[i].measure << endl;
                if(process.nom_calcul=="incr") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].Wchap);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t_post[j+1].Wchap);
                        Un=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t_post[j+1].Wchap);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t_post[j+1].Wchap);
                        Ut=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        Vec<double,TI::template SubType<0>::T::dim> a0;
                        a0.set(0.);
                        tmp0.set(1.);
                        for( unsigned k=0;k< Un.size()/TI::template SubType<0>::T::dim; k++) {
                            if (norm_2(Un[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]) >
                                    1e-6 or norm_2(Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)])
                                    < 1e-6) {
                                tmp0[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                                = a0;
                                Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                                =a0;
                            }
                        }
                        double a=dot(vecx,Ut);
                        double b=dot(vecy,Ut);
                        double c=dot(vecz,Ut);
                        dissi_inter[j+1]=std::sqrt((pow(a,2)+pow(b,2)+pow(c,2)))/Inter[i].measure;
                        surf[j+1]=dot(tmp0,Inter[i].side[0].M*tmp0)/TI::template SubType<0>
                        ::T::dim;
                    }
                } else if(process.nom_calcul=="latin") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].Wchap);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t[j+1].Wchap);
                        Un=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t[j+1].Wchap);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t[j+1].Wchap);
                        Ut=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        Vec<double,TI::template SubType<0>::T::dim> a0;
                        a0.set(0.);
                        tmp0.set(1.);
                        for( unsigned k=0;k< Un.size()/TI::template SubType<0>::T::dim; k++) {
                            if (norm_2(Un[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]) >
                                    1e-6 or norm_2(Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)])
                                    < 1e-6) {
                                tmp0[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                                = a0;
                                Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                                =a0;
                            }
                        }
                        double a=dot(vecx,Ut);
                        double b=dot(vecy,Ut);
                        double c=dot(vecz,Ut);
                        dissi_inter[j+1]=std::sqrt((pow(a,2)+pow(b,2)+pow(c,2)))/Inter[i].measure;
                        surf[j+1]=dot(tmp0,Inter[i].side[0].M*tmp0)/TI::template SubType<0>
                        ::T::dim;
                    }
                } else {
                    std::cout << "Nom de calcul nom pris en compte" << endl;
                    assert(0);
                }
                dissipation+=dissi_inter;
                //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter <<  endl;
                std::cout << "\t\t" << surf <<  endl;
            }
        }
    }
}

/**
Fonction permettant le calcul de la moyenne du saut de déplacement tangent à partir de la structure Inter et du Temps avec les quantités n*/
template <class TV1, class TI>
void calcul_Ut_lin(TV1 &S, TI &Inter,Vec<double> &dissipation,Param &process) {
    dissipation.set(0.);
    Vec<double> dissi_inter,surf;
    dissi_inter.resize(dissipation.size());
    surf.resize(dissipation.size());
    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned i=S[j].edge[e].internum;
            unsigned data=S[j].edge[e].datanum;
            if ((Inter[i].comp=="Contact" or Inter[i].comp=="Contact_jeu" or Inter[i].comp=="Contact_jeu_physique" or Inter[i].comp=="Contact_ep") and data == 0) {
                dissi_inter.set(0.);
                surf.set(0.);
                Vec<double> vecx,vecy,vecz,Un,Ut,tmp1,tmp0;
                vecx.resize(Inter[i].side[0].nodeeq.size());
                vecx.set(0.);
                vecy.resize(Inter[i].side[0].nodeeq.size());
                vecy.set(0.);
                vecz.resize(Inter[i].side[0].nodeeq.size());
                vecz.set(0.);
                Un.resize(Inter[i].side[0].nodeeq.size());
                Un.set(0.);
                Ut.resize(Inter[i].side[0].nodeeq.size());
                Ut.set(0.);
                tmp0.resize(Inter[i].side[0].nodeeq.size());
                tmp0.set(0.);
                tmp1.resize(Inter[i].side[0].nodeeq.size());
                tmp1.set(0.);
                for( unsigned ii=0;ii<Inter[i].side[0].nodeeq.size()/TI::template SubType<0>::T::dim ;ii++ ) {
                    vecx[3*ii]=1.;
                    vecy[3*ii+1]=1.;
                    vecz[3*ii+2]=1.;
                }
                std::cout << Inter[i].measure << endl;
                if(process.nom_calcul=="incr") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t_post[j+1].W);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t_post[j+1].W);
                        Un=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t_post[j+1].W);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t_post[j+1].W);
                        Ut=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        Vec<double,TI::template SubType<0>::T::dim> a0;
                        a0.set(0.);
                        tmp0.set(1.);
                        for( unsigned k=0;k< Un.size()/TI::template SubType<0>::T::dim; k++) {
                            if (norm_2(Un[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]) >
                                    1e-6 or norm_2(Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)])
                                    < 1e-6) {
                                tmp0[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                                = a0;
                                Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                                =a0;
                            }
                        }
                        double a=dot(vecx,Ut);
                        double b=dot(vecy,Ut);
                        double c=dot(vecz,Ut);
                        dissi_inter[j+1]=std::sqrt((pow(a,2)+pow(b,2)+pow(c,2)))/Inter[i].measure;
                        surf[j+1]=dot(tmp0,Inter[i].side[0].M*tmp0)/TI::template SubType<0>
                        ::T::dim;
                    }
                } else if(process.nom_calcul=="latin") {
                    for(unsigned j=0 ;j<dissi_inter.size()-1 ;j++ ) {
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pn(Inter[i].side[0].t[j+1].W);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pn(Inter[i].side[1].t[j+1].W);
                        Un=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        tmp0=Inter[i].side[0].M*Inter[i].side[0].Pt(Inter[i].side[0].t[j+1].W);
                        tmp1=Inter[i].side[1].M*Inter[i].side[1].Pt(Inter[i].side[1].t[j+1].W);
                        Ut=tmp0[Inter[i].side[0].ddlcorresp]-tmp1[Inter[i].side[1].ddlcorresp];
                        Vec<double,TI::template SubType<0>::T::dim> a0;
                        a0.set(0.);
                        tmp0.set(1.);
                        for( unsigned k=0;k< Un.size()/TI::template SubType<0>::T::dim; k++) {
                            if (norm_2(Un[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]) >
                                    1e-6 or norm_2(Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)])
                                    < 1e-6) {
                                tmp0[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                                = a0;
                                Ut[range(k*TI::template SubType<0>::T::dim,(k+1)*TI::template SubType<0>::T::dim)]
                                =a0;
                            }
                        }
                        double a=dot(vecx,Ut);
                        double b=dot(vecy,Ut);
                        double c=dot(vecz,Ut);
                        dissi_inter[j+1]=std::sqrt((pow(a,2)+pow(b,2)+pow(c,2)))/Inter[i].measure;
                        surf[j+1]=dot(tmp0,Inter[i].side[0].M*tmp0)/TI::template SubType<0>
                        ::T::dim;
                    }
                } else {
                    std::cout << "Nom de calcul nom pris en compte" << endl;
                    assert(0);
                }
                dissipation+=dissi_inter;
                //std::cout << "Contribution interface " << Inter[i].num << " entre les pieces " << Inter[i].vois << " : " << dissi_inter << endl;
                std::cout << "\t\t" << surf <<  endl;
            }
        }
    }
}



/*
./depouille_resu_bride 80 5 16 bride_90/bride_talon_90_f_0.*cisaillementR*resu
./depouille_resu_bride 100 5 20 bride_90/bride_talon_90_f_0.*traction*resu
./depouille_resu_bride 120 5 24 bride_90/bride_talon_90_f_0.*cisaillementT*resu
./depouille_resu_bride 100 5 20 bride_90/bride_talon_90_f_0.*Mt*resu
*/


