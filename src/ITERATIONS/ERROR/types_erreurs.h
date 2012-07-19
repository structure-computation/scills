/** \ingroup   calcul_erreur
 * \brief Calcul d'une erreur en determinant l'energie par interface : 
 * 
 * En bouclant sur les sous-structures et sur les interfaces entourant chaque sous-structure, on calcul les termes suivants :
 * - numérateur : \f$ (\int_{\Gamma} F*\dot{W} - Fchap*\dot{\hat{W}} )^2  \f$ en quasistatique
 * - dénominateur : \f$ (\int_{\Gamma} F*\dot{W} + Fchap*\dot{\hat{W}} )^2 \f$ en quasistatique
 * que l'on somme au numérateur et dénominateur précédents.
 */
struct calcerror_ener {
    template<class SST, class TV2>
    void  operator()(SST &S, TV2 &Inter, Vec<TYPEREEL> &frac, Process &process) const {
        
        unsigned data,q;
        
        unsigned imic = process.temps->pt;
        for(unsigned j=0;j<S.edge.size();++j) {
            q=S.edge[j].internum;
            data=S.edge[j].datanum;
            Vec<TYPEREEL> tempF,tempW,tempFchap, tempWchap;
            tempF=Inter[q].side[data].t[imic].F;
            tempW=Inter[q].side[data].t[imic].Wp;
            tempWchap=Inter[q].side[data].t[imic].Wpchap;
            tempFchap=Inter[q].side[data].t[imic].Fchap;
            frac[0] += std::pow(dot(tempF,Inter[q].side[data].M*tempW) - dot(tempFchap,Inter[q].side[data].M*tempWchap) ,2);
            frac[1] += std::pow(dot(tempF,Inter[q].side[data].M*tempW) + dot(tempFchap,Inter[q].side[data].M*tempWchap) ,2);
        }
    }
};

/** \ingroup   calcul_erreur
 * \brief Calcul d'une erreur en résidu sur le déplacement (ou vitesses) des interfaces : 
 * 
 * En bouclant sur les cotes des sous-structures, on calcule les termes suivants :
 * - \f$ \Vert \int_{\Gamma} \hat{Wp} - \dot{Wp}) \Vert_2  \f$ 
 */
struct calcerror_residu_depl_post {
    template<class SST, class TV2>
    void  operator()(SST &S, TV2 &Inter, Vec<TYPEREEL> &frac, Process &process) const {
        
        unsigned data,q;
        
        unsigned imic = process.temps->pt;
        for(unsigned j=0;j<S.edge.size();++j) {
            q=S.edge[j].internum;
            data=S.edge[j].datanum;
            if ((Inter[q].type == "Int" or Inter[q].comp == "periodique") and data == 0) {
                Vec<TYPEREEL> W1,W2;
                W1=Inter[q].side[data].t_post[imic].Wp;
                W2=Inter[q].side[1-data].t_post[imic].Wp;
                frac[0] += std::pow(dot(W1-W2,Inter[q].side[data].M*(W1-W2)) ,0.5);
            }
        }
    }
};
struct calcerror_residu_depl {
    template<class SST, class TV2>
    void  operator()(SST &S, TV2 &Inter, Vec<TYPEREEL> &frac, Process &process) const {
        
        unsigned data,q;
        
        unsigned imic = process.temps->pt;
        for(unsigned j=0;j<S.edge.size();++j) {
            q=S.edge[j].internum;
            data=S.edge[j].datanum;
            if ((Inter[q].type == "Int" or Inter[q].comp == "periodique") and data == 0) {
                Vec<TYPEREEL> W1,W2;
                W1=Inter[q].side[data].t[imic].Wp;
                W2=Inter[q].side[1-data].t[imic].Wp;
                frac[0] += std::pow(dot(W1-W2,Inter[q].side[data].M*(W1-W2)) ,0.5);
            }
        }
    }
};

/** \ingroup   calcul_erreur
 * \brief Calcul d'un indicateur d'erreur entre les solutions \f$ \mathbf{s} \f$ et \f$ \mathbf{\hat{s}} \f$ avec une norme sur les directions de recherche
 * 
 * En bouclant sur les sous-structures et sur les interfaces entourant chaque sous-structure, on calcule les termes suivants :
 * - numérateur : \f$ \int_{\Gamma} h(F-\hat{F})^2 + k(\dot{W} - \dot{\hat{W}})^2  \f$ en quasistatique
 * - dénominateur : \f$ \int_{\Gamma} h(F+\hat{F})^2 + k(\dot{W} + \dot{\hat{W}})^2  \f$ en quasistatique
 * que l'on somme au numérateur et dénominateur précédents.
 */
struct calcerror_ddr {
    template<class SST, class TV2>
    void  operator()(SST &S, TV2 &Inter, Vec<TYPEREEL> &frac, Process &process) const {
        
        Vec<TYPEREEL,2> errF, errW;
        errF.set(0.0);
        errW.set(0.0);
        
        unsigned data,q;
        unsigned imic = process.temps->pt;
        
        for(unsigned j=0;j<S.edge.size();++j) {
            q=S.edge[j].internum;
            data=S.edge[j].datanum;
            Vec<TYPEREEL> &tempF=Inter[q].side[data].t[imic].F;
            Vec<TYPEREEL> &tempW=Inter[q].side[data].t[imic].Wp;
            Vec<TYPEREEL> &tempWchap=Inter[q].side[data].t[imic].Wpchap;
            Vec<TYPEREEL> &tempFchap=Inter[q].side[data].t[imic].Fchap;
            const Vec<TYPEREEL> &JJ=Inter[q].jeu;
            Vec<TYPEREEL> temp=tempF-tempFchap;
            TYPEREEL temp0=dot(temp,Inter[q].side[data].M*Inter[q].side[data].hglo*(temp));
            temp=tempW-tempWchap;
            TYPEREEL temp1=dot(temp,Inter[q].side[data].M*Inter[q].side[data].kglo*(temp));
            temp=tempF+tempFchap;
            errF[1] +=dot(temp,Inter[q].side[data].M*Inter[q].side[data].hglo*(temp));
            //temp=tempW-JJ+tempWchap;  // TEST
            temp=tempW+tempWchap;  // TEST
            errW[1] +=dot(temp,Inter[q].side[data].M*Inter[q].side[data].kglo*(temp));
            frac[2+q]+=temp0+temp1;
            errF[0]+=temp0;
            errW[0]+=temp1;
            //cout << "sst " << S.num << " edge " << j << " interface " << q << " cote " << data << " type et comp " << Inter[q].type << " " << Inter[q].comp << " nombres : " << dot(tempF-tempFchap,Inter[q].side[data].M*Inter[q].side[data].hglo*(tempF-tempFchap)) << " " << dot(tempF+tempFchap,Inter[q].side[data].M*Inter[q].side[data].hglo*(tempF+tempFchap)) << " et " <<  dot(tempW-tempWchap,Inter[q].side[data].M*Inter[q].side[data].kglo*(tempW-tempWchap)) << " " << dot(tempW+tempWchap,Inter[q].side[data].M*Inter[q].side[data].kglo*(tempW+tempWchap)) << " et " << norm_2(tempF) << " " << norm_2(tempFchap) << " " <<norm_2(tempW) << " " <<norm_2(tempWchap) << endl;
        }
        frac[0]+=errF[0]+errW[0];
        frac[1]+=errF[1]+errW[1];
    }
    
};

/** \ingroup   calcul_erreur
 * \brief Calcul d'un indicateur d'erreur entre les solutions \f$ \mathbf{s} \f$ et \f$ \mathbf{\hat{s}} \f$ avec une norme sur les directions de recherche
 * 
 * En bouclant sur les sous-structures et sur les interfaces entourant chaque sous-structure, on calcule les termes suivants :
 * - numérateur : \f$ \int_{\Gamma} h(F-\hat{F})^2 + k(\dot{W} - \dot{\hat{W}})^2  \f$ en quasistatique
 * - dénominateur : \f$ \int_{\Gamma} h(F+\hat{F})^2 + k(\dot{W} + \dot{\hat{W}})^2  \f$ en quasistatique
 * que l'on somme au numérateur et dénominateur précédents.
 */
struct calcerror_ddr_post {
    template<class SST, class TV2>
    void  operator()(SST &S, TV2 &Inter, Vec<TYPEREEL> &frac, Process &process) const {
        
        Vec<TYPEREEL,2> errF, errW;
        errF.set(0.0);
        errW.set(0.0);
        
        unsigned data,q;
        unsigned imic = process.temps->pt;
        
        
        for(unsigned j=0;j<S.edge.size();++j) {
            q=S.edge[j].internum;
            data=S.edge[j].datanum;
            Vec<TYPEREEL> &tempF=Inter[q].side[data].t_post[imic].F;
            Vec<TYPEREEL> &tempW=Inter[q].side[data].t_post[imic].Wp;
            Vec<TYPEREEL> &tempWchap=Inter[q].side[data].t_post[imic].Wpchap;
            Vec<TYPEREEL> &tempFchap=Inter[q].side[data].t_post[imic].Fchap;
            Vec<TYPEREEL> temp=tempF-tempFchap;
            TYPEREEL temp0=dot(temp,Inter[q].side[data].M*Inter[q].side[data].hglo*(temp));
            temp=tempW-tempWchap;
            TYPEREEL temp1=dot(temp,Inter[q].side[data].M*Inter[q].side[data].kglo*(temp));
            temp=tempF+tempFchap;
            errF[1] +=dot(temp,Inter[q].side[data].M*Inter[q].side[data].hglo*(temp));
            temp=tempW+tempWchap;
            errW[1] +=dot(temp,Inter[q].side[data].M*Inter[q].side[data].kglo*(temp));
            frac[2+q]+=temp0+temp1;
            errF[0]+=temp0;
            errW[0]+=temp1;
        }
        frac[0]+=errF[0]+errW[0];
        frac[1]+=errF[1]+errW[1];
    }
    
};




/** \ingroup   calcul_erreur
 * \brief Calcul d'un indicateur d'erreur entre les solutions \f$ \mathbf{s} \f$ et \f$ \mathbf{\hat{s}} \f$ avec une norme sur la dissipation
 * 
 * On calcule l'energie dissipee sur les interfaces de contact a partir des quantites chapeaux et des quantites n, puis on effectue l erreur relation entre les dissipations
 */
struct calcerror_dissi {
    template<class SST, class TV2>
    void  operator()(SST &S, TV2 &Inter, Vec<TYPEREEL> &frac, Process &process) const {
        
        unsigned data,q;
        
        for(unsigned i=0;i<S.edge.size();++i) {
            q=S.edge[i].internum;
            data=S.edge[i].datanum;
            TYPEREEL temp1=0;
            TYPEREEL temp2=0;
            if (Inter[q].comp=="Contact_jeu_physique" or Inter[q].comp=="Contact_ep"){
                for(unsigned j=0 ;j<process.temps->nbpastemps ;j++ ) {
                    temp1+=process.temps->dt*(
                        dot((Inter[q].side[data].Pt(Inter[q].side[data].t[j+1].Fchap)+Inter[q].side[data].Pt(Inter[q].side[data].t[j].Fchap))/2.0,Inter[q].side[data].M*Inter[q].side[data].Pt(Inter[q].side[data].t[j+1].Wpchap)));
                    temp2+=process.temps->dt*(
                        dot((Inter[q].side[data].Pt(Inter[q].side[data].t[j+1].F)+Inter[q].side[data].Pt(Inter[q].side[data].t[j].F))/2.0,Inter[q].side[data].M*Inter[q].side[data].Pt(Inter[q].side[data].t[j+1].Wp)));
                }
                frac[0]+=temp1;
                frac[2+q]+=temp1-temp2;
                frac[1]+=temp2;
            }
        }
    }
};

struct calcerror_dissi_post {
    template<class SST, class TV2>
    void  operator()(SST &S, TV2 &Inter, Vec<TYPEREEL> &frac, Process &process) const {
        
        unsigned data,q;
        
        for(unsigned i=0;i<S.edge.size();++i) {
            q=S.edge[i].internum;
            data=S.edge[i].datanum;
            TYPEREEL temp1=0;
            TYPEREEL temp2=0;
            if (Inter[q].comp=="Contact_jeu_physique" or Inter[q].comp=="Contact_ep"){
                for(unsigned j=0 ;j<(unsigned)process.temps->pt_cur ;j++ ) {
                    temp1+=process.temps->dt*(
                        dot((Inter[q].side[data].Pt(Inter[q].side[data].t_post[j+1].Fchap)+Inter[q].side[data].Pt(Inter[q].side[data].t_post[j].Fchap))/2.0,Inter[q].side[data].M*Inter[q].side[data].Pt(Inter[q].side[data].t_post[j+1].Wpchap)));
                    temp2+=process.temps->dt*(
                        dot((Inter[q].side[data].Pt(Inter[q].side[data].t_post[j+1].F)+Inter[q].side[data].Pt(Inter[q].side[data].t_post[j].F))/2.0,Inter[q].side[data].M*Inter[q].side[data].Pt(Inter[q].side[data].t_post[j+1].Wp)));
                }
                frac[0]+=temp1;
                frac[2+q]+=temp1-temp2;
                frac[1]+=temp2;
            }
        }
    }
};
