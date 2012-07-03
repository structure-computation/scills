#ifndef PARAM_COMP_INTER_H
#define PARAM_COMP_INTER_H

#include "main_typedef.h"

struct PARAM_DAMAGE_INTER{
    PARAM_DAMAGE_INTER(){Yc=0.01; Yo=0.; n=0.5; alpha=1.5; gamma=0.3; kn=10000; knc=10000; kt=10000;}
    Scalar kn, kt,knc;
    Scalar gamma,alpha, Yc, Yo, n;
};
 
struct PARAM_COMP_INTER{
    Scalar coeffrottement;
    Vector f_coeffrottement;
    Vector jeu;
    Sc2String fcts_spatiales;
    Scalar Gcrit;
    bool degradable;
    Vector d;//,dmax;
    //Vector Y;//,Ymax;
    unsigned nbpastempsimpos;
    LMT::Vec<bool> comportement;
    int convergence; ///< =-1 si le calcul du pas de temps n'est pas convergence, >=0 sinon
                     ///< =0 après l'etape locale si aucun comportement d'element est mis à jour, >0 sinon

    PARAM_COMP_INTER(unsigned nbnodeeq){
        jeu.resize(nbnodeeq,0.); 
        //Ymax.resize(nbnodeeq,0.);
        //Y.resize(nbnodeeq,0.);
        d.resize(nbnodeeq,0.);
        //dmax.resize(nbnodeeq,0.);
        f_coeffrottement.resize(nbnodeeq,0.);
        comportement.resize(nbnodeeq,0);
        convergence = -1;
    }
    
    void free(){
        jeu.free();
        //dmax.free();
        d.free();
        //Y.free();
        //Ymax.free();
        t.free();
        comportement.free();
    }

    struct Time{
        Vector d;//,Y;//,dmax,Ymax;
        void allocate(unsigned nbnodeeq)
        {
            //Ymax.resize(nbnodeeq,0.);
            //Y.resize(nbnodeeq,0.);
            d.resize(nbnodeeq,0.);
            //dmax.resize(nbnodeeq,0.);
        }
    };
    Vec<Time> t;
    
    typedef PARAM_DAMAGE_INTER PARAM_DAMAGE;
    PARAM_DAMAGE param_damage;
};

#endif //PARAM_COMP_INTER_H

