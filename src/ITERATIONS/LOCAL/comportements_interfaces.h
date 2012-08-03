#ifndef COMPORTEMENT_INTERFACES_H
#define COMPORTEMENT_INTERFACES_H

#include "../../DEFINITIONS/SstCarac_InterCarac.h"
#include "../../DEFINITIONS/Interface.h"
#include "../../UTILITAIRES/utilitaires.h"
#include "../../DEFINITIONS/TimeData.h"


//interface exterieure de type deplacement impose pour tous les ddls
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procedure pour les interfaces exterieures à deplacement (ou vitesse) impose
*/
void compt_CL_vit (Interface &Inter,int &imic) {
    Inter.side[0].t[imic].Fchap=Inter.side[0].t[imic].F + Inter.side[0].kglo*(Inter.side[0].t[imic].Wpchap - Inter.side[0].t[imic].Wp);
}

//interface exterieure de type effort impose pour tous les ddls
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procedure pour les interfaces exterieures à effort impose
*/
void compt_CL_eff (Interface &Inter,int &imic) {
    Inter.side[0].t[imic].Wpchap=Inter.side[0].t[imic].Wp+Inter.side[0].hglo*(Inter.side[0].t[imic].Fchap - Inter.side[0].t[imic].F);
}

//interface exterieure de type : deplacement normal a l'interface plane impose + effort tangentiel nul
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procedure pour les interfaces exterieures de type symetrie ou ayant seulement le deplacement normal donne (ou vitesse)
*/
void compt_CL_sym (Interface &Inter,int &imic) {
    Vector &WWchap1 = Inter.side[0].t[imic].Wpchap;
    Vector &QQchap1 = Inter.side[0].t[imic].Fchap;
    const Vector &Q1 = Inter.side[0].t[imic].F;
    const Vector &WW1 = Inter.side[0].t[imic].Wp;

    Scalar ht = Inter.side[0].ht;
    Scalar kt = Inter.side[0].kt;
    Scalar kn = Inter.side[0].kn;

    Vector neq = Inter.side[0].neq;
    
//     PRINT(Inter.side[0].t[imic].Wp);

    for(unsigned i = 0; i < Inter.side[0].nodeeq.size(); ++i) {
        Vec<unsigned> rep = range(i*DIM,(i+1)*DIM);
        Point n = neq[rep];             /// normale de 1 vers 2
        Point F1 = Q1[rep];
        Point Fchap1 = QQchap1[rep];
        Point W1 = WW1[rep];
        Point Wchap1 = WWchap1[rep];
        Point PtWchap1;
        Point PnFchap1;

        Scalar Wchap1n = dot(Wchap1,n); ///donnee = 0 symetrie ou !=0 pour depl normal
        Wchap1.set(0.);
        Wchap1 = Wchap1n*n + ht * ( kt * ProjT(W1,n) - ProjT(F1,n) );
        Fchap1.set(0.);
        Fchap1 = dot(F1,n)*n + kn * (Wchap1n*n - dot(W1,n)*n);

        WWchap1[rep] = Wchap1;
        QQchap1[rep] = Fchap1;
    }

    Inter.side[0].t[imic].Wpchap=WWchap1;
    Inter.side[0].t[imic].Fchap=QQchap1;
    
//     PRINT(Inter.side[0].t[imic].Wpchap);
}

void comportement_local_interface(Interface &Inter, unsigned pt, Scalar dt);

/// Structure permettant de definir l'operateur de direction de recherche local a partir de direction normale et tangentielle
/// S'utilise ensuite comme une matrice
struct Kloc{
    /// Rigidites normale (kn) et tangentielle (kt)
    Scalar kn,kt;
    /// Vecteur normal
    Point n;
    /// Surcharge pour la multiplication d'un vecteur W
    template<class TV> Point operator*(TV &W){return kn*n*dot(n,W)+kt*(W-n*dot(n,W));}
    /// Surcharge pour la multiplication d'un vecteur W constant
    template<class TV> Point operator*(const TV &W){return kn*n*dot(n,W)+kt*(W-n*dot(n,W));}
};


//interface interieure de type parfaite
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procedure pour les interfaces interieures parfaites
*/
void compt_parfait (Interface &Inter,TimeData &temps) {
    const unsigned pt = temps.pt;
    const Scalar dt = temps.dt;
    
    LMT::Vec<unsigned> &list1 = Inter.side[0].ddlcorresp;
    LMT::Vec<unsigned> &list2 = Inter.side[1].ddlcorresp;
    const Vector &n1_ref = Inter.side[0].neq[list1];
    const Vector &F1_ref = Inter.side[0].t[pt].F[list1];
    const Vector &F2_ref = Inter.side[1].t[pt].F[list2];
    const Vector &Wp1_ref = Inter.side[0].t[pt].Wp[list1];
    const Vector &Wp2_ref = Inter.side[1].t[pt].Wp[list2];
    const Vector &Ep_imposee_ref = Inter.Ep_imposee[list1];
    const Vector &old_Ep_imposee_ref = Inter.old_Ep_imposee[list1];
    const Vector &Precharge_ref = Inter.precharge[list1];
    const Vector &old_Precharge_ref = Inter.old_precharge[list1];
    Vector Wpchap1 = Inter.side[0].t[pt].Wpchap[list1];
    Vector Wpchap2 = Inter.side[1].t[pt].Wpchap[list2];
    Vector Fchap1 = Inter.side[0].t[pt].Fchap[list1];
    Vector Fchap2 = Inter.side[1].t[pt].Fchap[list2];

    /// Creation des operateurs locaux de direction de recherche
    Kloc k1;
    k1.kn=Inter.side[0].kn;
    k1.kt=Inter.side[0].kt;
    
    Kloc k2;
    k2.kn=Inter.side[1].kn;
    k2.kt=Inter.side[1].kt;
    
    Kloc h;
    h.kn=1./(Inter.side[1].kn+Inter.side[0].kn);
    h.kt=1./(Inter.side[1].kt+Inter.side[0].kt);
    
    /// Travail point par point
    for(unsigned i = 0; i < Inter.side[0].nodeeq.size(); i++) {
        /// Creation du reperage du point
        LMT::Vec<unsigned> rep = range(i*DIM,(i+1)*DIM);
        
        /// Recuperation des vecteurs locaux
        Point n1 = n1_ref[rep];
        Point F1 = F1_ref[rep];
        Point F2 = F2_ref[rep];
        Point Wp1 = Wp1_ref[rep];
        Point Wp2 = Wp2_ref[rep];
        Point dEp_imposee = (Ep_imposee_ref[rep] - old_Ep_imposee_ref[rep])/dt;
        Point dPrecharge = Precharge_ref[rep] - old_Precharge_ref[rep];
        
        /// Creation des matrices elementaires de direction de recherche
        k1.n = n1;
        k2.n = n1;
        h.n = n1;
        
        /// Calcul des grandeurs
        Wpchap1[rep] = h*( (k1*Wp1 + k2*Wp2) - (F1 + F2) - k2*dEp_imposee);
        Wpchap2[rep] = Wpchap1[rep] + dEp_imposee;
        Fchap1[rep] = F1 + k1*(Wpchap1-Wp1) - dPrecharge;
    }

    Inter.side[0].t[pt].Wpchap[list1] = Wpchap1;
    Inter.side[1].t[pt].Wpchap[list2] = Wpchap2;
    Inter.side[0].t[pt].Fchap[list1] = Fchap1;
    Inter.side[1].t[pt].Fchap[list2] = -1.0*Fchap1;
}


//interface interieure de type contact avec frottement
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procedure pour les interfaces interieures de type contact avec frottement, avec jeu ou non
*/
void compt_contact (Interface &Inter,TimeData &temps) {/*
    const unsigned pt = temps.pt;
    const Scalar dt = temps.dt;
    
    LMT::Vec<unsigned> &list1 = Inter.side[0].ddlcorresp;
    LMT::Vec<unsigned> &list2 = Inter.side[1].ddlcorresp;
    const Vector &n1_ref = Inter.side[0].neq[list1];
    const Vector &F1_ref = Inter.side[0].t[pt].F[list1];
    const Vector &F2_ref = Inter.side[1].t[pt].F[list2];
    const Vector &Wp1_ref = Inter.side[0].t[pt].Wp[list1];
    const Vector &Wp2_ref = Inter.side[1].t[pt].Wp[list2];
    const Vector &Ep_imposee_ref = Inter.Ep_imposee[list1];
    const Vector &old_Ep_imposee_ref = Inter.old_Ep_imposee[list1];
    const Vector &old_Wchap1_ref = Inter.side[0].t[pt-1].Wp[list1];
    const Vector &old_Wchap2_ref = Inter.side[1].t[pt-1].Wp[list2];
    Vector Wpchap1 = Inter.side[0].t[pt].Wpchap[list1];
    Vector Wpchap2 = Inter.side[1].t[pt].Wpchap[list2];
    Vector Fchap1 = Inter.side[0].t[pt].Fchap[list1];
    Vector Fchap2 = Inter.side[1].t[pt].Fchap[list2];
    
    Scalar f=Inter.coeffrottement;

    Scalar h1n=Inter.side[0].hn;
    Scalar h2n=Inter.side[1].hn;
    Scalar h1t=Inter.side[0].ht;
    Scalar h2t=Inter.side[1].ht;

    // travail pt par pt
    const unsigned nbpts=Inter.side[0].nodeeq.size();
// #ifdef LOOK_CONTACT_ZONE
    //     Vec<unsigned> status;status.resize(nbpts*DIM);
// #endif

    for(unsigned i = 0; i < nbpts; ++i) {
        /// Recuperation des indices de stockage
        Vec<unsigned> rep=range(i*DIM,(i+1)*DIM);
        /// Recuperation des vecteurs locaux
        Point n = n1_ref[rep];                      /// Normale de 1 vers 2
        Point F1 = F1_ref[rep];                     /// Efforts du solide 1 sur l'interface calcules a l'etape lineaire
        Point F2 = F2_ref[rep];                     /// Efforts du solide 2 sur l'interface calcules a l'etape lineaire
        Point Wp1 = Wp1_ref[rep];                   /// Vitesse du solide 1 calculee a l'etape lineaire
        Point Wp2 = Wp2_ref[rep];                   /// Vitesse du solide 2 calculee a l'etape lineaire
        Point old_Wchap1 = old_Wchap1_ref[rep];     /// Deplacement du solide 1 calculee a la derniere etape locale du pas de temps precedent
        Point old_Wchap2 = old_Wchap2_ref[rep];     /// Deplacement du solide 2 calculee a la derniere etape locale du pas de temps precedent
        
        Point 
        
        /// Test de contact
        Scalar dWchap_n = dot(n1,(old_Wchap2-old_Wchap1)) + dt*(dot(n1,(Wp2-Wp1)) - (h2.kn*dot(n1,F2) - h1.kn*dot(n1,F1)));
        if (dWchap_n > Ep_n) {
            /// separation de l'interface
            Fchap1 = 0.0;
            Fchap2 = 0.0;
            Wpchap1[rep] = Wp1-h1n*dot(n,F1)*n - h1t* ProjT(F1,n);
            Wpchap2[rep] = Wp2-h2n*dot(n,F2)*n - h2t* ProjT(F2,n);
        } else { 
            /// contact
            Scalar Fchap1n = N;
            Scalar Fchap2n = -1.0*Fchap1n;
            Scalar Wpchap1n = dot(Wp1,n) + h1n*( Fchap1n - dot(F1,n) );
            Scalar Wpchap2n = dot(Wp2,n) + h2n*( Fchap2n - dot(F2,n) );

            /// test glissement adherence
            Point T = ((ProjT(Wp2,n) - ProjT(Wp1,n)) - (h2t*ProjT(F2,n) - h1t*ProjT(F1,n)) - dEp_t) / (h1t+h2t) - dPrecharge_t;
            Scalar g = f*std::abs(Fchap1n);

            Scalar normT = norm_2(T);
            Point Fchap1t,Fchap2t,Wpchap1t,Wpchap2t;

            //Attention au signe : modifie pour eviter la division par 0
            if (normT <= g) {
                /// adherence
                Fchap1t=T;
                Fchap2t=-1.0*Fchap1t;
                Wpchap1t= ProjT(Wp1,n) + h1t*(Fchap1t - ProjT(F1,n) );
                Wpchap2t= Wpchap1t;
// #ifdef LOOK_CONTACT_ZONE
//                 status[i*DIM]=1;
// endif
            } else if (normT > g) {
                /// glissement
                Fchap1t=T*g/normT;
                Fchap2t=-1.0*Fchap1t;
                Wpchap1t= ProjT(Wp1,n) + h1t*(Fchap1t - ProjT(F1,n) );
                Wpchap2t= ProjT(Wp2,n) + h2t*(Fchap2t - ProjT(F2,n) );
// #ifdef LOOK_CONTACT_ZONE
//                 status[i*DIM]=2;
// #endif
            }
            Wpchap1=Wpchap1n*n+Wpchap1t;
            Fchap1=Fchap1n*n+Fchap1t;
            Wpchap2=Wpchap2n*n+Wpchap2t;
            Fchap2=Fchap2n*n+Fchap2t;
        }

        ///integration
        Wchap1 = Wpchap1*dt + Wchap1old;
        Wchap2 = Wpchap2*dt + Wchap2old;
        WWchap1[rep]=Wchap1;
        WWchap2[rep]=Wchap2;

        Qchap1[rep]=Fchap1;
        Qchap2[rep]=Fchap2;
        WWpchap1[rep]=Wpchap1;
        WWpchap2[rep]=Wpchap2;

    }

    Inter.side[0].t[imic].Wpchap[list1]=WWpchap1;
    Inter.side[1].t[imic].Wpchap[list2]=WWpchap2;
    //if (temps.type_de_calcul=="Qstat") {
    Inter.side[0].t[imic].Wchap[list1]=WWchap1;
    Inter.side[1].t[imic].Wchap[list2]=WWchap2;
    //} else if(temps.type_de_calcul=="stat") {
    //    Inter.side[0].t[imic].Wchap[list1]=WWpchap1;
    //    Inter.side[1].t[imic].Wchap[list2]=WWpchap2;
    //} else {
    //    cout << "comportements contact : type de calcul non implemente " << endl;
    //    assert(0);
    //}
    Inter.side[0].t[imic].Fchap[list1]=Qchap1;
    Inter.side[1].t[imic].Fchap[list2]=Qchap2;
    //*/
}



//interface interieure de type defaut de forme
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procedure pour les interfaces interieures de type defaut de forme
*/
void compt_contact_ep (Interface &Inter,TimeData &temps) {/*

    const int imic = temps.pt;
    const Scalar dt=temps.dt;

    LMT::Vec<unsigned> &list1=Inter.side[0].ddlcorresp;
    LMT::Vec<unsigned> &list2=Inter.side[1].ddlcorresp;
    
    Vector WWpchap1 = Inter.side[0].t[imic].Wpchap[list1];
    Vector WWpchap2 = Inter.side[1].t[imic].Wpchap[list2];
    Vector Qchap1 = Inter.side[0].t[imic].Fchap[list1];
    Vector Qchap2 = Inter.side[1].t[imic].Fchap[list2];
    const Vector &Q1 = Inter.side[0].t[imic].F[list1];
    const Vector &Q2 = Inter.side[1].t[imic].F[list2];
    const Vector &WWp1 = Inter.side[0].t[imic].Wp[list1];
    const Vector &WWp2 = Inter.side[1].t[imic].Wp[list2];
    const Vector &WWchap1old = Inter.side[0].t[imic-1].Wchap[list1];
    const Vector &WWchap2old = Inter.side[1].t[imic-1].Wchap[list2];
    Vector WWchap1 = Inter.side[0].t[imic].Wchap[list1];
    Vector WWchap2 = Inter.side[1].t[imic].Wchap[list2];
    
    const Vector &JJ = Inter.jeu[list1];
    const Vector &neq = Inter.side[1].neq[list1];
    
    Scalar f = Inter.coeffrottement;
    const Vector &frott = Inter.coeffrottement_vec;
    
    Scalar h1n = Inter.side[0].hn;
    Scalar h2n = Inter.side[1].hn;
    Scalar h1t = Inter.side[0].ht;
    Scalar h2t = Inter.side[1].ht;

// #ifdef LOOK_CONTACT_ZONE
//     Vec<unsigned> status;status.resize(nbpts*INTER::dim);
// #endif

    /// travail pt par pt
    const unsigned nbpts=Inter.side[0].nodeeq.size();
    for(unsigned i=0;i<nbpts;++i) {
        /// Recuperation des indices de stockage
        Vec<unsigned> rep=range(i*DIM,(i+1)*DIM);
        /// Recuperation des vecteurs locaux
        Point n = -1.0*neq[rep];            /// Normale de 2 vers 1
        Point F1 = Q1[rep];                 /// Efforts du solide 1 sur l'interface calcules a l'etape lineaire
        Point F2 = Q2[rep];                 /// Efforts du solide 2 sur l'interface calcules a l'etape lineaire
        Point Fchap1 = Qchap1[rep];         /// Efforts du solide 1 sur l'interface calcules a l'etape locale precedente
        Point Fchap2 = Qchap2[rep];         /// Efforts du solide 2 sur l'interface calcules a l'etape locale precedente
        Point Wp1 = WWp1[rep];              /// Vitesse du solide 1 calculee a l'etape lineaire
        Point Wp2 = WWp2[rep];              /// Vitesse du solide 2 calculee a l'etape lineaire
        Point Wpchap1 = WWpchap1[rep];      /// Vitesse du solide 1 calculee a l'etape locale precedente
        Point Wpchap2 = WWpchap2[rep];      /// Vitesse du solide 2 calculee a l'etape locale precedente
        Point Wchap1old = WWchap1old[rep];  /// Deplacement du solide 1 calculee a la derniere etape locale du pas de temps precedent
        Point Wchap2old = WWchap2old[rep];  /// Deplacement du solide 2 calculee a la derniere etape locale du pas de temps precedent
        Point Wchap1 = WWchap1[rep];        /// Deplacement du solide 1 calculee a l'etape locale precedente
        Point Wchap2 = WWchap2[rep];        /// Deplacement du solide 2 calculee a l'etape locale precedente
        Point jeu=JJ[rep];                  /// Jeu entre les deux cotes de l'interface calcule a l'etape locale precedente

        /// test contact normal
        Scalar N=0;
        //if (temps.type_de_calcul=="stat")
        //    N = ( dot(n,(Wp2-Wp1)) -h2n*dot(n,F2) + h1n*dot(n,F1) )/(h2n+h1n);
        //else
        Wp1 += jeu;
        Wchap1old += jeu;
        N = ( dot(n,Wp2) - dot(n,Wp1) + (dot(n,Wchap2old) - dot(n,Wchap1old))/dt -h2n*dot(n,F2) + h1n*dot(n,F1) )/(h2n+h1n);

        if (N > 0.0) {
            /// separation de l'interface
            Fchap1 = 0.0;
            Fchap2 = 0.0;
            Wpchap1 = Wp1-h1n*dot(n,F1)*n - h1t* ProjT(F1,n);
            Wpchap2 = Wp2-h2n*dot(n,F2)*n - h2t* ProjT(F2,n);
//#ifdef LOOK_CONTACT_ZONE
//            status[i*DIM]=0;
//#endif
        } else {
            /// contact
            Scalar Fchap1n = N;
            Scalar Fchap2n = -1.0*Fchap1n;
            Scalar Wpchap1n = dot(Wp1,n) + h1n*( Fchap1n - dot(n,F1) );
            Scalar Wpchap2n = dot(Wp2,n) + h2n*( Fchap2n - dot(n,F2) );

            /// test glissement adherence
            Point Test = (ProjT(Wp2,n)-ProjT(Wp1,n) -h2t*ProjT(F2,n) + h1t*ProjT(F1,n))/(h2t+h1t);
            Scalar normT = norm_2(Test);
            
            Scalar g = frott[i]*std::abs(Fchap1n);
//            Scalar g = f*std::abs(Fchap1n);
            Point Fchap1t,Fchap2t,Wpchap1t,Wpchap2t;

            ///Attention au signe : modifie pour eviter la division par 0
            if (normT <= g) {
                /// adherence
                Fchap1t = Test;
                Fchap2t = -1.0*Fchap1t;
                Wpchap1t = ProjT(Wp1,n) + h1t*(Fchap1t - ProjT(F1,n) );
                Wpchap2t = Wpchap1t;
//#ifdef LOOK_CONTACT_ZONE
//                status[i*INTER::dim]=1;
//#endif
            } else if (normT > g) {
                /// glissement
                Fchap1t = Test*g/normT;
                Fchap2t = -1.0*Fchap1t;
                Wpchap1t = ProjT(Wp1,n) + h1t*(Fchap1t - ProjT(F1,n) );
                Wpchap2t = ProjT(Wp2,n) + h2t*(Fchap2t - ProjT(F2,n) );
//#ifdef LOOK_CONTACT_ZONE
//                status[i*INTER::dim]=2;
//#endif
            }

            Wpchap1 = Wpchap1n*n + Wpchap1t - jeu;
            Fchap1 = Fchap1n*n+Fchap1t;
            Wpchap2 = Wpchap2n*n+Wpchap2t;
            Fchap2 = Fchap2n*n+Fchap2t;
        }

        ///integration
        Wchap1 = Wpchap1*dt + (Wchap1old-jeu);
        Wchap2 = Wpchap2*dt + Wchap2old;
        WWchap1[rep] = Wchap1;
        WWchap2[rep] = Wchap2;

        Qchap1[rep] = Fchap1;
        Qchap2[rep] = Fchap2;
        WWpchap1[rep] = Wpchap1;
        WWpchap2[rep] = Wpchap2;

    }

    Inter.side[0].t[imic].Wpchap[list1] = WWpchap1;
    Inter.side[1].t[imic].Wpchap[list2] = WWpchap2;
    Inter.side[0].t[imic].Wchap[list1] = WWchap1;
    Inter.side[1].t[imic].Wchap[list2] = WWchap2;
    Inter.side[0].t[imic].Fchap[list1] = Qchap1;
    Inter.side[1].t[imic].Fchap[list2] = Qchap2;
*/}



//interface interieure de type parfait cassable => contact avec frottement
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procedure pour les interfaces interieures de type parfait cassable
*/
void compt_breakable (Interface &Inter,TimeData &temps) {

    const int imic = temps.pt;
    const Scalar dt = temps.dt;
    
    Interface::Side &side_0 = Inter.side[0];
    Interface::Side &side_1 = Inter.side[1];

    const LMT::Vec<unsigned> &list1 = side_0.ddlcorresp;
    const LMT::Vec<unsigned> &list2 = side_1.ddlcorresp;
    
    const Vector &Q1           = side_0.t[imic].F[list1];
    const Vector &Q2           = side_1.t[imic].F[list2];
    const Vector &WWp1         = side_0.t[imic].Wp[list1];
    const Vector &WWp2         = side_1.t[imic].Wp[list2];
    const Vector &WWchap1old   = side_0.t[imic-1].Wchap[list1];
    const Vector &WWchap2old   = side_1.t[imic-1].Wchap[list2];
    Vector WWchap1             = side_0.t[imic].Wchap[list1];
    Vector WWchap2             = side_1.t[imic].Wchap[list2];
    Vector WWpchap1            = side_0.t[imic].Wpchap[list1];
    Vector WWpchap2            = side_1.t[imic].Wpchap[list2];
    Vector Qchap1              = side_0.t[imic].Fchap[list1];
    Vector Qchap2              = side_1.t[imic].Fchap[list2];
    
    const Vector &neq = side_0.neq[list1];
    
    const Scalar f = Inter.coeffrottement;
    
    const Scalar h1n = side_0.hn;
    const Scalar h2n = side_1.hn;
    const Scalar h1t = side_0.ht;
    const Scalar h2t = side_1.ht;
    
    // travail pt par pt
    const unsigned nbpts = side_0.nodeeq.size();
    for(unsigned i=0;i<nbpts;++i) {
        const Vec<unsigned> rep=range(i*DIM,(i+1)*DIM);
        Point n         = neq[rep]; /// normale de 1 vers 2
        Point F1        = Q1[rep];
        Point F2        = Q2[rep];
        Point Fchap1    = Qchap1[rep];
        Point Fchap2    = Qchap2[rep];
        Point Wp1       = WWp1[rep];
        Point Wp2       = WWp2[rep];
        Point Wpchap1   = WWpchap1[rep];
        Point Wpchap2   = WWpchap2[rep];
        Point Wchap1old = WWchap1old[rep];
        Point Wchap2old = WWchap2old[rep];
        Point Wchap1    = WWchap1[rep];
        Point Wchap2    = WWchap2[rep];

        /// test contact normal
        Scalar N = ( dot(n,(Wp2-Wp1)) + dot(n,(Wchap2old-Wchap1old)/dt) -h2n*dot(n,F2) + h1n*dot(n,F1) )/(h2n+h1n);
        
        /// si la convergence du calcul iteratif est OK, on met à jour le comportement des elements qui ne sont pas deja casse
        if (Inter.convergence >= 0 and Inter.comportement[i] == false){
            if (N > 0.0 and std::abs(dot(n,F1)) > Inter.matprop->Fcr_n){
                ///David dit de mettre 10% de plus
                Inter.comportement[i] = true;
                Inter.convergence++;
            }
        }
        
        if (Inter.comportement[i] == 1){
            ///interface cassee...
            if(N > 0.0) {
                /// separation de l'interface...
                Fchap1  = 0.0;
                Fchap2  = 0.0;
                Wpchap1 = Wp1-h1n*dot(n,F1)*n - h1t* ProjT(F1,n);
                Wpchap2 = Wp2-h2n*dot(n,F2)*n - h2t* ProjT(F2,n);
            } else {
                /// contact...
                Scalar Fchap1n  = N;
                Scalar Fchap2n  = -1.0*Fchap1n;
                Scalar Wpchap1n = dot(Wp1,n) + h1n*( Fchap1n - dot(n,F1) );
                Scalar Wpchap2n = dot(Wp2,n) + h2n*( Fchap2n - dot(n,F2) );

                /// test de glissement adherence
                Point T;
                T=(ProjT(Wp2,n)-ProjT(Wp1,n) -h2t*ProjT(F2,n) + h1t*ProjT(F1,n))/(h2t+h1t);
                Scalar g = f*std::abs(Fchap1n);

                Scalar normT = norm_2(T);
                Point Fchap1t,Fchap2t,Wpchap1t,Wpchap2t;

                ///Attention au signe : modifie pour eviter la division par 0
                if (normT <= g) {
                    /// adherence...
                    Fchap1t  = T;
                    Fchap2t  = -1.0*Fchap1t;
                    Wpchap1t = ProjT(Wp1,n) + h1t*(Fchap1t - ProjT(F1,n) );
                    Wpchap2t = Wpchap1t;
                } else if (normT > g) {
                    /// glissement...
                    Fchap1t  = T*g/normT;
                    Fchap2t  = -1.0*Fchap1t;
                    Wpchap1t = ProjT(Wp1,n) + h1t*(Fchap1t - ProjT(F1,n) );
                    Wpchap2t = ProjT(Wp2,n) + h2t*(Fchap2t - ProjT(F2,n) );
                }
                /// assemblage des composantes normales et tangentielles
                Wpchap1 = Wpchap1n*n+Wpchap1t;
                Fchap1  = Fchap1n*n+Fchap1t;
                Wpchap2 = Wpchap2n*n+Wpchap2t;
                Fchap2  = Fchap2n*n+Fchap2t;
            }
            
            ///integration
            Wchap1       = Wpchap1*dt + Wchap1old;
            Wchap2       = Wpchap2*dt + Wchap2old;
            WWchap1[rep] = Wchap1;
            WWchap2[rep] = Wchap2;

            Qchap1[rep]   = Fchap1;
            Qchap2[rep]   = Fchap2;
            WWpchap1[rep] = Wpchap1;
            WWpchap2[rep] = Wpchap2;
            
        } else {
            /// interface parfaite...
            /// creation des operateurs locaux de direction de recherche
            Kloc kloc1;
            kloc1.kn = side_0.kn;
            kloc1.kt = side_0.kt;
            kloc1.n = n;
            
            Kloc kloc2;
            kloc2.kn = side_1.kn;
            kloc2.kt = side_1.kt;
            kloc2.n = n;
            
            Kloc hloc;
            hloc.kn = 1.0/(side_1.kn+side_0.kn);
            hloc.kt = 1.0/(side_1.kt+side_0.kt);
            hloc.n  = n;
            
            /// travail point par point
            WWpchap1[rep] = hloc*(kloc1*Wp1+kloc2*Wp2-(F1+F2));
            WWpchap2[rep] = hloc*(kloc1*Wp1+kloc2*Wp2-(F1+F2));
            Qchap1[rep]   = F1+kloc1*(Wpchap1[rep]-Wp1);
            Qchap2[rep]   = -1.0*Qchap1[rep];
            
            /// integration
            Wchap1       = WWpchap1[rep]*dt + Wchap1old;
            Wchap2       = WWpchap1[rep]*dt + Wchap2old;
            WWchap1[rep] = Wchap1;
            WWchap2[rep] = Wchap2;
        }
        
        /// Stockage des resultats
        side_0.t[imic].Wpchap[list1] = WWpchap1;
        side_1.t[imic].Wpchap[list2] = WWpchap2;
        side_0.t[imic].Wchap[list1]  = WWchap1;
        side_1.t[imic].Wchap[list2]  = WWchap2;
        side_0.t[imic].Fchap[list1]  = Qchap1;
        side_1.t[imic].Fchap[list2]  = Qchap2;
    }
}










struct assign_d_mesh {
    template<class TE> void operator() (TE &e,unsigned &incr, Vector &d) const {
        e.d=d[incr];
        incr+=1;
    }
};


//interface interieure avec loi cohesive (type mesomodele)
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procedure pour les interfaces interieures de type cohesif
*/
void compt_cohesif (Interface &Inter,TimeData &temps) {/*

    int imic = temps.pt;
    Scalar dt=temps.dt;

    /// Correspondance des noeuds sur l'interface
    const Vec<unsigned> &list1 = Inter.side[0].ddlcorresp;
    const Vec<unsigned> &list2 = Inter.side[1].ddlcorresp;
    /// Champs connus
    //const Vector &FF1=Inter.side[0].t[imic].F[list1];
    //const Vector &FF2=Inter.side[1].t[imic].F[list2];
    const Vector &WWp1 = Inter.side[0].t[imic].Wp[list1];
    const Vector &WWp2 = Inter.side[1].t[imic].Wp[list2];
    const Vector &WWchap1old = Inter.side[0].t[imic-1].Wchap[list1];
    const Vector &WWchap2old = Inter.side[1].t[imic-1].Wchap[list2];
    /// Champs recherches
    Vector WWchap1 = Inter.side[0].t[imic].Wchap[list1];
    Vector WWchap2 = Inter.side[1].t[imic].Wchap[list2];
    Vector FFchap1 = Inter.side[0].t[imic].Fchap[list1];
    Vector FFchap2 = Inter.side[1].t[imic].Fchap[list2];

    /// Vecteur des normales
    const Vector &neq = Inter.side[0].neq[list1];

    /// Matrice identite
    DenseMatrix Id;
    Id.set(0.);
    Id.diag() = Point(1);

    //Scalar eps=1e-3;

    /// parametres de direction de recherche
    //Scalar h1n=Inter.side[0].hn;
    //Scalar h2n=Inter.side[1].hn;
    //Scalar h1t=Inter.side[0].ht;
    //Scalar h2t=Inter.side[1].ht;
    
    bool contact = false;

    /// travail pt par pt
    const unsigned nbpts=Inter.side[0].nodeeq.size();
    for(unsigned i=0;i<nbpts;++i) {
        const Vec<unsigned> rep = range(i*DIM,(i+1)*DIM);   /// Vecteur des numeros des noeuds de l'element
        const Point n = neq[rep];               /// normale de 1 vers 2
        ///matrice normale
        DenseMatrix ntn;
        tens(n,n,ntn);
        ///matrices des directions de recherche
        //DenseMatrix H1 = h1n*ntn + h1t*(Id-ntn);
        //DenseMatrix H2 = h2n*ntn + h2t*(Id-ntn);
        
        ///variables connues
        //Point F1=FF1[rep];
        //Point F2=FF2[rep];
        const Point Wp1 = WWp1[rep];
        const Point Wp2 = WWp2[rep];
        const Point Wchap1old = WWchap1old[rep];
        const Point Wchap2old = WWchap2old[rep];
        const Scalar dold=Inter.t[imic-1].d[i];               // TODO verification
        
        ///variables inconnues
        //Point Fchap1(0.);
        //Point Fchap2(0.);
        //Point Fchap1t(0.);
        //Point Wchap1=WWchap1[rep];
        //Point Wchap2=WWchap2[rep];
        const Point Wchap1 = Wp1 * dt + Wchap1old;
        const Point Wchap2 = Wp2 * dt + Wchap2old;
        
        bool compression = dot(n,(Wchap2-Wchap1)) < 0 ;
        
        const Point jumpWchap = Wchap2 - Wchap1;
        const Scalar jumpn = dot(jumpWchap,n);
        const Point jumpt = jumpWchap - jumpn * n;
        const Scalar Yn = 0.5 * Inter.matprop->kn * jumpn * jumpn * (jumpn > 0);
        const Scalar Yt = 0.5 * Inter.matprop->kt * norm_2(jumpt) * norm_2(jumpt);
        const Scalar Y = std::pow( std::pow(Yn,Inter.matprop->alpha) + (Inter.matprop->gamma*std::pow(Yt,Inter.matprop->alpha)) , 1./Inter.matprop->alpha) ;
        const Scalar w = (Inter.matprop->n/(Inter.matprop->n+1))* std::pow( ((Y-Inter.matprop->Yo)>=0)*(Y-Inter.matprop->Yo)/(Inter.matprop->Yc-Inter.matprop->Yo), Inter.matprop->n ) ;
        const Scalar d = std::min( std::max(w,dold), 1.0);
        
        DenseMatrix Kc = Inter.matprop->knc*ntn + Inter.matprop->kt*(1.-d)*(Id - ntn);         /// comportement en compression
        DenseMatrix Kt = Inter.matprop->kn*(1.-d)*ntn +  Inter.matprop->kt*(1.-d)*(Id - ntn);  /// comportement en traction
        DenseMatrix Ktilde = Kt * (!compression) + Kc * compression ;
        
        Point Fchap1 = Ktilde * jumpWchap;
        Point Fchap2 = -Fchap1;
        
        Inter.t[imic].d[i] = d ;
        
        /// affectation des champs par elements aux champs par interface
        FFchap1[rep] = Fchap1 ;
        FFchap2[rep] = Fchap2 ;
        WWchap1[rep] = Wchap1 ;
        WWchap2[rep] = Wchap2 ;
    }
    
    /// affectation des champs chapeaux calcules
    Inter.side[0].t[imic].Wchap[list1] = WWchap1 ;
    Inter.side[1].t[imic].Wchap[list2] = WWchap2 ;
    Inter.side[0].t[imic].Fchap[list1] = FFchap1 ;
    Inter.side[1].t[imic].Fchap[list2] = FFchap2 ;
    
    /// calcul et affectation des champs de vitesses
    Inter.side[0].t[imic].Wpchap = (Inter.side[0].t[imic].Wchap - Inter.side[0].t[imic-1].Wchap)/dt;
    Inter.side[1].t[imic].Wpchap = (Inter.side[1].t[imic].Wchap - Inter.side[1].t[imic-1].Wchap)/dt;

        ///initialisation de la boucle iterative
        Scalar Ymax = Inter.param_comp->Y[i];
        Scalar d = dold;
        Scalar erreur = 1.;
        unsigned it=0;
        //boucle iterative
        while(erreur>eps) {
            it+=1;

            if(dold>=1.-eps) {
                ///etape necessaire lorsque d vaut 1 dans la boucle iterative : on ne gere pas le frottement tant que toute l'interface n'est pas cassee
                Scalar N=0;
                //if (temps.type_de_calcul=="stat")
                //    N = ( dot(n,(Wp2-Wp1)) -h2n*dot(n,F2) + h1n*dot(n,F1) )/(h2n+h1n);
                //else
                    N=( dot(n,(Wp2-Wp1)) + dot(n,(Wchap2old-Wchap1old)/dt) -h2n*dot(n,F2) + h1n*dot(n,F1) )/(h2n+h1n);

                if (N>0.0) { /// separation de l'interface
                    Fchap1=0.0;
                    Fchap2=0.0;
                    //Wchap1=(Wp1-h1n*dot(n,F1)*n - h1t* ProjT(F1,n))*dt+Wchap1old;
                    //Wchap2=(Wp2-h2n*dot(n,F2)*n - h2t* ProjT(F2,n))*dt+Wchap2old;
                } else { /// contact
                    Fchap1 = N*n;
                    Fchap2 = -1.*Fchap1;
                }
                
                //if(temps.type_de_calcul=="Qstat"){
                  Wchap1=(Wp1 + H1*(Fchap1 - F1) )*dt + Wchap1old;
                  Wchap2=(Wp2 + H2*(Fchap2 - F2) )*dt + Wchap2old;
                //}
                //else if(temps.type_de_calcul=="stat"){
                //  Wchap1=(Wp1 + H1*(Fchap1 - F1) );
                //  Wchap2=(Wp2 + H2*(Fchap2 - F2) );                
                //}
                //else{cout << "Mauvais type de calcul" << endl; assert(0);}
                
                erreur=0;///pour sortir de la boucle
            } else {///dans tous les autres cas
                ///matrice pour relation de comportement d'interface prenant en compte l'endommagement
                TMAT Kc;
                Kc = damage.knc*ntn +  damage.kt*(1.-dold)*(Id - ntn);
                TMAT Kt;
                Kt = damage.kn*(1.-dold)*ntn +  damage.kt*(1.-dold)*(Id - ntn);
                TMAT invKc = 1./damage.knc*ntn + 1./(damage.kt*(1.-dold))*(Id - ntn);
                TMAT invKt = 1./(damage.kn*(1.-dold))*ntn +  1./(damage.kt*(1.-dold))*(Id - ntn);

                ///1ere etape : calcul de Fchap1
                Point Wtemp;
                //if(temps.type_de_calcul=="Qstat")
                Wtemp= ( -1.*H2*F2 + H1*F1 + (Wp2-Wp1) + (Wchap2old-Wchap1old)/dt );
                //else if(temps.type_de_calcul=="stat")
                //  Wtemp= ( -1.*H2*F2 + H1*F1 + (Wp2-Wp1)  );
                //else{cout << "Mauvais type de calcul" << endl; assert(0);}
                
                Fchap1 = inv(H2+H1+invKc)*Wtemp;
                if (dot(Fchap1,n)>=0) {
                    Fchap1 = inv(H2+H1+invKt)*Wtemp;
                }
                ///calcul de Fchap2, Wchap1, Wchap2
                Fchap2 = -1.*Fchap1;
                //if(temps.type_de_calcul=="Qstat"){
                Wchap1=(Wp1 + H1*(Fchap1 - F1 ))*dt + Wchap1old;
                Wchap2=(Wp2 + H2*(Fchap2 - F2) )*dt + Wchap2old;
                //}
                //else if(temps.type_de_calcul=="stat"){
                //  Wchap1=(Wp1 + H1*(Fchap1 - F1) );
                //  Wchap2=(Wp2 + H2*(Fchap2 - F2) );                
                //}
                //else{cout << "Mauvais type de calcul" << endl; assert(0);}
                
                ///calcul des forces d'endommagement normales et tangentielles a partir de Fchap1
                Scalar Ydn = (dot(Fchap1,n)>=0)*( std::pow((dot(Fchap1,n)/(1.-d)),2.)/(2.*damage.kn) );
                Fchap1t = Fchap1 - dot(Fchap1,n)*n;
                Scalar Yt = ( std::pow(norm_2(Fchap1t)/(1.-d),2.)/(2.*damage.kt) );
                Scalar Ynew = std::pow( std::pow(Ydn,damage.alpha)+ (damage.gamma*std::pow(Yt,damage.alpha)) , 1./damage.alpha);
                Ymax = std::max(Ynew,Inter.param_comp->Ymax[i]);
                ///calcul du nouvel endommagement pour cette iteration
                Scalar Ypos = ((Ymax-damage.Yo)>=0)*(Ymax-damage.Yo);
                //        cout<< "ftemp " << Ftemp << " nodeeq "<<Inter.side[0].nodeeq[i]<< endl;
                //         cout<< " n " <<n <<"F1 " << F1 << " F2 " << F2 << " nodeeq "<<Inter.side[0].nodeeq[i]<< endl;
                //      cout << "Appuyez sur entrÃ©e pour continuer...";
                //        cin.ignore( numeric_limits<streamsize>::max(), '\n' );

                d = std::pow((damage.n/(damage.n+1))*(Ypos/(damage.Yc-damage.Yo)),damage.n);
                d = std::min(d,1.);
                if (d>=1.-eps) {
                    ///test contact (sans frottement)
                    Scalar N=0;
                    //if (temps.type_de_calcul=="stat")
                    //    N = ( dot(n,(Wp2-Wp1)) -h2n*dot(n,F2) + h1n*dot(n,F1) )/(h2n+h1n);
                    //else
                        N=( dot(n,(Wp2-Wp1)) + dot(n,(Wchap2old-Wchap1old)/dt) -h2n*dot(n,F2) + h1n*dot(n,F1) )/(h2n+h1n);

                    if (N>0.0) { // separation de l'interface
                        Fchap1=0.0;
                        Fchap2=0.0;
                        //Wchap1=(Wp1-h1n*dot(n,F1)*n - h1t* ProjT(F1,n))*dt+Wchap1old;
                        //Wchap2=(Wp2-h2n*dot(n,F2)*n - h2t* ProjT(F2,n))*dt+Wchap2old;
                    } else { // contact
                        Fchap1 = N*n;
                        Fchap2 = -1.*Fchap1;
                    }
                    
                     //if(temps.type_de_calcul=="Qstat"){
                        Wchap1=(Wp1 + H1*(Fchap1 - F1) )*dt + Wchap1old;
                        Wchap2=(Wp2 + H2*(Fchap2 - F2) )*dt + Wchap2old;
                     //}
                     //else if(temps.type_de_calcul=="stat"){
                     //   Wchap1=(Wp1 + H1*(Fchap1 - F1) );
                     //   Wchap2=(Wp2 + H2*(Fchap2 - F2) );                
                     //}
                     //else{cout << "Mauvais type de calcul" << endl; assert(0);}
                    //Wchap1=(Wp1 + H1*(Fchap1 - F1) )*dt + Wchap1old;
                    //Wchap2=(Wp2 + H2*(Fchap2 - F2) )*dt + Wchap2old;
                    
                    dold=d;
                } else {
                    erreur=std::abs(d-dold)/std::abs(dold+1e-10);
                    dold=d;
                }
            }
        }
        Inter.param_comp->t[imic].d[i]=d;
         
        //cout<< it << endl;
        ///modification dmax pour le pas de temps suivant si necessaire
        if(d>=Inter.param_comp->t[imic].dmax[i]) {
            Inter.param_comp->t[imic].dmax[i]=d;
            Inter.param_comp->t[imic].Ymax[i]=Ymax;
        }

        FFchap1[rep]=Fchap1;
        FFchap2[rep]=Fchap2;
        WWchap1[rep]=Wchap1;
        WWchap2[rep]=Wchap2;
    }

    Inter.side[0].t[imic].Wchap[list1]=WWchap1;
    Inter.side[1].t[imic].Wchap[list2]=WWchap2;
    Inter.side[0].t[imic].Fchap[list1]=FFchap1;
    Inter.side[1].t[imic].Fchap[list2]=FFchap2;
    
    //cout << Inter.param_comp->t[imic].d<<" d et Y " << Inter.param_comp->t[imic].Y <<" et F " << Inter.side[0].t[imic].Fchap << " sautW " << Inter.side[1].t[imic].Wchap - Inter.side[0].t[imic].Wchap<< endl;
    
    Inter.side[0].t[imic].Wpchap=(Inter.side[0].t[imic].Wchap - Inter.side[0].t[imic-1].Wchap)/dt;
    Inter.side[1].t[imic].Wpchap=(Inter.side[1].t[imic].Wchap - Inter.side[1].t[imic-1].Wchap)/dt;
    unsigned incr=0;
    apply(Inter.side[0].mesh.elem_list,assign_d_mesh(),incr,Inter.param_comp->t[imic].d);
    incr=0;
    apply(Inter.side[1].mesh.elem_list,assign_d_mesh(),incr,Inter.param_comp->t[imic].d);
   
    ///modification du type d'interface si dmax = 1 pour tous les points
    TT dmin;
    for(unsigned i=0;i<nbpts;i++) {
       dmin = min(dmin,Inter.param_comp->d[i]);
    }
    if (dmin>=1.-eps) {
       cout <<Inter.num <<  " Inter modif" << endl;
       Inter.comp="Contact";
    }  */
}

#endif // COMPORTEMENT_INTERFACES_H
