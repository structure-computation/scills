
//interface exterieure de type deplacement impose pour tous les ddls
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procédure pour les interfaces extérieures à déplacement imposé (ou vitesse)
*/
template<class INTER>
void compt_CL_depl (INTER &Inter,int &imic) {
    Inter.side[0].t[imic].Fchap=Inter.side[0].t[imic].F + Inter.side[0].kglo*(Inter.side[0].t[imic].Wpchap - Inter.side[0].t[imic].Wp);
}

//interface exterieure de type effort impose pour tous les ddls
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procédure pour les interfaces extérieures à effort imposé
*/
template<class INTER>
void compt_CL_eff (INTER &Inter,int &imic) {
    Inter.side[0].t[imic].Wpchap=Inter.side[0].t[imic].Wp+Inter.side[0].hglo*(Inter.side[0].t[imic].Fchap - Inter.side[0].t[imic].F);
}

//interface exterieure de type : deplacement normal a l'interface plane impose + effort tangentiel nul
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procédure pour les interfaces extérieures de type symétrie ou ayant seulement le déplacement normal donné (ou vitesse)
*/
template<class INTER>
void compt_CL_sym (INTER &Inter,int &imic) {
    typedef typename INTER::T TT;
    Vec<TT> &WWchap1=Inter.side[0].t[imic].Wpchap;
    Vec<TT> &QQchap1=Inter.side[0].t[imic].Fchap;
    const Vec<TT> &Q1=Inter.side[0].t[imic].F;
    const Vec<TT> &WW1=Inter.side[0].t[imic].Wp;

    TT ht = Inter.side[0].ht;
    TT kt = Inter.side[0].kt;
    TT kn = Inter.side[0].kn;

    Vec<TT> neq=Inter.side[0].neq;
    unsigned dim=INTER::dim;

    for(unsigned i=0;i<Inter.side[0].nodeeq.size();++i) {
        Vec<unsigned> rep=range(i*dim,(i+1)*dim);
        Vec<TT,INTER::dim> n=neq[rep]; // normale de 1 vers 2
        Vec<TT,INTER::dim> F1=Q1[rep], Fchap1=QQchap1[rep], W1=WW1[rep], Wchap1=WWchap1[rep], PtWchap1, PnFchap1 ;

        TT Wchap1n = dot(Wchap1,n);//donnee = 0 symetrie ou !=0 pour depl normal
        Wchap1.set(0.);
        Wchap1 = Wchap1n*n + ht * ( kt * ProjT(W1,n) - ProjT(F1,n) );
        Fchap1.set(0.);
        Fchap1 = dot(F1,n)*n + kn * (Wchap1n*n - dot(W1,n)*n);

        WWchap1[rep] = Wchap1;
        QQchap1[rep] = Fchap1;
    }

    Inter.side[0].t[imic].Wpchap=WWchap1;
    Inter.side[0].t[imic].Fchap=QQchap1;
}


//structure permettant de définir l'operateur de direction de recherche local a partir de direction normale et tangentielle
// s'utilise ensuite comme une matrice
template<class T,unsigned s>
struct Kloc{
   T kn,kt;
   Vec<T,s,void> n;
   template<class TV> Vec<T,s,void> operator*(TV &W){return kn*n*dot(n,W)+kt*(W-n*dot(n,W));}
   template<class TV> Vec<T,s,void> operator*(const TV &W){return kn*n*dot(n,W)+kt*(W-n*dot(n,W));}
};


//interface interieure de type parfaite
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procédure pour les interfaces intérieures parfaites
*/
template<class INTER>
void compt_parfait (INTER &Inter,int &imic) {
    typedef typename INTER::T T;
    typedef Mat <T , Gen<>, SparseLine<> > TMAT;
    Vec<unsigned> &list1=(Inter.side[0].ddlcorresp);
    Vec<unsigned> &list2=(Inter.side[1].ddlcorresp);
    Vec<T> Wchap1=Inter.side[0].t[imic].Wpchap[list1];
    Vec<T> Wchap2=Inter.side[1].t[imic].Wpchap[list2];
    Vec<T> Fchap1=Inter.side[0].t[imic].Fchap[list1];
    Vec<T> Fchap2=Inter.side[1].t[imic].Fchap[list2];
    const Vec<T> &Q1=Inter.side[0].t[imic].F[list1];
    const Vec<T> &Q2=Inter.side[1].t[imic].F[list2];
    const Vec<T> &WW1=Inter.side[0].t[imic].Wp[list1];
    const Vec<T> &WW2=Inter.side[1].t[imic].Wp[list2];
    const Vec<T> &neq1=(Inter.side[0].neq)[list1];
    const Vec<T> &JJ=Inter.param_comp->jeu[list1];
    //const Vec<T> &neq2=(Inter.side[1].neq)[list2];

    //creation des operateurs locaux de direction de recherche
    Kloc<T,INTER::dim> kloc1;
    kloc1.kn=Inter.side[0].kn;
    kloc1.kt=Inter.side[0].kt;
    
    Kloc<T,INTER::dim> kloc2;
    kloc2.kn=Inter.side[1].kn;
    kloc2.kt=Inter.side[1].kt;
    
    Kloc<T,INTER::dim> hloc;
    hloc.kn=1./(Inter.side[1].kn+Inter.side[0].kn);
    hloc.kt=1./(Inter.side[1].kt+Inter.side[0].kt);
    
    //travail point par point
    for(unsigned i=0;i<Inter.side[0].nodeeq.size();i++) {
        //creation du reperage du point
        Vec<unsigned> rep=range(i*INTER::dim,(i+1)*INTER::dim);
        //assignation de la normale aux operateurs elementaires
        Vec<T,INTER::dim> n1=neq1[rep];
        kloc1.n=n1;
        kloc2.n=n1;
        hloc.n=n1;
        // travail point par point
        Vec<T,INTER::dim> F1=Q1[rep],F2=Q2[rep], W1=WW1[rep],W2=WW2[rep] ,jeu=JJ[rep]  ;
        W1 += jeu;
        Wchap1[rep]=hloc*(kloc1*W1+kloc2*W2-(F1+F2));
        Fchap1[rep]=F1+kloc1*(Wchap1[rep]-W1);
    }

    Inter.side[0].t[imic].Wpchap[list1]=Wchap1;
    Inter.side[1].t[imic].Wpchap[list2]=Wchap1;
    Inter.side[0].t[imic].Fchap[list1]=Fchap1;
    Inter.side[1].t[imic].Fchap[list2]=-1.0*Fchap1;
    
    
/*    cout << WW1 << endl;
    cout << WW2 << endl;
    cout << Wchap1 << endl;
    cout << Q1 << endl;
    cout << Q2 << endl;
    cout << Fchap1 << endl;*/
   
    
}

//interface interieure de type jeu impose = idem interface parfaite avec prise en compte du jeu
//comme il y a des tests stat Qstat, on a préféré séparé dans un premier temps
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procédure pour les interfaces intérieures de type jeu imposé
*/
template<class INTER>
void compt_jeu_impose (INTER &Inter,TEMPS &temps) {
    int imic = temps.pt;
    unsigned pt_cur=temps.pt_cur;
    typedef typename INTER::T T;
    typedef Mat <T , Gen<>, SparseLine<> > TMAT;
    Vec<unsigned> &list1=(Inter.side[0].ddlcorresp);
    Vec<unsigned> &list2=(Inter.side[1].ddlcorresp);
    Vec<T> Wchap1=Inter.side[0].t[imic].Wpchap[list1];
    Vec<T> Wchap2=Inter.side[1].t[imic].Wpchap[list2];
    Vec<T> Fchap1=Inter.side[0].t[imic].Fchap[list1];
    Vec<T> Fchap2=Inter.side[1].t[imic].Fchap[list2];
    const Vec<T> &Q1=Inter.side[0].t[imic].F[list1];
    const Vec<T> &Q2=Inter.side[1].t[imic].F[list2];
    const Vec<T> &WW1=Inter.side[0].t[imic].Wp[list1];
    const Vec<T> &WW2=Inter.side[1].t[imic].Wp[list2];
    const Vec<T> &JJ=Inter.param_comp->jeu[list1];
    const Vec<T> &neq1=(Inter.side[0].neq)[list1];
    //const Vec<T> &neq2=(Inter.side[1].neq)[list2];

    //creation des operateurs locaux de direction de recherche
    Kloc<T,INTER::dim> kloc1;
    kloc1.kn=Inter.side[0].kn;
    kloc1.kt=Inter.side[0].kt;
    
    Kloc<T,INTER::dim> kloc2;
    kloc2.kn=Inter.side[1].kn;
    kloc2.kt=Inter.side[1].kt;
    
    Kloc<T,INTER::dim> hloc;
    hloc.kn=1./(Inter.side[1].kn+Inter.side[0].kn);
    hloc.kt=1./(Inter.side[1].kt+Inter.side[0].kt);

    //travail point par point    
    for(unsigned i=0;i<Inter.side[0].nodeeq.size();i++) {
        //creation du reperage du point
        Vec<unsigned> rep=range(i*INTER::dim,(i+1)*INTER::dim);
        //creation des matrices elementaires de direction de recherche
        Vec<T,INTER::dim> n1=neq1[rep];
        kloc1.n=n1;
        kloc2.n=n1;
        hloc.n=n1;

        // travail point par point
        Vec<T,INTER::dim> F1=Q1[rep],F2=Q2[rep], W1=WW1[rep],W2=WW2[rep],jeu=JJ[rep]  ;
        if (pt_cur==1)
            Wchap1[rep]=hloc * ( kloc1*W1+kloc2*W2 -(F1+F2) -kloc2*jeu/temps.dt);
        else
            Wchap1[rep]=hloc * ( kloc1*W1+kloc2*W2 -(F1+F2));

        Fchap1[rep]=F1+kloc1*(Wchap1[rep]-W1);

    }

    if(pt_cur<=Inter.param_comp->nbpastempsimpos)
        Wchap2=Wchap1+JJ/(temps.dt*Inter.param_comp->nbpastempsimpos);
    else
        Wchap2=Wchap1;

    Inter.side[0].t[imic].Wpchap[list1]=Wchap1;
    Inter.side[1].t[imic].Wpchap[list2]=Wchap2;
    Inter.side[0].t[imic].Fchap[list1]=Fchap1;
    Inter.side[1].t[imic].Fchap[list2]=-1.0*Fchap1;


    //     int imic = temps.pt;
    //     typedef typename INTER::T T;
    //     Vec<unsigned> &list1=Inter.side[0].ddlcorresp;
    //     Vec<unsigned> &list2=Inter.side[1].ddlcorresp;
    //     Vec<T> Wchap1=Inter.side[0].t[imic].Wpchap[list1];
    //     Vec<T> Wchap2=Inter.side[1].t[imic].Wpchap[list2];
    //     Vec<T> Qchap1=Inter.side[0].t[imic].Fchap[list1];
    //     Vec<T> Qchap2=Inter.side[1].t[imic].Fchap[list2];
    //     const Vec<T> &Q1=Inter.side[0].t[imic].F[list1];
    //     const Vec<T> &Q2=Inter.side[1].t[imic].F[list2];
    //     const Vec<T> &W1=Inter.side[0].t[imic].Wp[list1];
    //     const Vec<T> &W2=Inter.side[1].t[imic].Wp[list2];
    //
    //     T k1,k2;
    //     k1 = Inter.side[0].kn;
    //     k2 = Inter.side[1].kn;
    //     T ktot=k1+k2;
    //     T iktot=1/ktot;
    //
    //     if (temps.type_de_calcul=="Qstat") {
    //         if (imic==1) {
    //             Wchap1=iktot * ( k1*W1+k2*W2 -(Q1+Q2)-k2*(Inter.param_comp->jeu[list1]/temps.dt));
    //             Wchap2=Wchap1+(Inter.param_comp->jeu[list1]/temps.dt);
    //         } else {
    //             Wchap1=iktot * ( k1*W1+k2*W2 -(Q1+Q2));
    //             Wchap2=Wchap1;
    //         }
    //     } else if (temps.type_de_calcul=="stat") {
    //         Wchap1=iktot * ( k1*W1+k2*W2 -(Q1+Q2) - k2*Inter.param_comp->jeu[list1]);
    //         Wchap2=Wchap1+Inter.param_comp->jeu[list1];
    //     } else {
    //         cout << "comportements jeu impose : type de calcul non implemente " << endl;
    //         assert(0);
    //     }
    //
    //     Qchap1=Q1+k1*(Wchap1-W1);
    //     Qchap2=-1.0*Qchap1;
    //
    //     Inter.side[0].t[imic].Wpchap[list1]=Wchap1;
    //     Inter.side[1].t[imic].Wpchap[list2]=Wchap2;
    //     Inter.side[0].t[imic].Fchap[list1]=Qchap1;
    //     Inter.side[1].t[imic].Fchap[list2]=Qchap2;
}


struct apply_type_elem{
  template<class TE,class TV> void operator() ( TE &e, unsigned &compt, TV &status, unsigned &dim) const{
    e.type=status[compt*dim];
    compt+=1;
  }
};


//interface interieure de type contact avec frottement
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procédure pour les interfaces intérieures de type contact avec frottement, avec jeu ou non
*/
template<class INTER>
void compt_contact (INTER &Inter,TEMPS &temps) {
    typedef typename INTER::T TT;

    int imic = temps.pt;
    double dt=temps.dt;

    Vec<unsigned> &list1=Inter.side[0].ddlcorresp;
    Vec<unsigned> &list2=Inter.side[1].ddlcorresp;
    Vec<TT> WWpchap1=Inter.side[0].t[imic].Wpchap[list1];
    Vec<TT> WWpchap2=Inter.side[1].t[imic].Wpchap[list2];
    Vec<TT> Qchap1=Inter.side[0].t[imic].Fchap[list1];
    Vec<TT> Qchap2=Inter.side[1].t[imic].Fchap[list2];
    const Vec<TT> &Q1=Inter.side[0].t[imic].F[list1];
    const Vec<TT> &Q2=Inter.side[1].t[imic].F[list2];
    const Vec<TT> &WWp1=Inter.side[0].t[imic].Wp[list1];
    const Vec<TT> &WWp2=Inter.side[1].t[imic].Wp[list2];
    const Vec<TT> &WWchap1old=Inter.side[0].t[imic-1].Wchap[list1];
    const Vec<TT> &WWchap2old=Inter.side[1].t[imic-1].Wchap[list2];
    Vec<TT> WWchap1=Inter.side[0].t[imic].Wchap[list1];
    Vec<TT> WWchap2=Inter.side[1].t[imic].Wchap[list2];
    
    const Vec<TT> &neq=Inter.side[0].neq[list1];

//     if (Inter.num == 15) cout << list1 << endl;
//     if (Inter.num == 15) cout << list2 << endl;
//     if (Inter.num == 15) cout << Inter.side[0].neq << endl;
    
    TT f=Inter.param_comp->coeffrottement;

    TT h1n=Inter.side[0].hn;
    TT h2n=Inter.side[1].hn;
    TT h1t=Inter.side[0].ht;
    TT h2t=Inter.side[1].ht;
    unsigned dim=INTER::dim;

    // travail pt par pt
    unsigned nbpts=Inter.side[0].nodeeq.size();
// #ifdef LOOK_CONTACT_ZONE
//     Vec<unsigned> status;status.resize(nbpts*INTER::dim);
// #endif

    for(unsigned i=0;i<nbpts;++i) {
        Vec<unsigned> rep=range(i*dim,(i+1)*dim);
        Vec<TT,INTER::dim> n=neq[rep]; // normale de 1 vers 2
        Vec<TT,INTER::dim> F1=Q1[rep], F2=Q2[rep],Fchap1=Qchap1[rep], Fchap2=Qchap2[rep], Wp1=WWp1[rep],Wp2=WWp2[rep],Wpchap1=WWpchap1[rep],Wpchap2=WWpchap2[rep] , Wchap1old=WWchap1old[rep],Wchap2old=WWchap2old[rep] , Wchap1=WWchap1[rep],Wchap2=WWchap2[rep];

        // test contact normal
        TT N=0;
        //if (temps.type_de_calcul=="stat")
        //    N = ( dot(n,(Wp2-Wp1)) -h2n*dot(n,F2) + h1n*dot(n,F1) )/(h2n+h1n);
        //else
            N=( dot(n,(Wp2-Wp1)) + dot(n,(Wchap2old-Wchap1old)/dt) -h2n*dot(n,F2) + h1n*dot(n,F1) )/(h2n+h1n);

        if (N>0.0) { // separation de l'interface
            Fchap1=0.0;
            Fchap2=0.0;
            Wpchap1=Wp1-h1n*dot(n,F1)*n - h1t* ProjT(F1,n);
            Wpchap2=Wp2-h2n*dot(n,F2)*n - h2t* ProjT(F2,n);
/*#ifdef LOOK_CONTACT_ZONE
            status[i*INTER::dim]=0;
#endif*/
        } else { // contact
            TT Fchap1n=N;
            TT Fchap2n=-1.0*Fchap1n;
            TT Wpchap1n=dot(Wp1,n) + h1n*( Fchap1n - dot(n,F1) );
            TT Wpchap2n=dot(Wp2,n) + h2n*( Fchap2n - dot(n,F2) );

            // test glissement adherence
            Vec<TT,INTER::dim> T;
            T=(ProjT(Wp2,n)-ProjT(Wp1,n) -h2t*ProjT(F2,n) + h1t*ProjT(F1,n))/(h2t+h1t);
            TT g = f*std::abs(Fchap1n);

            TT normT=norm_2(T);
            Vec<TT,INTER::dim> Fchap1t,Fchap2t,Wpchap1t,Wpchap2t;

            //Attention au signe : modifie pour eviter la division par 0
            if (normT <= g) { // adherence
                Fchap1t=T;
                Fchap2t=-1.0*Fchap1t;
                Wpchap1t= ProjT(Wp1,n) + h1t*(Fchap1t - ProjT(F1,n) );
                Wpchap2t= Wpchap1t;
/*#ifdef LOOK_CONTACT_ZONE
                status[i*INTER::dim]=1;
#endif*/
            } else if (normT > g) { // glissement
                Fchap1t=T*g/normT;
                Fchap2t=-1.0*Fchap1t;
                Wpchap1t= ProjT(Wp1,n) + h1t*(Fchap1t - ProjT(F1,n) );
                Wpchap2t= ProjT(Wp2,n) + h2t*(Fchap2t - ProjT(F2,n) );
/*#ifdef LOOK_CONTACT_ZONE
                status[i*INTER::dim]=2;
#endif*/
            }

            Wpchap1=Wpchap1n*n+Wpchap1t;
            Fchap1=Fchap1n*n+Fchap1t;
            Wpchap2=Wpchap2n*n+Wpchap2t;
            Fchap2=Fchap2n*n+Fchap2t;
        }

        //integration
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
    
/*#ifdef LOOK_CONTACT_ZONE
    Vec<unsigned> status1,status2;status1.resize(nbpts*INTER::dim);status2.resize(nbpts*INTER::dim);
    status1[list1]=status;
    status2[list2]=status;
    unsigned compt=0;
    apply(Inter.side[0].mesh->elem_list,apply_type_elem(),compt,status1,dim);
    compt=0;
    apply(Inter.side[1].mesh->elem_list,apply_type_elem(),compt,status2,dim);
#endif*/
}



//interface interieure de type defaut de forme
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procédure pour les interfaces intérieures de type defaut de forme
*/
template<class INTER>
void compt_contact_ep (INTER &Inter,TEMPS &temps) {
    typedef typename INTER::T TT;

    int imic = temps.pt;
    double dt=temps.dt;

    Vec<unsigned> &list1=Inter.side[0].ddlcorresp;
    Vec<unsigned> &list2=Inter.side[1].ddlcorresp;
    Vec<TT> WWpchap1=Inter.side[0].t[imic].Wpchap[list1];
    Vec<TT> WWpchap2=Inter.side[1].t[imic].Wpchap[list2];
    Vec<TT> Qchap1=Inter.side[0].t[imic].Fchap[list1];
    Vec<TT> Qchap2=Inter.side[1].t[imic].Fchap[list2];
    const Vec<TT> &Q1=Inter.side[0].t[imic].F[list1];
    const Vec<TT> &Q2=Inter.side[1].t[imic].F[list2];
    const Vec<TT> &WWp1=Inter.side[0].t[imic].Wp[list1];
    const Vec<TT> &WWp2=Inter.side[1].t[imic].Wp[list2];
    const Vec<TT> &WWchap1old=Inter.side[0].t[imic-1].Wchap[list1];
    const Vec<TT> &WWchap2old=Inter.side[1].t[imic-1].Wchap[list2];
    Vec<TT> WWchap1=Inter.side[0].t[imic].Wchap[list1];
    Vec<TT> WWchap2=Inter.side[1].t[imic].Wchap[list2];
    
    const Vec<TT> &JJ=Inter.param_comp->jeu[list1];
    const Vec<TT> &neq=Inter.side[1].neq[list1];
//    const Vec<TT> &neq=Inter.side[0].neq[list1];

    TT f=Inter.param_comp->coeffrottement;
    const Vec<TT> &frott=Inter.param_comp->f_coeffrottement;

    TT h1n=Inter.side[0].hn;
    TT h2n=Inter.side[1].hn;
    TT h1t=Inter.side[0].ht;
    TT h2t=Inter.side[1].ht;
    unsigned dim=INTER::dim;

    // travail pt par pt
    unsigned nbpts=Inter.side[0].nodeeq.size();
// #ifdef LOOK_CONTACT_ZONE
//     Vec<unsigned> status;status.resize(nbpts*INTER::dim);
// #endif

    for(unsigned i=0;i<nbpts;++i) {
        Vec<unsigned> rep=range(i*dim,(i+1)*dim);
//        Vec<TT,INTER::dim> n=neq[rep]; // normale de 1 vers 2
        Vec<TT,INTER::dim> n=-1.*neq[rep]; // normale de 2 vers 1
        Vec<TT,INTER::dim> F1=Q1[rep], F2=Q2[rep],Fchap1=Qchap1[rep], Fchap2=Qchap2[rep], Wp1=WWp1[rep],Wp2=WWp2[rep],Wpchap1=WWpchap1[rep],Wpchap2=WWpchap2[rep] , Wchap1old=WWchap1old[rep],Wchap2old=WWchap2old[rep] , Wchap1=WWchap1[rep],Wchap2=WWchap2[rep],jeu=JJ[rep];

        // test contact normal
        TT N=0;
        //if (temps.type_de_calcul=="stat")
        //    N = ( dot(n,(Wp2-Wp1)) -h2n*dot(n,F2) + h1n*dot(n,F1) )/(h2n+h1n);
        //else
        Wp1 += jeu;
        Wchap1old += jeu;
        N=( dot(n,(Wp2-(Wp1))) + dot(n,(Wchap2old-(Wchap1old))/dt) -h2n*dot(n,F2) + h1n*dot(n,F1) )/(h2n+h1n);

        if (N>0.0) { // separation de l'interface
            Fchap1=0.0;
            Fchap2=0.0;
            Wpchap1=Wp1-h1n*dot(n,F1)*n - h1t* ProjT(F1,n);
            Wpchap2=Wp2-h2n*dot(n,F2)*n - h2t* ProjT(F2,n);
/*#ifdef LOOK_CONTACT_ZONE
            status[i*INTER::dim]=0;
#endif*/
        } else { // contact
            TT Fchap1n=N;
            TT Fchap2n=-1.0*Fchap1n;
            TT Wpchap1n=dot(Wp1,n) + h1n*( Fchap1n - dot(n,F1) );
            TT Wpchap2n=dot(Wp2,n) + h2n*( Fchap2n - dot(n,F2) );

            // test glissement adherence
            Vec<TT,INTER::dim> Test;
            Test=(ProjT(Wp2,n)-ProjT(Wp1,n) -h2t*ProjT(F2,n) + h1t*ProjT(F1,n))/(h2t+h1t);
            TT g = frott[i]*std::abs(Fchap1n);
//             TT g = f*std::abs(Fchap1n);

            TT normT=norm_2(Test);
            Vec<TT,INTER::dim> Fchap1t,Fchap2t,Wpchap1t,Wpchap2t;

            //Attention au signe : modifie pour eviter la division par 0
            if (normT <= g) { // adherence
                Fchap1t=Test;
                Fchap2t=-1.0*Fchap1t;
                Wpchap1t= ProjT(Wp1,n) + h1t*(Fchap1t - ProjT(F1,n) );
                Wpchap2t= Wpchap1t;
/*#ifdef LOOK_CONTACT_ZONE
                status[i*INTER::dim]=1;
#endif*/
            } else if (normT > g) { // glissement
                Fchap1t=Test*g/normT;
                Fchap2t=-1.0*Fchap1t;
                Wpchap1t= ProjT(Wp1,n) + h1t*(Fchap1t - ProjT(F1,n) );
                Wpchap2t= ProjT(Wp2,n) + h2t*(Fchap2t - ProjT(F2,n) );
/*#ifdef LOOK_CONTACT_ZONE
                status[i*INTER::dim]=2;
#endif*/
            }

            Wpchap1=Wpchap1n*n + Wpchap1t - jeu;
            Fchap1=Fchap1n*n+Fchap1t;
            Wpchap2=Wpchap2n*n+Wpchap2t;
            Fchap2=Fchap2n*n+Fchap2t;
        }

        //integration
        Wchap1 = Wpchap1*dt + (Wchap1old-jeu);
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
    Inter.side[0].t[imic].Wchap[list1]=WWchap1;
    Inter.side[1].t[imic].Wchap[list2]=WWchap2;
    Inter.side[0].t[imic].Fchap[list1]=Qchap1;
    Inter.side[1].t[imic].Fchap[list2]=Qchap2;
}





struct assign_d_mesh {
   template<class TE,class T> void operator() (TE &e,unsigned &incr, Vec<T> &d) const {
      e.d=d[incr];
      incr+=1;
   }
};


//interface interieure avec loi cohesive (type mesomodele)
/** \ingroup  etape_locale
\relates etape_locale_inter
\brief Procédure pour les interfaces intérieures de type cohésif
*/
template<class INTER>
void compt_cohesif (INTER &Inter,TEMPS &temps) {

    int imic = temps.pt;
    double dt=temps.dt;

    typedef typename INTER::T TT;
    typedef  Vec<TT,INTER::dim> TV;
    typedef Mat<TT,Gen<INTER::dim>,Dense<> > TMAT;

    Vec<unsigned> &list1=Inter.side[0].ddlcorresp;
    Vec<unsigned> &list2=Inter.side[1].ddlcorresp;
    Vec<TT> WWchap1=Inter.side[0].t[imic].Wchap[list1];
    Vec<TT> WWchap2=Inter.side[1].t[imic].Wchap[list2];
    Vec<TT> FFchap1=Inter.side[0].t[imic].Fchap[list1];
    Vec<TT> FFchap2=Inter.side[1].t[imic].Fchap[list2];
    const Vec<TT> &FF1=Inter.side[0].t[imic].F[list1];
    const Vec<TT> &FF2=Inter.side[1].t[imic].F[list2];
    const Vec<TT> &WWp1=Inter.side[0].t[imic].Wp[list1];
    const Vec<TT> &WWp2=Inter.side[1].t[imic].Wp[list2];
    const Vec<TT> &WWchap1old=Inter.side[0].t[imic-1].Wchap[list1];
    const Vec<TT> &WWchap2old=Inter.side[1].t[imic-1].Wchap[list2];

    const Vec<TT> &neq=Inter.side[0].neq[list1];

    TMAT Id;
    Id.set(0.);
    Id.diag()=Vec<TT,INTER::dim>(1);

    double eps=1e-3;

    typename INTER::PARAM_COMP::PARAM_DAMAGE damage=Inter.param_comp->param_damage;

    // parametres de direction de recherche
    TT h1n=Inter.side[0].hn;
    TT h2n=Inter.side[1].hn;
    TT h1t=Inter.side[0].ht;
    TT h2t=Inter.side[1].ht;

    // travail pt par pt
    unsigned nbpts=Inter.side[0].nodeeq.size();
    for(unsigned i=0;i<nbpts;++i) {
        Vec<unsigned> rep=range(i*INTER::dim,(i+1)*INTER::dim);
        TV n=neq[rep]; // normale de 1 vers 2
        //matrice normale nxn
        TMAT ntn;
        tens(n,n,ntn);
        //matrices des directions de recherche
        TMAT H1 = h1n*ntn + h1t*(Id-ntn);
        TMAT H2 = h2n*ntn + h2t*(Id-ntn);
        
        //variables connues
        TV F1=FF1[rep], F2=FF2[rep], Wp1=WWp1[rep],Wp2=WWp2[rep],Wchap1old=WWchap1old[rep],Wchap2old=WWchap2old[rep];
        //variables inconnues
        TV Fchap1(0.), Fchap2(0.),Wchap1=WWchap1[rep],Wchap2=WWchap2[rep],Fchap1t(0.);

        //initialisation de la boucle iterative
        TT dold=Inter.param_comp->dmax[i];
        TT d = dold;
        TT Ymax = Inter.param_comp->Ymax[i];
        TT erreur = 1.;
        unsigned it=0;
        //boucle iterative
        while(erreur>eps) {
            it+=1;

            if(dold>=1.-eps) {
                //etape necessaire lorsque d vaut 1 dans la boucle itérative : on ne gere pas le frottement tant que toute l'interface n'est pas cassee
                TT N=0;
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
                //  Wchap1=(Wp1 + H1*(Fchap1 - F1) );
                //  Wchap2=(Wp2 + H2*(Fchap2 - F2) );                
                //}
                //else{cout << "Mauvais type de calcul" << endl; assert(0);}
                
                erreur=0;//pour sortir de la boucle
            } else {//dans tous les autres cas
                //matrice pour relation de comportement d'interface prenant en compte l'endommagement
                TMAT Kc;
                Kc = damage.knc*ntn +  damage.kt*(1.-dold)*(Id - ntn);
                TMAT Kt;
                Kt = damage.kn*(1.-dold)*ntn +  damage.kt*(1.-dold)*(Id - ntn);
                TMAT invKc = 1./damage.knc*ntn + 1./(damage.kt*(1.-dold))*(Id - ntn);
                TMAT invKt = 1./(damage.kn*(1.-dold))*ntn +  1./(damage.kt*(1.-dold))*(Id - ntn);

                //1ere etape : calcul de Fchap1
                TV Wtemp;
                //if(temps.type_de_calcul=="Qstat")
                  Wtemp= ( -1.*H2*F2 + H1*F1 + (Wp2-Wp1) + (Wchap2old-Wchap1old)/dt );
                //else if(temps.type_de_calcul=="stat")
                //  Wtemp= ( -1.*H2*F2 + H1*F1 + (Wp2-Wp1)  );
                //else{cout << "Mauvais type de calcul" << endl; assert(0);}
                
                Fchap1 = inv(H2+H1+invKc)*Wtemp;
                if (dot(Fchap1,n)>=0) {
                    Fchap1 = inv(H2+H1+invKt)*Wtemp;
                }
                //calcul de Fchap2, Wchap1, Wchap2
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
                
                //calcul des forces d'endommagement normales et tangentielles a partir de Fchap1
                TT Ydn = (dot(Fchap1,n)>=0)*( std::pow((dot(Fchap1,n)/(1.-d)),2.)/(2.*damage.kn) );
                Fchap1t = Fchap1 - dot(Fchap1,n)*n;
                TT Yt = ( std::pow(norm_2(Fchap1t)/(1.-d),2.)/(2.*damage.kt) );
                TT Ynew = std::pow( std::pow(Ydn,damage.alpha)+ (damage.gamma*std::pow(Yt,damage.alpha)) , 1./damage.alpha);
                Ymax = std::max(Ynew,Inter.param_comp->Ymax[i]);
                //calcul du nouvel endommagement pour cette iteration
                TT Ypos = ((Ymax-damage.Yo)>=0)*(Ymax-damage.Yo);
                //        cout<< "ftemp " << Ftemp << " nodeeq "<<Inter.side[0].nodeeq[i]<< endl;
                //         cout<< " n " <<n <<"F1 " << F1 << " F2 " << F2 << " nodeeq "<<Inter.side[0].nodeeq[i]<< endl;
                //      cout << "Appuyez sur entrÃ©e pour continuer...";
                //        cin.ignore( numeric_limits<streamsize>::max(), '\n' );

                d = std::pow((damage.n/(damage.n+1))*(Ypos/(damage.Yc-damage.Yo)),damage.n);
                d = std::min(d,1.);
                if (d>=1.-eps) {
                    //test contact (sans frottement)
                    TT N=0;
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
        //modification dmax pour le pas de temps suivant si necessaire
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
/*      unsigned incr=0;
      apply(Inter.side[0].mesh.elem_list,assign_d_mesh(),incr,Inter.param_comp->t[imic].d);
      incr=0;
      apply(Inter.side[1].mesh.elem_list,assign_d_mesh(),incr,Inter.param_comp->t[imic].d);*/
   
      //modification du type d'interface si dmax = 1 pour tous les points
/*      TT dmin;
      for(unsigned i=0;i<nbpts;i++) {
         dmin = min(dmin,Inter.param_comp->d[i]);
      }
      if (dmin>=1.-eps) {
         cout <<Inter.num <<  " Inter modif" << endl;
         Inter.comp="Contact";
      }  */    

}


