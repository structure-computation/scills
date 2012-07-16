#ifndef MPISENDRECVALL
#define MPISENDRECVALL
//
// C++ Interface: mpi_sendrecvall
//
// Description:
//
//
// Author: Alain CAIGNOT <caignot@lmt.ens-cachan.fr>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "mpi_lmt_functions.h"

using namespace LMT;


template <class T1, class T2>
void SendInterAll(T1 &vecintertoexchange,T2 &Inter,MPI_Request &request,Process &process, Vec<TYPEREEL> &vectosend) {
    int nbqsend=6;
    vectosend.resize(0);
    for( unsigned i=0;i<vecintertoexchange.size() ;i++ ) {
        int sizebefore=vectosend.size();
        vectosend.resize(sizebefore+nbqsend*Inter[vecintertoexchange[i].num].side[0].t.size()*Inter[vecintertoexchange[i].num].side[0].t[0].Fchap.size()+2);
        vectosend[sizebefore]=vecintertoexchange[i].num;
        vectosend[sizebefore+1]=vecintertoexchange[i].sidetosend;
        for(unsigned pt=0;pt<Inter[vecintertoexchange[i].num].side[0].t.size();pt++) {
            vectosend[range(sizebefore+2+nbqsend*pt*Inter[vecintertoexchange[i].num].side[0].t[0].Fchap.size(),sizebefore+2+(nbqsend*pt+1)*Inter[vecintertoexchange[i].num].side[0].t[0].Fchap.size())]=Inter[vecintertoexchange[i].num].side[vecintertoexchange[i].sidetosend].t[pt].W;
            vectosend[range(sizebefore+2+(nbqsend*pt+1)*Inter[vecintertoexchange[i].num].side[0].t[0].Fchap.size(),sizebefore+2+(nbqsend*pt+2)*Inter[vecintertoexchange[i].num].side[0].t[0].Fchap.size())]=Inter[vecintertoexchange[i].num].side[vecintertoexchange[i].sidetosend].t[pt].Wchap;
            vectosend[range(sizebefore+2+(nbqsend*pt+2)*Inter[vecintertoexchange[i].num].side[0].t[0].Fchap.size(),sizebefore+2+(nbqsend*pt+3)*Inter[vecintertoexchange[i].num].side[0].t[0].Fchap.size())]=Inter[vecintertoexchange[i].num].side[vecintertoexchange[i].sidetosend].t[pt].F;
            vectosend[range(sizebefore+2+(nbqsend*pt+3)*Inter[vecintertoexchange[i].num].side[0].t[0].Fchap.size(),sizebefore+2+(nbqsend*pt+4)*Inter[vecintertoexchange[i].num].side[0].t[0].Fchap.size())]=Inter[vecintertoexchange[i].num].side[vecintertoexchange[i].sidetosend].t[pt].Fchap;
            vectosend[range(sizebefore+2+(nbqsend*pt+4)*Inter[vecintertoexchange[i].num].side[0].t[0].Fchap.size(),sizebefore+2+(nbqsend*pt+5)*Inter[vecintertoexchange[i].num].side[0].t[0].Fchap.size())]=Inter[vecintertoexchange[i].num].side[vecintertoexchange[i].sidetosend].t[pt].Wp;
            vectosend[range(sizebefore+2+(nbqsend*pt+5)*Inter[vecintertoexchange[i].num].side[0].t[0].Fchap.size(),sizebefore+2+(nbqsend*pt+6)*Inter[vecintertoexchange[i].num].side[0].t[0].Fchap.size())]=Inter[vecintertoexchange[i].num].side[vecintertoexchange[i].sidetosend].t[pt].Wpchap;
        }
    }
    MPI_Isend(vectosend,vecintertoexchange[0].to,request);
}

template <class T1, class T2>
void RecvInterAll(T1 &MPIsource,T2 &Inter, Process &process) {
    int nbqsend=6;
    Vec<TYPEREEL> vectorecv;
    MPI_Recv(vectorecv,MPIsource);
    int sizebefore=0,num,side;

    while( true ) {
        num=(int) vectorecv[sizebefore];
        side=(int) vectorecv[sizebefore+1];
        for(unsigned pt=0;pt<Inter[num].side[0].t.size();pt++) {
            Inter[num].side[side].t[pt].W=vectorecv[range(sizebefore+2+nbqsend*pt*Inter[num].side[0].t[0].Fchap.size(),sizebefore+2+(nbqsend*pt+1)*Inter[num].side[0].t[0].Fchap.size())];
            Inter[num].side[side].t[pt].Wchap=vectorecv[range(sizebefore+2+(nbqsend*pt+1)*Inter[num].side[0].t[0].Fchap.size(),sizebefore+2+(nbqsend*pt+2)*Inter[num].side[0].t[0].Fchap.size())];
            Inter[num].side[side].t[pt].F=vectorecv[range(sizebefore+2+(nbqsend*pt+2)*Inter[num].side[0].t[0].Fchap.size(),sizebefore+2+(nbqsend*pt+3)*Inter[num].side[0].t[0].Fchap.size())];
            Inter[num].side[side].t[pt].Fchap=vectorecv[range(sizebefore+2+(nbqsend*pt+3)*Inter[num].side[0].t[0].Fchap.size(),sizebefore+2+(nbqsend*pt+4)*Inter[num].side[0].t[0].Fchap.size())];
            Inter[num].side[side].t[pt].Wp=vectorecv[range(sizebefore+2+(nbqsend*pt+4)*Inter[num].side[0].t[0].Fchap.size(),sizebefore+2+(nbqsend*pt+5)*Inter[num].side[0].t[0].Fchap.size())];
            Inter[num].side[side].t[pt].Wpchap=vectorecv[range(sizebefore+2+(nbqsend*pt+5)*Inter[num].side[0].t[0].Fchap.size(),sizebefore+2+(nbqsend*pt+6)*Inter[num].side[0].t[0].Fchap.size())];
        }
        sizebefore+=nbqsend*Inter[num].side[0].t.size()*Inter[num].side[0].t[0].Fchap.size()+2;
        if (sizebefore >= (int)vectorecv.size())
            break;
    }
}

template <class T1, class T2>
void SendInterAllPost(T1 &vecintertoexchange,T2 &Inter,MPI_Request &request,Process &process, Vec<TYPEREEL> &vectosend) {
    int nbqsend=6;
    vectosend.resize(0);
    for( unsigned i=0;i<vecintertoexchange.size() ;i++ ) {
        int sizebefore=vectosend.size();
        vectosend.resize(sizebefore+nbqsend*Inter[vecintertoexchange[i].num].side[0].t_post.size()*Inter[vecintertoexchange[i].num].side[0].t_post[0].Fchap.size()+2);
        vectosend[sizebefore]=vecintertoexchange[i].num;
        vectosend[sizebefore+1]=vecintertoexchange[i].sidetosend;
        for(unsigned pt=0;pt<Inter[vecintertoexchange[i].num].side[0].t_post.size();pt++) {
            vectosend[range(sizebefore+2+nbqsend*pt*Inter[vecintertoexchange[i].num].side[0].t_post[0].Fchap.size(),sizebefore+2+(nbqsend*pt+1)*Inter[vecintertoexchange[i].num].side[0].t_post[0].Fchap.size())]=Inter[vecintertoexchange[i].num].side[vecintertoexchange[i].sidetosend].t_post[pt].W;
            vectosend[range(sizebefore+2+(nbqsend*pt+1)*Inter[vecintertoexchange[i].num].side[0].t_post[0].Fchap.size(),sizebefore+2+(nbqsend*pt+2)*Inter[vecintertoexchange[i].num].side[0].t_post[0].Fchap.size())]=Inter[vecintertoexchange[i].num].side[vecintertoexchange[i].sidetosend].t_post[pt].Wchap;
            vectosend[range(sizebefore+2+(nbqsend*pt+2)*Inter[vecintertoexchange[i].num].side[0].t_post[0].Fchap.size(),sizebefore+2+(nbqsend*pt+3)*Inter[vecintertoexchange[i].num].side[0].t_post[0].Fchap.size())]=Inter[vecintertoexchange[i].num].side[vecintertoexchange[i].sidetosend].t_post[pt].F;
            vectosend[range(sizebefore+2+(nbqsend*pt+3)*Inter[vecintertoexchange[i].num].side[0].t_post[0].Fchap.size(),sizebefore+2+(nbqsend*pt+4)*Inter[vecintertoexchange[i].num].side[0].t_post[0].Fchap.size())]=Inter[vecintertoexchange[i].num].side[vecintertoexchange[i].sidetosend].t_post[pt].Fchap;
            vectosend[range(sizebefore+2+(nbqsend*pt+4)*Inter[vecintertoexchange[i].num].side[0].t_post[0].Fchap.size(),sizebefore+2+(nbqsend*pt+5)*Inter[vecintertoexchange[i].num].side[0].t_post[0].Fchap.size())]=Inter[vecintertoexchange[i].num].side[vecintertoexchange[i].sidetosend].t_post[pt].Wp;
            vectosend[range(sizebefore+2+(nbqsend*pt+5)*Inter[vecintertoexchange[i].num].side[0].t_post[0].Fchap.size(),sizebefore+2+(nbqsend*pt+6)*Inter[vecintertoexchange[i].num].side[0].t_post[0].Fchap.size())]=Inter[vecintertoexchange[i].num].side[vecintertoexchange[i].sidetosend].t_post[pt].Wpchap;
        }
    }
    MPI_Isend(vectosend,vecintertoexchange[0].to,request);
}

template <class T1, class T2>
void RecvInterAllPost(T1 &MPIsource,T2 &Inter, Process &process) {
    int nbqsend=6;
    Vec<TYPEREEL> vectorecv;
    MPI_Recv(vectorecv,MPIsource);
    int sizebefore=0,num,side;

    while( true ) {
        num=(int) vectorecv[sizebefore];
        side=(int) vectorecv[sizebefore+1];
        for(unsigned pt=0;pt<Inter[num].side[0].t_post.size();pt++) {
            Inter[num].side[side].t_post[pt].W=vectorecv[range(sizebefore+2+nbqsend*pt*Inter[num].side[0].t_post[0].Fchap.size(),sizebefore+2+(nbqsend*pt+1)*Inter[num].side[0].t_post[0].Fchap.size())];
            Inter[num].side[side].t_post[pt].Wchap=vectorecv[range(sizebefore+2+(nbqsend*pt+1)*Inter[num].side[0].t_post[0].Fchap.size(),sizebefore+2+(nbqsend*pt+2)*Inter[num].side[0].t_post[0].Fchap.size())];
            Inter[num].side[side].t_post[pt].F=vectorecv[range(sizebefore+2+(nbqsend*pt+2)*Inter[num].side[0].t_post[0].Fchap.size(),sizebefore+2+(nbqsend*pt+3)*Inter[num].side[0].t_post[0].Fchap.size())];
            Inter[num].side[side].t_post[pt].Fchap=vectorecv[range(sizebefore+2+(nbqsend*pt+3)*Inter[num].side[0].t_post[0].Fchap.size(),sizebefore+2+(nbqsend*pt+4)*Inter[num].side[0].t_post[0].Fchap.size())];
            Inter[num].side[side].t_post[pt].Wp=vectorecv[range(sizebefore+2+(nbqsend*pt+4)*Inter[num].side[0].t_post[0].Fchap.size(),sizebefore+2+(nbqsend*pt+5)*Inter[num].side[0].t_post[0].Fchap.size())];
            Inter[num].side[side].t_post[pt].Wpchap=vectorecv[range(sizebefore+2+(nbqsend*pt+5)*Inter[num].side[0].t_post[0].Fchap.size(),sizebefore+2+(nbqsend*pt+6)*Inter[num].side[0].t_post[0].Fchap.size())];
        }
        sizebefore+=nbqsend*Inter[num].side[0].t_post.size()*Inter[num].side[0].t_post[0].Fchap.size()+2;
        if (sizebefore >= (int)vectorecv.size())
            break;
    }
}


template <class T1,class T2>
void SendRecvInterAll(T1 &vecintertoexchangebypro, T2 &Inter, Process &process) {
    Vec<Vec<TYPEREEL> > vectosend;
    vectosend.resize(vecintertoexchangebypro.size());

    MPI_Request request[vecintertoexchangebypro.size()];
    for( unsigned i=0;i< vecintertoexchangebypro.size() ;i++ ) {
      if (process.nom_calcul=="latin")
        SendInterAll(vecintertoexchangebypro[i].inter,Inter,request[i],process,vectosend[i]);
      if (process.nom_calcul=="incr")
        SendInterAllPost(vecintertoexchangebypro[i].inter,Inter,request[i],process,vectosend[i]);
    }
    int nbtorecv = vecintertoexchangebypro.size();
    while( nbtorecv != 0) {
        //On regarde si on a reçu un message de quelqu'un (MPI_ANY_SOURCE).
        int flag;
        MPI_Status status1;
        MPI_Iprobe(MPI_ANY_SOURCE,201,MPI_COMM_WORLD,&flag, &status1);
        //Si on a recu un message alors on va le traiter
        if (flag == 1) {
          if (process.nom_calcul=="latin")
            RecvInterAll(status1.MPI_SOURCE,Inter,process);
          if (process.nom_calcul=="incr")
            RecvInterAllPost(status1.MPI_SOURCE,Inter,process);
          nbtorecv--;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for( unsigned i=0;i< vecintertoexchangebypro.size() ;i++ )
        MPI_Request_free(&request[i]);
}


#endif //MPISENDRECVALL
