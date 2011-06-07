#ifndef MPI_TRANSACTION
#define MPI_TRANSACTION
//
// C++ Interface: mpi_transactions
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

using namespace std;
using namespace LMT;

template <class TV>
void SendbigF(Param &process,TV &bigF ) {
    Vec<TYPEREEL> temp;
    temp.resize(bigF.size());
    if (process.rank==0)
      bigF.set(0.0);
    MPI_Reduce(bigF.ptr(),temp.ptr(),bigF.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if (process.rank==0)
      bigF=temp;

}


template <class T1, class T2>
void RecvInter(T1 &MPIsource,T2 &Inter, Param &process) {
    Vec<TYPEREEL> vectorecv;
    MPI_Recv(vectorecv,MPIsource);
    int repere=0,num,side;
    while( true ) {
        num=(int) vectorecv[repere];
        side=(int) vectorecv[repere+1];
        for(unsigned pt=0;pt<Inter[num].side[0].t.size();pt++) {
            Inter[num].side[side].t[pt].F=vectorecv[range(repere+2+2*pt*Inter[num].side[0].t[0].F.size(),repere+2+(2*pt+1)*Inter[num].side[0].t[0].F.size())];
            Inter[num].side[side].t[pt].Wp=vectorecv[range(repere+2+(2*pt+1)*Inter[num].side[0].t[0].F.size(),repere+2+(2*pt+2)*Inter[num].side[0].t[0].F.size())];
        }
        repere+=2*Inter[num].side[0].t.size()*Inter[num].side[0].t[0].F.size()+2;
        if (repere >= (int)vectorecv.size())
            break;
    }

}

template <class T1, class T2>
void SendInter(T1 &vecintertoexchange,T2 &Inter,MPI_Request &request,Param &process, Vec<TYPEREEL> &vectosend) {
    vectosend.resize(0);
    for( unsigned i=0;i<vecintertoexchange.size() ;i++ ) {
        int sizebefore=vectosend.size();
        vectosend.resize(sizebefore+2*Inter[vecintertoexchange[i].num].side[0].t.size()*Inter[vecintertoexchange[i].num].side[0].t[0].F.size()+2);
        vectosend[sizebefore]=vecintertoexchange[i].num;
        vectosend[sizebefore+1]=vecintertoexchange[i].sidetosend;
        for(unsigned pt=0;pt<Inter[vecintertoexchange[i].num].side[0].t.size();pt++) {

            vectosend[range(sizebefore+2+2*pt*Inter[vecintertoexchange[i].num].side[0].t[0].F.size(),sizebefore+2+(2*pt+1)*Inter[vecintertoexchange[i].num].side[0].t[0].F.size())]=Inter[vecintertoexchange[i].num].side[vecintertoexchange[i].sidetosend].t[pt].F;
            vectosend[range(sizebefore+2+(2*pt+1)*Inter[vecintertoexchange[i].num].side[0].t[0].F.size(),sizebefore+2+(2*pt+2)*Inter[vecintertoexchange[i].num].side[0].t[0].F.size())]=Inter[vecintertoexchange[i].num].side[vecintertoexchange[i].sidetosend].t[pt].Wp;
        }
    }
    MPI_Isend(vectosend,vecintertoexchange[0].to,request);
}

template <class T1,class T2>
void SendRecvInter(T1 &vecintertoexchangebypro, T2 &Inter, Param &process) {
    Vec<Vec<TYPEREEL> > vectosend;
    vectosend.resize(vecintertoexchangebypro.size());

    MPI_Request request[vecintertoexchangebypro.size()];
    for( unsigned i=0;i< vecintertoexchangebypro.size() ;i++ ) {
        SendInter(vecintertoexchangebypro[i].inter,Inter,request[i],process,vectosend[i]);
    }
    int nbtorecv = vecintertoexchangebypro.size();
    while( nbtorecv != 0) {
        //On regarde si on a reçu un message de quelqu'un (MPI_ANY_SOURCE).
        int flag;
        MPI_Status status1;
        MPI_Iprobe(MPI_ANY_SOURCE,201,MPI_COMM_WORLD,&flag, &status1);
        //Si on a recu un message alors on va le traiter
        if (flag == 1) {
            RecvInter(status1.MPI_SOURCE,Inter,process);
            nbtorecv--;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for( unsigned i=0;i< vecintertoexchangebypro.size() ;i++ )
        MPI_Request_free(&request[i]);
}


#endif //MPI_TRANSACTION
