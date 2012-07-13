#ifndef SENDRECVTOMASTER_H
#define SENDRECVTOMASTER_H
//
// C++ Interface: mpi_sendrecvtomaster
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

template <class T1>
void SendSST(T1 &SubS,MPI_Request &request,Process &process, Vec<TYPEREEL> &vectosend) {
    vectosend.resize(0);
    for( unsigned i=0;i<SubS.size() ;i++ ) {
      int sizebefore=vectosend.size();
        if (process.nom_calcul=="latin") {
            vectosend.resize(sizebefore+SubS[i].t[1].q.size()*SubS[i].t.size());
            for(unsigned pt=1;pt<SubS[i].t.size();pt++) {
                vectosend[range(sizebefore+pt*SubS[i].t[pt].q.size(),sizebefore+(pt+1)*SubS[i].t[pt].q.size())]=SubS[i].t[pt].q;
            }
        } else {
            vectosend.resize(sizebefore+SubS[i].t_post[1].q.size()*SubS[i].t_post.size());
            for(unsigned pt=1;pt<SubS[i].t_post.size();pt++) {
              vectosend[range(sizebefore+pt*SubS[i].t_post[pt].q.size(),sizebefore+(pt+1)*SubS[i].t_post[pt].q.size())]=SubS[i].t_post[pt].q;
            }
        }
    }
    MPI_Isend(vectosend,0,request,204);
}
template <class T1>
void RecvSST(T1 &S,int &MPIsource,Process &process) {
    Vec<TYPEREEL> vectorecv;
    MPI_Recv(vectorecv,MPIsource,204);
    int repere=0;
    for(unsigned i=0;i<process.parallelisation->repartition_sst[MPIsource].size();i++) {
        int num = process.parallelisation->repartition_sst[MPIsource][i];
        if (process.nom_calcul=="latin") {
            for(unsigned pt=1;pt<S[num].t.size();pt++) {
                S[num].t[pt].q=vectorecv[range(repere+pt*S[num].t[pt].q.size(),repere+(pt+1)*S[num].t[pt].q.size())];
            }
            repere+=S[num].t.size()*S[num].t[1].q.size();
        } else {
            for(unsigned pt=1;pt<S[num].t_post.size();pt++) {
                S[num].t_post[pt].q=vectorecv[range(repere+pt*S[num].t_post[pt].q.size(),repere+(pt+1)*S[num].t_post[pt].q.size())];
            }
            repere+=S[num].t_post.size()*S[num].t_post[1].q.size();
        }
    }
}

template <class T1, class T2>
void SendInterMaster(T1 &vecintertoexchange,T2 &Inter,MPI_Request &request,Process &process, Vec<TYPEREEL> &vectosend) {
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
void RecvInterMaster(T1 &MPIsource,T2 &Inter, Process &process) {
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
void SendInterMasterPost(T1 &vecintertoexchange,T2 &Inter,MPI_Request &request,Process &process, Vec<TYPEREEL> &vectosend) {
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
void RecvInterMasterPost(T1 &MPIsource,T2 &Inter, Process &process) {
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

template <class T1,class T2,class T3,class T4>
void SendRecvInterMaster(T1 &intertoexchangeformaster, T2 &Inter,T3 &S,T4 &SubS, Process &process) {
    Vec<TYPEREEL> vectosend1,vectosend2,vectosend3,vectosend4;
    MPI_Request request,request2;

    if(not process.parallelisation->is_master_cpu()) {
        if (process.nom_calcul=="latin")
            SendInterMaster(intertoexchangeformaster,Inter,request,process,vectosend1);
        if (process.nom_calcul=="incr")
            SendInterMasterPost(intertoexchangeformaster,Inter,request,process,vectosend1);
        SendSST(SubS, request2,process, vectosend2);
    }
    if(process.parallelisation->is_master_cpu()) {
        int nbtorecv = (process.parallelisation->size-1)*2;
        while( nbtorecv != 0) {
            //On regarde si on a reçu un message de quelqu'un (MPI_ANY_SOURCE).
            int flag=0;
            MPI_Status status1;
            MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag, &status1);
            //Si on a recu un message alors on va le traiter
            if (flag == 1) {
                if (status1.MPI_TAG == 201) {
                    if (process.nom_calcul=="latin")
                        RecvInterMaster(status1.MPI_SOURCE,Inter,process);
                    if (process.nom_calcul=="incr")
                        RecvInterMasterPost(status1.MPI_SOURCE,Inter,process);
                }
                if (status1.MPI_TAG == 204) {
                    RecvSST(S,status1.MPI_SOURCE,process);
                }
                nbtorecv-=1;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for( unsigned i=0;i< vecintertoexchangebypro.size() ;i++ )
        MPI_Request_free(&request[i]);
}
#endif //SENDRECVTOMASTER_H
