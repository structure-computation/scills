// fichiers de definition des classes
#include "definition_PARAM.h"
#include "definition_GLOB.h"
#include "definition_SST_time.h"
#include "definition_INTER_time.h"
#include "definition_CL.h"


template< class TV1 > 
Sst<DIM,TYPEREEL>* find_sst(TV1 &S,int id_) {
    for (int i_group=0; i_group<S.size(); i_group++) {
        if (S[i_group].num == id_) {
            return &S[i_group];
            break;
        }
    }
}

template< class TV1 > 
int find_index_sst(TV1 &S, int id_) {
    for (int i_group=0; i_group<S.size(); i_group++) {
        if (S[i_group].num == id_) {
            return i_group;
            break;
        }
    }
}