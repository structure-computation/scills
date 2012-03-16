// fichiers de definition des classes
#include "definition_PARAM.h"
#include "definition_GLOB.h"
#include "definition_SST_time.h"
#include "definition_INTER_time.h"
#include "definition_CL.h"

// trouver une sst à paritr de son id----------------------------------------------
template< class TV1 > 
Sst* find_sst(TV1 &S,int id_) {
    for (int i_group=0; i_group<S.size(); i_group++) {
        if (S[i_group].id == id_) {
            return &S[i_group];
            break;
        }
    }
}

template< class TV1 > 
int find_index_sst(TV1 &S, int id_) {
    for (int i_group=0; i_group<S.size(); i_group++) {
        if (S[i_group].id == id_) {
            return i_group;
            break;
        }
    }
}

// trouver une interface à paritr de son id----------------------------------------------
template< class TV1 > 
Interface* find_inter(TV1 &Inter,int id_) {
    for (int i_group=0; i_group<Inter.size(); i_group++) {
        if (Inter[i_group].id == id_) {
            return &Inter[i_group];
            break;
        }
    }
}

template< class TV1 > 
int find_index_inter(TV1 &Inter,int id_) {
    for (int i_group=0; i_group<Inter.size(); i_group++) {
        if (Inter[i_group].id == id_) {
            return i_group;
            break;
        }
    }
}

// trouver un group d' interfaces à paritr de leur edge_id----------------------------------------------
template< class TV1 > 
Vec< Interface* > find_group_edges(TV1 &Inter, int edge_id_) {
    Vec< Interface* > group_edges_temp;
    for (int i_group=0; i_group<Inter.size(); i_group++) {
        if (Inter[i_group].edge_id == edge_id_) {
            group_edges_temp.push_back(&Inter[i_group]);
        }
    }
    return group_edges_temp;
}

template< class TV1 > 
Vec< int > find_index_group_edges(TV1 &Inter, int edge_id_) {
    Vec< int > group_edges_temp;
    for (int i_group=0; i_group<Inter.size(); i_group++) {
        if (Inter[i_group].edge_id == edge_id_) {
            group_edges_temp.push_back(i_group);
        }
    }
    return group_edges_temp;
}
