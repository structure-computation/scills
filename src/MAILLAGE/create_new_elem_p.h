#ifndef CREATE_NEW_ELEM_P_H
#define CREATE_NEW_ELEM_P_H
#include "mesh/tetra.h"
#include "mesh/tetra_10.h"
#include "mesh/wedge.h"
#include "mesh/wedge_15.h"
#include "mesh/hexa.h"
#include "mesh/hexa_20.h"

namespace LMT{
///Fonction generique
template<class TE, class TM> void create_new_elem_p(TE &e, TM &m, Vec<unsigned> &rep_nodes){}
///Fonctions specialisees
template<class TNB,class TN,class TD,unsigned NET, class TM> void create_new_elem_p(Element<Tetra,TNB,TN,TD,NET> &e, TM &m, Vec<unsigned> &rep_nodes){ m.add_element(Tetra_10(),DefaultBehavior(),rep_nodes.ptr() ); } 
template<class TNB,class TN,class TD,unsigned NET, class TM> void create_new_elem_p(Element<Wedge,TNB,TN,TD,NET> &e, TM &m, Vec<unsigned> &rep_nodes){ m.add_element(Wedge_15(),DefaultBehavior(),rep_nodes.ptr() ); } 
template<class TNB,class TN,class TD,unsigned NET, class TM> void create_new_elem_p(Element<Hexa,TNB,TN,TD,NET> &e, TM &m, Vec<unsigned> &rep_nodes){ m.add_element(Hexa_20(),DefaultBehavior(),rep_nodes.ptr() ); } 
template<class TNB,class TN,class TD,unsigned NET, class TM> void create_new_elem_p(Element<Triangle,TNB,TN,TD,NET> &e, TM &m, Vec<unsigned> &rep_nodes){ m.add_element(Triangle_6(),DefaultBehavior(),rep_nodes.ptr() ); } 
template<class TNB,class TN,class TD,unsigned NET, class TM> void create_new_elem_p(Element<Quad,TNB,TN,TD,NET> &e, TM &m, Vec<unsigned> &rep_nodes){ m.add_element(Quad_8(),DefaultBehavior(),rep_nodes.ptr() ); } 
template<class TNB,class TN,class TD,unsigned NET, class TM> void create_new_elem_p(Element<Bar,TNB,TN,TD,NET> &e, TM &m, Vec<unsigned> &rep_nodes){ m.add_element(Bar_3(),DefaultBehavior(),rep_nodes.ptr() ); } 
};

#endif //CREATE_NEW_ELEM_P_H

