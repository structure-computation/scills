#ifndef DEFINITION_FCTS
#define DEFINITION_FCTS
//fichiers contenant les differentes declarations des fonctions utilisees dans le programme principal
#include "GeometryUser.h"
#include "DataUser.h"
#include "FieldStructureUser.h"

#include "Param.h"
#include "definition_SST_time.h"
#include "definition_INTER_time.h"
#include "Boundary.h"
#include "containers/vec_mt.h"
#include "containers/vecpointedvalues.h"

using namespace LMT;

void multiscale_geometry_mesh(DataUser                          &data_user, 
                              GeometryUser                      &geometry_user,
                              Vec<Sst>                          &S,
                              Vec<Interface>                    &Inter, 
                              Param                             &process, 
                              Vec<Boundary>                     &CL,
                              Vec<VecPointedValues<Sst> >       &Stot,
                              Vec<VecPointedValues<Sst> >       &SubS,
                              Vec<VecPointedValues<Interface> > &SubI) __attribute__((noinline));

void assignation_materials_property_SST(DataUser &data_user, 
                                        Vec<SstCarac>      &matprops, 
                                        Vec<Sst>           &S, 
                                        Vec<Interface>     &Inter,
                                        Param              &process, 
                                        FieldStructureUser &field_structure_user) __attribute__((noinline));


void assignation_materials_property_INTER(DataUser           &data_user, 
                                          Vec<Interface>     &Inter, 
                                          Vec<Sst>           &S, 
                                          Param              &process, 
                                          FieldStructureUser &field_structure_user) __attribute__((noinline));


// template<class TV1, class TV2, class TV3, class TV4,class TV5> void multiscale_geometry_mesh(const XmlNode &n,TV1 &S,TV2 &Inter, Param &process, TV5 &CL,TV3 &Stot,TV3 &SubS,TV4 &SubI) ;
// template<class TV1, class TV2> void assignation_materials_property_SST(const XmlNode &n, TV1 &S, TV2 &Inter,Param &process) __attribute__((noinline));
// template<class TV1, class TV2> void assignation_materials_property_INTER(const XmlNode &n, TV2 &Inter, TV1 &S, Param &process) __attribute__((noinline));
template <class TV1, class TV2, class TV3, class TV4>  void multiscale_operateurs(TV1 &S, TV1 &SubS,TV2 &Inter, TV3 &SubI,Param &process,  TV4 &Global, DataUser &data_user) ;
template <class TV1,class TV2, class TV3, class TV4, class TV5, class TV6> void multiscale_iterate_latin(TV1 &S,TV2 &SubS, TV3 &Inter, TV4 &SubI, Param &process, TV5 &Global, TV6 &CL, DataUser &data_user);
template <class TV1,class TV2, class TV3, class TV4, class TV5, class TV6> void multiscale_iterate_incr(TV1 &S,TV2 &SubS, TV3 &Inter, TV4 &SubI, Param &process, TV5 &Global, TV6 &CL, DataUser &data_user, GeometryUser &geometry_user, FieldStructureUser &field_structure_user);

#endif DEFINITION_FCTS