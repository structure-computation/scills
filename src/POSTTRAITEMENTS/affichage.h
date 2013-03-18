#ifndef AFFICHAGE_H
#define AFFICHAGE_H


//#include "../DEFINITIONS/Process.h"
//#include "../DEFINITIONS/Sst.h"
//#include "../DEFINITIONS/Interface.h"
//#include "../COMPUTE/DataUser.h"

//#include "../../LMT/include/containers/vec.h"
//#include "../../LMT/include/containers/vecpointedvalues.h"
//using namespace LMT;

#include "../DEFINITIONS/structure_typedef.h"


template <class TV3,class TV4, class TV1> void affichage_maillage(TV3 &SubS, TV4 &Inter,TV1 & S, Process &process, DataUser &data_user) __attribute__((noinline));
void affichage_resultats(PointedSubstructures &S,  Process &process, DataUser &data_user) __attribute__((noinline));
template <class TV3> void affichage_depl_pt(TV3 &S, Process &process) __attribute__((noinline));
template <class TV3, class TV2> void affichage_var_inter(TV3 &S, TV2 &Inter, Process &process) __attribute__((noinline));
template <class TV1,class TV2> void affichage_inter_data(TV2 &Inter, TV1 &S, Process &process) __attribute__((noinline));
void affichage_resultats_inter(PointedInterfaces &Inter, Substructures &S , Process &process, DataUser &data_user) __attribute__((noinline));
template <class TV3,class TV2> void affichage_energie(TV3 &S,TV2 &Inter, Process &process, DataUser &data_user) __attribute__((noinline));

void create_pvd_geometry(PointedSubstructures &S, Substructures &SS, VecInterfaces &Inter,Process &process) __attribute__((noinline));
template <class TV1,class TV2, class TV3> void create_pvd_geometry_sst_by_group(TV1 &S, TV3 &SS, TV2 &Inter,Process &process) __attribute__((noinline));
void create_pvd_results(PointedSubstructures &S, Substructures &SS, VecInterfaces &Inter,Process &process) __attribute__((noinline));

void affichage_resultats_temps(Process &process);
void affichage_inter_temps(Process &process);

#endif //AFFICHAGE_H

