#ifndef AFFICHAGE_H
#define AFFICHAGE_H
#include "definition_PARAM_AFFICHAGE.h"
#include "containers/vecpointedvalues.h"

//permet un pont entre les differents cpp et les fonctions d'affichage sans avoir a tout recompiler
using namespace LMT;
using namespace std;

template <class TV3,class TV4, class TV1> void affichage_maillage(TV3 &SubS, TV4 &Inter,TV1 & S, Param &process) __attribute__((noinline));
template <class TV3> void affichage_resultats(TV3 &S,  Param &process, DataUser &data_user) __attribute__((noinline));
template <class TV3> void affichage_depl_pt(TV3 &S, Param &process) __attribute__((noinline));
template <class TV3, class TV2> void affichage_var_inter(TV3 &S, TV2 &Inter, Param &process) __attribute__((noinline));
template <class TV1,class TV2> void affichage_inter_data(TV2 &Inter, TV1 &S, Param &process) __attribute__((noinline));
template <class TV1,class TV4> void affichage_resultats_inter(TV4 &Inter, TV1 &S , Param &process) __attribute__((noinline));
template <class TV3,class TV2> void affichage_energie(TV3 &S,TV2 &Inter, Param &process, DataUser &data_user) __attribute__((noinline));

void affichage_resultats_temps(Param &process);
void affichage_inter_temps(Param &process);

#endif //AFFICHAGE_H

