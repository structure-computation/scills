#ifndef GLOB_H
#define GLOB_H
#include "util/solveLDL.h"

using namespace LMT;
using namespace std;

//********************************
// Operateurs probleme macro
//********************************
/** \ingroup Operateurs
\brief Operateurs macro parametrables

La classe est parametrable par :
dim_ dimension du probleme
TT_ : type de flottant
*/
template<unsigned dim_, class TT_> struct Glob
{
  Glob(){max_erreur=0;}
  static const unsigned dim = dim_; ///< dimension 2d ou 3d
  typedef TT_ T; ///< type connu de l'exterieur

  Vec<T> bigW; ///< vecteur solution du probleme macro
  Vec<T> bigF; ///< vecteur second membre macro
  //Inv<T, Sym<>, SparseLine<> > bigL; // factorisation de choleski de la matrice macro 
  LDL_solver l; ///< solver base sur la librairie LDL
  
  Vec<unsigned> repddlMbloq; ///< reperage des ddls macro bloques
  T coefpenalisation; ///< coef de penalisation de la matrice macro
  T max_erreur; /// max de au cours des iterations norm_2(Global.bigW)/Global.bigW.size()

  void allocations(Param &process); ///< defini dans allocate.h
  void resolmacro();///< résolution du problème macro 
};

template<unsigned dim_, class TT_> void Glob<dim_,TT_>::resolmacro(){ 
  //Mise de bigF dans bigW pour obtenir la solution
  bigW=bigF;
  //resolution
  l.solve( bigW );
   //mise a zero des ddls bloques
  bigW[repddlMbloq]=0.0;
}

#endif // GLOB_H
