#ifndef GLOB_H
#define GLOB_H
#include "util/solveLDL.h"
#include "definition_PARAM_MULTI.h"

using namespace LMT;

//********************************
// Operateurs probleme macro
//********************************
/** \ingroup Operateurs
\brief Operateurs macro parametrables

La classe est parametrable par :
dim_ dimension du probleme
TT_ : type de flottant
*/
struct Glob
{
  Glob(){max_erreur=0;}

  Vec<TYPEREEL> bigW; ///< vecteur solution du probleme macro
  Vec<TYPEREEL> bigF; ///< vecteur second membre macro
  //Inv<TYPEREEL, Sym<>, SparseLine<> > bigL; // factorisation de choleski de la matrice macro 
  LDL_solver l; ///< solver base sur la librairie LDL
  
  Vec<unsigned> repddlMbloq; ///< reperage des ddls macro bloques
  TYPEREEL coefpenalisation; ///< coef de penalisation de la matrice macro
  TYPEREEL max_erreur; /// max de au cours des iterations norm_2(Global.bigW)/Global.bigW.size()

  void allocations(Param &process)
  {
      bigW.resize(process.multiscale->sizeM);
      bigW.set(0.0);
      bigF.resize(process.multiscale->sizeM);
      bigF.set(0.0);
  }
  void resolmacro(){///< résolution du problème macro
    //Mise de bigF dans bigW pour obtenir la solution
    bigW=bigF;
    //resolution
    l.solve( bigW );
    //mise a zero des ddls bloques
    bigW[repddlMbloq]=0.0;
  }
};

#endif // GLOB_H
