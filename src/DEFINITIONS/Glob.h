#ifndef GLOB_H
#define GLOB_H
#include "../../LMT/include/util/solveLDL.h"
#include "Param.h"
#include "MULTI.h"

using namespace LMT;

//********************************
// Operateurs probleme macro
//********************************
/** \ingroup Operateurs
\brief Operateurs macro parametrables
*/
struct Glob {
    Vec<TYPEREEL> bigW; ///< vecteur solution du probleme macro
    Vec<TYPEREEL> bigF; ///< vecteur second membre macro
    LDL_solver l; ///< solver base sur la librairie LDL

    Vec<unsigned> repddlMbloq; ///< reperage des ddls macro bloques
    TYPEREEL coefpenalisation; ///< coef de penalisation de la matrice macro
    TYPEREEL max_erreur; /// max de au cours des iterations norm_2(Global.bigW)/Global.bigW.size()

    Glob(){
        max_erreur=0;
    }
    
    void allocations(Param &process){
        bigW.resize(process.multiscale->sizeM);
        bigW.set(0.0);
        bigF.resize(process.multiscale->sizeM);
        bigF.set(0.0);
    }
    
    void resolmacro(){ 
        //Mise de bigF dans bigW pour obtenir la solution
        bigW=bigF;
        //resolution
        l.solve( bigW );
        //mise a zero des ddls bloques
        bigW[repddlMbloq]=0.0;
    }
};

#endif // GLOB_H
