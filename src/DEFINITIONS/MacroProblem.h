#ifndef MACROPROBLEM_H
#define MACROPROBLEM_H
#include "../../LMT/include/util/solveLDL.h"
#include "Process.h"
#include "MultiScaleParameters.h"

using namespace LMT;

//********************************
// Operateurs probleme macro
//********************************
/** \ingroup Operateurs
\brief Operateurs macro parametrables
*/
struct MacroProblem {
    //attributs==============================================================================================
    Vec<TYPEREEL> bigW; ///< vecteur solution du probleme macro
    Vec<TYPEREEL> bigF; ///< vecteur second membre macro
    LDL_solver l;       ///< solver base sur la librairie LDL

    Vec<unsigned> repddlMbloq; ///< reperage des ddls macro bloques
    TYPEREEL coefpenalisation; ///< coef de penalisation de la matrice macro
    TYPEREEL max_erreur;       ///< max de au cours des iterations norm_2(Global.bigW)/Global.bigW.size()

    //methodes===============================================================================================
    MacroProblem();
    void allocations(Process &process);
    void resolmacro();
};

#endif // GLOB_H
