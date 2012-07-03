#ifndef MACROPROBLEM_H
#define MACROPROBLEM_H

#include "main_typedef.h"
#include "../../LMT/include/util/solveLDL.h"

//********************************
// Operateurs probleme macro
//********************************
/** \ingroup Operateurs
\brief Operateurs macro parametrables
*/
struct MacroProblem {
    //attributs==============================================================================================
    Vector bigW;    /// vecteur solution du probleme macro
    Vector bigF;    /// vecteur second membre macro
    LDL_solver l;   /// solver base sur la librairie LDL

    LMT::Vec<unsigned> repddlMbloq; /// reperage des ddls macro bloques
    Scalar coefpenalisation;        /// coef de penalisation de la matrice macro
    Scalar max_erreur;              /// max de au cours des iterations norm_2(Global.bigW)/Global.bigW.size()

    //methodes===============================================================================================
    MacroProblem();
    void allocations(int size_macro_pb);
    void resolmacro();
};

#endif // MACROPROBLEM_H
