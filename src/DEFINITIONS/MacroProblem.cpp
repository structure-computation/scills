#include "MacroProblem.h"
#include "MultiScaleData.h"

MacroProblem::MacroProblem(){
    max_erreur=0;
}

void MacroProblem::allocations(int size_macro_pb){
    bigW.resize(size_macro_pb);
    bigW.set(0.0);
    bigF.resize(size_macro_pb);
    bigF.set(0.0);
}

void MacroProblem::resolmacro(){ 
    //Mise de bigF dans bigW pour obtenir la solution
    bigW=bigF;
    //resolution
    l.solve( bigW );
    //mise a zero des ddls bloques
    bigW[repddlMbloq]=0.0;
}