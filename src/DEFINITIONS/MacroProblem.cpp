#include "MacroProblem.h"

MacroProblem::MacroProblem(){
    max_erreur=0;
}

void MacroProblem::allocations(Process &process){
    bigW.resize(process.multiscale->sizeM);
    bigW.set(0.0);
    bigF.resize(process.multiscale->sizeM);
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