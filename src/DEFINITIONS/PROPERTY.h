#ifndef PARAM_PROPERTY_H
#define PARAM_PROPERTY_H
using namespace LMT;

/** Parametres materiau definis de maniere globale sur la structure 
*/
struct PROPERTY{
  PROPERTY(){
     deltaT=0.; 
     //epshorsplan.set(0.);
  }
  TYPEREEL deltaT; /// variation de temperature
  //Vec<TYPEREEL,DIM> epshorsplan; // parametre pour les def planes generalisees

};

#endif // PARAM_PROPERTY_H

