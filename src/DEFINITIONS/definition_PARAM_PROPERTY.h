#ifndef PARAM_PROPERTY_H
#define PARAM_PROPERTY_H
using namespace LMT;
using namespace std;

/** Parametres materiau definis de maniere globale sur la structure 
*/
struct PROPERTY{
  PROPERTY(){
     deltaT=0.; 
     //epshorsplan.set(0.);
  }
  double deltaT; /// variation de temperature
  //Vec<double,3> epshorsplan; // parametre pour les def planes generalisees

};

#endif // PARAM_PROPERTY_H

