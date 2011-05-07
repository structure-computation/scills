#ifndef PARAM_TEMPS_H
#define PARAM_TEMPS_H
using namespace LMT;
using namespace std;

/** \ingroup Paramètres_calcul
Parametres temporels
*/
struct TEMPS{
  TEMPS(){
     pt=0;
     pt_cur=0;
     dt=1;
     theta = 1;
     nbpastemps=1;
     type_de_calcul="stat";
  }
  int pt;            ///< piquet de temps de calcul
  int pt_cur;   ///< piquet de temps (ou intervalle de temps) courant
  double dt;         ///< pas de temps
  double theta;   ///< paramètre de la theta_methode
  unsigned nbpastemps; ///< nb de pas de temps total
  string type_de_calcul; ///< type de calcul, choix entre : stat (statique (1 pas de temps)), qstat (plusieurs pas de temps)
};

#endif // PARAM_TEMPS_H

