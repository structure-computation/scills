#ifndef PARAM_TEMPS_H
#define PARAM_TEMPS_H
using namespace LMT;
using namespace std;

/** \ingroup Paramètres_calcul
Parametres temporels
*/

struct STEP{
    STEP(){
        dt=1;           ///< pas de temps dans le step
        t_ini=0;        ///< piquet de temps initial du step
        t_fin=1;        ///< piquet de temps final du step
        nb_time_step=1; ///< nombre de pas de temps dans le step
        pt_cur=0;       ///< pas de temps courant dans le step
    }
    double dt, t_ini, t_fin;
    int nb_time_step;
    int pt_cur;
};

struct TEMPS{
  TEMPS(){
     pt=0;
     pt_cur=0;
     dt=1;
     theta = 1;
     nbpastemps=1;
     type_de_calcul="stat";
     nb_step=1;
  }
  int pt;            ///< piquet de temps de calcul (0 ou 1)
  int pt_cur;   ///< piquet de temps (ou intervalle de temps) courant global
  int step_cur;         ///< step courant
  double dt;         ///< pas de temps
  double theta;   ///< paramètre de la theta_methode
  unsigned nbpastemps; ///< nb de pas de temps total
  string type_de_calcul; ///< type de calcul, choix entre : stat (statique (1 pas de temps)), qstat (plusieurs pas de temps)
  Vec<STEP> time_step;  ///< steps de calcul
  int nb_step;          ///< nombre de step de calcul
};

#endif // PARAM_TEMPS_H

