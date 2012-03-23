#include "../COMPUTE/DataUser.h"
#include "Boundary.h"

using namespace Metil;
using namespace Codegen;
using namespace LMT;

/** \ingroup Conditions_limites
\brief Lecture des conditions aux limites
 
 champ parametres : contient box, comportement
 
 champ fct_spatiale : contient la fct spatiale appliquee sur cette zone : entrer une fct de x ou y en 2d et x,y ou x,z ou y,z ou x,y,z  en 3d
 
 champ fct_temporelle : contient la fct temporelle appliquee sur l'intervalle de temps donne, rentrer une fonction de t
 
 On multiplie dans prelocalstage.h la valeur de la fonction spatiale en un noeud equivalent par la fonction temporel (en statique comme en quasistatique)
*/

void read_CL(DataUser &data_user, Vec<Boundary > &CL, Param &process) {
    const unsigned nbCL = data_user.behaviour_bc.size();
    const unsigned nbStep = process.temps->nb_step;
    unsigned nbfct_temporelle;
    CL.resize(nbCL);
    for(unsigned i=0;i<nbCL;++i) {
        CL[i].id = data_user.behaviour_bc[i].id;
        CL[i].comp = data_user.behaviour_bc[i].type;
        std::cout << "data_user.behaviour_bc[i].type = " << data_user.behaviour_bc[i].type << std::endl;
        if (CL[i].comp=="depl_nul") {
            CL[i].comp = "depl";
            data_user.behaviour_bc[i].type = "depl";
        }
        else if (CL[i].comp=="sym") {
            CL[i].fcts_spatiales.resize(nbStep);
            CL[i].fcts_temporelles.resize(nbStep);
            CL[i].intervalles_temps.resize(nbStep);
            for(unsigned i_step=0;i_step<nbStep;i_step++){
                CL[i].intervalles_temps[i_step][0]=0;
                CL[i].intervalles_temps[i_step][1]=100000;
                CL[i].fcts_spatiales[i_step].resize(data_user.dim,"0");
                CL[i].fcts_temporelles[i_step]="1";
            }
        }
        else{
            if(process.temps->type_de_calcul=="stat") {
                CL[i].fcts_spatiales.resize(nbStep);
                CL[i].fcts_temporelles.resize(nbStep);
                CL[i].intervalles_temps.resize(nbStep);                
                for(unsigned i_step=0;i_step<nbStep;i_step++){
                    CL[i].intervalles_temps[i_step][0]=0;
                    CL[i].intervalles_temps[i_step][1]=100000;
                    CL[i].fcts_spatiales[i_step].resize(data_user.dim,"0");
                    for(int d=0; d<data_user.dim; d++){
                        CL[i].fcts_spatiales[i_step][d] = data_user.behaviour_bc[i].step[i_step].CL_step_prop[d];
                    }
                    CL[i].fcts_temporelles[i_step]="1";
                }    
            }else if(process.temps->type_de_calcul=="Qstat") {
                CL[i].fcts_spatiales.resize(nbStep);
                CL[i].fcts_temporelles.resize(nbStep);
                CL[i].intervalles_temps.resize(nbStep);                 
                for(unsigned i_step=0;i_step<nbStep;i_step++){
                    PRINT(i_step);
                    CL[i].intervalles_temps[i_step][0]=process.temps->time_step[i_step].t_ini;
                    CL[i].intervalles_temps[i_step][1]=process.temps->time_step[i_step].t_fin;
                    CL[i].fcts_spatiales[i_step].resize(data_user.dim,"0");
                    for(int d=0; d<data_user.dim; d++){
                        CL[i].fcts_spatiales[i_step][d] = data_user.behaviour_bc[i].step[i_step].CL_step_prop[d];
                    }
                    CL[i].fcts_temporelles[i_step] = data_user.behaviour_bc[i].step[i_step].CL_step_prop[3];
                }
            }else {
                std::cout << "Mauvais type de calcul" << std::endl;
                assert(0);
            }
        }
    }
};

// modification des valeurs des boites de CL en fonction du parametre d'echelle
void modif_CL_scale(Vec<Boundary > &CL, TYPEREEL &scale) {
    for(unsigned i=0;i<CL.size();++i)
        CL[i].box=CL[i].box*scale;
};
