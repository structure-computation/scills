#include "DataUser.h"

using namespace Metil;
using namespace Codegen;
using namespace std;
using namespace LMT;

/** \ingroup Conditions_limites
\brief Lecture des conditions aux limites
 
 champ parametres : contient box, comportement
 
 champ fct_spatiale : contient la fct spatiale appliquee sur cette zone : entrer une fct de x ou y en 2d et x,y ou x,z ou y,z ou x,y,z  en 3d
 
 champ fct_temporelle : contient la fct temporelle appliquee sur l'intervalle de temps donne, rentrer une fonction de t
 
 On multiplie dans prelocalstage.h la valeur de la fonction spatiale en un noeud equivalent par la fonction temporel (en statique comme en quasistatique)
*/

template<class BOUNDARY>
void read_CL(DataUser &data_user, Vec<BOUNDARY > &CL, Param &process) {
};

template<class BOUNDARY>
void read_CL(const XmlNode &n, Vec<BOUNDARY > &CL, Param &process) {
    unsigned nbCL = n.nb_elements("CL");
    unsigned nbfct_temporelle;
    CL.resize(nbCL);
    string num_sst;
    for(unsigned i=0;i<nbCL;++i) {
        XmlNode ncl= n.get_element("CL",i);
        XmlNode npar= ncl.get_element("parametres");
        npar.get_attribute("comp",CL[i].comp);
        Vec<typename BOUNDARY::T> box;
        npar.get_attribute("box",box);
        npar.get_attribute("num_sst",num_sst,"");
        Vec<string> vec_sst_num;
        vec_sst_num=tokenize(num_sst,' ');
        CL[i].sst_num.resize(vec_sst_num.size());
        for( unsigned ii=0;ii<vec_sst_num.size() ;ii++ ){
            CL[i].sst_num[ii]=atoi(vec_sst_num[ii].c_str());
        }

        CL[i].box[0]=box[range((int)BOUNDARY::dim)];
        CL[i].box[1]=box[range((int)BOUNDARY::dim,(int)(2*BOUNDARY::dim))];

        //si le comportement est du type depl_normal : assignation de la valeur de la fct spatiale au deplacement donne
        //si le comportement est du type symetrie : assignation de 0
        if (CL[i].comp=="sym") {
            CL[i].fcts_spatiales.set("0");
            CL[i].fcts_temporelles.resize(1);
            CL[i].fcts_temporelles[0].resize(1);
            CL[i].intervalles_temps.resize(1);
            CL[i].intervalles_temps[0]=Vec<typename BOUNDARY::T,2>(0,100000);
            CL[i].fcts_temporelles[0]="1";
        } else if (CL[i].comp=="periodique") {
            CL[i].fcts_spatiales.set("0");
            CL[i].fcts_temporelles.resize(1);
            CL[i].fcts_temporelles[0].resize(1);
            CL[i].intervalles_temps.resize(1);
            CL[i].intervalles_temps[0]=Vec<typename BOUNDARY::T,2>(0,100000);
            CL[i].fcts_temporelles[0]="1";
            Vec<typename BOUNDARY::T> box1;
            npar.get_attribute("box1",box1);
            CL[i].box1[0]=box1[range((int)BOUNDARY::dim)];
            CL[i].box1[1]=box1[range((int)BOUNDARY::dim,(int)(2*BOUNDARY::dim))];
        } else if (CL[i].comp=="depl_normal") {

            cout << "ATTENTION : une condition limite en deplacement normal n'est valable que pour des surfaces planes" << endl;
//             XmlNode nfs= ncl.get_element("fct_spatiale");
//             string valfcts;
//             nfs.get_attribute("fonction",valfcts);
//             CL[i].fcts_spatiales.set("0");
//             CL[i].fcts_spatiales[0]=valfcts;
            //lecture des fcts temporelles definies pour un intervalle de temps donne
            if(process.temps->type_de_calcul=="Qstat") {
                CL[i].fcts_spatiales.resize(process.temps->nb_step);
                CL[i].fcts_temporelles.resize(process.temps->nb_step);
                CL[i].intervalles_temps.resize(process.temps->nb_step);                
                for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
                    CL[i].intervalles_temps[0]=process.temps->time_step[i_step].t_ini;
                    CL[i].intervalles_temps[1]=process.temps->time_step[i_step].t_fin;
                    XmlNode nfs= ncl.get_element("fct_step",i_step);
                    string valfcts;
                    nfs.get_attribute("fonction_spatiale",valfcts);
                    CL[i].fcts_spatiales[i_step]=tokenize(valfcts,';');
                    nfs.get_attribute("fonction_temporelle",valfcts);
                    CL[i].fcts_temporelles[i_step]=valfcts;
                }
//                 nbfct_temporelle= ncl.nb_elements("fct_temporelle");
//                 CL[i].fcts_temporelles.resize(nbfct_temporelle);
//                 CL[i].intervalles_temps.resize(nbfct_temporelle);
//                 for(unsigned j=0;j<nbfct_temporelle;++j) {
//                     XmlNode nft= ncl.get_element("fct_temporelle",j);
//                     nft.get_attribute("intervalle",CL[i].intervalles_temps[j]);
//                    ////modif DAVID 02-09-2007 
//                     string valfctstps;
//                     nft.get_attribute("fonction",valfctstps);
//                     CL[i].fcts_temporelles[j]=valfctstps;
//                     //nft.get_attribute("fonction",CL[i].fcts_temporelles[j]);
//                     ////fin modif DAVID 02-09-2007 
//                 }
            } else if(process.temps->type_de_calcul=="stat") {
                CL[i].fcts_spatiales.resize(process.temps->nb_step);
                CL[i].fcts_temporelles.resize(process.temps->nb_step);
                CL[i].intervalles_temps.resize(process.temps->nb_step);                
                for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
                    CL[i].intervalles_temps[0]=0;
                    CL[i].intervalles_temps[1]=100000;
                    XmlNode nfs= ncl.get_element("fct_step",i_step);
                    string valfcts;
                    nfs.get_attribute("fonction_spatiale",valfcts);
                    CL[i].fcts_spatiales[i_step]=tokenize(valfcts,';');
                    CL[i].fcts_temporelles[i_step]="1";
                }
//                 XmlNode nfs= ncl.get_element("fct_spatiale");
//                 string valfcts;
//                 nfs.get_attribute("fonction",valfcts);
//                 CL[i].fcts_spatiales.set("0");
//                 CL[i].fcts_spatiales[0]=valfcts;
//                 CL[i].fcts_temporelles.resize(1);
//                 CL[i].intervalles_temps.resize(1);
//                 CL[i].fcts_temporelles[0].resize(1);
//                 CL[i].intervalles_temps[0]=Vec<typename BOUNDARY::T,2>(0,100000);
//                 CL[i].fcts_temporelles[0]="1";
            } else {
                std::cout << "Mauvais type de calcul" << std::endl;
                assert(0);
            }
        } else { //sinon lecture des differentes fcts_spatiales
            //lecture de la fct spatiale sur la zone donnee : separation des valeurs selon x,y,z par des ; ne pas mettre d'espace du tout et stockage dans fcts_spatiales vecteur de string
//             XmlNode nfs= ncl.get_element("fct_spatiale");
//             string valfcts;
//             nfs.get_attribute("fonction",valfcts);
//             CL[i].fcts_spatiales=tokenize(valfcts,';');
            //lecture des fcts temporelles definies pour un intervalle de temps donne
            if(process.temps->type_de_calcul=="Qstat") {
                CL[i].fcts_spatiales.resize(process.temps->nb_step);
                CL[i].fcts_temporelles.resize(process.temps->nb_step);
                CL[i].intervalles_temps.resize(process.temps->nb_step);                
                for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
                    CL[i].intervalles_temps[0]=process.temps->time_step[i_step].t_ini;
                    CL[i].intervalles_temps[1]=process.temps->time_step[i_step].t_fin;
                    XmlNode nfs= ncl.get_element("fct_step",i_step);
                    string valfcts;
                    nfs.get_attribute("fonction_spatiale",valfcts);
                    CL[i].fcts_spatiales[i_step]=tokenize(valfcts,';');
                    nfs.get_attribute("fonction_temporelle",valfcts);
                    CL[i].fcts_temporelles[i_step]=valfcts;
                }
 /*               XmlNode nfs= ncl.get_element("fct_spatiale");
                string valfcts;
                nfs.get_attribute("fonction",valfcts);
                CL[i].fcts_spatiales=tokenize(valfcts,';');
                
                nbfct_temporelle= ncl.nb_elements("fct_temporelle");
                CL[i].fcts_temporelles.resize(nbfct_temporelle);
                
                for(unsigned j=0;j<nbfct_temporelle;++j) {
                    XmlNode nft= ncl.get_element("fct_temporelle",j);
                    nft.get_attribute("intervalle",CL[i].intervalles_temps[j]);
                    ////modif DAVID 02-09-2007 
                    string valfctstps;
                    nft.get_attribute("fonction",valfctstps);
                    Vec<string> vecvalfctstps;
                    vecvalfctstps=tokenize(valfctstps,';');
                    CL[i].fcts_temporelles[j].resize(vecvalfctstps.size());
                    CL[i].ft.resize(vecvalfctstps.size());
                    CL[i].fcts_temporelles[j]=vecvalfctstps;
                     //nft.get_attribute("fonction",CL[i].fcts_temporelles[j]);
                    ////fin modif DAVID 02-09-2007 
                }*/
            } else if(process.temps->type_de_calcul=="stat") {
             
                CL[i].fcts_spatiales.resize(process.temps->nb_step);
                CL[i].fcts_temporelles.resize(process.temps->nb_step);
                CL[i].intervalles_temps.resize(process.temps->nb_step);                
                for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
                    CL[i].intervalles_temps[0]=0;
                    CL[i].intervalles_temps[1]=100000;
                    XmlNode nfs= ncl.get_element("fct_step",i_step);
                    string valfcts;
                    nfs.get_attribute("fonction_spatiale",valfcts);
                    CL[i].fcts_spatiales[i_step]=tokenize(valfcts,';');
                    CL[i].fcts_temporelles[i_step]="1";
                }
//                 XmlNode nfs= ncl.get_element("fct_spatiale");
//                 string valfcts;
//                 nfs.get_attribute("fonction",valfcts);
//                 CL[i].fcts_spatiales.resize(1);
//                 CL[i].fcts_spatiales[0]=tokenize(valfcts,';');
//                  CL[i].fcts_temporelles.resize(1);
//                 CL[i].intervalles_temps.resize(1);
//                 CL[i].fcts_temporelles[0].resize(1);
//                 CL[i].intervalles_temps[0]=Vec<typename BOUNDARY::T,2>(0,100000);
//                 CL[i].fcts_temporelles[0]="1";
            } else {
                std::cout << "Mauvais type de calcul" << std::endl;
                assert(0);
            }
        }

    }
    
    
};

// modification des valeurs des boites de CL en fonction du parametre d'echelle
template<class BOUNDARY>
void modif_CL_scale(Vec<BOUNDARY > &CL, typename BOUNDARY::T &scale) {
    for(unsigned i=0;i<CL.size();++i)
        CL[i].box=CL[i].box*scale;
};
