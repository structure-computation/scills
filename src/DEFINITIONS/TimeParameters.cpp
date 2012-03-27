#include "TimeParameters.h"
#include "../UTILS/utils_2.h"

TimeParameters::TimeParameters(){
    pt=0;
    pt_cur=0;
    dt=1;
    theta = 1;
    nbpastemps=1;
    type_de_calcul="stat";
    nb_step=1;
    current_time=0;
}


void TimeParameters::read_data_user(DataUser &data_user){
    if(data_user.options.Temp_statique == "statique"){
        type_de_calcul= "stat";
    }else if(data_user.options.Temp_statique == "quasistatique"){
        type_de_calcul= "Qstat";
    }
    
    if (type_de_calcul=="stat") {
        nbpastemps=1;
        dt=0;
        nb_step = 1;
        time_step.resize(nb_step);
        nbpastemps=nb_step;
        for(unsigned i_step=0;i_step<nb_step;i_step++){
            time_step[i_step].dt=1;
            time_step[i_step].t_ini=0;
            time_step[i_step].t_fin=0;
            time_step[i_step].nb_time_step=1;
        }
        
    }else if(type_de_calcul=="Qstat") {
        nb_step = data_user.time_step.size();
        time_step.resize(nb_step);
        nbpastemps= 0;
        for(unsigned i_step=0;i_step<nb_step;i_step++){
            time_step[i_step].dt = data_user.time_step[i_step].dt;
            time_step[i_step].t_ini = data_user.time_step[i_step].ti;
            time_step[i_step].t_fin = data_user.time_step[i_step].tf;
            time_step[i_step].nb_time_step = data_user.time_step[i_step].nb_time_step;
            nbpastemps += time_step[i_step].nb_time_step;
        }
    }
}

void TimeParameters::display_all_data(){
    std::cout << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "******************* Debug TimeParameters : ******************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    debug("int pt                   ",pt);
    debug("int pt_cur               ",pt_cur);
    debug("int step_cur             ",step_cur);
    debug("double current_time      ",current_time);
    debug("double dt                ",dt);
    debug("double theta             ",theta);
    debug("unsigned nbpastemps      ",nbpastemps);
    debug("Sc2String type_de_calcul ",type_de_calcul);
    debug("int nb_step              ",nb_step);
    //debug("Vec<STEP> time_step",);
    std::cout << "*************************************************************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
}
