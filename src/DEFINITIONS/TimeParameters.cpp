#include "TimeParameters.h"
#include "../UTILS/utils_2.h"

TimeParameters::TimeParameters(){
    theta = 1;
    nbpastemps=1;
    type_de_calcul="stat";
    nb_step=1;
    pt=0;
    pt_cur=0;
    dt=1;
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
    debug("TYPEREEL current_time    ",current_time);
    debug("TYPEREEL dt              ",dt);
    debug("TYPEREEL theta           ",theta);
    debug("unsigned nbpastemps      ",nbpastemps);
    debug("Sc2String type_de_calcul ",type_de_calcul);
    debug("int nb_step              ",nb_step);
    //debug("Vec<STEP> time_step",);
    std::cout << "*************************************************************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
}

void TimeParameters::init(){
    step_cur = 0;
    pt = 0;
    pt_cur = 0;
    current_time = time_step[0].t_ini;
    dt = time_step[0].dt;
}

void TimeParameters::next(){
    pt_cur ++;
    time_step[step_cur].pt_cur ++;
    /// En cas de changement de step
    if(time_step[step_cur].pt_cur >= time_step[step_cur].nb_time_step){
        step_cur ++;
        time_step[step_cur].pt_cur = 0;
        dt = time_step[step_cur].dt;
        __step_changed = true;
    }
    current_time = time_step[step_cur].t_ini + (time_step[step_cur].pt_cur+1) *dt;
}

bool TimeParameters::has_next(){
    return (pt_cur <= nbpastemps);
}

bool TimeParameters::step_changed(){
    bool b = __step_changed;
    __step_changed = false;
    return b;
}
