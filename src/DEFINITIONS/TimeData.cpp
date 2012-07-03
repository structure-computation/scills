#include "TimeData.h"
#include "../UTILS/utils_2.h"

TimeData::TimeData():
parameters(),
t_cur("t")
{
    parameters.addParameter(&t_cur); /// Ajoute t_cur comme parametre (de controle) du groupe
    theta = 1;
    nbpastemps=1;
    type_de_calcul="stat";
    nb_step=1;
    pt=0;
    pt_cur=0;
    dt=1;
    step_old = -1;
    step_cur = 0;
}

TimeData::~TimeData(){
    for(unsigned i = 0; i < parameters.user_parameters.size(); i++){
        delete parameters.user_parameters[i];
    }
    time_step.free();
}


void TimeData::read_data_user(const DataUser &data_user){
    const DataUser::Json_time_steps &time_steps = data_user.time_steps;
    /// Lecture du type de schema temporel
    if(time_steps.time_scheme == "static"){
        type_de_calcul= "stat";
        nbpastemps=1;
        dt=0;
        nb_step = 1;
        time_step.resize(nb_step);
        nbpastemps=nb_step;
        time_step[0].dt=1;
        time_step[0].t_ini=0;
        time_step[0].t_fin=0;
        time_step[0].nb_time_step=1;
    }else if(time_steps.time_scheme == "quasistatic"){
        type_de_calcul= "Qstat";
        nb_step = time_steps.collection_vec.size();
        time_step.resize(nb_step);
        expressions.resize(time_steps.parameter_collection_vec.size());
        /// Initialisation des parametres
        for(unsigned i_param = 0; i_param < time_steps.parameter_collection_vec.size(); i_param++){
            std::cout << "*******************************************************************************************" << std::endl;
            std::cout << "time_steps.parameter_collection_vec[i_param].name : " << time_steps.parameter_collection_vec[i_param].name << std::endl;
            UserParameter *PT = new UserParameter(time_steps.parameter_collection_vec[i_param].name,&parameters);
            parameters.addParameter(PT);
            std::cout << "PT->symbol : " << PT->symbol << std::endl;
            std::cout << "PT->group  : " << PT->group << std::endl;
            std::cout << "*******************************************************************************************" << std::endl;
            expressions[i_param].resize(nb_step);
        }
        /// Recuperation des steps et des expressions de parametres associees
        nbpastemps = 0;
        for(unsigned i_step=0;i_step<nb_step;i_step++){
            time_step[i_step].id = time_steps.collection_vec[i_step].id_in_calcul;
            time_step[i_step].dt = time_steps.collection_vec[i_step].time_step;
            time_step[i_step].t_ini = time_steps.collection_vec[i_step].initial_time;
            time_step[i_step].t_fin = time_steps.collection_vec[i_step].final_time;
            time_step[i_step].nb_time_step = time_steps.collection_vec[i_step].nb_time_steps;
            nbpastemps += time_step[i_step].nb_time_step;
            for(unsigned i_param = 0; i_param < time_steps.parameter_collection_vec.size(); i_param++){
                expressions[i_param][i_step] = time_steps.parameter_collection_vec[i_param].stepFunctions_vec[i_step].temporal_function_t;
            }
        }
    }else{
        std::cerr << "Type de schema temporel inconnu : " << time_steps.time_scheme << std::endl;
        assert(0);
    }
    
    /*
    /// Lecture du type de schema temporel
    if(data_user.options.Temp_statique == "statique"){
        type_de_calcul= "stat";
    }else if(data_user.options.Temp_statique == "quasistatique"){
        type_de_calcul= "Qstat";
    }
    
    /// Lecture des steps
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
    }//*/
}

void TimeData::affiche(){
    std::cout << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "********************* Debug TimeData : **********************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    debug("int pt                   ",pt);
    debug("int pt_cur               ",pt_cur);
    debug("int step_cur             ",step_cur);
    debug("int step_old             ",step_old);
    debug("TYPEREEL t_cur           ",t_cur.value);
    debug("TYPEREEL dt              ",dt);
    debug("TYPEREEL theta           ",theta);
    debug("unsigned nbpastemps      ",nbpastemps);
    debug("Sc2String type_de_calcul ",type_de_calcul);
    debug("int nb_step              ",nb_step);
    std::cout << "expressions               : (" << expressions.size();
    if(expressions.size()>0) std::cout << "," << expressions[0].size() << ")";
    std::cout << std::endl;
    for(int i = 0; i < expressions.size(); i++){
        std::cout << "    ";
        for(int j = 0; j < expressions[i].size(); j++){
            std::cout << expressions[i][j] << " , ";
        }
        std::cout << std::endl;
    }
    std::cout << "parameters : " << std::endl; parameters.affiche();
    //debug("Vec<STEP> time_step",);
    std::cout << "*************************************************************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
}

void TimeData::init(){
    step_old = -1;
    step_cur = 0;
    pt = 1;
    pt_cur = 1;
    t_cur.value = time_step[0].t_ini;
    dt = time_step[0].dt;
}

void TimeData::next(){
    pt_cur ++;
    time_step[step_cur].pt_cur ++;
    step_old = step_cur; ///
    /// En cas de changement de step
    if(time_step[step_cur].pt_cur >= time_step[step_cur].nb_time_step){
        step_cur ++;
        time_step[step_cur].pt_cur = 0;
        dt = time_step[step_cur].dt;
    }
    t_cur.value = time_step[step_cur].t_ini + (time_step[step_cur].pt_cur+1) *dt;
}

bool TimeData::has_next(){
    return (pt_cur <= nbpastemps);
}

bool TimeData::step_changed(){
    return (step_cur == step_old);
}

void TimeData::prepareParameters(){
    /// Pas tres utile puisque 'step_changed()' devrait renvoyer 'false' (meme a la premiere iteration) dans updateParameters,
    /// mais comme ca les parametres multi-resolution auront egalement leurs expressions assignees et pretes (comme dans les autres classes)
    for(int i = 0; i < expressions.size(); i++){
        parameters.user_parameters[i]->setExpression(expressions[i][step_cur]);
    }
    parameters.prepareParameters();
}

void TimeData::updateParameters(){
    if(step_changed()){
        for(int i = 0; i < expressions.size(); i++){
            parameters.user_parameters[i]->setExpression(expressions[i][step_cur]);
        }
        parameters.prepareParameters();
    }
    parameters.updateParameters();
}
