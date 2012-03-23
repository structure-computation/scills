#include "Param.h"
#include "definition_PARAM_AFFICHAGE.h"
#include "definition_PARAM_STRUCTURE.h"
#include "definition_PARAM_MULTI.h"
#include "definition_PARAM_LATIN.h"
#include "definition_PARAM_TEMPS.h"
#include "definition_PARAM_PROPERTY.h"
#include "definition_PARAM_MPI.h"


Param::Param()
{
    // initialisation des valeurs
    nb_threads=1;
    sousint=1;
    type_sousint="p";
    dim=2;
    nom_calcul="latin";
    recopie_t_post=0;
    save_data=0;
    allocation_memoire();
}

Param::~Param(){
    desallocation_memoire();
}


void Param::read_data_user(DataUser &data_user) {
    sousint = false;
    type_sousint = "h";
    rbm.bloq = false;
    rbm.mvts_bloques.resize(3);
    rbm.mvts_bloques[0]= "Ty";
    rbm.mvts_bloques[1]= "Tx";
    rbm.mvts_bloques[2]= "Rz";
    nb_threads = 1;
    save_data = true;
    read_data = false;
    reprise_calcul = 0;
    properties->deltaT = 0;
    
    structure->read_data_user(data_user);
    multiscale->read_data_user(data_user);
    latin->read_data_user(data_user);
    affichage->read_data_user(data_user);
    
    //     std::cout<< "data_user.options.Temp_statique = " << data_user.options.Temp_statique << std::endl;
    
    if(data_user.options.Temp_statique == "statique"){
        temps->type_de_calcul= "stat";
    }else if(data_user.options.Temp_statique == "quasistatique"){
        temps->type_de_calcul= "Qstat";
    }
    
    if (temps->type_de_calcul=="stat") {
        if (rank==0) std::cout << "************************" << std::endl;
        if (rank==0) std::cout << "     STATIQUE     " << std::endl;
        if (rank==0) std::cout << "************************" << std::endl;
        temps->nbpastemps=1;
        temps->dt=0;
        nom_calcul="incr";
        temps->nb_step = 1;
        temps->time_step.resize(temps->nb_step);
        temps->nbpastemps=temps->nb_step;
        for(unsigned i_step=0;i_step<temps->nb_step;i_step++){
            temps->time_step[i_step].dt=1;
            temps->time_step[i_step].t_ini=0;
            temps->time_step[i_step].t_fin=0;
            temps->time_step[i_step].nb_time_step=1;
        }
        
    }else if(temps->type_de_calcul=="Qstat") {
        if (rank==0) std::cout << "************************" << std::endl;
        if (rank==0) std::cout << "     QUASISTATIQUE      " << std::endl;
        if (rank==0) std::cout << "************************" << std::endl;
        nom_calcul="incr";
        temps->nb_step = data_user.time_step.size();
        temps->time_step.resize(temps->nb_step);
        temps->nbpastemps= 0;
        for(unsigned i_step=0;i_step<temps->nb_step;i_step++){
            temps->time_step[i_step].dt = data_user.time_step[i_step].dt;
            temps->time_step[i_step].t_ini = data_user.time_step[i_step].ti;
            temps->time_step[i_step].t_fin = data_user.time_step[i_step].tf;
            temps->time_step[i_step].nb_time_step = data_user.time_step[i_step].nb_time_step;
            temps->nbpastemps += temps->time_step[i_step].nb_time_step;
        }
    }
};


void Param::allocation_memoire(){
    affichage  = new AFFICHAGE;
    structure  = new STRUCTURE;
    latin      = new LATIN;
    multiscale = new MULTI;
    temps      = new TEMPS;
    properties = new PROPERTY;
    multi_mpi  = new MULTI_MPI;
    #ifdef PRINT_ALLOC
    total_allocated[ typeid(AFFICHAGE).name() ] += sizeof(AFFICHAGE);
    total_allocated[ typeid(STRUCTURE).name() ] += sizeof(STRUCTURE);
    total_allocated[ typeid(LATIN).name() ]     += sizeof(LATIN);
    total_allocated[ typeid(MULTI).name() ]     += sizeof(MULTI);
    total_allocated[ typeid(TEMPS).name() ]     += sizeof(TEMPS);
    total_allocated[ typeid(PROPERTY).name() ]  += sizeof(PROPERTY);
    total_allocated[ typeid(MULTI_MPI).name() ] += sizeof(MULTI_MPI);
    #endif
}

void Param::desallocation_memoire(){
    if (affichage  != NULL) delete affichage;
    if (structure  != NULL) delete structure;
    if (latin      != NULL) delete latin;
    if (multiscale != NULL) delete multiscale;
    if (temps      != NULL) delete temps;
    if (properties != NULL) delete properties;
    if (multi_mpi  != NULL) delete multi_mpi;
}
