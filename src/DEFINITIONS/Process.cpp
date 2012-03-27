#include "Process.h"
#include "SaveParameters.h"
#include "GeneralParameters.h"
#include "MultiScaleParameters.h"
#include "LatinParameters.h"
#include "TimeParameters.h"
#include "PROPERTY.h"
#include "MPIParameters.h"


Process::Process()
{
    // initialisation des valeurs
    nb_threads=1;
    sousint=1;
    type_sousint="p";
    nom_calcul="latin";
    recopie_t_post=0;
    save_data=0;
    allocation_memoire();
}

Process::~Process(){
    desallocation_memoire();
}


void Process::read_data_user(DataUser &data_user) {
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
    temps->read_data_user(data_user);
    
    if(data_user.options.Temp_statique == "statique"){
        if (rank==0) std::cout << "************************" << std::endl;
        if (rank==0) std::cout << "     STATIQUE           " << std::endl;
        if (rank==0) std::cout << "************************" << std::endl;
    }else if(data_user.options.Temp_statique == "quasistatique"){
        if (rank==0) std::cout << "************************" << std::endl;
        if (rank==0) std::cout << "     QUASISTATIQUE      " << std::endl;
        if (rank==0) std::cout << "************************" << std::endl;
    }
    nom_calcul = "incr";
};


void Process::allocation_memoire(){
    affichage  = new SaveParameters;
    structure  = new GeneralParameters;
    latin      = new LatinParameters;
    multiscale = new MultiScaleParameters;
    temps      = new TimeParameters;
    properties = new PROPERTY;
    multi_mpi  = new MPIParameters;
    #ifdef PRINT_ALLOC
    total_allocated[ typeid(SaveParameters).name() ] += sizeof(SaveParameters);
    total_allocated[ typeid(GeneralParameters).name() ] += sizeof(GeneralParameters);
    total_allocated[ typeid(LatinParameters).name() ]     += sizeof(LatinParameters);
    total_allocated[ typeid(MultiScaleParameters).name() ]     += sizeof(MultiScaleParameters);
    total_allocated[ typeid(TimeParameters).name() ]     += sizeof(TimeParameters);
    total_allocated[ typeid(PROPERTY).name() ]  += sizeof(PROPERTY);
    total_allocated[ typeid(MPIParameters).name() ] += sizeof(MPIParameters);
    #endif
}

void Process::desallocation_memoire(){
    if (affichage  != NULL) delete affichage;
    if (structure  != NULL) delete structure;
    if (latin      != NULL) delete latin;
    if (multiscale != NULL) delete multiscale;
    if (temps      != NULL) delete temps;
    if (properties != NULL) delete properties;
    if (multi_mpi  != NULL) delete multi_mpi;
}
