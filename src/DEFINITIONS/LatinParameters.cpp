#include "LatinParameters.h"

LatinParameters::LatinParameters(){
    iter=0;
    nbitermax=100;
    mu=0.8;
    ktype="scalaire_auto_CL";
    kfact=200e3;
    copydirection=1;
    critere_erreur=1e-4;
    facteur_relaxation=mu;
    type_error="ddr";
    list_error=0;
    save_depl_SST=1;
    alloc_quantites=1;
}

void LatinParameters::read_data_user(DataUser &data_user){
    nbitermax = data_user.options.LATIN_nb_iter_max;
    facteur_relaxation = 0.8;
    critere_erreur = data_user.options.LATIN_crit_error;
    critere_erreur_diss = 0;
    critere_erreur_auto_stop = data_user.options.LATIN_crit_error_auto_stop;
    type_error = "ddr";
    if (type_error=="dissipation")
        critere_erreur_diss = 0;
    ktype = "scalaire_auto_CL";
    kfact = 1;
    copydirection = 0; 
    list_error= 1;
}
