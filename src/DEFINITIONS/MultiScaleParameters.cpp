#include "MultiScaleParameters.h"
#include "../UTILS/utils_2.h"


MultiScaleParameters::MultiScaleParameters(){
    multiechelle=1;
    type_base_macro=3;
    sizeM=0;
    erreur_macro=1e-8;
    opti_multi=1;
}

void MultiScaleParameters::read_data_user(DataUser &data_user) {
    multiechelle = data_user.options.multiechelle;
    type_base_macro = 3;
    opti_multi = 0;
    erreur_macro = 0.000001;
}

void MultiScaleParameters::display_all_data(){
    std::cout << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "**************** Debug MultiScaleParameters : ***************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    debug("int multiechelle         ",multiechelle);
    debug("unsigned type_base_macro ",type_base_macro);
    debug("unsigned sizeM           ",sizeM);
    debug("Vec<unsigned> nb_fcts_macro_par_inter",nb_fcts_macro_par_inter);
    debug("Vec<unsigned> inter_correspondantes ",inter_correspondantes);
    debug("bool nbmacro_identique ",nbmacro_identique);
    debug("bool opti_multi        ",opti_multi);
    debug("double erreur_macro    ",erreur_macro);
    std::cout << "*************************************************************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
}