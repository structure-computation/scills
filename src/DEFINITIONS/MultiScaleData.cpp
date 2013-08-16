#include "MultiScaleData.h"
#include "../UTILS/utils_2.h"


MultiScaleData::MultiScaleData(){
    multiechelle=1;
    type_base_macro=3;
    sizeM=0;
    erreur_macro=1e-6;
    opti_multi=1;
}

void MultiScaleData::read_data_user(const DataUser &data_user) {
    //multiechelle = (data_user.mesh.nb_groups_elem > 2); /// Uniquement si 3 pieces ou plus
    multiechelle = data_user.options.convergence_method_LATIN.multiscale;
    type_base_macro = 3;
    opti_multi = 0;
    erreur_macro = 0.000001;
}

void MultiScaleData::display_all_data(){
    std::cout << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "******************* Debug MultiScaleData : ******************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    debug("int multiechelle         ",multiechelle);
    debug("unsigned type_base_macro ",type_base_macro);
    debug("unsigned sizeM           ",sizeM);
    debug("Vec<unsigned> nb_fcts_macro_par_inter",nb_fcts_macro_par_inter);
    debug("Vec<unsigned> inter_correspondantes ",inter_correspondantes);
    debug("bool nbmacro_identique ",nbmacro_identique);
    debug("bool opti_multi        ",opti_multi);
    debug("TYPEREEL erreur_macro  ",erreur_macro);
    std::cout << "*************************************************************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
}