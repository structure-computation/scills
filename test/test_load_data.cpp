#ifdef METIL_COMP_DIRECTIVE
    #pragma src_file multiscale_geometry_mesh.cpp
    #pragma src_file assignation_material_properties_Sst.cpp
    #pragma src_file assignation_material_properties_Interface.cpp
    #pragma src_file multiscale_operateurs.cpp
    #if DIM == 2
        #pragma src_file formulation_2_double_elasticity_isotropy_stat_Qstat.cpp
        #pragma src_file formulation_2_double_elasticity_orthotropy_stat_Qstat.cpp
        #pragma src_file formulation_2_double_elasticity_damageable_isotropy_stat_Qstat.cpp
        #pragma src_file formulation_2_double_plasticity_isotropy_stat_Qstat.cpp
        #pragma src_file formulation_2_double_mesomodele.cpp
    #else
        #pragma src_file formulation_3_double_elasticity_isotropy_stat_Qstat.cpp
        #pragma src_file formulation_3_double_elasticity_orthotropy_stat_Qstat.cpp
        #pragma src_file formulation_3_double_elasticity_damageable_isotropy_stat_Qstat.cpp
        #pragma src_file formulation_3_double_plasticity_isotropy_stat_Qstat.cpp
        #pragma src_file formulation_3_double_mesomodele.cpp
    #endif
    #pragma src_file iterate_stat_Qstat.cpp
    #pragma src_file affichage.cpp
#endif



//lecture des parametres
#include "Process.h"

using namespace Metil;

#ifndef INFO_TIME
#define INFO_TIME
#endif 

#include "../LMT/include/io/ioexception.h"
using namespace LMT;
Crout crout;
#ifdef PRINT_ALLOC
namespace LMT {
    std::map<std::string,long long> total_allocated;
};
#endif


//*********************************************************************
// procedure multiechelle quasistatique ou statique (2 ou 3d)
//*********************************************************************
int main(int argc,char **argv) {
    /// lecture du fichier de donnees
    if ( argc!=3 and argc !=4 ) {
        std::cerr << "usage : nom_executable + dimension +  nom_complet_du_fichier_de_donnees. (+ mpi)" << std::endl;
        std::cerr << "ex : ./SC_multi_2 num_model num_calcul  (+ mpi)" << std::endl;
        return 1;
    }
    try {
        Sc2String id_model = argv[ 1 ];
        Sc2String id_calcul = argv[ 2 ];
        
        Process process;
        process.initialisation_MPI(argc, argv);
        process.lecture_fichiers(id_model, id_calcul);
        process.geometry_user->split_group_edges_within_geometry(*process.data_user);
        process.geometry_user->visualize_group_edges_within_geometry(*process.data_user);
        process.test_load_data();
        process.finalisation_MPI();
        if(process.parallelisation->is_master_cpu()) std::cout << "End of SC_multi_" << DIM << ".exe " << id_model << " " << id_calcul << std::endl;
        
    } catch ( const IoException &e ) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}
