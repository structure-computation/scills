// Author: D. Violeau , 24/03/2006

//librairies Hugo
#include "containers/mat.h"
#include "containers/vecpointedvalues.h"
#include "mesh/mesh.h"
#include <fstream>
#include <map>
#include "containers/vec_mt.h"

//fichiers de definition des variables
#include "Boundary.h"
#include "Process.h"
#include "GeneralParameters.h"
#include "Sst.h"
#include "Interface.h"


// fonction utilisees pour le maillage et la geometrie
#include "calculate_measure_G_SST.h"
#include "create_mesure_G_neq_INTER.h"
#include "create_SST_INTER.h"
#include "p_surdiscretisation.h"
#include "surdiscretisation.h"
#include "p_create_SST_edge.h"

#include "assignation_mpi.h"


#include "DataUser.h"
#include "GeometryUser.h"

using namespace LMT;


void multiscale_geometry_mesh(DataUser                          &data_user, 
                              GeometryUser                      &geometry_user,
                              Vec<Sst>                          &S,
                              Vec<Interface>                    &Inter, 
                              Process                           &process, 
                              Vec<Boundary>                     &CL,
                              Vec<VecPointedValues<Sst> >       &Stot,
                              Vec<VecPointedValues<Sst> >       &SubS,
                              Vec<VecPointedValues<Interface> > &SubI) {
    
    /// creation SST et INTERFACE - maillage
    create_SST_INTER(data_user, geometry_user,S,Inter,CL,process,Stot,SubS,SubI);
    
    ///creation des normales, centre de gravite et mesure des interfaces
    if (process.parallelisation->is_master_cpu()) std::cout << " - Centre de gravite des interfaces + mesure + normales equivalentes" << std::endl;
    apply_mt(SubI,process.parallelisation->nb_threads,calculate_measure_G_neq_INTER(),S);
    
    if (process.sousint==1) {
        if (process.parallelisation->is_master_cpu()) std::cout << " - Modification des SST : Sur-discretisation des elements" << std::endl;
        if (process.type_sousint=="h") {
            /// surdiscretisation des ssts
            apply_mt(SubS,process.parallelisation->nb_threads,surdiscretise_SST());
            /// correspondance interface - cote SST - creation des nodeeq
            apply_mt(SubI,process.parallelisation->nb_threads,decoup_INTER(),S);
        } else if(process.type_sousint=="p") {
            //apply_mt(S,process.parallelisation->nb_threads,p_surdiscretise_SST());
            for( unsigned i=0;i<SubS.size() ;i++ ){
                //               S[i].mesh.sousint(Number< TV1::template SubType<0>::T::dim-1 >() );
                SubS[i].mesh.sousintegration = "p";
            }
            apply_mt(SubI,process.parallelisation->nb_threads,p_decoup_INTER(),S);
        } else {if (process.parallelisation->is_master_cpu()) std::cout << "Mauvais choix de sous integration : h ou p " << std::endl; assert(0);}
        
    } else {
        /// correspondance interface - cote SST - creation des nodeeq
        apply_mt(SubI,process.parallelisation->nb_threads,copy_INTER(),S);
    }
    
    if (process.parallelisation->is_master_cpu()) std::cout << " - Centre de gravite des SST et cotes + mesure SST" << std::endl;
    apply_mt(SubS,process.parallelisation->nb_threads,calculate_measure_G_SST());
    
    ///mettre à jour la taille des sst sur tous les pros
    if (process.parallelisation->is_multi_cpu() and process.sousint==1)
        for (unsigned i=0;i<S.size();i++){
            for( unsigned k1=0;k1 <process.parallelisation->repartition_sst.size()  ;k1++ )
                if (find(process.parallelisation->repartition_sst[k1],LMT::_1==(int)i)) {
                    MPI_Bcast(&S[i].mesh.node_list_size,1,MPI_INT,k1,MPI_COMM_WORLD);
                    MPI_Bcast(&S[i].mesh.elem_list_size,1,MPI_INT,k1,MPI_COMM_WORLD);
                    MPI_Bcast(S[i].G.ptr(),S[i].G.size(),MPI_DOUBLE,k1,MPI_COMM_WORLD);
                    break;
                }
        }
        
        if (process.parallelisation->is_master_cpu())
            calcul_taille_probleme(S, Inter);
}


