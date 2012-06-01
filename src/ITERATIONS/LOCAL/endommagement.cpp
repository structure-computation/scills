#include "endommagement.h"
#include "../../MAILLAGE/elements_variables_accessor.h"
#include "../../../LMT/include/containers/mat.h"

using namespace LMT;

void calcul_endommagement(Sst &S, Process &process){
    /// Grandeurs utiles
    int pt = process.temps->pt;
    Sst::Time &t_old = S.t[pt-1];
    Sst::Time &t_new = S.t[pt];
    const TYPEREEL Yo = S.matprop.Yo;
    const TYPEREEL b_c = S.matprop.b_c;
    const TYPEREEL dmax = S.matprop.dmax;
    
    /// Reactualisation de la formulation
    S.f->set_mesh(S.mesh.m);
    upload_q(S,t_new);
    upload_d1(S,t_new);
    S.f->update_variables();
    S.f->call_after_solve();
    /// Recuperation des contraintes a corriger
    Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > all_sigma;
    all_sigma.resize(S.mesh.elem_list_size);
    download_sigma(S,all_sigma);
    Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > all_epsilon;
    all_epsilon.resize(S.mesh.elem_list_size);
    download_epsilon(S,all_epsilon);
    
    /// CALCUL DE L'ENDOMMAGEMENT
    for(unsigned i_elem = 0; i_elem < S.mesh.elem_list_size; i_elem++){
        Vec<TYPEREEL,DIM*(DIM+1)/2> &sigma = all_sigma[i_elem];
        Vec<TYPEREEL,DIM*(DIM+1)/2> &epsilon = all_epsilon[i_elem];
        TYPEREEL &d_old = t_old.d1[i_elem];
        
        /// Evaluation de la fonction seuil
        TYPEREEL energie,f;
#if DIM == 2
        energie = 0.5*(sigma[0]*epsilon[0]
                      +sigma[1]*epsilon[1]
                      +2*sigma[2]*epsilon[2]);  /// Energie de deformation en 2D
#else
        energie = 0.5*(sigma[0]*epsilon[0]
                      +sigma[1]*epsilon[1]
                      +sigma[2]*epsilon[2]
                      +2*sigma[3]*epsilon[3]
                      +2*sigma[4]*epsilon[4]
                      +2*sigma[5]*epsilon[5]);  /// Energie de deformation en 3D
#endif
        f = energie - Yo/(1-d_old/(1+b_c))/(1-d_old/(1+b_c));
        
        if(f>0){ /// Correction de l'endommagement
            t_new.d1[i_elem] = std::min((1+b_c)*(1+std::sqrt(Yo/energie)),dmax);  /// On impose f = 0 puis d <= dmax
            S.maj_operateurs = true;
        }else{ /// pas d'endommagement
            t_new.d1[i_elem] = d_old;
        }
    }
}
