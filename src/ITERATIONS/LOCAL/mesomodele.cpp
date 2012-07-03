#include "mesomodele.h"
#include "../../DEFINITIONS/TimeData.h"
#include "../../MAILLAGE/elements_variables_accessor.h"
#include "../../../LMT/include/containers/mat.h"

using namespace LMT;


// Reactualisation des matrices de rigidite
void reactualisation_rigidite::operator()(Sst &S, VecInterfaces &Inter, const Process &process, const DataUser &data_user) const {
    if (S.update_operator) {
        /// Chargement des variables intervenant dans la definition des operateurs
        if(S.endommageable){ 
            /// Chargement de l'endommagement
            upload_d1(S,process.temps->pt-1);
            upload_d2(S,process.temps->pt-1);
            upload_df(S,process.temps->pt-1);
        }
        /// Creation de la matrice de rigidite
        S.mesh.load();
        S.assign_material_on_element(data_user);
        S.f->set_mesh(S.mesh.m);
        S.f->update_variables();
        S.f->clean_mats();
        S.f->free_matrices();
        S.f->allocate_matrices();
        S.f->assemble(true,false);
#if LDL
        SymetricMatrix *Kl;
        S.f->get_mat(Kl);
        SymetricMatrix &K = *Kl;
#else
        S.f->get_mat( S.K );
        CholModMatrix &K = *S.K;
#endif
        /// ajout des directions de recherche sur chaque cote
        for(unsigned j=0;j<S.edge.size();++j) {
            unsigned q=S.edge[j].internum;
            unsigned data=S.edge[j].datanum;            
            K[S.edge[j].repddledge] += Inter[q].side[data].Nt*(Inter[q].side[data].M*(Inter[q].side[data].kglo*Inter[q].side[data].N))/process.temps->dt; // a optimiser
        }
        
        S.update_operator = false;
    }
}

/* // L'ANELASTICITE N'EST PAS CALCULEE
/// Retourne l'ecrouissage associee a une plasticite cumulee 'p' sur la sous-structure 'S'
inline Scalar fonction_ecrouissage(Sst &S, Scalar p) {
    return S.matprop->plast_ecrouissage_mult * std::pow( p , S.matprop->plast_ecrouissage_expo ) + S.matprop->plast_ecrouissage_init;
}

/// Retourne la contrainte equivalente du vecteur 'sigma' dans le solide 'S'
inline Scalar contrainte_equivalente(const Sst &S,const Vec<Scalar,DIM*(DIM+1)/2> &sigma){
    return std::sqrt((sigma[1]*sigma[1]+sigma[2]*sigma[2])*S.matprop->couplage
                     +sigma[3]*sigma[3]
                     +sigma[4]*sigma[4]
                     +sigma[5]*sigma[5]);
}
//*/

void calcul_mesomodele(Sst &S, Process &process){
    /// Grandeurs utiles
    int pt = process.temps->pt;
    Sst::Time &t_old = S.t[pt-1];
    Sst::Time &t_new = S.t[pt];
    const Scalar &E1 = S.matprop->elastic_modulus_1;
    const Scalar &E2 = S.matprop->elastic_modulus_2;
    const Scalar &E3 = S.matprop->elastic_modulus_3;
    const Scalar &nu12 = S.matprop->poisson_ratio_12;
    const Scalar &nu13 = S.matprop->poisson_ratio_13;
    const Scalar &nu23 = S.matprop->poisson_ratio_23;
    const Scalar &G12 = S.matprop->shear_modulus_12;
    const Scalar &G13 = S.matprop->shear_modulus_13;
    const Scalar &G23 = S.matprop->shear_modulus_23;
    const Scalar &sqrt_Y0 = std::sqrt(S.matprop->Yo);
    const Scalar &sqrt_Yc = std::sqrt(S.matprop->Yc);
    const Scalar &Ycf = S.matprop->Ycf;
    const Scalar &b_c = S.matprop->b_c;
    const Scalar &dmax = S.matprop->dmax;
    const Scalar &a = S.matprop->a;
    const Scalar &tau_c = S.matprop->tau_c;
    const Scalar &dt = process.temps->dt;
    
    /// Reactualisation de la formulation
    S.f->set_mesh(S.mesh.m);
    upload_q(S,t_new);
    upload_d1(S,t_old);
    upload_d2(S,t_old);
    upload_df(S,t_old);
    S.f->update_variables();
    S.f->call_after_solve();
    /// Recuperation des contraintes a corriger
    Vec<Vec<Scalar,DIM*(DIM+1)/2> > all_sigma;
    all_sigma.resize(S.mesh.elem_list_size);
    download_sigma(S,all_sigma);
    
    /// CALCUL DE L'ENDOMMAGEMENT
    for(unsigned i_elem = 0; i_elem < S.mesh.elem_list_size; i_elem++){
        Vec<Scalar,DIM*(DIM+1)/2> &sigma = all_sigma[i_elem];
        Scalar d1_old = t_old.d1[i_elem];
        Scalar d2_old = t_old.d2[i_elem];
        Scalar df_old = t_old.df[i_elem];
        Scalar &d1_new = t_new.d1[i_elem];
        Scalar &d2_new = t_new.d2[i_elem];
        Scalar &df_new = t_new.df[i_elem];
        
        /// Calcul des forces associees a l'endommagement
        Vec<Scalar,3> &Yd_old = t_old.Yd[i_elem];
        Vec<Scalar,3> &Yd_new = t_new.Yd[i_elem];
        /// Calcul de Yd1
        Yd_new[0] = std::max((sigma[3]*sigma[3]/G12 + sigma[4]*sigma[4]/G23 + sigma[5]*sigma[5]/G13)/(1.0-d1_old)/(1.0-d1_old),Yd_new[0]);
        /// Calcul de Yd2
        Yd_new[1] = std::max((sigma[2]*sigma[2]/E2 + sigma[3]*sigma[3]/E3)/(1.0-d2_old)/(1.0-d2_old),Yd_new[1]);
        /// Calcul de Ydf
        Yd_new[2] = std::max((sigma[0]*sigma[0]/E1 - sigma[0]*sigma[1]*(2*nu12/E1) - sigma[0]*sigma[2]*(2*nu13/E1) - sigma[1]*sigma[2]*(2*nu23/E2))/(1.0-df_old)/(1.0-df_old),Yd_new[2]);
        
        Scalar wd = std::max((std::sqrt(Yd_new[0])-sqrt_Y0)/(sqrt_Yc-sqrt_Y0),0.0);
        if(S.matprop->effet_retard){
            /// Calcul de d1 avec effet retard
            Scalar expo = -a*(wd > d1_new)*(wd - d1_new);
            Scalar erreur = 1.0 ;
            while (erreur > 0.001) {
                erreur = d1_new;
                d1_new = d1_old + (dt/tau_c)*(1-std::exp(expo)) ;
                erreur = std::abs(erreur-d1_new);
            }
            d1_new = std::min(d1_new,dmax);
            d1_new = std::max(d1_new,d1_old);
            /// Calcul de d2 avec effet retard
            expo = -a*(b_c*wd > d2_new)*(b_c*wd - d2_new);
            erreur = 1.0 ;
            while (erreur > 0.001) {
                erreur = d2_new;
                d2_new = d2_old + (dt/tau_c)*(1-std::exp(-a*expo)) ;
                erreur = std::abs(erreur-d2_new);
            }
            d2_new = std::min(d2_new,dmax);
            d2_new = std::max(d2_new,d2_old);
        } else {
            /// Calcul de d1 sans effet retard
            d1_new = std::min(wd,dmax);
            d1_new = std::max(d1_new,d1_old);
            /// Calcul de d2 sans effet retard
            d2_new = std::min(b_c*wd,dmax);
            d2_new = std::max(d2_new,d2_old);
        }
        
        /// Calcul de df
        if(Yd_new[2] >= Ycf){
            df_new = dmax;
        } else {
            df_new = df_old;
        }
        
        /// Verification des changements
        for(unsigned i = 0; i < 3; i++){
            if(d1_new != d1_old or d2_new != d2_old or df_new != df_old){
                S.update_operator = true;
                /* DEBUG
                std::cout << "Modification de l'endommagement sur l'element " << i_elem << std::endl;
                std::cout << "    Yd  = " << Yd_new[0] << std::endl;
                std::cout << "    Yd' = " << Yd_new[1] << std::endl;
                std::cout << "    Ydf = " << Yd_new[2] << std::endl;
                std::cout << "    delta d  = " << (dt/tau_c)*(1-std::exp(-a*(wd > D_new[0])*(wd - D_new[0]))) << std::endl;
                std::cout << "    delta d' = " << (dt/tau_c)*(1-std::exp(-a*(b_c*wd > D_new[1])*(b_c*wd - D_new[1]))) << std::endl;
                std::cout << "    new d  = " << D_new[0] << std::endl;
                std::cout << "    new d' = " << D_new[1] << std::endl;
                std::cout << "    new df = " << D_new[2] << std::endl;
                */
            }
        }
        
        /* DEBUG
        if(i_elem == 1000){
            std::cout << "Modification de l'endommagement sur l'element " << i_elem << std::endl;
            std::cout << "    sigma  = " << sigma << std::endl;
            std::cout << "    old d  = " << D_old[0] << std::endl;
            std::cout << "    old d' = " << D_old[1] << std::endl;
            std::cout << "    old df = " << D_old[2] << std::endl;
            std::cout << "    Yd  = " << Yd_new[0] << std::endl;
            std::cout << "    Yd' = " << Yd_new[1] << std::endl;
            std::cout << "    Ydf = " << Yd_new[2] << std::endl;
            std::cout << "    wd = " << wd << std::endl;
            std::cout << "    new d  = " << D_new[0] << std::endl;
            std::cout << "    new d' = " << D_new[1] << std::endl;
            std::cout << "    new df = " << D_new[2] << std::endl;
        }
        //*/
    }
    
    return;
    /* // L'ANELASTICITE N'EST PAS CALCULEE
    /// CALCUL DE L'ANELASTICITE ASSOCIEE
    /// Application de l'algorithme de retour radial
    for(unsigned i_elem = 0; i_elem < S.mesh.elem_list_size;i_elem++){
        
        /// Calcul de la fonction seuil f
        Scalar p_old = t_old.p[i_elem];
        Scalar R_p = t_new.R_p[i_elem] = fonction_ecrouissage(S,p_old);
        Vec<Scalar,DIM*(DIM+1)/2> sig_new = t_new.sigma[i_elem];
        Scalar seq_new = contrainte_equivalente(S,sig_new);
        Scalar f = seq_new - R_p;
        
        if (f > 0) {
            /// Calcul de la variation de contrainte
            Vec<Scalar,DIM*(DIM+1)/2> sig_old;
            Scalar seq_old, seq_dot;
            if (process.temps->pt_cur>1) {
                sig_old = t_old.sigma[i_elem];
                seq_old = contrainte_equivalente(S,sig_old);
                seq_dot = (seq_new - seq_old)/process.temps->dt;
            }
            else {
                seq_dot = seq_new/process.temps->dt;
            }
            Scalar K = seq_dot/(S.matprop->plast_cinematique_coef*S.matprop->plast_ecrouissage_expo) ;
            
            /// Calcul de la variation de "plasticite"
            Scalar erreur = 1 ;
            Scalar p_new = p_old ;
            while (erreur > 0.0001) {
                erreur = p_new ;
                if (S.matprop->plast_ecrouissage_expo == 1.0) {
                    p_new = p_old + dt*K ;
                    erreur = 0 ;
                } else {
                    p_new = p_old + dt*K*std::pow(p_new,1-S.matprop->plast_ecrouissage_expo);
                    if ((p_new+erreur) != 0.0) erreur = (p_new-erreur)/(p_new+erreur);
                    else erreur=(p_new-erreur);
                }
            }
            /// Maj des valeurs
            t_new.p[i_elem] = p_new ;
            for (unsigned j=1;j<3;j++) {
                if (seq_new != 0) t_new.epsilon_p[i_elem][j] = p_new*sig_new[j]/seq_new*S.matprop->couplage ;
                else t_new.epsilon_p[i_elem][j] = 0 ;
            }
            for (unsigned j=3;j<6;j++) {
                if (seq_new != 0) t_new.epsilon_p[i_elem][j] = p_new*sig_new[j]/seq_new ;
                else t_new.epsilon_p[i_elem][j] = 0 ;
            }
        } else { /// Pas de plastification
            t_new.p[i_elem] = p_old;
            t_new.epsilon_p[i_elem] = t_old.epsilon_p[i_elem];
        }
        
        if(i_elem == 1000){
            std::cout << "    old p  = " << t_old.p[i_elem] << std::endl;
            std::cout << "    old Rp = " << t_old.R_p[i_elem] << std::endl;
            std::cout << "    old S = " << contrainte_equivalente(S,t_old.sigma[i_elem]) << std::endl;
            std::cout << "    new S = " << contrainte_equivalente(S,t_new.sigma[i_elem]) << std::endl;
            std::cout << "    f = " << f << std::endl;
            std::cout << "    new p  = " << t_new.p[i_elem] << std::endl;
            std::cout << "    new Rp = " << fonction_ecrouissage(S,t_new.p[i_elem]) << std::endl;
            std::cout << "    epsilon p = " << t_new.epsilon_p[i_elem] << std::endl;/// Reactualisation de la formulation
        }
    }
    
    S.f->set_mesh(S.mesh.m);
    upload_q(S,t_new.q);
    upload_p(S,t_old.p);
    upload_epsilon_p(S,t_new.epsilon_p);
    S.f->update_variables();
    S.f->call_after_solve();
    /// Recuperation des contraintes a corriger
    Vec<Vec<Scalar,DIM*(DIM+1)/2> > new_sigma;
    new_sigma.resize(S.mesh.elem_list_size);
    download_sigma(S,new_sigma);
    std::cout << "    new sigma : " << new_sigma[1000] << std::endl;
    std::cout << "    new seq : " << contrainte_equivalente(S,new_sigma[1000]) << std::endl;
    //*/
}
