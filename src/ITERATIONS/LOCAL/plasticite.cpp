#include "plasticite.h"
#include "../../DEFINITIONS/TimeData.h"


/// Retourne l'ecrouissage associee a une plasticite cumulee 'p' sur la sous-structure 'S'
inline Scalar fonction_ecrouissage(Sst &S, Scalar p) {
    return S.matprop->plast_ecrouissage_mult * std::pow( p , S.matprop->plast_ecrouissage_expo ) + S.matprop->plast_ecrouissage_init;
}

/// Retourne la derivee de la fonction d'ecrouissage sur 'S' en 'p'
inline Scalar derivee_fonction_ecrouissage(Sst &S, Scalar p) {
    return S.matprop->plast_ecrouissage_mult * S.matprop->plast_ecrouissage_expo * std::pow( p , S.matprop->plast_ecrouissage_expo - 1 );
}


/// Retourne la contrainte equivalente du vecteur 'sigma' dans le solide 'S'
inline Scalar contrainte_equivalente(const Sst &S,const VoigtVector &sigma){
    return std::sqrt(sigma[0]*sigma[0]-sigma[0]*sigma[1]
                    +sigma[1]*sigma[1]-sigma[1]*sigma[2]
                    +sigma[2]*sigma[2]-sigma[2]*sigma[0]
                    +3*sigma[3]*sigma[3]
                    +3*sigma[4]*sigma[4]
                    +3*sigma[5]*sigma[5]);
}

/// Retourne la derivee (tensorielle) de la fonction de contrainte equivalente sur 'S' en 'sigma'
inline VoigtVector &derivee_contrainte_equivalente(const Sst &S,VoigtVector &sigma){
    VoigtVector *df_dsigma = new VoigtVector;
    (*df_dsigma) = sigma;
    Scalar spher_sigma = (sigma[0]+sigma[1]+sigma[2])/3;  /// Pression hydrostatique
    Scalar sigma_eq = contrainte_equivalente(S,sigma);    /// Contrainte equivalente
    (*df_dsigma)[0] -= spher_sigma;
    (*df_dsigma)[1] -= spher_sigma;
    (*df_dsigma)[2] -= spher_sigma;    /// 'sigma' est devenu sa partie deviatorique
    (*df_dsigma) *= 1.5/sigma_eq;      /// 'sigma' est devenu df/dsigma
    return (*df_dsigma);
}


/// Calcul l'eventuelle plastification de la sous-structure 'S'
void calcul_plasticite(Sst &S, Process &process) {
    /// Reconstruction et stockage de sigma sur les Sst
    Sst::Time &t_old = S.t[process.temps->pt-1];
    Sst::Time &t_new = S.t[process.temps->pt];
    const Scalar mu = S.matprop->elastic_modulus/(2.*(1.+S.matprop->poisson_ratio));    /// Second coefficient de Lame
    bool cinematique = (S.matprop->type_plast.find("cinematique") < S.matprop->type_plast.size());
    /// Reactualisation de la formulation
    S.f->set_mesh(S.mesh.m);
    upload_q(S,t_new);
    upload_p(S,t_new);
    upload_epsilon_p(S,t_old);
    S.f->update_variables();
    S.f->call_after_solve();
    /// Recuperation des contraintes a corriger
    Vec<VoigtVector> all_sigma;
    all_sigma.resize(S.mesh.elem_list_size);
    download_sigma(S,all_sigma);
    /// Application de l'algorithme de retour radial
    for(unsigned i_elem = 0; i_elem < S.mesh.elem_list_size;i_elem++){
    
        /// Calcul du predicateur alpha_dev et de la fonction seuil f
        Scalar p_old = t_old.p[i_elem];
        Scalar R_p = t_new.R_p[i_elem] = fonction_ecrouissage(S,p_old);
        VoigtVector alpha_elas = all_sigma[i_elem];
        if(cinematique) {alpha_elas -= t_old.X_p[i_elem];}
        Scalar alpha_eq = contrainte_equivalente(S,alpha_elas);
        Scalar f = alpha_eq - R_p;
        
        if(f>0){ /// Correction necessaire
            
            /// Calcul de la variation de p (Methode de Newton)
            Scalar dp = 0;//*
            Scalar erreur = f/R_p;   /// On adimensionne par la limite d'elasticite actuelle
            while(erreur > 1e-6){
                dp = dp + R_p*erreur/(3*mu+derivee_fonction_ecrouissage(S,p_old+dp)); /// Attention au signe
                erreur = (alpha_eq - 3*mu*dp - fonction_ecrouissage(S,p_old+dp))/R_p;
            }
            t_new.p[i_elem] = p_old + dp;
            
            /// Calcul de la variation de epsilon_p
            VoigtVector depsilon_p = dp*derivee_contrainte_equivalente(S,alpha_elas);
            t_new.epsilon_p[i_elem] = t_old.epsilon_p[i_elem] + depsilon_p;
            
            /// Calcul de la variation de X_p
            if(cinematique) {t_new.X_p[i_elem] = t_old.X_p[i_elem] + S.matprop->plast_cinematique_coef * depsilon_p;}
            
        } else { /// Pas de plastification
            t_new.p[i_elem] = p_old;
            t_new.epsilon_p[i_elem] = t_old.epsilon_p[i_elem];
            if(cinematique) {t_new.X_p[i_elem] = t_old.X_p[i_elem];}
        }
        /*  DEBUG
        //if(i_elem == 281) std::cout << "R_0      : " << S.matprop->plast_ecrouissage_init << endl;
        //if(i_elem == 281) std::cout << "k_p      : " << S.matprop->plast_ecrouissage_mult << endl;
        //if(i_elem == 281) std::cout << "m_p      : " << S.matprop->plast_ecrouissage_expo << endl;
        if(i_elem == 281) std::cout << "p_old    : " << p_old << std::endl << "epsilon_p_old : " << t_old.epsilon_p[i_elem] << std::endl;
        if(i_elem == 281) std::cout << "R_p      : " << R_p << std::endl;
        if(i_elem == 281) std::cout << "alpha_eq : " << alpha_eq << std::endl;
        if(i_elem == 281) std::cout << "f        : " << f << std::endl;
        if(i_elem == 281) std::cout << "p : " << t_new.p[i_elem] << std::endl << "epsilon_p : " << t_new.epsilon_p[i_elem] << std::endl;
        //*/
    }
}