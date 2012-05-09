#include "plasticity_functions.h"


/// Retourne l'ecrouissage associee a une plasticite cumulee 'p' sur la sous-structure 'S'
inline TYPEREEL fonction_ecrouissage(Sst &S, TYPEREEL p) {
    return S.matprop.plast_ecrouissage_mult * std::pow( p , S.matprop.plast_ecrouissage_expo ) + S.matprop.plast_ecrouissage_init;
}

/// Retourne la derivee de la fonction d'ecrouissage sur 'S' en 'p'
inline TYPEREEL derivee_fonction_ecrouissage(Sst &S, TYPEREEL p) {
    return S.matprop.plast_ecrouissage_mult * S.matprop.plast_ecrouissage_expo * std::pow( p , S.matprop.plast_ecrouissage_expo - 1 );
}


/// Retourne la contrainte equivalente du vecteur 'sigma' dans le solide 'S'
inline TYPEREEL contrainte_equivalente(const Sst &S,const Vec<TYPEREEL,DIM*(DIM+1)/2> &sigma){
    return std::sqrt(sigma[0]*sigma[0]-sigma[0]*sigma[1]
                    +sigma[1]*sigma[1]-sigma[1]*sigma[2]
                    +sigma[2]*sigma[2]-sigma[2]*sigma[0]
                    +3*sigma[3]*sigma[3]
                    +3*sigma[4]*sigma[4]
                    +3*sigma[5]*sigma[5]);
}

/// Retourne la derivee (tensorielle) de la fonction de contrainte equivalente sur 'S' en 'sigma'
inline Vec<TYPEREEL,DIM*(DIM+1)/2> &derivee_contrainte_equivalente(const Sst &S,Vec<TYPEREEL,DIM*(DIM+1)/2> &sigma){
    Vec<TYPEREEL,DIM*(DIM+1)/2> *df_dsigma = new Vec<TYPEREEL,DIM*(DIM+1)/2>;
    (*df_dsigma) = sigma;
    TYPEREEL spher_sigma = (sigma[0]+sigma[1]+sigma[2])/3;  /// Pression hydrostatique
    TYPEREEL sigma_eq = contrainte_equivalente(S,sigma);    /// Contrainte equivalente
    (*df_dsigma)[0] -= spher_sigma;
    (*df_dsigma)[1] -= spher_sigma;
    (*df_dsigma)[2] -= spher_sigma;    /// 'sigma' est devenu sa partie deviatorique
    (*df_dsigma) *= 1.5/sigma_eq;      /// 'sigma' est devenu df/dsigma
    return (*df_dsigma);
}


/// Calcul l'eventuelle plastification de la sous-structure 'S'
void calcul_plasticite(Sst &S, Process &process) {
    /// Reconstruction et stockage de sigma sur les Sst
    if (S.plastique) // or S.endommageable) ENDOMMAGEMENT NON IMPLEMENTE
    {
        Sst::Time &t_old = S.t[process.temps->pt-1];
        Sst::Time &t_cur = S.t[process.temps->pt];
        bool cinematique = (S.matprop.type_plast.find("cinematique") < S.matprop.type_plast.size());
        /// Reactualisation de la formulation
        S.f->set_mesh(S.mesh.m);
        upload_q(S,t_cur.q);
        upload_p(S,t_cur.p);
        upload_epsilon_p(S,t_old.epsilon_p);
        S.f->update_variables();
        S.f->call_after_solve();
        /// Recuperation des contraintes a corriger
        download_sigma(S,t_cur.sigma);
        /// Grandeurs utiles
        TYPEREEL mu = S.matprop.elastic_modulus/(2.*(1.+S.matprop.poisson_ratio));      /// Second coefficient de Lame
        if(S.mesh.type_formulation == "mesomodele") mu = (S.matprop.shear_modulus_12+S.matprop.shear_modulus_13+S.matprop.shear_modulus_23)/3;    // LE TEMPS D'IMPLEMENTER LA VERSION POUR LES MATERIAUX ORTHOTROPES
        /// Application de l'algorithme de retour radial
        for(unsigned i_elem = 0; i_elem < S.mesh.elem_list_size;i_elem++){
        
            /// Calcul du predicateur alpha_dev et de la fonction seuil f
            TYPEREEL p_old = t_old.p[i_elem];
            TYPEREEL R_p = t_cur.R_p[i_elem] = fonction_ecrouissage(S,p_old);
            Vec<TYPEREEL,DIM*(DIM+1)/2> alpha_elas = t_cur.sigma[i_elem];
            if(cinematique) {alpha_elas -= t_old.X_p[i_elem];}
            TYPEREEL alpha_eq = contrainte_equivalente(S,alpha_elas);
            TYPEREEL f = alpha_eq - R_p;
            
            if(f>0){ /// Correction necessaire
                
                /// Calcul de la variation de p (Methode de Newton)
                TYPEREEL dp = 0;//*
                TYPEREEL erreur = f/R_p;   /// On adimensionne par la limite d'elasticite actuelle
                while(erreur > 1e-6){
                    dp = dp + R_p*erreur/(3*mu+derivee_fonction_ecrouissage(S,p_old+dp)); /// Attention au signe
                    erreur = (alpha_eq - 3*mu*dp - fonction_ecrouissage(S,p_old+dp))/R_p;
                }
                t_cur.p[i_elem] = p_old + dp;
                
                /// Calcul de la variation de epsilon_p
                Vec<TYPEREEL,DIM*(DIM+1)/2> depsilon_p = dp*derivee_contrainte_equivalente(S,alpha_elas);
                t_cur.epsilon_p[i_elem] = t_old.epsilon_p[i_elem] + depsilon_p;
                
                /// Calcul de la variation de X_p
                if(cinematique) {t_cur.X_p[i_elem] = t_old.X_p[i_elem] + S.matprop.plast_cinematique_coef * depsilon_p;}
                
            } else { /// Pas de plastification
                t_cur.p[i_elem] = p_old;
                t_cur.epsilon_p[i_elem] = t_old.epsilon_p[i_elem];
                if(cinematique) {t_cur.X_p[i_elem] = t_old.X_p[i_elem];}
            }
            /*  AFFICHAGE
            //if(i_elem == 281) std::cout << "R_0      : " << S.matprop.plast_ecrouissage_init << endl;
            //if(i_elem == 281) std::cout << "k_p      : " << S.matprop.plast_ecrouissage_mult << endl;
            //if(i_elem == 281) std::cout << "m_p      : " << S.matprop.plast_ecrouissage_expo << endl;
            if(i_elem == 281) std::cout << "p_old    : " << p_old << std::endl << "epsilon_p_old : " << t_old.epsilon_p[i_elem] << std::endl;
            if(i_elem == 281) std::cout << "R_p      : " << R_p << std::endl;
            if(i_elem == 281) std::cout << "alpha_eq : " << alpha_eq << std::endl;
            if(i_elem == 281) std::cout << "f        : " << f << std::endl;
            if(i_elem == 281) std::cout << "p : " << t_cur.p[i_elem] << std::endl << "epsilon_p : " << t_cur.epsilon_p[i_elem] << std::endl;
            //*/
        }
    }
}