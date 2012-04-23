#include "plasticity_functions.h"

/*
void calcul_sigma_sst(Sst &S, Process &process) {
    if (S.plastique) // or S.endommageable) ENDOMMAGEMENT NON IMPLEMENTE
    {
        unsigned elem;
        S.f->get_result() = S.t[process.temps->pt].q;
        elem=0 ;
        apply(S.mesh->elem_list,recuperation_deformation_plastique(),S.t[process.temps->pt].epsilon_p,elem) ;
        S.f->set_mesh(S.mesh.m);
        S.f->update_variables();
        S.f->call_after_solve();
        S.t[process.temps->pt].sigma.resize(0);
        elem = 0 ;
        apply(S.mesh->elem_list,stockage_sigma(),S.t[process.temps->pt].sigma,elem);
        elem = 0 ;
        apply(S.mesh->elem_list,stockage_sigma_vm(),S.t[process.temps->pt].sigma,elem);
    }
}
*/

TYPEREEL fonction_ecrouissage(Sst &S, TYPEREEL p) {
    return S.matprop.k_p * std::pow( p , S.matprop.m_p ) + S.matprop.R0;
}
TYPEREEL derivee_fonction_ecrouissage(Sst &S, TYPEREEL p) {
    return S.matprop.k_p * S.matprop.m_p * std::pow( p , S.matprop.m_p - 1 );
}


TYPEREEL calcul_contrainte_equivalente(const Sst &S,const Vec<TYPEREEL> &sigma){
    TYPEREEL seq = 0.0;
    TYPEREEL trace_sigma = 0.0;
    for(unsigned i = 0; i < DIM; i++)
        trace_sigma += sigma[i];
    for(unsigned i = 0; i < DIM; i++)  /// Boucle sur les termes diagonaux
        seq += (sigma[i] - trace_sigma/3.0)*(sigma[i] - trace_sigma/3.0);
    for(unsigned i = DIM; i < DIM*(DIM+1)/2; i++)    /// Boucle sur les termes hors diagonale
        seq += 2.0*sigma[i]*sigma[i];
    return std::sqrt(1.5*seq);
}

TYPEREEL trace(Vec<TYPEREEL,DIM*(DIM+1)/2> voigt){
    return voigt[0]+voigt[1]+voigt[2];
}

void deviatorique(const Vec<TYPEREEL,DIM*(DIM+1)/2> &voigt,Vec<TYPEREEL,DIM*(DIM+1)/2> &dev){
    TYPEREEL _trace = trace(voigt);
    dev[0] = voigt[0]-_trace/3;
    dev[1] = voigt[1]-_trace/3;
    dev[2] = voigt[2]-_trace/3;
    dev[3] = voigt[3];
    dev[4] = voigt[4];
    dev[5] = voigt[5];
}

TYPEREEL produit_scalaire_voigt(const Vec<TYPEREEL,DIM*(DIM+1)/2> &voigt1,const Vec<TYPEREEL,DIM*(DIM+1)/2> &voigt2){
    TYPEREEL res = 0;
    for(int i = 0; i < 3; i++)
        res += voigt1[i]*voigt2[i];
    for(int i = 3; i < 6; i++)
        res += 2*voigt1[i]*voigt2[i];
    return res;
}

TYPEREEL von_mises(const Vec<TYPEREEL,DIM*(DIM+1)/2> &sigma){
    Vec<TYPEREEL,DIM*(DIM+1)/2> dev;
    deviatorique(sigma,dev);
    return std::sqrt(1.5*produit_scalaire_voigt(dev,dev));
}

TYPEREEL R_max = 0;

void calcul_plasticite(Sst &S, Process &process) {
    /// Reconstruction et stockage de sigma sur les Sst
    if (S.plastique) // or S.endommageable) ENDOMMAGEMENT NON IMPLEMENTE
    {
        Sst::Time &t_old = S.t[process.temps->pt-1];
        Sst::Time &t_cur = S.t[process.temps->pt];
        /// Chargememt de q (dep) sur les elements
        S.f->set_mesh(S.mesh.m);
        S.f->get_result() = t_cur.q;
        /// Chargement de epsilon_p sur les elements
        unsigned i_elem = 0;
        apply(S.mesh->elem_list,chargement_deformation_plastique(),t_old.epsilon_p,i_elem);
        /// Reactualisation de la formulation
        S.f->update_variables();
        S.f->call_after_solve();
        /// Recuperation de sigma
        i_elem = 0;
        apply(S.mesh->elem_list,stockage_sigma(),t_cur.sigma,i_elem);
        /// Recuperation de sigma von mises
        i_elem = 0;
        apply(S.mesh->elem_list,compare_vm(),S,t_cur.sigma,i_elem);
        for(unsigned i_elem = 0; i_elem < S.mesh->elem_list.size();i_elem++){
            /// Calcul et stockage de l'ecrouissage R sur les Sst
            t_cur.R[i_elem] = fonction_ecrouissage(S,t_old.p[i_elem]);
        
            /// Calcul de p et epsilon_p
            TYPEREEL sigma_eq = calcul_contrainte_equivalente(S,t_cur.sigma[i_elem]);
            TYPEREEL f = sigma_eq - t_cur.R[i_elem];
            
            if(f>0){
                /// Grandeurs utiles
                TYPEREEL p_old = t_old.p[i_elem];
                TYPEREEL R_p = t_cur.R[i_elem];
                TYPEREEL mu = S.matprop.elastic_modulus_1/(2.*(1.+S.matprop.poisson_ratio_12));      /// Second coefficient de Lame
                const Vec<TYPEREEL,DIM*(DIM+1)/2> &sigma_elas = t_cur.sigma[i_elem];
                
                /// Calcul de la variation de p (Methode de Newton)
                TYPEREEL dp = 0;//*
                TYPEREEL erreur = f/R_p;   /// On adimensionne par la limite d'elasticite actuelle
                while(erreur > 1e-6){
                    dp = dp + R_p*erreur/(3*mu+derivee_fonction_ecrouissage(S,p_old+dp)); /// Attention au signe
                    erreur = (sigma_eq - 3*mu*dp - fonction_ecrouissage(S,p_old+dp))/R_p;
                }//*/
                //dp = f/(3*mu+derivee_fonction_ecrouissage(S,p_old));
                t_cur.p[i_elem] = p_old + dp;
                
                /// Calcul de la variation de epsilon_p
                Vec<TYPEREEL,DIM*(DIM+1)/2> sigma_dev;  /// Partie deviatorique de sigma_elas
                deviatorique(sigma_elas,sigma_dev);
                Vec<TYPEREEL,DIM*(DIM+1)/2> depsilon_p = 1.5*dp*sigma_dev/sigma_eq;
                t_cur.epsilon_p[i_elem] = t_old.epsilon_p[i_elem] + depsilon_p;
                
                //*   TEST
                Vec<TYPEREEL,DIM*(DIM+1)/2> sigma_new = sigma_elas - 2*mu*depsilon_p;
                R_max = (R_max>fonction_ecrouissage(S,t_cur.p[i_elem]))? R_max : fonction_ecrouissage(S,t_cur.p[i_elem]);
                if(S.id == 0 and i_elem == 281){
                    //std::cout << "Numero de l'element : " << i_elem << std::endl;
                    std::cout << "f old         : " << f << std::endl;
                    std::cout << "dp            : " << dp << std::endl;
                    std::cout << "p_old         : " << t_old.p[i_elem] << std::endl;
                    std::cout << "p_new         : " << t_cur.p[i_elem] << std::endl;
                    std::cout << "old R(p)      : " << R_p << std::endl;
                    std::cout << "new R(p)      : " << fonction_ecrouissage(S,t_cur.p[i_elem]) << std::endl;
                    std::cout << "f new         : " << calcul_contrainte_equivalente(S,sigma_new)-fonction_ecrouissage(S,t_cur.p[i_elem]) << std::endl;
                    //std::cout << "sigma dev : " << sigma_dev << std::endl;
                    //std::cout << "dEpsPn    : " << dEpsPn << std::endl;
                    //std::cout << "epsilon_p : " << elem.epsilon_p << std::endl;
                    std::cout << "epsilon_p old : " << t_old.epsilon_p[i_elem] << std::endl;
                    std::cout << "epsilon_p new : " << t_cur.epsilon_p[i_elem] << std::endl;
                    std::cout << "sigma old : " << sigma_elas << std::endl;
                    std::cout << "sigma new : " << sigma_new << std::endl;
                    std::cout << "Test0 : " << std::endl 
                              << "    trace eps_p:" << trace(depsilon_p) << std::endl
                              << "    R_new      :" << fonction_ecrouissage(S,t_cur.p[i_elem]) << std::endl 
                              << "    sigma_eq_1 :" << sigma_eq - 3*mu*dp << std::endl 
                              << "    sigma_eq_2 :" << calcul_contrainte_equivalente(S,sigma_new) << std::endl
                              << "    sigma_eq_3 :" << von_mises(sigma_new) << std::endl;
                    //assert(0);
                }
                //*/
            }
        }
    }
    std::cout << "************************************************************************************ R_max : " << R_max << std::endl;
}