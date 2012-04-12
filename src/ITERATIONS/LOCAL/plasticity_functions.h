#ifndef PLASTICITY_FUNCTIONS_H
#define PLASTICITY_FUNCTIONS_H


#include "../../DEFINITIONS/Process.h"
#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/Interface.h"
#include "../../UTILS/utils_2.h"
//#include <boost/bind/bind_template.hpp>


TYPEREEL f_max = 0;
TYPEREEL p_max = 0;


/** Recuperation de la deformation plastique du "stock" vers la formulation
 * 
 */
struct recuperation_deformation_plastique {
    template<class TE, class ME>
    void operator()(TE &elem, ME &ep, unsigned &i_elem) const {
        elem.epsilon_p[0] = ep[i_elem] ;
        ++i_elem;    ///Sert d'index d'iteration sur les elements
    }
};


/**Stockage des deformations plastiques calculees
 * Les deformations plastiques calculees sont stockees dans le vecteur Time::epsilon_p du pas de temps courant
 */
struct stockage_deformation_plastique {
    template<class TE,class ME>
    void operator()(TE &elem, ME &ep) const {
        ep.push_back(elem.epsilon_p[0]);
    }
};


/**Stockage de la plasticite cumulee calculee
 * La plasticite cumulee calculee est stockee dans le vecteur Time::p du pas de temps courant
 */
struct stockage_plasticite_cumulee {
    template<class TE,class ME>
    void operator()(TE &elem, ME &p) const {
        p.push_back(elem.p);
    }
};


/**Stockage des contraintes calculees
 * Les contraintes calculee sont stockees dans le vecteur Time::sigma du pas de temps courant
 */
struct stockage_sigma {
    template<class TE,class ME>
    void operator()(TE &elem, ME &s, unsigned &i_elem) const {
        s.push_back(elem.sigma[0]);
        ++i_elem;    ///Sert d'index d'iteration
    }
};


void calcul_sigma_sst(Sst &S, Process &process) {
    if (S.plastique) // or S.endommageable) ENDOMMAGEMENT NON IMPLEMENTE
    {
        unsigned elem;
        S.f->get_result() = S.t[process.temps->pt_cur].q;
        elem=0 ;
        apply(S.mesh->elem_list,recuperation_deformation_plastique(),S.t[process.temps->pt_cur].epsilon_p,elem) ;
        S.f->set_mesh(S.mesh.m);
        S.f->update_variables();
        S.f->call_after_solve();
        S.t[process.temps->pt_cur].sigma.resize(0);
        elem = 0 ;
        apply(S.mesh->elem_list,stockage_sigma(),S.t[process.temps->pt_cur].sigma,elem);
    }
}


void calcul_ecrouissage_sst(Sst &S, Process &process) {
    if (S.plastique) {
        int elem = 0 ;
        for (unsigned i=0;i<S.mesh->elem_list.size();++i) {
            S.t[process.temps->pt_cur].R[i] = S.matprop.k_p * std::pow( S.t[process.temps->pt_cur].p[i] , S.matprop.m_p );
        }
    }
}


/** Calcul de la contrainte equivalente
 * Retourne la contrainte equivalente (type Von Mises) sur un element de S soumis a sigma
 */
TYPEREEL calcul_contrainte_equivalente(const Sst &S,const Vec<TYPEREEL> &sigma){
    TYPEREEL seq = 0;
    const Vec<Vec<TYPEREEL,DIM*(DIM+1)/2>,DIM*(DIM+1)/2> &coeff_seq = S.matprop.coeff_seq;
    for(int i = 0; i < DIM*(DIM+1)/2; i++){
        for(int j = 0; j < DIM*(DIM+1)/2; j++){
            seq+=coeff_seq[i][j]*sigma[i]*sigma[j];
        }
    }
    return std::sqrt(seq);
}


/** Recherche de la plasticite cumulee
 * Le premier argument servira de stockage pour le resultat
 */
void recherche_plasticite_cumulee(TYPEREEL &p, TYPEREEL K, TYPEREEL dt, TYPEREEL pold, TYPEREEL m) {
    TYPEREEL erreur = 1 ;
    TYPEREEL pnew = p ;
    while (erreur > 0.0001) {
        erreur = pnew ;
        if (m == 1.0) {
            pnew = pold + dt*K ;
            erreur = 0 ;
        }
        else {
            pnew = pold + dt*K*std::pow(pnew,1-m);
            if ((pnew+erreur) != 0.0) 
                erreur = (pnew-erreur)/(pnew+erreur);
            else 
                erreur=(pnew-erreur);
        }
    }
    p = pnew ;
}



/** Calcul de la plasticitee (non couplee a l'endommagement)
 * Calcul de la plasticite cumulee (p) et de la deformation plastique (epsilon_p) sur un element elem de S
 */
struct calcul_plasticite_elem {
    template<class TE>
    void operator()(TE &elem, Sst &S, Process &process, unsigned &i_elem) const {
        /// Test du signe de la fonction seuil
        TYPEREEL f,seq;
        Vec<TYPEREEL> sigma = S.t[process.temps->pt_cur].sigma[i_elem] ;
        TYPEREEL R_p = S.t[process.temps->pt_cur].R[i_elem];
        seq = calcul_contrainte_equivalente(S,sigma);
        //seq = std::pow(sigma[0],2) + std::pow(sigma[1],2) + std::pow(sigma[2],2) + std::pow(sigma[3],2) + std::pow(sigma[4],2) + std::pow(sigma[5],2) ;
        //seq = std::sqrt(seq);
        f = seq - R_p - S.matprop.R0;
        
        //std::cout << "Seuil :" << f << "    ";
        /// Si la fonction seuil est positive : correction de la plasticite cumulee
        if (f > 0) {
            /// Calcul de la variation de contrainte seq_dot
            Vec<TYPEREEL> sig_old;
            TYPEREEL seq_old, seq_dot;
            if (process.temps->pt_cur>1) {
                sig_old = S.t[process.temps->pt_cur-1].sigma[i_elem] ;
                seq_old = calcul_contrainte_equivalente(S,sig_old);
                //seq_old = std::pow(sig_old[0],2) + std::pow(sig_old[1],2) + std::pow(sig_old[2],2) + std::pow(sig_old[3],2) + std::pow(sig_old[4],2) + std::pow(sig_old[5],2) ;
                //seq_old = std::sqrt(seq_old);
                seq_dot = (seq - seq_old)/process.temps->dt;
            }
            else {
                seq_dot = seq/process.temps->dt;
            }
            ///Recherche de la nouvelle plasticite cumulee
            TYPEREEL K = seq_dot/(S.matprop.k_p*S.matprop.m_p);
            if (process.temps->pt_cur>1)
                recherche_plasticite_cumulee(elem.p,K,process.temps->dt,S.t[process.temps->pt_cur-1].p[i_elem],S.matprop.m_p);
            else {
                double zero = 0.0;
                recherche_plasticite_cumulee(elem.p,K,process.temps->dt,zero,S.matprop.m_p);
            }
            //*
            if(i_elem == 12415){
                std::cout << "sigma : " << sigma << std::endl;
                std::cout << "R0            : " << S.matprop.R0 << std::endl;
                std::cout << "R_p           : " << R_p << std::endl;
                std::cout << "f             : " << f << std::endl;
                std::cout << "sigma eq. new : " << seq << std::endl;
                std::cout << "sigma eq. old : " << seq_old << std::endl;
                std::cout << "sigma eq. dot : " << seq_dot << std::endl;
                std::cout << "p             : " << elem.p << std::endl;
            }//*/
            f_max = (f_max < f)? f : f_max ;
            p_max = (p_max < elem.p)? elem.p : p_max;
        } else {
            if (process.temps->pt_cur>1) 
                elem.p = S.t[process.temps->pt_cur-1].p[i_elem];
            else 
                elem.p = 0 ;
        }
        
        /// Calcul des deformations plastiques (epsilon_p)
        TYPEREEL sig_trace = (sigma[0]+sigma[1]+sigma[2])/3;
        for (int i = 0; i < DIM*(DIM+1)/2; i++) {
            if (seq != 0)
                elem.epsilon_p[0][i] = 3/2*elem.p*((i<3)? sigma[i]-sig_trace: sigma[i])/seq ;
            else 
                elem.epsilon_p[0][i] = 0 ;
        }
        //*
        if(i_elem == 12415){
            std::cout << "epsilon_p : " << elem.epsilon_p << std::endl;
            std::cout << "f new     : " << seq-(S.matprop.k_p * std::pow( elem.p , S.matprop.m_p ) + S.matprop.R0) << std::endl;
            std::cout << "sigma new : " << (sigma - 200000*elem.epsilon_p[0]) << std::endl;
        }//*/
        i_elem++;    ///Sert d'index d'iteration
    }
};


/* PRISE EN COMPTE ERREUR_PLASTICITE NON IMPLEMENTEE
// Calcul de l'erreur en plasticite
template<class TE, class TTE, class ME>
void calcul_erreur_plasticite(TE &erreur, TTE &p, ME &ep, TTE &pold, ME &epold) {
    erreur.set(0.0);
    for (unsigned i=0;i<p.size();++i) {
        erreur[0] += (p[i] - pold[i])*(p[i] - pold[i]);
        erreur[1] += (p[i] + pold[i])*(p[i] + pold[i]);
        for (unsigned j=0;j<ep[i].size();++j) {
            erreur[0] += (ep[i][j] - epold[i][j])*(ep[i][j] - epold[i][j]);
            erreur[1] += (ep[i][j] + epold[i][j])*(ep[i][j] + epold[i][j]);
        }
    }
}

// Calcul de la contribution de la plasticite a l'erreur
struct calcerror_plastique {
    template<class SST, class TE>
    void operator()(SST &S, Param &process, TE &erreur_plastiquet) const {
        if (S.plastique) {
            cout << "Contribution de la plasticite a l'erreur non code !" << endl ;
        }
    }
};//*/


// Etape lineaire ------------------------------------------------------------------------------------------
void calcul_plasticite(Sst &S, Process &process) {
    /// Reconstruction et stockage de sigma sur les Sst
    calcul_sigma_sst(S,process);
    /// Calcul et stockage de d et Yd
    //calcul_forces_thermo_endo(S,process) ;    ENDOMMAGEMENT NON IMPLEMENTE
    /// Calcul et stockage de l'ecrouissage R sur les Sst
    calcul_ecrouissage_sst(S,process);
    
    /// Calcul de p et epsilon_p sur chaque element
    unsigned i_elem = 0 ;
    /* ENDOMMAGEMENT NON IMPLEMENTE
    if (S.endommageable) 
        apply(S.mesh.elem_list,calcul_plasticite_couplee_elem(),S,process,i_elem);
    else*/ 
        apply(S.mesh->elem_list,calcul_plasticite_elem(),S,process,i_elem);
    //std::cout << "max(f) : " << f_max << "        max(p) : " << p_max << std::endl;
    
    /// On stocke ensuite les deux variables de plasticite dans un coin
    /* PRISE EN COMPTE ERREUR_PLASTICITE NON IMPLEMENTEE
    Vec<TYPEREEL> p_old ;
    Vec< Vec<TYPEREEL> > ep_old ;
    if (process.latin->contribution_error == "sst") {
        p_old = S.t[process.temps->pt_cur].p ;
        ep_old = S.t[process.temps->pt_cur].epsilon_p ;
    }
    //*/
    /// Stockage de p et epsilon_p
    S.t[process.temps->pt_cur].p.resize(0);
    apply(S.mesh->elem_list,stockage_plasticite_cumulee(),S.t[process.temps->pt_cur].p);
    S.t[process.temps->pt_cur].epsilon_p.resize(0);
    apply(S.mesh->elem_list,stockage_deformation_plastique(),S.t[process.temps->pt_cur].epsilon_p);
    /// Calcul et prise en compte
    /* PRISE EN COMPTE ERREUR_PLASTICITE NON IMPLEMENTEE
    if (process.latin->contribution_error == "sst")
        calcul_erreur_plasticite(S.t[process.temps->pt_cur].erreur_plasticite,S.t[process.temps->pt_cur].p,S.t[process.temps->pt_cur].epsilon_p,p_old,ep_old);
    //*/
}

#endif //PLASTICITY_FUNCTIONS_H
