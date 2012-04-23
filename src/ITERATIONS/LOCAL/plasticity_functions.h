#ifndef PLASTICITY_FUNCTIONS_H
#define PLASTICITY_FUNCTIONS_H


#include "../../DEFINITIONS/Process.h"
#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/Interface.h"
#include "../../UTILS/utils_2.h"
//#include <boost/bind/bind_template.hpp>



/** Recuperation de la deformation plastique du "stock" vers la formulation
 * 
 */
struct chargement_deformation_plastique {
    template<class TE, class ME>
    void operator()(TE &elem, ME &ep, unsigned &i_elem) const {
        elem.epsilon_p[0] = ep[i_elem] ;
        ++i_elem;    ///Sert d'index d'iteration sur les elements
    }
};


/**Stockage des deformations plastiques calculees
 * Les deformations plastiques calculees sont stockees dans le vecteur Time::epsilon_p du pas de temps courant
 *//*
struct stockage_deformation_plastique {
    template<class TE,class ME>
    void operator()(TE &elem, ME &ep) const {
        ep.push_back(elem.epsilon_p[0]);
    }
};
*/

/**Stockage de la plasticite cumulee calculee
 * La plasticite cumulee calculee est stockee dans le vecteur Time::p du pas de temps courant
 */
struct chargement_plasticite_cumulee {
    template<class TE,class ME>
    void operator()(TE &elem, ME &p,unsigned &i_elem) const {
        elem.p = p[i_elem];
        i_elem ++;
    }
};


/**Stockage de la plasticite cumulee calculee
 * La plasticite cumulee calculee est stockee dans le vecteur Time::p du pas de temps courant
 *//*
struct stockage_plasticite_cumulee {
    template<class TE,class ME>
    void operator()(TE &elem, ME &p) const {
        p.push_back(elem.p);
    }
};*/


/**Stockage des contraintes calculees
 * Les contraintes calculee sont stockees dans le vecteur Time::sigma du pas de temps courant
 */
struct stockage_sigma {
    template<class TE,class ME>
    void operator()(TE &elem, ME &s, unsigned &i_elem) const {
        s[i_elem] = elem.sigma[0];
        ++i_elem;    ///Sert d'index d'iteration
    }
};

/** Calcul de la contrainte equivalente
 * Retourne la contrainte equivalente (type Von Mises) sur un element de S soumis a sigma
 */
TYPEREEL calcul_contrainte_equivalente(const Sst &S,const Vec<TYPEREEL> &sigma);    //CALCULE DANS LA FORMULATION (elem.sigma_von_mises)

struct compare_vm {
    template<class TE,class ME>
    void operator()(TE &elem,Sst &S, ME &s, unsigned &i_elem) const {
        if(S.id == 0 and i_elem == 281) std::cout << "********************************************************" << std::endl;
        if(S.id == 0 and i_elem == 281) std::cout << "sigma                : " << elem.sigma[0] << std::endl;
        if(S.id == 0 and i_elem == 281) std::cout << "sigma_eq formulation : " << elem.sigma_von_mises << std::endl;
        if(S.id == 0 and i_elem == 281) std::cout << "sigma_eq calcule     : " << calcul_contrainte_equivalente(S,elem.sigma[0]) << std::endl;
        if(S.id == 0 and i_elem == 281) std::cout << "********************************************************" << std::endl << std::endl;
        ++i_elem;    ///Sert d'index d'iteration
    }
};

/*
void calcul_sigma_sst(Sst &S, Process &process);


void calcul_ecrouissage_sst(Sst &S, Process &process);
*/

/** Recherche de la plasticite cumulee
 * Le premier argument servira de stockage pour le resultat
 */
/*
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
*/

//TYPEREEL seq_max,p_max,f_max;

/** Calcul de la plasticitee (non couplee a l'endommagement)
 * Calcul de la plasticite cumulee (p) et de la deformation plastique (epsilon_p) sur un element elem de S
 *//*
struct calcul_plasticite_elem {
    template<class TE>
    void operator()(TE &elem, Sst &S, Process &process, unsigned &i_elem) const {
        /// Predicateur elastique
        Vec<TYPEREEL> sigma = S.t[process.temps->pt].sigma[i_elem] ;
        TYPEREEL sigma_eq = calcul_contrainte_equivalente(S,sigma);     /// Contrainte equivalente
        TYPEREEL R_p = S.t[process.temps->pt].R[i_elem];                /// Ecrouissage actuel
        //std::cout << R_p << " , " << elem.R_p << std::endl;
        //TYPEREEL sigma_eq = elem.sigma_von_mises;                       /// Contrainte equivalente
        TYPEREEL f = sigma_eq - R_p;                                    /// Valeur de la fonction seuil
        
        /// Si la fonction seuil est positive : correction de la plasticite cumulee - algorithme de retour radial
        if (f > 0) {
            /// Calcul de la variation de p (Methode de Newton)
            TYPEREEL mu = S.matprop.elastic_modulus_1/(2.*(1.+S.matprop.poisson_ratio_12));      /// Second coefficient de Lame
            TYPEREEL dPn = 0;
            TYPEREEL erreur = f/R_p;   /// On adimensionne par la limite d'elasticite actuelle
            while(erreur > 1e-6){
                dPn = dPn + R_p*erreur/(3*mu+S.matprop.m_p*S.matprop.k_p*std::pow(S.t[1].p[i_elem]+dPn,S.matprop.m_p-1.)); /// Attention au signe
                erreur = (sigma_eq - 3*mu*dPn - S.matprop.R0 - S.matprop.k_p*std::pow(S.t[1].p[i_elem],S.matprop.m_p))/R_p;
            }
            //dPn = f/(3*mu+S.matprop.m_p*S.matprop.k_p*std::pow(elem.p,S.matprop.m_p-1));
            S.t[1].p[i_elem] += dPn;
            
            /// Trace de sigma
            TYPEREEL trace_sigma = 0;
            for(unsigned i = 0; i < DIM; i++)
                trace_sigma += sigma[i];
            /// Partie deviatorique de sigma
            Vec<TYPEREEL,DIM*(DIM+1)/2> sigma_dev;
            for(unsigned i = 0; i < DIM; i++)               /// Boucle sur les termes diagonaux
                sigma_dev[i] = sigma[i] - trace_sigma/3;
            for(unsigned i = DIM; i < DIM*(DIM+1)/2; i++)   /// Boucle sur les termes hors diagonale
                sigma_dev[i] = sigma[i];
        
            /// Calcul de la variation de epsilon_p
            Vec<TYPEREEL,DIM*(DIM+1)/2> dEpsPn = 3/2*dPn*sigma_dev/sigma_eq;
            elem.epsilon_p[0] += dEpsPn;
            
            /// Stockage des resultats
            //S.t[process.temps->pt].p[i_elem] = elem.p;
            //S.t[process.temps->pt].epsilon_p[i_elem] = elem.epsilon_p[0];
            Vec<TYPEREEL,DIM*(DIM+1)/2> sigma_new = sigma-2*mu*dEpsPn;
            //elem.sigma_von_mises = calcul_contrainte_equivalente(S,S.t[process.temps->pt].sigma[i_elem]);
            
            //*   TEST
            if(S.id == 0 and i_elem == 281){
                //std::cout << "Numero de l'element : " << i_elem << std::endl;
                std::cout << "sigma eq. 1   : " << sigma_eq << std::endl;
                std::cout << "sigma eq. 2   : " << elem.sigma_von_mises << std::endl;
                std::cout << "f old         : " << f << std::endl;
                //std::cout << "Erreur init   : " << f/(S.matprop.R0+R_p) << std::endl;
                //std::cout << "Erreur fin    : " << (S.matprop.R0+R_p)*erreur << std::endl;
                std::cout << "dp              : " << dPn << std::endl;
                std::cout << "p               : " << S.t[1].p[i_elem] << std::endl;
                std::cout << "old R(p)        : " << R_p << std::endl;
                std::cout << "new R(p)        : " << S.matprop.k_p * std::pow( S.t[1].p[i_elem] , S.matprop.m_p ) + S.matprop.R0 << std::endl;
                std::cout << "f new           : " << calcul_contrainte_equivalente(S,sigma_new)-(S.matprop.k_p * std::pow( S.t[1].p[i_elem] , S.matprop.m_p ) + S.matprop.R0) << std::endl;
                //std::cout << "sigma dev : " << sigma_dev << std::endl;
                //std::cout << "dEpsPn    : " << dEpsPn << std::endl;
                //std::cout << "epsilon_p : " << elem.epsilon_p << std::endl;
                std::cout << "sigma old : " << sigma << std::endl;
                std::cout << "sigma new : " << sigma_new << std::endl;
                //std::cout << "Test0 : " << dPn-f/(3*mu+S.matprop.k_p) << std::endl;
                std::cout << "Test Newton      : " << sigma_eq - 3*mu*dPn- (S.matprop.k_p * std::pow( S.t[1].p[i_elem] , S.matprop.m_p ) + S.matprop.R0) << std::endl;
                std::cout << "Test Trace(dEps) : " << dEpsPn[0]+dEpsPn[1]+dEpsPn[2] << std::endl;
                //std::cout << "Test2 : " << sigma_eq-3*mu*dPn << " , " << calcul_contrainte_equivalente(S,S.t[process.temps->pt].sigma[i_elem]) << std::endl;
                //assert(0);
            }
            //*//*
        }
        i_elem++;    ///Sert d'index d'iteration sur les elements
    }
};
*/

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


// Etape locale ------------------------------------------------------------------------------------------
void calcul_plasticite(Sst &S, Process &process);

#endif //PLASTICITY_FUNCTIONS_H
