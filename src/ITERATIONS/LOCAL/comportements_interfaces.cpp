#include "../../DEFINITIONS/main_typedef.h"
#include "../../DEFINITIONS/Interface.h"
#include "../../UTILITAIRES/utilitaires.h"
#include <cmath>

// #define CHECK_COMPORTEMENT_INTERFACES   /// A definir pour que les comportements d'interface soit verifies
#define CHECK(err_name,err_def) Scalar err_name = norm_2(err_def); if(err_name > err_max){ PRINT(err_def); }

Scalar err_max = 1e-6;  /// Erreur maximale toleree sur la verification des equations


void Interface::NodalState::check_ddr(){
    /// Erreur en equilibre des efforts
    CHECK(err_equilibre,Fchap1+Fchap2);
    /// Erreur en direction de recherche sur le cote 1     <--- A MODIFIER POUR PRISE EN COMPTE PRECHARGE
    CHECK(err_ddr1,Fchap1-F1-h1*(Wpchap1 - Wp1));
    /// Erreur en direction de recherche sur le cote 2     <--- IDEM
    CHECK(err_ddr2,Fchap2-F2-h2*(Wpchap2 - Wp2));
}

void Interface::NodalState::comportement_parfait(){
    /// Derivee temporelle de l'epaisseur imposee
    Point dEp_imposee = (Ep_imposee - old_Ep_imposee)/dt;
    /// Construction de l'operateur (ne pas oublier que les normales sont identiques et que les matrices diagonales dans une base commune sont cools)
    Interface::LocalOperator h;
    h.kn = 1.0/(k1.kn + k2.kn);
    h.kt = 1.0/(k1.kt + k2.kt);
    h.n = n1;
    Wpchap1 = h*( (k1*Wp1 + k2*Wp2) - (F1 + F2) - k2*dEp_imposee);
    //Wpchap2 = Wpchap1 + dEp_imposee;
    Fchap1 = F1 + k1*(Wpchap1 - Wp1) - (Precharge - old_Precharge);
    Fchap2 = -1.0*Fchap1;
    Wpchap2 = Wp2 + h2*(Fchap2 - F2);
}

void Interface::NodalState::check_comportement_parfait(){
    /// Erreur en saut de vitesse
    CHECK(err_dWp,Wpchap2 - Wpchap1 - (Ep_imposee - old_Ep_imposee)/dt);
}


void Interface::NodalState::comportement_elastique(){
    Point dEp_imposee = (Ep_imposee - old_Ep_imposee)/dt;
    /// Construction de l'operateur (ne pas oublier que les normales sont identiques et que les matrices diagonales dans une base commune sont cools)
    Interface::LocalOperator Htilde;
    Htilde.kn = 1.0/(1.0/dt + K.kn*(h1.kn + h2.kn));
    Htilde.kt = 1.0/(1.0/dt + K.kt*(h1.kt + h2.kt));
    Htilde.n = n1;
    /// Calcul de la variation d'epaisseur due a l'elasticite
    Ep_elastique = Htilde*( (Wp2 - Wp1) - (h2*F2 - h1*F1) - dEp_imposee + old_Ep_elastique/dt);
    /// Calcul des autres valeurs
    Fchap1 = (1-d)*(K*Ep_elastique) - (Precharge - old_Precharge);
    Fchap2 = -1.0*Fchap1;
    Wpchap1 = h1*(Fchap1 - F1) + Wp1;
    Wpchap2 = h2*(Fchap2 - F2) + Wp2;
    /*if(id < 5){
        PRINT(F1);
        PRINT(F2);
        PRINT(Wp1);
        PRINT(Wp2);
        PRINT(Ep_elastique);
        PRINT(Fchap1);
        PRINT(Fchap2);
        PRINT(Wpchap1);
        PRINT(Wpchap2);
    }*/
}


void Interface::NodalState::check_comportement_elastique(){
    /// Erreur en saut de vitesse
    CHECK(err_dWpchap,Wpchap2 - Wpchap1 - (Ep_imposee - old_Ep_imposee)/dt - (Ep_elastique - old_Ep_elastique)/dt);
    /// Erreur en rdc elastique
    CHECK(err_K,Fchap1 - K*Ep_elastique);
}


void Interface::NodalState::comportement_cohesif(){
    Point dEp_imposee = (Ep_imposee - old_Ep_imposee)/dt;
    const unsigned max_iter = 1000;
    Scalar epsilon = 1.0e-6;
    Scalar dmax = 1.0-1.0e-3;
    
    Scalar Yn,Yt,Y,Fchap1n,Fchap1t;
    Scalar alpha = matprop->alpha;
    Scalar gamma = matprop->gamma;
    Scalar n = matprop->n;
    Scalar Yc = matprop->Yc;
    Scalar Yo = matprop->Yo;
    Scalar d_coeff = std::pow( n/(n+1) / (Yc-Yo) ,n);
    
    Interface::LocalOperator Ktilde;
    Ktilde.kn = 1.0/dt + (1-d)*K.kn*(h1.kn + h2.kn);
    Ktilde.kt = 1.0/dt + (1-d)*K.kt*(h1.kt + h2.kt);
    Ktilde.n = n1;
    Interface::LocalOperator Htilde;
    Htilde.kn = 1.0/Ktilde.kn;
    Htilde.kt = 1.0/Ktilde.kt;
    Htilde.n = n1;
    
    /// On realise un calcul purement elastique qui servira de predicateur
    comportement_elastique();
    Point F_elastique = Fchap1;
    Point delta_Ep;
    Point residu = Ktilde * delta_Ep - F_elastique;
    Scalar error = norm_2(residu);
    
    unsigned i_iter = 0;
    while(i_iter < max_iter and error > epsilon){
        /// Correction du jeu du a l'elasticite
        delta_Ep = delta_Ep - Htilde*residu;
        
        /// Calcul du coefficient de degradation associe
        Fchap1 = (1-d)*(K*delta_Ep);
        Fchap1n = dot(Fchap1,n1);
        Fchap1t = norm_2(Fchap1-Fchap1n*n1);
        Yn = 0.5/(K.kn*(1-d)*(1-d)) * (Fchap1n>0.0) * Fchap1n * Fchap1n;
        Yt = 0.5/(K.kt*(1-d)*(1-d)) * Fchap1t * Fchap1t;
        Y = std::pow( std::pow(Yn,alpha) + std::pow(gamma*Yt,alpha) ,1.0/alpha);
        d = d_coeff * std::pow( (Y > Yo) * std::abs(Y - Yo) , n);
        
        /// Mise a jour des operateurs
        Ktilde.kn = 1.0/dt + (1-d)*K.kn*(h1.kn + h2.kn);
        Ktilde.kt = 1.0/dt + (1-d)*K.kt*(h1.kt + h2.kt);
        Htilde.kn = 1.0/Ktilde.kn;
        Htilde.kt = 1.0/Ktilde.kt;
        
        /// Calculs des grandeurs de controle
        residu = Ktilde * delta_Ep - F_elastique;
        error = norm_2(residu);
        i_iter++;
    }
    
    d = std::max(d,old_d);
    if(d > dmax){
        /// Ruine du materiau
        comportement = true;
    }
    else{
        /// Le materiau survit au chargement
        Fchap1 = (1-d)*(K*delta_Ep);
        Fchap2 = -1.0*Fchap1;
        Wpchap1 = h1*(Fchap1 - F1) + Wp1;
        Wpchap2 = h2*(Fchap2 - F2) + Wp2;
    }
}


void Interface::NodalState::check_comportement_cohesif(){
    /// A definir
}


void Interface::NodalState::comportement_cassable(){
    /// !!! On suppose qu'un comportement parfait ou elastique a deja ete calcule
    /// si la convergence du calcul iteratif est OK, on met à jour le comportement des elements qui ne sont pas deja casse
    if (interface.convergence >= 0 and comportement == false){
        /// test contact normal : apres le calcul en supposant la cohesion des 2 cotes (c.f. comportement_local_interface, plus bas)
        /// on verifie si (< -Fchap1_n >+ / Fcr_n)^2 + (Fchap1_t / Fcr_t)^2 > 1 avec < >+ la partie positive
        /// que l'on reecrit k * N2 + T2 > R en multipliant par Fcr_t ^ 2
        Scalar N = dot(Fchap1,n1);  /// Effort normal (de 1 vers 2)
        Point T = ProjT(Fchap1,n1); /// Effort tangentiel
        Scalar N2 = N*N*(N < 0);    /// norme au carre de l'effort normal en traction
        Scalar T2 = dot(T,T);       /// norme au carre de l'effort tangentiel
        Scalar R = (matprop->Fcr_t*matprop->Fcr_t);
        Scalar k = R/(matprop->Fcr_n*matprop->Fcr_n);
        if (k * N2 + T2 > R){
            ///David dit de mettre 10% de plus (A QUOI ???)
            /// on a franchi la limite de rupture
            comportement = true;
            interface.convergence++;
        }
    }
}


void Interface::NodalState::check_comportement_cassable(){
    /// A definir
}


void Interface::NodalState::comportement_contact_parfait(){
    /// Quelques grandeurs utiles
    Scalar Ep_n = dot(Ep_imposee,n1);
    Point dEp = (Ep_imposee -old_Ep_imposee)/dt;
    Scalar dEp_n = dot(dEp,n1);
    Point dEp_t = ProjT(dEp,n1);
    Scalar dPrecharge_n = dot(Precharge,n1) - dot(old_Precharge,n1);
    Point dPrecharge_t = ProjT(Precharge,n1) - ProjT(old_Precharge,n1);
    
    /// Test de contact
    Scalar dWchap_n = dot(n1,(old_Wchap2-old_Wchap1)) + dt*(dot(n1,(Wp2-Wp1)) - (h2.kn*dot(n1,F2) - h1.kn*dot(n1,F1)));
    if (dWchap_n > Ep_n) {
        /// separation des bords
        Fchap1 = 0.0;
        Fchap2 = 0.0;
        Wpchap1 = Wp1 - h1*F1;
        Wpchap2 = Wp2 - h2*F2;
    }
    else{
        /// collision des bords
        /// Calcul des valeurs normales : on conserve la partie normale d'un calcul d'interface parfaite
        comportement_parfait();
        Scalar Wpchap1n = dot(n1,Wpchap1);
        Scalar Wpchap2n = dot(n1,Wpchap2);
        Scalar Fchap1n = dot(n1,Fchap1);
        Scalar Fchap2n = -1.0*Fchap1n;
        
        /// Test de glissement adherence
        /// Effort tangentiel
        Point T = ((ProjT(Wp2,n1) - ProjT(Wp1,n1)) - (h2.kt*ProjT(F2,n1) - h1.kt*ProjT(F1,n1)) - dEp_t) / (h1.kt+h2.kt) - dPrecharge_t;
        Scalar normT = norm_2(T);
        /// Limite d'adherence connaissant l'effort normal
        Scalar g = coeffrottement*std::abs(Fchap1n);
        
        /// Parties tangentielles de la solution
        Point Fchap1t,Fchap2t,Wpchap1t,Wpchap2t;
        
        if (normT <= g) {
            /// adherence
            Fchap1t = T;
        } else if (normT > g) {
            /// glissement
            Fchap1t = T*g/normT;
        }
        Fchap2t = -1.0*Fchap1t;
        Wpchap1t = ProjT(Wp1,n1)+h1.kt*(Fchap1t-ProjT(F1,n1));
        Wpchap2t = ProjT(Wp2,n1)+h2.kt*(Fchap2t-ProjT(F2,n1));
        
        /// Assemblage des parties normales et tangentielles
        Wpchap1 = Wpchap1n*n1 + Wpchap1t;
        Wpchap2 = Wpchap2n*n1 + Wpchap2t;
        Fchap1 = Fchap1n*n1 + Fchap1t;
        Fchap2 = Fchap2n*n1 + Fchap2t;
    }
}


void Interface::NodalState::check_comportement_contact_parfait(){
    /// A definir
}


void Interface::NodalState::comportement_contact_elastique(){
    /// Quelques grandeurs utiles
    Scalar Ep_n = dot(Ep_imposee,n1);
    Point dEp = (Ep_imposee -old_Ep_imposee)/dt;
    Scalar dEp_n = dot(dEp,n1);
    Point dEp_t = ProjT(dEp,n1);
    Scalar dPrecharge_n = dot(Precharge,n1) - dot(old_Precharge,n1);
    Point dPrecharge_t = ProjT(Precharge,n1) - ProjT(old_Precharge,n1);
    
    /// Test de contact
    Scalar dWchap_n = dot(n1,(old_Wchap2-old_Wchap1)) + dt*(dot(n1,(Wp2-Wp1)) - (h2.kn*dot(n1,F2) - h1.kn*dot(n1,F1)));
    if (dWchap_n > Ep_n) {
        /// separation des bords
        Fchap1 = 0.0;
        Fchap2 = 0.0;
        Wpchap1 = Wp1 - h1*F1;
        Wpchap2 = Wp2 - h2*F2;
        Ep_elastique.set(0.0);
    }
    else{
        /// collision des bords
        /// Calcul des valeurs normales : on conserve la partie normale d'un calcul d'interface elastique
        comportement_elastique();
        Scalar Wpchap1n = dot(n1,Wpchap1);
        Scalar Wpchap2n = dot(n1,Wpchap2);
        Scalar Fchap1n = dot(n1,Fchap1);
        Scalar Fchap2n = -1.0*Fchap1n;
        
        /// Test de glissement adherence
        /// Effort tangentiel
        Point T = ((ProjT(Wp2,n1) - ProjT(Wp1,n1)) - (h2.kt*ProjT(F2,n1) - h1.kt*ProjT(F1,n1)) - dEp_t + ProjT(old_Ep_elastique,n1)/dt ) / (h1.kt + h2.kt + 1.0/(K.kt*dt)) - dPrecharge_t;
        Scalar normT = norm_2(T);
        /// Limite d'adherence connaissant l'effort normal
        Scalar g = coeffrottement*std::abs(Fchap1n);
        
        /// Parties tangentielles de la solution
        Point Fchap1t,Fchap2t,Wpchap1t,Wpchap2t;
        
        if (normT <= g) {
            /// adherence
            Fchap1t = T;
        } 
        else if (normT > g) {
            /// glissement
            Fchap1t = T*g/normT;
        }
        Fchap2t = -1.0*Fchap1t;
        Wpchap1t = ProjT(Wp1,n1)+h1.kt*(Fchap1t-ProjT(F1,n1));
        Wpchap2t = ProjT(Wp2,n1)+h2.kt*(Fchap2t-ProjT(F2,n1));
        
        /// Assemblage des parties normales et tangentielles
        Wpchap1 = Wpchap1n*n1 + Wpchap1t;
        Wpchap2 = Wpchap2n*n1 + Wpchap2t;
        Fchap1 = Fchap1n*n1 + Fchap1t;
        Fchap2 = Fchap2n*n1 + Fchap2t;
        /// Calcul de la deformation elastique
        Interface::LocalOperator H; /// Operateur de souplesse
        H.kn = 1/K.kn;
        H.kt = 1/K.kt;
        Ep_elastique = H*(Fchap1 + (Precharge - old_Precharge));
    }
}


void Interface::NodalState::check_comportement_contact_elastique(){
    /// A definir
}



void comportement_local_interface(Interface &Inter, unsigned pt, Scalar dt){
    const unsigned nb_nodes = Inter.side[0].nodeeq.size();
    Interface::NodalState node(Inter,pt,dt);
    
    /// Interface parfaite
    if(Inter.comp == Interface::comp_parfait){
        for(unsigned i_node = 0; i_node < nb_nodes; i_node++){
            /// On recupere les valeurs au noeud
            node.set_node(i_node);
            /// On lui applique le comportement adequat
            node.comportement_parfait();
            /// On peut verifier si les equations sont respectees
            #ifdef CHECK_COMPORTEMENT_INTERFACES
            node.check_ddr();
            node.check_comportement_parfait();
            #endif
            /// On sauvegarde les resultats
            node.store_results();
        }
    }
    /// Interface elastique
    else if(Inter.comp == Interface::comp_elastique){
        for(unsigned i_node = 0; i_node < nb_nodes; i_node++){
            node.set_node(i_node);
            node.comportement_elastique();
            #ifdef CHECK_COMPORTEMENT_INTERFACES
            node.check_ddr();
            node.check_comportement_elastique();
            #endif
            node.store_results();
        }
    }
    /// Interface parfaite cassable
    else if(Inter.comp == Interface::comp_cassable_parfait){
        for(unsigned i_node = 0; i_node < nb_nodes; i_node++){
            node.set_node(i_node);
            /// On verifie si le noeud doit casser
            if(not node.comportement){
                node.comportement_parfait();
                #ifdef CHECK_COMPORTEMENT_INTERFACES
                node.check_ddr();
                node.check_comportement_parfait();
                #endif
                node.comportement_cassable();
                #ifdef CHECK_COMPORTEMENT_INTERFACES
                node.check_ddr();
                node.check_comportement_cassable();
                #endif
            }
            /// S'il a casse, on applique le contact
            if(node.comportement){
                node.comportement_contact_parfait();
                #ifdef CHECK_COMPORTEMENT_INTERFACES
                node.check_ddr();
                node.check_comportement_contact_parfait();
                #endif
            }
            node.store_results();
        }
    }
    /// Interface elastique cassable
    else if(Inter.comp == Interface::comp_cassable_elastique){
        for(unsigned i_node = 0; i_node < nb_nodes; i_node++){
            node.set_node(i_node);
            /// On verifie si le noeud doit casser
            if(not node.comportement){
                node.comportement_elastique();
                #ifdef CHECK_COMPORTEMENT_INTERFACES
                node.check_ddr();
                node.check_comportement_elastique();
                #endif
                node.comportement_cassable();
                #ifdef CHECK_COMPORTEMENT_INTERFACES
                node.check_ddr();
                node.check_comportement_cassable();
                #endif
            }
            /// S'il a casse, on applique le contact
            if(node.comportement){
                node.comportement_contact_elastique();
                #ifdef CHECK_COMPORTEMENT_INTERFACES
                node.check_ddr();
                node.check_comportement_contact_elastique();
                #endif
            }
            node.store_results();
        }
    }
    /// Interface contact parfait
    else if(Inter.comp == Interface::comp_contact_parfait){
        for(unsigned i_node = 0; i_node < nb_nodes; i_node++){
            node.set_node(i_node);
            node.comportement_contact_parfait();
            #ifdef CHECK_COMPORTEMENT_INTERFACES
            node.check_ddr();
            node.check_comportement_contact_parfait();
            #endif
            node.store_results();
        }
    }
    /// Interface contact elastique
    else if(Inter.comp == Interface::comp_contact_elastique){
        for(unsigned i_node = 0; i_node < nb_nodes; i_node++){
            node.set_node(i_node);
            node.comportement_contact_elastique();
            #ifdef CHECK_COMPORTEMENT_INTERFACES
            node.check_ddr();
            node.check_comportement_contact_elastique();
            #endif
            node.store_results();
        }
    }
    /// Interface cohesive
    else if(Inter.comp == Interface::comp_cohesive){
        for(unsigned i_node = 0; i_node < nb_nodes; i_node++){
            node.set_node(i_node);
            if(not node.comportement){
                node.comportement_cohesif();
            }
            if(node.comportement){
                node.comportement_contact_parfait();    /// L'interface est ruinee d'ou un contact parfait et non elastique
                #ifdef CHECK_COMPORTEMENT_INTERFACES
                node.check_ddr();
                node.check_comportement_contact_parfait();
                #endif
            }
            node.store_results();
        }
    }
}