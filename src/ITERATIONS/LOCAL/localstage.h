#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/Interface.h"
#include "comportements_interfaces.h"
#include "../LOCAL/plasticity_functions.h"

/** \defgroup etape_locale Etape locale
\ingroup LATIN
\brief Description des phases de l'étape locale
 
 Cette étape de la stratégie LATIN est constituée des phases suivantes :
 - pour chaque intervalle de temps on fait :
   - étape locale pour chaque interface : etape_locale_inter
   - étape locale pour chaque sous-structure : etape_locale_sst
La procédure qui gère la boucle locale est écrite dans etape_locale().
*/

/** \ingroup  etape_locale
\brief Etape locale pour les interfaces
 
Selon le type d'interface on appelle la procédure correspondante :
- interface de type vitesse (Interface::comp == "vit") : compt_CL_vit()
- interface de type effort (Interface::comp == "effort") : compt_CL_eff()
- interface de type symétrie (Interface::comp == "sym") : compt_CL_sym()
- interface de type vitesse normale donnée  (Interface::comp == "vit_normale") : compt_CL_sym()
- interface de type parfait (Interface::comp == "Parfait") : compt_parfait()
- interface de type contact avec ou sans frottement (Interface::comp == "Contact") : compt_contact()
- interface de type jeu contact avec ou sans frottement (Interface::comp == "Contact_jeu" ou "Contact_physique") : compt_contact()
- interface de type epaisse contact avec ou sans frottement (Interface::comp == "Contact_ep" ) : compt_contact()
- interface de type cohésive (Interface::comp == "Cohesive") : compt_cohesif()
 
Il suffit donc de rajouter un comportement dans cette procédure et programmer la fonction correspondante pour ajouter un comportement d'interface.
*/
struct etape_locale_inter {
    void operator()(Interface &Inter,Process &process) const {
//         std::cout << Inter.type << std::endl;
//         std::cout << Inter.comp << std::endl;
        if (Inter.type=="Ext") {
            if(Inter.comp=="vit" or Inter.comp=="depl") {
                compt_CL_vit(Inter,process.temps->pt);
            } else if (Inter.comp=="effort" or Inter.comp=="effort_normal") {
                compt_CL_eff(Inter,process.temps->pt);
/*            } else if (Inter.comp=="effort_normal") {
                compt_CL_eff_normal(Inter,process.temps->pt);*/
            } else if (Inter.comp=="sym") {
                compt_CL_sym(Inter,process.temps->pt);
            } else if (Inter.comp=="vit_normale" or Inter.comp == "depl_normal") {
                compt_CL_sym(Inter,process.temps->pt);
            } else if (Inter.comp=="periodique") {
                compt_parfait(Inter,process.temps->pt);
            }
        } else if (Inter.type=="Int") {
            if (Inter.comp=="Parfait") {
                compt_parfait(Inter,process.temps->pt);
            } else if (Inter.comp=="Contact" or Inter.comp=="Contact_jeu" or Inter.comp=="Contact_jeu_physique") {
                compt_contact(Inter,*process.temps);
            } else if (Inter.comp=="Contact_ep" ) {
                compt_contact_ep(Inter,*process.temps);
            } else if (Inter.comp=="Cohesive") {
                compt_cohesif(Inter,*process.temps);
            } else if (Inter.comp=="Jeu_impose") {
                compt_jeu_impose(Inter,*process.temps);
            } else if (Inter.comp=="Breakable") {
                compt_breakable(Inter,*process.temps);
            } else {
                std::cout<< "Erreur : comportement interface non reconnu" << endl;
                assert(0);
            }
        } else {
            std::cout<< "Erreur : comportement interface non reconnu" << endl;
            assert(0);
        }
    }
};


/** \ingroup   etape_locale
\brief Etape locale pour les sous-structures.
*/
struct etape_locale_sst {
    void operator()(Sst &S,Process &process) const {
        if (S.plastique) calcul_plasticite(S,process);
    }
};


/** \ingroup  etape_locale
\brief Procédure principale pour l'étape locale
*/
void etape_locale(Vec<VecPointedValues<Interface > > &Inter,Vec<VecPointedValues<Sst > > &S,Process &process) {
    apply(Inter,etape_locale_inter(),process);
    apply(S,etape_locale_sst(),process);
};

