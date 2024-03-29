#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/Interface.h"
#include "comportements_interfaces.h"
#include "plasticite.h"
#include "endommagement.h"
#include "mesomodele.h"

/** \defgroup etape_locale Etape locale
\ingroup LATIN
\brief Description des phases de l'�tape locale
 
 Cette �tape de la strat�gie LATIN est constitu�e des phases suivantes :
 - pour chaque intervalle de temps on fait :
   - �tape locale pour chaque interface : etape_locale_inter
   - �tape locale pour chaque sous-structure : etape_locale_sst
La proc�dure qui g�re la boucle locale est �crite dans etape_locale().
*/

/** \ingroup  etape_locale
\brief Etape locale pour les interfaces
 
Selon le type d'interface on appelle la proc�dure correspondante :
- interface de type vitesse (Interface::comp == "vit") : compt_CL_vit()
- interface de type effort (Interface::comp == "effort") : compt_CL_eff()
- interface de type sym�trie (Interface::comp == "sym") : compt_CL_sym()
- interface de type vitesse normale donn�e  (Interface::comp == "vit_normale") : compt_CL_sym()
- interface de type parfait (Interface::comp == "Parfait") : compt_parfait()
- interface de type contact avec ou sans frottement (Interface::comp == "Contact") : compt_contact()
- interface de type jeu contact avec ou sans frottement (Interface::comp == "Contact_jeu" ou "Contact_physique") : compt_contact()
- interface de type epaisse contact avec ou sans frottement (Interface::comp == "Contact_ep" ) : compt_contact()
- interface de type coh�sive (Interface::comp == "Cohesive") : compt_cohesif()
 
Il suffit donc de rajouter un comportement dans cette proc�dure et programmer la fonction correspondante pour ajouter un comportement d'interface.
*/
struct etape_locale_inter {
    void operator()(Interface &Inter,Process &process) const {
        //std::cout << Inter.type << std::endl;
        //std::cout << Inter.comp << std::endl;
        if (Inter.type == Interface::type_ext) {
            if(Inter.comp == Interface::comp_deplacement_nul
            or Inter.comp == Interface::comp_deplacement 
            or Inter.comp == Interface::comp_vitesse_nulle 
            or Inter.comp == Interface::comp_vitesse) {
                compt_CL_vit(Inter,process.temps->pt);
            } else if (Inter.comp == Interface::comp_effort or Inter.comp == Interface::comp_effort_normal) {
                compt_CL_eff(Inter,process.temps->pt);
            /*} else if (Inter.comp=="effort_normal") {
                compt_CL_eff_normal(Inter,process.temps->pt);*/
            } else if (Inter.comp == Interface::comp_symetrie) {
                compt_CL_sym(Inter,process.temps->pt);
            } else if (Inter.comp == Interface::comp_vitesse_normale or Inter.comp == Interface::comp_deplacement_normal) {
                compt_CL_sym(Inter,process.temps->pt);
            } else if (Inter.comp == Interface::comp_periodique) {
                compt_parfait(Inter,*process.temps);
            } else if (Inter.comp == Interface::comp_cinetic_torseur) {
                compt_CL_torseur_cinetic(Inter, process.temps->pt, *process.CL);
            }
        } else if (Inter.type == Interface::type_int) {
//             if (Inter.comp == Interface::comp_periodique) {
//                 compt_parfait(Inter,*process.temps);
//             } else
            comportement_local_interface(Inter,process.temps->pt,process.temps->dt);
            /*if (Inter.comp==Interface::comp_parfait) {
                compt_parfait(Inter,process.temps->pt);
                //compt_jeu_impose(Inter,*process.temps);
            } else if (Inter.comp=="Contact" or Inter.comp=="Contact_jeu" or Inter.comp=="Contact_jeu_physique") {
                compt_contact(Inter,*process.temps);
            } else if (Inter.comp==Interface::comp_contact_parfait ) {
                //compt_contact_ep(Inter,*process.temps);
                compt_contact(Inter,*process.temps);
            } else if (Inter.comp==Interface::comp_cohesive) {
                compt_cohesif(Inter,*process.temps);
            } else if (Inter.comp=="Jeu_impose") {
                compt_jeu_impose(Inter,*process.temps);
            } else if (Inter.comp==Interface::comp_cassable_parfait) {
                compt_breakable(Inter,*process.temps);
            } else {
                std::cout<< "Erreur : comportement interface non reconnu" << std::endl;
                assert(0);
            }*/
        } else {
            std::cout<< "Erreur : comportement interface non reconnu" << std::endl;
            assert(0);
        }
    }
};


/** \ingroup   etape_locale
\brief Etape locale pour les sous-structures.
*/
struct etape_locale_sst {
    void operator()(Sst &S,Process &process) const {
        if      (S.f == S.pb.formulation_plasticity_isotropy_stat_Qstat) {calcul_plasticite(S,process);}
        else if (S.f == S.pb.formulation_elasticity_damageable_isotropy_stat_Qstat) {calcul_endommagement(S,process);}
        else if (S.f == S.pb.formulation_mesomodele) {calcul_mesomodele(S,process);}
    }
};


/** \ingroup  etape_locale
\brief Proc�dure principale pour l'�tape locale
*/
void etape_locale(PointedInterfaces &Inter,PointedSubstructures &S,Process &process) {
    apply(Inter,etape_locale_inter(),process);
    apply(S,etape_locale_sst(),process);
};

