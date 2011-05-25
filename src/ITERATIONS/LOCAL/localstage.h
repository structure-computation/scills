#include "comportements_interfaces.h"

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
- interface de type deplacement (Interface::comp == "depl") : compt_CL_depl()
- interface de type effort (Interface::comp == "effort") : compt_CL_eff()
- interface de type sym�trie (Interface::comp == "sym") : compt_CL_sym()
- interface de type d�placement normal donn�  (Interface::comp == "depl_normal") : compt_CL_sym()
- interface de type parfait (Interface::comp == "Parfait") : compt_parfait()
- interface de type contact avec ou sans frottement (Interface::comp == "Contact") : compt_contact()
- interface de type jeu contact avec ou sans frottement (Interface::comp == "Contact_jeu" ou "Contact_physique") : compt_contact()
- interface de type epaisse contact avec ou sans frottement (Interface::comp == "Contact_ep" ) : compt_contact()
- interface de type coh�sive (Interface::comp == "Cohesive") : compt_cohesif()
 
Il suffit donc de rajouter un comportement dans cette proc�dure et programmer la fonction correspondante pour ajouter un comportement d'interface.
*/
struct etape_locale_inter {
    template<class INTER>
    void operator()(INTER &Inter,Param &process) const {
//         std::cout << Inter.type << std::endl;
//         std::cout << Inter.comp << std::endl;
        if (Inter.type=="Ext") {
            if(Inter.comp=="depl") {
                compt_CL_depl(Inter,process.temps->pt);
            } else if (Inter.comp=="effort") {
                compt_CL_eff(Inter,process.temps->pt);
            } else if (Inter.comp=="sym") {
                compt_CL_sym(Inter,process.temps->pt);
            } else if (Inter.comp=="depl_normal") {
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

#include "comportements_sst.h"
/** \ingroup   etape_locale
\brief Etape locale pour les sous-structures.
*/
struct etape_locale_sst {
    template<class SST>
    void operator()(SST &S) const {
#ifdef FORMUENDO
        if (S.f->get_name()=="elasticity_orthotropy_damage_stat_Qstat")
            compt_sst_damage(S);
#endif

    }
};


#include "containers/evaluate_nb_cycles.h"
/** \ingroup  etape_locale
\brief Proc�dure principale pour l'�tape locale
*/
template<class TV2,class TV1>
void etape_locale(TV2 &Inter,TV1 &S,Param &process) {
    apply(Inter,etape_locale_inter(),process);
//     apply(S,etape_locale_sst());
};

