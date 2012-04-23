#include "../../DEFINITIONS/Interface.h"
#include "../../DEFINITIONS/Sst.h"
#include "../../../LMT/include/containers/vec.h"
#include "../../../LMT/include/containers/vecpointedvalues.h"
#include "../../../LMT/include/containers/matcholamd.h"
using namespace LMT;

//**********************************************
// Optimisation des directions de recherche
//**********************************************

//Calcul les composantes normale et tangentielle du module d'elasticite
void calcul_En_Et(Sst &S, TYPEREEL &En, TYPEREEL &Et){
    //formulation isotrope ou viscoelastique 
    if (S.f->get_name()=="elasticity_isotropy_stat_Qstat" or S.f->get_name()=="elasticity_viscosity_Qstat" or S.f->get_name()=="plasticity_isotropy_stat_Qstat"){
        En = S.matprop.elastic_modulus;
        Et = En;
    }
    if(S.f->get_name()=="elasticity_orthotropy_stat_Qstat" or S.f->get_name()=="elasticity_orthotropy_damage_stat_Qstat" or S.f->get_name()=="mesomodele") {
        En = (S.matprop.elastic_modulus_1 + 2.*S.matprop.elastic_modulus_2)/3.;
        Et = En;
    }
}

// recherche des SST a prendre en compte pour la raideur du domaine complementaire
void find_SST_in_box(Vec<Sst> &S, Vec<TYPEREEL,DIM> &normale, Vec<TYPEREEL,DIM> &G, TYPEREEL &rayon, Vec<unsigned> &vois) {
    for(unsigned k=0;k<vois.size();++k) {
        unsigned ii=vois[k];
        for(unsigned j=0;j<S[ii].vois.size();++j) {
            int voisj=S[ii].vois[j];
            if (voisj!=-1){
                // verification si voisj est dans la boite
                if ( (length(S[(unsigned)voisj].G - G)<=rayon) and (dot(S[(unsigned)voisj].G - G, normale) >=0) and (find(vois,LMT::_1==(unsigned)voisj)==0 ) )
                    vois.push_back((unsigned)voisj);
            }
        }
    }
}

void modification_direction_CL(Interface &Inter, TYPEREEL &kn, TYPEREEL &kt, TYPEREEL &hn, TYPEREEL &ht) {
    TYPEREEL facteur = 1000.;
    if(Inter.type=="Ext") {
        if(Inter.comp=="depl" or Inter.comp=="depl_nul" or Inter.comp=="vit" or Inter.comp=="vit_nulle") {
            kn = kn * facteur;
            hn = hn / facteur;
            kt = kt * facteur;
            ht = ht / facteur;
        } else if(Inter.comp=="effort" or Inter.comp=="effort_normal") {
            kn = kn / facteur ;
            hn = hn * facteur ;
            kt = kt / facteur ;
            ht = ht * facteur ;
        } else if(Inter.comp=="sym" or Inter.comp=="depl_normal" or Inter.comp=="vit_normale") {
            kn = kn * facteur ;
            hn = hn / facteur ;
            kt = kt / facteur ;
            ht = ht * facteur ;
        } else if(Inter.comp=="periodique") {
            //ben on fait rien :)
        } else {
            std::cout << "optimisation direction : Type de condition limite non pris en compte : " << Inter.comp << endl;
            assert(0);
        }
    } else if(Inter.type=="Int") {
        if(Inter.comp!="Parfait" and Inter.comp!="Jeu_impose" and Inter.comp!="Cohesive" and Inter.comp!="Contact_ep") {
            TYPEREEL facteur_frottement;
            TYPEREEL eps=1e-6;
            if(Inter.param_comp->coeffrottement<=eps){facteur_frottement=1e-3;}
            else{facteur_frottement=Inter.param_comp->coeffrottement;}
            //facteur_frottement=1;
            kt = kt * facteur_frottement;
            ht = ht / facteur_frottement;
        }
    }
}


/** \ingroup  Interfaces
\brief Application des directions de recherche par coté d'interface
 
   Plusieurs choix de directions sont possibles selon le paramètre LatinParameters::ktype.
   Dans tous les cas une matrice globale sur l'interface Interface::Side::kglob et son inverse 
 Interface::Side::hglob est déterminée à partir des directions normale et tangentielle. 
   La matrice est constituée de bloc 3x3 ou 2x2 du type 
      \f$ k_n* \vec{n} x \vec{n} + k_t * (I -\vec{n} x \vec{n})  \f$ 
 où \f$ k_n \f$ et \f$ k_t \f$ sont déterminés automatiquement pour l'interface considérée et 
 modifiés selon les conditions aux limites. Ces matrices globales sur l'interface sont utilisées 
 pour la construction de la matrice de rigidité et le calcul des seconds membres dans les étapes 
 micro1 et micro2 ... Par contre pour l'étape locale, on effectue les opérations point par point 
 et on utilise donc la matrice locale. Ceci permet de travailler directement avec le vecteur 
 concaténé des quantités sur l'interface en chaque point.
*/
struct optimise_direction {
    void operator()(Interface &inter, Vec<VecPointedValues<Sst> > &S, LatinParameters &latin) const {
        typedef Interface::TMATS TMAT;
        for(unsigned q=0;q<inter.side.size();++q) {
            // recherche des ssts complementaires
            unsigned ii=0,jj=0;
            if (inter.side.size()==1) {
                // seule la sst voisine est prise en compte pour les interfaces de bord (plus de ssts n'est pas necessaire car la direction est deja approchee ...)
                ii=inter.side[0].vois[0];
                jj=inter.side[0].vois[1];
            } else {
                //on determine la sous-structure complementaire de l'interface consideree
                ii=inter.side[!q].vois[0];
                jj=inter.side[!q].vois[1];
            }

            TYPEREEL L=0, En=0, Et=0;
            //calcul de la longueur de l'interface à prendre en compte pour la determination automatique du scalaire
#if DIM==2
                L = inter.measure;
#elif DIM==3
                L= std::sqrt(inter.measure);
#endif
            //determination du module d'young correspondant
            //fonction generee a partir du __init__.py
            calcul_En_Et(S[ii],En,Et);
            //calcul des directions de recherche normale et tangentielle
            TYPEREEL kn=0, kt=0, hn=0, ht=0;


            TMAT k, h;
            k.resize(DIM);
            h.resize(DIM);

            TMAT Id;
            Id.resize(DIM);
            Id.diag()=Vec<TYPEREEL,DIM>(1);

            if(latin.ktype=="scalaire_donne") {
                kn = latin.kfact;
                kt = latin.kfact;
                hn = 1./kn;
                ht = 1./kt;
            } else if (latin.ktype=="scalaire_donne_CL") {
                kn = latin.kfact;
                kt = latin.kfact;
                hn = 1./kn;
                ht = 1./kt;
                modification_direction_CL(inter, kn, kt, hn, ht);
            } else if (latin.ktype=="scalaire_auto") {
                kn = En/L * latin.kfact;
                kt = Et/L * latin.kfact;
                hn = 1./kn;
                ht = 1./kt;
            } else if (latin.ktype=="scalaire_auto_CL") {
                kn = En/L * latin.kfact;
                kt = Et/L * latin.kfact;
                hn = 1./kn;
                ht = 1./kt;
                modification_direction_CL(inter, kn, kt, hn, ht);
            } else if (latin.ktype=="scalaire_eigenvalue_CL") {
                std::cout << "Pas encore bien implemente" << endl;
                assert(0);
            } else {
                std::cout << "type de direction de recherche non pris en compte " << endl;
                assert(0);
            }


            //creation de la matrice kglo et hglo et assignation des direction h, k , hn, ht, kn, kt
            inter.side[q].hn = hn;
            inter.side[q].ht = ht;
            inter.side[q].kn = kn;
            inter.side[q].kt = kt;
            inter.side[q].kglo.resize(inter.side[q].nodeeq.size()*DIM);
            inter.side[q].hglo.resize(inter.side[q].nodeeq.size()*DIM);

            for(unsigned i=0;i<inter.side[q].nodeeq.size();i++) {
                Vec<unsigned> rep=range(i*DIM,(i+1)*DIM);
                Vec<TYPEREEL,DIM> n=inter.side[q].neq[rep];
                TMAT nn;
                tens(n,n,nn);

                inter.side[q].kglo(rep,rep) = kn * nn + kt * (Id - nn);
                inter.side[q].hglo(rep,rep) = hn * nn + ht * (Id - nn);
            }
        }

        //copie des directions de recherche
        if(latin.copydirection==1) {
            TYPEREEL kn=0,kt=0,hn=0,ht=0;
            TMAT kglo, hglo;
            kglo.resize(inter.side[0].nodeeq.size()*DIM);
            hglo.resize(inter.side[0].nodeeq.size()*DIM);
            for(unsigned q=0;q<inter.side.size();q++) {
                kn+=inter.side[q].kn;
                kt+=inter.side[q].kt;
                hn+=inter.side[q].hn;
                ht+=inter.side[q].ht;
                kglo+=inter.side[q].kglo;
                hglo+=inter.side[q].hglo;
            }
            kn/=inter.side.size();
            kt/=inter.side.size();
            hn/=inter.side.size();
            ht/=inter.side.size();
            kglo=kglo/inter.side.size();
            hglo=hglo/inter.side.size();
            for(unsigned q=0;q<inter.side.size();q++) {
                inter.side[q].kn= kn ;
                inter.side[q].kt= kt ;
                inter.side[q].hn= hn ;
                inter.side[q].ht= ht ;
                inter.side[q].kglo = kglo ;
                inter.side[q].hglo = hglo ;
            }
        }
    }
};



