#include "../../DEFINITIONS/TimeParameters.h"
//fonctions utilisees pour le probleme macro

/** \ingroup etape_lineaire
 \relates macroassemble
 \brief Cr�ation du second membre macro pour les interfaces de type sym�trie ou d�placement normal impos�.
 
Deux proc�dures sont �crites pour le 2d et pour le 3d (en changeant l'argument Interface<>). On s�lectionne les composantes macro diff�rentes de la r�sultante et moment selon la normale des efforts Fchap pour le pas de temps donn�.
*/
// Ajout des conditions aux limites pour le cas de CL symetriques dans le probleme macro
void add_bigF_CLsym(Vec<TYPEREEL> &bigF, Interface &inter, unsigned &data,int &imic) {
    Vec<unsigned> rep=range(inter.side[data].MeM.nb_cols());
    Vec<unsigned> repnimp;
#if DIM == 2
    if(find(rep,LMT::_1==(unsigned)0)){repnimp.push_back(0);}
    if(find(rep,LMT::_1==(unsigned)3)){repnimp.push_back(3);}
    //repnimp = {0,3}
#elif DIM == 3
    if(find(rep,LMT::_1==(unsigned)0)){repnimp.push_back(0);}
    if(find(rep,LMT::_1==(unsigned)1)){repnimp.push_back(1);}
    for(unsigned i=5;i<inter.side[data].MeM.nb_cols();i++)
         if(find(rep,LMT::_1==i)){repnimp.push_back(i);}
    //repnimp = {0,1,5,6,7,8}
#endif
    //{std::cout << "id : " << inter.id << std::endl;}
    //{std::cout << "repnimp : " << repnimp << std::endl;}
    //{std::cout << "inter.repddl : " << inter.repddl << std::endl;}
    Vec<unsigned> repF=inter.repddl[repnimp];
    bigF[repF] += trans(inter.side[data].MeM(range(inter.side[data].MeM.nb_rows()),repnimp) ) * inter.side[data].t[imic].Fchap;
}


/** \ingroup etape_lineaire
 \relates macroassemble
 \brief Assemblage du second membre du probleme macro 
 
 En bouclant sur les interfaces voisines de la sous-structure consid�r�e, on ajoute des quantit�s uniquement pour les interfaces autres qu'� d�placement impos�.
 Ayant rep�r� les ddls macro de l'interface consid�r�e,
 - On ajoute dans tous les cas le terme provenant de l'�tape micro 1 (Sst::Time::Fadd).
 - Si l'interface est de type effort, on ajoute les efforts donn�s (partie macro).
 - Si l'interface est de type sym�trie ou d�placement normal donn�, on utilise la fonction add_bigF_CLsym().
 */
struct macroassemble {
    void operator()(Sst &S,Vec<Interface> &Inter,TimeParameters &temps, MacroProblem &Global) const {
        int imic=temps.pt;   
        for(unsigned j=0;j<S.edge.size();++j) {
            unsigned q=S.edge[j].internum;
            unsigned data=S.edge[j].datanum;
            ///reperage de l'interface dans bigF
            Vec<unsigned> repF=Inter[q].repddl,repFadd=S.edge[j].repLE;
            if (Inter[q].type!="Ext" or (Inter[q].comp!="vit" and Inter[q].comp!="depl")) {
                /// participation des efforts de l'etape micro 1
                Global.bigF[repF] -= S.Fadd[repFadd];
                /// participation des efforts imposes
                if ( Inter[q].comp=="effort" or Inter[q].comp=="effort_normal") {
                    Global.bigF[repF] += trans(Inter[q].side[data].MeM)*Inter[q].side[data].t[imic].Fchap;
                } else if ( Inter[q].comp=="sym" or Inter[q].comp=="depl_normal" or Inter[q].comp=="vit_normale" ) {
                    add_bigF_CLsym(Global.bigF,Inter[q],data,imic);
                }
            }
        }
    }
};



/** \ingroup etape_lineaire
\relates interextr
\brief Extraction de WtildeM pour les interfaces de type sym�trie ou d�placement normal
 
Deux proc�dures sont �crites pour le 2d et pour le 3d. \f$ \tilde{W}^M \f$ pour le pas de temps donn� est impos� �gal � 0 puis on s�lectionne les composantes macro diff�rentes de la translation et rotation autour de la normale � l'interface.
*/
// Modification des multiplicateurs pour le cas de CL symetriques a la sortie du probleme macro
void modif_WtildeM_CLsym(Interface &Inter,Vec<TYPEREEL> &bigW, unsigned &data,int &imic) {
    Vec<unsigned> rep=range(Inter.side[data].MeM.nb_cols());
    Vec<unsigned> repnimp;
#if DIM == 2
    if(find(rep,LMT::_1==(unsigned)0)){repnimp.push_back(0);}
    if(find(rep,LMT::_1==(unsigned)1)){repnimp.push_back(1);}
    for(unsigned i=5;i<Inter.side[data].MeM.nb_cols();i++)
        if(find(rep,LMT::_1==i)){repnimp.push_back(i);}
    //repnimp = {0,1,5,6,7,8};
#elif DIM == 3
    if(find(rep,LMT::_1==(unsigned)0)){repnimp.push_back(0);}
    if(find(rep,LMT::_1==(unsigned)3)){repnimp.push_back(3);}
    //repnimp = {0,3};
#endif
    Vec<unsigned> repW=Inter.repddl[repnimp];
    Inter.side[data].t[imic].WtildeM.set(0.0);
    Inter.side[data].t[imic].WtildeM=Inter.side[data].eM(range(Inter.side[data].eM.nb_rows()),repnimp) * bigW[repW];
}


/** \ingroup etape_lineaire
\relates interextrmacro
 \brief Extraction des multiplicateurs a partir du vecteur obtenu dans le probleme macro
 
 En bouclant sur les sous-structures puis sur les bords de celles ci, on rep�re le cot� des interfaces concern�es et on extrait du vecteur Glob::bigW le multiplicateur selon le type d'interface :
 - Pour les interfaces � d�placement impos�, \f$ \tilde{W}^M \f$ est impos� �gal � 0.
 - Pour les interfaces de type sym�trie, on applique la proc�dure modif_WtildeM_CLsym()
 - Pour les autres interfaces on extrait les composantes macro concern�es � partir de Interface::repddl et on en d�duit la r�partition macro sur l'interface (par Interface::Side::eM)
 
 */
struct interextrmacro {
    void operator()(Sst &S,Vec<Interface> &Inter,TimeParameters &temps, MacroProblem &Global) const {
        int imic=temps.pt;
        for(unsigned j=0;j<S.edge.size();++j) {
            unsigned data=S.edge[j].datanum;
            unsigned q=S.edge[j].internum;
            if (Inter[q].comp == "depl" or Inter[q].comp == "vit") {
                Inter[q].side[data].t[imic].WtildeM.set(0.0);
            } else if (Inter[q].comp == "sym" or Inter[q].comp == "depl_normal" or Inter[q].comp == "vit_normale") {
                modif_WtildeM_CLsym(Inter[q],Global.bigW,data,imic);
            } else {
                Inter[q].side[data].t[imic].WtildeM=Inter[q].side[data].eM * Global.bigW[Inter[q].repddl];
            }
        }
    }
};

