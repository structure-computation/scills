#ifndef MESH_DATA_ACCESSORS_SST_H
#define MESH_DATA_ACCESSORS_SST_H

/** \brief Accesseurs (et mutateur) pour les proprietes du maillage LMT::Mesh des Sst
 * 
 * 
 * Ce header regroupe l'ensemble des fonctions associes a la manipulation des variables d'elements de maillage definies par les fichiers de formulation.
 * Pour chaque 'Variable' xxx du fichiers de formulation sont definis :
 *   1) des fonctions pour prendre les valeurs de xxx dans le maillage de 'S' et les stocker dans un stockage externe:
 *     a) 'download_xxx(Sst &S, Vec< TypeXxx > &sto)' on stocke dans le vecteur 'sto' (attention a initialiser sa taille) - TypeXxx etant le type des valeurs de xxx
 *     b) 'download_xxx(Sst &S, Sst::Time &t)' prend les valeurs de xxx dans le maillage de 'S' et le stocke dans le vecteur correspondant de 't'
 *     c) 'download_xxx(Sst &S, unsigned pt)' prend les valeurs de xxx dans le maillage de 'S' et le stocke dans le vecteur correspondant de 'S.t[pt]' - 'pt' representant l'indice du pas de temps courant
 *   2) des fonctions pour assigner les valeurs de xxx dans le maillage de 'S' a celles stocker dans un stockage externe sous la forme 'upload_xxx(...)'
 *   3) un foncteur '__Download_xxx', utilise par les fonctions 'download_xxx' pour un apply sur les elements
 *   4) un foncteur '__Upload_xxx', equivalent de '__Download_xxx'
 * 
 * Un template est disponible a la fin du fichier (en commentaire). Il suffit de remplacer les '***' par les noms correspondant
 * 
 * Ces fonctions sont actuellement definis pour:
 *   - sigma (vecteur en notation de Voigt)
 *   - sigma_von_mises (scalaire) PAS UTILISEE
 *   - epsilon (vecteur en notation de Voigt)
 *   - epsilon_e (vecteur en notation de Voigt) PAS UTILISEE
 *   - epsilon_p (vecteur en notation de Voigt)
 *   - p (scalaire)
 *   - R_p (scalaire)
 *   - d1 (scalaire)
 *   - d2 (scalaire)
 *   - df (scalaire)
 *   - Yd (vecteur)
 * 
 * Remarque: Les fonctions 'upload_q' et 'upload_f_vol_e' - utilisant la meme nomenclature - ont egalement ete definies pour le deplacement et les efforts volumiques, par convenience.
 * 
 * Attention !!!
 * - Ne pas definir les surcharges utilisant 'pt' ou 't' si la variable n'a pas de vecteur de stockage dans la definition de 'Sst::Time'
 * - Ne pas appliquer les fonctions d'une variable a un element sur laquelle elle n'est pas definie, sous peine d'erreur a la compilation!
 */



#include "../DEFINITIONS/main_typedef.h"
#include "../DEFINITIONS/Sst.h"


//*
///----------------------------   Q   ----------------------------------
void upload_q(Sst &S,Vector &sto);
void upload_q(Sst &S,Sst::Time &t);
void upload_q(Sst &S,unsigned pt);
//*/
//*
///---------------------------- SIGMA ----------------------------------
/// Recuperation depuis le maillage
void download_sigma(Sst &S,Vec<VoigtVector> &sto);

struct __Download_sigma{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        sto[i_elem] = e.sigma[0];
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_sigma(Sst &S,Vec<VoigtVector> &sto);

struct __Upload_sigma{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        e.sigma[0] = sto[i_elem];
        i_elem++;
    }
};
//*/
/*
///---------------------------- SIGMA VON MISES ----------------------------
/// Recuperation depuis le maillage
void download_sigma_von_mises(Sst &S,Vector &sto);

struct __Download_sigma_von_mises{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        sto[i_elem] = e.sigma_von_mises;
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_sigma_von_mises(Sst &S,Vector &sto);

struct __Upload_sigma_von_mises{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        e.sigma_von_mises = sto[i_elem];
        i_elem++;
    }
};
//*/
//*
///---------------------------- EPSILON ---------------------------------
/// Recuperation depuis le maillage
void download_epsilon(Sst &S,Vec<VoigtVector> &sto);

struct __Download_epsilon{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        sto[i_elem] = e.epsilon[0];
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_epsilon(Sst &S,Vec<VoigtVector> &sto);

struct __Upload_epsilon{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        e.epsilon[0] = sto[i_elem];
        i_elem++;
    }
};
//*/
/*
///---------------------------- EPSILON E ---------------------------------
/// Recuperation depuis le maillage
void download_epsilon_e(Sst &S,Vec<VoigtVector> &sto);

struct __Download_epsilon_e{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        sto[i_elem] = e.epsilon_e[0];
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_epsilon_e(Sst &S,Vec<VoigtVector> &sto);

struct __Upload_epsilon_e{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        e.epsilon_e[0] = sto[i_elem];
        i_elem++;
    }
};
*/

///---------------------------- EPSILON P ---------------------------------
/// Recuperation depuis le maillage
void download_epsilon_p(Sst &S,Vec<VoigtVector> &sto);
void download_epsilon_p(Sst &S,Sst::Time &t);
void download_epsilon_p(Sst &S,unsigned pt);

struct __Download_epsilon_p{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        sto[i_elem] = e.epsilon_p[0];
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_epsilon_p(Sst &S,Vec<VoigtVector> &sto);
void upload_epsilon_p(Sst &S,Sst::Time &t);
void upload_epsilon_p(Sst &S,unsigned pt);

struct __Upload_epsilon_p{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        e.epsilon_p[0] = sto[i_elem];
        i_elem++;
    }
};

//*
///---------------------------     P     ----------------------------------
/// Recuperation depuis le maillage
void download_p(Sst &S,Vector &sto);
void download_p(Sst &S,Sst::Time &t);
void download_p(Sst &S,unsigned pt);

struct __Download_p{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        sto[i_elem] = e.plast_cumulee;
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_p(Sst &S,Vector &sto);
void upload_p(Sst &S,Sst::Time &t);
void upload_p(Sst &S,unsigned pt);

struct __Upload_p{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        e.plast_cumulee = sto[i_elem];
        i_elem++;
    }
};
//*/
//*
///-------------------------------  R_p  -------------------------------------------
/// Recuperation depuis le maillage
void download_R_p(Sst &S,Vector &sto);
void download_R_p(Sst &S,Sst::Time &t);
void download_R_p(Sst &S,unsigned pt);

struct __Download_R_p{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        sto[i_elem] = e.plast_ecrouissage;
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_R_p(Sst &S,Vector &sto);
void upload_R_p(Sst &S,Sst::Time &t);
void upload_R_p(Sst &S,unsigned pt);

struct __Upload_R_p{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        e.plast_ecrouissage = sto[i_elem];
        i_elem++;
    }
};
//*/
//*
///-------------------------------   d1   -------------------------------------------
/// Recuperation depuis le maillage
void download_d1(Sst &S,Vector &sto);
void download_d1(Sst &S,Sst::Time &t);
void download_d1(Sst &S,unsigned pt);

struct __Download_d1{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        sto[i_elem] = e.d1;
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_d1(Sst &S,Vector &sto);
void upload_d1(Sst &S,Sst::Time &t);
void upload_d1(Sst &S,unsigned pt);

struct __Upload_d1{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        e.d1 = sto[i_elem];
        i_elem++;
    }
};
//*/
//*
///-------------------------------   d2   -------------------------------------------
/// Recuperation depuis le maillage
void download_d2(Sst &S,Vector &sto);
void download_d2(Sst &S,Sst::Time &t);
void download_d2(Sst &S,unsigned pt);

struct __Download_d2{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        sto[i_elem] = e.d2;
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_d2(Sst &S,Vector &sto);
void upload_d2(Sst &S,Sst::Time &t);
void upload_d2(Sst &S,unsigned pt);

struct __Upload_d2{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        e.d2 = sto[i_elem];
        i_elem++;
    }
};
//*/
//*
///-------------------------------   df   -------------------------------------------
/// Recuperation depuis le maillage
void download_df(Sst &S,Vector &sto);
void download_df(Sst &S,Sst::Time &t);
void download_df(Sst &S,unsigned pt);

struct __Download_df{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        sto[i_elem] = e.df;
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_df(Sst &S,Vector &sto);
void upload_df(Sst &S,Sst::Time &t);
void upload_df(Sst &S,unsigned pt);

struct __Upload_df{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        e.df = sto[i_elem];
        i_elem++;
    }
};
//*/
//*
///-------------------------------   Yd   ------------------------------------------
/// Recuperation depuis le maillage
void download_Yd(Sst &S,Vec<Vec<TYPEREEL,3> > &sto);
void download_Yd(Sst &S,Sst::Time &t);
void download_Yd(Sst &S,unsigned pt);

struct __Download_Yd{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        sto[i_elem] = e.Yd;
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_Yd(Sst &S,Vec<Vec<TYPEREEL,3> > &sto);
void upload_Yd(Sst &S,Sst::Time &t);
void upload_Yd(Sst &S,unsigned pt);

struct __Upload_Yd{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        e.Yd = sto[i_elem];
        i_elem++;
    }
};
//*/

/*    TEMPLATE POUR LES FUTURES VARIABLES : 
 * Remplacer les '*xxx* par le nom de la variable et les *TypeXxx* par son type 
///-------------------------------------------------------------------------------
/// Recuperation depuis le maillage
void download_*xxx*(Sst &S,Vec<*TypeXxx*> &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_*xxx*(),sto,i_elem);
}

void download_*xxx*(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_*xxx*(),t.*xxx*,i_elem);
}

void download_*xxx*(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_*xxx*(),S.t[pt].*xxx*,i_elem);
}

struct __Download_*xxx*{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        sto[i_elem] = e.*xxx*;
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_*xxx*(Sst &S,Vec<*TypeXxx*> &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_*xxx*(),sto,i_elem);
}

void upload_*xxx*(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_*xxx*(),t.*xxx*,i_elem);
}

void upload_*xxx*(Sst &S,unsigned pt){
    apply(S.mesh->elem_list,__Upload_*xxx*(),S.t[pt].*xxx*,i_elem);
}

struct __Upload_*xxx*{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        e.*xxx* = sto[i_elem];
        i_elem++;
    }
};
*/


#endif //ELEMENTS_VARIABLES_ACCESSORS_H
