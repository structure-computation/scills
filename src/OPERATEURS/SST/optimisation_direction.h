#include "containers/matcholamd.h"
using namespace LMT;
using namespace std;

//**********************************************
// Optimisation des directions de recherche
//**********************************************

// recherche des SST a prendre en compte pour la raideur du domaine complementaire
template<class TV1,class TV,class T>
void find_SST_in_box(TV1 &S, TV &normale, TV &G, T &rayon, Vec<unsigned> &vois) {
    for(unsigned k=0;k<vois.size();++k) {
        unsigned ii=vois[k];
        for(unsigned j=0;j<S[ii].vois.size();++j) {
            int voisj=S[ii].vois[j];
            //std::cout << "voisj " <<voisj << endl;
            if (voisj!=-1) {
                // verification si voisj est dans la boite
                //std::cout << "length " << (length(S[(unsigned)voisj].G - G)<=rayon) << "norm " << (dot(S[(unsigned)voisj].G - G, normale) >=0) << "invois" << (find(vois,LMT::_1==(unsigned)voisj)==0 )<< endl;
                if ( (length(S[(unsigned)voisj].G - G)<=rayon) and (dot(S[(unsigned)voisj].G - G, normale) >=0) and (find(vois,LMT::_1==(unsigned)voisj)==0 ) )
                    vois.push_back((unsigned)voisj);
            }
        }
    }
}

// // creation de la raideur et maillage du complementaire
// template<class INTER,class TV1, class TF, class TM>
// void create_mesh_supp(INTER &inter,TV1 &S, unsigned ii, unsigned jj, TF &f, TM &m) {
//     typedef typename INTER::T TT;
//     typedef Vec<TT,INTER::dim> TV;
//     // creation de la boite semi circulaire dans laquelle rechercher les SSTS a assembler
//     TV box1,box2;
//     box1=inter.box[0];
//     box2=inter.box[1];
//     TT rayon = length(box1-box2);
//     Vec<TT, INTER::dim> centre = inter.G;
//     Vec<TT, INTER::dim> normale = (S[ii].G-inter.G);
//     normale=normale/norm_2(normale);
// 
//     // recherche si les cdg des ssts sont dans la zone
//     Vec<unsigned> vois(ii);
//     unsigned sizevois;
//     //std::cout << "ii " << ii << "ray " << rayon  << endl;
//     do {
//         sizevois=vois.size();
//         find_SST_in_box(S, normale,centre, rayon, vois);
//     } while(sizevois!=vois.size());
// 
//     //Creation du maillage globales et matrice de raideur
//     if (vois==Vec<TT>(ii)) {
//         m.append(S[ii].mesh);
//         f.allocate_matrices();
//         f.matrices(Number<0>())=S[ii].f.matrices(Number<0>());
//     } else {
//         for(unsigned i=0;i<vois.size();++i) {
//             m.append(S[vois[i]].mesh);
//             remove_doubles(m,1e-6);
//         }
//         std::cout << m.node_list.size() << endl;
//         //remove_doubles(m,1e-2);
//         //affichmesh(m);
//         //calcul raideur
//         f.allocate_matrices();
//         f.assemble();
//         std::cout << "K0 assemblee " << endl;
//         //    display_structure(f.matrices(Number<0>()));
//     }
// }

// // resolution du probleme au valeur propre pour le cas d'un parametre scalaire (interface parfaite)
// template<class TMAT, class TV, class INTER> void solve_eig_power_Pm(TMAT &K, TV &repddl, INTER &interside )
// {
//   typedef typename TMAT::T T;
//   T res;
//   T eps=1e-6;
//   T normX1X0=1;
//   //creation second membre
//   Vec<T> F,X0,X1,U,X2;
//   //X1=Vec<T>::random(interside.Pm.nb_rows());
//   X1.resize(interside.Pm.nb_rows(),1);
//   X1=X1/norm_2(X1);
//   F.resize(K.nb_rows() );
//   K.get_factorization();
//   //Inv<T, Sym<>, SparseLine<> > invK;
//   //get_factorization(K,invK.fact);
//   //Inv<T, Sym<>, SparseLine<> > invK=inv(K);
//   Mat<T,Gen<>, SparseLine<> > NtMP,PN;
//   NtMP=interside.Nt*(interside.M*interside.Pm);
//   PN=interside.Pm*interside.N;
//   unsigned it=0;
//   while (normX1X0 >eps && it<10000)
//   {
//     X0=X1;
//     F.set(0.0);
//     F[repddl]=NtMP*X0;
//     U=K.solve(F);
//     X1=PN*U[repddl];
//     X1=X1/norm_2(X1);
//     X2=X1-X0;
//     normX1X0=norm_2(X2);
//     it+=1;
//   }
//   X0=X1;
//   //std::cout << "X0 " <<X0 << endl;
//   F.set(0.0);
//   F[repddl]=NtMP*X0;
//   //std::cout << "F " << F << endl;
//   X1=PN*U[repddl];
//   res=dot(X1,X0);
//
//   // creation de la matrice kloc et cloc
//   interside.cloc.resize(interside.Pm.nb_rows());
//   interside.kloc.resize(interside.Pm.nb_rows());
//   for(unsigned i=0;i<interside.Pm.nb_rows();++i)
//   {
//     interside.cloc(i,i)=res;
//     interside.kloc(i,i)=1.0/res;
//   }
//   //cn et ct valent aussi res
//   interside.hn=res;
//   interside.ht=res;
//   interside.kn=1.0/res;
//   interside.kt=1.0/res;
//   interside.h=res;
//   interside.k=1.0/res;
//
// }
//
// // resolution du probleme au valeur propre pour le cas d'un parametre scalaire selon la normale et un pour la partie tangentielle (interface contact)
// template<class TMAT, class TV, class INTERSIDE> void solve_eig_power_Pnt(TMAT &K, TV &repddl, INTERSIDE &interside)
// {
//   typedef typename TMAT::T T;
//   T cn,ct,c;
//   T eps=1e-6;
//   T normX1X0=1;
//   Vec<T> F,X0,X1,U,X2;
//   //vecteur initial
//   X1=Vec<T>::random(interside.Pm.nb_rows());
//   X1=X1/norm_2(X1);
//   F.resize( K.nb_rows() );
//   //Inv<T, Sym<>, SparseLine<> > invK=inv(K);
//   Mat<T,Gen<>, SparseLine<> > NtMP,PN;
//   NtMP=interside.Nt*interside.M*interside.Pm*interside.Pn;
//   PN=interside.Pn*interside.Pm*interside.N;
//
//   // recherche valeur propre pour l'operateur normal
//   unsigned it=0;
//   while (normX1X0 >eps && it<10000)
//   {
//     X0=X1;
//     F.set(0.0);
//     F[repddl]=NtMP*X0;
//     U=invK*F;
//     X1=PN*U[repddl];
//     X1=X1/norm_2(X1);
//     X2=X1-X0;
//     normX1X0=norm_2(X2);
//     it+=1;
//   }
//   X0=X1;
//   F.set(0.0);
//   F[repddl]=NtMP*X0;
//   X1=PN*U[repddl];
//   cn=dot(X1,X0);
//
//   // recherche valeur propre pour l'operateur tangentiel
//   it=0;
//   normX1X0=1;
//   X1=Vec<T>::random(interside.Pm.nb_rows());
//   X1=X1/norm_2(X1);
//   NtMP=interside.Nt*interside.M*interside.Pm*interside.Pt;
//   PN=interside.Pt*interside.Pm*interside.N;
//   while (normX1X0 >eps && it<10000)
//   {
//     X0=X1;
//     F.set(0.0);
//     F[repddl]=NtMP*X0;
//     U=invK*F;
//     X1=PN*U[repddl];
//     X1=X1/norm_2(X1);
//     X2=X1-X0;
//     normX1X0=norm_2(X2);
//     it+=1;
//   }
//   X0=X1;
//   F.set(0.0);
//   F[repddl]=NtMP*X0;
//   X1=PN*U[repddl];
//   ct=dot(X1,X0);
//
//   interside.hn=cn;
//   interside.kn=1/cn;
//   interside.ht=ct;
//   interside.kt=1/ct;
//
//   //recherche pour les interfaces parfaites
//   it=0;
//   normX1X0=1;
//   X1=Vec<T>::random(interside.Pm.nb_rows());
//   X1=X1/norm_2(X1);
//   NtMP=interside.Nt*interside.M*interside.Pm;
//   PN=interside.Pm*interside.N;
//   while (normX1X0 >eps && it<10000)
//   {
//     X0=X1;
//     F.set(0.0);
//     F[repddl]=NtMP*X0;
//     U=invK*F;
//     X1=PN*U[repddl];
//     X1=X1/norm_2(X1);
//     X2=X1-X0;
//     normX1X0=norm_2(X2);
//     it+=1;
//   }
//   X0=X1;
//   F.set(0.0);
//   F[repddl]=NtMP*X0;
//   X1=PN*U[repddl];
//   c=dot(X1,X0);
//
//   // creation de la matrice kloc et cloc
//   interside.h=c;
//   interside.k=1.0/c;
//   interside.cloc.resize(interside.Pm.nb_rows());
//   interside.kloc.resize(interside.Pm.nb_rows());
//   for(unsigned i=0;i<interside.Pm.nb_rows();++i)
//   {
//     interside.cloc(i,i)=c;
//     interside.kloc(i,i)=1.0/c;
//   }
//
//
// }
//


template<class INTER, class T>
void modification_direction_CL(INTER &Inter, T &kn, T &kt, T &hn, T &ht) {
    T facteur = 1000.;
//     std::cout << Inter.type << std::endl;
//     std::cout << Inter.comp << std::endl;
    if(Inter.type=="Ext") {
        if(Inter.comp=="depl" or Inter.comp=="depl_nul") {
            kn = kn * facteur;
            hn = hn / facteur;
            kt = kt * facteur;
            ht = ht / facteur;
        } else if(Inter.comp=="effort" or Inter.comp=="effort_normal") {
            kn = kn / facteur ;
            hn = hn * facteur ;
            kt = kt / facteur ;
            ht = ht * facteur ;
        } else if(Inter.comp=="sym" or Inter.comp=="depl_normal") {
            kn = kn * facteur ;
            hn = hn / facteur ;
            kt = kt / facteur ;
            ht = ht * facteur ;
        } else if(Inter.comp=="periodique") {
            //ben on fait rien :)
        } else {
            std::cout << "optimisation direction : Type de condition limite non pris en compte" << endl;
            assert(0);
        }
    } else if(Inter.type=="Int") {
        if(Inter.comp!="Parfait" and Inter.comp!="Jeu_impose" and Inter.comp!="Cohesive" and Inter.comp!="Contact_ep") {
            T facteur_frottement;
            T eps=1e-6;
            if(Inter.param_comp->coeffrottement<=eps){facteur_frottement=1e-3;}
            else{facteur_frottement=Inter.param_comp->coeffrottement;}
            //facteur_frottement=1;
            kt = kt * facteur_frottement;
            ht = ht / facteur_frottement;
        }
    }
}

#include "calcul_En_Et.h"

/** \ingroup  Interfaces
\brief Application des directions de recherche par coté d'interface
 
Plusieurs choix de directions sont possibles selon le paramètre Param::ktype.
Dans tous les cas une matrice globale sur l'interface Interface::Side::kglob et son inverse Interface::Side::hglob est déterminée à partir des directions normale et tangentielle. La matrice est constituée de bloc 3x3 ou 2x2 du type \f$ k_n* \vec{n} x \vec{n} + k_t * (I -\vec{n} x \vec{n})  \f$ où \f$ k_n \f$ et \f$ k_t \f$ sont déterminés automatiquement pour l'interface considérée et modifiés selon les conditions aux limites. Ces matrices globales sur l'interface sont utilisées pour la construction de la matrice de rigidité et le calcul des seconds membres dans les étapes micro1 et micro2 ... Par contre pour l'étape locale, on effectue les opérations point par point et on utilise donc la matrice locale. Ceci permet de travailler directement avec le vecteur concaténé des quantités sur l'interface en chaque point.
*/
struct optimise_direction {
    template<class TV1, class INTER>
    void operator()(INTER &inter, TV1 &S, LATIN &process) const {
        typedef typename INTER::TMATS TMAT;
        typedef typename INTER::T T;
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

            T L=0, En=0, Et=0;
            //calcul de la longueur de l'interface à prendre en compte pour la determination automatique du scalaire
            if (INTER::dim==2)
                L = inter.measure;
            else if (INTER::dim==3)
                L= std::sqrt(inter.measure);

            //determination du module d'young correspondant
            //fonction generee a partir du __init__.py
            calcul_En_Et(S[ii],En,Et);
            //calcul des directions de recherche normale et tangentielle
            T kn=0, kt=0, hn=0, ht=0;


            TMAT k, h;
            k.resize(INTER::dim);
            h.resize(INTER::dim);

            TMAT Id;
            Id.resize(INTER::dim);
            Id.diag()=Vec<T,INTER::dim>(1);

            if(process.ktype=="scalaire_donne") {
                kn = process.kfact;
                kt = process.kfact;
                hn = 1./kn;
                ht = 1./kt;
            } else if (process.ktype=="scalaire_donne_CL") {
                kn = process.kfact;
                kt = process.kfact;
                hn = 1./kn;
                ht = 1./kt;
                modification_direction_CL(inter, kn, kt, hn, ht);
            } else if (process.ktype=="scalaire_auto") {
                kn = En/L * process.kfact;
                kt = Et/L * process.kfact;
                hn = 1./kn;
                ht = 1./kt;
            } else if (process.ktype=="scalaire_auto_CL") {
                kn = En/L * process.kfact;
                kt = Et/L * process.kfact;
                hn = 1./kn;
                ht = 1./kt;
                modification_direction_CL(inter, kn, kt, hn, ht);
            } else if (process.ktype=="scalaire_eigenvalue_CL") {
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
            inter.side[q].kglo.resize(inter.side[q].nodeeq.size()*INTER::dim);
            inter.side[q].hglo.resize(inter.side[q].nodeeq.size()*INTER::dim);

            for(unsigned i=0;i<inter.side[q].nodeeq.size();i++) {
                Vec<unsigned> rep=range(i*INTER::dim,(i+1)*INTER::dim);
                Vec<T,INTER::dim> n=inter.side[q].neq[rep];
                TMAT nn;
                tens(n,n,nn);

                inter.side[q].kglo(rep,rep) = kn * nn + kt * (Id - nn);
                inter.side[q].hglo(rep,rep) = hn * nn + ht * (Id - nn);
            }
        }

        //copie des directions de recherche
        if(process.copydirection==1) {
            T kn=0,kt=0,hn=0,ht=0;
            TMAT kglo, hglo;
            kglo.resize(inter.side[0].nodeeq.size()*INTER::dim);
            hglo.resize(inter.side[0].nodeeq.size()*INTER::dim);
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



