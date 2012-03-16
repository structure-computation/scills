#ifndef CREATE_MESURE_G_NEQ_INTER
#define CREATE_MESURE_G_NEQ_INTER

#include "algebre.h"
#include "utilitaires.h"

//********************************************
// calcul mesure cdg pour les interfaces
//*******************************************
/** \ingroup Operateurs_inter
\brief  Calcul du centre de gravité de l'interface et de la mesure de celle ci (Interface::G)
*/
struct CalcMeasure_G {
   template<class INTER > void operator()(INTER &Inter) const {
      Inter.G=barycenter_constant_rho(*Inter.side[0].mesh);
      Inter.measure=measure(*Inter.side[0].mesh);
   }
};

//*************************************************
// calcul des normales aux elements de l'interface :
// normales toujours de l'interieure de la Sst vers le cote correspondant
//*************************************************
struct Calculate_Normales {
   template<class TE, class TM,class TV, class SIDE> void operator() (TE &e, TM &m, TV &box, SIDE &side) const {
      double eps=1e-6;
      //verification si l'element est dans la boite
      typedef typename TM::TNode::Pvec Pvec;
      Pvec Gs = center(e);
      //test si l'element est dans la boite
      //if( bool_vec( Gs <= box[1]+eps and  Gs >= box[0]-eps) pt_in_box(center(*m2.elem_list[i]),box1)==1) {
      if( pt_in_box(Gs,box,eps)==1 ) {
         // on recherche donc l'element du maillage d'interface correspondant
         for(unsigned i=0;i<side.mesh->elem_list.size();i++) {
            Pvec Gi = center(*side.mesh->elem_list[i]);
            if( (length( Gi - Gs)<=eps) ) {
               //on calcule alors la normale exterieure a la structure
               Pvec Ge = center(*m.skin.get_parents_of(e)[0]);
               Pvec normext = Gs-Ge;
               //on determine la normale a l'element e
               Pvec n = sample_normal(e);
               //on adapte alors le signe de la normale selon la normale exterieure a la structure
               n = sign(dot(n,normext))*n;
               typename TE::T normn=norm_2(n);
               n/=normn;
               //on assigne alors la normale trouvee au vecteur neq
               side.neq[range(i*TM::dim,(i+1)*TM::dim)]=n;
               break;
            }
         }
      }
   }
};
/** \ingroup Operateurs_inter
\brief Création des normales à l'interface. (Interface::Side::neq)
 
Les normales sont créées pour les deux cotés de l'interface (elles sont opposées et toujours du centre de la sous-structure voisine vers l'interface). Elles correspondent donc bien à la normale sortante de la sous-structure (pour que F = sigma n ).
 
Les normales sont calculées par élément et donc aux noeuds équivalents de l'interface pour lesquels les grandeurs F, W sont définies. Celles ci sont normalisées.
*/
template<class TV1, class INTER > void calculate_normales(INTER &Inter, TV1 &S) {
   //le maillage de peau des SST a deja ete mis a jour dans calculate_measure_G_SST.h
   for(unsigned q=0;q<Inter.side.size();++q) {
      //selection des elements de peau de la sous-structure voisine correspondant
      unsigned ii = Inter.vois[2*q];
      //preselection des elements du maillage de peau qui peuvent convenir (pour accelerer les calculs)
      Vec<typename TV1::template SubType<0>::T::TMESH::TM::Pvec,2> box = create_box_mesh(*Inter.side[q].mesh);

      Inter.side[q].neq.resize(Inter.side[q].mesh->elem_list.size()*DIM,0.);
      S[ii].mesh->update_skin();
      apply(S[ii].mesh->skin.elem_list,Calculate_Normales(),*S[ii].mesh.m,box,Inter.side[q]);
      S[ii].mesh.unload();
   }
}


/** \ingroup  Operateurs_inter
\brief Calcul des quantités géométriques par interfaces : centre de gravité, mesure et normales par élément
 
Le centre de gravité et la mesure sont déterminées uniquement sur le coté 0 de l'interface. 
Les normales pour un coté donné de l'interface sont toujours dirigées de la sous-structure vers l'interface.
*/
struct calculate_measure_G_neq_INTER {
   template<class INTER,class TV1> void operator() (INTER &Inter, TV1 &S) const {
      Inter.G=barycenter_constant_rho(*Inter.side[0].mesh);
      Inter.measure=measure(*Inter.side[0].mesh);
      calculate_normales(Inter,S);
   }
};

#endif CREATE_MESURE_G_NEQ_INTER
