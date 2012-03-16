//#include "definition_PARAM_DAMAGE_SST.h"
/*
using namespace LMT;
using namespace std;

struct find_elem_vois {
   template<class TE,class TV> void operator() (TE &e , TV &G, unsigned &num, bool &flag) const{
      double eps=1e-5;
      typedef typename TE::Pvec Pvec;
      if (flag==0){
         Pvec G_i = center(e);
         Pvec normale_i=vect_prod(e.node(1)->pos - e.node(0)->pos,
                                 e.node(2)->pos - e.node(0)->pos);
         if (abs(G[0]-G_i[0])<=eps and abs(G[1]-G_i[1])<=eps){
         //if (abs(dot(G-G_i,normale_i)-1)<=eps or abs(dot(G-G_i,normale_i)+1)<=eps){
            e.elem.push_back(num);
            flag=1;
         }
      }
   }
};

struct affect_num {
   template<class TE,class T> void operator() (TE &e, T &num_ei) const{
      e.d[0]=num_ei;
   }
};

struct extract_vois {
   template<class TE,class T> void operator() (TE &e, T &voisins) const{
      voisins=e.elem;
   }
};

struct free_voisins {
   template<class TE> void operator() (TE &e) const{
      e.elem.free();
   }
};

template<class TMI, class TM> void reperage_elem(TMI &minter, TM &m){
   typedef typename TM::TNode::Pvec Pvec;
   //on vide la valeur 0 existant par defaut dans le champ voisin
   apply(minter.elem_list,free_voisins());
   //boucle sur les elements volumiques et recherche des elements ayant les memes abscisses dans le plan que la nappe consideree
   //cout << "Attention selection valide uniquement pour des plis plans !" << endl;
   for(unsigned i=0;i<m.elem_list.size();i++){
      Pvec G = center(*m.elem_list[i]);
      bool flag=0;
      apply(minter.elem_list,find_elem_vois(),G,i,flag);
   }
   //pour verification
//    for(unsigned i=0;i<minter.elem_list.size();i++){
//       Vec<unsigned> voisins;
//       TMI::TElemList::apply_on_down_cast(minter.elem_list[i],extract_vois(),voisins);
//       //cout << i << " vois "<<voisins << endl;
//       for(unsigned j=0;j<voisins.size();j++)
//          TM::TElemList::apply_on_down_cast(m.elem_list[voisins[j]],affect_num(),i);
//       
//    }
}



/** \ingroup  Mesomodele 
\brief Création du maillage surfacique représentatif de la sous-structure

A la fin de cette procédure, un maillage surfacique est créé et contient pour chaque élément les numéros des éléments volumiques de la sous-structure situés dans l'épaisseur.
- On repère tout d'abord l'interface entre deux sous-structures de materiaux differents. Le maillage de l'interface correspond alors au maillage recherché.
- On recherche ensuite les éléments situés dans un même pli (sous-structure) verticalement. Pour l'instant la procédure ne fonctionne que pour des plis plans. TODO : Utiliser les normales aux éléments pour les plis courbes.
*/
template<class TV1, class TV2> void create_mesh_meso(TV1 &S, TV2 &Inter){
   for(unsigned i=0;i<S.size();i++){
      bool flag=0;
      
   //1ere etape : creation du maillage surfacique représentatif
      for(unsigned j=0;j<S[i].vois.size();j++){
         int numvois = S[i].vois[j];
         if(numvois != -1)
            if(S[numvois].typmat != S[i].typmat){
               unsigned q = S[i].edge[j].internum;
               S[i].param_damage->mesh = Inter[q].side[0].mesh;
               //translation du maillage
               translate_mesh(S[i].param_damage->mesh,S[i].G - Inter[q].G);
               flag = 1;
            }
         if (flag==1) break;
      }
      
   //2eme etape : reperage des elements
   reperage_elem(S[i].param_damage->mesh, S[i].mesh);
      
   }
}

/** \ingroup  Mesomodele 
\brief Assignation des parametres meso aux maillages representatifs par sous-structure

Les differentes quantites necessaires au modele meso sont appliquees au maillage representatif de la sous-structure
*/
struct assignation_material_to_SST_meso{
   template<class SST, class TV3> void operator()(SST &S,TV3 &matprop) const{
         S.param_damage->Yo = matprop[S.typmat].param_damage.Yo;
         S.param_damage->Yop = matprop[S.typmat].param_damage.Yop;
         S.param_damage->Ysp = matprop[S.typmat].param_damage.Ysp;
         S.param_damage->Yc = matprop[S.typmat].param_damage.Yc;
         S.param_damage->Ycp = matprop[S.typmat].param_damage.Ycp;
         S.param_damage->b = matprop[S.typmat].param_damage.b;
   }
};

struct calc_measure {
   template<class TE> void operator() (TE &e, typename TE::T &mes) const {
      mes +=measure(e);
   }
};

struct calcul_epaisseur_pli {
   template<class TE,class TM> void operator() (TE &e, TM &m, typename TM::Tpos &epaisseur) const {
      typename TM::Tpos mesi=0.;
      typename TM::Tpos mese=0.;
      mese=measure(e);
      for(unsigned i=0;i<e.elem.size();i++){
         unsigned num=e.elem[i];
         TM::TElemList::apply_on_down_cast(m.elem_list[num],calc_measure(),mesi);
      }
      epaisseur = mesi/mese;
   }
};

/** \ingroup  Mesomodele 
\brief Procedure principale pour la creation et l'initialisation du maillage meso

On aloue tout d'abord la memoire pour la structure PARAM_DAMAGE_SST contenue dans chaque sous-structure.

On détermine ensuite le maillage représentatif stocké dans param_damage.mesh et contenant le repérage des différents éléments de la sous-structure.

On assigne enfin les propriétés matériaux pour le maillage représentatif et on calcule l'épaisseur du pli.
*/
template<class TV1, class TV2,class TV3> void initialisation_mesh_meso(TV1 &S, TV2 &Inter, TV3 &matprop){
   for(unsigned i=0;i<S.size();i++)
      S[i].param_damage ;//= new PARAM_DAMAGE_SST<TV1::template SubType<0>::T::dim,TYPEREEL>;
   
   create_mesh_meso(S,Inter);
   
   apply(S,assignation_material_to_SST_meso(),matprop);  
   
   //calcul de l'epaisseur du pli
   for(unsigned i=0;i<S.size();i++){
      apply(S[i].param_damage->mesh.elem_list,calcul_epaisseur_pli(),S[i].mesh,S[i].param_damage->epaisseur);
   cout << S[i].param_damage->epaisseur << endl;
      }
//    DisplayParaview dp;
//    for(unsigned i=0;i<S.size();i++){
//       dp.add_mesh(S[i].param_damage->mesh);
//       dp.add_mesh(S[i].mesh,"tmp/paraview",Vec<string>("d"));
//    }
//    dp.exec();
}
*/
