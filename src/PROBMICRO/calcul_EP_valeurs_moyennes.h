
//valeur moyenne des deformations et contrainte par elements
struct add_sigmoy_epsmoy_elem{
   template<class TE,class TV> void operator()(TE &e,TV &sigmoy, TV &epsmoy) const {
        sigmoy+=e.sigma[0]*measure(e);
        epsmoy+=e.epsilon[0]*measure(e);}
};

//calcul de sigma moyen et epsilon moyen par SST
// --------> non prise en compte des SST dont le comportement est egal a 2
struct calcul_sigmoy_epsmoy_SST{
   template<class SST> void operator()(SST &S, Vec<typename SST::T,SST::dim*(SST::dim+1)/2 >  &sigmoy, Vec<typename SST::T,SST::dim*(SST::dim+1)/2 >  &epsmoy) const{
   if (S.typmat!= 2) {
      apply(S.mesh.elem_list,add_sigmoy_epsmoy_elem(),sigmoy,epsmoy);
      }
   }
};

//ajout de l'energie par element
struct add_ener_elem{
   template<class TE> void operator()(TE &e, typename TE::T &ener) const {
      ener+=e.ener;
   }
};

// calcul d'energie par SST : (sigma*epsilon*volume_element)/2
struct calcul_energie_SST{
template<class SST, class TV> void operator()(SST &S, typename SST::T &ener, TV &sigmoy, TV &epsmoy) const{
   //mise a jour des deplacements contraintes deformations pour la SST
   assign_dep_contrainte calcdep_sig;
   calcdep_sig(S);
   //calcul de l'energie de cette SST et ajout a ener
   apply(S.mesh.elem_list,add_ener_elem(),ener);
   //calcul sigmoy, epsmoy
   apply(S.mesh.elem_list,add_sigmoy_epsmoy_elem(),sigmoy,epsmoy);
   }
};

//calcul de l'energie potentielle de la structure, de sigma et epsilon moyen
template<class TV1,class PARAMMICRO> void calcul_EP_sigmoy_epsmoy(typename PARAMMICRO::T &Ep, PARAMMICRO &parammicro, TV1 &S){
   parammicro.sigmoy.set(0.);
   parammicro.epsmoy.set(0.);
   for(unsigned i=0;i<S.size();++i){
      calcul_energie_SST calcener;
      calcener(S[i],Ep,parammicro.sigmoy,parammicro.epsmoy);
      }
   parammicro.sigmoy/=parammicro.volumetot;
   parammicro.epsmoy/=parammicro.volumetot;
   //calcul du terme d'effort impose intervenant dans l'Ep
   //TYPEREEL Fimp=0.0;
   //apply(Inter,calcul_eff_imposes(),Fimp);
   //parammicro.Edini=parammicro.Edini-Fimp;
}

//calcul du terme a effort impose pour la determination de l'energie potentielle
struct calcul_eff_imposes{
   template<class INTER> void operator()(INTER &Inter, typename INTER::T &Fimp) const{
      if(Inter.type=="Ext" && Inter.comp=="effort"){
         Fimp+=dot(Inter.side[0].W,(*Inter.side[0].M*Inter.side[0].F));
         }
   }
};
