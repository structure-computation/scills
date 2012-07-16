/*
struct calcul_Y_moy_sst{
   template<class TE, class TV> void operator()(TE &e, TV &Y) const{
      Y+=e.Yde*measure(e);
   }
};

//comportement endommageable orthotrope : determination de d
template<class SST> void compt_sst_damage_uniforme (SST &S){

   typedef typename SST::T T;
   //1ere etape : sauvegarde du deplacement pour cette iteration (utile pour ne pas stocker les contraintes chap et def chap), mais demande un peu plus de calculs
   // Rq : deja fait dans la relaxation
   //mise a jour des deplacements
   S.f->get_result()=S.q;
   S.f->update_variables();
   //2eme etape : boucle tant que pour determination de d
   S.mesh.d=S.mesh.dmax;
   T erreur_d=1.;
   //cout << "dmax " << S.mesh.dmax <<endl;
   while(erreur_d>1e-6){
   
      //calcul de Yde par element : cf formulation
      S.f->call_after_solve_2();
      
      //calcul de Yd moyen sur la sst
      Vec<T,SST::dim> Ydmoy;
      Ydmoy.set(0.);
      apply(S.mesh.elem_list,calcul_Y_moy_sst(),Ydmoy);
      Ydmoy=Ydmoy/S.measure;
      //cout << "YDMOY "<< Ydmoy <<endl;
      //calcul des endommagements d et dp
      T Ydb=Ydmoy[0]+ S.mesh.b*Ydmoy[1];
      T Ymax=max(Ydb,S.mesh.Ydmax);
      S.mesh.Yd=Ymax;
      T d=1./sqrt(S.mesh.Yc)* (sqrt(Ymax)-sqrt(S.mesh.Yo));
      if (d<=0.){d=0.;}
      else if (d>=1. or sqrt(Ydmoy[1])>=sqrt(S.mesh.Ysp)){d=1.;}
      
      T dp=1./sqrt(S.mesh.Ycp)* (sqrt(Ymax)-sqrt(S.mesh.Yop));
      if (dp<=0.){dp=0.;}
      else if (dp>=1. or sqrt(Ydmoy[1])>=sqrt(S.mesh.Ysp)){dp=1.;}
      
      Vec<T,SST::dim> dnew;
      dnew.set(0.);
      dnew[0]=d;
      dnew[1]=dp;
      
      erreur_d = norm_2(dnew-S.mesh.d)/(norm_2(S.mesh.d)+1e-10);
      S.mesh.d=dnew;
      //cout << "d "<<S.mesh.d <<endl;
   }
   
   //mise a jour de depold 
   for(unsigned i=0;i<S.mesh.node_list.size();i++){
      Vec<unsigned> rep=range(i*SST::dim,(i+1)*SST::dim);
      S.mesh.node_list[i].depold=S.q[rep];
   }

}


struct affectation_valeur_d {
   template<class TE, class TV> void operator() (TE &e , TV &d) const {
      e.d=d;
   }
};

struct calcul_Yd_moyen_pli {
   template<class TE,class TV> void operator() (TE &e, TV &Ydmoy_pli) const {
      Ydmoy_pli+=e.Ye;
   }
};

struct CallAfterSolve_2 {
   template<class TE,class TM,class NameFormulation,class NameVariant=DefaultBehavior,class ScalarType_=typename TM::Tpos> void operator()(TE &e,Formulation<TM,elasticity_orthotropy_damage_stat_Qstat> &f) const {
      unsigned in[ TE::nb_nodes+1 ];
      for(unsigned i=0;i<TE::nb_nodes;++i)
            in[i] = f.indice_noda[ f.m->node_list.number(*e.node(i)) ];
      in[ TE::nb_nodes ] = f.indice_elem[TE::num_in_elem_list][e.number];
      CaracFormulationForElement<NameFormulation,TE,NameVariant,ScalarType>::after_solve_2(e,f,in);
   }
};

struct calcul_endommagement_elem {
   template<class TE,class SST> void operator() (TE &e, SST &S) const {
   typedef typename SST::T T;
   
   T erreur_d=1.;
   //cout << "dmax " << S.mesh.dmax <<endl;
   e.d=e.dmax;
   while(erreur_d>1e-6){
      for(unsigned j=0;j<e.voisins.size();j++)
         TM::TElemList::apply_on_down_cast(S.mesh.elem_list[e.voisins[j]],affectation_valeur_d(),e.d); 
      //calcul de Yd integre par element : cf formulation
      Vec<void *> elem_list;
      elem_list.resize(e.voisins.size());
      for(unsigned j=0;j<elem_list.size();j++)
         elem_list[j]=S.mesh.elem_list[j];
      S.f->call_after_solve_2( elem_list );
        
      //calcul de Yd moyen dans le pli
      Vec<T,2> Ydmoy;
      Ydmoy.set(0.);
      for(unsigned j=0;j<e.voisins.size();j++)
         TM::TElemList::apply_on_down_cast(m.elem_list[e.voisins[j]],calcul_Ydmoy_pli(),Ydmoy);
      Ydmoy=Ydmoy/(measure(e)*S.param_damage->epaisseur);

      //calcul des endommagements d et dp
      T Ydb=Ydmoy[0]+ S.param_damage->b*Ydmoy[1];
      T Ymax=max(Ydb,S.param_damage->Ydmax);
      e.Y=Ymax;
      T d=1./sqrt(S.param_damage->Yc)* (sqrt(Ymax)-sqrt(S.param_damage->Yo));
      if (d<=0.){d=0.;}
      else if (d>=1. or sqrt(Ydmoy[1])>=sqrt(S.param_damage->Ysp)){d=1.;}
      
      T dp=1./sqrt(S.param_damage->Ycp)* (sqrt(Ymax)-sqrt(S.param_damage->Yop));
      if (dp<=0.){dp=0.;}
      else if (dp>=1. or sqrt(Ydmoy[1])>=sqrt(S.param_damage->Ysp)){dp=1.;}
      
      Vec<T,SST::dim> dnew;
      dnew.set(0.);
      dnew[0]=d;
      dnew[1]=dp;
      
      erreur_d = norm_2(dnew-e.d)/(norm_2(e.d)+1e-10);
      e.d=dnew;
      //cout << "d "<<S.mesh.d <<endl;
   }
   
   }	
};
           

//comportement endommageable orthotrope : determination de d par element
template<class SST> void compt_sst_damage(SST &S){

   typedef typename SST::T T;
   
   //mise a jour des deplacements
   S.f->get_result()=S.q;   
   S.f->update_variables();             
   //assignation depchap (= dep)
   for(unsigned i=0;i<m.node_list.size();i++)
       S.mesh.node_list[i].depchap=m.node_list[i].dep;
               
   //boucle sur chaque element du maillage de reference
   apply(S.param_damage->mesh.elem_list,calcul_endommagement_elem(),S);      
  
   
   //1ere etape : sauvegarde du deplacement pour cette iteration (utile pour ne pas stocker les contraintes chap et def chap), mais demande un peu plus de calculs
   // Rq : deja fait dans la relaxation
   //mise a jour des deplacements
   S.f->get_result()=S.q;
   S.f->update_variables();
   //2eme etape : boucle tant que pour determination de d
   S.mesh.d=S.mesh.dmax;
   T erreur_d=1.;
   //cout << "dmax " << S.mesh.dmax <<endl;
   while(erreur_d>1e-6){
   
      //calcul de Yde par element : cf formulation
      S.f->call_after_solve_2();
      
      //calcul de Yd moyen sur la sst
      Vec<T,SST::dim> Ydmoy;
      Ydmoy.set(0.);
      apply(S.mesh.elem_list,calcul_Y_moy_sst(),Ydmoy);
      Ydmoy=Ydmoy/S.measure;
      //cout << "YDMOY "<< Ydmoy <<endl;
      //calcul des endommagements d et dp
      T Ydb=Ydmoy[0]+ S.mesh.b*Ydmoy[1];
      T Ymax=max(Ydb,S.mesh.Ydmax);
      S.mesh.Yd=Ymax;
      T d=1./sqrt(S.mesh.Yc)* (sqrt(Ymax)-sqrt(S.mesh.Yo));
      if (d<=0.){d=0.;}
      else if (d>=1. or sqrt(Ydmoy[1])>=sqrt(S.mesh.Ysp)){d=1.;}
      
      T dp=1./sqrt(S.mesh.Ycp)* (sqrt(Ymax)-sqrt(S.mesh.Yop));
      if (dp<=0.){dp=0.;}
      else if (dp>=1. or sqrt(Ydmoy[1])>=sqrt(S.mesh.Ysp)){dp=1.;}
      
      Vec<T,SST::dim> dnew;
      dnew.set(0.);
      dnew[0]=d;
      dnew[1]=dp;
      
      erreur_d = norm_2(dnew-S.mesh.d)/(norm_2(S.mesh.d)+1e-10);
      S.mesh.d=dnew;
      //cout << "d "<<S.mesh.d <<endl;
   }
   
   //mise a jour de depold 
   for(unsigned i=0;i<S.mesh.node_list.size();i++){
      Vec<unsigned> rep=range(i*SST::dim,(i+1)*SST::dim);
      S.mesh.node_list[i].depold=S.q[rep];
   }

}*/
