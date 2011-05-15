using namespace LMT;
using namespace std;

//allocations des quantites par SST
template<unsigned dim_, class TT_> void Sst<dim_,TT_>::Time::allocations(unsigned nbnode,Param &process)
{
//    if (process.rank>0 or process.size==1) oldq.resize(nbnode);//otable
//    if (process.rank>0 or process.size==1) oldq.set(0.0);
   if (process.rank>0 or process.size==1) q.resize(nbnode);
   if (process.rank>0 or process.size==1) q.set(0.0);
//    if (process.rank>0 or process.size==1) q1.resize(nbnode);//otable
//    if (process.rank>0 or process.size==1) q1.set(0.0);
//    if (process.rank>0 or process.size==1) q2.resize(nbnode);//otable
//    if (process.rank>0 or process.size==1) q2.set(0.0);
   
}

//allocations des quantites par INTER
template<unsigned dim_, class TT_> void Interface<dim_,TT_>::Side::Time::allocations(unsigned nbnodeeq,Param &process)
{
//    if (process.rank>0 or process.size==1) Wd.resize(nbnodeeq);//otable
//    if (process.rank>0 or process.size==1) Wd.set(0.0);
//    if (process.rank>0 or process.size==1) Qd.resize(nbnodeeq);//otable
//    if (process.rank>0 or process.size==1) Qd.set(0.0);
//    if (process.rank>0 or process.size==1) Wtemp.resize(nbnodeeq);//otable
//    if (process.rank>0 or process.size==1) Wtemp.set(0.0);
   if (process.rank>0 or process.size==1) F.resize(nbnodeeq);
   if (process.rank>0 or process.size==1) F.set(0.0);
   if (process.rank>0 or process.size==1) Wp.resize(nbnodeeq);
   if (process.rank>0 or process.size==1) Wp.set(0.0);
   if (process.rank>0 or process.size==1) W.resize(nbnodeeq);
   if (process.rank>0 or process.size==1) W.set(0.0);
   if (process.rank>0 or process.size==1) oldW.resize(nbnodeeq);
   if (process.rank>0 or process.size==1) oldW.set(0.0);
   if (process.rank>0 or process.size==1) Fchap.resize(nbnodeeq);
   if (process.rank>0 or process.size==1) Fchap.set(0.0);
   if (process.rank>0 or process.size==1) Wchap.resize(nbnodeeq);
   if (process.rank>0 or process.size==1) Wchap.set(0.0);
   if (process.rank>0 or process.size==1) Wpchap.resize(nbnodeeq);
   if (process.rank>0 or process.size==1) Wpchap.set(0.0);
   if (process.rank>0 or process.size==1) WtildeM.resize(nbnodeeq);
   if (process.rank>0 or process.size==1) WtildeM.set(0.0);
   if (process.rank>0 or process.size==1) oldF.resize(nbnodeeq);
   if (process.rank>0 or process.size==1) oldF.set(0.0);
   if (process.rank>0 or process.size==1) oldWp.resize(nbnodeeq);
   if (process.rank>0 or process.size==1) oldWp.set(0.0);
//    if (process.rank>0 or process.size==1) F1.resize(nbnodeeq);//otable
//    if (process.rank>0 or process.size==1) F1.set(0.0);
//    if (process.rank>0 or process.size==1) W1.resize(nbnodeeq);//otable
//    if (process.rank>0 or process.size==1) W1.set(0.0);
//    if (process.rank>0 or process.size==1) F2.resize(nbnodeeq);//otable
//    if (process.rank>0 or process.size==1) F2.set(0.0);
//    if (process.rank>0 or process.size==1) W2.resize(nbnodeeq);//otable
//    if (process.rank>0 or process.size==1) W2.set(0.0);
//    if (process.rank>0 or process.size==1) Wp1.resize(nbnodeeq);//non utilise
//    if (process.rank>0 or process.size==1) Wp1.set(0.0);
//    if (process.rank>0 or process.size==1) Wp2.resize(nbnodeeq);//non utilise
//    if (process.rank>0 or process.size==1) Wp2.set(0.0);
}

//allocations pour le probleme macro
template<unsigned dim_, class TT_> void Glob<dim_,TT_>::allocations(Param &process)
{
   bigW.resize(process.multiscale->sizeM);
   bigW.set(0.0);
   bigF.resize(process.multiscale->sizeM);
   bigF.set(0.0);
}

//allocation des differentes quantites
/** \ingroup   Strategie_iterative
\brief Allocation de mémoire pour les sous-structures, interfaces et vecteurs macro. 

Pour chaque pas de temps les différents vecteurs sont initialisés à 0. On n'alloue pas la mémoire pour les différents vecteurs en chaque pas de temps des ssts car ils ne sont stockés qu'à convergence.
*/
template<class TV1, class TV2, class GLOB> void allocate_quantities(TV1 &S, TV2 &Inter,Param &process,GLOB &Global){
   unsigned nbpastemps;
   if(process.nom_calcul=="incr")
      nbpastemps=1;
   else if (process.nom_calcul=="latin")
      nbpastemps=process.temps->nbpastemps;
   else{std::cout << "Nom de clacul non implemente dans allocate.h" << std::endl;assert(0);}
   
   if (process.latin->save_depl_SST==1)
   for(unsigned i=0;i<S.size();++i){
      S[i].t.resize(nbpastemps+1);
      if(process.recopie_t_post==1)
         S[i].t=S[i].t_post;
      else
         S[i].t_post.resize(process.temps->nbpastemps+1);//utile pour le post traitement en incremental
      for(unsigned pt=1;pt<(nbpastemps+1);pt++)
         S[i].t[pt].allocations(S[i].mesh.node_list_size*TV2::template SubType<0>::T::dim,process);
      if (process.nom_calcul=="incr")
        for(unsigned pt=1;pt<(process.temps->nbpastemps+1);pt++)
          S[i].t_post[pt].allocations(S[i].mesh.node_list_size*TV2::template SubType<0>::T::dim,process);

   }
   if (process.rank>0 or process.size==1)
   for(unsigned i=0;i<Inter.size();++i){
      for(unsigned j=0;j<Inter[i].side.size();++j){
         Inter[i].side[j].t.resize(nbpastemps+1);
         if(process.recopie_t_post==1)
             Inter[i].side[j].t=Inter[i].side[j].t_post;
         else{
            Inter[i].side[j].t_post.resize(process.temps->nbpastemps+1); //utile pour le post traitement en incremental
            for(unsigned pt=0;pt<(nbpastemps+1);pt++)
               Inter[i].side[j].t[pt].allocations(Inter[i].side[j].nodeeq.size()*TV2::template SubType<0>::T::dim,process);
            for(unsigned pt=0;pt<(process.temps->nbpastemps+1);pt++)
              if(process.nom_calcul=="incr")
                  Inter[i].side[j].t_post[pt].allocations(Inter[i].side[j].nodeeq.size()*TV2::template SubType<0>::T::dim,process);                        

            }
      }
      //allocation des parametres de comportements dependant du temps (exemple endommagement)
      Inter[i].param_comp->t.resize(nbpastemps+1);
      for(unsigned j=0;j<Inter[i].param_comp->t.size();j++)
         Inter[i].param_comp->t[j].allocate(Inter[i].side[0].nodeeq.size());
   }

   if(process.recopie_t_post!=1){   
      Global.allocations(process);
      process.latin->error.resize(process.latin->nbitermax+1);
      process.latin->error.set(0.0);
   }

}

