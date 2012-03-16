using namespace LMT;
using namespace std;


//allocation de t_post uniquement pour les sorties de resultat
template<class TV1, class TV2> void allocate_quantities_post(TV1 &S, TV2 &Inter,Param &process){
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
         S[i].t[pt].allocations(S[i].mesh.node_list_size*DIM,process);
      if (process.nom_calcul=="incr")
        for(unsigned pt=1;pt<(process.temps->nbpastemps+1);pt++)
          S[i].t_post[pt].allocations(S[i].mesh.node_list_size*DIM,process);

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
               Inter[i].side[j].t[pt].allocations(Inter[i].side[j].nodeeq.size()*DIM,process);
            for(unsigned pt=0;pt<(process.temps->nbpastemps+1);pt++)
              if(process.nom_calcul=="incr")
                  Inter[i].side[j].t_post[pt].allocations(Inter[i].side[j].nodeeq.size()*DIM,process);                        

            }
      }
      //allocation des parametres de comportements dependant du temps (exemple endommagement)
      Inter[i].param_comp->t.resize(nbpastemps+1);
      for(unsigned j=0;j<Inter[i].param_comp->t.size();j++)
         Inter[i].param_comp->t[j].allocate(Inter[i].side[0].nodeeq.size());
   }
    
}


//allocation des differentes quantites
/** \ingroup   Strategie_iterative
\brief Allocation de m�moire pour les sous-structures, interfaces et vecteurs macro. 

Pour chaque pas de temps les diff�rents vecteurs sont initialis�s � 0. On n'alloue pas la m�moire pour les diff�rents vecteurs en chaque pas de temps des ssts car ils ne sont stock�s qu'� convergence.
*/
template<class TV1, class TV2, class GLOB> void allocate_quantities(TV1 &S, TV2 &Inter,Param &process,GLOB &Global){
   unsigned nbpastemps;
   if(process.nom_calcul=="incr")
      nbpastemps=1;
   else if (process.nom_calcul=="latin")
      nbpastemps=process.temps->nbpastemps;
   else{std::cout << "Nom de clacul non implemente dans allocate.h" << endl;assert(0);}
   
   if (process.latin->save_depl_SST==1)
   for(unsigned i=0;i<S.size();++i){
      S[i].t.resize(nbpastemps+1);
      if(process.recopie_t_post==1)
         S[i].t=S[i].t_post;
/*      else
         S[i].t_post.resize(process.temps->nbpastemps+1);//utile pour le post traitement en incremental*/
      for(unsigned pt=1;pt<(nbpastemps+1);pt++)
         S[i].t[pt].allocations(S[i].mesh.node_list_size*DIM,process);
//       if (process.nom_calcul=="incr")
//         for(unsigned pt=1;pt<(process.temps->nbpastemps+1);pt++)
         //           S[i].t_post[pt].allocations(S[i].mesh.node_list_size*DIM,process);

   }
   if (process.rank>0 or process.size==1)
   for(unsigned i=0;i<Inter.size();++i){
      for(unsigned j=0;j<Inter[i].side.size();++j){
         Inter[i].side[j].t.resize(nbpastemps+1);
         if(process.recopie_t_post==1)
             Inter[i].side[j].t=Inter[i].side[j].t_post;
         else{
/*            Inter[i].side[j].t_post.resize(process.temps->nbpastemps+1); //utile pour le post traitement en incremental*/
            for(unsigned pt=0;pt<(nbpastemps+1);pt++)
                Inter[i].side[j].t[pt].allocations(Inter[i].side[j].nodeeq.size()*DIM,process);
//             for(unsigned pt=0;pt<(process.temps->nbpastemps+1);pt++)
//               if(process.nom_calcul=="incr")
//                   Inter[i].side[j].t_post[pt].allocations(Inter[i].side[j].nodeeq.size()*DIM,process);                        

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

