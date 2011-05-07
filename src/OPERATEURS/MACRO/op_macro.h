using namespace LMT;
using namespace std;

template<class TV1> void renum_loc(TV1 &S, Vec<unsigned> &repSloc, Vec<unsigned> &repSglob, Vec<unsigned> &numI){

   Vec<unsigned > newrepSloc;
   repSglob.append(repSloc);
   
   for(unsigned i=0;i<repSloc.size();++i){
      unsigned k=repSloc[i];
      Vec<int> vois=S[k].vois;
      for(unsigned j=0;j<vois.size();++j){
         int jj=vois[j];
         if (jj!=-1){
            if ( (find(repSglob,_1==(unsigned)jj)==0) and (find(newrepSloc,_1==(unsigned)jj)==0) )
               newrepSloc.push_back((unsigned)jj); // ajout du voisin s'il n'est pas deja dans la liste repSglob et si ce n'est pas le bord
         }
         unsigned q=S[k].edge[j].internum;
         if ( find(numI,_1==q)==0 )
            numI.push_back(q); // ajout du numero d'interface s'il n'est pas deja dans numI
      }
   }
   repSloc=newrepSloc;
}

/** \ingroup Operateurs_macro
\brief Renumerotation des interfaces en utilisant une numerotation frontale
*/
template<class TV1> void renumerotation_interfaces(TV1 &S, Vec<unsigned > &renum){
   Vec<unsigned > repSloc(0);
   Vec<unsigned > repSglob;
   do{
      renum_loc(S,repSloc,repSglob,renum);
   }
   while(repSglob.size() !=S.size());
   Vec<unsigned> num;
   num.resize(renum.size());
   for(unsigned i=0;i<num.size();++i){
      unsigned q=renum[i];
      num[q]=i;
   }
   renum=num;
}

//*******************************************************
// reperage des ddls d'interface dans le probleme macro
//*******************************************************
/** \ingroup Operateurs_macro
 \brief Rep�rage des ddls d'interfaces dans le probl�me macro.
 
 Dans cette proc�dure, on d�termine la taille du probl�me macro et on assigne le champ Interface::repddl pour rep�rer les ddls macro d'une interface dans le probl�me macro. Par d�faut, on utilise la num�rotation correspondant � la position de l'interface dans le vecteur d'interface jusqu'au nombre de fonction macro sur cette interface.
 
 Rq : La renum�rotation des interfaces est non n�cessaire car une renum�rotation efficace est directement effectu�e lors de la factorisation.
 */
template<class TV1, class TV2> void Repere_ddl_Inter(TV1 &S, TV2 &Inter, Param &process){
   unsigned sizeM=0;
   for(unsigned q=0;q<Inter.size();++q){
      Inter[q].repddl=range(sizeM,sizeM+Inter[q].nb_macro_espace);
      sizeM+=Inter[q].nb_macro_espace; 
      }
   process.multiscale->sizeM=sizeM;  
}

//*******************************************************
// Assemblage du probleme macro global
//*******************************************************
/**\ingroup Operateurs_macro
 \brief Assemblage du probl�me macro.
 
 En bouclant sur les sous-structures et connaissant les interfaces adjacentes et leur position dans le probl�me macro, on ajoute l'op�rateur homog�n�is� par sous-structure � la bonne place dans la matrice macro. Celle ci est ensuite directement factoris� apr�s p�nalisation (on ne la stocke pas la matrice macro dans Glob)
 */
template<class TV1,class TV2> Mat<TYPEREEL, Sym<>, SparseLine<> > Assem_prob_macro(TV1 &S, TV2 &Inter, Param &process){
   Mat<TYPEREEL, Sym<>, SparseLine<> > bigK;
   bigK.resize(process.multiscale->sizeM);
//    cout << bigK.size() << endl;
   for(unsigned i=0;i<S.size();++i){
      Vec<unsigned> LErep,Krep;
      for(unsigned j=0;j<S[i].edge.size();++j){
         unsigned q=S[i].edge[j].internum;
         //reperage des inconnues dans LE 
         LErep.append(S[i].edge[j].repLE);
         //reperage des inconnues dans bigK
         Krep.append(Inter[q].repddl);
      }
      //TODO macro adaptation des signes de LE pour approche en effort
      bigK[Krep]+=S[i].LE(LErep,LErep);
      //cout << "SST " << i << endl;
      //cout << LErep << endl;
      //cout << Krep << endl;
      //cout << S[i].LE(LErep,LErep).diag() << endl;
   }
   return bigK;        
};

//************************************
// Reperage des ddls macro a bloquer
//************************************
//blocage mvt de corps rigide particulier
/** \ingroup Operateurs_macro
\brief Blocage des mouvements de corps rigides d�finis par l'utilisateur dans RBM::mvts_bloques

Lorsque seuls quelques mouvements de corps rigides ne sont pas bloqu�s (ex : utilisation d'une seule condition de sym�trie), le programme ne peut pas d�terminer tout seul les mouvements � bloquer. Cependant dans le cas o� toutes les conditions aux limites sont des conditions en effort, les 6 mouvements sont bloqu�s automatiquement.

Cette proc�dure est �crite d'une part pour le 2D et pour le 3D d'autre part. On rep�re parmi les interfaces celles qui poss�dent une direction principale (la 3�me) parall�le � un axe demand� par l'utilisateur (ex : blocage de la translation selon x : Tx, on recherche une direction colin�aire � x). On bloque ensuite la translation selon cet axe si celle ci est demand�e, ou la rotation autour de cet axe. On ajoute les ddl macro correspondant � la liste des ddl bloqu�s.
*/
template<class T> void bloqrbm(Vec<Interface<3,T> > &Inter, Param &process,Vec<unsigned> &repddlMbloq){
// definition des translations possibles
   Vec<string> translations("Tx","Ty","Tz");
   Vec<Vec<T> > vectbase(Vec<T>(1.,0.,0.),Vec<T>(0.,1.,0.),Vec<T>(0.,0.,1.));
   T eps=1e-10;

   for(unsigned l=0;l<translations.size();++l){
      if (find(process.rbm.mvts_bloques,_1==translations[l])==1){
         for(unsigned q=0;q<Inter.size();++q){
            if( length(Inter[q].BPI[2]-vectbase[l])<eps or length(-1.0*Inter[q].BPI[2]-vectbase[l])<eps ){
               repddlMbloq.push_back(Inter[q].repddl[2]);// blocage translation selon N3
               break;
            }
         }
      }
   }
   
// definition des rotations possibles
   Vec<string> rotations("Rx","Ry","Rz");
   Vec<Vec<T> > vectbase2(Vec<T>(1.,0.,0.),Vec<T>(0.,1.,0.),Vec<T>(0.,0.,1.));


   for(unsigned l=0;l<rotations.size();++l){
      if (find(process.rbm.mvts_bloques,_1==rotations[l])==1){
         for(unsigned q=0;q<Inter.size();++q){
            if( length(Inter[q].BPI[2]-vectbase2[l])<eps or length(-1.0*Inter[q].BPI[2]-vectbase2[l])<eps ){
               repddlMbloq.push_back(Inter[q].repddl[5]); // blocage rotation autour de N3
               break;
            }
         }
      }
   }

   remove_doubles(repddlMbloq);

}

/** \ingroup Operateurs_macro
\brief Blocage des ddl macro en fonction des conditions aux limites

On rep�re les interfaces � d�placement impos�s pour lesquels tous les ddls macro correspondant sont stock�s. Pour les interfaces de type sym�trie ou d�placement normal impos�,on ne rep�re que les ddls macro selont la direction 3 de la base principale d'inertie de l'interface (l'interface est n�cessairement plane pour avoir sym�trie).

Si aucun ddl macro n'est bloqu�, on rep�re les 6 ddls d'une interface permettant de bloquer les modes de corps rigides de la structure enti�re, de fa�on automatique. Enfin pour bloquer seulement quelques modes de corps rigides, l'utilisateur doit sp�cifier les translations et rotations � bloquer (Tx, Ty, Tz, Rx, Ry et Rz) et la fonction bloqrbm() d�tecte automatiquement les ddls macro � bloquer.
*/
template<class T> Vec<unsigned> macro_CL(Vec<Interface<3,T> > &Inter, Param &process){
   //creation du vecteur contenant les ddls a bloquer
   Vec<unsigned> repddlMbloq;
   bool bloq=0;
   for(unsigned q=0;q<Inter.size();++q){
      if (Inter[q].type=="Ext" and Inter[q].comp=="depl"){
         //ddl bloques
         repddlMbloq.append(Inter[q].repddl);
         bloq=1;
      }
      else if(Inter[q].type=="Ext" and ( Inter[q].comp=="sym" or Inter[q].comp=="depl_normal" )){
         Vec<unsigned,3> repimp(2,3,4);
         repddlMbloq.append(Inter[q].repddl[repimp]);
         bloq=1;
      }
   }
   
   if (bloq==0){ // blocage des mvts de corps rigide
      for(unsigned q=0;q<Inter.size();++q){
         if (Inter[q].type=="Ext"){
            cout << "Blocage mvts corps rigide : interface " << q << endl;
            repddlMbloq.append(Inter[q].repddl[range(6)]);
            break;
         }
      }
   }
   else if (process.rbm.bloq==1){
      cout << "\t Blocage mvts corps rigide selon mvts_bloques"  << endl;
      bloqrbm(Inter, process,repddlMbloq);
   }
  
   return repddlMbloq;
}

// idem en dimension 2
//blocage mvt de corps rigide definis par l'utilisateur dans rbm.mvts_bloques par Tx, Ty, Mz
/** \ingroup Operateurs_macro
\brief Blocage des mouvements de corps rigides d�finis par l'utilisateur dans RBM::mvts_bloques

Lorsque seuls quelques mouvements de corps rigides ne sont pas bloqu�s (ex : utilisation d'une seule condition de sym�trie), le programme ne peut pas d�terminer tout seul les mouvements � bloquer. Cependant dans le cas o� toutes les conditions aux limites sont des conditions en effort, les 6 mouvements sont bloqu�s automatiquement.

Cette proc�dure est �crite d'une part pour le 2D et pour le 3D d'autre part. On rep�re parmi les interfaces celles qui poss�dent une direction principale (la 3�me) parall�le � un axe demand� par l'utilisateur (ex : blocage de la translation selon x : Tx, on recherche une direction colin�aire � x). On bloque ensuite la translation selon cet axe si celle ci est demand�e, ou la rotation autour de cet axe. On ajoute les ddl macro correspondant � la liste des ddl bloqu�s.
*/
template<class T> void bloqrbm(Vec<Interface<2,T> > &Inter, Param &process,Vec<unsigned> &repddlMbloq){
// definition des translations possibles
   Vec<string> translations("Tx","Ty");
   Vec<Vec<T> > vectbase(Vec<T>(1.,0.),Vec<T>(0.,1.));
   T eps=1e-10;

   for(unsigned l=0;l<translations.size();++l){
      if (find(process.rbm.mvts_bloques,_1==translations[l])==1){
         for(unsigned q=0;q<Inter.size();++q){
            if( length(Inter[q].BPI[1]-vectbase[l])<eps or length(-1.0*Inter[q].BPI[1]-vectbase[l])<eps ){
               repddlMbloq.push_back(Inter[q].repddl[1]);// blocage translation selon N3
               break;
            }
         }
      }
   }
   
// definition des rotations possibles
   Vec<string> rotations("Rz");
   Vec<Vec<T,2>,1> vectbase2;
   vectbase2[0]=Vec<T>(0.,1.);


   for(unsigned l=0;l<rotations.size();++l){
      if (find(process.rbm.mvts_bloques,_1==rotations[l])==1){
         for(unsigned q=0;q<Inter.size();++q){
            if( length(Inter[q].BPI[1]-vectbase2[l])<eps or length(-1.0*Inter[q].BPI[1]-vectbase2[l])<eps ){
               repddlMbloq.push_back(Inter[q].repddl[2]); // blocage rotation autour de N3
               break;
            }
         }
      }
   }

   remove_doubles(repddlMbloq);
  
}

/** \ingroup Operateurs_macro
\brief Blocage des ddl macro en fonction des conditions aux limites

On rep�re les interfaces � d�placement impos�s pour lesquels tous les ddls macro correspondant sont stock�s. Pour les interfaces de type sym�trie ou d�placement normal impos�,on ne rep�re que les ddls macro selont la direction 3 de la base principale d'inertie de l'interface (l'interface est n�cessairement plane pour avoir sym�trie).

Si aucun ddl macro n'est bloqu�, on rep�re les 6 ddls d'une interface permettant de bloquer les modes de corps rigides de la structure enti�re, de fa�on automatique. Enfin pour bloquer seulement quelques modes de corps rigides, l'utilisateur doit sp�cifier les translations et rotations � bloquer (Tx, Ty, Tz, Rx, Ry et Rz) et la fonction bloqrbm() d�tecte automatiquement les ddls macro � bloquer.
*/
template<class T> Vec<unsigned> macro_CL(Vec<Interface<2,T> > &Inter, Param &process){
   //creation du vecteur contenant les ddls a bloquer
   Vec<unsigned> repddlMbloq;
   bool bloq=0;
   for(unsigned q=0;q<Inter.size();++q){
      if (Inter[q].type=="Ext" and Inter[q].comp=="depl"){
         //ddl bloques
         repddlMbloq.append(Inter[q].repddl);
         bloq=1;
      }
      else if(Inter[q].type=="Ext" and (Inter[q].comp=="sym" or Inter[q].comp=="depl_normal") ){
         Vec<unsigned,2> repimp(1,2);
         repddlMbloq.append(Inter[q].repddl[repimp]);
         bloq=1;
      }
   }

   if (bloq==0 && process.rbm.bloq==0){ // blocage des mvts de corps rigide
      for(unsigned q=0;q<Inter.size();++q){
         if (Inter[q].type=="Ext"){
            cout << "\t Blocage mvts corps rigide : interface " << q << endl;
            repddlMbloq.append(Inter[q].repddl[range(3)]);
            break;
         }
      }

   }
   else if (process.rbm.bloq==1){
      cout << "\t Blocage mvts corps rigide selon mvts_bloques"  << endl;
      bloqrbm(Inter, process,repddlMbloq);
   }

   return repddlMbloq;
}
