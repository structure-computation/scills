


//calcul de la somme de <sign>^2+sigt^2 pour selection d'une interface  
template<class INTER> typename INTER::T calc_sign_sigt(INTER &Inter){
   typedef typename INTER::T TT;
   
   Vec<TT> Fn1 = Inter.side[0].Pn*Inter.side[0].F;
   Vec<TT> Fn2 = Inter.side[1].Pn*Inter.side[1].F;
   Vec<TT> Ft1 = Inter.side[0].Pt*Inter.side[0].F;
   Vec<TT> Ft2 = Inter.side[1].Pt*Inter.side[1].F;
   
   // utilisation de la partie positive de la composante normale
   for(unsigned i=0;i<Fn1.size();++i){
      if (Fn1[i]<0.0) {Fn1[i]=0.0;}
      if (Fn2[i]>0.0) {Fn2[i]=0.0;}
   }
   
   // calcul d'un critere
   TT sn1=dot(Fn1,Inter.side[0].M*Fn1);
   TT sn2=dot(Fn2,Inter.side[0].M*Fn2);
   TT sn = (sn1+sn2)/2;
   
   TT st1=dot(Ft1,Inter.side[0].M*Ft1);
   TT st2=dot(Ft2,Inter.side[0].M*Ft2);
   TT st = (st1+st2)/2;
   
   TT res;
   res = st + sn ;
   return res;
   
      
      
}         

//********************************************************************************************
// Procedure permettant de selectionner a priori des interfaces susceptibles de se degrader pour ce chargement
//********************************************************************************************
// si l'interface est non degradable ou est deja degradee le critere est mis a zero
// sinon on fait appel a calc_sign_sigt defini ci dessus
struct calcul_sign_sigt{
   template<class INTER> void operator()(INTER &Inter, unsigned i, Vec<typename INTER::T> &valinter ) const{
      if (Inter.parammicro->degradable==0 or Inter.comp=="Contact") {valinter[i]=0;}
      else {valinter[i] = calc_sign_sigt(Inter);}
      }
};
