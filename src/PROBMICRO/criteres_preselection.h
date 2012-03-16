


//calcul de la somme de <sign>^2+sigt^2 pour selection d'une interface  
template<class INTER> TYPEREEL calc_sign_sigt(INTER &Inter){
   
    Vec<TYPEREEL> Fn1 = Inter.side[0].Pn*Inter.side[0].F;
    Vec<TYPEREEL> Fn2 = Inter.side[1].Pn*Inter.side[1].F;
    Vec<TYPEREEL> Ft1 = Inter.side[0].Pt*Inter.side[0].F;
    Vec<TYPEREEL> Ft2 = Inter.side[1].Pt*Inter.side[1].F;
   
   // utilisation de la partie positive de la composante normale
   for(unsigned i=0;i<Fn1.size();++i){
      if (Fn1[i]<0.0) {Fn1[i]=0.0;}
      if (Fn2[i]>0.0) {Fn2[i]=0.0;}
   }
   
   // calcul d'un critere
   TYPEREEL sn1=dot(Fn1,Inter.side[0].M*Fn1);
   TYPEREEL sn2=dot(Fn2,Inter.side[0].M*Fn2);
   TYPEREEL sn = (sn1+sn2)/2;
   
   TYPEREEL st1=dot(Ft1,Inter.side[0].M*Ft1);
   TYPEREEL st2=dot(Ft2,Inter.side[0].M*Ft2);
   TYPEREEL st = (st1+st2)/2;
   
   TYPEREEL res;
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
