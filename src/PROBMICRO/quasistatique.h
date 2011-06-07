

//********************************************************************************************
//Procedure permettant le calcul des contraintes, deformations moyennes sur toute la structure
//********************************************************************************************
// stockage du resultat dans le fichier de donnee f
template<class TV1> void calcul_sigmoy_epsmoy_structure(TV1 &S, Param &process,std::ofstream &f){
      cout << "Calcul contraintes deformations" << endl;
      apply(S,assign_dep_contrainte());
      cout << "Calcul contraintes deformations moyennes" << endl;
      Vec<TYPEREEL,TV1::template SubType<0>::T::dim*(TV1::template SubType<0>::dim+1)/2 > sigmoy,epsmoy;
      sigmoy.set(0.0);
      epsmoy.set(0.0);
      apply(S,calcul_sigmoy_epsmoy_SST(),sigmoy,epsmoy);
      sigmoy /= process.volumetot;
      epsmoy /= process.volumetot;
      f << " sigmoy " << sigmoy << " epsmoy " << epsmoy << endl;
      cout << "Sigma moyen " << sigmoy << " Epsilon moyen " << epsmoy << endl;
   }


