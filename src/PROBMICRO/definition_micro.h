
// classe parametres micro
template<unsigned dim_, class TT_> struct Param_Micro{
   Param_Micro() { Epini=0.0;sigmoy.set(0.);epsmoy.set(0.);num_Inter_degradee=-1;}
   static const unsigned dim = dim_;
   typedef TT_ T; // type connu de l'exterieur
   struct INTER{
      unsigned num;
      Vec<T> W;
      T Ep;
      T G;
   };
   Vec<INTER> Inter_select;
   Vec<unsigned> num_Inter_select;
   int num_Inter_degradee;
   T Epini;
   unsigned nbinterdeg;
   Vec<T,dim> sigmoy,epsmoy;
   
   T volumetot;
   Vec<T,dim_*2> boxdegrad;
   int signG;
   unsigned nb_inter_casse_ini;
};


// lecture des parametres micro 
template<class MICRO> void read_micro(MICRO &parammicro,const XmlNode &n){ 
   XmlNode np= n.get_element("parametres_micro");
   np.get_attribute("boxdegrad",parammicro.boxdegrad);
   Vec<typename MICRO::T> box
   nc.get_attribute( "box",box);
   parammicro.boxdegrad[0]=box[range(MICRO::dim)];
   propinter[i].boxparammicro.boxdegrad[1]=box[range(MICRO::dim,(unsigned)(2*MICRO::dim))];
   np.get_attribute("nbinterdeg",parammicro.nbinterdeg);
   np.get_attribute("nb_inter_casse_ini",parammicro.nb_inter_casse_ini);
//    string type_chargement;
//    np.get_attribute("type_chargement",type_chargement);
//    if (type_chargement=="effort")
//       parammicro.signG=1;
//    else 
//       parammicro.signG=-1;
}
