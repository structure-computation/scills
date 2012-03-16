
// classe parametres micro
struct Param_Micro{
   Param_Micro() { Epini=0.0;sigmoy.set(0.);epsmoy.set(0.);num_Inter_degradee=-1;}
   struct INTER{
      unsigned num;
      Vec<TYPEREEL> W;
      TYPEREEL Ep;
      TYPEREEL G;
   };
   Vec<INTER> Inter_select;
   Vec<unsigned> num_Inter_select;
   int num_Inter_degradee;
   TYPEREEL Epini;
   unsigned nbinterdeg;
   Vec<TYPEREEL,DIM> sigmoy,epsmoy;
   
   TYPEREEL volumetot;
   Vec<TYPEREEL,DIM*2> boxdegrad;
   int signG;
   unsigned nb_inter_casse_ini;
};


// lecture des parametres micro 
template<class MICRO> void read_micro(MICRO &parammicro,const XmlNode &n){ 
   XmlNode np= n.get_element("parametres_micro");
   np.get_attribute("boxdegrad",parammicro.boxdegrad);
   Vec<TYPEREEL> box
   nc.get_attribute( "box",box);
   parammicro.boxdegrad[0]=box[range(DIM)];
   propinter[i].boxparammicro.boxdegrad[1]=box[range(DIM,(unsigned)(2*DIM))];
   np.get_attribute("nbinterdeg",parammicro.nbinterdeg);
   np.get_attribute("nb_inter_casse_ini",parammicro.nb_inter_casse_ini);
//    string type_chargement;
//    np.get_attribute("type_chargement",type_chargement);
//    if (type_chargement=="effort")
//       parammicro.signG=1;
//    else 
//       parammicro.signG=-1;
}
