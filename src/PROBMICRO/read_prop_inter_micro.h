
// lecture des proprietes materiau d'interface
template<class TV4> void read_prop_inter_micro(TV4 &propintermicro,const XmlNode &n){        
 
   XmlNode nmat=n.get_element("proprietes_interfaces");
   unsigned nbmat = nmat.nb_elements("coefficients_micro");

   propintermicro.resize(nbmat);

   for(unsigned i=0;i<nbmat;++i){
        XmlNode nc = nmat.get_element("coefficients_micro",i);
        nc.get_attribute( "coeffrottement", propintermicro[i].coeffrottement);
        nc.get_attribute( "Gcrit", propintermicro[i].Gcrit);
        nc.get_attribute( "name", propintermicro[i].name);
        nc.get_attribute( "Yc", propintermicro[i].Yc);
        nc.get_attribute( "Yo", propintermicro[i].Yo);
        nc.get_attribute( "n", propintermicro[i].n);
        nc.get_attribute( "alpha", propintermicro[i].alpha);
        nc.get_attribute( "gamma", propintermicro[i].gamma);
        nc.get_attribute( "kn", propintermicro[i].kn);
        nc.get_attribute( "knc", propintermicro[i].knc);
        nc.get_attribute( "kt", propintermicro[i].kt);
   }
  
};

