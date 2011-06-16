

// assignation de ruptures initiales alleatoires
template<class TV2> Vec<typename TV2::template SubType<0>::T,TV2::template SubType<0>::T::dim> integre_bord(TV2 &Inter){
    typedef typename TV2::template SubType<0>::T TT;
    static const unsigned dim=TV2::template SubType<0>::T::dim;
    Vec<TT,TV2::template SubType<0>::T::dim> res;
    res.set(0.0);
    for(unsigned i=0;i < Inter.size();i++){
        if (Inter[i].type=="Ext"  and Inter[i].refCL==1){
            for(unsigned j=0;j<TV2::template SubType<0>::T::dim;j++){
                Vec<TT> unit;
                unit.resize(Inter[i].side[0].W.size());
                unit.set(0.0);
                Vec<int> rep=range(j,(unsigned)dim,Inter[i].side[0].W.size());
                for(unsigned k=0;k<rep.size();k++){
                    unit[rep[k]]=1.0;
                    }
                res[j]+=dot(unit,Inter[i].side[0].M*Inter[i].side[0].F);
                }
            }
      }
    return res;
}



// assignation de ruptures initiales aleatoires selon une loi uniforme
template<class TV2> void prerupture(TV2 &Inter,  unsigned nombre){
   typedef typename TV2::template SubType<0>::T TT;
   Vec<TT> liste_inter=Vec<TT>::random(Inter.size());
   Vec<unsigned> ind=sort_with_index(liste_inter);
   unsigned nb_inter_deg_ini=0;
   unsigned indfinal=0;
   while(nb_inter_deg_ini!=nombre){
      for(unsigned i=0;i<nombre-nb_inter_deg_ini;i++){
         if (Inter[ind[i+indfinal]].type=="Int" and Inter[ind[i+indfinal]].parammicro->degradable==1){
            Inter[ind[i+indfinal]].comp="Contact";
            nb_inter_deg_ini+=1;
            }
         }
         cout << "nbinterdeg " << nb_inter_deg_ini << endl;
         indfinal+=nombre-nb_inter_deg_ini;
         cout << "indfinal " <<indfinal << endl;
      }
}

