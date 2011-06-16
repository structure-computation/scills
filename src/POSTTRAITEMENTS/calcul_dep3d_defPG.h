
struct calcul_dep3d_orthodefPG{
template<class SST> void operator()(SST &S) const {
   if (S.f->get_name()=="elasticity_orthotropy_DPG"){
      cout << "S.mesh.epshorsplan " << S.mesh.epshorsplan<< endl;
      for(unsigned i=0;i<S.mesh.node_list.size();++i){
         S.mesh.node_list[i].depDPG[0]=S.mesh.node_list[i].dep[0];
         S.mesh.node_list[i].depDPG[1]=S.mesh.node_list[i].dep[1];
         
         S.mesh.node_list[i].depDPG[2]=S.mesh.node_list[i].pos[0]*S.mesh.epshorsplan[0]+S.mesh.node_list[i].pos[1]*S.mesh.epshorsplan[1];
      }
   }
}
};

