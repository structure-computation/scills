#ifndef CALCULATE_MEASURE_G_SST_H
#define CALCULATE_MEASURE_G_SST_H


/** \ingroup  Maillage_geometrie 
\brief Calcul de la mesure et centre de gravite des sst à partir du maillage
*/
struct calculate_measure_G_SST {
   template<class SST> void operator() (SST &S) const {
      S.mesh.load();
      S.G=barycenter_constant_rho(*S.mesh.m);
      S.measure=measure(*S.mesh.m);
      for(unsigned j=0;j<S.edge.size();++j) 
         S.edge[j].G = barycenter_constant_rho(*S.edge[j].mesh);

      S.mesh.unload();
   }
};



/** \ingroup  Maillage_geometrie
\brief Calcul de la taille du problème en nombre d'elements
*/
template<class TV1, class TV2>
void calcul_taille_probleme(TV1 &S, TV2 &Inter) {
    unsigned nbelem=0;
    unsigned nbnode=0;
    for(unsigned i=0;i<S.size();++i) {
        nbelem+=S[i].mesh.elem_list_size;
        nbnode+=S[i].mesh.node_list_size;
    }
    std::cout << "\t Taille du probleme : " << S.size() << " ssts, " << nbelem << " elements, " << nbnode << " noeuds"<< std::endl;
    std::cout << "\t Nbre d'interfaces : " << Inter.size()  << std::endl;
}

#endif //CALCULATE_MEASURE_G_SST_H
