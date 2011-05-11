//fichiers contenant les differentes declarations des fonctions utilisees dans le programme principal
#include "GeometryUser.h"
#include "DataUser.h"
using namespace LMT;
using namespace std;

template<class TV1, class TV2, class TV3, class TV4,class TV5> void multiscale_geometry_mesh(const XmlNode &n,TV1 &S,TV2 &Inter, Param &process, TV5 &CL,TV3 &Stot,TV3 &SubS,TV4 &SubI) ;
template<class TV1, class TV2, class TV3, class TV4,class TV5> void multiscale_geometry_mesh(DataUser &data_user, GeometryUser &geometry_user,TV1 &S,TV2 &Inter, Param &process, TV5 &CL,TV3 &Stot,TV3 &SubS,TV4 &SubI) ;
template<class TV1, class TV2> void assignation_materials_property_SST(const XmlNode &n, TV1 &S, TV2 &Inter,Param &process) __attribute__((noinline));
template<class TV1, class TV2> void assignation_materials_property_INTER(const XmlNode &n, TV2 &Inter, TV1 &S, Param &process) __attribute__((noinline));
template <class TV1, class TV2, class TV3, class TV4>  void multiscale_operateurs(const XmlNode &n,TV1 &S, TV1 &SubS,TV2 &Inter, TV3 &SubI,Param &process,  TV4 &Global) ;
template <class TV1,class TV2, class TV3, class TV4, class TV5, class TV6> void multiscale_iterate_latin(TV1 &S,TV2 &SubS, TV3 &Inter, TV4 &SubI, Param &process, TV5 &Global, TV6 &CL);
template <class TV1,class TV2, class TV3, class TV4, class TV5, class TV6> void multiscale_iterate_incr(TV1 &S,TV2 &SubS, TV3 &Inter, TV4 &SubI, Param &process, TV5 &Global, TV6 &CL,const XmlNode &n);

