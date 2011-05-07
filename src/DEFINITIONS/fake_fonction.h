#ifndef FAKE_FONCTION_H
#define FAKE_FONCTION_H

// #include "definition_GLOB.h"
//#include "definition_fcts.h"
//fichiers contenant les differentes declarations des fonctions utilisees dans le programme principal
using namespace LMT;
using namespace std;



void fake_fonction(){
//pour le 3d
    XmlNode n;
    Param process;
    
#ifdef DIMENSION3
Vec<Sst<3,TYPEREEL> > S;
Vec<Interface<3,TYPEREEL> > Inter;
Vec<Boundary<3,TYPEREEL> > CL;
// Vec<VecPointedValues<Sst<3,TYPEREEL> > > SubS,Stot;
// Vec<VecPointedValues<Interface<3,TYPEREEL> > > SubI;
// Glob<3,TYPEREEL> Global;

multiscale_geometry_mesh(n,S, Inter, process, CL) ;
/*assignation_materials_property_SST(n, S, Inter,process);
assignation_materials_property_INTER(n, Inter, S, process);
multiscale_operateurs(n,Stot, SubS,Inter,SubI, process,Global);
multiscale_iterate_latin(S,SubS, Inter, SubI, process,Global,CL) ;
multiscale_iterate_incr(S,SubS, Inter,SubI, process, Global, CL,n) ;*/
#endif

//pour le 2d
#ifdef DIMENSION2
Vec<Sst<2,TYPEREEL> > S2;
Vec<Interface<2,TYPEREEL> > Inter2;
Vec<Boundary<2,TYPEREEL> > CL2;
// Vec<VecPointedValues<Sst<2,TYPEREEL> > > SubS2,Stot2;
// Vec<VecPointedValues<Interface<2,TYPEREEL> > > SubI2;
// Glob<2,TYPEREEL> Global2;

multiscale_geometry_mesh(n,S2, Inter2, process, CL2) ;
// assignation_materials_property_SST(n, S2, Inter2,process);
// assignation_materials_property_INTER(n, Inter2, S2, process);
// multiscale_operateurs(n,Stot2, SubS2,Inter2,SubI2, process,Global2);
// multiscale_iterate_latin(S2,SubS2, Inter2, SubI2, process,Global2,CL2) ;
// multiscale_iterate_incr(S2,SubS2, Inter2,SubI2, process, Global2, CL2,n) ;

#endif

// void calcul_cohesif(const XmlNode &n,Vec<SST> &S, Vec<INTER > &Inter, Param &process,  Vec<Boundary<DIMENSION,TYPEREEL> > &CL,Glob<DIMENSION,TYPEREEL> &Global);
}

#endif //FAKE_FONCTION_H

