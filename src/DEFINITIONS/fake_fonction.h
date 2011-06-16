#ifndef FAKE_FONCTION_H
#define FAKE_FONCTION_H

// #include "definition_GLOB.h"
//#include "definition_fcts.h"
//fichiers contenant les differentes declarations des fonctions utilisees dans le programme principal
using namespace LMT;
using namespace std;



void fake_fonction(){
//pour le 3d
//     XmlNode n;
    Param process;
    
   
Vec<Sst<DIM,TYPEREEL> > S;
Vec<Interface<DIM,TYPEREEL> > Inter;
Vec<Boundary<DIM,TYPEREEL> > CL;
// Vec<VecPointedValues<Sst<3,TYPEREEL> > > SubS,Stot;
// Vec<VecPointedValues<Interface<3,TYPEREEL> > > SubI;
// Glob<3,TYPEREEL> Global;

// multiscale_geometry_mesh(n,S, Inter, process, CL) ;
/*assignation_materials_property_SST(n, S, Inter,process);
assignation_materials_property_INTER(n, Inter, S, process);
multiscale_operateurs(n,Stot, SubS,Inter,SubI, process,Global);
multiscale_iterate_latin(S,SubS, Inter, SubI, process,Global,CL) ;
multiscale_iterate_incr(S,SubS, Inter,SubI, process, Global, CL,n) ;*/

}

#endif //FAKE_FONCTION_H

