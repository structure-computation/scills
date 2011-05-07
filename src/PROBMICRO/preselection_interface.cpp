//librairie Hugo
#include "containers/mat.h"
using namespace LMT;
using namespace std;
#include "mesh/mesh.h"

// fichiers de definition des variables 

#include "definition_INTER.h"
#include "definition_PARAM_MICRO_INTER.h"
#include "definition_SST.h"
#include "definition_micro.h"

//fcts utilisees dans ce .cpp
#include "criteres_preselection.h"

void preselection_interface(Vec<Interface<DIMENSION,TYPEREEL> > &Inter, Param_Micro<DIMENSION,TYPEREEL> &parammicro, Vec<Sst<DIMENSION,TYPEREEL> > &S){

   //-----------------------------------
   // critere de selection d'interfaces
   //-----------------------------------
   Vec<TYPEREEL> valinter;
   valinter.resize(Inter.size());
   apply_wi(Inter,calcul_sign_sigt(),valinter);
   Vec<unsigned> ind;
   ind=sort_with_index(valinter);
      // selection des n derniers termes
      Vec<unsigned> repI;
      repI.resize(parammicro.nbinterdeg);
      repI=ind[range(ind.size()-parammicro.nbinterdeg,ind.size())];
      for(unsigned i=0;i<repI.size();++i){
         Param_Micro<DIMENSION,TYPEREEL>::INTER interloc;
         interloc.num=repI[i];
         parammicro.Inter_select.push_back(interloc);
      }
      parammicro.num_Inter_select=repI;
}

