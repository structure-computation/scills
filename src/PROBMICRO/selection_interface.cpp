
//librairie Hugo
#include "containers/mat.h"
using namespace LMT;
using namespace std;
#include "mesh/mesh.h"

// fichiers de definition des variables 

#include "definition_INTER.h"
#include "definition_PARAM_MICRO_INTER.h"
#include "definition_PARAM.h"
#include "definition_PARAM_AFFICHAGE.h"
#include "definition_micro.h"

//fcts utilisees dans ce .cpp
#include "selection_inter.h"

void selection_inter(Vec<Interface<DIMENSION,TYPEREEL> > &Inter, Param_Micro<DIMENSION,TYPEREEL> &parammicro, Param &process, unsigned &it, bool &increaseload,unsigned &numfile){
  //--------------------------------------
  //selection de l'interface a casser
  //--------------------------------------
            
  Vec<TYPEREEL> critere;
  critere.resize(parammicro.nbinterdeg);
  cout << "****************************************" << endl;
  cout << "Diagnostic sur chaque interface potentielle :" << endl;           
  cout << "****************************************"<< endl;

  Vec<unsigned> repI=parammicro.num_Inter_select;
   for(unsigned i=0;i<repI.size();++i){
      critere[i]=parammicro.Inter_select[i].G/Inter[repI[i]].parammicro->Gcrit;
      cout << "Interface numero :" << "----- \t \t " << i << endl;
      cout << "Energie potentielle :  \t \t  \t" << parammicro.Inter_select[i].G << endl;
      cout << "% de l'energie critique : \t \t " << 100*critere[i] << endl;
      }
   cout << "*" << endl;
            
   if (find(critere,_1>=1.0) ){
      //recherche de l'interface a critere max
      unsigned interface_select=0;
      TYPEREEL val=max(critere);
      for(unsigned k=0;k<critere.size();++k)
         if(critere[k]==val){interface_select=k;break;}
         
      // modification finale du comportement des interfaces selectionnees
      for(unsigned i=0;i<repI.size();++i){
         if (i==interface_select){Inter[repI[i]].comp="Contact"; 
            cout << "*****  Detection de rupture  *****";
            cout<< "Interface numero " << repI[i] << " degradee " << endl;
            cout << "Valeur du critere actif : " << 100 * critere[i];}
         else {   Inter[repI[i]].comp="Parfait";     }
      }
      //affichage et stockage des interfaces testees
      parammicro.num_Inter_degradee=repI[interface_select];
      stockage_interfaces_testees(Inter,*process.affichage,numfile-1,repI,parammicro.num_Inter_degradee);
      //augmentation du flag d'iterations
      it+=1;
      }
   else {
      // les interfaces restent parfaites pour ce chargement 
            cout << "*****  Aucune rupture detectee  *****   (" << 100*max(critere) << " % ) du critere" << endl;
      for(unsigned i=0;i<repI.size();++i)
         {
         Inter[repI[i]].comp="Parfait";} 
         //affichage et stockage des interfaces testees
         parammicro.num_Inter_degradee=-1;
         stockage_interfaces_testees(Inter,*process.affichage,numfile-1,repI,parammicro.num_Inter_degradee);
         // augmentation chargement
         increaseload=1;
         
         process.critere_max=max(critere)-0.02;
         process.modifCL=1;
      }

}

