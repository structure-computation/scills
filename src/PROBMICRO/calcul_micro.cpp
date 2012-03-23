//librairie Hugo
#include "containers/mat.h"
using namespace LMT;
using namespace std;

#include "mesh/mesh.h"
#include "mesh/problem.h"
#include "displayparaview2.h"
#include <fstream>
#include <map>
#include "containers/gnuplot.h"

// fichiers de definition des variables et fcts pont entre les fichiers cpp
#include "definition_SST.h"
#include "definition_INTER.h"
#include "definition_PARAM_MICRO_INTER.h"
#include "Param.h"
#include "definition_PARAM_STRUCTURE.h"
#include "definition_PARAM_MULTI.h"
#include "definition_PARAM_LATIN.h"
#include "definition_PARAM_AFFICHAGE.h"
#include "definition_materials_property.h"
#include "definition_GLOB.h"
#include "Boundary.h"
#include "definition_fcts.h"

// fonctions speciales
#include "utilitaires.h"

//fichiers pour la strategie multiechelle
#include "read_prop_inter_micro.h"
#include "reconstruction_dep_contrainte.h"
#include "allocate.h"
#include "prelocalstage.h"

// definitions pour le probleme de degradations discretes
#include "definition_micro.h"
//fonction de pont entre les differents cpp utilises dans calcul_micro.cpp
#include "definition_fcts_cpp_micro.h"

//procedures utilisees dans ce .cpp
#include "modification_prop_interfaces.h"
#include "prerupture.h"
#include "calcul_EP_valeurs_moyennes.h"
#include "sortie_fichiers_micro.h"


//procedures rajoutee pour faire une pause
#include <limits> 



void calcul_micro(const XmlNode &n,Vec<Sst<DIM,TYPEREEL> > &S, Vec<Interface<DIM,TYPEREEL> > &Inter, Param &process,  Vec<Boundary<DIM,TYPEREEL> > &CL,Glob<DIM,TYPEREEL> &Global)
{

  cout << "*************************" <<endl;
  cout << "  Pretraitements micro "<<endl;
  cout << "*************************" <<endl;

      //lecture des parametres micro et proprietes d'interface
  Param_Micro<DIM,TYPEREEL> parammicro;
      read_micro(parammicro,n);
      //modification du comportement des interfaces interieures si besoin : contact defini a priori
      for(unsigned q=0;q<Inter.size();++q)
         Inter[q].parammicro = new PARAM_MICRO_INTER;
      
      Vec<InterCarac<DIM,TYPEREEL> > propintermicro;
      read_prop_inter_micro(propintermicro,n);
      
      // selection des interfaces potentiellement degradables et assignation des parametres materiaux


      // ligne pour le scale (a priori plus utilisee))
      modif_boxdegrad_scale(parammicro,process.structure->scale);
      modif_inter_degrad(Inter,S, propintermicro, parammicro);
      
      prerupture(Inter,parammicro.nb_inter_casse_ini);
      
      //calcul du volume total de la structure
      parammicro.volumetot=0.0;
      for(unsigned i=0;i<S.size();++i)
         parammicro.volumetot+=S[i].measure;
      
      //initialisation fichier contenant les pas de temps + iterations + numero du fichier paraview correspondant + le numero de l'interface cassee
      unsigned numfile=0;
      std::ofstream output_ptit,output_valmoy;
      create_file_output(output_ptit,output_valmoy, *process.affichage);

  cout << "*******************************"<<endl;
  cout << "  Procedure incrementale micro "<<endl;
  cout << "*******************************"<<endl;

      //allocations des differentes quantites pour le calcul multiechelle
      cout << "Allocations des quantites d'interfaces et SST" << endl;
      allocate_quantities(S,Inter,process,Global);


      for(unsigned pt=1;pt<=process.temps->nbpastemps;++pt)
      {
         process.temps->pt=pt;
      
   cout << "Appuyez sur entrée pour continuer...";
   cout << process.temps->modifCL << endl;
    cin.ignore( numeric_limits<streamsize>::max(), '\n' );


         // application du chargement aux interfaces
         cout << "*****************************************" << endl;
         cout << "* Modification des CL " << endl;
         assign_CL_values_space_time(Inter, CL,process);
         cout << endl;
         cout << "******************************************" << endl;
      



 



         // flag pour augmentation du chargement
         bool increaseload=0; //1 pour augmentation du chargement (nouveau piquet de temps pt)
         unsigned it=0; //compteur des iterations pour un piquet de temps pt
      
         //iterations pour un piquet de temps
         while(increaseload!=1)
         {
      
            //strategie iterative
            cout << "Resolution etat donne" << endl;
            multiscale_iterate(S, Inter,process, Global);
      
            //calcul energie potentielle pour cette configuration + sigma moyen - epsilon moyen
            parammicro.Epini=0.0;
            calcul_EP_sigmoy_epsmoy(parammicro.Epini,parammicro,S);
            //stockage des def et contraintes moyennes dans le fichier output_valmoy
            output_valmoy << parammicro.sigmoy <<"\t" << parammicro.epsmoy <<endl;
      
            //sauvegarde deformee initiale pour cet increment donne
            //affichage_maillage(S,Inter,process);
            //modif_name_paraviewfile(*process.affichage,numfile);
      
            //sortie des iterations , piquet de temps et interface degradee dans le fichier fptit
            output_ptit << numfile << " " << pt << " " << it << " " << parammicro.num_Inter_degradee << endl;
            numfile +=1;
      
            cout << "*****************************************" << endl;
            cout << "* Cassures alternatives d'interfaces " << endl;
            cout << "*****************************************" << endl;
      
            //---------------------------------------------------------------------------------------
            // critere de preselection d'interfaces parmi les interfaces degradables
            //---------------------------------------------------------------------------------------
            preselection_interface(Inter,parammicro,S);
      
            //-----------------------------------------------------------------------------------------
            // nouveaux calculs sur les n interfaces selectionnees en les faisant casser alternativement
            //-----------------------------------------------------------------------------------------
            Vec<unsigned> repI = parammicro.num_Inter_select;
      
            for(unsigned i=0;i<repI.size();++i)
            {
            //modification du comportement de chacune des interfaces : 1 contact les autres parfaites
            for(unsigned j=0;j<repI.size();++j)
            {
               if (repI[i]==repI[j]){Inter[repI[j]].comp="Contact";}
               else {Inter[repI[j]].comp="Parfait";}
            }
            //strategie iterative en repartant du calcul precedent (pour converger plus vite, il faudrait repartir du calcul avant la premiere degradation)
            multiscale_iterate(S, Inter,process,Global);
      
            // calcul de la nouvelle Ep
            parammicro.Inter_select[i].Ep=0.0;
            calcul_EP_sigmoy_epsmoy(parammicro.Inter_select[i].Ep,parammicro,S);
      
            //calcul du taux de restitution discret et stockage pour chaque interface testee
            parammicro.Inter_select[i].G=-1.0*(parammicro.Inter_select[i].Ep -parammicro.Epini)/Inter[repI[i]].measure;
            }
      
            //----------------------------------------------------------
            //selection de l'interface a casser et nouveau calcul ou non
            //----------------------------------------------------------
            selection_inter(Inter,parammicro,process,it,increaseload,numfile);
      
         }
      }
      
      for(unsigned q=0;q<Inter.size();++q)
         delete Inter[q].parammicro;
      
}


