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
#include "definition_PARAM.h"
#include "definition_PARAM_STRUCTURE.h"
#include "definition_PARAM_MULTI.h"
#include "definition_PARAM_LATIN.h"
#include "definition_PARAM_AFFICHAGE.h"
#include "definition_PARAM_TEMPS.h"
#include "definition_materials_property.h"
#include "definition_GLOB.h"
#include "definition_CL.h"
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

//procedures utilisees dans ce .cpp
#include "modification_prop_interfaces.h"
#include "prerupture.h"
//#include "calcul_EP_valeurs_moyennes.h"
#include "sortie_fichiers_micro.h"
#include "modif_etat_cohesif.h"
#include "affichage.h"
#include "affichage_mesh_INTER.h"

#include "modification_endommagements_max.h"

//procedures rajoutee pour faire une pause
#include <limits> 

void calcul_cohesif(const XmlNode &n,Vec<Sst<DIMENSION,TYPEREEL> > &S, Vec<Interface<DIMENSION,TYPEREEL> > &Inter, Param &process,  Vec<Boundary<DIMENSION,TYPEREEL> > &CL,Glob<DIMENSION,TYPEREEL> &Global)
{

  cout << "*************************" <<endl;
  cout << "  Pretraitements micro "<<endl;
  cout << "*************************" <<endl;

/*      //lecture des parametres micro et proprietes d'interface
      Param_Micro<DIMENSION,TYPEREEL> parammicro;
      read_micro(parammicro,n);
      //modification du comportement des interfaces interieures si besoin : contact defini a priori
      for(unsigned q=0;q<Inter.size();++q)
         Inter[q].parammicro = new PARAM_MICRO_INTER(Inter[q].side[0].nodeeq.size());
      
      Vec<InterCarac<DIMENSION,TYPEREEL> > propintermicro;
      read_prop_inter_micro(propintermicro,n);
      
      // selection des interfaces potentiellement degradables et assignation des parametres materiaux
      // ligne pour le scale (a priori plus utilisee))
      modif_boxdegrad_scale(parammicro,process.structure->scale);
      modif_inter_degrad(Inter,S, propintermicro, parammicro);
      
      prerupture(Inter,parammicro.nb_inter_casse_ini);*/
            
      //initialisation fichier contenant les pas de temps + iterations + numero du fichier paraview correspondant + le numero de l'interface cassee
      //unsigned numfile=0;
      std::ofstream output_ptit,output_valmoy;
      create_file_output(output_ptit,output_valmoy, *process.affichage);

  cout << "*******************************"<<endl;
  cout << "  Procedure incrementale  "<<endl;
  cout << "*******************************"<<endl;

      //allocations des differentes quantites pour le calcul multiechelle
      cout << "Allocations des quantites d'interfaces et SST" << endl;
      allocate_quantities(S,Inter,process,Global);

ofstream file_sort;
string nom_sortie="endosortie";
file_sort.open(nom_sortie.c_str(),ofstream::out );
      for(unsigned pt=1;pt<=process.temps->nbpastemps;++pt)
      {
         process.temps->pt=pt;
      
//         cout << "Appuyez sur entrée pour continuer...";
//         cout << process.temps->modifCL << endl;
//         cin.ignore( numeric_limits<streamsize>::max(), '\n' );

         // application du chargement aux interfaces
         cout << "*****************************************" << endl;
         cout << "* Modification des CL " << endl;
         assign_CL_values_space_time(Inter, CL,process);
         cout << endl;
         cout << "******************************************" << endl;

         //strategie iterative
         cout << "Resolution etat donne" << endl;
         multiscale_iterate(S, Inter,process, Global);
        
        //modification du comportement des interfaces cohesives si necessaire

   //     cout << "JUste avant modif etat cohesif" << endl ;
       // modif_etat_cohesif(Inter);
   //     cout << "JUste apres modif etat cohesif" << endl ;

         //modification des valeurs des endommagements et force maximaux au cours du temps

       
modification_endommagements_max(S);
 
file_sort <<  pt << " " << S[0].mesh.d[1] << " " << S[0].mesh.Yd << endl;
 
   //sortie des resultats de calculs dans fichier paraview : SST et INTER
        //sauvegarde deformee  pour cet increment donne
        //process.affichage->type_affichage="Sinterieur";
        affichage_resultats(S,process);
        modif_name_paraviewfile(*process.affichage,pt);

/*        process.affichage->type_affichage="Inter";
        affich_INTER_endommagement(Inter,process);
        modif_name_paraviewfile_inter(*process.affichage,pt);*/
    }
    
/*    for(unsigned q=0;q<Inter.size();++q)
         delete Inter[q].parammicro;*/
      
}


