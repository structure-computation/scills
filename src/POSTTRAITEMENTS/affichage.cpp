//librairie Hugo
#include "containers/mat.h"

#include "mesh/mesh.h"
#include "displayparaview2.h"
#include <fstream>
#include <map>
#include "containers/gnuplot.h"

// fichiers de definition des variables
#include "definition_PARAM_MICRO_INTER.h"
#include "definition_PARAM.h"
#include "definition_PARAM_AFFICHAGE.h"
#include "definition_SST_time.h"
#include "definition_INTER_time.h"

#include "mise_a_jour_quantites.h"

#include "affichage_mesh_SST.h"
#include "affichage_mesh_INTER.h"

#include "definition_PARAM_TEMPS.h"
#include "affichage_resultats_time.h"
#include "create_file_pvd.h"

#include "calculs_energies.h"
#include "extraction_quantites.h"

#include "containers/vecpointedvalues.h"
#include "mpi.h"

//#include "calcul_dep3d_defPG.h"
using namespace LMT;
using namespace std;


#include "affichage.h"
void fake_affichage() {
    XmlNode n;
    Param process;

    Vec<Sst<DIM,TYPEREEL> > S3;
    Vec<VecPointedValues<Sst<DIM,TYPEREEL> > > SP3;
    Vec<Interface<DIM,TYPEREEL> > Inter3;
    Vec<VecPointedValues<Interface<DIM,TYPEREEL> > > InterP3;
    
    affichage_maillage(SP3, InterP3,S3, process);
    affichage_resultats(SP3, process);
    affichage_depl_pt(SP3, process);
    affichage_var_inter(SP3,Inter3, process);
    affichage_inter_data(Inter3, S3, process);
    affichage_resultats_inter(InterP3, S3 , process);
    affichage_energie(SP3,Inter3, process);

}



/** \defgroup Post_Traitement
\brief Outils pour le post traitement du programme.
 
Plusieurs fonctions sont accessibles dans tous les fichiers, en incluant le fichier affichage.h.
- affichage_maillage() : affichage du maillage des sous-structures, bord ou intérieur et interfaces
- affichage_resultats() : affichage du maillage des sous-structures avec champs solutions choisis
- 
 
*/
/** \ingroup Post_Traitement
\brief Creation d'un fichier pvd regroupant les differentes solutions pour chaque pas de temps et lancement de paraview avec ce fichier
*/
void affichage_resultats_temps(Param &process) {
      create_file_pvd(process,"sst_");
      string namepvd = process.affichage->repertoire_save+"sst_"+process.affichage->name_data+".pvd";
      cout << "nom pvd : " << namepvd << endl;
      //string cmd = process.affichage->repertoire_save+"paraview --data="+namepvd;
      string cmd = "paraview";
      if (process.affichage->command_file=="") system(cmd.c_str());
};

/** \ingroup Post_Traitement
\brief Affichage des champs après calcul sur les interfaces pour tous les pas de temps
 
Possibilité de choisir une interface donnée ou toutes les interfaces.
*/
void affichage_inter_temps(Param &process) {
    create_file_pvd(process,"inter_");
    string namepvd = process.affichage->repertoire_save+"inter_"+process.affichage->name_data+".pvd";
    cout << "nom pvd : " << namepvd << endl;
    //string cmd = "paraview --data="+namepvd;
    string cmd = "paraview";
    if (process.affichage->command_file=="") system(cmd.c_str());
}




/**\ingroup Post_Traitement
 \brief Affichage du maillage des sous-structures ou des interfaces avec éclaté
 
 Pour les sous-structures, 
 - si le champ type_affichage est "Sinterieur", on appelle affich_SST()
 - si le champ type_affichage est "Sbord", on appelle affich_SST_edge()
 Pour les interfaces,
 - si le champ type_affichage est "Inter", on appelle affich_INTER()
 */
template <class TV3,class TV4, class TV1> 
void affichage_maillage(TV3 &S, TV4 &Inter,TV1 &Stot, Param &process){
    if (process.affichage->affich_mesh==1) {
      if (process.size==1 or process.rank>0){
        if (process.affichage->type_affichage=="Sinterieur" or process.affichage->type_affichage=="Sbord" or process.affichage->type_affichage=="all")
            affich_SST(S,process);
        else if (process.affichage->type_affichage=="Inter" or process.affichage->type_affichage=="all") {
            affich_INTER(Inter,Stot, process);
        } else
            cout << "erreur d'affichage" << endl;
      }
      if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
      if (process.size>1 and process.rank==0){create_file_pvtu(process,process.affichage->type_affichage); string cmd = "paraview"; if (process.affichage->command_file=="") system(cmd.c_str());}

    }
}

/**\ingroup Post_Traitement
 \brief Affichage des champs après calcul sur les sous-structures
 
 On appelle affich_SST_resultat() pour créer le fichier paraview de résultat pour chaque pas de temps.
 */
template <class TV3> 
void affichage_resultats(TV3 &S,  Param &process) {
    if (process.affichage->affich_resultat==1)
      if (process.size == 1 or process.rank > 0) 
        affich_SST_resultat_latin(S,process);
};


/**\ingroup Post_Traitement
 \brief Affichage des champs après calcul sur les interfaces
 
 On appelle affich_resultats_inter() pour créer le fichier paraview de résultat pour chaque pas de temps.
 */
template <class TV1,class TV4> 
void affichage_resultats_inter(TV4 &Inter, TV1 &S , Param &process) {
  if (process.affichage->affich_resultat==1)
      if (process.size == 1 or process.rank > 0) 
        affich_INTER_resultat(Inter,S,process);
};

/** \ingroup Post_Traitement
\brief Affichage des champs après calcul sur les interfaces
 
Possibilité de choisir une interface donnée ou toutes les interfaces.
*/
template <class TV1,class TV2> 
void affichage_inter_data(TV2 &Inter, TV1 &S, Param &process){
    if (process.affichage->affich_inter_data==1) {
        affich_inter_data_time(Inter,S,process);
    }
}


/** \ingroup Post_Traitement
\brief Affichage de l'evolution du déplacement d'un point donné par ses coordonnées 
*/
template <class TV3> 
void affichage_depl_pt(TV3 &S, Param &process){
    if(process.affichage->affich_depl_pt==1) extraction_depl_pt(S, process);
}

/** \ingroup Post_Traitement
\brief Affichage de l'évolution de l'énergie dissipée ou de l'énergie imposée au cours du temps à partir des quantités chapeaux ou des quantités n de l'interface
 
Selon les parametres du champ AFFICHAGE::param_ener, on sélectionne le type d'énergie et les quantités retenues.
0 - 0 : energie dissipee sur les quantites chapeau
0 - 1 : energie dissipee sur les quantites n
1 - 0 : energie imposee sur les quantites chapeau
1 - 1 : energie imposee sur les quantites n

Faux en MPI pour certaine fonction qui necessite d avoir les W des deux cotes et ils sont pas transferes pour le posttraitement : est-ce utile de le faire ? non pour le moment
*/
template <class TV3,class TV2> 
void affichage_energie(TV3 &S,TV2 &Inter, Param &process){
    Vec<double> energie,temp;
    energie.resize(process.temps->nbpastemps + 1);temp.resize(process.temps->nbpastemps + 1);energie.set(0.);temp.set(0.);
    if(process.affichage->param_ener[0]==0 and process.affichage->param_ener[1]==0) {
        if (process.rank>0 or process.size==1) calcul_ener_dissi_chap(S,Inter,energie,process);
        if (process.size>1) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.rank==0) cout << "Dissipation (c) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==0 and process.affichage->param_ener[1]==1) {
      if (process.rank>0 or process.size==1) calcul_ener_dissi_lin(S,Inter,energie,process);
      if (process.size>1) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
      if (process.rank==0) cout << "Dissipation (n) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==1 and process.affichage->param_ener[1]==0) {
      if (process.rank>0 or process.size==1) calcul_ener_imp_chap(S,Inter,energie,process);
      if (process.size>1) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
      if (process.rank==0) cout << "Energie imposee (c) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==1 and process.affichage->param_ener[1]==1) {
      if (process.rank>0 or process.size==1) calcul_ener_imp_lin(S,Inter,energie,process);
      if (process.size>1) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
      if (process.rank==0) cout << "Energie imposee (n) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==2 and process.affichage->param_ener[1]==0) {
      if (process.rank>0 or process.size==1) calcul_Ft2_chap(S,Inter,energie,process);
      if (process.size>1) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
      if (process.rank==0) cout << "Effort tangent carre (c) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==2 and process.affichage->param_ener[1]==1) {
      if (process.rank>0 or process.size==1) calcul_Ft2_lin(S,Inter,energie,process);
      if (process.size>1) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
      if (process.rank==0) cout << "Effort tangent carre (n) faux en mpi : " << energie << endl;
    } else if(process.affichage->param_ener[0]==3 and process.affichage->param_ener[1]==0) {
      if (process.rank>0 or process.size==1) calcul_Fn_chap(S,Inter,energie,process);
      if (process.size>1) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
      if (process.rank==0) cout << "Effort normal (c) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==3 and process.affichage->param_ener[1]==1) {
      if (process.rank>0 or process.size==1) calcul_Fn_lin(S,Inter,energie,process);
      if (process.size>1) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
      if (process.rank==0) cout << "Effort normal (n) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==4 and process.affichage->param_ener[1]==0) {
      if (process.rank>0 or process.size==1) calcul_Un_chap(S,Inter,energie,process);
      if (process.size>1) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
      if (process.rank==0) cout << "Deplacement normal (c) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==4 and process.affichage->param_ener[1]==1) {
      if (process.rank>0 or process.size==1) calcul_Un_lin(S,Inter,energie,process);
      if (process.size>1) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
      if (process.rank==0) cout << "Deplacement normal (n) faux en mpi  : " << energie << endl;
    } else if(process.affichage->param_ener[0]==5 and process.affichage->param_ener[1]==0) {
      if (process.rank>0 or process.size==1) calcul_Ut_chap(S,Inter,energie,process);
      if (process.size>1) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
      if (process.rank==0) cout << "Deplacement tangent (c) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==5 and process.affichage->param_ener[1]==1) {
      if (process.rank>0 or process.size==1) calcul_Ut_lin(S,Inter,energie,process);
      if (process.size>1) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
      if (process.rank==0) cout << "Deplacement tanget (n) faux en mpi  : " << energie << endl;
    }  else if(process.affichage->param_ener[0]==6) {
        if (process.rank>0 or process.size==1) calcul_energie_elastique(S,Inter,energie,process);
        if (process.size>1) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.rank==0) cout << "Energie elastique : " << energie << endl;
    } else {
        cout << "Mauvais choix d'energie" << endl;
        assert(0);
    }
    if (process.affichage->command_file=="" and process.rank==0){
      GnuPlot gp;
      gp.plot(energie);
      gp.wait();
    }
}


