//librairie Hugo
#include "../../LMT/include/containers/mat.h"
#include "../../LMT/include/containers/vecpointedvalues.h"
#include "../../LMT/include/mesh/mesh.h"
#include "../../LMT/include/containers/gnuplot.h"

#include "displayparaview2.h"

#include <fstream>
#include <map>

// fichiers de definition des variables
#include "../DEFINITIONS/Process.h"
#include "../DEFINITIONS/SavingData.h"
#include "../DEFINITIONS/Sst.h"
#include "../DEFINITIONS/Interface.h"
#include "../DEFINITIONS/TimeData.h"

#include "../ITERATIONS/manipulate_quantities.h"

#include "affichage_mesh_SST.h"
#include "affichage_mesh_INTER.h"

#include "affichage_resultats_time.h"
#include "create_file_pvd.h"

#include "calculs_energies.h"
#include "extraction_quantites.h"

#include "mpi.h"

//#include "calcul_dep3d_defPG.h"
using namespace LMT;


#include "affichage.h"
void fake_affichage() {
    Process process;
    DataUser data_user;
    Vec<Sst> S3;
    Vec<VecPointedValues<Sst> > SP3;
    Vec<Interface> Inter3;
    Vec<VecPointedValues<Interface> > InterP3;
    
    affichage_maillage(SP3, InterP3,S3, process, data_user);
    affichage_resultats(SP3, process, data_user);
    affichage_depl_pt(SP3, process);
    affichage_var_inter(SP3,Inter3, process);
    //affichage_inter_data(Inter3, S3, process);
    affichage_resultats_inter(InterP3, S3 , process, data_user);
    affichage_energie(SP3,Inter3, process, data_user);
    create_pvd_geometry(SP3, S3, Inter3,process);
    create_pvd_results(SP3, S3, Inter3,process);
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
void affichage_resultats_temps(Process &process) {
      create_file_pvd(process,"sst_");
      Sc2String namepvd = process.affichage->repertoire_save+"sst_"+process.affichage->name_data+".pvd";
      std::cout << "nom pvd : " << namepvd << endl;
      //Sc2String cmd = process.affichage->repertoire_save+"paraview --data="+namepvd;
      Sc2String cmd = "paraview";
      if (process.affichage->command_file=="") int tmp=system(cmd.c_str());
};

/** \ingroup Post_Traitement
\brief Affichage des champs après calcul sur les interfaces pour tous les pas de temps
 
Possibilité de choisir une interface donnée ou toutes les interfaces.
*/
void affichage_inter_temps(Process &process) {
    create_file_pvd(process,"inter_");
    Sc2String namepvd = process.affichage->repertoire_save+"inter_"+process.affichage->name_data+".pvd";
    std::cout << "nom pvd : " << namepvd << endl;
    //Sc2String cmd = "paraview --data="+namepvd;
    Sc2String cmd = "paraview";
    if (process.affichage->command_file=="") int tmp=system(cmd.c_str());
}

template <class TV> 
void create_pvd_interfaces(TV &Inter, Process &process){
   //Sc2String save_directory=process.affichage->repertoire_save+"results/inter/";
    
   //creation du nom et du fichier pvd
   std::ostringstream spvd;
   
   spvd<<process.affichage->repertoire_save<<"results/Geometry_inter_proc"<< process.parallelisation->rank<<".pvd";
   Sc2String namepvd( spvd.str() );
   
   std::ofstream file_pvd;
   file_pvd.open( (namepvd).c_str(), std::ofstream::out);
   
   //ecriture entete
   file_pvd << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" \n\t byte_order=\"LittleEndian\" \n\t compressor=\"vtkZLibDataCompressor\" >\n <Collection> " <<std::endl;

  
   //ecriture des lignes correspondant aux fichier a lire
    for(unsigned i=0;i<Inter.size();++i){
        //nom du fichier a lire pour le piquet de temps et le processeur considere
        std::ostringstream ss;
        ss << process.affichage->repertoire_save<<"results/Geometry_inter/proc_"<< process.parallelisation->rank << "_Inter_"<< Inter[i].id <<  ".vtu";
        Sc2String nom1( ss.str() );
        //designation du dataset sous paraview (valeur du piquet de temps, processeur)    
        std::ostringstream ss2;
        ss2 << " <DataSet timestep=\"1\"" << " group=\" " <<i << "\" part=\" \"\n\t\t file=\"" << nom1 << "\"/>\n" ;
        Sc2String ligne1(ss2.str()); 
        file_pvd << ligne1;
    }          
   //fin du fichier pvd
   file_pvd << " </Collection> \n </VTKFile>";   
}


/**
Fonction permettant de sortir des fichiers pvd pour la mise en donnée
*/
template <class TV1,class TV2, class TV3>
void create_pvd_geometry(TV1 &S, TV3 &SS, TV2 &Inter,Process &process) {
  create_pvd_geometry_inter(S, SS, Inter,process);
  create_pvd_geometry_sst(S, SS, Inter,process);
  //create_pvd_geometry_sst_by_group(S, SS, Inter,process);
}
  
/**
Fonction permettant de sortir des fichiers pvd pour la mise en donnée
*/
template <class TV1,class TV2, class TV3>
void create_pvd_geometry_inter(TV1 &S, TV3 &SS, TV2 &Inter,Process &process) {

     ///preparation des noms et des repertoires pour ecriture des resultats
    Sc2String save_directory=process.affichage->repertoire_save+"Geometry/";
    Sc2String base_filename= save_directory;
    base_filename << "inter_proc_" ; 
    Sc2String namefile = base_filename;
    namefile << process.parallelisation->rank << ".pvd";
    
    ofstream os( namefile.c_str() );
    if (process.parallelisation->is_master_cpu())  os <<"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" \n\t byte_order=\"LittleEndian\" \n\t compressor=\"vtkZLibDataCompressor\" >\n <Collection> " <<std::endl;

    if (process.parallelisation->is_multi_cpu()) 
       process.parallelisation->synchronisation();
    
    unsigned data=process.affichage->side;
    for(unsigned i=0;i<Inter.size();++i) {
            if (SS[Inter[i].vois[data*2]].num_proc==process.parallelisation->rank){		
      ostringstream ss, ss2;
      ss<<"./inter/inter_id_"<<Inter[i].id<<".vtu";
      Sc2String namefile_entity(ss.str());
      ss2 << " <DataSet timestep=\"1\"" << " group=\" " <<Inter[i].id << "\" part=\" \"\n\t\t file=\"" << namefile_entity << "\"/>\n" ;
      Sc2String ligne1(ss2.str()); 
      os << ligne1;
            }
    }    
    os.close();
    
    ///Fin des écritures
    if (process.parallelisation->is_multi_cpu()) {
      process.parallelisation->synchronisation();
           
      if (process.parallelisation->is_master_cpu()){
	///concaténation des fichiers
	Sc2String cmd;
	cmd << "cat ";
	for (unsigned i=0;i<process.parallelisation->size;i++)
	  cmd << base_filename << i << ".pvd ";
	cmd << "> " << save_directory;
 	cmd << "geometry_inter.pvd" ; 
	int tmp=system(cmd.c_str());

	Sc2String namefile2 = process.affichage->repertoire_save+"Geometry/geometry_inter.pvd";
	fstream os2;
	os2.open ( namefile2.c_str(), fstream::in | fstream::out | fstream::ate );
	os2 << " </Collection> \n </VTKFile>";
	os2.close();
	Sc2String cmd2;
	cmd2 << "rm ";
	for (unsigned i=0;i<process.parallelisation->size;i++)
	  cmd2 << base_filename << i << ".pvd ";
	tmp=system(cmd2.c_str());
      }
      	process.parallelisation->synchronisation();

    } else
    {
        
	Sc2String cmd;
	cmd << "mv " << base_filename << "0" << ".pvd ";
	cmd << " " << save_directory;	
	cmd << "geometry_inter.pvd" ; 
	int tmp=system(cmd.c_str());  
	
	Sc2String namefile2 = process.affichage->repertoire_save+"Geometry/geometry_inter.pvd";
	fstream os2;
	os2.open ( namefile2.c_str(),fstream::in | fstream::out | fstream::ate  );
	os2 << " </Collection> \n </VTKFile>";
	os2.close();
    }
    
    
}


/**
Fonction permettant de sortir des fichiers pvd pour la mise en donnée
*/
template <class TV1,class TV2, class TV3>
void create_pvd_geometry_sst(TV1 &S, TV3 &SS, TV2 &Inter,Process &process) {

     ///preparation des noms et des repertoires pour ecriture des resultats
    Sc2String save_directory=process.affichage->repertoire_save+"Geometry/";
    Sc2String base_filename= save_directory;
    base_filename << "sst_proc_" ; 
    Sc2String namefile = base_filename;
    namefile << process.parallelisation->rank << ".pvd";
    
    ofstream os( namefile.c_str() );
    if (process.parallelisation->is_master_cpu())  os <<"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" \n\t byte_order=\"LittleEndian\" \n\t compressor=\"vtkZLibDataCompressor\" >\n <Collection> " <<std::endl;

    //     std::cout << S.size() << endl;
    if (process.parallelisation->is_multi_cpu()) 
       process.parallelisation->synchronisation();
    
    for(unsigned i=0;i<S.size();i++) {
      ostringstream ss, ss2;
      ss<<"./sst/sst_id_"<<S[i].id<<".vtu";
      Sc2String namefile_entity(ss.str());
      ss2 << " <DataSet timestep=\"1\"" << " group=\" " <<S[i].id << "\" part=\" \"\n\t\t file=\"" << namefile_entity << "\"/>\n" ;
      Sc2String ligne1(ss2.str()); 
      os << ligne1;
    }
    
    os.close();
    ///Fin des écritures
    if (process.parallelisation->is_multi_cpu()) {
      process.parallelisation->synchronisation();
           
      if (process.parallelisation->is_master_cpu()){
	///concaténation des fichiers
	Sc2String cmd;
	cmd << "cat ";
	for (unsigned i=0;i<process.parallelisation->size;i++)
	  cmd << base_filename << i << ".pvd ";
	cmd << "> " << save_directory;
 	cmd << "geometry_sst.pvd" ; 
	int tmp=system(cmd.c_str());

	Sc2String namefile2 = process.affichage->repertoire_save+"Geometry/geometry_sst.pvd";
	fstream os2;
	os2.open ( namefile2.c_str(), fstream::in | fstream::out | fstream::ate );
	os2 << " </Collection> \n </VTKFile>";
	os2.close();
	Sc2String cmd2;
	cmd2 << "rm ";
	for (unsigned i=0;i<process.parallelisation->size;i++)
	  cmd2 << base_filename << i << ".pvd ";
	tmp=system(cmd2.c_str());
// 	std::cout << cmd2 << endl;
      }
      	process.parallelisation->synchronisation();

    } else
    {
        
	Sc2String cmd;
	cmd << "mv " << base_filename << "0" << ".pvd ";
	cmd << " " << save_directory;	
	cmd << "geometry_sst.pvd" ; 
	int tmp=system(cmd.c_str());  
	
	Sc2String namefile2 = process.affichage->repertoire_save+"Geometry/geometry_sst.pvd";
	fstream os2;
	os2.open ( namefile2.c_str(),fstream::in | fstream::out | fstream::ate  );
	os2 << " </Collection> \n </VTKFile>";
	os2.close();
    }
    
    
}

/**
Fonction permettant de sortir des fichiers pvd pour la mise en donnée
*/
template <class TV1,class TV2, class TV3>
void create_pvd_geometry_sst_by_group(TV1 &S, TV3 &SS, TV2 &Inter,Process &process) {
/*
     ///preparation des noms et des repertoires pour ecriture des resultats
    Sc2String save_directory=process.affichage->repertoire_save+"Geometry/";
    Sc2String base_filename= save_directory;
    base_filename << "sst_proc_" ; 
    for(unsigned j=0;j<process.sst_materials->size();j++) {
      Sc2String namefile = base_filename;
      namefile << process.parallelisation->rank << "_group_"<<(*process.sst_materials[j]).id<<".pvd";
    
      ofstream os( namefile.c_str() );
      if (process.parallelisation->is_master_cpu()) os <<"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" \n\t byte_order=\"LittleEndian\" \n\t compressor=\"vtkZLibDataCompressor\" >\n <Collection> " <<std::endl;
    
      if (process.parallelisation->is_multi_cpu()) 
       process.parallelisation->synchronisation();
    
      for(unsigned i=0;i<S.size();i++) {
	ostringstream ss, ss2;
	ss<<"./sst/proc_"<< process.parallelisation->rank << "_sst_id_"<<S[i].id<<".vtu";
	Sc2String namefile_entity(ss.str());
	ss2 << " <DataSet timestep=\"1\"" << " group=\" " <<S[i].id << "\" part=\" \"\n\t\t file=\"" << namefile_entity << "\"/>\n" ;
	Sc2String ligne1(ss2.str()); 
	if(S[i].id_material==(*process.sst_materials[j]).id){
	  os << ligne1;
	}
	
      }

      os.close();
    }
    */
   
    
/*    
    
    ///Fin des écritures
    if (process.parallelisation->is_multi_cpu()) {
      process.parallelisation->synchronisation();
      
      
      for(unsigned j=0;j<process.sst_materials->size();j++) {
	Sc2String namefile = base_filename;
	namefile << process.parallelisation->rank << "_group_"<<(*process.sst_materials[j]).id<<".pvd";
      
	ofstream os( namefile.c_str() );
	
	if (process.parallelisation->is_master_cpu()){
	  
	  
	  ///concaténation des fichiers
	  Sc2String cmd;
	  cmd << "cat ";
	  for (unsigned i=0;i<process.parallelisation->size;i++)
	    cmd << base_filename << i << "_group_"<<(*process.sst_materials[j]).id << ".pvd ";
	  cmd << "> " << save_directory;
	  cmd << "geometry_group_ssts_"<< j<<".pvd" ; 
	  int tmp=system(cmd.c_str());

	  Sc2String namefile2 = process.affichage->repertoire_save+"Geometry/geometry_sst_group_"<<(*process.sst_materials[j]).id<<".pvd";
	  fstream os2;
	  os2.open ( namefile2.c_str(), fstream::in | fstream::out | fstream::ate );
	  os2 << " </Collection> \n </VTKFile>";
	  os2.close();
	  Sc2String cmd2;
	  cmd2 << "rm ";
	  for (unsigned i=0;i<process.parallelisation->size;i++)
	    cmd2 << base_filename << i << "_group_"<<(*process.sst_materials[j]).id << ".pvd ";
	  tmp=system(cmd2.c_str());
  // 	std::cout << cmd2 << endl;
	}
	  process.parallelisation->synchronisation();
      }

    } else
    {
       
       for(unsigned j=0;j<process.sst_materials->size();j++) { 
	Sc2String cmd;
	cmd << "mv " << base_filename << "0" << "_group_"<<(*process.sst_materials[j]).id << ".pvd ";
	cmd << " " << save_directory;	
	cmd << "geometry_sst-group_"<< (*process.sst_materials[j]).id << ".pvd" ; 
	int tmp=system(cmd.c_str());  
	
	Sc2String namefile2 = process.affichage->repertoire_save+"Geometry/geometry_sst_group_"<< (*process.sst_materials[j]).id <<".pvd";
	fstream os2;
	os2.open ( namefile2.c_str(),fstream::in | fstream::out | fstream::ate  );
	os2 << " </Collection> \n </VTKFile>";
	os2.close();
       }
    }
    */
    
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
void affichage_maillage(TV3 &S, TV4 &Inter,TV1 &Stot, Process &process, DataUser &data_user){
    PRINT(process.affichage->type_affichage);
    PRINT(process.affichage->affich_mesh);
    PRINT(process.parallelisation->is_local_cpu());
    process.affichage->affich_mesh=1;
    if (process.affichage->affich_mesh==1) {
        if (process.parallelisation->is_local_cpu()){
            std::cout << "type " << process.affichage->type_affichage << std::endl;
            if (process.affichage->type_affichage=="Sinterieur" or process.affichage->type_affichage=="Sbord" ){
                affich_SST(S,process);
            }else if (process.affichage->type_affichage=="Inter" ) {
                affich_INTER(Inter,Stot, process);
            }else if (process.affichage->type_affichage=="all") { 
                process.affichage->type_affichage="Sbord";
                PRINT(process.affichage->type_affichage);
                affich_SST(S,process);
		//create_pvd_by_material(S,process);
                process.affichage->type_affichage="Inter";
                PRINT(process.affichage->type_affichage);
                affich_INTER(Inter,Stot, process);
		 //create_pvd_inter_by_behaviour(Inter,S,process);
		 //create_pvd_inter_edge(Inter,S,process);
            } else {
                std::cout << "erreur d'affichage" << endl;
            }
        }
        process.parallelisation->synchronisation();
        if (not process.parallelisation->is_local_cpu()){
	  //create_file_pvtu(process,"sst_"); 
	//Sc2String cmd = "paraview"; if (process.affichage->command_file=="") int tmp=system(cmd.c_str());
	  
	}
        //if (not process.parallelisation->is_local_cpu()) create_file_pvd_geometry(process,data_user,"Geometry_sst");
        //if (not process.parallelisation->is_local_cpu()) create_file_pvd_geometry(process,data_user,"Geometry_inter");

    }
    //assert(0);
}


/**
Fonction permettant de sortir des fichiers pvd pour la mise en donnée
*/
template <class TV1,class TV2, class TV3>
void create_pvd_results(TV1 &S, TV3 &SS, TV2 &Inter,Process &process) {
  create_pvd_results_sst(S, SS, Inter,process);
  create_pvd_results_sst_skin(S, SS, Inter,process);
  create_pvd_results_inter(S, SS, Inter,process);
}
  

/**
Fonction permettant de sortir des fichiers pvd pour les résultats
*/
template <class TV1,class TV2, class TV3>
void create_pvd_results_sst(TV1 &S, TV3 &SS, TV2 &Inter,Process &process) {
    //for (unsigned j=0;j<2;j++) {
     ///preparation des noms et des repertoires pour ecriture des resultats
    Sc2String save_directory=process.affichage->repertoire_save+"results/";
    Sc2String base_filename= save_directory;
    if(process.multiresolution->nb_calculs>1)
        base_filename<<"resolution_"<<process.multiresolution->m<<"_";
    base_filename << "proc_" ; 
    Sc2String namefile = base_filename;
    namefile << process.parallelisation->rank << ".pvd";
        
    ofstream os( namefile.c_str() );
    if (process.parallelisation->is_master_cpu())  os <<"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" \n\t byte_order=\"LittleEndian\" \n\t compressor=\"vtkZLibDataCompressor\" >\n <Collection> " <<std::endl;

    //     std::cout << S.size() << endl;
    if (process.parallelisation->is_multi_cpu()) 
       process.parallelisation->synchronisation();
    
    Sc2String local_filename="./sst_bulk/result";
//     if(j==0) local_filename="./sst_bulk/result_";
//     else local_filename="./sst_skin/result_"; 
    if(process.multiresolution->nb_calculs>1)
        local_filename<<"_resolution_"<<process.multiresolution->m;
   
    
    for(unsigned i=0;i<S.size();i++) {
      for(unsigned k=1;k<S[i].t_post.size();k++){
      ostringstream ss, ss2;
      ss<<local_filename << "_sst_id_"<<S[i].id<<"_time_"<<k<<".vtu";
      Sc2String namefile_entity(ss.str());
      ss2 << " <DataSet timestep=\" "<< k<< " \" group=\" " <<S[i].id << "\" part=\" \"\n\t\t file=\"" << namefile_entity << "\"/>\n" ;
      Sc2String ligne1(ss2.str()); 
      os << ligne1;
      }
    }
    
    os.close();
    ///Fin des écritures
    if (process.parallelisation->is_multi_cpu()) {
      process.parallelisation->synchronisation();
           
      if (process.parallelisation->is_master_cpu()){
	///concaténation des fichiers
	Sc2String cmd;
	cmd << "cat ";
	for (unsigned i=0;i<process.parallelisation->size;i++)
	  cmd << base_filename << i << ".pvd ";
	cmd << "> " << save_directory;
	Sc2String name_file_output="results_sst";
	if(process.multiresolution->nb_calculs>1)
	  name_file_output<<"_resolution_"<<process.multiresolution->m;
	name_file_output<<".pvd";
	cmd << name_file_output ;
	int tmp=system(cmd.c_str());

	Sc2String namefile2 = process.affichage->repertoire_save+"results/"+name_file_output;
	fstream os2;
	os2.open ( namefile2.c_str(), fstream::in | fstream::out | fstream::ate );
	os2 << " </Collection> \n </VTKFile>";
	os2.close();
	Sc2String cmd2;
	cmd2 << "rm ";
	for (unsigned i=0;i<process.parallelisation->size;i++)
	  cmd2 << base_filename << i << ".pvd ";
	tmp=system(cmd2.c_str());
// 	std::cout << cmd2 << endl;
      }
      process.parallelisation->synchronisation();
    } else
    {
        
	Sc2String cmd;
	cmd << "mv " << base_filename << "0" << ".pvd ";
	cmd << " " << save_directory;
	Sc2String name_file_output="results_sst";
	if(process.multiresolution->nb_calculs>1)
	  name_file_output<<"_resolution_"<<process.multiresolution->m;
	name_file_output<<".pvd";
	cmd << name_file_output ;
	int tmp=system(cmd.c_str());

	Sc2String namefile2 = process.affichage->repertoire_save+"results/"+name_file_output;
	fstream os2;
	os2.open ( namefile2.c_str(),fstream::in | fstream::out | fstream::ate  );
	os2 << " </Collection> \n </VTKFile>";
	os2.close();
     }
//     }
    
}

/**
Fonction permettant de sortir des fichiers pvd pour les résultats
*/
template <class TV1,class TV2, class TV3>
void create_pvd_results_sst_skin(TV1 &S, TV3 &SS, TV2 &Inter,Process &process) {
    //for (unsigned j=0;j<2;j++) {
     ///preparation des noms et des repertoires pour ecriture des resultats
    Sc2String save_directory=process.affichage->repertoire_save+"results/";
    Sc2String base_filename= save_directory;
    if(process.multiresolution->nb_calculs>1)
        base_filename<<"resolution_"<<process.multiresolution->m<<"_";
    base_filename << "proc_" ; 
    Sc2String namefile = base_filename;
    namefile << process.parallelisation->rank << ".pvd";
        
    ofstream os( namefile.c_str() );
    if (process.parallelisation->is_master_cpu())  os <<"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" \n\t byte_order=\"LittleEndian\" \n\t compressor=\"vtkZLibDataCompressor\" >\n <Collection> " <<std::endl;

    //     std::cout << S.size() << endl;
    if (process.parallelisation->is_multi_cpu()) 
       process.parallelisation->synchronisation();
    
    Sc2String local_filename="./sst_skin/result";
//     if(j==0) local_filename="./sst_bulk/result_";
//     else local_filename="./sst_skin/result_";
    if(process.multiresolution->nb_calculs>1)
        local_filename<<"_resolution_"<<process.multiresolution->m;
   
    
    for(unsigned i=0;i<S.size();i++) {
      for(unsigned k=1;k<S[i].t_post.size();k++){
      ostringstream ss, ss2;
      ss<<local_filename<<"_sst_id_"<<S[i].id<<"_time_"<<k<<".vtu";
      Sc2String namefile_entity(ss.str());
      ss2 << " <DataSet timestep=\" "<< k<< " \" group=\" " <<S[i].id << "\" part=\" \"\n\t\t file=\"" << namefile_entity << "\"/>\n" ;
      Sc2String ligne1(ss2.str()); 
      os << ligne1;
      }
    }
    
    os.close();
    ///Fin des écritures
    if (process.parallelisation->is_multi_cpu()) {
      process.parallelisation->synchronisation();
           
      if (process.parallelisation->is_master_cpu()){
	///concaténation des fichiers
	Sc2String cmd;
	cmd << "cat ";
	for (unsigned i=0;i<process.parallelisation->size;i++)
	  cmd << base_filename << i << ".pvd ";
	cmd << "> " << save_directory;
	Sc2String name_file_output="results_sst_skin";
	if(process.multiresolution->nb_calculs>1)
	  name_file_output<<"_resolution_"<<process.multiresolution->m;
	name_file_output<<".pvd";
	cmd << name_file_output ;
	int tmp=system(cmd.c_str());

	Sc2String namefile2 = process.affichage->repertoire_save+"results/"+name_file_output;
	fstream os2;
	os2.open ( namefile2.c_str(), fstream::in | fstream::out | fstream::ate );
	os2 << " </Collection> \n </VTKFile>";
	os2.close();
	Sc2String cmd2;
	cmd2 << "rm ";
	for (unsigned i=0;i<process.parallelisation->size;i++)
	  cmd2 << base_filename << i << ".pvd ";
	tmp=system(cmd2.c_str());
// 	std::cout << cmd2 << endl;
      }
      process.parallelisation->synchronisation();
    } else
    {
        
	Sc2String cmd;
	cmd << "mv " << base_filename << "0" << ".pvd ";
	cmd << " " << save_directory;
	Sc2String name_file_output="results_sst_skin";
	if(process.multiresolution->nb_calculs>1)
	  name_file_output<<"_resolution_"<<process.multiresolution->m;
	name_file_output<<".pvd";
	cmd << name_file_output ;
	int tmp=system(cmd.c_str());

	Sc2String namefile2 = process.affichage->repertoire_save+"results/"+name_file_output;
// 	if(j==0) namefile2 = process.affichage->repertoire_save+"results/results_sst.pvd";
// 	else namefile2 = process.affichage->repertoire_save+"results/results_sst_skin.pvd";
	fstream os2;
	os2.open ( namefile2.c_str(),fstream::in | fstream::out | fstream::ate  );
	os2 << " </Collection> \n </VTKFile>";
	os2.close();
     }
//     }
    
}


/**
Fonction permettant de sortir des fichiers pvd pour les résultats
*/
template <class TV1,class TV2, class TV3>
void create_pvd_results_inter(TV1 &S, TV3 &SS, TV2 &Inter,Process &process) {
     ///preparation des noms et des repertoires pour ecriture des resultats
    Sc2String save_directory=process.affichage->repertoire_save+"results/";
    Sc2String base_filename= save_directory;
    if(process.multiresolution->nb_calculs>1)
        base_filename<<"resolution_"<<process.multiresolution->m<<"_";
    base_filename << "proc_" ; 
    Sc2String namefile = base_filename;
    namefile << process.parallelisation->rank << ".pvd";
        
    ofstream os( namefile.c_str() );
    if (process.parallelisation->is_master_cpu())  os <<"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" \n\t byte_order=\"LittleEndian\" \n\t compressor=\"vtkZLibDataCompressor\" >\n <Collection> " <<std::endl;

    //     std::cout << S.size() << endl;
    if (process.parallelisation->is_multi_cpu()) 
       process.parallelisation->synchronisation();
    
    Sc2String local_filename="./inter/result";
    if(process.multiresolution->nb_calculs>1)
        local_filename<<"resolution_"<<process.multiresolution->m;
   
    unsigned data=process.affichage->side;
    
    for(unsigned i=0;i<Inter.size();++i) {
            if (SS[Inter[i].vois[data*2]].num_proc==process.parallelisation->rank){	
	      unsigned side = (Inter[i].type == Interface::type_ext)? 0 : data;
	      for(unsigned k=1;k<Inter[i].side[side].t_post.size();k++){
		ostringstream ss, ss2;
		ss<<local_filename << "_inter_id_"<<Inter[i].id<<"_time_"<<k<<".vtu";
		Sc2String namefile_entity(ss.str());
		ss2 << " <DataSet timestep=\" "<< k<< " \" group=\" " <<Inter[i].id << "\" part=\" \"\n\t\t file=\"" << namefile_entity << "\"/>\n" ;
		Sc2String ligne1(ss2.str()); 
		os << ligne1;
	      }
            }
    }    
    
    os.close();
    ///Fin des écritures
    if (process.parallelisation->is_multi_cpu()) {
      process.parallelisation->synchronisation();
           
      if (process.parallelisation->is_master_cpu()){
	///concaténation des fichiers
	Sc2String cmd;
	cmd << "cat ";
	for (unsigned i=0;i<process.parallelisation->size;i++)
	  cmd << base_filename << i << ".pvd ";
	cmd << "> " << save_directory;
	Sc2String name_file_output="results_inter";
	if(process.multiresolution->nb_calculs>1)
	  name_file_output<<"_resolution_"<<process.multiresolution->m;
	name_file_output<<".pvd";
	cmd << name_file_output ;
	int tmp=system(cmd.c_str());

	Sc2String namefile2 = process.affichage->repertoire_save+"results/"+name_file_output;
	fstream os2;
	os2.open ( namefile2.c_str(), fstream::in | fstream::out | fstream::ate );
	os2 << " </Collection> \n </VTKFile>";
	os2.close();
	Sc2String cmd2;
	cmd2 << "rm ";
	for (unsigned i=0;i<process.parallelisation->size;i++)
	  cmd2 << base_filename << i << ".pvd ";
	tmp=system(cmd2.c_str());
	
      }
      process.parallelisation->synchronisation();
    } else
    {
        
	Sc2String cmd;
	cmd << "mv " << base_filename << "0" << ".pvd ";
	cmd << " " << save_directory;
	Sc2String name_file_output="results_inter";
	if(process.multiresolution->nb_calculs>1)
	  name_file_output<<"_resolution_"<<process.multiresolution->m;
	name_file_output<<".pvd";
	cmd << name_file_output ;
	int tmp=system(cmd.c_str());

	Sc2String namefile2 = process.affichage->repertoire_save+"results/"+name_file_output;
	fstream os2;
	os2.open ( namefile2.c_str(),fstream::in | fstream::out | fstream::ate  );
	os2 << " </Collection> \n </VTKFile>";
	os2.close();
     }    
}


/**\ingroup Post_Traitement
 \brief Affichage des champs après calcul sur les sous-structures
 
 On appelle affich_SST_resultat() pour créer le fichier paraview de résultat pour chaque pas de temps.
 */
void affichage_resultats(Vec<VecPointedValues<Sst> > &S,  Process &process, DataUser &data_user) {
    //PRINT(process.affichage->affich_resultat);
    
    if (process.affichage->affich_resultat==1)
        if (process.parallelisation->is_local_cpu()) {
            write_paraview_results(S,process, data_user);
            //create_file_pvd(process,data_user,"sst_bulk");    /// TMP, test sauvegarde a la fin de chaque pas de temps
            //create_file_pvd(process,data_user,"sst_skin");    /// TMP, test sauvegarde a la fin de chaque pas de temps
        }
};


/**\ingroup Post_Traitement
 \brief Affichage des champs après calcul sur les interfaces
 
 On appelle affich_resultats_inter() pour créer le fichier paraview de résultat pour chaque pas de temps.
 */
void affichage_resultats_inter(Vec<VecPointedValues<Interface> > &Inter, Vec<Sst> &S , Process &process, DataUser &data_user) {
  if (process.affichage->affich_resultat==1)
      if (process.parallelisation->is_local_cpu()) {
        affich_INTER_resultat(Inter,S,process);
        //create_file_pvd(process,data_user,"inter");   /// TMP, test sauvegarde a la fin de chaque pas de temps
      }
};

/** \ingroup Post_Traitement
\brief Affichage des champs après calcul sur les interfaces
 
Possibilité de choisir une interface donnée ou toutes les interfaces.
*/
template <class TV1,class TV2> 
void affichage_inter_data(TV2 &Inter, TV1 &S, Process &process){
    if (process.affichage->affich_inter_data==1) {
        affich_inter_data_time(Inter,S,process);
    }
}


/** \ingroup Post_Traitement
\brief Affichage de l'evolution du déplacement d'un point donné par ses coordonnées 
*/
template <class TV3> 
void affichage_depl_pt(TV3 &S, Process &process){
    if(process.affichage->affich_depl_pt==1) extraction_depl_pt(S, process);
}

/** \ingroup Post_Traitement
\brief Affichage de l'évolution de l'énergie dissipée ou de l'énergie imposée au cours du temps à partir des quantités chapeaux ou des quantités n de l'interface
 
 Selon les parametres du champ SavingData::param_ener, on sélectionne le type d'énergie et les quantités retenues.
0 - 0 : energie dissipee sur les quantites chapeau
0 - 1 : energie dissipee sur les quantites n
1 - 0 : energie imposee sur les quantites chapeau
1 - 1 : energie imposee sur les quantites n

Faux en MPI pour certaine fonction qui necessite d avoir les W des deux cotes et ils sont pas transferes pour le posttraitement : est-ce utile de le faire ? non pour le moment
*/
template <class TV3,class TV2> 
void affichage_energie(TV3 &S,TV2 &Inter, Process &process, DataUser &data_user){
    Vec<double> energie,temp;
    energie.resize(process.temps->nbpastemps + 1);temp.resize(process.temps->nbpastemps + 1);energie.set(0.);temp.set(0.);
    if(process.affichage->param_ener[0]==0 and process.affichage->param_ener[1]==0) {
        if (process.parallelisation->is_local_cpu()) calcul_ener_dissi_chap(S,Inter,energie,process);
        if (process.parallelisation->is_multi_cpu()) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.parallelisation->is_master_cpu()) std::cout << "Dissipation (c) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==0 and process.affichage->param_ener[1]==1) {
        if (process.parallelisation->is_local_cpu()) calcul_ener_dissi_lin(S,Inter,energie,process);
        if (process.parallelisation->is_multi_cpu()) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.parallelisation->is_master_cpu()) std::cout << "Dissipation (n) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==1 and process.affichage->param_ener[1]==0) {
        if (process.parallelisation->is_local_cpu()) calcul_ener_imp_chap(S,Inter,energie,process);
        if (process.parallelisation->is_multi_cpu()) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.parallelisation->is_master_cpu()) std::cout << "Energie imposee (c) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==1 and process.affichage->param_ener[1]==1) {
        if (process.parallelisation->is_local_cpu()) calcul_ener_imp_lin(S,Inter,energie,process);
        if (process.parallelisation->is_multi_cpu()) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.parallelisation->is_master_cpu()) std::cout << "Energie imposee (n) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==2 and process.affichage->param_ener[1]==0) {
        if (process.parallelisation->is_local_cpu()) calcul_Ft2_chap(S,Inter,energie,process);
        if (process.parallelisation->is_multi_cpu()) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.parallelisation->is_master_cpu()) std::cout << "Effort tangent carre (c) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==2 and process.affichage->param_ener[1]==1) {
        if (process.parallelisation->is_local_cpu()) calcul_Ft2_lin(S,Inter,energie,process);
        if (process.parallelisation->is_multi_cpu()) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.parallelisation->is_master_cpu()) std::cout << "Effort tangent carre (n) faux en mpi : " << energie << endl;
    } else if(process.affichage->param_ener[0]==3 and process.affichage->param_ener[1]==0) {
        if (process.parallelisation->is_local_cpu()) calcul_Fn_chap(S,Inter,energie,process);
        if (process.parallelisation->is_multi_cpu()) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.parallelisation->is_master_cpu()) std::cout << "Effort normal (c) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==3 and process.affichage->param_ener[1]==1) {
        if (process.parallelisation->is_local_cpu()) calcul_Fn_lin(S,Inter,energie,process);
        if (process.parallelisation->is_multi_cpu()) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.parallelisation->is_master_cpu()) std::cout << "Effort normal (n) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==4 and process.affichage->param_ener[1]==0) {
        if (process.parallelisation->is_local_cpu()) calcul_Un_chap(S,Inter,energie,process);
        if (process.parallelisation->is_multi_cpu()) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.parallelisation->is_master_cpu()) std::cout << "Deplacement normal (c) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==4 and process.affichage->param_ener[1]==1) {
        if (process.parallelisation->is_local_cpu()) calcul_Un_lin(S,Inter,energie,process);
        if (process.parallelisation->is_multi_cpu()) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.parallelisation->is_master_cpu()) std::cout << "Deplacement normal (n) faux en mpi  : " << energie << endl;
    } else if(process.affichage->param_ener[0]==5 and process.affichage->param_ener[1]==0) {
        if (process.parallelisation->is_local_cpu()) calcul_Ut_chap(S,Inter,energie,process);
        if (process.parallelisation->is_multi_cpu()) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.parallelisation->is_master_cpu()) std::cout << "Deplacement tangent (c) : " << energie << endl;
    } else if(process.affichage->param_ener[0]==5 and process.affichage->param_ener[1]==1) {
        if (process.parallelisation->is_local_cpu()) calcul_Ut_lin(S,Inter,energie,process);
        if (process.parallelisation->is_multi_cpu()) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.parallelisation->is_master_cpu()) std::cout << "Deplacement tanget (n) faux en mpi  : " << energie << endl;
    }  else if(process.affichage->param_ener[0]==6) {
        if (process.parallelisation->is_local_cpu()) calcul_energie_elastique(S,Inter,energie,process, data_user);
        if (process.parallelisation->is_multi_cpu()) {MPI_Reduce(energie.ptr(),temp.ptr(),temp.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);energie=temp;}
        if (process.parallelisation->is_master_cpu()) std::cout << "Energie elastique : " << energie << endl;
    } else {
        std::cout << "Mauvais choix d'energie" << endl;
        assert(0);
    }
    if (process.affichage->command_file=="" and process.parallelisation->is_master_cpu()){
        GnuPlot gp;
        gp.plot(energie);
        gp.wait();
    }
}


