#ifndef CALCULS_RESULTANTES_H
#define CALCULS_RESULTANTES_H
//
// C++ Implementation: calculs_energies
//
// Description:
//
//
// Author: Jeemie Bellec <bellec@structure-computation.com>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//

//librairie LMTpp
#include "../../LMT/include/containers/mat.h"
#include "../../LMT/include/containers/vec.h"
#include "../../LMT/include/mesh/mesh.h"

#include <fstream>
#include <map>

// fichiers de definition des variables
#include "../DEFINITIONS/Process.h"
#include "../DEFINITIONS/SavingData.h"
#include "../DEFINITIONS/LatinData.h"
#include "../DEFINITIONS/Sst.h"
#include "../DEFINITIONS/Interface.h"


/**
Fonction permettant le calcul des resultantes sur une interface de contact.
*/
template <class TV1,class TV2, class TV3>
void calcul_resultante(TV1 &S, TV3 &SS, TV2 &Inter,Process &process) {
    Vector resultante;
    resultante.set(0.);
    
     ///preparation des noms et des repertoires pour ecriture des resultats
    Sc2String save_directory=process.affichage->repertoire_save+"results/";
    Sc2String base_filename= save_directory;
    if(process.multiresolution->nb_calculs>1)
        base_filename<<"resolution_"<<process.multiresolution->m<<"_";
    base_filename << "proc_" ; 
    Sc2String namefile = base_filename;
    namefile << process.parallelisation->rank << ".csv";
    
    ofstream os( namefile.c_str() );
    if (process.parallelisation->is_master_cpu())  os << "Numero interface ; Numero resolution ; Numero pas de temps ; Comportement interface ; Numero SST cote ; Numero SST cote 2 ; Materiau 1 ; Materiau 2 ; Fx ; Fy ; Fz ; Ux ; Uy ; Uz;" << endl;
//     std::cout << S.size() << endl;
    if (process.parallelisation->is_multi_cpu()) 
       process.parallelisation->synchronisation();
    
    for(unsigned j=0;j<S.size();j++) {
        for(unsigned e=0;e<S[j].edge.size();++e) {
            unsigned i=S[j].edge[e].internum;
            unsigned data=S[j].edge[e].datanum;
	    
// 	    std::cout << "interface : " << i << " cote " << data << " comp : " << Inter[i].comp << endl;

	    ///résultantes
	    Scalar resx,resy,resz;
		    
// 	    if ((Inter[i].comp == Interface::comp_cassable_parfait or
// 		Inter[i].comp == Interface::comp_cassable_elastique or
// 		Inter[i].comp == Interface::comp_contact_parfait or
//  		Inter[i].comp == Interface::comp_contact_elastique or
// 		Inter[i].comp == Interface::comp_cohesive) and data==0) {
	    if (data==0) {
// 	        std::cout << "type calcul " << process.nom_calcul << endl;
		if(process.nom_calcul=="incr") {
		    Vector vecx,vecy,vecz;
		    ///vecteur de projections
		    vecx.resize(Inter[i].side[0].nodeeq.size()*DIM);
		    vecx.set(0.);
		    vecy.resize(Inter[i].side[0].nodeeq.size()*DIM);
		    vecy.set(0.);
		    vecz.resize(Inter[i].side[0].nodeeq.size()*DIM);
		    vecz.set(0.);
		    for( unsigned ii=0;ii<Inter[i].side[0].nodeeq.size();ii++ ) {
			vecx[DIM*ii]=1.;
			vecy[DIM*ii+1]=1.;
			if(DIM == 3){
			    vecz[DIM*ii+2]=1.;
			}
		    }
		    
		    ///NumSst
		    int nums1 = Inter[i].side[0].vois[0];
		    int nums2 = -1;
		    if (Inter[i].side.size() > 1)
		      nums2 = Inter[i].side[1].vois[0];
		    
		    ///Pour chaque pas de temps
		    for(unsigned k=0;k<Inter[i].side[0].t_post.size();k++){
			resx = 0;
			resy = 0;
			resz = 0;
		        Vector vectemp;
			vectemp.resize(vecx.size());
			///Résultante effort
			vectemp=Inter[i].side[0].M*Inter[i].side[0].t_post[k].Fchap;
			resx = dot(vecx,vectemp);
			resy = dot(vecy,vectemp);
			if(DIM == 3){
			    resz = dot(vecz,vectemp);
			}
			///Ecriture dans le fichiers
			os << i << ";" << ";" << k << ";" << Inter[i].comp << ";" << nums1 << ";" << nums2 << ";" << SS[nums1].typmat << ";" ;
			if (nums2 == -1) 
			  os <<  ";" ;
			else
			  os << SS[nums2].typmat << ";";
			os << resx << ";" << resy << ";" << resz << ";" ;

			///Résultante déplacement
			vectemp=Inter[i].side[0].M*Inter[i].side[0].t_post[k].W;
// 			std::cout << Inter[i].side[0].t_post[k].Wchap << endl;
			resx = dot(vecx,vectemp)/Inter[i].measure;
			resy = dot(vecy,vectemp)/Inter[i].measure;
			if(DIM == 3){
			    resz = dot(vecz,vectemp)/Inter[i].measure;
			}
			os << resx << ";" << resy << ";" << resz << ";" << endl;
		    }
		} else {
		    std::cout << "Nom de calcul nom pris en compte" << endl;
		    assert(0);
		}
// 		std::cout << "Interface " << Inter[i].num << ", entre les pieces : " << Inter[i].vois[0] << " " << Inter[i].vois[2] << ", resx " << resx << ", resy " << resy << ", resz " << resz <<  endl;
	    }
	}
    }
    
    ///Fin des écritures
    if (process.parallelisation->is_multi_cpu()) {
      process.parallelisation->synchronisation();
      
      if (process.parallelisation->is_master_cpu()){
	///concaténation des fichiers
	Sc2String cmd;
	cmd << "cat ";
	for (unsigned i=0;i<process.parallelisation->size;i++)
	  cmd << base_filename << i << ".csv ";
	cmd << "> " << save_directory;
	if(process.multiresolution->nb_calculs>1)
	  cmd<<"resolution_"<<process.multiresolution->m<<"_";
	cmd << "resultante.csv" ; 
	system(cmd.c_str());
// 	std::cout << cmd << endl;
	
	Sc2String cmd2;
	cmd2 << "rm ";
	for (unsigned i=0;i<process.parallelisation->size;i++)
	  cmd2 << base_filename << i << ".csv ";
	system(cmd2.c_str());
// 	std::cout << cmd2 << endl;
      }
    } else
    {
	Sc2String cmd;
	cmd << "mv " << base_filename << "0" << ".csv ";
	cmd << " " << save_directory;
	if(process.multiresolution->nb_calculs>1)
	  cmd<<"resolution_"<<process.multiresolution->m<<"_";
	cmd << "resultante.csv" ; 
	system(cmd.c_str());
      
      
    }
    
    os.close();
}
#endif //CALCULS_RESULTANTES_H



