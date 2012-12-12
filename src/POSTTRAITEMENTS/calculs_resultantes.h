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
template <class TV1,class TV2>
void calcul_resultante(TV1 &S, TV2 &Inter,Process &process) {
    Vector resultante;
    resultante.set(0.);
    
     ///preparation des noms et des repertoires pour ecriture des resultats
    Sc2String save_directory=process.affichage->repertoire_save+"results/";
    Sc2String base_filename= save_directory;
    if(process.multiresolution->nb_calculs>1)
        base_filename<<"resolution_"<<process.multiresolution->m<<"_";
    base_filename << "proc_" << process.parallelisation->rank << "_time_";
    Sc2String namefile = base_filename;
    namefile << process.temps->pt_cur << ".csv";
    
    
    for(unsigned i=0;i<Inter.size();i++) {
        ///résultantes
        Scalar resx,resy,resz;
        resx = 0;
        resy = 0;
        resz = 0;
                
        if (Inter[i].comp == Interface::comp_cassable_parfait or
            Inter[i].comp == Interface::comp_cassable_elastique or
            Inter[i].comp == Interface::comp_contact_parfait or
            Inter[i].comp == Interface::comp_contact_elastique or
            Inter[i].comp == Interface::comp_cohesive) {
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
                resx = dot(vecx,Inter[i].side[0].t_post[1].F);
                resy = dot(vecy,Inter[i].side[0].t_post[1].F);
                if(DIM == 3){
                    resz = dot(vecz,Inter[i].side[0].t_post[1].F);
                }  
            } else {
                std::cout << "Nom de calcul nom pris en compte" << endl;
                assert(0);
            }
            std::cout << "Interface " << Inter[i].num << ", entre les pieces : " << Inter[i].vois[0] << " " << Inter[i].vois[2] << ", resx " << resx << ", resy " << resy << ", resz " << resz <<  endl;
        }
    }
}
#endif //CALCULS_RESULTANTES_H



