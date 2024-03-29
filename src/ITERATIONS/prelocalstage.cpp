#include "prelocalstage.h"
#include "../COMPUTE/DataUser.h"
#include "Process.h"
#include "../DEFINITIONS/TimeData.h"
#include "../DEFINITIONS/Sst.h"
#include "../DEFINITIONS/Interface.h"
#include "../DEFINITIONS/Boundary.h"
#include "../UTILITAIRES/utilitaires.h"
#include "../../LMT/include/codegen/codegen.h"
#include "../../LMT/include/containers/basicops.h"
#include "../../UTILITAIRES/algebre.h"

using namespace Codegen;



/// Calcul des valeurs numeriques des fonctions spatiales de la CL et assignation dans le vecteur V associe au chargement defini (Fchap ou Wpchap)
void assign_CL_spatial_temporel(Vector &V, Vec<Point > &nodeeq, Boundary &CL,Scalar dt) {
    Ex::MapExNum values = Boundary::CL_parameters.getParentsValues();               /// Recuperation des parametres temporels et de multiresolution
    for(unsigned i=0;i<nodeeq.size();++i){
        for(unsigned i_dir=0;i_dir<DIM;++i_dir){
            values[Boundary::CL_parameters.main_parameters[i_dir]->self_ex] = nodeeq[i][i_dir];  /// Chargement des coordonnees du point (main_parameters pointe vers x, y et z)
        }
        Point data;
        for(unsigned i_dir=0;i_dir<DIM;++i_dir){
            data[i_dir] = CL.fcts_spatiales[i_dir].updateValue(values);     /// Evaluation des composantes de la CL
        }
        V[range(i*DIM,(i+1)*DIM)]=data;                                     /// Stockage du resultat
    }
}


/// Calcul des valeurs numeriques des fonctions spatiales de la CL en deplacement normal et assignation dans le vecteur V associe (Wpchap)
void assign_CL_spatial_temporel_normale(Vector &V, Vec<Point > &nodeeq, Vector &neqs, Boundary &CL,Scalar dt) {
    Ex::MapExNum values = Boundary::CL_parameters.getParentsValues();   /// Recuperation des parametres temporels et de multiresolution
    for(unsigned i=0;i<nodeeq.size();++i){
        for(unsigned i_dir=0;i_dir<DIM;++i_dir){
            values[Boundary::CL_parameters.main_parameters[i_dir]->self_ex] = nodeeq[i][i_dir]; /// Chargement des coordonnees du point
        }
        //PRINT(CL.fcts_spatiales);
        Scalar data = CL.fcts_spatiales[0].updateValue(values); /// Evaluation de la composante normale
        //PRINT(data);
        Point temp=V[range(i*DIM,(i+1)*DIM)];                   /// Recuperation de la valeur actuelle sur l'interface
        Point neq = neqs[range(i*DIM,(i+1)*DIM)];               /// Recuperation de la normale de l'element
        V[range(i*DIM,(i+1)*DIM)]=ProjT(temp,neq)+data*neq;     /// Calcul et stockage du resultat
    }
    PRINT(V);
}


void initialise_CL_values(PointedInterfaces &Inter, Boundaries &CL){
    Ex::MapExNum values = Boundary::CL_parameters.getParentsValues();
    for(unsigned i_inter = 0; i_inter < Inter.size(); i_inter++){
        /// Teste si la CL doit etre initialisee
        if(Inter[i_inter].type != Interface::type_ext or (Inter[i_inter].comp != Interface::comp_deplacement and Inter[i_inter].comp != Interface::comp_deplacement_normal)){
            continue;   /// Si inutile, passer a l'interface suivante
        }
        bool depl = Inter[i_inter].comp == Interface::comp_deplacement;
        bool depl_normal = Inter[i_inter].comp == Interface::comp_deplacement_normal;
        for(unsigned i=0;i<Inter[i_inter].side[0].nodeeq.size();++i){
            for(unsigned i_dir = 0; i_dir < DIM; ++i_dir){
                values[Boundary::CL_parameters.main_parameters[i_dir]->self_ex] = Inter[i_inter].side[0].nodeeq[i][i_dir];
            }
//             if(depl){
//                 /// Cas d'un deplacement impose
//                 Point data;
//                 for(unsigned i_dir = 0; i_dir < DIM; i_dir++){
//                     data[i_dir] = (Scalar)CL[Inter[i_inter].refCL].fcts_spatiales[i_dir].updateValue(values);
//                 }
//                 Inter[i_inter].side[0].t[0].W[range(i*DIM,(i+1)*DIM)]=data;
//             }else if(depl_normal){
//                 /// Cas d'un deplacement normal
//                 Scalar data = (Scalar)CL[Inter[i_inter].refCL].fcts_spatiales[0].updateValue(values);
//                 Point neq = Inter[i_inter].side[0].neq[range(i*DIM,(i+1)*DIM)];
//                 Point temp=Inter[i_inter].side[0].t[0].W[range(i*DIM,(i+1)*DIM)];
//                 Inter[i_inter].side[0].t[0].W[range(i*DIM,(i+1)*DIM)]=ProjT(temp,neq)+data*neq;
//             }else{
//                 /// Rien a faire pour les autres type de CL
//             }
        }
    }
}


/// Mise a jour des CL et assignation des nouvelles valeurs a Fchap et Wpchap (pour une iteration incrementale)
void update_CL_values(PointedInterfaces &Inter, Boundaries &CL, Process &process, DataUser &data_user ) {
    Scalar dt = process.temps->dt;
    std::cout << "Mise a jour des Conditions aux Limites sur le processeur " << process.parallelisation->rank << " : " << std::endl;
    for(unsigned i_inter=0;i_inter<Inter.size();++i_inter) {
        std::cout << "\tid : " << Inter[i_inter].id  << "\ttype : " << Inter[i_inter].type << "\tcomportement : " << Inter[i_inter].comp << "\tvaleur : " << Inter[i_inter].comp;
        if (Inter[i_inter].type == Interface::type_ext and Inter[i_inter].comp != Interface::comp_periodique) {
            if(Inter[i_inter].comp == Interface::comp_effort) {
                assign_CL_spatial_temporel(Inter[i_inter].side[0].t[1].Fchap,Inter[i_inter].side[0].nodeeq,CL[Inter[i_inter].refCL],dt);
                //std::cout << Inter[i_inter].side[0].t[1].Fchap[0] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Fchap[1] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Fchap[2] << " ; ";
            } else if(Inter[i_inter].comp == Interface::comp_effort_normal) {
                assign_CL_spatial_temporel_normale(Inter[i_inter].side[0].t[1].Fchap,Inter[i_inter].side[0].nodeeq,Inter[i_inter].side[0].neq,CL[Inter[i_inter].refCL],dt);
                //std::cout << Inter[i_inter].side[0].t[1].Fchap[0] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Fchap[1] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Fchap[2] << " ; ";
            } else if(Inter[i_inter].comp == Interface::comp_deplacement or Inter[i_inter].comp == Interface::comp_deplacement_nul) {
                assign_CL_spatial_temporel(Inter[i_inter].side[0].t[1].Wchap,Inter[i_inter].side[0].nodeeq,CL[Inter[i_inter].refCL],dt);
                Inter[i_inter].side[0].t[1].Wpchap = (Inter[i_inter].side[0].t[1].Wchap - Inter[i_inter].side[0].t[0].Wchap)/dt;
                PRINT(Inter[i_inter].side[0].t[1].Wpchap);
                PRINT(Inter[i_inter].side[0].t[0].W);
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[0] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[1] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[2] << " ; ";
            } else if(Inter[i_inter].comp == Interface::comp_vitesse or Inter[i_inter].comp == Interface::comp_vitesse_nulle) {
                assign_CL_spatial_temporel(Inter[i_inter].side[0].t[1].Wpchap,Inter[i_inter].side[0].nodeeq,CL[Inter[i_inter].refCL],dt);
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[0] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[1] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[2] << " ; ";
            } else if(Inter[i_inter].comp == Interface::comp_cinetic_torseur) {
                Inter[i_inter].assign_W_CL_torseur_cinetic_temporel(Inter[i_inter].side[0].t[1].Wchap,Inter[i_inter].side[0].nodeeq,CL[Inter[i_inter].refCL],true);
                //std::cout << Inter[i_inter].side[0].t[1].Fchap[0] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Fchap[1] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Fchap[2] << " ; ";
            } else if(Inter[i_inter].comp == Interface::comp_symetrie) {
                if(process.temps->pt_cur==1) {
                    if(process.reprise_calcul==0){
                        assign_CL_spatial_temporel(Inter[i_inter].side[0].t[1].Wpchap,Inter[i_inter].side[0].nodeeq,CL[Inter[i_inter].refCL],dt);
                    }
                    if(process.reprise_calcul==0){
                        Inter[i_inter].side[0].t[1].Fchap.set(0.0);
                    }
                }
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[0] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[1] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[2] << " ; ";
            } else if(Inter[i_inter].comp == Interface::comp_vitesse_normale) {
                if(process.temps->pt_cur==1) {
                    assign_CL_spatial_temporel_normale(Inter[i_inter].side[0].t[1].Wpchap,Inter[i_inter].side[0].nodeeq,Inter[i_inter].side[0].neq,CL[Inter[i_inter].refCL],dt);//le Wpchap evolue au cours des iterations donc si on reprend on initialise avec le resultat du calcul precedent donc on fait rien...
                    if(process.reprise_calcul==0)
                        Inter[i_inter].side[0].t[1].Fchap.set(0.0);
                } else {
                    Vector Wpchapnormal;
                    Wpchapnormal.resize(Inter[i_inter].side[0].t[1].Wpchap.size(),0.0);
                    assign_CL_spatial_temporel_normale(Wpchapnormal,Inter[i_inter].side[0].nodeeq,Inter[i_inter].side[0].neq,CL[Inter[i_inter].refCL],dt);
                    Inter[i_inter].side[0].t[1].Wpchap = Inter[i_inter].side[0].Pt(Inter[i_inter].side[0].t[1].Wpchap)+Wpchapnormal;
                }
                //PRINT(Inter[i_inter].side[0].t[1].Wpchap);
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[0] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[1] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[2] << " ; ";
            } else if(Inter[i_inter].comp == Interface::comp_deplacement_normal) {
                if(process.temps->pt_cur==1) {
                    assign_CL_spatial_temporel_normale(Inter[i_inter].side[0].t[1].Wchap,Inter[i_inter].side[0].nodeeq,Inter[i_inter].side[0].neq,CL[Inter[i_inter].refCL],dt);//le Wpchap evolue au cours des iterations donc si on reprend on initialise avec le resultat du calcul precedent donc on fait rien...
                    Inter[i_inter].side[0].t[1].Wpchap = (Inter[i_inter].side[0].t[1].Wchap - Inter[i_inter].side[0].t[0].Wchap)/dt;
                    if(process.reprise_calcul==0)
                        Inter[i_inter].side[0].t[1].Fchap.set(0.0);
                } else {
                    Vector Wchapnormal;
                    Wchapnormal.resize(Inter[i_inter].side[0].t[1].Wpchap.size(),0.0);
                    assign_CL_spatial_temporel_normale(Wchapnormal,Inter[i_inter].side[0].nodeeq,Inter[i_inter].side[0].neq,CL[Inter[i_inter].refCL],dt);
                    Inter[i_inter].side[0].t[1].Wpchap = Inter[i_inter].side[0].Pt(Inter[i_inter].side[0].t[1].Wpchap) + (Wchapnormal - Inter[i_inter].side[0].Pn(Inter[i_inter].side[0].t[0].Wchap))/dt;
                }
                //PRINT(Inter[i_inter].side[0].t[1].Wpchap);
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[0] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[1] << " ; ";
                //std::cout << Inter[i_inter].side[0].t[1].Wpchap[2] << " ; ";
            } else {
                std::cout << "Erreur d'interface ext - prelocalstage " << std::endl;
                assert(0);
            }
        } else if( Inter[i_inter].comp==Interface::comp_contact_parfait or Inter[i_inter].comp==Interface::comp_parfait or Inter[i_inter].comp==Interface::comp_elastique or Inter[i_inter].comp==Interface::comp_cassable_elastique or Inter[i_inter].comp==Interface::comp_cassable_parfait or Inter[i_inter].comp==Interface::comp_cohesive) {
          //else if(Inter[i_inter].comp=="Contact_jeu" or Inter[i_inter].comp=="Contact_jeu_physique" or Inter[i_inter].comp==Interface::comp_contact_ep or Inter[i_inter].comp==Interface::comp_parfait) {
            ///le jeu est reparti en moyenne sur chacun des deplacements des cotes 1 et 2
            //if(process.temps->pt_cur==1) {
                Vector dep_jeu = Inter[i_inter].Ep_imposee - Inter[i_inter].old_Ep_imposee ;
                Vector dep_precharge = ( Inter[i_inter].precharge - Inter[i_inter].old_precharge );
                
                Scalar R0 = Inter[i_inter].side[1].kn/(Inter[i_inter].side[1].kn+Inter[i_inter].side[0].kn);
                Scalar R1 = Inter[i_inter].side[0].kn/(Inter[i_inter].side[1].kn+Inter[i_inter].side[0].kn);
                
                Inter[i_inter].side[1].t[process.temps->pt-1].W[Inter[i_inter].side[1].ddlcorresp] += R1 * dep_jeu - dt * Inter[i_inter].side[1].hglo * dep_precharge ;
                Inter[i_inter].side[0].t[process.temps->pt-1].W += - 1. * R0 * dep_jeu + dt * Inter[i_inter].side[0].hglo * dep_precharge ;
                
//                 Inter[i_inter].side[1].t[process.temps->pt-1].F[Inter[i_inter].side[1].ddlcorresp] += - 1. * dep_precharge;
//                 Inter[i_inter].side[0].t[process.temps->pt-1].F +=  dep_precharge;
 
//                 if(Inter[i_inter].id==5){
//                     PRINT("  ");
//                     PRINT(R0);
//                     PRINT(R1);
//                     PRINT(dep_jeu[LMT::range(0,DIM*1)]);
//                     PRINT("on est dans prélocal stage");
//                     PRINT(Inter[i_inter].id);
//                     //PRINT(Inter[i_inter].side[0].t[process.temps->pt].dEp_imposee[LMT::range(0,DIM*1)]);
// //                     Vector testt1 ;
// //                     testt1.resize(dep_precharge.size(),1);
// //                     Vector testt = Inter[i_inter].side[0].kglo * ( (dt * Inter[i_inter].side[0].hglo * testt1)/process.temps->dt ) ;
// //                     PRINT(testt);
// //                     
// //                     
// //                     Scalar test;
// //                     test = 0;
// //                     for(unsigned i=0;i<Inter[i_inter].side[0].nodeeq.size()*DIM;++i) {
// //                         for(unsigned j=0;j<Inter[i_inter].side[0].nodeeq.size()*DIM;++j) {
// //                             test += Inter[i_inter].side[0].kglo(i,j) * Inter[i_inter].side[0].hglo(i,j);
// //                         }
// //                     }
// //                     PRINT(test);
// //                     PRINT((1./test));
//                     PRINT(Inter[i_inter].measure);
//                     PRINT(dep_precharge[LMT::range(0,DIM*1)]);
//                     PRINT(dt);
//                     PRINT("  ");
//                 }
                //if (Inter[i_inter].num == 15 ) std::cout << "Jeu (cote 0) : " << Inter[i_inter].side[0].t[0].Wchap << endl;
                //if (Inter[i_inter].num == 15 ) std::cout << "Jeu (cote 1) : " << Inter[i_inter].side[1].t[0].Wchap << endl;
            //}
        }
     
        std::cout << std::endl;
    }
    std::cout << "Fin de la mise a jour des Conditions aux Limites sur le processeur " << process.parallelisation->rank << " : " << std::endl;
}
