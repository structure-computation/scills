#include "Sst.h"



Sst::Sst() : pb(*mesh.m,true),plastique(false) {}     ///< constructeur de la formulation pour la sous-structure

Sst::~Sst() {free();}    ///destructeur de la SST

void Sst::free(){
    vois.free();
    edge.free();
    LE.free();
    ///K est libéré apres factorisation de toute facon et non calculer où y a pas besoin
    fvol.free();
    t.free();
    t_post.free();
}

/// Trouver une sst à partir de son id----------------------------------------------
Sst* Sst::find_sst(Vec<Sst> &S,int id_) {
    for (int i_group=0; i_group<S.size(); i_group++) {
        if (S[i_group].id == id_) {
            return &S[i_group];
            break;
        }
    }
}

/// Trouver l'index d'une sst à partir de son id -----------------------------------
int Sst::find_index_sst(Vec<Sst> &S, int id_) {
    for (int i_group=0; i_group<S.size(); i_group++) {
        if (S[i_group].id == id_) {
            return i_group;
            break;
        }
    }
}

void Sst::assign_material_on_element(DataUser &data_user){
    if (mesh.type_formulation=="isotrope") {
        ///formulation isotrope 
        mesh->elastic_modulus = matprop.elastic_modulus ; 
        mesh->poisson_ratio   = matprop.poisson_ratio   ; 
        mesh->deltaT          = matprop.deltaT          ; 
        mesh->resolution      = matprop.resolution      ; 
        mesh->alpha           = matprop.alpha           ; 
        mesh->f_vol           = matprop.f_vol           ; 
        mesh->density         = matprop.density         ; 
        mesh.load_f_vol_e(matprop.f_vol_e,data_user);
    } else if (mesh.type_formulation=="orthotrope") {
        ///formulation orthotrope 
        mesh->elastic_modulus_1 = matprop.elastic_modulus_1;
        mesh->elastic_modulus_2 = matprop.elastic_modulus_2;
        mesh->elastic_modulus_3 = matprop.elastic_modulus_3;
        mesh->poisson_ratio_12  = matprop.poisson_ratio_12 ;
        mesh->poisson_ratio_13  = matprop.poisson_ratio_13 ;
        mesh->poisson_ratio_23  = matprop.poisson_ratio_23 ;
        mesh->shear_modulus_12  = matprop.shear_modulus_12 ;
        mesh->shear_modulus_13  = matprop.shear_modulus_13 ;
        mesh->shear_modulus_23  = matprop.shear_modulus_23 ;
        mesh->v1                = matprop.v1               ;
        mesh->v2                = matprop.v2               ;
        mesh->deltaT            = matprop.deltaT           ;
        mesh->resolution        = matprop.resolution       ;
        mesh->alpha_1           = matprop.alpha_1          ;
        mesh->alpha_2           = matprop.alpha_2          ;
        mesh->alpha_3           = matprop.alpha_3          ;
        mesh->f_vol             = matprop.f_vol            ;
        mesh->density           = matprop.density          ;
        
        if(matprop.comp == "mesomodele"){
            mesh->k_p      = matprop.k_p;
            mesh->m_p      = matprop.m_p;
            mesh->R0       = matprop.R0;
            mesh->couplage = matprop.coefvm_composite;
        }
        
        if(matprop.comp == "mesomodele"){
            mesh->Yo           = matprop.Yo;
            mesh->Yc           = matprop.Yc;
            mesh->Ycf          = matprop.Ycf;
            mesh->dmax         = matprop.dmax;
            mesh->b_c          = matprop.b_c;
            mesh->effet_retard = matprop.effet_retard;
            mesh->a            = matprop.a;
            mesh->tau_c        = matprop.tau_c;
        }
        
        mesh.load_f_vol_e(matprop.f_vol_e,data_user);
    }
}


void Sst::Time::allocations(unsigned nbddl,unsigned nbelem,const Process &process,bool plastique){
    if (process.rank>0 or process.size==1){
        /// Deplacements
        q.resize(nbddl);
        q.set(0.0);
        /// Comportement plastique
        if(plastique){
            /// Plasticite cumulee
            p.resize(nbelem);
            p.set(0.0);
            /// Ecrouissage isotrope
            R.resize(nbelem);
            R.set(0.0);
            /// Deformations plastiques
            epsilon_p.resize(nbelem);
            for( unsigned i = 0; i < nbelem; ++i ) {
                epsilon_p[i].resize( DIM*(DIM+1)/2);
                epsilon_p[i].set(0.0);
            }
            /// Contraintes
            sigma.resize(nbelem);
            for( unsigned i = 0; i < nbelem; ++i ) {
                sigma[i].resize( DIM*(DIM+1)/2);
                sigma[i].set(0.0);
            }
        }
        /* Debug
        if(plastique){
            std::cout << "--------------TEST-------------" << std::endl;
            std::cout << this << "   " << epsilon_p.size() << std::endl;
            for(int i = 0; i < nbelem; i++){
                std::cout << i << ": ";
                for(int j = 0; j < DIM*(DIM+1)/2; j++)
                    std::cout << epsilon_p[0][i] << " , ";
                std::cout << endl;
            }
            std::cout << endl;
        }
        //*/
    }
}

void Sst::Time::display_all_data(){
    std::cout << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "********************** Debug Sst::Time : ********************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "q         : " << q.size()         << std::endl;
    std::cout << "p         : " << p.size()         << std::endl;
    std::cout << "R         : " << R.size()         << std::endl;
    std::cout << "epsilon_p : " << epsilon_p.size() << std::endl;
    std::cout << "sigma     : " << sigma.size()     << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << std::endl;
}
