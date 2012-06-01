#include "Sst.h"



Sst::Sst() : pb(*mesh.m,true),plastique(false),endommageable(false),maj_operateurs(false) {}     /// constructeur de la formulation pour la sous-structure

Sst::~Sst() {free();}    ///destructeur de la SST

void Sst::free(){
    vois.free();
    edge.free();
    LE.free();
    ///K est liberee apres factorisation de toute facon et non calculer où y a pas besoin
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

void Sst::assign_material_on_element(const DataUser &data_user){
    if (matprop.type.find("isotrope")<matprop.type.size()) {
        ///formulation isotrope 
        mesh->elastic_modulus = matprop.elastic_modulus ; 
        mesh->poisson_ratio   = matprop.poisson_ratio   ; 
        mesh->deltaT          = matprop.deltaT          ; 
        mesh->resolution      = matprop.resolution      ; 
        mesh->alpha           = matprop.alpha           ; 
        mesh->f_vol           = matprop.f_vol           ; 
        mesh->density         = matprop.density         ; 
        mesh.load_f_vol_e(matprop.f_vol_e,data_user);
    } else if (matprop.type.find("orthotrope")<matprop.type.size()) {
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
        mesh.load_f_vol_e(matprop.f_vol_e,data_user);
    }
        
    if(matprop.comp.find("plastique")<matprop.comp.size()){
        mesh->plast_ecrouissage_init = matprop.plast_ecrouissage_init;
        mesh->plast_ecrouissage_mult = matprop.plast_ecrouissage_mult;
        mesh->plast_ecrouissage_expo = matprop.plast_ecrouissage_expo;
        mesh->plast_cinematique_coef = matprop.plast_cinematique_coef;
    }
    
    if(matprop.comp.find("endommageable")<matprop.comp.size()){
        mesh->Yo           = matprop.Yo;
        mesh->b_c          = matprop.b_c;
    }
    
    if(matprop.comp.find("mesomodele")<matprop.comp.size()){
        mesh->Yo           = matprop.Yo;
        mesh->Yc           = matprop.Yc;
        mesh->Ycf          = matprop.Ycf;
        mesh->dmax         = matprop.dmax;
        mesh->b_c          = matprop.b_c;
        mesh->effet_retard = matprop.effet_retard;
        mesh->a            = matprop.a;
        mesh->tau_c        = matprop.tau_c;
    }
}


void Sst::Time::allocations(unsigned nbddl,unsigned nbelem,const Process &process,Sst &S){
    if (process.parallelisation->is_local_cpu()){
        /// Deplacements
        q.resize(nbddl);
        q.set(0.0);
        /// Comportement plastique
        if(S.f == S.pb.formulation_plasticity_isotropy_stat_Qstat){
            /// Plasticite cumulee
            p.resize(nbelem);
            p.set(0.0);
            /// Ecrouissage isotrope
            R_p.resize(nbelem);
            R_p.set(0.0);
            /// Deformations plastiques
            epsilon_p.resize(nbelem);
            for( unsigned i = 0; i < nbelem; ++i ) {
                epsilon_p[i].set(0.0);
            }
            /// Origines des domaines elastiques
            if(S.matprop.type_plast.find("cinematique") < S.matprop.type_plast.size()){
                X_p.resize(nbelem);
                for( unsigned i = 0; i < nbelem; ++i ) {
                    X_p[i].set(0.0);
                }
            }
        }
        /// Comportement endommageable
        if(S.f == S.pb.formulation_elasticity_damageable_isotropy_stat_Qstat or S.f == S.pb.formulation_mesomodele){
            /// Endommagement standart ou micro-fissuration du mesomodele
            d1.resize(nbelem);
            d1.set(0.0);
        }
        /// Comportement mesomodele
        if(S.f == S.pb.formulation_mesomodele){
            /// Decohesion fibres/matrice du mesomodele
            d2.resize(nbelem);
            d2.set(0.0);
            /// Rupture fragile des fibres du mesomodele
            df.resize(nbelem);
            df.set(0.0);
            /// Forces associees
            Yd.resize(nbelem);
            for( unsigned i = 0; i < nbelem; ++i )
                Yd[i].set(0.0);
        }
    }
}

void Sst::Time::display_all_data(){
    std::cout << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "********************** Debug Sst::Time : ********************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "q         : " << q.size()         << std::endl;
    std::cout << "p         : " << p.size()         << std::endl;
    std::cout << "R_p       : " << R_p.size()       << std::endl;
    std::cout << "epsilon_p : " << epsilon_p.size() << std::endl;
    std::cout << "d1        : " << d1.size()        << std::endl;
    std::cout << "d2        : " << d2.size()        << std::endl;
    std::cout << "df        : " << df.size()        << std::endl;
    std::cout << "Yd        : " << Yd.size()        << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << std::endl;
}
