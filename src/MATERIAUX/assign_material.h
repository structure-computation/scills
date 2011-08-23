#ifndef ASSIGN_MATERIALS_SST_H
#define ASSIGN_MATERIALS_SST_H

using namespace LMT;
using namespace std;

//Attention fichier generer dans : generation_auto.py
template<class SST, class SSTCARAC> void assign_material(SST &S, Vec<SSTCARAC> &matprop, Param &process, DataUser &data_user){
//formulation isotrope 
if (matprop[S.typmat].type=="isotrope") {
	 S.f = S.pb.formulation_elasticity_isotropy_stat_Qstat; 
	 S.mesh.elastic_modulus=matprop[S.typmat].coef[0]; 
	 S.mesh.poisson_ratio=matprop[S.typmat].coef[1]; 
	 S.mesh.deltaT=matprop[S.typmat].coefth[1];  
	 S.mesh.resolution=matprop[S.typmat].resolution; 
	 S.mesh.alpha=matprop[S.typmat].coefth[0];  
	 S.mesh.f_vol=matprop[S.typmat].f_vol;  
	 S.mesh.f_vol_e=matprop[S.typmat].f_vol_e;  
	 S.mesh.density=matprop[S.typmat].density; 
	 S.mesh.type_formulation="isotrope"; 
}
//formulation orthotrope 
if (matprop[S.typmat].type=="orthotrope") {
	 S.f = S.pb.formulation_elasticity_orthotropy_stat_Qstat; 
	 S.mesh.elastic_modulus_1=matprop[S.typmat].coef[0];  
	 S.mesh.elastic_modulus_2=matprop[S.typmat].coef[1];  
	 S.mesh.elastic_modulus_3=matprop[S.typmat].coef[2];  
	 S.mesh.poisson_ratio_12=matprop[S.typmat].coef[3];   
	 S.mesh.poisson_ratio_13=matprop[S.typmat].coef[4];   
	 S.mesh.poisson_ratio_23=matprop[S.typmat].coef[5];   
	 S.mesh.shear_modulus_12=matprop[S.typmat].coef[6];   
	 S.mesh.shear_modulus_13=matprop[S.typmat].coef[7];   
	 S.mesh.shear_modulus_23=matprop[S.typmat].coef[8];   
	 S.mesh.v1=matprop[S.typmat].direction[0];            
	 S.mesh.v2=matprop[S.typmat].direction[1];            
	 S.mesh.deltaT=matprop[S.typmat].coefth[3];           
	 S.mesh.resolution=matprop[S.typmat].resolution;      
	 S.mesh.alpha_1=matprop[S.typmat].coefth[0];          
	 S.mesh.alpha_2=matprop[S.typmat].coefth[1];          
	 S.mesh.alpha_3=matprop[S.typmat].coefth[2];          
	 S.mesh.f_vol=matprop[S.typmat].f_vol;  
	 S.mesh.density=matprop[S.typmat].density; 
	 S.mesh.f_vol_e=matprop[S.typmat].f_vol_e;  
	 S.mesh.type_formulation="orthotrope"; 
}
}

template<class T1> void assign_material_on_element(T1 &S, DataUser &data_user){
//formulation isotrope 
if (S.mesh.type_formulation=="isotrope") {
	 S.mesh->elastic_modulus= S.mesh.elastic_modulus ; 
	 S.mesh->poisson_ratio  = S.mesh.poisson_ratio   ; 
	 S.mesh->deltaT         = S.mesh.deltaT          ; 
	 S.mesh->resolution     = S.mesh.resolution      ; 
	 S.mesh->alpha          = S.mesh.alpha           ; 
	 S.mesh->f_vol          = S.mesh.f_vol           ; 
	 S.mesh->density        = S.mesh.density           ; 
	 S.mesh.load_f_vol_e(data_user);
}
//formulation orthotrope 
if (S.mesh.type_formulation=="orthotrope") {
	 S.mesh->elastic_modulus_1=S.mesh.elastic_modulus_1;   
	 S.mesh->elastic_modulus_2=S.mesh.elastic_modulus_2;   
	 S.mesh->elastic_modulus_3=S.mesh.elastic_modulus_3;   
	 S.mesh->poisson_ratio_12 =S.mesh.poisson_ratio_12 ;   
	 S.mesh->poisson_ratio_13 =S.mesh.poisson_ratio_13 ;   
	 S.mesh->poisson_ratio_23 =S.mesh.poisson_ratio_23 ;   
	 S.mesh->shear_modulus_12 =S.mesh.shear_modulus_12 ;   
	 S.mesh->shear_modulus_13 =S.mesh.shear_modulus_13 ;   
	 S.mesh->shear_modulus_23 =S.mesh.shear_modulus_23 ;   
	 S.mesh->v1               =S.mesh.v1               ;   
	 S.mesh->v2               =S.mesh.v2               ;   
	 S.mesh->deltaT           =S.mesh.deltaT           ;   
	 S.mesh->resolution       =S.mesh.resolution       ;   
	 S.mesh->alpha_1          =S.mesh.alpha_1          ;   
	 S.mesh->alpha_2          =S.mesh.alpha_2          ;   
	 S.mesh->alpha_3          =S.mesh.alpha_3          ;   
	 S.mesh->f_vol            =S.mesh.f_vol            ;   
	 S.mesh->density        = S.mesh.density           ; 
	 S.mesh.load_f_vol_e(data_user);
}

}
#endif //ASSIGN_MATERIALS_SST_H 
