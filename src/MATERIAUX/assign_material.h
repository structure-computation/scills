#ifndef ASSIGN_MATERIALS_SST_H
#define ASSIGN_MATERIALS_SST_H

using namespace LMT;
using namespace std;

//Attention fichier generer dans : generation_auto.py
template<class SST, class SSTCARAC> void assign_material(SST &S, Vec<SSTCARAC> &matprop, Param &process){
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
//          S.mesh.load_f_vol_e();
         
	 S.mesh.type_formulation="isotrope"; 
}
}

template<class T1> void assign_material_on_element(T1 &S){
//formulation isotrope 
if (S.mesh.type_formulation=="isotrope") {
	 S.mesh->elastic_modulus= S.mesh.elastic_modulus ; 
	 S.mesh->poisson_ratio  = S.mesh.poisson_ratio   ; 
	 S.mesh->deltaT         = S.mesh.deltaT          ; 
	 S.mesh->resolution     = S.mesh.resolution      ; 
	 S.mesh->alpha          = S.mesh.alpha           ; 
	 S.mesh->f_vol          = S.mesh.f_vol           ; 
         S.mesh->density        = S.mesh.density         ;
         S.mesh.load_f_vol_e();
}

}
#endif //ASSIGN_MATERIALS_SST_H 
