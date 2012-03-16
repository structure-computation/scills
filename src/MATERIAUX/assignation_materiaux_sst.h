using namespace LMT;
using namespace std;


//********************************************
// Affectation de materiaux aux SST
//******************************************** 
/**\ingroup Materiaux
 \brief Structure contenant un operateur() (utilisable avec apply) permettant l'affectation des valeurs des proprietes materiaux au maillage associé à la sous-structure.
 
 Chaque sous-structure comporte un champ typmat (Sst::typmat) dont la valeur correspond à une composante du vecteur des propriétés matériaux SstCaract. 
 Selon le nom correspondant extrait du vecteur de propriétés matériaux, on déclare une formulation du type donné. 
 Si certaines formulations ne veulent pas être inclus dans la base des formulations disponibles, il suffit de commenter chaque bloc if et d'enlever la formulation correspondante dans le fichier SConstruct.  
 */
struct assignation_material_to_SST{
    template<class SST, class MatProps> void operator()(SST &S,MatProps &matprops) const{ 
        S.matprop = matprops[S.typmat]; 
        //formulation isotrope 
        if (matprops[S.typmat].type=="isotrope") {
            S.f = S.pb.formulation_elasticity_isotropy_stat_Qstat;
            S.mesh.type_formulation="isotrope"; 
        }
        //formulation orthotrope 
        if (matprops[S.typmat].type=="orthotrope") {
            S.f = S.pb.formulation_elasticity_orthotropy_stat_Qstat;
            S.mesh.type_formulation="orthotrope"; 
        }
        //formulation mesomodele 
        if (matprops[S.typmat].type=="mesomodele") {
            S.f = S.pb.formulation_mesomodele;
            S.mesh.type_formulation="mesomodele"; 
        }
   }
};


template<class T1> void assign_material_on_element(T1 &S, DataUser &data_user){
    //formulation isotrope 
    if (S.mesh.type_formulation=="isotrope") {
        S.mesh->elastic_modulus = S.matprop.elastic_modulus ; 
        S.mesh->poisson_ratio   = S.matprop.poisson_ratio   ; 
        S.mesh->deltaT          = S.matprop.deltaT          ; 
        S.mesh->resolution      = S.matprop.resolution      ; 
        S.mesh->alpha           = S.matprop.alpha           ; 
        S.mesh->f_vol           = S.matprop.f_vol           ; 
        S.mesh->density         = S.matprop.density         ; 
        S.mesh.load_f_vol_e(S.matprop.f_vol_e,data_user);
    }
    //formulation orthotrope 
    if (S.mesh.type_formulation=="orthotrope") {
        S.mesh->elastic_modulus_1 = S.matprop.elastic_modulus_1;
        S.mesh->elastic_modulus_2 = S.matprop.elastic_modulus_2;
        S.mesh->elastic_modulus_3 = S.matprop.elastic_modulus_3;
        S.mesh->poisson_ratio_12  = S.matprop.poisson_ratio_12 ;
        S.mesh->poisson_ratio_13  = S.matprop.poisson_ratio_13 ;
        S.mesh->poisson_ratio_23  = S.matprop.poisson_ratio_23 ;
        S.mesh->shear_modulus_12  = S.matprop.shear_modulus_12 ;
        S.mesh->shear_modulus_13  = S.matprop.shear_modulus_13 ;
        S.mesh->shear_modulus_23  = S.matprop.shear_modulus_23 ;
        S.mesh->v1                = S.matprop.v1               ;
        S.mesh->v2                = S.matprop.v2               ;
        S.mesh->deltaT            = S.matprop.deltaT           ;
        S.mesh->resolution        = S.matprop.resolution       ;
        S.mesh->alpha_1           = S.matprop.alpha_1          ;
        S.mesh->alpha_2           = S.matprop.alpha_2          ;
        S.mesh->alpha_3           = S.matprop.alpha_3          ;
        S.mesh->f_vol             = S.matprop.f_vol            ;
        S.mesh->density           = S.matprop.density          ;
        S.mesh.load_f_vol_e(S.matprop.f_vol_e,data_user);
    }
    //formulation mesomodele 
    if (S.mesh.type_formulation=="mesomodele") {
        S.mesh->elastic_modulus_1 = S.matprop.elastic_modulus_1;
        S.mesh->elastic_modulus_2 = S.matprop.elastic_modulus_2;
        S.mesh->elastic_modulus_3 = S.matprop.elastic_modulus_3;
        S.mesh->poisson_ratio_12  = S.matprop.poisson_ratio_12 ;
        S.mesh->poisson_ratio_13  = S.matprop.poisson_ratio_13 ;
        S.mesh->poisson_ratio_23  = S.matprop.poisson_ratio_23 ;
        S.mesh->shear_modulus_12  = S.matprop.shear_modulus_12 ;
        S.mesh->shear_modulus_13  = S.matprop.shear_modulus_13 ;
        S.mesh->shear_modulus_23  = S.matprop.shear_modulus_23 ;
        S.mesh->v1                = S.matprop.v1               ;
        S.mesh->v2                = S.matprop.v2               ;
        S.mesh->deltaT            = S.matprop.deltaT           ;
        S.mesh->resolution        = S.matprop.resolution       ;
        S.mesh->alpha_1           = S.matprop.alpha_1          ;
        S.mesh->alpha_2           = S.matprop.alpha_2          ;
        S.mesh->alpha_3           = S.matprop.alpha_3          ;
        S.mesh->f_vol             = S.matprop.f_vol            ;
        S.mesh->density           = S.matprop.density          ;
        
        S.mesh->k_p      = S.matprop.k_p;
        S.mesh->m_p      = S.matprop.m_p;
        S.mesh->R0       = S.matprop.R0;
        S.mesh->couplage = S.matprop.coefvm_composite;
        
        S.mesh->Yo           = S.matprop.Yo;
        S.mesh->Yc           = S.matprop.Yc;
        S.mesh->Ycf          = S.matprop.Ycf;
        S.mesh->dmax         = S.matprop.dmax;
        S.mesh->b_c          = S.matprop.b_c;
        S.mesh->effet_retard = S.matprop.effet_retard;
        S.mesh->a            = S.matprop.a;
        S.mesh->tau_c        = S.matprop.tau_c;
        
        S.mesh.load_f_vol_e(S.matprop.f_vol_e,data_user);
    }
}
