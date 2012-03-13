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
   template<class SST, class TV3> void operator()(SST &S,TV3 &matprop) const{
       //formulation isotrope 
       if (matprop[S.typmat].type=="isotrope") {
           S.f = S.pb.formulation_elasticity_isotropy_stat_Qstat; 
           S.mesh.elastic_modulus = matprop[S.typmat].coef[0]; 
           S.mesh.poisson_ratio   = matprop[S.typmat].coef[1]; 
           S.mesh.deltaT          = matprop[S.typmat].coefth[1];  
           S.mesh.resolution      = matprop[S.typmat].resolution; 
           S.mesh.alpha           = matprop[S.typmat].coefth[0];  
           S.mesh.f_vol           = matprop[S.typmat].f_vol;  
           S.mesh.f_vol_e         = matprop[S.typmat].f_vol_e;  
           S.mesh.density         = matprop[S.typmat].density; 
           S.mesh.type_formulation="isotrope"; 
       }
       //formulation orthotrope 
       if (matprop[S.typmat].type=="orthotrope") {
           S.f = S.pb.formulation_elasticity_orthotropy_stat_Qstat;
           S.mesh.elastic_modulus_1 = matprop[S.typmat].coef[0];
           S.mesh.elastic_modulus_2 = matprop[S.typmat].coef[1];
           S.mesh.elastic_modulus_3 = matprop[S.typmat].coef[2];
           S.mesh.poisson_ratio_12  = matprop[S.typmat].coef[3];
           S.mesh.poisson_ratio_13  = matprop[S.typmat].coef[4];
           S.mesh.poisson_ratio_23  = matprop[S.typmat].coef[5];
           S.mesh.shear_modulus_12  = matprop[S.typmat].coef[6];
           S.mesh.shear_modulus_13  = matprop[S.typmat].coef[7];
           S.mesh.shear_modulus_23  = matprop[S.typmat].coef[8];
           S.mesh.v1                = matprop[S.typmat].direction[0];  
           S.mesh.v2                = matprop[S.typmat].direction[1];
           S.mesh.deltaT            = matprop[S.typmat].coefth[3];
           S.mesh.resolution        = matprop[S.typmat].resolution;
           S.mesh.alpha_1           = matprop[S.typmat].coefth[0];
           S.mesh.alpha_2           = matprop[S.typmat].coefth[1];
           S.mesh.alpha_3           = matprop[S.typmat].coefth[2];
           S.mesh.f_vol             = matprop[S.typmat].f_vol;
           S.mesh.density           = matprop[S.typmat].density;
           S.mesh.f_vol_e           = matprop[S.typmat].f_vol_e;
           S.mesh.type_formulation="orthotrope"; 
       }
       //formulation mesomodele 
       if (matprop[S.typmat].type=="mesomodele") {
           S.f = S.pb.formulation_elasticity_orthotropy_stat_Qstat;
           S.mesh.elastic_modulus_1 = matprop[S.typmat].coef[0];
           S.mesh.elastic_modulus_2 = matprop[S.typmat].coef[1];
           S.mesh.elastic_modulus_3 = matprop[S.typmat].coef[2];
           S.mesh.poisson_ratio_12  = matprop[S.typmat].coef[3];
           S.mesh.poisson_ratio_13  = matprop[S.typmat].coef[4];
           S.mesh.poisson_ratio_23  = matprop[S.typmat].coef[5];
           S.mesh.shear_modulus_12  = matprop[S.typmat].coef[6];
           S.mesh.shear_modulus_13  = matprop[S.typmat].coef[7];
           S.mesh.shear_modulus_23  = matprop[S.typmat].coef[8];
           S.mesh.v1                = matprop[S.typmat].direction[0];  
           S.mesh.v2                = matprop[S.typmat].direction[1];
           S.mesh.deltaT            = matprop[S.typmat].coefth[3];
           S.mesh.resolution        = matprop[S.typmat].resolution;
           S.mesh.alpha_1           = matprop[S.typmat].coefth[0];
           S.mesh.alpha_2           = matprop[S.typmat].coefth[1];
           S.mesh.alpha_3           = matprop[S.typmat].coefth[2];
           S.mesh.f_vol             = matprop[S.typmat].f_vol;
           S.mesh.density           = matprop[S.typmat].density;
           S.mesh.f_vol_e           = matprop[S.typmat].f_vol_e;
           
           S.mesh.k_p      = matprop[S.typmat].coefp[0];
           S.mesh.m_p      = matprop[S.typmat].coefp[1];
           S.mesh.R0       = matprop[S.typmat].coefp[2];
           S.mesh.couplage = matprop[S.typmat].coefp[3];
           
           S.mesh.Yo           = matprop[S.typmat].coefendom[0];
           S.mesh.Yc           = matprop[S.typmat].coefendom[1];
           S.mesh.Ycf          = matprop[S.typmat].coefendom[2];
           S.mesh.dmax         = matprop[S.typmat].coefendom[3];
           S.mesh.b_c          = matprop[S.typmat].coefendom[4];
           S.mesh.effet_retard = matprop[S.typmat].coefendom[5];
           S.mesh.a            = matprop[S.typmat].coefendom[6];
           S.mesh.tau_c        = matprop[S.typmat].coefendom[7];
           
           S.mesh.type_formulation="mesomodele"; 
       }
   }
};


template<class T1> void assign_material_on_element(T1 &S, DataUser &data_user){
    //formulation isotrope 
    if (S.mesh.type_formulation=="isotrope") {
        S.mesh->elastic_modulus = S.mesh.elastic_modulus ; 
        S.mesh->poisson_ratio   = S.mesh.poisson_ratio   ; 
        S.mesh->deltaT          = S.mesh.deltaT          ; 
        S.mesh->resolution      = S.mesh.resolution      ; 
        S.mesh->alpha           = S.mesh.alpha           ; 
        S.mesh->f_vol           = S.mesh.f_vol           ; 
        S.mesh->density         = S.mesh.density           ; 
        S.mesh.load_f_vol_e(data_user);
    }
    //formulation orthotrope 
    if (S.mesh.type_formulation=="orthotrope") {
        S.mesh->elastic_modulus_1 = S.mesh.elastic_modulus_1;
        S.mesh->elastic_modulus_2 = S.mesh.elastic_modulus_2;
        S.mesh->elastic_modulus_3 = S.mesh.elastic_modulus_3;
        S.mesh->poisson_ratio_12  = S.mesh.poisson_ratio_12 ;
        S.mesh->poisson_ratio_13  = S.mesh.poisson_ratio_13 ;
        S.mesh->poisson_ratio_23  = S.mesh.poisson_ratio_23 ;
        S.mesh->shear_modulus_12  = S.mesh.shear_modulus_12 ;
        S.mesh->shear_modulus_13  = S.mesh.shear_modulus_13 ;
        S.mesh->shear_modulus_23  = S.mesh.shear_modulus_23 ;
        S.mesh->v1                = S.mesh.v1               ;
        S.mesh->v2                = S.mesh.v2               ;
        S.mesh->deltaT            = S.mesh.deltaT           ;
        S.mesh->resolution        = S.mesh.resolution       ;
        S.mesh->alpha_1           = S.mesh.alpha_1          ;
        S.mesh->alpha_2           = S.mesh.alpha_2          ;
        S.mesh->alpha_3           = S.mesh.alpha_3          ;
        S.mesh->f_vol             = S.mesh.f_vol            ;
        S.mesh->density           = S.mesh.density          ;
        S.mesh.load_f_vol_e(data_user);
    }
    //formulation mesomodele 
    if (S.mesh.type_formulation=="mesomodele") {
        S.mesh->elastic_modulus_1 = S.mesh.elastic_modulus_1;
        S.mesh->elastic_modulus_2 = S.mesh.elastic_modulus_2;
        S.mesh->elastic_modulus_3 = S.mesh.elastic_modulus_3;
        S.mesh->poisson_ratio_12  = S.mesh.poisson_ratio_12 ;
        S.mesh->poisson_ratio_13  = S.mesh.poisson_ratio_13 ;
        S.mesh->poisson_ratio_23  = S.mesh.poisson_ratio_23 ;
        S.mesh->shear_modulus_12  = S.mesh.shear_modulus_12 ;
        S.mesh->shear_modulus_13  = S.mesh.shear_modulus_13 ;
        S.mesh->shear_modulus_23  = S.mesh.shear_modulus_23 ;
        S.mesh->v1                = S.mesh.v1               ;
        S.mesh->v2                = S.mesh.v2               ;
        S.mesh->deltaT            = S.mesh.deltaT           ;
        S.mesh->resolution        = S.mesh.resolution       ;
        S.mesh->alpha_1           = S.mesh.alpha_1          ;
        S.mesh->alpha_2           = S.mesh.alpha_2          ;
        S.mesh->alpha_3           = S.mesh.alpha_3          ;
        S.mesh->f_vol             = S.mesh.f_vol            ;
        S.mesh->density           = S.mesh.density          ;
        
        /*S.mesh->k_p    = S.mesh.k_p;
        S.mesh->m_p      = S.mesh.m_p;
        S.mesh->R0       = S.mesh.R0;
        S.mesh->couplage = S.mesh.couplage;
        
        S.mesh->Yo           = S.mesh.Yo;
        S.mesh->Yc           = S.mesh.Yc;
        S.mesh->Ycf          = S.mesh.Ycf;
        S.mesh->dmax         = S.mesh.dmax;
        S.mesh->b_c          = S.mesh.b_c;
        S.mesh->effet_retard = S.mesh.effet_retard;
        S.mesh->a            = S.mesh.a;
        S.mesh->tau_c        = S.mesh.tau_c;
        
        S.mesh.load_f_vol_e(data_user);*/
    }
}
