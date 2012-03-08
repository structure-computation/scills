###################################################################################################
# Generation automatique des fichiers 
#     MATERIAUX / assign_materials.h : assignation des parametres materiaux
#     OPERATEURS/SST/calcul_En_Et.h : calcul des modules d'young moyen pour direction de recherche
###################################################################################################

def assign_material(formulations,target):
   output = file( target, 'w' )
   #ecriture de l'entete
   output.write('#ifndef ASSIGN_MATERIALS_SST_H\n')
   output.write('#define ASSIGN_MATERIALS_SST_H\n\n')
   output.write('using namespace LMT;\nusing namespace std;\n\n')
   
   output.write('//Attention fichier generer dans : generation_auto.py\n')
   output.write('template<class SST, class SSTCARAC> void assign_material(SST &S, Vec<SSTCARAC> &matprop, Param &process, DataUser &data_user){\n')

   for name in formulations:      
      if (name=="elasticity_isotropy_stat_Qstat"):
         #ecriture formulation isotropy
         output.write('//formulation isotrope \n')
         output.write('if (matprop[S.typmat].type=="isotrope") { \n')
         output.write('\t S.f = S.pb.formulation_elasticity_isotropy_stat_Qstat; \n')
#          output.write('\t S.f->set_mesh(S.mesh.m); \n')
         output.write('\t S.mesh.elastic_modulus = matprop[S.typmat].coef[0];    \n')
         output.write('\t S.mesh.poisson_ratio   = matprop[S.typmat].coef[1];    \n')
         output.write('\t S.mesh.deltaT          = matprop[S.typmat].coefth[1];  \n')
         output.write('\t S.mesh.resolution      = matprop[S.typmat].resolution; \n')
         output.write('\t S.mesh.alpha           = matprop[S.typmat].coefth[0];  \n')
         output.write('\t S.mesh.f_vol           = matprop[S.typmat].f_vol;      \n')
         output.write('\t S.mesh.f_vol_e         = matprop[S.typmat].f_vol_e;    \n')
         output.write('\t S.mesh.density         = matprop[S.typmat].density;    \n')
         output.write('\t S.mesh.type_formulation = "isotrope"; \n')
         output.write('}\n')
      elif (name=="elasticity_orthotropy_stat_Qstat"):
         #ecriture formulation orthotropy
         output.write('//formulation orthotrope \n')
         output.write('if (matprop[S.typmat].type=="orthotrope") {\n')
         output.write('\t S.f = S.pb.formulation_elasticity_orthotropy_stat_Qstat; \n')      
#          output.write('\t S.f->set_mesh(S.mesh.m); \n')
         output.write('\t S.mesh.elastic_modulus_1 = matprop[S.typmat].coef[0]     ; \n')
         output.write('\t S.mesh.elastic_modulus_2 = matprop[S.typmat].coef[1]     ; \n')
         output.write('\t S.mesh.elastic_modulus_3 = matprop[S.typmat].coef[2]     ; \n')
         output.write('\t S.mesh.poisson_ratio_12  = matprop[S.typmat].coef[3]     ; \n')
         output.write('\t S.mesh.poisson_ratio_13  = matprop[S.typmat].coef[4]     ; \n')
         output.write('\t S.mesh.poisson_ratio_23  = matprop[S.typmat].coef[5]     ; \n')
         output.write('\t S.mesh.shear_modulus_12  = matprop[S.typmat].coef[6]     ; \n')
         output.write('\t S.mesh.shear_modulus_13  = matprop[S.typmat].coef[7]     ; \n')
         output.write('\t S.mesh.shear_modulus_23  = matprop[S.typmat].coef[8]     ; \n')
         output.write('\t S.mesh.v1                = matprop[S.typmat].direction[0]; \n')
         output.write('\t S.mesh.v2                = matprop[S.typmat].direction[1]; \n')
         output.write('\t S.mesh.deltaT            = matprop[S.typmat].coefth[3]   ; \n')
         output.write('\t S.mesh.resolution        = matprop[S.typmat].resolution  ; \n')
         output.write('\t S.mesh.alpha_1           = matprop[S.typmat].coefth[0]   ; \n')
         output.write('\t S.mesh.alpha_2           = matprop[S.typmat].coefth[1]   ; \n')
         output.write('\t S.mesh.alpha_3           = matprop[S.typmat].coefth[2]   ; \n')
         output.write('\t S.mesh.f_vol             = matprop[S.typmat].f_vol       ; \n')
         output.write('\t S.mesh.density           = matprop[S.typmat].density     ; \n')
         output.write('\t S.mesh.f_vol_e           = matprop[S.typmat].f_vol_e     ; \n')
         output.write('\t S.mesh.type_formulation = "orthotrope"; \n')
         output.write('}\n')
      elif (name=="mesomodele"):
         #ecriture formulation orthotropy with plasticity
         output.write("//formulation du mesomodele d'endommagement\n")
         output.write('if (matprop[S.typmat].type == "mesomodele") {\n')
         output.write('\t S.f = S.pb.formulation_mesomodele; \n')
         output.write('\t S.plastique = 1; \n')
         output.write('\t S.endommageable = 3; \n')
         output.write('\t S.mesh.elastic_modulus_1 = matprop[S.typmat].coef[0]                  ; \n')
         output.write('\t S.mesh.elastic_modulus_2 = matprop[S.typmat].coef[1]                  ; \n')
         output.write('\t S.mesh.elastic_modulus_3 = matprop[S.typmat].coef[2]                  ; \n')
         output.write('\t S.mesh.poisson_ratio_12  = matprop[S.typmat].coef[3]                  ; \n')
         output.write('\t S.mesh.poisson_ratio_13  = matprop[S.typmat].coef[4]                  ; \n')
         output.write('\t S.mesh.poisson_ratio_23  = matprop[S.typmat].coef[5]                  ; \n')
         output.write('\t S.mesh.shear_modulus_12  = matprop[S.typmat].coef[6]                  ; \n')
         output.write('\t S.mesh.shear_modulus_13  = matprop[S.typmat].coef[7]                  ; \n')
         output.write('\t S.mesh.shear_modulus_23  = matprop[S.typmat].coef[8]                  ; \n')
         output.write('\t S.mesh.v1                = matprop[S.typmat].direction[0]             ; \n')
         output.write('\t S.mesh.v2                = matprop[S.typmat].direction[1]             ; \n')
         output.write('\t S.mesh.deltaT            = matprop[S.typmat].coefth[3]                ; \n')
         output.write('\t S.mesh.resolution        = matprop[S.typmat].resolution               ; \n')
         output.write('\t S.mesh.alpha_1           = matprop[S.typmat].coefth[0]                ; \n')
         output.write('\t S.mesh.alpha_2           = matprop[S.typmat].coefth[1]                ; \n')
         output.write('\t S.mesh.alpha_3           = matprop[S.typmat].coefth[2]                ; \n')
         output.write('\t S.mesh.k_p               = matprop[S.typmat].coefp[0]                 ; \n')
         output.write('\t S.mesh.m                 = matprop[S.typmat].coefp[1]                 ; \n')
         output.write('\t S.mesh.R0                = matprop[S.typmat].coefp[2]                 ; \n')
         output.write('\t S.mesh.couplage          = matprop[S.typmat].coefp[3]                 ; \n')
         output.write('\t S.mesh.Yo                = matprop[S.typmat].param_damage.Yo          ; \n')
         output.write('\t S.mesh.Yc                = matprop[S.typmat].param_damage.Yc          ; \n')
         output.write('\t S.mesh.Ycf               = matprop[S.typmat].param_damage.Ycf         ; \n')
         output.write('\t S.mesh.dmax              = matprop[S.typmat].param_damage.dmax        ; \n')
         output.write('\t S.mesh.effet_retard      = matprop[S.typmat].param_damage.effet_retard; \n')
         output.write('\t if (matprop[S.typmat].param_damage.effet_retard) {\n')
         output.write('\t \t S.mesh.tau_c = matprop[S.typmat].param_damage.tau; \n')
         output.write('\t \t S.mesh.a     = matprop[S.typmat].param_damage.a  ; \n')
         output.write('\t \t S.mesh.b_c   = matprop[S.typmat].param_damage.b_c; \n')
         output.write('\t }\n')
         output.write('\t S.mesh.type_formulation = "mesomodele"; \n')
         output.write('}\n')
      else : print "formulation non definie"
                      
#    output.write('assign_material_on_element(S);\n')
   output.write('}\n')
   output.write('\n')
   
   output.write('template<class T1> void assign_material_on_element(T1 &S, DataUser &data_user){\n')

   for name in formulations:      
      if (name=="elasticity_isotropy_stat_Qstat"):
         #ecriture formulation isotropy
         output.write('//formulation isotrope \n')
         output.write('if (S.mesh.type_formulation=="isotrope") {\n')
         output.write('\t S.mesh->elastic_modulus= S.mesh.elastic_modulus ; \n')
         output.write('\t S.mesh->poisson_ratio  = S.mesh.poisson_ratio   ; \n')
         output.write('\t S.mesh->deltaT         = S.mesh.deltaT          ; \n')
         output.write('\t S.mesh->resolution     = S.mesh.resolution      ; \n')
         output.write('\t S.mesh->alpha          = S.mesh.alpha           ; \n')
         output.write('\t S.mesh->f_vol          = S.mesh.f_vol           ; \n')
         output.write('\t S.mesh->density        = S.mesh.density         ; \n')
         output.write('\t S.mesh.load_f_vol_e(data_user); \n')
#          output.write('\t cout <<   S.mesh.elastic_modulus << " " <<   S.mesh.poisson_ratio   << " " <<   S.mesh.deltaT  << " " << S.mesh.resolution  << " " << S.mesh.alpha  << " " << S.mesh.f_vol ; \n')    
         output.write('}\n')                      
      elif (name=="elasticity_orthotropy_stat_Qstat"):
         #ecriture formulation orthotropy         
         output.write('//formulation orthotrope \n')
         output.write('if (S.mesh.type_formulation=="orthotrope") {\n')
         output.write('\t S.mesh->elastic_modulus_1 = S.mesh.elastic_modulus_1; \n')
         output.write('\t S.mesh->elastic_modulus_2 = S.mesh.elastic_modulus_2; \n')
         output.write('\t S.mesh->elastic_modulus_3 = S.mesh.elastic_modulus_3; \n')
         output.write('\t S.mesh->poisson_ratio_12  = S.mesh.poisson_ratio_12 ; \n')
         output.write('\t S.mesh->poisson_ratio_13  = S.mesh.poisson_ratio_13 ; \n')
         output.write('\t S.mesh->poisson_ratio_23  = S.mesh.poisson_ratio_23 ; \n')
         output.write('\t S.mesh->shear_modulus_12  = S.mesh.shear_modulus_12 ; \n')
         output.write('\t S.mesh->shear_modulus_13  = S.mesh.shear_modulus_13 ; \n')
         output.write('\t S.mesh->shear_modulus_23  = S.mesh.shear_modulus_23 ; \n')
         output.write('\t S.mesh->v1                = S.mesh.v1               ; \n')
         output.write('\t S.mesh->v2                = S.mesh.v2               ; \n')
         output.write('\t S.mesh->deltaT            = S.mesh.deltaT           ; \n')
         output.write('\t S.mesh->resolution        = S.mesh.resolution       ; \n')
         output.write('\t S.mesh->alpha_1           = S.mesh.alpha_1          ; \n')
         output.write('\t S.mesh->alpha_2           = S.mesh.alpha_2          ; \n')
         output.write('\t S.mesh->alpha_3           = S.mesh.alpha_3          ; \n')
         output.write('\t S.mesh->f_vol             = S.mesh.f_vol            ; \n')
         output.write('\t S.mesh->density           = S.mesh.density          ; \n')
         output.write('\t S.mesh.load_f_vol_e(data_user); \n')
         output.write('}\n')
      elif (name=="mesomodele"):
         #ecriture formulation orthotropy with plasticity
         output.write("//formulation du mesomodele d'endommagement\n")
         output.write('if (S.mesh.type_formulation=="mesomodele") {\n')
         output.write('\t S.mesh->elastic_modulus_1 = S.mesh.elastic_modulus_1; \n')
         output.write('\t S.mesh->elastic_modulus_2 = S.mesh.elastic_modulus_2; \n')
         output.write('\t S.mesh->elastic_modulus_3 = S.mesh.elastic_modulus_3; \n')
         output.write('\t S.mesh->poisson_ratio_12 = S.mesh.poisson_ratio_12  ; \n')
         output.write('\t S.mesh->poisson_ratio_13 = S.mesh.poisson_ratio_13  ; \n')
         output.write('\t S.mesh->poisson_ratio_23 = S.mesh.poisson_ratio_23  ; \n')
         output.write('\t S.mesh->shear_modulus_12 = S.mesh.shear_modulus_12  ; \n')
         output.write('\t S.mesh->shear_modulus_13 = S.mesh.shear_modulus_13  ; \n')
         output.write('\t S.mesh->shear_modulus_23 = S.mesh.shear_modulus_23  ; \n')
         output.write('\t S.mesh->v1               = S.mesh.v1                ; \n')
         output.write('\t S.mesh->v2               = S.mesh.v2                ; \n')
         output.write('\t S.mesh->deltaT           = S.mesh.deltaT            ; \n')
         output.write('\t S.mesh->resolution       = S.mesh.resolution        ; \n')
         output.write('\t S.mesh->alpha_1          = S.mesh.alpha_1           ; \n')
         output.write('\t S.mesh->alpha_2          = S.mesh.alpha_2           ; \n')
         output.write('\t S.mesh->alpha_3          = S.mesh.alpha_3           ; \n')
         output.write('\t S.mesh->k_p              = S.mesh.k_p               ; \n')
         output.write('\t S.mesh->m                = S.mesh.m                 ; \n')
         output.write('\t S.mesh->R0               = S.mesh.R0                ; \n')
         output.write('\t S.mesh->couplage         = S.mesh.couplage          ; \n')
         output.write('\t S.mesh->Yo               = S.mesh.Yo                ; \n')
         output.write('\t S.mesh->Yc               = S.mesh.Yc                ; \n')
         output.write('\t S.mesh->Ycf              = S.mesh.Ycf               ; \n')
         output.write('\t S.mesh->dmax             = S.mesh.dmax              ; \n')
         output.write('\t S.mesh->effet_retard     = S.mesh.effet_retard      ; \n')
         output.write('\t if (S.mesh.effet_retard) {\n')
         output.write('\t \t S.mesh->tau_c = S.mesh.tau_c; \n')
         output.write('\t \t S.mesh->a     = S.mesh.a    ; \n')
         output.write('\t \t S.mesh->b_c   = S.mesh.b_c  ; \n')
         output.write('\t }\n')
         output.write('\t S.mesh.load_f_vol_e(data_user); \n')
         output.write('}\n')
      else : print "formulation non definie"
                      
   output.write('\n')
   output.write('}\n')

   
   output.write('#endif //ASSIGN_MATERIALS_SST_H \n')

def create_new_elem_p(self,target):
   output = file( target, 'w' )
   #ecriture de l'entete
   output.write('#include "mesh/triangle.h"\n')
   output.write('#include "mesh/quad.h"\n')
   output.write('#include "mesh/quad_8.h"\n')
   output.write('#include "mesh/triangle_6.h"\n')
   output.write('#include "mesh/tetra.h"\n')
   output.write('#include "mesh/tetra_10.h"\n')
   
   output.write('namespace LMT{\n');
   output.write('//Attention fichier generer dans : generation_auto.py\n')
   output.write('///Fonction generique\n')
   output.write('template<class TE, class TM> void create_new_elem_p(TE &e, TM &m, Vec<unsigned> &rep_nodes){}\n')
   output.write('///Fonctions specialisees\n')
   
   for name in self.elements:
      if (name=="Triangle"):
         output.write('template<class TNB,class TN,class TD,unsigned NET, class TM> void create_new_elem_p(Element<Triangle,TNB,TN,TD,NET> &e, TM &m, Vec<unsigned> &rep_nodes){ m.add_element(Triangle_6(),DefaultBehavior(),rep_nodes.ptr() ); } \n')
      elif(name=="Quad"):
         output.write('template<class TNB,class TN,class TD,unsigned NET, class TM> void create_new_elem_p(Element<Quad,TNB,TN,TD,NET> &e, TM &m, Vec<unsigned> &rep_nodes){ m.add_element(Quad_8(),DefaultBehavior(),rep_nodes.ptr() ); } \n')
      elif(name=="Tetra"):
         output.write('template<class TNB,class TN,class TD,unsigned NET, class TM> void create_new_elem_p(Element<Tetra,TNB,TN,TD,NET> &e, TM &m, Vec<unsigned> &rep_nodes){ m.add_element(Tetra_10(),DefaultBehavior(),rep_nodes.ptr() ); } \n')
      elif(name=="Wedge"):
         output.write('template<class TNB,class TN,class TD,unsigned NET, class TM> void create_new_elem_p(Element<Wedge,TNB,TN,TD,NET> &e, TM &m, Vec<unsigned> &rep_nodes){ m.add_element(Wedge_15(),DefaultBehavior(),rep_nodes.ptr() ); } \n')
      elif(name=="Hexa"):
         output.write('template<class TNB,class TN,class TD,unsigned NET, class TM> void create_new_elem_p(Element<Hexa,TNB,TN,TD,NET> &e, TM &m, Vec<unsigned> &rep_nodes){ m.add_element(Hexa_20(),DefaultBehavior(),rep_nodes.ptr() ); } \n')
   
   output.write('};\n')   
   
def calcul_En_Et(formulations,target):
   output = file( target, 'w' )
   #ecriture de l'entete
   output.write('#ifndef CALCUL_EN_ET_H\n')
   output.write('#define CALCUL_EN_ET_H\n\n')
   output.write('using namespace LMT;\nusing namespace std;\n\n')
   
   output.write('//Attention fichier generer dans : generation_auto.py\n')
   output.write('template<class SST> void calcul_En_Et(SST &S, typename SST::T &En, typename SST::T &Et){\n')

   for name in formulations:   
      #ecriture directions isotropy             
      if (name=="elasticity_isotropy_stat_Qstat" or name=="elasticity_viscosity_Qstat"):
         output.write('//formulation isotrope ou viscoelastique \n')
         output.write('if (S.f->get_name()=="elasticity_isotropy_stat_Qstat" or S.f->get_name()=="elasticity_viscosity_Qstat"){\n')
         output.write('\t En = S.mesh.elastic_modulus;\n')
         output.write('\t Et = S.mesh.elastic_modulus;\n')
         output.write('}\n')
      #ecriture directions orthotropy  
      elif (name=="elasticity_orthotropy_stat_Qstat" or name=="elasticity_orthotropy_damage_stat_Qstat"):
         output.write('if(S.f->get_name()=="elasticity_orthotropy_stat_Qstat" or S.f->get_name()=="elasticity_orthotropy_damage_stat_Qstat") {\n')
         output.write('\t En = (S.mesh.elastic_modulus_1 + 2.*S.mesh.elastic_modulus_2)/3.;\n')
         output.write('\t Et = En;\n')
         output.write('}\n')
      else : print "Formulation non definie"
   
   output.write('}\n')
   output.write('#endif //CALCUL_EN_ET_H \n')     

def generation_auto(formulations):
   assign_material(formulations,'src/MATERIAUX/assign_material.h')
   calcul_En_Et(formulations,'src/OPERATEURS/SST/calcul_En_Et.h')
   print "Generation auto : ok" 
   