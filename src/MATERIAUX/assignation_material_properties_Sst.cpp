#include "assignation_material_properties_Sst.h"

#include "../../LMT/include/containers/mat.h"
#include "../../LMT/include/codegen/codegen.h"
#include "../../LMT/include/containers/basicops.h"

using namespace Codegen;
using namespace LMT;

void assignation_materials_property_SST(DataUser &data_user, Vec<SstCarac> &matprops, Vec<Sst> &S, Process &process, FieldStructureUser &field_structure_user){
    if (process.rank == 0) std::cout << "\t Assignation du materiau aux SST" << std::endl;
    // lecture des proprietes materiau des ssts
    
    BasicVec<BasicVec<TYPEREEL > > mat_prop_temp;
    read_matprop(matprops,process,data_user, mat_prop_temp); 
    //assignation des proprietes aux SST
    apply_mt(S,process.nb_threads,assignation_material_to_SST(),matprops,process.plasticite);      
    //assignation des proprietes aux group_elements (pour GPU)
    field_structure_user.assign_material_id_to_group_elements(data_user);
    field_structure_user.assign_material_properties_to_group_elements(data_user,mat_prop_temp);
}

void read_matprop(Vec<SstCarac> &matprops, Process &process, DataUser &data_user , BasicVec<BasicVec<TYPEREEL> > &mat_prop_temp) {
    unsigned nbmat = data_user.behaviour_materials.size();
    matprops.resize(nbmat);
    mat_prop_temp.resize(nbmat);
    for(unsigned i=0;i<nbmat;++i) {
        matprops[i].id = data_user.behaviour_materials[i].id;
        matprops[i].type_num = data_user.behaviour_materials[i].type_num;
        matprops[i].type = data_user.behaviour_materials[i].type;
        matprops[i].type_plast = data_user.behaviour_materials[i].type_plast;
        matprops[i].type_endo = data_user.behaviour_materials[i].type_endo;
        matprops[i].comp = data_user.behaviour_materials[i].comp;
        if(data_user.dim == 2){
            if (data_user.behaviour_materials[i].resolution =="CP")
                matprops[i].resolution=1;
            else if (data_user.behaviour_materials[i].resolution =="DP")
                matprops[i].resolution=0;
            else {
                std::cout << "type de resolution non implemente : choix contrainte_plane ou deformation_plane" << std::endl;
                assert(0);
            }
        }else{
            matprops[i].resolution=0;
        }
        
        std::vector<Ex> symbols;
        if (DIM==2) {
            symbols.push_back("x");
            symbols.push_back("y");
        }
        else if (DIM==3) {
            symbols.push_back("x");
            symbols.push_back("y");
            symbols.push_back("z");
        }
        
        if(data_user.options.Multiresolution_on==1){
            //ajout des variables de multiresolution aux symboles
            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
                Sc2String var="V"; var<<i_par; //nom de la variable de multiresolution
                symbols.push_back(var.c_str());
            }
        }
        
        Ex Pi = symbol("PI");
        symbols.push_back(Pi);
        
        Vec<Sc2String> vstr;
        vstr.resize(data_user.dim, "0");
        
        for(int i_fvol=0; i_fvol<data_user.behaviour_bc_volume.size(); i_fvol++){
            if(data_user.behaviour_bc_volume[i_fvol].select){
                for(int d=0; d<data_user.dim; d++){
                    vstr[d] += " + " + data_user.behaviour_bc_volume[i_fvol].step[0].CLv_step_prop[d] + " * " + data_user.behaviour_bc_volume[i_fvol].step[0].CLv_step_prop[6];
                }
            }
        }
        if(process.rank == 0){
            std::cout << "force volumique 0 = " << vstr[0] << std::endl;
            std::cout << "force volumique 1 = " << vstr[1] << std::endl;
#if DIM == 3
            std::cout << "force volumique 2 = " << vstr[2] << std::endl;
#endif
        }
        
        Vec<Ex> expr;
        expr.resize(DIM);
        for(unsigned d2=0;d2<DIM;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
            expr[d2] = read_ex(vstr[d2],symbols);
        }
        
        Vec<double,DIM> data;
        Ex::MapExNum var;
        for(unsigned d2=0;d2<DIM;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
            var[symbols[d2]]= 0.;
        }
        if(data_user.options.Multiresolution_on==1){
            //evaluation des variables de multiresolution aux symboles
            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                var[symbols[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
        }
        var[symbols[DIM+data_user.Multiresolution_parameters.size()]]=M_PI;
        
        for(unsigned d2=0;d2<DIM;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
            data[d2] = (double)expr[d2].subs_numerical(var);
        
        matprops[i].f_vol_e=vstr;
        matprops[i].f_vol=data;
        
        
        mat_prop_temp[i].resize(data_user.behaviour_materials[i].mat_prop.size());
        for(int i_prop=0; i_prop<data_user.behaviour_materials[i].mat_prop.size(); i_prop++){
            if(data_user.behaviour_materials[i].mat_prop[i_prop] == ""){
                data_user.behaviour_materials[i].mat_prop[i_prop] = "0";
            }
            Ex expr_temp;
            expr_temp = read_ex(data_user.behaviour_materials[i].mat_prop[i_prop],symbols);
            Ex::MapExNum var_temp;
            for(unsigned d2=0;d2<DIM;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                var_temp[symbols[d2]]= 0.;
            }
            if(data_user.options.Multiresolution_on==1){
                //evaluation des variables de multiresolution aux symboles
                for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                    var_temp[symbols[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
            }
            var_temp[symbols[DIM+data_user.Multiresolution_parameters.size()]]=M_PI;
            
            mat_prop_temp[i][i_prop] = (TYPEREEL) expr_temp.subs_numerical(var_temp);
        }
        matprops[i].density = mat_prop_temp[i][3];
        
        std::cout << "id : " << matprops[i].id << "    type : " << matprops[i].type << "    comp : " << matprops[i].comp << std::endl;
        if(matprops[i].type.find("isotrope")<matprops[i].type.size()) {/// comportement thermo-elastique isotrope
            matprops[i].elastic_modulus = mat_prop_temp[i][0];   /// E
            matprops[i].poisson_ratio   = mat_prop_temp[i][1];   /// nu
            matprops[i].alpha           = mat_prop_temp[i][2];   /// alpha
            matprops[i].deltaT          = 0;                     /// deltaT
        } else if (matprops[i].type.find("orthotrope")<matprops[i].type.size()) {/// comportement thermo-elastique orthotrope
            /// coefficients d'elasticite
            matprops[i].elastic_modulus_1 = mat_prop_temp[i][14];   /// E1
            matprops[i].elastic_modulus_2 = mat_prop_temp[i][15];   /// E2
            matprops[i].elastic_modulus_3 = mat_prop_temp[i][16];   /// E3
            
            matprops[i].poisson_ratio_12 = mat_prop_temp[i][20];   /// nu12
            matprops[i].poisson_ratio_13 = mat_prop_temp[i][22];   /// nu13
            matprops[i].poisson_ratio_23 = mat_prop_temp[i][21];   /// nu23
            
            matprops[i].shear_modulus_12 = mat_prop_temp[i][17];   /// G12
            matprops[i].shear_modulus_13 = mat_prop_temp[i][19];   /// G13
            matprops[i].shear_modulus_23 = mat_prop_temp[i][18];   /// G23
            
            /// directions d'orthotropie
            for(int d=0; d<data_user.dim; d++){
                matprops[i].v1[d]=mat_prop_temp[i][d+5];
                matprops[i].v2[d]=mat_prop_temp[i][d+8];
            }
            matprops[i].v1=matprops[i].v1/norm_2(matprops[i].v1);
            matprops[i].v2=matprops[i].v2/norm_2(matprops[i].v2);
            
            /// coefficients thermiques
            matprops[i].alpha_1 = mat_prop_temp[i][23];      ///alpha_1
            matprops[i].alpha_2 = mat_prop_temp[i][24];      ///alpha_2
            matprops[i].alpha_3 = mat_prop_temp[i][25];      ///alpha_3
            matprops[i].deltaT  = 0;                         /// deltaT
        }
        
        if(matprops[i].comp.find("visqueux")<matprops[i].comp.size()){/// Comportement visqueux
            matprops[i].viscosite       = mat_prop_temp[i][4];   /// viscosite
        }
        if (matprops[i].comp.find("plastique")<matprops[i].comp.size() or matprops[i].comp.find("mesomodele")<matprops[i].comp.size()) {
            /// parametres de plasticite
            matprops[i].plast_ecrouissage_init = mat_prop_temp[i][26];     /// R0
            matprops[i].plast_ecrouissage_mult = mat_prop_temp[i][27];     /// k_p
            matprops[i].plast_ecrouissage_expo = mat_prop_temp[i][28];     /// m_p
            matprops[i].plast_cinematique_coef = mat_prop_temp[i][29];     /// couplage
        }
        if (matprops[i].comp.find("mesomodele")<matprops[i].comp.size()) {
            /// parametres d'endommagement
            matprops[i].Yo           = mat_prop_temp[i][30];     /// Yo
            matprops[i].Yc           = mat_prop_temp[i][31];     /// Yc
            matprops[i].Ycf          = mat_prop_temp[i][32];     /// Ycf
            matprops[i].dmax         = mat_prop_temp[i][33];     /// dmax
            matprops[i].b_c          = mat_prop_temp[i][34];     /// b_c
            matprops[i].effet_retard = mat_prop_temp[i][35];     /// effet_retard
            matprops[i].a            = mat_prop_temp[i][36];     /// a
            matprops[i].tau_c        = mat_prop_temp[i][37];     /// tau_c
        }
    }
};

void assignation_material_to_SST::operator()(Sst &S,Vec<SstCarac> &matprops,bool &plasticite) const{ 
    S.matprop = matprops[S.typmat]; 
    //formulation isotrope 
    if (matprops[S.typmat].type.find("isotrope")<matprops[S.typmat].type.size() and matprops[S.typmat].comp.find("elastique")<matprops[S.typmat].comp.size()) {
        S.f = S.pb.formulation_elasticity_isotropy_stat_Qstat;
        S.mesh.type_formulation="isotrope"; 
    } else
    //formulation orthotrope 
    if (matprops[S.typmat].type.find("orthotrope")<matprops[S.typmat].type.size() and matprops[S.typmat].comp.find("elastique")<matprops[S.typmat].comp.size()) {
        S.f = S.pb.formulation_elasticity_orthotropy_stat_Qstat;
        S.mesh.type_formulation="orthotrope"; 
    } else
    //*//formulation plastique isotrope 
    if (matprops[S.typmat].type.find("isotrope")<matprops[S.typmat].type.size() and matprops[S.typmat].comp.find("plastique")<matprops[S.typmat].comp.size()) {
        plasticite = true;   /// Informe process qu'au moins un des materiaux est plastifiable
        S.plastique = true;
        S.f = S.pb.formulation_plasticity_isotropy_stat_Qstat;
        S.mesh.type_formulation="plastique"; 
    }//*/
    //formulation mesomodele 
    if ((matprops[S.typmat].type.find("orthotrope")<matprops[S.typmat].type.size()  and matprops[S.typmat].comp.find("plastique")<matprops[S.typmat].comp.size()) or matprops[S.typmat].comp.find("mesomodele")<matprops[S.typmat].comp.size()) {
        plasticite = true;   /// Informe process qu'au moins un des materiaux est plastifiable
        S.plastique = true;
        S.f = S.pb.formulation_mesomodele;  // TODO creer une formulation_plasticity_isotropy_stat_Qstat
        S.mesh.type_formulation="mesomodele"; 
    }
    if(not S.f){
        std::cerr << "Aucune formulation ne correspond aux type et comportement de materiau demandes!" << std::endl;
        assert(0);
    }
}
