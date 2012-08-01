#include "SstCarac_InterCarac.h"

ParameterGroup SstCarac::sst_materials_parameters;

/// Constructeur. Les symboles des parametres materiaux ne sont pas important
SstCarac::SstCarac():
density("rho",&sst_materials_parameters),
elastic_modulus("E",&sst_materials_parameters),
poisson_ratio("nu",&sst_materials_parameters),
elastic_modulus_1("E1",&sst_materials_parameters),
elastic_modulus_2("E2",&sst_materials_parameters),
elastic_modulus_3("E3",&sst_materials_parameters),
poisson_ratio_12("nu12",&sst_materials_parameters),
poisson_ratio_13("nu13",&sst_materials_parameters),
poisson_ratio_23("nu23",&sst_materials_parameters),
shear_modulus_12("G12",&sst_materials_parameters),
shear_modulus_13("G13",&sst_materials_parameters),
shear_modulus_23("G23",&sst_materials_parameters),
alpha("alpha",&sst_materials_parameters),
alpha_1("alpha1",&sst_materials_parameters),
alpha_2("alpha2",&sst_materials_parameters),
alpha_3("alpha3",&sst_materials_parameters),
plast_ecrouissage_init("R0",&sst_materials_parameters),
plast_ecrouissage_mult("k_p",&sst_materials_parameters),
plast_ecrouissage_expo("m_p",&sst_materials_parameters),
plast_cinematique_coef("c_p",&sst_materials_parameters),
Yo("Yo",&sst_materials_parameters),
Yc("Yc",&sst_materials_parameters),
Ycf("Ycf",&sst_materials_parameters),
viscosite("viscosite",&sst_materials_parameters)
{
    dt=1.0;
    resolution=1;
    /// Initialisation des grandeurs a -1 pour detecter les erreurs de chargement
    deltaT = -1;
    dmax = -1;
    b_c = -1;
    a = -1;
    tau_c = -1;
    couplage = -1;
    if(DIM == 2){
        v1[0] = 1; v2[0] = 0;
        v1[1] = 0; v2[1] = 1;
    }else{
        v1[0] = 1; v2[0] = 0;
        v1[1] = 0; v2[1] = 1;
        v1[2] = 0; v2[2] = 0;
    }
    ///Initialisation des coefficients de la contrainte equivalente a Von Mises standart
    /*
    for(int i = 0; i < DIM*(DIM+1)/2; i++)
        coeff_seq[i].set(0.0);
    if(DIM == 2){
        ;
    } else if(DIM == 3){
        coeff_seq[0][0] = 1;
        coeff_seq[1][1] = 1;
        coeff_seq[2][2] = 1;
        coeff_seq[0][1] = -1;
        coeff_seq[0][2] = -1;
        coeff_seq[1][2] = -1;
        coeff_seq[3][3] = 3;
        coeff_seq[4][4] = 3;
        coeff_seq[5][5] = 3;
    }
    //*/
}

void SstCarac::read_data_user(int index, Metil::DataUser& data_user){
    const DataUser::Json_materials &new_material = data_user.materials_vec[index];
    id = new_material.id_in_calcul;
    type_num = new_material.type_num;
    type = new_material.mtype;
    type_plast = new_material.type_plast;
    type_endo = new_material.type_endo;
    comp = new_material.comp;
    std::cout << "id : " << id << "    type : " << type << "    comp : " << comp << std::endl;
    if(data_user.dim == 2){/*   2D EN STAND-BY
        if (new_material.resolution =="CP"){
            resolution=1;
        } else if (new_material.resolution =="DP"){
            resolution=0;
        } else {
            std::cout << "type de resolution non implemente : choix contrainte_plane ou deformation_plane" << std::endl;
            assert(0);
        }*/
    }else{
        resolution=0;
    }
    
    density.setExpression(new_material.rho);
    
    if(type.find("isotrope")<type.size()) {/// comportement thermo-elastique isotrope
            elastic_modulus = new_material.E_1;         /// E
            poisson_ratio   = new_material.nu_12;       /// nu
            alpha           = new_material.alpha_1;     /// alpha
            deltaT          = 0;                        /// deltaT
    } else if (type.find("orthotrope")<type.size()) {/// comportement thermo-elastique orthotrope
            /// coefficients d'elasticite
            elastic_modulus_1 = new_material.E_1;      /// E1
            elastic_modulus_2 = new_material.E_2;      /// E2
            elastic_modulus_3 = new_material.E_3;      /// E3
            
            poisson_ratio_12 = new_material.nu_12;     /// nu12
            poisson_ratio_13 = new_material.nu_13;     /// nu13
            poisson_ratio_23 = new_material.nu_23;     /// nu23
            
            shear_modulus_12.setExpression(new_material.cis_1);    /// G12
            shear_modulus_13.setExpression(new_material.cis_2);    /// G13
            shear_modulus_23.setExpression(new_material.cis_3);    /// G23
            
            /// directions d'orthotropie TODO version parametrique
            v1[0]=new_material.dir_1_x;
            v2[0]=new_material.dir_2_x;
            v1[1]=new_material.dir_1_y;
            v2[1]=new_material.dir_2_y;
            if(DIM == 3){
                v1[2]=new_material.dir_1_z;
                v2[2]=new_material.dir_2_z;
            }
            v1=v1/norm_2(v1);
            v2=v2/norm_2(v2);
            
            /// coefficients thermiques
            alpha_1.setExpression(new_material.alpha_1);           ///alpha_1
            alpha_2.setExpression(new_material.alpha_2);           ///alpha_2
            alpha_3.setExpression(new_material.alpha_3);           ///alpha_3
            deltaT  = 0;                                           /// deltaT
    }
    
    if(comp.find("visqueux")<comp.size()){/// Comportement visqueux
            viscosite.setExpression(new_material.viscosite);   /// viscosite
    }
    if (comp.find("pl")<comp.size() or comp.find("mesomodele")<comp.size()) {
        /// parametres de plasticite
        plast_ecrouissage_init.setExpression(new_material.R0);                         /// R0
        plast_ecrouissage_mult.setExpression(new_material.k_p);                        /// k_p
        plast_ecrouissage_expo.setExpression(new_material.m_p);                        /// m_p
        plast_cinematique_coef.setExpression(new_material.coeff_plast_cinematique);    /// couplage
    }
    if (comp.find("en")<comp.size()){
        /// parametres d'endommagement
        Yo.setExpression(new_material.Yo); /// Yo (limite initiale d'endommagement)
        dmax = new_material.dmax;           /// dmax (valeur maximale de l'endommagement - pour eviter les divisions par 0)
        b_c  = new_material.b_c;            /// b_c (voir endommagement.cpp::calcul_endommagement)
    }
    if (comp.find("mesomodele")<comp.size()) {
        /// parametres du mesomodele
        Yo.setExpression(new_material.Yo);         /// Yo
        Yc.setExpression(new_material.Yc);         /// Yc
        Ycf.setExpression(new_material.Ycf);       /// Ycf
        dmax         = new_material.dmax;           /// dmax
        b_c          = new_material.b_c;            /// b_c
        effet_retard = new_material.effet_retard;   /// effet_retard
        a            = new_material.a;              /// a
        tau_c        = new_material.tau_c;          /// tau_c
    }//*/
}

void SstCarac::prepareParameters(){
    sst_materials_parameters.prepareParameters();
}

void SstCarac::updateParameters(){
    sst_materials_parameters.updateParameters();
}

void SstCarac::affiche(){
    std::cout << std::endl << std::endl;
    std::cout << "******************SSTCARAC_DEBUG*********************" << std::endl;
    std::cout << "id : " << id << std::endl;
    std::cout << "type : " << type << std::endl;
    std::cout << "comp : " << comp << std::endl;
    std::cout << "density : " << density << std::endl;
    if (type == "isotrope"){
        std::cout << "elastic_modulus : " << elastic_modulus << std::endl;
        std::cout << "poisson_ratio : " << poisson_ratio << std::endl;
        std::cout << "alpha : " << alpha << std::endl;
    } else if (type == "orthotrope"){
        std::cout << "v1 : " << v1[0] << "," << v1[1] << "," << v1[2] << std::endl;
        std::cout << "v2 : " << v2[0] << "," << v2[1] << "," << v2[2] << std::endl;
        std::cout << "elastic_modulus_1 : " << elastic_modulus_1 << std::endl;
        std::cout << "elastic_modulus_2 : " << elastic_modulus_2 << std::endl;
        std::cout << "elastic_modulus_3 : " << elastic_modulus_3 << std::endl;
        std::cout << "poisson_ratio_12 : " << poisson_ratio_12 << std::endl;
        std::cout << "poisson_ratio_13 : " << poisson_ratio_13 << std::endl;
        std::cout << "poisson_ratio_23 : " << poisson_ratio_23 << std::endl;
        std::cout << "shear_modulus_12 : " << shear_modulus_12 << std::endl;
        std::cout << "shear_modulus_13 : " << shear_modulus_13 << std::endl;
        std::cout << "shear_modulus_23 : " << shear_modulus_23 << std::endl;
        std::cout << "alpha_1 : " << alpha_1 << std::endl;
        std::cout << "alpha_2 : " << alpha_2 << std::endl;
        std::cout << "alpha_3 : " << alpha_3 << std::endl;
    }
    std::cout << "deltaT  : " << deltaT << std::endl;
    if (comp.find("pl")<comp.size() or comp == "mesomodele"){
        std::cout << "plast_ecrouissage_init : " << plast_ecrouissage_init << std::endl;
        std::cout << "plast_ecrouissage_mult : " << plast_ecrouissage_mult << std::endl;
        std::cout << "plast_ecrouissage_expo : " << plast_ecrouissage_expo << std::endl;
        std::cout << "plast_cinematique_coef : " << plast_cinematique_coef << std::endl;
    }
    if (comp == "mesomodele"){
        std::cout << "Yo           : " << Yo << std::endl;
        std::cout << "Yc           : " << Yc << std::endl;
        std::cout << "Ycf          : " << Ycf << std::endl;
        std::cout << "dmax         : " << dmax << std::endl;
        std::cout << "b_c          : " << b_c << std::endl;
        std::cout << "effet_retard : " << effet_retard << std::endl;
        std::cout << "a            : " << a << std::endl;
        std::cout << "tau_c        : " << tau_c << std::endl;
    }
    std::cout << "*****************************************************" << std::endl;
}


ParameterGroup InterCarac::inter_materials_parameters;

/// Constructeur. Les symboles des parametres materiaux ne sont pas important
InterCarac::InterCarac():
Ep_n("Ep_n",&inter_materials_parameters),
Ep_x("Ep_x",&inter_materials_parameters),
Ep_y("Ep_y",&inter_materials_parameters),
Ep_z("Ep_z",&inter_materials_parameters),
Preload_n("Preload_n",&inter_materials_parameters),
Preload_x("Preload_x",&inter_materials_parameters),
Preload_y("Preload_y",&inter_materials_parameters),
Preload_z("Preload_z",&inter_materials_parameters),
f("f",&inter_materials_parameters),
Fcr_n("Fcr_n",&inter_materials_parameters),
Fcr_t("Fcr_t",&inter_materials_parameters),
Kn("Kn",&inter_materials_parameters),
Kt("Kt",&inter_materials_parameters),
Knc("Knc",&inter_materials_parameters),
Rop("Rop",&inter_materials_parameters),
kp("kp",&inter_materials_parameters),
mp("mp",&inter_materials_parameters),
n("n",&inter_materials_parameters),
Yo("Yo",&inter_materials_parameters),
Yc("Yc",&inter_materials_parameters),
alpha("alpha",&inter_materials_parameters),
gamma("gamma",&inter_materials_parameters)
{
}

void InterCarac::allocate(){
}

void InterCarac::free(){
}


void InterCarac::read_data_user(int index,const DataUser &data_user){
    const DataUser::Json_links &link = data_user.links_vec[index];
    id = link.id_in_calcul;
    name = link.name;
    type_num = link.type_num;
    degradable = (type_num == 5);
    name = link.name;
    switch(type_num){
        case 0:
            type = "parfait";
            comp = "Parfait";
            break;
        case 1:
            type = "elastique";
            comp = "Elastique";
            break;
        case 2:
            type = "contact parfait";
            comp="Contact Parfait";
            break;
        case 3:
            type = "parfait cassable";
            comp = "Parfait Cassable";
            break;
        case 4:
            type = "elastique cassable";
            comp = "Elastique Cassable";
            break;
        case 5:
            type = "cohesive";
            comp = "Cohesive";
            break;
        default:
            std::cerr  << "comportement d'interface non reconnu : type_num = " << link.type_num << std::endl;
            assert(0);
    }
    /// Caracteristique du prechargement
    Ep_Type   = link.Ep_type;     /// Type de condition interne (epaisseur/precharge impose(e)/normal(e))
    Ep_n      = link.Ep_n;        /// jeux ou epaisseur normal
    Ep_x      = link.Ep_x;        /// jeux ou epaisseur selon x
    Ep_y      = link.Ep_y;        /// jeux ou epaisseur selon y
    Ep_z      = link.Ep_z;        /// jeux ou epaisseur selon z
    Preload_n = link.Preload_n;   /// precharge normale
    Preload_x = link.Preload_x;   /// precharge selon x
    Preload_y = link.Preload_y;   /// precharge selon y
    Preload_z = link.Preload_z;   /// precharge selon z
    
    /// Caracteristiques pour le contact avec frottement
    f         = link.f;           /// coefficient de frottement
    
    /// Caracteristiques pour le comportement cassable
    Fcr_n     = link.Fcr_n;       /// limite en rupture normale
    Fcr_t     = link.Fcr_t;       /// limite en rupture tangentielle
    
    /// Caracteristiques pour le comportement elastique
    Kn        = link.Kn;          /// Raideur normale (globale ou en traction) de l'interface
    Kt        = link.Kt;          /// Raideur tangentielle de l'interface
    Knc       = link.Knc;         /// Raideur normale en compression de l'interface
    
    /// Caracteristiques pour le comportement plastique
    Rop       = link.Rop;         /// Contrainte initiale de plastification
    kp        = link.kp;          /// Multiplicateur de la loi d'ecrouissage
    mp        = link.np;          /// Exposant de la loi d'ecrouissage
    
    /// Caracteristiques pour les interfaces cohesives (type mesomodele)
    Yo        = link.Yo;          /// Valeur minimale d'effort thermodynamique equivalent pour amorcer l'endommagement
    Yc        = link.Yc;          /// Valeur critique d'effort thermodynamique equivalent pour l'endommagement
    alpha     = link.alpha;       /// exposant dans le calcul de l'effort thermodynamique equivalent
    gamma     = link.gamma;       /// couplage normal/tangentiel dans le calcul de l'effort thermodynamique equivalent
    n         = link.n;           /// coefficient dans le calcul de l'endommagement associe a un effort thermodynamique
}

void InterCarac::prepareParameters(){
    inter_materials_parameters.prepareParameters();
}

void InterCarac::updateParameters(){
    inter_materials_parameters.updateParameters();
}

void InterCarac::affiche(){
    std::cout << "----------------------------------------------------------------------"<< std::endl;
    PRINT(id);
    PRINT(name);
    PRINT(type);
    PRINT(comp);
    PRINT(Ep_Type);
    switch(Ep_Type){
        case 0:
            PRINT(Ep_n.str_expr);
            break;
        case 1:
            PRINT(Ep_x.str_expr);
            PRINT(Ep_y.str_expr);
            #if DIM == 3
            PRINT(Ep_z.str_expr);
            #endif
            break;
        case 2:
            PRINT(Preload_n.str_expr);
            break;
        case 3:
            PRINT(Preload_x.str_expr);
            PRINT(Preload_y.str_expr);
            #if DIM == 3
            PRINT(Preload_z.str_expr);
            #endif
            break;
            break;
    }
    PRINT(f.str_expr);
    PRINT(Fcr_n.str_expr);
    PRINT(Fcr_t.str_expr);
    PRINT(Kn.str_expr);
    PRINT(Knc.str_expr);
    PRINT(Kt.str_expr);
    PRINT(Rop.str_expr);
    PRINT(kp.str_expr);
    PRINT(mp.str_expr);
    PRINT(Yo.str_expr);
    PRINT(Yc.str_expr);
    PRINT(alpha.str_expr);
    PRINT(gamma.str_expr);
    PRINT(n);
    std::cout << "----------------------------------------------------------------------"<< std::endl;
}
