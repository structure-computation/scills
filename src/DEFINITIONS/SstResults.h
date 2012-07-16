#include "../../LMT/include/containers/vec.h"
#include "../COMPUTE/DataUser.h"
#include "Parameter.h"
using namespace LMT;



struct SstState{
    Vec<TYPEREEL> q;                              /// vecteur des deplacements aux noeuds
    Vec<TYPEREEL> p;                              /// vecteur de la plasticite cumulee aux points de Gauss
    Vec<TYPEREEL> R_p;                            /// vecteur de l'ecrouissage aux points de Gauss
    Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > epsilon_p;  /// vecteur des deformations plastiques aux points de Gauss
    Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > X_p;        /// vecteur des centres du domaine elastique de chaque elements
    //Vec<TYPEREEL,2> erreur_plasticite;          /// vecteur des erreurs en plasticite cumulee calculees aux points de Gauss PRISE EN COMPTE ERREUR_PLASTICITE NON IMPLEMENTEE
    Vec<TYPEREEL> d1;                             /// vecteur des endommagements (modele standart ou micro-fissuration du mesomodele) aux points de Gauss
    Vec<TYPEREEL> d2;                             /// vecteur des endommagements (decohesion fibres/matrice du mesomodele) aux points de Gauss
    Vec<TYPEREEL> df;                             /// vecteur des endommagements (rupture fragile des fibres du mesomodele) aux points de Gauss
    Vec<Vec<TYPEREEL,3> > Yd;                     /// vecteur des forces thermodynamyques d'endommagement (type mesomodele) aux points de Gauss
    
    SstState(const Process& process, Sst &S){
        allocate(process,S);
    }
    
    void clear(){
        q.set(0.0);
        p.set(0.0);
        R_p.set(0.0);
        for( unsigned i = 0; i < nbelem; ++i ) {
            epsilon_p[i].set(0.0);
        }
        for( unsigned i = 0; i < nbelem; ++i ) {
            X_p[i].set(0.0);
        }
        d1.set(0.0);
        d2.set(0.0);
        df.set(0.0);
        for( unsigned i = 0; i < nbelem; ++i ){
            Yd[i].set(0.0);
        }
    }
    
    ~SstState(){
        q.free();
        p.free();
        R_p.free();
        epsilon_p.free();
        X_p.free();
        d1.free();
        d2.free();
        df.free();
        Yd.free();
    }
    
    void allocate(const Process &process,Sst &S){
        unsigned nbddl = S.mesh.node_list_size*DIM;
        unsigned nbelem = S.mesh.elem_list_size;
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
};

struct SstResults{
    int nb_result_stored;
    Vec<SstState *> states;
    SstResults(){}
    ~SstResults(){
        free();
    }
    void allocate(const Process& process,Sst& S){
        states.resize(nb_result_stored);
        for(int i = 0; i < nb_result_stored; i++){
            states[i] = new SstState(process,S);
        }
    }
    
    void free(){
        for(int i = 0; i < nb_result_stored; i++){
            states[i].free();
        }
        states.free();
    }
    
    void next_pt(){
        
    }
    
    SstState &operator[](int pt){
        return *(states[pt]);
    }
};

struct SstMaterial{
    /// Attributs communs
    int id;                         /// identite du materiaux dans data_user
    int type_num;                   /// numero d'identite du comportement materiaux : 0=isotrope elastique
    Sc2String type;                 /// type de formulation : isotrope, orthotrope, orthotrope endommageable
    Sc2String comp;                 /// type de comportement : elastique, endommageable, plastique, mesomodele...
    bool resolution;                /// type de resolution contrainte_plane (1) ou deformation_plane (0) : utilise en 2d
    Vec<TYPEREEL,DIM> v1,v2;        /// direction pour les materiaux orthotropes
    TYPEREEL density;               /// densite du materiaux
    Vec<Parameter,DIM> f_vol;        /// Champs de force volumique constant
    TYPEREEL dt;                    /// pas de temps lu uniquement pour la quasistatique (obtenu a partir de process.temps->dt)
    
    /// Comportement elastique
    /// Isotrope
    Parameter elastic_modulus;
    Parameter poisson_ratio;
    /// Orthotrope
    Parameter elastic_modulus_1;
    Parameter elastic_modulus_2;
    Parameter elastic_modulus_3;
    Parameter poisson_ratio_12;
    Parameter poisson_ratio_13;
    Parameter poisson_ratio_23;
    Parameter shear_modulus_12;
    Parameter shear_modulus_13;
    Parameter shear_modulus_23;
    
    /// Comportement thermique
    /// Isotrope
    Parameter alpha;
    /// Orthotrope
    Parameter alpha_1;
    Parameter alpha_2;
    Parameter alpha_3;
    //TYPEREEL deltaT;                    /// Variable pour le stockage de la variation de temperature PAS UTILISE POUR LE MOMENT
    
    /// Compoprtement plastique
    Sc2String type_plast;   /// Options decrivant l'evolution du domaine elastique (parfaite, isotrope, cinematique) et la contrainte equivalente (von_mises)
    /// Coefficients de la loi d'ecrouissage : R(p) = R0 + k_p * p ^ m_p
    Parameter plast_ecrouissage_init;    /// Limite d'elasticite initiale
    Parameter plast_ecrouissage_mult;    /// Coefficient multiplicateur de la loi d'ecrouissage
    Parameter plast_ecrouissage_expo;    /// Exposant de la loi d'ecrouissage
    Parameter plast_cinematique_coef;    /// Coefficient multiplicateur du modele cinematique : dX = C*depsilon_p
    Vec<Vec<TYPEREEL,DIM*(DIM+1)/2>,DIM*(DIM+1)/2> coeff_seq;     /// Coefficients multiplicateurs pour la fonction seuil (matrice de la forme bilineaire associee a la norme energetique)
    /// COEFF_SEQ N'EST PAS UTILISE POUR L'INSTANT
    
    /// Comportement endommageable de type mesomodele de composite
    Sc2String type_endo;    /// Options decrivant l'endommagement A SPECIFIER!!!
    Parameter Yo,Yc,Ycf;    /// Efforts limites de l'endommagement de la rupture
    Parameter dmax;         /// Maximum possible de l'endommagement avant rupture (pour eviter les divisions par 0)
    Parameter couplage;     /// Coefficient pour le calcul de la contrainte equivalente
    Parameter b_c;          /// Couplage entre micro-fissuration et decohesion fibres/matrice
    bool effet_retard;      /// Indique si un effet retard doit etre applique a l'endommagement
    Parameter a,tau_c;      /// Coefficients de l'effet_retard de l'endommagement
    
    /// Comportement visqueux
    Parameter viscosite;     ///< Viscosite du materiau
    
    SstMaterial():
    #if DIM == 2
    f_vol(Vec<Parameter,2>(Parameter("f_vol_x",Parameter::Material),Parameter("f_vol_y",Parameter::Material))),
    #else
    f_vol(Vec<Parameter,3>(Parameter("f_vol_x",Parameter::Material),Parameter("f_vol_y",Parameter::Material),Parameter("f_vol_z",Parameter::Material))),
    #endif
    elastic_modulus(Parameter("E",Parameter::Material)),
    poisson_ratio(Parameter("nu",Parameter::Material)),
    elastic_modulus_1(Parameter("E1",Parameter::Material)),
    elastic_modulus_2(Parameter("E2",Parameter::Material)),
    elastic_modulus_3(Parameter("E3",Parameter::Material)),
    poisson_ratio_12(Parameter("nu12",Parameter::Material)),
    poisson_ratio_13(Parameter("nu13",Parameter::Material)),
    poisson_ratio_23(Parameter("nu23",Parameter::Material)),
    shear_modulus_12(Parameter("G12",Parameter::Material)),
    shear_modulus_13(Parameter("G13",Parameter::Material)),
    shear_modulus_23(Parameter("G23",Parameter::Material)),
    alpha(Parameter("alpha",Parameter::Material)),
    alpha_1(Parameter("alpha1",Parameter::Material)),
    alpha_2(Parameter("alpha2",Parameter::Material)),
    alpha_3(Parameter("alpha3",Parameter::Material)),
    plast_ecrouissage_init(Parameter("R0",Parameter::Material)),
    plast_ecrouissage_mult(Parameter("k_p",Parameter::Material)),
    plast_ecrouissage_expo(Parameter("m_p",Parameter::Material)),
    plast_cinematique_coef(Parameter("c_p",Parameter::Material)),
    Yo(Parameter("Yo",Parameter::Material)),
    Yc(Parameter("Yc",Parameter::Material)),
    Ycf(Parameter("Yf"),Parameter::Material),
    dmax(Parameter("dmax"),Parameter::Material),
    couplage(Parameter("couplage"),Parameter::Material),
    b_c(Parameter("b_c"),Parameter::Material),
    a(Parameter("a"),Parameter::Material),
    tau_c(Parameter("tau_c"),Parameter::Material)
    {
    }
    
    void read_data_user(ScJsonReader::Json_materials &material,Process& process){
        id = material.id;
        type_num = material.type_num;
        type = material.mtype;
        type_plast = material.type_plast;
        type_endo = material.type_endo;
        comp = material.comp;
        if(data_user.dim == 2){
            if (material.resolution =="CP")
                resolution=true;
            else if (material.resolution =="DP")
                resolution=false;
            else {
                std::cout << "type de resolution non implemente : choix contrainte_plane ou deformation_plane" << std::endl;
                assert(false);
            }
        }else{
            resolution=false;
        }
        
        for(int i_fvol=0; i_fvol<data_user.behaviour_bc_volume.size(); i_fvol++){
            if(data_user.behaviour_bc_volume[i_fvol].select){
                for(int d=0; d<data_user.dim; d++){
                    vstr[d] += " + " + data_user.behaviour_bc_volume[i_fvol].step[0].CLv_step_prop[d] + " * " + data_user.behaviour_bc_volume[i_fvol].step[0].CLv_step_prop[6];
                }
            }
        }
        
        //f_vol_e=vstr;
        f_vol=data;
        
        std::cout << "id : " << id << "    type : " << type << "    comp : " << comp << std::endl;
        if(type.find("isotrope")<type.size()) {/// comportement thermo-elastique isotrope
            elastic_modulus = material.E;       /// E
            poisson_ratio   = material.nu;      /// nu
            alpha           = material.alpha;   /// alpha
            //deltaT          = 0;                     /// deltaT
        } else if (type.find("orthotrope")<type.size()) {/// comportement thermo-elastique orthotrope
            /// coefficients d'elasticite
            elastic_modulus_1 = material.E1;   /// E1
            elastic_modulus_2 = material.E2;   /// E2
            elastic_modulus_3 = material.E3;   /// E3
            
            poisson_ratio_12 = material.nu12;   /// nu12
            poisson_ratio_13 = material.nu13;   /// nu13
            poisson_ratio_23 = material.nu23;   /// nu23
            
            shear_modulus_12 = material.G12;   /// G12
            shear_modulus_13 = material.G13;   /// G13
            shear_modulus_23 = material.G23;   /// G23
            
            /// directions d'orthotropie
            v1[0] = material.dir_1_x;
            v1[1] = material.dir_1_y;
            v1[2] = material.dir_1_z;
            v2[0] = material.dir_2_x;
            v2[1] = material.dir_2_y;
            v2[2] = material.dir_2_z;
            
            /// coefficients thermiques
            alpha_1 = material.alpha_1;      ///alpha_1
            alpha_2 = material.alpha_2;      ///alpha_2
            alpha_3 = material.alpha_3;      ///alpha_3
            //deltaT  = 0;                         /// deltaT
        }
        
        if(comp.find("visqueux")<comp.size()){/// Comportement visqueux
            viscosite       = material.viscosite;   /// viscosite
        }
        if (comp.find("plastique")<comp.size() or comp.find("mesomodele")<comp.size()) {
            /// parametres de plasticite
            plast_ecrouissage_init = material.R0;   /// R0
            plast_ecrouissage_mult = material.k_p;  /// k_p
            plast_ecrouissage_expo = material.m_p;  /// m_p
            plast_cinematique_coef = material.c_p;  /// c_p
        }
        if (comp.find("endommageable")<comp.size()){
            /// parametres d'endommagement
            Yo           = material.Yo;     /// Yo (limite initiale d'endommagement)
            dmax         = material.dmax;   /// dmax (valeur maximale de l'endommagement - pour eviter les divisions par 0)
            b_c          = material.b_c;    /// b_c (voir endommagement.cpp::calcul_endommagement)
        }
        if (comp.find("mesomodele")<comp.size()) {
            /// parametres du mesomodele
            Yo           = material.Yo;             /// Yo
            Yc           = material.Yc;             /// Yc
            Ycf          = material.Ycf;            /// Ycf
            dmax         = material.dmax;           /// dmax
            b_c          = material.b_c;            /// b_c
            effet_retard = material.effet_retard;   /// effet_retard
            a            = material.a;              /// a
            tau_c        = material.tau_c;          /// tau_c
        }
    }
};