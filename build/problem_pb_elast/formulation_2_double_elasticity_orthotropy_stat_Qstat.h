
#include "formulation/formulation.h"
namespace LMT {
#ifndef ELASTICITY_ORTHOTROPY_STAT_QSTAT_2_NUM_0_CARAC_H
#define ELASTICITY_ORTHOTROPY_STAT_QSTAT_2_NUM_0_CARAC_H
#ifndef STRUCT_ELASTICITY_ORTHOTROPY_STAT_QSTAT
#define STRUCT_ELASTICITY_ORTHOTROPY_STAT_QSTAT
struct elasticity_orthotropy_stat_Qstat {
  static const char *name() { return "elasticity_orthotropy_stat_Qstat"; }
};
#endif // ELASTICITY_ORTHOTROPY_STAT_QSTAT

template<class P_T>
class CaracFormulation<elasticity_orthotropy_stat_Qstat,2,P_T>  {
public:
  typedef P_T T;
  static const char *name() { return "elasticity_orthotropy_stat_Qstat"; }
  static const bool matrix_will_be_definite_positive=true;
  static const bool has_nodal_matrix = false;
  static const bool has_IS_contact_matrix=false;
  static const bool need_skin_assembly=false;
  typedef Norm1_is_inf Name_convergence_criterium;
  static const unsigned nb_vectors = 4;
  static const unsigned nb_matrices = 4;
  static const unsigned auto_contact = false;
  static const bool friction_coeff_is_a_nodal_variable = 0;
  static const unsigned offset_of_pos_unknown=3;
  static const unsigned pos_is_an_unknown = false;
  static const unsigned nb_der_var = 0;
  template<class TF> static void add_to_der_vars( TF &f, const Vec<T> &v ) {
  }
  static bool is_unknown(const std::string &s) { return (s=="dep"); }
  static unsigned num_in_vec_unknown(const std::string &s) { if ( s=="dep" )return 0; return 0;  }
  template<unsigned num_mat,unsigned inner=0> struct NodalMatricesCarac {
      static const bool symm = 1;
      static const bool herm = false;
      static const bool diag = false;
  };
  template<unsigned num_mat,unsigned inner=0> struct GlobalMatricesCarac {
      static const bool symm = 1;
      static const bool herm = false;
      static const bool diag = false;
  };
  
  static const unsigned nb_nodal_unknowns = 2;
  template<class TE,class TTs,class Tvec>
  inline static void set_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice) {
    node.dep[1]=vecs[0][indice+1]; node.dep[0]=vecs[0][indice+0];
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
    node.dep[1]=vecs[0][indice+1]; node.dep[0]=vecs[0][indice+0];
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_nodal_initial_conditions(const TE &node,const TTs &f,Tvec &vecs,unsigned indice) {
    vecs[0][indice+0]=node.dep[0]; vecs[1][indice+0]=node.dep[0]; vecs[2][indice+0]=node.dep[0]; vecs[3][indice+0]=node.dep[0]; vecs[4][indice+0]=node.dep[0];
    vecs[0][indice+1]=node.dep[1]; vecs[1][indice+1]=node.dep[1]; vecs[2][indice+1]=node.dep[1]; vecs[3][indice+1]=node.dep[1]; vecs[4][indice+1]=node.dep[1];
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_nodal_error(const TE &node,const TTs &f,const Tvec &vecs,int indice) {
    T reg0=vecs[1][indice+1]-vecs[0][indice+1]; T reg1=vecs[1][indice+0]-vecs[0][indice+0]; reg0=abs(reg0); reg1=abs(reg1); return max(reg0,reg1);
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_nodal(const TE &node,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
    old_vec[indice+1]=vecs[1][indice+1]; old_vec[indice+0]=vecs[1][indice+0];
  }
  
  static const unsigned nb_global_unknowns = 0;
  template<class TE,class TTs,class Tvec>
  inline static void set_global_unknowns(TE &mesh,const TTs &f,const Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_global_unknowns(TE &mesh,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_global_initial_conditions(const TE &mesh,const TTs &f,Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_global_error(const TE &mesh,const TTs &f,const Tvec &vecs,int indice) {
    return 0;
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_global(const TE &mesh,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
  }
  
  static const unsigned nb_skin_nodal_unknowns = 0;
  template<class TE,class TTs,class Tvec>
  inline static void set_skin_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_skin_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_skin_nodal_initial_conditions(const TE &node,const TTs &f,Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_skin_nodal_error(const TE &node,const TTs &f,const Tvec &vecs,int indice) {
    return 0;
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_skin_nodal(const TE &node,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
  }
  
  static const unsigned nb_skin_global_unknowns = 0;
  template<class TE,class TTs,class Tvec>
  inline static void set_skin_global_unknowns(TE &mesh,const TTs &f,const Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_skin_global_unknowns(TE &mesh,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_skin_global_initial_conditions(const TE &mesh,const TTs &f,Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_skin_global_error(const TE &mesh,const TTs &f,const Tvec &vecs,int indice) {
    return 0;
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_skin_global(const TE &mesh,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
  }
};
#endif // ELASTICITY_ORTHOTROPY_STAT_QSTAT_2_NUM_0_CARAC_H
extern "C" void apply_on_elements_after_solve_0_elasticity_orthotropy_stat_Qstat_Triangle_0( double * );
extern "C" void apply_on_elements_after_solve_1_elasticity_orthotropy_stat_Qstat_Triangle_1( double * );
extern "C" void apply_on_elements_after_solve_2_elasticity_orthotropy_stat_Qstat_Triangle_2( double * );
extern "C" void apply_on_elements_after_solve_3_elasticity_orthotropy_stat_Qstat_Triangle_3( double * );
extern "C" void apply_on_elements_after_solve_4_elasticity_orthotropy_stat_Qstat_Triangle_4( double * );
extern "C" void apply_on_elements_after_solve_5_elasticity_orthotropy_stat_Qstat_Triangle_5( double * );
extern "C" void apply_on_elements_after_solve_6_elasticity_orthotropy_stat_Qstat_Triangle_6( double * );
extern "C" void apply_on_elements_after_solve_7_elasticity_orthotropy_stat_Qstat_Triangle_7( double * );
extern "C" void apply_on_elements_after_solve_8_elasticity_orthotropy_stat_Qstat_Triangle_8( double * );
extern "C" void apply_on_elements_after_solve_9_elasticity_orthotropy_stat_Qstat_Triangle_9( double * );
extern "C" void apply_on_elements_after_solve_10_elasticity_orthotropy_stat_Qstat_Triangle_10( double * );
extern "C" void apply_on_elements_after_solve_11_elasticity_orthotropy_stat_Qstat_Triangle_11( double * );
extern "C" void apply_on_elements_after_solve_12_elasticity_orthotropy_stat_Qstat_Triangle_12( double * );
extern "C" void apply_on_elements_after_solve_13_elasticity_orthotropy_stat_Qstat_Triangle_13( double * );
extern "C" void apply_on_elements_after_solve_14_elasticity_orthotropy_stat_Qstat_Triangle_14( double * );
class Triangle;
template<unsigned A,class B,class C> class Node;
template<class A,class B,class C,class D,unsigned E> class Element;

// Carac for ...
template<class P_T_pos,class P_ND,class TED,unsigned nim,class TM,class T>
class CaracFormulationForElement<elasticity_orthotropy_stat_Qstat,Element<Triangle,DefaultBehavior,Node<2,P_T_pos,P_ND>,TED,nim>,TM,T> {
public:
    template<unsigned num_mat,unsigned inner=0> struct ElemMatricesCarac {
        static const bool symm = true;
        static const bool herm = false;
        static const bool diag = false;
        static const bool linear = true;
    };
    static const unsigned order_integration = 0;
    static const bool has_elementary_matrix = true;
    static const bool has_skin_elementary_matrix = false;
    template<class TE,class TF, class TVEVE> static void after_solve(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; T reg2=2*(*f.m).shear_modulus_12; reg1=1.0/reg1; reg0=1.0/reg0;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=reg2*reg3; T reg5=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg6=1.0/(*f.m).elastic_modulus_3;
    T reg7=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg8=pow((*f.m).v1[0],2); T reg9=pow((*f.m).v1[1],2); T reg10=pow((*f.m).v2[0],2); T reg11=pow((*f.m).v2[1],2);
    T reg12=1.0/(*f.m).elastic_modulus_2; T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; reg9=reg8+reg9; reg8=pow((*f.m).v1[2],2); reg11=reg10+reg11;
    reg10=pow((*f.m).v2[2],2); T reg14=reg6*reg4; T reg15=reg7*reg4; T reg16=reg5*reg4; reg10=reg11+reg10;
    reg11=reg5*reg15; reg8=reg9+reg8; reg9=reg13*reg14; T reg17=reg5*reg16; T reg18=reg12*reg14;
    T reg19=reg13*reg16; T reg20=reg12*reg15; reg8=pow(reg8,0.5); reg11=reg9+reg11; reg17=reg18-reg17;
    reg10=pow(reg10,0.5); reg18=1.0/(*f.m).elastic_modulus_1; T reg21=(*f.m).v2[2]/reg10; T reg22=reg19+reg20; T reg23=reg18*reg17;
    T reg24=(*f.m).v1[1]/reg8; T reg25=reg13*reg11; T reg26=(*f.m).v2[1]/reg10; T reg27=(*f.m).v1[2]/reg8; T reg28=reg7*reg16;
    reg25=reg23-reg25; reg23=reg2*reg0; T reg29=reg5*reg3; T reg30=reg7*reg3; reg3=reg6*reg3;
    T reg31=reg7*reg22; reg14=reg18*reg14; T reg32=reg7*reg15; T reg33=reg27*reg26; T reg34=reg12*reg4;
    reg4=reg13*reg4; T reg35=reg24*reg21; reg10=(*f.m).v2[0]/reg10; reg8=(*f.m).v1[0]/reg8; reg28=reg9+reg28;
    reg32=reg14-reg32; reg16=reg18*reg16; reg9=reg4*reg7; reg14=reg23*reg5; reg15=reg13*reg15;
    T reg36=reg27*reg10; T reg37=reg35-reg33; T reg38=reg2*reg1; T reg39=reg23*reg7; T reg40=reg3*reg12;
    T reg41=reg8*reg21; T reg42=reg7*reg34; reg31=reg25-reg31; reg23=reg23*reg6; reg25=2*reg10;
    T reg43=reg30*reg5; T reg44=2*reg8; reg3=reg3*reg13; T reg45=reg29*reg5; reg32=reg32/reg31;
    T reg46=2*reg37; reg11=reg11/reg31; reg42=reg19+reg42; reg34=reg18*reg34; reg28=reg28/reg31;
    reg17=reg17/reg31; T reg47=reg25*reg26; T reg48=reg38*reg6; T reg49=reg39*reg5; T reg50=pow(reg26,2);
    T reg51=reg24*reg10; T reg52=reg8*reg26; T reg53=pow(reg10,2); T reg54=pow(reg8,2); T reg55=reg44*reg24;
    T reg56=reg36-reg41; T reg57=pow(reg24,2); reg15=reg16+reg15; reg4=reg4*reg13; reg43=reg3+reg43;
    reg29=reg29*reg13; reg3=reg38*reg5; T reg58=reg23*reg12; reg9=reg16+reg9; reg30=reg30*reg12;
    reg23=reg23*reg13; reg45=reg40-reg45; reg16=reg14*reg5; reg38=reg38*reg7; reg14=reg14*reg13;
    reg45=reg45*reg18; reg40=reg52-reg51; reg43=reg43*reg13; reg4=reg34-reg4; reg39=reg39*reg12;
    reg23=reg49+reg23; reg15=reg15/reg31; reg34=reg48*reg12; reg48=reg48*reg13; reg16=reg58-reg16;
    reg42=reg42/reg31; reg49=reg54*reg17; reg58=reg53*reg11; reg22=reg22/reg31; reg9=reg9/reg31;
    T reg59=reg29+reg30; T reg60=reg47*reg32; T reg61=reg55*reg28; T reg62=reg50*reg32; T reg63=reg57*reg28;
    T reg64=reg53*reg32; T reg65=reg54*reg28; T reg66=reg47*reg11; T reg67=reg55*reg17; T reg68=reg50*reg11;
    T reg69=reg57*reg17; T reg70=reg38*reg5; T reg71=reg3*reg5; T reg72=reg46*reg56; T reg73=pow(reg56,2);
    T reg74=pow(reg37,2); T reg75=pow(reg27,2); T reg76=pow(reg21,2); reg43=reg45-reg43; reg45=reg57*reg42;
    T reg77=reg50*reg9; reg62=reg63+reg62; reg68=reg69+reg68; reg63=reg73*reg22; reg4=reg4/reg31;
    reg38=reg38*reg12; reg69=reg74*reg15; T reg78=reg75*reg17; T reg79=reg76*reg11; reg64=reg65+reg64;
    reg3=reg3*reg13; reg65=reg55*reg42; reg66=reg67+reg66; reg67=reg72*reg22; T reg80=reg47*reg9;
    T reg81=reg14+reg39; T reg82=reg72*reg15; reg23=reg23*reg13; T reg83=reg76*reg32; reg71=reg34-reg71;
    reg16=reg16*reg18; reg34=reg75*reg28; T reg84=reg54*reg42; T reg85=reg53*reg9; reg60=reg61+reg60;
    reg61=reg73*reg15; T reg86=reg74*reg22; reg58=reg49+reg58; reg49=pow(reg40,2); reg59=reg59*reg7;
    reg70=reg48+reg70; reg48=reg72*reg4; T reg87=reg49*reg15; reg80=reg65+reg80; reg65=reg76*reg9;
    T reg88=reg75*reg42; T reg89=reg73*reg4; reg77=reg45+reg77; reg45=reg74*reg4; reg85=reg84+reg85;
    reg59=reg43-reg59; reg23=reg16-reg23; reg82=reg60+reg82; reg86=reg58+reg86; reg70=reg70*reg13;
    reg16=2*reg26; reg43=reg25*reg21; reg58=reg3+reg38; reg60=2*reg24; reg84=reg10*reg26;
    T reg90=reg8*reg24; reg63=reg68+reg63; reg79=reg78+reg79; reg68=reg49*reg22; reg78=reg44*reg27;
    reg71=reg71*reg18; reg67=reg66+reg67; reg81=reg81*reg7; reg69=reg64+reg69; reg61=reg62+reg61;
    reg83=reg34+reg83; reg34=reg53*reg61; reg62=reg54*reg63; reg70=reg71-reg70; reg81=reg23-reg81;
    reg58=reg58*reg7; reg59=reg59/reg31; reg23=reg61*reg84; reg64=reg63*reg90; reg66=reg54*reg67;
    reg71=reg53*reg82; reg68=reg79+reg68; reg79=reg78*reg17; T reg91=reg43*reg11; reg87=reg83+reg87;
    reg48=reg80+reg48; reg80=reg53*reg69; reg83=reg54*reg86; T reg92=reg57*reg86; T reg93=reg50*reg69;
    T reg94=reg57*reg63; T reg95=reg50*reg61; T reg96=reg57*reg67; T reg97=reg50*reg82; T reg98=reg56*reg37;
    reg52=reg51+reg52; reg51=reg24*reg26; T reg99=reg8*reg10; T reg100=2*reg56; reg46=reg46*reg40;
    T reg101=reg60*reg27; T reg102=reg16*reg21; T reg103=reg82*reg84; T reg104=reg67*reg90; reg45=reg85+reg45;
    reg89=reg77+reg89; reg77=reg8*reg56; reg85=reg24*reg37; reg65=reg88+reg65; reg88=reg49*reg4;
    T reg105=reg78*reg28; T reg106=reg43*reg32; T reg107=reg54*reg68; reg97=reg96+reg97; reg96=reg73*reg48;
    T reg108=2*reg27; T reg109=reg89*reg98; T reg110=reg27*reg21; reg32=reg102*reg32; reg58=reg70-reg58;
    reg70=reg74*reg89; T reg111=reg60*reg26; T reg112=reg44*reg10; reg100=reg100*reg40; reg34=reg62+reg34;
    reg28=reg101*reg28; reg23=reg64+reg23; reg62=reg46*reg15; reg106=reg105+reg106; reg64=reg74*reg48;
    reg105=reg90*reg86; T reg113=reg69*reg84; T reg114=reg78*reg42; T reg115=reg43*reg9; T reg116=reg99*reg59;
    T reg117=reg51*reg59; T reg118=reg52*reg59; T reg119=reg74*reg45; reg80=reg83+reg80; reg71=reg66+reg71;
    reg88=reg65+reg88; reg93=reg92+reg93; reg11=reg102*reg11; reg17=reg101*reg17; reg65=reg46*reg22;
    reg91=reg79+reg91; reg66=reg73*reg45; reg79=reg53*reg87; reg95=reg94+reg95; reg83=reg73*reg89;
    reg92=reg57*reg68; reg94=reg50*reg87; T reg120=reg48*reg98; T reg121=reg26*reg37; reg103=reg104+reg103;
    reg81=reg81/reg31; reg77=reg85+reg77; reg85=reg24*reg56; reg104=reg8*reg37; T reg122=reg10*reg56;
    T reg123=reg45*reg98; reg62=reg106+reg62; reg106=reg13*reg112; reg32=reg28+reg32; reg15=reg100*reg15;
    reg28=reg12*reg111; T reg124=reg116*reg111; T reg125=reg27*reg40; reg79=reg107+reg79; reg66=reg93+reg66;
    reg93=reg10*reg37; reg107=reg26*reg56; T reg126=reg74*reg88; reg83=reg95+reg83; reg95=reg117*reg111;
    T reg127=reg53*reg18; T reg128=reg50*reg13; reg122=reg121+reg122; reg94=reg92+reg94; reg92=reg18*reg112;
    reg121=reg13*reg111; T reg129=reg73*reg88; T reg130=reg50*reg12; reg96=reg97+reg96; reg113=reg105+reg113;
    reg70=reg34+reg70; reg34=reg117*reg112; reg97=reg118*reg111; reg105=reg53*reg13; T reg131=reg52*reg117;
    reg109=reg23+reg109; reg31=reg58/reg31; reg23=reg68*reg90; reg58=reg87*reg84; T reg132=reg81*reg77;
    T reg133=reg81*reg85; T reg134=reg81*reg104; T reg135=reg21*reg108; T reg136=reg116*reg112; reg60=reg60*reg56;
    reg119=reg80+reg119; reg44=reg44*reg37; reg120=reg103+reg120; reg80=reg52*reg118; reg103=reg110*reg59;
    T reg137=reg118*reg112; reg64=reg71+reg64; reg22=reg100*reg22; reg11=reg17+reg11; reg9=reg102*reg9;
    reg17=reg46*reg4; reg42=reg101*reg42; reg65=reg91+reg65; reg115=reg114+reg115; reg71=reg54*reg65;
    reg97=reg96+reg97; reg91=reg53*reg62; reg96=reg53*reg7; reg114=reg50*reg5; T reg138=reg132*reg44;
    T reg139=reg50*reg62; reg137=reg64+reg137; reg126=reg79+reg126; reg64=reg133*reg44; reg79=reg132*reg60;
    T reg140=reg103*reg112; T reg141=reg57*reg65; T reg142=reg31*reg122; T reg143=reg31*reg107; T reg144=reg31*reg93;
    T reg145=reg134*reg44; reg136=reg119+reg136; reg119=reg81*reg125; reg4=reg100*reg4; reg9=reg42+reg9;
    reg17=reg115+reg17; reg42=reg21*reg40; reg115=reg8*reg40; T reg146=reg27*reg37; reg15=reg32+reg15;
    reg124=reg66+reg124; reg32=reg134*reg60; reg95=reg83+reg95; reg66=reg133*reg60; reg129=reg94+reg129;
    reg83=reg103*reg111; reg34=reg70+reg34; reg70=reg52*reg116; reg94=reg53*(*f.m).alpha_2; reg36=reg41+reg36;
    reg41=reg5*reg135; reg106=reg28-reg106; reg105=reg130-reg105; reg28=reg76*reg5; reg130=(*f.m).alpha_1*reg54;
    reg22=reg11+reg22; reg11=(*f.m).alpha_1*reg57; reg123=reg113+reg123; reg113=reg50*(*f.m).alpha_2; T reg147=reg132*reg77;
    reg80=reg120+reg80; reg25=reg25*reg37; reg128=reg127-reg128; reg120=reg5*reg111; reg16=reg16*reg56;
    reg127=reg7*reg135; reg131=reg109+reg131; reg121=reg92-reg121; reg92=reg133*reg77; reg109=reg88*reg98;
    reg58=reg23+reg58; reg23=reg7*reg112; T reg148=reg76*reg7; reg94=reg130+reg94; reg130=reg143*reg16;
    T reg149=reg74*(*f.m).alpha_3; reg66=reg95+reg66; reg32=reg124+reg32; reg95=reg144*reg16; reg124=reg73*(*f.m).alpha_3;
    T reg150=reg57*reg13; reg113=reg11+reg113; reg145=reg136+reg145; reg11=reg144*reg25; reg136=reg31*reg42;
    reg92=reg131+reg92; reg131=reg143*reg122; reg109=reg58+reg109; reg58=reg52*reg103; T reg151=reg36*reg59;
    T reg152=reg54*reg18; reg4=reg9+reg4; reg147=reg80+reg147; reg9=reg142*reg122; reg80=reg65*reg90;
    T reg153=reg62*reg84; T reg154=reg134*reg77; T reg155=reg24*reg40; T reg156=reg27*reg56; reg115=reg146+reg115;
    reg70=reg123+reg70; reg123=reg143*reg25; reg28=reg105-reg28; reg79=reg97+reg79; reg97=reg142*reg16;
    reg41=reg106-reg41; reg105=reg76*(*f.m).alpha_2; reg139=reg141+reg139; reg120=reg23+reg120; reg140=reg126+reg140;
    reg23=reg119*reg44; reg106=reg73*reg17; reg126=reg57*reg22; reg141=reg50*reg15; reg35=reg33+reg35;
    reg33=reg6*reg135; reg146=reg53*reg15; T reg157=reg54*reg22; reg138=reg137+reg138; reg137=reg142*reg25;
    reg91=reg71+reg91; reg71=reg74*reg17; reg114=reg96+reg114; reg96=reg76*reg6; reg148=reg128-reg148;
    reg64=reg34+reg64; reg34=(*f.m).alpha_1*reg75; reg127=reg121-reg127; reg83=reg129+reg83; reg121=reg54*reg13;
    reg128=reg119*reg60; reg129=reg57*reg12; T reg158=reg10*reg40; T reg159=reg21*reg37; reg95=reg32+reg95;
    reg106=reg139+reg106; reg32=reg151*reg111; reg139=reg51*reg28; T reg160=reg99*reg148; reg141=reg126+reg141;
    reg126=reg136*reg16; reg128=reg83+reg128; reg83=reg54*reg7; T reg161=reg57*reg28; T reg162=reg73*reg4;
    T reg163=reg57*reg5; reg120=reg33-reg120; reg130=reg66+reg130; reg114=reg96-reg114; reg158=reg159+reg158;
    reg10=reg10*reg21; reg33=reg15*reg84; reg97=reg79+reg97; reg8=reg8*reg27; reg149=reg94+reg149;
    reg66=reg22*reg90; reg124=reg113+reg124; reg79=reg17*reg98; reg105=reg34+reg105; reg34=reg49*(*f.m).alpha_3;
    reg94=(*f.m).alpha_1*reg90; reg96=reg84*(*f.m).alpha_2; reg113=reg54*reg127; reg159=reg57*reg41; T reg164=reg53*reg148;
    T reg165=reg50*reg28; reg153=reg80+reg153; reg80=reg54*reg148; reg9=reg147+reg9; reg147=reg26*reg40;
    reg154=reg70+reg154; reg70=reg75*reg7; reg150=reg152-reg150; reg152=reg75*reg5; reg121=reg129-reg121;
    reg129=reg21*reg56; T reg166=reg74*reg4; reg146=reg157+reg146; reg157=reg151*reg112; reg71=reg91+reg71;
    reg137=reg138+reg137; reg91=reg136*reg25; reg23=reg140+reg23; reg123=reg64+reg123; reg11=reg145+reg11;
    reg64=reg53*reg127; reg138=reg50*reg41; reg155=reg156+reg155; reg140=reg119*reg77; reg58=reg109+reg58;
    reg109=reg51*reg41; reg145=reg144*reg122; reg131=reg92+reg131; reg92=reg99*reg127; reg59=reg35*reg59;
    reg156=reg81*reg115; reg32=reg106+reg32; reg91=reg23+reg91; reg111=reg59*reg111; reg23=reg11*reg149;
    reg106=reg123*reg124; reg162=reg141+reg162; reg126=reg128+reg126; reg138=reg64+reg138; reg64=reg76*reg120;
    reg128=reg95*reg149; reg141=reg130*reg124; reg139=reg160+reg139; reg160=reg110*reg114; reg112=reg59*reg112;
    reg166=reg146+reg166; reg34=reg105+reg34; reg96=reg94+reg96; reg94=reg98*(*f.m).alpha_3; reg105=(*f.m).alpha_1*reg8;
    reg146=reg10*(*f.m).alpha_2; reg147=reg129+reg147; reg79=reg153+reg79; reg129=reg52*reg151; reg33=reg66+reg33;
    reg26=reg26*reg21; reg27=reg24*reg27; reg37=reg40*reg37; reg24=reg137*reg131; reg66=reg75*reg114;
    reg161=reg80+reg161; reg70=reg150-reg70; reg145=reg154+reg145; reg80=reg131*reg97; reg152=reg121-reg152;
    reg121=reg52*reg2; reg84=reg84*reg2; reg150=reg123*reg9; reg153=reg9*reg130; reg154=reg136*reg122;
    T reg167=reg156*reg44; reg157=reg71+reg157; reg140=reg58+reg140; reg163=reg83+reg163; reg58=reg75*reg6;
    reg71=reg156*reg60; reg81=reg81*reg155; reg83=reg31*reg158; reg109=reg92+reg109; reg92=reg4*reg98;
    reg159=reg113+reg159; reg113=reg75*reg120; T reg168=reg76*reg114; reg165=reg164+reg165; reg164=reg110*reg120;
    T reg169=reg83*reg25; reg167=reg157+reg167; reg157=reg47*reg84; reg168=reg165+reg168; reg164=reg109+reg164;
    reg163=reg58-reg163; reg58=reg83*reg16; reg71=reg32+reg71; reg32=reg50*reg152; reg106=reg23+reg106;
    reg91=reg91*reg34; reg23=reg53*reg70; reg154=reg140+reg154; reg64=reg138+reg64; reg109=reg47*reg121;
    reg111=reg162+reg111; reg138=reg81*reg60; reg31=reg31*reg147; reg80=reg153-reg80; reg24=reg150-reg24;
    reg140=reg123*reg97; reg150=reg137*reg130; reg153=reg145*reg149; reg162=reg131*reg124; reg129=reg79+reg129;
    reg79=reg156*reg77; reg56=reg40*reg56; reg165=reg55*reg84; reg66=reg161+reg66; reg161=reg52*reg59;
    reg92=reg33+reg92; reg33=reg57*reg152; T reg170=reg54*reg70; T reg171=reg36*reg1; T reg172=reg10*reg1;
    reg141=reg128+reg141; reg126=reg126*reg34; reg113=reg159+reg113; reg128=reg55*reg121; reg112=reg166+reg112;
    reg94=reg96+reg94; reg96=reg81*reg44; reg159=reg26*(*f.m).alpha_2; reg146=reg105+reg146; reg105=(*f.m).alpha_1*reg27;
    reg166=reg37*(*f.m).alpha_3; T reg173=reg52*reg84; T reg174=reg52*reg121; reg160=reg139+reg160; reg139=reg81*reg77;
    T reg175=reg78*reg172; reg165=reg66+reg165; reg157=reg168+reg157; reg161=reg92+reg161; reg174=reg164+reg174;
    reg66=reg75*reg163; reg33=reg170+reg33; reg92=reg43*reg172; reg90=reg90*reg2; reg96=reg112+reg96;
    reg112=reg31*reg25; reg164=reg35*reg0; reg168=reg26*reg0; reg166=reg146+reg166; reg138=reg111+reg138;
    reg111=reg36*reg172; reg146=reg11*reg80; reg170=reg95*reg24; reg159=reg105+reg159; reg105=reg56*(*f.m).alpha_3;
    reg150=reg140-reg150; reg140=reg36*reg171; T reg176=reg31*reg16; reg162=reg153+reg162; reg79=reg129+reg79;
    reg129=reg83*reg122; reg154=reg154*reg34; reg153=reg99*reg70; T reg177=reg51*reg152; T reg178=reg43*reg171;
    reg109=reg64+reg109; reg128=reg113+reg128; reg64=reg78*reg171; reg113=reg137*reg94; reg91=reg106+reg91;
    reg32=reg23+reg32; reg23=reg76*reg163; reg58=reg71+reg58; reg126=reg141+reg126; reg71=reg97*reg94;
    reg169=reg167+reg169; reg173=reg160+reg173; reg169=reg169*reg166; reg106=reg110*reg163; reg177=reg153+reg177;
    reg141=reg47*reg90; reg153=reg31*reg122; reg160=reg9*reg94; reg154=reg162+reg154; reg162=reg35*reg164;
    reg167=reg11*reg9; reg139=reg161+reg139; reg23=reg32+reg23; reg32=reg145*reg97; reg170=reg146-reg170;
    reg113=reg91+reg113; reg91=reg9*reg95; reg146=reg145*reg150; reg176=reg138+reg176; reg129=reg79+reg129;
    reg79=reg8*reg1; reg58=reg58*reg166; reg138=reg101*reg168; reg140=reg174+reg140; reg105=reg159+reg105;
    reg159=elem.pos(1)[0]-elem.pos(0)[0]; reg161=reg102*reg168; reg92=reg157+reg92; reg64=reg128+reg64; reg128=reg101*reg164;
    reg111=reg173+reg111; reg157=reg102*reg164; reg173=elem.pos(2)[1]-elem.pos(0)[1]; reg112=reg96+reg112; reg96=elem.pos(1)[1]-elem.pos(0)[1];
    reg174=reg55*reg90; reg66=reg33+reg66; reg178=reg109+reg178; reg33=reg35*reg168; reg71=reg126+reg71;
    reg175=reg165+reg175; reg109=elem.pos(2)[0]-elem.pos(0)[0]; reg126=reg137*reg145; reg176=reg176*reg105; reg33=reg111+reg33;
    reg58=reg71+reg58; reg128=reg64+reg128; reg106=reg177+reg106; reg157=reg178+reg157; reg112=reg112*reg105;
    reg64=reg159*reg173; reg71=reg52*reg90; reg146=reg170+reg146; reg169=reg113+reg169; reg153=reg139+reg153;
    reg111=reg27*reg0; reg126=reg167-reg126; reg113=reg123*reg145; reg139=reg11*reg97; reg174=reg66+reg174;
    reg66=reg78*reg79; reg165=reg137*reg95; reg138=reg175+reg138; reg167=reg109*reg96; reg162=reg140+reg162;
    reg140=reg43*reg79; reg141=reg23+reg141; reg129=reg129*reg166; reg161=reg92+reg161; reg23=reg11*reg131;
    reg92=reg131*reg95; reg170=reg145*reg130; reg160=reg154+reg160; reg32=reg91-reg32; reg91=reg102*reg111;
    reg154=reg123*reg95; reg80=reg80/reg146; reg113=reg23-reg113; reg32=reg32/reg146; reg165=reg139-reg165;
    reg126=reg126/reg146; reg23=reg11*reg130; reg112=reg169+reg112; reg170=reg92-reg170; reg140=reg141+reg140;
    reg24=reg24/reg146; reg92=reg101*reg111; reg66=reg174+reg66; reg139=reg33*reg157; reg141=reg33*reg128;
    reg169=reg162*reg161; reg174=reg36*reg79; reg71=reg106+reg71; reg106=reg162*reg138; reg167=reg64-reg167;
    reg58=reg176+reg58; reg129=reg160+reg129; reg153=reg105*reg153; reg96=reg96/reg167; reg150=reg150/reg146;
    reg109=reg109/reg167; reg165=reg165/reg146; reg92=reg66+reg92; reg113=reg113/reg146; reg170=reg170/reg146;
    reg24=reg58*reg24; reg173=reg173/reg167; reg80=reg80*reg112; reg91=reg140+reg91; reg153=reg129+reg153;
    reg64=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; reg66=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; reg126=reg58*reg126; reg129=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg174=reg71+reg174;
    reg71=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg140=reg35*reg111; reg160=reg128*reg161; reg175=reg157*reg138; reg141=reg106-reg141;
    reg32=reg32*reg112; reg154=reg23-reg154; reg139=reg169-reg139; reg159=reg159/reg167; reg165=reg153*reg165;
    reg24=reg80-reg24; reg150=reg153*reg150; reg32=reg126-reg32; reg112=reg170*reg112; reg113=reg58*reg113;
    reg23=reg141*reg91; reg146=reg154/reg146; reg58=reg96*reg129; reg80=reg173*reg71; reg106=reg109*reg64;
    reg126=reg159*reg66; reg140=reg174+reg140; reg154=reg92*reg139; reg160=reg175-reg160; reg169=1-(*f.m).resolution;
    reg106=reg126-reg106; reg58=reg80-reg58; reg71=reg109*reg71; reg129=reg159*reg129; reg66=reg96*reg66;
    reg80=PNODE(2).dep[1]-PNODE(0).dep[1]; reg24=reg150+reg24; reg23=reg154-reg23; reg126=reg140*reg160; reg150=PNODE(1).dep[0]-PNODE(0).dep[0];
    reg165=reg32-reg165; reg32=PNODE(2).dep[0]-PNODE(0).dep[0]; reg154=reg140*reg157; reg170=reg33*reg91; reg174=reg162*reg92;
    reg175=reg140*reg161; reg176=PNODE(1).dep[1]-PNODE(0).dep[1]; reg64=reg173*reg64; reg177=reg140*reg138; reg178=reg162*reg91;
    reg146=reg153*reg146; reg113=reg112-reg113; reg112=reg33*reg92; reg153=reg140*reg128; T reg179=reg92*reg161;
    reg177=reg112-reg177; reg112=reg138*reg91; T reg180=reg150*reg109; T reg181=(*f.m).resolution*reg124; T reg182=(*f.m).resolution*reg149;
    reg58=reg106+reg58; reg106=reg157*reg92; T reg183=reg32*reg159; T reg184=reg80*reg96; reg24=reg169*reg24;
    reg153=reg174-reg153; reg174=reg128*reg91; reg154=reg178-reg154; reg175=reg170-reg175; reg170=reg173*reg176;
    reg113=reg146+reg113; reg126=reg23+reg126; reg66=reg64-reg66; elem.epsilon[0][0]=reg66; reg71=reg129-reg71;
    elem.epsilon[0][1]=reg71; reg165=reg169*reg165; reg149=(*f.m).deltaT*reg149; reg124=(*f.m).deltaT*reg124; reg58=0.5*reg58;
    elem.epsilon[0][2]=reg58; reg23=reg66-reg149; reg112=reg179-reg112; reg64=reg71-reg124; reg174=reg106-reg174;
    reg106=(*f.m).deltaT*reg94; reg80=reg80*reg159; reg180=reg183-reg180; reg177=reg177/reg126; reg176=reg109*reg176;
    reg32=reg32*reg96; reg184=reg170-reg184; reg94=(*f.m).resolution*reg94; reg154=reg154/reg126; reg175=reg175/reg126;
    reg141=reg141/reg126; reg153=reg153/reg126; reg150=reg150*reg173; reg139=reg139/reg126; reg165=reg181+reg165;
    reg182=reg24+reg182; reg113=reg169*reg113; reg24=(*f.m).resolution*reg153; reg129=(*f.m).resolution*reg154; reg146=(*f.m).resolution*reg139;
    reg170=(*f.m).resolution*reg141; reg112=reg112/reg126; reg174=reg174/reg126; reg131=reg131*reg169; reg126=reg160/reg126;
    reg145=reg145*reg169; reg160=(*f.m).resolution*reg175; reg178=(*f.m).resolution*reg177; reg179=reg23*reg154; reg181=reg64*reg153;
    reg183=reg64*reg141; T reg185=reg23*reg139; T reg186=reg58-reg106; reg113=reg94+reg113; reg165=(*f.m).deltaT*reg165;
    reg182=(*f.m).deltaT*reg182; reg32=reg150-reg32; reg176=reg80-reg176; reg11=reg11*reg169; reg184=reg180+reg184;
    reg123=reg123*reg169; reg95=reg169*reg95; reg130=reg169*reg130; reg184=0.5*reg184; reg137=reg137*reg169;
    reg80=reg186*reg126; reg94=(*f.m).resolution*reg174; reg183=reg185-reg183; reg113=(*f.m).deltaT*reg113; reg150=(*f.m).resolution*reg112;
    reg180=reg32-reg182; reg129=reg95-reg129; reg24=reg130+reg24; reg95=reg176-reg165; reg160=reg145+reg160;
    reg170=reg123-reg170; reg146=reg11+reg146; reg11=(*f.m).resolution*reg126; reg9=reg9*reg169; reg178=reg131-reg178;
    reg123=reg186*reg174; reg179=reg181-reg179; reg97=reg169*reg97; reg130=reg178*reg95; reg131=reg24*reg95;
    reg150=reg9+reg150; reg183=reg80+reg183; reg11=reg137+reg11; reg94=reg97-reg94; reg9=reg184-reg113;
    reg123=reg179-reg123; reg80=reg160*reg180; reg97=reg146*reg180; reg137=reg170*reg95; reg145=reg129*reg180;
    reg169=reg11*reg9; reg179=reg94*reg9; reg130=reg80+reg130; reg131=reg145+reg131; reg137=reg97+reg137;
    reg80=reg183+reg123; reg97=reg150*reg9; reg23=reg23*reg175; reg64=reg64*reg177; reg80=reg80/3;
    reg169=reg137+reg169; reg179=reg131+reg179; reg130=reg97+reg130; reg165=reg71-reg165; reg180=reg180*reg169;
    reg186=reg186*reg112; reg123=reg123-reg80; reg95=reg95*reg179; reg183=reg183-reg80; reg182=reg66-reg182;
    reg64=reg23-reg64; reg130=2*reg130; reg113=reg58-reg113; reg183=pow(reg183,2); reg129=reg129*reg182;
    reg123=pow(reg123,2); reg9=reg130*reg9; reg24=reg165*reg24; reg180=reg95+reg180; reg170=reg165*reg170;
    reg146=reg146*reg182; reg64=reg186+reg64; reg9=reg180+reg9; reg80=pow(reg80,2); reg178=reg165*reg178;
    reg182=reg160*reg182; reg123=reg183+reg123; reg94=reg113*reg94; reg146=reg170+reg146; reg24=reg129+reg24;
    reg23=2*reg64; reg11=reg113*reg11; reg9=reg167*reg9; reg80=reg123+reg80; reg178=reg182+reg178;
    reg150=reg113*reg150; reg94=reg24+reg94; elem.sigma[0][1]=reg94; reg11=reg146+reg11; elem.sigma[0][0]=reg11;
    reg23=reg64*reg23; reg24=0.33333333333333331483*reg9; reg9=0.16666666666666665741*reg9; reg23=reg80+reg23; reg178=reg150+reg178;
    elem.sigma[0][2]=reg178; reg58=reg54*reg11; reg64=reg57*reg94; reg66=reg51*reg94; reg71=reg99*reg11;
    reg11=reg53*reg11; reg94=reg50*reg94; reg80=reg52*reg178; reg66=reg71+reg66; reg71=reg47*reg178;
    reg94=reg11+reg94; reg178=reg55*reg178; reg64=reg58+reg64; reg23=1.5*reg23; reg24=reg9+reg24;
    elem.sigma_von_mises=pow(reg23,0.5); elem.sigma_local[0][2]=reg66+reg80; elem.sigma_local[0][1]=reg94+reg71; elem.sigma_local[0][0]=reg64+reg178; elem.ener=reg24/2;
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_2(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_3(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_4(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_5(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_6(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_7(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_8(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_9(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_10(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_11(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_12(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_13(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_14(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_15(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
  
  static const unsigned nb_elementary_unknowns = 0;
  template<class TE,class TTs,class Tvec>
  inline static void set_elementary_unknowns(TE &elem,const TTs &f,const Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_elementary_unknowns(TE &elem,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_elementary_initial_conditions(const TE &elem,const TTs &f,Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_elementary_error(const TE &elem,const TTs &f,const Tvec &vecs,int indice) {
    return 0;
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_elementary(const TE &elem,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
  }
  
  static const unsigned nb_skin_elementary_unknowns = 0;
  template<class TE,class TTs,class Tvec>
  inline static void set_skin_elementary_unknowns(TE &elem,const TTs &f,const Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_skin_elementary_unknowns(TE &elem,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_skin_elementary_initial_conditions(const TE &elem,const TTs &f,Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_skin_elementary_error(const TE &elem,const TTs &f,const Tvec &vecs,int indice) {
    return 0;
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_skin_elementary(const TE &elem,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
  }
};

// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; reg1=1.0/reg1; T reg2=2*(*f.m).shear_modulus_12;
    T reg3=reg1*reg0; reg2=1.0/reg2; T reg4=reg2*reg3; T reg5=pow((*f.m).v1[0],2); T reg6=pow((*f.m).v2[1],2);
    T reg7=pow((*f.m).v2[0],2); T reg8=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg11=pow((*f.m).v1[1],2);
    T reg12=reg9*reg4; T reg13=reg8*reg4; T reg14=reg10*reg4; reg11=reg5+reg11; reg5=pow((*f.m).v1[2],2);
    T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; reg6=reg7+reg6; reg7=pow((*f.m).v2[2],2); T reg17=reg8*reg14;
    T reg18=reg16*reg12; T reg19=reg8*reg13; reg7=reg6+reg7; reg6=reg15*reg12; reg5=reg11+reg5;
    reg11=reg15*reg14; reg7=pow(reg7,0.5); reg5=pow(reg5,0.5); T reg20=reg16*reg13; reg17=reg18+reg17;
    T reg21=1.0/(*f.m).elastic_modulus_1; reg19=reg6-reg19; reg6=reg16*reg17; T reg22=reg20+reg11; T reg23=(*f.m).v1[1]/reg5;
    T reg24=(*f.m).v1[2]/reg5; T reg25=reg21*reg19; T reg26=(*f.m).v2[2]/reg7; T reg27=(*f.m).v2[1]/reg7; reg12=reg21*reg12;
    T reg28=reg2*reg0; T reg29=reg23*reg26; T reg30=reg24*reg27; T reg31=reg16*reg4; T reg32=reg10*reg14;
    reg5=(*f.m).v1[0]/reg5; T reg33=reg9*reg3; reg6=reg25-reg6; reg25=reg10*reg22; reg4=reg15*reg4;
    T reg34=reg10*reg13; reg7=(*f.m).v2[0]/reg7; T reg35=reg8*reg3; reg3=reg10*reg3; T reg36=reg28*reg10;
    T reg37=2*reg7; T reg38=reg33*reg15; reg32=reg12-reg32; reg12=reg28*reg9; reg13=reg21*reg13;
    T reg39=reg31*reg10; T reg40=reg29-reg30; T reg41=reg5*reg26; T reg42=reg3*reg8; T reg43=reg24*reg7;
    T reg44=reg35*reg8; reg33=reg33*reg16; reg28=reg28*reg8; reg14=reg16*reg14; T reg45=reg2*reg1;
    T reg46=2*reg5; reg25=reg6-reg25; reg34=reg18+reg34; reg6=reg10*reg4; reg18=pow(reg23,2);
    T reg47=pow(reg5,2); T reg48=reg43-reg41; T reg49=reg5*reg27; T reg50=reg23*reg7; T reg51=reg45*reg9;
    T reg52=pow(reg7,2); T reg53=pow(reg27,2); reg14=reg13+reg14; T reg54=reg37*reg27; reg31=reg31*reg16;
    T reg55=reg36*reg8; reg34=reg34/reg25; reg35=reg35*reg16; T reg56=2*reg40; T reg57=reg45*reg10;
    reg45=reg45*reg8; reg19=reg19/reg25; reg3=reg3*reg15; reg4=reg21*reg4; reg32=reg32/reg25;
    T reg58=reg46*reg23; reg17=reg17/reg25; reg6=reg20+reg6; T reg59=reg28*reg8; T reg60=reg12*reg16;
    reg12=reg12*reg15; reg42=reg33+reg42; reg39=reg13+reg39; reg44=reg38-reg44; reg13=reg53*reg32;
    reg33=reg52*reg17; reg31=reg4-reg31; reg4=reg47*reg19; reg38=reg58*reg34; T reg61=reg54*reg32;
    reg22=reg22/reg25; reg60=reg55+reg60; reg55=reg18*reg34; reg39=reg39/reg25; T reg62=reg51*reg15;
    reg51=reg51*reg16; T reg63=reg45*reg8; T reg64=reg52*reg32; T reg65=reg47*reg34; reg14=reg14/reg25;
    T reg66=reg57*reg8; T reg67=reg18*reg19; T reg68=reg54*reg17; T reg69=reg53*reg17; T reg70=reg58*reg19;
    reg59=reg12-reg59; reg12=reg35+reg3; reg42=reg42*reg16; reg6=reg6/reg25; reg44=reg44*reg21;
    T reg71=pow(reg24,2); T reg72=reg49-reg50; T reg73=pow(reg26,2); T reg74=pow(reg40,2); T reg75=pow(reg48,2);
    T reg76=reg56*reg48; reg36=reg36*reg15; reg28=reg28*reg16; T reg77=reg71*reg34; T reg78=reg73*reg32;
    T reg79=reg28+reg36; T reg80=reg75*reg14; reg13=reg55+reg13; reg55=reg58*reg6; T reg81=reg18*reg6;
    T reg82=reg74*reg14; reg64=reg65+reg64; reg65=reg53*reg39; reg63=reg62-reg63; reg66=reg51+reg66;
    reg51=reg76*reg22; reg68=reg70+reg68; reg69=reg67+reg69; reg62=reg75*reg22; reg67=reg73*reg17;
    reg70=reg71*reg19; T reg83=pow(reg72,2); reg31=reg31/reg25; reg45=reg45*reg16; reg57=reg57*reg15;
    reg33=reg4+reg33; reg42=reg44-reg42; reg12=reg12*reg10; reg4=reg74*reg22; reg59=reg59*reg21;
    reg44=reg76*reg14; reg60=reg60*reg16; T reg84=reg54*reg39; reg61=reg38+reg61; reg38=reg52*reg39;
    T reg85=reg47*reg6; T reg86=reg73*reg39; reg66=reg66*reg16; reg84=reg55+reg84; reg55=reg45+reg57;
    reg76=reg76*reg31; reg4=reg33+reg4; reg33=reg75*reg31; reg79=reg79*reg10; reg65=reg81+reg65;
    reg63=reg63*reg21; reg81=reg7*reg27; T reg87=reg5*reg23; T reg88=reg71*reg6; reg44=reg61+reg44;
    reg60=reg59-reg60; reg12=reg42-reg12; reg38=reg85+reg38; reg42=reg74*reg31; reg59=reg83*reg14;
    reg78=reg77+reg78; reg61=reg37*reg26; reg80=reg13+reg80; reg82=reg64+reg82; reg13=2*reg23;
    reg64=2*reg27; reg51=reg68+reg51; reg68=reg46*reg24; reg77=reg83*reg22; reg67=reg70+reg67;
    reg62=reg69+reg62; reg69=reg18*reg62; reg70=reg53*reg80; reg85=reg53*reg82; T reg89=reg18*reg51;
    T reg90=reg18*reg4; T reg91=reg53*reg44; T reg92=reg48*reg40; reg12=reg12/reg25; T reg93=reg44*reg81;
    reg76=reg84+reg76; reg84=reg68*reg34; T reg94=reg61*reg32; reg59=reg78+reg59; reg79=reg60-reg79;
    reg60=reg83*reg31; reg86=reg88+reg86; reg42=reg38+reg42; reg66=reg63-reg66; reg33=reg65+reg33;
    reg38=reg61*reg17; reg63=reg68*reg19; reg55=reg55*reg10; reg77=reg67+reg77; reg65=reg80*reg81;
    reg67=reg62*reg87; reg44=reg52*reg44; reg78=reg47*reg51; reg80=reg52*reg80; reg62=reg47*reg62;
    reg88=reg52*reg82; T reg95=reg47*reg4; reg49=reg50+reg49; reg50=reg23*reg27; T reg96=reg5*reg7;
    T reg97=2*reg48; reg56=reg56*reg72; T reg98=2*reg24; T reg99=reg13*reg24; T reg100=reg64*reg26;
    reg51=reg51*reg87; T reg101=reg5*reg48; T reg102=reg23*reg40; T reg103=reg46*reg7; T reg104=reg13*reg27;
    T reg105=reg47*reg77; T reg106=reg15*reg104; T reg107=reg16*reg103; T reg108=reg53*reg15; T reg109=reg52*reg16;
    T reg110=reg16*reg104; T reg111=reg21*reg103; T reg112=reg53*reg16; T reg113=reg52*reg21; T reg114=reg52*reg59;
    reg101=reg102+reg101; reg102=reg27*reg40; reg44=reg78+reg44; reg78=reg74*reg76; reg17=reg100*reg17;
    reg19=reg99*reg19; T reg115=reg56*reg22; reg38=reg63+reg38; reg55=reg66-reg55; reg4=reg87*reg4;
    reg82=reg82*reg81; reg63=reg7*reg48; reg65=reg67+reg65; reg66=reg33*reg92; reg79=reg79/reg25;
    reg67=reg68*reg6; T reg116=reg61*reg39; T reg117=reg18*reg77; T reg118=reg53*reg59; reg32=reg100*reg32;
    reg34=reg99*reg34; T reg119=reg56*reg14; reg94=reg84+reg94; reg91=reg89+reg91; reg84=reg75*reg76;
    reg89=reg96*reg12; reg98=reg26*reg98; T reg120=reg50*reg12; T reg121=reg49*reg12; reg97=reg97*reg72;
    T reg122=reg75*reg42; reg85=reg90+reg85; reg90=reg24*reg26; T reg123=reg75*reg33; T reg124=reg5*reg40;
    T reg125=reg23*reg48; reg88=reg95+reg88; reg95=reg74*reg42; reg60=reg86+reg60; reg70=reg69+reg70;
    reg80=reg62+reg80; reg33=reg74*reg33; reg93=reg51+reg93; reg76=reg76*reg92; reg123=reg70+reg123;
    reg51=reg73*reg10; reg112=reg113-reg112; reg62=reg120*reg104; reg69=reg52*reg10; reg70=reg10*reg98;
    reg110=reg111-reg110; reg59=reg59*reg81; reg77=reg77*reg87; reg82=reg4+reg82; reg42=reg42*reg92;
    reg25=reg55/reg25; reg118=reg117+reg118; reg4=reg73*reg8; reg55=reg90*reg12; reg39=reg100*reg39;
    reg86=reg89*reg103; reg95=reg88+reg95; reg78=reg44+reg78; reg44=reg121*reg103; reg6=reg99*reg6;
    reg22=reg97*reg22; reg17=reg19+reg17; reg124=reg79*reg124; reg125=reg79*reg125; reg115=reg38+reg115;
    reg19=reg74*reg60; reg114=reg105+reg114; reg38=reg79*reg101; reg56=reg56*reg31; reg116=reg67+reg116;
    reg32=reg34+reg32; reg107=reg106-reg107; reg34=reg10*reg103; reg119=reg94+reg119; reg67=reg8*reg104;
    reg88=reg24*reg72; reg94=reg75*reg60; reg105=reg49*reg121; reg76=reg93+reg76; reg93=reg120*reg103;
    reg33=reg80+reg33; reg63=reg102+reg63; reg121=reg121*reg104; reg84=reg91+reg84; reg13=reg13*reg48;
    reg80=reg7*reg40; reg91=reg27*reg48; reg109=reg108-reg109; reg66=reg65+reg66; reg120=reg49*reg120;
    reg14=reg97*reg14; reg46=reg46*reg40; reg122=reg85+reg122; reg65=reg89*reg104; reg85=reg8*reg98;
    reg102=reg53*reg8; reg14=reg32+reg14; reg121=reg84+reg121; reg64=reg64*reg48; reg32=reg124*reg46;
    reg84=reg55*reg103; reg19=reg114+reg19; reg39=reg6+reg39; reg15=reg18*reg15; reg93=reg33+reg93;
    reg6=reg125*reg46; reg33=reg55*reg104; reg94=reg118+reg94; reg56=reg116+reg56; reg106=reg18*reg16;
    reg4=reg109-reg4; reg37=reg37*reg40; reg108=reg125*reg13; reg16=reg47*reg16; reg70=reg110-reg70;
    reg85=reg107-reg85; reg51=reg112-reg51; reg62=reg123+reg62; reg107=reg18*reg115; reg109=reg53*reg119;
    reg59=reg77+reg59; reg77=reg52*(*f.m).alpha_2; reg42=reg82+reg42; reg89=reg49*reg89; reg43=reg41+reg43;
    reg41=reg124*reg13; reg65=reg122+reg65; reg120=reg66+reg120; reg66=reg73*reg9; reg102=reg69+reg102;
    reg69=reg38*reg101; reg105=reg76+reg105; reg76=reg24*reg40; reg82=reg5*reg72; reg98=reg9*reg98;
    reg67=reg34+reg67; reg125=reg125*reg101; reg34=reg26*reg72; reg31=reg97*reg31; reg86=reg95+reg86;
    reg95=(*f.m).alpha_1*reg47; reg44=reg78+reg44; reg78=reg38*reg46; reg22=reg17+reg22; reg88=reg79*reg88;
    reg60=reg60*reg92; reg80=reg25*reg80; reg91=reg25*reg91; reg17=reg47*reg115; reg97=reg52*reg119;
    reg110=(*f.m).alpha_1*reg18; reg111=reg53*(*f.m).alpha_2; reg21=reg47*reg21; reg38=reg38*reg13; reg112=reg25*reg63;
    reg113=reg18*reg8; reg114=reg18*reg85; reg116=reg47*reg70; reg117=reg47*reg10; reg118=reg53*reg14;
    reg122=reg112*reg64; reg38=reg121+reg38; reg109=reg107+reg109; reg107=reg75*reg56; reg121=reg53*reg85;
    reg123=reg18*reg22; T reg126=reg52*reg70; T reg127=reg23*reg72; T reg128=reg24*reg48; reg82=reg76+reg82;
    reg69=reg105+reg69; reg76=reg26*reg40; reg105=reg7*reg72; T reg129=reg112*reg63; reg115=reg115*reg87;
    reg41=reg65+reg41; reg77=reg95+reg77; reg65=reg74*(*f.m).alpha_3; reg111=reg110+reg111; reg95=reg75*(*f.m).alpha_3;
    reg110=(*f.m).alpha_1*reg71; T reg130=reg73*(*f.m).alpha_2; T reg131=reg80*reg64; reg34=reg25*reg34; T reg132=reg43*reg12;
    reg119=reg119*reg81; reg31=reg39+reg31; reg108=reg62+reg108; reg39=reg91*reg64; reg33=reg94+reg33;
    reg62=reg88*reg13; reg89=reg42+reg89; reg42=reg52*reg14; reg94=reg47*reg22; reg29=reg30+reg29;
    reg30=reg74*reg56; reg97=reg17+reg97; reg112=reg112*reg37; reg78=reg44+reg78; reg60=reg59+reg60;
    reg55=reg49*reg55; reg17=reg96*reg51; reg44=reg50*reg4; reg32=reg86+reg32; reg59=reg80*reg37;
    reg86=reg88*reg46; reg84=reg19+reg84; reg10=reg71*reg10; reg106=reg21-reg106; reg8=reg71*reg8;
    reg6=reg93+reg6; reg19=reg91*reg37; reg16=reg15-reg16; reg85=reg50*reg85; reg70=reg96*reg70;
    reg15=reg47*reg51; reg124=reg124*reg101; reg125=reg120+reg125; reg102=reg66-reg102; reg21=reg18*reg4;
    reg67=reg98-reg67; reg4=reg53*reg4; reg91=reg91*reg63; reg51=reg52*reg51; reg66=reg26*reg48;
    reg56=reg56*reg92; reg93=reg34*reg37; reg86=reg84+reg86; reg105=reg76+reg105; reg131=reg41+reg131;
    reg10=reg106-reg10; reg41=reg90*reg67; reg85=reg70+reg85; reg4=reg51+reg4; reg51=reg73*reg102;
    reg22=reg22*reg87; reg127=reg128+reg127; reg70=reg81*reg2; reg39=reg108+reg39; reg14=reg14*reg81;
    reg8=reg16-reg8; reg16=reg49*reg2; reg74=reg74*reg31; reg42=reg94+reg42; reg65=reg77+reg65;
    reg124=reg89+reg124; reg95=reg111+reg95; reg76=reg132*reg103; reg130=reg110+reg130; reg83=reg83*(*f.m).alpha_3;
    reg77=(*f.m).alpha_1*reg87; reg81=reg81*(*f.m).alpha_2; reg30=reg97+reg30; reg75=reg75*reg31; reg80=reg80*reg63;
    reg91=reg125+reg91; reg112=reg78+reg112; reg82=reg79*reg82; reg129=reg69+reg129; reg55=reg60+reg55;
    reg60=reg27*reg72; reg5=reg5*reg24; reg12=reg29*reg12; reg119=reg115+reg119; reg7=reg7*reg26;
    reg107=reg109+reg107; reg59=reg32+reg59; reg32=reg34*reg64; reg69=reg90*reg102; reg62=reg33+reg62;
    reg102=reg71*reg102; reg33=reg73*reg67; reg121=reg126+reg121; reg78=reg132*reg104; reg44=reg17+reg44;
    reg114=reg116+reg114; reg67=reg71*reg67; reg122=reg38+reg122; reg118=reg123+reg118; reg9=reg71*reg9;
    reg21=reg15+reg21; reg113=reg117+reg113; reg88=reg88*reg101; reg19=reg6+reg19; reg6=reg52*reg10;
    reg81=reg77+reg81; reg105=reg25*reg105; reg15=reg53*reg8; reg17=reg7*(*f.m).alpha_2; reg38=reg92*(*f.m).alpha_3;
    reg77=(*f.m).alpha_1*reg5; reg24=reg23*reg24; reg83=reg130+reg83; reg26=reg27*reg26; reg23=reg43*reg1;
    reg27=reg131*reg65; reg84=reg39*reg95; reg89=reg19*reg95; reg80=reg124+reg80; reg94=reg82*reg13;
    reg75=reg118+reg75; reg104=reg12*reg104; reg103=reg12*reg103; reg78=reg107+reg78; reg33=reg121+reg33;
    reg97=reg54*reg16; reg74=reg42+reg74; reg60=reg66+reg60; reg42=reg82*reg46; reg66=reg58*reg70;
    reg76=reg30+reg76; reg40=reg72*reg40; reg30=reg47*reg10; reg98=reg18*reg8; reg14=reg22+reg14;
    reg34=reg34*reg63; reg88=reg55+reg88; reg22=reg129*reg39; reg41=reg85+reg41; reg92=reg31*reg92;
    reg31=reg54*reg70; reg51=reg4+reg51; reg4=reg59*reg65; reg55=reg49*reg16; reg132=reg49*reg132;
    reg56=reg119+reg56; reg32=reg62+reg32; reg62=reg19*reg129; reg70=reg49*reg70; reg69=reg44+reg69;
    reg102=reg21+reg102; reg21=reg91*reg122; reg113=reg9-reg113; reg93=reg86+reg93; reg7=reg7*reg1;
    reg9=reg112*reg91; reg16=reg58*reg16; reg67=reg114+reg67; reg127=reg79*reg127; reg92=reg14+reg92;
    reg12=reg49*reg12; reg84=reg27+reg84; reg48=reg72*reg48; reg103=reg74+reg103; reg9=reg62-reg9;
    reg46=reg127*reg46; reg32=reg32*reg83; reg73=reg73*reg113; reg2=reg87*reg2; reg14=reg80*reg65;
    reg27=reg61*reg7; reg31=reg51+reg31; reg44=reg91*reg95; reg21=reg22-reg21; reg22=reg19*reg122;
    reg51=reg112*reg39; reg17=reg77+reg17; reg40=reg40*(*f.m).alpha_3; reg62=reg68*reg23; reg16=reg67+reg16;
    reg60=reg25*reg60; reg25=reg26*reg0; reg67=reg29*reg0; reg72=reg43*reg23; reg8=reg50*reg8;
    reg10=reg96*reg10; reg55=reg41+reg55; reg132=reg56+reg132; reg82=reg82*reg101; reg26=reg26*(*f.m).alpha_2;
    reg41=reg68*reg7; reg70=reg69+reg70; reg66=reg102+reg66; reg34=reg88+reg34; reg98=reg30+reg98;
    reg30=reg105*reg64; reg93=reg93*reg83; reg89=reg4+reg89; reg94=reg78+reg94; reg104=reg75+reg104;
    reg13=reg127*reg13; reg7=reg43*reg7; reg97=reg33+reg97; reg4=reg105*reg37; reg42=reg76+reg42;
    reg23=reg61*reg23; reg33=(*f.m).alpha_1*reg24; reg71=reg71*reg113; reg15=reg6+reg15; reg38=reg81+reg38;
    reg1=reg5*reg1; reg41=reg66+reg41; reg71=reg98+reg71; reg5=reg58*reg2; reg30=reg94+reg30;
    reg37=reg60*reg37; reg46=reg103+reg46; reg6=reg112*reg38; reg93=reg89+reg93; reg13=reg104+reg13;
    reg64=reg60*reg64; reg40=reg17+reg40; reg26=reg33+reg26; reg48=reg48*(*f.m).alpha_3; reg7=reg70+reg7;
    reg17=reg99*reg25; reg62=reg16+reg62; reg16=reg99*reg67; reg33=reg100*reg67; reg23=reg97+reg23;
    reg51=reg22-reg51; reg22=reg131*reg9; reg56=reg59*reg21; reg66=reg29*reg25; reg4=reg42+reg4;
    reg73=reg15+reg73; reg15=reg54*reg2; reg27=reg31+reg27; reg25=reg100*reg25; reg12=reg92+reg12;
    reg101=reg127*reg101; reg67=reg29*reg67; reg72=reg55+reg72; reg83=reg34*reg83; reg31=reg122*reg38;
    reg113=reg90*reg113; reg44=reg14+reg44; reg32=reg84+reg32; reg8=reg10+reg8; reg82=reg132+reg82;
    reg105=reg105*reg63; reg48=reg26+reg48; reg68=reg68*reg1; reg5=reg71+reg5; reg10=reg129*reg38;
    reg17=reg41+reg17; reg16=reg62+reg16; reg33=reg23+reg33; reg113=reg8+reg113; reg8=reg59*reg129;
    reg14=reg80*reg122; reg23=reg129*reg131; reg26=reg80*reg51; reg22=reg56-reg22; reg31=reg32+reg31;
    reg30=reg30*reg40; reg66=reg7+reg66; reg67=reg72+reg67; reg15=reg73+reg15; reg0=reg24*reg0;
    reg61=reg61*reg1; reg63=reg60*reg63; reg101=reg12+reg101; reg105=reg82+reg105; reg25=reg27+reg25;
    reg4=reg4*reg40; reg37=reg46+reg37; reg6=reg93+reg6; reg2=reg49*reg2; reg64=reg13+reg64;
    reg7=reg112*reg80; reg83=reg44+reg83; reg12=reg67*reg17; reg13=reg59*reg122; reg26=reg22+reg26;
    reg1=reg43*reg1; reg22=reg112*reg131; reg24=reg66*reg16; reg27=reg67*reg25; reg30=reg31+reg30;
    reg2=reg113+reg2; reg68=reg5+reg68; reg61=reg15+reg61; reg100=reg100*reg0; reg37=reg37*reg48;
    reg5=reg66*reg33; reg64=reg64*reg48; reg63=reg101+reg63; reg4=reg6+reg4; reg6=reg80*reg39;
    reg40=reg105*reg40; reg15=reg59*reg91; reg31=reg19*reg80; reg32=reg91*reg131; reg14=reg23-reg14;
    reg99=reg99*reg0; reg10=reg83+reg10; reg7=reg8-reg7; reg37=reg4+reg37; reg4=reg33*reg17;
    reg5=reg27-reg5; reg9=reg9/reg26; reg30=reg64+reg30; reg0=reg29*reg0; reg24=reg12-reg24;
    reg1=reg2+reg1; reg6=reg32-reg6; reg100=reg61+reg100; reg40=reg10+reg40; reg2=reg16*reg25;
    reg7=reg7/reg26; reg63=reg48*reg63; reg8=reg19*reg131; reg22=reg13-reg22; reg10=reg59*reg39;
    reg14=reg14/reg26; reg31=reg15-reg31; reg99=reg68+reg99; reg21=reg21/reg26; reg0=reg1+reg0;
    reg21=reg21*reg37; reg9=reg30*reg9; reg1=reg99*reg5; reg14=reg14*reg37; reg7=reg30*reg7;
    reg63=reg40+reg63; reg8=reg10-reg8; reg22=reg22/reg26; reg51=reg51/reg26; reg31=reg31/reg26;
    reg6=reg6/reg26; reg2=reg4-reg2; reg4=reg24*reg100; reg4=reg1-reg4; reg26=reg8/reg26;
    reg1=reg0*reg2; reg8=reg66*reg100; reg10=reg0*reg25; reg66=reg66*reg99; reg12=reg0*reg17;
    reg9=reg21-reg9; reg51=reg63*reg51; reg22=reg63*reg22; reg14=reg7-reg14; reg37=reg6*reg37;
    reg31=reg30*reg31; reg10=reg8-reg10; reg6=reg67*reg99; reg7=reg0*reg33; reg8=1-(*f.m).resolution;
    reg67=reg67*reg100; reg0=reg0*reg16; reg1=reg4+reg1; reg17=reg17*reg100; reg25=reg99*reg25;
    reg12=reg66-reg12; reg31=reg37-reg31; reg26=reg63*reg26; reg9=reg51+reg9; reg22=reg14-reg22;
    reg99=reg33*reg99; reg100=reg16*reg100; reg12=reg12/reg1; reg0=reg6-reg0; reg4=elem.pos(2)[1]-elem.pos(0)[1];
    reg17=reg25-reg17; reg6=elem.pos(1)[1]-elem.pos(0)[1]; reg10=reg10/reg1; reg13=(*f.m).resolution*reg95; reg7=reg67-reg7;
    reg22=reg8*reg22; reg31=reg26+reg31; reg9=reg8*reg9; reg14=elem.pos(1)[0]-elem.pos(0)[0]; reg15=(*f.m).resolution*reg65;
    reg16=elem.pos(2)[0]-elem.pos(0)[0]; reg5=reg5/reg1; reg31=reg8*reg31; reg100=reg99-reg100; reg22=reg13+reg22;
    reg13=reg14*reg4; reg21=reg16*reg6; reg17=reg17/reg1; reg7=reg7/reg1; reg0=reg0/reg1;
    reg80=reg80*reg8; reg91=reg91*reg8; reg24=reg24/reg1; reg15=reg9+reg15; reg9=(*f.m).resolution*reg12;
    reg23=(*f.m).resolution*reg10; reg25=(*f.m).resolution*reg38; reg2=reg2/reg1; reg22=(*f.m).deltaT*reg22; reg15=(*f.m).deltaT*reg15;
    reg31=reg25+reg31; reg9=reg91-reg9; reg23=reg80+reg23; reg25=(*f.m).resolution*reg17; reg26=(*f.m).resolution*reg0;
    reg27=(*f.m).resolution*reg7; reg29=(*f.m).resolution*reg5; reg30=(*f.m).resolution*reg24; reg129=reg129*reg8; reg39=reg8*reg39;
    reg131=reg8*reg131; reg19=reg19*reg8; reg59=reg59*reg8; reg21=reg13-reg21; reg1=reg100/reg1;
    reg13=(*f.m).resolution*reg2; reg4=reg4/reg21; reg32=reg9*reg22; reg33=reg23*reg15; reg112=reg112*reg8;
    reg122=reg8*reg122; reg29=reg59+reg29; reg30=reg19-reg30; reg8=(*f.m).resolution*reg1; reg27=reg131-reg27;
    reg26=reg39+reg26; reg16=reg16/reg21; reg6=reg6/reg21; reg25=reg129+reg25; reg31=(*f.m).deltaT*reg31;
    reg14=reg14/reg21; reg19=reg33+reg32; reg34=reg25*reg31; reg37=reg30*reg22; reg39=reg29*reg15;
    reg40=reg16-reg14; reg8=reg122-reg8; reg13=reg112+reg13; reg41=reg27*reg15; reg42=reg6-reg4;
    reg43=reg26*reg22; reg44=0.5*reg6; reg46=0.5*reg42; reg48=0.5*reg14; reg51=0.5*reg4;
    reg55=reg13*reg31; reg56=reg39+reg37; reg59=0.5*reg40; reg60=reg34+reg19; reg61=0.5*reg16;
    reg62=reg41+reg43; reg63=reg8*reg31; reg64=reg4*reg23; reg66=reg62+reg63; reg67=reg25*reg51;
    reg68=reg9*reg16; reg69=reg25*reg48; reg70=reg23*reg6; reg71=reg9*reg14; reg72=reg44*reg25;
    reg73=reg59*reg25; reg74=reg42*reg23; reg75=reg55+reg56; reg76=reg25*reg46; reg77=reg40*reg9;
    reg78=reg25*reg61; reg79=2*reg60; reg80=reg8*reg48; reg81=reg27*reg6; reg69=reg69-reg70;
    reg82=reg59*reg8; reg83=reg42*reg27; reg84=reg8*reg46; reg85=reg44*reg13; reg71=reg71-reg72;
    reg86=reg59*reg13; reg87=reg30*reg14; reg88=reg4*reg27; reg89=reg42*reg29; reg90=reg8*reg61;
    reg91=reg4*reg29; reg92=reg13*reg61; reg93=reg4*reg75; reg74=reg73+reg74; reg73=reg13*reg46;
    reg94=reg40*reg30; reg77=reg76+reg77; reg76=reg26*reg16; reg97=reg8*reg51; reg98=reg6*reg75;
    reg99=1-var_inter[0]; reg100=reg79*reg61; reg101=reg79*reg48; reg102=reg14*reg66; reg103=reg16*reg66;
    reg104=reg51*reg79; reg105=reg13*reg51; reg106=reg30*reg16; reg107=reg44*reg79; reg64=reg64-reg78;
    reg67=reg67-reg68; reg108=reg44*reg8; reg109=reg40*reg26; reg110=reg13*reg48; reg111=reg29*reg6;
    reg112=reg26*reg14; reg113=reg107-reg102; reg114=reg100-reg93; reg115=reg40*reg66; reg94=reg73+reg94;
    reg74=2*reg74; reg64=2*reg64; reg73=reg98-reg101; reg89=reg86+reg89; reg86=reg103-reg104;
    reg116=var_inter[1]*elem.f_vol_e[1]; reg117=var_inter[0]*elem.f_vol_e[0]; reg105=reg105-reg106; reg118=var_inter[1]*elem.f_vol_e[0]; reg119=var_inter[0]*elem.f_vol_e[1];
    reg120=reg79*reg46; reg109=reg84+reg109; reg69=2*reg69; reg67=2*reg67; reg110=reg110-reg111;
    reg80=reg80-reg81; reg91=reg91-reg92; reg88=reg88-reg90; reg99=reg99-var_inter[1]; reg77=2*reg77;
    reg84=reg42*reg75; reg87=reg87-reg85; reg121=reg59*reg79; reg112=reg112-reg108; reg97=reg97-reg76;
    reg82=reg83+reg82; reg71=2*reg71; reg113=reg113-reg116; reg83=reg69*reg48; reg122=reg67*reg48;
    reg123=reg91*reg6; reg124=reg105*reg6; reg125=reg64*reg48; reg73=reg73-reg118; reg126=reg89*reg4;
    reg127=reg74*reg61; reg128=reg94*reg4; reg129=reg94*reg6; reg130=reg77*reg48; reg86=reg86-reg119;
    reg131=reg105*reg4; reg132=reg67*reg61; T reg133=reg89*reg6; reg114=reg114-reg117; T reg134=reg77*reg61;
    T reg135=reg74*reg48; T reg136=reg44*reg74; T reg137=reg44*reg69; T reg138=reg14*reg80; T reg139=reg109*reg14;
    T reg140=reg82*reg14; T reg141=reg40*reg97; T reg142=reg88*reg14; T reg143=reg44*reg64; T reg144=reg67*reg46;
    T reg145=reg87*reg6; T reg146=reg14*reg97; T reg147=reg44*reg67; T reg148=elem.f_vol_e[1]*reg99; T reg149=reg40*reg80;
    T reg150=reg69*reg46; T reg151=reg71*reg48; T reg152=reg42*reg91; T reg153=reg40*reg112; T reg154=reg71*reg46;
    T reg155=reg59*reg64; T reg156=reg110*reg6; T reg157=reg84+reg121; T reg158=reg82*reg16; T reg159=reg59*reg69;
    T reg160=reg77*reg51; T reg161=reg59*reg71; T reg162=reg99*elem.f_vol_e[0]; T reg163=reg16*reg112; T reg164=reg42*reg87;
    T reg165=reg109*reg16; reg89=reg42*reg89; T reg166=reg40*reg88; T reg167=reg64*reg51; reg88=reg88*reg16;
    T reg168=reg59*reg74; T reg169=reg71*reg51; T reg170=reg59*reg77; reg94=reg42*reg94; reg112=reg14*reg112;
    T reg171=reg44*reg71; reg80=reg16*reg80; T reg172=reg69*reg51; T reg173=reg64*reg46; reg97=reg16*reg97;
    T reg174=reg67*reg51; reg91=reg91*reg4; reg64=reg64*reg61; T reg175=reg110*reg4; reg105=reg42*reg105;
    reg69=reg69*reg61; reg87=reg4*reg87; T reg176=reg77*reg46; reg67=reg59*reg67; reg77=reg44*reg77;
    T reg177=reg115+reg120; reg109=reg109*reg40; reg82=reg40*reg82; T reg178=reg74*reg46; reg110=reg42*reg110;
    reg71=reg71*reg61; reg74=reg74*reg51; T reg179=reg162+reg157; reg124=reg122-reg124; reg147=reg146-reg147;
    reg97=reg174-reg97; reg158=reg74-reg158; reg133=reg135-reg133; reg132=reg131-reg132; reg71=reg87-reg71;
    reg69=reg175-reg69; reg156=reg83-reg156; reg145=reg151-reg145; reg80=reg172-reg80; reg74=reg148+reg177;
    reg88=reg167-reg88; reg163=reg169-reg163; reg123=reg125-reg123; reg155=reg152+reg155; reg154=reg153+reg154;
    reg113=reg21*reg113; reg73=reg21*reg73; reg127=reg126-reg127; reg86=reg21*reg86; reg134=reg128-reg134;
    reg67=reg105+reg67; reg159=reg110+reg159; reg114=reg21*reg114; reg168=reg89+reg168; reg64=reg91-reg64;
    reg171=reg112-reg171; reg173=reg166+reg173; reg170=reg94+reg170; reg176=reg109+reg176; reg165=reg160-reg165;
    reg161=reg164+reg161; reg178=reg82+reg178; reg137=reg138-reg137; reg129=reg130-reg129; reg143=reg142-reg143;
    reg77=reg139-reg77; reg136=reg140-reg136; reg144=reg141+reg144; reg150=reg149+reg150; reg132=reg132*reg21;
    reg163=reg163*reg21; reg77=reg77*reg21; reg168=reg168*reg21; reg80=reg80*reg21; reg64=reg21*reg64;
    reg176=reg176*reg21; reg171=reg21*reg171; reg97=reg97*reg21; reg69=reg69*reg21; reg173=reg173*reg21;
    reg147=reg21*reg147; reg88=reg88*reg21; reg170=reg170*reg21; reg136=reg136*reg21; reg165=reg165*reg21;
    reg137=reg21*reg137; reg178=reg178*reg21; reg158=reg158*reg21; reg161=reg161*reg21; reg71=reg71*reg21;
    reg156=reg156*reg21; reg155=reg155*reg21; reg124=reg124*reg21; reg154=reg21*reg154; reg113=ponderation*reg113;
    reg150=reg21*reg150; reg82=reg21*reg179; reg123=reg123*reg21; reg73=ponderation*reg73; reg143=reg143*reg21;
    reg83=reg21*reg74; reg129=reg129*reg21; reg159=reg159*reg21; reg144=reg21*reg144; reg134=reg21*reg134;
    reg67=reg67*reg21; reg114=ponderation*reg114; reg133=reg133*reg21; reg127=reg21*reg127; reg86=ponderation*reg86;
    reg145=reg145*reg21; T tmp_2_4=ponderation*reg69; T tmp_2_0=ponderation*reg127; T tmp_1_1=ponderation*reg176; reg69=ponderation*reg82;
    T vec_0=reg69; T tmp_1_0=ponderation*reg178; T tmp_1_3=ponderation*reg144; T tmp_1_5=ponderation*reg154; T tmp_2_3=ponderation*reg132;
    T tmp_2_2=ponderation*reg64; T tmp_2_1=ponderation*reg134; T tmp_1_4=ponderation*reg150; T tmp_5_3=ponderation*reg147; T tmp_5_0=ponderation*reg136;
    T tmp_1_2=ponderation*reg173; T tmp_5_1=ponderation*reg77; reg64=ponderation*reg83; T vec_1=reg64; T tmp_4_5=ponderation*reg145;
    T tmp_5_2=ponderation*reg143; T tmp_4_4=ponderation*reg156; T tmp_4_3=ponderation*reg124; T vec_5=-reg113; T tmp_4_2=ponderation*reg123;
    T vec_4=-reg73; T tmp_4_1=ponderation*reg129; T vec_3=-reg86; T tmp_4_0=ponderation*reg133; T vec_2=-reg114;
    T tmp_0_2=ponderation*reg155; T tmp_0_3=ponderation*reg67; T tmp_0_4=ponderation*reg159; T tmp_3_5=ponderation*reg163; T tmp_5_4=ponderation*reg137;
    T tmp_3_4=ponderation*reg80; T tmp_0_0=ponderation*reg168; T tmp_3_3=ponderation*reg97; T tmp_5_5=ponderation*reg171; T tmp_3_2=ponderation*reg88;
    T tmp_0_1=ponderation*reg170; T tmp_3_1=ponderation*reg165; T tmp_3_0=ponderation*reg158; T tmp_0_5=ponderation*reg161; T tmp_2_5=ponderation*reg71;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_3;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_4;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_5;
    matrix(indices[0]+1,indices[0]+0) += tmp_1_0;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_3;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_4;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_5;
    matrix(indices[1]+0,indices[0]+0) += tmp_2_0;
    matrix(indices[1]+0,indices[0]+1) += tmp_2_1;
    matrix(indices[1]+0,indices[1]+0) += tmp_2_2;
    matrix(indices[1]+0,indices[1]+1) += tmp_2_3;
    matrix(indices[1]+0,indices[2]+0) += tmp_2_4;
    matrix(indices[1]+0,indices[2]+1) += tmp_2_5;
    matrix(indices[1]+1,indices[0]+0) += tmp_3_0;
    matrix(indices[1]+1,indices[0]+1) += tmp_3_1;
    matrix(indices[1]+1,indices[1]+0) += tmp_3_2;
    matrix(indices[1]+1,indices[1]+1) += tmp_3_3;
    matrix(indices[1]+1,indices[2]+0) += tmp_3_4;
    matrix(indices[1]+1,indices[2]+1) += tmp_3_5;
    matrix(indices[2]+0,indices[0]+0) += tmp_4_0;
    matrix(indices[2]+0,indices[0]+1) += tmp_4_1;
    matrix(indices[2]+0,indices[1]+0) += tmp_4_2;
    matrix(indices[2]+0,indices[1]+1) += tmp_4_3;
    matrix(indices[2]+0,indices[2]+0) += tmp_4_4;
    matrix(indices[2]+0,indices[2]+1) += tmp_4_5;
    matrix(indices[2]+1,indices[0]+0) += tmp_5_0;
    matrix(indices[2]+1,indices[0]+1) += tmp_5_1;
    matrix(indices[2]+1,indices[1]+0) += tmp_5_2;
    matrix(indices[2]+1,indices[1]+1) += tmp_5_3;
    matrix(indices[2]+1,indices[2]+0) += tmp_5_4;
    matrix(indices[2]+1,indices[2]+1) += tmp_5_5;
sollicitation[indices[0]+0] += vec_0;
sollicitation[indices[0]+1] += vec_1;
sollicitation[indices[1]+0] += vec_2;
sollicitation[indices[1]+1] += vec_3;
sollicitation[indices[2]+0] += vec_4;
sollicitation[indices[2]+1] += vec_5;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; reg1=1.0/reg1; T reg2=2*(*f.m).shear_modulus_12;
    T reg3=reg1*reg0; reg2=1.0/reg2; T reg4=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg5=1.0/(*f.m).elastic_modulus_3; T reg6=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg7=reg2*reg3; T reg8=pow((*f.m).v1[0],2); T reg9=pow((*f.m).v1[1],2); T reg10=pow((*f.m).v2[0],2); T reg11=pow((*f.m).v2[1],2);
    reg9=reg8+reg9; reg8=pow((*f.m).v1[2],2); T reg12=1.0/(*f.m).elastic_modulus_2; reg11=reg10+reg11; reg10=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1;
    T reg13=pow((*f.m).v2[2],2); T reg14=reg6*reg7; T reg15=reg4*reg7; T reg16=reg5*reg7; T reg17=reg4*reg14;
    reg13=reg11+reg13; reg11=reg10*reg16; reg8=reg9+reg8; reg9=reg4*reg15; T reg18=reg12*reg16;
    reg17=reg11+reg17; reg9=reg18-reg9; reg13=pow(reg13,0.5); reg8=pow(reg8,0.5); reg18=reg10*reg15;
    T reg19=1.0/(*f.m).elastic_modulus_1; T reg20=reg12*reg14; T reg21=reg19*reg9; T reg22=(*f.m).v2[2]/reg13; T reg23=reg18+reg20;
    T reg24=(*f.m).v2[1]/reg13; T reg25=(*f.m).v1[1]/reg8; T reg26=(*f.m).v1[2]/reg8; T reg27=reg10*reg17; reg13=(*f.m).v2[0]/reg13;
    T reg28=reg10*reg7; T reg29=reg2*reg0; reg16=reg19*reg16; T reg30=reg6*reg14; T reg31=reg6*reg23;
    T reg32=reg5*reg3; reg27=reg21-reg27; reg21=reg26*reg24; T reg33=reg6*reg15; reg8=(*f.m).v1[0]/reg8;
    reg7=reg12*reg7; T reg34=reg4*reg3; T reg35=reg25*reg22; reg3=reg6*reg3; reg30=reg16-reg30;
    reg15=reg19*reg15; reg16=reg28*reg6; T reg36=reg3*reg4; T reg37=reg29*reg4; T reg38=reg34*reg4;
    T reg39=reg32*reg10; reg14=reg10*reg14; T reg40=reg29*reg6; reg32=reg32*reg12; T reg41=reg2*reg1;
    reg29=reg29*reg5; T reg42=2*reg13; reg31=reg27-reg31; reg27=reg35-reg21; T reg43=2*reg8;
    reg33=reg11+reg33; reg11=reg8*reg22; T reg44=reg26*reg13; T reg45=reg6*reg7; reg33=reg33/reg31;
    reg28=reg28*reg10; T reg46=reg44-reg11; reg14=reg15+reg14; T reg47=reg40*reg4; T reg48=reg41*reg5;
    reg45=reg18+reg45; reg16=reg15+reg16; reg17=reg17/reg31; reg30=reg30/reg31; reg9=reg9/reg31;
    reg7=reg19*reg7; reg15=pow(reg8,2); T reg49=pow(reg25,2); T reg50=reg43*reg25; T reg51=reg25*reg13;
    T reg52=reg8*reg24; T reg53=pow(reg13,2); T reg54=pow(reg24,2); T reg55=reg42*reg24; T reg56=reg37*reg4;
    T reg57=reg29*reg10; reg29=reg29*reg12; T reg58=2*reg27; reg34=reg34*reg10; T reg59=reg41*reg4;
    reg36=reg39+reg36; reg3=reg3*reg12; reg38=reg32-reg38; reg41=reg41*reg6; reg32=reg54*reg17;
    reg39=reg49*reg9; T reg60=reg41*reg4; T reg61=reg59*reg4; T reg62=reg48*reg10; reg48=reg48*reg12;
    T reg63=reg52-reg51; reg57=reg47+reg57; reg14=reg14/reg31; reg40=reg40*reg12; reg37=reg37*reg10;
    reg28=reg7-reg28; reg45=reg45/reg31; reg7=reg50*reg9; reg47=reg55*reg17; T reg64=reg15*reg33;
    T reg65=reg53*reg30; T reg66=reg49*reg33; T reg67=reg54*reg30; T reg68=reg50*reg33; T reg69=reg55*reg30;
    reg38=reg38*reg19; reg56=reg29-reg56; reg23=reg23/reg31; reg36=reg36*reg10; reg16=reg16/reg31;
    reg29=reg34+reg3; T reg70=reg53*reg17; T reg71=reg15*reg9; T reg72=pow(reg46,2); T reg73=pow(reg27,2);
    T reg74=pow(reg26,2); T reg75=pow(reg22,2); T reg76=reg58*reg46; T reg77=reg73*reg14; reg65=reg64+reg65;
    reg64=reg76*reg23; reg47=reg7+reg47; reg7=reg75*reg17; T reg78=reg74*reg9; T reg79=reg72*reg23;
    reg32=reg39+reg32; reg60=reg62+reg60; reg61=reg48-reg61; reg56=reg56*reg19; reg39=reg37+reg40;
    reg57=reg57*reg10; reg48=reg76*reg14; reg62=reg15*reg45; T reg80=reg53*reg16; T reg81=reg49*reg45;
    T reg82=reg54*reg16; T reg83=reg55*reg16; T reg84=reg50*reg45; reg70=reg71+reg70; reg71=reg73*reg23;
    reg29=reg29*reg6; reg69=reg68+reg69; reg36=reg38-reg36; reg38=pow(reg63,2); reg41=reg41*reg12;
    reg59=reg59*reg10; reg68=reg75*reg30; reg67=reg66+reg67; reg66=reg72*reg14; T reg85=reg74*reg33;
    reg28=reg28/reg31; reg83=reg84+reg83; reg61=reg61*reg19; reg57=reg56-reg57; reg76=reg76*reg28;
    reg56=reg75*reg16; reg39=reg39*reg6; reg84=reg74*reg45; T reg86=2*reg25; T reg87=reg43*reg26;
    reg48=reg69+reg48; reg69=reg72*reg28; reg82=reg81+reg82; reg80=reg62+reg80; reg62=reg8*reg25;
    reg81=reg13*reg24; T reg88=reg73*reg28; reg71=reg70+reg71; reg68=reg85+reg68; reg64=reg47+reg64;
    reg47=reg38*reg23; reg7=reg78+reg7; reg70=reg38*reg14; reg78=reg42*reg22; reg79=reg32+reg79;
    reg29=reg36-reg29; reg32=reg59+reg41; reg60=reg60*reg10; reg36=2*reg24; reg66=reg67+reg66;
    reg77=reg65+reg77; reg65=reg53*reg77; reg67=reg15*reg64; reg85=reg53*reg66; T reg89=reg78*reg17;
    T reg90=reg8*reg46; T reg91=reg25*reg27; reg69=reg82+reg69; reg82=reg87*reg9; T reg92=reg15*reg79;
    T reg93=reg15*reg71; reg52=reg51+reg52; reg60=reg61-reg60; reg58=reg58*reg63; reg51=2*reg46;
    reg39=reg57-reg39; reg32=reg32*reg6; reg29=reg29/reg31; reg57=reg86*reg26; reg61=reg8*reg13;
    T reg94=reg36*reg22; T reg95=reg25*reg24; reg70=reg68+reg70; reg68=reg87*reg33; T reg96=reg78*reg30;
    reg47=reg7+reg47; reg7=reg53*reg48; T reg97=reg46*reg27; reg88=reg80+reg88; reg80=reg66*reg81;
    T reg98=reg49*reg64; T reg99=reg54*reg48; T reg100=reg79*reg62; reg79=reg49*reg79; reg48=reg48*reg81;
    reg64=reg64*reg62; T reg101=2*reg26; T reg102=reg86*reg24; reg66=reg54*reg66; T reg103=reg43*reg13;
    T reg104=reg54*reg77; T reg105=reg49*reg71; reg56=reg84+reg56; reg84=reg38*reg28; reg76=reg83+reg76;
    reg83=reg53*reg10; T reg106=reg26*reg22; reg39=reg39/reg31; reg77=reg77*reg81; T reg107=reg78*reg16;
    reg51=reg51*reg63; T reg108=reg8*reg27; T reg109=reg69*reg97; T reg110=reg10*reg102; T reg111=reg19*reg103;
    reg104=reg105+reg104; reg105=reg25*reg46; T reg112=reg54*reg10; T reg113=reg72*reg88; T reg114=reg53*reg19;
    reg7=reg67+reg7; reg67=reg87*reg45; T reg115=reg72*reg76; reg32=reg60-reg32; reg90=reg91+reg90;
    reg60=reg73*reg76; reg99=reg98+reg99; reg84=reg56+reg84; reg80=reg100+reg80; reg76=reg76*reg97;
    reg48=reg64+reg48; reg85=reg92+reg85; reg56=reg73*reg69; reg64=reg54*reg70; reg91=reg73*reg88;
    reg65=reg93+reg65; reg92=reg24*reg27; reg93=reg13*reg46; reg98=reg15*reg47; reg100=reg53*reg70;
    reg69=reg72*reg69; reg17=reg94*reg17; reg30=reg94*reg30; reg33=reg57*reg33; reg101=reg22*reg101;
    T reg116=reg58*reg14; reg96=reg68+reg96; reg71=reg62*reg71; reg68=reg61*reg29; T reg117=reg95*reg29;
    reg9=reg57*reg9; T reg118=reg52*reg29; T reg119=reg12*reg102; reg66=reg79+reg66; reg89=reg82+reg89;
    reg79=reg54*reg12; reg82=reg10*reg103; T reg120=reg49*reg47; T reg121=reg58*reg23; reg69=reg66+reg69;
    reg66=reg117*reg102; T reg122=reg52*reg117; reg109=reg80+reg109; reg47=reg47*reg62; reg80=reg68*reg103;
    reg91=reg65+reg91; reg56=reg85+reg56; reg117=reg117*reg103; reg100=reg98+reg100; reg65=reg73*reg84;
    reg82=reg119-reg82; reg85=reg4*reg101; reg83=reg79-reg83; reg79=reg75*reg4; reg77=reg71+reg77;
    reg88=reg88*reg97; reg110=reg111-reg110; reg71=reg6*reg101; reg112=reg114-reg112; reg98=reg75*reg6;
    reg113=reg104+reg113; reg104=reg68*reg102; reg86=reg86*reg46; reg43=reg43*reg27; reg111=reg4*reg102;
    reg114=reg6*reg103; reg119=reg54*reg4; T reg123=reg39*reg90; reg105=reg39*reg105; T reg124=reg53*reg6;
    reg108=reg39*reg108; T reg125=reg106*reg29; reg116=reg96+reg116; reg30=reg33+reg30; reg14=reg51*reg14;
    reg93=reg92+reg93; reg33=reg24*reg46; reg92=reg13*reg27; reg16=reg94*reg16; reg45=reg57*reg45;
    reg58=reg58*reg28; reg107=reg67+reg107; reg67=reg26*reg63; reg70=reg70*reg81; reg23=reg51*reg23;
    reg17=reg9+reg17; reg121=reg89+reg121; reg64=reg120+reg64; reg9=reg72*reg84; reg76=reg48+reg76;
    reg48=reg52*reg118; reg60=reg7+reg60; reg115=reg99+reg115; reg7=reg118*reg102; reg118=reg118*reg103;
    reg31=reg32/reg31; reg19=reg15*reg19; reg32=reg54*(*f.m).alpha_2; reg89=(*f.m).alpha_1*reg49; reg28=reg51*reg28;
    reg111=reg114+reg111; reg51=reg49*reg121; reg16=reg45+reg16; reg45=reg31*reg93; reg58=reg107+reg58;
    reg96=reg53*(*f.m).alpha_2; reg99=reg105*reg86; reg66=reg69+reg66; reg80=reg91+reg80; reg69=reg54*reg116;
    reg33=reg31*reg33; reg14=reg30+reg14; reg9=reg64+reg9; reg30=reg125*reg102; reg64=(*f.m).alpha_1*reg15;
    reg44=reg11+reg44; reg101=reg5*reg101; reg92=reg31*reg92; reg11=reg123*reg86; reg91=reg108*reg86;
    reg104=reg113+reg104; reg67=reg39*reg67; reg7=reg115+reg7; reg107=reg75*reg5; reg119=reg124+reg119;
    reg12=reg49*reg12; reg113=reg15*reg10; reg71=reg110-reg71; reg98=reg112-reg98; reg110=reg22*reg63;
    reg122=reg109+reg122; reg109=reg105*reg90; reg10=reg49*reg10; reg70=reg47+reg70; reg84=reg84*reg97;
    reg23=reg17+reg23; reg118=reg60+reg118; reg17=reg123*reg43; reg42=reg42*reg27; reg36=reg36*reg46;
    reg48=reg76+reg48; reg123=reg123*reg90; reg47=reg53*reg116; reg60=reg15*reg121; reg76=reg108*reg43;
    reg117=reg56+reg117; reg105=reg105*reg43; reg56=reg26*reg27; reg112=reg8*reg63; reg65=reg100+reg65;
    reg100=reg125*reg103; reg85=reg82-reg85; reg79=reg83-reg79; reg88=reg77+reg88; reg68=reg52*reg68;
    reg77=reg61*reg98; reg82=reg45*reg36; reg69=reg51+reg69; reg51=reg15*reg98; reg83=reg49*reg79;
    reg114=reg95*reg85; reg98=reg53*reg98; reg115=reg54*reg79; reg120=reg72*reg58; reg124=reg25*reg63;
    reg110=reg31*reg110; reg79=reg95*reg79; T reg126=reg44*reg29; T reg127=reg61*reg71; T reg128=reg13*reg63;
    reg112=reg56+reg112; reg56=reg22*reg27; reg28=reg16+reg28; reg16=reg26*reg46; reg96=reg64+reg96;
    reg64=reg73*(*f.m).alpha_3; T reg129=reg15*reg23; T reg130=reg53*reg14; reg32=reg89+reg32; reg89=reg72*(*f.m).alpha_3;
    T reg131=(*f.m).alpha_1*reg74; T reg132=reg75*(*f.m).alpha_2; reg68=reg88+reg68; reg108=reg108*reg90; reg91=reg104+reg91;
    reg88=reg92*reg36; reg109=reg122+reg109; reg104=reg33*reg93; reg99=reg66+reg99; reg66=reg33*reg36;
    reg84=reg70+reg84; reg125=reg52*reg125; reg30=reg9+reg30; reg9=reg67*reg86; reg70=reg15*reg71;
    reg122=reg49*reg85; reg123=reg48+reg123; reg48=reg45*reg93; reg121=reg121*reg62; reg116=reg116*reg81;
    reg11=reg7+reg11; reg100=reg65+reg100; reg7=reg74*reg6; reg10=reg19-reg10; reg19=reg67*reg43;
    reg119=reg107-reg119; reg65=reg74*reg4; reg113=reg12-reg113; reg45=reg45*reg42; reg6=reg15*reg6;
    reg17=reg118+reg17; reg4=reg49*reg4; reg111=reg101-reg111; reg47=reg60+reg47; reg35=reg21+reg35;
    reg71=reg53*reg71; reg12=reg54*reg14; reg21=reg49*reg23; reg33=reg33*reg42; reg105=reg117+reg105;
    reg60=reg73*reg58; reg101=reg92*reg42; reg76=reg80+reg76; reg85=reg54*reg85; reg80=reg75*reg111;
    reg9=reg30+reg9; reg122=reg70+reg122; reg30=reg74*reg111; reg85=reg71+reg85; reg70=reg110*reg36;
    reg71=reg52*reg2; reg107=reg81*reg2; reg48=reg123+reg48; reg116=reg121+reg116; reg45=reg17+reg45;
    reg82=reg11+reg82; reg58=reg58*reg97; reg115=reg98+reg115; reg11=reg75*reg119; reg120=reg69+reg120;
    reg17=reg126*reg102; reg72=reg72*reg28; reg64=reg96+reg64; reg101=reg76+reg101; reg33=reg105+reg33;
    reg89=reg32+reg89; reg132=reg131+reg132; reg38=reg38*(*f.m).alpha_3; reg32=(*f.m).alpha_1*reg62; reg69=reg81*(*f.m).alpha_2;
    reg130=reg129+reg130; reg73=reg73*reg28; reg65=reg113-reg65; reg108=reg68+reg108; reg92=reg92*reg93;
    reg7=reg10-reg7; reg88=reg91+reg88; reg104=reg109+reg104; reg19=reg100+reg19; reg10=reg110*reg42;
    reg66=reg99+reg66; reg125=reg84+reg125; reg67=reg67*reg90; reg68=reg22*reg46; reg8=reg8*reg26;
    reg29=reg35*reg29; reg76=reg24*reg63; reg4=reg6+reg4; reg111=reg106*reg111; reg124=reg16+reg124;
    reg81=reg14*reg81; reg114=reg127+reg114; reg112=reg39*reg112; reg12=reg21+reg12; reg5=reg74*reg5;
    reg23=reg23*reg62; reg128=reg56+reg128; reg13=reg13*reg22; reg6=reg74*reg119; reg14=reg126*reg103;
    reg60=reg47+reg60; reg83=reg51+reg83; reg79=reg77+reg79; reg119=reg106*reg119; reg22=reg24*reg22;
    reg26=reg25*reg26; reg92=reg108+reg92; reg16=reg52*reg71; reg111=reg114+reg111; reg10=reg19+reg10;
    reg4=reg5-reg4; reg5=reg52*reg107; reg27=reg63*reg27; reg72=reg12+reg72; reg12=reg13*(*f.m).alpha_2;
    reg19=(*f.m).alpha_1*reg8; reg102=reg29*reg102; reg21=reg97*(*f.m).alpha_3; reg69=reg32+reg69; reg38=reg132+reg38;
    reg24=reg50*reg107; reg103=reg29*reg103; reg6=reg83+reg6; reg73=reg130+reg73; reg128=reg31*reg128;
    reg25=reg45*reg104; reg32=reg49*reg65; reg97=reg28*reg97; reg28=reg66*reg89; reg47=reg112*reg43;
    reg51=reg53*reg7; reg56=reg54*reg65; reg76=reg68+reg76; reg68=reg101*reg64; reg81=reg23+reg81;
    reg58=reg116+reg58; reg124=reg39*reg124; reg23=reg33*reg89; reg14=reg60+reg14; reg11=reg115+reg11;
    reg107=reg55*reg107; reg126=reg52*reg126; reg17=reg120+reg17; reg39=reg112*reg86; reg60=reg88*reg64;
    reg77=reg104*reg82; reg83=reg33*reg48; reg80=reg85+reg80; reg84=reg48*reg66; reg30=reg122+reg30;
    reg110=reg110*reg93; reg67=reg125+reg67; reg85=reg50*reg71; reg119=reg79+reg119; reg79=reg44*reg1;
    reg91=reg15*reg7; reg70=reg9+reg70; reg71=reg55*reg71; reg13=reg13*reg1; reg32=reg91+reg32;
    reg74=reg74*reg4; reg86=reg124*reg86; reg28=reg60+reg28; reg39=reg17+reg39; reg9=reg128*reg36;
    reg24=reg6+reg24; reg103=reg73+reg103; reg43=reg124*reg43; reg6=reg87*reg13; reg102=reg72+reg102;
    reg70=reg70*reg38; reg2=reg62*reg2; reg46=reg63*reg46; reg17=reg35*reg0; reg60=reg22*reg0;
    reg62=reg104*reg89; reg63=reg92*reg64; reg85=reg30+reg85; reg30=reg87*reg79; reg56=reg51+reg56;
    reg75=reg75*reg4; reg107=reg11+reg107; reg11=reg78*reg13; reg97=reg81+reg97; reg29=reg52*reg29;
    reg10=reg10*reg38; reg23=reg68+reg23; reg51=reg78*reg79; reg71=reg80+reg71; reg76=reg31*reg76;
    reg77=reg84-reg77; reg25=reg83-reg25; reg31=reg33*reg82; reg68=reg45*reg66; reg79=reg44*reg79;
    reg16=reg111+reg16; reg13=reg44*reg13; reg5=reg119+reg5; reg65=reg95*reg65; reg7=reg61*reg7;
    reg112=reg112*reg90; reg126=reg58+reg126; reg110=reg67+reg110; reg58=reg128*reg42; reg47=reg14+reg47;
    reg27=reg27*(*f.m).alpha_3; reg12=reg19+reg12; reg21=reg69+reg21; reg14=(*f.m).alpha_1*reg26; reg22=reg22*(*f.m).alpha_2;
    reg79=reg16+reg79; reg16=reg35*reg17; reg51=reg71+reg51; reg19=reg94*reg17; reg11=reg107+reg11;
    reg67=reg94*reg60; reg69=reg88*reg25; reg71=reg82*reg21; reg70=reg28+reg70; reg9=reg39+reg9;
    reg29=reg97+reg29; reg90=reg124*reg90; reg74=reg32+reg74; reg28=reg50*reg2; reg6=reg24+reg6;
    reg68=reg31-reg68; reg27=reg12+reg27; reg36=reg76*reg36; reg12=reg45*reg21; reg10=reg23+reg10;
    reg62=reg63+reg62; reg86=reg102+reg86; reg112=reg126+reg112; reg128=reg128*reg93; reg23=reg101*reg77;
    reg24=reg57*reg60; reg65=reg7+reg65; reg4=reg106*reg4; reg38=reg110*reg38; reg22=reg14+reg22;
    reg43=reg103+reg43; reg30=reg85+reg30; reg17=reg57*reg17; reg42=reg76*reg42; reg13=reg5+reg13;
    reg60=reg35*reg60; reg75=reg56+reg75; reg5=reg55*reg2; reg1=reg8*reg1; reg46=reg46*(*f.m).alpha_3;
    reg58=reg47+reg58; reg28=reg74+reg28; reg87=reg87*reg1; reg69=reg23-reg69; reg7=reg92*reg68;
    reg8=reg48*reg88; reg14=reg92*reg82; reg23=reg101*reg48; reg31=reg45*reg92; reg71=reg70+reg71;
    reg9=reg9*reg27; reg16=reg79+reg16; reg42=reg43+reg42; reg60=reg13+reg60; reg2=reg52*reg2;
    reg4=reg65+reg4; reg128=reg112+reg128; reg36=reg86+reg36; reg58=reg58*reg27; reg12=reg10+reg12;
    reg17=reg30+reg17; reg46=reg22+reg46; reg24=reg6+reg24; reg5=reg75+reg5; reg78=reg78*reg1;
    reg67=reg11+reg67; reg90=reg29+reg90; reg38=reg62+reg38; reg19=reg51+reg19; reg0=reg26*reg0;
    reg6=reg48*reg21; reg93=reg76*reg93; reg10=reg101*reg104; reg1=reg44*reg1; reg11=reg16*reg24;
    reg13=reg92*reg66; reg14=reg8-reg14; reg8=reg16*reg67; reg78=reg5+reg78; reg9=reg71+reg9;
    reg5=reg60*reg19; reg31=reg23-reg31; reg94=reg94*reg0; reg22=reg33*reg92; reg23=reg104*reg88;
    reg93=reg90+reg93; reg2=reg4+reg2; reg36=reg36*reg46; reg6=reg38+reg6; reg42=reg42*reg46;
    reg4=reg60*reg17; reg58=reg12+reg58; reg12=reg45*reg88; reg7=reg69+reg7; reg27=reg128*reg27;
    reg57=reg57*reg0; reg87=reg28+reg87; reg26=reg101*reg82; reg94=reg78+reg94; reg27=reg6+reg27;
    reg5=reg8-reg5; reg22=reg10-reg22; reg93=reg46*reg93; reg9=reg36+reg9; reg0=reg35*reg0;
    reg1=reg2+reg1; reg2=reg101*reg66; reg6=reg17*reg67; reg12=reg26-reg12; reg4=reg11-reg4;
    reg8=reg19*reg24; reg77=reg77/reg7; reg10=reg33*reg88; reg31=reg31/reg7; reg42=reg58+reg42;
    reg57=reg87+reg57; reg25=reg25/reg7; reg13=reg23-reg13; reg14=reg14/reg7; reg25=reg9*reg25;
    reg68=reg68/reg7; reg12=reg12/reg7; reg14=reg14*reg42; reg31=reg9*reg31; reg11=reg4*reg94;
    reg77=reg77*reg42; reg10=reg2-reg10; reg13=reg13/reg7; reg0=reg1+reg0; reg93=reg27+reg93;
    reg6=reg8-reg6; reg1=reg57*reg5; reg22=reg22/reg7; reg7=reg10/reg7; reg2=reg60*reg57;
    reg8=reg0*reg67; reg60=reg60*reg94; reg11=reg1-reg11; reg1=reg0*reg24; reg10=reg0*reg6;
    reg68=reg93*reg68; reg25=reg77-reg25; reg12=reg93*reg12; reg14=reg31-reg14; reg42=reg13*reg42;
    reg22=reg9*reg22; reg9=reg0*reg17; reg13=reg16*reg57; reg0=reg0*reg19; reg8=reg60-reg8;
    reg16=reg16*reg94; reg10=reg11+reg10; reg11=1-(*f.m).resolution; reg22=reg42-reg22; reg7=reg93*reg7;
    reg24=reg24*reg94; reg25=reg68+reg25; reg12=reg14-reg12; reg1=reg2-reg1; reg67=reg57*reg67;
    reg22=reg7+reg22; reg8=reg8/reg10; reg0=reg16-reg0; reg12=reg11*reg12; reg2=elem.pos(2)[1]-elem.pos(0)[1];
    reg7=elem.pos(1)[1]-elem.pos(0)[1]; reg14=(*f.m).resolution*reg89; reg25=reg11*reg25; reg16=(*f.m).resolution*reg64; reg9=reg13-reg9;
    reg1=reg1/reg10; reg24=reg67-reg24; reg13=elem.pos(2)[0]-elem.pos(0)[0]; reg23=elem.pos(1)[0]-elem.pos(0)[0]; reg94=reg17*reg94;
    reg57=reg19*reg57; reg12=reg14+reg12; reg16=reg25+reg16; reg22=reg11*reg22; reg14=reg23*reg2;
    reg17=(*f.m).resolution*reg1; reg19=(*f.m).resolution*reg8; reg25=(*f.m).resolution*reg21; reg94=reg57-reg94; reg5=reg5/reg10;
    reg24=reg24/reg10; reg26=reg13*reg7; reg92=reg92*reg11; reg0=reg0/reg10; reg104=reg104*reg11;
    reg9=reg9/reg10; reg4=reg4/reg10; reg94=reg94/reg10; reg10=reg6/reg10; reg6=(*f.m).resolution*reg4;
    reg27=(*f.m).resolution*reg5; reg22=reg25+reg22; reg17=reg104-reg17; reg19=reg92+reg19; reg25=(*f.m).resolution*reg24;
    reg26=reg14-reg26; reg101=reg101*reg11; reg14=(*f.m).resolution*reg9; reg33=reg33*reg11; reg88=reg11*reg88;
    reg66=reg11*reg66; reg28=(*f.m).resolution*reg0; reg48=reg48*reg11; reg16=(*f.m).deltaT*reg16; reg12=(*f.m).deltaT*reg12;
    reg29=reg19*reg16; reg2=reg2/reg26; reg30=reg17*reg12; reg45=reg45*reg11; reg82=reg11*reg82;
    reg27=reg101+reg27; reg6=reg33-reg6; reg23=reg23/reg26; reg7=reg7/reg26; reg13=reg13/reg26;
    reg11=(*f.m).resolution*reg94; reg28=reg88-reg28; reg14=reg66+reg14; reg25=reg48+reg25; reg31=(*f.m).resolution*reg10;
    reg22=(*f.m).deltaT*reg22; reg32=reg7-reg2; reg33=reg29+reg30; reg11=reg82-reg11; reg35=reg28*reg16;
    reg31=reg45+reg31; reg36=reg25*reg22; reg38=reg14*reg12; reg39=reg6*reg12; reg42=reg13-reg23;
    reg43=reg27*reg16; reg44=0.5*reg13; reg45=reg36+reg33; reg46=0.5*reg2; reg47=0.5*reg7;
    reg48=reg11*reg22; reg51=0.5*reg23; reg56=reg35+reg38; reg57=reg31*reg22; reg58=reg43+reg39;
    reg60=0.5*reg42; reg62=0.5*reg32; reg63=reg19*reg7; reg65=reg25*reg51; reg66=reg56+reg48;
    reg67=reg32*reg19; reg68=2*reg45; reg69=reg25*reg44; reg70=reg60*reg25; reg71=reg25*reg62;
    reg72=reg25*reg46; reg73=reg47*reg25; reg74=reg2*reg19; reg75=reg17*reg13; reg76=reg17*reg23;
    reg77=reg42*reg17; reg78=reg57+reg58; reg67=reg70+reg67; reg70=reg6*reg23; reg79=1-var_inter[0];
    reg76=reg76-reg73; reg80=reg2*reg27; reg81=reg31*reg44; reg82=reg31*reg62; reg83=reg11*reg62;
    reg84=reg46*reg68; reg85=reg13*reg66; reg86=reg42*reg6; reg87=reg68*reg51; reg88=reg7*reg78;
    reg90=reg47*reg31; reg91=reg23*reg66; reg92=reg47*reg68; reg93=reg60*reg31; reg77=reg71+reg77;
    reg71=reg32*reg27; reg96=reg14*reg23; reg97=reg27*reg7; reg98=reg31*reg51; reg99=reg47*reg11;
    reg100=reg68*reg44; reg101=reg2*reg78; reg102=reg11*reg46; reg72=reg72-reg75; reg74=reg74-reg69;
    reg103=reg14*reg13; reg104=reg6*reg13; reg65=reg65-reg63; reg105=reg31*reg46; reg106=reg11*reg44;
    reg107=reg2*reg28; reg108=reg28*reg7; reg109=reg42*reg14; reg110=reg11*reg51; reg111=var_inter[1]*elem.f_vol_e[1];
    reg112=reg92-reg91; reg72=2*reg72; reg110=reg110-reg108; reg71=reg93+reg71; reg77=2*reg77;
    reg96=reg96-reg99; reg105=reg105-reg104; reg86=reg82+reg86; reg82=var_inter[0]*elem.f_vol_e[0]; reg93=var_inter[0]*elem.f_vol_e[1];
    reg74=2*reg74; reg67=2*reg67; reg65=2*reg65; reg70=reg70-reg90; reg80=reg80-reg81;
    reg113=reg42*reg66; reg107=reg107-reg106; reg114=reg68*reg62; reg98=reg98-reg97; reg76=2*reg76;
    reg109=reg83+reg109; reg83=reg100-reg101; reg115=reg85-reg84; reg102=reg102-reg103; reg116=var_inter[1]*elem.f_vol_e[0];
    reg79=reg79-var_inter[1]; reg117=reg32*reg78; reg118=reg88-reg87; reg119=reg60*reg68; reg120=reg60*reg65;
    reg121=reg60*reg74; reg122=reg32*reg98; reg123=reg32*reg105; reg124=reg32*reg80; reg125=reg32*reg71;
    reg126=reg60*reg72; reg127=reg79*elem.f_vol_e[0]; reg128=reg13*reg102; reg129=reg72*reg46; reg130=reg42*reg107;
    reg131=reg74*reg62; reg132=reg76*reg44; T reg133=reg2*reg70; T reg134=reg113+reg114; reg83=reg83-reg82;
    T reg135=reg47*reg76; T reg136=reg23*reg96; T reg137=reg42*reg102; T reg138=reg72*reg62; T reg139=reg42*reg110;
    T reg140=reg65*reg62; T reg141=reg42*reg96; T reg142=reg76*reg62; T reg143=reg105*reg2; T reg144=reg72*reg44;
    T reg145=reg65*reg51; T reg146=reg80*reg2; T reg147=reg74*reg44; T reg148=reg98*reg2; T reg149=reg65*reg44;
    T reg150=reg60*reg67; T reg151=reg32*reg86; T reg152=reg60*reg77; T reg153=elem.f_vol_e[1]*reg79; T reg154=reg117+reg119;
    reg112=reg112-reg111; T reg155=reg70*reg7; T reg156=reg76*reg51; reg118=reg118-reg116; reg98=reg98*reg7;
    T reg157=reg109*reg42; reg115=reg115-reg93; reg96=reg13*reg96; T reg158=reg76*reg46; reg70=reg32*reg70;
    reg76=reg60*reg76; T reg159=reg13*reg110; T reg160=reg65*reg46; T reg161=reg77*reg62; reg96=reg158-reg96;
    reg98=reg145-reg98; reg132=reg133-reg132; reg159=reg160-reg159; reg144=reg143-reg144; reg155=reg156-reg155;
    reg128=reg129-reg128; reg149=reg148-reg149; reg135=reg136-reg135; reg129=reg127+reg154; reg133=reg153+reg134;
    reg83=reg26*reg83; reg138=reg137+reg138; reg140=reg139+reg140; reg142=reg141+reg142; reg147=reg146-reg147;
    reg131=reg130+reg131; reg161=reg157+reg161; reg76=reg70+reg76; reg152=reg151+reg152; reg115=reg26*reg115;
    reg118=reg26*reg118; reg112=reg26*reg112; reg121=reg124+reg121; reg150=reg125+reg150; reg126=reg123+reg126;
    reg120=reg122+reg120; reg128=reg128*reg26; reg161=reg161*reg26; reg132=reg132*reg26; reg131=reg131*reg26;
    reg149=reg149*reg26; reg147=reg26*reg147; reg120=reg120*reg26; reg144=reg144*reg26; reg142=reg26*reg142;
    reg140=reg26*reg140; reg138=reg26*reg138; reg126=reg126*reg26; reg70=reg26*reg133; reg122=reg26*reg129;
    reg121=reg121*reg26; reg135=reg26*reg135; reg112=ponderation*reg112; reg155=reg155*reg26; reg118=ponderation*reg118;
    reg98=reg98*reg26; reg115=ponderation*reg115; reg83=ponderation*reg83; reg96=reg96*reg26; reg152=reg152*reg26;
    reg76=reg76*reg26; reg150=reg150*reg26; reg159=reg159*reg26; T tmp_0_1=ponderation*reg152; T tmp_0_5=ponderation*reg76;
    T tmp_5_5=ponderation*reg135; T vec_4=-reg118; T tmp_3_3=ponderation*reg128; reg76=ponderation*reg122; T vec_0=reg76;
    T tmp_0_2=ponderation*reg121; T tmp_4_4=ponderation*reg98; reg98=ponderation*reg70; T vec_1=reg98; T tmp_0_4=ponderation*reg120;
    T tmp_4_5=ponderation*reg155; T vec_3=-reg115; T tmp_1_3=ponderation*reg138; T tmp_0_3=ponderation*reg126; T tmp_2_5=ponderation*reg132;
    T tmp_1_4=ponderation*reg140; T vec_5=-reg112; T vec_2=-reg83; T tmp_1_5=ponderation*reg142; T tmp_3_4=ponderation*reg159;
    T tmp_2_3=ponderation*reg144; T tmp_3_5=ponderation*reg96; T tmp_1_1=ponderation*reg161; T tmp_2_2=ponderation*reg147; T tmp_0_0=ponderation*reg150;
    T tmp_2_4=ponderation*reg149; T tmp_1_2=ponderation*reg131;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_3;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_4;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_5;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_3;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_4;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_5;
    matrix(indices[1]+0,indices[1]+0) += tmp_2_2;
    matrix(indices[1]+0,indices[1]+1) += tmp_2_3;
    matrix(indices[1]+0,indices[2]+0) += tmp_2_4;
    matrix(indices[1]+0,indices[2]+1) += tmp_2_5;
    matrix(indices[1]+1,indices[1]+1) += tmp_3_3;
    matrix(indices[1]+1,indices[2]+0) += tmp_3_4;
    matrix(indices[1]+1,indices[2]+1) += tmp_3_5;
    matrix(indices[2]+0,indices[2]+0) += tmp_4_4;
    matrix(indices[2]+0,indices[2]+1) += tmp_4_5;
    matrix(indices[2]+1,indices[2]+1) += tmp_5_5;
sollicitation[indices[0]+0] += vec_0;
sollicitation[indices[0]+1] += vec_1;
sollicitation[indices[1]+0] += vec_2;
sollicitation[indices[1]+1] += vec_3;
sollicitation[indices[2]+0] += vec_4;
sollicitation[indices[2]+1] += vec_5;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=pow((*f.m).v1[0],2); T reg1=pow((*f.m).v1[1],2); reg1=reg0+reg1; reg0=pow((*f.m).v1[2],2); T reg2=pow((*f.m).v2[0],2);
    T reg3=pow((*f.m).v2[1],2); T reg4=pow((*f.m).v2[2],2); reg0=reg1+reg0; reg3=reg2+reg3; reg1=2*(*f.m).shear_modulus_23;
    reg2=2*(*f.m).shear_modulus_13; reg0=pow(reg0,0.5); T reg5=2*(*f.m).shear_modulus_12; reg4=reg3+reg4; reg1=1.0/reg1;
    reg2=1.0/reg2; reg5=1.0/reg5; reg3=reg2*reg1; reg4=pow(reg4,0.5); T reg6=(*f.m).v1[0]/reg0;
    T reg7=(*f.m).v1[1]/reg0; T reg8=(*f.m).v2[0]/reg4; T reg9=(*f.m).v2[1]/reg4; reg0=(*f.m).v1[2]/reg0; T reg10=2*reg6;
    T reg11=2*reg7; T reg12=reg5*reg3; T reg13=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg14=1.0/(*f.m).elastic_modulus_3; T reg15=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg16=reg11*reg9; T reg17=reg13*reg12; reg4=(*f.m).v2[2]/reg4; T reg18=reg10*reg8; T reg19=reg14*reg12;
    T reg20=pow(reg8,2); T reg21=2*reg0; T reg22=1.0/(*f.m).elastic_modulus_1; T reg23=reg15*reg12; T reg24=1.0/(*f.m).elastic_modulus_2;
    T reg25=pow(reg9,2); T reg26=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg27=reg22*reg18; T reg28=reg20*reg26; T reg29=reg25*reg24;
    T reg30=reg26*reg18; T reg31=reg25*reg26; T reg32=reg24*reg16; T reg33=reg26*reg16; T reg34=reg13*reg23;
    T reg35=pow(reg4,2); T reg36=reg20*reg22; reg21=reg4*reg21; T reg37=reg24*reg19; T reg38=reg26*reg19;
    T reg39=reg13*reg17; T reg40=reg35*reg13; reg28=reg29-reg28; reg31=reg36-reg31; reg29=reg35*reg15;
    reg36=reg15*reg21; reg33=reg27-reg33; reg27=reg20*reg15; T reg41=reg25*reg13; reg30=reg32-reg30;
    reg39=reg37-reg39; reg34=reg38+reg34; reg32=reg26*reg17; reg37=pow(reg6,2); T reg42=pow(reg7,2);
    T reg43=reg13*reg21; T reg44=reg15*reg18; T reg45=reg24*reg23; T reg46=reg13*reg16; T reg47=reg37*reg26;
    reg43=reg30-reg43; reg36=reg33-reg36; reg30=reg22*reg39; reg33=reg26*reg34; T reg48=reg42*reg24;
    T reg49=reg6*reg8; T reg50=reg7*reg9; T reg51=reg37*reg22; T reg52=reg42*reg26; reg40=reg28-reg40;
    reg28=reg32+reg45; reg46=reg44+reg46; reg21=reg14*reg21; reg44=pow(reg0,2); reg41=reg27+reg41;
    reg27=reg35*reg14; reg29=reg31-reg29; reg31=reg7*reg8; T reg53=reg6*reg9; T reg54=reg15*reg23;
    T reg55=reg7*reg4; T reg56=reg0*reg4; T reg57=reg31+reg53; T reg58=reg13*reg3; T reg59=reg0*reg9;
    T reg60=reg20*reg29; reg19=reg22*reg19; T reg61=reg15*reg28; reg52=reg51-reg52; reg51=reg44*reg15;
    T reg62=reg25*reg40; T reg63=reg50*reg40; T reg64=reg49*reg29; reg46=reg21-reg46; reg21=reg25*reg43;
    T reg65=reg20*reg36; T reg66=reg24*reg12; T reg67=reg6*reg4; T reg68=reg15*reg17; T reg69=reg5*reg1;
    T reg70=reg0*reg8; T reg71=reg15*reg3; T reg72=reg42*reg43; T reg73=reg37*reg36; T reg74=reg42*reg13;
    T reg75=reg8*reg9; T reg76=reg44*reg13; reg33=reg30-reg33; reg29=reg37*reg29; reg40=reg42*reg40;
    reg30=reg37*reg15; reg47=reg48-reg47; reg12=reg26*reg12; reg36=reg49*reg36; reg43=reg50*reg43;
    reg3=reg14*reg3; reg48=2*reg8; reg41=reg27-reg41; reg68=reg38+reg68; reg27=reg10*reg7;
    reg74=reg30+reg74; reg30=reg71*reg13; reg38=reg58*reg13; T reg77=reg3*reg26; T reg78=reg44*reg14;
    reg51=reg52-reg51; reg52=reg5*reg2; T reg79=reg15*reg66; reg3=reg3*reg24; reg23=reg26*reg23;
    T reg80=reg55-reg59; T reg81=reg69*reg15; T reg82=reg44*reg46; reg72=reg73+reg72; reg73=reg8*reg4;
    T reg83=reg56*reg46; reg43=reg36+reg43; reg36=reg57*reg5; T reg84=reg48*reg9; T reg85=reg75*reg5;
    reg54=reg19-reg54; reg17=reg22*reg17; reg19=reg67+reg70; reg40=reg29+reg40; reg29=reg44*reg41;
    reg62=reg60+reg62; reg60=reg69*reg13; reg69=reg69*reg14; T reg86=reg12*reg15; reg61=reg33-reg61;
    reg21=reg65+reg21; reg46=reg35*reg46; reg63=reg64+reg63; reg76=reg47-reg76; reg33=reg56*reg41;
    reg41=reg35*reg41; reg74=reg78-reg74; reg54=reg54/reg61; reg33=reg63+reg33; reg23=reg17+reg23;
    reg67=reg70-reg67; reg47=reg57*reg36; reg83=reg43+reg83; reg86=reg17+reg86; reg66=reg22*reg66;
    reg39=reg39/reg61; reg17=reg57*reg85; reg34=reg34/reg61; reg79=reg32+reg79; reg41=reg62+reg41;
    reg43=reg27*reg85; reg29=reg40+reg29; reg55=reg59+reg55; reg14=reg52*reg14; reg40=2*reg80;
    reg59=reg9*reg4; reg62=2*reg9; reg63=reg48*reg4; reg64=reg81*reg13; reg65=reg42*reg76;
    reg70=reg37*reg51; reg78=reg25*reg76; reg82=reg72+reg82; reg72=reg27*reg36; T reg87=reg6*reg7;
    T reg88=reg20*reg51; T reg89=reg73*reg2; T reg90=reg19*reg2; reg12=reg12*reg26; reg58=reg58*reg26;
    T reg91=reg52*reg13; reg71=reg71*reg24; reg52=reg52*reg15; reg38=reg3-reg38; reg3=reg10*reg0;
    reg30=reg77+reg30; reg77=reg69*reg24; reg68=reg68/reg61; reg69=reg69*reg26; T reg92=reg60*reg13;
    reg46=reg21+reg46; reg36=reg84*reg36; reg85=reg84*reg85; reg69=reg64+reg69; reg21=reg91*reg13;
    reg64=reg6*reg0; T reg93=reg14*reg26; reg14=reg14*reg24; reg13=reg52*reg13; reg28=reg28/reg61;
    reg43=reg29+reg43; reg29=reg40*reg67; T reg94=pow(reg67,2); T reg95=pow(reg80,2); T reg96=reg27*reg68;
    T reg97=reg44*reg74; reg65=reg70+reg65; reg70=reg55*reg1; reg72=reg82+reg72; reg82=reg3*reg90;
    T reg98=reg59*reg1; reg78=reg88+reg78; reg88=reg35*reg74; T reg99=reg27*reg39; T reg100=reg84*reg34;
    T reg101=reg25*reg54; T reg102=reg42*reg68; reg85=reg41+reg85; reg41=reg63*reg89; T reg103=reg37*reg68;
    T reg104=reg20*reg54; reg5=reg87*reg5; T reg105=reg11*reg0; T reg106=reg62*reg4; T reg107=reg19*reg90;
    reg47=reg83+reg47; reg83=reg19*reg89; reg17=reg33+reg17; reg33=reg37*reg39; T reg108=reg20*reg34;
    reg86=reg86/reg61; reg23=reg23/reg61; reg76=reg50*reg76; reg51=reg49*reg51; reg12=reg66-reg12;
    reg60=reg60*reg26; reg81=reg81*reg24; reg38=reg38*reg22; reg79=reg79/reg61; reg30=reg30*reg26;
    reg66=reg58+reg71; T reg109=reg42*reg39; reg92=reg77-reg92; reg36=reg46+reg36; reg90=reg63*reg90;
    reg46=reg25*reg34; reg77=reg84*reg54; reg89=reg3*reg89; reg77=reg96+reg77; reg101=reg102+reg101;
    reg96=reg94*reg23; reg102=reg95*reg23; T reg110=reg84*reg86; T reg111=reg27*reg79; T reg112=reg25*reg86;
    T reg113=reg42*reg79; T reg114=reg20*reg86; T reg115=reg37*reg79; T reg116=reg29*reg23; reg2=reg64*reg2;
    reg97=reg65+reg97; reg65=reg27*reg5; T reg117=reg106*reg70; reg90=reg36+reg90; reg92=reg92*reg22;
    reg66=reg66*reg15; reg30=reg38-reg30; reg24=reg52*reg24; reg91=reg91*reg26; reg12=reg12/reg61;
    reg36=reg105*reg98; reg38=reg55*reg70; reg52=reg95*reg28; reg108=reg33+reg108; reg107=reg47+reg107;
    reg33=reg55*reg98; reg83=reg17+reg83; reg82=reg72+reg82; reg74=reg56*reg74; reg76=reg51+reg76;
    reg17=reg7*reg0; reg89=reg43+reg89; reg69=reg69*reg26; reg98=reg106*reg98; reg43=reg60+reg81;
    reg21=reg14-reg21; reg41=reg85+reg41; reg104=reg103+reg104; reg70=reg105*reg70; reg13=reg93+reg13;
    reg88=reg78+reg88; reg14=reg29*reg28; reg100=reg99+reg100; reg47=reg84*reg5; reg46=reg109+reg46;
    reg51=reg94*reg28; reg38=reg107+reg38; reg70=reg82+reg70; reg33=reg83+reg33; reg47=reg88+reg47;
    reg66=reg30-reg66; reg1=reg17*reg1; reg5=reg57*reg5; reg74=reg76+reg74; reg30=reg63*reg2;
    reg116=reg77+reg116; reg29=reg29*reg12; reg110=reg111+reg110; reg114=reg115+reg114; reg98=reg41+reg98;
    reg41=reg95*reg12; reg72=reg94*reg12; reg112=reg113+reg112; reg102=reg104+reg102; reg14=reg100+reg14;
    reg51=reg46+reg51; reg96=reg101+reg96; reg46=reg91+reg24; reg26=reg13*reg26; reg117=reg90+reg117;
    reg22=reg21*reg22; reg43=reg43*reg15; reg36=reg89+reg36; reg52=reg108+reg52; reg13=reg3*reg2;
    reg65=reg97+reg65; reg69=reg92-reg69; reg41=reg114+reg41; reg21=reg96*reg75; reg76=reg51*reg87;
    reg15=reg46*reg15; reg46=reg106*reg1; reg77=reg33*reg70; reg30=reg47+reg30; reg47=reg7*reg80;
    reg78=reg6*reg67; reg66=reg66/reg61; reg82=reg102*reg75; reg83=reg87*reg52; reg2=reg19*reg2;
    reg5=reg74+reg5; reg74=reg105*reg1; reg85=reg38*reg98; reg88=reg38*reg36; reg29=reg110+reg29;
    reg89=reg67*reg80; reg72=reg112+reg72; reg90=reg116*reg75; reg92=reg14*reg87; reg43=reg69-reg43;
    reg13=reg65+reg13; reg65=reg33*reg117; reg26=reg22-reg26; reg22=reg20*reg102; reg69=reg37*reg52;
    reg93=reg70*reg98; reg74=reg13+reg74; reg46=reg30+reg46; reg13=reg49*reg66; reg30=reg50*reg66;
    reg97=reg57*reg66; reg99=reg8*reg67; reg100=reg9*reg80; reg78=reg47+reg78; reg47=reg7*reg67;
    reg101=reg6*reg80; reg103=reg25*reg116; reg104=reg42*reg14; reg107=reg25*reg96; reg108=reg42*reg51;
    reg102=reg25*reg102; reg52=reg42*reg52; reg109=reg72*reg89; reg21=reg76+reg21; reg90=reg92+reg90;
    reg76=reg29*reg89; reg2=reg5+reg2; reg1=reg55*reg1; reg5=reg41*reg89; reg82=reg83+reg82;
    reg116=reg20*reg116; reg14=reg37*reg14; reg65=reg85-reg65; reg77=reg88-reg77; reg83=reg117*reg36;
    reg51=reg37*reg51; reg15=reg26-reg15; reg43=reg43/reg61; reg96=reg20*reg96; reg26=elem.pos(1)[1]-elem.pos(0)[1];
    reg85=reg57*reg30; reg88=elem.pos(2)[1]-elem.pos(0)[1]; reg76=reg90+reg76; reg90=reg57*reg97; reg61=reg15/reg61;
    reg102=reg52+reg102; reg15=reg94*reg41; reg1=reg2+reg1; reg101=reg43*reg101; reg47=reg43*reg47;
    reg2=reg43*reg78; reg99=reg100+reg99; reg52=reg9*reg67; reg92=reg8*reg80; reg107=reg108+reg107;
    reg100=reg94*reg72; reg108=reg74*reg65; reg103=reg104+reg103; reg104=reg94*reg29; reg72=reg95*reg72;
    reg96=reg51+reg96; reg22=reg69+reg22; reg41=reg95*reg41; reg116=reg14+reg116; reg29=reg95*reg29;
    reg14=elem.pos(2)[0]-elem.pos(0)[0]; reg51=elem.pos(1)[0]-elem.pos(0)[0]; reg5=reg82+reg5; reg69=reg57*reg13; reg93=reg83-reg93;
    reg82=reg77*reg46; reg109=reg21+reg109; reg21=reg36*reg46; reg83=reg74*reg98; reg110=reg97*reg16;
    reg104=reg103+reg104; reg103=reg30*reg18; reg72=reg96+reg72; reg96=reg33*reg46; reg30=reg30*reg16;
    reg100=reg107+reg100; reg82=reg108-reg82; reg107=reg1*reg93; reg108=reg14*reg26; reg36=reg1*reg36;
    reg111=reg13*reg16; reg15=reg102+reg15; reg41=reg22+reg41; reg13=reg13*reg18; reg10=reg10*reg80;
    reg22=reg61*reg99; reg33=reg33*reg74; reg52=reg61*reg52; reg92=reg61*reg92; reg29=reg116+reg29;
    reg11=reg11*reg67; reg97=reg97*reg18; reg102=reg51*reg88; reg112=reg2*reg78; reg90=reg76+reg90;
    reg69=reg5+reg69; reg5=reg47*reg78; reg85=reg109+reg85; reg98=reg1*reg98; reg76=reg101*reg78;
    reg109=reg38*reg46; reg111=reg15+reg111; reg15=reg101*reg11; reg107=reg82+reg107; reg48=reg48*reg80;
    reg30=reg100+reg30; reg82=reg47*reg11; reg100=reg1*reg117; reg1=reg1*reg70; reg98=reg96-reg98;
    reg62=reg62*reg67; reg110=reg104+reg110; reg96=reg2*reg11; reg38=reg38*reg74; reg36=reg33-reg36;
    reg74=reg117*reg74; reg46=reg70*reg46; reg21=reg83-reg21; reg33=reg22*reg99; reg112=reg90+reg112;
    reg70=reg52*reg99; reg5=reg85+reg5; reg83=reg92*reg99; reg76=reg69+reg76; reg2=reg2*reg10;
    reg97=reg29+reg97; reg108=reg102-reg108; reg47=reg47*reg10; reg103=reg72+reg103; reg101=reg101*reg10;
    reg13=reg41+reg13; reg29=reg22*reg62; reg96=reg110+reg96; reg88=reg88/reg108; reg14=reg14/reg108;
    reg26=reg26/reg108; reg51=reg51/reg108; reg33=reg112+reg33; reg21=reg21/reg107; reg41=reg52*reg62;
    reg82=reg30+reg82; reg70=reg5+reg70; reg46=reg74-reg46; reg83=reg76+reg83; reg36=reg36/reg107;
    reg22=reg22*reg48; reg100=reg109-reg100; reg101=reg13+reg101; reg5=reg92*reg48; reg98=reg98/reg107;
    reg47=reg103+reg47; reg52=reg52*reg48; reg13=1-(*f.m).resolution; reg1=reg38-reg1; reg15=reg111+reg15;
    reg92=reg92*reg62; reg2=reg97+reg2; reg1=reg1/reg107; reg46=reg46/reg107; reg52=reg47+reg52;
    reg93=reg93/reg107; reg30=reg26-reg88; reg77=reg77/reg107; reg5=reg101+reg5; reg22=reg2+reg22;
    reg2=reg14-reg51; reg65=reg65/reg107; reg107=reg100/reg107; reg38=reg33*reg13; reg47=(*f.m).resolution*reg98;
    reg69=(*f.m).resolution*reg36; reg72=(*f.m).resolution*reg21; reg74=reg70*reg13; reg76=reg83*reg13; reg92=reg15+reg92;
    reg41=reg82+reg41; reg29=reg96+reg29; reg15=0.5*reg88; reg82=(*f.m).resolution*reg46; reg85=0.5*reg26;
    reg90=0.5*reg30; reg96=(*f.m).resolution*reg1; reg97=(*f.m).resolution*reg107; reg100=0.5*reg51; reg101=reg5*reg13;
    reg102=reg52*reg13; reg103=(*f.m).resolution*reg77; reg104=(*f.m).resolution*reg65; reg109=0.5*reg14; reg110=(*f.m).resolution*reg93;
    reg47=reg76+reg47; reg69=reg74-reg69; reg74=reg22*reg13; reg76=0.5*reg2; reg111=reg13*reg92;
    reg112=reg13*reg41; reg113=reg13*reg29; reg72=reg38+reg72; reg38=reg88*reg47; reg114=reg72*reg109;
    reg115=reg69*reg51; reg116=reg85*reg72; reg104=reg101+reg104; reg101=reg76*reg72; reg117=reg30*reg47;
    reg82=reg113-reg82; reg113=reg69*reg14; reg103=reg102-reg103; reg110=reg74+reg110; reg74=reg2*reg69;
    reg102=reg47*reg26; reg97=reg111-reg97; reg96=reg112+reg96; reg111=reg72*reg100; reg112=reg72*reg15;
    T reg118=reg72*reg90; T reg119=reg103*reg14; T reg120=reg110*reg15; T reg121=reg82*reg90; reg112=reg112-reg113;
    reg38=reg38-reg114; T reg122=reg85*reg82; T reg123=reg96*reg51; T reg124=reg88*reg97; T reg125=reg82*reg109;
    T reg126=reg96*reg14; T reg127=reg30*reg97; T reg128=reg76*reg82; T reg129=reg82*reg100; T reg130=reg82*reg15;
    T reg131=reg97*reg26; T reg132=reg30*reg104; T reg133=reg76*reg110; T reg134=reg2*reg96; reg117=reg101+reg117;
    reg101=reg85*reg110; T reg135=reg110*reg90; T reg136=reg2*reg103; reg74=reg118+reg74; reg111=reg111-reg102;
    reg118=reg110*reg109; T reg137=reg110*reg100; T reg138=reg88*reg104; T reg139=reg104*reg26; T reg140=reg103*reg51;
    reg115=reg115-reg116; reg115=2*reg115; reg132=reg133+reg132; reg123=reg123-reg122; reg124=reg124-reg125;
    reg140=reg140-reg101; reg117=2*reg117; reg128=reg127+reg128; reg136=reg135+reg136; reg138=reg138-reg118;
    reg130=reg130-reg126; reg74=2*reg74; reg129=reg129-reg131; reg134=reg121+reg134; reg38=2*reg38;
    reg111=2*reg111; reg112=2*reg112; reg120=reg120-reg119; reg137=reg137-reg139; reg121=reg38*reg15;
    reg127=reg14*reg123; reg133=reg85*reg112; reg135=reg115*reg109; T reg141=reg88*reg140; T reg142=reg111*reg109;
    T reg143=reg137*reg88; T reg144=reg115*reg15; T reg145=reg124*reg14; T reg146=reg112*reg109; T reg147=reg120*reg88;
    T reg148=reg14*reg129; T reg149=reg111*reg15; T reg150=reg14*reg130; T reg151=reg112*reg15; T reg152=reg38*reg109;
    T reg153=reg138*reg88; T reg154=reg136*reg88; T reg155=reg117*reg109; T reg156=reg132*reg88; T reg157=reg132*reg26;
    T reg158=reg74*reg109; T reg159=reg74*reg15; T reg160=reg117*reg100; T reg161=reg74*reg100; T reg162=reg136*reg26;
    T reg163=reg38*reg100; T reg164=reg138*reg26; T reg165=reg112*reg100; T reg166=reg120*reg26; T reg167=reg128*reg14;
    T reg168=reg85*reg74; T reg169=reg134*reg14; T reg170=reg137*reg26; T reg171=reg115*reg100; T reg172=reg140*reg26;
    T reg173=reg128*reg51; T reg174=reg85*reg117; T reg175=reg134*reg51; T reg176=reg124*reg51; T reg177=reg85*reg38;
    reg128=reg2*reg128; T reg178=reg117*reg15; reg134=reg134*reg2; T reg179=reg51*reg130; T reg180=reg111*reg90;
    T reg181=reg2*reg129; T reg182=reg74*reg90; T reg183=reg112*reg90; reg130=reg2*reg130; T reg184=reg76*reg111;
    reg137=reg30*reg137; reg112=reg76*reg112; T reg185=reg85*reg115; T reg186=reg51*reg123; reg120=reg30*reg120;
    T reg187=reg76*reg38; reg138=reg30*reg138; T reg188=reg85*reg111; reg129=reg51*reg129; reg124=reg2*reg124;
    reg38=reg38*reg90; reg132=reg30*reg132; T reg189=reg76*reg117; reg117=reg117*reg90; reg123=reg2*reg123;
    reg136=reg30*reg136; T reg190=reg115*reg90; reg111=reg111*reg100; reg140=reg30*reg140; reg115=reg76*reg115;
    reg74=reg76*reg74; reg169=reg159-reg169; reg162=reg161-reg162; reg182=reg134+reg182; reg157=reg160-reg157;
    reg115=reg140+reg115; reg38=reg124+reg38; reg164=reg163-reg164; reg166=reg165-reg166; reg133=reg179-reg133;
    reg170=reg111-reg170; reg187=reg138+reg187; reg127=reg144-reg127; reg146=reg147-reg146; reg148=reg149-reg148;
    reg142=reg143-reg142; reg74=reg136+reg74; reg112=reg120+reg112; reg150=reg151-reg150; reg145=reg121-reg145;
    reg135=reg141-reg135; reg189=reg132+reg189; reg184=reg137+reg184; reg167=reg178-reg167; reg117=reg128+reg117;
    reg174=reg173-reg174; reg152=reg153-reg152; reg158=reg154-reg158; reg172=reg171-reg172; reg185=reg186-reg185;
    reg177=reg176-reg177; reg190=reg123+reg190; reg183=reg130+reg183; reg188=reg129-reg188; reg155=reg156-reg155;
    reg168=reg175-reg168; reg180=reg181+reg180; reg148=reg148*reg108; reg188=reg108*reg188; reg142=reg142*reg108;
    reg150=reg150*reg108; reg177=reg177*reg108; reg112=reg112*reg108; reg135=reg135*reg108; reg145=reg145*reg108;
    reg152=reg108*reg152; reg185=reg108*reg185; reg189=reg189*reg108; reg169=reg169*reg108; reg184=reg184*reg108;
    reg183=reg108*reg183; reg180=reg108*reg180; reg117=reg117*reg108; reg167=reg167*reg108; reg170=reg170*reg108;
    reg166=reg166*reg108; reg133=reg108*reg133; reg155=reg108*reg155; reg164=reg164*reg108; reg190=reg108*reg190;
    reg38=reg38*reg108; reg162=reg162*reg108; reg172=reg172*reg108; reg146=reg146*reg108; reg127=reg127*reg108;
    reg158=reg108*reg158; reg187=reg187*reg108; reg168=reg168*reg108; reg74=reg74*reg108; reg174=reg174*reg108;
    reg182=reg182*reg108; reg115=reg115*reg108; reg157=reg157*reg108; T tmp_0_0=ponderation*reg189; T tmp_2_5=ponderation*reg135;
    T tmp_2_0=ponderation*reg155; T tmp_2_1=ponderation*reg158; T tmp_2_2=ponderation*reg152; T tmp_0_1=ponderation*reg74; T tmp_1_4=ponderation*reg180;
    T tmp_2_4=ponderation*reg142; T tmp_2_3=ponderation*reg146; T tmp_1_5=ponderation*reg190; T tmp_4_4=ponderation*reg170; T tmp_4_3=ponderation*reg166;
    T tmp_1_2=ponderation*reg38; T tmp_0_5=ponderation*reg115; T tmp_4_2=ponderation*reg164; T tmp_5_3=ponderation*reg133; T tmp_4_1=ponderation*reg162;
    T tmp_4_5=ponderation*reg172; T tmp_4_0=ponderation*reg157; T tmp_1_1=ponderation*reg182; T tmp_5_0=ponderation*reg174; T tmp_3_5=ponderation*reg127;
    T tmp_5_1=ponderation*reg168; T tmp_0_2=ponderation*reg187; T tmp_3_4=ponderation*reg148; T tmp_5_4=ponderation*reg188; T tmp_3_3=ponderation*reg150;
    T tmp_0_3=ponderation*reg112; T tmp_3_2=ponderation*reg145; T tmp_5_2=ponderation*reg177; T tmp_5_5=ponderation*reg185; T tmp_3_1=ponderation*reg169;
    T tmp_0_4=ponderation*reg184; T tmp_1_0=ponderation*reg117; T tmp_3_0=ponderation*reg167; T tmp_1_3=ponderation*reg183;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_3;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_4;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_5;
    matrix(indices[0]+1,indices[0]+0) += tmp_1_0;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_3;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_4;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_5;
    matrix(indices[1]+0,indices[0]+0) += tmp_2_0;
    matrix(indices[1]+0,indices[0]+1) += tmp_2_1;
    matrix(indices[1]+0,indices[1]+0) += tmp_2_2;
    matrix(indices[1]+0,indices[1]+1) += tmp_2_3;
    matrix(indices[1]+0,indices[2]+0) += tmp_2_4;
    matrix(indices[1]+0,indices[2]+1) += tmp_2_5;
    matrix(indices[1]+1,indices[0]+0) += tmp_3_0;
    matrix(indices[1]+1,indices[0]+1) += tmp_3_1;
    matrix(indices[1]+1,indices[1]+0) += tmp_3_2;
    matrix(indices[1]+1,indices[1]+1) += tmp_3_3;
    matrix(indices[1]+1,indices[2]+0) += tmp_3_4;
    matrix(indices[1]+1,indices[2]+1) += tmp_3_5;
    matrix(indices[2]+0,indices[0]+0) += tmp_4_0;
    matrix(indices[2]+0,indices[0]+1) += tmp_4_1;
    matrix(indices[2]+0,indices[1]+0) += tmp_4_2;
    matrix(indices[2]+0,indices[1]+1) += tmp_4_3;
    matrix(indices[2]+0,indices[2]+0) += tmp_4_4;
    matrix(indices[2]+0,indices[2]+1) += tmp_4_5;
    matrix(indices[2]+1,indices[0]+0) += tmp_5_0;
    matrix(indices[2]+1,indices[0]+1) += tmp_5_1;
    matrix(indices[2]+1,indices[1]+0) += tmp_5_2;
    matrix(indices[2]+1,indices[1]+1) += tmp_5_3;
    matrix(indices[2]+1,indices[2]+0) += tmp_5_4;
    matrix(indices[2]+1,indices[2]+1) += tmp_5_5;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=pow((*f.m).v1[0],2); T reg1=pow((*f.m).v1[1],2); reg1=reg0+reg1; reg0=pow((*f.m).v1[2],2); T reg2=pow((*f.m).v2[0],2);
    T reg3=pow((*f.m).v2[1],2); reg0=reg1+reg0; reg1=pow((*f.m).v2[2],2); reg3=reg2+reg3; reg2=2*(*f.m).shear_modulus_23;
    T reg4=2*(*f.m).shear_modulus_13; reg2=1.0/reg2; reg4=1.0/reg4; reg0=pow(reg0,0.5); reg1=reg3+reg1;
    reg3=2*(*f.m).shear_modulus_12; reg1=pow(reg1,0.5); T reg5=(*f.m).v1[1]/reg0; T reg6=(*f.m).v1[0]/reg0; T reg7=reg4*reg2;
    reg3=1.0/reg3; reg0=(*f.m).v1[2]/reg0; T reg8=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg11=reg3*reg7; T reg12=(*f.m).v2[0]/reg1; T reg13=2*reg5; T reg14=2*reg6; T reg15=(*f.m).v2[1]/reg1;
    T reg16=pow(reg12,2); T reg17=pow(reg15,2); T reg18=reg14*reg12; T reg19=reg13*reg15; T reg20=2*reg0;
    T reg21=reg8*reg11; reg1=(*f.m).v2[2]/reg1; T reg22=1.0/(*f.m).elastic_modulus_1; T reg23=1.0/(*f.m).elastic_modulus_2; T reg24=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1;
    T reg25=reg10*reg11; T reg26=reg9*reg11; reg20=reg1*reg20; T reg27=pow(reg1,2); T reg28=reg16*reg22;
    T reg29=reg17*reg24; T reg30=reg22*reg18; T reg31=reg23*reg26; T reg32=reg24*reg19; T reg33=reg16*reg24;
    T reg34=reg17*reg23; T reg35=reg24*reg18; T reg36=reg23*reg19; T reg37=reg10*reg25; T reg38=reg24*reg26;
    T reg39=reg10*reg21; reg33=reg34-reg33; reg34=reg27*reg10; T reg40=reg10*reg20; T reg41=reg16*reg8;
    T reg42=reg17*reg10; T reg43=reg23*reg21; reg35=reg36-reg35; reg36=reg24*reg25; reg32=reg30-reg32;
    reg39=reg38+reg39; reg30=reg10*reg19; T reg44=reg8*reg18; T reg45=pow(reg5,2); T reg46=reg8*reg20;
    reg37=reg31-reg37; reg29=reg28-reg29; reg28=reg27*reg8; reg31=pow(reg6,2); T reg47=reg24*reg39;
    T reg48=reg36+reg43; T reg49=reg27*reg9; reg20=reg9*reg20; reg42=reg41+reg42; reg30=reg44+reg30;
    reg40=reg35-reg40; reg35=reg6*reg15; reg41=reg22*reg37; reg44=reg5*reg12; T reg50=reg5*reg15;
    T reg51=reg6*reg12; reg28=reg29-reg28; reg46=reg32-reg46; reg29=reg31*reg22; reg32=reg31*reg24;
    T reg52=pow(reg0,2); T reg53=reg45*reg23; T reg54=reg45*reg24; reg34=reg33-reg34; reg33=reg50*reg40;
    T reg55=reg51*reg46; T reg56=reg24*reg11; T reg57=reg9*reg7; T reg58=reg17*reg34; T reg59=reg51*reg28;
    T reg60=reg50*reg34; reg47=reg41-reg47; reg41=reg8*reg7; reg7=reg10*reg7; T reg61=reg8*reg48;
    T reg62=reg3*reg2; T reg63=reg16*reg46; T reg64=reg17*reg40; reg40=reg45*reg40; reg46=reg31*reg46;
    T reg65=reg16*reg28; reg11=reg23*reg11; T reg66=reg31*reg8; T reg67=2*reg12; T reg68=reg52*reg8;
    reg54=reg29-reg54; reg29=reg0*reg1; T reg69=reg44+reg35; T reg70=reg12*reg15; T reg71=reg45*reg10;
    T reg72=reg52*reg10; reg32=reg53-reg32; reg53=reg5*reg1; T reg73=reg0*reg15; T reg74=reg8*reg25;
    T reg75=reg6*reg1; T reg76=reg0*reg12; reg28=reg31*reg28; reg26=reg22*reg26; reg34=reg45*reg34;
    reg42=reg49-reg42; reg49=reg8*reg21; reg30=reg20-reg30; reg20=reg57*reg23; reg49=reg26-reg49;
    reg68=reg54-reg68; reg26=reg29*reg42; reg33=reg55+reg33; reg54=reg29*reg30; reg55=reg62*reg8;
    T reg77=reg3*reg4; T reg78=reg8*reg11; reg21=reg24*reg21; T reg79=reg62*reg10; reg74=reg38+reg74;
    reg38=reg52*reg30; reg40=reg46+reg40; reg72=reg32-reg72; reg32=reg69*reg3; reg46=reg56*reg8;
    T reg80=reg14*reg5; reg25=reg22*reg25; T reg81=reg53-reg73; reg71=reg66+reg71; reg66=reg52*reg42;
    reg34=reg28+reg34; reg61=reg47-reg61; reg28=reg52*reg9; reg42=reg27*reg42; reg58=reg65+reg58;
    reg47=reg75+reg76; reg65=reg12*reg1; reg30=reg27*reg30; reg64=reg63+reg64; reg63=reg67*reg15;
    reg62=reg62*reg9; reg60=reg59+reg60; reg59=reg70*reg3; T reg82=reg7*reg10; reg57=reg57*reg24;
    T reg83=reg41*reg10; T reg84=reg47*reg4; T reg85=reg15*reg1; reg53=reg73+reg53; reg73=reg31*reg68;
    T reg86=reg45*reg72; reg37=reg37/reg61; reg49=reg49/reg61; reg42=reg58+reg42; reg58=reg6*reg5;
    reg74=reg74/reg61; reg39=reg39/reg61; T reg87=reg63*reg59; reg78=reg36+reg78; T reg88=reg55*reg10;
    reg75=reg76-reg75; reg76=reg77*reg8; reg11=reg22*reg11; T reg89=reg69*reg59; reg26=reg60+reg26;
    reg9=reg77*reg9; reg60=reg14*reg0; reg66=reg34+reg66; reg59=reg80*reg59; reg34=reg69*reg32;
    reg83=reg57+reg83; reg41=reg41*reg23; reg77=reg77*reg10; reg71=reg28-reg71; reg7=reg7*reg24;
    reg54=reg33+reg54; reg28=reg62*reg23; reg56=reg56*reg24; reg33=reg17*reg72; reg57=reg16*reg68;
    reg21=reg25+reg21; T reg90=reg67*reg1; T reg91=reg65*reg4; T reg92=2*reg81; reg46=reg25+reg46;
    reg25=reg63*reg32; reg30=reg64+reg30; reg82=reg20-reg82; reg38=reg40+reg38; reg32=reg80*reg32;
    reg20=2*reg15; reg40=reg79*reg10; reg62=reg62*reg24; reg64=reg63*reg49; T reg93=reg45*reg74;
    T reg94=reg6*reg0; T reg95=reg80*reg74; T reg96=reg17*reg49; T reg97=reg60*reg91; reg55=reg55*reg23;
    reg79=reg79*reg24; reg82=reg82*reg22; reg83=reg83*reg24; T reg98=reg7+reg41; reg56=reg11-reg56;
    reg21=reg21/reg61; reg40=reg28-reg40; reg46=reg46/reg61; reg25=reg30+reg25; reg11=reg90*reg84;
    reg28=reg16*reg39; reg30=reg31*reg37; reg78=reg78/reg61; reg48=reg48/reg61; T reg99=reg20*reg1;
    T reg100=reg13*reg0; reg59=reg66+reg59; reg66=reg76*reg10; reg10=reg77*reg10; T reg101=reg9*reg24;
    reg9=reg9*reg23; reg62=reg88+reg62; reg88=reg52*reg71; reg86=reg73+reg86; reg73=reg53*reg2;
    T reg102=reg85*reg2; reg32=reg38+reg32; reg38=reg60*reg84; reg33=reg57+reg33; reg57=reg27*reg71;
    reg84=reg47*reg84; reg34=reg54+reg34; reg54=reg47*reg91; reg89=reg26+reg89; reg72=reg50*reg72;
    reg68=reg51*reg68; reg3=reg58*reg3; reg87=reg42+reg87; reg91=reg90*reg91; reg26=reg92*reg75;
    reg42=pow(reg75,2); T reg103=pow(reg81,2); T reg104=reg80*reg37; T reg105=reg63*reg39; T reg106=reg17*reg39;
    T reg107=reg45*reg37; T reg108=reg16*reg49; T reg109=reg31*reg74; reg57=reg33+reg57; reg33=reg63*reg3;
    T reg110=reg53*reg73; reg84=reg34+reg84; reg64=reg95+reg64; reg105=reg104+reg105; reg34=reg53*reg102;
    reg54=reg89+reg54; reg89=reg26*reg48; reg77=reg77*reg24; reg95=reg63*reg46; reg104=reg80*reg78;
    reg91=reg87+reg91; reg23=reg76*reg23; reg76=reg103*reg48; reg28=reg30+reg28; reg30=reg99*reg73;
    reg11=reg25+reg11; reg25=reg26*reg21; reg87=reg31*reg78; T reg111=reg16*reg46; reg40=reg40*reg22;
    reg96=reg93+reg96; reg93=reg42*reg21; reg72=reg68+reg72; reg68=reg99*reg102; reg98=reg98*reg8;
    reg71=reg29*reg71; reg83=reg82-reg83; reg82=reg45*reg78; T reg112=reg17*reg46; T reg113=reg79+reg55;
    reg62=reg62*reg24; T reg114=reg80*reg3; reg88=reg86+reg88; reg97=reg59+reg97; reg59=reg103*reg21;
    reg10=reg9-reg10; reg106=reg107+reg106; reg9=reg42*reg48; reg86=reg5*reg0; reg4=reg94*reg4;
    reg108=reg109+reg108; reg102=reg100*reg102; reg56=reg56/reg61; reg66=reg101+reg66; reg73=reg100*reg73;
    reg38=reg32+reg38; reg32=reg90*reg4; reg62=reg40-reg62; reg89=reg105+reg89; reg40=reg103*reg56;
    reg111=reg87+reg111; reg87=reg60*reg4; reg59=reg108+reg59; reg24=reg66*reg24; reg25=reg64+reg25;
    reg22=reg10*reg22; reg68=reg91+reg68; reg10=reg77+reg23; reg76=reg28+reg76; reg114=reg88+reg114;
    reg30=reg11+reg30; reg2=reg86*reg2; reg93=reg96+reg93; reg26=reg26*reg56; reg34=reg54+reg34;
    reg102=reg97+reg102; reg95=reg104+reg95; reg9=reg106+reg9; reg11=reg42*reg56; reg112=reg82+reg112;
    reg3=reg69*reg3; reg71=reg72+reg71; reg110=reg84+reg110; reg113=reg113*reg8; reg98=reg83-reg98;
    reg73=reg38+reg73; reg33=reg57+reg33; reg28=reg58*reg76; reg38=reg59*reg70; reg26=reg95+reg26;
    reg54=reg34*reg30; reg57=reg6*reg75; reg11=reg112+reg11; reg64=reg110*reg68; reg40=reg111+reg40;
    reg66=reg110*reg102; reg72=reg75*reg81; reg8=reg10*reg8; reg24=reg22-reg24; reg87=reg114+reg87;
    reg113=reg62-reg113; reg10=reg5*reg81; reg22=reg100*reg2; reg62=reg34*reg73; reg4=reg47*reg4;
    reg3=reg71+reg3; reg98=reg98/reg61; reg32=reg33+reg32; reg33=reg25*reg70; reg71=reg89*reg58;
    reg82=reg99*reg2; reg83=reg9*reg58; reg84=reg93*reg70; reg88=reg69*reg98; reg91=reg45*reg89;
    reg95=reg16*reg93; reg96=reg31*reg9; reg97=reg17*reg25; reg57=reg10+reg57; reg10=reg30*reg102;
    reg101=reg12*reg75; reg104=reg15*reg81; reg62=reg66-reg62; reg113=reg113/reg61; reg66=reg16*reg59;
    reg105=reg31*reg76; reg106=reg6*reg81; reg8=reg24-reg8; reg24=reg40*reg72; reg38=reg28+reg38;
    reg59=reg17*reg59; reg82=reg32+reg82; reg54=reg64-reg54; reg89=reg31*reg89; reg22=reg87+reg22;
    reg33=reg71+reg33; reg28=reg26*reg72; reg76=reg45*reg76; reg25=reg16*reg25; reg93=reg17*reg93;
    reg32=reg50*reg98; reg64=reg51*reg98; reg71=reg73*reg68; reg87=reg5*reg75; reg84=reg83+reg84;
    reg83=reg11*reg72; reg9=reg45*reg9; reg2=reg53*reg2; reg4=reg3+reg4; reg97=reg91+reg97;
    reg3=reg42*reg26; reg91=elem.pos(2)[1]-elem.pos(0)[1]; reg107=elem.pos(1)[1]-elem.pos(0)[1]; reg108=reg42*reg40; reg109=reg42*reg11;
    reg93=reg9+reg93; reg59=reg76+reg59; reg95=reg96+reg95; reg11=reg103*reg11; reg25=reg89+reg25;
    reg26=reg103*reg26; reg24=reg38+reg24; reg9=reg69*reg64; reg83=reg84+reg83; reg38=reg69*reg32;
    reg28=reg33+reg28; reg33=reg69*reg88; reg2=reg4+reg2; reg101=reg104+reg101; reg71=reg10-reg71;
    reg40=reg103*reg40; reg66=reg105+reg66; reg61=reg8/reg61; reg4=reg12*reg81; reg8=reg15*reg75;
    reg10=reg62*reg82; reg76=elem.pos(2)[0]-elem.pos(0)[0]; reg84=reg113*reg57; reg87=reg113*reg87; reg89=reg22*reg54;
    reg106=reg113*reg106; reg96=elem.pos(1)[0]-elem.pos(0)[0]; reg104=reg64*reg19; reg8=reg61*reg8; reg105=reg106*reg57;
    reg9=reg24+reg9; reg24=reg61*reg101; reg111=reg88*reg18; reg26=reg25+reg26; reg25=reg22*reg68;
    reg112=reg2*reg102; reg109=reg93+reg109; reg93=reg32*reg19; reg13=reg13*reg75; reg14=reg14*reg81;
    reg32=reg32*reg18; reg11=reg95+reg11; reg102=reg102*reg82; reg95=reg96*reg91; reg10=reg89-reg10;
    reg89=reg2*reg71; reg40=reg66+reg40; reg88=reg88*reg19; reg66=reg34*reg82; reg3=reg97+reg3;
    reg68=reg2*reg68; reg64=reg64*reg18; reg34=reg34*reg22; reg97=reg84*reg57; reg33=reg28+reg33;
    reg28=reg76*reg107; reg114=reg87*reg57; reg38=reg83+reg38; reg4=reg61*reg4; reg108=reg59+reg108;
    reg59=reg73*reg82; reg67=reg67*reg81; reg73=reg2*reg73; reg68=reg66-reg68; reg66=reg110*reg22;
    reg2=reg2*reg30; reg82=reg110*reg82; reg89=reg10+reg89; reg20=reg20*reg75; reg102=reg25-reg102;
    reg112=reg34-reg112; reg22=reg30*reg22; reg10=reg106*reg14; reg32=reg11+reg32; reg64=reg40+reg64;
    reg11=reg87*reg14; reg28=reg95-reg28; reg111=reg26+reg111; reg25=reg84*reg14; reg105=reg9+reg105;
    reg9=reg4*reg101; reg114=reg38+reg114; reg104=reg108+reg104; reg106=reg106*reg13; reg26=reg8*reg101;
    reg97=reg33+reg97; reg93=reg109+reg93; reg87=reg87*reg13; reg30=reg24*reg101; reg88=reg3+reg88;
    reg84=reg84*reg13; reg10=reg64+reg10; reg3=reg4*reg67; reg91=reg91/reg28; reg2=reg82-reg2;
    reg11=reg32+reg11; reg76=reg76/reg28; reg107=reg107/reg28; reg96=reg96/reg28; reg32=reg8*reg67;
    reg68=reg68/reg89; reg25=reg111+reg25; reg33=reg24*reg67; reg30=reg97+reg30; reg9=reg105+reg9;
    reg26=reg114+reg26; reg106=reg104+reg106; reg4=reg4*reg20; reg24=reg24*reg20; reg73=reg66-reg73;
    reg112=reg112/reg89; reg102=reg102/reg89; reg84=reg88+reg84; reg59=reg22-reg59; reg22=1-(*f.m).resolution;
    reg8=reg8*reg20; reg87=reg93+reg87; reg34=reg107-reg91; reg4=reg106+reg4; reg73=reg73/reg89;
    reg62=reg62/reg89; reg8=reg87+reg8; reg2=reg2/reg89; reg54=reg54/reg89; reg24=reg84+reg24;
    reg38=reg76-reg96; reg40=reg9*reg22; reg64=reg26*reg22; reg3=reg10+reg3; reg10=reg30*reg22;
    reg71=reg71/reg89; reg89=reg59/reg89; reg59=(*f.m).resolution*reg68; reg66=(*f.m).resolution*reg112; reg82=(*f.m).resolution*reg102;
    reg32=reg11+reg32; reg33=reg25+reg33; reg11=0.5*reg107; reg25=(*f.m).resolution*reg62; reg83=0.5*reg96;
    reg84=0.5*reg34; reg87=0.5*reg38; reg88=0.5*reg76; reg93=(*f.m).resolution*reg54; reg95=(*f.m).resolution*reg71;
    reg97=reg3*reg22; reg66=reg64-reg66; reg64=(*f.m).resolution*reg89; reg104=reg32*reg22; reg105=reg33*reg22;
    reg106=reg22*reg4; reg108=reg22*reg8; reg109=reg22*reg24; reg82=reg10+reg82; reg10=(*f.m).resolution*reg73;
    reg110=(*f.m).resolution*reg2; reg59=reg40+reg59; reg40=0.5*reg91; reg111=reg66*reg76; reg114=reg82*reg40;
    reg10=reg108+reg10; reg64=reg109-reg64; reg108=reg66*reg96; reg109=reg91*reg59; T reg115=reg82*reg88;
    reg110=reg106-reg110; reg106=reg11*reg82; reg95=reg105+reg95; reg105=reg82*reg83; T reg116=reg59*reg107;
    reg93=reg97+reg93; reg97=reg34*reg59; T reg117=reg87*reg82; T reg118=reg82*reg84; T reg119=reg38*reg66;
    reg25=reg104-reg25; reg104=reg95*reg83; reg119=reg118+reg119; reg118=reg93*reg107; T reg120=reg25*reg96;
    reg114=reg114-reg111; T reg121=reg38*reg25; T reg122=reg64*reg84; T reg123=reg95*reg84; T reg124=reg25*reg76;
    T reg125=reg11*reg64; T reg126=reg10*reg96; reg105=reg105-reg116; T reg127=reg95*reg40; T reg128=reg11*reg95;
    T reg129=reg95*reg88; T reg130=reg10*reg76; T reg131=reg87*reg95; T reg132=reg64*reg88; T reg133=reg91*reg93;
    T reg134=reg64*reg40; T reg135=reg91*reg110; T reg136=reg34*reg93; T reg137=reg38*reg10; reg108=reg108-reg106;
    T reg138=reg110*reg107; T reg139=reg64*reg83; reg97=reg117+reg97; reg109=reg109-reg115; reg114=2*reg114;
    reg108=2*reg108; reg127=reg127-reg124; reg137=reg122+reg137; reg135=reg135-reg132; reg120=reg120-reg128;
    reg104=reg104-reg118; reg133=reg133-reg129; reg109=2*reg109; reg139=reg139-reg138; reg134=reg134-reg130;
    reg126=reg126-reg125; reg136=reg131+reg136; reg97=2*reg97; reg121=reg123+reg121; reg105=2*reg105;
    reg119=2*reg119; reg117=reg87*reg108; reg122=reg76*reg134; reg123=reg114*reg84; reg131=reg119*reg84;
    T reg140=reg38*reg134; T reg141=reg108*reg40; T reg142=reg114*reg88; T reg143=reg127*reg91; T reg144=reg76*reg126;
    T reg145=reg38*reg135; T reg146=reg109*reg84; T reg147=reg11*reg108; T reg148=reg96*reg126; T reg149=reg87*reg119;
    T reg150=reg34*reg136; T reg151=reg87*reg97; T reg152=reg105*reg40; T reg153=reg137*reg38; T reg154=reg120*reg107;
    T reg155=reg108*reg83; T reg156=reg76*reg139; T reg157=reg34*reg121; T reg158=reg104*reg107; T reg159=reg34*reg104;
    T reg160=reg108*reg88; T reg161=reg91*reg120; reg120=reg34*reg120; T reg162=reg109*reg88; T reg163=reg133*reg91;
    T reg164=reg105*reg83; T reg165=reg87*reg114; reg108=reg108*reg84; T reg166=reg87*reg105; T reg167=reg34*reg127;
    reg126=reg38*reg126; T reg168=reg105*reg88; reg104=reg104*reg91; T reg169=reg114*reg40; T reg170=reg34*reg133;
    T reg171=reg87*reg109; T reg172=reg38*reg139; T reg173=reg105*reg84; reg149=reg157+reg149; reg156=reg152-reg156;
    reg160=reg161-reg160; reg158=reg164-reg158; reg142=reg143-reg142; reg162=reg163-reg162; reg117=reg120+reg117;
    reg154=reg155-reg154; reg123=reg140+reg123; reg144=reg141-reg144; reg151=reg150+reg151; reg165=reg167+reg165;
    reg173=reg172+reg173; reg166=reg159+reg166; reg168=reg104-reg168; reg108=reg126+reg108; reg147=reg148-reg147;
    reg131=reg153+reg131; reg146=reg145+reg146; reg171=reg170+reg171; reg122=reg169-reg122; reg156=reg156*reg28;
    reg160=reg160*reg28; reg168=reg168*reg28; reg122=reg122*reg28; reg142=reg142*reg28; reg162=reg28*reg162;
    reg108=reg28*reg108; reg166=reg166*reg28; reg173=reg28*reg173; reg165=reg165*reg28; reg171=reg171*reg28;
    reg117=reg117*reg28; reg131=reg131*reg28; reg123=reg28*reg123; reg146=reg146*reg28; reg147=reg28*reg147;
    reg154=reg154*reg28; reg151=reg151*reg28; reg158=reg158*reg28; reg144=reg144*reg28; reg149=reg149*reg28;
    T tmp_5_5=ponderation*reg147; T tmp_1_3=ponderation*reg123; T tmp_0_1=ponderation*reg149; T tmp_2_3=ponderation*reg142; T tmp_1_4=ponderation*reg173;
    T tmp_0_0=ponderation*reg151; T tmp_1_5=ponderation*reg108; T tmp_2_4=ponderation*reg168; T tmp_2_2=ponderation*reg162; T tmp_0_4=ponderation*reg166;
    T tmp_2_5=ponderation*reg160; T tmp_0_3=ponderation*reg165; T tmp_0_2=ponderation*reg171; T tmp_0_5=ponderation*reg117; T tmp_1_1=ponderation*reg131;
    T tmp_1_2=ponderation*reg146; T tmp_3_3=ponderation*reg122; T tmp_4_5=ponderation*reg154; T tmp_4_4=ponderation*reg158; T tmp_3_5=ponderation*reg144;
    T tmp_3_4=ponderation*reg156;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_3;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_4;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_5;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_3;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_4;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_5;
    matrix(indices[1]+0,indices[1]+0) += tmp_2_2;
    matrix(indices[1]+0,indices[1]+1) += tmp_2_3;
    matrix(indices[1]+0,indices[2]+0) += tmp_2_4;
    matrix(indices[1]+0,indices[2]+1) += tmp_2_5;
    matrix(indices[1]+1,indices[1]+1) += tmp_3_3;
    matrix(indices[1]+1,indices[2]+0) += tmp_3_4;
    matrix(indices[1]+1,indices[2]+1) += tmp_3_5;
    matrix(indices[2]+0,indices[2]+0) += tmp_4_4;
    matrix(indices[2]+0,indices[2]+1) += tmp_4_5;
    matrix(indices[2]+1,indices[2]+1) += tmp_5_5;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim,unsigned symmetric_version>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=2*(*f.m).shear_modulus_13; T reg1=2*(*f.m).shear_modulus_23; reg1=1.0/reg1; T reg2=2*(*f.m).shear_modulus_12; reg0=1.0/reg0;
    T reg3=reg0*reg1; reg2=1.0/reg2; T reg4=pow((*f.m).v1[0],2); T reg5=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg6=1.0/(*f.m).elastic_modulus_3;
    T reg7=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg8=reg2*reg3; T reg9=pow((*f.m).v2[1],2); T reg10=pow((*f.m).v2[0],2); T reg11=pow((*f.m).v1[1],2);
    T reg12=reg5*reg8; T reg13=reg7*reg8; T reg14=reg6*reg8; reg9=reg10+reg9; reg11=reg4+reg11;
    reg4=pow((*f.m).v1[2],2); reg10=pow((*f.m).v2[2],2); T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg17=reg16*reg14;
    T reg18=reg5*reg12; T reg19=reg5*reg13; reg4=reg11+reg4; reg10=reg9+reg10; reg9=reg15*reg14;
    reg11=reg16*reg12; T reg20=reg15*reg13; reg19=reg17+reg19; reg4=pow(reg4,0.5); T reg21=1.0/(*f.m).elastic_modulus_1;
    reg10=pow(reg10,0.5); reg18=reg9-reg18; reg9=(*f.m).v1[2]/reg4; T reg22=(*f.m).v1[1]/reg4; T reg23=reg16*reg19;
    T reg24=(*f.m).v2[2]/reg10; T reg25=(*f.m).v2[1]/reg10; T reg26=reg21*reg18; T reg27=reg11+reg20; reg14=reg21*reg14;
    T reg28=reg7*reg13; T reg29=reg16*reg8; T reg30=reg2*reg1; reg4=(*f.m).v1[0]/reg4; T reg31=reg6*reg3;
    T reg32=reg7*reg3; reg3=reg5*reg3; reg10=(*f.m).v2[0]/reg10; reg23=reg26-reg23; reg26=reg7*reg27;
    reg8=reg15*reg8; T reg33=reg7*reg12; T reg34=reg9*reg25; T reg35=reg22*reg24; T reg36=2*reg10;
    T reg37=reg30*reg6; T reg38=reg32*reg5; T reg39=reg3*reg5; T reg40=reg31*reg16; reg31=reg31*reg15;
    T reg41=reg30*reg7; T reg42=reg2*reg0; reg13=reg16*reg13; reg30=reg30*reg5; T reg43=2*reg4;
    T reg44=reg29*reg7; reg12=reg21*reg12; reg28=reg14-reg28; reg33=reg17+reg33; reg14=reg9*reg10;
    reg26=reg23-reg26; reg17=reg35-reg34; reg23=reg7*reg8; T reg45=reg4*reg24; reg18=reg18/reg26;
    T reg46=2*reg17; T reg47=reg42*reg7; T reg48=reg4*reg25; T reg49=reg14-reg45; reg32=reg32*reg15;
    T reg50=reg42*reg5; reg3=reg3*reg16; T reg51=reg22*reg10; T reg52=reg43*reg22; reg29=reg29*reg16;
    reg13=reg12+reg13; reg8=reg21*reg8; T reg53=reg41*reg5; reg44=reg12+reg44; reg33=reg33/reg26;
    reg42=reg42*reg6; reg12=pow(reg22,2); reg28=reg28/reg26; T reg54=pow(reg4,2); reg19=reg19/reg26;
    reg23=reg11+reg23; reg39=reg31-reg39; reg31=reg36*reg25; T reg55=pow(reg25,2); T reg56=pow(reg10,2);
    reg38=reg40+reg38; reg40=reg30*reg5; T reg57=reg37*reg16; reg37=reg37*reg15; T reg58=pow(reg9,2);
    reg27=reg27/reg26; reg23=reg23/reg26; T reg59=reg31*reg28; T reg60=reg52*reg33; T reg61=reg55*reg28;
    T reg62=reg12*reg33; T reg63=reg56*reg28; T reg64=reg54*reg33; T reg65=reg31*reg19; T reg66=reg52*reg18;
    T reg67=reg55*reg19; T reg68=reg12*reg18; T reg69=reg47*reg5; T reg70=reg50*reg5; T reg71=reg42*reg16;
    reg42=reg42*reg15; T reg72=pow(reg24,2); reg57=reg53+reg57; reg53=reg48-reg51; reg39=reg39*reg21;
    reg41=reg41*reg15; T reg73=reg46*reg49; T reg74=pow(reg49,2); reg30=reg30*reg16; reg38=reg38*reg16;
    T reg75=reg3+reg32; T reg76=pow(reg17,2); reg29=reg8-reg29; reg13=reg13/reg26; reg8=reg54*reg18;
    T reg77=reg56*reg19; reg44=reg44/reg26; reg40=reg37-reg40; reg47=reg47*reg15; reg57=reg57*reg16;
    reg50=reg50*reg16; reg37=reg30+reg41; reg38=reg39-reg38; reg39=reg54*reg23; T reg78=reg56*reg44;
    reg59=reg60+reg59; reg70=reg42-reg70; reg42=reg12*reg23; reg69=reg71+reg69; reg60=reg55*reg44;
    reg75=reg75*reg7; reg67=reg68+reg67; reg68=reg74*reg27; reg71=reg76*reg27; T reg79=reg72*reg28;
    T reg80=reg58*reg33; T reg81=reg74*reg13; reg61=reg62+reg61; reg40=reg40*reg21; reg29=reg29/reg26;
    reg62=reg76*reg13; reg63=reg64+reg63; reg64=reg31*reg44; T reg82=reg73*reg13; T reg83=reg73*reg27;
    reg65=reg66+reg65; reg66=reg52*reg23; reg77=reg8+reg77; reg8=reg72*reg19; T reg84=reg58*reg18;
    T reg85=pow(reg53,2); T reg86=reg72*reg44; T reg87=reg58*reg23; reg37=reg37*reg7; T reg88=reg74*reg29;
    reg60=reg42+reg60; reg70=reg70*reg21; reg42=reg76*reg29; reg78=reg39+reg78; reg69=reg69*reg16;
    reg39=reg50+reg47; reg68=reg67+reg68; reg82=reg59+reg82; reg8=reg84+reg8; reg59=reg85*reg27;
    reg83=reg65+reg83; reg57=reg40-reg57; reg75=reg38-reg75; reg71=reg77+reg71; reg62=reg63+reg62;
    reg38=2*reg25; reg81=reg61+reg81; reg40=reg36*reg24; reg79=reg80+reg79; reg61=reg85*reg13;
    reg63=reg43*reg9; reg65=2*reg22; reg67=reg10*reg25; reg64=reg66+reg64; reg66=reg4*reg22;
    reg73=reg73*reg29; reg46=reg46*reg53; reg77=reg63*reg18; reg80=reg40*reg19; reg84=reg49*reg17;
    T reg89=reg55*reg81; reg75=reg75/reg26; T reg90=reg43*reg10; T reg91=reg56*reg81; T reg92=reg54*reg68;
    T reg93=reg54*reg83; T reg94=reg12*reg68; T reg95=reg4*reg49; T reg96=reg55*reg82; T reg97=reg54*reg71;
    T reg98=reg56*reg62; T reg99=reg56*reg82; T reg100=reg22*reg17; reg61=reg79+reg61; reg79=reg55*reg62;
    reg68=reg68*reg66; reg81=reg81*reg67; T reg101=reg65*reg9; T reg102=reg38*reg24; T reg103=reg12*reg83;
    reg83=reg83*reg66; reg82=reg82*reg67; T reg104=reg12*reg71; T reg105=reg22*reg25; reg86=reg87+reg86;
    reg87=reg85*reg29; T reg106=reg4*reg10; reg88=reg60+reg88; reg37=reg57-reg37; reg57=2*reg9;
    reg42=reg78+reg42; reg60=reg65*reg25; reg73=reg64+reg73; reg69=reg70-reg69; reg39=reg39*reg7;
    reg64=reg40*reg28; reg59=reg8+reg59; reg48=reg51+reg48; reg8=2*reg49; reg51=reg63*reg33;
    reg81=reg68+reg81; reg68=reg88*reg84; reg70=reg76*reg73; reg87=reg86+reg87; reg89=reg94+reg89;
    reg78=reg22*reg49; reg86=reg54*reg59; reg94=reg56*reg61; T reg107=reg76*reg88; reg95=reg100+reg95;
    reg82=reg83+reg82; reg83=reg73*reg84; reg100=reg15*reg60; T reg108=reg16*reg90; T reg109=reg55*reg15;
    T reg110=reg56*reg16; T reg111=reg4*reg17; T reg112=reg55*reg61; T reg113=reg12*reg59; reg64=reg51+reg64;
    reg51=reg46*reg13; reg71=reg66*reg71; reg33=reg101*reg33; reg62=reg62*reg67; reg28=reg102*reg28;
    T reg114=reg40*reg44; T reg115=reg63*reg23; T reg116=reg74*reg42; reg91=reg92+reg91; reg79=reg104+reg79;
    reg88=reg74*reg88; reg92=reg106*reg75; reg99=reg93+reg99; reg93=reg48*reg75; reg104=reg105*reg75;
    reg73=reg74*reg73; reg96=reg103+reg96; reg103=reg9*reg24; reg19=reg102*reg19; reg18=reg101*reg18;
    T reg117=reg10*reg49; T reg118=reg46*reg27; reg57=reg24*reg57; reg80=reg77+reg80; reg77=reg56*reg21;
    T reg119=reg55*reg16; reg37=reg37/reg26; reg8=reg8*reg53; T reg120=reg76*reg42; reg39=reg69-reg39;
    reg98=reg97+reg98; reg69=reg25*reg17; reg97=reg21*reg90; T reg121=reg16*reg60; T reg122=reg74*reg87;
    T reg123=reg10*reg17; reg112=reg113+reg112; reg23=reg101*reg23; reg27=reg8*reg27; reg46=reg46*reg29;
    reg26=reg39/reg26; reg19=reg18+reg19; reg65=reg65*reg49; reg114=reg115+reg114; reg43=reg43*reg17;
    reg13=reg8*reg13; reg18=reg56*reg7; reg28=reg33+reg28; reg116=reg79+reg116; reg33=reg92*reg60;
    reg39=reg55*reg5; reg118=reg80+reg118; reg79=reg25*reg49; reg80=reg76*reg87; reg113=reg104*reg60;
    reg88=reg89+reg88; reg89=reg7*reg90; reg115=reg104*reg90; reg51=reg64+reg51; reg64=reg93*reg90;
    reg73=reg96+reg73; reg70=reg99+reg70; reg96=reg93*reg60; reg99=reg37*reg95; T reg124=reg92*reg90;
    reg107=reg91+reg107; reg117=reg69+reg117; reg59=reg59*reg66; reg61=reg61*reg67; reg69=reg5*reg60;
    reg91=reg72*reg7; reg119=reg77-reg119; reg77=reg7*reg57; reg108=reg100-reg108; reg121=reg97-reg121;
    reg97=reg5*reg57; reg100=reg9*reg53; reg110=reg109-reg110; reg109=reg72*reg5; reg120=reg98+reg120;
    reg98=reg103*reg75; reg111=reg37*reg111; reg93=reg48*reg93; reg42=reg42*reg84; reg83=reg82+reg83;
    reg78=reg37*reg78; reg68=reg81+reg68; reg104=reg48*reg104; reg62=reg71+reg62; reg44=reg102*reg44;
    reg94=reg86+reg94; reg71=(*f.m).alpha_1*reg12; reg46=reg114+reg46; reg81=reg55*(*f.m).alpha_2; reg44=reg23+reg44;
    reg29=reg8*reg29; reg14=reg45+reg14; reg21=reg54*reg21; reg38=reg38*reg49; reg80=reg94+reg80;
    reg8=reg98*reg90; reg124=reg120+reg124; reg23=reg26*reg117; reg79=reg26*reg79; reg123=reg26*reg123;
    reg39=reg18+reg39; reg18=reg55*reg51; reg45=reg72*reg6; reg82=reg12*reg118; reg100=reg37*reg100;
    reg27=reg19+reg27; reg69=reg89+reg69; reg42=reg62+reg42; reg92=reg48*reg92; reg19=reg111*reg43;
    reg62=reg99*reg95; reg122=reg112+reg122; reg86=reg98*reg60; reg93=reg83+reg93; reg83=reg12*reg16;
    reg104=reg68+reg104; reg68=reg78*reg95; reg61=reg59+reg61; reg87=reg87*reg84; reg57=reg6*reg57;
    reg91=reg119-reg91; reg97=reg108-reg97; reg77=reg121-reg77; reg109=reg110-reg109; reg16=reg54*reg16;
    reg15=reg12*reg15; reg64=reg70+reg64; reg59=reg99*reg43; reg13=reg28+reg13; reg28=reg24*reg53;
    reg99=reg99*reg65; reg70=reg56*(*f.m).alpha_2; reg96=reg73+reg96; reg73=reg9*reg17; reg89=reg4*reg53;
    reg94=(*f.m).alpha_1*reg54; reg113=reg88+reg113; reg88=reg111*reg65; reg108=reg78*reg43; reg33=reg116+reg33;
    reg115=reg107+reg115; reg78=reg78*reg65; reg36=reg36*reg17; reg107=reg56*reg51; reg110=reg54*reg118;
    reg69=reg57-reg69; reg57=reg55*reg109; reg112=reg12*reg5; reg114=reg56*reg91; reg116=reg123*reg38;
    reg70=reg94+reg70; reg39=reg45-reg39; reg45=reg24*reg17; reg94=reg10*reg53; reg88=reg33+reg88;
    reg35=reg34+reg35; reg33=reg76*(*f.m).alpha_3; reg34=reg54*reg7; reg16=reg15-reg16; reg5=reg58*reg5;
    reg83=reg21-reg83; reg7=reg58*reg7; reg62=reg93+reg62; reg15=reg23*reg117; reg118=reg118*reg66;
    reg51=reg51*reg67; reg21=reg54*reg91; reg93=reg12*reg109; reg99=reg96+reg99; reg96=reg23*reg38;
    reg8=reg80+reg8; reg80=reg100*reg43; reg89=reg73+reg89; reg73=reg9*reg49; reg119=reg22*reg53;
    reg120=reg54*reg77; reg121=reg12*reg97; reg18=reg82+reg18; reg82=reg74*reg46; T reg125=reg105*reg97;
    T reg126=reg106*reg77; reg19=reg124+reg19; reg124=reg123*reg36; T reg127=reg12*reg27; T reg128=reg55*reg13;
    reg109=reg105*reg109; reg91=reg106*reg91; T reg129=reg14*reg75; reg29=reg44+reg29; reg59=reg64+reg59;
    reg44=(*f.m).alpha_1*reg58; reg28=reg26*reg28; reg64=reg72*(*f.m).alpha_2; reg97=reg55*reg97; reg77=reg56*reg77;
    reg23=reg23*reg36; reg107=reg110+reg107; reg110=reg76*reg46; reg78=reg113+reg78; reg113=reg79*reg38;
    T reg130=reg74*(*f.m).alpha_3; T reg131=reg54*reg27; T reg132=reg56*reg13; reg108=reg115+reg108; reg92=reg42+reg92;
    reg111=reg111*reg95; reg42=reg79*reg36; reg68=reg104+reg68; reg79=reg79*reg117; reg86=reg122+reg86;
    reg87=reg61+reg87; reg98=reg48*reg98; reg61=reg100*reg65; reg81=reg71+reg81; reg116=reg88+reg116;
    reg124=reg19+reg124; reg128=reg127+reg128; reg19=reg129*reg60; reg82=reg18+reg82; reg96=reg99+reg96;
    reg113=reg78+reg113; reg61=reg86+reg61; reg18=reg28*reg38; reg71=reg67*(*f.m).alpha_2; reg132=reg131+reg132;
    reg76=reg76*reg29; reg111=reg92+reg111; reg123=reg123*reg117; reg79=reg68+reg79; reg98=reg87+reg98;
    reg100=reg100*reg95; reg42=reg108+reg42; reg5=reg16-reg5; reg7=reg83-reg7; reg15=reg62+reg15;
    reg51=reg118+reg51; reg46=reg46*reg84; reg93=reg21+reg93; reg16=reg58*reg39; reg27=reg27*reg66;
    reg13=reg13*reg67; reg119=reg73+reg119; reg121=reg120+reg121; reg21=reg58*reg69; reg80=reg8+reg80;
    reg8=reg103*reg69; reg125=reg126+reg125; reg4=reg4*reg9; reg10=reg10*reg24; reg6=reg58*reg6;
    reg112=reg34+reg112; reg57=reg114+reg57; reg34=reg72*reg39; reg94=reg45+reg94; reg45=reg24*reg49;
    reg62=reg25*reg53; reg74=reg74*reg29; reg67=reg67*reg2; reg68=reg48*reg2; reg109=reg91+reg109;
    reg39=reg103*reg39; reg33=reg70+reg33; reg130=reg81+reg130; reg89=reg37*reg89; reg75=reg35*reg75;
    reg64=reg44+reg64; reg85=reg85*(*f.m).alpha_3; reg44=reg28*reg36; reg69=reg72*reg69; reg97=reg77+reg97;
    reg23=reg59+reg23; reg110=reg107+reg110; reg59=reg129*reg90; reg70=(*f.m).alpha_1*reg66; reg9=reg22*reg9;
    reg17=reg53*reg17; reg24=reg25*reg24; reg44=reg80+reg44; reg22=reg52*reg67; reg16=reg93+reg16;
    reg25=reg12*reg5; reg73=reg54*reg7; reg77=reg14*reg0; reg78=reg10*reg0; reg112=reg6-reg112;
    reg85=reg64+reg85; reg6=reg31*reg68; reg69=reg97+reg69; reg59=reg110+reg59; reg64=reg89*reg43;
    reg71=reg70+reg71; reg70=reg84*(*f.m).alpha_3; reg123=reg111+reg123; reg100=reg98+reg100; reg28=reg28*reg117;
    reg80=(*f.m).alpha_1*reg4; reg46=reg51+reg46; reg129=reg48*reg129; reg10=reg10*(*f.m).alpha_2; reg13=reg27+reg13;
    reg84=reg29*reg84; reg21=reg121+reg21; reg27=reg48*reg68; reg8=reg125+reg8; reg68=reg52*reg68;
    reg29=reg56*reg7; reg51=reg55*reg5; reg81=reg48*reg67; reg39=reg109+reg39; reg83=reg23*reg79;
    reg94=reg26*reg94; reg76=reg132+reg76; reg90=reg75*reg90; reg18=reg61+reg18; reg61=reg79*reg96;
    reg86=reg42*reg15; reg87=reg15*reg113; reg88=reg124*reg33; reg34=reg57+reg34; reg57=reg42*reg130;
    reg67=reg31*reg67; reg19=reg82+reg19; reg62=reg45+reg62; reg74=reg128+reg74; reg60=reg75*reg60;
    reg45=reg89*reg65; reg82=reg116*reg33; reg91=reg113*reg130; reg119=reg37*reg119; reg51=reg29+reg51;
    reg29=reg63*reg77; reg68=reg21+reg68; reg62=reg26*reg62; reg83=reg86-reg83; reg84=reg13+reg84;
    reg75=reg48*reg75; reg72=reg72*reg112; reg61=reg87-reg61; reg13=reg24*(*f.m).alpha_2; reg21=(*f.m).alpha_1*reg9;
    reg17=reg17*(*f.m).alpha_3; reg10=reg80+reg10; reg70=reg71+reg70; reg60=reg74+reg60; reg65=reg119*reg65;
    reg2=reg66*reg2; reg26=reg94*reg38; reg45=reg19+reg45; reg43=reg119*reg43; reg90=reg76+reg90;
    reg27=reg8+reg27; reg8=reg14*reg77; reg57=reg88+reg57; reg44=reg44*reg85; reg5=reg105*reg5;
    reg7=reg106*reg7; reg19=reg79*reg130; reg37=reg123*reg33; reg91=reg82+reg91; reg18=reg18*reg85;
    reg77=reg40*reg77; reg6=reg69+reg6; reg64=reg59+reg64; reg59=reg94*reg36; reg67=reg34+reg67;
    reg28=reg100+reg28; reg129=reg46+reg129; reg89=reg89*reg95; reg22=reg16+reg22; reg16=reg63*reg78;
    reg34=reg14*reg78; reg46=reg23*reg113; reg58=reg58*reg112; reg25=reg73+reg25; reg78=reg40*reg78;
    reg49=reg53*reg49; reg24=reg24*reg1; reg53=reg42*reg96; reg81=reg39+reg81; reg39=reg35*reg1;
    reg34=reg81+reg34; reg44=reg57+reg44; reg57=reg23*reg70; reg112=reg103*reg112; reg5=reg7+reg5;
    reg29=reg68+reg29; reg7=reg101*reg39; reg49=reg49*(*f.m).alpha_3; reg95=reg119*reg95; reg75=reg84+reg75;
    reg13=reg21+reg13; reg85=reg28*reg85; reg19=reg37+reg19; reg65=reg60+reg65; reg38=reg62*reg38;
    reg18=reg91+reg18; reg21=reg96*reg70; reg94=reg94*reg117; reg17=reg10+reg17; reg10=reg102*reg39;
    reg77=reg6+reg77; reg16=reg22+reg16; reg6=reg31*reg2; reg59=reg64+reg59; reg89=reg129+reg89;
    reg22=reg116*reg83; reg28=reg101*reg24; reg37=reg124*reg61; reg60=reg52*reg2; reg78=reg67+reg78;
    reg72=reg51+reg72; reg51=reg35*reg24; reg46=reg53-reg46; reg24=reg102*reg24; reg58=reg25+reg58;
    reg43=reg90+reg43; reg39=reg35*reg39; reg36=reg62*reg36; reg26=reg45+reg26; reg0=reg4*reg0;
    reg8=reg27+reg8; reg60=reg58+reg60; reg21=reg18+reg21; reg26=reg26*reg17; reg4=reg124*reg15;
    reg94=reg89+reg94; reg18=reg123*reg96; reg51=reg34+reg51; reg10=reg77+reg10; reg25=reg15*reg116;
    reg22=reg37-reg22; reg36=reg43+reg36; reg27=reg123*reg46; reg1=reg9*reg1; reg2=reg48*reg2;
    reg95=reg75+reg95; reg117=reg62*reg117; reg112=reg5+reg112; reg63=reg63*reg0; reg24=reg78+reg24;
    reg57=reg44+reg57; reg59=reg59*reg17; reg5=reg23*reg123; reg40=reg40*reg0; reg49=reg13+reg49;
    reg9=reg15*reg70; reg85=reg19+reg85; reg28=reg16+reg28; reg7=reg29+reg7; reg39=reg8+reg39;
    reg6=reg72+reg6; reg38=reg65+reg38; reg0=reg14*reg0; reg117=reg95+reg117; reg2=reg112+reg2;
    reg8=reg39*reg28; reg13=reg39*reg24; reg14=reg51*reg10; reg101=reg101*reg1; reg63=reg60+reg63;
    reg16=reg79*reg116; reg18=reg25-reg18; reg19=reg123*reg113; reg27=reg22+reg27; reg22=reg124*reg79;
    reg40=reg6+reg40; reg5=reg4-reg5; reg102=reg102*reg1; reg4=reg42*reg123; reg6=reg124*reg96;
    reg25=reg23*reg116; reg29=reg51*reg7; reg38=reg38*reg49; reg59=reg57+reg59; reg26=reg21+reg26;
    reg36=reg36*reg49; reg9=reg85+reg9; reg17=reg94*reg17; reg26=reg38+reg26; reg101=reg63+reg101;
    reg14=reg13-reg14; reg61=reg61/reg27; reg18=reg18/reg27; reg17=reg9+reg17; reg117=reg49*reg117;
    reg19=reg16-reg19; reg29=reg8-reg29; reg83=reg83/reg27; reg8=reg10*reg28; reg36=reg59+reg36;
    reg5=reg5/reg27; reg9=reg7*reg24; reg1=reg35*reg1; reg0=reg2+reg0; reg102=reg40+reg102;
    reg4=reg22-reg4; reg2=reg42*reg116; reg25=reg6-reg25; reg6=reg124*reg113; reg9=reg8-reg9;
    reg8=reg101*reg14; reg19=reg19/reg27; reg5=reg26*reg5; reg18=reg18*reg36; reg4=reg4/reg27;
    reg83=reg26*reg83; reg61=reg61*reg36; reg46=reg46/reg27; reg25=reg25/reg27; reg2=reg6-reg2;
    reg117=reg17+reg117; reg6=reg29*reg102; reg1=reg0+reg1; reg0=reg1*reg28; reg27=reg2/reg27;
    reg46=reg117*reg46; reg83=reg61-reg83; reg25=reg117*reg25; reg18=reg5-reg18; reg36=reg19*reg36;
    reg4=reg26*reg4; reg6=reg8-reg6; reg2=reg1*reg9; reg5=reg51*reg101; reg8=reg1*reg24;
    reg51=reg51*reg102; reg13=reg1*reg7; reg16=1-(*f.m).resolution; reg83=reg46+reg83; reg8=reg51-reg8;
    reg25=reg18-reg25; reg17=reg39*reg101; reg27=reg117*reg27; reg4=reg36-reg4; reg28=reg28*reg102;
    reg24=reg101*reg24; reg0=reg5-reg0; reg2=reg6+reg2; reg1=reg1*reg10; reg39=reg39*reg102;
    reg28=reg24-reg28; reg102=reg7*reg102; reg101=reg10*reg101; reg0=reg0/reg2; reg13=reg17-reg13;
    reg4=reg27+reg4; reg8=reg8/reg2; reg25=reg16*reg25; reg5=(*f.m).resolution*reg33; reg83=reg16*reg83;
    reg1=reg39-reg1; reg6=(*f.m).resolution*reg130; reg28=reg28/reg2; reg102=reg101-reg102; reg123=reg123*reg16;
    reg79=reg79*reg16; reg1=reg1/reg2; reg29=reg29/reg2; reg13=reg13/reg2; reg25=reg6+reg25;
    reg5=reg83+reg5; reg4=reg16*reg4; reg14=reg14/reg2; reg6=(*f.m).resolution*reg70; reg7=(*f.m).resolution*reg8;
    reg10=(*f.m).resolution*reg0; reg10=reg79-reg10; reg113=reg16*reg113; reg17=(*f.m).resolution*reg13; reg18=elem.pos(2)[1]-elem.pos(0)[1];
    reg7=reg123+reg7; reg19=elem.pos(1)[1]-elem.pos(0)[1]; reg25=(*f.m).deltaT*reg25; reg5=(*f.m).deltaT*reg5; reg15=reg15*reg16;
    reg21=(*f.m).resolution*reg1; reg22=(*f.m).resolution*reg29; reg124=reg124*reg16; reg42=reg42*reg16; reg24=(*f.m).resolution*reg14;
    reg4=reg6+reg4; reg116=reg16*reg116; reg6=(*f.m).resolution*reg28; reg26=elem.pos(2)[0]-elem.pos(0)[0]; reg9=reg9/reg2;
    reg27=elem.pos(1)[0]-elem.pos(0)[0]; reg2=reg102/reg2; reg6=reg15+reg6; reg15=(*f.m).resolution*reg9; reg17=reg113+reg17;
    reg21=reg116-reg21; reg23=reg23*reg16; reg34=reg27*reg18; reg35=(*f.m).resolution*reg2; reg22=reg42-reg22;
    reg4=(*f.m).deltaT*reg4; reg36=reg26*reg19; reg24=reg124+reg24; reg37=reg10*reg25; reg38=reg7*reg5;
    reg96=reg16*reg96; reg35=reg96-reg35; reg16=reg17*reg25; reg39=reg21*reg5; reg15=reg23+reg15;
    reg23=reg6*reg4; reg36=reg34-reg36; reg34=reg22*reg25; reg40=reg38+reg37; reg42=reg24*reg5;
    reg43=reg42+reg34; reg44=reg23+reg40; reg45=reg15*reg4; reg18=reg18/reg36; reg27=reg27/reg36;
    reg46=reg39+reg16; reg49=reg35*reg4; reg19=reg19/reg36; reg26=reg26/reg36; reg51=0.5*reg26;
    reg53=reg19-reg18; reg57=reg45+reg43; reg58=reg46+reg49; reg59=0.5*reg27; reg60=0.5*reg19;
    reg61=reg26-reg27; reg62=0.5*reg18; reg63=2*reg44; reg64=1-var_inter[0]; reg65=reg27*reg58;
    reg66=reg60*reg63; reg67=0.5*reg53; reg68=reg19*reg57; reg69=reg63*reg59; reg71=reg26*reg58;
    reg72=reg62*reg63; reg73=reg63*reg51; reg74=reg18*reg57; reg75=0.5*reg61; reg76=var_inter[1]*elem.f_vol_e[1];
    reg64=reg64-var_inter[1]; reg77=var_inter[1]*elem.f_vol_e[0]; reg78=var_inter[0]*elem.f_vol_e[1]; reg79=var_inter[0]*elem.f_vol_e[0]; reg80=reg68-reg69;
    reg81=reg71-reg72; reg82=reg66-reg65; reg83=reg53*reg57; reg84=reg73-reg74; reg85=reg75*reg63;
    reg86=reg63*reg67; reg87=reg61*reg58; reg80=reg80-reg77; reg82=reg82-reg76; reg84=reg84-reg79;
    reg88=elem.f_vol_e[1]*reg64; reg89=reg87+reg86; reg81=reg81-reg78; reg90=reg64*elem.f_vol_e[0]; reg91=reg83+reg85;
    reg81=reg36*reg81; reg82=reg36*reg82; reg84=reg36*reg84; reg80=reg36*reg80; reg92=reg88+reg89;
    reg93=reg90+reg91; reg81=ponderation*reg81; reg94=reg36*reg92; reg80=ponderation*reg80; reg82=ponderation*reg82;
    reg95=reg36*reg93; reg84=ponderation*reg84; T vec_2=-reg84; T vec_4=-reg80; T vec_5=-reg82;
    T vec_3=-reg81; reg80=ponderation*reg94; T vec_1=reg80; reg81=ponderation*reg95; T vec_0=reg81;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
sollicitation[indices[0]+0] += vec_0;
sollicitation[indices[0]+1] += vec_1;
sollicitation[indices[1]+0] += vec_2;
sollicitation[indices[1]+1] += vec_3;
sollicitation[indices[2]+0] += vec_4;
sollicitation[indices[2]+1] += vec_5;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim,unsigned symmetric_version>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_residual( TD ponderation, const TD *var_inter,
      const Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices ) { 
  #define PNODE(N) (*elem.node(N))
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; T reg2=2*(*f.m).shear_modulus_12; reg1=1.0/reg1; reg0=1.0/reg0;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=pow((*f.m).v2[1],2); T reg5=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg6=pow((*f.m).v2[0],2);
    T reg7=1.0/(*f.m).elastic_modulus_3; T reg8=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg9=pow((*f.m).v1[0],2); T reg10=pow((*f.m).v1[1],2); T reg11=reg2*reg3;
    T reg12=1.0/(*f.m).elastic_modulus_2; T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; reg4=reg6+reg4; reg6=reg5*reg11; T reg14=pow((*f.m).v2[2],2);
    reg10=reg9+reg10; reg9=pow((*f.m).v1[2],2); T reg15=reg8*reg11; T reg16=reg7*reg11; reg14=reg4+reg14;
    reg4=reg8*reg6; T reg17=reg8*reg15; T reg18=reg13*reg16; reg9=reg10+reg9; reg10=reg12*reg16;
    reg9=pow(reg9,0.5); T reg19=1.0/(*f.m).elastic_modulus_1; T reg20=reg12*reg6; T reg21=reg13*reg15; reg4=reg18+reg4;
    reg14=pow(reg14,0.5); reg17=reg10-reg17; reg10=(*f.m).v1[2]/reg9; T reg22=(*f.m).v2[1]/reg14; T reg23=(*f.m).v2[2]/reg14;
    T reg24=(*f.m).v1[1]/reg9; T reg25=reg19*reg17; T reg26=reg21+reg20; T reg27=reg13*reg4; T reg28=reg7*reg3;
    T reg29=reg5*reg15; T reg30=reg5*reg3; reg3=reg8*reg3; T reg31=reg24*reg23; T reg32=reg10*reg22;
    T reg33=reg2*reg0; reg9=(*f.m).v1[0]/reg9; reg16=reg19*reg16; T reg34=reg12*reg11; T reg35=reg5*reg6;
    T reg36=reg5*reg26; reg27=reg25-reg27; reg14=(*f.m).v2[0]/reg14; reg11=reg13*reg11; reg6=reg13*reg6;
    reg25=reg33*reg8; T reg37=reg11*reg5; reg15=reg19*reg15; reg35=reg16-reg35; reg29=reg18+reg29;
    reg16=reg2*reg1; reg18=reg33*reg5; T reg38=reg28*reg12; reg28=reg28*reg13; T reg39=reg3*reg8;
    T reg40=reg9*reg23; T reg41=reg30*reg8; T reg42=reg10*reg14; reg33=reg33*reg7; T reg43=2*reg14;
    reg36=reg27-reg36; reg27=reg5*reg34; T reg44=reg31-reg32; T reg45=2*reg9; T reg46=reg9*reg22;
    T reg47=reg24*reg14; reg34=reg19*reg34; reg35=reg35/reg36; reg4=reg4/reg36; reg27=reg21+reg27;
    T reg48=reg16*reg7; T reg49=reg18*reg8; reg29=reg29/reg36; reg37=reg15+reg37; T reg50=reg42-reg40;
    reg6=reg15+reg6; reg11=reg11*reg13; reg17=reg17/reg36; reg3=reg3*reg13; reg15=reg25*reg8;
    T reg51=reg16*reg8; reg30=reg30*reg12; T reg52=reg33*reg13; reg33=reg33*reg12; reg16=reg16*reg5;
    reg41=reg28+reg41; reg39=reg38-reg39; reg28=pow(reg24,2); reg38=reg45*reg24; T reg53=pow(reg9,2);
    T reg54=2*reg44; T reg55=pow(reg14,2); T reg56=pow(reg22,2); T reg57=reg43*reg22; reg39=reg39*reg19;
    T reg58=reg57*reg4; reg18=reg18*reg12; reg25=reg25*reg13; reg11=reg34-reg11; reg34=pow(reg23,2);
    reg6=reg6/reg36; reg37=reg37/reg36; T reg59=reg55*reg4; T reg60=reg53*reg17; reg26=reg26/reg36;
    T reg61=reg46-reg47; T reg62=reg38*reg17; T reg63=reg57*reg35; T reg64=reg38*reg29; T reg65=pow(reg10,2);
    T reg66=pow(reg44,2); T reg67=pow(reg50,2); T reg68=reg54*reg50; T reg69=reg28*reg29; reg52=reg49+reg52;
    reg49=reg56*reg35; T reg70=reg48*reg12; reg48=reg48*reg13; T reg71=reg51*reg8; T reg72=reg55*reg35;
    reg27=reg27/reg36; T reg73=reg16*reg8; reg15=reg33-reg15; reg33=reg53*reg29; T reg74=reg28*reg17;
    T reg75=reg3+reg30; T reg76=reg56*reg4; reg41=reg41*reg13; reg58=reg62+reg58; reg72=reg33+reg72;
    reg33=reg68*reg26; reg62=reg34*reg4; T reg77=reg65*reg17; T reg78=reg67*reg26; reg76=reg74+reg76;
    reg73=reg48+reg73; reg71=reg70-reg71; reg48=reg25+reg18; reg52=reg52*reg13; reg70=pow(reg61,2);
    reg74=reg56*reg37; T reg79=reg28*reg27; T reg80=reg55*reg37; T reg81=reg53*reg27; T reg82=reg68*reg6;
    T reg83=reg38*reg27; T reg84=reg57*reg37; reg15=reg15*reg19; reg75=reg75*reg5; reg41=reg39-reg41;
    reg16=reg16*reg12; reg51=reg51*reg13; reg11=reg11/reg36; reg39=reg66*reg6; reg63=reg64+reg63;
    reg59=reg60+reg59; reg60=reg66*reg26; reg64=reg34*reg35; T reg85=reg65*reg29; reg49=reg69+reg49;
    reg69=reg67*reg6; reg48=reg48*reg5; T reg86=reg65*reg27; T reg87=reg67*reg11; reg33=reg58+reg33;
    reg82=reg63+reg82; reg80=reg81+reg80; reg58=2*reg22; reg63=reg43*reg23; reg81=reg66*reg11;
    reg39=reg72+reg39; reg74=reg79+reg74; reg72=reg34*reg37; reg84=reg83+reg84; reg71=reg71*reg19;
    reg68=reg68*reg11; reg52=reg15-reg52; reg73=reg73*reg13; reg69=reg49+reg69; reg15=reg51+reg16;
    reg75=reg41-reg75; reg64=reg85+reg64; reg78=reg76+reg78; reg41=reg70*reg6; reg62=reg77+reg62;
    reg49=reg70*reg26; reg76=reg45*reg10; reg77=2*reg24; reg60=reg59+reg60; reg59=reg14*reg22;
    reg79=reg9*reg24; reg49=reg62+reg49; reg62=reg77*reg10; reg83=reg50*reg44; reg15=reg15*reg5;
    reg85=reg9*reg50; T reg88=reg24*reg44; reg73=reg71-reg73; reg48=reg52-reg48; reg52=2*reg50;
    reg54=reg54*reg61; reg68=reg84+reg68; reg71=reg70*reg11; reg72=reg86+reg72; reg87=reg74+reg87;
    reg74=reg69*reg59; reg84=reg78*reg79; reg86=reg53*reg33; T reg89=reg55*reg82; T reg90=reg55*reg69;
    T reg91=reg53*reg78; T reg92=reg55*reg39; T reg93=reg53*reg60; reg75=reg75/reg36; reg46=reg47+reg46;
    reg47=reg24*reg22; T reg94=reg9*reg14; reg41=reg64+reg41; reg64=reg33*reg79; T reg95=reg82*reg59;
    T reg96=reg63*reg4; T reg97=reg76*reg17; T reg98=reg58*reg23; reg81=reg80+reg81; reg33=reg28*reg33;
    reg69=reg56*reg69; reg78=reg28*reg78; reg82=reg56*reg82; reg80=reg56*reg39; T reg99=reg28*reg60;
    T reg100=reg63*reg35; T reg101=reg76*reg29; T reg102=reg77*reg22; reg60=reg79*reg60; reg39=reg39*reg59;
    T reg103=reg54*reg6; T reg104=reg66*reg87; reg29=reg62*reg29; reg74=reg84+reg74; reg84=reg53*reg49;
    reg69=reg78+reg69; reg78=reg67*reg87; reg89=reg86+reg89; reg86=reg66*reg68; T reg105=reg66*reg81;
    reg90=reg91+reg90; reg35=reg98*reg35; reg96=reg97+reg96; reg91=reg54*reg26; reg97=reg56*reg41;
    reg17=reg62*reg17; reg4=reg98*reg4; T reg106=reg94*reg75; T reg107=reg47*reg75; T reg108=reg45*reg14;
    T reg109=reg68*reg83; reg95=reg64+reg95; reg87=reg87*reg83; reg52=reg52*reg61; reg64=2*reg10;
    T reg110=reg55*reg41; reg71=reg72+reg71; reg48=reg48/reg36; reg72=reg76*reg27; T reg111=reg63*reg37;
    T reg112=reg14*reg50; T reg113=reg46*reg75; T reg114=reg22*reg44; reg80=reg99+reg80; reg99=reg10*reg23;
    T reg115=reg28*reg49; T reg116=reg9*reg44; T reg117=reg24*reg50; reg92=reg93+reg92; reg15=reg73-reg15;
    reg73=reg67*reg81; reg68=reg67*reg68; reg100=reg101+reg100; reg85=reg88+reg85; reg82=reg33+reg82;
    reg33=reg55*reg19; reg88=reg55*reg13; reg45=reg45*reg44; reg93=reg67*reg71; reg101=reg19*reg108;
    reg112=reg114+reg112; reg114=reg56*reg13; T reg118=reg13*reg102; T reg119=reg56*reg12; T reg120=reg13*reg108;
    T reg121=reg99*reg75; T reg122=reg12*reg102; reg86=reg89+reg86; reg104=reg90+reg104; reg89=reg48*reg85;
    reg90=reg107*reg108; T reg123=reg46*reg107; reg117=reg48*reg117; reg87=reg74+reg87; reg74=reg66*reg71;
    reg110=reg84+reg110; reg116=reg48*reg116; reg68=reg82+reg68; reg26=reg52*reg26; reg4=reg17+reg4;
    reg91=reg96+reg91; reg107=reg107*reg102; reg105=reg92+reg105; reg78=reg69+reg78; reg103=reg100+reg103;
    reg17=reg106*reg102; reg73=reg80+reg73; reg36=reg15/reg36; reg15=reg10*reg61; reg37=reg98*reg37;
    reg27=reg62*reg27; reg69=reg113*reg108; reg54=reg54*reg11; reg111=reg72+reg111; reg77=reg77*reg50;
    reg64=reg23*reg64; reg72=reg113*reg102; reg97=reg115+reg97; reg49=reg49*reg79; reg41=reg41*reg59;
    reg81=reg81*reg83; reg39=reg60+reg39; reg60=reg22*reg50; reg80=reg14*reg44; reg82=reg106*reg108;
    reg113=reg46*reg113; reg109=reg95+reg109; reg6=reg52*reg6; reg35=reg29+reg35; reg29=reg55*(*f.m).alpha_2;
    reg84=(*f.m).alpha_1*reg28; reg58=reg58*reg50; reg92=reg116*reg45; reg106=reg46*reg106; reg81=reg39+reg81;
    reg39=reg117*reg45; reg54=reg111+reg54; reg90=reg104+reg90; reg11=reg52*reg11; reg43=reg43*reg44;
    reg6=reg35+reg6; reg35=reg56*(*f.m).alpha_2; reg52=reg53*reg91; reg37=reg27+reg37; reg82=reg105+reg82;
    reg27=reg55*reg103; reg114=reg33-reg114; reg33=reg34*reg5; reg41=reg49+reg41; reg71=reg71*reg83;
    reg49=reg23*reg61; reg93=reg97+reg93; reg26=reg4+reg26; reg113=reg109+reg113; reg4=reg89*reg85;
    reg95=(*f.m).alpha_1*reg53; reg96=reg28*reg91; reg97=reg56*reg103; reg100=reg8*reg102; reg104=reg5*reg108;
    reg105=reg117*reg77; reg107=reg78+reg107; reg78=reg56*reg8; reg109=reg55*reg5; reg69=reg86+reg69;
    reg17=reg73+reg17; reg73=reg116*reg77; reg74=reg110+reg74; reg86=reg121*reg108; reg123=reg87+reg123;
    reg15=reg48*reg15; reg80=reg36*reg80; reg60=reg36*reg60; reg117=reg117*reg85; reg87=reg89*reg45;
    reg72=reg68+reg72; reg120=reg122-reg120; reg89=reg89*reg77; reg68=reg8*reg64; reg110=reg36*reg112;
    reg88=reg119-reg88; reg111=reg34*reg8; reg115=reg10*reg44; reg119=reg121*reg102; reg122=reg9*reg61;
    reg42=reg40+reg42; reg118=reg101-reg118; reg40=reg5*reg64; reg122=reg115+reg122; reg101=reg66*(*f.m).alpha_3;
    reg73=reg17+reg73; reg17=reg80*reg58; reg105=reg107+reg105; reg107=reg60*reg58; reg35=reg84+reg35;
    reg119=reg93+reg119; reg84=reg15*reg77; reg93=reg67*(*f.m).alpha_3; reg115=(*f.m).alpha_1*reg65; T reg124=reg34*(*f.m).alpha_2;
    T reg125=reg42*reg75; reg49=reg36*reg49; T reg126=reg23*reg44; T reg127=reg14*reg61; reg19=reg53*reg19;
    reg31=reg32+reg31; reg11=reg37+reg11; reg32=reg24*reg61; reg100=reg104+reg100; reg64=reg7*reg64;
    reg78=reg109+reg78; reg37=reg34*reg7; reg104=reg10*reg50; reg109=reg56*reg6; T reg128=reg28*reg26;
    T reg129=reg67*reg54; reg97=reg96+reg97; reg96=reg110*reg58; reg89=reg72+reg89; reg29=reg95+reg29;
    reg72=reg80*reg43; reg39=reg90+reg39; reg90=reg60*reg43; reg86=reg74+reg86; reg74=reg15*reg45;
    reg68=reg120-reg68; reg111=reg88-reg111; reg116=reg116*reg85; reg106=reg81+reg106; reg81=reg55*reg6;
    reg88=reg53*reg26; reg95=reg66*reg54; reg27=reg52+reg27; reg52=reg110*reg43; reg87=reg69+reg87;
    reg69=reg28*reg13; reg12=reg28*reg12; reg33=reg114-reg33; reg40=reg118-reg40; reg13=reg53*reg13;
    reg110=reg110*reg112; reg4=reg113+reg4; reg91=reg91*reg79; reg103=reg103*reg59; reg121=reg46*reg121;
    reg71=reg41+reg71; reg92=reg82+reg92; reg117=reg123+reg117; reg60=reg60*reg112; reg67=reg67*reg11;
    reg17=reg73+reg17; reg41=reg94*reg40; reg73=reg47*reg68; reg82=reg125*reg108; reg101=reg29+reg101;
    reg95=reg27+reg95; reg14=reg14*reg23; reg32=reg104+reg32; reg52=reg87+reg52; reg93=reg35+reg93;
    reg54=reg54*reg83; reg103=reg91+reg103; reg124=reg115+reg124; reg70=reg70*(*f.m).alpha_3; reg27=(*f.m).alpha_1*reg79;
    reg29=reg59*(*f.m).alpha_2; reg6=reg6*reg59; reg35=reg65*reg5; reg26=reg26*reg79; reg69=reg19-reg69;
    reg110=reg4+reg110; reg80=reg80*reg112; reg4=reg23*reg50; reg60=reg117+reg60; reg127=reg126+reg127;
    reg19=reg22*reg61; reg78=reg37-reg78; reg37=reg28*reg8; reg5=reg53*reg5; reg109=reg128+reg109;
    reg100=reg64-reg100; reg121=reg71+reg121; reg64=reg125*reg102; reg129=reg97+reg129; reg116=reg106+reg116;
    reg15=reg15*reg85; reg96=reg89+reg96; reg71=reg53*reg33; reg66=reg66*reg11; reg81=reg88+reg81;
    reg87=reg49*reg58; reg84=reg119+reg84; reg88=reg28*reg111; reg9=reg9*reg10; reg107=reg105+reg107;
    reg89=reg94*reg33; reg91=reg47*reg111; reg97=reg28*reg68; reg104=reg53*reg40; reg68=reg56*reg68;
    reg40=reg55*reg40; reg72=reg92+reg72; reg75=reg31*reg75; reg90=reg39+reg90; reg74=reg86+reg74;
    reg122=reg48*reg122; reg39=reg49*reg43; reg33=reg55*reg33; reg111=reg56*reg111; reg13=reg12-reg13;
    reg8=reg65*reg8; reg12=reg14*(*f.m).alpha_2; reg86=(*f.m).alpha_1*reg9; reg37=reg5+reg37; reg7=reg65*reg7;
    reg32=reg48*reg32; reg5=reg83*(*f.m).alpha_3; reg39=reg74+reg39; reg48=reg107*reg93; reg74=reg122*reg77;
    reg64=reg129+reg64; reg8=reg13-reg8; reg13=reg17*reg101; reg19=reg4+reg19; reg29=reg27+reg29;
    reg4=reg90*reg110; reg111=reg33+reg111; reg27=reg90*reg93; reg33=reg34*reg78; reg87=reg84+reg87;
    reg91=reg89+reg91; reg84=reg65*reg100; reg54=reg103+reg54; reg125=reg46*reg125; reg97=reg104+reg97;
    reg23=reg22*reg23; reg10=reg24*reg10; reg44=reg61*reg44; reg127=reg36*reg127; reg83=reg11*reg83;
    reg49=reg49*reg112; reg11=reg65*reg78; reg88=reg71+reg88; reg15=reg121+reg15; reg22=reg34*reg100;
    reg68=reg40+reg68; reg35=reg69-reg35; reg80=reg116+reg80; reg24=reg46*reg2; reg59=reg59*reg2;
    reg6=reg26+reg6; reg26=reg110*reg107; reg40=reg60*reg96; reg70=reg124+reg70; reg69=reg72*reg101;
    reg82=reg95+reg82; reg71=reg122*reg45; reg66=reg81+reg66; reg108=reg75*reg108; reg100=reg99*reg100;
    reg73=reg41+reg73; reg102=reg75*reg102; reg67=reg109+reg67; reg41=reg52*reg60; reg78=reg99*reg78;
    reg37=reg7-reg37; reg102=reg67+reg102; reg5=reg29+reg5; reg45=reg32*reg45; reg75=reg46*reg75;
    reg108=reg66+reg108; reg7=reg56*reg8; reg29=reg90*reg96; reg77=reg32*reg77; reg19=reg36*reg19;
    reg14=reg14*reg1; reg36=reg42*reg1; reg66=reg53*reg35; reg67=reg28*reg8; reg41=reg4-reg41;
    reg22=reg68+reg22; reg4=reg57*reg24; reg68=reg55*reg35; reg49=reg15+reg49; reg11=reg88+reg11;
    reg15=reg38*reg59; reg83=reg6+reg83; reg122=reg122*reg85; reg125=reg54+reg125; reg6=reg23*(*f.m).alpha_2;
    reg54=(*f.m).alpha_1*reg10; reg44=reg44*(*f.m).alpha_3; reg12=reg86+reg12; reg50=reg61*reg50; reg40=reg26-reg40;
    reg78=reg91+reg78; reg26=reg60*reg93; reg61=reg80*reg101; reg39=reg39*reg70; reg71=reg82+reg71;
    reg27=reg69+reg27; reg74=reg64+reg74; reg64=reg127*reg58; reg84=reg97+reg84; reg100=reg73+reg100;
    reg69=reg38*reg24; reg48=reg13+reg48; reg13=reg127*reg43; reg73=reg46*reg59; reg24=reg46*reg24;
    reg33=reg111+reg33; reg87=reg87*reg70; reg81=reg52*reg107; reg59=reg57*reg59; reg13=reg71+reg13;
    reg70=reg49*reg70; reg49=reg17*reg41; reg50=reg50*(*f.m).alpha_3; reg26=reg61+reg26; reg61=reg52*reg5;
    reg39=reg27+reg39; reg15=reg11+reg15; reg11=reg76*reg14; reg6=reg54+reg6; reg127=reg127*reg112;
    reg122=reg125+reg122; reg27=reg42*reg14; reg73=reg78+reg73; reg35=reg94*reg35; reg44=reg12+reg44;
    reg8=reg47*reg8; reg81=reg29-reg81; reg12=reg72*reg40; reg24=reg100+reg24; reg85=reg32*reg85;
    reg29=reg42*reg36; reg43=reg19*reg43; reg87=reg48+reg87; reg32=reg96*reg5; reg45=reg108+reg45;
    reg75=reg83+reg75; reg64=reg74+reg64; reg7=reg68+reg7; reg14=reg63*reg14; reg2=reg79*reg2;
    reg48=reg76*reg36; reg77=reg102+reg77; reg69=reg84+reg69; reg34=reg34*reg37; reg58=reg19*reg58;
    reg23=reg23*reg0; reg54=reg31*reg0; reg67=reg66+reg67; reg65=reg65*reg37; reg59=reg33+reg59;
    reg4=reg22+reg4; reg36=reg63*reg36; reg48=reg69+reg48; reg22=reg62*reg54; reg33=reg62*reg23;
    reg34=reg7+reg34; reg7=reg57*reg2; reg66=reg98*reg23; reg14=reg59+reg14; reg50=reg6+reg50;
    reg58=reg77+reg58; reg61=reg39+reg61; reg85=reg75+reg85; reg112=reg19*reg112; reg29=reg24+reg29;
    reg6=reg31*reg54; reg8=reg35+reg8; reg13=reg13*reg44; reg37=reg99*reg37; reg19=reg80*reg96;
    reg43=reg45+reg43; reg24=reg110*reg17; reg32=reg87+reg32; reg64=reg64*reg44; reg54=reg98*reg54;
    reg1=reg9*reg1; reg9=reg80*reg81; reg127=reg122+reg127; reg65=reg67+reg65; reg35=reg38*reg2;
    reg49=reg12-reg49; reg23=reg31*reg23; reg70=reg26+reg70; reg36=reg4+reg36; reg4=reg72*reg110;
    reg11=reg15+reg11; reg12=reg110*reg5; reg27=reg73+reg27; reg15=reg52*reg80; reg23=reg27+reg23;
    reg7=reg34+reg7; reg112=reg85+reg112; reg6=reg29+reg6; reg63=reg63*reg1; reg58=reg58*reg50;
    reg26=reg52*reg17; reg27=reg72*reg96; reg9=reg49+reg9; reg29=reg90*reg80; reg15=reg4-reg15;
    reg12=reg70+reg12; reg44=reg127*reg44; reg76=reg76*reg1; reg35=reg65+reg35; reg0=reg10*reg0;
    reg64=reg32+reg64; reg4=reg72*reg60; reg43=reg43*reg50; reg13=reg61+reg13; reg10=reg80*reg107;
    reg19=reg24-reg19; reg54=reg36+reg54; reg33=reg11+reg33; reg37=reg8+reg37; reg2=reg46*reg2;
    reg22=reg48+reg22; reg66=reg14+reg66; reg8=reg60*reg17; reg11=elem.pos(2)[0]-elem.pos(0)[0]; reg10=reg8-reg10;
    reg41=reg41/reg9; reg8=elem.pos(1)[0]-elem.pos(0)[0]; reg40=reg40/reg9; reg64=reg58+reg64; reg98=reg98*reg0;
    reg63=reg7+reg63; reg19=reg19/reg9; reg7=reg23*reg22; reg14=reg23*reg54; reg24=reg6*reg66;
    reg32=reg6*reg33; reg62=reg62*reg0; reg76=reg35+reg76; reg43=reg13+reg43; reg2=reg37+reg2;
    reg1=reg42*reg1; reg13=elem.pos(2)[1]-elem.pos(0)[1]; reg34=elem.pos(1)[1]-elem.pos(0)[1]; reg35=reg90*reg17; reg26=reg27-reg26;
    reg27=reg72*reg107; reg44=reg12+reg44; reg112=reg50*reg112; reg29=reg4-reg29; reg15=reg15/reg9;
    reg112=reg44+reg112; reg62=reg76+reg62; reg14=reg24-reg14; reg10=reg10/reg9; reg29=reg29/reg9;
    reg81=reg81/reg9; reg4=reg8*reg13; reg35=reg27-reg35; reg40=reg40*reg43; reg7=reg32-reg7;
    reg1=reg2+reg1; reg0=reg31*reg0; reg2=reg11*reg34; reg12=reg54*reg33; reg41=reg64*reg41;
    reg24=reg22*reg66; reg26=reg26/reg9; reg98=reg63+reg98; reg15=reg64*reg15; reg19=reg19*reg43;
    reg9=reg35/reg9; reg24=reg12-reg24; reg12=reg62*reg14; reg29=reg64*reg29; reg43=reg10*reg43;
    reg19=reg15-reg19; reg26=reg112*reg26; reg41=reg40-reg41; reg81=reg112*reg81; reg0=reg1+reg0;
    reg2=reg4-reg2; reg1=reg7*reg98; reg4=reg0*reg24; reg1=reg12-reg1; reg10=vectors[0][indices[1]+0]-vectors[0][indices[0]+0];
    reg13=reg13/reg2; reg8=reg8/reg2; reg34=reg34/reg2; reg11=reg11/reg2; reg12=reg0*reg33;
    reg29=reg43-reg29; reg9=reg112*reg9; reg15=1-(*f.m).resolution; reg26=reg19-reg26; reg41=reg81+reg41;
    reg19=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg27=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; reg31=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg32=reg0*reg66; reg35=reg23*reg98;
    reg23=reg23*reg62; reg33=reg33*reg98; reg36=(*f.m).resolution*reg93; reg37=reg8*reg27; reg39=reg11*reg10;
    reg40=reg13*reg19; reg29=reg9+reg29; reg26=reg15*reg26; reg66=reg62*reg66; reg9=(*f.m).resolution*reg101;
    reg42=reg6*reg98; reg43=reg0*reg54; reg41=reg15*reg41; reg6=reg6*reg62; reg4=reg1+reg4;
    reg32=reg35-reg32; reg1=reg34*reg31; reg12=reg23-reg12; reg0=reg0*reg22; reg12=reg12/reg4;
    reg98=reg22*reg98; reg29=reg15*reg29; reg31=reg8*reg31; reg0=reg6-reg0; reg9=reg41+reg9;
    reg26=reg36+reg26; reg62=reg54*reg62; reg6=(*f.m).resolution*reg5; reg19=reg11*reg19; reg33=reg66-reg33;
    reg39=reg37-reg39; reg27=reg34*reg27; reg32=reg32/reg4; reg1=reg40-reg1; reg43=reg42-reg43;
    reg10=reg13*reg10; reg33=reg33/reg4; reg60=reg60*reg15; reg98=reg62-reg98; reg80=reg80*reg15;
    reg27=reg10-reg27; reg39=reg1+reg39; reg1=(*f.m).resolution*reg32; reg14=reg14/reg4; reg43=reg43/reg4;
    reg7=reg7/reg4; reg0=reg0/reg4; reg26=(*f.m).deltaT*reg26; reg9=(*f.m).deltaT*reg9; reg19=reg31-reg19;
    reg29=reg6+reg29; reg6=(*f.m).resolution*reg12; reg10=(*f.m).resolution*reg33; reg39=0.5*reg39; reg29=(*f.m).deltaT*reg29;
    reg24=reg24/reg4; reg22=(*f.m).resolution*reg14; reg23=(*f.m).resolution*reg7; reg27=reg27-reg9; reg72=reg72*reg15;
    reg90=reg90*reg15; reg17=reg15*reg17; reg107=reg15*reg107; reg110=reg110*reg15; reg1=reg80+reg1;
    reg6=reg60-reg6; reg31=(*f.m).resolution*reg43; reg35=(*f.m).resolution*reg0; reg4=reg98/reg4; reg19=reg19-reg26;
    reg10=reg110+reg10; reg36=reg1*reg27; reg37=reg6*reg19; reg31=reg17-reg31; reg39=reg39-reg29;
    reg17=(*f.m).resolution*reg24; reg35=reg107+reg35; reg22=reg72+reg22; reg96=reg15*reg96; reg23=reg90-reg23;
    reg40=(*f.m).resolution*reg4; reg15=reg52*reg15; reg17=reg15+reg17; reg15=reg31*reg27; reg41=reg35*reg19;
    reg42=reg39*reg10; reg37=reg36+reg37; reg27=reg22*reg27; reg40=reg96-reg40; reg19=reg23*reg19;
    reg36=reg39*reg17; reg19=reg27+reg19; reg27=reg11-reg8; reg37=reg42+reg37; reg42=reg34-reg13;
    reg39=reg39*reg40; reg41=reg15+reg41; reg19=reg36+reg19; reg39=reg41+reg39; reg15=0.5*reg8;
    reg36=0.5*reg13; reg41=0.5*reg34; reg44=0.5*reg11; reg45=0.5*reg42; reg37=2*reg37;
    reg48=0.5*reg27; reg49=1-var_inter[0]; reg50=reg41*reg37; reg52=reg39*reg8; reg49=reg49-var_inter[1];
    reg54=reg19*reg13; reg58=reg37*reg44; reg59=reg39*reg11; reg60=reg48*reg37; reg61=reg42*reg19;
    reg62=reg37*reg36; reg63=reg39*reg27; reg64=reg37*reg45; reg65=reg37*reg15; reg66=reg19*reg34;
    reg67=reg49*elem.f_vol_e[0]; reg60=reg61+reg60; reg61=var_inter[1]*elem.f_vol_e[1]; reg68=var_inter[1]*elem.f_vol_e[0]; reg52=reg52-reg50;
    reg65=reg65-reg66; reg69=var_inter[0]*elem.f_vol_e[1]; reg64=reg63+reg64; reg63=elem.f_vol_e[1]*reg49; reg70=var_inter[0]*elem.f_vol_e[0];
    reg54=reg54-reg58; reg62=reg62-reg59; reg62=reg62-reg69; reg54=reg54-reg70; reg64=reg64-reg63;
    reg65=reg65-reg68; reg52=reg52-reg61; reg60=reg60-reg67; reg54=reg54*reg2; reg52=reg52*reg2;
    reg64=reg64*reg2; reg65=reg65*reg2; reg60=reg60*reg2; reg62=reg62*reg2; T vec_2=ponderation*reg54;
    T vec_0=ponderation*reg60; T vec_1=ponderation*reg64; T vec_4=ponderation*reg65; T vec_5=ponderation*reg52; T vec_3=ponderation*reg62;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
sollicitation[indices[0]+0] += vec_0;
sollicitation[indices[0]+1] += vec_1;
sollicitation[indices[1]+0] += vec_2;
sollicitation[indices[1]+1] += vec_3;
sollicitation[indices[2]+0] += vec_4;
sollicitation[indices[2]+1] += vec_5;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
#ifndef ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_false_true_true
#define ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_false_true_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const typename TM::TNode &node,
      const unsigned *indices){ 
  #define PNODE(N) node
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
#endif
// 
#ifndef ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_true_true_true
#define ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_true_true_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const typename TM::TNode &node,
      const unsigned *indices){ 
  #define PNODE(N) node
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
#endif
// 
#ifndef ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_false_true_false
#define ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_false_true_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const typename TM::TNode &node,
      const unsigned *indices){ 
  #define PNODE(N) node
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
#endif
// 
#ifndef ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_true_true_false
#define ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_true_true_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const typename TM::TNode &node,
      const unsigned *indices){ 
  #define PNODE(N) node
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
#endif
// 
#ifndef ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_symmetric_version_false_true
#define ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_symmetric_version_false_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE ,unsigned symmetric_version>
void add_nodal_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const typename TM::TNode &node,
      const unsigned *indices){ 
  #define PNODE(N) node
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
#endif
// 
#ifndef ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_symmetric_version_false_false
#define ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_symmetric_version_false_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE ,unsigned symmetric_version>
void add_nodal_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const typename TM::TNode &node,
      const unsigned *indices){ 
  #define PNODE(N) node
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
#endif
// 
#ifndef ADD_NODAL_RESIDUAL_elasticity_orthotropy_stat_Qstat
#define ADD_NODAL_RESIDUAL_elasticity_orthotropy_stat_Qstat
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE>
void add_nodal_residual(
      const Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const typename TM::TNode &node,
      const unsigned *indices ) { 
  #define PNODE(N) node
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
#endif
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}

#ifndef elasticity_orthotropy_stat_Qstat_read_material_to_mesh
#define elasticity_orthotropy_stat_Qstat_read_material_to_mesh
template<class TM, class T, bool wont_add_nz>
void read_material_to_mesh_(const XmlNode &n, Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f){ 
    if(n.has_attribute("elastic_modulus_1"))  
        n.get_attribute("elastic_modulus_1", f.m->elastic_modulus_1 ); 
    else  
        std::cerr << "Warning using default value of elastic_modulus_1 : " << f.m->elastic_modulus_1 << std::endl; 

    if(n.has_attribute("alpha_2"))  
        n.get_attribute("alpha_2", f.m->alpha_2 ); 
    else  
        std::cerr << "Warning using default value of alpha_2 : " << f.m->alpha_2 << std::endl; 

    if(n.has_attribute("alpha_3"))  
        n.get_attribute("alpha_3", f.m->alpha_3 ); 
    else  
        std::cerr << "Warning using default value of alpha_3 : " << f.m->alpha_3 << std::endl; 

    if(n.has_attribute("alpha_1"))  
        n.get_attribute("alpha_1", f.m->alpha_1 ); 
    else  
        std::cerr << "Warning using default value of alpha_1 : " << f.m->alpha_1 << std::endl; 

    if(n.has_attribute("f_vol"))  
        n.get_attribute("f_vol", f.m->f_vol ); 
    else  
        std::cerr << "Warning using default value of f_vol : " << f.m->f_vol << std::endl; 

    if(n.has_attribute("shear_modulus_13"))  
        n.get_attribute("shear_modulus_13", f.m->shear_modulus_13 ); 
    else  
        std::cerr << "Warning using default value of shear_modulus_13 : " << f.m->shear_modulus_13 << std::endl; 

    if(n.has_attribute("shear_modulus_12"))  
        n.get_attribute("shear_modulus_12", f.m->shear_modulus_12 ); 
    else  
        std::cerr << "Warning using default value of shear_modulus_12 : " << f.m->shear_modulus_12 << std::endl; 

    if(n.has_attribute("deltaT"))  
        n.get_attribute("deltaT", f.m->deltaT ); 
    else  
        std::cerr << "Warning using default value of deltaT : " << f.m->deltaT << std::endl; 

    if(n.has_attribute("density"))  
        n.get_attribute("density", f.m->density ); 
    else  
        std::cerr << "Warning using default value of density : " << f.m->density << std::endl; 

    if(n.has_attribute("poisson_ratio_13"))  
        n.get_attribute("poisson_ratio_13", f.m->poisson_ratio_13 ); 
    else  
        std::cerr << "Warning using default value of poisson_ratio_13 : " << f.m->poisson_ratio_13 << std::endl; 

    if(n.has_attribute("poisson_ratio_12"))  
        n.get_attribute("poisson_ratio_12", f.m->poisson_ratio_12 ); 
    else  
        std::cerr << "Warning using default value of poisson_ratio_12 : " << f.m->poisson_ratio_12 << std::endl; 

    if(n.has_attribute("v2"))  
        n.get_attribute("v2", f.m->v2 ); 
    else  
        std::cerr << "Warning using default value of v2 : " << f.m->v2 << std::endl; 

    if(n.has_attribute("elastic_modulus_3"))  
        n.get_attribute("elastic_modulus_3", f.m->elastic_modulus_3 ); 
    else  
        std::cerr << "Warning using default value of elastic_modulus_3 : " << f.m->elastic_modulus_3 << std::endl; 

    if(n.has_attribute("elastic_modulus_2"))  
        n.get_attribute("elastic_modulus_2", f.m->elastic_modulus_2 ); 
    else  
        std::cerr << "Warning using default value of elastic_modulus_2 : " << f.m->elastic_modulus_2 << std::endl; 

    if(n.has_attribute("shear_modulus_23"))  
        n.get_attribute("shear_modulus_23", f.m->shear_modulus_23 ); 
    else  
        std::cerr << "Warning using default value of shear_modulus_23 : " << f.m->shear_modulus_23 << std::endl; 

    if(n.has_attribute("v1"))  
        n.get_attribute("v1", f.m->v1 ); 
    else  
        std::cerr << "Warning using default value of v1 : " << f.m->v1 << std::endl; 

    if(n.has_attribute("resolution"))  
        n.get_attribute("resolution", f.m->resolution ); 
    else  
        std::cerr << "Warning using default value of resolution : " << f.m->resolution << std::endl; 

    if(n.has_attribute("poisson_ratio_23"))  
        n.get_attribute("poisson_ratio_23", f.m->poisson_ratio_23 ); 
    else  
        std::cerr << "Warning using default value of poisson_ratio_23 : " << f.m->poisson_ratio_23 << std::endl; 

  };
#endif // elasticity_orthotropy_stat_Qstat_read_material_to_mesh
} // namespace LMT


#include "formulation/formulation.h"
namespace LMT {
#ifndef ELASTICITY_ORTHOTROPY_STAT_QSTAT_2_NUM_0_CARAC_H
#define ELASTICITY_ORTHOTROPY_STAT_QSTAT_2_NUM_0_CARAC_H
#ifndef STRUCT_ELASTICITY_ORTHOTROPY_STAT_QSTAT
#define STRUCT_ELASTICITY_ORTHOTROPY_STAT_QSTAT
struct elasticity_orthotropy_stat_Qstat {
  static const char *name() { return "elasticity_orthotropy_stat_Qstat"; }
};
#endif // ELASTICITY_ORTHOTROPY_STAT_QSTAT

template<class P_T>
class CaracFormulation<elasticity_orthotropy_stat_Qstat,2,P_T>  {
public:
  typedef P_T T;
  static const char *name() { return "elasticity_orthotropy_stat_Qstat"; }
  static const bool matrix_will_be_definite_positive=true;
  static const bool has_nodal_matrix = false;
  static const bool has_IS_contact_matrix=false;
  static const bool need_skin_assembly=false;
  typedef Norm1_is_inf Name_convergence_criterium;
  static const unsigned nb_vectors = 4;
  static const unsigned nb_matrices = 4;
  static const unsigned auto_contact = false;
  static const bool friction_coeff_is_a_nodal_variable = 0;
  static const unsigned offset_of_pos_unknown=3;
  static const unsigned pos_is_an_unknown = false;
  static const unsigned nb_der_var = 0;
  template<class TF> static void add_to_der_vars( TF &f, const Vec<T> &v ) {
  }
  static bool is_unknown(const std::string &s) { return (s=="dep"); }
  static unsigned num_in_vec_unknown(const std::string &s) { if ( s=="dep" )return 0; return 0;  }
  template<unsigned num_mat,unsigned inner=0> struct NodalMatricesCarac {
      static const bool symm = 1;
      static const bool herm = false;
      static const bool diag = false;
  };
  template<unsigned num_mat,unsigned inner=0> struct GlobalMatricesCarac {
      static const bool symm = 1;
      static const bool herm = false;
      static const bool diag = false;
  };
  
  static const unsigned nb_nodal_unknowns = 2;
  template<class TE,class TTs,class Tvec>
  inline static void set_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice) {
    node.dep[0]=vecs[0][indice+0]; node.dep[1]=vecs[0][indice+1];
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
    node.dep[0]=vecs[0][indice+0]; node.dep[1]=vecs[0][indice+1];
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_nodal_initial_conditions(const TE &node,const TTs &f,Tvec &vecs,unsigned indice) {
    vecs[0][indice+0]=node.dep[0]; vecs[1][indice+0]=node.dep[0]; vecs[2][indice+0]=node.dep[0]; vecs[3][indice+0]=node.dep[0]; vecs[4][indice+0]=node.dep[0];
    vecs[0][indice+1]=node.dep[1]; vecs[1][indice+1]=node.dep[1]; vecs[2][indice+1]=node.dep[1]; vecs[3][indice+1]=node.dep[1]; vecs[4][indice+1]=node.dep[1];
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_nodal_error(const TE &node,const TTs &f,const Tvec &vecs,int indice) {
    T reg0=vecs[1][indice+1]-vecs[0][indice+1]; T reg1=vecs[1][indice+0]-vecs[0][indice+0]; reg0=abs(reg0); reg1=abs(reg1); return max(reg0,reg1);
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_nodal(const TE &node,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
    old_vec[indice+1]=vecs[1][indice+1]; old_vec[indice+0]=vecs[1][indice+0];
  }
  
  static const unsigned nb_global_unknowns = 0;
  template<class TE,class TTs,class Tvec>
  inline static void set_global_unknowns(TE &mesh,const TTs &f,const Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_global_unknowns(TE &mesh,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_global_initial_conditions(const TE &mesh,const TTs &f,Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_global_error(const TE &mesh,const TTs &f,const Tvec &vecs,int indice) {
    return 0;
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_global(const TE &mesh,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
  }
  
  static const unsigned nb_skin_nodal_unknowns = 0;
  template<class TE,class TTs,class Tvec>
  inline static void set_skin_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_skin_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_skin_nodal_initial_conditions(const TE &node,const TTs &f,Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_skin_nodal_error(const TE &node,const TTs &f,const Tvec &vecs,int indice) {
    return 0;
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_skin_nodal(const TE &node,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
  }
  
  static const unsigned nb_skin_global_unknowns = 0;
  template<class TE,class TTs,class Tvec>
  inline static void set_skin_global_unknowns(TE &mesh,const TTs &f,const Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_skin_global_unknowns(TE &mesh,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_skin_global_initial_conditions(const TE &mesh,const TTs &f,Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_skin_global_error(const TE &mesh,const TTs &f,const Tvec &vecs,int indice) {
    return 0;
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_skin_global(const TE &mesh,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
  }
};
#endif // ELASTICITY_ORTHOTROPY_STAT_QSTAT_2_NUM_0_CARAC_H
extern "C" void apply_on_elements_after_solve_0_elasticity_orthotropy_stat_Qstat_Quad_0( double * );
extern "C" void apply_on_elements_after_solve_1_elasticity_orthotropy_stat_Qstat_Quad_1( double * );
extern "C" void apply_on_elements_after_solve_2_elasticity_orthotropy_stat_Qstat_Quad_2( double * );
extern "C" void apply_on_elements_after_solve_3_elasticity_orthotropy_stat_Qstat_Quad_3( double * );
extern "C" void apply_on_elements_after_solve_4_elasticity_orthotropy_stat_Qstat_Quad_4( double * );
extern "C" void apply_on_elements_after_solve_5_elasticity_orthotropy_stat_Qstat_Quad_5( double * );
extern "C" void apply_on_elements_after_solve_6_elasticity_orthotropy_stat_Qstat_Quad_6( double * );
extern "C" void apply_on_elements_after_solve_7_elasticity_orthotropy_stat_Qstat_Quad_7( double * );
extern "C" void apply_on_elements_after_solve_8_elasticity_orthotropy_stat_Qstat_Quad_8( double * );
extern "C" void apply_on_elements_after_solve_9_elasticity_orthotropy_stat_Qstat_Quad_9( double * );
extern "C" void apply_on_elements_after_solve_10_elasticity_orthotropy_stat_Qstat_Quad_10( double * );
extern "C" void apply_on_elements_after_solve_11_elasticity_orthotropy_stat_Qstat_Quad_11( double * );
extern "C" void apply_on_elements_after_solve_12_elasticity_orthotropy_stat_Qstat_Quad_12( double * );
extern "C" void apply_on_elements_after_solve_13_elasticity_orthotropy_stat_Qstat_Quad_13( double * );
extern "C" void apply_on_elements_after_solve_14_elasticity_orthotropy_stat_Qstat_Quad_14( double * );
class Quad;
template<unsigned A,class B,class C> class Node;
template<class A,class B,class C,class D,unsigned E> class Element;

// Carac for ...
template<class P_T_pos,class P_ND,class TED,unsigned nim,class TM,class T>
class CaracFormulationForElement<elasticity_orthotropy_stat_Qstat,Element<Quad,DefaultBehavior,Node<2,P_T_pos,P_ND>,TED,nim>,TM,T> {
public:
    template<unsigned num_mat,unsigned inner=0> struct ElemMatricesCarac {
        static const bool symm = true;
        static const bool herm = false;
        static const bool diag = false;
        static const bool linear = true;
    };
    static const unsigned order_integration = 2;
    static const bool has_elementary_matrix = true;
    static const bool has_skin_elementary_matrix = false;
    template<class TE,class TF, class TVEVE> static void after_solve(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
    T reg0=2*(*f.m).shear_modulus_13; T reg1=2*(*f.m).shear_modulus_23; T reg2=2*(*f.m).shear_modulus_12; reg0=1.0/reg0; reg1=1.0/reg1;
    reg2=1.0/reg2; T reg3=reg0*reg1; T reg4=pow((*f.m).v2[0],2); T reg5=pow((*f.m).v2[1],2); T reg6=pow((*f.m).v1[1],2);
    T reg7=pow((*f.m).v1[0],2); T reg8=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=reg2*reg3; T reg11=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg12=1.0/(*f.m).elastic_modulus_2; T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; reg6=reg7+reg6; reg7=pow((*f.m).v1[2],2); T reg14=reg8*reg10;
    T reg15=reg11*reg10; reg5=reg4+reg5; reg4=pow((*f.m).v2[2],2); T reg16=reg9*reg10; reg7=reg6+reg7;
    reg4=reg5+reg4; reg5=reg12*reg16; reg6=reg11*reg14; T reg17=reg11*reg15; T reg18=reg13*reg16;
    T reg19=reg12*reg14; reg6=reg18+reg6; T reg20=reg13*reg15; T reg21=1.0/(*f.m).elastic_modulus_1; reg17=reg5-reg17;
    reg4=pow(reg4,0.5); reg7=pow(reg7,0.5); reg5=reg13*reg6; T reg22=reg21*reg17; T reg23=reg20+reg19;
    T reg24=(*f.m).v2[2]/reg4; T reg25=(*f.m).v2[1]/reg4; T reg26=(*f.m).v1[2]/reg7; T reg27=(*f.m).v1[1]/reg7; T reg28=reg8*reg23;
    reg5=reg22-reg5; reg22=reg12*reg10; T reg29=reg8*reg15; reg16=reg21*reg16; T reg30=reg8*reg14;
    reg10=reg13*reg10; reg7=(*f.m).v1[0]/reg7; reg4=(*f.m).v2[0]/reg4; T reg31=reg8*reg3; T reg32=reg11*reg3;
    T reg33=reg26*reg25; T reg34=reg27*reg24; reg3=reg9*reg3; T reg35=reg2*reg1; T reg36=reg11*reg32;
    T reg37=reg11*reg31; reg28=reg5-reg28; reg5=reg34-reg33; T reg38=2*reg7; T reg39=reg13*reg3;
    T reg40=reg26*reg4; T reg41=reg7*reg24; T reg42=reg9*reg35; reg3=reg12*reg3; T reg43=reg8*reg35;
    T reg44=2*reg4; T reg45=reg8*reg10; reg15=reg21*reg15; reg30=reg16-reg30; reg35=reg11*reg35;
    reg16=reg2*reg0; T reg46=reg8*reg22; reg14=reg13*reg14; reg29=reg18+reg29; reg31=reg12*reg31;
    reg45=reg15+reg45; reg18=reg8*reg16; T reg47=reg11*reg16; reg30=reg30/reg28; reg36=reg3-reg36;
    reg14=reg15+reg14; reg6=reg6/reg28; reg46=reg20+reg46; reg37=reg39+reg37; reg29=reg29/reg28;
    reg10=reg13*reg10; reg3=reg12*reg42; reg42=reg13*reg42; reg32=reg13*reg32; reg15=reg11*reg35;
    reg39=reg11*reg43; reg16=reg9*reg16; T reg48=2*reg5; T reg49=reg44*reg25; T reg50=reg27*reg4;
    T reg51=reg7*reg25; T reg52=reg40-reg41; T reg53=pow(reg7,2); T reg54=pow(reg27,2); T reg55=reg27*reg38;
    T reg56=pow(reg4,2); T reg57=pow(reg25,2); reg17=reg17/reg28; reg22=reg21*reg22; T reg58=reg17*reg54;
    T reg59=reg6*reg57; T reg60=reg17*reg55; T reg61=reg6*reg49; T reg62=reg29*reg53; T reg63=reg30*reg56;
    T reg64=reg29*reg54; T reg65=reg48*reg52; T reg66=reg30*reg57; T reg67=reg29*reg55; T reg68=reg30*reg49;
    T reg69=pow(reg52,2); T reg70=pow(reg5,2); T reg71=reg11*reg18; T reg72=reg11*reg47; T reg73=reg13*reg16;
    reg16=reg12*reg16; reg39=reg42+reg39; reg15=reg3-reg15; reg3=reg51-reg50; reg42=reg32+reg31;
    reg37=reg13*reg37; reg36=reg21*reg36; T reg74=pow(reg26,2); reg43=reg12*reg43; reg35=reg13*reg35;
    T reg75=pow(reg24,2); reg45=reg45/reg28; T reg76=reg17*reg53; T reg77=reg6*reg56; reg46=reg46/reg28;
    reg10=reg22-reg10; reg14=reg14/reg28; reg23=reg23/reg28; reg22=reg14*reg65; reg68=reg67+reg68;
    reg10=reg10/reg28; reg67=reg46*reg54; reg72=reg16-reg72; reg16=reg30*reg75; T reg78=reg29*reg74;
    T reg79=reg14*reg69; reg66=reg64+reg66; reg71=reg73+reg71; reg64=reg14*reg70; reg63=reg62+reg63;
    reg62=reg45*reg57; reg59=reg58+reg59; reg58=reg23*reg69; reg73=reg23*reg65; reg61=reg60+reg61;
    reg60=reg17*reg74; T reg80=reg6*reg75; T reg81=reg46*reg55; T reg82=reg45*reg49; T reg83=reg23*reg70;
    reg77=reg76+reg77; reg47=reg13*reg47; reg18=reg12*reg18; reg37=reg36-reg37; reg42=reg8*reg42;
    reg36=pow(reg3,2); reg15=reg21*reg15; reg76=reg35+reg43; reg39=reg13*reg39; T reg84=reg45*reg56;
    T reg85=reg46*reg53; reg80=reg60+reg80; reg76=reg8*reg76; reg83=reg77+reg83; reg84=reg85+reg84;
    reg60=reg10*reg70; reg58=reg59+reg58; reg59=2*reg27; reg77=reg26*reg38; reg72=reg21*reg72;
    reg85=reg47+reg18; reg39=reg15-reg39; reg42=reg37-reg42; reg15=reg4*reg25; reg37=reg7*reg27;
    reg71=reg13*reg71; reg82=reg81+reg82; reg81=reg10*reg65; T reg86=2*reg25; T reg87=reg44*reg24;
    reg22=reg68+reg22; reg68=reg45*reg75; T reg88=reg14*reg36; reg16=reg78+reg16; reg79=reg66+reg79;
    reg66=reg46*reg74; reg78=reg23*reg36; reg64=reg63+reg64; reg63=reg10*reg69; reg62=reg67+reg62;
    reg73=reg61+reg73; reg81=reg82+reg81; reg61=reg15*reg79; reg67=reg37*reg58; reg82=reg22*reg57;
    reg48=reg48*reg3; reg42=reg42/reg28; T reg89=reg73*reg54; T reg90=reg79*reg57; T reg91=2*reg52;
    T reg92=reg58*reg54; T reg93=reg7*reg4; T reg94=reg27*reg25; reg51=reg50+reg51; reg50=reg64*reg57;
    T reg95=reg27*reg5; T reg96=reg26*reg59; T reg97=reg83*reg54; T reg98=reg7*reg52; reg78=reg80+reg78;
    reg85=reg8*reg85; reg63=reg62+reg63; reg71=reg72-reg71; reg88=reg16+reg88; reg16=reg6*reg87;
    reg60=reg84+reg60; reg62=reg58*reg53; reg72=reg79*reg56; reg76=reg39-reg76; reg39=reg17*reg77;
    reg80=reg29*reg77; reg84=reg30*reg87; T reg99=reg64*reg56; T reg100=reg73*reg53; T reg101=reg83*reg53;
    T reg102=reg22*reg56; T reg103=reg86*reg24; reg68=reg66+reg68; reg66=reg15*reg22; T reg104=reg10*reg36;
    T reg105=reg37*reg73; T reg106=reg52*reg5; reg76=reg76/reg28; T reg107=reg106*reg63; reg61=reg67+reg61;
    reg67=2*reg26; T reg108=reg25*reg59; reg50=reg97+reg50; reg97=reg69*reg60; reg90=reg92+reg90;
    reg92=reg69*reg63; T reg109=reg4*reg38; T reg110=reg78*reg54; T reg111=reg88*reg57; reg104=reg68+reg104;
    reg68=reg46*reg77; T reg112=reg45*reg87; T reg113=reg42*reg93; T reg114=reg42*reg94; T reg115=reg42*reg51;
    reg82=reg89+reg82; reg99=reg101+reg99; reg89=reg70*reg60; reg101=reg69*reg81; reg72=reg62+reg72;
    reg62=reg70*reg63; T reg116=reg78*reg53; T reg117=reg88*reg56; reg102=reg100+reg102; reg100=reg70*reg81;
    T reg118=reg15*reg64; T reg119=reg37*reg83; reg30=reg30*reg103; reg29=reg29*reg96; T reg120=reg4*reg52;
    T reg121=reg14*reg48; reg84=reg80+reg84; reg80=reg25*reg5; reg6=reg6*reg103; reg17=reg17*reg96;
    reg91=reg91*reg3; T reg122=reg23*reg48; T reg123=reg26*reg24; reg16=reg39+reg16; reg39=reg7*reg5;
    T reg124=reg27*reg52; T reg125=reg106*reg81; reg98=reg95+reg98; reg66=reg105+reg66; reg85=reg71-reg85;
    reg71=reg76*reg39; reg95=reg76*reg124; reg105=reg26*reg3; T reg126=reg76*reg98; reg14=reg14*reg91;
    reg120=reg80+reg120; reg118=reg119+reg118; reg89=reg99+reg89; reg80=reg106*reg60; reg99=reg109*reg115;
    reg101=reg82+reg101; reg82=reg108*reg115; reg119=reg109*reg113; T reg127=reg25*reg52; T reg128=reg4*reg5;
    reg100=reg102+reg100; reg62=reg72+reg62; reg72=reg109*reg114; reg102=reg70*reg104; reg117=reg116+reg117;
    reg116=reg15*reg88; reg28=reg85/reg28; reg85=reg37*reg78; T reg129=reg51*reg114; reg107=reg61+reg107;
    reg61=reg67*reg24; reg125=reg66+reg125; reg59=reg52*reg59; reg66=reg51*reg115; reg38=reg5*reg38;
    reg97=reg50+reg97; reg50=reg108*reg113; reg122=reg16+reg122; reg92=reg90+reg92; reg16=reg108*reg114;
    reg6=reg17+reg6; reg23=reg23*reg91; reg111=reg110+reg111; reg17=reg69*reg104; reg121=reg84+reg121;
    reg112=reg68+reg112; reg68=reg10*reg48; reg46=reg46*reg96; reg45=reg45*reg103; reg84=reg42*reg123;
    reg30=reg29+reg30; reg29=reg13*reg57; reg90=reg21*reg56; reg110=reg108*reg12; T reg130=reg13*reg56;
    T reg131=reg109*reg13; T reg132=reg109*reg21; T reg133=reg108*reg13; T reg134=reg12*reg57; T reg135=reg108*reg84;
    T reg136=reg108*reg11; reg17=reg111+reg17; reg102=reg117+reg102; reg111=reg109*reg84; reg133=reg132-reg133;
    reg117=reg61*reg8; reg132=reg59*reg95; reg16=reg92+reg16; reg99=reg100+reg99; reg92=reg38*reg126;
    reg29=reg90-reg29; reg90=reg122*reg54; reg100=reg122*reg53; T reg137=reg121*reg56; T reg138=reg59*reg71;
    T reg139=reg28*reg120; T reg140=reg28*reg127; T reg141=reg28*reg128; T reg142=reg76*reg105; reg131=reg110-reg131;
    reg110=reg61*reg11; T reg143=(*f.m).alpha_2*reg56; reg130=reg134-reg130; reg134=reg8*reg56; reg10=reg10*reg91;
    reg119=reg89+reg119; reg89=reg38*reg71; reg45=reg46+reg45; reg46=reg11*reg75; reg82=reg101+reg82;
    reg101=reg59*reg126; T reg144=reg11*reg57; reg68=reg112+reg68; reg112=reg109*reg8; reg72=reg62+reg72;
    reg62=reg38*reg95; T reg145=reg26*reg5; T reg146=(*f.m).alpha_1*reg53; reg40=reg41+reg40; reg41=reg7*reg3;
    T reg147=reg24*reg3; reg14=reg30+reg14; reg23=reg6+reg23; reg6=reg98*reg126; reg66=reg125+reg66;
    reg30=reg106*reg104; reg116=reg85+reg116; reg85=reg98*reg95; reg129=reg107+reg129; reg107=(*f.m).alpha_2*reg57;
    reg86=reg86*reg52; reg44=reg44*reg5; reg125=(*f.m).alpha_1*reg54; reg80=reg118+reg80; reg118=reg51*reg113;
    T reg148=reg8*reg75; T reg149=reg121*reg57; reg50=reg97+reg50; reg97=reg37*reg122; T reg150=reg15*reg121;
    reg143=reg146+reg143; reg146=reg86*reg141; T reg151=reg4*reg3; reg34=reg33+reg34; reg10=reg45+reg10;
    reg33=reg24*reg5; reg46=reg130-reg46; reg110=reg131-reg110; reg45=reg27*reg3; reg130=reg26*reg52;
    reg131=reg42*reg40; reg41=reg145+reg41; reg138=reg50+reg138; reg50=reg98*reg71; reg132=reg16+reg132;
    reg16=reg86*reg140; reg107=reg125+reg107; reg118=reg80+reg118; reg117=reg133-reg117; reg85=reg129+reg85;
    reg80=reg13*reg53; reg125=reg120*reg140; reg129=reg69*(*f.m).alpha_3; reg133=(*f.m).alpha_1*reg74; reg135=reg17+reg135;
    reg17=reg59*reg142; reg145=reg12*reg54; reg148=reg29-reg148; reg30=reg116+reg30; reg29=reg51*reg84;
    reg116=(*f.m).alpha_2*reg75; T reg152=reg70*(*f.m).alpha_3; reg6=reg66+reg6; reg66=reg120*reg139; reg111=reg102+reg111;
    reg136=reg112+reg136; reg102=reg61*reg9; reg112=reg38*reg142; T reg153=reg23*reg53; T reg154=reg44*reg140;
    T reg155=reg14*reg56; reg101=reg82+reg101; reg82=reg86*reg139; T reg156=reg23*reg54; T reg157=reg14*reg57;
    reg62=reg72+reg62; reg144=reg134+reg144; reg72=reg9*reg75; reg134=reg44*reg141; reg89=reg119+reg89;
    reg119=reg28*reg147; T reg158=reg13*reg54; reg137=reg100+reg137; reg149=reg90+reg149; reg90=reg70*reg68;
    reg100=reg21*reg53; T reg159=reg44*reg139; reg92=reg99+reg92; reg99=reg69*reg68; reg90=reg137+reg90;
    reg137=reg24*reg52; reg151=reg33+reg151; reg45=reg130+reg45; reg33=reg25*reg3; reg130=reg117*reg53;
    reg152=reg143+reg152; reg50=reg118+reg50; reg99=reg149+reg99; reg159=reg92+reg159; reg92=reg120*reg141;
    reg82=reg101+reg82; reg101=reg108*reg131; reg146=reg138+reg146; reg118=reg109*reg131; reg138=reg37*reg23;
    reg143=reg148*reg53; reg149=reg15*reg14; T reg160=reg15*(*f.m).alpha_2; T reg161=reg106*reg68; reg150=reg97+reg150;
    reg97=reg46*reg54; reg66=reg6+reg66; reg6=reg37*(*f.m).alpha_1; T reg162=reg36*(*f.m).alpha_3; reg7=reg7*reg26;
    reg116=reg133+reg116; reg4=reg4*reg24; reg133=reg69*reg10; T reg163=reg98*reg142; reg29=reg30+reg29;
    reg157=reg156+reg157; reg129=reg107+reg129; reg30=reg70*reg10; reg125=reg85+reg125; reg155=reg153+reg155;
    reg85=reg110*reg54; reg107=reg148*reg56; reg153=reg46*reg57; reg156=reg44*reg119; reg16=reg132+reg16;
    reg112=reg111+reg112; reg111=reg11*reg54; reg132=reg148*reg93; T reg164=reg46*reg94; reg144=reg72-reg144;
    reg72=reg117*reg56; reg136=reg102-reg136; reg134=reg89+reg134; reg89=reg110*reg57; reg154=reg62+reg154;
    reg42=reg42*reg34; reg62=reg8*reg53; reg102=reg76*reg41; T reg165=reg8*reg74; T reg166=reg110*reg94;
    T reg167=reg117*reg93; reg158=reg100-reg158; reg100=reg11*reg74; T reg168=reg86*reg119; reg80=reg145-reg80;
    reg17=reg135+reg17; reg163=reg29+reg163; reg29=reg125*reg159; reg135=reg136*reg75; reg89=reg72+reg89;
    reg100=reg80-reg100; reg72=reg129*reg154; reg80=reg120*reg119; reg145=reg146*reg152; T reg169=reg16*reg129;
    T reg170=reg28*reg151; reg76=reg76*reg45; T reg171=reg7*(*f.m).alpha_1; reg149=reg138+reg149; reg138=reg106*(*f.m).alpha_3;
    reg160=reg6+reg160; reg6=reg106*reg10; T reg172=reg51*reg131; reg161=reg150+reg161; reg25=reg25*reg24;
    reg97=reg143+reg97; reg143=reg144*reg74; reg150=reg16*reg66; reg133=reg157+reg133; reg157=reg66*reg154;
    reg162=reg116+reg162; reg116=reg9*reg74; T reg173=reg152*reg134; reg111=reg62+reg111; reg62=reg82*reg125;
    T reg174=reg144*reg123; reg164=reg132+reg164; reg132=reg144*reg75; reg101=reg99+reg101; reg118=reg90+reg118;
    reg109=reg109*reg42; reg90=reg38*reg102; reg156=reg112+reg156; reg108=reg108*reg42; reg15=reg15*reg2;
    reg99=reg4*(*f.m).alpha_2; reg112=reg59*reg102; reg92=reg50+reg92; reg30=reg155+reg30; reg50=reg136*reg123;
    reg166=reg167+reg166; reg168=reg17+reg168; reg165=reg158-reg165; reg17=reg136*reg74; reg85=reg130+reg85;
    reg33=reg137+reg33; reg153=reg107+reg153; reg5=reg3*reg5; reg26=reg27*reg26; reg27=reg2*reg51;
    reg107=reg25*(*f.m).alpha_2; reg130=reg0*reg40; reg6=reg149+reg6; reg137=reg4*reg0; reg149=reg26*(*f.m).alpha_1;
    reg155=reg51*reg42; reg158=reg5*(*f.m).alpha_3; reg167=reg59*reg76; T reg175=reg165*reg53; reg108=reg133+reg108;
    reg99=reg171+reg99; reg62=reg150-reg62; reg17=reg85+reg17; reg29=reg157-reg29; reg85=reg100*reg54;
    reg133=reg82*reg154; reg150=reg16*reg159; reg168=reg168*reg162; reg157=reg27*reg49; reg171=reg86*reg170;
    reg112=reg101+reg112; reg101=reg152*reg92; T reg176=reg129*reg125; reg174=reg164+reg174; reg164=reg15*reg51;
    reg50=reg166+reg50; reg166=reg44*reg170; reg90=reg118+reg90; reg30=reg109+reg30; reg109=reg27*reg51;
    reg118=reg38*reg76; T reg177=reg15*reg49; reg132=reg153+reg132; reg153=reg98*reg102; T reg178=reg27*reg55;
    reg138=reg160+reg138; reg160=reg165*reg56; reg172=reg161+reg172; reg143=reg97+reg143; reg28=reg33*reg28;
    reg97=reg15*reg55; reg161=reg100*reg57; reg80=reg163+reg80; reg52=reg3*reg52; reg72=reg173+reg72;
    reg111=reg116-reg111; reg156=reg162*reg156; reg135=reg89+reg135; reg169=reg145+reg169; reg89=0.78867513459481286553*elem.pos(0)[0];
    reg116=reg44*reg28; reg145=0.78867513459481286553*elem.pos(0)[1]; reg163=reg52*(*f.m).alpha_3; reg107=reg149+reg107; reg149=reg86*reg28;
    reg118=reg30+reg118; reg158=reg99+reg158; reg167=reg108+reg167; reg30=0.21132486540518713447*elem.pos(0)[0]; reg99=0.21132486540518713447*elem.pos(1)[0];
    reg108=0.21132486540518713447*elem.pos(0)[1]; reg173=0.78867513459481286553*elem.pos(1)[1]; reg150=reg133-reg150; reg171=reg112+reg171; reg112=reg146*reg29;
    reg133=reg62*reg134; T reg179=reg98*reg76; reg155=reg6+reg155; reg6=reg120*reg170; reg153=reg172+reg153;
    reg97=reg143+reg97; reg143=reg137*reg77; reg172=reg82*reg138; reg168=reg169+reg168; reg166=reg90+reg166;
    reg90=reg25*reg1; reg109=reg50+reg109; reg50=reg130*reg40; reg169=reg130*reg87; reg37=reg37*reg2;
    T reg180=reg137*reg87; T reg181=reg138*reg159; reg156=reg72+reg156; reg157=reg135+reg157; reg177=reg132+reg177;
    reg72=0.5*elem.pos(1)[0]; reg132=0.5*elem.pos(1)[1]; reg135=0.5*elem.pos(0)[1]; T reg182=0.78867513459481286553*elem.pos(1)[0]; reg164=reg174+reg164;
    reg174=reg111*reg74; T reg183=reg137*reg40; T reg184=reg1*reg34; T reg185=0.21132486540518713447*elem.pos(1)[1]; T reg186=reg130*reg77;
    reg178=reg17+reg178; reg17=reg100*reg94; reg161=reg160+reg161; reg160=reg111*reg75; T reg187=0.5*elem.pos(0)[0];
    reg85=reg175+reg85; reg80=reg162*reg80; reg176=reg101+reg176; reg101=reg165*reg93; reg175=reg187+reg72;
    T reg188=reg132-reg135; T reg189=reg92*reg159; reg183=reg164+reg183; reg164=reg90*reg34; T reg190=reg173+reg108;
    T reg191=reg184*reg96; reg186=reg178+reg186; reg178=0.78867513459481286553*elem.pos(2)[1]; T reg192=0.78867513459481286553*elem.pos(2)[0]; T reg193=reg99-reg30;
    T reg194=reg138*reg66; reg80=reg176+reg80; reg160=reg161+reg160; reg161=reg37*reg49; reg171=reg171*reg158;
    reg172=reg168+reg172; reg116=reg118+reg116; reg166=reg158*reg166; reg181=reg156+reg181; reg118=reg184*reg34;
    reg50=reg109+reg50; reg163=reg107+reg163; reg108=reg185-reg108; reg185=reg185+reg145; reg174=reg85+reg174;
    reg85=reg90*reg103; reg180=reg177+reg180; reg107=reg7*reg0; reg109=0.21132486540518713447*elem.pos(2)[1]; reg112=reg133-reg112;
    reg133=reg146*reg66; reg156=reg37*reg55; reg168=reg92*reg150; reg176=reg66*reg134; reg177=reg90*reg96;
    reg143=reg97+reg143; reg97=reg82*reg92; reg149=reg167+reg149; reg99=reg99+reg89; reg6=reg153+reg6;
    reg153=reg111*reg123; reg17=reg101+reg17; reg187=reg72-reg187; reg72=0.5*elem.pos(2)[0]; reg30=reg182+reg30;
    reg179=reg155+reg179; reg101=0.21132486540518713447*elem.pos(2)[0]; reg155=reg120*reg28; reg167=reg184*reg103; reg169=reg157+reg169;
    reg135=reg132+reg135; reg132=0.5*elem.pos(2)[1]; reg185=reg109-reg185; reg157=0.21132486540518713447*elem.pos(3)[1]; reg190=reg178-reg190;
    reg156=reg174+reg156; reg174=reg107*reg77; reg161=reg160+reg161; reg160=reg107*reg87; T reg195=reg26*reg1;
    reg145=reg173-reg145; reg85=reg180+reg85; reg89=reg182-reg89; reg118=reg50+reg118; reg164=reg183+reg164;
    reg167=reg169+reg167; reg50=reg37*reg51; reg153=reg17+reg153; reg17=reg82*reg134; reg169=reg92*reg154;
    reg189=reg176-reg189; reg188=reg132+reg188; reg175=reg72-reg175; reg173=0.5*elem.pos(3)[0]; reg187=reg72+reg187;
    reg72=0.5*elem.pos(3)[1]; reg155=reg179+reg155; reg135=reg132-reg135; reg132=reg125*reg134; reg168=reg112+reg168;
    reg112=reg16*reg92; reg193=reg192+reg193; reg97=reg133-reg97; reg133=reg146*reg125; reg99=reg101-reg99;
    reg176=0.78867513459481286553*elem.pos(3)[1]; reg178=reg108+reg178; reg108=0.21132486540518713447*elem.pos(3)[0]; reg30=reg192-reg30; reg179=0.78867513459481286553*elem.pos(3)[0];
    reg191=reg186+reg191; reg6=reg158*reg6; reg194=reg80+reg194; reg149=reg149*reg163; reg171=reg172+reg171;
    reg116=reg116*reg163; reg166=reg181+reg166; reg177=reg143+reg177; reg80=reg146*reg159; reg112=reg133-reg112;
    reg133=0.78867513459481286553*PNODE(1).dep[1]; reg143=0.21132486540518713447*PNODE(0).dep[1]; reg172=reg16*reg134; reg180=0.78867513459481286553*PNODE(1).dep[0]; reg181=0.21132486540518713447*PNODE(0).dep[0];
    reg62=reg62/reg168; reg182=0.5*vectors[0][indices[0]+0]; reg183=0.5*vectors[0][indices[1]+0]; reg186=reg195*reg96; reg174=reg156+reg174;
    reg99=reg179+reg99; reg156=reg177*reg118; reg149=reg171+reg149; reg171=reg85*reg118; reg97=reg97/reg168;
    reg169=reg132-reg169; reg178=reg178-reg176; reg6=reg194+reg6; reg155=reg163*reg155; reg185=reg176+reg185;
    reg30=reg30+reg108; reg190=reg190+reg157; reg132=reg191*reg164; reg179=reg193-reg179; reg188=reg188-reg72;
    reg176=reg146*reg154; reg175=reg173+reg175; reg192=reg107*reg40; reg50=reg153+reg50; reg101=reg89+reg101;
    reg173=reg187-reg173; reg72=reg135+reg72; reg109=reg145+reg109; reg89=reg167*reg164; reg189=reg189/reg168;
    reg135=reg195*reg103; reg160=reg161+reg160; reg145=0.21132486540518713447*PNODE(1).dep[0]; reg153=0.21132486540518713447*PNODE(1).dep[1]; reg116=reg166+reg116;
    reg29=reg29/reg168; reg80=reg17-reg80; reg17=0.5*vectors[0][indices[0]+1]; reg161=0.5*vectors[0][indices[1]+1]; reg166=0.78867513459481286553*PNODE(0).dep[0];
    reg187=0.78867513459481286553*PNODE(0).dep[1]; reg132=reg156-reg132; reg135=reg160+reg135; reg156=reg178*reg99; reg160=reg191*reg85;
    reg193=reg177*reg167; reg194=0.5*vectors[0][indices[2]+0]; T reg196=reg30*reg178; reg189=reg189*reg149; reg89=reg171-reg89;
    reg97=reg97*reg116; reg171=reg183-reg182; reg182=reg183+reg182; reg29=reg29*reg149; reg62=reg62*reg116;
    reg155=reg6+reg155; reg176=reg172-reg176; reg112=reg112/reg168; reg80=reg80/reg168; reg150=reg150/reg168;
    reg169=reg169/reg168; reg6=reg175*reg188; reg172=reg72*reg173; reg183=reg161-reg17; reg161=reg17+reg161;
    reg17=0.5*vectors[0][indices[2]+1]; T reg197=reg179*reg185; reg186=reg174+reg186; reg174=reg145-reg181; reg181=reg181+reg180;
    T reg198=0.78867513459481286553*PNODE(2).dep[1]; T reg199=reg143+reg133; T reg200=0.78867513459481286553*PNODE(2).dep[0]; reg143=reg153-reg143; T reg201=0.21132486540518713447*PNODE(2).dep[0];
    reg145=reg145+reg166; T reg202=0.21132486540518713447*PNODE(2).dep[1]; reg153=reg153+reg187; reg192=reg50+reg192; reg50=reg195*reg34;
    reg108=reg101-reg108; reg157=reg109-reg157; reg101=reg179*reg190; reg156=reg197-reg156; reg168=reg176/reg168;
    reg150=reg150*reg155; reg29=reg62-reg29; reg80=reg80*reg155; reg97=reg189-reg97; reg116=reg112*reg116;
    reg149=reg169*reg149; reg183=reg17+reg183; reg181=reg200-reg181; reg187=reg133-reg187; reg62=reg89*reg186;
    reg109=0.5*vectors[0][indices[3]+1]; reg112=0.21132486540518713447*PNODE(3).dep[0]; reg153=reg202-reg153; reg161=reg17-reg161; reg199=reg198-reg199;
    reg17=0.21132486540518713447*PNODE(3).dep[1]; reg196=reg101-reg196; reg143=reg198+reg143; reg101=0.78867513459481286553*PNODE(3).dep[1]; reg145=reg201-reg145;
    reg171=reg194+reg171; reg133=0.5*vectors[0][indices[3]+0]; reg182=reg194-reg182; reg6=reg172-reg6; reg166=reg180-reg166;
    reg50=reg192+reg50; reg169=reg108*reg190; reg172=reg30*reg157; reg160=reg193-reg160; reg176=0.78867513459481286553*PNODE(3).dep[0];
    reg200=reg174+reg200; reg174=reg135*reg132; reg149=reg116-reg149; reg72=reg72/reg6; reg155=reg168*reg155;
    reg166=reg201+reg166; reg116=reg177*reg50; reg80=reg97-reg80; reg29=reg150+reg29; reg97=reg191*reg50;
    reg145=reg176+reg145; reg150=reg186*reg164; reg172=reg169-reg172; reg153=reg101+reg153; reg161=reg161+reg109;
    reg202=reg187+reg202; reg168=1-(*f.m).resolution; reg109=reg183-reg109; reg169=reg167*reg50; reg180=reg135*reg164;
    reg183=reg157*reg99; reg187=reg186*reg118; reg189=reg85*reg50; reg192=reg108*reg185; reg193=reg135*reg118;
    reg194=reg50*reg160; reg174=reg62-reg174; reg62=reg190/reg196; reg197=reg179/reg156; reg198=reg99/reg156;
    reg201=reg185/reg156; T reg203=reg178/reg156; reg176=reg200-reg176; reg200=reg30/reg196; reg178=reg178/reg196;
    reg175=reg175/reg6; reg173=reg173/reg6; reg6=reg188/reg6; reg181=reg181+reg112; reg179=reg179/reg196;
    reg199=reg199+reg17; reg101=reg143-reg101; reg171=reg171-reg133; reg182=reg133+reg182; reg133=(*f.m).resolution*reg129;
    reg29=reg168*reg29; reg97=reg187-reg97; reg116=reg150-reg116; reg143=reg167*reg186; reg17=reg202-reg17;
    reg150=reg191*reg135; reg187=reg85*reg186; reg188=(*f.m).resolution*reg152; reg202=reg177*reg135; reg112=reg166-reg112;
    reg183=reg192-reg183; reg169=reg193-reg169; reg189=reg180-reg189; reg194=reg174+reg194; reg166=reg179*reg181;
    reg174=reg200*reg176; reg180=reg62*reg101; reg192=reg199*reg178; reg193=reg145*reg197; T reg204=reg176*reg198;
    T reg205=reg101*reg201; T reg206=reg153*reg203; T reg207=reg108/reg172; reg30=reg30/reg172; reg190=reg190/reg172;
    T reg208=reg182*reg173; T reg209=reg171*reg175; reg80=reg168*reg80; T reg210=reg72*reg109; T reg211=reg161*reg6;
    reg149=reg155+reg149; reg155=reg157/reg172; reg99=reg99/reg183; T reg212=(*f.m).resolution*reg138; reg80=reg133+reg80;
    reg133=reg62*reg176; T reg213=reg181*reg178; T reg214=reg179*reg199; T reg215=reg200*reg101; T reg216=reg181*reg207;
    T reg217=reg199*reg155; reg185=reg185/reg183; reg174=reg166-reg174; reg29=reg188+reg29; reg157=reg157/reg183;
    reg206=reg205-reg206; reg192=reg180-reg192; reg176=reg176*reg201; reg166=reg145*reg203; reg180=reg153*reg197;
    reg101=reg101*reg198; reg149=reg168*reg149; reg188=reg17*reg190; reg204=reg193-reg204; reg193=reg112*reg30;
    reg116=reg116/reg194; reg211=reg210-reg211; reg97=reg97/reg194; reg109=reg109*reg175; reg209=reg208-reg209;
    reg150=reg143-reg150; reg202=reg187-reg202; reg89=reg89/reg194; reg132=reg132/reg194; reg189=reg189/reg194;
    reg161=reg161*reg173; reg169=reg169/reg194; reg182=reg182*reg6; reg171=reg72*reg171; reg108=reg108/reg183;
    reg92=reg92*reg168; reg16=reg16*reg168; reg125=reg125*reg168; reg143=(*f.m).resolution*reg169; reg146=reg146*reg168;
    reg187=(*f.m).resolution*reg97; reg154=reg168*reg154; reg205=(*f.m).resolution*reg189; reg208=(*f.m).resolution*reg116; reg134=reg168*reg134;
    reg109=reg161-reg109; elem.epsilon[0][1]=reg109; reg211=reg209+reg211; reg149=reg212+reg149; reg29=reg29*(*f.m).deltaT;
    reg80=reg80*(*f.m).deltaT; reg206=reg204+reg206; reg101=reg180-reg101; reg166=reg176-reg166; reg192=reg174+reg192;
    reg215=reg214-reg215; reg213=reg133-reg213; reg129=reg129*(*f.m).deltaT; reg152=reg152*(*f.m).deltaT; reg133=reg157*reg153;
    reg161=reg185*reg17; reg174=reg99*reg112; reg176=reg108*reg145; reg217=reg188-reg217; reg193=reg216-reg193;
    reg160=reg160/reg194; reg180=reg112*reg190; reg181=reg181*reg155; reg182=reg171-reg182; elem.epsilon[0][0]=reg182;
    reg171=reg17*reg30; reg188=(*f.m).resolution*reg132; reg204=(*f.m).resolution*reg89; reg202=reg202/reg194; reg199=reg199*reg207;
    reg194=reg150/reg194; reg143=reg146-reg143; reg187=reg16+reg187; reg205=reg92+reg205; reg208=reg125-reg208;
    reg166=reg166-reg29; reg112=reg185*reg112; reg145=reg157*reg145; reg192=0.5*reg192; reg211=0.5*reg211;
    elem.epsilon[0][2]=reg211; reg153=reg108*reg153; reg17=reg99*reg17; reg149=reg149*(*f.m).deltaT; reg213=reg213-reg29;
    reg215=reg215-reg80; reg206=0.5*reg206; reg16=reg109-reg129; reg171=reg199-reg171; reg138=reg138*(*f.m).deltaT;
    reg181=reg180-reg181; reg92=(*f.m).resolution*reg160; reg101=reg101-reg80; reg125=reg182-reg152; reg174=reg176-reg174;
    reg82=reg82*reg168; reg133=reg161-reg133; reg146=(*f.m).resolution*reg194; reg66=reg66*reg168; reg217=reg193+reg217;
    reg159=reg168*reg159; reg150=(*f.m).resolution*reg202; reg134=reg204+reg134; reg188=reg154-reg188; reg154=reg16*reg132;
    reg161=reg16*reg97; reg168=reg215*reg188; reg176=reg125*reg169; reg180=reg213*reg134; reg133=reg174+reg133;
    reg171=reg171-reg80; reg174=reg166*reg134; reg193=reg101*reg208; reg199=reg166*reg205; reg204=reg125*reg89;
    reg209=reg101*reg187; reg210=reg166*reg143; reg206=reg206-reg149; reg212=reg101*reg188; reg181=reg181-reg29;
    reg214=reg211-reg138; reg216=reg213*reg143; T reg218=reg215*reg187; reg17=reg153-reg17; reg145=reg112-reg145;
    reg217=0.5*reg217; reg112=reg213*reg205; reg153=reg215*reg208; reg150=reg66+reg150; reg192=reg192-reg149;
    reg146=reg82-reg146; reg159=reg92+reg159; reg66=reg206*reg159; reg82=reg181*reg205; reg212=reg174+reg212;
    reg92=reg192*reg159; reg217=reg217-reg149; reg209=reg210+reg209; reg174=reg206*reg146; reg210=reg214*reg160;
    T reg219=reg181*reg134; reg145=reg145-reg29; T reg220=reg171*reg208; T reg221=reg171*reg188; T reg222=reg192*reg150;
    reg133=0.5*reg133; reg168=reg180+reg168; reg153=reg112+reg153; reg154=reg204-reg154; reg112=reg171*reg187;
    reg176=reg161-reg176; reg161=reg214*reg194; reg17=reg17-reg80; reg180=reg181*reg143; reg204=reg192*reg146;
    reg218=reg216+reg218; reg216=reg206*reg150; reg193=reg199+reg193; reg204=reg218+reg204; reg199=reg217*reg146;
    reg112=reg180+reg112; reg92=reg168+reg92; reg222=reg153+reg222; reg153=reg17*reg188; reg221=reg219+reg221;
    reg168=reg17*reg208; reg180=reg145*reg205; reg210=reg154+reg210; reg161=reg176-reg161; reg154=reg17*reg187;
    reg133=reg133-reg149; reg176=reg217*reg159; reg216=reg193+reg216; reg193=reg145*reg134; reg174=reg209+reg174;
    reg66=reg212+reg66; reg209=reg145*reg143; reg212=reg217*reg150; reg220=reg82+reg220; reg92=reg213*reg92;
    reg82=reg133*reg159; reg154=reg209+reg154; reg204=reg215*reg204; reg209=reg133*reg146; reg199=reg112+reg199;
    reg212=reg220+reg212; reg66=reg166*reg66; reg176=reg221+reg176; reg174=reg101*reg174; reg216=2*reg216;
    reg101=reg210+reg161; reg112=reg133*reg150; reg168=reg180+reg168; reg153=reg193+reg153; reg222=2*reg222;
    reg16=reg16*reg116; reg125=reg125*reg189; reg101=reg101/3; reg112=reg168+reg112; reg209=reg154+reg209;
    reg82=reg153+reg82; reg199=reg171*reg199; reg176=reg181*reg176; reg212=2*reg212; reg216=reg206*reg216;
    reg174=reg66+reg174; reg222=reg192*reg222; reg204=reg92+reg204; reg216=reg174+reg216; reg209=reg17*reg209;
    reg82=reg145*reg82; reg109=reg109-reg80; reg112=2*reg112; reg212=reg217*reg212; reg182=reg182-reg29;
    reg161=reg161-reg101; reg16=reg125-reg16; reg199=reg176+reg199; reg210=reg210-reg101; reg214=reg214*reg202;
    reg222=reg204+reg222; reg216=reg216*reg156; reg161=pow(reg161,2); reg17=reg182*reg134; reg209=reg82+reg209;
    reg66=reg182*reg143; reg82=reg109*reg187; reg211=reg211-reg149; reg222=reg196*reg222; reg92=reg109*reg188;
    reg112=reg133*reg112; reg212=reg199+reg212; reg210=pow(reg210,2); reg16=reg214+reg16; reg125=reg211*reg146;
    reg82=reg66+reg82; reg182=reg182*reg205; reg109=reg109*reg208; reg112=reg209+reg112; reg216=0.25*reg216;
    reg161=reg210+reg161; reg66=reg211*reg159; reg92=reg17+reg92; reg212=reg212*reg172; reg17=2*reg16;
    reg222=0.25*reg222; reg101=pow(reg101,2); reg66=reg92+reg66; elem.sigma[0][0]=reg66; reg125=reg82+reg125;
    elem.sigma[0][1]=reg125; reg212=0.25*reg212; reg101=reg161+reg101; reg17=reg16*reg17; reg216=reg222+reg216;
    reg211=reg211*reg150; reg112=reg183*reg112; reg109=reg182+reg109; reg112=0.25*reg112; reg17=reg101+reg17;
    reg16=reg94*reg125; reg82=reg93*reg66; reg216=reg212+reg216; reg92=reg125*reg57; reg211=reg109+reg211;
    elem.sigma[0][2]=reg211; reg101=reg66*reg53; reg125=reg125*reg54; reg66=reg66*reg56; reg17=1.5*reg17;
    reg109=reg211*reg55; reg125=reg101+reg125; reg112=reg216+reg112; reg16=reg82+reg16; reg82=reg51*reg211;
    reg211=reg49*reg211; reg92=reg66+reg92; elem.sigma_local[0][1]=reg92+reg211; elem.sigma_local[0][0]=reg125+reg109; elem.sigma_von_mises=pow(reg17,0.5);
    elem.ener=reg112/2; elem.sigma_local[0][2]=reg16+reg82;
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_2(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_3(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_4(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_5(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_6(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_7(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_8(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_9(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_10(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_11(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_12(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_13(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_14(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_15(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
      #undef PNODE
    }
  
  static const unsigned nb_elementary_unknowns = 0;
  template<class TE,class TTs,class Tvec>
  inline static void set_elementary_unknowns(TE &elem,const TTs &f,const Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_elementary_unknowns(TE &elem,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_elementary_initial_conditions(const TE &elem,const TTs &f,Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_elementary_error(const TE &elem,const TTs &f,const Tvec &vecs,int indice) {
    return 0;
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_elementary(const TE &elem,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
  }
  
  static const unsigned nb_skin_elementary_unknowns = 0;
  template<class TE,class TTs,class Tvec>
  inline static void set_skin_elementary_unknowns(TE &elem,const TTs &f,const Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_skin_elementary_unknowns(TE &elem,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_skin_elementary_initial_conditions(const TE &elem,const TTs &f,Tvec &vecs,unsigned indice) {
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_skin_elementary_error(const TE &elem,const TTs &f,const Tvec &vecs,int indice) {
    return 0;
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_skin_elementary(const TE &elem,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
  }
};

// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=2*(*f.m).shear_modulus_13; T reg1=2*(*f.m).shear_modulus_23; T reg2=2*(*f.m).shear_modulus_12; reg0=1.0/reg0; reg1=1.0/reg1;
    reg2=1.0/reg2; T reg3=reg0*reg1; T reg4=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg5=1.0/(*f.m).elastic_modulus_3; T reg6=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg7=reg2*reg3; T reg8=pow((*f.m).v1[0],2); T reg9=pow((*f.m).v1[1],2); T reg10=pow((*f.m).v2[0],2); T reg11=pow((*f.m).v2[1],2);
    reg9=reg8+reg9; reg8=pow((*f.m).v1[2],2); reg11=reg10+reg11; reg10=pow((*f.m).v2[2],2); T reg12=1.0/(*f.m).elastic_modulus_2;
    T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg14=reg5*reg7; T reg15=reg4*reg7; T reg16=reg6*reg7; reg10=reg11+reg10;
    reg11=reg12*reg14; reg8=reg9+reg8; reg9=reg4*reg15; T reg17=reg4*reg16; T reg18=reg13*reg14;
    T reg19=1.0/(*f.m).elastic_modulus_1; reg9=reg11-reg9; reg10=pow(reg10,0.5); reg17=reg18+reg17; reg8=pow(reg8,0.5);
    reg11=reg12*reg16; T reg20=reg13*reg15; T reg21=(*f.m).v1[2]/reg8; T reg22=(*f.m).v1[1]/reg8; T reg23=reg20+reg11;
    T reg24=reg19*reg9; T reg25=(*f.m).v2[1]/reg10; T reg26=(*f.m).v2[2]/reg10; T reg27=reg13*reg17; T reg28=reg6*reg23;
    T reg29=reg12*reg7; reg27=reg24-reg27; reg24=reg6*reg15; reg14=reg19*reg14; T reg30=reg6*reg16;
    reg7=reg13*reg7; T reg31=reg4*reg3; T reg32=reg6*reg3; T reg33=reg2*reg1; T reg34=reg21*reg25;
    reg8=(*f.m).v1[0]/reg8; reg10=(*f.m).v2[0]/reg10; reg3=reg5*reg3; T reg35=reg22*reg26; reg16=reg13*reg16;
    T reg36=reg12*reg3; T reg37=2*reg10; T reg38=reg2*reg0; T reg39=reg6*reg7; reg3=reg13*reg3;
    reg15=reg19*reg15; reg30=reg14-reg30; reg14=reg4*reg31; T reg40=reg4*reg32; T reg41=reg6*reg29;
    reg24=reg18+reg24; reg18=reg5*reg33; T reg42=reg6*reg33; T reg43=2*reg8; T reg44=reg35-reg34;
    reg33=reg4*reg33; reg28=reg27-reg28; reg27=reg8*reg26; T reg45=reg21*reg10; T reg46=pow(reg25,2);
    T reg47=pow(reg10,2); reg32=reg12*reg32; reg40=reg3+reg40; reg14=reg36-reg14; reg3=reg22*reg43;
    reg36=pow(reg22,2); T reg48=pow(reg8,2); T reg49=reg8*reg25; T reg50=reg45-reg27; T reg51=reg6*reg38;
    reg29=reg19*reg29; reg9=reg9/reg28; T reg52=reg22*reg10; T reg53=reg5*reg38; T reg54=2*reg44;
    reg31=reg13*reg31; T reg55=reg4*reg42; reg7=reg13*reg7; reg16=reg15+reg16; reg38=reg4*reg38;
    T reg56=reg37*reg25; reg39=reg15+reg39; reg15=reg4*reg33; reg24=reg24/reg28; reg41=reg20+reg41;
    reg17=reg17/reg28; T reg57=reg12*reg18; reg18=reg13*reg18; reg30=reg30/reg28; reg39=reg39/reg28;
    T reg58=reg12*reg53; reg53=reg13*reg53; T reg59=reg4*reg38; T reg60=reg4*reg51; reg33=reg13*reg33;
    T reg61=reg9*reg36; T reg62=reg17*reg46; T reg63=reg9*reg3; T reg64=reg17*reg56; reg55=reg18+reg55;
    reg42=reg12*reg42; reg18=reg31+reg32; reg41=reg41/reg28; reg7=reg29-reg7; reg15=reg57-reg15;
    reg29=reg9*reg48; reg40=reg13*reg40; reg16=reg16/reg28; reg23=reg23/reg28; reg57=reg17*reg47;
    reg14=reg19*reg14; T reg65=reg24*reg3; T reg66=pow(reg50,2); T reg67=reg30*reg56; T reg68=reg30*reg47;
    T reg69=reg24*reg48; T reg70=pow(reg44,2); T reg71=reg24*reg36; T reg72=reg54*reg50; T reg73=reg49-reg52;
    T reg74=reg30*reg46; T reg75=pow(reg26,2); T reg76=pow(reg21,2); T reg77=reg24*reg76; T reg78=reg30*reg75;
    reg60=reg53+reg60; reg38=reg13*reg38; reg53=reg16*reg66; T reg79=reg41*reg36; T reg80=pow(reg73,2);
    reg67=reg65+reg67; reg65=reg16*reg72; reg7=reg7/reg28; reg55=reg13*reg55; T reg81=reg33+reg42;
    T reg82=reg39*reg47; T reg83=reg41*reg48; reg59=reg58-reg59; reg58=reg9*reg76; T reg84=reg17*reg75;
    T reg85=reg23*reg66; reg62=reg61+reg62; reg15=reg19*reg15; reg61=reg41*reg3; reg40=reg14-reg40;
    reg18=reg6*reg18; reg14=reg39*reg56; reg64=reg63+reg64; reg57=reg29+reg57; reg29=reg23*reg70;
    reg63=reg23*reg72; reg68=reg69+reg68; reg69=reg16*reg70; reg51=reg12*reg51; T reg86=reg39*reg46;
    reg74=reg71+reg74; reg60=reg13*reg60; reg63=reg64+reg63; reg18=reg40-reg18; reg40=reg23*reg80;
    reg64=2*reg25; reg71=reg37*reg26; T reg87=reg38+reg51; reg84=reg58+reg84; reg85=reg62+reg85;
    reg58=reg7*reg66; reg62=reg41*reg76; T reg88=reg39*reg75; reg86=reg79+reg86; reg79=reg10*reg25;
    T reg89=reg8*reg22; reg14=reg61+reg14; reg29=reg57+reg29; reg72=reg7*reg72; reg69=reg68+reg69;
    reg53=reg74+reg53; reg78=reg77+reg78; reg57=reg16*reg80; reg61=reg21*reg43; reg55=reg15-reg55;
    reg15=2*reg22; reg65=reg67+reg65; reg81=reg6*reg81; reg59=reg19*reg59; reg82=reg83+reg82;
    reg67=reg7*reg70; reg58=reg86+reg58; reg88=reg62+reg88; reg62=reg7*reg80; reg60=reg59-reg60;
    reg59=reg85*reg36; reg68=reg8*reg10; reg74=reg53*reg46; reg49=reg52+reg49; reg81=reg55-reg81;
    reg52=reg22*reg25; reg72=reg14+reg72; reg87=reg6*reg87; reg14=reg22*reg44; reg55=reg8*reg50;
    reg67=reg82+reg67; reg77=reg79*reg65; reg82=reg53*reg47; reg83=reg85*reg48; reg86=reg89*reg63;
    T reg90=reg30*reg71; T reg91=reg24*reg61; T reg92=reg69*reg47; T reg93=reg29*reg48; reg57=reg78+reg57;
    reg78=reg10*reg43; T reg94=reg25*reg15; T reg95=2*reg21; T reg96=reg63*reg36; T reg97=reg65*reg46;
    T reg98=reg50*reg44; reg85=reg89*reg85; reg53=reg79*reg53; T reg99=2*reg50; reg54=reg54*reg73;
    T reg100=reg69*reg46; T reg101=reg9*reg61; T reg102=reg17*reg71; reg40=reg84+reg40; reg65=reg65*reg47;
    reg63=reg63*reg48; reg84=reg64*reg26; T reg103=reg29*reg36; reg18=reg18/reg28; T reg104=reg21*reg15;
    reg69=reg79*reg69; reg29=reg89*reg29; T reg105=reg66*reg72; reg97=reg96+reg97; reg96=reg66*reg67;
    reg53=reg85+reg53; reg85=reg98*reg58; reg100=reg103+reg100; reg74=reg59+reg74; reg59=reg66*reg58;
    reg103=reg57*reg46; T reg106=reg40*reg36; reg55=reg14+reg55; reg14=reg25*reg44; T reg107=reg10*reg50;
    reg77=reg86+reg77; reg86=reg98*reg72; T reg108=reg78*reg19; T reg109=reg13*reg46; T reg110=reg19*reg47;
    reg102=reg101+reg102; reg72=reg70*reg72; reg65=reg63+reg65; reg63=reg57*reg47; reg101=reg18*reg52;
    T reg111=reg40*reg48; reg30=reg30*reg84; reg24=reg24*reg104; T reg112=reg16*reg54; reg58=reg70*reg58;
    reg82=reg83+reg82; reg90=reg91+reg90; reg83=reg70*reg67; reg92=reg93+reg92; reg91=reg18*reg68;
    reg93=reg23*reg54; reg95=reg95*reg26; reg62=reg88+reg62; reg88=reg41*reg61; reg81=reg81/reg28;
    T reg113=reg39*reg71; T reg114=reg18*reg49; T reg115=reg22*reg50; T reg116=reg8*reg44; T reg117=reg21*reg26;
    reg9=reg9*reg104; T reg118=reg94*reg13; reg99=reg99*reg73; reg17=reg17*reg84; T reg119=reg13*reg47;
    reg87=reg60-reg87; reg60=reg12*reg46; T reg120=reg78*reg13; T reg121=reg94*reg12; T reg122=reg94*reg114;
    reg105=reg97+reg105; reg69=reg29+reg69; reg85=reg53+reg85; reg29=reg49*reg101; reg17=reg9+reg17;
    reg23=reg23*reg99; reg93=reg102+reg93; reg67=reg98*reg67; reg9=reg66*reg62; reg103=reg106+reg103;
    reg53=reg94*reg101; reg112=reg90+reg112; reg30=reg24+reg30; reg16=reg16*reg99; reg59=reg74+reg59;
    reg96=reg100+reg96; reg24=reg94*reg91; reg107=reg14+reg107; reg28=reg87/reg28; reg14=reg25*reg50;
    reg74=reg10*reg44; reg87=reg21*reg73; reg120=reg121-reg120; reg90=reg81*reg55; reg97=reg94*reg4;
    reg100=reg78*reg6; reg109=reg110-reg109; reg43=reg44*reg43; reg102=reg6*reg75; reg15=reg50*reg15;
    reg106=reg6*reg47; reg101=reg78*reg101; reg54=reg7*reg54; reg110=reg49*reg114; reg115=reg81*reg115;
    reg116=reg81*reg116; reg86=reg77+reg86; reg113=reg88+reg113; reg41=reg41*reg104; reg77=reg4*reg46;
    reg72=reg65+reg72; reg114=reg78*reg114; reg65=reg78*reg91; reg88=reg4*reg75; reg83=reg92+reg83;
    reg39=reg39*reg84; reg119=reg60-reg119; reg60=reg70*reg62; reg63=reg111+reg63; reg92=reg18*reg117;
    reg118=reg108-reg118; reg108=reg95*reg4; reg111=reg95*reg6; reg58=reg82+reg58; reg40=reg89*reg40;
    reg57=reg79*reg57; reg82=reg94*reg92; reg9=reg103+reg9; reg108=reg120-reg108; reg88=reg119-reg88;
    reg62=reg98*reg62; reg122=reg105+reg122; reg103=reg15*reg90; reg105=reg78*reg92; reg60=reg63+reg60;
    reg114=reg72+reg114; reg63=reg43*reg90; reg72=(*f.m).alpha_1*reg48; reg119=reg93*reg48; reg120=reg112*reg47;
    reg45=reg27+reg45; reg97=reg100+reg97; reg95=reg95*reg5; reg77=reg106+reg77; reg27=reg5*reg75;
    reg90=reg55*reg90; reg110=reg86+reg110; reg53=reg59+reg53; reg59=reg15*reg115; reg16=reg30+reg16;
    reg30=reg21*reg44; reg86=reg8*reg73; reg100=(*f.m).alpha_1*reg36; reg106=(*f.m).alpha_2*reg46; reg121=reg43*reg115;
    reg23=reg17+reg23; reg101=reg58+reg101; reg64=reg64*reg50; reg37=reg37*reg44; reg87=reg81*reg87;
    reg17=reg15*reg116; reg24=reg96+reg24; reg57=reg40+reg57; reg67=reg69+reg67; reg91=reg49*reg91;
    reg74=reg28*reg74; reg14=reg28*reg14; reg40=reg28*reg107; reg19=reg19*reg48; reg115=reg55*reg115;
    reg65=reg83+reg65; reg58=reg43*reg116; reg29=reg85+reg29; reg12=reg12*reg36; reg69=reg13*reg48;
    reg83=reg93*reg36; reg85=reg112*reg46; reg111=reg118-reg111; reg96=reg26*reg73; reg102=reg109-reg102;
    reg54=reg113+reg54; reg13=reg13*reg36; reg39=reg41+reg39; reg99=reg7*reg99; reg7=(*f.m).alpha_2*reg47;
    reg58=reg65+reg58; reg41=reg88*reg36; reg65=reg102*reg48; reg109=reg4*reg76; reg69=reg12-reg69;
    reg12=reg37*reg14; reg77=reg27-reg77; reg27=reg37*reg74; reg113=reg18*reg45; reg96=reg28*reg96;
    reg99=reg39+reg99; reg4=reg4*reg36; reg121=reg101+reg121; reg39=reg6*reg48; reg97=reg95-reg97;
    reg95=reg64*reg40; reg101=reg16*reg47; reg85=reg83+reg85; reg83=reg66*reg54; reg118=reg23*reg36;
    T reg123=reg16*reg46; T reg124=reg26*reg44; T reg125=reg10*reg73; reg6=reg6*reg76; reg13=reg19-reg13;
    reg7=reg72+reg7; reg19=reg70*(*f.m).alpha_3; reg106=reg100+reg106; reg72=reg66*(*f.m).alpha_3; reg100=(*f.m).alpha_1*reg76;
    T reg126=(*f.m).alpha_2*reg75; T reg127=reg108*reg46; T reg128=reg111*reg47; T reg129=reg64*reg74; reg17=reg24+reg17;
    reg24=reg88*reg46; reg91=reg67+reg91; reg116=reg55*reg116; reg67=reg102*reg47; T reg130=reg107*reg14;
    reg115=reg29+reg115; reg29=reg111*reg68; T reg131=reg108*reg52; reg108=reg108*reg36; reg111=reg111*reg48;
    reg105=reg60+reg105; reg60=reg43*reg87; reg63=reg114+reg63; reg114=reg37*reg40; reg35=reg34+reg35;
    reg120=reg119+reg120; reg112=reg79*reg112; reg93=reg89*reg93; reg34=reg70*reg54; reg40=reg107*reg40;
    reg90=reg110+reg90; reg59=reg53+reg59; reg14=reg64*reg14; reg86=reg30+reg86; reg82=reg9+reg82;
    reg103=reg122+reg103; reg9=reg23*reg48; reg88=reg88*reg52; reg102=reg102*reg68; reg62=reg57+reg62;
    reg92=reg49*reg92; reg30=reg22*reg73; reg53=reg21*reg50; reg57=reg15*reg87; reg110=reg25*reg73;
    reg60=reg105+reg60; reg105=reg37*reg96; reg119=reg26*reg50; reg122=reg78*reg113; reg108=reg111+reg108;
    reg111=reg97*reg76; reg114=reg63+reg114; reg12=reg121+reg12; reg101=reg9+reg101; reg27=reg58+reg27;
    reg70=reg70*reg99; reg24=reg67+reg24; reg9=reg97*reg75; reg34=reg120+reg34; reg129=reg17+reg129;
    reg17=reg77*reg75; reg30=reg53+reg30; reg127=reg128+reg127; reg125=reg124+reg125; reg53=reg79*(*f.m).alpha_2;
    reg58=reg89*(*f.m).alpha_1; reg16=reg79*reg16; reg23=reg89*reg23; reg54=reg98*reg54; reg112=reg93+reg112;
    reg40=reg90+reg40; reg14=reg59+reg14; reg80=reg80*(*f.m).alpha_3; reg126=reg100+reg126; reg72=reg106+reg72;
    reg19=reg7+reg19; reg87=reg55*reg87; reg92=reg62+reg92; reg57=reg82+reg57; reg7=reg64*reg96;
    reg59=reg77*reg117; reg88=reg102+reg88; reg95=reg103+reg95; reg83=reg85+reg83; reg62=reg94*reg113;
    reg123=reg118+reg123; reg66=reg66*reg99; reg97=reg97*reg117; reg86=reg81*reg86; reg8=reg8*reg21;
    reg10=reg10*reg26; reg18=reg18*reg35; reg6=reg13-reg6; reg109=reg69-reg109; reg5=reg5*reg76;
    reg4=reg39+reg4; reg79=reg79*reg2; reg13=reg2*reg49; reg41=reg65+reg41; reg77=reg77*reg76;
    reg131=reg29+reg131; reg74=reg107*reg74; reg116=reg91+reg116; reg130=reg115+reg130; reg74=reg116+reg74;
    reg70=reg101+reg70; reg80=reg126+reg80; reg29=reg10*(*f.m).alpha_2; reg39=reg8*(*f.m).alpha_1; reg63=reg98*(*f.m).alpha_3;
    reg94=reg94*reg18; reg66=reg123+reg66; reg78=reg78*reg18; reg53=reg58+reg53; reg125=reg28*reg125;
    reg30=reg81*reg30; reg58=reg14*reg72; reg65=reg129*reg19; reg67=reg72*reg12; reg69=reg19*reg27;
    reg105=reg60+reg105; reg60=reg130*reg114; reg81=reg95*reg130; reg82=reg40*reg12; reg85=reg14*reg40;
    reg99=reg98*reg99; reg16=reg23+reg16; reg113=reg49*reg113; reg54=reg112+reg54; reg96=reg107*reg96;
    reg122=reg34+reg122; reg23=reg43*reg86; reg87=reg92+reg87; reg7=reg57+reg7; reg62=reg83+reg62;
    reg34=reg15*reg86; reg57=reg6*reg47; reg83=reg109*reg46; reg90=reg79*reg3; reg77=reg41+reg77;
    reg17=reg24+reg17; reg24=reg79*reg56; reg41=reg109*reg36; reg91=reg6*reg48; reg92=reg0*reg45;
    reg9=reg127+reg9; reg93=reg13*reg56; reg10=reg10*reg0; reg4=reg5-reg4; reg26=reg25*reg26;
    reg21=reg22*reg21; reg110=reg119+reg110; reg44=reg73*reg44; reg59=reg88+reg59; reg79=reg79*reg49;
    reg5=reg13*reg3; reg13=reg13*reg49; reg97=reg131+reg97; reg111=reg108+reg111; reg13=reg97+reg13;
    reg22=reg1*reg35; reg25=reg26*reg1; reg6=reg6*reg68; reg109=reg109*reg52; reg2=reg89*reg2;
    reg96=reg87+reg96; reg87=reg14*reg114; reg88=reg95*reg12; reg7=reg7*reg80; reg50=reg73*reg50;
    reg73=reg10*reg45; reg79=reg59+reg79; reg18=reg49*reg18; reg59=reg64*reg125; reg89=reg19*reg74;
    reg97=reg92*reg45; reg99=reg16+reg99; reg43=reg43*reg30; reg70=reg78+reg70; reg16=reg10*reg61;
    reg90=reg77+reg90; reg34=reg62+reg34; reg94=reg66+reg94; reg15=reg15*reg30; reg76=reg4*reg76;
    reg62=reg72*reg130; reg41=reg91+reg41; reg86=reg55*reg86; reg113=reg54+reg113; reg44=reg44*(*f.m).alpha_3;
    reg29=reg39+reg29; reg63=reg53+reg63; reg39=reg21*(*f.m).alpha_1; reg81=reg85-reg81; reg67=reg69+reg67;
    reg5=reg111+reg5; reg53=reg92*reg61; reg105=reg80*reg105; reg83=reg57+reg83; reg75=reg4*reg75;
    reg26=reg26*(*f.m).alpha_2; reg60=reg82-reg60; reg24=reg17+reg24; reg10=reg10*reg71; reg17=reg37*reg125;
    reg28=reg110*reg28; reg58=reg65+reg58; reg93=reg9+reg93; reg92=reg92*reg71; reg23=reg122+reg23;
    reg9=reg22*reg104; reg59=reg34+reg59; reg53=reg5+reg53; reg30=reg55*reg30; reg18=reg99+reg18;
    reg5=reg95*reg63; reg34=reg129*reg60; reg125=reg107*reg125; reg62=reg89+reg62; reg86=reg113+reg86;
    reg7=reg58+reg7; reg54=reg81*reg27; reg96=reg80*reg96; reg92=reg93+reg92; reg117=reg4*reg117;
    reg4=reg22*reg84; reg44=reg29+reg44; reg109=reg6+reg109; reg6=reg25*reg84; reg10=reg24+reg10;
    reg50=reg50*(*f.m).alpha_3; reg16=reg90+reg16; reg24=reg25*reg104; reg76=reg41+reg76; reg29=reg2*reg3;
    reg26=reg39+reg26; reg64=reg64*reg28; reg15=reg94+reg15; reg0=reg8*reg0; reg73=reg79+reg73;
    reg8=reg63*reg114; reg105=reg67+reg105; reg39=reg2*reg56; reg75=reg83+reg75; reg17=reg23+reg17;
    reg43=reg70+reg43; reg37=reg37*reg28; reg25=reg25*reg35; reg87=reg88-reg87; reg22=reg22*reg35;
    reg97=reg13+reg97; reg50=reg26+reg50; reg25=reg73+reg25; reg13=reg74*reg114; reg23=reg63*reg40;
    reg96=reg62+reg96; reg22=reg97+reg22; reg8=reg105+reg8; reg17=reg44*reg17; reg64=reg15+reg64;
    reg59=reg59*reg44; reg5=reg7+reg5; reg4=reg92+reg4; reg6=reg10+reg6; reg34=reg54-reg34;
    reg1=reg21*reg1; reg29=reg76+reg29; reg61=reg0*reg61; reg7=reg74*reg87; reg37=reg43+reg37;
    reg10=reg129*reg40; reg71=reg0*reg71; reg28=reg107*reg28; reg30=reg18+reg30; reg39=reg75+reg39;
    reg125=reg86+reg125; reg9=reg53+reg9; reg24=reg16+reg24; reg15=reg95*reg74; reg16=reg40*reg27;
    reg2=reg2*reg49; reg117=reg109+reg117; reg7=reg34+reg7; reg18=reg74*reg12; reg13=reg16-reg13;
    reg16=1-var_inter[0]; reg37=reg37*reg50; reg17=reg8+reg17; reg84=reg1*reg84; reg71=reg39+reg71;
    reg8=reg130*reg27; reg21=reg129*reg114; reg45=reg0*reg45; reg2=reg117+reg2; reg0=reg14*reg74;
    reg15=reg10-reg15; reg125=reg44*reg125; reg23=reg96+reg23; reg10=reg129*reg130; reg26=reg95*reg27;
    reg34=reg9*reg25; reg28=reg30+reg28; reg30=1-var_inter[1]; reg64=reg64*reg50; reg104=reg1*reg104;
    reg61=reg29+reg61; reg59=reg5+reg59; reg5=reg4*reg25; reg29=reg6*reg22; reg39=reg24*reg22;
    reg18=reg8-reg18; reg8=reg9*reg6; reg13=reg13/reg7; reg41=reg16*elem.pos(0)[0]; reg60=reg60/reg7;
    reg0=reg10-reg0; reg81=reg81/reg7; reg15=reg15/reg7; reg10=var_inter[0]*elem.pos(1)[1]; reg43=reg16*elem.pos(0)[1];
    reg45=reg2+reg45; reg35=reg1*reg35; reg1=reg129*reg12; reg2=reg30*elem.pos(0)[0]; reg21=reg26-reg21;
    reg26=reg30*elem.pos(1)[0]; reg44=reg14*reg27; reg37=reg17+reg37; reg34=reg39-reg34; reg17=var_inter[0]*elem.pos(1)[0];
    reg84=reg71+reg84; reg5=reg29-reg5; reg29=reg30*elem.pos(1)[1]; reg39=reg30*elem.pos(0)[1]; reg64=reg59+reg64;
    reg53=reg24*reg4; reg28=reg50*reg28; reg104=reg61+reg104; reg125=reg23+reg125; reg28=reg125+reg28;
    reg23=var_inter[1]*elem.pos(2)[0]; reg26=reg26-reg2; reg50=var_inter[0]*elem.pos(2)[0]; reg54=reg5*reg104; reg35=reg45+reg35;
    reg29=reg29-reg39; reg45=var_inter[1]*elem.pos(2)[1]; reg13=reg13*reg64; reg8=reg53-reg8; reg15=reg15*reg37;
    reg21=reg21/reg7; reg53=reg41+reg17; reg55=reg84*reg34; reg60=reg60*reg64; reg87=reg87/reg7;
    reg18=reg18/reg7; reg81=reg81*reg37; reg1=reg44-reg1; reg0=reg0/reg7; reg44=reg10+reg43;
    reg57=var_inter[0]*elem.pos(2)[1]; reg37=reg0*reg37; reg64=reg18*reg64; reg15=reg13-reg15; reg0=reg16*elem.pos(3)[1];
    reg57=reg57-reg44; reg87=reg87*reg28; reg60=reg81-reg60; reg21=reg21*reg28; reg50=reg50-reg53;
    reg13=reg16*elem.pos(3)[0]; reg55=reg54-reg55; reg18=reg35*reg8; reg54=var_inter[1]*elem.pos(3)[1]; reg45=reg29+reg45;
    reg29=var_inter[1]*elem.pos(3)[0]; reg26=reg23+reg26; reg23=reg104*reg25; reg25=reg84*reg25; reg7=reg1/reg7;
    reg1=reg24*reg35; reg58=reg6*reg35; reg24=reg24*reg84; reg6=reg6*reg104; reg1=reg23-reg1;
    reg23=reg9*reg35; reg21=reg15-reg21; reg64=reg37-reg64; reg28=reg7*reg28; reg7=reg84*reg22;
    reg35=reg4*reg35; reg22=reg104*reg22; reg58=reg25-reg58; reg18=reg55+reg18; reg15=1-(*f.m).resolution;
    reg26=reg26-reg29; reg45=reg45-reg54; reg13=reg50+reg13; reg60=reg87+reg60; reg0=reg57+reg0;
    reg24=reg6-reg24; reg35=reg7-reg35; reg84=reg9*reg84; reg104=reg4*reg104; reg58=reg58/reg18;
    reg4=reg0*reg26; reg1=reg1/reg18; reg23=reg22-reg23; reg6=(*f.m).resolution*reg19; reg7=reg13*reg45;
    reg21=reg15*reg21; reg60=reg15*reg60; reg9=(*f.m).resolution*reg72; reg64=reg28+reg64; reg74=reg74*reg15;
    reg130=reg130*reg15; reg64=reg15*reg64; reg22=(*f.m).resolution*reg63; reg25=(*f.m).resolution*reg58; reg21=reg9+reg21;
    reg60=reg6+reg60; reg6=(*f.m).resolution*reg1; reg7=reg4-reg7; reg34=reg34/reg18; reg35=reg35/reg18;
    reg5=reg5/reg18; reg23=reg23/reg18; reg84=reg104-reg84; reg24=reg24/reg18; reg129=reg129*reg15;
    reg14=reg14*reg15; reg12=reg15*reg12; reg27=reg15*reg27; reg40=reg40*reg15; reg4=(*f.m).resolution*reg35;
    reg9=(*f.m).resolution*reg23; reg28=(*f.m).resolution*reg24; reg37=(*f.m).resolution*reg5; reg45=reg45/reg7; reg13=reg13/reg7;
    reg50=(*f.m).resolution*reg34; reg26=reg26/reg7; reg0=reg0/reg7; reg25=reg74+reg25; reg6=reg130-reg6;
    reg84=reg84/reg18; reg18=reg8/reg18; reg21=reg21*(*f.m).deltaT; reg60=reg60*(*f.m).deltaT; reg64=reg22+reg64;
    reg8=reg30*reg0; reg22=reg25*reg60; reg114=reg15*reg114; reg55=reg6*reg21; reg15=reg95*reg15;
    reg64=reg64*(*f.m).deltaT; reg57=var_inter[1]*reg13; reg59=reg16*reg45; reg61=var_inter[1]*reg0; reg62=(*f.m).resolution*reg84;
    reg65=reg16*reg26; reg66=reg30*reg13; reg67=var_inter[0]*reg26; reg27=reg37+reg27; reg37=var_inter[0]*reg45;
    reg50=reg12-reg50; reg4=reg129-reg4; reg9=reg14+reg9; reg28=reg40+reg28; reg12=(*f.m).resolution*reg18;
    reg14=reg4*reg60; reg40=reg9*reg21; reg69=reg61-reg37; reg70=reg27*reg60; reg71=reg50*reg21;
    reg73=reg67-reg57; reg74=reg59+reg61; reg75=reg65+reg57; reg76=reg28*reg64; reg114=reg12+reg114;
    reg12=reg8+reg37; reg62=reg15-reg62; reg15=reg66-reg65; reg77=reg66+reg67; reg78=reg22+reg55;
    reg79=reg59-reg8; reg80=0.5*reg79; reg81=0.5*reg15; reg82=0.5*reg69; reg83=0.5*reg73;
    reg85=reg14+reg40; reg86=0.5*reg74; reg87=0.5*reg75; reg88=reg71+reg70; reg89=reg78+reg76;
    reg90=reg114*reg64; reg91=0.5*reg12; reg92=0.5*reg77; reg93=reg62*reg64; reg94=reg73*reg6;
    reg95=reg88+reg90; reg96=reg80*reg28; reg97=reg15*reg6; reg98=2*reg89; reg99=reg82*reg28;
    reg100=reg85+reg93; reg101=reg79*reg25; reg102=reg81*reg28; reg103=reg12*reg25; reg104=reg92*reg28;
    reg105=reg28*reg86; reg106=reg75*reg6; reg107=reg83*reg28; reg108=reg69*reg25; reg109=reg87*reg28;
    reg110=reg74*reg25; reg111=reg91*reg28; reg112=reg77*reg6; reg113=reg91*reg62; reg115=reg69*reg4;
    reg116=reg77*reg9; reg117=reg92*reg62; reg118=reg12*reg4; reg119=reg16*var_inter[1]; reg120=reg30*var_inter[0];
    reg121=reg80*reg62; reg122=reg12*reg27; reg123=reg92*reg114; reg111=reg111-reg112; reg124=reg77*reg50;
    reg125=reg91*reg114; reg107=reg108+reg107; reg108=reg74*reg95; reg126=reg75*reg100; reg127=reg98*reg86;
    reg128=reg15*reg9; reg129=reg87*reg98; reg130=reg77*reg100; reg131=reg91*reg98; reg103=reg103-reg104;
    T reg132=reg80*reg114; T reg133=reg15*reg50; reg96=reg97+reg96; reg97=reg92*reg98; T reg134=reg12*reg95;
    T reg135=reg81*reg114; T reg136=reg62*reg86; T reg137=reg75*reg9; T reg138=reg87*reg62; T reg139=reg74*reg4;
    T reg140=reg82*reg62; T reg141=reg73*reg9; T reg142=reg83*reg62; T reg143=reg75*reg50; T reg144=reg74*reg27;
    T reg145=reg87*reg114; reg99=reg94+reg99; reg94=reg73*reg50; reg102=reg101+reg102; reg101=reg83*reg114;
    T reg146=reg69*reg27; reg106=reg106-reg105; reg109=reg109-reg110; T reg147=reg114*reg86; T reg148=reg82*reg114;
    T reg149=reg79*reg27; T reg150=reg81*reg62; T reg151=reg79*reg4; reg103=2*reg103; T reg152=reg15*reg100;
    reg132=reg133+reg132; reg143=reg143-reg147; reg99=2*reg99; reg133=reg83*reg98; reg148=reg94+reg148;
    reg96=2*reg96; reg94=reg69*reg95; reg149=reg135+reg149; reg101=reg146+reg101; reg135=reg80*reg98;
    reg138=reg138-reg139; reg146=reg82*reg98; T reg153=reg81*reg98; T reg154=reg73*reg100; reg150=reg151+reg150;
    reg151=reg79*reg95; reg137=reg137-reg136; reg106=2*reg106; T reg155=reg130-reg131; reg102=2*reg102;
    reg118=reg118-reg117; T reg156=reg108-reg129; T reg157=reg127-reg126; reg145=reg145-reg144; T reg158=reg119*elem.f_vol_e[1];
    T reg159=reg119*elem.f_vol_e[0]; reg113=reg113-reg116; T reg160=reg97-reg134; reg142=reg115+reg142; reg125=reg125-reg124;
    reg109=2*reg109; reg115=reg120*elem.f_vol_e[1]; reg107=2*reg107; T reg161=var_inter[0]*var_inter[1]; reg140=reg141+reg140;
    reg111=2*reg111; reg141=reg30*reg16; reg121=reg128+reg121; reg122=reg122-reg123; reg128=reg120*elem.f_vol_e[0];
    T reg162=reg107*reg91; T reg163=reg91*reg102; T reg164=reg150*reg77; T reg165=reg96*reg91; T reg166=reg77*reg121;
    T reg167=reg75*reg113; T reg168=reg77*reg113; T reg169=reg77*reg118; T reg170=reg111*reg91; reg157=reg157-reg158;
    T reg171=reg145*reg79; reg156=reg156-reg159; T reg172=reg109*reg81; T reg173=reg125*reg79; T reg174=reg111*reg81;
    T reg175=reg122*reg79; T reg176=reg148*reg79; T reg177=reg103*reg81; T reg178=reg141*elem.f_vol_e[0]; T reg179=reg132*reg79;
    T reg180=reg99*reg81; T reg181=reg96*reg81; T reg182=reg149*reg79; T reg183=reg101*reg79; T reg184=reg106*reg80;
    T reg185=reg15*reg138; T reg186=reg109*reg80; T reg187=reg145*reg12; T reg188=reg109*reg92; T reg189=reg143*reg12;
    T reg190=reg106*reg92; T reg191=reg107*reg82; T reg192=reg73*reg142; T reg193=reg99*reg82; T reg194=reg73*reg140;
    T reg195=reg109*reg82; T reg196=reg73*reg138; T reg197=reg106*reg82; T reg198=reg73*reg137; T reg199=reg87*reg102;
    T reg200=reg149*reg74; T reg201=reg96*reg87; T reg202=reg132*reg74; T reg203=reg103*reg87; T reg204=reg122*reg74;
    T reg205=reg111*reg87; T reg206=reg125*reg74; T reg207=reg107*reg81; T reg208=reg15*reg140; T reg209=reg99*reg80;
    T reg210=reg15*reg142; T reg211=reg161*elem.f_vol_e[0]; T reg212=reg161*elem.f_vol_e[1]; T reg213=reg107*reg80; T reg214=reg15*reg113;
    T reg215=reg111*reg80; T reg216=reg15*reg118; T reg217=reg103*reg80; T reg218=reg103*reg91; T reg219=reg15*reg121;
    T reg220=reg81*reg102; T reg221=reg141*elem.f_vol_e[1]; T reg222=reg77*reg142; T reg223=reg99*reg91; T reg224=reg77*reg140;
    T reg225=reg109*reg91; T reg226=reg77*reg138; T reg227=reg106*reg91; T reg228=reg77*reg137; T reg229=reg83*reg102;
    T reg230=reg149*reg69; T reg231=reg96*reg83; T reg232=reg132*reg69; T reg233=reg103*reg83; T reg234=reg122*reg69;
    T reg235=reg111*reg83; T reg236=reg125*reg69; T reg237=reg107*reg83; T reg238=reg101*reg69; T reg239=reg99*reg83;
    T reg240=reg109*reg83; T reg241=reg145*reg69; T reg242=reg75*reg118; T reg243=reg103*reg86; T reg244=reg106*reg83;
    T reg245=reg143*reg69; T reg246=reg82*reg102; T reg247=reg150*reg73; T reg248=reg96*reg82; T reg249=reg73*reg121;
    T reg250=reg103*reg82; reg118=reg73*reg118; T reg251=reg111*reg82; reg113=reg73*reg113; reg121=reg75*reg121;
    T reg252=reg96*reg86; T reg253=reg111*reg86; T reg254=reg135+reg152; T reg255=reg143*reg79; reg142=reg75*reg142;
    reg122=reg122*reg12; reg103=reg103*reg92; T reg256=reg107*reg86; T reg257=reg146+reg154; reg140=reg75*reg140;
    reg125=reg125*reg12; reg111=reg111*reg92; T reg258=reg99*reg86; T reg259=reg153+reg151; reg138=reg75*reg138;
    T reg260=reg109*reg86; T reg261=reg80*reg102; T reg262=reg101*reg12; T reg263=reg107*reg92; T reg264=reg150*reg15;
    T reg265=reg75*reg137; T reg266=reg106*reg86; T reg267=reg96*reg80; T reg268=reg148*reg12; T reg269=reg99*reg92;
    reg107=reg107*reg87; reg101=reg101*reg74; reg99=reg99*reg87; T reg270=reg148*reg74; reg137=reg15*reg137;
    reg155=reg155-reg115; reg109=reg109*reg87; reg145=reg145*reg74; T reg271=reg106*reg81; reg149=reg149*reg12;
    T reg272=reg92*reg102; reg106=reg106*reg87; reg143=reg143*reg74; reg160=reg160-reg128; T reg273=reg133+reg94;
    reg132=reg132*reg12; reg150=reg150*reg75; reg102=reg102*reg86; reg96=reg96*reg92; reg148=reg148*reg69;
    reg226=reg225-reg226; reg206=reg205-reg206; reg155=reg7*reg155; reg241=reg240+reg241; reg205=reg221+reg254;
    reg264=reg261+reg264; reg228=reg227-reg228; reg168=reg170-reg168; reg238=reg237+reg238; reg230=reg229+reg230;
    reg236=reg235+reg236; reg170=reg178+reg259; reg222=reg162-reg222; reg160=reg7*reg160; reg234=reg233+reg234;
    reg224=reg223-reg224; reg255=reg271+reg255; reg232=reg231+reg232; reg101=reg107-reg101; reg270=reg99-reg270;
    reg219=reg267+reg219; reg145=reg109-reg145; reg216=reg217+reg216; reg143=reg106-reg143; reg214=reg215+reg214;
    reg102=reg150-reg102; reg210=reg213+reg210; reg252=reg121-reg252; reg208=reg209+reg208; reg253=reg167-reg253;
    reg204=reg203-reg204; reg202=reg201-reg202; reg200=reg199-reg200; reg256=reg142-reg256; reg198=reg197+reg198;
    reg196=reg195+reg196; reg258=reg140-reg258; reg194=reg193+reg194; reg192=reg191+reg192; reg113=reg251+reg113;
    reg260=reg138-reg260; reg118=reg250+reg118; reg249=reg248+reg249; reg266=reg265-reg266; reg247=reg246+reg247;
    reg245=reg244+reg245; reg243=reg242-reg243; reg164=reg163-reg164; reg190=reg189-reg190; reg99=reg212+reg257;
    reg188=reg187-reg188; reg269=reg268-reg269; reg185=reg186+reg185; reg263=reg262-reg263; reg183=reg207+reg183;
    reg220=reg182+reg220; reg111=reg125-reg111; reg179=reg181+reg179; reg103=reg122-reg103; reg96=reg132-reg96;
    reg239=reg148+reg239; reg175=reg177+reg175; reg176=reg180+reg176; reg173=reg174+reg173; reg272=reg149-reg272;
    reg137=reg184+reg137; reg156=reg156*reg7; reg171=reg172+reg171; reg157=reg157*reg7; reg106=reg211+reg273;
    reg166=reg165-reg166; reg169=reg218-reg169; reg253=reg7*reg253; reg202=reg7*reg202; reg96=reg7*reg96;
    reg200=reg7*reg200; reg179=reg179*reg7; reg256=reg7*reg256; reg198=reg7*reg198; reg222=reg7*reg222;
    reg103=reg7*reg103; reg196=reg7*reg196; reg111=reg7*reg111; reg194=reg7*reg194; reg258=reg7*reg258;
    reg192=reg7*reg192; reg183=reg7*reg183; reg113=reg7*reg113; reg220=reg220*reg7; reg160=ponderation*reg160;
    reg260=reg7*reg260; reg118=reg7*reg118; reg166=reg7*reg166; reg263=reg7*reg263; reg249=reg7*reg249;
    reg157=ponderation*reg157; reg169=reg7*reg169; reg101=reg7*reg101; reg156=ponderation*reg156; reg219=reg7*reg219;
    reg270=reg7*reg270; reg155=ponderation*reg155; reg137=reg7*reg137; reg145=reg7*reg145; reg216=reg7*reg216;
    reg173=reg173*reg7; reg272=reg7*reg272; reg143=reg7*reg143; reg214=reg7*reg214; reg171=reg7*reg171;
    reg176=reg7*reg176; reg102=reg7*reg102; reg168=reg7*reg168; reg210=reg7*reg210; reg175=reg175*reg7;
    reg252=reg7*reg252; reg243=reg7*reg243; reg208=reg7*reg208; reg239=reg239*reg7; reg204=reg7*reg204;
    reg232=reg7*reg232; reg269=reg7*reg269; reg264=reg7*reg264; reg255=reg7*reg255; reg241=reg7*reg241;
    reg190=reg7*reg190; reg238=reg7*reg238; reg230=reg7*reg230; reg107=reg7*reg99; reg234=reg7*reg234;
    reg206=reg7*reg206; reg109=reg7*reg106; reg236=reg7*reg236; reg121=reg7*reg170; reg122=reg7*reg205;
    reg164=reg7*reg164; reg226=reg7*reg226; reg247=reg7*reg247; reg185=reg7*reg185; reg266=reg7*reg266;
    reg188=reg7*reg188; reg228=reg7*reg228; reg245=reg7*reg245; reg224=reg7*reg224; T tmp_4_2=ponderation*reg234;
    reg125=ponderation*reg122; T vec_1=reg125; reg132=ponderation*reg121; T vec_0=reg132; T tmp_0_7=ponderation*reg255;
    T tmp_1_4=ponderation*reg210; T tmp_0_5=ponderation*reg176; T tmp_4_1=ponderation*reg232; T tmp_3_3=ponderation*reg168; T tmp_6_7=ponderation*reg143;
    T tmp_7_0=ponderation*reg102; T tmp_1_3=ponderation*reg214; T tmp_0_2=ponderation*reg175; T tmp_3_6=ponderation*reg226; T tmp_6_4=ponderation*reg101;
    T vec_7=-reg157; T vec_6=-reg156; T tmp_3_2=ponderation*reg169; T tmp_3_7=ponderation*reg228; T tmp_1_1=ponderation*reg219;
    T tmp_6_5=ponderation*reg270; reg101=ponderation*reg109; T vec_4=reg101; T tmp_2_7=ponderation*reg190; T tmp_1_7=ponderation*reg137;
    T vec_3=-reg155; T tmp_0_3=ponderation*reg173; T tmp_4_0=ponderation*reg230; T tmp_1_2=ponderation*reg216; T tmp_6_6=ponderation*reg145;
    T tmp_7_2=ponderation*reg243; T tmp_3_0=ponderation*reg164; T tmp_2_0=ponderation*reg272; T tmp_3_4=ponderation*reg222; T tmp_1_0=ponderation*reg264;
    T tmp_5_5=ponderation*reg194; T tmp_0_6=ponderation*reg171; T tmp_2_3=ponderation*reg111; T tmp_3_5=ponderation*reg224; T tmp_7_5=ponderation*reg258;
    T tmp_5_4=ponderation*reg192; T tmp_0_4=ponderation*reg183; T tmp_6_3=ponderation*reg206; T tmp_0_0=ponderation*reg220; T tmp_5_3=ponderation*reg113;
    T tmp_4_7=ponderation*reg245; T tmp_7_7=ponderation*reg266; T tmp_5_2=ponderation*reg118; T tmp_7_6=ponderation*reg260; T tmp_5_1=ponderation*reg249;
    T tmp_5_0=ponderation*reg247; T tmp_1_6=ponderation*reg185; T tmp_2_4=ponderation*reg263; T tmp_7_1=ponderation*reg252; T tmp_2_6=ponderation*reg188;
    T tmp_4_5=ponderation*reg239; T tmp_6_2=ponderation*reg204; T tmp_4_3=ponderation*reg236; T vec_2=-reg160; T tmp_7_3=ponderation*reg253;
    T tmp_6_1=ponderation*reg202; reg102=ponderation*reg107; T vec_5=reg102; T tmp_2_1=ponderation*reg96; T tmp_3_1=ponderation*reg166;
    T tmp_6_0=ponderation*reg200; T tmp_4_4=ponderation*reg238; T tmp_0_1=ponderation*reg179; T tmp_2_5=ponderation*reg269; T tmp_5_7=ponderation*reg198;
    T tmp_1_5=ponderation*reg208; T tmp_7_4=ponderation*reg256; T tmp_5_6=ponderation*reg196; T tmp_4_6=ponderation*reg241; T tmp_2_2=ponderation*reg103;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_3;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_4;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_5;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_7;
    matrix(indices[0]+1,indices[0]+0) += tmp_1_0;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_3;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_4;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_5;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_7;
    matrix(indices[1]+0,indices[0]+0) += tmp_2_0;
    matrix(indices[1]+0,indices[0]+1) += tmp_2_1;
    matrix(indices[1]+0,indices[1]+0) += tmp_2_2;
    matrix(indices[1]+0,indices[1]+1) += tmp_2_3;
    matrix(indices[1]+0,indices[2]+0) += tmp_2_4;
    matrix(indices[1]+0,indices[2]+1) += tmp_2_5;
    matrix(indices[1]+0,indices[3]+0) += tmp_2_6;
    matrix(indices[1]+0,indices[3]+1) += tmp_2_7;
    matrix(indices[1]+1,indices[0]+0) += tmp_3_0;
    matrix(indices[1]+1,indices[0]+1) += tmp_3_1;
    matrix(indices[1]+1,indices[1]+0) += tmp_3_2;
    matrix(indices[1]+1,indices[1]+1) += tmp_3_3;
    matrix(indices[1]+1,indices[2]+0) += tmp_3_4;
    matrix(indices[1]+1,indices[2]+1) += tmp_3_5;
    matrix(indices[1]+1,indices[3]+0) += tmp_3_6;
    matrix(indices[1]+1,indices[3]+1) += tmp_3_7;
    matrix(indices[2]+0,indices[0]+0) += tmp_4_0;
    matrix(indices[2]+0,indices[0]+1) += tmp_4_1;
    matrix(indices[2]+0,indices[1]+0) += tmp_4_2;
    matrix(indices[2]+0,indices[1]+1) += tmp_4_3;
    matrix(indices[2]+0,indices[2]+0) += tmp_4_4;
    matrix(indices[2]+0,indices[2]+1) += tmp_4_5;
    matrix(indices[2]+0,indices[3]+0) += tmp_4_6;
    matrix(indices[2]+0,indices[3]+1) += tmp_4_7;
    matrix(indices[2]+1,indices[0]+0) += tmp_5_0;
    matrix(indices[2]+1,indices[0]+1) += tmp_5_1;
    matrix(indices[2]+1,indices[1]+0) += tmp_5_2;
    matrix(indices[2]+1,indices[1]+1) += tmp_5_3;
    matrix(indices[2]+1,indices[2]+0) += tmp_5_4;
    matrix(indices[2]+1,indices[2]+1) += tmp_5_5;
    matrix(indices[2]+1,indices[3]+0) += tmp_5_6;
    matrix(indices[2]+1,indices[3]+1) += tmp_5_7;
    matrix(indices[3]+0,indices[0]+0) += tmp_6_0;
    matrix(indices[3]+0,indices[0]+1) += tmp_6_1;
    matrix(indices[3]+0,indices[1]+0) += tmp_6_2;
    matrix(indices[3]+0,indices[1]+1) += tmp_6_3;
    matrix(indices[3]+0,indices[2]+0) += tmp_6_4;
    matrix(indices[3]+0,indices[2]+1) += tmp_6_5;
    matrix(indices[3]+0,indices[3]+0) += tmp_6_6;
    matrix(indices[3]+0,indices[3]+1) += tmp_6_7;
    matrix(indices[3]+1,indices[0]+0) += tmp_7_0;
    matrix(indices[3]+1,indices[0]+1) += tmp_7_1;
    matrix(indices[3]+1,indices[1]+0) += tmp_7_2;
    matrix(indices[3]+1,indices[1]+1) += tmp_7_3;
    matrix(indices[3]+1,indices[2]+0) += tmp_7_4;
    matrix(indices[3]+1,indices[2]+1) += tmp_7_5;
    matrix(indices[3]+1,indices[3]+0) += tmp_7_6;
    matrix(indices[3]+1,indices[3]+1) += tmp_7_7;
sollicitation[indices[0]+0] += vec_0;
sollicitation[indices[0]+1] += vec_1;
sollicitation[indices[1]+0] += vec_2;
sollicitation[indices[1]+1] += vec_3;
sollicitation[indices[2]+0] += vec_4;
sollicitation[indices[2]+1] += vec_5;
sollicitation[indices[3]+0] += vec_6;
sollicitation[indices[3]+1] += vec_7;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=2*(*f.m).shear_modulus_13; T reg1=2*(*f.m).shear_modulus_23; T reg2=2*(*f.m).shear_modulus_12; reg0=1.0/reg0; reg1=1.0/reg1;
    T reg3=reg0*reg1; reg2=1.0/reg2; T reg4=reg2*reg3; T reg5=pow((*f.m).v1[0],2); T reg6=pow((*f.m).v1[1],2);
    T reg7=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg8=1.0/(*f.m).elastic_modulus_3; T reg9=pow((*f.m).v2[1],2); T reg10=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg11=pow((*f.m).v2[0],2);
    T reg12=pow((*f.m).v1[2],2); reg6=reg5+reg6; reg5=1.0/(*f.m).elastic_modulus_2; T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg14=reg8*reg4;
    T reg15=reg7*reg4; T reg16=pow((*f.m).v2[2],2); reg9=reg11+reg9; reg11=reg10*reg4; reg12=reg6+reg12;
    reg16=reg9+reg16; reg6=reg5*reg14; reg9=reg7*reg15; T reg17=reg13*reg14; T reg18=reg7*reg11;
    T reg19=1.0/(*f.m).elastic_modulus_1; reg16=pow(reg16,0.5); T reg20=reg5*reg11; reg12=pow(reg12,0.5); reg9=reg6-reg9;
    reg6=reg13*reg15; reg18=reg17+reg18; T reg21=reg6+reg20; T reg22=(*f.m).v1[1]/reg12; T reg23=(*f.m).v1[2]/reg12;
    T reg24=reg13*reg18; T reg25=reg19*reg9; T reg26=(*f.m).v2[2]/reg16; T reg27=(*f.m).v2[1]/reg16; T reg28=reg2*reg1;
    T reg29=reg8*reg3; reg12=(*f.m).v1[0]/reg12; T reg30=reg10*reg3; reg3=reg7*reg3; reg24=reg25-reg24;
    reg25=reg10*reg15; T reg31=reg5*reg4; T reg32=reg23*reg27; T reg33=reg22*reg26; reg16=(*f.m).v2[0]/reg16;
    reg14=reg19*reg14; T reg34=reg10*reg11; reg4=reg13*reg4; T reg35=reg10*reg21; T reg36=2*reg12;
    T reg37=reg7*reg28; T reg38=reg33-reg32; reg11=reg13*reg11; T reg39=reg10*reg28; T reg40=reg2*reg0;
    T reg41=reg7*reg30; reg28=reg8*reg28; T reg42=reg23*reg16; T reg43=reg12*reg26; T reg44=reg10*reg31;
    T reg45=2*reg16; reg35=reg24-reg35; reg24=reg5*reg29; reg25=reg17+reg25; reg17=reg7*reg3;
    reg29=reg13*reg29; reg34=reg14-reg34; reg14=reg10*reg4; reg15=reg19*reg15; reg30=reg5*reg30;
    T reg46=reg10*reg40; reg3=reg13*reg3; reg4=reg13*reg4; reg11=reg15+reg11; T reg47=reg7*reg40;
    reg14=reg15+reg14; reg9=reg9/reg35; reg31=reg19*reg31; reg25=reg25/reg35; reg44=reg6+reg44;
    reg18=reg18/reg35; reg34=reg34/reg35; reg15=reg7*reg39; T reg48=reg42-reg43; T reg49=reg45*reg27;
    T reg50=pow(reg12,2); T reg51=reg7*reg37; T reg52=reg13*reg28; T reg53=pow(reg22,2); T reg54=reg22*reg36;
    reg28=reg5*reg28; T reg55=reg12*reg27; reg41=reg29+reg41; reg40=reg8*reg40; reg29=2*reg38;
    T reg56=pow(reg16,2); reg17=reg24-reg17; reg24=reg22*reg16; T reg57=pow(reg27,2); T reg58=reg9*reg53;
    T reg59=reg18*reg57; T reg60=pow(reg48,2); T reg61=reg7*reg46; T reg62=reg9*reg50; T reg63=pow(reg38,2);
    T reg64=reg9*reg54; T reg65=reg18*reg49; reg21=reg21/reg35; T reg66=reg25*reg50; T reg67=reg34*reg56;
    T reg68=reg25*reg53; T reg69=reg34*reg57; T reg70=reg34*reg49; T reg71=reg25*reg54; T reg72=pow(reg26,2);
    reg39=reg5*reg39; reg17=reg19*reg17; reg37=reg13*reg37; reg41=reg13*reg41; reg44=reg44/reg35;
    T reg73=reg3+reg30; T reg74=reg29*reg48; T reg75=pow(reg23,2); reg51=reg28-reg51; reg4=reg31-reg4;
    reg28=reg55-reg24; reg15=reg52+reg15; reg11=reg11/reg35; reg31=reg5*reg40; reg40=reg13*reg40;
    reg14=reg14/reg35; reg52=reg7*reg47; T reg76=reg18*reg56; T reg77=pow(reg28,2); T reg78=reg14*reg56;
    T reg79=reg44*reg50; T reg80=reg11*reg74; reg70=reg71+reg70; reg71=reg34*reg72; reg47=reg13*reg47;
    reg46=reg5*reg46; reg41=reg17-reg41; reg73=reg10*reg73; reg17=reg44*reg53; T reg81=reg14*reg57;
    T reg82=reg44*reg54; T reg83=reg14*reg49; reg51=reg19*reg51; reg15=reg13*reg15; reg4=reg4/reg35;
    T reg84=reg37+reg39; reg52=reg31-reg52; reg61=reg40+reg61; reg59=reg58+reg59; reg31=reg21*reg60;
    reg40=reg9*reg75; reg58=reg18*reg72; reg65=reg64+reg65; reg64=reg21*reg74; reg67=reg66+reg67;
    reg66=reg11*reg63; reg69=reg68+reg69; reg68=reg11*reg60; T reg85=reg25*reg75; T reg86=reg21*reg63;
    reg76=reg62+reg76; reg62=reg45*reg26; reg80=reg70+reg80; reg70=reg11*reg77; reg71=reg85+reg71;
    reg68=reg69+reg68; reg66=reg67+reg66; reg64=reg65+reg64; reg65=reg21*reg77; reg58=reg40+reg58;
    reg31=reg59+reg31; reg40=reg47+reg46; reg61=reg13*reg61; reg52=reg19*reg52; reg84=reg10*reg84;
    reg15=reg51-reg15; reg73=reg41-reg73; reg78=reg79+reg78; reg41=reg4*reg63; reg81=reg17+reg81;
    reg17=reg4*reg60; reg51=reg44*reg75; reg59=reg14*reg72; reg83=reg82+reg83; reg74=reg4*reg74;
    reg67=2*reg27; reg69=reg23*reg36; reg86=reg76+reg86; reg76=2*reg22; reg79=reg16*reg27;
    reg82=reg12*reg22; reg85=reg79*reg68; reg73=reg73/reg35; T reg87=reg80*reg56; T reg88=reg64*reg50;
    T reg89=reg82*reg31; T reg90=reg86*reg53; reg41=reg78+reg41; reg78=reg66*reg57; T reg91=reg31*reg53;
    reg17=reg81+reg17; reg81=reg68*reg57; T reg92=reg64*reg53; reg59=reg51+reg59; reg51=reg4*reg77;
    T reg93=reg80*reg57; reg74=reg83+reg74; reg83=reg86*reg50; T reg94=reg66*reg56; T reg95=reg67*reg26;
    T reg96=reg34*reg62; T reg97=reg25*reg69; reg70=reg71+reg70; reg29=reg29*reg28; reg71=2*reg48;
    T reg98=reg12*reg16; T reg99=reg22*reg27; T reg100=reg18*reg62; T reg101=reg9*reg69; reg55=reg24+reg55;
    reg65=reg58+reg65; reg24=reg22*reg38; reg58=reg12*reg48; reg40=reg10*reg40; reg61=reg52-reg61;
    reg52=reg23*reg76; reg84=reg15-reg84; reg80=reg79*reg80; reg64=reg82*reg64; reg15=reg48*reg38;
    reg31=reg31*reg50; T reg102=reg27*reg76; reg68=reg68*reg56; T reg103=reg16*reg36; T reg104=2*reg23;
    T reg105=reg73*reg55; T reg106=reg22*reg48; T reg107=reg65*reg50; reg51=reg59+reg51; reg59=reg12*reg38;
    reg100=reg101+reg100; reg101=reg21*reg29; T reg108=reg23*reg26; reg66=reg79*reg66; reg86=reg82*reg86;
    T reg109=reg13*reg56; reg84=reg84/reg35; T reg110=reg102*reg5; T reg111=reg70*reg56; T reg112=reg73*reg99;
    T reg113=reg73*reg98; reg85=reg89+reg85; reg89=reg15*reg17; reg104=reg104*reg26; reg58=reg24+reg58;
    reg24=reg13*reg57; T reg114=reg27*reg38; T reg115=reg14*reg62; T reg116=reg44*reg69; T reg117=reg16*reg48;
    reg40=reg61-reg40; reg80=reg64+reg80; reg61=reg15*reg74; reg64=reg102*reg13; T reg118=reg103*reg19;
    reg25=reg25*reg52; T reg119=reg63*reg17; T reg120=reg70*reg57; T reg121=reg65*reg53; reg34=reg34*reg95;
    reg68=reg31+reg68; reg93=reg92+reg93; reg31=reg60*reg74; reg18=reg18*reg95; reg9=reg9*reg52;
    reg92=reg63*reg41; reg94=reg83+reg94; reg17=reg60*reg17; reg81=reg91+reg81; reg96=reg97+reg96;
    reg87=reg88+reg87; reg71=reg71*reg28; reg83=reg5*reg57; reg88=reg11*reg29; reg91=reg103*reg13;
    reg74=reg63*reg74; reg97=reg60*reg41; reg78=reg90+reg78; reg90=reg19*reg56; T reg122=reg27*reg48;
    reg24=reg90-reg24; reg117=reg114+reg117; reg18=reg9+reg18; reg21=reg21*reg71; reg9=reg16*reg38;
    reg90=reg10*reg72; reg88=reg96+reg88; reg96=reg10*reg56; reg34=reg25+reg34; reg25=reg23*reg28;
    reg11=reg11*reg71; reg36=reg38*reg36; reg101=reg100+reg101; reg76=reg48*reg76; reg100=reg63*reg51;
    reg91=reg110-reg91; reg111=reg107+reg111; reg107=reg104*reg7; reg66=reg86+reg66; reg86=reg102*reg7;
    reg110=reg102*reg105; reg114=reg7*reg57; reg31=reg93+reg31; reg97=reg78+reg97; reg78=reg102*reg113;
    reg93=reg7*reg72; reg74=reg87+reg74; reg87=reg103*reg112; reg119=reg68+reg119; reg92=reg94+reg92;
    reg68=reg103*reg113; reg17=reg81+reg17; reg81=reg102*reg112; reg94=reg60*reg51; T reg123=reg103*reg105;
    reg120=reg121+reg120; reg105=reg55*reg105; reg61=reg80+reg61; reg80=reg104*reg10; reg115=reg116+reg115;
    reg29=reg4*reg29; reg70=reg79*reg70; reg44=reg44*reg52; reg14=reg14*reg95; reg65=reg82*reg65;
    reg64=reg118-reg64; reg116=reg103*reg10; reg112=reg55*reg112; reg89=reg85+reg89; reg85=reg73*reg108;
    reg109=reg83-reg109; reg35=reg40/reg35; reg59=reg84*reg59; reg106=reg84*reg106; reg41=reg15*reg41;
    reg40=reg84*reg58; reg83=reg103*reg85; reg123=reg74+reg123; reg74=reg36*reg40; reg118=reg101*reg50;
    reg121=reg88*reg56; reg5=reg5*reg53; T reg124=reg13*reg50; reg80=reg64-reg80; reg90=reg24-reg90;
    reg113=reg55*reg113; reg41=reg66+reg41; reg24=(*f.m).alpha_2*reg57; reg64=(*f.m).alpha_1*reg53; reg66=(*f.m).alpha_2*reg56;
    reg13=reg13*reg53; reg87=reg119+reg87; reg78=reg97+reg78; reg97=reg76*reg59; reg81=reg17+reg81;
    reg17=reg76*reg106; reg94=reg120+reg94; reg119=reg102*reg85; reg110=reg31+reg110; reg31=reg76*reg40;
    reg120=reg101*reg53; T reg125=reg88*reg57; reg21=reg18+reg21; reg18=reg36*reg59; reg11=reg34+reg11;
    reg34=reg36*reg106; reg45=reg45*reg38; reg67=reg67*reg48; reg42=reg43+reg42; reg43=reg23*reg38;
    T reg126=reg12*reg28; T reg127=reg26*reg28; reg40=reg58*reg40; reg105=reg61+reg105; reg51=reg15*reg51;
    reg70=reg65+reg70; reg106=reg58*reg106; reg112=reg89+reg112; reg104=reg104*reg8; reg61=reg35*reg117;
    reg86=reg116+reg86; reg122=reg35*reg122; reg14=reg44+reg14; reg19=reg19*reg50; reg114=reg96+reg114;
    reg107=reg91-reg107; reg29=reg115+reg29; reg93=reg109-reg93; reg100=reg111+reg100; reg71=reg4*reg71;
    reg25=reg84*reg25; reg68=reg92+reg68; reg4=reg8*reg72; reg9=reg35*reg9; reg44=(*f.m).alpha_1*reg50;
    reg65=reg67*reg122; reg17=reg81+reg17; reg33=reg32+reg33; reg13=reg19-reg13; reg19=reg10*reg75;
    reg106=reg112+reg106; reg32=reg117*reg122; reg97=reg78+reg97; reg78=reg67*reg9; reg126=reg43+reg126;
    reg43=reg90*reg98; reg81=reg93*reg53; reg89=reg90*reg50; reg91=reg80*reg56; reg92=reg107*reg57;
    reg96=reg73*reg42; reg113=reg41+reg113; reg59=reg58*reg59; reg34=reg87+reg34; reg122=reg45*reg122;
    reg41=reg107*reg53; reg87=reg80*reg50; reg71=reg14+reg71; reg14=reg45*reg9; reg18=reg68+reg18;
    reg68=(*f.m).alpha_2*reg72; reg109=(*f.m).alpha_1*reg75; reg111=reg60*(*f.m).alpha_3; reg90=reg90*reg56; reg112=reg93*reg57;
    reg24=reg64+reg24; reg10=reg10*reg50; reg127=reg35*reg127; reg64=reg7*reg53; reg115=reg63*(*f.m).alpha_3;
    reg66=reg44+reg66; reg44=reg16*reg28; reg31=reg110+reg31; reg110=reg67*reg61; reg116=reg63*reg29;
    reg51=reg70+reg51; reg121=reg118+reg121; reg85=reg55*reg85; reg86=reg104-reg86; reg70=reg45*reg61;
    reg74=reg123+reg74; reg125=reg120+reg125; reg104=reg60*reg29; reg88=reg79*reg88; reg118=reg21*reg53;
    reg101=reg82*reg101; reg120=reg11*reg57; reg61=reg117*reg61; reg40=reg105+reg40; reg105=reg36*reg25;
    reg83=reg100+reg83; reg93=reg93*reg99; reg7=reg7*reg75; reg100=reg76*reg25; reg119=reg94+reg119;
    reg107=reg107*reg99; reg80=reg80*reg98; reg114=reg4-reg114; reg4=reg11*reg56; reg124=reg5-reg124;
    reg5=reg22*reg28; reg94=reg23*reg48; reg123=reg21*reg50; T reg128=reg26*reg38; reg59=reg113+reg59;
    reg9=reg117*reg9; reg25=reg58*reg25; reg113=reg114*reg108; reg5=reg94+reg5; reg94=reg86*reg75;
    reg41=reg87+reg41; reg93=reg43+reg93; reg16=reg16*reg26; reg61=reg40+reg61; reg40=reg86*reg108;
    reg88=reg101+reg88; reg29=reg15*reg29; reg85=reg51+reg85; reg107=reg80+reg107; reg32=reg106+reg32;
    reg12=reg12*reg23; reg64=reg10+reg64; reg122=reg34+reg122; reg10=reg27*reg28; reg34=reg26*reg48;
    reg8=reg8*reg75; reg44=reg128+reg44; reg19=reg13-reg19; reg11=reg79*reg11; reg13=reg114*reg75;
    reg78=reg97+reg78; reg81=reg89+reg81; reg92=reg91+reg92; reg86=reg86*reg72; reg65=reg17+reg65;
    reg17=reg102*reg96; reg104=reg125+reg104; reg100=reg119+reg100; reg43=reg67*reg127; reg70=reg74+reg70;
    reg63=reg63*reg71; reg4=reg123+reg4; reg110=reg31+reg110; reg7=reg124-reg7; reg31=reg2*reg55;
    reg51=reg79*reg2; reg74=reg103*reg96; reg116=reg121+reg116; reg73=reg73*reg33; reg126=reg84*reg126;
    reg60=reg60*reg71; reg105=reg83+reg105; reg79=reg79*(*f.m).alpha_2; reg80=reg82*(*f.m).alpha_1; reg83=reg45*reg127;
    reg120=reg118+reg120; reg77=reg77*(*f.m).alpha_3; reg14=reg18+reg14; reg68=reg109+reg68; reg111=reg24+reg111;
    reg112=reg90+reg112; reg114=reg114*reg72; reg115=reg66+reg115; reg21=reg82*reg21; reg17=reg104+reg17;
    reg25=reg85+reg25; reg127=reg117*reg127; reg38=reg28*reg38; reg18=reg76*reg126; reg60=reg120+reg60;
    reg102=reg102*reg73; reg5=reg84*reg5; reg24=reg16*(*f.m).alpha_2; reg66=reg12*(*f.m).alpha_1; reg9=reg59+reg9;
    reg26=reg27*reg26; reg27=reg15*(*f.m).alpha_3; reg79=reg80+reg79; reg23=reg22*reg23; reg77=reg68+reg77;
    reg44=reg35*reg44; reg103=reg103*reg73; reg22=reg51*reg54; reg13=reg81+reg13; reg59=reg7*reg53;
    reg68=reg19*reg50; reg43=reg100+reg43; reg80=reg0*reg42; reg16=reg16*reg0; reg10=reg34+reg10;
    reg11=reg21+reg11; reg94=reg41+reg94; reg21=reg31*reg54; reg71=reg15*reg71; reg74=reg116+reg74;
    reg15=reg36*reg126; reg34=reg65*reg111; reg41=reg78*reg115; reg81=reg65*reg61; reg84=reg19*reg56;
    reg85=reg7*reg57; reg87=reg61*reg122; reg114=reg112+reg114; reg89=reg51*reg49; reg90=reg110*reg32;
    reg63=reg4+reg63; reg4=reg111*reg122; reg91=reg115*reg14; reg97=reg32*reg70; reg100=reg31*reg49;
    reg86=reg92+reg86; reg113=reg93+reg113; reg51=reg51*reg55; reg29=reg88+reg29; reg40=reg107+reg40;
    reg31=reg31*reg55; reg83=reg105+reg83; reg64=reg8-reg64; reg96=reg55*reg96; reg126=reg58*reg126;
    reg8=reg67*reg44; reg83=reg77*reg83; reg4=reg91+reg4; reg88=reg16*reg42; reg51=reg113+reg51;
    reg35=reg10*reg35; reg90=reg81-reg90; reg76=reg76*reg5; reg10=reg80*reg62; reg81=reg26*reg1;
    reg63=reg103+reg63; reg36=reg36*reg5; reg91=reg1*reg33; reg102=reg60+reg102; reg48=reg28*reg48;
    reg100=reg86+reg100; reg18=reg17+reg18; reg17=reg16*reg69; reg22=reg13+reg22; reg97=reg87-reg97;
    reg89=reg114+reg89; reg16=reg16*reg62; reg13=reg110*reg122; reg127=reg25+reg127; reg59=reg68+reg59;
    reg75=reg64*reg75; reg25=reg65*reg70; reg28=reg111*reg32; reg60=reg115*reg9; reg19=reg19*reg98;
    reg71=reg11+reg71; reg73=reg55*reg73; reg7=reg7*reg99; reg26=reg26*(*f.m).alpha_2; reg11=reg23*(*f.m).alpha_1;
    reg15=reg74+reg15; reg38=reg38*(*f.m).alpha_3; reg24=reg66+reg24; reg43=reg43*reg77; reg34=reg41+reg34;
    reg27=reg79+reg27; reg96=reg29+reg96; reg29=reg80*reg69; reg80=reg80*reg42; reg2=reg82*reg2;
    reg31=reg40+reg31; reg40=reg45*reg44; reg85=reg84+reg85; reg21=reg94+reg21; reg72=reg64*reg72;
    reg40=reg15+reg40; reg75=reg59+reg75; reg126=reg96+reg126; reg0=reg12*reg0; reg73=reg71+reg73;
    reg5=reg58*reg5; reg8=reg18+reg8; reg12=reg91*reg95; reg10=reg100+reg10; reg15=reg90*reg14;
    reg44=reg117*reg44; reg18=reg78*reg97; reg41=reg27*reg70; reg58=reg2*reg49; reg72=reg85+reg72;
    reg59=reg91*reg52; reg29=reg21+reg29; reg38=reg24+reg38; reg26=reg11+reg26; reg48=reg48*(*f.m).alpha_3;
    reg43=reg34+reg43; reg11=reg110*reg27; reg28=reg60+reg28; reg127=reg77*reg127; reg7=reg19+reg7;
    reg91=reg91*reg33; reg80=reg31+reg80; reg108=reg64*reg108; reg19=reg81*reg33; reg88=reg51+reg88;
    reg21=reg2*reg54; reg25=reg13-reg25; reg13=reg81*reg95; reg16=reg89+reg16; reg17=reg22+reg17;
    reg81=reg81*reg52; reg45=reg45*reg35; reg36=reg63+reg36; reg83=reg4+reg83; reg76=reg102+reg76;
    reg67=reg67*reg35; reg18=reg15-reg18; reg4=reg9*reg25; reg41=reg83+reg41; reg44=reg126+reg44;
    reg15=reg9*reg70; reg67=reg76+reg67; reg22=reg27*reg61; reg127=reg28+reg127; reg24=reg78*reg61;
    reg40=reg38*reg40; reg8=reg8*reg38; reg11=reg43+reg11; reg48=reg26+reg48; reg5=reg73+reg5;
    reg35=reg117*reg35; reg26=reg61*reg14; reg28=reg110*reg9; reg81=reg17+reg81; reg59=reg29+reg59;
    reg58=reg72+reg58; reg62=reg0*reg62; reg13=reg16+reg13; reg12=reg10+reg12; reg1=reg23*reg1;
    reg21=reg75+reg21; reg69=reg0*reg69; reg108=reg7+reg108; reg2=reg2*reg55; reg19=reg88+reg19;
    reg91=reg80+reg91; reg45=reg36+reg45; reg52=reg1*reg52; reg2=reg108+reg2; reg42=reg0*reg42;
    reg0=reg59*reg19; reg35=reg5+reg35; reg4=reg18+reg4; reg44=reg38*reg44; reg22=reg127+reg22;
    reg5=reg78*reg32; reg28=reg24-reg28; reg7=reg65*reg9; reg10=reg32*reg14; reg15=reg26-reg15;
    reg16=reg9*reg122; reg17=reg110*reg14; reg18=reg78*reg70; reg40=reg41+reg40; reg45=reg45*reg48;
    reg67=reg67*reg48; reg8=reg11+reg8; reg11=1-var_inter[1]; reg23=1-var_inter[0]; reg62=reg58+reg62;
    reg95=reg1*reg95; reg24=reg81*reg91; reg26=reg13*reg91; reg29=reg12*reg19; reg69=reg21+reg69;
    reg18=reg17-reg18; reg17=reg65*reg14; reg21=reg59*reg13; reg16=reg10-reg16; reg29=reg26-reg29;
    reg42=reg2+reg42; reg33=reg1*reg33; reg15=reg15/reg4; reg97=reg97/reg4; reg7=reg5-reg7;
    reg95=reg62+reg95; reg28=reg28/reg4; reg90=reg90/reg4; reg44=reg22+reg44; reg1=reg11*elem.pos(1)[0];
    reg2=reg11*elem.pos(0)[0]; reg5=var_inter[0]*elem.pos(1)[0]; reg10=reg11*elem.pos(1)[1]; reg22=reg11*elem.pos(0)[1]; reg35=reg48*reg35;
    reg26=reg81*reg12; reg67=reg8+reg67; reg0=reg24-reg0; reg52=reg69+reg52; reg45=reg40+reg45;
    reg8=reg23*elem.pos(0)[0]; reg24=var_inter[0]*elem.pos(1)[1]; reg31=reg23*elem.pos(0)[1]; reg34=reg78*reg122; reg36=var_inter[1]*elem.pos(2)[0];
    reg35=reg44+reg35; reg90=reg90*reg45; reg38=reg95*reg0; reg40=reg29*reg52; reg97=reg97*reg67;
    reg1=reg1-reg2; reg28=reg28*reg45; reg15=reg15*reg67; reg41=var_inter[0]*elem.pos(2)[0]; reg18=reg18/reg4;
    reg21=reg26-reg21; reg25=reg25/reg4; reg26=reg24+reg31; reg43=var_inter[0]*elem.pos(2)[1]; reg33=reg42+reg33;
    reg34=reg17-reg34; reg16=reg16/reg4; reg17=var_inter[1]*elem.pos(2)[1]; reg42=reg8+reg5; reg7=reg7/reg4;
    reg10=reg10-reg22; reg44=var_inter[1]*elem.pos(3)[0]; reg48=reg52*reg19; reg1=reg36+reg1; reg17=reg10+reg17;
    reg10=reg13*reg33; reg36=reg81*reg33; reg51=var_inter[1]*elem.pos(3)[1]; reg19=reg95*reg19; reg4=reg34/reg4;
    reg34=reg23*elem.pos(3)[1]; reg38=reg40-reg38; reg43=reg43-reg26; reg40=reg23*elem.pos(3)[0]; reg41=reg41-reg42;
    reg58=reg33*reg21; reg67=reg16*reg67; reg45=reg7*reg45; reg28=reg15-reg28; reg18=reg18*reg35;
    reg97=reg90-reg97; reg25=reg25*reg35; reg7=1-(*f.m).resolution; reg58=reg38+reg58; reg81=reg81*reg95;
    reg1=reg1-reg44; reg13=reg13*reg52; reg36=reg48-reg36; reg15=reg59*reg33; reg97=reg25+reg97;
    reg34=reg43+reg34; reg40=reg41+reg40; reg18=reg28-reg18; reg35=reg4*reg35; reg67=reg45-reg67;
    reg17=reg17-reg51; reg4=reg95*reg91; reg10=reg19-reg10; reg33=reg12*reg33; reg91=reg52*reg91;
    reg52=reg12*reg52; reg18=reg7*reg18; reg12=reg34*reg1; reg16=(*f.m).resolution*reg111; reg15=reg91-reg15;
    reg36=reg36/reg58; reg33=reg4-reg33; reg4=(*f.m).resolution*reg115; reg67=reg35+reg67; reg95=reg59*reg95;
    reg10=reg10/reg58; reg97=reg7*reg97; reg81=reg13-reg81; reg13=reg40*reg17; reg97=reg4+reg97;
    reg18=reg16+reg18; reg33=reg33/reg58; reg13=reg12-reg13; reg32=reg32*reg7; reg9=reg9*reg7;
    reg15=reg15/reg58; reg4=(*f.m).resolution*reg27; reg95=reg52-reg95; reg12=(*f.m).resolution*reg10; reg16=(*f.m).resolution*reg36;
    reg81=reg81/reg58; reg0=reg0/reg58; reg29=reg29/reg58; reg67=reg7*reg67; reg14=reg7*reg14;
    reg34=reg34/reg13; reg122=reg7*reg122; reg1=reg1/reg13; reg78=reg78*reg7; reg65=reg65*reg7;
    reg40=reg40/reg13; reg17=reg17/reg13; reg19=(*f.m).resolution*reg33; reg25=(*f.m).resolution*reg15; reg21=reg21/reg58;
    reg58=reg95/reg58; reg28=(*f.m).resolution*reg81; reg12=reg9+reg12; reg16=reg32-reg16; reg9=(*f.m).resolution*reg0;
    reg61=reg61*reg7; reg67=reg4+reg67; reg4=(*f.m).resolution*reg29; reg18=reg18*(*f.m).deltaT; reg97=reg97*(*f.m).deltaT;
    reg32=reg16*reg18; reg35=reg12*reg97; reg67=reg67*(*f.m).deltaT; reg38=(*f.m).resolution*reg58; reg41=reg11*reg34;
    reg14=reg4+reg14; reg9=reg122-reg9; reg70=reg7*reg70; reg19=reg78-reg19; reg25=reg65+reg25;
    reg28=reg61+reg28; reg7=reg110*reg7; reg4=reg11*reg40; reg43=reg23*reg1; reg45=var_inter[1]*reg40;
    reg48=var_inter[1]*reg34; reg52=reg23*reg17; reg59=var_inter[0]*reg1; reg60=(*f.m).resolution*reg21; reg61=var_inter[0]*reg17;
    reg62=reg48-reg61; reg63=reg59-reg45; reg64=reg52+reg48; reg65=reg43+reg45; reg66=reg4+reg59;
    reg68=reg41+reg61; reg69=reg25*reg18; reg71=reg19*reg97; reg72=reg4-reg43; reg38=reg7-reg38;
    reg7=reg14*reg97; reg73=reg9*reg18; reg74=reg52-reg41; reg70=reg60+reg70; reg60=reg28*reg67;
    reg75=reg35+reg32; reg76=0.5*reg74; reg77=0.5*reg72; reg78=0.5*reg66; reg79=0.5*reg68;
    reg80=0.5*reg64; reg82=0.5*reg63; reg83=0.5*reg62; reg84=0.5*reg65; reg85=reg38*reg67;
    reg86=reg71+reg69; reg87=reg70*reg67; reg88=reg73+reg7; reg89=reg75+reg60; reg90=reg86+reg85;
    reg91=reg66*reg16; reg92=reg79*reg28; reg93=reg88+reg87; reg94=2*reg89; reg95=reg62*reg12;
    reg96=reg82*reg28; reg100=reg63*reg16; reg101=reg83*reg28; reg102=reg28*reg80; reg103=reg84*reg28;
    reg104=reg64*reg12; reg105=reg65*reg16; reg106=reg78*reg28; reg107=reg68*reg12; reg108=reg76*reg28;
    reg109=reg77*reg28; reg110=reg74*reg12; reg112=reg72*reg16; reg113=reg79*reg38; reg114=reg66*reg25;
    reg116=reg63*reg25; reg117=reg83*reg38; reg118=reg78*reg38; reg119=reg68*reg19; reg120=reg64*reg93;
    reg121=reg77*reg70; reg122=reg62*reg14; reg123=reg82*reg70; reg124=reg82*reg38; reg125=reg62*reg19;
    reg101=reg100+reg101; reg100=reg63*reg9; reg126=reg83*reg70; reg127=reg84*reg94; reg128=reg64*reg19;
    T reg129=reg84*reg38; T reg130=reg76*reg70; T reg131=reg72*reg9; T reg132=reg66*reg90; T reg133=reg79*reg94;
    reg108=reg112+reg108; reg112=reg78*reg94; T reg134=reg68*reg93; T reg135=reg94*reg80; T reg136=reg65*reg90;
    reg92=reg92-reg91; T reg137=reg76*reg38; T reg138=reg72*reg25; T reg139=reg84*reg70; T reg140=reg64*reg14;
    reg109=reg110+reg109; reg110=reg78*reg70; T reg141=reg74*reg14; T reg142=reg68*reg14; T reg143=reg23*var_inter[1];
    T reg144=reg11*var_inter[0]; reg103=reg103-reg104; T reg145=reg70*reg80; T reg146=reg65*reg9; reg107=reg107-reg106;
    reg105=reg105-reg102; reg96=reg95+reg96; reg95=reg66*reg9; T reg147=reg65*reg25; T reg148=reg38*reg80;
    T reg149=reg79*reg70; T reg150=reg82*reg94; T reg151=reg62*reg93; reg129=reg129-reg128; reg130=reg131+reg130;
    reg131=reg83*reg94; T reg152=reg63*reg90; T reg153=var_inter[0]*var_inter[1]; reg103=2*reg103; reg149=reg149-reg95;
    reg124=reg125+reg124; reg126=reg100+reg126; reg107=2*reg107; reg146=reg146-reg145; reg119=reg119-reg118;
    reg101=2*reg101; reg139=reg139-reg140; reg100=reg120-reg127; reg125=reg11*reg23; reg117=reg116+reg117;
    reg96=2*reg96; reg123=reg122+reg123; reg105=2*reg105; reg141=reg121+reg141; reg116=reg76*reg94;
    reg113=reg113-reg114; reg137=reg138+reg137; reg121=reg74*reg93; reg92=2*reg92; reg122=reg72*reg90;
    reg138=reg144*elem.f_vol_e[1]; reg142=reg142-reg110; T reg154=reg143*elem.f_vol_e[0]; T reg155=reg143*elem.f_vol_e[1]; T reg156=reg112-reg134;
    T reg157=reg135-reg136; reg147=reg147-reg148; T reg158=reg77*reg94; reg109=2*reg109; T reg159=reg144*elem.f_vol_e[0];
    reg108=2*reg108; T reg160=reg132-reg133; T reg161=reg103*reg82; T reg162=reg103*reg78; T reg163=reg72*reg117;
    T reg164=reg125*elem.f_vol_e[1]; T reg165=reg139*reg62; T reg166=reg139*reg68; T reg167=reg101*reg76; T reg168=reg72*reg124;
    T reg169=reg101*reg82; T reg170=reg125*elem.f_vol_e[0]; T reg171=reg66*reg113; T reg172=reg92*reg79; T reg173=reg149*reg74;
    T reg174=reg92*reg77; reg100=reg100-reg154; T reg175=reg96*reg79; T reg176=reg105*reg78; T reg177=reg146*reg68;
    T reg178=reg66*reg124; T reg179=reg142*reg74; reg157=reg157-reg155; T reg180=reg101*reg79; T reg181=reg66*reg117;
    T reg182=reg107*reg77; T reg183=reg103*reg76; T reg184=reg103*reg79; T reg185=reg66*reg129; T reg186=reg126*reg62;
    T reg187=reg105*reg79; T reg188=reg66*reg147; T reg189=reg96*reg82; T reg190=reg123*reg62; T reg191=reg108*reg77;
    reg160=reg160-reg138; reg156=reg156-reg159; T reg192=reg141*reg74; T reg193=reg116+reg122; T reg194=reg123*reg68;
    T reg195=reg158+reg121; T reg196=reg107*reg76; T reg197=reg96*reg78; T reg198=reg105*reg80; T reg199=reg65*reg147;
    T reg200=reg126*reg68; T reg201=reg105*reg76; T reg202=reg146*reg64; T reg203=reg105*reg84; T reg204=reg139*reg64;
    T reg205=reg103*reg84; T reg206=reg72*reg147; T reg207=reg101*reg78; T reg208=reg72*reg137; T reg209=reg108*reg76;
    T reg210=reg153*elem.f_vol_e[0]; T reg211=reg153*elem.f_vol_e[1]; T reg212=reg77*reg109; T reg213=reg142*reg68; T reg214=reg146*reg74;
    T reg215=reg105*reg77; T reg216=reg107*reg78; T reg217=reg103*reg77; reg139=reg139*reg74; T reg218=reg96*reg76;
    T reg219=reg105*reg82; reg146=reg146*reg62; T reg220=reg72*reg113; T reg221=reg92*reg76; T reg222=reg92*reg78;
    T reg223=reg149*reg68; T reg224=reg101*reg83; T reg225=reg63*reg117; T reg226=reg72*reg119; T reg227=reg103*reg83;
    T reg228=reg96*reg77; T reg229=reg123*reg74; T reg230=reg63*reg129; T reg231=reg72*reg129; T reg232=reg150+reg151;
    T reg233=reg131+reg152; T reg234=reg130*reg74; reg147=reg63*reg147; T reg235=reg126*reg74; T reg236=reg101*reg77;
    reg105=reg105*reg83; reg185=reg184-reg185; reg165=reg161+reg165; reg162=reg166-reg162; reg216=reg213-reg216;
    reg181=reg180-reg181; reg171=reg172-reg171; reg230=reg227+reg230; reg212=reg192+reg212; reg207=reg200-reg207;
    reg146=reg219+reg146; reg197=reg194-reg197; reg178=reg175-reg178; reg147=reg105+reg147; reg176=reg177-reg176;
    reg188=reg187-reg188; reg206=reg201+reg206; reg225=reg224+reg225; reg231=reg183+reg231; reg190=reg189+reg190;
    reg220=reg221+reg220; reg226=reg196+reg226; reg168=reg218+reg168; reg229=reg228+reg229; reg163=reg167+reg163;
    reg222=reg223-reg222; reg105=reg211+reg233; reg161=reg210+reg232; reg160=reg13*reg160; reg169=reg186+reg169;
    reg156=reg13*reg156; reg166=reg164+reg193; reg157=reg157*reg13; reg167=reg170+reg195; reg198=reg199-reg198;
    reg100=reg100*reg13; reg139=reg217+reg139; reg234=reg191+reg234; reg235=reg236+reg235; reg214=reg215+reg214;
    reg179=reg182+reg179; reg208=reg209+reg208; reg204=reg205-reg204; reg173=reg174+reg173; reg202=reg203-reg202;
    reg197=reg13*reg197; reg146=reg13*reg146; reg226=reg13*reg226; reg139=reg13*reg139; reg202=reg13*reg202;
    reg225=reg13*reg225; reg216=reg13*reg216; reg235=reg13*reg235; reg229=reg13*reg229; reg231=reg13*reg231;
    reg230=reg13*reg230; reg172=reg13*reg105; reg147=reg13*reg147; reg214=reg13*reg214; reg174=reg13*reg161;
    reg198=reg13*reg198; reg160=ponderation*reg160; reg204=reg13*reg204; reg206=reg13*reg206; reg156=ponderation*reg156;
    reg208=reg13*reg208; reg212=reg212*reg13; reg175=reg13*reg166; reg177=reg13*reg167; reg207=reg13*reg207;
    reg234=reg234*reg13; reg162=reg13*reg162; reg179=reg179*reg13; reg176=reg13*reg176; reg173=reg173*reg13;
    reg171=reg13*reg171; reg100=ponderation*reg100; reg178=reg13*reg178; reg157=ponderation*reg157; reg181=reg13*reg181;
    reg163=reg13*reg163; reg185=reg13*reg185; reg169=reg169*reg13; reg188=reg13*reg188; reg168=reg13*reg168;
    reg165=reg13*reg165; reg190=reg13*reg190; reg222=reg13*reg222; reg220=reg13*reg220; T tmp_3_3=ponderation*reg171;
    T tmp_2_4=ponderation*reg197; T tmp_6_7=ponderation*reg202; T tmp_2_2=ponderation*reg216; T tmp_0_3=ponderation*reg173; T tmp_6_6=ponderation*reg204;
    T tmp_1_3=ponderation*reg220; T tmp_2_7=ponderation*reg176; T tmp_1_1=ponderation*reg208; T tmp_2_5=ponderation*reg207; T tmp_0_2=ponderation*reg179;
    T tmp_0_6=ponderation*reg139; T tmp_1_7=ponderation*reg206; T tmp_0_5=ponderation*reg235; T tmp_4_7=ponderation*reg146; T tmp_0_7=ponderation*reg214;
    T tmp_2_6=ponderation*reg162; T tmp_4_6=ponderation*reg165; T tmp_0_1=ponderation*reg234; T tmp_1_2=ponderation*reg226; reg139=ponderation*reg172;
    T vec_5=reg139; T tmp_3_7=ponderation*reg188; T tmp_0_4=ponderation*reg229; T tmp_5_7=ponderation*reg147; T tmp_2_3=ponderation*reg222;
    reg146=ponderation*reg174; T vec_4=reg146; T tmp_4_5=ponderation*reg169; T tmp_5_6=ponderation*reg230; T vec_3=-reg160;
    T tmp_3_6=ponderation*reg185; T vec_2=-reg156; T tmp_1_5=ponderation*reg163; reg147=ponderation*reg175; T vec_1=reg147;
    T tmp_0_0=ponderation*reg212; T tmp_3_5=ponderation*reg181; T tmp_4_4=ponderation*reg190; T vec_7=-reg157; reg156=ponderation*reg177;
    T vec_0=reg156; T tmp_1_4=ponderation*reg168; T tmp_3_4=ponderation*reg178; T tmp_7_7=ponderation*reg198; T tmp_5_5=ponderation*reg225;
    T vec_6=-reg100; T tmp_1_6=ponderation*reg231;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_3;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_4;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_5;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_7;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_3;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_4;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_5;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_7;
    matrix(indices[1]+0,indices[1]+0) += tmp_2_2;
    matrix(indices[1]+0,indices[1]+1) += tmp_2_3;
    matrix(indices[1]+0,indices[2]+0) += tmp_2_4;
    matrix(indices[1]+0,indices[2]+1) += tmp_2_5;
    matrix(indices[1]+0,indices[3]+0) += tmp_2_6;
    matrix(indices[1]+0,indices[3]+1) += tmp_2_7;
    matrix(indices[1]+1,indices[1]+1) += tmp_3_3;
    matrix(indices[1]+1,indices[2]+0) += tmp_3_4;
    matrix(indices[1]+1,indices[2]+1) += tmp_3_5;
    matrix(indices[1]+1,indices[3]+0) += tmp_3_6;
    matrix(indices[1]+1,indices[3]+1) += tmp_3_7;
    matrix(indices[2]+0,indices[2]+0) += tmp_4_4;
    matrix(indices[2]+0,indices[2]+1) += tmp_4_5;
    matrix(indices[2]+0,indices[3]+0) += tmp_4_6;
    matrix(indices[2]+0,indices[3]+1) += tmp_4_7;
    matrix(indices[2]+1,indices[2]+1) += tmp_5_5;
    matrix(indices[2]+1,indices[3]+0) += tmp_5_6;
    matrix(indices[2]+1,indices[3]+1) += tmp_5_7;
    matrix(indices[3]+0,indices[3]+0) += tmp_6_6;
    matrix(indices[3]+0,indices[3]+1) += tmp_6_7;
    matrix(indices[3]+1,indices[3]+1) += tmp_7_7;
sollicitation[indices[0]+0] += vec_0;
sollicitation[indices[0]+1] += vec_1;
sollicitation[indices[1]+0] += vec_2;
sollicitation[indices[1]+1] += vec_3;
sollicitation[indices[2]+0] += vec_4;
sollicitation[indices[2]+1] += vec_5;
sollicitation[indices[3]+0] += vec_6;
sollicitation[indices[3]+1] += vec_7;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=pow((*f.m).v1[0],2); T reg1=pow((*f.m).v1[1],2); reg1=reg0+reg1; reg0=pow((*f.m).v1[2],2); T reg2=pow((*f.m).v2[0],2);
    T reg3=pow((*f.m).v2[1],2); reg0=reg1+reg0; reg3=reg2+reg3; reg1=pow((*f.m).v2[2],2); reg2=2*(*f.m).shear_modulus_23;
    T reg4=2*(*f.m).shear_modulus_13; reg0=pow(reg0,0.5); reg1=reg3+reg1; reg2=1.0/reg2; reg4=1.0/reg4;
    reg3=2*(*f.m).shear_modulus_12; reg1=pow(reg1,0.5); T reg5=reg4*reg2; T reg6=(*f.m).v1[1]/reg0; T reg7=(*f.m).v1[0]/reg0;
    reg3=1.0/reg3; T reg8=2*reg6; T reg9=reg3*reg5; T reg10=2*reg7; T reg11=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg12=1.0/(*f.m).elastic_modulus_3; T reg13=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg14=(*f.m).v2[1]/reg1; T reg15=(*f.m).v2[0]/reg1; reg0=(*f.m).v1[2]/reg0;
    T reg16=reg14*reg8; T reg17=reg15*reg10; T reg18=pow(reg15,2); reg1=(*f.m).v2[2]/reg1; T reg19=2*reg0;
    T reg20=1.0/(*f.m).elastic_modulus_1; T reg21=reg12*reg9; T reg22=reg11*reg9; T reg23=reg13*reg9; T reg24=1.0/(*f.m).elastic_modulus_2;
    T reg25=pow(reg14,2); T reg26=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; reg19=reg19*reg1; T reg27=reg24*reg21; T reg28=pow(reg1,2);
    T reg29=reg16*reg26; T reg30=reg17*reg20; T reg31=reg26*reg25; T reg32=reg20*reg18; T reg33=reg26*reg18;
    T reg34=reg24*reg25; T reg35=reg17*reg26; T reg36=reg16*reg24; T reg37=reg13*reg23; T reg38=reg26*reg21;
    T reg39=reg13*reg22; reg29=reg30-reg29; reg30=reg19*reg11; reg37=reg27-reg37; reg31=reg32-reg31;
    reg27=reg11*reg28; reg39=reg38+reg39; reg32=reg26*reg23; T reg40=reg13*reg28; T reg41=reg16*reg13;
    T reg42=reg24*reg22; reg33=reg34-reg33; reg34=reg11*reg18; T reg43=reg19*reg13; T reg44=reg13*reg25;
    reg35=reg36-reg35; reg36=reg17*reg11; T reg45=pow(reg7,2); T reg46=pow(reg6,2); reg30=reg29-reg30;
    reg29=reg26*reg45; T reg47=reg24*reg46; reg41=reg36+reg41; reg19=reg19*reg12; reg36=reg7*reg15;
    T reg48=reg6*reg14; reg44=reg34+reg44; reg34=reg12*reg28; T reg49=reg26*reg46; T reg50=reg20*reg37;
    T reg51=reg20*reg45; T reg52=reg26*reg39; T reg53=reg32+reg42; reg43=reg35-reg43; reg35=reg6*reg15;
    T reg54=reg7*reg14; reg40=reg33-reg40; reg27=reg31-reg27; reg31=pow(reg0,2); reg33=2*reg15;
    T reg55=reg13*reg46; T reg56=reg13*reg5; T reg57=reg30*reg45; reg41=reg19-reg41; reg19=reg0*reg1;
    T reg58=reg30*reg36; T reg59=reg43*reg48; T reg60=reg35+reg54; T reg61=reg43*reg46; T reg62=reg11*reg23;
    T reg63=reg11*reg31; reg49=reg51-reg49; reg51=reg24*reg9; T reg64=reg27*reg18; reg44=reg34-reg44;
    reg43=reg43*reg25; reg34=reg3*reg2; reg30=reg30*reg18; reg52=reg50-reg52; reg50=reg40*reg25;
    T reg65=reg11*reg53; T reg66=reg0*reg14; T reg67=reg6*reg1; T reg68=reg27*reg45; T reg69=reg15*reg14;
    T reg70=reg12*reg5; T reg71=reg13*reg31; T reg72=reg40*reg46; T reg73=reg7*reg1; reg27=reg27*reg36;
    reg9=reg26*reg9; reg29=reg47-reg29; reg47=reg11*reg22; reg40=reg40*reg48; T reg74=reg11*reg45;
    reg21=reg20*reg21; T reg75=reg0*reg15; reg5=reg11*reg5; T reg76=reg6*reg10; T reg77=reg3*reg4;
    T reg78=reg24*reg70; T reg79=reg41*reg31; reg22=reg26*reg22; T reg80=reg11*reg9; T reg81=reg11*reg51;
    T reg82=reg11*reg34; reg55=reg74+reg55; reg62=reg38+reg62; reg23=reg20*reg23; reg47=reg21-reg47;
    reg21=reg12*reg31; reg61=reg57+reg61; reg70=reg26*reg70; reg38=reg13*reg56; reg57=reg44*reg31;
    reg72=reg68+reg72; reg63=reg49-reg63; reg49=reg13*reg5; reg68=reg67-reg66; reg74=reg15*reg1;
    T reg83=reg3*reg60; T reg84=reg69*reg3; T reg85=reg12*reg34; reg34=reg13*reg34; reg40=reg27+reg40;
    reg27=reg33*reg14; reg71=reg29-reg71; reg29=reg44*reg19; reg65=reg52-reg65; reg52=reg73+reg75;
    reg59=reg58+reg59; reg58=reg41*reg19; reg41=reg41*reg28; reg43=reg30+reg43; reg44=reg44*reg28;
    reg50=reg64+reg50; reg30=reg0*reg10; reg73=reg75-reg73; reg64=reg11*reg77; reg5=reg24*reg5;
    reg75=reg83*reg60; reg58=reg59+reg58; reg59=reg84*reg60; reg29=reg40+reg29; reg40=reg84*reg76;
    reg57=reg72+reg57; reg38=reg78-reg38; reg49=reg70+reg49; reg70=reg71*reg46; reg72=reg63*reg45;
    reg78=reg4*reg52; T reg86=reg74*reg4; T reg87=reg24*reg85; reg85=reg26*reg85; T reg88=reg13*reg34;
    T reg89=reg13*reg82; reg12=reg12*reg77; T reg90=reg33*reg1; T reg91=2*reg14; T reg92=2*reg68;
    reg67=reg66+reg67; reg55=reg21-reg55; reg21=reg63*reg18; reg66=reg14*reg1; T reg93=reg71*reg25;
    T reg94=reg7*reg6; T reg95=reg83*reg76; reg79=reg61+reg79; reg51=reg20*reg51; reg37=reg37/reg65;
    reg62=reg62/reg65; reg81=reg32+reg81; reg39=reg39/reg65; reg47=reg47/reg65; reg56=reg26*reg56;
    reg80=reg23+reg80; reg77=reg13*reg77; reg83=reg83*reg27; reg22=reg23+reg22; reg44=reg50+reg44;
    reg41=reg43+reg41; reg9=reg26*reg9; reg84=reg84*reg27; reg23=reg78*reg30; reg40=reg57+reg40;
    reg34=reg26*reg34; reg43=reg78*reg90; reg63=reg63*reg36; reg71=reg71*reg48; reg95=reg79+reg95;
    reg50=reg47*reg18; reg57=reg62*reg45; reg61=reg39*reg27; reg79=reg37*reg76; T reg96=reg39*reg25;
    T reg97=reg37*reg46; T reg98=reg86*reg30; T reg99=reg86*reg90; reg84=reg44+reg84; reg78=reg78*reg52;
    reg75=reg58+reg75; reg44=pow(reg73,2); reg58=pow(reg68,2); T reg100=reg91*reg1; T reg101=reg47*reg27;
    T reg102=reg62*reg76; reg86=reg86*reg52; reg59=reg29+reg59; reg83=reg41+reg83; reg29=reg55*reg28;
    reg93=reg21+reg93; reg21=reg47*reg25; reg41=reg92*reg73; T reg103=reg62*reg46; T reg104=reg55*reg31;
    reg70=reg72+reg70; reg72=reg2*reg67; T reg105=reg66*reg2; T reg106=reg13*reg77; T reg107=reg7*reg0;
    reg49=reg26*reg49; T reg108=reg26*reg12; reg3=reg94*reg3; T reg109=reg56+reg5; reg12=reg24*reg12;
    reg88=reg87-reg88; reg89=reg85+reg89; reg81=reg81/reg65; reg53=reg53/reg65; reg85=reg37*reg45;
    reg87=reg39*reg18; reg80=reg80/reg65; reg22=reg22/reg65; reg9=reg51-reg9; reg51=reg0*reg8;
    T reg110=1-var_inter[0]; T reg111=1-var_inter[1]; reg82=reg24*reg82; reg13=reg13*reg64; reg38=reg20*reg38;
    T reg112=reg111*elem.pos(0)[1]; reg101=reg102+reg101; reg102=reg22*reg41; T reg113=reg81*reg45; T reg114=reg80*reg18;
    reg50=reg57+reg50; reg57=var_inter[0]*elem.pos(1)[0]; T reg115=reg53*reg44; reg96=reg97+reg96; reg97=reg72*reg100;
    reg43=reg83+reg43; reg89=reg26*reg89; reg83=reg105*reg100; reg49=reg38-reg49; reg106=reg12-reg106;
    reg13=reg108+reg13; reg12=reg80*reg25; reg38=reg81*reg76; reg108=reg80*reg27; T reg116=reg22*reg58;
    reg64=reg24*reg64; reg24=reg81*reg46; T reg117=reg53*reg41; reg21=reg103+reg21; reg103=reg22*reg44;
    reg61=reg79+reg61; reg109=reg11*reg109; reg79=reg111*elem.pos(1)[0]; T reg118=reg111*elem.pos(0)[0]; reg77=reg26*reg77;
    T reg119=reg34+reg82; reg88=reg20*reg88; T reg120=reg111*elem.pos(1)[1]; T reg121=reg72*reg51; reg55=reg55*reg19;
    reg71=reg63+reg71; reg23=reg95+reg23; reg63=reg110*elem.pos(0)[1]; reg95=var_inter[0]*elem.pos(1)[1]; T reg122=reg110*elem.pos(0)[0];
    T reg123=reg105*reg51; reg98=reg40+reg98; reg40=reg3*reg76; reg104=reg70+reg104; reg70=reg6*reg0;
    reg4=reg107*reg4; reg87=reg85+reg87; reg85=reg53*reg58; reg9=reg9/reg65; T reg124=reg3*reg27;
    reg105=reg105*reg67; reg86=reg59+reg86; reg78=reg75+reg78; reg72=reg72*reg67; reg99=reg84+reg99;
    reg29=reg93+reg29; reg105=reg86+reg105; reg79=reg79-reg118; reg120=reg120-reg112; reg59=var_inter[1]*elem.pos(2)[1];
    reg75=var_inter[0]*elem.pos(2)[0]; reg84=reg122+reg57; reg115=reg96+reg115; reg102=reg101+reg102; reg124=reg29+reg124;
    reg29=reg4*reg90; reg86=reg77+reg64; reg13=reg26*reg13; reg85=reg87+reg85; reg97=reg43+reg97;
    reg26=reg4*reg30; reg40=reg104+reg40; reg106=reg20*reg106; reg2=reg70*reg2; reg89=reg88-reg89;
    reg119=reg11*reg119; reg109=reg49-reg109; reg83=reg99+reg83; reg72=reg78+reg72; reg121=reg23+reg121;
    reg116=reg50+reg116; reg41=reg9*reg41; reg108=reg38+reg108; reg3=reg3*reg60; reg55=reg71+reg55;
    reg20=reg9*reg44; reg12=reg24+reg12; reg103=reg21+reg103; reg21=reg9*reg58; reg114=reg113+reg114;
    reg123=reg98+reg123; reg23=reg95+reg63; reg24=var_inter[1]*elem.pos(2)[0]; reg38=var_inter[0]*elem.pos(2)[1]; reg117=reg61+reg117;
    reg13=reg106-reg13; reg29=reg124+reg29; reg43=reg2*reg100; reg20=reg12+reg20; reg109=reg109/reg65;
    reg12=reg73*reg68; reg49=reg6*reg68; reg50=reg7*reg73; reg61=var_inter[1]*elem.pos(3)[0]; reg71=reg123*reg72;
    reg78=reg83*reg72; reg79=reg24+reg79; reg24=reg97*reg105; reg119=reg89-reg119; reg41=reg108+reg41;
    reg3=reg55+reg3; reg55=reg110*elem.pos(3)[1]; reg4=reg4*reg52; reg87=reg69*reg103; reg88=reg94*reg115;
    reg89=reg69*reg102; reg75=reg75-reg84; reg93=var_inter[1]*elem.pos(3)[1]; reg96=reg69*reg116; reg98=reg94*reg85;
    reg59=reg120+reg59; reg86=reg11*reg86; reg21=reg114+reg21; reg11=reg2*reg51; reg99=reg94*reg117;
    reg26=reg40+reg26; reg38=reg38-reg23; reg40=reg121*reg105; reg101=reg110*elem.pos(3)[0]; reg104=reg15*reg73;
    reg106=reg14*reg68; reg108=reg6*reg73; reg113=reg85*reg46; reg86=reg13-reg86; reg13=reg7*reg68;
    reg50=reg49+reg50; reg49=reg12*reg20; reg59=reg59-reg93; reg79=reg79-reg61; reg89=reg99+reg89;
    reg99=reg12*reg41; reg114=reg109*reg36; reg120=reg109*reg48; reg124=reg109*reg60; reg85=reg85*reg45;
    T reg125=reg116*reg18; T reg126=reg115*reg45; T reg127=reg103*reg18; T reg128=reg117*reg45; T reg129=reg102*reg18;
    reg43=reg29+reg43; reg2=reg2*reg67; reg4=reg3+reg4; reg55=reg38+reg55; reg101=reg75+reg101;
    reg11=reg26+reg11; reg24=reg78-reg24; reg116=reg116*reg25; reg119=reg119/reg65; reg115=reg115*reg46;
    reg103=reg103*reg25; reg117=reg117*reg46; reg102=reg102*reg25; reg40=reg71-reg40; reg3=reg123*reg97;
    reg26=reg121*reg83; reg87=reg88+reg87; reg96=reg98+reg96; reg29=reg12*reg21; reg38=reg60*reg114;
    reg116=reg113+reg116; reg71=reg44*reg21; reg29=reg96+reg29; reg103=reg115+reg103; reg75=reg58*reg41;
    reg129=reg128+reg129; reg78=reg15*reg68; reg88=reg14*reg73; reg96=reg44*reg20; reg49=reg87+reg49;
    reg20=reg58*reg20; reg127=reg126+reg127; reg87=reg60*reg120; reg104=reg106+reg104; reg21=reg58*reg21;
    reg125=reg85+reg125; reg102=reg117+reg102; reg41=reg44*reg41; reg85=reg119*reg50; reg108=reg119*reg108;
    reg13=reg119*reg13; reg26=reg3-reg26; reg99=reg89+reg99; reg3=reg60*reg124; reg2=reg4+reg2;
    reg4=reg101*reg59; reg89=reg43*reg40; reg98=reg24*reg11; reg65=reg86/reg65; reg86=reg55*reg79;
    reg3=reg99+reg3; reg99=reg65*reg104; reg106=reg50*reg85; reg113=reg16*reg120; reg96=reg103+reg96;
    reg88=reg65*reg88; reg78=reg65*reg78; reg20=reg127+reg20; reg103=reg16*reg114; reg71=reg116+reg71;
    reg38=reg29+reg38; reg120=reg17*reg120; reg29=reg11*reg105; reg115=reg50*reg13; reg8=reg73*reg8;
    reg116=reg123*reg2; reg123=reg123*reg43; reg117=reg83*reg11; reg10=reg68*reg10; reg75=reg129+reg75;
    reg87=reg49+reg87; reg41=reg102+reg41; reg49=reg16*reg124; reg105=reg43*reg105; reg102=reg2*reg26;
    reg83=reg83*reg2; reg89=reg98-reg89; reg98=reg50*reg108; reg21=reg125+reg21; reg4=reg86-reg4;
    reg114=reg17*reg114; reg124=reg17*reg124; reg120=reg20+reg120; reg20=reg10*reg108; reg86=reg10*reg13;
    reg101=reg101/reg4; reg116=reg29-reg116; reg59=reg59/reg4; reg29=reg97*reg11; reg125=reg121*reg2;
    reg114=reg21+reg114; reg79=reg79/reg4; reg98=reg87+reg98; reg21=reg104*reg88; reg87=reg104*reg99;
    reg106=reg3+reg106; reg55=reg55/reg4; reg3=reg10*reg85; reg49=reg41+reg49; reg108=reg8*reg108;
    reg113=reg96+reg113; reg85=reg8*reg85; reg13=reg8*reg13; reg91=reg91*reg73; reg83=reg105-reg83;
    reg103=reg71+reg103; reg124=reg75+reg124; reg102=reg89+reg102; reg33=reg33*reg68; reg115=reg38+reg115;
    reg11=reg11*reg72; reg123=reg117-reg123; reg38=reg104*reg78; reg72=reg43*reg72; reg2=reg97*reg2;
    reg43=reg121*reg43; reg83=reg83/reg102; reg41=reg91*reg99; reg85=reg49+reg85; reg49=var_inter[0]*reg79;
    reg2=reg72-reg2; reg71=reg91*reg88; reg108=reg113+reg108; reg72=var_inter[0]*reg59; reg75=reg91*reg78;
    reg13=reg103+reg13; reg89=reg110*reg59; reg3=reg124+reg3; reg99=reg33*reg99; reg96=reg110*reg79;
    reg97=reg111*reg101; reg20=reg120+reg20; reg38=reg115+reg38; reg88=reg33*reg88; reg103=1-(*f.m).resolution;
    reg123=reg123/reg102; reg43=reg29-reg43; reg29=var_inter[1]*reg55; reg116=reg116/reg102; reg125=reg11-reg125;
    reg11=reg111*reg55; reg105=var_inter[1]*reg101; reg78=reg33*reg78; reg87=reg106+reg87; reg21=reg98+reg21;
    reg86=reg114+reg86; reg40=reg40/reg102; reg41=reg85+reg41; reg85=reg97+reg49; reg98=reg38*reg103;
    reg75=reg13+reg75; reg13=reg97-reg96; reg106=reg11+reg72; reg113=reg29-reg72; reg71=reg108+reg71;
    reg108=reg89+reg29; reg114=reg21*reg103; reg115=reg87*reg103; reg117=reg96+reg105; reg120=reg49-reg105;
    reg78=reg86+reg78; reg86=reg89-reg11; reg2=reg2/reg102; reg24=reg24/reg102; reg125=reg125/reg102;
    reg26=reg26/reg102; reg102=reg43/reg102; reg43=(*f.m).resolution*reg123; reg121=(*f.m).resolution*reg116; reg88=reg20+reg88;
    reg99=reg3+reg99; reg3=(*f.m).resolution*reg83; reg20=0.5*reg106; reg124=0.5*reg85; reg126=0.5*reg120;
    reg3=reg98+reg3; reg121=reg114-reg121; reg43=reg115+reg43; reg98=(*f.m).resolution*reg102; reg114=0.5*reg108;
    reg115=(*f.m).resolution*reg40; reg127=(*f.m).resolution*reg24; reg128=(*f.m).resolution*reg26; reg129=reg103*reg78; T reg130=reg103*reg88;
    T reg131=reg103*reg99; T reg132=reg75*reg103; T reg133=reg71*reg103; T reg134=reg41*reg103; T reg135=0.5*reg13;
    T reg136=0.5*reg86; T reg137=(*f.m).resolution*reg2; T reg138=(*f.m).resolution*reg125; T reg139=0.5*reg117; T reg140=0.5*reg113;
    T reg141=reg120*reg121; T reg142=reg140*reg43; T reg143=reg86*reg3; reg98=reg134-reg98; reg134=reg135*reg43;
    reg138=reg133+reg138; reg137=reg132-reg137; reg131=reg128+reg131; reg115=reg130-reg115; reg129=reg127+reg129;
    reg127=reg108*reg3; reg128=reg139*reg43; reg130=reg126*reg43; reg132=reg113*reg3; reg133=reg117*reg121;
    T reg144=reg43*reg114; T reg145=reg106*reg3; T reg146=reg124*reg43; T reg147=reg85*reg121; T reg148=reg20*reg43;
    T reg149=reg13*reg121; T reg150=reg136*reg43; T reg151=reg108*reg137; reg128=reg128-reg127; reg134=reg143+reg134;
    reg143=reg86*reg129; T reg152=reg140*reg98; T reg153=reg120*reg138; T reg154=reg108*reg129; T reg155=reg126*reg98;
    T reg156=reg113*reg137; T reg157=reg139*reg131; reg130=reg132+reg130; reg132=reg20*reg98; T reg158=reg85*reg138;
    reg150=reg149+reg150; reg149=reg13*reg115; T reg159=reg124*reg98; T reg160=reg106*reg137; reg133=reg133-reg144;
    T reg161=reg117*reg115; T reg162=reg136*reg98; T reg163=reg136*reg131; T reg164=reg131*reg114; T reg165=reg13*reg138;
    T reg166=reg140*reg131; T reg167=reg113*reg129; reg145=reg145-reg146; T reg168=reg98*reg114; T reg169=reg126*reg131;
    T reg170=reg117*reg138; T reg171=reg120*reg115; T reg172=reg85*reg115; T reg173=reg135*reg98; T reg174=reg86*reg137;
    T reg175=reg106*reg129; T reg176=reg124*reg131; T reg177=reg20*reg131; T reg178=reg135*reg131; reg142=reg141+reg142;
    reg148=reg148-reg147; reg141=reg139*reg98; reg160=reg160-reg159; reg133=2*reg133; reg162=reg165+reg162;
    reg148=2*reg148; reg145=2*reg145; reg161=reg161-reg164; reg143=reg178+reg143; reg163=reg149+reg163;
    reg150=2*reg150; reg175=reg175-reg176; reg141=reg141-reg151; reg130=2*reg130; reg152=reg153+reg152;
    reg134=2*reg134; reg142=2*reg142; reg173=reg174+reg173; reg128=2*reg128; reg155=reg156+reg155;
    reg170=reg170-reg168; reg166=reg171+reg166; reg157=reg157-reg154; reg177=reg177-reg172; reg169=reg167+reg169;
    reg132=reg132-reg158; reg149=reg139*reg134; reg153=reg117*reg170; reg156=reg133*reg114; reg165=reg133*reg140;
    reg167=reg120*reg141; reg171=reg120*reg170; reg174=reg142*reg135; reg178=reg157*reg86; T reg179=reg161*reg86;
    T reg180=reg133*reg135; T reg181=reg128*reg135; T reg182=reg166*reg86; T reg183=reg120*reg162; T reg184=reg150*reg140;
    T reg185=reg128*reg20; T reg186=reg85*reg141; T reg187=reg133*reg20; T reg188=reg85*reg170; T reg189=reg173*reg120;
    T reg190=reg140*reg134; T reg191=reg126*reg134; T reg192=reg143*reg113; T reg193=reg150*reg126; T reg194=reg163*reg113;
    T reg195=reg161*reg113; T reg196=reg133*reg126; T reg197=reg145*reg126; T reg198=reg175*reg113; T reg199=reg145*reg114;
    T reg200=reg148*reg126; T reg201=reg177*reg113; T reg202=reg117*reg160; T reg203=reg130*reg126; T reg204=reg169*reg113;
    T reg205=reg157*reg113; T reg206=reg142*reg126; T reg207=reg128*reg126; T reg208=reg128*reg140; T reg209=reg157*reg106;
    T reg210=reg128*reg124; T reg211=reg120*reg152; T reg212=reg142*reg140; T reg213=reg161*reg106; T reg214=reg133*reg124;
    T reg215=reg20*reg134; T reg216=reg173*reg85; T reg217=reg120*reg155; T reg218=reg130*reg140; T reg219=reg150*reg20;
    T reg220=reg85*reg162; T reg221=reg120*reg132; T reg222=reg148*reg140; T reg223=reg117*reg132; T reg224=reg85*reg160;
    T reg225=reg148*reg20; T reg226=reg85*reg132; T reg227=reg120*reg160; T reg228=reg145*reg140; T reg229=reg130*reg20;
    T reg230=reg85*reg155; T reg231=reg142*reg20; T reg232=reg85*reg152; T reg233=reg145*reg135; T reg234=reg175*reg86;
    T reg235=reg13*reg152; T reg236=reg142*reg136; T reg237=reg173*reg13; T reg238=reg136*reg134; T reg239=reg130*reg135;
    T reg240=reg135*reg134; T reg241=reg130*reg139; T reg242=reg169*reg108; T reg243=reg169*reg86; T reg244=reg142*reg139;
    T reg245=reg166*reg108; T reg246=reg13*reg155; T reg247=reg130*reg136; T reg248=reg128*reg139; reg157=reg157*reg108;
    T reg249=reg177*reg108; T reg250=reg133*reg139; T reg251=reg150*reg124; T reg252=reg175*reg106; T reg253=reg145*reg124;
    T reg254=reg163*reg106; T reg255=reg124*reg134; T reg256=reg143*reg106; T reg257=reg177*reg106; reg170=reg13*reg170;
    T reg258=reg148*reg124; reg133=reg133*reg136; reg177=reg177*reg86; reg169=reg169*reg106; T reg259=reg130*reg124;
    T reg260=reg13*reg141; T reg261=reg166*reg106; T reg262=reg142*reg124; T reg263=reg128*reg136; T reg264=reg163*reg86;
    T reg265=reg148*reg135; reg166=reg166*reg113; reg142=reg142*reg114; reg152=reg117*reg152; T reg266=reg150*reg139;
    T reg267=reg150*reg136; reg130=reg130*reg114; reg155=reg117*reg155; T reg268=reg143*reg86; T reg269=reg150*reg135;
    reg163=reg163*reg108; reg143=reg143*reg108; T reg270=reg13*reg162; T reg271=reg148*reg114; T reg272=reg145*reg20;
    T reg273=reg145*reg139; reg145=reg145*reg136; reg161=reg161*reg108; T reg274=reg148*reg139; reg132=reg13*reg132;
    reg148=reg148*reg136; reg173=reg173*reg117; reg128=reg128*reg114; reg134=reg134*reg114; reg141=reg117*reg141;
    reg160=reg13*reg160; reg175=reg175*reg108; reg162=reg117*reg162; reg150=reg150*reg114; reg170=reg133+reg170;
    reg204=reg203+reg204; reg205=reg207+reg205; reg255=reg256-reg255; reg206=reg166+reg206; reg227=reg228+reg227;
    reg230=reg229-reg230; reg183=reg184+reg183; reg226=reg225-reg226; reg240=reg268+reg240; reg260=reg263+reg260;
    reg201=reg200+reg201; reg232=reg231-reg232; reg235=reg236+reg235; reg186=reg185-reg186; reg189=reg190+reg189;
    reg270=reg267+reg270; reg199=reg202-reg199; reg198=reg197+reg198; reg188=reg187-reg188; reg246=reg247+reg246;
    reg160=reg145+reg160; reg194=reg193+reg194; reg132=reg148+reg132; reg192=reg191+reg192; reg195=reg196+reg195;
    reg234=reg233+reg234; reg237=reg238+reg237; reg249=reg274-reg249; reg242=reg241-reg242; reg245=reg244-reg245;
    reg243=reg239+reg243; reg157=reg248-reg157; reg161=reg250-reg161; reg182=reg174+reg182; reg167=reg208+reg167;
    reg134=reg173-reg134; reg175=reg273-reg175; reg150=reg162-reg150; reg271=reg223-reg271; reg163=reg266-reg163;
    reg178=reg181+reg178; reg130=reg155-reg130; reg142=reg152-reg142; reg143=reg149-reg143; reg179=reg180+reg179;
    reg171=reg165+reg171; reg156=reg153-reg156; reg128=reg141-reg128; reg253=reg252-reg253; reg177=reg265+reg177;
    reg216=reg215-reg216; reg220=reg219-reg220; reg258=reg257-reg258; reg217=reg218+reg217; reg214=reg213-reg214;
    reg221=reg222+reg221; reg259=reg169-reg259; reg251=reg254-reg251; reg210=reg209-reg210; reg211=reg212+reg211;
    reg224=reg272-reg224; reg262=reg261-reg262; reg264=reg269+reg264; reg227=reg4*reg227; reg221=reg4*reg221;
    reg205=reg4*reg205; reg171=reg4*reg171; reg206=reg206*reg4; reg143=reg4*reg143; reg167=reg4*reg167;
    reg177=reg177*reg4; reg163=reg4*reg163; reg183=reg4*reg183; reg217=reg4*reg217; reg175=reg4*reg175;
    reg195=reg4*reg195; reg243=reg4*reg243; reg211=reg4*reg211; reg189=reg4*reg189; reg170=reg4*reg170;
    reg220=reg4*reg220; reg255=reg4*reg255; reg216=reg4*reg216; reg214=reg4*reg214; reg251=reg4*reg251;
    reg210=reg4*reg210; reg262=reg4*reg262; reg253=reg4*reg253; reg182=reg4*reg182; reg258=reg4*reg258;
    reg178=reg4*reg178; reg259=reg4*reg259; reg156=reg4*reg156; reg264=reg264*reg4; reg128=reg4*reg128;
    reg234=reg234*reg4; reg142=reg4*reg142; reg130=reg4*reg130; reg179=reg4*reg179; reg271=reg4*reg271;
    reg249=reg4*reg249; reg150=reg4*reg150; reg242=reg4*reg242; reg134=reg4*reg134; reg245=reg4*reg245;
    reg161=reg4*reg161; reg157=reg4*reg157; reg192=reg4*reg192; reg160=reg4*reg160; reg270=reg4*reg270;
    reg188=reg4*reg188; reg194=reg4*reg194; reg186=reg4*reg186; reg198=reg4*reg198; reg132=reg4*reg132;
    reg232=reg4*reg232; reg237=reg4*reg237; reg246=reg4*reg246; reg230=reg4*reg230; reg201=reg4*reg201;
    reg235=reg4*reg235; reg226=reg4*reg226; reg240=reg240*reg4; reg204=reg4*reg204; reg224=reg4*reg224;
    reg260=reg4*reg260; reg199=reg4*reg199; T tmp_0_0=ponderation*reg240; T tmp_6_4=ponderation*reg242; T tmp_6_3=ponderation*reg249;
    T tmp_1_0=ponderation*reg237; T tmp_6_5=ponderation*reg245; T tmp_4_5=ponderation*reg206; T tmp_0_7=ponderation*reg179; T tmp_0_2=ponderation*reg234;
    T tmp_1_1=ponderation*reg270; T tmp_0_1=ponderation*reg264; T tmp_1_2=ponderation*reg160; T tmp_2_4=ponderation*reg259; T tmp_2_3=ponderation*reg258;
    T tmp_1_3=ponderation*reg132; T tmp_2_2=ponderation*reg253; T tmp_0_3=ponderation*reg177; T tmp_1_4=ponderation*reg246; T tmp_2_1=ponderation*reg251;
    T tmp_2_0=ponderation*reg255; T tmp_1_5=ponderation*reg235; T tmp_1_7=ponderation*reg170; T tmp_1_6=ponderation*reg260; T tmp_2_7=ponderation*reg214;
    T tmp_3_0=ponderation*reg216; T tmp_5_3=ponderation*reg221; T tmp_3_1=ponderation*reg220; T tmp_7_2=ponderation*reg199; T tmp_5_2=ponderation*reg227;
    T tmp_3_2=ponderation*reg224; T tmp_3_3=ponderation*reg226; T tmp_5_1=ponderation*reg183; T tmp_3_4=ponderation*reg230; T tmp_3_5=ponderation*reg232;
    T tmp_5_0=ponderation*reg189; T tmp_3_6=ponderation*reg186; T tmp_3_7=ponderation*reg188; T tmp_4_7=ponderation*reg195; T tmp_4_0=ponderation*reg192;
    T tmp_4_1=ponderation*reg194; T tmp_4_2=ponderation*reg198; T tmp_4_6=ponderation*reg205; T tmp_4_3=ponderation*reg201; T tmp_4_4=ponderation*reg204;
    T tmp_6_6=ponderation*reg157; T tmp_6_2=ponderation*reg175; T tmp_6_7=ponderation*reg161; T tmp_7_0=ponderation*reg134; T tmp_6_1=ponderation*reg163;
    T tmp_7_1=ponderation*reg150; T tmp_7_3=ponderation*reg271; T tmp_6_0=ponderation*reg143; T tmp_7_4=ponderation*reg130; T tmp_7_5=ponderation*reg142;
    T tmp_5_7=ponderation*reg171; T tmp_7_6=ponderation*reg128; T tmp_7_7=ponderation*reg156; T tmp_0_6=ponderation*reg178; T tmp_5_6=ponderation*reg167;
    T tmp_0_5=ponderation*reg182; T tmp_0_4=ponderation*reg243; T tmp_5_5=ponderation*reg211; T tmp_2_5=ponderation*reg262; T tmp_2_6=ponderation*reg210;
    T tmp_5_4=ponderation*reg217;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_3;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_4;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_5;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_7;
    matrix(indices[0]+1,indices[0]+0) += tmp_1_0;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_3;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_4;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_5;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_7;
    matrix(indices[1]+0,indices[0]+0) += tmp_2_0;
    matrix(indices[1]+0,indices[0]+1) += tmp_2_1;
    matrix(indices[1]+0,indices[1]+0) += tmp_2_2;
    matrix(indices[1]+0,indices[1]+1) += tmp_2_3;
    matrix(indices[1]+0,indices[2]+0) += tmp_2_4;
    matrix(indices[1]+0,indices[2]+1) += tmp_2_5;
    matrix(indices[1]+0,indices[3]+0) += tmp_2_6;
    matrix(indices[1]+0,indices[3]+1) += tmp_2_7;
    matrix(indices[1]+1,indices[0]+0) += tmp_3_0;
    matrix(indices[1]+1,indices[0]+1) += tmp_3_1;
    matrix(indices[1]+1,indices[1]+0) += tmp_3_2;
    matrix(indices[1]+1,indices[1]+1) += tmp_3_3;
    matrix(indices[1]+1,indices[2]+0) += tmp_3_4;
    matrix(indices[1]+1,indices[2]+1) += tmp_3_5;
    matrix(indices[1]+1,indices[3]+0) += tmp_3_6;
    matrix(indices[1]+1,indices[3]+1) += tmp_3_7;
    matrix(indices[2]+0,indices[0]+0) += tmp_4_0;
    matrix(indices[2]+0,indices[0]+1) += tmp_4_1;
    matrix(indices[2]+0,indices[1]+0) += tmp_4_2;
    matrix(indices[2]+0,indices[1]+1) += tmp_4_3;
    matrix(indices[2]+0,indices[2]+0) += tmp_4_4;
    matrix(indices[2]+0,indices[2]+1) += tmp_4_5;
    matrix(indices[2]+0,indices[3]+0) += tmp_4_6;
    matrix(indices[2]+0,indices[3]+1) += tmp_4_7;
    matrix(indices[2]+1,indices[0]+0) += tmp_5_0;
    matrix(indices[2]+1,indices[0]+1) += tmp_5_1;
    matrix(indices[2]+1,indices[1]+0) += tmp_5_2;
    matrix(indices[2]+1,indices[1]+1) += tmp_5_3;
    matrix(indices[2]+1,indices[2]+0) += tmp_5_4;
    matrix(indices[2]+1,indices[2]+1) += tmp_5_5;
    matrix(indices[2]+1,indices[3]+0) += tmp_5_6;
    matrix(indices[2]+1,indices[3]+1) += tmp_5_7;
    matrix(indices[3]+0,indices[0]+0) += tmp_6_0;
    matrix(indices[3]+0,indices[0]+1) += tmp_6_1;
    matrix(indices[3]+0,indices[1]+0) += tmp_6_2;
    matrix(indices[3]+0,indices[1]+1) += tmp_6_3;
    matrix(indices[3]+0,indices[2]+0) += tmp_6_4;
    matrix(indices[3]+0,indices[2]+1) += tmp_6_5;
    matrix(indices[3]+0,indices[3]+0) += tmp_6_6;
    matrix(indices[3]+0,indices[3]+1) += tmp_6_7;
    matrix(indices[3]+1,indices[0]+0) += tmp_7_0;
    matrix(indices[3]+1,indices[0]+1) += tmp_7_1;
    matrix(indices[3]+1,indices[1]+0) += tmp_7_2;
    matrix(indices[3]+1,indices[1]+1) += tmp_7_3;
    matrix(indices[3]+1,indices[2]+0) += tmp_7_4;
    matrix(indices[3]+1,indices[2]+1) += tmp_7_5;
    matrix(indices[3]+1,indices[3]+0) += tmp_7_6;
    matrix(indices[3]+1,indices[3]+1) += tmp_7_7;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=pow((*f.m).v1[0],2); T reg1=pow((*f.m).v1[1],2); reg1=reg0+reg1; reg0=pow((*f.m).v1[2],2); T reg2=pow((*f.m).v2[0],2);
    T reg3=pow((*f.m).v2[1],2); reg0=reg1+reg0; reg3=reg2+reg3; reg1=pow((*f.m).v2[2],2); reg2=2*(*f.m).shear_modulus_23;
    T reg4=2*(*f.m).shear_modulus_13; reg0=pow(reg0,0.5); reg1=reg3+reg1; reg2=1.0/reg2; reg4=1.0/reg4;
    reg3=2*(*f.m).shear_modulus_12; reg1=pow(reg1,0.5); T reg5=reg4*reg2; T reg6=(*f.m).v1[1]/reg0; T reg7=(*f.m).v1[0]/reg0;
    reg3=1.0/reg3; T reg8=2*reg6; T reg9=reg3*reg5; T reg10=2*reg7; T reg11=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg12=1.0/(*f.m).elastic_modulus_3; T reg13=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg14=(*f.m).v2[1]/reg1; T reg15=(*f.m).v2[0]/reg1; reg0=(*f.m).v1[2]/reg0;
    T reg16=reg14*reg8; T reg17=reg15*reg10; T reg18=pow(reg15,2); reg1=(*f.m).v2[2]/reg1; T reg19=2*reg0;
    T reg20=1.0/(*f.m).elastic_modulus_1; T reg21=reg12*reg9; T reg22=reg11*reg9; T reg23=reg13*reg9; T reg24=1.0/(*f.m).elastic_modulus_2;
    T reg25=pow(reg14,2); T reg26=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; reg19=reg19*reg1; T reg27=reg24*reg21; T reg28=pow(reg1,2);
    T reg29=reg16*reg26; T reg30=reg17*reg20; T reg31=reg26*reg25; T reg32=reg20*reg18; T reg33=reg26*reg18;
    T reg34=reg24*reg25; T reg35=reg17*reg26; T reg36=reg16*reg24; T reg37=reg13*reg23; T reg38=reg26*reg21;
    T reg39=reg13*reg22; reg29=reg30-reg29; reg30=reg19*reg11; reg37=reg27-reg37; reg31=reg32-reg31;
    reg27=reg11*reg28; reg39=reg38+reg39; reg32=reg26*reg23; T reg40=reg13*reg28; T reg41=reg16*reg13;
    T reg42=reg24*reg22; reg33=reg34-reg33; reg34=reg11*reg18; T reg43=reg19*reg13; T reg44=reg13*reg25;
    reg35=reg36-reg35; reg36=reg17*reg11; T reg45=pow(reg7,2); T reg46=pow(reg6,2); reg30=reg29-reg30;
    reg29=reg26*reg45; T reg47=reg24*reg46; reg41=reg36+reg41; reg19=reg19*reg12; reg36=reg7*reg15;
    T reg48=reg6*reg14; reg44=reg34+reg44; reg34=reg12*reg28; T reg49=reg26*reg46; T reg50=reg20*reg37;
    T reg51=reg20*reg45; T reg52=reg26*reg39; T reg53=reg32+reg42; reg43=reg35-reg43; reg35=reg6*reg15;
    T reg54=reg7*reg14; reg40=reg33-reg40; reg27=reg31-reg27; reg31=pow(reg0,2); reg33=2*reg15;
    T reg55=reg13*reg46; T reg56=reg13*reg5; T reg57=reg30*reg45; reg41=reg19-reg41; reg19=reg0*reg1;
    T reg58=reg30*reg36; T reg59=reg43*reg48; T reg60=reg35+reg54; T reg61=reg43*reg46; T reg62=reg11*reg23;
    T reg63=reg11*reg31; reg49=reg51-reg49; reg51=reg24*reg9; T reg64=reg27*reg18; reg44=reg34-reg44;
    reg43=reg43*reg25; reg34=reg3*reg2; reg30=reg30*reg18; reg52=reg50-reg52; reg50=reg40*reg25;
    T reg65=reg11*reg53; T reg66=reg0*reg14; T reg67=reg6*reg1; T reg68=reg27*reg45; T reg69=reg15*reg14;
    T reg70=reg12*reg5; T reg71=reg13*reg31; T reg72=reg40*reg46; T reg73=reg7*reg1; reg27=reg27*reg36;
    reg9=reg26*reg9; reg29=reg47-reg29; reg47=reg11*reg22; reg40=reg40*reg48; T reg74=reg11*reg45;
    reg21=reg20*reg21; T reg75=reg0*reg15; reg5=reg11*reg5; T reg76=reg6*reg10; T reg77=reg3*reg4;
    T reg78=reg24*reg70; T reg79=reg41*reg31; reg22=reg26*reg22; T reg80=reg11*reg9; T reg81=reg11*reg51;
    T reg82=reg11*reg34; reg55=reg74+reg55; reg62=reg38+reg62; reg23=reg20*reg23; reg47=reg21-reg47;
    reg21=reg12*reg31; reg61=reg57+reg61; reg70=reg26*reg70; reg38=reg13*reg56; reg57=reg44*reg31;
    reg72=reg68+reg72; reg63=reg49-reg63; reg49=reg13*reg5; reg68=reg67-reg66; reg74=reg15*reg1;
    T reg83=reg3*reg60; T reg84=reg69*reg3; T reg85=reg12*reg34; reg34=reg13*reg34; reg40=reg27+reg40;
    reg27=reg33*reg14; reg71=reg29-reg71; reg29=reg44*reg19; reg65=reg52-reg65; reg52=reg73+reg75;
    reg59=reg58+reg59; reg58=reg41*reg19; reg41=reg41*reg28; reg43=reg30+reg43; reg44=reg44*reg28;
    reg50=reg64+reg50; reg30=reg0*reg10; reg73=reg75-reg73; reg64=reg11*reg77; reg5=reg24*reg5;
    reg75=reg83*reg60; reg58=reg59+reg58; reg59=reg84*reg60; reg29=reg40+reg29; reg40=reg84*reg76;
    reg57=reg72+reg57; reg38=reg78-reg38; reg49=reg70+reg49; reg70=reg71*reg46; reg72=reg63*reg45;
    reg78=reg4*reg52; T reg86=reg74*reg4; T reg87=reg24*reg85; reg85=reg26*reg85; T reg88=reg13*reg34;
    T reg89=reg13*reg82; reg12=reg12*reg77; T reg90=reg33*reg1; T reg91=2*reg14; T reg92=2*reg68;
    reg67=reg66+reg67; reg55=reg21-reg55; reg21=reg7*reg6; reg66=reg14*reg1; reg79=reg61+reg79;
    reg61=reg83*reg76; reg83=reg83*reg27; reg51=reg20*reg51; reg37=reg37/reg65; reg62=reg62/reg65;
    reg81=reg32+reg81; reg39=reg39/reg65; T reg93=reg63*reg18; reg47=reg47/reg65; T reg94=reg71*reg25;
    reg80=reg23+reg80; reg77=reg13*reg77; reg44=reg50+reg44; reg22=reg23+reg22; reg84=reg84*reg27;
    reg9=reg26*reg9; reg56=reg26*reg56; reg41=reg43+reg41; reg34=reg26*reg34; reg82=reg24*reg82;
    reg38=reg20*reg38; reg49=reg26*reg49; reg23=reg56+reg5; reg88=reg87-reg88; reg89=reg85+reg89;
    reg43=reg24*reg12; reg12=reg26*reg12; reg63=reg63*reg36; reg71=reg71*reg48; reg50=1-var_inter[0];
    reg85=reg0*reg8; reg87=reg78*reg52; reg75=reg58+reg75; reg59=reg29+reg59; reg29=reg86*reg52;
    reg58=reg78*reg90; reg83=reg41+reg83; reg41=reg86*reg90; reg84=reg44+reg84; reg44=reg55*reg28;
    reg94=reg93+reg94; reg78=reg78*reg30; reg61=reg79+reg61; reg79=reg92*reg73; reg93=pow(reg73,2);
    T reg95=pow(reg68,2); T reg96=reg91*reg1; T reg97=reg47*reg27; T reg98=reg62*reg76; T reg99=reg47*reg25;
    T reg100=reg62*reg46; T reg101=reg47*reg18; T reg102=reg62*reg45; T reg103=reg39*reg27; T reg104=reg37*reg76;
    T reg105=reg39*reg25; T reg106=reg37*reg46; T reg107=reg13*reg64; reg13=reg13*reg77; T reg108=reg37*reg45;
    T reg109=reg7*reg0; T reg110=reg39*reg18; reg53=reg53/reg65; T reg111=reg66*reg2; reg80=reg80/reg65;
    reg86=reg86*reg30; T reg112=1-var_inter[1]; reg40=reg57+reg40; reg9=reg51-reg9; reg51=reg2*reg67;
    reg81=reg81/reg65; reg3=reg21*reg3; reg70=reg72+reg70; reg57=reg55*reg31; reg22=reg22/reg65;
    reg4=reg109*reg4; reg57=reg70+reg57; reg97=reg98+reg97; reg70=reg22*reg79; reg72=reg112*elem.pos(1)[0];
    reg98=reg81*reg45; T reg113=reg80*reg18; T reg114=reg51*reg67; reg87=reg75+reg87; reg75=reg111*reg85;
    reg86=reg40+reg86; reg40=reg22*reg93; reg99=reg100+reg99; reg100=reg6*reg0; T reg115=reg22*reg95;
    reg101=reg102+reg101; reg88=reg20*reg88; reg102=reg53*reg79; reg103=reg104+reg103; reg104=reg53*reg93;
    reg105=reg106+reg105; reg89=reg26*reg89; reg107=reg12+reg107; reg12=var_inter[0]*elem.pos(1)[0]; reg106=reg80*reg27;
    T reg116=reg81*reg76; T reg117=reg3*reg76; T reg118=reg34+reg82; T reg119=reg80*reg25; T reg120=reg81*reg46;
    reg77=reg26*reg77; T reg121=reg112*elem.pos(0)[1]; reg58=reg83+reg58; reg83=reg50*elem.pos(0)[1]; T reg122=var_inter[0]*elem.pos(1)[1];
    reg71=reg63+reg71; reg55=reg55*reg19; reg63=reg112*elem.pos(1)[1]; T reg123=reg111*reg96; reg110=reg108+reg110;
    reg108=reg50*elem.pos(0)[0]; T reg124=reg53*reg95; T reg125=reg112*elem.pos(0)[0]; reg78=reg61+reg78; reg61=reg51*reg85;
    reg13=reg43-reg13; reg111=reg111*reg67; reg29=reg59+reg29; reg44=reg94+reg44; reg43=reg3*reg27;
    reg23=reg11*reg23; reg49=reg38-reg49; reg9=reg9/reg65; reg64=reg24*reg64; reg51=reg51*reg96;
    reg41=reg84+reg41; reg118=reg11*reg118; reg24=reg4*reg30; reg89=reg88-reg89; reg75=reg86+reg75;
    reg23=reg49-reg23; reg38=var_inter[0]*elem.pos(2)[0]; reg79=reg9*reg79; reg106=reg116+reg106; reg49=reg9*reg93;
    reg119=reg120+reg119; reg59=reg9*reg95; reg113=reg98+reg113; reg123=reg41+reg123; reg51=reg58+reg51;
    reg41=reg4*reg90; reg43=reg44+reg43; reg61=reg78+reg61; reg124=reg110+reg124; reg70=reg97+reg70;
    reg40=reg99+reg40; reg115=reg101+reg115; reg2=reg100*reg2; reg102=reg103+reg102; reg104=reg105+reg104;
    reg44=reg77+reg64; reg107=reg26*reg107; reg13=reg20*reg13; reg117=reg57+reg117; reg20=var_inter[1]*elem.pos(2)[0];
    reg111=reg29+reg111; reg3=reg3*reg60; reg26=var_inter[1]*elem.pos(2)[1]; reg55=reg71+reg55; reg29=reg108+reg12;
    reg63=reg63-reg121; reg72=reg72-reg125; reg114=reg87+reg114; reg57=var_inter[0]*elem.pos(2)[1]; reg58=reg122+reg83;
    reg71=reg69*reg115; reg72=reg20+reg72; reg20=reg73*reg68; reg78=reg21*reg124; reg84=reg21*reg104;
    reg86=var_inter[1]*elem.pos(3)[0]; reg23=reg23/reg65; reg87=reg2*reg85; reg4=reg4*reg52; reg38=reg38-reg29;
    reg24=reg117+reg24; reg44=reg11*reg44; reg11=reg50*elem.pos(3)[0]; reg118=reg89-reg118; reg41=reg43+reg41;
    reg107=reg13-reg107; reg13=reg2*reg96; reg43=reg75*reg114; reg88=var_inter[1]*elem.pos(3)[1]; reg79=reg106+reg79;
    reg57=reg57-reg58; reg89=reg50*elem.pos(3)[1]; reg49=reg119+reg49; reg59=reg113+reg59; reg94=reg123*reg114;
    reg97=reg51*reg111; reg26=reg63+reg26; reg3=reg55+reg3; reg55=reg6*reg68; reg63=reg61*reg111;
    reg98=reg69*reg70; reg99=reg21*reg102; reg101=reg7*reg73; reg103=reg69*reg40; reg97=reg94-reg97;
    reg94=reg70*reg18; reg105=reg7*reg68; reg106=reg15*reg73; reg110=reg6*reg73; reg113=reg102*reg45;
    reg116=reg14*reg68; reg11=reg38+reg11; reg101=reg55+reg101; reg38=reg40*reg18; reg55=reg23*reg36;
    reg89=reg57+reg89; reg57=reg23*reg48; reg117=reg23*reg60; reg119=reg124*reg45; reg120=reg115*reg18;
    T reg126=reg20*reg79; reg98=reg99+reg98; reg99=reg20*reg49; reg103=reg84+reg103; reg84=reg20*reg59;
    reg71=reg78+reg71; reg78=reg61*reg123; T reg127=reg75*reg51; reg63=reg43-reg63; reg4=reg3+reg4;
    reg70=reg70*reg25; reg102=reg102*reg46; reg40=reg40*reg25; reg3=reg104*reg46; reg115=reg115*reg25;
    reg124=reg124*reg46; reg2=reg2*reg67; reg13=reg41+reg13; reg104=reg104*reg45; reg26=reg26-reg88;
    reg72=reg72-reg86; reg44=reg107-reg44; reg118=reg118/reg65; reg87=reg24+reg87; reg115=reg124+reg115;
    reg120=reg119+reg120; reg24=reg95*reg59; reg2=reg4+reg2; reg65=reg44/reg65; reg38=reg104+reg38;
    reg4=reg95*reg49; reg94=reg113+reg94; reg41=reg95*reg79; reg106=reg116+reg106; reg43=reg14*reg73;
    reg44=reg15*reg68; reg59=reg93*reg59; reg40=reg3+reg40; reg49=reg93*reg49; reg70=reg102+reg70;
    reg79=reg93*reg79; reg78=reg127-reg78; reg84=reg71+reg84; reg3=reg60*reg55; reg99=reg103+reg99;
    reg71=reg60*reg57; reg102=reg13*reg63; reg126=reg98+reg126; reg98=reg60*reg117; reg103=reg118*reg101;
    reg110=reg118*reg110; reg104=reg11*reg26; reg107=reg97*reg87; reg105=reg118*reg105; reg113=reg89*reg72;
    reg104=reg113-reg104; reg113=reg16*reg55; reg59=reg115+reg59; reg115=reg75*reg2; reg24=reg120+reg24;
    reg55=reg17*reg55; reg4=reg38+reg4; reg38=reg123*reg87; reg116=reg17*reg57; reg75=reg75*reg13;
    reg41=reg94+reg41; reg94=reg17*reg117; reg49=reg40+reg49; reg57=reg16*reg57; reg40=reg87*reg111;
    reg79=reg70+reg79; reg117=reg16*reg117; reg102=reg107-reg102; reg70=reg2*reg78; reg111=reg13*reg111;
    reg123=reg123*reg2; reg3=reg84+reg3; reg84=reg101*reg105; reg8=reg73*reg8; reg71=reg99+reg71;
    reg99=reg101*reg110; reg107=reg101*reg103; reg98=reg126+reg98; reg119=reg65*reg106; reg44=reg65*reg44;
    reg43=reg65*reg43; reg10=reg68*reg10; reg120=reg10*reg110; reg116=reg4+reg116; reg4=reg61*reg13;
    reg124=reg87*reg114; reg126=reg51*reg2; reg87=reg51*reg87; reg114=reg13*reg114; reg115=reg40-reg115;
    reg70=reg102+reg70; reg89=reg89/reg104; reg13=reg10*reg105; reg55=reg24+reg55; reg33=reg33*reg68;
    reg26=reg26/reg104; reg24=reg8*reg103; reg117=reg79+reg117; reg11=reg11/reg104; reg113=reg59+reg113;
    reg105=reg8*reg105; reg2=reg61*reg2; reg72=reg72/reg104; reg57=reg49+reg57; reg110=reg8*reg110;
    reg40=reg106*reg119; reg107=reg98+reg107; reg49=reg106*reg44; reg84=reg3+reg84; reg103=reg10*reg103;
    reg99=reg71+reg99; reg3=reg106*reg43; reg91=reg91*reg73; reg75=reg38-reg75; reg123=reg111-reg123;
    reg94=reg41+reg94; reg38=reg50*reg26; reg41=reg91*reg119; reg24=reg117+reg24; reg51=reg112*reg11;
    reg59=reg50*reg72; reg110=reg57+reg110; reg115=reg115/reg70; reg103=reg94+reg103; reg57=reg91*reg43;
    reg61=var_inter[1]*reg11; reg3=reg99+reg3; reg105=reg113+reg105; reg71=reg91*reg44; reg2=reg124-reg2;
    reg79=var_inter[0]*reg26; reg75=reg75/reg70; reg120=reg116+reg120; reg123=reg123/reg70; reg4=reg87-reg4;
    reg119=reg33*reg119; reg126=reg114-reg126; reg43=reg33*reg43; reg87=1-(*f.m).resolution; reg40=reg107+reg40;
    reg49=reg84+reg49; reg84=var_inter[0]*reg72; reg13=reg55+reg13; reg44=reg33*reg44; reg55=reg112*reg89;
    reg94=var_inter[1]*reg89; reg126=reg126/reg70; reg98=reg51-reg59; reg63=reg63/reg70; reg99=reg84-reg61;
    reg102=reg51+reg84; reg57=reg110+reg57; reg107=reg94-reg79; reg97=reg97/reg70; reg41=reg24+reg41;
    reg24=reg38-reg55; reg110=reg3*reg87; reg111=reg49*reg87; reg113=reg40*reg87; reg4=reg4/reg70;
    reg114=(*f.m).resolution*reg123; reg116=(*f.m).resolution*reg115; reg44=reg13+reg44; reg13=(*f.m).resolution*reg75; reg119=reg103+reg119;
    reg78=reg78/reg70; reg103=reg55+reg79; reg117=reg38+reg94; reg124=reg59+reg61; reg70=reg2/reg70;
    reg43=reg120+reg43; reg71=reg105+reg71; reg2=reg71*reg87; reg105=0.5*reg117; reg120=reg41*reg87;
    reg127=(*f.m).resolution*reg63; T reg128=reg87*reg119; T reg129=(*f.m).resolution*reg97; T reg130=reg87*reg43; T reg131=reg87*reg44;
    T reg132=0.5*reg24; T reg133=0.5*reg98; T reg134=reg57*reg87; T reg135=(*f.m).resolution*reg78; T reg136=0.5*reg124;
    T reg137=0.5*reg107; T reg138=0.5*reg99; reg13=reg113+reg13; reg116=reg110-reg116; reg114=reg111+reg114;
    reg110=(*f.m).resolution*reg4; reg111=(*f.m).resolution*reg70; reg113=0.5*reg103; T reg139=(*f.m).resolution*reg126; T reg140=0.5*reg102;
    T reg141=reg117*reg114; reg131=reg129+reg131; reg129=reg136*reg13; T reg142=reg99*reg116; T reg143=reg137*reg13;
    T reg144=reg102*reg116; reg110=reg120-reg110; reg120=reg113*reg13; T reg145=reg107*reg114; T reg146=reg133*reg13;
    T reg147=reg13*reg105; T reg148=reg103*reg114; T reg149=reg124*reg116; T reg150=reg132*reg13; T reg151=reg98*reg116;
    T reg152=reg138*reg13; reg111=reg134+reg111; reg139=reg2-reg139; reg128=reg135+reg128; reg127=reg130-reg127;
    reg2=reg140*reg13; reg130=reg24*reg114; reg146=reg130+reg146; reg130=reg24*reg131; reg134=reg136*reg110;
    reg135=reg117*reg139; T reg153=reg137*reg110; T reg154=reg99*reg111; T reg155=reg136*reg128; T reg156=reg138*reg110;
    T reg157=reg107*reg139; T reg158=reg117*reg131; T reg159=reg113*reg128; T reg160=reg113*reg110; T reg161=reg102*reg111;
    T reg162=reg102*reg127; T reg163=reg140*reg110; T reg164=reg103*reg139; reg120=reg120-reg144; T reg165=reg132*reg110;
    T reg166=reg98*reg111; reg129=reg129-reg141; reg150=reg151+reg150; reg151=reg137*reg128; T reg167=reg99*reg127;
    T reg168=reg98*reg127; T reg169=reg132*reg128; T reg170=reg107*reg131; T reg171=reg138*reg128; reg143=reg142+reg143;
    reg142=reg128*reg105; T reg172=reg140*reg128; T reg173=reg103*reg131; T reg174=reg133*reg128; T reg175=reg124*reg111;
    T reg176=reg124*reg127; T reg177=reg110*reg105; reg148=reg148-reg2; reg152=reg145+reg152; reg149=reg149-reg147;
    reg120=2*reg120; reg165=reg166+reg165; reg130=reg174+reg130; reg129=2*reg129; reg146=2*reg146;
    reg151=reg167+reg151; reg150=2*reg150; reg148=2*reg148; reg169=reg168+reg169; reg152=2*reg152;
    reg171=reg170+reg171; reg143=2*reg143; reg173=reg173-reg172; reg134=reg134-reg135; reg159=reg159-reg162;
    reg153=reg154+reg153; reg155=reg155-reg158; reg156=reg157+reg156; reg175=reg175-reg177; reg149=2*reg149;
    reg160=reg160-reg161; reg176=reg176-reg142; reg164=reg164-reg163; reg145=reg149*reg113; reg154=reg102*reg175;
    reg157=reg155*reg107; reg166=reg171*reg107; reg167=reg149*reg140; reg168=reg149*reg138; reg170=reg143*reg113;
    reg174=reg102*reg153; T reg178=reg120*reg113; T reg179=reg152*reg138; T reg180=reg149*reg137; T reg181=reg176*reg107;
    T reg182=reg102*reg156; T reg183=reg99*reg175; T reg184=reg99*reg134; T reg185=reg129*reg138; T reg186=reg152*reg113;
    T reg187=reg102*reg160; T reg188=reg129*reg140; T reg189=reg102*reg134; T reg190=reg129*reg113; T reg191=reg143*reg138;
    T reg192=reg99*reg153; T reg193=reg143*reg137; T reg194=reg155*reg103; T reg195=reg129*reg137; T reg196=reg176*reg103;
    T reg197=reg98*reg175; T reg198=reg149*reg132; T reg199=reg171*reg24; T reg200=reg98*reg134; T reg201=reg129*reg132;
    T reg202=reg98*reg153; T reg203=reg143*reg132; T reg204=reg143*reg133; T reg205=reg98*reg156; T reg206=reg152*reg132;
    T reg207=reg149*reg133; T reg208=reg98*reg160; T reg209=reg120*reg132; T reg210=reg98*reg164; T reg211=reg148*reg132;
    T reg212=reg151*reg24; T reg213=reg98*reg165; T reg214=reg169*reg24; T reg215=reg150*reg132; T reg216=reg130*reg24;
    T reg217=reg129*reg133; T reg218=reg155*reg24; T reg219=reg150*reg133; reg175=reg124*reg175; T reg220=reg176*reg24;
    reg176=reg176*reg117; T reg221=reg149*reg136; T reg222=reg173*reg24; reg155=reg155*reg117; T reg223=reg129*reg136;
    reg149=reg149*reg105; T reg224=reg133*reg146; T reg225=reg148*reg133; T reg226=reg143*reg140; T reg227=reg151*reg103;
    T reg228=reg159*reg24; T reg229=reg152*reg140; T reg230=reg171*reg103; T reg231=reg120*reg133; T reg232=reg151*reg107;
    T reg233=reg120*reg140; T reg234=reg159*reg103; T reg235=reg152*reg133; T reg236=reg148*reg140; T reg237=reg173*reg103;
    reg220=reg207+reg220; reg218=reg217+reg218; reg192=reg193+reg192; reg182=reg186-reg182; reg174=reg170-reg174;
    reg212=reg204+reg212; reg189=reg190-reg189; reg181=reg168+reg181; reg199=reg235+reg199; reg154=reg145-reg154;
    reg166=reg179+reg166; reg157=reg185+reg157; reg228=reg231+reg228; reg191=reg232+reg191; reg229=reg230-reg229;
    reg233=reg234-reg233; reg226=reg227-reg226; reg236=reg237-reg236; reg197=reg198+reg197; reg200=reg201+reg200;
    reg202=reg203+reg202; reg205=reg206+reg205; reg208=reg209+reg208; reg210=reg211+reg210; reg213=reg215+reg213;
    reg183=reg180+reg183; reg155=reg223-reg155; reg224=reg216+reg224; reg188=reg194-reg188; reg176=reg221-reg176;
    reg214=reg219+reg214; reg167=reg196-reg167; reg184=reg195+reg184; reg149=reg175-reg149; reg187=reg178-reg187;
    reg222=reg225+reg222; reg183=reg104*reg183; reg157=reg104*reg157; reg181=reg104*reg181; reg184=reg104*reg184;
    reg192=reg104*reg192; reg188=reg104*reg188; reg226=reg104*reg226; reg222=reg222*reg104; reg149=reg104*reg149;
    reg214=reg214*reg104; reg176=reg104*reg176; reg155=reg104*reg155; reg224=reg224*reg104; reg229=reg104*reg229;
    reg213=reg104*reg213; reg233=reg104*reg233; reg236=reg104*reg236; reg210=reg104*reg210; reg197=reg104*reg197;
    reg200=reg104*reg200; reg208=reg104*reg208; reg202=reg104*reg202; reg205=reg104*reg205; reg189=reg104*reg189;
    reg220=reg104*reg220; reg228=reg228*reg104; reg167=reg104*reg167; reg166=reg104*reg166; reg191=reg191*reg104;
    reg187=reg104*reg187; reg218=reg104*reg218; reg199=reg104*reg199; reg154=reg104*reg154; reg182=reg104*reg182;
    reg212=reg104*reg212; reg174=reg104*reg174; T tmp_1_2=ponderation*reg210; T tmp_0_3=ponderation*reg228; T tmp_1_3=ponderation*reg208;
    T tmp_4_5=ponderation*reg191; T tmp_1_1=ponderation*reg213; T tmp_0_4=ponderation*reg199; T tmp_0_0=ponderation*reg224; T tmp_0_5=ponderation*reg212;
    T tmp_0_1=ponderation*reg214; T tmp_0_6=ponderation*reg218; T tmp_0_2=ponderation*reg222; T tmp_0_7=ponderation*reg220; T tmp_4_4=ponderation*reg166;
    T tmp_3_7=ponderation*reg154; T tmp_4_6=ponderation*reg157; T tmp_3_6=ponderation*reg189; T tmp_3_5=ponderation*reg174; T tmp_4_7=ponderation*reg181;
    T tmp_3_4=ponderation*reg182; T tmp_3_3=ponderation*reg187; T tmp_5_5=ponderation*reg192; T tmp_2_7=ponderation*reg167; T tmp_2_6=ponderation*reg188;
    T tmp_5_6=ponderation*reg184; T tmp_2_5=ponderation*reg226; T tmp_7_7=ponderation*reg149; T tmp_5_7=ponderation*reg183; T tmp_6_7=ponderation*reg176;
    T tmp_6_6=ponderation*reg155; T tmp_2_4=ponderation*reg229; T tmp_2_3=ponderation*reg233; T tmp_2_2=ponderation*reg236; T tmp_1_7=ponderation*reg197;
    T tmp_1_6=ponderation*reg200; T tmp_1_5=ponderation*reg202; T tmp_1_4=ponderation*reg205;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_3;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_4;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_5;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_7;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_3;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_4;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_5;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_7;
    matrix(indices[1]+0,indices[1]+0) += tmp_2_2;
    matrix(indices[1]+0,indices[1]+1) += tmp_2_3;
    matrix(indices[1]+0,indices[2]+0) += tmp_2_4;
    matrix(indices[1]+0,indices[2]+1) += tmp_2_5;
    matrix(indices[1]+0,indices[3]+0) += tmp_2_6;
    matrix(indices[1]+0,indices[3]+1) += tmp_2_7;
    matrix(indices[1]+1,indices[1]+1) += tmp_3_3;
    matrix(indices[1]+1,indices[2]+0) += tmp_3_4;
    matrix(indices[1]+1,indices[2]+1) += tmp_3_5;
    matrix(indices[1]+1,indices[3]+0) += tmp_3_6;
    matrix(indices[1]+1,indices[3]+1) += tmp_3_7;
    matrix(indices[2]+0,indices[2]+0) += tmp_4_4;
    matrix(indices[2]+0,indices[2]+1) += tmp_4_5;
    matrix(indices[2]+0,indices[3]+0) += tmp_4_6;
    matrix(indices[2]+0,indices[3]+1) += tmp_4_7;
    matrix(indices[2]+1,indices[2]+1) += tmp_5_5;
    matrix(indices[2]+1,indices[3]+0) += tmp_5_6;
    matrix(indices[2]+1,indices[3]+1) += tmp_5_7;
    matrix(indices[3]+0,indices[3]+0) += tmp_6_6;
    matrix(indices[3]+0,indices[3]+1) += tmp_6_7;
    matrix(indices[3]+1,indices[3]+1) += tmp_7_7;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim,unsigned symmetric_version>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; T reg2=2*(*f.m).shear_modulus_12; reg0=1.0/reg0; reg1=1.0/reg1;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=pow((*f.m).v1[0],2); T reg5=pow((*f.m).v1[1],2); T reg6=reg2*reg3;
    T reg7=pow((*f.m).v2[0],2); T reg8=pow((*f.m).v2[1],2); T reg9=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg10=1.0/(*f.m).elastic_modulus_3; T reg11=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg12=1.0/(*f.m).elastic_modulus_2; T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg14=reg10*reg6; T reg15=reg9*reg6; reg5=reg4+reg5;
    reg4=pow((*f.m).v1[2],2); T reg16=reg11*reg6; reg8=reg7+reg8; reg7=pow((*f.m).v2[2],2); reg4=reg5+reg4;
    reg5=reg12*reg14; T reg17=reg11*reg16; T reg18=reg13*reg14; T reg19=reg11*reg15; reg7=reg8+reg7;
    reg8=1.0/(*f.m).elastic_modulus_1; T reg20=reg12*reg15; reg7=pow(reg7,0.5); T reg21=reg13*reg16; reg19=reg18+reg19;
    reg17=reg5-reg17; reg4=pow(reg4,0.5); reg5=(*f.m).v1[1]/reg4; T reg22=(*f.m).v1[2]/reg4; T reg23=(*f.m).v2[1]/reg7;
    T reg24=(*f.m).v2[2]/reg7; T reg25=reg8*reg17; T reg26=reg13*reg19; T reg27=reg21+reg20; T reg28=reg10*reg3;
    T reg29=reg2*reg0; T reg30=reg13*reg6; T reg31=reg9*reg15; reg14=reg8*reg14; T reg32=reg9*reg16;
    T reg33=reg9*reg3; reg3=reg11*reg3; reg6=reg12*reg6; T reg34=reg9*reg27; reg26=reg25-reg26;
    reg25=reg22*reg23; T reg35=reg5*reg24; reg7=(*f.m).v2[0]/reg7; reg4=(*f.m).v1[0]/reg4; T reg36=reg4*reg24;
    T reg37=reg22*reg7; reg15=reg13*reg15; T reg38=reg35-reg25; T reg39=reg11*reg29; T reg40=2*reg4;
    reg16=reg8*reg16; T reg41=reg9*reg29; T reg42=reg9*reg6; T reg43=reg12*reg28; reg32=reg18+reg32;
    reg28=reg13*reg28; reg18=reg11*reg3; T reg44=reg2*reg1; T reg45=2*reg7; T reg46=reg11*reg33;
    reg31=reg14-reg31; reg34=reg26-reg34; reg29=reg10*reg29; reg14=reg9*reg30; reg26=reg45*reg23;
    T reg47=reg5*reg7; T reg48=reg4*reg23; reg6=reg8*reg6; reg17=reg17/reg34; T reg49=reg37-reg36;
    T reg50=pow(reg4,2); T reg51=pow(reg5,2); T reg52=reg5*reg40; T reg53=2*reg38; reg32=reg32/reg34;
    reg42=reg21+reg42; reg19=reg19/reg34; reg31=reg31/reg34; T reg54=pow(reg7,2); T reg55=pow(reg23,2);
    T reg56=reg12*reg29; reg46=reg28+reg46; reg18=reg43-reg18; reg28=reg9*reg44; reg29=reg13*reg29;
    reg33=reg12*reg33; reg43=reg11*reg39; T reg57=reg11*reg41; reg3=reg13*reg3; T reg58=reg10*reg44;
    reg30=reg13*reg30; reg15=reg16+reg15; reg44=reg11*reg44; reg14=reg16+reg14; reg16=reg32*reg50;
    T reg59=reg3+reg33; T reg60=reg31*reg54; T reg61=reg53*reg49; reg46=reg13*reg46; T reg62=reg32*reg51;
    reg18=reg8*reg18; T reg63=reg31*reg55; reg41=reg12*reg41; T reg64=reg32*reg52; reg39=reg13*reg39;
    reg14=reg14/reg34; T reg65=reg19*reg26; reg43=reg56-reg43; reg56=pow(reg49,2); T reg66=reg19*reg54;
    reg57=reg29+reg57; reg29=reg17*reg50; reg30=reg6-reg30; reg6=reg12*reg58; T reg67=reg17*reg52;
    T reg68=reg31*reg26; reg15=reg15/reg34; T reg69=pow(reg38,2); reg58=reg13*reg58; T reg70=reg11*reg44;
    T reg71=reg11*reg28; T reg72=reg19*reg55; T reg73=reg17*reg51; T reg74=pow(reg24,2); T reg75=pow(reg22,2);
    T reg76=reg48-reg47; reg27=reg27/reg34; reg42=reg42/reg34; T reg77=reg27*reg56; reg72=reg73+reg72;
    reg71=reg58+reg71; reg70=reg6-reg70; reg6=reg39+reg41; reg57=reg13*reg57; reg43=reg8*reg43;
    reg59=reg9*reg59; reg46=reg18-reg46; reg28=reg12*reg28; reg44=reg13*reg44; reg18=reg42*reg51;
    reg58=reg14*reg55; reg73=reg42*reg52; T reg78=reg14*reg26; T reg79=reg17*reg75; T reg80=reg19*reg74;
    reg65=reg67+reg65; reg67=reg27*reg61; reg60=reg16+reg60; reg16=reg15*reg69; reg63=reg62+reg63;
    reg62=reg15*reg56; T reg81=reg32*reg75; T reg82=reg31*reg74; reg68=reg64+reg68; reg64=reg15*reg61;
    T reg83=reg42*reg50; T reg84=reg14*reg54; T reg85=pow(reg76,2); T reg86=reg27*reg69; reg30=reg30/reg34;
    reg66=reg29+reg66; reg62=reg63+reg62; reg16=reg60+reg16; reg67=reg65+reg67; reg29=reg27*reg85;
    reg80=reg79+reg80; reg60=reg7*reg23; reg63=reg4*reg5; reg77=reg72+reg77; reg65=reg44+reg28;
    reg61=reg30*reg61; reg78=reg73+reg78; reg71=reg13*reg71; reg70=reg8*reg70; reg72=reg14*reg74;
    reg73=reg42*reg75; reg79=reg30*reg56; reg58=reg18+reg58; reg6=reg9*reg6; reg57=reg43-reg57;
    reg18=reg30*reg69; reg84=reg83+reg84; reg59=reg46-reg59; reg43=2*reg5; reg46=reg22*reg40;
    reg86=reg66+reg86; reg66=2*reg23; reg83=reg45*reg24; reg82=reg81+reg82; reg64=reg68+reg64;
    reg68=reg15*reg85; reg6=reg57-reg6; reg57=reg5*reg38; reg81=reg4*reg49; T reg87=reg16*reg55;
    T reg88=reg86*reg51; T reg89=reg22*reg43; reg59=reg59/reg34; T reg90=reg77*reg51; T reg91=reg63*reg67;
    T reg92=reg62*reg54; T reg93=reg60*reg62; T reg94=reg77*reg50; T reg95=reg7*reg40; T reg96=reg23*reg43;
    reg62=reg62*reg55; T reg97=reg67*reg51; T reg98=reg64*reg55; reg77=reg63*reg77; reg68=reg82+reg68;
    reg82=reg19*reg83; T reg99=reg17*reg46; T reg100=reg32*reg46; T reg101=reg49*reg38; reg29=reg80+reg29;
    reg80=reg31*reg83; T reg102=reg16*reg54; T reg103=reg86*reg50; T reg104=reg66*reg24; reg65=reg9*reg65;
    reg71=reg70-reg71; reg53=reg53*reg76; reg70=reg64*reg54; T reg105=2*reg49; T reg106=reg4*reg7;
    reg67=reg67*reg50; T reg107=reg5*reg23; reg48=reg47+reg48; reg64=reg60*reg64; reg61=reg78+reg61;
    reg47=2*reg22; reg78=reg30*reg85; reg72=reg73+reg72; reg79=reg58+reg79; reg18=reg84+reg18;
    reg58=reg14*reg83; reg73=reg42*reg46; reg65=reg71-reg65; reg71=reg8*reg54; reg84=reg69*reg61;
    reg78=reg72+reg78; reg70=reg67+reg70; reg67=reg13*reg55; reg6=reg6/reg34; reg72=reg68*reg54;
    T reg108=reg29*reg50; reg47=reg47*reg24; T reg109=reg95*reg8; T reg110=reg69*reg79; reg92=reg94+reg92;
    reg94=reg96*reg13; reg62=reg90+reg62; reg90=reg56*reg79; T reg111=reg29*reg51; T reg112=reg68*reg55;
    T reg113=reg13*reg54; reg87=reg88+reg87; reg88=reg56*reg18; T reg114=reg7*reg49; T reg115=reg23*reg38;
    reg81=reg57+reg81; reg57=reg12*reg55; T reg116=reg5*reg49; T reg117=reg4*reg38; T reg118=reg22*reg24;
    T reg119=reg95*reg13; reg105=reg105*reg76; T reg120=reg96*reg12; T reg121=reg59*reg48; T reg122=reg59*reg107;
    reg31=reg31*reg104; reg32=reg32*reg89; T reg123=reg15*reg53; reg80=reg100+reg80; reg102=reg103+reg102;
    reg100=reg69*reg18; reg103=reg59*reg106; reg19=reg19*reg104; reg17=reg17*reg89; T reg124=reg27*reg53;
    reg82=reg99+reg82; reg99=reg101*reg61; reg16=reg60*reg16; reg86=reg63*reg86; reg64=reg91+reg64;
    reg79=reg101*reg79; reg98=reg97+reg98; reg61=reg56*reg61; reg93=reg77+reg93; reg34=reg65/reg34;
    reg15=reg15*reg105; reg65=reg9*reg54; reg77=reg59*reg118; reg117=reg6*reg117; reg94=reg109-reg94;
    reg91=reg9*reg74; reg116=reg6*reg116; reg119=reg120-reg119; reg67=reg71-reg67; reg71=reg69*reg78;
    reg72=reg108+reg72; reg97=reg47*reg11; reg79=reg93+reg79; reg27=reg27*reg105; reg19=reg17+reg19;
    reg124=reg82+reg124; reg14=reg14*reg104; reg42=reg42*reg89; reg17=reg95*reg103; reg100=reg102+reg100;
    reg53=reg30*reg53; reg58=reg73+reg58; reg73=reg96*reg11; reg16=reg86+reg16; reg82=reg95*reg9;
    reg18=reg101*reg18; reg86=reg95*reg121; reg84=reg70+reg84; reg70=reg11*reg55; reg123=reg80+reg123;
    reg31=reg32+reg31; reg32=reg96*reg103; reg88=reg87+reg88; reg40=reg38*reg40; reg80=reg7*reg38;
    reg87=reg23*reg49; reg43=reg49*reg43; reg68=reg60*reg68; reg61=reg98+reg61; reg90=reg62+reg90;
    reg62=reg96*reg122; reg114=reg115+reg114; reg99=reg64+reg99; reg112=reg111+reg112; reg64=reg56*reg78;
    reg93=reg6*reg81; reg98=reg48*reg121; reg102=reg48*reg122; reg113=reg57-reg113; reg57=reg22*reg76;
    reg108=reg11*reg74; reg122=reg95*reg122; reg110=reg92+reg110; reg29=reg63*reg29; reg92=reg47*reg9;
    reg121=reg96*reg121; reg12=reg12*reg51; reg108=reg113-reg108; reg102=reg79+reg102; reg79=reg40*reg117;
    reg109=reg124*reg50; reg111=reg123*reg54; reg17=reg100+reg17; reg15=reg31+reg15; reg73=reg82+reg73;
    reg37=reg36+reg37; reg47=reg47*reg10; reg31=reg81*reg116; reg57=reg6*reg57; reg78=reg101*reg78;
    reg32=reg88+reg32; reg68=reg29+reg68; reg29=reg43*reg117; reg97=reg119-reg97; reg36=reg4*reg76;
    reg18=reg16+reg18; reg16=reg24*reg76; reg70=reg65+reg70; reg65=reg10*reg74; reg8=reg8*reg50;
    reg103=reg48*reg103; reg82=reg22*reg38; reg88=(*f.m).alpha_2*reg55; reg100=(*f.m).alpha_1*reg51; reg80=reg34*reg80;
    reg113=reg81*reg93; reg45=reg45*reg38; reg115=reg13*reg50; reg119=(*f.m).alpha_2*reg54; reg87=reg34*reg87;
    reg13=reg13*reg51; reg121=reg61+reg121; reg61=reg43*reg93; reg92=reg94-reg92; reg94=reg95*reg77;
    reg122=reg110+reg122; reg71=reg72+reg71; reg72=reg123*reg55; reg110=reg124*reg51; reg120=reg40*reg116;
    reg91=reg67-reg91; reg27=reg19+reg27; reg19=reg96*reg77; reg64=reg112+reg64; reg105=reg30*reg105;
    reg14=reg42+reg14; reg93=reg40*reg93; reg86=reg84+reg86; reg53=reg58+reg53; reg66=reg66*reg49;
    reg30=reg34*reg114; reg42=(*f.m).alpha_1*reg50; reg98=reg99+reg98; reg116=reg43*reg116; reg62=reg90+reg62;
    reg58=reg43*reg57; reg19=reg64+reg19; reg29=reg32+reg29; reg77=reg48*reg77; reg31=reg102+reg31;
    reg32=reg66*reg87; reg64=reg66*reg30; reg67=reg114*reg87; reg61=reg121+reg61; reg78=reg68+reg78;
    reg68=reg66*reg80; reg123=reg60*reg123; reg124=reg63*reg124; reg116=reg62+reg116; reg113=reg98+reg113;
    reg62=reg114*reg30; reg93=reg86+reg93; reg30=reg45*reg30; reg79=reg17+reg79; reg17=(*f.m).alpha_2*reg74;
    reg105=reg14+reg105; reg14=(*f.m).alpha_1*reg75; reg84=reg45*reg80; reg86=reg56*(*f.m).alpha_3; reg88=reg100+reg88;
    reg111=reg109+reg111; reg90=reg69*reg53; reg16=reg34*reg16; reg98=reg69*(*f.m).alpha_3; reg103=reg18+reg103;
    reg117=reg81*reg117; reg119=reg42+reg119; reg13=reg8-reg13; reg8=reg9*reg75; reg18=reg27*reg50;
    reg42=reg15*reg54; reg99=reg40*reg57; reg100=reg59*reg37; reg102=reg15*reg55; reg94=reg71+reg94;
    reg71=reg27*reg51; reg109=reg56*reg53; reg72=reg110+reg72; reg87=reg45*reg87; reg120=reg122+reg120;
    reg110=reg5*reg76; reg112=reg22*reg49; reg121=reg24*reg38; reg122=reg7*reg76; reg36=reg82+reg36;
    reg82=reg11*reg75; reg35=reg25+reg35; reg115=reg12-reg115; reg12=reg92*reg54; reg9=reg9*reg50;
    reg25=reg92*reg106; T reg125=reg97*reg107; reg11=reg11*reg51; T reg126=reg108*reg55; T reg127=reg91*reg54;
    T reg128=reg91*reg50; reg91=reg91*reg106; T reg129=reg108*reg51; reg108=reg108*reg107; T reg130=reg97*reg51;
    reg92=reg92*reg50; reg97=reg97*reg55; reg73=reg47-reg73; reg70=reg65-reg70; reg85=reg85*(*f.m).alpha_3;
    reg47=reg63*(*f.m).alpha_1; reg65=reg60*(*f.m).alpha_2; reg30=reg93+reg30; reg97=reg12+reg97; reg17=reg14+reg17;
    reg108=reg91+reg108; reg12=reg70*reg118; reg86=reg88+reg86; reg98=reg119+reg98; reg125=reg25+reg125;
    reg14=reg73*reg118; reg56=reg56*reg105; reg102=reg71+reg102; reg25=reg96*reg100; reg109=reg72+reg109;
    reg11=reg9+reg11; reg9=reg60*reg2; reg15=reg60*reg15; reg27=reg63*reg27; reg60=reg2*reg48;
    reg53=reg101*reg53; reg123=reg124+reg123; reg10=reg10*reg75; reg62=reg113+reg62; reg59=reg59*reg35;
    reg82=reg115-reg82; reg36=reg6*reg36; reg8=reg13-reg8; reg57=reg81*reg57; reg77=reg78+reg77;
    reg67=reg31+reg67; reg7=reg7*reg24; reg4=reg4*reg22; reg84=reg79+reg84; reg80=reg114*reg80;
    reg117=reg103+reg117; reg87=reg120+reg87; reg99=reg94+reg99; reg13=reg45*reg16; reg31=reg73*reg74;
    reg32=reg116+reg32; reg71=reg70*reg75; reg122=reg121+reg122; reg72=reg24*reg49; reg129=reg128+reg129;
    reg78=reg95*reg100; reg70=reg70*reg74; reg90=reg111+reg90; reg73=reg73*reg75; reg130=reg92+reg130;
    reg110=reg112+reg110; reg79=reg66*reg16; reg58=reg19+reg58; reg64=reg61+reg64; reg126=reg127+reg126;
    reg42=reg18+reg42; reg69=reg69*reg105; reg68=reg29+reg68; reg18=reg23*reg76; reg19=reg32*reg62;
    reg110=reg6*reg110; reg6=reg9*reg26; reg29=reg68*reg98; reg61=reg32*reg86; reg88=reg62*reg87;
    reg69=reg42+reg69; reg71=reg129+reg71; reg11=reg10-reg11; reg38=reg76*reg38; reg22=reg5*reg22;
    reg5=reg8*reg54; reg24=reg23*reg24; reg10=reg1*reg37; reg80=reg117+reg80; reg23=reg82*reg55;
    reg42=reg86*reg87; reg91=reg8*reg50; reg18=reg72+reg18; reg72=reg82*reg51; reg92=reg98*reg84;
    reg57=reg77+reg57; reg16=reg114*reg16; reg122=reg34*reg122; reg77=reg9*reg52; reg13=reg99+reg13;
    reg93=reg60*reg26; reg31=reg97+reg31; reg25=reg109+reg25; reg94=reg43*reg36; reg105=reg101*reg105;
    reg15=reg27+reg15; reg27=reg60*reg48; reg14=reg125+reg14; reg70=reg126+reg70; reg56=reg102+reg56;
    reg96=reg96*reg59; reg73=reg130+reg73; reg97=reg67*reg30; reg78=reg90+reg78; reg90=reg40*reg36;
    reg99=reg7*reg1; reg100=reg48*reg100; reg53=reg123+reg53; reg95=reg95*reg59; reg7=reg7*(*f.m).alpha_2;
    reg102=reg4*(*f.m).alpha_1; reg101=reg101*(*f.m).alpha_3; reg103=reg64*reg67; reg65=reg47+reg65; reg85=reg17+reg85;
    reg79=reg58+reg79; reg12=reg108+reg12; reg9=reg9*reg48; reg60=reg60*reg52; reg59=reg48*reg59;
    reg49=reg76*reg49; reg17=reg98*reg80; reg69=reg95+reg69; reg13=reg85*reg13; reg2=reg63*reg2;
    reg79=reg79*reg85; reg42=reg92+reg42; reg105=reg15+reg105; reg15=reg10*reg46; reg34=reg18*reg34;
    reg40=reg40*reg110; reg18=reg86*reg67; reg16=reg57+reg16; reg36=reg81*reg36; reg100=reg53+reg100;
    reg60=reg73+reg60; reg6=reg70+reg6; reg47=reg24*(*f.m).alpha_2; reg53=reg22*(*f.m).alpha_1; reg38=reg38*(*f.m).alpha_3;
    reg7=reg102+reg7; reg101=reg65+reg101; reg9=reg12+reg9; reg12=reg99*reg37; reg57=reg45*reg122;
    reg103=reg19-reg103; reg90=reg78+reg90; reg97=reg88-reg97; reg43=reg43*reg110; reg96=reg56+reg96;
    reg19=reg64*reg87; reg56=reg66*reg122; reg94=reg25+reg94; reg25=reg32*reg30; reg27=reg14+reg27;
    reg14=reg10*reg37; reg58=reg0*reg35; reg23=reg5+reg23; reg74=reg11*reg74; reg72=reg91+reg72;
    reg24=reg24*reg0; reg75=reg11*reg75; reg82=reg82*reg107; reg8=reg8*reg106; reg5=reg99*reg46;
    reg10=reg10*reg83; reg77=reg71+reg77; reg93=reg31+reg93; reg61=reg29+reg61; reg99=reg99*reg83;
    reg25=reg19-reg25; reg19=reg24*reg89; reg79=reg61+reg79; reg29=reg64*reg101; reg1=reg4*reg1;
    reg4=reg68*reg97; reg5=reg77+reg5; reg31=reg103*reg84; reg75=reg72+reg75; reg61=reg2*reg52;
    reg16=reg85*reg16; reg18=reg17+reg18; reg56=reg94+reg56; reg43=reg96+reg43; reg66=reg66*reg34;
    reg17=reg58*reg89; reg15=reg60+reg15; reg60=reg24*reg35; reg12=reg9+reg12; reg38=reg7+reg38;
    reg47=reg53+reg47; reg49=reg49*(*f.m).alpha_3; reg14=reg27+reg14; reg74=reg23+reg74; reg99=reg6+reg99;
    reg6=reg58*reg35; reg24=reg24*reg104; reg7=reg2*reg26; reg10=reg93+reg10; reg58=reg58*reg104;
    reg45=reg45*reg34; reg40=reg69+reg40; reg118=reg11*reg118; reg57=reg90+reg57; reg82=reg8+reg82;
    reg36=reg100+reg36; reg122=reg114*reg122; reg59=reg105+reg59; reg110=reg81*reg110; reg13=reg42+reg13;
    reg8=reg101*reg30; reg7=reg74+reg7; reg9=reg62*reg84; reg17=reg15+reg17; reg11=reg80*reg30;
    reg15=reg101*reg62; reg16=reg18+reg16; reg45=reg40+reg45; reg19=reg5+reg19; reg8=reg13+reg8;
    reg57=reg38*reg57; reg56=reg56*reg38; reg29=reg79+reg29; reg34=reg114*reg34; reg110=reg59+reg110;
    reg122=reg36+reg122; reg0=reg22*reg0; reg118=reg82+reg118; reg2=reg2*reg48; reg58=reg10+reg58;
    reg24=reg99+reg24; reg61=reg75+reg61; reg46=reg1*reg46; reg49=reg47+reg49; reg60=reg12+reg60;
    reg66=reg43+reg66; reg4=reg31-reg4; reg5=reg80*reg25; reg6=reg14+reg6; reg10=reg68*reg62;
    reg12=reg64*reg80; reg83=reg1*reg83; reg45=reg45*reg49; reg56=reg29+reg56; reg66=reg66*reg49;
    reg15=reg16+reg15; reg122=reg38*reg122; reg89=reg0*reg89; reg46=reg61+reg46; reg13=reg58*reg60;
    reg14=reg24*reg6; reg16=reg19*reg6; reg104=reg0*reg104; reg83=reg7+reg83; reg12=reg10-reg12;
    reg7=reg68*reg67; reg10=reg32*reg80; reg5=reg4+reg5; reg34=reg110+reg34; reg2=reg118+reg2;
    reg37=reg1*reg37; reg1=reg17*reg60; reg4=reg67*reg84; reg11=reg9-reg11; reg9=reg80*reg87;
    reg57=reg8+reg57; reg8=reg68*reg30; reg18=reg64*reg84; reg122=reg15+reg122; reg34=reg49*reg34;
    reg104=reg83+reg104; reg1=reg16-reg1; reg15=reg19*reg58; reg16=reg17*reg24; reg35=reg0*reg35;
    reg37=reg2+reg37; reg8=reg18-reg8; reg13=reg14-reg13; reg0=reg68*reg87; reg89=reg46+reg89;
    reg103=reg103/reg5; reg12=reg12/reg5; reg97=reg97/reg5; reg10=reg7-reg10; reg45=reg57+reg45;
    reg2=1-var_inter[1]; reg66=reg56+reg66; reg7=1-var_inter[0]; reg9=reg4-reg9; reg4=reg32*reg84;
    reg11=reg11/reg5; reg8=reg8/reg5; reg16=reg15-reg16; reg25=reg25/reg5; reg9=reg9/reg5;
    reg14=reg2*elem.pos(1)[0]; reg35=reg37+reg35; reg15=reg13*reg89; reg11=reg11*reg66; reg12=reg12*reg45;
    reg18=reg7*elem.pos(0)[1]; reg22=var_inter[0]*elem.pos(1)[1]; reg97=reg97*reg66; reg0=reg4-reg0; reg10=reg10/reg5;
    reg34=reg122+reg34; reg4=reg104*reg1; reg23=reg2*elem.pos(0)[1]; reg27=reg2*elem.pos(1)[1]; reg29=var_inter[0]*elem.pos(1)[0];
    reg31=reg7*elem.pos(0)[0]; reg36=reg2*elem.pos(0)[0]; reg103=reg103*reg45; reg37=var_inter[0]*elem.pos(2)[1]; reg38=reg89*reg60;
    reg40=reg19*reg35; reg42=reg22+reg18; reg5=reg0/reg5; reg0=reg24*reg35; reg60=reg104*reg60;
    reg43=reg31+reg29; reg46=var_inter[0]*elem.pos(2)[0]; reg47=var_inter[1]*elem.pos(2)[1]; reg27=reg27-reg23; reg14=reg14-reg36;
    reg25=reg25*reg34; reg97=reg103-reg97; reg8=reg8*reg34; reg12=reg11-reg12; reg45=reg10*reg45;
    reg66=reg9*reg66; reg9=reg35*reg16; reg4=reg15-reg4; reg10=var_inter[1]*elem.pos(2)[0]; reg97=reg25+reg97;
    reg11=reg58*reg35; reg15=reg7*elem.pos(3)[0]; reg46=reg46-reg43; reg25=reg104*reg6; reg14=reg10+reg14;
    reg10=var_inter[1]*elem.pos(3)[0]; reg35=reg17*reg35; reg40=reg38-reg40; reg24=reg24*reg89; reg19=reg19*reg104;
    reg47=reg27+reg47; reg27=1-(*f.m).resolution; reg38=var_inter[1]*elem.pos(3)[1]; reg9=reg4+reg9; reg6=reg89*reg6;
    reg37=reg37-reg42; reg0=reg60-reg0; reg8=reg12-reg8; reg4=reg7*elem.pos(3)[1]; reg66=reg45-reg66;
    reg34=reg5*reg34; reg0=reg0/reg9; reg47=reg47-reg38; reg97=reg27*reg97; reg66=reg34+reg66;
    reg5=(*f.m).resolution*reg86; reg19=reg24-reg19; reg12=(*f.m).resolution*reg98; reg104=reg17*reg104; reg11=reg25-reg11;
    reg15=reg46+reg15; reg8=reg27*reg8; reg4=reg37+reg4; reg35=reg6-reg35; reg40=reg40/reg9;
    reg14=reg14-reg10; reg89=reg58*reg89; reg19=reg19/reg9; reg13=reg13/reg9; reg35=reg35/reg9;
    reg6=reg4*reg14; reg104=reg89-reg104; reg1=reg1/reg9; reg17=reg15*reg47; reg11=reg11/reg9;
    reg24=(*f.m).resolution*reg101; reg8=reg5+reg8; reg97=reg12+reg97; reg66=reg27*reg66; reg5=(*f.m).resolution*reg0;
    reg12=(*f.m).resolution*reg40; reg80=reg80*reg27; reg67=reg67*reg27; reg25=(*f.m).resolution*reg13; reg34=(*f.m).resolution*reg11;
    reg37=(*f.m).resolution*reg35; reg45=(*f.m).resolution*reg19; reg5=reg80+reg5; reg12=reg67-reg12; reg66=reg24+reg66;
    reg97=reg97*(*f.m).deltaT; reg8=reg8*(*f.m).deltaT; reg17=reg6-reg17; reg84=reg27*reg84; reg87=reg27*reg87;
    reg68=reg68*reg27; reg32=reg32*reg27; reg6=(*f.m).resolution*reg1; reg16=reg16/reg9; reg62=reg62*reg27;
    reg9=reg104/reg9; reg24=(*f.m).resolution*reg16; reg30=reg27*reg30; reg45=reg62+reg45; reg27=reg64*reg27;
    reg37=reg32+reg37; reg34=reg68-reg34; reg4=reg4/reg17; reg6=reg87-reg6; reg84=reg25+reg84;
    reg25=(*f.m).resolution*reg9; reg14=reg14/reg17; reg15=reg15/reg17; reg47=reg47/reg17; reg66=reg66*(*f.m).deltaT;
    reg32=reg5*reg97; reg46=reg12*reg8; reg49=reg37*reg8; reg53=reg6*reg8; reg56=reg34*reg97;
    reg57=reg32+reg46; reg58=reg45*reg66; reg25=reg27-reg25; reg30=reg24+reg30; reg24=var_inter[1]*reg4;
    reg27=reg2*reg4; reg59=reg7*reg47; reg60=var_inter[1]*reg15; reg61=reg7*reg14; reg62=reg2*reg15;
    reg63=var_inter[0]*reg47; reg64=reg84*reg97; reg65=var_inter[0]*reg14; reg67=reg53+reg64; reg68=reg30*reg66;
    reg69=reg57+reg58; reg70=reg62+reg65; reg71=reg27+reg63; reg72=reg59+reg24; reg73=reg61+reg60;
    reg74=reg56+reg49; reg75=reg25*reg66; reg76=reg74+reg75; reg77=reg24-reg63; reg78=reg65-reg60;
    reg79=2*reg69; reg80=0.5*reg71; reg81=0.5*reg70; reg82=reg67+reg68; reg83=0.5*reg73;
    reg85=0.5*reg72; reg87=reg59-reg27; reg88=reg62-reg61; reg89=0.5*reg77; reg90=0.5*reg78;
    reg91=reg80*reg79; reg92=reg81*reg79; reg93=reg71*reg82; reg94=reg70*reg76; reg95=0.5*reg87;
    reg96=0.5*reg88; reg99=reg83*reg79; reg100=reg72*reg82; reg102=reg2*var_inter[0]; reg103=reg7*var_inter[1];
    reg104=reg73*reg76; reg105=reg79*reg85; reg108=reg2*reg7; reg109=reg100-reg99; reg110=reg103*elem.f_vol_e[0];
    reg111=reg89*reg79; reg112=reg102*elem.f_vol_e[1]; reg113=reg78*reg76; reg114=reg103*elem.f_vol_e[1]; reg115=var_inter[0]*var_inter[1];
    reg116=reg94-reg91; reg117=reg102*elem.f_vol_e[0]; reg118=reg88*reg76; reg119=reg77*reg82; reg120=reg87*reg82;
    reg121=reg95*reg79; reg122=reg92-reg93; reg123=reg96*reg79; reg124=reg90*reg79; reg125=reg105-reg104;
    reg116=reg116-reg112; reg126=reg111+reg113; reg127=reg124+reg119; reg128=reg115*elem.f_vol_e[1]; reg129=reg115*elem.f_vol_e[0];
    reg130=reg108*elem.f_vol_e[0]; T reg131=reg108*elem.f_vol_e[1]; reg109=reg109-reg110; T reg132=reg123+reg120; reg125=reg125-reg114;
    reg122=reg122-reg117; T reg133=reg121+reg118; reg109=reg109*reg17; T reg134=reg130+reg132; T reg135=reg128+reg126;
    reg125=reg125*reg17; T reg136=reg129+reg127; T reg137=reg131+reg133; reg116=reg17*reg116; reg122=reg17*reg122;
    T reg138=reg17*reg137; T reg139=reg17*reg136; reg125=ponderation*reg125; T reg140=reg17*reg134; reg109=ponderation*reg109;
    T reg141=reg17*reg135; reg116=ponderation*reg116; reg122=ponderation*reg122; T vec_2=-reg122; reg122=ponderation*reg141;
    T vec_5=reg122; T reg142=ponderation*reg140; T vec_0=reg142; T vec_6=-reg109; T vec_3=-reg116;
    reg109=ponderation*reg138; T vec_1=reg109; reg116=ponderation*reg139; T vec_4=reg116; T vec_7=-reg125;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
sollicitation[indices[0]+0] += vec_0;
sollicitation[indices[0]+1] += vec_1;
sollicitation[indices[1]+0] += vec_2;
sollicitation[indices[1]+1] += vec_3;
sollicitation[indices[2]+0] += vec_4;
sollicitation[indices[2]+1] += vec_5;
sollicitation[indices[3]+0] += vec_6;
sollicitation[indices[3]+1] += vec_7;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim,unsigned symmetric_version>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_residual( TD ponderation, const TD *var_inter,
      const Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices ) { 
  #define PNODE(N) (*elem.node(N))
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; T reg2=2*(*f.m).shear_modulus_12; reg1=1.0/reg1;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg5=reg2*reg3; T reg6=1.0/(*f.m).elastic_modulus_3;
    T reg7=pow((*f.m).v1[0],2); T reg8=pow((*f.m).v1[1],2); T reg9=pow((*f.m).v2[0],2); T reg10=pow((*f.m).v2[1],2); T reg11=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg12=pow((*f.m).v2[2],2); reg10=reg9+reg10; reg9=reg11*reg5; T reg13=reg4*reg5; T reg14=pow((*f.m).v1[2],2);
    reg8=reg7+reg8; reg7=reg6*reg5; T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg17=reg4*reg9;
    T reg18=reg16*reg7; T reg19=reg4*reg13; T reg20=reg15*reg7; reg12=reg10+reg12; reg14=reg8+reg14;
    reg17=reg18+reg17; reg8=reg16*reg13; reg14=pow(reg14,0.5); reg12=pow(reg12,0.5); reg10=reg15*reg9;
    T reg21=1.0/(*f.m).elastic_modulus_1; reg19=reg20-reg19; reg20=reg16*reg17; T reg22=reg8+reg10; T reg23=(*f.m).v1[2]/reg14;
    T reg24=(*f.m).v1[1]/reg14; T reg25=reg21*reg19; T reg26=(*f.m).v2[1]/reg12; T reg27=(*f.m).v2[2]/reg12; T reg28=reg2*reg0;
    T reg29=reg6*reg3; T reg30=reg11*reg3; T reg31=reg24*reg27; T reg32=reg23*reg26; reg14=(*f.m).v1[0]/reg14;
    T reg33=reg11*reg9; reg12=(*f.m).v2[0]/reg12; T reg34=reg16*reg5; reg7=reg21*reg7; T reg35=reg11*reg13;
    reg5=reg15*reg5; reg3=reg4*reg3; T reg36=reg11*reg22; reg20=reg25-reg20; reg25=reg2*reg1;
    T reg37=2*reg12; T reg38=reg14*reg27; T reg39=reg4*reg28; reg9=reg16*reg9; T reg40=reg11*reg34;
    T reg41=reg23*reg12; T reg42=reg11*reg28; T reg43=reg15*reg29; reg36=reg20-reg36; reg13=reg21*reg13;
    reg29=reg16*reg29; reg20=reg4*reg3; T reg44=reg31-reg32; T reg45=reg4*reg30; reg33=reg7-reg33;
    reg35=reg18+reg35; reg28=reg6*reg28; reg7=2*reg14; reg18=reg11*reg5; reg5=reg21*reg5;
    reg19=reg19/reg36; reg35=reg35/reg36; reg18=reg8+reg18; reg17=reg17/reg36; reg34=reg16*reg34;
    reg33=reg33/reg36; reg9=reg13+reg9; reg40=reg13+reg40; reg13=reg4*reg25; T reg46=reg24*reg12;
    T reg47=reg14*reg26; T reg48=reg41-reg38; T reg49=pow(reg14,2); T reg50=pow(reg24,2); T reg51=reg37*reg26;
    T reg52=reg6*reg25; T reg53=reg24*reg7; T reg54=reg4*reg42; T reg55=reg4*reg39; T reg56=reg16*reg28;
    reg28=reg15*reg28; reg45=reg29+reg45; reg20=reg43-reg20; reg29=2*reg44; reg25=reg11*reg25;
    reg43=pow(reg12,2); T reg57=pow(reg26,2); reg30=reg15*reg30; reg3=reg16*reg3; T reg58=reg35*reg53;
    T reg59=reg33*reg51; T reg60=reg33*reg57; T reg61=reg35*reg50; T reg62=pow(reg44,2); reg9=reg9/reg36;
    T reg63=pow(reg48,2); T reg64=reg29*reg48; reg34=reg5-reg34; reg39=reg16*reg39; reg42=reg15*reg42;
    reg20=reg21*reg20; reg45=reg16*reg45; reg5=reg3+reg30; reg55=reg28-reg55; reg54=reg56+reg54;
    reg28=reg15*reg52; reg52=reg16*reg52; reg56=reg4*reg13; T reg65=reg4*reg25; T reg66=reg19*reg50;
    T reg67=reg17*reg57; T reg68=reg19*reg53; T reg69=reg17*reg51; T reg70=reg35*reg49; T reg71=reg33*reg43;
    reg40=reg40/reg36; T reg72=reg17*reg43; T reg73=reg19*reg49; reg22=reg22/reg36; reg18=reg18/reg36;
    T reg74=pow(reg23,2); T reg75=reg47-reg46; T reg76=pow(reg27,2); T reg77=reg22*reg63; T reg78=reg19*reg74;
    T reg79=reg17*reg76; T reg80=reg40*reg43; T reg81=reg18*reg49; T reg82=reg18*reg50; T reg83=reg40*reg57;
    T reg84=reg9*reg64; reg59=reg58+reg59; reg58=reg18*reg53; T reg85=reg33*reg76; T reg86=reg35*reg74;
    T reg87=reg40*reg51; T reg88=reg9*reg63; reg60=reg61+reg60; reg61=reg9*reg62; reg71=reg70+reg71;
    reg69=reg68+reg69; reg68=reg22*reg64; reg34=reg34/reg36; reg70=reg22*reg62; reg72=reg73+reg72;
    reg13=reg16*reg13; reg25=reg15*reg25; reg45=reg20-reg45; reg5=reg11*reg5; reg55=reg21*reg55;
    reg54=reg16*reg54; reg67=reg66+reg67; reg20=reg39+reg42; reg65=reg52+reg65; reg52=pow(reg75,2);
    reg56=reg28-reg56; reg68=reg69+reg68; reg28=reg12*reg26; reg66=reg22*reg52; reg64=reg34*reg64;
    reg69=reg13+reg25; reg77=reg67+reg77; reg87=reg58+reg87; reg58=reg14*reg24; reg65=reg16*reg65;
    reg5=reg45-reg5; reg79=reg78+reg79; reg45=reg40*reg76; reg56=reg21*reg56; reg54=reg55-reg54;
    reg55=reg18*reg74; reg67=reg34*reg63; reg83=reg82+reg83; reg20=reg11*reg20; reg73=reg34*reg62;
    reg80=reg81+reg80; reg78=2*reg24; reg81=reg23*reg7; reg82=2*reg26; reg70=reg72+reg70;
    reg61=reg71+reg61; reg88=reg60+reg88; reg85=reg86+reg85; reg60=reg9*reg52; reg84=reg59+reg84;
    reg59=reg37*reg27; reg71=reg24*reg26; reg72=reg14*reg12; reg86=reg68*reg49; T reg89=2*reg48;
    reg20=reg54-reg20; reg54=reg84*reg43; reg29=reg29*reg75; T reg90=reg58*reg68; T reg91=reg28*reg84;
    T reg92=reg14*reg48; T reg93=reg24*reg44; T reg94=reg28*reg88; T reg95=reg58*reg77; reg47=reg46+reg47;
    reg5=reg5/reg36; reg46=reg77*reg49; T reg96=reg88*reg43; T reg97=reg48*reg44; T reg98=reg70*reg50;
    T reg99=reg61*reg57; reg77=reg77*reg50; reg88=reg88*reg57; reg68=reg68*reg50; reg84=reg84*reg57;
    reg64=reg87+reg64; reg87=reg34*reg52; reg45=reg55+reg45; reg60=reg85+reg60; reg66=reg79+reg66;
    reg67=reg83+reg67; reg55=reg35*reg81; reg79=reg33*reg59; reg73=reg80+reg73; reg65=reg56-reg65;
    reg56=reg23*reg78; reg69=reg11*reg69; reg80=reg82*reg27; reg83=reg19*reg81; reg85=reg61*reg43;
    T reg100=reg70*reg49; T reg101=reg17*reg59; reg101=reg83+reg101; reg83=reg22*reg29; T reg102=reg63*reg73;
    reg88=reg77+reg88; reg77=reg63*reg67; reg84=reg68+reg84; reg99=reg98+reg99; reg68=reg63*reg64;
    reg19=reg19*reg56; reg98=reg66*reg50; T reg103=reg12*reg48; T reg104=reg26*reg44; T reg105=reg60*reg57;
    T reg106=reg66*reg49; reg69=reg65-reg69; reg65=reg62*reg64; reg54=reg86+reg54; reg33=reg33*reg80;
    reg89=reg89*reg75; reg35=reg35*reg56; reg86=reg9*reg29; reg20=reg20/reg36; T reg107=reg23*reg27;
    reg79=reg55+reg79; reg55=reg12*reg7; T reg108=reg26*reg78; T reg109=reg14*reg44; T reg110=2*reg23;
    reg92=reg93+reg92; reg93=reg24*reg48; reg17=reg17*reg80; T reg111=reg60*reg43; reg70=reg58*reg70;
    reg61=reg28*reg61; reg94=reg95+reg94; reg95=reg97*reg67; reg91=reg90+reg91; reg64=reg97*reg64;
    reg67=reg62*reg67; reg87=reg45+reg87; reg45=reg18*reg81; reg90=reg40*reg59; T reg112=reg5*reg72;
    T reg113=reg5*reg71; T reg114=reg5*reg47; reg96=reg46+reg96; reg46=reg62*reg73; reg85=reg100+reg85;
    reg64=reg91+reg64; reg91=reg55*reg114; reg100=reg108*reg16; reg36=reg69/reg36; reg69=reg55*reg21;
    reg65=reg54+reg65; reg60=reg28*reg60; reg66=reg58*reg66; reg54=reg47*reg113; reg95=reg94+reg95;
    reg94=reg62*reg87; reg111=reg106+reg111; reg106=reg16*reg57; T reg115=reg21*reg43; reg7=reg44*reg7;
    reg78=reg48*reg78; reg73=reg97*reg73; reg110=reg110*reg27; reg61=reg70+reg61; reg70=reg23*reg75;
    reg9=reg9*reg89; reg33=reg35+reg33; reg35=reg108*reg15; reg86=reg79+reg86; reg79=reg20*reg92;
    reg93=reg20*reg93; reg109=reg20*reg109; T reg116=reg55*reg16; T reg117=reg5*reg107; T reg118=reg15*reg57;
    T reg119=reg16*reg43; reg22=reg22*reg89; reg17=reg19+reg17; reg83=reg101+reg83; reg40=reg40*reg80;
    reg18=reg18*reg56; reg29=reg34*reg29; reg90=reg45+reg90; reg46=reg85+reg46; reg19=reg55*reg112;
    reg67=reg96+reg67; reg45=reg47*reg114; reg85=reg26*reg48; reg96=reg108*reg113; reg77=reg88+reg77;
    reg88=reg12*reg44; reg103=reg104+reg103; reg105=reg98+reg105; reg98=reg63*reg87; reg101=reg108*reg112;
    reg114=reg108*reg114; reg102=reg99+reg102; reg68=reg84+reg68; reg113=reg55*reg113; reg84=reg11*reg43;
    reg70=reg20*reg70; reg99=reg4*reg57; reg104=(*f.m).alpha_2*reg43; T reg120=reg55*reg11; T reg121=reg108*reg4;
    T reg122=reg86*reg43; T reg123=reg83*reg49; reg22=reg17+reg22; reg89=reg34*reg89; reg17=(*f.m).alpha_1*reg49;
    reg34=reg78*reg109; reg40=reg18+reg40; reg101=reg102+reg101; reg29=reg90+reg29; reg41=reg38+reg41;
    reg114=reg68+reg114; reg18=reg78*reg79; reg38=reg11*reg76; reg106=reg115-reg106; reg68=reg7*reg109;
    reg90=reg110*reg11; reg100=reg69-reg100; reg19=reg46+reg19; reg46=reg108*reg117; reg98=reg105+reg98;
    reg69=reg4*reg76; reg119=reg118-reg119; reg102=reg83*reg50; reg105=reg110*reg4; reg9=reg33+reg9;
    reg33=reg86*reg57; reg116=reg35-reg116; reg35=reg78*reg93; reg115=reg36*reg103; reg85=reg36*reg85;
    reg88=reg36*reg88; reg96=reg77+reg96; reg77=reg14*reg75; reg118=reg7*reg93; T reg124=(*f.m).alpha_1*reg50;
    reg37=reg37*reg44; reg87=reg97*reg87; reg60=reg66+reg60; reg112=reg47*reg112; reg73=reg61+reg73;
    reg94=reg111+reg94; reg61=reg27*reg75; reg66=reg55*reg117; reg82=reg82*reg48; reg113=reg67+reg113;
    reg45=reg64+reg45; reg64=reg92*reg79; reg67=(*f.m).alpha_2*reg57; reg91=reg65+reg91; reg79=reg7*reg79;
    reg93=reg92*reg93; reg65=reg23*reg44; reg54=reg95+reg54; reg95=reg24*reg75; reg111=reg23*reg48;
    T reg125=reg82*reg85; reg35=reg96+reg35; reg77=reg65+reg77; reg66=reg94+reg66; reg109=reg92*reg109;
    reg118=reg113+reg118; reg46=reg98+reg46; reg65=reg78*reg70; reg112=reg73+reg112; reg73=reg82*reg115;
    reg18=reg114+reg18; reg94=reg37*reg85; reg67=reg124+reg67; reg68=reg19+reg68; reg19=reg63*(*f.m).alpha_3;
    reg96=(*f.m).alpha_1*reg74; reg98=reg37*reg88; reg113=(*f.m).alpha_2*reg76; reg104=reg17+reg104; reg86=reg28*reg86;
    reg79=reg91+reg79; reg17=reg37*reg115; reg91=reg16*reg50; reg83=reg58*reg83; reg115=reg103*reg115;
    reg64=reg45+reg64; reg89=reg40+reg89; reg34=reg101+reg34; reg40=reg82*reg88; reg45=reg12*reg75;
    reg101=reg27*reg44; reg114=reg9*reg57; reg124=reg22*reg50; T reg126=reg63*reg29; reg122=reg123+reg122;
    reg123=reg62*reg29; reg33=reg102+reg33; reg102=reg62*(*f.m).alpha_3; T reg127=reg5*reg41; reg117=reg47*reg117;
    reg87=reg60+reg87; reg60=reg22*reg49; T reg128=reg9*reg43; reg85=reg103*reg85; reg93=reg54+reg93;
    reg61=reg36*reg61; reg54=reg7*reg70; reg21=reg21*reg49; reg121=reg120+reg121; reg110=reg110*reg6;
    reg99=reg84+reg99; reg84=reg6*reg76; reg38=reg106-reg38; reg90=reg100-reg90; reg16=reg16*reg49;
    reg15=reg15*reg50; reg69=reg119-reg69; reg105=reg116-reg105; reg31=reg32+reg31; reg32=reg82*reg61;
    reg65=reg46+reg65; reg46=reg90*reg72; reg100=reg105*reg71; reg125=reg35+reg125; reg35=reg69*reg50;
    reg106=reg38*reg49; reg116=reg11*reg49; reg121=reg110-reg121; reg40=reg34+reg40; reg5=reg5*reg31;
    reg34=reg4*reg50; reg77=reg20*reg77; reg99=reg84-reg99; reg11=reg11*reg74; reg91=reg21-reg91;
    reg102=reg104+reg102; reg19=reg67+reg19; reg113=reg96+reg113; reg52=reg52*(*f.m).alpha_3; reg21=reg58*(*f.m).alpha_1;
    reg67=reg28*(*f.m).alpha_2; reg63=reg63*reg89; reg114=reg124+reg114; reg109=reg112+reg109; reg88=reg103*reg88;
    reg85=reg93+reg85; reg117=reg87+reg117; reg84=reg108*reg127; reg126=reg33+reg126; reg70=reg92*reg70;
    reg115=reg64+reg115; reg73=reg18+reg73; reg86=reg83+reg86; reg29=reg97*reg29; reg22=reg58*reg22;
    reg18=reg38*reg72; reg33=reg69*reg71; reg9=reg28*reg9; reg95=reg111+reg95; reg62=reg62*reg89;
    reg64=reg37*reg61; reg54=reg66+reg54; reg66=reg26*reg75; reg83=reg27*reg48; reg69=reg69*reg57;
    reg14=reg14*reg23; reg128=reg60+reg128; reg12=reg12*reg27; reg94=reg118+reg94; reg38=reg38*reg43;
    reg45=reg101+reg45; reg60=reg55*reg127; reg98=reg68+reg98; reg123=reg122+reg123; reg4=reg4*reg74;
    reg68=reg105*reg50; reg16=reg15-reg16; reg15=reg90*reg49; reg105=reg105*reg57; reg17=reg79+reg17;
    reg90=reg90*reg43; reg79=reg99*reg76; reg87=reg7*reg77; reg60=reg123+reg60; reg70=reg117+reg70;
    reg61=reg103*reg61; reg69=reg38+reg69; reg88=reg109+reg88; reg38=reg12*(*f.m).alpha_2; reg93=reg14*(*f.m).alpha_1;
    reg62=reg128+reg62; reg96=reg97*(*f.m).alpha_3; reg67=reg21+reg67; reg21=reg121*reg74; reg52=reg113+reg52;
    reg68=reg15+reg68; reg45=reg36*reg45; reg4=reg16-reg4; reg6=reg6*reg74; reg95=reg20*reg95;
    reg34=reg116+reg34; reg11=reg91-reg11; reg27=reg26*reg27; reg23=reg24*reg23; reg44=reg75*reg44;
    reg28=reg28*reg2; reg15=reg2*reg47; reg66=reg83+reg66; reg35=reg106+reg35; reg16=reg99*reg74;
    reg20=reg125*reg19; reg64=reg54+reg64; reg24=reg40*reg102; reg26=reg19*reg94; reg54=reg102*reg98;
    reg83=reg85*reg17; reg91=reg73*reg85; reg101=reg115*reg94; reg104=reg125*reg115; reg89=reg97*reg89;
    reg9=reg22+reg9; reg22=reg121*reg76; reg105=reg90+reg105; reg127=reg47*reg127; reg29=reg86+reg29;
    reg86=reg78*reg77; reg32=reg65+reg32; reg100=reg46+reg100; reg84=reg126+reg84; reg55=reg55*reg5;
    reg99=reg99*reg107; reg33=reg18+reg33; reg108=reg108*reg5; reg63=reg114+reg63; reg121=reg121*reg107;
    reg18=reg28*reg53; reg16=reg35+reg16; reg48=reg75*reg48; reg20=reg24+reg20; reg24=reg37*reg45;
    reg87=reg60+reg87; reg32=reg32*reg52; reg79=reg69+reg79; reg35=reg28*reg51; reg36=reg66*reg36;
    reg7=reg7*reg95; reg62=reg55+reg62; reg121=reg100+reg121; reg86=reg84+reg86; reg46=reg15*reg47;
    reg55=reg15*reg51; reg78=reg78*reg95; reg108=reg63+reg108; reg22=reg105+reg22; reg60=reg4*reg50;
    reg96=reg67+reg96; reg63=reg82*reg45; reg38=reg93+reg38; reg44=reg44*(*f.m).alpha_3; reg65=reg23*(*f.m).alpha_1;
    reg66=reg27*(*f.m).alpha_2; reg67=reg11*reg49; reg69=reg102*reg88; reg75=reg19*reg85; reg12=reg12*reg1;
    reg84=reg1*reg41; reg26=reg54+reg26; reg54=1-var_inter[1]; reg90=1-var_inter[0]; reg64=reg52*reg64;
    reg127=reg29+reg127; reg77=reg92*reg77; reg21=reg68+reg21; reg15=reg15*reg53; reg61=reg70+reg61;
    reg99=reg33+reg99; reg29=reg73*reg94; reg5=reg47*reg5; reg83=reg101-reg83; reg89=reg9+reg89;
    reg91=reg104-reg91; reg34=reg6-reg34; reg6=reg4*reg57; reg9=reg11*reg43; reg33=reg125*reg17;
    reg28=reg28*reg47; reg5=reg89+reg5; reg95=reg92*reg95; reg68=reg54*elem.pos(0)[1]; reg46=reg121+reg46;
    reg70=reg84*reg41; reg89=reg54*elem.pos(1)[1]; reg74=reg34*reg74; reg60=reg67+reg60; reg67=reg90*elem.pos(0)[1];
    reg92=var_inter[0]*elem.pos(1)[1]; reg2=reg58*reg2; reg37=reg37*reg36; reg58=reg0*reg31; reg77=reg127+reg77;
    reg45=reg103*reg45; reg63=reg86+reg63; reg86=reg54*elem.pos(1)[0]; reg93=reg54*elem.pos(0)[0]; reg27=reg27*reg0;
    reg7=reg62+reg7; reg62=reg90*elem.pos(0)[0]; reg82=reg82*reg36; reg97=reg73*reg96; reg32=reg20+reg32;
    reg20=reg12*reg59; reg35=reg79+reg35; reg24=reg87+reg24; reg28=reg99+reg28; reg79=reg91*reg98;
    reg87=reg40*reg83; reg99=reg12*reg41; reg76=reg34*reg76; reg6=reg9+reg6; reg33=reg29-reg33;
    reg9=reg84*reg81; reg15=reg21+reg15; reg21=reg96*reg17; reg64=reg26+reg64; reg78=reg108+reg78;
    reg55=reg22+reg55; reg18=reg16+reg18; reg11=reg11*reg72; reg4=reg4*reg71; reg44=reg38+reg44;
    reg66=reg65+reg66; reg48=reg48*(*f.m).alpha_3; reg12=reg12*reg81; reg75=reg69+reg75; reg61=reg52*reg61;
    reg16=var_inter[0]*elem.pos(1)[0]; reg84=reg84*reg59; reg21=reg64+reg21; reg22=reg73*reg88; reg45=reg77+reg45;
    reg26=reg40*reg115; reg29=reg115*reg98; reg12=reg18+reg12; reg18=reg88*reg17; reg38=reg27*reg31;
    reg99=reg28+reg99; reg28=reg27*reg56; reg52=reg96*reg115; reg64=reg88*reg33; reg65=reg2*reg53;
    reg74=reg60+reg74; reg87=reg79-reg87; reg61=reg75+reg61; reg24=reg44*reg24; reg97=reg32+reg97;
    reg63=reg63*reg44; reg95=reg5+reg95; reg36=reg103*reg36; reg5=reg62+reg16; reg32=reg92+reg67;
    reg60=var_inter[0]*elem.pos(2)[1]; reg69=var_inter[1]*elem.pos(2)[1]; reg75=reg58*reg80; reg84=reg55+reg84; reg55=var_inter[0]*elem.pos(2)[0];
    reg27=reg27*reg80; reg48=reg66+reg48; reg82=reg78+reg82; reg20=reg35+reg20; reg35=reg2*reg51;
    reg76=reg6+reg76; reg6=reg58*reg56; reg9=reg15+reg9; reg86=reg86-reg93; reg15=var_inter[1]*elem.pos(2)[0];
    reg58=reg58*reg31; reg107=reg34*reg107; reg70=reg46+reg70; reg4=reg11+reg4; reg1=reg14*reg1;
    reg37=reg7+reg37; reg89=reg89-reg68; reg86=reg15+reg86; reg27=reg20+reg27; reg82=reg82*reg48;
    reg63=reg97+reg63; reg58=reg70+reg58; reg7=reg40*reg17; reg38=reg99+reg38; reg37=reg37*reg48;
    reg24=reg21+reg24; reg11=reg73*reg98; reg59=reg1*reg59; reg35=reg76+reg35; reg64=reg87+reg64;
    reg14=reg88*reg94; reg0=reg23*reg0; reg6=reg9+reg6; reg2=reg2*reg47; reg9=reg40*reg85;
    reg22=reg26-reg22; reg28=reg12+reg28; reg12=reg125*reg88; reg18=reg29-reg18; reg15=reg85*reg98;
    reg69=reg89+reg69; reg107=reg4+reg107; reg36=reg95+reg36; reg4=var_inter[1]*elem.pos(3)[1]; reg81=reg1*reg81;
    reg55=reg55-reg5; reg20=reg90*elem.pos(3)[0]; reg65=reg74+reg65; reg21=reg90*elem.pos(3)[1]; reg60=reg60-reg32;
    reg23=var_inter[1]*elem.pos(3)[0]; reg75=reg84+reg75; reg45=reg44*reg45; reg52=reg61+reg52; reg91=reg91/reg64;
    reg41=reg1*reg41; reg1=reg54*vectors[0][indices[0]+1]; reg22=reg22/reg64; reg14=reg15-reg14; reg15=reg54*vectors[0][indices[1]+0];
    reg26=reg6*reg38; reg56=reg0*reg56; reg81=reg65+reg81; reg2=reg107+reg2; reg12=reg9-reg12;
    reg83=reg83/reg64; reg18=reg18/reg64; reg86=reg86-reg23; reg20=reg55+reg20; reg69=reg69-reg4;
    reg9=reg54*vectors[0][indices[0]+0]; reg45=reg52+reg45; reg36=reg48*reg36; reg82=reg63+reg82; reg29=reg40*reg94;
    reg7=reg11-reg7; reg21=reg60+reg21; reg11=reg125*reg98; reg34=var_inter[0]*vectors[0][indices[1]+1]; reg37=reg24+reg37;
    reg24=reg90*vectors[0][indices[0]+1]; reg80=reg0*reg80; reg59=reg35+reg59; reg35=var_inter[0]*vectors[0][indices[1]+0]; reg44=reg90*vectors[0][indices[0]+0];
    reg46=reg54*vectors[0][indices[1]+1]; reg48=reg28*reg58; reg52=reg27*reg58; reg55=reg75*reg38; reg29=reg11-reg29;
    reg11=var_inter[0]*vectors[0][indices[2]+0]; reg56=reg81+reg56; reg55=reg52-reg55; reg80=reg59+reg80; reg18=reg18*reg82;
    reg35=reg44+reg35; reg22=reg22*reg37; reg83=reg83*reg82; reg91=reg91*reg37; reg36=reg45+reg36;
    reg44=reg6*reg27; reg45=var_inter[1]*vectors[0][indices[2]+0]; reg52=reg28*reg75; reg26=reg48-reg26; reg48=reg21*reg86;
    reg59=reg20*reg69; reg60=var_inter[0]*vectors[0][indices[2]+1]; reg34=reg24+reg34; reg24=var_inter[1]*vectors[0][indices[2]+1]; reg12=reg12/reg64;
    reg1=reg46-reg1; reg14=reg14/reg64; reg33=reg33/reg64; reg41=reg2+reg41; reg7=reg7/reg64;
    reg31=reg0*reg31; reg9=reg15-reg9; reg35=reg11-reg35; reg59=reg48-reg59; reg44=reg52-reg44;
    reg45=reg9+reg45; reg0=reg55*reg56; reg2=reg90*vectors[0][indices[3]+0]; reg33=reg33*reg36; reg83=reg91-reg83;
    reg7=reg7*reg36; reg34=reg60-reg34; reg22=reg18-reg22; reg37=reg12*reg37; reg82=reg14*reg82;
    reg9=var_inter[1]*vectors[0][indices[3]+0]; reg31=reg41+reg31; reg64=reg29/reg64; reg24=reg1+reg24; reg1=reg90*vectors[0][indices[3]+1];
    reg11=var_inter[1]*vectors[0][indices[3]+1]; reg12=reg80*reg26; reg2=reg35+reg2; reg14=reg28*reg31; reg83=reg33+reg83;
    reg11=reg24-reg11; reg15=reg27*reg31; reg18=reg80*reg38; reg7=reg22-reg7; reg22=reg31*reg44;
    reg36=reg64*reg36; reg69=reg69/reg59; reg38=reg56*reg38; reg20=reg20/reg59; reg9=reg45-reg9;
    reg86=reg86/reg59; reg21=reg21/reg59; reg1=reg34+reg1; reg24=1-(*f.m).resolution; reg82=reg37-reg82;
    reg12=reg0-reg12; reg28=reg28*reg80; reg82=reg36+reg82; reg27=reg27*reg56; reg0=reg86*reg2;
    reg29=(*f.m).resolution*reg102; reg33=reg20*reg9; reg14=reg38-reg14; reg34=reg11*reg21; reg22=reg12+reg22;
    reg12=(*f.m).resolution*reg19; reg35=reg6*reg31; reg36=reg80*reg58; reg83=reg24*reg83; reg31=reg75*reg31;
    reg58=reg56*reg58; reg15=reg18-reg15; reg18=reg69*reg1; reg7=reg24*reg7; reg31=reg36-reg31;
    reg11=reg11*reg20; reg15=reg15/reg22; reg18=reg34-reg18; reg33=reg0-reg33; reg9=reg21*reg9;
    reg1=reg86*reg1; reg2=reg69*reg2; reg82=reg24*reg82; reg35=reg58-reg35; reg28=reg27-reg28;
    reg14=reg14/reg22; reg0=(*f.m).resolution*reg96; reg83=reg29+reg83; reg7=reg12+reg7; reg80=reg6*reg80;
    reg56=reg75*reg56; reg2=reg9-reg2; reg83=reg83*(*f.m).deltaT; reg82=reg0+reg82; reg88=reg88*reg24;
    reg85=reg85*reg24; reg0=(*f.m).resolution*reg15; reg55=reg55/reg22; reg6=(*f.m).resolution*reg14; reg31=reg31/reg22;
    reg26=reg26/reg22; reg7=reg7*(*f.m).deltaT; reg35=reg35/reg22; reg18=reg33+reg18; reg28=reg28/reg22;
    reg80=reg56-reg80; reg11=reg1-reg11; reg6=reg85-reg6; reg0=reg88+reg0; reg82=reg82*(*f.m).deltaT;
    reg1=(*f.m).resolution*reg28; reg9=(*f.m).resolution*reg55; reg12=(*f.m).resolution*reg26; reg27=(*f.m).resolution*reg35; reg29=(*f.m).resolution*reg31;
    reg115=reg115*reg24; reg11=reg11-reg7; reg80=reg80/reg22; reg125=reg125*reg24; reg2=reg2-reg83;
    reg40=reg40*reg24; reg22=reg44/reg22; reg18=0.5*reg18; reg94=reg24*reg94; reg98=reg24*reg98;
    reg18=reg18-reg82; reg1=reg115+reg1; reg33=reg2*reg0; reg34=reg11*reg6; reg27=reg125+reg27;
    reg29=reg40-reg29; reg12=reg94-reg12; reg98=reg9+reg98; reg9=(*f.m).resolution*reg22; reg17=reg24*reg17;
    reg36=(*f.m).resolution*reg80; reg24=reg73*reg24; reg37=reg18*reg1; reg34=reg33+reg34; reg33=reg2*reg98;
    reg38=reg90*reg69; reg40=reg11*reg27; reg2=reg2*reg29; reg41=var_inter[0]*reg86; reg44=reg90*reg86;
    reg17=reg9+reg17; reg11=reg11*reg12; reg9=var_inter[1]*reg21; reg45=var_inter[1]*reg20; reg36=reg24-reg36;
    reg24=reg54*reg20; reg46=var_inter[0]*reg69; reg48=reg54*reg21; reg40=reg2+reg40; reg2=reg41-reg45;
    reg52=reg9-reg46; reg56=reg38-reg48; reg58=reg18*reg36; reg60=reg44+reg45; reg61=reg24+reg41;
    reg33=reg11+reg33; reg11=reg48+reg46; reg63=reg24-reg44; reg37=reg34+reg37; reg34=reg38+reg9;
    reg18=reg18*reg17; reg58=reg40+reg58; reg40=0.5*reg34; reg64=0.5*reg60; reg65=0.5*reg63;
    reg66=0.5*reg52; reg18=reg33+reg18; reg33=0.5*reg2; reg70=0.5*reg56; reg73=0.5*reg61;
    reg74=0.5*reg11; reg37=2*reg37; reg75=reg56*reg18; reg76=reg37*reg74; reg77=reg61*reg58;
    reg78=reg37*reg33; reg79=reg52*reg18; reg81=reg37*reg66; reg84=reg2*reg58; reg85=reg90*var_inter[1];
    reg87=var_inter[0]*var_inter[1]; reg88=reg54*var_inter[0]; reg89=reg54*reg90; reg91=reg37*reg65; reg94=reg37*reg70;
    reg95=reg63*reg58; reg97=reg34*reg18; reg99=reg37*reg64; reg100=reg60*reg58; reg101=reg37*reg40;
    reg103=reg11*reg18; reg104=reg37*reg73; reg95=reg94+reg95; reg94=reg89*elem.f_vol_e[1]; reg79=reg78+reg79;
    reg78=reg87*elem.f_vol_e[0]; reg103=reg103-reg104; reg105=reg88*elem.f_vol_e[0]; reg84=reg81+reg84; reg81=reg87*elem.f_vol_e[1];
    reg106=reg85*elem.f_vol_e[0]; reg99=reg99-reg97; reg75=reg91+reg75; reg91=reg89*elem.f_vol_e[0]; reg107=reg85*elem.f_vol_e[1];
    reg100=reg100-reg101; reg108=reg88*elem.f_vol_e[1]; reg76=reg76-reg77; reg100=reg100-reg107; reg103=reg103-reg105;
    reg84=reg84-reg81; reg99=reg99-reg106; reg79=reg79-reg78; reg75=reg75-reg91; reg76=reg76-reg108;
    reg95=reg95-reg94; reg75=reg75*reg59; reg100=reg59*reg100; reg103=reg103*reg59; reg84=reg59*reg84;
    reg76=reg59*reg76; reg99=reg59*reg99; reg95=reg95*reg59; reg79=reg59*reg79; T vec_6=ponderation*reg99;
    T vec_3=ponderation*reg76; T vec_2=ponderation*reg103; T vec_5=ponderation*reg84; T vec_4=ponderation*reg79; T vec_1=ponderation*reg95;
    T vec_7=ponderation*reg100; T vec_0=ponderation*reg75;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
sollicitation[indices[0]+0] += vec_0;
sollicitation[indices[0]+1] += vec_1;
sollicitation[indices[1]+0] += vec_2;
sollicitation[indices[1]+1] += vec_3;
sollicitation[indices[2]+0] += vec_4;
sollicitation[indices[2]+1] += vec_5;
sollicitation[indices[3]+0] += vec_6;
sollicitation[indices[3]+1] += vec_7;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
#ifndef ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_false_true_true
#define ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_false_true_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const typename TM::TNode &node,
      const unsigned *indices){ 
  #define PNODE(N) node
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
#endif
// 
#ifndef ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_true_true_true
#define ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_true_true_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const typename TM::TNode &node,
      const unsigned *indices){ 
  #define PNODE(N) node
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
#endif
// 
#ifndef ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_false_true_false
#define ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_false_true_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const typename TM::TNode &node,
      const unsigned *indices){ 
  #define PNODE(N) node
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
#endif
// 
#ifndef ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_true_true_false
#define ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_true_true_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const typename TM::TNode &node,
      const unsigned *indices){ 
  #define PNODE(N) node
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
#endif
// 
#ifndef ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_symmetric_version_false_true
#define ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_symmetric_version_false_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE ,unsigned symmetric_version>
void add_nodal_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const typename TM::TNode &node,
      const unsigned *indices){ 
  #define PNODE(N) node
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
#endif
// 
#ifndef ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_symmetric_version_false_false
#define ADD_NODAL_MATRIX_elasticity_orthotropy_stat_Qstat_symmetric_version_false_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE ,unsigned symmetric_version>
void add_nodal_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const typename TM::TNode &node,
      const unsigned *indices){ 
  #define PNODE(N) node
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
#endif
// 
#ifndef ADD_NODAL_RESIDUAL_elasticity_orthotropy_stat_Qstat
#define ADD_NODAL_RESIDUAL_elasticity_orthotropy_stat_Qstat
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE>
void add_nodal_residual(
      const Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const typename TM::TNode &node,
      const unsigned *indices ) { 
  #define PNODE(N) node
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
#endif
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const Element<Bar,DefaultBehavior,Node<2,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}

#ifndef elasticity_orthotropy_stat_Qstat_read_material_to_mesh
#define elasticity_orthotropy_stat_Qstat_read_material_to_mesh
template<class TM, class T, bool wont_add_nz>
void read_material_to_mesh_(const XmlNode &n, Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f){ 
    if(n.has_attribute("elastic_modulus_1"))  
        n.get_attribute("elastic_modulus_1", f.m->elastic_modulus_1 ); 
    else  
        std::cerr << "Warning using default value of elastic_modulus_1 : " << f.m->elastic_modulus_1 << std::endl; 

    if(n.has_attribute("alpha_2"))  
        n.get_attribute("alpha_2", f.m->alpha_2 ); 
    else  
        std::cerr << "Warning using default value of alpha_2 : " << f.m->alpha_2 << std::endl; 

    if(n.has_attribute("alpha_3"))  
        n.get_attribute("alpha_3", f.m->alpha_3 ); 
    else  
        std::cerr << "Warning using default value of alpha_3 : " << f.m->alpha_3 << std::endl; 

    if(n.has_attribute("alpha_1"))  
        n.get_attribute("alpha_1", f.m->alpha_1 ); 
    else  
        std::cerr << "Warning using default value of alpha_1 : " << f.m->alpha_1 << std::endl; 

    if(n.has_attribute("f_vol"))  
        n.get_attribute("f_vol", f.m->f_vol ); 
    else  
        std::cerr << "Warning using default value of f_vol : " << f.m->f_vol << std::endl; 

    if(n.has_attribute("shear_modulus_13"))  
        n.get_attribute("shear_modulus_13", f.m->shear_modulus_13 ); 
    else  
        std::cerr << "Warning using default value of shear_modulus_13 : " << f.m->shear_modulus_13 << std::endl; 

    if(n.has_attribute("shear_modulus_12"))  
        n.get_attribute("shear_modulus_12", f.m->shear_modulus_12 ); 
    else  
        std::cerr << "Warning using default value of shear_modulus_12 : " << f.m->shear_modulus_12 << std::endl; 

    if(n.has_attribute("deltaT"))  
        n.get_attribute("deltaT", f.m->deltaT ); 
    else  
        std::cerr << "Warning using default value of deltaT : " << f.m->deltaT << std::endl; 

    if(n.has_attribute("density"))  
        n.get_attribute("density", f.m->density ); 
    else  
        std::cerr << "Warning using default value of density : " << f.m->density << std::endl; 

    if(n.has_attribute("poisson_ratio_13"))  
        n.get_attribute("poisson_ratio_13", f.m->poisson_ratio_13 ); 
    else  
        std::cerr << "Warning using default value of poisson_ratio_13 : " << f.m->poisson_ratio_13 << std::endl; 

    if(n.has_attribute("poisson_ratio_12"))  
        n.get_attribute("poisson_ratio_12", f.m->poisson_ratio_12 ); 
    else  
        std::cerr << "Warning using default value of poisson_ratio_12 : " << f.m->poisson_ratio_12 << std::endl; 

    if(n.has_attribute("v2"))  
        n.get_attribute("v2", f.m->v2 ); 
    else  
        std::cerr << "Warning using default value of v2 : " << f.m->v2 << std::endl; 

    if(n.has_attribute("elastic_modulus_3"))  
        n.get_attribute("elastic_modulus_3", f.m->elastic_modulus_3 ); 
    else  
        std::cerr << "Warning using default value of elastic_modulus_3 : " << f.m->elastic_modulus_3 << std::endl; 

    if(n.has_attribute("elastic_modulus_2"))  
        n.get_attribute("elastic_modulus_2", f.m->elastic_modulus_2 ); 
    else  
        std::cerr << "Warning using default value of elastic_modulus_2 : " << f.m->elastic_modulus_2 << std::endl; 

    if(n.has_attribute("shear_modulus_23"))  
        n.get_attribute("shear_modulus_23", f.m->shear_modulus_23 ); 
    else  
        std::cerr << "Warning using default value of shear_modulus_23 : " << f.m->shear_modulus_23 << std::endl; 

    if(n.has_attribute("v1"))  
        n.get_attribute("v1", f.m->v1 ); 
    else  
        std::cerr << "Warning using default value of v1 : " << f.m->v1 << std::endl; 

    if(n.has_attribute("resolution"))  
        n.get_attribute("resolution", f.m->resolution ); 
    else  
        std::cerr << "Warning using default value of resolution : " << f.m->resolution << std::endl; 

    if(n.has_attribute("poisson_ratio_23"))  
        n.get_attribute("poisson_ratio_23", f.m->poisson_ratio_23 ); 
    else  
        std::cerr << "Warning using default value of poisson_ratio_23 : " << f.m->poisson_ratio_23 << std::endl; 

  };
#endif // elasticity_orthotropy_stat_Qstat_read_material_to_mesh
} // namespace LMT

