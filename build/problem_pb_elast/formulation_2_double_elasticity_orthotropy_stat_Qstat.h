
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
  static const unsigned nb_matrices = 1;
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
    T reg0=vecs[1][indice+0]-vecs[0][indice+0]; T reg1=vecs[1][indice+1]-vecs[0][indice+1]; reg0=abs(reg0); reg1=abs(reg1); return max(reg1,reg0);
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; reg1=1.0/reg1; T reg2=2*(*f.m).shear_modulus_12;
    T reg3=reg1*reg0; reg2=1.0/reg2; T reg4=reg2*reg3; T reg5=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg6=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg7=pow((*f.m).v1[0],2); T reg8=pow((*f.m).v1[1],2); T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=pow((*f.m).v2[0],2); T reg11=pow((*f.m).v2[1],2);
    T reg12=pow((*f.m).v1[2],2); T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg14=1.0/(*f.m).elastic_modulus_2; reg8=reg7+reg8; reg7=reg5*reg4;
    T reg15=pow((*f.m).v2[2],2); reg11=reg10+reg11; reg10=reg9*reg4; T reg16=reg6*reg4; T reg17=reg6*reg7;
    T reg18=reg13*reg10; T reg19=reg6*reg16; T reg20=reg14*reg10; reg15=reg11+reg15; reg12=reg8+reg12;
    reg12=pow(reg12,0.5); reg8=reg14*reg7; reg11=reg13*reg16; reg17=reg18+reg17; reg15=pow(reg15,0.5);
    reg19=reg20-reg19; reg20=1.0/(*f.m).elastic_modulus_1; T reg21=(*f.m).v2[2]/reg15; T reg22=(*f.m).v2[1]/reg15; T reg23=(*f.m).v1[1]/reg12;
    T reg24=reg20*reg19; T reg25=(*f.m).v1[2]/reg12; T reg26=reg11+reg8; T reg27=reg13*reg17; reg10=reg20*reg10;
    T reg28=reg6*reg3; T reg29=reg5*reg7; T reg30=reg13*reg4; T reg31=reg25*reg22; T reg32=reg5*reg26;
    T reg33=reg23*reg21; reg12=(*f.m).v1[0]/reg12; reg27=reg24-reg27; reg15=(*f.m).v2[0]/reg15; reg24=reg9*reg3;
    reg3=reg5*reg3; reg4=reg14*reg4; T reg34=reg2*reg0; T reg35=reg5*reg16; T reg36=reg34*reg5;
    T reg37=reg14*reg24; reg16=reg20*reg16; reg29=reg10-reg29; reg7=reg13*reg7; reg10=reg34*reg6;
    T reg38=reg5*reg30; T reg39=reg2*reg1; T reg40=2*reg15; reg34=reg34*reg9; T reg41=2*reg12;
    reg35=reg18+reg35; reg18=reg3*reg6; T reg42=reg33-reg31; T reg43=reg6*reg28; T reg44=reg25*reg15;
    T reg45=reg12*reg21; reg24=reg13*reg24; T reg46=reg5*reg4; reg32=reg27-reg32; reg27=reg36*reg6;
    T reg47=reg39*reg6; T reg48=reg10*reg6; T reg49=reg34*reg13; T reg50=pow(reg22,2); T reg51=reg39*reg9;
    reg34=reg34*reg14; T reg52=reg22*reg40; T reg53=pow(reg15,2); reg18=reg24+reg18; reg43=reg37-reg43;
    reg4=reg20*reg4; reg39=reg39*reg5; reg7=reg16+reg7; reg30=reg13*reg30; reg19=reg19/reg32;
    reg35=reg35/reg32; reg46=reg11+reg46; reg24=2*reg42; reg37=reg23*reg41; T reg54=pow(reg23,2);
    T reg55=pow(reg12,2); T reg56=reg44-reg45; reg17=reg17/reg32; reg29=reg29/reg32; reg28=reg13*reg28;
    reg38=reg16+reg38; reg16=reg12*reg22; reg3=reg3*reg14; T reg57=reg23*reg15; T reg58=pow(reg25,2);
    reg36=reg36*reg14; reg7=reg7/reg32; T reg59=pow(reg56,2); reg48=reg34-reg48; reg34=reg50*reg17;
    T reg60=pow(reg42,2); T reg61=pow(reg21,2); reg10=reg10*reg13; T reg62=reg54*reg19; T reg63=reg55*reg35;
    reg27=reg49+reg27; reg46=reg46/reg32; reg49=reg51*reg14; T reg64=reg53*reg29; reg51=reg51*reg13;
    T reg65=reg47*reg6; T reg66=reg50*reg29; T reg67=reg54*reg35; T reg68=reg39*reg6; T reg69=reg16-reg57;
    T reg70=reg37*reg35; T reg71=reg53*reg17; reg30=reg4-reg30; reg26=reg26/reg32; reg4=reg37*reg19;
    T reg72=reg52*reg17; reg43=reg43*reg20; T reg73=reg55*reg19; reg38=reg38/reg32; T reg74=reg28+reg3;
    T reg75=reg56*reg24; reg18=reg18*reg13; T reg76=reg52*reg29; T reg77=reg52*reg38; T reg78=reg75*reg7;
    reg76=reg70+reg76; reg70=reg58*reg35; T reg79=reg61*reg29; T reg80=reg59*reg7; T reg81=reg50*reg38;
    T reg82=reg54*reg46; T reg83=reg37*reg46; T reg84=reg53*reg38; T reg85=reg55*reg46; reg68=reg51+reg68;
    reg66=reg67+reg66; reg34=reg62+reg34; reg51=reg59*reg26; reg62=pow(reg69,2); reg67=reg58*reg19;
    T reg86=reg61*reg17; reg47=reg47*reg13; reg73=reg71+reg73; reg30=reg30/reg32; reg71=reg60*reg26;
    reg39=reg39*reg14; reg72=reg4+reg72; reg4=reg75*reg26; reg18=reg43-reg18; reg74=reg74*reg5;
    reg43=reg60*reg7; reg65=reg49-reg65; reg64=reg63+reg64; reg49=reg10+reg36; reg27=reg27*reg13;
    reg48=reg48*reg20; reg74=reg18-reg74; reg51=reg34+reg51; reg65=reg65*reg20; reg68=reg68*reg13;
    reg18=reg25*reg41; reg34=2*reg23; reg63=reg21*reg40; T reg87=reg47+reg39; reg84=reg85+reg84;
    reg85=reg60*reg30; T reg88=2*reg22; reg81=reg82+reg81; reg82=reg59*reg30; T reg89=reg58*reg46;
    reg71=reg73+reg71; reg49=reg49*reg5; reg73=reg61*reg38; reg27=reg48-reg27; reg79=reg70+reg79;
    reg80=reg66+reg80; reg48=reg62*reg7; reg66=reg75*reg30; reg77=reg83+reg77; reg43=reg64+reg43;
    reg64=reg12*reg23; reg70=reg15*reg22; reg4=reg72+reg4; reg78=reg76+reg78; reg72=reg62*reg26;
    reg86=reg67+reg86; reg73=reg89+reg73; reg67=reg62*reg30; reg76=reg54*reg51; reg24=reg69*reg24;
    reg83=2*reg56; reg89=reg80*reg50; T reg90=reg12*reg15; T reg91=reg21*reg88; reg66=reg77+reg66;
    reg77=reg43*reg50; T reg92=reg23*reg22; T reg93=reg51*reg64; reg87=reg87*reg5; T reg94=reg80*reg70;
    reg68=reg65-reg68; reg49=reg27-reg49; reg72=reg86+reg72; reg27=reg56*reg42; reg65=reg54*reg71;
    reg86=reg12*reg56; T reg95=reg23*reg42; T reg96=reg63*reg17; T reg97=reg18*reg19; reg16=reg57+reg16;
    reg74=reg74/reg32; reg57=reg78*reg70; T reg98=reg4*reg64; T reg99=reg18*reg35; T reg100=reg25*reg34;
    reg48=reg79+reg48; reg79=reg63*reg29; T reg101=reg55*reg71; T reg102=reg43*reg53; T reg103=reg55*reg51;
    T reg104=reg80*reg53; T reg105=reg78*reg50; T reg106=reg4*reg54; T reg107=reg4*reg55; reg85=reg84+reg85;
    reg82=reg81+reg82; reg81=reg78*reg53; reg17=reg91*reg17; reg29=reg91*reg29; reg19=reg100*reg19;
    reg84=reg23*reg56; T reg108=reg24*reg26; reg96=reg97+reg96; reg35=reg100*reg35; reg97=reg66*reg59;
    T reg109=reg12*reg42; T reg110=reg24*reg7; reg79=reg99+reg79; reg99=reg85*reg60; reg102=reg101+reg102;
    reg57=reg98+reg57; reg67=reg73+reg67; reg73=reg66*reg27; reg98=reg48*reg50; reg101=reg92*reg74;
    T reg111=reg25*reg21; T reg112=reg72*reg54; T reg113=reg90*reg74; T reg114=reg82*reg59; reg89=reg76+reg89;
    reg76=reg66*reg60; T reg115=reg85*reg59; reg77=reg65+reg77; reg65=reg63*reg38; T reg116=reg18*reg46;
    T reg117=reg48*reg53; T reg118=reg72*reg55; T reg119=reg43*reg70; reg87=reg68-reg87; reg68=reg22*reg42;
    T reg120=reg71*reg64; T reg121=reg15*reg56; reg83=reg69*reg83; reg94=reg93+reg94; reg93=reg15*reg41;
    T reg122=reg82*reg27; T reg123=reg22*reg34; reg49=reg49/reg32; T reg124=2*reg25; T reg125=reg16*reg74;
    reg81=reg107+reg81; reg105=reg106+reg105; reg106=reg82*reg60; reg104=reg103+reg104; reg86=reg95+reg86;
    reg95=reg85*reg27; reg103=reg72*reg64; reg119=reg120+reg119; reg107=reg48*reg70; reg120=reg125*reg93;
    T reg126=reg125*reg16; T reg127=reg25*reg69; reg73=reg57+reg73; reg76=reg81+reg76; reg57=reg101*reg16;
    reg81=reg125*reg123; reg97=reg105+reg97; reg122=reg94+reg122; reg34=reg56*reg34; reg41=reg42*reg41;
    reg32=reg87/reg32; reg121=reg68+reg121; reg68=reg22*reg56; reg65=reg116+reg65; reg87=reg24*reg30;
    reg46=reg100*reg46; reg38=reg91*reg38; reg94=reg111*reg74; reg105=reg15*reg42; reg116=reg67*reg60;
    reg117=reg118+reg117; reg118=reg109*reg49; T reg128=reg84*reg49; reg7=reg83*reg7; reg29=reg35+reg29;
    reg35=reg101*reg93; reg106=reg104+reg106; reg110=reg79+reg110; reg79=reg113*reg93; reg99=reg102+reg99;
    reg102=reg86*reg49; reg115=reg77+reg115; reg77=reg113*reg123; reg114=reg89+reg114; reg89=reg101*reg123;
    reg98=reg112+reg98; reg104=reg67*reg59; reg112=reg53*reg20; T reg129=reg50*reg13; T reg130=reg93*reg20;
    T reg131=reg123*reg13; T reg132=reg21*reg124; T reg133=reg123*reg14; reg26=reg83*reg26; T reg134=reg93*reg13;
    reg17=reg19+reg17; reg19=reg50*reg14; T reg135=reg53*reg13; reg108=reg96+reg108; reg96=reg110*reg50;
    T reg136=reg108*reg54; reg7=reg29+reg7; reg44=reg45+reg44; reg26=reg17+reg26; reg17=reg25*reg42;
    reg29=reg12*reg69; reg81=reg97+reg81; reg45=reg102*reg34; reg97=reg61*reg6; reg135=reg19-reg135;
    reg19=reg132*reg6; reg88=reg56*reg88; reg40=reg42*reg40; T reg137=reg21*reg69; reg87=reg65+reg87;
    reg38=reg46+reg38; reg30=reg83*reg30; reg46=reg127*reg49; reg65=reg32*reg105; T reg138=reg32*reg68;
    T reg139=reg32*reg121; reg79=reg99+reg79; reg99=reg118*reg41; reg35=reg106+reg35; reg106=reg128*reg41;
    reg116=reg117+reg116; reg117=(*f.m).alpha_1*reg55; T reg140=reg50*(*f.m).alpha_2; T reg141=(*f.m).alpha_1*reg54; T reg142=reg53*(*f.m).alpha_2;
    T reg143=reg94*reg93; reg120=reg76+reg120; reg76=reg102*reg41; T reg144=reg108*reg55; T reg145=reg110*reg53;
    reg77=reg115+reg77; reg115=reg118*reg34; reg89=reg114+reg89; reg114=reg128*reg34; reg104=reg98+reg104;
    reg98=reg94*reg123; T reg146=reg102*reg86; reg126=reg73+reg126; reg73=reg61*reg5; reg129=reg112-reg129;
    reg112=reg132*reg5; reg131=reg130-reg131; reg130=reg53*reg5; T reg147=reg67*reg27; T reg148=reg123*reg6;
    reg95=reg119+reg95; reg119=reg113*reg16; reg134=reg133-reg134; reg133=reg93*reg5; T reg149=reg50*reg6;
    reg57=reg122+reg57; reg122=reg128*reg86; reg107=reg103+reg107; reg103=reg108*reg64; T reg150=reg65*reg88;
    reg115=reg77+reg115; reg77=reg110*reg70; T reg151=reg7*reg53; T reg152=reg26*reg55; reg122=reg57+reg122;
    reg57=reg138*reg121; reg114=reg89+reg114; reg89=reg138*reg88; reg98=reg104+reg98; reg104=reg46*reg34;
    T reg153=reg139*reg121; T reg154=reg54*reg13; reg146=reg126+reg146; reg147=reg107+reg147; reg107=(*f.m).alpha_3*reg59;
    reg140=reg141+reg140; reg126=(*f.m).alpha_3*reg60; reg142=reg117+reg142; reg119=reg95+reg119; reg45=reg81+reg45;
    reg81=reg139*reg88; reg95=reg118*reg86; reg116=reg143+reg116; reg117=reg46*reg41; reg96=reg136+reg96;
    reg136=reg87*reg59; reg141=reg26*reg54; reg143=reg7*reg50; reg76=reg120+reg76; reg120=reg139*reg40;
    T reg155=reg25*reg56; reg29=reg17+reg29; reg145=reg144+reg145; reg17=reg87*reg60; reg148=reg133+reg148;
    reg133=reg132*reg9; reg144=reg138*reg40; reg106=reg35+reg106; reg35=reg65*reg40; reg99=reg79+reg99;
    reg149=reg130+reg149; reg79=reg61*reg9; reg130=reg55*reg20; T reg156=reg32*reg137; T reg157=reg23*reg69;
    T reg158=reg44*reg74; reg30=reg38+reg30; reg19=reg134-reg19; reg38=reg21*reg42; reg134=reg15*reg69;
    T reg159=reg94*reg16; reg33=reg31+reg33; reg73=reg129-reg73; reg31=(*f.m).alpha_1*reg58; reg112=reg131-reg112;
    reg97=reg135-reg97; reg129=reg55*reg13; reg131=reg54*reg14; reg135=reg61*(*f.m).alpha_2; reg153=reg146+reg153;
    reg146=reg158*reg123; reg136=reg96+reg136; reg149=reg79-reg149; reg57=reg122+reg57; reg81=reg45+reg81;
    reg45=reg46*reg86; reg79=reg97*reg54; reg148=reg133-reg148; reg159=reg147+reg159; reg126=reg142+reg126;
    reg96=reg73*reg55; reg122=reg97*reg92; reg133=reg112*reg90; reg142=reg22*reg69; reg147=reg21*reg56;
    reg134=reg38+reg134; reg77=reg103+reg77; reg38=reg87*reg27; reg103=reg26*reg64; T reg160=reg7*reg70;
    reg12=reg12*reg25; reg15=reg15*reg21; T reg161=reg92*reg19; T reg162=reg73*reg90; T reg163=reg55*reg5;
    T reg164=reg54*reg6; T reg165=reg65*reg121; reg95=reg119+reg95; reg119=reg30*reg59; reg143=reg141+reg143;
    reg157=reg155+reg157; reg154=reg130-reg154; reg130=reg58*reg5; reg104=reg98+reg104; reg117=reg116+reg117;
    reg98=reg50*reg19; reg116=reg112*reg53; reg141=reg156*reg40; reg89=reg114+reg89; reg150=reg115+reg150;
    reg120=reg76+reg120; reg17=reg145+reg17; reg76=reg30*reg60; reg151=reg152+reg151; reg114=reg158*reg93;
    reg107=reg140+reg107; reg135=reg31+reg135; reg31=(*f.m).alpha_3*reg62; reg115=(*f.m).alpha_1*reg64; reg140=reg70*(*f.m).alpha_2;
    reg145=reg58*reg6; reg129=reg131-reg129; reg74=reg33*reg74; reg131=reg29*reg49; reg152=reg73*reg53;
    reg155=reg97*reg50; reg35=reg99+reg35; reg144=reg106+reg144; reg99=reg54*reg19; reg106=reg112*reg55;
    T reg166=reg156*reg88; reg99=reg106+reg99; reg106=(*f.m).alpha_3*reg27; T reg167=reg35*reg126; T reg168=(*f.m).alpha_1*reg12;
    reg140=reg115+reg140; reg115=reg144*reg107; T reg169=reg15*(*f.m).alpha_2; T reg170=reg156*reg121; reg45=reg159+reg45;
    reg165=reg95+reg165; reg95=reg58*reg148; reg141=reg117+reg141; reg155=reg152+reg155; reg117=reg89*reg153;
    reg152=reg61*reg149; reg159=reg120*reg57; T reg171=reg81*reg57; T reg172=reg144*reg153; reg42=reg69*reg42;
    reg25=reg23*reg25; reg22=reg22*reg21; reg142=reg147+reg142; reg145=reg129-reg145; reg49=reg49*reg157;
    reg23=reg32*reg134; reg79=reg96+reg79; reg96=reg16*reg2; reg70=reg70*reg2; reg130=reg154-reg130;
    reg164=reg163+reg164; reg129=reg58*reg9; reg147=reg111*reg148; reg161=reg133+reg161; reg133=reg111*reg149;
    reg122=reg162+reg122; reg38=reg77+reg38; reg77=reg158*reg16; reg160=reg103+reg160; reg103=reg30*reg27;
    reg93=reg74*reg93; reg76=reg151+reg76; reg151=reg131*reg41; reg114=reg17+reg114; reg31=reg135+reg31;
    reg17=reg150*reg126; reg135=reg89*reg107; reg154=reg58*reg149; reg123=reg74*reg123; reg119=reg143+reg119;
    reg98=reg116+reg98; reg116=reg61*reg148; reg166=reg104+reg166; reg104=reg131*reg34; reg146=reg136+reg146;
    reg136=reg107*reg57; reg143=reg165*reg126; reg106=reg140+reg106; reg169=reg168+reg169; reg140=(*f.m).alpha_3*reg42;
    reg162=(*f.m).alpha_1*reg25; reg163=reg22*(*f.m).alpha_2; reg170=reg45+reg170; reg79=reg154+reg79; reg171=reg117-reg171;
    reg159=reg172-reg159; reg45=reg144*reg81; reg117=reg120*reg89; reg154=reg52*reg70; reg152=reg155+reg152;
    reg166=reg31*reg166; reg135=reg17+reg135; reg56=reg69*reg56; reg17=reg145*reg50; reg155=reg130*reg53;
    reg168=reg37*reg96; reg95=reg99+reg95; reg99=reg52*reg96; reg172=reg145*reg54; T reg173=reg130*reg55;
    T reg174=reg44*reg1; T reg175=reg15*reg1; reg116=reg98+reg116; reg98=reg23*reg88; reg32=reg32*reg142;
    reg104=reg146+reg104; reg164=reg129-reg164; reg123=reg119+reg123; reg119=reg49*reg34; reg129=reg16*reg96;
    reg147=reg161+reg147; reg146=reg16*reg70; reg133=reg122+reg133; reg77=reg38+reg77; reg38=reg131*reg86;
    reg122=reg49*reg41; reg93=reg76+reg93; reg103=reg160+reg103; reg76=reg74*reg16; reg160=reg23*reg40;
    reg151=reg114+reg151; reg115=reg167+reg115; reg141=reg141*reg31; reg114=reg37*reg70; reg161=reg130*reg90;
    reg167=reg145*reg92; reg146=reg133+reg146; reg133=reg44*reg175; T reg176=reg32*reg88; reg129=reg147+reg129;
    reg147=reg18*reg175; T reg177=reg44*reg174; reg119=reg123+reg119; reg64=reg64*reg2; reg123=reg35*reg171;
    T reg178=reg150*reg159; T reg179=reg106*reg81; T reg180=reg63*reg175; reg154=reg152+reg154; reg117=reg45-reg117;
    reg166=reg135+reg166; reg45=reg22*reg0; reg160=reg151+reg160; reg114=reg79+reg114; reg79=reg49*reg86;
    reg76=reg103+reg76; reg170=reg31*reg170; reg136=reg143+reg136; reg103=reg23*reg121; reg140=reg169+reg140;
    reg38=reg77+reg38; reg163=reg162+reg163; reg77=(*f.m).alpha_3*reg56; reg122=reg93+reg122; reg93=reg32*reg40;
    reg135=reg120*reg106; reg168=reg95+reg168; reg95=reg18*reg174; reg141=reg115+reg141; reg115=reg63*reg174;
    reg99=reg116+reg99; reg98=reg104+reg98; reg17=reg155+reg17; reg104=reg61*reg164; reg116=reg33*reg0;
    reg172=reg173+reg172; reg143=reg58*reg164; reg151=reg100*reg116; reg176=reg119+reg176; reg119=reg33*reg45;
    reg133=reg146+reg133; reg95=reg168+reg95; reg146=reg120*reg165; reg152=reg111*reg164; reg167=reg161+reg167;
    reg155=elem.pos(1)[0]-elem.pos(0)[0]; reg93=reg122+reg93; reg135=reg141+reg135; reg77=reg163+reg77; reg115=reg99+reg115;
    reg99=reg91*reg116; reg122=reg150*reg153; reg103=reg38+reg103; reg38=reg81*reg165; reg160=reg160*reg140;
    reg170=reg136+reg170; reg136=reg106*reg153; reg141=elem.pos(2)[1]-elem.pos(0)[1]; reg161=elem.pos(2)[0]-elem.pos(0)[0]; reg162=reg35*reg153;
    reg79=reg76+reg79; reg76=reg32*reg121; reg163=reg91*reg45; reg98=reg140*reg98; reg179=reg166+reg179;
    reg104=reg17+reg104; reg17=reg52*reg64; reg147=reg114+reg147; reg114=elem.pos(1)[1]-elem.pos(0)[1]; reg180=reg154+reg180;
    reg154=reg100*reg45; reg178=reg123-reg178; reg123=reg117*reg165; reg166=reg37*reg64; reg168=reg12*reg1;
    reg169=reg33*reg116; reg177=reg129+reg177; reg143=reg172+reg143; reg93=reg93*reg77; reg160=reg135+reg160;
    reg99=reg115+reg99; reg115=reg150*reg57; reg129=reg120*reg150; reg17=reg104+reg17; reg123=reg178+reg123;
    reg76=reg79+reg76; reg79=reg89*reg165; reg104=reg25*reg0; reg154=reg147+reg154; reg136=reg170+reg136;
    reg103=reg140*reg103; reg38=reg122-reg38; reg151=reg95+reg151; reg146=reg162-reg146; reg119=reg133+reg119;
    reg169=reg177+reg169; reg166=reg143+reg166; reg176=reg77*reg176; reg95=reg35*reg57; reg122=reg16*reg64;
    reg152=reg167+reg152; reg133=reg18*reg168; reg135=reg114*reg161; reg163=reg180+reg163; reg143=reg141*reg155;
    reg147=reg144*reg165; reg162=reg63*reg168; reg167=reg35*reg81; reg98=reg179+reg98; reg129=reg167-reg129;
    reg167=reg99*reg119; reg146=reg146/reg123; reg159=reg159/reg123; reg170=reg163*reg169; reg172=reg89*reg35;
    reg173=reg154*reg169; reg177=reg44*reg168; reg122=reg152+reg122; reg79=reg115-reg79; reg147=reg95-reg147;
    reg133=reg166+reg133; reg95=reg100*reg104; reg115=reg150*reg144; reg38=reg38/reg123; reg171=reg171/reg123;
    reg152=reg91*reg104; reg162=reg17+reg162; reg135=reg143-reg135; reg103=reg136+reg103; reg17=reg151*reg119;
    reg176=reg98+reg176; reg76=reg77*reg76; reg93=reg160+reg93; reg98=reg33*reg104; reg177=reg122+reg177;
    reg147=reg147/reg123; reg122=reg151*reg163; reg117=reg117/reg123; reg136=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg143=reg154*reg99;
    reg160=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg166=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; reg129=reg129/reg123; reg178=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; reg152=reg162+reg152;
    reg95=reg133+reg95; reg115=reg172-reg115; reg146=reg146*reg176; reg38=reg38*reg93; reg141=reg141/reg135;
    reg159=reg159*reg176; reg171=reg171*reg93; reg17=reg173-reg17; reg161=reg161/reg135; reg79=reg79/reg123;
    reg167=reg170-reg167; reg155=reg155/reg135; reg114=reg114/reg135; reg76=reg103+reg76; reg103=reg166*reg155;
    reg133=reg178*reg161; reg162=reg141*reg136; reg98=reg177+reg98; reg117=reg117*reg76; reg159=reg171-reg159;
    reg170=reg160*reg114; reg171=reg95*reg167; reg129=reg129*reg76; reg122=reg143-reg122; reg143=reg152*reg17;
    reg176=reg147*reg176; reg123=reg115/reg123; reg93=reg79*reg93; reg38=reg146-reg38; reg79=PNODE(1).dep[0]-PNODE(0).dep[0];
    reg129=reg38-reg129; reg176=reg93-reg176; reg178=reg141*reg178; reg38=PNODE(1).dep[1]-PNODE(0).dep[1]; reg166=reg166*reg114;
    reg93=PNODE(2).dep[0]-PNODE(0).dep[0]; reg160=reg160*reg155; reg115=PNODE(2).dep[1]-PNODE(0).dep[1]; reg136=reg136*reg161; reg159=reg117+reg159;
    reg76=reg123*reg76; reg117=reg163*reg98; reg123=reg169*reg95; reg133=reg103-reg133; reg103=reg119*reg95;
    reg146=reg151*reg98; reg147=reg152*reg119; reg172=reg99*reg98; reg173=reg152*reg169; reg177=reg122*reg98;
    reg179=reg154*reg98; reg180=1-(*f.m).resolution; reg143=reg171-reg143; reg170=reg162-reg170; reg162=(*f.m).deltaT*reg126;
    reg171=reg107*(*f.m).deltaT; reg159=reg180*reg159; reg172=reg173-reg172; reg170=reg133+reg170; reg177=reg143+reg177;
    reg107=(*f.m).resolution*reg107; reg176=reg76+reg176; reg76=reg114*reg115; reg129=reg180*reg129; reg126=(*f.m).resolution*reg126;
    reg117=reg147-reg117; reg133=reg151*reg152; reg143=reg141*reg38; reg166=reg178-reg166; elem.epsilon[0][0]=reg166;
    reg147=reg161*reg79; reg173=reg154*reg152; reg136=reg160-reg136; elem.epsilon[0][1]=reg136; reg146=reg123-reg146;
    reg123=reg155*reg93; reg160=reg163*reg95; reg178=reg99*reg95; reg179=reg103-reg179; reg103=reg106*(*f.m).deltaT;
    reg133=reg178-reg133; reg179=reg179/reg177; reg170=0.5*reg170; elem.epsilon[0][2]=reg170; reg178=reg136-reg171;
    T reg181=reg166-reg162; reg176=reg180*reg176; reg106=(*f.m).resolution*reg106; reg167=reg167/reg177; reg172=reg172/reg177;
    reg146=reg146/reg177; reg173=reg160-reg173; reg115=reg155*reg115; reg76=reg143-reg76; reg147=reg123-reg147;
    reg159=reg126+reg159; reg129=reg107+reg129; reg38=reg161*reg38; reg93=reg114*reg93; reg17=reg17/reg177;
    reg79=reg141*reg79; reg117=reg117/reg177; reg93=reg79-reg93; reg144=reg144*reg180; reg35=reg35*reg180;
    reg79=(*f.m).resolution*reg167; reg107=(*f.m).resolution*reg17; reg76=reg147+reg76; reg38=reg115-reg38; reg176=reg106+reg176;
    reg159=reg159*(*f.m).deltaT; reg129=reg129*(*f.m).deltaT; reg106=reg170-reg103; reg115=reg181*reg167; reg123=reg178*reg17;
    reg126=reg178*reg146; reg143=reg181*reg172; reg173=reg173/reg177; reg122=reg122/reg177; reg177=reg133/reg177;
    reg133=(*f.m).resolution*reg146; reg150=reg150*reg180; reg147=(*f.m).resolution*reg117; reg89=reg89*reg180; reg165=reg180*reg165;
    reg57=reg180*reg57; reg160=(*f.m).resolution*reg172; T reg182=(*f.m).resolution*reg179; T reg183=reg106*reg122; T reg184=(*f.m).resolution*reg177;
    reg123=reg115-reg123; reg115=reg38-reg129; T reg185=reg93-reg159; reg176=reg176*(*f.m).deltaT; reg153=reg180*reg153;
    reg182=reg57-reg182; reg147=reg165+reg147; reg133=reg89+reg133; reg160=reg150-reg160; reg143=reg126-reg143;
    reg107=reg144-reg107; reg35=reg79+reg35; reg76=0.5*reg76; reg81=reg180*reg81; reg57=(*f.m).resolution*reg122;
    reg180=reg120*reg180; reg79=reg106*reg177; reg89=(*f.m).resolution*reg173; reg79=reg143-reg79; reg123=reg183+reg123;
    reg120=reg182*reg115; reg126=reg147*reg185; reg143=reg133*reg115; reg180=reg57+reg180; reg184=reg81-reg184;
    reg57=reg160*reg185; reg89=reg153+reg89; reg81=reg76-reg176; reg144=reg35*reg185; reg150=reg107*reg115;
    reg153=reg79+reg123; reg150=reg144+reg150; reg144=reg180*reg81; reg165=reg89*reg81; reg120=reg126+reg120;
    reg126=reg184*reg81; reg143=reg57+reg143; reg165=reg120+reg165; reg144=reg150+reg144; reg126=reg143+reg126;
    reg178=reg178*reg179; reg153=reg153/3; reg181=reg181*reg117; reg165=2*reg165; reg129=reg136-reg129;
    reg123=reg123-reg153; reg159=reg166-reg159; reg185=reg185*reg144; reg178=reg181-reg178; reg106=reg106*reg173;
    reg79=reg79-reg153; reg115=reg115*reg126; reg35=reg159*reg35; reg178=reg106+reg178; reg107=reg129*reg107;
    reg115=reg185+reg115; reg81=reg81*reg165; reg79=pow(reg79,2); reg123=pow(reg123,2); reg133=reg129*reg133;
    reg160=reg159*reg160; reg176=reg170-reg176; reg184=reg176*reg184; reg147=reg159*reg147; reg182=reg129*reg182;
    reg133=reg160+reg133; reg180=reg176*reg180; reg107=reg35+reg107; reg81=reg115+reg81; reg153=pow(reg153,2);
    reg79=reg123+reg79; reg35=2*reg178; reg153=reg79+reg153; reg89=reg176*reg89; reg182=reg147+reg182;
    reg184=reg133+reg184; elem.sigma[0][1]=reg184; reg81=reg135*reg81; reg180=reg107+reg180; elem.sigma[0][0]=reg180;
    reg35=reg178*reg35; reg57=0.16666666666666665741*reg81; reg81=0.33333333333333331483*reg81; reg89=reg182+reg89; elem.sigma[0][2]=reg89;
    reg79=reg180*reg55; reg106=reg184*reg54; reg35=reg153+reg35; reg107=reg180*reg53; reg115=reg184*reg50;
    reg184=reg184*reg92; reg180=reg180*reg90; reg81=reg57+reg81; reg57=reg89*reg16; reg184=reg180+reg184;
    reg35=1.5*reg35; reg120=reg89*reg52; reg115=reg107+reg115; reg89=reg89*reg37; reg106=reg79+reg106;
    elem.sigma_von_mises=pow(reg35,0.5); elem.ener=reg81/2; elem.sigma_local[0][2]=reg184+reg57; elem.sigma_local[0][1]=reg115+reg120; elem.sigma_local[0][0]=reg106+reg89;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; T reg2=2*(*f.m).shear_modulus_12; reg1=1.0/reg1; reg0=1.0/reg0;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=1.0/(*f.m).elastic_modulus_3; T reg5=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg6=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg7=reg2*reg3; T reg8=pow((*f.m).v2[1],2); T reg9=pow((*f.m).v2[0],2); T reg10=pow((*f.m).v1[0],2); T reg11=pow((*f.m).v1[1],2);
    T reg12=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg13=pow((*f.m).v1[2],2); reg8=reg9+reg8; reg9=pow((*f.m).v2[2],2); T reg14=reg6*reg7;
    T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=reg5*reg7; reg11=reg10+reg11; reg10=reg4*reg7; T reg17=reg12*reg10;
    T reg18=reg5*reg16; T reg19=reg5*reg14; T reg20=reg15*reg10; reg9=reg8+reg9; reg13=reg11+reg13;
    reg8=1.0/(*f.m).elastic_modulus_1; reg13=pow(reg13,0.5); reg9=pow(reg9,0.5); reg11=reg15*reg14; T reg21=reg12*reg16;
    reg19=reg17+reg19; reg18=reg20-reg18; reg20=(*f.m).v2[1]/reg9; T reg22=reg12*reg19; T reg23=reg8*reg18;
    T reg24=reg21+reg11; T reg25=(*f.m).v2[2]/reg9; T reg26=(*f.m).v1[2]/reg13; T reg27=(*f.m).v1[1]/reg13; T reg28=reg26*reg20;
    T reg29=reg27*reg25; reg9=(*f.m).v2[0]/reg9; T reg30=reg6*reg3; T reg31=reg2*reg0; reg13=(*f.m).v1[0]/reg13;
    T reg32=reg4*reg3; T reg33=reg12*reg7; T reg34=reg6*reg14; reg10=reg8*reg10; reg3=reg5*reg3;
    reg22=reg23-reg22; reg23=reg6*reg24; T reg35=reg6*reg16; reg7=reg15*reg7; T reg36=reg15*reg32;
    reg23=reg22-reg23; reg22=reg29-reg28; T reg37=2*reg13; T reg38=2*reg9; T reg39=reg6*reg33;
    T reg40=reg31*reg4; T reg41=reg31*reg5; reg31=reg31*reg6; T reg42=reg2*reg1; T reg43=reg30*reg5;
    reg32=reg12*reg32; T reg44=reg5*reg3; reg16=reg8*reg16; T reg45=reg26*reg9; T reg46=reg13*reg25;
    reg14=reg12*reg14; reg34=reg10-reg34; reg35=reg17+reg35; reg10=reg6*reg7; reg17=reg42*reg6;
    T reg47=reg42*reg5; reg34=reg34/reg23; reg10=reg21+reg10; reg30=reg30*reg15; reg3=reg12*reg3;
    reg35=reg35/reg23; T reg48=2*reg22; T reg49=reg20*reg38; reg18=reg18/reg23; T reg50=pow(reg20,2);
    T reg51=pow(reg9,2); T reg52=reg27*reg37; T reg53=reg27*reg9; T reg54=pow(reg27,2); T reg55=pow(reg13,2);
    reg19=reg19/reg23; T reg56=reg13*reg20; T reg57=reg45-reg46; reg42=reg42*reg4; T reg58=reg40*reg12;
    T reg59=reg41*reg5; reg40=reg40*reg15; reg43=reg32+reg43; reg32=reg31*reg5; reg44=reg36-reg44;
    reg33=reg12*reg33; reg39=reg16+reg39; reg7=reg8*reg7; reg14=reg16+reg14; reg16=reg57*reg48;
    reg36=reg49*reg34; T reg60=reg52*reg35; T reg61=reg50*reg34; T reg62=reg54*reg35; T reg63=reg51*reg34;
    T reg64=reg55*reg35; reg32=reg58+reg32; reg58=reg49*reg19; T reg65=reg52*reg18; T reg66=reg42*reg15;
    reg42=reg42*reg12; T reg67=reg47*reg5; T reg68=reg56-reg53; T reg69=reg17*reg5; reg10=reg10/reg23;
    T reg70=reg54*reg18; T reg71=reg50*reg19; reg33=reg7-reg33; reg7=reg55*reg18; reg24=reg24/reg23;
    reg14=reg14/reg23; reg44=reg44*reg8; reg31=reg31*reg15; reg43=reg43*reg12; reg41=reg41*reg12;
    T reg72=reg3+reg30; reg39=reg39/reg23; T reg73=reg51*reg19; T reg74=pow(reg57,2); T reg75=pow(reg22,2);
    T reg76=pow(reg25,2); reg59=reg40-reg59; reg40=pow(reg26,2); T reg77=reg50*reg39; reg59=reg59*reg8;
    T reg78=reg16*reg24; reg58=reg65+reg58; reg65=reg76*reg19; T reg79=reg40*reg18; T reg80=pow(reg68,2);
    reg67=reg66-reg67; reg69=reg42+reg69; reg42=reg52*reg10; reg72=reg72*reg6; reg66=reg49*reg39;
    reg71=reg70+reg71; reg70=reg74*reg24; reg32=reg32*reg12; reg43=reg44-reg43; reg33=reg33/reg23;
    reg47=reg47*reg12; reg7=reg73+reg7; reg44=reg75*reg24; reg17=reg17*reg15; reg73=reg51*reg39;
    T reg81=reg55*reg10; T reg82=reg16*reg14; reg36=reg60+reg36; reg60=reg41+reg31; T reg83=reg76*reg34;
    T reg84=reg40*reg35; T reg85=reg74*reg14; reg61=reg62+reg61; reg63=reg64+reg63; reg62=reg75*reg14;
    reg64=reg54*reg10; T reg86=reg75*reg33; reg73=reg81+reg73; reg81=reg9*reg20; T reg87=reg13*reg27;
    reg82=reg36+reg82; reg60=reg60*reg6; reg78=reg58+reg78; reg62=reg63+reg62; reg36=reg26*reg37;
    reg70=reg71+reg70; reg58=2*reg27; reg63=reg80*reg14; reg83=reg84+reg83; reg71=reg25*reg38;
    reg44=reg7+reg44; reg7=reg47+reg17; reg84=2*reg20; reg69=reg69*reg12; T reg88=reg76*reg39;
    reg72=reg43-reg72; reg85=reg61+reg85; reg67=reg67*reg8; reg43=reg40*reg10; reg61=reg80*reg24;
    T reg89=reg74*reg33; reg77=reg64+reg77; reg65=reg79+reg65; reg32=reg59-reg32; reg66=reg42+reg66;
    reg16=reg16*reg33; reg42=reg26*reg58; reg59=reg36*reg18; reg64=reg71*reg19; reg79=reg57*reg22;
    T reg90=reg78*reg55; T reg91=reg85*reg50; T reg92=2*reg26; T reg93=reg54*reg70; T reg94=reg20*reg58;
    reg72=reg72/reg23; T reg95=reg9*reg37; T reg96=reg62*reg50; T reg97=reg54*reg44; reg60=reg32-reg60;
    reg32=reg25*reg84; reg69=reg67-reg69; reg67=reg13*reg57; reg7=reg7*reg6; T reg98=reg27*reg22;
    T reg99=reg82*reg51; reg86=reg73+reg86; reg56=reg53+reg56; reg53=reg27*reg20; reg73=reg13*reg9;
    T reg100=2*reg57; reg48=reg68*reg48; reg89=reg77+reg89; reg88=reg43+reg88; reg43=reg80*reg33;
    reg77=reg85*reg81; T reg101=reg82*reg50; T reg102=reg78*reg54; reg16=reg66+reg16; reg63=reg83+reg63;
    reg61=reg65+reg61; reg85=reg85*reg51; reg65=reg55*reg70; reg78=reg78*reg87; reg70=reg70*reg87;
    reg66=reg36*reg35; reg82=reg82*reg81; reg83=reg71*reg34; T reg103=reg55*reg44; T reg104=reg62*reg51;
    T reg105=reg86*reg75; T reg106=reg53*reg72; reg104=reg103+reg104; reg7=reg69-reg7; reg69=reg27*reg57;
    reg43=reg88+reg43; reg88=reg13*reg22; reg103=reg89*reg75; T reg107=reg16*reg75; reg85=reg65+reg85;
    reg65=reg73*reg72; T reg108=reg26*reg25; reg100=reg68*reg100; T reg109=reg51*reg8; T reg110=reg89*reg74;
    reg91=reg93+reg91; reg93=reg61*reg54; T reg111=reg63*reg50; T reg112=reg16*reg79; reg82=reg78+reg82;
    reg78=reg50*reg12; reg44=reg44*reg87; T reg113=reg86*reg74; reg62=reg62*reg81; reg96=reg97+reg96;
    reg60=reg60/reg23; reg97=reg95*reg8; T reg114=reg94*reg12; reg89=reg89*reg79; reg77=reg70+reg77;
    reg67=reg98+reg67; reg70=reg71*reg39; reg98=reg9*reg57; T reg115=reg20*reg22; reg83=reg66+reg83;
    reg66=reg94*reg15; T reg116=reg48*reg14; T reg117=reg95*reg12; reg18=reg42*reg18; reg99=reg90+reg99;
    reg90=reg48*reg24; reg19=reg32*reg19; T reg118=reg50*reg15; reg64=reg59+reg64; reg59=reg36*reg10;
    reg101=reg102+reg101; reg92=reg25*reg92; reg35=reg42*reg35; reg102=reg61*reg55; reg16=reg16*reg74;
    reg34=reg32*reg34; T reg119=reg51*reg12; T reg120=reg63*reg51; T reg121=reg56*reg72; T reg122=reg108*reg72;
    reg112=reg82+reg112; reg110=reg91+reg110; reg82=reg121*reg56; reg91=reg106*reg94; T reg123=reg43*reg74;
    reg34=reg35+reg34; reg70=reg59+reg70; reg111=reg93+reg111; reg116=reg83+reg116; reg105=reg104+reg105;
    reg35=reg65*reg95; reg120=reg102+reg120; reg59=reg43*reg75; reg83=reg121*reg94; reg16=reg101+reg16;
    reg103=reg85+reg103; reg85=reg106*reg95; reg93=reg67*reg60; reg69=reg69*reg60; reg89=reg77+reg89;
    reg106=reg106*reg56; reg61=reg61*reg87; reg63=reg63*reg81; reg86=reg86*reg79; reg62=reg44+reg62;
    reg88=reg88*reg60; reg44=reg92*reg5; reg39=reg32*reg39; reg23=reg7/reg23; reg7=reg20*reg57;
    reg77=reg9*reg22; reg119=reg118-reg119; reg101=reg26*reg68; reg121=reg121*reg95; reg107=reg99+reg107;
    reg19=reg18+reg19; reg24=reg100*reg24; reg18=reg76*reg5; reg99=reg94*reg5; reg102=reg95*reg6;
    reg114=reg97-reg114; reg117=reg66-reg117; reg66=reg50*reg5; reg97=reg51*reg6; reg104=reg92*reg6;
    reg14=reg100*reg14; reg48=reg48*reg33; reg98=reg115+reg98; reg115=reg65*reg94; reg10=reg42*reg10;
    reg58=reg57*reg58; reg90=reg64+reg90; reg113=reg96+reg113; reg64=reg76*reg6; reg78=reg109-reg78;
    reg37=reg22*reg37; reg96=reg122*reg94; reg109=reg54*reg12; reg118=reg93*reg58; reg64=reg78-reg64;
    reg15=reg54*reg15; reg8=reg55*reg8; reg104=reg114-reg104; reg12=reg55*reg12; reg83=reg16+reg83;
    reg16=reg90*reg54; reg14=reg34+reg14; reg45=reg46+reg45; reg34=reg23*reg98; reg7=reg23*reg7;
    reg77=reg23*reg77; reg65=reg65*reg56; reg86=reg62+reg86; reg46=reg50*(*f.m).alpha_2; reg62=reg26*reg22;
    reg18=reg119-reg18; reg78=reg13*reg68; reg123=reg111+reg123; reg111=(*f.m).alpha_1*reg54; reg114=reg93*reg67;
    reg82=reg112+reg82; reg112=reg122*reg95; reg119=reg69*reg58; reg91=reg110+reg91; reg84=reg57*reg84;
    reg110=reg88*reg58; reg115=reg113+reg115; reg48=reg70+reg48; reg38=reg22*reg38; reg39=reg10+reg39;
    reg10=reg25*reg68; reg70=reg116*reg51; reg113=reg90*reg55; reg33=reg100*reg33; reg93=reg93*reg37;
    reg121=reg107+reg121; reg100=(*f.m).alpha_1*reg55; reg99=reg102+reg99; reg92=reg92*reg4; reg66=reg97+reg66;
    reg97=reg76*reg4; reg44=reg117-reg44; reg24=reg19+reg24; reg101=reg101*reg60; reg19=reg69*reg67;
    reg106=reg89+reg106; reg69=reg69*reg37; reg89=reg116*reg50; reg63=reg61+reg63; reg43=reg43*reg79;
    reg61=reg51*(*f.m).alpha_2; reg85=reg103+reg85; reg59=reg120+reg59; reg102=reg88*reg37; reg35=reg105+reg35;
    reg99=reg92-reg99; reg12=reg15-reg12; reg15=reg104*reg55; reg92=reg40*reg5; reg103=reg54*reg44;
    reg105=reg64*reg55; reg66=reg97-reg66; reg97=reg18*reg54; reg107=reg27*reg68; reg117=reg34*reg98;
    reg114=reg82+reg114; reg82=reg45*reg72; reg90=reg90*reg87; reg120=(*f.m).alpha_3*reg75; reg116=reg116*reg81;
    reg59=reg112+reg59; reg61=reg100+reg61; reg100=reg101*reg37; reg112=reg76*(*f.m).alpha_2; T reg124=reg9*reg68;
    T reg125=(*f.m).alpha_1*reg40; reg46=reg111+reg46; reg111=(*f.m).alpha_3*reg74; T reg126=reg25*reg22; T reg127=reg18*reg50;
    T reg128=reg7*reg38; T reg129=reg64*reg51; reg118=reg83+reg118; reg109=reg8-reg109; reg8=reg40*reg6;
    reg119=reg91+reg119; reg83=reg34*reg84; reg91=reg101*reg58; reg96=reg123+reg96; reg123=reg7*reg84;
    T reg130=reg48*reg74; reg89=reg16+reg89; reg122=reg122*reg56; reg43=reg63+reg43; reg7=reg7*reg98;
    reg65=reg86+reg65; reg16=reg50*reg44; reg63=reg104*reg51; reg88=reg88*reg67; reg19=reg106+reg19;
    reg10=reg23*reg10; reg86=reg77*reg38; reg102=reg35+reg102; reg5=reg54*reg5; reg6=reg55*reg6;
    reg33=reg39+reg33; reg44=reg53*reg44; reg69=reg85+reg69; reg104=reg104*reg73; reg35=reg14*reg50;
    reg39=reg24*reg54; reg18=reg18*reg53; reg64=reg64*reg73; reg29=reg28+reg29; reg93=reg121+reg93;
    reg34=reg34*reg38; reg70=reg113+reg70; reg28=reg48*reg75; reg78=reg62+reg78; reg62=reg26*reg57;
    reg85=reg24*reg55; reg106=reg14*reg51; reg110=reg115+reg110; reg113=reg77*reg84; reg72=reg29*reg72;
    reg78=reg78*reg60; reg115=reg81*reg2; reg121=reg56*reg2; reg97=reg105+reg97; reg107=reg62+reg107;
    reg124=reg126+reg124; reg62=reg25*reg57; reg105=reg20*reg68; reg13=reg13*reg26; reg9=reg9*reg25;
    reg126=reg10*reg38; reg100=reg59+reg100; reg117=reg114+reg117; reg116=reg90+reg116; reg48=reg48*reg79;
    reg24=reg24*reg87; reg14=reg14*reg81; reg59=reg76*reg66; reg127=reg129+reg127; reg90=reg40*reg99;
    reg103=reg15+reg103; reg15=reg40*reg66; reg92=reg12-reg92; reg8=reg109-reg8; reg91=reg96+reg91;
    reg88=reg65+reg88; reg77=reg77*reg98; reg7=reg19+reg7; reg122=reg43+reg122; reg101=reg101*reg67;
    reg123=reg119+reg123; reg113=reg110+reg113; reg75=reg33*reg75; reg106=reg85+reg106; reg12=reg82*reg95;
    reg28=reg70+reg28; reg34=reg93+reg34; reg18=reg64+reg18; reg66=reg108*reg66; reg44=reg104+reg44;
    reg19=reg108*reg99; reg4=reg40*reg4; reg5=reg6+reg5; reg81=reg81*(*f.m).alpha_2; reg6=(*f.m).alpha_1*reg87;
    reg80=(*f.m).alpha_3*reg80; reg112=reg125+reg112; reg111=reg46+reg111; reg99=reg76*reg99; reg120=reg61+reg120;
    reg16=reg63+reg16; reg43=reg10*reg84; reg83=reg118+reg83; reg130=reg89+reg130; reg46=reg82*reg94;
    reg86=reg102+reg86; reg35=reg39+reg35; reg74=reg33*reg74; reg128=reg69+reg128; reg39=reg83*reg7;
    reg97=reg15+reg97; reg80=reg112+reg80; reg14=reg24+reg14; reg33=reg33*reg79; reg15=reg34*reg7;
    reg5=reg4-reg5; reg4=reg56*reg121; reg19=reg44+reg19; reg24=reg56*reg115; reg66=reg18+reg66;
    reg12=reg28+reg12; reg18=reg78*reg37; reg28=reg49*reg115; reg59=reg127+reg59; reg115=reg52*reg115;
    reg25=reg20*reg25; reg26=reg27*reg26; reg22=reg68*reg22; reg126=reg100+reg126; reg105=reg62+reg105;
    reg20=reg92*reg54; reg27=reg8*reg55; reg44=reg9*(*f.m).alpha_2; reg61=(*f.m).alpha_1*reg13; reg62=reg45*reg1;
    reg79=(*f.m).alpha_3*reg79; reg81=reg6+reg81; reg48=reg116+reg48; reg82=reg82*reg56; reg9=reg9*reg1;
    reg6=reg123*reg117; reg63=reg128*reg117; reg64=reg92*reg50; reg65=reg86*reg120; reg95=reg72*reg95;
    reg75=reg106+reg75; reg74=reg35+reg74; reg124=reg23*reg124; reg35=reg8*reg51; reg77=reg88+reg77;
    reg107=reg60*reg107; reg94=reg72*reg94; reg101=reg122+reg101; reg60=reg123*reg111; reg69=reg113*reg120;
    reg43=reg91+reg43; reg10=reg10*reg98; reg99=reg16+reg99; reg16=reg49*reg121; reg70=reg78*reg58;
    reg46=reg130+reg46; reg90=reg103+reg90; reg85=reg128*reg111; reg121=reg52*reg121; reg95=reg75+reg95;
    reg75=reg124*reg38; reg88=reg25*reg0; reg89=reg29*reg0; reg20=reg27+reg20; reg40=reg40*reg5;
    reg10=reg101+reg10; reg27=reg111*reg7; reg91=reg77*reg120; reg57=reg68*reg57; reg126=reg126*reg80;
    reg85=reg65+reg85; reg18=reg12+reg18; reg12=reg71*reg62; reg16=reg99+reg16; reg8=reg8*reg73;
    reg92=reg92*reg53; reg24=reg66+reg24; reg65=reg45*reg9; reg60=reg69+reg60; reg4=reg19+reg4;
    reg19=reg45*reg62; reg66=reg34*reg123; reg68=reg128*reg83; reg15=reg63-reg15; reg105=reg23*reg105;
    reg39=reg6-reg39; reg43=reg80*reg43; reg37=reg107*reg37; reg2=reg87*reg2; reg6=reg71*reg9;
    reg78=reg78*reg67; reg28=reg59+reg28; reg121=reg90+reg121; reg62=reg36*reg62; reg9=reg36*reg9;
    reg115=reg97+reg115; reg33=reg14+reg33; reg82=reg48+reg82; reg79=reg81+reg79; reg44=reg61+reg44;
    reg22=(*f.m).alpha_3*reg22; reg76=reg76*reg5; reg64=reg35+reg64; reg70=reg46+reg70; reg14=reg124*reg84;
    reg25=reg25*(*f.m).alpha_2; reg23=(*f.m).alpha_1*reg26; reg72=reg72*reg56; reg94=reg74+reg94; reg58=reg107*reg58;
    reg35=reg113*reg15; reg72=reg33+reg72; reg67=reg107*reg67; reg33=reg86*reg39; reg124=reg124*reg98;
    reg78=reg82+reg78; reg1=reg13*reg1; reg22=reg44+reg22; reg13=reg79*reg83; reg40=reg20+reg40;
    reg20=reg52*reg2; reg43=reg60+reg43; reg10=reg80*reg10; reg58=reg94+reg58; reg14=reg70+reg14;
    reg75=reg18+reg75; reg27=reg91+reg27; reg84=reg105*reg84; reg76=reg64+reg76; reg18=reg49*reg2;
    reg12=reg16+reg12; reg16=reg42*reg89; reg62=reg121+reg62; reg92=reg8+reg92; reg5=reg108*reg5;
    reg8=reg32*reg89; reg6=reg28+reg6; reg28=reg32*reg88; reg65=reg24+reg65; reg24=reg29*reg88;
    reg19=reg4+reg19; reg89=reg29*reg89; reg66=reg68-reg66; reg4=reg34*reg79; reg126=reg85+reg126;
    reg25=reg23+reg25; reg57=(*f.m).alpha_3*reg57; reg9=reg115+reg9; reg88=reg42*reg88; reg37=reg95+reg37;
    reg38=reg105*reg38; reg35=reg33-reg35; reg67=reg72+reg67; reg98=reg105*reg98; reg23=reg66*reg77;
    reg89=reg19+reg89; reg28=reg6+reg28; reg6=reg113*reg117; reg19=reg83*reg77; reg24=reg65+reg24;
    reg33=reg86*reg117; reg2=reg56*reg2; reg5=reg92+reg5; reg10=reg27+reg10; reg27=reg79*reg117;
    reg84=reg58+reg84; reg57=reg25+reg57; reg16=reg62+reg16; reg38=reg37+reg38; reg25=reg34*reg77;
    reg71=reg71*reg1; reg18=reg76+reg18; reg88=reg9+reg88; reg20=reg40+reg20; reg36=reg36*reg1;
    reg13=reg43+reg13; reg0=reg26*reg0; reg8=reg12+reg8; reg14=reg22*reg14; reg75=reg75*reg22;
    reg4=reg126+reg4; reg124=reg78+reg124; reg19=reg6-reg19; reg6=reg123*reg77; reg1=reg45*reg1;
    reg2=reg5+reg2; reg38=reg38*reg57; reg36=reg20+reg36; reg5=reg88*reg89; reg23=reg35+reg23;
    reg9=reg86*reg7; reg25=reg33-reg25; reg12=reg128*reg77; reg20=reg86*reg83; reg26=reg34*reg113;
    reg33=reg16*reg24; reg32=reg32*reg0; reg71=reg18+reg71; reg14=reg13+reg14; reg27=reg10+reg27;
    reg42=reg42*reg0; reg98=reg67+reg98; reg10=reg8*reg24; reg84=reg57*reg84; reg124=reg22*reg124;
    reg75=reg4+reg75; reg4=reg113*reg7; reg13=reg28*reg89; reg18=reg88*reg8; reg32=reg71+reg32;
    reg39=reg39/reg23; reg26=reg20-reg26; reg20=reg123*reg86; reg22=reg113*reg128; reg12=reg9-reg12;
    reg84=reg14+reg84; reg25=reg25/reg23; reg9=reg16*reg28; reg15=reg15/reg23; reg6=reg4-reg6;
    reg124=reg27+reg124; reg19=reg19/reg23; reg0=reg29*reg0; reg1=reg2+reg1; reg98=reg57*reg98;
    reg33=reg5-reg33; reg42=reg36+reg42; reg10=reg13-reg10; reg38=reg75+reg38; reg98=reg124+reg98;
    reg15=reg15*reg84; reg2=reg32*reg33; reg0=reg1+reg0; reg9=reg18-reg9; reg26=reg26/reg23;
    reg66=reg66/reg23; reg39=reg39*reg38; reg12=reg12/reg23; reg19=reg19*reg38; reg1=reg42*reg10;
    reg6=reg6/reg23; reg25=reg25*reg84; reg22=reg20-reg22; reg15=reg39-reg15; reg66=reg66*reg98;
    reg38=reg6*reg38; reg84=reg12*reg84; reg19=reg25-reg19; reg23=reg22/reg23; reg26=reg26*reg98;
    reg4=reg32*reg24; reg5=reg28*reg0; reg6=reg88*reg0; reg24=reg24*reg42; reg12=reg9*reg0;
    reg2=reg1-reg2; reg28=reg28*reg42; reg1=reg8*reg0; reg13=reg89*reg42; reg89=reg32*reg89;
    reg88=reg88*reg32; reg5=reg4-reg5; reg0=reg16*reg0; reg12=reg2+reg12; reg15=reg66+reg15;
    reg2=1-(*f.m).resolution; reg6=reg24-reg6; reg26=reg19-reg26; reg84=reg38-reg84; reg98=reg23*reg98;
    reg4=elem.pos(2)[0]-elem.pos(0)[0]; reg26=reg2*reg26; reg84=reg98+reg84; reg1=reg89-reg1; reg5=reg5/reg12;
    reg0=reg13-reg0; reg6=reg6/reg12; reg42=reg8*reg42; reg32=reg16*reg32; reg88=reg28-reg88;
    reg8=elem.pos(1)[1]-elem.pos(0)[1]; reg13=elem.pos(1)[0]-elem.pos(0)[0]; reg14=(*f.m).resolution*reg111; reg16=elem.pos(2)[1]-elem.pos(0)[1]; reg18=(*f.m).resolution*reg120;
    reg15=reg2*reg15; reg84=reg2*reg84; reg19=(*f.m).resolution*reg79; reg15=reg18+reg15; reg26=reg14+reg26;
    reg10=reg10/reg12; reg1=reg1/reg12; reg14=(*f.m).resolution*reg5; reg18=(*f.m).resolution*reg6; reg33=reg33/reg12;
    reg0=reg0/reg12; reg7=reg2*reg7; reg77=reg2*reg77; reg32=reg42-reg32; reg88=reg88/reg12;
    reg20=reg8*reg4; reg22=reg16*reg13; reg9=reg9/reg12; reg12=reg32/reg12; reg20=reg22-reg20;
    reg22=(*f.m).resolution*reg10; reg23=(*f.m).resolution*reg33; reg15=reg15*(*f.m).deltaT; reg84=reg19+reg84; reg26=reg26*(*f.m).deltaT;
    reg18=reg7-reg18; reg14=reg77+reg14; reg7=(*f.m).resolution*reg88; reg19=(*f.m).resolution*reg0; reg86=reg86*reg2;
    reg128=reg128*reg2; reg113=reg113*reg2; reg123=reg123*reg2; reg24=(*f.m).resolution*reg1; reg117=reg2*reg117;
    reg34=reg34*reg2; reg84=reg84*(*f.m).deltaT; reg7=reg117+reg7; reg4=reg4/reg20; reg25=(*f.m).resolution*reg12;
    reg16=reg16/reg20; reg83=reg2*reg83; reg19=reg123+reg19; reg24=reg113-reg24; reg2=reg14*reg15;
    reg27=reg18*reg26; reg23=reg128-reg23; reg28=(*f.m).resolution*reg9; reg13=reg13/reg20; reg8=reg8/reg20;
    reg86=reg22+reg86; reg22=reg4-reg13; reg29=reg19*reg26; reg32=reg24*reg15; reg35=reg8-reg16;
    reg36=reg86*reg15; reg37=reg23*reg26; reg25=reg83-reg25; reg34=reg28+reg34; reg28=reg7*reg84;
    reg38=reg2+reg27; reg39=reg25*reg84; reg40=reg32+reg29; reg42=0.5*reg4; reg43=0.5*reg22;
    reg44=0.5*reg16; reg45=0.5*reg35; reg46=reg36+reg37; reg48=0.5*reg13; reg57=reg34*reg84;
    reg58=reg38+reg28; reg59=0.5*reg8; reg60=reg4*reg18; reg61=reg45*reg7; reg62=reg7*reg44;
    reg63=reg22*reg18; reg64=reg42*reg7; reg65=reg16*reg14; reg66=2*reg58; reg67=reg8*reg14;
    reg68=reg48*reg7; reg69=reg46+reg57; reg70=reg59*reg7; reg71=reg13*reg18; reg72=reg40+reg39;
    reg74=reg7*reg43; reg75=reg35*reg14; reg61=reg63+reg61; reg63=reg42*reg66; reg76=reg16*reg69;
    reg77=reg22*reg19; reg78=reg13*reg19; reg80=reg44*reg66; reg81=reg4*reg72; reg82=reg16*reg24;
    reg83=reg59*reg66; reg85=reg13*reg72; reg87=reg42*reg25; reg89=reg8*reg69; reg90=reg4*reg19;
    reg91=reg48*reg25; reg92=reg25*reg44; reg93=reg48*reg66; reg94=reg8*reg24; reg95=1-var_inter[0];
    reg96=reg34*reg43; reg97=reg35*reg86; reg74=reg75+reg74; reg75=reg45*reg25; reg98=reg59*reg34;
    reg99=reg13*reg23; reg100=reg59*reg25; reg71=reg71-reg70; reg68=reg68-reg67; reg101=reg48*reg34;
    reg102=reg8*reg86; reg103=reg22*reg23; reg104=reg35*reg24; reg62=reg62-reg60; reg105=reg45*reg34;
    reg106=reg34*reg44; reg107=reg4*reg23; reg108=reg25*reg43; reg65=reg65-reg64; reg109=reg42*reg34;
    reg110=reg16*reg86; reg112=reg89-reg93; reg113=reg83-reg85; reg74=2*reg74; reg78=reg78-reg100;
    reg105=reg103+reg105; reg96=reg97+reg96; reg97=var_inter[1]*elem.f_vol_e[1]; reg99=reg99-reg98; reg68=2*reg68;
    reg101=reg101-reg102; reg62=2*reg62; reg103=reg81-reg80; reg114=reg63-reg76; reg115=reg22*reg72;
    reg116=reg45*reg66; reg91=reg91-reg94; reg71=2*reg71; reg95=reg95-var_inter[1]; reg61=2*reg61;
    reg117=reg66*reg43; reg77=reg75+reg77; reg75=reg35*reg69; reg108=reg104+reg108; reg92=reg92-reg90;
    reg110=reg110-reg109; reg106=reg106-reg107; reg104=var_inter[0]*elem.f_vol_e[0]; reg82=reg82-reg87; reg118=var_inter[0]*elem.f_vol_e[1];
    reg119=var_inter[1]*elem.f_vol_e[0]; reg65=2*reg65; reg121=reg35*reg106; reg122=reg48*reg71; reg123=reg62*reg43;
    reg124=reg99*reg8; reg125=reg71*reg43; reg126=reg22*reg108; reg127=reg101*reg8; reg128=reg65*reg43;
    reg129=reg35*reg110; reg130=reg8*reg106; T reg131=reg48*reg62; T reg132=reg45*reg74; T reg133=reg8*reg110;
    T reg134=reg48*reg65; T reg135=reg115+reg116; T reg136=reg105*reg8; T reg137=reg35*reg105; T reg138=reg95*elem.f_vol_e[0];
    T reg139=reg75+reg117; T reg140=reg74*reg43; T reg141=reg59*reg71; T reg142=reg13*reg78; T reg143=reg59*reg68;
    T reg144=reg13*reg91; T reg145=reg35*reg96; T reg146=reg59*reg62; T reg147=reg13*reg92; T reg148=reg95*elem.f_vol_e[1];
    T reg149=reg35*reg99; T reg150=reg42*reg61; T reg151=reg59*reg65; T reg152=reg13*reg82; T reg153=reg22*reg77;
    T reg154=reg68*reg43; T reg155=reg13*reg77; T reg156=reg59*reg74; T reg157=reg13*reg108; T reg158=reg35*reg101;
    T reg159=reg62*reg42; reg106=reg16*reg106; T reg160=reg22*reg92; T reg161=reg42*reg65; reg110=reg16*reg110;
    T reg162=reg45*reg62; reg105=reg16*reg105; reg92=reg4*reg92; reg62=reg62*reg44; T reg163=reg22*reg91;
    reg112=reg112-reg119; T reg164=reg4*reg82; T reg165=reg65*reg44; T reg166=reg42*reg74; T reg167=reg16*reg96;
    reg77=reg4*reg77; T reg168=reg61*reg44; T reg169=reg59*reg61; T reg170=reg48*reg68; reg113=reg113-reg97;
    reg108=reg4*reg108; T reg171=reg45*reg68; T reg172=reg22*reg78; T reg173=reg45*reg71; T reg174=reg48*reg74;
    reg114=reg114-reg104; reg96=reg96*reg8; T reg175=reg45*reg61; reg78=reg4*reg78; T reg176=reg44*reg71;
    reg74=reg74*reg44; reg71=reg42*reg71; T reg177=reg61*reg43; reg103=reg103-reg118; reg91=reg4*reg91;
    T reg178=reg68*reg44; reg99=reg16*reg99; reg82=reg22*reg82; reg68=reg68*reg42; reg61=reg48*reg61;
    reg101=reg16*reg101; reg65=reg45*reg65; T reg179=reg148+reg135; reg113=reg20*reg113; reg171=reg163+reg171;
    reg132=reg126+reg132; reg114=reg20*reg114; reg175=reg153+reg175; reg112=reg20*reg112; reg77=reg168-reg77;
    reg162=reg160+reg162; reg103=reg20*reg103; reg125=reg149+reg125; reg65=reg82+reg65; reg154=reg158+reg154;
    reg156=reg157-reg156; reg96=reg174-reg96; reg140=reg145+reg140; reg124=reg122-reg124; reg173=reg172+reg173;
    reg164=reg165-reg164; reg128=reg129+reg128; reg166=reg167-reg166; reg92=reg62-reg92; reg177=reg137+reg177;
    reg150=reg105-reg150; reg78=reg176-reg78; reg71=reg99-reg71; reg161=reg110-reg161; reg68=reg101-reg68;
    reg159=reg106-reg159; reg91=reg178-reg91; reg136=reg61-reg136; reg61=reg139+reg138; reg133=reg134-reg133;
    reg141=reg142-reg141; reg108=reg74-reg108; reg130=reg131-reg130; reg143=reg144-reg143; reg127=reg170-reg127;
    reg169=reg155-reg169; reg151=reg152-reg151; reg121=reg123+reg121; reg146=reg147-reg146; reg136=reg20*reg136;
    reg91=reg20*reg91; reg161=reg20*reg161; reg62=reg20*reg179; reg96=reg20*reg96; reg150=reg20*reg150;
    reg154=reg20*reg154; reg162=reg20*reg162; reg133=reg20*reg133; reg121=reg20*reg121; reg166=reg20*reg166;
    reg140=reg20*reg140; reg78=reg20*reg78; reg127=reg20*reg127; reg114=ponderation*reg114; reg173=reg20*reg173;
    reg130=reg20*reg130; reg171=reg20*reg171; reg74=reg20*reg61; reg113=ponderation*reg113; reg141=reg20*reg141;
    reg143=reg20*reg143; reg108=reg20*reg108; reg146=reg20*reg146; reg151=reg20*reg151; reg77=reg20*reg77;
    reg169=reg20*reg169; reg156=reg20*reg156; reg124=reg20*reg124; reg125=reg20*reg125; reg112=ponderation*reg112;
    reg164=reg20*reg164; reg128=reg20*reg128; reg177=reg20*reg177; reg65=reg20*reg65; reg132=reg20*reg132;
    reg92=reg20*reg92; reg159=reg20*reg159; reg71=reg20*reg71; reg175=reg20*reg175; reg68=reg20*reg68;
    reg103=ponderation*reg103; sollicitation[indices[2]+1]+=-reg113; T tmp_3_3=ponderation*reg92; reg82=ponderation*reg62; sollicitation[indices[0]+1]+=reg82;
    T tmp_3_4=ponderation*reg91; T tmp_4_0=ponderation*reg96; T tmp_3_0=ponderation*reg108; sollicitation[indices[1]+1]+=-reg103; T tmp_3_1=ponderation*reg77;
    sollicitation[indices[1]+0]+=-reg114; sollicitation[indices[2]+0]+=-reg112; T tmp_3_5=ponderation*reg78; T tmp_3_2=ponderation*reg164; reg77=ponderation*reg74;
    sollicitation[indices[0]+0]+=reg77; T tmp_5_5=ponderation*reg141; T tmp_5_4=ponderation*reg143; T tmp_5_3=ponderation*reg146; T tmp_5_2=ponderation*reg151;
    T tmp_5_1=ponderation*reg169; T tmp_5_0=ponderation*reg156; T tmp_4_5=ponderation*reg124; T tmp_0_5=ponderation*reg125; T tmp_4_4=ponderation*reg127;
    T tmp_0_2=ponderation*reg128; T tmp_0_1=ponderation*reg177; T tmp_1_0=ponderation*reg132; T tmp_2_5=ponderation*reg71; T tmp_1_1=ponderation*reg175;
    T tmp_2_4=ponderation*reg68; T tmp_2_3=ponderation*reg159; T tmp_1_2=ponderation*reg65; T tmp_2_2=ponderation*reg161; T tmp_2_1=ponderation*reg150;
    T tmp_1_3=ponderation*reg162; T tmp_2_0=ponderation*reg166; T tmp_1_5=ponderation*reg173; T tmp_1_4=ponderation*reg171; T tmp_0_0=ponderation*reg140;
    T tmp_0_4=ponderation*reg154; T tmp_0_3=ponderation*reg121; T tmp_4_3=ponderation*reg130; T tmp_4_2=ponderation*reg133; T tmp_4_1=ponderation*reg136;
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
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; T reg2=2*(*f.m).shear_modulus_12; reg1=1.0/reg1; reg0=1.0/reg0;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=1.0/(*f.m).elastic_modulus_3; T reg5=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg6=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg7=reg2*reg3; T reg8=pow((*f.m).v2[1],2); T reg9=pow((*f.m).v2[0],2); T reg10=pow((*f.m).v1[0],2); T reg11=pow((*f.m).v1[1],2);
    T reg12=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg13=pow((*f.m).v1[2],2); reg8=reg9+reg8; reg9=pow((*f.m).v2[2],2); T reg14=reg6*reg7;
    T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=reg5*reg7; reg11=reg10+reg11; reg10=reg4*reg7; T reg17=reg12*reg10;
    T reg18=reg5*reg16; T reg19=reg5*reg14; T reg20=reg15*reg10; reg9=reg8+reg9; reg13=reg11+reg13;
    reg8=1.0/(*f.m).elastic_modulus_1; reg13=pow(reg13,0.5); reg9=pow(reg9,0.5); reg11=reg15*reg14; T reg21=reg12*reg16;
    reg19=reg17+reg19; reg18=reg20-reg18; reg20=(*f.m).v2[1]/reg9; T reg22=reg12*reg19; T reg23=reg8*reg18;
    T reg24=reg21+reg11; T reg25=(*f.m).v2[2]/reg9; T reg26=(*f.m).v1[2]/reg13; T reg27=(*f.m).v1[1]/reg13; T reg28=reg26*reg20;
    T reg29=reg27*reg25; reg9=(*f.m).v2[0]/reg9; T reg30=reg6*reg3; T reg31=reg2*reg0; reg13=(*f.m).v1[0]/reg13;
    T reg32=reg4*reg3; T reg33=reg12*reg7; T reg34=reg6*reg14; reg10=reg8*reg10; reg3=reg5*reg3;
    reg22=reg23-reg22; reg23=reg6*reg24; T reg35=reg6*reg16; reg7=reg15*reg7; T reg36=reg15*reg32;
    reg23=reg22-reg23; reg22=reg29-reg28; T reg37=2*reg13; T reg38=2*reg9; T reg39=reg6*reg33;
    T reg40=reg31*reg4; T reg41=reg31*reg5; reg31=reg31*reg6; T reg42=reg2*reg1; T reg43=reg30*reg5;
    reg32=reg12*reg32; T reg44=reg5*reg3; reg16=reg8*reg16; T reg45=reg26*reg9; T reg46=reg13*reg25;
    reg14=reg12*reg14; reg34=reg10-reg34; reg35=reg17+reg35; reg10=reg6*reg7; reg17=reg42*reg6;
    T reg47=reg42*reg5; reg34=reg34/reg23; reg10=reg21+reg10; reg30=reg30*reg15; reg3=reg12*reg3;
    reg35=reg35/reg23; T reg48=2*reg22; T reg49=reg20*reg38; reg18=reg18/reg23; T reg50=pow(reg20,2);
    T reg51=pow(reg9,2); T reg52=reg27*reg37; T reg53=reg27*reg9; T reg54=pow(reg27,2); T reg55=pow(reg13,2);
    reg19=reg19/reg23; T reg56=reg13*reg20; T reg57=reg45-reg46; reg42=reg42*reg4; T reg58=reg40*reg12;
    T reg59=reg41*reg5; reg40=reg40*reg15; reg43=reg32+reg43; reg32=reg31*reg5; reg44=reg36-reg44;
    reg33=reg12*reg33; reg39=reg16+reg39; reg7=reg8*reg7; reg14=reg16+reg14; reg16=reg57*reg48;
    reg36=reg49*reg34; T reg60=reg52*reg35; T reg61=reg50*reg34; T reg62=reg54*reg35; T reg63=reg51*reg34;
    T reg64=reg55*reg35; reg32=reg58+reg32; reg58=reg49*reg19; T reg65=reg52*reg18; T reg66=reg42*reg15;
    reg42=reg42*reg12; T reg67=reg47*reg5; T reg68=reg56-reg53; T reg69=reg17*reg5; reg10=reg10/reg23;
    T reg70=reg54*reg18; T reg71=reg50*reg19; reg33=reg7-reg33; reg7=reg55*reg18; reg24=reg24/reg23;
    reg14=reg14/reg23; reg44=reg44*reg8; reg31=reg31*reg15; reg43=reg43*reg12; reg41=reg41*reg12;
    T reg72=reg3+reg30; reg39=reg39/reg23; T reg73=reg51*reg19; T reg74=pow(reg57,2); T reg75=pow(reg22,2);
    T reg76=pow(reg25,2); reg59=reg40-reg59; reg40=pow(reg26,2); T reg77=reg50*reg39; reg59=reg59*reg8;
    T reg78=reg16*reg24; reg58=reg65+reg58; reg65=reg76*reg19; T reg79=reg40*reg18; T reg80=pow(reg68,2);
    reg67=reg66-reg67; reg69=reg42+reg69; reg42=reg52*reg10; reg72=reg72*reg6; reg66=reg49*reg39;
    reg71=reg70+reg71; reg70=reg74*reg24; reg32=reg32*reg12; reg43=reg44-reg43; reg33=reg33/reg23;
    reg47=reg47*reg12; reg7=reg73+reg7; reg44=reg75*reg24; reg17=reg17*reg15; reg73=reg51*reg39;
    T reg81=reg55*reg10; T reg82=reg16*reg14; T reg83=reg41+reg31; reg36=reg60+reg36; reg60=reg76*reg34;
    T reg84=reg40*reg35; T reg85=reg74*reg14; reg61=reg62+reg61; reg63=reg64+reg63; reg62=reg75*reg14;
    reg64=reg54*reg10; T reg86=reg75*reg33; reg73=reg81+reg73; reg81=reg9*reg20; T reg87=reg13*reg27;
    reg82=reg36+reg82; reg83=reg83*reg6; reg78=reg58+reg78; reg62=reg63+reg62; reg70=reg71+reg70;
    reg36=reg26*reg37; reg58=2*reg27; reg63=reg80*reg14; reg60=reg84+reg60; reg71=reg25*reg38;
    reg44=reg7+reg44; reg7=reg47+reg17; reg84=2*reg20; reg69=reg69*reg12; T reg88=reg76*reg39;
    reg72=reg43-reg72; reg85=reg61+reg85; reg67=reg67*reg8; reg43=reg40*reg10; reg61=reg80*reg24;
    T reg89=reg74*reg33; reg77=reg64+reg77; reg65=reg79+reg65; reg32=reg59-reg32; reg66=reg42+reg66;
    reg16=reg16*reg33; reg16=reg66+reg16; reg63=reg60+reg63; reg42=reg85*reg51; reg59=reg55*reg70;
    reg60=reg82*reg51; reg64=reg71*reg19; reg66=reg36*reg18; reg61=reg65+reg61; reg65=reg62*reg51;
    reg79=reg25*reg84; T reg90=reg55*reg44; T reg91=reg26*reg58; T reg92=reg57*reg22; T reg93=reg78*reg55;
    T reg94=2*reg26; T reg95=reg82*reg81; reg82=reg82*reg50; T reg96=reg20*reg58; T reg97=reg78*reg54;
    reg78=reg78*reg87; T reg98=reg9*reg37; T reg99=reg70*reg87; T reg100=reg85*reg81; T reg101=reg27*reg22;
    reg56=reg53+reg56; reg7=reg7*reg6; reg53=reg27*reg20; T reg102=reg13*reg9; reg69=reg67-reg69;
    reg67=2*reg57; reg48=reg68*reg48; reg83=reg32-reg83; reg32=reg80*reg33; reg88=reg43+reg88;
    reg89=reg77+reg89; reg43=reg54*reg44; reg86=reg73+reg86; reg73=reg13*reg57; reg85=reg85*reg50;
    reg77=reg36*reg35; T reg103=reg71*reg34; reg70=reg54*reg70; reg72=reg72/reg23; T reg104=reg62*reg50;
    T reg105=reg16*reg75; T reg106=reg89*reg92; reg7=reg69-reg7; reg65=reg90+reg65; reg104=reg43+reg104;
    reg83=reg83/reg23; reg43=reg86*reg75; reg85=reg70+reg85; reg69=reg89*reg74; reg95=reg78+reg95;
    reg89=reg89*reg75; reg70=reg16*reg92; reg42=reg59+reg42; reg59=reg63*reg50; reg78=reg86*reg74;
    reg90=reg61*reg54; T reg107=reg27*reg57; T reg108=reg13*reg22; T reg109=reg26*reg25; reg67=reg68*reg67;
    reg32=reg88+reg32; reg88=reg51*reg12; T reg110=reg50*reg15; T reg111=reg98*reg12; reg34=reg79*reg34;
    reg35=reg91*reg35; T reg112=reg48*reg14; reg103=reg77+reg103; reg77=reg96*reg15; T reg113=reg36*reg10;
    T reg114=reg71*reg39; T reg115=reg102*reg72; T reg116=reg53*reg72; T reg117=reg56*reg72; reg19=reg79*reg19;
    reg18=reg91*reg18; T reg118=reg48*reg24; reg64=reg66+reg64; reg73=reg101+reg73; reg66=reg96*reg12;
    reg101=reg98*reg8; T reg119=reg50*reg12; T reg120=reg51*reg8; reg16=reg16*reg74; reg82=reg97+reg82;
    reg44=reg44*reg87; reg62=reg62*reg81; reg100=reg99+reg100; reg97=reg63*reg51; reg94=reg25*reg94;
    reg99=reg61*reg55; reg60=reg93+reg60; reg93=reg9*reg57; T reg121=reg20*reg22; reg16=reg82+reg16;
    reg58=reg57*reg58; reg37=reg22*reg37; reg62=reg44+reg62; reg86=reg86*reg92; reg106=reg100+reg106;
    reg44=reg116*reg56; reg61=reg61*reg87; reg63=reg63*reg81; reg93=reg121+reg93; reg70=reg95+reg70;
    reg82=reg117*reg56; reg95=reg32*reg74; reg59=reg90+reg59; reg90=reg20*reg57; reg100=reg9*reg22;
    reg43=reg65+reg43; reg65=reg115*reg98; reg89=reg42+reg89; reg42=reg116*reg98; reg116=reg116*reg96;
    reg69=reg85+reg69; reg85=reg26*reg68; reg97=reg99+reg97; reg99=reg32*reg75; reg121=reg76*reg5;
    reg88=reg110-reg88; reg110=reg94*reg5; reg14=reg67*reg14; reg34=reg35+reg34; reg112=reg103+reg112;
    reg114=reg113+reg114; reg48=reg48*reg33; reg10=reg91*reg10; reg39=reg79*reg39; reg35=reg109*reg72;
    reg24=reg67*reg24; reg19=reg18+reg19; reg118=reg64+reg118; reg108=reg108*reg83; reg107=reg107*reg83;
    reg18=reg73*reg83; reg66=reg101-reg66; reg64=reg94*reg6; reg119=reg120-reg119; reg101=reg76*reg6;
    reg103=reg117*reg96; reg23=reg7/reg23; reg7=reg50*reg5; reg113=reg96*reg5; reg120=reg51*reg6;
    reg105=reg60+reg105; reg111=reg77-reg111; reg60=reg115*reg96; reg78=reg104+reg78; reg117=reg117*reg98;
    reg77=reg98*reg6; reg104=reg107*reg73; reg14=reg34+reg14; reg95=reg59+reg95; reg63=reg61+reg63;
    reg32=reg32*reg92; reg34=reg35*reg96; reg121=reg88-reg121; reg110=reg111-reg110; reg24=reg19+reg24;
    reg19=reg54*reg12; reg59=reg18*reg73; reg61=reg25*reg68; reg15=reg54*reg15; reg88=reg18*reg37;
    reg101=reg119-reg101; reg111=reg50*(*f.m).alpha_2; reg48=reg114+reg48; reg38=reg22*reg38; reg114=reg118*reg55;
    reg84=reg57*reg84; reg7=reg120+reg7; reg119=reg112*reg51; reg39=reg10+reg39; reg10=reg76*reg4;
    reg86=reg62+reg86; reg115=reg115*reg56; reg103=reg16+reg103; reg18=reg18*reg58; reg33=reg67*reg33;
    reg16=reg118*reg54; reg62=reg112*reg50; reg44=reg106+reg44; reg65=reg43+reg65; reg43=reg108*reg37;
    reg67=reg35*reg98; reg64=reg66-reg64; reg66=reg51*(*f.m).alpha_2; reg99=reg97+reg99; reg45=reg46+reg45;
    reg42=reg89+reg42; reg46=reg107*reg37; reg89=reg13*reg68; reg97=reg26*reg22; reg100=reg23*reg100;
    reg107=reg107*reg58; reg116=reg69+reg116; reg69=(*f.m).alpha_1*reg55; reg85=reg85*reg83; reg106=reg108*reg58;
    reg60=reg78+reg60; reg82=reg70+reg82; reg117=reg105+reg117; reg8=reg55*reg8; reg70=reg23*reg93;
    reg78=(*f.m).alpha_1*reg54; reg90=reg23*reg90; reg12=reg55*reg12; reg94=reg94*reg4; reg113=reg77+reg113;
    reg77=reg64*reg102; reg105=reg53*reg110; reg120=reg45*reg72; reg88=reg117+reg88; reg117=reg85*reg37;
    T reg122=reg121*reg53; reg59=reg82+reg59; reg19=reg8-reg19; reg8=reg40*reg6; reg82=reg101*reg102;
    reg89=reg97+reg89; reg18=reg103+reg18; reg97=reg70*reg84; reg103=reg70*reg93; T reg123=reg14*reg50;
    reg62=reg16+reg62; reg16=reg48*reg74; T reg124=reg26*reg57; reg99=reg67+reg99; reg33=reg39+reg33;
    reg39=reg24*reg54; reg66=reg69+reg66; reg118=reg118*reg87; reg70=reg70*reg38; reg67=reg85*reg58;
    reg69=(*f.m).alpha_3*reg75; reg34=reg95+reg34; reg6=reg55*reg6; reg113=reg94-reg113; reg94=reg24*reg55;
    reg95=reg14*reg51; T reg125=reg76*(*f.m).alpha_2; reg12=reg15-reg12; reg15=reg50*reg110; T reg126=reg64*reg51;
    reg61=reg23*reg61; T reg127=(*f.m).alpha_1*reg40; T reg128=reg40*reg5; T reg129=reg27*reg68; reg43=reg65+reg43;
    reg65=reg100*reg38; T reg130=reg101*reg51; T reg131=reg121*reg54; reg101=reg101*reg55; reg5=reg54*reg5;
    reg46=reg42+reg46; reg42=reg90*reg38; T reg132=reg90*reg84; reg107=reg116+reg107; reg29=reg28+reg29;
    reg28=reg100*reg84; reg106=reg60+reg106; reg7=reg10-reg7; reg108=reg108*reg73; reg115=reg86+reg115;
    reg119=reg114+reg119; reg112=reg112*reg81; reg10=(*f.m).alpha_3*reg74; reg60=reg9*reg68; reg86=reg48*reg75;
    reg104=reg44+reg104; reg90=reg90*reg93; reg111=reg78+reg111; reg32=reg63+reg32; reg35=reg35*reg56;
    reg121=reg121*reg50; reg110=reg54*reg110; reg64=reg64*reg55; reg44=reg25*reg22; reg4=reg40*reg4;
    reg5=reg6+reg5; reg72=reg29*reg72; reg112=reg118+reg112; reg48=reg48*reg92; reg105=reg77+reg105;
    reg24=reg24*reg87; reg14=reg14*reg81; reg6=reg76*reg7; reg63=reg109*reg113; reg121=reg130+reg121;
    reg10=reg111+reg10; reg77=reg81*(*f.m).alpha_2; reg78=(*f.m).alpha_1*reg87; reg125=reg127+reg125; reg122=reg82+reg122;
    reg82=reg109*reg7; reg128=reg12-reg128; reg80=(*f.m).alpha_3*reg80; reg103=reg59+reg103; reg110=reg64+reg110;
    reg12=reg40*reg113; reg108=reg115+reg108; reg100=reg100*reg93; reg59=reg20*reg68; reg64=reg25*reg57;
    reg60=reg44+reg60; reg90=reg104+reg90; reg86=reg119+reg86; reg44=reg120*reg98; reg35=reg32+reg35;
    reg85=reg85*reg73; reg95=reg94+reg95; reg75=reg33*reg75; reg113=reg76*reg113; reg15=reg126+reg15;
    reg132=reg107+reg132; reg129=reg124+reg129; reg65=reg43+reg65; reg131=reg101+reg131; reg81=reg81*reg2;
    reg42=reg46+reg42; reg28=reg106+reg28; reg32=reg56*reg2; reg13=reg13*reg26; reg8=reg19-reg8;
    reg9=reg9*reg25; reg67=reg34+reg67; reg74=reg33*reg74; reg123=reg39+reg123; reg19=reg120*reg96;
    reg16=reg62+reg16; reg69=reg66+reg69; reg97=reg18+reg97; reg70=reg88+reg70; reg117=reg99+reg117;
    reg89=reg89*reg83; reg18=reg61*reg38; reg7=reg40*reg7; reg34=reg61*reg84; reg39=reg9*reg1;
    reg43=reg45*reg1; reg75=reg95+reg75; reg5=reg4-reg5; reg22=reg68*reg22; reg98=reg72*reg98;
    reg4=reg8*reg55; reg26=reg27*reg26; reg25=reg20*reg25; reg20=reg128*reg50; reg27=reg52*reg32;
    reg46=reg128*reg54; reg62=reg8*reg51; reg12=reg110+reg12; reg66=reg52*reg81; reg131=reg7+reg131;
    reg59=reg64+reg59; reg48=reg112+reg48; reg120=reg120*reg56; reg14=reg24+reg14; reg33=reg33*reg92;
    reg7=reg49*reg81; reg6=reg121+reg6; reg100=reg108+reg100; reg81=reg56*reg81; reg63=reg105+reg63;
    reg24=reg56*reg32; reg82=reg122+reg82; reg85=reg35+reg85; reg61=reg61*reg93; reg44=reg86+reg44;
    reg35=reg89*reg37; reg18=reg117+reg18; reg113=reg15+reg113; reg60=reg23*reg60; reg32=reg49*reg32;
    reg15=reg132*reg103; reg34=reg67+reg34; reg80=reg125+reg80; reg64=reg42*reg103; reg19=reg16+reg19;
    reg16=reg97*reg90; reg67=reg70*reg90; reg96=reg72*reg96; reg74=reg123+reg74; reg86=(*f.m).alpha_1*reg13;
    reg88=reg132*reg10; reg94=reg28*reg69; reg9=reg9*(*f.m).alpha_2; reg95=reg65*reg69; reg99=reg42*reg10;
    reg101=reg89*reg58; reg92=(*f.m).alpha_3*reg92; reg129=reg83*reg129; reg77=reg78+reg77; reg67=reg64-reg67;
    reg16=reg15-reg16; reg15=reg42*reg97; reg64=reg36*reg39; reg66=reg131+reg66; reg57=reg68*reg57;
    reg68=reg25*(*f.m).alpha_2; reg78=(*f.m).alpha_1*reg26; reg120=reg48+reg120; reg89=reg89*reg73; reg22=(*f.m).alpha_3*reg22;
    reg9=reg86+reg9; reg33=reg14+reg33; reg72=reg72*reg56; reg99=reg95+reg99; reg128=reg128*reg53;
    reg8=reg8*reg102; reg18=reg18*reg80; reg81=reg82+reg81; reg88=reg94+reg88; reg34=reg80*reg34;
    reg35=reg44+reg35; reg14=reg45*reg39; reg44=reg60*reg38; reg98=reg75+reg98; reg37=reg129*reg37;
    reg48=reg100*reg69; reg75=reg10*reg90; reg24=reg63+reg24; reg63=reg45*reg43; reg59=reg23*reg59;
    reg32=reg113+reg32; reg23=reg71*reg43; reg61=reg85+reg61; reg2=reg87*reg2; reg25=reg25*reg0;
    reg82=reg29*reg0; reg46=reg4+reg46; reg40=reg40*reg5; reg4=reg70*reg132; reg76=reg76*reg5;
    reg43=reg36*reg43; reg27=reg12+reg27; reg20=reg62+reg20; reg101=reg19+reg101; reg12=reg60*reg84;
    reg58=reg129*reg58; reg39=reg71*reg39; reg92=reg77+reg92; reg7=reg6+reg7; reg96=reg74+reg96;
    reg76=reg20+reg76; reg6=reg79*reg82; reg19=reg79*reg25; reg64=reg66+reg64; reg20=reg70*reg92;
    reg18=reg99+reg18; reg62=reg49*reg2; reg34=reg88+reg34; reg66=reg92*reg97; reg23=reg32+reg23;
    reg37=reg98+reg37; reg32=reg29*reg82; reg63=reg24+reg63; reg24=reg65*reg16; reg44=reg35+reg44;
    reg14=reg81+reg14; reg35=reg29*reg25; reg74=reg28*reg67; reg39=reg7+reg39; reg73=reg129*reg73;
    reg72=reg33+reg72; reg1=reg13*reg1; reg61=reg80*reg61; reg12=reg101+reg12; reg22=reg9+reg22;
    reg75=reg48+reg75; reg60=reg60*reg93; reg89=reg120+reg89; reg58=reg96+reg58; reg84=reg59*reg84;
    reg68=reg78+reg68; reg57=(*f.m).alpha_3*reg57; reg4=reg15-reg4; reg38=reg59*reg38; reg25=reg91*reg25;
    reg5=reg109*reg5; reg40=reg46+reg40; reg128=reg8+reg128; reg43=reg27+reg43; reg82=reg91*reg82;
    reg7=reg52*reg2; reg8=reg92*reg103; reg25=reg64+reg25; reg74=reg24-reg74; reg32=reg63+reg32;
    reg6=reg23+reg6; reg36=reg36*reg1; reg7=reg40+reg7; reg9=reg70*reg100; reg13=reg4*reg100;
    reg15=reg65*reg103; reg0=reg26*reg0; reg23=reg28*reg103; reg24=reg97*reg100; reg73=reg72+reg73;
    reg93=reg59*reg93; reg60=reg89+reg60; reg57=reg68+reg57; reg2=reg56*reg2; reg5=reg128+reg5;
    reg84=reg58+reg84; reg82=reg43+reg82; reg20=reg18+reg20; reg44=reg44*reg22; reg66=reg34+reg66;
    reg12=reg22*reg12; reg61=reg75+reg61; reg38=reg37+reg38; reg62=reg76+reg62; reg71=reg71*reg1;
    reg19=reg39+reg19; reg35=reg14+reg35; reg24=reg23-reg24; reg1=reg45*reg1; reg8=reg61+reg8;
    reg14=reg132*reg100; reg18=reg65*reg90; reg2=reg5+reg2; reg79=reg79*reg0; reg5=reg70*reg28;
    reg9=reg15-reg9; reg60=reg22*reg60; reg71=reg62+reg71; reg15=reg65*reg97; reg22=reg42*reg100;
    reg44=reg20+reg44; reg38=reg38*reg57; reg12=reg66+reg12; reg84=reg57*reg84; reg20=reg28*reg90;
    reg23=reg82*reg35; reg13=reg74+reg13; reg91=reg91*reg0; reg26=reg6*reg35; reg27=reg25*reg32;
    reg33=reg19*reg32; reg36=reg7+reg36; reg93=reg73+reg93; reg7=reg82*reg19; reg22=reg18-reg22;
    reg93=reg57*reg93; reg60=reg8+reg60; reg26=reg33-reg26; reg8=reg28*reg42; reg38=reg44+reg38;
    reg84=reg12+reg84; reg5=reg15-reg5; reg79=reg71+reg79; reg12=reg132*reg65; reg24=reg24/reg13;
    reg0=reg29*reg0; reg14=reg20-reg14; reg67=reg67/reg13; reg16=reg16/reg13; reg15=reg25*reg6;
    reg91=reg36+reg91; reg1=reg2+reg1; reg9=reg9/reg13; reg23=reg27-reg23; reg0=reg1+reg0;
    reg5=reg5/reg13; reg14=reg14/reg13; reg22=reg22/reg13; reg93=reg60+reg93; reg1=reg79*reg23;
    reg16=reg16*reg38; reg67=reg67*reg84; reg8=reg12-reg8; reg2=reg91*reg26; reg9=reg9*reg84;
    reg24=reg24*reg38; reg4=reg4/reg13; reg7=reg15-reg7; reg1=reg2-reg1; reg2=reg7*reg0;
    reg12=reg35*reg91; reg15=reg19*reg0; reg35=reg79*reg35; reg13=reg8/reg13; reg38=reg14*reg38;
    reg8=reg25*reg0; reg84=reg22*reg84; reg4=reg4*reg93; reg67=reg16-reg67; reg5=reg5*reg93;
    reg24=reg9-reg24; reg9=reg79*reg32; reg14=reg6*reg0; reg25=reg25*reg79; reg19=reg19*reg91;
    reg32=reg32*reg91; reg67=reg4+reg67; reg15=reg35-reg15; reg2=reg1+reg2; reg1=1-(*f.m).resolution;
    reg0=reg82*reg0; reg84=reg38-reg84; reg93=reg13*reg93; reg5=reg24-reg5; reg8=reg12-reg8;
    reg84=reg93+reg84; reg67=reg1*reg67; reg4=(*f.m).resolution*reg10; reg25=reg19-reg25; reg12=elem.pos(1)[1]-elem.pos(0)[1];
    reg13=(*f.m).resolution*reg69; reg8=reg8/reg2; reg0=reg32-reg0; reg91=reg6*reg91; reg79=reg82*reg79;
    reg15=reg15/reg2; reg6=elem.pos(2)[1]-elem.pos(0)[1]; reg14=reg9-reg14; reg5=reg1*reg5; reg9=elem.pos(2)[0]-elem.pos(0)[0];
    reg16=elem.pos(1)[0]-elem.pos(0)[0]; reg23=reg23/reg2; reg79=reg91-reg79; reg14=reg14/reg2; reg0=reg0/reg2;
    reg26=reg26/reg2; reg25=reg25/reg2; reg18=(*f.m).resolution*reg92; reg19=(*f.m).resolution*reg15; reg20=(*f.m).resolution*reg8;
    reg22=reg12*reg9; reg24=reg6*reg16; reg84=reg1*reg84; reg67=reg13+reg67; reg5=reg4+reg5;
    reg90=reg1*reg90; reg100=reg1*reg100; reg5=reg5*(*f.m).deltaT; reg67=reg67*(*f.m).deltaT; reg4=(*f.m).resolution*reg25;
    reg84=reg18+reg84; reg13=(*f.m).resolution*reg23; reg19=reg100+reg19; reg18=(*f.m).resolution*reg26; reg20=reg90-reg20;
    reg27=(*f.m).resolution*reg0; reg29=(*f.m).resolution*reg14; reg103=reg1*reg103; reg22=reg24-reg22; reg132=reg132*reg1;
    reg28=reg28*reg1; reg7=reg7/reg2; reg42=reg42*reg1; reg65=reg65*reg1; reg2=reg79/reg2;
    reg84=reg84*(*f.m).deltaT; reg4=reg103+reg4; reg24=reg19*reg67; reg32=reg20*reg5; reg70=reg70*reg1;
    reg9=reg9/reg22; reg16=reg16/reg22; reg6=reg6/reg22; reg97=reg1*reg97; reg1=(*f.m).resolution*reg2;
    reg12=reg12/reg22; reg33=(*f.m).resolution*reg7; reg65=reg18+reg65; reg13=reg42-reg13; reg29=reg28-reg29;
    reg27=reg132+reg27; reg18=reg4*reg84; reg28=reg24+reg32; reg34=reg29*reg67; reg35=reg13*reg5;
    reg36=reg27*reg5; reg37=reg65*reg67; reg38=reg9-reg16; reg39=reg12-reg6; reg70=reg33+reg70;
    reg1=reg97-reg1; reg33=0.5*reg9; reg40=0.5*reg12; reg42=reg37+reg35; reg43=reg70*reg84;
    reg44=0.5*reg6; reg45=0.5*reg39; reg46=reg28+reg18; reg48=0.5*reg16; reg57=reg1*reg84;
    reg58=reg34+reg36; reg59=0.5*reg38; reg60=reg33*reg4; reg61=reg6*reg19; reg62=reg9*reg20;
    reg63=reg58+reg57; reg64=reg16*reg20; reg66=reg4*reg59; reg68=reg42+reg43; reg71=reg40*reg4;
    reg72=reg4*reg44; reg73=reg12*reg19; reg74=reg48*reg4; reg75=reg39*reg19; reg76=2*reg46;
    reg77=reg38*reg20; reg78=reg45*reg4; reg64=reg64-reg71; reg79=reg48*reg76; reg80=reg12*reg68;
    reg81=reg6*reg68; reg78=reg77+reg78; reg77=reg38*reg27; reg82=reg9*reg63; reg83=reg6*reg29;
    reg85=reg33*reg1; reg86=reg44*reg76; reg87=reg33*reg76; reg88=reg9*reg27; reg89=reg1*reg44;
    reg66=reg75+reg66; reg75=reg12*reg29; reg90=reg48*reg1; reg91=reg70*reg59; reg93=reg16*reg27;
    reg94=reg39*reg65; reg95=reg48*reg70; reg96=reg12*reg65; reg97=reg45*reg70; reg98=reg40*reg1;
    reg74=reg74-reg73; reg72=reg72-reg62; reg99=reg70*reg44; reg100=reg9*reg13; reg101=reg40*reg70;
    reg103=1-var_inter[0]; reg104=reg38*reg13; reg61=reg61-reg60; reg105=reg33*reg70; reg106=reg6*reg65;
    reg107=reg16*reg13; reg108=reg40*reg76; reg109=reg45*reg1; reg110=reg16*reg63; reg91=reg94+reg91;
    reg97=reg104+reg97; reg94=var_inter[1]*elem.f_vol_e[0]; reg66=2*reg66; reg104=var_inter[0]*elem.f_vol_e[1]; reg111=reg87-reg81;
    reg112=var_inter[0]*elem.f_vol_e[0]; reg113=reg82-reg86; reg114=reg108-reg110; reg115=reg80-reg79; reg74=2*reg74;
    reg116=reg76*reg59; reg95=reg95-reg96; reg117=reg39*reg68; reg107=reg107-reg101; reg118=reg38*reg63;
    reg72=2*reg72; reg103=reg103-var_inter[1]; reg99=reg99-reg100; reg61=2*reg61; reg106=reg106-reg105;
    reg119=reg45*reg76; reg77=reg109+reg77; reg78=2*reg78; reg64=2*reg64; reg109=var_inter[1]*elem.f_vol_e[1];
    reg93=reg93-reg98; reg90=reg90-reg75; reg83=reg83-reg85; reg89=reg89-reg88; reg120=reg45*reg74;
    reg121=reg38*reg90; reg122=reg45*reg72; reg123=reg38*reg89; reg124=reg45*reg61; reg125=reg38*reg83;
    reg111=reg111-reg112; reg114=reg114-reg109; reg126=reg118+reg119; reg127=reg45*reg78; reg115=reg115-reg94;
    reg128=reg64*reg59; reg129=reg38*reg77; reg113=reg113-reg104; reg130=reg78*reg59; reg131=reg33*reg64;
    reg132=reg6*reg107; T reg133=reg74*reg33; T reg134=reg6*reg95; T reg135=reg72*reg33; T reg136=reg6*reg99;
    T reg137=reg74*reg44; T reg138=reg9*reg90; T reg139=reg33*reg61; T reg140=reg6*reg106; T reg141=reg48*reg74;
    T reg142=reg45*reg64; T reg143=reg38*reg93; T reg144=reg39*reg97; T reg145=reg66*reg59; T reg146=reg39*reg91;
    T reg147=reg39*reg107; T reg148=reg74*reg59; T reg149=reg39*reg95; T reg150=reg72*reg59; T reg151=reg44*reg64;
    reg95=reg95*reg12; T reg152=reg9*reg93; T reg153=reg9*reg89; T reg154=reg103*elem.f_vol_e[0]; T reg155=reg117+reg116;
    T reg156=reg103*elem.f_vol_e[1]; reg107=reg107*reg12; T reg157=reg61*reg59; T reg158=reg48*reg64; reg93=reg16*reg93;
    reg64=reg40*reg64; T reg159=reg39*reg106; T reg160=reg39*reg99; T reg161=reg72*reg44; reg139=reg140-reg139;
    reg107=reg158-reg107; reg122=reg123+reg122; reg160=reg150+reg160; reg123=reg155+reg154; reg138=reg137-reg138;
    reg148=reg149+reg148; reg142=reg143+reg142; reg114=reg22*reg114; reg120=reg121+reg120; reg113=reg22*reg113;
    reg145=reg146+reg145; reg64=reg93-reg64; reg153=reg161-reg153; reg130=reg144+reg130; reg111=reg22*reg111;
    reg152=reg151-reg152; reg131=reg132-reg131; reg157=reg159+reg157; reg133=reg134-reg133; reg128=reg147+reg128;
    reg93=reg156+reg126; reg127=reg129+reg127; reg95=reg141-reg95; reg115=reg22*reg115; reg135=reg136-reg135;
    reg124=reg125+reg124; reg121=reg22*reg93; reg148=reg22*reg148; reg152=reg22*reg152; reg95=reg22*reg95;
    reg138=reg22*reg138; reg114=ponderation*reg114; reg160=reg22*reg160; reg125=reg22*reg123; reg64=reg22*reg64;
    reg111=ponderation*reg111; reg130=reg22*reg130; reg153=reg22*reg153; reg131=reg22*reg131; reg133=reg22*reg133;
    reg128=reg22*reg128; reg115=ponderation*reg115; reg157=reg22*reg157; reg127=reg22*reg127; reg135=reg22*reg135;
    reg124=reg22*reg124; reg139=reg22*reg139; reg122=reg22*reg122; reg107=reg22*reg107; reg142=reg22*reg142;
    reg120=reg22*reg120; reg113=ponderation*reg113; reg145=reg22*reg145; T tmp_3_3=ponderation*reg153; T tmp_3_4=ponderation*reg138;
    sollicitation[indices[2]+1]+=-reg114; sollicitation[indices[1]+1]+=-reg113; sollicitation[indices[2]+0]+=-reg115; reg113=ponderation*reg125; sollicitation[indices[0]+0]+=reg113;
    T tmp_5_5=ponderation*reg64; T tmp_4_5=ponderation*reg107; T tmp_4_4=ponderation*reg95; T tmp_0_2=ponderation*reg157; T tmp_0_1=ponderation*reg130;
    T tmp_2_5=ponderation*reg131; T tmp_2_4=ponderation*reg133; T tmp_0_5=ponderation*reg128; T tmp_1_1=ponderation*reg127; T tmp_2_3=ponderation*reg135;
    T tmp_1_2=ponderation*reg124; T tmp_2_2=ponderation*reg139; T tmp_1_3=ponderation*reg122; T tmp_1_5=ponderation*reg142; T tmp_1_4=ponderation*reg120;
    T tmp_0_0=ponderation*reg145; T tmp_0_4=ponderation*reg148; T tmp_0_3=ponderation*reg160; reg64=ponderation*reg121; sollicitation[indices[0]+1]+=reg64;
    T tmp_3_5=ponderation*reg152; sollicitation[indices[1]+0]+=-reg111;
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
    T reg0=pow((*f.m).v1[1],2); T reg1=pow((*f.m).v1[0],2); T reg2=pow((*f.m).v2[1],2); reg0=reg1+reg0; reg1=pow((*f.m).v1[2],2);
    T reg3=pow((*f.m).v2[0],2); reg2=reg3+reg2; reg3=pow((*f.m).v2[2],2); T reg4=2*(*f.m).shear_modulus_23; T reg5=2*(*f.m).shear_modulus_13;
    reg1=reg0+reg1; reg0=2*(*f.m).shear_modulus_12; reg1=pow(reg1,0.5); reg5=1.0/reg5; reg4=1.0/reg4;
    reg3=reg2+reg3; reg0=1.0/reg0; reg2=(*f.m).v1[0]/reg1; T reg6=(*f.m).v1[1]/reg1; T reg7=reg5*reg4;
    reg3=pow(reg3,0.5); T reg8=2*reg2; T reg9=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg10=1.0/(*f.m).elastic_modulus_3; T reg11=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg12=2*reg6; T reg13=(*f.m).v2[0]/reg3; T reg14=(*f.m).v2[1]/reg3; reg1=(*f.m).v1[2]/reg1; T reg15=reg0*reg7;
    T reg16=reg9*reg15; T reg17=pow(reg13,2); T reg18=pow(reg14,2); T reg19=reg10*reg15; T reg20=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1;
    T reg21=1.0/(*f.m).elastic_modulus_2; T reg22=reg11*reg15; T reg23=1.0/(*f.m).elastic_modulus_1; T reg24=2*reg1; T reg25=reg14*reg12;
    reg3=(*f.m).v2[2]/reg3; T reg26=reg13*reg8; reg24=reg3*reg24; T reg27=reg17*reg23; T reg28=reg18*reg20;
    T reg29=reg26*reg23; T reg30=reg25*reg20; T reg31=reg17*reg20; T reg32=pow(reg3,2); T reg33=reg18*reg21;
    T reg34=reg26*reg20; T reg35=reg25*reg21; T reg36=reg21*reg19; T reg37=reg11*reg22; T reg38=reg20*reg19;
    T reg39=reg11*reg16; T reg40=reg24*reg9; reg28=reg27-reg28; reg27=reg32*reg9; reg30=reg29-reg30;
    reg29=pow(reg2,2); T reg41=reg32*reg11; reg31=reg33-reg31; reg33=pow(reg6,2); T reg42=reg24*reg11;
    T reg43=reg21*reg16; T reg44=reg20*reg22; reg34=reg35-reg34; reg35=reg17*reg9; T reg45=reg18*reg11;
    reg37=reg36-reg37; reg36=reg26*reg9; T reg46=reg25*reg11; reg39=reg38+reg39; T reg47=reg6*reg14;
    T reg48=reg2*reg14; T reg49=reg2*reg13; reg41=reg31-reg41; reg45=reg35+reg45; reg31=reg6*reg13;
    reg35=reg32*reg10; T reg50=reg33*reg21; T reg51=reg29*reg20; reg40=reg30-reg40; reg24=reg24*reg10;
    reg30=reg29*reg23; reg46=reg36+reg46; reg42=reg34-reg42; reg27=reg28-reg27; reg28=reg44+reg43;
    reg34=reg23*reg37; reg36=reg20*reg39; T reg52=reg33*reg20; T reg53=pow(reg1,2); T reg54=reg18*reg42;
    T reg55=reg9*reg7; T reg56=reg0*reg4; T reg57=reg9*reg28; T reg58=reg11*reg7; reg36=reg34-reg36;
    reg46=reg24-reg46; reg24=reg33*reg42; reg7=reg10*reg7; reg45=reg35-reg45; reg34=reg20*reg15;
    reg35=reg1*reg14; T reg59=reg33*reg11; T reg60=reg27*reg49; T reg61=reg29*reg9; T reg62=reg40*reg29;
    T reg63=reg41*reg47; T reg64=reg40*reg49; reg42=reg47*reg42; reg19=reg23*reg19; T reg65=reg6*reg3;
    reg15=reg21*reg15; T reg66=reg27*reg17; reg52=reg30-reg52; reg30=reg53*reg11; T reg67=reg53*reg9;
    reg51=reg50-reg51; reg50=reg9*reg22; T reg68=reg41*reg18; T reg69=reg9*reg16; T reg70=reg13*reg14;
    T reg71=reg2*reg3; T reg72=reg1*reg13; T reg73=reg31+reg48; reg27=reg27*reg29; T reg74=reg1*reg3;
    reg41=reg41*reg33; reg40=reg40*reg17; T reg75=2*reg13; T reg76=reg14*reg75; T reg77=reg74*reg45;
    reg63=reg60+reg63; reg60=reg56*reg11; T reg78=reg56*reg9; T reg79=reg0*reg5; T reg80=reg20*reg7;
    T reg81=reg71+reg72; reg56=reg56*reg10; T reg82=reg11*reg58; T reg83=reg65-reg35; T reg84=reg32*reg46;
    T reg85=reg55*reg11; reg54=reg40+reg54; reg40=reg6*reg8; reg59=reg61+reg59; reg61=reg70*reg0;
    T reg86=reg73*reg0; reg41=reg27+reg41; reg30=reg51-reg30; reg67=reg52-reg67; reg27=reg13*reg3;
    reg57=reg36-reg57; reg36=reg32*reg45; reg68=reg66+reg68; reg51=reg9*reg15; reg50=reg38+reg50;
    reg38=reg53*reg46; reg24=reg62+reg24; reg69=reg19-reg69; reg22=reg23*reg22; reg7=reg21*reg7;
    reg45=reg53*reg45; reg19=reg9*reg34; reg16=reg20*reg16; reg52=reg53*reg10; reg42=reg64+reg42;
    reg46=reg74*reg46; reg69=reg69/reg57; reg71=reg72-reg71; reg62=reg67*reg17; reg64=reg30*reg18;
    reg51=reg44+reg51; reg50=reg50/reg57; reg59=reg52-reg59; reg36=reg68+reg36; reg52=reg76*reg61;
    reg37=reg37/reg57; reg39=reg39/reg57; reg66=reg14*reg3; reg68=reg2*reg6; reg72=reg3*reg75;
    T reg87=reg40*reg86; reg38=reg24+reg38; reg24=reg73*reg86; reg77=reg63+reg77; reg63=reg73*reg61;
    reg84=reg54+reg84; reg86=reg76*reg86; reg61=reg40*reg61; reg41=reg45+reg41; reg45=reg1*reg8;
    reg46=reg42+reg46; reg19=reg22+reg19; reg65=reg35+reg65; reg15=reg23*reg15; reg16=reg22+reg16;
    reg34=reg20*reg34; reg22=reg56*reg20; reg35=reg60*reg11; reg55=reg55*reg21; reg42=reg81*reg5;
    reg54=reg78*reg11; T reg88=reg79*reg11; reg10=reg79*reg10; T reg89=reg67*reg29; reg58=reg20*reg58;
    T reg90=reg30*reg33; reg56=reg56*reg21; reg79=reg79*reg9; T reg91=2*reg83; T reg92=reg27*reg5;
    reg85=reg80+reg85; reg82=reg7-reg82; reg7=2*reg14; reg19=reg19/reg57; reg80=reg10*reg21;
    T reg93=reg58+reg55; T reg94=reg1*reg12; reg85=reg85*reg20; reg16=reg16/reg57; T reg95=reg3*reg7;
    reg35=reg56-reg35; reg54=reg22+reg54; reg34=reg15-reg34; reg15=reg71*reg91; reg22=reg76*reg69;
    reg56=reg40*reg50; T reg96=reg18*reg69; T reg97=reg33*reg50; T reg98=reg17*reg69; T reg99=reg29*reg50;
    T reg100=reg76*reg39; T reg101=reg40*reg37; T reg102=reg53*reg59; T reg103=reg79*reg11; T reg104=reg33*reg37;
    T reg105=reg18*reg39; T reg106=pow(reg83,2); T reg107=reg2*reg1; T reg108=pow(reg71,2); T reg109=reg17*reg39;
    reg90=reg89+reg90; reg67=reg67*reg49; reg30=reg30*reg47; reg89=reg65*reg4; reg51=reg51/reg57;
    T reg110=reg81*reg42; reg24=reg46+reg24; reg46=reg66*reg4; T reg111=reg29*reg37; reg60=reg60*reg20;
    reg28=reg28/reg57; reg78=reg78*reg21; reg11=reg88*reg11; reg10=reg10*reg20; reg82=reg82*reg23;
    reg0=reg68*reg0; reg63=reg77+reg63; reg77=reg81*reg92; reg86=reg84+reg86; reg84=reg72*reg42;
    reg64=reg62+reg64; reg42=reg45*reg42; reg87=reg38+reg87; reg38=reg32*reg59; reg62=reg72*reg92;
    reg52=reg36+reg52; reg92=reg45*reg92; reg61=reg41+reg61; reg36=reg40*reg0; reg102=reg90+reg102;
    reg5=reg107*reg5; reg41=reg65*reg89; reg110=reg24+reg110; reg24=reg65*reg46; reg77=reg63+reg77;
    reg59=reg74*reg59; reg30=reg67+reg30; reg63=reg108*reg28; reg105=reg104+reg105; reg103=reg10+reg103;
    reg11=reg80-reg11; reg10=reg18*reg19; reg67=reg33*reg51; reg80=reg17*reg19; reg90=reg29*reg51;
    reg104=reg15*reg16; reg22=reg56+reg22; reg56=reg108*reg16; reg96=reg97+reg96; reg97=reg40*reg51;
    T reg112=reg76*reg19; T reg113=reg106*reg16; reg98=reg99+reg98; reg99=reg6*reg1; T reg114=reg15*reg28;
    reg100=reg101+reg100; reg101=reg95*reg89; reg84=reg86+reg84; reg93=reg93*reg9; reg85=reg82-reg85;
    reg21=reg79*reg21; reg79=reg106*reg28; reg34=reg34/reg57; reg111=reg109+reg111; reg88=reg88*reg20;
    reg35=reg35*reg23; reg92=reg61+reg92; reg61=reg94*reg46; reg42=reg87+reg42; reg89=reg94*reg89;
    reg38=reg64+reg38; reg64=reg76*reg0; reg82=reg60+reg78; reg54=reg54*reg20; reg46=reg95*reg46;
    reg62=reg52+reg62; reg54=reg35-reg54; reg61=reg92+reg61; reg101=reg84+reg101; reg104=reg22+reg104;
    reg15=reg15*reg34; reg112=reg97+reg112; reg82=reg82*reg9; reg79=reg111+reg79; reg22=reg88+reg21;
    reg114=reg100+reg114; reg35=reg108*reg34; reg56=reg96+reg56; reg20=reg103*reg20; reg113=reg98+reg113;
    reg59=reg30+reg59; reg93=reg85-reg93; reg10=reg67+reg10; reg30=reg106*reg34; reg24=reg77+reg24;
    reg80=reg90+reg80; reg52=reg72*reg5; reg41=reg110+reg41; reg64=reg38+reg64; reg4=reg99*reg4;
    reg0=reg73*reg0; reg63=reg105+reg63; reg36=reg102+reg36; reg38=reg45*reg5; reg23=reg11*reg23;
    reg89=reg42+reg89; reg46=reg62+reg46; reg20=reg23-reg20; reg35=reg10+reg35; reg93=reg93/reg57;
    reg82=reg54-reg82; reg9=reg22*reg9; reg10=reg104*reg70; reg11=reg114*reg68; reg22=reg56*reg70;
    reg23=reg71*reg83; reg42=reg6*reg83; reg54=reg2*reg71; reg15=reg112+reg15; reg62=reg94*reg4;
    reg38=reg36+reg38; reg36=reg63*reg68; reg67=reg113*reg70; reg77=reg79*reg68; reg5=reg81*reg5;
    reg0=reg59+reg0; reg59=reg101*reg24; reg84=reg89*reg24; reg85=reg95*reg4; reg52=reg64+reg52;
    reg64=reg46*reg41; reg86=reg61*reg41; reg30=reg80+reg30; reg80=reg113*reg18; reg87=reg114*reg33;
    reg90=reg104*reg18; reg5=reg0+reg5; reg4=reg65*reg4; reg0=reg89*reg46; reg92=reg61*reg101;
    reg84=reg86-reg84; reg86=reg56*reg17; reg96=reg29*reg63; reg113=reg113*reg17; reg97=reg29*reg79;
    reg59=reg64-reg59; reg64=reg73*reg93; reg98=reg47*reg93; reg10=reg11+reg10; reg11=reg15*reg23;
    reg100=reg35*reg23; reg22=reg36+reg22; reg82=reg82/reg57; reg36=reg30*reg23; reg67=reg77+reg67;
    reg85=reg52+reg85; reg114=reg114*reg29; reg104=reg104*reg17; reg9=reg20-reg9; reg79=reg33*reg79;
    reg20=reg14*reg83; reg56=reg56*reg18; reg62=reg38+reg62; reg38=reg2*reg83; reg52=reg6*reg71;
    reg63=reg33*reg63; reg77=reg13*reg71; reg102=reg49*reg93; reg54=reg42+reg54; reg100=reg22+reg100;
    reg22=reg98*reg73; reg80=reg79+reg80; reg86=reg96+reg86; reg42=reg35*reg106; reg79=reg15*reg108;
    reg96=elem.pos(1)[0]-elem.pos(0)[0]; reg103=elem.pos(2)[1]-elem.pos(0)[1]; reg105=elem.pos(2)[0]-elem.pos(0)[0]; reg109=reg64*reg73; reg11=reg10+reg11;
    reg56=reg63+reg56; reg90=reg87+reg90; reg77=reg20+reg77; reg0=reg92-reg0; reg10=reg14*reg71;
    reg20=reg102*reg73; reg36=reg67+reg36; reg63=reg62*reg59; reg67=reg30*reg108; reg35=reg35*reg108;
    reg87=elem.pos(1)[1]-elem.pos(0)[1]; reg38=reg38*reg82; reg52=reg52*reg82; reg92=reg54*reg82; reg104=reg114+reg104;
    reg4=reg5+reg4; reg15=reg15*reg106; reg5=reg85*reg84; reg110=reg13*reg83; reg113=reg97+reg113;
    reg30=reg30*reg106; reg57=reg9/reg57; reg79=reg90+reg79; reg9=reg64*reg25; reg90=reg52*reg54;
    reg97=reg85*reg24; reg111=reg87*reg105; reg22=reg100+reg22; reg100=reg103*reg96; reg8=reg83*reg8;
    reg112=reg46*reg4; reg12=reg71*reg12; reg64=reg64*reg26; reg15=reg104+reg15; reg24=reg24*reg62;
    reg20=reg36+reg20; reg36=reg38*reg54; reg104=reg102*reg26; reg30=reg113+reg30; reg113=reg92*reg54;
    reg114=reg57*reg77; reg10=reg57*reg10; T reg115=reg61*reg85; reg109=reg11+reg109; reg42=reg86+reg42;
    reg11=reg98*reg26; reg110=reg57*reg110; reg98=reg98*reg25; reg86=reg0*reg4; reg61=reg61*reg4;
    reg5=reg63-reg5; reg35=reg56+reg35; reg46=reg46*reg62; reg102=reg102*reg25; reg67=reg80+reg67;
    reg56=reg101*reg62; reg102=reg67+reg102; reg63=reg38*reg12; reg113=reg109+reg113; reg67=reg89*reg4;
    reg111=reg100-reg111; reg89=reg89*reg85; reg112=reg97-reg112; reg80=reg114*reg77; reg36=reg20+reg36;
    reg20=reg110*reg77; reg97=reg10*reg77; reg85=reg85*reg41; reg90=reg22+reg90; reg9=reg79+reg9;
    reg22=reg92*reg12; reg86=reg5+reg86; reg98=reg35+reg98; reg115=reg46-reg115; reg4=reg101*reg4;
    reg5=reg52*reg8; reg11=reg42+reg11; reg75=reg83*reg75; reg62=reg41*reg62; reg61=reg24-reg61;
    reg7=reg71*reg7; reg38=reg38*reg8; reg104=reg30+reg104; reg92=reg92*reg8; reg64=reg15+reg64;
    reg52=reg52*reg12; reg52=reg98+reg52; reg15=reg10*reg7; reg24=1-(*f.m).resolution; reg97=reg90+reg97;
    reg61=reg61/reg86; reg63=reg102+reg63; reg89=reg56-reg89; reg80=reg113+reg80; reg38=reg104+reg38;
    reg30=reg110*reg75; reg5=reg11+reg5; reg110=reg110*reg7; reg115=reg115/reg86; reg11=reg114*reg7;
    reg10=reg10*reg75; reg105=reg105/reg111; reg22=reg9+reg22; reg103=reg103/reg111; reg4=reg85-reg4;
    reg114=reg114*reg75; reg92=reg64+reg92; reg112=reg112/reg86; reg87=reg87/reg111; reg96=reg96/reg111;
    reg67=reg62-reg67; reg20=reg36+reg20; reg89=reg89/reg86; reg9=reg105-reg96; reg0=reg0/reg86;
    reg35=reg87-reg103; reg15=reg52+reg15; reg110=reg63+reg110; reg36=(*f.m).resolution*reg115; reg41=(*f.m).resolution*reg61;
    reg42=(*f.m).resolution*reg112; reg30=reg38+reg30; reg11=reg22+reg11; reg10=reg5+reg10; reg59=reg59/reg86;
    reg4=reg4/reg86; reg114=reg92+reg114; reg84=reg84/reg86; reg86=reg67/reg86; reg5=reg24*reg80;
    reg22=reg24*reg97; reg38=reg24*reg20; reg46=0.5*reg35; reg52=0.5*reg87; reg56=0.5*reg9;
    reg62=0.5*reg105; reg63=0.5*reg103; reg42=reg38+reg42; reg41=reg22-reg41; reg36=reg5+reg36;
    reg5=(*f.m).resolution*reg4; reg22=(*f.m).resolution*reg86; reg38=(*f.m).resolution*reg89; reg64=0.5*reg96; reg67=reg24*reg11;
    reg79=reg15*reg24; reg85=reg110*reg24; reg90=reg114*reg24; reg92=reg10*reg24; reg98=reg30*reg24;
    reg100=(*f.m).resolution*reg0; reg101=(*f.m).resolution*reg59; reg102=(*f.m).resolution*reg84; reg104=reg64*reg36; reg109=reg62*reg36;
    reg113=reg103*reg42; T reg116=reg35*reg42; T reg117=reg36*reg56; T reg118=reg46*reg36; T reg119=reg96*reg41;
    T reg120=reg9*reg41; T reg121=reg52*reg36; reg98=reg101+reg98; reg102=reg92-reg102; reg90=reg100+reg90;
    reg5=reg85-reg5; reg22=reg79+reg22; reg38=reg67-reg38; reg67=reg87*reg42; reg79=reg105*reg41;
    reg85=reg36*reg63; reg92=reg62*reg90; reg100=reg103*reg98; reg101=reg105*reg22; T reg122=reg64*reg90;
    T reg123=reg87*reg98; reg118=reg120+reg118; reg120=reg46*reg38; reg117=reg116+reg117; reg116=reg9*reg102;
    T reg124=reg46*reg90; reg85=reg85-reg79; T reg125=reg62*reg38; T reg126=reg103*reg5; reg119=reg119-reg121;
    T reg127=reg9*reg22; T reg128=reg35*reg5; T reg129=reg96*reg22; T reg130=reg52*reg38; T reg131=reg38*reg56;
    T reg132=reg64*reg38; T reg133=reg52*reg90; T reg134=reg90*reg63; T reg135=reg105*reg102; T reg136=reg96*reg102;
    T reg137=reg38*reg63; T reg138=reg35*reg98; T reg139=reg90*reg56; reg113=reg113-reg109; reg104=reg104-reg67;
    T reg140=reg87*reg5; reg127=reg120+reg127; reg104=2*reg104; reg119=2*reg119; reg129=reg129-reg130;
    reg136=reg136-reg133; reg85=2*reg85; reg126=reg126-reg125; reg131=reg128+reg131; reg132=reg132-reg140;
    reg122=reg122-reg123; reg139=reg138+reg139; reg100=reg100-reg92; reg113=2*reg113; reg118=2*reg118;
    reg137=reg137-reg101; reg117=2*reg117; reg124=reg116+reg124; reg134=reg134-reg135; reg116=reg118*reg63;
    reg120=reg124*reg87; reg128=reg64*reg118; reg138=reg64*reg113; T reg141=reg87*reg100; T reg142=reg35*reg124;
    T reg143=reg62*reg118; T reg144=reg105*reg127; T reg145=reg64*reg85; T reg146=reg85*reg56; T reg147=reg35*reg122;
    T reg148=reg87*reg134; T reg149=reg139*reg87; T reg150=reg64*reg117; T reg151=reg122*reg87; T reg152=reg52*reg118;
    T reg153=reg104*reg56; T reg154=reg35*reg139; T reg155=reg9*reg127; T reg156=reg35*reg136; T reg157=reg113*reg63;
    T reg158=reg105*reg126; T reg159=reg119*reg56; T reg160=reg85*reg63; T reg161=reg105*reg137; T reg162=reg117*reg56;
    T reg163=reg105*reg132; T reg164=reg104*reg63; T reg165=reg103*reg134; T reg166=reg85*reg62; T reg167=reg52*reg85;
    T reg168=reg96*reg137; reg122=reg103*reg122; T reg169=reg104*reg62; T reg170=reg46*reg104; T reg171=reg9*reg126;
    T reg172=reg103*reg136; T reg173=reg62*reg119; T reg174=reg9*reg132; T reg175=reg46*reg113; T reg176=reg117*reg63;
    T reg177=reg118*reg56; T reg178=reg52*reg113; reg126=reg96*reg126; T reg179=reg35*reg100; T reg180=reg113*reg56;
    reg134=reg35*reg134; reg137=reg9*reg137; T reg181=reg64*reg119; reg136=reg136*reg87; reg127=reg96*reg127;
    reg85=reg46*reg85; T reg182=reg52*reg117; T reg183=reg96*reg131; T reg184=reg64*reg104; reg139=reg103*reg139;
    T reg185=reg62*reg117; T reg186=reg52*reg119; reg117=reg46*reg117; T reg187=reg9*reg131; reg124=reg103*reg124;
    T reg188=reg96*reg129; T reg189=reg46*reg119; T reg190=reg9*reg129; reg118=reg46*reg118; reg100=reg103*reg100;
    reg131=reg105*reg131; reg104=reg52*reg104; reg113=reg62*reg113; reg132=reg96*reg132; reg129=reg105*reg129;
    reg119=reg63*reg119; reg186=reg188-reg186; reg152=reg127-reg152; reg182=reg183-reg182; reg141=reg138-reg141;
    reg149=reg150-reg149; reg178=reg126-reg178; reg144=reg116-reg144; reg134=reg146+reg134; reg104=reg132-reg104;
    reg148=reg145-reg148; reg151=reg184-reg151; reg170=reg174+reg170; reg167=reg168-reg167; reg131=reg176-reg131;
    reg120=reg128-reg120; reg189=reg190+reg189; reg185=reg139-reg185; reg117=reg187+reg117; reg143=reg124-reg143;
    reg158=reg157-reg158; reg159=reg156+reg159; reg113=reg100-reg113; reg118=reg155+reg118; reg161=reg160-reg161;
    reg166=reg165-reg166; reg153=reg147+reg153; reg85=reg137+reg85; reg136=reg181-reg136; reg180=reg179+reg180;
    reg177=reg142+reg177; reg175=reg171+reg175; reg162=reg154+reg162; reg129=reg119-reg129; reg173=reg172-reg173;
    reg169=reg122-reg169; reg163=reg164-reg163; reg170=reg111*reg170; reg129=reg111*reg129; reg159=reg111*reg159;
    reg175=reg111*reg175; reg117=reg111*reg117; reg149=reg111*reg149; reg118=reg111*reg118; reg85=reg111*reg85;
    reg158=reg111*reg158; reg144=reg111*reg144; reg131=reg111*reg131; reg189=reg111*reg189; reg185=reg111*reg185;
    reg143=reg111*reg143; reg113=reg111*reg113; reg166=reg111*reg166; reg169=reg111*reg169; reg173=reg111*reg173;
    reg177=reg111*reg177; reg180=reg111*reg180; reg136=reg111*reg136; reg182=reg111*reg182; reg152=reg111*reg152;
    reg178=reg111*reg178; reg167=reg111*reg167; reg104=reg111*reg104; reg186=reg111*reg186; reg120=reg111*reg120;
    reg161=reg111*reg161; reg151=reg111*reg151; reg163=reg111*reg163; reg134=reg111*reg134; reg148=reg111*reg148;
    reg162=reg111*reg162; reg141=reg111*reg141; reg153=reg111*reg153; T tmp_0_2=ponderation*reg180; T tmp_4_4=ponderation*reg151;
    T tmp_4_3=ponderation*reg148; T tmp_4_5=ponderation*reg136; T tmp_1_3=ponderation*reg85; T tmp_5_0=ponderation*reg182; T tmp_4_2=ponderation*reg141;
    T tmp_5_1=ponderation*reg152; T tmp_5_2=ponderation*reg178; T tmp_4_0=ponderation*reg149; T tmp_5_3=ponderation*reg167; T tmp_1_4=ponderation*reg170;
    T tmp_4_1=ponderation*reg120; T tmp_5_4=ponderation*reg104; T tmp_5_5=ponderation*reg186; T tmp_3_4=ponderation*reg163; T tmp_3_2=ponderation*reg158;
    T tmp_3_1=ponderation*reg144; T tmp_0_5=ponderation*reg159; T tmp_3_0=ponderation*reg131; T tmp_3_3=ponderation*reg161; T tmp_1_5=ponderation*reg189;
    T tmp_0_0=ponderation*reg162; T tmp_2_0=ponderation*reg185; T tmp_1_0=ponderation*reg117; T tmp_2_1=ponderation*reg143; T tmp_0_4=ponderation*reg153;
    T tmp_2_2=ponderation*reg113; T tmp_1_1=ponderation*reg118; T tmp_2_3=ponderation*reg166; T tmp_0_3=ponderation*reg134; T tmp_2_4=ponderation*reg169;
    T tmp_2_5=ponderation*reg173; T tmp_1_2=ponderation*reg175; T tmp_0_1=ponderation*reg177; T tmp_3_5=ponderation*reg129;
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
    T reg0=pow((*f.m).v1[0],2); T reg1=pow((*f.m).v1[1],2); T reg2=pow((*f.m).v2[1],2); T reg3=pow((*f.m).v2[0],2); T reg4=pow((*f.m).v1[2],2);
    reg1=reg0+reg1; reg0=2*(*f.m).shear_modulus_13; T reg5=2*(*f.m).shear_modulus_23; T reg6=pow((*f.m).v2[2],2); reg2=reg3+reg2;
    reg4=reg1+reg4; reg6=reg2+reg6; reg5=1.0/reg5; reg0=1.0/reg0; reg4=pow(reg4,0.5);
    reg1=2*(*f.m).shear_modulus_12; reg1=1.0/reg1; reg6=pow(reg6,0.5); reg2=reg0*reg5; reg3=(*f.m).v1[0]/reg4;
    T reg7=(*f.m).v1[1]/reg4; T reg8=2*reg7; T reg9=2*reg3; T reg10=(*f.m).v2[1]/reg6; T reg11=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg12=1.0/(*f.m).elastic_modulus_3; T reg13=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg14=(*f.m).v2[0]/reg6; T reg15=reg1*reg2; reg4=(*f.m).v1[2]/reg4;
    T reg16=pow(reg10,2); T reg17=pow(reg14,2); reg6=(*f.m).v2[2]/reg6; T reg18=reg11*reg15; T reg19=reg13*reg15;
    T reg20=reg12*reg15; T reg21=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg22=1.0/(*f.m).elastic_modulus_2; T reg23=1.0/(*f.m).elastic_modulus_1; T reg24=reg14*reg9;
    T reg25=2*reg4; T reg26=reg10*reg8; T reg27=reg22*reg20; T reg28=reg13*reg19; reg25=reg6*reg25;
    T reg29=reg26*reg22; T reg30=reg24*reg21; T reg31=reg13*reg18; T reg32=reg26*reg21; T reg33=reg24*reg23;
    T reg34=reg16*reg21; T reg35=pow(reg6,2); T reg36=reg17*reg23; T reg37=reg21*reg20; T reg38=reg16*reg22;
    T reg39=reg17*reg21; T reg40=reg22*reg18; T reg41=reg21*reg19; T reg42=reg25*reg13; reg39=reg38-reg39;
    reg38=reg35*reg13; T reg43=reg35*reg11; reg34=reg36-reg34; reg36=reg25*reg11; reg32=reg33-reg32;
    reg33=reg26*reg13; T reg44=reg24*reg11; T reg45=reg16*reg13; T reg46=reg17*reg11; reg30=reg29-reg30;
    reg29=pow(reg7,2); T reg47=pow(reg3,2); reg31=reg37+reg31; reg28=reg27-reg28; reg27=reg23*reg28;
    reg33=reg44+reg33; reg25=reg25*reg12; reg44=reg47*reg21; T reg48=reg21*reg31; reg45=reg46+reg45;
    reg46=reg35*reg12; T reg49=reg29*reg22; T reg50=reg41+reg40; reg42=reg30-reg42; reg30=reg7*reg14;
    T reg51=pow(reg4,2); reg38=reg39-reg38; reg39=reg3*reg10; T reg52=reg47*reg23; T reg53=reg7*reg10;
    T reg54=reg3*reg14; reg36=reg32-reg36; reg43=reg34-reg43; reg32=reg29*reg21; reg34=reg29*reg42;
    T reg55=reg4*reg6; T reg56=reg30+reg39; T reg57=2*reg14; T reg58=reg14*reg10; T reg59=reg38*reg29;
    T reg60=reg43*reg47; reg33=reg25-reg33; reg25=reg43*reg17; T reg61=reg38*reg16; reg45=reg46-reg45;
    reg46=reg29*reg13; T reg62=reg47*reg11; T reg63=reg11*reg2; T reg64=reg1*reg5; T reg65=reg36*reg47;
    reg48=reg27-reg48; reg27=reg11*reg50; T reg66=reg22*reg15; T reg67=reg11*reg19; T reg68=reg13*reg2;
    reg20=reg23*reg20; T reg69=reg11*reg18; reg15=reg21*reg15; reg2=reg12*reg2; T reg70=reg16*reg42;
    T reg71=reg36*reg17; T reg72=reg7*reg6; T reg73=reg4*reg10; T reg74=reg4*reg14; T reg75=reg3*reg6;
    T reg76=reg51*reg13; reg38=reg38*reg53; reg44=reg49-reg44; reg32=reg52-reg32; reg43=reg43*reg54;
    reg36=reg36*reg54; reg42=reg53*reg42; reg49=reg51*reg11; reg59=reg60+reg59; reg19=reg23*reg19;
    reg69=reg20-reg69; reg18=reg21*reg18; reg20=reg55*reg33; reg42=reg36+reg42; reg36=reg10*reg57;
    reg52=reg7*reg9; reg60=reg14*reg6; reg67=reg37+reg67; reg37=reg75+reg74; T reg77=reg11*reg66;
    T reg78=reg51*reg33; reg34=reg65+reg34; reg38=reg43+reg38; reg27=reg48-reg27; reg43=reg55*reg45;
    reg48=reg51*reg45; reg65=reg1*reg0; T reg79=reg64*reg11; reg49=reg32-reg49; reg32=reg64*reg13;
    T reg80=reg21*reg2; reg45=reg35*reg45; T reg81=reg13*reg68; T reg82=reg72-reg73; reg76=reg44-reg76;
    reg44=reg51*reg12; reg46=reg62+reg46; reg61=reg25+reg61; reg25=reg63*reg13; reg2=reg22*reg2;
    reg64=reg64*reg12; reg70=reg71+reg70; reg62=reg58*reg1; reg71=reg56*reg1; reg33=reg35*reg33;
    T reg83=reg11*reg15; reg28=reg28/reg27; T reg84=reg56*reg71; reg43=reg38+reg43; reg38=reg36*reg71;
    reg75=reg74-reg75; reg33=reg70+reg33; reg20=reg42+reg20; reg42=reg56*reg62; reg69=reg69/reg27;
    reg67=reg67/reg27; reg77=reg41+reg77; reg70=reg65*reg13; reg63=reg63*reg22; reg74=reg65*reg11;
    reg68=reg21*reg68; reg45=reg61+reg45; reg46=reg44-reg46; reg81=reg2-reg81; reg25=reg80+reg25;
    reg2=reg64*reg22; reg44=reg60*reg0; reg61=reg37*reg0; reg80=reg49*reg47; T reg85=reg76*reg29;
    reg64=reg64*reg21; T reg86=reg32*reg13; T reg87=reg76*reg16; T reg88=reg49*reg17; T reg89=reg3*reg7;
    T reg90=reg79*reg13; T reg91=2*reg82; T reg92=2*reg10; T reg93=reg6*reg57; reg12=reg65*reg12;
    reg65=reg4*reg9; T reg94=reg10*reg6; reg72=reg73+reg72; reg73=reg36*reg62; reg71=reg52*reg71;
    reg78=reg34+reg78; reg31=reg31/reg27; reg62=reg52*reg62; reg59=reg48+reg59; reg66=reg23*reg66;
    reg18=reg19+reg18; reg15=reg21*reg15; reg83=reg19+reg83; reg19=reg65*reg44; reg62=reg59+reg62;
    reg49=reg49*reg54; reg34=reg29*reg67; reg76=reg76*reg53; reg48=reg16*reg69; reg59=reg52*reg67;
    T reg95=reg36*reg69; T reg96=reg68+reg63; reg77=reg77/reg27; reg90=reg64+reg90; reg15=reg66-reg15;
    reg64=reg47*reg28; reg66=reg29*reg28; T reg97=reg16*reg31; T reg98=reg35*reg46; reg86=reg2-reg86;
    reg87=reg88+reg87; reg71=reg78+reg71; reg2=reg65*reg61; reg78=reg52*reg28; reg88=reg36*reg31;
    T reg99=reg51*reg46; reg85=reg80+reg85; reg80=reg72*reg5; T reg100=reg94*reg5; T reg101=reg17*reg31;
    T reg102=reg47*reg67; reg1=reg89*reg1; T reg103=reg17*reg69; T reg104=reg12*reg22; T reg105=reg37*reg44;
    reg42=reg43+reg42; reg43=reg4*reg8; reg73=reg45+reg73; reg81=reg81*reg23; reg45=reg6*reg92;
    reg44=reg93*reg44; reg32=reg32*reg21; T reg106=reg74*reg13; reg38=reg33+reg38; reg33=reg93*reg61;
    T reg107=pow(reg82,2); T reg108=reg75*reg91; reg50=reg50/reg27; reg61=reg37*reg61; reg84=reg20+reg84;
    reg20=pow(reg75,2); reg18=reg18/reg27; T reg109=reg3*reg4; reg25=reg25*reg21; reg79=reg79*reg22;
    reg13=reg70*reg13; reg12=reg12*reg21; reg83=reg83/reg27; reg98=reg87+reg98; reg106=reg12+reg106;
    reg12=reg7*reg4; reg90=reg90*reg21; reg87=reg36*reg1; reg2=reg71+reg2; reg71=reg43*reg80;
    reg13=reg104-reg13; reg104=reg16*reg83; T reg110=reg52*reg77; T reg111=reg29*reg77; T reg112=reg36*reg83;
    T reg113=reg17*reg83; T reg114=reg47*reg77; T reg115=reg108*reg18; reg95=reg59+reg95; reg59=reg20*reg18;
    reg48=reg34+reg48; reg34=reg32+reg79; T reg116=reg107*reg18; reg44=reg73+reg44; reg103=reg102+reg103;
    reg73=reg45*reg100; reg102=reg108*reg50; reg88=reg78+reg88; reg78=reg72*reg80; reg61=reg84+reg61;
    reg15=reg15/reg27; reg96=reg96*reg11; reg46=reg55*reg46; reg80=reg45*reg80; reg76=reg49+reg76;
    reg33=reg38+reg33; reg25=reg81-reg25; reg19=reg62+reg19; reg38=reg43*reg100; reg0=reg109*reg0;
    reg86=reg86*reg23; reg99=reg85+reg99; reg49=reg52*reg1; reg64=reg101+reg64; reg70=reg70*reg21;
    reg22=reg74*reg22; reg62=reg107*reg50; reg105=reg42+reg105; reg100=reg72*reg100; reg97=reg66+reg97;
    reg42=reg20*reg50; reg66=reg107*reg15; reg100=reg105+reg100; reg104=reg111+reg104; reg74=reg20*reg15;
    reg80=reg33+reg80; reg71=reg2+reg71; reg34=reg34*reg11; reg23=reg13*reg23; reg73=reg44+reg73;
    reg108=reg108*reg15; reg112=reg110+reg112; reg62=reg64+reg62; reg2=reg70+reg22; reg38=reg19+reg38;
    reg42=reg97+reg42; reg21=reg106*reg21; reg90=reg86-reg90; reg13=reg65*reg0; reg49=reg99+reg49;
    reg5=reg12*reg5; reg102=reg88+reg102; reg87=reg98+reg87; reg116=reg103+reg116; reg113=reg114+reg113;
    reg1=reg56*reg1; reg78=reg61+reg78; reg115=reg95+reg115; reg46=reg76+reg46; reg96=reg25-reg96;
    reg59=reg48+reg59; reg19=reg93*reg0; reg96=reg96/reg27; reg25=reg62*reg89; reg33=reg102*reg89;
    reg44=reg115*reg58; reg48=reg116*reg58; reg61=reg42*reg89; reg13=reg49+reg13; reg49=reg45*reg5;
    reg34=reg90-reg34; reg64=reg59*reg58; reg19=reg87+reg19; reg76=reg43*reg5; reg81=reg75*reg82;
    reg0=reg37*reg0; reg1=reg46+reg1; reg11=reg2*reg11; reg2=reg73*reg78; reg66=reg113+reg66;
    reg74=reg104+reg74; reg108=reg112+reg108; reg46=reg80*reg100; reg21=reg23-reg21; reg23=reg7*reg82;
    reg84=reg3*reg75; reg85=reg71*reg100; reg86=reg38*reg78; reg87=reg115*reg16; reg88=reg38*reg80;
    reg90=reg102*reg29; reg46=reg2-reg46; reg85=reg86-reg85; reg49=reg19+reg49; reg2=reg56*reg96;
    reg19=reg29*reg62; reg86=reg116*reg16; reg95=reg47*reg42; reg97=reg59*reg17; reg98=reg108*reg81;
    reg44=reg33+reg44; reg116=reg116*reg17; reg42=reg29*reg42; reg33=reg74*reg81; reg64=reg61+reg64;
    reg59=reg59*reg16; reg62=reg47*reg62; reg61=reg54*reg96; reg99=reg66*reg81; reg48=reg25+reg48;
    reg25=reg53*reg96; reg5=reg72*reg5; reg0=reg1+reg0; reg1=reg14*reg75; reg101=reg10*reg82;
    reg102=reg102*reg47; reg11=reg21-reg11; reg84=reg23+reg84; reg34=reg34/reg27; reg115=reg115*reg17;
    reg21=reg71*reg73; reg76=reg13+reg76; reg13=reg7*reg75; reg23=reg3*reg82; reg103=reg25*reg56;
    reg33=reg64+reg33; reg64=reg66*reg20; reg104=reg2*reg56; reg27=reg11/reg27; reg11=reg74*reg107;
    reg97=reg95+reg97; reg1=reg101+reg1; reg115=reg102+reg115; reg98=reg44+reg98; reg44=reg108*reg107;
    reg23=reg23*reg34; reg5=reg0+reg5; reg0=reg76*reg46; reg13=reg13*reg34; reg95=reg61*reg56;
    reg99=reg48+reg99; reg108=reg108*reg20; reg87=reg90+reg87; reg48=elem.pos(1)[0]-elem.pos(0)[0]; reg90=reg49*reg85;
    reg74=reg74*reg20; reg101=elem.pos(2)[1]-elem.pos(0)[1]; reg66=reg66*reg107; reg116=reg62+reg116; reg59=reg42+reg59;
    reg21=reg88-reg21; reg42=elem.pos(2)[0]-elem.pos(0)[0]; reg62=reg14*reg82; reg88=reg84*reg34; reg102=reg10*reg75;
    reg105=elem.pos(1)[1]-elem.pos(0)[1]; reg86=reg19+reg86; reg19=reg100*reg76; reg103=reg33+reg103; reg33=reg13*reg84;
    reg66=reg116+reg66; reg106=reg61*reg24; reg110=reg25*reg26; reg111=reg2*reg24; reg61=reg61*reg26;
    reg44=reg115+reg44; reg64=reg86+reg64; reg108=reg87+reg108; reg2=reg2*reg26; reg90=reg0-reg90;
    reg0=reg101*reg48; reg86=reg73*reg76; reg87=reg21*reg5; reg112=reg38*reg49; reg113=reg105*reg42;
    reg62=reg27*reg62; reg102=reg27*reg102; reg114=reg27*reg1; reg25=reg25*reg24; reg11=reg97+reg11;
    reg9=reg82*reg9; reg97=reg88*reg84; reg104=reg98+reg104; reg8=reg75*reg8; reg95=reg99+reg95;
    reg98=reg23*reg84; reg74=reg59+reg74; reg100=reg49*reg100; reg38=reg38*reg5; reg73=reg73*reg5;
    reg110=reg74+reg110; reg59=reg13*reg8; reg87=reg90+reg87; reg57=reg82*reg57; reg92=reg75*reg92;
    reg74=reg49*reg78; reg90=reg80*reg5; reg78=reg78*reg76; reg73=reg100-reg73; reg112=reg86-reg112;
    reg113=reg0-reg113; reg5=reg71*reg5; reg2=reg108+reg2; reg0=reg88*reg8; reg97=reg104+reg97;
    reg86=reg114*reg1; reg49=reg71*reg49; reg38=reg19-reg38; reg76=reg80*reg76; reg106=reg66+reg106;
    reg19=reg23*reg9; reg66=reg102*reg1; reg33=reg103+reg33; reg71=reg62*reg1; reg98=reg95+reg98;
    reg25=reg11+reg25; reg13=reg13*reg9; reg23=reg23*reg8; reg61=reg64+reg61; reg111=reg44+reg111;
    reg88=reg88*reg9; reg88=reg111+reg88; reg23=reg61+reg23; reg112=reg112/reg87; reg49=reg76-reg49;
    reg11=reg62*reg92; reg101=reg101/reg113; reg0=reg2+reg0; reg2=reg114*reg92; reg114=reg114*reg57;
    reg59=reg110+reg59; reg44=reg102*reg92; reg61=1-(*f.m).resolution; reg42=reg42/reg113; reg86=reg97+reg86;
    reg19=reg106+reg19; reg62=reg62*reg57; reg38=reg38/reg87; reg5=reg78-reg5; reg48=reg48/reg113;
    reg105=reg105/reg113; reg90=reg74-reg90; reg66=reg33+reg66; reg71=reg98+reg71; reg73=reg73/reg87;
    reg13=reg25+reg13; reg102=reg102*reg57; reg46=reg46/reg87; reg90=reg90/reg87; reg11=reg23+reg11;
    reg23=(*f.m).resolution*reg112; reg102=reg13+reg102; reg62=reg19+reg62; reg114=reg88+reg114; reg44=reg59+reg44;
    reg13=(*f.m).resolution*reg38; reg19=(*f.m).resolution*reg73; reg25=reg61*reg86; reg33=reg61*reg66; reg59=reg61*reg71;
    reg64=reg42-reg48; reg74=reg105-reg101; reg85=reg85/reg87; reg5=reg5/reg87; reg21=reg21/reg87;
    reg2=reg0+reg2; reg87=reg49/reg87; reg0=0.5*reg74; reg49=0.5*reg64; reg76=0.5*reg105;
    reg23=reg25+reg23; reg13=reg33-reg13; reg19=reg59+reg19; reg25=0.5*reg42; reg33=0.5*reg48;
    reg59=(*f.m).resolution*reg87; reg78=(*f.m).resolution*reg5; reg80=(*f.m).resolution*reg90; reg88=reg61*reg2; reg95=reg44*reg61;
    reg97=reg11*reg61; reg98=reg114*reg61; reg99=reg102*reg61; reg100=(*f.m).resolution*reg46; reg103=(*f.m).resolution*reg21;
    reg104=0.5*reg101; reg106=reg62*reg61; reg108=(*f.m).resolution*reg85; reg59=reg88-reg59; reg78=reg95+reg78;
    reg80=reg97-reg80; reg98=reg103+reg98; reg108=reg99-reg108; reg106=reg100+reg106; reg88=reg64*reg13;
    reg95=reg23*reg49; reg97=reg74*reg19; reg99=reg0*reg23; reg100=reg33*reg23; reg103=reg105*reg19;
    reg110=reg48*reg13; reg111=reg76*reg23; reg115=reg25*reg23; reg116=reg101*reg19; T reg117=reg23*reg104;
    T reg118=reg42*reg13; T reg119=reg0*reg59; T reg120=reg59*reg104; T reg121=reg42*reg78; reg110=reg110-reg111;
    T reg122=reg64*reg78; T reg123=reg101*reg80; T reg124=reg0*reg98; T reg125=reg64*reg108; reg95=reg97+reg95;
    reg97=reg98*reg49; T reg126=reg74*reg106; T reg127=reg25*reg59; T reg128=reg76*reg98; T reg129=reg48*reg108;
    T reg130=reg33*reg59; reg117=reg117-reg118; reg100=reg100-reg103; T reg131=reg105*reg80; T reg132=reg105*reg106;
    T reg133=reg33*reg98; T reg134=reg98*reg104; T reg135=reg48*reg78; reg99=reg88+reg99; reg88=reg42*reg108;
    T reg136=reg76*reg59; T reg137=reg101*reg106; T reg138=reg25*reg98; reg116=reg116-reg115; reg130=reg130-reg131;
    reg95=2*reg95; reg116=2*reg116; reg97=reg126+reg97; reg129=reg129-reg128; reg137=reg137-reg138;
    reg117=2*reg117; reg100=2*reg100; reg99=2*reg99; reg133=reg133-reg132; reg120=reg120-reg121;
    reg123=reg123-reg127; reg135=reg135-reg136; reg110=2*reg110; reg134=reg134-reg88; reg124=reg125+reg124;
    reg122=reg119+reg122; reg119=reg0*reg116; reg125=reg117*reg49; reg126=reg64*reg123; T reg139=reg76*reg110;
    T reg140=reg48*reg135; T reg141=reg64*reg120; T reg142=reg0*reg99; T reg143=reg133*reg105; T reg144=reg42*reg135;
    T reg145=reg0*reg117; T reg146=reg104*reg110; T reg147=reg42*reg130; T reg148=reg100*reg104; T reg149=reg42*reg120;
    T reg150=reg110*reg49; T reg151=reg117*reg104; T reg152=reg74*reg137; T reg153=reg64*reg122; T reg154=reg129*reg105;
    T reg155=reg33*reg110; T reg156=reg116*reg49; T reg157=reg74*reg134; T reg158=reg74*reg97; T reg159=reg25*reg116;
    T reg160=reg100*reg25; T reg161=reg101*reg129; T reg162=reg25*reg110; T reg163=reg95*reg49; T reg164=reg101*reg137;
    T reg165=reg33*reg100; reg129=reg74*reg129; reg135=reg64*reg135; reg110=reg0*reg110; T reg166=reg101*reg133;
    T reg167=reg74*reg124; T reg168=reg117*reg25; T reg169=reg100*reg49; T reg170=reg99*reg49; T reg171=reg101*reg134;
    T reg172=reg64*reg130; reg133=reg74*reg133; T reg173=reg0*reg100; reg173=reg172+reg173; reg159=reg164-reg159;
    reg168=reg171-reg168; reg163=reg158+reg163; reg156=reg152+reg156; reg154=reg155-reg154; reg160=reg166-reg160;
    reg150=reg129+reg150; reg170=reg167+reg170; reg139=reg140-reg139; reg157=reg125+reg157; reg169=reg133+reg169;
    reg119=reg126+reg119; reg143=reg165-reg143; reg110=reg135+reg110; reg144=reg146-reg144; reg162=reg161-reg162;
    reg147=reg148-reg147; reg145=reg141+reg145; reg149=reg151-reg149; reg142=reg153+reg142; reg154=reg113*reg154;
    reg110=reg113*reg110; reg150=reg113*reg150; reg119=reg113*reg119; reg142=reg113*reg142; reg139=reg113*reg139;
    reg170=reg113*reg170; reg157=reg113*reg157; reg169=reg113*reg169; reg143=reg113*reg143; reg162=reg113*reg162;
    reg144=reg113*reg144; reg145=reg113*reg145; reg147=reg113*reg147; reg149=reg113*reg149; reg173=reg113*reg173;
    reg160=reg113*reg160; reg156=reg113*reg156; reg163=reg113*reg163; reg168=reg113*reg168; reg159=reg113*reg159;
    T tmp_1_3=ponderation*reg145; T tmp_1_4=ponderation*reg173; T tmp_5_5=ponderation*reg139; T tmp_4_5=ponderation*reg154; T tmp_4_4=ponderation*reg143;
    T tmp_0_2=ponderation*reg156; T tmp_0_1=ponderation*reg170; T tmp_2_5=ponderation*reg162; T tmp_2_4=ponderation*reg160; T tmp_2_3=ponderation*reg168;
    T tmp_2_2=ponderation*reg159; T tmp_1_5=ponderation*reg110; T tmp_0_0=ponderation*reg163; T tmp_0_4=ponderation*reg169; T tmp_0_3=ponderation*reg157;
    T tmp_3_5=ponderation*reg144; T tmp_3_4=ponderation*reg147; T tmp_3_3=ponderation*reg149; T tmp_0_5=ponderation*reg150; T tmp_1_1=ponderation*reg142;
    T tmp_1_2=ponderation*reg119;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; T reg2=2*(*f.m).shear_modulus_12; reg1=1.0/reg1; reg0=1.0/reg0;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=1.0/(*f.m).elastic_modulus_3; T reg5=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg6=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg7=reg2*reg3; T reg8=pow((*f.m).v2[1],2); T reg9=pow((*f.m).v2[0],2); T reg10=pow((*f.m).v1[0],2); T reg11=pow((*f.m).v1[1],2);
    T reg12=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg13=pow((*f.m).v1[2],2); reg8=reg9+reg8; reg9=pow((*f.m).v2[2],2); T reg14=reg6*reg7;
    T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=reg5*reg7; reg11=reg10+reg11; reg10=reg4*reg7; T reg17=reg12*reg10;
    T reg18=reg5*reg16; T reg19=reg5*reg14; T reg20=reg15*reg10; reg9=reg8+reg9; reg13=reg11+reg13;
    reg8=1.0/(*f.m).elastic_modulus_1; reg13=pow(reg13,0.5); reg9=pow(reg9,0.5); reg11=reg15*reg14; T reg21=reg12*reg16;
    reg19=reg17+reg19; reg18=reg20-reg18; reg20=(*f.m).v2[1]/reg9; T reg22=reg12*reg19; T reg23=reg8*reg18;
    T reg24=reg21+reg11; T reg25=(*f.m).v2[2]/reg9; T reg26=(*f.m).v1[2]/reg13; T reg27=(*f.m).v1[1]/reg13; T reg28=reg26*reg20;
    reg9=(*f.m).v2[0]/reg9; T reg29=reg27*reg25; T reg30=reg6*reg3; T reg31=reg2*reg0; reg13=(*f.m).v1[0]/reg13;
    T reg32=reg4*reg3; T reg33=reg12*reg7; T reg34=reg6*reg14; reg10=reg8*reg10; reg3=reg5*reg3;
    reg22=reg23-reg22; reg23=reg6*reg24; T reg35=reg6*reg16; reg7=reg15*reg7; T reg36=reg15*reg32;
    reg23=reg22-reg23; reg22=reg29-reg28; T reg37=reg31*reg4; T reg38=2*reg13; T reg39=reg30*reg5;
    T reg40=2*reg9; T reg41=reg5*reg3; reg32=reg12*reg32; T reg42=reg6*reg33; T reg43=reg31*reg5;
    reg31=reg31*reg6; T reg44=reg2*reg1; reg16=reg8*reg16; T reg45=reg26*reg9; T reg46=reg13*reg25;
    reg14=reg12*reg14; reg34=reg10-reg34; reg35=reg17+reg35; reg10=reg6*reg7; reg30=reg30*reg15;
    reg3=reg12*reg3; reg34=reg34/reg23; reg17=2*reg22; reg10=reg21+reg10; T reg47=reg20*reg40;
    reg35=reg35/reg23; T reg48=pow(reg20,2); T reg49=pow(reg9,2); T reg50=reg27*reg38; reg18=reg18/reg23;
    T reg51=pow(reg27,2); T reg52=reg27*reg9; T reg53=pow(reg13,2); T reg54=reg13*reg20; reg19=reg19/reg23;
    T reg55=reg45-reg46; reg14=reg16+reg14; reg7=reg8*reg7; reg42=reg16+reg42; reg33=reg12*reg33;
    reg16=reg44*reg4; T reg56=reg31*reg5; T reg57=reg43*reg5; T reg58=reg37*reg12; reg37=reg37*reg15;
    reg39=reg32+reg39; reg41=reg36-reg41; reg32=reg44*reg5; reg44=reg44*reg6; reg36=reg48*reg34;
    T reg59=reg51*reg35; T reg60=reg49*reg34; T reg61=reg53*reg35; T reg62=reg47*reg19; T reg63=reg50*reg18;
    reg57=reg37-reg57; reg42=reg42/reg23; reg56=reg58+reg56; reg37=reg54-reg52; reg58=reg16*reg15;
    reg16=reg16*reg12; T reg64=reg32*reg5; reg10=reg10/reg23; T reg65=reg44*reg5; T reg66=reg51*reg18;
    T reg67=reg48*reg19; T reg68=reg53*reg18; reg24=reg24/reg23; reg33=reg7-reg33; reg31=reg31*reg15;
    reg14=reg14/reg23; reg43=reg43*reg12; reg7=pow(reg55,2); T reg69=pow(reg22,2); reg41=reg41*reg8;
    T reg70=pow(reg25,2); T reg71=reg50*reg35; T reg72=reg47*reg34; reg39=reg39*reg12; T reg73=reg55*reg17;
    T reg74=reg49*reg19; T reg75=reg3+reg30; T reg76=pow(reg26,2); reg32=reg32*reg12; T reg77=reg50*reg10;
    T reg78=reg47*reg42; reg64=reg58-reg64; reg65=reg16+reg65; reg75=reg75*reg6; reg67=reg66+reg67;
    reg44=reg44*reg15; reg16=reg7*reg24; reg39=reg41-reg39; reg33=reg33/reg23; reg41=pow(reg37,2);
    reg58=reg7*reg14; reg36=reg59+reg36; reg59=reg76*reg35; reg66=reg70*reg34; T reg79=reg69*reg14;
    reg60=reg61+reg60; reg61=reg73*reg24; reg62=reg63+reg62; reg72=reg71+reg72; reg63=reg73*reg14;
    reg71=reg70*reg19; T reg80=reg76*reg18; T reg81=reg69*reg24; reg68=reg74+reg68; reg57=reg57*reg8;
    reg74=reg53*reg10; T reg82=reg49*reg42; reg56=reg56*reg12; T reg83=reg51*reg10; T reg84=reg48*reg42;
    T reg85=reg43+reg31; T reg86=reg69*reg33; reg65=reg65*reg12; T reg87=reg32+reg44; reg81=reg68+reg81;
    reg16=reg67+reg16; reg82=reg74+reg82; reg78=reg77+reg78; reg73=reg73*reg33; reg61=reg62+reg61;
    reg62=2*reg20; reg67=reg25*reg40; reg63=reg72+reg63; reg68=reg41*reg24; reg71=reg80+reg71;
    reg72=2*reg27; reg74=reg9*reg20; reg77=reg13*reg27; reg80=reg26*reg38; T reg88=reg41*reg14;
    reg66=reg59+reg66; reg59=reg7*reg33; reg84=reg83+reg84; reg79=reg60+reg79; reg60=reg70*reg42;
    reg58=reg36+reg58; reg75=reg39-reg75; reg56=reg57-reg56; reg36=reg76*reg10; reg85=reg85*reg6;
    reg64=reg64*reg8; reg39=reg51*reg16; reg57=reg41*reg33; reg83=reg9*reg38; reg88=reg66+reg88;
    reg86=reg82+reg86; reg65=reg64-reg65; reg64=reg61*reg51; reg66=reg58*reg48; reg82=reg26*reg72;
    T reg89=reg61*reg53; T reg90=reg79*reg48; T reg91=reg51*reg81; reg75=reg75/reg23; T reg92=reg80*reg35;
    T reg93=reg67*reg34; T reg94=reg63*reg48; T reg95=2*reg55; T reg96=reg16*reg77; reg73=reg78+reg73;
    reg87=reg87*reg6; reg85=reg56-reg85; reg17=reg37*reg17; reg56=reg58*reg74; reg78=reg63*reg49;
    T reg97=reg13*reg9; T reg98=reg27*reg22; reg59=reg84+reg59; reg84=reg13*reg55; T reg99=reg27*reg20;
    reg58=reg58*reg49; reg16=reg53*reg16; reg61=reg61*reg77; T reg100=reg67*reg19; T reg101=reg80*reg18;
    reg63=reg63*reg74; T reg102=reg25*reg62; reg68=reg71+reg68; reg60=reg36+reg60; reg36=reg20*reg72;
    reg71=2*reg26; reg54=reg52+reg54; reg52=reg53*reg81; T reg103=reg55*reg22; T reg104=reg79*reg49;
    T reg105=reg86*reg69; T reg106=reg54*reg75; T reg107=reg36*reg15; reg18=reg82*reg18; reg19=reg102*reg19;
    T reg108=reg27*reg55; reg95=reg37*reg95; T reg109=reg26*reg25; T reg110=reg13*reg22; reg104=reg52+reg104;
    reg52=reg49*reg8; T reg111=reg59*reg69; reg58=reg16+reg58; reg16=reg49*reg12; T reg112=reg59*reg103;
    T reg113=reg97*reg75; T reg114=reg86*reg7; reg56=reg96+reg56; reg94=reg64+reg94; reg90=reg91+reg90;
    reg64=reg73*reg7; reg91=reg99*reg75; reg85=reg85/reg23; reg96=reg67*reg42; T reg115=reg36*reg12;
    reg71=reg25*reg71; T reg116=reg9*reg55; T reg117=reg20*reg22; T reg118=reg83*reg8; T reg119=reg73*reg103;
    reg57=reg60+reg57; reg60=reg83*reg12; T reg120=reg88*reg48; T reg121=reg68*reg51; reg63=reg61+reg63;
    reg61=reg80*reg10; reg93=reg92+reg93; reg92=reg48*reg15; reg87=reg65-reg87; reg65=reg68*reg53;
    T reg122=reg88*reg49; reg66=reg39+reg66; reg73=reg73*reg69; reg78=reg89+reg78; reg84=reg98+reg84;
    reg59=reg59*reg7; reg39=reg17*reg24; reg89=reg48*reg12; reg79=reg79*reg74; reg81=reg81*reg77;
    reg100=reg101+reg100; reg34=reg102*reg34; reg35=reg82*reg35; reg98=reg17*reg14; reg68=reg68*reg77;
    reg101=reg91*reg54; reg112=reg56+reg112; reg86=reg86*reg103; reg79=reg81+reg79; reg56=reg91*reg36;
    reg120=reg121+reg120; reg81=reg57*reg7; reg96=reg61+reg96; reg17=reg17*reg33; reg10=reg82*reg10;
    reg42=reg102*reg42; reg61=reg26*reg37; reg121=reg84*reg85; reg108=reg108*reg85; reg110=reg110*reg85;
    T reg123=reg109*reg75; T reg124=reg36*reg5; T reg125=reg70*reg6; reg89=reg52-reg89; reg52=reg9*reg22;
    T reg126=reg20*reg55; reg116=reg117+reg116; reg38=reg22*reg38; reg72=reg55*reg72; reg117=reg71*reg5;
    reg16=reg92-reg16; reg92=reg70*reg5; T reg127=reg57*reg69; reg122=reg65+reg122; reg39=reg100+reg39;
    reg19=reg18+reg19; reg91=reg91*reg83; reg111=reg58+reg111; reg24=reg95*reg24; reg18=reg71*reg6;
    reg58=reg113*reg83; reg105=reg104+reg105; reg115=reg118-reg115; reg65=reg113*reg36; reg114=reg90+reg114;
    reg64=reg94+reg64; reg90=reg106*reg36; reg59=reg66+reg59; reg23=reg87/reg23; reg66=reg106*reg54;
    reg106=reg106*reg83; reg119=reg63+reg119; reg73=reg78+reg73; reg88=reg88*reg74; reg60=reg107-reg60;
    reg14=reg95*reg14; reg34=reg35+reg34; reg35=reg49*reg6; reg63=reg48*reg5; reg98=reg93+reg98;
    reg78=reg83*reg6; reg15=reg51*reg15; reg87=reg53*reg12; reg12=reg51*reg12; reg24=reg19+reg24;
    reg18=reg115-reg18; reg90=reg64+reg90; reg19=reg121*reg72; reg33=reg95*reg33; reg42=reg10+reg42;
    reg10=reg39*reg51; reg17=reg96+reg17; reg92=reg16-reg92; reg16=reg98*reg48; reg14=reg34+reg14;
    reg125=reg89-reg125; reg88=reg68+reg88; reg57=reg57*reg103; reg34=reg49*(*f.m).alpha_2; reg127=reg122+reg127;
    reg64=reg110*reg72; reg65=reg114+reg65; reg52=reg23*reg52; reg126=reg23*reg126; reg68=reg98*reg49;
    reg89=reg39*reg53; reg93=reg121*reg38; reg106=reg73+reg106; reg121=reg121*reg84; reg66=reg119+reg66;
    reg73=reg123*reg83; reg117=reg60-reg117; reg62=reg55*reg62; reg40=reg22*reg40; reg60=reg70*reg4;
    reg63=reg35+reg63; reg71=reg71*reg4; reg124=reg78+reg124; reg35=(*f.m).alpha_1*reg53; reg78=reg25*reg37;
    reg61=reg61*reg85; reg45=reg46+reg45; reg46=reg26*reg22; reg94=reg13*reg37; reg95=reg23*reg116;
    reg96=reg48*(*f.m).alpha_2; reg8=reg53*reg8; reg100=reg123*reg36; reg81=reg120+reg81; reg58=reg105+reg58;
    reg104=reg110*reg38; reg105=reg108*reg72; reg56=reg59+reg56; reg86=reg79+reg86; reg59=reg108*reg38;
    reg91=reg111+reg91; reg108=reg108*reg84; reg101=reg112+reg101; reg113=reg113*reg54; reg79=(*f.m).alpha_1*reg51;
    reg107=reg25*reg22; reg39=reg39*reg77; reg111=(*f.m).alpha_3*reg7; reg34=reg35+reg34; reg35=reg9*reg37;
    reg112=reg18*reg53; reg114=reg95*reg116; reg96=reg79+reg96; reg79=reg51*reg117; reg115=reg61*reg38;
    reg121=reg66+reg121; reg127=reg73+reg127; reg66=reg125*reg49; reg73=(*f.m).alpha_3*reg69; reg118=reg92*reg48;
    reg98=reg98*reg74; reg12=reg8-reg12; reg64=reg65+reg64; reg8=reg52*reg62; reg123=reg123*reg54;
    reg57=reg88+reg57; reg65=reg126*reg116; reg108=reg101+reg108; reg88=reg76*reg6; reg110=reg110*reg84;
    reg113=reg86+reg113; reg105=reg56+reg105; reg56=reg126*reg62; reg100=reg81+reg100; reg81=reg61*reg72;
    reg86=reg26*reg55; reg94=reg46+reg94; reg29=reg28+reg29; reg28=reg27*reg37; reg46=reg92*reg51;
    reg101=reg125*reg53; reg124=reg71-reg124; reg63=reg60-reg63; reg60=reg51*reg5; reg6=reg53*reg6;
    reg71=reg99*reg117; reg119=reg18*reg97; reg92=reg92*reg99; reg125=reg125*reg97; reg93=reg106+reg93;
    reg106=reg95*reg40; reg68=reg89+reg68; reg89=reg17*reg69; reg120=reg24*reg53; reg122=reg14*reg49;
    reg126=reg126*reg40; reg59=reg91+reg59; reg91=reg52*reg40; T reg128=(*f.m).alpha_1*reg76; T reg129=reg70*(*f.m).alpha_2;
    T reg130=reg45*reg75; reg33=reg42+reg33; reg18=reg18*reg49; reg117=reg48*reg117; reg104=reg58+reg104;
    reg19=reg90+reg19; reg95=reg95*reg62; reg16=reg10+reg16; reg10=reg17*reg7; reg78=reg23*reg78;
    reg87=reg15-reg87; reg5=reg76*reg5; reg15=reg14*reg48; reg42=reg24*reg51; reg58=(*f.m).alpha_1*reg77;
    reg41=(*f.m).alpha_3*reg41; reg89=reg68+reg89; reg68=reg109*reg124; reg129=reg128+reg129; reg71=reg119+reg71;
    reg106=reg93+reg106; reg90=reg109*reg63; reg111=reg96+reg111; reg92=reg125+reg92; reg114=reg121+reg114;
    reg79=reg112+reg79; reg93=reg76*reg124; reg98=reg39+reg98; reg17=reg17*reg103; reg118=reg66+reg118;
    reg56=reg105+reg56; reg24=reg24*reg77; reg14=reg14*reg74; reg39=reg70*reg63; reg60=reg6+reg60;
    reg6=reg74*reg2; reg66=reg54*reg2; reg8=reg64+reg8; reg46=reg101+reg46; reg61=reg61*reg84;
    reg123=reg57+reg123; reg28=reg86+reg28; reg4=reg76*reg4; reg65=reg108+reg65; reg69=reg33*reg69;
    reg122=reg120+reg122; reg35=reg107+reg35; reg57=reg25*reg55; reg52=reg52*reg116; reg64=reg130*reg83;
    reg110=reg113+reg110; reg86=reg20*reg37; reg74=reg74*(*f.m).alpha_2; reg13=reg13*reg26; reg9=reg9*reg25;
    reg96=reg78*reg40; reg115=reg127+reg115; reg124=reg70*reg124; reg81=reg100+reg81; reg117=reg18+reg117;
    reg73=reg34+reg73; reg7=reg33*reg7; reg75=reg29*reg75; reg94=reg94*reg85; reg91=reg104+reg91;
    reg63=reg76*reg63; reg126=reg59+reg126; reg88=reg12-reg88; reg95=reg19+reg95; reg10=reg16+reg10;
    reg5=reg87-reg5; reg12=reg130*reg36; reg15=reg42+reg15; reg16=reg78*reg62; reg28=reg85*reg28;
    reg46=reg63+reg46; reg41=reg129+reg41; reg17=reg98+reg17; reg124=reg117+reg124; reg25=reg20*reg25;
    reg26=reg27*reg26; reg22=reg37*reg22; reg18=reg47*reg66; reg19=reg126*reg111; reg61=reg123+reg61;
    reg16=reg81+reg16; reg20=reg8*reg73; reg27=reg56*reg111; reg78=reg78*reg116; reg96=reg115+reg96;
    reg52=reg110+reg52; reg12=reg10+reg12; reg10=reg91*reg73; reg86=reg57+reg86; reg34=reg5*reg51;
    reg42=reg88*reg53; reg57=reg94*reg72; reg59=reg9*(*f.m).alpha_2; reg63=(*f.m).alpha_1*reg13; reg81=reg45*reg1;
    reg85=(*f.m).alpha_3*reg103; reg74=reg58+reg74; reg9=reg9*reg1; reg58=reg106*reg65; reg87=reg95*reg65;
    reg39=reg118+reg39; reg98=reg47*reg6; reg100=reg126*reg114; reg101=reg54*reg6; reg68=reg71+reg68;
    reg71=reg54*reg66; reg104=reg5*reg48; reg105=reg88*reg49; reg36=reg75*reg36; reg107=reg56*reg114;
    reg90=reg92+reg90; reg64=reg89+reg64; reg130=reg130*reg54; reg6=reg50*reg6; reg60=reg4-reg60;
    reg93=reg79+reg93; reg83=reg75*reg83; reg69=reg122+reg69; reg7=reg15+reg7; reg66=reg50*reg66;
    reg14=reg24+reg14; reg103=reg33*reg103; reg4=reg94*reg38; reg35=reg23*reg35; reg5=reg5*reg99;
    reg101=reg90+reg101; reg78=reg61+reg78; reg15=reg45*reg9; reg88=reg88*reg97; reg19=reg10+reg19;
    reg87=reg107-reg87; reg2=reg77*reg2; reg38=reg28*reg38; reg83=reg69+reg83; reg10=reg35*reg40;
    reg4=reg64+reg4; reg24=reg35*reg62; reg33=reg45*reg81; reg71=reg68+reg71; reg57=reg12+reg57;
    reg36=reg7+reg36; reg72=reg28*reg72; reg7=reg106*reg56; reg12=reg126*reg95; reg61=reg25*reg0;
    reg64=reg29*reg0; reg58=reg100-reg58; reg34=reg42+reg34; reg76=reg76*reg60; reg42=(*f.m).alpha_1*reg26;
    reg22=(*f.m).alpha_3*reg22; reg59=reg63+reg59; reg85=reg74+reg85; reg63=reg67*reg81; reg98=reg39+reg98;
    reg18=reg124+reg18; reg55=reg37*reg55; reg130=reg17+reg130; reg94=reg94*reg84; reg17=reg67*reg9;
    reg37=reg52*reg73; reg39=reg111*reg65; reg75=reg75*reg54; reg103=reg14+reg103; reg9=reg80*reg9;
    reg6=reg46+reg6; reg96=reg96*reg41; reg66=reg93+reg66; reg81=reg80*reg81; reg27=reg20+reg27;
    reg16=reg41*reg16; reg86=reg23*reg86; reg104=reg105+reg104; reg70=reg70*reg60; reg25=reg25*(*f.m).alpha_2;
    reg14=reg82*reg61; reg5=reg88+reg5; reg60=reg109*reg60; reg75=reg103+reg75; reg84=reg28*reg84;
    reg20=reg102*reg61; reg17=reg98+reg17; reg72=reg36+reg72; reg15=reg101+reg15; reg81=reg66+reg81;
    reg7=reg12-reg7; reg62=reg86*reg62; reg12=reg8*reg58; reg23=reg91*reg87; reg61=reg29*reg61;
    reg28=reg82*reg64; reg9=reg6+reg9; reg6=reg47*reg2; reg70=reg104+reg70; reg96=reg19+reg96;
    reg19=reg106*reg85; reg55=(*f.m).alpha_3*reg55; reg25=reg42+reg25; reg16=reg27+reg16; reg27=reg102*reg64;
    reg36=reg50*reg2; reg76=reg34+reg76; reg34=reg85*reg95; reg63=reg18+reg63; reg22=reg59+reg22;
    reg1=reg13*reg1; reg24=reg57+reg24; reg40=reg86*reg40; reg38=reg83+reg38; reg94=reg130+reg94;
    reg10=reg4+reg10; reg39=reg37+reg39; reg78=reg41*reg78; reg35=reg35*reg116; reg33=reg71+reg33;
    reg64=reg29*reg64; reg14=reg9+reg14; reg55=reg25+reg55; reg60=reg5+reg60; reg2=reg54*reg2;
    reg64=reg33+reg64; reg35=reg94+reg35; reg28=reg81+reg28; reg84=reg75+reg84; reg80=reg80*reg1;
    reg36=reg76+reg36; reg116=reg86*reg116; reg6=reg70+reg6; reg20=reg17+reg20; reg67=reg67*reg1;
    reg0=reg26*reg0; reg61=reg15+reg61; reg4=reg91*reg114; reg40=reg38+reg40; reg5=reg106*reg52;
    reg9=reg95*reg52; reg13=reg8*reg114; reg15=reg7*reg52; reg12=reg23-reg12; reg19=reg96+reg19;
    reg10=reg10*reg22; reg34=reg16+reg34; reg24=reg22*reg24; reg62=reg72+reg62; reg27=reg63+reg27;
    reg16=reg85*reg114; reg78=reg39+reg78; reg16=reg78+reg16; reg35=reg22*reg35; reg17=reg91*reg65;
    reg18=reg20*reg64; reg5=reg4-reg5; reg4=reg14*reg64; reg22=reg126*reg52; reg23=reg28*reg61;
    reg25=reg91*reg95; reg116=reg84+reg116; reg26=reg106*reg8; reg80=reg36+reg80; reg82=reg82*reg0;
    reg33=reg27*reg61; reg10=reg19+reg10; reg40=reg40*reg55; reg62=reg55*reg62; reg24=reg34+reg24;
    reg67=reg6+reg67; reg15=reg12+reg15; reg102=reg102*reg0; reg6=reg8*reg65; reg9=reg13-reg9;
    reg1=reg45*reg1; reg2=reg60+reg2; reg12=reg56*reg52; reg82=reg80+reg82; reg26=reg25-reg26;
    reg13=reg56*reg91; reg19=reg8*reg126; reg35=reg16+reg35; reg22=reg17-reg22; reg87=reg87/reg15;
    reg16=reg28*reg20; reg9=reg9/reg15; reg5=reg5/reg15; reg116=reg55*reg116; reg17=reg14*reg27;
    reg102=reg67+reg102; reg12=reg6-reg12; reg1=reg2+reg1; reg33=reg18-reg33; reg62=reg24+reg62;
    reg0=reg29*reg0; reg40=reg10+reg40; reg23=reg4-reg23; reg58=reg58/reg15; reg0=reg1+reg0;
    reg19=reg13-reg19; reg16=reg17-reg16; reg9=reg9*reg40; reg1=reg102*reg23; reg2=reg82*reg33;
    reg58=reg58*reg62; reg116=reg35+reg116; reg87=reg87*reg40; reg12=reg12/reg15; reg26=reg26/reg15;
    reg5=reg5*reg62; reg7=reg7/reg15; reg22=reg22/reg15; reg26=reg26*reg116; reg7=reg7*reg116;
    reg4=reg61*reg82; reg62=reg22*reg62; reg6=reg20*reg0; reg58=reg87-reg58; reg10=reg14*reg0;
    reg61=reg102*reg61; reg9=reg5-reg9; reg1=reg2-reg1; reg40=reg12*reg40; reg2=reg16*reg0;
    reg15=reg19/reg15; reg5=reg102*reg64; reg12=reg28*reg0; reg10=reg4-reg10; reg62=reg40-reg62;
    reg116=reg15*reg116; reg26=reg9-reg26; reg6=reg61-reg6; reg4=1-(*f.m).resolution; reg64=reg64*reg82;
    reg0=reg27*reg0; reg58=reg7+reg58; reg20=reg20*reg82; reg2=reg1+reg2; reg14=reg14*reg102;
    reg58=reg4*reg58; reg10=reg10/reg2; reg12=reg64-reg12; reg26=reg4*reg26; reg102=reg28*reg102;
    reg1=(*f.m).resolution*reg111; reg62=reg116+reg62; reg6=reg6/reg2; reg82=reg27*reg82; reg0=reg5-reg0;
    reg5=(*f.m).resolution*reg73; reg14=reg20-reg14; reg65=reg4*reg65; reg52=reg4*reg52; reg26=reg1+reg26;
    reg58=reg5+reg58; reg62=reg4*reg62; reg1=(*f.m).resolution*reg85; reg5=(*f.m).resolution*reg6; reg102=reg82-reg102;
    reg14=reg14/reg2; reg33=reg33/reg2; reg12=reg12/reg2; reg7=(*f.m).resolution*reg10; reg23=reg23/reg2;
    reg0=reg0/reg2; reg91=reg91*reg4; reg9=(*f.m).resolution*reg33; reg5=reg52+reg5; reg13=(*f.m).resolution*reg23;
    reg16=reg16/reg2; reg2=reg102/reg2; reg15=elem.pos(1)[0]-elem.pos(0)[0]; reg7=reg65-reg7; reg17=elem.pos(2)[1]-elem.pos(0)[1];
    reg18=elem.pos(2)[0]-elem.pos(0)[0]; reg8=reg8*reg4; reg56=reg56*reg4; reg126=reg126*reg4; reg114=reg4*reg114;
    reg26=reg26*(*f.m).deltaT; reg19=elem.pos(1)[1]-elem.pos(0)[1]; reg20=(*f.m).resolution*reg12; reg22=(*f.m).resolution*reg14; reg24=(*f.m).resolution*reg0;
    reg58=reg58*(*f.m).deltaT; reg62=reg1+reg62; reg22=reg114+reg22; reg1=(*f.m).resolution*reg2; reg106=reg106*reg4;
    reg25=(*f.m).resolution*reg16; reg27=reg19*reg18; reg95=reg4*reg95; reg62=reg62*(*f.m).deltaT; reg4=reg17*reg15;
    reg20=reg56+reg20; reg28=reg7*reg26; reg24=reg8-reg24; reg13=reg126-reg13; reg91=reg9+reg91;
    reg8=reg5*reg58; reg106=reg25+reg106; reg1=reg95-reg1; reg9=reg20*reg26; reg25=reg24*reg58;
    reg29=reg91*reg58; reg34=reg13*reg26; reg27=reg4-reg27; reg4=reg8+reg28; reg35=reg22*reg62;
    reg17=reg17/reg27; reg36=reg25+reg9; reg37=reg4+reg35; reg15=reg15/reg27; reg19=reg19/reg27;
    reg38=reg106*reg62; reg39=reg29+reg34; reg40=reg1*reg62; reg18=reg18/reg27; reg41=reg36+reg40;
    reg42=0.5*reg18; reg45=0.5*reg15; reg46=reg19-reg17; reg52=reg18-reg15; reg55=0.5*reg19;
    reg56=0.5*reg17; reg57=reg39+reg38; reg59=2*reg37; reg60=0.5*reg52; reg61=reg19*reg57;
    reg63=reg45*reg59; reg64=reg15*reg41; reg65=reg18*reg41; reg66=reg56*reg59; reg67=reg55*reg59;
    reg68=0.5*reg46; reg69=1-var_inter[0]; reg70=reg42*reg59; reg71=reg17*reg57; reg72=reg67-reg64;
    reg74=reg61-reg63; reg75=reg52*reg41; reg76=reg65-reg66; reg77=reg70-reg71; reg78=reg68*reg59;
    reg79=reg59*reg60; reg80=var_inter[0]*elem.f_vol_e[0]; reg81=reg46*reg57; reg82=var_inter[0]*elem.f_vol_e[1]; reg83=var_inter[1]*elem.f_vol_e[0];
    reg69=reg69-var_inter[1]; reg84=var_inter[1]*elem.f_vol_e[1]; reg72=reg72-reg84; reg86=reg69*elem.f_vol_e[1]; reg87=reg75+reg78;
    reg77=reg77-reg80; reg74=reg74-reg83; reg88=reg81+reg79; reg89=reg69*elem.f_vol_e[0]; reg76=reg76-reg82;
    reg72=reg27*reg72; reg74=reg27*reg74; reg76=reg27*reg76; reg90=reg86+reg87; reg92=reg88+reg89;
    reg77=reg27*reg77; reg72=ponderation*reg72; reg93=reg27*reg90; reg74=ponderation*reg74; reg77=ponderation*reg77;
    reg94=reg27*reg92; reg76=ponderation*reg76; sollicitation[indices[2]+0]+=-reg74; sollicitation[indices[1]+0]+=-reg77; reg74=ponderation*reg94;
    sollicitation[indices[0]+0]+=reg74; reg77=ponderation*reg93; sollicitation[indices[0]+1]+=reg77; sollicitation[indices[1]+1]+=-reg76; sollicitation[indices[2]+1]+=-reg72;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; reg1=1.0/reg1; T reg2=2*(*f.m).shear_modulus_12;
    T reg3=reg1*reg0; reg2=1.0/reg2; T reg4=1.0/(*f.m).elastic_modulus_3; T reg5=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg6=reg2*reg3;
    T reg7=pow((*f.m).v1[0],2); T reg8=pow((*f.m).v1[1],2); T reg9=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg10=pow((*f.m).v2[0],2); T reg11=pow((*f.m).v2[1],2);
    T reg12=reg4*reg6; reg8=reg7+reg8; reg7=pow((*f.m).v1[2],2); reg11=reg10+reg11; reg10=pow((*f.m).v2[2],2);
    T reg13=reg9*reg6; T reg14=reg5*reg6; T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; reg7=reg8+reg7;
    reg10=reg11+reg10; reg8=reg5*reg13; reg11=reg16*reg12; T reg17=reg5*reg14; T reg18=reg15*reg12;
    T reg19=reg15*reg13; T reg20=reg16*reg14; reg8=reg11+reg8; T reg21=1.0/(*f.m).elastic_modulus_1; reg10=pow(reg10,0.5);
    reg7=pow(reg7,0.5); reg17=reg18-reg17; reg18=(*f.m).v1[2]/reg7; T reg22=(*f.m).v2[1]/reg10; T reg23=(*f.m).v1[1]/reg7;
    T reg24=(*f.m).v2[2]/reg10; T reg25=reg16*reg8; T reg26=reg20+reg19; T reg27=reg21*reg17; T reg28=reg9*reg13;
    T reg29=reg16*reg6; reg12=reg21*reg12; T reg30=reg5*reg3; reg25=reg27-reg25; reg27=reg4*reg3;
    T reg31=reg9*reg26; reg3=reg9*reg3; T reg32=reg23*reg24; T reg33=reg18*reg22; reg6=reg15*reg6;
    T reg34=reg9*reg14; reg10=(*f.m).v2[0]/reg10; reg7=(*f.m).v1[0]/reg7; T reg35=reg2*reg0; T reg36=reg9*reg6;
    T reg37=reg35*reg9; reg34=reg11+reg34; reg11=reg2*reg1; T reg38=reg35*reg5; reg31=reg25-reg31;
    reg25=reg7*reg24; T reg39=reg18*reg10; T reg40=reg32-reg33; reg13=reg16*reg13; T reg41=reg9*reg29;
    T reg42=reg15*reg27; reg14=reg21*reg14; reg28=reg12-reg28; reg35=reg35*reg4; reg12=2*reg7;
    T reg43=reg3*reg5; T reg44=reg5*reg30; reg27=reg16*reg27; T reg45=2*reg10; reg43=reg27+reg43;
    reg44=reg42-reg44; reg27=reg35*reg15; reg42=reg11*reg5; reg35=reg35*reg16; T reg46=reg11*reg9;
    T reg47=reg38*reg5; T reg48=reg37*reg5; reg11=reg11*reg4; reg3=reg3*reg15; reg30=reg16*reg30;
    reg28=reg28/reg31; T reg49=2*reg40; T reg50=reg22*reg45; T reg51=pow(reg22,2); T reg52=pow(reg10,2);
    reg13=reg14+reg13; T reg53=reg23*reg12; T reg54=pow(reg23,2); T reg55=pow(reg7,2); reg41=reg14+reg41;
    reg14=reg39-reg25; T reg56=reg7*reg22; T reg57=reg23*reg10; reg6=reg21*reg6; reg8=reg8/reg31;
    reg36=reg20+reg36; reg34=reg34/reg31; reg17=reg17/reg31; reg29=reg16*reg29; reg47=reg27-reg47;
    reg27=reg52*reg28; T reg58=reg55*reg34; T reg59=pow(reg24,2); T reg60=pow(reg40,2); T reg61=pow(reg14,2);
    T reg62=reg30+reg3; reg26=reg26/reg31; reg43=reg43*reg16; T reg63=reg50*reg8; T reg64=reg53*reg17;
    T reg65=reg56-reg57; reg44=reg44*reg21; reg41=reg41/reg31; reg13=reg13/reg31; reg36=reg36/reg31;
    T reg66=reg55*reg17; T reg67=reg54*reg17; reg37=reg37*reg15; T reg68=reg51*reg8; reg38=reg38*reg16;
    T reg69=reg46*reg5; T reg70=reg50*reg28; T reg71=reg53*reg34; reg29=reg6-reg29; reg6=reg42*reg5;
    T reg72=reg11*reg16; reg11=reg11*reg15; T reg73=reg14*reg49; T reg74=reg52*reg8; T reg75=reg54*reg34;
    T reg76=pow(reg18,2); reg48=reg35+reg48; reg35=reg51*reg28; T reg77=reg61*reg26; reg70=reg71+reg70;
    reg71=reg73*reg13; T reg78=reg60*reg13; T reg79=reg55*reg36; T reg80=reg52*reg41; T reg81=reg54*reg36;
    T reg82=reg51*reg41; reg27=reg58+reg27; reg58=reg76*reg17; T reg83=reg59*reg8; T reg84=reg53*reg36;
    reg35=reg75+reg35; reg63=reg64+reg63; reg64=reg73*reg26; reg75=reg50*reg41; T reg85=reg61*reg13;
    T reg86=reg59*reg28; T reg87=pow(reg65,2); T reg88=reg76*reg34; reg43=reg44-reg43; reg69=reg72+reg69;
    reg62=reg62*reg9; reg47=reg47*reg21; reg48=reg48*reg16; reg44=reg38+reg37; reg29=reg29/reg31;
    reg46=reg46*reg15; reg72=reg60*reg26; reg66=reg74+reg66; reg68=reg67+reg68; reg6=reg11-reg6;
    reg42=reg42*reg16; reg11=reg87*reg26; reg67=reg42+reg46; reg74=reg10*reg22; reg64=reg63+reg64;
    reg71=reg70+reg71; reg62=reg43-reg62; reg69=reg69*reg16; reg43=2*reg22; reg63=reg24*reg45;
    reg6=reg6*reg21; reg70=2*reg23; T reg89=reg18*reg12; reg78=reg27+reg78; reg48=reg47-reg48;
    reg27=reg87*reg13; reg44=reg44*reg9; reg85=reg35+reg85; reg86=reg88+reg86; reg80=reg79+reg80;
    reg35=reg59*reg41; reg47=reg76*reg36; reg79=reg61*reg29; reg82=reg81+reg82; reg81=reg60*reg29;
    reg77=reg68+reg77; reg68=reg7*reg23; reg75=reg84+reg75; reg73=reg73*reg29; reg72=reg66+reg72;
    reg83=reg58+reg83; reg58=reg64*reg55; reg66=reg54*reg72; reg84=reg78*reg51; reg88=reg18*reg70;
    reg73=reg75+reg73; reg75=reg77*reg68; T reg90=reg85*reg74; T reg91=reg24*reg43; T reg92=reg14*reg40;
    T reg93=reg7*reg14; T reg94=reg23*reg40; reg56=reg57+reg56; reg57=reg23*reg22; T reg95=reg7*reg10;
    T reg96=2*reg14; reg49=reg65*reg49; reg44=reg48-reg44; reg48=reg64*reg68; T reg97=reg87*reg29;
    reg35=reg47+reg35; reg47=reg71*reg74; reg79=reg82+reg79; reg81=reg80+reg81; reg80=reg71*reg51;
    reg64=reg64*reg54; reg82=reg85*reg51; T reg98=reg54*reg77; reg27=reg86+reg27; reg86=reg63*reg28;
    T reg99=reg89*reg34; reg69=reg6-reg69; reg85=reg85*reg52; reg6=reg89*reg17; reg67=reg67*reg9;
    reg77=reg55*reg77; T reg100=reg63*reg8; T reg101=reg55*reg72; T reg102=reg78*reg52; reg62=reg62/reg31;
    reg71=reg71*reg52; reg11=reg83+reg11; reg97=reg35+reg97; reg102=reg101+reg102; reg35=reg81*reg60;
    reg47=reg48+reg47; reg48=reg11*reg54; reg100=reg6+reg100; reg6=reg49*reg26; reg83=reg73*reg92;
    reg80=reg64+reg80; reg64=reg79*reg61; reg82=reg98+reg82; reg98=reg57*reg62; reg72=reg72*reg68;
    reg78=reg78*reg74; reg67=reg69-reg67; reg69=reg22*reg70; reg101=reg79*reg92; reg28=reg91*reg28;
    reg34=reg88*reg34; T reg103=2*reg18; T reg104=reg49*reg13; reg86=reg99+reg86; reg90=reg75+reg90;
    reg75=reg95*reg62; reg99=reg63*reg41; reg8=reg91*reg8; T reg105=reg11*reg55; T reg106=reg27*reg52;
    T reg107=reg10*reg12; reg84=reg66+reg84; reg66=reg89*reg36; reg44=reg44/reg31; reg17=reg88*reg17;
    reg79=reg79*reg60; reg85=reg77+reg85; reg77=reg56*reg62; T reg108=reg81*reg61; reg96=reg65*reg96;
    T reg109=reg73*reg60; reg71=reg58+reg71; reg58=reg18*reg24; reg73=reg73*reg61; T reg110=reg7*reg40;
    T reg111=reg23*reg14; T reg112=reg27*reg51; T reg113=reg22*reg40; T reg114=reg10*reg14; reg93=reg94+reg93;
    reg6=reg100+reg6; reg114=reg113+reg114; reg26=reg96*reg26; reg8=reg17+reg8; reg73=reg80+reg73;
    reg17=reg77*reg69; reg12=reg40*reg12; reg70=reg14*reg70; reg103=reg24*reg103; reg80=reg93*reg44;
    reg111=reg111*reg44; reg110=reg110*reg44; reg94=reg58*reg62; reg101=reg90+reg101; reg41=reg91*reg41;
    reg36=reg88*reg36; reg90=reg98*reg56; reg49=reg49*reg29; reg99=reg66+reg99; reg11=reg11*reg68;
    reg27=reg27*reg74; reg83=reg47+reg83; reg47=reg77*reg56; reg66=reg69*reg15; reg100=reg107*reg16;
    reg113=reg51*reg15; T reg115=reg52*reg16; T reg116=reg69*reg16; T reg117=reg107*reg21; T reg118=reg51*reg16;
    T reg119=reg52*reg21; T reg120=reg97*reg61; reg112=reg48+reg112; reg108=reg84+reg108; reg48=reg98*reg69;
    reg64=reg82+reg64; reg31=reg67/reg31; reg67=reg75*reg69; reg82=reg10*reg40; reg84=reg22*reg14;
    T reg121=reg18*reg65; T reg122=reg97*reg60; reg106=reg105+reg106; reg13=reg96*reg13; reg98=reg98*reg107;
    reg79=reg85+reg79; reg85=reg75*reg107; reg35=reg102+reg35; reg109=reg71+reg109; reg77=reg77*reg107;
    reg78=reg72+reg78; reg104=reg86+reg104; reg81=reg81*reg92; reg28=reg34+reg28; reg34=reg69*reg5;
    reg43=reg14*reg43; reg71=(*f.m).alpha_1*reg55; reg72=reg24*reg65; reg86=reg7*reg65; reg102=reg18*reg40;
    reg105=reg51*(*f.m).alpha_2; reg47=reg83+reg47; reg83=reg80*reg93; reg39=reg25+reg39; reg13=reg28+reg13;
    reg81=reg78+reg81; reg75=reg75*reg56; reg100=reg66-reg100; reg90=reg101+reg90; reg25=reg111*reg93;
    reg28=reg104*reg51; reg66=reg6*reg54; reg78=reg80*reg70; reg17=reg73+reg17; reg26=reg8+reg26;
    reg27=reg11+reg27; reg97=reg97*reg92; reg8=reg52*reg9; reg11=reg51*reg5; reg45=reg40*reg45;
    reg73=reg107*reg9; reg49=reg99+reg49; reg41=reg36+reg41; reg29=reg96*reg29; reg36=(*f.m).alpha_1*reg54;
    reg96=reg103*reg9; reg120=reg112+reg120; reg99=reg94*reg69; reg116=reg117-reg116; reg77=reg109+reg77;
    reg121=reg121*reg44; reg82=reg31*reg82; reg84=reg31*reg84; reg101=reg31*reg114; reg122=reg106+reg122;
    reg85=reg35+reg85; reg35=reg110*reg12; reg106=reg111*reg12; reg98=reg79+reg98; reg79=reg94*reg107;
    reg109=reg59*reg9; reg112=reg110*reg70; reg67=reg108+reg67; reg108=reg103*reg5; reg111=reg111*reg70;
    reg115=reg113-reg115; reg118=reg119-reg118; reg80=reg80*reg12; reg113=reg6*reg55; reg117=reg59*reg5;
    reg119=reg104*reg52; T reg123=reg52*(*f.m).alpha_2; reg48=reg64+reg48; reg64=reg23*reg65; T reg124=reg59*reg4;
    reg122=reg79+reg122; reg79=reg49*reg61; reg28=reg66+reg28; reg66=reg10*reg65; T reg125=reg24*reg40;
    T reg126=reg54*reg16; reg11=reg8+reg11; reg78=reg17+reg78; reg8=reg101*reg43; reg17=reg82*reg43;
    reg112=reg67+reg112; reg111=reg48+reg111; reg104=reg104*reg74; reg6=reg6*reg68; reg48=reg101*reg114;
    reg83=reg47+reg83; reg47=(*f.m).alpha_3*reg61; reg67=reg13*reg52; T reg127=reg26*reg55; T reg128=reg84*reg43;
    reg105=reg36+reg105; reg99=reg120+reg99; reg36=reg121*reg70; reg120=reg49*reg60; reg108=reg100-reg108;
    reg119=reg113+reg119; reg32=reg33+reg32; reg86=reg102+reg86; reg33=reg18*reg14; reg100=reg121*reg12;
    reg102=reg13*reg51; reg113=reg26*reg54; reg106=reg98+reg106; reg98=reg82*reg45; reg35=reg85+reg35;
    reg15=reg54*reg15; reg21=reg55*reg21; reg72=reg31*reg72; reg75=reg81+reg75; reg110=reg110*reg93;
    reg81=reg39*reg62; reg29=reg41+reg29; reg25=reg90+reg25; reg41=reg84*reg114; reg97=reg27+reg97;
    reg94=reg94*reg56; reg34=reg73+reg34; reg103=reg103*reg4; reg27=(*f.m).alpha_3*reg60; reg123=reg71+reg123;
    reg117=reg115-reg117; reg84=reg84*reg45; reg71=reg59*(*f.m).alpha_2; reg73=(*f.m).alpha_1*reg76; reg101=reg101*reg45;
    reg80=reg77+reg80; reg109=reg118-reg109; reg96=reg116-reg96; reg16=reg55*reg16; reg62=reg32*reg62;
    reg128=reg111+reg128; reg71=reg73+reg71; reg27=reg123+reg27; reg87=(*f.m).alpha_3*reg87; reg41=reg25+reg41;
    reg25=reg76*reg9; reg48=reg83+reg48; reg126=reg21-reg126; reg94=reg97+reg94; reg104=reg6+reg104;
    reg49=reg49*reg92; reg121=reg121*reg93; reg64=reg33+reg64; reg34=reg103-reg34; reg6=reg117*reg54;
    reg21=reg109*reg55; reg33=reg24*reg14; reg73=reg22*reg65; reg77=reg117*reg51; reg83=reg109*reg52;
    reg11=reg124-reg11; reg98=reg35+reg98; reg84=reg106+reg84; reg35=reg76*reg5; reg66=reg125+reg66;
    reg16=reg15-reg16; reg15=(*f.m).alpha_1*reg68; reg85=reg74*(*f.m).alpha_2; reg5=reg54*reg5; reg9=reg55*reg9;
    reg86=reg86*reg44; reg47=reg105+reg47; reg7=reg7*reg18; reg36=reg99+reg36; reg10=reg10*reg24;
    reg110=reg75+reg110; reg82=reg82*reg114; reg26=reg26*reg68; reg60=reg29*reg60; reg120=reg119+reg120;
    reg75=reg72*reg43; reg101=reg80+reg101; reg80=reg96*reg95; reg90=reg81*reg107; reg67=reg127+reg67;
    reg97=reg51*reg108; reg99=reg96*reg52; reg100=reg122+reg100; reg102=reg113+reg102; reg103=reg57*reg108;
    reg105=reg72*reg45; reg106=reg81*reg69; reg96=reg96*reg55; reg108=reg54*reg108; reg61=reg29*reg61;
    reg79=reg28+reg79; reg117=reg117*reg57; reg8=reg78+reg8; reg109=reg109*reg95; reg13=reg13*reg74;
    reg17=reg112+reg17; reg106=reg79+reg106; reg28=reg86*reg70; reg29=reg29*reg92; reg73=reg33+reg73;
    reg33=reg59*reg11; reg77=reg83+reg77; reg78=reg58*reg34; reg87=reg71+reg87; reg5=reg9+reg5;
    reg4=reg76*reg4; reg61=reg102+reg61; reg69=reg62*reg69; reg85=reg15+reg85; reg92=(*f.m).alpha_3*reg92;
    reg64=reg44*reg64; reg66=reg31*reg66; reg82=reg110+reg82; reg9=reg76*reg34; reg108=reg96+reg108;
    reg117=reg109+reg117; reg15=reg58*reg11; reg44=reg128*reg47; reg71=reg17*reg27; reg105=reg100+reg105;
    reg97=reg99+reg97; reg34=reg59*reg34; reg35=reg16-reg35; reg6=reg21+reg6; reg11=reg76*reg11;
    reg75=reg36+reg75; reg25=reg126-reg25; reg16=reg84*reg47; reg103=reg80+reg103; reg21=reg98*reg27;
    reg36=reg101*reg41; reg79=reg8*reg41; reg80=reg84*reg48; reg83=reg128*reg48; reg81=reg81*reg56;
    reg49=reg104+reg49; reg74=reg74*reg2; reg107=reg62*reg107; reg60=reg67+reg60; reg67=reg56*reg2;
    reg96=reg86*reg12; reg24=reg22*reg24; reg18=reg23*reg18; reg40=reg65*reg40; reg90=reg120+reg90;
    reg13=reg26+reg13; reg72=reg72*reg114; reg121=reg94+reg121; reg22=(*f.m).alpha_1*reg7; reg23=reg10*(*f.m).alpha_2;
    reg26=reg50*reg67; reg34=reg97+reg34; reg70=reg64*reg70; reg69=reg61+reg69; reg81=reg49+reg81;
    reg23=reg22+reg23; reg86=reg86*reg93; reg44=reg71+reg44; reg75=reg87*reg75; reg22=reg56*reg74;
    reg15=reg117+reg15; reg92=reg85+reg92; reg49=reg47*reg41; reg79=reg83-reg79; reg61=reg82*reg27;
    reg9=reg108+reg9; reg36=reg80-reg36; reg71=reg84*reg8; reg80=reg35*reg54; reg83=reg25*reg55;
    reg73=reg31*reg73; reg31=reg101*reg128; reg62=reg62*reg56; reg85=reg53*reg67; reg94=reg25*reg52;
    reg72=reg121+reg72; reg33=reg77+reg33; reg77=reg50*reg74; reg5=reg4-reg5; reg28=reg106+reg28;
    reg4=reg66*reg43; reg14=reg65*reg14; reg67=reg56*reg67; reg96=reg90+reg96; reg29=reg13+reg29;
    reg10=reg10*reg1; reg13=reg66*reg45; reg65=reg39*reg1; reg16=reg21+reg16; reg105=reg105*reg87;
    reg40=(*f.m).alpha_3*reg40; reg74=reg53*reg74; reg12=reg64*reg12; reg21=reg24*(*f.m).alpha_2; reg90=reg35*reg51;
    reg107=reg60+reg107; reg6=reg11+reg6; reg78=reg103+reg78; reg11=(*f.m).alpha_1*reg18; reg85=reg9+reg85;
    reg9=reg89*reg65; reg60=reg32*reg0; reg97=reg39*reg65; reg99=reg39*reg10; reg2=reg68*reg2;
    reg49=reg61+reg49; reg67=reg78+reg67; reg22=reg15+reg22; reg24=reg24*reg0; reg72=reg87*reg72;
    reg105=reg16+reg105; reg15=reg63*reg10; reg77=reg33+reg77; reg4=reg28+reg4; reg40=reg23+reg40;
    reg21=reg11+reg21; reg14=(*f.m).alpha_3*reg14; reg70=reg69+reg70; reg43=reg73*reg43; reg59=reg59*reg5;
    reg13=reg96+reg13; reg90=reg94+reg90; reg12=reg107+reg12; reg86=reg81+reg86; reg66=reg66*reg114;
    reg11=reg98*reg79; reg16=reg17*reg36; reg45=reg73*reg45; reg31=reg71-reg31; reg93=reg64*reg93;
    reg62=reg29+reg62; reg35=reg35*reg57; reg80=reg83+reg80; reg76=reg76*reg5; reg25=reg25*reg95;
    reg74=reg6+reg74; reg10=reg89*reg10; reg6=reg92*reg8; reg75=reg44+reg75; reg26=reg34+reg26;
    reg65=reg63*reg65; reg23=reg101*reg92; reg10=reg74+reg10; reg99=reg22+reg99; reg22=reg32*reg24;
    reg28=reg32*reg60; reg5=reg58*reg5; reg35=reg25+reg35; reg25=reg88*reg24; reg45=reg12+reg45;
    reg97=reg67+reg97; reg12=reg88*reg60; reg9=reg85+reg9; reg4=reg40*reg4; reg6=reg75+reg6;
    reg14=reg21+reg14; reg65=reg26+reg65; reg60=reg91*reg60; reg21=reg50*reg2; reg59=reg90+reg59;
    reg13=reg13*reg40; reg23=reg105+reg23; reg43=reg70+reg43; reg16=reg11-reg16; reg11=reg31*reg82;
    reg26=reg17*reg48; reg29=reg8*reg82; reg66=reg86+reg66; reg33=reg98*reg48; reg34=reg101*reg82;
    reg1=reg7*reg1; reg93=reg62+reg93; reg72=reg49+reg72; reg7=reg92*reg48; reg114=reg73*reg114;
    reg15=reg77+reg15; reg24=reg91*reg24; reg76=reg80+reg76; reg44=reg53*reg2; reg49=reg101*reg17;
    reg66=reg40*reg66; reg7=reg72+reg7; reg40=reg98*reg8; reg28=reg97+reg28; reg25=reg10+reg25;
    reg10=reg84*reg82; reg34=reg33-reg34; reg33=reg98*reg41; reg22=reg99+reg22; reg60=reg65+reg60;
    reg114=reg93+reg114; reg13=reg23+reg13; reg45=reg45*reg14; reg4=reg6+reg4; reg43=reg14*reg43;
    reg2=reg56*reg2; reg5=reg35+reg5; reg11=reg16+reg11; reg24=reg15+reg24; reg12=reg9+reg12;
    reg6=reg17*reg41; reg29=reg26-reg29; reg63=reg63*reg1; reg89=reg89*reg1; reg9=reg128*reg82;
    reg44=reg76+reg44; reg0=reg18*reg0; reg21=reg59+reg21; reg45=reg13+reg45; reg63=reg21+reg63;
    reg91=reg91*reg0; reg43=reg4+reg43; reg2=reg5+reg2; reg1=reg39*reg1; reg4=reg60*reg22;
    reg5=elem.pos(1)[1]-elem.pos(0)[1]; reg88=reg88*reg0; reg13=elem.pos(2)[1]-elem.pos(0)[1]; reg15=elem.pos(2)[0]-elem.pos(0)[0]; reg89=reg44+reg89;
    reg114=reg14*reg114; reg66=reg7+reg66; reg7=reg12*reg22; reg79=reg79/reg11; reg10=reg33-reg10;
    reg14=reg128*reg98; reg49=reg40-reg49; reg29=reg29/reg11; reg16=reg17*reg84; reg18=elem.pos(1)[0]-elem.pos(0)[0];
    reg34=reg34/reg11; reg21=reg24*reg28; reg23=reg25*reg28; reg36=reg36/reg11; reg9=reg6-reg9;
    reg6=reg13*reg18; reg10=reg10/reg11; reg114=reg66+reg114; reg4=reg21-reg4; reg9=reg9/reg11;
    reg21=reg5*reg15; reg31=reg31/reg11; reg88=reg89+reg88; reg26=reg25*reg60; reg7=reg23-reg7;
    reg23=reg12*reg24; reg16=reg14-reg16; reg91=reg63+reg91; reg79=reg79*reg45; reg36=reg36*reg43;
    reg29=reg29*reg45; reg1=reg2+reg1; reg0=reg32*reg0; reg34=reg34*reg43; reg49=reg49/reg11;
    reg2=reg88*reg4; reg21=reg6-reg21; reg6=reg91*reg7; reg11=reg16/reg11; reg23=reg26-reg23;
    reg43=reg10*reg43; reg31=reg31*reg114; reg36=reg79-reg36; reg49=reg49*reg114; reg29=reg34-reg29;
    reg45=reg9*reg45; reg0=reg1+reg0; reg36=reg31+reg36; reg114=reg11*reg114; reg1=1-(*f.m).resolution;
    reg5=reg5/reg21; reg49=reg29-reg49; reg18=reg18/reg21; reg13=reg13/reg21; reg6=reg2-reg6;
    reg2=reg23*reg0; reg9=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg10=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg11=reg91*reg22; reg15=reg15/reg21;
    reg14=reg24*reg0; reg16=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; reg26=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; reg22=reg22*reg88; reg43=reg45-reg43;
    reg29=reg25*reg0; reg43=reg114+reg43; reg49=reg1*reg49; reg25=reg25*reg91; reg24=reg24*reg88;
    reg2=reg6+reg2; reg6=reg91*reg28; reg31=reg60*reg0; reg28=reg28*reg88; reg14=reg11-reg14;
    reg0=reg12*reg0; reg29=reg22-reg29; reg11=reg9*reg13; reg22=reg10*reg5; reg32=(*f.m).resolution*reg27;
    reg33=reg15*reg16; reg34=reg26*reg18; reg35=(*f.m).resolution*reg47; reg36=reg1*reg36; reg25=reg24-reg25;
    reg33=reg34-reg33; reg91=reg12*reg91; reg9=reg9*reg15; reg10=reg10*reg18; reg31=reg6-reg31;
    reg14=reg14/reg2; reg26=reg26*reg5; reg16=reg13*reg16; reg0=reg28-reg0; reg29=reg29/reg2;
    reg88=reg60*reg88; reg49=reg35+reg49; reg36=reg32+reg36; reg22=reg11-reg22; reg43=reg1*reg43;
    reg6=(*f.m).resolution*reg92; reg36=reg36*(*f.m).deltaT; reg43=reg6+reg43; reg0=reg0/reg2; reg49=reg49*(*f.m).deltaT;
    reg7=reg7/reg2; reg26=reg16-reg26; reg82=reg1*reg82; reg41=reg1*reg41; reg6=(*f.m).resolution*reg29;
    reg31=reg31/reg2; reg11=(*f.m).resolution*reg14; reg4=reg4/reg2; reg9=reg10-reg9; reg25=reg25/reg2;
    reg91=reg88-reg91; reg22=reg33+reg22; reg10=(*f.m).resolution*reg7; reg91=reg91/reg2; reg17=reg17*reg1;
    reg12=(*f.m).resolution*reg4; reg128=reg128*reg1; reg16=(*f.m).resolution*reg0; reg24=(*f.m).resolution*reg31; reg98=reg98*reg1;
    reg84=reg84*reg1; reg48=reg1*reg48; reg28=(*f.m).resolution*reg25; reg26=reg26-reg36; reg9=reg9-reg49;
    reg22=0.5*reg22; reg43=reg43*(*f.m).deltaT; reg6=reg41-reg6; reg11=reg82+reg11; reg2=reg23/reg2;
    reg16=reg128+reg16; reg24=reg17-reg24; reg17=(*f.m).resolution*reg2; reg28=reg48+reg28; reg98=reg12+reg98;
    reg12=reg26*reg11; reg23=reg9*reg6; reg101=reg101*reg1; reg8=reg1*reg8; reg22=reg22-reg43;
    reg10=reg84-reg10; reg1=(*f.m).resolution*reg91; reg32=reg26*reg98; reg33=reg9*reg10; reg26=reg26*reg24;
    reg9=reg9*reg16; reg34=reg28*reg22; reg23=reg12+reg23; reg101=reg17+reg101; reg1=reg8-reg1;
    reg8=reg5-reg13; reg12=reg101*reg22; reg17=reg15-reg18; reg33=reg32+reg33; reg9=reg26+reg9;
    reg34=reg23+reg34; reg22=reg1*reg22; reg23=0.5*reg17; reg26=0.5*reg8; reg32=1-var_inter[0];
    reg35=0.5*reg5; reg39=0.5*reg18; reg40=0.5*reg15; reg41=0.5*reg13; reg22=reg9+reg22;
    reg34=2*reg34; reg12=reg33+reg12; reg9=reg8*reg12; reg33=reg39*reg34; reg44=reg5*reg12;
    reg45=reg34*reg23; reg48=reg35*reg34; reg58=reg18*reg22; reg59=reg13*reg12; reg60=reg40*reg34;
    reg61=reg26*reg34; reg62=reg15*reg22; reg63=reg41*reg34; reg64=reg17*reg22; reg32=reg32-var_inter[1];
    reg65=var_inter[1]*elem.f_vol_e[0]; reg66=reg32*elem.f_vol_e[1]; reg61=reg64+reg61; reg33=reg33-reg44; reg58=reg58-reg48;
    reg64=var_inter[1]*elem.f_vol_e[1]; reg67=var_inter[0]*elem.f_vol_e[1]; reg45=reg9+reg45; reg9=reg32*elem.f_vol_e[0]; reg63=reg63-reg62;
    reg59=reg59-reg60; reg68=var_inter[0]*elem.f_vol_e[0]; reg63=reg63-reg67; reg33=reg33-reg65; reg45=reg45-reg9;
    reg59=reg59-reg68; reg61=reg61-reg66; reg58=reg58-reg64; reg61=reg61*reg21; reg59=reg21*reg59;
    reg63=reg21*reg63; reg45=reg45*reg21; reg33=reg21*reg33; reg58=reg21*reg58; sollicitation[indices[0]+1]+=ponderation*reg61;
    sollicitation[indices[1]+0]+=ponderation*reg59; sollicitation[indices[0]+0]+=ponderation*reg45; sollicitation[indices[1]+1]+=ponderation*reg63; sollicitation[indices[2]+1]+=ponderation*reg58; sollicitation[indices[2]+0]+=ponderation*reg33;
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
  static const unsigned nb_matrices = 1;
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
    node.dep[0]=vecs[0][indice+0]; node.dep[1]=vecs[0][indice+1];
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_nodal_initial_conditions(const TE &node,const TTs &f,Tvec &vecs,unsigned indice) {
    vecs[0][indice+0]=node.dep[0]; vecs[1][indice+0]=node.dep[0]; vecs[2][indice+0]=node.dep[0]; vecs[3][indice+0]=node.dep[0]; vecs[4][indice+0]=node.dep[0];
    vecs[0][indice+1]=node.dep[1]; vecs[1][indice+1]=node.dep[1]; vecs[2][indice+1]=node.dep[1]; vecs[3][indice+1]=node.dep[1]; vecs[4][indice+1]=node.dep[1];
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_nodal_error(const TE &node,const TTs &f,const Tvec &vecs,int indice) {
    T reg0=vecs[1][indice+0]-vecs[0][indice+0]; T reg1=vecs[1][indice+1]-vecs[0][indice+1]; reg1=abs(reg1); reg0=abs(reg0); return max(reg1,reg0);
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_nodal(const TE &node,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
    old_vec[indice+0]=vecs[1][indice+0]; old_vec[indice+1]=vecs[1][indice+1];
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; reg1=1.0/reg1; T reg2=2*(*f.m).shear_modulus_12;
    T reg3=reg1*reg0; reg2=1.0/reg2; T reg4=reg2*reg3; T reg5=pow((*f.m).v2[1],2); T reg6=pow((*f.m).v2[0],2);
    T reg7=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg8=1.0/(*f.m).elastic_modulus_3; T reg9=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg10=pow((*f.m).v1[1],2); T reg11=pow((*f.m).v1[0],2);
    T reg12=reg7*reg4; T reg13=reg9*reg4; T reg14=reg8*reg4; reg10=reg11+reg10; reg11=pow((*f.m).v1[2],2);
    T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; reg5=reg6+reg5; reg6=pow((*f.m).v2[2],2); T reg17=reg7*reg12;
    T reg18=reg16*reg14; T reg19=reg7*reg13; T reg20=reg15*reg14; reg11=reg10+reg11; reg6=reg5+reg6;
    reg17=reg20-reg17; reg11=pow(reg11,0.5); reg19=reg18+reg19; reg5=reg16*reg12; reg10=reg15*reg13;
    reg20=1.0/(*f.m).elastic_modulus_1; reg6=pow(reg6,0.5); T reg21=(*f.m).v1[2]/reg11; T reg22=(*f.m).v1[1]/reg11; T reg23=reg20*reg17;
    T reg24=reg16*reg19; T reg25=reg5+reg10; T reg26=(*f.m).v2[1]/reg6; T reg27=(*f.m).v2[2]/reg6; T reg28=reg7*reg3;
    T reg29=reg2*reg0; reg6=(*f.m).v2[0]/reg6; T reg30=reg16*reg4; T reg31=reg9*reg13; reg14=reg20*reg14;
    T reg32=reg9*reg12; reg4=reg15*reg4; T reg33=reg9*reg25; T reg34=reg21*reg26; T reg35=reg22*reg27;
    reg24=reg23-reg24; reg23=reg9*reg3; reg11=(*f.m).v1[0]/reg11; reg3=reg8*reg3; T reg36=reg9*reg30;
    T reg37=reg35-reg34; T reg38=reg21*reg6; T reg39=reg11*reg27; T reg40=reg7*reg29; reg13=reg16*reg13;
    T reg41=reg9*reg29; T reg42=2*reg11; reg12=reg20*reg12; T reg43=reg9*reg4; T reg44=reg15*reg3;
    reg29=reg8*reg29; reg31=reg14-reg31; reg32=reg18+reg32; reg14=2*reg6; reg18=reg7*reg23;
    T reg45=reg2*reg1; reg33=reg24-reg33; reg24=reg7*reg28; reg3=reg16*reg3; T reg46=reg14*reg26;
    T reg47=2*reg37; reg31=reg31/reg33; T reg48=reg8*reg45; T reg49=reg11*reg26; T reg50=reg22*reg6;
    T reg51=reg9*reg45; reg45=reg7*reg45; reg19=reg19/reg33; reg43=reg5+reg43; reg32=reg32/reg33;
    reg4=reg20*reg4; reg24=reg44-reg24; reg23=reg15*reg23; reg30=reg16*reg30; reg44=reg7*reg41;
    reg17=reg17/reg33; reg18=reg3+reg18; reg3=pow(reg26,2); T reg52=pow(reg6,2); T reg53=reg15*reg29;
    reg29=reg16*reg29; T reg54=reg7*reg40; T reg55=reg22*reg42; reg13=reg12+reg13; T reg56=pow(reg22,2);
    reg28=reg16*reg28; reg36=reg12+reg36; reg12=pow(reg11,2); T reg57=reg38-reg39; T reg58=reg3*reg31;
    T reg59=reg7*reg45; reg44=reg29+reg44; reg29=reg46*reg31; T reg60=reg16*reg48; T reg61=reg55*reg32;
    reg48=reg15*reg48; T reg62=reg56*reg32; reg40=reg16*reg40; T reg63=reg52*reg31; T reg64=reg12*reg32;
    reg25=reg25/reg33; reg41=reg15*reg41; reg43=reg43/reg33; T reg65=reg49-reg50; T reg66=pow(reg21,2);
    T reg67=reg7*reg51; T reg68=reg56*reg17; T reg69=reg3*reg19; T reg70=reg55*reg17; T reg71=reg46*reg19;
    T reg72=reg12*reg17; T reg73=reg52*reg19; reg24=reg20*reg24; reg18=reg16*reg18; reg36=reg36/reg33;
    reg54=reg53-reg54; reg30=reg4-reg30; reg4=pow(reg27,2); reg13=reg13/reg33; reg53=pow(reg37,2);
    T reg74=reg28+reg23; T reg75=pow(reg57,2); T reg76=reg47*reg57; reg73=reg72+reg73; reg69=reg68+reg69;
    reg68=reg4*reg19; reg54=reg20*reg54; reg72=reg75*reg25; T reg77=reg53*reg25; reg59=reg48-reg59;
    reg48=reg4*reg31; T reg78=reg66*reg32; T reg79=pow(reg65,2); T reg80=reg75*reg13; reg63=reg64+reg63;
    reg58=reg62+reg58; reg62=reg66*reg17; reg51=reg15*reg51; reg64=reg53*reg13; reg74=reg9*reg74;
    T reg81=reg76*reg25; T reg82=reg3*reg36; T reg83=reg56*reg43; reg30=reg30/reg33; T reg84=reg52*reg36;
    T reg85=reg12*reg43; reg71=reg70+reg71; reg67=reg60+reg67; reg44=reg16*reg44; reg60=reg55*reg43;
    reg70=reg40+reg41; reg18=reg24-reg18; reg24=reg76*reg13; reg29=reg61+reg29; reg61=reg46*reg36;
    reg45=reg16*reg45; reg59=reg20*reg59; reg72=reg69+reg72; reg69=reg45+reg51; reg77=reg73+reg77;
    reg67=reg16*reg67; reg73=reg21*reg42; T reg86=2*reg22; reg74=reg18-reg74; reg44=reg54-reg44;
    reg70=reg9*reg70; reg84=reg85+reg84; reg18=reg53*reg30; reg54=reg6*reg26; reg24=reg29+reg24;
    reg61=reg60+reg61; reg29=reg11*reg22; reg82=reg83+reg82; reg60=reg75*reg30; reg83=reg76*reg30;
    reg85=reg66*reg43; T reg87=reg4*reg36; T reg88=reg79*reg13; reg48=reg78+reg48; reg80=reg58+reg80;
    reg58=reg14*reg27; reg78=2*reg26; T reg89=reg79*reg25; reg81=reg71+reg81; reg64=reg63+reg64;
    reg68=reg62+reg68; reg62=reg58*reg19; reg63=reg73*reg17; reg67=reg59-reg67; reg69=reg9*reg69;
    reg59=reg57*reg37; reg74=reg74/reg33; reg71=reg29*reg72; T reg90=reg54*reg80; T reg91=reg52*reg24;
    reg89=reg68+reg89; reg68=reg56*reg72; T reg92=reg12*reg81; T reg93=reg3*reg80; T reg94=reg22*reg37;
    T reg95=reg3*reg24; T reg96=reg56*reg81; reg49=reg50+reg49; reg50=reg22*reg26; T reg97=reg11*reg6;
    T reg98=2*reg57; reg47=reg47*reg65; T reg99=reg78*reg27; reg83=reg61+reg83; reg70=reg44-reg70;
    reg88=reg48+reg88; reg44=reg73*reg32; reg48=reg58*reg31; reg61=reg79*reg30; reg18=reg84+reg18;
    reg60=reg82+reg60; reg87=reg85+reg87; reg82=reg3*reg64; reg84=reg21*reg86; reg85=reg52*reg80;
    T reg100=reg12*reg72; T reg101=reg54*reg24; T reg102=reg29*reg81; T reg103=reg56*reg77; T reg104=reg11*reg57;
    T reg105=reg52*reg64; T reg106=reg12*reg77; reg101=reg102+reg101; reg102=reg22*reg57; reg70=reg70/reg33;
    reg104=reg94+reg104; reg94=reg6*reg57; T reg107=reg11*reg37; reg69=reg67-reg69; reg82=reg103+reg82;
    reg67=reg59*reg83; reg103=reg75*reg18; T reg108=reg75*reg60; reg93=reg68+reg93; reg98=reg98*reg65;
    reg68=reg26*reg37; T reg109=reg21*reg27; T reg110=reg75*reg83; T reg111=reg53*reg60; reg85=reg100+reg85;
    reg100=reg50*reg74; T reg112=reg6*reg42; T reg113=reg49*reg74; reg32=reg84*reg32; reg31=reg99*reg31;
    T reg114=reg54*reg64; T reg115=reg29*reg77; T reg116=reg53*reg18; reg105=reg106+reg105; reg61=reg87+reg61;
    reg87=reg59*reg60; reg106=reg73*reg43; T reg117=reg58*reg36; reg90=reg71+reg90; reg71=reg56*reg89;
    T reg118=reg97*reg74; T reg119=reg3*reg88; reg48=reg44+reg48; reg19=reg99*reg19; reg44=reg53*reg83;
    reg91=reg92+reg91; reg17=reg84*reg17; reg92=reg47*reg13; T reg120=reg47*reg25; T reg121=reg52*reg88;
    T reg122=2*reg21; T reg123=reg12*reg89; reg62=reg63+reg62; reg63=reg26*reg86; reg95=reg96+reg95;
    reg96=reg47*reg30; reg43=reg84*reg43; reg117=reg106+reg117; reg13=reg98*reg13; reg31=reg32+reg31;
    reg92=reg48+reg92; reg32=reg3*reg16; reg48=reg52*reg20; reg106=reg49*reg113; reg67=reg101+reg67;
    reg101=reg112*reg113; reg44=reg91+reg44; reg91=reg63*reg15; T reg124=reg53*reg61; reg121=reg123+reg121;
    reg123=reg112*reg16; T reg125=reg112*reg100; reg111=reg85+reg111; reg85=reg3*reg15; T reg126=reg52*reg16;
    reg94=reg68+reg94; reg68=reg112*reg118; reg116=reg105+reg116; reg105=reg26*reg57; T reg127=reg6*reg37;
    T reg128=reg21*reg65; T reg129=reg63*reg16; T reg130=reg104*reg70; T reg131=reg102*reg70; T reg132=reg107*reg70;
    T reg133=reg112*reg20; T reg134=reg109*reg74; reg36=reg99*reg36; T reg135=reg29*reg89; T reg136=reg54*reg88;
    reg42=reg37*reg42; reg103=reg82+reg103; reg82=reg63*reg118; T reg137=reg75*reg61; T reg138=reg59*reg18;
    reg114=reg115+reg114; reg86=reg57*reg86; reg115=reg27*reg122; reg119=reg71+reg119; reg71=reg49*reg100;
    T reg139=reg63*reg100; reg25=reg98*reg25; reg19=reg17+reg19; reg120=reg62+reg120; reg33=reg69/reg33;
    reg110=reg95+reg110; reg17=reg63*reg113; reg108=reg93+reg108; reg87=reg90+reg87; reg71=reg87+reg71;
    reg62=reg104*reg131; reg69=reg27*reg65; reg68=reg116+reg68; reg87=reg63*reg7; reg90=reg42*reg132;
    reg139=reg108+reg139; reg93=reg112*reg9; reg95=reg86*reg131; reg108=reg4*reg7; reg136=reg135+reg136;
    reg116=reg59*reg61; reg14=reg14*reg37; reg125=reg111+reg125; reg111=reg42*reg131; reg126=reg85-reg126;
    reg85=reg52*reg92; reg78=reg78*reg57; reg124=reg121+reg124; reg121=reg112*reg134; reg135=reg115*reg7;
    reg123=reg91-reg123; reg91=reg12*reg120; T reg140=reg52*reg9; T reg141=reg3*reg7; reg30=reg98*reg30;
    reg36=reg43+reg36; reg138=reg114+reg138; reg96=reg117+reg96; reg43=(*f.m).alpha_1*reg56; reg106=reg67+reg106;
    reg67=reg104*reg130; reg114=(*f.m).alpha_2*reg52; reg13=reg31+reg13; reg31=reg3*reg92; reg117=reg56*reg120;
    T reg142=reg86*reg130; reg17=reg110+reg17; reg32=reg48-reg32; reg48=reg4*reg9; reg25=reg19+reg25;
    reg19=reg11*reg65; reg110=reg21*reg37; reg137=reg119+reg137; reg119=reg94*reg33; reg129=reg133-reg129;
    reg133=reg86*reg132; T reg143=reg105*reg33; T reg144=reg127*reg33; T reg145=reg115*reg9; reg38=reg39+reg38;
    reg82=reg103+reg82; reg39=reg128*reg70; reg103=reg63*reg134; T reg146=(*f.m).alpha_1*reg12; reg101=reg44+reg101;
    reg44=reg42*reg130; T reg147=(*f.m).alpha_2*reg3; T reg148=reg49*reg118; reg141=reg140+reg141; reg140=reg115*reg8;
    reg87=reg93+reg87; reg103=reg137+reg103; reg93=reg86*reg39; reg95=reg139+reg95; reg137=reg78*reg143;
    reg67=reg106+reg67; reg106=reg4*reg8; reg135=reg123-reg135; reg108=reg126-reg108; reg123=reg78*reg144;
    reg133=reg82+reg133; reg82=reg56*reg15; reg126=reg12*reg16; reg145=reg129-reg145; reg129=(*f.m).alpha_2*reg4;
    reg139=(*f.m).alpha_1*reg66; T reg149=(*f.m).alpha_3*reg75; reg147=reg43+reg147; reg48=reg32-reg48; reg32=(*f.m).alpha_3*reg53;
    reg114=reg146+reg114; reg43=reg56*reg16; reg146=reg6*reg65; T reg150=reg27*reg37; T reg151=reg22*reg65;
    T reg152=reg21*reg57; reg19=reg110+reg19; reg35=reg34+reg35; reg34=reg42*reg39; reg121=reg124+reg121;
    reg110=reg14*reg143; reg111=reg125+reg111; reg124=reg14*reg144; reg90=reg68+reg90; reg68=reg12*reg20;
    reg125=reg69*reg33; T reg153=reg38*reg74; reg30=reg36+reg30; reg142=reg17+reg142; reg17=reg78*reg119;
    reg31=reg117+reg31; reg36=reg75*reg96; reg117=reg56*reg25; T reg154=reg3*reg13; T reg155=reg49*reg134;
    reg116=reg136+reg116; reg136=reg94*reg143; reg148=reg138+reg148; reg138=reg104*reg132; reg62=reg71+reg62;
    reg71=reg12*reg25; T reg156=reg94*reg119; T reg157=reg52*reg13; reg44=reg101+reg44; reg101=reg14*reg119;
    T reg158=reg53*reg96; T reg159=reg29*reg120; T reg160=reg54*reg92; reg85=reg91+reg85; reg11=reg11*reg21;
    reg6=reg6*reg27; reg91=reg97*reg145; T reg161=reg50*reg135; T reg162=reg75*reg30; reg137=reg95+reg137;
    reg154=reg117+reg154; reg95=reg78*reg125; reg117=reg63*reg153; T reg163=reg48*reg12; T reg164=reg108*reg56;
    reg36=reg31+reg36; reg93=reg103+reg93; reg17=reg142+reg17; reg31=reg145*reg12; reg103=(*f.m).alpha_2*reg54;
    reg142=(*f.m).alpha_1*reg29; T reg165=(*f.m).alpha_3*reg79; reg129=reg139+reg129; reg139=reg94*reg144; reg149=reg147+reg149;
    reg138=reg148+reg138; reg156=reg67+reg156; reg67=reg97*reg48; reg147=reg50*reg108; reg136=reg62+reg136;
    reg160=reg159+reg160; reg62=reg59*reg96; reg155=reg116+reg155; reg116=reg104*reg39; reg32=reg114+reg32;
    reg157=reg71+reg157; reg71=reg53*reg30; reg146=reg150+reg146; reg114=reg26*reg65; reg148=reg108*reg3;
    reg126=reg82-reg126; reg151=reg152+reg151; reg82=reg112*reg153; reg124=reg90+reg124; reg141=reg106-reg141;
    reg101=reg44+reg101; reg158=reg85+reg158; reg110=reg111+reg110; reg44=reg145*reg52; reg85=reg135*reg3;
    reg87=reg140-reg87; reg34=reg121+reg34; reg90=reg14*reg125; reg106=reg12*reg9; reg111=reg56*reg7;
    reg121=reg27*reg57; reg140=reg48*reg52; reg150=reg66*reg7; reg123=reg133+reg123; reg133=reg135*reg56;
    reg152=reg54*reg13; reg159=reg66*reg9; reg43=reg68-reg43; reg68=reg19*reg70; reg74=reg35*reg74;
    T reg166=reg29*reg25; reg111=reg106+reg111; reg106=reg136*reg101; T reg167=reg17*reg136; T reg168=reg156*reg110;
    T reg169=reg156*reg137; T reg170=(*f.m).alpha_2*reg6; T reg171=(*f.m).alpha_1*reg11; reg139=reg138+reg139; reg138=(*f.m).alpha_3*reg59;
    reg90=reg34+reg90; reg103=reg142+reg103; reg165=reg129+reg165; reg34=reg42*reg68; reg159=reg43-reg159;
    reg82=reg158+reg82; reg70=reg151*reg70; reg117=reg36+reg117; reg36=reg86*reg68; reg95=reg93+reg95;
    reg43=reg146*reg33; reg162=reg154+reg162; reg63=reg63*reg74; reg93=reg109*reg87; reg161=reg91+reg161;
    reg91=reg49*reg2; reg54=reg54*reg2; reg150=reg126-reg150; reg126=reg109*reg141; reg147=reg67+reg147;
    reg67=reg66*reg8; reg129=reg87*reg66; reg133=reg31+reg133; reg31=reg149*reg137; reg142=reg32*reg123;
    reg148=reg140+reg148; reg140=reg4*reg141; reg154=reg141*reg66; reg164=reg163+reg164; reg85=reg44+reg85;
    reg44=reg4*reg87; reg116=reg155+reg116; reg112=reg112*reg74; reg155=reg94*reg125; reg114=reg121+reg114;
    reg121=reg149*reg110; reg62=reg160+reg62; reg158=reg49*reg153; reg160=reg59*reg30; reg152=reg166+reg152;
    reg163=reg32*reg124; reg26=reg26*reg27; reg71=reg157+reg71; reg21=reg22*reg21; reg37=reg65*reg37;
    reg129=reg133+reg129; reg22=reg54*reg55; reg133=reg49*reg91; reg154=reg164+reg154; reg157=reg86*reg70;
    reg63=reg162+reg63; reg162=reg6*reg1; reg164=reg150*reg56; reg36=reg117+reg36; reg117=reg38*reg1;
    reg166=reg78*reg43; T reg172=reg159*reg12; reg158=reg62+reg158; reg62=reg104*reg68; T reg173=reg14*reg43;
    reg34=reg82+reg34; reg160=reg152+reg160; reg82=reg49*reg74; reg152=reg42*reg70; reg112=reg71+reg112;
    reg111=reg67-reg111; reg67=reg46*reg91; reg44=reg85+reg44; reg71=reg46*reg54; reg140=reg148+reg140;
    reg33=reg114*reg33; reg85=reg150*reg3; reg148=reg159*reg52; T reg174=reg91*reg55; T reg175=reg17*reg110;
    reg106=reg168-reg106; reg168=reg101*reg137; reg167=reg169-reg167; reg121=reg163+reg121; reg90=reg165*reg90;
    reg163=(*f.m).alpha_2*reg26; reg169=(*f.m).alpha_1*reg21; T reg176=(*f.m).alpha_3*reg37; reg170=reg171+reg170; reg138=reg103+reg138;
    reg31=reg142+reg31; reg95=reg165*reg95; reg103=reg139*reg32; reg142=reg136*reg149; reg126=reg147+reg126;
    reg155=reg116+reg155; reg116=reg49*reg54; reg57=reg65*reg57; reg93=reg161+reg93; reg147=0.78867513459481286553*elem.pos(0)[1];
    reg161=reg138*reg101; reg171=reg14*reg33; reg90=reg121+reg90; reg71=reg140+reg71; reg121=reg162*reg58;
    reg140=0.5*elem.pos(0)[1]; T reg177=0.5*elem.pos(1)[1]; reg176=reg170+reg176; reg170=reg162*reg38; reg67=reg44+reg67;
    reg44=reg78*reg33; T reg178=reg117*reg58; reg163=reg169+reg163; reg169=(*f.m).alpha_3*reg57; reg62=reg158+reg62;
    reg158=0.5*elem.pos(0)[0]; reg152=reg112+reg152; reg173=reg34+reg173; reg34=reg104*reg70; reg82=reg160+reg82;
    reg157=reg63+reg157; reg63=reg167*reg124; reg29=reg29*reg2; reg112=reg117*reg38; reg160=reg106*reg123;
    reg133=reg93+reg133; reg93=reg94*reg43; reg166=reg36+reg166; reg164=reg172+reg164; reg36=reg111*reg66;
    reg155=reg155*reg165; reg142=reg103+reg142; reg103=0.5*elem.pos(1)[0]; reg22=reg154+reg22; reg154=reg162*reg73;
    reg172=0.78867513459481286553*elem.pos(1)[0]; T reg179=0.21132486540518713447*elem.pos(1)[1]; T reg180=0.21132486540518713447*elem.pos(0)[1]; T reg181=0.78867513459481286553*elem.pos(1)[1]; T reg182=reg17*reg138;
    T reg183=0.21132486540518713447*elem.pos(0)[0]; reg116=reg126+reg116; reg168=reg175-reg168; reg126=0.78867513459481286553*elem.pos(0)[0]; reg175=reg4*reg111;
    reg85=reg148+reg85; reg148=reg26*reg0; T reg184=reg97*reg159; T reg185=reg117*reg73; reg174=reg129+reg174;
    reg129=reg50*reg150; T reg186=reg35*reg0; reg95=reg31+reg95; reg31=0.21132486540518713447*elem.pos(1)[0]; T reg187=reg156*reg124;
    T reg188=reg17*reg139; T reg189=reg156*reg123; T reg190=reg139*reg101; T reg191=reg139*reg168; reg160=reg63-reg160;
    reg93=reg62+reg93; reg34=reg82+reg34; reg62=reg94*reg33; reg169=reg163+reg169; reg171=reg152+reg171;
    reg161=reg90+reg161; reg173=reg176*reg173; reg182=reg95+reg182; reg166=reg166*reg176; reg155=reg142+reg155;
    reg63=reg138*reg156; reg82=reg148*reg35; reg170=reg116+reg170; reg90=reg109*reg111; reg129=reg184+reg129;
    reg95=reg31-reg183; reg116=0.78867513459481286553*elem.pos(2)[0]; reg142=0.78867513459481286553*elem.pos(2)[1]; reg152=reg181+reg180; reg163=reg186*reg99;
    reg178=reg67+reg178; reg67=reg148*reg99; reg180=reg179-reg180; reg121=reg71+reg121; reg183=reg183+reg172;
    reg71=reg46*reg29; reg175=reg85+reg175; reg85=reg186*reg84; reg44=reg157+reg44; reg157=reg11*reg1;
    reg185=reg174+reg185; reg36=reg164+reg36; reg164=reg29*reg55; reg174=reg148*reg84; reg154=reg22+reg154;
    reg22=reg177+reg140; reg140=reg177-reg140; reg177=0.21132486540518713447*elem.pos(2)[0]; reg184=0.5*elem.pos(2)[1]; T reg192=reg158+reg103;
    reg31=reg31+reg126; reg179=reg179+reg147; T reg193=0.21132486540518713447*elem.pos(2)[1]; T reg194=0.5*elem.pos(2)[0]; T reg195=reg186*reg35;
    reg112=reg133+reg112; reg158=reg103-reg158; reg171=reg169*reg171; reg67=reg121+reg67; reg103=reg157*reg58;
    reg71=reg175+reg71; reg179=reg193-reg179; reg85=reg185+reg85; reg31=reg177-reg31; reg166=reg182+reg166;
    reg44=reg44*reg169; reg174=reg154+reg174; reg63=reg155+reg63; reg93=reg176*reg93; reg121=reg139*reg110;
    reg190=reg187-reg190; reg133=reg136*reg124; reg62=reg34+reg62; reg34=reg139*reg137; reg188=reg189-reg188;
    reg154=reg136*reg123; reg147=reg181-reg147; reg191=reg160+reg191; reg155=reg17*reg124; reg126=reg172-reg126;
    reg163=reg178+reg163; reg160=reg101*reg123; reg173=reg161+reg173; reg161=0.5*elem.pos(3)[1]; reg95=reg95+reg116;
    reg22=reg184-reg22; reg172=0.78867513459481286553*elem.pos(3)[0]; reg158=reg158+reg194; reg152=reg142-reg152; reg175=0.21132486540518713447*elem.pos(3)[1];
    reg178=0.5*elem.pos(3)[0]; reg180=reg142+reg180; reg195=reg112+reg195; reg140=reg184+reg140; reg90=reg129+reg90;
    reg112=0.78867513459481286553*elem.pos(3)[1]; reg82=reg170+reg82; reg183=reg116-reg183; reg116=0.21132486540518713447*elem.pos(3)[0]; reg192=reg194-reg192;
    reg129=reg21*reg0; reg164=reg36+reg164; reg36=reg157*reg73; reg142=reg49*reg29; reg170=0.21132486540518713447*PNODE(0).dep[0];
    reg160=reg155-reg160; reg155=reg124*reg137; reg181=reg110*reg123; reg182=0.21132486540518713447*PNODE(1).dep[0]; reg121=reg133-reg121;
    reg190=reg190/reg191; reg133=reg157*reg38; reg184=0.5*vectors[0][indices[0]+0]; reg185=0.78867513459481286553*PNODE(1).dep[0]; reg187=0.5*vectors[0][indices[1]+0];
    reg189=0.21132486540518713447*PNODE(0).dep[1]; reg194=0.78867513459481286553*PNODE(1).dep[1]; T reg196=0.21132486540518713447*PNODE(1).dep[1]; reg31=reg172+reg31; reg22=reg22+reg161;
    reg171=reg173+reg171; reg158=reg158-reg178; reg179=reg112+reg179; reg44=reg166+reg44; reg161=reg140-reg161;
    reg192=reg178+reg192; reg140=reg129*reg84; reg36=reg164+reg36; reg62=reg169*reg62; reg183=reg183+reg116;
    reg112=reg180-reg112; reg103=reg71+reg103; reg71=reg129*reg99; reg93=reg63+reg93; reg152=reg152+reg175;
    reg172=reg95-reg172; reg177=reg126+reg177; reg63=0.78867513459481286553*PNODE(0).dep[1]; reg193=reg147+reg193; reg95=0.5*vectors[0][indices[1]+1];
    reg142=reg90+reg142; reg106=reg106/reg191; reg34=reg154-reg34; reg90=0.78867513459481286553*PNODE(0).dep[0]; reg126=reg163*reg82;
    reg147=reg67*reg195; reg154=reg174*reg195; reg164=reg85*reg82; reg188=reg188/reg191; reg167=reg167/reg191;
    reg166=0.5*vectors[0][indices[0]+1]; reg173=reg112*reg31; reg178=reg166+reg95; reg126=reg147-reg126; reg147=0.5*vectors[0][indices[2]+1];
    reg180=reg161*reg192; T reg197=reg187-reg184; T reg198=0.5*vectors[0][indices[2]+0]; reg184=reg187+reg184; reg166=reg95-reg166;
    reg140=reg36+reg140; reg36=reg112*reg183; reg71=reg103+reg71; reg95=reg172*reg152; reg103=reg63+reg196;
    reg116=reg177-reg116; reg175=reg193-reg175; reg167=reg167*reg171; reg106=reg106*reg44; reg177=0.21132486540518713447*PNODE(2).dep[1];
    reg187=reg90+reg182; reg188=reg188*reg171; reg193=0.21132486540518713447*PNODE(2).dep[0]; reg190=reg190*reg44; reg182=reg182-reg170;
    T reg199=0.78867513459481286553*PNODE(2).dep[0]; reg170=reg170+reg185; T reg200=0.78867513459481286553*PNODE(2).dep[1]; T reg201=reg189+reg194; reg189=reg196-reg189;
    reg196=reg172*reg179; reg34=reg34/reg191; reg121=reg121/reg191; reg168=reg168/reg191; reg160=reg160/reg191;
    reg181=reg155-reg181; reg155=reg85*reg67; T reg202=reg174*reg163; reg164=reg154-reg164; reg62=reg93+reg62;
    reg133=reg142+reg133; reg93=reg22*reg158; reg142=reg129*reg35; reg168=reg168*reg62; reg106=reg167-reg106;
    reg154=reg140*reg126; reg184=reg198-reg184; reg90=reg185-reg90; reg167=0.5*vectors[0][indices[3]+0]; reg198=reg197+reg198;
    reg201=reg200-reg201; reg185=reg71*reg164; reg197=0.21132486540518713447*PNODE(3).dep[1]; reg180=reg93-reg180; reg63=reg194-reg63;
    reg173=reg196-reg173; reg189=reg200+reg189; reg93=0.78867513459481286553*PNODE(3).dep[1]; reg160=reg160*reg62; reg188=reg190-reg188;
    reg171=reg34*reg171; reg44=reg121*reg44; reg142=reg133+reg142; reg155=reg202-reg155; reg36=reg95-reg36;
    reg34=reg183*reg175; reg95=0.21132486540518713447*PNODE(3).dep[0]; reg170=reg199-reg170; reg191=reg181/reg191; reg199=reg182+reg199;
    reg166=reg147+reg166; reg103=reg177-reg103; reg121=reg152*reg116; reg133=0.5*vectors[0][indices[3]+1]; reg181=0.78867513459481286553*PNODE(3).dep[0];
    reg187=reg193-reg187; reg178=reg147-reg178; reg63=reg177+reg63; reg147=1-(*f.m).resolution; reg201=reg201+reg197;
    reg177=reg71*reg82; reg182=reg163*reg142; reg190=reg174*reg142; reg194=reg172/reg36; reg196=reg85*reg142;
    reg200=reg71*reg195; reg202=reg142*reg155; reg103=reg103+reg93; reg160=reg188-reg160; reg187=reg187+reg181;
    reg185=reg154-reg185; reg62=reg191*reg62; reg44=reg171-reg44; reg181=reg199-reg181; reg93=reg189-reg93;
    reg154=reg175*reg31; reg171=reg116*reg179; reg188=reg67*reg142; reg189=reg140*reg195; reg191=reg183/reg36;
    reg199=reg152/reg36; reg184=reg167+reg184; reg106=reg168+reg106; reg167=reg198-reg167; reg178=reg178+reg133;
    reg168=reg112/reg36; reg90=reg193+reg90; reg133=reg166-reg133; reg166=reg140*reg82; reg192=reg192/reg180;
    reg158=reg158/reg180; reg161=reg161/reg180; reg180=reg22/reg180; reg170=reg170+reg95; reg34=reg121-reg34;
    reg172=reg172/reg173; reg22=reg31/reg173; reg112=reg112/reg173; reg121=reg179/reg173; reg193=reg181*reg191;
    reg198=reg194*reg170; reg188=reg177-reg188; reg196=reg189-reg196; reg182=reg200-reg182; reg177=reg32*(*f.m).resolution;
    reg160=reg160*reg147; reg154=reg171-reg154; reg171=reg158*reg184; reg189=reg187*reg172; reg200=reg181*reg22;
    T reg203=reg116/reg34; reg183=reg183/reg34; T reg204=reg192*reg167; reg152=reg152/reg34; T reg205=reg180*reg133;
    T reg206=reg161*reg178; T reg207=reg93*reg121; T reg208=reg103*reg112; T reg209=reg175/reg34; reg95=reg90-reg95;
    reg197=reg63-reg197; reg63=reg174*reg71; reg90=reg140*reg67; T reg210=reg85*reg71; T reg211=reg149*(*f.m).resolution;
    reg44=reg62+reg44; reg202=reg185+reg202; reg106=reg106*reg147; reg190=reg166-reg190; reg62=reg140*reg163;
    reg166=reg201*reg168; reg185=reg93*reg199; reg175=reg175/reg154; reg179=reg179/reg154; reg31=reg31/reg154;
    reg116=reg116/reg154; reg196=reg196/reg202; reg160=reg211+reg160; reg190=reg190/reg202; reg210=reg62-reg210;
    reg63=reg90-reg63; reg62=reg209*reg201; reg90=reg197*reg152; reg211=reg95*reg183; T reg212=reg203*reg170;
    reg167=reg180*reg167; reg184=reg161*reg184; reg208=reg207-reg208; reg200=reg189-reg200; reg206=reg205-reg206;
    reg204=reg171-reg204; reg178=reg158*reg178; reg133=reg192*reg133; reg171=reg93*reg22; reg106=reg177+reg106;
    reg44=reg44*reg147; reg177=reg103*reg172; reg189=reg187*reg112; reg205=reg181*reg121; reg181=reg181*reg199;
    reg207=reg138*(*f.m).resolution; T reg213=reg170*reg168; reg164=reg164/reg202; reg188=reg188/reg202; T reg214=reg194*reg201;
    reg93=reg93*reg191; reg182=reg182/reg202; reg193=reg198-reg193; reg166=reg185-reg166; reg126=reg126/reg202;
    reg185=reg187*reg116; reg210=reg210/reg202; reg189=reg205-reg189; reg62=reg90-reg62; reg90=reg95*reg31;
    reg198=reg197*reg179; reg211=reg212-reg211; reg171=reg177-reg171; reg177=reg103*reg175; reg32=reg32*(*f.m).deltaT;
    reg155=reg155/reg202; reg205=reg197*reg183; reg201=reg203*reg201; reg202=reg63/reg202; reg170=reg209*reg170;
    reg208=reg200+reg208; reg149=reg149*(*f.m).deltaT; reg63=reg95*reg152; reg124=reg147*reg124; reg110=reg147*reg110;
    reg123=reg147*reg123; reg137=reg147*reg137; reg166=reg193+reg166; reg139=reg139*reg147; reg93=reg214-reg93;
    reg136=reg136*reg147; reg213=reg181-reg213; reg206=reg204+reg206; reg184=reg167-reg184; elem.epsilon[0][0]=reg184;
    reg133=reg178-reg133; elem.epsilon[0][1]=reg133; reg160=reg160*(*f.m).deltaT; reg106=(*f.m).deltaT*reg106; reg44=reg207+reg44;
    reg167=(*f.m).resolution*reg126; reg178=(*f.m).resolution*reg182; reg181=(*f.m).resolution*reg164; reg193=(*f.m).resolution*reg196; reg200=(*f.m).resolution*reg190;
    reg204=(*f.m).resolution*reg188; reg197=reg197*reg31; reg181=reg110-reg181; reg124=reg167+reg124; reg93=reg93-reg160;
    reg110=(*f.m).resolution*reg202; reg17=reg17*reg147; reg166=0.5*reg166; reg167=(*f.m).resolution*reg210; reg170=reg63-reg170;
    reg103=reg103*reg116; reg101=reg147*reg101; reg206=0.5*reg206; elem.epsilon[0][2]=reg206; reg205=reg201-reg205;
    reg62=reg211+reg62; reg187=reg187*reg175; reg95=reg95*reg179; reg204=reg139+reg204; reg200=reg136-reg200;
    reg177=reg198-reg177; reg171=reg171-reg160; reg63=(*f.m).resolution*reg155; reg189=reg189-reg106; reg44=reg44*(*f.m).deltaT;
    reg136=reg184-reg32; reg139=reg133-reg149; reg193=reg137+reg193; reg178=reg123-reg178; reg213=reg213-reg106;
    reg156=reg147*reg156; reg208=0.5*reg208; reg90=reg185-reg90; reg138=reg138*(*f.m).deltaT; reg123=reg139*reg164;
    reg177=reg90+reg177; reg90=reg213*reg204; reg187=reg95-reg187; reg95=reg93*reg200; reg197=reg103-reg197;
    reg62=0.5*reg62; reg205=reg205-reg160; reg170=reg170-reg106; reg103=reg171*reg200; reg137=reg189*reg204;
    reg147=reg171*reg193; reg185=reg189*reg178; reg208=reg208-reg44; reg198=reg171*reg181; reg110=reg156+reg110;
    reg156=reg189*reg124; reg167=reg17-reg167; reg101=reg63+reg101; reg17=reg93*reg193; reg63=reg213*reg178;
    reg166=reg166-reg44; reg201=reg93*reg181; reg207=reg136*reg126; reg211=reg206-reg138; reg212=reg213*reg124;
    reg214=reg139*reg196; T reg215=reg136*reg182; T reg216=reg211*reg155; T reg217=reg205*reg181; reg177=0.5*reg177;
    reg215=reg214-reg215; reg214=reg211*reg210; reg62=reg62-reg44; T reg218=reg170*reg178; T reg219=reg205*reg193;
    reg95=reg90+reg95; reg197=reg197-reg160; reg90=reg170*reg204; T reg220=reg205*reg200; reg187=reg187-reg106;
    reg17=reg63+reg17; reg63=reg166*reg167; T reg221=reg166*reg110; reg147=reg185+reg147; reg185=reg208*reg167;
    reg198=reg156+reg198; reg103=reg137+reg103; reg137=reg208*reg110; reg156=reg208*reg101; reg201=reg212+reg201;
    reg212=reg166*reg101; T reg222=reg170*reg124; reg123=reg207-reg123; reg207=reg62*reg101; reg137=reg103+reg137;
    reg217=reg222+reg217; reg219=reg218+reg219; reg103=reg62*reg167; reg156=reg198+reg156; reg185=reg147+reg185;
    reg220=reg90+reg220; reg90=reg62*reg110; reg212=reg201+reg212; reg177=reg177-reg44; reg147=reg187*reg178;
    reg198=reg197*reg193; reg214=reg215-reg214; reg63=reg17+reg63; reg17=reg197*reg181; reg201=reg187*reg204;
    reg215=reg197*reg200; reg216=reg123+reg216; reg95=reg221+reg95; reg123=reg187*reg124; reg218=reg177*reg101;
    reg212=reg213*reg212; reg156=reg189*reg156; reg198=reg147+reg198; reg147=reg216+reg214; reg189=reg177*reg167;
    reg213=reg177*reg110; reg63=reg93*reg63; reg207=reg217+reg207; reg95=2*reg95; reg185=reg171*reg185;
    reg17=reg123+reg17; reg103=reg219+reg103; reg137=2*reg137; reg90=reg220+reg90; reg215=reg201+reg215;
    reg137=reg208*reg137; reg189=reg198+reg189; reg185=reg156+reg185; reg213=reg215+reg213; reg218=reg17+reg218;
    reg136=reg136*reg188; reg166=reg95*reg166; reg90=2*reg90; reg139=reg139*reg190; reg103=reg205*reg103;
    reg207=reg170*reg207; reg147=reg147/3; reg63=reg212+reg63; reg90=reg62*reg90; reg213=2*reg213;
    reg139=reg136-reg139; reg103=reg207+reg103; reg214=reg214-reg147; reg137=reg185+reg137; reg133=reg133-reg160;
    reg63=reg166+reg63; reg189=reg197*reg189; reg211=reg211*reg202; reg218=reg187*reg218; reg184=reg184-reg106;
    reg216=reg216-reg147; reg17=reg184*reg124; reg63=reg63*reg36; reg213=reg177*reg213; reg216=pow(reg216,2);
    reg90=reg103+reg90; reg206=reg206-reg44; reg214=pow(reg214,2); reg139=reg211+reg139; reg137=reg137*reg173;
    reg189=reg218+reg189; reg62=reg133*reg181; reg93=reg133*reg193; reg95=reg184*reg178; reg103=2*reg139;
    reg62=reg17+reg62; reg63=0.25*reg63; reg17=reg206*reg167; reg93=reg95+reg93; reg95=reg206*reg101;
    reg147=pow(reg147,2); reg214=reg216+reg214; reg137=0.25*reg137; reg90=reg90*reg34; reg184=reg184*reg204;
    reg133=reg133*reg200; reg213=reg189+reg213; reg90=0.25*reg90; reg137=reg63+reg137; reg17=reg93+reg17;
    elem.sigma[0][1]=reg17; reg133=reg184+reg133; reg206=reg206*reg110; reg147=reg214+reg147; reg103=reg139*reg103;
    reg213=reg213*reg154; reg95=reg62+reg95; elem.sigma[0][0]=reg95; reg213=0.25*reg213; reg62=reg95*reg52;
    reg63=reg17*reg3; reg90=reg137+reg90; reg103=reg147+reg103; reg93=reg97*reg95; reg123=reg50*reg17;
    reg17=reg17*reg56; reg95=reg95*reg12; reg206=reg133+reg206; elem.sigma[0][2]=reg206; reg63=reg62+reg63;
    reg213=reg90+reg213; reg62=reg46*reg206; reg123=reg93+reg123; reg103=1.5*reg103; reg90=reg49*reg206;
    reg206=reg206*reg55; reg17=reg95+reg17; elem.sigma_von_mises=pow(reg103,0.5); elem.sigma_local[0][0]=reg17+reg206; elem.sigma_local[0][2]=reg123+reg90;
    elem.sigma_local[0][1]=reg63+reg62; elem.ener=reg213/2;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; reg1=1.0/reg1; T reg2=2*(*f.m).shear_modulus_12;
    T reg3=reg1*reg0; reg2=1.0/reg2; T reg4=reg2*reg3; T reg5=pow((*f.m).v2[1],2); T reg6=pow((*f.m).v2[0],2);
    T reg7=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg8=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=pow((*f.m).v1[1],2); T reg11=pow((*f.m).v1[0],2);
    reg10=reg11+reg10; reg11=pow((*f.m).v1[2],2); T reg12=reg9*reg4; T reg13=reg7*reg4; T reg14=reg8*reg4;
    T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=pow((*f.m).v2[2],2); reg5=reg6+reg5; reg6=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg17=reg6*reg12;
    T reg18=reg7*reg13; T reg19=reg15*reg12; reg16=reg5+reg16; reg11=reg10+reg11; reg5=reg7*reg14;
    reg10=reg6*reg13; T reg20=1.0/(*f.m).elastic_modulus_1; reg11=pow(reg11,0.5); reg5=reg17+reg5; T reg21=reg15*reg14;
    reg18=reg19-reg18; reg16=pow(reg16,0.5); reg19=reg10+reg21; T reg22=reg20*reg18; T reg23=reg6*reg5;
    T reg24=(*f.m).v1[1]/reg11; T reg25=(*f.m).v1[2]/reg11; T reg26=(*f.m).v2[1]/reg16; T reg27=(*f.m).v2[2]/reg16; reg11=(*f.m).v1[0]/reg11;
    reg23=reg22-reg23; reg16=(*f.m).v2[0]/reg16; reg22=reg8*reg19; T reg28=reg8*reg3; T reg29=reg9*reg3;
    T reg30=reg15*reg4; reg3=reg7*reg3; T reg31=reg25*reg26; T reg32=reg8*reg13; T reg33=reg24*reg27;
    reg4=reg6*reg4; T reg34=reg2*reg0; reg12=reg20*reg12; T reg35=reg8*reg14; T reg36=reg8*reg30;
    T reg37=reg15*reg29; T reg38=2*reg16; T reg39=reg8*reg34; reg32=reg17+reg32; reg35=reg12-reg35;
    reg12=reg2*reg1; reg13=reg20*reg13; reg29=reg6*reg29; reg17=reg7*reg3; T reg40=reg7*reg28;
    T reg41=reg11*reg27; T reg42=reg25*reg16; reg22=reg23-reg22; reg23=reg9*reg34; T reg43=reg33-reg31;
    T reg44=reg8*reg4; T reg45=2*reg11; reg34=reg7*reg34; reg14=reg6*reg14; reg4=reg6*reg4;
    reg14=reg13+reg14; reg3=reg6*reg3; T reg46=reg7*reg12; reg28=reg15*reg28; T reg47=reg8*reg12;
    reg17=reg37-reg17; reg40=reg29+reg40; reg29=reg15*reg23; reg23=reg6*reg23; reg37=2*reg43;
    T reg48=reg7*reg34; T reg49=reg7*reg39; reg12=reg9*reg12; T reg50=reg11*reg26; T reg51=reg24*reg16;
    reg18=reg18/reg22; T reg52=reg42-reg41; T reg53=reg38*reg26; T reg54=pow(reg26,2); reg30=reg20*reg30;
    reg32=reg32/reg22; reg36=reg10+reg36; reg5=reg5/reg22; T reg55=pow(reg11,2); reg35=reg35/reg22;
    T reg56=pow(reg24,2); T reg57=pow(reg16,2); T reg58=reg24*reg45; reg44=reg13+reg44; reg49=reg23+reg49;
    reg13=reg15*reg12; reg12=reg6*reg12; reg48=reg29-reg48; reg23=reg7*reg46; reg29=reg7*reg47;
    T reg59=reg56*reg18; T reg60=reg54*reg5; T reg61=reg58*reg18; T reg62=reg53*reg5; T reg63=reg55*reg32;
    T reg64=reg57*reg35; T reg65=reg56*reg32; T reg66=reg54*reg35; T reg67=reg58*reg32; T reg68=reg53*reg35;
    reg14=reg14/reg22; reg44=reg44/reg22; T reg69=reg57*reg5; reg4=reg30-reg4; reg30=reg55*reg18;
    reg19=reg19/reg22; reg34=reg6*reg34; reg36=reg36/reg22; reg39=reg15*reg39; T reg70=pow(reg27,2);
    T reg71=pow(reg43,2); reg17=reg20*reg17; T reg72=reg50-reg51; reg40=reg6*reg40; T reg73=reg3+reg28;
    T reg74=pow(reg52,2); T reg75=pow(reg25,2); T reg76=reg37*reg52; reg64=reg63+reg64; reg63=reg53*reg44;
    T reg77=reg58*reg36; T reg78=reg71*reg14; reg47=reg15*reg47; reg66=reg65+reg66; reg65=reg74*reg14;
    T reg79=reg75*reg32; T reg80=reg54*reg44; T reg81=reg56*reg36; reg40=reg17-reg40; reg17=reg70*reg35;
    T reg82=reg57*reg44; T reg83=reg55*reg36; reg73=reg8*reg73; T reg84=reg34+reg39; reg48=reg20*reg48;
    reg68=reg67+reg68; reg67=reg76*reg14; reg49=reg6*reg49; reg4=reg4/reg22; T reg85=reg75*reg18;
    T reg86=reg70*reg5; reg46=reg6*reg46; T reg87=reg74*reg19; reg69=reg30+reg69; reg60=reg59+reg60;
    reg62=reg61+reg62; reg30=reg76*reg19; reg59=reg71*reg19; reg23=reg13-reg23; reg29=reg12+reg29;
    reg12=pow(reg72,2); reg23=reg20*reg23; reg67=reg68+reg67; reg87=reg60+reg87; reg59=reg69+reg59;
    reg13=reg25*reg45; reg84=reg8*reg84; reg86=reg85+reg86; reg60=reg12*reg14; reg61=reg12*reg19;
    reg29=reg6*reg29; reg17=reg79+reg17; reg78=reg64+reg78; reg64=2*reg24; reg68=reg46+reg47;
    reg30=reg62+reg30; reg65=reg66+reg65; reg63=reg77+reg63; reg62=reg11*reg24; reg66=reg70*reg44;
    reg69=reg16*reg26; reg77=reg75*reg36; reg76=reg76*reg4; reg79=reg74*reg4; reg80=reg81+reg80;
    reg81=reg71*reg4; reg82=reg83+reg82; reg73=reg40-reg73; reg40=reg38*reg27; reg83=2*reg26;
    reg49=reg48-reg49; reg48=reg57*reg67; reg85=reg55*reg59; reg68=reg8*reg68; T reg88=reg55*reg30;
    T reg89=reg69*reg67; T reg90=reg62*reg30; reg67=reg54*reg67; T reg91=reg69*reg65; T reg92=reg62*reg87;
    T reg93=reg57*reg65; reg30=reg56*reg30; T reg94=reg83*reg27; reg37=reg37*reg72; T reg95=reg57*reg78;
    reg29=reg23-reg29; reg23=2*reg52; T reg96=reg11*reg16; T reg97=reg25*reg64; T reg98=reg24*reg26;
    reg50=reg51+reg50; reg84=reg49-reg84; reg73=reg73/reg22; reg49=reg55*reg87; reg51=reg24*reg43;
    T reg99=reg11*reg52; reg65=reg54*reg65; reg87=reg56*reg87; T reg100=reg26*reg64; T reg101=2*reg25;
    T reg102=reg52*reg43; T reg103=reg16*reg45; reg76=reg63+reg76; reg63=reg12*reg4; reg66=reg77+reg66;
    reg79=reg80+reg79; reg81=reg82+reg81; reg77=reg40*reg35; reg80=reg13*reg32; reg60=reg17+reg60;
    reg17=reg54*reg78; reg82=reg56*reg59; reg61=reg86+reg61; reg86=reg40*reg5; T reg104=reg13*reg18;
    reg68=reg29-reg68; reg23=reg23*reg72; reg29=reg13*reg36; T reg105=reg40*reg44; T reg106=reg71*reg76;
    T reg107=reg54*reg6; T reg108=reg57*reg20; reg48=reg88+reg48; reg88=reg54*reg15; T reg109=reg54*reg60;
    T reg110=reg74*reg76; reg67=reg30+reg67; reg76=reg102*reg76; reg89=reg90+reg89; reg30=reg96*reg73;
    reg90=reg98*reg73; reg101=reg27*reg101; T reg111=reg50*reg73; T reg112=reg56*reg61; T reg113=reg57*reg60;
    T reg114=reg55*reg61; T reg115=reg102*reg79; reg91=reg92+reg91; reg92=reg57*reg6; T reg116=reg26*reg43;
    T reg117=reg71*reg79; reg93=reg49+reg93; reg59=reg62*reg59; reg78=reg69*reg78; reg95=reg85+reg95;
    reg49=reg71*reg81; reg85=reg103*reg20; T reg118=reg16*reg52; T reg119=reg100*reg15; T reg120=reg100*reg6;
    reg17=reg82+reg17; reg82=reg74*reg81; reg5=reg94*reg5; reg77=reg80+reg77; reg84=reg84/reg22;
    reg80=reg37*reg14; reg99=reg51+reg99; reg32=reg97*reg32; reg35=reg94*reg35; reg63=reg66+reg63;
    reg86=reg104+reg86; reg51=reg25*reg27; reg66=reg37*reg19; reg18=reg97*reg18; reg79=reg74*reg79;
    reg104=reg11*reg43; reg65=reg87+reg65; reg87=reg103*reg6; T reg121=reg24*reg52; reg107=reg108-reg107;
    reg106=reg48+reg106; reg48=reg103*reg111; reg117=reg93+reg117; reg93=reg103*reg90; reg108=reg70*reg8;
    T reg122=reg100*reg90; reg79=reg65+reg79; reg64=reg52*reg64; reg45=reg43*reg45; reg65=reg100*reg30;
    reg110=reg67+reg110; reg120=reg85-reg120; reg67=reg101*reg8; reg85=reg100*reg111; T reg123=reg71*reg63;
    reg113=reg114+reg113; reg82=reg17+reg82; reg60=reg69*reg60; reg17=reg25*reg72; reg114=reg16*reg43;
    T reg124=reg103*reg30; reg49=reg95+reg49; reg22=reg68/reg22; reg68=reg26*reg52; reg66=reg86+reg66;
    reg5=reg18+reg5; reg19=reg23*reg19; reg118=reg116+reg118; reg87=reg119-reg87; reg18=reg101*reg7;
    reg86=reg99*reg84; reg80=reg77+reg80; reg35=reg32+reg35; reg14=reg23*reg14; reg105=reg29+reg105;
    reg37=reg37*reg4; reg36=reg97*reg36; reg44=reg94*reg44; reg92=reg88-reg92; reg29=reg70*reg7;
    reg32=reg51*reg73; reg104=reg104*reg84; reg121=reg121*reg84; reg90=reg50*reg90; reg77=reg54*reg7;
    reg115=reg91+reg115; reg76=reg89+reg76; reg111=reg50*reg111; reg88=reg57*reg8; reg89=reg103*reg8;
    reg61=reg62*reg61; reg78=reg59+reg78; reg81=reg102*reg81; reg59=reg100*reg7; reg91=reg74*reg63;
    reg109=reg112+reg109; reg29=reg92-reg29; reg77=reg88+reg77; reg88=reg70*reg9; reg108=reg107-reg108;
    reg92=reg56*reg6; reg101=reg101*reg9; reg59=reg89+reg59; reg15=reg56*reg15; reg122=reg79+reg122;
    reg79=reg64*reg121; reg6=reg55*reg6; reg18=reg87-reg18; reg67=reg120-reg67; reg87=reg27*reg72;
    reg85=reg110+reg85; reg89=reg64*reg86; reg95=(*f.m).alpha_2*reg54; reg107=(*f.m).alpha_1*reg56; reg110=(*f.m).alpha_1*reg55;
    reg112=(*f.m).alpha_2*reg57; reg91=reg109+reg91; reg109=reg11*reg72; reg116=reg25*reg43; reg119=reg100*reg32;
    reg42=reg41+reg42; reg41=reg54*reg80; reg120=reg56*reg66; T reg125=reg57*reg80; T reg126=reg55*reg66;
    T reg127=reg45*reg86; reg48=reg106+reg48; reg106=reg103*reg32; reg123=reg113+reg123; reg113=reg45*reg121;
    reg93=reg117+reg93; reg81=reg78+reg81; reg30=reg50*reg30; reg90=reg115+reg90; reg121=reg99*reg121;
    reg111=reg76+reg111; reg86=reg99*reg86; reg60=reg61+reg60; reg63=reg102*reg63; reg61=reg45*reg104;
    reg124=reg49+reg124; reg20=reg55*reg20; reg49=reg118*reg22; reg19=reg5+reg19; reg68=reg68*reg22;
    reg114=reg114*reg22; reg14=reg35+reg14; reg37=reg105+reg37; reg44=reg36+reg44; reg4=reg23*reg4;
    reg17=reg17*reg84; reg65=reg82+reg65; reg5=reg64*reg104; reg83=reg83*reg52; reg38=reg38*reg43;
    reg23=reg98*reg18; reg119=reg91+reg119; reg80=reg69*reg80; reg66=reg62*reg66; reg35=reg118*reg49;
    reg86=reg111+reg86; reg36=reg56*reg7; reg76=reg118*reg68; reg121=reg90+reg121; reg78=reg55*reg19;
    reg77=reg88-reg77; reg82=reg108*reg57; reg104=reg99*reg104; reg59=reg101-reg59; reg30=reg81+reg30;
    reg81=reg18*reg54; reg88=reg67*reg57; reg90=(*f.m).alpha_2*reg70; reg91=reg29*reg54; reg101=reg83*reg68;
    reg79=reg122+reg79; reg7=reg75*reg7; reg6=reg15-reg6; reg15=reg42*reg73; reg4=reg44+reg4;
    reg44=reg75*reg8; reg92=reg20-reg92; reg8=reg55*reg8; reg20=reg16*reg72; reg105=reg27*reg43;
    reg87=reg87*reg22; reg111=reg24*reg72; reg115=reg25*reg52; reg109=reg116+reg109; reg61=reg124+reg61;
    reg116=reg38*reg114; reg32=reg50*reg32; reg33=reg31+reg33; reg63=reg60+reg63; reg31=reg57*reg14;
    reg60=reg96*reg67; reg117=reg64*reg17; reg89=reg85+reg89; reg85=reg83*reg49; reg122=reg96*reg108;
    reg124=reg98*reg29; reg41=reg120+reg41; reg18=reg18*reg56; reg67=reg67*reg55; reg120=reg74*reg37;
    reg127=reg48+reg127; reg49=reg38*reg49; reg29=reg29*reg56; reg108=reg108*reg55; reg48=reg56*reg19;
    T reg128=reg54*reg14; T reg129=(*f.m).alpha_3*reg71; reg125=reg126+reg125; reg126=reg71*reg37; T reg130=reg83*reg114;
    reg5=reg65+reg5; reg112=reg110+reg112; reg95=reg107+reg95; reg68=reg38*reg68; reg106=reg123+reg106;
    reg65=reg45*reg17; reg113=reg93+reg113; reg93=(*f.m).alpha_3*reg74; reg107=(*f.m).alpha_1*reg75; reg68=reg113+reg68;
    reg110=reg69*reg2; reg111=reg115+reg111; reg130=reg5+reg130; reg29=reg108+reg29; reg5=reg77*reg75;
    reg108=reg51*reg59; reg129=reg112+reg129; reg71=reg71*reg4; reg112=reg51*reg77; reg44=reg92-reg44;
    reg92=reg50*reg2; reg101=reg79+reg101; reg128=reg48+reg128; reg20=reg105+reg20; reg48=reg27*reg52;
    reg79=reg26*reg72; reg74=reg74*reg4; reg105=reg83*reg87; reg109=reg109*reg84; reg16=reg16*reg27;
    reg11=reg11*reg25; reg77=reg70*reg77; reg91=reg82+reg91; reg7=reg6-reg7; reg93=reg95+reg93;
    reg126=reg125+reg126; reg73=reg33*reg73; reg6=reg103*reg15; reg12=(*f.m).alpha_3*reg12; reg82=reg38*reg87;
    reg124=reg122+reg124; reg114=reg118*reg114; reg104=reg30+reg104; reg30=(*f.m).alpha_1*reg62; reg36=reg8+reg36;
    reg9=reg75*reg9; reg85=reg89+reg85; reg8=(*f.m).alpha_2*reg69; reg65=reg106+reg65; reg35=reg86+reg35;
    reg86=reg59*reg75; reg80=reg66+reg80; reg37=reg102*reg37; reg23=reg60+reg23; reg49=reg127+reg49;
    reg116=reg61+reg116; reg17=reg99*reg17; reg32=reg63+reg32; reg90=reg107+reg90; reg60=reg100*reg15;
    reg31=reg78+reg31; reg120=reg41+reg120; reg76=reg121+reg76; reg14=reg69*reg14; reg19=reg62*reg19;
    reg117=reg119+reg117; reg81=reg88+reg81; reg59=reg70*reg59; reg18=reg67+reg18; reg41=reg50*reg110;
    reg108=reg23+reg108; reg112=reg124+reg112; reg23=reg50*reg92; reg105=reg117+reg105; reg103=reg103*reg73;
    reg71=reg31+reg71; reg31=reg45*reg109; reg6=reg126+reg6; reg82=reg65+reg82; reg20=reg20*reg22;
    reg84=reg111*reg84; reg4=reg102*reg4; reg14=reg19+reg14; reg15=reg50*reg15; reg37=reg80+reg37;
    reg27=reg26*reg27; reg25=reg24*reg25; reg43=reg72*reg43; reg19=reg93*reg101; reg24=reg129*reg130;
    reg26=reg93*reg68; reg61=reg129*reg116; reg63=reg76*reg49; reg65=reg85*reg76; reg66=reg35*reg68;
    reg67=reg35*reg101; reg69=(*f.m).alpha_2*reg16; reg78=(*f.m).alpha_1*reg11; reg102=(*f.m).alpha_3*reg102; reg8=reg30+reg8;
    reg12=reg90+reg12; reg36=reg9-reg36; reg77=reg91+reg77; reg79=reg48+reg79; reg74=reg128+reg74;
    reg100=reg100*reg73; reg9=reg53*reg110; reg30=reg7*reg54; reg48=reg44*reg57; reg87=reg118*reg87;
    reg17=reg32+reg17; reg32=reg92*reg58; reg86=reg18+reg86; reg16=reg16*reg1; reg18=reg42*reg1;
    reg114=reg104+reg114; reg59=reg81+reg59; reg92=reg53*reg92; reg110=reg110*reg58; reg5=reg29+reg5;
    reg29=reg44*reg55; reg60=reg120+reg60; reg80=reg64*reg109; reg81=reg7*reg56; reg88=reg27*reg0;
    reg63=reg66-reg63; reg66=reg85*reg68; reg80=reg60+reg80; reg60=reg49*reg101; reg89=reg18*reg40;
    reg90=reg83*reg20; reg100=reg74+reg100; reg41=reg112+reg41; reg74=reg16*reg42; reg102=reg8+reg102;
    reg87=reg17+reg87; reg64=reg64*reg84; reg69=reg78+reg69; reg43=(*f.m).alpha_3*reg43; reg73=reg50*reg73;
    reg105=reg12*reg105; reg4=reg14+reg4; reg23=reg108+reg23; reg8=reg18*reg42; reg19=reg24+reg19;
    reg22=reg79*reg22; reg81=reg29+reg81; reg75=reg36*reg75; reg31=reg6+reg31; reg110=reg5+reg110;
    reg5=reg16*reg13; reg27=(*f.m).alpha_2*reg27; reg6=(*f.m).alpha_1*reg25; reg14=reg38*reg20; reg82=reg12*reg82;
    reg32=reg86+reg32; reg18=reg18*reg13; reg26=reg61+reg26; reg30=reg48+reg30; reg70=reg70*reg36;
    reg7=reg98*reg7; reg44=reg96*reg44; reg92=reg59+reg92; reg109=reg99*reg109; reg2=reg62*reg2;
    reg65=reg67-reg65; reg45=reg45*reg84; reg16=reg16*reg40; reg17=reg114*reg129; reg24=reg76*reg93;
    reg9=reg77+reg9; reg29=reg33*reg0; reg103=reg71+reg103; reg52=reg72*reg52; reg15=reg37+reg15;
    reg60=reg66-reg60; reg82=reg26+reg82; reg75=reg81+reg75; reg26=reg2*reg58; reg5=reg110+reg5;
    reg37=reg88*reg97; reg18=reg32+reg18; reg32=reg29*reg97; reg70=reg30+reg70; reg30=reg53*reg2;
    reg16=reg9+reg16; reg9=reg88*reg94; reg89=reg92+reg89; reg48=reg29*reg94; reg1=reg11*reg1;
    reg84=reg99*reg84; reg73=reg4+reg73; reg20=reg118*reg20; reg109=reg15+reg109; reg90=reg80+reg90;
    reg74=reg41+reg74; reg88=reg88*reg33; reg64=reg100+reg64; reg83=reg83*reg22; reg105=reg19+reg105;
    reg36=reg51*reg36; reg8=reg23+reg8; reg7=reg44+reg7; reg4=reg85*reg102; reg29=reg29*reg33;
    reg11=reg63*reg130; reg15=reg65*reg116; reg24=reg17+reg24; reg38=reg38*reg22; reg17=reg102*reg49;
    reg14=reg31+reg14; reg52=(*f.m).alpha_3*reg52; reg27=reg6+reg27; reg43=reg69+reg43; reg45=reg103+reg45;
    reg12=reg87*reg12; reg12=reg24+reg12; reg0=reg25*reg0; reg2=reg50*reg2; reg32=reg18+reg32;
    reg48=reg89+reg48; reg52=reg27+reg52; reg26=reg75+reg26; reg13=reg1*reg13; reg36=reg7+reg36;
    reg6=reg102*reg35; reg30=reg70+reg30; reg40=reg1*reg40; reg37=reg5+reg37; reg4=reg105+reg4;
    reg90=reg90*reg43; reg9=reg16+reg9; reg5=reg85*reg114; reg7=reg35*reg116; reg38=reg45+reg38;
    reg16=reg114*reg49; reg11=reg15-reg11; reg20=reg109+reg20; reg83=reg64+reg83; reg88=reg74+reg88;
    reg84=reg73+reg84; reg22=reg118*reg22; reg15=reg35*reg130; reg29=reg8+reg29; reg14=reg43*reg14;
    reg8=reg114*reg60; reg17=reg82+reg17; reg6=reg12+reg6; reg12=reg85*reg116; reg18=reg76*reg130;
    reg20=reg43*reg20; reg19=reg49*reg130; reg38=reg52*reg38; reg14=reg17+reg14; reg17=1-var_inter[0];
    reg97=reg0*reg97; reg13=reg26+reg13; reg23=1-var_inter[1]; reg8=reg11+reg8; reg42=reg1*reg42;
    reg2=reg36+reg2; reg1=reg76*reg116; reg22=reg84+reg22; reg16=reg7-reg16; reg7=reg114*reg101;
    reg90=reg4+reg90; reg83=reg83*reg52; reg4=reg114*reg68; reg11=reg48*reg88; reg24=reg9*reg29;
    reg94=reg0*reg94; reg40=reg30+reg40; reg5=reg15-reg5; reg15=reg37*reg29; reg25=reg32*reg88;
    reg25=reg15-reg25; reg15=reg23*elem.pos(0)[1]; reg26=reg23*elem.pos(1)[1]; reg65=reg65/reg8; reg27=reg37*reg48;
    reg5=reg5/reg8; reg63=reg63/reg8; reg7=reg18-reg7; reg11=reg24-reg11; reg18=reg32*reg9;
    reg20=reg6+reg20; reg22=reg52*reg22; reg16=reg16/reg8; reg6=reg23*elem.pos(1)[0]; reg24=reg23*elem.pos(0)[0];
    reg30=elem.pos(1)[0]*var_inter[0]; reg31=reg17*elem.pos(0)[0]; reg4=reg1-reg4; reg1=reg116*reg101; reg19=reg12-reg19;
    reg12=elem.pos(1)[1]*var_inter[0]; reg83=reg90+reg83; reg97=reg13+reg97; reg94=reg40+reg94; reg42=reg2+reg42;
    reg33=reg0*reg33; reg38=reg14+reg38; reg0=reg17*elem.pos(0)[1]; reg2=reg68*reg130; reg19=reg19/reg8;
    reg33=reg42+reg33; reg13=elem.pos(2)[1]*var_inter[0]; reg14=reg12+reg0; reg2=reg1-reg2; reg60=reg60/reg8;
    reg1=reg97*reg11; reg7=reg7/reg8; reg4=reg4/reg8; reg63=reg63*reg83; reg26=reg26-reg15;
    reg65=reg65*reg38; reg36=reg94*reg25; reg40=elem.pos(2)[1]*var_inter[1]; reg18=reg27-reg18; reg5=reg5*reg38;
    reg6=reg6-reg24; reg16=reg16*reg83; reg27=reg31+reg30; reg41=elem.pos(2)[0]*var_inter[0]; reg42=elem.pos(2)[0]*var_inter[1];
    reg22=reg20+reg22; reg13=reg13-reg14; reg20=elem.pos(3)[1]*reg17; reg36=reg1-reg36; reg1=reg33*reg18;
    reg43=reg94*reg88; reg41=reg41-reg27; reg44=elem.pos(3)[0]*reg17; reg45=reg9*reg33; reg51=reg37*reg33;
    reg88=reg97*reg88; reg63=reg65-reg63; reg19=reg19*reg22; reg60=reg60*reg22; reg52=elem.pos(3)[1]*var_inter[1];
    reg8=reg2/reg8; reg26=reg40+reg26; reg5=reg16-reg5; reg42=reg6+reg42; reg2=elem.pos(3)[0]*var_inter[1];
    reg83=reg4*reg83; reg38=reg7*reg38; reg4=reg97*reg29; reg83=reg38-reg83; reg22=reg8*reg22;
    reg26=reg26-reg52; reg45=reg43-reg45; reg19=reg5-reg19; reg5=reg32*reg33; reg33=reg48*reg33;
    reg29=reg94*reg29; reg51=reg88-reg51; reg1=reg36+reg1; reg63=reg60+reg63; reg20=reg13+reg20;
    reg42=reg42-reg2; reg37=reg37*reg94; reg9=reg97*reg9; reg44=reg41+reg44; reg6=1-(*f.m).resolution;
    reg7=reg26*reg44; reg8=reg42*reg20; reg37=reg9-reg37; reg9=reg93*(*f.m).resolution; reg83=reg22+reg83;
    reg45=reg45/reg1; reg13=reg129*(*f.m).resolution; reg63=reg63*reg6; reg33=reg29-reg33; reg94=reg32*reg94;
    reg48=reg97*reg48; reg5=reg4-reg5; reg51=reg51/reg1; reg19=reg19*reg6; reg33=reg33/reg1;
    reg11=reg11/reg1; reg4=reg102*(*f.m).resolution; reg16=(*f.m).resolution*reg45; reg22=(*f.m).resolution*reg51; reg25=reg25/reg1;
    reg19=reg9+reg19; reg5=reg5/reg1; reg83=reg83*reg6; reg63=reg13+reg63; reg7=reg8-reg7;
    reg76=reg76*reg6; reg94=reg48-reg94; reg114=reg114*reg6; reg37=reg37/reg1; reg130=reg6*reg130;
    reg68=reg6*reg68; reg35=reg6*reg35; reg101=reg6*reg101; reg116=reg6*reg116; reg8=(*f.m).resolution*reg25;
    reg94=reg94/reg1; reg1=reg18/reg1; reg9=(*f.m).resolution*reg11; reg26=reg26/reg7; reg22=reg76-reg22;
    reg16=reg114+reg16; reg19=reg19*(*f.m).deltaT; reg13=(*f.m).resolution*reg37; reg63=(*f.m).deltaT*reg63; reg18=(*f.m).resolution*reg5;
    reg29=(*f.m).resolution*reg33; reg83=reg4+reg83; reg44=reg44/reg7; reg42=reg42/reg7; reg20=reg20/reg7;
    reg4=reg23*reg44; reg32=reg17*reg42; reg36=var_inter[0]*reg26; reg38=var_inter[0]*reg42; reg83=reg83*(*f.m).deltaT;
    reg40=reg23*reg20; reg41=reg17*reg26; reg43=(*f.m).resolution*reg1; reg85=reg85*reg6; reg48=reg63*reg16;
    reg59=reg19*reg22; reg13=reg35+reg13; reg18=reg101+reg18; reg29=reg130-reg29; reg8=reg68-reg8;
    reg116=reg9+reg116; reg9=(*f.m).resolution*reg94; reg49=reg6*reg49; reg6=var_inter[1]*reg20; reg35=var_inter[1]*reg44;
    reg60=reg83*reg13; reg61=reg48+reg59; reg62=reg35+reg32; reg64=reg4+reg38; reg65=reg63*reg29;
    reg66=reg19*reg18; reg67=reg19*reg8; reg68=reg41+reg6; reg69=reg40+reg36; reg70=reg38-reg35;
    reg71=reg63*reg116; reg72=reg4-reg32; reg73=reg6-reg36; reg74=reg41-reg40; reg49=reg43+reg49;
    reg9=reg85-reg9; reg43=reg71+reg67; reg75=reg83*reg49; reg76=0.5*reg70; reg77=0.5*reg73;
    reg78=0.5*reg62; reg79=0.5*reg72; reg80=0.5*reg74; reg81=reg65+reg66; reg82=reg83*reg9;
    reg84=0.5*reg64; reg85=0.5*reg69; reg86=0.5*reg68; reg87=reg61+reg60; reg88=reg86*reg13;
    reg89=reg62*reg22; reg90=reg13*reg80; reg91=reg22*reg72; reg92=reg78*reg13; reg95=reg68*reg16;
    reg97=reg13*reg79; reg99=reg16*reg74; reg100=reg22*reg64; reg101=reg13*reg85; reg103=reg73*reg16;
    reg104=reg76*reg13; reg105=reg70*reg22; reg106=reg77*reg13; reg107=2*reg87; reg108=reg43+reg75;
    reg109=reg82+reg81; reg110=reg16*reg69; reg111=reg13*reg84; reg112=reg23*var_inter[0]; reg113=reg17*var_inter[1];
    reg114=reg9*reg85; reg115=reg18*reg64; reg117=reg76*reg9; reg118=reg73*reg29; reg119=reg68*reg29;
    reg120=reg77*reg9; reg121=reg70*reg18; reg122=reg68*reg116; reg123=reg78*reg49; reg124=reg78*reg9;
    reg125=reg29*reg69; reg126=reg9*reg84; reg92=reg92-reg95; reg127=reg62*reg8; reg128=reg86*reg49;
    reg130=reg18*reg72; T reg131=reg9*reg80; T reg132=reg29*reg74; T reg133=reg9*reg79; reg89=reg89-reg88;
    T reg134=reg108*reg68; T reg135=reg107*reg78; T reg136=reg109*reg64; T reg137=reg107*reg85; T reg138=reg107*reg84;
    T reg139=reg108*reg69; reg110=reg110-reg111; T reg140=reg8*reg64; T reg141=reg49*reg85; T reg142=reg49*reg84;
    T reg143=reg116*reg69; reg90=reg91+reg90; reg91=reg49*reg80; T reg144=reg8*reg72; reg97=reg99+reg97;
    reg99=reg49*reg79; T reg145=reg116*reg74; reg101=reg101-reg100; T reg146=reg73*reg116; T reg147=reg76*reg49;
    reg106=reg105+reg106; reg105=reg77*reg49; T reg148=reg70*reg8; reg104=reg103+reg104; reg103=reg86*reg9;
    T reg149=reg109*reg62; T reg150=reg62*reg18; T reg151=reg107*reg86; reg90=2*reg90; reg89=2*reg89;
    reg132=reg133+reg132; reg130=reg131+reg130; reg143=reg143-reg142; reg125=reg125-reg126; reg141=reg141-reg140;
    reg110=2*reg110; reg131=reg151-reg149; reg133=reg108*reg74; T reg152=reg107*reg79; T reg153=reg109*reg72;
    T reg154=reg107*reg80; T reg155=reg138-reg139; reg114=reg114-reg115; T reg156=var_inter[0]*var_inter[1]; reg104=2*reg104;
    T reg157=reg23*reg17; reg150=reg150-reg103; reg105=reg148+reg105; reg106=2*reg106; reg147=reg146+reg147;
    reg101=2*reg101; reg146=reg113*elem.f_vol_e[1]; reg148=reg113*elem.f_vol_e[0]; reg99=reg145+reg99; reg97=2*reg97;
    reg123=reg123-reg122; reg92=2*reg92; reg145=reg112*elem.f_vol_e[1]; reg91=reg144+reg91; reg127=reg127-reg128;
    reg144=reg112*elem.f_vol_e[0]; T reg158=reg134-reg135; T reg159=reg107*reg77; T reg160=reg109*reg70; reg124=reg124-reg119;
    T reg161=reg107*reg76; reg121=reg120+reg121; reg120=reg108*reg73; T reg162=reg136-reg137; reg118=reg117+reg118;
    reg117=reg110*reg86; T reg163=reg90*reg86; T reg164=reg130*reg62; T reg165=reg97*reg86; T reg166=reg132*reg62;
    T reg167=reg127*reg68; T reg168=reg89*reg78; T reg169=reg123*reg68; T reg170=reg92*reg78; T reg171=reg105*reg68;
    T reg172=reg106*reg78; T reg173=reg147*reg68; T reg174=reg104*reg78; T reg175=reg141*reg68; T reg176=reg101*reg78;
    T reg177=reg143*reg68; T reg178=reg110*reg78; T reg179=reg91*reg68; T reg180=reg90*reg78; T reg181=reg99*reg68;
    T reg182=reg97*reg78; T reg183=reg89*reg77; T reg184=reg150*reg70; T reg185=reg92*reg77; T reg186=reg124*reg70;
    T reg187=reg105*reg69; T reg188=reg99*reg69; T reg189=reg97*reg84; T reg190=reg91*reg69; T reg191=reg90*reg84;
    T reg192=reg104*reg84; T reg193=reg147*reg69; T reg194=reg101*reg84; T reg195=reg141*reg69; T reg196=reg143*reg69;
    T reg197=reg110*reg84; T reg198=reg106*reg76; T reg199=reg125*reg62; reg131=reg131-reg146; T reg200=reg110*reg79;
    T reg201=reg89*reg86; T reg202=reg150*reg62; T reg203=reg92*reg86; T reg204=reg124*reg62; T reg205=reg106*reg86;
    T reg206=reg121*reg62; T reg207=reg104*reg86; T reg208=reg118*reg62; T reg209=reg101*reg86; T reg210=reg114*reg62;
    T reg211=reg132*reg72; T reg212=reg89*reg79; T reg213=reg90*reg79; T reg214=reg127*reg74; T reg215=reg120+reg161;
    T reg216=reg91*reg74; T reg217=reg92*reg79; T reg218=reg124*reg72; T reg219=reg156*elem.f_vol_e[0]; T reg220=reg92*reg80;
    T reg221=reg97*reg79; T reg222=reg123*reg74; T reg223=reg156*elem.f_vol_e[1]; T reg224=reg99*reg74; T reg225=reg141*reg74;
    T reg226=reg101*reg79; T reg227=reg89*reg80; T reg228=reg157*elem.f_vol_e[0]; T reg229=reg147*reg74; T reg230=reg106*reg79;
    T reg231=reg150*reg72; T reg232=reg105*reg74; T reg233=reg160+reg159; T reg234=reg104*reg79; reg158=reg158-reg148;
    T reg235=reg106*reg77; T reg236=reg121*reg70; T reg237=reg104*reg77; T reg238=reg118*reg70; T reg239=reg101*reg77;
    T reg240=reg114*reg70; T reg241=reg110*reg77; T reg242=reg125*reg70; T reg243=reg90*reg77; T reg244=reg130*reg70;
    T reg245=reg97*reg77; T reg246=reg132*reg70; T reg247=reg89*reg76; T reg248=reg127*reg73; T reg249=reg92*reg76;
    T reg250=reg123*reg73; T reg251=reg143*reg74; T reg252=reg110*reg80; T reg253=reg125*reg72; T reg254=reg121*reg72;
    T reg255=reg110*reg85; T reg256=reg90*reg80; T reg257=reg130*reg72; T reg258=reg106*reg80; T reg259=reg157*elem.f_vol_e[1];
    T reg260=reg97*reg80; reg110=reg110*reg76; reg143=reg143*reg73; T reg261=reg90*reg76; reg127=reg127*reg69;
    T reg262=reg89*reg84; reg91=reg91*reg73; T reg263=reg97*reg76; reg99=reg99*reg73; reg97=reg97*reg85;
    reg150=reg150*reg64; reg132=reg132*reg64; reg89=reg89*reg85; reg124=reg124*reg64; T reg264=reg92*reg85;
    reg90=reg90*reg85; reg121=reg121*reg64; T reg265=reg106*reg85; T reg266=reg118*reg64; T reg267=reg104*reg85;
    T reg268=reg114*reg64; T reg269=reg101*reg85; T reg270=reg153+reg154; reg125=reg125*reg64; reg130=reg130*reg64;
    reg162=reg162-reg145; T reg271=reg101*reg80; T reg272=reg104*reg80; reg106=reg106*reg84; reg114=reg114*reg72;
    reg118=reg118*reg72; reg155=reg155-reg144; reg105=reg105*reg73; reg104=reg104*reg76; T reg273=reg133+reg152;
    reg147=reg147*reg73; reg123=reg123*reg69; reg92=reg92*reg84; reg141=reg141*reg73; reg101=reg101*reg76;
    reg239=reg240+reg239; reg237=reg238+reg237; reg266=reg267-reg266; reg230=reg232+reg230; reg235=reg236+reg235;
    reg185=reg186+reg185; reg268=reg269-reg268; reg183=reg184+reg183; reg213=reg216+reg213; reg155=reg7*reg155;
    reg184=reg233+reg223; reg181=reg182-reg181; reg179=reg180-reg179; reg125=reg255-reg125; reg234=reg229+reg234;
    reg177=reg178-reg177; reg175=reg176-reg175; reg130=reg90-reg130; reg212=reg214+reg212; reg173=reg174-reg173;
    reg158=reg7*reg158; reg110=reg143+reg110; reg263=reg99+reg263; reg256=reg257+reg256; reg260=reg211+reg260;
    reg101=reg141+reg101; reg252=reg253+reg252; reg150=reg89-reg150; reg220=reg218+reg220; reg271=reg114+reg271;
    reg131=reg7*reg131; reg221=reg224+reg221; reg89=reg215+reg219; reg217=reg222+reg217; reg249=reg250+reg249;
    reg247=reg248+reg247; reg124=reg264-reg124; reg261=reg91+reg261; reg245=reg246+reg245; reg104=reg147+reg104;
    reg243=reg244+reg243; reg121=reg265-reg121; reg226=reg225+reg226; reg241=reg242+reg241; reg227=reg231+reg227;
    reg258=reg254+reg258; reg198=reg105+reg198; reg163=reg164-reg163; reg197=reg196-reg197; reg92=reg123-reg92;
    reg251=reg200+reg251; reg191=reg190-reg191; reg262=reg127-reg262; reg117=reg199-reg117; reg194=reg195-reg194;
    reg165=reg166-reg165; reg272=reg118+reg272; reg90=reg228+reg273; reg192=reg193-reg192; reg209=reg210-reg209;
    reg162=reg7*reg162; reg167=reg168-reg167; reg171=reg172-reg171; reg106=reg187-reg106; reg201=reg202-reg201;
    reg207=reg208-reg207; reg169=reg170-reg169; reg189=reg188-reg189; reg91=reg270+reg259; reg203=reg204-reg203;
    reg205=reg206-reg205; reg132=reg97-reg132; reg251=reg251*reg7; reg213=reg7*reg213; reg198=reg198*reg7;
    reg271=reg271*reg7; reg131=ponderation*reg131; reg260=reg260*reg7; reg150=reg150*reg7; reg263=reg263*reg7;
    reg247=reg247*reg7; reg205=reg205*reg7; reg207=reg207*reg7; reg262=reg262*reg7; reg124=reg124*reg7;
    reg249=reg249*reg7; reg203=reg203*reg7; reg256=reg256*reg7; reg155=ponderation*reg155; reg162=ponderation*reg162;
    reg106=reg106*reg7; reg234=reg234*reg7; reg97=reg7*reg184; reg201=reg7*reg201; reg189=reg189*reg7;
    reg158=ponderation*reg158; reg252=reg252*reg7; reg230=reg230*reg7; reg104=reg104*reg7; reg226=reg226*reg7;
    reg192=reg192*reg7; reg227=reg227*reg7; reg99=reg7*reg89; reg194=reg194*reg7; reg220=reg220*reg7;
    reg101=reg101*reg7; reg221=reg7*reg221; reg197=reg197*reg7; reg110=reg110*reg7; reg191=reg191*reg7;
    reg217=reg217*reg7; reg258=reg258*reg7; reg261=reg261*reg7; reg92=reg92*reg7; reg212=reg212*reg7;
    reg177=reg177*reg7; reg117=reg117*reg7; reg239=reg239*reg7; reg105=reg7*reg91; reg266=reg266*reg7;
    reg179=reg179*reg7; reg167=reg167*reg7; reg114=reg7*reg90; reg237=reg237*reg7; reg125=reg125*reg7;
    reg163=reg163*reg7; reg235=reg235*reg7; reg181=reg181*reg7; reg272=reg272*reg7; reg185=reg185*reg7;
    reg183=reg183*reg7; reg268=reg268*reg7; reg165=reg165*reg7; reg171=reg171*reg7; reg132=reg132*reg7;
    reg245=reg245*reg7; reg173=reg173*reg7; reg169=reg169*reg7; reg209=reg209*reg7; reg243=reg243*reg7;
    reg175=reg175*reg7; reg241=reg241*reg7; reg130=reg130*reg7; reg121=reg121*reg7; T tmp_6_0=ponderation*reg181;
    T tmp_1_7=ponderation*reg227; T tmp_2_3=ponderation*reg194; T tmp_1_6=ponderation*reg220; reg118=ponderation*reg97; sollicitation[indices[2]+1]+=reg118;
    T tmp_4_3=ponderation*reg101; T tmp_6_7=ponderation*reg167; T tmp_6_4=ponderation*reg173; T tmp_3_1=ponderation*reg130; T tmp_5_7=ponderation*reg183;
    T tmp_6_5=ponderation*reg171; T tmp_2_0=ponderation*reg189; T tmp_2_5=ponderation*reg106; T tmp_2_2=ponderation*reg197; T tmp_0_5=ponderation*reg230;
    T tmp_1_4=ponderation*reg272; T tmp_6_2=ponderation*reg177; reg101=ponderation*reg105; sollicitation[indices[0]+1]+=reg101; T tmp_3_2=ponderation*reg125;
    reg106=ponderation*reg99; sollicitation[indices[2]+0]+=reg106; T tmp_2_4=ponderation*reg192; T tmp_0_4=ponderation*reg234; T tmp_6_3=ponderation*reg175;
    T tmp_1_2=ponderation*reg252; T tmp_0_3=ponderation*reg226; T tmp_6_1=ponderation*reg179; T tmp_6_6=ponderation*reg169; sollicitation[indices[3]+0]+=-reg158;
    T tmp_4_4=ponderation*reg104; T tmp_7_7=ponderation*reg201; T tmp_3_4=ponderation*reg266; T tmp_4_0=ponderation*reg263; T tmp_5_3=ponderation*reg239;
    T tmp_7_2=ponderation*reg117; T tmp_1_1=ponderation*reg256; T tmp_7_6=ponderation*reg203; sollicitation[indices[3]+1]+=-reg131; T tmp_5_2=ponderation*reg241;
    T tmp_3_7=ponderation*reg150; T tmp_5_1=ponderation*reg243; T tmp_2_7=ponderation*reg262; T tmp_7_5=ponderation*reg205; T tmp_0_2=ponderation*reg251;
    T tmp_3_5=ponderation*reg121; T tmp_7_3=ponderation*reg209; T tmp_7_4=ponderation*reg207; T tmp_5_0=ponderation*reg245; T tmp_4_6=ponderation*reg249;
    T tmp_3_6=ponderation*reg124; T tmp_3_0=ponderation*reg132; T tmp_4_7=ponderation*reg247; T tmp_7_0=ponderation*reg165; T tmp_0_0=ponderation*reg221;
    T tmp_4_2=ponderation*reg110; T tmp_5_6=ponderation*reg185; T tmp_0_6=ponderation*reg217; sollicitation[indices[1]+0]+=-reg155; T tmp_2_1=ponderation*reg191;
    T tmp_3_3=ponderation*reg268; reg104=ponderation*reg114; sollicitation[indices[0]+0]+=reg104; T tmp_1_3=ponderation*reg271; T tmp_1_5=ponderation*reg258;
    T tmp_4_1=ponderation*reg261; T tmp_5_5=ponderation*reg235; T tmp_7_1=ponderation*reg163; T tmp_0_7=ponderation*reg212; T tmp_5_4=ponderation*reg237;
    T tmp_4_5=ponderation*reg198; T tmp_2_6=ponderation*reg92; sollicitation[indices[1]+1]+=-reg162; T tmp_1_0=ponderation*reg260; T tmp_0_1=ponderation*reg213;
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
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=2*(*f.m).shear_modulus_13; T reg1=2*(*f.m).shear_modulus_23; T reg2=2*(*f.m).shear_modulus_12; reg0=1.0/reg0; reg1=1.0/reg1;
    T reg3=reg0*reg1; reg2=1.0/reg2; T reg4=pow((*f.m).v2[1],2); T reg5=pow((*f.m).v2[0],2); T reg6=pow((*f.m).v1[1],2);
    T reg7=pow((*f.m).v1[0],2); T reg8=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg11=reg2*reg3;
    reg6=reg7+reg6; reg7=pow((*f.m).v1[2],2); T reg12=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg13=reg9*reg11; T reg14=1.0/(*f.m).elastic_modulus_2;
    T reg15=reg10*reg11; reg4=reg5+reg4; reg5=pow((*f.m).v2[2],2); T reg16=reg8*reg11; T reg17=reg10*reg15;
    T reg18=reg14*reg13; reg7=reg6+reg7; reg5=reg4+reg5; reg4=reg10*reg16; reg6=reg12*reg13;
    reg5=pow(reg5,0.5); reg4=reg6+reg4; T reg19=reg12*reg15; T reg20=reg14*reg16; T reg21=1.0/(*f.m).elastic_modulus_1;
    reg17=reg18-reg17; reg7=pow(reg7,0.5); reg18=reg12*reg4; T reg22=(*f.m).v2[1]/reg5; T reg23=(*f.m).v2[2]/reg5;
    T reg24=reg21*reg17; T reg25=(*f.m).v1[2]/reg7; T reg26=reg19+reg20; T reg27=(*f.m).v1[1]/reg7; reg13=reg21*reg13;
    T reg28=reg8*reg16; T reg29=reg2*reg1; T reg30=reg27*reg23; T reg31=reg25*reg22; T reg32=reg8*reg3;
    T reg33=reg9*reg3; T reg34=reg12*reg11; reg3=reg10*reg3; reg5=(*f.m).v2[0]/reg5; reg7=(*f.m).v1[0]/reg7;
    T reg35=reg8*reg15; reg11=reg14*reg11; T reg36=reg8*reg26; reg18=reg24-reg18; reg24=reg10*reg32;
    T reg37=reg10*reg29; reg16=reg12*reg16; T reg38=2*reg7; T reg39=reg10*reg3; T reg40=reg12*reg33;
    reg36=reg18-reg36; reg35=reg6+reg35; reg6=reg2*reg0; reg18=reg30-reg31; reg33=reg14*reg33;
    T reg41=reg8*reg29; T reg42=reg8*reg11; T reg43=reg8*reg34; T reg44=reg25*reg5; T reg45=reg7*reg23;
    T reg46=2*reg5; reg15=reg21*reg15; reg28=reg13-reg28; reg29=reg9*reg29; reg17=reg17/reg36;
    reg13=pow(reg22,2); T reg47=reg8*reg6; T reg48=pow(reg5,2); T reg49=reg9*reg6; T reg50=pow(reg7,2);
    T reg51=pow(reg27,2); T reg52=reg10*reg41; T reg53=reg10*reg37; T reg54=reg12*reg29; reg29=reg14*reg29;
    reg39=reg33-reg39; reg11=reg21*reg11; reg24=reg40+reg24; reg33=reg27*reg38; reg40=reg46*reg22;
    T reg55=reg27*reg5; T reg56=reg7*reg22; T reg57=2*reg18; reg28=reg28/reg36; reg4=reg4/reg36;
    reg43=reg15+reg43; reg16=reg15+reg16; reg34=reg12*reg34; reg42=reg19+reg42; reg32=reg14*reg32;
    reg3=reg12*reg3; reg6=reg10*reg6; reg35=reg35/reg36; reg15=reg44-reg45; T reg58=pow(reg23,2);
    T reg59=reg14*reg49; reg49=reg12*reg49; T reg60=reg10*reg6; T reg61=reg10*reg47; T reg62=reg51*reg17;
    T reg63=reg13*reg4; T reg64=reg33*reg17; T reg65=reg40*reg4; T reg66=reg50*reg35; T reg67=reg48*reg28;
    T reg68=reg51*reg35; T reg69=reg13*reg28; T reg70=reg33*reg35; T reg71=reg40*reg28; reg37=reg12*reg37;
    T reg72=pow(reg18,2); reg52=reg54+reg52; reg41=reg14*reg41; reg53=reg29-reg53; reg29=pow(reg15,2);
    reg54=reg57*reg15; reg26=reg26/reg36; T reg73=reg3+reg32; T reg74=reg50*reg17; T reg75=reg48*reg4;
    reg24=reg12*reg24; reg43=reg43/reg36; reg16=reg16/reg36; reg42=reg42/reg36; reg39=reg21*reg39;
    reg34=reg11-reg34; reg11=pow(reg25,2); T reg76=reg56-reg55; T reg77=reg13*reg43; T reg78=reg51*reg42;
    T reg79=reg48*reg43; T reg80=reg50*reg42; T reg81=reg72*reg26; T reg82=reg54*reg16; reg71=reg70+reg71;
    reg70=reg58*reg28; T reg83=reg11*reg35; T reg84=reg29*reg16; reg69=reg68+reg69; reg68=reg72*reg16;
    reg67=reg66+reg67; reg34=reg34/reg36; reg75=reg74+reg75; reg66=reg54*reg26; reg65=reg64+reg65;
    reg64=reg58*reg4; reg74=reg11*reg17; T reg85=reg29*reg26; reg63=reg62+reg63; reg61=reg49+reg61;
    reg6=reg12*reg6; reg60=reg59-reg60; reg49=reg37+reg41; reg52=reg12*reg52; reg53=reg21*reg53;
    reg73=reg8*reg73; reg24=reg39-reg24; reg47=reg14*reg47; reg39=reg33*reg42; reg59=pow(reg76,2);
    reg62=reg40*reg43; T reg86=reg59*reg16; reg70=reg83+reg70; reg62=reg39+reg62; reg54=reg54*reg34;
    reg84=reg69+reg84; reg39=reg7*reg27; reg60=reg21*reg60; reg68=reg67+reg68; reg67=reg46*reg23;
    reg61=reg12*reg61; reg69=2*reg22; reg83=reg6+reg47; reg81=reg75+reg81; reg85=reg63+reg85;
    reg64=reg74+reg64; reg66=reg65+reg66; reg63=reg59*reg26; reg65=reg58*reg43; reg79=reg80+reg79;
    reg74=reg72*reg34; reg75=reg25*reg38; reg80=2*reg27; reg73=reg24-reg73; reg52=reg53-reg52;
    reg24=reg11*reg42; reg53=reg5*reg22; reg77=reg78+reg77; reg78=reg29*reg34; reg49=reg8*reg49;
    reg82=reg71+reg82; reg71=reg13*reg84; T reg87=reg13*reg68; reg57=reg57*reg76; T reg88=reg51*reg85;
    T reg89=2*reg15; reg63=reg64+reg63; reg83=reg8*reg83; reg64=reg7*reg5; reg73=reg73/reg36;
    reg61=reg60-reg61; reg60=reg7*reg15; T reg90=reg27*reg22; T reg91=reg39*reg66; T reg92=reg27*reg18;
    T reg93=reg53*reg82; T reg94=reg48*reg84; T reg95=reg50*reg85; T reg96=reg25*reg80; reg49=reg52-reg49;
    reg56=reg55+reg56; reg78=reg77+reg78; reg74=reg79+reg74; reg65=reg24+reg65; reg24=reg51*reg66;
    reg52=reg15*reg18; reg55=reg48*reg82; reg66=reg50*reg66; reg77=reg59*reg34; reg79=reg67*reg28;
    T reg97=reg75*reg35; reg82=reg13*reg82; T reg98=2*reg25; reg86=reg70+reg86; reg54=reg62+reg54;
    reg62=reg22*reg80; reg70=reg48*reg68; T reg99=reg50*reg81; T reg100=reg51*reg81; T reg101=reg75*reg17;
    T reg102=reg67*reg4; reg85=reg39*reg85; T reg103=reg5*reg38; T reg104=reg69*reg23; reg84=reg53*reg84;
    reg87=reg100+reg87; reg100=reg29*reg74; T reg105=reg62*reg12; T reg106=reg103*reg21; T reg107=reg13*reg12;
    reg55=reg66+reg55; reg66=reg48*reg21; T reg108=reg72*reg54; reg71=reg88+reg71; reg88=reg13*reg14;
    reg60=reg92+reg60; reg81=reg39*reg81; reg92=reg48*reg12; T reg109=reg50*reg63; T reg110=reg62*reg14;
    T reg111=reg25*reg23; T reg112=reg48*reg86; reg68=reg53*reg68; T reg113=reg22*reg18; reg93=reg91+reg93;
    reg91=reg52*reg78; reg84=reg85+reg84; reg89=reg89*reg76; reg98=reg23*reg98; reg85=reg5*reg15;
    T reg114=reg52*reg54; T reg115=reg7*reg18; T reg116=reg27*reg15; T reg117=reg72*reg78; reg94=reg95+reg94;
    reg49=reg49/reg36; reg83=reg61-reg83; reg61=reg72*reg74; reg70=reg99+reg70; reg102=reg101+reg102;
    reg95=reg57*reg26; reg17=reg96*reg17; reg4=reg104*reg4; reg99=reg103*reg12; reg54=reg29*reg54;
    reg82=reg24+reg82; reg79=reg97+reg79; reg24=reg57*reg16; reg35=reg96*reg35; reg28=reg104*reg28;
    reg77=reg65+reg77; reg65=reg75*reg42; reg97=reg67*reg43; reg101=reg64*reg73; T reg118=reg90*reg73;
    T reg119=reg56*reg73; T reg120=reg51*reg63; reg78=reg29*reg78; T reg121=reg13*reg86; reg63=reg39*reg63;
    T reg122=reg58*reg10; reg16=reg89*reg16; T reg123=reg103*reg119; reg108=reg55+reg108; reg28=reg35+reg28;
    reg35=reg48*reg8; reg24=reg79+reg24; reg55=reg98*reg10; reg79=reg29*reg77; T reg124=reg5*reg18;
    T reg125=reg22*reg15; reg80=reg15*reg80; reg97=reg65+reg97; reg57=reg57*reg34; reg78=reg71+reg78;
    reg99=reg110-reg99; reg42=reg96*reg42; reg43=reg104*reg43; reg65=reg62*reg118; reg92=reg88-reg92;
    reg71=reg111*reg73; reg86=reg53*reg86; reg115=reg115*reg49; reg116=reg116*reg49; reg88=reg25*reg76;
    reg110=reg13*reg10; T reg126=reg60*reg49; reg74=reg52*reg74; reg95=reg102+reg95; reg85=reg113+reg85;
    reg121=reg120+reg121; reg102=reg103*reg118; reg113=reg98*reg8; reg120=reg62*reg101; reg91=reg84+reg91;
    reg118=reg56*reg118; reg100=reg87+reg100; reg61=reg70+reg61; reg36=reg83/reg36; reg70=reg103*reg101;
    reg83=reg103*reg8; reg84=reg56*reg119; reg114=reg93+reg114; reg105=reg106-reg105; reg117=reg94+reg117;
    reg87=reg58*reg8; reg54=reg82+reg54; reg119=reg62*reg119; reg112=reg109+reg112; reg26=reg89*reg26;
    reg4=reg17+reg4; reg38=reg18*reg38; reg107=reg66-reg107; reg17=reg72*reg77; reg68=reg81+reg68;
    reg66=reg62*reg10; reg81=reg51*reg95; reg82=reg23*reg76; reg93=reg80*reg126; reg94=reg58*reg9;
    reg106=reg25*reg18; reg110=reg35+reg110; reg35=reg80*reg115; reg66=reg83+reg66; reg83=reg80*reg116;
    reg65=reg78+reg65; reg119=reg54+reg119; reg54=reg13*reg24; reg78=reg7*reg76; reg120=reg100+reg120;
    reg98=reg98*reg9; reg14=reg51*reg14; reg100=(*f.m).alpha_1*reg50; reg86=reg63+reg86; reg77=reg52*reg77;
    reg63=reg51*reg12; reg109=reg38*reg116; reg102=reg117+reg102; reg117=reg38*reg115; reg70=reg61+reg70;
    reg26=reg4+reg26; reg21=reg50*reg21; reg4=reg85*reg36; reg61=(*f.m).alpha_2*reg48; reg125=reg125*reg36;
    reg124=reg124*reg36; reg16=reg28+reg16; reg57=reg97+reg57; reg43=reg42+reg43; reg34=reg89*reg34;
    reg122=reg92-reg122; reg28=(*f.m).alpha_1*reg51; reg88=reg88*reg49; reg44=reg45+reg44; reg55=reg99-reg55;
    reg42=reg48*reg24; reg45=reg50*reg95; reg89=reg38*reg126; reg123=reg108+reg123; reg92=reg103*reg71;
    reg17=reg112+reg17; reg69=reg69*reg15; reg46=reg46*reg18; reg74=reg68+reg74; reg101=reg56*reg101;
    reg118=reg91+reg118; reg116=reg60*reg116; reg12=reg50*reg12; reg113=reg105-reg113; reg84=reg114+reg84;
    reg126=reg60*reg126; reg87=reg107-reg87; reg79=reg121+reg79; reg68=reg62*reg71; reg91=(*f.m).alpha_2*reg13;
    reg110=reg94-reg110; reg94=reg11*reg8; reg63=reg21-reg63; reg66=reg98-reg66; reg12=reg14-reg12;
    reg14=reg90*reg55; reg21=reg11*reg10; reg97=reg64*reg113; reg98=(*f.m).alpha_2*reg58; reg99=(*f.m).alpha_1*reg11;
    reg105=(*f.m).alpha_3*reg29; reg91=reg28+reg91; reg28=reg90*reg122; reg107=reg64*reg87; reg108=(*f.m).alpha_3*reg72;
    reg61=reg100+reg61; reg100=reg87*reg50; reg95=reg39*reg95; reg24=reg53*reg24; reg112=reg51*reg26;
    reg114=reg13*reg16; reg121=reg69*reg124; reg35=reg120+reg35; reg34=reg43+reg34; reg43=reg44*reg73;
    reg8=reg50*reg8; reg82=reg82*reg36; reg117=reg70+reg117; reg70=reg46*reg124; reg120=reg48*reg16;
    reg109=reg102+reg109; reg102=reg46*reg125; T reg127=reg50*reg26; reg71=reg56*reg71; reg77=reg86+reg77;
    reg86=reg72*reg57; reg42=reg45+reg42; reg45=reg85*reg125; reg116=reg118+reg116; reg118=reg46*reg4;
    reg89=reg123+reg89; reg115=reg60*reg115; reg101=reg74+reg101; reg92=reg17+reg92; reg17=reg38*reg88;
    reg10=reg51*reg10; reg74=reg5*reg76; reg123=reg23*reg18; T reg128=reg80*reg88; T reg129=reg27*reg76;
    T reg130=reg25*reg15; reg78=reg106+reg78; reg30=reg31+reg30; reg68=reg79+reg68; reg31=reg55*reg13;
    reg79=reg113*reg48; reg125=reg69*reg125; reg83=reg65+reg83; reg65=reg122*reg13; reg87=reg87*reg48;
    reg93=reg119+reg93; reg122=reg122*reg51; reg106=reg29*reg57; reg54=reg81+reg54; reg81=reg85*reg4;
    reg126=reg84+reg126; reg4=reg69*reg4; reg113=reg113*reg50; reg55=reg55*reg51; reg17=reg92+reg17;
    reg84=reg46*reg82; reg92=reg69*reg82; reg128=reg68+reg128; reg68=reg53*reg2; reg9=reg11*reg9;
    reg118=reg89+reg118; reg10=reg8+reg10; reg8=reg56*reg2; reg106=reg54+reg106; reg72=reg72*reg34;
    reg86=reg42+reg86; reg42=reg62*reg43; reg121=reg35+reg121; reg35=reg103*reg43; reg4=reg93+reg4;
    reg125=reg83+reg125; reg120=reg127+reg120; reg114=reg112+reg114; reg29=reg29*reg34; reg98=reg99+reg98;
    reg59=(*f.m).alpha_3*reg59; reg54=(*f.m).alpha_1*reg39; reg83=(*f.m).alpha_2*reg53; reg94=reg63-reg94; reg122=reg100+reg122;
    reg63=reg110*reg11; reg81=reg126+reg81; reg55=reg113+reg55; reg89=reg66*reg11; reg65=reg87+reg65;
    reg87=reg58*reg110; reg31=reg79+reg31; reg79=reg58*reg66; reg7=reg7*reg25; reg5=reg5*reg23;
    reg129=reg130+reg129; reg74=reg123+reg74; reg93=reg23*reg15; reg99=reg22*reg76; reg66=reg111*reg66;
    reg14=reg97+reg14; reg115=reg101+reg115; reg124=reg85*reg124; reg102=reg109+reg102; reg45=reg116+reg45;
    reg110=reg111*reg110; reg28=reg107+reg28; reg71=reg77+reg71; reg88=reg60*reg88; reg70=reg117+reg70;
    reg78=reg78*reg49; reg73=reg30*reg73; reg108=reg61+reg108; reg16=reg53*reg16; reg26=reg39*reg26;
    reg21=reg12-reg21; reg57=reg52*reg57; reg105=reg91+reg105; reg24=reg95+reg24; reg92=reg128+reg92;
    reg103=reg103*reg73; reg72=reg120+reg72; reg12=reg38*reg78; reg35=reg86+reg35; reg84=reg17+reg84;
    reg74=reg74*reg36; reg49=reg129*reg49; reg34=reg52*reg34; reg16=reg26+reg16; reg43=reg56*reg43;
    reg57=reg24+reg57; reg23=reg22*reg23; reg25=reg27*reg25; reg18=reg76*reg18; reg17=reg105*reg125;
    reg22=reg108*reg121; reg24=reg105*reg102; reg26=reg108*reg70; reg27=reg45*reg118; reg53=reg4*reg45;
    reg61=reg81*reg102; reg77=reg81*reg125; reg86=(*f.m).alpha_2*reg5; reg91=(*f.m).alpha_1*reg7; reg52=(*f.m).alpha_3*reg52;
    reg83=reg54+reg83; reg59=reg98+reg59; reg89=reg55+reg89; reg88=reg71+reg88; reg82=reg85*reg82;
    reg54=reg68*reg33; reg42=reg106+reg42; reg55=reg80*reg78; reg63=reg122+reg63; reg5=reg5*reg0;
    reg71=reg44*reg0; reg29=reg114+reg29; reg62=reg62*reg73; reg95=reg21*reg51; reg110=reg28+reg110;
    reg28=reg56*reg8; reg66=reg14+reg66; reg14=reg94*reg50; reg97=reg56*reg68; reg98=reg8*reg33;
    reg124=reg115+reg124; reg99=reg93+reg99; reg10=reg9-reg10; reg9=reg94*reg48; reg8=reg40*reg8;
    reg93=reg21*reg13; reg79=reg31+reg79; reg68=reg40*reg68; reg87=reg65+reg87; reg31=reg4*reg102;
    reg52=reg83+reg52; reg86=reg91+reg86; reg27=reg61-reg27; reg53=reg77-reg53; reg61=(*f.m).alpha_2*reg23;
    reg65=(*f.m).alpha_1*reg25; reg77=reg30*reg1; reg23=reg23*reg1; reg18=(*f.m).alpha_3*reg18; reg38=reg38*reg49;
    reg36=reg99*reg36; reg82=reg88+reg82; reg103=reg72+reg103; reg72=reg46*reg74; reg12=reg35+reg12;
    reg93=reg9+reg93; reg58=reg58*reg10; reg68=reg87+reg68; reg9=reg5*reg67; reg8=reg79+reg8;
    reg35=reg71*reg67; reg79=reg71*reg75; reg98=reg89+reg98; reg83=reg5*reg75; reg54=reg63+reg54;
    reg43=reg57+reg43; reg78=reg60*reg78; reg15=reg76*reg15; reg11=reg10*reg11; reg95=reg14+reg95;
    reg34=reg16+reg34; reg73=reg56*reg73; reg14=reg45*reg105; reg16=reg124*reg108; reg92=reg59*reg92;
    reg17=reg22+reg17; reg84=reg59*reg84; reg24=reg26+reg24; reg22=reg118*reg125; reg80=reg80*reg49;
    reg5=reg5*reg44; reg97=reg110+reg97; reg62=reg29+reg62; reg94=reg64*reg94; reg21=reg90*reg21;
    reg55=reg42+reg55; reg28=reg66+reg28; reg71=reg71*reg44; reg2=reg39*reg2; reg26=reg69*reg74;
    reg11=reg95+reg11; reg10=reg111*reg10; reg0=reg7*reg0; reg74=reg85*reg74; reg78=reg43+reg78;
    reg26=reg55+reg26; reg7=reg2*reg33; reg21=reg94+reg21; reg69=reg69*reg36; reg83=reg54+reg83;
    reg29=reg23*reg96; reg22=reg31-reg22; reg31=reg27*reg121; reg39=reg53*reg70; reg5=reg97+reg5;
    reg42=reg23*reg30; reg84=reg24+reg84; reg24=reg52*reg118; reg92=reg17+reg92; reg17=reg4*reg52;
    reg71=reg28+reg71; reg28=reg77*reg30; reg14=reg16+reg14; reg59=reg82*reg59; reg15=(*f.m).alpha_3*reg15;
    reg61=reg65+reg61; reg18=reg86+reg18; reg49=reg60*reg49; reg73=reg34+reg73; reg80=reg62+reg80;
    reg35=reg8+reg35; reg38=reg103+reg38; reg8=reg77*reg104; reg77=reg77*reg96; reg72=reg12+reg72;
    reg79=reg98+reg79; reg23=reg23*reg104; reg58=reg93+reg58; reg9=reg68+reg9; reg12=reg40*reg2;
    reg46=reg46*reg36; reg7=reg11+reg7; reg11=reg81*reg121; reg77=reg79+reg77; reg8=reg35+reg8;
    reg15=reg61+reg15; reg75=reg0*reg75; reg74=reg78+reg74; reg10=reg21+reg10; reg1=reg25*reg1;
    reg2=reg56*reg2; reg31=reg39-reg31; reg16=reg124*reg22; reg23=reg9+reg23; reg28=reg71+reg28;
    reg49=reg73+reg49; reg36=reg85*reg36; reg26=reg26*reg18; reg9=reg124*reg118; reg17=reg92+reg17;
    reg67=reg0*reg67; reg72=reg18*reg72; reg12=reg58+reg12; reg29=reg83+reg29; reg46=reg38+reg46;
    reg69=reg80+reg69; reg42=reg5+reg42; reg5=reg81*reg70; reg24=reg84+reg24; reg59=reg14+reg59;
    reg14=reg52*reg81; reg21=reg4*reg124; reg72=reg24+reg72; reg36=reg49+reg36; reg46=reg15*reg46;
    reg26=reg17+reg26; reg69=reg69*reg15; reg44=reg0*reg44; reg67=reg12+reg67; reg16=reg31+reg16;
    reg0=1-var_inter[1]; reg104=reg1*reg104; reg14=reg59+reg14; reg74=reg18*reg74; reg12=reg118*reg121;
    reg17=reg4*reg70; reg18=reg124*reg102; reg9=reg5-reg9; reg5=reg45*reg70; reg24=reg124*reg125;
    reg96=reg1*reg96; reg75=reg7+reg75; reg21=reg11-reg21; reg7=reg45*reg121; reg11=1-var_inter[0];
    reg25=reg77*reg42; reg2=reg10+reg2; reg10=reg29*reg28; reg31=reg23*reg28; reg34=reg8*reg42;
    reg35=reg102*reg121; reg24=reg7-reg24; reg27=reg27/reg16; reg7=elem.pos(1)[0]*var_inter[0]; reg30=reg1*reg30;
    reg9=reg9/reg16; reg1=reg11*elem.pos(0)[1]; reg38=elem.pos(1)[1]*var_inter[0]; reg39=reg11*elem.pos(0)[0]; reg18=reg5-reg18;
    reg53=reg53/reg16; reg12=reg17-reg12; reg21=reg21/reg16; reg5=reg70*reg125; reg17=reg0*elem.pos(1)[1];
    reg43=reg0*elem.pos(0)[1]; reg96=reg75+reg96; reg49=reg0*elem.pos(1)[0]; reg54=reg0*elem.pos(0)[0]; reg34=reg31-reg34;
    reg44=reg2+reg44; reg36=reg15*reg36; reg25=reg10-reg25; reg46=reg72+reg46; reg69=reg26+reg69;
    reg104=reg67+reg104; reg2=reg77*reg23; reg74=reg14+reg74; reg10=reg29*reg8; reg21=reg21*reg46;
    reg27=reg27*reg69; reg12=reg12/reg16; reg24=reg24/reg16; reg49=reg49-reg54; reg14=elem.pos(2)[0]*var_inter[1];
    reg9=reg9*reg69; reg15=elem.pos(2)[1]*var_inter[1]; reg22=reg22/reg16; reg18=reg18/reg16; reg26=elem.pos(2)[1]*var_inter[0];
    reg35=reg5-reg35; reg2=reg10-reg2; reg5=elem.pos(2)[0]*var_inter[0]; reg10=reg39+reg7; reg30=reg44+reg30;
    reg31=reg38+reg1; reg36=reg74+reg36; reg44=reg104*reg25; reg55=reg96*reg34; reg53=reg53*reg46;
    reg17=reg17-reg43; reg21=reg9-reg21; reg46=reg24*reg46; reg69=reg18*reg69; reg16=reg35/reg16;
    reg12=reg12*reg36; reg27=reg53-reg27; reg22=reg22*reg36; reg9=reg23*reg30; reg18=reg104*reg42;
    reg42=reg96*reg42; reg24=reg30*reg2; reg44=reg55-reg44; reg35=reg29*reg30; reg53=elem.pos(3)[1]*reg11;
    reg26=reg26-reg31; reg55=elem.pos(3)[1]*var_inter[1]; reg17=reg15+reg17; reg15=elem.pos(3)[0]*var_inter[1]; reg5=reg5-reg10;
    reg14=reg49+reg14; reg49=elem.pos(3)[0]*reg11; reg36=reg16*reg36; reg69=reg46-reg69; reg24=reg44+reg24;
    reg53=reg26+reg53; reg49=reg5+reg49; reg5=reg104*reg28; reg16=reg8*reg30; reg14=reg14-reg15;
    reg27=reg22+reg27; reg28=reg96*reg28; reg22=1-(*f.m).resolution; reg12=reg21-reg12; reg23=reg96*reg23;
    reg29=reg29*reg104; reg35=reg42-reg35; reg9=reg18-reg9; reg17=reg17-reg55; reg30=reg77*reg30;
    reg18=reg17*reg49; reg12=reg12*reg22; reg104=reg77*reg104; reg21=reg14*reg53; reg8=reg96*reg8;
    reg69=reg36+reg69; reg27=reg27*reg22; reg26=reg105*(*f.m).resolution; reg35=reg35/reg24; reg30=reg28-reg30;
    reg9=reg9/reg24; reg28=reg108*(*f.m).resolution; reg16=reg5-reg16; reg29=reg23-reg29; reg27=reg28+reg27;
    reg69=reg69*reg22; reg18=reg21-reg18; reg104=reg8-reg104; reg29=reg29/reg24; reg124=reg124*reg22;
    reg45=reg45*reg22; reg5=(*f.m).resolution*reg35; reg8=(*f.m).resolution*reg9; reg21=reg52*(*f.m).resolution; reg12=reg26+reg12;
    reg25=reg25/reg24; reg16=reg16/reg24; reg30=reg30/reg24; reg34=reg34/reg24; reg12=reg12*(*f.m).deltaT;
    reg27=(*f.m).deltaT*reg27; reg17=reg17/reg18; reg23=(*f.m).resolution*reg25; reg69=reg21+reg69; reg49=reg49/reg18;
    reg21=(*f.m).resolution*reg34; reg70=reg22*reg70; reg102=reg22*reg102; reg121=reg22*reg121; reg125=reg22*reg125;
    reg81=reg22*reg81; reg2=reg2/reg24; reg26=(*f.m).resolution*reg16; reg28=(*f.m).resolution*reg30; reg36=(*f.m).resolution*reg29;
    reg24=reg104/reg24; reg53=reg53/reg18; reg8=reg124+reg8; reg5=reg45-reg5; reg14=reg14/reg18;
    reg42=reg27*reg8; reg44=reg11*reg14; reg69=reg69*(*f.m).deltaT; reg45=reg11*reg17; reg46=reg0*reg53;
    reg57=reg12*reg5; reg58=reg0*reg49; reg59=var_inter[0]*reg14; reg60=var_inter[0]*reg17; reg61=var_inter[1]*reg53;
    reg62=var_inter[1]*reg49; reg36=reg81+reg36; reg28=reg125+reg28; reg26=reg121-reg26; reg23=reg102-reg23;
    reg63=(*f.m).resolution*reg2; reg118=reg22*reg118; reg22=reg4*reg22; reg70=reg21+reg70; reg4=(*f.m).resolution*reg24;
    reg21=reg61-reg60; reg65=reg58-reg44; reg66=reg59-reg62; reg67=reg12*reg28; reg68=reg27*reg26;
    reg71=reg58+reg59; reg72=reg45+reg61; reg73=reg46+reg60; reg74=reg62+reg44; reg75=reg69*reg36;
    reg76=reg42+reg57; reg4=reg22-reg4; reg22=reg12*reg23; reg77=reg27*reg70; reg118=reg63+reg118;
    reg63=reg45-reg46; reg78=0.5*reg74; reg79=0.5*reg63; reg80=reg76+reg75; reg81=0.5*reg66;
    reg82=0.5*reg65; reg83=0.5*reg72; reg84=reg77+reg22; reg85=reg69*reg118; reg86=0.5*reg21;
    reg87=0.5*reg71; reg88=reg68+reg67; reg89=reg69*reg4; reg91=0.5*reg73; reg92=2*reg80;
    reg93=reg21*reg8; reg94=reg36*reg87; reg95=reg84+reg85; reg96=reg8*reg73; reg97=reg36*reg82;
    reg98=reg5*reg71; reg99=reg81*reg36; reg100=reg36*reg91; reg101=reg72*reg8; reg102=reg78*reg36;
    reg103=reg83*reg36; reg104=reg5*reg65; reg106=reg8*reg63; reg107=reg36*reg79; reg109=reg89+reg88;
    reg110=reg86*reg36; reg111=reg66*reg5; reg112=reg74*reg5; reg113=reg28*reg71; reg114=reg0*var_inter[0];
    reg115=reg4*reg91; reg116=reg81*reg4; reg117=reg21*reg26; reg119=reg66*reg28; reg120=reg86*reg4;
    reg121=reg78*reg4; reg122=reg72*reg70; reg123=reg78*reg118; reg124=reg11*var_inter[1]; reg125=reg26*reg73;
    reg126=reg4*reg87; reg127=reg28*reg65; reg128=reg4*reg79; reg112=reg112-reg103; reg129=reg72*reg26;
    reg102=reg102-reg101; reg130=reg74*reg23; T reg131=reg83*reg118; T reg132=reg83*reg4; T reg133=reg74*reg28;
    reg110=reg111+reg110; reg111=reg118*reg82; T reg134=reg86*reg118; T reg135=reg66*reg23; T reg136=reg70*reg63;
    T reg137=reg92*reg83; reg97=reg106+reg97; reg106=reg95*reg72; T reg138=reg92*reg78; reg99=reg93+reg99;
    reg93=reg81*reg118; T reg139=reg21*reg70; T reg140=reg23*reg65; T reg141=reg118*reg79; T reg142=reg92*reg87;
    T reg143=reg109*reg74; reg96=reg96-reg94; T reg144=reg118*reg87; T reg145=reg70*reg73; T reg146=reg95*reg73;
    T reg147=reg23*reg71; T reg148=reg118*reg91; T reg149=reg109*reg71; reg107=reg104+reg107; reg100=reg100-reg98;
    reg104=reg92*reg91; T reg150=reg142-reg146; reg112=2*reg112; reg133=reg133-reg132; reg130=reg130-reg131;
    T reg151=reg92*reg79; reg127=reg128+reg127; reg102=2*reg102; reg128=reg109*reg65; reg121=reg121-reg129;
    T reg152=reg137-reg143; reg125=reg125-reg126; T reg153=reg92*reg86; reg119=reg120+reg119; reg120=reg149-reg104;
    T reg154=reg109*reg66; T reg155=reg92*reg81; reg117=reg116+reg117; reg116=reg95*reg21; reg115=reg115-reg113;
    T reg156=reg106-reg138; T reg157=reg0*reg11; reg99=2*reg99; reg97=2*reg97; reg145=reg145-reg144;
    reg111=reg136+reg111; reg134=reg135+reg134; reg110=2*reg110; reg96=2*reg96; reg93=reg139+reg93;
    reg135=reg124*elem.f_vol_e[1]; reg141=reg140+reg141; reg136=reg124*elem.f_vol_e[0]; reg139=reg95*reg63; reg140=reg114*elem.f_vol_e[1];
    T reg158=reg114*elem.f_vol_e[0]; reg100=2*reg100; reg148=reg148-reg147; T reg159=reg92*reg82; T reg160=var_inter[0]*var_inter[1];
    reg123=reg123-reg122; reg107=2*reg107; T reg161=reg116+reg155; T reg162=reg112*reg83; T reg163=reg133*reg74;
    T reg164=reg130*reg72; T reg165=reg112*reg78; T reg166=reg123*reg72; T reg167=reg102*reg78; T reg168=reg99*reg91;
    T reg169=reg115*reg71; T reg170=reg100*reg91; T reg171=reg154+reg153; T reg172=reg112*reg87; T reg173=reg130*reg73;
    T reg174=reg148*reg73; T reg175=reg100*reg87; T reg176=reg102*reg87; T reg177=reg93*reg73; T reg178=reg99*reg87;
    T reg179=reg134*reg73; T reg180=reg123*reg73; T reg181=reg110*reg87; T reg182=reg145*reg63; T reg183=reg107*reg82;
    T reg184=reg141*reg63; T reg185=reg97*reg82; T reg186=reg111*reg63; T reg187=reg157*elem.f_vol_e[0]; T reg188=reg139+reg159;
    T reg189=reg123*reg63; T reg190=reg128+reg151; T reg191=reg102*reg82; T reg192=reg130*reg63; T reg193=reg112*reg82;
    reg150=reg150-reg158; T reg194=reg127*reg65; T reg195=reg107*reg79; T reg196=reg125*reg65; T reg197=reg96*reg79;
    reg120=reg120-reg140; reg123=reg123*reg21; T reg198=reg102*reg81; reg130=reg130*reg21; T reg199=reg112*reg81;
    T reg200=reg119*reg66; T reg201=reg110*reg86; T reg202=reg121*reg66; T reg203=reg102*reg86; T reg204=reg133*reg66;
    T reg205=reg112*reg86; T reg206=reg112*reg91; T reg207=reg133*reg71; T reg208=reg93*reg21; T reg209=reg99*reg81;
    T reg210=reg134*reg21; reg156=reg156-reg136; T reg211=reg115*reg65; T reg212=reg100*reg79; T reg213=reg117*reg65;
    T reg214=reg99*reg79; T reg215=reg119*reg65; T reg216=reg110*reg79; T reg217=reg121*reg65; T reg218=reg102*reg79;
    reg152=reg152-reg135; reg133=reg133*reg65; reg112=reg112*reg79; T reg219=reg145*reg73; T reg220=reg96*reg87;
    T reg221=reg157*elem.f_vol_e[1]; T reg222=reg160*elem.f_vol_e[0]; T reg223=reg160*elem.f_vol_e[1]; T reg224=reg110*reg82; T reg225=reg134*reg63;
    T reg226=reg99*reg82; T reg227=reg93*reg63; T reg228=reg100*reg82; T reg229=reg148*reg63; T reg230=reg96*reg82;
    T reg231=reg110*reg81; T reg232=reg119*reg71; T reg233=reg121*reg71; T reg234=reg110*reg91; T reg235=reg102*reg91;
    T reg236=reg117*reg71; reg152=reg18*reg152; reg232=reg234-reg232; reg218=reg217+reg218; reg212=reg211+reg212;
    reg193=reg192+reg193; reg150=reg18*reg150; reg201=reg200+reg201; reg195=reg194+reg195; reg216=reg215+reg216;
    reg181=reg179-reg181; reg197=reg196+reg197; reg214=reg213+reg214; reg120=reg18*reg120; reg199=reg130+reg199;
    reg198=reg123+reg198; reg156=reg18*reg156; reg231=reg210+reg231; reg182=reg230+reg182; reg183=reg184+reg183;
    reg169=reg170-reg169; reg228=reg229+reg228; reg185=reg186+reg185; reg236=reg168-reg236; reg226=reg227+reg226;
    reg224=reg225+reg224; reg172=reg173-reg172; reg162=reg163-reg162; reg123=reg187+reg188; reg220=reg219-reg220;
    reg130=reg190+reg221; reg112=reg133+reg112; reg191=reg189+reg191; reg176=reg180-reg176; reg203=reg202+reg203;
    reg164=reg165-reg164; reg133=reg161+reg222; reg205=reg204+reg205; reg233=reg235-reg233; reg166=reg167-reg166;
    reg209=reg208+reg209; reg163=reg171+reg223; reg175=reg174-reg175; reg207=reg206-reg207; reg178=reg177-reg178;
    reg162=reg18*reg162; reg152=ponderation*reg152; reg175=reg175*reg18; reg165=reg18*reg123; reg176=reg176*reg18;
    reg112=reg112*reg18; reg191=reg191*reg18; reg167=reg18*reg130; reg209=reg209*reg18; reg164=reg164*reg18;
    reg218=reg218*reg18; reg168=reg18*reg133; reg224=reg224*reg18; reg207=reg207*reg18; reg172=reg172*reg18;
    reg226=reg226*reg18; reg236=reg236*reg18; reg185=reg18*reg185; reg228=reg228*reg18; reg170=reg18*reg163;
    reg169=reg169*reg18; reg182=reg182*reg18; reg183=reg18*reg183; reg220=reg220*reg18; reg231=reg231*reg18;
    reg233=reg233*reg18; reg178=reg178*reg18; reg201=reg201*reg18; reg120=ponderation*reg120; reg212=reg212*reg18;
    reg199=reg199*reg18; reg203=reg203*reg18; reg156=ponderation*reg156; reg198=reg198*reg18; reg205=reg205*reg18;
    reg214=reg214*reg18; reg232=reg232*reg18; reg197=reg197*reg18; reg150=ponderation*reg150; reg193=reg193*reg18;
    reg166=reg166*reg18; reg181=reg181*reg18; reg216=reg216*reg18; reg195=reg195*reg18; T tmp_6_7=ponderation*reg164;
    T tmp_4_6=ponderation*reg198; sollicitation[indices[1]+1]+=-reg120; T tmp_5_6=ponderation*reg203; T tmp_0_0=ponderation*reg185; T tmp_0_3=ponderation*reg228;
    T tmp_2_5=ponderation*reg181; T tmp_1_6=ponderation*reg218; reg120=ponderation*reg170; sollicitation[indices[2]+1]+=reg120; T tmp_4_7=ponderation*reg199;
    T tmp_1_3=ponderation*reg212; T tmp_0_2=ponderation*reg182; T tmp_0_7=ponderation*reg193; T tmp_3_6=ponderation*reg233; T tmp_4_4=ponderation*reg209;
    T tmp_0_1=ponderation*reg183; T tmp_6_6=ponderation*reg166; T tmp_3_3=ponderation*reg169; T tmp_4_5=ponderation*reg231; T tmp_5_5=ponderation*reg201;
    T tmp_2_2=ponderation*reg220; T tmp_1_7=ponderation*reg112; T tmp_1_1=ponderation*reg195; T tmp_2_6=ponderation*reg176; reg112=ponderation*reg165;
    sollicitation[indices[0]+0]+=reg112; sollicitation[indices[1]+0]+=-reg150; T tmp_0_6=ponderation*reg191; T tmp_1_2=ponderation*reg197; T tmp_3_7=ponderation*reg207;
    sollicitation[indices[3]+1]+=-reg152; T tmp_1_5=ponderation*reg216; T tmp_3_4=ponderation*reg236; T tmp_7_7=ponderation*reg162; T tmp_5_7=ponderation*reg205;
    sollicitation[indices[3]+0]+=-reg156; T tmp_0_5=ponderation*reg224; reg150=ponderation*reg168; sollicitation[indices[2]+0]+=reg150; T tmp_1_4=ponderation*reg214;
    reg152=ponderation*reg167; sollicitation[indices[0]+1]+=reg152; T tmp_2_3=ponderation*reg175; T tmp_2_7=ponderation*reg172; T tmp_0_4=ponderation*reg226;
    T tmp_3_5=ponderation*reg232; T tmp_2_4=ponderation*reg178;
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
    T reg3=pow((*f.m).v2[1],2); reg3=reg2+reg3; reg2=pow((*f.m).v2[2],2); reg0=reg1+reg0; reg1=2*(*f.m).shear_modulus_13;
    T reg4=2*(*f.m).shear_modulus_23; reg0=pow(reg0,0.5); reg2=reg3+reg2; reg3=2*(*f.m).shear_modulus_12; reg1=1.0/reg1;
    reg4=1.0/reg4; reg2=pow(reg2,0.5); T reg5=reg1*reg4; T reg6=(*f.m).v1[1]/reg0; T reg7=(*f.m).v1[0]/reg0;
    reg3=1.0/reg3; T reg8=(*f.m).v2[0]/reg2; T reg9=(*f.m).v2[1]/reg2; T reg10=2*reg6; T reg11=2*reg7;
    T reg12=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg13=1.0/(*f.m).elastic_modulus_3; T reg14=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg15=reg3*reg5; reg0=(*f.m).v1[2]/reg0;
    reg2=(*f.m).v2[2]/reg2; T reg16=1.0/(*f.m).elastic_modulus_1; T reg17=2*reg0; T reg18=reg9*reg10; T reg19=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1;
    T reg20=1.0/(*f.m).elastic_modulus_2; T reg21=pow(reg8,2); T reg22=reg8*reg11; T reg23=pow(reg9,2); T reg24=reg12*reg15;
    T reg25=reg13*reg15; T reg26=reg14*reg15; T reg27=reg18*reg20; reg17=reg2*reg17; T reg28=reg14*reg24;
    T reg29=reg21*reg16; T reg30=reg19*reg25; T reg31=reg14*reg26; T reg32=reg22*reg19; T reg33=reg18*reg19;
    T reg34=pow(reg2,2); T reg35=reg20*reg25; T reg36=reg23*reg20; T reg37=reg22*reg16; T reg38=reg21*reg19;
    T reg39=reg23*reg19; reg31=reg35-reg31; reg28=reg30+reg28; reg35=reg19*reg26; T reg40=reg20*reg24;
    T reg41=reg34*reg12; reg39=reg29-reg39; reg29=reg17*reg14; T reg42=reg18*reg14; reg38=reg36-reg38;
    reg36=reg34*reg14; T reg43=reg22*reg12; reg33=reg37-reg33; reg37=reg17*reg12; T reg44=reg23*reg14;
    T reg45=reg21*reg12; T reg46=pow(reg7,2); T reg47=pow(reg6,2); reg32=reg27-reg32; reg27=reg34*reg13;
    T reg48=reg47*reg19; reg44=reg45+reg44; reg29=reg32-reg29; reg32=reg35+reg40; reg45=reg19*reg28;
    reg17=reg17*reg13; reg42=reg43+reg42; reg43=reg16*reg31; reg36=reg38-reg36; reg38=reg47*reg20;
    T reg49=reg46*reg19; reg37=reg33-reg37; reg41=reg39-reg41; reg33=reg6*reg8; reg39=reg7*reg9;
    T reg50=reg6*reg9; T reg51=reg7*reg8; T reg52=pow(reg0,2); T reg53=reg46*reg16; T reg54=reg12*reg24;
    T reg55=reg47*reg14; reg25=reg16*reg25; T reg56=reg36*reg47; reg44=reg27-reg44; reg42=reg17-reg42;
    reg17=reg12*reg26; reg27=reg37*reg21; T reg57=reg7*reg2; T reg58=reg20*reg15; T reg59=reg12*reg32;
    T reg60=reg33+reg39; T reg61=reg0*reg2; T reg62=reg29*reg23; T reg63=reg41*reg46; reg45=reg43-reg45;
    reg43=2*reg8; T reg64=reg0*reg8; T reg65=reg0*reg9; T reg66=reg6*reg2; T reg67=reg8*reg9;
    T reg68=reg52*reg12; reg48=reg53-reg48; reg53=reg13*reg5; T reg69=reg52*reg14; reg49=reg38-reg49;
    reg38=reg29*reg47; T reg70=reg12*reg5; reg5=reg14*reg5; T reg71=reg41*reg21; T reg72=reg36*reg23;
    T reg73=reg3*reg4; reg29=reg50*reg29; T reg74=reg51*reg37; reg37=reg37*reg46; reg15=reg19*reg15;
    reg36=reg50*reg36; reg41=reg51*reg41; T reg75=reg46*reg12; T reg76=reg44*reg52; reg56=reg63+reg56;
    reg38=reg37+reg38; reg37=reg6*reg11; reg63=reg13*reg73; T reg77=reg42*reg52; T reg78=reg8*reg2;
    T reg79=reg14*reg70; reg59=reg45-reg59; reg45=reg14*reg5; T reg80=reg19*reg53; T reg81=reg12*reg58;
    reg17=reg30+reg17; reg54=reg25-reg54; reg26=reg16*reg26; reg25=reg12*reg15; reg30=reg14*reg73;
    reg24=reg19*reg24; T reg82=reg3*reg1; reg53=reg20*reg53; reg73=reg12*reg73; reg55=reg75+reg55;
    reg75=reg67*reg3; T reg83=reg60*reg3; T reg84=reg57+reg64; T reg85=reg52*reg13; T reg86=reg43*reg9;
    reg36=reg41+reg36; reg41=reg61*reg44; T reg87=reg66-reg65; reg29=reg74+reg29; reg74=reg61*reg42;
    reg42=reg34*reg42; reg68=reg48-reg68; reg62=reg27+reg62; reg69=reg49-reg69; reg72=reg71+reg72;
    reg44=reg34*reg44; reg70=reg20*reg70; reg81=reg35+reg81; reg28=reg28/reg59; reg54=reg54/reg59;
    reg55=reg85-reg55; reg27=reg14*reg82; reg41=reg36+reg41; reg5=reg19*reg5; reg36=reg60*reg75;
    reg25=reg26+reg25; reg48=reg60*reg83; reg74=reg29+reg74; reg24=reg26+reg24; reg15=reg19*reg15;
    reg77=reg38+reg77; reg26=reg83*reg37; reg29=reg68*reg21; reg38=reg69*reg23; reg49=reg75*reg37;
    reg76=reg56+reg76; reg44=reg72+reg44; reg75=reg86*reg75; reg56=reg69*reg47; reg71=reg68*reg46;
    reg72=reg84*reg1; reg85=reg78*reg1; reg42=reg62+reg42; reg62=reg0*reg11; T reg88=reg7*reg6;
    reg83=reg86*reg83; T reg89=reg9*reg2; reg57=reg64-reg57; reg58=reg16*reg58; reg31=reg31/reg59;
    reg64=reg43*reg2; T reg90=2*reg9; T reg91=2*reg87; reg66=reg65+reg66; reg17=reg17/reg59;
    reg13=reg13*reg82; reg45=reg53-reg45; reg79=reg80+reg79; reg82=reg12*reg82; reg53=reg14*reg73;
    reg65=reg20*reg63; reg63=reg19*reg63; reg80=reg14*reg30; T reg92=reg20*reg13; reg13=reg19*reg13;
    T reg93=reg14*reg27; reg26=reg77+reg26; reg77=reg72*reg62; reg38=reg29+reg38; reg29=reg34*reg55;
    reg14=reg14*reg82; reg75=reg44+reg75; reg44=reg85*reg64; T reg94=reg85*reg62; reg49=reg76+reg49;
    reg53=reg63+reg53; reg63=reg55*reg52; reg56=reg71+reg56; reg80=reg65-reg80; reg65=reg66*reg4;
    reg71=reg89*reg4; reg76=reg0*reg10; T reg95=reg5+reg70; T reg96=reg7*reg0; reg79=reg19*reg79;
    T reg97=1-var_inter[1]; T reg98=1-var_inter[0]; T reg99=reg72*reg84; reg48=reg74+reg48; reg85=reg85*reg84;
    reg36=reg41+reg36; reg69=reg50*reg69; reg68=reg51*reg68; reg41=reg86*reg54; reg3=reg88*reg3;
    reg74=reg37*reg17; T reg100=reg23*reg54; T reg101=reg47*reg17; T reg102=reg91*reg57; T reg103=pow(reg57,2);
    T reg104=pow(reg87,2); T reg105=reg90*reg2; T reg106=reg21*reg54; T reg107=reg46*reg17; T reg108=reg86*reg28;
    T reg109=reg37*reg31; reg72=reg72*reg64; reg83=reg42+reg83; reg42=reg23*reg28; T reg110=reg47*reg31;
    reg25=reg25/reg59; reg30=reg19*reg30; T reg111=reg21*reg28; reg45=reg16*reg45; reg81=reg81/reg59;
    T reg112=reg46*reg31; reg24=reg24/reg59; reg73=reg20*reg73; reg15=reg58-reg15; reg32=reg32/reg59;
    reg58=reg3*reg37; reg53=reg19*reg53; reg63=reg56+reg63; reg100=reg101+reg100; reg94=reg49+reg94;
    reg49=reg71*reg76; reg56=reg30+reg73; reg80=reg16*reg80; reg82=reg20*reg82; reg1=reg96*reg1;
    reg20=elem.pos(1)[1]*var_inter[0]; reg101=reg98*elem.pos(0)[1]; T reg113=reg97*elem.pos(0)[1]; T reg114=reg97*elem.pos(1)[1]; reg69=reg68+reg69;
    reg55=reg61*reg55; reg68=reg47*reg81; T reg115=reg23*reg25; reg85=reg36+reg85; reg36=reg71*reg66;
    T reg116=reg37*reg81; T reg117=reg86*reg25; T reg118=reg104*reg24; reg106=reg107+reg106; reg15=reg15/reg59;
    reg41=reg74+reg41; reg74=reg102*reg24; reg107=reg102*reg32; reg108=reg109+reg108; reg109=reg46*reg81;
    T reg119=reg21*reg25; T reg120=reg65*reg105; reg72=reg83+reg72; reg83=reg103*reg32; reg42=reg110+reg42;
    reg110=reg98*elem.pos(0)[0]; reg71=reg71*reg105; reg44=reg75+reg44; reg14=reg13+reg14; reg13=reg86*reg3;
    reg29=reg38+reg29; reg93=reg92-reg93; reg38=reg65*reg76; reg77=reg26+reg77; reg27=reg19*reg27;
    reg26=elem.pos(1)[0]*var_inter[0]; reg75=reg97*elem.pos(0)[0]; reg92=reg6*reg0; T reg121=reg97*elem.pos(1)[0]; T reg122=reg103*reg24;
    reg65=reg65*reg66; reg99=reg48+reg99; reg95=reg12*reg95; reg111=reg112+reg111; reg48=reg104*reg32;
    reg79=reg45-reg79; reg53=reg80-reg53; reg45=reg1*reg62; reg58=reg63+reg58; reg56=reg12*reg56;
    reg63=elem.pos(2)[0]*var_inter[0]; reg121=reg121-reg75; reg49=reg94+reg49; reg80=elem.pos(2)[0]*var_inter[1]; reg94=reg110+reg26;
    reg4=reg92*reg4; reg112=elem.pos(2)[1]*var_inter[0]; reg114=reg114-reg113; reg107=reg108+reg107; reg74=reg41+reg74;
    reg83=reg42+reg83; reg48=reg111+reg48; reg120=reg72+reg120; reg119=reg109+reg119; reg122=reg100+reg122;
    reg41=reg27+reg82; reg71=reg44+reg71; reg14=reg19*reg14; reg118=reg106+reg118; reg19=reg1*reg64;
    reg13=reg29+reg13; reg93=reg16*reg93; reg38=reg77+reg38; reg115=reg68+reg115; reg16=reg103*reg15;
    reg29=elem.pos(2)[1]*var_inter[1]; reg42=reg104*reg15; reg44=reg20+reg101; reg36=reg85+reg36; reg65=reg99+reg65;
    reg55=reg69+reg55; reg95=reg79-reg95; reg117=reg116+reg117; reg102=reg102*reg15; reg3=reg60*reg3;
    reg68=reg57*reg87; reg69=elem.pos(3)[1]*var_inter[1]; reg42=reg119+reg42; reg41=reg12*reg41; reg56=reg53-reg56;
    reg14=reg93-reg14; reg1=reg1*reg84; reg12=reg49*reg65; reg114=reg29+reg114; reg29=reg88*reg107;
    reg53=reg4*reg105; reg72=reg88*reg48; reg77=reg67*reg118; reg19=reg13+reg19; reg3=reg55+reg3;
    reg13=reg67*reg74; reg55=reg71*reg65; reg79=reg120*reg36; reg85=reg7*reg57; reg93=reg6*reg87;
    reg102=reg117+reg102; reg99=reg67*reg122; reg100=reg88*reg83; reg95=reg95/reg59; reg16=reg115+reg16;
    reg106=reg38*reg36; reg108=elem.pos(3)[0]*var_inter[1]; reg80=reg121+reg80; reg112=reg112-reg44; reg63=reg63-reg94;
    reg109=elem.pos(3)[0]*reg98; reg111=elem.pos(3)[1]*reg98; reg115=reg4*reg76; reg45=reg58+reg45; reg58=reg8*reg57;
    reg80=reg80-reg108; reg4=reg4*reg66; reg41=reg14-reg41; reg114=reg114-reg69; reg14=reg47*reg83;
    reg99=reg100+reg99; reg100=reg68*reg16; reg116=reg47*reg48; reg117=reg23*reg118; reg1=reg3+reg1;
    reg79=reg55-reg79; reg3=reg47*reg107; reg55=reg21*reg122; reg109=reg63+reg109; reg63=reg51*reg95;
    reg119=reg6*reg57; reg121=reg50*reg95; T reg123=reg23*reg74; T reg124=reg7*reg87; reg118=reg21*reg118;
    T reg125=reg60*reg95; reg83=reg46*reg83; reg74=reg21*reg74; reg115=reg45+reg115; reg56=reg56/reg59;
    reg48=reg46*reg48; reg111=reg112+reg111; reg85=reg93+reg85; reg106=reg12-reg106; reg12=reg49*reg120;
    reg45=reg68*reg42; reg107=reg46*reg107; reg93=reg38*reg71; reg53=reg19+reg53; reg19=reg9*reg87;
    reg112=reg68*reg102; reg13=reg29+reg13; reg77=reg72+reg77; reg122=reg23*reg122; reg29=reg103*reg42;
    reg117=reg116+reg117; reg122=reg14+reg122; reg58=reg19+reg58; reg74=reg107+reg74; reg14=reg104*reg102;
    reg19=reg9*reg57; reg72=reg8*reg87; reg107=reg85*reg56; reg116=reg60*reg63; reg45=reg77+reg45;
    reg77=reg115*reg79; reg119=reg119*reg56; reg124=reg124*reg56; reg123=reg3+reg123; reg102=reg103*reg102;
    reg100=reg99+reg100; reg4=reg1+reg4; reg93=reg12-reg93; reg1=reg60*reg121; reg3=reg114*reg109;
    reg12=reg80*reg111; reg55=reg83+reg55; reg83=reg104*reg16; reg42=reg104*reg42; reg118=reg48+reg118;
    reg48=reg53*reg106; reg99=reg60*reg125; reg112=reg13+reg112; reg59=reg41/reg59; reg16=reg103*reg16;
    reg29=reg117+reg29; reg19=reg19*reg59; reg13=reg58*reg59; reg41=reg18*reg63; reg72=reg72*reg59;
    reg102=reg123+reg102; reg117=reg18*reg125; reg3=reg12-reg3; reg12=reg18*reg121; reg16=reg122+reg16;
    reg116=reg45+reg116; reg99=reg112+reg99; reg45=reg53*reg36; reg112=reg85*reg124; reg122=reg85*reg107;
    reg123=reg71*reg4; T reg126=reg4*reg93; reg48=reg77-reg48; reg36=reg115*reg36; reg121=reg22*reg121;
    reg77=reg49*reg4; reg83=reg55+reg83; reg42=reg118+reg42; reg63=reg22*reg63; reg10=reg57*reg10;
    reg14=reg74+reg14; reg125=reg22*reg125; reg71=reg115*reg71; reg49=reg49*reg53; reg1=reg100+reg1;
    reg11=reg87*reg11; reg55=reg85*reg119; reg114=reg114/reg3; reg126=reg48+reg126; reg109=reg109/reg3;
    reg80=reg80/reg3; reg48=reg58*reg19; reg111=reg111/reg3; reg74=reg53*reg65; reg100=reg120*reg4;
    reg65=reg115*reg65; reg123=reg45-reg123; reg63=reg42+reg63; reg4=reg38*reg4; reg42=reg11*reg124;
    reg90=reg90*reg57; reg77=reg36-reg77; reg120=reg115*reg120; reg53=reg38*reg53; reg49=reg71-reg49;
    reg43=reg43*reg87; reg36=reg10*reg107; reg117=reg102+reg117; reg41=reg29+reg41; reg124=reg10*reg124;
    reg107=reg11*reg107; reg125=reg14+reg125; reg14=reg58*reg72; reg12=reg16+reg12; reg112=reg116+reg112;
    reg16=reg10*reg119; reg119=reg11*reg119; reg121=reg83+reg121; reg122=reg99+reg122; reg55=reg1+reg55;
    reg1=reg58*reg13; reg29=var_inter[1]*reg111; reg38=reg43*reg19; reg14=reg112+reg14; reg48=reg55+reg48;
    reg45=reg97*reg111; reg53=reg120-reg53; reg107=reg125+reg107; reg49=reg49/reg126; reg55=reg90*reg72;
    reg71=1-(*f.m).resolution; reg42=reg63+reg42; reg124=reg41+reg124; reg41=var_inter[1]*reg109; reg63=reg43*reg13;
    reg72=reg43*reg72; reg83=var_inter[0]*reg80; reg99=reg98*reg114; reg19=reg90*reg19; reg16=reg12+reg16;
    reg100=reg74-reg100; reg12=var_inter[0]*reg114; reg123=reg123/reg126; reg1=reg122+reg1; reg119=reg121+reg119;
    reg36=reg117+reg36; reg13=reg90*reg13; reg4=reg65-reg4; reg77=reg77/reg126; reg65=reg97*reg109;
    reg74=reg98*reg80; reg19=reg16+reg19; reg16=(*f.m).resolution*reg49; reg63=reg107+reg63; reg55=reg124+reg55;
    reg38=reg119+reg38; reg102=reg99-reg45; reg107=reg65-reg74; reg112=reg45+reg12; reg115=reg65+reg83;
    reg72=reg42+reg72; reg42=reg29-reg12; reg116=reg83-reg41; reg117=reg99+reg29; reg118=reg41+reg74;
    reg79=reg79/reg126; reg100=reg100/reg126; reg106=reg106/reg126; reg4=reg4/reg126; reg93=reg93/reg126;
    reg126=reg53/reg126; reg53=reg14*reg71; reg119=reg48*reg71; reg120=reg71*reg1; reg121=(*f.m).resolution*reg123;
    reg122=(*f.m).resolution*reg77; reg13=reg36+reg13; reg36=0.5*reg116; reg124=0.5*reg42; reg125=0.5*reg118;
    T reg127=0.5*reg117; T reg128=(*f.m).resolution*reg106; T reg129=(*f.m).resolution*reg79; T reg130=(*f.m).resolution*reg93; T reg131=reg71*reg72;
    T reg132=reg71*reg38; T reg133=reg71*reg63; T reg134=reg71*reg55; T reg135=reg71*reg19; T reg136=reg13*reg71;
    T reg137=(*f.m).resolution*reg100; T reg138=0.5*reg112; T reg139=0.5*reg115; T reg140=0.5*reg102; T reg141=0.5*reg107;
    reg16=reg120+reg16; reg122=reg119-reg122; reg121=reg53+reg121; reg53=(*f.m).resolution*reg126; reg119=(*f.m).resolution*reg4;
    reg53=reg136-reg53; reg119=reg135+reg119; reg137=reg134-reg137; reg133=reg130+reg133; reg128=reg132-reg128;
    reg131=reg129+reg131; reg120=reg127*reg16; reg129=reg118*reg122; reg130=reg122*reg115; reg132=reg16*reg138;
    reg134=reg125*reg16; reg135=reg42*reg121; reg136=reg36*reg16; T reg142=reg117*reg121; T reg143=reg16*reg139;
    T reg144=reg16*reg141; T reg145=reg116*reg122; T reg146=reg124*reg16; T reg147=reg122*reg107; T reg148=reg121*reg112;
    T reg149=reg16*reg140; T reg150=reg121*reg102; T reg151=reg128*reg107; T reg152=reg117*reg131; T reg153=reg125*reg133;
    T reg154=reg118*reg119; T reg155=reg119*reg115; reg148=reg148-reg143; reg144=reg150+reg144; reg129=reg129-reg120;
    reg150=reg36*reg53; T reg156=reg42*reg137; T reg157=reg127*reg53; T reg158=reg118*reg128; T reg159=reg127*reg133;
    T reg160=reg137*reg102; T reg161=reg124*reg53; T reg162=reg116*reg119; T reg163=reg53*reg141; T reg164=reg125*reg53;
    T reg165=reg117*reg137; T reg166=reg124*reg133; T reg167=reg116*reg128; T reg168=reg137*reg112; T reg169=reg53*reg139;
    T reg170=reg133*reg139; T reg171=reg131*reg112; reg136=reg135+reg136; reg135=reg133*reg141; reg149=reg147+reg149;
    reg147=reg36*reg133; T reg172=reg42*reg131; reg146=reg145+reg146; reg134=reg134-reg142; reg145=reg119*reg107;
    reg132=reg132-reg130; T reg173=reg53*reg140; T reg174=reg131*reg102; T reg175=reg133*reg140; T reg176=reg53*reg138;
    T reg177=reg133*reg138; T reg178=reg128*reg115; reg171=reg171-reg170; reg154=reg154-reg157; reg164=reg164-reg165;
    reg162=reg161+reg162; reg160=reg163+reg160; reg156=reg150+reg156; reg176=reg176-reg155; reg144=2*reg144;
    reg175=reg151+reg175; reg149=2*reg149; reg145=reg173+reg145; reg168=reg168-reg169; reg166=reg167+reg166;
    reg136=2*reg136; reg147=reg172+reg147; reg132=2*reg132; reg146=2*reg146; reg177=reg177-reg178;
    reg134=2*reg134; reg148=2*reg148; reg158=reg158-reg159; reg135=reg174+reg135; reg129=2*reg129;
    reg153=reg153-reg152; reg150=reg145*reg115; reg151=reg132*reg125; reg161=reg149*reg138; reg163=reg160*reg115;
    reg167=reg144*reg138; reg172=reg177*reg117; reg173=reg146*reg124; reg174=reg145*reg116; T reg179=reg149*reg124;
    T reg180=reg129*reg139; T reg181=reg136*reg125; T reg182=reg162*reg116; T reg183=reg147*reg117; T reg184=reg175*reg117;
    T reg185=reg136*reg127; T reg186=reg158*reg112; T reg187=reg134*reg139; T reg188=reg149*reg125; T reg189=reg154*reg115;
    T reg190=reg129*reg138; T reg191=reg148*reg140; T reg192=reg164*reg115; T reg193=reg134*reg138; T reg194=reg171*reg117;
    T reg195=reg148*reg125; T reg196=reg153*reg42; T reg197=reg134*reg36; T reg198=reg162*reg115; T reg199=reg146*reg138;
    T reg200=reg156*reg115; T reg201=reg136*reg138; T reg202=reg158*reg42; T reg203=reg129*reg36; T reg204=reg176*reg115;
    T reg205=reg132*reg138; T reg206=reg168*reg115; T reg207=reg160*reg116; T reg208=reg144*reg124; T reg209=reg144*reg127;
    T reg210=reg129*reg124; T reg211=reg132*reg139; T reg212=reg177*reg112; T reg213=reg129*reg127; T reg214=reg154*reg116;
    T reg215=reg145*reg118; T reg216=reg149*reg127; T reg217=reg148*reg127; T reg218=reg162*reg118; T reg219=reg154*reg118;
    T reg220=reg156*reg116; T reg221=reg136*reg124; T reg222=reg134*reg124; T reg223=reg134*reg127; T reg224=reg164*reg116;
    T reg225=reg176*reg118; T reg226=reg132*reg127; T reg227=reg164*reg118; T reg228=reg146*reg127; T reg229=reg146*reg125;
    T reg230=reg153*reg112; T reg231=reg168*reg116; T reg232=reg156*reg118; T reg233=reg171*reg102; T reg234=reg166*reg117;
    T reg235=reg134*reg125; T reg236=reg148*reg124; T reg237=reg146*reg139; T reg238=reg153*reg117; T reg239=reg135*reg117;
    T reg240=reg144*reg125; T reg241=reg129*reg125; T reg242=reg158*reg117; T reg243=reg166*reg112; T reg244=reg136*reg139;
    T reg245=reg147*reg112; T reg246=reg176*reg116; T reg247=reg132*reg124; T reg248=reg160*reg118; T reg249=reg147*reg42;
    T reg250=reg168*reg118; T reg251=reg136*reg36; T reg252=reg146*reg36; T reg253=reg149*reg141; T reg254=reg166*reg42;
    T reg255=reg148*reg139; T reg256=reg175*reg102; T reg257=reg171*reg112; T reg258=reg149*reg139; T reg259=reg144*reg141;
    reg176=reg176*reg107; T reg260=reg175*reg112; T reg261=reg132*reg140; reg158=reg158*reg102; T reg262=reg144*reg140;
    reg160=reg160*reg107; reg156=reg156*reg107; T reg263=reg136*reg140; reg162=reg162*reg107; T reg264=reg146*reg140;
    T reg265=reg144*reg139; T reg266=reg135*reg112; reg164=reg164*reg107; T reg267=reg129*reg140; T reg268=reg134*reg140;
    T reg269=reg135*reg102; reg154=reg154*reg107; reg129=reg129*reg141; T reg270=reg148*reg138; reg147=reg147*reg102;
    reg171=reg171*reg42; reg136=reg136*reg141; T reg271=reg148*reg36; T reg272=reg132*reg141; T reg273=reg149*reg36;
    reg146=reg146*reg141; T reg274=reg177*reg102; reg149=reg149*reg140; reg134=reg134*reg141; reg175=reg175*reg42;
    reg145=reg145*reg107; reg148=reg148*reg141; reg168=reg168*reg107; reg135=reg135*reg42; reg132=reg132*reg36;
    reg153=reg153*reg102; reg166=reg166*reg102; reg177=reg177*reg42; reg144=reg144*reg36; reg185=reg232-reg185;
    reg173=reg182+reg173; reg267=reg154+reg267; reg194=reg195-reg194; reg226=reg225-reg226; reg146=reg166+reg146;
    reg213=reg219-reg213; reg172=reg151-reg172; reg136=reg147+reg136; reg238=reg235-reg238; reg233=reg148+reg233;
    reg134=reg153+reg134; reg252=reg254+reg252; reg242=reg241-reg242; reg210=reg214+reg210; reg239=reg240-reg239;
    reg234=reg229-reg234; reg255=reg257-reg255; reg209=reg248-reg209; reg258=reg260-reg258; reg216=reg215-reg216;
    reg222=reg224+reg222; reg217=reg250-reg217; reg183=reg181-reg183; reg272=reg274+reg272; reg184=reg188-reg184;
    reg265=reg266-reg265; reg187=reg230-reg187; reg179=reg174+reg179; reg180=reg186-reg180; reg253=reg256+reg253;
    reg251=reg249+reg251; reg163=reg167-reg163; reg208=reg207+reg208; reg150=reg161-reg150; reg132=reg177+reg132;
    reg206=reg270-reg206; reg203=reg202+reg203; reg149=reg145+reg149; reg204=reg205-reg204; reg271=reg171+reg271;
    reg200=reg201-reg200; reg197=reg196+reg197; reg273=reg175+reg273; reg198=reg199-reg198; reg144=reg135+reg144;
    reg192=reg193-reg192; reg191=reg168+reg191; reg189=reg190-reg189; reg263=reg156+reg263; reg247=reg246+reg247;
    reg264=reg162+reg264; reg211=reg212-reg211; reg261=reg176+reg261; reg244=reg245-reg244; reg223=reg227-reg223;
    reg236=reg231+reg236; reg221=reg220+reg221; reg228=reg218-reg228; reg237=reg243-reg237; reg259=reg269+reg259;
    reg262=reg160+reg262; reg268=reg164+reg268; reg129=reg158+reg129; reg129=reg129*reg3; reg184=reg184*reg3;
    reg221=reg221*reg3; reg197=reg197*reg3; reg149=reg149*reg3; reg173=reg173*reg3; reg247=reg247*reg3;
    reg203=reg203*reg3; reg239=reg239*reg3; reg236=reg236*reg3; reg208=reg208*reg3; reg134=reg134*reg3;
    reg222=reg222*reg3; reg210=reg210*reg3; reg179=reg179*reg3; reg262=reg262*reg3; reg223=reg223*reg3;
    reg265=reg265*reg3; reg213=reg3*reg213; reg255=reg255*reg3; reg211=reg211*reg3; reg267=reg267*reg3;
    reg244=reg244*reg3; reg237=reg237*reg3; reg268=reg268*reg3; reg187=reg187*reg3; reg180=reg180*reg3;
    reg264=reg264*reg3; reg163=reg163*reg3; reg263=reg263*reg3; reg150=reg150*reg3; reg261=reg261*reg3;
    reg206=reg206*reg3; reg204=reg204*reg3; reg191=reg191*reg3; reg200=reg200*reg3; reg198=reg198*reg3;
    reg251=reg251*reg3; reg192=reg192*reg3; reg189=reg189*reg3; reg132=reg132*reg3; reg144=reg144*reg3;
    reg273=reg273*reg3; reg271=reg271*reg3; reg194=reg194*reg3; reg253=reg3*reg253; reg172=reg172*reg3;
    reg183=reg183*reg3; reg259=reg3*reg259; reg234=reg234*reg3; reg146=reg146*reg3; reg238=reg238*reg3;
    reg136=reg136*reg3; reg242=reg242*reg3; reg209=reg209*reg3; reg272=reg272*reg3; reg216=reg216*reg3;
    reg228=reg228*reg3; reg258=reg258*reg3; reg185=reg185*reg3; reg226=reg226*reg3; reg252=reg252*reg3;
    reg233=reg233*reg3; reg217=reg217*reg3; T tmp_0_1=ponderation*reg253; T tmp_4_2=ponderation*reg271; T tmp_2_0=ponderation*reg265;
    T tmp_1_5=ponderation*reg264; T tmp_0_5=ponderation*reg146; T tmp_4_3=ponderation*reg132; T tmp_2_1=ponderation*reg258; T tmp_0_0=ponderation*reg259;
    T tmp_4_4=ponderation*reg251; T tmp_0_6=ponderation*reg134; T tmp_1_0=ponderation*reg262; T tmp_7_7=ponderation*reg213; T tmp_0_7=ponderation*reg129;
    T tmp_1_2=ponderation*reg191; T tmp_0_4=ponderation*reg136; T tmp_1_7=ponderation*reg267; T tmp_1_3=ponderation*reg261; T tmp_0_2=ponderation*reg233;
    T tmp_0_3=ponderation*reg272; T tmp_4_5=ponderation*reg252; T tmp_1_4=ponderation*reg263; T tmp_1_6=ponderation*reg268; T tmp_2_2=ponderation*reg255;
    T tmp_5_3=ponderation*reg247; T tmp_7_6=ponderation*reg223; T tmp_7_5=ponderation*reg228; T tmp_5_4=ponderation*reg221; T tmp_7_4=ponderation*reg185;
    T tmp_7_3=ponderation*reg226; T tmp_5_5=ponderation*reg173; T tmp_7_2=ponderation*reg217; T tmp_7_1=ponderation*reg216; T tmp_5_6=ponderation*reg222;
    T tmp_7_0=ponderation*reg209; T tmp_6_7=ponderation*reg242; T tmp_5_7=ponderation*reg210; T tmp_6_6=ponderation*reg238; T tmp_6_5=ponderation*reg234;
    T tmp_6_0=ponderation*reg239; T tmp_6_4=ponderation*reg183; T tmp_6_3=ponderation*reg172; T tmp_6_1=ponderation*reg184; T tmp_6_2=ponderation*reg194;
    T tmp_4_1=ponderation*reg273; T tmp_1_1=ponderation*reg149; T tmp_4_0=ponderation*reg144; T tmp_3_7=ponderation*reg189; T tmp_3_6=ponderation*reg192;
    T tmp_3_5=ponderation*reg198; T tmp_3_4=ponderation*reg200; T tmp_4_6=ponderation*reg197; T tmp_3_3=ponderation*reg204; T tmp_3_2=ponderation*reg206;
    T tmp_4_7=ponderation*reg203; T tmp_3_1=ponderation*reg150; T tmp_3_0=ponderation*reg163; T tmp_5_0=ponderation*reg208; T tmp_2_7=ponderation*reg180;
    T tmp_2_6=ponderation*reg187; T tmp_5_1=ponderation*reg179; T tmp_2_5=ponderation*reg237; T tmp_2_4=ponderation*reg244; T tmp_5_2=ponderation*reg236;
    T tmp_2_3=ponderation*reg211;
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
    T reg0=pow((*f.m).v1[0],2); T reg1=pow((*f.m).v1[1],2); T reg2=pow((*f.m).v2[1],2); T reg3=pow((*f.m).v2[0],2); T reg4=pow((*f.m).v1[2],2);
    reg1=reg0+reg1; reg0=2*(*f.m).shear_modulus_23; reg4=reg1+reg4; reg1=2*(*f.m).shear_modulus_13; reg2=reg3+reg2;
    reg3=pow((*f.m).v2[2],2); reg0=1.0/reg0; reg1=1.0/reg1; reg4=pow(reg4,0.5); T reg5=2*(*f.m).shear_modulus_12;
    reg3=reg2+reg3; reg2=reg1*reg0; reg5=1.0/reg5; T reg6=(*f.m).v1[0]/reg4; T reg7=(*f.m).v1[1]/reg4;
    reg3=pow(reg3,0.5); T reg8=reg5*reg2; T reg9=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg10=1.0/(*f.m).elastic_modulus_3; T reg11=2*reg7;
    T reg12=2*reg6; T reg13=(*f.m).v2[1]/reg3; T reg14=(*f.m).v2[0]/reg3; reg4=(*f.m).v1[2]/reg4; T reg15=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg16=reg13*reg11; T reg17=2*reg4; T reg18=reg9*reg8; T reg19=reg14*reg12; T reg20=reg15*reg8;
    T reg21=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg22=pow(reg13,2); T reg23=pow(reg14,2); T reg24=1.0/(*f.m).elastic_modulus_2; reg3=(*f.m).v2[2]/reg3;
    T reg25=reg10*reg8; T reg26=1.0/(*f.m).elastic_modulus_1; T reg27=reg23*reg21; T reg28=reg19*reg26; T reg29=reg16*reg24;
    T reg30=reg16*reg21; T reg31=reg19*reg21; T reg32=reg22*reg24; T reg33=reg9*reg20; reg17=reg3*reg17;
    T reg34=reg21*reg25; T reg35=reg22*reg21; T reg36=reg23*reg26; T reg37=pow(reg3,2); T reg38=reg24*reg25;
    T reg39=reg9*reg18; T reg40=reg37*reg15; reg35=reg36-reg35; reg36=reg24*reg20; T reg41=reg21*reg18;
    T reg42=reg17*reg15; reg30=reg28-reg30; reg33=reg34+reg33; reg28=reg19*reg15; T reg43=reg16*reg9;
    T reg44=reg22*reg9; T reg45=reg23*reg15; reg39=reg38-reg39; reg38=pow(reg7,2); reg31=reg29-reg31;
    reg29=reg17*reg9; T reg46=pow(reg6,2); reg27=reg32-reg27; reg32=reg37*reg9; T reg47=pow(reg4,2);
    T reg48=reg26*reg39; T reg49=reg7*reg14; reg43=reg28+reg43; reg28=reg46*reg26; reg17=reg17*reg10;
    T reg50=reg6*reg13; T reg51=reg38*reg21; T reg52=reg6*reg14; T reg53=reg41+reg36; reg44=reg45+reg44;
    reg45=reg37*reg10; T reg54=reg7*reg13; T reg55=reg21*reg33; reg29=reg31-reg29; reg40=reg35-reg40;
    reg32=reg27-reg32; reg42=reg30-reg42; reg27=reg46*reg21; reg30=reg38*reg24; reg31=reg52*reg40;
    reg35=reg54*reg32; T reg56=reg40*reg46; T reg57=reg32*reg38; T reg58=reg42*reg46; T reg59=reg29*reg38;
    T reg60=reg15*reg2; T reg61=reg10*reg2; reg40=reg40*reg23; reg32=reg32*reg22; T reg62=reg42*reg23;
    T reg63=reg29*reg22; reg2=reg9*reg2; T reg64=2*reg14; T reg65=reg4*reg3; T reg66=reg49+reg50;
    T reg67=reg5*reg0; T reg68=reg21*reg8; T reg69=reg15*reg20; reg25=reg26*reg25; T reg70=reg15*reg18;
    reg8=reg24*reg8; T reg71=reg15*reg53; reg55=reg48-reg55; reg51=reg28-reg51; reg28=reg47*reg15;
    reg48=reg14*reg13; T reg72=reg6*reg3; T reg73=reg4*reg14; T reg74=reg47*reg9; T reg75=reg4*reg13;
    T reg76=reg7*reg3; reg43=reg17-reg43; reg44=reg45-reg44; reg29=reg54*reg29; reg42=reg52*reg42;
    reg27=reg30-reg27; reg17=reg38*reg9; reg30=reg46*reg15; reg17=reg30+reg17; reg30=reg47*reg10;
    reg45=reg10*reg67; T reg77=reg72+reg73; T reg78=reg64*reg13; T reg79=reg9*reg60; T reg80=reg48*reg5;
    reg20=reg21*reg20; T reg81=reg9*reg67; T reg82=reg15*reg68; T reg83=reg66*reg5; T reg84=reg9*reg2;
    T reg85=reg21*reg61; T reg86=reg37*reg43; reg63=reg62+reg63; reg62=reg37*reg44; reg32=reg40+reg32;
    reg40=reg76-reg75; reg61=reg24*reg61; reg67=reg15*reg67; T reg87=reg5*reg1; T reg88=reg7*reg12;
    T reg89=reg43*reg47; reg59=reg58+reg59; reg58=reg44*reg47; reg57=reg56+reg57; reg71=reg55-reg71;
    reg35=reg31+reg35; reg44=reg65*reg44; reg18=reg26*reg18; reg69=reg25-reg69; reg74=reg27-reg74;
    reg25=reg14*reg3; reg29=reg42+reg29; reg27=reg15*reg8; reg43=reg65*reg43; reg28=reg51-reg28;
    reg70=reg34+reg70; reg31=reg28*reg23; reg34=reg74*reg22; reg43=reg29+reg43; reg62=reg32+reg62;
    reg29=reg78*reg80; reg32=reg66*reg83; reg39=reg39/reg71; reg86=reg63+reg86; reg42=reg78*reg83;
    reg72=reg73-reg72; reg51=reg6*reg7; reg58=reg57+reg58; reg55=reg80*reg88; reg56=reg74*reg38;
    reg57=reg28*reg46; reg63=reg4*reg12; reg44=reg35+reg44; reg35=reg77*reg1; reg73=reg25*reg1;
    reg89=reg59+reg89; reg83=reg83*reg88; reg2=reg21*reg2; reg80=reg66*reg80; reg59=reg9*reg87;
    reg60=reg24*reg60; T reg90=reg15*reg87; T reg91=2*reg40; T reg92=reg9*reg67; reg27=reg41+reg27;
    reg33=reg33/reg71; T reg93=reg9*reg81; reg87=reg10*reg87; reg10=reg21*reg45; reg45=reg24*reg45;
    reg70=reg70/reg71; reg69=reg69/reg71; reg76=reg75+reg76; reg75=2*reg13; T reg94=reg64*reg3;
    reg79=reg85+reg79; reg8=reg26*reg8; reg17=reg30-reg17; reg84=reg61-reg84; reg82=reg18+reg82;
    reg20=reg18+reg20; reg68=reg21*reg68; reg18=reg13*reg3; reg30=reg23*reg33; reg61=reg91*reg72;
    reg85=reg46*reg39; reg74=reg54*reg74; reg53=reg53/reg71; reg82=reg82/reg71; reg20=reg20/reg71;
    reg55=reg58+reg55; reg58=reg73*reg63; T reg95=reg4*reg11; T reg96=reg35*reg94; reg42=reg86+reg42;
    reg86=reg35*reg77; reg32=reg43+reg32; reg43=reg73*reg94; reg29=reg62+reg29; reg62=reg6*reg4;
    reg27=reg27/reg71; reg5=reg51*reg5; T reg97=reg37*reg17; reg34=reg31+reg34; reg35=reg35*reg63;
    reg73=reg73*reg77; reg83=reg89+reg83; reg80=reg44+reg80; reg31=reg75*reg3; reg44=pow(reg40,2);
    reg89=pow(reg72,2); reg28=reg52*reg28; reg67=reg24*reg67; T reg98=reg38*reg39; T reg99=reg22*reg69;
    reg93=reg45-reg93; reg45=reg9*reg90; T reg100=1-var_inter[1]; T reg101=reg88*reg70; T reg102=reg78*reg69;
    reg92=reg10+reg92; reg9=reg9*reg59; reg10=1-var_inter[0]; T reg103=reg24*reg87; reg87=reg21*reg87;
    reg84=reg26*reg84; T reg104=reg23*reg69; T reg105=reg46*reg70; reg79=reg21*reg79; T reg106=reg78*reg33;
    T reg107=reg17*reg47; T reg108=reg88*reg39; T reg109=reg38*reg70; T reg110=reg2+reg60; reg56=reg57+reg56;
    reg57=reg76*reg0; T reg111=reg18*reg0; T reg112=reg22*reg33; reg68=reg8-reg68; reg81=reg21*reg81;
    reg92=reg21*reg92; reg8=elem.pos(1)[1]*var_inter[0]; T reg113=reg7*reg4; reg93=reg26*reg93; reg110=reg15*reg110;
    T reg114=reg81+reg67; T reg115=reg10*elem.pos(0)[1]; reg79=reg84-reg79; reg84=reg44*reg20; reg104=reg105+reg104;
    reg105=reg61*reg53; reg106=reg108+reg106; reg108=reg89*reg53; reg112=reg98+reg112; reg99=reg109+reg99;
    reg98=reg89*reg20; reg45=reg87+reg45; reg9=reg103-reg9; reg102=reg101+reg102; reg87=reg61*reg20;
    reg101=reg46*reg27; reg103=reg23*reg82; reg109=reg100*elem.pos(0)[1]; T reg116=reg100*elem.pos(1)[1]; T reg117=reg10*elem.pos(0)[0];
    T reg118=elem.pos(1)[0]*var_inter[0]; reg43=reg29+reg43; reg29=reg111*reg95; T reg119=reg100*elem.pos(1)[0]; T reg120=reg100*elem.pos(0)[0];
    T reg121=reg78*reg82; T reg122=reg88*reg27; reg1=reg62*reg1; reg35=reg83+reg35; reg73=reg80+reg73;
    reg80=reg111*reg76; reg83=reg57*reg95; T reg123=reg38*reg27; T reg124=reg78*reg5; reg97=reg34+reg97;
    reg34=reg22*reg82; reg59=reg21*reg59; reg68=reg68/reg71; reg90=reg24*reg90; reg30=reg85+reg30;
    reg24=reg44*reg53; reg85=reg57*reg31; reg96=reg42+reg96; reg74=reg28+reg74; reg17=reg65*reg17;
    reg57=reg57*reg76; reg86=reg32+reg86; reg28=reg5*reg88; reg107=reg56+reg107; reg58=reg55+reg58;
    reg111=reg111*reg31; reg84=reg104+reg84; reg121=reg122+reg121; reg57=reg86+reg57; reg61=reg61*reg68;
    reg32=reg89*reg68; reg42=reg8+reg115; reg55=elem.pos(2)[1]*var_inter[1]; reg17=reg74+reg17; reg56=elem.pos(2)[1]*var_inter[0];
    reg5=reg66*reg5; reg24=reg30+reg24; reg30=elem.pos(2)[0]*var_inter[1]; reg119=reg119-reg120; reg74=elem.pos(2)[0]*var_inter[0];
    reg86=reg117+reg118; reg103=reg101+reg103; reg116=reg116-reg109; reg101=reg44*reg68; reg87=reg102+reg87;
    reg80=reg73+reg80; reg98=reg99+reg98; reg34=reg123+reg34; reg111=reg43+reg111; reg43=reg1*reg94;
    reg85=reg96+reg85; reg124=reg97+reg124; reg83=reg35+reg83; reg110=reg79-reg110; reg0=reg113*reg0;
    reg92=reg93-reg92; reg114=reg15*reg114; reg29=reg58+reg29; reg28=reg107+reg28; reg105=reg106+reg105;
    reg108=reg112+reg108; reg45=reg21*reg45; reg9=reg26*reg9; reg21=reg59+reg90; reg26=reg1*reg63;
    reg32=reg34+reg32; reg45=reg9-reg45; reg9=reg85*reg80; reg34=reg83*reg80; reg21=reg15*reg21;
    reg15=elem.pos(3)[0]*var_inter[1]; reg30=reg119+reg30; reg43=reg124+reg43; reg35=reg0*reg31; reg101=reg103+reg101;
    reg5=reg17+reg5; reg1=reg1*reg77; reg61=reg121+reg61; reg17=reg51*reg105; reg58=reg48*reg87;
    reg73=elem.pos(3)[1]*var_inter[1]; reg116=reg55+reg116; reg55=reg51*reg108; reg79=reg48*reg98; reg93=reg6*reg72;
    reg96=reg72*reg40; reg97=reg7*reg40; reg99=reg48*reg84; reg102=reg51*reg24; reg114=reg92-reg114;
    reg74=reg74-reg86; reg92=elem.pos(3)[0]*reg10; reg26=reg28+reg26; reg28=reg0*reg95; reg110=reg110/reg71;
    reg56=reg56-reg42; reg103=elem.pos(3)[1]*reg10; reg104=reg29*reg57; reg106=reg111*reg57; reg107=reg83*reg111;
    reg34=reg104-reg34; reg104=reg29*reg85; reg112=reg46*reg108; reg119=reg23*reg98; reg121=reg96*reg32;
    reg79=reg55+reg79; reg55=reg96*reg101; reg99=reg102+reg99; reg102=reg46*reg105; reg122=reg23*reg84;
    reg123=reg46*reg24; reg124=reg66*reg110; reg9=reg106-reg9; reg106=reg23*reg87; reg24=reg38*reg24;
    reg84=reg22*reg84; reg108=reg38*reg108; reg98=reg22*reg98; reg87=reg22*reg87; T reg125=reg54*reg110;
    T reg126=reg52*reg110; reg105=reg38*reg105; reg58=reg17+reg58; reg21=reg45-reg21; reg17=reg14*reg72;
    reg28=reg26+reg28; reg26=reg13*reg40; reg116=reg116-reg73; reg45=reg96*reg61; reg92=reg74+reg92;
    reg103=reg56+reg103; reg93=reg97+reg93; reg56=reg7*reg72; reg74=reg6*reg40; reg114=reg114/reg71;
    reg0=reg0*reg76; reg1=reg5+reg1; reg30=reg30-reg15; reg35=reg43+reg35; reg5=reg89*reg32;
    reg43=reg66*reg125; reg121=reg79+reg121; reg17=reg26+reg17; reg98=reg108+reg98; reg26=reg89*reg101;
    reg79=reg13*reg72; reg84=reg24+reg84; reg45=reg58+reg45; reg24=reg44*reg61; reg58=reg66*reg126;
    reg55=reg99+reg55; reg97=reg14*reg40; reg106=reg102+reg106; reg99=reg35*reg34; reg32=reg44*reg32;
    reg102=reg66*reg124; reg87=reg105+reg87; reg61=reg89*reg61; reg107=reg104-reg107; reg0=reg1+reg0;
    reg1=reg93*reg114; reg119=reg112+reg119; reg71=reg21/reg71; reg56=reg56*reg114; reg74=reg74*reg114;
    reg21=reg116*reg92; reg104=reg28*reg9; reg105=reg30*reg103; reg122=reg123+reg122; reg101=reg44*reg101;
    reg108=reg29*reg0; reg97=reg97*reg71; reg24=reg106+reg24; reg106=reg19*reg124; reg61=reg87+reg61;
    reg124=reg16*reg124; reg87=reg28*reg80; reg80=reg35*reg80; reg21=reg105-reg21; reg26=reg84+reg26;
    reg84=reg16*reg126; reg105=reg16*reg125; reg112=reg111*reg0; reg5=reg98+reg5; reg98=reg93*reg1;
    reg12=reg40*reg12; reg123=reg0*reg107; T reg127=reg93*reg56; reg43=reg121+reg43; reg111=reg28*reg111;
    reg29=reg29*reg35; reg99=reg104-reg99; reg126=reg19*reg126; reg104=reg93*reg74; reg79=reg79*reg71;
    reg11=reg72*reg11; reg102=reg45+reg102; reg45=reg17*reg71; reg125=reg19*reg125; reg32=reg119+reg32;
    reg101=reg122+reg101; reg58=reg55+reg58; reg123=reg99+reg123; reg64=reg64*reg40; reg55=reg35*reg57;
    reg75=reg75*reg72; reg57=reg28*reg57; reg99=reg11*reg56; reg105=reg5+reg105; reg5=reg85*reg0;
    reg119=reg11*reg1; reg124=reg61+reg124; reg116=reg116/reg21; reg92=reg92/reg21; reg30=reg30/reg21;
    reg103=reg103/reg21; reg126=reg101+reg126; reg61=reg12*reg74; reg35=reg83*reg35; reg29=reg111-reg29;
    reg101=reg17*reg79; reg127=reg43+reg127; reg85=reg28*reg85; reg28=reg17*reg97; reg104=reg58+reg104;
    reg125=reg32+reg125; reg56=reg12*reg56; reg108=reg87-reg108; reg106=reg24+reg106; reg1=reg12*reg1;
    reg98=reg102+reg98; reg0=reg83*reg0; reg24=reg17*reg45; reg112=reg80-reg112; reg74=reg11*reg74;
    reg84=reg26+reg84; reg29=reg29/reg123; reg112=reg112/reg123; reg26=reg100*reg92; reg32=reg10*reg30;
    reg5=reg55-reg5; reg35=reg85-reg35; reg0=reg57-reg0; reg108=reg108/reg123; reg43=1-(*f.m).resolution;
    reg55=reg100*reg103; reg61=reg126+reg61; reg57=reg64*reg97; reg58=reg10*reg116; reg101=reg127+reg101;
    reg28=reg104+reg28; reg56=reg125+reg56; reg80=reg64*reg79; reg1=reg106+reg1; reg83=reg64*reg45;
    reg24=reg98+reg24; reg74=reg84+reg74; reg97=reg75*reg97; reg99=reg105+reg99; reg79=reg75*reg79;
    reg84=var_inter[1]*reg92; reg85=var_inter[1]*reg103; reg87=var_inter[0]*reg116; reg98=var_inter[0]*reg30; reg45=reg75*reg45;
    reg119=reg124+reg119; reg102=(*f.m).resolution*reg29; reg104=(*f.m).resolution*reg108; reg57=reg61+reg57; reg61=reg26-reg32;
    reg105=reg58-reg55; reg97=reg74+reg97; reg74=reg55+reg87; reg80=reg56+reg80; reg56=reg43*reg24;
    reg106=reg26+reg98; reg111=reg28*reg43; reg83=reg1+reg83; reg1=(*f.m).resolution*reg112; reg121=reg101*reg43;
    reg79=reg99+reg79; reg35=reg35/reg123; reg45=reg119+reg45; reg99=reg85-reg87; reg119=reg98-reg84;
    reg107=reg107/reg123; reg122=reg58+reg85; reg124=reg84+reg32; reg9=reg9/reg123; reg0=reg0/reg123;
    reg34=reg34/reg123; reg123=reg5/reg123; reg5=0.5*reg61; reg125=0.5*reg105; reg126=0.5*reg106;
    reg127=0.5*reg74; reg1=reg111+reg1; reg104=reg121-reg104; reg102=reg56+reg102; reg56=0.5*reg99;
    reg111=0.5*reg119; reg121=(*f.m).resolution*reg34; T reg128=0.5*reg124; T reg129=0.5*reg122; T reg130=(*f.m).resolution*reg35;
    T reg131=(*f.m).resolution*reg0; T reg132=(*f.m).resolution*reg123; T reg133=reg45*reg43; T reg134=reg43*reg79; T reg135=reg43*reg97;
    T reg136=reg43*reg83; T reg137=reg43*reg80; T reg138=reg43*reg57; T reg139=(*f.m).resolution*reg107; T reg140=(*f.m).resolution*reg9;
    T reg141=reg129*reg102; T reg142=reg102*reg5; T reg143=reg1*reg105; T reg144=reg124*reg104; T reg145=reg104*reg61;
    T reg146=reg102*reg125; T reg147=reg56*reg102; T reg148=reg119*reg104; T reg149=reg111*reg102; T reg150=reg99*reg1;
    T reg151=reg102*reg127; T reg152=reg104*reg106; reg136=reg139+reg136; reg132=reg135-reg132; reg135=reg1*reg74;
    reg139=reg102*reg126; T reg153=reg128*reg102; reg131=reg134+reg131; reg130=reg133-reg130; reg121=reg137-reg121;
    reg138=reg140+reg138; reg133=reg122*reg1; reg147=reg148+reg147; reg134=reg129*reg130; reg137=reg111*reg130;
    reg140=reg56*reg136; reg148=reg119*reg121; T reg154=reg99*reg132; T reg155=reg56*reg130; T reg156=reg119*reg131;
    reg149=reg150+reg149; reg150=reg128*reg136; T reg157=reg122*reg138; T reg158=reg111*reg136; T reg159=reg99*reg138;
    T reg160=reg128*reg130; T reg161=reg122*reg132; reg151=reg151-reg152; T reg162=reg124*reg121; T reg163=reg129*reg136;
    reg146=reg145+reg146; reg153=reg153-reg133; reg135=reg135-reg139; reg145=reg138*reg74; T reg164=reg136*reg125;
    T reg165=reg121*reg61; T reg166=reg136*reg126; reg144=reg144-reg141; T reg167=reg130*reg125; T reg168=reg131*reg61;
    reg142=reg143+reg142; reg143=reg121*reg106; T reg169=reg136*reg127; T reg170=reg130*reg126; T reg171=reg136*reg5;
    T reg172=reg138*reg105; T reg173=reg132*reg74; T reg174=reg124*reg131; T reg175=reg130*reg127; T reg176=reg131*reg106;
    reg135=2*reg135; reg174=reg174-reg134; reg160=reg160-reg161; reg156=reg155+reg156; reg154=reg137+reg154;
    reg175=reg175-reg176; reg173=reg173-reg170; reg168=reg167+reg168; reg150=reg150-reg157; reg144=2*reg144;
    reg162=reg162-reg163; reg153=2*reg153; reg140=reg148+reg140; reg151=2*reg151; reg169=reg169-reg143;
    reg171=reg172+reg171; reg147=2*reg147; reg142=2*reg142; reg164=reg165+reg164; reg146=2*reg146;
    reg158=reg159+reg158; reg145=reg145-reg166; reg149=2*reg149; reg137=reg151*reg126; reg148=reg169*reg74;
    reg155=reg158*reg74; reg159=reg149*reg126; reg165=reg140*reg74; reg167=reg147*reg126; reg172=reg144*reg129;
    T reg177=reg174*reg124; T reg178=reg150*reg74; T reg179=reg153*reg126; T reg180=reg162*reg74; T reg181=reg144*reg126;
    T reg182=reg151*reg127; T reg183=reg175*reg106; T reg184=reg149*reg127; T reg185=reg154*reg106; T reg186=reg147*reg127;
    T reg187=reg156*reg106; T reg188=reg153*reg127; T reg189=reg160*reg106; T reg190=reg144*reg127; T reg191=reg174*reg106;
    T reg192=reg158*reg99; T reg193=reg145*reg105; T reg194=reg150*reg105; T reg195=reg153*reg5; T reg196=reg146*reg5;
    T reg197=reg164*reg105; T reg198=reg162*reg105; T reg199=reg144*reg5; T reg200=reg142*reg5; T reg201=reg171*reg105;
    T reg202=reg168*reg61; T reg203=reg146*reg125; T reg204=reg173*reg61; T reg205=reg135*reg125; T reg206=reg150*reg99;
    T reg207=reg153*reg111; T reg208=reg162*reg99; T reg209=reg144*reg111; T reg210=reg156*reg119; T reg211=reg147*reg56;
    T reg212=reg160*reg119; T reg213=reg153*reg56; T reg214=reg174*reg119; T reg215=reg144*reg56; T reg216=reg153*reg128;
    reg150=reg150*reg122; T reg217=reg144*reg128; reg162=reg162*reg122; T reg218=reg147*reg5; T reg219=reg140*reg105;
    T reg220=reg147*reg111; T reg221=reg175*reg61; T reg222=reg151*reg125; T reg223=reg135*reg126; T reg224=reg145*reg74;
    reg144=reg144*reg125; reg174=reg174*reg61; T reg225=reg154*reg61; T reg226=reg149*reg125; T reg227=reg169*reg105;
    T reg228=reg149*reg5; T reg229=reg156*reg61; T reg230=reg147*reg125; T reg231=reg153*reg125; T reg232=reg160*reg61;
    T reg233=reg158*reg105; T reg234=reg151*reg5; T reg235=reg140*reg99; T reg236=reg149*reg111; T reg237=reg135*reg5;
    reg203=reg202+reg203; reg159=reg155-reg159; reg189=reg188-reg189; reg226=reg225+reg226; reg209=reg208+reg209;
    reg218=reg219+reg218; reg137=reg148-reg137; reg193=reg237+reg193; reg205=reg204+reg205; reg231=reg232+reg231;
    reg228=reg233+reg228; reg207=reg206+reg207; reg230=reg229+reg230; reg162=reg217-reg162; reg172=reg177-reg172;
    reg236=reg192+reg236; reg191=reg190-reg191; reg150=reg216-reg150; reg185=reg184-reg185; reg220=reg235+reg220;
    reg223=reg224-reg223; reg183=reg182-reg183; reg213=reg212+reg213; reg195=reg194+reg195; reg187=reg186-reg187;
    reg181=reg180-reg181; reg222=reg221+reg222; reg200=reg201+reg200; reg196=reg197+reg196; reg199=reg198+reg199;
    reg179=reg178-reg179; reg211=reg210+reg211; reg234=reg227+reg234; reg144=reg174+reg144; reg167=reg165-reg167;
    reg215=reg214+reg215; reg211=reg211*reg21; reg209=reg209*reg21; reg234=reg234*reg21; reg172=reg21*reg172;
    reg207=reg207*reg21; reg231=reg231*reg21; reg218=reg218*reg21; reg203=reg203*reg21; reg193=reg193*reg21;
    reg196=reg21*reg196; reg144=reg144*reg21; reg199=reg199*reg21; reg200=reg21*reg200; reg195=reg195*reg21;
    reg220=reg220*reg21; reg159=reg159*reg21; reg228=reg228*reg21; reg167=reg167*reg21; reg205=reg205*reg21;
    reg137=reg137*reg21; reg226=reg226*reg21; reg185=reg185*reg21; reg187=reg187*reg21; reg223=reg223*reg21;
    reg183=reg183*reg21; reg162=reg162*reg21; reg189=reg189*reg21; reg150=reg150*reg21; reg215=reg215*reg21;
    reg236=reg236*reg21; reg230=reg230*reg21; reg191=reg191*reg21; reg179=reg179*reg21; reg213=reg213*reg21;
    reg181=reg181*reg21; reg222=reg222*reg21; T tmp_2_7=ponderation*reg181; T tmp_0_0=ponderation*reg200; T tmp_3_3=ponderation*reg183;
    T tmp_3_4=ponderation*reg185; T tmp_0_6=ponderation*reg195; T tmp_1_2=ponderation*reg205; T tmp_4_5=ponderation*reg220; T tmp_3_5=ponderation*reg187;
    T tmp_0_1=ponderation*reg196; T tmp_4_4=ponderation*reg236; T tmp_3_6=ponderation*reg189; T tmp_0_5=ponderation*reg218; T tmp_7_7=ponderation*reg172;
    T tmp_3_7=ponderation*reg191; T tmp_0_3=ponderation*reg234; T tmp_5_5=ponderation*reg211; T tmp_5_6=ponderation*reg213; T tmp_1_5=ponderation*reg230;
    T tmp_4_7=ponderation*reg209; T tmp_5_7=ponderation*reg215; T tmp_6_6=ponderation*reg150; T tmp_4_6=ponderation*reg207; T tmp_6_7=ponderation*reg162;
    T tmp_2_2=ponderation*reg223; T tmp_1_4=ponderation*reg226; T tmp_1_6=ponderation*reg231; T tmp_2_3=ponderation*reg137; T tmp_0_2=ponderation*reg193;
    T tmp_1_1=ponderation*reg203; T tmp_2_4=ponderation*reg159; T tmp_2_5=ponderation*reg167; T tmp_0_4=ponderation*reg228; T tmp_1_7=ponderation*reg144;
    T tmp_0_7=ponderation*reg199; T tmp_1_3=ponderation*reg222; T tmp_2_6=ponderation*reg179;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; T reg2=2*(*f.m).shear_modulus_12; reg1=1.0/reg1; reg0=1.0/reg0;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=pow((*f.m).v1[0],2); T reg5=pow((*f.m).v1[1],2); T reg6=pow((*f.m).v2[0],2);
    T reg7=pow((*f.m).v2[1],2); T reg8=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg11=reg2*reg3;
    reg5=reg4+reg5; reg4=reg8*reg11; T reg12=pow((*f.m).v1[2],2); T reg13=reg10*reg11; T reg14=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1;
    T reg15=1.0/(*f.m).elastic_modulus_2; reg7=reg6+reg7; reg6=reg9*reg11; T reg16=pow((*f.m).v2[2],2); T reg17=reg15*reg6;
    reg12=reg5+reg12; reg5=reg10*reg13; T reg18=reg14*reg6; T reg19=reg10*reg4; reg16=reg7+reg16;
    reg7=reg15*reg4; T reg20=reg14*reg13; reg19=reg18+reg19; T reg21=1.0/(*f.m).elastic_modulus_1; reg12=pow(reg12,0.5);
    reg5=reg17-reg5; reg16=pow(reg16,0.5); reg17=reg20+reg7; T reg22=reg21*reg5; T reg23=(*f.m).v2[2]/reg16;
    T reg24=reg14*reg19; T reg25=(*f.m).v2[1]/reg16; T reg26=(*f.m).v1[1]/reg12; T reg27=(*f.m).v1[2]/reg12; T reg28=reg10*reg3;
    T reg29=reg14*reg11; T reg30=reg8*reg4; T reg31=reg27*reg25; reg6=reg21*reg6; T reg32=reg26*reg23;
    reg16=(*f.m).v2[0]/reg16; T reg33=reg8*reg3; T reg34=reg8*reg13; reg3=reg9*reg3; reg11=reg15*reg11;
    reg12=(*f.m).v1[0]/reg12; T reg35=reg8*reg17; reg24=reg22-reg24; reg22=reg2*reg0; T reg36=reg27*reg16;
    T reg37=reg12*reg23; T reg38=2*reg16; T reg39=reg32-reg31; T reg40=reg2*reg1; T reg41=reg8*reg22;
    T reg42=reg15*reg3; T reg43=2*reg12; reg35=reg24-reg35; reg24=reg9*reg22; T reg44=reg10*reg33;
    T reg45=reg8*reg11; T reg46=reg10*reg28; reg34=reg18+reg34; reg3=reg14*reg3; reg30=reg6-reg30;
    reg13=reg21*reg13; reg6=reg8*reg29; reg4=reg14*reg4; reg22=reg10*reg22; reg18=pow(reg25,2);
    reg5=reg5/reg35; T reg47=reg26*reg43; reg11=reg21*reg11; T reg48=reg8*reg40; reg4=reg13+reg4;
    T reg49=pow(reg16,2); reg33=reg15*reg33; T reg50=reg10*reg40; reg34=reg34/reg35; reg28=reg14*reg28;
    reg6=reg13+reg6; reg45=reg20+reg45; reg19=reg19/reg35; reg30=reg30/reg35; reg29=reg14*reg29;
    reg13=reg36-reg37; reg46=reg42-reg46; reg42=pow(reg12,2); T reg51=pow(reg26,2); T reg52=reg26*reg16;
    reg44=reg3+reg44; reg3=reg12*reg25; T reg53=2*reg39; T reg54=reg15*reg24; reg24=reg14*reg24;
    T reg55=reg10*reg22; T reg56=reg10*reg41; reg40=reg9*reg40; T reg57=reg38*reg25; T reg58=reg18*reg30;
    reg17=reg17/reg35; T reg59=reg42*reg5; T reg60=reg47*reg34; T reg61=reg49*reg19; T reg62=reg57*reg30;
    T reg63=reg51*reg34; reg46=reg21*reg46; reg45=reg45/reg35; reg44=reg14*reg44; T reg64=reg28+reg33;
    reg55=reg54-reg55; reg56=reg24+reg56; reg24=reg15*reg40; reg40=reg14*reg40; reg54=pow(reg23,2);
    T reg65=reg57*reg19; T reg66=reg47*reg5; T reg67=pow(reg39,2); T reg68=pow(reg13,2); reg41=reg15*reg41;
    reg22=reg14*reg22; T reg69=reg53*reg13; T reg70=pow(reg27,2); T reg71=reg18*reg19; T reg72=reg51*reg5;
    reg29=reg11-reg29; reg11=reg10*reg48; reg6=reg6/reg35; T reg73=reg3-reg52; T reg74=reg10*reg50;
    T reg75=reg49*reg30; T reg76=reg42*reg34; reg4=reg4/reg35; T reg77=reg67*reg4; reg61=reg59+reg61;
    reg59=reg67*reg17; reg75=reg76+reg75; reg58=reg63+reg58; reg63=pow(reg73,2); reg76=reg68*reg4;
    T reg78=reg70*reg34; T reg79=reg54*reg30; T reg80=reg57*reg6; T reg81=reg47*reg45; T reg82=reg18*reg6;
    T reg83=reg51*reg45; T reg84=reg49*reg6; T reg85=reg42*reg45; reg62=reg60+reg62; reg60=reg69*reg4;
    reg71=reg72+reg71; reg72=reg68*reg17; T reg86=reg70*reg5; T reg87=reg54*reg19; reg64=reg8*reg64;
    reg11=reg40+reg11; reg55=reg21*reg55; reg74=reg24-reg74; reg48=reg15*reg48; reg56=reg14*reg56;
    reg65=reg66+reg65; reg29=reg29/reg35; reg50=reg14*reg50; reg44=reg46-reg44; reg24=reg22+reg41;
    reg40=reg69*reg17; reg40=reg65+reg40; reg72=reg71+reg72; reg46=reg27*reg43; reg65=2*reg26;
    reg66=reg50+reg48; reg59=reg61+reg59; reg64=reg44-reg64; reg44=reg38*reg23; reg69=reg69*reg29;
    reg80=reg81+reg80; reg84=reg85+reg84; reg11=reg14*reg11; reg77=reg75+reg77; reg61=reg67*reg29;
    reg60=reg62+reg60; reg62=reg12*reg26; reg71=reg63*reg17; reg87=reg86+reg87; reg75=reg16*reg25;
    reg82=reg83+reg82; reg56=reg55-reg56; reg55=reg63*reg4; reg79=reg78+reg79; reg78=reg68*reg29;
    reg74=reg21*reg74; reg76=reg58+reg76; reg58=reg70*reg45; reg81=reg54*reg6; reg24=reg8*reg24;
    reg83=2*reg25; reg85=reg18*reg76; reg24=reg56-reg24; reg56=reg51*reg72; reg86=reg16*reg43;
    reg69=reg80+reg69; reg80=reg18*reg77; T reg88=reg51*reg59; T reg89=reg63*reg29; reg81=reg58+reg81;
    reg78=reg82+reg78; reg55=reg79+reg55; reg61=reg84+reg61; reg58=reg75*reg76; reg79=reg62*reg72;
    reg82=reg46*reg34; reg84=reg44*reg30; reg64=reg64/reg35; T reg90=reg42*reg59; T reg91=reg12*reg13;
    T reg92=reg26*reg39; T reg93=reg49*reg77; reg72=reg42*reg72; reg76=reg49*reg76; reg3=reg52+reg3;
    reg52=reg26*reg25; T reg94=reg12*reg16; T reg95=2*reg13; reg53=reg53*reg73; T reg96=reg62*reg40;
    reg71=reg87+reg71; reg87=reg83*reg23; T reg97=reg46*reg5; T reg98=reg44*reg19; T reg99=reg49*reg60;
    T reg100=reg42*reg40; T reg101=reg25*reg65; T reg102=reg18*reg60; T reg103=2*reg27; reg40=reg51*reg40;
    T reg104=reg13*reg39; reg11=reg74-reg11; reg66=reg8*reg66; reg74=reg27*reg65; reg60=reg75*reg60;
    T reg105=reg52*reg64; T reg106=reg68*reg69; reg99=reg100+reg99; reg102=reg40+reg102; reg40=reg67*reg69;
    reg100=reg49*reg55; T reg107=reg101*reg15; T reg108=reg86*reg14; T reg109=reg18*reg15; T reg110=reg3*reg64;
    T reg111=reg42*reg71; reg24=reg24/reg35; T reg112=reg49*reg14; T reg113=reg67*reg78; T reg114=reg46*reg45;
    reg66=reg11-reg66; reg76=reg72+reg76; reg11=reg101*reg14; reg72=reg86*reg21; T reg115=reg94*reg64;
    T reg116=reg67*reg61; reg93=reg90+reg93; reg89=reg81+reg89; reg59=reg62*reg59; reg77=reg75*reg77;
    reg81=reg44*reg6; reg90=reg49*reg21; T reg117=reg18*reg14; reg80=reg88+reg80; reg88=reg68*reg61;
    T reg118=reg51*reg71; reg95=reg95*reg73; T reg119=reg27*reg23; T reg120=reg18*reg55; reg19=reg87*reg19;
    reg5=reg74*reg5; reg85=reg56+reg85; reg56=reg68*reg78; T reg121=reg12*reg39; T reg122=reg26*reg13;
    reg91=reg92+reg91; reg103=reg23*reg103; reg92=reg16*reg13; T reg123=reg25*reg39; reg69=reg104*reg69;
    reg60=reg96+reg60; reg30=reg87*reg30; reg98=reg97+reg98; reg96=reg53*reg4; reg97=reg53*reg17;
    reg34=reg74*reg34; reg84=reg82+reg84; reg78=reg104*reg78; reg58=reg79+reg58; reg79=reg27*reg73;
    reg43=reg39*reg43; reg82=reg101*reg105; reg56=reg85+reg56; reg116=reg93+reg116; reg97=reg98+reg97;
    reg85=reg67*reg89; reg93=reg86*reg115; reg6=reg87*reg6; reg78=reg58+reg78; reg58=reg101*reg10;
    reg45=reg74*reg45; reg120=reg118+reg120; reg30=reg34+reg30; reg34=reg16*reg39; reg98=reg25*reg13;
    reg118=reg119*reg64; reg4=reg95*reg4; reg92=reg123+reg92; reg69=reg60+reg69; reg60=reg68*reg89;
    reg121=reg121*reg24; reg122=reg122*reg24; reg123=reg3*reg110; reg35=reg66/reg35; reg66=reg91*reg24;
    reg65=reg13*reg65; reg106=reg102+reg106; reg102=reg101*reg110; reg96=reg84+reg96; reg112=reg109-reg112;
    reg84=reg54*reg10; reg100=reg111+reg100; reg77=reg59+reg77; reg110=reg86*reg110; reg11=reg72-reg11;
    reg59=reg103*reg8; reg72=reg101*reg115; reg88=reg80+reg88; reg40=reg99+reg40; reg117=reg90-reg117;
    reg80=reg54*reg8; reg71=reg62*reg71; reg61=reg104*reg61; reg55=reg75*reg55; reg90=reg86*reg8;
    reg99=reg3*reg105; reg109=reg18*reg10; reg111=reg49*reg8; reg53=reg53*reg29; reg81=reg114+reg81;
    reg19=reg5+reg19; reg17=reg95*reg17; reg113=reg76+reg113; reg105=reg86*reg105; reg108=reg107-reg108;
    reg5=reg103*reg10; reg76=reg43*reg122; reg107=(*f.m).alpha_1*reg51; reg114=(*f.m).alpha_2*reg18; reg36=reg37+reg36;
    reg37=reg49*reg96; T reg124=reg43*reg121; T reg125=reg42*reg97; reg93=reg116+reg93; reg116=reg86*reg118;
    reg21=reg42*reg21; T reg126=reg92*reg35; reg105=reg113+reg105; reg113=reg27*reg39; T reg127=(*f.m).alpha_1*reg42;
    reg98=reg98*reg35; reg110=reg40+reg110; reg85=reg100+reg85; reg79=reg79*reg24; reg40=reg43*reg66;
    reg34=reg34*reg35; reg58=reg90+reg58; reg103=reg103*reg9; reg109=reg111+reg109; reg90=reg54*reg9;
    reg53=reg81+reg53; reg17=reg19+reg17; reg5=reg108-reg5; reg19=reg65*reg121; reg72=reg88+reg72;
    reg84=reg112-reg84; reg89=reg104*reg89; reg15=reg51*reg15; reg81=reg42*reg14; reg59=reg11-reg59;
    reg80=reg117-reg80; reg55=reg71+reg55; reg14=reg51*reg14; reg11=reg91*reg122; reg99=reg78+reg99;
    reg61=reg77+reg61; reg115=reg3*reg115; reg4=reg30+reg4; reg123=reg69+reg123; reg30=reg91*reg66;
    reg83=reg83*reg13; reg38=reg38*reg39; reg69=reg101*reg118; reg60=reg120+reg60; reg71=reg23*reg73;
    reg77=(*f.m).alpha_2*reg49; reg102=reg106+reg102; reg66=reg65*reg66; reg78=reg51*reg97; reg88=reg18*reg96;
    reg122=reg65*reg122; reg82=reg56+reg82; reg29=reg95*reg29; reg56=reg12*reg73; reg6=reg45+reg6;
    reg45=reg36*reg64; reg29=reg6+reg29; reg116=reg85+reg116; reg96=reg75*reg96; reg97=reg62*reg97;
    reg6=reg92*reg126; reg71=reg71*reg35; reg85=reg38*reg98; reg76=reg105+reg76; reg30=reg123+reg30;
    reg95=reg38*reg34; reg124=reg93+reg124; reg56=reg113+reg56; reg32=reg31+reg32; reg31=reg5*reg18;
    reg93=reg59*reg49; reg100=(*f.m).alpha_2*reg54; reg105=(*f.m).alpha_1*reg70; reg106=reg84*reg18; reg108=(*f.m).alpha_3*reg68;
    reg114=reg107+reg114; reg107=reg80*reg49; reg111=(*f.m).alpha_3*reg67; reg77=reg127+reg77; reg112=reg5*reg51;
    reg58=reg103-reg58; reg109=reg90-reg109; reg90=reg59*reg42; reg89=reg55+reg89; reg55=reg70*reg8;
    reg14=reg21-reg14; reg118=reg3*reg118; reg21=reg70*reg10; reg81=reg15-reg81; reg8=reg42*reg8;
    reg10=reg51*reg10; reg15=reg83*reg98; reg5=reg52*reg5; reg59=reg94*reg59; reg98=reg92*reg98;
    reg103=reg80*reg42; reg11=reg99+reg11; reg99=reg84*reg51; reg121=reg91*reg121; reg115=reg61+reg115;
    reg61=reg43*reg79; reg84=reg52*reg84; reg40=reg110+reg40; reg110=reg38*reg126; reg80=reg94*reg80;
    reg19=reg72+reg19; reg37=reg125+reg37; reg72=reg67*reg53; reg113=reg83*reg34; reg122=reg82+reg122;
    reg82=reg65*reg79; reg69=reg60+reg69; reg66=reg102+reg66; reg60=reg42*reg17; reg102=reg49*reg4;
    reg117=reg27*reg13; reg120=reg26*reg73; reg123=reg18*reg4; reg126=reg83*reg126; reg125=reg23*reg39;
    reg88=reg78+reg88; reg78=reg68*reg53; reg127=reg51*reg17; T reg128=reg16*reg73; T reg129=reg54*reg109;
    reg15=reg122+reg15; reg79=reg91*reg79; reg118=reg89+reg118; reg31=reg93+reg31; reg89=reg54*reg58;
    reg6=reg30+reg6; reg84=reg80+reg84; reg30=reg119*reg109; reg80=reg119*reg58; reg5=reg59+reg5;
    reg120=reg117+reg120; reg113=reg19+reg113; reg106=reg107+reg106; reg12=reg12*reg27; reg99=reg103+reg99;
    reg109=reg109*reg70; reg128=reg125+reg128; reg16=reg16*reg23; reg19=reg25*reg73; reg112=reg90+reg112;
    reg58=reg58*reg70; reg59=reg23*reg13; reg82=reg69+reg82; reg95=reg124+reg95; reg10=reg8+reg10;
    reg8=reg75*reg2; reg69=reg3*reg2; reg102=reg60+reg102; reg111=reg77+reg111; reg108=reg114+reg108;
    reg64=reg32*reg64; reg100=reg105+reg100; reg67=reg67*reg29; reg126=reg66+reg126; reg63=(*f.m).alpha_3*reg63;
    reg56=reg56*reg24; reg78=reg88+reg78; reg68=reg68*reg29; reg60=reg101*reg45; reg66=(*f.m).alpha_1*reg62;
    reg77=(*f.m).alpha_2*reg75; reg123=reg127+reg123; reg96=reg97+reg96; reg55=reg14-reg55; reg53=reg104*reg53;
    reg21=reg81-reg21; reg61=reg116+reg61; reg17=reg62*reg17; reg4=reg75*reg4; reg14=reg38*reg71;
    reg110=reg40+reg110; reg98=reg11+reg98; reg85=reg76+reg85; reg72=reg37+reg72; reg11=reg86*reg45;
    reg9=reg70*reg9; reg37=reg83*reg71; reg34=reg92*reg34; reg121=reg115+reg121; reg19=reg59+reg19;
    reg23=reg25*reg23; reg37=reg82+reg37; reg27=reg26*reg27; reg39=reg73*reg39; reg30=reg84+reg30;
    reg25=reg108*reg15; reg26=reg111*reg113; reg60=reg78+reg60; reg53=reg96+reg53; reg79=reg118+reg79;
    reg71=reg92*reg71; reg45=reg3*reg45; reg4=reg17+reg4; reg29=reg104*reg29; reg17=reg16*reg1;
    reg40=reg36*reg1; reg59=reg55*reg42; reg75=reg21*reg51; reg34=reg121+reg34; reg109=reg99+reg109;
    reg76=reg8*reg47; reg101=reg101*reg64; reg58=reg112+reg58; reg78=reg69*reg47; reg81=reg55*reg49;
    reg82=reg21*reg18; reg129=reg106+reg129; reg84=reg57*reg8; reg68=reg123+reg68; reg89=reg31+reg89;
    reg31=reg57*reg69; reg88=reg65*reg56; reg90=reg111*reg95; reg8=reg3*reg8; reg93=reg108*reg85;
    reg96=reg98*reg110; reg97=reg126*reg98; reg99=reg6*reg85; reg103=reg6*reg15; reg80=reg5+reg80;
    reg69=reg3*reg69; reg86=reg86*reg64; reg10=reg9-reg10; reg24=reg120*reg24; reg63=reg100+reg63;
    reg77=reg66+reg77; reg128=reg128*reg35; reg14=reg61+reg14; reg11=reg72+reg11; reg104=(*f.m).alpha_3*reg104;
    reg67=reg102+reg67; reg5=reg43*reg56; reg16=(*f.m).alpha_2*reg16; reg9=(*f.m).alpha_1*reg12; reg61=reg38*reg128;
    reg97=reg103-reg97; reg5=reg11+reg5; reg76=reg109+reg76; reg11=reg17*reg46; reg96=reg99-reg96;
    reg66=reg126*reg85; reg72=reg110*reg15; reg78=reg58+reg78; reg58=reg40*reg46; reg82=reg81+reg82;
    reg54=reg54*reg10; reg84=reg129+reg84; reg81=reg17*reg44; reg17=reg17*reg36; reg31=reg89+reg31;
    reg89=reg40*reg44; reg8=reg30+reg8; reg13=reg73*reg13; reg35=reg19*reg35; reg19=reg83*reg128;
    reg88=reg60+reg88; reg2=reg62*reg2; reg101=reg68+reg101; reg65=reg65*reg24; reg30=(*f.m).alpha_1*reg27;
    reg60=(*f.m).alpha_2*reg23; reg64=reg3*reg64; reg29=reg4+reg29; reg56=reg91*reg56; reg45=reg53+reg45;
    reg104=reg77+reg104; reg40=reg40*reg36; reg69=reg80+reg69; reg71=reg79+reg71; reg16=reg9+reg16;
    reg39=(*f.m).alpha_3*reg39; reg23=reg23*reg0; reg4=reg32*reg0; reg75=reg59+reg75; reg70=reg10*reg70;
    reg21=reg52*reg21; reg43=reg43*reg24; reg9=reg34*reg111; reg55=reg94*reg55; reg25=reg26+reg25;
    reg26=reg98*reg108; reg14=reg63*reg14; reg93=reg90+reg93; reg86=reg67+reg86; reg37=reg63*reg37;
    reg21=reg55+reg21; reg54=reg82+reg54; reg53=reg57*reg2; reg14=reg93+reg14; reg55=reg104*reg110;
    reg10=reg119*reg10; reg1=reg12*reg1; reg17=reg8+reg17; reg19=reg88+reg19; reg81=reg84+reg81;
    reg8=reg23*reg87; reg12=reg4*reg74; reg58=reg78+reg58; reg59=reg23*reg32; reg72=reg66-reg72;
    reg65=reg101+reg65; reg83=reg83*reg35; reg23=reg23*reg74; reg11=reg76+reg11; reg39=reg16+reg39;
    reg16=reg96*reg113; reg38=reg38*reg35; reg62=reg97*reg95; reg43=reg86+reg43; reg70=reg75+reg70;
    reg66=reg2*reg47; reg67=reg126*reg104; reg40=reg69+reg40; reg68=reg4*reg87; reg4=reg4*reg32;
    reg64=reg29+reg64; reg37=reg25+reg37; reg89=reg31+reg89; reg26=reg9+reg26; reg63=reg71*reg63;
    reg24=reg91*reg24; reg13=(*f.m).alpha_3*reg13; reg60=reg30+reg60; reg56=reg45+reg56; reg128=reg92*reg128;
    reg61=reg5+reg61; reg16=reg62-reg16; reg23=reg11+reg23; reg24=reg64+reg24; reg5=reg34*reg72;
    reg83=reg65+reg83; reg4=reg40+reg4; reg38=reg43+reg38; reg9=reg104*reg6; reg63=reg26+reg63;
    reg46=reg1*reg46; reg66=reg70+reg66; reg128=reg56+reg128; reg8=reg81+reg8; reg10=reg21+reg10;
    reg61=reg39*reg61; reg55=reg14+reg55; reg2=reg3*reg2; reg11=reg6*reg95; reg19=reg19*reg39;
    reg14=reg126*reg34; reg67=reg37+reg67; reg59=reg17+reg59; reg0=reg27*reg0; reg44=reg1*reg44;
    reg53=reg54+reg53; reg68=reg89+reg68; reg17=reg34*reg110; reg21=reg6*reg113; reg12=reg58+reg12;
    reg13=reg60+reg13; reg35=reg92*reg35; reg128=reg39*reg128; reg9=reg63+reg9; reg25=reg12*reg59;
    reg35=reg24+reg35; reg24=reg23*reg4; reg26=reg110*reg113; reg27=reg8*reg4; reg29=reg126*reg95;
    reg30=reg68*reg59; reg31=reg34*reg85; reg61=reg55+reg61; reg17=reg11-reg17; reg11=reg98*reg95;
    reg19=reg67+reg19; reg83=reg83*reg13; reg46=reg66+reg46; reg74=reg0*reg74; reg5=reg16+reg5;
    reg2=reg10+reg2; reg36=reg1*reg36; reg38=reg13*reg38; reg44=reg53+reg44; reg87=reg0*reg87;
    reg1=reg98*reg113; reg14=reg21-reg14; reg10=reg34*reg15; reg16=1-var_inter[1]; reg38=reg61+reg38;
    reg21=1-var_inter[0]; reg37=reg12*reg8; reg39=reg23*reg68; reg25=reg24-reg25; reg83=reg19+reg83;
    reg32=reg0*reg32; reg36=reg2+reg36; reg96=reg96/reg5; reg30=reg27-reg30; reg17=reg17/reg5;
    reg74=reg46+reg74; reg10=reg1-reg10; reg31=reg11-reg31; reg14=reg14/reg5; reg35=reg13*reg35;
    reg128=reg9+reg128; reg87=reg44+reg87; reg0=reg95*reg15; reg97=reg97/reg5; reg26=reg29-reg26;
    reg1=reg85*reg113; reg1=reg0-reg1; reg10=reg10/reg5; reg0=reg21*elem.pos(0)[0]; reg32=reg36+reg32;
    reg17=reg17*reg83; reg31=reg31/reg5; reg72=reg72/reg5; reg2=reg74*reg30; reg26=reg26/reg5;
    reg97=reg97*reg38; reg14=reg14*reg38; reg96=reg96*reg83; reg37=reg39-reg37; reg35=reg128+reg35;
    reg9=reg16*elem.pos(0)[0]; reg11=reg16*elem.pos(1)[0]; reg13=reg87*reg25; reg19=reg16*elem.pos(1)[1]; reg24=reg16*elem.pos(0)[1];
    reg27=elem.pos(1)[0]*var_inter[0]; reg29=reg21*elem.pos(0)[1]; reg36=elem.pos(1)[1]*var_inter[0]; reg72=reg72*reg35; reg39=reg32*reg37;
    reg96=reg97-reg96; reg5=reg1/reg5; reg1=reg87*reg59; reg13=reg2-reg13; reg2=reg8*reg32;
    reg59=reg74*reg59; reg40=reg23*reg32; reg43=reg0+reg27; reg44=elem.pos(2)[0]*var_inter[0]; reg45=elem.pos(2)[1]*var_inter[0];
    reg26=reg26*reg35; reg11=reg11-reg9; reg46=elem.pos(2)[0]*var_inter[1]; reg53=elem.pos(2)[1]*var_inter[1]; reg83=reg31*reg83;
    reg14=reg17-reg14; reg38=reg10*reg38; reg19=reg19-reg24; reg10=reg36+reg29; reg17=reg12*reg32;
    reg2=reg1-reg2; reg19=reg53+reg19; reg1=1-(*f.m).resolution; reg31=reg74*reg4; reg32=reg68*reg32;
    reg53=elem.pos(3)[1]*var_inter[1]; reg54=elem.pos(3)[0]*var_inter[1]; reg96=reg72+reg96; reg26=reg14-reg26; reg35=reg5*reg35;
    reg83=reg38-reg83; reg4=reg87*reg4; reg23=reg23*reg87; reg39=reg13+reg39; reg40=reg59-reg40;
    reg5=elem.pos(3)[1]*reg21; reg45=reg45-reg10; reg13=elem.pos(3)[0]*reg21; reg44=reg44-reg43; reg46=reg11+reg46;
    reg8=reg74*reg8; reg26=reg26*reg1; reg46=reg46-reg54; reg68=reg74*reg68; reg2=reg2/reg39;
    reg5=reg45+reg5; reg11=reg111*(*f.m).resolution; reg32=reg4-reg32; reg4=reg108*(*f.m).resolution; reg13=reg44+reg13;
    reg83=reg35+reg83; reg17=reg31-reg17; reg96=reg96*reg1; reg23=reg8-reg23; reg19=reg19-reg53;
    reg40=reg40/reg39; reg87=reg12*reg87; reg23=reg23/reg39; reg17=reg17/reg39; reg87=reg68-reg87;
    reg25=reg25/reg39; reg34=reg34*reg1; reg98=reg98*reg1; reg8=reg104*(*f.m).resolution; reg12=(*f.m).resolution*reg2;
    reg14=(*f.m).resolution*reg40; reg32=reg32/reg39; reg26=reg4+reg26; reg30=reg30/reg39; reg83=reg83*reg1;
    reg96=reg11+reg96; reg4=reg46*reg5; reg11=reg19*reg13; reg11=reg4-reg11; reg26=reg26*(*f.m).deltaT;
    reg96=(*f.m).deltaT*reg96; reg83=reg8+reg83; reg14=reg98-reg14; reg12=reg34+reg12; reg4=(*f.m).resolution*reg23;
    reg8=(*f.m).resolution*reg17; reg31=(*f.m).resolution*reg32; reg6=reg1*reg6; reg15=reg1*reg15; reg113=reg1*reg113;
    reg85=reg1*reg85; reg95=reg1*reg95; reg34=(*f.m).resolution*reg30; reg35=(*f.m).resolution*reg25; reg87=reg87/reg39;
    reg39=reg37/reg39; reg37=reg96*reg12; reg38=reg26*reg14; reg83=reg83*(*f.m).deltaT; reg5=reg5/reg11;
    reg13=reg13/reg11; reg4=reg6+reg4; reg8=reg15+reg8; reg31=reg113-reg31; reg19=reg19/reg11;
    reg35=reg85-reg35; reg95=reg34+reg95; reg6=(*f.m).resolution*reg87; reg126=reg126*reg1; reg110=reg1*reg110;
    reg1=(*f.m).resolution*reg39; reg46=reg46/reg11; reg15=reg96*reg31; reg34=reg26*reg8; reg44=var_inter[0]*reg19;
    reg45=reg96*reg95; reg55=reg16*reg13; reg56=var_inter[0]*reg46; reg58=reg26*reg35; reg59=reg21*reg46;
    reg6=reg126-reg6; reg110=reg1+reg110; reg1=reg83*reg4; reg60=reg16*reg5; reg61=var_inter[1]*reg5;
    reg62=reg21*reg19; reg63=var_inter[1]*reg13; reg64=reg37+reg38; reg65=reg45+reg58; reg66=reg83*reg110;
    reg67=reg62+reg61; reg68=reg63+reg59; reg69=reg64+reg1; reg70=reg83*reg6; reg71=reg60+reg44;
    reg72=reg55+reg56; reg73=reg15+reg34; reg74=2*reg69; reg75=reg65+reg66; reg76=reg70+reg73;
    reg77=reg61-reg44; reg78=reg56-reg63; reg79=0.5*reg68; reg80=0.5*reg67; reg81=reg62-reg60;
    reg82=reg55-reg59; reg84=0.5*reg72; reg85=0.5*reg71; reg86=reg21*var_inter[1]; reg88=reg16*var_inter[0];
    reg89=0.5*reg78; reg90=reg76*reg72; reg91=reg74*reg85; reg92=0.5*reg77; reg93=0.5*reg82;
    reg97=reg74*reg79; reg98=reg74*reg84; reg99=reg75*reg71; reg100=reg75*reg67; reg101=reg76*reg68;
    reg102=reg74*reg80; reg103=0.5*reg81; reg105=reg88*elem.f_vol_e[0]; reg106=reg88*elem.f_vol_e[1]; reg107=reg100-reg97;
    reg109=reg74*reg92; reg112=reg76*reg78; reg113=reg86*elem.f_vol_e[1]; reg114=reg74*reg89; reg115=reg75*reg77;
    reg116=reg102-reg101; reg117=reg86*elem.f_vol_e[0]; reg118=reg76*reg82; reg119=reg75*reg81; reg120=reg90-reg91;
    reg121=var_inter[0]*var_inter[1]; reg122=reg98-reg99; reg123=reg74*reg93; reg124=reg16*reg21; reg125=reg74*reg103;
    reg107=reg107-reg117; reg122=reg122-reg105; reg126=reg124*elem.f_vol_e[1]; reg127=reg112+reg109; reg128=reg119+reg123;
    reg129=reg121*elem.f_vol_e[0]; T reg130=reg121*elem.f_vol_e[1]; T reg131=reg118+reg125; T reg132=reg115+reg114; T reg133=reg124*elem.f_vol_e[0];
    reg116=reg116-reg113; reg120=reg120-reg106; T reg134=reg131+reg126; reg116=reg11*reg116; reg107=reg11*reg107;
    T reg135=reg127+reg130; T reg136=reg133+reg128; T reg137=reg132+reg129; reg122=reg11*reg122; reg120=reg11*reg120;
    reg116=ponderation*reg116; reg107=ponderation*reg107; T reg138=reg11*reg135; T reg139=reg11*reg137; reg120=ponderation*reg120;
    reg122=ponderation*reg122; T reg140=reg11*reg136; T reg141=reg11*reg134; sollicitation[indices[3]+0]+=-reg107; sollicitation[indices[3]+1]+=-reg116;
    reg107=ponderation*reg138; sollicitation[indices[2]+1]+=reg107; reg116=ponderation*reg139; sollicitation[indices[2]+0]+=reg116; sollicitation[indices[1]+1]+=-reg120;
    sollicitation[indices[1]+0]+=-reg122; reg120=ponderation*reg140; sollicitation[indices[0]+0]+=reg120; reg122=ponderation*reg141; sollicitation[indices[0]+1]+=reg122;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; T reg2=2*(*f.m).shear_modulus_12; reg0=1.0/reg0; reg1=1.0/reg1;
    T reg3=reg1*reg0; reg2=1.0/reg2; T reg4=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg5=pow((*f.m).v1[0],2); T reg6=1.0/(*f.m).elastic_modulus_3;
    T reg7=pow((*f.m).v1[1],2); T reg8=pow((*f.m).v2[1],2); T reg9=pow((*f.m).v2[0],2); T reg10=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg11=reg2*reg3;
    T reg12=reg10*reg11; T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg14=1.0/(*f.m).elastic_modulus_2; T reg15=reg4*reg11; T reg16=reg6*reg11;
    reg7=reg5+reg7; reg5=pow((*f.m).v1[2],2); reg8=reg9+reg8; reg9=pow((*f.m).v2[2],2); reg5=reg7+reg5;
    reg9=reg8+reg9; reg7=reg4*reg15; reg8=reg13*reg16; T reg17=reg4*reg12; T reg18=reg14*reg16;
    reg9=pow(reg9,0.5); reg5=pow(reg5,0.5); T reg19=reg14*reg12; T reg20=reg13*reg15; reg17=reg8+reg17;
    reg7=reg18-reg7; reg18=1.0/(*f.m).elastic_modulus_1; T reg21=(*f.m).v1[2]/reg5; T reg22=reg13*reg17; T reg23=(*f.m).v1[1]/reg5;
    T reg24=(*f.m).v2[1]/reg9; T reg25=reg20+reg19; T reg26=(*f.m).v2[2]/reg9; T reg27=reg18*reg7; reg5=(*f.m).v1[0]/reg5;
    T reg28=reg10*reg3; T reg29=reg10*reg15; T reg30=reg14*reg11; T reg31=reg4*reg3; T reg32=reg10*reg25;
    reg3=reg6*reg3; T reg33=reg2*reg0; reg22=reg27-reg22; reg27=reg10*reg12; reg11=reg13*reg11;
    reg9=(*f.m).v2[0]/reg9; reg16=reg18*reg16; T reg34=reg23*reg26; T reg35=reg21*reg24; T reg36=reg4*reg33;
    T reg37=2*reg5; T reg38=reg14*reg3; T reg39=2*reg9; reg27=reg16-reg27; reg15=reg18*reg15;
    reg16=reg34-reg35; T reg40=reg10*reg11; reg3=reg13*reg3; T reg41=reg4*reg31; reg32=reg22-reg32;
    reg22=reg21*reg9; T reg42=reg10*reg33; T reg43=reg5*reg26; T reg44=reg10*reg30; reg12=reg13*reg12;
    T reg45=reg2*reg1; reg33=reg6*reg33; T reg46=reg4*reg28; reg29=reg8+reg29; reg8=reg10*reg45;
    reg30=reg18*reg30; reg28=reg14*reg28; T reg47=reg39*reg24; T reg48=2*reg16; T reg49=reg6*reg45;
    T reg50=reg4*reg42; reg41=reg38-reg41; reg38=reg4*reg36; T reg51=reg13*reg33; reg33=reg14*reg33;
    reg46=reg3+reg46; reg3=pow(reg23,2); T reg52=pow(reg5,2); T reg53=reg22-reg43; T reg54=reg5*reg24;
    T reg55=reg23*reg9; reg17=reg17/reg32; reg44=reg20+reg44; reg29=reg29/reg32; reg27=reg27/reg32;
    T reg56=reg23*reg37; T reg57=pow(reg9,2); T reg58=pow(reg24,2); reg40=reg15+reg40; reg45=reg4*reg45;
    reg12=reg15+reg12; reg31=reg13*reg31; reg7=reg7/reg32; reg11=reg13*reg11; reg15=reg4*reg8;
    T reg59=reg3*reg7; T reg60=reg57*reg27; T reg61=reg48*reg53; T reg62=reg58*reg17; T reg63=reg52*reg29;
    T reg64=pow(reg21,2); T reg65=pow(reg53,2); T reg66=reg47*reg17; T reg67=pow(reg16,2); T reg68=reg56*reg7;
    T reg69=pow(reg26,2); reg11=reg30-reg11; reg36=reg13*reg36; reg42=reg14*reg42; reg12=reg12/reg32;
    reg40=reg40/reg32; reg30=reg57*reg17; T reg70=reg52*reg7; reg41=reg18*reg41; reg44=reg44/reg32;
    reg46=reg13*reg46; T reg71=reg31+reg28; reg25=reg25/reg32; T reg72=reg54-reg55; reg38=reg33-reg38;
    reg50=reg51+reg50; reg33=reg14*reg49; reg49=reg13*reg49; reg51=reg4*reg45; T reg73=reg3*reg29;
    T reg74=reg58*reg27; T reg75=reg56*reg29; T reg76=reg47*reg27; T reg77=reg65*reg25; reg62=reg59+reg62;
    reg59=reg65*reg12; T reg78=reg64*reg29; T reg79=reg69*reg27; T reg80=pow(reg72,2); reg11=reg11/reg32;
    T reg81=reg47*reg40; T reg82=reg56*reg44; reg45=reg13*reg45; reg8=reg14*reg8; T reg83=reg58*reg40;
    T reg84=reg3*reg44; reg46=reg41-reg46; reg41=reg57*reg40; T reg85=reg52*reg44; reg71=reg10*reg71;
    reg38=reg18*reg38; reg50=reg13*reg50; T reg86=reg36+reg42; reg15=reg49+reg15; reg76=reg75+reg76;
    reg49=reg61*reg12; reg51=reg33-reg51; reg33=reg64*reg7; reg75=reg69*reg17; T reg87=reg67*reg12;
    reg60=reg63+reg60; reg63=reg61*reg25; reg30=reg70+reg30; reg70=reg67*reg25; reg66=reg68+reg66;
    reg74=reg73+reg74; reg68=reg69*reg40; reg73=reg64*reg44; T reg88=2*reg23; T reg89=reg65*reg11;
    T reg90=reg45+reg8; reg83=reg84+reg83; reg84=reg67*reg11; reg41=reg85+reg41; reg71=reg46-reg71;
    reg70=reg30+reg70; reg50=reg38-reg50; reg86=reg10*reg86; reg87=reg60+reg87; reg15=reg13*reg15;
    reg51=reg18*reg51; reg49=reg76+reg49; reg30=reg80*reg12; reg38=2*reg24; reg46=reg39*reg26;
    reg79=reg78+reg79; reg60=reg5*reg23; reg76=reg80*reg25; reg78=reg21*reg37; reg59=reg74+reg59;
    reg61=reg61*reg11; reg74=reg9*reg24; reg81=reg82+reg81; reg63=reg66+reg63; reg75=reg33+reg75;
    reg77=reg62+reg77; reg33=reg46*reg27; reg62=reg78*reg29; reg66=reg52*reg77; reg86=reg50-reg86;
    reg50=reg60*reg77; reg82=reg74*reg59; reg77=reg3*reg77; reg85=reg53*reg16; reg68=reg73+reg68;
    reg73=reg80*reg11; reg89=reg83+reg89; reg83=reg57*reg59; T reg91=reg60*reg63; T reg92=reg74*reg49;
    reg61=reg81+reg61; reg81=reg21*reg88; T reg93=reg58*reg49; T reg94=reg3*reg63; reg84=reg41+reg84;
    reg71=reg71/reg32; reg41=reg52*reg70; reg59=reg58*reg59; T reg95=reg57*reg87; T reg96=reg58*reg87;
    T reg97=reg3*reg70; reg90=reg10*reg90; T reg98=reg23*reg24; reg63=reg52*reg63; reg76=reg75+reg76;
    reg54=reg55+reg54; reg55=reg5*reg9; reg30=reg79+reg30; reg75=reg23*reg16; reg79=reg38*reg26;
    T reg99=2*reg53; reg49=reg57*reg49; reg48=reg48*reg72; T reg100=reg78*reg7; T reg101=reg5*reg53;
    reg15=reg51-reg15; reg51=reg46*reg17; reg99=reg99*reg72; T reg102=reg52*reg76; reg70=reg60*reg70;
    reg87=reg74*reg87; reg73=reg68+reg73; reg68=reg58*reg30; T reg103=reg3*reg76; T reg104=reg57*reg30;
    T reg105=reg54*reg71; reg90=reg15-reg90; reg82=reg50+reg82; reg15=reg85*reg89; reg50=reg65*reg89;
    T reg106=reg21*reg26; T reg107=reg85*reg61; reg92=reg91+reg92; reg91=reg23*reg53; T reg108=reg5*reg16;
    reg83=reg66+reg83; reg89=reg67*reg89; reg66=reg98*reg71; T reg109=reg55*reg71; reg33=reg62+reg33;
    reg62=reg67*reg61; reg49=reg63+reg49; reg63=reg48*reg12; T reg110=reg9*reg37; reg86=reg86/reg32;
    reg29=reg81*reg29; reg17=reg79*reg17; T reg111=reg9*reg53; reg59=reg77+reg59; reg27=reg79*reg27;
    reg96=reg97+reg96; reg77=reg24*reg16; reg97=reg65*reg84; T reg112=reg67*reg84; T reg113=2*reg21;
    T reg114=reg78*reg44; T reg115=reg46*reg40; reg61=reg65*reg61; reg93=reg94+reg93; reg101=reg75+reg101;
    reg51=reg100+reg51; reg75=reg48*reg25; reg94=reg24*reg88; reg7=reg81*reg7; reg95=reg41+reg95;
    reg41=reg57*reg18; reg100=reg94*reg13; T reg116=reg57*reg13; T reg117=reg110*reg18; reg44=reg81*reg44;
    T reg118=reg58*reg14; reg40=reg79*reg40; T reg119=reg58*reg13; T reg120=reg110*reg13; T reg121=reg94*reg14;
    T reg122=reg94*reg109; reg25=reg99*reg25; reg17=reg7+reg17; reg37=reg16*reg37; reg7=reg110*reg109;
    reg112=reg95+reg112; reg75=reg51+reg75; reg97=reg96+reg97; reg61=reg93+reg61; reg51=reg94*reg105;
    reg88=reg53*reg88; reg113=reg26*reg113; reg93=reg101*reg86; reg91=reg91*reg86; reg108=reg108*reg86;
    reg87=reg70+reg87; reg84=reg85*reg84; reg104=reg102+reg104; reg15=reg82+reg15; reg70=reg54*reg66;
    reg82=reg67*reg73; reg76=reg60*reg76; reg30=reg74*reg30; reg95=reg106*reg71; reg96=reg65*reg73;
    reg68=reg103+reg68; reg102=reg110*reg105; reg62=reg49+reg62; reg49=reg110*reg66; reg103=reg9*reg16;
    T reg123=reg24*reg53; reg66=reg94*reg66; reg50=reg59+reg50; reg12=reg99*reg12; reg27=reg29+reg27;
    reg89=reg83+reg89; reg29=reg21*reg72; reg32=reg90/reg32; reg63=reg33+reg63; reg107=reg92+reg107;
    reg105=reg54*reg105; reg48=reg48*reg11; reg115=reg114+reg115; reg111=reg77+reg111; reg33=reg101*reg91;
    reg66=reg50+reg66; reg50=reg88*reg91; reg70=reg15+reg70; reg105=reg107+reg105; reg15=reg101*reg93;
    reg109=reg54*reg109; reg84=reg87+reg84; reg59=reg57*reg63; reg96=reg68+reg96; reg29=reg29*reg86;
    reg68=reg52*reg75; reg103=reg103*reg32; reg25=reg17+reg25; reg39=reg39*reg16; reg17=reg37*reg108;
    reg77=(*f.m).alpha_2*reg58; reg83=(*f.m).alpha_1*reg3; reg7=reg112+reg7; reg87=(*f.m).alpha_2*reg57; reg90=reg110*reg95;
    reg12=reg27+reg12; reg27=reg26*reg72; reg92=reg5*reg72; reg107=reg21*reg16; reg112=(*f.m).alpha_1*reg52;
    reg102=reg62+reg102; reg51=reg61+reg51; reg61=reg88*reg93; reg93=reg37*reg93; reg62=reg3*reg75;
    reg114=reg58*reg63; reg38=reg38*reg53; T reg124=reg111*reg32; reg123=reg123*reg32; T reg125=reg57*reg10;
    reg49=reg89+reg49; reg11=reg99*reg11; reg120=reg121-reg120; reg89=reg113*reg4; reg91=reg37*reg91;
    reg116=reg118-reg116; reg99=reg69*reg4; reg118=reg94*reg95; reg100=reg117-reg100; reg117=reg113*reg10;
    reg22=reg43+reg22; reg119=reg41-reg119; reg41=reg69*reg10; reg40=reg44+reg40; reg48=reg115+reg48;
    reg122=reg97+reg122; reg43=reg88*reg108; reg44=reg58*reg4; reg97=reg94*reg4; reg30=reg76+reg30;
    reg73=reg85*reg73; reg76=reg110*reg10; reg82=reg104+reg82; reg104=reg57*reg12; reg95=reg54*reg95;
    reg115=reg111*reg123; reg33=reg70+reg33; reg41=reg119-reg41; reg118=reg96+reg118; reg97=reg76+reg97;
    reg70=reg23*reg72; reg18=reg52*reg18; reg76=reg21*reg53; reg92=reg107+reg92; reg113=reg113*reg6;
    reg117=reg100-reg117; reg96=reg52*reg13; reg100=reg38*reg103; reg43=reg122+reg43; reg107=reg9*reg72;
    reg119=reg37*reg29; reg121=reg26*reg16; reg122=reg39*reg103; T reg126=reg88*reg29; reg17=reg7+reg17;
    reg7=(*f.m).alpha_2*reg69; T reg127=(*f.m).alpha_1*reg64; reg73=reg30+reg73; reg30=(*f.m).alpha_3*reg65; reg77=reg83+reg77;
    reg90=reg82+reg90; reg13=reg3*reg13; reg82=reg52*reg25; reg83=(*f.m).alpha_3*reg67; reg87=reg112+reg87;
    reg11=reg40+reg11; reg40=reg38*reg123; reg89=reg120-reg89; reg63=reg74*reg63; reg75=reg60*reg75;
    reg123=reg39*reg123; reg112=reg58*reg12; reg120=reg3*reg25; reg50=reg66+reg50; reg66=reg22*reg71;
    reg59=reg68+reg59; reg68=reg39*reg124; T reg128=reg111*reg124; reg15=reg105+reg15; reg27=reg27*reg32;
    reg105=reg69*reg6; reg93=reg102+reg93; reg102=reg67*reg48; reg99=reg116-reg99; reg116=reg65*reg48;
    reg114=reg62+reg114; reg44=reg125+reg44; reg91=reg49+reg91; reg124=reg38*reg124; reg61=reg51+reg61;
    reg34=reg35+reg34; reg108=reg101*reg108; reg14=reg3*reg14; reg109=reg84+reg109; reg119=reg90+reg119;
    reg123=reg91+reg123; reg9=reg9*reg26; reg35=reg39*reg27; reg68=reg93+reg68; reg102=reg59+reg102;
    reg5=reg5*reg21; reg49=reg110*reg66; reg67=reg67*reg11; reg104=reg82+reg104; reg12=reg74*reg12;
    reg51=reg24*reg72; reg59=reg26*reg53; reg107=reg121+reg107; reg70=reg76+reg70; reg62=reg89*reg58;
    reg76=reg117*reg57; reg82=reg99*reg58; reg84=reg41*reg57; reg90=reg89*reg3; reg91=reg117*reg52;
    reg93=reg99*reg3; reg121=reg41*reg52; reg29=reg101*reg29; reg95=reg73+reg95; reg115=reg33+reg115;
    reg103=reg111*reg103; reg108=reg109+reg108; reg71=reg34*reg71; reg92=reg92*reg86; reg65=reg65*reg11;
    reg112=reg120+reg112; reg33=reg94*reg66; reg116=reg114+reg116; reg124=reg61+reg124; reg126=reg118+reg126;
    reg122=reg17+reg122; reg17=reg38*reg27; reg61=(*f.m).alpha_2*reg74; reg73=(*f.m).alpha_1*reg60; reg80=(*f.m).alpha_3*reg80;
    reg7=reg127+reg7; reg30=reg77+reg30; reg83=reg87+reg83; reg97=reg113-reg97; reg44=reg105-reg44;
    reg77=reg3*reg4; reg87=reg52*reg10; reg96=reg14-reg96; reg4=reg64*reg4; reg13=reg18-reg13;
    reg10=reg64*reg10; reg100=reg43+reg100; reg25=reg60*reg25; reg48=reg85*reg48; reg63=reg75+reg63;
    reg128=reg15+reg128; reg40=reg50+reg40; reg41=reg55*reg41; reg99=reg98*reg99; reg117=reg55*reg117;
    reg89=reg98*reg89; reg14=reg124*reg115; reg15=reg128*reg123; reg18=reg115*reg68; reg43=reg128*reg40;
    reg50=reg83*reg122; reg75=reg30*reg123; reg105=reg83*reg100; reg109=reg30*reg40; reg66=reg54*reg66;
    reg12=reg25+reg12; reg107=reg107*reg32; reg11=reg85*reg11; reg48=reg63+reg48; reg86=reg70*reg86;
    reg16=reg72*reg16; reg21=reg23*reg21; reg26=reg24*reg26; reg29=reg95+reg29; reg27=reg111*reg27;
    reg85=(*f.m).alpha_3*reg85; reg61=reg73+reg61; reg80=reg7+reg80; reg7=reg54*reg2; reg74=reg74*reg2;
    reg77=reg87+reg77; reg6=reg64*reg6; reg4=reg96-reg4; reg10=reg13-reg10; reg13=reg106*reg97;
    reg89=reg117+reg89; reg93=reg121+reg93; reg23=reg44*reg64; reg24=reg106*reg44; reg99=reg41+reg99;
    reg90=reg91+reg90; reg25=reg97*reg64; reg51=reg59+reg51; reg82=reg84+reg82; reg44=reg69*reg44;
    reg97=reg69*reg97; reg62=reg76+reg62; reg17=reg126+reg17; reg35=reg119+reg35; reg49=reg102+reg49;
    reg41=reg37*reg92; reg67=reg104+reg67; reg33=reg116+reg33; reg59=reg88*reg92; reg103=reg108+reg103;
    reg63=(*f.m).alpha_1*reg5; reg70=(*f.m).alpha_2*reg9; reg94=reg94*reg71; reg65=reg112+reg65; reg110=reg110*reg71;
    reg77=reg6-reg77; reg97=reg62+reg97; reg6=reg47*reg7; reg75=reg50+reg75; reg50=(*f.m).alpha_2*reg26;
    reg62=(*f.m).alpha_1*reg21; reg35=reg80*reg35; reg110=reg67+reg110; reg109=reg105+reg109; reg17=reg80*reg17;
    reg71=reg54*reg71; reg11=reg12+reg11; reg92=reg101*reg92; reg66=reg48+reg66; reg85=reg61+reg85;
    reg70=reg63+reg70; reg41=reg49+reg41; reg53=reg72*reg53; reg24=reg99+reg24; reg12=reg54*reg74;
    reg48=reg39*reg107; reg16=(*f.m).alpha_3*reg16; reg37=reg37*reg86; reg49=reg54*reg7; reg13=reg89+reg13;
    reg61=reg103*reg83; reg63=reg115*reg30; reg67=reg68*reg40; reg9=reg9*reg1; reg72=reg22*reg1;
    reg73=reg124*reg123; reg18=reg15-reg18; reg15=reg10*reg52; reg76=reg4*reg3; reg27=reg29+reg27;
    reg14=reg43-reg14; reg23=reg93+reg23; reg29=reg74*reg56; reg88=reg88*reg86; reg94=reg65+reg94;
    reg25=reg90+reg25; reg7=reg7*reg56; reg43=reg38*reg107; reg65=1-var_inter[1]; reg74=reg47*reg74;
    reg44=reg82+reg44; reg82=1-var_inter[0]; reg84=reg10*reg57; reg87=reg4*reg58; reg59=reg33+reg59;
    reg32=reg51*reg32; reg48=reg41+reg48; reg33=reg65*elem.pos(0)[0]; reg41=reg65*elem.pos(1)[0]; reg38=reg38*reg32;
    reg39=reg39*reg32; reg16=reg70+reg16; reg37=reg110+reg37; reg2=reg60*reg2; reg51=reg65*elem.pos(1)[1];
    reg43=reg59+reg43; reg88=reg94+reg88; reg59=reg9*reg46; reg74=reg44+reg74; reg69=reg69*reg77;
    reg6=reg97+reg6; reg44=reg72*reg46; reg87=reg84+reg87; reg60=reg82*elem.pos(0)[0]; reg70=elem.pos(1)[0]*var_inter[0];
    reg84=elem.pos(1)[1]*var_inter[0]; reg89=reg82*elem.pos(0)[1]; reg86=reg101*reg86; reg71=reg11+reg71; reg11=reg72*reg78;
    reg7=reg25+reg7; reg107=reg111*reg107; reg92=reg66+reg92; reg10=reg55*reg10; reg4=reg98*reg4;
    reg25=reg9*reg78; reg29=reg23+reg29; reg12=reg24+reg12; reg9=reg9*reg22; reg64=reg77*reg64;
    reg76=reg15+reg76; reg15=reg34*reg0; reg26=reg26*reg0; reg49=reg13+reg49; reg72=reg72*reg22;
    reg13=reg65*elem.pos(0)[1]; reg67=reg73-reg67; reg23=reg124*reg85; reg17=reg109+reg17; reg35=reg75+reg35;
    reg24=reg85*reg68; reg63=reg61+reg63; reg80=reg27*reg80; reg27=reg18*reg100; reg61=reg14*reg122;
    reg53=(*f.m).alpha_3*reg53; reg50=reg62+reg50; reg23=reg17+reg23; reg11=reg7+reg11; reg43=reg43*reg16;
    reg7=reg103*reg68; reg17=reg15*reg34; reg62=reg15*reg81; reg66=reg85*reg128; reg107=reg92+reg107;
    reg69=reg87+reg69; reg73=reg26*reg81; reg75=reg47*reg2; reg86=reg71+reg86; reg59=reg74+reg59;
    reg71=reg26*reg79; reg32=reg111*reg32; reg44=reg6+reg44; reg15=reg15*reg79; reg24=reg35+reg24;
    reg48=reg16*reg48; reg6=reg60+reg70; reg35=elem.pos(2)[0]*var_inter[0]; reg74=elem.pos(2)[1]*var_inter[0]; reg87=reg84+reg89;
    reg53=reg50+reg53; reg38=reg88+reg38; reg41=reg41-reg33; reg50=elem.pos(2)[0]*var_inter[1]; reg9=reg12+reg9;
    reg26=reg26*reg34; reg12=elem.pos(2)[1]*var_inter[1]; reg27=reg61-reg27; reg61=reg103*reg67; reg51=reg51-reg13;
    reg77=reg106*reg77; reg25=reg29+reg25; reg72=reg49+reg72; reg29=reg128*reg122; reg80=reg63+reg80;
    reg49=reg2*reg56; reg64=reg76+reg64; reg39=reg37+reg39; reg37=reg124*reg103; reg63=reg128*reg100;
    reg4=reg10+reg4; reg1=reg5*reg1; reg38=reg38*reg53; reg32=reg86+reg32; reg2=reg54*reg2;
    reg43=reg23+reg43; reg5=elem.pos(3)[1]*reg82; reg74=reg74-reg87; reg77=reg4+reg77; reg61=reg27+reg61;
    reg50=reg41+reg50; reg4=elem.pos(3)[0]*var_inter[1]; reg51=reg12+reg51; reg10=elem.pos(3)[1]*var_inter[1]; reg0=reg21*reg0;
    reg12=reg115*reg100; reg37=reg63-reg37; reg49=reg64+reg49; reg78=reg1*reg78; reg21=reg103*reg40;
    reg73=reg25+reg73; reg23=reg115*reg122; reg7=reg29-reg7; reg25=reg103*reg123; reg62=reg11+reg62;
    reg11=reg124*reg122; reg75=reg69+reg75; reg46=reg1*reg46; reg27=reg68*reg100; reg71=reg59+reg71;
    reg15=reg44+reg15; reg48=reg24+reg48; reg35=reg35-reg6; reg24=elem.pos(3)[0]*reg82; reg39=reg53*reg39;
    reg26=reg9+reg26; reg107=reg16*reg107; reg66=reg80+reg66; reg17=reg72+reg17; reg18=reg18/reg61;
    reg21=reg12-reg21; reg107=reg66+reg107; reg7=reg7/reg61; reg32=reg53*reg32; reg9=var_inter[0]*vectors[0][indices[1]+1];
    reg12=reg82*vectors[0][indices[0]+0]; reg16=reg62*reg26; reg5=reg74+reg5; reg29=var_inter[0]*vectors[0][indices[1]+0]; reg41=reg65*vectors[0][indices[0]+0];
    reg39=reg48+reg39; reg24=reg35+reg24; reg25=reg23-reg25; reg46=reg75+reg46; reg23=reg82*vectors[0][indices[0]+1];
    reg79=reg0*reg79; reg51=reg51-reg10; reg35=reg122*reg40; reg27=reg11-reg27; reg11=reg123*reg100;
    reg44=reg65*vectors[0][indices[0]+1]; reg48=reg65*vectors[0][indices[1]+1]; reg22=reg1*reg22; reg2=reg77+reg2; reg1=reg73*reg17;
    reg53=reg71*reg17; reg50=reg50-reg4; reg14=reg14/reg61; reg38=reg43+reg38; reg43=reg65*vectors[0][indices[1]+0];
    reg59=reg15*reg26; reg37=reg37/reg61; reg81=reg0*reg81; reg78=reg49+reg78; reg41=reg43-reg41;
    reg44=reg48-reg44; reg43=var_inter[1]*vectors[0][indices[2]+0]; reg22=reg2+reg22; reg2=var_inter[1]*vectors[0][indices[2]+1]; reg34=reg0*reg34;
    reg0=var_inter[0]*vectors[0][indices[2]+1]; reg48=reg73*reg15; reg14=reg14*reg39; reg79=reg46+reg79; reg27=reg27/reg61;
    reg18=reg18*reg38; reg16=reg1-reg16; reg11=reg35-reg11; reg37=reg37*reg39; reg7=reg7*reg38;
    reg1=reg51*reg24; reg35=var_inter[0]*vectors[0][indices[2]+0]; reg59=reg53-reg59; reg46=reg62*reg71; reg32=reg107+reg32;
    reg25=reg25/reg61; reg67=reg67/reg61; reg9=reg23+reg9; reg12=reg29+reg12; reg23=reg50*reg5;
    reg21=reg21/reg61; reg81=reg78+reg81; reg43=reg41+reg43; reg12=reg35-reg12; reg29=var_inter[1]*vectors[0][indices[3]+0];
    reg1=reg23-reg1; reg67=reg67*reg32; reg18=reg14-reg18; reg2=reg44+reg2; reg27=reg27*reg32;
    reg37=reg7-reg37; reg39=reg21*reg39; reg38=reg25*reg38; reg61=reg11/reg61; reg7=var_inter[1]*vectors[0][indices[3]+1];
    reg46=reg48-reg46; reg9=reg0-reg9; reg0=reg82*vectors[0][indices[3]+0]; reg11=reg82*vectors[0][indices[3]+1]; reg34=reg22+reg34;
    reg14=reg79*reg16; reg21=reg81*reg59; reg22=reg79*reg26; reg23=reg71*reg34; reg25=reg73*reg34;
    reg26=reg81*reg26; reg35=reg34*reg46; reg11=reg9+reg11; reg12=reg0+reg12; reg51=reg51/reg1;
    reg0=1-(*f.m).resolution; reg24=reg24/reg1; reg7=reg2-reg7; reg14=reg21-reg14; reg38=reg39-reg38;
    reg32=reg61*reg32; reg27=reg37-reg27; reg18=reg67+reg18; reg29=reg43-reg29; reg50=reg50/reg1;
    reg5=reg5/reg1; reg2=reg62*reg34; reg25=reg26-reg25; reg23=reg22-reg23; reg9=reg51*reg11;
    reg21=reg30*(*f.m).resolution; reg18=reg18*reg0; reg35=reg14+reg35; reg14=reg5*reg7; reg22=reg29*reg24;
    reg26=reg12*reg50; reg71=reg81*reg71; reg37=reg81*reg17; reg34=reg15*reg34; reg73=reg73*reg79;
    reg27=reg27*reg0; reg39=reg83*(*f.m).resolution; reg17=reg79*reg17; reg38=reg32+reg38; reg32=reg85*(*f.m).resolution;
    reg38=reg38*reg0; reg18=reg39+reg18; reg29=reg29*reg5; reg15=reg81*reg15; reg73=reg71-reg73;
    reg79=reg62*reg79; reg27=reg21+reg27; reg7=reg24*reg7; reg34=reg17-reg34; reg22=reg26-reg22;
    reg11=reg50*reg11; reg23=reg23/reg35; reg9=reg14-reg9; reg12=reg12*reg51; reg2=reg37-reg2;
    reg25=reg25/reg35; reg12=reg29-reg12; reg103=reg103*reg0; reg115=reg115*reg0; reg7=reg11-reg7;
    reg38=reg32+reg38; reg34=reg34/reg35; reg18=(*f.m).deltaT*reg18; reg27=reg27*(*f.m).deltaT; reg73=reg73/reg35;
    reg59=reg59/reg35; reg79=reg15-reg79; reg2=reg2/reg35; reg11=(*f.m).resolution*reg23; reg14=(*f.m).resolution*reg25;
    reg9=reg22+reg9; reg16=reg16/reg35; reg122=reg0*reg122; reg40=reg0*reg40; reg100=reg0*reg100;
    reg79=reg79/reg35; reg35=reg46/reg35; reg123=reg0*reg123; reg7=reg7-reg27; reg38=reg38*(*f.m).deltaT;
    reg12=reg12-reg18; reg9=0.5*reg9; reg15=(*f.m).resolution*reg34; reg17=(*f.m).resolution*reg2; reg21=(*f.m).resolution*reg59;
    reg128=reg0*reg128; reg22=(*f.m).resolution*reg73; reg26=(*f.m).resolution*reg16; reg11=reg103+reg11; reg14=reg115-reg14;
    reg9=reg9-reg38; reg68=reg0*reg68; reg29=(*f.m).resolution*reg79; reg22=reg128+reg22; reg32=reg12*reg11;
    reg37=reg7*reg14; reg0=reg124*reg0; reg122=reg21+reg122; reg21=(*f.m).resolution*reg35; reg26=reg123-reg26;
    reg15=reg100-reg15; reg17=reg40+reg17; reg39=var_inter[1]*reg24; reg40=reg7*reg26; reg41=reg12*reg122;
    reg43=var_inter[0]*reg50; reg44=var_inter[0]*reg51; reg12=reg12*reg15; reg7=reg7*reg17; reg46=reg65*reg5;
    reg48=reg65*reg24; reg49=reg82*reg51; reg53=reg82*reg50; reg37=reg32+reg37; reg68=reg21+reg68;
    reg21=reg9*reg22; reg32=var_inter[1]*reg5; reg29=reg0-reg29; reg40=reg41+reg40; reg0=reg9*reg68;
    reg41=reg48+reg43; reg61=reg46+reg44; reg7=reg12+reg7; reg12=reg48-reg53; reg62=reg49-reg46;
    reg21=reg37+reg21; reg37=reg39+reg53; reg63=reg49+reg32; reg9=reg9*reg29; reg64=reg43-reg39;
    reg66=reg32-reg44; reg0=reg40+reg0; reg21=2*reg21; reg7=reg9+reg7; reg9=0.5*reg64;
    reg40=0.5*reg66; reg67=0.5*reg37; reg69=0.5*reg63; reg71=0.5*reg12; reg72=0.5*reg62;
    reg74=0.5*reg41; reg75=0.5*reg61; reg76=reg21*reg75; reg77=reg7*reg41; reg78=reg21*reg74;
    reg80=reg0*reg61; reg81=reg21*reg69; reg86=reg7*reg37; reg88=reg21*reg72; reg90=reg7*reg12;
    reg91=reg21*reg71; reg92=reg0*reg62; reg93=reg82*var_inter[1]; reg94=reg65*reg82; reg95=var_inter[0]*var_inter[1];
    reg96=reg65*var_inter[0]; reg97=reg21*reg67; reg99=reg0*reg63; reg100=reg21*reg40; reg101=reg7*reg64;
    reg102=reg21*reg9; reg103=reg0*reg66; reg104=reg93*elem.f_vol_e[0]; reg97=reg97-reg99; reg105=reg95*elem.f_vol_e[1];
    reg100=reg101+reg100; reg91=reg92+reg91; reg92=reg94*elem.f_vol_e[0]; reg86=reg86-reg81; reg101=reg93*elem.f_vol_e[1];
    reg106=reg95*elem.f_vol_e[0]; reg102=reg103+reg102; reg88=reg90+reg88; reg90=reg94*elem.f_vol_e[1]; reg76=reg76-reg77;
    reg103=reg96*elem.f_vol_e[1]; reg80=reg80-reg78; reg107=reg96*elem.f_vol_e[0]; reg91=reg91-reg92; reg76=reg76-reg103;
    reg80=reg80-reg107; reg100=reg100-reg105; reg102=reg102-reg106; reg88=reg88-reg90; reg86=reg86-reg101;
    reg97=reg97-reg104; reg100=reg1*reg100; reg80=reg1*reg80; reg76=reg1*reg76; reg91=reg1*reg91;
    reg97=reg1*reg97; reg102=reg1*reg102; reg88=reg1*reg88; reg86=reg1*reg86; sollicitation[indices[0]+0]+=ponderation*reg91;
    sollicitation[indices[0]+1]+=ponderation*reg88; sollicitation[indices[2]+0]+=ponderation*reg102; sollicitation[indices[3]+0]+=ponderation*reg97; sollicitation[indices[2]+1]+=ponderation*reg100; sollicitation[indices[1]+0]+=ponderation*reg80;
    sollicitation[indices[3]+1]+=ponderation*reg86; sollicitation[indices[1]+1]+=ponderation*reg76;
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

