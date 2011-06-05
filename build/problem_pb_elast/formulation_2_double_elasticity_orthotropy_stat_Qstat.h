
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
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=pow((*f.m).v1[1],2); T reg5=pow((*f.m).v1[0],2); T reg6=pow((*f.m).v2[0],2);
    T reg7=pow((*f.m).v2[1],2); T reg8=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=reg2*reg3; reg4=reg5+reg4;
    reg5=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg11=pow((*f.m).v1[2],2); T reg12=reg10*reg8; T reg13=reg5*reg10; T reg14=reg10*reg9;
    reg11=reg4+reg11; reg4=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=pow((*f.m).v2[2],2); reg7=reg6+reg7;
    reg6=reg5*reg12; reg16=reg7+reg16; reg7=reg4*reg14; T reg17=reg5*reg13; reg11=pow(reg11,0.5);
    T reg18=reg15*reg14; T reg19=(*f.m).v1[0]/reg11; T reg20=1.0/(*f.m).elastic_modulus_1; T reg21=reg15*reg12; T reg22=reg4*reg13;
    reg6=reg7+reg6; T reg23=(*f.m).v1[1]/reg11; reg17=reg18-reg17; reg16=pow(reg16,0.5); reg18=(*f.m).v2[2]/reg16;
    reg11=(*f.m).v1[2]/reg11; T reg24=(*f.m).v2[1]/reg16; reg16=(*f.m).v2[0]/reg16; T reg25=reg22+reg21; T reg26=2*reg19;
    T reg27=2*reg23; T reg28=reg20*reg17; T reg29=reg4*reg6; T reg30=reg27*reg24; T reg31=reg26*reg16;
    T reg32=reg5*reg3; T reg33=reg3*reg8; T reg34=2*reg11; T reg35=reg23*reg18; T reg36=reg11*reg24;
    T reg37=reg2*reg0; T reg38=reg4*reg10; T reg39=reg12*reg8; reg14=reg20*reg14; reg3=reg3*reg9;
    T reg40=reg13*reg8; reg10=reg15*reg10; T reg41=reg25*reg8; reg29=reg28-reg29; reg28=pow(reg24,2);
    T reg42=pow(reg16,2); T reg43=reg4*reg3; T reg44=reg2*reg1; T reg45=reg5*reg32; T reg46=pow(reg18,2);
    T reg47=reg38*reg8; T reg48=reg37*reg9; T reg49=2*reg16; T reg50=reg31*reg20; T reg51=reg37*reg8;
    reg3=reg15*reg3; T reg52=reg19*reg18; T reg53=reg11*reg16; reg13=reg20*reg13; T reg54=reg30*reg4;
    T reg55=reg28*reg4; T reg56=reg28*reg15; reg39=reg14-reg39; reg14=reg31*reg4; reg41=reg29-reg41;
    reg29=reg30*reg15; T reg57=reg42*reg4; reg37=reg5*reg37; reg12=reg4*reg12; T reg58=reg5*reg33;
    reg40=reg7+reg40; reg7=reg42*reg20; T reg59=reg34*reg18; T reg60=reg35-reg36; T reg61=reg10*reg8;
    T reg62=reg49*reg24; reg33=reg15*reg33; reg47=reg13+reg47; T reg63=reg5*reg51; reg17=reg17/reg41;
    T reg64=reg44*reg9; T reg65=reg23*reg26; reg55=reg7-reg55; reg7=reg46*reg8; reg40=reg40/reg41;
    reg61=reg22+reg61; reg6=reg6/reg41; T reg66=reg30*reg5; reg39=reg39/reg41; T reg67=reg46*reg5;
    reg57=reg56-reg57; reg56=reg59*reg5; reg14=reg29-reg14; reg45=reg3-reg45; reg32=reg4*reg32;
    reg58=reg43+reg58; reg54=reg50-reg54; reg3=reg59*reg8; reg29=reg53-reg52; reg43=reg19*reg24;
    reg50=reg23*reg16; T reg68=reg42*reg8; T reg69=2*reg60; reg10=reg20*reg10; T reg70=pow(reg19,2);
    T reg71=pow(reg23,2); T reg72=reg44*reg8; reg44=reg5*reg44; T reg73=reg31*reg8; reg12=reg13+reg12;
    reg13=reg28*reg5; T reg74=reg5*reg37; T reg75=reg4*reg48; reg48=reg15*reg48; reg38=reg4*reg38;
    T reg76=reg17*reg71; reg37=reg4*reg37; T reg77=reg65*reg17; T reg78=reg28*reg6; T reg79=reg4*reg71;
    reg58=reg4*reg58; reg45=reg20*reg45; T reg80=reg32+reg33; T reg81=pow(reg11,2); reg74=reg48-reg74;
    reg13=reg68+reg13; reg51=reg15*reg51; reg63=reg75+reg63; reg48=reg42*reg39; reg68=reg40*reg70;
    reg66=reg73+reg66; reg73=reg15*reg64; reg64=reg4*reg64; reg75=reg5*reg44; T reg82=reg59*reg9;
    T reg83=reg5*reg72; T reg84=reg62*reg6; reg7=reg55-reg7; reg55=reg62*reg39; T reg85=reg65*reg40;
    T reg86=pow(reg60,2); T reg87=pow(reg29,2); T reg88=reg69*reg29; T reg89=reg19*reg16; T reg90=reg23*reg24;
    T reg91=reg28*reg39; T reg92=reg40*reg71; reg3=reg54-reg3; reg54=reg4*reg70; T reg93=reg15*reg71;
    reg67=reg57-reg67; reg56=reg14-reg56; reg14=reg43-reg50; reg57=reg20*reg70; reg12=reg12/reg41;
    reg47=reg47/reg41; T reg94=reg42*reg6; T reg95=reg17*reg70; reg38=reg10-reg38; reg10=reg46*reg9;
    reg61=reg61/reg41; reg25=reg25/reg41; T reg96=reg88*reg12; T reg97=reg61*reg70; T reg98=reg42*reg47;
    T reg99=reg90*reg56; T reg100=reg71*reg56; T reg101=reg88*reg25; reg84=reg77+reg84; reg77=reg61*reg71;
    T reg102=reg28*reg47; T reg103=reg70*reg3; T reg104=reg46*reg6; T reg105=reg17*reg81; T reg106=reg87*reg25;
    T reg107=reg62*reg47; T reg108=reg65*reg61; reg55=reg85+reg55; reg85=reg71*reg67; T reg109=reg89*reg7;
    T reg110=reg46*reg39; T reg111=reg40*reg81; T reg112=reg87*reg12; reg91=reg92+reg91; reg92=reg28*reg56;
    T reg113=reg86*reg12; reg48=reg68+reg48; reg68=reg70*reg7; T reg114=reg42*reg3; T reg115=reg28*reg67;
    T reg116=reg42*reg7; T reg117=pow(reg14,2); T reg118=reg89*reg3; T reg119=reg37+reg51; reg63=reg4*reg63;
    reg74=reg20*reg74; reg80=reg80*reg8; T reg120=reg11*reg18; reg43=reg50+reg43; reg50=reg5*reg81;
    reg54=reg93-reg54; reg58=reg45-reg58; reg72=reg15*reg72; reg45=reg70*reg8; reg44=reg4*reg44;
    reg38=reg38/reg41; reg93=reg16*reg24; T reg121=reg5*reg71; T reg122=reg86*reg25; reg94=reg95+reg94;
    reg13=reg10-reg13; reg10=reg90*reg67; reg66=reg82-reg66; reg83=reg64+reg83; reg79=reg57-reg79;
    reg57=reg81*reg8; reg75=reg73-reg75; reg78=reg76+reg78; reg64=reg81*reg9; reg121=reg45+reg121;
    reg115=reg116+reg115; reg45=reg16*reg18; reg73=reg46*reg13; reg92=reg114+reg92; reg99=reg118+reg99;
    reg76=reg19*reg23; reg82=reg46*reg66; reg95=reg120*reg66; reg50=reg54-reg50; reg54=reg81*reg66;
    reg100=reg103+reg100; reg113=reg48+reg113; reg102=reg77+reg102; reg57=reg79-reg57; reg48=reg120*reg13;
    reg77=reg86*reg38; reg79=reg87*reg38; reg98=reg97+reg98; reg10=reg109+reg10; reg97=reg26*reg11;
    reg96=reg55+reg96; reg55=reg61*reg81; reg103=reg49*reg18; reg109=2*reg24; reg114=reg46*reg47;
    reg116=reg12*reg117; reg110=reg111+reg110; reg112=reg91+reg112; reg107=reg108+reg107; reg91=reg88*reg38;
    reg53=reg52+reg53; reg106=reg78+reg106; reg52=reg44+reg72; reg104=reg105+reg104; reg78=reg25*reg117;
    reg83=reg4*reg83; reg101=reg84+reg101; reg75=reg20*reg75; reg119=reg119*reg8; reg84=reg93*reg2;
    reg63=reg74-reg63; reg74=reg43*reg2; reg80=reg58-reg80; reg122=reg94+reg122; reg85=reg68+reg85;
    reg58=reg81*reg13; reg68=reg42*reg96; reg52=reg52*reg8; reg91=reg107+reg91; reg94=reg97*reg17;
    reg105=reg103*reg6; reg107=reg43*reg84; reg48=reg10+reg48; reg10=reg24*reg18; reg108=reg38*reg117;
    reg114=reg55+reg114; reg73=reg115+reg73; reg55=reg62*reg84; reg79=reg102+reg79; reg102=reg29*reg60;
    reg119=reg63-reg119; reg77=reg98+reg77; reg63=reg28*reg50; reg98=reg42*reg57; reg111=reg106*reg70;
    reg115=reg42*reg113; reg118=reg122*reg70; T reg123=reg42*reg112; T reg124=reg93*reg96; T reg125=reg65*reg74;
    T reg126=reg76*reg101; reg54=reg100+reg54; reg78=reg104+reg78; reg95=reg99+reg95; reg99=reg43*reg74;
    reg100=reg93*reg112; reg104=reg76*reg106; reg83=reg75-reg83; reg75=reg101*reg70; T reg127=reg101*reg71;
    reg80=reg80/reg41; T reg128=reg109*reg18; reg116=reg110+reg116; reg110=reg45*reg1; reg69=reg69*reg14;
    T reg129=2*reg29; T reg130=reg53*reg1; T reg131=reg122*reg71; T reg132=reg28*reg113; T reg133=reg70*reg57;
    T reg134=reg71*reg50; reg35=reg36+reg35; reg36=reg23*reg60; T reg135=reg19*reg29; T reg136=reg28*reg112;
    T reg137=reg106*reg71; reg121=reg64-reg121; reg82=reg92+reg82; reg64=reg65*reg84; reg92=reg62*reg74;
    reg58=reg85+reg58; reg85=reg27*reg11; T reg138=reg28*reg96; T reg139=reg97*reg40; T reg140=reg103*reg39;
    reg52=reg83-reg52; reg63=reg98+reg63; reg83=reg102*reg91; reg136=reg137+reg136; reg98=reg87*reg79;
    reg137=reg78*reg71; reg123=reg111+reg123; reg111=reg86*reg79; reg124=reg126+reg124; reg126=reg87*reg77;
    reg132=reg131+reg132; reg131=reg86*reg91; T reg141=reg78*reg70; reg68=reg75+reg68; reg75=reg28*reg116;
    T reg142=reg103*reg110; T reg143=reg42*reg116; reg55=reg73+reg55; reg73=reg103*reg130; T reg144=reg102*reg79;
    reg100=reg104+reg100; reg92=reg82+reg92; reg82=reg93*reg113; reg104=reg76*reg122; T reg145=reg46*reg121;
    T reg146=reg87*reg91; reg138=reg127+reg138; reg107=reg48+reg107; reg48=reg103*reg47; reg127=reg97*reg61;
    reg108=reg114+reg108; reg105=reg94+reg105; reg94=reg69*reg25; reg17=reg85*reg17; reg6=reg128*reg6;
    reg114=reg76*reg2; T reg147=reg90*reg50; T reg148=reg97*reg110; reg64=reg58+reg64; reg39=reg128*reg39;
    reg40=reg85*reg40; reg58=reg69*reg12; reg140=reg139+reg140; reg119=reg119/reg41; reg129=reg129*reg14;
    reg139=reg10*reg0; T reg149=reg35*reg0; T reg150=reg19*reg60; T reg151=reg23*reg29; reg135=reg36+reg135;
    reg134=reg133+reg134; reg36=reg81*reg121; reg133=reg60*reg24; T reg152=reg29*reg16; T reg153=reg89*reg57;
    reg125=reg54+reg125; reg99=reg95+reg99; reg54=reg43*reg80; reg95=reg97*reg130; T reg154=reg53*reg130;
    T reg155=reg90*reg80; reg115=reg118+reg115; reg118=reg89*reg80; T reg156=reg86*reg77; T reg157=reg19*reg11;
    T reg158=reg53*reg110; reg36=reg134+reg36; reg134=reg87*reg108; reg75=reg137+reg75; reg137=reg157*reg1;
    T reg159=elem.pos(1)[0]-elem.pos(0)[0]; T reg160=reg85*reg149; T reg161=reg43*reg54; reg83=reg124+reg83; reg98=reg136+reg98;
    reg124=reg30*reg155; reg95=reg125+reg95; reg125=reg85*reg139; reg145=reg63+reg145; reg63=reg62*reg114;
    reg148=reg64+reg148; reg64=reg35*reg149; reg154=reg99+reg154; reg82=reg104+reg82; reg99=reg102*reg77;
    reg142=reg55+reg142; reg55=reg128*reg139; reg104=reg65*reg114; reg144=reg100+reg144; reg100=reg43*reg155;
    reg136=reg30*reg54; reg146=reg138+reg146; reg138=reg76*reg78; T reg162=reg93*reg116; T reg163=reg86*reg108;
    T reg164=reg35*reg139; reg158=reg107+reg158; reg61=reg85*reg61; reg47=reg128*reg47; reg143=reg141+reg143;
    reg107=reg120*reg80; reg141=reg150*reg119; T reg165=reg151*reg119; T reg166=reg135*reg119; reg25=reg129*reg25;
    reg6=reg17+reg6; reg17=reg14*reg11; T reg167=reg31*reg155; reg111=reg123+reg111; reg156=reg115+reg156;
    reg41=reg52/reg41; reg52=reg31*reg118; reg115=reg60*reg16; reg94=reg105+reg94; reg105=reg29*reg24;
    reg152=reg133+reg152; reg147=reg153+reg147; reg123=reg31*reg54; reg126=reg132+reg126; reg73=reg92+reg73;
    reg92=reg30*reg118; reg132=reg128*reg149; reg26=reg26*reg60; reg133=reg120*reg121; reg153=reg23*reg11;
    T reg168=elem.pos(1)[1]-elem.pos(0)[1]; reg27=reg27*reg29; T reg169=elem.pos(2)[1]-elem.pos(0)[1]; reg58=reg140+reg58; reg39=reg40+reg39;
    reg12=reg129*reg12; reg40=elem.pos(2)[0]-elem.pos(0)[0]; reg48=reg127+reg48; reg127=reg69*reg38; reg131=reg68+reg131;
    reg161=reg83+reg161; reg68=reg135*reg166; reg83=reg43*reg118; reg99=reg82+reg99; reg167=reg111+reg167;
    reg82=reg102*reg108; reg162=reg138+reg162; reg109=reg109*reg29; reg111=reg26*reg165; reg138=reg135*reg165;
    reg100=reg144+reg100; reg140=reg31*reg107; reg163=reg143+reg163; reg143=reg153*reg0; reg144=(*f.m).alpha_1*reg70;
    T reg170=reg169*reg159; T reg171=reg43*reg114; reg133=reg147+reg133; reg104=reg36+reg104; reg36=reg97*reg137;
    reg147=reg40*reg168; reg12=reg39+reg12; reg125=reg148+reg125; reg127=reg48+reg127; reg47=reg61+reg47;
    reg64=reg154+reg64; reg38=reg129*reg38; reg25=reg6+reg25; reg6=reg17*reg119; reg39=reg115*reg41;
    reg48=reg105*reg41; reg61=reg152*reg41; reg148=reg42*(*f.m).alpha_2; reg160=reg95+reg160; reg95=(*f.m).alpha_1*reg71;
    reg154=reg28*(*f.m).alpha_2; reg52=reg156+reg52; reg156=reg26*reg141; reg55=reg142+reg55; reg142=reg28*reg58;
    T reg172=reg94*reg71; T reg173=reg42*reg58; T reg174=reg27*reg166; reg136=reg146+reg136; reg164=reg158+reg164;
    reg146=reg60*reg11; reg158=reg94*reg70; T reg175=reg30*reg107; reg134=reg75+reg134; reg75=reg14*reg18;
    T reg176=reg26*reg166; reg49=reg49*reg60; T reg177=reg27*reg141; reg92=reg126+reg92; reg63=reg145+reg63;
    reg126=reg103*reg137; reg19=reg19*reg14; reg124=reg98+reg124; reg98=reg27*reg165; reg123=reg131+reg123;
    reg132=reg73+reg132; reg73=reg25*reg71; reg131=reg75*reg41; reg145=reg25*reg70; reg147=reg170-reg147;
    reg170=(*f.m).alpha_3*reg86; T reg178=reg28*reg12; reg144=reg148+reg144; reg148=reg87*reg127; reg142=reg172+reg142;
    reg19=reg146+reg19; reg11=reg29*reg11; reg146=reg109*reg39; reg177=reg92+reg177; reg171=reg133+reg171;
    reg98=reg124+reg98; reg92=reg109*reg48; reg36=reg104+reg36; reg104=reg85*reg143; reg124=reg53*reg137;
    reg175=reg134+reg175; reg133=reg27*reg6; reg134=reg42*reg12; reg174=reg136+reg174; reg136=reg109*reg61;
    reg38=reg47+reg38; reg47=reg53*reg80; reg68=reg161+reg68; reg16=reg14*reg16; reg156=reg52+reg156;
    reg52=reg49*reg39; reg161=reg64*reg55; reg172=reg86*reg127; reg173=reg158+reg173; reg158=reg160*reg164;
    T reg179=reg49*reg48; reg176=reg123+reg176; reg123=reg43*reg107; reg82=reg162+reg82; reg138=reg100+reg138;
    reg100=reg152*reg48; reg162=reg49*reg61; reg111=reg167+reg111; reg167=reg132*reg164; reg23=reg23*reg14;
    reg154=reg95+reg154; reg95=(*f.m).alpha_3*reg87; reg140=reg163+reg140; reg163=reg26*reg6; T reg180=(*f.m).alpha_1*reg81;
    T reg181=reg46*(*f.m).alpha_2; reg83=reg99+reg83; reg99=reg135*reg141; T reg182=reg64*reg125; reg126=reg63+reg126;
    reg63=reg93*reg58; T reg183=reg128*reg143; T reg184=reg76*reg94; T reg185=reg60*reg18; T reg186=reg152*reg61;
    T reg187=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; T reg188=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg24=reg14*reg24; T reg189=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; reg124=reg171+reg124;
    reg170=reg144+reg170; reg95=reg154+reg95; reg23=reg11+reg23; reg181=reg180+reg181; reg11=(*f.m).alpha_3*reg117;
    reg144=(*f.m).alpha_1*reg76; reg154=reg93*(*f.m).alpha_2; reg168=reg168/reg147; reg40=reg40/reg147; reg16=reg185+reg16;
    reg167=reg161-reg167; reg158=reg182-reg158; reg161=reg125*reg132; reg171=reg160*reg55; reg180=reg29*reg18;
    reg182=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg185=reg35*reg143; reg63=reg184+reg63; reg184=reg86*reg38; reg134=reg145+reg134;
    reg136=reg174+reg136; reg148=reg142+reg148; reg142=reg30*reg47; reg186=reg68+reg186; reg178=reg73+reg178;
    reg68=reg87*reg38; reg52=reg156+reg52; reg73=reg135*reg6; reg123=reg82+reg123; reg100=reg138+reg100;
    reg179=reg111+reg179; reg82=reg31*reg47; reg172=reg173+reg172; reg111=reg152*reg39; reg162=reg176+reg162;
    reg99=reg83+reg99; reg183=reg126+reg183; reg163=reg140+reg163; reg83=reg49*reg131; reg104=reg36+reg104;
    reg169=reg169/reg147; reg159=reg159/reg147; reg80=reg35*reg80; reg36=reg19*reg119; reg146=reg177+reg146;
    reg92=reg98+reg92; reg133=reg175+reg133; reg93=reg93*reg12; reg98=reg102*reg127; reg126=reg109*reg131;
    reg76=reg76*reg25; reg24=reg180+reg24; reg184=reg134+reg184; reg31=reg31*reg80; reg134=reg92*reg95;
    reg138=reg146*reg170; reg126=reg133+reg126; reg133=reg179*reg95; reg82=reg172+reg82; reg140=reg168*reg188;
    reg145=reg26*reg36; reg156=reg92*reg186; reg172=reg52*reg170; reg173=reg179*reg186; reg142=reg148+reg142;
    reg148=reg136*reg100; reg30=reg30*reg80; reg174=reg162*reg100; reg68=reg178+reg68; reg175=reg159*reg187;
    reg176=reg40*reg189; reg177=reg169*reg182; reg178=reg27*reg36; reg180=reg102*reg38; T reg190=reg45*(*f.m).alpha_2;
    T reg191=(*f.m).alpha_1*reg157; reg93=reg76+reg93; reg76=reg158*reg183; T reg192=reg152*reg131; T reg193=(*f.m).alpha_3*reg102;
    T reg194=reg16*reg41; reg73=reg123+reg73; reg119=reg23*reg119; reg154=reg144+reg154; reg11=reg181+reg11;
    reg98=reg63+reg98; reg171=reg161-reg171; reg185=reg124+reg185; reg111=reg99+reg111; reg63=reg167*reg104;
    reg60=reg14*reg60; reg83=reg163+reg83; reg99=reg43*reg47; reg123=reg135*reg36; reg140=reg177-reg140;
    reg176=reg175-reg176; reg180=reg93+reg180; reg182=reg40*reg182; reg93=reg43*reg80; reg188=reg159*reg188;
    reg124=reg26*reg119; reg187=reg168*reg187; reg144=reg132*reg185; reg41=reg24*reg41; reg189=reg169*reg189;
    reg133=reg172+reg133; reg83=reg83*reg11; reg161=reg100*reg95; reg163=reg111*reg170; reg134=reg138+reg134;
    reg126=reg126*reg11; reg138=reg64*reg183; reg172=reg160*reg185; reg29=reg14*reg29; reg175=reg64*reg104;
    reg177=reg171*reg185; reg145=reg82+reg145; reg82=reg49*reg194; reg181=reg27*reg119; reg76=reg63-reg76;
    reg30=reg68+reg30; reg148=reg156-reg148; reg192=reg73+reg192; reg174=reg173-reg174; reg63=reg179*reg136;
    reg99=reg98+reg99; reg68=reg10*(*f.m).alpha_2; reg73=reg162*reg92; reg98=(*f.m).alpha_1*reg153; reg156=reg109*reg194;
    reg178=reg142+reg178; reg142=(*f.m).alpha_3*reg60; reg31=reg184+reg31; reg193=reg154+reg193; reg190=reg191+reg190;
    reg140=reg176+reg140; reg154=reg136*reg193; reg161=reg163+reg161; reg192=reg192*reg11; reg163=reg170*(*f.m).deltaT;
    reg173=reg95*(*f.m).deltaT; reg182=reg188-reg182; elem.epsilon[0][1]=reg182; reg187=reg189-reg187; elem.epsilon[0][0]=reg187;
    reg176=reg86*reg57; reg184=reg109*reg41; reg181=reg30+reg181; reg156=reg178+reg156; reg30=reg49*reg41;
    reg178=reg146*reg174; reg188=reg52*reg148; reg124=reg31+reg124; reg82=reg145+reg82; reg142=reg190+reg142;
    reg31=reg160*reg183; reg68=reg98+reg68; reg98=reg104*reg132; reg145=(*f.m).alpha_3*reg29; reg189=reg125*reg185;
    reg177=reg76+reg177; reg172=reg175-reg172; reg76=reg183*reg164; reg175=reg104*reg164; reg144=reg138-reg144;
    reg138=reg55*reg185; reg126=reg134+reg126; reg83=reg133+reg83; reg133=reg162*reg193; reg123=reg99+reg123;
    reg99=reg152*reg194; reg134=reg87*reg67; reg190=reg86*reg7; reg93=reg180+reg93; reg180=reg87*reg50;
    reg73=reg63-reg73; reg63=reg135*reg119; reg191=reg193*(*f.m).deltaT; reg184=reg181+reg184; reg180=reg176+reg180;
    reg176=reg117*reg121; reg181=reg152*reg41; T reg195=reg182-reg173; reg145=reg68+reg145; reg68=reg187-reg163;
    reg99=reg123+reg99; reg123=reg117*reg13; reg134=reg190+reg134; reg30=reg124+reg30; reg124=reg87*reg56;
    reg190=reg86*reg3; reg140=0.5*reg140; elem.epsilon[0][2]=reg140; T reg196=reg162*reg111; reg63=reg93+reg63;
    reg82=reg82*reg142; reg133=reg83+reg133; reg192=reg161+reg192; reg83=reg186*reg193; reg93=reg146*reg186;
    reg161=reg111*reg73; reg178=reg188-reg178; reg188=reg136*reg111; T reg197=reg52*reg186; T reg198=reg125*reg183;
    reg144=reg144/reg177; reg31=reg98-reg31; reg98=reg104*reg55; reg154=reg126+reg154; reg138=reg76-reg138;
    reg158=reg158/reg177; reg172=reg172/reg177; reg156=reg156*reg142; reg167=reg167/reg177; reg189=reg175-reg189;
    reg76=reg88*reg114; reg30=reg30*reg145; reg126=reg52*reg136; reg175=reg88*reg84; T reg199=reg162*reg146;
    reg156=reg154+reg156; reg154=reg146*reg100; reg188=reg93-reg188; reg123=reg134+reg123; reg93=reg140-reg191;
    reg134=reg144*reg68; T reg200=reg172*reg195; T reg201=reg92*reg111; reg184=reg184*reg145; T reg202=reg158*reg195;
    T reg203=reg167*reg68; T reg204=reg52*reg100; T reg205=reg117*reg66; reg124=reg190+reg124; reg31=reg31/reg177;
    reg181=reg63+reg181; reg196=reg197-reg196; reg63=reg179*reg111; reg171=reg171/reg177; reg198=reg98-reg198;
    reg189=reg189/reg177; reg161=reg178+reg161; reg99=reg99*reg142; reg138=reg138/reg177; reg83=reg192+reg83;
    reg82=reg133+reg82; reg176=reg180+reg176; reg68=reg138*reg68; reg195=reg189*reg195; reg98=reg179*reg146;
    reg133=reg69*reg110; reg175=reg123+reg175; reg174=reg174/reg161; reg199=reg126-reg199; reg123=reg88*reg74;
    reg205=reg124+reg205; reg201=reg154-reg201; reg196=reg196/reg161; reg63=reg204-reg63; reg177=reg198/reg177;
    reg124=reg52*reg92; reg188=reg188/reg161; reg148=reg148/reg161; reg76=reg176+reg76; reg126=reg69*reg137;
    reg154=reg171*reg93; reg202=reg203-reg202; reg134=reg200-reg134; reg30=reg82+reg30; reg184=reg156+reg184;
    reg82=reg31*reg93; reg99=reg83+reg99; reg181=reg181*reg145; reg148=reg148*reg30; reg174=reg174*reg184;
    reg188=reg188*reg30; reg196=reg196*reg184; reg201=reg201/reg161; reg181=reg99+reg181; reg195=reg68-reg195;
    reg93=reg177*reg93; reg68=reg69*reg130; reg123=reg205+reg123; reg63=reg63/reg161; reg73=reg73/reg161;
    reg199=reg199/reg161; reg82=reg134-reg82; reg133=reg175+reg133; reg83=reg129*reg143; reg126=reg76+reg126;
    reg76=reg129*reg139; reg98=reg124-reg98; reg202=reg154+reg202; reg174=reg148-reg174; reg199=reg199*reg181;
    reg188=reg196-reg188; reg30=reg201*reg30; reg99=reg82*reg55; reg124=reg202*reg183; reg195=reg93+reg195;
    reg83=reg126+reg83; reg184=reg63*reg184; reg63=reg151*reg67; reg73=reg73*reg181; reg93=reg202*reg104;
    reg161=reg98/reg161; reg76=reg133+reg76; reg98=reg82*reg125; reg126=reg150*reg7; reg133=reg151*reg50;
    reg134=reg150*reg57; reg148=reg129*reg149; reg68=reg123+reg68; reg133=reg134+reg133; reg99=reg124+reg99;
    reg123=reg195*reg132; reg148=reg68+reg148; reg174=reg73+reg174; reg181=reg161*reg181; reg184=reg30-reg184;
    reg30=PNODE(2).dep[0]-PNODE(0).dep[0]; reg67=reg105*reg67; reg7=reg115*reg7; reg68=PNODE(1).dep[1]-PNODE(0).dep[1]; reg73=PNODE(2).dep[1]-PNODE(0).dep[1];
    reg124=PNODE(1).dep[0]-PNODE(0).dep[0]; reg50=reg105*reg50; reg57=reg115*reg57; reg134=1-(*f.m).resolution; reg63=reg126+reg63;
    reg126=reg17*reg13; reg98=reg93+reg98; reg199=reg188-reg199; reg93=reg195*reg160; reg154=reg83*reg202;
    reg156=reg76*reg82; reg161=reg150*reg3; reg175=reg17*reg121; reg176=reg151*reg56; reg178=reg68*reg169;
    reg121=reg75*reg121; reg180=reg135*reg84; reg95=(*f.m).resolution*reg95; reg50=reg57+reg50; reg126=reg63+reg126;
    reg56=reg105*reg56; reg174=reg134*reg174; reg57=reg73*reg168; reg63=reg124*reg40; reg123=reg99+reg123;
    reg67=reg7+reg67; reg7=reg30*reg159; reg99=reg135*reg114; reg199=reg134*reg199; reg13=reg75*reg13;
    reg184=reg181+reg184; reg175=reg133+reg175; reg176=reg161+reg176; reg93=reg98+reg93; reg3=reg115*reg3;
    reg170=(*f.m).resolution*reg170; reg98=reg17*reg66; reg156=reg154+reg156; reg133=reg148*reg195; reg30=reg168*reg30;
    reg73=reg73*reg159; reg154=reg135*reg74; reg180=reg126+reg180; reg63=reg7-reg63; reg121=reg50+reg121;
    reg114=reg152*reg114; reg7=reg19*reg137; reg99=reg175+reg99; reg57=reg178-reg57; reg68=reg68*reg40;
    reg124=reg124*reg169; reg84=reg152*reg84; reg13=reg67+reg13; reg50=reg19*reg110; reg98=reg176+reg98;
    reg123=reg173+reg123; reg133=reg156+reg133; reg56=reg3+reg56; reg93=reg163+reg93; reg199=reg95+reg199;
    reg174=reg170+reg174; reg184=reg134*reg184; reg193=(*f.m).resolution*reg193; reg11=reg11*(*f.m).deltaT; reg66=reg75*reg66;
    reg199=reg199*(*f.m).deltaT; reg3=reg19*reg130; reg154=reg98+reg154; reg66=reg56+reg66; reg74=reg152*reg74;
    reg30=reg124-reg30; reg7=reg99+reg7; reg56=reg23*reg143; reg50=reg180+reg50; reg67=reg23*reg139;
    reg133=reg11+reg133; reg95=reg93+reg123; reg52=reg52*reg134; reg179=reg179*reg134; reg57=reg63+reg57;
    reg68=reg73-reg68; reg146=reg146*reg134; reg92=reg92*reg134; reg111=reg111*reg134; reg100=reg100*reg134;
    reg63=(*f.m).resolution*reg144; reg73=(*f.m).resolution*reg172; reg98=(*f.m).resolution*reg138; reg99=(*f.m).resolution*reg189; reg110=reg16*reg110;
    reg84=reg13+reg84; reg184=reg193+reg184; reg137=reg16*reg137; reg114=reg121+reg114; reg174=reg174*(*f.m).deltaT;
    reg13=(*f.m).resolution*reg167; reg121=(*f.m).resolution*reg158; reg74=reg66+reg74; reg63=reg146-reg63; reg73=reg92+reg73;
    reg66=reg23*reg149; reg56=reg7+reg56; reg98=reg111+reg98; reg139=reg24*reg139; reg110=reg84+reg110;
    reg99=reg100-reg99; reg7=reg68-reg199; reg84=(*f.m).resolution*reg171; reg92=reg30-reg174; reg143=reg24*reg143;
    reg137=reg114+reg137; reg184=reg184*(*f.m).deltaT; reg162=reg162*reg134; reg57=0.5*reg57; reg95=reg133+reg95;
    reg136=reg136*reg134; reg134=reg186*reg134; reg100=(*f.m).resolution*reg31; reg111=reg202*reg185; reg67=reg50+reg67;
    reg50=reg82*reg164; reg3=reg154+reg3; reg130=reg16*reg130; reg114=(*f.m).resolution*reg177; reg52=reg13+reg52;
    reg121=reg179-reg121; reg13=reg121*reg7; reg124=reg63*reg92; reg126=reg52*reg92; reg146=reg73*reg7;
    reg154=reg98*reg92; reg156=reg99*reg7; reg161=reg67*reg82; reg170=reg56*reg202; reg130=reg74+reg130;
    reg149=reg24*reg149; reg95=reg95/3; reg50=reg111+reg50; reg74=reg64*reg195; reg162=reg84+reg162;
    reg100=reg136-reg100; reg139=reg110+reg139; reg114=reg134+reg114; reg143=reg137+reg143; reg84=reg57-reg184;
    reg66=reg3+reg66; reg3=reg100*reg84; reg202=reg202*reg143; reg82=reg82*reg139; reg146=reg124+reg146;
    reg156=reg154+reg156; reg110=reg114*reg84; reg161=reg170+reg161; reg111=reg162*reg84; reg13=reg126+reg13;
    reg124=reg195*reg66; reg149=reg130+reg149; reg123=reg123-reg95; reg93=reg93-reg95; reg74=reg50+reg74;
    reg124=reg161+reg124; reg111=reg13+reg111; reg82=reg202+reg82; reg195=reg149*reg195; reg3=reg146+reg3;
    reg74=reg191+reg74; reg142=reg142*(*f.m).deltaT; reg110=reg156+reg110; reg95=reg133-reg95; reg93=pow(reg93,2);
    reg123=pow(reg123,2); reg195=reg82+reg195; reg95=pow(reg95,2); reg7=reg7*reg3; reg92=reg92*reg111;
    reg124=reg142+reg124; reg123=reg93+reg123; reg174=reg187-reg174; reg145=reg145*(*f.m).deltaT; reg199=reg182-reg199;
    reg110=2*reg110; reg13=2*reg74; reg52=reg52*reg174; reg121=reg121*reg199; reg195=reg145+reg195;
    reg95=reg123+reg95; reg184=reg140-reg184; reg13=reg74*reg13; reg50=2*reg124; reg84=reg84*reg110;
    reg63=reg63*reg174; reg73=reg73*reg199; reg92=reg7+reg92; reg162=reg162*reg184; reg73=reg63+reg73;
    reg100=reg100*reg184; reg7=2*reg195; reg121=reg52+reg121; reg174=reg98*reg174; reg199=reg99*reg199;
    reg13=reg95+reg13; reg50=reg124*reg50; reg84=reg92+reg84; reg7=reg195*reg7; reg84=reg147*reg84;
    reg162=reg121+reg162; elem.sigma[0][0]=reg162; reg100=reg73+reg100; elem.sigma[0][1]=reg100; reg199=reg174+reg199;
    reg184=reg114*reg184; reg50=reg13+reg50; reg13=0.16666666666666665741*reg84; reg84=0.33333333333333331483*reg84; reg184=reg199+reg184;
    elem.sigma[0][2]=reg184; reg52=reg70*reg162; reg7=reg50+reg7; reg50=reg71*reg100; reg63=reg42*reg162;
    reg73=reg28*reg100; reg162=reg89*reg162; reg100=reg90*reg100; reg74=reg62*reg184; reg82=reg65*reg184;
    reg50=reg52+reg50; reg73=reg63+reg73; reg7=1.5*reg7; reg84=reg13+reg84; reg184=reg43*reg184;
    reg100=reg162+reg100; elem.sigma_local[0][1]=reg73+reg74; elem.ener=reg84/2; elem.sigma_local[0][2]=reg100+reg184; elem.sigma_von_mises=pow(reg7,0.5);
    elem.sigma_local[0][0]=reg50+reg82;
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
    T reg0=2*(*f.m).shear_modulus_13; T reg1=2*(*f.m).shear_modulus_23; reg1=1.0/reg1; reg0=1.0/reg0; T reg2=2*(*f.m).shear_modulus_12;
    T reg3=reg0*reg1; reg2=1.0/reg2; T reg4=pow((*f.m).v2[1],2); T reg5=pow((*f.m).v2[0],2); T reg6=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg7=pow((*f.m).v1[1],2); T reg8=pow((*f.m).v1[0],2); T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg11=reg2*reg3;
    T reg12=reg10*reg11; T reg13=pow((*f.m).v2[2],2); T reg14=reg11*reg9; reg4=reg5+reg4; reg5=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1;
    T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=reg11*reg6; reg7=reg8+reg7; reg8=pow((*f.m).v1[2],2); reg13=reg4+reg13;
    reg8=reg7+reg8; reg4=reg10*reg16; reg7=reg5*reg14; T reg17=reg15*reg14; T reg18=reg10*reg12;
    T reg19=1.0/(*f.m).elastic_modulus_1; reg18=reg17-reg18; reg4=reg7+reg4; reg17=reg5*reg12; T reg20=reg15*reg16;
    reg13=pow(reg13,0.5); reg8=pow(reg8,0.5); T reg21=(*f.m).v2[2]/reg13; T reg22=(*f.m).v2[1]/reg13; T reg23=reg19*reg18;
    T reg24=reg17+reg20; T reg25=(*f.m).v1[2]/reg8; T reg26=reg5*reg4; T reg27=(*f.m).v1[1]/reg8; T reg28=reg2*reg1;
    T reg29=reg10*reg3; T reg30=reg3*reg6; T reg31=reg5*reg11; T reg32=reg16*reg6; reg14=reg19*reg14;
    reg26=reg23-reg26; reg23=reg25*reg22; T reg33=reg27*reg21; reg3=reg3*reg9; reg8=(*f.m).v1[0]/reg8;
    T reg34=reg12*reg6; reg11=reg15*reg11; T reg35=reg24*reg6; reg13=(*f.m).v2[0]/reg13; reg34=reg7+reg34;
    reg32=reg14-reg32; reg7=reg11*reg6; reg12=reg19*reg12; reg14=reg31*reg6; T reg36=reg10*reg29;
    T reg37=reg5*reg3; T reg38=reg10*reg30; reg3=reg15*reg3; T reg39=reg33-reg23; T reg40=reg28*reg6;
    T reg41=reg25*reg13; T reg42=reg8*reg21; T reg43=reg28*reg9; T reg44=reg2*reg0; reg35=reg26-reg35;
    reg26=2*reg8; T reg45=2*reg13; reg16=reg5*reg16; reg28=reg10*reg28; reg11=reg19*reg11;
    T reg46=reg15*reg43; reg43=reg5*reg43; reg18=reg18/reg35; T reg47=reg10*reg28; reg34=reg34/reg35;
    reg7=reg17+reg7; reg4=reg4/reg35; reg32=reg32/reg35; reg14=reg12+reg14; reg16=reg12+reg16;
    reg31=reg5*reg31; reg29=reg5*reg29; reg12=reg10*reg44; reg30=reg15*reg30; T reg48=reg44*reg6;
    reg38=reg37+reg38; reg36=reg3-reg36; reg3=reg41-reg42; reg37=reg8*reg22; T reg49=reg27*reg13;
    T reg50=pow(reg13,2); T reg51=2*reg39; reg44=reg44*reg9; T reg52=pow(reg8,2); T reg53=reg45*reg22;
    T reg54=reg10*reg40; T reg55=reg27*reg26; T reg56=pow(reg22,2); T reg57=pow(reg27,2); T reg58=reg18*reg52;
    T reg59=reg50*reg4; T reg60=reg50*reg32; T reg61=reg34*reg52; T reg62=pow(reg21,2); reg31=reg11-reg31;
    reg11=reg53*reg4; T reg63=reg55*reg18; T reg64=pow(reg25,2); reg24=reg24/reg35; reg14=reg14/reg35;
    reg16=reg16/reg35; T reg65=reg37-reg49; T reg66=reg56*reg4; T reg67=reg18*reg57; reg47=reg46-reg47;
    reg36=reg19*reg36; reg54=reg43+reg54; reg43=reg29+reg30; reg46=reg53*reg32; T reg68=reg34*reg57;
    T reg69=reg56*reg32; reg38=reg5*reg38; T reg70=reg51*reg3; T reg71=reg55*reg34; reg40=reg15*reg40;
    T reg72=pow(reg3,2); T reg73=reg15*reg44; reg44=reg5*reg44; T reg74=reg10*reg12; T reg75=pow(reg39,2);
    reg28=reg5*reg28; T reg76=reg10*reg48; reg7=reg7/reg35; T reg77=reg72*reg16; reg69=reg68+reg69;
    reg38=reg36-reg38; reg36=reg75*reg16; reg68=reg34*reg64; reg60=reg61+reg60; reg47=reg19*reg47;
    reg54=reg5*reg54; reg61=reg28+reg40; reg46=reg71+reg46; reg71=reg70*reg16; T reg78=reg7*reg52;
    T reg79=reg50*reg14; T reg80=reg7*reg57; T reg81=reg56*reg14; T reg82=reg55*reg7; T reg83=reg62*reg32;
    reg74=reg73-reg74; reg73=reg53*reg14; reg76=reg44+reg76; reg66=reg67+reg66; reg44=reg72*reg24;
    reg67=reg18*reg64; T reg84=reg62*reg4; reg11=reg63+reg11; reg63=reg70*reg24; T reg85=pow(reg65,2);
    reg43=reg43*reg6; reg48=reg15*reg48; reg12=reg5*reg12; reg31=reg31/reg35; T reg86=reg75*reg24;
    reg59=reg58+reg59; reg58=reg45*reg21; reg84=reg67+reg84; reg83=reg68+reg83; reg67=reg70*reg31;
    reg73=reg82+reg73; reg68=reg16*reg85; reg74=reg19*reg74; reg82=reg72*reg31; T reg87=reg62*reg14;
    T reg88=2*reg27; T reg89=reg26*reg25; reg76=reg5*reg76; T reg90=reg24*reg85; reg71=reg46+reg71;
    reg46=reg12+reg48; reg61=reg61*reg6; reg77=reg69+reg77; reg69=2*reg22; reg63=reg11+reg63;
    reg54=reg47-reg54; reg11=reg13*reg22; reg44=reg66+reg44; reg36=reg60+reg36; reg81=reg80+reg81;
    reg86=reg59+reg86; reg43=reg38-reg43; reg38=reg8*reg27; reg47=reg75*reg31; reg79=reg78+reg79;
    reg59=reg7*reg64; reg60=reg56*reg36; reg66=reg11*reg71; reg78=reg38*reg63; reg80=reg26*reg13;
    T reg91=reg50*reg77; T reg92=reg44*reg52; T reg93=reg44*reg57; T reg94=reg56*reg77; T reg95=reg86*reg57;
    T reg96=reg58*reg32; T reg97=reg89*reg34; T reg98=reg50*reg71; T reg99=reg63*reg52; reg68=reg83+reg68;
    reg83=reg88*reg25; T reg100=reg69*reg21; reg51=reg51*reg65; T reg101=2*reg3; T reg102=reg8*reg13;
    T reg103=reg27*reg22; reg37=reg49+reg37; reg49=reg8*reg3; T reg104=reg27*reg39; reg71=reg56*reg71;
    reg63=reg63*reg57; T reg105=reg50*reg36; reg61=reg54-reg61; reg67=reg73+reg67; reg54=reg86*reg52;
    reg76=reg74-reg76; reg44=reg38*reg44; reg77=reg11*reg77; reg46=reg46*reg6; reg73=reg31*reg85;
    reg87=reg59+reg87; reg59=reg3*reg39; reg43=reg43/reg35; reg47=reg79+reg47; reg74=reg58*reg4;
    reg79=reg89*reg18; T reg106=reg88*reg22; reg82=reg81+reg82; reg81=2*reg25; reg90=reg84+reg90;
    reg84=reg56*reg15; T reg107=reg27*reg3; T reg108=reg102*reg43; T reg109=reg8*reg39; reg101=reg101*reg65;
    T reg110=reg50*reg5; T reg111=reg25*reg21; T reg112=reg103*reg43; T reg113=reg37*reg43; reg105=reg54+reg105;
    reg54=reg75*reg47; T reg114=reg58*reg14; T reg115=reg89*reg7; T reg116=reg50*reg19; T reg117=reg56*reg5;
    reg73=reg87+reg73; reg96=reg97+reg96; reg87=reg51*reg16; reg97=reg75*reg82; reg34=reg83*reg34;
    reg32=reg100*reg32; reg91=reg92+reg91; reg92=reg80*reg19; T reg118=reg106*reg5; T reg119=reg72*reg67;
    T reg120=reg59*reg67; reg86=reg38*reg86; reg36=reg11*reg36; reg98=reg99+reg98; reg67=reg75*reg67;
    reg60=reg95+reg60; reg66=reg78+reg66; reg78=reg72*reg47; reg46=reg76-reg46; reg77=reg44+reg77;
    reg44=reg59*reg82; reg94=reg93+reg94; reg4=reg100*reg4; reg18=reg83*reg18; reg82=reg72*reg82;
    reg76=reg51*reg24; reg81=reg81*reg21; reg74=reg79+reg74; reg49=reg104+reg49; reg79=reg80*reg5;
    reg93=reg90*reg57; reg95=reg106*reg15; reg99=reg56*reg68; reg104=reg39*reg22; T reg121=reg3*reg13;
    T reg122=reg50*reg68; T reg123=reg90*reg52; reg71=reg63+reg71; reg61=reg61/reg35; reg47=reg59*reg47;
    reg54=reg105+reg54; reg63=reg72*reg73; reg82=reg94+reg82; reg99=reg93+reg99; reg93=reg49*reg61;
    reg94=reg107*reg61; reg105=reg80*reg108; reg36=reg86+reg36; reg86=reg106*reg113; T reg124=reg109*reg61;
    reg119=reg71+reg119; reg71=reg106*reg112; reg44=reg77+reg44; reg110=reg84-reg110; reg77=reg62*reg10;
    reg84=reg81*reg10; reg79=reg95-reg79; reg95=reg50*reg6; T reg125=reg56*reg10; T reg126=reg80*reg6;
    T reg127=reg106*reg10; T reg128=reg37*reg113; reg97=reg91+reg97; reg91=reg80*reg112; reg122=reg123+reg122;
    reg123=reg75*reg73; reg120=reg66+reg120; reg118=reg92-reg118; reg66=reg81*reg6; reg67=reg98+reg67;
    reg113=reg80*reg113; reg78=reg60+reg78; reg60=reg106*reg108; reg68=reg11*reg68; reg90=reg38*reg90;
    reg117=reg116-reg117; reg92=reg62*reg6; reg112=reg37*reg112; reg114=reg115+reg114; reg98=reg51*reg31;
    reg7=reg83*reg7; reg14=reg100*reg14; reg24=reg101*reg24; reg4=reg18+reg4; reg76=reg74+reg76;
    reg35=reg46/reg35; reg26=reg26*reg39; reg18=reg3*reg22; reg46=reg65*reg25; reg16=reg101*reg16;
    reg88=reg88*reg3; reg32=reg34+reg32; reg87=reg96+reg87; reg121=reg104+reg121; reg34=reg111*reg43;
    reg74=reg39*reg13; reg125=reg95+reg125; reg95=reg62*reg9; reg96=reg56*(*f.m).alpha_2; reg104=reg65*reg21;
    reg81=reg81*reg9; reg127=reg126+reg127; reg47=reg36+reg47; reg108=reg37*reg108; reg36=reg50*reg87;
    reg115=reg76*reg52; reg116=(*f.m).alpha_1*reg52; reg126=reg26*reg93; reg113=reg67+reg113; reg67=reg5*reg52;
    reg84=reg79-reg84; reg79=reg106*reg34; reg63=reg99+reg63; reg86=reg119+reg86; reg99=reg88*reg93;
    reg119=reg8*reg65; T reg129=reg76*reg57; T reg130=reg56*reg87; T reg131=reg88*reg94; reg71=reg82+reg71;
    reg77=reg110-reg77; reg24=reg4+reg24; reg4=reg50*(*f.m).alpha_2; reg82=reg88*reg124; reg60=reg78+reg60;
    reg78=reg39*reg25; reg110=(*f.m).alpha_1*reg57; reg41=reg42+reg41; reg15=reg15*reg57; reg68=reg90+reg68;
    reg73=reg59*reg73; reg42=reg26*reg94; reg91=reg97+reg91; reg16=reg32+reg16; reg128=reg120+reg128;
    reg93=reg49*reg93; reg32=reg26*reg124; reg105=reg54+reg105; reg98=reg114+reg98; reg14=reg7+reg14;
    reg5=reg5*reg57; reg19=reg19*reg52; reg7=reg121*reg35; reg54=reg18*reg35; reg90=reg74*reg35;
    reg31=reg101*reg31; reg97=reg46*reg61; reg112=reg44+reg112; reg45=reg45*reg39; reg69=reg69*reg3;
    reg94=reg49*reg94; reg123=reg122+reg123; reg66=reg118-reg66; reg44=reg80*reg34; reg92=reg117-reg92;
    reg124=reg49*reg124; reg114=reg121*reg7; reg76=reg38*reg76; reg87=reg11*reg87; reg117=reg50*reg92;
    reg118=reg69*reg7; reg120=reg103*reg84; reg99=reg86+reg99; reg86=reg102*reg66; reg122=reg56*reg77;
    T reg132=reg39*reg21; T reg133=reg10*reg57; reg116=reg4+reg116; reg4=reg52*reg6; T reg134=reg50*reg66;
    reg127=reg81-reg127; reg81=reg121*reg54; T reg135=reg102*reg92; T reg136=reg103*reg77; reg73=reg68+reg73;
    reg125=reg95-reg125; reg68=reg52*reg66; reg34=reg37*reg34; reg95=reg57*reg77; T reg137=reg65*reg13;
    T reg138=reg52*reg92; T reg139=reg57*reg84; T reg140=reg56*reg16; T reg141=reg24*reg57; T reg142=reg56*reg84;
    reg93=reg128+reg93; reg94=reg112+reg94; reg112=reg72*reg98; reg130=reg129+reg130; reg108=reg47+reg108;
    reg47=reg41*reg43; reg33=reg23+reg33; reg6=reg64*reg6; reg67=reg15-reg67; reg42=reg91+reg42;
    reg15=reg45*reg54; reg96=reg110+reg96; reg23=(*f.m).alpha_3*reg72; reg5=reg19-reg5; reg19=reg50*reg16;
    reg91=reg24*reg52; reg10=reg10*reg64; reg110=(*f.m).alpha_1*reg64; reg128=reg62*(*f.m).alpha_2; reg44=reg123+reg44;
    reg123=reg26*reg97; reg31=reg14+reg31; reg126=reg113+reg126; reg7=reg45*reg7; reg36=reg115+reg36;
    reg14=reg75*reg98; reg113=(*f.m).alpha_3*reg75; reg32=reg105+reg32; reg105=reg88*reg97; reg79=reg63+reg79;
    reg63=reg45*reg90; reg115=reg27*reg65; reg129=reg3*reg25; reg119=reg78+reg119; reg54=reg69*reg54;
    reg131=reg71+reg131; reg71=reg69*reg90; reg78=reg104*reg35; reg82=reg60+reg82; reg98=reg59*reg98;
    reg114=reg93+reg114; reg7=reg126+reg7; reg60=reg62*reg125; reg124=reg108+reg124; reg90=reg121*reg90;
    reg122=reg117+reg122; reg93=(*f.m).alpha_1*reg38; reg13=reg13*reg21; reg8=reg8*reg25; reg108=reg11*(*f.m).alpha_2;
    reg117=reg45*reg78; reg123=reg44+reg123; reg87=reg76+reg87; reg97=reg49*reg97; reg63=reg32+reg63;
    reg32=reg111*reg125; reg81=reg94+reg81; reg34=reg73+reg34; reg136=reg135+reg136; reg137=reg132+reg137;
    reg6=reg5-reg6; reg5=reg65*reg22; reg15=reg42+reg15; reg42=reg3*reg21; reg19=reg91+reg19;
    reg44=reg75*reg31; reg23=reg96+reg23; reg73=reg64*reg125; reg95=reg138+reg95; reg10=reg67-reg10;
    reg67=reg72*reg31; reg140=reg141+reg140; reg113=reg116+reg113; reg115=reg129+reg115; reg76=reg106*reg47;
    reg112=reg130+reg112; reg71=reg82+reg71; reg118=reg99+reg118; reg54=reg131+reg54; reg133=reg4+reg133;
    reg9=reg64*reg9; reg105=reg79+reg105; reg4=reg69*reg78; reg120=reg86+reg120; reg14=reg36+reg14;
    reg142=reg134+reg142; reg36=reg62*reg127; reg79=reg111*reg127; reg82=reg119*reg61; reg86=reg64*reg127;
    reg139=reg68+reg139; reg68=reg80*reg47; reg91=(*f.m).alpha_3*reg85; reg24=reg38*reg24; reg16=reg11*reg16;
    reg94=reg37*reg2; reg11=reg11*reg2; reg43=reg33*reg43; reg128=reg110+reg128; reg79=reg120+reg79;
    reg96=(*f.m).alpha_3*reg59; reg99=reg37*reg94; reg110=(*f.m).alpha_1*reg8; reg116=reg54*reg23; reg120=reg71*reg113;
    reg126=reg15*reg23; reg129=reg52*reg6; reg130=reg57*reg10; reg131=reg13*(*f.m).alpha_2; reg73=reg95+reg73;
    reg95=reg55*reg11; reg132=reg54*reg114; reg134=reg15*reg114; reg86=reg139+reg86; reg135=reg55*reg94;
    reg138=reg118*reg81; reg139=reg50*reg6; reg141=reg56*reg10; T reg143=reg7*reg81; reg91=reg128+reg91;
    reg60=reg122+reg60; reg122=reg53*reg11; reg128=reg63*reg113; reg108=reg93+reg108; reg78=reg121*reg78;
    reg97=reg34+reg97; reg39=reg65*reg39; reg25=reg27*reg25; reg90=reg124+reg90; reg21=reg22*reg21;
    reg36=reg142+reg36; reg22=reg53*reg94; reg27=reg41*reg0; reg13=reg13*reg0; reg106=reg106*reg43;
    reg67=reg140+reg67; reg34=reg88*reg82; reg76=reg112+reg76; reg133=reg9-reg133; reg4=reg105+reg4;
    reg80=reg80*reg43; reg44=reg19+reg44; reg9=reg26*reg82; reg68=reg14+reg68; reg14=reg37*reg11;
    reg32=reg136+reg32; reg117=reg123+reg117; reg19=reg137*reg35; reg61=reg115*reg61; reg16=reg24+reg16;
    reg31=reg59*reg31; reg5=reg42+reg5; reg98=reg87+reg98; reg47=reg37*reg47; reg34=reg76+reg34;
    reg130=reg129+reg130; reg64=reg64*reg133; reg24=reg69*reg19; reg131=reg110+reg131; reg39=(*f.m).alpha_3*reg39;
    reg42=(*f.m).alpha_1*reg25; reg59=reg21*(*f.m).alpha_2; reg76=reg90*reg113; reg82=reg49*reg82; reg106=reg67+reg106;
    reg88=reg88*reg61; reg2=reg38*reg2; reg21=reg21*reg1; reg38=reg33*reg1; reg67=reg41*reg13;
    reg14=reg32+reg14; reg26=reg26*reg61; reg80=reg44+reg80; reg117=reg117*reg91; reg32=reg45*reg19;
    reg9=reg68+reg9; reg31=reg16+reg31; reg43=reg37*reg43; reg116=reg120+reg116; reg35=reg5*reg35;
    reg4=reg4*reg91; reg96=reg108+reg96; reg16=reg15*reg118; reg143=reg134-reg143; reg126=reg128+reg126;
    reg44=reg7*reg54; reg62=reg62*reg133; reg141=reg139+reg141; reg3=reg65*reg3; reg138=reg132-reg138;
    reg65=reg58*reg13; reg68=reg89*reg27; reg135=reg86+reg135; reg86=reg41*reg27; reg99=reg79+reg99;
    reg22=reg36+reg22; reg36=reg58*reg27; reg122=reg60+reg122; reg47=reg98+reg47; reg60=reg81*reg23;
    reg78=reg97+reg78; reg79=reg102*reg6; reg87=reg89*reg13; reg95=reg73+reg95; reg73=reg103*reg10;
    reg93=reg118*reg96; reg44=reg16-reg44; reg4=reg116+reg4; reg16=reg7*reg96; reg117=reg126+reg117;
    reg97=reg71*reg143; reg60=reg76+reg60; reg76=reg63*reg138; reg78=reg78*reg91; reg39=reg131+reg39;
    reg0=reg8*reg0; reg59=reg42+reg59; reg3=(*f.m).alpha_3*reg3; reg32=reg9+reg32; reg26=reg80+reg26;
    reg45=reg45*reg35; reg82=reg47+reg82; reg67=reg14+reg67; reg19=reg121*reg19; reg8=reg83*reg21;
    reg69=reg69*reg35; reg68=reg135+reg68; reg9=reg83*reg38; reg62=reg141+reg62; reg14=reg53*reg2;
    reg88=reg106+reg88; reg24=reg34+reg24; reg43=reg31+reg43; reg61=reg49*reg61; reg31=reg33*reg21;
    reg34=reg33*reg38; reg86=reg99+reg86; reg65=reg122+reg65; reg42=reg100*reg21; reg64=reg130+reg64;
    reg47=reg55*reg2; reg111=reg111*reg133; reg73=reg79+reg73; reg87=reg95+reg87; reg79=reg100*reg38;
    reg36=reg22+reg36; reg58=reg58*reg0; reg97=reg76-reg97; reg22=reg90*reg44; reg3=reg59+reg3;
    reg59=reg71*reg114; reg76=reg118*reg90; reg80=reg63*reg114; reg24=reg24*reg39; reg95=reg7*reg90;
    reg93=reg4+reg93; reg4=reg37*reg2; reg61=reg43+reg61; reg35=reg121*reg35; reg34=reg86+reg34;
    reg32=reg32*reg39; reg89=reg89*reg0; reg47=reg64+reg47; reg42=reg65+reg42; reg16=reg117+reg16;
    reg31=reg67+reg31; reg78=reg60+reg78; reg43=reg114*reg96; reg45=reg26+reg45; reg1=reg25*reg1;
    reg69=reg88+reg69; reg19=reg82+reg19; reg14=reg62+reg14; reg111=reg73+reg111; reg9=reg68+reg9;
    reg8=reg87+reg8; reg79=reg36+reg79; reg22=reg97+reg22; reg4=reg111+reg4; reg25=reg34*reg8;
    reg43=reg78+reg43; reg19=reg19*reg39; reg32=reg16+reg32; reg45=reg45*reg3; reg41=reg41*reg0;
    reg89=reg47+reg89; reg16=reg79*reg31; reg26=reg7*reg71; reg35=reg61+reg35; reg36=reg63*reg118;
    reg58=reg14+reg58; reg14=reg71*reg81; reg83=reg83*reg1; reg47=reg15*reg90; reg76=reg59-reg76;
    reg95=reg80-reg95; reg59=reg34*reg42; reg100=reg100*reg1; reg60=reg9*reg31; reg61=reg63*reg81;
    reg62=reg54*reg90; reg24=reg93+reg24; reg69=reg69*reg3; reg69=reg24+reg69; reg33=reg33*reg1;
    reg45=reg32+reg45; reg41=reg4+reg41; reg35=reg35*reg3; reg19=reg43+reg19; reg138=reg138/reg22;
    reg76=reg76/reg22; reg4=reg9*reg42; reg62=reg14-reg62; reg143=reg143/reg22; reg14=reg8*reg79;
    reg95=reg95/reg22; reg60=reg25-reg60; reg16=reg59-reg16; reg47=reg61-reg47; reg100=reg58+reg100;
    reg83=reg89+reg83; reg24=reg63*reg54; reg26=reg36-reg26; reg25=reg15*reg71; reg35=reg19+reg35;
    reg19=reg16*reg83; reg32=reg60*reg100; reg4=reg14-reg4; reg95=reg95*reg69; reg76=reg76*reg45;
    reg138=reg138*reg45; reg143=reg143*reg69; reg25=reg24-reg25; reg26=reg26/reg22; reg44=reg44/reg22;
    reg47=reg47/reg22; reg62=reg62/reg22; reg33=reg41+reg33; reg22=reg25/reg22; reg44=reg44*reg35;
    reg143=reg138-reg143; reg14=reg8*reg33; reg32=reg19-reg32; reg19=reg4*reg33; reg26=reg26*reg35;
    reg76=reg95-reg76; reg24=reg83*reg31; reg25=reg42*reg33; reg36=reg100*reg31; reg69=reg47*reg69;
    reg45=reg62*reg45; reg69=reg45-reg69; reg41=elem.pos(2)[0]-elem.pos(0)[0]; reg35=reg22*reg35; reg22=reg79*reg33;
    reg43=elem.pos(1)[1]-elem.pos(0)[1]; reg143=reg44+reg143; reg26=reg76-reg26; reg44=elem.pos(2)[1]-elem.pos(0)[1]; reg45=elem.pos(1)[0]-elem.pos(0)[0];
    reg47=1-(*f.m).resolution; reg19=reg32+reg19; reg32=reg9*reg33; reg25=reg36-reg25; reg14=reg24-reg14;
    reg24=reg34*reg100; reg36=reg8*reg100; reg58=reg34*reg83; reg59=reg83*reg42; reg61=(*f.m).resolution*reg113;
    reg36=reg59-reg36; reg59=reg44*reg45; reg25=reg25/reg19; reg62=reg41*reg43; reg69=reg35+reg69;
    reg35=reg9*reg100; reg64=(*f.m).resolution*reg23; reg65=reg83*reg79; reg143=reg47*reg143; reg22=reg24-reg22;
    reg32=reg58-reg32; reg26=reg47*reg26; reg14=reg14/reg19; reg35=reg65-reg35; reg32=reg32/reg19;
    reg36=reg36/reg19; reg16=reg16/reg19; reg60=reg60/reg19; reg62=reg59-reg62; reg90=reg90*reg47;
    reg81=reg81*reg47; reg22=reg22/reg19; reg26=reg64+reg26; reg143=reg61+reg143; reg24=(*f.m).resolution*reg96;
    reg58=(*f.m).resolution*reg25; reg59=(*f.m).resolution*reg14; reg69=reg47*reg69; reg43=reg43/reg62; reg41=reg41/reg62;
    reg44=reg44/reg62; reg45=reg45/reg62; reg59=reg81-reg59; reg58=reg90+reg58; reg61=(*f.m).resolution*reg36;
    reg64=(*f.m).resolution*reg32; reg65=(*f.m).resolution*reg22; reg69=reg24+reg69; reg143=reg143*(*f.m).deltaT; reg26=reg26*(*f.m).deltaT;
    reg114=reg114*reg47; reg54=reg54*reg47; reg71=reg71*reg47; reg15=reg15*reg47; reg63=reg63*reg47;
    reg24=(*f.m).resolution*reg16; reg67=(*f.m).resolution*reg60; reg35=reg35/reg19; reg19=reg4/reg19; reg4=reg59*reg26;
    reg68=reg58*reg143; reg73=reg41-reg45; reg76=reg43-reg44; reg69=reg69*(*f.m).deltaT; reg61=reg114+reg61;
    reg64=reg54+reg64; reg65=reg71-reg65; reg67=reg15-reg67; reg15=(*f.m).resolution*reg19; reg7=reg7*reg47;
    reg47=reg118*reg47; reg54=(*f.m).resolution*reg35; reg63=reg24+reg63; reg24=reg63*reg143; reg71=reg67*reg26;
    reg78=reg68+reg4; reg80=reg61*reg69; reg81=0.5*reg43; reg82=0.5*reg45; reg86=0.5*reg44;
    reg87=0.5*reg41; reg88=0.5*reg73; reg89=0.5*reg76; reg90=reg65*reg143; reg93=reg64*reg26;
    reg54=reg47-reg54; reg7=reg15+reg7; reg15=reg89*reg61; reg47=reg73*reg59; reg95=reg78+reg80;
    reg97=reg90+reg93; reg98=reg54*reg69; reg99=reg7*reg69; reg105=reg24+reg71; reg106=reg76*reg58;
    reg108=reg88*reg61; reg110=reg44*reg58; reg111=reg43*reg58; reg112=reg82*reg61; reg114=reg81*reg61;
    reg116=reg41*reg59; reg117=reg86*reg61; reg118=reg87*reg61; reg120=reg45*reg59; reg122=2*reg95;
    reg15=reg47+reg15; reg47=reg41*reg67; reg123=reg86*reg7; reg124=reg81*reg54; reg126=reg97+reg98;
    reg128=reg45*reg64; reg129=reg89*reg7; reg130=reg73*reg67; reg117=reg117-reg116; reg131=reg44*reg65;
    reg132=reg87*reg54; reg112=reg112-reg111; reg134=reg43*reg63; reg135=reg41*reg64; reg136=reg86*reg54;
    reg138=reg82*reg7; reg139=reg43*reg65; reg140=reg82*reg54; reg141=reg87*reg7; reg142=reg73*reg64;
    T reg144=reg89*reg54; T reg145=reg88*reg7; T reg146=reg81*reg7; T reg147=reg105+reg99; reg120=reg120-reg114;
    reg108=reg106+reg108; reg106=reg44*reg63; T reg148=reg88*reg54; T reg149=reg45*reg67; T reg150=reg76*reg63;
    reg110=reg110-reg118; T reg151=reg76*reg65; T reg152=reg82*reg122; T reg153=reg43*reg147; reg145=reg150+reg145;
    reg136=reg136-reg135; reg138=reg138-reg134; reg150=reg45*reg126; T reg154=reg81*reg122; reg140=reg140-reg139;
    reg108=2*reg108; reg120=2*reg120; reg110=2*reg110; reg15=2*reg15; reg128=reg128-reg124;
    reg123=reg123-reg47; reg144=reg142+reg144; reg149=reg149-reg146; reg129=reg130+reg129; reg112=2*reg112;
    reg130=reg44*reg147; reg142=reg87*reg122; reg106=reg106-reg141; reg117=2*reg117; reg148=reg151+reg148;
    reg151=reg86*reg122; T reg155=reg41*reg126; reg131=reg131-reg132; T reg156=reg41*reg140; T reg157=reg86*reg112;
    T reg158=reg41*reg136; T reg159=reg88*reg108; T reg160=reg76*reg129; T reg161=reg88*reg110; T reg162=reg88*reg15;
    T reg163=reg76*reg106; T reg164=reg88*reg112; T reg165=reg76*reg138; T reg166=reg88*reg117; T reg167=reg76*reg123;
    T reg168=reg81*reg15; T reg169=reg87*reg15; T reg170=reg76*reg149; T reg171=reg88*reg120; T reg172=reg76*reg145;
    T reg173=reg73*reg148; T reg174=reg89*reg108; T reg175=reg73*reg144; T reg176=reg44*reg138; T reg177=reg87*reg112;
    T reg178=reg44*reg149; T reg179=reg87*reg120; T reg180=reg86*reg108; T reg181=reg41*reg148; T reg182=reg86*reg15;
    T reg183=reg41*reg144; T reg184=reg86*reg110; T reg185=reg41*reg131; T reg186=reg86*reg117; T reg187=reg89*reg15;
    T reg188=reg73*reg131; T reg189=reg89*reg110; T reg190=reg73*reg136; T reg191=reg89*reg117; T reg192=reg73*reg140;
    T reg193=reg44*reg129; T reg194=reg43*reg106; reg131=reg45*reg131; T reg195=reg89*reg122; T reg196=reg81*reg110;
    T reg197=reg87*reg108; T reg198=reg44*reg145; T reg199=reg82*reg112; T reg200=reg81*reg120; reg138=reg43*reg138;
    T reg201=reg89*reg120; T reg202=reg73*reg128; reg148=reg45*reg148; T reg203=reg45*reg128; T reg204=reg89*reg112;
    reg136=reg45*reg136; T reg205=reg142-reg130; T reg206=reg154-reg150; T reg207=reg81*reg117; T reg208=reg155-reg151;
    T reg209=reg153-reg152; T reg210=reg82*reg117; T reg211=reg43*reg123; reg140=reg45*reg140; reg112=reg81*reg112;
    T reg212=reg82*reg120; T reg213=reg76*reg147; reg129=reg43*reg129; T reg214=reg82*reg108; reg145=reg43*reg145;
    reg149=reg43*reg149; reg117=reg87*reg117; reg123=reg44*reg123; reg15=reg82*reg15; reg144=reg45*reg144;
    reg128=reg41*reg128; reg120=reg86*reg120; T reg215=reg88*reg122; T reg216=reg73*reg126; T reg217=reg87*reg110;
    reg106=reg44*reg106; reg110=reg82*reg110; reg108=reg81*reg108; reg211=reg210-reg211; reg117=reg123-reg117;
    reg208=reg62*reg208; reg112=reg140-reg112; reg187=reg175+reg187; reg205=reg62*reg205; reg161=reg163+reg161;
    reg177=reg176-reg177; reg123=reg213+reg215; reg149=reg212-reg149; reg185=reg184-reg185; reg179=reg178-reg179;
    reg138=reg199-reg138; reg183=reg182-reg183; reg129=reg15-reg129; reg200=reg203-reg200; reg15=reg216+reg195;
    reg181=reg180-reg181; reg174=reg173+reg174; reg164=reg165+reg164; reg156=reg157-reg156; reg158=reg186-reg158;
    reg128=reg120-reg128; reg145=reg214-reg145; reg108=reg148-reg108; reg168=reg144-reg168; reg162=reg160+reg162;
    reg217=reg106-reg217; reg169=reg193-reg169; reg194=reg110-reg194; reg172=reg159+reg172; reg197=reg198-reg197;
    reg189=reg188+reg189; reg209=reg62*reg209; reg166=reg167+reg166; reg191=reg190+reg191; reg207=reg136-reg207;
    reg206=reg62*reg206; reg204=reg192+reg204; reg171=reg170+reg171; reg201=reg202+reg201; reg196=reg131-reg196;
    reg171=reg62*reg171; reg174=reg62*reg174; reg168=reg62*reg168; reg200=reg62*reg200; reg112=reg62*reg112;
    reg164=reg62*reg164; reg149=reg62*reg149; reg207=reg62*reg207; reg196=reg62*reg196; reg108=reg62*reg108;
    reg189=reg62*reg189; reg209=ponderation*reg209; reg191=reg62*reg191; reg166=reg62*reg166; reg206=ponderation*reg206;
    reg204=reg62*reg204; reg201=reg62*reg201; reg194=reg62*reg194; reg172=reg62*reg172; reg197=reg62*reg197;
    reg169=reg62*reg169; reg162=reg62*reg162; reg217=reg62*reg217; reg129=reg129*reg62; reg128=reg62*reg128;
    reg158=reg62*reg158; reg156=reg62*reg156; reg117=reg62*reg117; reg145=reg62*reg145; reg177=reg62*reg177;
    reg161=reg62*reg161; reg106=reg62*reg123; reg179=reg62*reg179; reg181=reg62*reg181; reg138=reg62*reg138;
    reg183=reg62*reg183; reg110=reg62*reg15; reg185=reg62*reg185; reg187=reg62*reg187; reg211=reg62*reg211;
    reg205=ponderation*reg205; reg208=ponderation*reg208; T tmp_4_2=ponderation*reg194; T tmp_0_3=ponderation*reg166; T tmp_4_4=ponderation*reg138;
    T tmp_0_2=ponderation*reg161; T tmp_4_3=ponderation*reg211; T tmp_4_1=ponderation*reg129; T tmp_4_5=ponderation*reg149; T tmp_1_1=ponderation*reg187;
    sollicitation[indices[1]+0]+=-reg205; T tmp_1_2=ponderation*reg189; sollicitation[indices[1]+1]+=-reg208; sollicitation[indices[2]+0]+=-reg209; T tmp_1_3=ponderation*reg191;
    sollicitation[indices[2]+1]+=-reg206; T tmp_1_4=ponderation*reg204; T tmp_1_5=ponderation*reg201; T tmp_0_0=ponderation*reg172; T tmp_2_0=ponderation*reg197;
    T tmp_2_1=ponderation*reg169; T tmp_0_1=ponderation*reg162; T tmp_2_2=ponderation*reg217; T tmp_3_5=ponderation*reg128; T tmp_3_3=ponderation*reg158;
    T tmp_3_4=ponderation*reg156; T tmp_5_0=ponderation*reg108; T tmp_0_4=ponderation*reg164; T tmp_5_1=ponderation*reg168; T tmp_5_2=ponderation*reg196;
    T tmp_5_3=ponderation*reg207; T tmp_0_5=ponderation*reg171; T tmp_5_4=ponderation*reg112; T tmp_1_0=ponderation*reg174; T tmp_5_5=ponderation*reg200;
    T tmp_4_0=ponderation*reg145; T tmp_2_3=ponderation*reg117; T tmp_2_4=ponderation*reg177; T tmp_2_5=ponderation*reg179; reg108=ponderation*reg106;
    sollicitation[indices[0]+0]+=reg108; T tmp_3_0=ponderation*reg181; T tmp_3_1=ponderation*reg183; reg112=ponderation*reg110; sollicitation[indices[0]+1]+=reg112;
    T tmp_3_2=ponderation*reg185;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; reg1=1.0/reg1; T reg2=2*(*f.m).shear_modulus_12;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=reg2*reg3; T reg5=pow((*f.m).v1[0],2); T reg6=pow((*f.m).v1[1],2);
    T reg7=1.0/(*f.m).elastic_modulus_3; T reg8=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg9=pow((*f.m).v2[0],2); T reg10=pow((*f.m).v2[1],2); T reg11=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg12=reg4*reg8; T reg13=1.0/(*f.m).elastic_modulus_2; reg6=reg5+reg6; reg5=pow((*f.m).v1[2],2); T reg14=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1;
    T reg15=reg11*reg4; T reg16=reg4*reg7; T reg17=pow((*f.m).v2[2],2); reg10=reg9+reg10; reg9=reg11*reg15;
    reg17=reg10+reg17; reg10=reg14*reg16; T reg18=reg11*reg12; T reg19=reg13*reg16; reg5=reg6+reg5;
    reg9=reg19-reg9; reg6=1.0/(*f.m).elastic_modulus_1; reg18=reg10+reg18; reg19=reg14*reg15; T reg20=reg13*reg12;
    reg17=pow(reg17,0.5); reg5=pow(reg5,0.5); T reg21=reg6*reg9; T reg22=(*f.m).v1[2]/reg5; T reg23=(*f.m).v2[1]/reg17;
    T reg24=reg14*reg18; T reg25=(*f.m).v1[1]/reg5; T reg26=(*f.m).v2[2]/reg17; T reg27=reg19+reg20; T reg28=reg14*reg4;
    T reg29=reg12*reg8; reg16=reg6*reg16; T reg30=reg2*reg0; T reg31=reg15*reg8; reg4=reg13*reg4;
    T reg32=reg11*reg3; T reg33=reg3*reg8; T reg34=reg22*reg23; T reg35=reg27*reg8; T reg36=reg25*reg26;
    reg24=reg21-reg24; reg3=reg3*reg7; reg5=(*f.m).v1[0]/reg5; reg17=(*f.m).v2[0]/reg17; reg21=reg4*reg8;
    reg35=reg24-reg35; reg31=reg10+reg31; reg10=reg36-reg34; reg24=reg22*reg17; T reg37=2*reg5;
    T reg38=reg5*reg26; T reg39=2*reg17; T reg40=reg30*reg7; T reg41=reg11*reg33; T reg42=reg11*reg32;
    T reg43=reg14*reg3; reg3=reg13*reg3; T reg44=reg30*reg8; T reg45=reg2*reg1; reg29=reg16-reg29;
    reg15=reg6*reg15; reg16=reg28*reg8; reg12=reg14*reg12; reg30=reg11*reg30; T reg46=reg24-reg38;
    reg29=reg29/reg35; T reg47=reg5*reg23; T reg48=reg25*reg17; reg18=reg18/reg35; T reg49=pow(reg25,2);
    T reg50=pow(reg5,2); T reg51=reg45*reg7; T reg52=reg11*reg44; T reg53=reg11*reg30; T reg54=reg14*reg40;
    reg40=reg13*reg40; reg21=reg19+reg21; reg41=reg43+reg41; reg16=reg15+reg16; reg42=reg3-reg42;
    reg31=reg31/reg35; reg9=reg9/reg35; reg3=reg45*reg8; reg12=reg15+reg12; reg33=reg13*reg33;
    reg45=reg11*reg45; reg28=reg14*reg28; reg32=reg14*reg32; reg4=reg6*reg4; reg15=reg25*reg37;
    reg43=pow(reg17,2); T reg55=pow(reg23,2); T reg56=reg39*reg23; T reg57=2*reg10; T reg58=reg11*reg3;
    T reg59=reg11*reg45; T reg60=reg14*reg51; reg51=reg13*reg51; T reg61=reg57*reg46; reg52=reg54+reg52;
    reg53=reg40-reg53; reg40=reg32+reg33; reg41=reg14*reg41; reg42=reg6*reg42; reg44=reg13*reg44;
    reg30=reg14*reg30; reg28=reg4-reg28; reg12=reg12/reg35; reg16=reg16/reg35; reg4=reg9*reg49;
    reg54=reg55*reg18; T reg62=reg15*reg9; T reg63=reg56*reg18; T reg64=reg31*reg50; T reg65=reg43*reg29;
    T reg66=pow(reg46,2); T reg67=pow(reg22,2); T reg68=reg56*reg29; T reg69=reg15*reg31; T reg70=pow(reg10,2);
    T reg71=reg55*reg29; T reg72=reg31*reg49; T reg73=pow(reg26,2); T reg74=reg47-reg48; T reg75=reg43*reg18;
    reg21=reg21/reg35; T reg76=reg9*reg50; reg27=reg27/reg35; reg3=reg13*reg3; reg68=reg69+reg68;
    reg69=reg15*reg21; reg41=reg42-reg41; reg42=reg61*reg12; reg40=reg40*reg8; T reg77=reg61*reg27;
    reg63=reg62+reg63; reg62=reg56*reg16; reg53=reg6*reg53; T reg78=reg73*reg18; T reg79=reg9*reg67;
    reg52=reg14*reg52; T reg80=reg66*reg27; reg54=reg4+reg54; reg4=reg30+reg44; T reg81=reg21*reg50;
    reg58=reg60+reg58; reg59=reg51-reg59; reg51=reg43*reg16; reg60=reg55*reg16; T reg82=reg21*reg49;
    reg75=reg76+reg75; reg71=reg72+reg71; reg72=reg66*reg12; reg28=reg28/reg35; reg76=pow(reg74,2);
    T reg83=reg70*reg27; reg65=reg64+reg65; reg64=reg31*reg67; T reg84=reg73*reg29; reg45=reg14*reg45;
    T reg85=reg70*reg12; T reg86=reg45+reg3; reg83=reg75+reg83; reg75=reg37*reg22; reg58=reg14*reg58;
    T reg87=2*reg25; reg62=reg69+reg62; reg80=reg54+reg80; reg85=reg65+reg85; reg60=reg82+reg60;
    reg54=reg66*reg28; reg78=reg79+reg78; reg65=reg27*reg76; reg69=reg5*reg25; reg79=reg39*reg26;
    reg77=reg63+reg77; reg42=reg68+reg42; reg63=2*reg23; reg72=reg71+reg72; reg68=reg17*reg23;
    reg71=reg73*reg16; reg82=reg21*reg67; T reg88=reg12*reg76; reg84=reg64+reg84; reg40=reg41-reg40;
    reg51=reg81+reg51; reg52=reg53-reg52; reg41=reg70*reg28; reg4=reg4*reg8; reg53=reg61*reg28;
    reg59=reg6*reg59; reg71=reg82+reg71; reg64=reg69*reg77; reg81=reg63*reg26; reg88=reg84+reg88;
    reg82=reg28*reg76; reg84=reg46*reg10; reg41=reg51+reg41; reg51=reg68*reg42; T reg89=reg25*reg10;
    T reg90=reg83*reg50; reg40=reg40/reg35; T reg91=reg5*reg46; reg54=reg60+reg54; reg60=reg43*reg42;
    T reg92=reg77*reg50; T reg93=reg43*reg85; T reg94=reg87*reg22; T reg95=reg43*reg72; T reg96=reg83*reg49;
    T reg97=reg55*reg85; T reg98=reg80*reg49; T reg99=reg55*reg72; reg77=reg77*reg49; reg42=reg55*reg42;
    T reg100=reg80*reg50; reg58=reg59-reg58; reg59=reg87*reg23; reg86=reg86*reg8; reg53=reg62+reg53;
    reg62=2*reg22; reg4=reg52-reg4; reg52=reg37*reg17; reg57=reg57*reg74; T reg101=2*reg46;
    reg80=reg69*reg80; reg65=reg78+reg65; reg78=reg5*reg17; T reg102=reg25*reg23; reg72=reg68*reg72;
    T reg103=reg75*reg31; T reg104=reg75*reg9; reg47=reg48+reg47; reg48=reg79*reg18; T reg105=reg79*reg29;
    T reg106=reg66*reg41; reg85=reg68*reg85; T reg107=reg75*reg21; reg82=reg71+reg82; reg99=reg98+reg99;
    reg71=reg66*reg54; reg83=reg69*reg83; reg98=reg65*reg49; T reg108=reg55*reg88; reg72=reg80+reg72;
    reg80=reg78*reg40; reg42=reg77+reg42; reg77=reg66*reg53; T reg109=reg70*reg41; reg93=reg90+reg93;
    reg90=reg47*reg40; reg51=reg64+reg51; reg64=reg79*reg16; T reg110=reg70*reg53; reg60=reg92+reg60;
    reg29=reg81*reg29; reg31=reg94*reg31; reg92=reg57*reg12; reg53=reg84*reg53; reg105=reg103+reg105;
    reg103=reg43*reg88; T reg111=reg84*reg54; T reg112=reg65*reg50; reg54=reg70*reg54; reg95=reg100+reg95;
    reg100=reg102*reg40; reg97=reg96+reg97; reg96=reg59*reg13; reg62=reg62*reg26; T reg113=reg52*reg14;
    T reg114=reg55*reg13; T reg115=reg43*reg14; T reg116=reg25*reg46; T reg117=reg5*reg10; reg4=reg4/reg35;
    reg101=reg101*reg74; reg86=reg58-reg86; reg91=reg89+reg91; reg58=reg22*reg26; reg48=reg104+reg48;
    reg89=reg57*reg27; reg9=reg94*reg9; reg18=reg81*reg18; reg104=reg10*reg23; T reg118=reg43*reg6;
    T reg119=reg55*reg14; T reg120=reg59*reg14; T reg121=reg52*reg6; T reg122=reg46*reg17; T reg123=reg117*reg4;
    T reg124=reg46*reg23; reg106=reg97+reg106; reg97=reg58*reg40; T reg125=reg59*reg80; reg113=reg96-reg113;
    reg96=reg55*reg11; T reg126=reg52*reg90; reg41=reg84*reg41; reg16=reg81*reg16; reg21=reg94*reg21;
    T reg127=reg52*reg80; T reg128=reg116*reg4; T reg129=reg74*reg22; T reg130=reg91*reg4; T reg131=reg59*reg11;
    reg109=reg93+reg109; reg122=reg104+reg122; reg88=reg68*reg88; reg93=reg52*reg100; reg104=reg43*reg8;
    reg65=reg69*reg65; reg111=reg72+reg111; reg72=reg47*reg100; T reg132=reg10*reg17; reg54=reg95+reg54;
    reg95=reg66*reg82; reg120=reg121-reg120; reg35=reg86/reg35; reg12=reg101*reg12; reg29=reg31+reg29;
    reg92=reg105+reg92; reg89=reg48+reg89; reg18=reg9+reg18; reg27=reg101*reg27; reg9=reg62*reg8;
    reg119=reg118-reg119; reg31=reg73*reg8; reg48=reg47*reg90; reg77=reg42+reg77; reg90=reg59*reg90;
    reg53=reg51+reg53; reg42=reg52*reg8; reg85=reg83+reg85; reg51=reg57*reg28; reg64=reg107+reg64;
    reg37=reg37*reg10; reg87=reg87*reg46; reg103=reg112+reg103; reg83=reg62*reg11; reg71=reg99+reg71;
    reg115=reg114-reg115; reg86=reg73*reg11; reg100=reg59*reg100; reg99=reg70*reg82; reg110=reg60+reg110;
    reg108=reg98+reg108; reg60=reg10*reg22; reg98=reg5*reg74; reg126=reg110+reg126; reg105=reg37*reg130;
    reg107=reg89*reg50; reg110=reg87*reg130; reg90=reg77+reg90; reg77=reg89*reg49; reg112=reg55*reg92;
    reg114=reg59*reg97; reg95=reg108+reg95; reg108=reg87*reg128; reg100=reg71+reg100; reg71=(*f.m).alpha_1*reg50;
    reg118=reg87*reg123; reg125=reg106+reg125; reg131=reg42+reg131; reg62=reg62*reg7; reg41=reg85+reg41;
    reg80=reg47*reg80; reg24=reg38+reg24; reg72=reg111+reg72; reg38=reg91*reg128; reg88=reg65+reg88;
    reg82=reg84*reg82; reg42=reg43*reg92; reg31=reg119-reg31; reg6=reg6*reg50; reg65=reg43*(*f.m).alpha_2;
    reg85=reg122*reg35; reg27=reg18+reg27; reg18=reg124*reg35; reg12=reg29+reg12; reg9=reg120-reg9;
    reg29=reg132*reg35; reg106=reg14*reg50; reg111=(*f.m).alpha_1*reg49; reg119=reg129*reg4; reg120=reg55*(*f.m).alpha_2;
    reg83=reg113-reg83; reg28=reg101*reg28; reg16=reg21+reg16; reg51=reg64+reg51; reg39=reg39*reg10;
    reg63=reg63*reg46; reg86=reg115-reg86; reg13=reg13*reg49; reg96=reg104+reg96; reg48=reg53+reg48;
    reg130=reg91*reg130; reg21=reg52*reg97; reg99=reg103+reg99; reg53=reg73*reg7; reg14=reg14*reg49;
    reg128=reg37*reg128; reg93=reg54+reg93; reg54=reg74*reg26; reg64=reg37*reg123; reg127=reg109+reg127;
    reg14=reg6-reg14; reg6=reg63*reg85; reg103=reg11*reg67; reg104=reg49*reg83; reg109=reg55*reg12;
    reg112=reg77+reg112; reg77=reg43*reg9; reg113=reg66*reg51; reg115=reg55*reg83; reg121=(*f.m).alpha_1*reg67;
    T reg133=reg27*reg49; T reg134=reg73*(*f.m).alpha_2; T reg135=reg43*reg31; T reg136=reg55*reg86; T reg137=reg78*reg31;
    T reg138=reg102*reg86; T reg139=reg67*reg8; reg130=reg48+reg130; reg48=reg122*reg85; reg89=reg69*reg89;
    reg92=reg68*reg92; reg97=reg47*reg97; reg82=reg88+reg82; reg96=reg53-reg96; reg11=reg11*reg49;
    reg53=reg122*reg18; reg38=reg72+reg38; reg72=reg78*reg9; reg88=reg102*reg83; reg131=reg62-reg131;
    reg71=reg65+reg71; reg123=reg91*reg123; reg80=reg41+reg80; reg41=(*f.m).alpha_3*reg70; reg120=reg111+reg120;
    reg8=reg50*reg8; reg62=(*f.m).alpha_3*reg66; reg65=reg50*reg31; reg111=reg49*reg86; T reg140=reg50*reg9;
    reg106=reg13-reg106; reg28=reg16+reg28; reg13=reg24*reg40; reg16=reg54*reg35; reg108=reg100+reg108;
    reg100=reg63*reg29; reg118=reg125+reg118; reg64=reg127+reg64; reg125=reg39*reg29; reg128=reg93+reg128;
    reg93=reg39*reg18; reg127=reg43*reg12; T reg141=reg27*reg50; reg21=reg99+reg21; reg99=reg37*reg119;
    reg105=reg126+reg105; reg85=reg39*reg85; reg126=reg70*reg51; reg42=reg107+reg42; reg110=reg90+reg110;
    reg90=reg87*reg119; reg114=reg95+reg114; reg36=reg34+reg36; reg98=reg60+reg98; reg34=reg46*reg22;
    reg60=reg25*reg74; reg95=reg10*reg26; reg18=reg63*reg18; reg107=reg74*reg17; reg27=reg69*reg27;
    reg51=reg84*reg51; reg92=reg89+reg92; reg125=reg64+reg125; reg64=reg74*reg23; reg48=reg130+reg48;
    reg89=reg68*reg2; reg93=reg128+reg93; reg128=reg47*reg2; reg130=reg46*reg26; reg107=reg95+reg107;
    reg111=reg65+reg111; reg99=reg21+reg99; reg21=reg39*reg16; reg65=reg67*reg96; reg95=reg67*reg131;
    reg85=reg105+reg85; reg119=reg91*reg119; reg104=reg140+reg104; reg7=reg67*reg7; reg62=reg120+reg62;
    reg11=reg8+reg11; reg134=reg121+reg134; reg40=reg36*reg40; reg8=reg98*reg4; reg41=reg71+reg41;
    reg71=(*f.m).alpha_3*reg76; reg103=reg106-reg103; reg136=reg135+reg136; reg105=reg73*reg96; reg115=reg77+reg115;
    reg77=reg73*reg131; reg106=(*f.m).alpha_1*reg69; reg120=reg58*reg131; reg60=reg34+reg60; reg34=reg68*(*f.m).alpha_2;
    reg88=reg72+reg88; reg17=reg17*reg26; reg5=reg5*reg22; reg12=reg68*reg12; reg100=reg118+reg100;
    reg113=reg112+reg113; reg68=reg59*reg13; reg18=reg108+reg18; reg109=reg133+reg109; reg90=reg114+reg90;
    reg72=reg66*reg28; reg108=reg63*reg16; reg123=reg80+reg123; reg138=reg137+reg138; reg29=reg122*reg29;
    reg80=reg58*reg96; reg112=reg70*reg28; reg139=reg14-reg139; reg14=reg52*reg13; reg53=reg38+reg53;
    reg126=reg42+reg126; reg6=reg110+reg6; reg97=reg82+reg97; reg127=reg141+reg127; reg38=reg18*reg48;
    reg42=reg93*reg48; reg82=reg24*reg1; reg110=reg17*reg1; reg114=reg6*reg53; reg26=reg23*reg26;
    reg68=reg113+reg68; reg23=reg56*reg128; reg77=reg115+reg77; reg113=reg87*reg8; reg22=reg25*reg22;
    reg71=reg134+reg71; reg25=reg85*reg53; reg72=reg109+reg72; reg59=reg59*reg40; reg11=reg7-reg11;
    reg34=reg106+reg34; reg7=reg47*reg89; reg106=reg56*reg89; reg105=reg136+reg105; reg64=reg130+reg64;
    reg10=reg74*reg10; reg17=reg17*(*f.m).alpha_2; reg109=(*f.m).alpha_1*reg5; reg115=(*f.m).alpha_3*reg84; reg118=reg55*reg103;
    reg121=reg107*reg35; reg80=reg138+reg80; reg4=reg60*reg4; reg130=reg18*reg62; reg28=reg84*reg28;
    reg12=reg27+reg12; reg65=reg111+reg65; reg27=reg15*reg89; reg13=reg47*reg13; reg51=reg92+reg51;
    reg29=reg123+reg29; reg52=reg52*reg40; reg112=reg127+reg112; reg84=reg15*reg128; reg95=reg104+reg95;
    reg21=reg99+reg21; reg92=reg37*reg8; reg14=reg126+reg14; reg108=reg90+reg108; reg16=reg122*reg16;
    reg119=reg97+reg119; reg90=reg47*reg128; reg120=reg88+reg120; reg88=reg50*reg139; reg97=reg49*reg103;
    reg99=reg43*reg139; reg104=reg125*reg41; reg111=reg93*reg62; reg123=reg100*reg41; reg126=reg85*reg18;
    reg127=reg93*reg6; reg25=reg42-reg25; reg21=reg21*reg71; reg114=reg38-reg114; reg113=reg68+reg113;
    reg38=reg63*reg121; reg42=reg36*reg0; reg68=reg26*reg0; reg84=reg95+reg84; reg130=reg123+reg130;
    reg115=reg34+reg115; reg97=reg88+reg97; reg108=reg108*reg71; reg17=reg109+reg17; reg10=(*f.m).alpha_3*reg10;
    reg34=(*f.m).alpha_1*reg22; reg26=reg26*(*f.m).alpha_2; reg88=reg75*reg82; reg95=reg75*reg110; reg27=reg65+reg27;
    reg87=reg87*reg4; reg59=reg72+reg59; reg67=reg67*reg11; reg65=reg29*reg41; reg2=reg69*reg2;
    reg111=reg104+reg111; reg69=reg53*reg62; reg72=reg79*reg82; reg23=reg77+reg23; reg77=reg78*reg139;
    reg104=reg102*reg103; reg40=reg47*reg40; reg109=reg79*reg110; reg106=reg105+reg106; reg28=reg12+reg28;
    reg46=reg74*reg46; reg12=reg24*reg82; reg13=reg51+reg13; reg37=reg37*reg4; reg90=reg120+reg90;
    reg7=reg80+reg7; reg92=reg14+reg92; reg35=reg64*reg35; reg118=reg99+reg118; reg52=reg112+reg52;
    reg14=reg39*reg121; reg8=reg91*reg8; reg16=reg119+reg16; reg73=reg73*reg11; reg51=reg24*reg110;
    reg67=reg97+reg67; reg58=reg58*reg11; reg104=reg77+reg104; reg74=reg15*reg2; reg14=reg92+reg14;
    reg1=reg5*reg1; reg5=reg81*reg42; reg77=reg125*reg114; reg80=reg100*reg25; reg92=reg36*reg42;
    reg97=reg6*reg115; reg12=reg90+reg12; reg39=reg39*reg35; reg37=reg52+reg37; reg108=reg130+reg108;
    reg4=reg91*reg4; reg126=reg127-reg126; reg52=reg85*reg115; reg21=reg111+reg21; reg90=reg36*reg68;
    reg99=reg94*reg68; reg51=reg7+reg51; reg95=reg27+reg95; reg73=reg118+reg73; reg7=reg56*reg2;
    reg8=reg13+reg8; reg10=reg17+reg10; reg88=reg84+reg88; reg121=reg122*reg121; reg26=reg34+reg26;
    reg46=(*f.m).alpha_3*reg46; reg63=reg63*reg35; reg72=reg23+reg72; reg13=reg94*reg42; reg38=reg113+reg38;
    reg40=reg28+reg40; reg17=reg81*reg68; reg109=reg106+reg109; reg69=reg65+reg69; reg16=reg16*reg71;
    reg87=reg59+reg87; reg38=reg38*reg10; reg4=reg40+reg4; reg35=reg122*reg35; reg14=reg14*reg10;
    reg52=reg21+reg52; reg97=reg108+reg97; reg121=reg8+reg121; reg8=reg48*reg115; reg16=reg69+reg16;
    reg75=reg75*reg1; reg74=reg67+reg74; reg13=reg88+reg13; reg99=reg95+reg99; reg90=reg51+reg90;
    reg63=reg87+reg63; reg7=reg73+reg7; reg79=reg79*reg1; reg46=reg26+reg46; reg17=reg109+reg17;
    reg5=reg72+reg5; reg0=reg22*reg0; reg58=reg104+reg58; reg21=reg47*reg2; reg92=reg12+reg92;
    reg80=reg77-reg80; reg12=reg29*reg126; reg22=reg100*reg48; reg23=reg6*reg29; reg39=reg37+reg39;
    reg26=reg85*reg29; reg27=reg125*reg48; reg38=reg97+reg38; reg63=reg63*reg46; reg28=reg125*reg6;
    reg34=reg5*reg90; reg75=reg74+reg75; reg37=reg92*reg17; reg94=reg94*reg0; reg40=reg92*reg99;
    reg51=reg13*reg90; reg8=reg16+reg8; reg121=reg121*reg10; reg16=reg93*reg29; reg81=reg81*reg0;
    reg79=reg7+reg79; reg7=reg125*reg53; reg26=reg27-reg26; reg21=reg58+reg21; reg24=reg24*reg1;
    reg23=reg22-reg23; reg22=reg100*reg53; reg14=reg52+reg14; reg39=reg39*reg46; reg35=reg4+reg35;
    reg4=reg18*reg29; reg27=reg85*reg100; reg12=reg80+reg12; reg26=reg26/reg12; reg52=reg93*reg100;
    reg4=reg22-reg4; reg24=reg21+reg24; reg34=reg37-reg34; reg81=reg79+reg81; reg51=reg40-reg51;
    reg35=reg35*reg46; reg121=reg8+reg121; reg8=reg125*reg18; reg23=reg23/reg12; reg21=reg99*reg5;
    reg22=reg13*reg17; reg63=reg38+reg63; reg25=reg25/reg12; reg27=reg28-reg27; reg94=reg75+reg94;
    reg39=reg14+reg39; reg114=reg114/reg12; reg16=reg7-reg16; reg36=reg36*reg0; reg25=reg25*reg63;
    reg114=reg114*reg39; reg26=reg26*reg63; reg23=reg23*reg39; reg35=reg121+reg35; reg52=reg8-reg52;
    reg27=reg27/reg12; reg126=reg126/reg12; reg16=reg16/reg12; reg4=reg4/reg12; reg36=reg24+reg36;
    reg7=reg51*reg81; reg22=reg21-reg22; reg8=reg34*reg94; reg63=reg16*reg63; reg39=reg4*reg39;
    reg23=reg26-reg23; reg27=reg27*reg35; reg4=reg22*reg36; reg7=reg8-reg7; reg8=reg81*reg90;
    reg14=reg17*reg36; reg16=reg94*reg90; reg12=reg52/reg12; reg21=reg99*reg36; reg126=reg126*reg35;
    reg25=reg114-reg25; reg24=reg13*reg36; reg26=1-(*f.m).resolution; reg14=reg8-reg14; reg21=reg16-reg21;
    reg8=reg99*reg81; reg16=reg92*reg94; reg28=reg94*reg17; reg25=reg126+reg25; reg37=elem.pos(2)[1]-elem.pos(0)[1];
    reg27=reg23-reg27; reg35=reg12*reg35; reg63=reg39-reg63; reg12=elem.pos(2)[0]-elem.pos(0)[0]; reg23=elem.pos(1)[0]-elem.pos(0)[0];
    reg38=reg5*reg36; reg39=elem.pos(1)[1]-elem.pos(0)[1]; reg4=reg7+reg4; reg7=reg92*reg81; reg14=reg14/reg4;
    reg38=reg7-reg38; reg24=reg16-reg24; reg21=reg21/reg4; reg7=reg37*reg23; reg16=reg94*reg5;
    reg40=reg13*reg81; reg52=reg12*reg39; reg25=reg26*reg25; reg8=reg28-reg8; reg28=(*f.m).resolution*reg62;
    reg27=reg26*reg27; reg63=reg35+reg63; reg35=(*f.m).resolution*reg41; reg58=(*f.m).resolution*reg21; reg59=(*f.m).resolution*reg14;
    reg65=(*f.m).resolution*reg115; reg38=reg38/reg4; reg53=reg53*reg26; reg29=reg29*reg26; reg27=reg28+reg27;
    reg51=reg51/reg4; reg25=reg35+reg25; reg24=reg24/reg4; reg34=reg34/reg4; reg52=reg7-reg52;
    reg63=reg26*reg63; reg8=reg8/reg4; reg40=reg16-reg40; reg63=reg65+reg63; reg59=reg29+reg59;
    reg25=reg25*(*f.m).deltaT; reg27=reg27*(*f.m).deltaT; reg58=reg53-reg58; reg7=(*f.m).resolution*reg8; reg37=reg37/reg52;
    reg23=reg23/reg52; reg16=(*f.m).resolution*reg24; reg28=(*f.m).resolution*reg38; reg48=reg48*reg26; reg18=reg18*reg26;
    reg100=reg100*reg26; reg93=reg93*reg26; reg125=reg125*reg26; reg29=(*f.m).resolution*reg34; reg35=(*f.m).resolution*reg51;
    reg22=reg22/reg4; reg39=reg39/reg52; reg12=reg12/reg52; reg4=reg40/reg4; reg40=reg39-reg37;
    reg53=reg12-reg23; reg65=reg58*reg27; reg67=reg59*reg25; reg7=reg48+reg7; reg16=reg18+reg16;
    reg28=reg100-reg28; reg35=reg93-reg35; reg125=reg29+reg125; reg63=reg63*(*f.m).deltaT; reg18=(*f.m).resolution*reg4;
    reg6=reg6*reg26; reg26=reg85*reg26; reg29=(*f.m).resolution*reg22; reg26=reg29+reg26; reg29=reg28*reg25;
    reg18=reg6-reg18; reg6=0.5*reg39; reg48=0.5*reg23; reg69=0.5*reg37; reg72=0.5*reg12;
    reg73=0.5*reg40; reg74=reg7*reg63; reg75=reg16*reg27; reg77=reg67+reg65; reg79=0.5*reg53;
    reg80=reg125*reg25; reg84=reg35*reg27; reg85=reg29+reg75; reg87=reg18*reg63; reg88=reg53*reg58;
    reg93=reg48*reg7; reg95=reg73*reg7; reg97=reg39*reg59; reg100=reg80+reg84; reg104=reg72*reg7;
    reg105=reg26*reg63; reg106=reg69*reg7; reg108=reg40*reg59; reg109=reg79*reg7; reg111=reg12*reg58;
    reg112=reg37*reg59; reg113=reg77+reg74; reg114=reg6*reg7; reg118=reg23*reg58; reg119=2*reg113;
    reg120=reg53*reg35; reg121=reg73*reg26; reg123=reg100+reg105; reg126=reg85+reg87; reg93=reg93-reg97;
    reg112=reg112-reg104; reg127=reg48*reg26; reg130=reg39*reg125; reg106=reg106-reg111; reg133=reg69*reg26;
    reg134=reg12*reg35; reg135=reg23*reg35; reg136=reg6*reg26; reg118=reg118-reg114; reg137=reg72*reg26;
    reg138=reg53*reg16; reg140=reg73*reg18; reg141=reg37*reg125; reg109=reg108+reg109; reg95=reg88+reg95;
    reg88=reg79*reg26; reg108=reg40*reg125; T reg142=reg37*reg28; T reg143=reg72*reg18; T reg144=reg23*reg16;
    T reg145=reg6*reg18; T reg146=reg12*reg16; T reg147=reg69*reg18; T reg148=reg39*reg28; T reg149=reg48*reg18;
    reg95=2*reg95; reg141=reg141-reg137; reg112=2*reg112; reg121=reg120+reg121; reg149=reg149-reg148;
    reg147=reg147-reg146; reg142=reg142-reg143; reg88=reg108+reg88; reg109=2*reg109; reg140=reg138+reg140;
    reg118=2*reg118; reg135=reg135-reg136; reg133=reg133-reg134; reg106=2*reg106; reg127=reg127-reg130;
    reg93=2*reg93; reg108=reg69*reg119; reg120=reg6*reg119; reg138=reg12*reg126; T reg150=reg72*reg119;
    T reg151=reg37*reg123; T reg152=reg23*reg126; reg144=reg144-reg145; T reg153=reg48*reg119; T reg154=reg39*reg123;
    T reg155=reg79*reg119; T reg156=reg79*reg118; T reg157=reg40*reg123; T reg158=reg72*reg106; T reg159=reg79*reg109;
    T reg160=reg37*reg133; T reg161=reg53*reg140; T reg162=reg37*reg127; T reg163=reg154-reg153; T reg164=reg72*reg112;
    T reg165=reg37*reg141; T reg166=reg53*reg149; T reg167=reg48*reg93; T reg168=reg72*reg93; T reg169=reg120-reg152;
    T reg170=reg73*reg106; T reg171=reg53*reg147; T reg172=reg40*reg88; T reg173=reg73*reg118; T reg174=reg53*reg144;
    T reg175=reg73*reg93; T reg176=reg37*reg135; T reg177=reg73*reg112; T reg178=reg53*reg142; T reg179=reg72*reg118;
    T reg180=reg69*reg106; T reg181=reg73*reg95; T reg182=reg40*reg127; T reg183=reg12*reg147; T reg184=reg79*reg112;
    T reg185=reg69*reg93; reg127=reg39*reg127; T reg186=reg79*reg106; T reg187=reg12*reg149; T reg188=reg40*reg141;
    T reg189=reg150-reg151; T reg190=reg73*reg119; T reg191=reg48*reg118; T reg192=reg40*reg133; T reg193=reg79*reg95;
    T reg194=reg53*reg126; T reg195=reg39*reg135; T reg196=reg138-reg108; T reg197=reg79*reg93; T reg198=reg69*reg118;
    T reg199=reg23*reg144; T reg200=reg40*reg121; reg118=reg6*reg118; reg144=reg12*reg144; reg135=reg40*reg135;
    reg179=reg176-reg179; reg175=reg166+reg175; reg197=reg182+reg197; reg166=reg194+reg190; reg186=reg192+reg186;
    reg169=reg52*reg169; reg168=reg162-reg168; reg189=reg52*reg189; reg158=reg160-reg158; reg184=reg188+reg184;
    reg156=reg135+reg156; reg163=reg52*reg163; reg196=reg52*reg196; reg144=reg198-reg144; reg187=reg185-reg187;
    reg183=reg180-reg183; reg127=reg167-reg127; reg193=reg200+reg193; reg195=reg191-reg195; reg118=reg199-reg118;
    reg172=reg159+reg172; reg135=reg157+reg155; reg164=reg165-reg164; reg170=reg171+reg170; reg173=reg174+reg173;
    reg181=reg161+reg181; reg177=reg178+reg177; reg156=reg52*reg156; reg159=reg52*reg135; reg144=reg52*reg144;
    reg172=reg52*reg172; reg186=reg52*reg186; reg118=reg52*reg118; reg196=ponderation*reg196; reg127=reg52*reg127;
    reg197=reg52*reg197; reg195=reg52*reg195; reg184=reg52*reg184; reg189=ponderation*reg189; reg175=reg52*reg175;
    reg179=reg52*reg179; reg193=reg52*reg193; reg163=ponderation*reg163; reg169=ponderation*reg169; reg183=reg52*reg183;
    reg160=reg52*reg166; reg164=reg52*reg164; reg158=reg52*reg158; reg177=reg52*reg177; reg170=reg52*reg170;
    reg187=reg52*reg187; reg168=reg52*reg168; reg173=reg52*reg173; reg181=reg52*reg181; sollicitation[indices[1]+0]+=-reg189;
    T tmp_3_4=ponderation*reg187; T tmp_4_4=ponderation*reg127; T tmp_2_5=ponderation*reg179; T tmp_3_3=ponderation*reg183; T tmp_3_5=ponderation*reg144;
    T tmp_0_3=ponderation*reg186; reg127=ponderation*reg160; sollicitation[indices[0]+1]+=reg127; T tmp_1_1=ponderation*reg181; sollicitation[indices[2]+1]+=-reg169;
    T tmp_1_2=ponderation*reg177; T tmp_1_3=ponderation*reg170; T tmp_1_5=ponderation*reg173; T tmp_2_3=ponderation*reg158; T tmp_2_2=ponderation*reg164;
    reg144=ponderation*reg159; sollicitation[indices[0]+0]+=reg144; sollicitation[indices[2]+0]+=-reg163; T tmp_0_5=ponderation*reg156; T tmp_2_4=ponderation*reg168;
    T tmp_5_5=ponderation*reg118; T tmp_0_0=ponderation*reg172; sollicitation[indices[1]+1]+=-reg196; T tmp_1_4=ponderation*reg175; T tmp_4_5=ponderation*reg195;
    T tmp_0_4=ponderation*reg197; T tmp_0_2=ponderation*reg184; T tmp_0_1=ponderation*reg193;
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
    T reg0=pow((*f.m).v1[0],2); T reg1=pow((*f.m).v1[1],2); T reg2=pow((*f.m).v2[1],2); T reg3=pow((*f.m).v2[0],2); T reg4=pow((*f.m).v1[2],2);
    reg1=reg0+reg1; reg0=2*(*f.m).shear_modulus_23; reg4=reg1+reg4; reg2=reg3+reg2; reg1=pow((*f.m).v2[2],2);
    reg3=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; reg4=pow(reg4,0.5); reg3=1.0/reg3; T reg5=2*(*f.m).shear_modulus_12;
    reg1=reg2+reg1; reg1=pow(reg1,0.5); reg5=1.0/reg5; reg2=reg3*reg0; T reg6=(*f.m).v1[0]/reg4;
    T reg7=(*f.m).v1[1]/reg4; T reg8=reg5*reg2; T reg9=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; reg4=(*f.m).v1[2]/reg4; T reg10=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg11=(*f.m).v2[1]/reg1; T reg12=(*f.m).v2[0]/reg1; T reg13=2*reg7; T reg14=1.0/(*f.m).elastic_modulus_3; T reg15=2*reg6;
    reg1=(*f.m).v2[2]/reg1; T reg16=2*reg4; T reg17=reg13*reg11; T reg18=reg15*reg12; T reg19=1.0/(*f.m).elastic_modulus_1;
    T reg20=1.0/(*f.m).elastic_modulus_2; T reg21=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg22=pow(reg11,2); T reg23=pow(reg12,2); T reg24=reg8*reg9;
    T reg25=reg10*reg8; T reg26=reg8*reg14; reg16=reg16*reg1; T reg27=reg10*reg24; T reg28=reg21*reg26;
    T reg29=reg10*reg25; T reg30=reg20*reg26; T reg31=reg22*reg21; T reg32=reg23*reg19; T reg33=reg17*reg20;
    T reg34=reg18*reg21; T reg35=reg18*reg19; T reg36=reg17*reg21; T reg37=reg23*reg21; T reg38=reg22*reg20;
    T reg39=pow(reg1,2); reg34=reg33-reg34; reg33=reg16*reg10; reg37=reg38-reg37; reg38=reg39*reg10;
    reg29=reg30-reg29; reg30=pow(reg7,2); reg31=reg32-reg31; reg32=reg39*reg9; T reg40=pow(reg6,2);
    T reg41=reg16*reg9; T reg42=reg23*reg9; reg36=reg35-reg36; reg35=reg22*reg10; T reg43=reg18*reg9;
    T reg44=reg17*reg10; T reg45=reg20*reg24; T reg46=reg21*reg25; reg27=reg28+reg27; reg44=reg43+reg44;
    reg16=reg16*reg14; reg43=reg21*reg27; reg38=reg37-reg38; reg37=reg7*reg12; T reg47=reg6*reg11;
    reg33=reg34-reg33; reg35=reg42+reg35; reg34=pow(reg4,2); reg42=reg39*reg14; T reg48=reg19*reg29;
    T reg49=reg19*reg40; T reg50=reg46+reg45; T reg51=reg21*reg30; T reg52=reg6*reg12; reg32=reg31-reg32;
    reg31=reg7*reg11; reg41=reg36-reg41; reg36=reg21*reg40; T reg53=reg20*reg30; T reg54=reg22*reg33;
    T reg55=reg40*reg32; T reg56=reg4*reg1; T reg57=reg37+reg47; T reg58=reg5*reg0; T reg59=reg6*reg1;
    T reg60=reg4*reg12; T reg61=reg4*reg11; T reg62=reg7*reg1; T reg63=reg10*reg2; T reg64=reg2*reg9;
    reg44=reg16-reg44; reg16=reg12*reg11; T reg65=reg10*reg34; reg2=reg2*reg14; reg35=reg42-reg35;
    reg42=reg52*reg41; T reg66=reg40*reg9; T reg67=reg31*reg33; T reg68=reg10*reg30; T reg69=reg50*reg9;
    T reg70=reg23*reg32; T reg71=reg52*reg32; reg43=reg48-reg43; reg48=reg31*reg38; T reg72=reg40*reg41;
    reg36=reg53-reg36; reg53=reg20*reg8; T reg73=reg22*reg38; T reg74=reg25*reg9; T reg75=2*reg12;
    T reg76=reg30*reg33; reg51=reg49-reg51; reg49=reg34*reg9; T reg77=reg23*reg41; reg26=reg19*reg26;
    T reg78=reg30*reg38; T reg79=reg24*reg9; reg8=reg21*reg8; T reg80=reg7*reg15; T reg81=reg10*reg64;
    reg49=reg51-reg49; reg67=reg42+reg67; reg42=reg56*reg44; reg51=reg58*reg14; reg65=reg36-reg65;
    reg36=reg39*reg44; T reg82=reg59+reg60; T reg83=reg62-reg61; T reg84=reg75*reg11; reg73=reg70+reg73;
    reg70=reg39*reg35; reg54=reg77+reg54; reg69=reg43-reg69; reg43=reg56*reg35; reg48=reg71+reg48;
    reg76=reg72+reg76; reg71=reg34*reg44; reg72=reg53*reg9; reg74=reg28+reg74; reg28=reg34*reg35;
    reg78=reg55+reg78; reg79=reg26-reg79; reg25=reg19*reg25; reg26=reg8*reg9; reg55=reg10*reg58;
    reg24=reg21*reg24; reg77=reg57*reg5; T reg85=reg16*reg5; T reg86=reg10*reg63; T reg87=reg21*reg2;
    T reg88=reg34*reg14; reg68=reg66+reg68; reg2=reg20*reg2; reg58=reg58*reg9; reg66=reg5*reg3;
    T reg89=reg12*reg1; T reg90=reg57*reg85; reg43=reg48+reg43; reg71=reg76+reg71; reg36=reg54+reg36;
    reg68=reg88-reg68; reg48=reg89*reg3; reg59=reg60-reg59; reg54=reg80*reg85; reg60=reg84*reg77;
    reg28=reg78+reg28; reg76=reg22*reg65; reg70=reg73+reg70; reg73=reg84*reg85; reg42=reg67+reg42;
    reg67=reg57*reg77; reg78=reg23*reg49; reg88=reg30*reg65; T reg91=reg80*reg77; T reg92=reg40*reg49;
    T reg93=reg82*reg3; reg14=reg66*reg14; T reg94=reg10*reg58; T reg95=reg6*reg7; T reg96=reg10*reg55;
    T reg97=reg21*reg51; reg51=reg20*reg51; reg81=reg87+reg81; reg87=reg11*reg1; reg53=reg19*reg53;
    reg29=reg29/reg69; reg74=reg74/reg69; reg72=reg46+reg72; reg27=reg27/reg69; reg79=reg79/reg69;
    reg26=reg25+reg26; reg86=reg2-reg86; reg24=reg25+reg24; reg8=reg21*reg8; reg63=reg21*reg63;
    reg2=reg66*reg9; reg66=reg10*reg66; reg64=reg20*reg64; reg62=reg61+reg62; reg25=reg75*reg1;
    reg61=2*reg83; T reg98=reg15*reg4; T reg99=2*reg11; T reg100=reg34*reg68; T reg101=reg29*reg40;
    reg50=reg50/reg69; T reg102=reg99*reg1; T reg103=reg22*reg27; reg54=reg28+reg54; reg28=reg98*reg48;
    reg72=reg72/reg69; T reg104=reg29*reg30; T reg105=pow(reg83,2); T reg106=reg82*reg48; reg90=reg43+reg90;
    reg43=reg84*reg79; T reg107=reg80*reg74; T reg108=reg80*reg29; T reg109=reg84*reg27; reg55=reg21*reg55;
    reg58=reg20*reg58; reg8=reg53-reg8; reg5=reg95*reg5; reg53=reg74*reg40; reg24=reg24/reg69;
    reg86=reg19*reg86; T reg110=reg23*reg79; T reg111=reg13*reg4; T reg112=reg87*reg0; reg81=reg21*reg81;
    T reg113=reg62*reg0; reg26=reg26/reg69; T reg114=reg23*reg27; reg88=reg92+reg88; reg60=reg36+reg60;
    reg36=reg25*reg93; reg92=reg39*reg68; T reg115=reg10*reg66; T reg116=reg21*reg14; reg76=reg78+reg76;
    reg78=reg52*reg49; reg14=reg20*reg14; reg10=reg10*reg2; T reg117=reg31*reg65; reg94=reg97+reg94;
    reg97=reg61*reg59; T reg118=reg82*reg93; reg96=reg51-reg96; reg67=reg42+reg67; reg42=pow(reg59,2);
    reg73=reg70+reg73; reg51=reg22*reg79; reg70=reg74*reg30; T reg119=reg6*reg4; reg91=reg71+reg91;
    reg71=reg25*reg48; T reg120=reg98*reg93; T reg121=reg63+reg64; reg71=reg73+reg71; reg114=reg101+reg114;
    reg36=reg60+reg36; reg60=reg102*reg113; reg92=reg76+reg92; reg73=reg84*reg5; reg76=reg105*reg50;
    reg101=reg102*reg112; reg120=reg91+reg120; reg8=reg8/reg69; reg91=reg7*reg4; reg3=reg119*reg3;
    T reg122=reg111*reg113; reg106=reg90+reg106; reg90=reg62*reg112; reg117=reg78+reg117; reg78=reg56*reg68;
    reg28=reg54+reg28; reg54=reg62*reg113; T reg123=reg80*reg5; reg118=reg67+reg118; reg67=reg111*reg112;
    reg100=reg88+reg100; reg88=reg42*reg50; reg103=reg104+reg103; reg104=reg84*reg26; T reg124=reg80*reg72;
    T reg125=reg22*reg26; T reg126=reg72*reg30; reg10=reg116+reg10; reg116=reg23*reg26; T reg127=reg72*reg40;
    T reg128=reg97*reg24; reg43=reg107+reg43; reg107=reg42*reg24; reg115=reg14-reg115; reg51=reg70+reg51;
    reg14=reg105*reg24; reg110=reg53+reg110; reg121=reg121*reg9; reg96=reg19*reg96; reg94=reg21*reg94;
    reg53=reg55+reg58; reg109=reg108+reg109; reg81=reg86-reg81; reg2=reg20*reg2; reg20=reg97*reg50;
    reg66=reg21*reg66; reg70=reg42*reg8; reg125=reg126+reg125; reg122=reg120+reg122; reg67=reg28+reg67;
    reg28=reg105*reg8; reg116=reg127+reg116; reg86=reg57*reg5; reg101=reg71+reg101; reg10=reg21*reg10;
    reg76=reg114+reg76; reg90=reg106+reg90; reg128=reg43+reg128; reg107=reg51+reg107; reg21=reg25*reg3;
    reg121=reg81-reg121; reg115=reg19*reg115; reg60=reg36+reg60; reg73=reg92+reg73; reg14=reg110+reg14;
    reg53=reg53*reg9; reg88=reg103+reg88; reg20=reg109+reg20; reg78=reg117+reg78; reg54=reg118+reg54;
    reg94=reg96-reg94; reg123=reg100+reg123; reg19=reg98*reg3; reg36=reg97*reg8; reg104=reg124+reg104;
    reg43=reg66+reg2; reg0=reg91*reg0; reg51=reg54*reg101; reg71=reg54*reg67; reg36=reg104+reg36;
    reg81=reg95*reg76; reg92=reg16*reg14; reg96=reg59*reg83; reg100=reg7*reg83; reg10=reg115-reg10;
    reg70=reg125+reg70; reg103=reg95*reg88; reg104=reg16*reg107; reg106=reg16*reg128; reg108=reg95*reg20;
    reg121=reg121/reg69; reg9=reg43*reg9; reg43=reg122*reg90; reg28=reg116+reg28; reg109=reg111*reg0;
    reg19=reg123+reg19; reg86=reg78+reg86; reg78=reg82*reg3; reg110=reg102*reg0; reg53=reg94-reg53;
    reg21=reg73+reg21; reg73=reg6*reg59; reg94=reg60*reg90; reg114=reg23*reg107; reg73=reg100+reg73;
    reg100=reg96*reg36; reg115=reg59*reg12; reg106=reg108+reg106; reg108=reg83*reg11; reg78=reg86+reg78;
    reg86=reg88*reg40; reg53=reg53/reg69; reg116=reg76*reg30; reg117=reg22*reg14; reg118=reg57*reg121;
    reg120=reg31*reg121; reg123=reg52*reg121; reg88=reg88*reg30; reg107=reg22*reg107; reg9=reg10-reg9;
    reg109=reg19+reg109; reg10=reg20*reg30; reg19=reg22*reg128; reg76=reg76*reg40; reg14=reg23*reg14;
    reg110=reg21+reg110; reg43=reg71-reg43; reg21=reg67*reg60; reg71=reg96*reg70; reg104=reg103+reg104;
    reg103=reg122*reg101; reg124=reg7*reg59; reg125=reg6*reg83; reg126=reg96*reg28; reg92=reg81+reg92;
    reg81=reg62*reg0; reg94=reg51-reg94; reg20=reg20*reg40; reg128=reg23*reg128; reg69=reg9/reg69;
    reg9=elem.pos(2)[0]-elem.pos(0)[0]; reg51=elem.pos(1)[1]-elem.pos(0)[1]; reg127=reg94*reg109; T reg129=elem.pos(2)[1]-elem.pos(0)[1]; reg114=reg86+reg114;
    reg86=reg105*reg70; T reg130=reg43*reg110; T reg131=reg105*reg28; reg128=reg20+reg128; reg20=reg105*reg36;
    reg14=reg76+reg14; reg76=reg73*reg53; reg117=reg116+reg117; reg28=reg42*reg28; reg116=reg124*reg53;
    T reg132=reg125*reg53; reg107=reg88+reg107; reg70=reg42*reg70; reg19=reg10+reg19; reg36=reg42*reg36;
    reg10=reg83*reg12; reg88=reg59*reg11; reg115=reg108+reg115; reg108=elem.pos(1)[0]-elem.pos(0)[0]; T reg133=reg57*reg118;
    reg81=reg78+reg81; reg126=reg92+reg126; reg78=reg57*reg123; reg100=reg106+reg100; reg103=reg21-reg103;
    reg21=reg57*reg120; reg71=reg104+reg71; reg130=reg127-reg130; reg92=reg103*reg81; reg104=reg110*reg90;
    reg106=reg101*reg81; reg127=reg67*reg110; T reg134=reg109*reg101; T reg135=reg67*reg81; T reg136=reg109*reg90;
    T reg137=reg18*reg118; T reg138=reg9*reg51; reg20=reg128+reg20; reg128=reg129*reg108; T reg139=reg10*reg69;
    T reg140=reg88*reg69; T reg141=reg115*reg69; reg78=reg126+reg78; reg126=reg73*reg132; reg131=reg14+reg131;
    reg14=reg18*reg123; reg15=reg15*reg83; reg21=reg71+reg21; reg71=reg73*reg116; reg123=reg17*reg123;
    reg70=reg107+reg70; reg86=reg114+reg86; reg107=reg17*reg120; reg133=reg100+reg133; reg100=reg73*reg76;
    reg28=reg117+reg28; reg120=reg18*reg120; reg13=reg13*reg59; reg118=reg17*reg118; reg36=reg19+reg36;
    reg123=reg28+reg123; reg99=reg99*reg59; reg127=reg134-reg127; reg75=reg75*reg83; reg19=reg13*reg132;
    reg28=reg15*reg116; reg114=reg60*reg81; reg14=reg131+reg14; reg132=reg15*reg132; reg117=reg13*reg76;
    reg118=reg36+reg118; reg120=reg86+reg120; reg107=reg70+reg107; reg116=reg13*reg116; reg36=reg115*reg139;
    reg126=reg78+reg126; reg92=reg130+reg92; reg70=reg54*reg109; reg71=reg21+reg71; reg106=reg104-reg106;
    reg21=reg115*reg140; reg78=reg54*reg110; reg100=reg133+reg100; reg86=reg122*reg81; reg104=reg115*reg141;
    reg137=reg20+reg137; reg135=reg136-reg135; reg138=reg128-reg138; reg20=reg122*reg110; reg128=reg109*reg60;
    reg76=reg15*reg76; reg130=reg75*reg140; reg116=reg107+reg116; reg132=reg14+reg132; reg14=reg75*reg139;
    reg140=reg99*reg140; reg21=reg71+reg21; reg71=reg99*reg141; reg117=reg118+reg117; reg9=reg9/reg138;
    reg104=reg100+reg104; reg20=reg128-reg20; reg127=reg127/reg92; reg76=reg137+reg76; reg100=1-(*f.m).resolution;
    reg141=reg75*reg141; reg135=reg135/reg92; reg86=reg70-reg86; reg19=reg123+reg19; reg139=reg99*reg139;
    reg129=reg129/reg138; reg108=reg108/reg138; reg28=reg120+reg28; reg106=reg106/reg92; reg51=reg51/reg138;
    reg114=reg78-reg114; reg36=reg126+reg36; reg141=reg76+reg141; reg139=reg19+reg139; reg130=reg28+reg130;
    reg19=reg51-reg129; reg28=reg9-reg108; reg43=reg43/reg92; reg86=reg86/reg92; reg103=reg103/reg92;
    reg20=reg20/reg92; reg114=reg114/reg92; reg92=reg94/reg92; reg70=reg36*reg100; reg76=reg21*reg100;
    reg78=reg104*reg100; reg94=(*f.m).resolution*reg106; reg107=(*f.m).resolution*reg135; reg118=(*f.m).resolution*reg127; reg14=reg132+reg14;
    reg71=reg117+reg71; reg140=reg116+reg140; reg116=0.5*reg51; reg117=0.5*reg108; reg120=0.5*reg129;
    reg123=0.5*reg9; reg126=reg139*reg100; reg128=reg140*reg100; reg131=reg71*reg100; reg132=reg141*reg100;
    reg133=(*f.m).resolution*reg114; reg134=(*f.m).resolution*reg86; reg136=(*f.m).resolution*reg20; reg137=reg130*reg100; T reg142=reg14*reg100;
    T reg143=(*f.m).resolution*reg103; T reg144=(*f.m).resolution*reg92; T reg145=(*f.m).resolution*reg43; reg94=reg70+reg94; reg107=reg76-reg107;
    reg118=reg78+reg118; reg70=0.5*reg28; reg76=0.5*reg19; reg78=reg51*reg94; T reg146=reg9*reg107;
    T reg147=reg108*reg107; T reg148=reg116*reg118; T reg149=reg70*reg118; T reg150=reg19*reg94; T reg151=reg117*reg118;
    reg142=reg144+reg142; reg145=reg137-reg145; reg132=reg143+reg132; reg133=reg126-reg133; reg134=reg128+reg134;
    reg136=reg131-reg136; reg126=reg120*reg118; reg128=reg28*reg107; reg131=reg76*reg118; reg137=reg129*reg94;
    reg143=reg123*reg118; reg151=reg151-reg78; reg144=reg116*reg136; T reg152=reg108*reg134; T reg153=reg117*reg132;
    T reg154=reg51*reg142; reg126=reg126-reg146; T reg155=reg120*reg132; T reg156=reg9*reg145; reg149=reg150+reg149;
    reg150=reg70*reg132; T reg157=reg19*reg142; T reg158=reg123*reg136; T reg159=reg129*reg133; reg147=reg147-reg148;
    T reg160=reg116*reg132; T reg161=reg19*reg133; T reg162=reg70*reg136; T reg163=reg108*reg145; T reg164=reg28*reg134;
    T reg165=reg76*reg136; reg137=reg137-reg143; T reg166=reg117*reg136; T reg167=reg76*reg132; T reg168=reg28*reg145;
    T reg169=reg123*reg132; T reg170=reg129*reg142; reg131=reg128+reg131; reg128=reg51*reg133; T reg171=reg9*reg134;
    T reg172=reg120*reg136; reg150=reg157+reg150; reg163=reg163-reg160; reg172=reg172-reg171; reg166=reg166-reg128;
    reg149=2*reg149; reg147=2*reg147; reg167=reg168+reg167; reg162=reg161+reg162; reg165=reg164+reg165;
    reg126=2*reg126; reg137=2*reg137; reg153=reg153-reg154; reg152=reg152-reg144; reg170=reg170-reg169;
    reg151=2*reg151; reg159=reg159-reg158; reg155=reg155-reg156; reg131=2*reg131; reg157=reg70*reg131;
    reg161=reg51*reg155; reg164=reg117*reg126; reg168=reg120*reg149; T reg173=reg123*reg147; T reg174=reg9*reg162;
    T reg175=reg51*reg170; T reg176=reg117*reg137; T reg177=reg19*reg170; T reg178=reg129*reg163; T reg179=reg120*reg131;
    T reg180=reg51*reg167; T reg181=reg9*reg166; T reg182=reg120*reg151; T reg183=reg19*reg167; T reg184=reg28*reg159;
    T reg185=reg120*reg147; T reg186=reg9*reg152; T reg187=reg76*reg131; T reg188=reg120*reg126; T reg189=reg76*reg137;
    T reg190=reg9*reg172; T reg191=reg117*reg149; T reg192=reg70*reg137; T reg193=reg51*reg150; T reg194=reg9*reg159;
    T reg195=reg120*reg137; T reg196=reg9*reg165; T reg197=reg117*reg131; T reg198=reg70*reg147; reg167=reg129*reg167;
    T reg199=reg19*reg163; T reg200=reg123*reg149; T reg201=reg129*reg150; T reg202=reg117*reg151; T reg203=reg19*reg155;
    T reg204=reg76*reg147; T reg205=reg28*reg152; T reg206=reg108*reg172; T reg207=reg116*reg126; T reg208=reg70*reg126;
    T reg209=reg76*reg151; T reg210=reg123*reg131; T reg211=reg116*reg147; reg152=reg108*reg152; T reg212=reg19*reg153;
    reg131=reg116*reg131; T reg213=reg116*reg151; T reg214=reg108*reg166; T reg215=reg70*reg151; T reg216=reg51*reg153;
    reg150=reg19*reg150; reg147=reg117*reg147; reg163=reg51*reg163; reg151=reg123*reg151; reg153=reg129*reg153;
    reg172=reg28*reg172; T reg217=reg108*reg162; T reg218=reg116*reg149; T reg219=reg28*reg165; reg165=reg108*reg165;
    T reg220=reg76*reg126; T reg221=reg70*reg149; reg159=reg108*reg159; T reg222=reg116*reg137; reg149=reg76*reg149;
    reg162=reg28*reg162; reg126=reg123*reg126; reg155=reg129*reg155; reg166=reg28*reg166; reg137=reg123*reg137;
    reg170=reg129*reg170; reg220=reg172+reg220; reg189=reg184+reg189; reg187=reg219+reg187; reg194=reg195-reg194;
    reg174=reg168-reg174; reg173=reg178-reg173; reg151=reg153-reg151; reg126=reg155-reg126; reg196=reg179-reg196;
    reg149=reg162+reg149; reg198=reg199+reg198; reg192=reg177+reg192; reg208=reg203+reg208; reg180=reg197-reg180;
    reg215=reg212+reg215; reg216=reg202-reg216; reg163=reg147-reg163; reg161=reg164-reg161; reg218=reg217-reg218;
    reg175=reg176-reg175; reg131=reg165-reg131; reg222=reg159-reg222; reg137=reg170-reg137; reg193=reg191-reg193;
    reg210=reg167-reg210; reg186=reg185-reg186; reg200=reg201-reg200; reg204=reg205+reg204; reg181=reg182-reg181;
    reg209=reg166+reg209; reg207=reg206-reg207; reg190=reg188-reg190; reg150=reg221+reg150; reg211=reg152-reg211;
    reg213=reg214-reg213; reg157=reg183+reg157; reg151=reg138*reg151; reg216=reg138*reg216; reg207=reg138*reg207;
    reg163=reg138*reg163; reg126=reg138*reg126; reg218=reg138*reg218; reg131=reg138*reg131; reg149=reg138*reg149;
    reg213=reg138*reg213; reg222=reg138*reg222; reg137=reg138*reg137; reg198=reg138*reg198; reg211=reg138*reg211;
    reg210=reg138*reg210; reg215=reg138*reg215; reg192=reg138*reg192; reg200=reg138*reg200; reg208=reg138*reg208;
    reg204=reg138*reg204; reg209=reg138*reg209; reg220=reg138*reg220; reg150=reg138*reg150; reg157=reg138*reg157;
    reg189=reg138*reg189; reg190=reg138*reg190; reg181=reg138*reg181; reg187=reg138*reg187; reg186=reg138*reg186;
    reg194=reg138*reg194; reg193=reg138*reg193; reg196=reg138*reg196; reg180=reg180*reg138; reg173=reg138*reg173;
    reg175=reg138*reg175; reg161=reg138*reg161; reg174=reg138*reg174; T tmp_1_1=ponderation*reg187; T tmp_1_4=ponderation*reg209;
    T tmp_0_4=ponderation*reg215; T tmp_2_3=ponderation*reg126; T tmp_3_3=ponderation*reg190; T tmp_5_5=ponderation*reg211; T tmp_4_2=ponderation*reg175;
    T tmp_0_3=ponderation*reg208; T tmp_2_4=ponderation*reg151; T tmp_0_1=ponderation*reg157; T tmp_1_3=ponderation*reg220; T tmp_1_2=ponderation*reg189;
    T tmp_5_4=ponderation*reg213; T tmp_4_4=ponderation*reg216; T tmp_0_0=ponderation*reg150; T tmp_4_3=ponderation*reg161; T tmp_5_3=ponderation*reg207;
    T tmp_1_0=ponderation*reg149; T tmp_3_0=ponderation*reg174; T tmp_5_1=ponderation*reg131; T tmp_4_1=ponderation*reg180; T tmp_4_0=ponderation*reg193;
    T tmp_5_2=ponderation*reg222; T tmp_2_2=ponderation*reg137; T tmp_5_0=ponderation*reg218; T tmp_0_5=ponderation*reg198; T tmp_3_1=ponderation*reg196;
    T tmp_2_1=ponderation*reg210; T tmp_3_5=ponderation*reg186; T tmp_2_5=ponderation*reg173; T tmp_2_0=ponderation*reg200; T tmp_3_2=ponderation*reg194;
    T tmp_0_2=ponderation*reg192; T tmp_4_5=ponderation*reg163; T tmp_3_4=ponderation*reg181; T tmp_1_5=ponderation*reg204;
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
    T reg0=pow((*f.m).v1[1],2); T reg1=pow((*f.m).v1[0],2); T reg2=pow((*f.m).v2[1],2); T reg3=pow((*f.m).v2[0],2); T reg4=pow((*f.m).v1[2],2);
    reg0=reg1+reg0; reg2=reg3+reg2; reg4=reg0+reg4; reg0=pow((*f.m).v2[2],2); reg1=2*(*f.m).shear_modulus_23;
    reg3=2*(*f.m).shear_modulus_13; T reg5=2*(*f.m).shear_modulus_12; reg4=pow(reg4,0.5); reg3=1.0/reg3; reg1=1.0/reg1;
    reg0=reg2+reg0; reg5=1.0/reg5; reg2=reg3*reg1; reg0=pow(reg0,0.5); T reg6=(*f.m).v1[1]/reg4;
    T reg7=(*f.m).v1[0]/reg4; T reg8=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; reg4=(*f.m).v1[2]/reg4; T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=reg5*reg2;
    T reg11=(*f.m).v2[0]/reg0; T reg12=(*f.m).v2[1]/reg0; T reg13=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg14=2*reg6; T reg15=2*reg7;
    T reg16=reg10*reg8; T reg17=reg14*reg12; T reg18=reg15*reg11; T reg19=2*reg4; T reg20=reg13*reg10;
    T reg21=reg10*reg9; reg0=(*f.m).v2[2]/reg0; T reg22=pow(reg12,2); T reg23=1.0/(*f.m).elastic_modulus_1; T reg24=1.0/(*f.m).elastic_modulus_2;
    T reg25=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg26=pow(reg11,2); T reg27=pow(reg0,2); T reg28=reg22*reg24; T reg29=reg24*reg21;
    T reg30=reg13*reg20; T reg31=reg25*reg21; reg19=reg19*reg0; T reg32=reg17*reg25; T reg33=reg26*reg25;
    T reg34=reg17*reg24; T reg35=reg18*reg25; T reg36=reg26*reg23; T reg37=reg22*reg25; T reg38=reg18*reg23;
    T reg39=reg13*reg16; T reg40=reg25*reg20; T reg41=reg24*reg16; T reg42=reg17*reg13; T reg43=pow(reg7,2);
    reg30=reg29-reg30; reg29=pow(reg6,2); reg39=reg31+reg39; reg33=reg28-reg33; reg28=reg19*reg13;
    reg35=reg34-reg35; reg34=reg26*reg8; T reg44=reg27*reg13; T reg45=reg22*reg13; reg37=reg36-reg37;
    reg36=reg18*reg8; T reg46=reg27*reg8; T reg47=reg19*reg8; reg32=reg38-reg32; reg38=reg40+reg41;
    reg42=reg36+reg42; reg46=reg37-reg46; reg36=reg7*reg11; reg37=reg25*reg39; reg44=reg33-reg44;
    reg33=pow(reg4,2); reg19=reg19*reg9; T reg48=reg23*reg30; T reg49=reg6*reg12; reg28=reg35-reg28;
    reg47=reg32-reg47; reg32=reg23*reg43; reg35=reg7*reg12; T reg50=reg6*reg11; T reg51=reg25*reg43;
    T reg52=reg24*reg29; reg45=reg34+reg45; reg34=reg27*reg9; T reg53=reg25*reg29; T reg54=reg25*reg10;
    T reg55=reg4*reg11; T reg56=reg5*reg1; T reg57=reg13*reg29; T reg58=reg29*reg28; T reg59=reg7*reg0;
    T reg60=reg43*reg47; reg53=reg32-reg53; reg32=reg33*reg8; T reg61=reg13*reg2; T reg62=reg2*reg8;
    T reg63=reg43*reg8; T reg64=reg29*reg44; T reg65=reg43*reg46; T reg66=reg4*reg12; T reg67=reg6*reg0;
    reg2=reg2*reg9; T reg68=reg4*reg0; T reg69=reg50+reg35; reg51=reg52-reg51; reg42=reg19-reg42;
    reg19=reg36*reg47; reg52=2*reg11; T reg70=reg13*reg33; T reg71=reg22*reg44; reg45=reg34-reg45;
    reg34=reg26*reg47; T reg72=reg11*reg12; T reg73=reg22*reg28; T reg74=reg26*reg46; T reg75=reg36*reg46;
    T reg76=reg16*reg8; reg21=reg23*reg21; T reg77=reg20*reg8; reg10=reg24*reg10; T reg78=reg49*reg28;
    T reg79=reg38*reg8; reg37=reg48-reg37; reg48=reg49*reg44; T reg80=reg13*reg62; reg57=reg63+reg57;
    reg63=reg27*reg42; T reg81=reg33*reg9; T reg82=reg72*reg5; T reg83=reg69*reg5; T reg84=reg6*reg15;
    T reg85=reg13*reg61; T reg86=reg25*reg2; T reg87=reg56*reg9; T reg88=reg27*reg45; reg71=reg74+reg71;
    reg74=reg52*reg12; reg70=reg51-reg70; reg79=reg37-reg79; reg37=reg11*reg0; reg48=reg75+reg48;
    reg51=reg10*reg8; reg78=reg19+reg78; reg77=reg31+reg77; reg19=reg68*reg42; reg31=reg33*reg42;
    reg58=reg60+reg58; reg76=reg21-reg76; reg20=reg23*reg20; reg21=reg54*reg8; reg32=reg53-reg32;
    reg53=reg13*reg56; reg16=reg25*reg16; reg60=reg33*reg45; reg64=reg65+reg64; reg65=reg5*reg3;
    reg56=reg56*reg8; reg75=reg68*reg45; reg2=reg24*reg2; T reg89=reg67-reg66; reg73=reg34+reg73;
    reg34=reg59+reg55; T reg90=reg37*reg3; reg63=reg73+reg63; reg73=reg74*reg83; T reg91=reg74*reg82;
    reg88=reg71+reg88; reg71=reg22*reg70; T reg92=reg26*reg32; T reg93=reg84*reg83; reg31=reg58+reg31;
    reg58=reg84*reg82; reg60=reg64+reg60; reg64=reg69*reg82; reg75=reg48+reg75; reg59=reg55-reg59;
    reg48=reg15*reg4; reg55=reg29*reg70; T reg94=reg43*reg32; T reg95=reg34*reg3; reg80=reg86+reg80;
    reg57=reg81-reg57; reg81=reg24*reg87; reg87=reg25*reg87; reg76=reg76/reg79; reg85=reg2-reg85;
    reg2=reg13*reg53; reg39=reg39/reg79; reg51=reg40+reg51; reg86=reg65*reg8; reg77=reg77/reg79;
    T reg96=reg13*reg56; reg62=reg24*reg62; reg30=reg30/reg79; reg9=reg65*reg9; reg10=reg23*reg10;
    reg21=reg20+reg21; T reg97=2*reg89; reg67=reg66+reg67; reg19=reg78+reg19; reg66=reg69*reg83;
    reg78=reg7*reg6; T reg98=reg12*reg0; reg16=reg20+reg16; reg54=reg25*reg54; reg61=reg25*reg61;
    reg65=reg13*reg65; reg20=2*reg12; T reg99=reg52*reg0; T reg100=reg84*reg77; T reg101=reg74*reg76;
    T reg102=reg77*reg29; reg58=reg60+reg58; reg60=reg48*reg90; T reg103=reg22*reg76; reg21=reg21/reg79;
    T reg104=reg26*reg39; T reg105=reg30*reg43; reg38=reg38/reg79; reg93=reg31+reg93; reg31=reg48*reg95;
    reg51=reg51/reg79; reg71=reg92+reg71; reg92=reg27*reg57; T reg106=reg20*reg0; T reg107=pow(reg89,2);
    reg91=reg88+reg91; reg88=reg99*reg90; T reg108=reg99*reg95; reg73=reg63+reg73; reg63=reg7*reg4;
    T reg109=reg97*reg59; T reg110=pow(reg59,2); T reg111=reg49*reg70; T reg112=reg14*reg4; reg96=reg87+reg96;
    reg87=reg67*reg1; T reg113=reg98*reg1; T reg114=reg24*reg9; reg9=reg25*reg9; T reg115=reg13*reg65;
    reg5=reg78*reg5; T reg116=reg34*reg95; reg66=reg19+reg66; reg13=reg13*reg86; reg19=reg26*reg76;
    T reg117=reg77*reg43; T reg118=reg30*reg29; T reg119=reg74*reg39; T reg120=reg84*reg30; T reg121=reg22*reg39;
    reg16=reg16/reg79; reg64=reg75+reg64; reg75=reg34*reg90; reg54=reg10-reg54; reg53=reg25*reg53;
    reg56=reg24*reg56; reg85=reg23*reg85; reg80=reg25*reg80; reg10=reg61+reg62; T reg122=reg36*reg32;
    reg55=reg94+reg55; reg2=reg81-reg2; reg81=reg33*reg57; reg94=reg112*reg87; reg31=reg93+reg31;
    reg60=reg58+reg60; reg58=reg112*reg113; reg92=reg71+reg92; reg71=reg74*reg5; reg93=reg106*reg113;
    reg88=reg91+reg88; reg116=reg66+reg116; reg66=reg67*reg87; reg81=reg55+reg81; reg55=reg84*reg5;
    reg3=reg63*reg3; reg104=reg105+reg104; reg115=reg114-reg115; reg75=reg64+reg75; reg19=reg117+reg19;
    reg64=reg107*reg16; reg91=reg107*reg38; reg13=reg9+reg13; reg103=reg102+reg103; reg9=reg110*reg16;
    reg102=reg109*reg38; reg121=reg118+reg121; reg105=reg110*reg38; reg114=reg51*reg43; reg101=reg100+reg101;
    reg100=reg109*reg16; reg117=reg74*reg21; reg118=reg84*reg51; reg119=reg120+reg119; reg65=reg25*reg65;
    reg54=reg54/reg79; reg86=reg24*reg86; reg24=reg6*reg4; reg120=reg26*reg21; reg80=reg85-reg80;
    reg10=reg10*reg8; reg111=reg122+reg111; reg108=reg73+reg108; reg73=reg106*reg87; reg2=reg23*reg2;
    reg85=reg51*reg29; reg96=reg25*reg96; reg122=reg22*reg21; T reg123=reg53+reg56; T reg124=reg67*reg113;
    T reg125=reg68*reg57; reg120=reg114+reg120; reg100=reg101+reg100; reg58=reg60+reg58; reg9=reg103+reg9;
    reg64=reg19+reg64; reg19=reg69*reg5; reg125=reg111+reg125; reg1=reg24*reg1; reg60=reg48*reg3;
    reg55=reg81+reg55; reg124=reg75+reg124; reg73=reg108+reg73; reg123=reg123*reg8; reg96=reg2-reg96;
    reg66=reg116+reg66; reg115=reg23*reg115; reg91=reg104+reg91; reg2=reg109*reg54; reg13=reg25*reg13;
    reg10=reg80-reg10; reg23=reg99*reg3; reg25=reg65+reg86; reg105=reg121+reg105; reg71=reg92+reg71;
    reg117=reg118+reg117; reg75=reg110*reg54; reg122=reg85+reg122; reg102=reg119+reg102; reg94=reg31+reg94;
    reg93=reg88+reg93; reg31=reg107*reg54; reg10=reg10/reg79; reg80=reg66*reg58; reg81=reg78*reg105;
    reg85=reg72*reg9; reg88=reg94*reg124; reg92=reg66*reg93; reg101=reg73*reg124; reg123=reg96-reg123;
    reg96=reg78*reg102; reg13=reg115-reg13; reg8=reg25*reg8; reg25=reg72*reg100; reg19=reg125+reg19;
    reg103=reg34*reg3; reg104=reg7*reg59; reg108=reg6*reg89; reg2=reg117+reg2; reg111=reg112*reg1;
    reg60=reg55+reg60; reg31=reg120+reg31; reg55=reg106*reg1; reg23=reg71+reg23; reg71=reg78*reg91;
    reg114=reg72*reg64; reg75=reg122+reg75; reg115=reg59*reg89; reg116=reg91*reg29; reg55=reg23+reg55;
    reg8=reg13-reg8; reg13=reg36*reg10; reg23=reg22*reg64; reg104=reg108+reg104; reg108=reg26*reg100;
    reg117=reg102*reg43; reg118=reg49*reg10; reg119=reg69*reg10; reg120=reg6*reg59; reg123=reg123/reg79;
    reg121=reg7*reg89; reg88=reg80-reg88; reg80=reg58*reg73; reg122=reg26*reg9; reg125=reg105*reg43;
    reg91=reg91*reg43; reg64=reg26*reg64; T reg126=reg94*reg93; reg85=reg81+reg85; reg81=reg115*reg75;
    reg103=reg19+reg103; reg111=reg60+reg111; reg114=reg71+reg114; reg25=reg96+reg25; reg19=reg115*reg2;
    reg60=reg115*reg31; reg101=reg92-reg101; reg100=reg22*reg100; reg102=reg102*reg29; reg71=reg59*reg11;
    reg92=reg89*reg12; reg9=reg22*reg9; reg105=reg105*reg29; reg96=reg67*reg1; reg126=reg80-reg126;
    reg96=reg103+reg96; reg71=reg92+reg71; reg80=reg59*reg12; reg92=reg89*reg11; reg103=reg107*reg31;
    reg64=reg91+reg64; reg122=reg125+reg122; reg91=reg107*reg75; reg125=reg104*reg123; T reg127=reg120*reg123;
    T reg128=reg121*reg123; reg108=reg117+reg108; reg117=reg107*reg2; reg23=reg116+reg23; reg31=reg110*reg31;
    reg9=reg105+reg9; reg75=reg110*reg75; reg100=reg102+reg100; reg102=reg88*reg55; reg79=reg8/reg79;
    reg8=elem.pos(1)[0]-elem.pos(0)[0]; reg105=reg101*reg111; reg116=reg69*reg119; reg19=reg25+reg19; reg25=elem.pos(2)[0]-elem.pos(0)[0];
    T reg129=reg69*reg118; reg81=reg85+reg81; reg85=elem.pos(2)[1]-elem.pos(0)[1]; T reg130=elem.pos(1)[1]-elem.pos(0)[1]; T reg131=reg69*reg13;
    reg60=reg114+reg60; reg2=reg110*reg2; reg114=reg58*reg55; T reg132=reg111*reg93; T reg133=reg58*reg96;
    T reg134=reg111*reg124; T reg135=reg93*reg96; T reg136=reg71*reg79; T reg137=reg80*reg79; T reg138=reg92*reg79;
    T reg139=reg17*reg119; reg2=reg100+reg2; reg100=reg25*reg130; T reg140=reg85*reg8; T reg141=reg17*reg118;
    reg75=reg9+reg75; reg131=reg60+reg131; reg14=reg14*reg59; reg102=reg105-reg102; reg9=reg17*reg13;
    reg60=reg126*reg96; reg31=reg23+reg31; reg23=reg104*reg128; reg105=reg104*reg125; reg119=reg18*reg119;
    reg117=reg108+reg117; reg108=reg55*reg124; reg15=reg15*reg89; reg116=reg19+reg116; reg118=reg18*reg118;
    reg91=reg122+reg91; reg129=reg81+reg129; reg19=reg104*reg127; reg13=reg18*reg13; reg103=reg64+reg103;
    reg64=reg71*reg137; reg105=reg116+reg105; reg81=reg71*reg136; reg20=reg20*reg59; reg60=reg102+reg60;
    reg52=reg52*reg89; reg102=reg66*reg55; reg100=reg140-reg100; reg19=reg129+reg19; reg116=reg71*reg138;
    reg23=reg131+reg23; reg122=reg14*reg125; reg139=reg2+reg139; reg2=reg14*reg127; reg141=reg75+reg141;
    reg75=reg14*reg128; reg9=reg31+reg9; reg125=reg15*reg125; reg119=reg117+reg119; reg127=reg15*reg127;
    reg118=reg91+reg118; reg128=reg15*reg128; reg13=reg103+reg13; reg31=reg73*reg96; reg91=reg66*reg111;
    reg133=reg134-reg133; reg103=reg111*reg73; reg135=reg108-reg135; reg108=reg94*reg55; reg117=reg94*reg96;
    reg114=reg132-reg114; reg129=reg52*reg137; reg125=reg119+reg125; reg119=reg52*reg136; reg108=reg103-reg108;
    reg114=reg114/reg60; reg75=reg9+reg75; reg9=reg20*reg138; reg25=reg25/reg100; reg130=reg130/reg100;
    reg103=1-(*f.m).resolution; reg2=reg141+reg2; reg137=reg20*reg137; reg122=reg139+reg122; reg136=reg20*reg136;
    reg85=reg85/reg100; reg8=reg8/reg100; reg81=reg105+reg81; reg64=reg19+reg64; reg116=reg23+reg116;
    reg138=reg52*reg138; reg117=reg91-reg117; reg135=reg135/reg60; reg133=reg133/reg60; reg127=reg118+reg127;
    reg31=reg102-reg31; reg128=reg13+reg128; reg136=reg122+reg136; reg9=reg75+reg9; reg101=reg101/reg60;
    reg31=reg31/reg60; reg13=reg25-reg8; reg19=reg130-reg85; reg88=reg88/reg60; reg23=reg116*reg103;
    reg75=reg64*reg103; reg91=reg81*reg103; reg138=reg128+reg138; reg102=(*f.m).resolution*reg114; reg105=(*f.m).resolution*reg133;
    reg117=reg117/reg60; reg119=reg125+reg119; reg118=(*f.m).resolution*reg135; reg137=reg2+reg137; reg108=reg108/reg60;
    reg60=reg126/reg60; reg129=reg127+reg129; reg105=reg75-reg105; reg118=reg23+reg118; reg2=reg9*reg103;
    reg102=reg91+reg102; reg23=reg137*reg103; reg75=0.5*reg25; reg91=0.5*reg85; reg122=reg136*reg103;
    reg125=0.5*reg8; reg126=0.5*reg130; reg127=(*f.m).resolution*reg108; reg128=(*f.m).resolution*reg117; reg131=(*f.m).resolution*reg31;
    reg132=0.5*reg13; reg134=0.5*reg19; reg139=(*f.m).resolution*reg88; reg140=(*f.m).resolution*reg101; reg141=(*f.m).resolution*reg60;
    T reg142=reg119*reg103; T reg143=reg129*reg103; T reg144=reg138*reg103; T reg145=reg126*reg102; reg128=reg23+reg128;
    reg23=reg8*reg105; reg127=reg122-reg127; reg122=reg75*reg102; T reg146=reg134*reg102; reg131=reg2-reg131;
    reg142=reg141+reg142; reg139=reg143-reg139; reg144=reg140+reg144; reg2=reg13*reg105; reg140=reg85*reg118;
    reg141=reg132*reg102; reg143=reg125*reg102; T reg147=reg19*reg118; T reg148=reg25*reg105; T reg149=reg91*reg102;
    T reg150=reg130*reg118; T reg151=reg132*reg142; T reg152=reg134*reg142; T reg153=reg13*reg128; T reg154=reg19*reg144;
    reg140=reg140-reg122; T reg155=reg134*reg127; T reg156=reg25*reg139; T reg157=reg91*reg142; T reg158=reg125*reg142;
    reg149=reg149-reg148; T reg159=reg130*reg144; T reg160=reg130*reg131; T reg161=reg85*reg144; T reg162=reg75*reg142;
    reg143=reg143-reg150; T reg163=reg125*reg127; T reg164=reg91*reg127; T reg165=reg25*reg128; T reg166=reg126*reg127;
    T reg167=reg8*reg128; T reg168=reg8*reg139; T reg169=reg126*reg142; T reg170=reg75*reg127; T reg171=reg85*reg131;
    reg23=reg23-reg145; reg141=reg147+reg141; reg147=reg13*reg139; reg146=reg2+reg146; reg168=reg168-reg169;
    reg23=2*reg23; reg152=reg147+reg152; reg146=2*reg146; reg161=reg161-reg162; reg140=2*reg140;
    reg158=reg158-reg159; reg149=2*reg149; reg157=reg157-reg156; reg167=reg167-reg166; reg143=2*reg143;
    reg163=reg163-reg160; reg164=reg164-reg165; reg171=reg171-reg170; reg155=reg153+reg155; reg151=reg154+reg151;
    reg141=2*reg141; reg2=reg130*reg168; reg147=reg125*reg23; reg153=reg132*reg23; reg154=reg19*reg157;
    T reg172=reg130*reg158; T reg173=reg75*reg143; T reg174=reg25*reg167; T reg175=reg91*reg23; T reg176=reg132*reg149;
    T reg177=reg25*reg163; T reg178=reg91*reg143; T reg179=reg85*reg158; reg158=reg19*reg158; T reg180=reg25*reg164;
    T reg181=reg132*reg140; T reg182=reg19*reg151; T reg183=reg19*reg161; T reg184=reg132*reg146; T reg185=reg13*reg155;
    T reg186=reg19*reg152; T reg187=reg134*reg146; T reg188=reg19*reg168; T reg189=reg13*reg171; T reg190=reg134*reg140;
    T reg191=reg91*reg149; T reg192=reg13*reg164; T reg193=reg134*reg149; T reg194=reg13*reg163; T reg195=reg75*reg149;
    T reg196=reg85*reg157; T reg197=reg132*reg141; T reg198=reg132*reg143; T reg199=reg75*reg140; T reg200=reg85*reg161;
    T reg201=reg125*reg143; T reg202=reg75*reg23; T reg203=reg134*reg23; T reg204=reg13*reg167; reg168=reg85*reg168;
    T reg205=reg134*reg143; reg23=reg126*reg23; reg167=reg8*reg167; reg176=reg154+reg176; reg198=reg158+reg198;
    reg181=reg183+reg181; reg193=reg192+reg193; reg2=reg147-reg2; reg23=reg167-reg23; reg190=reg189+reg190;
    reg172=reg201-reg172; reg187=reg185+reg187; reg174=reg175-reg174; reg199=reg200-reg199; reg202=reg168-reg202;
    reg177=reg178-reg177; reg173=reg179-reg173; reg205=reg194+reg205; reg180=reg191-reg180; reg195=reg196-reg195;
    reg203=reg204+reg203; reg184=reg186+reg184; reg153=reg188+reg153; reg182=reg197+reg182; reg23=reg100*reg23;
    reg205=reg100*reg205; reg181=reg100*reg181; reg2=reg100*reg2; reg172=reg100*reg172; reg174=reg100*reg174;
    reg177=reg100*reg177; reg176=reg100*reg176; reg180=reg100*reg180; reg184=reg100*reg184; reg182=reg100*reg182;
    reg198=reg100*reg198; reg153=reg100*reg153; reg195=reg100*reg195; reg173=reg100*reg173; reg202=reg100*reg202;
    reg187=reg100*reg187; reg190=reg100*reg190; reg193=reg100*reg193; reg203=reg100*reg203; reg199=reg100*reg199;
    T tmp_1_4=ponderation*reg205; T tmp_0_4=ponderation*reg198; T tmp_0_0=ponderation*reg182; T tmp_0_5=ponderation*reg153; T tmp_0_1=ponderation*reg184;
    T tmp_5_5=ponderation*reg23; T tmp_2_3=ponderation*reg195; T tmp_3_3=ponderation*reg180; T tmp_1_5=ponderation*reg203; T tmp_0_3=ponderation*reg176;
    T tmp_2_4=ponderation*reg173; T tmp_3_4=ponderation*reg177; T tmp_2_5=ponderation*reg202; T tmp_2_2=ponderation*reg199; T tmp_3_5=ponderation*reg174;
    T tmp_1_1=ponderation*reg187; T tmp_4_5=ponderation*reg2; T tmp_4_4=ponderation*reg172; T tmp_1_2=ponderation*reg190; T tmp_1_3=ponderation*reg193;
    T tmp_0_2=ponderation*reg181;
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
    T reg0=2*(*f.m).shear_modulus_13; T reg1=2*(*f.m).shear_modulus_23; T reg2=2*(*f.m).shear_modulus_12; reg0=1.0/reg0; reg1=1.0/reg1;
    reg2=1.0/reg2; T reg3=reg0*reg1; T reg4=reg2*reg3; T reg5=1.0/(*f.m).elastic_modulus_3; T reg6=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg7=pow((*f.m).v2[0],2); T reg8=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg9=pow((*f.m).v2[1],2); T reg10=pow((*f.m).v1[0],2); T reg11=pow((*f.m).v1[1],2);
    T reg12=1.0/(*f.m).elastic_modulus_2; T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg14=reg4*reg6; T reg15=pow((*f.m).v2[2],2); reg9=reg7+reg9;
    reg11=reg10+reg11; reg7=reg4*reg5; reg10=pow((*f.m).v1[2],2); T reg16=reg8*reg4; T reg17=reg8*reg14;
    T reg18=reg13*reg7; reg10=reg11+reg10; reg11=reg8*reg16; T reg19=reg12*reg7; reg15=reg9+reg15;
    reg15=pow(reg15,0.5); reg11=reg19-reg11; reg10=pow(reg10,0.5); reg9=1.0/(*f.m).elastic_modulus_1; reg17=reg18+reg17;
    reg19=reg13*reg16; T reg20=reg12*reg14; T reg21=(*f.m).v1[1]/reg10; T reg22=(*f.m).v2[1]/reg15; T reg23=(*f.m).v2[2]/reg15;
    T reg24=(*f.m).v1[2]/reg10; T reg25=reg9*reg11; T reg26=reg19+reg20; T reg27=reg13*reg17; T reg28=reg3*reg5;
    reg10=(*f.m).v1[0]/reg10; T reg29=reg21*reg23; reg27=reg25-reg27; reg25=reg3*reg6; reg15=(*f.m).v2[0]/reg15;
    reg3=reg8*reg3; T reg30=reg24*reg22; T reg31=reg16*reg6; T reg32=reg26*reg6; T reg33=reg2*reg1;
    T reg34=reg12*reg4; reg7=reg9*reg7; reg4=reg13*reg4; T reg35=reg14*reg6; T reg36=2*reg10;
    T reg37=reg4*reg6; reg14=reg13*reg14; T reg38=reg8*reg33; T reg39=reg2*reg0; T reg40=reg33*reg6;
    T reg41=reg10*reg23; T reg42=reg24*reg15; T reg43=reg34*reg6; T reg44=2*reg15; reg33=reg33*reg5;
    reg32=reg27-reg32; reg27=reg8*reg25; reg35=reg7-reg35; reg16=reg9*reg16; reg7=reg8*reg3;
    T reg45=reg13*reg28; T reg46=reg29-reg30; reg31=reg18+reg31; reg28=reg12*reg28; reg34=reg9*reg34;
    reg31=reg31/reg32; reg17=reg17/reg32; reg43=reg19+reg43; reg11=reg11/reg32; reg35=reg35/reg32;
    reg18=reg10*reg22; T reg47=reg21*reg15; T reg48=reg39*reg5; T reg49=2*reg46; T reg50=reg8*reg40;
    T reg51=reg42-reg41; T reg52=reg8*reg38; T reg53=reg13*reg33; reg33=reg12*reg33; reg27=reg45+reg27;
    reg45=reg44*reg22; reg7=reg28-reg7; reg28=pow(reg22,2); T reg54=pow(reg15,2); T reg55=reg39*reg6;
    T reg56=pow(reg21,2); reg25=reg12*reg25; reg39=reg8*reg39; reg3=reg13*reg3; T reg57=reg21*reg36;
    reg4=reg13*reg4; reg14=reg16+reg14; T reg58=pow(reg10,2); reg37=reg16+reg37; reg16=reg31*reg58;
    T reg59=reg54*reg35; T reg60=reg11*reg56; T reg61=reg28*reg17; T reg62=reg57*reg11; T reg63=reg45*reg17;
    T reg64=pow(reg24,2); T reg65=reg54*reg17; reg37=reg37/reg32; reg14=reg14/reg32; reg4=reg34-reg4;
    reg38=reg13*reg38; reg40=reg12*reg40; reg34=reg11*reg58; reg7=reg9*reg7; reg27=reg13*reg27;
    T reg66=reg3+reg25; reg52=reg33-reg52; reg50=reg53+reg50; reg33=reg12*reg48; reg48=reg13*reg48;
    reg53=reg8*reg39; T reg67=reg18-reg47; T reg68=reg8*reg55; T reg69=reg57*reg31; T reg70=reg45*reg35;
    T reg71=reg49*reg51; T reg72=reg28*reg35; T reg73=reg31*reg56; T reg74=pow(reg51,2); T reg75=pow(reg46,2);
    T reg76=pow(reg23,2); reg43=reg43/reg32; reg26=reg26/reg32; T reg77=reg43*reg56; reg55=reg12*reg55;
    T reg78=reg28*reg37; reg27=reg7-reg27; reg66=reg66*reg6; reg7=reg57*reg43; T reg79=reg54*reg37;
    T reg80=reg45*reg37; reg52=reg9*reg52; reg50=reg13*reg50; T reg81=reg38+reg40; reg53=reg33-reg53;
    reg68=reg48+reg68; reg61=reg60+reg61; reg33=reg74*reg26; reg48=reg11*reg64; reg60=reg76*reg17;
    reg63=reg62+reg63; reg62=reg71*reg26; T reg82=pow(reg67,2); T reg83=reg74*reg14; reg70=reg69+reg70;
    reg72=reg73+reg72; reg69=reg31*reg64; reg4=reg4/reg32; reg73=reg71*reg14; T reg84=reg75*reg14;
    reg59=reg16+reg59; reg65=reg34+reg65; reg16=reg75*reg26; reg34=reg76*reg35; T reg85=reg43*reg58;
    reg39=reg13*reg39; T reg86=reg15*reg22; T reg87=reg39+reg55; reg68=reg13*reg68; T reg88=reg10*reg21;
    reg83=reg72+reg83; reg53=reg9*reg53; reg33=reg61+reg33; reg60=reg48+reg60; reg48=reg26*reg82;
    reg84=reg59+reg84; reg59=2*reg22; reg62=reg63+reg62; reg61=reg44*reg23; reg16=reg65+reg16;
    reg63=2*reg21; reg65=reg36*reg24; reg66=reg27-reg66; reg27=reg76*reg37; reg72=reg43*reg64;
    T reg89=reg14*reg82; reg80=reg7+reg80; reg7=reg71*reg4; reg50=reg52-reg50; reg52=reg74*reg4;
    reg73=reg70+reg73; reg78=reg77+reg78; reg34=reg69+reg34; reg81=reg81*reg6; reg69=reg75*reg4;
    reg79=reg85+reg79; reg70=reg65*reg31; reg77=reg61*reg35; reg69=reg79+reg69; reg52=reg78+reg52;
    reg78=reg88*reg33; reg79=reg86*reg83; reg85=reg86*reg73; T reg90=reg88*reg62; reg27=reg72+reg27;
    reg72=reg4*reg82; reg89=reg34+reg89; reg34=reg28*reg73; T reg91=reg62*reg56; reg7=reg80+reg7;
    reg80=reg28*reg83; T reg92=reg16*reg58; T reg93=reg33*reg56; T reg94=reg28*reg84; T reg95=reg16*reg56;
    T reg96=reg54*reg84; reg73=reg54*reg73; reg62=reg62*reg58; reg83=reg54*reg83; reg33=reg33*reg58;
    T reg97=reg21*reg46; T reg98=reg10*reg51; T reg99=reg61*reg17; T reg100=reg65*reg11; reg48=reg60+reg48;
    reg87=reg87*reg6; reg68=reg53-reg68; reg81=reg50-reg81; reg66=reg66/reg32; reg50=2*reg51;
    reg49=reg49*reg67; reg53=reg63*reg24; reg60=reg36*reg15; T reg101=reg63*reg22; T reg102=2*reg24;
    T reg103=reg51*reg46; T reg104=reg59*reg23; T reg105=reg21*reg22; T reg106=reg10*reg15; reg18=reg47+reg18;
    reg47=reg51*reg15; T reg107=reg24*reg23; T reg108=reg21*reg51; reg34=reg91+reg34; reg91=reg74*reg7;
    reg96=reg92+reg96; reg35=reg104*reg35; reg31=reg53*reg31; reg92=reg54*reg89; T reg109=reg49*reg14;
    reg77=reg70+reg77; reg70=reg75*reg69; T reg110=reg10*reg46; T reg111=reg48*reg58; T reg112=reg54*reg9;
    reg102=reg102*reg23; T reg113=reg28*reg13; T reg114=reg18*reg66; T reg115=reg105*reg66; T reg116=reg106*reg66;
    T reg117=reg75*reg52; T reg118=reg61*reg37; T reg119=reg65*reg43; T reg120=reg74*reg52; reg80=reg93+reg80;
    reg50=reg50*reg67; reg93=reg54*reg13; reg98=reg97+reg98; reg97=reg48*reg56; reg83=reg33+reg83;
    reg33=reg74*reg69; reg72=reg27+reg72; reg94=reg95+reg94; reg27=reg28*reg89; reg95=reg46*reg22;
    T reg121=reg28*reg12; T reg122=reg60*reg9; reg99=reg100+reg99; reg81=reg81/reg32; reg100=reg75*reg7;
    T reg123=reg49*reg26; reg7=reg103*reg7; reg79=reg78+reg79; reg85=reg90+reg85; reg87=reg68-reg87;
    reg68=reg101*reg13; reg11=reg53*reg11; reg78=reg60*reg13; reg73=reg62+reg73; reg17=reg104*reg17;
    reg62=reg101*reg12; reg84=reg86*reg84; reg16=reg88*reg16; reg52=reg103*reg52; reg7=reg85+reg7;
    reg48=reg88*reg48; reg89=reg86*reg89; reg85=reg107*reg66; reg90=reg102*reg6; reg120=reg80+reg120;
    reg80=reg101*reg115; reg37=reg104*reg37; reg43=reg53*reg43; reg117=reg83+reg117; reg113=reg112-reg113;
    reg83=reg76*reg6; reg26=reg50*reg26; reg17=reg11+reg17; reg11=reg60*reg116; reg70=reg96+reg70;
    reg123=reg99+reg123; reg96=reg101*reg8; reg99=reg67*reg24; reg112=reg60*reg6; reg33=reg94+reg33;
    reg94=reg101*reg116; T reg124=reg28*reg8; T reg125=reg98*reg81; T reg126=reg54*reg6; T reg127=reg108*reg81;
    T reg128=reg110*reg81; reg14=reg50*reg14; reg35=reg31+reg35; reg31=reg18*reg114; reg109=reg77+reg109;
    reg100=reg73+reg100; reg91=reg34+reg91; reg34=reg101*reg114; reg114=reg60*reg114; reg36=reg36*reg46;
    reg68=reg122-reg68; reg73=reg75*reg72; reg63=reg63*reg51; reg69=reg103*reg69; reg84=reg16+reg84;
    reg92=reg111+reg92; reg16=reg49*reg4; reg118=reg119+reg118; reg78=reg62-reg78; reg62=reg102*reg8;
    reg77=reg18*reg115; reg115=reg60*reg115; reg32=reg87/reg32; reg27=reg97+reg27; reg87=reg74*reg72;
    reg97=reg46*reg15; reg111=reg51*reg22; reg52=reg79+reg52; reg47=reg95+reg47; reg93=reg121-reg93;
    reg79=reg76*reg8; reg95=reg123*reg58; reg114=reg100+reg114; reg83=reg113-reg83; reg100=reg36*reg125;
    reg113=reg54*reg109; reg72=reg103*reg72; reg89=reg48+reg89; reg48=reg98*reg127; reg77=reg52+reg77;
    reg31=reg7+reg31; reg7=reg98*reg125; reg116=reg18*reg116; reg69=reg84+reg69; reg52=reg54*(*f.m).alpha_2;
    reg84=reg28*reg109; reg119=reg123*reg56; reg125=reg63*reg125; reg34=reg91+reg34; reg91=(*f.m).alpha_1*reg56;
    reg121=reg101*reg85; reg87=reg27+reg87; reg27=reg28*(*f.m).alpha_2; reg122=reg63*reg127; reg80=reg120+reg80;
    reg120=reg63*reg128; reg94=reg33+reg94; reg33=reg13*reg56; reg13=reg13*reg58; reg9=reg9*reg58;
    reg12=reg12*reg56; T reg129=reg47*reg32; T reg130=reg111*reg32; T reg131=reg97*reg32; T reg132=reg99*reg81;
    reg4=reg50*reg4; reg37=reg43+reg37; reg16=reg118+reg16; reg43=reg46*reg24; reg118=reg10*reg67;
    T reg133=reg67*reg23; reg14=reg35+reg14; reg44=reg44*reg46; reg59=reg59*reg51; reg35=(*f.m).alpha_1*reg58;
    reg96=reg112+reg96; reg102=reg102*reg5; reg79=reg93-reg79; reg62=reg78-reg62; reg78=reg76*reg5;
    reg124=reg126+reg124; reg26=reg17+reg26; reg42=reg41+reg42; reg115=reg117+reg115; reg90=reg68-reg90;
    reg127=reg36*reg127; reg11=reg70+reg11; reg17=reg60*reg85; reg73=reg92+reg73; reg41=reg36*reg128;
    reg68=reg105*reg62; reg70=reg58*reg6; reg48=reg77+reg48; reg77=reg47*reg130; reg92=reg74*reg16;
    reg84=reg119+reg84; reg93=reg8*reg56; reg112=reg44*reg129; reg35=reg52+reg35; reg52=reg59*reg129;
    reg125=reg34+reg125; reg34=reg54*reg90; reg117=reg28*reg62; reg119=reg106*reg83; reg126=reg67*reg15;
    T reg134=reg56*reg79; reg128=reg98*reg128; reg116=reg69+reg116; reg69=reg58*reg83; reg100=reg114+reg100;
    reg114=reg58*reg90; T reg135=reg56*reg62; reg7=reg31+reg7; reg129=reg47*reg129; reg123=reg88*reg123;
    reg109=reg86*reg109; reg127=reg115+reg127; reg31=reg44*reg130; reg115=reg28*reg14; T reg136=reg26*reg56;
    T reg137=reg54*reg83; T reg138=reg28*reg79; T reg139=reg106*reg90; T reg140=reg42*reg66; reg17=reg73+reg17;
    reg72=reg89+reg72; reg13=reg12-reg13; reg27=reg91+reg27; reg12=reg59*reg131; reg120=reg94+reg120;
    reg33=reg9-reg33; reg9=reg133*reg32; reg6=reg64*reg6; reg73=reg26*reg58; reg8=reg8*reg64;
    reg89=reg54*reg14; reg91=(*f.m).alpha_3*reg74; reg85=reg18*reg85; reg94=(*f.m).alpha_1*reg64; T reg141=reg76*(*f.m).alpha_2;
    T reg142=reg36*reg132; reg96=reg102-reg96; reg102=reg46*reg23; T reg143=reg63*reg132; reg124=reg78-reg124;
    reg78=(*f.m).alpha_3*reg75; reg121=reg87+reg121; reg29=reg30+reg29; reg30=reg105*reg79; reg87=reg21*reg67;
    T reg144=reg51*reg24; reg118=reg43+reg118; reg113=reg95+reg113; reg43=reg75*reg16; reg130=reg59*reg130;
    reg122=reg80+reg122; reg80=reg44*reg131; reg41=reg11+reg41; reg4=reg37+reg4; reg11=reg107*reg124;
    reg93=reg70+reg93; reg5=reg64*reg5; reg37=reg18*reg2; reg30=reg119+reg30; reg70=reg86*reg2;
    reg95=reg86*(*f.m).alpha_2; reg119=(*f.m).alpha_1*reg88; T reg145=(*f.m).alpha_3*reg82; reg141=reg94+reg141; reg91=reg27+reg91;
    reg6=reg33-reg6; reg8=reg13-reg8; reg78=reg35+reg78; reg13=reg76*reg96; reg117=reg34+reg117;
    reg27=reg76*reg124; reg138=reg137+reg138; reg33=reg107*reg96; reg68=reg139+reg68; reg14=reg86*reg14;
    reg26=reg88*reg26; reg16=reg103*reg16; reg109=reg123+reg109; reg129=reg7+reg129; reg7=reg64*reg96;
    reg135=reg114+reg135; reg34=reg64*reg124; reg134=reg69+reg134; reg43=reg113+reg43; reg112=reg100+reg112;
    reg35=reg44*reg9; reg142=reg17+reg142; reg17=reg74*reg4; reg115=reg136+reg115; reg10=reg10*reg24;
    reg52=reg125+reg52; reg31=reg127+reg31; reg80=reg41+reg80; reg92=reg84+reg92; reg126=reg102+reg126;
    reg41=reg118*reg81; reg66=reg29*reg66; reg87=reg144+reg87; reg15=reg15*reg23; reg69=reg101*reg140;
    reg132=reg98*reg132; reg85=reg72+reg85; reg77=reg48+reg77; reg131=reg47*reg131; reg128=reg116+reg128;
    reg48=reg51*reg23; reg72=reg75*reg4; reg89=reg73+reg89; reg12=reg120+reg12; reg73=reg60*reg140;
    reg84=reg67*reg22; reg130=reg122+reg130; reg143=reg121+reg143; reg86=reg59*reg9; reg94=reg54*reg6;
    reg100=reg57*reg37; reg7=reg135+reg7; reg60=reg60*reg66; reg23=reg22*reg23; reg22=reg15*(*f.m).alpha_2;
    reg145=reg141+reg145; reg102=(*f.m).alpha_1*reg10; reg72=reg89+reg72; reg95=reg119+reg95; reg89=(*f.m).alpha_3*reg103;
    reg113=reg112*reg77; reg114=reg52*reg77; reg116=reg31*reg129; reg13=reg117+reg13; reg117=reg45*reg37;
    reg119=reg130*reg129; reg84=reg48+reg84; reg48=reg45*reg70; reg120=reg80*reg78; reg27=reg138+reg27;
    reg121=reg31*reg91; reg122=reg12*reg78; reg123=reg130*reg91; reg11=reg30+reg11; reg30=reg18*reg70;
    reg46=reg67*reg46; reg24=reg21*reg24; reg21=reg28*reg8; reg69=reg92+reg69; reg92=reg63*reg41;
    reg125=reg42*reg0; reg15=reg15*reg0; reg93=reg5-reg93; reg81=reg87*reg81; reg5=reg126*reg32;
    reg4=reg103*reg4; reg14=reg26+reg14; reg35=reg142+reg35; reg140=reg18*reg140; reg16=reg109+reg16;
    reg17=reg115+reg17; reg101=reg101*reg66; reg73=reg43+reg73; reg26=reg36*reg41; reg131=reg128+reg131;
    reg132=reg85+reg132; reg9=reg47*reg9; reg43=reg56*reg8; reg85=reg58*reg6; reg103=reg57*reg70;
    reg33=reg68+reg33; reg68=reg18*reg37; reg34=reg134+reg34; reg86=reg143+reg86; reg48=reg27+reg48;
    reg27=reg61*reg15; reg51=reg67*reg51; reg67=reg112*reg130; reg109=reg31*reg52; reg64=reg64*reg93;
    reg26=reg73+reg26; reg73=reg44*reg5; reg103=reg34+reg103; reg113=reg116-reg113; reg34=reg65*reg15;
    reg43=reg85+reg43; reg85=reg29*reg1; reg115=reg23*reg1; reg121=reg120+reg121; reg35=reg35*reg145;
    reg32=reg84*reg32; reg100=reg7+reg100; reg123=reg122+reg123; reg7=reg65*reg125; reg2=reg88*reg2;
    reg86=reg86*reg145; reg21=reg94+reg21; reg76=reg76*reg93; reg88=reg77*reg91; reg94=reg131*reg78;
    reg30=reg11+reg30; reg11=reg42*reg15; reg23=reg23*(*f.m).alpha_2; reg116=(*f.m).alpha_1*reg24; reg46=(*f.m).alpha_3*reg46;
    reg120=reg105*reg8; reg22=reg102+reg22; reg89=reg95+reg89; reg95=reg106*reg6; reg60=reg72+reg60;
    reg36=reg36*reg81; reg72=reg42*reg125; reg68=reg33+reg68; reg92=reg69+reg92; reg33=reg59*reg5;
    reg9=reg132+reg9; reg114=reg119-reg114; reg117=reg13+reg117; reg13=reg61*reg125; reg63=reg63*reg81;
    reg101=reg17+reg101; reg140=reg16+reg140; reg41=reg98*reg41; reg66=reg18*reg66; reg4=reg14+reg4;
    reg73=reg26+reg73; reg7=reg100+reg7; reg14=reg53*reg85; reg23=reg116+reg23; reg51=(*f.m).alpha_3*reg51;
    reg16=reg29*reg115; reg9=reg9*reg145; reg88=reg94+reg88; reg76=reg21+reg76; reg17=reg45*reg2;
    reg11=reg30+reg11; reg81=reg98*reg81; reg21=reg52*reg89; reg86=reg123+reg86; reg66=reg4+reg66;
    reg59=reg59*reg32; reg46=reg22+reg46; reg33=reg92+reg33; reg63=reg101+reg63; reg0=reg10*reg0;
    reg4=reg53*reg115; reg34=reg103+reg34; reg5=reg47*reg5; reg41=reg140+reg41; reg36=reg60+reg36;
    reg44=reg44*reg32; reg10=reg57*reg2; reg64=reg43+reg64; reg22=reg29*reg85; reg72=reg68+reg72;
    reg26=reg12*reg113; reg30=reg80*reg114; reg67=reg109-reg67; reg13=reg117+reg13; reg43=reg104*reg85;
    reg60=reg104*reg115; reg27=reg48+reg27; reg120=reg95+reg120; reg107=reg107*reg93; reg35=reg121+reg35;
    reg48=reg112*reg89; reg51=reg23+reg51; reg23=reg80*reg129; reg81=reg66+reg81; reg14=reg7+reg14;
    reg32=reg47*reg32; reg7=reg18*reg2; reg107=reg120+reg107; reg16=reg11+reg16; reg1=reg24*reg1;
    reg4=reg34+reg4; reg60=reg27+reg60; reg43=reg13+reg43; reg11=reg112*reg131; reg22=reg72+reg22;
    reg44=reg36+reg44; reg10=reg64+reg10; reg65=reg65*reg0; reg59=reg63+reg59; reg61=reg61*reg0;
    reg17=reg76+reg17; reg26=reg30-reg26; reg13=reg131*reg67; reg21=reg86+reg21; reg24=reg12*reg129;
    reg73=reg73*reg46; reg48=reg35+reg48; reg9=reg88+reg9; reg27=reg129*reg89; reg30=reg52*reg131;
    reg5=reg41+reg5; reg33=reg33*reg46; reg34=reg14*reg16; reg35=reg130*reg131; reg33=reg21+reg33;
    reg32=reg81+reg32; reg21=reg80*reg77; reg36=reg12*reg77; reg30=reg24-reg30; reg11=reg23-reg11;
    reg59=reg59*reg51; reg13=reg26+reg13; reg42=reg42*reg0; reg7=reg107+reg7; reg5=reg5*reg46;
    reg27=reg9+reg27; reg9=reg22*reg4; reg23=reg112*reg12; reg24=reg80*reg52; reg73=reg48+reg73;
    reg44=reg44*reg51; reg61=reg17+reg61; reg17=reg31*reg131; reg26=reg22*reg60; reg41=reg43*reg16;
    reg104=reg104*reg1; reg53=reg53*reg1; reg65=reg10+reg65; reg29=reg29*reg1; reg10=reg4*reg43;
    reg114=reg114/reg13; reg34=reg9-reg34; reg30=reg30/reg13; reg5=reg27+reg5; reg32=reg32*reg51;
    reg104=reg61+reg104; reg44=reg73+reg44; reg9=reg14*reg60; reg27=reg31*reg12; reg42=reg7+reg42;
    reg23=reg24-reg23; reg7=reg80*reg130; reg35=reg36-reg35; reg113=reg113/reg13; reg41=reg26-reg41;
    reg17=reg21-reg17; reg11=reg11/reg13; reg59=reg33+reg59; reg53=reg65+reg53; reg21=reg41*reg53;
    reg24=reg34*reg104; reg29=reg42+reg29; reg113=reg113*reg59; reg114=reg114*reg44; reg32=reg5+reg32;
    reg27=reg7-reg27; reg23=reg23/reg13; reg67=reg67/reg13; reg17=reg17/reg13; reg35=reg35/reg13;
    reg30=reg30*reg44; reg11=reg11*reg59; reg9=reg10-reg9; reg23=reg23*reg32; reg30=reg11-reg30;
    reg113=reg114-reg113; reg44=reg35*reg44; reg59=reg17*reg59; reg67=reg67*reg32; reg13=reg27/reg13;
    reg5=reg4*reg29; reg7=reg53*reg16; reg10=reg60*reg29; reg11=reg104*reg16; reg17=reg9*reg29;
    reg24=reg21-reg24; reg21=1-(*f.m).resolution; reg32=reg13*reg32; reg59=reg44-reg59; reg23=reg30-reg23;
    reg113=reg67+reg113; reg13=reg22*reg104; reg10=reg11-reg10; reg11=reg22*reg53; reg26=reg4*reg104;
    reg17=reg24+reg17; reg24=reg53*reg60; reg27=reg43*reg29; reg30=reg14*reg29; reg5=reg7-reg5;
    reg30=reg11-reg30; reg59=reg32+reg59; reg5=reg5/reg17; reg10=reg10/reg17; reg7=(*f.m).resolution*reg78;
    reg27=reg13-reg27; reg11=(*f.m).resolution*reg91; reg26=reg24-reg26; reg23=reg21*reg23; reg13=reg53*reg43;
    reg113=reg21*reg113; reg24=reg14*reg104; reg30=reg30/reg17; reg32=(*f.m).resolution*reg5; reg113=reg7+reg113;
    reg23=reg11+reg23; reg27=reg27/reg17; reg34=reg34/reg17; reg131=reg131*reg21; reg7=(*f.m).resolution*reg10;
    reg41=reg41/reg17; reg77=reg77*reg21; reg24=reg13-reg24; reg26=reg26/reg17; reg59=reg21*reg59;
    reg11=(*f.m).resolution*reg89; reg12=reg12*reg21; reg130=reg130*reg21; reg7=reg131+reg7; reg31=reg31*reg21;
    reg80=reg80*reg21; reg13=(*f.m).resolution*reg26; reg33=(*f.m).resolution*reg27; reg32=reg77-reg32; reg35=(*f.m).resolution*reg30;
    reg129=reg129*reg21; reg9=reg9/reg17; reg36=elem.pos(1)[0]-elem.pos(0)[0]; reg17=reg24/reg17; reg24=elem.pos(1)[1]-elem.pos(0)[1];
    reg42=elem.pos(2)[1]-elem.pos(0)[1]; reg44=elem.pos(2)[0]-elem.pos(0)[0]; reg48=(*f.m).resolution*reg34; reg61=(*f.m).resolution*reg41; reg59=reg11+reg59;
    reg113=reg113*(*f.m).deltaT; reg23=reg23*(*f.m).deltaT; reg11=(*f.m).resolution*reg17; reg80=reg61+reg80; reg48=reg31-reg48;
    reg33=reg12-reg33; reg12=reg32*reg23; reg35=reg130+reg35; reg52=reg52*reg21; reg31=(*f.m).resolution*reg9;
    reg59=reg59*(*f.m).deltaT; reg61=reg42*reg36; reg63=reg44*reg24; reg13=reg129+reg13; reg21=reg112*reg21;
    reg64=reg7*reg113; reg21=reg31+reg21; reg31=reg64+reg12; reg65=reg13*reg59; reg66=reg80*reg113;
    reg67=reg48*reg23; reg63=reg61-reg63; reg11=reg52-reg11; reg52=reg33*reg113; reg61=reg35*reg23;
    reg68=reg11*reg59; reg69=reg66+reg67; reg72=reg52+reg61; reg24=reg24/reg63; reg44=reg44/reg63;
    reg73=reg31+reg65; reg36=reg36/reg63; reg42=reg42/reg63; reg76=reg21*reg59; reg77=0.5*reg44;
    reg81=0.5*reg42; reg86=0.5*reg36; reg88=0.5*reg24; reg92=reg69+reg76; reg94=2*reg73;
    reg95=reg24-reg42; reg100=reg44-reg36; reg101=reg72+reg68; reg102=reg42*reg92; reg103=reg86*reg94;
    reg107=reg24*reg92; reg109=reg36*reg101; reg112=reg88*reg94; reg114=reg77*reg94; reg116=0.5*reg95;
    reg117=0.5*reg100; reg119=reg44*reg101; reg120=reg81*reg94; reg121=reg114-reg102; reg122=reg112-reg109;
    reg123=reg95*reg92; reg127=reg119-reg120; reg128=reg117*reg94; reg129=reg107-reg103; reg130=reg100*reg101;
    reg131=reg116*reg94; reg132=reg130+reg131; reg134=reg123+reg128; reg122=reg63*reg122; reg129=reg63*reg129;
    reg127=reg63*reg127; reg121=reg63*reg121; reg127=ponderation*reg127; reg129=ponderation*reg129; reg122=ponderation*reg122;
    reg135=reg63*reg132; reg121=ponderation*reg121; reg136=reg63*reg134; sollicitation[indices[1]+1]+=-reg127; sollicitation[indices[2]+0]+=-reg129;
    sollicitation[indices[1]+0]+=-reg121; reg121=ponderation*reg135; sollicitation[indices[0]+1]+=reg121; reg127=ponderation*reg136; sollicitation[indices[0]+0]+=reg127;
    sollicitation[indices[2]+1]+=-reg122;
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
    T reg0=2*(*f.m).shear_modulus_13; T reg1=2*(*f.m).shear_modulus_23; T reg2=2*(*f.m).shear_modulus_12; reg0=1.0/reg0; reg1=1.0/reg1;
    reg2=1.0/reg2; T reg3=reg0*reg1; T reg4=reg2*reg3; T reg5=1.0/(*f.m).elastic_modulus_3; T reg6=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg7=pow((*f.m).v2[0],2); T reg8=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg9=pow((*f.m).v2[1],2); T reg10=pow((*f.m).v1[0],2); T reg11=pow((*f.m).v1[1],2);
    T reg12=1.0/(*f.m).elastic_modulus_2; T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg14=reg4*reg5; T reg15=reg8*reg4; reg11=reg10+reg11;
    reg10=pow((*f.m).v1[2],2); T reg16=reg4*reg6; reg9=reg7+reg9; reg7=pow((*f.m).v2[2],2); reg10=reg11+reg10;
    reg7=reg9+reg7; reg9=reg12*reg14; reg11=reg13*reg14; T reg17=reg8*reg16; T reg18=reg8*reg15;
    reg18=reg9-reg18; reg17=reg11+reg17; reg9=reg13*reg15; reg10=pow(reg10,0.5); T reg19=reg12*reg16;
    T reg20=1.0/(*f.m).elastic_modulus_1; reg7=pow(reg7,0.5); T reg21=reg20*reg18; T reg22=(*f.m).v1[1]/reg10; T reg23=(*f.m).v1[2]/reg10;
    T reg24=reg13*reg17; T reg25=(*f.m).v2[1]/reg7; T reg26=(*f.m).v2[2]/reg7; T reg27=reg9+reg19; T reg28=reg2*reg1;
    reg10=(*f.m).v1[0]/reg10; T reg29=reg8*reg3; T reg30=reg16*reg6; T reg31=reg13*reg4; T reg32=reg3*reg6;
    reg3=reg3*reg5; reg14=reg20*reg14; reg4=reg12*reg4; T reg33=reg15*reg6; reg24=reg21-reg24;
    reg21=reg23*reg25; T reg34=reg22*reg26; T reg35=reg27*reg6; reg7=(*f.m).v2[0]/reg7; T reg36=reg13*reg3;
    reg3=reg12*reg3; T reg37=reg28*reg6; reg15=reg20*reg15; T reg38=reg2*reg0; T reg39=reg4*reg6;
    reg30=reg14-reg30; reg33=reg11+reg33; reg16=reg13*reg16; reg11=reg8*reg28; reg35=reg24-reg35;
    reg14=reg31*reg6; reg24=2*reg7; reg28=reg28*reg5; T reg40=reg34-reg21; T reg41=reg23*reg7;
    T reg42=2*reg10; T reg43=reg8*reg32; T reg44=reg10*reg26; T reg45=reg8*reg29; reg39=reg9+reg39;
    reg16=reg15+reg16; reg33=reg33/reg35; T reg46=pow(reg7,2); T reg47=pow(reg25,2); reg14=reg15+reg14;
    reg4=reg20*reg4; reg15=reg24*reg25; T reg48=reg10*reg25; T reg49=2*reg40; T reg50=reg22*reg7;
    T reg51=reg41-reg44; reg18=reg18/reg35; T reg52=pow(reg10,2); T reg53=pow(reg22,2); T reg54=reg38*reg6;
    reg45=reg3-reg45; reg32=reg12*reg32; reg43=reg36+reg43; reg3=reg8*reg38; reg30=reg30/reg35;
    reg29=reg13*reg29; reg36=reg12*reg28; reg28=reg13*reg28; T reg55=reg8*reg11; T reg56=reg22*reg42;
    T reg57=reg8*reg37; reg31=reg13*reg31; reg17=reg17/reg35; reg38=reg38*reg5; T reg58=reg46*reg17;
    T reg59=reg18*reg52; T reg60=reg33*reg53; T reg61=reg47*reg30; T reg62=reg56*reg33; T reg63=reg15*reg30;
    reg45=reg20*reg45; reg37=reg12*reg37; reg11=reg13*reg11; reg43=reg13*reg43; T reg64=reg29+reg32;
    reg31=reg4-reg31; reg55=reg36-reg55; reg57=reg28+reg57; reg4=reg12*reg38; reg16=reg16/reg35;
    reg38=reg13*reg38; reg28=reg8*reg3; reg36=reg8*reg54; T reg65=reg18*reg53; reg14=reg14/reg35;
    T reg66=reg47*reg17; T reg67=reg56*reg18; T reg68=reg15*reg17; T reg69=reg33*reg52; T reg70=reg46*reg30;
    T reg71=pow(reg23,2); T reg72=reg48-reg50; reg27=reg27/reg35; reg39=reg39/reg35; T reg73=reg49*reg51;
    T reg74=pow(reg51,2); T reg75=pow(reg40,2); T reg76=pow(reg26,2); T reg77=reg33*reg71; T reg78=reg76*reg30;
    T reg79=reg11+reg37; reg63=reg62+reg63; reg57=reg13*reg57; reg62=reg73*reg16; reg55=reg20*reg55;
    T reg80=reg39*reg52; T reg81=reg46*reg14; reg64=reg64*reg6; T reg82=reg39*reg53; reg43=reg45-reg43;
    reg45=reg47*reg14; reg28=reg4-reg28; reg4=reg74*reg16; reg61=reg60+reg61; reg36=reg38+reg36;
    reg38=reg75*reg16; reg70=reg69+reg70; reg66=reg65+reg66; reg60=reg74*reg27; reg65=reg18*reg71;
    reg69=reg76*reg17; reg68=reg67+reg68; reg67=reg73*reg27; T reg83=pow(reg72,2); T reg84=reg15*reg14;
    T reg85=reg56*reg39; reg58=reg59+reg58; reg31=reg31/reg35; reg54=reg12*reg54; reg59=reg75*reg27;
    reg3=reg13*reg3; reg78=reg77+reg78; reg84=reg85+reg84; reg28=reg20*reg28; reg77=reg73*reg31;
    reg4=reg61+reg4; reg59=reg58+reg59; reg36=reg13*reg36; reg58=reg42*reg23; reg61=reg3+reg54;
    reg38=reg70+reg38; reg70=2*reg25; reg85=reg24*reg26; reg67=reg68+reg67; reg60=reg66+reg60;
    reg69=reg65+reg69; reg65=reg27*reg83; reg66=2*reg22; reg68=reg16*reg83; reg79=reg79*reg6;
    T reg86=reg7*reg25; T reg87=reg10*reg22; reg45=reg82+reg45; reg82=reg74*reg31; reg57=reg55-reg57;
    reg62=reg63+reg62; reg55=reg39*reg71; reg63=reg75*reg31; reg81=reg80+reg81; reg80=reg76*reg14;
    reg64=reg43-reg64; reg43=reg59*reg53; T reg88=reg86*reg4; T reg89=reg87*reg60; T reg90=reg70*reg26;
    T reg91=reg67*reg52; reg48=reg50+reg48; reg50=reg22*reg25; T reg92=reg10*reg7; T reg93=reg46*reg62;
    T reg94=2*reg51; reg49=reg49*reg72; reg82=reg45+reg82; reg63=reg81+reg63; reg80=reg55+reg80;
    reg45=reg47*reg62; reg55=reg85*reg30; reg81=reg58*reg33; T reg95=reg31*reg83; reg68=reg78+reg68;
    reg77=reg84+reg77; reg78=reg67*reg53; reg84=reg47*reg4; T reg96=reg60*reg53; T reg97=reg59*reg52;
    T reg98=reg46*reg38; T reg99=reg47*reg38; reg60=reg60*reg52; T reg100=reg66*reg23; reg4=reg46*reg4;
    T reg101=reg51*reg40; reg64=reg64/reg35; reg79=reg57-reg79; reg36=reg28-reg36; reg61=reg61*reg6;
    reg65=reg69+reg65; reg28=reg58*reg18; reg57=reg85*reg17; reg62=reg86*reg62; reg67=reg87*reg67;
    reg69=reg10*reg51; T reg102=reg22*reg40; T reg103=reg49*reg27; reg57=reg28+reg57; reg88=reg89+reg88;
    reg28=reg75*reg63; reg98=reg97+reg98; reg89=reg101*reg77; reg97=reg42*reg7; reg61=reg36-reg61;
    reg36=reg101*reg82; reg84=reg96+reg84; reg96=reg74*reg82; T reg104=reg92*reg64; T reg105=reg40*reg25;
    T reg106=reg51*reg7; T reg107=reg75*reg77; T reg108=reg50*reg64; T reg109=reg65*reg53; reg30=reg90*reg30;
    reg79=reg79/reg35; T reg110=reg47*reg68; reg33=reg100*reg33; reg55=reg81+reg55; reg81=reg49*reg16;
    T reg111=reg48*reg64; reg95=reg80+reg95; reg80=reg46*reg68; T reg112=reg65*reg52; reg38=reg86*reg38;
    reg59=reg87*reg59; T reg113=2*reg23; reg82=reg75*reg82; reg4=reg60+reg4; reg60=reg22*reg51;
    T reg114=reg10*reg40; reg69=reg102+reg69; reg62=reg67+reg62; reg67=reg58*reg39; reg45=reg78+reg45;
    reg77=reg74*reg77; reg93=reg91+reg93; reg18=reg100*reg18; reg17=reg90*reg17; reg78=reg66*reg25;
    reg94=reg94*reg72; reg91=reg85*reg14; reg102=reg74*reg63; T reg115=reg23*reg26; reg99=reg43+reg99;
    reg43=reg49*reg31; T reg116=reg115*reg64; T reg117=reg74*reg95; reg39=reg100*reg39; reg91=reg67+reg91;
    reg14=reg90*reg14; reg80=reg112+reg80; reg67=reg75*reg95; reg112=reg97*reg108; reg82=reg4+reg82;
    reg102=reg99+reg102; reg4=reg78*reg104; reg99=reg97*reg104; reg28=reg98+reg28; reg98=reg97*reg111;
    reg89=reg62+reg89; reg62=reg48*reg111; reg107=reg93+reg107; reg96=reg84+reg96; reg84=reg78*reg108;
    reg93=reg69*reg79; T reg118=reg60*reg79; T reg119=reg114*reg79; reg110=reg109+reg110; reg30=reg33+reg30;
    reg81=reg55+reg81; reg63=reg101*reg63; reg38=reg59+reg38; reg33=reg46*reg13; reg106=reg105+reg106;
    reg55=reg51*reg25; reg35=reg61/reg35; reg59=reg40*reg7; reg61=reg47*reg12; reg105=reg97*reg13;
    reg111=reg78*reg111; reg103=reg57+reg103; reg77=reg45+reg77; reg17=reg18+reg17; reg27=reg94*reg27;
    reg18=reg78*reg12; reg113=reg113*reg26; reg68=reg86*reg68; reg65=reg87*reg65; reg45=reg46*reg20;
    reg66=reg66*reg51; reg57=reg47*reg13; reg42=reg42*reg40; reg108=reg48*reg108; reg16=reg94*reg16;
    reg109=reg78*reg13; reg36=reg88+reg36; reg88=reg72*reg23; T reg120=reg97*reg20; T reg121=reg47*reg8;
    reg68=reg65+reg68; reg95=reg101*reg95; reg108=reg36+reg108; reg63=reg38+reg63; reg36=reg46*reg6;
    reg104=reg48*reg104; reg111=reg77+reg111; reg38=reg66*reg93; reg65=reg103*reg53; reg77=reg47*reg81;
    reg98=reg107+reg98; reg105=reg18-reg105; reg18=reg46*reg81; reg107=reg103*reg52; T reg122=reg42*reg93;
    T reg123=reg47*(*f.m).alpha_2; T reg124=(*f.m).alpha_1*reg53; T reg125=reg46*(*f.m).alpha_2; reg93=reg69*reg93; T reg126=reg76*reg6;
    reg57=reg45-reg57; reg45=reg78*reg116; reg117=reg110+reg117; reg110=reg113*reg6; reg109=reg120-reg109;
    reg120=(*f.m).alpha_1*reg52; T reg127=reg78*reg8; T reg128=reg66*reg118; reg84=reg96+reg84; reg96=reg76*reg8;
    reg62=reg89+reg62; reg89=reg97*reg6; T reg129=reg66*reg119; reg4=reg102+reg4; reg33=reg61-reg33;
    reg61=reg113*reg8; reg102=reg69*reg118; reg99=reg28+reg99; reg24=reg24*reg40; reg28=reg10*reg72;
    reg70=reg70*reg51; T reg130=reg40*reg23; T reg131=reg42*reg119; reg43=reg91+reg43; reg14=reg39+reg14;
    reg31=reg94*reg31; reg118=reg42*reg118; reg16=reg30+reg16; reg41=reg44+reg41; reg112=reg82+reg112;
    reg30=reg106*reg35; reg39=reg55*reg35; reg67=reg80+reg67; reg44=reg97*reg116; reg80=reg59*reg35;
    reg82=reg88*reg79; reg27=reg17+reg27; reg17=reg72*reg26; reg91=reg40*reg26; reg93=reg62+reg93;
    reg62=reg106*reg30; T reg132=reg66*reg82; reg45=reg117+reg45; reg31=reg14+reg31; reg14=reg41*reg64;
    reg103=reg87*reg103; reg81=reg86*reg81; reg117=reg70*reg39; reg119=reg69*reg119; reg104=reg63+reg104;
    reg102=reg108+reg102; reg127=reg89+reg127; reg113=reg113*reg5; reg63=reg106*reg39; reg121=reg36+reg121;
    reg36=reg76*reg5; reg34=reg21+reg34; reg95=reg68+reg95; reg116=reg48*reg116; reg21=reg47*reg16;
    reg68=reg27*reg53; reg61=reg105-reg61; reg89=reg74*reg43; reg77=reg65+reg77; reg28=reg130+reg28;
    reg65=reg51*reg23; reg105=reg70*reg30; reg38=reg111+reg38; reg96=reg33-reg96; reg33=reg22*reg72;
    reg108=reg75*reg43; reg18=reg107+reg18; reg118=reg112+reg118; reg107=reg13*reg53; reg123=reg124+reg123;
    reg111=(*f.m).alpha_3*reg74; reg112=reg27*reg52; reg124=reg46*reg16; reg30=reg24*reg30; reg130=(*f.m).alpha_1*reg71;
    reg122=reg98+reg122; reg98=reg76*(*f.m).alpha_2; reg126=reg57-reg126; reg39=reg24*reg39; reg44=reg67+reg44;
    reg57=reg42*reg82; reg128=reg84+reg128; reg67=reg17*reg35; reg12=reg12*reg53; reg120=reg125+reg120;
    reg84=(*f.m).alpha_3*reg75; reg125=reg70*reg80; reg129=reg4+reg129; reg20=reg20*reg52; reg13=reg13*reg52;
    reg131=reg99+reg131; reg4=reg24*reg80; reg110=reg109-reg110; reg99=reg72*reg7; reg109=reg53*reg61;
    reg84=reg120+reg84; reg63=reg102+reg63; reg102=reg92*reg110; reg120=reg50*reg61; T reg133=reg52*reg110;
    T reg134=reg50*reg96; T reg135=reg53*reg96; reg81=reg103+reg81; reg43=reg101*reg43; reg103=reg86*(*f.m).alpha_2;
    reg27=reg87*reg27; reg62=reg93+reg62; reg16=reg86*reg16; reg93=(*f.m).alpha_1*reg87; T reg136=reg47*reg61;
    T reg137=(*f.m).alpha_3*reg83; reg98=reg130+reg98; reg130=reg46*reg110; reg111=reg123+reg111; reg123=reg47*reg96;
    T reg138=reg46*reg126; T reg139=reg92*reg126; reg82=reg69*reg82; reg10=reg10*reg23; reg7=reg7*reg26;
    reg116=reg95+reg116; reg89=reg77+reg89; reg105=reg38+reg105; reg33=reg65+reg33; reg38=reg70*reg67;
    reg132=reg45+reg132; reg64=reg34*reg64; reg117=reg128+reg117; reg45=reg28*reg79; reg13=reg12-reg13;
    reg12=reg8*reg71; reg125=reg129+reg125; reg4=reg131+reg4; reg99=reg91+reg99; reg65=reg51*reg26;
    reg77=reg72*reg25; reg91=reg75*reg31; reg124=reg112+reg124; reg107=reg20-reg107; reg20=reg97*reg14;
    reg108=reg18+reg108; reg30=reg122+reg30; reg39=reg118+reg39; reg57=reg44+reg57; reg18=reg24*reg67;
    reg44=reg71*reg6; reg121=reg36-reg121; reg127=reg113-reg127; reg36=reg74*reg31; reg21=reg68+reg21;
    reg8=reg8*reg53; reg68=reg52*reg126; reg119=reg104+reg119; reg80=reg106*reg80; reg6=reg52*reg6;
    reg95=reg78*reg14; reg104=reg71*reg121; reg109=reg133+reg109; reg135=reg68+reg135; reg40=reg72*reg40;
    reg68=reg71*reg127; reg112=reg117*reg111; reg44=reg107-reg44; reg107=reg48*reg2; reg86=reg86*reg2;
    reg77=reg65+reg77; reg137=reg98+reg137; reg65=reg125*reg84; reg123=reg138+reg123; reg98=reg76*reg121;
    reg113=reg117*reg62; reg118=reg39*reg62; reg12=reg13-reg12; reg13=reg105*reg63; reg122=reg30*reg63;
    reg5=reg71*reg5; reg136=reg130+reg136; reg128=reg76*reg127; reg129=reg7*(*f.m).alpha_2; reg130=(*f.m).alpha_1*reg10;
    reg131=(*f.m).alpha_3*reg101; reg103=reg93+reg103; reg8=reg6+reg8; reg6=reg4*reg84; reg93=reg39*reg111;
    reg31=reg101*reg31; reg16=reg27+reg16; reg14=reg48*reg14; reg43=reg81+reg43; reg67=reg106*reg67;
    reg82=reg116+reg82; reg134=reg139+reg134; reg27=reg115*reg121; reg80=reg119+reg80; reg78=reg78*reg64;
    reg36=reg21+reg36; reg21=reg66*reg45; reg95=reg89+reg95; reg38=reg132+reg38; reg79=reg33*reg79;
    reg81=reg99*reg35; reg97=reg97*reg64; reg91=reg124+reg91; reg89=reg42*reg45; reg20=reg108+reg20;
    reg18=reg57+reg18; reg23=reg22*reg23; reg26=reg25*reg26; reg22=reg115*reg127; reg120=reg102+reg120;
    reg51=reg72*reg51; reg25=reg24*reg81; reg57=reg63*reg111; reg72=reg80*reg84; reg101=reg48*reg86;
    reg68=reg109+reg68; reg102=reg56*reg107; reg108=reg15*reg86; reg98=reg123+reg98; reg97=reg91+reg97;
    reg42=reg42*reg79; reg91=reg46*reg44; reg66=reg66*reg79; reg35=reg77*reg35; reg38=reg38*reg137;
    reg67=reg82+reg67; reg22=reg120+reg22; reg112=reg65+reg112; reg13=reg113-reg13; reg27=reg134+reg27;
    reg78=reg36+reg78; reg122=reg118-reg122; reg36=reg39*reg105; reg65=reg30*reg117; reg18=reg18*reg137;
    reg93=reg6+reg93; reg6=reg47*reg12; reg8=reg5-reg8; reg5=reg70*reg81; reg21=reg95+reg21;
    reg82=reg48*reg107; reg95=reg52*reg44; reg109=reg41*reg0; reg113=reg53*reg12; reg64=reg48*reg64;
    reg31=reg16+reg31; reg7=reg7*reg0; reg131=reg103+reg131; reg45=reg69*reg45; reg14=reg43+reg14;
    reg129=reg130+reg129; reg40=(*f.m).alpha_3*reg40; reg16=(*f.m).alpha_1*reg23; reg43=reg26*(*f.m).alpha_2; reg104=reg135+reg104;
    reg103=reg56*reg86; reg128=reg136+reg128; reg116=reg15*reg107; reg89=reg20+reg89; reg76=reg76*reg8;
    reg20=reg34*reg1; reg26=reg26*reg1; reg6=reg91+reg6; reg91=reg92*reg44; reg5=reg21+reg5;
    reg21=reg50*reg12; reg57=reg72+reg57; reg67=reg67*reg137; reg18=reg93+reg18; reg101=reg27+reg101;
    reg113=reg95+reg113; reg27=reg30*reg131; reg71=reg71*reg8; reg66=reg78+reg66; reg70=reg70*reg35;
    reg72=reg41*reg7; reg78=reg58*reg7; reg103=reg104+reg103; reg102=reg68+reg102; reg38=reg112+reg38;
    reg2=reg87*reg2; reg68=reg105*reg131; reg87=reg58*reg109; reg93=reg4*reg13; reg95=reg125*reg122;
    reg51=(*f.m).alpha_3*reg51; reg43=reg16+reg43; reg24=reg24*reg35; reg65=reg36-reg65; reg42=reg97+reg42;
    reg40=reg129+reg40; reg108=reg98+reg108; reg16=reg85*reg7; reg25=reg89+reg25; reg36=reg41*reg109;
    reg82=reg22+reg82; reg79=reg69*reg79; reg64=reg31+reg64; reg22=reg85*reg109; reg116=reg128+reg116;
    reg45=reg14+reg45; reg81=reg106*reg81; reg22=reg116+reg22; reg70=reg66+reg70; reg14=reg90*reg20;
    reg16=reg108+reg16; reg24=reg42+reg24; reg68=reg38+reg68; reg5=reg5*reg40; reg78=reg103+reg78;
    reg31=reg100*reg26; reg38=reg62*reg131; reg67=reg57+reg67; reg72=reg101+reg72; reg42=reg100*reg20;
    reg87=reg102+reg87; reg57=reg34*reg20; reg36=reg82+reg36; reg66=reg90*reg26; reg81=reg45+reg81;
    reg0=reg10*reg0; reg71=reg113+reg71; reg10=reg56*reg2; reg45=reg4*reg62; reg82=reg30*reg80;
    reg89=reg15*reg2; reg79=reg64+reg79; reg64=reg105*reg80; reg76=reg6+reg76; reg6=reg125*reg62;
    reg35=reg106*reg35; reg21=reg91+reg21; reg91=reg80*reg65; reg95=reg93-reg95; reg27=reg18+reg27;
    reg25=reg25*reg40; reg18=reg34*reg26; reg115=reg115*reg8; reg51=reg43+reg51; reg1=reg23*reg1;
    reg35=reg79+reg35; reg10=reg71+reg10; reg58=reg58*reg0; reg31=reg78+reg31; reg18=reg72+reg18;
    reg23=reg4*reg63; reg82=reg45-reg82; reg43=reg39*reg80; reg45=reg48*reg2; reg71=reg117*reg80;
    reg72=reg4*reg105; reg78=reg30*reg125; reg64=reg6-reg64; reg6=reg125*reg63; reg81=reg81*reg40;
    reg38=reg67+reg38; reg66=reg16+reg66; reg42=reg87+reg42; reg70=reg70*reg51; reg5=reg68+reg5;
    reg85=reg85*reg0; reg89=reg76+reg89; reg14=reg22+reg14; reg24=reg24*reg51; reg25=reg27+reg25;
    reg91=reg95+reg91; reg57=reg36+reg57; reg115=reg21+reg115; reg45=reg115+reg45; reg16=elem.pos(1)[1]-elem.pos(0)[1];
    reg21=elem.pos(2)[1]-elem.pos(0)[1]; reg41=reg41*reg0; reg122=reg122/reg91; reg82=reg82/reg91; reg71=reg6-reg71;
    reg64=reg64/reg91; reg43=reg23-reg43; reg13=reg13/reg91; reg6=reg4*reg117; reg78=reg72-reg78;
    reg22=reg39*reg125; reg85=reg89+reg85; reg24=reg25+reg24; reg23=reg42*reg18; reg70=reg5+reg70;
    reg81=reg38+reg81; reg35=reg35*reg51; reg5=elem.pos(2)[0]-elem.pos(0)[0]; reg90=reg90*reg1; reg100=reg100*reg1;
    reg58=reg10+reg58; reg10=reg57*reg31; reg25=reg14*reg18; reg27=reg57*reg66; reg36=elem.pos(1)[0]-elem.pos(0)[0];
    reg38=reg42*reg66; reg25=reg27-reg25; reg90=reg85+reg90; reg41=reg45+reg41; reg71=reg71/reg91;
    reg34=reg34*reg1; reg27=reg31*reg14; reg43=reg43/reg91; reg23=reg10-reg23; reg65=reg65/reg91;
    reg78=reg78/reg91; reg22=reg6-reg22; reg35=reg81+reg35; reg6=reg5*reg16; reg100=reg58+reg100;
    reg10=reg21*reg36; reg13=reg13*reg24; reg122=reg122*reg70; reg64=reg64*reg24; reg82=reg82*reg70;
    reg34=reg41+reg34; reg6=reg10-reg6; reg38=reg27-reg38; reg64=reg82-reg64; reg78=reg78*reg35;
    reg24=reg71*reg24; reg122=reg13-reg122; reg10=reg23*reg90; reg13=reg25*reg100; reg70=reg43*reg70;
    reg65=reg65*reg35; reg91=reg22/reg91; reg35=reg91*reg35; reg22=reg100*reg18; reg70=reg24-reg70;
    reg24=reg66*reg34; reg78=reg64-reg78; reg27=reg31*reg34; reg16=reg16/reg6; reg41=reg90*reg18;
    reg5=reg5/reg6; reg43=reg38*reg34; reg10=reg13-reg10; reg13=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg122=reg65+reg122;
    reg45=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg21=reg21/reg6; reg58=1-(*f.m).resolution; reg36=reg36/reg6; reg64=vectors[0][indices[1]+0]-vectors[0][indices[0]+0];
    reg65=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; reg122=reg58*reg122; reg67=reg100*reg66; reg68=reg31*reg90; reg71=(*f.m).resolution*reg84;
    reg72=(*f.m).resolution*reg111; reg70=reg35+reg70; reg78=reg58*reg78; reg35=reg14*reg34; reg27=reg22-reg27;
    reg22=reg57*reg90; reg76=reg5*reg64; reg79=reg42*reg34; reg24=reg41-reg24; reg41=reg36*reg65;
    reg81=reg57*reg100; reg82=reg21*reg45; reg43=reg10+reg43; reg10=reg16*reg13; reg85=(*f.m).resolution*reg131;
    reg10=reg82-reg10; reg35=reg22-reg35; reg24=reg24/reg43; reg79=reg81-reg79; reg27=reg27/reg43;
    reg68=reg67-reg68; reg22=reg100*reg14; reg67=reg42*reg90; reg76=reg41-reg76; reg70=reg58*reg70;
    reg122=reg71+reg122; reg78=reg72+reg78; reg45=reg5*reg45; reg13=reg36*reg13; reg64=reg21*reg64;
    reg65=reg16*reg65; reg76=reg10+reg76; reg79=reg79/reg43; reg78=reg78*(*f.m).deltaT; reg67=reg22-reg67;
    reg122=reg122*(*f.m).deltaT; reg68=reg68/reg43; reg70=reg85+reg70; reg23=reg23/reg43; reg65=reg64-reg65;
    reg35=reg35/reg43; reg45=reg13-reg45; reg25=reg25/reg43; reg10=(*f.m).resolution*reg27; reg13=(*f.m).resolution*reg24;
    reg80=reg80*reg58; reg63=reg63*reg58; reg65=reg65-reg122; reg76=0.5*reg76; reg38=reg38/reg43;
    reg22=(*f.m).resolution*reg68; reg13=reg80+reg13; reg43=reg67/reg43; reg70=reg70*(*f.m).deltaT; reg10=reg63-reg10;
    reg41=(*f.m).resolution*reg79; reg63=(*f.m).resolution*reg23; reg64=(*f.m).resolution*reg25; reg125=reg125*reg58; reg67=(*f.m).resolution*reg35;
    reg62=reg62*reg58; reg45=reg45-reg78; reg117=reg117*reg58; reg39=reg39*reg58; reg4=reg4*reg58;
    reg41=reg117+reg41; reg67=reg125-reg67; reg63=reg39-reg63; reg30=reg30*reg58; reg4=reg64+reg4;
    reg39=reg13*reg65; reg64=reg10*reg45; reg76=reg76-reg70; reg58=reg105*reg58; reg71=(*f.m).resolution*reg43;
    reg72=(*f.m).resolution*reg38; reg22=reg62+reg22; reg62=reg4*reg65; reg80=reg63*reg45; reg65=reg67*reg65;
    reg45=reg41*reg45; reg71=reg58-reg71; reg30=reg72+reg30; reg58=reg22*reg76; reg64=reg39+reg64;
    reg39=reg16-reg21; reg72=reg5-reg36; reg45=reg65+reg45; reg65=reg71*reg76; reg58=reg64+reg58;
    reg80=reg62+reg80; reg76=reg30*reg76; reg76=reg80+reg76; reg62=0.5*reg39; reg64=0.5*reg72;
    reg80=0.5*reg5; reg81=0.5*reg21; reg58=2*reg58; reg82=0.5*reg36; reg85=0.5*reg16;
    reg65=reg45+reg65; reg45=reg85*reg58; reg87=reg80*reg58; reg89=reg36*reg65; reg91=reg81*reg58;
    reg93=reg21*reg76; reg95=reg64*reg58; reg97=reg62*reg58; reg98=reg72*reg65; reg101=reg82*reg58;
    reg102=reg39*reg76; reg103=reg16*reg76; reg104=reg5*reg65; reg95=reg102+reg95; reg93=reg93-reg87;
    reg97=reg98+reg97; reg91=reg91-reg104; reg101=reg101-reg103; reg89=reg89-reg45; reg93=reg6*reg93;
    reg91=reg6*reg91; reg101=reg6*reg101; reg89=reg6*reg89; reg97=reg6*reg97; reg95=reg6*reg95;
    sollicitation[indices[2]+0]+=ponderation*reg101; sollicitation[indices[1]+1]+=ponderation*reg91; sollicitation[indices[2]+1]+=ponderation*reg89; sollicitation[indices[0]+1]+=ponderation*reg97; sollicitation[indices[1]+0]+=ponderation*reg93;
    sollicitation[indices[0]+0]+=ponderation*reg95;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; T reg2=2*(*f.m).shear_modulus_12; reg1=1.0/reg1; reg0=1.0/reg0;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=reg2*reg3; T reg5=pow((*f.m).v1[0],2); T reg6=pow((*f.m).v1[1],2);
    T reg7=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg8=1.0/(*f.m).elastic_modulus_3; T reg9=pow((*f.m).v2[0],2); T reg10=pow((*f.m).v2[1],2); T reg11=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg12=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg13=1.0/(*f.m).elastic_modulus_2; T reg14=reg8*reg4; reg6=reg5+reg6; reg5=pow((*f.m).v1[2],2);
    reg10=reg9+reg10; reg9=pow((*f.m).v2[2],2); T reg15=reg7*reg4; T reg16=reg11*reg4; T reg17=reg7*reg16;
    T reg18=reg14*reg13; T reg19=reg7*reg15; reg5=reg6+reg5; reg6=reg14*reg12; reg9=reg10+reg9;
    reg10=reg12*reg15; T reg20=reg13*reg16; reg17=reg6+reg17; reg5=pow(reg5,0.5); reg9=pow(reg9,0.5);
    reg19=reg18-reg19; reg18=1.0/(*f.m).elastic_modulus_1; T reg21=(*f.m).v2[2]/reg9; T reg22=(*f.m).v2[1]/reg9; T reg23=(*f.m).v1[2]/reg5;
    T reg24=(*f.m).v1[1]/reg5; T reg25=reg18*reg19; T reg26=reg10+reg20; T reg27=reg12*reg17; reg14=reg14*reg18;
    T reg28=reg11*reg16; T reg29=reg12*reg4; T reg30=reg8*reg3; reg27=reg25-reg27; reg25=reg2*reg0;
    reg5=(*f.m).v1[0]/reg5; T reg31=reg7*reg3; T reg32=reg11*reg15; reg3=reg11*reg3; reg4=reg13*reg4;
    reg9=(*f.m).v2[0]/reg9; T reg33=reg11*reg26; T reg34=reg23*reg22; T reg35=reg24*reg21; reg16=reg12*reg16;
    T reg36=reg7*reg25; T reg37=reg23*reg9; T reg38=2*reg9; reg33=reg27-reg33; reg27=reg11*reg29;
    T reg39=reg13*reg30; reg15=reg18*reg15; T reg40=reg5*reg21; reg28=reg14-reg28; reg14=reg2*reg1;
    T reg41=reg35-reg34; T reg42=reg11*reg25; T reg43=2*reg24; T reg44=reg7*reg3; T reg45=2*reg5;
    reg25=reg8*reg25; T reg46=reg7*reg31; reg30=reg12*reg30; T reg47=reg11*reg4; reg32=reg6+reg32;
    reg6=pow(reg24,2); T reg48=reg8*reg14; reg19=reg19/reg33; reg32=reg32/reg33; reg28=reg28/reg33;
    reg47=reg10+reg47; T reg49=reg9*reg45; T reg50=reg37-reg40; T reg51=reg22*reg38; reg27=reg15+reg27;
    reg17=reg17/reg33; T reg52=reg24*reg45; T reg53=pow(reg5,2); reg44=reg30+reg44; reg30=reg24*reg9;
    T reg54=reg5*reg22; reg46=reg39-reg46; reg39=reg13*reg25; T reg55=pow(reg9,2); T reg56=reg11*reg14;
    T reg57=pow(reg22,2); reg25=reg12*reg25; reg3=reg13*reg3; reg14=reg7*reg14; reg31=reg12*reg31;
    T reg58=reg22*reg43; T reg59=2*reg41; T reg60=2*reg23; reg16=reg15+reg16; reg4=reg18*reg4;
    reg15=reg7*reg42; reg29=reg12*reg29; T reg61=reg7*reg36; T reg62=reg57*reg12; T reg63=pow(reg41,2);
    T reg64=pow(reg21,2); T reg65=reg55*reg18; T reg66=reg49*reg18; T reg67=pow(reg23,2); T reg68=pow(reg50,2);
    T reg69=reg50*reg59; T reg70=reg54-reg30; T reg71=reg57*reg13; reg44=reg12*reg44; T reg72=reg31+reg3;
    T reg73=reg49*reg12; reg61=reg39-reg61; reg15=reg25+reg15; reg25=reg58*reg13; reg39=reg13*reg48;
    reg48=reg12*reg48; T reg74=reg7*reg14; T reg75=reg7*reg56; T reg76=reg6*reg19; T reg77=reg57*reg17;
    T reg78=reg21*reg60; T reg79=reg52*reg19; T reg80=reg6*reg32; T reg81=reg51*reg17; T reg82=reg53*reg32;
    T reg83=reg55*reg28; T reg84=reg51*reg28; T reg85=reg52*reg32; reg16=reg16/reg33; T reg86=reg58*reg12;
    reg27=reg27/reg33; reg29=reg4-reg29; reg36=reg12*reg36; reg4=reg53*reg19; reg26=reg26/reg33;
    T reg87=reg55*reg17; reg42=reg13*reg42; T reg88=reg55*reg12; reg47=reg47/reg33; T reg89=reg57*reg28;
    reg46=reg18*reg46; T reg90=reg58*reg7; T reg91=reg78*reg11; T reg92=reg49*reg11; T reg93=reg57*reg7;
    T reg94=reg55*reg11; reg86=reg66-reg86; reg73=reg25-reg73; reg25=reg63*reg26; reg87=reg4+reg87;
    reg4=reg78*reg7; reg88=reg71-reg88; reg66=reg64*reg7; reg71=reg63*reg16; reg83=reg82+reg83;
    reg82=reg69*reg26; reg81=reg79+reg81; reg79=reg64*reg17; T reg95=reg67*reg19; T reg96=reg68*reg26;
    reg77=reg76+reg77; reg75=reg48+reg75; reg74=reg39-reg74; reg39=reg36+reg42; reg15=reg12*reg15;
    reg61=reg18*reg61; reg72=reg11*reg72; reg44=reg46-reg44; reg56=reg13*reg56; reg14=reg12*reg14;
    reg29=reg29/reg33; reg89=reg80+reg89; reg46=reg68*reg16; reg48=reg67*reg32; reg76=reg64*reg28;
    reg84=reg85+reg84; reg80=reg69*reg16; reg85=reg53*reg47; T reg97=reg55*reg27; T reg98=reg6*reg47;
    T reg99=reg57*reg27; T reg100=reg52*reg47; T reg101=reg51*reg27; T reg102=reg64*reg11; reg62=reg65-reg62;
    reg65=pow(reg70,2); reg102=reg62-reg102; reg91=reg86-reg91; reg62=reg53*reg12; reg86=reg6*reg13;
    reg75=reg12*reg75; reg66=reg88-reg66; reg82=reg81+reg82; reg25=reg87+reg25; reg4=reg73-reg4;
    reg73=reg24*reg22; reg81=reg5*reg9; reg87=reg14+reg56; reg88=2*reg22; T reg103=reg21*reg38;
    T reg104=reg23*reg45; reg80=reg84+reg80; reg84=reg64*reg27; T reg105=reg65*reg26; reg79=reg95+reg79;
    reg95=reg64*reg8; T reg106=reg67*reg47; reg93=reg94+reg93; reg94=reg68*reg29; reg99=reg98+reg99;
    reg96=reg77+reg96; reg90=reg92+reg90; reg77=reg78*reg8; reg15=reg61-reg15; reg46=reg89+reg46;
    reg39=reg11*reg39; reg76=reg48+reg76; reg48=reg65*reg16; reg71=reg83+reg71; reg72=reg44-reg72;
    reg44=reg53*reg18; reg61=reg5*reg24; reg83=reg9*reg22; reg74=reg18*reg74; reg89=reg69*reg29;
    reg101=reg100+reg101; reg97=reg85+reg97; reg85=reg63*reg29; reg92=reg6*reg12; reg92=reg44-reg92;
    reg44=reg6*reg82; reg39=reg15-reg39; reg15=reg55*reg71; reg93=reg95-reg93; reg95=reg67*reg11;
    reg98=reg57*reg80; reg100=reg67*reg7; T reg107=reg61*reg96; T reg108=reg83*reg46; T reg109=reg6*reg96;
    T reg110=reg53*reg25; T reg111=reg57*reg46; reg72=reg72/reg33; reg62=reg86-reg62; reg86=reg6*reg7;
    reg85=reg97+reg85; reg75=reg74-reg75; reg74=reg104*reg32; reg97=reg53*reg11; T reg112=reg103*reg28;
    T reg113=reg61*reg82; reg87=reg11*reg87; T reg114=reg83*reg80; T reg115=reg104*reg19; T reg116=reg73*reg4;
    T reg117=reg24*reg41; T reg118=reg5*reg50; T reg119=reg103*reg17; T reg120=reg81*reg91; T reg121=reg65*reg29;
    T reg122=reg55*reg102; T reg123=reg57*reg66; T reg124=reg73*reg66; T reg125=reg55*reg91; T reg126=reg57*reg4;
    T reg127=reg55*reg80; T reg128=reg53*reg82; T reg129=reg53*reg96; T reg130=reg50*reg41; reg89=reg101+reg89;
    reg48=reg76+reg48; reg76=reg81*reg102; reg101=reg6*reg66; reg84=reg106+reg84; reg106=reg23*reg43;
    T reg131=reg21*reg88; T reg132=reg53*reg102; reg94=reg99+reg94; reg59=reg70*reg59; reg99=reg55*reg46;
    reg90=reg77-reg90; reg77=2*reg50; reg105=reg79+reg105; reg79=reg23*reg21; reg54=reg30+reg54;
    reg30=reg57*reg71; T reg133=reg53*reg91; T reg134=reg6*reg4; T reg135=reg6*reg25; T reg136=reg57*reg48;
    reg98=reg44+reg98; reg99=reg129+reg99; reg44=reg55*reg48; reg129=0.5*elem.pos(1)[0]; T reg137=reg53*reg105;
    T reg138=reg63*reg94; reg112=reg74+reg112; reg74=reg59*reg16; T reg139=reg6*reg105; reg127=reg128+reg127;
    reg17=reg131*reg17; reg19=reg106*reg19; reg128=reg68*reg94; T reg140=reg63*reg89; T reg141=reg59*reg26;
    reg119=reg115+reg119; reg111=reg109+reg111; reg87=reg75-reg87; reg75=reg63*reg85; reg30=reg135+reg30;
    reg109=reg68*reg85; reg101=reg132+reg101; reg115=reg67*reg93; reg77=reg70*reg77; reg37=reg40+reg37;
    reg40=reg5*reg41; reg132=reg24*reg50; reg134=reg133+reg134; reg118=reg117+reg118; reg117=reg67*reg90;
    reg121=reg84+reg121; reg84=reg73*reg72; reg123=reg122+reg123; reg122=reg64*reg93; reg126=reg125+reg126;
    reg125=reg9*reg21; reg133=reg64*reg90; reg135=reg104*reg47; reg124=reg76+reg124; reg76=reg79*reg93;
    T reg142=reg103*reg27; T reg143=reg22*reg41; T reg144=reg9*reg50; T reg145=reg81*reg72; reg116=reg120+reg116;
    reg120=reg79*reg90; T reg146=reg68*reg89; reg39=reg39/reg33; reg15=reg110+reg15; reg32=reg106*reg32;
    reg110=reg61*reg25; T reg147=reg83*reg71; reg28=reg131*reg28; reg108=reg107+reg108; reg107=reg130*reg94;
    reg114=reg113+reg114; reg113=reg130*reg89; reg95=reg92-reg95; reg100=reg62-reg100; reg62=reg67*reg8;
    reg86=reg97+reg86; reg92=reg54*reg2; reg97=reg83*reg2; T reg148=0.5*elem.pos(0)[0]; T reg149=reg54*reg72;
    T reg150=0.5*elem.pos(0)[1]; T reg151=0.5*elem.pos(1)[1]; reg47=reg106*reg47; reg27=reg131*reg27; reg75=reg15+reg75;
    reg15=reg49*reg145; T reg152=reg40*reg39; T reg153=reg118*reg39; T reg154=reg132*reg39; T reg155=reg79*reg72;
    reg86=reg62-reg86; reg62=reg125*reg1; T reg156=reg37*reg1; T reg157=reg53*reg95; T reg158=reg6*reg100;
    reg115=reg101+reg115; reg101=reg52*reg97; reg117=reg134+reg117; reg134=reg52*reg92; T reg159=reg55*reg95;
    T reg160=reg57*reg100; reg122=reg123+reg122; reg123=reg51*reg97; reg133=reg126+reg133; reg126=reg51*reg92;
    reg76=reg124+reg76; reg124=reg54*reg97; reg120=reg116+reg120; reg116=reg54*reg92; reg138=reg99+reg138;
    reg99=reg49*reg84; reg44=reg137+reg44; reg137=reg63*reg121; reg140=reg127+reg140; reg127=reg49*reg149;
    reg109=reg30+reg109; reg30=reg58*reg145; reg128=reg111+reg128; reg111=reg58*reg84; reg136=reg139+reg136;
    reg139=reg68*reg121; reg146=reg98+reg146; reg98=reg58*reg149; reg147=reg110+reg147; reg110=reg130*reg85;
    reg107=reg108+reg107; reg108=reg54*reg84; T reg161=reg61*reg105; T reg162=reg83*reg48; reg113=reg114+reg113;
    reg114=reg54*reg149; T reg163=reg22*reg50; T reg164=reg9*reg41; T reg165=reg22*reg21; reg16=reg77*reg16;
    reg28=reg32+reg28; reg74=reg112+reg74; reg32=reg23*reg70; reg35=reg34+reg35; reg141=reg119+reg141;
    reg17=reg19+reg17; reg26=reg77*reg26; reg19=reg129+reg148; reg34=0.5*elem.pos(2)[1]; reg112=reg150+reg151;
    reg150=reg151-reg150; reg119=0.5*elem.pos(2)[0]; reg148=reg129-reg148; reg129=reg59*reg29; reg33=reg87/reg33;
    reg142=reg135+reg142; reg43=reg50*reg43; reg144=reg143+reg144; reg45=reg41*reg45; reg30=reg109+reg30;
    reg87=0.5*elem.pos(3)[1]; reg109=reg43*reg152; reg112=reg34-reg112; reg111=reg128+reg111; reg128=reg43*reg154;
    reg139=reg136+reg139; reg135=reg58*reg155; reg98=reg146+reg98; reg136=reg43*reg153; reg143=reg118*reg153;
    reg114=reg113+reg114; reg113=reg130*reg121; reg162=reg161+reg162; reg146=reg118*reg154; reg151=reg6*reg141;
    reg161=reg57*reg74; reg108=reg107+reg108; reg107=reg54*reg145; reg110=reg147+reg110; reg123=reg122+reg123;
    reg122=reg103*reg62; reg147=reg5*reg23; reg126=reg133+reg126; reg133=reg103*reg156; T reg166=reg81*reg95;
    T reg167=reg73*reg100; reg88=reg50*reg88; T reg168=reg21*reg70; reg38=reg41*reg38; reg124=reg76+reg124;
    reg76=reg37*reg62; reg116=reg120+reg116; reg120=reg37*reg156; T reg169=(*f.m).alpha_2*reg55; T reg170=(*f.m).alpha_1*reg6;
    T reg171=(*f.m).alpha_2*reg57; reg148=reg119+reg148; T reg172=0.5*elem.pos(3)[0]; T reg173=reg61*reg2; reg150=reg34+reg150;
    reg34=reg165*reg0; T reg174=reg35*reg0; T reg175=(*f.m).alpha_1*reg53; reg158=reg157+reg158; reg157=reg67*reg86;
    reg19=reg119-reg19; reg101=reg115+reg101; reg115=reg104*reg62; reg119=reg23*reg41; reg5=reg5*reg70;
    reg134=reg117+reg134; reg117=reg104*reg156; reg160=reg159+reg160; reg159=reg64*reg86; T reg176=reg53*reg141;
    T reg177=reg55*reg74; reg27=reg47+reg27; reg129=reg142+reg129; reg47=reg32*reg39; reg16=reg28+reg16;
    reg29=reg77*reg29; reg99=reg138+reg99; reg28=reg45*reg154; reg138=reg45*reg153; reg142=reg144*reg33;
    T reg178=reg163*reg33; reg26=reg17+reg26; reg17=reg164*reg33; reg127=reg140+reg127; reg15=reg75+reg15;
    reg137=reg44+reg137; reg44=reg49*reg155; reg75=reg45*reg152; reg120=reg116+reg120; reg29=reg27+reg29;
    reg27=reg51*reg173; reg116=reg37*reg72; reg76=reg124+reg76; reg177=reg176+reg177; reg124=reg63*reg129;
    reg140=reg35*reg34; reg159=reg160+reg159; reg160=reg6*reg26; reg5=reg119+reg5; reg9=reg9*reg70;
    reg119=reg21*reg41; reg176=reg57*reg16; T reg179=reg106*reg174; T reg180=reg144*reg142; reg143=reg114+reg143;
    reg133=reg126+reg133; reg114=reg131*reg174; reg126=reg61*reg141; T reg181=reg83*reg74; T reg182=reg168*reg33;
    T reg183=reg24*reg23; T reg184=reg54*reg155; reg113=reg162+reg113; reg138=reg127+reg138; reg127=reg144*reg178;
    reg146=reg108+reg146; reg108=reg38*reg142; reg162=reg131*reg34; reg122=reg123+reg122; reg167=reg166+reg167;
    reg123=reg79*reg86; reg24=reg24*reg70; reg23=reg23*reg50; reg166=reg118*reg152; reg107=reg110+reg107;
    reg110=(*f.m).alpha_3*reg63; T reg185=0.5*vectors[0][indices[0]+0]; reg148=reg148-reg172; T reg186=reg43*reg47; reg135=reg139+reg135;
    reg139=reg52*reg173; reg171=reg170+reg171; reg170=(*f.m).alpha_3*reg68; reg157=reg158+reg157; reg150=reg150-reg87;
    reg158=reg88*reg178; reg128=reg111+reg128; reg111=reg53*reg26; T reg187=reg55*reg16; T reg188=(*f.m).alpha_1*reg67;
    T reg189=(*f.m).alpha_2*reg64; T reg190=reg147*reg1; reg28=reg99+reg28; reg109=reg30+reg109; reg30=reg38*reg178;
    reg99=reg88*reg17; T reg191=reg35*reg174; reg117=reg134+reg117; reg134=0.5*vectors[0][indices[0]+1]; T reg192=reg68*reg129;
    reg161=reg151+reg161; reg75=reg15+reg75; reg15=0.5*vectors[0][indices[1]+1]; reg151=reg38*reg17; T reg193=reg88*reg142;
    reg136=reg98+reg136; reg98=0.5*vectors[0][indices[1]+0]; T reg194=reg45*reg47; reg44=reg137+reg44; reg137=reg106*reg34;
    reg115=reg101+reg115; reg87=reg112+reg87; reg19=reg172+reg19; reg169=reg175+reg169; reg194=reg44+reg194;
    reg44=reg130*reg129; reg101=reg38*reg182; reg181=reg126+reg181; reg112=reg87*reg148; reg30=reg28+reg30;
    reg180=reg143+reg180; reg28=reg183*reg0; reg162=reg122+reg162; reg122=reg61*reg26; reg126=reg83*reg16;
    reg24=reg23+reg24; reg23=reg103*reg190; reg139=reg157+reg139; reg143=reg104*reg190; reg27=reg159+reg27;
    reg157=reg150*reg19; reg137=reg115+reg137; reg151=reg75+reg151; reg179=reg117+reg179; reg124=reg177+reg124;
    reg140=reg76+reg140; reg22=reg22*reg70; reg75=reg21*reg50; reg9=reg119+reg9; reg76=reg68*reg29;
    reg176=reg160+reg176; reg191=reg120+reg191; reg115=reg58*reg116; reg192=reg161+reg192; reg117=reg49*reg116;
    reg193=reg136+reg193; reg119=reg15-reg134; reg120=reg88*reg182; reg186=reg135+reg186; reg110=reg169+reg110;
    reg170=reg171+reg170; reg158=reg128+reg158; reg187=reg111+reg187; reg189=reg188+reg189; reg111=(*f.m).alpha_3*reg65;
    reg99=reg109+reg99; reg109=reg63*reg29; reg61=(*f.m).alpha_1*reg61; reg83=(*f.m).alpha_2*reg83; reg15=reg134+reg15;
    reg128=reg118*reg47; reg184=reg113+reg184; reg113=0.5*vectors[0][indices[2]+1]; reg134=reg98+reg185; reg135=0.5*vectors[0][indices[2]+0];
    reg185=reg98-reg185; reg127=reg146+reg127; reg108=reg138+reg108; reg98=reg5*reg39; reg123=reg167+reg123;
    reg72=reg35*reg72; reg166=reg107+reg166; reg107=reg144*reg17; reg114=reg133+reg114; reg133=reg54*reg173;
    reg136=reg193*reg127; reg138=(*f.m).alpha_1*reg147; reg146=(*f.m).alpha_3*reg130; reg83=reg61+reg83; reg61=reg108*reg127;
    reg159=(*f.m).alpha_2*reg125; reg160=reg151*reg110; reg185=reg185+reg135; reg161=0.5*vectors[0][indices[3]+0]; reg111=reg189+reg111;
    reg167=reg30*reg180; reg169=reg30*reg170; reg134=reg135-reg134; reg135=reg158*reg180; reg22=reg75+reg22;
    reg23=reg27+reg23; reg27=reg131*reg28; reg75=reg162*reg191; reg171=reg37*reg190; reg133=reg123+reg133;
    reg123=reg137*reg191; reg172=reg114*reg140; reg175=reg158*reg170; reg119=reg113+reg119; reg177=reg99*reg110;
    reg157=reg112-reg157; reg112=0.5*vectors[0][indices[3]+1]; reg15=reg113-reg15; reg113=reg106*reg28; reg143=reg139+reg143;
    reg39=reg24*reg39; reg41=reg70*reg41; reg139=reg9*reg33; reg188=reg179*reg140; reg44=reg181+reg44;
    reg181=reg54*reg116; reg189=reg144*reg182; reg128=reg184+reg128; reg126=reg122+reg126; reg122=reg130*reg29;
    reg184=reg45*reg98; reg117=reg124+reg117; reg107=reg166+reg107; reg115=reg192+reg115; reg101=reg194+reg101;
    reg124=reg43*reg98; reg76=reg176+reg76; reg58=reg58*reg72; reg120=reg186+reg120; reg109=reg187+reg109;
    reg49=reg49*reg72; reg50=reg70*reg50; reg122=reg126+reg122; reg126=reg38*reg139; reg33=reg22*reg33;
    reg166=reg54*reg72; reg101=reg101*reg111; reg176=reg118*reg98; reg181=reg44+reg181; reg172=reg75-reg172;
    reg27=reg23+reg27; reg15=reg15+reg112; reg124=reg115+reg124; reg23=reg88*reg139; reg58=reg76+reg58;
    reg44=reg43*reg39; reg75=reg127*reg170; reg76=reg107*reg110; reg115=reg35*reg28; reg171=reg133+reg171;
    reg120=reg120*reg111; reg175=reg177+reg175; reg136=reg135-reg136; reg112=reg119-reg112; reg134=reg161+reg134;
    reg161=reg185-reg161; reg189=reg128+reg189; reg184=reg117+reg184; reg61=reg167-reg61; reg117=reg30*reg193;
    reg119=reg108*reg158; reg19=reg19/reg157; reg49=reg109+reg49; reg109=(*f.m).alpha_2*reg165; reg146=reg83+reg146;
    reg83=reg45*reg39; reg148=reg148/reg157; reg169=reg160+reg169; reg188=reg123-reg188; reg123=reg137*reg114;
    reg159=reg138+reg159; reg113=reg143+reg113; reg87=reg87/reg157; reg128=(*f.m).alpha_3*reg41; reg133=reg179*reg162;
    reg135=(*f.m).alpha_1*reg183; reg157=reg150/reg157; reg119=reg117-reg119; reg117=reg27*reg188; reg115=reg171+reg115;
    reg138=0.78867513459481286553*elem.pos(0)[1]; reg128=reg159+reg128; reg133=reg123-reg133; reg123=0.78867513459481286553*elem.pos(1)[0]; reg83=reg49+reg83;
    reg49=reg99*reg61; reg109=reg135+reg109; reg135=(*f.m).alpha_3*reg50; reg143=reg193*reg146; reg150=reg151*reg136;
    reg120=reg175+reg120; reg159=reg19*reg161; reg160=reg148*reg134; reg167=reg87*reg112; reg171=reg157*reg15;
    reg175=0.78867513459481286553*elem.pos(1)[1]; reg177=reg108*reg146; reg185=0.21132486540518713447*elem.pos(0)[1]; reg189=reg189*reg111; reg126=reg184+reg126;
    reg184=0.21132486540518713447*elem.pos(1)[0]; reg186=reg172*reg113; reg187=0.21132486540518713447*elem.pos(0)[0]; reg23=reg124+reg23; reg101=reg169+reg101;
    reg75=reg76+reg75; reg44=reg58+reg44; reg58=reg88*reg33; reg176=reg181+reg176; reg76=reg144*reg139;
    reg124=reg118*reg39; reg166=reg122+reg166; reg122=0.78867513459481286553*elem.pos(0)[0]; reg169=reg38*reg33; reg181=0.21132486540518713447*elem.pos(1)[1];
    reg192=reg180*reg146; reg159=reg160-reg159; reg189=reg75+reg189; reg75=reg122+reg184; reg160=reg181-reg185;
    reg194=reg179*reg115; T reg195=reg114*reg115; reg143=reg120+reg143; reg23=reg23*reg128; reg120=reg187+reg123;
    T reg196=reg27*reg191; T reg197=0.21132486540518713447*elem.pos(2)[0]; T reg198=reg113*reg191; reg169=reg83+reg169; reg58=reg44+reg58;
    reg44=reg108*reg107; reg185=reg185+reg175; reg83=0.78867513459481286553*elem.pos(2)[0]; T reg199=reg151*reg180; reg187=reg184-reg187;
    reg184=0.78867513459481286553*elem.pos(2)[1]; reg15=reg148*reg15; T reg200=reg193*reg107; T reg201=reg144*reg33; reg124=reg166+reg124;
    reg177=reg101+reg177; reg76=reg176+reg76; reg181=reg138+reg181; reg126=reg126*reg128; reg101=reg99*reg180;
    reg171=reg167-reg171; reg166=reg115*reg133; reg167=reg107*reg119; reg49=reg150-reg49; reg150=0.21132486540518713447*elem.pos(2)[1];
    reg117=reg186-reg117; reg161=reg87*reg161; reg112=reg19*reg112; reg135=reg109+reg135; reg134=reg157*reg134;
    reg185=reg184-reg185; reg109=0.21132486540518713447*elem.pos(3)[1]; reg176=reg137*reg115; reg194=reg198-reg194; reg112=reg15-reg112;
    elem.epsilon[0][1]=reg112; reg187=reg83+reg187; reg15=reg113*reg140; reg186=reg162*reg115; reg166=reg117+reg166;
    reg117=0.78867513459481286553*elem.pos(3)[0]; reg23=reg143+reg23; reg143=reg179*reg27; reg198=reg113*reg114; reg58=reg58*reg135;
    reg195=reg196-reg195; reg196=reg27*reg140; reg126=reg177+reg126; reg169=reg169*reg135; reg76=reg76*reg128;
    reg134=reg161-reg134; elem.epsilon[0][0]=reg134; reg167=reg49+reg167; reg49=reg108*reg99; reg138=reg175-reg138;
    reg161=reg151*reg193; reg201=reg124+reg201; reg124=reg99*reg127; reg200=reg101-reg200; reg101=reg63*reg95;
    reg175=reg68*reg100; reg181=reg150-reg181; reg177=reg158*reg107; reg122=reg123-reg122; reg123=reg30*reg107;
    T reg202=reg151*reg127; T reg203=reg63*reg102; T reg204=reg68*reg66; reg44=reg199-reg44; reg75=reg197-reg75;
    reg160=reg184+reg160; reg184=(*f.m).deltaT*reg170; reg199=0.78867513459481286553*elem.pos(3)[1]; T reg205=(*f.m).deltaT*reg110; reg120=reg83-reg120;
    reg192=reg189+reg192; reg83=0.21132486540518713447*elem.pos(3)[0]; reg171=reg159+reg171; reg159=0.78867513459481286553*PNODE(0).dep[1]; reg194=reg194/reg166;
    reg189=0.78867513459481286553*PNODE(1).dep[1]; reg49=reg161-reg49; reg188=reg188/reg166; reg122=reg197+reg122; reg181=reg181+reg199;
    reg186=reg196-reg186; reg123=reg202-reg123; reg176=reg15-reg176; reg15=(*f.m).deltaT*reg146; reg161=0.78867513459481286553*PNODE(0).dep[0];
    reg196=reg151*reg158; reg195=reg195/reg166; reg197=0.78867513459481286553*PNODE(1).dep[0]; reg202=reg30*reg99; T reg206=0.21132486540518713447*PNODE(0).dep[0];
    reg138=reg150+reg138; reg187=reg187-reg117; reg171=0.5*reg171; elem.epsilon[0][2]=reg171; reg201=reg201*reg135;
    reg76=reg192+reg76; reg150=reg134-reg205; reg192=reg112-reg184; reg120=reg120+reg83; reg199=reg160-reg199;
    reg58=reg23+reg58; reg117=reg75+reg117; reg23=reg113*reg162; reg143=reg198-reg143; reg75=reg68*reg4;
    reg160=reg63*reg91; reg204=reg203+reg204; reg198=reg65*reg93; reg203=reg137*reg27; reg169=reg126+reg169;
    reg136=reg136/reg167; reg126=0.21132486540518713447*PNODE(0).dep[1]; reg200=reg200/reg167; reg172=reg172/reg166; T reg207=0.21132486540518713447*PNODE(1).dep[1];
    reg175=reg101+reg175; reg101=reg65*reg86; reg44=reg44/reg167; reg177=reg124-reg177; reg185=reg185+reg109;
    reg61=reg61/reg167; reg124=0.21132486540518713447*PNODE(1).dep[0]; reg202=reg196-reg202; reg133=reg133/reg166; reg123=reg123/reg167;
    reg203=reg23-reg203; reg143=reg143/reg166; reg119=reg119/reg167; reg177=reg177/reg167; reg49=reg49/reg167;
    reg200=reg200*reg169; reg44=reg44*reg58; reg23=reg192*reg188; reg196=0.78867513459481286553*PNODE(2).dep[0]; T reg208=0.78867513459481286553*PNODE(2).dep[1];
    T reg209=reg117*reg199; reg101=reg175+reg101; reg175=reg69*reg173; T reg210=reg124-reg206; T reg211=reg171-reg15;
    T reg212=reg207+reg159; reg206=reg206+reg197; T reg213=0.21132486540518713447*PNODE(2).dep[1]; reg83=reg122-reg83; reg109=reg138-reg109;
    reg122=reg192*reg194; reg138=reg150*reg195; reg124=reg124+reg161; T reg214=0.21132486540518713447*PNODE(2).dep[0]; reg198=reg204+reg198;
    reg204=reg69*reg97; T reg215=reg181*reg187; T reg216=reg172*reg150; T reg217=reg185*reg187; reg75=reg160+reg75;
    reg201=reg76+reg201; reg76=reg65*reg90; reg160=reg199*reg120; T reg218=reg126+reg189; reg126=reg207-reg126;
    reg136=reg136*reg169; reg61=reg61*reg58; reg186=reg186/reg166; reg176=reg176/reg166; reg218=reg208-reg218;
    reg124=reg214-reg124; reg207=0.21132486540518713447*PNODE(3).dep[1]; reg58=reg123*reg58; reg138=reg122-reg138; reg23=reg216-reg23;
    reg122=reg211*reg143; reg169=reg177*reg169; reg206=reg196-reg206; reg210=reg196+reg210; reg123=0.78867513459481286553*PNODE(3).dep[0];
    reg212=reg213-reg212; reg161=reg197-reg161; reg177=0.78867513459481286553*PNODE(3).dep[1]; reg126=reg208+reg126; reg159=reg189-reg159;
    reg200=reg44-reg200; reg49=reg49*reg201; reg119=reg119*reg201; reg61=reg136-reg61; reg150=reg150*reg186;
    reg192=reg192*reg176; reg44=reg211*reg133; reg204=reg198+reg204; reg136=reg59*reg62; reg76=reg75+reg76;
    reg75=reg120*reg109; reg189=reg69*reg92; reg196=reg185*reg83; reg197=reg59*reg190; reg175=reg101+reg175;
    reg101=0.21132486540518713447*PNODE(3).dep[0]; reg160=reg217-reg160; reg167=reg202/reg167; reg209=reg215-reg209; reg166=reg203/reg166;
    reg211=reg211*reg166; reg192=reg150-reg192; reg61=reg119+reg61; reg161=reg214+reg161; reg49=reg200-reg49;
    reg212=reg177+reg212; reg201=reg167*reg201; reg58=reg169-reg58; reg119=reg181/reg209; reg150=reg199/reg209;
    reg122=reg138-reg122; reg75=reg196-reg75; reg138=reg187/reg209; reg167=reg117/reg209; reg124=reg123+reg124;
    reg169=reg181*reg83; reg196=reg120/reg160; reg206=reg101+reg206; reg198=reg77*reg28; reg197=reg175+reg197;
    reg175=reg117*reg109; reg159=reg213+reg159; reg123=reg210-reg123; reg23=reg44+reg23; reg187=reg187/reg160;
    reg199=reg199/reg160; reg44=reg185/reg160; reg218=reg218+reg207; reg177=reg126-reg177; reg189=reg76+reg189;
    reg76=reg59*reg156; reg126=1-(*f.m).resolution; reg136=reg204+reg136; reg200=reg77*reg34; reg202=reg150*reg212;
    reg58=reg201+reg58; reg201=reg77*reg174; reg76=reg189+reg76; reg189=reg27*reg23; reg203=reg196*reg123;
    reg175=reg169-reg175; reg169=reg119*reg177; reg170=(*f.m).resolution*reg170; reg120=reg120/reg75; reg204=reg83/reg75;
    reg208=reg109/reg75; reg210=reg167*reg123; reg185=reg185/reg75; reg200=reg136+reg200; reg136=reg138*reg124;
    reg61=reg126*reg61; reg213=reg132*reg66; reg198=reg197+reg198; reg197=reg40*reg102; reg214=reg187*reg206;
    reg110=(*f.m).resolution*reg110; reg207=reg159-reg207; reg101=reg161-reg101; reg192=reg211+reg192; reg49=reg126*reg49;
    reg159=reg113*reg23; reg161=reg199*reg218; reg211=reg162*reg122; reg215=reg40*reg95; reg216=reg132*reg100;
    reg217=reg137*reg122; T reg219=reg44*reg177; T reg220=reg208*reg218; T reg221=reg185*reg207; T reg222=reg196*reg177;
    T reg223=reg187*reg218; reg83=reg83/reg175; reg117=reg117/reg175; T reg224=reg199*reg206; reg201=reg76+reg201;
    reg76=reg44*reg123; T reg225=reg32*reg86; reg181=reg181/reg175; reg109=reg109/reg175; reg161=reg219-reg161;
    reg219=reg114*reg192; reg216=reg215+reg216; reg211=reg189+reg211; reg203=reg214-reg203; reg189=reg204*reg206;
    reg177=reg167*reg177; reg217=reg159+reg217; reg159=reg138*reg212; reg100=reg163*reg100; reg95=reg164*reg95;
    reg214=reg179*reg192; reg215=reg150*reg124; T reg226=reg198*reg23; T reg227=reg200*reg122; reg123=reg119*reg123;
    T reg228=reg120*reg101; T reg229=reg132*reg4; T reg230=reg40*reg91; reg146=(*f.m).resolution*reg146; reg210=reg136-reg210;
    reg213=reg197+reg213; reg136=reg32*reg93; reg49=reg170+reg49; reg61=reg110+reg61; reg102=reg164*reg102;
    reg66=reg163*reg66; reg58=reg126*reg58; reg202=reg169-reg202; reg110=reg83*reg124; reg169=reg120*reg207;
    reg218=reg204*reg218; reg206=reg208*reg206; reg170=reg185*reg101; reg220=reg221-reg220; reg177=reg159-reg177;
    reg159=reg118*reg173; reg202=reg210+reg202; reg197=(*f.m).resolution*reg172; reg228=reg189-reg228; reg225=reg216+reg225;
    reg136=reg213+reg136; reg189=reg118*reg97; reg215=reg123-reg215; reg123=(*f.m).resolution*reg188; reg210=reg201*reg192;
    reg4=reg163*reg4; reg91=reg164*reg91; reg49=(*f.m).deltaT*reg49; reg61=(*f.m).deltaT*reg61; reg58=reg146+reg58;
    reg93=reg168*reg93; reg66=reg102+reg66; reg229=reg230+reg229; reg102=reg32*reg90; reg151=reg151*reg126;
    reg227=reg226+reg227; reg214=reg217+reg214; reg100=reg95+reg100; reg86=reg168*reg86; reg30=reg30*reg126;
    reg99=reg99*reg126; reg127=reg127*reg126; reg107=reg107*reg126; reg158=reg158*reg126; reg95=reg117*reg101;
    reg146=(*f.m).resolution*reg176; reg213=(*f.m).resolution*reg186; reg216=(*f.m).resolution*reg194; reg217=(*f.m).resolution*reg195; reg219=reg211+reg219;
    reg161=reg203+reg161; reg224=reg76-reg224; reg222=reg223-reg222; reg76=reg181*reg207; reg203=reg109*reg212;
    reg95=reg110-reg95; reg169=reg218-reg169; reg206=reg170-reg206; reg203=reg76-reg203; reg101=reg181*reg101;
    reg214=reg205+reg214; reg220=reg228+reg220; reg124=reg109*reg124; reg212=reg83*reg212; reg207=reg117*reg207;
    reg219=reg184+reg219; reg111=(*f.m).deltaT*reg111; reg210=reg227+reg210; reg193=reg193*reg126; reg108=reg108*reg126;
    reg102=reg229+reg102; reg76=(*f.m).resolution*reg166; reg177=reg177-reg49; reg215=reg215-reg61; reg222=reg222-reg49;
    reg202=0.5*reg202; reg110=reg118*reg92; reg170=(*f.m).resolution*reg143; reg211=(*f.m).resolution*reg133; reg218=reg5*reg62;
    reg189=reg136+reg189; reg86=reg100+reg86; reg173=reg144*reg173; reg100=reg5*reg190; reg159=reg225+reg159;
    reg97=reg144*reg97; reg93=reg66+reg93; reg161=0.5*reg161; reg224=reg224-reg61; reg90=reg168*reg90;
    reg4=reg91+reg4; reg58=(*f.m).deltaT*reg58; reg146=reg127-reg146; reg107=reg213+reg107; reg126=reg180*reg126;
    reg158=reg216+reg158; reg151=reg197+reg151; reg217=reg99-reg217; reg123=reg30-reg123; reg161=reg161-reg58;
    reg97=reg93+reg97; reg62=reg9*reg62; reg30=reg107*reg224; reg92=reg144*reg92; reg90=reg4+reg90;
    reg100=reg159+reg100; reg4=reg24*reg28; reg108=reg211+reg108; reg66=reg158*reg222; reg91=reg217*reg224;
    reg218=reg189+reg218; reg93=reg24*reg34; reg203=reg95+reg203; reg124=reg101-reg124; reg126=reg76+reg126;
    reg76=reg146*reg222; reg207=reg212-reg207; reg170=reg193-reg170; reg95=reg5*reg156; reg110=reg102+reg110;
    reg202=reg202-reg58; reg99=reg107*reg215; reg101=reg146*reg177; reg210=reg111+reg210; reg220=0.5*reg220;
    reg102=reg214+reg219; reg127=reg123*reg222; reg136=reg151*reg224; reg206=reg206-reg61; reg159=reg123*reg177;
    reg180=reg217*reg215; reg190=reg9*reg190; reg173=reg86+reg173; reg86=reg158*reg177; reg189=reg151*reg215;
    reg169=reg169-reg49; reg193=reg108*reg202; reg197=reg170*reg202; reg86=reg180+reg86; reg180=reg115*reg23;
    reg102=reg210+reg102; reg211=reg140*reg122; reg159=reg189+reg159; reg207=reg207-reg49; reg92=reg90+reg92;
    reg156=reg9*reg156; reg90=reg126*reg202; reg101=reg99+reg101; reg28=reg22*reg28; reg99=reg217*reg206;
    reg189=reg158*reg169; reg190=reg173+reg190; reg173=reg146*reg169; reg212=reg151*reg206; reg213=reg123*reg169;
    reg62=reg97+reg62; reg34=reg22*reg34; reg97=reg170*reg161; reg216=reg107*reg206; reg4=reg100+reg4;
    reg124=reg124-reg61; reg95=reg110+reg95; reg76=reg30+reg76; reg30=reg126*reg161; reg100=reg24*reg174;
    reg203=0.5*reg203; reg93=reg218+reg93; reg110=reg108*reg161; reg220=reg220-reg58; reg127=reg136+reg127;
    reg66=reg91+reg66; reg211=reg180+reg211; reg100=reg95+reg100; reg34=reg62+reg34; reg62=reg191*reg192;
    reg102=reg102/3; reg156=reg92+reg156; reg174=reg22*reg174; reg91=reg4*reg23; reg92=reg93*reg122;
    reg28=reg190+reg28; reg189=reg99+reg189; reg95=reg170*reg220; reg99=reg217*reg124; reg136=reg158*reg207;
    reg90=reg101+reg90; reg213=reg212+reg213; reg101=reg108*reg220; reg197=reg86+reg197; reg86=reg151*reg124;
    reg180=reg123*reg207; reg190=reg146*reg207; reg193=reg159+reg193; reg159=reg107*reg124; reg203=reg203-reg58;
    reg30=reg76+reg30; reg76=reg126*reg220; reg173=reg216+reg173; reg97=reg66+reg97; reg110=reg127+reg110;
    reg101=reg213+reg101; reg23=reg28*reg23; reg95=reg189+reg95; reg122=reg34*reg122; reg66=reg100*reg192;
    reg92=reg91+reg92; reg90=2*reg90; reg62=reg211+reg62; reg197=reg177*reg197; reg219=reg219-reg102;
    reg214=reg214-reg102; reg193=reg215*reg193; reg76=reg173+reg76; reg97=reg222*reg97; reg110=reg224*reg110;
    reg136=reg99+reg136; reg174=reg156+reg174; reg91=reg170*reg203; reg99=reg108*reg203; reg30=2*reg30;
    reg190=reg159+reg190; reg180=reg86+reg180; reg86=reg126*reg203; reg66=reg92+reg66; reg62=reg15+reg62;
    reg91=reg136+reg91; reg90=reg202*reg90; reg102=reg210-reg102; reg219=pow(reg219,2); reg128=(*f.m).deltaT*reg128;
    reg214=pow(reg214,2); reg99=reg180+reg99; reg76=2*reg76; reg192=reg174*reg192; reg110=reg97+reg110;
    reg193=reg197+reg193; reg101=reg206*reg101; reg122=reg23+reg122; reg95=reg169*reg95; reg30=reg161*reg30;
    reg86=reg190+reg86; reg99=reg124*reg99; reg112=reg112-reg49; reg101=reg95+reg101; reg192=reg122+reg192;
    reg219=reg214+reg219; reg134=reg134-reg61; reg102=pow(reg102,2); reg135=(*f.m).deltaT*reg135; reg91=reg207*reg91;
    reg76=reg220*reg76; reg110=reg30+reg110; reg23=2*reg62; reg66=reg128+reg66; reg193=reg90+reg193;
    reg86=2*reg86; reg101=reg76+reg101; reg110=reg160*reg110; reg192=reg135+reg192; reg193=reg209*reg193;
    reg99=reg91+reg99; reg30=reg123*reg112; reg76=reg151*reg134; reg102=reg219+reg102; reg171=reg171-reg58;
    reg90=2*reg66; reg23=reg62*reg23; reg62=reg217*reg134; reg86=reg203*reg86; reg91=reg158*reg112;
    reg90=reg66*reg90; reg66=2*reg192; reg23=reg102+reg23; reg92=reg170*reg171; reg134=reg107*reg134;
    reg112=reg146*reg112; reg91=reg62+reg91; reg62=reg108*reg171; reg30=reg76+reg30; reg101=reg75*reg101;
    reg99=reg86+reg99; reg193=0.25*reg193; reg110=0.25*reg110; reg92=reg91+reg92; elem.sigma[0][1]=reg92;
    reg66=reg192*reg66; reg112=reg134+reg112; reg171=reg126*reg171; reg62=reg30+reg62; elem.sigma[0][0]=reg62;
    reg99=reg175*reg99; reg90=reg23+reg90; reg101=0.25*reg101; reg193=reg110+reg193; reg23=reg57*reg92;
    reg30=reg55*reg62; reg99=0.25*reg99; reg90=reg66+reg90; reg66=reg6*reg92; reg76=reg53*reg62;
    reg171=reg112+reg171; elem.sigma[0][2]=reg171; reg62=reg81*reg62; reg92=reg73*reg92; reg101=reg193+reg101;
    reg99=reg101+reg99; reg66=reg76+reg66; reg76=reg52*reg171; reg23=reg30+reg23; reg30=reg51*reg171;
    reg92=reg62+reg92; reg171=reg54*reg171; reg90=1.5*reg90; elem.sigma_local[0][0]=reg66+reg76; elem.sigma_local[0][1]=reg23+reg30;
    elem.sigma_local[0][2]=reg92+reg171; elem.sigma_von_mises=pow(reg90,0.5); elem.ener=reg99/2;
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
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg5=1.0/(*f.m).elastic_modulus_3; T reg6=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg7=pow((*f.m).v2[1],2); T reg8=pow((*f.m).v2[0],2); T reg9=pow((*f.m).v1[1],2); T reg10=pow((*f.m).v1[0],2); T reg11=reg2*reg3;
    reg9=reg10+reg9; reg10=pow((*f.m).v1[2],2); reg7=reg8+reg7; reg8=pow((*f.m).v2[2],2); T reg12=reg4*reg11;
    T reg13=reg5*reg11; T reg14=1.0/(*f.m).elastic_modulus_2; T reg15=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg16=reg6*reg11; reg10=reg9+reg10;
    reg8=reg7+reg8; reg7=reg6*reg12; reg9=reg13*reg14; T reg17=reg13*reg15; T reg18=reg6*reg16;
    reg10=pow(reg10,0.5); reg8=pow(reg8,0.5); T reg19=reg14*reg12; T reg20=reg15*reg16; reg7=reg17+reg7;
    T reg21=1.0/(*f.m).elastic_modulus_1; reg18=reg9-reg18; reg9=reg20+reg19; T reg22=(*f.m).v2[2]/reg8; T reg23=(*f.m).v2[1]/reg8;
    T reg24=reg21*reg18; T reg25=(*f.m).v1[2]/reg10; T reg26=(*f.m).v1[1]/reg10; T reg27=reg15*reg7; reg27=reg24-reg27;
    reg24=reg2*reg0; T reg28=reg4*reg16; T reg29=reg15*reg11; T reg30=reg4*reg3; T reg31=reg4*reg12;
    T reg32=reg5*reg3; T reg33=reg26*reg22; T reg34=reg25*reg23; T reg35=reg4*reg9; reg11=reg14*reg11;
    reg10=(*f.m).v1[0]/reg10; reg13=reg13*reg21; reg3=reg6*reg3; reg8=(*f.m).v2[0]/reg8; T reg36=2*reg8;
    T reg37=reg15*reg32; T reg38=2*reg10; reg28=reg17+reg28; reg17=reg33-reg34; T reg39=reg6*reg30;
    T reg40=reg25*reg8; reg35=reg27-reg35; reg27=reg4*reg29; T reg41=reg10*reg22; T reg42=reg2*reg1;
    T reg43=reg5*reg24; T reg44=reg4*reg24; reg32=reg14*reg32; reg16=reg21*reg16; reg31=reg13-reg31;
    reg12=reg15*reg12; reg24=reg6*reg24; reg13=reg4*reg11; T reg45=reg6*reg3; T reg46=reg6*reg24;
    reg27=reg16+reg27; T reg47=reg26*reg8; T reg48=reg10*reg23; T reg49=reg15*reg43; reg43=reg14*reg43;
    reg12=reg16+reg12; reg16=reg6*reg44; reg18=reg18/reg35; T reg50=reg40-reg41; T reg51=reg5*reg42;
    reg39=reg37+reg39; reg37=reg23*reg36; T reg52=2*reg17; T reg53=pow(reg10,2); T reg54=pow(reg23,2);
    T reg55=pow(reg26,2); T reg56=pow(reg8,2); T reg57=reg26*reg38; reg45=reg32-reg45; reg11=reg21*reg11;
    reg30=reg14*reg30; reg32=reg6*reg42; reg7=reg7/reg35; reg42=reg4*reg42; reg29=reg15*reg29;
    reg13=reg20+reg13; reg31=reg31/reg35; reg3=reg15*reg3; reg28=reg28/reg35; T reg58=reg55*reg18;
    T reg59=reg53*reg28; T reg60=reg54*reg7; reg24=reg15*reg24; T reg61=reg37*reg31; T reg62=reg53*reg18;
    reg13=reg13/reg35; T reg63=reg50*reg52; T reg64=reg56*reg31; T reg65=reg6*reg42; T reg66=reg3+reg30;
    T reg67=reg6*reg32; T reg68=reg15*reg51; reg51=reg14*reg51; T reg69=reg55*reg28; T reg70=reg54*reg31;
    reg16=reg49+reg16; reg12=reg12/reg35; reg46=reg43-reg46; reg43=pow(reg25,2); reg49=reg57*reg18;
    reg44=reg14*reg44; T reg71=reg48-reg47; reg39=reg15*reg39; T reg72=reg37*reg7; T reg73=reg57*reg28;
    T reg74=pow(reg22,2); reg27=reg27/reg35; T reg75=pow(reg50,2); T reg76=pow(reg17,2); reg45=reg21*reg45;
    reg9=reg9/reg35; T reg77=reg56*reg7; reg29=reg11-reg29; reg11=reg43*reg28; T reg78=reg75*reg12;
    reg70=reg69+reg70; reg29=reg29/reg35; reg69=reg74*reg31; reg16=reg15*reg16; T reg79=reg53*reg13;
    T reg80=reg24+reg44; T reg81=reg55*reg13; T reg82=reg56*reg27; T reg83=reg37*reg27; T reg84=reg57*reg13;
    T reg85=reg54*reg27; reg46=reg21*reg46; T reg86=pow(reg71,2); reg39=reg45-reg39; reg77=reg62+reg77;
    reg61=reg73+reg61; reg45=reg74*reg7; reg42=reg14*reg42; reg62=reg43*reg18; reg73=reg75*reg9;
    reg60=reg58+reg60; reg58=reg63*reg9; T reg87=reg63*reg12; reg72=reg49+reg72; reg64=reg59+reg64;
    reg49=reg76*reg12; reg65=reg68+reg65; reg66=reg4*reg66; reg59=reg76*reg9; reg67=reg51-reg67;
    reg32=reg15*reg32; reg59=reg77+reg59; reg66=reg39-reg66; reg82=reg79+reg82; reg87=reg61+reg87;
    reg39=reg86*reg12; reg51=reg76*reg29; reg69=reg11+reg69; reg78=reg70+reg78; reg49=reg64+reg49;
    reg85=reg81+reg85; reg11=reg75*reg29; reg61=reg43*reg13; reg64=reg74*reg27; reg58=reg72+reg58;
    reg68=reg86*reg9; reg83=reg84+reg83; reg70=reg8*reg23; reg72=reg10*reg26; reg16=reg46-reg16;
    reg80=reg4*reg80; reg67=reg21*reg67; reg65=reg15*reg65; reg46=reg32+reg42; reg73=reg60+reg73;
    reg45=reg62+reg45; reg60=2*reg23; reg62=reg22*reg36; reg77=2*reg26; reg79=reg25*reg38;
    reg81=reg63*reg29; reg84=reg55*reg58; T reg88=reg54*reg87; reg66=reg66/reg35; T reg89=reg62*reg7;
    T reg90=reg54*reg78; T reg91=reg55*reg73; T reg92=reg72*reg58; T reg93=reg70*reg87; reg58=reg53*reg58;
    reg87=reg56*reg87; reg39=reg69+reg39; reg69=reg79*reg28; T reg94=reg62*reg31; reg51=reg82+reg51;
    reg11=reg85+reg11; reg64=reg61+reg64; reg61=reg86*reg29; reg81=reg83+reg81; reg82=reg53*reg59;
    reg83=reg56*reg49; reg85=reg53*reg73; T reg95=reg56*reg78; T reg96=reg54*reg49; T reg97=reg55*reg59;
    reg65=reg67-reg65; reg80=reg16-reg80; reg46=reg4*reg46; reg16=reg8*reg38; reg67=reg10*reg50;
    T reg98=reg26*reg17; T reg99=reg23*reg77; reg48=reg47+reg48; reg47=reg26*reg23; T reg100=reg10*reg8;
    T reg101=2*reg50; T reg102=reg79*reg18; reg68=reg45+reg68; reg45=reg50*reg17; T reg103=2*reg25;
    reg78=reg70*reg78; reg73=reg72*reg73; T reg104=reg25*reg77; T reg105=reg22*reg60; reg52=reg71*reg52;
    reg80=reg80/reg35; T reg106=reg76*reg51; reg83=reg82+reg83; reg95=reg85+reg95; reg82=reg76*reg11;
    reg85=reg55*reg68; T reg107=reg54*reg39; T reg108=reg48*reg66; T reg109=reg47*reg66; T reg110=reg100*reg66;
    T reg111=reg62*reg27; T reg112=reg79*reg13; reg61=reg64+reg61; reg46=reg65-reg46; reg31=reg105*reg31;
    reg28=reg104*reg28; reg64=reg52*reg12; reg94=reg69+reg94; reg7=reg105*reg7; reg18=reg104*reg18;
    reg65=reg52*reg9; reg89=reg102+reg89; reg69=reg8*reg50; reg102=reg23*reg17; T reg113=reg99*reg14;
    reg67=reg98+reg67; reg98=reg16*reg15; T reg114=reg26*reg50; T reg115=reg10*reg17; T reg116=reg54*reg14;
    T reg117=reg25*reg22; reg101=reg71*reg101; reg78=reg73+reg78; reg73=reg45*reg11; T reg118=reg56*reg15;
    T reg119=reg99*reg15; reg90=reg91+reg90; reg91=reg16*reg21; T reg120=reg54*reg15; T reg121=reg56*reg21;
    T reg122=reg76*reg81; reg87=reg58+reg87; reg93=reg92+reg93; reg58=reg45*reg81; reg88=reg84+reg88;
    reg81=reg75*reg81; reg96=reg97+reg96; reg84=reg75*reg51; reg11=reg75*reg11; reg103=reg22*reg103;
    reg92=reg53*reg68; reg97=reg56*reg39; reg59=reg72*reg59; reg49=reg70*reg49; reg31=reg28+reg31;
    reg28=reg103*reg4; T reg123=reg16*reg109; reg64=reg94+reg64; reg84=reg96+reg84; reg106=reg83+reg106;
    reg120=reg121-reg120; reg83=reg74*reg4; reg94=reg16*reg110; reg82=reg95+reg82; reg95=reg16*reg108;
    reg122=reg87+reg122; reg27=reg105*reg27; reg87=reg48*reg108; reg58=reg93+reg58; reg93=reg117*reg66;
    reg96=reg16*reg4; reg121=reg56*reg4; reg98=reg113-reg98; reg113=reg99*reg110; T reg124=reg103*reg6;
    T reg125=reg115*reg80; T reg126=reg54*reg6; reg111=reg112+reg111; reg112=reg114*reg80; reg118=reg116-reg118;
    reg116=reg67*reg80; T reg127=reg75*reg61; reg107=reg85+reg107; reg85=reg52*reg29; reg13=reg104*reg13;
    T reg128=reg74*reg6; T reg129=reg99*reg6; reg119=reg91-reg119; reg12=reg101*reg12; reg73=reg78+reg73;
    reg78=reg48*reg109; reg77=reg50*reg77; reg91=reg25*reg71; reg38=reg17*reg38; reg97=reg92+reg97;
    reg92=reg76*reg61; reg49=reg59+reg49; reg65=reg89+reg65; reg69=reg102+reg69; reg51=reg45*reg51;
    reg68=reg72*reg68; reg7=reg18+reg7; reg9=reg101*reg9; reg11=reg90+reg11; reg39=reg70*reg39;
    reg109=reg99*reg109; reg81=reg88+reg81; reg108=reg99*reg108; reg18=reg8*reg17; reg35=reg46/reg35;
    reg46=reg23*reg50; reg36=reg17*reg36; reg29=reg101*reg29; reg59=reg91*reg80; reg88=reg38*reg112;
    reg89=reg99*reg93; reg123=reg82+reg123; reg82=reg53*reg65; reg126=reg121+reg126; reg109=reg11+reg109;
    reg11=reg77*reg112; reg108=reg81+reg108; reg81=reg38*reg125; reg94=reg106+reg94; reg90=reg77*reg116;
    reg102=reg56*reg64; reg106=(*f.m).alpha_1*reg53; reg129=reg96+reg129; reg96=reg55*reg65; reg21=reg53*reg21;
    reg103=reg103*reg5; reg121=reg54*reg64; T reg130=reg69*reg35; reg60=reg50*reg60; T reg131=reg46*reg35;
    T reg132=reg18*reg35; reg127=reg107+reg127; reg107=reg53*reg15; reg28=reg119-reg28; reg78=reg73+reg78;
    reg12=reg31+reg12; reg112=reg67*reg112; reg31=reg16*reg93; reg83=reg120-reg83; reg73=(*f.m).alpha_2*reg56;
    reg9=reg7+reg9; reg15=reg55*reg15; reg39=reg68+reg39; reg7=(*f.m).alpha_1*reg55; reg68=reg67*reg116;
    reg87=reg58+reg87; reg58=(*f.m).alpha_2*reg54; reg61=reg45*reg61; reg27=reg13+reg27; reg13=reg22*reg71;
    reg85=reg111+reg85; reg111=reg74*reg5; reg51=reg49+reg51; reg110=reg48*reg110; reg113=reg84+reg113;
    reg49=reg77*reg125; reg124=reg98-reg124; reg84=reg10*reg71; reg98=reg25*reg17; reg128=reg118-reg128;
    reg92=reg97+reg92; reg40=reg41+reg40; reg116=reg38*reg116; reg14=reg55*reg14; reg95=reg122+reg95;
    reg49=reg113+reg49; reg41=reg60*reg132; reg97=reg56*reg12; reg113=reg53*reg9; reg64=reg70*reg64;
    reg118=reg47*reg128; reg119=reg100*reg83; reg120=reg43*reg4; reg15=reg21-reg15; reg21=reg43*reg6;
    reg107=reg14-reg107; reg4=reg53*reg4; reg6=reg55*reg6; reg126=reg111-reg126; reg129=reg103-reg129;
    reg14=reg53*reg83; reg103=reg55*reg128; reg111=reg53*reg28; reg122=reg55*reg124; T reg133=reg56*reg83;
    T reg134=reg54*reg128; T reg135=reg54*reg124; T reg136=reg56*reg28; reg11=reg109+reg11; reg109=reg60*reg131;
    reg89=reg127+reg89; reg127=reg77*reg59; reg90=reg108+reg90; reg108=reg60*reg130; reg121=reg96+reg121;
    reg96=reg75*reg85; T reg137=reg55*reg9; T reg138=reg54*reg12; reg110=reg51+reg110; reg125=reg67*reg125;
    reg112=reg78+reg112; reg51=reg69*reg131; reg61=reg39+reg61; reg93=reg48*reg93; reg68=reg87+reg68;
    reg39=reg69*reg130; reg65=reg72*reg65; reg29=reg27+reg29; reg31=reg92+reg31; reg27=reg38*reg59;
    reg131=reg36*reg131; reg78=reg47*reg124; reg87=(*f.m).alpha_1*reg43; reg92=(*f.m).alpha_2*reg74; reg116=reg95+reg116;
    reg130=reg36*reg130; reg88=reg123+reg88; reg102=reg82+reg102; reg82=reg76*reg85; reg95=(*f.m).alpha_3*reg75;
    reg33=reg34+reg33; reg84=reg98+reg84; reg34=reg25*reg50; reg98=reg26*reg71; reg123=reg36*reg132;
    T reg139=reg100*reg28; T reg140=reg22*reg17; T reg141=reg8*reg71; reg81=reg94+reg81; reg58=reg7+reg58;
    reg7=reg40*reg66; reg94=(*f.m).alpha_3*reg76; reg73=reg106+reg73; reg106=reg13*reg35; T reg142=reg117*reg129;
    reg78=reg139+reg78; reg132=reg69*reg132; reg125=reg110+reg125; reg98=reg34+reg98; reg141=reg140+reg141;
    reg118=reg119+reg118; reg34=reg43*reg129; reg110=reg75*reg29; reg138=reg137+reg138; reg119=reg22*reg50;
    reg137=reg23*reg71; reg21=reg107-reg21; reg5=reg43*reg5; reg6=reg4+reg6; reg120=reg15-reg120;
    reg4=reg70*reg2; reg15=reg48*reg2; reg12=reg70*reg12; reg9=reg72*reg9; reg103=reg14+reg103;
    reg85=reg45*reg85; reg64=reg65+reg64; reg39=reg68+reg39; reg14=reg43*reg126; reg70=(*f.m).alpha_2*reg70;
    reg65=(*f.m).alpha_1*reg72; reg68=(*f.m).alpha_3*reg86; reg92=reg87+reg92; reg59=reg67*reg59; reg93=reg61+reg93;
    reg95=reg58+reg95; reg94=reg73+reg94; reg51=reg112+reg51; reg122=reg111+reg122; reg58=reg60*reg106;
    reg127=reg89+reg127; reg66=reg33*reg66; reg61=reg74*reg126; reg73=reg84*reg80; reg134=reg133+reg134;
    reg135=reg136+reg135; reg87=reg117*reg126; reg89=reg74*reg129; reg109=reg11+reg109; reg123=reg81+reg123;
    reg131=reg88+reg131; reg27=reg31+reg27; reg41=reg49+reg41; reg11=reg16*reg7; reg31=reg36*reg106;
    reg82=reg102+reg82; reg130=reg116+reg130; reg49=reg99*reg7; reg8=reg8*reg22; reg10=reg10*reg25;
    reg97=reg113+reg97; reg96=reg121+reg96; reg108=reg90+reg108; reg81=reg76*reg29; reg88=reg131*reg39;
    reg90=reg108*reg51; reg102=(*f.m).alpha_3*reg45; reg70=reg65+reg70; reg65=reg41*reg94; reg107=reg37*reg4;
    reg111=reg130*reg51; reg11=reg82+reg11; reg6=reg5-reg6; reg68=reg92+reg68; reg5=reg131*reg95;
    reg82=reg123*reg94; reg92=(*f.m).alpha_1*reg10; reg112=(*f.m).alpha_2*reg8; reg113=reg109*reg39; reg116=reg109*reg95;
    reg121=reg38*reg73; reg85=reg64+reg85; reg7=reg48*reg7; reg61=reg134+reg61; reg16=reg16*reg66;
    reg81=reg97+reg81; reg12=reg9+reg12; reg29=reg45*reg29; reg9=reg55*reg21; reg45=reg53*reg120;
    reg64=reg40*reg1; reg8=reg8*reg1; reg34=reg122+reg34; reg97=reg57*reg15; reg17=reg71*reg17;
    reg25=reg26*reg25; reg22=reg23*reg22; reg23=reg54*reg21; reg99=reg99*reg66; reg110=reg138+reg110;
    reg137=reg119+reg137; reg26=reg37*reg15; reg89=reg135+reg89; reg119=reg48*reg4; reg122=reg56*reg120;
    reg58=reg127+reg58; reg31=reg27+reg31; reg14=reg103+reg14; reg106=reg69*reg106; reg59=reg93+reg59;
    reg87=reg118+reg87; reg27=reg141*reg35; reg80=reg98*reg80; reg93=reg57*reg4; reg49=reg96+reg49;
    reg96=reg77*reg73; reg103=reg48*reg15; reg142=reg78+reg142; reg132=reg125+reg132; reg107=reg61+reg107;
    reg61=reg79*reg64; reg78=reg62*reg8; reg97=reg34+reg97; reg34=reg62*reg64; reg31=reg31*reg68;
    reg5=reg82+reg5; reg116=reg65+reg116; reg58=reg58*reg68; reg2=reg72*reg2; reg65=reg22*reg0;
    reg72=reg33*reg0; reg23=reg122+reg23; reg74=reg74*reg6; reg9=reg45+reg9; reg26=reg89+reg26;
    reg43=reg43*reg6; reg45=reg47*reg21; reg82=reg100*reg120; reg89=reg132*reg94; reg118=reg51*reg95;
    reg122=reg79*reg8; reg93=reg14+reg93; reg17=(*f.m).alpha_3*reg17; reg14=reg40*reg64; reg119=reg87+reg119;
    reg87=(*f.m).alpha_1*reg25; reg22=(*f.m).alpha_2*reg22; reg103=reg142+reg103; reg66=reg48*reg66; reg29=reg12+reg29;
    reg77=reg77*reg80; reg99=reg110+reg99; reg112=reg92+reg112; reg38=reg38*reg80; reg106=reg59+reg106;
    reg16=reg81+reg16; reg35=reg137*reg35; reg73=reg67*reg73; reg7=reg85+reg7; reg12=reg40*reg8;
    reg59=reg36*reg27; reg81=reg130*reg109; reg50=reg71*reg50; reg71=reg131*reg108; reg96=reg49+reg96;
    reg111=reg88-reg111; reg49=reg60*reg27; reg102=reg70+reg102; reg90=reg113-reg90; reg121=reg11+reg121;
    reg11=reg33*reg72; reg14=reg103+reg14; reg70=reg105*reg72; reg122=reg93+reg122; reg85=reg104*reg65;
    reg34=reg26+reg34; reg118=reg89+reg118; reg60=reg60*reg35; reg77=reg99+reg77; reg49=reg96+reg49;
    reg61=reg97+reg61; reg26=reg104*reg72; reg88=reg33*reg65; reg74=reg23+reg74; reg23=reg37*reg2;
    reg12=reg119+reg12; reg81=reg71-reg81; reg71=reg105*reg65; reg78=reg107+reg78; reg59=reg121+reg59;
    reg31=reg5+reg31; reg5=reg130*reg102; reg89=reg41*reg111; reg92=reg123*reg90; reg50=(*f.m).alpha_3*reg50;
    reg22=reg87+reg22; reg58=reg116+reg58; reg87=reg108*reg102; reg17=reg112+reg17; reg1=reg10*reg1;
    reg38=reg16+reg38; reg117=reg117*reg6; reg80=reg67*reg80; reg106=reg106*reg68; reg66=reg29+reg66;
    reg43=reg9+reg43; reg9=reg57*reg2; reg45=reg82+reg45; reg27=reg69*reg27; reg73=reg7+reg73;
    reg36=reg36*reg35; reg106=reg118+reg106; reg7=reg39*reg102; reg88=reg12+reg88; reg11=reg14+reg11;
    reg49=reg49*reg17; reg70=reg34+reg70; reg87=reg58+reg87; reg10=reg48*reg2; reg117=reg45+reg117;
    reg50=reg22+reg50; reg59=reg59*reg17; reg5=reg31+reg5; reg89=reg92-reg89; reg12=reg132*reg81;
    reg14=reg130*reg132; reg16=reg123*reg39; reg22=reg108*reg132; reg29=reg41*reg39; reg36=reg38+reg36;
    reg60=reg77+reg60; reg27=reg73+reg27; reg80=reg66+reg80; reg35=reg69*reg35; reg0=reg25*reg0;
    reg9=reg43+reg9; reg79=reg79*reg1; reg25=1-var_inter[0]; reg85=reg122+reg85; reg26=reg61+reg26;
    reg23=reg74+reg23; reg62=reg62*reg1; reg31=1-var_inter[1]; reg71=reg78+reg71; reg34=reg31*elem.pos(1)[0];
    reg38=reg123*reg51; reg43=reg109*reg132; reg22=reg29-reg22; reg29=reg41*reg51; reg45=reg31*elem.pos(0)[1];
    reg58=reg31*elem.pos(1)[1]; reg12=reg89+reg12; reg61=elem.pos(0)[0]*reg25; reg66=elem.pos(1)[0]*var_inter[0]; reg73=reg25*elem.pos(0)[1];
    reg74=elem.pos(1)[1]*var_inter[0]; reg14=reg16-reg14; reg16=reg131*reg132; reg77=reg123*reg108; reg78=reg130*reg41;
    reg82=reg31*elem.pos(0)[0]; reg59=reg5+reg59; reg36=reg36*reg50; reg27=reg27*reg17; reg7=reg106+reg7;
    reg49=reg87+reg49; reg60=reg60*reg50; reg10=reg117+reg10; reg35=reg80+reg35; reg40=reg40*reg1;
    reg79=reg9+reg79; reg5=reg71*reg11; reg105=reg105*reg0; reg62=reg23+reg62; reg9=reg85*reg11;
    reg23=reg26*reg88; reg104=reg104*reg0; reg80=reg70*reg88; reg87=reg61+reg66; reg23=reg9-reg23;
    reg9=reg85*reg70; reg34=reg34-reg82; reg89=reg26*reg71; reg92=reg131*reg41; reg60=reg49+reg60;
    reg105=reg62+reg105; reg111=reg111/reg12; reg14=reg14/reg12; reg43=reg29-reg43; reg104=reg79+reg104;
    reg78=reg77-reg78; reg29=reg123*reg109; reg80=reg5-reg80; reg5=elem.pos(2)[0]*var_inter[1]; reg22=reg22/reg12;
    reg90=reg90/reg12; reg49=elem.pos(2)[1]*var_inter[1]; reg36=reg59+reg36; reg58=reg58-reg45; reg16=reg38-reg16;
    reg33=reg33*reg0; reg38=elem.pos(2)[0]*var_inter[0]; reg40=reg10+reg40; reg10=reg73+reg74; reg35=reg35*reg50;
    reg27=reg7+reg27; reg7=elem.pos(2)[1]*var_inter[0]; reg58=reg49+reg58; reg38=reg38-reg87; reg33=reg40+reg33;
    reg43=reg43/reg12; reg40=elem.pos(3)[1]*var_inter[1]; reg49=elem.pos(3)[0]*var_inter[1]; reg34=reg5+reg34; reg92=reg29-reg92;
    reg111=reg111*reg60; reg90=reg90*reg36; reg35=reg27+reg35; reg81=reg81/reg12; reg5=elem.pos(3)[1]*reg25;
    reg7=reg7-reg10; reg27=reg105*reg23; reg29=elem.pos(3)[0]*reg25; reg59=reg80*reg104; reg22=reg22*reg36;
    reg89=reg9-reg89; reg16=reg16/reg12; reg14=reg14*reg60; reg78=reg78/reg12; reg34=reg34-reg49;
    reg78=reg78*reg35; reg111=reg90-reg111; reg81=reg81*reg35; reg22=reg14-reg22; reg36=reg43*reg36;
    reg60=reg16*reg60; reg58=reg58-reg40; reg5=reg7+reg5; reg29=reg38+reg29; reg27=reg59-reg27;
    reg7=reg33*reg89; reg12=reg92/reg12; reg9=reg105*reg88; reg14=reg71*reg33; reg16=reg104*reg88;
    reg38=reg85*reg33; reg60=reg36-reg60; reg36=reg58*reg29; reg35=reg12*reg35; reg78=reg22-reg78;
    reg111=reg81+reg111; reg12=reg34*reg5; reg7=reg27+reg7; reg22=reg105*reg11; reg27=reg70*reg33;
    reg43=reg104*reg11; reg14=reg9-reg14; reg9=reg26*reg33; reg59=1-(*f.m).resolution; reg38=reg16-reg38;
    reg16=reg104*reg71; reg62=reg85*reg105; reg78=reg59*reg78; reg111=reg59*reg111; reg77=(*f.m).resolution*reg95;
    reg60=reg35+reg60; reg35=reg26*reg105; reg62=reg16-reg62; reg16=reg104*reg70; reg38=reg38/reg7;
    reg9=reg43-reg9; reg14=reg14/reg7; reg27=reg22-reg27; reg22=(*f.m).resolution*reg94; reg36=reg12-reg36;
    reg5=reg5/reg36; reg12=(*f.m).resolution*reg102; reg43=(*f.m).resolution*reg14; reg79=(*f.m).resolution*reg38; reg34=reg34/reg36;
    reg60=reg59*reg60; reg9=reg9/reg7; reg23=reg23/reg7; reg111=reg22+reg111; reg78=reg77+reg78;
    reg27=reg27/reg7; reg58=reg58/reg36; reg29=reg29/reg36; reg132=reg132*reg59; reg51=reg51*reg59;
    reg35=reg16-reg35; reg62=reg62/reg7; reg80=reg80/reg7; reg111=(*f.m).deltaT*reg111; reg123=reg123*reg59;
    reg131=reg131*reg59; reg41=reg41*reg59; reg109=reg109*reg59; reg16=reg31*reg29; reg39=reg39*reg59;
    reg60=reg12+reg60; reg132=reg43+reg132; reg79=reg51-reg79; reg12=var_inter[0]*reg58; reg22=var_inter[0]*reg34;
    reg43=(*f.m).resolution*reg23; reg51=var_inter[1]*reg5; reg77=var_inter[1]*reg29; reg81=reg31*reg5; reg90=(*f.m).resolution*reg80;
    reg92=(*f.m).resolution*reg27; reg35=reg35/reg7; reg93=(*f.m).resolution*reg9; reg96=(*f.m).resolution*reg62; reg7=reg89/reg7;
    reg78=(*f.m).deltaT*reg78; reg89=reg25*reg34; reg97=reg25*reg58; reg99=reg51-reg12; reg103=(*f.m).resolution*reg35;
    reg60=(*f.m).deltaT*reg60; reg106=reg132*reg111; reg107=reg16+reg22; reg110=reg12+reg81; reg112=reg22-reg77;
    reg39=reg96+reg39; reg96=reg97-reg81; reg130=reg130*reg59; reg113=reg79*reg78; reg109=reg93+reg109;
    reg92=reg41-reg92; reg41=reg89+reg77; reg43=reg131-reg43; reg123=reg90+reg123; reg90=reg16-reg89;
    reg93=reg97+reg51; reg59=reg108*reg59; reg108=(*f.m).resolution*reg7; reg116=reg109*reg78; reg117=reg92*reg111;
    reg118=reg106+reg113; reg119=reg39*reg60; reg121=reg43*reg78; reg122=reg123*reg111; reg125=0.5*reg99;
    reg127=0.5*reg112; reg131=0.5*reg110; reg133=0.5*reg107; reg103=reg59-reg103; reg130=reg108+reg130;
    reg59=0.5*reg96; reg108=0.5*reg90; reg134=0.5*reg93; reg135=0.5*reg41; reg136=reg110*reg132;
    reg138=reg93*reg132; reg139=reg135*reg39; reg140=reg59*reg39; reg142=reg103*reg60; T reg143=reg117+reg116;
    T reg144=reg112*reg79; T reg145=reg125*reg39; T reg146=reg130*reg60; T reg147=reg122+reg121; T reg148=reg41*reg79;
    T reg149=reg134*reg39; T reg150=reg118+reg119; T reg151=reg96*reg132; T reg152=reg108*reg39; T reg153=reg127*reg39;
    T reg154=reg99*reg132; T reg155=reg131*reg39; T reg156=reg107*reg79; T reg157=reg133*reg39; T reg158=reg90*reg79;
    T reg159=reg112*reg43; T reg160=reg99*reg123; T reg161=reg125*reg130; T reg162=reg134*reg103; T reg163=reg41*reg109;
    T reg164=reg135*reg103; T reg165=reg93*reg92; T reg166=reg125*reg103; T reg167=reg112*reg109; T reg168=reg127*reg103;
    T reg169=reg131*reg103; T reg170=reg99*reg92; T reg171=reg107*reg109; T reg172=reg133*reg103; T reg173=reg59*reg103;
    reg145=reg144+reg145; reg144=reg127*reg130; T reg174=2*reg150; reg136=reg136-reg157; reg153=reg154+reg153;
    reg154=reg59*reg130; T reg175=reg90*reg43; T reg176=reg131*reg130; T reg177=reg147+reg146; T reg178=reg107*reg43;
    reg155=reg155-reg156; T reg179=reg110*reg123; T reg180=reg143+reg142; T reg181=reg133*reg130; reg140=reg158+reg140;
    reg158=reg134*reg130; T reg182=reg41*reg43; reg148=reg148-reg149; T reg183=reg135*reg130; T reg184=reg93*reg123;
    T reg185=reg96*reg92; T reg186=reg108*reg103; reg139=reg139-reg138; T reg187=reg108*reg130; T reg188=reg96*reg123;
    reg152=reg151+reg152; reg151=reg90*reg109; T reg189=reg110*reg92; reg166=reg167+reg166; reg154=reg175+reg154;
    reg140=2*reg140; reg144=reg160+reg144; reg164=reg164-reg165; reg139=2*reg139; reg160=reg134*reg174;
    reg167=reg41*reg180; reg175=reg93*reg177; T reg190=reg135*reg174; reg163=reg163-reg162; reg182=reg182-reg158;
    reg173=reg151+reg173; reg151=reg107*reg180; T reg191=reg131*reg174; T reg192=reg133*reg174; T reg193=reg110*reg177;
    reg148=2*reg148; reg145=2*reg145; reg183=reg183-reg184; reg153=2*reg153; reg152=2*reg152;
    reg176=reg176-reg178; reg169=reg169-reg171; reg189=reg189-reg172; reg155=2*reg155; reg187=reg188+reg187;
    reg168=reg170+reg168; reg179=reg179-reg181; reg186=reg185+reg186; reg136=2*reg136; reg161=reg159+reg161;
    reg159=reg135*reg148; reg170=reg112*reg166; reg185=reg125*reg145; reg188=reg112*reg164; T reg194=reg93*reg183;
    T reg195=reg93*reg179; T reg196=reg135*reg139; T reg197=reg125*reg139; T reg198=reg93*reg161; T reg199=reg135*reg155;
    T reg200=reg112*reg163; T reg201=reg135*reg136; T reg202=reg135*reg145; T reg203=reg125*reg148; T reg204=reg93*reg154;
    T reg205=reg93*reg144; T reg206=reg135*reg152; T reg207=reg135*reg140; T reg208=reg135*reg153; T reg209=reg93*reg187;
    T reg210=reg93*reg176; T reg211=reg99*reg187; T reg212=reg127*reg152; T reg213=reg127*reg136; T reg214=reg99*reg154;
    T reg215=reg127*reg140; T reg216=reg99*reg179; T reg217=reg99*reg176; T reg218=reg99*reg144; T reg219=reg127*reg155;
    T reg220=reg127*reg153; T reg221=reg99*reg161; T reg222=reg127*reg145; T reg223=reg99*reg182; T reg224=reg127*reg148;
    T reg225=reg112*reg186; T reg226=reg125*reg152; T reg227=reg112*reg173; T reg228=reg125*reg140; T reg229=reg112*reg189;
    T reg230=reg125*reg136; T reg231=reg112*reg169; T reg232=reg125*reg155; T reg233=reg112*reg168; T reg234=reg125*reg153;
    T reg235=reg134*reg136; T reg236=reg90*reg180; T reg237=reg59*reg174; T reg238=reg192-reg193; T reg239=reg151-reg191;
    T reg240=reg99*reg177; T reg241=reg127*reg174; T reg242=reg112*reg180; T reg243=reg125*reg174; T reg244=reg175-reg190;
    T reg245=reg160-reg167; T reg246=reg127*reg139; T reg247=reg96*reg144; T reg248=reg108*reg140; T reg249=reg96*reg154;
    T reg250=reg96*reg161; T reg251=reg108*reg136; T reg252=reg96*reg179; T reg253=reg96*reg187; T reg254=reg108*reg155;
    T reg255=reg96*reg176; T reg256=reg108*reg152; T reg257=reg99*reg183; T reg258=reg108*reg153; T reg259=reg93*reg182;
    T reg260=reg41*reg186; T reg261=reg134*reg152; T reg262=reg41*reg173; T reg263=reg134*reg140; T reg264=reg41*reg169;
    T reg265=reg134*reg155; T reg266=reg41*reg168; T reg267=reg134*reg153; T reg268=reg41*reg166; T reg269=reg134*reg145;
    T reg270=reg41*reg164; T reg271=reg134*reg139; T reg272=reg41*reg163; T reg273=reg134*reg148; T reg274=reg96*reg182;
    T reg275=reg108*reg148; T reg276=reg108*reg174; T reg277=reg41*reg189; T reg278=reg96*reg177; T reg279=reg96*reg183;
    T reg280=reg108*reg145; T reg281=reg108*reg139; T reg282=reg107*reg166; T reg283=reg90*reg164; T reg284=reg59*reg139;
    T reg285=reg131*reg145; T reg286=reg90*reg163; T reg287=reg59*reg148; reg187=reg110*reg187; T reg288=reg107*reg168;
    T reg289=reg133*reg152; T reg290=reg133*reg136; reg154=reg110*reg154; T reg291=reg133*reg140; reg179=reg110*reg179;
    T reg292=reg107*reg169; T reg293=reg131*reg153; T reg294=reg131*reg155; reg176=reg110*reg176; T reg295=reg133*reg155;
    reg144=reg110*reg144; T reg296=reg107*reg189; T reg297=reg133*reg153; reg161=reg110*reg161; T reg298=reg133*reg145;
    reg183=reg110*reg183; T reg299=reg133*reg139; reg182=reg110*reg182; T reg300=reg133*reg148; T reg301=reg107*reg173;
    T reg302=reg131*reg152; T reg303=reg107*reg186; T reg304=reg131*reg140; T reg305=reg59*reg136; reg169=reg90*reg169;
    reg155=reg59*reg155; reg189=reg90*reg189; reg164=reg107*reg164; reg152=reg59*reg152; reg148=reg131*reg148;
    reg139=reg131*reg139; reg163=reg107*reg163; reg168=reg90*reg168; reg153=reg59*reg153; reg136=reg131*reg136;
    reg173=reg90*reg173; reg186=reg90*reg186; reg145=reg59*reg145; reg140=reg59*reg140; reg166=reg90*reg166;
    reg305=reg189+reg305; reg297=reg144-reg297; reg271=reg270-reg271; reg269=reg268-reg269; reg267=reg266-reg267;
    reg298=reg161-reg298; reg205=reg208-reg205; reg140=reg173+reg140; reg265=reg264-reg265; reg198=reg202-reg198;
    reg263=reg262-reg263; reg300=reg182-reg300; reg246=reg257+reg246; reg299=reg183-reg299; reg261=reg260-reg261;
    reg194=reg196-reg194; reg259=reg159-reg259; reg247=reg258+reg247; reg145=reg166+reg145; reg245=reg36*reg245;
    reg249=reg248+reg249; reg284=reg283+reg284; reg244=reg36*reg244; reg144=reg242+reg243; reg280=reg250+reg280;
    reg159=reg240+reg241; reg287=reg286+reg287; reg239=reg36*reg239; reg153=reg168+reg153; reg238=reg36*reg238;
    reg289=reg187-reg289; reg161=reg236+reg237; reg252=reg251+reg252; reg253=reg256+reg253; reg291=reg154-reg291;
    reg155=reg169+reg155; reg154=reg276+reg278; reg290=reg179-reg290; reg279=reg281+reg279; reg152=reg186+reg152;
    reg295=reg176-reg295; reg255=reg254+reg255; reg274=reg275+reg274; reg273=reg272-reg273; reg301=reg304-reg301;
    reg219=reg217+reg219; reg220=reg218+reg220; reg209=reg206-reg209; reg222=reg221+reg222; reg235=reg277-reg235;
    reg282=reg285-reg282; reg224=reg223+reg224; reg203=reg200+reg203; reg226=reg225+reg226; reg288=reg293-reg288;
    reg228=reg227+reg228; reg197=reg188+reg197; reg296=reg136-reg296; reg230=reg229+reg230; reg185=reg170+reg185;
    reg292=reg294-reg292; reg232=reg231+reg232; reg234=reg233+reg234; reg210=reg199-reg210; reg212=reg211+reg212;
    reg303=reg302-reg303; reg163=reg148-reg163; reg215=reg214+reg215; reg195=reg201-reg195; reg216=reg213+reg216;
    reg204=reg207-reg204; reg164=reg139-reg164; reg288=reg36*reg288; reg136=reg36*reg159; reg164=reg36*reg164;
    reg247=reg36*reg247; reg287=reg36*reg287; reg239=ponderation*reg239; reg140=reg36*reg140; reg212=reg36*reg212;
    reg228=reg36*reg228; reg238=ponderation*reg238; reg246=reg36*reg246; reg230=reg36*reg230; reg289=reg36*reg289;
    reg139=reg36*reg161; reg155=reg36*reg155; reg292=reg36*reg292; reg163=reg36*reg163; reg232=reg36*reg232;
    reg291=reg36*reg291; reg148=reg36*reg154; reg253=reg36*reg253; reg220=reg36*reg220; reg153=reg36*reg153;
    reg282=reg36*reg282; reg152=reg36*reg152; reg249=reg36*reg249; reg216=reg36*reg216; reg222=reg36*reg222;
    reg280=reg36*reg280; reg255=reg36*reg255; reg305=reg36*reg305; reg145=reg36*reg145; reg215=reg36*reg215;
    reg245=ponderation*reg245; reg224=reg36*reg224; reg284=reg36*reg284; reg244=ponderation*reg244; reg252=reg36*reg252;
    reg219=reg36*reg219; reg274=reg36*reg274; reg166=reg36*reg144; reg226=reg36*reg226; reg197=reg36*reg197;
    reg271=reg36*reg271; reg297=reg36*reg297; reg269=reg36*reg269; reg235=reg36*reg235; reg267=reg36*reg267;
    reg203=reg36*reg203; reg298=reg36*reg298; reg265=reg36*reg265; reg263=reg36*reg263; reg209=reg36*reg209;
    reg261=reg36*reg261; reg301=reg36*reg301; reg299=reg36*reg299; reg259=reg36*reg259; reg204=reg36*reg204;
    reg194=reg36*reg194; reg300=reg36*reg300; reg198=reg36*reg198; reg195=reg36*reg195; reg205=reg36*reg205;
    reg303=reg36*reg303; reg210=reg36*reg210; reg185=reg36*reg185; reg290=reg36*reg290; reg279=reg36*reg279;
    reg295=reg36*reg295; reg273=reg36*reg273; reg234=reg36*reg234; reg296=reg36*reg296; T tmp_3_5=ponderation*reg282;
    T tmp_5_2=ponderation*reg230; T tmp_4_3=ponderation*reg219; T tmp_3_2=ponderation*reg296; T tmp_0_2=ponderation*reg252; T tmp_2_6=ponderation*reg299;
    T tmp_6_7=ponderation*reg259; T tmp_1_3=ponderation*reg155; T tmp_2_2=ponderation*reg290; T tmp_0_0=ponderation*reg253; T tmp_6_6=ponderation*reg194;
    T tmp_4_2=ponderation*reg216; T tmp_0_6=ponderation*reg279; T tmp_6_1=ponderation*reg204; T tmp_3_0=ponderation*reg303; T tmp_3_7=ponderation*reg163;
    T tmp_4_6=ponderation*reg246; T tmp_6_3=ponderation*reg210; T tmp_5_4=ponderation*reg234; T tmp_6_2=ponderation*reg195; T tmp_0_4=ponderation*reg247;
    T tmp_6_4=ponderation*reg205; T tmp_4_0=ponderation*reg212; T tmp_1_1=ponderation*reg140; reg140=ponderation*reg148; sollicitation[indices[0]+0]+=reg140;
    T tmp_2_7=ponderation*reg300; T tmp_5_3=ponderation*reg232; T tmp_4_1=ponderation*reg215; T tmp_3_6=ponderation*reg164; T tmp_6_5=ponderation*reg198;
    T tmp_1_2=ponderation*reg305; T tmp_1_0=ponderation*reg152; T tmp_0_3=ponderation*reg255; T tmp_2_1=ponderation*reg291; T tmp_5_5=ponderation*reg185;
    sollicitation[indices[3]+0]+=-reg244; T tmp_7_4=ponderation*reg267; T tmp_7_7=ponderation*reg273; T tmp_1_6=ponderation*reg284; T tmp_0_7=ponderation*reg274;
    T tmp_7_2=ponderation*reg235; reg152=ponderation*reg166; sollicitation[indices[2]+1]+=reg152; sollicitation[indices[1]+0]+=-reg238; T tmp_7_5=ponderation*reg269;
    T tmp_5_0=ponderation*reg226; T tmp_2_4=ponderation*reg297; reg155=ponderation*reg136; sollicitation[indices[2]+0]+=reg155; T tmp_3_3=ponderation*reg292;
    T tmp_5_6=ponderation*reg197; T tmp_5_1=ponderation*reg228; T tmp_7_6=ponderation*reg271; sollicitation[indices[1]+1]+=-reg239; T tmp_1_7=ponderation*reg287;
    T tmp_7_0=ponderation*reg261; T tmp_6_0=ponderation*reg209; T tmp_1_4=ponderation*reg153; T tmp_4_4=ponderation*reg220; T tmp_0_1=ponderation*reg249;
    reg153=ponderation*reg139; sollicitation[indices[0]+1]+=reg153; T tmp_7_1=ponderation*reg263; T tmp_3_1=ponderation*reg301; T tmp_4_5=ponderation*reg222;
    T tmp_2_5=ponderation*reg298; T tmp_0_5=ponderation*reg280; T tmp_2_0=ponderation*reg289; T tmp_7_3=ponderation*reg265; T tmp_1_5=ponderation*reg145;
    T tmp_2_3=ponderation*reg295; sollicitation[indices[3]+1]+=-reg245; T tmp_5_7=ponderation*reg203; T tmp_4_7=ponderation*reg224; T tmp_3_4=ponderation*reg288;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; reg1=1.0/reg1; T reg2=2*(*f.m).shear_modulus_12;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=pow((*f.m).v2[0],2); T reg5=pow((*f.m).v2[1],2); T reg6=reg2*reg3;
    T reg7=pow((*f.m).v1[1],2); T reg8=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg11=pow((*f.m).v1[0],2);
    reg7=reg11+reg7; reg11=pow((*f.m).v1[2],2); T reg12=reg8*reg6; T reg13=reg10*reg6; T reg14=reg9*reg6;
    reg5=reg4+reg5; reg4=pow((*f.m).v2[2],2); T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; reg4=reg5+reg4;
    reg11=reg7+reg11; reg5=reg14*reg15; reg7=reg8*reg12; T reg17=reg14*reg16; T reg18=reg8*reg13;
    reg7=reg5-reg7; reg5=reg15*reg13; reg4=pow(reg4,0.5); reg11=pow(reg11,0.5); T reg19=1.0/(*f.m).elastic_modulus_1;
    reg18=reg17+reg18; T reg20=reg16*reg12; T reg21=(*f.m).v1[1]/reg11; T reg22=(*f.m).v1[2]/reg11; T reg23=reg19*reg7;
    T reg24=reg16*reg18; T reg25=(*f.m).v2[2]/reg4; T reg26=reg20+reg5; T reg27=(*f.m).v2[1]/reg4; T reg28=reg22*reg27;
    T reg29=reg10*reg3; T reg30=reg8*reg3; reg4=(*f.m).v2[0]/reg4; T reg31=reg21*reg25; reg11=(*f.m).v1[0]/reg11;
    T reg32=reg2*reg0; reg3=reg9*reg3; T reg33=reg16*reg6; T reg34=reg10*reg13; reg14=reg14*reg19;
    T reg35=reg10*reg12; reg6=reg15*reg6; T reg36=reg10*reg26; reg24=reg23-reg24; reg23=reg11*reg25;
    T reg37=reg31-reg28; T reg38=reg9*reg32; T reg39=2*reg11; T reg40=reg22*reg4; T reg41=2*reg4;
    T reg42=reg8*reg29; T reg43=reg8*reg30; T reg44=reg16*reg3; T reg45=reg10*reg32; reg36=reg24-reg36;
    reg24=reg10*reg6; reg35=reg17+reg35; reg34=reg14-reg34; reg12=reg19*reg12; reg3=reg15*reg3;
    reg14=reg10*reg33; reg32=reg8*reg32; reg13=reg16*reg13; reg17=reg2*reg1; T reg46=reg11*reg27;
    T reg47=reg8*reg17; reg30=reg16*reg30; reg33=reg16*reg33; reg13=reg12+reg13; reg14=reg12+reg14;
    reg34=reg34/reg36; reg18=reg18/reg36; reg24=reg20+reg24; reg35=reg35/reg36; reg7=reg7/reg36;
    reg6=reg19*reg6; reg29=reg15*reg29; reg12=reg10*reg17; reg43=reg3-reg43; reg3=2*reg37;
    T reg48=reg27*reg41; T reg49=pow(reg27,2); T reg50=pow(reg4,2); T reg51=reg21*reg39; T reg52=pow(reg21,2);
    T reg53=pow(reg11,2); T reg54=reg40-reg23; T reg55=reg21*reg4; reg42=reg44+reg42; reg44=reg15*reg38;
    reg38=reg16*reg38; T reg56=reg8*reg32; reg17=reg9*reg17; T reg57=reg8*reg45; reg24=reg24/reg36;
    T reg58=reg54*reg3; T reg59=reg8*reg12; reg42=reg16*reg42; T reg60=reg8*reg47; T reg61=reg16*reg17;
    reg17=reg15*reg17; reg57=reg38+reg57; reg43=reg19*reg43; reg32=reg16*reg32; reg38=reg53*reg7;
    T reg62=reg48*reg34; T reg63=reg51*reg35; T reg64=reg50*reg18; reg45=reg15*reg45; T reg65=reg53*reg35;
    T reg66=reg50*reg34; T reg67=reg52*reg35; T reg68=reg49*reg34; T reg69=reg46-reg55; reg33=reg6-reg33;
    reg6=reg30+reg29; T reg70=pow(reg22,2); reg13=reg13/reg36; T reg71=reg48*reg18; reg14=reg14/reg36;
    T reg72=reg51*reg7; reg56=reg44-reg56; reg44=pow(reg25,2); T reg73=reg52*reg7; T reg74=pow(reg54,2);
    reg26=reg26/reg36; T reg75=pow(reg37,2); T reg76=reg49*reg18; reg68=reg67+reg68; reg67=reg75*reg13;
    reg47=reg16*reg47; reg66=reg65+reg66; reg33=reg33/reg36; reg65=reg58*reg26; reg59=reg61+reg59;
    reg71=reg72+reg71; reg76=reg73+reg76; reg60=reg17-reg60; reg17=reg74*reg26; reg57=reg16*reg57;
    reg61=reg44*reg18; reg64=reg38+reg64; reg38=reg75*reg26; reg72=reg51*reg24; reg73=reg48*reg14;
    T reg77=reg32+reg45; T reg78=reg70*reg7; T reg79=reg49*reg14; reg56=reg19*reg56; T reg80=reg53*reg24;
    T reg81=reg44*reg34; reg6=reg10*reg6; T reg82=reg50*reg14; reg62=reg63+reg62; reg42=reg43-reg42;
    reg43=reg70*reg35; reg63=reg58*reg13; T reg83=reg52*reg24; reg12=reg15*reg12; T reg84=pow(reg69,2);
    T reg85=reg74*reg13; reg61=reg78+reg61; reg78=reg84*reg26; reg63=reg62+reg63; reg62=2*reg21;
    reg38=reg64+reg38; reg82=reg80+reg82; reg77=reg10*reg77; reg17=reg76+reg17; reg73=reg72+reg73;
    reg64=reg58*reg33; reg72=reg75*reg33; reg76=reg25*reg41; reg80=2*reg27; reg60=reg19*reg60;
    T reg86=reg47+reg12; reg59=reg16*reg59; reg85=reg68+reg85; reg67=reg66+reg67; reg65=reg71+reg65;
    reg66=reg11*reg21; reg68=reg4*reg27; reg6=reg42-reg6; reg79=reg83+reg79; reg81=reg43+reg81;
    reg42=reg84*reg13; reg57=reg56-reg57; reg43=reg74*reg33; reg56=reg22*reg39; reg71=reg70*reg24;
    reg83=reg44*reg14; T reg87=reg76*reg18; reg3=reg69*reg3; reg72=reg82+reg72; reg86=reg10*reg86;
    reg82=reg52*reg38; T reg88=reg49*reg67; T reg89=reg25*reg80; T reg90=reg68*reg85; T reg91=reg22*reg62;
    T reg92=reg56*reg7; reg78=reg61+reg78; reg61=reg50*reg67; T reg93=reg53*reg38; T reg94=reg50*reg85;
    T reg95=reg53*reg17; T reg96=reg49*reg63; T reg97=reg52*reg17; T reg98=reg4*reg39; T reg99=reg27*reg62;
    T reg100=2*reg22; T reg101=reg54*reg37; T reg102=reg66*reg65; reg85=reg49*reg85; reg42=reg81+reg42;
    reg81=reg68*reg63; reg43=reg79+reg43; reg83=reg71+reg83; reg71=reg53*reg65; reg63=reg50*reg63;
    reg79=reg56*reg35; T reg103=2*reg54; reg59=reg60-reg59; reg65=reg52*reg65; reg60=reg11*reg4;
    T reg104=reg21*reg27; reg46=reg55+reg46; reg77=reg57-reg77; reg64=reg73+reg64; reg17=reg66*reg17;
    reg6=reg6/reg36; reg55=reg76*reg34; reg57=reg21*reg37; reg73=reg11*reg54; T reg105=reg84*reg33;
    reg77=reg77/reg36; reg96=reg65+reg96; reg86=reg59-reg86; reg59=reg3*reg26; reg87=reg92+reg87;
    reg85=reg97+reg85; reg81=reg102+reg81; reg65=reg101*reg64; reg63=reg71+reg63; reg71=reg75*reg64;
    reg105=reg83+reg105; reg83=reg50*reg19; reg92=reg49*reg16; reg97=reg98*reg19; reg102=reg99*reg16;
    T reg106=reg56*reg24; T reg107=reg76*reg14; T reg108=reg50*reg16; T reg109=reg60*reg6; T reg110=reg49*reg15;
    T reg111=reg104*reg6; T reg112=reg98*reg16; T reg113=reg46*reg6; T reg114=reg99*reg15; T reg115=reg49*reg42;
    T reg116=reg52*reg78; reg61=reg93+reg61; reg93=reg75*reg72; reg94=reg95+reg94; reg95=reg74*reg72;
    reg88=reg82+reg88; reg82=reg75*reg43; T reg117=reg50*reg42; T reg118=reg53*reg78; reg100=reg25*reg100;
    reg38=reg66*reg38; reg67=reg68*reg67; T reg119=reg4*reg54; T reg120=reg27*reg37; reg73=reg57+reg73;
    reg57=reg21*reg54; reg18=reg89*reg18; reg7=reg91*reg7; T reg121=reg11*reg37; reg55=reg79+reg55;
    reg79=reg22*reg25; reg103=reg69*reg103; reg64=reg74*reg64; reg90=reg17+reg90; reg17=reg101*reg43;
    T reg122=reg3*reg13; reg35=reg91*reg35; reg34=reg89*reg34; reg43=reg74*reg43; reg82=reg94+reg82;
    reg26=reg103*reg26; reg18=reg7+reg18; reg93=reg61+reg93; reg59=reg87+reg59; reg7=reg98*reg109;
    reg61=reg98*reg111; reg115=reg116+reg115; reg87=reg74*reg105; reg94=reg73*reg77; reg116=reg57*reg77;
    T reg123=reg121*reg77; T reg124=reg79*reg6; reg14=reg89*reg14; reg24=reg91*reg24; T reg125=reg3*reg33;
    reg107=reg106+reg107; reg122=reg55+reg122; reg34=reg35+reg34; reg13=reg103*reg13; reg35=reg4*reg37;
    reg55=reg100*reg10; reg106=reg99*reg8; T reg126=reg98*reg10; reg92=reg83-reg92; reg83=reg44*reg10;
    T reg127=reg22*reg69; T reg128=reg49*reg8; reg17=reg90+reg17; reg90=reg46*reg111; T reg129=reg98*reg113;
    T reg130=reg50*reg10; reg112=reg114-reg112; reg71=reg63+reg71; reg78=reg66*reg78; reg42=reg68*reg42;
    reg63=reg46*reg113; reg65=reg81+reg65; reg36=reg86/reg36; reg64=reg96+reg64; reg113=reg99*reg113;
    reg95=reg88+reg95; reg81=reg99*reg109; reg111=reg99*reg111; reg43=reg85+reg43; reg85=reg44*reg8;
    reg108=reg110-reg108; reg62=reg54*reg62; reg86=reg100*reg8; reg39=reg37*reg39; reg117=reg118+reg117;
    reg88=reg75*reg105; reg67=reg38+reg67; reg72=reg101*reg72; reg119=reg120+reg119; reg102=reg97-reg102;
    reg38=reg27*reg54; reg85=reg108-reg85; reg96=reg39*reg116; reg33=reg103*reg33; reg111=reg43+reg111;
    reg86=reg112-reg86; reg43=reg127*reg77; reg97=reg62*reg123; reg108=(*f.m).alpha_1*reg53; reg110=reg62*reg116;
    reg106=reg126+reg106; reg61=reg82+reg61; reg82=reg39*reg123; reg7=reg93+reg7; reg100=reg100*reg9;
    reg93=reg50*reg122; reg128=reg130+reg128; reg112=reg44*reg9; reg19=reg53*reg19; reg114=reg53*reg59;
    reg118=reg119*reg36; reg81=reg95+reg81; reg95=reg38*reg36; reg120=reg35*reg36; reg87=reg115+reg87;
    reg115=reg99*reg124; reg116=reg73*reg116; reg90=reg17+reg90; reg40=reg23+reg40; reg17=reg22*reg37;
    reg23=reg11*reg69; reg126=reg98*reg124; reg88=reg117+reg88; reg109=reg46*reg109; reg72=reg67+reg72;
    reg67=reg25*reg69; reg41=reg37*reg41; reg80=reg54*reg80; reg117=reg49*reg122; reg130=reg52*reg59;
    reg26=reg18+reg26; reg18=reg62*reg94; reg113=reg64+reg113; reg14=reg24+reg14; reg15=reg52*reg15;
    reg24=reg39*reg94; reg125=reg107+reg125; reg64=reg53*reg16; reg55=reg102-reg55; reg83=reg92-reg83;
    reg129=reg71+reg129; reg16=reg52*reg16; reg94=reg73*reg94; reg63=reg65+reg63; reg105=reg101*reg105;
    reg65=(*f.m).alpha_2*reg49; reg71=(*f.m).alpha_1*reg52; reg42=reg78+reg42; reg78=(*f.m).alpha_2*reg50; reg13=reg34+reg13;
    reg34=reg80*reg120; reg97=reg81+reg97; reg81=reg50*reg13; reg92=reg53*reg26; reg122=reg68*reg122;
    reg102=reg104*reg85; reg107=reg60*reg83; T reg131=reg70*reg10; reg16=reg19-reg16; reg19=reg70*reg8;
    reg64=reg15-reg64; reg10=reg53*reg10; reg8=reg52*reg8; reg128=reg112-reg128; reg106=reg100-reg106;
    reg15=reg53*reg83; reg100=reg52*reg85; reg112=reg53*reg55; T reg132=reg52*reg86; T reg133=reg50*reg83;
    T reg134=reg49*reg85; T reg135=reg50*reg55; T reg136=reg49*reg86; reg110=reg111+reg110; reg111=reg80*reg95;
    reg115=reg87+reg115; reg87=reg62*reg43; reg18=reg113+reg18; reg113=reg80*reg118; reg117=reg130+reg117;
    reg130=reg74*reg125; T reg137=reg52*reg26; T reg138=reg49*reg13; reg109=reg72+reg109; reg123=reg73*reg123;
    reg116=reg90+reg116; reg72=reg119*reg95; reg105=reg42+reg105; reg124=reg46*reg124; reg94=reg63+reg94;
    reg42=reg119*reg118; reg59=reg66*reg59; reg93=reg114+reg93; reg118=reg41*reg118; reg24=reg129+reg24;
    reg23=reg17+reg23; reg17=reg22*reg54; reg63=reg21*reg69; reg90=reg39*reg43; reg126=reg88+reg126;
    reg95=reg41*reg95; reg96=reg61+reg96; reg61=reg25*reg37; reg88=reg41*reg120; reg82=reg7+reg82;
    reg7=reg4*reg69; reg33=reg14+reg33; reg14=reg67*reg36; reg114=reg40*reg6; reg129=reg75*reg125;
    T reg139=(*f.m).alpha_2*reg44; reg31=reg28+reg31; reg28=(*f.m).alpha_1*reg70; T reg140=reg60*reg55; T reg141=reg104*reg86;
    T reg142=(*f.m).alpha_3*reg74; reg65=reg71+reg65; reg78=reg108+reg78; reg71=(*f.m).alpha_3*reg75; reg4=reg4*reg25;
    reg108=reg46*reg2; reg11=reg11*reg22; T reg143=(*f.m).alpha_3*reg84; T reg144=reg74*reg33; reg130=reg117+reg130;
    reg117=reg99*reg114; reg42=reg94+reg42; reg138=reg137+reg138; reg131=reg16-reg131; reg16=(*f.m).alpha_2*reg68;
    reg94=reg27*reg69; reg137=(*f.m).alpha_1*reg66; reg124=reg105+reg124; reg142=reg65+reg142; reg71=reg78+reg71;
    reg72=reg116+reg72; reg65=reg44*reg106; reg43=reg73*reg43; reg134=reg133+reg134; reg78=reg79*reg106;
    reg141=reg140+reg141; reg120=reg119*reg120; reg123=reg109+reg123; reg63=reg17+reg63; reg7=reg61+reg7;
    reg17=reg25*reg54; reg61=reg79*reg128; reg105=reg44*reg128; reg139=reg28+reg139; reg28=reg68*reg2;
    reg8=reg10+reg8; reg10=reg98*reg114; reg129=reg93+reg129; reg118=reg24+reg118; reg34=reg97+reg34;
    reg24=reg41*reg14; reg90=reg126+reg90; reg95=reg96+reg95; reg88=reg82+reg88; reg111=reg110+reg111;
    reg13=reg68*reg13; reg26=reg66*reg26; reg102=reg107+reg102; reg9=reg70*reg9; reg68=reg23*reg77;
    reg6=reg31*reg6; reg19=reg64-reg19; reg113=reg18+reg113; reg122=reg59+reg122; reg125=reg101*reg125;
    reg81=reg92+reg81; reg18=reg75*reg33; reg59=reg70*reg106; reg132=reg112+reg132; reg64=reg70*reg128;
    reg100=reg15+reg100; reg136=reg135+reg136; reg15=reg80*reg14; reg87=reg115+reg87; reg82=reg88*reg71;
    reg8=reg9-reg8; reg43=reg124+reg43; reg14=reg119*reg14; reg9=reg51*reg108; reg92=reg95*reg142;
    reg59=reg132+reg59; reg33=reg101*reg33; reg13=reg26+reg13; reg26=reg48*reg108; reg143=reg139+reg143;
    reg65=reg136+reg65; reg93=reg34*reg71; reg96=reg111*reg142; reg114=reg46*reg114; reg16=reg137+reg16;
    reg101=(*f.m).alpha_3*reg101; reg125=reg122+reg125; reg97=(*f.m).alpha_2*reg4; reg107=(*f.m).alpha_1*reg11; reg109=reg40*reg1;
    reg64=reg100+reg64; reg100=reg113*reg72; reg110=reg50*reg131; reg112=reg53*reg131; reg15=reg87+reg15;
    reg77=reg63*reg77; reg87=reg7*reg36; reg115=reg52*reg19; reg116=reg46*reg28; reg24=reg90+reg24;
    reg10=reg129+reg10; reg90=reg39*reg68; reg61=reg102+reg61; reg98=reg98*reg6; reg18=reg81+reg18;
    reg81=reg118*reg72; reg102=reg49*reg19; reg122=reg111*reg42; reg124=reg95*reg42; reg126=reg46*reg108;
    reg78=reg141+reg78; reg120=reg123+reg120; reg105=reg134+reg105; reg99=reg99*reg6; reg144=reg138+reg144;
    reg94=reg17+reg94; reg17=reg48*reg28; reg123=reg51*reg28; reg4=reg4*reg1; reg129=reg62*reg68;
    reg117=reg130+reg117; reg37=reg69*reg37; reg22=reg21*reg22; reg25=reg27*reg25; reg9=reg59+reg9;
    reg21=reg56*reg109; reg27=reg118*reg111; reg92=reg82+reg92; reg24=reg24*reg143; reg59=reg72*reg142;
    reg96=reg93+reg96; reg15=reg15*reg143; reg82=reg56*reg4; reg123=reg64+reg123; reg2=reg66*reg2;
    reg64=reg25*reg0; reg66=reg31*reg0; reg115=reg112+reg115; reg70=reg70*reg8; reg93=reg120*reg71;
    reg112=reg60*reg131; reg130=reg104*reg19; reg33=reg13+reg33; reg36=reg94*reg36; reg68=reg73*reg68;
    reg114=reg125+reg114; reg101=reg16+reg101; reg116=reg61+reg116; reg13=reg40*reg4; reg14=reg43+reg14;
    reg16=reg40*reg109; reg126=reg78+reg126; reg62=reg62*reg77; reg99=reg144+reg99; reg43=reg80*reg87;
    reg129=reg117+reg129; reg54=reg69*reg54; reg17=reg105+reg17; reg61=reg76*reg4; reg69=reg41*reg87;
    reg81=reg124-reg81; reg90=reg10+reg90; reg100=reg122-reg100; reg102=reg110+reg102; reg44=reg44*reg8;
    reg25=(*f.m).alpha_2*reg25; reg10=reg95*reg113; reg78=(*f.m).alpha_1*reg22; reg98=reg18+reg98; reg39=reg39*reg77;
    reg37=(*f.m).alpha_3*reg37; reg97=reg107+reg97; reg18=reg76*reg109; reg26=reg65+reg26; reg6=reg46*reg6;
    reg61=reg17+reg61; reg17=reg89*reg64; reg14=reg14*reg143; reg82=reg123+reg82; reg65=reg91*reg64;
    reg105=reg89*reg66; reg59=reg93+reg59; reg93=reg48*reg2; reg44=reg102+reg44; reg18=reg26+reg18;
    reg21=reg9+reg21; reg9=reg91*reg66; reg27=reg10-reg27; reg10=reg34*reg81; reg26=reg88*reg100;
    reg54=(*f.m).alpha_3*reg54; reg25=reg78+reg25; reg37=reg97+reg37; reg77=reg73*reg77; reg6=reg33+reg6;
    reg87=reg119*reg87; reg68=reg114+reg68; reg33=reg31*reg66; reg16=reg126+reg16; reg80=reg80*reg36;
    reg62=reg99+reg62; reg43=reg129+reg43; reg78=reg31*reg64; reg13=reg116+reg13; reg41=reg41*reg36;
    reg39=reg98+reg39; reg69=reg90+reg69; reg90=reg113*reg101; reg1=reg11*reg1; reg15=reg96+reg15;
    reg70=reg115+reg70; reg11=reg118*reg101; reg24=reg92+reg24; reg92=reg51*reg2; reg79=reg79*reg8;
    reg130=reg112+reg130; reg78=reg13+reg78; reg90=reg15+reg90; reg79=reg130+reg79; reg43=reg43*reg37;
    reg105=reg18+reg105; reg13=reg42*reg101; reg14=reg59+reg14; reg33=reg16+reg33; reg69=reg69*reg37;
    reg15=reg46*reg2; reg11=reg24+reg11; reg54=reg25+reg54; reg16=reg118*reg120; reg18=reg88*reg42;
    reg10=reg26-reg10; reg24=reg120*reg27; reg25=reg113*reg120; reg26=reg34*reg42; reg41=reg39+reg41;
    reg80=reg62+reg80; reg87=reg68+reg87; reg77=reg6+reg77; reg36=reg119*reg36; reg0=reg22*reg0;
    reg92=reg70+reg92; reg56=reg56*reg1; reg6=1-var_inter[0]; reg65=reg82+reg65; reg9=reg21+reg9;
    reg93=reg44+reg93; reg76=reg76*reg1; reg21=1-var_inter[1]; reg17=reg61+reg17; reg22=reg21*elem.pos(1)[0];
    reg39=reg88*reg72; reg44=reg111*reg120; reg25=reg26-reg25; reg26=reg34*reg72; reg59=reg21*elem.pos(0)[1];
    reg61=reg21*elem.pos(1)[1]; reg24=reg10+reg24; reg10=elem.pos(0)[0]*reg6; reg62=elem.pos(1)[0]*var_inter[0]; reg68=reg6*elem.pos(0)[1];
    reg70=elem.pos(1)[1]*var_inter[0]; reg16=reg18-reg16; reg18=reg95*reg120; reg82=reg88*reg113; reg96=reg118*reg34;
    reg97=reg21*elem.pos(0)[0]; reg69=reg11+reg69; reg41=reg41*reg54; reg87=reg87*reg37; reg13=reg14+reg13;
    reg43=reg90+reg43; reg80=reg80*reg54; reg11=reg105*reg78; reg91=reg91*reg0; reg76=reg93+reg76;
    reg89=reg89*reg0; reg56=reg92+reg56; reg36=reg77+reg36; reg15=reg79+reg15; reg14=reg9*reg78;
    reg40=reg40*reg1; reg77=reg17*reg33; reg79=reg65*reg33; reg81=reg81/reg24; reg90=reg9*reg17;
    reg44=reg26-reg44; reg22=reg22-reg97; reg26=reg65*reg105; reg25=reg25/reg24; reg100=reg100/reg24;
    reg92=elem.pos(2)[1]*var_inter[1]; reg91=reg56+reg91; reg61=reg61-reg59; reg56=elem.pos(2)[0]*var_inter[0]; reg89=reg76+reg89;
    reg76=reg10+reg62; reg87=reg13+reg87; reg36=reg36*reg54; reg31=reg31*reg0; reg40=reg15+reg40;
    reg41=reg69+reg41; reg80=reg43+reg80; reg13=elem.pos(2)[0]*var_inter[1]; reg15=reg95*reg34; reg96=reg82-reg96;
    reg43=reg88*reg111; reg11=reg77-reg11; reg14=reg79-reg14; reg18=reg39-reg18; reg16=reg16/reg24;
    reg39=reg68+reg70; reg69=elem.pos(2)[1]*var_inter[0]; reg31=reg40+reg31; reg56=reg56-reg76; reg18=reg18/reg24;
    reg40=reg11*reg91; reg77=elem.pos(3)[1]*var_inter[1]; reg15=reg43-reg15; reg61=reg92+reg61; reg44=reg44/reg24;
    reg43=reg89*reg14; reg27=reg27/reg24; reg79=elem.pos(3)[0]*var_inter[1]; reg96=reg96/reg24; reg22=reg13+reg22;
    reg13=elem.pos(3)[0]*reg6; reg36=reg87+reg36; reg90=reg26-reg90; reg25=reg25*reg41; reg81=reg81*reg80;
    reg69=reg69-reg39; reg16=reg16*reg80; reg100=reg100*reg41; reg26=elem.pos(3)[1]*reg6; reg22=reg22-reg79;
    reg82=reg17*reg31; reg26=reg69+reg26; reg69=reg89*reg78; reg13=reg56+reg13; reg56=reg31*reg90;
    reg43=reg40-reg43; reg27=reg27*reg36; reg81=reg100-reg81; reg40=reg91*reg78; reg96=reg96*reg36;
    reg24=reg15/reg24; reg61=reg61-reg77; reg15=reg65*reg31; reg25=reg16-reg25; reg41=reg44*reg41;
    reg80=reg18*reg80; reg81=reg27+reg81; reg16=reg22*reg26; reg96=reg25-reg96; reg36=reg24*reg36;
    reg80=reg41-reg80; reg18=reg61*reg13; reg56=reg43+reg56; reg24=reg65*reg89; reg25=reg89*reg33;
    reg27=1-(*f.m).resolution; reg41=reg105*reg31; reg43=reg91*reg17; reg44=reg91*reg33; reg15=reg40-reg15;
    reg82=reg69-reg82; reg40=reg9*reg31; reg69=reg9*reg89; reg87=reg91*reg105; reg24=reg43-reg24;
    reg43=(*f.m).resolution*reg142; reg80=reg36+reg80; reg15=reg15/reg56; reg40=reg44-reg40; reg36=(*f.m).resolution*reg71;
    reg96=reg27*reg96; reg82=reg82/reg56; reg41=reg25-reg41; reg18=reg16-reg18; reg81=reg27*reg81;
    reg81=reg36+reg81; reg80=reg27*reg80; reg96=reg43+reg96; reg16=(*f.m).resolution*reg101; reg120=reg120*reg27;
    reg72=reg72*reg27; reg25=(*f.m).resolution*reg82; reg36=(*f.m).resolution*reg15; reg61=reg61/reg18; reg14=reg14/reg56;
    reg11=reg11/reg56; reg41=reg41/reg56; reg13=reg13/reg18; reg24=reg24/reg56; reg69=reg87-reg69;
    reg40=reg40/reg56; reg22=reg22/reg18; reg26=reg26/reg18; reg69=reg69/reg56; reg88=reg88*reg27;
    reg95=reg95*reg27; reg43=var_inter[0]*reg22; reg34=reg34*reg27; reg111=reg111*reg27; reg56=reg90/reg56;
    reg42=reg42*reg27; reg44=var_inter[0]*reg61; reg36=reg72-reg36; reg120=reg25+reg120; reg25=reg6*reg61;
    reg96=(*f.m).deltaT*reg96; reg72=reg6*reg22; reg87=reg21*reg13; reg90=(*f.m).resolution*reg41; reg81=(*f.m).deltaT*reg81;
    reg92=(*f.m).resolution*reg40; reg80=reg16+reg80; reg16=(*f.m).resolution*reg11; reg93=reg21*reg26; reg98=(*f.m).resolution*reg24;
    reg99=(*f.m).resolution*reg14; reg100=var_inter[1]*reg13; reg102=var_inter[1]*reg26; reg107=reg44+reg93; reg110=reg25+reg102;
    reg112=reg43-reg100; reg114=reg102-reg44; reg115=reg87+reg43; reg80=(*f.m).deltaT*reg80; reg42=reg98+reg42;
    reg98=reg36*reg96; reg116=reg120*reg81; reg117=(*f.m).resolution*reg56; reg111=reg92+reg111; reg90=reg34-reg90;
    reg99=reg95-reg99; reg34=(*f.m).resolution*reg69; reg88=reg16+reg88; reg113=reg113*reg27; reg27=reg118*reg27;
    reg16=reg87-reg72; reg92=reg25-reg93; reg95=reg72+reg100; reg118=0.5*reg16; reg122=0.5*reg110;
    reg123=0.5*reg95; reg124=reg116+reg98; reg125=reg42*reg80; reg126=reg111*reg96; reg129=reg90*reg81;
    reg130=reg88*reg81; reg132=reg99*reg96; reg34=reg113-reg34; reg113=0.5*reg114; reg133=0.5*reg112;
    reg134=0.5*reg107; reg27=reg117+reg27; reg117=0.5*reg92; reg135=0.5*reg115; reg136=reg135*reg42;
    reg137=reg107*reg120; reg138=reg110*reg120; reg139=reg123*reg42; reg140=reg124+reg125; reg141=reg117*reg42;
    reg144=reg34*reg80; T reg145=reg129+reg126; T reg146=reg112*reg36; T reg147=reg27*reg80; T reg148=reg95*reg36;
    T reg149=reg122*reg42; T reg150=reg130+reg132; T reg151=reg113*reg42; T reg152=reg16*reg36; T reg153=reg115*reg36;
    T reg154=reg134*reg42; T reg155=reg92*reg120; T reg156=reg118*reg42; T reg157=reg114*reg120; T reg158=reg133*reg42;
    T reg159=reg133*reg34; T reg160=reg134*reg34; T reg161=reg114*reg90; T reg162=reg112*reg111; T reg163=reg113*reg34;
    T reg164=reg110*reg90; T reg165=2*reg140; T reg166=reg123*reg34; T reg167=reg95*reg111; T reg168=reg122*reg34;
    T reg169=reg115*reg111; reg151=reg146+reg151; reg146=reg112*reg99; T reg170=reg133*reg27; T reg171=reg113*reg27;
    T reg172=reg114*reg88; reg158=reg157+reg158; reg157=reg134*reg27; T reg173=reg115*reg99; reg154=reg154-reg153;
    T reg174=reg135*reg27; T reg175=reg107*reg88; reg137=reg137-reg136; T reg176=reg117*reg27; T reg177=reg16*reg99;
    reg141=reg152+reg141; reg152=reg145+reg144; T reg178=reg150+reg147; T reg179=reg95*reg99; T reg180=reg117*reg34;
    T reg181=reg16*reg111; reg148=reg148-reg149; T reg182=reg122*reg27; reg156=reg155+reg156; reg155=reg107*reg90;
    T reg183=reg123*reg27; T reg184=reg110*reg88; T reg185=reg135*reg34; T reg186=reg118*reg27; T reg187=reg92*reg88;
    reg139=reg139-reg138; reg170=reg172+reg170; reg141=2*reg141; reg157=reg157-reg173; reg139=2*reg139;
    reg166=reg166-reg164; reg179=reg179-reg182; reg172=reg110*reg178; reg186=reg187+reg186; reg176=reg177+reg176;
    reg167=reg167-reg168; reg154=2*reg154; reg177=reg95*reg152; reg158=2*reg158; reg156=2*reg156;
    reg187=reg122*reg165; T reg188=reg135*reg165; reg160=reg160-reg169; T reg189=reg107*reg178; reg151=2*reg151;
    reg171=reg146+reg171; reg146=reg134*reg165; reg137=2*reg137; T reg190=reg115*reg152; reg159=reg161+reg159;
    reg183=reg183-reg184; reg175=reg175-reg174; reg161=reg123*reg165; reg180=reg181+reg180; reg163=reg162+reg163;
    reg148=2*reg148; reg155=reg155-reg185; reg162=reg112*reg163; reg181=reg113*reg151; T reg191=reg118*reg141;
    T reg192=reg115*reg163; T reg193=reg134*reg139; T reg194=reg133*reg148; T reg195=reg114*reg179; T reg196=reg115*reg166;
    T reg197=reg118*reg137; T reg198=reg92*reg176; T reg199=reg133*reg151; T reg200=reg134*reg148; T reg201=reg114*reg171;
    T reg202=reg115*reg167; T reg203=reg92*reg171; T reg204=reg114*reg170; T reg205=reg133*reg158; T reg206=reg114*reg178;
    T reg207=reg190-reg146; T reg208=reg133*reg165; T reg209=reg188-reg189; T reg210=reg112*reg152; T reg211=reg117*reg165;
    T reg212=reg16*reg152; T reg213=reg92*reg183; T reg214=reg113*reg165; T reg215=reg92*reg178; T reg216=reg118*reg148;
    T reg217=reg118*reg165; T reg218=reg92*reg179; T reg219=reg172-reg161; T reg220=reg122*reg148; T reg221=reg95*reg167;
    T reg222=reg110*reg179; T reg223=reg123*reg148; T reg224=reg187-reg177; T reg225=reg110*reg183; T reg226=reg123*reg139;
    T reg227=reg113*reg148; T reg228=reg112*reg167; T reg229=reg118*reg139; T reg230=reg113*reg139; T reg231=reg112*reg166;
    T reg232=reg135*reg151; T reg233=reg107*reg171; T reg234=reg92*reg157; T reg235=reg135*reg158; T reg236=reg107*reg170;
    T reg237=reg135*reg154; T reg238=reg107*reg157; T reg239=reg118*reg156; T reg240=reg107*reg175; T reg241=reg135*reg137;
    T reg242=reg117*reg148; reg167=reg16*reg167; T reg243=reg118*reg158; T reg244=reg117*reg139; T reg245=reg16*reg166;
    T reg246=reg92*reg170; T reg247=reg117*reg151; T reg248=reg16*reg163; T reg249=reg16*reg180; T reg250=reg117*reg158;
    T reg251=reg16*reg159; T reg252=reg117*reg141; T reg253=reg114*reg183; T reg254=reg133*reg139; T reg255=reg117*reg154;
    T reg256=reg16*reg160; T reg257=reg118*reg151; T reg258=reg117*reg137; T reg259=reg16*reg155; T reg260=reg134*reg151;
    T reg261=reg92*reg175; T reg262=reg115*reg159; T reg263=reg115*reg160; T reg264=reg134*reg158; T reg265=reg134*reg154;
    T reg266=reg92*reg186; reg148=reg135*reg148; reg179=reg107*reg179; reg183=reg107*reg183; T reg267=reg118*reg154;
    T reg268=reg135*reg139; reg268=reg183-reg268; reg247=reg248+reg247; reg202=reg200-reg202; reg194=reg195+reg194;
    reg181=reg162+reg181; reg246=reg243+reg246; reg162=reg217+reg215; reg183=reg210+reg214; reg250=reg251+reg250;
    reg213=reg229+reg213; reg263=reg265-reg263; reg196=reg193-reg196; reg255=reg256+reg255; reg193=reg212+reg211;
    reg262=reg264-reg262; reg261=reg197+reg261; reg252=reg249+reg252; reg254=reg253+reg254; reg209=reg18*reg209;
    reg195=reg206+reg208; reg192=reg260-reg192; reg258=reg259+reg258; reg232=reg233-reg232; reg207=reg18*reg207;
    reg230=reg231+reg230; reg234=reg267+reg234; reg237=reg238-reg237; reg227=reg228+reg227; reg224=reg18*reg224;
    reg198=reg191+reg198; reg199=reg201+reg199; reg225=reg226-reg225; reg241=reg240-reg241; reg222=reg223-reg222;
    reg242=reg167+reg242; reg257=reg203+reg257; reg266=reg239+reg266; reg244=reg245+reg244; reg218=reg216+reg218;
    reg235=reg236-reg235; reg219=reg18*reg219; reg220=reg221-reg220; reg148=reg179-reg148; reg205=reg204+reg205;
    reg198=reg18*reg198; reg257=reg18*reg257; reg167=reg18*reg195; reg246=reg18*reg246; reg179=reg18*reg183;
    reg261=reg18*reg261; reg224=ponderation*reg224; reg254=reg18*reg254; reg266=reg18*reg266; reg234=reg18*reg234;
    reg219=ponderation*reg219; reg227=reg18*reg227; reg244=reg18*reg244; reg230=reg18*reg230; reg181=reg18*reg181;
    reg242=reg18*reg242; reg194=reg18*reg194; reg241=reg18*reg241; reg199=reg18*reg199; reg205=reg18*reg205;
    reg237=reg18*reg237; reg202=reg18*reg202; reg196=reg18*reg196; reg235=reg18*reg235; reg192=reg18*reg192;
    reg232=reg18*reg232; reg262=reg18*reg262; reg268=reg18*reg268; reg263=reg18*reg263; reg148=reg18*reg148;
    reg207=ponderation*reg207; reg258=reg18*reg258; reg209=ponderation*reg209; reg191=reg18*reg193; reg252=reg18*reg252;
    reg197=reg18*reg162; reg213=reg18*reg213; reg255=reg18*reg255; reg218=reg18*reg218; reg220=reg18*reg220;
    reg225=reg18*reg225; reg250=reg18*reg250; reg247=reg18*reg247; reg222=reg18*reg222; T tmp_4_6=ponderation*reg254;
    T tmp_2_6=ponderation*reg268; T tmp_1_5=ponderation*reg247; T tmp_2_5=ponderation*reg232; T tmp_1_2=ponderation*reg258; T tmp_1_1=ponderation*reg252;
    T tmp_2_4=ponderation*reg235; T tmp_1_6=ponderation*reg244; T tmp_2_3=ponderation*reg237; T tmp_1_4=ponderation*reg250; T tmp_0_7=ponderation*reg218;
    T tmp_0_3=ponderation*reg234; T tmp_2_2=ponderation*reg241; T tmp_1_3=ponderation*reg255; T tmp_0_4=ponderation*reg246; T tmp_1_7=ponderation*reg242;
    sollicitation[indices[1]+1]+=-reg207; sollicitation[indices[1]+0]+=-reg209; reg200=ponderation*reg191; sollicitation[indices[0]+1]+=reg200; reg201=ponderation*reg167;
    sollicitation[indices[2]+0]+=reg201; reg203=ponderation*reg197; sollicitation[indices[0]+0]+=reg203; reg204=ponderation*reg179; sollicitation[indices[2]+1]+=reg204;
    T tmp_0_6=ponderation*reg213; T tmp_7_7=ponderation*reg220; T tmp_6_7=ponderation*reg222; sollicitation[indices[3]+0]+=-reg219; T tmp_6_6=ponderation*reg225;
    T tmp_5_7=ponderation*reg227; sollicitation[indices[3]+1]+=-reg224; T tmp_5_6=ponderation*reg230; T tmp_5_5=ponderation*reg181; T tmp_4_7=ponderation*reg194;
    T tmp_0_5=ponderation*reg257; T tmp_4_5=ponderation*reg199; T tmp_4_4=ponderation*reg205; T tmp_0_1=ponderation*reg198; T tmp_3_7=ponderation*reg202;
    T tmp_3_6=ponderation*reg196; T tmp_3_5=ponderation*reg192; T tmp_0_0=ponderation*reg266; T tmp_3_4=ponderation*reg262; T tmp_3_3=ponderation*reg263;
    T tmp_0_2=ponderation*reg261; T tmp_2_7=ponderation*reg148;
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
    T reg0=pow((*f.m).v1[1],2); T reg1=pow((*f.m).v1[0],2); T reg2=pow((*f.m).v2[0],2); T reg3=pow((*f.m).v2[1],2); T reg4=pow((*f.m).v1[2],2);
    reg0=reg1+reg0; reg1=2*(*f.m).shear_modulus_23; reg3=reg2+reg3; reg2=pow((*f.m).v2[2],2); reg4=reg0+reg4;
    reg0=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; reg2=reg3+reg2; reg4=pow(reg4,0.5); reg3=2*(*f.m).shear_modulus_12;
    reg1=1.0/reg1; reg3=1.0/reg3; T reg5=(*f.m).v1[0]/reg4; T reg6=(*f.m).v1[1]/reg4; reg2=pow(reg2,0.5);
    T reg7=reg0*reg1; reg4=(*f.m).v1[2]/reg4; T reg8=(*f.m).v2[0]/reg2; T reg9=(*f.m).v2[1]/reg2; T reg10=2*reg5;
    T reg11=2*reg6; T reg12=reg3*reg7; T reg13=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg14=1.0/(*f.m).elastic_modulus_3; T reg15=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg16=reg15*reg12; T reg17=reg13*reg12; T reg18=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg19=1.0/(*f.m).elastic_modulus_2; T reg20=1.0/(*f.m).elastic_modulus_1;
    T reg21=2*reg4; T reg22=reg9*reg11; T reg23=reg8*reg10; T reg24=reg14*reg12; T reg25=pow(reg9,2);
    T reg26=pow(reg8,2); reg2=(*f.m).v2[2]/reg2; T reg27=reg25*reg18; T reg28=reg23*reg18; T reg29=reg26*reg20;
    T reg30=reg23*reg20; T reg31=reg25*reg19; T reg32=reg13*reg16; T reg33=pow(reg2,2); T reg34=reg22*reg18;
    T reg35=reg24*reg19; T reg36=reg26*reg18; reg21=reg2*reg21; T reg37=reg13*reg17; T reg38=reg22*reg19;
    T reg39=reg24*reg18; T reg40=reg26*reg15; reg28=reg38-reg28; reg38=reg33*reg15; T reg41=reg25*reg13;
    T reg42=reg23*reg15; T reg43=reg18*reg17; reg27=reg29-reg27; reg29=reg19*reg16; T reg44=reg21*reg15;
    T reg45=pow(reg6,2); reg34=reg30-reg34; reg30=pow(reg5,2); reg37=reg35-reg37; reg32=reg39+reg32;
    reg35=reg33*reg13; T reg46=reg21*reg13; T reg47=reg22*reg13; reg36=reg31-reg36; reg31=reg20*reg37;
    T reg48=reg6*reg8; T reg49=reg5*reg9; reg44=reg34-reg44; reg34=reg45*reg18; T reg50=reg18*reg32;
    T reg51=reg30*reg18; T reg52=reg43+reg29; T reg53=pow(reg4,2); reg38=reg27-reg38; reg46=reg28-reg46;
    reg35=reg36-reg35; reg27=reg30*reg20; reg47=reg42+reg47; reg21=reg21*reg14; reg28=reg5*reg8;
    reg36=reg33*reg14; reg41=reg40+reg41; reg40=reg45*reg19; reg42=reg6*reg9; T reg54=reg8*reg9;
    reg50=reg31-reg50; reg31=reg53*reg15; T reg55=reg28*reg38; reg51=reg40-reg51; reg40=reg53*reg13;
    T reg56=reg45*reg13; T reg57=2*reg8; T reg58=reg15*reg52; T reg59=reg30*reg15; T reg60=reg42*reg35;
    T reg61=reg4*reg2; T reg62=reg48+reg49; T reg63=reg18*reg12; T reg64=reg15*reg16; T reg65=reg45*reg35;
    T reg66=reg30*reg38; reg24=reg24*reg20; reg34=reg27-reg34; reg27=reg15*reg17; reg47=reg21-reg47;
    reg21=reg30*reg44; T reg67=reg45*reg46; T reg68=reg14*reg7; T reg69=reg42*reg46; reg41=reg36-reg41;
    reg36=reg28*reg44; T reg70=reg3*reg1; T reg71=reg6*reg2; T reg72=reg26*reg38; T reg73=reg4*reg9;
    T reg74=reg25*reg35; reg12=reg19*reg12; T reg75=reg4*reg8; T reg76=reg13*reg7; T reg77=reg5*reg2;
    reg7=reg15*reg7; T reg78=reg26*reg44; T reg79=reg25*reg46; T reg80=reg14*reg70; reg58=reg50-reg58;
    reg50=reg8*reg2; reg31=reg34-reg31; reg60=reg55+reg60; reg34=reg61*reg47; reg69=reg36+reg69;
    reg36=reg61*reg41; reg55=reg18*reg68; T reg81=reg13*reg7; T reg82=reg13*reg76; T reg83=reg62*reg3;
    T reg84=reg54*reg3; T reg85=reg15*reg12; reg27=reg39+reg27; reg67=reg21+reg67; reg21=reg53*reg47;
    reg68=reg19*reg68; reg56=reg59+reg56; reg39=reg53*reg14; reg59=reg15*reg63; T reg86=reg13*reg70;
    reg16=reg18*reg16; reg74=reg72+reg74; reg72=reg33*reg41; reg65=reg66+reg65; reg79=reg78+reg79;
    reg66=reg33*reg47; reg78=reg53*reg41; reg40=reg51-reg40; reg51=reg77+reg75; T reg87=reg71-reg73;
    T reg88=reg9*reg57; T reg89=reg3*reg0; reg70=reg15*reg70; T reg90=reg6*reg10; reg17=reg20*reg17;
    reg64=reg24-reg64; reg24=reg13*reg86; reg85=reg43+reg85; reg32=reg32/reg58; T reg91=reg19*reg80;
    reg64=reg64/reg58; reg80=reg18*reg80; reg27=reg27/reg58; T reg92=reg13*reg70; reg81=reg55+reg81;
    reg14=reg14*reg89; reg59=reg17+reg59; reg82=reg68-reg82; reg16=reg17+reg16; reg63=reg18*reg63;
    reg76=reg18*reg76; reg17=reg13*reg89; reg7=reg19*reg7; reg89=reg15*reg89; reg55=reg45*reg40;
    reg68=reg62*reg83; reg78=reg65+reg78; reg65=reg90*reg84; reg34=reg69+reg34; reg21=reg67+reg21;
    reg67=reg90*reg83; reg69=reg26*reg31; T reg93=reg25*reg40; T reg94=reg62*reg84; T reg95=reg88*reg83;
    reg66=reg79+reg66; reg72=reg74+reg72; reg74=reg88*reg84; reg36=reg60+reg36; reg37=reg37/reg58;
    reg12=reg20*reg12; reg60=reg9*reg2; reg79=reg5*reg6; reg71=reg73+reg71; reg73=2*reg87;
    reg56=reg39-reg56; reg39=2*reg9; T reg96=reg2*reg57; T reg97=reg4*reg10; reg77=reg75-reg77;
    reg75=reg50*reg0; T reg98=reg51*reg0; T reg99=reg30*reg31; T reg100=reg13*reg89; T reg101=reg45*reg37;
    T reg102=reg88*reg32; T reg103=reg25*reg32; T reg104=reg90*reg37; T reg105=reg42*reg40; reg13=reg13*reg17;
    T reg106=reg18*reg14; reg14=reg19*reg14; reg94=reg36+reg94; reg36=reg51*reg75; reg92=reg80+reg92;
    reg24=reg91-reg24; reg68=reg34+reg68; reg34=reg51*reg98; reg3=reg79*reg3; reg80=reg60*reg1;
    reg91=reg71*reg1; reg55=reg99+reg55; reg99=reg53*reg56; reg65=reg78+reg65; reg78=reg97*reg75;
    reg67=reg21+reg67; reg21=reg97*reg98; reg93=reg69+reg93; reg69=reg33*reg56; T reg107=reg88*reg64;
    T reg108=reg90*reg27; reg74=reg72+reg74; reg72=reg25*reg64; T reg109=reg96*reg75; reg95=reg66+reg95;
    reg66=reg96*reg98; T reg110=reg28*reg31; T reg111=reg45*reg27; T reg112=reg26*reg64; T reg113=reg30*reg27;
    T reg114=reg26*reg32; T reg115=reg30*reg37; reg70=reg19*reg70; reg86=reg18*reg86; reg63=reg12-reg63;
    reg12=reg4*reg11; T reg116=reg2*reg39; reg16=reg16/reg58; T reg117=pow(reg87,2); T reg118=pow(reg77,2);
    T reg119=reg77*reg73; T reg120=reg5*reg4; reg82=reg20*reg82; T reg121=1-var_inter[1]; T reg122=reg76+reg7;
    T reg123=1-var_inter[0]; reg59=reg59/reg58; reg52=reg52/reg58; reg81=reg18*reg81; reg85=reg85/reg58;
    T reg124=reg6*reg4; reg0=reg120*reg0; reg99=reg55+reg99; reg114=reg115+reg114; reg55=reg117*reg52;
    reg115=reg90*reg3; reg78=reg65+reg78; reg65=reg12*reg80; T reg125=reg88*reg59; T reg126=reg90*reg85;
    T reg127=reg71*reg91; reg34=reg68+reg34; reg68=reg71*reg80; reg36=reg94+reg36; reg94=reg61*reg56;
    reg105=reg110+reg105; reg110=reg116*reg91; reg66=reg95+reg66; reg95=reg121*elem.pos(0)[0]; T reg128=reg121*elem.pos(1)[0];
    T reg129=reg116*reg80; reg109=reg74+reg109; reg74=reg121*elem.pos(0)[1]; T reg130=reg121*elem.pos(1)[1]; T reg131=reg88*reg3;
    reg69=reg93+reg69; reg93=elem.pos(0)[0]*reg123; T reg132=elem.pos(1)[0]*var_inter[0]; T reg133=reg12*reg91; reg21=reg67+reg21;
    reg67=reg123*elem.pos(0)[1]; T reg134=elem.pos(1)[1]*var_inter[0]; reg100=reg106+reg100; reg106=reg119*reg52; reg89=reg19*reg89;
    reg102=reg104+reg102; reg112=reg113+reg112; reg13=reg14-reg13; reg14=reg117*reg16; reg17=reg18*reg17;
    reg72=reg111+reg72; reg81=reg82-reg81; reg19=reg118*reg16; reg82=reg25*reg59; reg104=reg86+reg70;
    reg111=reg45*reg85; reg103=reg101+reg103; reg101=reg26*reg59; reg24=reg20*reg24; reg113=reg30*reg85;
    T reg135=reg119*reg16; reg107=reg108+reg107; reg92=reg18*reg92; reg122=reg15*reg122; reg63=reg63/reg58;
    reg108=reg118*reg52; reg108=reg103+reg108; reg103=reg93+reg132; reg133=reg21+reg133; reg21=elem.pos(2)[0]*var_inter[0];
    reg131=reg69+reg131; reg69=elem.pos(2)[1]*var_inter[0]; T reg136=reg96*reg0; reg106=reg102+reg106; reg102=reg67+reg134;
    reg65=reg78+reg65; reg78=reg97*reg0; reg115=reg99+reg115; reg14=reg112+reg14; reg1=reg124*reg1;
    reg19=reg72+reg19; reg135=reg107+reg135; reg101=reg113+reg101; reg72=reg117*reg63; reg55=reg114+reg55;
    reg104=reg15*reg104; reg68=reg36+reg68; reg13=reg20*reg13; reg92=reg24-reg92; reg20=reg62*reg3;
    reg122=reg81-reg122; reg94=reg105+reg94; reg110=reg66+reg110; reg24=reg119*reg63; reg36=elem.pos(2)[0]*var_inter[1];
    reg100=reg18*reg100; reg128=reg128-reg95; reg125=reg126+reg125; reg82=reg111+reg82; reg129=reg109+reg129;
    reg18=reg118*reg63; reg66=reg17+reg89; reg127=reg34+reg127; reg130=reg130-reg74; reg34=elem.pos(2)[1]*var_inter[1];
    reg81=reg54*reg19; reg99=reg79*reg106; reg105=reg54*reg135; reg72=reg101+reg72; reg101=reg6*reg87;
    reg107=reg5*reg77; reg109=reg77*reg87; reg111=reg65*reg127; reg112=reg110*reg68; reg113=reg129*reg127;
    reg136=reg131+reg136; reg114=reg116*reg1; reg126=elem.pos(3)[1]*var_inter[1]; reg130=reg34+reg130; reg34=reg133*reg68;
    reg66=reg15*reg66; reg15=elem.pos(3)[0]*var_inter[1]; reg21=reg21-reg103; reg131=elem.pos(3)[0]*reg123; reg128=reg36+reg128;
    reg69=reg69-reg102; reg100=reg13-reg100; reg13=elem.pos(3)[1]*reg123; reg36=reg12*reg1; reg78=reg115+reg78;
    reg20=reg94+reg20; reg94=reg51*reg0; reg122=reg122/reg58; reg104=reg92-reg104; reg24=reg125+reg24;
    reg18=reg82+reg18; reg82=reg79*reg55; reg92=reg54*reg14; reg115=reg79*reg108; reg125=reg26*reg14;
    T reg137=reg30*reg55; T reg138=reg45*reg106; T reg139=reg62*reg122; T reg140=reg42*reg122; reg104=reg104/reg58;
    T reg141=reg28*reg122; reg34=reg111-reg34; reg111=reg65*reg110; T reg142=reg133*reg129; T reg143=reg25*reg135;
    T reg144=reg5*reg87; T reg145=reg6*reg77; reg107=reg101+reg107; reg101=reg9*reg87; T reg146=reg8*reg77;
    reg135=reg26*reg135; reg130=reg130-reg126; reg114=reg136+reg114; reg106=reg30*reg106; reg66=reg100-reg66;
    reg128=reg128-reg15; reg55=reg45*reg55; reg14=reg25*reg14; reg131=reg21+reg131; reg21=reg26*reg19;
    reg100=reg45*reg108; reg19=reg25*reg19; reg13=reg69+reg13; reg108=reg30*reg108; reg36=reg78+reg36;
    reg94=reg20+reg94; reg20=reg71*reg1; reg92=reg82+reg92; reg112=reg113-reg112; reg69=reg109*reg72;
    reg78=reg109*reg18; reg81=reg115+reg81; reg82=reg109*reg24; reg105=reg99+reg105; reg69=reg92+reg69;
    reg92=reg130*reg131; reg99=reg62*reg141; reg113=reg144*reg104; reg115=reg145*reg104; reg136=reg107*reg104;
    T reg147=reg117*reg18; reg21=reg108+reg21; reg125=reg137+reg125; reg20=reg94+reg20; reg82=reg105+reg82;
    reg19=reg100+reg19; reg94=reg117*reg72; reg18=reg118*reg18; reg100=reg62*reg139; reg105=reg128*reg13;
    reg58=reg66/reg58; reg78=reg81+reg78; reg72=reg118*reg72; reg66=reg112*reg36; reg14=reg55+reg14;
    reg146=reg101+reg146; reg135=reg106+reg135; reg55=reg62*reg140; reg143=reg138+reg143; reg81=reg118*reg24;
    reg101=reg9*reg77; reg106=reg114*reg34; reg108=reg8*reg87; reg24=reg117*reg24; reg142=reg111-reg142;
    reg100=reg82+reg100; reg55=reg78+reg55; reg24=reg135+reg24; reg78=reg23*reg139; reg82=reg107*reg115;
    reg147=reg21+reg147; reg21=reg23*reg140; reg111=reg65*reg114; reg135=reg36*reg129; reg137=reg22*reg141;
    reg72=reg14+reg72; reg14=reg20*reg142; reg81=reg143+reg81; reg139=reg22*reg139; reg138=reg129*reg20;
    reg106=reg66-reg106; reg92=reg105-reg92; reg66=reg114*reg68; reg99=reg69+reg99; reg69=reg108*reg58;
    reg105=reg101*reg58; reg143=reg146*reg58; T reg148=reg107*reg113; T reg149=reg36*reg68; reg11=reg77*reg11;
    T reg150=reg65*reg20; reg94=reg125+reg94; reg140=reg22*reg140; reg18=reg19+reg18; reg19=reg107*reg136;
    reg141=reg23*reg141; reg10=reg87*reg10; reg125=reg114*reg127; T reg151=reg36*reg127; T reg152=reg110*reg20;
    reg138=reg66-reg138; reg66=reg10*reg115; reg21=reg147+reg21; reg147=reg133*reg114; T reg153=reg36*reg110;
    reg14=reg106+reg14; reg111=reg135-reg111; reg150=reg149-reg150; reg106=reg133*reg20; reg135=reg10*reg113;
    reg141=reg94+reg141; reg140=reg18+reg140; reg115=reg11*reg115; reg82=reg55+reg82; reg18=reg146*reg105;
    reg55=reg146*reg69; reg148=reg99+reg148; reg19=reg100+reg19; reg94=reg146*reg143; reg13=reg13/reg92;
    reg128=reg128/reg92; reg131=reg131/reg92; reg130=reg130/reg92; reg139=reg81+reg139; reg81=reg11*reg136;
    reg39=reg77*reg39; reg57=reg87*reg57; reg78=reg24+reg78; reg136=reg10*reg136; reg113=reg11*reg113;
    reg137=reg72+reg137; reg24=reg39*reg69; reg72=var_inter[1]*reg131; reg113=reg137+reg113; reg150=reg150/reg14;
    reg106=reg151-reg106; reg94=reg19+reg94; reg19=reg57*reg143; reg66=reg21+reg66; reg138=reg138/reg14;
    reg143=reg39*reg143; reg81=reg139+reg81; reg152=reg125-reg152; reg21=reg121*reg13; reg136=reg78+reg136;
    reg78=reg57*reg105; reg69=reg57*reg69; reg135=reg141+reg135; reg99=reg121*reg131; reg100=reg123*reg128;
    reg125=var_inter[0]*reg130; reg55=reg148+reg55; reg137=var_inter[0]*reg128; reg111=reg111/reg14; reg115=reg140+reg115;
    reg105=reg39*reg105; reg147=reg153-reg147; reg139=reg123*reg130; reg140=1-(*f.m).resolution; reg18=reg82+reg18;
    reg82=var_inter[1]*reg13; reg141=(*f.m).resolution*reg111; reg148=reg18*reg140; reg149=(*f.m).resolution*reg150; reg34=reg34/reg14;
    reg143=reg81+reg143; reg81=(*f.m).resolution*reg138; reg106=reg106/reg14; reg147=reg147/reg14; reg142=reg142/reg14;
    reg151=reg55*reg140; reg153=reg94*reg140; reg69=reg135+reg69; reg135=reg99-reg100; T reg154=reg125+reg21;
    T reg155=reg99+reg137; T reg156=reg139-reg21; T reg157=reg100+reg72; T reg158=reg82-reg125; reg112=reg112/reg14;
    T reg159=reg137-reg72; reg105=reg115+reg105; reg14=reg152/reg14; reg115=reg139+reg82; reg24=reg113+reg24;
    reg19=reg136+reg19; reg78=reg66+reg78; reg66=(*f.m).resolution*reg147; reg153=reg141+reg153; reg113=0.5*reg135;
    reg149=reg148-reg149; reg151=reg81+reg151; reg81=0.5*reg156; reg136=reg69*reg140; reg141=reg78*reg140;
    reg148=reg19*reg140; reg152=reg143*reg140; T reg160=reg105*reg140; T reg161=reg24*reg140; T reg162=0.5*reg158;
    T reg163=0.5*reg159; T reg164=0.5*reg154; T reg165=0.5*reg157; T reg166=0.5*reg115; T reg167=0.5*reg155;
    T reg168=(*f.m).resolution*reg34; T reg169=(*f.m).resolution*reg112; T reg170=(*f.m).resolution*reg142; T reg171=(*f.m).resolution*reg14; T reg172=(*f.m).resolution*reg106;
    T reg173=reg115*reg151; T reg174=reg156*reg151; T reg175=reg113*reg153; T reg176=reg163*reg153; T reg177=reg158*reg151;
    T reg178=reg154*reg151; T reg179=reg166*reg153; T reg180=reg164*reg153; T reg181=reg155*reg149; T reg182=reg157*reg149;
    T reg183=reg167*reg153; T reg184=reg135*reg149; T reg185=reg165*reg153; reg66=reg152-reg66; reg160=reg172+reg160;
    reg171=reg161-reg171; reg148=reg170+reg148; reg168=reg141-reg168; reg136=reg169+reg136; reg141=reg162*reg153;
    reg152=reg159*reg149; reg161=reg81*reg153; reg169=reg115*reg171; reg170=reg162*reg66; reg172=reg159*reg160;
    T reg186=reg165*reg66; T reg187=reg163*reg66; T reg188=reg164*reg66; T reg189=reg158*reg171; T reg190=reg155*reg160;
    T reg191=reg167*reg66; T reg192=reg156*reg171; T reg193=reg113*reg66; T reg194=reg154*reg171; T reg195=reg135*reg160;
    T reg196=reg81*reg66; T reg197=reg162*reg148; T reg198=reg159*reg168; reg141=reg152+reg141; reg152=reg163*reg148;
    T reg199=reg158*reg136; reg176=reg177+reg176; reg177=reg164*reg148; T reg200=reg155*reg168; reg180=reg180-reg181;
    T reg201=reg167*reg148; T reg202=reg154*reg136; reg178=reg178-reg183; T reg203=reg81*reg148; T reg204=reg135*reg168;
    reg161=reg184+reg161; reg184=reg166*reg66; T reg205=reg157*reg160; T reg206=reg156*reg136; reg185=reg185-reg173;
    T reg207=reg165*reg148; T reg208=reg113*reg148; reg175=reg174+reg175; reg174=reg115*reg136; reg182=reg182-reg179;
    T reg209=reg157*reg168; T reg210=reg166*reg148; reg185=2*reg185; reg186=reg186-reg169; reg182=2*reg182;
    reg207=reg207-reg174; reg188=reg188-reg190; reg175=2*reg175; reg152=reg199+reg152; reg187=reg189+reg187;
    reg170=reg172+reg170; reg208=reg206+reg208; reg180=2*reg180; reg176=2*reg176; reg197=reg198+reg197;
    reg178=2*reg178; reg193=reg192+reg193; reg203=reg204+reg203; reg205=reg205-reg184; reg141=2*reg141;
    reg209=reg209-reg210; reg196=reg195+reg196; reg177=reg177-reg200; reg161=2*reg161; reg202=reg202-reg201;
    reg194=reg194-reg191; reg172=reg158*reg177; reg189=reg158*reg152; reg192=reg166*reg141; reg195=reg158*reg202;
    reg198=reg163*reg161; reg199=reg163*reg180; reg204=reg157*reg170; reg206=reg166*reg176; T reg211=reg158*reg203;
    T reg212=reg163*reg176; T reg213=reg156*reg197; T reg214=reg155*reg196; T reg215=reg156*reg203; T reg216=reg155*reg194;
    T reg217=reg113*reg161; T reg218=reg164*reg180; T reg219=reg164*reg176; T reg220=reg155*reg188; T reg221=reg155*reg187;
    T reg222=reg164*reg141; T reg223=reg155*reg170; T reg224=reg166*reg182; T reg225=reg164*reg185; T reg226=reg155*reg186;
    T reg227=reg157*reg205; T reg228=reg164*reg182; T reg229=reg155*reg205; T reg230=reg166*reg185; T reg231=reg158*reg208;
    T reg232=reg163*reg175; T reg233=reg157*reg186; T reg234=reg163*reg178; T reg235=reg166*reg178; T reg236=reg159*reg170;
    T reg237=reg162*reg141; T reg238=reg165*reg182; T reg239=reg159*reg186; T reg240=reg162*reg185; T reg241=reg115*reg207;
    T reg242=reg159*reg205; T reg243=reg162*reg182; T reg244=reg165*reg185; T reg245=reg165*reg175; T reg246=reg115*reg208;
    T reg247=reg115*reg197; T reg248=reg165*reg161; T reg249=reg115*reg203; T reg250=reg165*reg141; T reg251=reg165*reg178;
    T reg252=reg115*reg202; T reg253=reg115*reg152; T reg254=reg165*reg180; T reg255=reg115*reg177; T reg256=reg165*reg176;
    T reg257=reg158*reg197; T reg258=reg157*reg187; T reg259=reg163*reg141; T reg260=reg166*reg180; T reg261=reg158*reg209;
    T reg262=reg163*reg182; T reg263=reg157*reg188; T reg264=reg159*reg193; T reg265=reg162*reg175; T reg266=reg166*reg161;
    T reg267=reg159*reg196; T reg268=reg162*reg161; T reg269=reg157*reg196; T reg270=reg159*reg194; T reg271=reg162*reg178;
    T reg272=reg166*reg175; T reg273=reg159*reg188; T reg274=reg162*reg180; T reg275=reg157*reg193; T reg276=reg159*reg187;
    T reg277=reg162*reg176; T reg278=reg115*reg209; T reg279=reg156*reg207; reg170=reg135*reg170; T reg280=reg81*reg141;
    reg186=reg135*reg186; T reg281=reg81*reg185; T reg282=reg113*reg175; T reg283=reg113*reg185; reg205=reg135*reg205;
    T reg284=reg81*reg182; T reg285=reg156*reg177; T reg286=reg154*reg208; T reg287=reg167*reg175; T reg288=reg167*reg178;
    reg203=reg154*reg203; T reg289=reg167*reg161; T reg290=reg113*reg180; T reg291=reg154*reg202; reg177=reg154*reg177;
    T reg292=reg167*reg180; T reg293=reg157*reg194; T reg294=reg113*reg141; T reg295=reg156*reg209; T reg296=reg135*reg193;
    T reg297=reg81*reg175; reg196=reg135*reg196; T reg298=reg81*reg161; T reg299=reg163*reg185; T reg300=reg164*reg178;
    T reg301=reg113*reg182; reg194=reg135*reg194; T reg302=reg81*reg178; T reg303=reg156*reg152; reg188=reg135*reg188;
    reg180=reg81*reg180; T reg304=reg158*reg207; reg187=reg135*reg187; T reg305=reg81*reg176; T reg306=reg113*reg176;
    reg207=reg154*reg207; reg202=reg156*reg202; reg185=reg167*reg185; reg141=reg167*reg141; reg178=reg113*reg178;
    reg197=reg154*reg197; reg209=reg154*reg209; reg182=reg167*reg182; reg176=reg167*reg176; reg152=reg154*reg152;
    reg161=reg164*reg161; reg175=reg164*reg175; reg193=reg155*reg193; reg208=reg156*reg208; reg303=reg306+reg303;
    reg216=reg300-reg216; reg295=reg301+reg295; reg278=reg238-reg278; reg277=reg276+reg277; reg180=reg188+reg180;
    reg274=reg273+reg274; reg185=reg207-reg185; reg193=reg175-reg193; reg305=reg187+reg305; reg272=reg275-reg272;
    reg271=reg270+reg271; reg255=reg254-reg255; reg220=reg218-reg220; reg279=reg283+reg279; reg297=reg296+reg297;
    reg249=reg248-reg249; reg214=reg161-reg214; reg247=reg250-reg247; reg294=reg213+reg294; reg246=reg245-reg246;
    reg298=reg196+reg298; reg299=reg304+reg299; reg215=reg217+reg215; reg243=reg242+reg243; reg235=reg293-reg235;
    reg182=reg209-reg182; reg252=reg251-reg252; reg241=reg244-reg241; reg240=reg239+reg240; reg253=reg256-reg253;
    reg302=reg194+reg302; reg237=reg236+reg237; reg223=reg222-reg223; reg212=reg189+reg212; reg284=reg205+reg284;
    reg285=reg290+reg285; reg224=reg227-reg224; reg206=reg258-reg206; reg199=reg172+reg199; reg176=reg152-reg176;
    reg287=reg286-reg287; reg195=reg234+reg195; reg226=reg225-reg226; reg192=reg204-reg192; reg198=reg211+reg198;
    reg289=reg203-reg289; reg208=reg282+reg208; reg288=reg291-reg288; reg232=reg231+reg232; reg229=reg228-reg229;
    reg292=reg177-reg292; reg230=reg233-reg230; reg268=reg267+reg268; reg141=reg197-reg141; reg202=reg178+reg202;
    reg259=reg257+reg259; reg260=reg263-reg260; reg280=reg170+reg280; reg266=reg269-reg266; reg265=reg264+reg265;
    reg262=reg261+reg262; reg281=reg186+reg281; reg221=reg219-reg221; reg266=reg92*reg266; reg253=reg92*reg253;
    reg230=reg92*reg230; reg192=reg92*reg192; reg272=reg92*reg272; reg299=reg92*reg299; reg278=reg92*reg278;
    reg215=reg92*reg215; reg208=reg92*reg208; reg206=reg92*reg206; reg247=reg92*reg247; reg202=reg92*reg202;
    reg224=reg92*reg224; reg285=reg92*reg285; reg294=reg92*reg294; reg303=reg92*reg303; reg260=reg92*reg260;
    reg241=reg92*reg241; reg265=reg92*reg265; reg281=reg92*reg281; reg262=reg92*reg262; reg259=reg92*reg259;
    reg284=reg92*reg284; reg212=reg92*reg212; reg199=reg92*reg199; reg287=reg92*reg287; reg195=reg92*reg195;
    reg289=reg92*reg289; reg198=reg92*reg198; reg288=reg92*reg288; reg232=reg92*reg232; reg229=reg92*reg229;
    reg292=reg92*reg292; reg226=reg92*reg226; reg223=reg92*reg223; reg176=reg92*reg176; reg221=reg92*reg221;
    reg141=reg92*reg141; reg220=reg92*reg220; reg185=reg92*reg185; reg216=reg92*reg216; reg235=reg92*reg235;
    reg214=reg92*reg214; reg182=reg92*reg182; reg193=reg92*reg193; reg295=reg92*reg295; reg255=reg92*reg255;
    reg252=reg92*reg252; reg297=reg92*reg297; reg249=reg92*reg249; reg246=reg92*reg246; reg298=reg92*reg298;
    reg243=reg92*reg243; reg240=reg92*reg240; reg302=reg92*reg302; reg237=reg92*reg237; reg277=reg92*reg277;
    reg280=reg92*reg280; reg268=reg92*reg268; reg305=reg92*reg305; reg271=reg92*reg271; reg279=reg92*reg279;
    reg274=reg92*reg274; reg180=reg92*reg180; T tmp_2_7=ponderation*reg182; T tmp_2_0=ponderation*reg287; T tmp_4_6=ponderation*reg299;
    T tmp_2_6=ponderation*reg185; T tmp_1_5=ponderation*reg280; T tmp_0_0=ponderation*reg208; T tmp_0_7=ponderation*reg295; T tmp_1_4=ponderation*reg305;
    T tmp_2_5=ponderation*reg141; T tmp_1_3=ponderation*reg180; T tmp_1_0=ponderation*reg297; T tmp_0_3=ponderation*reg285; T tmp_2_4=ponderation*reg176;
    T tmp_0_6=ponderation*reg279; T tmp_1_1=ponderation*reg298; T tmp_0_2=ponderation*reg202; T tmp_1_6=ponderation*reg281; T tmp_2_3=ponderation*reg292;
    T tmp_0_4=ponderation*reg303; T tmp_1_7=ponderation*reg284; T tmp_2_2=ponderation*reg288; T tmp_2_1=ponderation*reg289; T tmp_1_2=ponderation*reg302;
    T tmp_4_5=ponderation*reg259; T tmp_7_1=ponderation*reg266; T tmp_4_7=ponderation*reg262; T tmp_5_0=ponderation*reg265; T tmp_7_0=ponderation*reg272;
    T tmp_5_1=ponderation*reg268; T tmp_5_2=ponderation*reg271; T tmp_6_7=ponderation*reg278; T tmp_5_3=ponderation*reg274; T tmp_5_4=ponderation*reg277;
    T tmp_6_6=ponderation*reg241; T tmp_5_5=ponderation*reg237; T tmp_5_6=ponderation*reg240; T tmp_6_5=ponderation*reg247; T tmp_5_7=ponderation*reg243;
    T tmp_6_0=ponderation*reg246; T tmp_6_4=ponderation*reg253; T tmp_6_1=ponderation*reg249; T tmp_6_2=ponderation*reg252; T tmp_6_3=ponderation*reg255;
    T tmp_3_0=ponderation*reg193; T tmp_0_1=ponderation*reg215; T tmp_3_1=ponderation*reg214; T tmp_7_2=ponderation*reg235; T tmp_0_5=ponderation*reg294;
    T tmp_3_2=ponderation*reg216; T tmp_3_3=ponderation*reg220; T tmp_7_7=ponderation*reg224; T tmp_3_4=ponderation*reg221; T tmp_3_5=ponderation*reg223;
    T tmp_7_6=ponderation*reg230; T tmp_3_6=ponderation*reg226; T tmp_3_7=ponderation*reg229; T tmp_7_5=ponderation*reg192; T tmp_4_0=ponderation*reg232;
    T tmp_4_1=ponderation*reg198; T tmp_7_4=ponderation*reg206; T tmp_4_2=ponderation*reg195; T tmp_4_3=ponderation*reg199; T tmp_7_3=ponderation*reg260;
    T tmp_4_4=ponderation*reg212;
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
    T reg0=pow((*f.m).v1[1],2); T reg1=pow((*f.m).v1[0],2); T reg2=pow((*f.m).v2[1],2); T reg3=pow((*f.m).v2[0],2); reg0=reg1+reg0;
    reg1=pow((*f.m).v1[2],2); T reg4=2*(*f.m).shear_modulus_23; T reg5=pow((*f.m).v2[2],2); T reg6=2*(*f.m).shear_modulus_13; reg2=reg3+reg2;
    reg1=reg0+reg1; reg5=reg2+reg5; reg6=1.0/reg6; reg1=pow(reg1,0.5); reg0=2*(*f.m).shear_modulus_12;
    reg4=1.0/reg4; reg2=(*f.m).v1[1]/reg1; reg3=reg6*reg4; reg0=1.0/reg0; T reg7=(*f.m).v1[0]/reg1;
    reg5=pow(reg5,0.5); T reg8=2*reg7; T reg9=2*reg2; T reg10=reg0*reg3; T reg11=(*f.m).v2[1]/reg5;
    T reg12=(*f.m).v2[0]/reg5; T reg13=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; reg1=(*f.m).v1[2]/reg1; T reg14=1.0/(*f.m).elastic_modulus_3; T reg15=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    reg5=(*f.m).v2[2]/reg5; T reg16=2*reg1; T reg17=reg11*reg9; T reg18=reg13*reg10; T reg19=reg12*reg8;
    T reg20=1.0/(*f.m).elastic_modulus_2; T reg21=1.0/(*f.m).elastic_modulus_1; T reg22=pow(reg11,2); T reg23=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg24=reg15*reg10;
    T reg25=reg14*reg10; T reg26=pow(reg12,2); T reg27=pow(reg5,2); T reg28=reg15*reg18; T reg29=reg19*reg21;
    T reg30=reg25*reg23; T reg31=reg15*reg24; T reg32=reg22*reg20; T reg33=reg17*reg20; T reg34=reg22*reg23;
    T reg35=reg19*reg23; T reg36=reg26*reg21; reg16=reg5*reg16; T reg37=reg17*reg23; T reg38=reg26*reg23;
    T reg39=reg25*reg20; T reg40=reg22*reg15; T reg41=reg16*reg15; T reg42=reg27*reg13; T reg43=reg19*reg13;
    reg35=reg33-reg35; reg33=reg26*reg13; T reg44=pow(reg7,2); reg34=reg36-reg34; reg36=pow(reg2,2);
    T reg45=reg17*reg15; reg37=reg29-reg37; reg29=reg16*reg13; reg31=reg39-reg31; reg39=reg20*reg18;
    T reg46=reg27*reg15; reg28=reg30+reg28; T reg47=reg23*reg24; reg38=reg32-reg38; reg46=reg38-reg46;
    reg32=pow(reg1,2); reg38=reg2*reg11; T reg48=reg7*reg12; reg41=reg35-reg41; reg35=reg21*reg31;
    T reg49=reg7*reg11; T reg50=reg2*reg12; T reg51=reg47+reg39; reg42=reg34-reg42; reg34=reg44*reg23;
    reg29=reg37-reg29; reg37=reg36*reg23; T reg52=reg36*reg20; T reg53=reg23*reg28; reg45=reg43+reg45;
    reg16=reg16*reg14; reg40=reg33+reg40; reg33=reg27*reg14; reg43=reg44*reg21; T reg54=reg32*reg13;
    T reg55=reg22*reg46; T reg56=2*reg12; reg25=reg25*reg21; T reg57=reg26*reg29; T reg58=reg22*reg41;
    T reg59=reg36*reg15; T reg60=reg26*reg42; T reg61=reg44*reg13; T reg62=reg13*reg18; T reg63=reg23*reg10;
    T reg64=reg48*reg42; T reg65=reg38*reg46; T reg66=reg1*reg5; T reg67=reg50+reg49; T reg68=reg14*reg3;
    T reg69=reg0*reg4; T reg70=reg48*reg29; T reg71=reg38*reg41; reg40=reg33-reg40; reg45=reg16-reg45;
    reg53=reg35-reg53; reg16=reg15*reg3; reg33=reg13*reg51; reg3=reg13*reg3; reg35=reg32*reg15;
    T reg72=reg44*reg42; T reg73=reg36*reg46; T reg74=reg7*reg5; T reg75=reg1*reg12; T reg76=reg1*reg11;
    T reg77=reg2*reg5; reg37=reg43-reg37; reg43=reg44*reg29; T reg78=reg36*reg41; T reg79=reg12*reg11;
    reg10=reg20*reg10; T reg80=reg13*reg24; reg34=reg52-reg34; reg33=reg53-reg33; reg54=reg37-reg54;
    reg37=reg12*reg5; reg35=reg34-reg35; reg80=reg30+reg80; reg30=reg15*reg3; reg55=reg60+reg55;
    reg34=reg27*reg40; reg58=reg57+reg58; reg52=reg27*reg45; reg62=reg25-reg62; reg25=reg13*reg69;
    reg24=reg21*reg24; reg53=reg14*reg69; reg65=reg64+reg65; reg57=reg66*reg40; reg60=reg20*reg68;
    reg64=reg13*reg63; reg71=reg70+reg71; reg70=reg66*reg45; reg69=reg15*reg69; reg18=reg23*reg18;
    T reg81=reg32*reg14; T reg82=reg2*reg8; T reg83=reg77-reg76; reg59=reg61+reg59; reg61=reg79*reg0;
    T reg84=reg67*reg0; reg68=reg23*reg68; reg73=reg72+reg73; reg72=reg32*reg40; T reg85=reg15*reg16;
    reg78=reg43+reg78; reg43=reg11*reg56; T reg86=reg74+reg75; T reg87=reg32*reg45; T reg88=reg13*reg10;
    T reg89=reg0*reg6; reg18=reg24+reg18; T reg90=reg13*reg89; reg31=reg31/reg33; reg64=reg24+reg64;
    reg62=reg62/reg33; reg28=reg28/reg33; reg88=reg47+reg88; reg80=reg80/reg33; reg24=reg26*reg54;
    T reg91=reg22*reg35; T reg92=reg15*reg89; reg34=reg55+reg34; reg55=reg43*reg61; reg3=reg20*reg3;
    reg30=reg68+reg30; reg52=reg58+reg52; reg58=reg43*reg84; reg89=reg14*reg89; reg14=reg15*reg25;
    reg57=reg65+reg57; reg65=reg67*reg61; reg68=reg20*reg53; reg53=reg23*reg53; reg70=reg71+reg70;
    reg71=reg67*reg84; T reg93=reg15*reg69; reg10=reg21*reg10; T reg94=reg11*reg5; T reg95=reg7*reg2;
    reg77=reg76+reg77; reg76=2*reg83; reg59=reg81-reg59; reg81=2*reg11; T reg96=reg5*reg56;
    T reg97=reg1*reg8; reg74=reg75-reg74; reg63=reg23*reg63; reg75=reg37*reg6; T reg98=reg86*reg6;
    T reg99=reg44*reg54; T reg100=reg36*reg35; reg72=reg73+reg72; reg73=reg82*reg61; reg16=reg23*reg16;
    reg87=reg78+reg87; reg78=reg82*reg84; reg85=reg60-reg85; reg60=reg44*reg80; T reg101=reg26*reg62;
    T reg102=reg43*reg28; T reg103=reg82*reg31; T reg104=reg22*reg28; T reg105=reg36*reg31; T reg106=reg15*reg90;
    reg15=reg15*reg92; T reg107=reg23*reg89; reg89=reg20*reg89; reg14=reg53+reg14; reg93=reg68-reg93;
    reg53=reg86*reg98; reg71=reg70+reg71; reg68=reg86*reg75; reg65=reg57+reg65; reg57=reg38*reg35;
    reg70=reg48*reg54; T reg108=reg96*reg98; reg58=reg52+reg58; reg52=reg96*reg75; reg55=reg34+reg55;
    reg34=reg27*reg59; reg91=reg24+reg91; reg24=reg97*reg98; reg78=reg87+reg78; reg87=reg97*reg75;
    reg73=reg72+reg73; reg72=reg32*reg59; reg100=reg99+reg100; reg99=reg77*reg4; T reg109=reg94*reg4;
    reg0=reg95*reg0; T reg110=reg43*reg62; T reg111=reg82*reg80; T reg112=reg22*reg62; T reg113=reg36*reg80;
    T reg114=reg26*reg28; reg25=reg20*reg25; reg18=reg18/reg33; T reg115=reg44*reg31; reg69=reg23*reg69;
    reg63=reg10-reg63; reg10=reg1*reg9; reg85=reg21*reg85; T reg116=1-var_inter[0]; T reg117=reg5*reg81;
    reg64=reg64/reg33; reg51=reg51/reg33; reg88=reg88/reg33; T reg118=reg16+reg3; T reg119=1-var_inter[1];
    T reg120=pow(reg83,2); T reg121=pow(reg74,2); T reg122=reg74*reg76; T reg123=reg7*reg1; reg30=reg23*reg30;
    T reg124=reg43*reg0; T reg125=reg119*elem.pos(1)[1]; reg108=reg58+reg108; reg58=reg36*reg88; T reg126=reg119*elem.pos(0)[1];
    reg52=reg55+reg52; reg55=reg117*reg109; T reg127=reg43*reg64; T reg128=reg82*reg88; T reg129=reg119*elem.pos(1)[0];
    T reg130=reg119*elem.pos(0)[0]; T reg131=reg22*reg64; reg34=reg91+reg34; reg91=elem.pos(0)[0]*reg116; T reg132=elem.pos(1)[0]*var_inter[0];
    T reg133=reg10*reg99; reg24=reg78+reg24; reg78=reg116*elem.pos(0)[1]; T reg134=elem.pos(1)[1]*var_inter[0]; T reg135=reg10*reg109;
    reg87=reg73+reg87; reg73=reg82*reg0; reg72=reg100+reg72; reg6=reg123*reg6; reg100=reg2*reg1;
    reg114=reg115+reg114; reg115=reg120*reg51; reg101=reg60+reg101; reg60=reg122*reg51; reg102=reg103+reg102;
    reg92=reg23*reg92; reg103=reg121*reg51; reg104=reg105+reg104; reg90=reg20*reg90; reg68=reg65+reg68;
    reg20=reg77*reg109; reg106=reg107+reg106; reg15=reg89-reg15; reg30=reg85-reg30; reg118=reg13*reg118;
    reg65=reg69+reg25; reg14=reg23*reg14; reg53=reg71+reg53; reg71=reg77*reg99; reg93=reg21*reg93;
    reg85=reg117*reg99; reg89=reg26*reg64; reg105=reg44*reg88; reg57=reg70+reg57; reg70=reg122*reg18;
    reg110=reg111+reg110; reg63=reg63/reg33; reg107=reg66*reg59; reg111=reg121*reg18; reg112=reg113+reg112;
    reg113=reg120*reg18; reg71=reg53+reg71; reg85=reg108+reg85; reg115=reg114+reg115; reg53=reg67*reg0;
    reg108=elem.pos(2)[0]*var_inter[0]; reg114=elem.pos(2)[0]*var_inter[1]; reg124=reg34+reg124; reg34=reg96*reg6; T reg136=reg91+reg132;
    reg133=reg24+reg133; reg20=reg68+reg20; reg4=reg100*reg4; reg129=reg129-reg130; reg55=reg52+reg55;
    reg73=reg72+reg73; reg24=reg97*reg6; reg125=reg125-reg126; reg107=reg57+reg107; reg135=reg87+reg135;
    reg52=reg78+reg134; reg57=elem.pos(2)[1]*var_inter[1]; reg68=elem.pos(2)[1]*var_inter[0]; reg103=reg104+reg103; reg72=reg122*reg63;
    reg87=reg92+reg90; reg127=reg128+reg127; reg106=reg23*reg106; reg23=reg121*reg63; reg15=reg21*reg15;
    reg131=reg58+reg131; reg65=reg13*reg65; reg118=reg30-reg118; reg14=reg93-reg14; reg21=reg120*reg63;
    reg89=reg105+reg89; reg70=reg110+reg70; reg113=reg101+reg113; reg60=reg102+reg60; reg111=reg112+reg111;
    reg30=reg79*reg113; reg58=reg135*reg71; reg93=reg55*reg71; reg101=elem.pos(3)[1]*var_inter[1]; reg125=reg57+reg125;
    reg57=reg85*reg20; reg102=reg10*reg4; reg106=reg15-reg106; reg87=reg13*reg87; reg13=reg95*reg103;
    reg15=reg79*reg111; reg24=reg73+reg24; reg129=reg114+reg129; reg73=elem.pos(3)[0]*var_inter[1]; reg104=reg133*reg20;
    reg105=reg95*reg60; reg110=reg79*reg70; reg72=reg127+reg72; reg34=reg124+reg34; reg68=reg68-reg52;
    reg112=reg117*reg4; reg114=reg7*reg74; reg124=reg2*reg83; reg21=reg89+reg21; reg53=reg107+reg53;
    reg89=elem.pos(3)[0]*reg116; reg107=elem.pos(3)[1]*reg116; reg127=reg86*reg6; reg128=reg74*reg83; reg65=reg14-reg65;
    reg14=reg95*reg115; reg23=reg131+reg23; reg118=reg118/reg33; reg108=reg108-reg136; reg129=reg129-reg73;
    reg131=reg11*reg83; T reg137=reg48*reg118; T reg138=reg7*reg83; T reg139=reg12*reg74; reg114=reg124+reg114;
    reg124=reg77*reg4; reg127=reg53+reg127; reg53=reg2*reg74; T reg140=reg44*reg60; reg57=reg93-reg57;
    reg93=reg26*reg70; T reg141=reg44*reg115; reg115=reg36*reg115; T reg142=reg22*reg113; reg113=reg26*reg113;
    T reg143=reg133*reg55; reg89=reg108+reg89; reg108=reg135*reg85; reg104=reg58-reg104; reg58=reg36*reg103;
    T reg144=reg22*reg111; reg103=reg44*reg103; reg65=reg65/reg33; reg60=reg36*reg60; reg70=reg22*reg70;
    reg111=reg26*reg111; reg107=reg68+reg107; reg68=reg38*reg118; T reg145=reg67*reg118; T reg146=reg128*reg72;
    reg110=reg105+reg110; reg87=reg106-reg87; reg105=reg128*reg23; reg15=reg13+reg15; reg102=reg24+reg102;
    reg13=reg128*reg21; reg30=reg14+reg30; reg125=reg125-reg101; reg112=reg34+reg112; reg14=reg114*reg65;
    reg24=reg125*reg89; reg111=reg103+reg111; reg34=reg120*reg23; reg103=reg53*reg65; reg106=reg120*reg21;
    reg113=reg141+reg113; reg141=reg138*reg65; T reg147=reg129*reg107; reg23=reg121*reg23; reg144=reg58+reg144;
    reg146=reg110+reg146; reg58=reg67*reg145; reg110=reg67*reg68; T reg148=reg112*reg104; reg105=reg15+reg105;
    reg70=reg60+reg70; reg15=reg121*reg72; reg124=reg127+reg124; reg60=reg67*reg137; reg127=reg57*reg102;
    reg72=reg120*reg72; reg13=reg30+reg13; reg93=reg140+reg93; reg21=reg121*reg21; reg142=reg115+reg142;
    reg143=reg108-reg143; reg139=reg131+reg139; reg33=reg87/reg33; reg30=reg11*reg74; reg87=reg12*reg83;
    reg15=reg70+reg15; reg70=reg17*reg68; reg108=reg17*reg145; reg23=reg144+reg23; reg115=reg87*reg33;
    reg131=reg30*reg33; reg148=reg127-reg148; reg127=reg124*reg143; reg140=reg17*reg137; reg137=reg19*reg137;
    reg106=reg113+reg106; reg21=reg142+reg21; reg113=reg55*reg124; reg142=reg139*reg33; reg9=reg74*reg9;
    reg24=reg147-reg24; reg72=reg93+reg72; reg145=reg19*reg145; reg93=reg112*reg20; reg144=reg135*reg124;
    reg8=reg83*reg8; reg147=reg114*reg14; reg58=reg146+reg58; reg146=reg114*reg103; reg68=reg19*reg68;
    reg34=reg111+reg34; reg110=reg105+reg110; reg105=reg102*reg55; reg111=reg114*reg141; reg60=reg13+reg60;
    reg13=reg102*reg20; T reg149=reg135*reg112; T reg150=reg102*reg71; T reg151=reg85*reg124; reg127=reg148+reg127;
    reg148=reg112*reg71; reg125=reg125/reg24; T reg152=reg102*reg85; T reg153=reg133*reg112; reg56=reg83*reg56;
    reg81=reg74*reg81; reg145=reg72+reg145; reg72=reg8*reg14; T reg154=reg139*reg142; reg147=reg58+reg147;
    reg58=reg139*reg131; reg146=reg110+reg146; reg110=reg139*reg115; reg111=reg60+reg111; reg14=reg9*reg14;
    reg108=reg15+reg108; reg149=reg105-reg149; reg15=reg9*reg103; reg70=reg23+reg70; reg23=reg9*reg141;
    reg140=reg21+reg140; reg89=reg89/reg24; reg129=reg129/reg24; reg107=reg107/reg24; reg137=reg106+reg137;
    reg141=reg8*reg141; reg144=reg13-reg144; reg13=reg133*reg124; reg103=reg8*reg103; reg68=reg34+reg68;
    reg113=reg93-reg113; reg58=reg146+reg58; reg110=reg111+reg110; reg21=var_inter[0]*reg129; reg149=reg149/reg127;
    reg34=reg81*reg142; reg14=reg108+reg14; reg60=var_inter[0]*reg125; reg113=reg113/reg127; reg93=reg81*reg131;
    reg15=reg70+reg15; reg70=reg56*reg115; reg141=reg137+reg141; reg115=reg81*reg115; reg23=reg140+reg23;
    reg151=reg148-reg151; reg153=reg152-reg153; reg105=reg119*reg89; reg106=reg116*reg129; reg108=reg119*reg107;
    reg103=reg68+reg103; reg131=reg56*reg131; reg144=reg144/reg127; reg68=var_inter[1]*reg89; reg72=reg145+reg72;
    reg142=reg56*reg142; reg111=reg116*reg125; reg154=reg147+reg154; reg137=1-(*f.m).resolution; reg13=reg150-reg13;
    reg140=var_inter[1]*reg107; reg143=reg143/reg127; reg104=reg104/reg127; reg151=reg151/reg127; reg145=(*f.m).resolution*reg144;
    reg13=reg13/reg127; reg153=reg153/reg127; reg146=(*f.m).resolution*reg113; reg131=reg103+reg131; reg103=reg154*reg137;
    reg147=reg58*reg137; reg148=reg111+reg140; reg150=reg110*reg137; reg34=reg14+reg34; reg14=reg21-reg68;
    reg152=reg140-reg60; reg93=reg15+reg93; reg115=reg23+reg115; reg15=reg105+reg21; reg23=reg60+reg108;
    reg142=reg72+reg142; reg70=reg141+reg70; reg127=reg57/reg127; reg57=reg105-reg106; reg72=reg111-reg108;
    reg141=reg106+reg68; T reg155=(*f.m).resolution*reg149; T reg156=(*f.m).resolution*reg153; T reg157=(*f.m).resolution*reg13; reg103=reg155+reg103;
    reg145=reg147-reg145; reg150=reg146+reg150; reg146=reg34*reg137; reg147=reg93*reg137; reg155=reg115*reg137;
    T reg158=reg142*reg137; T reg159=reg131*reg137; T reg160=reg70*reg137; T reg161=0.5*reg152; T reg162=0.5*reg14;
    T reg163=0.5*reg23; T reg164=0.5*reg15; T reg165=0.5*reg72; T reg166=0.5*reg57; T reg167=0.5*reg148;
    T reg168=0.5*reg141; T reg169=(*f.m).resolution*reg104; T reg170=(*f.m).resolution*reg127; T reg171=(*f.m).resolution*reg143; T reg172=(*f.m).resolution*reg151;
    T reg173=reg166*reg103; T reg174=reg162*reg103; T reg175=reg152*reg150; T reg176=reg167*reg103; T reg177=reg141*reg145;
    T reg178=reg163*reg103; T reg179=reg15*reg145; T reg180=reg57*reg145; T reg181=reg164*reg103; T reg182=reg23*reg150;
    T reg183=reg148*reg150; T reg184=reg168*reg103; T reg185=reg165*reg103; reg160=reg170+reg160; reg169=reg159-reg169;
    reg158=reg171+reg158; reg172=reg155-reg172; reg147=reg157+reg147; reg156=reg146-reg156; reg146=reg161*reg103;
    reg155=reg14*reg145; reg157=reg72*reg150; reg159=reg163*reg156; reg170=reg164*reg156; reg171=reg162*reg156;
    T reg186=reg152*reg172; T reg187=reg14*reg147; T reg188=reg161*reg156; T reg189=reg15*reg147; T reg190=reg148*reg172;
    T reg191=reg165*reg156; T reg192=reg57*reg147; T reg193=reg168*reg156; T reg194=reg141*reg147; T reg195=reg23*reg172;
    T reg196=reg167*reg156; reg185=reg180+reg185; reg180=reg57*reg169; T reg197=reg165*reg158; reg182=reg182-reg181;
    T reg198=reg23*reg160; T reg199=reg164*reg158; reg178=reg178-reg179; T reg200=reg15*reg169; T reg201=reg163*reg158;
    reg174=reg175+reg174; reg175=reg152*reg160; T reg202=reg162*reg158; reg146=reg155+reg146; reg155=reg14*reg169;
    T reg203=reg161*reg158; reg177=reg177-reg176; reg173=reg157+reg173; reg157=reg141*reg169; T reg204=reg148*reg160;
    T reg205=reg168*reg158; T reg206=reg72*reg160; T reg207=reg166*reg158; T reg208=reg167*reg158; reg184=reg184-reg183;
    reg188=reg187+reg188; reg201=reg201-reg200; reg207=reg206+reg207; reg178=2*reg178; reg171=reg186+reg171;
    reg194=reg194-reg196; reg173=2*reg173; reg174=2*reg174; reg159=reg159-reg189; reg193=reg193-reg190;
    reg205=reg205-reg204; reg177=2*reg177; reg198=reg198-reg199; reg202=reg175+reg202; reg197=reg180+reg197;
    reg146=2*reg146; reg203=reg155+reg203; reg184=2*reg184; reg195=reg195-reg170; reg185=2*reg185;
    reg182=2*reg182; reg157=reg157-reg208; reg191=reg192+reg191; reg155=reg164*reg174; reg175=reg168*reg177;
    reg180=reg148*reg157; reg186=reg164*reg177; reg187=reg166*reg185; reg192=reg23*reg203; reg206=reg23*reg157;
    T reg209=reg141*reg194; T reg210=reg164*reg146; T reg211=reg166*reg184; T reg212=reg167*reg177; T reg213=reg164*reg184;
    T reg214=reg23*reg205; T reg215=reg162*reg146; T reg216=reg152*reg203; T reg217=reg162*reg174; T reg218=reg152*reg202;
    T reg219=reg152*reg157; T reg220=reg162*reg177; T reg221=reg15*reg194; T reg222=reg163*reg177; T reg223=reg14*reg188;
    T reg224=reg161*reg146; T reg225=reg15*reg193; T reg226=reg163*reg184; T reg227=reg14*reg193; T reg228=reg161*reg184;
    T reg229=reg15*reg188; T reg230=reg163*reg146; T reg231=reg14*reg194; T reg232=reg161*reg177; T reg233=reg15*reg171;
    T reg234=reg168*reg184; T reg235=reg148*reg205; T reg236=reg15*reg159; T reg237=reg163*reg174; T reg238=reg163*reg178;
    T reg239=reg166*reg178; T reg240=reg165*reg146; T reg241=reg57*reg188; T reg242=reg72*reg201; T reg243=reg165*reg174;
    T reg244=reg57*reg171; T reg245=reg166*reg173; T reg246=reg166*reg174; T reg247=reg165*reg178; T reg248=reg57*reg159;
    T reg249=reg166*reg177; T reg250=reg72*reg202; T reg251=reg165*reg182; T reg252=reg57*reg195; T reg253=reg162*reg184;
    T reg254=reg165*reg185; T reg255=reg57*reg191; reg157=reg72*reg157; T reg256=reg166*reg146; T reg257=reg152*reg205;
    T reg258=reg23*reg202; T reg259=reg72*reg197; T reg260=reg164*reg178; T reg261=reg23*reg201; T reg262=reg72*reg203;
    T reg263=reg166*reg182; T reg264=reg23*reg198; T reg265=reg164*reg182; T reg266=reg72*reg198; reg177=reg165*reg177;
    reg205=reg72*reg205; T reg267=reg57*reg193; T reg268=reg165*reg184; T reg269=reg72*reg207; reg194=reg57*reg194;
    reg269=reg245+reg269; reg232=reg231+reg232; reg242=reg239+reg242; reg224=reg223+reg224; reg266=reg263+reg266;
    reg180=reg175-reg180; reg250=reg246+reg250; reg253=reg257+reg253; reg220=reg219+reg220; reg235=reg234-reg235;
    reg228=reg227+reg228; reg256=reg262+reg256; reg259=reg187+reg259; reg212=reg209-reg212; reg215=reg216+reg215;
    reg251=reg252+reg251; reg155=reg258-reg155; reg233=reg237-reg233; reg236=reg238-reg236; reg243=reg244+reg243;
    reg260=reg261-reg260; reg229=reg230-reg229; reg210=reg192-reg210; reg265=reg264-reg265; reg186=reg206-reg186;
    reg225=reg226-reg225; reg177=reg194+reg177; reg254=reg255+reg254; reg205=reg211+reg205; reg221=reg222-reg221;
    reg240=reg241+reg240; reg157=reg249+reg157; reg247=reg248+reg247; reg217=reg218+reg217; reg268=reg267+reg268;
    reg213=reg214-reg213; reg210=reg24*reg210; reg247=reg24*reg247; reg155=reg24*reg155; reg212=reg24*reg212;
    reg260=reg24*reg260; reg259=reg24*reg259; reg242=reg24*reg242; reg256=reg24*reg256; reg243=reg24*reg243;
    reg265=reg24*reg265; reg205=reg24*reg205; reg177=reg24*reg177; reg266=reg24*reg266; reg269=reg24*reg269;
    reg268=reg24*reg268; reg240=reg24*reg240; reg217=reg24*reg217; reg215=reg24*reg215; reg157=reg24*reg157;
    reg221=reg24*reg221; reg220=reg24*reg220; reg225=reg24*reg225; reg224=reg24*reg224; reg254=reg24*reg254;
    reg253=reg24*reg253; reg229=reg24*reg229; reg228=reg24*reg228; reg250=reg24*reg250; reg213=reg24*reg213;
    reg180=reg24*reg180; reg186=reg24*reg186; reg251=reg24*reg251; reg235=reg24*reg235; reg236=reg24*reg236;
    reg232=reg24*reg232; reg233=reg24*reg233; T tmp_1_4=ponderation*reg243; T tmp_0_3=ponderation*reg242; T tmp_1_1=ponderation*reg254;
    T tmp_0_6=ponderation*reg205; T tmp_0_7=ponderation*reg157; T tmp_4_6=ponderation*reg253; T tmp_0_4=ponderation*reg250; T tmp_1_3=ponderation*reg247;
    T tmp_1_2=ponderation*reg251; T tmp_4_4=ponderation*reg217; T tmp_4_5=ponderation*reg215; T tmp_3_7=ponderation*reg221; T tmp_4_7=ponderation*reg220;
    T tmp_3_6=ponderation*reg225; T tmp_5_5=ponderation*reg224; T tmp_3_5=ponderation*reg229; T tmp_5_6=ponderation*reg228; T tmp_3_4=ponderation*reg233;
    T tmp_5_7=ponderation*reg232; T tmp_3_3=ponderation*reg236; T tmp_6_6=ponderation*reg235; T tmp_2_7=ponderation*reg186; T tmp_6_7=ponderation*reg180;
    T tmp_2_6=ponderation*reg213; T tmp_7_7=ponderation*reg212; T tmp_2_5=ponderation*reg210; T tmp_2_4=ponderation*reg155; T tmp_0_5=ponderation*reg256;
    T tmp_2_3=ponderation*reg260; T tmp_0_1=ponderation*reg259; T tmp_2_2=ponderation*reg265; T tmp_1_7=ponderation*reg177; T tmp_0_0=ponderation*reg269;
    T tmp_1_6=ponderation*reg268; T tmp_0_2=ponderation*reg266; T tmp_1_5=ponderation*reg240;
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
    T reg0=2*(*f.m).shear_modulus_13; T reg1=2*(*f.m).shear_modulus_23; reg1=1.0/reg1; reg0=1.0/reg0; T reg2=2*(*f.m).shear_modulus_12;
    reg2=1.0/reg2; T reg3=reg0*reg1; T reg4=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg5=1.0/(*f.m).elastic_modulus_3; T reg6=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg7=pow((*f.m).v2[0],2); T reg8=pow((*f.m).v2[1],2); T reg9=pow((*f.m).v1[1],2); T reg10=pow((*f.m).v1[0],2); T reg11=reg2*reg3;
    reg8=reg7+reg8; reg7=pow((*f.m).v2[2],2); T reg12=pow((*f.m).v1[2],2); T reg13=reg4*reg11; T reg14=reg5*reg11;
    T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg17=reg6*reg11; reg9=reg10+reg9; reg12=reg9+reg12;
    reg9=reg4*reg13; reg10=reg14*reg16; reg7=reg8+reg7; reg8=reg14*reg15; T reg18=reg4*reg17;
    reg9=reg8-reg9; reg12=pow(reg12,0.5); reg8=1.0/(*f.m).elastic_modulus_1; reg7=pow(reg7,0.5); T reg19=reg15*reg17;
    T reg20=reg16*reg13; reg18=reg10+reg18; T reg21=reg8*reg9; T reg22=(*f.m).v1[1]/reg12; T reg23=(*f.m).v1[2]/reg12;
    T reg24=(*f.m).v2[1]/reg7; T reg25=(*f.m).v2[2]/reg7; T reg26=reg16*reg18; T reg27=reg20+reg19; T reg28=reg4*reg3;
    reg12=(*f.m).v1[0]/reg12; reg7=(*f.m).v2[0]/reg7; T reg29=reg23*reg24; T reg30=reg22*reg25; T reg31=reg2*reg1;
    T reg32=reg5*reg3; T reg33=reg16*reg11; T reg34=reg6*reg17; reg14=reg14*reg8; T reg35=reg6*reg13;
    reg11=reg15*reg11; T reg36=reg6*reg27; reg26=reg21-reg26; reg3=reg6*reg3; reg35=reg10+reg35;
    reg10=reg6*reg11; reg21=reg12*reg25; T reg37=reg4*reg3; T reg38=reg23*reg7; T reg39=reg5*reg31;
    T reg40=reg4*reg28; T reg41=reg16*reg32; T reg42=reg30-reg29; T reg43=2*reg12; T reg44=reg4*reg31;
    reg17=reg16*reg17; reg36=reg26-reg36; reg31=reg6*reg31; reg26=2*reg7; reg13=reg8*reg13;
    reg34=reg14-reg34; reg32=reg15*reg32; reg14=reg6*reg33; T reg45=reg2*reg0; T reg46=2*reg42;
    T reg47=reg22*reg7; T reg48=reg4*reg44; T reg49=reg12*reg24; reg18=reg18/reg36; T reg50=reg38-reg21;
    reg34=reg34/reg36; T reg51=reg15*reg39; T reg52=reg24*reg26; T reg53=pow(reg12,2); T reg54=pow(reg22,2);
    T reg55=reg22*reg43; reg39=reg16*reg39; reg14=reg13+reg14; T reg56=pow(reg7,2); T reg57=pow(reg24,2);
    reg17=reg13+reg17; reg9=reg9/reg36; reg3=reg15*reg3; reg13=reg4*reg45; T reg58=reg5*reg45;
    reg11=reg8*reg11; reg45=reg6*reg45; reg33=reg16*reg33; reg40=reg32-reg40; reg32=reg4*reg31;
    reg28=reg16*reg28; reg37=reg41+reg37; reg35=reg35/reg36; reg10=reg20+reg10; reg41=reg57*reg18;
    T reg59=reg55*reg9; T reg60=reg54*reg9; T reg61=reg52*reg18; reg48=reg51-reg48; reg51=reg4*reg45;
    T reg62=reg4*reg13; T reg63=reg16*reg58; reg58=reg15*reg58; reg27=reg27/reg36; T reg64=reg53*reg35;
    T reg65=reg56*reg34; reg32=reg39+reg32; reg39=reg54*reg35; reg10=reg10/reg36; T reg66=reg57*reg34;
    reg33=reg11-reg33; reg11=reg52*reg34; T reg67=reg55*reg35; reg44=reg16*reg44; T reg68=reg53*reg9;
    reg17=reg17/reg36; T reg69=reg56*reg18; reg31=reg15*reg31; reg40=reg8*reg40; reg37=reg16*reg37;
    T reg70=reg28+reg3; T reg71=reg49-reg47; T reg72=reg50*reg46; reg14=reg14/reg36; T reg73=pow(reg25,2);
    T reg74=pow(reg23,2); T reg75=pow(reg42,2); T reg76=pow(reg50,2); T reg77=reg54*reg10; reg62=reg58-reg62;
    reg58=reg55*reg10; T reg78=reg52*reg14; T reg79=reg72*reg27; reg51=reg63+reg51; reg33=reg33/reg36;
    reg63=reg75*reg27; reg69=reg68+reg69; reg68=pow(reg71,2); reg11=reg67+reg11; reg41=reg60+reg41;
    reg60=reg76*reg27; reg67=reg74*reg9; T reg80=reg73*reg18; T reg81=reg72*reg17; reg61=reg59+reg61;
    reg59=reg76*reg17; reg48=reg8*reg48; T reg82=reg73*reg34; T reg83=reg44+reg31; reg13=reg16*reg13;
    T reg84=reg56*reg14; reg65=reg64+reg65; reg64=reg57*reg14; reg66=reg39+reg66; reg70=reg6*reg70;
    reg45=reg15*reg45; reg39=reg75*reg17; T reg85=reg53*reg10; reg32=reg16*reg32; T reg86=reg74*reg35;
    reg37=reg40-reg37; reg80=reg67+reg80; reg40=reg68*reg27; reg70=reg37-reg70; reg60=reg41+reg60;
    reg32=reg48-reg32; reg63=reg69+reg63; reg37=reg13+reg45; reg83=reg6*reg83; reg51=reg16*reg51;
    reg79=reg61+reg79; reg62=reg8*reg62; reg39=reg65+reg39; reg84=reg85+reg84; reg41=reg75*reg33;
    reg64=reg77+reg64; reg48=reg76*reg33; reg59=reg66+reg59; reg61=reg74*reg10; reg65=reg73*reg14;
    reg66=reg12*reg22; reg67=reg7*reg24; reg81=reg11+reg81; reg11=reg68*reg17; reg82=reg86+reg82;
    reg69=reg23*reg43; reg77=2*reg22; reg78=reg58+reg78; reg58=reg72*reg33; reg85=2*reg24;
    reg86=reg25*reg26; T reg87=reg25*reg85; T reg88=reg57*reg81; T reg89=reg54*reg79; T reg90=reg86*reg34;
    T reg91=reg23*reg77; T reg92=reg69*reg35; T reg93=reg56*reg39; T reg94=reg53*reg63; reg11=reg82+reg11;
    reg40=reg80+reg40; reg80=reg67*reg81; reg82=reg66*reg79; reg58=reg78+reg58; reg37=reg6*reg37;
    reg81=reg56*reg81; reg51=reg62-reg51; reg79=reg53*reg79; reg62=reg68*reg33; reg65=reg61+reg65;
    reg48=reg64+reg48; reg83=reg32-reg83; reg32=reg67*reg59; reg61=reg66*reg60; reg41=reg84+reg41;
    reg64=reg56*reg59; reg70=reg70/reg36; reg78=reg53*reg60; reg84=reg50*reg42; T reg95=reg57*reg39;
    T reg96=2*reg23; T reg97=reg54*reg63; T reg98=reg24*reg77; T reg99=reg7*reg43; T reg100=reg12*reg50;
    T reg101=reg22*reg42; T reg102=reg86*reg18; reg46=reg71*reg46; reg60=reg54*reg60; reg49=reg47+reg49;
    reg59=reg57*reg59; reg47=2*reg50; T reg103=reg69*reg9; T reg104=reg12*reg7; T reg105=reg22*reg24;
    reg96=reg25*reg96; reg18=reg87*reg18; T reg106=reg84*reg48; reg32=reg61+reg32; reg61=reg104*reg70;
    reg39=reg67*reg39; reg83=reg83/reg36; T reg107=reg56*reg11; T reg108=reg53*reg40; reg9=reg91*reg9;
    reg62=reg65+reg62; reg95=reg97+reg95; reg65=reg49*reg70; reg81=reg79+reg81; reg79=reg75*reg58;
    reg97=reg99*reg16; reg37=reg51-reg37; reg51=reg46*reg27; reg102=reg103+reg102; reg103=reg69*reg10;
    T reg109=reg98*reg15; reg80=reg82+reg80; reg82=reg84*reg58; T reg110=reg86*reg14; T reg111=reg105*reg70;
    reg58=reg76*reg58; reg88=reg89+reg88; reg89=reg57*reg15; T reg112=reg46*reg17; reg90=reg92+reg90;
    reg92=reg57*reg16; T reg113=reg75*reg41; reg93=reg94+reg93; reg100=reg101+reg100; reg94=reg22*reg50;
    reg101=reg56*reg16; T reg114=reg12*reg42; T reg115=reg54*reg40; T reg116=reg57*reg11; T reg117=reg23*reg25;
    reg47=reg71*reg47; reg63=reg66*reg63; T reg118=reg98*reg16; reg59=reg60+reg59; reg60=reg76*reg41;
    T reg119=reg56*reg8; T reg120=reg75*reg48; reg64=reg78+reg64; reg48=reg76*reg48; reg78=reg7*reg50;
    T reg121=reg24*reg42; T reg122=reg99*reg8; reg34=reg87*reg34; reg35=reg91*reg35; reg51=reg102+reg51;
    reg27=reg47*reg27; reg18=reg9+reg18; reg17=reg47*reg17; reg34=reg35+reg34; reg112=reg90+reg112;
    reg110=reg103+reg110; reg9=reg46*reg33; reg10=reg91*reg10; reg14=reg87*reg14; reg35=reg99*reg6;
    reg116=reg115+reg116; reg90=reg98*reg4; reg102=reg76*reg62; reg113=reg93+reg113; reg93=reg99*reg61;
    reg58=reg88+reg58; reg88=reg98*reg65; reg92=reg119-reg92; reg120=reg64+reg120; reg64=reg99*reg111;
    reg39=reg63+reg39; reg41=reg84*reg41; reg63=reg96*reg6; reg107=reg108+reg107; reg103=reg75*reg62;
    reg106=reg32+reg106; reg32=reg49*reg111; reg118=reg122-reg118; reg40=reg66*reg40; reg11=reg67*reg11;
    reg79=reg81+reg79; reg81=reg99*reg65; reg82=reg80+reg82; reg65=reg49*reg65; reg36=reg37/reg36;
    reg101=reg89-reg101; reg37=reg96*reg4; reg97=reg109-reg97; reg80=reg117*reg70; reg60=reg95+reg60;
    reg89=reg98*reg61; reg95=reg114*reg83; reg108=reg94*reg83; reg77=reg50*reg77; reg43=reg42*reg43;
    reg48=reg59+reg48; reg111=reg98*reg111; reg78=reg121+reg78; reg59=reg100*reg83; reg109=reg24*reg50;
    reg115=reg7*reg42; reg119=reg73*reg6; reg121=reg57*reg4; reg122=reg23*reg71; T reg123=reg56*reg6;
    T reg124=reg73*reg4; T reg125=reg56*reg112; T reg126=reg53*reg51; reg33=reg47*reg33; reg15=reg54*reg15;
    T reg127=reg43*reg59; reg81=reg79+reg81; reg8=reg53*reg8; reg79=reg53*reg16; reg63=reg118-reg63;
    reg118=reg99*reg80; reg103=reg107+reg103; reg107=reg78*reg36; T reg128=reg109*reg36; reg93=reg113+reg93;
    reg113=reg122*reg83; T reg129=reg115*reg36; T reg130=reg43*reg95; T reg131=reg43*reg108; reg64=reg120+reg64;
    reg119=reg92-reg119; reg121=reg123+reg121; reg96=reg96*reg5; reg90=reg35+reg90; reg102=reg116+reg102;
    reg35=reg98*reg80; reg92=(*f.m).alpha_1*reg53; reg88=reg58+reg88; reg58=reg77*reg59; reg116=reg54*reg51;
    reg120=reg57*reg112; reg16=reg54*reg16; reg41=reg39+reg41; reg61=reg49*reg61; reg32=reg106+reg32;
    reg39=reg100*reg108; reg106=(*f.m).alpha_2*reg56; reg11=reg40+reg11; reg62=reg84*reg62; reg40=(*f.m).alpha_1*reg54;
    reg123=(*f.m).alpha_2*reg57; reg65=reg82+reg65; reg59=reg100*reg59; reg82=reg77*reg95; reg89=reg60+reg89;
    reg37=reg97-reg37; reg85=reg50*reg85; reg124=reg101-reg124; reg26=reg42*reg26; reg27=reg18+reg27;
    reg17=reg34+reg17; reg18=reg25*reg71; reg111=reg48+reg111; reg34=reg73*reg5; reg38=reg21+reg38;
    reg21=reg23*reg42; reg48=reg12*reg71; reg14=reg10+reg14; reg108=reg77*reg108; reg9=reg110+reg9;
    reg125=reg126+reg125; reg80=reg49*reg80; reg59=reg65+reg59; reg62=reg11+reg62; reg51=reg66*reg51;
    reg10=reg75*reg9; reg11=reg74*reg4; reg60=reg78*reg107; reg65=reg78*reg128; reg39=reg32+reg39;
    reg32=reg53*reg27; reg35=reg102+reg35; reg97=reg85*reg128; reg101=reg77*reg113; reg108=reg111+reg108;
    reg58=reg88+reg58; reg88=reg85*reg107; reg120=reg116+reg120; reg102=reg76*reg9; reg16=reg8-reg16;
    reg8=reg74*reg6; reg110=reg54*reg27; reg111=reg57*reg17; reg116=reg56*reg17; reg126=reg85*reg129;
    reg61=reg41+reg61; reg95=reg100*reg95; reg112=reg67*reg112; reg82=reg89+reg82; reg123=reg40+reg123;
    reg40=reg38*reg70; reg33=reg14+reg33; reg6=reg53*reg6; reg14=reg104*reg119; reg4=reg54*reg4;
    reg41=(*f.m).alpha_3*reg75; reg106=reg92+reg106; reg90=reg96-reg90; reg121=reg34-reg121; reg30=reg29+reg30;
    reg29=reg105*reg37; reg34=reg104*reg63; reg89=reg54*reg124; reg92=reg7*reg71; reg96=reg25*reg42;
    T reg132=reg53*reg119; T reg133=reg22*reg71; T reg134=reg23*reg50; reg48=reg21+reg48; reg21=reg105*reg124;
    T reg135=reg57*reg37; reg107=reg26*reg107; T reg136=reg56*reg63; T reg137=reg56*reg119; reg127=reg81+reg127;
    reg79=reg15-reg79; reg15=reg57*reg124; reg81=(*f.m).alpha_2*reg73; T reg138=reg43*reg113; reg118=reg103+reg118;
    reg103=(*f.m).alpha_1*reg74; reg128=reg26*reg128; reg131=reg64+reg131; reg64=(*f.m).alpha_3*reg76; T reg139=reg26*reg129;
    reg130=reg93+reg130; reg93=reg54*reg37; T reg140=reg18*reg36; T reg141=reg53*reg63; reg112=reg51+reg112;
    reg9=reg84*reg9; reg27=reg66*reg27; reg17=reg67*reg17; reg51=reg76*reg33; reg111=reg110+reg111;
    reg110=reg73*reg121; reg15=reg137+reg15; reg137=reg98*reg40; reg102=reg120+reg102; reg120=(*f.m).alpha_2*reg67;
    T reg142=(*f.m).alpha_1*reg66; T reg143=(*f.m).alpha_3*reg68; reg81=reg103+reg81; reg64=reg123+reg64; reg60=reg59+reg60;
    reg113=reg100*reg113; reg80=reg62+reg80; reg41=reg106+reg41; reg65=reg39+reg65; reg39=reg117*reg90;
    reg29=reg34+reg29; reg34=reg117*reg121; reg21=reg14+reg21; reg129=reg78*reg129; reg95=reg61+reg95;
    reg14=reg73*reg90; reg135=reg136+reg135; reg11=reg79-reg11; reg59=reg85*reg140; reg101=reg35+reg101;
    reg35=reg48*reg83; reg70=reg30*reg70; reg5=reg74*reg5; reg7=reg7*reg25; reg61=reg49*reg2;
    reg67=reg67*reg2; reg97=reg108+reg97; reg133=reg134+reg133; reg92=reg96+reg92; reg12=reg12*reg23;
    reg126=reg82+reg126; reg62=reg25*reg50; reg79=reg24*reg71; reg4=reg6+reg4; reg107=reg127+reg107;
    reg10=reg125+reg10; reg6=reg74*reg90; reg88=reg58+reg88; reg93=reg141+reg93; reg58=reg99*reg40;
    reg82=reg26*reg140; reg138=reg118+reg138; reg128=reg131+reg128; reg96=reg74*reg121; reg89=reg132+reg89;
    reg116=reg32+reg116; reg32=reg75*reg33; reg139=reg130+reg139; reg8=reg16-reg8; reg16=reg49*reg61;
    reg4=reg5-reg4; reg39=reg29+reg39; reg143=reg81+reg143; reg120=reg142+reg120; reg5=(*f.m).alpha_3*reg84;
    reg29=reg56*reg8; reg81=reg55*reg61; reg103=reg57*reg11; reg6=reg93+reg6; reg33=reg84*reg33;
    reg110=reg15+reg110; reg15=reg55*reg67; reg96=reg89+reg96; reg17=reg27+reg17; reg40=reg49*reg40;
    reg9=reg112+reg9; reg27=reg54*reg11; reg84=reg52*reg67; reg89=reg53*reg8; reg93=reg38*reg0;
    reg14=reg135+reg14; reg106=reg52*reg61; reg108=(*f.m).alpha_2*reg7; reg7=reg7*reg0; reg112=(*f.m).alpha_1*reg12;
    reg34=reg21+reg34; reg21=reg49*reg67; reg118=reg107*reg65; reg123=reg126*reg41; reg125=reg97*reg64;
    reg59=reg101+reg59; reg137=reg102+reg137; reg101=reg77*reg35; reg99=reg99*reg70; reg32=reg116+reg32;
    reg102=reg128*reg64; reg116=reg88*reg65; reg127=reg128*reg60; reg130=reg139*reg41; reg131=reg43*reg35;
    reg58=reg10+reg58; reg51=reg111+reg51; reg140=reg78*reg140; reg113=reg80+reg113; reg129=reg95+reg129;
    reg25=reg24*reg25; reg23=reg22*reg23; reg42=reg71*reg42; reg79=reg62+reg79; reg83=reg133*reg83;
    reg10=reg92*reg36; reg22=reg97*reg60; reg82=reg138+reg82; reg98=reg98*reg70; reg24=reg105*reg11;
    reg116=reg22-reg116; reg2=reg66*reg2; reg22=reg104*reg8; reg50=reg71*reg50; reg81=reg6+reg81;
    reg6=reg69*reg7; reg15=reg96+reg15; reg62=reg69*reg93; reg103=reg29+reg103; reg73=reg73*reg4;
    reg42=(*f.m).alpha_3*reg42; reg108=reg112+reg108; reg84=reg110+reg84; reg29=reg86*reg7; reg74=reg74*reg4;
    reg27=reg89+reg27; reg66=reg107*reg97; reg71=reg128*reg88; reg80=reg30*reg1; reg89=reg25*reg1;
    reg118=reg127-reg118; reg102=reg130+reg102; reg106=reg14+reg106; reg14=reg86*reg93; reg140=reg113+reg140;
    reg95=reg65*reg64; reg77=reg77*reg83; reg40=reg9+reg40; reg35=reg100*reg35; reg9=reg129*reg41;
    reg98=reg51+reg98; reg51=reg38*reg93; reg33=reg17+reg33; reg70=reg49*reg70; reg16=reg39+reg16;
    reg59=reg59*reg143; reg125=reg123+reg125; reg17=reg85*reg10; reg101=reg137+reg101; reg43=reg43*reg83;
    reg82=reg82*reg143; reg5=reg120+reg5; reg25=(*f.m).alpha_2*reg25; reg39=(*f.m).alpha_1*reg23; reg21=reg34+reg21;
    reg36=reg79*reg36; reg99=reg32+reg99; reg32=reg26*reg10; reg131=reg58+reg131; reg34=reg38*reg7;
    reg140=reg140*reg143; reg95=reg9+reg95; reg82=reg102+reg82; reg62=reg81+reg62; reg9=reg91*reg80;
    reg34=reg21+reg34; reg21=reg107*reg5; reg73=reg103+reg73; reg58=reg30*reg89; reg81=reg52*reg2;
    reg14=reg106+reg14; reg96=reg30*reg80; reg51=reg16+reg51; reg29=reg84+reg29; reg16=reg87*reg89;
    reg117=reg117*reg4; reg24=reg22+reg24; reg59=reg125+reg59; reg22=reg88*reg5; reg84=reg87*reg80;
    reg25=reg39+reg25; reg42=reg108+reg42; reg32=reg131+reg32; reg43=reg99+reg43; reg26=reg26*reg36;
    reg50=(*f.m).alpha_3*reg50; reg66=reg71-reg66; reg83=reg100*reg83; reg70=reg33+reg70; reg17=reg101+reg17;
    reg33=reg55*reg2; reg10=reg78*reg10; reg74=reg27+reg74; reg6=reg15+reg6; reg15=reg139*reg116;
    reg35=reg40+reg35; reg77=reg98+reg77; reg85=reg85*reg36; reg27=reg126*reg118; reg0=reg12*reg0;
    reg12=reg91*reg89; reg84=reg14+reg84; reg50=reg25+reg50; reg117=reg24+reg117; reg14=reg49*reg2;
    reg1=reg23*reg1; reg21=reg82+reg21; reg32=reg32*reg42; reg58=reg34+reg58; reg26=reg43+reg26;
    reg36=reg78*reg36; reg83=reg70+reg83; reg10=reg35+reg10; reg96=reg51+reg96; reg22=reg59+reg22;
    reg17=reg17*reg42; reg85=reg77+reg85; reg140=reg95+reg140; reg23=reg60*reg5; reg24=reg129*reg66;
    reg81=reg73+reg81; reg25=reg107*reg129; reg86=reg86*reg0; reg9=reg62+reg9; reg69=reg69*reg0;
    reg33=reg74+reg33; reg12=reg6+reg12; reg6=reg126*reg60; reg16=reg29+reg16; reg29=reg88*reg129;
    reg27=reg15-reg27; reg15=reg139*reg60; reg86=reg81+reg86; reg87=reg87*reg1; reg69=reg33+reg69;
    reg36=reg83+reg36; reg91=reg91*reg1; reg33=reg97*reg129; reg34=reg139*reg65; reg17=reg22+reg17;
    reg85=reg85*reg50; reg24=reg27+reg24; reg10=reg10*reg42; reg23=reg140+reg23; reg22=reg107*reg126;
    reg27=reg139*reg88; reg14=reg117+reg14; reg38=reg38*reg0; reg35=reg12*reg96; reg39=reg16*reg96;
    reg40=reg126*reg65; reg43=reg84*reg58; reg51=reg128*reg129; reg25=reg15-reg25; reg32=reg21+reg32;
    reg26=reg26*reg50; reg29=reg6-reg29; reg6=reg9*reg58; reg33=reg40-reg33; reg85=reg17+reg85;
    reg10=reg23+reg10; reg36=reg36*reg50; reg118=reg118/reg24; reg29=reg29/reg24; reg116=reg116/reg24;
    reg25=reg25/reg24; reg26=reg32+reg26; reg51=reg34-reg51; reg6=reg35-reg6; reg15=reg12*reg84;
    reg17=reg9*reg16; reg21=reg139*reg97; reg22=reg27-reg22; reg23=reg128*reg126; reg30=reg30*reg1;
    reg38=reg14+reg38; reg87=reg86+reg87; reg91=reg69+reg91; reg43=reg39-reg43; reg14=1-var_inter[1];
    reg27=1-var_inter[0]; reg33=reg33/reg24; reg36=reg10+reg36; reg116=reg116*reg26; reg118=reg118*reg85;
    reg29=reg29*reg26; reg25=reg25*reg85; reg10=elem.pos(1)[1]*var_inter[0]; reg32=reg27*elem.pos(0)[1]; reg51=reg51/reg24;
    reg66=reg66/reg24; reg34=elem.pos(1)[0]*var_inter[0]; reg22=reg22/reg24; reg35=elem.pos(0)[0]*reg27; reg23=reg21-reg23;
    reg21=reg14*elem.pos(1)[1]; reg39=reg14*elem.pos(0)[1]; reg40=reg14*elem.pos(1)[0]; reg59=reg14*elem.pos(0)[0]; reg30=reg38+reg30;
    reg38=reg87*reg6; reg62=reg43*reg91; reg17=reg15-reg17; reg15=reg16*reg30; reg69=reg87*reg58;
    reg70=reg32+reg10; reg71=reg30*reg17; reg73=elem.pos(2)[1]*var_inter[0]; reg38=reg62-reg38; reg62=reg35+reg34;
    reg74=elem.pos(2)[0]*var_inter[0]; reg24=reg23/reg24; reg21=reg21-reg39; reg23=elem.pos(2)[1]*var_inter[1]; reg40=reg40-reg59;
    reg77=elem.pos(2)[0]*var_inter[1]; reg118=reg116-reg118; reg81=reg12*reg30; reg85=reg51*reg85; reg51=reg91*reg58;
    reg26=reg33*reg26; reg66=reg66*reg36; reg29=reg25-reg29; reg22=reg22*reg36; reg25=reg9*reg30;
    reg81=reg51-reg81; reg33=elem.pos(3)[0]*reg27; reg85=reg26-reg85; reg74=reg74-reg62; reg26=reg12*reg87;
    reg71=reg38+reg71; reg22=reg29-reg22; reg29=reg91*reg16; reg21=reg23+reg21; reg23=elem.pos(3)[1]*var_inter[1];
    reg15=reg69-reg15; reg38=reg91*reg96; reg51=reg84*reg30; reg36=reg24*reg36; reg40=reg77+reg40;
    reg24=1-(*f.m).resolution; reg69=elem.pos(3)[1]*reg27; reg73=reg73-reg70; reg77=reg87*reg96; reg118=reg66+reg118;
    reg66=elem.pos(3)[0]*var_inter[1]; reg82=reg91*reg84; reg83=reg9*reg87; reg40=reg40-reg66; reg81=reg81/reg71;
    reg86=(*f.m).resolution*reg64; reg25=reg38-reg25; reg15=reg15/reg71; reg22=reg24*reg22; reg51=reg77-reg51;
    reg69=reg73+reg69; reg118=reg24*reg118; reg26=reg29-reg26; reg33=reg74+reg33; reg85=reg36+reg85;
    reg21=reg21-reg23; reg29=(*f.m).resolution*reg41; reg36=(*f.m).resolution*reg81; reg85=reg24*reg85; reg118=reg29+reg118;
    reg22=reg86+reg22; reg29=(*f.m).resolution*reg15; reg38=(*f.m).resolution*reg5; reg129=reg129*reg24; reg65=reg65*reg24;
    reg51=reg51/reg71; reg73=reg21*reg33; reg74=reg40*reg69; reg26=reg26/reg71; reg6=reg6/reg71;
    reg25=reg25/reg71; reg83=reg82-reg83; reg43=reg43/reg71; reg17=reg17/reg71; reg71=reg83/reg71;
    reg73=reg74-reg73; reg139=reg139*reg24; reg128=reg128*reg24; reg126=reg126*reg24; reg97=reg97*reg24;
    reg60=reg60*reg24; reg129=reg29+reg129; reg36=reg65-reg36; reg29=(*f.m).resolution*reg6; reg65=(*f.m).resolution*reg26;
    reg74=(*f.m).resolution*reg43; reg77=(*f.m).resolution*reg25; reg85=reg38+reg85; reg38=(*f.m).resolution*reg51; reg118=(*f.m).deltaT*reg118;
    reg22=(*f.m).deltaT*reg22; reg82=reg36*reg22; reg69=reg69/reg73; reg97=reg77+reg97; reg38=reg126-reg38;
    reg77=reg129*reg118; reg29=reg128-reg29; reg139=reg74+reg139; reg85=(*f.m).deltaT*reg85; reg40=reg40/reg73;
    reg88=reg88*reg24; reg24=reg107*reg24; reg74=(*f.m).resolution*reg71; reg60=reg65+reg60; reg65=(*f.m).resolution*reg17;
    reg21=reg21/reg73; reg33=reg33/reg73; reg83=reg60*reg85; reg86=reg14*reg69; reg95=reg77+reg82;
    reg98=reg27*reg21; reg74=reg88-reg74; reg24=reg65+reg24; reg65=var_inter[1]*reg33; reg88=reg139*reg118;
    reg99=var_inter[1]*reg69; reg101=reg29*reg22; reg102=reg27*reg40; reg103=reg14*reg33; reg106=var_inter[0]*reg40;
    reg107=reg97*reg22; reg108=var_inter[0]*reg21; reg110=reg38*reg118; reg111=reg74*reg85; reg112=reg110+reg107;
    reg113=reg102+reg65; reg116=reg108+reg86; reg117=reg24*reg85; reg120=reg103+reg106; reg123=reg88+reg101;
    reg125=reg98+reg99; reg126=reg95+reg83; reg127=reg123+reg117; reg128=2*reg126; reg130=reg112+reg111;
    reg131=reg106-reg65; reg132=reg99-reg108; reg134=0.5*reg116; reg135=0.5*reg120; reg136=reg98-reg86;
    reg137=0.5*reg125; reg138=0.5*reg113; reg140=reg103-reg102; reg141=reg134*reg128; reg142=reg120*reg130;
    T reg144=reg135*reg128; T reg145=reg116*reg127; T reg146=reg137*reg128; T reg147=0.5*reg140; T reg148=reg113*reg130;
    T reg149=reg138*reg128; T reg150=0.5*reg136; T reg151=reg125*reg127; T reg152=0.5*reg131; T reg153=0.5*reg132;
    T reg154=reg131*reg130; T reg155=reg153*reg128; T reg156=reg152*reg128; T reg157=reg132*reg127; T reg158=reg146-reg148;
    T reg159=reg142-reg141; T reg160=reg151-reg149; T reg161=reg144-reg145; T reg162=reg147*reg128; T reg163=reg150*reg128;
    T reg164=reg136*reg127; T reg165=reg140*reg130; reg158=reg73*reg158; reg161=reg73*reg161; reg160=reg73*reg160;
    reg159=reg73*reg159; T reg166=reg165+reg163; T reg167=reg157+reg156; T reg168=reg162+reg164; T reg169=reg154+reg155;
    reg160=ponderation*reg160; T reg170=reg73*reg169; T reg171=reg73*reg168; T reg172=reg73*reg167; T reg173=reg73*reg166;
    reg159=ponderation*reg159; reg158=ponderation*reg158; reg161=ponderation*reg161; sollicitation[indices[3]+1]+=-reg158; sollicitation[indices[3]+0]+=-reg160;
    reg158=ponderation*reg170; sollicitation[indices[2]+1]+=reg158; reg160=ponderation*reg172; sollicitation[indices[2]+0]+=reg160; sollicitation[indices[1]+1]+=-reg159;
    sollicitation[indices[1]+0]+=-reg161; reg159=ponderation*reg173; sollicitation[indices[0]+1]+=reg159; reg161=ponderation*reg171; sollicitation[indices[0]+0]+=reg161;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; T reg2=2*(*f.m).shear_modulus_12; reg1=1.0/reg1;
    T reg3=reg1*reg0; reg2=1.0/reg2; T reg4=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg5=1.0/(*f.m).elastic_modulus_3; T reg6=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg7=pow((*f.m).v1[0],2); T reg8=pow((*f.m).v1[1],2); T reg9=reg2*reg3; T reg10=pow((*f.m).v2[1],2); T reg11=pow((*f.m).v2[0],2);
    reg8=reg7+reg8; reg7=reg6*reg9; T reg12=1.0/(*f.m).elastic_modulus_2; T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg14=reg5*reg9;
    T reg15=reg4*reg9; T reg16=pow((*f.m).v2[2],2); reg10=reg11+reg10; reg11=pow((*f.m).v1[2],2); reg16=reg10+reg16;
    reg11=reg8+reg11; reg8=reg4*reg7; reg10=reg14*reg12; T reg17=reg4*reg15; T reg18=reg14*reg13;
    T reg19=reg12*reg7; reg11=pow(reg11,0.5); T reg20=reg13*reg15; reg8=reg18+reg8; reg17=reg10-reg17;
    reg10=1.0/(*f.m).elastic_modulus_1; reg16=pow(reg16,0.5); T reg21=reg10*reg17; T reg22=reg13*reg8; T reg23=reg20+reg19;
    T reg24=(*f.m).v1[1]/reg11; T reg25=(*f.m).v1[2]/reg11; T reg26=(*f.m).v2[1]/reg16; T reg27=(*f.m).v2[2]/reg16; T reg28=reg25*reg26;
    T reg29=reg24*reg27; reg11=(*f.m).v1[0]/reg11; reg16=(*f.m).v2[0]/reg16; T reg30=reg6*reg3; T reg31=reg4*reg3;
    T reg32=reg2*reg0; reg3=reg5*reg3; reg22=reg21-reg22; reg21=reg6*reg23; T reg33=reg12*reg9;
    T reg34=reg6*reg15; reg14=reg14*reg10; T reg35=reg6*reg7; reg9=reg13*reg9; reg34=reg18+reg34;
    reg18=reg6*reg33; reg21=reg22-reg21; reg22=reg5*reg32; T reg36=reg11*reg27; T reg37=reg25*reg16;
    T reg38=reg4*reg30; reg35=reg14-reg35; reg14=reg12*reg3; T reg39=reg6*reg9; reg15=reg10*reg15;
    T reg40=reg4*reg32; T reg41=reg29-reg28; reg7=reg13*reg7; T reg42=reg4*reg31; reg3=reg13*reg3;
    T reg43=2*reg11; T reg44=reg2*reg1; reg32=reg6*reg32; T reg45=2*reg16; T reg46=2*reg41;
    reg34=reg34/reg21; reg30=reg12*reg30; reg18=reg20+reg18; reg8=reg8/reg21; T reg47=reg4*reg44;
    reg35=reg35/reg21; reg31=reg13*reg31; reg7=reg15+reg7; reg9=reg13*reg9; T reg48=reg26*reg45;
    reg17=reg17/reg21; reg33=reg10*reg33; T reg49=reg24*reg16; T reg50=pow(reg26,2); T reg51=reg11*reg26;
    T reg52=pow(reg16,2); T reg53=reg24*reg43; T reg54=reg37-reg36; T reg55=pow(reg24,2); T reg56=pow(reg11,2);
    reg39=reg15+reg39; reg42=reg14-reg42; reg38=reg3+reg38; reg3=reg6*reg44; reg14=reg12*reg22;
    reg44=reg5*reg44; reg22=reg13*reg22; reg15=reg4*reg32; T reg57=reg4*reg40; reg39=reg39/reg21;
    reg7=reg7/reg21; reg9=reg33-reg9; reg33=reg56*reg34; T reg58=reg52*reg35; T reg59=reg55*reg34;
    T reg60=reg50*reg35; reg40=reg13*reg40; T reg61=reg54*reg46; T reg62=pow(reg54,2); T reg63=reg48*reg8;
    T reg64=reg53*reg17; reg38=reg13*reg38; T reg65=reg31+reg30; reg57=reg14-reg57; reg15=reg22+reg15;
    reg14=reg12*reg44; reg44=reg13*reg44; reg18=reg18/reg21; reg22=reg4*reg47; reg23=reg23/reg21;
    T reg66=reg4*reg3; T reg67=reg52*reg8; reg32=reg12*reg32; T reg68=pow(reg27,2); T reg69=reg55*reg17;
    T reg70=reg50*reg8; T reg71=reg56*reg17; reg42=reg10*reg42; T reg72=reg51-reg49; T reg73=pow(reg25,2);
    T reg74=pow(reg41,2); T reg75=reg48*reg35; T reg76=reg53*reg34; T reg77=reg74*reg7; reg58=reg33+reg58;
    reg47=reg13*reg47; reg63=reg64+reg63; reg33=reg68*reg8; reg64=reg48*reg39; reg57=reg10*reg57;
    T reg78=reg61*reg23; T reg79=reg53*reg18; reg9=reg9/reg21; T reg80=reg74*reg23; reg67=reg71+reg67;
    reg70=reg69+reg70; reg65=reg6*reg65; reg69=reg62*reg23; reg38=reg42-reg38; reg42=reg73*reg17;
    reg71=reg73*reg34; T reg81=reg52*reg39; reg75=reg76+reg75; reg66=reg44+reg66; reg44=reg56*reg18;
    reg76=reg62*reg7; reg60=reg59+reg60; reg59=reg61*reg7; T reg82=reg68*reg35; reg22=reg14-reg22;
    reg14=reg40+reg32; T reg83=pow(reg72,2); reg15=reg13*reg15; T reg84=reg55*reg18; reg3=reg12*reg3;
    T reg85=reg50*reg39; reg66=reg13*reg66; reg64=reg79+reg64; reg33=reg42+reg33; reg42=reg83*reg23;
    reg79=reg47+reg3; T reg86=reg61*reg9; reg22=reg10*reg22; reg69=reg70+reg69; reg65=reg38-reg65;
    reg38=reg83*reg7; reg80=reg67+reg80; reg59=reg75+reg59; reg14=reg6*reg14; reg15=reg57-reg15;
    reg81=reg44+reg81; reg44=reg73*reg18; reg57=reg68*reg39; reg76=reg60+reg76; reg60=reg74*reg9;
    reg67=2*reg24; reg70=reg62*reg9; reg85=reg84+reg85; reg75=reg25*reg43; reg84=reg27*reg45;
    reg77=reg58+reg77; reg78=reg63+reg78; reg58=reg11*reg24; reg63=reg16*reg26; T reg87=2*reg26;
    reg82=reg71+reg82; reg71=reg75*reg34; T reg88=reg52*reg59; T reg89=reg56*reg80; reg14=reg15-reg14;
    reg15=reg52*reg77; T reg90=reg24*reg41; reg60=reg81+reg60; reg51=reg49+reg51; reg49=reg24*reg26;
    reg81=reg56*reg78; T reg91=reg56*reg69; T reg92=reg52*reg76; T reg93=reg11*reg16; reg66=reg22-reg66;
    reg22=2*reg54; reg46=reg72*reg46; reg79=reg6*reg79; T reg94=reg27*reg87; reg42=reg33+reg42;
    reg33=reg55*reg80; T reg95=reg54*reg41; T reg96=reg50*reg77; T reg97=reg75*reg17; reg86=reg64+reg86;
    reg64=reg84*reg8; T reg98=reg50*reg59; reg65=reg65/reg21; T reg99=reg58*reg69; T reg100=reg63*reg76;
    T reg101=reg55*reg78; T reg102=reg83*reg9; reg57=reg44+reg57; reg44=reg25*reg67; reg70=reg85+reg70;
    reg38=reg82+reg38; reg69=reg55*reg69; reg76=reg50*reg76; reg82=reg11*reg54; reg85=reg84*reg35;
    reg78=reg58*reg78; reg59=reg63*reg59; reg79=reg66-reg79; reg14=reg14/reg21; reg98=reg101+reg98;
    reg66=reg62*reg86; reg34=reg44*reg34; reg35=reg94*reg35; reg80=reg58*reg80; reg77=reg63*reg77;
    reg102=reg57+reg102; reg57=reg75*reg18; reg101=reg84*reg39; reg100=reg99+reg100; reg99=reg95*reg70;
    T reg103=reg93*reg65; T reg104=reg49*reg65; T reg105=reg51*reg65; reg59=reg78+reg59; reg78=reg95*reg86;
    reg86=reg74*reg86; reg88=reg81+reg88; reg15=reg89+reg15; reg81=reg74*reg60; reg92=reg91+reg92;
    reg89=reg74*reg70; reg91=reg56*reg42; T reg106=reg52*reg38; reg8=reg94*reg8; reg17=reg44*reg17;
    T reg107=reg46*reg23; T reg108=reg62*reg60; reg96=reg33+reg96; reg33=reg16*reg54; T reg109=reg26*reg41;
    T reg110=reg16*reg43; reg64=reg97+reg64; reg76=reg69+reg76; reg82=reg90+reg82; reg69=reg46*reg7;
    reg85=reg71+reg85; reg71=reg50*reg38; reg90=reg55*reg42; reg22=reg72*reg22; reg97=2*reg25;
    T reg111=reg26*reg67; T reg112=reg25*reg27; T reg113=reg11*reg41; T reg114=reg24*reg54; reg70=reg62*reg70;
    reg89=reg92+reg89; reg92=reg110*reg104; T reg115=reg110*reg103; reg81=reg15+reg81; reg106=reg91+reg106;
    reg15=reg74*reg102; reg107=reg64+reg107; reg86=reg88+reg86; reg64=reg110*reg105; reg88=reg82*reg14;
    reg91=reg114*reg14; T reg116=reg113*reg14; reg8=reg17+reg8; reg17=reg112*reg65; reg39=reg94*reg39;
    reg18=reg44*reg18; T reg117=reg46*reg9; reg23=reg22*reg23; reg101=reg57+reg101; reg7=reg22*reg7;
    reg35=reg34+reg35; reg69=reg85+reg69; reg33=reg109+reg33; reg42=reg58*reg42; reg34=reg26*reg54;
    reg57=reg16*reg41; reg85=reg51*reg104; reg109=reg25*reg72; reg70=reg76+reg70; reg104=reg111*reg104;
    reg99=reg100+reg99; reg71=reg90+reg71; reg76=reg62*reg102; reg60=reg95*reg60; reg77=reg80+reg77;
    reg80=reg110*reg10; reg66=reg98+reg66; reg90=reg111*reg105; reg98=reg52*reg13; reg100=reg111*reg13;
    T reg118=reg50*reg12; T reg119=reg110*reg13; reg21=reg79/reg21; reg79=reg111*reg103; reg108=reg96+reg108;
    reg43=reg41*reg43; reg96=reg111*reg12; reg78=reg59+reg78; reg59=reg50*reg13; reg38=reg63*reg38;
    reg67=reg54*reg67; reg97=reg27*reg97; reg105=reg51*reg105; T reg120=reg52*reg10; reg64=reg86+reg64;
    reg105=reg78+reg105; reg115=reg81+reg115; reg78=reg57*reg21; reg81=reg34*reg21; reg86=reg33*reg21;
    T reg121=reg68*reg4; T reg122=reg52*reg6; T reg123=reg43*reg116; reg98=reg118-reg98; reg60=reg77+reg60;
    reg102=reg95*reg102; reg77=(*f.m).alpha_2*reg52; reg103=reg51*reg103; reg118=(*f.m).alpha_1*reg55; T reg124=(*f.m).alpha_2*reg50;
    T reg125=reg109*reg14; reg117=reg101+reg117; reg38=reg42+reg38; reg85=reg99+reg85; reg39=reg18+reg39;
    reg9=reg22*reg9; reg18=reg82*reg91; reg42=reg43*reg88; reg99=reg82*reg88; reg37=reg36+reg37;
    reg36=reg25*reg41; reg101=reg11*reg72; T reg126=reg43*reg91; T reg127=reg110*reg17; reg23=reg8+reg23;
    reg8=reg27*reg72; T reg128=reg67*reg116; reg79=reg108+reg79; reg45=reg41*reg45; reg87=reg54*reg87;
    reg108=reg97*reg4; reg119=reg96-reg119; reg15=reg106+reg15; reg96=reg68*reg6; reg106=reg111*reg4;
    T reg129=reg110*reg6; T reg130=reg52*reg69; T reg131=reg50*reg4; T reg132=reg56*reg107; reg100=reg80-reg100;
    reg80=reg50*reg69; T reg133=reg55*reg107; T reg134=reg97*reg6; reg7=reg35+reg7; reg88=reg67*reg88;
    reg90=reg66+reg90; reg35=reg111*reg17; reg76=reg71+reg76; reg59=reg120-reg59; reg92=reg89+reg92;
    reg91=reg67*reg91; reg104=reg70+reg104; reg66=(*f.m).alpha_1*reg56; reg99=reg105+reg99; reg70=reg33*reg86;
    reg96=reg59-reg96; reg107=reg58*reg107; reg59=reg56*reg13; reg69=reg63*reg69; reg134=reg100-reg134;
    reg71=reg43*reg125; reg13=reg55*reg13; reg127=reg15+reg127; reg130=reg132+reg130; reg15=reg74*reg117;
    reg89=reg56*reg23; reg100=reg52*reg7; reg108=reg119-reg108; reg128=reg79+reg128; reg79=reg87*reg78;
    reg91=reg104+reg91; reg104=reg87*reg81; reg35=reg76+reg35; reg76=reg67*reg125; reg88=reg90+reg88;
    reg90=reg87*reg86; reg80=reg133+reg80; reg105=reg62*reg117; reg121=reg98-reg121; reg98=reg55*reg23;
    reg119=reg50*reg7; reg12=reg55*reg12; reg103=reg60+reg103; reg116=reg82*reg116; reg86=reg45*reg86;
    reg18=reg85+reg18; reg60=reg33*reg81; reg42=reg64+reg42; reg102=reg38+reg102; reg17=reg51*reg17;
    reg106=reg129+reg106; reg97=reg97*reg5; reg38=(*f.m).alpha_2*reg68; reg64=(*f.m).alpha_1*reg73; reg85=reg16*reg72;
    reg120=reg27*reg41; reg131=reg122+reg131; reg122=reg68*reg5; reg129=(*f.m).alpha_3*reg62; reg9=reg39+reg9;
    reg39=reg37*reg65; reg124=reg118+reg124; reg118=reg24*reg72; reg81=reg45*reg81; reg77=reg66+reg77;
    reg126=reg92+reg126; reg66=reg45*reg78; reg29=reg28+reg29; reg28=(*f.m).alpha_3*reg74; reg123=reg115+reg123;
    reg101=reg36+reg101; reg10=reg56*reg10; reg36=reg25*reg54; reg92=reg8*reg21; reg115=reg50*reg121;
    reg132=reg52*reg96; reg17=reg102+reg17; reg125=reg82*reg125; reg60=reg18+reg60; reg18=reg87*reg92;
    reg76=reg35+reg76; reg35=reg55*reg108; reg102=reg56*reg134; reg70=reg99+reg70; reg99=reg73*reg4;
    reg133=reg56*reg96; reg104=reg91+reg104; reg91=reg55*reg121; T reg135=reg49*reg108; T reg136=reg93*reg134;
    reg59=reg12-reg59; reg28=reg77+reg28; reg116=reg103+reg116; reg78=reg33*reg78; reg12=reg62*reg9;
    reg119=reg98+reg119; reg77=reg49*reg121; reg98=reg93*reg96; reg129=reg124+reg129; reg38=reg64+reg38;
    reg64=reg111*reg39; reg105=reg80+reg105; reg80=(*f.m).alpha_3*reg83; reg103=reg50*reg108; reg124=reg52*reg134;
    reg90=reg88+reg90; reg88=(*f.m).alpha_1*reg58; T reg137=(*f.m).alpha_2*reg63; reg16=reg16*reg27; T reg138=reg56*reg6;
    reg23=reg58*reg23; reg7=reg63*reg7; T reg139=reg74*reg9; reg100=reg89+reg100; reg4=reg55*reg4;
    reg89=reg110*reg39; reg15=reg130+reg15; reg106=reg97-reg106; reg131=reg122-reg131; reg86=reg42+reg86;
    reg65=reg29*reg65; reg42=reg101*reg14; reg13=reg10-reg13; reg6=reg73*reg6; reg10=reg45*reg92;
    reg71=reg127+reg71; reg66=reg123+reg66; reg81=reg126+reg81; reg97=reg27*reg54; reg85=reg120+reg85;
    reg69=reg107+reg69; reg79=reg128+reg79; reg118=reg36+reg118; reg36=reg26*reg72; reg11=reg11*reg25;
    reg117=reg95*reg117; reg107=reg104*reg70; reg120=reg81*reg70; reg78=reg116+reg78; reg6=reg13-reg6;
    reg13=reg90*reg60; reg116=reg86*reg60; reg92=reg33*reg92; reg125=reg17+reg125; reg17=reg104*reg129;
    reg122=reg79*reg28; reg9=reg95*reg9; reg7=reg23+reg7; reg80=reg38+reg80; reg117=reg69+reg117;
    reg23=(*f.m).alpha_2*reg16; reg38=(*f.m).alpha_1*reg11; reg69=reg81*reg129; reg123=reg66*reg28; reg137=reg88+reg137;
    reg39=reg51*reg39; reg95=(*f.m).alpha_3*reg95; reg88=reg68*reg131; reg115=reg132+reg115; reg89=reg15+reg89;
    reg15=reg43*reg42; reg5=reg73*reg5; reg139=reg100+reg139; reg110=reg110*reg65; reg18=reg76+reg18;
    reg76=reg73*reg106; reg35=reg102+reg35; reg27=reg26*reg27; reg25=reg24*reg25; reg24=reg73*reg131;
    reg91=reg133+reg91; reg41=reg72*reg41; reg63=reg63*reg2; reg36=reg97+reg36; reg26=reg51*reg2;
    reg99=reg59-reg99; reg59=reg85*reg21; reg97=reg112*reg131; reg4=reg138+reg4; reg135=reg136+reg135;
    reg111=reg111*reg65; reg12=reg119+reg12; reg77=reg98+reg77; reg14=reg118*reg14; reg10=reg71+reg10;
    reg71=reg67*reg42; reg64=reg105+reg64; reg98=reg68*reg106; reg103=reg124+reg103; reg100=reg112*reg106;
    reg102=reg60*reg129; reg16=reg16*reg1; reg105=reg86*reg104; reg69=reg123+reg69; reg10=reg10*reg80;
    reg119=reg78*reg28; reg4=reg5-reg4; reg18=reg18*reg80; reg17=reg122+reg17; reg100=reg135+reg100;
    reg5=reg51*reg26; reg122=reg51*reg63; reg97=reg77+reg97; reg77=reg48*reg26; reg98=reg103+reg98;
    reg103=reg48*reg63; reg88=reg115+reg88; reg115=reg50*reg99; reg95=reg137+reg95; reg123=reg52*reg6;
    reg124=reg53*reg26; reg23=reg38+reg23; reg76=reg35+reg76; reg35=reg53*reg63; reg41=(*f.m).alpha_3*reg41;
    reg38=(*f.m).alpha_1*reg25; reg126=(*f.m).alpha_2*reg27; reg24=reg91+reg24; reg13=reg107-reg13; reg91=reg55*reg99;
    reg116=reg120-reg116; reg107=reg81*reg90; reg120=reg56*reg6; reg127=reg37*reg1; reg65=reg51*reg65;
    reg9=reg7+reg9; reg42=reg82*reg42; reg39=reg117+reg39; reg92=reg125+reg92; reg7=1-var_inter[1];
    reg67=reg67*reg14; reg111=reg12+reg111; reg12=1-var_inter[0]; reg21=reg36*reg21; reg54=reg72*reg54;
    reg15=reg89+reg15; reg72=reg45*reg59; reg110=reg139+reg110; reg43=reg43*reg14; reg71=reg64+reg71;
    reg64=reg87*reg59; reg10=reg69+reg10; reg69=reg86*reg95; reg27=reg27*reg0; reg2=reg58*reg2;
    reg58=reg29*reg0; reg91=reg120+reg91; reg35=reg24+reg35; reg24=reg75*reg16; reg73=reg73*reg4;
    reg89=reg7*elem.pos(0)[1]; reg117=reg7*elem.pos(1)[1]; reg120=elem.pos(0)[0]*reg12; reg125=elem.pos(1)[0]*var_inter[0]; reg128=reg7*elem.pos(1)[0];
    reg130=reg7*elem.pos(0)[0]; reg41=reg23+reg41; reg126=reg38+reg126; reg54=(*f.m).alpha_3*reg54; reg23=reg66*reg13;
    reg38=reg79*reg116; reg132=reg37*reg127; reg5=reg100+reg5; reg100=reg12*elem.pos(0)[1]; reg133=elem.pos(1)[1]*var_inter[0];
    reg135=reg37*reg16; reg122=reg97+reg122; reg105=reg107-reg105; reg97=reg49*reg99; reg107=reg93*reg6;
    reg136=reg84*reg127; reg77=reg98+reg77; reg98=reg84*reg16; reg103=reg88+reg103; reg68=reg68*reg4;
    reg115=reg123+reg115; reg88=reg75*reg127; reg124=reg76+reg124; reg14=reg82*reg14; reg65=reg9+reg65;
    reg59=reg33*reg59; reg42=reg39+reg42; reg18=reg17+reg18; reg87=reg87*reg21; reg67=reg111+reg67;
    reg64=reg71+reg64; reg9=reg90*reg95; reg45=reg45*reg21; reg43=reg110+reg43; reg92=reg92*reg80;
    reg102=reg119+reg102; reg72=reg15+reg72; reg135=reg122+reg135; reg14=reg65+reg14; reg88=reg124+reg88;
    reg15=reg44*reg58; reg117=reg117-reg89; reg21=reg33*reg21; reg17=reg29*reg27; reg68=reg115+reg68;
    reg38=reg23-reg38; reg23=reg78*reg105; reg39=reg48*reg2; reg65=reg79*reg70; reg71=reg100+reg133;
    reg76=reg29*reg58; reg132=reg5+reg132; reg5=reg86*reg78; reg110=reg70*reg95; reg98=reg103+reg98;
    reg103=reg94*reg27; reg111=elem.pos(2)[1]*var_inter[0]; reg112=reg112*reg4; reg115=reg66*reg70; reg92=reg102+reg92;
    reg136=reg77+reg136; reg97=reg107+reg97; reg64=reg64*reg41; reg77=reg94*reg58; reg102=reg90*reg78;
    reg9=reg18+reg9; reg128=reg128-reg130; reg18=elem.pos(2)[1]*var_inter[1]; reg1=reg11*reg1; reg11=elem.pos(2)[0]*var_inter[1];
    reg72=reg72*reg41; reg45=reg43+reg45; reg73=reg91+reg73; reg43=reg53*reg2; reg91=reg120+reg125;
    reg69=reg10+reg69; reg87=reg67+reg87; reg54=reg126+reg54; reg59=reg42+reg59; reg24=reg35+reg24;
    reg10=reg44*reg27; reg35=elem.pos(2)[0]*var_inter[0]; reg35=reg35-reg91; reg72=reg69+reg72; reg0=reg25*reg0;
    reg77=reg136+reg77; reg45=reg45*reg54; reg10=reg24+reg10; reg24=elem.pos(3)[0]*reg12; reg112=reg97+reg112;
    reg25=reg51*reg2; reg111=reg111-reg71; reg17=reg135+reg17; reg64=reg9+reg64; reg9=reg86*reg79;
    reg87=reg87*reg54; reg128=reg11+reg128; reg23=reg38+reg23; reg11=elem.pos(3)[0]*var_inter[1]; reg38=elem.pos(3)[1]*reg12;
    reg21=reg14+reg21; reg39=reg68+reg39; reg84=reg84*reg1; reg76=reg132+reg76; reg59=reg59*reg41;
    reg14=reg81*reg78; reg117=reg18+reg117; reg5=reg115-reg5; reg18=reg66*reg60; reg110=reg92+reg110;
    reg75=reg75*reg1; reg43=reg73+reg43; reg42=reg66*reg90; reg103=reg98+reg103; reg67=elem.pos(3)[1]*var_inter[1];
    reg68=reg104*reg78; reg15=reg88+reg15; reg102=reg65-reg102; reg65=reg79*reg60; reg24=reg35+reg24;
    reg35=reg12*vectors[0][indices[0]+1]; reg21=reg21*reg54; reg59=reg110+reg59; reg69=reg12*vectors[0][indices[0]+0]; reg128=reg128-reg11;
    reg117=reg117-reg67; reg73=reg15*reg17; reg88=reg77*reg17; reg45=reg72+reg45; reg72=reg103*reg76;
    reg92=reg10*reg76; reg97=reg7*vectors[0][indices[0]+0]; reg98=reg7*vectors[0][indices[1]+0]; reg75=reg43+reg75; reg44=reg44*reg0;
    reg43=reg81*reg79; reg9=reg42-reg9; reg42=reg66*reg104; reg107=var_inter[0]*vectors[0][indices[1]+1]; reg110=var_inter[0]*vectors[0][indices[1]+0];
    reg87=reg64+reg87; reg38=reg111+reg38; reg13=reg13/reg23; reg37=reg37*reg1; reg25=reg112+reg25;
    reg102=reg102/reg23; reg64=reg7*vectors[0][indices[0]+1]; reg68=reg65-reg68; reg65=reg7*vectors[0][indices[1]+1]; reg116=reg116/reg23;
    reg5=reg5/reg23; reg94=reg94*reg0; reg84=reg39+reg84; reg14=reg18-reg14; reg5=reg5*reg87;
    reg18=var_inter[1]*vectors[0][indices[2]+1]; reg102=reg102*reg45; reg116=reg116*reg87; reg39=var_inter[1]*vectors[0][indices[2]+0]; reg13=reg13*reg45;
    reg21=reg59+reg21; reg73=reg92-reg73; reg59=reg10*reg77; reg92=reg15*reg103; reg68=reg68/reg23;
    reg14=reg14/reg23; reg105=reg105/reg23; reg9=reg9/reg23; reg43=reg42-reg43; reg29=reg29*reg0;
    reg37=reg25+reg37; reg25=reg128*reg38; reg42=var_inter[0]*vectors[0][indices[2]+0]; reg88=reg72-reg88; reg94=reg84+reg94;
    reg72=reg117*reg24; reg69=reg110+reg69; reg107=reg35+reg107; reg97=reg98-reg97; reg44=reg75+reg44;
    reg35=var_inter[0]*vectors[0][indices[2]+1]; reg64=reg65-reg64; reg105=reg105*reg21; reg116=reg13-reg116; reg9=reg9*reg21;
    reg102=reg5-reg102; reg45=reg68*reg45; reg29=reg37+reg29; reg69=reg42-reg69; reg5=reg94*reg73;
    reg13=reg88*reg44; reg72=reg25-reg72; reg25=var_inter[1]*vectors[0][indices[3]+0]; reg23=reg43/reg23; reg37=reg12*vectors[0][indices[3]+0];
    reg42=reg12*vectors[0][indices[3]+1]; reg107=reg35-reg107; reg18=reg64+reg18; reg92=reg59-reg92; reg35=var_inter[1]*vectors[0][indices[3]+1];
    reg87=reg14*reg87; reg97=reg39+reg97; reg14=reg29*reg92; reg69=reg37+reg69; reg37=1-(*f.m).resolution;
    reg21=reg23*reg21; reg9=reg102-reg9; reg25=reg97-reg25; reg116=reg105+reg116; reg87=reg45-reg87;
    reg5=reg13-reg5; reg38=reg38/reg72; reg128=reg128/reg72; reg35=reg18-reg35; reg24=reg24/reg72;
    reg13=reg10*reg29; reg117=reg117/reg72; reg18=reg94*reg17; reg23=reg44*reg17; reg39=reg103*reg29;
    reg42=reg107+reg42; reg43=reg117*reg42; reg14=reg5+reg14; reg116=reg37*reg116; reg5=(*f.m).resolution*reg129;
    reg87=reg21+reg87; reg21=(*f.m).resolution*reg28; reg45=reg77*reg29; reg59=reg44*reg103; reg13=reg23-reg13;
    reg23=reg10*reg94; reg64=reg15*reg29; reg65=reg94*reg76; reg9=reg37*reg9; reg39=reg18-reg39;
    reg18=reg44*reg76; reg68=reg128*reg69; reg75=reg24*reg25; reg84=reg35*reg38; reg25=reg38*reg25;
    reg45=reg65-reg45; reg39=reg39/reg14; reg64=reg18-reg64; reg13=reg13/reg14; reg35=reg35*reg24;
    reg42=reg128*reg42; reg18=reg44*reg77; reg9=reg5+reg9; reg75=reg68-reg75; reg69=reg117*reg69;
    reg43=reg84-reg43; reg5=(*f.m).resolution*reg95; reg23=reg59-reg23; reg87=reg37*reg87; reg59=reg15*reg94;
    reg116=reg21+reg116; reg88=reg88/reg14; reg35=reg42-reg35; reg21=(*f.m).resolution*reg39; reg116=(*f.m).deltaT*reg116;
    reg9=(*f.m).deltaT*reg9; reg42=(*f.m).resolution*reg13; reg43=reg75+reg43; reg87=reg5+reg87; reg60=reg60*reg37;
    reg78=reg78*reg37; reg23=reg23/reg14; reg59=reg18-reg59; reg64=reg64/reg14; reg73=reg73/reg14;
    reg45=reg45/reg14; reg69=reg25-reg69; reg5=(*f.m).resolution*reg73; reg42=reg60-reg42; reg43=0.5*reg43;
    reg18=(*f.m).resolution*reg23; reg78=reg21+reg78; reg81=reg81*reg37; reg35=reg35-reg9; reg79=reg79*reg37;
    reg66=reg66*reg37; reg21=(*f.m).resolution*reg64; reg70=reg70*reg37; reg87=(*f.m).deltaT*reg87; reg69=reg69-reg116;
    reg25=(*f.m).resolution*reg88; reg92=reg92/reg14; reg14=reg59/reg14; reg104=reg104*reg37; reg59=(*f.m).resolution*reg45;
    reg86=reg86*reg37; reg70=reg18+reg70; reg18=(*f.m).resolution*reg14; reg60=reg69*reg78; reg65=reg35*reg42;
    reg43=reg43-reg87; reg37=reg90*reg37; reg68=(*f.m).resolution*reg92; reg104=reg21+reg104; reg59=reg79-reg59;
    reg5=reg81-reg5; reg66=reg25+reg66; reg21=reg12*reg117; reg86=reg68+reg86; reg18=reg37-reg18;
    reg25=reg35*reg5; reg37=reg12*reg128; reg68=reg7*reg24; reg65=reg60+reg65; reg60=reg43*reg70;
    reg75=reg69*reg66; reg79=var_inter[0]*reg117; reg35=reg35*reg104; reg81=var_inter[0]*reg128; reg69=reg69*reg59;
    reg84=reg7*reg38; reg90=var_inter[1]*reg24; reg97=var_inter[1]*reg38; reg98=reg81-reg90; reg102=reg21+reg97;
    reg105=reg43*reg18; reg35=reg69+reg35; reg69=reg97-reg79; reg107=reg68-reg37; reg110=reg21-reg84;
    reg43=reg43*reg86; reg111=reg68+reg81; reg112=reg79+reg84; reg60=reg65+reg60; reg65=reg37+reg90;
    reg25=reg75+reg25; reg75=0.5*reg107; reg60=2*reg60; reg115=0.5*reg69; reg119=0.5*reg98;
    reg105=reg35+reg105; reg43=reg25+reg43; reg25=0.5*reg112; reg35=0.5*reg111; reg122=0.5*reg110;
    reg123=0.5*reg65; reg124=0.5*reg102; reg126=reg69*reg43; reg132=reg119*reg60; reg135=reg98*reg105;
    reg136=reg115*reg60; reg137=reg102*reg43; reg138=reg123*reg60; reg139=reg65*reg105; T reg140=reg124*reg60;
    T reg141=reg107*reg105; T reg142=reg122*reg60; T reg143=reg75*reg60; T reg144=reg110*reg43; T reg145=reg112*reg43;
    T reg146=reg35*reg60; T reg147=reg111*reg105; T reg148=reg25*reg60; reg132=reg126+reg132; reg142=reg141+reg142;
    reg145=reg145-reg146; reg136=reg135+reg136; reg144=reg143+reg144; reg139=reg139-reg140; reg138=reg138-reg137;
    reg148=reg148-reg147; reg142=reg72*reg142; reg144=reg72*reg144; reg139=reg72*reg139; reg138=reg72*reg138;
    reg148=reg72*reg148; reg136=reg72*reg136; reg145=reg72*reg145; reg132=reg72*reg132; sollicitation[indices[3]+0]+=ponderation*reg138;
    sollicitation[indices[0]+1]+=ponderation*reg142; sollicitation[indices[1]+1]+=ponderation*reg148; sollicitation[indices[3]+1]+=ponderation*reg139; sollicitation[indices[2]+1]+=ponderation*reg136; sollicitation[indices[0]+0]+=ponderation*reg144;
    sollicitation[indices[2]+0]+=ponderation*reg132; sollicitation[indices[1]+0]+=ponderation*reg145;
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

