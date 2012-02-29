
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
    T reg3=pow((*f.m).v1[1],2); T reg4=pow((*f.m).v1[0],2); reg2=1.0/reg2; T reg5=reg1*reg0; T reg6=pow((*f.m).v2[0],2);
    T reg7=reg2*reg5; T reg8=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg11=pow((*f.m).v2[1],2);
    reg3=reg4+reg3; reg4=pow((*f.m).v1[2],2); T reg12=reg10*reg7; T reg13=reg8*reg7; T reg14=reg9*reg7;
    reg4=reg3+reg4; reg3=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=pow((*f.m).v2[2],2); reg11=reg6+reg11;
    reg6=reg3*reg14; T reg17=reg8*reg13; reg16=reg11+reg16; reg11=reg15*reg14; T reg18=reg8*reg12;
    reg4=pow(reg4,0.5); T reg19=(*f.m).v1[0]/reg4; reg18=reg18+reg6; T reg20=1.0/(*f.m).elastic_modulus_1; T reg21=(*f.m).v1[1]/reg4;
    T reg22=reg15*reg12; reg16=pow(reg16,0.5); reg17=reg11-reg17; reg11=reg3*reg13; reg4=(*f.m).v1[2]/reg4;
    T reg23=reg20*reg17; T reg24=2*reg21; T reg25=2*reg19; T reg26=reg18*reg3; T reg27=reg11+reg22;
    T reg28=(*f.m).v2[0]/reg16; T reg29=(*f.m).v2[1]/reg16; reg16=(*f.m).v2[2]/reg16; T reg30=reg4*reg29; T reg31=reg21*reg16;
    reg26=reg23-reg26; reg23=2*reg4; T reg32=reg24*reg29; T reg33=reg9*reg5; T reg34=reg10*reg5;
    T reg35=reg27*reg10; T reg36=reg15*reg7; T reg37=reg10*reg13; reg14=reg20*reg14; T reg38=reg10*reg12;
    reg7=reg3*reg7; T reg39=reg2*reg0; reg5=reg8*reg5; T reg40=reg25*reg28; T reg41=pow(reg29,2);
    T reg42=pow(reg28,2); reg13=reg20*reg13; reg38=reg14-reg38; reg14=reg34*reg8; reg6=reg37+reg6;
    reg37=reg36*reg10; T reg43=reg5*reg8; T reg44=reg31-reg30; T reg45=reg33*reg3; T reg46=reg32*reg15;
    T reg47=reg2*reg1; T reg48=reg3*reg40; T reg49=reg39*reg10; T reg50=reg41*reg15; T reg51=reg42*reg3;
    reg33=reg33*reg15; T reg52=reg32*reg3; T reg53=reg20*reg40; T reg54=reg41*reg3; T reg55=reg42*reg20;
    T reg56=reg19*reg16; reg35=reg26-reg35; reg26=reg4*reg28; T reg57=reg39*reg8; reg12=reg3*reg12;
    T reg58=reg7*reg10; T reg59=pow(reg16,2); T reg60=reg23*reg16; T reg61=2*reg28; reg39=reg39*reg9;
    T reg62=reg60*reg10; reg54=reg55-reg54; reg52=reg53-reg52; reg53=reg59*reg10; reg55=reg59*reg8;
    reg36=reg36*reg20; reg51=reg50-reg51; reg50=reg61*reg29; reg17=reg17/reg35; T reg63=reg60*reg8;
    reg48=reg46-reg48; reg46=reg25*reg21; T reg64=2*reg44; T reg65=reg47*reg9; T reg66=reg49*reg8;
    T reg67=reg57*reg8; reg43=reg33-reg43; reg33=reg39*reg3; T reg68=reg32*reg8; T reg69=reg10*reg40;
    T reg70=reg41*reg8; reg14=reg45+reg14; reg45=reg42*reg10; reg39=reg39*reg15; reg38=reg38/reg35;
    reg18=reg18/reg35; reg58=reg13+reg58; reg37=reg11+reg37; T reg71=reg21*reg28; reg6=reg6/reg35;
    T reg72=pow(reg19,2); T reg73=pow(reg21,2); reg5=reg5*reg3; T reg74=reg47*reg8; reg34=reg34*reg15;
    reg47=reg47*reg10; T reg75=reg19*reg29; reg7=reg7*reg3; reg12=reg13+reg12; reg13=reg26-reg56;
    T reg76=reg42*reg38; T reg77=reg5+reg34; T reg78=reg72*reg17; T reg79=reg42*reg18; reg58=reg58/reg35;
    T reg80=reg72*reg6; reg67=reg39-reg67; reg12=reg12/reg35; reg39=reg21*reg29; T reg81=reg19*reg28;
    T reg82=reg50*reg18; reg66=reg33+reg66; reg33=reg46*reg17; T reg83=reg41*reg38; T reg84=reg64*reg13;
    T reg85=reg73*reg6; T reg86=reg41*reg18; T reg87=reg73*reg17; T reg88=pow(reg13,2); T reg89=reg47*reg8;
    T reg90=pow(reg44,2); T reg91=reg65*reg15; T reg92=reg75-reg71; reg65=reg65*reg3; T reg93=reg74*reg8;
    reg53=reg54-reg53; reg62=reg52-reg62; reg52=reg72*reg20; reg54=reg72*reg3; T reg94=reg73*reg15;
    reg49=reg49*reg15; reg57=reg57*reg3; reg55=reg51-reg55; reg51=pow(reg4,2); reg7=reg36-reg7;
    reg63=reg48-reg63; reg37=reg37/reg35; reg36=reg46*reg6; reg48=reg50*reg38; T reg95=reg59*reg9;
    reg70=reg45+reg70; reg27=reg27/reg35; reg14=reg14*reg3; reg45=reg60*reg9; reg68=reg69+reg68;
    reg69=reg73*reg3; reg43=reg43*reg20; T reg96=reg51*reg10; T reg97=reg28*reg29; T reg98=reg73*reg63;
    T reg99=reg41*reg63; T reg100=reg72*reg62; T reg101=reg73*reg55; reg7=reg7/reg35; reg69=reg52-reg69;
    reg52=reg72*reg53; T reg102=reg27*reg84; reg82=reg33+reg82; reg33=reg42*reg62; T reg103=reg41*reg55;
    T reg104=reg42*reg53; T reg105=reg59*reg18; T reg106=reg51*reg17; T reg107=reg27*reg88; reg86=reg87+reg86;
    reg74=reg74*reg3; reg89=reg65+reg89; reg93=reg91-reg93; reg65=reg51*reg8; reg54=reg94-reg54;
    reg87=reg57+reg49; reg47=reg47*reg15; reg66=reg66*reg3; reg14=reg43-reg14; reg67=reg67*reg20;
    reg77=reg77*reg10; reg43=reg63*reg39; reg91=reg62*reg81; reg94=reg55*reg39; T reg108=reg53*reg81;
    T reg109=reg72*reg37; T reg110=pow(reg92,2); reg68=reg45-reg68; reg79=reg78+reg79; reg45=reg27*reg90;
    reg70=reg95-reg70; reg78=reg73*reg8; reg95=reg72*reg10; reg75=reg71+reg75; reg71=reg4*reg16;
    T reg111=reg42*reg58; T reg112=reg41*reg58; T reg113=reg73*reg37; T reg114=reg51*reg6; T reg115=reg59*reg38;
    reg76=reg80+reg76; reg80=reg50*reg58; T reg116=reg46*reg37; T reg117=reg12*reg88; reg83=reg85+reg83;
    reg85=reg12*reg84; T reg118=reg12*reg90; reg48=reg36+reg48; reg36=reg51*reg9; reg77=reg14-reg77;
    reg26=reg56+reg26; reg66=reg67-reg66; reg118=reg76+reg118; reg87=reg87*reg10; reg14=2*reg29;
    reg56=reg61*reg16; reg93=reg93*reg20; reg111=reg109+reg111; reg89=reg89*reg3; reg67=reg7*reg90;
    reg43=reg91+reg43; reg76=reg71*reg68; reg91=reg7*reg84; reg103=reg104+reg103; reg104=reg59*reg70;
    reg109=reg71*reg70; reg94=reg108+reg94; reg99=reg33+reg99; reg33=reg59*reg68; reg85=reg48+reg85;
    reg96=reg69-reg96; reg48=reg12*reg110; reg65=reg54-reg65; reg54=reg2*reg75; reg69=reg97*reg2;
    reg108=reg25*reg4; reg115=reg114+reg115; reg117=reg83+reg117; reg78=reg95+reg78; reg83=reg51*reg70;
    reg101=reg52+reg101; reg102=reg82+reg102; reg45=reg79+reg45; reg52=reg59*reg58; reg79=reg51*reg68;
    reg82=reg51*reg37; reg98=reg100+reg98; reg95=reg27*reg110; reg105=reg106+reg105; reg80=reg116+reg80;
    reg100=reg74+reg47; reg106=reg19*reg21; reg114=reg7*reg88; reg107=reg86+reg107; reg112=reg113+reg112;
    reg86=reg28*reg16; reg113=reg42*reg96; reg104=reg103+reg104; reg103=reg50*reg69; reg116=reg73*reg45;
    T reg119=reg46*reg54; reg79=reg98+reg79; reg33=reg99+reg33; reg98=reg50*reg54; reg99=reg41*reg118;
    reg52=reg82+reg52; reg82=reg42*reg118; T reg120=reg41*reg117; T reg121=reg73*reg107; T reg122=reg72*reg45;
    T reg123=reg29*reg16; T reg124=reg44*reg13; T reg125=reg7*reg110; reg48=reg115+reg48; reg115=reg19*reg13;
    T reg126=reg97*reg117; T reg127=reg108*reg6; T reg128=reg38*reg56; T reg129=reg44*reg21; reg91=reg80+reg91;
    reg80=reg107*reg106; T reg130=reg97*reg85; T reg131=reg102*reg106; T reg132=reg75*reg54; reg76=reg43+reg76;
    reg43=reg75*reg69; reg109=reg94+reg109; reg94=reg41*reg65; T reg133=reg73*reg65; T reg134=reg72*reg96;
    reg64=reg64*reg92; reg87=reg66-reg87; reg66=reg1*reg26; T reg135=reg86*reg1; T reg136=2*reg13;
    T reg137=reg46*reg69; reg77=reg77/reg35; reg83=reg101+reg83; reg101=reg24*reg4; reg95=reg105+reg95;
    reg105=reg41*reg85; T reg138=reg18*reg56; T reg139=reg108*reg17; reg31=reg30+reg31; reg30=reg73*reg102;
    reg78=reg36-reg78; reg36=reg72*reg102; reg67=reg111+reg67; reg111=reg72*reg107; T reg140=reg42*reg85;
    reg89=reg93-reg89; reg93=reg42*reg117; reg114=reg112+reg114; reg100=reg100*reg10; reg112=reg14*reg16;
    reg17=reg101*reg17; T reg141=reg44*reg19; reg18=reg18*reg112; T reg142=reg27*reg64; T reg143=reg108*reg37;
    T reg144=reg58*reg56; T reg145=reg12*reg64; T reg146=reg21*reg13; reg128=reg127+reg128; reg6=reg101*reg6;
    reg87=reg87/reg35; reg138=reg139+reg138; reg100=reg89-reg100; reg89=reg28*reg13; reg136=reg136*reg92;
    reg127=reg19*reg4; reg126=reg80+reg126; reg80=reg124*reg114; reg125=reg52+reg125; reg115=reg129+reg115;
    reg52=reg44*reg29; reg98=reg33+reg98; reg33=reg56*reg66; reg38=reg38*reg112; reg129=reg65*reg39;
    reg139=reg96*reg81; T reg147=reg108*reg135; T reg148=reg73*reg95; T reg149=reg88*reg114; reg120=reg121+reg120;
    reg133=reg134+reg133; reg121=reg88*reg67; reg134=reg51*reg78; reg99=reg116+reg99; reg140=reg36+reg140;
    reg137=reg83+reg137; reg36=reg90*reg91; reg83=reg124*reg91; reg130=reg131+reg130; reg116=reg26*reg66;
    reg132=reg76+reg132; reg76=reg97*reg118; reg131=reg106*reg45; T reg150=reg59*reg78; reg94=reg113+reg94;
    reg113=reg26*reg135; reg43=reg109+reg43; reg109=reg106*reg2; reg103=reg104+reg103; reg104=reg56*reg135;
    T reg151=reg108*reg66; reg119=reg79+reg119; reg79=reg123*reg0; T reg152=reg0*reg31; T reg153=reg72*reg95;
    T reg154=reg42*reg48; reg82=reg122+reg82; reg122=reg90*reg67; T reg155=reg90*reg114; reg93=reg111+reg93;
    reg111=reg77*reg81; T reg156=reg77*reg39; T reg157=reg41*reg48; T reg158=reg88*reg91; reg105=reg30+reg105;
    reg30=reg77*reg75; T reg159=elem.pos(2)[0]-elem.pos(0)[0]; reg121=reg99+reg121; reg99=reg30*reg75; reg134=reg133+reg134;
    reg122=reg82+reg122; reg35=reg100/reg35; reg12=reg12*reg136; reg38=reg6+reg38; reg145=reg128+reg145;
    reg6=reg30*reg40; reg36=reg140+reg36; reg82=reg111*reg40; reg100=reg46*reg109; reg128=reg4*reg92;
    reg83=reg130+reg83; reg130=reg156*reg75; reg80=reg126+reg80; reg126=reg97*reg48; reg24=reg24*reg13;
    reg133=reg95*reg106; reg140=reg21*reg4; T reg160=reg50*reg109; reg150=reg94+reg150; reg94=reg87*reg146;
    T reg161=reg101*reg152; reg151=reg119+reg151; reg119=reg87*reg141; T reg162=reg71*reg78; reg129=reg139+reg129;
    reg139=reg87*reg115; T reg163=reg77*reg71; T reg164=reg101*reg79; reg147=reg137+reg147; reg137=elem.pos(1)[1]-elem.pos(0)[1];
    T reg165=elem.pos(1)[0]-elem.pos(0)[0]; reg76=reg131+reg76; reg131=reg124*reg67; reg113=reg43+reg113; reg43=reg31*reg79;
    reg58=reg58*reg112; reg37=reg101*reg37; T reg166=reg7*reg64; reg144=reg143+reg144; reg116=reg132+reg116;
    reg132=reg31*reg152; reg143=elem.pos(2)[1]-elem.pos(0)[1]; T reg167=reg156*reg32; reg149=reg120+reg149; reg120=reg111*reg32;
    T reg168=reg29*reg13; T reg169=reg112*reg79; reg104=reg103+reg104; reg142=reg138+reg142; reg148=reg157+reg148;
    reg103=reg88*reg125; reg138=reg44*reg28; reg18=reg17+reg18; reg27=reg27*reg136; reg89=reg52+reg89;
    reg17=reg90*reg125; reg155=reg93+reg155; reg52=reg156*reg40; reg93=reg127*reg1; reg157=reg30*reg32;
    reg25=reg44*reg25; T reg170=reg112*reg152; reg158=reg105+reg158; reg33=reg98+reg33; reg154=reg153+reg154;
    reg17=reg154+reg17; reg167=reg149+reg167; reg98=reg94*reg24; reg105=reg35*reg138; reg120=reg121+reg120;
    reg121=reg119*reg24; reg149=reg159*reg137; reg153=reg35*reg168; reg132=reg116+reg132; reg170=reg33+reg170;
    reg12=reg38+reg12; reg33=reg35*reg89; reg38=reg56*reg93; reg160=reg150+reg160; reg61=reg44*reg61;
    reg116=reg75*reg109; reg161=reg151+reg161; reg150=(*f.m).alpha_1*reg73; reg151=(*f.m).alpha_2*reg41; reg162=reg129+reg162;
    reg129=reg87*reg128; reg164=reg147+reg164; reg147=(*f.m).alpha_2*reg42; reg154=reg143*reg165; reg27=reg18+reg27;
    reg18=reg16*reg92; reg131=reg76+reg131; reg43=reg113+reg43; reg7=reg7*reg136; reg58=reg37+reg58;
    reg166=reg144+reg166; reg37=reg111*reg75; reg169=reg104+reg169; reg76=reg163*reg40; reg104=reg108*reg93;
    reg100=reg134+reg100; reg113=(*f.m).alpha_1*reg72; reg14=reg14*reg13; reg130=reg80+reg130; reg80=reg94*reg115;
    reg103=reg148+reg103; reg134=reg163*reg32; reg126=reg133+reg126; reg133=reg124*reg125; reg144=reg41*reg145;
    reg148=reg73*reg142; reg99=reg83+reg99; reg83=reg139*reg115; T reg171=reg139*reg24; T reg172=reg44*reg4;
    reg19=reg19*reg92; reg157=reg158+reg157; reg158=reg94*reg25; T reg173=reg42*reg145; reg52=reg155+reg52;
    reg155=reg72*reg142; T reg174=reg139*reg25; reg6=reg36+reg6; reg36=reg140*reg0; T reg175=reg119*reg25;
    reg82=reg122+reg82; reg122=reg166*reg88; reg144=reg148+reg144; reg113=reg147+reg113; reg147=reg129*reg25;
    reg148=reg90*(*f.m).alpha_3; reg151=reg150+reg151; reg150=reg169*reg132; reg175=reg82+reg175; reg82=reg33*reg14;
    T reg176=reg170*reg43; reg171=reg157+reg171; reg157=reg105*reg61; reg158=reg52+reg158; reg116=reg162+reg116;
    reg52=reg42*reg12; reg38=reg160+reg38; reg160=reg72*reg27; reg173=reg155+reg173; reg155=reg166*reg90;
    reg162=reg101*reg36; T reg177=reg119*reg115; reg37=reg131+reg37; reg104=reg100+reg104; reg100=reg35*reg18;
    reg131=(*f.m).alpha_2*reg59; reg76=reg17+reg76; reg17=reg41*reg12; T reg178=reg105*reg14; T reg179=reg73*reg27;
    reg121=reg120+reg121; reg120=reg153*reg14; reg174=reg6+reg174; reg6=reg33*reg61; T reg180=reg164*reg132;
    reg134=reg103+reg134; reg103=reg129*reg24; T reg181=reg153*reg61; T reg182=(*f.m).alpha_1*reg51; T reg183=reg88*(*f.m).alpha_3;
    reg98=reg167+reg98; reg80=reg130+reg80; reg149=reg154-reg149; reg130=reg153*reg89; reg154=reg97*reg145;
    reg7=reg58+reg7; reg21=reg21*reg92; reg58=reg106*reg142; reg133=reg126+reg133; reg126=reg112*reg36;
    reg167=reg44*reg16; reg28=reg28*reg92; T reg184=reg161*reg43; T reg185=reg33*reg89; T reg186=reg77*reg26;
    reg4=reg4*reg13; reg19=reg172+reg19; reg83=reg99+reg83; reg99=reg163*reg75; reg172=reg26*reg93;
    reg130=reg80+reg130; reg122=reg144+reg122; reg99=reg133+reg99; reg80=reg186*reg32; reg159=reg159/reg149;
    reg143=reg143/reg149; reg133=reg169*reg161; reg144=reg170*reg164; reg147=reg76+reg147; reg76=reg100*reg61;
    reg184=reg180-reg184; reg6=reg174+reg6; reg155=reg173+reg155; reg173=reg186*reg40; reg52=reg160+reg52;
    reg160=reg7*reg90; reg137=reg137/reg149; reg165=reg165/reg149; reg29=reg29*reg92; reg148=reg113+reg148;
    reg113=reg16*reg13; reg28=reg167+reg28; reg183=reg151+reg183; reg38=reg126+reg38; reg21=reg4+reg21;
    reg181=reg158+reg181; reg4=reg97*reg12; reg126=reg106*reg27; reg157=reg175+reg157; reg151=reg166*reg124;
    reg154=reg58+reg154; reg103=reg134+reg103; reg58=reg100*reg14; reg176=reg150-reg176; reg185=reg83+reg185;
    reg82=reg171+reg82; reg83=reg129*reg115; reg134=reg87*reg19; reg150=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; reg17=reg179+reg17;
    reg158=reg7*reg88; reg120=reg98+reg120; reg98=reg105*reg89; reg177=reg37+reg177; reg106=(*f.m).alpha_1*reg106;
    reg37=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; reg167=reg110*(*f.m).alpha_3; reg162=reg104+reg162; reg104=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg171=vectors[0][indices[1]+1]-vectors[0][indices[0]+1];
    reg131=reg182+reg131; reg178=reg121+reg178; reg77=reg77*reg31; reg97=(*f.m).alpha_2*reg97; reg121=reg31*reg36;
    reg172=reg116+reg172; reg44=reg44*reg92; reg80=reg122+reg80; reg116=reg134*reg24; reg167=reg131+reg167;
    reg83=reg99+reg83; reg99=reg162*reg176; reg122=reg100*reg89; reg97=reg106+reg97; reg133=reg144-reg133;
    reg87=reg87*reg21; reg158=reg17+reg158; reg32=reg77*reg32; reg58=reg103+reg58; reg98=reg177+reg98;
    reg151=reg154+reg151; reg17=reg186*reg75; reg4=reg126+reg4; reg103=reg7*reg124; reg106=reg6*reg130;
    reg126=reg82*reg130; reg131=reg184*reg38; reg144=reg181*reg185; reg154=reg120*reg185; reg121=reg172+reg121;
    reg29=reg113+reg29; reg113=reg150*reg165; reg172=reg148*reg157; reg174=reg37*reg159; reg173=reg155+reg173;
    reg155=reg183*reg181; reg175=reg134*reg25; reg177=(*f.m).alpha_2*reg86; reg179=reg104*reg137; reg180=reg148*reg178;
    reg76=reg147+reg76; reg147=reg171*reg143; reg182=(*f.m).alpha_1*reg127; T reg187=reg124*(*f.m).alpha_3; T reg188=reg183*reg120;
    T reg189=reg35*reg28; reg40=reg77*reg40; reg160=reg52+reg160; reg187=reg97+reg187; reg52=reg87*reg25;
    reg40=reg160+reg40; reg17=reg151+reg17; reg97=reg134*reg115; reg175=reg173+reg175; reg151=reg189*reg61;
    reg103=reg4+reg103; reg4=reg77*reg75; reg160=reg170*reg121; reg173=reg6*reg120; reg32=reg158+reg32;
    reg158=reg87*reg24; T reg190=reg161*reg121; reg174=reg113-reg174; reg126=reg154-reg126; reg113=reg162*reg132;
    reg35=reg35*reg29; reg154=reg181*reg82; reg106=reg144-reg106; reg104=reg104*reg165; reg13=reg92*reg13;
    reg171=reg171*reg159; reg188=reg180+reg188; reg58=reg167*reg58; reg116=reg80+reg116; reg80=reg189*reg14;
    reg131=reg99-reg131; reg122=reg83+reg122; reg83=reg148*reg98; reg99=reg183*reg130; reg144=(*f.m).alpha_2*reg123;
    reg180=(*f.m).alpha_1*reg140; T reg191=reg38*reg132; reg150=reg150*reg137; reg76=reg167*reg76; T reg192=reg133*reg121;
    reg155=reg172+reg155; reg179=reg147-reg179; reg147=reg44*(*f.m).alpha_3; reg177=reg182+reg177; reg37=reg37*reg143;
    reg52=reg40+reg52; reg40=reg35*reg61; reg174=reg179+reg174; reg160=reg191-reg160; reg171=reg104-reg171;
    elem.epsilon[0][1]=reg171; reg104=reg169*reg121; reg97=reg17+reg97; reg17=reg189*reg89; reg147=reg177+reg147;
    reg151=reg175+reg151; reg172=(*f.m).deltaT*reg183; reg4=reg103+reg4; reg103=reg87*reg115; reg175=(*f.m).deltaT*reg148;
    reg177=reg38*reg43; reg179=reg38*reg161; reg144=reg180+reg144; reg180=reg170*reg162; reg182=reg164*reg121;
    reg190=reg113-reg190; reg150=reg37-reg150; elem.epsilon[0][0]=reg150; reg37=reg13*(*f.m).alpha_3; reg192=reg131+reg192;
    reg113=reg162*reg43; reg76=reg155+reg76; reg131=reg187*reg6; reg80=reg116+reg80; reg99=reg83+reg99;
    reg83=reg157*reg126; reg122=reg167*reg122; reg116=reg96*reg90; reg173=reg154-reg173; reg154=reg65*reg88;
    reg155=reg55*reg88; reg191=reg35*reg14; reg158=reg32+reg158; reg32=reg53*reg90; T reg193=reg106*reg178;
    reg58=reg188+reg58; reg188=reg187*reg82; reg80=reg147*reg80; reg103=reg4+reg103; reg4=reg169*reg162;
    T reg194=reg78*reg110; T reg195=(*f.m).deltaT*reg187; T reg196=reg35*reg89; T reg197=reg98*reg6; reg179=reg180-reg179;
    reg180=reg70*reg110; reg155=reg32+reg155; reg32=reg187*reg185; reg193=reg83-reg193; reg174=0.5*reg174;
    elem.epsilon[0][2]=reg174; reg83=reg157*reg185; reg122=reg99+reg122; reg154=reg116+reg154; reg99=reg63*reg88;
    reg176=reg176/reg192; reg17=reg97+reg17; reg188=reg58+reg188; reg58=reg62*reg90; reg104=reg177-reg104;
    reg184=reg184/reg192; reg190=reg190/reg192; reg97=reg173*reg98; reg182=reg113-reg182; reg191=reg158+reg191;
    reg37=reg144+reg37; reg113=reg178*reg185; reg40=reg52+reg40; reg52=reg150-reg175; reg116=reg82*reg98;
    reg144=reg171-reg172; reg151=reg147*reg151; reg160=reg160/reg192; reg131=reg76+reg131; reg76=reg38*reg164;
    reg158=reg68*reg110; reg99=reg58+reg99; reg97=reg193+reg97; reg17=reg147*reg17; reg58=reg98*reg120;
    reg177=reg84*reg109; reg194=reg154+reg194; reg154=reg178*reg130; reg116=reg113-reg116; reg191=reg37*reg191;
    reg80=reg188+reg80; reg104=reg104/reg192; reg113=reg84*reg69; reg182=reg182/reg192; reg133=reg133/reg192;
    reg179=reg179/reg192; reg40=reg37*reg40; reg188=reg174-reg195; reg180=reg155+reg180; reg155=reg176*reg52;
    reg151=reg131+reg151; reg131=reg184*reg144; reg193=reg190*reg144; T reg198=reg160*reg52; reg76=reg4-reg76;
    reg4=reg6*reg178; T reg199=reg157*reg82; T reg200=reg181*reg98; reg197=reg83-reg197; reg83=reg157*reg130;
    reg196=reg103+reg196; reg32=reg122+reg32; reg192=reg76/reg192; reg76=reg84*reg54; reg113=reg180+reg113;
    reg103=reg64*reg135; reg191=reg80+reg191; reg52=reg104*reg52; reg40=reg151+reg40; reg80=reg133*reg188;
    reg131=reg155-reg131; reg198=reg193-reg198; reg122=reg179*reg188; reg151=reg181*reg178; reg4=reg199-reg4;
    reg155=reg157*reg120; reg144=reg182*reg144; reg200=reg83-reg200; reg197=reg197/reg97; reg106=reg106/reg97;
    reg58=reg154-reg58; reg17=reg32+reg17; reg177=reg194+reg177; reg116=reg116/reg97; reg126=reg126/reg97;
    reg32=reg64*reg93; reg196=reg37*reg196; reg158=reg99+reg158; reg103=reg113+reg103; reg83=reg136*reg79;
    reg196=reg17+reg196; reg188=reg188*reg192; reg131=reg80+reg131; reg151=reg155-reg151; reg122=reg198-reg122;
    reg4=reg4/reg97; reg173=reg173/reg97; reg200=reg200/reg97; reg144=reg52-reg144; reg58=reg58/reg97;
    reg32=reg177+reg32; reg17=reg136*reg36; reg76=reg158+reg76; reg52=reg64*reg66; reg197=reg197*reg191;
    reg116=reg116*reg40; reg106=reg106*reg191; reg126=reg126*reg40; reg144=reg188+reg144; reg52=reg76+reg52;
    reg97=reg151/reg97; reg76=reg136*reg152; reg80=reg131*reg38; reg83=reg103+reg83; reg173=reg173*reg196;
    reg106=reg126-reg106; reg4=reg4*reg196; reg116=reg197-reg116; reg40=reg58*reg40; reg191=reg200*reg191;
    reg17=reg32+reg17; reg32=reg96*reg141; reg58=reg65*reg146; reg99=reg53*reg141; reg103=reg55*reg146;
    reg113=reg122*reg164; reg126=reg131*reg162; reg151=reg122*reg169; reg55=reg55*reg168; reg4=reg116-reg4;
    reg116=reg131*reg17; reg196=reg97*reg196; reg191=reg40-reg191; reg40=1-(*f.m).resolution; reg97=reg122*reg83;
    reg154=PNODE(1).dep[0]-PNODE(0).dep[0]; reg155=reg144*reg161; reg158=PNODE(2).dep[1]-PNODE(0).dep[1]; reg177=PNODE(2).dep[0]-PNODE(0).dep[0]; reg96=reg96*reg138;
    reg180=reg70*reg128; reg103=reg99+reg103; reg58=reg32+reg58; reg32=reg78*reg128; reg113=reg126+reg113;
    reg99=reg62*reg141; reg65=reg65*reg168; reg151=reg80+reg151; reg80=PNODE(1).dep[1]-PNODE(0).dep[1]; reg126=reg63*reg146;
    reg53=reg53*reg138; reg106=reg173+reg106; reg173=reg144*reg170; reg76=reg52+reg76; reg52=reg109*reg115;
    reg78=reg78*reg18; reg65=reg96+reg65; reg148=reg148*(*f.m).resolution; reg180=reg103+reg180; reg96=reg69*reg115;
    reg70=reg70*reg18; reg103=reg177*reg165; reg188=reg158*reg137; reg173=reg151+reg173; reg62=reg62*reg138;
    reg63=reg63*reg168; reg55=reg53+reg55; reg155=reg113+reg155; reg183=reg183*(*f.m).resolution; reg106=reg106*reg40;
    reg4=reg4*reg40; reg126=reg99+reg126; reg53=reg80*reg143; reg99=reg68*reg128; reg191=reg196+reg191;
    reg113=reg154*reg159; reg32=reg58+reg32; reg97=reg116+reg97; reg58=reg144*reg76; reg155=reg175+reg155;
    reg109=reg109*reg89; reg78=reg65+reg78; reg52=reg32+reg52; reg188=reg53-reg188; reg113=reg103-reg113;
    reg187=reg187*(*f.m).resolution; reg32=reg93*reg19; reg58=reg97+reg58; reg158=reg158*reg165; reg191=reg191*reg40;
    reg106=reg148+reg106; reg4=reg183+reg4; reg167=(*f.m).deltaT*reg167; reg63=reg62+reg63; reg68=reg68*reg18;
    reg99=reg126+reg99; reg53=reg54*reg115; reg154=reg154*reg143; reg177=reg177*reg137; reg80=reg80*reg159;
    reg70=reg55+reg70; reg69=reg69*reg89; reg173=reg172+reg173; reg55=reg135*reg19; reg96=reg180+reg96;
    reg58=reg167+reg58; reg62=reg155+reg173; reg54=reg54*reg89; reg177=reg154-reg177; reg65=reg160*(*f.m).resolution;
    reg69=reg70+reg69; reg55=reg96+reg55; reg70=reg190*(*f.m).resolution; reg96=reg36*reg21; reg32=reg52+reg32;
    reg52=reg104*(*f.m).resolution; reg97=reg182*(*f.m).resolution; reg135=reg135*reg28; reg120=reg40*reg120; reg191=reg187+reg191;
    reg106=(*f.m).deltaT*reg106; reg4=(*f.m).deltaT*reg4; reg178=reg40*reg178; reg157=reg157*reg40; reg181=reg181*reg40;
    reg113=reg188+reg113; reg109=reg78+reg109; reg78=reg176*(*f.m).resolution; reg103=reg184*(*f.m).resolution; reg80=reg158-reg80;
    reg116=reg79*reg21; reg53=reg99+reg53; reg99=reg66*reg19; reg93=reg93*reg28; reg98=reg98*reg40;
    reg130=reg40*reg130; reg68=reg63+reg68; reg96=reg32+reg96; reg32=reg122*reg43; reg63=reg80-reg4;
    reg113=0.5*reg113; reg135=reg69+reg135; reg79=reg79*reg29; reg69=reg177-reg106; reg93=reg109+reg93;
    reg36=reg36*reg29; reg66=reg66*reg28; reg54=reg68+reg54; reg185=reg40*reg185; reg82=reg82*reg40;
    reg157=reg78+reg157; reg103=reg181-reg103; reg65=reg178-reg65; reg120=reg70+reg120; reg68=reg152*reg21;
    reg62=reg58+reg62; reg99=reg53+reg99; reg6=reg40*reg6; reg191=(*f.m).deltaT*reg191; reg40=reg133*(*f.m).resolution;
    reg55=reg116+reg55; reg97=reg130-reg97; reg52=reg98+reg52; reg53=reg192*(*f.m).resolution; reg70=reg179*(*f.m).resolution;
    reg78=reg131*reg121; reg98=reg131*reg96; reg109=reg103*reg63; reg116=reg144*reg132; reg126=reg122*reg55;
    reg68=reg99+reg68; reg99=reg157*reg69; reg32=reg78+reg32; reg78=reg65*reg69; reg130=reg113-reg191;
    reg148=reg120*reg63; reg151=reg52*reg69; reg154=reg97*reg63; reg53=reg185+reg53; reg79=reg135+reg79;
    reg6=reg40+reg6; reg66=reg54+reg66; reg62=reg62/3; reg152=reg152*reg29; reg70=reg82-reg70;
    reg36=reg93+reg36; reg122=reg122*reg79; reg131=reg131*reg36; reg173=reg173-reg62; reg155=reg155-reg62;
    reg40=reg144*reg68; reg126=reg98+reg126; reg152=reg66+reg152; reg116=reg32+reg116; reg109=reg99+reg109;
    reg32=reg6*reg130; reg54=reg53*reg130; reg154=reg151+reg154; reg148=reg78+reg148; reg66=reg70*reg130;
    reg54=reg154+reg54; reg147=(*f.m).deltaT*reg147; reg144=reg144*reg152; reg122=reg131+reg122; reg66=reg148+reg66;
    reg32=reg109+reg32; reg155=pow(reg155,2); reg40=reg126+reg40; reg173=pow(reg173,2); reg62=reg58-reg62;
    reg116=reg195+reg116; reg54=2*reg54; reg63=reg63*reg66; reg173=reg155+reg173; reg62=pow(reg62,2);
    reg4=reg171-reg4; reg58=2*reg116; reg106=reg150-reg106; reg40=reg147+reg40; reg37=(*f.m).deltaT*reg37;
    reg69=reg69*reg32; reg144=reg122+reg144; reg157=reg157*reg106; reg103=reg103*reg4; reg69=reg63+reg69;
    reg130=reg130*reg54; reg191=reg174-reg191; reg65=reg65*reg106; reg120=reg120*reg4; reg62=reg173+reg62;
    reg58=reg116*reg58; reg63=2*reg40; reg144=reg37+reg144; reg58=reg62+reg58; reg130=reg69+reg130;
    reg63=reg40*reg63; reg40=2*reg144; reg103=reg157+reg103; reg4=reg97*reg4; reg106=reg52*reg106;
    reg6=reg6*reg191; reg120=reg65+reg120; reg70=reg70*reg191; reg63=reg58+reg63; reg191=reg53*reg191;
    reg4=reg106+reg4; reg70=reg120+reg70; elem.sigma[0][1]=reg70; reg40=reg144*reg40; reg6=reg103+reg6;
    elem.sigma[0][0]=reg6; reg130=reg130*reg149; reg52=reg81*reg6; reg53=reg39*reg70; reg58=0.33333333333333331483*reg130;
    reg130=0.16666666666666665741*reg130; reg191=reg4+reg191; elem.sigma[0][2]=reg191; reg4=reg72*reg6; reg62=reg73*reg70;
    reg40=reg63+reg40; reg6=reg42*reg6; reg70=reg41*reg70; reg40=1.5*reg40; reg63=reg50*reg191;
    reg70=reg6+reg70; reg6=reg46*reg191; reg62=reg4+reg62; reg130=reg58+reg130; reg191=reg75*reg191;
    reg52=reg53+reg52; elem.sigma_von_mises=pow(reg40,0.5); elem.sigma_local[0][2]=reg52+reg191; elem.sigma_local[0][0]=reg62+reg6; elem.sigma_local[0][1]=reg70+reg63;
    elem.ener=reg130/2;
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
    T reg3=reg1*reg0; reg2=1.0/reg2; T reg4=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg5=1.0/(*f.m).elastic_modulus_3; T reg6=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg7=reg2*reg3; T reg8=pow((*f.m).v1[0],2); T reg9=pow((*f.m).v2[1],2); T reg10=pow((*f.m).v2[0],2); T reg11=pow((*f.m).v1[1],2);
    T reg12=pow((*f.m).v2[2],2); reg9=reg10+reg9; reg10=1.0/(*f.m).elastic_modulus_2; T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg14=pow((*f.m).v1[2],2);
    reg11=reg8+reg11; reg8=reg5*reg7; T reg15=reg6*reg7; T reg16=reg4*reg7; reg14=reg11+reg14;
    reg11=reg10*reg8; T reg17=reg6*reg16; T reg18=reg6*reg15; T reg19=reg13*reg8; reg12=reg9+reg12;
    reg9=reg10*reg16; T reg20=reg13*reg15; reg17=reg17+reg19; reg14=pow(reg14,0.5); reg12=pow(reg12,0.5);
    T reg21=1.0/(*f.m).elastic_modulus_1; reg18=reg11-reg18; reg11=reg20+reg9; T reg22=reg17*reg13; T reg23=(*f.m).v1[1]/reg14;
    T reg24=(*f.m).v1[2]/reg14; T reg25=(*f.m).v2[1]/reg12; T reg26=reg21*reg18; T reg27=(*f.m).v2[2]/reg12; T reg28=reg10*reg7;
    T reg29=reg23*reg27; T reg30=reg11*reg4; T reg31=reg24*reg25; T reg32=reg4*reg16; reg7=reg13*reg7;
    reg8=reg21*reg8; T reg33=reg4*reg15; reg14=(*f.m).v1[0]/reg14; T reg34=reg4*reg3; T reg35=reg5*reg3;
    reg3=reg6*reg3; reg12=(*f.m).v2[0]/reg12; reg22=reg26-reg22; reg26=reg2*reg0; T reg36=reg35*reg10;
    reg35=reg35*reg13; T reg37=reg3*reg6; T reg38=reg26*reg4; reg32=reg8-reg32; reg8=reg2*reg1;
    reg15=reg21*reg15; T reg39=reg7*reg4; T reg40=reg26*reg6; reg16=reg13*reg16; T reg41=reg34*reg6;
    reg19=reg33+reg19; reg33=2*reg12; reg26=reg26*reg5; T reg42=reg28*reg4; T reg43=reg14*reg27;
    reg30=reg22-reg30; reg22=2*reg14; T reg44=reg29-reg31; T reg45=reg24*reg12; T reg46=reg14*reg25;
    T reg47=reg23*reg12; reg39=reg15+reg39; T reg48=pow(reg14,2); T reg49=pow(reg23,2); T reg50=reg45-reg43;
    reg32=reg32/reg30; reg19=reg19/reg30; reg42=reg20+reg42; reg17=reg17/reg30; T reg51=reg40*reg6;
    T reg52=pow(reg12,2); T reg53=reg26*reg13; reg26=reg26*reg10; reg18=reg18/reg30; reg41=reg35+reg41;
    reg35=reg38*reg6; reg37=reg36-reg37; reg36=2*reg44; T reg54=pow(reg25,2); T reg55=reg8*reg4;
    T reg56=reg8*reg5; reg34=reg34*reg10; reg8=reg8*reg6; reg3=reg3*reg13; T reg57=reg33*reg25;
    reg28=reg28*reg21; T reg58=reg22*reg23; reg7=reg7*reg13; reg16=reg15+reg16; reg15=pow(reg27,2);
    T reg59=reg57*reg17; T reg60=reg58*reg18; T reg61=reg46-reg47; T reg62=reg36*reg50; T reg63=pow(reg50,2);
    T reg64=pow(reg44,2); T reg65=reg48*reg19; T reg66=reg52*reg32; T reg67=reg49*reg19; T reg68=reg54*reg32;
    T reg69=reg58*reg19; T reg70=reg57*reg32; reg16=reg16/reg30; reg39=reg39/reg30; reg7=reg28-reg7;
    reg28=reg52*reg17; T reg71=reg48*reg18; reg40=reg40*reg13; reg38=reg38*reg10; reg11=reg11/reg30;
    reg37=reg37*reg21; reg41=reg41*reg13; T reg72=reg3+reg34; reg42=reg42/reg30; reg51=reg26-reg51;
    reg35=reg53+reg35; reg26=reg56*reg10; reg53=pow(reg24,2); reg56=reg56*reg13; T reg73=reg8*reg6;
    T reg74=reg55*reg6; T reg75=reg49*reg18; T reg76=reg54*reg17; T reg77=reg53*reg19; reg74=reg56+reg74;
    reg51=reg51*reg21; reg56=reg15*reg32; reg73=reg26-reg73; reg55=reg55*reg10; reg26=reg52*reg39;
    T reg78=reg48*reg42; reg59=reg60+reg59; reg60=reg58*reg42; T reg79=pow(reg61,2); T reg80=reg40+reg38;
    T reg81=reg11*reg62; reg8=reg8*reg13; reg35=reg35*reg13; T reg82=reg16*reg62; reg70=reg69+reg70;
    reg72=reg72*reg4; reg69=reg15*reg17; T reg83=reg53*reg18; T reg84=reg57*reg39; T reg85=reg11*reg63;
    reg76=reg75+reg76; reg7=reg7/reg30; reg66=reg65+reg66; reg65=reg16*reg64; reg28=reg71+reg28;
    reg71=reg11*reg64; reg68=reg67+reg68; reg67=reg16*reg63; reg75=reg54*reg39; T reg86=reg49*reg42;
    reg41=reg37-reg41; reg35=reg51-reg35; reg37=reg12*reg25; reg84=reg60+reg84; reg51=reg14*reg23;
    reg72=reg41-reg72; reg81=reg59+reg81; reg80=reg80*reg4; reg41=reg16*reg79; reg56=reg77+reg56;
    reg73=reg73*reg21; reg74=reg74*reg13; reg67=reg68+reg67; reg59=reg8+reg55; reg60=reg7*reg62;
    reg65=reg66+reg65; reg85=reg76+reg85; reg69=reg83+reg69; reg66=reg11*reg79; reg68=reg7*reg64;
    reg26=reg78+reg26; reg76=reg22*reg24; reg75=reg86+reg75; reg77=2*reg23; reg78=reg7*reg63;
    reg82=reg70+reg82; reg70=2*reg25; reg83=reg33*reg27; reg86=reg53*reg42; T reg87=reg15*reg39;
    reg71=reg28+reg71; reg28=reg48*reg71; reg74=reg73-reg74; reg73=reg52*reg65; T reg88=reg49*reg71;
    T reg89=reg54*reg65; T reg90=reg14*reg50; T reg91=2*reg50; reg36=reg36*reg61; T reg92=reg22*reg12;
    reg59=reg59*reg4; reg68=reg26+reg68; reg26=reg7*reg79; reg87=reg86+reg87; reg86=reg52*reg67;
    T reg93=reg48*reg85; T reg94=reg49*reg85; reg66=reg69+reg66; reg78=reg75+reg78; reg69=reg70*reg27;
    reg85=reg85*reg51; reg75=reg37*reg67; reg46=reg47+reg46; reg47=reg23*reg25; reg72=reg72/reg30;
    T reg95=reg17*reg83; T reg96=reg44*reg50; T reg97=reg54*reg82; T reg98=reg49*reg81; T reg99=reg81*reg51;
    T reg100=2*reg24; T reg101=reg37*reg82; T reg102=reg76*reg18; T reg103=reg77*reg25; reg60=reg84+reg60;
    reg84=reg76*reg19; T reg104=reg32*reg83; T reg105=reg77*reg24; reg81=reg48*reg81; reg82=reg52*reg82;
    reg67=reg54*reg67; reg41=reg56+reg41; reg80=reg35-reg80; reg35=reg44*reg23; reg56=reg14*reg12;
    T reg106=reg72*reg47; T reg107=reg13*reg92; T reg108=reg72*reg46; T reg109=reg72*reg56; T reg110=reg49*reg66;
    reg90=reg35+reg90; reg59=reg74-reg59; reg35=reg23*reg50; reg74=reg39*reg83; T reg111=reg76*reg42;
    T reg112=reg96*reg60; reg101=reg99+reg101; reg80=reg80/reg30; reg99=reg103*reg10; reg104=reg84+reg104;
    reg84=reg16*reg36; reg19=reg105*reg19; reg32=reg32*reg69; reg26=reg87+reg26; reg18=reg105*reg18;
    reg17=reg17*reg69; reg87=reg54*reg41; reg75=reg85+reg75; reg85=reg96*reg78; T reg113=reg52*reg13;
    T reg114=reg63*reg68; reg89=reg88+reg89; reg88=reg12*reg50; T reg115=reg44*reg25; T reg116=reg54*reg13;
    reg97=reg98+reg97; reg98=reg63*reg60; reg60=reg64*reg60; reg82=reg81+reg82; reg81=reg103*reg13;
    T reg117=reg21*reg92; T reg118=reg24*reg27; T reg119=reg44*reg14; reg100=reg100*reg27; reg71=reg51*reg71;
    reg67=reg94+reg67; reg73=reg28+reg73; reg28=reg64*reg68; reg94=reg52*reg41; T reg120=reg48*reg66;
    T reg121=reg63*reg78; T reg122=reg11*reg36; reg95=reg102+reg95; reg102=reg54*reg10; reg91=reg91*reg61;
    reg65=reg37*reg65; reg86=reg93+reg86; reg78=reg64*reg78; reg93=reg52*reg21; reg66=reg66*reg51;
    T reg123=reg109*reg92; reg41=reg37*reg41; T reg124=reg108*reg103; T reg125=reg44*reg12; reg28=reg73+reg28;
    reg73=reg15*reg4; T reg126=reg52*reg4; T reg127=reg54*reg6; reg98=reg97+reg98; reg107=reg99-reg107;
    reg22=reg44*reg22; reg97=reg25*reg50; reg78=reg86+reg78; reg86=reg106*reg92; reg88=reg115+reg88;
    reg99=reg63*reg26; reg110=reg87+reg110; reg84=reg104+reg84; reg87=reg64*reg26; reg32=reg19+reg32;
    reg16=reg16*reg91; reg19=reg100*reg6; reg94=reg120+reg94; reg104=reg15*reg6; reg11=reg11*reg91;
    reg17=reg18+reg17; reg74=reg111+reg74; reg18=reg7*reg36; reg122=reg95+reg122; reg42=reg105*reg42;
    reg39=reg39*reg69; reg68=reg96*reg68; reg95=reg80*reg90; reg65=reg71+reg65; reg71=reg80*reg35;
    reg111=reg80*reg119; reg115=reg72*reg118; reg121=reg67+reg121; reg67=reg106*reg103; reg113=reg102-reg113;
    reg85=reg75+reg85; reg77=reg77*reg50; reg75=reg4*reg92; reg112=reg101+reg112; reg101=reg108*reg46;
    reg102=reg103*reg6; reg120=reg100*reg4; reg60=reg82+reg60; reg108=reg108*reg92; reg82=reg24*reg61;
    reg106=reg106*reg46; reg81=reg117-reg81; reg30=reg59/reg30; reg116=reg93-reg116; reg114=reg89+reg114;
    reg59=reg109*reg103; reg104=reg113-reg104; reg89=reg49*reg13; reg73=reg116-reg73; reg120=reg81-reg120;
    reg13=reg48*reg13; reg10=reg49*reg10; reg81=(*f.m).alpha_2*reg52; reg33=reg44*reg33; reg93=reg27*reg61;
    reg45=reg43+reg45; reg43=reg15*reg5; reg127=reg126+reg127; reg70=reg70*reg50; reg100=reg100*reg5;
    reg102=reg75+reg102; reg75=(*f.m).alpha_1*reg48; reg113=reg115*reg92; reg108=reg60+reg108; reg60=reg95*reg22;
    reg116=reg48*reg122; reg117=reg52*reg84; reg59=reg114+reg59; reg114=reg111*reg77; reg109=reg109*reg46;
    reg68=reg65+reg68; reg67=reg121+reg67; reg65=reg71*reg77; reg121=reg14*reg61; reg126=reg44*reg24;
    T reg128=(*f.m).alpha_1*reg49; T reg129=(*f.m).alpha_2*reg54; reg19=reg107-reg19; reg107=reg80*reg82; reg11=reg17+reg11;
    reg17=reg115*reg103; T reg130=reg30*reg125; T reg131=reg30*reg97; reg7=reg7*reg91; reg106=reg85+reg106;
    reg124=reg98+reg124; reg85=reg95*reg77; reg39=reg42+reg39; reg18=reg74+reg18; reg42=reg30*reg88;
    reg74=reg49*reg122; reg98=reg54*reg84; reg16=reg32+reg16; reg32=reg71*reg90; reg41=reg66+reg41;
    reg26=reg96*reg26; reg95=reg95*reg90; reg101=reg112+reg101; reg99=reg110+reg99; reg66=reg111*reg22;
    reg21=reg48*reg21; reg87=reg94+reg87; reg123=reg28+reg123; reg86=reg78+reg86; reg71=reg71*reg22;
    reg28=reg49*reg104; reg117=reg116+reg117; reg78=reg42*reg33; reg60=reg108+reg60; reg94=reg107*reg22;
    reg113=reg87+reg113; reg84=reg37*reg84; reg122=reg51*reg122; reg102=reg100-reg102; reg87=reg49*reg19;
    reg100=reg42*reg88; reg95=reg101+reg95; reg101=reg48*reg120; reg108=reg54*reg16; reg110=reg49*reg11;
    reg129=reg128+reg129; reg112=reg19*reg47; reg109=reg68+reg109; reg68=reg120*reg56; reg111=reg111*reg90;
    reg116=reg72*reg45; reg128=reg30*reg93; reg121=reg126+reg121; reg7=reg39+reg7; reg65=reg67+reg65;
    reg39=reg131*reg70; reg67=reg23*reg61; reg126=reg24*reg50; T reg132=reg130*reg70; reg114=reg59+reg114;
    reg59=reg73*reg56; T reg133=reg104*reg47; T reg134=reg48*reg73; T reg135=reg52*reg16; T reg136=reg48*reg11;
    T reg137=reg18*reg64; T reg138=reg18*reg63; reg29=reg31+reg29; reg98=reg74+reg98; reg42=reg42*reg70;
    reg13=reg10-reg13; reg85=reg124+reg85; reg10=reg53*reg6; reg31=reg107*reg77; reg17=reg99+reg17;
    reg74=reg44*reg27; reg99=reg12*reg61; reg124=reg52*reg73; T reg139=reg54*reg104; reg89=reg21-reg89;
    reg21=reg53*reg4; T reg140=reg131*reg33; reg71=reg86+reg71; reg86=reg54*reg19; T reg141=reg52*reg120;
    T reg142=reg130*reg33; reg66=reg123+reg66; reg123=reg64*(*f.m).alpha_3; reg75=reg81+reg75; reg81=reg63*(*f.m).alpha_3;
    reg115=reg115*reg46; reg26=reg41+reg26; reg41=(*f.m).alpha_1*reg53; T reg143=(*f.m).alpha_2*reg15; reg127=reg43-reg127;
    reg6=reg49*reg6; reg32=reg106+reg32; reg131=reg131*reg88; reg4=reg48*reg4; reg99=reg74+reg99;
    reg43=reg27*reg50; reg72=reg72*reg29; reg74=reg25*reg61; reg140=reg71+reg140; reg137=reg117+reg137;
    reg39=reg65+reg39; reg65=reg7*reg63; reg108=reg110+reg108; reg67=reg126+reg67; reg142=reg66+reg142;
    reg78=reg60+reg78; reg14=reg14*reg24; reg60=reg53*reg102; reg87=reg101+reg87; reg12=reg12*reg27;
    reg66=reg7*reg64; reg42=reg85+reg42; reg71=reg128*reg33; reg94=reg113+reg94; reg135=reg136+reg135;
    reg5=reg53*reg5; reg85=reg128*reg70; reg31=reg17+reg31; reg17=reg2*reg46; reg132=reg114+reg132;
    reg101=reg37*reg2; reg6=reg4+reg6; reg138=reg98+reg138; reg130=reg130*reg88; reg111=reg109+reg111;
    reg4=reg116*reg103; reg98=reg80*reg121; reg106=reg116*reg92; reg18=reg18*reg96; reg84=reg122+reg84;
    reg100=reg95+reg100; reg123=reg75+reg123; reg107=reg107*reg90; reg115=reg26+reg115; reg112=reg68+reg112;
    reg81=reg129+reg81; reg26=reg118*reg102; reg143=reg41+reg143; reg41=reg79*(*f.m).alpha_3; reg68=(*f.m).alpha_1*reg51;
    reg75=(*f.m).alpha_2*reg37; reg131=reg32+reg131; reg10=reg13-reg10; reg13=reg118*reg127; reg32=reg15*reg127;
    reg86=reg141+reg86; reg95=reg15*reg102; reg109=reg53*reg127; reg28=reg134+reg28; reg139=reg124+reg139;
    reg21=reg89-reg21; reg133=reg59+reg133; reg11=reg51*reg11; reg16=reg37*reg16; reg16=reg11+reg16;
    reg92=reg72*reg92; reg66=reg135+reg66; reg11=reg1*reg45; reg37=reg48*reg21; reg59=reg49*reg10;
    reg89=reg78*reg131; reg18=reg84+reg18; reg84=reg98*reg22; reg110=reg39*reg100; reg116=reg116*reg46;
    reg71=reg94+reg71; reg94=reg30*reg99; reg106=reg137+reg106; reg113=reg42*reg131; reg114=reg140*reg100;
    reg7=reg7*reg96; reg32=reg139+reg32; reg117=reg57*reg101; reg74=reg43+reg74; reg95=reg86+reg95;
    reg43=reg57*reg17; reg96=reg96*(*f.m).alpha_3; reg86=(*f.m).alpha_1*reg14; reg85=reg31+reg85; reg31=(*f.m).alpha_2*reg12;
    reg75=reg68+reg75; reg6=reg5-reg6; reg41=reg143+reg41; reg4=reg138+reg4; reg5=reg98*reg77;
    reg27=reg25*reg27; reg24=reg23*reg24; reg44=reg44*reg61; reg107=reg115+reg107; reg128=reg128*reg88;
    reg65=reg108+reg65; reg103=reg72*reg103; reg12=reg12*reg1; reg23=reg123*reg132; reg25=reg46*reg101;
    reg68=reg81*reg39; reg108=reg58*reg17; reg115=reg81*reg140; reg13=reg133+reg13; reg60=reg87+reg60;
    reg87=reg123*reg142; reg109=reg28+reg109; reg80=reg80*reg67; reg28=reg58*reg101; reg130=reg111+reg130;
    reg26=reg112+reg26; reg111=reg46*reg17; reg112=reg54*reg10; reg122=reg52*reg21; reg30=reg30*reg74;
    reg98=reg98*reg90; reg116=reg18+reg116; reg113=reg110-reg113; reg18=reg10*reg47; reg110=reg21*reg56;
    reg124=(*f.m).alpha_2*reg27; reg126=(*f.m).alpha_1*reg24; reg71=reg41*reg71; reg128=reg107+reg128; reg112=reg122+reg112;
    reg103=reg65+reg103; reg77=reg80*reg77; reg2=reg51*reg2; reg15=reg15*reg6; reg51=reg94*reg70;
    reg5=reg4+reg5; reg4=reg45*reg12; reg27=reg27*reg0; reg65=reg0*reg29; reg25=reg13+reg25;
    reg50=reg61*reg50; reg115=reg87+reg115; reg59=reg37+reg59; reg53=reg53*reg6; reg28=reg109+reg28;
    reg13=reg76*reg12; reg117=reg32+reg117; reg32=reg83*reg12; reg37=reg78*reg39; reg61=reg140*reg42;
    reg87=reg123*reg130; reg107=reg81*reg131; reg43=reg95+reg43; reg95=reg83*reg11; reg108=reg60+reg108;
    reg22=reg80*reg22; reg92=reg66+reg92; reg60=reg76*reg11; reg96=reg75+reg96; reg72=reg72*reg46;
    reg66=reg45*reg11; reg44=reg44*(*f.m).alpha_3; reg7=reg16+reg7; reg68=reg23+reg68; reg111=reg26+reg111;
    reg89=reg114-reg89; reg16=reg94*reg33; reg84=reg106+reg84; reg31=reg86+reg31; reg85=reg41*reg85;
    reg124=reg126+reg124; reg23=reg105*reg27; reg26=reg96*reg78; reg75=reg69*reg27; reg32=reg117+reg32;
    reg13=reg28+reg13; reg44=reg31+reg44; reg51=reg5+reg51; reg4=reg25+reg4; reg5=reg105*reg65;
    reg95=reg43+reg95; reg25=reg69*reg65; reg50=reg50*(*f.m).alpha_3; reg66=reg111+reg66; reg28=reg29*reg65;
    reg60=reg108+reg60; reg118=reg118*reg6; reg18=reg110+reg18; reg31=reg96*reg42; reg80=reg80*reg90;
    reg72=reg7+reg72; reg16=reg84+reg16; reg85=reg68+reg85; reg22=reg92+reg22; reg33=reg30*reg33;
    reg94=reg94*reg88; reg98=reg116+reg98; reg7=reg89*reg132; reg43=reg142*reg113; reg128=reg41*reg128;
    reg68=reg29*reg27; reg107=reg87+reg107; reg37=reg61-reg37; reg71=reg115+reg71; reg77=reg103+reg77;
    reg15=reg112+reg15; reg1=reg14*reg1; reg14=reg57*reg2; reg70=reg30*reg70; reg53=reg59+reg53;
    reg59=reg58*reg2; reg128=reg107+reg128; reg61=reg42*reg130; reg80=reg72+reg80; reg72=reg132*reg100;
    reg84=reg96*reg100; reg30=reg30*reg88; reg25=reg95+reg25; reg31=reg85+reg31; reg85=reg142*reg100;
    reg33=reg22+reg33; reg7=reg43-reg7; reg22=reg37*reg130; reg23=reg13+reg23; reg51=reg44*reg51;
    reg75=reg32+reg75; reg26=reg71+reg26; reg14=reg15+reg14; reg0=reg24*reg0; reg16=reg44*reg16;
    reg83=reg83*reg1; reg59=reg53+reg59; reg76=reg76*reg1; reg70=reg77+reg70; reg118=reg18+reg118;
    reg13=reg46*reg2; reg50=reg124+reg50; reg68=reg4+reg68; reg94=reg98+reg94; reg4=reg130*reg78;
    reg5=reg60+reg5; reg28=reg66+reg28; reg13=reg118+reg13; reg105=reg105*reg0; reg76=reg59+reg76;
    reg45=reg45*reg1; reg30=reg80+reg30; reg4=reg85-reg4; reg70=reg50*reg70; reg15=reg5*reg68;
    reg18=reg140*reg130; reg24=reg78*reg132; reg32=reg142*reg131; reg43=reg142*reg42; reg51=reg31+reg51;
    reg83=reg14+reg83; reg14=reg130*reg39; reg61=reg72-reg61; reg31=reg132*reg131; reg53=reg25*reg68;
    reg16=reg26+reg16; reg84=reg128+reg84; reg33=reg50*reg33; reg22=reg7+reg22; reg94=reg44*reg94;
    reg7=reg75*reg28; reg69=reg69*reg0; reg26=reg23*reg28; reg18=reg32-reg18; reg32=reg25*reg23;
    reg59=reg75*reg5; reg4=reg4/reg22; reg45=reg13+reg45; reg15=reg26-reg15; reg89=reg89/reg22;
    reg14=reg31-reg14; reg83=reg69+reg83; reg61=reg61/reg22; reg113=reg113/reg22; reg33=reg16+reg33;
    reg94=reg84+reg94; reg13=reg140*reg132; reg30=reg50*reg30; reg70=reg51+reg70; reg24=reg43-reg24;
    reg16=reg142*reg39; reg29=reg29*reg0; reg105=reg76+reg105; reg53=reg7-reg53; reg7=reg105*reg53;
    reg4=reg4*reg70; reg13=reg16-reg13; reg61=reg61*reg33; reg29=reg45+reg29; reg14=reg14/reg22;
    reg37=reg37/reg22; reg89=reg89*reg70; reg24=reg24/reg22; reg16=reg15*reg83; reg113=reg113*reg33;
    reg59=reg32-reg59; reg30=reg94+reg30; reg18=reg18/reg22; reg89=reg113-reg89; reg24=reg24*reg30;
    reg61=reg4-reg61; reg33=reg14*reg33; reg16=reg7-reg16; reg70=reg18*reg70; reg4=reg59*reg29;
    reg7=reg105*reg68; reg22=reg13/reg22; reg13=reg23*reg29; reg14=reg75*reg29; reg18=reg83*reg68;
    reg37=reg37*reg30; reg26=reg5*reg29; reg31=elem.pos(2)[0]-elem.pos(0)[0]; reg24=reg61-reg24; reg14=reg18-reg14;
    reg30=reg22*reg30; reg70=reg33-reg70; reg18=elem.pos(1)[0]-elem.pos(0)[0]; reg22=1-(*f.m).resolution; reg32=elem.pos(2)[1]-elem.pos(0)[1];
    reg33=reg105*reg28; reg43=elem.pos(1)[1]-elem.pos(0)[1]; reg45=reg83*reg23; reg4=reg16+reg4; reg16=reg83*reg28;
    reg89=reg37+reg89; reg13=reg7-reg13; reg7=reg25*reg29; reg37=reg75*reg105; reg24=reg24*reg22;
    reg70=reg30+reg70; reg30=reg83*reg5; reg51=reg31*reg43; reg45=reg37-reg45; reg37=reg32*reg18;
    reg60=reg81*(*f.m).resolution; reg61=reg25*reg105; reg7=reg16-reg7; reg14=reg14/reg4; reg16=reg123*(*f.m).resolution;
    reg89=reg89*reg22; reg13=reg13/reg4; reg26=reg33-reg26; reg33=reg96*(*f.m).resolution; reg131=reg22*reg131;
    reg130=reg130*reg22; reg53=reg53/reg4; reg7=reg7/reg4; reg15=reg15/reg4; reg26=reg26/reg4;
    reg51=reg37-reg51; reg45=reg45/reg4; reg30=reg61-reg30; reg24=reg60+reg24; reg89=reg16+reg89;
    reg70=reg70*reg22; reg16=reg14*(*f.m).resolution; reg37=reg13*(*f.m).resolution; reg31=reg31/reg51; reg59=reg59/reg4;
    reg32=reg32/reg51; reg60=reg53*(*f.m).resolution; reg61=reg15*(*f.m).resolution; reg43=reg43/reg51; reg18=reg18/reg51;
    reg4=reg30/reg4; reg30=reg26*(*f.m).resolution; reg39=reg22*reg39; reg132=reg22*reg132; reg140=reg140*reg22;
    reg142=reg142*reg22; reg37=reg131-reg37; reg16=reg130+reg16; reg66=reg7*(*f.m).resolution; reg100=reg22*reg100;
    reg69=reg45*(*f.m).resolution; reg70=reg33+reg70; reg24=(*f.m).deltaT*reg24; reg89=(*f.m).deltaT*reg89; reg33=reg43-reg32;
    reg71=reg16*reg89; reg72=reg37*reg24; reg69=reg100+reg69; reg70=(*f.m).deltaT*reg70; reg76=reg59*(*f.m).resolution;
    reg39=reg30+reg39; reg66=reg132-reg66; reg30=reg4*(*f.m).resolution; reg78=reg22*reg78; reg22=reg42*reg22;
    reg61=reg140-reg61; reg142=reg60+reg142; reg42=reg31-reg18; reg60=reg61*reg24; reg77=reg142*reg89;
    reg80=reg71+reg72; reg84=reg69*reg70; reg85=0.5*reg43; reg86=0.5*reg18; reg87=0.5*reg42;
    reg30=reg22-reg30; reg22=reg39*reg24; reg92=0.5*reg33; reg94=reg66*reg89; reg95=0.5*reg31;
    reg98=0.5*reg32; reg78=reg76+reg78; reg76=reg69*reg95; reg100=reg87*reg69; reg103=reg33*reg16;
    reg106=reg16*reg32; reg107=reg69*reg85; reg108=reg37*reg18; reg109=reg60+reg77; reg110=reg37*reg31;
    reg111=reg98*reg69; reg112=reg69*reg86; reg113=reg16*reg43; reg114=reg22+reg94; reg115=reg78*reg70;
    reg116=reg30*reg70; reg117=reg80+reg84; reg118=reg42*reg37; reg122=reg69*reg92; reg124=reg109+reg115;
    reg126=reg30*reg86; reg128=reg39*reg18; reg112=reg112-reg113; reg100=reg103+reg100; reg122=reg118+reg122;
    reg103=2*reg117; reg118=reg114+reg116; reg129=reg30*reg92; reg130=reg42*reg39; reg131=reg87*reg78;
    reg132=reg61*reg18; reg133=reg78*reg85; reg134=reg33*reg142; reg135=reg98*reg30; reg136=reg39*reg31;
    reg137=reg66*reg32; reg138=reg30*reg95; reg139=reg87*reg30; reg140=reg142*reg32; reg141=reg78*reg95;
    reg143=reg33*reg66; reg111=reg111-reg110; T reg144=reg142*reg43; reg106=reg106-reg76; T reg145=reg78*reg86;
    reg108=reg108-reg107; T reg146=reg61*reg31; T reg147=reg98*reg78; T reg148=reg66*reg43; T reg149=reg30*reg85;
    T reg150=reg78*reg92; T reg151=reg42*reg61; T reg152=reg103*reg95; reg128=reg128-reg149; T reg153=reg103*reg85;
    T reg154=reg124*reg32; reg112=2*reg112; reg126=reg126-reg148; reg132=reg132-reg133; reg145=reg145-reg144;
    reg108=2*reg108; reg139=reg143+reg139; reg129=reg130+reg129; reg111=2*reg111; reg140=reg140-reg141;
    reg106=2*reg106; reg147=reg147-reg146; reg150=reg151+reg150; reg100=2*reg100; reg122=2*reg122;
    reg131=reg134+reg131; reg130=reg118*reg31; reg134=reg103*reg86; reg135=reg135-reg136; reg143=reg124*reg43;
    reg137=reg137-reg138; reg151=reg118*reg18; T reg155=reg103*reg98; T reg156=reg135*reg18; T reg157=reg112*reg92;
    T reg158=reg100*reg86; T reg159=reg131*reg43; T reg160=reg132*reg43; T reg161=reg103*reg92; T reg162=reg42*reg128;
    T reg163=reg86*reg108; T reg164=reg122*reg86; T reg165=reg150*reg43; T reg166=reg145*reg43; T reg167=reg86*reg106;
    T reg168=reg42*reg129; T reg169=reg140*reg43; T reg170=reg147*reg43; T reg171=reg86*reg111; T reg172=reg92*reg108;
    T reg173=reg85*reg106; T reg174=reg137*reg18; T reg175=reg85*reg112; T reg176=reg100*reg95; T reg177=reg122*reg92;
    T reg178=reg126*reg18; T reg179=reg150*reg32; T reg180=reg131*reg32; T reg181=reg86*reg112; T reg182=reg103*reg87;
    T reg183=reg42*reg126; T reg184=reg129*reg18; T reg185=reg33*reg124; T reg186=reg85*reg111; T reg187=reg100*reg85;
    T reg188=reg139*reg18; T reg189=reg92*reg106; T reg190=reg137*reg42; T reg191=reg143-reg134; T reg192=reg111*reg92;
    T reg193=reg135*reg42; T reg194=reg130-reg155; reg131=reg33*reg131; T reg195=reg122*reg87; T reg196=reg100*reg87;
    T reg197=reg122*reg95; reg150=reg33*reg150; T reg198=reg122*reg85; T reg199=reg33*reg147; reg122=reg98*reg122;
    reg129=reg129*reg31; T reg200=reg98*reg106; T reg201=reg87*reg106; reg137=reg137*reg31; T reg202=reg98*reg111;
    T reg203=reg33*reg140; reg135=reg135*reg31; T reg204=reg98*reg112; T reg205=reg128*reg31; T reg206=reg98*reg108;
    reg126=reg126*reg31; T reg207=reg87*reg111; T reg208=reg42*reg139; T reg209=reg100*reg92; T reg210=reg87*reg108;
    T reg211=reg33*reg145; T reg212=reg33*reg132; reg139=reg139*reg31; reg100=reg98*reg100; T reg213=reg87*reg112;
    T reg214=reg95*reg108; reg132=reg132*reg32; T reg215=reg118*reg42; reg112=reg112*reg95; reg145=reg145*reg32;
    reg111=reg111*reg95; reg147=reg147*reg32; reg128=reg128*reg18; reg108=reg85*reg108; reg106=reg106*reg95;
    reg140=reg140*reg32; T reg216=reg152-reg154; T reg217=reg153-reg151; T reg218=reg161+reg215; reg129=reg122-reg129;
    reg122=reg185+reg182; reg159=reg158-reg159; reg169=reg167-reg169; reg216=reg216*reg51; reg217=reg217*reg51;
    reg186=reg156-reg186; reg191=reg191*reg51; reg175=reg178-reg175; reg194=reg194*reg51; reg168=reg177+reg168;
    reg195=reg150+reg195; reg139=reg100-reg139; reg208=reg209+reg208; reg210=reg212+reg210; reg165=reg164-reg165;
    reg137=reg200-reg137; reg214=reg132-reg214; reg135=reg202-reg135; reg112=reg145-reg112; reg126=reg204-reg126;
    reg111=reg147-reg111; reg205=reg206-reg205; reg106=reg140-reg106; reg201=reg203+reg201; reg190=reg189+reg190;
    reg193=reg192+reg193; reg196=reg131+reg196; reg162=reg172+reg162; reg170=reg171-reg170; reg176=reg180-reg176;
    reg197=reg179-reg197; reg166=reg181-reg166; reg183=reg157+reg183; reg160=reg163-reg160; reg199=reg207+reg199;
    reg213=reg211+reg213; reg187=reg188-reg187; reg173=reg174-reg173; reg108=reg128-reg108; reg198=reg184-reg198;
    reg162=reg162*reg51; reg196=reg196*reg51; reg201=reg201*reg51; reg217=ponderation*reg217; reg216=ponderation*reg216;
    reg176=reg176*reg51; reg205=reg205*reg51; reg197=reg197*reg51; reg126=reg126*reg51; reg191=ponderation*reg191;
    reg183=reg183*reg51; reg135=reg135*reg51; reg100=reg122*reg51; reg137=reg137*reg51; reg175=reg175*reg51;
    reg194=ponderation*reg194; reg128=reg218*reg51; reg129=reg129*reg51; reg139=reg139*reg51; reg199=reg199*reg51;
    reg195=reg195*reg51; reg213=reg213*reg51; reg198=reg198*reg51; reg173=reg173*reg51; reg187=reg187*reg51;
    reg186=reg186*reg51; reg108=reg108*reg51; reg160=reg160*reg51; reg168=reg168*reg51; reg208=reg208*reg51;
    reg210=reg210*reg51; reg159=reg159*reg51; reg214=reg214*reg51; reg166=reg166*reg51; reg193=reg193*reg51;
    reg169=reg169*reg51; reg190=reg190*reg51; reg170=reg170*reg51; reg106=reg106*reg51; reg112=reg112*reg51;
    reg165=reg165*reg51; reg111=reg111*reg51; T tmp_5_0=ponderation*reg187; reg131=ponderation*reg128; T vec_1=reg131;
    T tmp_4_2=ponderation*reg169; T tmp_5_5=ponderation*reg108; T tmp_0_4=ponderation*reg213; T tmp_1_4=ponderation*reg183; T tmp_1_5=ponderation*reg162;
    T tmp_4_5=ponderation*reg160; T tmp_0_3=ponderation*reg199; T tmp_2_0=ponderation*reg176; T tmp_3_5=ponderation*reg205; T tmp_4_1=ponderation*reg165;
    T tmp_4_0=ponderation*reg159; T tmp_4_4=ponderation*reg166; T tmp_4_3=ponderation*reg170; T tmp_1_3=ponderation*reg193; reg108=ponderation*reg100;
    T vec_0=reg108; T tmp_1_2=ponderation*reg190; T tmp_0_0=ponderation*reg196; T tmp_0_2=ponderation*reg201; T tmp_2_1=ponderation*reg197;
    T tmp_2_2=ponderation*reg106; T tmp_3_4=ponderation*reg126; T tmp_2_3=ponderation*reg111; T tmp_3_3=ponderation*reg135; T tmp_2_4=ponderation*reg112;
    T tmp_3_2=ponderation*reg137; T tmp_2_5=ponderation*reg214; T tmp_3_1=ponderation*reg129; T tmp_3_0=ponderation*reg139; T tmp_0_5=ponderation*reg210;
    T tmp_1_0=ponderation*reg208; T tmp_0_1=ponderation*reg195; T tmp_1_1=ponderation*reg168; T vec_3=-reg194; T tmp_5_3=ponderation*reg186;
    T vec_4=-reg191; T vec_2=-reg216; T tmp_5_2=ponderation*reg173; T tmp_5_1=ponderation*reg198; T vec_5=-reg217;
    T tmp_5_4=ponderation*reg175;
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
    T reg3=reg1*reg0; reg2=1.0/reg2; T reg4=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg5=1.0/(*f.m).elastic_modulus_3; T reg6=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg7=reg2*reg3; T reg8=pow((*f.m).v1[0],2); T reg9=pow((*f.m).v2[1],2); T reg10=pow((*f.m).v2[0],2); T reg11=pow((*f.m).v1[1],2);
    T reg12=pow((*f.m).v2[2],2); reg9=reg10+reg9; reg10=1.0/(*f.m).elastic_modulus_2; T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg14=pow((*f.m).v1[2],2);
    reg11=reg8+reg11; reg8=reg5*reg7; T reg15=reg6*reg7; T reg16=reg4*reg7; reg14=reg11+reg14;
    reg11=reg10*reg8; T reg17=reg6*reg16; T reg18=reg6*reg15; T reg19=reg13*reg8; reg12=reg9+reg12;
    reg14=pow(reg14,0.5); reg17=reg17+reg19; reg9=reg13*reg15; T reg20=1.0/(*f.m).elastic_modulus_1; T reg21=reg10*reg16;
    reg18=reg11-reg18; reg12=pow(reg12,0.5); reg11=(*f.m).v1[1]/reg14; T reg22=(*f.m).v1[2]/reg14; T reg23=(*f.m).v2[1]/reg12;
    T reg24=reg20*reg18; T reg25=(*f.m).v2[2]/reg12; T reg26=reg9+reg21; T reg27=reg17*reg13; reg12=(*f.m).v2[0]/reg12;
    T reg28=reg6*reg3; reg14=(*f.m).v1[0]/reg14; T reg29=reg4*reg16; T reg30=reg4*reg3; reg8=reg20*reg8;
    T reg31=reg13*reg7; reg3=reg5*reg3; reg27=reg24-reg27; reg24=reg2*reg0; T reg32=reg22*reg23;
    T reg33=reg4*reg15; T reg34=reg11*reg25; reg7=reg10*reg7; T reg35=reg26*reg4; reg16=reg13*reg16;
    T reg36=reg24*reg6; T reg37=reg31*reg4; T reg38=reg2*reg1; T reg39=reg24*reg4; T reg40=reg3*reg10;
    reg3=reg3*reg13; reg24=reg24*reg5; T reg41=reg28*reg6; T reg42=reg30*reg6; reg35=reg27-reg35;
    reg27=reg7*reg4; T reg43=2*reg12; reg19=reg33+reg19; reg33=2*reg14; T reg44=reg34-reg32;
    T reg45=reg14*reg25; T reg46=reg22*reg12; reg29=reg8-reg29; reg15=reg20*reg15; reg31=reg31*reg13;
    reg29=reg29/reg35; reg16=reg15+reg16; reg28=reg28*reg13; reg17=reg17/reg35; reg8=reg38*reg6;
    reg30=reg30*reg10; reg27=reg9+reg27; T reg47=reg38*reg4; reg19=reg19/reg35; reg37=reg15+reg37;
    reg41=reg40-reg41; reg15=pow(reg11,2); reg7=reg7*reg20; reg40=reg43*reg23; T reg48=pow(reg23,2);
    T reg49=pow(reg12,2); T reg50=pow(reg14,2); reg18=reg18/reg35; reg38=reg38*reg5; T reg51=reg11*reg12;
    T reg52=reg39*reg6; T reg53=reg14*reg23; T reg54=2*reg44; T reg55=reg36*reg6; T reg56=reg24*reg13;
    reg24=reg24*reg10; T reg57=reg33*reg11; T reg58=reg46-reg45; reg42=reg3+reg42; reg3=reg53-reg51;
    reg26=reg26/reg35; T reg59=reg50*reg19; T reg60=pow(reg22,2); T reg61=pow(reg25,2); reg27=reg27/reg35;
    T reg62=pow(reg44,2); T reg63=pow(reg58,2); T reg64=reg54*reg58; T reg65=reg57*reg19; T reg66=reg40*reg29;
    reg42=reg42*reg13; reg41=reg41*reg20; T reg67=reg28+reg30; reg55=reg24-reg55; reg39=reg39*reg10;
    reg36=reg36*reg13; reg52=reg56+reg52; reg24=reg38*reg10; reg38=reg38*reg13; reg56=reg8*reg6;
    reg31=reg7-reg31; reg7=reg47*reg6; T reg68=reg15*reg18; T reg69=reg48*reg17; reg16=reg16/reg35;
    T reg70=reg57*reg18; T reg71=reg40*reg17; T reg72=reg48*reg29; T reg73=reg15*reg19; reg37=reg37/reg35;
    T reg74=reg49*reg29; T reg75=reg49*reg17; T reg76=reg50*reg18; reg69=reg68+reg69; reg7=reg38+reg7;
    reg56=reg24-reg56; reg24=reg50*reg27; reg38=reg36+reg39; reg52=reg52*reg13; reg55=reg55*reg20;
    reg67=reg67*reg4; reg68=reg16*reg64; reg75=reg76+reg75; reg42=reg41-reg42; reg47=reg47*reg10;
    reg41=reg26*reg62; reg66=reg65+reg66; reg65=reg48*reg37; reg76=reg15*reg27; reg8=reg8*reg13;
    T reg77=reg40*reg37; T reg78=reg57*reg27; T reg79=pow(reg3,2); reg31=reg31/reg35; T reg80=reg49*reg37;
    reg74=reg59+reg74; reg59=reg16*reg62; reg72=reg73+reg72; reg73=reg16*reg63; T reg81=reg60*reg19;
    T reg82=reg61*reg29; T reg83=reg26*reg64; reg71=reg70+reg71; reg70=reg61*reg17; T reg84=reg26*reg63;
    T reg85=reg60*reg18; reg41=reg75+reg41; reg75=reg16*reg79; reg65=reg76+reg65; reg7=reg7*reg13;
    reg84=reg69+reg84; reg69=reg12*reg23; reg76=2*reg11; T reg86=reg14*reg11; reg82=reg81+reg82;
    reg81=reg31*reg64; reg73=reg72+reg73; reg77=reg78+reg77; reg72=reg8+reg47; reg59=reg74+reg59;
    reg74=reg33*reg22; reg78=reg61*reg37; T reg87=reg31*reg62; T reg88=reg60*reg27; reg80=reg24+reg80;
    reg24=2*reg23; T reg89=reg43*reg25; T reg90=reg31*reg63; reg67=reg42-reg67; reg56=reg56*reg20;
    reg52=reg55-reg52; reg42=reg26*reg79; reg68=reg66+reg68; reg38=reg38*reg4; reg70=reg85+reg70;
    reg83=reg71+reg83; reg42=reg70+reg42; reg55=reg74*reg19; reg72=reg72*reg4; reg87=reg80+reg87;
    reg66=reg31*reg79; reg70=reg76*reg23; reg78=reg88+reg78; reg38=reg52-reg38; reg53=reg51+reg53;
    reg51=reg29*reg89; reg52=reg11*reg23; reg71=reg14*reg12; reg80=2*reg58; reg54=reg54*reg3;
    reg90=reg65+reg90; reg65=reg33*reg12; reg85=reg48*reg68; reg88=reg24*reg25; T reg91=reg15*reg83;
    T reg92=reg15*reg84; T reg93=reg48*reg73; T reg94=reg44*reg11; reg81=reg77+reg81; reg77=reg14*reg58;
    T reg95=reg17*reg89; T reg96=reg49*reg73; reg67=reg67/reg35; T reg97=reg50*reg84; T reg98=reg49*reg59;
    reg75=reg82+reg75; reg82=reg74*reg18; T reg99=reg50*reg41; T reg100=reg49*reg68; T reg101=reg15*reg41;
    T reg102=reg83*reg86; reg68=reg69*reg68; T reg103=2*reg22; reg7=reg56-reg7; reg73=reg69*reg73;
    reg84=reg84*reg86; reg56=reg76*reg22; T reg104=reg48*reg59; T reg105=reg44*reg58; reg83=reg50*reg83;
    reg51=reg55+reg51; reg55=reg67*reg53; reg98=reg99+reg98; reg99=reg16*reg54; reg103=reg103*reg25;
    reg19=reg56*reg19; T reg106=reg67*reg52; reg72=reg7-reg72; reg17=reg17*reg88; reg18=reg56*reg18;
    reg7=reg26*reg54; reg95=reg82+reg95; reg38=reg38/reg35; reg82=reg37*reg89; T reg107=reg11*reg58;
    T reg108=reg74*reg27; reg41=reg86*reg41; T reg109=reg70*reg10; T reg110=reg13*reg65; T reg111=reg48*reg10;
    reg59=reg69*reg59; T reg112=reg49*reg13; T reg113=reg70*reg13; T reg114=reg20*reg65; T reg115=reg48*reg13;
    T reg116=reg49*reg20; reg73=reg84+reg73; reg84=reg105*reg90; reg77=reg94+reg77; reg94=reg105*reg81;
    reg68=reg102+reg68; reg66=reg78+reg66; reg78=reg67*reg71; reg80=reg80*reg3; reg102=reg22*reg25;
    T reg117=reg44*reg14; T reg118=reg15*reg42; reg100=reg83+reg100; reg83=reg62*reg81; reg104=reg101+reg104;
    reg101=reg63*reg87; reg81=reg63*reg81; reg85=reg91+reg85; reg93=reg92+reg93; reg91=reg63*reg90;
    reg92=reg48*reg75; T reg119=reg49*reg75; T reg120=reg50*reg42; reg90=reg62*reg90; reg96=reg97+reg96;
    reg97=reg62*reg87; reg29=reg29*reg88; T reg121=reg44*reg23; T reg122=reg12*reg58; T reg123=reg78*reg70;
    reg101=reg104+reg101; reg104=reg22*reg3; reg110=reg109-reg110; reg109=reg103*reg6; reg112=reg111-reg112;
    reg111=reg61*reg6; T reg124=reg55*reg65; reg83=reg100+reg83; reg59=reg41+reg59; reg87=reg105*reg87;
    reg113=reg114-reg113; reg41=reg103*reg4; reg122=reg121+reg122; reg115=reg116-reg115; reg100=reg61*reg4;
    reg114=reg70*reg6; reg116=reg55*reg53; reg94=reg68+reg94; reg68=reg4*reg65; reg97=reg98+reg97;
    reg98=reg78*reg65; reg26=reg26*reg80; reg17=reg18+reg17; reg99=reg51+reg99; reg29=reg19+reg29;
    reg90=reg96+reg90; reg18=reg106*reg65; reg19=reg106*reg70; reg91=reg93+reg91; reg16=reg16*reg80;
    reg119=reg120+reg119; reg51=reg62*reg66; reg118=reg92+reg118; reg92=reg63*reg66; reg106=reg106*reg53;
    reg35=reg72/reg35; reg84=reg73+reg84; reg72=reg44*reg12; reg7=reg95+reg7; reg73=reg23*reg58;
    reg81=reg85+reg81; reg55=reg55*reg70; reg85=reg31*reg54; reg27=reg56*reg27; reg37=reg37*reg88;
    reg82=reg108+reg82; reg93=reg49*reg4; reg76=reg76*reg58; reg75=reg69*reg75; reg95=reg67*reg102;
    reg96=reg48*reg6; reg42=reg42*reg86; reg108=reg38*reg77; reg120=reg38*reg107; reg121=reg38*reg117;
    reg33=reg44*reg33; T reg125=reg15*reg7; T reg126=(*f.m).alpha_2*reg48; T reg127=(*f.m).alpha_1*reg15; T reg128=reg48*reg99;
    reg123=reg101+reg123; reg101=reg121*reg76; T reg129=reg35*reg122; T reg130=reg108*reg76; reg55=reg81+reg55;
    reg46=reg45+reg46; reg45=reg25*reg3; reg81=reg49*reg99; T reg131=reg50*reg7; reg109=reg110-reg109;
    reg110=reg35*reg73; T reg132=reg35*reg72; reg26=reg17+reg26; reg98=reg97+reg98; reg17=reg121*reg33;
    reg97=reg120*reg76; T reg133=reg38*reg104; reg18=reg90+reg18; reg90=reg120*reg33; reg19=reg91+reg19;
    reg24=reg24*reg58; reg16=reg29+reg16; reg29=(*f.m).alpha_2*reg49; reg51=reg119+reg51; reg106=reg84+reg106;
    reg92=reg118+reg92; reg84=reg95*reg70; reg31=reg31*reg80; reg37=reg27+reg37; reg85=reg82+reg85;
    reg27=(*f.m).alpha_1*reg50; reg82=reg95*reg65; reg114=reg68+reg114; reg103=reg103*reg5; reg68=reg15*reg13;
    reg91=reg44*reg22; reg78=reg78*reg53; reg87=reg59+reg87; reg59=reg14*reg3; reg66=reg105*reg66;
    reg100=reg115-reg100; reg75=reg42+reg75; reg41=reg113-reg41; reg13=reg50*reg13; reg10=reg15*reg10;
    reg42=reg108*reg77; reg111=reg112-reg111; reg112=reg61*reg5; reg108=reg108*reg33; reg120=reg120*reg77;
    reg96=reg93+reg96; reg124=reg83+reg124; reg43=reg44*reg43; reg20=reg50*reg20; reg116=reg94+reg116;
    reg83=reg12*reg3; reg93=reg50*reg41; reg114=reg103-reg114; reg94=reg100*reg71; reg103=reg109*reg52;
    reg95=reg95*reg53; reg113=reg111*reg52; reg66=reg75+reg66; reg59=reg91+reg59; reg75=reg110*reg24;
    reg91=reg22*reg58; reg115=reg15*reg111; reg118=reg50*reg100; reg99=reg69*reg99; reg7=reg86*reg7;
    reg97=reg19+reg97; reg19=reg129*reg122; reg42=reg116+reg42; reg116=reg41*reg71; reg108=reg124+reg108;
    reg119=reg129*reg43; reg81=reg131+reg81; reg124=reg85*reg62; reg131=reg11*reg3; T reg134=reg50*reg26;
    T reg135=reg49*reg16; reg34=reg32+reg34; reg32=reg44*reg25; T reg136=reg50*reg4; reg120=reg106+reg120;
    reg106=reg110*reg122; reg126=reg127+reg126; reg127=reg63*(*f.m).alpha_3; T reg137=(*f.m).alpha_1*reg60; reg101=reg123+reg101;
    reg123=(*f.m).alpha_2*reg61; T reg138=reg15*reg6; T reg139=reg62*(*f.m).alpha_3; reg96=reg112-reg96; reg112=reg133*reg33;
    T reg140=reg132*reg24; reg82=reg51+reg82; reg78=reg87+reg78; reg121=reg121*reg77; reg27=reg29+reg27;
    reg29=reg15*reg109; reg51=reg35*reg45; reg87=reg67*reg46; reg31=reg37+reg31; reg37=reg85*reg63;
    reg128=reg125+reg128; reg4=reg60*reg4; reg68=reg20-reg68; reg20=reg48*reg16; reg129=reg129*reg24;
    reg130=reg55+reg130; reg6=reg60*reg6; reg13=reg10-reg13; reg10=reg49*reg100; reg55=reg48*reg111;
    reg125=reg133*reg76; reg84=reg92+reg84; reg92=reg15*reg26; T reg141=reg48*reg109; reg110=reg110*reg43;
    reg90=reg18+reg90; reg18=reg132*reg43; reg17=reg98+reg17; reg98=reg49*reg41; T reg142=reg102*reg96;
    reg103=reg116+reg103; reg113=reg94+reg113; reg94=reg23*reg3; reg116=reg102*reg114; reg132=reg132*reg122;
    reg121=reg78+reg121; reg4=reg68-reg4; reg6=reg13-reg6; reg5=reg60*reg5; reg138=reg136+reg138;
    reg13=reg69*reg2; reg68=reg2*reg53; reg19=reg42+reg19; reg99=reg7+reg99; reg85=reg85*reg105;
    reg7=reg61*reg96; reg55=reg10+reg55; reg133=reg133*reg77; reg26=reg86*reg26; reg16=reg69*reg16;
    reg95=reg66+reg95; reg106=reg120+reg106; reg67=reg67*reg34; reg10=reg38*reg59; reg14=reg14*reg22;
    reg12=reg12*reg25; reg131=reg91+reg131; reg83=reg32+reg83; reg32=reg25*reg58; reg141=reg98+reg141;
    reg42=reg61*reg114; reg66=reg31*reg62; reg135=reg134+reg135; reg78=reg87*reg70; reg91=(*f.m).alpha_1*reg86;
    reg98=reg87*reg65; reg120=reg79*(*f.m).alpha_3; reg37=reg128+reg37; reg123=reg137+reg123; reg124=reg81+reg124;
    reg81=reg60*reg114; reg29=reg93+reg29; reg129=reg130+reg129; reg127=reg126+reg127; reg20=reg92+reg20;
    reg119=reg108+reg119; reg140=reg101+reg140; reg75=reg97+reg75; reg69=(*f.m).alpha_2*reg69; reg18=reg17+reg18;
    reg110=reg90+reg110; reg115=reg118+reg115; reg112=reg82+reg112; reg17=reg51*reg43; reg82=reg60*reg96;
    reg125=reg84+reg125; reg84=reg51*reg24; reg139=reg27+reg139; reg27=reg31*reg63; reg51=reg51*reg122;
    reg133=reg95+reg133; reg84=reg125+reg84; reg7=reg55+reg7; reg78=reg37+reg78; reg37=reg10*reg76;
    reg38=reg38*reg131; reg65=reg67*reg65; reg55=reg35*reg83; reg66=reg135+reg66; reg90=reg40*reg13;
    reg92=reg1*reg46; reg93=reg12*reg1; reg95=reg50*reg4; reg97=reg15*reg6; reg101=reg119*reg106;
    reg108=reg129*reg106; reg118=reg110*reg19; reg17=reg112+reg17; reg112=reg75*reg19; reg70=reg67*reg70;
    reg27=reg20+reg27; reg138=reg5-reg138; reg94=reg32+reg94; reg5=reg139*reg18; reg142=reg113+reg142;
    reg20=reg53*reg13; reg32=reg127*reg110; reg116=reg103+reg116; reg103=reg53*reg68; reg98=reg124+reg98;
    reg113=reg10*reg33; reg124=reg139*reg140; reg125=reg127*reg75; reg25=reg23*reg25; reg22=reg11*reg22;
    reg44=reg44*reg3; reg16=reg26+reg16; reg31=reg31*reg105; reg87=reg87*reg53; reg85=reg99+reg85;
    reg132=reg121+reg132; reg11=reg48*reg6; reg23=reg57*reg68; reg81=reg29+reg81; reg42=reg141+reg42;
    reg26=reg40*reg68; reg29=reg49*reg4; reg69=reg91+reg69; reg120=reg123+reg120; reg91=reg57*reg13;
    reg82=reg115+reg82; reg12=(*f.m).alpha_2*reg12; reg99=(*f.m).alpha_1*reg14; reg105=reg105*(*f.m).alpha_3; reg108=reg112-reg108;
    reg2=reg86*reg2; reg86=reg25*reg0; reg112=reg0*reg34; reg61=reg61*reg138; reg97=reg95+reg97;
    reg76=reg38*reg76; reg60=reg60*reg138; reg11=reg29+reg11; reg35=reg35*reg94; reg58=reg3*reg58;
    reg3=reg127*reg106; reg84=reg120*reg84; reg125=reg124+reg125; reg29=reg139*reg132; reg101=reg118-reg101;
    reg95=reg110*reg129; reg115=reg46*reg92; reg118=reg119*reg75; reg121=reg4*reg71; reg123=reg6*reg52;
    reg17=reg120*reg17; reg32=reg5+reg32; reg103=reg116+reg103; reg5=reg89*reg92; reg26=reg42+reg26;
    reg20=reg142+reg20; reg42=reg46*reg93; reg105=reg69+reg105; reg91=reg82+reg91; reg69=reg74*reg93;
    reg12=reg99+reg12; reg44=reg44*(*f.m).alpha_3; reg82=reg89*reg93; reg23=reg81+reg23; reg90=reg7+reg90;
    reg7=reg74*reg92; reg87=reg85+reg87; reg10=reg10*reg77; reg51=reg133+reg51; reg33=reg38*reg33;
    reg65=reg66+reg65; reg66=(*f.m).alpha_1*reg22; reg25=(*f.m).alpha_2*reg25; reg113=reg98+reg113; reg70=reg27+reg70;
    reg37=reg78+reg37; reg31=reg16+reg31; reg16=reg55*reg43; reg27=reg55*reg24; reg67=reg67*reg53;
    reg102=reg102*reg138; reg78=reg88*reg112; reg25=reg66+reg25; reg5=reg26+reg5; reg58=reg58*(*f.m).alpha_3;
    reg24=reg35*reg24; reg26=reg56*reg86; reg69=reg91+reg69; reg76=reg70+reg76; reg27=reg37+reg27;
    reg42=reg20+reg42; reg123=reg121+reg123; reg7=reg23+reg7; reg118=reg95-reg118; reg20=reg57*reg2;
    reg60=reg97+reg60; reg23=reg56*reg112; reg43=reg35*reg43; reg33=reg65+reg33; reg16=reg113+reg16;
    reg37=reg101*reg140; reg44=reg12+reg44; reg1=reg14*reg1; reg12=reg40*reg2; reg61=reg11+reg61;
    reg11=reg18*reg108; reg82=reg90+reg82; reg14=reg88*reg86; reg10=reg87+reg10; reg55=reg55*reg122;
    reg51=reg120*reg51; reg3=reg29+reg3; reg29=reg105*reg129; reg84=reg125+reg84; reg65=reg34*reg112;
    reg115=reg103+reg115; reg67=reg31+reg67; reg38=reg38*reg77; reg31=reg105*reg119; reg17=reg32+reg17;
    reg32=reg34*reg86; reg66=reg18*reg19; reg70=reg53*reg2; reg102=reg123+reg102; reg29=reg84+reg29;
    reg31=reg17+reg31; reg16=reg44*reg16; reg74=reg74*reg1; reg20=reg60+reg20; reg65=reg115+reg65;
    reg0=reg22*reg0; reg17=reg129*reg132; reg12=reg61+reg12; reg22=reg140*reg19; reg78=reg5+reg78;
    reg89=reg89*reg1; reg35=reg35*reg122; reg5=reg118*reg132; reg37=reg11-reg37; reg26=reg69+reg26;
    reg32=reg42+reg32; reg14=reg82+reg14; reg38=reg67+reg38; reg43=reg33+reg43; reg55=reg10+reg55;
    reg23=reg7+reg23; reg58=reg25+reg58; reg7=reg132*reg119; reg27=reg44*reg27; reg24=reg76+reg24;
    reg51=reg3+reg51; reg3=reg105*reg19; reg35=reg38+reg35; reg10=reg23*reg32; reg11=reg78*reg32;
    reg3=reg51+reg3; reg88=reg88*reg0; reg25=reg14*reg65; reg24=reg58*reg24; reg55=reg44*reg55;
    reg27=reg29+reg27; reg89=reg12+reg89; reg5=reg37+reg5; reg12=reg26*reg65; reg56=reg56*reg0;
    reg74=reg20+reg74; reg16=reg31+reg16; reg70=reg102+reg70; reg46=reg46*reg1; reg17=reg22-reg17;
    reg20=reg132*reg75; reg22=reg140*reg106; reg29=reg18*reg106; reg7=reg66-reg7; reg31=reg110*reg132;
    reg43=reg58*reg43; reg33=reg119*reg140; reg37=reg18*reg129; reg34=reg34*reg0; reg38=reg18*reg75;
    reg56=reg74+reg56; reg46=reg70+reg46; reg20=reg22-reg20; reg101=reg101/reg5; reg22=reg14*reg23;
    reg42=reg78*reg26; reg10=reg12-reg10; reg31=reg29-reg31; reg89=reg88+reg89; reg7=reg7/reg5;
    reg55=reg3+reg55; reg35=reg58*reg35; reg43=reg16+reg43; reg3=reg110*reg140; reg24=reg27+reg24;
    reg33=reg37-reg33; reg108=reg108/reg5; reg17=reg17/reg5; reg11=reg25-reg11; reg12=reg10*reg89;
    reg31=reg31/reg5; reg35=reg55+reg35; reg118=reg118/reg5; reg34=reg46+reg34; reg3=reg38-reg3;
    reg33=reg33/reg5; reg108=reg108*reg43; reg101=reg101*reg24; reg22=reg42-reg22; reg17=reg17*reg43;
    reg16=reg56*reg11; reg20=reg20/reg5; reg7=reg7*reg24; reg33=reg33*reg35; reg25=reg89*reg32;
    reg101=reg108-reg101; reg24=reg31*reg24; reg118=reg118*reg35; reg43=reg20*reg43; reg17=reg7-reg17;
    reg7=reg22*reg34; reg5=reg3/reg5; reg12=reg16-reg12; reg3=reg26*reg34; reg16=reg56*reg32;
    reg20=reg14*reg34; reg27=reg78*reg34; reg29=reg89*reg26; reg31=elem.pos(1)[0]-elem.pos(0)[0]; reg37=reg89*reg65;
    reg38=reg14*reg56; reg101=reg118+reg101; reg42=elem.pos(1)[1]-elem.pos(0)[1]; reg46=elem.pos(2)[0]-elem.pos(0)[0]; reg51=reg56*reg65;
    reg33=reg17-reg33; reg17=elem.pos(2)[1]-elem.pos(0)[1]; reg35=reg5*reg35; reg24=reg43-reg24; reg5=1-(*f.m).resolution;
    reg7=reg12+reg7; reg3=reg16-reg3; reg12=reg23*reg34; reg20=reg25-reg20; reg16=reg46*reg42;
    reg25=reg17*reg31; reg29=reg38-reg29; reg38=reg139*(*f.m).resolution; reg43=reg89*reg23; reg55=reg78*reg56;
    reg33=reg33*reg5; reg3=reg3/reg7; reg12=reg51-reg12; reg101=reg101*reg5; reg20=reg20/reg7;
    reg27=reg37-reg27; reg24=reg35+reg24; reg35=reg127*(*f.m).resolution; reg37=reg3*(*f.m).resolution; reg132=reg132*reg5;
    reg106=reg5*reg106; reg51=reg20*(*f.m).resolution; reg60=reg105*(*f.m).resolution; reg16=reg25-reg16; reg29=reg29/reg7;
    reg43=reg55-reg43; reg12=reg12/reg7; reg10=reg10/reg7; reg27=reg27/reg7; reg11=reg11/reg7;
    reg24=reg24*reg5; reg101=reg38+reg101; reg33=reg35+reg33; reg37=reg106-reg37; reg51=reg132+reg51;
    reg25=reg29*(*f.m).resolution; reg101=(*f.m).deltaT*reg101; reg33=(*f.m).deltaT*reg33; reg35=reg12*(*f.m).resolution; reg24=reg60+reg24;
    reg38=reg27*(*f.m).resolution; reg18=reg18*reg5; reg110=reg110*reg5; reg19=reg5*reg19; reg140=reg5*reg140;
    reg75=reg5*reg75; reg55=reg11*(*f.m).resolution; reg60=reg10*(*f.m).resolution; reg43=reg43/reg7; reg7=reg22/reg7;
    reg31=reg31/reg16; reg42=reg42/reg16; reg46=reg46/reg16; reg17=reg17/reg16; reg24=(*f.m).deltaT*reg24;
    reg119=reg5*reg119; reg25=reg19+reg25; reg19=reg43*(*f.m).resolution; reg75=reg35+reg75; reg38=reg140-reg38;
    reg60=reg110-reg60; reg18=reg55+reg18; reg5=reg129*reg5; reg22=reg46-reg31; reg35=reg42-reg17;
    reg55=reg7*(*f.m).resolution; reg61=reg37*reg33; reg66=reg51*reg101; reg67=0.5*reg17; reg69=reg18*reg101;
    reg70=0.5*reg46; reg74=0.5*reg22; reg119=reg55+reg119; reg19=reg5-reg19; reg5=reg38*reg101;
    reg55=0.5*reg42; reg76=0.5*reg35; reg81=reg66+reg61; reg82=reg25*reg24; reg84=reg60*reg33;
    reg85=0.5*reg31; reg87=reg75*reg33; reg88=reg25*reg70; reg90=reg51*reg17; reg91=reg25*reg55;
    reg95=reg37*reg31; reg97=reg25*reg76; reg98=reg37*reg46; reg99=reg22*reg37; reg102=reg74*reg25;
    reg103=reg35*reg51; reg106=reg67*reg25; reg108=reg25*reg85; reg110=reg81+reg82; reg113=reg51*reg42;
    reg115=reg119*reg24; reg116=reg87+reg5; reg118=reg19*reg24; reg121=reg84+reg69; reg108=reg108-reg113;
    reg90=reg90-reg88; reg123=reg60*reg31; reg124=reg116+reg118; reg125=reg60*reg46; reg126=reg67*reg119;
    reg128=reg119*reg85; reg129=reg18*reg42; reg106=reg106-reg98; reg130=2*reg110; reg132=reg119*reg55;
    reg133=reg119*reg70; reg134=reg18*reg17; reg135=reg121+reg115; reg97=reg99+reg97; reg99=reg35*reg18;
    reg136=reg74*reg119; reg137=reg119*reg76; reg140=reg22*reg60; reg102=reg103+reg102; reg103=reg22*reg75;
    reg141=reg19*reg76; reg142=reg67*reg19; T reg143=reg75*reg46; T reg144=reg19*reg70; T reg145=reg38*reg17;
    T reg146=reg75*reg31; reg95=reg95-reg91; T reg147=reg19*reg55; T reg148=reg19*reg85; T reg149=reg38*reg42;
    T reg150=reg130*reg67; reg136=reg99+reg136; reg108=2*reg108; reg142=reg142-reg143; reg106=2*reg106;
    reg95=2*reg95; reg134=reg134-reg133; reg146=reg146-reg147; reg128=reg128-reg129; reg145=reg145-reg144;
    reg126=reg126-reg125; reg141=reg103+reg141; reg90=2*reg90; reg99=reg130*reg70; reg103=reg124*reg46;
    T reg151=reg135*reg17; reg137=reg140+reg137; reg97=2*reg97; reg148=reg148-reg149; reg140=reg130*reg85;
    T reg152=reg135*reg42; reg102=2*reg102; reg123=reg123-reg132; T reg153=reg124*reg31; T reg154=reg130*reg55;
    T reg155=reg106*reg70; T reg156=reg76*reg90; T reg157=reg102*reg74; T reg158=reg128*reg17; T reg159=reg108*reg70;
    T reg160=reg145*reg22; T reg161=reg106*reg76; T reg162=reg70*reg95; T reg163=reg123*reg17; T reg164=reg142*reg22;
    T reg165=reg35*reg136; T reg166=reg76*reg95; T reg167=reg22*reg148; T reg168=reg130*reg74; T reg169=reg130*reg76;
    T reg170=reg108*reg76; T reg171=reg146*reg31; T reg172=reg124*reg22; T reg173=reg99-reg151; T reg174=reg74*reg108;
    T reg175=reg55*reg95; T reg176=reg134*reg17; T reg177=reg67*reg106; T reg178=reg142*reg46; T reg179=reg35*reg135;
    T reg180=reg90*reg70; T reg181=reg67*reg108; T reg182=reg148*reg46; T reg183=reg67*reg95; T reg184=reg146*reg46;
    T reg185=reg128*reg42; T reg186=reg85*reg95; T reg187=reg123*reg42; T reg188=reg126*reg17; T reg189=reg97*reg76;
    reg128=reg35*reg128; T reg190=reg74*reg106; T reg191=reg35*reg126; T reg192=reg85*reg108; T reg193=reg22*reg141;
    reg146=reg22*reg146; T reg194=reg103-reg150; T reg195=reg74*reg90; T reg196=reg35*reg134; reg95=reg74*reg95;
    T reg197=reg152-reg140; T reg198=reg35*reg137; reg123=reg35*reg123; T reg199=reg154-reg153; T reg200=reg97*reg74;
    reg167=reg170+reg167; reg195=reg196+reg195; reg146=reg166+reg146; reg173=reg173*reg16; reg187=reg186-reg187;
    reg193=reg189+reg193; reg174=reg128+reg174; reg185=reg192-reg185; reg175=reg171-reg175; reg180=reg176-reg180;
    reg162=reg163-reg162; reg178=reg177-reg178; reg182=reg181-reg182; reg191=reg190+reg191; reg128=reg179+reg168;
    reg200=reg198+reg200; reg164=reg161+reg164; reg197=reg197*reg16; reg159=reg158-reg159; reg184=reg183-reg184;
    reg194=reg194*reg16; reg160=reg156+reg160; reg95=reg123+reg95; reg123=reg169+reg172; reg199=reg199*reg16;
    reg157=reg165+reg157; reg155=reg188-reg155; reg146=reg146*reg16; reg174=reg174*reg16; reg194=ponderation*reg194;
    reg173=ponderation*reg173; reg167=reg167*reg16; reg175=reg175*reg16; reg156=reg123*reg16; reg197=ponderation*reg197;
    reg199=ponderation*reg199; reg178=reg178*reg16; reg158=reg128*reg16; reg162=reg162*reg16; reg157=reg157*reg16;
    reg184=reg184*reg16; reg164=reg164*reg16; reg200=reg200*reg16; reg159=reg159*reg16; reg160=reg160*reg16;
    reg95=reg95*reg16; reg155=reg155*reg16; reg195=reg195*reg16; reg187=reg187*reg16; reg185=reg185*reg16;
    reg191=reg191*reg16; reg193=reg193*reg16; reg182=reg182*reg16; reg180=reg180*reg16; T tmp_1_4=ponderation*reg167;
    T tmp_0_0=ponderation*reg157; T tmp_0_3=ponderation*reg191; T tmp_2_5=ponderation*reg162; T tmp_2_4=ponderation*reg159; T vec_4=-reg197;
    T tmp_3_3=ponderation*reg178; reg157=ponderation*reg158; T vec_0=reg157; T tmp_3_4=ponderation*reg182; T tmp_0_1=ponderation*reg200;
    T tmp_1_1=ponderation*reg193; T vec_3=-reg194; T tmp_1_3=ponderation*reg164; T tmp_1_2=ponderation*reg160; T tmp_5_5=ponderation*reg175;
    T tmp_0_5=ponderation*reg95; reg95=ponderation*reg156; T vec_1=reg95; T tmp_2_3=ponderation*reg155; T tmp_3_5=ponderation*reg184;
    T vec_5=-reg199; T tmp_2_2=ponderation*reg180; T tmp_4_5=ponderation*reg187; T vec_2=-reg173; T tmp_0_2=ponderation*reg195;
    T tmp_0_4=ponderation*reg174; T tmp_4_4=ponderation*reg185; T tmp_1_5=ponderation*reg146;
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
    T reg0=pow((*f.m).v1[1],2); T reg1=pow((*f.m).v1[0],2); T reg2=pow((*f.m).v2[1],2); T reg3=pow((*f.m).v2[0],2); T reg4=pow((*f.m).v1[2],2);
    reg0=reg1+reg0; reg1=2*(*f.m).shear_modulus_23; T reg5=pow((*f.m).v2[2],2); reg2=reg3+reg2; reg3=2*(*f.m).shear_modulus_13;
    reg4=reg0+reg4; reg1=1.0/reg1; reg5=reg2+reg5; reg3=1.0/reg3; reg4=pow(reg4,0.5);
    reg0=2*(*f.m).shear_modulus_12; reg2=reg3*reg1; reg0=1.0/reg0; T reg6=(*f.m).v1[0]/reg4; T reg7=(*f.m).v1[1]/reg4;
    reg5=pow(reg5,0.5); T reg8=reg0*reg2; T reg9=2*reg7; T reg10=2*reg6; reg4=(*f.m).v1[2]/reg4;
    T reg11=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg12=1.0/(*f.m).elastic_modulus_3; T reg13=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg14=(*f.m).v2[0]/reg5; T reg15=(*f.m).v2[1]/reg5;
    T reg16=pow(reg14,2); T reg17=reg13*reg8; T reg18=reg9*reg15; T reg19=reg10*reg14; T reg20=reg12*reg8;
    T reg21=reg11*reg8; reg5=(*f.m).v2[2]/reg5; T reg22=pow(reg15,2); T reg23=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg24=1.0/(*f.m).elastic_modulus_2;
    T reg25=1.0/(*f.m).elastic_modulus_1; T reg26=2*reg4; T reg27=reg23*reg20; T reg28=reg11*reg21; T reg29=reg11*reg17;
    T reg30=reg16*reg25; T reg31=reg22*reg23; T reg32=reg25*reg19; T reg33=reg18*reg23; T reg34=reg16*reg23;
    T reg35=reg22*reg24; T reg36=reg24*reg20; T reg37=pow(reg5,2); reg26=reg26*reg5; T reg38=reg23*reg19;
    T reg39=reg18*reg24; T reg40=reg22*reg11; T reg41=reg16*reg13; T reg42=pow(reg6,2); T reg43=pow(reg7,2);
    reg28=reg36-reg28; reg36=reg13*reg19; T reg44=reg18*reg11; reg29=reg29+reg27; T reg45=reg23*reg21;
    T reg46=reg24*reg17; reg34=reg35-reg34; reg35=reg37*reg13; reg31=reg30-reg31; reg30=reg37*reg11;
    T reg47=reg26*reg11; reg38=reg39-reg38; reg39=reg26*reg13; reg33=reg32-reg33; reg32=reg6*reg15;
    T reg48=reg7*reg14; T reg49=reg6*reg14; reg30=reg34-reg30; reg34=reg7*reg15; T reg50=reg43*reg24;
    T reg51=reg37*reg12; T reg52=reg42*reg23; reg39=reg33-reg39; reg40=reg41+reg40; reg26=reg26*reg12;
    reg44=reg36+reg44; reg33=reg29*reg23; reg35=reg31-reg35; reg31=reg45+reg46; reg36=reg43*reg23;
    reg41=reg25*reg28; reg47=reg38-reg47; reg38=reg42*reg25; T reg53=pow(reg4,2); T reg54=reg42*reg13;
    T reg55=reg14*reg15; T reg56=reg43*reg11; T reg57=reg13*reg21; reg40=reg51-reg40; reg51=reg24*reg8;
    T reg58=reg16*reg35; T reg59=reg22*reg30; T reg60=reg31*reg13; T reg61=2*reg14; T reg62=reg16*reg39;
    reg33=reg41-reg33; reg41=reg35*reg49; T reg63=reg30*reg34; T reg64=reg12*reg2; T reg65=reg22*reg47;
    reg52=reg50-reg52; reg50=reg53*reg11; T reg66=reg53*reg13; reg36=reg38-reg36; reg44=reg26-reg44;
    reg26=reg39*reg49; reg38=reg47*reg34; T reg67=reg42*reg39; T reg68=reg43*reg47; T reg69=reg4*reg14;
    T reg70=reg6*reg5; T reg71=reg4*reg15; T reg72=reg7*reg5; T reg73=reg0*reg1; T reg74=reg4*reg5;
    T reg75=reg43*reg30; T reg76=reg42*reg35; reg8=reg23*reg8; T reg77=reg13*reg17; T reg78=reg13*reg2;
    T reg79=reg48+reg32; reg20=reg25*reg20; reg2=reg11*reg2; T reg80=reg73*reg11; T reg81=reg0*reg3;
    T reg82=reg37*reg44; reg17=reg23*reg17; reg65=reg62+reg65; reg62=reg74*reg40; reg63=reg41+reg63;
    reg41=reg8*reg13; T reg83=reg73*reg13; reg21=reg25*reg21; T reg84=reg37*reg40; reg77=reg20-reg77;
    reg59=reg58+reg59; reg20=reg51*reg13; reg58=reg14*reg5; reg27=reg57+reg27; reg57=reg64*reg24;
    T reg85=reg53*reg44; reg68=reg67+reg68; reg67=reg53*reg40; reg75=reg76+reg75; reg76=reg70+reg69;
    T reg86=reg53*reg12; reg56=reg54+reg56; reg60=reg33-reg60; reg33=reg72-reg71; reg73=reg73*reg12;
    reg54=reg10*reg7; T reg87=reg78*reg11; T reg88=reg61*reg15; T reg89=reg0*reg79; T reg90=reg55*reg0;
    T reg91=reg2*reg11; reg66=reg36-reg66; reg64=reg64*reg23; reg50=reg52-reg50; reg38=reg26+reg38;
    reg26=reg74*reg44; reg36=reg83*reg11; reg12=reg81*reg12; reg67=reg75+reg67; reg78=reg78*reg24;
    reg52=reg54*reg90; reg75=reg81*reg11; reg85=reg68+reg85; reg68=reg54*reg89; reg2=reg2*reg23;
    T reg92=reg22*reg50; T reg93=reg16*reg66; T reg94=reg43*reg50; reg81=reg81*reg13; T reg95=reg80*reg11;
    T reg96=reg73*reg23; reg73=reg73*reg24; T reg97=reg6*reg7; T reg98=reg88*reg89; reg82=reg65+reg82;
    reg65=reg58*reg3; T reg99=reg3*reg76; reg87=reg64+reg87; reg84=reg59+reg84; reg59=reg88*reg90;
    reg91=reg57-reg91; reg57=reg42*reg66; reg64=reg61*reg5; T reg100=2*reg15; reg70=reg69-reg70;
    reg69=2*reg33; reg72=reg71+reg72; reg56=reg86-reg56; reg28=reg28/reg60; reg71=reg10*reg4;
    reg51=reg51*reg25; reg86=reg79*reg89; reg26=reg38+reg26; reg38=reg79*reg90; reg62=reg63+reg62;
    reg63=reg15*reg5; reg27=reg27/reg60; reg20=reg45+reg20; reg8=reg8*reg23; reg29=reg29/reg60;
    reg77=reg77/reg60; reg41=reg21+reg41; reg17=reg21+reg17; reg21=reg81*reg11; T reg101=reg53*reg56;
    reg11=reg75*reg11; T reg102=reg12*reg23; reg12=reg12*reg24; T reg103=reg22*reg77; reg94=reg57+reg94;
    reg36=reg96+reg36; reg95=reg73-reg95; reg17=reg17/reg60; reg41=reg41/reg60; reg57=reg2+reg78;
    reg73=reg16*reg29; reg87=reg87*reg23; reg96=reg9*reg4; T reg104=reg42*reg28; reg91=reg91*reg25;
    T reg105=reg1*reg72; T reg106=reg63*reg1; reg92=reg93+reg92; reg93=reg37*reg56; T reg107=reg71*reg99;
    reg68=reg85+reg68; reg85=reg16*reg77; T reg108=reg42*reg27; T reg109=reg100*reg5; reg8=reg51-reg8;
    reg51=reg71*reg65; reg52=reg67+reg52; reg67=pow(reg33,2); T reg110=pow(reg70,2); T reg111=reg69*reg70;
    T reg112=reg43*reg27; T reg113=reg66*reg49; T reg114=reg50*reg34; T reg115=reg88*reg29; T reg116=reg54*reg28;
    T reg117=reg22*reg29; T reg118=reg43*reg28; T reg119=reg64*reg99; reg98=reg82+reg98; reg80=reg80*reg23;
    reg86=reg26+reg86; reg26=reg76*reg65; reg38=reg62+reg38; reg83=reg83*reg24; reg62=reg54*reg27;
    reg0=reg97*reg0; reg82=reg76*reg99; T reg120=reg64*reg65; reg59=reg84+reg59; reg84=reg88*reg77;
    T reg121=reg6*reg4; reg31=reg31/reg60; reg20=reg20/reg60; reg11=reg12-reg11; reg101=reg94+reg101;
    reg12=reg109*reg105; reg21=reg102+reg21; reg94=reg54*reg0; reg119=reg98+reg119; reg117=reg118+reg117;
    reg98=reg31*reg110; reg102=reg17*reg67; reg26=reg38+reg26; reg115=reg116+reg115; reg38=reg31*reg111;
    reg116=reg22*reg41; reg118=reg74*reg56; reg75=reg75*reg23; reg93=reg92+reg93; reg114=reg113+reg114;
    reg92=reg88*reg0; reg113=reg96*reg105; reg107=reg68+reg107; reg85=reg108+reg85; reg68=reg109*reg106;
    reg120=reg59+reg120; reg59=reg17*reg111; reg8=reg8/reg60; reg51=reg52+reg51; reg52=reg96*reg106;
    reg84=reg62+reg84; reg62=reg72*reg105; reg108=reg31*reg67; reg3=reg121*reg3; T reg122=reg88*reg41;
    reg73=reg104+reg73; reg104=reg17*reg110; reg24=reg81*reg24; reg103=reg112+reg103; reg36=reg36*reg23;
    reg81=reg43*reg20; reg87=reg91-reg87; reg95=reg95*reg25; reg91=reg80+reg83; reg112=reg7*reg4;
    reg82=reg86+reg82; reg86=reg16*reg41; T reg123=reg54*reg20; reg57=reg57*reg13; T reg124=reg72*reg106;
    T reg125=reg42*reg20; reg68=reg120+reg68; reg102=reg85+reg102; reg85=reg8*reg110; reg108=reg73+reg108;
    reg52=reg51+reg52; reg62=reg82+reg62; reg51=reg8*reg67; reg86=reg125+reg86; reg113=reg107+reg113;
    reg73=reg64*reg3; reg92=reg93+reg92; reg59=reg84+reg59; reg1=reg112*reg1; reg91=reg91*reg13;
    reg94=reg101+reg94; reg122=reg123+reg122; reg82=reg75+reg24; reg98=reg117+reg98; reg84=reg8*reg111;
    reg36=reg95-reg36; reg38=reg115+reg38; reg93=reg71*reg3; reg124=reg26+reg124; reg116=reg81+reg116;
    reg23=reg21*reg23; reg21=reg79*reg0; reg118=reg114+reg118; reg104=reg103+reg104; reg57=reg87-reg57;
    reg12=reg119+reg12; reg25=reg11*reg25; reg84=reg122+reg84; reg11=reg33*reg70; reg26=reg68*reg62;
    reg81=reg6*reg70; reg87=reg33*reg7; reg95=reg52*reg62; reg101=reg12*reg124; reg93=reg94+reg93;
    reg94=reg55*reg104; reg103=reg98*reg97; reg107=reg96*reg1; reg85=reg116+reg85; reg73=reg92+reg73;
    reg92=reg97*reg108; reg114=reg55*reg102; reg21=reg118+reg21; reg115=reg76*reg3; reg13=reg82*reg13;
    reg23=reg25-reg23; reg91=reg36-reg91; reg57=reg57/reg60; reg51=reg86+reg51; reg25=reg109*reg1;
    reg36=reg113*reg124; reg82=reg38*reg97; reg86=reg55*reg59; reg94=reg103+reg94; reg103=reg11*reg85;
    reg101=reg26-reg101; reg115=reg21+reg115; reg13=reg23-reg13; reg21=reg33*reg15; reg81=reg87+reg81;
    reg36=reg95-reg36; reg91=reg91/reg60; reg23=reg12*reg52; reg26=reg57*reg49; reg87=reg57*reg34;
    reg95=reg57*reg79; reg116=reg42*reg108; reg117=reg16*reg102; reg118=reg7*reg70; reg73=reg25+reg73;
    reg25=reg14*reg70; reg119=reg72*reg1; reg120=reg68*reg113; reg122=reg11*reg84; reg86=reg82+reg86;
    reg82=reg22*reg104; reg123=reg43*reg98; reg102=reg22*reg102; reg108=reg43*reg108; reg125=reg22*reg59;
    T reg126=reg43*reg38; reg107=reg93+reg107; reg93=reg33*reg6; reg114=reg92+reg114; reg98=reg42*reg98;
    reg38=reg42*reg38; reg59=reg16*reg59; reg104=reg16*reg104; reg92=reg11*reg51; T reg127=elem.pos(2)[0]-elem.pos(0)[0];
    reg59=reg38+reg59; reg119=reg115+reg119; reg125=reg126+reg125; reg38=reg36*reg73; reg115=reg110*reg84;
    reg92=reg114+reg92; reg103=reg94+reg103; reg120=reg23-reg120; reg23=reg87*reg79; reg25=reg21+reg25;
    reg21=reg107*reg101; reg94=elem.pos(1)[0]-elem.pos(0)[0]; reg102=reg108+reg102; reg108=reg110*reg51; reg114=reg33*reg14;
    reg126=reg15*reg70; T reg128=elem.pos(1)[1]-elem.pos(0)[1]; reg84=reg67*reg84; reg60=reg13/reg60; reg13=elem.pos(2)[1]-elem.pos(0)[1];
    T reg129=reg91*reg93; T reg130=reg91*reg118; T reg131=reg91*reg81; T reg132=reg26*reg79; T reg133=reg67*reg85;
    reg104=reg98+reg104; reg117=reg116+reg117; reg98=reg95*reg79; reg122=reg86+reg122; reg51=reg67*reg51;
    reg85=reg110*reg85; reg82=reg123+reg82; reg98=reg122+reg98; reg86=reg87*reg19; reg116=reg131*reg81;
    reg133=reg104+reg133; reg10=reg33*reg10; reg115=reg125+reg115; reg104=reg73*reg52; reg122=reg68*reg107;
    reg123=reg95*reg18; reg125=reg127*reg128; reg132=reg92+reg132; reg92=reg129*reg81; T reg134=reg13*reg94;
    T reg135=reg52*reg119; T reg136=reg107*reg124; reg95=reg95*reg19; T reg137=reg60*reg114; T reg138=reg60*reg126;
    reg84=reg59+reg84; reg108=reg102+reg108; reg59=reg60*reg25; reg102=reg68*reg119; reg9=reg9*reg70;
    T reg139=reg73*reg124; reg51=reg117+reg51; reg117=reg26*reg19; reg26=reg26*reg18; T reg140=reg120*reg119;
    reg38=reg21-reg38; reg23=reg103+reg23; reg21=reg130*reg81; reg87=reg87*reg18; reg85=reg82+reg85;
    reg95=reg84+reg95; reg82=reg131*reg10; reg26=reg108+reg26; reg84=reg129*reg9; reg87=reg85+reg87;
    reg85=reg130*reg9; reg61=reg33*reg61; reg103=reg59*reg25; reg116=reg98+reg116; reg125=reg134-reg125;
    reg98=reg138*reg25; reg21=reg23+reg21; reg100=reg100*reg70; reg131=reg131*reg9; reg123=reg115+reg123;
    reg23=reg137*reg25; reg92=reg132+reg92; reg130=reg130*reg10; reg86=reg133+reg86; reg129=reg129*reg10;
    reg117=reg51+reg117; reg51=reg73*reg113; reg104=reg122-reg104; reg108=reg12*reg107; reg135=reg136-reg135;
    reg115=reg113*reg119; reg102=reg139-reg102; reg140=reg38+reg140; reg38=reg73*reg62; reg122=reg12*reg119;
    reg132=reg107*reg62; reg82=reg95+reg82; reg51=reg108-reg51; reg129=reg117+reg129; reg95=reg137*reg61;
    reg85=reg87+reg85; reg87=reg138*reg100; reg108=reg59*reg100; reg131=reg123+reg131; reg130=reg86+reg130;
    reg138=reg138*reg61; reg13=reg13/reg125; reg23=reg92+reg23; reg127=reg127/reg125; reg98=reg21+reg98;
    reg103=reg116+reg103; reg115=reg132-reg115; reg135=reg135/reg140; reg137=reg137*reg100; reg84=reg26+reg84;
    reg102=reg102/reg140; reg122=reg38-reg122; reg104=reg104/reg140; reg21=1-(*f.m).resolution; reg94=reg94/reg125;
    reg128=reg128/reg125; reg59=reg59*reg61; reg26=reg102*(*f.m).resolution; reg38=reg21*reg103; reg86=reg104*(*f.m).resolution;
    reg92=reg135*(*f.m).resolution; reg116=reg21*reg98; reg117=reg23*reg21; reg123=reg127-reg94; reg108=reg131+reg108;
    reg59=reg82+reg59; reg137=reg84+reg137; reg87=reg85+reg87; reg51=reg51/reg140; reg120=reg120/reg140;
    reg115=reg115/reg140; reg36=reg36/reg140; reg82=reg128-reg13; reg95=reg129+reg95; reg101=reg101/reg140;
    reg138=reg130+reg138; reg140=reg122/reg140; reg84=reg21*reg87; reg85=0.5*reg128; reg122=0.5*reg94;
    reg129=0.5*reg123; reg130=0.5*reg13; reg131=reg120*(*f.m).resolution; reg132=reg101*(*f.m).resolution; reg133=reg36*(*f.m).resolution;
    reg26=reg117+reg26; reg92=reg116-reg92; reg86=reg38+reg86; reg38=reg95*reg21; reg116=reg51*(*f.m).resolution;
    reg117=reg115*(*f.m).resolution; reg134=reg140*(*f.m).resolution; reg136=0.5*reg127; reg139=0.5*reg82; T reg141=reg138*reg21;
    T reg142=reg21*reg59; T reg143=reg21*reg137; T reg144=reg108*reg21; T reg145=reg92*reg94; T reg146=reg86*reg85;
    T reg147=reg129*reg86; T reg148=reg82*reg26; T reg149=reg130*reg86; reg116=reg144-reg116; reg84=reg117+reg84;
    reg134=reg143-reg134; reg142=reg131+reg142; reg133=reg141-reg133; reg38=reg132+reg38; reg117=reg26*reg13;
    reg131=reg86*reg136; reg132=reg86*reg122; reg141=reg26*reg128; reg143=reg92*reg127; reg144=reg86*reg139;
    T reg150=reg123*reg92; reg117=reg117-reg131; T reg151=reg134*reg13; T reg152=reg133*reg127; T reg153=reg142*reg136;
    T reg154=reg38*reg13; T reg155=reg130*reg142; T reg156=reg116*reg136; T reg157=reg116*reg139; T reg158=reg123*reg84;
    T reg159=reg129*reg142; T reg160=reg82*reg38; T reg161=reg84*reg127; T reg162=reg134*reg128; T reg163=reg116*reg122;
    T reg164=reg142*reg85; T reg165=reg133*reg94; reg132=reg132-reg141; T reg166=reg142*reg122; T reg167=reg38*reg128;
    reg144=reg150+reg144; reg150=reg142*reg139; T reg168=reg123*reg133; reg149=reg149-reg143; T reg169=reg84*reg94;
    T reg170=reg116*reg85; reg147=reg148+reg147; reg148=reg129*reg116; T reg171=reg82*reg134; reg145=reg145-reg146;
    T reg172=reg130*reg116; reg151=reg151-reg156; reg157=reg158+reg157; reg155=reg155-reg152; reg148=reg171+reg148;
    reg117=2*reg117; reg145=2*reg145; reg165=reg165-reg164; reg154=reg154-reg153; reg172=reg172-reg161;
    reg147=2*reg147; reg159=reg160+reg159; reg163=reg163-reg162; reg149=2*reg149; reg150=reg168+reg150;
    reg169=reg169-reg170; reg132=2*reg132; reg144=2*reg144; reg166=reg166-reg167; reg158=reg155*reg13;
    reg160=reg129*reg149; reg168=reg144*reg129; reg171=reg117*reg136; T reg173=reg154*reg13; T reg174=reg144*reg122;
    T reg175=reg157*reg94; T reg176=reg139*reg145; T reg177=reg151*reg94; T reg178=reg85*reg117; T reg179=reg122*reg149;
    T reg180=reg123*reg169; T reg181=reg147*reg85; T reg182=reg172*reg94; T reg183=reg172*reg123; T reg184=reg132*reg139;
    T reg185=reg82*reg165; T reg186=reg129*reg132; T reg187=reg130*reg144; T reg188=reg123*reg157; T reg189=reg148*reg127;
    T reg190=reg130*reg147; T reg191=reg82*reg166; T reg192=reg129*reg145; T reg193=reg136*reg145; T reg194=reg165*reg13;
    T reg195=reg144*reg136; T reg196=reg155*reg128; T reg197=reg132*reg136; T reg198=reg166*reg13; T reg199=reg82*reg159;
    T reg200=reg149*reg136; T reg201=reg149*reg139; reg172=reg172*reg127; T reg202=reg147*reg139; T reg203=reg122*reg117;
    reg155=reg82*reg155; T reg204=reg147*reg122; T reg205=reg130*reg132; T reg206=reg163*reg127; T reg207=reg85*reg145;
    T reg208=reg169*reg94; T reg209=reg123*reg148; T reg210=reg129*reg117; reg166=reg166*reg128; T reg211=reg150*reg128;
    T reg212=reg159*reg128; T reg213=reg130*reg145; T reg214=reg82*reg154; reg169=reg169*reg127; T reg215=reg123*reg163;
    T reg216=reg144*reg85; T reg217=reg85*reg149; reg148=reg148*reg94; T reg218=reg82*reg150; T reg219=reg122*reg132;
    reg163=reg163*reg94; reg132=reg85*reg132; reg157=reg157*reg127; reg159=reg159*reg13; reg144=reg144*reg139;
    T reg220=reg147*reg136; T reg221=reg130*reg117; T reg222=reg151*reg127; reg147=reg147*reg129; reg117=reg139*reg117;
    reg154=reg154*reg128; reg149=reg130*reg149; reg151=reg151*reg123; reg165=reg165*reg128; reg145=reg122*reg145;
    reg150=reg150*reg13; reg209=reg202+reg209; reg188=reg144+reg188; reg168=reg218+reg168; reg181=reg148-reg181;
    reg165=reg145-reg165; reg147=reg199+reg147; reg166=reg219-reg166; reg212=reg204-reg212; reg196=reg179-reg196;
    reg169=reg213-reg169; reg157=reg187-reg157; reg206=reg205-reg206; reg210=reg214+reg210; reg207=reg208-reg207;
    reg172=reg149-reg172; reg195=reg150-reg195; reg151=reg117+reg151; reg220=reg159-reg220; reg222=reg221-reg222;
    reg132=reg163-reg132; reg154=reg203-reg154; reg217=reg182-reg217; reg180=reg176+reg180; reg178=reg177-reg178;
    reg216=reg175-reg216; reg171=reg173-reg171; reg200=reg158-reg200; reg155=reg160+reg155; reg197=reg198-reg197;
    reg193=reg194-reg193; reg211=reg174-reg211; reg186=reg191+reg186; reg189=reg190-reg189; reg192=reg185+reg192;
    reg215=reg184+reg215; reg183=reg201+reg183; reg212=reg212*reg125; reg147=reg147*reg125; reg168=reg168*reg125;
    reg186=reg186*reg125; reg155=reg155*reg125; reg215=reg215*reg125; reg180=reg180*reg125; reg195=reg195*reg125;
    reg207=reg207*reg125; reg132=reg132*reg125; reg220=reg220*reg125; reg209=reg209*reg125; reg183=reg183*reg125;
    reg192=reg192*reg125; reg189=reg189*reg125; reg193=reg193*reg125; reg197=reg197*reg125; reg200=reg200*reg125;
    reg171=reg171*reg125; reg216=reg216*reg125; reg178=reg178*reg125; reg217=reg217*reg125; reg157=reg157*reg125;
    reg222=reg222*reg125; reg151=reg151*reg125; reg172=reg172*reg125; reg210=reg210*reg125; reg206=reg206*reg125;
    reg169=reg169*reg125; reg211=reg211*reg125; reg166=reg166*reg125; reg188=reg188*reg125; reg154=reg154*reg125;
    reg196=reg196*reg125; reg181=reg181*reg125; reg165=reg165*reg125; T tmp_5_2=ponderation*reg178; T tmp_1_5=ponderation*reg180;
    T tmp_5_3=ponderation*reg217; T tmp_4_3=ponderation*reg196; T tmp_3_1=ponderation*reg157; T tmp_3_2=ponderation*reg222; T tmp_4_2=ponderation*reg154;
    T tmp_1_2=ponderation*reg151; T tmp_2_0=ponderation*reg220; T tmp_5_5=ponderation*reg207; T tmp_4_0=ponderation*reg212; T tmp_3_3=ponderation*reg172;
    T tmp_0_2=ponderation*reg210; T tmp_5_4=ponderation*reg132; T tmp_4_1=ponderation*reg211; T tmp_3_4=ponderation*reg206; T tmp_1_3=ponderation*reg183;
    T tmp_1_0=ponderation*reg209; T tmp_0_1=ponderation*reg168; T tmp_0_5=ponderation*reg192; T tmp_0_4=ponderation*reg186; T tmp_1_1=ponderation*reg188;
    T tmp_5_0=ponderation*reg181; T tmp_3_0=ponderation*reg189; T tmp_2_5=ponderation*reg193; T tmp_2_4=ponderation*reg197; T tmp_4_5=ponderation*reg165;
    T tmp_0_3=ponderation*reg155; T tmp_2_3=ponderation*reg200; T tmp_2_2=ponderation*reg171; T tmp_0_0=ponderation*reg147; T tmp_2_1=ponderation*reg195;
    T tmp_1_4=ponderation*reg215; T tmp_4_4=ponderation*reg166; T tmp_5_1=ponderation*reg216; T tmp_3_5=ponderation*reg169;
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
    reg0=reg1+reg0; reg4=reg0+reg4; reg0=2*(*f.m).shear_modulus_13; reg1=2*(*f.m).shear_modulus_23; T reg5=pow((*f.m).v2[2],2);
    reg2=reg3+reg2; reg4=pow(reg4,0.5); reg3=2*(*f.m).shear_modulus_12; reg0=1.0/reg0; reg1=1.0/reg1;
    reg5=reg2+reg5; reg5=pow(reg5,0.5); reg2=(*f.m).v1[1]/reg4; T reg6=(*f.m).v1[0]/reg4; T reg7=reg0*reg1;
    reg3=1.0/reg3; T reg8=2*reg6; T reg9=2*reg2; T reg10=reg3*reg7; T reg11=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg12=1.0/(*f.m).elastic_modulus_3; reg4=(*f.m).v1[2]/reg4; T reg13=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg14=(*f.m).v2[1]/reg5; T reg15=(*f.m).v2[0]/reg5;
    T reg16=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg17=1.0/(*f.m).elastic_modulus_2; T reg18=1.0/(*f.m).elastic_modulus_1; T reg19=reg9*reg14; T reg20=2*reg4;
    T reg21=reg12*reg10; T reg22=reg11*reg10; reg5=(*f.m).v2[2]/reg5; T reg23=pow(reg15,2); T reg24=pow(reg14,2);
    T reg25=reg13*reg10; T reg26=reg8*reg15; T reg27=reg18*reg26; T reg28=reg19*reg16; T reg29=reg24*reg16;
    T reg30=reg23*reg18; T reg31=reg23*reg16; T reg32=reg24*reg17; T reg33=reg16*reg26; T reg34=reg19*reg17;
    T reg35=pow(reg5,2); T reg36=reg16*reg21; T reg37=reg11*reg22; T reg38=reg17*reg21; T reg39=reg11*reg25;
    reg20=reg20*reg5; reg28=reg27-reg28; reg27=reg20*reg13; T reg40=reg24*reg11; T reg41=reg23*reg13;
    T reg42=reg35*reg11; T reg43=reg17*reg25; reg31=reg32-reg31; reg29=reg30-reg29; reg30=reg35*reg13;
    reg32=reg20*reg11; T reg44=reg19*reg11; reg33=reg34-reg33; reg34=pow(reg2,2); T reg45=reg16*reg22;
    reg39=reg39+reg36; T reg46=pow(reg6,2); reg37=reg38-reg37; reg38=reg13*reg26; T reg47=pow(reg4,2);
    T reg48=reg2*reg15; T reg49=reg6*reg14; T reg50=reg34*reg16; reg30=reg29-reg30; reg29=reg35*reg12;
    reg40=reg41+reg40; reg41=reg18*reg37; T reg51=reg39*reg16; T reg52=reg6*reg15; T reg53=reg2*reg14;
    T reg54=reg46*reg18; reg44=reg38+reg44; reg20=reg20*reg12; reg32=reg33-reg32; reg42=reg31-reg42;
    reg31=reg34*reg17; reg33=reg45+reg43; reg38=reg46*reg16; reg27=reg28-reg27; reg28=reg13*reg22;
    reg21=reg18*reg21; T reg55=reg13*reg25; T reg56=reg17*reg10; reg10=reg16*reg10; T reg57=reg3*reg1;
    reg44=reg20-reg44; reg40=reg29-reg40; reg20=reg34*reg11; reg29=reg46*reg13; T reg58=reg48+reg49;
    T reg59=reg4*reg5; reg51=reg41-reg51; reg41=reg15*reg14; T reg60=reg30*reg52; T reg61=reg42*reg53;
    T reg62=reg27*reg52; T reg63=reg32*reg53; T reg64=reg34*reg32; T reg65=reg46*reg27; T reg66=reg34*reg42;
    T reg67=reg46*reg30; T reg68=reg12*reg7; reg38=reg31-reg38; reg31=reg47*reg11; reg50=reg54-reg50;
    reg54=reg47*reg13; T reg69=reg24*reg32; T reg70=reg23*reg27; T reg71=reg24*reg42; T reg72=reg23*reg30;
    T reg73=reg13*reg7; reg7=reg11*reg7; T reg74=reg4*reg14; T reg75=reg33*reg13; T reg76=2*reg15;
    T reg77=reg2*reg5; T reg78=reg6*reg5; T reg79=reg4*reg15; T reg80=reg76*reg14; reg63=reg62+reg63;
    reg62=reg47*reg12; T reg81=reg59*reg40; reg61=reg60+reg61; reg60=reg35*reg40; reg71=reg72+reg71;
    reg20=reg29+reg20; reg29=reg68*reg17; reg72=reg8*reg2; T reg82=reg77-reg74; T reg83=reg57*reg13;
    reg69=reg70+reg69; reg70=reg35*reg44; T reg84=reg59*reg44; reg75=reg51-reg75; reg54=reg50-reg54;
    reg50=reg41*reg3; reg51=reg3*reg58; reg64=reg65+reg64; reg31=reg38-reg31; reg38=reg47*reg44;
    reg68=reg68*reg16; reg65=reg7*reg11; T reg85=reg78+reg79; T reg86=reg73*reg11; T reg87=reg57*reg12;
    T reg88=reg15*reg5; T reg89=reg56*reg13; reg36=reg28+reg36; reg55=reg21-reg55; reg22=reg18*reg22;
    reg21=reg10*reg13; reg57=reg57*reg11; reg25=reg16*reg25; reg66=reg67+reg66; reg28=reg47*reg40;
    reg67=reg3*reg0; reg78=reg79-reg78; reg79=reg76*reg5; T reg90=reg88*reg0; T reg91=reg0*reg85;
    reg55=reg55/reg75; T reg92=2*reg14; reg12=reg67*reg12; T reg93=reg46*reg54; T reg94=reg34*reg31;
    reg39=reg39/reg75; T reg95=reg24*reg31; reg38=reg64+reg38; reg64=reg83*reg11; T reg96=reg58*reg50;
    T reg97=2*reg82; reg77=reg74+reg77; reg89=reg45+reg89; reg74=reg67*reg11; T reg98=reg72*reg51;
    reg65=reg29-reg65; reg36=reg36/reg75; reg29=reg57*reg11; T reg99=reg72*reg50; T reg100=reg87*reg16;
    reg86=reg68+reg86; reg87=reg87*reg17; reg73=reg73*reg17; reg68=reg23*reg54; reg7=reg7*reg16;
    reg21=reg22+reg21; reg25=reg22+reg25; reg67=reg67*reg13; reg22=reg80*reg50; reg70=reg69+reg70;
    reg60=reg71+reg60; reg69=reg80*reg51; reg56=reg56*reg18; reg81=reg61+reg81; reg10=reg10*reg16;
    reg61=reg8*reg4; reg71=reg6*reg2; reg84=reg63+reg84; reg20=reg62-reg20; reg37=reg37/reg75;
    reg28=reg66+reg28; reg62=reg58*reg51; reg63=reg14*reg5; reg66=reg72*reg37; T reg101=reg24*reg39;
    T reg102=reg34*reg37; T reg103=reg31*reg53; reg29=reg87-reg29; reg87=reg61*reg90; T reg104=reg80*reg39;
    T reg105=pow(reg78,2); T reg106=reg54*reg52; reg64=reg100+reg64; reg100=pow(reg82,2); T reg107=reg6*reg4;
    T reg108=reg92*reg5; T reg109=reg12*reg17; reg12=reg12*reg16; T reg110=reg74*reg11; reg11=reg67*reg11;
    reg99=reg28+reg99; reg3=reg71*reg3; reg28=reg24*reg55; T reg111=reg34*reg36; T reg112=reg79*reg91;
    reg69=reg70+reg69; reg96=reg81+reg96; reg70=reg79*reg90; reg22=reg60+reg22; reg60=reg85*reg90;
    reg81=reg9*reg4; T reg113=reg23*reg55; reg83=reg83*reg17; reg57=reg57*reg16; T reg114=reg46*reg36;
    reg10=reg56-reg10; reg25=reg25/reg75; reg21=reg21/reg75; reg56=reg23*reg39; T reg115=reg46*reg37;
    reg33=reg33/reg75; reg89=reg89/reg75; T reg116=reg97*reg78; T reg117=reg7+reg73; reg98=reg38+reg98;
    reg86=reg86*reg16; reg38=reg61*reg91; reg65=reg65*reg18; T reg118=reg35*reg20; reg95=reg68+reg95;
    reg68=reg47*reg20; reg62=reg84+reg62; reg84=reg85*reg91; T reg119=reg72*reg36; T reg120=reg80*reg55;
    T reg121=reg63*reg1; T reg122=reg1*reg77; reg94=reg93+reg94; reg93=reg77*reg121; reg118=reg95+reg118;
    reg95=reg80*reg3; T reg123=reg81*reg122; reg38=reg98+reg38; reg60=reg96+reg60; reg96=reg59*reg20;
    reg103=reg106+reg103; reg98=reg77*reg122; reg87=reg99+reg87; reg84=reg62+reg84; reg62=reg81*reg121;
    reg64=reg64*reg16; reg29=reg29*reg18; reg117=reg117*reg13; reg86=reg65-reg86; reg17=reg67*reg17;
    reg65=reg72*reg3; reg67=reg23*reg21; reg99=reg46*reg89; reg68=reg94+reg68; reg94=reg25*reg116;
    reg120=reg119+reg120; reg106=reg25*reg105; reg0=reg107*reg0; reg28=reg111+reg28; reg111=reg108*reg122;
    reg112=reg69+reg112; reg69=reg25*reg100; reg119=reg108*reg121; reg70=reg22+reg70; reg113=reg114+reg113;
    reg74=reg74*reg16; reg10=reg10/reg75; reg22=reg2*reg4; reg114=reg34*reg89; T reg124=reg24*reg21;
    T reg125=reg72*reg89; T reg126=reg80*reg21; T reg127=reg33*reg116; reg104=reg66+reg104; reg56=reg115+reg56;
    reg66=reg33*reg105; reg115=reg57+reg83; reg110=reg109-reg110; reg11=reg12+reg11; reg12=reg33*reg100;
    reg101=reg102+reg101; reg102=reg61*reg0; reg109=reg10*reg116; reg1=reg22*reg1; reg126=reg125+reg126;
    reg65=reg68+reg65; reg68=reg58*reg3; reg96=reg103+reg96; reg12=reg56+reg12; reg56=reg10*reg105;
    reg124=reg114+reg124; reg103=reg79*reg0; reg114=reg10*reg100; reg95=reg118+reg95; reg119=reg70+reg119;
    reg69=reg113+reg69; reg111=reg112+reg111; reg106=reg28+reg106; reg94=reg120+reg94; reg67=reg99+reg67;
    reg117=reg86-reg117; reg64=reg29-reg64; reg115=reg115*reg13; reg18=reg110*reg18; reg16=reg11*reg16;
    reg11=reg74+reg17; reg66=reg101+reg66; reg127=reg104+reg127; reg62=reg87+reg62; reg93=reg60+reg93;
    reg98=reg84+reg98; reg123=reg38+reg123; reg115=reg64-reg115; reg102=reg65+reg102; reg28=reg81*reg1;
    reg29=reg127*reg71; reg38=reg41*reg94; reg16=reg18-reg16; reg103=reg95+reg103; reg18=reg6*reg78;
    reg60=reg82*reg2; reg13=reg11*reg13; reg109=reg126+reg109; reg114=reg67+reg114; reg11=reg82*reg78;
    reg64=reg71*reg12; reg65=reg111*reg93; reg67=reg62*reg98; reg56=reg124+reg56; reg70=reg119*reg98;
    reg84=reg123*reg93; reg68=reg96+reg68; reg86=reg85*reg0; reg87=reg108*reg1; reg95=reg41*reg106;
    reg117=reg117/reg75; reg96=reg66*reg71; reg99=reg41*reg69; reg101=reg23*reg106; reg13=reg16-reg13;
    reg65=reg70-reg65; reg16=reg34*reg127; reg70=reg24*reg94; reg104=reg15*reg78; reg103=reg87+reg103;
    reg87=reg2*reg78; reg110=reg46*reg66; reg112=reg82*reg6; reg18=reg60+reg18; reg60=reg119*reg123;
    reg113=reg82*reg14; reg106=reg24*reg106; reg66=reg34*reg66; reg118=reg24*reg69; reg120=reg34*reg12;
    reg94=reg23*reg94; reg127=reg46*reg127; reg84=reg67-reg84; reg67=reg111*reg62; reg115=reg115/reg75;
    reg124=reg117*reg53; reg28=reg102+reg28; reg69=reg23*reg69; reg12=reg46*reg12; reg102=reg117*reg58;
    reg86=reg68+reg86; reg68=reg77*reg1; reg125=reg11*reg56; reg126=reg117*reg52; reg95=reg96+reg95;
    reg99=reg64+reg99; reg64=reg11*reg114; reg38=reg29+reg38; reg29=reg11*reg109; reg96=elem.pos(1)[0]-elem.pos(0)[0];
    T reg128=reg100*reg109; reg94=reg127+reg94; reg127=elem.pos(1)[1]-elem.pos(0)[1]; T reg129=elem.pos(2)[0]-elem.pos(0)[0]; reg75=reg13/reg75;
    reg13=reg84*reg103; T reg130=reg115*reg87; T reg131=reg124*reg58; reg64=reg99+reg64; reg99=reg126*reg58;
    T reg132=reg115*reg18; reg68=reg86+reg68; reg86=reg100*reg114; T reg133=reg115*reg112; reg125=reg95+reg125;
    reg104=reg113+reg104; reg109=reg105*reg109; reg70=reg16+reg70; reg101=reg110+reg101; reg16=reg100*reg56;
    reg56=reg105*reg56; reg106=reg66+reg106; reg66=reg82*reg15; reg95=reg14*reg78; reg60=reg67-reg60;
    reg29=reg38+reg29; reg69=reg12+reg69; reg12=elem.pos(2)[1]-elem.pos(0)[1]; reg118=reg120+reg118; reg114=reg105*reg114;
    reg38=reg28*reg65; reg67=reg102*reg58; reg110=reg119*reg68; reg9=reg9*reg78; reg113=reg75*reg66;
    reg16=reg101+reg16; reg101=reg124*reg26; reg120=reg62*reg68; reg8=reg82*reg8; T reg134=reg103*reg62;
    T reg135=reg75*reg95; T reg136=reg75*reg104; T reg137=reg133*reg18; T reg138=reg132*reg18; T reg139=reg119*reg28;
    reg109=reg70+reg109; reg67=reg29+reg67; reg29=reg102*reg19; reg99=reg64+reg99; reg64=reg28*reg93;
    reg70=reg126*reg19; reg114=reg118+reg114; reg118=reg12*reg96; reg128=reg94+reg128; reg102=reg102*reg26;
    reg94=reg129*reg127; T reg140=reg130*reg18; reg131=reg125+reg131; reg13=reg38-reg13; reg38=reg60*reg68;
    reg124=reg124*reg19; reg56=reg106+reg56; reg86=reg69+reg86; reg126=reg126*reg26; reg69=reg103*reg93;
    reg102=reg128+reg102; reg106=reg136*reg104; reg138=reg67+reg138; reg67=reg132*reg8; reg70=reg114+reg70;
    reg76=reg82*reg76; reg114=reg135*reg104; reg125=reg133*reg9; reg140=reg131+reg140; reg128=reg130*reg9;
    reg124=reg56+reg124; reg94=reg118-reg94; reg130=reg130*reg8; reg101=reg16+reg101; reg29=reg109+reg29;
    reg132=reg132*reg9; reg133=reg133*reg8; reg126=reg86+reg126; reg16=reg103*reg123; reg56=reg111*reg28;
    reg120=reg64-reg120; reg64=reg113*reg104; reg137=reg99+reg137; reg92=reg92*reg78; reg86=reg123*reg68;
    reg38=reg13+reg38; reg13=reg103*reg98; reg99=reg111*reg68; reg109=reg28*reg98; reg134=reg139-reg134;
    reg110=reg69-reg110; reg69=reg135*reg92; reg130=reg101+reg130; reg135=reg135*reg76; reg106=reg138+reg106;
    reg101=1-(*f.m).resolution; reg125=reg70+reg125; reg70=reg113*reg92; reg114=reg140+reg114; reg128=reg124+reg128;
    reg129=reg129/reg94; reg132=reg29+reg132; reg64=reg137+reg64; reg12=reg12/reg94; reg16=reg56-reg16;
    reg134=reg134/reg38; reg29=reg136*reg92; reg110=reg110/reg38; reg96=reg96/reg94; reg86=reg109-reg86;
    reg120=reg120/reg38; reg127=reg127/reg94; reg99=reg13-reg99; reg133=reg126+reg133; reg113=reg113*reg76;
    reg67=reg102+reg67; reg136=reg136*reg76; reg13=reg129-reg96; reg99=reg99/reg38; reg65=reg65/reg38;
    reg29=reg132+reg29; reg56=reg110*(*f.m).resolution; reg102=reg120*(*f.m).resolution; reg109=reg134*(*f.m).resolution; reg118=reg127-reg12;
    reg113=reg133+reg113; reg136=reg67+reg136; reg135=reg130+reg135; reg70=reg125+reg70; reg16=reg16/reg38;
    reg60=reg60/reg38; reg67=reg64*reg101; reg124=reg101*reg114; reg125=reg101*reg106; reg86=reg86/reg38;
    reg38=reg84/reg38; reg69=reg128+reg69; reg84=reg29*reg101; reg126=0.5*reg129; reg128=0.5*reg13;
    reg130=reg16*(*f.m).resolution; reg131=reg86*(*f.m).resolution; reg132=reg99*(*f.m).resolution; reg133=0.5*reg12; reg56=reg67+reg56;
    reg102=reg124-reg102; reg109=reg125+reg109; reg67=reg101*reg69; reg124=reg38*(*f.m).resolution; reg125=reg65*(*f.m).resolution;
    reg137=0.5*reg127; reg138=0.5*reg118; reg139=reg60*(*f.m).resolution; reg140=0.5*reg96; T reg141=reg113*reg101;
    T reg142=reg135*reg101; T reg143=reg101*reg136; T reg144=reg101*reg70; reg67=reg131+reg67; reg132=reg144-reg132;
    reg130=reg84-reg130; reg143=reg139+reg143; reg124=reg142-reg124; reg141=reg125+reg141; reg84=reg133*reg109;
    reg125=reg102*reg129; reg131=reg102*reg96; reg139=reg109*reg137; reg142=reg128*reg109; reg144=reg118*reg56;
    T reg145=reg109*reg126; T reg146=reg56*reg12; T reg147=reg56*reg127; T reg148=reg109*reg140; T reg149=reg13*reg102;
    T reg150=reg109*reg138; T reg151=reg13*reg124; T reg152=reg143*reg138; reg150=reg149+reg150; reg131=reg131-reg139;
    reg142=reg144+reg142; reg144=reg133*reg130; reg149=reg67*reg129; T reg153=reg118*reg141; T reg154=reg130*reg126;
    T reg155=reg132*reg12; T reg156=reg128*reg143; T reg157=reg13*reg67; T reg158=reg130*reg138; T reg159=reg143*reg137;
    T reg160=reg124*reg96; reg148=reg148-reg147; T reg161=reg143*reg140; T reg162=reg141*reg127; reg84=reg84-reg125;
    T reg163=reg130*reg137; T reg164=reg67*reg96; T reg165=reg130*reg140; T reg166=reg132*reg127; T reg167=reg141*reg12;
    T reg168=reg143*reg126; reg146=reg146-reg145; T reg169=reg124*reg129; T reg170=reg133*reg143; reg144=reg144-reg149;
    reg167=reg167-reg168; reg156=reg153+reg156; reg146=2*reg146; reg155=reg155-reg154; reg142=2*reg142;
    reg158=reg157+reg158; reg131=2*reg131; reg170=reg170-reg169; reg165=reg165-reg166; reg160=reg160-reg159;
    reg164=reg164-reg163; reg150=2*reg150; reg84=2*reg84; reg152=reg151+reg152; reg148=2*reg148;
    reg161=reg161-reg162; reg151=reg150*reg138; reg153=reg133*reg131; reg157=reg164*reg129; T reg171=reg138*reg146;
    T reg172=reg142*reg128; T reg173=reg128*reg148; T reg174=reg150*reg128; T reg175=reg167*reg12; T reg176=reg137*reg131;
    T reg177=reg155*reg13; T reg178=reg84*reg138; T reg179=reg164*reg96; T reg180=reg118*reg167; T reg181=reg160*reg127;
    T reg182=reg140*reg131; T reg183=reg118*reg161; T reg184=reg161*reg127; T reg185=reg165*reg129; T reg186=reg133*reg148;
    T reg187=reg144*reg129; T reg188=reg133*reg84; T reg189=reg148*reg138; T reg190=reg118*reg152; T reg191=reg128*reg84;
    T reg192=reg13*reg165; T reg193=reg138*reg131; reg164=reg13*reg164; T reg194=reg140*reg148; T reg195=reg118*reg170;
    T reg196=reg118*reg160; T reg197=reg126*reg131; reg160=reg160*reg12; T reg198=reg148*reg126; reg161=reg161*reg12;
    reg131=reg128*reg131; T reg199=reg128*reg146; T reg200=reg118*reg156; T reg201=reg13*reg158; T reg202=reg84*reg126;
    T reg203=reg146*reg126; T reg204=reg170*reg12; T reg205=reg144*reg13; reg187=reg188-reg187; reg131=reg196+reg131;
    reg203=reg175-reg203; reg172=reg200+reg172; reg192=reg189+reg192; reg197=reg160-reg197; reg164=reg193+reg164;
    reg195=reg191+reg195; reg198=reg161-reg198; reg174=reg190+reg174; reg173=reg183+reg173; reg202=reg204-reg202;
    reg176=reg179-reg176; reg184=reg194-reg184; reg205=reg178+reg205; reg177=reg171+reg177; reg201=reg151+reg201;
    reg181=reg182-reg181; reg199=reg180+reg199; reg185=reg186-reg185; reg157=reg153-reg157; reg197=reg197*reg94;
    reg201=reg201*reg94; reg164=reg164*reg94; reg174=reg174*reg94; reg195=reg195*reg94; reg198=reg198*reg94;
    reg199=reg199*reg94; reg157=reg157*reg94; reg202=reg202*reg94; reg205=reg205*reg94; reg173=reg173*reg94;
    reg176=reg176*reg94; reg177=reg177*reg94; reg203=reg203*reg94; reg181=reg181*reg94; reg172=reg172*reg94;
    reg131=reg131*reg94; reg192=reg192*reg94; reg184=reg184*reg94; reg187=reg187*reg94; reg185=reg185*reg94;
    T tmp_2_3=ponderation*reg202; T tmp_4_4=ponderation*reg184; T tmp_3_3=ponderation*reg187; T tmp_0_4=ponderation*reg173; T tmp_5_5=ponderation*reg176;
    T tmp_3_5=ponderation*reg157; T tmp_1_2=ponderation*reg177; T tmp_0_5=ponderation*reg131; T tmp_2_2=ponderation*reg203; T tmp_1_3=ponderation*reg205;
    T tmp_4_5=ponderation*reg181; T tmp_0_3=ponderation*reg195; T tmp_0_2=ponderation*reg199; T tmp_2_4=ponderation*reg198; T tmp_1_1=ponderation*reg201;
    T tmp_3_4=ponderation*reg185; T tmp_1_5=ponderation*reg164; T tmp_1_4=ponderation*reg192; T tmp_2_5=ponderation*reg197; T tmp_0_0=ponderation*reg172;
    T tmp_0_1=ponderation*reg174;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; reg1=1.0/reg1; T reg2=2*(*f.m).shear_modulus_12;
    T reg3=reg1*reg0; reg2=1.0/reg2; T reg4=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg5=1.0/(*f.m).elastic_modulus_3; T reg6=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg7=reg2*reg3; T reg8=pow((*f.m).v1[0],2); T reg9=pow((*f.m).v2[1],2); T reg10=pow((*f.m).v2[0],2); T reg11=pow((*f.m).v1[1],2);
    T reg12=pow((*f.m).v2[2],2); reg9=reg10+reg9; reg10=1.0/(*f.m).elastic_modulus_2; T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg14=pow((*f.m).v1[2],2);
    reg11=reg8+reg11; reg8=reg5*reg7; T reg15=reg6*reg7; T reg16=reg4*reg7; reg14=reg11+reg14;
    reg11=reg10*reg8; T reg17=reg6*reg16; T reg18=reg6*reg15; T reg19=reg13*reg8; reg12=reg9+reg12;
    reg18=reg11-reg18; reg9=reg10*reg16; reg11=1.0/(*f.m).elastic_modulus_1; reg17=reg17+reg19; T reg20=reg13*reg15;
    reg14=pow(reg14,0.5); reg12=pow(reg12,0.5); T reg21=(*f.m).v2[2]/reg12; T reg22=(*f.m).v2[1]/reg12; T reg23=reg20+reg9;
    T reg24=(*f.m).v1[2]/reg14; T reg25=(*f.m).v1[1]/reg14; T reg26=reg17*reg13; T reg27=reg11*reg18; T reg28=reg23*reg4;
    T reg29=reg24*reg22; T reg30=reg10*reg7; T reg31=reg4*reg15; reg8=reg11*reg8; T reg32=reg4*reg16;
    T reg33=reg25*reg21; reg7=reg13*reg7; T reg34=reg2*reg0; T reg35=reg6*reg3; T reg36=reg4*reg3;
    reg26=reg27-reg26; reg14=(*f.m).v1[0]/reg14; reg12=(*f.m).v2[0]/reg12; reg3=reg5*reg3; reg28=reg26-reg28;
    reg26=2*reg12; reg27=reg24*reg12; T reg37=reg14*reg21; T reg38=2*reg14; T reg39=reg33-reg29;
    T reg40=reg34*reg5; T reg41=reg36*reg6; T reg42=reg35*reg6; T reg43=reg30*reg4; T reg44=reg3*reg13;
    reg19=reg31+reg19; reg3=reg3*reg10; reg32=reg8-reg32; reg15=reg11*reg15; reg8=reg7*reg4;
    reg31=reg34*reg4; reg34=reg34*reg6; reg16=reg13*reg16; T reg45=reg2*reg1; reg43=reg20+reg43;
    reg19=reg19/reg28; reg35=reg35*reg13; T reg46=reg27-reg37; T reg47=reg45*reg6; reg36=reg36*reg10;
    T reg48=2*reg39; reg18=reg18/reg28; reg17=reg17/reg28; T reg49=pow(reg12,2); T reg50=reg14*reg22;
    T reg51=pow(reg22,2); T reg52=reg45*reg4; T reg53=pow(reg14,2); T reg54=reg26*reg22; T reg55=reg25*reg12;
    reg8=reg15+reg8; T reg56=reg40*reg10; reg40=reg40*reg13; T reg57=reg34*reg6; reg41=reg44+reg41;
    reg44=reg31*reg6; T reg58=reg38*reg25; reg16=reg15+reg16; reg45=reg45*reg5; reg7=reg7*reg13;
    reg15=pow(reg25,2); reg42=reg3-reg42; reg32=reg32/reg28; reg30=reg30*reg11; reg7=reg30-reg7;
    reg3=pow(reg24,2); reg16=reg16/reg28; reg34=reg34*reg13; reg30=reg54*reg17; T reg59=reg58*reg18;
    T reg60=reg50-reg55; reg31=reg31*reg10; T reg61=reg51*reg17; T reg62=reg15*reg18; T reg63=reg52*reg6;
    reg42=reg42*reg11; T reg64=reg47*reg6; T reg65=reg45*reg13; reg45=reg45*reg10; reg41=reg41*reg13;
    reg44=reg40+reg44; reg57=reg56-reg57; reg40=reg35+reg36; reg56=reg49*reg32; T reg66=reg53*reg19;
    T reg67=reg54*reg32; T reg68=reg58*reg19; T reg69=reg48*reg46; T reg70=pow(reg46,2); T reg71=pow(reg39,2);
    reg43=reg43/reg28; reg23=reg23/reg28; T reg72=pow(reg21,2); T reg73=reg51*reg32; T reg74=reg15*reg19;
    T reg75=reg53*reg18; T reg76=reg49*reg17; reg8=reg8/reg28; reg64=reg45-reg64; reg52=reg52*reg10;
    reg73=reg74+reg73; reg45=reg16*reg70; reg74=reg16*reg69; T reg77=reg58*reg43; reg41=reg42-reg41;
    reg42=reg34+reg31; T reg78=reg54*reg8; reg44=reg44*reg13; reg67=reg68+reg67; reg57=reg57*reg11;
    reg68=reg3*reg19; T reg79=reg72*reg32; reg40=reg40*reg4; reg7=reg7/reg28; T reg80=pow(reg60,2);
    reg76=reg75+reg76; reg75=reg23*reg71; T reg81=reg23*reg69; reg30=reg59+reg30; reg59=reg53*reg43;
    T reg82=reg49*reg8; reg56=reg66+reg56; reg66=reg72*reg17; T reg83=reg3*reg18; T reg84=reg15*reg43;
    T reg85=reg23*reg70; reg61=reg62+reg61; reg62=reg51*reg8; T reg86=reg16*reg71; reg47=reg47*reg13;
    reg63=reg65+reg63; reg45=reg73+reg45; reg74=reg67+reg74; reg86=reg56+reg86; reg56=reg16*reg80;
    reg79=reg68+reg79; reg65=reg7*reg69; reg78=reg77+reg78; reg67=reg72*reg8; reg68=reg3*reg43;
    reg73=reg7*reg70; reg62=reg84+reg62; reg77=reg7*reg71; reg82=reg59+reg82; reg59=2*reg22;
    reg84=reg26*reg21; T reg87=2*reg25; reg75=reg76+reg75; reg76=reg38*reg24; T reg88=reg12*reg22;
    T reg89=reg14*reg25; reg64=reg64*reg11; reg63=reg63*reg13; T reg90=reg47+reg52; reg85=reg61+reg85;
    reg42=reg42*reg4; reg66=reg83+reg66; reg61=reg23*reg80; reg81=reg30+reg81; reg44=reg57-reg44;
    reg40=reg41-reg40; reg50=reg55+reg50; reg30=reg25*reg22; reg41=reg53*reg81; reg55=reg14*reg12;
    reg57=reg49*reg74; reg83=reg15*reg75; T reg91=reg51*reg86; reg42=reg44-reg42; reg44=2*reg46;
    reg40=reg40/reg28; T reg92=reg15*reg85; T reg93=reg51*reg45; T reg94=reg76*reg18; T reg95=reg17*reg84;
    T reg96=reg85*reg89; T reg97=reg88*reg45; T reg98=reg88*reg74; T reg99=reg81*reg89; T reg100=reg87*reg22;
    T reg101=2*reg24; reg61=reg66+reg61; reg66=reg39*reg46; T reg102=reg38*reg12; T reg103=reg87*reg24;
    reg90=reg90*reg4; T reg104=reg59*reg21; T reg105=reg14*reg46; reg63=reg64-reg63; reg48=reg48*reg60;
    reg81=reg15*reg81; reg74=reg51*reg74; reg77=reg82+reg77; reg73=reg62+reg73; reg67=reg68+reg67;
    reg62=reg7*reg80; reg65=reg78+reg65; reg64=reg76*reg19; reg68=reg32*reg84; reg78=reg53*reg75;
    reg45=reg49*reg45; reg85=reg53*reg85; reg82=reg49*reg86; T reg106=reg39*reg25; reg56=reg79+reg56;
    reg45=reg85+reg45; reg79=reg71*reg73; reg90=reg63-reg90; reg63=reg53*reg61; reg85=reg71*reg77;
    T reg107=reg49*reg56; T reg108=reg39*reg22; T reg109=reg12*reg46; T reg110=reg51*reg56; reg42=reg42/reg28;
    reg74=reg81+reg74; reg81=reg70*reg65; T reg111=reg100*reg10; T reg112=reg13*reg102; T reg113=reg51*reg10;
    T reg114=reg49*reg13; T reg115=reg100*reg13; T reg116=reg11*reg102; T reg117=reg51*reg13; T reg118=reg49*reg11;
    reg19=reg103*reg19; T reg119=reg70*reg73; reg93=reg92+reg93; reg92=reg70*reg77; reg91=reg83+reg91;
    reg83=reg25*reg46; T reg120=reg8*reg84; T reg121=reg71*reg65; reg57=reg41+reg57; reg41=reg76*reg43;
    T reg122=reg39*reg14; reg62=reg67+reg62; reg67=reg24*reg21; reg44=reg44*reg60; reg32=reg32*reg104;
    T reg123=reg40*reg55; T reg124=reg40*reg30; T reg125=reg15*reg61; T reg126=reg40*reg50; reg101=reg101*reg21;
    reg98=reg99+reg98; reg65=reg66*reg65; reg73=reg66*reg73; reg97=reg96+reg97; reg82=reg78+reg82;
    reg86=reg88*reg86; reg75=reg89*reg75; reg105=reg106+reg105; reg78=reg23*reg48; reg68=reg64+reg68;
    reg95=reg94+reg95; reg18=reg103*reg18; reg17=reg17*reg104; reg64=reg16*reg48; reg94=reg4*reg102;
    reg115=reg116-reg115; reg96=reg101*reg4; reg99=reg51*reg6; reg106=reg49*reg4; reg117=reg118-reg117;
    reg116=reg72*reg4; reg23=reg23*reg44; reg17=reg18+reg17; reg78=reg95+reg78; reg86=reg75+reg86;
    reg77=reg66*reg77; reg16=reg16*reg44; reg32=reg19+reg32; reg18=reg70*reg62; reg121=reg57+reg121;
    reg81=reg74+reg81; reg19=reg126*reg100; reg125=reg110+reg125; reg57=reg22*reg46; reg28=reg90/reg28;
    reg74=reg71*reg62; reg107=reg63+reg107; reg63=reg39*reg12; reg75=reg124*reg102; reg79=reg45+reg79;
    reg45=reg123*reg102; reg112=reg111-reg112; reg90=reg126*reg102; reg95=reg101*reg6; reg85=reg82+reg85;
    reg114=reg113-reg114; reg82=reg72*reg6; reg110=reg100*reg6; reg111=reg40*reg67; reg113=reg42*reg122;
    reg118=reg42*reg83; reg92=reg91+reg92; reg91=reg123*reg100; reg126=reg126*reg50; reg65=reg98+reg65;
    reg98=reg42*reg105; reg56=reg88*reg56; reg61=reg61*reg89; T reg127=reg24*reg60; reg119=reg93+reg119;
    reg93=reg124*reg50; reg73=reg97+reg73; reg124=reg124*reg100; reg64=reg68+reg64; reg87=reg87*reg46;
    reg109=reg108+reg109; reg120=reg41+reg120; reg38=reg39*reg38; reg41=reg7*reg48; reg43=reg103*reg43;
    reg8=reg8*reg104; reg68=(*f.m).alpha_1*reg53; reg97=reg118*reg105; reg93=reg73+reg93; reg75=reg79+reg75;
    reg73=reg118*reg38; reg118=reg118*reg87; reg79=reg111*reg102; reg77=reg86+reg77; reg123=reg123*reg50;
    reg11=reg53*reg11; reg74=reg107+reg74; reg86=reg28*reg63; reg16=reg32+reg16; reg32=reg28*reg57;
    reg59=reg59*reg46; reg18=reg125+reg18; reg107=reg111*reg100; reg7=reg7*reg44; reg23=reg17+reg23;
    reg27=reg37+reg27; reg17=reg42*reg127; reg8=reg43+reg8; reg37=reg72*reg5; reg99=reg106+reg99;
    reg43=reg98*reg105; reg101=reg101*reg5; reg110=reg94+reg110; reg126=reg65+reg126; reg26=reg39*reg26;
    reg65=reg21*reg60; reg94=reg28*reg109; reg62=reg66*reg62; reg56=reg61+reg56; reg61=(*f.m).alpha_2*reg49;
    reg41=reg120+reg41; reg45=reg85+reg45; reg85=reg113*reg38; reg106=reg113*reg87; reg96=reg115-reg96;
    reg108=reg53*reg13; reg116=reg117-reg116; reg10=reg15*reg10; reg115=reg98*reg38; reg13=reg15*reg13;
    reg117=reg53*reg78; reg120=reg49*reg64; reg82=reg114-reg82; reg90=reg121+reg90; reg95=reg112-reg95;
    reg112=(*f.m).alpha_2*reg51; reg114=(*f.m).alpha_1*reg15; reg91=reg92+reg91; reg92=reg51*reg64; reg19=reg81+reg19;
    reg98=reg98*reg87; reg81=reg14*reg60; reg121=reg15*reg78; reg124=reg119+reg124; reg119=reg39*reg24;
    reg64=reg88*reg64; reg125=reg51*reg82; T reg128=reg40*reg27; T reg129=reg82*reg30; T reg130=reg94*reg26;
    T reg131=reg49*reg116; T reg132=reg3*reg4; reg13=reg11-reg13; reg81=reg119+reg81; reg115=reg90+reg115;
    reg7=reg8+reg7; reg33=reg29+reg33; reg4=reg53*reg4; reg8=reg15*reg6; reg111=reg111*reg50;
    reg62=reg56+reg62; reg11=reg49*reg16; reg29=reg94*reg109; reg56=reg53*reg23; reg43=reg126+reg43;
    reg90=reg32*reg59; reg119=reg12*reg60; reg126=reg39*reg21; T reg133=reg15*reg23; T reg134=reg51*reg16;
    T reg135=reg49*reg96; T reg136=reg51*reg95; T reg137=reg71*(*f.m).alpha_3; T reg138=reg116*reg55; reg68=reg61+reg68;
    reg113=reg113*reg105; reg123=reg77+reg123; reg61=reg41*reg71; reg78=reg89*reg78; reg120=reg117+reg120;
    reg77=reg24*reg46; reg85=reg45+reg85; reg45=reg86*reg26; reg117=(*f.m).alpha_1*reg3; T reg139=reg25*reg60;
    T reg140=reg41*reg70; reg92=reg121+reg92; reg73=reg75+reg73; reg118=reg124+reg118; reg75=reg32*reg26;
    reg121=reg95*reg30; reg94=reg94*reg59; reg98=reg19+reg98; reg79=reg74+reg79; reg19=reg96*reg55;
    reg74=(*f.m).alpha_2*reg72; reg107=reg18+reg107; reg18=reg17*reg87; reg124=reg17*reg38; T reg141=reg28*reg65;
    reg110=reg101-reg110; reg6=reg3*reg6; reg108=reg10-reg108; reg10=reg70*(*f.m).alpha_3; reg32=reg32*reg109;
    reg101=reg53*reg96; T reg142=reg15*reg95; reg106=reg91+reg106; reg97=reg93+reg97; reg112=reg114+reg112;
    reg91=reg15*reg82; reg93=reg53*reg116; reg114=reg86*reg59; reg99=reg37-reg99; reg137=reg68+reg137;
    reg74=reg117+reg74; reg37=reg80*(*f.m).alpha_3; reg41=reg41*reg66; reg64=reg78+reg64; reg111=reg62+reg111;
    reg40=reg40*reg33; reg17=reg17*reg105; reg29=reg43+reg29; reg43=reg21*reg46; reg121=reg19+reg121;
    reg32=reg97+reg32; reg129=reg138+reg129; reg19=reg67*reg110; reg62=reg42*reg81; reg14=reg14*reg24;
    reg12=reg12*reg21; reg68=reg67*reg99; reg134=reg133+reg134; reg78=reg7*reg70; reg10=reg112+reg10;
    reg114=reg106+reg114; reg132=reg13-reg132; reg13=(*f.m).alpha_1*reg89; reg5=reg3*reg5; reg8=reg4+reg8;
    reg4=(*f.m).alpha_2*reg88; reg97=reg3*reg110; reg142=reg101+reg142; reg101=reg88*reg2; reg106=reg2*reg50;
    reg112=reg3*reg99; reg91=reg93+reg91; reg6=reg108-reg6; reg93=reg128*reg100; reg45=reg85+reg45;
    reg139=reg77+reg139; reg140=reg92+reg140; reg75=reg73+reg75; reg94=reg98+reg94; reg73=reg141*reg26;
    reg18=reg107+reg18; reg77=reg141*reg59; reg124=reg79+reg124; reg79=reg72*reg99; reg125=reg131+reg125;
    reg85=reg7*reg71; reg11=reg56+reg11; reg90=reg118+reg90; reg119=reg126+reg119; reg56=reg128*reg102;
    reg16=reg88*reg16; reg136=reg135+reg136; reg88=reg72*reg110; reg61=reg120+reg61; reg86=reg86*reg109;
    reg23=reg89*reg23; reg113=reg123+reg113; reg130=reg115+reg130; reg92=reg22*reg60; reg85=reg11+reg85;
    reg73=reg124+reg73; reg42=reg42*reg139; reg11=reg130*reg32; reg92=reg43+reg92; reg102=reg40*reg102;
    reg43=reg15*reg6; reg98=reg53*reg132; reg39=reg39*reg60; reg24=reg25*reg24; reg21=reg22*reg21;
    reg22=reg1*reg27; reg25=reg12*reg1; reg107=reg62*reg38; reg56=reg61+reg56; reg61=reg10*reg75;
    reg8=reg5-reg8; reg5=reg137*reg45; reg128=reg128*reg50; reg16=reg23+reg16; reg7=reg7*reg66;
    reg23=reg137*reg114; reg108=reg10*reg90; reg86=reg113+reg86; reg113=reg51*reg6; reg115=reg49*reg132;
    reg117=reg58*reg106; reg97=reg142+reg97; reg118=reg58*reg101; reg112=reg91+reg112; reg77=reg18+reg77;
    reg12=(*f.m).alpha_2*reg12; reg18=(*f.m).alpha_1*reg14; reg91=reg94*reg32; reg66=reg66*(*f.m).alpha_3; reg120=reg75*reg29;
    reg4=reg13+reg4; reg37=reg74+reg37; reg13=reg54*reg106; reg88=reg136+reg88; reg74=reg54*reg101;
    reg79=reg125+reg79; reg93=reg140+reg93; reg123=reg62*reg87; reg78=reg134+reg78; reg100=reg40*reg100;
    reg124=reg90*reg29; reg68=reg129+reg68; reg125=reg50*reg101; reg19=reg121+reg19; reg121=reg50*reg106;
    reg17=reg111+reg17; reg141=reg141*reg109; reg111=reg28*reg119; reg41=reg64+reg41; reg64=reg75*reg94;
    reg11=reg120-reg11; reg66=reg4+reg66; reg46=reg60*reg46; reg12=reg18+reg12; reg39=reg39*(*f.m).alpha_3;
    reg4=(*f.m).alpha_1*reg24; reg3=reg3*reg8; reg43=reg98+reg43; reg18=(*f.m).alpha_2*reg21; reg60=reg0*reg33;
    reg21=reg21*reg0; reg91=reg124-reg91; reg141=reg17+reg141; reg2=reg89*reg2; reg118=reg112+reg118;
    reg17=reg76*reg25; reg117=reg97+reg117; reg89=reg76*reg22; reg87=reg42*reg87; reg100=reg78+reg100;
    reg78=reg132*reg55; reg125=reg68+reg125; reg68=reg111*reg59; reg123=reg93+reg123; reg93=reg130*reg90;
    reg38=reg42*reg38; reg102=reg85+reg102; reg85=reg27*reg25; reg74=reg79+reg74; reg79=reg84*reg25;
    reg97=reg111*reg26; reg107=reg56+reg107; reg128=reg41+reg128; reg13=reg88+reg13; reg41=reg84*reg22;
    reg62=reg62*reg105; reg121=reg19+reg121; reg19=reg27*reg22; reg40=reg40*reg50; reg56=reg137*reg86;
    reg88=reg10*reg32; reg28=reg28*reg92; reg72=reg72*reg8; reg113=reg115+reg113; reg7=reg16+reg7;
    reg61=reg5+reg61; reg77=reg37*reg77; reg73=reg37*reg73; reg5=reg6*reg30; reg108=reg23+reg108;
    reg16=reg11*reg114; reg19=reg121+reg19; reg23=reg104*reg60; reg41=reg13+reg41; reg42=reg42*reg105;
    reg13=reg104*reg21; reg79=reg74+reg79; reg74=reg33*reg60; reg67=reg67*reg8; reg40=reg7+reg40;
    reg97=reg107+reg97; reg7=reg33*reg21; reg85=reg125+reg85; reg88=reg56+reg88; reg141=reg37*reg141;
    reg38=reg102+reg38; reg26=reg28*reg26; reg93=reg64-reg93; reg68=reg123+reg68; reg5=reg78+reg5;
    reg73=reg61+reg73; reg87=reg100+reg87; reg59=reg28*reg59; reg72=reg113+reg72; reg1=reg14*reg1;
    reg111=reg111*reg109; reg14=reg54*reg2; reg46=reg46*(*f.m).alpha_3; reg18=reg4+reg18; reg4=reg66*reg94;
    reg3=reg43+reg3; reg17=reg118+reg17; reg43=reg58*reg2; reg56=reg103*reg21; reg61=reg66*reg130;
    reg64=reg103*reg60; reg78=reg45*reg91; reg77=reg108+reg77; reg39=reg12+reg39; reg89=reg117+reg89;
    reg62=reg128+reg62; reg111=reg62+reg111; reg59=reg87+reg59; reg26=reg38+reg26; reg12=reg66*reg29;
    reg38=reg45*reg29; reg62=reg93*reg86; reg87=reg114*reg29; reg98=reg94*reg86; reg97=reg39*reg97;
    reg100=reg86*reg130; reg56=reg17+reg56; reg64=reg89+reg64; reg74=reg19+reg74; reg67=reg5+reg67;
    reg68=reg39*reg68; reg16=reg78-reg16; reg5=reg50*reg2; reg23=reg41+reg23; reg76=reg76*reg1;
    reg43=reg3+reg43; reg28=reg28*reg109; reg42=reg40+reg42; reg7=reg85+reg7; reg13=reg79+reg13;
    reg4=reg77+reg4; reg84=reg84*reg1; reg14=reg72+reg14; reg61=reg73+reg61; reg0=reg24*reg0;
    reg46=reg18+reg46; reg141=reg88+reg141; reg27=reg27*reg1; reg5=reg67+reg5; reg28=reg42+reg28;
    reg3=reg13*reg74; reg17=reg56*reg74; reg111=reg39*reg111; reg12=reg141+reg12; reg59=reg46*reg59;
    reg68=reg4+reg68; reg84=reg14+reg84; reg104=reg104*reg0; reg4=reg64*reg7; reg26=reg46*reg26;
    reg97=reg61+reg97; reg14=reg130*reg114; reg18=reg45*reg94; reg19=reg75*reg86; reg100=reg38-reg100;
    reg24=reg45*reg32; reg38=reg86*reg90; reg98=reg87-reg98; reg40=reg114*reg32; reg62=reg16+reg62;
    reg76=reg43+reg76; reg103=reg103*reg0; reg16=reg23*reg7; reg28=reg46*reg28; reg111=reg12+reg111;
    reg38=reg40-reg38; reg11=reg11/reg62; reg59=reg68+reg59; reg100=reg100/reg62; reg103=reg76+reg103;
    reg84=reg104+reg84; reg16=reg3-reg16; reg26=reg97+reg26; reg3=reg75*reg114; reg14=reg18-reg14;
    reg12=reg45*reg90; reg19=reg24-reg19; reg91=reg91/reg62; reg33=reg33*reg0; reg98=reg98/reg62;
    reg18=reg13*reg64; reg27=reg5+reg27; reg5=reg23*reg56; reg4=reg17-reg4; reg100=reg100*reg59;
    reg98=reg98*reg26; reg11=reg11*reg59; reg91=reg91*reg26; reg19=reg19/reg62; reg17=reg4*reg84;
    reg93=reg93/reg62; reg3=reg12-reg3; reg33=reg27+reg33; reg12=reg103*reg16; reg38=reg38/reg62;
    reg14=reg14/reg62; reg18=reg5-reg18; reg28=reg111+reg28; reg17=reg12-reg17; reg5=reg18*reg33;
    reg93=reg93*reg28; reg11=reg91-reg11; reg14=reg14*reg28; reg98=reg100-reg98; reg26=reg38*reg26;
    reg12=reg84*reg7; reg59=reg19*reg59; reg19=reg13*reg33; reg24=reg103*reg7; reg27=reg56*reg33;
    reg62=reg3/reg62; reg3=reg23*reg33; reg38=reg103*reg74; reg19=reg12-reg19; reg27=reg24-reg27;
    reg12=reg64*reg33; reg5=reg17+reg5; reg17=reg84*reg56; reg24=reg84*reg74; reg40=reg13*reg103;
    reg41=1-(*f.m).resolution; reg59=reg26-reg59; reg28=reg62*reg28; reg11=reg93+reg11; reg14=reg98-reg14;
    reg26=reg137*(*f.m).resolution; reg42=reg84*reg64; reg43=reg23*reg103; reg61=reg10*(*f.m).resolution; reg14=reg14*reg41;
    reg17=reg40-reg17; reg3=reg24-reg3; reg19=reg19/reg5; reg12=reg38-reg12; reg27=reg27/reg5;
    reg11=reg11*reg41; reg59=reg28+reg59; reg59=reg59*reg41; reg11=reg26+reg11; reg14=reg61+reg14;
    reg16=reg16/reg5; reg12=reg12/reg5; reg4=reg4/reg5; reg42=reg43-reg42; reg3=reg3/reg5;
    reg24=reg66*(*f.m).resolution; reg17=reg17/reg5; reg26=reg19*(*f.m).resolution; reg86=reg86*reg41; reg32=reg41*reg32;
    reg28=reg27*(*f.m).resolution; reg75=reg75*reg41; reg38=reg17*(*f.m).resolution; reg114=reg41*reg114; reg18=reg18/reg5;
    reg90=reg41*reg90; reg40=elem.pos(2)[1]-elem.pos(0)[1]; reg43=reg3*(*f.m).resolution; reg61=reg12*(*f.m).resolution; reg45=reg45*reg41;
    reg26=reg86+reg26; reg62=reg16*(*f.m).resolution; reg14=(*f.m).deltaT*reg14; reg11=(*f.m).deltaT*reg11; reg59=reg24+reg59;
    reg29=reg41*reg29; reg24=elem.pos(1)[0]-elem.pos(0)[0]; reg67=elem.pos(2)[0]-elem.pos(0)[0]; reg68=elem.pos(1)[1]-elem.pos(0)[1]; reg5=reg42/reg5;
    reg42=reg4*(*f.m).resolution; reg28=reg32-reg28; reg130=reg41*reg130; reg59=(*f.m).deltaT*reg59; reg38=reg29+reg38;
    reg29=reg5*(*f.m).resolution; reg90=reg61+reg90; reg43=reg114-reg43; reg42=reg75-reg42; reg45=reg62+reg45;
    reg41=reg94*reg41; reg32=reg40*reg24; reg61=reg26*reg11; reg62=reg28*reg14; reg72=reg67*reg68;
    reg73=reg18*(*f.m).resolution; reg75=reg38*reg59; reg76=reg43*reg11; reg77=reg90*reg14; reg78=reg45*reg11;
    reg79=reg42*reg14; reg85=reg61+reg62; reg72=reg32-reg72; reg29=reg41-reg29; reg130=reg73+reg130;
    reg24=reg24/reg72; reg68=reg68/reg72; reg40=reg40/reg72; reg67=reg67/reg72; reg32=reg130*reg59;
    reg41=reg85+reg75; reg73=reg29*reg59; reg86=reg77+reg76; reg87=reg79+reg78; reg88=reg87+reg32;
    reg89=reg86+reg73; reg91=2*reg41; reg93=0.5*reg24; reg94=0.5*reg68; reg97=0.5*reg40;
    reg98=reg68-reg40; reg100=reg67-reg24; reg102=0.5*reg67; reg104=0.5*reg98; reg107=reg88*reg40;
    reg108=reg91*reg102; reg111=reg91*reg97; reg112=0.5*reg100; reg113=reg89*reg67; reg114=reg91*reg94;
    reg115=reg89*reg24; reg117=reg91*reg93; reg118=reg88*reg68; reg120=reg113-reg111; reg121=reg91*reg112;
    reg123=reg91*reg104; reg124=reg89*reg100; reg125=reg118-reg117; reg126=reg108-reg107; reg128=reg114-reg115;
    reg129=reg98*reg88; reg126=reg126*reg72; reg125=reg125*reg72; reg128=reg128*reg72; reg131=reg129+reg121;
    reg133=reg123+reg124; reg120=reg120*reg72; reg126=ponderation*reg126; reg125=ponderation*reg125; reg134=reg133*reg72;
    reg135=reg131*reg72; reg120=ponderation*reg120; reg128=ponderation*reg128; reg136=ponderation*reg135; T vec_0=reg136;
    T vec_3=-reg120; T vec_5=-reg128; reg120=ponderation*reg134; T vec_1=reg120; T vec_4=-reg125;
    T vec_2=-reg126;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; reg1=1.0/reg1; T reg2=2*(*f.m).shear_modulus_12;
    T reg3=reg1*reg0; reg2=1.0/reg2; T reg4=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg5=1.0/(*f.m).elastic_modulus_3; T reg6=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg7=reg2*reg3; T reg8=pow((*f.m).v1[0],2); T reg9=pow((*f.m).v2[1],2); T reg10=pow((*f.m).v2[0],2); T reg11=pow((*f.m).v1[1],2);
    T reg12=pow((*f.m).v2[2],2); reg9=reg10+reg9; reg10=1.0/(*f.m).elastic_modulus_2; T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg14=pow((*f.m).v1[2],2);
    reg11=reg8+reg11; reg8=reg5*reg7; T reg15=reg6*reg7; T reg16=reg4*reg7; reg14=reg11+reg14;
    reg11=reg10*reg8; T reg17=reg6*reg16; T reg18=reg6*reg15; T reg19=reg13*reg8; reg12=reg9+reg12;
    reg17=reg17+reg19; reg9=reg13*reg15; T reg20=reg10*reg16; T reg21=1.0/(*f.m).elastic_modulus_1; reg18=reg11-reg18;
    reg14=pow(reg14,0.5); reg12=pow(reg12,0.5); reg11=reg17*reg13; T reg22=reg9+reg20; T reg23=reg21*reg18;
    T reg24=(*f.m).v1[1]/reg14; T reg25=(*f.m).v1[2]/reg14; T reg26=(*f.m).v2[1]/reg12; T reg27=(*f.m).v2[2]/reg12; reg12=(*f.m).v2[0]/reg12;
    reg14=(*f.m).v1[0]/reg14; T reg28=reg4*reg15; T reg29=reg25*reg26; T reg30=reg24*reg27; T reg31=reg10*reg7;
    reg8=reg21*reg8; T reg32=reg4*reg16; reg7=reg13*reg7; T reg33=reg5*reg3; T reg34=reg22*reg4;
    T reg35=reg6*reg3; T reg36=reg2*reg0; reg3=reg4*reg3; reg11=reg23-reg11; reg19=reg28+reg19;
    reg23=reg3*reg6; reg28=2*reg14; T reg37=reg2*reg1; T reg38=reg31*reg4; T reg39=reg25*reg12;
    T reg40=reg14*reg27; T reg41=reg35*reg6; reg32=reg8-reg32; reg15=reg21*reg15; reg8=reg30-reg29;
    T reg42=reg36*reg6; T reg43=reg33*reg13; reg34=reg11-reg34; reg33=reg33*reg10; reg11=reg36*reg4;
    reg36=reg36*reg5; T reg44=reg7*reg4; T reg45=2*reg12; reg16=reg13*reg16; T reg46=pow(reg14,2);
    T reg47=pow(reg24,2); reg23=reg43+reg23; reg19=reg19/reg34; reg43=reg45*reg26; T reg48=pow(reg26,2);
    T reg49=pow(reg12,2); reg18=reg18/reg34; T reg50=2*reg8; T reg51=reg37*reg5; reg16=reg15+reg16;
    reg31=reg31*reg21; T reg52=reg11*reg6; reg32=reg32/reg34; reg44=reg15+reg44; reg17=reg17/reg34;
    reg38=reg9+reg38; reg15=reg28*reg24; T reg53=reg42*reg6; T reg54=reg36*reg13; reg36=reg36*reg10;
    T reg55=reg37*reg4; reg3=reg3*reg10; reg37=reg37*reg6; T reg56=reg39-reg40; reg35=reg35*reg13;
    T reg57=reg14*reg26; T reg58=reg24*reg12; reg7=reg7*reg13; reg41=reg33-reg41; reg33=reg43*reg32;
    T reg59=reg15*reg19; T reg60=reg51*reg10; T reg61=reg46*reg18; reg51=reg51*reg13; reg52=reg54+reg52;
    reg54=reg37*reg6; T reg62=reg50*reg56; T reg63=pow(reg56,2); T reg64=pow(reg8,2); T reg65=reg55*reg6;
    T reg66=reg49*reg17; T reg67=pow(reg27,2); reg11=reg11*reg10; reg42=reg42*reg13; T reg68=reg47*reg18;
    T reg69=reg48*reg17; reg44=reg44/reg34; reg7=reg31-reg7; reg16=reg16/reg34; reg31=reg15*reg18;
    T reg70=reg43*reg17; T reg71=reg47*reg19; T reg72=reg48*reg32; T reg73=reg49*reg32; T reg74=reg46*reg19;
    reg53=reg36-reg53; reg23=reg23*reg13; reg36=pow(reg25,2); reg38=reg38/reg34; reg22=reg22/reg34;
    reg41=reg41*reg21; T reg75=reg35+reg3; T reg76=reg57-reg58; T reg77=reg67*reg17; T reg78=reg15*reg38;
    T reg79=reg22*reg63; reg69=reg68+reg69; reg68=reg36*reg18; reg66=reg61+reg66; reg75=reg75*reg4;
    reg61=reg22*reg64; T reg80=reg22*reg62; reg70=reg31+reg70; reg31=reg43*reg44; T reg81=reg48*reg44;
    T reg82=reg47*reg38; reg73=reg74+reg73; reg74=reg16*reg64; reg52=reg52*reg13; T reg83=reg67*reg32;
    T reg84=reg42+reg11; T reg85=reg36*reg19; reg33=reg59+reg33; reg23=reg41-reg23; reg7=reg7/reg34;
    reg41=reg16*reg63; reg59=reg16*reg62; reg72=reg71+reg72; reg71=reg49*reg44; T reg86=reg46*reg38;
    reg53=reg53*reg21; reg65=reg51+reg65; reg37=reg37*reg13; reg55=reg55*reg10; reg54=reg60-reg54;
    reg51=pow(reg76,2); reg65=reg65*reg13; reg75=reg23-reg75; reg52=reg53-reg52; reg23=reg37+reg55;
    reg77=reg68+reg77; reg53=reg22*reg51; reg79=reg69+reg79; reg54=reg54*reg21; reg84=reg84*reg4;
    reg74=reg73+reg74; reg41=reg72+reg41; reg31=reg78+reg31; reg83=reg85+reg83; reg60=reg16*reg51;
    reg68=reg7*reg62; reg59=reg33+reg59; reg33=2*reg26; reg69=reg45*reg27; reg72=2*reg24;
    reg71=reg86+reg71; reg73=reg7*reg64; reg61=reg66+reg61; reg81=reg82+reg81; reg66=reg7*reg63;
    reg78=reg12*reg26; reg82=reg14*reg24; reg85=reg28*reg25; reg86=reg36*reg38; T reg87=reg67*reg44;
    reg80=reg70+reg80; reg70=reg49*reg74; T reg88=reg46*reg61; reg68=reg31+reg68; reg31=reg78*reg59;
    T reg89=reg46*reg79; T reg90=reg14*reg56; reg75=reg75/reg34; reg60=reg83+reg60; reg83=reg85*reg19;
    T reg91=reg32*reg69; reg73=reg71+reg73; reg66=reg81+reg66; reg87=reg86+reg87; reg71=reg7*reg51;
    reg53=reg77+reg53; reg77=reg85*reg18; reg81=reg17*reg69; reg86=reg79*reg82; T reg92=reg78*reg41;
    T reg93=reg8*reg56; T reg94=reg8*reg24; T reg95=reg48*reg59; T reg96=reg47*reg80; reg23=reg23*reg4;
    T reg97=reg72*reg25; reg65=reg54-reg65; reg54=reg33*reg27; T reg98=reg49*reg41; T reg99=reg80*reg82;
    reg41=reg48*reg41; reg79=reg47*reg79; T reg100=reg48*reg74; T reg101=reg47*reg61; reg59=reg49*reg59;
    reg80=reg46*reg80; reg57=reg58+reg57; reg84=reg52-reg84; reg52=reg24*reg26; reg58=reg14*reg12;
    T reg102=2*reg56; reg50=reg50*reg76; T reg103=reg12*reg56; T reg104=reg8*reg26; reg90=reg94+reg90;
    reg94=reg28*reg12; reg71=reg87+reg71; reg32=reg32*reg54; reg19=reg97*reg19; reg87=reg16*reg50;
    reg91=reg83+reg91; reg17=reg17*reg54; reg18=reg97*reg18; reg83=reg22*reg50; reg81=reg77+reg81;
    reg77=reg85*reg38; T reg105=reg93*reg68; reg31=reg99+reg31; reg99=reg47*reg53; T reg106=reg63*reg66;
    reg41=reg79+reg41; reg79=reg63*reg73; reg100=reg101+reg100; reg101=reg64*reg68; reg59=reg80+reg59;
    reg80=reg8*reg14; T reg107=reg25*reg27; reg102=reg102*reg76; T reg108=reg44*reg69; T reg109=reg24*reg56;
    T reg110=2*reg25; T reg111=reg72*reg26; T reg112=reg93*reg66; reg92=reg86+reg92; reg74=reg78*reg74;
    reg61=reg82*reg61; reg23=reg65-reg23; reg84=reg84/reg34; reg70=reg88+reg70; reg65=reg64*reg73;
    reg86=reg48*reg60; reg98=reg89+reg98; reg66=reg64*reg66; reg88=reg46*reg53; reg95=reg96+reg95;
    reg89=reg75*reg57; reg96=reg75*reg52; T reg113=reg75*reg58; reg68=reg63*reg68; T reg114=reg49*reg60;
    T reg115=reg84*reg90; reg101=reg59+reg101; reg59=reg84*reg109; reg110=reg110*reg27; reg112=reg92+reg112;
    reg92=reg84*reg80; T reg116=reg96*reg57; T reg117=reg89*reg94; reg34=reg23/reg34; reg72=reg72*reg56;
    reg53=reg53*reg82; reg23=reg25*reg76; T reg118=reg49*reg13; T reg119=reg48*reg10; T reg120=reg111*reg13;
    T reg121=reg21*reg94; T reg122=reg13*reg94; T reg123=reg48*reg13; T reg124=reg111*reg10; T reg125=reg49*reg21;
    T reg126=reg89*reg57; reg105=reg31+reg105; reg31=reg8*reg12; T reg127=reg26*reg56; reg99=reg86+reg99;
    reg86=reg63*reg71; T reg128=reg113*reg94; reg65=reg70+reg65; reg60=reg78*reg60; reg70=reg7*reg50;
    reg108=reg77+reg108; reg87=reg91+reg87; reg77=reg96*reg111; reg106=reg41+reg106; reg32=reg19+reg32;
    reg16=reg16*reg102; reg73=reg93*reg73; reg74=reg61+reg74; reg103=reg104+reg103; reg19=reg113*reg111;
    reg79=reg100+reg79; reg68=reg95+reg68; reg89=reg89*reg111; reg28=reg8*reg28; reg114=reg88+reg114;
    reg41=reg64*reg71; reg66=reg98+reg66; reg22=reg22*reg102; reg44=reg44*reg54; reg17=reg18+reg17;
    reg38=reg97*reg38; reg83=reg81+reg83; reg96=reg96*reg94; reg18=reg75*reg107; reg126=reg105+reg126;
    reg61=reg115*reg90; reg81=reg48*reg6; reg88=reg49*reg4; reg122=reg124-reg122; reg19=reg79+reg19;
    reg79=reg49*reg87; reg91=reg46*reg83; reg95=reg110*reg6; reg118=reg119-reg118; reg39=reg40+reg39;
    reg40=reg67*reg6; reg117=reg101+reg117; reg98=reg115*reg28; reg100=reg59*reg90; reg96=reg66+reg96;
    reg66=reg59*reg72; reg60=reg53+reg60; reg71=reg93*reg71; reg77=reg106+reg77; reg116=reg112+reg116;
    reg53=reg18*reg94; reg101=(*f.m).alpha_1*reg46; reg104=reg111*reg6; reg41=reg114+reg41; reg59=reg59*reg28;
    reg105=reg92*reg72; reg113=reg113*reg57; reg106=reg4*reg94; reg73=reg74+reg73; reg74=reg48*reg87;
    reg112=reg92*reg28; reg7=reg7*reg102; reg44=reg38+reg44; reg86=reg99+reg86; reg70=reg108+reg70;
    reg38=(*f.m).alpha_2*reg49; reg128=reg65+reg128; reg89=reg68+reg89; reg16=reg32+reg16; reg32=reg18*reg111;
    reg65=reg47*reg83; reg115=reg115*reg72; reg22=reg17+reg22; reg33=reg33*reg56; reg17=reg34*reg103;
    reg68=reg34*reg31; reg99=reg84*reg23; reg108=reg27*reg76; reg45=reg8*reg45; reg120=reg121-reg120;
    reg114=reg110*reg4; reg119=(*f.m).alpha_1*reg47; reg123=reg125-reg123; reg121=reg67*reg4; reg124=reg8*reg25;
    reg125=reg14*reg76; T reg129=reg34*reg127; T reg130=(*f.m).alpha_2*reg48; T reg131=reg48*reg16; reg79=reg91+reg79;
    reg91=reg47*reg22; T reg132=reg70*reg64; reg7=reg44+reg7; reg115=reg89+reg115; reg44=reg75*reg39;
    reg89=reg17*reg45; reg98=reg117+reg98; reg117=reg49*reg16; T reg133=reg46*reg22; T reg134=(*f.m).alpha_1*reg36;
    T reg135=reg12*reg76; T reg136=reg8*reg27; T reg137=(*f.m).alpha_2*reg67; reg130=reg119+reg130; reg119=reg63*(*f.m).alpha_3;
    reg105=reg19+reg105; reg19=reg68*reg33; reg113=reg73+reg113; reg92=reg92*reg90; reg73=reg129*reg45;
    reg59=reg96+reg59; reg96=reg70*reg63; T reg138=reg99*reg72; reg66=reg77+reg66; reg77=reg129*reg33;
    reg74=reg65+reg74; reg121=reg123-reg121; reg40=reg118-reg40; reg95=reg122-reg95; reg81=reg88+reg81;
    reg101=reg38+reg101; reg38=reg67*reg5; reg65=reg68*reg45; reg112=reg128+reg112; reg88=reg64*(*f.m).alpha_3;
    reg125=reg124+reg125; reg118=reg24*reg76; reg122=reg25*reg56; reg61=reg126+reg61; reg123=reg17*reg103;
    reg124=reg47*reg13; reg83=reg82*reg83; reg87=reg78*reg87; reg17=reg17*reg33; reg10=reg47*reg10;
    reg126=reg34*reg108; reg13=reg46*reg13; reg100=reg116+reg100; reg129=reg129*reg103; reg116=reg99*reg28;
    reg32=reg86+reg32; reg53=reg41+reg53; reg21=reg46*reg21; reg71=reg60+reg71; reg114=reg120-reg114;
    reg30=reg29+reg30; reg18=reg18*reg57; reg104=reg106+reg104; reg110=reg110*reg5; reg29=reg46*reg4;
    reg41=reg47*reg6; reg88=reg101+reg88; reg81=reg38-reg81; reg38=reg126*reg33; reg60=reg44*reg111;
    reg117=reg133+reg117; reg92=reg113+reg92; reg68=reg68*reg103; reg96=reg74+reg96; reg118=reg122+reg118;
    reg119=reg130+reg119; reg138=reg32+reg138; reg12=reg12*reg27; reg104=reg110-reg104; reg14=reg14*reg25;
    reg32=reg44*reg94; reg132=reg79+reg132; reg89=reg98+reg89; reg17=reg115+reg17; reg116=reg53+reg116;
    reg53=reg126*reg45; reg74=reg114*reg58; reg79=reg95*reg52; reg129=reg100+reg129; reg18=reg71+reg18;
    reg99=reg99*reg90; reg65=reg112+reg65; reg123=reg61+reg123; reg87=reg83+reg87; reg70=reg70*reg93;
    reg13=reg10-reg13; reg6=reg36*reg6; reg124=reg21-reg124; reg4=reg36*reg4; reg22=reg82*reg22;
    reg16=reg78*reg16; reg10=reg48*reg95; reg21=reg49*reg114; reg61=reg48*reg40; reg71=reg49*reg121;
    reg83=reg84*reg125; reg75=reg75*reg30; reg86=reg51*(*f.m).alpha_3; reg137=reg134+reg137; reg135=reg136+reg135;
    reg98=reg27*reg56; reg100=reg26*reg76; reg101=(*f.m).alpha_2*reg78; reg106=(*f.m).alpha_1*reg82; reg110=reg40*reg52;
    reg112=reg121*reg58; reg113=reg47*reg95; reg77=reg66+reg77; reg66=reg46*reg121; reg115=reg47*reg40;
    reg120=reg46*reg114; reg19=reg105+reg19; reg131=reg91+reg131; reg91=reg7*reg63; reg73=reg59+reg73;
    reg59=reg7*reg64; reg16=reg22+reg16; reg7=reg7*reg93; reg22=reg67*reg104; reg10=reg21+reg10;
    reg21=reg36*reg81; reg4=reg124-reg4; reg113=reg120+reg113; reg105=reg34*reg135; reg120=reg67*reg81;
    reg61=reg71+reg61; reg71=(*f.m).alpha_1*reg14; reg59=reg117+reg59; reg100=reg98+reg100; reg94=reg75*reg94;
    reg38=reg138+reg38; reg101=reg106+reg101; reg93=reg93*(*f.m).alpha_3; reg98=reg119*reg77; reg106=reg88*reg19;
    reg86=reg137+reg86; reg117=reg83*reg28; reg27=reg26*reg27; reg25=reg24*reg25; reg8=reg8*reg76;
    reg24=reg119*reg73; reg26=reg88*reg65; reg32=reg132+reg32; reg60=reg96+reg60; reg84=reg84*reg118;
    reg96=reg36*reg104; reg91=reg131+reg91; reg111=reg75*reg111; reg122=reg107*reg81; reg79=reg74+reg79;
    reg74=reg107*reg104; reg53=reg116+reg53; reg110=reg112+reg110; reg112=reg77*reg123; reg116=reg73*reg123;
    reg124=reg17*reg129; reg128=reg2*reg57; reg78=reg78*reg2; reg99=reg18+reg99; reg18=reg89*reg129;
    reg126=reg126*reg103; reg130=(*f.m).alpha_2*reg12; reg115=reg66+reg115; reg68=reg92+reg68; reg6=reg13-reg6;
    reg44=reg44*reg57; reg70=reg87+reg70; reg5=reg36*reg5; reg41=reg29+reg41; reg13=reg83*reg72;
    reg126=reg99+reg126; reg93=reg101+reg93; reg75=reg75*reg57; reg122=reg110+reg122; reg94=reg59+reg94;
    reg29=reg57*reg78; reg59=reg105*reg45; reg130=reg71+reg130; reg34=reg34*reg100; reg117=reg32+reg117;
    reg32=(*f.m).alpha_1*reg25; reg66=(*f.m).alpha_2*reg27; reg44=reg70+reg44; reg74=reg79+reg74; reg70=reg57*reg128;
    reg83=reg83*reg90; reg71=reg47*reg6; reg8=reg8*(*f.m).alpha_3; reg79=reg46*reg4; reg41=reg5-reg41;
    reg5=reg1*reg39; reg12=reg12*reg1; reg28=reg84*reg28; reg7=reg16+reg7; reg16=reg15*reg128;
    reg111=reg91+reg111; reg96=reg113+reg96; reg87=reg43*reg78; reg120=reg61+reg120; reg13=reg60+reg13;
    reg72=reg84*reg72; reg60=reg105*reg33; reg61=reg89*reg77; reg21=reg115+reg21; reg91=reg73*reg17;
    reg24=reg26+reg24; reg53=reg86*reg53; reg18=reg116-reg18; reg26=reg88*reg68; reg92=reg48*reg6;
    reg99=reg49*reg4; reg101=reg15*reg78; reg56=reg76*reg56; reg124=reg112-reg124; reg98=reg106+reg98;
    reg38=reg86*reg38; reg76=reg119*reg129; reg22=reg10+reg22; reg10=reg43*reg128; reg101=reg21+reg101;
    reg83=reg44+reg83; reg105=reg105*reg103; reg21=reg85*reg12; reg44=reg4*reg58; reg106=reg6*reg52;
    reg33=reg34*reg33; reg61=reg91-reg61; reg60=reg13+reg60; reg13=reg18*reg19; reg91=reg65*reg124;
    reg76=reg26+reg76; reg126=reg86*reg126; reg2=reg82*reg2; reg26=reg39*reg5; reg70=reg74+reg70;
    reg29=reg122+reg29; reg74=reg93*reg17; reg66=reg32+reg66; reg56=reg56*(*f.m).alpha_3; reg36=reg36*reg41;
    reg38=reg98+reg38; reg45=reg34*reg45; reg28=reg94+reg28; reg8=reg130+reg8; reg32=reg93*reg89;
    reg53=reg24+reg53; reg92=reg99+reg92; reg67=reg67*reg41; reg71=reg79+reg71; reg24=reg39*reg12;
    reg72=reg111+reg72; reg79=reg85*reg5; reg87=reg120+reg87; reg82=reg69*reg12; reg59=reg117+reg59;
    reg16=reg96+reg16; reg94=reg0*reg30; reg84=reg84*reg90; reg75=reg7+reg75; reg10=reg22+reg10;
    reg7=reg69*reg5; reg27=reg27*reg0; reg22=reg93*reg123; reg126=reg76+reg126; reg107=reg107*reg41;
    reg106=reg44+reg106; reg21=reg101+reg21; reg44=reg97*reg27; reg60=reg8*reg60; reg79=reg16+reg79;
    reg16=reg43*reg2; reg67=reg92+reg67; reg45=reg28+reg45; reg28=reg97*reg94; reg74=reg38+reg74;
    reg56=reg66+reg56; reg59=reg8*reg59; reg32=reg53+reg32; reg38=reg68*reg89; reg82=reg87+reg82;
    reg53=reg54*reg27; reg34=reg34*reg103; reg84=reg75+reg84; reg7=reg10+reg7; reg105=reg83+reg105;
    reg10=reg54*reg94; reg66=reg65*reg123; reg75=reg17*reg68; reg33=reg72+reg33; reg1=reg14*reg1;
    reg26=reg70+reg26; reg14=reg30*reg94; reg13=reg91-reg13; reg70=reg30*reg27; reg72=reg61*reg68;
    reg24=reg29+reg24; reg29=reg19*reg123; reg36=reg71+reg36; reg71=reg15*reg2; reg76=reg89*reg19;
    reg59=reg32+reg59; reg32=reg65*reg17; reg83=reg73*reg68; reg45=reg56*reg45; reg38=reg66-reg38;
    reg70=reg24+reg70; reg24=reg65*reg129; reg53=reg82+reg53; reg34=reg84+reg34; reg44=reg21+reg44;
    reg14=reg26+reg14; reg72=reg13+reg72; reg10=reg7+reg10; reg7=reg68*reg77; reg75=reg29-reg75;
    reg13=reg19*reg129; reg107=reg106+reg107; reg21=reg57*reg2; reg0=reg25*reg0; reg105=reg8*reg105;
    reg22=reg126+reg22; reg71=reg36+reg71; reg85=reg85*reg1; reg60=reg74+reg60; reg33=reg56*reg33;
    reg69=reg69*reg1; reg16=reg67+reg16; reg28=reg79+reg28; reg21=reg107+reg21; reg83=reg24-reg83;
    reg39=reg39*reg1; reg97=reg97*reg0; reg85=reg71+reg85; reg33=reg60+reg33; reg105=reg22+reg105;
    reg34=reg56*reg34; reg22=reg10*reg70; reg45=reg59+reg45; reg24=reg65*reg77; reg76=reg32-reg76;
    reg25=reg73*reg19; reg18=reg18/reg72; reg7=reg13-reg7; reg75=reg75/reg72; reg124=reg124/reg72;
    reg13=reg44*reg14; reg26=reg53*reg14; reg69=reg16+reg69; reg38=reg38/reg72; reg16=elem.pos(2)[1]-elem.pos(0)[1];
    reg29=reg28*reg70; reg32=elem.pos(2)[0]-elem.pos(0)[0]; reg36=elem.pos(1)[0]-elem.pos(0)[0]; reg59=elem.pos(1)[1]-elem.pos(0)[1]; reg54=reg54*reg0;
    reg60=reg32*reg59; reg34=reg105+reg34; reg38=reg38*reg33; reg75=reg75*reg45; reg18=reg18*reg33;
    reg124=reg124*reg45; reg66=reg10*reg44; reg29=reg13-reg29; reg7=reg7/reg72; reg25=reg24-reg25;
    reg97=reg85+reg97; reg76=reg76/reg72; reg83=reg83/reg72; reg69=reg54+reg69; reg22=reg26-reg22;
    reg61=reg61/reg72; reg30=reg30*reg0; reg39=reg21+reg39; reg13=reg16*reg36; reg21=reg53*reg28;
    reg30=reg39+reg30; reg60=reg13-reg60; reg33=reg83*reg33; reg72=reg25/reg72; reg21=reg66-reg21;
    reg13=reg29*reg69; reg24=reg97*reg22; reg75=reg38-reg75; reg76=reg76*reg34; reg45=reg7*reg45;
    reg18=reg124-reg18; reg61=reg61*reg34; reg18=reg61+reg18; reg36=reg36/reg60; reg59=reg59/reg60;
    reg33=reg45-reg33; reg7=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; reg34=reg72*reg34; reg25=1-(*f.m).resolution; reg76=reg75-reg76;
    reg26=reg53*reg30; reg38=reg69*reg70; reg16=reg16/reg60; reg32=reg32/reg60; reg39=reg21*reg30;
    reg45=reg97*reg70; reg54=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; reg61=reg44*reg30; reg13=reg24-reg13; reg24=vectors[0][indices[2]+1]-vectors[0][indices[0]+1];
    reg66=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg67=reg66*reg16; reg71=reg24*reg59; reg72=reg97*reg14; reg74=reg119*(*f.m).resolution;
    reg26=reg38-reg26; reg18=reg18*reg25; reg38=reg10*reg30; reg76=reg76*reg25; reg75=reg69*reg14;
    reg79=reg54*reg36; reg82=reg7*reg32; reg83=reg88*(*f.m).resolution; reg84=reg28*reg30; reg39=reg13+reg39;
    reg13=reg53*reg97; reg85=reg69*reg44; reg61=reg45-reg61; reg33=reg34+reg33; reg24=reg24*reg36;
    reg38=reg75-reg38; reg54=reg54*reg59; reg66=reg66*reg32; reg34=reg69*reg28; reg33=reg33*reg25;
    reg18=reg83+reg18; reg45=reg10*reg97; reg76=reg74+reg76; reg85=reg13-reg85; reg61=reg61/reg39;
    reg84=reg72-reg84; reg82=reg79-reg82; reg71=reg67-reg71; reg13=reg93*(*f.m).resolution; reg7=reg7*reg16;
    reg26=reg26/reg39; reg84=reg84/reg39; reg22=reg22/reg39; reg38=reg38/reg39; reg29=reg29/reg39;
    reg82=reg71+reg82; reg68=reg68*reg25; reg129=reg25*reg129; reg34=reg45-reg34; reg66=reg24-reg66;
    reg85=reg85/reg39; reg54=reg7-reg54; reg7=reg26*(*f.m).resolution; reg24=reg61*(*f.m).resolution; reg33=reg13+reg33;
    reg18=(*f.m).deltaT*reg18; reg76=(*f.m).deltaT*reg76; reg73=reg73*reg25; reg13=reg84*(*f.m).resolution; reg19=reg25*reg19;
    reg77=reg25*reg77; reg24=reg129-reg24; reg7=reg68+reg7; reg45=reg85*(*f.m).resolution; reg34=reg34/reg39;
    reg39=reg21/reg39; reg82=0.5*reg82; reg65=reg65*reg25; reg21=reg38*(*f.m).resolution; reg67=reg29*(*f.m).resolution;
    reg33=(*f.m).deltaT*reg33; reg123=reg25*reg123; reg54=reg54-reg18; reg66=reg66-reg76; reg68=reg22*(*f.m).resolution;
    reg71=reg66*reg24; reg72=reg54*reg7; reg77=reg13+reg77; reg21=reg19-reg21; reg67=reg73-reg67;
    reg65=reg68+reg65; reg13=reg34*(*f.m).resolution; reg17=reg17*reg25; reg82=reg82-reg33; reg19=reg39*(*f.m).resolution;
    reg89=reg25*reg89; reg45=reg123+reg45; reg25=reg54*reg65; reg68=reg66*reg77; reg89=reg19+reg89;
    reg54=reg54*reg21; reg13=reg17-reg13; reg17=reg82*reg45; reg71=reg72+reg71; reg66=reg66*reg67;
    reg17=reg71+reg17; reg19=reg32-reg36; reg25=reg66+reg25; reg66=reg82*reg89; reg71=reg59-reg16;
    reg82=reg82*reg13; reg54=reg68+reg54; reg82=reg54+reg82; reg54=0.5*reg16; reg68=0.5*reg59;
    reg72=0.5*reg32; reg73=0.5*reg71; reg74=0.5*reg36; reg66=reg25+reg66; reg17=2*reg17;
    reg25=0.5*reg19; reg75=reg17*reg68; reg79=reg82*reg36; reg83=reg17*reg25; reg87=reg71*reg66;
    reg91=reg17*reg73; reg92=reg19*reg82; reg96=reg66*reg59; reg98=reg17*reg74; reg99=reg17*reg72;
    reg101=reg66*reg16; reg105=reg82*reg32; reg106=reg17*reg54; reg79=reg79-reg75; reg98=reg98-reg96;
    reg106=reg106-reg105; reg101=reg101-reg99; reg92=reg91+reg92; reg83=reg87+reg83; reg101=reg101*reg60;
    reg106=reg106*reg60; reg92=reg92*reg60; reg98=reg98*reg60; reg79=reg79*reg60; reg83=reg83*reg60;
    T vec_2=ponderation*reg101; T vec_3=ponderation*reg106; T vec_1=ponderation*reg92; T vec_5=ponderation*reg79; T vec_0=ponderation*reg83;
    T vec_4=ponderation*reg98;
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
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=pow((*f.m).v2[1],2); T reg5=pow((*f.m).v2[0],2); T reg6=1.0/(*f.m).elastic_modulus_3;
    T reg7=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg8=reg2*reg3; T reg9=pow((*f.m).v1[0],2); T reg10=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg11=pow((*f.m).v1[1],2);
    reg4=reg5+reg4; reg5=reg10*reg8; reg11=reg9+reg11; reg9=reg7*reg8; T reg12=pow((*f.m).v1[2],2);
    T reg13=reg6*reg8; T reg14=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=pow((*f.m).v2[2],2); reg16=reg4+reg16;
    reg4=reg7*reg5; reg12=reg11+reg12; reg11=reg14*reg13; T reg17=reg7*reg9; T reg18=reg15*reg13;
    reg12=pow(reg12,0.5); T reg19=reg15*reg5; T reg20=reg14*reg9; reg4=reg11+reg4; reg17=reg18-reg17;
    reg18=1.0/(*f.m).elastic_modulus_1; reg16=pow(reg16,0.5); T reg21=(*f.m).v1[1]/reg12; T reg22=(*f.m).v1[2]/reg12; T reg23=(*f.m).v2[1]/reg16;
    T reg24=(*f.m).v2[2]/reg16; T reg25=reg18*reg17; T reg26=reg14*reg4; T reg27=reg20+reg19; T reg28=reg2*reg0;
    T reg29=reg7*reg3; T reg30=reg10*reg3; reg3=reg6*reg3; reg13=reg18*reg13; T reg31=reg21*reg24;
    reg26=reg25-reg26; reg25=reg22*reg23; reg12=(*f.m).v1[0]/reg12; reg16=(*f.m).v2[0]/reg16; T reg32=reg14*reg8;
    T reg33=reg10*reg5; T reg34=reg10*reg9; T reg35=reg10*reg27; reg8=reg15*reg8; T reg36=reg10*reg8;
    T reg37=reg7*reg28; reg5=reg14*reg5; T reg38=reg22*reg16; reg33=reg13-reg33; reg13=reg2*reg1;
    T reg39=reg10*reg28; T reg40=reg12*reg24; T reg41=2*reg21; T reg42=2*reg16; T reg43=reg15*reg3;
    T reg44=reg31-reg25; reg35=reg26-reg35; reg3=reg14*reg3; reg26=reg7*reg29; T reg45=reg7*reg30;
    reg28=reg6*reg28; reg34=reg11+reg34; reg11=reg10*reg32; reg9=reg18*reg9; T reg46=2*reg12;
    reg32=reg14*reg32; T reg47=2*reg22; T reg48=pow(reg16,2); T reg49=reg38-reg40; T reg50=reg21*reg46;
    reg29=reg14*reg29; T reg51=pow(reg12,2); T reg52=pow(reg21,2); T reg53=reg7*reg13; reg30=reg15*reg30;
    reg34=reg34/reg35; T reg54=reg7*reg39; reg26=reg43-reg26; reg36=reg20+reg36; reg17=reg17/reg35;
    reg43=reg6*reg13; reg13=reg10*reg13; reg33=reg33/reg35; reg45=reg3+reg45; reg11=reg9+reg11;
    reg3=reg41*reg23; T reg55=reg15*reg28; reg28=reg14*reg28; T reg56=2*reg44; T reg57=reg42*reg23;
    T reg58=reg16*reg46; reg4=reg4/reg35; reg8=reg18*reg8; T reg59=reg7*reg37; reg5=reg9+reg5;
    reg9=pow(reg23,2); T reg60=reg21*reg16; T reg61=reg12*reg23; reg39=reg15*reg39; T reg62=reg14*reg43;
    T reg63=reg58*reg18; reg54=reg28+reg54; reg27=reg27/reg35; reg59=reg55-reg59; reg28=reg9*reg14;
    reg55=reg48*reg18; T reg64=reg7*reg53; T reg65=reg3*reg14; T reg66=reg48*reg14; T reg67=reg9*reg15;
    T reg68=reg48*reg4; T reg69=reg58*reg14; T reg70=reg3*reg15; T reg71=reg50*reg17; T reg72=reg57*reg4;
    T reg73=pow(reg22,2); reg43=reg15*reg43; T reg74=reg51*reg17; T reg75=reg61-reg60; T reg76=pow(reg24,2);
    T reg77=pow(reg44,2); T reg78=pow(reg49,2); T reg79=reg56*reg49; T reg80=reg29+reg30; T reg81=reg51*reg34;
    reg37=reg14*reg37; T reg82=reg48*reg33; reg26=reg18*reg26; T reg83=reg24*reg47; T reg84=reg9*reg4;
    T reg85=reg52*reg17; T reg86=reg52*reg34; reg32=reg8-reg32; reg11=reg11/reg35; reg36=reg36/reg35;
    reg45=reg14*reg45; reg5=reg5/reg35; reg8=reg7*reg13; T reg87=reg9*reg33; T reg88=reg57*reg33;
    T reg89=reg50*reg34; reg87=reg86+reg87; reg86=reg76*reg10; reg28=reg55-reg28; reg54=reg14*reg54;
    reg55=reg77*reg5; reg82=reg81+reg82; reg81=reg78*reg5; T reg90=reg37+reg39; T reg91=reg10*reg83;
    reg65=reg63-reg65; reg63=reg76*reg7; reg66=reg67-reg66; reg84=reg85+reg84; reg53=reg14*reg53;
    reg68=reg74+reg68; reg13=reg15*reg13; reg67=reg83*reg7; reg74=reg77*reg27; reg69=reg70-reg69;
    reg64=reg43-reg64; reg45=reg26-reg45; reg26=reg79*reg27; reg72=reg71+reg72; reg32=reg32/reg35;
    reg80=reg10*reg80; reg8=reg62+reg8; reg43=pow(reg75,2); reg62=reg79*reg5; reg70=reg73*reg17;
    reg71=reg76*reg4; reg88=reg89+reg88; reg85=reg78*reg27; reg89=reg51*reg36; T reg92=reg48*reg11;
    T reg93=reg52*reg36; T reg94=reg9*reg11; T reg95=reg50*reg36; T reg96=reg57*reg11; T reg97=reg76*reg33;
    T reg98=reg3*reg7; T reg99=reg9*reg7; T reg100=reg73*reg34; reg59=reg18*reg59; T reg101=reg10*reg58;
    T reg102=reg48*reg10; T reg103=reg73*reg36; reg68=reg74+reg68; reg74=reg52*reg14; T reg104=2*reg23;
    T reg105=reg78*reg32; reg71=reg70+reg71; reg94=reg93+reg94; reg70=reg43*reg27; reg93=reg42*reg24;
    reg64=reg18*reg64; T reg106=reg53+reg13; reg67=reg69-reg67; reg55=reg82+reg55; reg80=reg45-reg80;
    reg45=reg76*reg6; reg69=reg43*reg5; reg82=reg77*reg32; reg26=reg72+reg26; reg92=reg89+reg92;
    reg72=reg22*reg46; reg90=reg10*reg90; reg89=reg16*reg23; reg86=reg28-reg86; reg28=reg12*reg21;
    reg62=reg88+reg62; reg85=reg84+reg85; reg54=reg59-reg54; reg59=reg83*reg6; reg98=reg101+reg98;
    reg8=reg14*reg8; reg91=reg65-reg91; reg65=reg51*reg14; reg84=reg52*reg15; reg88=reg76*reg11;
    reg101=reg21*reg23; T reg107=reg12*reg16; reg97=reg100+reg97; reg81=reg87+reg81; reg63=reg66-reg63;
    reg66=reg51*reg18; reg96=reg95+reg96; reg87=reg79*reg32; reg99=reg102+reg99; reg95=reg63*reg101;
    reg100=reg63*reg9; reg106=reg10*reg106; reg99=reg45-reg99; reg45=reg86*reg48; reg102=reg104*reg24;
    T reg108=reg86*reg107; T reg109=reg51*reg10; reg98=reg59-reg98; reg59=reg91*reg107; T reg110=reg67*reg101;
    T reg111=reg62*reg89; T reg112=reg73*reg7; reg65=reg84-reg65; reg84=reg52*reg7; T reg113=reg28*reg26;
    reg74=reg66-reg74; reg66=reg73*reg10; T reg114=reg86*reg51; reg70=reg71+reg70; reg71=reg63*reg52;
    T reg115=reg51*reg26; T reg116=reg62*reg48; reg8=reg64-reg8; reg64=reg91*reg51; T reg117=reg67*reg52;
    T reg118=reg41*reg22; T reg119=reg55*reg48; T reg120=reg91*reg48; T reg121=reg93*reg33; T reg122=reg72*reg34;
    T reg123=reg67*reg9; T reg124=reg81*reg9; T reg125=reg52*reg85; T reg126=reg72*reg17; reg80=reg80/reg35;
    reg69=reg97+reg69; reg97=reg43*reg32; T reg127=reg62*reg9; T reg128=reg93*reg4; reg88=reg103+reg88;
    reg103=reg52*reg26; reg105=reg94+reg105; reg56=reg56*reg75; reg94=2*reg49; T reg129=reg51*reg85;
    T reg130=reg81*reg48; T reg131=reg55*reg9; T reg132=reg22*reg24; reg61=reg60+reg61; reg60=reg52*reg68;
    reg82=reg92+reg82; reg90=reg54-reg90; reg54=reg49*reg44; reg92=reg28*reg85; reg87=reg96+reg87;
    reg96=reg21*reg44; T reg133=reg12*reg49; T reg134=reg81*reg89; T reg135=reg51*reg68; reg33=reg102*reg33;
    reg34=reg118*reg34; reg97=reg88+reg97; reg66=reg74-reg66; reg74=reg61*reg2; reg88=reg89*reg2;
    T reg136=reg72*reg36; T reg137=reg93*reg11; T reg138=reg28*reg68; T reg139=reg55*reg89; reg90=reg90/reg35;
    T reg140=reg82*reg78; reg131=reg60+reg131; reg60=reg23*reg44; T reg141=reg16*reg49; T reg142=reg98*reg73;
    reg117=reg64+reg117; reg64=reg16*reg24; reg106=reg8-reg106; reg4=reg102*reg4; reg17=reg118*reg17;
    reg8=reg56*reg27; reg128=reg126+reg128; reg126=reg99*reg73; reg71=reg114+reg71; reg127=reg103+reg127;
    reg103=reg87*reg78; reg121=reg122+reg121; reg114=reg56*reg5; reg122=reg99*reg132; reg95=reg108+reg95;
    reg108=reg105*reg77; reg112=reg65-reg112; reg65=reg52*reg70; T reg143=reg51*reg70; T reg144=reg69*reg48;
    T reg145=reg98*reg76; reg123=reg120+reg123; reg116=reg115+reg116; reg115=reg87*reg77; reg120=reg105*reg78;
    reg124=reg125+reg124; reg133=reg96+reg133; reg96=0.5*elem.pos(0)[1]; reg125=0.5*elem.pos(1)[1]; T reg146=reg21*reg49;
    T reg147=reg12*reg44; reg38=reg40+reg38; reg94=reg94*reg75; reg40=0.5*elem.pos(1)[0]; T reg148=reg99*reg76;
    reg100=reg45+reg100; reg45=0.5*elem.pos(0)[0]; T reg149=reg107*reg80; T reg150=reg101*reg80; T reg151=reg61*reg80;
    reg84=reg109+reg84; reg109=reg73*reg6; reg119=reg135+reg119; reg135=reg82*reg77; reg134=reg92+reg134;
    reg110=reg59+reg110; reg59=reg98*reg132; reg92=reg105*reg54; reg130=reg129+reg130; reg129=reg69*reg9;
    reg111=reg113+reg111; reg113=reg87*reg54; T reg152=reg82*reg54; reg41=reg41*reg49; reg139=reg138+reg139;
    reg122=reg95+reg122; reg95=reg150*reg61; reg31=reg25+reg31; reg25=reg28*reg70; reg138=reg88*reg61;
    reg103=reg127+reg103; reg127=reg22*reg75; T reg153=reg69*reg89; T reg154=reg97*reg78; T reg155=reg23*reg49;
    T reg156=reg16*reg44; reg92=reg134+reg92; reg35=reg106/reg35; reg129=reg65+reg129; reg120=reg124+reg120;
    reg65=reg66*reg48; reg59=reg110+reg59; reg106=reg74*reg61; reg137=reg136+reg137; reg110=reg56*reg32;
    reg124=reg125-reg96; reg36=reg118*reg36; reg134=reg40+reg45; reg11=reg102*reg11; reg144=reg143+reg144;
    reg136=reg97*reg77; reg143=reg150*reg58; reg115=reg116+reg115; reg116=reg112*reg9; reg108=reg130+reg108;
    reg130=reg132*reg80; T reg157=reg147*reg90; T reg158=0.5*elem.pos(2)[1]; T reg159=reg151*reg58; T reg160=reg146*reg90;
    T reg161=reg133*reg90; reg125=reg96+reg125; reg46=reg44*reg46; reg96=reg149*reg58; reg45=reg40-reg45;
    reg40=0.5*elem.pos(2)[0]; reg113=reg111+reg113; reg111=reg88*reg57; T reg162=reg151*reg61; reg135=reg119+reg135;
    reg148=reg100+reg148; reg100=reg74*reg57; reg119=reg88*reg50; reg126=reg71+reg126; reg71=reg151*reg3;
    reg141=reg60+reg141; reg142=reg117+reg142; reg60=reg74*reg50; reg117=reg112*reg52; reg140=reg131+reg140;
    reg131=reg66*reg51; reg8=reg128+reg8; reg128=reg23*reg24; reg4=reg17+reg4; reg17=reg149*reg3;
    reg27=reg94*reg27; T reg163=reg38*reg1; T reg164=reg64*reg1; reg145=reg123+reg145; reg114=reg121+reg114;
    reg84=reg109-reg84; reg33=reg34+reg33; reg5=reg94*reg5; reg34=reg150*reg3; reg109=reg163*reg93;
    reg153=reg25+reg153; reg96=reg135+reg96; reg25=reg157*reg46; reg121=reg160*reg46; reg162=reg113+reg162;
    reg113=reg161*reg133; reg95=reg92+reg95; reg92=reg160*reg133; reg123=reg112*reg101; reg143=reg108+reg143;
    reg108=reg66*reg107; reg135=reg157*reg41; reg17=reg140+reg17; reg104=reg104*reg49; reg140=reg12*reg22;
    reg27=reg4+reg27; reg5=reg33+reg5; reg110=reg137+reg110; reg152=reg139+reg152; reg4=reg149*reg61;
    reg11=reg36+reg11; reg32=reg94*reg32; reg33=(*f.m).alpha_2*reg9; reg36=(*f.m).alpha_1*reg52; reg137=reg127*reg90;
    reg139=(*f.m).alpha_2*reg48; T reg165=reg156*reg35; T reg166=reg155*reg35; T reg167=reg141*reg35; T reg168=reg130*reg3;
    reg154=reg129+reg154; reg129=reg161*reg46; reg159=reg115+reg159; reg115=reg160*reg41; T reg169=reg84*reg76;
    reg116=reg65+reg116; reg65=reg97*reg54; T reg170=reg130*reg58; reg136=reg144+reg136; reg134=reg40-reg134;
    reg124=reg124+reg158; reg144=reg163*reg38; reg106=reg59+reg106; reg59=0.5*elem.pos(3)[0]; reg34=reg120+reg34;
    reg120=reg28*reg2; reg71=reg103+reg71; reg103=reg128*reg0; T reg171=reg31*reg0; T reg172=(*f.m).alpha_1*reg51;
    reg117=reg131+reg117; reg131=reg84*reg73; reg119=reg126+reg119; reg126=reg164*reg72; reg138=reg122+reg138;
    reg122=reg164*reg38; reg60=reg142+reg60; reg142=reg163*reg72; T reg173=reg164*reg93; T reg174=reg114*reg48;
    T reg175=reg51*reg8; T reg176=reg22*reg44; reg12=reg12*reg75; T reg177=reg114*reg9; reg111=reg148+reg111;
    reg148=reg52*reg8; reg42=reg42*reg44; T reg178=reg161*reg41; reg40=reg45+reg40; reg100=reg145+reg100;
    reg45=0.5*elem.pos(3)[1]; reg125=reg158-reg125; reg145=reg24*reg75; reg158=reg5*reg9; T reg179=reg167*reg141;
    T reg180=(*f.m).alpha_2*reg76; T reg181=(*f.m).alpha_1*reg73; T reg182=(*f.m).alpha_3*reg78; T reg183=reg84*reg132; reg33=reg36+reg33;
    reg36=reg27*reg52; reg123=reg108+reg123; reg108=reg114*reg89; reg113=reg162+reg113; reg162=reg110*reg78;
    reg40=reg40-reg59; reg177=reg148+reg177; reg124=reg124-reg45; reg134=reg59+reg134; reg59=reg167*reg104;
    reg178=reg71+reg178; reg45=reg125+reg45; reg71=(*f.m).alpha_3*reg77; reg125=reg28*reg8; reg139=reg172+reg139;
    reg148=reg22*reg49; reg12=reg176+reg12; reg115=reg34+reg115; reg65=reg153+reg65; reg109=reg100+reg109;
    reg34=reg130*reg61; reg122=reg138+reg122; reg100=reg103*reg31; reg138=reg103*reg102; reg173=reg111+reg173;
    reg111=reg171*reg118; reg142=reg60+reg142; reg60=reg110*reg77; reg153=reg103*reg118; reg126=reg119+reg126;
    reg174=reg175+reg174; reg169=reg116+reg169; reg116=reg120*reg57; reg119=reg137*reg46; reg170=reg136+reg170;
    reg136=reg166*reg42; reg121=reg143+reg121; reg143=0.5*vectors[0][indices[1]+1]; reg172=reg171*reg31; reg144=reg106+reg144;
    reg106=0.5*vectors[0][indices[0]+1]; reg129=reg159+reg129; reg159=reg140*reg1; reg175=reg167*reg42; reg131=reg117+reg131;
    reg117=reg120*reg50; reg176=0.5*vectors[0][indices[0]+0]; T reg184=0.5*vectors[0][indices[1]+0]; T reg185=reg165*reg104; reg135=reg17+reg135;
    reg17=reg24*reg44; reg16=reg16*reg75; reg22=reg21*reg22; reg4=reg152+reg4; reg152=reg157*reg133;
    T reg186=reg5*reg48; T reg187=reg27*reg51; reg32=reg11+reg32; reg11=reg38*reg80; T reg188=reg145*reg35;
    reg25=reg96+reg25; reg96=reg165*reg42; T reg189=reg166*reg141; reg92=reg95+reg92; reg95=reg166*reg104;
    reg21=reg21*reg75; T reg190=reg137*reg41; T reg191=reg171*reg102; reg168=reg154+reg168; reg179=reg113+reg179;
    reg113=reg159*reg93; reg154=reg143-reg106; reg143=reg106+reg143; reg191=reg109+reg191; reg106=(*f.m).alpha_2*reg89;
    reg109=(*f.m).alpha_1*reg28; T reg192=(*f.m).alpha_3*reg43; reg180=reg181+reg180; reg116=reg169+reg116; reg169=reg32*reg77;
    reg181=0.5*vectors[0][indices[2]+1]; reg186=reg187+reg186; reg182=reg33+reg182; reg96=reg25+reg96; reg60=reg174+reg60;
    reg25=reg11*reg58; reg33=reg12*reg90; reg80=reg31*reg80; reg71=reg139+reg71; reg138=reg173+reg138;
    reg175=reg129+reg175; reg129=0.5*vectors[0][indices[2]+0]; reg139=reg11*reg3; reg162=reg177+reg162; reg158=reg36+reg158;
    reg36=reg24*reg49; reg173=reg32*reg78; reg23=reg23*reg75; reg174=reg184-reg176; reg176=reg184+reg176;
    reg59=reg178+reg59; reg185=reg135+reg185; reg135=reg110*reg54; reg117=reg131+reg117; reg131=reg159*reg72;
    reg28=reg27*reg28; reg177=reg188*reg104; reg89=reg5*reg89; reg152=reg4+reg152; reg153=reg126+reg153;
    reg111=reg142+reg111; reg4=reg165*reg141; reg126=reg137*reg133; reg34=reg65+reg34; reg21=reg148+reg21;
    reg95=reg115+reg95; reg189=reg92+reg189; reg190=reg168+reg190; reg183=reg123+reg183; reg65=reg188*reg42;
    reg119=reg170+reg119; reg92=reg120*reg61; reg16=reg17+reg16; reg17=reg124*reg134; reg136=reg121+reg136;
    reg172=reg144+reg172; reg115=reg40*reg45; reg108=reg125+reg108; reg100=reg122+reg100; reg121=reg22*reg0;
    reg169=reg186+reg169; reg122=(*f.m).alpha_3*reg54; reg106=reg109+reg106; reg109=reg96*reg71; reg123=reg159*reg38;
    reg173=reg158+reg173; reg125=reg16*reg35; reg44=reg75*reg44; reg177=reg190+reg177; reg142=reg179*reg136;
    reg90=reg21*reg90; reg144=reg189*reg59; reg4=reg152+reg4; reg148=reg182*reg136; reg152=reg182*reg95;
    reg58=reg80*reg58; reg158=(*f.m).alpha_2*reg64; reg168=reg71*reg185; reg192=reg180+reg192; reg170=reg189*reg175;
    reg178=(*f.m).alpha_1*reg140; reg92=reg183+reg92; reg3=reg80*reg3; reg65=reg119+reg65; reg17=reg115-reg17;
    reg113=reg116+reg113; reg115=reg121*reg102; reg154=reg181+reg154; reg116=0.5*vectors[0][indices[3]+1]; reg143=reg181-reg143;
    reg119=reg138*reg172; reg176=reg129-reg176; reg180=0.5*vectors[0][indices[3]+0]; reg174=reg129+reg174; reg135=reg108+reg135;
    reg23=reg36+reg23; reg36=reg11*reg61; reg131=reg117+reg131; reg108=reg179*reg95; reg117=reg153*reg172;
    reg129=reg111*reg100; reg181=reg33*reg46; reg183=reg191*reg100; reg25=reg60+reg25; reg126=reg34+reg126;
    reg34=reg188*reg141; reg60=reg32*reg54; reg89=reg28+reg89; reg28=reg33*reg41; reg139=reg162+reg139;
    reg162=reg121*reg118; reg34=reg126+reg34; reg115=reg113+reg115; reg113=reg90*reg41; reg3=reg173+reg3;
    reg126=reg111*reg138; reg173=reg153*reg191; reg124=reg124/reg17; reg45=reg45/reg17; reg36=reg135+reg36;
    reg135=reg33*reg133; reg184=reg125*reg42; reg181=reg25+reg181; reg60=reg89+reg60; reg25=reg80*reg61;
    reg144=reg108-reg144; reg174=reg174-reg180; reg176=reg180+reg176; reg143=reg143+reg116; reg116=reg154-reg116;
    reg152=reg168+reg152; reg177=reg192*reg177; reg134=reg134/reg17; reg123=reg92+reg123; reg162=reg131+reg162;
    reg89=reg182*reg189; reg92=reg71*reg4; reg183=reg119-reg183; reg108=reg121*reg31; reg129=reg117-reg129;
    reg117=reg90*reg46; reg58=reg169+reg58; reg148=reg109+reg148; reg49=reg75*reg49; reg65=reg192*reg65;
    reg35=reg23*reg35; reg109=(*f.m).alpha_2*reg128; reg119=(*f.m).alpha_1*reg22; reg131=(*f.m).alpha_3*reg44; reg158=reg178+reg158;
    reg122=reg106+reg122; reg106=reg175*reg95; reg154=reg125*reg104; reg28=reg139+reg28; reg139=reg136*reg59;
    reg170=reg142-reg170; reg17=reg40/reg17; reg106=reg139-reg106; reg40=reg96*reg144; reg139=reg170*reg185;
    reg142=reg35*reg104; reg113=reg3+reg113; reg3=reg35*reg42; reg117=reg58+reg117; reg58=0.78867513459481286553*elem.pos(0)[1];
    reg184=reg181+reg184; reg168=0.21132486540518713447*elem.pos(0)[1]; reg169=0.78867513459481286553*elem.pos(1)[1]; reg178=0.78867513459481286553*elem.pos(1)[0]; reg180=0.21132486540518713447*elem.pos(1)[0];
    reg181=0.21132486540518713447*elem.pos(0)[0]; reg186=reg162*reg183; reg187=reg116*reg45; reg190=reg143*reg124; T reg193=reg176*reg17;
    T reg194=reg174*reg134; T reg195=0.78867513459481286553*elem.pos(0)[0]; T reg196=reg129*reg115; reg108=reg123+reg108; reg25=reg60+reg25;
    reg60=reg90*reg133; reg123=reg125*reg141; reg34=reg192*reg34; reg89=reg92+reg89; reg135=reg36+reg135;
    reg36=reg122*reg59; reg177=reg152+reg177; reg154=reg28+reg154; reg28=0.21132486540518713447*elem.pos(1)[1]; reg131=reg158+reg131;
    reg109=reg119+reg109; reg92=(*f.m).alpha_3*reg49; reg119=reg122*reg175; reg65=reg148+reg65; reg126=reg173-reg126;
    reg36=reg177+reg36; reg154=reg131*reg154; reg148=reg168+reg169; reg152=reg126*reg108; reg158=reg180-reg181;
    reg173=reg4*reg175; reg177=0.78867513459481286553*elem.pos(2)[0]; reg180=reg180+reg195; reg92=reg109+reg92; reg194=reg193-reg194;
    reg34=reg89+reg34; reg89=reg122*reg179; reg190=reg187-reg190; reg196=reg186-reg196; reg116=reg116*reg134;
    reg123=reg135+reg123; reg143=reg143*reg17; reg109=reg35*reg141; reg176=reg176*reg124; reg60=reg25+reg60;
    reg174=reg174*reg45; reg3=reg117+reg3; reg25=reg179*reg185; reg117=reg106*reg4; reg135=reg4*reg59;
    reg139=reg40-reg139; reg142=reg113+reg142; reg119=reg65+reg119; reg184=reg131*reg184; reg40=reg115*reg172;
    reg181=reg181+reg178; reg65=0.21132486540518713447*elem.pos(2)[0]; reg113=0.78867513459481286553*elem.pos(2)[1]; reg186=0.21132486540518713447*elem.pos(2)[1]; reg187=reg191*reg108;
    reg193=reg96*reg179; T reg197=reg58+reg28; T reg198=reg111*reg108; reg168=reg28-reg168; reg28=reg162*reg172;
    reg158=reg158+reg177; T reg199=0.21132486540518713447*elem.pos(3)[0]; T reg200=0.78867513459481286553*elem.pos(3)[0]; reg195=reg178-reg195; reg178=0.21132486540518713447*elem.pos(3)[1];
    reg148=reg113-reg148; reg181=reg177-reg181; reg177=(*f.m).deltaT*reg71; reg197=reg186-reg197; reg58=reg169-reg58;
    reg152=reg196+reg152; reg169=reg162*reg100; reg198=reg28-reg198; reg28=reg153*reg108; reg196=reg162*reg191;
    T reg201=reg111*reg115; reg3=reg92*reg3; reg184=reg119+reg184; reg117=reg139+reg117; reg168=reg113+reg168;
    reg113=reg175*reg185; reg119=reg189*reg185; reg135=reg25-reg135; reg25=reg4*reg95; reg139=reg96*reg189;
    reg173=reg193-reg173; reg193=reg4*reg136; T reg202=0.78867513459481286553*elem.pos(3)[1]; T reg203=reg96*reg59; T reg204=(*f.m).deltaT*reg182;
    reg176=reg174-reg176; elem.epsilon[0][0]=reg176; reg109=reg60+reg109; reg116=reg143-reg116; elem.epsilon[0][1]=reg116;
    reg123=reg131*reg123; reg89=reg34+reg89; reg194=reg190+reg194; reg180=reg65-reg180; reg34=reg66*reg77;
    reg60=reg112*reg78; reg143=reg86*reg77; reg174=reg63*reg78; reg142=reg92*reg142; reg154=reg36+reg154;
    reg36=reg115*reg100; reg187=reg40-reg187; reg40=reg138*reg108; reg190=0.78867513459481286553*PNODE(0).dep[0]; T reg205=0.78867513459481286553*PNODE(0).dep[1];
    T reg206=reg116-reg204; reg60=reg34+reg60; reg173=reg173/reg117; reg34=0.21132486540518713447*PNODE(0).dep[1]; reg170=reg170/reg117;
    T reg207=0.78867513459481286553*PNODE(1).dep[1]; reg25=reg119-reg25; reg119=reg84*reg43; reg135=reg135/reg117; reg144=reg144/reg117;
    T reg208=0.21132486540518713447*PNODE(1).dep[1]; T reg209=reg153*reg115; reg201=reg196-reg201; reg196=reg162*reg138; reg28=reg169-reg28;
    reg197=reg197+reg202; reg198=reg198/reg152; reg129=reg129/reg152; reg40=reg36-reg40; reg58=reg186+reg58;
    reg195=reg65+reg195; reg174=reg143+reg174; reg36=reg99*reg43; reg187=reg187/reg152; reg183=reg183/reg152;
    reg65=(*f.m).deltaT*reg122; reg143=reg67*reg78; reg169=reg91*reg77; reg194=0.5*reg194; elem.epsilon[0][2]=reg194;
    reg180=reg200+reg180; reg181=reg181+reg199; reg202=reg168-reg202; reg109=reg92*reg109; reg123=reg89+reg123;
    reg89=reg176-reg177; reg142=reg154+reg142; reg3=reg184+reg3; reg154=0.21132486540518713447*PNODE(1).dep[0]; reg168=reg136*reg185;
    reg184=0.21132486540518713447*PNODE(0).dep[0]; reg113=reg203-reg113; reg200=reg158-reg200; reg158=reg96*reg95; reg186=0.78867513459481286553*PNODE(1).dep[0];
    reg193=reg139-reg193; reg148=reg148+reg178; reg139=reg148*reg200; reg178=reg58-reg178; reg58=reg198*reg206;
    reg199=reg195-reg199; reg195=reg197*reg200; reg203=reg208-reg34; reg34=reg34+reg207; T reg210=0.78867513459481286553*PNODE(2).dep[1];
    T reg211=reg184+reg186; reg184=reg154-reg184; T reg212=0.78867513459481286553*PNODE(2).dep[0]; T reg213=reg202*reg181; reg173=reg173*reg142;
    reg135=reg135*reg3; T reg214=reg129*reg206; T reg215=reg98*reg43; reg143=reg169+reg143; reg170=reg170*reg142;
    reg144=reg144*reg3; reg109=reg123+reg109; reg123=reg89*reg183; reg168=reg158-reg168; reg158=0.21132486540518713447*PNODE(2).dep[0];
    reg154=reg154+reg190; reg113=reg113/reg117; reg106=reg106/reg117; reg193=reg193/reg117; reg169=0.21132486540518713447*PNODE(2).dep[1];
    T reg216=reg89*reg187; T reg217=reg180*reg202; T reg218=reg194-reg65; reg119=reg60+reg119; reg60=reg120*reg79;
    reg36=reg174+reg36; reg174=reg88*reg79; reg40=reg40/reg152; reg28=reg28/reg152; reg126=reg126/reg152;
    reg201=reg201/reg152; reg209=reg196-reg209; reg25=reg25/reg117; reg208=reg208+reg205; reg196=reg218*reg126;
    reg213=reg139-reg213; reg214=reg123-reg214; reg34=reg210-reg34; reg123=0.21132486540518713447*PNODE(3).dep[1]; reg203=reg210+reg203;
    reg139=0.78867513459481286553*PNODE(3).dep[1]; reg210=reg218*reg201; reg205=reg207-reg205; reg216=reg58-reg216; reg174=reg36+reg174;
    reg36=reg74*reg79; reg215=reg143+reg215; reg106=reg106*reg109; reg58=reg164*reg56; reg143=reg199*reg148;
    reg207=reg178*reg181; reg117=reg168/reg117; reg152=reg209/reg152; reg170=reg144-reg170; reg113=reg113*reg109;
    reg208=reg169-reg208; reg154=reg158-reg154; reg135=reg173-reg135; reg3=reg25*reg3; reg190=reg186-reg190;
    reg142=reg193*reg142; reg217=reg195-reg217; reg60=reg119+reg60; reg206=reg28*reg206; reg184=reg212+reg184;
    reg25=0.78867513459481286553*PNODE(3).dep[0]; reg89=reg89*reg40; reg119=0.21132486540518713447*PNODE(3).dep[0]; reg211=reg212-reg211; reg144=reg159*reg56;
    reg168=reg121*reg94; reg144=reg60+reg144; reg60=reg180/reg217; reg173=reg200/reg217; reg58=reg174+reg58;
    reg174=reg103*reg94; reg186=1-(*f.m).resolution; reg36=reg215+reg36; reg193=reg163*reg56; reg205=reg169+reg205;
    reg210=reg216-reg210; reg218=reg218*reg152; reg206=reg89-reg206; reg190=reg158+reg190; reg154=reg25+reg154;
    reg208=reg139+reg208; reg139=reg203-reg139; reg34=reg34+reg123; reg211=reg119+reg211; reg89=reg197*reg199;
    reg158=reg178*reg180; reg25=reg184-reg25; reg142=reg3-reg142; reg109=reg117*reg109; reg3=reg181/reg213;
    reg207=reg143-reg207; reg113=reg135-reg113; reg117=reg202/reg217; reg135=reg197/reg217; reg143=reg148/reg213;
    reg202=reg202/reg213; reg196=reg214+reg196; reg200=reg200/reg213; reg106=reg170+reg106; reg169=reg139*reg143;
    reg170=reg153*reg210; reg184=reg25*reg60; reg195=reg173*reg154; reg203=reg117*reg208; reg209=reg139*reg135;
    reg142=reg109+reg142; reg206=reg218+reg206; reg109=reg196*reg115; reg212=reg210*reg138; reg214=reg66*reg147;
    reg215=reg112*reg146; reg216=reg86*reg147; reg218=reg63*reg146; reg168=reg144+reg168; reg144=reg196*reg162;
    reg181=reg181/reg207; T reg219=reg199/reg207; reg174=reg58+reg174; reg58=reg178/reg207; reg148=reg148/reg207;
    reg71=(*f.m).resolution*reg71; reg182=(*f.m).resolution*reg182; reg106=reg106*reg186; reg113=reg113*reg186; reg158=reg89-reg158;
    reg89=reg25*reg3; T reg220=reg211*reg200; reg193=reg36+reg193; reg36=reg171*reg94; T reg221=reg34*reg202;
    reg119=reg190-reg119; reg123=reg205-reg123; reg190=reg173*reg208; reg170=reg144+reg170; reg144=reg84*reg127;
    reg215=reg214+reg215; reg197=reg197/reg158; reg205=reg211*reg202; reg214=reg117*reg154; T reg222=reg25*reg143;
    reg25=reg25*reg135; reg89=reg220-reg89; reg36=reg193+reg36; reg221=reg169-reg221; reg106=reg71+reg106;
    reg113=reg182+reg113; reg71=reg148*reg123; reg142=reg142*reg186; reg169=reg34*reg58; reg86=reg86*reg156;
    reg182=reg181*reg119; reg63=reg63*reg155; reg193=reg67*reg146; reg220=reg91*reg147; reg112=reg112*reg155;
    T reg223=reg139*reg3; reg212=reg109+reg212; reg109=reg206*reg191; reg66=reg66*reg156; T reg224=reg211*reg219;
    reg180=reg180/reg158; reg199=reg199/reg158; reg178=reg178/reg158; T reg225=reg174*reg210; T reg226=reg196*reg168;
    reg122=(*f.m).resolution*reg122; reg184=reg195-reg184; reg195=reg34*reg200; T reg227=reg111*reg206; reg203=reg209-reg203;
    reg209=reg99*reg127; reg218=reg216+reg218; reg139=reg139*reg60; reg216=reg98*reg127; reg193=reg220+reg193;
    reg220=(*f.m).resolution*reg198; T reg228=(*f.m).resolution*reg40; T reg229=(*f.m).resolution*reg28; reg169=reg71-reg169; reg142=reg122+reg142;
    reg106=reg106*(*f.m).deltaT; reg113=reg113*(*f.m).deltaT; reg91=reg91*reg156; reg67=reg67*reg155; reg71=reg181*reg123;
    reg34=reg34*reg219; reg211=reg211*reg58; reg122=reg148*reg119; reg205=reg222-reg205; reg227=reg170+reg227;
    reg225=reg226+reg225; reg170=reg36*reg206; reg109=reg212+reg109; reg223=reg195-reg223; reg195=reg180*reg119;
    reg212=reg199*reg154; reg222=reg178*reg208; reg226=reg197*reg123; reg184=reg203+reg184; reg203=reg88*reg133;
    reg209=reg218+reg209; reg218=reg120*reg133; reg139=reg190-reg139; reg144=reg215+reg144; reg214=reg25-reg214;
    reg89=reg221+reg89; reg25=(*f.m).resolution*reg129; reg190=(*f.m).resolution*reg183; reg182=reg224-reg182; reg96=reg96*reg186;
    reg136=reg186*reg136; reg112=reg66+reg112; reg84=reg84*reg145; reg4=reg186*reg4; reg95=reg186*reg95;
    reg63=reg86+reg63; reg66=(*f.m).resolution*reg187; reg189=reg186*reg189; reg99=reg99*reg145; reg185=reg186*reg185;
    reg86=reg159*reg12; reg223=reg223-reg113; reg179=reg186*reg179; reg227=reg177+reg227; reg203=reg209+reg203;
    reg209=reg164*reg12; reg215=(*f.m).resolution*reg201; reg170=reg225+reg170; reg222=reg226-reg222; reg195=reg212-reg195;
    reg192=(*f.m).deltaT*reg192; reg96=reg190+reg96; reg25=reg136-reg25; reg66=reg185-reg66; reg220=reg95+reg220;
    reg228=reg4+reg228; reg229=reg189-reg229; reg182=reg169+reg182; reg4=(*f.m).resolution*reg126; reg175=reg186*reg175;
    reg71=reg34-reg71; reg208=reg199*reg208; reg123=reg180*reg123; reg89=0.5*reg89; reg211=reg122-reg211;
    reg59=reg186*reg59; reg34=(*f.m).resolution*reg152; reg218=reg144+reg218; reg139=reg139-reg113; reg84=reg112+reg84;
    reg120=reg120*reg141; reg184=0.5*reg184; reg214=reg214-reg106; reg99=reg63+reg99; reg88=reg88*reg141;
    reg63=reg74*reg133; reg205=reg205-reg106; reg216=reg193+reg216; reg119=reg197*reg119; reg154=reg178*reg154;
    reg142=reg142*(*f.m).deltaT; reg67=reg91+reg67; reg98=reg98*reg145; reg109=reg204+reg109; reg63=reg216+reg63;
    reg175=reg4+reg175; reg4=reg163*reg12; reg91=reg25*reg139; reg98=reg67+reg98; reg215=reg59-reg215;
    reg120=reg84+reg120; reg159=reg159*reg16; reg59=reg96*reg214; reg88=reg99+reg88; reg164=reg164*reg16;
    reg67=reg205*reg228; reg84=reg223*reg229; reg95=reg223*reg220; reg99=reg205*reg66; reg182=0.5*reg182;
    reg89=reg89-reg142; reg123=reg208-reg123; reg170=reg192+reg170; reg112=reg227+reg109; reg195=reg222+reg195;
    reg154=reg119-reg154; reg119=reg205*reg96; reg122=reg220*reg139; reg136=reg66*reg214; reg211=reg211-reg106;
    reg144=reg223*reg25; reg169=reg103*reg21; reg209=reg203+reg209; reg184=reg184-reg142; reg71=reg71-reg113;
    reg74=reg74*reg141; reg185=reg228*reg214; reg86=reg218+reg86; reg186=reg229*reg139; reg34=reg179+reg34;
    reg179=reg121*reg21; reg189=reg229*reg71; reg190=reg228*reg211; reg193=reg89*reg175; reg195=0.5*reg195;
    reg203=reg66*reg211; reg95=reg99+reg95; reg99=reg89*reg215; reg179=reg86+reg179; reg86=reg215*reg184;
    reg122=reg136+reg122; reg136=reg34*reg184; reg186=reg185+reg186; reg123=reg123-reg113; reg185=reg220*reg71;
    reg208=reg175*reg184; reg169=reg209+reg169; reg91=reg59+reg91; reg159=reg120+reg159; reg121=reg121*reg23;
    reg59=reg25*reg71; reg164=reg88+reg164; reg144=reg119+reg144; reg88=reg171*reg21; reg103=reg103*reg23;
    reg119=reg96*reg211; reg84=reg67+reg84; reg67=reg89*reg34; reg4=reg63+reg4; reg74=reg98+reg74;
    reg63=reg210*reg100; reg98=reg196*reg108; reg163=reg163*reg16; reg154=reg154-reg106; reg182=reg182-reg142;
    reg112=reg170+reg112; reg171=reg171*reg23; reg120=reg96*reg154; reg59=reg119+reg59; reg119=reg66*reg154;
    reg209=reg210*reg169; reg212=reg196*reg179; reg216=reg34*reg182; reg208=reg91+reg208; reg91=reg206*reg172;
    reg63=reg98+reg63; reg86=reg122+reg86; reg112=reg112/3; reg195=reg195-reg142; reg189=reg190+reg189;
    reg98=reg175*reg182; reg122=reg215*reg182; reg67=reg84+reg67; reg193=reg144+reg193; reg84=reg229*reg123;
    reg144=reg228*reg154; reg103=reg164+reg103; reg99=reg95+reg99; reg121=reg159+reg121; reg95=reg25*reg123;
    reg185=reg203+reg185; reg136=reg186+reg136; reg88=reg4+reg88; reg163=reg74+reg163; reg4=reg220*reg123;
    reg99=reg223*reg99; reg216=reg189+reg216; reg171=reg163+reg171; reg193=reg205*reg193; reg95=reg120+reg95;
    reg122=reg185+reg122; reg74=reg175*reg195; reg98=reg59+reg98; reg210=reg210*reg103; reg196=reg196*reg121;
    reg59=reg206*reg88; reg209=reg212+reg209; reg4=reg119+reg4; reg119=reg215*reg195; reg67=2*reg67;
    reg91=reg63+reg91; reg84=reg144+reg84; reg63=reg34*reg195; reg109=reg109-reg112; reg136=2*reg136;
    reg208=reg214*reg208; reg86=reg139*reg86; reg227=reg227-reg112; reg112=reg170-reg112; reg119=reg4+reg119;
    reg63=reg84+reg63; reg99=reg193+reg99; reg91=reg65+reg91; reg216=2*reg216; reg109=pow(reg109,2);
    reg59=reg209+reg59; reg227=pow(reg227,2); reg122=reg71*reg122; reg131=(*f.m).deltaT*reg131; reg67=reg89*reg67;
    reg74=reg95+reg74; reg136=reg184*reg136; reg206=reg171*reg206; reg98=reg211*reg98; reg86=reg208+reg86;
    reg210=reg196+reg210; reg136=reg86+reg136; reg109=reg227+reg109; reg63=2*reg63; reg116=reg116-reg113;
    reg176=reg176-reg106; reg122=reg98+reg122; reg119=reg123*reg119; reg112=pow(reg112,2); reg67=reg99+reg67;
    reg92=(*f.m).deltaT*reg92; reg216=reg182*reg216; reg74=reg154*reg74; reg4=2*reg91; reg206=reg210+reg206;
    reg59=reg131+reg59; reg112=reg109+reg112; reg67=reg67*reg213; reg216=reg122+reg216; reg4=reg91*reg4;
    reg206=reg92+reg206; reg71=2*reg59; reg136=reg217*reg136; reg84=reg176*reg96; reg86=reg116*reg25;
    reg119=reg74+reg119; reg63=reg195*reg63; reg194=reg194-reg142; reg74=reg176*reg66; reg89=reg116*reg220;
    reg71=reg59*reg71; reg116=reg116*reg229; reg4=reg112+reg4; reg176=reg176*reg228; reg59=reg194*reg215;
    reg91=2*reg206; reg89=reg74+reg89; reg216=reg207*reg216; reg67=0.25*reg67; reg86=reg84+reg86;
    reg63=reg119+reg63; reg136=0.25*reg136; reg74=reg194*reg175; reg74=reg86+reg74; elem.sigma[0][0]=reg74;
    reg216=0.25*reg216; reg91=reg206*reg91; reg136=reg67+reg136; reg59=reg89+reg59; elem.sigma[0][1]=reg59;
    reg63=reg158*reg63; reg116=reg176+reg116; reg194=reg194*reg34; reg71=reg4+reg71; reg71=reg91+reg71;
    reg63=0.25*reg63; reg216=reg136+reg216; reg4=reg59*reg101; reg67=reg74*reg107; reg84=reg59*reg9;
    reg86=reg74*reg48; reg59=reg59*reg52; reg74=reg74*reg51; reg194=reg116+reg194; elem.sigma[0][2]=reg194;
    reg84=reg86+reg84; reg86=reg194*reg57; reg71=1.5*reg71; reg89=reg194*reg50; reg59=reg74+reg59;
    reg4=reg67+reg4; reg63=reg216+reg63; reg194=reg194*reg61; elem.sigma_local[0][0]=reg59+reg89; elem.ener=reg63/2;
    elem.sigma_local[0][1]=reg84+reg86; elem.sigma_von_mises=pow(reg71,0.5); elem.sigma_local[0][2]=reg4+reg194;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; T reg2=2*(*f.m).shear_modulus_12; reg1=1.0/reg1; reg0=1.0/reg0;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg5=reg2*reg3; T reg6=1.0/(*f.m).elastic_modulus_3;
    T reg7=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg8=pow((*f.m).v2[1],2); T reg9=pow((*f.m).v2[0],2); T reg10=pow((*f.m).v1[1],2); T reg11=pow((*f.m).v1[0],2);
    T reg12=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg13=1.0/(*f.m).elastic_modulus_2; T reg14=pow((*f.m).v1[2],2); T reg15=reg4*reg5; reg10=reg11+reg10;
    reg11=reg6*reg5; reg8=reg9+reg8; reg9=reg7*reg5; T reg16=pow((*f.m).v2[2],2); T reg17=reg7*reg9;
    T reg18=reg7*reg15; T reg19=reg12*reg11; reg16=reg8+reg16; reg8=reg13*reg11; reg14=reg10+reg14;
    reg10=reg12*reg9; T reg20=reg13*reg15; reg18=reg19+reg18; reg14=pow(reg14,0.5); reg17=reg8-reg17;
    reg16=pow(reg16,0.5); reg8=1.0/(*f.m).elastic_modulus_1; T reg21=(*f.m).v2[2]/reg16; T reg22=reg12*reg18; T reg23=(*f.m).v2[1]/reg16;
    T reg24=reg8*reg17; T reg25=(*f.m).v1[1]/reg14; T reg26=reg10+reg20; T reg27=(*f.m).v1[2]/reg14; T reg28=reg4*reg15;
    T reg29=reg6*reg3; T reg30=reg4*reg9; reg11=reg8*reg11; T reg31=reg13*reg5; reg14=(*f.m).v1[0]/reg14;
    T reg32=reg4*reg26; reg16=(*f.m).v2[0]/reg16; reg22=reg24-reg22; reg24=reg7*reg3; T reg33=reg2*reg0;
    T reg34=reg27*reg23; T reg35=reg25*reg21; reg3=reg4*reg3; reg5=reg12*reg5; T reg36=reg7*reg33;
    reg15=reg12*reg15; T reg37=reg12*reg29; T reg38=reg7*reg24; T reg39=reg2*reg1; reg32=reg22-reg32;
    reg22=2*reg16; T reg40=reg7*reg3; T reg41=reg4*reg33; T reg42=reg35-reg34; T reg43=reg4*reg5;
    T reg44=reg27*reg16; T reg45=reg14*reg21; reg30=reg19+reg30; reg33=reg6*reg33; reg28=reg11-reg28;
    reg9=reg8*reg9; reg29=reg13*reg29; reg11=2*reg14; reg19=reg4*reg31; reg17=reg17/reg32;
    reg15=reg9+reg15; reg5=reg12*reg5; T reg46=2*reg42; reg3=reg13*reg3; T reg47=reg14*reg23;
    T reg48=reg7*reg36; T reg49=reg44-reg45; T reg50=reg12*reg33; reg24=reg12*reg24; reg33=reg13*reg33;
    T reg51=reg25*reg11; reg28=reg28/reg32; reg18=reg18/reg32; reg31=reg8*reg31; T reg52=reg6*reg39;
    reg19=reg10+reg19; reg30=reg30/reg32; T reg53=reg7*reg41; T reg54=reg25*reg16; reg38=reg29-reg38;
    reg29=pow(reg25,2); T reg55=pow(reg14,2); T reg56=pow(reg16,2); T reg57=pow(reg23,2); reg43=reg9+reg43;
    reg9=reg4*reg39; reg39=reg7*reg39; T reg58=reg22*reg23; reg40=reg37+reg40; reg37=pow(reg27,2);
    T reg59=reg47-reg54; reg5=reg31-reg5; reg31=reg55*reg30; T reg60=reg56*reg28; T reg61=reg29*reg30;
    T reg62=reg57*reg28; T reg63=reg51*reg30; T reg64=reg58*reg28; reg15=reg15/reg32; T reg65=reg56*reg18;
    reg43=reg43/reg32; reg36=reg12*reg36; reg19=reg19/reg32; reg41=reg13*reg41; T reg66=reg55*reg17;
    reg53=reg50+reg53; reg50=reg13*reg52; reg52=reg12*reg52; T reg67=reg7*reg39; T reg68=reg51*reg17;
    T reg69=reg58*reg18; T reg70=reg7*reg9; T reg71=pow(reg21,2); reg48=reg33-reg48; reg33=reg29*reg17;
    T reg72=pow(reg42,2); reg26=reg26/reg32; reg40=reg12*reg40; T reg73=reg24+reg3; reg38=reg8*reg38;
    T reg74=reg46*reg49; T reg75=pow(reg49,2); T reg76=reg57*reg18; T reg77=reg74*reg15; T reg78=reg36+reg41;
    T reg79=reg58*reg43; T reg80=reg71*reg18; T reg81=reg51*reg19; reg9=reg13*reg9; reg65=reg66+reg65;
    reg66=reg37*reg17; reg67=reg50-reg67; reg76=reg33+reg76; reg33=reg75*reg26; reg50=reg55*reg19;
    T reg82=reg56*reg43; reg39=reg12*reg39; reg70=reg52+reg70; reg52=reg57*reg43; T reg83=reg29*reg19;
    reg40=reg38-reg40; reg73=reg4*reg73; reg48=reg8*reg48; reg5=reg5/reg32; reg38=pow(reg59,2);
    T reg84=reg74*reg26; reg69=reg68+reg69; reg60=reg31+reg60; reg31=reg72*reg15; reg62=reg61+reg62;
    reg64=reg63+reg64; reg53=reg12*reg53; reg61=reg72*reg26; reg63=reg71*reg28; reg68=reg37*reg30;
    T reg85=reg75*reg15; reg70=reg12*reg70; T reg86=reg27*reg11; T reg87=reg22*reg21; reg53=reg48-reg53;
    reg48=reg39+reg9; T reg88=2*reg23; reg65=reg61+reg65; reg61=2*reg25; reg67=reg8*reg67;
    T reg89=reg74*reg5; reg79=reg81+reg79; reg73=reg40-reg73; reg78=reg4*reg78; reg40=reg38*reg26;
    reg77=reg64+reg77; reg80=reg66+reg80; reg64=reg38*reg15; reg63=reg68+reg63; reg85=reg62+reg85;
    reg33=reg76+reg33; reg82=reg50+reg82; reg50=reg72*reg5; reg62=reg14*reg25; reg52=reg83+reg52;
    reg66=reg16*reg23; reg68=reg75*reg5; reg31=reg60+reg31; reg60=reg37*reg19; reg76=reg71*reg43;
    reg84=reg69+reg84; reg69=reg86*reg17; reg81=reg55*reg33; reg83=reg85*reg56; T reg90=reg88*reg21;
    T reg91=reg87*reg18; T reg92=reg61*reg23; T reg93=reg29*reg84; T reg94=reg16*reg11; T reg95=reg77*reg57;
    reg46=reg46*reg59; T reg96=2*reg49; T reg97=reg14*reg16; T reg98=reg25*reg23; reg47=reg54+reg47;
    reg54=2*reg27; T reg99=reg25*reg42; T reg100=reg14*reg49; T reg101=reg85*reg66; T reg102=reg62*reg33;
    T reg103=reg49*reg42; reg73=reg73/reg32; T reg104=reg38*reg5; reg76=reg60+reg76; reg89=reg79+reg89;
    reg68=reg52+reg68; reg48=reg4*reg48; reg50=reg82+reg50; reg70=reg67-reg70; reg52=reg62*reg84;
    reg60=reg77*reg66; reg67=reg29*reg65; reg79=reg31*reg57; reg82=reg87*reg28; reg78=reg53-reg78;
    reg53=reg86*reg30; reg40=reg80+reg40; reg33=reg29*reg33; reg85=reg85*reg57; reg64=reg63+reg64;
    reg63=reg55*reg65; reg80=reg31*reg56; T reg105=reg61*reg27; reg77=reg77*reg56; reg84=reg55*reg84;
    reg17=reg105*reg17; T reg106=reg46*reg26; reg91=reg69+reg91; reg18=reg90*reg18; reg48=reg70-reg48;
    reg78=reg78/reg32; reg69=reg16*reg49; reg70=reg23*reg42; T reg107=reg89*reg103; reg60=reg52+reg60;
    reg65=reg62*reg65; reg52=reg14*reg42; T reg108=reg27*reg21; T reg109=reg56*reg8; T reg110=reg57*reg12;
    reg96=reg96*reg59; T reg111=reg89*reg75; T reg112=reg94*reg8; T reg113=reg92*reg12; reg95=reg93+reg95;
    reg93=reg64*reg57; T reg114=reg56*reg12; T reg115=reg29*reg40; T reg116=reg57*reg13; T reg117=reg68*reg75;
    T reg118=reg94*reg12; reg85=reg33+reg85; reg79=reg67+reg79; reg33=reg50*reg75; reg67=reg92*reg13;
    reg101=reg102+reg101; reg102=reg68*reg103; reg68=reg68*reg72; reg83=reg81+reg83; reg81=reg50*reg72;
    reg80=reg63+reg80; reg63=reg47*reg73; T reg119=reg98*reg73; T reg120=reg97*reg73; T reg121=reg87*reg43;
    T reg122=reg86*reg19; reg104=reg76+reg104; reg28=reg90*reg28; reg30=reg105*reg30; reg76=reg46*reg15;
    reg82=reg53+reg82; reg54=reg21*reg54; reg53=reg55*reg40; T reg123=reg64*reg56; reg31=reg31*reg66;
    reg100=reg99+reg100; reg77=reg84+reg77; reg89=reg89*reg72; reg84=reg25*reg49; reg31=reg65+reg31;
    reg50=reg50*reg103; reg65=reg23*reg49; reg99=reg16*reg42; T reg124=reg27*reg59; reg61=reg61*reg49;
    T reg125=reg56*reg4; T reg126=reg57*reg7; T reg127=reg4*reg94; T reg128=reg92*reg7; reg11=reg42*reg11;
    reg69=reg70+reg69; reg64=reg64*reg66; reg40=reg62*reg40; reg70=reg71*reg4; reg110=reg109-reg110;
    reg109=reg4*reg54; reg113=reg112-reg113; reg112=reg71*reg7; reg114=reg116-reg114; reg116=reg54*reg7;
    reg118=reg67-reg118; reg67=reg119*reg47; reg102=reg101+reg102; reg101=reg119*reg92; T reg129=reg120*reg94;
    reg81=reg80+reg81; reg93=reg115+reg93; reg80=reg104*reg75; reg115=reg63*reg94; reg89=reg77+reg89;
    reg111=reg95+reg111; reg77=reg100*reg78; reg95=reg63*reg92; T reg130=reg84*reg78; T reg131=reg52*reg78;
    T reg132=reg108*reg73; reg43=reg90*reg43; reg19=reg105*reg19; T reg133=reg46*reg5; reg121=reg122+reg121;
    reg15=reg96*reg15; reg28=reg30+reg28; reg76=reg82+reg76; reg26=reg96*reg26; reg18=reg17+reg18;
    reg17=reg104*reg72; reg123=reg53+reg123; reg117=reg85+reg117; reg68=reg83+reg68; reg119=reg119*reg94;
    reg32=reg48/reg32; reg30=reg120*reg92; reg106=reg91+reg106; reg33=reg79+reg33; reg63=reg63*reg47;
    reg107=reg60+reg107; reg13=reg29*reg13; reg48=reg27*reg42; reg116=reg118-reg116; reg53=reg14*reg59;
    reg112=reg114-reg112; reg60=(*f.m).alpha_2*reg57; reg128=reg127+reg128; reg54=reg54*reg6; reg126=reg125+reg126;
    reg79=reg71*reg6; reg82=reg132*reg94; reg17=reg123+reg17; reg22=reg22*reg42; reg88=reg88*reg49;
    reg83=reg21*reg59; reg85=(*f.m).alpha_1*reg29; reg101=reg117+reg101; reg91=reg130*reg61; reg114=reg131*reg61;
    reg30=reg33+reg30; reg80=reg93+reg80; reg33=reg132*reg92; reg115=reg89+reg115; reg89=reg77*reg11;
    reg93=reg55*reg106; reg67=reg102+reg67; reg102=reg130*reg100; reg104=reg104*reg103; reg64=reg40+reg64;
    reg40=reg76*reg56; reg95=reg111+reg95; reg111=reg77*reg61; reg117=reg29*reg12; reg44=reg45+reg44;
    reg45=reg29*reg106; reg118=reg76*reg57; reg70=reg110-reg70; reg109=reg113-reg109; reg110=(*f.m).alpha_1*reg55;
    reg12=reg55*reg12; reg5=reg96*reg5; reg43=reg19+reg43; reg133=reg121+reg133; reg50=reg31+reg50;
    reg15=reg28+reg15; reg26=reg18+reg26; reg18=reg124*reg78; reg120=reg120*reg47; reg19=(*f.m).alpha_2*reg56;
    reg28=reg99*reg32; reg31=reg65*reg32; reg130=reg130*reg11; reg119=reg68+reg119; reg129=reg81+reg129;
    reg77=reg77*reg100; reg8=reg55*reg8; reg68=reg131*reg11; reg63=reg107+reg63; reg81=reg69*reg32;
    reg107=reg109*reg55; reg113=reg16*reg59; reg121=(*f.m).alpha_1*reg37; reg104=reg64+reg104; reg64=reg44*reg73;
    reg122=reg112*reg57; reg123=reg70*reg56; reg111=reg95+reg111; reg95=reg81*reg88; reg5=reg43+reg5;
    reg43=reg112*reg29; reg125=reg37*reg4; reg127=reg70*reg55; reg40=reg93+reg40; reg68=reg129+reg68;
    reg93=reg109*reg56; reg129=reg116*reg57; reg114=reg30+reg114; reg30=reg83*reg32; reg33=reg80+reg33;
    reg80=reg18*reg61; reg131=reg131*reg100; reg89=reg115+reg89; reg115=reg81*reg22; T reg134=reg28*reg88;
    T reg135=reg31*reg88; reg91=reg101+reg91; reg101=reg21*reg42; reg120=reg50+reg120; reg50=(*f.m).alpha_2*reg71;
    T reg136=reg28*reg22; reg132=reg132*reg47; T reg137=reg116*reg29; T reg138=reg26*reg55; reg82=reg17+reg82;
    reg126=reg79-reg126; reg17=reg15*reg56; reg79=reg29*reg7; reg53=reg48+reg53; reg4=reg55*reg4;
    reg76=reg76*reg66; reg106=reg62*reg106; reg48=reg109*reg97; reg19=reg110+reg19; reg110=(*f.m).alpha_3*reg72;
    reg7=reg37*reg7; T reg139=reg31*reg22; reg130=reg119+reg130; reg119=reg27*reg49; T reg140=reg25*reg59;
    reg12=reg13-reg12; reg81=reg81*reg69; reg77=reg63+reg77; reg13=reg116*reg98; reg63=reg26*reg29;
    T reg141=reg15*reg57; T reg142=reg133*reg72; reg60=reg85+reg60; reg117=reg8-reg117; reg8=reg133*reg75;
    reg85=reg70*reg97; T reg143=reg112*reg98; reg102=reg67+reg102; reg128=reg54-reg128; reg118=reg45+reg118;
    reg45=reg18*reg11; reg54=(*f.m).alpha_3*reg75; reg31=reg31*reg69; reg35=reg34+reg35; reg134=reg114+reg134;
    reg31=reg102+reg31; reg133=reg133*reg103; reg7=reg12-reg7; reg12=reg5*reg72; reg34=reg5*reg75;
    reg81=reg77+reg81; reg125=reg117-reg125; reg142=reg40+reg142; reg136=reg68+reg136; reg26=reg26*reg62;
    reg40=reg64*reg94; reg76=reg106+reg76; reg17=reg138+reg17; reg132=reg104+reg132; reg115=reg89+reg115;
    reg18=reg18*reg100; reg15=reg15*reg66; reg67=reg23*reg59; reg95=reg111+reg95; reg8=reg118+reg8;
    reg68=reg30*reg22; reg77=reg64*reg92; reg89=reg47*reg2; reg102=reg66*reg2; reg143=reg85+reg143;
    reg45=reg82+reg45; reg82=reg126*reg108; reg131=reg120+reg131; reg79=reg4+reg79; reg6=reg37*reg6;
    reg28=reg28*reg69; reg139=reg130+reg139; reg110=reg19+reg110; reg140=reg119+reg140; reg14=reg14*reg27;
    reg16=reg16*reg21; reg13=reg48+reg13; reg4=reg128*reg108; reg122=reg123+reg122; reg19=reg126*reg71;
    reg135=reg91+reg135; reg129=reg93+reg129; reg48=reg128*reg71; reg80=reg33+reg80; reg66=(*f.m).alpha_2*reg66;
    reg33=(*f.m).alpha_1*reg62; reg85=reg30*reg88; reg141=reg63+reg141; reg63=(*f.m).alpha_3*reg38; reg50=reg121+reg50;
    reg91=reg53*reg78; reg93=reg128*reg37; reg137=reg107+reg137; reg54=reg60+reg54; reg73=reg35*reg73;
    reg60=reg126*reg37; reg43=reg127+reg43; reg113=reg101+reg113; reg101=reg21*reg49; reg85=reg80+reg85;
    reg28=reg131+reg28; reg77=reg8+reg77; reg8=reg91*reg61; reg80=reg89*reg51; reg93=reg137+reg93;
    reg104=reg81*reg139; reg106=reg31*reg95; reg107=reg102*reg51; reg60=reg43+reg60; reg67=reg101+reg67;
    reg68=reg45+reg68; reg43=reg31*reg115; reg45=reg7*reg29; reg101=reg125*reg55; reg111=reg44*reg1;
    reg114=reg16*reg1; reg79=reg6-reg79; reg6=reg54*reg139; reg42=reg59*reg42; reg27=reg25*reg27;
    reg21=reg23*reg21; reg23=reg89*reg47; reg4=reg13+reg4; reg13=reg136*reg110; reg25=reg102*reg47;
    reg82=reg143+reg82; reg78=reg140*reg78; reg117=reg113*reg32; reg5=reg5*reg103; reg15=reg26+reg15;
    reg64=reg64*reg47; reg133=reg76+reg133; reg26=reg7*reg57; reg92=reg73*reg92; reg34=reg141+reg34;
    reg76=reg125*reg56; reg30=reg30*reg69; reg18=reg132+reg18; reg40=reg142+reg40; reg118=reg91*reg11;
    reg19=reg122+reg19; reg12=reg17+reg12; reg17=reg102*reg58; reg94=reg73*reg94; reg119=reg110*reg134;
    reg120=reg54*reg135; reg16=(*f.m).alpha_2*reg16; reg121=(*f.m).alpha_1*reg14; reg48=reg129+reg48; reg122=reg89*reg58;
    reg103=(*f.m).alpha_3*reg103; reg66=reg33+reg66; reg63=reg50+reg63; reg33=reg81*reg135; reg120=reg119+reg120;
    reg50=reg110*reg28; reg119=reg54*reg31; reg85=reg63*reg85; reg61=reg78*reg61; reg92=reg34+reg92;
    reg6=reg13+reg6; reg68=reg63*reg68; reg13=(*f.m).alpha_2*reg21; reg34=reg115*reg135; reg123=(*f.m).alpha_1*reg27;
    reg127=reg139*reg95; reg42=(*f.m).alpha_3*reg42; reg43=reg104-reg43; reg16=reg121+reg16; reg106=reg33-reg106;
    reg103=reg66+reg103; reg33=reg111*reg44; reg23=reg4+reg23; reg4=reg111*reg86; reg80=reg93+reg80;
    reg66=reg114*reg86; reg107=reg60+reg107; reg60=reg114*reg44; reg25=reg82+reg25; reg37=reg79*reg37;
    reg45=reg101+reg45; reg82=reg35*reg0; reg21=reg21*reg0; reg93=reg7*reg98; reg2=reg62*reg2;
    reg62=reg125*reg97; reg49=reg59*reg49; reg32=reg67*reg32; reg73=reg73*reg47; reg64=reg133+reg64;
    reg91=reg91*reg100; reg5=reg15+reg5; reg30=reg18+reg30; reg15=reg117*reg88; reg26=reg76+reg26;
    reg71=reg79*reg71; reg8=reg77+reg8; reg11=reg78*reg11; reg94=reg12+reg94; reg17=reg19+reg17;
    reg12=reg114*reg87; reg18=reg111*reg87; reg19=reg117*reg22; reg118=reg40+reg118; reg122=reg48+reg122;
    reg15=reg8+reg15; reg119=reg50+reg119; reg30=reg63*reg30; reg73=reg5+reg73; reg78=reg78*reg100;
    reg68=reg6+reg68; reg5=reg82*reg35; reg33=reg23+reg33; reg6=reg103*reg95; reg85=reg120+reg85;
    reg22=reg32*reg22; reg11=reg94+reg11; reg8=reg103*reg115; reg93=reg62+reg93; reg108=reg79*reg108;
    reg19=reg118+reg19; reg23=reg21*reg35; reg60=reg25+reg60; reg42=reg16+reg42; reg12=reg17+reg12;
    reg16=reg21*reg90; reg13=reg123+reg13; reg18=reg122+reg18; reg49=(*f.m).alpha_3*reg49; reg17=reg82*reg90;
    reg25=reg82*reg105; reg4=reg80+reg4; reg40=reg21*reg105; reg66=reg107+reg66; reg48=reg136*reg106;
    reg50=reg2*reg51; reg59=reg2*reg58; reg71=reg26+reg71; reg26=reg43*reg134; reg37=reg45+reg37;
    reg1=reg14*reg1; reg34=reg127-reg34; reg61=reg92+reg61; reg117=reg117*reg69; reg91=reg64+reg91;
    reg88=reg32*reg88; reg59=reg71+reg59; reg14=1-var_inter[1]; reg88=reg61+reg88; reg22=reg11+reg22;
    reg30=reg119+reg30; reg11=reg103*reg81; reg49=reg13+reg49; reg87=reg1*reg87; reg16=reg12+reg16;
    reg5=reg33+reg5; reg17=reg18+reg17; reg23=reg60+reg23; reg12=reg2*reg47; reg108=reg93+reg108;
    reg13=1-var_inter[0]; reg32=reg32*reg69; reg78=reg73+reg78; reg117=reg91+reg117; reg6=reg85+reg6;
    reg19=reg42*reg19; reg8=reg68+reg8; reg18=reg28*reg115; reg15=reg42*reg15; reg33=reg28*reg95;
    reg45=reg136*reg81; reg60=reg81*reg134; reg61=reg34*reg28; reg0=reg27*reg0; reg50=reg37+reg50;
    reg86=reg1*reg86; reg26=reg48-reg26; reg25=reg4+reg25; reg40=reg66+reg40; reg4=reg13*elem.pos(0)[0];
    reg27=reg17*reg23; reg37=reg16*reg5; reg32=reg78+reg32; reg48=reg14*elem.pos(1)[1]; reg62=reg14*elem.pos(0)[1];
    reg64=reg25*reg23; reg66=reg31*reg134; reg68=reg14*elem.pos(0)[0]; reg71=reg14*elem.pos(1)[0]; reg90=reg0*reg90;
    reg105=reg0*reg105; reg86=reg50+reg86; reg87=reg59+reg87; reg33=reg60-reg33; reg50=reg28*reg135;
    reg59=reg136*reg31; reg18=reg45-reg18; reg45=reg28*reg139; reg60=reg136*reg95; reg73=reg115*reg134;
    reg19=reg8+reg19; reg22=reg49*reg22; reg8=elem.pos(1)[1]*var_inter[0]; reg76=elem.pos(0)[1]*reg13; reg77=elem.pos(1)[0]*var_inter[0];
    reg78=reg40*reg5; reg44=reg1*reg44; reg12=reg108+reg12; reg88=reg49*reg88; reg15=reg6+reg15;
    reg61=reg26+reg61; reg11=reg30+reg11; reg117=reg42*reg117; reg6=reg136*reg135; reg73=reg60-reg73;
    reg26=reg139*reg134; reg30=reg25*reg16; reg60=reg40*reg17; reg64=reg78-reg64; reg22=reg19+reg22;
    reg19=reg76+reg8; reg78=elem.pos(2)[0]*var_inter[0]; reg80=reg4+reg77; reg32=reg49*reg32; reg85=elem.pos(2)[1]*var_inter[0];
    reg45=reg59-reg45; reg117=reg11+reg117; reg18=reg18/reg61; reg43=reg43/reg61; reg50=reg66-reg50;
    reg88=reg15+reg88; reg105=reg86+reg105; reg90=reg87+reg90; reg35=reg0*reg35; reg33=reg33/reg61;
    reg44=reg12+reg44; reg27=reg37-reg27; reg71=reg71-reg68; reg11=var_inter[1]*elem.pos(2)[0]; reg12=var_inter[1]*elem.pos(2)[1];
    reg48=reg48-reg62; reg106=reg106/reg61; reg32=reg117+reg32; reg73=reg73/reg61; reg15=reg64*reg90;
    reg50=reg50/reg61; reg34=reg34/reg61; reg37=var_inter[1]*elem.pos(3)[0]; reg45=reg45/reg61; reg11=reg71+reg11;
    reg106=reg106*reg22; reg85=reg85-reg19; reg43=reg43*reg88; reg59=reg105*reg27; reg66=var_inter[1]*elem.pos(3)[1];
    reg33=reg33*reg22; reg12=reg48+reg12; reg78=reg78-reg80; reg48=reg13*elem.pos(3)[0]; reg18=reg18*reg88;
    reg30=reg60-reg30; reg35=reg44+reg35; reg26=reg6-reg26; reg6=reg13*elem.pos(3)[1]; reg22=reg50*reg22;
    reg33=reg18-reg33; reg43=reg106-reg43; reg6=reg85+reg6; reg18=reg30*reg35; reg73=reg73*reg32;
    reg34=reg34*reg32; reg44=reg90*reg23; reg50=reg40*reg35; reg11=reg11-reg37; reg88=reg45*reg88;
    reg61=reg26/reg61; reg15=reg59-reg15; reg12=reg12-reg66; reg26=reg105*reg23; reg48=reg78+reg48;
    reg45=reg16*reg35; reg59=reg25*reg35; reg45=reg44-reg45; reg44=reg17*reg35; reg60=reg90*reg5;
    reg71=reg105*reg5; reg50=reg26-reg50; reg18=reg15+reg18; reg15=reg105*reg16; reg26=reg40*reg90;
    reg78=reg11*reg6; reg34=reg43+reg34; reg73=reg33-reg73; reg33=1-(*f.m).resolution; reg43=reg12*reg48;
    reg32=reg61*reg32; reg88=reg22-reg88; reg59=reg71-reg59; reg50=reg50/reg18; reg88=reg32+reg88;
    reg45=reg45/reg18; reg22=reg105*reg17; reg34=reg34*reg33; reg32=reg25*reg90; reg73=reg73*reg33;
    reg26=reg15-reg26; reg15=(*f.m).resolution*reg54; reg43=reg78-reg43; reg61=(*f.m).resolution*reg110; reg44=reg60-reg44;
    reg60=(*f.m).resolution*reg103; reg71=(*f.m).resolution*reg45; reg78=(*f.m).resolution*reg50; reg27=reg27/reg18; reg31=reg33*reg31;
    reg44=reg44/reg18; reg64=reg64/reg18; reg28=reg33*reg28; reg34=reg61+reg34; reg11=reg11/reg43;
    reg12=reg12/reg43; reg48=reg48/reg43; reg6=reg6/reg43; reg88=reg88*reg33; reg26=reg26/reg18;
    reg32=reg22-reg32; reg59=reg59/reg18; reg73=reg15+reg73; reg73=reg73*(*f.m).deltaT; reg15=var_inter[0]*reg12;
    reg22=reg13*reg11; reg34=reg34*(*f.m).deltaT; reg61=var_inter[1]*reg48; reg88=reg60+reg88; reg60=var_inter[0]*reg11;
    reg85=reg13*reg12; reg86=var_inter[1]*reg6; reg87=reg14*reg6; reg91=reg14*reg48; reg32=reg32/reg18;
    reg92=(*f.m).resolution*reg64; reg71=reg28+reg71; reg28=(*f.m).resolution*reg26; reg78=reg31-reg78; reg31=(*f.m).resolution*reg59;
    reg93=(*f.m).resolution*reg44; reg136=reg136*reg33; reg81=reg33*reg81; reg135=reg33*reg135; reg18=reg30/reg18;
    reg30=(*f.m).resolution*reg27; reg134=reg33*reg134; reg139=reg33*reg139; reg28=reg81+reg28; reg81=(*f.m).resolution*reg18;
    reg94=reg91-reg22; reg101=reg85-reg87; reg104=reg22+reg61; reg106=reg85+reg86; reg107=reg73*reg78;
    reg108=reg34*reg71; reg117=reg91+reg60; reg115=reg33*reg115; reg118=reg86-reg15; reg88=reg88*(*f.m).deltaT;
    reg31=reg135+reg31; reg95=reg33*reg95; reg93=reg134-reg93; reg33=(*f.m).resolution*reg32; reg119=reg87+reg15;
    reg136=reg30+reg136; reg92=reg139-reg92; reg30=reg60-reg61; reg120=0.5*reg119; reg121=0.5*reg94;
    reg122=0.5*reg30; reg123=0.5*reg101; reg127=0.5*reg104; reg129=0.5*reg106; reg130=0.5*reg117;
    reg131=0.5*reg118; reg132=reg108+reg107; reg133=reg88*reg28; reg134=reg34*reg136; reg135=reg34*reg93;
    reg137=reg73*reg31; reg138=reg73*reg92; reg115=reg81+reg115; reg33=reg95-reg33; reg81=reg30*reg78;
    reg95=reg131*reg28; reg139=reg123*reg28; reg141=reg88*reg33; reg142=reg135+reg137; reg143=reg122*reg28;
    T reg144=reg118*reg71; T reg145=reg132+reg133; T reg146=reg94*reg78; T reg147=reg120*reg28; T reg148=reg117*reg78;
    T reg149=reg88*reg115; T reg150=reg138+reg134; T reg151=reg130*reg28; T reg152=reg119*reg71; T reg153=reg129*reg28;
    T reg154=reg121*reg28; T reg155=reg101*reg71; T reg156=reg104*reg78; T reg157=reg127*reg28; T reg158=reg106*reg71;
    reg139=reg146+reg139; reg146=reg150+reg149; T reg159=2*reg145; T reg160=reg142+reg141; T reg161=reg101*reg136;
    T reg162=reg121*reg115; reg154=reg155+reg154; reg155=reg119*reg136; T reg163=reg130*reg115; reg152=reg152-reg151;
    T reg164=reg117*reg92; T reg165=reg120*reg115; reg147=reg147-reg148; T reg166=reg94*reg92; T reg167=reg123*reg115;
    T reg168=reg118*reg136; T reg169=reg122*reg115; reg143=reg144+reg143; reg144=reg30*reg92; T reg170=reg131*reg115;
    reg95=reg81+reg95; reg81=reg119*reg93; T reg171=reg30*reg31; T reg172=reg131*reg33; T reg173=reg122*reg33;
    T reg174=reg118*reg93; reg156=reg156-reg153; T reg175=reg123*reg33; T reg176=reg129*reg115; T reg177=reg104*reg92;
    T reg178=reg120*reg33; T reg179=reg94*reg31; T reg180=reg106*reg93; T reg181=reg127*reg33; T reg182=reg127*reg115;
    T reg183=reg106*reg136; T reg184=reg129*reg33; T reg185=reg104*reg31; reg157=reg157-reg158; T reg186=reg130*reg33;
    T reg187=reg101*reg93; T reg188=reg117*reg31; T reg189=reg121*reg33; reg147=2*reg147; reg173=reg174+reg173;
    reg189=reg187+reg189; reg81=reg81-reg186; reg162=reg161+reg162; reg154=2*reg154; reg165=reg165-reg164;
    reg161=reg159*reg127; reg174=reg146*reg106; reg152=2*reg152; reg175=reg179+reg175; reg179=reg159*reg129;
    reg155=reg155-reg163; reg187=reg160*reg104; reg178=reg178-reg188; reg182=reg182-reg183; reg139=2*reg139;
    reg157=2*reg157; T reg190=reg146*reg119; reg181=reg181-reg180; T reg191=reg159*reg130; reg167=reg166+reg167;
    reg95=2*reg95; reg170=reg144+reg170; reg144=reg159*reg120; reg172=reg171+reg172; reg177=reg177-reg176;
    reg185=reg185-reg184; reg166=reg160*reg117; reg169=reg168+reg169; reg143=2*reg143; reg156=2*reg156;
    reg168=reg152*reg129; reg171=reg147*reg131; T reg192=reg30*reg175; T reg193=reg131*reg139; T reg194=reg95*reg129;
    T reg195=reg147*reg129; T reg196=reg104*reg178; T reg197=reg122*reg157; T reg198=reg118*reg182; T reg199=reg104*reg172;
    T reg200=reg147*reg130; T reg201=reg169*reg119; T reg202=reg30*reg81; T reg203=reg104*reg81; T reg204=reg118*reg177;
    T reg205=reg30*reg178; T reg206=reg104*reg173; T reg207=reg143*reg129; T reg208=reg122*reg156; T reg209=reg129*reg156;
    T reg210=reg104*reg181; T reg211=reg152*reg131; T reg212=reg185*reg104; T reg213=reg165*reg119; T reg214=reg30*reg189;
    T reg215=reg131*reg154; T reg216=reg129*reg157; T reg217=reg118*reg162; T reg218=reg160*reg94; T reg219=reg159*reg123;
    T reg220=reg185*reg117; T reg221=reg120*reg156; T reg222=reg147*reg120; T reg223=reg95*reg122; T reg224=reg117*reg178;
    T reg225=reg117*reg181; T reg226=reg120*reg157; T reg227=reg185*reg94; T reg228=reg123*reg156; T reg229=reg117*reg172;
    T reg230=reg95*reg120; T reg231=reg130*reg139; T reg232=reg119*reg162; T reg233=reg130*reg154; T reg234=reg117*reg173;
    T reg235=reg143*reg120; T reg236=reg119*reg167; T reg237=reg170*reg118; T reg238=reg146*reg101; T reg239=reg94*reg189;
    T reg240=reg143*reg122; T reg241=reg169*reg118; T reg242=reg123*reg154; T reg243=reg117*reg81; T reg244=reg147*reg122;
    T reg245=reg165*reg118; T reg246=reg159*reg121; T reg247=reg152*reg130; T reg248=reg152*reg122; T reg249=reg155*reg118;
    T reg250=reg94*reg175; T reg251=reg123*reg139; T reg252=reg122*reg139; T reg253=reg118*reg167; T reg254=reg155*reg119;
    T reg255=reg152*reg120; T reg256=reg122*reg154; T reg257=reg30*reg181; T reg258=reg147*reg121; T reg259=reg165*reg101;
    T reg260=reg131*reg157; T reg261=reg130*reg156; reg185=reg185*reg30; T reg262=reg131*reg156; T reg263=reg119*reg182;
    T reg264=reg152*reg121; T reg265=reg127*reg154; T reg266=reg147*reg123; reg178=reg94*reg178; T reg267=reg155*reg101;
    T reg268=reg106*reg162; T reg269=reg127*reg139; T reg270=reg121*reg139; T reg271=reg106*reg167; T reg272=reg95*reg130;
    T reg273=reg121*reg154; T reg274=reg101*reg182; T reg275=reg121*reg157; T reg276=reg123*reg157; reg181=reg94*reg181;
    reg167=reg101*reg167; T reg277=reg95*reg121; T reg278=reg95*reg123; T reg279=reg170*reg101; T reg280=reg94*reg172;
    T reg281=reg101*reg177; T reg282=reg143*reg121; T reg283=reg169*reg101; T reg284=reg121*reg156; T reg285=reg130*reg157;
    reg172=reg30*reg172; T reg286=reg95*reg131; T reg287=reg119*reg177; T reg288=reg143*reg123; T reg289=reg94*reg173;
    T reg290=reg159*reg131; T reg291=reg160*reg30; reg157=reg127*reg157; reg182=reg106*reg182; T reg292=reg143*reg130;
    T reg293=reg159*reg122; T reg294=reg146*reg118; reg156=reg127*reg156; reg177=reg106*reg177; T reg295=reg120*reg139;
    T reg296=reg166-reg144; T reg297=reg104*reg189; T reg298=reg129*reg154; T reg299=reg117*reg175; T reg300=reg191-reg190;
    reg175=reg104*reg175; reg139=reg129*reg139; T reg301=reg143*reg131; reg173=reg30*reg173; T reg302=reg152*reg127;
    reg155=reg155*reg106; reg162=reg101*reg162; reg154=reg120*reg154; reg152=reg152*reg123; reg81=reg94*reg81;
    reg147=reg147*reg127; reg165=reg165*reg106; reg189=reg117*reg189; T reg303=reg179-reg187; reg143=reg143*reg127;
    reg169=reg169*reg106; T reg304=reg170*reg119; T reg305=reg174-reg161; reg95=reg95*reg127; reg170=reg170*reg106;
    reg247=reg254-reg247; reg200=reg213-reg200; reg272=reg304-reg272; reg285=reg263-reg285; reg231=reg236-reg231;
    reg292=reg201-reg292; reg275=reg274+reg275; reg284=reg281+reg284; reg286=reg172+reg286; reg260=reg257+reg260;
    reg262=reg185+reg262; reg268=reg265-reg268; reg271=reg269-reg271; reg155=reg302-reg155; reg165=reg147-reg165;
    reg169=reg143-reg169; reg170=reg95-reg170; reg182=reg157-reg182; reg177=reg156-reg177; reg298=reg297-reg298;
    reg139=reg175-reg139; reg195=reg196-reg195; reg207=reg206-reg207; reg194=reg199-reg194; reg216=reg210-reg216;
    reg209=reg212-reg209; reg242=reg239+reg242; reg95=reg238+reg246; reg251=reg250+reg251; reg143=reg218+reg219;
    reg223=reg237+reg223; reg228=reg227+reg228; reg233=reg232-reg233; reg229=reg230-reg229; reg225=reg226-reg225;
    reg220=reg221-reg220; reg256=reg217+reg256; reg252=reg253+reg252; reg248=reg249+reg248; reg244=reg245+reg244;
    reg240=reg241+reg240; reg197=reg198+reg197; reg208=reg204+reg208; reg215=reg214+reg215; reg193=reg192+reg193;
    reg211=reg202+reg211; reg171=reg205+reg171; reg301=reg173+reg301; reg300=reg43*reg300; reg296=reg43*reg296;
    reg147=reg294+reg293; reg156=reg291+reg290; reg305=reg43*reg305; reg303=reg43*reg303; reg152=reg81+reg152;
    reg273=reg162+reg273; reg167=reg270+reg167; reg266=reg178+reg266; reg264=reg267+reg264; reg258=reg259+reg258;
    reg288=reg289+reg288; reg282=reg283+reg282; reg278=reg280+reg278; reg276=reg181+reg276; reg277=reg279+reg277;
    reg299=reg295-reg299; reg243=reg255-reg243; reg234=reg235-reg234; reg224=reg222-reg224; reg261=reg287-reg261;
    reg168=reg203-reg168; reg189=reg154-reg189; reg165=reg165*reg43; reg305=ponderation*reg305; reg169=reg169*reg43;
    reg81=reg43*reg156; reg170=reg170*reg43; reg154=reg43*reg147; reg200=reg200*reg43; reg182=reg182*reg43;
    reg296=ponderation*reg296; reg299=reg299*reg43; reg177=reg177*reg43; reg300=ponderation*reg300; reg298=reg298*reg43;
    reg139=reg139*reg43; reg171=reg171*reg43; reg195=reg195*reg43; reg211=reg211*reg43; reg247=reg247*reg43;
    reg207=reg207*reg43; reg193=reg193*reg43; reg168=reg168*reg43; reg285=reg285*reg43; reg277=reg43*reg277;
    reg276=reg43*reg276; reg275=reg43*reg275; reg278=reg43*reg278; reg282=reg282*reg43; reg288=reg43*reg288;
    reg272=reg272*reg43; reg301=reg301*reg43; reg258=reg258*reg43; reg261=reg261*reg43; reg286=reg286*reg43;
    reg266=reg43*reg266; reg260=reg260*reg43; reg264=reg264*reg43; reg262=reg262*reg43; reg167=reg167*reg43;
    reg268=reg268*reg43; reg273=reg43*reg273; reg292=reg292*reg43; reg271=reg271*reg43; reg152=reg43*reg152;
    reg189=reg189*reg43; reg155=reg155*reg43; reg303=ponderation*reg303; reg209=reg209*reg43; reg225=reg225*reg43;
    reg284=reg43*reg284; reg197=reg197*reg43; reg157=reg143*reg43; reg233=reg233*reg43; reg240=reg240*reg43;
    reg231=reg231*reg43; reg220=reg220*reg43; reg244=reg244*reg43; reg242=reg43*reg242; reg256=reg256*reg43;
    reg251=reg43*reg251; reg248=reg248*reg43; reg162=reg95*reg43; reg224=reg224*reg43; reg252=reg252*reg43;
    reg234=reg234*reg43; reg194=reg194*reg43; reg215=reg215*reg43; reg228=reg228*reg43; reg216=reg216*reg43;
    reg208=reg208*reg43; reg229=reg229*reg43; reg223=reg223*reg43; reg243=reg243*reg43; T tmp_2_5=ponderation*reg272;
    T tmp_2_4=ponderation*reg292; T tmp_6_0=ponderation*reg268; T tmp_1_7=ponderation*reg228; T tmp_0_5=ponderation*reg277; T tmp_4_0=ponderation*reg256;
    T tmp_1_4=ponderation*reg288; T tmp_0_1=ponderation*reg167; T tmp_1_1=ponderation*reg251; T tmp_5_7=ponderation*reg262; T tmp_1_6=ponderation*reg276;
    T tmp_2_7=ponderation*reg261; T tmp_2_0=ponderation*reg233; T tmp_0_4=ponderation*reg282; T tmp_3_3=ponderation*reg224; T tmp_0_2=ponderation*reg264;
    T tmp_5_6=ponderation*reg260; T tmp_3_7=ponderation*reg220; T tmp_2_6=ponderation*reg285; T tmp_1_5=ponderation*reg278; T tmp_1_3=ponderation*reg266;
    T tmp_0_0=ponderation*reg273; T tmp_3_4=ponderation*reg234; T tmp_5_5=ponderation*reg286; T tmp_0_6=ponderation*reg275; T tmp_3_6=ponderation*reg225;
    reg167=ponderation*reg157; T vec_1=reg167; T tmp_0_3=ponderation*reg258; T tmp_3_5=ponderation*reg229; T tmp_5_4=ponderation*reg301;
    T vec_3=-reg296; T tmp_0_7=ponderation*reg284; T tmp_6_7=ponderation*reg177; T vec_2=-reg300; T tmp_7_7=ponderation*reg209;
    T tmp_3_1=ponderation*reg299; T tmp_7_0=ponderation*reg298; T tmp_4_5=ponderation*reg223; T tmp_4_6=ponderation*reg197; T tmp_5_3=ponderation*reg171;
    T tmp_7_1=ponderation*reg139; T tmp_7_2=ponderation*reg168; T tmp_2_2=ponderation*reg247; T tmp_5_2=ponderation*reg211; T tmp_7_6=ponderation*reg216;
    T tmp_7_3=ponderation*reg195; T tmp_4_7=ponderation*reg208; T tmp_5_1=ponderation*reg193; T tmp_7_4=ponderation*reg207; T tmp_7_5=ponderation*reg194;
    T tmp_5_0=ponderation*reg215; T tmp_4_1=ponderation*reg252; T tmp_6_1=ponderation*reg271; T tmp_1_2=ponderation*reg152; reg139=ponderation*reg162;
    T vec_0=reg139; T tmp_6_2=ponderation*reg155; T vec_7=-reg303; T tmp_4_2=ponderation*reg248; T tmp_3_0=ponderation*reg189;
    T tmp_3_2=ponderation*reg243; T tmp_6_3=ponderation*reg165; T vec_6=-reg305; T tmp_1_0=ponderation*reg242; T tmp_6_4=ponderation*reg169;
    reg152=ponderation*reg81; T vec_5=reg152; T tmp_4_3=ponderation*reg244; T tmp_2_3=ponderation*reg200; T tmp_6_5=ponderation*reg170;
    reg155=ponderation*reg154; T vec_4=reg155; T tmp_4_4=ponderation*reg240; T tmp_2_1=ponderation*reg231; T tmp_6_6=ponderation*reg182;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; T reg2=2*(*f.m).shear_modulus_12; reg1=1.0/reg1; reg0=1.0/reg0;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg5=pow((*f.m).v2[1],2); T reg6=pow((*f.m).v2[0],2);
    T reg7=pow((*f.m).v1[1],2); T reg8=pow((*f.m).v1[0],2); T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg11=reg2*reg3;
    reg5=reg6+reg5; reg6=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg12=reg9*reg11; T reg13=1.0/(*f.m).elastic_modulus_2; T reg14=pow((*f.m).v1[2],2);
    T reg15=reg10*reg11; reg7=reg8+reg7; reg8=reg4*reg11; T reg16=pow((*f.m).v2[2],2); reg16=reg5+reg16;
    reg5=reg13*reg12; T reg17=reg10*reg15; T reg18=reg6*reg12; T reg19=reg10*reg8; reg14=reg7+reg14;
    reg17=reg5-reg17; reg5=1.0/(*f.m).elastic_modulus_1; reg19=reg18+reg19; reg7=reg6*reg15; T reg20=reg13*reg8;
    reg14=pow(reg14,0.5); reg16=pow(reg16,0.5); T reg21=(*f.m).v1[2]/reg14; T reg22=(*f.m).v1[1]/reg14; T reg23=(*f.m).v2[1]/reg16;
    T reg24=(*f.m).v2[2]/reg16; T reg25=reg5*reg17; T reg26=reg6*reg19; T reg27=reg7+reg20; reg16=(*f.m).v2[0]/reg16;
    T reg28=reg6*reg11; T reg29=reg4*reg8; T reg30=reg21*reg23; T reg31=reg22*reg24; T reg32=reg2*reg0;
    T reg33=reg10*reg3; T reg34=reg4*reg3; reg14=(*f.m).v1[0]/reg14; reg3=reg9*reg3; reg11=reg13*reg11;
    reg12=reg5*reg12; T reg35=reg4*reg27; T reg36=reg4*reg15; reg26=reg25-reg26; reg25=reg4*reg28;
    T reg37=reg2*reg1; T reg38=reg4*reg32; reg8=reg6*reg8; T reg39=reg10*reg32; reg15=reg5*reg15;
    reg29=reg12-reg29; reg12=2*reg14; T reg40=reg31-reg30; T reg41=reg6*reg3; T reg42=reg10*reg33;
    T reg43=reg21*reg16; T reg44=reg14*reg24; T reg45=reg10*reg34; reg32=reg9*reg32; T reg46=reg4*reg11;
    reg35=reg26-reg35; reg26=2*reg16; reg36=reg18+reg36; reg3=reg13*reg3; reg18=2*reg40;
    T reg47=reg4*reg37; T reg48=reg26*reg23; T reg49=pow(reg23,2); T reg50=pow(reg16,2); reg17=reg17/reg35;
    T reg51=reg22*reg12; reg36=reg36/reg35; reg46=reg7+reg46; T reg52=reg9*reg37; reg11=reg5*reg11;
    reg19=reg19/reg35; reg29=reg29/reg35; reg34=reg13*reg34; T reg53=reg10*reg38; T reg54=reg10*reg39;
    reg33=reg6*reg33; reg28=reg6*reg28; reg8=reg15+reg8; reg45=reg41+reg45; reg37=reg10*reg37;
    reg25=reg15+reg25; reg42=reg3-reg42; reg3=reg13*reg32; reg15=pow(reg22,2); reg41=pow(reg14,2);
    T reg55=reg43-reg44; reg32=reg6*reg32; T reg56=reg14*reg23; T reg57=reg22*reg16; T reg58=reg33+reg34;
    reg53=reg32+reg53; reg32=reg13*reg52; reg54=reg3-reg54; reg52=reg6*reg52; reg3=reg10*reg37;
    reg45=reg6*reg45; T reg59=reg18*reg55; T reg60=pow(reg55,2); T reg61=pow(reg40,2); T reg62=pow(reg24,2);
    T reg63=reg10*reg47; reg27=reg27/reg35; T reg64=reg15*reg17; T reg65=reg49*reg19; T reg66=reg51*reg17;
    T reg67=reg48*reg19; T reg68=reg56-reg57; T reg69=pow(reg21,2); reg25=reg25/reg35; reg39=reg6*reg39;
    reg8=reg8/reg35; reg28=reg11-reg28; reg11=reg50*reg19; T reg70=reg48*reg29; T reg71=reg41*reg17;
    T reg72=reg51*reg36; reg46=reg46/reg35; reg38=reg13*reg38; T reg73=reg49*reg29; T reg74=reg15*reg36;
    reg42=reg5*reg42; T reg75=reg50*reg29; T reg76=reg41*reg36; T reg77=reg41*reg46; T reg78=reg50*reg25;
    T reg79=reg15*reg46; reg3=reg32-reg3; reg32=pow(reg68,2); reg63=reg52+reg63; reg52=reg61*reg27;
    T reg80=reg59*reg8; reg65=reg64+reg65; reg64=reg60*reg27; T reg81=reg69*reg17; T reg82=reg62*reg19;
    reg70=reg72+reg70; reg67=reg66+reg67; reg66=reg59*reg27; reg72=reg62*reg29; T reg83=reg69*reg36;
    T reg84=reg60*reg8; reg73=reg74+reg73; reg74=reg61*reg8; reg28=reg28/reg35; reg75=reg76+reg75;
    reg45=reg42-reg45; reg58=reg4*reg58; reg42=reg48*reg25; reg54=reg5*reg54; reg76=reg51*reg46;
    reg47=reg13*reg47; reg53=reg6*reg53; reg11=reg71+reg11; reg71=reg39+reg38; reg37=reg6*reg37;
    T reg85=reg49*reg25; T reg86=reg14*reg22; reg82=reg81+reg82; reg81=reg32*reg27; T reg87=reg59*reg28;
    reg71=reg4*reg71; T reg88=reg16*reg23; reg78=reg77+reg78; reg66=reg67+reg66; reg42=reg76+reg42;
    reg67=reg61*reg28; reg76=reg26*reg24; reg11=reg52+reg11; reg52=reg32*reg8; reg72=reg83+reg72;
    reg53=reg54-reg53; reg54=2*reg23; reg84=reg73+reg84; reg73=2*reg22; reg77=reg21*reg12;
    reg3=reg5*reg3; reg85=reg79+reg85; reg58=reg45-reg58; reg63=reg6*reg63; reg74=reg75+reg74;
    reg45=reg60*reg28; reg75=reg37+reg47; reg80=reg70+reg80; reg70=reg69*reg46; reg79=reg62*reg25;
    reg64=reg65+reg64; reg65=reg16*reg12; reg83=reg55*reg40; T reg89=reg73*reg23; T reg90=reg54*reg24;
    T reg91=reg80*reg49; T reg92=reg84*reg88; T reg93=reg86*reg64; reg58=reg58/reg35; T reg94=2*reg21;
    T reg95=reg73*reg21; T reg96=reg22*reg23; reg56=reg57+reg56; reg63=reg3-reg63; reg3=reg22*reg40;
    reg57=reg14*reg55; reg75=reg4*reg75; T reg97=reg14*reg16; T reg98=reg15*reg66; reg71=reg53-reg71;
    reg81=reg82+reg81; reg53=2*reg55; reg18=reg18*reg68; reg82=reg77*reg17; T reg99=reg76*reg19;
    T reg100=reg84*reg49; T reg101=reg15*reg64; T reg102=reg74*reg49; T reg103=reg15*reg11; reg87=reg42+reg87;
    reg42=reg86*reg66; T reg104=reg80*reg88; T reg105=reg32*reg28; reg79=reg70+reg79; reg66=reg41*reg66;
    reg45=reg85+reg45; reg67=reg78+reg67; reg80=reg80*reg50; reg70=reg41*reg11; reg78=reg74*reg50;
    reg52=reg72+reg52; reg72=reg77*reg36; reg85=reg76*reg29; reg64=reg41*reg64; reg84=reg84*reg50;
    T reg106=reg89*reg13; T reg107=reg45*reg60; reg100=reg101+reg100; reg84=reg64+reg84; reg71=reg71/reg35;
    reg64=reg56*reg58; reg101=reg96*reg58; T reg108=reg97*reg58; T reg109=reg67*reg60; reg102=reg103+reg102;
    reg103=reg52*reg50; T reg110=reg41*reg81; reg91=reg98+reg91; reg98=reg87*reg60; T reg111=reg45*reg61;
    T reg112=reg76*reg25; T reg113=reg77*reg46; T reg114=reg89*reg6; T reg115=reg65*reg5; T reg116=reg49*reg6;
    T reg117=reg50*reg5; T reg118=reg67*reg61; reg53=reg53*reg68; reg75=reg63-reg75; reg78=reg70+reg78;
    reg63=reg50*reg6; reg70=reg21*reg24; T reg119=reg49*reg13; T reg120=reg14*reg40; T reg121=reg22*reg55;
    T reg122=reg87*reg61; reg80=reg66+reg80; reg66=reg65*reg6; reg99=reg82+reg99; reg82=reg18*reg27;
    reg17=reg95*reg17; reg19=reg90*reg19; reg57=reg3+reg57; reg3=reg52*reg49; T reg123=reg15*reg81;
    reg104=reg42+reg104; reg36=reg95*reg36; reg87=reg87*reg83; reg74=reg74*reg88; reg29=reg90*reg29;
    reg42=reg23*reg40; T reg124=reg16*reg55; reg85=reg72+reg85; reg72=reg18*reg8; reg92=reg93+reg92;
    reg45=reg45*reg83; reg94=reg24*reg94; reg105=reg79+reg105; reg11=reg86*reg11; reg79=reg108*reg65;
    reg114=reg115-reg114; reg93=reg57*reg71; reg115=reg4*reg94; T reg125=reg121*reg71; T reg126=reg62*reg10;
    reg107=reg100+reg107; reg100=reg101*reg89; reg72=reg85+reg72; reg116=reg117-reg116; reg85=reg64*reg65;
    reg122=reg80+reg122; reg45=reg92+reg45; reg80=reg101*reg56; reg82=reg99+reg82; reg8=reg53*reg8;
    reg92=reg21*reg68; reg87=reg104+reg87; reg99=reg64*reg56; reg29=reg36+reg29; reg124=reg42+reg124;
    reg36=reg105*reg60; reg3=reg123+reg3; reg118=reg78+reg118; reg42=reg62*reg4; reg78=reg16*reg40;
    reg104=reg23*reg55; reg81=reg86*reg81; reg67=reg67*reg83; reg73=reg73*reg55; reg35=reg75/reg35;
    reg74=reg11+reg74; reg25=reg90*reg25; reg46=reg95*reg46; reg98=reg91+reg98; reg11=reg89*reg10;
    reg75=reg18*reg28; reg91=reg4*reg65; reg101=reg101*reg65; reg117=reg49*reg10; reg64=reg64*reg89;
    reg123=reg50*reg4; reg112=reg113+reg112; reg111=reg84+reg111; reg27=reg53*reg27; reg84=reg120*reg71;
    reg63=reg119-reg63; reg113=reg70*reg58; reg19=reg17+reg19; reg17=reg108*reg89; reg109=reg102+reg109;
    reg102=reg94*reg10; reg52=reg52*reg88; reg103=reg110+reg103; reg110=reg105*reg61; reg12=reg40*reg12;
    reg66=reg106-reg66; reg42=reg116-reg42; reg106=reg93*reg12; reg126=reg63-reg126; reg108=reg108*reg56;
    reg63=reg84*reg12; reg115=reg114-reg115; reg114=reg41*reg82; reg67=reg74+reg67; reg79=reg118+reg79;
    reg101=reg111+reg101; reg13=reg15*reg13; reg74=reg72*reg50; reg43=reg44+reg43; reg102=reg66-reg102;
    reg44=reg41*reg6; reg66=reg125*reg12; reg105=reg105*reg83; reg52=reg81+reg52; reg6=reg15*reg6;
    reg27=reg19+reg27; reg19=(*f.m).alpha_2*reg49; reg81=(*f.m).alpha_1*reg15; reg111=reg124*reg35; reg116=reg104*reg35;
    reg118=reg78*reg35; reg80=reg45+reg80; reg45=reg92*reg71; reg119=reg84*reg73; reg17=reg109+reg17;
    reg109=(*f.m).alpha_2*reg50; T reg127=reg93*reg73; T reg128=reg113*reg65; reg110=reg103+reg110; reg26=reg26*reg40;
    reg103=(*f.m).alpha_1*reg41; reg28=reg53*reg28; reg25=reg46+reg25; reg54=reg54*reg55; reg75=reg112+reg75;
    reg11=reg91+reg11; reg94=reg94*reg9; reg117=reg123+reg117; reg46=reg62*reg9; reg64=reg98+reg64;
    reg85=reg122+reg85; reg99=reg87+reg99; reg93=reg93*reg57; reg5=reg41*reg5; reg87=reg72*reg49;
    reg91=reg21*reg40; reg98=reg14*reg68; reg112=reg113*reg89; reg36=reg3+reg36; reg100=reg107+reg100;
    reg3=reg125*reg73; reg125=reg125*reg57; reg107=reg24*reg68; reg122=reg15*reg82; reg8=reg29+reg8;
    reg108=reg67+reg108; reg84=reg84*reg57; reg29=reg75*reg60; reg31=reg30+reg31; reg30=reg45*reg73;
    reg87=reg122+reg87; reg67=reg21*reg55; reg122=reg111*reg54; reg123=reg27*reg15; T reg129=reg118*reg26;
    reg63=reg79+reg63; reg79=reg8*reg49; reg125=reg80+reg125; reg80=reg43*reg58; T reg130=reg107*reg35;
    T reg131=reg22*reg68; reg28=reg25+reg28; reg127=reg64+reg127; reg98=reg91+reg98; reg25=reg116*reg54;
    reg3=reg100+reg3; reg113=reg113*reg56; reg64=reg102*reg15; reg91=reg115*reg41; reg100=reg118*reg54;
    reg119=reg17+reg119; reg17=reg126*reg15; T reg132=reg42*reg41; T reg133=reg45*reg12; reg128=reg110+reg128;
    reg110=reg116*reg26; reg66=reg101+reg66; reg11=reg94-reg11; reg117=reg46-reg117; reg46=reg15*reg10;
    reg94=reg42*reg50; reg101=reg126*reg49; T reg134=reg126*reg96; T reg135=reg42*reg97; T reg136=reg115*reg50;
    T reg137=reg102*reg49; T reg138=reg115*reg97; T reg139=reg24*reg40; T reg140=reg16*reg68; T reg141=reg102*reg96;
    T reg142=(*f.m).alpha_2*reg62; T reg143=(*f.m).alpha_1*reg69; T reg144=(*f.m).alpha_3*reg60; reg19=reg81+reg19; reg81=(*f.m).alpha_3*reg61;
    reg109=reg103+reg109; reg103=reg41*reg4; reg105=reg52+reg105; reg4=reg69*reg4; reg6=reg5-reg6;
    reg5=reg75*reg61; reg74=reg114+reg74; reg52=reg111*reg26; reg106=reg85+reg106; reg112=reg36+reg112;
    reg10=reg69*reg10; reg93=reg99+reg93; reg111=reg111*reg124; reg82=reg86*reg82; reg44=reg13-reg44;
    reg13=reg27*reg41; reg116=reg116*reg124; reg72=reg72*reg88; reg36=reg8*reg50; reg116=reg125+reg116;
    reg144=reg19+reg144; reg122=reg127+reg122; reg101=reg94+reg101; reg4=reg6-reg4; reg6=reg117*reg62;
    reg142=reg143+reg142; reg19=reg117*reg70; reg85=(*f.m).alpha_3*reg32; reg94=(*f.m).alpha_1*reg86; reg29=reg87+reg29;
    reg134=reg135+reg134; reg87=reg80*reg89; reg118=reg118*reg124; reg84=reg108+reg84; reg99=(*f.m).alpha_2*reg88;
    reg141=reg138+reg141; reg108=reg28*reg61; reg36=reg13+reg36; reg129=reg63+reg129; reg137=reg136+reg137;
    reg140=reg139+reg140; reg13=reg11*reg62; reg10=reg44-reg10; reg79=reg123+reg79; reg44=reg28*reg60;
    reg17=reg132+reg17; reg63=reg117*reg69; reg133=reg128+reg133; reg58=reg31*reg58; reg114=reg130*reg26;
    reg72=reg82+reg72; reg100=reg119+reg100; reg82=reg98*reg71; reg75=reg75*reg83; reg25=reg3+reg25;
    reg64=reg91+reg64; reg3=reg11*reg69; reg131=reg67+reg131; reg27=reg27*reg86; reg67=reg130*reg54;
    reg30=reg112+reg30; reg8=reg8*reg88; reg81=reg109+reg81; reg91=reg80*reg65; reg5=reg74+reg5;
    reg74=reg24*reg55; reg109=reg23*reg68; reg112=reg11*reg70; reg45=reg45*reg57; reg9=reg69*reg9;
    reg46=reg103+reg46; reg52=reg106+reg52; reg113=reg105+reg113; reg88=reg88*reg2; reg103=reg56*reg2;
    reg14=reg14*reg21; reg16=reg16*reg24; reg110=reg66+reg110; reg111=reg93+reg111; reg66=reg144*reg25;
    reg93=reg81*reg100; reg89=reg58*reg89; reg44=reg79+reg44; reg118=reg84+reg118; reg45=reg113+reg45;
    reg130=reg130*reg124; reg79=reg4*reg50; reg84=reg10*reg49; reg6=reg101+reg6; reg101=reg88*reg48;
    reg13=reg137+reg13; reg105=reg103*reg48; reg106=(*f.m).alpha_2*reg16; reg113=(*f.m).alpha_1*reg14; reg119=(*f.m).alpha_3*reg83;
    reg99=reg94+reg99; reg85=reg142+reg85; reg109=reg74+reg109; reg40=reg68*reg40; reg21=reg22*reg21;
    reg24=reg23*reg24; reg22=reg129*reg81; reg3=reg64+reg3; reg23=reg103*reg51; reg83=reg28*reg83;
    reg8=reg27+reg8; reg65=reg58*reg65; reg71=reg131*reg71; reg27=reg140*reg35; reg108=reg36+reg108;
    reg28=reg82*reg12; reg91=reg5+reg91; reg114=reg133+reg114; reg80=reg80*reg56; reg75=reg72+reg75;
    reg5=reg103*reg56; reg112=reg141+reg112; reg36=reg111*reg25; reg64=reg111*reg110; reg72=reg116*reg122;
    reg74=reg116*reg52; reg19=reg134+reg19; reg94=reg88*reg56; reg46=reg9-reg46; reg16=reg16*reg1;
    reg9=reg43*reg1; reg123=reg4*reg41; reg125=reg10*reg15; reg127=reg82*reg73; reg63=reg17+reg63;
    reg17=reg88*reg51; reg128=reg144*reg110; reg87=reg29+reg87; reg67=reg30+reg67; reg105=reg13+reg105;
    reg13=reg31*reg0; reg55=reg68*reg55; reg29=reg24*reg0; reg58=reg58*reg56; reg72=reg36-reg72;
    reg40=(*f.m).alpha_3*reg40; reg30=reg16*reg76; reg74=reg64-reg74; reg36=reg110*reg122; reg101=reg6+reg101;
    reg6=reg52*reg25; reg2=reg86*reg2; reg64=reg16*reg77; reg68=reg9*reg76; reg86=reg10*reg96;
    reg132=reg4*reg97; reg23=reg3+reg23; reg125=reg123+reg125; reg69=reg46*reg69; reg119=reg99+reg119;
    reg80=reg75+reg80; reg82=reg82*reg57; reg24=(*f.m).alpha_2*reg24; reg3=reg9*reg77; reg75=(*f.m).alpha_1*reg21;
    reg106=reg113+reg106; reg83=reg8+reg83; reg17=reg63+reg17; reg35=reg109*reg35; reg65=reg108+reg65;
    reg12=reg71*reg12; reg8=reg27*reg26; reg28=reg91+reg28; reg63=reg144*reg116; reg91=reg81*reg118;
    reg99=reg9*reg43; reg128=reg22+reg128; reg114=reg85*reg114; reg130=reg45+reg130; reg5=reg112+reg5;
    reg127=reg87+reg127; reg22=reg27*reg54; reg89=reg44+reg89; reg94=reg19+reg94; reg19=reg16*reg43;
    reg73=reg71*reg73; reg84=reg79+reg84; reg67=reg85*reg67; reg62=reg46*reg62; reg66=reg93+reg66;
    reg12=reg65+reg12; reg1=reg14*reg1; reg40=reg106+reg40; reg14=reg13*reg95; reg3=reg23+reg3;
    reg26=reg35*reg26; reg130=reg85*reg130; reg63=reg91+reg63; reg73=reg89+reg73; reg23=reg29*reg95;
    reg64=reg17+reg64; reg114=reg128+reg114; reg67=reg66+reg67; reg17=reg119*reg122; reg44=reg119*reg52;
    reg54=reg35*reg54; reg71=reg71*reg57; reg58=reg83+reg58; reg45=reg2*reg51; reg69=reg125+reg69;
    reg22=reg127+reg22; reg65=reg29*reg31; reg30=reg101+reg30; reg66=reg29*reg90; reg19=reg94+reg19;
    reg79=reg129*reg72; reg83=reg74*reg100; reg70=reg46*reg70; reg62=reg84+reg62; reg68=reg105+reg68;
    reg84=reg13*reg90; reg87=reg2*reg48; reg86=reg132+reg86; reg99=reg5+reg99; reg6=reg36-reg6;
    reg5=reg13*reg31; reg55=(*f.m).alpha_3*reg55; reg24=reg75+reg24; reg8=reg28+reg8; reg82=reg80+reg82;
    reg27=reg27*reg124; reg87=reg62+reg87; reg76=reg1*reg76; reg28=1-var_inter[1]; reg54=reg73+reg54;
    reg65=reg19+reg65; reg66=reg30+reg66; reg19=reg2*reg56; reg70=reg86+reg70; reg71=reg58+reg71;
    reg84=reg68+reg84; reg26=reg12+reg26; reg35=reg35*reg124; reg55=reg24+reg55; reg5=reg99+reg5;
    reg12=reg6*reg118; reg24=reg111*reg100; reg83=reg79-reg83; reg30=reg118*reg122; reg77=reg1*reg77;
    reg36=reg129*reg111; reg23=reg64+reg23; reg14=reg3+reg14; reg3=reg118*reg52; reg27=reg82+reg27;
    reg0=reg21*reg0; reg45=reg69+reg45; reg21=1-var_inter[0]; reg8=reg40*reg8; reg17=reg67+reg17;
    reg22=reg40*reg22; reg44=reg114+reg44; reg130=reg63+reg130; reg58=reg119*reg111; reg19=reg70+reg19;
    reg43=reg1*reg43; reg12=reg83+reg12; reg35=reg71+reg35; reg62=reg14*reg65; reg90=reg0*reg90;
    reg76=reg87+reg76; reg63=reg116*reg100; reg30=reg24-reg30; reg24=reg118*reg25; reg64=reg129*reg116;
    reg3=reg36-reg3; reg36=reg118*reg110; reg67=reg129*reg122; reg68=reg52*reg100; reg27=reg40*reg27;
    reg58=reg130+reg58; reg54=reg55*reg54; reg22=reg17+reg22; reg8=reg44+reg8; reg26=reg55*reg26;
    reg17=reg21*elem.pos(0)[0]; reg44=elem.pos(1)[0]*var_inter[0]; reg69=reg28*elem.pos(1)[1]; reg70=reg28*elem.pos(0)[1]; reg71=reg28*elem.pos(1)[0];
    reg77=reg45+reg77; reg45=reg28*elem.pos(0)[0]; reg95=reg0*reg95; reg73=reg23*reg5; reg75=reg84*reg65;
    reg79=elem.pos(0)[1]*reg21; reg80=reg66*reg5; reg82=elem.pos(1)[1]*var_inter[0]; reg90=reg76+reg90; reg76=var_inter[1]*elem.pos(2)[0];
    reg83=reg17+reg44; reg27=reg58+reg27; reg35=reg55*reg35; reg68=reg67-reg68; reg58=reg129*reg25;
    reg67=reg14*reg66; reg69=reg69-reg70; reg54=reg22+reg54; reg24=reg63-reg24; reg74=reg74/reg12;
    reg22=var_inter[1]*elem.pos(2)[1]; reg3=reg3/reg12; reg36=reg64-reg36; reg63=reg110*reg100; reg26=reg8+reg26;
    reg75=reg80-reg75; reg95=reg77+reg95; reg43=reg19+reg43; reg8=elem.pos(2)[1]*var_inter[0]; reg19=reg79+reg82;
    reg62=reg73-reg62; reg31=reg0*reg31; reg64=elem.pos(2)[0]*var_inter[0]; reg73=reg23*reg84; reg72=reg72/reg12;
    reg71=reg71-reg45; reg30=reg30/reg12; reg35=reg27+reg35; reg64=reg64-reg83; reg27=reg21*elem.pos(3)[0];
    reg68=reg68/reg12; reg63=reg58-reg63; reg6=reg6/reg12; reg58=reg62*reg90; reg77=reg95*reg75;
    reg80=reg21*elem.pos(3)[1]; reg8=reg8-reg19; reg31=reg43+reg31; reg76=reg71+reg76; reg43=var_inter[1]*elem.pos(3)[0];
    reg30=reg30*reg26; reg74=reg74*reg54; reg72=reg72*reg26; reg3=reg3*reg54; reg67=reg73-reg67;
    reg24=reg24/reg12; reg22=reg69+reg22; reg69=var_inter[1]*elem.pos(3)[1]; reg36=reg36/reg12; reg80=reg8+reg80;
    reg54=reg36*reg54; reg68=reg68*reg35; reg26=reg24*reg26; reg30=reg3-reg30; reg27=reg64+reg27;
    reg12=reg63/reg12; reg22=reg22-reg69; reg76=reg76-reg43; reg6=reg6*reg35; reg3=reg23*reg31;
    reg8=reg95*reg65; reg24=reg66*reg31; reg36=reg90*reg65; reg74=reg72-reg74; reg58=reg77-reg58;
    reg63=reg67*reg31; reg64=reg84*reg31; reg71=reg95*reg5; reg24=reg36-reg24; reg36=reg14*reg31;
    reg3=reg8-reg3; reg8=reg95*reg66; reg72=reg76*reg80; reg73=reg23*reg90; reg77=1-(*f.m).resolution;
    reg6=reg74+reg6; reg74=reg90*reg5; reg68=reg30-reg68; reg63=reg58+reg63; reg35=reg12*reg35;
    reg12=reg22*reg27; reg54=reg26-reg54; reg64=reg74-reg64; reg24=reg24/reg63; reg36=reg71-reg36;
    reg3=reg3/reg63; reg12=reg72-reg12; reg26=reg95*reg84; reg30=reg14*reg90; reg68=reg68*reg77;
    reg73=reg8-reg73; reg8=(*f.m).resolution*reg81; reg58=(*f.m).resolution*reg144; reg6=reg6*reg77; reg54=reg35+reg54;
    reg35=(*f.m).resolution*reg119; reg73=reg73/reg63; reg30=reg26-reg30; reg76=reg76/reg12; reg26=(*f.m).resolution*reg24;
    reg71=(*f.m).resolution*reg3; reg64=reg64/reg63; reg116=reg77*reg116; reg75=reg75/reg63; reg36=reg36/reg63;
    reg118=reg77*reg118; reg62=reg62/reg63; reg68=reg58+reg68; reg54=reg54*reg77; reg6=reg8+reg6;
    reg80=reg80/reg12; reg27=reg27/reg12; reg22=reg22/reg12; reg25=reg77*reg25; reg111=reg77*reg111;
    reg8=(*f.m).resolution*reg64; reg100=reg77*reg100; reg58=(*f.m).resolution*reg36; reg68=reg68*(*f.m).deltaT; reg72=(*f.m).resolution*reg73;
    reg6=reg6*(*f.m).deltaT; reg54=reg35+reg54; reg26=reg118+reg26; reg71=reg116-reg71; reg35=reg21*reg22;
    reg67=reg67/reg63; reg74=var_inter[0]*reg76; reg86=reg28*reg80; reg63=reg30/reg63; reg30=reg28*reg27;
    reg87=reg21*reg76; reg89=(*f.m).resolution*reg62; reg91=var_inter[0]*reg22; reg93=(*f.m).resolution*reg75; reg94=var_inter[1]*reg80;
    reg129=reg129*reg77; reg110=reg77*reg110; reg99=var_inter[1]*reg27; reg101=reg30+reg74; reg105=reg86+reg91;
    reg106=reg94-reg91; reg72=reg111+reg72; reg108=reg74-reg99; reg111=reg35+reg94; reg58=reg25+reg58;
    reg8=reg100-reg8; reg25=reg87+reg99; reg100=reg35-reg86; reg89=reg110-reg89; reg129=reg93+reg129;
    reg93=reg30-reg87; reg110=(*f.m).resolution*reg63; reg54=reg54*(*f.m).deltaT; reg112=(*f.m).resolution*reg67; reg113=reg68*reg71;
    reg114=reg6*reg26; reg122=reg77*reg122; reg52=reg77*reg52; reg77=reg6*reg8; reg116=0.5*reg25;
    reg118=reg54*reg72; reg123=reg114+reg113; reg125=reg68*reg58; reg127=0.5*reg100; reg128=reg6*reg129;
    reg130=reg68*reg89; reg132=0.5*reg93; reg133=0.5*reg106; reg134=0.5*reg111; reg135=0.5*reg108;
    reg52=reg112+reg52; reg110=reg122-reg110; reg112=0.5*reg105; reg122=0.5*reg101; reg136=reg116*reg72;
    reg137=reg122*reg72; reg138=reg105*reg26; reg139=reg111*reg26; reg141=reg101*reg71; reg142=reg112*reg72;
    reg143=reg106*reg26; T reg145=reg135*reg72; T reg146=reg127*reg72; T reg147=reg93*reg71; T reg148=reg133*reg72;
    T reg149=reg108*reg71; T reg150=reg132*reg72; T reg151=reg100*reg26; T reg152=reg54*reg110; T reg153=reg77+reg125;
    T reg154=reg25*reg71; T reg155=reg123+reg118; T reg156=reg134*reg72; T reg157=reg54*reg52; T reg158=reg130+reg128;
    T reg159=reg105*reg8; T reg160=reg122*reg110; T reg161=reg127*reg110; reg146=reg147+reg146; reg148=reg149+reg148;
    reg147=reg93*reg58; reg149=reg101*reg58; T reg162=reg112*reg110; T reg163=reg108*reg58; reg154=reg154-reg156;
    T reg164=reg111*reg129; T reg165=reg116*reg52; reg136=reg136-reg139; T reg166=reg134*reg52; T reg167=reg25*reg89;
    T reg168=reg135*reg110; T reg169=reg106*reg8; T reg170=reg134*reg110; reg150=reg151+reg150; reg151=reg25*reg58;
    T reg171=reg132*reg52; T reg172=reg100*reg129; T reg173=reg116*reg110; T reg174=reg111*reg8; T reg175=reg158+reg157;
    T reg176=2*reg155; T reg177=reg153+reg152; T reg178=reg133*reg110; T reg179=reg105*reg129; T reg180=reg122*reg52;
    reg138=reg138-reg137; T reg181=reg101*reg89; T reg182=reg127*reg52; T reg183=reg133*reg52; T reg184=reg108*reg89;
    reg145=reg143+reg145; reg143=reg93*reg89; T reg185=reg135*reg52; T reg186=reg106*reg129; reg142=reg142-reg141;
    T reg187=reg112*reg52; reg142=2*reg142; reg178=reg163+reg178; reg151=reg151-reg170; reg185=reg186+reg185;
    reg148=2*reg148; reg179=reg179-reg180; reg162=reg162-reg149; reg187=reg187-reg181; reg145=2*reg145;
    reg168=reg169+reg168; reg159=reg159-reg160; reg173=reg173-reg174; reg138=2*reg138; reg183=reg184+reg183;
    reg182=reg143+reg182; reg146=2*reg146; reg165=reg165-reg164; reg143=reg175*reg105; reg163=reg176*reg116;
    reg136=2*reg136; reg169=reg175*reg111; reg184=reg177*reg25; reg186=reg176*reg134; reg167=reg167-reg166;
    T reg188=reg176*reg122; reg154=2*reg154; reg171=reg172+reg171; reg161=reg147+reg161; reg147=reg176*reg112;
    reg150=2*reg150; reg172=reg177*reg101; T reg189=reg101*reg178; T reg190=reg148*reg112; T reg191=reg101*reg168;
    T reg192=reg145*reg112; T reg193=reg101*reg162; T reg194=reg142*reg112; T reg195=reg122*reg154; T reg196=reg105*reg167;
    T reg197=reg179*reg100; T reg198=reg122*reg136; T reg199=reg105*reg165; T reg200=reg148*reg122; T reg201=reg183*reg105;
    T reg202=reg175*reg106; T reg203=reg176*reg135; T reg204=reg112*reg136; T reg205=reg132*reg146; T reg206=reg101*reg173;
    T reg207=reg112*reg154; T reg208=reg93*reg178; T reg209=reg135*reg154; T reg210=reg106*reg167; T reg211=reg151*reg101;
    T reg212=reg185*reg106; T reg213=reg145*reg135; T reg214=reg183*reg106; T reg215=reg106*reg165; T reg216=reg135*reg136;
    T reg217=reg100*reg171; T reg218=reg127*reg136; T reg219=reg93*reg173; T reg220=reg108*reg178; T reg221=reg148*reg133;
    T reg222=reg108*reg173; T reg223=reg133*reg136; T reg224=reg151*reg108; T reg225=reg133*reg154; T reg226=reg116*reg136;
    T reg227=reg111*reg165; T reg228=reg116*reg154; T reg229=reg111*reg167; T reg230=reg151*reg25; T reg231=reg134*reg154;
    T reg232=reg186-reg184; T reg233=reg175*reg100; T reg234=reg169-reg163; T reg235=reg176*reg132; T reg236=reg176*reg133;
    T reg237=reg177*reg93; T reg238=reg176*reg127; T reg239=reg148*reg127; T reg240=reg177*reg108; reg151=reg151*reg93;
    T reg241=reg127*reg154; T reg242=reg179*reg105; T reg243=reg138*reg122; T reg244=reg187*reg105; T reg245=reg142*reg122;
    T reg246=reg185*reg105; T reg247=reg145*reg122; T reg248=reg138*reg127; T reg249=reg145*reg132; T reg250=reg185*reg100;
    T reg251=reg93*reg159; T reg252=reg132*reg150; T reg253=reg148*reg135; T reg254=reg183*reg100; T reg255=reg100*reg182;
    T reg256=reg142*reg132; T reg257=reg148*reg132; T reg258=reg188-reg143; T reg259=reg93*reg161; T reg260=reg145*reg127;
    T reg261=reg187*reg100; T reg262=reg93*reg168; reg165=reg100*reg165; T reg263=reg93*reg162; T reg264=reg138*reg132;
    reg154=reg132*reg154; T reg265=reg142*reg127; T reg266=reg172-reg147; reg167=reg100*reg167; T reg267=reg132*reg136;
    T reg268=reg127*reg146; reg198=reg199-reg198; reg249=reg250+reg249; reg223=reg222+reg223; reg248=reg251+reg248;
    reg229=reg228-reg229; reg200=reg201-reg200; reg253=reg214+reg253; reg247=reg246-reg247; reg195=reg196-reg195;
    reg231=reg230-reg231; reg245=reg244-reg245; reg232=reg12*reg232; reg243=reg242-reg243; reg257=reg254+reg257;
    reg265=reg263+reg265; reg225=reg224+reg225; reg241=reg151+reg241; reg151=reg237+reg238; reg260=reg262+reg260;
    reg239=reg208+reg239; reg234=reg12*reg234; reg267=reg165+reg267; reg165=reg240+reg236; reg154=reg167+reg154;
    reg167=reg233+reg235; reg227=reg226-reg227; reg206=reg204-reg206; reg218=reg219+reg218; reg209=reg210+reg209;
    reg256=reg261+reg256; reg211=reg207-reg211; reg264=reg197+reg264; reg216=reg215+reg216; reg189=reg190-reg189;
    reg252=reg217+reg252; reg258=reg12*reg258; reg255=reg205+reg255; reg191=reg192-reg191; reg213=reg212+reg213;
    reg221=reg220+reg221; reg268=reg259+reg268; reg190=reg202+reg203; reg266=reg12*reg266; reg193=reg194-reg193;
    reg213=reg213*reg12; reg258=ponderation*reg258; reg192=reg12*reg165; reg216=reg216*reg12; reg268=reg12*reg268;
    reg194=reg151*reg12; reg154=reg12*reg154; reg196=reg167*reg12; reg264=reg264*reg12; reg209=reg209*reg12;
    reg267=reg12*reg267; reg218=reg12*reg218; reg234=ponderation*reg234; reg239=reg12*reg239; reg227=reg227*reg12;
    reg265=reg12*reg265; reg197=reg12*reg190; reg255=reg255*reg12; reg231=reg231*reg12; reg266=ponderation*reg266;
    reg229=reg229*reg12; reg221=reg221*reg12; reg195=reg195*reg12; reg193=reg193*reg12; reg253=reg253*reg12;
    reg252=reg12*reg252; reg198=reg198*reg12; reg191=reg191*reg12; reg249=reg249*reg12; reg200=reg200*reg12;
    reg223=reg223*reg12; reg189=reg189*reg12; reg247=reg247*reg12; reg256=reg256*reg12; reg260=reg12*reg260;
    reg225=reg225*reg12; reg257=reg12*reg257; reg241=reg241*reg12; reg211=reg211*reg12; reg243=reg243*reg12;
    reg232=ponderation*reg232; reg206=reg206*reg12; reg245=reg245*reg12; reg248=reg12*reg248; T tmp_1_2=ponderation*reg248;
    T tmp_6_6=ponderation*reg227; T tmp_1_6=ponderation*reg218; T tmp_5_6=ponderation*reg223; T vec_7=-reg232; T tmp_1_1=ponderation*reg268;
    T tmp_0_7=ponderation*reg154; T tmp_1_5=ponderation*reg239; T tmp_5_7=ponderation*reg225; T tmp_5_5=ponderation*reg221; T tmp_2_7=ponderation*reg195;
    T tmp_4_5=ponderation*reg253; T tmp_0_0=ponderation*reg252; T tmp_3_3=ponderation*reg193; T tmp_3_4=ponderation*reg191; T tmp_3_5=ponderation*reg189;
    T tmp_0_3=ponderation*reg256; T tmp_3_6=ponderation*reg206; T tmp_3_7=ponderation*reg211; T tmp_4_4=ponderation*reg213; T vec_2=-reg258;
    T tmp_4_6=ponderation*reg216; T tmp_0_2=ponderation*reg264; T tmp_4_7=ponderation*reg209; reg154=ponderation*reg197; T vec_4=reg154;
    T tmp_0_1=ponderation*reg255; T vec_3=-reg266; T tmp_6_7=ponderation*reg229; T tmp_7_7=ponderation*reg231; T tmp_1_3=ponderation*reg265;
    T vec_6=-reg234; T tmp_0_6=ponderation*reg267; reg189=ponderation*reg196; T vec_0=reg189; reg191=ponderation*reg192;
    T vec_5=reg191; reg193=ponderation*reg194; T vec_1=reg193; T tmp_0_5=ponderation*reg257; T tmp_1_4=ponderation*reg260;
    T tmp_1_7=ponderation*reg241; T tmp_2_2=ponderation*reg243; T tmp_2_3=ponderation*reg245; T tmp_2_4=ponderation*reg247; T tmp_0_4=ponderation*reg249;
    T tmp_2_5=ponderation*reg200; T tmp_2_6=ponderation*reg198;
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
    T reg3=pow((*f.m).v2[1],2); reg0=reg1+reg0; reg3=reg2+reg3; reg1=pow((*f.m).v2[2],2); reg2=2*(*f.m).shear_modulus_13;
    T reg4=2*(*f.m).shear_modulus_23; reg1=reg3+reg1; reg4=1.0/reg4; reg0=pow(reg0,0.5); reg3=2*(*f.m).shear_modulus_12;
    reg2=1.0/reg2; T reg5=reg2*reg4; T reg6=(*f.m).v1[0]/reg0; T reg7=(*f.m).v1[1]/reg0; reg1=pow(reg1,0.5);
    reg3=1.0/reg3; reg0=(*f.m).v1[2]/reg0; T reg8=(*f.m).v2[0]/reg1; T reg9=(*f.m).v2[1]/reg1; T reg10=2*reg7;
    T reg11=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg12=2*reg6; T reg13=reg3*reg5; T reg14=1.0/(*f.m).elastic_modulus_3; T reg15=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg16=reg11*reg13; T reg17=reg15*reg13; T reg18=pow(reg9,2); T reg19=pow(reg8,2); T reg20=reg14*reg13;
    T reg21=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg22=1.0/(*f.m).elastic_modulus_2; T reg23=1.0/(*f.m).elastic_modulus_1; reg1=(*f.m).v2[2]/reg1; T reg24=2*reg0;
    T reg25=reg10*reg9; T reg26=reg8*reg12; T reg27=reg26*reg21; T reg28=reg25*reg22; T reg29=reg18*reg22;
    T reg30=reg19*reg21; T reg31=reg22*reg20; T reg32=pow(reg1,2); T reg33=reg15*reg17; T reg34=reg21*reg20;
    T reg35=reg25*reg21; T reg36=reg15*reg16; T reg37=reg26*reg23; reg24=reg1*reg24; T reg38=reg19*reg23;
    T reg39=reg18*reg21; T reg40=reg18*reg15; T reg41=reg19*reg11; T reg42=reg32*reg15; T reg43=reg25*reg15;
    T reg44=reg11*reg26; reg30=reg29-reg30; reg29=pow(reg6,2); T reg45=reg24*reg15; reg33=reg31-reg33;
    reg31=reg11*reg24; reg39=reg38-reg39; reg36=reg34+reg36; reg38=reg21*reg17; T reg46=reg32*reg11;
    T reg47=reg22*reg16; reg35=reg37-reg35; reg27=reg28-reg27; reg28=pow(reg7,2); reg37=reg29*reg21;
    reg31=reg35-reg31; reg35=reg28*reg21; reg24=reg24*reg14; reg43=reg44+reg43; reg40=reg41+reg40;
    reg41=reg29*reg23; reg46=reg39-reg46; reg39=reg28*reg22; reg44=reg32*reg14; T reg48=reg7*reg8;
    reg42=reg30-reg42; reg45=reg27-reg45; reg27=reg23*reg33; reg30=reg21*reg36; T reg49=reg6*reg9;
    T reg50=reg6*reg8; T reg51=reg7*reg9; T reg52=reg38+reg47; T reg53=pow(reg0,2); T reg54=reg11*reg17;
    T reg55=reg46*reg29; T reg56=reg42*reg28; T reg57=reg31*reg29; T reg58=reg45*reg28; T reg59=reg3*reg4;
    T reg60=reg21*reg13; T reg61=reg11*reg16; T reg62=reg45*reg18; T reg63=reg15*reg5; T reg64=reg11*reg5;
    reg13=reg22*reg13; reg43=reg24-reg43; reg5=reg14*reg5; reg24=reg11*reg52; reg30=reg27-reg30;
    reg20=reg23*reg20; reg27=reg8*reg9; T reg65=reg29*reg11; T reg66=reg28*reg15; reg40=reg44-reg40;
    reg44=reg7*reg1; T reg67=reg0*reg9; T reg68=reg53*reg11; reg35=reg41-reg35; reg41=reg53*reg15;
    T reg69=reg0*reg8; T reg70=reg46*reg50; T reg71=reg42*reg51; T reg72=reg6*reg1; T reg73=reg46*reg19;
    T reg74=reg42*reg18; T reg75=reg48+reg49; reg37=reg39-reg37; reg39=reg31*reg19; T reg76=reg0*reg1;
    T reg77=reg31*reg50; T reg78=reg45*reg51; T reg79=2*reg8; T reg80=reg79*reg9; T reg81=reg27*reg3;
    T reg82=reg22*reg5; T reg83=reg44-reg67; T reg84=reg7*reg12; reg68=reg35-reg68; reg24=reg30-reg24;
    reg30=reg14*reg59; reg35=reg15*reg64; T reg85=reg43*reg76; reg71=reg70+reg71; reg78=reg77+reg78;
    reg70=reg40*reg76; reg66=reg65+reg66; reg65=reg53*reg14; reg77=reg72+reg69; T reg86=reg8*reg1;
    reg5=reg21*reg5; reg41=reg37-reg41; reg37=reg15*reg63; reg62=reg39+reg62; reg39=reg43*reg32;
    reg61=reg20-reg61; reg17=reg23*reg17; reg20=reg40*reg32; T reg87=reg11*reg60; reg54=reg34+reg54;
    reg34=reg15*reg59; T reg88=reg43*reg53; reg58=reg57+reg58; reg74=reg73+reg74; reg16=reg21*reg16;
    reg57=reg75*reg3; reg73=reg40*reg53; reg56=reg55+reg56; reg55=reg11*reg13; reg59=reg11*reg59;
    T reg89=reg3*reg2; T reg90=2*reg83; reg35=reg5+reg35; reg33=reg33/reg24; reg61=reg61/reg24;
    reg36=reg36/reg24; reg5=reg22*reg30; reg30=reg21*reg30; T reg91=reg15*reg34; T reg92=2*reg9;
    reg13=reg23*reg13; reg37=reg82-reg37; reg64=reg22*reg64; reg82=reg11*reg89; reg55=reg38+reg55;
    T reg93=reg15*reg89; reg63=reg21*reg63; reg87=reg17+reg87; reg72=reg69-reg72; reg60=reg21*reg60;
    reg44=reg67+reg44; reg54=reg54/reg24; reg16=reg17+reg16; reg17=reg81*reg75; reg66=reg65-reg66;
    reg70=reg71+reg70; reg65=reg68*reg19; reg67=reg41*reg18; reg69=reg86*reg2; reg71=reg77*reg2;
    T reg94=reg68*reg29; T reg95=reg41*reg28; reg73=reg56+reg73; reg56=reg81*reg84; reg88=reg58+reg88;
    reg58=reg57*reg84; reg20=reg74+reg20; reg74=reg81*reg80; T reg96=reg57*reg80; reg39=reg62+reg39;
    reg62=reg79*reg1; T reg97=reg15*reg59; T reg98=reg0*reg12; reg89=reg14*reg89; reg14=reg57*reg75;
    T reg99=reg6*reg7; T reg100=reg9*reg1; reg85=reg78+reg85; reg78=reg69*reg62; T reg101=reg6*reg0;
    reg96=reg39+reg96; reg39=reg71*reg62; reg74=reg20+reg74; reg20=reg10*reg0; T reg102=reg29*reg33;
    T reg103=reg71*reg98; reg58=reg88+reg58; reg88=reg69*reg98; reg56=reg73+reg56; reg73=reg92*reg1;
    T reg104=reg90*reg72; T reg105=reg66*reg32; reg67=reg65+reg67; reg65=reg66*reg53; reg95=reg94+reg95;
    reg94=pow(reg83,2); reg52=reg52/reg24; T reg106=reg44*reg4; T reg107=reg100*reg4; T reg108=pow(reg72,2);
    T reg109=1-var_inter[0]; reg3=reg99*reg3; reg55=reg55/reg24; T reg110=reg71*reg77; reg97=reg30+reg97;
    reg14=reg85+reg14; reg30=reg22*reg89; reg89=reg21*reg89; reg85=reg80*reg36; T reg111=reg69*reg77;
    T reg112=reg15*reg93; T reg113=reg84*reg33; reg17=reg70+reg17; reg70=reg18*reg36; T reg114=reg29*reg54;
    T reg115=reg19*reg61; reg15=reg15*reg82; T reg116=reg28*reg54; T reg117=reg18*reg61; T reg118=reg28*reg33;
    T reg119=reg84*reg54; T reg120=reg80*reg61; T reg121=reg68*reg50; T reg122=1-var_inter[1]; reg87=reg87/reg24;
    reg16=reg16/reg24; reg60=reg13-reg60; reg34=reg21*reg34; reg59=reg22*reg59; reg13=reg19*reg36;
    T reg123=reg41*reg51; reg37=reg23*reg37; reg91=reg5-reg91; reg35=reg21*reg35; reg5=reg63+reg64;
    reg65=reg95+reg65; reg95=elem.pos(0)[1]*reg109; reg78=reg74+reg78; reg39=reg96+reg39; reg74=reg28*reg55;
    reg96=reg19*reg87; T reg124=elem.pos(1)[0]*var_inter[0]; T reg125=reg106*reg20; T reg126=reg109*elem.pos(0)[0]; T reg127=reg29*reg55;
    T reg128=reg106*reg73; T reg129=reg104*reg16; T reg130=reg122*elem.pos(1)[1]; T reg131=reg122*elem.pos(0)[1]; T reg132=elem.pos(1)[1]*var_inter[0];
    reg103=reg58+reg103; reg58=reg107*reg73; T reg133=reg18*reg87; T reg134=reg107*reg20; T reg135=reg84*reg55;
    T reg136=reg80*reg87; T reg137=reg3*reg84; reg70=reg118+reg70; reg118=reg104*reg52; T reg138=reg108*reg52;
    reg123=reg121+reg123; reg88=reg56+reg88; reg2=reg101*reg2; reg85=reg113+reg85; reg112=reg30-reg112;
    reg30=reg66*reg76; reg111=reg17+reg111; reg17=reg107*reg44; reg56=reg34+reg59; reg97=reg21*reg97;
    reg110=reg14+reg110; reg91=reg23*reg91; reg14=reg106*reg44; reg5=reg11*reg5; reg35=reg37-reg35;
    reg82=reg22*reg82; reg13=reg102+reg13; reg93=reg21*reg93; reg60=reg60/reg24; reg22=reg94*reg52;
    reg37=reg3*reg80; reg105=reg67+reg105; reg120=reg119+reg120; reg67=reg122*elem.pos(1)[0]; reg102=reg122*elem.pos(0)[0];
    reg113=reg108*reg16; reg117=reg116+reg117; reg116=reg7*reg0; reg15=reg89+reg15; reg115=reg114+reg115;
    reg89=reg94*reg16; reg136=reg135+reg136; reg133=reg74+reg133; reg37=reg105+reg37; reg137=reg65+reg137;
    reg65=reg2*reg98; reg74=reg2*reg62; reg105=reg108*reg60; reg14=reg110+reg14; reg129=reg120+reg129;
    reg58=reg78+reg58; reg134=reg88+reg134; reg89=reg115+reg89; reg113=reg117+reg113; reg125=reg103+reg125;
    reg128=reg39+reg128; reg30=reg123+reg30; reg96=reg127+reg96; reg39=reg94*reg60; reg78=reg104*reg60;
    reg17=reg111+reg17; reg88=reg3*reg75; reg103=elem.pos(2)[1]*var_inter[0]; reg110=reg95+reg132; reg15=reg21*reg15;
    reg5=reg35-reg5; reg21=elem.pos(2)[0]*var_inter[0]; reg35=reg126+reg124; reg97=reg91-reg97; reg56=reg11*reg56;
    reg91=reg93+reg82; reg138=reg70+reg138; reg118=reg85+reg118; reg70=var_inter[1]*elem.pos(2)[1]; reg130=reg130-reg131;
    reg112=reg23*reg112; reg4=reg116*reg4; reg23=var_inter[1]*elem.pos(2)[0]; reg67=reg67-reg102; reg13=reg22+reg13;
    reg22=reg128*reg17; reg85=reg58*reg14; reg15=reg112-reg15; reg105=reg133+reg105; reg88=reg30+reg88;
    reg30=reg109*elem.pos(3)[1]; reg103=reg103-reg110; reg111=reg134*reg14; reg112=reg109*elem.pos(3)[0]; reg39=reg96+reg39;
    reg21=reg21-reg35; reg91=reg11*reg91; reg11=var_inter[1]*elem.pos(3)[1]; reg70=reg130+reg70; reg96=reg4*reg73;
    reg74=reg37+reg74; reg37=reg72*reg83; reg114=var_inter[1]*elem.pos(3)[0]; reg23=reg67+reg23; reg67=reg99*reg118;
    reg115=reg129*reg27; reg117=reg6*reg72; reg119=reg7*reg83; reg120=reg113*reg27; reg121=reg99*reg138;
    reg123=reg89*reg27; reg127=reg99*reg13; reg130=reg4*reg20; reg65=reg137+reg65; reg5=reg5/reg24;
    reg56=reg97-reg56; reg97=reg125*reg17; reg78=reg136+reg78; reg133=reg2*reg77; reg135=reg134*reg128;
    reg136=reg125*reg58; reg97=reg111-reg97; reg123=reg127+reg123; reg111=reg39*reg37; reg127=reg28*reg118;
    reg117=reg119+reg117; reg119=reg129*reg18; reg137=reg7*reg72; T reg139=reg6*reg83; T reg140=reg113*reg18;
    T reg141=reg28*reg138; T reg142=reg89*reg18; T reg143=reg28*reg13; reg129=reg129*reg19; reg118=reg29*reg118;
    reg120=reg121+reg120; reg121=reg105*reg37; reg56=reg56/reg24; reg22=reg85-reg22; reg85=reg4*reg44;
    reg133=reg88+reg133; reg13=reg29*reg13; reg89=reg89*reg19; reg88=reg75*reg5; T reg144=reg51*reg5;
    T reg145=reg50*reg5; reg130=reg65+reg130; reg65=reg8*reg72; T reg146=reg9*reg83; reg113=reg113*reg19;
    reg91=reg15-reg91; reg15=reg78*reg37; reg23=reg23-reg114; reg115=reg67+reg115; reg96=reg74+reg96;
    reg70=reg70-reg11; reg112=reg21+reg112; reg30=reg103+reg30; reg138=reg29*reg138; reg21=reg78*reg94;
    reg129=reg118+reg129; reg142=reg143+reg142; reg67=reg117*reg56; reg74=reg39*reg108; reg103=reg137*reg56;
    reg118=reg139*reg56; reg89=reg13+reg89; reg39=reg39*reg94; reg13=reg70*reg112; reg143=reg8*reg83;
    T reg147=reg9*reg72; T reg148=reg145*reg75; reg111=reg123+reg111; reg136=reg135-reg136; reg123=reg105*reg94;
    reg140=reg141+reg140; reg105=reg105*reg108; reg113=reg138+reg113; reg119=reg127+reg119; reg78=reg78*reg108;
    reg127=reg23*reg30; reg135=reg130*reg22; reg24=reg91/reg24; reg15=reg115+reg15; reg91=reg88*reg75;
    reg65=reg146+reg65; reg85=reg133+reg85; reg115=reg97*reg96; reg121=reg120+reg121; reg120=reg144*reg75;
    reg10=reg10*reg72; reg39=reg89+reg39; reg89=reg145*reg26; reg133=reg118*reg117; reg148=reg111+reg148;
    reg111=reg88*reg26; reg138=reg134*reg85; reg141=reg144*reg25; reg105=reg140+reg105; reg144=reg144*reg26;
    reg21=reg129+reg21; reg129=reg136*reg85; reg12=reg83*reg12; reg115=reg135-reg115; reg123=reg113+reg123;
    reg113=reg130*reg17; reg135=reg103*reg117; reg140=reg96*reg17; reg91=reg15+reg91; reg15=reg58*reg85;
    reg145=reg145*reg25; reg146=reg67*reg117; T reg149=reg130*reg58; reg13=reg127-reg13; reg74=reg142+reg74;
    reg127=reg134*reg96; reg142=reg143*reg24; reg88=reg88*reg25; reg120=reg121+reg120; reg121=reg147*reg24;
    T reg150=reg65*reg24; reg78=reg119+reg78; reg119=reg67*reg10; reg70=reg70/reg13; reg112=reg112/reg13;
    reg15=reg140-reg15; reg30=reg30/reg13; reg88=reg78+reg88; reg78=reg125*reg85; reg23=reg23/reg13;
    reg141=reg105+reg141; reg127=reg149-reg127; reg105=reg103*reg10; reg138=reg113-reg138; reg113=reg130*reg128;
    reg140=reg125*reg96; reg145=reg74+reg145; reg74=reg118*reg10; reg149=reg121*reg65; T reg151=reg150*reg65;
    reg103=reg103*reg12; reg67=reg67*reg12; reg146=reg91+reg146; reg92=reg92*reg72; reg111=reg21+reg111;
    reg21=reg142*reg65; reg133=reg148+reg133; reg79=reg79*reg83; reg89=reg39+reg89; reg118=reg118*reg12;
    reg144=reg123+reg144; reg39=reg96*reg14; reg135=reg120+reg135; reg91=reg128*reg85; reg129=reg115+reg129;
    reg115=reg130*reg14; reg149=reg135+reg149; reg67=reg111+reg67; reg111=reg150*reg79; reg150=reg150*reg92;
    reg119=reg88+reg119; reg88=reg121*reg79; reg74=reg145+reg74; reg120=reg142*reg92; reg103=reg144+reg103;
    reg151=reg146+reg151; reg121=reg121*reg92; reg105=reg141+reg105; reg21=reg133+reg21; reg142=reg142*reg79;
    reg123=reg109*reg23; reg133=var_inter[1]*reg112; reg135=var_inter[1]*reg30; reg118=reg89+reg118; reg91=reg39-reg91;
    reg39=reg109*reg70; reg89=reg122*reg30; reg141=reg122*reg112; reg15=reg15/reg129; reg144=1-(*f.m).resolution;
    reg145=var_inter[0]*reg23; reg78=reg115-reg78; reg127=reg127/reg129; reg115=var_inter[0]*reg70; reg140=reg113-reg140;
    reg138=reg138/reg129; reg97=reg97/reg129; reg78=reg78/reg129; reg113=reg141-reg123; reg146=reg135-reg115;
    reg148=reg145-reg133; reg140=reg140/reg129; reg22=reg22/reg129; reg136=reg136/reg129; T reg152=reg89+reg115;
    T reg153=(*f.m).resolution*reg15; T reg154=(*f.m).resolution*reg138; T reg155=(*f.m).resolution*reg127; T reg156=reg141+reg145; T reg157=reg39+reg135;
    T reg158=reg123+reg133; reg150=reg119+reg150; reg121=reg105+reg121; reg120=reg74+reg120; reg111=reg67+reg111;
    reg88=reg103+reg88; reg67=reg39-reg89; reg129=reg91/reg129; reg74=reg144*reg151; reg91=reg144*reg149;
    reg103=reg144*reg21; reg142=reg118+reg142; reg105=0.5*reg146; reg118=0.5*reg113; reg119=0.5*reg67;
    T reg159=0.5*reg148; T reg160=0.5*reg158; T reg161=0.5*reg157; T reg162=0.5*reg156; T reg163=0.5*reg152;
    T reg164=reg144*reg150; T reg165=reg144*reg121; T reg166=reg144*reg120; T reg167=(*f.m).resolution*reg129; T reg168=reg144*reg111;
    T reg169=reg144*reg88; T reg170=(*f.m).resolution*reg78; reg155=reg74+reg155; reg154=reg91-reg154; reg153=reg103+reg153;
    reg74=(*f.m).resolution*reg140; reg91=reg142*reg144; reg103=(*f.m).resolution*reg136; T reg171=(*f.m).resolution*reg22; T reg172=(*f.m).resolution*reg97;
    reg167=reg166-reg167; reg170=reg165+reg170; reg74=reg164-reg74; reg164=reg157*reg153; reg165=reg160*reg155;
    reg168=reg103+reg168; reg172=reg169-reg172; reg91=reg171+reg91; reg103=reg152*reg153; reg166=reg162*reg155;
    reg169=reg118*reg155; reg171=reg159*reg155; T reg173=reg161*reg155; T reg174=reg158*reg154; T reg175=reg67*reg153;
    T reg176=reg148*reg154; T reg177=reg105*reg155; T reg178=reg146*reg153; T reg179=reg163*reg155; T reg180=reg113*reg154;
    T reg181=reg119*reg155; T reg182=reg156*reg154; T reg183=reg159*reg74; T reg184=reg146*reg167; T reg185=reg163*reg74;
    T reg186=reg157*reg167; T reg187=reg160*reg74; T reg188=reg148*reg170; T reg189=reg105*reg74; T reg190=reg113*reg170;
    T reg191=reg119*reg74; T reg192=reg152*reg167; T reg193=reg162*reg74; T reg194=reg156*reg170; T reg195=reg158*reg170;
    T reg196=reg161*reg74; reg177=reg176+reg177; reg176=reg105*reg168; T reg197=reg148*reg172; reg171=reg178+reg171;
    reg178=reg159*reg168; T reg198=reg146*reg91; reg179=reg179-reg182; T reg199=reg163*reg168; T reg200=reg156*reg172;
    reg103=reg103-reg166; T reg201=reg162*reg168; T reg202=reg152*reg91; reg181=reg180+reg181; reg174=reg174-reg173;
    reg180=reg67*reg91; T reg203=reg157*reg91; reg169=reg175+reg169; reg175=reg160*reg168; T reg204=reg119*reg168;
    T reg205=reg67*reg167; T reg206=reg118*reg74; T reg207=reg161*reg168; T reg208=reg158*reg172; T reg209=reg113*reg172;
    T reg210=reg118*reg168; reg165=reg165-reg164; reg204=reg209+reg204; reg174=2*reg174; reg183=reg184+reg183;
    reg103=2*reg103; reg181=2*reg181; reg195=reg195-reg196; reg208=reg208-reg207; reg189=reg188+reg189;
    reg202=reg202-reg201; reg175=reg175-reg203; reg187=reg187-reg186; reg165=2*reg165; reg169=2*reg169;
    reg179=2*reg179; reg171=2*reg171; reg191=reg190+reg191; reg178=reg198+reg178; reg210=reg180+reg210;
    reg192=reg192-reg193; reg199=reg199-reg200; reg176=reg197+reg176; reg206=reg205+reg206; reg185=reg185-reg194;
    reg177=2*reg177; reg180=reg119*reg174; reg184=reg146*reg204; reg188=reg148*reg187; reg190=reg161*reg165;
    reg197=reg179*reg118; reg198=reg152*reg210; reg205=reg118*reg181; reg209=reg195*reg113; T reg211=reg179*reg163;
    T reg212=reg195*reg158; T reg213=reg202*reg67; T reg214=reg199*reg67; T reg215=reg157*reg210; T reg216=reg177*reg105;
    T reg217=reg161*reg174; T reg218=reg148*reg189; T reg219=reg156*reg192; T reg220=reg105*reg165; T reg221=reg103*reg118;
    T reg222=reg176*reg152; T reg223=reg177*reg162; T reg224=reg159*reg174; T reg225=reg146*reg208; T reg226=reg158*reg192;
    T reg227=reg103*reg159; T reg228=reg152*reg175; T reg229=reg162*reg165; T reg230=reg159*reg165; T reg231=reg146*reg175;
    T reg232=reg176*reg146; T reg233=reg156*reg191; T reg234=reg152*reg208; T reg235=reg162*reg174; T reg236=reg171*reg159;
    T reg237=reg178*reg146; T reg238=reg163*reg181; T reg239=reg163*reg169; T reg240=reg156*reg206; T reg241=reg179*reg159;
    T reg242=reg199*reg146; T reg243=reg162*reg169; T reg244=reg171*reg105; T reg245=reg148*reg183; T reg246=reg159*reg181;
    T reg247=reg152*reg204; T reg248=reg162*reg181; T reg249=reg179*reg105; T reg250=reg148*reg185; T reg251=reg202*reg152;
    T reg252=reg103*reg162; T reg253=reg103*reg105; T reg254=reg148*reg192; T reg255=reg103*reg161; T reg256=reg199*reg152;
    T reg257=reg179*reg162; T reg258=reg105*reg181; T reg259=reg148*reg191; T reg260=reg178*reg152; T reg261=reg171*reg162;
    T reg262=reg105*reg169; T reg263=reg148*reg206; T reg264=reg202*reg146; T reg265=reg171*reg160; T reg266=reg118*reg174;
    T reg267=reg113*reg185; T reg268=reg179*reg119; T reg269=reg178*reg157; T reg270=reg177*reg163; T reg271=reg177*reg159;
    T reg272=reg113*reg183; T reg273=reg171*reg119; T reg274=reg177*reg160; T reg275=reg176*reg157; T reg276=reg67*reg208;
    T reg277=reg113*reg189; T reg278=reg177*reg119; T reg279=reg160*reg165; T reg280=reg157*reg175; T reg281=reg156*reg183;
    T reg282=reg113*reg206; T reg283=reg119*reg169; T reg284=reg195*reg148; T reg285=reg160*reg181; T reg286=reg113*reg191;
    T reg287=reg119*reg181; T reg288=reg157*reg204; T reg289=reg156*reg189; T reg290=reg103*reg160; T reg291=reg103*reg163;
    reg202=reg202*reg157; T reg292=reg163*reg165; reg192=reg113*reg192; reg103=reg103*reg119; T reg293=reg179*reg160;
    reg199=reg199*reg157; T reg294=reg156*reg187; T reg295=reg158*reg185; reg179=reg179*reg161; reg204=reg67*reg204;
    T reg296=reg160*reg169; T reg297=reg177*reg118; reg185=reg156*reg185; reg183=reg158*reg183; reg176=reg176*reg67;
    T reg298=reg171*reg161; T reg299=reg118*reg169; T reg300=reg146*reg210; reg189=reg158*reg189; reg177=reg177*reg161;
    T reg301=reg171*reg118; reg210=reg67*reg210; reg178=reg178*reg67; T reg302=reg159*reg169; T reg303=reg158*reg187;
    reg187=reg113*reg187; T reg304=reg119*reg165; T reg305=reg163*reg174; T reg306=reg160*reg174; reg208=reg157*reg208;
    reg195=reg195*reg156; reg206=reg158*reg206; reg169=reg161*reg169; reg165=reg118*reg165; reg175=reg67*reg175;
    reg174=reg105*reg174; reg171=reg171*reg163; reg191=reg158*reg191; reg181=reg161*reg181; reg219=reg291-reg219;
    reg233=reg238-reg233; reg281=reg171-reg281; reg289=reg270-reg289; reg255=reg226-reg255; reg185=reg211-reg185;
    reg215=reg296-reg215; reg288=reg285-reg288; reg202=reg290-reg202; reg199=reg293-reg199; reg269=reg265-reg269;
    reg275=reg274-reg275; reg280=reg279-reg280; reg271=reg232+reg271; reg208=reg306-reg208; reg169=reg206-reg169;
    reg181=reg191-reg181; reg179=reg295-reg179; reg298=reg183-reg298; reg177=reg189-reg177; reg190=reg303-reg190;
    reg217=reg212-reg217; reg220=reg188+reg220; reg180=reg209+reg180; reg243=reg198-reg243; reg248=reg247-reg248;
    reg252=reg251-reg252; reg257=reg256-reg257; reg261=reg260-reg261; reg223=reg222-reg223; reg229=reg228-reg229;
    reg235=reg234-reg235; reg240=reg239-reg240; reg278=reg277+reg278; reg227=reg264+reg227; reg304=reg187+reg304;
    reg241=reg242+reg241; reg165=reg175+reg165; reg174=reg284+reg174; reg236=reg237+reg236; reg230=reg231+reg230;
    reg297=reg176+reg297; reg224=reg225+reg224; reg299=reg210+reg299; reg301=reg178+reg301; reg262=reg263+reg262;
    reg197=reg214+reg197; reg258=reg259+reg258; reg253=reg254+reg253; reg221=reg213+reg221; reg216=reg218+reg216;
    reg249=reg250+reg249; reg204=reg205+reg204; reg244=reg245+reg244; reg283=reg282+reg283; reg294=reg292-reg294;
    reg287=reg286+reg287; reg195=reg305-reg195; reg103=reg192+reg103; reg302=reg300+reg302; reg266=reg276+reg266;
    reg268=reg267+reg268; reg273=reg272+reg273; reg246=reg184+reg246; reg177=reg177*reg13; reg197=reg197*reg13;
    reg278=reg13*reg278; reg288=reg288*reg13; reg190=reg190*reg13; reg244=reg244*reg13; reg217=reg217*reg13;
    reg220=reg220*reg13; reg287=reg13*reg287; reg221=reg221*reg13; reg283=reg13*reg283; reg216=reg216*reg13;
    reg273=reg13*reg273; reg204=reg204*reg13; reg180=reg180*reg13; reg215=reg215*reg13; reg275=reg275*reg13;
    reg208=reg208*reg13; reg165=reg13*reg165; reg268=reg13*reg268; reg269=reg269*reg13; reg304=reg13*reg304;
    reg169=reg169*reg13; reg271=reg271*reg13; reg174=reg174*reg13; reg181=reg181*reg13; reg199=reg199*reg13;
    reg297=reg13*reg297; reg179=reg179*reg13; reg299=reg13*reg299; reg103=reg13*reg103; reg280=reg280*reg13;
    reg298=reg298*reg13; reg301=reg301*reg13; reg202=reg202*reg13; reg266=reg13*reg266; reg262=reg262*reg13;
    reg261=reg261*reg13; reg224=reg224*reg13; reg223=reg223*reg13; reg230=reg230*reg13; reg229=reg229*reg13;
    reg236=reg236*reg13; reg235=reg235*reg13; reg241=reg241*reg13; reg240=reg240*reg13; reg227=reg227*reg13;
    reg233=reg233*reg13; reg246=reg246*reg13; reg255=reg255*reg13; reg302=reg302*reg13; reg219=reg219*reg13;
    reg195=reg195*reg13; reg185=reg185*reg13; reg294=reg294*reg13; reg281=reg281*reg13; reg289=reg289*reg13;
    reg249=reg249*reg13; reg252=reg252*reg13; reg248=reg248*reg13; reg243=reg243*reg13; reg258=reg258*reg13;
    reg253=reg253*reg13; reg257=reg257*reg13; T tmp_2_1=ponderation*reg248; T tmp_6_5=ponderation*reg275; T tmp_4_1=ponderation*reg246;
    T tmp_1_4=ponderation*reg273; T tmp_3_1=ponderation*reg233; T tmp_6_4=ponderation*reg269; T tmp_1_3=ponderation*reg268; T tmp_5_2=ponderation*reg253;
    T tmp_1_5=ponderation*reg278; T tmp_4_0=ponderation*reg302; T tmp_5_5=ponderation*reg216; T tmp_0_1=ponderation*reg204; T tmp_6_3=ponderation*reg199;
    T tmp_5_3=ponderation*reg249; T tmp_3_5=ponderation*reg289; T tmp_3_4=ponderation*reg281; T tmp_6_0=ponderation*reg215; T tmp_1_0=ponderation*reg283;
    T tmp_1_7=ponderation*reg180; T tmp_3_6=ponderation*reg294; T tmp_6_1=ponderation*reg288; T tmp_1_1=ponderation*reg287; T tmp_0_7=ponderation*reg266;
    T tmp_3_3=ponderation*reg185; T tmp_2_0=ponderation*reg243; T tmp_5_6=ponderation*reg220; T tmp_6_2=ponderation*reg202; T tmp_3_7=ponderation*reg195;
    T tmp_1_2=ponderation*reg103; T tmp_3_2=ponderation*reg219; T tmp_7_1=ponderation*reg181; T tmp_7_2=ponderation*reg255; T tmp_5_4=ponderation*reg244;
    T tmp_2_5=ponderation*reg223; T tmp_0_0=ponderation*reg299; T tmp_2_2=ponderation*reg252; T tmp_7_3=ponderation*reg179; T tmp_4_6=ponderation*reg230;
    T tmp_7_6=ponderation*reg190; T tmp_0_4=ponderation*reg301; T tmp_2_4=ponderation*reg261; T tmp_7_4=ponderation*reg298; T tmp_4_7=ponderation*reg224;
    T tmp_0_3=ponderation*reg197; T tmp_5_0=ponderation*reg262; T tmp_7_5=ponderation*reg177; T tmp_2_3=ponderation*reg257; T tmp_3_0=ponderation*reg240;
    T tmp_5_1=ponderation*reg258; T tmp_6_6=ponderation*reg280; T tmp_4_2=ponderation*reg227; T tmp_1_6=ponderation*reg304; T tmp_0_6=ponderation*reg165;
    T tmp_4_5=ponderation*reg271; T tmp_6_7=ponderation*reg208; T tmp_2_7=ponderation*reg235; T tmp_7_7=ponderation*reg217; T tmp_4_3=ponderation*reg241;
    T tmp_2_6=ponderation*reg229; T tmp_7_0=ponderation*reg169; T tmp_0_5=ponderation*reg297; T tmp_0_2=ponderation*reg221; T tmp_5_7=ponderation*reg174;
    T tmp_4_4=ponderation*reg236;
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
    T reg3=pow((*f.m).v2[1],2); reg0=reg1+reg0; reg3=reg2+reg3; reg1=pow((*f.m).v2[2],2); reg2=2*(*f.m).shear_modulus_13;
    T reg4=2*(*f.m).shear_modulus_23; reg1=reg3+reg1; reg4=1.0/reg4; reg0=pow(reg0,0.5); reg3=2*(*f.m).shear_modulus_12;
    reg2=1.0/reg2; T reg5=reg2*reg4; T reg6=(*f.m).v1[0]/reg0; T reg7=(*f.m).v1[1]/reg0; reg1=pow(reg1,0.5);
    reg3=1.0/reg3; reg0=(*f.m).v1[2]/reg0; T reg8=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=(*f.m).v2[0]/reg1;
    T reg11=(*f.m).v2[1]/reg1; T reg12=2*reg7; T reg13=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg14=reg3*reg5; T reg15=2*reg6;
    T reg16=reg13*reg14; T reg17=reg8*reg14; T reg18=reg10*reg15; T reg19=reg12*reg11; T reg20=2*reg0;
    T reg21=1.0/(*f.m).elastic_modulus_1; T reg22=1.0/(*f.m).elastic_modulus_2; T reg23=reg9*reg14; T reg24=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg25=pow(reg10,2);
    reg1=(*f.m).v2[2]/reg1; T reg26=pow(reg11,2); T reg27=reg19*reg22; T reg28=reg18*reg24; T reg29=reg26*reg22;
    T reg30=reg25*reg24; T reg31=reg19*reg24; T reg32=pow(reg1,2); T reg33=reg18*reg21; T reg34=reg26*reg24;
    T reg35=reg25*reg21; T reg36=reg8*reg16; T reg37=reg22*reg23; reg20=reg1*reg20; T reg38=reg24*reg23;
    T reg39=reg8*reg17; T reg40=reg19*reg8; T reg41=reg13*reg18; T reg42=reg26*reg8; T reg43=reg25*reg13;
    reg28=reg27-reg28; reg27=reg20*reg8; reg30=reg29-reg30; reg29=reg32*reg8; reg31=reg33-reg31;
    reg33=reg13*reg20; reg39=reg37-reg39; reg34=reg35-reg34; reg35=reg32*reg13; reg36=reg38+reg36;
    reg37=reg24*reg17; T reg44=reg22*reg16; T reg45=pow(reg7,2); T reg46=pow(reg6,2); T reg47=reg24*reg36;
    T reg48=reg45*reg22; T reg49=reg45*reg24; T reg50=reg21*reg39; T reg51=reg37+reg44; T reg52=pow(reg0,2);
    reg29=reg30-reg29; reg30=reg7*reg11; T reg53=reg46*reg24; T reg54=reg6*reg10; reg33=reg31-reg33;
    reg27=reg28-reg27; reg35=reg34-reg35; reg28=reg6*reg11; reg31=reg7*reg10; reg34=reg32*reg9;
    reg42=reg43+reg42; reg20=reg20*reg9; reg43=reg46*reg21; reg40=reg41+reg40; reg41=reg13*reg5;
    reg47=reg50-reg47; reg50=reg8*reg5; T reg55=reg13*reg51; T reg56=reg35*reg54; T reg57=reg29*reg30;
    T reg58=reg52*reg13; T reg59=reg3*reg4; T reg60=reg22*reg14; reg14=reg24*reg14; T reg61=reg13*reg17;
    reg49=reg43-reg49; reg43=reg13*reg16; T reg62=reg33*reg46; T reg63=reg29*reg45; T reg64=reg35*reg46;
    T reg65=reg27*reg45; reg23=reg21*reg23; reg40=reg20-reg40; reg42=reg34-reg42; reg20=reg45*reg8;
    reg34=reg46*reg13; T reg66=reg6*reg1; T reg67=reg0*reg10; T reg68=reg0*reg11; T reg69=reg7*reg1;
    T reg70=reg10*reg11; T reg71=reg27*reg26; T reg72=reg33*reg25; T reg73=reg27*reg30; T reg74=reg33*reg54;
    T reg75=reg29*reg26; T reg76=reg31+reg28; T reg77=reg35*reg25; T reg78=reg52*reg8; reg5=reg9*reg5;
    reg53=reg48-reg53; reg48=reg0*reg1; T reg79=2*reg10; T reg80=reg24*reg5; T reg81=reg8*reg50;
    reg20=reg34+reg20; reg34=reg52*reg9; T reg82=reg42*reg52; T reg83=reg8*reg41; reg75=reg77+reg75;
    reg77=reg42*reg32; T reg84=reg40*reg52; reg73=reg74+reg73; reg74=reg40*reg48; reg78=reg53-reg78;
    reg63=reg64+reg63; reg53=reg69-reg68; reg64=reg66+reg67; T reg85=reg9*reg59; reg71=reg72+reg71;
    reg72=reg40*reg32; T reg86=reg10*reg1; reg65=reg62+reg65; reg55=reg47-reg55; reg47=reg79*reg11;
    reg62=reg76*reg3; reg16=reg24*reg16; T reg87=reg8*reg59; reg5=reg22*reg5; T reg88=reg70*reg3;
    reg57=reg56+reg57; reg56=reg13*reg14; reg58=reg49-reg58; reg49=reg42*reg48; reg17=reg21*reg17;
    reg43=reg23-reg43; reg59=reg13*reg59; reg23=reg3*reg2; T reg89=reg7*reg15; reg61=reg38+reg61;
    reg38=reg13*reg60; T reg90=reg64*reg2; T reg91=reg78*reg45; T reg92=reg58*reg46; T reg93=reg8*reg87;
    T reg94=reg11*reg1; T reg95=reg86*reg2; T reg96=reg6*reg7; T reg97=reg8*reg59; reg66=reg67-reg66;
    reg67=reg62*reg89; reg84=reg65+reg84; reg20=reg34-reg20; reg9=reg9*reg23; reg34=2*reg53;
    reg65=reg0*reg15; reg50=reg24*reg50; T reg98=reg8*reg23; reg41=reg22*reg41; reg14=reg24*reg14;
    T reg99=reg58*reg25; reg16=reg17+reg16; reg56=reg17+reg56; reg17=reg78*reg26; reg23=reg13*reg23;
    reg43=reg43/reg55; reg36=reg36/reg55; reg60=reg21*reg60; reg38=reg37+reg38; reg61=reg61/reg55;
    T reg100=reg88*reg89; reg49=reg57+reg49; reg57=reg24*reg85; reg85=reg22*reg85; T reg101=reg62*reg76;
    reg74=reg73+reg74; reg73=reg79*reg1; T reg102=reg62*reg47; reg72=reg71+reg72; reg71=2*reg11;
    T reg103=reg88*reg47; reg77=reg75+reg77; reg83=reg80+reg83; reg39=reg39/reg55; reg69=reg68+reg69;
    reg82=reg63+reg82; reg63=reg88*reg76; reg81=reg5-reg81; reg67=reg84+reg67; reg5=reg90*reg65;
    reg68=reg71*reg1; reg17=reg99+reg17; reg75=reg20*reg32; reg63=reg49+reg63; reg49=reg95*reg64;
    reg103=reg77+reg103; reg77=reg95*reg73; reg102=reg72+reg102; reg72=reg90*reg73; reg101=reg74+reg101;
    reg74=reg6*reg0; reg80=reg34*reg66; reg84=reg12*reg0; reg99=reg90*reg64; T reg104=pow(reg66,2);
    T reg105=pow(reg53,2); T reg106=1-var_inter[0]; reg51=reg51/reg55; T reg107=reg45*reg61; T reg108=reg26*reg43;
    T reg109=reg89*reg61; T reg110=reg47*reg43; reg38=reg38/reg55; T reg111=reg46*reg39; reg56=reg56/reg55;
    reg16=reg16/reg55; reg14=reg60-reg14; reg87=reg24*reg87; reg59=reg22*reg59; reg60=reg25*reg36;
    reg81=reg21*reg81; reg83=reg24*reg83; T reg112=reg50+reg41; reg93=reg85-reg93; reg97=reg57+reg97;
    reg57=reg22*reg9; reg9=reg24*reg9; reg85=reg8*reg98; T reg113=reg47*reg36; reg8=reg8*reg23;
    T reg114=reg89*reg39; T reg115=reg45*reg39; T reg116=reg26*reg36; T reg117=reg95*reg65; reg100=reg82+reg100;
    reg82=1-var_inter[1]; T reg118=reg20*reg52; reg91=reg92+reg91; reg92=reg69*reg4; T reg119=reg94*reg4;
    T reg120=reg78*reg30; T reg121=reg58*reg54; T reg122=reg46*reg61; T reg123=reg25*reg43; reg3=reg96*reg3;
    reg14=reg14/reg55; T reg124=reg20*reg48; reg117=reg100+reg117; reg98=reg24*reg98; reg60=reg111+reg60;
    reg23=reg22*reg23; reg22=reg92*reg84; reg100=reg89*reg38; reg83=reg81-reg83; reg5=reg67+reg5;
    reg112=reg13*reg112; reg67=reg105*reg16; reg81=reg47*reg56; reg93=reg21*reg93; reg118=reg91+reg118;
    reg97=reg24*reg97; reg91=reg3*reg89; reg75=reg17+reg75; reg17=reg87+reg59; reg111=reg3*reg47;
    reg2=reg74*reg2; reg99=reg101+reg99; reg101=reg80*reg51; T reg125=reg92*reg69; reg113=reg114+reg113;
    reg85=reg57-reg85; reg123=reg122+reg123; reg57=reg119*reg84; reg8=reg9+reg8; reg9=reg104*reg51;
    reg116=reg115+reg116; reg114=elem.pos(1)[1]*var_inter[0]; reg115=reg104*reg16; reg122=elem.pos(0)[1]*reg106; reg108=reg107+reg108;
    reg107=reg7*reg0; reg110=reg109+reg110; reg109=elem.pos(1)[0]*var_inter[0]; T reg126=reg80*reg16; T reg127=reg46*reg38;
    T reg128=reg92*reg68; reg72=reg102+reg72; reg102=reg106*elem.pos(0)[0]; T reg129=reg25*reg56; T reg130=reg119*reg69;
    reg49=reg63+reg49; reg63=reg82*elem.pos(1)[1]; reg120=reg121+reg120; reg121=reg105*reg51; T reg131=reg82*elem.pos(0)[0];
    T reg132=reg82*elem.pos(1)[0]; T reg133=reg26*reg56; T reg134=reg45*reg38; T reg135=reg82*elem.pos(0)[1]; T reg136=reg119*reg68;
    reg77=reg103+reg77; reg103=reg104*reg14; reg129=reg127+reg129; reg115=reg108+reg115; reg9=reg116+reg9;
    reg81=reg100+reg81; reg67=reg123+reg67; reg125=reg99+reg125; reg99=reg80*reg14; reg124=reg120+reg124;
    reg100=reg105*reg14; reg126=reg110+reg126; reg108=reg3*reg76; reg101=reg113+reg101; reg133=reg134+reg133;
    reg111=reg75+reg111; reg130=reg49+reg130; reg60=reg121+reg60; reg132=reg132-reg131; reg49=var_inter[1]*elem.pos(2)[0];
    reg136=reg77+reg136; reg63=reg63-reg135; reg75=var_inter[1]*elem.pos(2)[1]; reg128=reg72+reg128; reg72=reg102+reg109;
    reg77=elem.pos(2)[0]*var_inter[0]; reg110=reg122+reg114; reg113=elem.pos(2)[1]*var_inter[0]; reg4=reg107*reg4; reg22=reg5+reg22;
    reg91=reg118+reg91; reg5=reg2*reg65; reg57=reg117+reg57; reg17=reg13*reg17; reg8=reg24*reg8;
    reg97=reg93-reg97; reg24=reg98+reg23; reg112=reg83-reg112; reg85=reg21*reg85; reg21=reg2*reg73;
    reg75=reg63+reg75; reg63=var_inter[1]*elem.pos(3)[1]; reg17=reg97-reg17; reg83=reg22*reg130; reg8=reg85-reg8;
    reg99=reg81+reg99; reg77=reg77-reg72; reg81=reg106*elem.pos(3)[0]; reg85=reg66*reg53; reg93=reg4*reg68;
    reg21=reg111+reg21; reg97=reg57*reg125; reg113=reg113-reg110; reg111=reg106*elem.pos(3)[1]; reg116=reg128*reg130;
    reg117=reg136*reg125; reg112=reg112/reg55; reg118=reg4*reg84; reg120=reg126*reg70; reg121=reg96*reg101;
    reg5=reg91+reg5; reg91=reg115*reg70; reg123=reg96*reg9; reg127=reg2*reg64; reg134=reg67*reg70;
    T reg137=reg96*reg60; reg108=reg124+reg108; reg103=reg133+reg103; reg24=reg13*reg24; reg49=reg132+reg49;
    reg100=reg129+reg100; reg13=var_inter[1]*elem.pos(3)[0]; reg124=reg6*reg66; reg129=reg7*reg53; reg127=reg108+reg127;
    reg116=reg117-reg116; reg24=reg8-reg24; reg83=reg97-reg83; reg8=reg57*reg128; reg97=reg22*reg136;
    reg93=reg21+reg93; reg21=reg10*reg66; reg108=reg11*reg53; reg17=reg17/reg55; reg117=reg4*reg69;
    reg134=reg137+reg134; reg132=reg100*reg85; reg91=reg123+reg91; reg123=reg103*reg85; reg120=reg121+reg120;
    reg121=reg99*reg85; reg133=reg6*reg53; reg137=reg7*reg66; T reg138=reg67*reg26; T reg139=reg46*reg9;
    reg81=reg77+reg81; reg111=reg113+reg111; reg77=reg45*reg60; reg113=reg45*reg101; T reg140=reg76*reg112;
    T reg141=reg30*reg112; T reg142=reg126*reg26; T reg143=reg115*reg25; reg49=reg49-reg13; reg60=reg46*reg60;
    reg75=reg75-reg63; reg115=reg115*reg26; T reg144=reg54*reg112; reg101=reg46*reg101; reg126=reg126*reg25;
    reg118=reg5+reg118; reg67=reg67*reg25; reg9=reg45*reg9; reg124=reg129+reg124; reg97=reg8-reg97;
    reg126=reg101+reg126; reg5=reg99*reg105; reg21=reg108+reg21; reg55=reg24/reg55; reg8=reg49*reg111;
    reg143=reg139+reg143; reg24=reg103*reg105; reg142=reg113+reg142; reg99=reg99*reg104; reg101=reg118*reg116;
    reg108=reg83*reg93; reg103=reg103*reg104; reg113=reg140*reg76; reg121=reg120+reg121; reg120=reg133*reg17;
    reg115=reg9+reg115; reg9=reg144*reg76; reg123=reg91+reg123; reg91=reg75*reg81; reg132=reg134+reg132;
    reg129=reg141*reg76; reg134=reg100*reg104; reg138=reg77+reg138; reg77=reg137*reg17; reg117=reg127+reg117;
    reg127=reg124*reg17; reg139=reg11*reg66; T reg145=reg10*reg53; reg67=reg60+reg67; reg100=reg100*reg105;
    reg60=reg118*reg130; T reg146=reg140*reg19; reg99=reg142+reg99; reg103=reg115+reg103; reg115=reg141*reg19;
    reg142=reg136*reg117; T reg147=reg144*reg19; T reg148=reg57*reg117; reg108=reg101-reg108; reg101=reg97*reg117;
    reg134=reg138+reg134; reg138=reg77*reg124; reg129=reg123+reg129; reg123=reg120*reg124; T reg149=reg127*reg124;
    reg113=reg121+reg113; reg9=reg132+reg9; reg91=reg8-reg91; reg8=reg145*reg55; reg121=reg139*reg55;
    reg132=reg21*reg55; reg100=reg67+reg100; reg144=reg144*reg18; reg15=reg53*reg15; reg67=reg93*reg130;
    T reg150=reg57*reg93; reg5=reg126+reg5; reg126=reg118*reg136; reg12=reg12*reg66; reg140=reg140*reg18;
    reg24=reg143+reg24; reg141=reg141*reg18; reg143=reg120*reg12; reg147=reg134+reg147; reg49=reg49/reg91;
    reg134=reg8*reg21; reg123=reg9+reg123; reg9=reg118*reg125; T reg151=reg127*reg15; reg140=reg5+reg140;
    reg138=reg129+reg138; reg5=reg121*reg21; reg71=reg71*reg66; reg142=reg67-reg142; reg79=reg79*reg53;
    reg67=reg128*reg117; reg127=reg127*reg12; reg120=reg120*reg15; reg144=reg100+reg144; reg148=reg60-reg148;
    reg146=reg99+reg146; reg60=reg77*reg15; reg101=reg108+reg101; reg141=reg24+reg141; reg75=reg75/reg91;
    reg24=reg118*reg128; reg99=reg93*reg125; reg81=reg81/reg91; reg111=reg111/reg91; reg149=reg113+reg149;
    reg100=reg132*reg21; reg108=reg22*reg117; reg113=reg22*reg93; reg77=reg77*reg12; reg115=reg103+reg115;
    reg150=reg126-reg150; reg103=var_inter[1]*reg111; reg126=reg82*reg111; reg129=reg82*reg81; reg120=reg144+reg120;
    reg144=reg8*reg79; T reg152=reg106*reg75; reg60=reg141+reg60; reg141=var_inter[1]*reg81; reg5=reg138+reg5;
    reg138=var_inter[0]*reg75; T reg153=reg121*reg79; reg100=reg149+reg100; reg134=reg123+reg134; reg123=var_inter[0]*reg49;
    reg127=reg146+reg127; reg146=reg132*reg71; reg108=reg9-reg108; reg148=reg148/reg101; reg142=reg142/reg101;
    reg121=reg121*reg71; reg77=reg115+reg77; reg132=reg132*reg79; reg151=reg140+reg151; reg143=reg147+reg143;
    reg67=reg99-reg67; reg9=1-(*f.m).resolution; reg150=reg150/reg101; reg113=reg24-reg113; reg8=reg8*reg71;
    reg24=reg106*reg49; reg99=(*f.m).resolution*reg150; reg115=(*f.m).resolution*reg148; reg140=reg129+reg123; reg8=reg143+reg8;
    reg143=(*f.m).resolution*reg142; reg108=reg108/reg101; reg83=reg83/reg101; reg153=reg60+reg153; reg67=reg67/reg101;
    reg132=reg151+reg132; reg116=reg116/reg101; reg146=reg127+reg146; reg144=reg120+reg144; reg60=reg129-reg24;
    reg120=reg103-reg138; reg113=reg113/reg101; reg127=reg152-reg126; reg147=reg123-reg141; reg149=reg9*reg134;
    reg151=reg9*reg5; T reg154=reg9*reg100; T reg155=reg24+reg141; T reg156=reg152+reg103; reg101=reg97/reg101;
    reg121=reg77+reg121; reg77=reg126+reg138; reg97=0.5*reg155; T reg157=0.5*reg127; T reg158=0.5*reg60;
    T reg159=0.5*reg156; T reg160=(*f.m).resolution*reg83; T reg161=(*f.m).resolution*reg116; T reg162=(*f.m).resolution*reg101; T reg163=reg144*reg9;
    T reg164=reg9*reg153; T reg165=reg9*reg132; T reg166=reg9*reg8; T reg167=reg9*reg121; T reg168=reg9*reg146;
    T reg169=(*f.m).resolution*reg67; T reg170=(*f.m).resolution*reg108; T reg171=0.5*reg120; reg99=reg154+reg99; reg115=reg151-reg115;
    reg143=reg149+reg143; reg149=(*f.m).resolution*reg113; reg151=0.5*reg77; reg154=0.5*reg140; T reg172=0.5*reg147;
    T reg173=reg154*reg99; T reg174=reg77*reg143; T reg175=reg140*reg115; T reg176=reg151*reg99; T reg177=reg127*reg143;
    T reg178=reg158*reg99; T reg179=reg159*reg99; T reg180=reg120*reg143; T reg181=reg172*reg99; T reg182=reg155*reg115;
    T reg183=reg60*reg115; reg149=reg168-reg149; reg170=reg167+reg170; reg160=reg164-reg160; reg169=reg166-reg169;
    reg164=reg97*reg99; reg166=reg171*reg99; reg167=reg147*reg115; reg163=reg161+reg163; reg161=reg156*reg143;
    reg168=reg157*reg99; reg165=reg162+reg165; reg162=reg154*reg165; T reg184=reg97*reg165; T reg185=reg77*reg163;
    T reg186=reg155*reg170; reg174=reg174-reg173; T reg187=reg159*reg149; T reg188=reg154*reg149; T reg189=reg77*reg169;
    T reg190=reg127*reg163; T reg191=reg158*reg165; T reg192=reg157*reg149; T reg193=reg60*reg170; reg164=reg164-reg161;
    reg178=reg177+reg178; reg182=reg182-reg179; reg177=reg156*reg163; T reg194=reg60*reg160; T reg195=reg157*reg165;
    T reg196=reg159*reg165; reg168=reg183+reg168; reg183=reg155*reg160; T reg197=reg147*reg170; T reg198=reg97*reg149;
    T reg199=reg156*reg169; reg181=reg180+reg181; reg176=reg176-reg175; reg180=reg120*reg169; T reg200=reg171*reg149;
    T reg201=reg172*reg165; reg166=reg167+reg166; reg167=reg171*reg165; T reg202=reg151*reg165; T reg203=reg140*reg160;
    T reg204=reg151*reg149; T reg205=reg120*reg163; T reg206=reg140*reg170; T reg207=reg147*reg160; T reg208=reg172*reg149;
    reg182=2*reg182; reg178=2*reg178; reg186=reg186-reg187; reg201=reg205+reg201; reg164=2*reg164;
    reg168=2*reg168; reg181=2*reg181; reg167=reg207+reg167; reg195=reg194+reg195; reg183=reg183-reg196;
    reg208=reg180+reg208; reg184=reg184-reg177; reg174=2*reg174; reg200=reg197+reg200; reg185=reg185-reg162;
    reg202=reg202-reg203; reg189=reg189-reg188; reg176=2*reg176; reg198=reg198-reg199; reg204=reg204-reg206;
    reg166=2*reg166; reg192=reg193+reg192; reg191=reg190+reg191; reg180=reg77*reg184; reg190=reg167*reg77;
    reg193=reg154*reg164; reg194=reg181*reg154; reg197=reg201*reg77; reg205=reg166*reg154; reg207=reg77*reg183;
    T reg209=reg159*reg182; T reg210=reg186*reg155; T reg211=reg156*reg183; T reg212=reg147*reg200; T reg213=reg97*reg182;
    T reg214=reg166*reg171; T reg215=reg147*reg198; T reg216=reg171*reg164; T reg217=reg127*reg191; T reg218=reg158*reg178;
    T reg219=reg186*reg60; T reg220=reg186*reg147; T reg221=reg157*reg182; T reg222=reg166*reg172; T reg223=reg185*reg77;
    T reg224=reg174*reg154; T reg225=reg127*reg195; T reg226=reg171*reg182; T reg227=reg202*reg77; T reg228=reg176*reg154;
    T reg229=reg127*reg183; T reg230=reg201*reg120; T reg231=reg181*reg172; T reg232=reg166*reg158; T reg233=reg167*reg127;
    T reg234=reg167*reg120; T reg235=reg120*reg184; T reg236=reg60*reg208; T reg237=reg172*reg164; T reg238=reg158*reg182;
    reg183=reg120*reg183; T reg239=reg172*reg182; T reg240=reg181*reg158; T reg241=reg60*reg192; T reg242=reg157*reg168;
    T reg243=reg158*reg168; T reg244=reg201*reg127; T reg245=reg60*reg189; T reg246=reg174*reg157; T reg247=reg185*reg127;
    T reg248=reg174*reg158; T reg249=reg202*reg127; T reg250=reg176*reg158; T reg251=reg176*reg157; T reg252=reg60*reg204;
    T reg253=reg154*reg182; T reg254=reg127*reg184; T reg255=reg176*reg151; T reg256=reg157*reg164; T reg257=reg140*reg204;
    T reg258=reg181*reg151; T reg259=reg60*reg198; T reg260=reg97*reg164; T reg261=reg140*reg208; reg184=reg156*reg184;
    T reg262=reg158*reg164; T reg263=reg166*reg151; T reg264=reg140*reg200; T reg265=reg181*reg157; reg186=reg186*reg140;
    T reg266=reg166*reg157; T reg267=reg151*reg164; T reg268=reg60*reg200; T reg269=reg140*reg198; reg182=reg151*reg182;
    reg184=reg260-reg184; reg256=reg259+reg256; reg216=reg215+reg216; reg251=reg252+reg251; reg265=reg236+reg265;
    reg211=reg213-reg211; reg222=reg234+reg222; reg266=reg268+reg266; reg226=reg220+reg226; reg240=reg244+reg240;
    reg209=reg210-reg209; reg253=reg207-reg253; reg193=reg180-reg193; reg257=reg255-reg257; reg232=reg233+reg232;
    reg205=reg190-reg205; reg261=reg258-reg261; reg194=reg197-reg194; reg264=reg263-reg264; reg262=reg254+reg262;
    reg228=reg227-reg228; reg269=reg267-reg269; reg224=reg223-reg224; reg186=reg182-reg186; reg221=reg219+reg221;
    reg231=reg230+reg231; reg218=reg217+reg218; reg237=reg235+reg237; reg239=reg183+reg239; reg238=reg229+reg238;
    reg242=reg241+reg242; reg225=reg243+reg225; reg246=reg245+reg246; reg248=reg247+reg248; reg214=reg212+reg214;
    reg250=reg249+reg250; reg218=reg91*reg218; reg214=reg214*reg91; reg228=reg228*reg91; reg216=reg216*reg91;
    reg194=reg194*reg91; reg222=reg222*reg91; reg226=reg226*reg91; reg231=reg231*reg91; reg205=reg205*reg91;
    reg240=reg240*reg91; reg257=reg257*reg91; reg250=reg250*reg91; reg193=reg193*reg91; reg248=reg248*reg91;
    reg237=reg237*reg91; reg225=reg225*reg91; reg253=reg253*reg91; reg239=reg239*reg91; reg264=reg264*reg91;
    reg256=reg91*reg256; reg211=reg211*reg91; reg266=reg91*reg266; reg265=reg91*reg265; reg209=reg209*reg91;
    reg251=reg91*reg251; reg269=reg269*reg91; reg246=reg91*reg246; reg184=reg184*reg91; reg242=reg91*reg242;
    reg221=reg221*reg91; reg238=reg91*reg238; reg261=reg261*reg91; reg262=reg91*reg262; reg224=reg224*reg91;
    reg232=reg91*reg232; reg186=reg186*reg91; T tmp_3_3=ponderation*reg257; T tmp_3_4=ponderation*reg261; T tmp_2_7=ponderation*reg253;
    T tmp_2_6=ponderation*reg193; T tmp_2_5=ponderation*reg205; T tmp_2_4=ponderation*reg194; T tmp_2_3=ponderation*reg228; T tmp_2_2=ponderation*reg224;
    T tmp_1_7=ponderation*reg221; T tmp_7_7=ponderation*reg209; T tmp_6_7=ponderation*reg211; T tmp_6_6=ponderation*reg184; T tmp_1_6=ponderation*reg256;
    T tmp_1_5=ponderation*reg266; T tmp_1_4=ponderation*reg265; T tmp_1_3=ponderation*reg251; T tmp_1_2=ponderation*reg246; T tmp_1_1=ponderation*reg242;
    T tmp_0_7=ponderation*reg238; T tmp_0_6=ponderation*reg262; T tmp_0_5=ponderation*reg232; T tmp_0_0=ponderation*reg218; T tmp_5_5=ponderation*reg214;
    T tmp_5_6=ponderation*reg216; T tmp_4_5=ponderation*reg222; T tmp_5_7=ponderation*reg226; T tmp_0_4=ponderation*reg240; T tmp_0_3=ponderation*reg250;
    T tmp_0_2=ponderation*reg248; T tmp_0_1=ponderation*reg225; T tmp_4_7=ponderation*reg239; T tmp_4_6=ponderation*reg237; T tmp_4_4=ponderation*reg231;
    T tmp_3_7=ponderation*reg186; T tmp_3_6=ponderation*reg269; T tmp_3_5=ponderation*reg264;
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
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=pow((*f.m).v2[1],2); T reg5=reg2*reg3; T reg6=pow((*f.m).v2[0],2);
    T reg7=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg8=1.0/(*f.m).elastic_modulus_3; T reg9=pow((*f.m).v1[1],2); T reg10=pow((*f.m).v1[0],2); T reg11=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg12=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg13=1.0/(*f.m).elastic_modulus_2; T reg14=pow((*f.m).v1[2],2); reg9=reg10+reg9; reg10=pow((*f.m).v2[2],2);
    T reg15=reg7*reg5; reg4=reg6+reg4; reg6=reg8*reg5; T reg16=reg11*reg5; reg10=reg4+reg10;
    reg14=reg9+reg14; reg4=reg11*reg15; reg9=reg12*reg6; T reg17=reg11*reg16; T reg18=reg13*reg6;
    reg17=reg18-reg17; reg10=pow(reg10,0.5); reg4=reg9+reg4; reg18=reg12*reg16; T reg19=1.0/(*f.m).elastic_modulus_1;
    reg14=pow(reg14,0.5); T reg20=reg13*reg15; T reg21=(*f.m).v1[2]/reg14; T reg22=(*f.m).v2[1]/reg10; T reg23=(*f.m).v2[2]/reg10;
    T reg24=(*f.m).v1[1]/reg14; T reg25=reg18+reg20; T reg26=reg12*reg4; T reg27=reg19*reg17; T reg28=reg7*reg3;
    T reg29=reg11*reg3; T reg30=reg7*reg25; reg14=(*f.m).v1[0]/reg14; T reg31=reg7*reg16; reg3=reg8*reg3;
    reg26=reg27-reg26; reg6=reg19*reg6; reg10=(*f.m).v2[0]/reg10; reg27=reg13*reg5; T reg32=reg21*reg22;
    T reg33=reg24*reg23; reg5=reg12*reg5; T reg34=reg7*reg15; T reg35=reg2*reg0; T reg36=reg14*reg23;
    T reg37=reg21*reg10; reg30=reg26-reg30; reg26=reg33-reg32; T reg38=2*reg14; reg34=reg6-reg34;
    reg6=reg7*reg27; reg15=reg12*reg15; T reg39=reg11*reg35; T reg40=reg7*reg5; reg16=reg19*reg16;
    T reg41=reg2*reg1; T reg42=2*reg10; T reg43=reg7*reg35; T reg44=reg13*reg3; reg31=reg9+reg31;
    reg3=reg12*reg3; reg9=reg11*reg29; reg35=reg8*reg35; T reg45=reg11*reg28; T reg46=reg37-reg36;
    T reg47=reg11*reg43; T reg48=reg24*reg38; T reg49=pow(reg10,2); T reg50=pow(reg14,2); T reg51=pow(reg24,2);
    T reg52=pow(reg22,2); T reg53=reg42*reg22; reg40=reg16+reg40; T reg54=2*reg26; T reg55=reg8*reg41;
    T reg56=reg11*reg39; T reg57=reg12*reg35; T reg58=reg14*reg22; T reg59=reg24*reg10; reg35=reg13*reg35;
    reg45=reg3+reg45; reg9=reg44-reg9; reg17=reg17/reg30; reg3=reg7*reg41; reg28=reg13*reg28;
    reg41=reg11*reg41; reg15=reg16+reg15; reg5=reg12*reg5; reg29=reg12*reg29; reg27=reg19*reg27;
    reg4=reg4/reg30; reg31=reg31/reg30; reg34=reg34/reg30; reg6=reg18+reg6; reg43=reg13*reg43;
    reg16=reg48*reg31; reg39=reg12*reg39; reg44=reg53*reg34; T reg60=reg11*reg3; reg5=reg27-reg5;
    reg6=reg6/reg30; reg25=reg25/reg30; reg27=reg58-reg59; T reg61=pow(reg21,2); T reg62=reg51*reg17;
    T reg63=reg52*reg4; reg15=reg15/reg30; T reg64=reg48*reg17; T reg65=reg53*reg4; reg40=reg40/reg30;
    T reg66=reg50*reg17; T reg67=reg11*reg41; T reg68=reg12*reg55; reg55=reg13*reg55; T reg69=reg54*reg46;
    T reg70=pow(reg46,2); T reg71=pow(reg26,2); T reg72=pow(reg23,2); T reg73=reg50*reg31; T reg74=reg49*reg34;
    reg47=reg57+reg47; reg56=reg35-reg56; reg35=reg51*reg31; reg57=reg52*reg34; T reg75=reg29+reg28;
    reg45=reg12*reg45; reg9=reg19*reg9; T reg76=reg49*reg4; reg3=reg13*reg3; T reg77=reg70*reg25;
    reg63=reg62+reg63; reg76=reg66+reg76; reg74=reg73+reg74; reg41=reg12*reg41; reg62=reg71*reg25;
    reg60=reg68+reg60; reg44=reg16+reg44; reg16=reg69*reg15; reg66=reg53*reg40; reg68=reg48*reg6;
    reg5=reg5/reg30; reg73=pow(reg27,2); T reg78=reg69*reg25; T reg79=reg50*reg6; T reg80=reg49*reg40;
    T reg81=reg51*reg6; T reg82=reg52*reg40; reg56=reg19*reg56; T reg83=reg39+reg43; reg65=reg64+reg65;
    reg57=reg35+reg57; reg35=reg71*reg15; reg47=reg12*reg47; reg75=reg7*reg75; reg64=reg70*reg15;
    reg67=reg55-reg67; reg45=reg9-reg45; reg9=reg72*reg4; reg55=reg61*reg17; T reg84=reg61*reg31;
    T reg85=reg72*reg34; T reg86=reg71*reg5; reg80=reg79+reg80; reg47=reg56-reg47; reg78=reg65+reg78;
    reg83=reg7*reg83; reg76=reg62+reg76; reg56=2*reg22; reg62=reg42*reg23; reg60=reg12*reg60;
    reg16=reg44+reg16; reg35=reg74+reg35; reg75=reg45-reg75; reg44=reg14*reg24; reg45=reg73*reg15;
    reg65=2*reg24; reg74=reg21*reg38; reg64=reg57+reg64; reg85=reg84+reg85; reg57=reg41+reg3;
    reg79=reg69*reg5; reg67=reg19*reg67; reg66=reg68+reg66; reg77=reg63+reg77; reg9=reg55+reg9;
    reg55=reg72*reg40; reg63=reg10*reg22; reg68=reg61*reg6; reg84=reg73*reg25; reg82=reg81+reg82;
    reg81=reg70*reg5; reg45=reg85+reg45; reg57=reg7*reg57; reg85=reg65*reg21; reg58=reg59+reg58;
    reg59=reg16*reg49; reg83=reg47-reg83; reg47=reg16*reg52; T reg87=2*reg21; T reg88=reg64*reg49;
    T reg89=reg24*reg26; T reg90=reg44*reg78; T reg91=reg14*reg46; T reg92=reg50*reg77; reg79=reg66+reg79;
    reg66=reg35*reg52; T reg93=reg46*reg26; reg75=reg75/reg30; reg84=reg9+reg84; reg9=reg50*reg78;
    T reg94=reg65*reg22; reg16=reg16*reg63; T reg95=reg44*reg77; reg54=reg54*reg27; reg86=reg80+reg86;
    reg80=reg64*reg63; reg81=reg82+reg81; reg82=2*reg46; reg55=reg68+reg55; reg68=reg74*reg17;
    T reg96=reg62*reg4; T reg97=reg73*reg5; T reg98=reg10*reg38; T reg99=reg51*reg76; T reg100=reg56*reg23;
    T reg101=reg14*reg10; T reg102=reg62*reg34; T reg103=reg24*reg22; reg60=reg67-reg60; reg67=reg74*reg31;
    reg77=reg51*reg77; reg64=reg64*reg52; T reg104=reg35*reg49; T reg105=reg50*reg76; reg78=reg51*reg78;
    T reg106=reg94*reg13; reg87=reg23*reg87; T reg107=reg54*reg25; reg59=reg9+reg59; reg9=reg79*reg71;
    T reg108=reg86*reg70; reg66=reg99+reg66; reg91=reg89+reg91; reg89=reg24*reg46; reg99=reg14*reg26;
    T reg109=reg21*reg23; T reg110=reg22*reg26; T reg111=reg10*reg46; reg82=reg82*reg27; T reg112=reg58*reg75;
    T reg113=reg103*reg75; T reg114=reg101*reg75; reg102=reg67+reg102; reg67=reg54*reg15; reg80=reg95+reg80;
    reg96=reg68+reg96; reg31=reg85*reg31; reg34=reg100*reg34; reg68=reg79*reg70; reg97=reg55+reg97;
    reg47=reg78+reg47; reg57=reg60-reg57; reg64=reg77+reg64; reg55=reg81*reg70; reg35=reg35*reg63;
    reg76=reg44*reg76; reg60=reg74*reg6; reg77=reg62*reg40; reg78=reg51*reg84; reg88=reg92+reg88;
    reg92=reg45*reg52; reg95=reg81*reg71; reg4=reg100*reg4; reg17=reg85*reg17; T reg115=reg98*reg12;
    T reg116=reg52*reg13; reg83=reg83/reg30; T reg117=reg49*reg12; reg79=reg79*reg93; reg16=reg90+reg16;
    reg90=reg94*reg12; T reg118=reg98*reg19; T reg119=reg52*reg12; T reg120=reg49*reg19; T reg121=reg86*reg71;
    reg104=reg105+reg104; reg81=reg81*reg93; reg105=reg50*reg84; T reg122=reg45*reg49; reg111=reg110+reg111;
    reg110=reg112*reg58; reg81=reg80+reg81; reg65=reg65*reg46; reg79=reg16+reg79; reg86=reg86*reg93;
    reg38=reg26*reg38; reg45=reg45*reg63; reg122=reg105+reg122; reg16=reg112*reg98; reg9=reg59+reg9;
    reg84=reg44*reg84; reg59=reg113*reg58; reg35=reg76+reg35; reg76=reg97*reg71; reg40=reg100*reg40;
    reg90=reg118-reg90; reg80=reg94*reg11; reg105=reg7*reg98; reg118=reg109*reg75; T reg123=reg99*reg83;
    T reg124=reg89*reg83; T reg125=reg21*reg27; T reg126=reg91*reg83; T reg127=reg52*reg11; T reg128=reg49*reg7;
    reg108=reg66+reg108; reg66=reg114*reg94; T reg129=reg72*reg11; reg117=reg116-reg117; reg55=reg64+reg55;
    reg121=reg104+reg121; reg64=reg114*reg98; reg104=reg87*reg11; reg115=reg106-reg115; reg30=reg57/reg30;
    reg95=reg88+reg95; reg57=reg113*reg98; reg113=reg113*reg94; reg112=reg112*reg94; reg68=reg47+reg68;
    reg47=reg72*reg7; reg119=reg120-reg119; reg88=reg97*reg70; reg92=reg78+reg92; reg4=reg17+reg4;
    reg25=reg82*reg25; reg17=reg10*reg26; reg78=reg22*reg46; reg67=reg102+reg67; reg34=reg31+reg34;
    reg15=reg82*reg15; reg107=reg96+reg107; reg77=reg60+reg77; reg31=reg54*reg5; reg60=reg7*reg87;
    reg6=reg85*reg6; reg96=(*f.m).alpha_2*reg52; reg102=reg51*reg12; reg113=reg55+reg113; reg76=reg122+reg76;
    reg55=(*f.m).alpha_1*reg51; reg106=reg126*reg91; reg116=(*f.m).alpha_2*reg49; reg110=reg79+reg110; reg79=reg72*reg8;
    reg120=(*f.m).alpha_1*reg50; reg127=reg128+reg127; reg80=reg105+reg80; reg87=reg87*reg8; reg105=reg67*reg52;
    reg122=reg51*reg107; reg128=reg126*reg65; reg112=reg68+reg112; reg68=reg118*reg94; reg88=reg92+reg88;
    reg25=reg4+reg25; reg15=reg34+reg15; reg31=reg77+reg31; reg56=reg56*reg46; reg40=reg6+reg40;
    reg5=reg82*reg5; reg42=reg42*reg26; reg4=reg125*reg83; reg6=reg17*reg30; reg34=reg78*reg30;
    reg97=reg97*reg93; reg45=reg84+reg45; reg77=reg111*reg30; reg84=reg124*reg91; reg59=reg81+reg59;
    reg19=reg50*reg19; reg64=reg121+reg64; reg114=reg114*reg58; reg86=reg35+reg86; reg35=reg123*reg38;
    reg57=reg95+reg57; reg81=reg124*reg38; reg124=reg124*reg65; reg92=reg50*reg107; reg95=reg21*reg26;
    reg121=reg67*reg49; T reg130=reg14*reg27; reg126=reg126*reg38; reg16=reg9+reg16; reg37=reg36+reg37;
    reg66=reg108+reg66; reg47=reg119-reg47; reg60=reg90-reg60; reg12=reg50*reg12; reg13=reg51*reg13;
    reg9=reg123*reg65; reg129=reg117-reg129; reg36=reg23*reg27; reg90=reg118*reg98; reg104=reg115-reg104;
    reg106=reg110+reg106; reg108=reg77*reg111; reg107=reg44*reg107; reg67=reg67*reg63; reg33=reg32+reg33;
    reg32=reg6*reg42; reg126=reg16+reg126; reg16=reg34*reg56; reg110=reg77*reg42; reg115=reg47*reg49;
    reg124=reg113+reg124; reg113=reg129*reg52; reg117=reg60*reg49; reg119=reg104*reg52; T reg131=reg37*reg75;
    T reg132=reg23*reg26; T reg133=reg47*reg50; T reg134=reg129*reg51; T reg135=reg10*reg27; T reg136=reg47*reg101;
    T reg137=reg129*reg103; reg116=reg120+reg116; reg120=(*f.m).alpha_3*reg71; T reg138=reg6*reg56; reg9=reg66+reg9;
    reg96=reg55+reg96; reg55=reg50*reg7; reg84=reg59+reg84; reg59=reg51*reg11; reg12=reg13-reg12;
    reg127=reg79-reg127; reg13=reg4*reg38; reg66=(*f.m).alpha_3*reg70; reg79=reg34*reg111; reg11=reg61*reg11;
    reg80=reg87-reg80; reg87=(*f.m).alpha_1*reg61; T reg139=(*f.m).alpha_2*reg72; reg123=reg123*reg91; T reg140=reg36*reg30;
    reg114=reg86+reg114; reg35=reg64+reg35; reg97=reg45+reg97; reg118=reg118*reg58; reg90=reg76+reg90;
    reg102=reg19-reg102; reg19=reg31*reg70; reg105=reg122+reg105; reg121=reg92+reg121; reg45=reg31*reg71;
    reg77=reg77*reg56; reg128=reg112+reg128; reg64=reg24*reg27; reg76=reg21*reg46; reg130=reg95+reg130;
    reg86=reg4*reg65; reg68=reg88+reg68; reg7=reg61*reg7; reg88=reg104*reg51; reg92=reg60*reg50;
    reg34=reg34*reg42; reg95=reg25*reg50; reg112=reg15*reg49; reg81=reg57+reg81; reg5=reg40+reg5;
    reg40=reg104*reg103; reg57=reg60*reg101; reg122=reg15*reg52; T reg141=reg25*reg51; reg8=reg61*reg8;
    reg59=reg55+reg59; reg123=reg114+reg123; reg55=reg80*reg61; reg114=reg131*reg94; reg88=reg92+reg88;
    reg32=reg35+reg32; reg34=reg81+reg34; reg79=reg84+reg79; reg35=reg63*reg2; reg16=reg124+reg16;
    reg122=reg141+reg122; reg81=reg58*reg2; reg84=reg5*reg70; reg86=reg68+reg86; reg68=reg140*reg56;
    reg11=reg12-reg11; reg6=reg6*reg111; reg12=(*f.m).alpha_2*reg63; reg138=reg9+reg138; reg13=reg90+reg13;
    reg9=(*f.m).alpha_1*reg44; reg90=reg140*reg42; reg77=reg128+reg77; reg92=(*f.m).alpha_3*reg73; reg139=reg87+reg139;
    reg120=reg116+reg120; reg64=reg76+reg64; reg66=reg96+reg66; reg19=reg105+reg19; reg76=reg130*reg83;
    reg25=reg25*reg44; reg63=reg15*reg63; reg15=reg131*reg98; reg45=reg121+reg45; reg87=reg80*reg109;
    reg113=reg115+reg113; reg96=reg127*reg72; reg40=reg57+reg40; reg119=reg117+reg119; reg57=reg80*reg72;
    reg110=reg126+reg110; reg75=reg33*reg75; reg105=reg22*reg27; reg134=reg133+reg134; reg115=reg127*reg61;
    reg116=reg23*reg46; reg135=reg132+reg135; reg117=reg127*reg109; reg137=reg136+reg137; reg121=reg5*reg71;
    reg112=reg95+reg112; reg118=reg97+reg118; reg31=reg31*reg93; reg67=reg107+reg67; reg14=reg14*reg21;
    reg4=reg4*reg91; reg108=reg106+reg108; reg10=reg10*reg23; reg7=reg102-reg7; reg95=reg108*reg34;
    reg97=reg79*reg77; reg6=reg123+reg6; reg121=reg112+reg121; reg105=reg116+reg105; reg114=reg19+reg114;
    reg19=reg79*reg110; reg102=reg76*reg65; reg94=reg75*reg94; reg84=reg122+reg84; reg106=reg76*reg38;
    reg107=reg108*reg16; reg83=reg64*reg83; reg68=reg86+reg68; reg15=reg45+reg15; reg45=reg7*reg49;
    reg140=reg140*reg111; reg4=reg118+reg4; reg86=reg11*reg52; reg90=reg13+reg90; reg13=reg135*reg30;
    reg98=reg75*reg98; reg96=reg113+reg96; reg112=reg35*reg53; reg26=reg27*reg26; reg21=reg24*reg21;
    reg57=reg119+reg57; reg24=reg81*reg53; reg23=reg22*reg23; reg22=(*f.m).alpha_3*reg93; reg12=reg9+reg12;
    reg92=reg139+reg92; reg9=reg120*reg138; reg113=reg66*reg16; reg59=reg8-reg59; reg8=reg10*reg1;
    reg116=reg37*reg1; reg118=reg7*reg50; reg119=reg11*reg51; reg31=reg67+reg31; reg131=reg131*reg58;
    reg63=reg25+reg63; reg25=reg81*reg48; reg55=reg88+reg55; reg93=reg5*reg93; reg5=reg35*reg48;
    reg67=reg81*reg58; reg115=reg134+reg115; reg117=reg137+reg117; reg88=reg35*reg58; reg87=reg40+reg87;
    reg40=reg32*reg120; reg122=reg66*reg34; reg10=(*f.m).alpha_2*reg10; reg123=(*f.m).alpha_1*reg14; reg38=reg83*reg38;
    reg98=reg121+reg98; reg122=reg40+reg122; reg93=reg63+reg93; reg75=reg75*reg58; reg40=reg120*reg6;
    reg63=reg66*reg79; reg119=reg118+reg119; reg61=reg59*reg61; reg112=reg96+reg112; reg96=reg8*reg62;
    reg118=reg13*reg42; reg106=reg15+reg106; reg15=reg7*reg101; reg121=reg11*reg103; reg46=reg27*reg46;
    reg24=reg57+reg24; reg27=reg116*reg62; reg57=reg8*reg74; reg5=reg115+reg5; reg19=reg95-reg19;
    reg95=reg116*reg37; reg67=reg87+reg67; reg87=reg34*reg77; reg22=reg12+reg22; reg12=reg110*reg16;
    reg88=reg117+reg88; reg115=reg8*reg37; reg117=reg116*reg74; reg10=reg123+reg10; reg26=(*f.m).alpha_3*reg26;
    reg30=reg105*reg30; reg25=reg55+reg25; reg97=reg107-reg97; reg2=reg44*reg2; reg44=(*f.m).alpha_1*reg21;
    reg55=(*f.m).alpha_2*reg23; reg113=reg9+reg113; reg23=reg23*reg0; reg9=reg33*reg0; reg90=reg92*reg90;
    reg102=reg114+reg102; reg94=reg84+reg94; reg72=reg59*reg72; reg86=reg45+reg86; reg65=reg83*reg65;
    reg68=reg92*reg68; reg76=reg76*reg91; reg45=reg13*reg56; reg140=reg4+reg140; reg131=reg31+reg131;
    reg90=reg122+reg90; reg4=reg9*reg85; reg45=reg102+reg45; reg63=reg40+reg63; reg31=reg22*reg77;
    reg117=reg25+reg117; reg25=reg22*reg110; reg95=reg67+reg95; reg42=reg30*reg42; reg68=reg113+reg68;
    reg38=reg98+reg38; reg40=reg9*reg33; reg57=reg5+reg57; reg5=reg23*reg85; reg140=reg92*reg140;
    reg118=reg106+reg118; reg67=reg23*reg100; reg96=reg112+reg96; reg83=reg83*reg91; reg75=reg93+reg75;
    reg84=reg2*reg53; reg72=reg86+reg72; reg13=reg13*reg111; reg76=reg131+reg76; reg1=reg14*reg1;
    reg14=reg32*reg97; reg86=reg19*reg138; reg26=reg10+reg26; reg55=reg44+reg55; reg56=reg30*reg56;
    reg65=reg94+reg65; reg46=(*f.m).alpha_3*reg46; reg12=reg87-reg12; reg121=reg15+reg121; reg10=reg23*reg33;
    reg109=reg59*reg109; reg15=reg9*reg100; reg27=reg24+reg27; reg115=reg88+reg115; reg61=reg119+reg61;
    reg24=reg2*reg48; reg44=reg108*reg138; reg42=reg38+reg42; reg15=reg27+reg15; reg118=reg26*reg118;
    reg25=reg90+reg25; reg24=reg61+reg24; reg56=reg65+reg56; reg74=reg1*reg74; reg10=reg115+reg10;
    reg46=reg55+reg46; reg86=reg14-reg86; reg4=reg117+reg4; reg14=reg12*reg6; reg67=reg96+reg67;
    reg30=reg30*reg111; reg83=reg75+reg83; reg140=reg63+reg140; reg40=reg95+reg40; reg5=reg57+reg5;
    reg62=reg1*reg62; reg84=reg72+reg84; reg27=reg22*reg108; reg109=reg121+reg109; reg45=reg26*reg45;
    reg31=reg68+reg31; reg13=reg76+reg13; reg38=reg2*reg58; reg0=reg21*reg0; reg21=reg6*reg110;
    reg55=reg32*reg108; reg57=reg6*reg77; reg61=reg6*reg16; reg63=reg32*reg79; reg21=reg55-reg21;
    reg30=reg83+reg30; reg38=reg109+reg38; reg37=reg1*reg37; reg55=reg6*reg34; reg65=reg32*reg77;
    reg68=reg110*reg138; reg118=reg25+reg118; reg42=reg46*reg42; reg13=reg26*reg13; reg27=reg140+reg27;
    reg45=reg31+reg45; reg56=reg46*reg56; reg62=reg84+reg62; reg100=reg0*reg100; reg25=reg4*reg10;
    reg57=reg44-reg57; reg74=reg24+reg74; reg24=reg5*reg40; reg31=reg15*reg10; reg44=reg79*reg138;
    reg72=reg67*reg40; reg85=reg0*reg85; reg14=reg86+reg14; reg75=1-var_inter[1]; reg42=reg118+reg42;
    reg85=reg74+reg85; reg74=reg34*reg138; reg68=reg65-reg68; reg65=reg32*reg16; reg55=reg63-reg55;
    reg97=reg97/reg14; reg57=reg57/reg14; reg21=reg21/reg14; reg19=reg19/reg14; reg61=reg44-reg61;
    reg31=reg72-reg31; reg44=1-var_inter[0]; reg100=reg62+reg100; reg62=reg4*reg67; reg63=reg5*reg15;
    reg25=reg24-reg25; reg37=reg38+reg37; reg33=reg0*reg33; reg56=reg45+reg56; reg13=reg27+reg13;
    reg30=reg46*reg30; reg24=elem.pos(1)[1]*var_inter[0]; reg27=reg25*reg100; reg68=reg68/reg14; reg12=reg12/reg14;
    reg38=reg75*elem.pos(1)[0]; reg30=reg13+reg30; reg13=reg75*elem.pos(0)[0]; reg45=elem.pos(0)[1]*reg44; reg72=reg75*elem.pos(0)[1];
    reg76=reg75*elem.pos(1)[1]; reg74=reg65-reg74; reg65=elem.pos(1)[0]*var_inter[0]; reg83=reg85*reg31; reg21=reg21*reg56;
    reg62=reg63-reg62; reg61=reg61/reg14; reg63=reg44*elem.pos(0)[0]; reg97=reg97*reg42; reg57=reg57*reg42;
    reg19=reg19*reg56; reg55=reg55/reg14; reg33=reg37+reg33; reg76=reg76-reg72; reg37=var_inter[1]*elem.pos(2)[0];
    reg38=reg38-reg13; reg56=reg55*reg56; reg42=reg61*reg42; reg57=reg21-reg57; reg68=reg68*reg30;
    reg21=var_inter[1]*elem.pos(2)[1]; reg55=reg63+reg65; reg61=elem.pos(2)[0]*var_inter[0]; reg84=reg45+reg24; reg86=elem.pos(2)[1]*var_inter[0];
    reg19=reg97-reg19; reg87=reg100*reg10; reg88=reg67*reg33; reg90=reg85*reg10; reg93=reg5*reg33;
    reg14=reg74/reg14; reg74=reg62*reg33; reg27=reg83-reg27; reg12=reg12*reg30; reg83=1-(*f.m).resolution;
    reg94=reg5*reg100; reg95=reg85*reg67; reg93=reg90-reg93; reg90=reg4*reg33; reg88=reg87-reg88;
    reg87=reg85*reg40; reg96=reg15*reg33; reg97=reg44*elem.pos(3)[1]; reg86=reg86-reg84; reg98=reg44*elem.pos(3)[0];
    reg61=reg61-reg55; reg102=var_inter[1]*elem.pos(3)[1]; reg30=reg14*reg30; reg56=reg42-reg56; reg74=reg27+reg74;
    reg68=reg57-reg68; reg14=reg100*reg40; reg37=reg38+reg37; reg12=reg19+reg12; reg19=var_inter[1]*elem.pos(3)[0];
    reg21=reg76+reg21; reg27=(*f.m).resolution*reg120; reg94=reg95-reg94; reg68=reg68*reg83; reg38=reg4*reg100;
    reg42=reg85*reg15; reg57=(*f.m).resolution*reg66; reg12=reg12*reg83; reg93=reg93/reg74; reg90=reg87-reg90;
    reg37=reg37-reg19; reg88=reg88/reg74; reg56=reg30+reg56; reg96=reg14-reg96; reg21=reg21-reg102;
    reg97=reg86+reg97; reg98=reg61+reg98; reg14=(*f.m).resolution*reg93; reg30=(*f.m).resolution*reg88; reg61=(*f.m).resolution*reg22;
    reg79=reg83*reg79; reg6=reg83*reg6; reg68=reg57+reg68; reg94=reg94/reg74; reg57=reg37*reg97;
    reg38=reg42-reg38; reg12=reg27+reg12; reg25=reg25/reg74; reg90=reg90/reg74; reg27=reg21*reg98;
    reg56=reg56*reg83; reg96=reg96/reg74; reg31=reg31/reg74; reg27=reg57-reg27; reg42=(*f.m).resolution*reg25;
    reg57=(*f.m).resolution*reg31; reg32=reg32*reg83; reg34=reg83*reg34; reg38=reg38/reg74; reg138=reg83*reg138;
    reg16=reg83*reg16; reg74=reg62/reg74; reg108=reg83*reg108; reg62=(*f.m).resolution*reg96; reg76=(*f.m).resolution*reg90;
    reg14=reg79-reg14; reg30=reg6+reg30; reg68=reg68*(*f.m).deltaT; reg6=(*f.m).resolution*reg94; reg12=reg12*(*f.m).deltaT;
    reg56=reg61+reg56; reg37=reg37/reg27; reg61=reg68*reg14; reg79=reg12*reg30; reg56=reg56*(*f.m).deltaT;
    reg21=reg21/reg27; reg98=reg98/reg27; reg62=reg138-reg62; reg42=reg34-reg42; reg32=reg57+reg32;
    reg76=reg16+reg76; reg16=(*f.m).resolution*reg38; reg6=reg108+reg6; reg77=reg83*reg77; reg110=reg83*reg110;
    reg97=reg97/reg27; reg34=(*f.m).resolution*reg74; reg57=var_inter[0]*reg37; reg83=var_inter[0]*reg21; reg86=reg56*reg6;
    reg87=reg44*reg21; reg95=var_inter[1]*reg98; reg106=reg75*reg97; reg107=reg79+reg61; reg108=var_inter[1]*reg97;
    reg109=reg75*reg98; reg112=reg12*reg62; reg113=reg68*reg76; reg114=reg12*reg32; reg115=reg68*reg42;
    reg16=reg77-reg16; reg110=reg34+reg110; reg34=reg44*reg37; reg77=reg107+reg86; reg117=reg112+reg113;
    reg118=reg56*reg16; reg119=reg34+reg95; reg121=reg87+reg108; reg122=reg109+reg57; reg123=reg115+reg114;
    reg124=reg56*reg110; reg126=reg106+reg83; reg128=reg117+reg118; reg131=2*reg77; reg132=reg123+reg124;
    reg133=reg109-reg34; reg134=reg108-reg83; reg136=reg57-reg95; reg137=0.5*reg126; reg138=reg87-reg106;
    reg139=0.5*reg122; reg140=0.5*reg119; reg141=0.5*reg121; T reg142=reg131*reg139; T reg143=reg132*reg126;
    T reg144=reg131*reg141; T reg145=reg128*reg119; T reg146=reg132*reg121; T reg147=reg131*reg137; T reg148=reg128*reg122;
    T reg149=0.5*reg134; T reg150=0.5*reg136; T reg151=0.5*reg138; T reg152=0.5*reg133; T reg153=reg131*reg140;
    T reg154=reg132*reg134; T reg155=reg131*reg150; T reg156=reg132*reg138; T reg157=reg142-reg143; T reg158=reg148-reg147;
    T reg159=reg131*reg152; T reg160=reg146-reg153; T reg161=reg144-reg145; T reg162=reg128*reg136; T reg163=reg131*reg149;
    T reg164=reg128*reg133; T reg165=reg131*reg151; reg161=reg27*reg161; T reg166=reg162+reg163; T reg167=reg156+reg159;
    reg158=reg27*reg158; reg160=reg27*reg160; T reg168=reg164+reg165; T reg169=reg154+reg155; reg157=reg27*reg157;
    reg158=ponderation*reg158; reg161=ponderation*reg161; reg157=ponderation*reg157; reg160=ponderation*reg160; T reg170=reg167*reg27;
    T reg171=reg27*reg169; T reg172=reg168*reg27; T reg173=reg27*reg166; T vec_3=-reg158; T vec_6=-reg160;
    reg158=ponderation*reg173; T vec_5=reg158; reg160=ponderation*reg170; T vec_0=reg160; T vec_2=-reg157;
    reg157=ponderation*reg171; T vec_4=reg157; T vec_7=-reg161; reg161=ponderation*reg172; T vec_1=reg161;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; T reg2=2*(*f.m).shear_modulus_12; reg1=1.0/reg1; reg0=1.0/reg0;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=pow((*f.m).v2[1],2); T reg5=pow((*f.m).v2[0],2); T reg6=reg2*reg3;
    T reg7=1.0/(*f.m).elastic_modulus_3; T reg8=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg9=pow((*f.m).v1[0],2); T reg10=pow((*f.m).v1[1],2); T reg11=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg12=pow((*f.m).v2[2],2); reg4=reg5+reg4; reg5=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; reg10=reg9+reg10; reg9=reg11*reg6;
    T reg13=reg7*reg6; T reg14=1.0/(*f.m).elastic_modulus_2; T reg15=reg8*reg6; T reg16=pow((*f.m).v1[2],2); reg16=reg10+reg16;
    reg10=reg14*reg13; reg12=reg4+reg12; reg4=reg11*reg15; T reg17=reg5*reg13; T reg18=reg11*reg9;
    T reg19=reg5*reg9; T reg20=reg14*reg15; T reg21=1.0/(*f.m).elastic_modulus_1; reg16=pow(reg16,0.5); reg12=pow(reg12,0.5);
    reg18=reg10-reg18; reg4=reg17+reg4; reg10=(*f.m).v2[1]/reg12; T reg22=(*f.m).v2[2]/reg12; T reg23=reg21*reg18;
    T reg24=reg19+reg20; T reg25=(*f.m).v1[2]/reg16; T reg26=(*f.m).v1[1]/reg16; T reg27=reg5*reg4; T reg28=reg8*reg3;
    T reg29=reg7*reg3; reg3=reg11*reg3; reg12=(*f.m).v2[0]/reg12; reg16=(*f.m).v1[0]/reg16; T reg30=reg8*reg15;
    T reg31=reg25*reg10; T reg32=reg26*reg22; T reg33=reg5*reg6; reg27=reg23-reg27; reg23=reg8*reg24;
    T reg34=reg8*reg9; T reg35=reg2*reg0; reg6=reg14*reg6; reg13=reg21*reg13; reg34=reg17+reg34;
    reg30=reg13-reg30; reg13=reg8*reg6; reg9=reg21*reg9; reg17=reg8*reg33; T reg36=reg11*reg35;
    T reg37=2*reg16; reg15=reg5*reg15; reg23=reg27-reg23; reg27=reg32-reg31; T reg38=reg25*reg12;
    T reg39=reg16*reg22; T reg40=reg7*reg35; T reg41=reg11*reg28; T reg42=2*reg12; T reg43=reg11*reg3;
    T reg44=reg14*reg29; reg29=reg5*reg29; T reg45=reg2*reg1; reg35=reg8*reg35; reg15=reg9+reg15;
    T reg46=pow(reg16,2); T reg47=pow(reg26,2); T reg48=reg7*reg45; T reg49=reg11*reg35; T reg50=reg11*reg36;
    reg17=reg9+reg17; reg9=reg5*reg40; reg40=reg14*reg40; reg41=reg29+reg41; reg18=reg18/reg23;
    reg6=reg21*reg6; reg43=reg44-reg43; reg34=reg34/reg23; reg30=reg30/reg23; reg4=reg4/reg23;
    reg13=reg19+reg13; reg29=reg26*reg12; reg3=reg5*reg3; reg44=reg11*reg45; reg28=reg14*reg28;
    T reg51=reg16*reg10; reg33=reg5*reg33; reg45=reg8*reg45; T reg52=reg26*reg37; T reg53=pow(reg12,2);
    T reg54=pow(reg10,2); T reg55=2*reg27; T reg56=reg42*reg10; T reg57=reg38-reg39; reg24=reg24/reg23;
    T reg58=pow(reg57,2); T reg59=reg3+reg28; T reg60=pow(reg27,2); reg41=reg5*reg41; T reg61=pow(reg22,2);
    reg43=reg21*reg43; reg36=reg5*reg36; T reg62=reg46*reg18; reg33=reg6-reg33; reg6=reg53*reg4;
    T reg63=reg56*reg4; T reg64=reg52*reg18; T reg65=reg54*reg4; T reg66=reg47*reg18; reg35=reg14*reg35;
    reg15=reg15/reg23; T reg67=reg11*reg45; T reg68=pow(reg25,2); T reg69=reg11*reg44; T reg70=reg5*reg48;
    reg48=reg14*reg48; reg49=reg9+reg49; reg9=reg55*reg57; reg50=reg40-reg50; reg40=reg51-reg29;
    reg17=reg17/reg23; T reg71=reg56*reg30; T reg72=reg47*reg34; T reg73=reg54*reg30; T reg74=reg53*reg30;
    T reg75=reg46*reg34; reg13=reg13/reg23; T reg76=reg52*reg34; reg44=reg5*reg44; T reg77=reg54*reg17;
    T reg78=reg58*reg15; T reg79=reg53*reg17; T reg80=reg68*reg34; T reg81=reg61*reg30; reg73=reg72+reg73;
    reg72=reg9*reg24; reg63=reg64+reg63; reg64=reg47*reg13; T reg82=pow(reg40,2); reg71=reg76+reg71;
    reg76=reg61*reg4; T reg83=reg68*reg18; T reg84=reg9*reg15; reg6=reg62+reg6; reg45=reg14*reg45;
    reg41=reg43-reg41; reg59=reg8*reg59; reg50=reg21*reg50; reg49=reg5*reg49; reg33=reg33/reg23;
    reg43=reg36+reg35; reg62=reg60*reg24; reg69=reg48-reg69; reg74=reg75+reg74; reg48=reg46*reg13;
    reg67=reg70+reg67; reg70=reg56*reg17; reg75=reg60*reg15; T reg85=reg52*reg13; reg65=reg66+reg65;
    reg66=reg58*reg24; T reg86=reg25*reg37; reg79=reg48+reg79; reg48=reg60*reg33; reg77=reg64+reg77;
    reg64=reg58*reg33; T reg87=reg68*reg13; T reg88=reg61*reg17; reg70=reg85+reg70; reg85=reg9*reg33;
    T reg89=reg16*reg26; reg59=reg41-reg59; reg49=reg50-reg49; reg43=reg8*reg43; reg69=reg21*reg69;
    reg67=reg5*reg67; reg41=reg44+reg45; reg75=reg74+reg75; reg66=reg65+reg66; reg76=reg83+reg76;
    reg50=reg82*reg24; reg72=reg63+reg72; reg78=reg73+reg78; reg81=reg80+reg81; reg63=reg82*reg15;
    reg65=2*reg10; reg73=reg42*reg22; reg84=reg71+reg84; reg71=2*reg26; reg6=reg62+reg6;
    reg62=reg12*reg10; reg74=reg16*reg12; reg80=2*reg57; reg55=reg55*reg40; reg83=reg89*reg72;
    T reg90=reg84*reg62; reg63=reg81+reg63; reg81=reg84*reg54; T reg91=reg65*reg22; T reg92=reg47*reg72;
    T reg93=reg46*reg66; T reg94=reg86*reg34; reg85=reg70+reg85; reg70=reg71*reg25; T reg95=reg73*reg30;
    reg43=reg49-reg43; reg49=reg47*reg6; T reg96=reg82*reg33; T reg97=reg89*reg66; reg88=reg87+reg88;
    reg64=reg77+reg64; reg77=reg78*reg62; reg50=reg76+reg50; reg76=reg75*reg54; reg72=reg46*reg72;
    reg48=reg79+reg48; reg84=reg84*reg53; reg66=reg47*reg66; reg79=reg57*reg27; reg87=reg86*reg18;
    T reg98=reg73*reg4; T reg99=reg75*reg53; T reg100=reg78*reg54; T reg101=reg46*reg6; reg41=reg8*reg41;
    reg59=reg59/reg23; reg78=reg78*reg53; reg67=reg69-reg67; reg69=reg26*reg10; T reg102=reg16*reg57;
    T reg103=reg26*reg27; reg51=reg29+reg51; reg29=reg64*reg60; reg84=reg72+reg84; reg72=reg85*reg60;
    T reg104=reg71*reg10; T reg105=reg25*reg22; reg80=reg80*reg40; reg98=reg87+reg98; reg87=reg64*reg58;
    reg30=reg91*reg30; reg34=reg70*reg34; T reg106=reg55*reg24; T reg107=reg55*reg15; reg95=reg94+reg95;
    reg94=reg12*reg37; reg78=reg93+reg78; reg93=reg16*reg27; reg100=reg66+reg100; reg66=reg26*reg57;
    reg6=reg89*reg6; reg90=reg83+reg90; reg18=reg70*reg18; reg75=reg75*reg62; reg83=reg85*reg79;
    reg4=reg91*reg4; reg102=reg103+reg102; reg103=reg47*reg50; reg96=reg88+reg96; reg64=reg64*reg79;
    reg77=reg97+reg77; reg88=reg86*reg13; reg97=reg73*reg17; T reg108=2*reg25; T reg109=reg74*reg59;
    reg81=reg92+reg81; reg92=reg69*reg59; T reg110=reg51*reg59; reg85=reg85*reg58; T reg111=reg63*reg54;
    T reg112=reg46*reg50; T reg113=reg63*reg53; T reg114=reg12*reg57; T reg115=reg10*reg27; reg99=reg101+reg99;
    reg101=reg48*reg60; reg41=reg67-reg41; reg67=reg48*reg58; reg43=reg43/reg23; reg76=reg49+reg76;
    reg85=reg81+reg85; reg48=reg48*reg79; reg111=reg103+reg111; reg108=reg22*reg108; reg49=reg25*reg40;
    reg81=reg110*reg104; reg67=reg76+reg67; reg76=reg94*reg21; reg103=reg110*reg94; reg23=reg41/reg23;
    reg41=reg109*reg104; reg75=reg6+reg75; reg113=reg112+reg113; reg6=reg96*reg60; reg112=reg96*reg58;
    reg37=reg27*reg37; reg106=reg98+reg106; reg98=reg12*reg27; T reg116=reg10*reg57; reg114=reg115+reg114;
    reg115=reg53*reg5; reg29=reg78+reg29; reg78=reg92*reg51; reg72=reg84+reg72; reg84=reg92*reg94;
    T reg117=reg104*reg5; reg63=reg63*reg62; reg107=reg95+reg107; reg95=reg53*reg21; reg110=reg110*reg51;
    T reg118=reg54*reg5; reg24=reg80*reg24; T reg119=reg109*reg94; reg101=reg99+reg101; reg30=reg34+reg30;
    reg15=reg80*reg15; reg97=reg88+reg97; reg92=reg92*reg104; reg87=reg100+reg87; reg50=reg89*reg50;
    reg34=reg55*reg33; reg13=reg70*reg13; reg17=reg91*reg17; reg71=reg71*reg57; reg88=reg104*reg14;
    reg4=reg18+reg4; reg83=reg90+reg83; reg64=reg77+reg64; reg18=reg94*reg5; reg77=reg105*reg59;
    reg90=reg93*reg43; reg99=reg66*reg43; reg100=reg102*reg43; T reg120=reg54*reg14; reg63=reg50+reg63;
    reg92=reg87+reg92; reg24=reg4+reg24; reg4=reg99*reg71; reg84=reg29+reg84; reg29=reg99*reg37;
    reg50=(*f.m).alpha_1*reg47; reg65=reg65*reg57; reg81=reg85+reg81; reg85=reg100*reg71; reg87=(*f.m).alpha_2*reg54;
    T reg121=reg47*reg106; T reg122=reg107*reg54; reg112=reg111+reg112; reg111=reg77*reg104; T reg123=reg61*reg8;
    reg119=reg101+reg119; reg34=reg97+reg34; reg17=reg13+reg17; reg33=reg80*reg33; reg99=reg99*reg102;
    reg13=reg49*reg43; reg97=(*f.m).alpha_2*reg53; reg109=reg109*reg51; reg48=reg75+reg48; reg75=reg98*reg23;
    reg101=reg116*reg23; T reg124=(*f.m).alpha_1*reg46; reg6=reg113+reg6; reg113=reg77*reg94; T reg125=reg114*reg23;
    T reg126=reg53*reg8; reg78=reg64+reg78; reg64=reg104*reg11; T reg127=reg8*reg94; T reg128=reg54*reg11;
    T reg129=reg90*reg71; reg41=reg67+reg41; reg96=reg96*reg79; reg42=reg42*reg27; reg67=reg22*reg40;
    T reg130=reg16*reg40; T reg131=reg25*reg27; reg38=reg39+reg38; reg39=reg107*reg53; T reg132=reg46*reg106;
    T reg133=reg90*reg37; T reg134=reg100*reg37; reg103=reg72+reg103; reg110=reg83+reg110; reg100=reg100*reg102;
    reg18=reg88-reg18; reg72=reg108*reg11; reg15=reg30+reg15; reg115=reg120-reg115; reg30=reg61*reg11;
    reg117=reg76-reg117; reg76=reg8*reg108; reg118=reg95-reg118; reg83=reg34*reg60; reg88=reg24*reg46;
    reg95=reg15*reg53; reg77=reg77*reg51; reg96=reg63+reg96; reg129=reg41+reg129; reg41=reg75*reg65;
    reg4=reg92+reg4; reg63=reg101*reg65; reg111=reg112+reg111; reg92=reg13*reg71; reg112=reg75*reg42;
    reg133=reg119+reg133; reg21=reg46*reg21; reg119=reg101*reg114; reg120=reg67*reg23; reg97=reg124+reg97;
    reg124=(*f.m).alpha_3*reg60; T reg135=reg38*reg59; reg33=reg17+reg33; reg85=reg81+reg85; reg17=reg125*reg65;
    reg122=reg121+reg122; reg81=reg34*reg58; reg87=reg50+reg87; reg50=(*f.m).alpha_3*reg58; reg121=(*f.m).alpha_1*reg68;
    T reg136=(*f.m).alpha_2*reg61; T reg137=reg24*reg47; T reg138=reg15*reg54; T reg139=reg26*reg40; T reg140=reg25*reg57;
    reg130=reg131+reg130; reg32=reg31+reg32; reg100=reg110+reg100; reg31=reg125*reg114; reg106=reg89*reg106;
    reg107=reg107*reg62; reg72=reg18-reg72; reg30=reg115-reg30; reg14=reg47*reg14; reg18=reg46*reg5;
    reg76=reg117-reg76; reg123=reg118-reg123; reg5=reg47*reg5; reg39=reg132+reg39; reg125=reg125*reg42;
    reg134=reg103+reg134; reg99=reg78+reg99; reg78=reg61*reg7; reg128=reg126+reg128; reg108=reg108*reg7;
    reg64=reg127+reg64; reg103=reg13*reg37; reg113=reg6+reg113; reg101=reg101*reg42; reg29=reg84+reg29;
    reg6=reg22*reg27; reg84=reg12*reg40; reg109=reg48+reg109; reg90=reg90*reg102; reg16=reg16*reg25;
    reg12=reg12*reg22; reg48=reg10*reg40; reg110=reg22*reg57; reg84=reg6+reg84; reg6=reg72*reg54;
    reg115=reg76*reg53; reg117=reg30*reg54; reg118=reg123*reg53; reg13=reg13*reg102; reg77=reg96+reg77;
    reg119=reg99+reg119; reg124=reg97+reg124; reg50=reg87+reg50; reg136=reg121+reg136; reg87=(*f.m).alpha_3*reg82;
    reg96=(*f.m).alpha_1*reg89; reg97=(*f.m).alpha_2*reg62; reg63=reg4+reg63; reg24=reg24*reg89; reg15=reg15*reg62;
    reg41=reg129+reg41; reg75=reg75*reg114; reg90=reg109+reg90; reg4=reg33*reg60; reg99=reg72*reg47;
    reg109=reg76*reg46; reg95=reg88+reg95; reg88=reg123*reg74; reg121=reg30*reg69; reg126=reg30*reg47;
    reg127=reg123*reg46; reg129=reg76*reg74; reg131=reg72*reg69; reg132=reg135*reg94; reg101=reg29+reg101;
    reg83=reg39+reg83; reg125=reg134+reg125; reg103=reg113+reg103; reg64=reg108-reg64; reg128=reg78-reg128;
    reg29=reg46*reg8; reg39=reg120*reg42; reg78=reg47*reg11; reg108=reg33*reg58; reg138=reg137+reg138;
    reg113=reg135*reg104; reg81=reg122+reg81; reg17=reg85+reg17; reg85=reg120*reg65; reg139=reg140+reg139;
    reg59=reg32*reg59; reg122=reg130*reg43; reg31=reg100+reg31; reg107=reg106+reg107; reg112=reg133+reg112;
    reg92=reg111+reg92; reg8=reg68*reg8; reg5=reg21-reg5; reg11=reg68*reg11; reg34=reg34*reg79;
    reg18=reg14-reg18; reg14=reg128*reg61; reg117=reg118+reg117; reg87=reg136+reg87; reg43=reg139*reg43;
    reg39=reg103+reg39; reg21=(*f.m).alpha_3*reg79; reg104=reg59*reg104; reg108=reg138+reg108; reg100=reg122*reg71;
    reg113=reg81+reg113; reg81=reg122*reg37; reg94=reg59*reg94; reg132=reg83+reg132; reg85=reg92+reg85;
    reg83=reg84*reg23; reg13=reg77+reg13; reg4=reg95+reg4; reg97=reg96+reg97; reg120=reg120*reg114;
    reg75=reg90+reg75; reg135=reg135*reg51; reg34=reg107+reg34; reg8=reg5-reg8; reg11=reg18-reg11;
    reg48=reg110+reg48; reg5=reg31*reg63; reg18=reg31*reg101; reg77=reg119*reg17; reg90=reg119*reg125;
    reg92=reg112*reg124; reg95=reg50*reg101; reg96=reg124*reg41; reg103=reg50*reg63; reg22=reg10*reg22;
    reg25=reg26*reg25; reg27=reg40*reg27; reg10=(*f.m).alpha_1*reg16; reg26=(*f.m).alpha_2*reg12; reg7=reg68*reg7;
    reg6=reg115+reg6; reg106=reg64*reg61; reg78=reg29+reg78; reg62=reg62*reg2; reg29=reg51*reg2;
    reg107=reg64*reg105; reg131=reg129+reg131; reg126=reg127+reg126; reg110=reg128*reg68; reg111=reg128*reg105;
    reg121=reg88+reg121; reg99=reg109+reg99; reg88=reg64*reg68; reg79=reg33*reg79; reg15=reg24+reg15;
    reg71=reg43*reg71; reg24=1-var_inter[1]; reg33=reg83*reg65; reg100=reg113+reg100; reg37=reg43*reg37;
    reg94=reg4+reg94; reg4=reg83*reg42; reg81=reg132+reg81; reg109=reg29*reg51; reg107=reg131+reg107;
    reg113=reg62*reg51; reg111=reg121+reg111; reg59=reg59*reg51; reg79=reg15+reg79; reg122=reg122*reg102;
    reg135=reg34+reg135; reg15=1-var_inter[0]; reg57=reg40*reg57; reg34=reg29*reg56; reg106=reg6+reg106;
    reg6=reg62*reg56; reg14=reg117+reg14; reg40=reg11*reg54; reg115=reg8*reg53; reg120=reg13+reg120;
    reg77=reg5-reg77; reg90=reg18-reg90; reg5=reg101*reg17; reg13=reg125*reg63; reg95=reg92+reg95;
    reg18=(*f.m).alpha_2*reg22; reg92=(*f.m).alpha_1*reg25; reg39=reg87*reg39; reg27=(*f.m).alpha_3*reg27; reg26=reg10+reg26;
    reg21=reg97+reg21; reg103=reg96+reg103; reg85=reg87*reg85; reg10=reg8*reg46; reg96=reg11*reg47;
    reg97=reg124*reg75; reg23=reg48*reg23; reg117=reg50*reg119; reg104=reg108+reg104; reg108=reg29*reg52;
    reg88=reg99+reg88; reg78=reg7-reg78; reg7=reg38*reg1; reg12=reg12*reg1; reg99=reg62*reg52;
    reg110=reg126+reg110; reg6=reg14+reg6; reg117=reg97+reg117; reg120=reg87*reg120; reg96=reg10+reg96;
    reg10=reg12*reg73; reg14=reg24*elem.pos(1)[1]; reg34=reg106+reg34; reg97=reg32*reg0; reg22=reg22*reg0;
    reg2=reg89*reg2; reg89=reg7*reg73; reg106=reg15*elem.pos(0)[0]; reg118=elem.pos(1)[0]*var_inter[0]; reg121=elem.pos(0)[1]*reg15;
    reg126=elem.pos(1)[1]*var_inter[0]; reg122=reg135+reg122; reg83=reg83*reg114; reg127=reg7*reg86; reg108=reg88+reg108;
    reg59=reg79+reg59; reg43=reg43*reg102; reg79=reg8*reg74; reg88=reg11*reg69; reg129=reg12*reg86;
    reg99=reg110+reg99; reg113=reg111+reg113; reg110=reg12*reg38; reg111=reg112*reg77; reg131=reg90*reg41;
    reg13=reg5-reg13; reg39=reg95+reg39; reg5=reg21*reg125; reg85=reg103+reg85; reg95=reg21*reg17;
    reg109=reg107+reg109; reg103=reg7*reg38; reg68=reg78*reg68; reg33=reg100+reg33; reg57=(*f.m).alpha_3*reg57;
    reg61=reg78*reg61; reg40=reg115+reg40; reg18=reg92+reg18; reg27=reg26+reg27; reg65=reg23*reg65;
    reg42=reg23*reg42; reg37=reg94+reg37; reg71=reg104+reg71; reg26=reg24*elem.pos(1)[0]; reg92=reg24*elem.pos(0)[0];
    reg4=reg81+reg4; reg81=reg24*elem.pos(0)[1]; reg94=reg22*reg32; reg131=reg111-reg131; reg57=reg18+reg57;
    reg18=reg13*reg75; reg100=reg31*reg41; reg104=reg75*reg17; reg1=reg16*reg1; reg16=reg2*reg56;
    reg61=reg40+reg61; reg40=reg112*reg31; reg107=reg2*reg52; reg103=reg109+reg103; reg109=reg75*reg125;
    reg111=var_inter[1]*elem.pos(2)[0]; reg26=reg26-reg92; reg115=reg97*reg32; reg68=reg96+reg68; reg5=reg39+reg5;
    reg4=reg27*reg4; reg33=reg27*reg33; reg39=reg21*reg31; reg120=reg117+reg120; reg95=reg85+reg95;
    reg10=reg6+reg10; reg6=reg22*reg91; reg14=reg14-reg81; reg85=var_inter[1]*elem.pos(2)[1]; reg89=reg34+reg89;
    reg34=reg97*reg91; reg65=reg71+reg65; reg42=reg37+reg42; reg83=reg122+reg83; reg37=reg97*reg70;
    reg127=reg108+reg127; reg43=reg59+reg43; reg23=reg23*reg114; reg59=reg106+reg118; reg71=elem.pos(2)[0]*var_inter[0];
    reg96=elem.pos(2)[1]*var_inter[0]; reg108=reg121+reg126; reg110=reg113+reg110; reg129=reg99+reg129; reg99=reg22*reg70;
    reg88=reg79+reg88; reg105=reg78*reg105; reg99=reg129+reg99; reg86=reg1*reg86; reg65=reg57*reg65;
    reg33=reg95+reg33; reg18=reg131+reg18; reg79=reg75*reg63; reg42=reg57*reg42; reg37=reg127+reg37;
    reg4=reg5+reg4; reg23=reg43+reg23; reg5=reg2*reg51; reg43=reg112*reg119; reg95=reg119*reg41;
    reg105=reg88+reg105; reg104=reg100-reg104; reg88=reg125*reg41; reg100=reg112*reg17; reg113=reg75*reg101;
    reg109=reg40-reg109; reg94=reg110+reg94; reg34=reg89+reg34; reg40=var_inter[1]*elem.pos(3)[0]; reg111=reg26+reg111;
    reg0=reg25*reg0; reg16=reg61+reg16; reg25=var_inter[1]*elem.pos(3)[1]; reg73=reg1*reg73; reg85=reg14+reg85;
    reg83=reg27*reg83; reg39=reg120+reg39; reg71=reg71-reg59; reg107=reg68+reg107; reg115=reg103+reg115;
    reg96=reg96-reg108; reg6=reg10+reg6; reg10=reg15*elem.pos(3)[1]; reg14=reg15*elem.pos(3)[0]; reg26=reg15*vectors[0][indices[0]+0];
    reg109=reg109/reg18; reg91=reg0*reg91; reg61=reg24*vectors[0][indices[1]+0]; reg90=reg90/reg18; reg79=reg95-reg79;
    reg38=reg1*reg38; reg73=reg16+reg73; reg16=reg6*reg115; reg68=reg34*reg94; reg10=reg96+reg10;
    reg5=reg105+reg5; reg104=reg104/reg18; reg89=reg99*reg115; reg14=reg71+reg14; reg77=reg77/reg18;
    reg71=reg24*vectors[0][indices[0]+0]; reg95=var_inter[0]*vectors[0][indices[1]+1]; reg42=reg4+reg42; reg111=reg111-reg40; reg83=reg39+reg83;
    reg4=reg15*vectors[0][indices[0]+1]; reg23=reg57*reg23; reg39=reg24*vectors[0][indices[0]+1]; reg96=var_inter[0]*vectors[0][indices[1]+0]; reg103=reg101*reg41;
    reg86=reg107+reg86; reg65=reg33+reg65; reg88=reg100-reg88; reg33=reg112*reg63; reg85=reg85-reg25;
    reg100=reg37*reg94; reg70=reg0*reg70; reg113=reg43-reg113; reg43=reg24*vectors[0][indices[1]+1]; reg38=reg5+reg38;
    reg95=reg4+reg95; reg4=var_inter[1]*vectors[0][indices[2]+0]; reg32=reg0*reg32; reg90=reg90*reg65; reg5=reg85*reg14;
    reg109=reg109*reg65; reg26=reg96+reg26; reg91=reg73+reg91; reg71=reg61-reg71; reg104=reg104*reg42;
    reg61=var_inter[0]*vectors[0][indices[2]+0]; reg39=reg43-reg39; reg43=var_inter[1]*vectors[0][indices[2]+1]; reg73=var_inter[0]*vectors[0][indices[2]+1]; reg23=reg83+reg23;
    reg103=reg33-reg103; reg88=reg88/reg18; reg13=reg13/reg18; reg113=reg113/reg18; reg68=reg16-reg68;
    reg70=reg86+reg70; reg79=reg79/reg18; reg16=reg111*reg10; reg77=reg77*reg42; reg100=reg89-reg100;
    reg33=reg99*reg34; reg83=reg37*reg6; reg90=reg77-reg90; reg42=reg79*reg42; reg77=reg15*vectors[0][indices[3]+1];
    reg95=reg73-reg95; reg13=reg13*reg23; reg26=reg61-reg26; reg32=reg38+reg32; reg18=reg103/reg18;
    reg83=reg33-reg83; reg65=reg113*reg65; reg33=reg100*reg91; reg38=var_inter[1]*vectors[0][indices[3]+1]; reg61=var_inter[1]*vectors[0][indices[3]+0];
    reg43=reg39+reg43; reg71=reg4+reg71; reg88=reg88*reg23; reg4=reg70*reg68; reg39=reg15*vectors[0][indices[3]+0];
    reg5=reg16-reg5; reg104=reg109-reg104; reg23=reg18*reg23; reg88=reg104-reg88; reg13=reg90+reg13;
    reg111=reg111/reg5; reg77=reg95+reg77; reg65=reg42-reg65; reg16=reg6*reg32; reg18=reg70*reg94;
    reg42=reg91*reg94; reg61=reg71-reg61; reg71=reg99*reg32; reg73=1-(*f.m).resolution; reg38=reg43-reg38;
    reg10=reg10/reg5; reg14=reg14/reg5; reg85=reg85/reg5; reg39=reg26+reg39; reg26=reg83*reg32;
    reg33=reg4-reg33; reg16=reg42-reg16; reg4=reg91*reg115; reg42=reg37*reg32; reg26=reg33+reg26;
    reg71=reg18-reg71; reg18=reg70*reg6; reg33=reg99*reg91; reg43=(*f.m).resolution*reg124; reg65=reg23+reg65;
    reg23=(*f.m).resolution*reg50; reg13=reg13*reg73; reg88=reg88*reg73; reg79=reg61*reg14; reg86=reg111*reg39;
    reg89=reg77*reg85; reg90=reg38*reg10; reg95=reg70*reg115; reg96=reg34*reg32; reg89=reg90-reg89;
    reg79=reg86-reg79; reg38=reg38*reg14; reg61=reg61*reg10; reg65=reg65*reg73; reg77=reg77*reg111;
    reg86=(*f.m).resolution*reg21; reg88=reg23+reg88; reg13=reg43+reg13; reg33=reg18-reg33; reg16=reg16/reg26;
    reg96=reg4-reg96; reg42=reg95-reg42; reg71=reg71/reg26; reg39=reg85*reg39; reg4=reg70*reg34;
    reg18=reg37*reg91; reg100=reg100/reg26; reg39=reg61-reg39; reg38=reg77-reg38; reg42=reg42/reg26;
    reg18=reg4-reg18; reg79=reg89+reg79; reg88=reg88*(*f.m).deltaT; reg4=(*f.m).resolution*reg71; reg23=(*f.m).resolution*reg16;
    reg13=reg13*(*f.m).deltaT; reg119=reg73*reg119; reg65=reg86+reg65; reg33=reg33/reg26; reg68=reg68/reg26;
    reg75=reg73*reg75; reg96=reg96/reg26; reg23=reg75+reg23; reg4=reg119-reg4; reg39=reg39-reg13;
    reg65=reg65*(*f.m).deltaT; reg79=0.5*reg79; reg43=(*f.m).resolution*reg33; reg101=reg73*reg101; reg41=reg73*reg41;
    reg63=reg73*reg63; reg83=reg83/reg26; reg112=reg112*reg73; reg61=(*f.m).resolution*reg42; reg75=(*f.m).resolution*reg96;
    reg31=reg73*reg31; reg38=reg38-reg88; reg26=reg18/reg26; reg18=(*f.m).resolution*reg68; reg77=(*f.m).resolution*reg100;
    reg79=reg79-reg65; reg86=(*f.m).resolution*reg83; reg43=reg31+reg43; reg61=reg63+reg61; reg75=reg41-reg75;
    reg77=reg101-reg77; reg112=reg18+reg112; reg125=reg73*reg125; reg18=(*f.m).resolution*reg26; reg17=reg73*reg17;
    reg31=reg23*reg39; reg41=reg4*reg38; reg63=reg15*reg85; reg73=reg15*reg111; reg89=reg61*reg38;
    reg90=reg75*reg39; reg39=reg112*reg39; reg95=reg24*reg14; reg101=reg24*reg10; reg103=var_inter[0]*reg85;
    reg104=var_inter[1]*reg10; reg18=reg17-reg18; reg125=reg86+reg125; reg17=var_inter[1]*reg14; reg41=reg31+reg41;
    reg38=reg77*reg38; reg31=reg43*reg79; reg86=var_inter[0]*reg111; reg105=reg95+reg86; reg107=reg125*reg79;
    reg89=reg90+reg89; reg79=reg18*reg79; reg90=reg101+reg103; reg109=reg63-reg101; reg31=reg41+reg31;
    reg41=reg95-reg73; reg110=reg86-reg17; reg39=reg38+reg39; reg38=reg104-reg103; reg113=reg73+reg17;
    reg117=reg63+reg104; reg31=2*reg31; reg119=0.5*reg90; reg120=0.5*reg41; reg122=0.5*reg110;
    reg79=reg89+reg79; reg89=0.5*reg105; reg107=reg39+reg107; reg39=0.5*reg38; reg127=0.5*reg113;
    reg129=0.5*reg109; reg131=0.5*reg117; reg132=reg39*reg31; reg133=reg110*reg79; reg134=reg122*reg31;
    reg135=reg38*reg107; reg136=reg120*reg31; reg137=reg119*reg31; reg138=reg109*reg107; reg140=reg131*reg31;
    T reg141=reg113*reg79; T reg142=reg105*reg79; T reg143=reg127*reg31; T reg144=reg129*reg31; T reg145=reg41*reg79;
    T reg146=reg90*reg107; T reg147=reg89*reg31; T reg148=reg117*reg107; reg144=reg145+reg144; reg146=reg146-reg147;
    reg141=reg141-reg140; reg134=reg135+reg134; reg132=reg133+reg132; reg137=reg137-reg142; reg143=reg143-reg148;
    reg136=reg138+reg136; reg136=reg5*reg136; reg134=reg134*reg5; reg144=reg144*reg5; reg137=reg137*reg5;
    reg143=reg5*reg143; reg132=reg132*reg5; reg141=reg5*reg141; reg146=reg146*reg5; T vec_0=ponderation*reg136;
    T vec_5=ponderation*reg132; T vec_4=ponderation*reg134; T vec_1=ponderation*reg144; T vec_3=ponderation*reg137; T vec_7=ponderation*reg141;
    T vec_6=ponderation*reg143; T vec_2=ponderation*reg146;
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

