/// @author hugo LECLERC

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
  static const unsigned order_integration = 2;
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
    static const bool has_elementary_matrix = true;
    static const bool has_skin_elementary_matrix = false;
    template<class TE,class TF, class TVEVE> static void after_solve(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
    T reg0=2*(*f.m).shear_modulus_13; T reg1=2*(*f.m).shear_modulus_23; reg0=1.0/reg0; T reg2=2*(*f.m).shear_modulus_12; reg1=1.0/reg1;
    T reg3=reg0*reg1; T reg4=pow((*f.m).v1[1],2); T reg5=pow((*f.m).v1[0],2); reg2=1.0/reg2; T reg6=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg7=pow((*f.m).v1[2],2); T reg8=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg9=pow((*f.m).v2[0],2); T reg10=pow((*f.m).v2[1],2); T reg11=1.0/(*f.m).elastic_modulus_3;
    T reg12=reg2*reg3; reg4=reg5+reg4; reg5=reg11*reg12; T reg13=reg8*reg12; T reg14=reg6*reg12;
    reg10=reg9+reg10; reg9=pow((*f.m).v2[2],2); T reg15=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg16=1.0/(*f.m).elastic_modulus_2; reg4=reg7+reg4;
    reg7=reg16*reg5; T reg17=reg8*reg13; T reg18=reg15*reg5; reg4=pow(reg4,0.5); T reg19=reg8*reg14;
    reg9=reg10+reg9; reg10=(*f.m).v1[0]/reg4; T reg20=reg16*reg14; T reg21=(*f.m).v1[1]/reg4; T reg22=reg15*reg13;
    reg19=reg18+reg19; reg9=pow(reg9,0.5); reg17=reg7-reg17; reg7=1.0/(*f.m).elastic_modulus_1; T reg23=2*reg21;
    reg4=(*f.m).v1[2]/reg4; T reg24=(*f.m).v2[0]/reg9; T reg25=reg7*reg17; T reg26=(*f.m).v2[1]/reg9; reg9=(*f.m).v2[2]/reg9;
    T reg27=reg15*reg19; T reg28=2*reg10; T reg29=reg22+reg20; T reg30=reg23*reg26; T reg31=2*reg4;
    T reg32=reg24*reg28; T reg33=reg4*reg26; T reg34=reg21*reg9; T reg35=reg6*reg3; T reg36=reg2*reg1;
    reg27=reg25-reg27; reg25=pow(reg24,2); T reg37=pow(reg26,2); T reg38=reg6*reg13; T reg39=reg8*reg3;
    reg3=reg11*reg3; T reg40=reg6*reg29; T reg41=reg16*reg12; reg5=reg7*reg5; reg12=reg15*reg12;
    T reg42=reg6*reg14; T reg43=reg34-reg33; T reg44=reg8*reg36; T reg45=reg6*reg36; T reg46=reg2*reg0;
    T reg47=reg15*reg3; T reg48=reg8*reg39; T reg49=reg8*reg35; reg36=reg11*reg36; T reg50=reg30*reg15;
    T reg51=reg32*reg7; T reg52=reg37*reg15; T reg53=reg25*reg7; T reg54=reg4*reg24; T reg55=reg10*reg9;
    reg38=reg18+reg38; reg18=reg30*reg16; T reg56=reg6*reg41; T reg57=reg32*reg15; T reg58=reg31*reg9;
    T reg59=reg37*reg16; T reg60=reg25*reg15; reg42=reg5-reg42; reg13=reg7*reg13; reg40=reg27-reg40;
    reg3=reg16*reg3; reg5=reg6*reg12; reg27=pow(reg9,2); T reg61=2*reg24; reg14=reg15*reg14;
    reg12=reg15*reg12; T reg62=reg6*reg46; reg56=reg22+reg56; T reg63=reg8*reg46; reg35=reg16*reg35;
    reg60=reg59-reg60; reg42=reg42/reg40; reg39=reg15*reg39; reg59=reg27*reg8; T reg64=reg21*reg28;
    T reg65=reg27*reg6; T reg66=pow(reg21,2); T reg67=pow(reg10,2); T reg68=reg54-reg55; T reg69=reg61*reg26;
    T reg70=reg10*reg26; T reg71=reg21*reg24; reg52=reg53-reg52; reg50=reg51-reg50; reg46=reg11*reg46;
    reg51=reg37*reg8; reg53=reg6*reg32; reg5=reg13+reg5; T reg72=reg30*reg8; reg19=reg19/reg40;
    T reg73=reg6*reg58; T reg74=reg8*reg45; T reg75=reg25*reg6; reg13=reg14+reg13; reg14=reg8*reg44;
    T reg76=reg15*reg36; reg36=reg16*reg36; reg17=reg17/reg40; reg38=reg38/reg40; reg57=reg18-reg57;
    reg41=reg7*reg41; reg49=reg47+reg49; reg18=2*reg43; reg47=reg58*reg8; reg48=reg3-reg48;
    reg5=reg5/reg40; reg3=reg15*reg66; T reg77=reg38*reg67; T reg78=pow(reg43,2); T reg79=reg17*reg67;
    T reg80=pow(reg68,2); T reg81=reg18*reg68; T reg82=reg25*reg42; T reg83=reg10*reg24; T reg84=reg21*reg26;
    T reg85=reg16*reg46; reg46=reg15*reg46; T reg86=reg8*reg63; reg74=reg76+reg74; reg14=reg36-reg14;
    reg36=reg8*reg62; reg76=reg39+reg35; reg49=reg15*reg49; reg48=reg7*reg48; reg45=reg16*reg45;
    reg65=reg52-reg65; reg44=reg15*reg44; reg52=reg17*reg66; T reg87=reg37*reg19; T reg88=pow(reg4,2);
    T reg89=reg70-reg71; T reg90=reg17*reg64; T reg91=reg69*reg19; reg12=reg41-reg12; reg13=reg13/reg40;
    reg41=reg15*reg67; reg29=reg29/reg40; reg59=reg60-reg59; reg73=reg50-reg73; reg50=reg25*reg19;
    reg47=reg57-reg47; reg57=reg69*reg42; reg60=reg38*reg64; reg56=reg56/reg40; T reg92=reg7*reg67;
    T reg93=reg27*reg11; reg51=reg75+reg51; reg75=reg58*reg11; T reg94=reg38*reg66; T reg95=reg37*reg42;
    reg72=reg53+reg72; reg53=reg16*reg66; T reg96=reg24*reg26; T reg97=reg80*reg29; T reg98=reg56*reg66;
    T reg99=reg25*reg5; reg87=reg52+reg87; reg52=reg56*reg67; T reg100=reg17*reg88; T reg101=reg37*reg5;
    T reg102=reg37*reg47; reg63=reg15*reg63; reg79=reg50+reg79; reg50=reg78*reg29; reg41=reg53-reg41;
    reg62=reg16*reg62; reg36=reg46+reg36; reg49=reg48-reg49; reg46=reg56*reg64; reg48=reg69*reg5;
    reg76=reg6*reg76; reg53=reg6*reg67; T reg103=reg83*reg65; reg86=reg85-reg86; reg85=reg8*reg66;
    reg14=reg7*reg14; T reg104=reg84*reg59; reg51=reg93-reg51; reg74=reg15*reg74; reg93=reg44+reg45;
    reg72=reg75-reg72; reg75=reg81*reg29; T reg105=reg66*reg59; T reg106=pow(reg89,2); reg12=reg12/reg40;
    T reg107=reg81*reg13; reg57=reg60+reg57; reg3=reg92-reg3; reg60=reg6*reg88; reg92=reg67*reg73;
    T reg108=reg66*reg47; T reg109=reg27*reg42; T reg110=reg38*reg88; T reg111=reg80*reg13; reg95=reg94+reg95;
    reg82=reg77+reg82; reg77=reg78*reg13; reg70=reg71+reg70; reg71=reg4*reg9; reg94=reg25*reg73;
    T reg112=reg84*reg47; T reg113=reg83*reg73; T reg114=reg27*reg19; T reg115=reg8*reg88; T reg116=reg67*reg65;
    reg91=reg90+reg91; reg90=reg25*reg65; T reg117=reg37*reg59; reg112=reg113+reg112; reg113=reg80*reg12;
    reg97=reg87+reg97; reg77=reg82+reg77; reg75=reg91+reg75; reg101=reg98+reg101; reg111=reg95+reg111;
    reg86=reg7*reg86; reg82=reg71*reg72; reg114=reg100+reg114; reg109=reg110+reg109; reg87=reg81*reg12;
    reg91=reg29*reg106; reg48=reg46+reg48; reg46=reg106*reg13; reg99=reg52+reg99; reg52=reg78*reg12;
    reg95=reg63+reg62; reg50=reg79+reg50; reg60=reg3-reg60; reg36=reg15*reg36; reg107=reg57+reg107;
    reg3=reg27*reg5; reg57=reg56*reg88; reg79=reg4*reg28; reg98=reg70*reg2; reg100=reg61*reg9;
    reg110=2*reg26; T reg118=reg96*reg2; T reg119=reg11*reg88; T reg120=reg10*reg21; reg117=reg90+reg117;
    reg90=reg27*reg51; reg54=reg55+reg54; reg55=reg24*reg9; reg102=reg94+reg102; reg94=reg27*reg72;
    reg115=reg41-reg115; reg41=reg88*reg51; reg93=reg6*reg93; T reg121=reg88*reg72; reg108=reg92+reg108;
    reg74=reg14-reg74; reg104=reg103+reg104; reg85=reg53+reg85; reg76=reg49-reg76; reg105=reg116+reg105;
    reg14=reg71*reg51; reg113=reg101+reg113; reg49=reg54*reg0; reg53=reg79*reg38; reg92=reg100*reg42;
    reg101=reg68*reg43; reg52=reg99+reg52; reg99=reg66*reg115; reg103=reg26*reg9; reg116=reg67*reg60;
    T reg122=reg70*reg98; reg82=reg112+reg82; reg112=reg67*reg97; T reg123=reg25*reg111; T reg124=reg25*reg77;
    T reg125=reg67*reg50; T reg126=reg67*reg75; T reg127=reg25*reg107; T reg128=reg66*reg50; T reg129=reg37*reg77;
    T reg130=reg66*reg97; T reg131=reg37*reg111; T reg132=reg66*reg75; T reg133=reg37*reg107; T reg134=reg120*reg97;
    T reg135=reg96*reg111; T reg136=reg120*reg75; T reg137=reg96*reg107; reg87=reg48+reg87; reg85=reg119-reg85;
    reg48=reg106*reg12; reg3=reg57+reg3; reg57=reg55*reg0; reg18=reg18*reg89; reg119=reg110*reg9;
    T reg138=reg23*reg4; T reg139=reg25*reg60; T reg140=reg37*reg115; T reg141=reg100*reg19; T reg142=reg79*reg17;
    reg90=reg117+reg90; reg91=reg114+reg91; reg114=reg69*reg118; reg94=reg102+reg94; reg102=reg69*reg98;
    reg117=reg70*reg118; reg95=reg6*reg95; reg14=reg104+reg14; reg76=reg76/reg40; reg36=reg86-reg36;
    reg93=reg74-reg93; reg41=reg105+reg41; reg74=reg64*reg118; reg46=reg109+reg46; reg86=reg10*reg68;
    reg104=reg21*reg43; reg121=reg108+reg121; reg105=reg64*reg98; reg108=2*reg68; reg34=reg33+reg34;
    reg33=reg80*reg113; reg131=reg130+reg131; reg109=reg79*reg57; reg102=reg94+reg102; reg94=reg80*reg52;
    reg129=reg128+reg129; reg128=reg100*reg49; reg130=reg54*reg57; T reg143=reg78*reg87; reg127=reg126+reg127;
    reg117=reg14+reg117; reg99=reg116+reg99; reg14=reg25*reg46; reg116=reg88*reg85; reg126=reg79*reg49;
    reg74=reg41+reg74; reg41=reg67*reg91; T reg144=reg78*reg113; reg123=reg112+reg123; reg112=reg103*reg1;
    reg105=reg121+reg105; reg121=reg101*reg87; T reg145=reg34*reg1; T reg146=reg120*reg2; T reg147=reg54*reg49;
    reg140=reg139+reg140; reg139=reg27*reg85; reg137=reg136+reg137; reg122=reg82+reg122; reg82=reg84*reg115;
    reg136=reg101*reg113; reg135=reg134+reg135; reg134=reg83*reg60; reg114=reg90+reg114; reg90=reg96*reg77;
    T reg148=reg120*reg50; T reg149=reg80*reg87; reg133=reg132+reg133; reg132=reg100*reg57; T reg150=reg37*reg46;
    T reg151=reg66*reg91; reg19=reg119*reg19; reg42=reg119*reg42; reg38=reg138*reg38; T reg152=reg10*reg4;
    T reg153=reg70*reg76; T reg154=reg83*reg76; T reg155=reg84*reg76; T reg156=reg18*reg13; reg48=reg3+reg48;
    reg3=reg79*reg56; T reg157=reg24*reg68; T reg158=reg100*reg5; T reg159=reg26*reg43; reg92=reg53+reg92;
    reg141=reg142+reg141; reg93=reg93/reg40; reg53=reg18*reg29; reg86=reg104+reg86; reg124=reg125+reg124;
    reg17=reg138*reg17; reg95=reg36-reg95; reg36=reg21*reg68; reg108=reg108*reg89; reg104=reg10*reg43;
    reg125=reg78*reg52; reg142=reg120*reg91; reg53=reg141+reg53; reg141=reg26*reg68; T reg160=reg36*reg93;
    reg150=reg151+reg150; reg151=reg104*reg93; T reg161=reg70*reg155; reg29=reg108*reg29; reg19=reg17+reg19;
    reg33=reg131+reg33; reg17=reg30*reg155; reg131=reg138*reg112; reg109=reg74+reg109; reg74=reg34*reg145;
    T reg162=reg4*reg89; reg147=reg122+reg147; reg90=reg148+reg90; reg126=reg105+reg126; reg105=reg138*reg145;
    reg122=reg101*reg52; reg139=reg140+reg139; reg140=reg69*reg146; reg13=reg108*reg13; reg42=reg38+reg42;
    reg38=reg24*reg43; reg148=reg30*reg153; reg149=reg133+reg149; reg156=reg92+reg156; reg92=reg71*reg76;
    reg5=reg119*reg5; reg56=reg138*reg56; reg133=reg18*reg12; reg136=reg135+reg136; reg135=reg80*reg48;
    reg158=reg3+reg158; reg3=reg32*reg153; reg143=reg127+reg143; reg121=reg137+reg121; reg127=reg21*reg4;
    reg137=reg70*reg153; reg157=reg159+reg157; reg23=reg23*reg68; reg159=reg78*reg48; reg14=reg41+reg14;
    reg41=reg71*reg85; reg82=reg134+reg82; reg28=reg43*reg28; reg134=reg64*reg146; reg116=reg99+reg116;
    reg99=reg32*reg155; reg144=reg123+reg144; reg125=reg124+reg125; reg123=reg32*reg154; reg124=elem.pos(2)[0]-elem.pos(0)[0];
    T reg163=elem.pos(2)[1]-elem.pos(0)[1]; reg128=reg102+reg128; reg102=reg152*reg0; reg40=reg95/reg40; reg95=reg119*reg112;
    reg132=reg114+reg132; reg114=reg119*reg145; T reg164=reg86*reg93; reg94=reg129+reg94; reg129=reg30*reg154;
    T reg165=reg34*reg112; T reg166=reg96*reg46; reg130=reg117+reg130; reg117=elem.pos(1)[0]-elem.pos(0)[0]; T reg167=elem.pos(1)[1]-elem.pos(0)[1];
    reg12=reg108*reg12; reg5=reg56+reg5; reg166=reg142+reg166; reg56=reg101*reg48; reg142=reg117*reg163;
    reg134=reg116+reg134; reg116=(*f.m).alpha_1*reg67; T reg168=reg86*reg160; reg161=reg136+reg161; reg136=reg79*reg102;
    reg122=reg90+reg122; reg90=reg70*reg154; reg74=reg147+reg74; reg147=reg167*reg124; T reg169=reg25*(*f.m).alpha_2;
    T reg170=reg127*reg1; T reg171=reg86*reg164; reg137=reg121+reg137; reg121=(*f.m).alpha_1*reg66; T reg172=reg9*reg89;
    T reg173=reg37*(*f.m).alpha_2; reg133=reg158+reg133; reg158=reg162*reg93; T reg174=reg23*reg151; reg129=reg94+reg129;
    reg94=reg38*reg40; reg95=reg132+reg95; reg114=reg128+reg114; reg128=reg25*reg156; reg132=reg67*reg53;
    T reg175=reg28*reg164; reg3=reg143+reg3; reg143=reg141*reg40; T reg176=reg157*reg40; T reg177=reg70*reg146;
    reg41=reg82+reg41; reg110=reg110*reg68; reg82=reg32*reg92; reg159=reg14+reg159; reg61=reg61*reg43;
    reg14=reg28*reg160; reg99=reg144+reg99; reg123=reg125+reg123; reg125=reg28*reg151; reg131=reg109+reg131;
    reg10=reg10*reg89; reg109=reg4*reg43; reg13=reg42+reg13; reg105=reg126+reg105; reg42=reg37*reg156;
    reg126=reg66*reg53; reg140=reg139+reg140; reg139=reg23*reg164; reg148=reg149+reg148; reg144=reg100*reg102;
    reg165=reg130+reg165; reg130=reg30*reg92; reg17=reg33+reg17; reg33=reg23*reg160; reg29=reg19+reg29;
    reg135=reg150+reg135; reg24=reg24*reg89; reg19=reg9*reg43; reg149=reg78*(*f.m).alpha_3; reg150=(*f.m).alpha_1*reg88;
    reg147=reg142-reg147; reg142=reg80*(*f.m).alpha_3; T reg178=reg105*reg165; reg169=reg116+reg169; reg116=reg95*reg74;
    T reg179=reg131*reg74; T reg180=reg114*reg165; reg173=reg121+reg173; reg139=reg148+reg139; reg12=reg5+reg12;
    reg144=reg140+reg144; reg5=reg119*reg170; reg121=reg23*reg158; reg130=reg135+reg130; reg135=reg54*reg76;
    reg140=reg110*reg143; reg33=reg17+reg33; reg17=reg110*reg94; reg174=reg129+reg174; reg129=reg25*reg13;
    reg148=reg67*reg29; T reg181=reg78*reg133; reg128=reg132+reg128; reg132=reg61*reg176; reg175=reg3+reg175;
    reg3=reg172*reg40; T reg182=reg54*reg102; reg177=reg41+reg177; reg41=reg28*reg158; reg82=reg159+reg82;
    reg159=reg61*reg143; reg14=reg99+reg14; reg99=reg61*reg94; reg125=reg123+reg125; reg123=reg70*reg92;
    T reg183=reg138*reg170; reg56=reg166+reg56; reg166=reg157*reg143; reg171=reg137+reg171; reg168=reg161+reg168;
    reg137=reg157*reg176; reg136=reg134+reg136; reg134=reg120*reg53; reg161=reg96*reg156; reg21=reg21*reg89;
    T reg184=reg86*reg151; T reg185=reg27*(*f.m).alpha_2; reg90=reg122+reg90; reg4=reg4*reg68; reg10=reg109+reg10;
    reg109=reg37*reg13; reg122=reg66*reg29; T reg186=reg80*reg133; reg42=reg126+reg42; reg126=reg110*reg176;
    reg163=reg163/reg147; reg117=reg117/reg147; reg24=reg19+reg24; reg183=reg136+reg183; reg182=reg177+reg182;
    reg19=reg34*reg170; reg21=reg4+reg21; reg167=reg167/reg147; reg124=reg124/reg147; reg4=vectors[0][indices[1]+1]-vectors[0][indices[0]+1];
    reg136=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg177=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; T reg187=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; reg5=reg144+reg5; reg144=reg106*(*f.m).alpha_3;
    T reg188=(*f.m).alpha_1*reg120; T reg189=reg105*reg95; T reg190=reg101*reg133; reg161=reg134+reg161; reg137=reg171+reg137;
    reg134=(*f.m).alpha_2*reg96; reg171=reg86*reg158; reg123=reg56+reg123; reg166=reg168+reg166; reg56=reg131*reg114;
    reg168=reg157*reg94; reg184=reg90+reg184; reg178=reg179-reg178; reg90=reg10*reg93; reg17=reg174+reg17;
    reg76=reg34*reg76; reg140=reg33+reg140; reg121=reg130+reg121; reg33=reg110*reg3; reg130=reg80*reg12;
    reg109=reg122+reg109; reg126=reg139+reg126; reg186=reg42+reg186; reg42=reg30*reg135; reg122=reg9*reg68;
    reg26=reg26*reg89; reg99=reg125+reg99; reg159=reg14+reg159; reg41=reg82+reg41; reg14=reg61*reg3;
    reg132=reg175+reg132; reg180=reg116-reg180; reg181=reg128+reg181; reg82=reg32*reg135; reg149=reg169+reg149;
    reg185=reg150+reg185; reg120=reg120*reg29; reg96=reg96*reg13; reg116=reg78*reg12; reg142=reg173+reg142;
    reg129=reg148+reg129; reg125=reg17*reg149; reg128=reg140*reg142; reg116=reg129+reg116; reg33=reg121+reg33;
    reg32=reg32*reg76; reg14=reg41+reg14; reg41=reg24*reg40; reg189=reg56-reg189; reg93=reg21*reg93;
    reg19=reg182+reg19; reg56=reg124*reg187; reg82=reg181+reg82; reg121=reg28*reg90; reg129=reg167*reg136;
    reg139=reg163*reg4; reg148=reg117*reg177; reg150=reg99*reg149; reg169=reg159*reg142; reg173=reg101*reg12;
    reg96=reg120+reg96; reg144=reg185+reg144; reg120=reg70*reg135; reg190=reg161+reg190; reg161=reg157*reg3;
    reg171=reg123+reg171; reg168=reg184+reg168; reg134=reg188+reg134; reg123=(*f.m).alpha_3*reg101; reg174=reg5*reg178;
    reg175=(*f.m).alpha_1*reg152; reg179=reg132*reg166; reg26=reg122+reg26; reg42=reg186+reg42; reg122=reg23*reg90;
    reg43=reg89*reg43; reg181=reg140*reg137; reg130=reg109+reg130; reg30=reg30*reg76; reg109=reg180*reg183;
    reg182=reg159*reg137; reg184=reg126*reg166; reg185=(*f.m).alpha_2*reg55; reg186=reg159*reg126; reg188=reg114*reg19;
    T reg191=reg5*reg74; T reg192=reg105*reg19; T reg193=reg132*reg140; T reg194=reg183*reg74; reg179=reg182-reg179;
    reg184=reg181-reg184; reg129=reg139-reg129; reg14=reg14*reg144; reg169=reg150+reg169; reg139=(*f.m).alpha_2*reg103;
    reg150=(*f.m).alpha_1*reg127; reg181=(*f.m).alpha_3*reg43; reg185=reg175+reg185; reg128=reg125+reg128; reg33=reg33*reg144;
    reg125=reg168*reg149; reg175=reg166*reg142; reg187=reg163*reg187; reg177=reg167*reg177; reg123=reg134+reg123;
    reg136=reg117*reg136; reg174=reg109-reg174; reg4=reg124*reg4; reg109=reg19*reg189; reg56=reg148-reg56;
    reg122=reg42+reg122; reg42=reg110*reg41; reg134=reg28*reg93; reg32=reg116+reg32; reg30=reg130+reg30;
    reg116=reg61*reg41; reg121=reg82+reg121; reg82=reg23*reg93; reg40=reg26*reg40; reg161=reg171+reg161;
    reg120=reg190+reg120; reg130=reg86*reg90; reg68=reg89*reg68; reg173=reg96+reg173; reg96=reg70*reg76;
    reg130=reg120+reg130; reg120=reg157*reg41; reg148=reg110*reg40; reg82=reg30+reg82; reg96=reg173+reg96;
    reg30=reg86*reg93; reg175=reg125+reg175; reg161=reg161*reg144; reg42=reg122+reg42; reg181=reg185+reg181;
    reg4=reg136-reg4; elem.epsilon[0][1]=reg4; reg177=reg187-reg177; elem.epsilon[0][0]=reg177; reg33=reg128+reg33;
    reg129=reg56+reg129; reg56=reg132*reg123; reg14=reg169+reg14; reg122=reg126*reg123; reg193=reg186-reg193;
    reg125=reg17*reg179; reg128=reg99*reg184; reg136=reg105*reg5; reg169=reg183*reg114; reg171=reg131*reg19;
    reg192=reg194-reg192; reg173=reg183*reg165; reg182=reg95*reg19; reg188=reg191-reg188; reg185=reg5*reg165;
    reg116=reg121+reg116; reg134=reg32+reg134; reg32=reg61*reg40; reg109=reg174+reg109; reg121=(*f.m).alpha_3*reg68;
    reg139=reg150+reg139; reg150=reg80*reg59; reg174=reg78*reg60; reg186=reg80*reg115; reg187=(*f.m).deltaT*reg149;
    reg190=(*f.m).deltaT*reg142; reg191=reg78*reg65; reg161=reg175+reg161; reg175=reg137*reg123; reg194=reg183*reg95;
    T reg195=reg80*reg47; reg136=reg169-reg136; reg169=reg131*reg5; T reg196=reg157*reg40; reg30=reg96+reg30;
    reg96=reg78*reg73; reg171=reg173-reg171; reg120=reg130+reg120; reg192=reg192/reg109; reg180=reg180/reg109;
    reg178=reg178/reg109; reg130=reg106*reg51; reg182=reg185-reg182; reg148=reg82+reg148; reg150=reg191+reg150;
    reg186=reg174+reg186; reg82=reg106*reg85; reg173=(*f.m).deltaT*reg123; reg121=reg139+reg121; reg188=reg188/reg109;
    reg32=reg134+reg32; reg134=reg132*reg168; reg139=reg4-reg190; reg125=reg128-reg125; reg128=reg168*reg193;
    reg116=reg116*reg181; reg56=reg14+reg56; reg14=reg177-reg187; reg174=reg17*reg137; reg185=reg126*reg168;
    reg122=reg33+reg122; reg129=0.5*reg129; elem.epsilon[0][2]=reg129; reg42=reg42*reg181; reg33=reg99*reg137;
    reg148=reg148*reg121; reg191=reg99*reg126; reg182=reg182/reg109; T reg197=reg132*reg17; reg42=reg122+reg42;
    reg116=reg56+reg116; reg32=reg32*reg121; reg56=reg129-reg173; reg122=reg159*reg168; reg134=reg33-reg134;
    reg33=reg99*reg166; T reg198=reg140*reg168; reg171=reg171/reg109; reg185=reg174-reg185; reg174=reg188*reg14;
    reg189=reg189/reg109; T reg199=reg17*reg166; T reg200=reg180*reg14; T reg201=reg192*reg139; reg136=reg136/reg109;
    reg128=reg125+reg128; reg125=reg178*reg139; reg169=reg194-reg169; reg194=reg106*reg72; T reg202=reg81*reg118;
    reg130=reg150+reg130; reg82=reg186+reg82; reg150=reg81*reg146; reg196=reg30+reg196; reg195=reg96+reg195;
    reg175=reg161+reg175; reg120=reg120*reg181; reg125=reg200-reg125; reg120=reg175+reg120; reg196=reg196*reg121;
    reg14=reg182*reg14; reg139=reg171*reg139; reg184=reg184/reg128; reg185=reg185/reg128; reg109=reg169/reg109;
    reg194=reg195+reg194; reg198=reg199-reg198; reg179=reg179/reg128; reg148=reg42+reg148; reg134=reg134/reg128;
    reg122=reg33-reg122; reg30=reg99*reg140; reg197=reg191-reg197; reg33=reg159*reg17; reg42=reg189*reg56;
    reg32=reg116+reg32; reg96=reg81*reg98; reg202=reg130+reg202; reg116=reg136*reg56; reg174=reg201-reg174;
    reg130=reg18*reg57; reg161=reg18*reg102; reg150=reg82+reg150; reg198=reg198/reg128; reg116=reg174-reg116;
    reg184=reg184*reg32; reg179=reg179*reg148; reg122=reg122/reg128; reg193=reg193/reg128; reg197=reg197/reg128;
    reg185=reg185*reg32; reg134=reg134*reg148; reg33=reg30-reg33; reg161=reg150+reg161; reg30=reg18*reg49;
    reg96=reg194+reg96; reg82=reg108*reg170; reg125=reg42+reg125; reg56=reg109*reg56; reg196=reg120+reg196;
    reg42=reg108*reg112; reg139=reg14-reg139; reg130=reg202+reg130; reg42=reg130+reg42; reg197=reg197*reg196;
    reg128=reg33/reg128; reg185=reg134-reg185; reg32=reg198*reg32; reg14=reg5*reg125; reg33=reg108*reg145;
    reg30=reg96+reg30; reg96=reg183*reg125; reg82=reg161+reg82; reg120=reg95*reg116; reg179=reg184-reg179;
    reg193=reg193*reg196; reg130=reg131*reg116; reg148=reg122*reg148; reg139=reg56+reg139; reg56=PNODE(1).dep[0]-PNODE(0).dep[0];
    reg122=PNODE(2).dep[0]-PNODE(0).dep[0]; reg134=PNODE(1).dep[1]-PNODE(0).dep[1]; reg150=PNODE(2).dep[1]-PNODE(0).dep[1]; reg197=reg185-reg197; reg161=1-(*f.m).resolution;
    reg196=reg128*reg196; reg179=reg193+reg179; reg148=reg32-reg148; reg32=reg42*reg116; reg128=reg82*reg125;
    reg33=reg30+reg33; reg30=reg104*reg60; reg169=reg36*reg115; reg130=reg96+reg130; reg96=reg114*reg139;
    reg174=reg104*reg65; reg175=reg105*reg139; reg120=reg14+reg120; reg14=reg36*reg59; reg184=reg134*reg163;
    reg185=reg36*reg47; reg186=reg104*reg73; reg14=reg174+reg14; reg32=reg128+reg32; reg65=reg38*reg65;
    reg59=reg141*reg59; reg60=reg38*reg60; reg115=reg141*reg115; reg148=reg196+reg148; reg128=reg162*reg51;
    reg149=(*f.m).resolution*reg149; reg175=reg130+reg175; reg130=reg167*reg150; reg174=reg124*reg56; reg96=reg120+reg96;
    reg169=reg30+reg169; reg179=reg161*reg179; reg30=reg162*reg85; reg120=reg122*reg117; reg191=reg33*reg139;
    reg197=reg161*reg197; reg142=(*f.m).resolution*reg142; reg175=reg187+reg175; reg148=reg161*reg148; reg193=reg162*reg72;
    reg185=reg186+reg185; reg179=reg149+reg179; reg197=reg142+reg197; reg144=(*f.m).deltaT*reg144; reg30=reg169+reg30;
    reg142=reg86*reg146; reg191=reg32+reg191; reg128=reg14+reg128; reg14=reg86*reg118; reg47=reg141*reg47;
    reg73=reg38*reg73; reg150=reg150*reg117; reg122=reg167*reg122; reg56=reg56*reg163; reg123=(*f.m).resolution*reg123;
    reg96=reg190+reg96; reg130=reg184-reg130; reg134=reg124*reg134; reg51=reg172*reg51; reg59=reg65+reg59;
    reg85=reg172*reg85; reg115=reg60+reg115; reg174=reg120-reg174; reg47=reg73+reg47; reg72=reg172*reg72;
    reg134=reg150-reg134; reg32=reg10*reg57; reg14=reg128+reg14; reg122=reg56-reg122; reg130=reg174+reg130;
    reg56=reg175+reg96; reg60=(*f.m).resolution*reg180; reg65=(*f.m).resolution*reg188; reg191=reg144+reg191; reg73=(*f.m).resolution*reg192;
    reg120=(*f.m).resolution*reg178; reg128=reg10*reg102; reg142=reg30+reg142; reg140=reg140*reg161; reg17=reg17*reg161;
    reg159=reg159*reg161; reg146=reg157*reg146; reg85=reg115+reg85; reg193=reg185+reg193; reg51=reg59+reg51;
    reg118=reg157*reg118; reg30=reg86*reg98; reg148=reg123+reg148; reg179=(*f.m).deltaT*reg179; reg197=(*f.m).deltaT*reg197;
    reg168=reg168*reg161; reg99=reg99*reg161; reg166=reg166*reg161; reg59=(*f.m).resolution*reg171; reg115=(*f.m).resolution*reg182;
    reg56=reg191+reg56; reg98=reg157*reg98; reg72=reg47+reg72; reg120=reg159-reg120; reg99=reg60+reg99;
    reg137=reg137*reg161; reg126=reg126*reg161; reg161=reg132*reg161; reg47=(*f.m).resolution*reg109; reg60=(*f.m).resolution*reg136;
    reg130=0.5*reg130; reg65=reg17-reg65; reg140=reg73+reg140; reg17=(*f.m).resolution*reg189; reg73=reg134-reg197;
    reg102=reg24*reg102; reg146=reg85+reg146; reg85=reg122-reg179; reg168=reg115+reg168; reg115=reg10*reg49;
    reg30=reg193+reg30; reg128=reg142+reg128; reg123=reg21*reg170; reg59=reg166-reg59; reg118=reg51+reg118;
    reg57=reg24*reg57; reg32=reg14+reg32; reg14=reg21*reg112; reg148=(*f.m).deltaT*reg148; reg60=reg126-reg60;
    reg161=reg17+reg161; reg17=reg165*reg116; reg51=reg19*reg125; reg56=reg56/3; reg170=reg26*reg170;
    reg102=reg146+reg102; reg126=reg21*reg145; reg115=reg30+reg115; reg123=reg128+reg123; reg14=reg32+reg14;
    reg112=reg26*reg112; reg57=reg118+reg57; reg98=reg72+reg98; reg49=reg24*reg49; reg30=reg59*reg73;
    reg32=reg168*reg85; reg72=reg140*reg73; reg118=reg65*reg85; reg137=reg47+reg137; reg47=reg130-reg148;
    reg128=reg120*reg73; reg132=reg99*reg85; reg112=reg57+reg112; reg49=reg98+reg49; reg145=reg26*reg145;
    reg126=reg115+reg126; reg170=reg102+reg170; reg57=reg137*reg47; reg175=reg175-reg56; reg96=reg96-reg56;
    reg30=reg32+reg30; reg17=reg51+reg17; reg32=reg74*reg139; reg51=reg123*reg125; reg128=reg132+reg128;
    reg98=reg161*reg47; reg72=reg118+reg72; reg102=reg60*reg47; reg115=reg14*reg116; reg57=reg30+reg57;
    reg115=reg51+reg115; reg175=pow(reg175,2); reg102=reg72+reg102; reg96=pow(reg96,2); reg56=reg191-reg56;
    reg30=reg126*reg139; reg125=reg170*reg125; reg116=reg112*reg116; reg145=reg49+reg145; reg32=reg17+reg32;
    reg98=reg128+reg98; reg30=reg115+reg30; reg116=reg125+reg116; reg197=reg4-reg197; reg57=2*reg57;
    reg73=reg73*reg102; reg181=(*f.m).deltaT*reg181; reg85=reg85*reg98; reg139=reg145*reg139; reg32=reg173+reg32;
    reg96=reg175+reg96; reg56=pow(reg56,2); reg179=reg177-reg179; reg65=reg65*reg179; reg121=(*f.m).deltaT*reg121;
    reg32=pow(reg32,2); reg148=reg129-reg148; reg139=reg116+reg139; reg120=reg120*reg197; reg99=reg99*reg179;
    reg56=reg96+reg56; reg30=reg181+reg30; reg140=reg140*reg197; reg47=reg47*reg57; reg73=reg85+reg73;
    reg30=pow(reg30,2); reg139=reg121+reg139; reg32=reg56+reg32; reg47=reg73+reg47; reg120=reg99+reg120;
    reg161=reg161*reg148; reg140=reg65+reg140; reg60=reg60*reg148; reg179=reg168*reg179; reg197=reg59*reg197;
    reg47=reg147*reg47; reg30=reg32+reg30; reg139=pow(reg139,2); reg161=reg120+reg161; elem.sigma[0][0]=reg161;
    reg60=reg140+reg60; elem.sigma[0][1]=reg60; reg197=reg179+reg197; reg148=reg137*reg148; reg148=reg197+reg148;
    elem.sigma[0][2]=reg148; reg4=reg67*reg161; reg17=reg66*reg60; reg139=reg30+reg139; reg30=reg25*reg161;
    reg32=reg37*reg60; reg60=reg84*reg60; reg49=0.16666666666666665741*reg47; reg161=reg83*reg161; reg47=0.33333333333333331483*reg47;
    reg51=reg70*reg148; reg60=reg161+reg60; reg47=reg49+reg47; reg139=pow(reg139,0.5); reg49=reg69*reg148;
    reg32=reg30+reg32; reg148=reg64*reg148; reg17=reg4+reg17; elem.sigma_von_mises=0.86602540378443859659*reg139; elem.ener=reg47/2;
    elem.sigma_local[0][0]=reg17+reg148; elem.sigma_local[0][1]=reg32+reg49; elem.sigma_local[0][2]=reg60+reg51;

      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_2(TE &elem,TF &f, TVEVE &vectors,const unsigned *indices) {
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
    T reg0=2*(*f.m).shear_modulus_13; T reg1=2*(*f.m).shear_modulus_23; reg0=1.0/reg0; T reg2=2*(*f.m).shear_modulus_12; reg1=1.0/reg1;
    T reg3=reg0*reg1; reg2=1.0/reg2; T reg4=pow((*f.m).v1[0],2); T reg5=pow((*f.m).v1[1],2); T reg6=pow((*f.m).v2[0],2);
    T reg7=pow((*f.m).v2[1],2); T reg8=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg9=reg2*reg3; T reg10=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg11=1.0/(*f.m).elastic_modulus_3;
    T reg12=1.0/(*f.m).elastic_modulus_2; T reg13=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg14=reg10*reg9; reg7=reg6+reg7; reg6=pow((*f.m).v1[2],2);
    T reg15=reg8*reg9; T reg16=pow((*f.m).v2[2],2); reg5=reg4+reg5; reg4=reg11*reg9; T reg17=reg12*reg4;
    reg5=reg6+reg5; reg6=reg13*reg4; T reg18=reg10*reg14; T reg19=reg10*reg15; reg16=reg7+reg16;
    reg7=1.0/(*f.m).elastic_modulus_1; reg16=pow(reg16,0.5); reg18=reg17-reg18; reg5=pow(reg5,0.5); reg19=reg6+reg19;
    reg17=reg13*reg14; T reg20=reg12*reg15; T reg21=(*f.m).v2[1]/reg16; T reg22=(*f.m).v2[2]/reg16; T reg23=(*f.m).v1[2]/reg5;
    T reg24=(*f.m).v1[1]/reg5; T reg25=reg7*reg18; T reg26=reg13*reg19; T reg27=reg17+reg20; T reg28=reg10*reg3;
    reg26=reg25-reg26; reg25=reg8*reg27; T reg29=reg12*reg9; T reg30=reg8*reg14; reg4=reg7*reg4;
    T reg31=reg8*reg15; reg9=reg13*reg9; T reg32=reg11*reg3; T reg33=reg2*reg1; reg3=reg8*reg3;
    reg5=(*f.m).v1[0]/reg5; T reg34=reg23*reg21; reg16=(*f.m).v2[0]/reg16; T reg35=reg24*reg22; T reg36=reg10*reg28;
    T reg37=reg13*reg32; T reg38=reg11*reg33; T reg39=reg2*reg0; T reg40=reg8*reg33; reg32=reg12*reg32;
    T reg41=reg10*reg3; reg33=reg10*reg33; T reg42=reg5*reg22; T reg43=reg8*reg9; reg14=reg7*reg14;
    reg31=reg4-reg31; reg4=2*reg16; reg15=reg13*reg15; reg30=reg6+reg30; reg6=reg35-reg34;
    T reg44=reg23*reg16; T reg45=2*reg5; reg25=reg26-reg25; reg26=reg8*reg29; reg30=reg30/reg25;
    reg19=reg19/reg25; reg41=reg37+reg41; reg29=reg7*reg29; reg37=2*reg6; T reg46=reg5*reg21;
    T reg47=reg12*reg38; reg38=reg13*reg38; T reg48=reg10*reg33; reg18=reg18/reg25; T reg49=reg4*reg21;
    T reg50=reg10*reg40; T reg51=pow(reg21,2); T reg52=pow(reg16,2); T reg53=reg11*reg39; T reg54=reg24*reg16;
    T reg55=pow(reg24,2); T reg56=pow(reg5,2); reg28=reg13*reg28; reg9=reg13*reg9; reg15=reg15+reg14;
    reg43=reg14+reg43; reg14=reg24*reg45; reg3=reg12*reg3; T reg57=reg10*reg39; reg26=reg17+reg26;
    reg31=reg31/reg25; T reg58=reg44-reg42; reg36=reg32-reg36; reg39=reg8*reg39; reg32=pow(reg23,2);
    reg26=reg26/reg25; T reg59=reg52*reg31; T reg60=reg30*reg56; T reg61=reg30*reg55; T reg62=reg51*reg31;
    T reg63=reg49*reg19; T reg64=reg18*reg14; T reg65=reg51*reg19; T reg66=reg18*reg55; T reg67=reg10*reg39;
    T reg68=reg10*reg57; T reg69=reg13*reg53; reg53=reg12*reg53; T reg70=reg18*reg56; reg50=reg38+reg50;
    reg40=reg12*reg40; reg9=reg29-reg9; reg29=reg28+reg3; reg41=reg13*reg41; reg15=reg15/reg25;
    reg48=reg47-reg48; reg38=reg46-reg54; reg43=reg43/reg25; reg33=reg13*reg33; reg36=reg7*reg36;
    reg47=pow(reg22,2); T reg71=pow(reg6,2); T reg72=pow(reg58,2); T reg73=reg37*reg58; reg27=reg27/reg25;
    T reg74=reg49*reg31; T reg75=reg30*reg14; T reg76=reg52*reg19; T reg77=pow(reg38,2); reg9=reg9/reg25;
    reg57=reg13*reg57; reg70=reg76+reg70; reg76=reg71*reg27; reg39=reg12*reg39; reg41=reg36-reg41;
    reg29=reg8*reg29; reg36=reg49*reg43; reg48=reg7*reg48; reg50=reg13*reg50; T reg78=reg33+reg40;
    reg68=reg53-reg68; reg67=reg69+reg67; reg53=reg47*reg31; reg65=reg66+reg65; reg66=reg72*reg27;
    reg69=reg18*reg32; T reg79=reg47*reg19; reg63=reg64+reg63; reg64=reg73*reg27; T reg80=reg30*reg32;
    T reg81=reg72*reg15; reg62=reg61+reg62; reg61=reg71*reg15; reg59=reg60+reg59; reg74=reg75+reg74;
    reg60=reg73*reg15; reg75=reg26*reg56; T reg82=reg52*reg43; T reg83=reg26*reg55; T reg84=reg51*reg43;
    T reg85=reg26*reg14; T reg86=reg5*reg24; reg78=reg8*reg78; T reg87=2*reg24; T reg88=reg23*reg45;
    T reg89=2*reg21; T reg90=reg47*reg43; reg68=reg7*reg68; T reg91=reg16*reg21; reg53=reg80+reg53;
    reg80=reg26*reg32; reg67=reg13*reg67; reg82=reg75+reg82; reg75=reg72*reg9; T reg92=reg57+reg39;
    T reg93=reg73*reg9; reg76=reg70+reg76; reg36=reg85+reg36; reg66=reg65+reg66; reg79=reg69+reg79;
    reg65=reg27*reg77; reg84=reg83+reg84; reg29=reg41-reg29; reg64=reg63+reg64; reg41=reg71*reg9;
    reg50=reg48-reg50; reg48=reg77*reg15; reg81=reg62+reg81; reg60=reg74+reg60; reg61=reg59+reg61;
    reg59=reg4*reg22; reg62=reg89*reg22; reg41=reg82+reg41; reg63=reg87*reg21; reg29=reg29/reg25;
    reg93=reg36+reg93; reg37=reg37*reg38; reg36=2*reg58; reg69=2*reg23; reg70=reg87*reg23;
    reg78=reg50-reg78; reg50=reg91*reg60; reg67=reg68-reg67; reg68=reg58*reg6; reg92=reg8*reg92;
    reg65=reg79+reg65; reg74=reg59*reg31; reg79=reg88*reg18; reg82=reg59*reg19; reg83=reg88*reg30;
    reg48=reg53+reg48; reg53=reg52*reg81; reg85=reg56*reg64; T reg94=reg52*reg60; T reg95=reg55*reg76;
    T reg96=reg51*reg61; T reg97=reg56*reg66; reg90=reg80+reg90; reg80=reg55*reg66; T reg98=reg52*reg61;
    T reg99=reg51*reg81; T reg100=reg56*reg76; T reg101=reg77*reg9; T reg102=reg55*reg64; reg60=reg51*reg60;
    T reg103=reg24*reg21; reg46=reg54+reg46; reg75=reg84+reg75; reg54=reg5*reg16; reg81=reg91*reg81;
    reg66=reg86*reg66; reg84=reg16*reg45; reg64=reg86*reg64; T reg104=reg24*reg6; T reg105=reg5*reg58;
    T reg106=reg51*reg12; reg53=reg97+reg53; reg97=reg71*reg75; T reg107=reg37*reg15; T reg108=reg52*reg7;
    reg74=reg83+reg74; reg83=reg59*reg43; T reg109=reg52*reg13; T reg110=reg71*reg41; reg98=reg100+reg98;
    reg100=reg63*reg13; T reg111=reg46*reg29; T reg112=reg54*reg29; reg101=reg90+reg101; reg90=reg103*reg29;
    T reg113=reg51*reg13; T reg114=reg84*reg7; T reg115=reg88*reg26; T reg116=reg68*reg75; reg81=reg66+reg81;
    reg61=reg91*reg61; reg76=reg86*reg76; reg66=reg72*reg93; reg60=reg102+reg60; reg102=reg51*reg48;
    T reg117=reg55*reg65; reg75=reg72*reg75; reg99=reg80+reg99; reg50=reg64+reg50; reg64=reg72*reg41;
    reg96=reg95+reg96; reg80=reg68*reg93; reg93=reg71*reg93; reg94=reg85+reg94; reg85=reg63*reg12;
    reg95=reg84*reg13; T reg118=reg52*reg48; reg31=reg62*reg31; T reg119=reg56*reg65; reg30=reg70*reg30;
    T reg120=reg5*reg6; reg36=reg36*reg38; T reg121=reg23*reg22; T reg122=reg24*reg58; reg69=reg69*reg22;
    reg78=reg78/reg25; reg92=reg67-reg92; reg82=reg79+reg82; reg67=reg37*reg27; reg79=reg16*reg58;
    reg18=reg70*reg18; reg19=reg62*reg19; T reg123=reg21*reg6; reg105=reg104+reg105; reg65=reg86*reg65;
    reg48=reg91*reg48; reg75=reg99+reg75; reg99=reg63*reg90; reg104=reg16*reg6; reg43=reg62*reg43;
    reg26=reg70*reg26; T reg124=reg37*reg9; reg83=reg115+reg83; reg80=reg50+reg80; reg50=reg46*reg111;
    reg115=reg21*reg58; reg66=reg60+reg66; reg60=reg84*reg112; reg110=reg98+reg110; reg98=reg63*reg111;
    T reg125=reg72*reg101; reg102=reg117+reg102; reg117=reg23*reg38; reg61=reg76+reg61; reg41=reg68*reg41;
    reg97=reg53+reg97; reg53=reg84*reg90; reg76=reg105*reg78; T reg126=reg122*reg78; T reg127=reg120*reg78;
    reg116=reg81+reg116; reg90=reg46*reg90; reg81=reg121*reg29; reg111=reg84*reg111; reg100=reg114-reg100;
    reg114=reg8*reg69; reg113=reg108-reg113; reg108=reg47*reg8; reg25=reg92/reg25; reg92=reg47*reg10;
    reg109=reg106-reg109; reg107=reg74+reg107; reg67=reg82+reg67; reg74=reg69*reg10; reg31=reg30+reg31;
    reg15=reg36*reg15; reg19=reg18+reg19; reg95=reg85-reg95; reg79=reg123+reg79; reg93=reg94+reg93;
    reg118=reg119+reg118; reg18=reg71*reg101; reg30=reg63*reg10; reg82=reg8*reg84; reg27=reg36*reg27;
    reg45=reg6*reg45; reg85=reg63*reg112; reg87=reg87*reg58; reg94=reg51*reg10; reg106=reg52*reg8;
    reg64=reg96+reg64; reg96=reg63*reg81; reg125=reg102+reg125; reg102=reg5*reg38; reg111=reg93+reg111;
    reg93=reg45*reg76; reg119=reg87*reg126; reg123=reg56*reg67; T reg128=reg52*reg107; reg99=reg75+reg99;
    reg85=reg64+reg85; reg64=reg87*reg127; reg75=reg22*reg38; reg12=reg12*reg55; T reg129=reg13*reg56;
    reg114=reg100-reg114; reg108=reg113-reg108; reg13=reg13*reg55; reg92=reg109-reg92; reg74=reg95-reg74;
    reg95=reg47*reg11; reg94=reg106+reg94; reg69=reg69*reg11; reg30=reg82+reg30; reg82=(*f.m).alpha_1*reg56;
    reg100=reg105*reg76; reg50=reg80+reg50; reg101=reg68*reg101; reg48=reg65+reg48; reg65=reg105*reg126;
    reg90=reg116+reg90; reg44=reg42+reg44; reg112=reg46*reg112; reg41=reg61+reg41; reg42=reg51*reg107;
    reg61=reg55*reg67; reg76=reg87*reg76; reg98=reg66+reg98; reg66=reg23*reg6; reg18=reg118+reg18;
    reg80=reg84*reg81; reg4=reg4*reg6; reg89=reg89*reg58; reg126=reg45*reg126; reg53=reg97+reg53;
    reg97=reg79*reg25; reg7=reg7*reg56; reg106=reg51*(*f.m).alpha_2; reg109=(*f.m).alpha_1*reg55; reg113=reg115*reg25;
    reg116=reg104*reg25; reg27=reg19+reg27; reg15=reg31+reg15; reg19=reg117*reg78; reg124=reg83+reg124;
    reg60=reg110+reg60; reg31=reg45*reg127; reg9=reg36*reg9; reg43=reg26+reg43; reg26=reg52*(*f.m).alpha_2;
    reg102=reg66+reg102; reg66=reg51*reg74; reg83=reg16*reg38; reg30=reg69-reg30; reg69=reg54*reg108;
    reg110=reg103*reg92; reg118=reg4*reg116; T reg130=reg51*reg15; T reg131=reg55*reg27; reg94=reg95-reg94;
    reg31=reg60+reg31; reg60=reg72*reg124; reg42=reg61+reg42; reg61=reg10*reg55; reg95=reg89*reg97;
    reg76=reg98+reg76; reg81=reg46*reg81; reg101=reg48+reg101; reg9=reg43+reg9; reg43=reg79*reg113;
    reg65=reg90+reg65; reg48=reg44*reg29; reg35=reg34+reg35; reg127=reg105*reg127; reg112=reg41+reg112;
    reg34=reg75*reg25; reg41=reg47*(*f.m).alpha_2; reg90=(*f.m).alpha_1*reg32; reg100=reg50+reg100; reg50=reg56*reg108;
    reg98=reg55*reg92; T reg132=reg79*reg97; reg67=reg86*reg67; reg107=reg91*reg107; T reg133=reg56*reg114;
    T reg134=reg55*reg74; T reg135=reg22*reg6; T reg136=reg52*reg108; T reg137=reg51*reg92; T reg138=reg52*reg114;
    T reg139=reg89*reg116; T reg140=reg103*reg74; T reg141=reg54*reg114; reg129=reg12-reg129; reg64=reg85+reg64;
    reg13=reg7-reg13; reg7=reg8*reg32; reg106=reg109+reg106; reg10=reg10*reg32; reg12=reg4*reg113;
    reg80=reg18+reg80; reg126=reg53+reg126; reg18=reg52*reg15; reg53=reg56*reg27; reg26=reg82+reg26;
    reg82=reg71*(*f.m).alpha_3; reg128=reg123+reg128; reg97=reg4*reg97; reg85=reg87*reg19; reg96=reg125+reg96;
    reg109=reg71*reg124; reg8=reg8*reg56; reg123=reg45*reg19; reg125=reg72*(*f.m).alpha_3; reg93=reg111+reg93;
    reg111=reg23*reg58; T reg142=reg24*reg38; reg113=reg89*reg113; reg119=reg99+reg119; reg134=reg133+reg134;
    reg97=reg93+reg97; reg124=reg68*reg124; reg107=reg67+reg107; reg67=reg32*reg30; reg10=reg129-reg10;
    reg29=reg35*reg29; reg93=reg102*reg78; reg123=reg80+reg123; reg116=reg79*reg116; reg127=reg112+reg127;
    reg109=reg128+reg109; reg80=reg84*reg48; reg82=reg26+reg82; reg26=reg121*reg30; reg18=reg53+reg18;
    reg125=reg106+reg125; reg53=reg71*reg9; reg7=reg13-reg7; reg13=reg77*(*f.m).alpha_3; reg140=reg141+reg140;
    reg98=reg50+reg98; reg50=reg32*reg94; reg132=reg100+reg132; reg99=reg72*reg9; reg130=reg131+reg130;
    reg142=reg111+reg142; reg100=reg121*reg94; reg19=reg105*reg19; reg5=reg5*reg23; reg61=reg8+reg61;
    reg8=reg63*reg48; reg60=reg42+reg60; reg113=reg119+reg113; reg42=(*f.m).alpha_1*reg86; reg95=reg76+reg95;
    reg76=(*f.m).alpha_2*reg91; reg11=reg11*reg32; reg16=reg16*reg22; reg118=reg31+reg118; reg85=reg96+reg85;
    reg31=reg89*reg34; reg43=reg65+reg43; reg41=reg90+reg41; reg27=reg86*reg27; reg137=reg136+reg137;
    reg65=reg47*reg94; reg15=reg91*reg15; reg90=reg4*reg34; reg96=reg46*reg2; reg66=reg138+reg66;
    reg106=reg47*reg30; reg91=reg91*reg2; reg83=reg135+reg83; reg111=reg22*reg58; reg112=reg21*reg38;
    reg139=reg64+reg139; reg110=reg69+reg110; reg81=reg101+reg81; reg12=reg126+reg12; reg64=reg56*reg7;
    reg69=reg118*reg82; reg101=reg12*reg125; reg119=reg113*reg125; reg124=reg107+reg124; reg48=reg46*reg48;
    reg107=reg44*reg0; reg126=reg16*reg0; reg128=reg139*reg82; reg15=reg27+reg15; reg9=reg68*reg9;
    reg22=reg21*reg22; reg23=reg24*reg23; reg112=reg111+reg112; reg6=reg38*reg6; reg61=reg11-reg61;
    reg11=reg49*reg91; reg106=reg66+reg106; reg21=reg49*reg96; reg63=reg63*reg29; reg99=reg130+reg99;
    reg100=reg110+reg100; reg24=reg46*reg91; reg27=reg87*reg93; reg8=reg60+reg8; reg76=reg42+reg76;
    reg68=(*f.m).alpha_3*reg68; reg31=reg85+reg31; reg42=(*f.m).alpha_1*reg5; reg16=(*f.m).alpha_2*reg16; reg60=reg97*reg43;
    reg84=reg84*reg29; reg53=reg18+reg53; reg26=reg140+reg26; reg18=reg46*reg96; reg66=reg45*reg93;
    reg80=reg109+reg80; reg90=reg123+reg90; reg85=reg113*reg132; reg109=reg95*reg43; reg110=reg12*reg132;
    reg50=reg98+reg50; reg98=reg14*reg91; reg111=reg55*reg10; reg67=reg134+reg67; reg123=reg14*reg96;
    reg129=reg83*reg25; reg130=reg52*reg7; reg116=reg127+reg116; reg78=reg142*reg78; reg127=reg51*reg10;
    reg65=reg137+reg65; reg13=reg41+reg13; reg19=reg81+reg19; reg34=reg79*reg34; reg60=reg110-reg60;
    reg32=reg32*reg61; reg109=reg85-reg109; reg24=reg100+reg24; reg41=reg44*reg126; reg81=reg43*reg125;
    reg111=reg64+reg111; reg119=reg128+reg119; reg18=reg26+reg18; reg26=reg44*reg107; reg64=(*f.m).alpha_2*reg22;
    reg31=reg31*reg13; reg85=reg116*reg82; reg100=reg97*reg113; reg47=reg47*reg61; reg11=reg65+reg11;
    reg65=reg59*reg126; reg127=reg130+reg127; reg2=reg86*reg2; reg21=reg106+reg21; reg86=reg59*reg107;
    reg106=reg54*reg7; reg110=reg103*reg10; reg22=reg22*reg1; reg128=reg35*reg1; reg130=reg88*reg107;
    reg123=reg67+reg123; reg67=reg88*reg126; reg98=reg50+reg98; reg68=reg76+reg68; reg101=reg69+reg101;
    reg16=reg42+reg16; reg90=reg90*reg13; reg6=(*f.m).alpha_3*reg6; reg42=(*f.m).alpha_1*reg23; reg50=reg12*reg95;
    reg69=reg89*reg129; reg9=reg15+reg9; reg29=reg46*reg29; reg27=reg8+reg27; reg63=reg99+reg63;
    reg87=reg87*reg78; reg93=reg105*reg93; reg45=reg45*reg78; reg84=reg53+reg84; reg48=reg124+reg48;
    reg8=reg4*reg129; reg66=reg80+reg66; reg25=reg112*reg25; reg34=reg19+reg34; reg58=reg38*reg58;
    reg121=reg121*reg61; reg110=reg106+reg110; reg89=reg89*reg25; reg15=reg62*reg128; reg86=reg21+reg86;
    reg19=reg62*reg22; reg65=reg11+reg65; reg87=reg63+reg87; reg11=reg49*reg2; reg47=reg127+reg47;
    reg21=reg70*reg128; reg130=reg123+reg130; reg38=reg70*reg22; reg67=reg98+reg67; reg53=reg95*reg68;
    reg31=reg119+reg31; reg63=reg97*reg68; reg90=reg101+reg90; reg100=reg50-reg100; reg50=reg139*reg60;
    reg76=reg118*reg109; reg58=(*f.m).alpha_3*reg58; reg64=reg42+reg64; reg6=reg16+reg6; reg81=reg85+reg81;
    reg34=reg34*reg13; reg16=reg35*reg128; reg26=reg18+reg26; reg8=reg66+reg8; reg45=reg84+reg45;
    reg18=reg35*reg22; reg41=reg24+reg41; reg4=reg4*reg25; reg69=reg27+reg69; reg29=reg9+reg29;
    reg78=reg105*reg78; reg129=reg79*reg129; reg0=reg5*reg0; reg93=reg48+reg93; reg5=reg14*reg2;
    reg32=reg111+reg32; reg11=reg47+reg11; reg129=reg93+reg129; reg19=reg65+reg19; reg50=reg76-reg50;
    reg9=reg116*reg100; reg5=reg32+reg5; reg16=reg26+reg16; reg88=reg88*reg0; reg24=reg139*reg132;
    reg1=reg23*reg1; reg18=reg41+reg18; reg25=reg79*reg25; reg78=reg29+reg78; reg4=reg45+reg4;
    reg23=reg118*reg132; reg89=reg87+reg89; reg58=reg64+reg58; reg69=reg69*reg6; reg38=reg67+reg38;
    reg26=reg95*reg116; reg8=reg8*reg6; reg63=reg90+reg63; reg15=reg86+reg15; reg53=reg31+reg53;
    reg121=reg110+reg121; reg59=reg59*reg0; reg27=reg97*reg116; reg21=reg130+reg21; reg29=reg46*reg2;
    reg34=reg81+reg34; reg31=reg132*reg68; reg32=reg118*reg95; reg41=reg12*reg116; reg27=reg23-reg27;
    reg44=reg44*reg0; reg29=reg121+reg29; reg23=reg118*reg43; reg42=reg139*reg43; reg26=reg24-reg26;
    reg24=reg113*reg116; reg9=reg50+reg9; reg4=reg4*reg58; reg8=reg63+reg8; reg69=reg53+reg69;
    reg89=reg89*reg58; reg45=reg38*reg16; reg47=reg15*reg18; reg48=reg19*reg16; reg31=reg34+reg31;
    reg129=reg129*reg6; reg62=reg62*reg1; reg34=reg21*reg18; reg59=reg11+reg59; reg88=reg5+reg88;
    reg70=reg70*reg1; reg25=reg78+reg25; reg5=reg97*reg139; reg11=reg21*reg19; reg50=reg38*reg15;
    reg34=reg45-reg34; reg89=reg69+reg89; reg26=reg26/reg9; reg35=reg35*reg1; reg47=reg48-reg47;
    reg109=reg109/reg9; reg129=reg31+reg129; reg44=reg29+reg44; reg25=reg25*reg58; reg24=reg42-reg24;
    reg70=reg88+reg70; reg60=reg60/reg9; reg4=reg8+reg4; reg8=reg12*reg139; reg27=reg27/reg9;
    reg5=reg32-reg5; reg29=reg118*reg113; reg62=reg59+reg62; reg41=reg23-reg41; reg109=reg109*reg4;
    reg24=reg24/reg9; reg41=reg41/reg9; reg35=reg44+reg35; reg100=reg100/reg9; reg5=reg5/reg9;
    reg11=reg50-reg11; reg23=reg47*reg70; reg31=reg62*reg34; reg25=reg129+reg25; reg8=reg29-reg8;
    reg27=reg27*reg89; reg26=reg26*reg4; reg60=reg60*reg89; reg26=reg27-reg26; reg27=reg19*reg35;
    reg5=reg5*reg25; reg4=reg24*reg4; reg89=reg41*reg89; reg24=reg70*reg18; reg29=reg62*reg18;
    reg60=reg109-reg60; reg100=reg100*reg25; reg32=reg38*reg35; reg31=reg23-reg31; reg23=reg35*reg11;
    reg9=reg8/reg9; reg8=reg15*reg35; reg41=reg70*reg16; reg89=reg4-reg89; reg4=reg62*reg16;
    reg5=reg26-reg5; reg23=reg31+reg23; reg25=reg9*reg25; reg9=elem.pos(1)[1]-elem.pos(0)[1]; reg27=reg29-reg27;
    reg60=reg100+reg60; reg26=elem.pos(2)[0]-elem.pos(0)[0]; reg29=reg21*reg35; reg31=1-(*f.m).resolution; reg42=reg38*reg62;
    reg32=reg24-reg32; reg24=elem.pos(2)[1]-elem.pos(0)[1]; reg44=elem.pos(1)[0]-elem.pos(0)[0]; reg45=reg70*reg19; reg48=reg44*reg24;
    reg50=(*f.m).resolution*reg125; reg42=reg45-reg42; reg89=reg25+reg89; reg25=reg9*reg26; reg5=reg31*reg5;
    reg45=reg21*reg62; reg53=reg70*reg15; reg32=reg32/reg23; reg29=reg41-reg29; reg8=reg4-reg8;
    reg4=(*f.m).resolution*reg82; reg60=reg31*reg60; reg27=reg27/reg23; reg43=reg43*reg31; reg41=(*f.m).resolution*reg32;
    reg59=(*f.m).resolution*reg27; reg116=reg116*reg31; reg63=(*f.m).resolution*reg68; reg42=reg42/reg23; reg45=reg53-reg45;
    reg89=reg31*reg89; reg29=reg29/reg23; reg34=reg34/reg23; reg8=reg8/reg23; reg25=reg48-reg25;
    reg47=reg47/reg23; reg60=reg4+reg60; reg5=reg50+reg5; reg44=reg44/reg25; reg4=(*f.m).resolution*reg42;
    reg48=(*f.m).resolution*reg29; reg11=reg11/reg23; reg118=reg118*reg31; reg23=reg45/reg23; reg24=reg24/reg25;
    reg45=(*f.m).resolution*reg8; reg5=(*f.m).deltaT*reg5; reg116=reg59+reg116; reg41=reg43-reg41; reg43=(*f.m).resolution*reg34;
    reg60=(*f.m).deltaT*reg60; reg12=reg12*reg31; reg9=reg9/reg25; reg139=reg139*reg31; reg113=reg113*reg31;
    reg89=reg63+reg89; reg50=(*f.m).resolution*reg47; reg26=reg26/reg25; reg132=reg132*reg31; reg89=(*f.m).deltaT*reg89;
    reg113=reg48+reg113; reg48=reg116*reg60; reg43=reg12-reg43; reg45=reg139-reg45; reg132=reg4+reg132;
    reg4=reg9-reg24; reg12=reg26-reg44; reg53=reg41*reg5; reg59=(*f.m).resolution*reg11; reg63=(*f.m).resolution*reg23;
    reg118=reg50+reg118; reg95=reg95*reg31; reg31=reg97*reg31; reg50=reg48+reg53; reg63=reg95-reg63;
    reg31=reg59+reg31; reg59=reg113*reg5; reg64=reg132*reg89; reg65=reg45*reg60; reg66=reg43*reg5;
    reg67=reg118*reg60; reg69=0.5*reg12; reg76=0.5*reg4; reg78=0.5*reg26; reg80=0.5*reg24;
    reg81=0.5*reg44; reg84=0.5*reg9; reg85=reg81*reg132; reg86=reg67+reg66; reg87=reg31*reg89;
    reg88=reg12*reg41; reg90=reg9*reg116; reg93=reg63*reg89; reg95=reg76*reg132; reg97=reg26*reg41;
    reg98=reg80*reg132; reg99=reg84*reg132; reg100=reg65+reg59; reg101=reg4*reg116; reg106=reg78*reg132;
    reg109=reg24*reg116; reg110=reg69*reg132; reg111=reg44*reg41; reg119=reg50+reg64; reg121=reg81*reg63;
    reg123=reg26*reg43; reg124=reg9*reg45; reg127=reg76*reg31; reg129=reg12*reg43; reg110=reg101+reg110;
    reg101=reg69*reg31; reg130=reg4*reg118; reg131=reg26*reg113; reg133=reg80*reg63; reg134=2*reg119;
    reg135=reg76*reg63; reg85=reg85-reg90; reg136=reg12*reg113; reg137=reg81*reg31; reg138=reg86+reg87;
    reg139=reg69*reg63; reg140=reg4*reg45; reg111=reg111-reg99; reg109=reg109-reg106; reg141=reg9*reg118;
    T reg143=reg100+reg93; T reg144=reg78*reg31; T reg145=reg24*reg118; T reg146=reg84*reg63; T reg147=reg44*reg113;
    reg95=reg88+reg95; reg88=reg78*reg63; T reg148=reg24*reg45; T reg149=reg44*reg43; T reg150=reg84*reg31;
    reg98=reg98-reg97; T reg151=reg80*reg31; reg133=reg133-reg131; T reg152=reg24*reg138; T reg153=reg78*reg134;
    reg98=2*reg98; T reg154=reg80*reg134; reg151=reg151-reg123; T reg155=reg26*reg143; reg147=reg147-reg146;
    reg137=reg137-reg141; reg85=2*reg85; reg121=reg121-reg124; T reg156=reg81*reg134; reg101=reg130+reg101;
    reg110=2*reg110; reg127=reg129+reg127; reg129=reg9*reg138; reg130=reg44*reg143; reg135=reg136+reg135;
    reg136=reg84*reg134; reg95=2*reg95; reg145=reg145-reg144; reg109=2*reg109; reg111=2*reg111;
    reg139=reg140+reg139; reg148=reg148-reg88; reg149=reg149-reg150; reg140=reg76*reg134; T reg157=reg84*reg98;
    T reg158=reg12*reg143; T reg159=reg153-reg152; T reg160=reg155-reg154; T reg161=reg44*reg133; T reg162=reg136-reg130;
    T reg163=reg129-reg156; T reg164=reg69*reg134; T reg165=reg4*reg138; T reg166=reg84*reg111; T reg167=reg44*reg147;
    T reg168=reg84*reg85; T reg169=reg44*reg121; T reg170=reg76*reg85; T reg171=reg12*reg121; T reg172=reg76*reg111;
    T reg173=reg12*reg147; T reg174=reg9*reg145; T reg175=reg81*reg98; T reg176=reg9*reg151; T reg177=reg81*reg85;
    T reg178=reg9*reg137; T reg179=reg9*reg149; T reg180=reg44*reg139; T reg181=reg84*reg110; T reg182=reg44*reg135;
    T reg183=reg44*reg148; T reg184=reg84*reg109; T reg185=reg12*reg139; T reg186=reg76*reg110; T reg187=reg69*reg111;
    T reg188=reg4*reg145; T reg189=reg4*reg149; T reg190=reg69*reg95; T reg191=reg4*reg127; T reg192=reg12*reg148;
    T reg193=reg24*reg137; T reg194=reg69*reg110; T reg195=reg4*reg101; T reg196=reg81*reg109; T reg197=reg76*reg109;
    T reg198=reg78*reg110; T reg199=reg9*reg127; T reg200=reg81*reg95; T reg201=reg24*reg101; T reg202=reg78*reg85;
    reg101=reg9*reg101; T reg203=reg81*reg110; reg147=reg26*reg147; T reg204=reg80*reg111; T reg205=reg69*reg98;
    T reg206=reg4*reg151; reg121=reg26*reg121; T reg207=reg80*reg85; T reg208=reg12*reg135; T reg209=reg26*reg133;
    T reg210=reg80*reg98; T reg211=reg69*reg109; reg151=reg24*reg151; T reg212=reg78*reg95; T reg213=reg81*reg111;
    T reg214=reg84*reg95; T reg215=reg78*reg98; reg148=reg26*reg148; T reg216=reg80*reg109; reg109=reg78*reg109;
    reg135=reg26*reg135; T reg217=reg80*reg95; reg145=reg24*reg145; reg139=reg26*reg139; reg110=reg80*reg110;
    reg133=reg12*reg133; reg137=reg4*reg137; reg95=reg76*reg95; reg111=reg78*reg111; reg149=reg24*reg149;
    reg85=reg69*reg85; reg98=reg76*reg98; reg127=reg24*reg127; reg160=reg25*reg160; reg215=reg151-reg215;
    reg184=reg183-reg184; reg151=reg165+reg164; reg202=reg193-reg202; reg174=reg196-reg174; reg157=reg161-reg157;
    reg163=reg25*reg163; reg166=reg167-reg166; reg173=reg172+reg173; reg188=reg211+reg188; reg85=reg137+reg85;
    reg176=reg175-reg176; reg205=reg206+reg205; reg178=reg177-reg178; reg168=reg169-reg168; reg109=reg145-reg109;
    reg179=reg213-reg179; reg214=reg182-reg214; reg198=reg201-reg198; reg162=reg25*reg162; reg159=reg25*reg159;
    reg137=reg140+reg158; reg212=reg127-reg212; reg181=reg180-reg181; reg135=reg217-reg135; reg199=reg200-reg199;
    reg148=reg216-reg148; reg133=reg98+reg133; reg139=reg110-reg139; reg111=reg149-reg111; reg185=reg186+reg185;
    reg187=reg189+reg187; reg192=reg197+reg192; reg190=reg191+reg190; reg194=reg195+reg194; reg101=reg203-reg101;
    reg147=reg204-reg147; reg121=reg207-reg121; reg171=reg170+reg171; reg95=reg208+reg95; reg209=reg210-reg209;
    reg85=reg25*reg85; reg174=reg25*reg174; reg185=reg25*reg185; reg215=reg25*reg215; reg202=reg25*reg202;
    reg209=reg25*reg209; reg98=reg25*reg137; reg111=reg25*reg111; reg212=reg25*reg212; reg162=ponderation*reg162;
    reg159=ponderation*reg159; reg139=reg25*reg139; reg135=reg25*reg135; reg133=reg25*reg133; reg173=reg25*reg173;
    reg160=ponderation*reg160; reg148=reg25*reg148; reg109=reg25*reg109; reg171=reg25*reg171; reg188=reg25*reg188;
    reg163=ponderation*reg163; reg166=reg25*reg166; reg95=reg25*reg95; reg176=reg25*reg176; reg179=reg25*reg179;
    reg101=reg25*reg101; reg168=reg25*reg168; reg199=reg25*reg199; reg181=reg25*reg181; reg194=reg25*reg194;
    reg121=reg25*reg121; reg190=reg25*reg190; reg214=reg25*reg214; reg178=reg25*reg178; reg198=reg25*reg198;
    reg147=reg25*reg147; reg184=reg25*reg184; reg157=reg25*reg157; reg205=reg25*reg205; reg187=reg25*reg187;
    reg110=reg25*reg151; reg192=reg25*reg192; matrix(indices[0]+1,indices[1]+1)+=ponderation*reg133; matrix(indices[2]+1,indices[2]+0)+=ponderation*reg168; matrix(indices[0]+1,indices[0]+1)+=ponderation*reg95;
    matrix(indices[2]+1,indices[1]+1)+=ponderation*reg157; sollicitation[indices[2]+1]+=-reg162; matrix(indices[2]+1,indices[2]+1)+=ponderation*reg166; matrix(indices[0]+1,indices[1]+0)+=ponderation*reg192; matrix(indices[1]+0,indices[1]+1)+=ponderation*reg215;
    matrix(indices[0]+0,indices[1]+0)+=ponderation*reg188; matrix(indices[0]+1,indices[2]+1)+=ponderation*reg173; matrix(indices[0]+1,indices[2]+0)+=ponderation*reg171; matrix(indices[1]+1,indices[1]+1)+=ponderation*reg209; matrix(indices[1]+1,indices[2]+0)+=ponderation*reg121;
    matrix(indices[1]+1,indices[2]+1)+=ponderation*reg147; matrix(indices[2]+0,indices[0]+0)+=ponderation*reg101; matrix(indices[2]+0,indices[0]+1)+=ponderation*reg199; matrix(indices[0]+0,indices[0]+0)+=ponderation*reg194; matrix(indices[0]+0,indices[0]+1)+=ponderation*reg190;
    matrix(indices[0]+0,indices[2]+1)+=ponderation*reg187; matrix(indices[0]+1,indices[0]+0)+=ponderation*reg185; matrix(indices[1]+0,indices[2]+0)+=ponderation*reg202; matrix(indices[1]+0,indices[2]+1)+=ponderation*reg111; matrix(indices[1]+1,indices[0]+0)+=ponderation*reg139;
    matrix(indices[1]+1,indices[0]+1)+=ponderation*reg135; matrix(indices[1]+1,indices[1]+0)+=ponderation*reg148; sollicitation[indices[2]+0]+=-reg163; matrix(indices[1]+0,indices[1]+0)+=ponderation*reg109; sollicitation[indices[1]+1]+=-reg160;
    sollicitation[indices[1]+0]+=-reg159; matrix(indices[1]+0,indices[0]+1)+=ponderation*reg212; reg95=ponderation*reg98; sollicitation[indices[0]+1]+=reg95; matrix(indices[1]+0,indices[0]+0)+=ponderation*reg198;
    matrix(indices[0]+0,indices[2]+0)+=ponderation*reg85; reg85=ponderation*reg110; sollicitation[indices[0]+0]+=reg85; matrix(indices[2]+1,indices[1]+0)+=ponderation*reg184; matrix(indices[2]+1,indices[0]+1)+=ponderation*reg214;
    matrix(indices[2]+1,indices[0]+0)+=ponderation*reg181; matrix(indices[2]+0,indices[2]+1)+=ponderation*reg179; matrix(indices[0]+0,indices[1]+1)+=ponderation*reg205; matrix(indices[2]+0,indices[2]+0)+=ponderation*reg178; matrix(indices[2]+0,indices[1]+1)+=ponderation*reg176;
    matrix(indices[2]+0,indices[1]+0)+=ponderation*reg174;
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
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg5=pow((*f.m).v2[0],2); T reg6=pow((*f.m).v2[1],2);
    T reg7=1.0/(*f.m).elastic_modulus_3; T reg8=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg9=reg2*reg3; T reg10=pow((*f.m).v1[1],2); T reg11=pow((*f.m).v1[0],2);
    T reg12=pow((*f.m).v1[2],2); T reg13=1.0/(*f.m).elastic_modulus_2; T reg14=reg7*reg9; reg10=reg11+reg10; reg11=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1;
    T reg15=pow((*f.m).v2[2],2); reg6=reg5+reg6; reg5=reg8*reg9; T reg16=reg4*reg9; T reg17=reg8*reg5;
    T reg18=reg11*reg14; T reg19=reg8*reg16; reg15=reg6+reg15; reg6=reg13*reg14; reg10=reg12+reg10;
    reg19=reg18+reg19; reg12=reg11*reg5; T reg20=reg13*reg16; reg17=reg6-reg17; reg6=1.0/(*f.m).elastic_modulus_1;
    reg10=pow(reg10,0.5); reg15=pow(reg15,0.5); T reg21=(*f.m).v2[2]/reg15; T reg22=reg6*reg17; T reg23=(*f.m).v2[1]/reg15;
    T reg24=reg12+reg20; T reg25=(*f.m).v1[1]/reg10; T reg26=(*f.m).v1[2]/reg10; T reg27=reg11*reg19; reg10=(*f.m).v1[0]/reg10;
    T reg28=reg4*reg24; T reg29=reg4*reg3; reg27=reg22-reg27; reg22=reg13*reg9; T reg30=reg4*reg5;
    T reg31=reg8*reg3; T reg32=reg26*reg23; T reg33=reg25*reg21; T reg34=reg2*reg0; reg14=reg6*reg14;
    reg15=(*f.m).v2[0]/reg15; T reg35=reg4*reg16; reg9=reg11*reg9; reg3=reg7*reg3; T reg36=reg8*reg31;
    T reg37=reg11*reg3; T reg38=reg10*reg21; T reg39=reg26*reg15; T reg40=reg33-reg32; reg28=reg27-reg28;
    reg27=2*reg10; T reg41=2*reg15; T reg42=reg8*reg34; T reg43=reg2*reg1; T reg44=reg4*reg34;
    reg16=reg11*reg16; T reg45=reg8*reg29; reg5=reg6*reg5; reg34=reg7*reg34; reg3=reg13*reg3;
    reg35=reg14-reg35; reg14=reg4*reg9; reg30=reg18+reg30; reg18=reg4*reg22; T reg46=reg7*reg43;
    T reg47=reg8*reg44; T reg48=reg8*reg43; reg43=reg4*reg43; T reg49=reg8*reg42; T reg50=reg11*reg34;
    reg34=reg13*reg34; reg45=reg37+reg45; reg36=reg3-reg36; reg29=reg13*reg29; reg3=reg25*reg27;
    reg37=pow(reg25,2); T reg51=pow(reg10,2); T reg52=reg39-reg38; T reg53=reg10*reg23; T reg54=reg25*reg15;
    reg17=reg17/reg28; reg31=reg11*reg31; reg19=reg19/reg28; reg9=reg11*reg9; reg16=reg16+reg5;
    reg14=reg5+reg14; reg30=reg30/reg28; reg18=reg12+reg18; reg35=reg35/reg28; reg22=reg6*reg22;
    reg5=pow(reg15,2); T reg55=pow(reg23,2); T reg56=reg41*reg23; T reg57=2*reg40; T reg58=reg5*reg19;
    reg18=reg18/reg28; T reg59=reg31+reg29; T reg60=reg56*reg19; T reg61=reg17*reg3; T reg62=reg11*reg46;
    reg46=reg13*reg46; T reg63=reg30*reg51; T reg64=reg5*reg35; reg9=reg22-reg9; reg47=reg50+reg47;
    reg16=reg16/reg28; reg22=reg30*reg37; reg14=reg14/reg28; reg50=reg55*reg35; reg49=reg34-reg49;
    reg24=reg24/reg28; reg34=reg56*reg35; T reg65=reg30*reg3; T reg66=reg17*reg51; T reg67=pow(reg26,2);
    T reg68=pow(reg21,2); T reg69=pow(reg40,2); T reg70=pow(reg52,2); reg42=reg11*reg42; T reg71=reg8*reg43;
    T reg72=reg57*reg52; T reg73=reg17*reg37; T reg74=reg53-reg54; reg36=reg6*reg36; T reg75=reg55*reg19;
    T reg76=reg8*reg48; reg45=reg11*reg45; reg44=reg13*reg44; T reg77=pow(reg74,2); reg50=reg22+reg50;
    reg48=reg11*reg48; reg22=reg70*reg16; T reg78=reg30*reg67; reg66=reg58+reg66; reg58=reg69*reg16;
    reg9=reg9/reg28; reg64=reg63+reg64; reg63=reg69*reg24; reg43=reg13*reg43; reg45=reg36-reg45;
    reg36=reg72*reg24; reg59=reg4*reg59; reg60=reg61+reg60; reg61=reg55*reg14; T reg79=reg18*reg37;
    reg76=reg46-reg76; reg46=reg5*reg14; T reg80=reg18*reg51; T reg81=reg18*reg3; T reg82=reg56*reg14;
    reg71=reg62+reg71; reg62=reg17*reg67; T reg83=reg42+reg44; reg47=reg11*reg47; T reg84=reg68*reg19;
    T reg85=reg68*reg35; reg75=reg73+reg75; reg73=reg70*reg24; T reg86=reg72*reg16; reg34=reg65+reg34;
    reg49=reg6*reg49; reg65=reg26*reg27; reg73=reg75+reg73; reg61=reg79+reg61; reg75=2*reg25;
    reg79=reg68*reg14; T reg87=reg70*reg9; reg71=reg11*reg71; reg63=reg66+reg63; reg66=reg10*reg25;
    T reg88=reg15*reg23; T reg89=2*reg23; T reg90=reg41*reg21; T reg91=reg18*reg67; T reg92=reg48+reg43;
    reg85=reg78+reg85; reg78=reg77*reg16; T reg93=reg24*reg77; reg22=reg50+reg22; reg86=reg34+reg86;
    reg84=reg62+reg84; reg47=reg49-reg47; reg58=reg64+reg58; reg83=reg4*reg83; reg34=reg72*reg9;
    reg82=reg81+reg82; reg59=reg45-reg59; reg46=reg80+reg46; reg45=reg69*reg9; reg76=reg6*reg76;
    reg36=reg60+reg36; reg83=reg47-reg83; reg59=reg59/reg28; reg92=reg4*reg92; reg71=reg76-reg71;
    reg47=reg55*reg86; reg93=reg84+reg93; reg49=reg65*reg17; reg50=reg37*reg73; reg60=reg51*reg73;
    reg62=reg5*reg22; reg64=reg65*reg30; reg76=reg90*reg35; reg80=reg88*reg86; reg81=reg66*reg36;
    reg45=reg46+reg45; reg46=reg15*reg27; reg84=reg75*reg23; T reg94=reg55*reg58; T reg95=reg37*reg63;
    T reg96=2*reg26; reg87=reg61+reg87; reg61=reg51*reg36; T reg97=reg52*reg40; reg86=reg5*reg86;
    reg73=reg66*reg73; T reg98=reg88*reg22; reg79=reg91+reg79; reg91=reg77*reg9; reg53=reg54+reg53;
    reg36=reg37*reg36; reg54=reg25*reg23; T reg99=reg10*reg15; T reg100=reg90*reg19; T reg101=reg25*reg40;
    T reg102=reg10*reg52; T reg103=2*reg52; reg57=reg57*reg74; T reg104=reg89*reg21; T reg105=reg75*reg26;
    reg34=reg82+reg34; reg82=reg51*reg63; T reg106=reg5*reg58; reg78=reg85+reg78; reg22=reg55*reg22;
    reg91=reg79+reg91; reg100=reg49+reg100; reg49=reg57*reg24; reg17=reg105*reg17; reg19=reg104*reg19;
    reg79=reg90*reg14; reg35=reg104*reg35; reg30=reg105*reg30; reg85=reg57*reg16; reg76=reg64+reg76;
    reg64=reg65*reg18; T reg107=reg55*reg13; T reg108=reg5*reg11; reg106=reg82+reg106; reg82=reg69*reg45;
    T reg109=reg84*reg11; T reg110=reg46*reg6; T reg111=reg70*reg87; reg22=reg50+reg22; reg50=reg55*reg11;
    T reg112=reg5*reg6; reg62=reg60+reg62; reg60=reg97*reg34; reg80=reg81+reg80; reg81=reg69*reg87;
    T reg113=reg23*reg40; T reg114=reg15*reg52; T reg115=reg70*reg45; reg94=reg95+reg94; reg95=reg51*reg93;
    T reg116=reg5*reg78; reg63=reg66*reg63; reg58=reg88*reg58; reg96=reg96*reg21; reg86=reg61+reg86;
    reg98=reg73+reg98; reg87=reg97*reg87; reg61=reg69*reg34; reg92=reg71-reg92; reg83=reg83/reg28;
    reg34=reg70*reg34; reg47=reg36+reg47; reg36=reg99*reg59; reg71=reg54*reg59; reg73=reg53*reg59;
    reg103=reg103*reg74; T reg117=reg26*reg21; T reg118=reg10*reg40; T reg119=reg25*reg52; T reg120=reg37*reg93;
    T reg121=reg55*reg78; T reg122=reg46*reg11; T reg123=reg84*reg13; reg102=reg101+reg102; reg61=reg86+reg61;
    reg86=reg57*reg9; reg101=reg84*reg73; reg18=reg105*reg18; T reg124=reg46*reg73; T reg125=reg84*reg71;
    reg14=reg104*reg14; T reg126=reg69*reg91; reg116=reg95+reg116; reg34=reg47+reg34; reg47=reg119*reg83;
    reg111=reg22+reg111; reg121=reg120+reg121; reg22=reg46*reg36; reg115=reg94+reg115; reg94=reg84*reg36;
    reg95=reg46*reg71; reg120=reg117*reg59; T reg127=reg70*reg91; reg81=reg62+reg81; reg82=reg106+reg82;
    reg62=reg102*reg83; reg106=reg118*reg83; reg79=reg64+reg79; reg64=reg26*reg74; reg75=reg75*reg52;
    reg122=reg123-reg122; reg123=reg96*reg8; reg27=reg40*reg27; reg108=reg107-reg108; reg107=reg15*reg40;
    T reg128=reg23*reg52; T reg129=reg68*reg8; reg93=reg66*reg93; reg78=reg88*reg78; reg109=reg110-reg109;
    reg110=reg4*reg96; reg50=reg112-reg50; reg112=reg68*reg4; reg114=reg113+reg114; reg73=reg53*reg73;
    reg60=reg80+reg60; reg71=reg53*reg71; reg16=reg103*reg16; reg35=reg30+reg35; reg85=reg76+reg85;
    reg87=reg98+reg87; reg24=reg103*reg24; reg19=reg17+reg19; reg49=reg100+reg49; reg28=reg92/reg28;
    reg45=reg97*reg45; reg58=reg63+reg58; reg17=reg5*reg4; reg30=reg55*reg8; reg63=reg84*reg8;
    reg76=reg4*reg46; reg71=reg87+reg71; reg80=reg5*reg85; reg87=reg51*reg49; reg94=reg115+reg94;
    reg92=reg102*reg47; reg91=reg97*reg91; reg98=reg27*reg62; reg78=reg93+reg78; reg36=reg53*reg36;
    reg45=reg58+reg45; reg124=reg61+reg124; reg58=reg55*reg85; reg61=reg37*reg49; reg101=reg34+reg101;
    reg34=reg75*reg62; reg93=reg55*(*f.m).alpha_2; reg100=(*f.m).alpha_1*reg37; reg113=reg5*(*f.m).alpha_2; reg115=(*f.m).alpha_1*reg51;
    reg63=reg76+reg63; reg96=reg96*reg7; reg30=reg17+reg30; reg17=reg68*reg7; reg123=reg122-reg123;
    reg76=reg84*reg120; reg127=reg121+reg127; reg129=reg108-reg129; reg108=reg75*reg47; reg125=reg111+reg125;
    reg13=reg13*reg37; reg111=reg11*reg51; reg110=reg109-reg110; reg112=reg50-reg112; reg11=reg11*reg37;
    reg62=reg102*reg62; reg50=reg75*reg106; reg73=reg60+reg73; reg47=reg27*reg47; reg60=reg26*reg40;
    reg109=reg10*reg74; reg95=reg81+reg95; reg126=reg116+reg126; reg81=reg46*reg120; reg116=reg27*reg106;
    reg89=reg89*reg52; reg41=reg41*reg40; reg22=reg82+reg22; reg6=reg6*reg51; reg24=reg19+reg24;
    reg39=reg38+reg39; reg19=reg114*reg28; reg16=reg35+reg16; reg86=reg79+reg86; reg14=reg18+reg14;
    reg9=reg103*reg9; reg18=reg128*reg28; reg35=reg107*reg28; reg38=reg21*reg74; reg79=reg64*reg83;
    reg36=reg45+reg36; reg113=reg115+reg113; reg30=reg17-reg30; reg17=reg39*reg59; reg45=reg70*reg86;
    reg58=reg61+reg58; reg11=reg6-reg11; reg63=reg96-reg63; reg6=reg89*reg19; reg34=reg101+reg34;
    reg9=reg14+reg9; reg14=reg25*reg74; reg62=reg73+reg62; reg61=reg5*reg112; reg73=reg55*reg129;
    reg82=reg26*reg52; reg96=reg114*reg19; reg101=reg5*reg110; reg115=reg55*reg123; reg121=reg37*reg123;
    reg122=reg51*reg110; reg49=reg66*reg49; reg85=reg88*reg85; reg109=reg60+reg109; reg60=reg99*reg112;
    T reg130=reg54*reg129; T reg131=reg37*reg129; reg111=reg13-reg111; reg13=reg51*reg112; T reg132=reg21*reg40;
    T reg133=reg55*reg16; T reg134=reg37*reg24; T reg135=reg8*reg67; T reg136=reg99*reg110; T reg137=reg54*reg123;
    T reg138=reg15*reg74; reg47=reg95+reg47; reg95=reg41*reg18; reg8=reg8*reg37; T reg139=reg70*(*f.m).alpha_3;
    reg33=reg32+reg33; reg32=(*f.m).alpha_1*reg67; T reg140=reg68*(*f.m).alpha_2; T reg141=reg27*reg79; T reg142=reg4*reg67;
    T reg143=reg5*reg16; T reg144=reg51*reg24; reg92=reg71+reg92; reg71=reg114*reg18; reg81=reg126+reg81;
    reg80=reg87+reg80; reg4=reg4*reg51; reg87=reg69*reg86; reg106=reg102*reg106; reg126=reg75*reg79;
    reg76=reg127+reg76; reg127=reg69*(*f.m).alpha_3; T reg145=reg38*reg28; reg18=reg89*reg18; reg108=reg125+reg108;
    reg93=reg100+reg93; reg116=reg22+reg116; reg22=reg41*reg35; reg120=reg53*reg120; reg100=reg89*reg35;
    reg98=reg124+reg98; reg50=reg94+reg50; reg19=reg41*reg19; reg91=reg78+reg91; reg78=reg67*reg63;
    reg120=reg91+reg120; reg121=reg122+reg121; reg79=reg102*reg79; reg138=reg132+reg138; reg91=reg21*reg52;
    reg94=reg23*reg74; reg71=reg92+reg71; reg92=reg67*reg30; reg131=reg13+reg131; reg10=reg10*reg26;
    reg15=reg15*reg21; reg13=reg53*reg2; reg7=reg7*reg67; reg122=reg88*reg2; reg8=reg4+reg8;
    reg35=reg114*reg35; reg106=reg36+reg106; reg126=reg76+reg126; reg4=reg109*reg83; reg18=reg108+reg18;
    reg127=reg113+reg127; reg100=reg50+reg100; reg22=reg116+reg22; reg95=reg47+reg95; reg139=reg93+reg139;
    reg140=reg32+reg140; reg32=reg77*(*f.m).alpha_3; reg36=reg69*reg9; reg143=reg144+reg143; reg47=(*f.m).alpha_1*reg66;
    reg50=(*f.m).alpha_2*reg88; reg24=reg66*reg24; reg16=reg88*reg16; reg76=reg46*reg17; reg87=reg80+reg87;
    reg19=reg98+reg19; reg141=reg81+reg141; reg80=reg41*reg145; reg14=reg82+reg14; reg73=reg61+reg73;
    reg61=reg68*reg30; reg96=reg62+reg96; reg115=reg101+reg115; reg62=reg68*reg63; reg85=reg49+reg85;
    reg86=reg97*reg86; reg135=reg111-reg135; reg130=reg60+reg130; reg49=reg117*reg30; reg60=reg70*reg9;
    reg133=reg134+reg133; reg137=reg136+reg137; reg81=reg117*reg63; reg82=reg84*reg17; reg45=reg58+reg45;
    reg142=reg11-reg142; reg6=reg34+reg6; reg59=reg33*reg59; reg11=reg89*reg145; reg34=reg22*reg127;
    reg58=reg95*reg139; reg88=(*f.m).alpha_2*reg15; reg93=(*f.m).alpha_1*reg10; reg98=(*f.m).alpha_3*reg97; reg50=reg47+reg50;
    reg32=reg140+reg32; reg47=reg53*reg13; reg81=reg137+reg81; reg101=reg53*reg122; reg49=reg130+reg49;
    reg108=reg18*reg96; reg111=reg95*reg96; reg113=reg6*reg71; reg116=reg56*reg13; reg62=reg115+reg62;
    reg115=reg56*reg122; reg61=reg73+reg61; reg73=reg55*reg135; reg124=reg5*reg142; reg125=reg19*reg71;
    reg130=reg3*reg13; reg78=reg121+reg78; reg121=reg3*reg122; reg92=reg131+reg92; reg131=reg37*reg135;
    reg132=reg51*reg142; reg134=reg39*reg1; reg15=reg15*reg1; reg8=reg7-reg8; reg86=reg85+reg86;
    reg84=reg84*reg59; reg60=reg133+reg60; reg17=reg53*reg17; reg7=reg100*reg127; reg85=reg18*reg139;
    reg133=reg75*reg4; reg82=reg45+reg82; reg145=reg114*reg145; reg11=reg126+reg11; reg83=reg14*reg83;
    reg16=reg24+reg16; reg24=reg138*reg28; reg79=reg120+reg79; reg9=reg97*reg9; reg35=reg106+reg35;
    reg80=reg141+reg80; reg76=reg87+reg76; reg45=reg27*reg4; reg21=reg23*reg21; reg26=reg25*reg26;
    reg36=reg143+reg36; reg46=reg46*reg59; reg40=reg74*reg40; reg94=reg91+reg94; reg125=reg111-reg125;
    reg23=reg95*reg6; reg73=reg124+reg73; reg25=reg19*reg18; reg68=reg68*reg8; reg87=reg35*reg127;
    reg91=reg65*reg134; reg130=reg78+reg130; reg78=reg65*reg15; reg121=reg92+reg121; reg67=reg67*reg8;
    reg131=reg132+reg131; reg92=reg33*reg0; reg97=reg21*reg0; reg2=reg66*reg2; reg66=(*f.m).alpha_1*reg26;
    reg40=(*f.m).alpha_3*reg40; reg88=reg93+reg88; reg98=reg50+reg98; reg45=reg76+reg45; reg50=reg41*reg24;
    reg46=reg36+reg46; reg27=reg27*reg83; reg28=reg94*reg28; reg21=(*f.m).alpha_2*reg21; reg36=reg39*reg134;
    reg47=reg81+reg47; reg133=reg82+reg133; reg76=reg89*reg24; reg81=reg39*reg15; reg101=reg49+reg101;
    reg84=reg60+reg84; reg75=reg75*reg83; reg49=reg54*reg135; reg60=reg99*reg142; reg113=reg108-reg113;
    reg82=reg90*reg134; reg116=reg62+reg116; reg62=reg90*reg15; reg115=reg61+reg115; reg61=reg71*reg139;
    reg80=reg80*reg32; reg85=reg7+reg85; reg58=reg34+reg58; reg17=reg86+reg17; reg11=reg11*reg32;
    reg145=reg79+reg145; reg52=reg74*reg52; reg9=reg16+reg9; reg59=reg53*reg59; reg4=reg102*reg4;
    reg27=reg46+reg27; reg41=reg41*reg28; reg89=reg89*reg28; reg21=reg66+reg21; reg52=(*f.m).alpha_3*reg52;
    reg91=reg130+reg91; reg7=reg105*reg92; reg25=reg23-reg25; reg16=reg33*reg92; reg36=reg47+reg36;
    reg68=reg73+reg68; reg23=reg56*reg2; reg61=reg87+reg61; reg145=reg145*reg32; reg49=reg60+reg49;
    reg117=reg117*reg8; reg34=reg33*reg97; reg62=reg115+reg62; reg46=reg104*reg97; reg75=reg84+reg75;
    reg81=reg101+reg81; reg76=reg133+reg76; reg82=reg116+reg82; reg47=reg104*reg92; reg60=reg100*reg125;
    reg66=reg22*reg113; reg24=reg114*reg24; reg4=reg17+reg4; reg17=reg6*reg98; reg11=reg85+reg11;
    reg50=reg45+reg50; reg67=reg131+reg67; reg45=reg3*reg2; reg1=reg10*reg1; reg83=reg102*reg83;
    reg10=reg19*reg98; reg80=reg58+reg80; reg59=reg9+reg59; reg78=reg121+reg78; reg9=reg105*reg97;
    reg40=reg88+reg40; reg24=reg4+reg24; reg65=reg65*reg1; reg45=reg67+reg45; reg23=reg68+reg23;
    reg90=reg90*reg1; reg145=reg61+reg145; reg4=reg96*reg98; reg34=reg81+reg34; reg46=reg62+reg46;
    reg0=reg26*reg0; reg26=reg19*reg35; reg47=reg82+reg47; reg89=reg75+reg89; reg28=reg114*reg28;
    reg58=reg53*reg2; reg117=reg49+reg117; reg83=reg59+reg83; reg52=reg21+reg52; reg21=reg22*reg96;
    reg49=reg6*reg35; reg59=reg100*reg96; reg17=reg11+reg17; reg11=reg35*reg25; reg60=reg66-reg60;
    reg7=reg91+reg7; reg76=reg76*reg40; reg9=reg78+reg9; reg50=reg50*reg40; reg16=reg36+reg16;
    reg41=reg27+reg41; reg10=reg80+reg10; reg27=reg95*reg35; reg50=reg10+reg50; reg10=reg46*reg16;
    reg36=reg18*reg35; reg26=reg21-reg26; reg89=reg89*reg52; reg41=reg41*reg52; reg76=reg17+reg76;
    reg17=reg47*reg34; reg58=reg117+reg58; reg28=reg83+reg28; reg39=reg39*reg1; reg65=reg45+reg65;
    reg21=reg7*reg34; reg105=reg105*reg0; reg11=reg60+reg11; reg90=reg23+reg90; reg104=reg104*reg0;
    reg23=reg19*reg100; reg45=reg22*reg71; reg60=reg100*reg71; reg4=reg145+reg4; reg24=reg24*reg40;
    reg61=reg9*reg16; reg62=reg22*reg6; reg49=reg59-reg49; reg59=reg9*reg47; reg66=reg7*reg46;
    reg21=reg61-reg21; reg61=reg22*reg18; reg33=reg33*reg0; reg39=reg58+reg39; reg23=reg62-reg23;
    reg58=reg95*reg100; reg27=reg45-reg27; reg113=reg113/reg11; reg105=reg65+reg105; reg89=reg76+reg89;
    reg104=reg90+reg104; reg17=reg10-reg17; reg41=reg50+reg41; reg125=reg125/reg11; reg49=reg49/reg11;
    reg24=reg4+reg24; reg28=reg28*reg52; reg26=reg26/reg11; reg36=reg60-reg36; reg23=reg23/reg11;
    reg33=reg39+reg33; reg58=reg61-reg58; reg25=reg25/reg11; reg28=reg24+reg28; reg4=reg104*reg21;
    reg26=reg26*reg89; reg27=reg27/reg11; reg125=reg125*reg89; reg10=reg17*reg105; reg66=reg59-reg66;
    reg49=reg49*reg41; reg36=reg36/reg11; reg113=reg113*reg41; reg4=reg10-reg4; reg10=reg33*reg66;
    reg24=reg9*reg33; reg39=reg104*reg34; reg45=reg46*reg33; reg50=reg105*reg34; reg125=reg113-reg125;
    reg25=reg25*reg28; reg49=reg26-reg49; reg41=reg36*reg41; reg89=reg27*reg89; reg11=reg58/reg11;
    reg23=reg23*reg28; reg45=reg39-reg45; reg26=reg105*reg16; reg27=reg47*reg33; reg36=reg104*reg16;
    reg39=elem.pos(2)[0]-elem.pos(0)[0]; reg10=reg4+reg10; reg4=elem.pos(1)[1]-elem.pos(0)[1]; reg58=elem.pos(1)[0]-elem.pos(0)[0]; reg59=reg9*reg104;
    reg60=reg105*reg46; reg125=reg25+reg125; reg24=reg50-reg24; reg23=reg49-reg23; reg28=reg11*reg28;
    reg89=reg41-reg89; reg11=1-(*f.m).resolution; reg25=elem.pos(2)[1]-elem.pos(0)[1]; reg41=reg7*reg33; reg49=reg7*reg104;
    reg50=reg105*reg47; reg23=reg11*reg23; reg59=reg60-reg59; reg45=reg45/reg10; reg60=(*f.m).resolution*reg139;
    reg41=reg26-reg41; reg89=reg28+reg89; reg125=reg11*reg125; reg24=reg24/reg10; reg27=reg36-reg27;
    reg26=(*f.m).resolution*reg127; reg28=reg58*reg25; reg36=reg4*reg39; reg59=reg59/reg10; reg49=reg50-reg49;
    reg17=reg17/reg10; reg41=reg41/reg10; reg27=reg27/reg10; reg21=reg21/reg10; reg71=reg71*reg11;
    reg89=reg11*reg89; reg35=reg35*reg11; reg50=(*f.m).resolution*reg98; reg61=(*f.m).resolution*reg45; reg62=(*f.m).resolution*reg24;
    reg125=reg26+reg125; reg23=reg60+reg23; reg36=reg28-reg36; reg26=(*f.m).resolution*reg41; reg28=(*f.m).resolution*reg27;
    reg22=reg22*reg11; reg66=reg66/reg10; reg60=(*f.m).resolution*reg17; reg65=(*f.m).resolution*reg59; reg10=reg49/reg10;
    reg95=reg95*reg11; reg100=reg100*reg11; reg49=(*f.m).resolution*reg21; reg58=reg58/reg36; reg25=reg25/reg36;
    reg39=reg39/reg36; reg4=reg4/reg36; reg23=(*f.m).deltaT*reg23; reg125=(*f.m).deltaT*reg125; reg89=reg50+reg89;
    reg62=reg71-reg62; reg35=reg61+reg35; reg96=reg96*reg11; reg18=reg18*reg11; reg49=reg95-reg49;
    reg89=(*f.m).deltaT*reg89; reg22=reg60+reg22; reg6=reg6*reg11; reg50=reg4-reg25; reg60=reg39-reg58;
    reg61=(*f.m).resolution*reg66; reg96=reg65+reg96; reg18=reg26+reg18; reg11=reg19*reg11; reg28=reg100-reg28;
    reg19=(*f.m).resolution*reg10; reg26=reg62*reg23; reg65=reg35*reg125; reg67=reg96*reg89; reg68=reg65+reg26;
    reg71=reg49*reg23; reg73=0.5*reg60; reg74=reg22*reg125; reg75=reg28*reg125; reg76=reg18*reg23;
    reg78=0.5*reg50; reg79=0.5*reg39; reg80=0.5*reg25; reg81=0.5*reg58; reg82=0.5*reg4;
    reg11=reg61+reg11; reg19=reg6-reg19; reg6=reg82*reg96; reg61=reg60*reg62; reg83=reg73*reg96;
    reg84=reg50*reg35; reg85=reg80*reg96; reg86=reg11*reg89; reg87=reg74+reg71; reg88=reg39*reg62;
    reg90=reg4*reg35; reg91=reg81*reg96; reg93=reg78*reg96; reg95=reg25*reg35; reg100=reg79*reg96;
    reg101=reg58*reg62; reg106=reg68+reg67; reg108=reg19*reg89; reg111=reg75+reg76; reg113=reg111+reg108;
    reg95=reg95-reg100; reg115=reg80*reg19; reg116=reg79*reg19; reg117=reg39*reg18; reg120=reg25*reg28;
    reg93=reg61+reg93; reg61=reg87+reg86; reg121=reg25*reg22; reg124=reg79*reg11; reg126=reg78*reg11;
    reg130=reg50*reg22; reg131=reg60*reg18; reg132=reg78*reg19; reg133=reg73*reg11; reg136=2*reg106;
    reg137=reg4*reg28; reg140=reg82*reg11; reg141=reg58*reg49; reg143=reg81*reg19; reg144=reg39*reg49;
    reg145=reg80*reg11; T reg146=reg60*reg49; T reg147=reg82*reg19; T reg148=reg58*reg18; reg85=reg85-reg88;
    reg83=reg84+reg83; reg84=reg4*reg22; T reg149=reg81*reg11; reg101=reg101-reg6; reg91=reg91-reg90;
    reg101=2*reg101; reg95=2*reg95; T reg150=reg82*reg136; T reg151=reg58*reg113; T reg152=reg4*reg61;
    T reg153=reg81*reg136; T reg154=reg39*reg113; reg85=2*reg85; reg145=reg145-reg144; reg126=reg146+reg126;
    reg83=2*reg83; reg149=reg149-reg84; reg148=reg148-reg147; reg91=2*reg91; reg115=reg115-reg117;
    reg120=reg120-reg116; reg146=reg25*reg61; reg143=reg143-reg137; reg141=reg141-reg140; T reg155=reg79*reg136;
    reg133=reg130+reg133; reg130=reg80*reg136; reg132=reg131+reg132; reg93=2*reg93; reg121=reg121-reg124;
    reg131=reg60*reg143; T reg156=reg78*reg101; T reg157=reg78*reg91; T reg158=reg73*reg136; T reg159=reg60*reg148;
    T reg160=reg81*reg91; T reg161=reg4*reg149; T reg162=reg50*reg61; T reg163=reg4*reg141; T reg164=reg78*reg136;
    T reg165=reg60*reg113; T reg166=reg82*reg101; T reg167=reg58*reg148; T reg168=reg150-reg151; T reg169=reg155-reg146;
    T reg170=reg152-reg153; T reg171=reg154-reg130; T reg172=reg25*reg145; T reg173=reg79*reg85; T reg174=reg73*reg83;
    T reg175=reg50*reg133; T reg176=reg60*reg115; T reg177=reg25*reg149; T reg178=reg78*reg85; reg148=reg39*reg148;
    T reg179=reg80*reg101; T reg180=reg79*reg91; T reg181=reg60*reg132; T reg182=reg39*reg143; T reg183=reg80*reg91;
    T reg184=reg60*reg120; T reg185=reg78*reg95; T reg186=reg50*reg141; T reg187=reg39*reg115; T reg188=reg80*reg85;
    T reg189=reg73*reg95; T reg190=reg50*reg145; T reg191=reg79*reg101; reg141=reg25*reg141; T reg192=reg73*reg85;
    T reg193=reg78*reg93; T reg194=reg25*reg121; T reg195=reg79*reg95; T reg196=reg73*reg101; T reg197=reg50*reg121;
    reg149=reg50*reg149; reg101=reg81*reg101; T reg198=reg73*reg91; T reg199=reg73*reg93; T reg200=reg50*reg126;
    reg163=reg101-reg163; reg184=reg185+reg184; reg176=reg178+reg176; reg101=reg162+reg158; reg161=reg160-reg161;
    reg160=reg164+reg165; reg193=reg181+reg193; reg169=reg36*reg169; reg180=reg177-reg180; reg171=reg36*reg171;
    reg170=reg36*reg170; reg173=reg172-reg173; reg168=reg36*reg168; reg166=reg167-reg166; reg195=reg194-reg195;
    reg192=reg190+reg192; reg197=reg189+reg197; reg174=reg175+reg174; reg187=reg188-reg187; reg198=reg149+reg198;
    reg148=reg179-reg148; reg159=reg156+reg159; reg199=reg200+reg199; reg191=reg141-reg191; reg182=reg183-reg182;
    reg196=reg186+reg196; reg131=reg157+reg131; reg193=reg36*reg193; reg187=reg36*reg187; reg196=reg36*reg196;
    reg169=ponderation*reg169; reg182=reg36*reg182; reg195=reg36*reg195; reg171=ponderation*reg171; reg148=reg36*reg148;
    reg166=reg36*reg166; reg173=reg36*reg173; reg170=ponderation*reg170; reg168=ponderation*reg168; reg191=reg36*reg191;
    reg176=reg36*reg176; reg159=reg36*reg159; reg161=reg36*reg161; reg197=reg36*reg197; reg163=reg36*reg163;
    reg141=reg36*reg101; reg199=reg36*reg199; reg131=reg36*reg131; reg180=reg36*reg180; reg184=reg36*reg184;
    reg198=reg36*reg198; reg192=reg36*reg192; reg149=reg36*reg160; reg174=reg36*reg174; matrix(indices[1]+0,indices[1]+0)+=ponderation*reg195;
    matrix(indices[0]+0,indices[0]+1)+=ponderation*reg199; matrix(indices[2]+1,indices[2]+1)+=ponderation*reg166; matrix(indices[0]+0,indices[2]+1)+=ponderation*reg196; matrix(indices[0]+0,indices[1]+1)+=ponderation*reg192; matrix(indices[1]+0,indices[2]+0)+=ponderation*reg180;
    matrix(indices[0]+0,indices[1]+0)+=ponderation*reg197; matrix(indices[1]+0,indices[2]+1)+=ponderation*reg191; matrix(indices[0]+1,indices[2]+1)+=ponderation*reg159; matrix(indices[0]+1,indices[1]+1)+=ponderation*reg176; matrix(indices[2]+0,indices[2]+0)+=ponderation*reg161;
    matrix(indices[2]+0,indices[2]+1)+=ponderation*reg163; matrix(indices[0]+1,indices[2]+0)+=ponderation*reg131; reg131=ponderation*reg141; sollicitation[indices[0]+0]+=reg131; matrix(indices[0]+1,indices[1]+0)+=ponderation*reg184;
    matrix(indices[0]+0,indices[2]+0)+=ponderation*reg198; reg156=ponderation*reg149; sollicitation[indices[0]+1]+=reg156; matrix(indices[0]+1,indices[0]+1)+=ponderation*reg193; matrix(indices[1]+1,indices[1]+1)+=ponderation*reg187;
    sollicitation[indices[1]+0]+=-reg169; matrix(indices[1]+1,indices[2]+0)+=ponderation*reg182; sollicitation[indices[1]+1]+=-reg171; matrix(indices[1]+0,indices[1]+1)+=ponderation*reg173; matrix(indices[1]+1,indices[2]+1)+=ponderation*reg148;
    sollicitation[indices[2]+0]+=-reg170; matrix(indices[0]+0,indices[0]+0)+=ponderation*reg174; sollicitation[indices[2]+1]+=-reg168;
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
    T reg0=pow((*f.m).v1[1],2); T reg1=pow((*f.m).v1[0],2); reg0=reg1+reg0; reg1=pow((*f.m).v1[2],2); T reg2=pow((*f.m).v2[0],2);
    T reg3=pow((*f.m).v2[1],2); T reg4=pow((*f.m).v2[2],2); reg3=reg2+reg3; reg2=2*(*f.m).shear_modulus_13; reg0=reg1+reg0;
    reg1=2*(*f.m).shear_modulus_23; reg4=reg3+reg4; reg0=pow(reg0,0.5); reg3=2*(*f.m).shear_modulus_12; reg1=1.0/reg1;
    reg2=1.0/reg2; reg3=1.0/reg3; T reg5=(*f.m).v1[1]/reg0; T reg6=(*f.m).v1[0]/reg0; T reg7=reg2*reg1;
    reg4=pow(reg4,0.5); T reg8=(*f.m).v2[1]/reg4; T reg9=(*f.m).v2[0]/reg4; T reg10=2*reg5; T reg11=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg12=1.0/(*f.m).elastic_modulus_3; T reg13=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; reg0=(*f.m).v1[2]/reg0; T reg14=2*reg6; T reg15=reg3*reg7;
    T reg16=reg13*reg15; T reg17=2*reg0; T reg18=reg10*reg8; T reg19=reg11*reg15; T reg20=reg12*reg15;
    T reg21=reg9*reg14; T reg22=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg23=1.0/(*f.m).elastic_modulus_2; T reg24=1.0/(*f.m).elastic_modulus_1; reg4=(*f.m).v2[2]/reg4;
    T reg25=pow(reg9,2); T reg26=pow(reg8,2); T reg27=reg18*reg22; T reg28=pow(reg4,2); reg17=reg17*reg4;
    T reg29=reg25*reg22; T reg30=reg26*reg23; T reg31=reg21*reg22; T reg32=reg18*reg23; T reg33=reg23*reg20;
    T reg34=reg25*reg24; T reg35=reg26*reg22; T reg36=reg11*reg19; T reg37=reg21*reg24; T reg38=reg11*reg16;
    T reg39=reg22*reg20; T reg40=reg13*reg17; reg27=reg37-reg27; reg37=pow(reg6,2); T reg41=pow(reg5,2);
    T reg42=reg28*reg11; T reg43=reg23*reg16; reg29=reg30-reg29; reg30=reg22*reg19; T reg44=reg17*reg11;
    reg38=reg39+reg38; reg36=reg33-reg36; reg31=reg32-reg31; reg32=reg25*reg13; reg33=reg26*reg11;
    T reg45=reg28*reg13; T reg46=reg13*reg21; reg35=reg34-reg35; reg34=reg18*reg11; T reg47=reg24*reg37;
    reg34=reg46+reg34; reg17=reg17*reg12; reg40=reg27-reg40; reg27=reg22*reg37; reg46=reg23*reg41;
    T reg48=reg5*reg9; T reg49=reg6*reg8; reg45=reg35-reg45; reg35=reg30+reg43; reg33=reg32+reg33;
    reg42=reg29-reg42; reg29=reg28*reg12; reg32=pow(reg0,2); T reg50=reg24*reg36; reg44=reg31-reg44;
    reg31=reg22*reg41; T reg51=reg22*reg38; T reg52=reg5*reg8; T reg53=reg6*reg9; T reg54=reg53*reg45;
    T reg55=reg52*reg42; reg51=reg50-reg51; reg50=reg13*reg35; T reg56=reg23*reg15; T reg57=reg13*reg19;
    reg33=reg29-reg33; reg20=reg24*reg20; reg29=reg13*reg16; reg15=reg22*reg15; T reg58=reg12*reg7;
    T reg59=reg0*reg4; T reg60=reg48+reg49; T reg61=reg13*reg7; T reg62=reg41*reg44; reg34=reg17-reg34;
    reg17=reg37*reg40; T reg63=reg53*reg40; T reg64=reg52*reg44; T reg65=reg6*reg4; T reg66=reg0*reg9;
    T reg67=reg0*reg8; T reg68=reg5*reg4; T reg69=reg41*reg42; T reg70=reg3*reg1; T reg71=reg37*reg45;
    reg7=reg11*reg7; reg27=reg46-reg27; reg46=reg11*reg32; T reg72=reg25*reg40; T reg73=2*reg9;
    T reg74=reg26*reg44; T reg75=reg26*reg42; T reg76=reg9*reg8; reg31=reg47-reg31; reg47=reg13*reg37;
    T reg77=reg25*reg45; T reg78=reg11*reg41; T reg79=reg13*reg32; reg74=reg72+reg74; reg72=reg76*reg3;
    T reg80=reg60*reg3; reg46=reg27-reg46; reg27=reg28*reg34; reg64=reg63+reg64; reg63=reg59*reg34;
    reg62=reg17+reg62; reg17=reg11*reg61; T reg81=reg11*reg7; reg16=reg22*reg16; T reg82=reg9*reg4;
    T reg83=reg22*reg58; reg79=reg31-reg79; reg31=reg68-reg67; T reg84=reg32*reg33; reg69=reg71+reg69;
    reg71=reg5*reg14; T reg85=reg11*reg70; T reg86=reg13*reg70; T reg87=reg3*reg2; T reg88=reg12*reg32;
    reg78=reg47+reg78; reg55=reg54+reg55; reg47=reg59*reg33; reg75=reg77+reg75; reg54=reg28*reg33;
    reg57=reg39+reg57; reg39=reg13*reg56; reg50=reg51-reg50; reg51=reg65+reg66; reg77=reg32*reg34;
    reg70=reg12*reg70; T reg89=reg73*reg8; T reg90=reg13*reg15; reg58=reg23*reg58; reg19=reg24*reg19;
    reg29=reg20-reg29; reg20=reg89*reg72; T reg91=reg26*reg46; T reg92=reg13*reg87; reg12=reg12*reg87;
    reg54=reg75+reg54; reg75=reg22*reg70; T reg93=reg89*reg80; T reg94=reg25*reg79; T reg95=reg11*reg86;
    reg84=reg69+reg84; reg69=reg71*reg80; reg70=reg23*reg70; reg27=reg74+reg27; reg17=reg83+reg17;
    reg77=reg62+reg77; reg81=reg58-reg81; reg58=reg71*reg72; reg62=reg11*reg85; reg15=reg22*reg15;
    reg16=reg16+reg19; reg90=reg19+reg90; reg29=reg29/reg50; reg39=reg30+reg39; reg57=reg57/reg50;
    reg38=reg38/reg50; reg19=reg60*reg72; reg47=reg55+reg47; reg36=reg36/reg50; reg78=reg88-reg78;
    reg56=reg24*reg56; reg55=2*reg31; reg74=2*reg8; reg83=reg73*reg4; reg88=reg6*reg5;
    T reg96=reg0*reg14; T reg97=reg8*reg4; T reg98=reg41*reg46; reg61=reg23*reg61; T reg99=reg37*reg79;
    reg87=reg11*reg87; reg65=reg66-reg65; reg66=reg51*reg2; T reg100=reg60*reg80; reg63=reg64+reg63;
    reg64=reg82*reg2; reg7=reg22*reg7; reg68=reg67+reg68; reg67=reg26*reg29; T reg101=reg57*reg71;
    T reg102=reg57*reg41; T reg103=reg89*reg29; T reg104=reg83*reg66; reg93=reg27+reg93; reg27=reg25*reg29;
    T reg105=reg32*reg78; T reg106=reg57*reg37; T reg107=reg10*reg0; reg98=reg99+reg98; reg99=reg83*reg64;
    reg20=reg54+reg20; reg54=reg74*reg4; T reg108=pow(reg31,2); T reg109=pow(reg65,2); T reg110=reg55*reg65;
    T reg111=reg68*reg1; T reg112=reg28*reg78; reg91=reg94+reg91; reg94=reg97*reg1; T reg113=reg96*reg66;
    reg69=reg77+reg69; reg58=reg84+reg58; reg77=reg96*reg64; reg3=reg88*reg3; reg84=reg25*reg38;
    reg39=reg39/reg50; reg86=reg23*reg86; T reg114=reg36*reg37; reg90=reg90/reg50; reg16=reg16/reg50;
    T reg115=reg23*reg12; reg15=reg56-reg15; reg12=reg22*reg12; reg56=reg36*reg41; reg19=reg47+reg19;
    reg47=reg51*reg64; T reg116=reg11*reg87; reg85=reg22*reg85; reg100=reg63+reg100; reg11=reg11*reg92;
    reg63=reg51*reg66; reg35=reg35/reg50; T reg117=reg7+reg61; T reg118=reg52*reg46; T reg119=reg89*reg38;
    T reg120=reg36*reg71; reg62=reg70-reg62; reg70=reg53*reg79; reg17=reg22*reg17; T reg121=reg26*reg38;
    reg81=reg24*reg81; reg95=reg75+reg95; reg75=reg6*reg0; reg117=reg13*reg117; reg92=reg23*reg92;
    reg23=reg54*reg111; reg116=reg115-reg116; reg115=reg85+reg86; reg63=reg100+reg63; reg100=reg68*reg111;
    reg112=reg91+reg112; reg62=reg24*reg62; reg91=reg89*reg3; T reg122=reg107*reg94; T reg123=reg108*reg35;
    reg77=reg58+reg77; reg113=reg69+reg113; reg114=reg84+reg114; reg17=reg81-reg17; reg58=reg107*reg111;
    reg105=reg98+reg105; reg69=reg71*reg3; reg95=reg22*reg95; reg87=reg22*reg87; reg119=reg120+reg119;
    reg81=reg110*reg35; reg27=reg106+reg27; reg84=reg108*reg16; reg98=reg109*reg35; reg67=reg102+reg67;
    reg102=reg109*reg16; reg106=reg5*reg0; reg121=reg56+reg121; reg103=reg101+reg103; reg56=reg110*reg16;
    reg101=reg39*reg37; reg120=reg25*reg90; reg118=reg70+reg118; reg70=reg59*reg78; T reg124=reg39*reg41;
    reg99=reg20+reg99; reg20=reg54*reg94; reg11=reg12+reg11; reg2=reg75*reg2; reg12=reg68*reg94;
    reg47=reg19+reg47; reg104=reg93+reg104; reg15=reg15/reg50; reg19=reg89*reg90; reg93=reg39*reg71;
    T reg125=reg26*reg90; T reg126=reg60*reg3; reg69=reg105+reg69; reg105=reg96*reg2; reg12=reg47+reg12;
    reg100=reg63+reg100; reg23=reg104+reg23; reg1=reg106*reg1; reg122=reg77+reg122; reg70=reg118+reg70;
    reg47=reg87+reg92; reg123=reg114+reg123; reg11=reg22*reg11; reg20=reg99+reg20; reg84=reg27+reg84;
    reg98=reg121+reg98; reg22=reg110*reg15; reg116=reg24*reg116; reg102=reg67+reg102; reg56=reg103+reg56;
    reg19=reg93+reg19; reg24=reg83*reg2; reg27=reg109*reg15; reg81=reg119+reg81; reg117=reg17-reg117;
    reg120=reg101+reg120; reg95=reg62-reg95; reg17=reg108*reg15; reg58=reg113+reg58; reg115=reg13*reg115;
    reg91=reg112+reg91; reg125=reg124+reg125; reg62=reg51*reg2; reg63=reg20*reg100; reg126=reg70+reg126;
    reg67=reg76*reg84; reg22=reg19+reg22; reg27=reg125+reg27; reg17=reg120+reg17; reg19=reg58*reg12;
    reg70=reg65*reg31; reg77=reg88*reg98; reg93=reg5*reg31; reg99=reg76*reg102; reg101=reg88*reg123;
    reg117=reg117/reg50; reg103=reg122*reg100; reg104=reg76*reg56; reg112=reg88*reg81; reg113=reg23*reg12;
    reg114=reg107*reg1; reg105=reg69+reg105; reg24=reg91+reg24; reg69=reg54*reg1; reg115=reg95-reg115;
    reg47=reg13*reg47; reg13=reg6*reg65; reg11=reg116-reg11; reg91=reg70*reg22; reg104=reg112+reg104;
    reg99=reg77+reg99; reg77=reg122*reg23; reg95=reg70*reg27; reg19=reg103-reg19; reg67=reg101+reg67;
    reg101=reg70*reg17; reg103=reg58*reg20; reg112=reg41*reg98; reg116=reg26*reg102; reg115=reg115/reg50;
    reg118=reg26*reg84; reg119=reg41*reg123; reg120=reg25*reg56; reg121=reg41*reg81; reg56=reg26*reg56;
    reg81=reg37*reg81; reg114=reg105+reg114; reg69=reg24+reg69; reg13=reg93+reg13; reg102=reg25*reg102;
    reg98=reg37*reg98; reg84=reg25*reg84; reg123=reg37*reg123; reg47=reg11-reg47; reg11=reg5*reg65;
    reg24=reg6*reg31; reg93=reg60*reg117; reg105=reg52*reg117; reg124=reg53*reg117; reg125=reg68*reg1;
    reg62=reg126+reg62; reg126=reg8*reg31; T reg127=reg9*reg65; reg113=reg63-reg113; reg125=reg62+reg125;
    reg103=reg77-reg103; reg62=reg8*reg65; reg63=elem.pos(1)[0]-elem.pos(0)[0]; reg77=elem.pos(1)[1]-elem.pos(0)[1]; T reg128=reg60*reg93;
    reg91=reg104+reg91; reg104=elem.pos(2)[0]-elem.pos(0)[0]; T reg129=elem.pos(2)[1]-elem.pos(0)[1]; T reg130=reg108*reg17; reg84=reg123+reg84;
    reg123=reg60*reg105; reg95=reg99+reg95; reg127=reg126+reg127; reg99=reg13*reg115; reg126=reg11*reg115;
    T reg131=reg60*reg124; reg101=reg67+reg101; reg67=reg24*reg115; reg50=reg47/reg50; reg47=reg69*reg19;
    T reg132=reg113*reg114; T reg133=reg109*reg22; reg56=reg121+reg56; reg121=reg109*reg27; reg116=reg112+reg116;
    reg17=reg109*reg17; reg118=reg119+reg118; reg22=reg108*reg22; reg112=reg9*reg31; reg102=reg98+reg102;
    reg27=reg108*reg27; reg120=reg81+reg120; reg81=reg21*reg93; reg130=reg84+reg130; reg84=reg127*reg50;
    reg98=reg21*reg124; reg17=reg118+reg17; reg118=reg62*reg50; reg119=reg112*reg50; reg124=reg18*reg124;
    reg27=reg102+reg27; reg22=reg120+reg22; reg121=reg116+reg121; reg102=reg18*reg105; reg105=reg21*reg105;
    reg93=reg18*reg93; reg133=reg56+reg133; reg56=reg122*reg125; reg116=reg125*reg103; reg47=reg132-reg47;
    reg120=reg69*reg12; reg132=reg114*reg12; T reg134=reg13*reg99; reg128=reg91+reg128; reg91=reg63*reg129;
    T reg135=reg13*reg126; reg123=reg95+reg123; reg95=reg77*reg104; T reg136=reg114*reg20; T reg137=reg13*reg67;
    reg131=reg101+reg131; reg14=reg31*reg14; reg101=reg122*reg69; reg10=reg10*reg65; T reg138=reg20*reg125;
    reg138=reg120-reg138; reg124=reg17+reg124; reg17=reg10*reg67; reg120=reg58*reg125; reg56=reg132-reg56;
    reg132=reg114*reg23; reg102=reg121+reg102; reg121=reg10*reg126; T reg139=reg58*reg69; reg101=reg136-reg101;
    reg136=reg114*reg100; T reg140=reg14*reg99; reg81=reg22+reg81; reg22=reg23*reg125; T reg141=reg69*reg100;
    reg126=reg14*reg126; reg105=reg27+reg105; reg116=reg47+reg116; reg67=reg14*reg67; reg98=reg130+reg98;
    reg95=reg91-reg95; reg73=reg73*reg31; reg74=reg74*reg65; reg137=reg131+reg137; reg27=reg127*reg119;
    reg135=reg123+reg135; reg47=reg127*reg118; reg91=reg127*reg84; reg134=reg128+reg134; reg93=reg133+reg93;
    reg99=reg10*reg99; reg123=1-(*f.m).resolution; reg140=reg81+reg140; reg99=reg93+reg99; reg101=reg101/reg116;
    reg81=reg73*reg118; reg126=reg105+reg126; reg139=reg132-reg139; reg91=reg134+reg91; reg56=reg56/reg116;
    reg120=reg136-reg120; reg93=reg73*reg119; reg67=reg98+reg67; reg138=reg138/reg116; reg22=reg141-reg22;
    reg77=reg77/reg95; reg104=reg104/reg95; reg47=reg135+reg47; reg121=reg102+reg121; reg118=reg74*reg118;
    reg129=reg129/reg95; reg27=reg137+reg27; reg63=reg63/reg95; reg98=reg74*reg84; reg84=reg73*reg84;
    reg119=reg74*reg119; reg17=reg124+reg17; reg22=reg22/reg116; reg113=reg113/reg116; reg19=reg19/reg116;
    reg120=reg120/reg116; reg93=reg67+reg93; reg84=reg140+reg84; reg103=reg103/reg116; reg116=reg139/reg116;
    reg67=reg91*reg123; reg102=reg47*reg123; reg105=reg27*reg123; reg124=reg104-reg63; reg128=(*f.m).resolution*reg101;
    reg130=(*f.m).resolution*reg56; reg131=(*f.m).resolution*reg138; reg132=reg77-reg129; reg118=reg121+reg118; reg98=reg99+reg98;
    reg119=reg17+reg119; reg81=reg126+reg81; reg67=reg128+reg67; reg130=reg102-reg130; reg17=0.5*reg104;
    reg99=0.5*reg129; reg102=0.5*reg63; reg105=reg131+reg105; reg121=0.5*reg132; reg126=0.5*reg124;
    reg128=0.5*reg77; reg131=reg98*reg123; reg133=reg118*reg123; reg134=reg119*reg123; reg135=(*f.m).resolution*reg19;
    reg136=(*f.m).resolution*reg113; reg137=(*f.m).resolution*reg103; reg139=reg84*reg123; reg140=(*f.m).resolution*reg22; reg141=(*f.m).resolution*reg120;
    T reg142=(*f.m).resolution*reg116; T reg143=reg81*reg123; T reg144=reg93*reg123; T reg145=reg124*reg130; T reg146=reg126*reg67;
    T reg147=reg132*reg105; T reg148=reg121*reg67; T reg149=reg129*reg105; T reg150=reg17*reg67; T reg151=reg128*reg67;
    reg142=reg131-reg142; reg133=reg141+reg133; reg140=reg134-reg140; reg139=reg137+reg139; reg135=reg143-reg135;
    reg144=reg136+reg144; reg131=reg99*reg67; reg134=reg104*reg130; reg136=reg63*reg130; reg137=reg77*reg105;
    reg141=reg102*reg67; reg143=reg99*reg139; T reg152=reg104*reg135; T reg153=reg129*reg140; T reg154=reg121*reg142;
    T reg155=reg124*reg133; T reg156=reg77*reg140; T reg157=reg126*reg142; T reg158=reg132*reg140; T reg159=reg102*reg142;
    reg136=reg136-reg151; reg149=reg149-reg150; T reg160=reg17*reg139; T reg161=reg129*reg144; T reg162=reg128*reg142;
    T reg163=reg63*reg133; reg131=reg131-reg134; reg148=reg145+reg148; reg145=reg17*reg142; T reg164=reg128*reg139;
    T reg165=reg63*reg135; T reg166=reg104*reg133; T reg167=reg99*reg142; T reg168=reg132*reg144; T reg169=reg126*reg139;
    reg141=reg141-reg137; T reg170=reg121*reg139; T reg171=reg124*reg135; T reg172=reg77*reg144; T reg173=reg102*reg139;
    reg146=reg147+reg146; reg167=reg167-reg166; reg149=2*reg149; reg157=reg158+reg157; reg161=reg161-reg160;
    reg173=reg173-reg172; reg136=2*reg136; reg163=reg163-reg162; reg148=2*reg148; reg159=reg159-reg156;
    reg143=reg143-reg152; reg141=2*reg141; reg165=reg165-reg164; reg153=reg153-reg145; reg146=2*reg146;
    reg131=2*reg131; reg169=reg168+reg169; reg170=reg171+reg170; reg154=reg155+reg154; reg147=reg104*reg154;
    reg155=reg99*reg148; reg158=reg126*reg149; reg168=reg126*reg146; reg171=reg63*reg159; T reg174=reg124*reg154;
    T reg175=reg121*reg148; T reg176=reg104*reg157; T reg177=reg99*reg146; T reg178=reg132*reg143; T reg179=reg128*reg131;
    T reg180=reg128*reg149; T reg181=reg17*reg136; T reg182=reg129*reg165; T reg183=reg132*reg170; T reg184=reg63*reg167;
    T reg185=reg128*reg136; T reg186=reg132*reg165; T reg187=reg124*reg153; T reg188=reg63*reg163; T reg189=reg102*reg148;
    T reg190=reg77*reg170; T reg191=reg129*reg169; T reg192=reg121*reg131; T reg193=reg126*reg141; T reg194=reg128*reg141;
    T reg195=reg102*reg149; T reg196=reg129*reg161; T reg197=reg132*reg169; reg169=reg77*reg169; T reg198=reg102*reg146;
    T reg199=reg104*reg153; T reg200=reg99*reg149; T reg201=reg132*reg173; T reg202=reg104*reg159; T reg203=reg99*reg141;
    T reg204=reg128*reg148; T reg205=reg124*reg163; T reg206=reg121*reg149; T reg207=reg126*reg148; reg165=reg77*reg165;
    reg148=reg17*reg148; T reg208=reg77*reg161; T reg209=reg102*reg136; T reg210=reg124*reg167; T reg211=reg77*reg173;
    T reg212=reg102*reg131; T reg213=reg77*reg143; reg143=reg129*reg143; T reg214=reg17*reg131; T reg215=reg102*reg141;
    T reg216=reg99*reg131; reg167=reg104*reg167; T reg217=reg121*reg141; reg153=reg63*reg153; T reg218=reg124*reg157;
    reg173=reg129*reg173; reg141=reg17*reg141; T reg219=reg17*reg146; T reg220=reg121*reg146; reg154=reg63*reg154;
    reg163=reg104*reg163; reg159=reg124*reg159; T reg221=reg126*reg136; T reg222=reg99*reg136; reg131=reg126*reg131;
    reg149=reg17*reg149; reg146=reg128*reg146; reg157=reg63*reg157; reg170=reg129*reg170; reg136=reg121*reg136;
    reg161=reg132*reg161; reg131=reg178+reg131; reg193=reg201+reg193; reg161=reg158+reg161; reg159=reg217+reg159;
    reg205=reg136+reg205; reg208=reg195-reg208; reg213=reg212-reg213; reg211=reg215-reg211; reg165=reg209-reg165;
    reg146=reg157-reg146; reg204=reg154-reg204; reg180=reg153-reg180; reg179=reg184-reg179; reg194=reg171-reg194;
    reg185=reg188-reg185; reg221=reg186+reg221; reg149=reg196-reg149; reg218=reg220+reg218; reg207=reg183+reg207;
    reg181=reg182-reg181; reg168=reg197+reg168; reg176=reg177-reg176; reg210=reg192+reg210; reg199=reg200-reg199;
    reg214=reg143-reg214; reg167=reg216-reg167; reg147=reg155-reg147; reg190=reg189-reg190; reg175=reg174+reg175;
    reg187=reg206+reg187; reg141=reg173-reg141; reg148=reg170-reg148; reg202=reg203-reg202; reg169=reg198-reg169;
    reg163=reg222-reg163; reg219=reg191-reg219; reg146=reg95*reg146; reg202=reg95*reg202; reg159=reg95*reg159;
    reg149=reg95*reg149; reg207=reg95*reg207; reg205=reg95*reg205; reg180=reg95*reg180; reg214=reg95*reg214;
    reg165=reg95*reg165; reg208=reg95*reg208; reg163=reg95*reg163; reg210=reg95*reg210; reg213=reg95*reg213;
    reg167=reg95*reg167; reg211=reg95*reg211; reg204=reg95*reg204; reg131=reg95*reg131; reg175=reg95*reg175;
    reg219=reg95*reg219; reg185=reg95*reg185; reg148=reg95*reg148; reg169=reg95*reg169; reg187=reg95*reg187;
    reg193=reg95*reg193; reg194=reg95*reg194; reg199=reg95*reg199; reg190=reg95*reg190; reg147=reg95*reg147;
    reg176=reg95*reg176; reg181=reg95*reg181; reg221=reg95*reg221; reg179=reg95*reg179; reg218=reg95*reg218;
    reg161=reg95*reg161; reg141=reg95*reg141; reg168=reg95*reg168; matrix(indices[2]+1,indices[2]+1)+=ponderation*reg185; matrix(indices[2]+1,indices[0]+0)+=ponderation*reg146;
    matrix(indices[0]+0,indices[0]+0)+=ponderation*reg168; matrix(indices[1]+1,indices[2]+1)+=ponderation*reg163; matrix(indices[2]+1,indices[1]+0)+=ponderation*reg180; matrix(indices[2]+0,indices[0]+0)+=ponderation*reg169; matrix(indices[2]+1,indices[2]+0)+=ponderation*reg194;
    matrix(indices[2]+1,indices[0]+1)+=ponderation*reg204; matrix(indices[2]+1,indices[1]+1)+=ponderation*reg179; matrix(indices[1]+0,indices[1]+1)+=ponderation*reg214; matrix(indices[2]+0,indices[0]+1)+=ponderation*reg190; matrix(indices[0]+0,indices[1]+1)+=ponderation*reg131;
    matrix(indices[1]+0,indices[0]+0)+=ponderation*reg219; matrix(indices[0]+1,indices[0]+1)+=ponderation*reg175; matrix(indices[1]+0,indices[0]+1)+=ponderation*reg148; matrix(indices[0]+0,indices[2]+0)+=ponderation*reg193; matrix(indices[0]+1,indices[1]+0)+=ponderation*reg187;
    matrix(indices[1]+1,indices[1]+0)+=ponderation*reg199; matrix(indices[1]+1,indices[0]+1)+=ponderation*reg147; matrix(indices[1]+1,indices[0]+0)+=ponderation*reg176; matrix(indices[1]+0,indices[2]+1)+=ponderation*reg181; matrix(indices[1]+0,indices[2]+0)+=ponderation*reg141;
    matrix(indices[0]+0,indices[1]+0)+=ponderation*reg161; matrix(indices[0]+1,indices[0]+0)+=ponderation*reg218; matrix(indices[0]+0,indices[2]+1)+=ponderation*reg221; matrix(indices[0]+1,indices[2]+0)+=ponderation*reg159; matrix(indices[0]+0,indices[0]+1)+=ponderation*reg207;
    matrix(indices[0]+1,indices[2]+1)+=ponderation*reg205; matrix(indices[1]+0,indices[1]+0)+=ponderation*reg149; matrix(indices[2]+0,indices[1]+0)+=ponderation*reg208; matrix(indices[0]+1,indices[1]+1)+=ponderation*reg210; matrix(indices[2]+0,indices[1]+1)+=ponderation*reg213;
    matrix(indices[2]+0,indices[2]+0)+=ponderation*reg211; matrix(indices[1]+1,indices[1]+1)+=ponderation*reg167; matrix(indices[2]+0,indices[2]+1)+=ponderation*reg165; matrix(indices[1]+1,indices[2]+0)+=ponderation*reg202;
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
    T reg0=pow((*f.m).v1[1],2); T reg1=pow((*f.m).v1[0],2); reg0=reg1+reg0; reg1=pow((*f.m).v1[2],2); T reg2=pow((*f.m).v2[0],2);
    T reg3=pow((*f.m).v2[1],2); T reg4=2*(*f.m).shear_modulus_13; T reg5=2*(*f.m).shear_modulus_23; reg3=reg2+reg3; reg2=pow((*f.m).v2[2],2);
    reg0=reg1+reg0; reg5=1.0/reg5; reg1=2*(*f.m).shear_modulus_12; reg4=1.0/reg4; reg0=pow(reg0,0.5);
    reg2=reg3+reg2; reg3=(*f.m).v1[1]/reg0; T reg6=(*f.m).v1[0]/reg0; T reg7=reg4*reg5; reg1=1.0/reg1;
    reg2=pow(reg2,0.5); T reg8=2*reg6; T reg9=(*f.m).v2[1]/reg2; T reg10=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg11=1.0/(*f.m).elastic_modulus_3;
    T reg12=(*f.m).v2[0]/reg2; T reg13=2*reg3; reg0=(*f.m).v1[2]/reg0; T reg14=reg1*reg7; T reg15=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    reg2=(*f.m).v2[2]/reg2; T reg16=pow(reg12,2); T reg17=pow(reg9,2); T reg18=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg19=1.0/(*f.m).elastic_modulus_2;
    T reg20=1.0/(*f.m).elastic_modulus_1; T reg21=reg11*reg14; T reg22=reg13*reg9; T reg23=reg10*reg14; T reg24=reg12*reg8;
    T reg25=reg15*reg14; T reg26=2*reg0; T reg27=reg24*reg18; T reg28=reg22*reg19; T reg29=reg18*reg21;
    T reg30=reg10*reg25; T reg31=pow(reg2,2); T reg32=reg24*reg20; T reg33=reg17*reg18; T reg34=reg17*reg19;
    reg26=reg26*reg2; T reg35=reg16*reg18; T reg36=reg22*reg18; T reg37=reg10*reg23; T reg38=reg16*reg20;
    T reg39=reg19*reg21; T reg40=pow(reg3,2); T reg41=pow(reg6,2); T reg42=reg17*reg10; T reg43=reg15*reg24;
    T reg44=reg19*reg25; T reg45=reg18*reg23; reg30=reg29+reg30; T reg46=reg16*reg15; T reg47=reg22*reg10;
    reg37=reg39-reg37; reg27=reg28-reg27; reg36=reg32-reg36; reg28=reg15*reg26; reg32=reg31*reg10;
    reg35=reg34-reg35; reg34=reg26*reg10; reg33=reg38-reg33; reg38=reg31*reg15; reg42=reg46+reg42;
    reg39=reg31*reg11; reg38=reg33-reg38; reg33=reg20*reg37; reg46=reg3*reg12; T reg48=reg3*reg9;
    T reg49=reg6*reg12; T reg50=reg6*reg9; reg28=reg36-reg28; reg36=reg18*reg40; T reg51=reg18*reg41;
    T reg52=reg18*reg30; T reg53=reg19*reg40; reg47=reg43+reg47; reg43=pow(reg0,2); reg32=reg35-reg32;
    reg34=reg27-reg34; reg27=reg45+reg44; reg26=reg26*reg11; reg35=reg20*reg41; T reg54=reg16*reg38;
    T reg55=reg17*reg32; reg21=reg20*reg21; T reg56=reg46+reg50; T reg57=reg0*reg2; T reg58=reg15*reg23;
    T reg59=reg19*reg14; reg36=reg35-reg36; reg35=reg48*reg34; T reg60=reg15*reg27; T reg61=reg15*reg43;
    T reg62=reg15*reg25; T reg63=reg48*reg32; reg52=reg33-reg52; reg33=reg49*reg38; T reg64=reg40*reg34;
    T reg65=reg41*reg28; T reg66=reg15*reg41; T reg67=reg15*reg7; T reg68=reg10*reg40; T reg69=reg12*reg9;
    T reg70=reg1*reg5; T reg71=reg11*reg7; reg51=reg53-reg51; reg42=reg39-reg42; reg39=2*reg12;
    reg53=reg49*reg28; T reg72=reg10*reg43; T reg73=reg40*reg32; T reg74=reg41*reg38; T reg75=reg3*reg2;
    T reg76=reg0*reg9; reg14=reg18*reg14; T reg77=reg17*reg34; T reg78=reg16*reg28; T reg79=reg0*reg12;
    reg7=reg10*reg7; T reg80=reg6*reg2; reg47=reg26-reg47; reg26=reg39*reg9; T reg81=reg10*reg67;
    T reg82=reg56*reg1; reg55=reg54+reg55; reg54=reg31*reg42; reg35=reg53+reg35; reg53=reg10*reg7;
    T reg83=reg75-reg76; T reg84=reg57*reg47; T reg85=reg3*reg8; T reg86=reg10*reg70; T reg87=reg18*reg71;
    T reg88=reg15*reg70; T reg89=reg1*reg4; reg64=reg65+reg64; reg65=reg43*reg47; reg72=reg51-reg72;
    reg51=reg11*reg43; reg68=reg66+reg68; reg66=reg12*reg2; T reg90=reg43*reg42; reg73=reg74+reg73;
    reg74=reg31*reg47; reg77=reg78+reg77; reg25=reg18*reg25; reg60=reg52-reg60; reg61=reg36-reg61;
    reg36=reg80+reg79; reg63=reg33+reg63; reg33=reg57*reg42; reg58=reg29+reg58; reg29=reg15*reg59;
    reg52=reg69*reg1; reg62=reg21-reg62; reg70=reg11*reg70; reg23=reg20*reg23; reg71=reg19*reg71;
    reg21=reg15*reg14; reg78=reg15*reg89; T reg91=reg41*reg61; T reg92=reg36*reg4; reg68=reg51-reg68;
    reg81=reg87+reg81; reg51=reg19*reg70; reg87=reg66*reg4; T reg93=reg56*reg82; reg84=reg35+reg84;
    reg53=reg71-reg53; reg33=reg63+reg33; reg35=reg56*reg52; reg11=reg11*reg89; reg70=reg18*reg70;
    reg63=reg10*reg86; reg71=reg10*reg88; T reg94=reg39*reg2; reg62=reg62/reg60; T reg95=2*reg9;
    reg29=reg45+reg29; reg58=reg58/reg60; reg30=reg30/reg60; reg37=reg37/reg60; T reg96=2*reg83;
    reg75=reg76+reg75; reg76=reg6*reg3; T reg97=reg17*reg72; T reg98=reg16*reg61; T reg99=reg85*reg52;
    reg74=reg77+reg74; reg77=reg9*reg2; reg90=reg73+reg90; reg73=reg26*reg82; reg89=reg10*reg89;
    T reg100=reg85*reg82; reg67=reg19*reg67; reg80=reg79-reg80; reg79=reg26*reg52; reg7=reg18*reg7;
    reg54=reg55+reg54; reg55=reg40*reg72; reg65=reg64+reg65; reg14=reg18*reg14; reg59=reg20*reg59;
    reg64=reg0*reg8; reg25=reg25+reg23; reg21=reg23+reg21; reg23=reg58*reg85; T reg101=reg26*reg62;
    T reg102=reg49*reg61; T reg103=reg94*reg92; T reg104=reg48*reg72; T reg105=reg58*reg40; T reg106=reg16*reg62;
    T reg107=reg17*reg62; T reg108=reg58*reg41; reg73=reg74+reg73; reg27=reg27/reg60; reg97=reg98+reg97;
    reg74=reg31*reg68; reg98=reg96*reg80; T reg109=pow(reg80,2); T reg110=pow(reg83,2); T reg111=reg95*reg2;
    T reg112=reg13*reg0; reg35=reg33+reg35; reg33=reg36*reg87; reg79=reg54+reg79; reg54=reg94*reg87;
    reg93=reg84+reg93; reg84=reg36*reg92; T reg113=reg10*reg89; T reg114=reg18*reg11; reg11=reg19*reg11;
    T reg115=reg43*reg68; reg55=reg91+reg55; reg91=reg16*reg30; reg29=reg29/reg60; reg71=reg70+reg71;
    reg63=reg51-reg63; reg1=reg76*reg1; reg51=reg37*reg41; reg21=reg21/reg60; reg25=reg25/reg60;
    reg70=reg7+reg67; reg81=reg18*reg81; reg14=reg59-reg14; reg53=reg20*reg53; reg86=reg18*reg86;
    reg59=reg77*reg5; T reg116=reg75*reg5; reg88=reg19*reg88; T reg117=reg64*reg92; reg10=reg10*reg78;
    T reg118=reg6*reg0; reg99=reg90+reg99; reg100=reg65+reg100; reg65=reg37*reg85; reg90=reg64*reg87;
    T reg119=reg26*reg30; T reg120=reg37*reg40; T reg121=reg17*reg30; reg74=reg97+reg74; reg97=reg26*reg1;
    T reg122=reg57*reg68; T reg123=reg112*reg116; reg117=reg100+reg117; reg54=reg79+reg54; reg84=reg93+reg84;
    reg79=reg75*reg116; reg93=reg75*reg59; reg33=reg35+reg33; reg4=reg118*reg4; reg35=reg112*reg59;
    reg104=reg102+reg104; reg90=reg99+reg90; reg99=reg111*reg116; reg103=reg73+reg103; reg73=reg85*reg1;
    reg115=reg55+reg115; reg55=reg111*reg59; reg100=reg86+reg88; reg102=reg29*reg85; T reg124=reg17*reg21;
    reg113=reg11-reg113; reg11=reg29*reg40; reg10=reg114+reg10; reg114=reg16*reg21; T reg125=reg29*reg41;
    T reg126=reg98*reg25; reg101=reg23+reg101; reg121=reg120+reg121; reg23=reg109*reg27; reg119=reg65+reg119;
    reg65=reg98*reg27; reg120=reg109*reg25; reg107=reg105+reg107; reg106=reg108+reg106; reg105=reg110*reg25;
    reg108=reg3*reg0; reg14=reg14/reg60; reg89=reg18*reg89; reg51=reg91+reg51; reg91=reg110*reg27;
    reg78=reg19*reg78; reg81=reg53-reg81; reg70=reg15*reg70; reg19=reg26*reg21; reg71=reg18*reg71;
    reg63=reg20*reg63; reg53=reg64*reg4; reg93=reg33+reg93; reg79=reg84+reg79; reg19=reg102+reg19;
    reg120=reg107+reg120; reg33=reg109*reg14; reg123=reg117+reg123; reg84=reg98*reg14; reg124=reg11+reg124;
    reg126=reg101+reg126; reg122=reg104+reg122; reg114=reg125+reg114; reg11=reg110*reg14; reg35=reg90+reg35;
    reg90=reg56*reg1; reg99=reg103+reg99; reg5=reg108*reg5; reg70=reg81-reg70; reg71=reg63-reg71;
    reg100=reg15*reg100; reg63=reg94*reg4; reg113=reg20*reg113; reg10=reg18*reg10; reg73=reg115+reg73;
    reg97=reg74+reg97; reg105=reg106+reg105; reg55=reg54+reg55; reg65=reg119+reg65; reg23=reg121+reg23;
    reg91=reg51+reg91; reg18=reg89+reg78; reg20=reg69*reg120; reg51=reg76*reg65; reg54=reg69*reg126;
    reg74=reg69*reg105; reg81=reg76*reg91; reg101=reg112*reg5; reg102=reg76*reg23; reg103=reg99*reg93;
    reg53=reg73+reg53; reg73=reg35*reg79; reg104=reg55*reg79; reg106=reg80*reg83; reg107=reg123*reg93;
    reg84=reg19+reg84; reg19=reg6*reg80; reg115=reg3*reg83; reg18=reg15*reg18; reg11=reg114+reg11;
    reg10=reg113-reg10; reg63=reg97+reg63; reg15=reg111*reg5; reg33=reg124+reg33; reg100=reg71-reg100;
    reg90=reg122+reg90; reg70=reg70/reg60; reg71=reg36*reg4; reg97=reg16*reg105; reg100=reg100/reg60;
    reg113=reg41*reg91; reg103=reg104-reg103; reg104=reg106*reg84; reg54=reg51+reg54; reg51=reg56*reg70;
    reg101=reg53+reg101; reg53=reg48*reg70; reg114=reg49*reg70; reg18=reg10-reg18; reg10=reg106*reg33;
    reg71=reg90+reg71; reg90=reg75*reg5; reg117=reg123*reg55; reg119=reg35*reg99; reg107=reg73-reg107;
    reg19=reg115+reg19; reg73=reg40*reg65; reg115=reg17*reg126; reg121=reg17*reg120; reg122=reg40*reg23;
    reg124=reg3*reg80; reg125=reg6*reg83; T reg127=reg12*reg80; reg105=reg17*reg105; reg91=reg40*reg91;
    reg15=reg63+reg15; reg63=reg9*reg83; reg20=reg102+reg20; reg23=reg41*reg23; reg120=reg16*reg120;
    reg102=reg106*reg11; reg74=reg81+reg74; reg65=reg41*reg65; reg126=reg16*reg126; reg115=reg73+reg115;
    reg73=reg109*reg84; reg81=elem.pos(1)[1]-elem.pos(0)[1]; T reg128=elem.pos(2)[0]-elem.pos(0)[0]; T reg129=elem.pos(2)[1]-elem.pos(0)[1]; T reg130=elem.pos(1)[0]-elem.pos(0)[0];
    T reg131=reg9*reg80; T reg132=reg12*reg83; T reg133=reg109*reg33; reg121=reg122+reg121; reg122=reg109*reg11;
    reg105=reg91+reg105; reg60=reg18/reg60; reg84=reg110*reg84; reg126=reg65+reg126; reg33=reg110*reg33;
    reg120=reg23+reg120; reg18=reg15*reg107; reg11=reg110*reg11; reg23=reg19*reg100; reg65=reg124*reg100;
    reg91=reg125*reg100; reg97=reg113+reg97; reg90=reg71+reg90; reg71=reg56*reg51; reg104=reg54+reg104;
    reg117=reg119-reg117; reg127=reg63+reg127; reg102=reg74+reg102; reg54=reg56*reg114; reg63=reg103*reg101;
    reg10=reg20+reg10; reg20=reg56*reg53; reg74=reg127*reg60; reg113=reg24*reg114; reg11=reg97+reg11;
    reg20=reg10+reg20; reg10=reg19*reg65; reg97=reg24*reg51; reg84=reg126+reg84; reg119=reg19*reg91;
    reg126=reg19*reg23; reg71=reg104+reg71; reg33=reg120+reg33; reg104=reg24*reg53; reg54=reg102+reg54;
    reg102=reg101*reg93; reg120=reg81*reg128; T reg134=reg130*reg129; T reg135=reg35*reg15; T reg136=reg55*reg90;
    T reg137=reg101*reg55; reg51=reg22*reg51; reg73=reg115+reg73; reg115=reg15*reg93; T reg138=reg90*reg117;
    T reg139=reg35*reg90; reg18=reg63-reg18; reg8=reg83*reg8; reg53=reg22*reg53; reg133=reg121+reg133;
    reg63=reg132*reg60; reg121=reg131*reg60; reg114=reg22*reg114; reg122=reg105+reg122; reg13=reg13*reg80;
    reg39=reg39*reg83; reg135=reg137-reg135; reg95=reg95*reg80; reg120=reg134-reg120; reg105=reg123*reg15;
    reg134=reg13*reg23; reg51=reg73+reg51; reg73=reg13*reg65; reg53=reg133+reg53; reg133=reg123*reg90;
    reg137=reg13*reg91; reg114=reg122+reg114; reg139=reg102-reg139; reg136=reg115-reg136; reg23=reg8*reg23;
    reg102=reg101*reg79; reg97=reg84+reg97; reg84=reg99*reg90; reg115=reg101*reg99; reg122=reg15*reg79;
    reg138=reg18+reg138; reg18=reg127*reg74; reg91=reg8*reg91; reg113=reg11+reg113; reg119=reg54+reg119;
    reg11=reg127*reg63; reg65=reg8*reg65; reg126=reg71+reg126; reg104=reg33+reg104; reg10=reg20+reg10;
    reg20=reg127*reg121; reg33=reg39*reg63; reg54=reg95*reg121; reg73=reg53+reg73; reg63=reg95*reg63;
    reg137=reg114+reg137; reg128=reg128/reg120; reg81=reg81/reg120; reg53=reg39*reg74; reg23=reg97+reg23;
    reg65=reg104+reg65; reg121=reg39*reg121; reg18=reg126+reg18; reg20=reg10+reg20; reg11=reg119+reg11;
    reg91=reg113+reg91; reg84=reg122-reg84; reg105=reg115-reg105; reg136=reg136/reg138; reg135=reg135/reg138;
    reg10=1-(*f.m).resolution; reg134=reg51+reg134; reg74=reg95*reg74; reg139=reg139/reg138; reg129=reg129/reg120;
    reg133=reg102-reg133; reg130=reg130/reg120; reg117=reg117/reg138; reg133=reg133/reg138; reg51=(*f.m).resolution*reg136;
    reg84=reg84/reg138; reg71=reg20*reg10; reg97=reg18*reg10; reg102=(*f.m).resolution*reg139; reg107=reg107/reg138;
    reg104=(*f.m).resolution*reg135; reg105=reg105/reg138; reg113=reg11*reg10; reg53=reg23+reg53; reg63=reg137+reg63;
    reg121=reg65+reg121; reg33=reg91+reg33; reg138=reg103/reg138; reg54=reg73+reg54; reg74=reg134+reg74;
    reg23=reg128-reg130; reg65=reg81-reg129; reg73=(*f.m).resolution*reg117; reg91=(*f.m).resolution*reg138; reg103=(*f.m).resolution*reg107;
    reg114=0.5*reg65; reg115=0.5*reg128; reg119=0.5*reg129; reg122=0.5*reg130; reg126=0.5*reg81;
    reg97=reg104+reg97; reg102=reg71-reg102; reg113=reg51+reg113; reg51=reg74*reg10; reg71=reg54*reg10;
    reg104=reg63*reg10; reg134=reg53*reg10; reg137=reg121*reg10; T reg140=reg33*reg10; T reg141=(*f.m).resolution*reg84;
    T reg142=(*f.m).resolution*reg133; T reg143=(*f.m).resolution*reg105; T reg144=0.5*reg23; T reg145=reg115*reg97; reg71=reg142+reg71;
    reg142=reg128*reg102; reg141=reg104-reg141; reg104=reg119*reg97; reg134=reg73+reg134; reg140=reg91+reg140;
    reg73=reg129*reg113; reg143=reg51-reg143; reg51=reg126*reg97; reg91=reg23*reg102; T reg146=reg114*reg97;
    reg103=reg137-reg103; reg137=reg81*reg113; T reg147=reg122*reg97; T reg148=reg65*reg113; T reg149=reg144*reg97;
    T reg150=reg130*reg102; T reg151=reg126*reg134; T reg152=reg129*reg140; T reg153=reg115*reg134; T reg154=reg130*reg103;
    reg147=reg147-reg137; T reg155=reg130*reg71; reg73=reg73-reg145; reg146=reg91+reg146; reg150=reg150-reg51;
    reg91=reg126*reg143; T reg156=reg23*reg71; T reg157=reg129*reg141; T reg158=reg115*reg143; T reg159=reg114*reg134;
    T reg160=reg23*reg103; reg149=reg148+reg149; reg148=reg122*reg134; T reg161=reg81*reg140; T reg162=reg128*reg103;
    T reg163=reg119*reg134; T reg164=reg144*reg134; T reg165=reg119*reg143; T reg166=reg128*reg71; reg104=reg104-reg142;
    T reg167=reg65*reg140; T reg168=reg81*reg141; T reg169=reg122*reg143; T reg170=reg114*reg143; reg73=2*reg73;
    reg163=reg163-reg162; reg150=2*reg150; reg152=reg152-reg153; reg104=2*reg104; reg157=reg157-reg158;
    reg159=reg160+reg159; reg149=2*reg149; reg164=reg167+reg164; reg148=reg148-reg161; reg170=reg156+reg170;
    reg169=reg169-reg168; reg155=reg155-reg91; reg165=reg165-reg166; reg147=2*reg147; reg146=2*reg146;
    reg154=reg154-reg151; reg156=reg65*reg163; reg160=reg23*reg155; reg167=reg23*reg165; T reg171=reg126*reg150;
    T reg172=reg144*reg73; T reg173=reg114*reg104; T reg174=reg144*reg104; T reg175=reg65*reg148; T reg176=reg65*reg159;
    T reg177=reg130*reg155; T reg178=reg144*reg149; T reg179=reg65*reg164; T reg180=reg114*reg147; reg155=reg128*reg155;
    T reg181=reg119*reg150; T reg182=reg128*reg169; T reg183=reg119*reg147; T reg184=reg23*reg169; T reg185=reg128*reg165;
    T reg186=reg119*reg104; T reg187=reg144*reg147; T reg188=reg114*reg150; T reg189=reg65*reg154; T reg190=reg129*reg163;
    T reg191=reg115*reg104; T reg192=reg65*reg152; T reg193=reg129*reg148; T reg194=reg115*reg147; T reg195=reg23*reg170;
    T reg196=reg144*reg150; reg148=reg81*reg148; T reg197=reg122*reg147; T reg198=reg114*reg146; T reg199=reg114*reg73;
    T reg200=reg23*reg157; T reg201=reg129*reg154; T reg202=reg115*reg150; T reg203=reg144*reg146; T reg204=reg129*reg152;
    reg150=reg122*reg150; reg154=reg81*reg154; T reg205=reg115*reg73; reg187=reg175+reg187; reg174=reg156+reg174;
    reg191=reg190-reg191; reg203=reg176+reg203; reg185=reg186-reg185; reg184=reg180+reg184; reg148=reg197-reg148;
    reg160=reg188+reg160; reg182=reg183-reg182; reg194=reg193-reg194; reg205=reg204-reg205; reg198=reg195+reg198;
    reg155=reg181-reg155; reg154=reg150-reg154; reg196=reg189+reg196; reg167=reg173+reg167; reg178=reg179+reg178;
    reg202=reg201-reg202; reg200=reg199+reg200; reg171=reg177-reg171; reg192=reg172+reg192; reg171=reg120*reg171;
    reg184=reg120*reg184; reg148=reg120*reg148; reg160=reg120*reg160; reg154=reg120*reg154; reg187=reg120*reg187;
    reg191=reg120*reg191; reg185=reg120*reg185; reg205=reg120*reg205; reg182=reg120*reg182; reg167=reg120*reg167;
    reg155=reg120*reg155; reg178=reg120*reg178; reg174=reg120*reg174; reg203=reg120*reg203; reg200=reg120*reg200;
    reg196=reg120*reg196; reg192=reg120*reg192; reg194=reg120*reg194; reg202=reg120*reg202; reg198=reg120*reg198;
    matrix(indices[1]+0,indices[1]+1)+=ponderation*reg191; matrix(indices[2]+0,indices[2]+0)+=ponderation*reg148; matrix(indices[1]+0,indices[1]+0)+=ponderation*reg205; matrix(indices[0]+0,indices[1]+1)+=ponderation*reg174; matrix(indices[2]+0,indices[2]+1)+=ponderation*reg154;
    matrix(indices[0]+0,indices[1]+0)+=ponderation*reg192; matrix(indices[2]+1,indices[2]+1)+=ponderation*reg171; matrix(indices[0]+1,indices[0]+1)+=ponderation*reg198; matrix(indices[0]+1,indices[2]+1)+=ponderation*reg160; matrix(indices[0]+1,indices[1]+0)+=ponderation*reg200;
    matrix(indices[0]+1,indices[1]+1)+=ponderation*reg167; matrix(indices[0]+1,indices[2]+0)+=ponderation*reg184; matrix(indices[0]+0,indices[2]+0)+=ponderation*reg187; matrix(indices[1]+1,indices[1]+1)+=ponderation*reg185; matrix(indices[1]+1,indices[2]+0)+=ponderation*reg182;
    matrix(indices[1]+1,indices[2]+1)+=ponderation*reg155; matrix(indices[0]+0,indices[0]+0)+=ponderation*reg178; matrix(indices[0]+0,indices[0]+1)+=ponderation*reg203; matrix(indices[0]+0,indices[2]+1)+=ponderation*reg196; matrix(indices[1]+0,indices[2]+0)+=ponderation*reg194;
    matrix(indices[1]+0,indices[2]+1)+=ponderation*reg202;
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
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=pow((*f.m).v1[0],2); T reg5=pow((*f.m).v1[1],2); T reg6=reg2*reg3;
    T reg7=pow((*f.m).v2[1],2); T reg8=pow((*f.m).v2[0],2); T reg9=1.0/(*f.m).elastic_modulus_3; T reg10=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg11=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg12=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg13=1.0/(*f.m).elastic_modulus_2; reg7=reg8+reg7; reg8=pow((*f.m).v2[2],2); reg5=reg4+reg5;
    reg4=reg9*reg6; T reg14=pow((*f.m).v1[2],2); T reg15=reg10*reg6; T reg16=reg11*reg6; reg8=reg7+reg8;
    reg5=reg14+reg5; reg7=reg13*reg4; reg14=reg10*reg15; T reg17=reg12*reg4; T reg18=reg10*reg16;
    reg5=pow(reg5,0.5); reg14=reg7-reg14; reg7=1.0/(*f.m).elastic_modulus_1; reg18=reg17+reg18; T reg19=reg12*reg15;
    T reg20=reg13*reg16; reg8=pow(reg8,0.5); T reg21=reg7*reg14; T reg22=(*f.m).v1[1]/reg5; T reg23=reg12*reg18;
    T reg24=(*f.m).v2[2]/reg8; T reg25=reg19+reg20; T reg26=(*f.m).v2[1]/reg8; T reg27=(*f.m).v1[2]/reg5; T reg28=reg2*reg0;
    reg5=(*f.m).v1[0]/reg5; T reg29=reg10*reg3; T reg30=reg22*reg24; T reg31=reg27*reg26; T reg32=reg9*reg3;
    reg3=reg11*reg3; T reg33=reg12*reg6; T reg34=reg11*reg16; reg4=reg7*reg4; reg8=(*f.m).v2[0]/reg8;
    reg6=reg13*reg6; T reg35=reg11*reg15; T reg36=reg11*reg25; reg23=reg21-reg23; reg21=reg10*reg28;
    T reg37=2*reg5; T reg38=reg11*reg28; T reg39=reg2*reg1; reg28=reg9*reg28; reg36=reg23-reg36;
    reg23=reg30-reg31; reg35=reg17+reg35; reg17=reg11*reg33; T reg40=reg10*reg29; T reg41=reg13*reg32;
    reg32=reg12*reg32; T reg42=2*reg8; T reg43=reg11*reg6; T reg44=reg27*reg8; reg15=reg7*reg15;
    T reg45=reg5*reg24; reg34=reg4-reg34; reg16=reg12*reg16; reg4=reg10*reg3; reg33=reg12*reg33;
    reg6=reg7*reg6; T reg46=reg22*reg8; T reg47=pow(reg22,2); T reg48=pow(reg5,2); reg17=reg15+reg17;
    T reg49=reg5*reg26; T reg50=reg44-reg45; reg15=reg16+reg15; reg16=reg9*reg39; reg43=reg19+reg43;
    T reg51=reg10*reg38; reg35=reg35/reg36; reg18=reg18/reg36; T reg52=reg10*reg21; T reg53=reg12*reg28;
    reg28=reg13*reg28; reg4=reg32+reg4; reg32=2*reg23; reg40=reg41-reg40; reg34=reg34/reg36;
    reg14=reg14/reg36; reg41=reg11*reg39; T reg54=reg42*reg26; T reg55=pow(reg26,2); T reg56=pow(reg8,2);
    reg39=reg10*reg39; reg29=reg12*reg29; reg3=reg13*reg3; T reg57=reg22*reg37; reg43=reg43/reg36;
    T reg58=reg56*reg18; reg17=reg17/reg36; reg15=reg15/reg36; T reg59=reg14*reg48; T reg60=reg54*reg34;
    reg33=reg6-reg33; reg6=pow(reg23,2); reg40=reg7*reg40; T reg61=pow(reg50,2); T reg62=reg32*reg50;
    reg4=reg12*reg4; T reg63=reg29+reg3; T reg64=reg56*reg34; reg52=reg28-reg52; reg28=reg35*reg48;
    reg25=reg25/reg36; reg51=reg53+reg51; reg53=reg13*reg16; reg16=reg12*reg16; T reg65=reg10*reg39;
    T reg66=reg10*reg41; T reg67=reg14*reg47; T reg68=reg55*reg18; T reg69=reg14*reg57; T reg70=reg54*reg18;
    T reg71=reg35*reg57; T reg72=reg49-reg46; T reg73=pow(reg27,2); reg21=reg12*reg21; reg38=reg13*reg38;
    T reg74=pow(reg24,2); T reg75=reg35*reg47; T reg76=reg55*reg34; T reg77=reg55*reg17; T reg78=reg43*reg47;
    T reg79=reg56*reg17; reg65=reg53-reg65; reg53=reg43*reg48; reg66=reg16+reg66; reg16=reg62*reg15;
    reg60=reg71+reg60; reg68=reg67+reg68; reg67=reg61*reg25; reg71=reg14*reg73; T reg80=reg74*reg18;
    T reg81=reg74*reg34; reg70=reg69+reg70; reg69=reg62*reg25; T reg82=reg21+reg38; T reg83=reg6*reg15;
    reg51=reg12*reg51; reg39=reg12*reg39; T reg84=reg35*reg73; T reg85=reg61*reg15; reg52=reg7*reg52;
    reg76=reg75+reg76; reg63=reg11*reg63; reg59=reg58+reg59; reg64=reg28+reg64; reg28=reg43*reg57;
    reg58=reg54*reg17; reg4=reg40-reg4; reg40=reg6*reg25; reg41=reg13*reg41; reg33=reg33/reg36;
    reg75=pow(reg72,2); reg63=reg4-reg63; reg51=reg52-reg51; reg82=reg11*reg82; reg65=reg7*reg65;
    reg66=reg12*reg66; reg4=reg39+reg41; reg40=reg59+reg40; reg67=reg68+reg67; reg80=reg71+reg80;
    reg52=reg25*reg75; reg69=reg70+reg69; reg83=reg64+reg83; reg59=2*reg26; reg64=reg42*reg24;
    reg68=2*reg22; reg85=reg76+reg85; reg70=reg27*reg37; reg81=reg84+reg81; reg71=reg75*reg15;
    reg16=reg60+reg16; reg79=reg53+reg79; reg53=reg6*reg33; reg77=reg78+reg77; reg60=reg61*reg33;
    reg76=reg43*reg73; reg78=reg74*reg17; reg58=reg28+reg58; reg28=reg62*reg33; reg84=reg8*reg26;
    T reg86=reg5*reg22; T reg87=reg50*reg23; T reg88=reg86*reg67; T reg89=reg84*reg85; T reg90=reg59*reg24;
    reg82=reg51-reg82; reg32=reg32*reg72; reg51=reg75*reg33; reg78=reg76+reg78; reg76=2*reg50;
    T reg91=reg8*reg37; reg60=reg77+reg60; reg66=reg65-reg66; reg65=reg5*reg50; reg4=reg11*reg4;
    reg77=reg68*reg26; T reg92=reg22*reg23; T reg93=reg56*reg16; T reg94=reg48*reg69; T reg95=reg86*reg69;
    T reg96=reg84*reg16; T reg97=reg64*reg18; T reg98=2*reg27; reg52=reg80+reg52; reg49=reg46+reg49;
    reg46=reg70*reg14; reg80=reg22*reg26; T reg99=reg5*reg8; T reg100=reg56*reg83; T reg101=reg48*reg40;
    reg16=reg55*reg16; reg69=reg47*reg69; T reg102=reg64*reg34; T reg103=reg70*reg35; reg28=reg58+reg28;
    reg58=reg48*reg67; T reg104=reg56*reg85; reg85=reg55*reg85; reg67=reg47*reg67; T reg105=reg68*reg27;
    reg71=reg81+reg71; reg63=reg63/reg36; reg53=reg79+reg53; reg79=reg47*reg40; reg81=reg55*reg83;
    reg100=reg101+reg100; reg101=reg6*reg53; T reg106=reg56*reg71; reg34=reg90*reg34; reg35=reg105*reg35;
    T reg107=reg55*reg13; T reg108=reg56*reg7; T reg109=reg32*reg15; reg97=reg46+reg97; reg46=reg32*reg25;
    reg102=reg103+reg102; reg103=reg48*reg52; T reg110=reg87*reg28; reg14=reg105*reg14; T reg111=reg56*reg12;
    reg18=reg90*reg18; T reg112=reg55*reg12; reg76=reg76*reg72; T reg113=reg6*reg60; reg104=reg58+reg104;
    reg58=reg80*reg63; reg83=reg84*reg83; reg81=reg79+reg81; reg79=reg61*reg53; T reg114=reg77*reg13;
    T reg115=reg91*reg7; T reg116=reg99*reg63; reg40=reg86*reg40; reg85=reg67+reg85; reg67=reg61*reg60;
    T reg117=reg70*reg43; T reg118=reg47*reg52; T reg119=reg55*reg71; T reg120=reg64*reg17; reg16=reg69+reg16;
    reg69=reg26*reg23; T reg121=reg8*reg50; T reg122=reg61*reg28; T reg123=reg5*reg23; T reg124=reg22*reg50;
    reg4=reg66-reg4; reg93=reg94+reg93; reg28=reg6*reg28; reg96=reg95+reg96; reg65=reg92+reg65;
    reg66=reg27*reg24; reg92=reg49*reg63; reg98=reg98*reg24; reg94=reg77*reg12; reg51=reg78+reg51;
    reg60=reg87*reg60; reg89=reg88+reg89; reg78=reg91*reg12; reg82=reg82/reg36; reg17=reg90*reg17;
    reg43=reg105*reg43; reg112=reg108-reg112; reg88=reg74*reg11; reg95=reg123*reg82; reg108=reg65*reg82;
    T reg125=reg124*reg82; T reg126=reg66*reg63; T reg127=reg77*reg92; reg122=reg16+reg122; reg16=reg61*reg51;
    reg119=reg118+reg119; reg118=reg77*reg58; reg67=reg85+reg67; reg85=reg77*reg116; reg79=reg81+reg79;
    reg83=reg40+reg83; reg53=reg87*reg53; reg60=reg89+reg60; reg40=reg49*reg58; reg52=reg86*reg52;
    reg71=reg84*reg71; reg81=reg91*reg92; reg28=reg93+reg28; reg110=reg96+reg110; reg89=reg6*reg51;
    reg106=reg103+reg106; reg92=reg49*reg92; reg58=reg91*reg58; reg113=reg104+reg113; reg93=reg91*reg116;
    reg101=reg100+reg101; reg96=reg55*reg10; reg100=reg11*reg91; reg103=reg77*reg10; reg15=reg76*reg15;
    reg34=reg35+reg34; reg109=reg102+reg109; reg35=reg98*reg10; reg102=reg26*reg50; reg104=reg8*reg23;
    T reg128=reg74*reg10; reg111=reg107-reg111; reg25=reg76*reg25; reg18=reg14+reg18; reg46=reg97+reg46;
    reg36=reg4/reg36; reg37=reg23*reg37; reg68=reg68*reg50; reg4=reg27*reg72; reg94=reg115-reg94;
    reg14=reg11*reg98; reg120=reg117+reg120; reg97=reg56*reg11; reg121=reg69+reg121; reg69=reg32*reg33;
    reg78=reg114-reg78; reg35=reg78-reg35; reg78=reg27*reg23; reg81=reg28+reg81; reg42=reg42*reg23;
    reg28=reg37*reg108; reg51=reg87*reg51; reg107=reg48*reg46; reg114=reg56*reg109; reg71=reg52+reg71;
    reg52=reg5*reg72; reg115=reg65*reg125; reg118=reg67+reg118; reg67=reg68*reg125; reg117=reg74*reg9;
    reg96=reg97+reg96; reg97=reg55*reg109; reg98=reg98*reg9; T reg129=reg68*reg95; reg85=reg79+reg85;
    reg16=reg119+reg16; reg79=reg77*reg126; reg119=reg47*reg46; reg103=reg100+reg103; reg100=(*f.m).alpha_1*reg48;
    T reg130=reg68*reg108; T reg131=reg24*reg72; reg53=reg83+reg53; reg116=reg49*reg116; reg127=reg122+reg127;
    reg128=reg111-reg128; reg40=reg60+reg40; reg60=reg56*(*f.m).alpha_2; reg13=reg13*reg47; reg83=reg37*reg95;
    reg93=reg101+reg93; reg15=reg34+reg15; reg34=(*f.m).alpha_1*reg47; reg101=reg55*(*f.m).alpha_2; reg7=reg7*reg48;
    reg111=reg12*reg47; reg122=reg121*reg36; T reg132=reg102*reg36; T reg133=reg104*reg36; reg12=reg12*reg48;
    T reg134=reg4*reg82; reg14=reg94-reg14; reg88=reg112-reg88; reg33=reg76*reg33; reg17=reg43+reg17;
    reg69=reg120+reg69; reg108=reg65*reg108; reg44=reg45+reg44; reg59=reg59*reg50; reg43=reg91*reg126;
    reg92=reg110+reg92; reg89=reg106+reg89; reg125=reg37*reg125; reg58=reg113+reg58; reg25=reg18+reg25;
    reg18=reg11*reg73; reg45=reg99*reg14; reg51=reg71+reg51; reg71=reg61*reg69; reg97=reg119+reg97;
    reg111=reg7-reg111; reg96=reg117-reg96; reg7=reg10*reg73; reg10=reg10*reg47; reg11=reg11*reg48;
    reg126=reg49*reg126; reg101=reg34+reg101; reg34=reg59*reg122; reg130=reg127+reg130; reg94=reg61*(*f.m).alpha_3;
    reg106=(*f.m).alpha_1*reg73; reg110=reg74*(*f.m).alpha_2; reg112=reg80*reg128; reg113=reg56*reg14; reg117=reg55*reg128;
    reg119=reg55*reg35; reg60=reg100+reg60; reg100=reg6*(*f.m).alpha_3; reg95=reg65*reg95; reg116=reg53+reg116;
    reg53=reg80*reg35; reg108=reg92+reg108; reg92=reg121*reg122; reg46=reg86*reg46; reg120=reg56*reg88;
    reg127=reg47*reg35; T reg135=reg48*reg14; T reg136=reg47*reg128; T reg137=reg48*reg88; reg12=reg13-reg12;
    reg115=reg40+reg115; reg13=reg121*reg132; reg40=reg99*reg88; T reg138=reg55*reg15; T reg139=reg47*reg25;
    reg103=reg98-reg103; reg109=reg84*reg109; reg83=reg93+reg83; reg52=reg78+reg52; reg30=reg31+reg30;
    reg129=reg85+reg129; reg31=reg42*reg133; reg125=reg58+reg125; reg58=reg42*reg132; reg78=reg8*reg72;
    reg43=reg89+reg43; reg85=reg37*reg134; reg89=reg24*reg23; reg93=reg56*reg15; reg98=reg48*reg25;
    T reg140=reg131*reg36; reg28=reg81+reg28; reg122=reg42*reg122; reg81=reg6*reg69; reg114=reg107+reg114;
    reg132=reg59*reg132; reg67=reg118+reg67; reg79=reg16+reg79; reg33=reg17+reg33; reg16=reg44*reg63;
    reg17=reg27*reg50; reg107=reg22*reg72; reg118=reg59*reg133; T reg141=reg68*reg134; reg134=reg65*reg134;
    reg126=reg51+reg126; reg9=reg9*reg73; reg63=reg30*reg63; reg10=reg11+reg10; reg122=reg28+reg122;
    reg15=reg84*reg15; reg11=reg73*reg96; reg18=reg111-reg18; reg28=reg42*reg140; reg51=reg52*reg82;
    reg85=reg43+reg85; reg25=reg86*reg25; reg43=(*f.m).alpha_2*reg84; reg58=reg125+reg58; reg111=(*f.m).alpha_1*reg86;
    reg31=reg83+reg31; reg83=reg75*(*f.m).alpha_3; reg110=reg106+reg110; reg92=reg108+reg92; reg136=reg137+reg136;
    reg109=reg46+reg109; reg69=reg87*reg69; reg100=reg60+reg100; reg94=reg101+reg94; reg107=reg17+reg107;
    reg117=reg120+reg117; reg17=reg74*reg96; reg7=reg12-reg7; reg12=reg61*reg33; reg138=reg139+reg138;
    reg118=reg129+reg118; reg119=reg113+reg119; reg46=reg74*reg103; reg5=reg5*reg27; reg60=reg77*reg16;
    reg71=reg97+reg71; reg132=reg67+reg132; reg34=reg130+reg34; reg141=reg79+reg141; reg67=reg59*reg140;
    reg8=reg8*reg24; reg81=reg114+reg81; reg13=reg115+reg13; reg79=reg91*reg16; reg97=reg49*reg2;
    reg101=reg66*reg103; reg127=reg135+reg127; reg106=reg73*reg103; reg84=reg84*reg2; reg93=reg98+reg93;
    reg98=reg6*reg33; reg133=reg121*reg133; reg95=reg116+reg95; reg53=reg45+reg53; reg78=reg89+reg78;
    reg45=reg24*reg50; reg89=reg26*reg72; reg108=reg66*reg96; reg112=reg40+reg112; reg40=(*f.m).alpha_2*reg8;
    reg113=reg132*reg92; reg114=reg58*reg92; reg10=reg9-reg10; reg9=reg34*reg13; reg89=reg45+reg89;
    reg45=reg54*reg97; reg46=reg119+reg46; reg115=reg54*reg84; reg17=reg117+reg17; reg108=reg112+reg108;
    reg112=reg49*reg84; reg116=reg55*reg7; reg117=reg56*reg18; reg119=reg57*reg97; reg106=reg127+reg106;
    reg101=reg53+reg101; reg53=reg49*reg97; reg120=reg57*reg84; reg11=reg136+reg11; reg83=reg110+reg83;
    reg110=reg47*reg7; reg125=reg48*reg18; reg43=reg111+reg43; reg111=(*f.m).alpha_3*reg87; reg127=reg44*reg1;
    reg8=reg8*reg1; reg129=(*f.m).alpha_1*reg5; reg130=reg31*reg100; reg135=reg58*reg94; reg33=reg87*reg33;
    reg15=reg25+reg15; reg25=reg118*reg100; reg87=reg132*reg94; reg16=reg49*reg16; reg69=reg109+reg69;
    reg28=reg85+reg28; reg140=reg121*reg140; reg134=reg126+reg134; reg79=reg81+reg79; reg81=reg37*reg51;
    reg133=reg95+reg133; reg98=reg93+reg98; reg91=reg91*reg63; reg77=reg77*reg63; reg12=reg138+reg12;
    reg85=reg68*reg51; reg60=reg71+reg60; reg67=reg141+reg67; reg71=reg78*reg36; reg23=reg72*reg23;
    reg27=reg22*reg27; reg24=reg26*reg24; reg22=reg122*reg13; reg82=reg107*reg82; reg112=reg108+reg112;
    reg26=reg44*reg8; reg74=reg74*reg10; reg116=reg117+reg116; reg91=reg98+reg91; reg9=reg113-reg9;
    reg93=reg70*reg127; reg119=reg106+reg119; reg87=reg25+reg87; reg25=reg42*reg71; reg81=reg79+reg81;
    reg22=reg114-reg22; reg37=reg37*reg82; reg67=reg67*reg83; reg79=reg133*reg100; reg95=reg13*reg94;
    reg68=reg68*reg82; reg77=reg12+reg77; reg50=reg72*reg50; reg115=reg17+reg115; reg12=reg64*reg8;
    reg17=reg80*reg7; reg72=reg59*reg71; reg85=reg60+reg85; reg45=reg46+reg45; reg46=reg64*reg127;
    reg60=reg99*reg18; reg2=reg86*reg2; reg40=reg129+reg40; reg135=reg130+reg135; reg111=reg43+reg111;
    reg28=reg28*reg83; reg63=reg49*reg63; reg43=reg24*reg0; reg86=reg30*reg0; reg33=reg15+reg33;
    reg23=(*f.m).alpha_3*reg23; reg110=reg125+reg110; reg73=reg73*reg10; reg15=(*f.m).alpha_1*reg27; reg24=(*f.m).alpha_2*reg24;
    reg51=reg65*reg51; reg16=reg69+reg16; reg69=reg58*reg34; reg98=reg122*reg132; reg53=reg101+reg53;
    reg101=reg44*reg127; reg106=reg70*reg8; reg120=reg11+reg120; reg36=reg89*reg36; reg140=reg134+reg140;
    reg28=reg135+reg28; reg98=reg69-reg98; reg140=reg140*reg83; reg17=reg60+reg17; reg66=reg66*reg10;
    reg11=reg118*reg22; reg60=reg30*reg43; reg23=reg40+reg23; reg40=reg31*reg9; reg26=reg112+reg26;
    reg95=reg79+reg95; reg69=reg122*reg111; reg101=reg53+reg101; reg53=reg34*reg111; reg67=reg87+reg67;
    reg79=reg30*reg86; reg50=(*f.m).alpha_3*reg50; reg24=reg15+reg24; reg71=reg121*reg71; reg51=reg16+reg51;
    reg15=reg54*reg2; reg74=reg116+reg74; reg25=reg81+reg25; reg16=reg105*reg86; reg37=reg91+reg37;
    reg42=reg42*reg36; reg93=reg119+reg93; reg81=reg105*reg43; reg106=reg120+reg106; reg59=reg59*reg36;
    reg68=reg77+reg68; reg72=reg85+reg72; reg1=reg5*reg1; reg73=reg110+reg73; reg5=reg57*reg2;
    reg77=reg90*reg86; reg46=reg45+reg46; reg63=reg33+reg63; reg33=reg90*reg43; reg12=reg115+reg12;
    reg82=reg65*reg82; reg70=reg70*reg1; reg5=reg73+reg5; reg71=reg51+reg71; reg0=reg27*reg0;
    reg27=reg92*reg111; reg140=reg95+reg140; reg50=reg24+reg50; reg82=reg63+reg82; reg59=reg68+reg59;
    reg72=reg72*reg23; reg53=reg67+reg53; reg36=reg121*reg36; reg42=reg37+reg42; reg69=reg28+reg69;
    reg25=reg25*reg23; reg24=reg122*reg133; reg28=reg31*reg92; reg37=reg34*reg133; reg45=reg118*reg92;
    reg51=reg133*reg98; reg11=reg40-reg11; reg77=reg46+reg77; reg33=reg12+reg33; reg66=reg17+reg66;
    reg12=reg49*reg2; reg64=reg64*reg1; reg15=reg74+reg15; reg16=reg93+reg16; reg60=reg26+reg60;
    reg81=reg106+reg81; reg79=reg101+reg79; reg59=reg59*reg50; reg25=reg69+reg25; reg17=reg31*reg13;
    reg72=reg53+reg72; reg24=reg28-reg24; reg26=reg132*reg133; reg28=reg58*reg133; reg42=reg42*reg50;
    reg37=reg45-reg37; reg40=reg118*reg13; reg45=reg31*reg34; reg46=reg122*reg118; reg44=reg44*reg1;
    reg53=reg16*reg60; reg64=reg15+reg64; reg70=reg5+reg70; reg5=reg33*reg79; reg71=reg71*reg23;
    reg27=reg140+reg27; reg105=reg105*reg0; reg12=reg66+reg12; reg51=reg11+reg51; reg90=reg90*reg0;
    reg11=reg77*reg60; reg15=reg81*reg79; reg36=reg82+reg36; reg37=reg37/reg51; reg42=reg25+reg42;
    reg90=reg64+reg90; reg9=reg9/reg51; reg26=reg40-reg26; reg11=reg5-reg11; reg36=reg36*reg50;
    reg71=reg27+reg71; reg105=reg70+reg105; reg5=reg58*reg118; reg53=reg15-reg53; reg15=reg81*reg77;
    reg25=reg16*reg33; reg46=reg45-reg46; reg27=reg31*reg132; reg44=reg12+reg44; reg28=reg17-reg28;
    reg24=reg24/reg51; reg30=reg30*reg0; reg59=reg72+reg59; reg22=reg22/reg51; reg46=reg46/reg51;
    reg98=reg98/reg51; reg30=reg44+reg30; reg28=reg28/reg51; reg26=reg26/reg51; reg5=reg27-reg5;
    reg22=reg22*reg59; reg12=reg90*reg53; reg17=reg11*reg105; reg9=reg9*reg42; reg36=reg71+reg36;
    reg37=reg37*reg42; reg24=reg24*reg59; reg25=reg15-reg25; reg42=reg26*reg42; reg15=reg90*reg60;
    reg26=reg33*reg30; reg37=reg24-reg37; reg24=reg105*reg60; reg27=reg30*reg25; reg40=reg81*reg30;
    reg12=reg17-reg12; reg59=reg28*reg59; reg51=reg5/reg51; reg46=reg46*reg36; reg22=reg9-reg22;
    reg98=reg98*reg36; reg27=reg12+reg27; reg5=reg90*reg79; reg9=reg77*reg30; reg12=reg105*reg79;
    reg26=reg15-reg26; reg15=1-(*f.m).resolution; reg17=reg16*reg30; reg40=reg24-reg40; reg46=reg37-reg46;
    reg24=reg105*reg33; reg36=reg51*reg36; reg59=reg42-reg59; reg28=reg81*reg90; reg22=reg98+reg22;
    reg22=reg15*reg22; reg37=(*f.m).resolution*reg94; reg59=reg36+reg59; reg46=reg15*reg46; reg26=reg26/reg27;
    reg9=reg5-reg9; reg17=reg12-reg17; reg40=reg40/reg27; reg5=reg105*reg77; reg12=reg16*reg90;
    reg28=reg24-reg28; reg24=(*f.m).resolution*reg100; reg59=reg15*reg59; reg13=reg13*reg15; reg133=reg133*reg15;
    reg22=reg24+reg22; reg9=reg9/reg27; reg46=reg37+reg46; reg53=reg53/reg27; reg17=reg17/reg27;
    reg24=(*f.m).resolution*reg40; reg36=(*f.m).resolution*reg26; reg37=(*f.m).resolution*reg111; reg12=reg5-reg12; reg28=reg28/reg27;
    reg11=reg11/reg27; reg25=reg25/reg27; reg27=reg12/reg27; reg92=reg92*reg15; reg133=reg36+reg133;
    reg24=reg13-reg24; reg5=(*f.m).resolution*reg53; reg12=elem.pos(2)[1]-elem.pos(0)[1]; reg13=elem.pos(2)[0]-elem.pos(0)[0]; reg36=elem.pos(1)[1]-elem.pos(0)[1];
    reg42=elem.pos(1)[0]-elem.pos(0)[0]; reg44=(*f.m).resolution*reg11; reg59=reg37+reg59; reg22=(*f.m).deltaT*reg22; reg46=(*f.m).deltaT*reg46;
    reg132=reg132*reg15; reg118=reg118*reg15; reg37=(*f.m).resolution*reg28; reg31=reg31*reg15; reg45=(*f.m).resolution*reg17;
    reg51=(*f.m).resolution*reg9; reg58=reg58*reg15; reg59=(*f.m).deltaT*reg59; reg63=(*f.m).resolution*reg27; reg64=(*f.m).resolution*reg25;
    reg66=reg42*reg12; reg67=reg133*reg22; reg68=reg24*reg46; reg34=reg34*reg15; reg31=reg44+reg31;
    reg5=reg58-reg5; reg44=reg36*reg13; reg92=reg37+reg92; reg15=reg122*reg15; reg51=reg118-reg51;
    reg132=reg45+reg132; reg15=reg64+reg15; reg37=reg51*reg22; reg45=reg132*reg46; reg58=reg31*reg22;
    reg64=reg5*reg46; reg63=reg34-reg63; reg44=reg66-reg44; reg34=reg67+reg68; reg66=reg92*reg59;
    reg69=reg63*reg59; reg70=reg37+reg45; reg71=reg58+reg64; reg72=reg15*reg59; reg73=reg34+reg66;
    reg42=reg42/reg44; reg12=reg12/reg44; reg13=reg13/reg44; reg36=reg36/reg44; reg74=2*reg73;
    reg82=0.5*reg13; reg85=reg13-reg42; reg87=reg70+reg69; reg91=0.5*reg12; reg93=0.5*reg42;
    reg95=0.5*reg36; reg98=reg36-reg12; reg101=reg71+reg72; reg106=reg91*reg74; reg108=reg13*reg87;
    reg109=reg95*reg74; reg110=reg42*reg87; reg112=reg93*reg74; reg113=reg36*reg101; reg114=0.5*reg85;
    reg115=reg82*reg74; reg116=reg12*reg101; reg117=0.5*reg98; reg118=reg109-reg110; reg119=reg98*reg101;
    reg120=reg113-reg112; reg122=reg114*reg74; reg125=reg117*reg74; reg126=reg85*reg87; reg129=reg108-reg106;
    reg130=reg115-reg116; reg134=reg119+reg122; reg120=reg44*reg120; reg129=reg44*reg129; reg118=reg44*reg118;
    reg135=reg125+reg126; reg130=reg44*reg130; reg136=reg44*reg134; reg137=reg44*reg135; reg130=ponderation*reg130;
    reg118=ponderation*reg118; reg129=ponderation*reg129; reg120=ponderation*reg120; sollicitation[indices[1]+0]+=-reg130; reg130=ponderation*reg137;
    sollicitation[indices[0]+1]+=reg130; sollicitation[indices[2]+1]+=-reg118; sollicitation[indices[2]+0]+=-reg120; reg118=ponderation*reg136; sollicitation[indices[0]+0]+=reg118;
    sollicitation[indices[1]+1]+=-reg129;
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

#ifndef elasticity_orthotropy_stat_Qstat_read_material_to_mesh
#define elasticity_orthotropy_stat_Qstat_read_material_to_mesh
template<class TM, class T, bool wont_add_nz>
void read_material_to_mesh_(const XmlNode &n, Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f){ 
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

    if(n.has_attribute("elastic_modulus_1"))  
        n.get_attribute("elastic_modulus_1", f.m->elastic_modulus_1 ); 
    else  
        std::cerr << "Warning using default value of elastic_modulus_1 : " << f.m->elastic_modulus_1 << std::endl; 

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

/// @author hugo LECLERC

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
  static const unsigned order_integration = 2;
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
    static const bool has_elementary_matrix = true;
    static const bool has_skin_elementary_matrix = false;
    template<class TE,class TF, class TVEVE> static void after_solve(TE &elem,TF &f,TVEVE &vectors,const unsigned *indices) {
      #define PNODE(N) (*elem.node(N))
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; reg1=1.0/reg1; reg0=1.0/reg0; T reg2=2*(*f.m).shear_modulus_12;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=pow((*f.m).v1[1],2); T reg5=pow((*f.m).v1[0],2); T reg6=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2;
    T reg7=1.0/(*f.m).elastic_modulus_3; T reg8=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg9=reg2*reg3; T reg10=pow((*f.m).v2[1],2); T reg11=pow((*f.m).v2[0],2);
    T reg12=reg7*reg9; T reg13=reg6*reg9; T reg14=reg8*reg9; T reg15=pow((*f.m).v2[2],2); reg10=reg11+reg10;
    reg11=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg16=1.0/(*f.m).elastic_modulus_2; reg4=reg5+reg4; reg5=pow((*f.m).v1[2],2); reg15=reg10+reg15;
    reg5=reg4+reg5; reg4=reg16*reg12; reg10=reg6*reg13; T reg17=reg11*reg12; T reg18=reg6*reg14;
    T reg19=reg16*reg14; reg5=pow(reg5,0.5); T reg20=1.0/(*f.m).elastic_modulus_1; T reg21=reg11*reg13; reg15=pow(reg15,0.5);
    reg18=reg17+reg18; reg10=reg4-reg10; reg4=reg21+reg19; T reg22=(*f.m).v2[1]/reg15; T reg23=reg11*reg18;
    T reg24=(*f.m).v1[1]/reg5; T reg25=(*f.m).v2[2]/reg15; T reg26=reg20*reg10; T reg27=(*f.m).v1[2]/reg5; T reg28=reg8*reg14;
    reg15=(*f.m).v2[0]/reg15; reg12=reg20*reg12; reg23=reg26-reg23; reg26=reg8*reg13; T reg29=reg16*reg9;
    T reg30=reg8*reg4; reg9=reg11*reg9; T reg31=reg7*reg3; T reg32=reg2*reg0; T reg33=reg27*reg22;
    T reg34=reg24*reg25; T reg35=reg6*reg3; reg3=reg8*reg3; reg5=(*f.m).v1[0]/reg5; T reg36=reg16*reg31;
    reg31=reg11*reg31; T reg37=reg6*reg35; T reg38=reg6*reg3; reg30=reg23-reg30; reg23=2*reg15;
    T reg39=reg8*reg32; T reg40=2*reg5; T reg41=2*reg24; T reg42=reg2*reg1; T reg43=reg34-reg33;
    reg14=reg11*reg14; T reg44=reg6*reg32; T reg45=reg27*reg15; T reg46=reg8*reg9; T reg47=reg5*reg25;
    reg13=reg20*reg13; reg28=reg12-reg28; reg32=reg7*reg32; reg26=reg17+reg26; reg12=reg8*reg29;
    reg17=reg22*reg41; T reg48=2*reg27; T reg49=reg15*reg40; T reg50=2*reg43; T reg51=reg7*reg42;
    T reg52=reg22*reg23; T reg53=pow(reg22,2); T reg54=reg24*reg15; T reg55=reg5*reg22; T reg56=pow(reg15,2);
    T reg57=reg45-reg47; T reg58=pow(reg5,2); T reg59=pow(reg24,2); T reg60=reg24*reg40; reg18=reg18/reg30;
    reg38=reg31+reg38; reg12=reg21+reg12; reg10=reg10/reg30; reg28=reg28/reg30; reg29=reg20*reg29;
    reg26=reg26/reg30; reg31=reg16*reg32; reg46=reg13+reg46; T reg61=reg8*reg42; reg14=reg13+reg14;
    reg32=reg11*reg32; reg13=reg6*reg44; reg9=reg11*reg9; reg37=reg36-reg37; reg35=reg11*reg35;
    reg42=reg6*reg42; reg3=reg16*reg3; reg36=reg6*reg39; T reg62=reg53*reg16; T reg63=reg56*reg11;
    T reg64=reg17*reg11; reg37=reg20*reg37; reg38=reg11*reg38; T reg65=reg49*reg11; T reg66=reg17*reg16;
    reg39=reg16*reg39; reg12=reg12/reg30; T reg67=reg58*reg10; reg46=reg46/reg30; reg14=reg14/reg30;
    reg9=reg29-reg9; reg44=reg11*reg44; reg36=reg32+reg36; reg29=reg16*reg51; reg13=reg31-reg13;
    reg31=reg25*reg48; reg51=reg11*reg51; reg32=reg6*reg42; T reg68=reg6*reg61; T reg69=reg59*reg10;
    T reg70=reg53*reg18; T reg71=reg60*reg10; T reg72=reg52*reg18; T reg73=reg58*reg26; T reg74=reg56*reg28;
    T reg75=reg35+reg3; T reg76=reg59*reg26; T reg77=reg57*reg50; T reg78=pow(reg57,2); T reg79=pow(reg43,2);
    reg4=reg4/reg30; T reg80=reg56*reg18; T reg81=pow(reg25,2); T reg82=reg53*reg28; T reg83=reg60*reg26;
    T reg84=reg52*reg28; T reg85=reg56*reg20; T reg86=reg53*reg11; T reg87=reg49*reg20; T reg88=pow(reg27,2);
    T reg89=reg55-reg54; T reg90=reg77*reg4; reg38=reg37-reg38; reg61=reg16*reg61; reg42=reg11*reg42;
    reg37=reg79*reg4; reg67=reg80+reg67; reg74=reg73+reg74; reg73=reg79*reg14; reg82=reg76+reg82;
    reg76=reg78*reg14; reg80=reg81*reg8; reg86=reg85-reg86; reg85=reg88*reg26; T reg91=reg81*reg28;
    reg84=reg83+reg84; reg83=reg77*reg14; reg9=reg9/reg30; T reg92=reg31*reg8; reg64=reg87-reg64;
    reg87=reg58*reg12; T reg93=reg56*reg46; T reg94=reg59*reg12; T reg95=reg53*reg46; T reg96=reg60*reg12;
    T reg97=reg52*reg46; reg65=reg66-reg65; reg66=reg31*reg6; reg63=reg62-reg63; reg62=reg81*reg6;
    T reg98=reg78*reg4; reg70=reg69+reg70; reg69=reg53*reg6; reg68=reg51+reg68; reg51=reg88*reg10;
    reg32=reg29-reg32; reg29=pow(reg89,2); T reg99=reg49*reg8; T reg100=reg17*reg6; T reg101=reg44+reg39;
    reg36=reg11*reg36; reg75=reg8*reg75; reg72=reg71+reg72; reg71=reg81*reg18; T reg102=reg56*reg8;
    reg13=reg20*reg13; reg68=reg11*reg68; reg75=reg38-reg75; reg95=reg94+reg95; reg38=reg78*reg9;
    reg94=reg88*reg12; T reg103=reg81*reg46; T reg104=reg42+reg61; reg98=reg70+reg98; reg97=reg96+reg97;
    reg70=reg77*reg9; reg66=reg65-reg66; reg37=reg67+reg37; reg90=reg72+reg90; reg62=reg63-reg62;
    reg63=reg58*reg20; reg73=reg74+reg73; reg100=reg99+reg100; reg36=reg13-reg36; reg76=reg82+reg76;
    reg13=reg31*reg7; reg101=reg8*reg101; reg71=reg51+reg71; reg91=reg85+reg91; reg51=reg29*reg14;
    reg69=reg102+reg69; reg32=reg20*reg32; reg65=reg29*reg4; reg83=reg84+reg83; reg67=reg58*reg11;
    reg72=reg81*reg7; reg93=reg87+reg93; reg74=reg79*reg9; reg82=reg59*reg16; reg92=reg64-reg92;
    reg64=reg15*reg22; reg84=reg24*reg22; reg85=reg5*reg15; reg87=reg59*reg11; reg96=2*reg22;
    reg99=reg25*reg23; reg102=reg27*reg40; T reg105=reg5*reg24; reg80=reg86-reg80; reg86=reg105*reg90;
    T reg106=reg64*reg83; T reg107=reg56*reg83; T reg108=reg102*reg10; T reg109=reg99*reg18; T reg110=reg53*reg73;
    T reg111=reg59*reg37; T reg112=reg58*reg80; T reg113=reg59*reg66; T reg114=reg58*reg92; T reg115=reg59*reg62;
    T reg116=reg88*reg8; reg87=reg63-reg87; reg63=reg59*reg6; T reg117=reg27*reg41; reg69=reg72-reg69;
    reg72=reg105*reg98; T reg118=reg64*reg76; reg104=reg8*reg104; T reg119=reg25*reg96; reg68=reg32-reg68;
    reg50=reg89*reg50; reg100=reg13-reg100; reg13=2*reg57; reg101=reg36-reg101; reg32=reg53*reg83;
    reg36=reg27*reg25; reg55=reg54+reg55; reg54=reg59*reg90; reg65=reg71+reg65; reg71=reg24*reg43;
    T reg120=reg5*reg57; reg75=reg75/reg30; T reg121=reg53*reg76; T reg122=reg59*reg98; T reg123=reg57*reg43;
    reg74=reg93+reg74; reg93=reg88*reg6; T reg124=reg56*reg92; T reg125=reg53*reg66; T reg126=reg84*reg66;
    reg38=reg95+reg38; reg95=reg85*reg92; reg103=reg94+reg103; reg94=reg29*reg9; T reg127=reg56*reg76;
    T reg128=reg58*reg98; reg67=reg82-reg67; reg70=reg97+reg70; reg82=reg84*reg62; reg97=reg56*reg73;
    T reg129=reg58*reg37; T reg130=reg85*reg80; T reg131=reg102*reg26; T reg132=reg99*reg28; T reg133=reg58*reg90;
    T reg134=reg56*reg80; T reg135=reg53*reg62; reg51=reg91+reg51; reg91=reg58*reg8; reg116=reg87-reg116;
    reg87=reg102*reg12; T reg136=reg99*reg46; reg107=reg133+reg107; reg133=reg105*reg37; reg32=reg54+reg32;
    reg54=reg78*reg70; T reg137=reg123*reg38; reg118=reg72+reg118; reg72=reg64*reg73; T reg138=reg85*reg75;
    T reg139=reg84*reg75; T reg140=reg79*reg70; T reg141=reg55*reg75; reg93=reg67-reg93; reg67=reg79*reg74;
    reg97=reg129+reg97; reg28=reg119*reg28; reg129=reg88*reg7; reg18=reg119*reg18; T reg142=reg56*reg51;
    reg10=reg117*reg10; T reg143=reg50*reg4; reg109=reg108+reg109; reg108=reg58*reg65; reg110=reg111+reg110;
    reg111=reg78*reg74; reg106=reg86+reg106; reg86=reg123*reg70; reg26=reg117*reg26; T reg144=reg50*reg14;
    reg132=reg131+reg132; reg63=reg91+reg63; reg91=reg79*reg38; reg121=reg122+reg121; reg122=reg78*reg38;
    reg131=reg59*reg65; T reg145=reg53*reg51; reg94=reg103+reg94; reg127=reg128+reg127; reg45=reg47+reg45;
    reg13=reg89*reg13; reg115=reg112+reg115; reg47=reg88*reg69; reg103=reg22*reg43; reg112=reg15*reg57;
    reg128=reg81*reg100; reg125=reg124+reg125; reg82=reg130+reg82; reg124=reg36*reg69; reg113=reg114+reg113;
    reg114=reg88*reg100; reg130=reg81*reg69; T reg146=0.5*elem.pos(0)[0]; reg135=reg134+reg135; reg134=reg15*reg25;
    T reg147=reg36*reg100; reg126=reg95+reg126; reg95=reg55*reg2; T reg148=reg5*reg43; T reg149=0.5*elem.pos(1)[0];
    T reg150=reg64*reg2; T reg151=reg24*reg57; reg120=reg71+reg120; reg71=0.5*elem.pos(0)[1]; T reg152=0.5*elem.pos(1)[1];
    reg101=reg101/reg30; reg104=reg68-reg104; reg46=reg119*reg46; reg68=reg49*reg139; reg91=reg127+reg91;
    reg12=reg117*reg12; reg127=reg55*reg139; T reg153=reg55*reg150; reg124=reg82+reg124; reg82=reg49*reg138;
    reg67=reg97+reg67; reg97=reg36*reg75; T reg154=reg148*reg101; T reg155=reg151*reg101; T reg156=reg120*reg101;
    reg137=reg118+reg137; reg118=reg149+reg146; T reg157=reg152-reg71; T reg158=reg64*reg51; T reg159=reg78*reg94;
    reg145=reg131+reg145; reg34=reg33+reg34; reg33=reg105*reg65; reg131=reg17*reg139; reg122=reg121+reg122;
    reg121=reg27*reg89; T reg160=reg17*reg138; reg111=reg110+reg111; reg110=reg15*reg43; T reg161=reg22*reg57;
    reg86=reg106+reg86; reg106=reg55*reg141; reg112=reg103+reg112; reg54=reg32+reg54; reg32=reg17*reg141;
    reg103=reg55*reg95; T reg162=reg49*reg141; reg140=reg107+reg140; reg147=reg126+reg147; reg40=reg43*reg40;
    reg41=reg57*reg41; reg107=reg79*reg94; reg142=reg108+reg142; reg108=reg22*reg25; reg128=reg125+reg128;
    reg125=reg52*reg150; reg130=reg135+reg130; reg14=reg13*reg14; reg28=reg26+reg28; reg144=reg132+reg144;
    reg26=reg53*reg93; reg126=reg56*reg116; reg63=reg129-reg63; reg129=reg60*reg95; reg114=reg113+reg114;
    reg113=reg123*reg74; reg132=reg60*reg150; reg47=reg115+reg47; reg72=reg133+reg72; reg4=reg13*reg4;
    reg18=reg10+reg18; reg143=reg109+reg143; reg10=reg59*reg93; reg109=reg58*reg116; reg115=reg45*reg1;
    reg133=reg134*reg1; reg30=reg104/reg30; reg104=reg50*reg9; reg136=reg87+reg136; reg87=0.5*elem.pos(2)[0];
    reg146=reg149-reg146; reg135=reg52*reg95; reg71=reg152+reg71; reg149=0.5*elem.pos(2)[1]; reg152=reg88*reg63;
    T reg163=(*f.m).alpha_1*reg59; T reg164=reg120*reg155; reg135=reg128+reg135; reg128=reg25*reg89; T reg165=reg99*reg115;
    reg4=reg18+reg4; reg106=reg86+reg106; reg18=reg120*reg156; reg146=reg146+reg87; reg86=reg53*reg144;
    reg132=reg47+reg132; reg47=(*f.m).alpha_2*reg56; T reg166=reg102*reg133; T reg167=reg85*reg116; T reg168=reg45*reg115;
    T reg169=reg17*reg97; reg159=reg145+reg159; reg145=reg105*reg2; reg104=reg136+reg104; reg136=0.5*elem.pos(3)[1];
    reg158=reg33+reg158; reg33=reg123*reg94; T reg170=reg108*reg0; T reg171=reg41*reg155; reg131=reg122+reg131;
    reg122=reg34*reg0; T reg172=(*f.m).alpha_1*reg58; T reg173=0.5*elem.pos(3)[0]; reg46=reg12+reg46; reg10=reg109+reg10;
    reg12=reg27*reg43; reg109=reg5*reg89; T reg174=(*f.m).alpha_2*reg53; T reg175=reg41*reg154; reg160=reg111+reg160;
    reg9=reg13*reg9; reg111=reg49*reg97; reg107=reg142+reg107; reg118=reg87-reg118; reg26=reg126+reg26;
    reg5=reg5*reg27; reg87=reg81*reg63; reg71=reg149-reg71; reg126=reg40*reg155; reg68=reg91+reg68;
    reg91=reg59*reg143; reg142=reg110*reg30; reg127=reg137+reg127; reg14=reg28+reg14; reg28=reg45*reg133;
    reg153=reg124+reg153; reg124=reg161*reg30; reg137=reg40*reg154; reg82=reg67+reg82; reg125=reg130+reg125;
    reg67=reg99*reg133; reg162=reg140+reg162; reg130=reg84*reg93; reg157=reg149+reg157; reg129=reg114+reg129;
    reg114=reg112*reg30; reg140=reg40*reg156; reg113=reg72+reg113; reg23=reg43*reg23; reg32=reg54+reg32;
    reg103=reg147+reg103; reg54=reg102*reg115; reg72=reg41*reg156; reg147=reg55*reg138; reg96=reg57*reg96;
    reg149=reg58*reg143; T reg176=reg56*reg144; T reg177=reg121*reg101; T reg178=0.5*vectors[0][indices[0]+1]; T reg179=(*f.m).alpha_2*reg81;
    T reg180=0.5*vectors[0][indices[1]+0]; T reg181=(*f.m).alpha_3*reg78; reg174=reg163+reg174; reg163=0.5*vectors[0][indices[0]+0]; T reg182=reg112*reg124;
    T reg183=(*f.m).alpha_1*reg88; T reg184=reg120*reg154; reg147=reg113+reg147; reg164=reg127+reg164; reg113=0.5*vectors[0][indices[1]+1];
    reg127=reg5*reg1; reg152=reg10+reg152; reg10=reg60*reg145; reg166=reg132+reg166; reg132=reg117*reg170;
    reg54=reg129+reg54; reg129=reg117*reg122; reg87=reg26+reg87; reg26=reg52*reg145; reg67=reg125+reg67;
    reg125=reg119*reg170; reg165=reg135+reg165; reg135=reg119*reg122; reg130=reg167+reg130; reg167=reg36*reg63;
    reg28=reg153+reg28; reg153=reg34*reg170; reg168=reg103+reg168; reg103=reg34*reg122; T reg185=reg64*reg144;
    T reg186=reg105*reg143; T reg187=reg112*reg114; reg18=reg106+reg18; reg47=reg172+reg47; reg106=(*f.m).alpha_3*reg79;
    reg172=reg55*reg97; reg33=reg158+reg33; reg157=reg157-reg136; reg118=reg173+reg118; reg175=reg160+reg175;
    reg158=reg96*reg142; reg171=reg131+reg171; reg131=reg96*reg124; reg169=reg159+reg169; reg159=reg41*reg177;
    reg72=reg32+reg72; reg32=reg96*reg114; reg86=reg91+reg86; reg91=reg78*reg104; reg160=reg24*reg27;
    T reg188=reg59*reg4; T reg189=reg53*reg14; reg15=reg15*reg89; T reg190=reg25*reg43; reg126=reg68+reg126;
    reg68=reg23*reg124; reg111=reg107+reg111; reg107=reg40*reg177; T reg191=reg23*reg142; reg140=reg162+reg140;
    reg162=reg23*reg114; reg176=reg149+reg176; reg149=reg79*reg104; reg137=reg82+reg137; reg82=reg128*reg30;
    T reg192=reg45*reg75; reg9=reg46+reg9; reg46=reg58*reg4; reg71=reg136+reg71; reg173=reg146-reg173;
    reg136=reg56*reg14; reg24=reg24*reg89; reg27=reg27*reg57; reg109=reg12+reg109; reg10=reg152+reg10;
    reg12=reg102*reg127; reg132=reg166+reg132; reg129=reg54+reg129; reg26=reg87+reg26; reg54=reg99*reg127;
    reg87=reg64*reg14; reg146=reg105*reg4; reg125=reg67+reg125; reg135=reg165+reg135; reg67=reg71*reg173;
    reg152=reg123*reg104; reg185=reg186+reg185; reg136=reg46+reg136; reg187=reg18+reg187; reg18=reg79*reg9;
    reg167=reg130+reg167; reg46=reg55*reg145; reg130=reg78*reg9; reg68=reg126+reg68; reg191=reg137+reg191;
    reg126=reg113-reg178; reg107=reg111+reg107; reg111=reg23*reg82; reg178=reg113+reg178; reg113=0.5*vectors[0][indices[2]+1];
    reg162=reg140+reg162; reg137=reg180+reg163; reg149=reg176+reg149; reg140=reg49*reg192; reg165=reg109*reg101;
    reg75=reg34*reg75; reg166=0.5*vectors[0][indices[2]+0]; reg163=reg180-reg163; reg176=reg160*reg0; reg182=reg164+reg182;
    reg159=reg169+reg159; reg164=reg96*reg82; reg64=(*f.m).alpha_2*reg64; reg105=(*f.m).alpha_1*reg105; reg32=reg72+reg32;
    reg72=(*f.m).alpha_3*reg29; reg91=reg86+reg91; reg86=reg17*reg192; reg153=reg28+reg153; reg28=reg112*reg142;
    reg184=reg147+reg184; reg179=reg183+reg179; reg22=reg22*reg89; reg181=reg174+reg181; reg189=reg188+reg189;
    reg147=reg25*reg57; reg15=reg190+reg15; reg103=reg168+reg103; reg106=reg47+reg106; reg24=reg27+reg24;
    reg27=reg120*reg177; reg172=reg33+reg172; reg158=reg175+reg158; reg131=reg171+reg131; reg33=reg157*reg118;
    reg33=reg67-reg33; reg47=reg132*reg103; reg67=reg129*reg153; reg137=reg166-reg137; reg43=reg89*reg43;
    reg168=0.5*vectors[0][indices[3]+1]; reg22=reg147+reg22; reg147=reg135*reg153; reg46=reg167+reg46; reg167=reg45*reg127;
    reg178=reg113-reg178; reg169=reg162*reg182; reg171=reg32*reg182; reg174=reg68*reg187; reg126=reg113+reg126;
    reg113=reg131*reg187; reg54=reg26+reg54; reg64=reg105+reg64; reg26=reg119*reg176; reg105=(*f.m).alpha_3*reg123;
    reg175=reg158*reg106; reg180=reg125*reg103; reg183=reg68*reg181; reg186=reg117*reg176; reg12=reg10+reg12;
    reg10=(*f.m).alpha_1*reg5; reg72=reg179+reg72; reg179=reg131*reg181; reg188=(*f.m).alpha_2*reg134; reg190=reg191*reg106;
    reg101=reg24*reg101; reg166=reg163+reg166; reg163=reg15*reg30; T reg193=0.5*vectors[0][indices[3]+0]; reg49=reg49*reg75;
    reg164=reg159+reg164; reg159=reg55*reg192; T reg194=reg123*reg9; reg18=reg136+reg18; reg136=reg41*reg165;
    reg152=reg185+reg152; reg86=reg91+reg86; reg27=reg172+reg27; reg91=reg112*reg82; reg87=reg146+reg87;
    reg17=reg17*reg75; reg146=reg40*reg165; reg28=reg184+reg28; reg140=reg149+reg140; reg111=reg107+reg111;
    reg130=reg189+reg130; reg167=reg46+reg167; reg126=reg126-reg168; reg46=reg68*reg32; reg107=reg28*reg106;
    reg71=reg71/reg33; reg149=reg162*reg131; reg172=reg34*reg176; reg184=(*f.m).alpha_1*reg160; reg157=reg157/reg33;
    reg186=reg12+reg186; reg12=reg96*reg163; reg185=(*f.m).alpha_2*reg108; reg136=reg86+reg136; reg173=reg173/reg33;
    reg179=reg175+reg179; reg57=reg89*reg57; reg26=reg54+reg26; reg164=reg164*reg72; reg54=reg55*reg75;
    reg33=reg118/reg33; reg91=reg27+reg91; reg183=reg190+reg183; reg30=reg22*reg30; reg146=reg140+reg146;
    reg27=reg23*reg163; reg111=reg111*reg72; reg105=reg64+reg105; reg64=reg129*reg125; reg86=reg132*reg135;
    reg67=reg47-reg67; reg171=reg113-reg171; reg166=reg166-reg193; reg147=reg180-reg147; reg159=reg152+reg159;
    reg137=reg193+reg137; reg169=reg174-reg169; reg47=reg120*reg165; reg49=reg18+reg49; reg18=reg41*reg101;
    reg17=reg130+reg17; reg188=reg10+reg188; reg10=(*f.m).alpha_3*reg43; reg113=reg40*reg101; reg178=reg168+reg178;
    reg194=reg87+reg194; reg87=reg182*reg181; reg118=reg112*reg163; reg130=(*f.m).alpha_3*reg57; reg185=reg184+reg185;
    reg140=reg23*reg30; reg172=reg167+reg172; reg10=reg188+reg10; reg47=reg159+reg47; reg113=reg49+reg113;
    reg49=reg26*reg67; reg152=reg33*reg166; reg159=reg71*reg126; reg149=reg46-reg149; reg46=reg157*reg178;
    reg27=reg146+reg27; reg64=reg86-reg64; reg86=reg191*reg171; reg146=reg158*reg169; reg167=0.21132486540518713447*elem.pos(1)[0];
    reg168=0.21132486540518713447*elem.pos(0)[0]; reg174=reg96*reg30; reg18=reg17+reg18; reg17=0.78867513459481286553*elem.pos(1)[1]; reg175=0.21132486540518713447*elem.pos(0)[1];
    reg91=reg91*reg72; reg87=reg107+reg87; reg107=0.21132486540518713447*elem.pos(1)[1]; reg180=0.78867513459481286553*elem.pos(1)[0]; reg184=reg173*reg137;
    reg111=reg183+reg111; reg12=reg136+reg12; reg136=reg120*reg101; reg54=reg194+reg54; reg183=reg162*reg105;
    reg188=reg32*reg105; reg189=0.78867513459481286553*elem.pos(0)[1]; reg190=reg147*reg186; reg193=0.78867513459481286553*elem.pos(0)[0]; reg164=reg179+reg164;
    reg179=reg129*reg172; reg194=reg186*reg103; reg152=reg184-reg152; reg184=reg28*reg149; reg166=reg71*reg166;
    reg174=reg18+reg174; reg137=reg157*reg137; reg140=reg113+reg140; reg178=reg173*reg178; reg126=reg33*reg126;
    reg146=reg86-reg146; reg183=reg111+reg183; reg27=reg27*reg10; reg49=reg190-reg49; reg130=reg185+reg130;
    reg18=reg167+reg193; reg86=0.21132486540518713447*elem.pos(2)[0]; reg111=reg107+reg189; reg113=0.21132486540518713447*elem.pos(2)[1]; reg188=reg164+reg188;
    reg12=reg12*reg10; reg136=reg54+reg136; reg54=reg112*reg30; reg164=reg187*reg105; reg167=reg167-reg168;
    reg185=0.78867513459481286553*elem.pos(2)[0]; reg168=reg168+reg180; reg91=reg87+reg91; reg107=reg107-reg175; reg175=reg17+reg175;
    reg87=0.78867513459481286553*elem.pos(2)[1]; reg46=reg159-reg46; reg159=reg26*reg103; reg190=reg158*reg187; T reg195=reg135*reg172;
    T reg196=reg32*reg28; T reg197=reg191*reg187; T reg198=reg162*reg28; reg118=reg47+reg118; reg47=reg172*reg64;
    T reg199=reg186*reg153; T reg200=0.21132486540518713447*elem.pos(3)[1]; T reg201=reg79*reg80; reg137=reg166-reg137; elem.epsilon[0][0]=reg137;
    reg166=0.78867513459481286553*elem.pos(3)[0]; reg167=reg167+reg185; T reg202=reg79*reg116; reg179=reg194-reg179; reg194=reg132*reg172;
    T reg203=reg78*reg62; T reg204=reg78*reg93; T reg205=reg125*reg172; T reg206=reg186*reg135; reg118=reg118*reg10;
    reg164=reg91+reg164; reg91=reg129*reg26; reg195=reg159-reg195; reg46=reg152+reg46; reg152=reg26*reg153;
    reg159=reg158*reg182; reg196=reg190-reg196; reg190=(*f.m).deltaT*reg181; T reg207=reg131*reg28; T reg208=(*f.m).deltaT*reg106;
    T reg209=reg191*reg182; reg198=reg197-reg198; reg197=reg68*reg28; T reg210=reg191*reg32; T reg211=reg162*reg158;
    reg189=reg17-reg189; reg193=reg180-reg193; reg184=reg146+reg184; reg175=reg87-reg175; reg107=reg87+reg107;
    reg54=reg136+reg54; reg17=0.78867513459481286553*elem.pos(3)[1]; reg168=reg185-reg168; reg126=reg178-reg126; elem.epsilon[0][1]=reg126;
    reg87=0.21132486540518713447*elem.pos(3)[0]; reg47=reg49+reg47; reg174=reg174*reg130; reg12=reg188+reg12; reg111=reg113-reg111;
    reg18=reg86-reg18; reg140=reg140*reg130; reg27=reg183+reg27; reg195=reg195/reg47; reg205=reg152-reg205;
    reg49=0.21132486540518713447*PNODE(1).dep[0]; reg136=0.21132486540518713447*PNODE(1).dep[1]; reg146=0.21132486540518713447*PNODE(0).dep[1]; reg152=0.78867513459481286553*PNODE(1).dep[1]; reg178=0.78867513459481286553*PNODE(0).dep[0];
    reg180=0.78867513459481286553*PNODE(0).dep[1]; reg54=reg54*reg130; reg118=reg164+reg118; reg174=reg12+reg174; reg140=reg27+reg140;
    reg12=reg68*reg158; reg211=reg210-reg211; reg27=reg191*reg131; reg197=reg209-reg197; reg198=reg198/reg184;
    reg169=reg169/reg184; reg207=reg159-reg207; reg196=reg196/reg184; reg171=reg171/reg184; reg159=reg132*reg26;
    reg91=reg206-reg91; reg164=reg186*reg125; reg194=reg199-reg194; reg179=reg179/reg47; reg67=reg67/reg47;
    reg189=reg113+reg189; reg193=reg86+reg193; reg18=reg166+reg18; reg111=reg17+reg111; reg46=0.5*reg46;
    elem.epsilon[0][2]=reg46; reg86=reg126-reg190; reg168=reg168+reg87; reg17=reg107-reg17; reg204=reg202+reg204;
    reg107=reg29*reg63; reg175=reg200+reg175; reg113=(*f.m).deltaT*reg105; reg166=reg167-reg166; reg167=0.78867513459481286553*PNODE(1).dep[0];
    reg183=reg78*reg66; reg185=0.21132486540518713447*PNODE(0).dep[0]; reg188=reg137-reg208; reg147=reg147/reg47; reg199=reg79*reg92;
    reg203=reg201+reg203; reg201=reg29*reg69; reg202=reg188*reg195; reg206=0.78867513459481286553*PNODE(2).dep[0]; reg209=reg17*reg168;
    reg210=reg49-reg185; T reg212=reg166*reg111; T reg213=reg77*reg145; T reg214=reg17*reg18; reg87=reg193-reg87;
    reg12=reg27-reg12; reg200=reg189-reg200; reg211=reg211/reg184; reg149=reg149/reg184; reg197=reg197/reg184;
    reg27=0.21132486540518713447*PNODE(2).dep[1]; reg198=reg198*reg174; reg196=reg196*reg140; reg169=reg169*reg174; reg171=reg171*reg140;
    reg189=reg136+reg180; reg54=reg118+reg54; reg118=reg188*reg147; reg49=reg49+reg178; reg185=reg167+reg185;
    reg193=reg166*reg175; T reg215=0.21132486540518713447*PNODE(2).dep[0]; T reg216=reg146+reg152; T reg217=0.78867513459481286553*PNODE(2).dep[1]; reg146=reg136-reg146;
    reg136=reg46-reg113; T reg218=reg86*reg179; reg194=reg194/reg47; reg64=reg64/reg47; reg183=reg199+reg183;
    reg91=reg91/reg47; reg207=reg207/reg184; reg107=reg204+reg107; reg199=reg29*reg100; reg159=reg164-reg159;
    reg205=reg205/reg47; reg164=reg86*reg67; reg204=reg77*reg150; reg201=reg203+reg201; reg180=reg152-reg180;
    reg178=reg167-reg178; reg49=reg215-reg49; reg189=reg27-reg189; reg209=reg193-reg209; reg174=reg197*reg174;
    reg202=reg218-reg202; reg140=reg207*reg140; reg196=reg198-reg196; reg211=reg211*reg54; reg152=reg64*reg136;
    reg169=reg171-reg169; reg164=reg118-reg164; reg149=reg149*reg54; reg47=reg159/reg47; reg118=reg77*reg95;
    reg199=reg183+reg199; reg159=0.78867513459481286553*PNODE(3).dep[0]; reg210=reg206+reg210; reg167=0.21132486540518713447*PNODE(3).dep[0]; reg171=reg91*reg136;
    reg185=reg206-reg185; reg204=reg201+reg204; reg183=reg168*reg200; reg193=reg50*reg133; reg146=reg146+reg217;
    reg184=reg12/reg184; reg12=0.78867513459481286553*PNODE(3).dep[1]; reg197=reg50*reg127; reg198=reg175*reg87; reg201=0.21132486540518713447*PNODE(3).dep[1];
    reg188=reg188*reg205; reg216=reg217-reg216; reg214=reg212-reg214; reg86=reg86*reg194; reg213=reg107+reg213;
    reg210=reg210-reg159; reg178=reg215+reg178; reg189=reg12+reg189; reg49=reg159+reg49; reg185=reg185+reg167;
    reg174=reg140-reg174; reg54=reg184*reg54; reg216=reg201+reg216; reg193=reg204+reg193; reg211=reg196-reg211;
    reg107=reg13*reg170; reg169=reg149+reg169; reg12=reg146-reg12; reg140=reg17/reg209; reg146=reg175/reg209;
    reg136=reg47*reg136; reg149=reg168/reg209; reg159=reg166/reg209; reg171=reg202-reg171; reg166=reg166/reg214;
    reg184=reg18/reg214; reg196=reg111/reg214; reg152=reg164+reg152; reg17=reg17/reg214; reg164=reg50*reg115;
    reg202=1-(*f.m).resolution; reg203=reg13*reg176; reg86=reg188-reg86; reg180=reg27+reg180; reg118=reg199+reg118;
    reg197=reg213+reg197; reg183=reg198-reg183; reg27=reg18*reg200; reg188=reg111*reg87; reg198=reg166*reg49;
    reg199=reg132*reg171; reg204=reg149*reg210; reg107=reg193+reg107; reg193=reg200/reg183; reg27=reg188-reg27;
    reg188=reg140*reg216; reg175=reg175/reg183; reg206=reg159*reg185; reg207=reg26*reg152; reg168=reg168/reg183;
    reg212=reg87/reg183; reg106=(*f.m).resolution*reg106; reg213=reg146*reg12; reg86=reg136+reg86; reg167=reg178-reg167;
    reg201=reg180-reg201; reg136=reg17*reg189; reg211=reg202*reg211; reg178=reg125*reg171; reg203=reg197+reg203;
    reg169=reg202*reg169; reg180=reg186*reg152; reg164=reg118+reg164; reg118=reg13*reg122; reg181=(*f.m).resolution*reg181;
    reg197=reg184*reg210; reg174=reg54+reg174; reg54=reg196*reg12; reg215=reg168*reg167; reg197=reg198-reg197;
    reg198=reg175*reg201; reg217=reg212*reg185; reg218=reg184*reg12; T reg219=reg166*reg189; T reg220=reg17*reg49;
    T reg221=reg196*reg210; reg136=reg54-reg136; reg54=reg107*reg171; T reg222=reg203*reg152; T reg223=reg193*reg216;
    reg199=reg180+reg199; reg180=reg129*reg86; reg178=reg207+reg178; reg207=reg135*reg86; reg12=reg149*reg12;
    T reg224=reg159*reg216; T reg225=reg140*reg185; reg210=reg146*reg210; reg188=reg213-reg188; reg118=reg164+reg118;
    reg164=reg151*reg62; reg213=reg148*reg80; reg211=reg181+reg211; reg169=reg106+reg169; reg174=reg202*reg174;
    reg106=reg151*reg93; reg181=reg148*reg116; reg200=reg200/reg27; reg111=reg111/reg27; reg18=reg18/reg27;
    reg204=reg206-reg204; reg87=reg87/reg27; reg105=(*f.m).resolution*reg105; reg206=(*f.m).resolution*reg205; T reg226=(*f.m).resolution*reg194;
    reg218=reg219-reg218; reg220=reg221-reg220; reg136=reg197+reg136; reg191=reg191*reg202; reg197=reg87*reg49;
    reg219=reg18*reg167; reg174=reg105+reg174; reg169=(*f.m).deltaT*reg169; reg211=(*f.m).deltaT*reg211; reg105=reg121*reg69;
    reg80=reg110*reg80; reg62=reg161*reg62; reg221=reg111*reg201; T reg227=reg200*reg189; T reg228=(*f.m).resolution*reg195;
    T reg229=(*f.m).resolution*reg147; T reg230=(*f.m).resolution*reg67; T reg231=reg168*reg201; T reg232=(*f.m).resolution*reg179; reg216=reg212*reg216;
    reg185=reg193*reg185; T reg233=reg175*reg167; reg223=reg198-reg223; reg68=reg68*reg202; reg158=reg158*reg202;
    reg131=reg131*reg202; reg215=reg217-reg215; reg28=reg28*reg202; reg182=reg182*reg202; reg93=reg161*reg93;
    reg116=reg110*reg116; reg198=reg151*reg66; reg217=reg148*reg92; reg180=reg199+reg180; reg12=reg224-reg12;
    reg225=reg210-reg225; reg188=reg204+reg188; reg54=reg222+reg54; reg106=reg181+reg106; reg181=reg118*reg86;
    reg199=reg121*reg63; reg164=reg213+reg164; reg207=reg178+reg207; reg178=reg121*reg100; reg12=reg12-reg211;
    reg93=reg116+reg93; reg63=reg128*reg63; reg230=reg68-reg230; reg191=reg229+reg191; reg187=reg187*reg202;
    reg219=reg197-reg219; reg69=reg128*reg69; reg32=reg32*reg202; reg225=reg225-reg169; reg202=reg162*reg202;
    reg207=reg190+reg207; reg223=reg215+reg223; reg188=0.5*reg188; reg185=reg233-reg185; reg72=(*f.m).deltaT*reg72;
    reg227=reg221-reg227; reg231=reg216-reg231; reg167=reg111*reg167; reg49=reg200*reg49; reg181=reg54+reg181;
    reg54=(*f.m).resolution*reg64; reg189=reg87*reg189; reg92=reg110*reg92; reg66=reg161*reg66; reg68=(*f.m).resolution*reg91;
    reg116=(*f.m).resolution*reg47; reg201=reg18*reg201; reg62=reg80+reg62; reg105=reg164+reg105; reg80=reg120*reg150;
    reg174=(*f.m).deltaT*reg174; reg162=reg120*reg145; reg199=reg106+reg199; reg226=reg182-reg226; reg28=reg206+reg28;
    reg136=0.5*reg136; reg131=reg232+reg131; reg220=reg220-reg169; reg228=reg158-reg228; reg218=reg218-reg211;
    reg180=reg208+reg180; reg198=reg217+reg198; reg136=reg136-reg174; reg231=reg231-reg211; reg106=reg226*reg12;
    reg187=reg116+reg187; reg116=reg109*reg127; reg158=reg120*reg95; reg164=reg228*reg225; reg66=reg92+reg66;
    reg100=reg128*reg100; reg181=reg72+reg181; reg92=reg131*reg12; reg227=reg219+reg227; reg201=reg189-reg201;
    reg162=reg199+reg162; reg182=reg131*reg218; reg178=reg198+reg178; reg49=reg167-reg49; reg167=reg228*reg220;
    reg189=reg180+reg207; reg197=reg109*reg133; reg80=reg105+reg80; reg105=reg226*reg218; reg69=reg62+reg69;
    reg145=reg112*reg145; reg62=reg28*reg225; reg198=reg191*reg225; reg202=reg54+reg202; reg63=reg93+reg63;
    reg150=reg112*reg150; reg54=reg28*reg220; reg223=0.5*reg223; reg188=reg188-reg174; reg93=reg191*reg220;
    reg199=reg230*reg12; reg185=reg185-reg169; reg204=reg230*reg218; reg68=reg32-reg68; reg116=reg162+reg116;
    reg227=0.5*reg227; reg32=reg24*reg176; reg189=reg181+reg189; reg49=reg49-reg169; reg201=reg201-reg211;
    reg158=reg178+reg158; reg162=reg109*reg115; reg127=reg15*reg127; reg145=reg63+reg145; reg197=reg80+reg197;
    reg63=reg24*reg170; reg80=reg202*reg188; reg92=reg164+reg92; reg223=reg223-reg174; reg164=reg28*reg185;
    reg178=reg226*reg231; reg133=reg15*reg133; reg100=reg66+reg100; reg95=reg112*reg95; reg150=reg69+reg150;
    reg66=reg68*reg136; reg182=reg167+reg182; reg106=reg62+reg106; reg62=reg228*reg185; reg69=reg131*reg231;
    reg167=reg202*reg136; reg204=reg93+reg204; reg93=reg187*reg188; reg206=reg187*reg136; reg105=reg54+reg105;
    reg199=reg198+reg199; reg54=reg191*reg185; reg198=reg68*reg188; reg210=reg230*reg231; reg213=reg187*reg223;
    reg178=reg164+reg178; reg210=reg54+reg210; reg80=reg199+reg80; reg54=reg172*reg152; reg164=reg153*reg171;
    reg199=reg24*reg122; reg162=reg158+reg162; reg66=reg182+reg66; reg167=reg204+reg167; reg63=reg197+reg63;
    reg206=reg105+reg206; reg127=reg145+reg127; reg176=reg22*reg176; reg198=reg92+reg198; reg32=reg116+reg32;
    reg189=reg189/3; reg227=reg227-reg174; reg115=reg15*reg115; reg92=reg28*reg49; reg105=reg68*reg223;
    reg69=reg62+reg69; reg62=reg226*reg201; reg95=reg100+reg95; reg100=reg191*reg49; reg116=reg230*reg201;
    reg170=reg22*reg170; reg145=reg202*reg223; reg93=reg106+reg93; reg133=reg150+reg133; reg106=reg228*reg49;
    reg150=reg131*reg201; reg198=reg12*reg198; reg122=reg22*reg122; reg12=reg63*reg171; reg176=reg127+reg176;
    reg145=reg210+reg145; reg115=reg95+reg115; reg93=2*reg93; reg180=reg180-reg189; reg105=reg69+reg105;
    reg206=2*reg206; reg69=reg32*reg152; reg213=reg178+reg213; reg80=reg225*reg80; reg167=reg220*reg167;
    reg207=reg207-reg189; reg170=reg133+reg170; reg95=reg202*reg227; reg116=reg100+reg116; reg164=reg54+reg164;
    reg66=reg218*reg66; reg199=reg162+reg199; reg150=reg106+reg150; reg54=reg68*reg227; reg100=reg187*reg227;
    reg62=reg92+reg62; reg92=reg103*reg86; reg93=reg188*reg93; reg180=pow(reg180,2); reg189=reg181-reg189;
    reg12=reg69+reg12; reg69=reg199*reg86; reg207=pow(reg207,2); reg122=reg115+reg122; reg92=reg164+reg92;
    reg66=reg167+reg66; reg206=reg136*reg206; reg152=reg176*reg152; reg171=reg170*reg171; reg213=2*reg213;
    reg145=reg185*reg145; reg54=reg150+reg54; reg95=reg116+reg95; reg100=reg62+reg100; reg105=reg231*reg105;
    reg198=reg80+reg198; reg95=reg49*reg95; reg54=reg201*reg54; reg86=reg122*reg86; reg171=reg152+reg171;
    reg213=reg223*reg213; reg100=2*reg100; reg10=(*f.m).deltaT*reg10; reg105=reg145+reg105; reg207=reg180+reg207;
    reg189=pow(reg189,2); reg66=reg206+reg66; reg137=reg137-reg169; reg92=reg113+reg92; reg198=reg93+reg198;
    reg126=reg126-reg211; reg69=reg12+reg69; reg69=reg10+reg69; reg12=reg228*reg137; reg49=reg131*reg126;
    reg54=reg95+reg54; reg86=reg171+reg86; reg100=reg227*reg100; reg46=reg46-reg174; reg198=reg209*reg198;
    reg130=(*f.m).deltaT*reg130; reg62=reg230*reg126; reg105=reg213+reg105; reg92=pow(reg92,2); reg189=reg207+reg189;
    reg66=reg214*reg66; reg80=reg191*reg137; reg198=0.25*reg198; reg105=reg183*reg105; reg86=reg130+reg86;
    reg93=reg202*reg46; reg54=reg100+reg54; reg49=reg12+reg49; reg12=reg68*reg46; reg137=reg28*reg137;
    reg126=reg226*reg126; reg62=reg80+reg62; reg69=pow(reg69,2); reg92=reg189+reg92; reg66=0.25*reg66;
    reg69=reg92+reg69; reg86=pow(reg86,2); reg105=0.25*reg105; reg54=reg27*reg54; reg66=reg198+reg66;
    reg93=reg62+reg93; elem.sigma[0][0]=reg93; reg12=reg49+reg12; elem.sigma[0][1]=reg12; reg126=reg137+reg126;
    reg46=reg187*reg46; reg105=reg66+reg105; reg54=0.25*reg54; reg86=reg69+reg86; reg46=reg126+reg46;
    elem.sigma[0][2]=reg46; reg49=reg58*reg93; reg62=reg59*reg12; reg66=reg56*reg93; reg69=reg53*reg12;
    reg12=reg84*reg12; reg93=reg85*reg93; reg80=reg55*reg46; reg54=reg105+reg54; reg12=reg93+reg12;
    reg86=pow(reg86,0.5); reg92=reg52*reg46; reg69=reg66+reg69; reg62=reg49+reg62; reg46=reg60*reg46;
    elem.sigma_local[0][1]=reg69+reg92; elem.sigma_von_mises=0.86602540378443859659*reg86; elem.sigma_local[0][0]=reg62+reg46; elem.sigma_local[0][2]=reg12+reg80; elem.ener=reg54/2;

      #undef PNODE
    }
    template<class TE,class TF, class TVEVE> static void after_solve_2(TE &elem,TF &f, TVEVE &vectors,const unsigned *indices) {
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
    reg2=1.0/reg2; T reg3=reg0*reg1; T reg4=pow((*f.m).v2[0],2); T reg5=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg6=pow((*f.m).v2[1],2);
    T reg7=1.0/(*f.m).elastic_modulus_3; T reg8=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg9=pow((*f.m).v1[1],2); T reg10=pow((*f.m).v1[0],2); T reg11=reg2*reg3;
    T reg12=reg5*reg11; T reg13=reg7*reg11; T reg14=pow((*f.m).v2[2],2); reg6=reg4+reg6; reg4=pow((*f.m).v1[2],2);
    reg9=reg10+reg9; reg10=reg8*reg11; T reg15=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg16=1.0/(*f.m).elastic_modulus_2; reg4=reg9+reg4;
    reg9=reg8*reg12; reg14=reg6+reg14; reg6=reg16*reg13; T reg17=reg8*reg10; T reg18=reg15*reg13;
    T reg19=1.0/(*f.m).elastic_modulus_1; T reg20=reg15*reg10; reg9=reg18+reg9; reg4=pow(reg4,0.5); reg17=reg6-reg17;
    reg6=reg16*reg12; reg14=pow(reg14,0.5); T reg21=(*f.m).v1[1]/reg4; T reg22=(*f.m).v1[2]/reg4; T reg23=reg20+reg6;
    T reg24=reg19*reg17; T reg25=reg15*reg9; T reg26=(*f.m).v2[1]/reg14; T reg27=(*f.m).v2[2]/reg14; reg4=(*f.m).v1[0]/reg4;
    T reg28=reg22*reg26; reg14=(*f.m).v2[0]/reg14; T reg29=reg21*reg27; T reg30=reg7*reg3; T reg31=reg5*reg10;
    T reg32=reg16*reg11; T reg33=reg5*reg3; reg3=reg8*reg3; reg13=reg19*reg13; reg25=reg24-reg25;
    reg24=reg5*reg23; T reg34=reg2*reg1; T reg35=reg5*reg12; reg11=reg15*reg11; T reg36=reg8*reg34;
    reg12=reg15*reg12; T reg37=reg5*reg11; T reg38=reg2*reg0; reg31=reg18+reg31; reg18=reg29-reg28;
    reg10=reg19*reg10; reg35=reg13-reg35; reg13=reg4*reg27; T reg39=reg16*reg30; reg30=reg15*reg30;
    T reg40=reg8*reg3; T reg41=reg22*reg14; T reg42=reg8*reg33; T reg43=reg7*reg34; T reg44=2*reg14;
    reg24=reg25-reg24; reg25=2*reg4; T reg45=reg5*reg32; reg34=reg5*reg34; reg31=reg31/reg24;
    reg35=reg35/reg24; reg17=reg17/reg24; T reg46=reg41-reg13; reg45=reg20+reg45; reg9=reg9/reg24;
    T reg47=reg5*reg38; reg40=reg39-reg40; reg42=reg30+reg42; reg30=2*reg18; reg39=reg26*reg44;
    T reg48=pow(reg26,2); T reg49=pow(reg14,2); T reg50=reg16*reg43; T reg51=reg21*reg25; T reg52=pow(reg21,2);
    T reg53=pow(reg4,2); reg33=reg16*reg33; T reg54=reg8*reg38; reg3=reg15*reg3; reg11=reg15*reg11;
    reg12=reg10+reg12; reg37=reg10+reg37; reg32=reg19*reg32; reg10=reg4*reg26; T reg55=reg21*reg14;
    reg43=reg15*reg43; T reg56=reg8*reg36; T reg57=reg8*reg34; reg38=reg7*reg38; T reg58=reg15*reg38;
    T reg59=reg46*reg30; T reg60=reg51*reg31; T reg61=reg3+reg33; T reg62=reg39*reg35; T reg63=reg8*reg54;
    reg42=reg15*reg42; reg40=reg19*reg40; T reg64=reg39*reg9; T reg65=reg8*reg47; reg37=reg37/reg24;
    T reg66=reg51*reg17; reg12=reg12/reg24; reg34=reg16*reg34; reg36=reg15*reg36; reg23=reg23/reg24;
    T reg67=reg49*reg9; T reg68=reg52*reg17; T reg69=reg48*reg9; reg11=reg32-reg11; reg32=reg10-reg55;
    reg56=reg50-reg56; reg50=pow(reg22,2); T reg70=reg48*reg35; T reg71=reg52*reg31; reg57=reg43+reg57;
    reg45=reg45/reg24; reg43=pow(reg27,2); T reg72=pow(reg46,2); T reg73=pow(reg18,2); T reg74=reg53*reg17;
    reg38=reg16*reg38; T reg75=reg53*reg31; T reg76=reg49*reg35; T reg77=reg43*reg35; T reg78=reg72*reg12;
    reg70=reg71+reg70; reg71=reg50*reg17; T reg79=reg43*reg9; T reg80=reg50*reg31; T reg81=reg59*reg12;
    reg62=reg60+reg62; reg60=reg73*reg12; T reg82=reg59*reg23; reg76=reg75+reg76; reg11=reg11/reg24;
    reg64=reg66+reg64; reg56=reg19*reg56; reg57=reg15*reg57; reg66=reg39*reg37; reg75=reg51*reg45;
    T reg83=reg36+reg34; T reg84=pow(reg32,2); T reg85=reg48*reg37; reg61=reg5*reg61; T reg86=reg52*reg45;
    reg42=reg40-reg42; reg63=reg38-reg63; reg38=reg72*reg23; reg47=reg16*reg47; reg65=reg58+reg65;
    reg54=reg15*reg54; reg69=reg68+reg69; reg40=reg53*reg45; reg58=reg49*reg37; reg74=reg67+reg74;
    reg67=reg73*reg23; reg68=reg59*reg11; reg66=reg75+reg66; reg75=reg43*reg37; T reg87=reg50*reg45;
    T reg88=reg72*reg11; reg85=reg86+reg85; reg77=reg80+reg77; reg80=reg84*reg12; reg81=reg62+reg81;
    reg58=reg40+reg58; reg40=reg73*reg11; reg57=reg56-reg57; reg56=reg22*reg25; reg62=2*reg21;
    reg83=reg5*reg83; reg86=reg27*reg44; T reg89=2*reg26; reg61=reg42-reg61; reg63=reg19*reg63;
    reg65=reg15*reg65; reg42=reg54+reg47; reg67=reg74+reg67; reg38=reg69+reg38; reg78=reg70+reg78;
    reg79=reg71+reg79; reg69=reg84*reg23; reg60=reg76+reg60; reg82=reg64+reg82; reg64=reg4*reg21;
    reg70=reg14*reg26; reg68=reg66+reg68; reg66=reg86*reg9; reg71=reg56*reg17; reg74=reg84*reg11;
    reg75=reg87+reg75; reg83=reg57-reg83; reg88=reg85+reg88; reg80=reg77+reg80; reg61=reg61/reg24;
    reg40=reg58+reg40; reg65=reg63-reg65; reg69=reg79+reg69; reg42=reg5*reg42; reg57=reg56*reg31;
    reg58=reg86*reg35; reg63=reg52*reg38; reg76=reg48*reg78; reg77=reg4*reg46; reg79=reg21*reg18;
    reg10=reg55+reg10; reg55=reg21*reg26; reg85=reg4*reg14; reg87=2*reg46; reg30=reg32*reg30;
    T reg90=reg27*reg89; T reg91=reg22*reg62; T reg92=reg70*reg81; T reg93=reg52*reg82; T reg94=reg48*reg81;
    T reg95=reg53*reg67; T reg96=reg49*reg60; T reg97=reg53*reg38; T reg98=reg49*reg78; T reg99=reg64*reg82;
    reg82=reg53*reg82; reg81=reg49*reg81; reg38=reg64*reg38; reg78=reg70*reg78; T reg100=reg48*reg60;
    T reg101=reg26*reg62; T reg102=2*reg22; T reg103=reg46*reg18; T reg104=reg14*reg25; T reg105=reg52*reg67;
    reg83=reg83/reg24; T reg106=reg30*reg12; T reg107=reg49*reg15; reg92=reg99+reg92; reg99=reg52*reg69;
    T reg108=reg48*reg80; reg58=reg57+reg58; reg31=reg91*reg31; reg57=reg48*reg16; reg60=reg70*reg60;
    T reg109=reg104*reg15; reg67=reg64*reg67; reg94=reg93+reg94; reg42=reg65-reg42; reg17=reg91*reg17;
    reg65=reg30*reg23; reg66=reg71+reg66; reg98=reg97+reg98; reg9=reg90*reg9; reg71=reg86*reg37;
    reg93=reg73*reg88; reg97=reg49*reg19; T reg110=reg73*reg40; reg96=reg95+reg96; reg95=reg72*reg40;
    T reg111=reg104*reg19; T reg112=reg72*reg68; T reg113=reg53*reg69; T reg114=reg101*reg15; T reg115=reg49*reg80;
    T reg116=reg101*reg16; T reg117=reg85*reg61; T reg118=reg55*reg61; T reg119=reg14*reg46; T reg120=reg26*reg18;
    T reg121=reg73*reg68; reg74=reg75+reg74; reg76=reg63+reg76; reg63=reg72*reg88; reg88=reg103*reg88;
    reg75=reg10*reg61; T reg122=reg48*reg15; reg77=reg79+reg77; reg78=reg38+reg78; reg100=reg105+reg100;
    reg38=reg21*reg46; reg79=reg4*reg18; reg35=reg90*reg35; reg81=reg82+reg81; reg68=reg103*reg68;
    reg82=reg56*reg45; reg102=reg27*reg102; reg105=reg22*reg27; reg87=reg32*reg87; reg35=reg31+reg35;
    reg31=reg104*reg75; reg121=reg81+reg121; reg12=reg87*reg12; reg9=reg17+reg9; reg23=reg87*reg23;
    reg106=reg58+reg106; reg17=reg10*reg118; reg40=reg103*reg40; reg93=reg98+reg93; reg58=reg104*reg118;
    reg88=reg78+reg88; reg78=reg102*reg5; reg81=reg101*reg8; reg60=reg67+reg60; reg68=reg92+reg68;
    reg67=reg49*reg5; reg92=reg104*reg5; reg115=reg113+reg115; reg98=reg48*reg8; reg80=reg70*reg80;
    reg69=reg64*reg69; reg113=reg73*reg74; reg107=reg57-reg107; reg57=reg72*reg74; reg108=reg99+reg108;
    reg99=reg43*reg8; T reg123=reg43*reg5; T reg124=reg77*reg83; T reg125=reg101*reg117; reg37=reg90*reg37;
    reg62=reg46*reg62; reg71=reg82+reg71; reg82=reg38*reg83; T reg126=reg79*reg83; T reg127=reg22*reg32;
    reg25=reg18*reg25; reg118=reg101*reg118; reg63=reg76+reg63; reg76=reg14*reg18; T reg128=reg26*reg46;
    T reg129=reg105*reg61; reg119=reg120+reg119; reg120=reg101*reg75; reg75=reg10*reg75; reg110=reg96+reg110;
    reg112=reg94+reg112; reg94=reg104*reg117; reg96=reg30*reg11; reg109=reg116-reg109; reg116=reg102*reg8;
    reg114=reg111-reg114; reg45=reg91*reg45; reg95=reg100+reg95; reg65=reg66+reg65; reg24=reg42/reg24;
    reg122=reg97-reg122; reg42=reg25*reg126; reg94=reg110+reg94; reg11=reg87*reg11; reg37=reg45+reg37;
    reg31=reg121+reg31; reg45=reg25*reg124; reg58=reg93+reg58; reg66=reg25*reg82; reg123=reg122-reg123;
    reg93=reg53*reg65; reg97=reg49*reg106; reg100=reg104*reg129; reg110=reg127*reg83; reg19=reg53*reg19;
    reg96=reg71+reg96; reg71=reg119*reg24; reg111=reg76*reg24; reg121=reg128*reg24; reg113=reg115+reg113;
    reg99=reg107-reg99; reg107=reg101*reg129; reg116=reg109-reg116; reg109=(*f.m).alpha_2*reg49; reg115=(*f.m).alpha_1*reg52;
    reg122=(*f.m).alpha_2*reg48; reg120=reg112+reg120; reg112=reg62*reg124; T reg130=reg52*reg65; T reg131=reg48*reg106;
    reg74=reg103*reg74; reg80=reg69+reg80; reg69=reg52*reg15; T reg132=reg43*reg7; reg98=reg67+reg98;
    reg102=reg102*reg7; reg81=reg92+reg81; reg40=reg60+reg40; reg117=reg10*reg117; reg60=(*f.m).alpha_1*reg53;
    reg17=reg88+reg17; reg67=reg77*reg82; reg88=reg27*reg32; reg92=reg62*reg126; reg125=reg95+reg125;
    reg44=reg18*reg44; reg89=reg46*reg89; reg124=reg77*reg124; reg75=reg68+reg75; reg78=reg114-reg78;
    reg68=reg4*reg32; reg95=reg22*reg18; reg57=reg108+reg57; reg15=reg53*reg15; reg118=reg63+reg118;
    reg41=reg13+reg41; reg16=reg52*reg16; reg82=reg62*reg82; reg23=reg9+reg23; reg12=reg35+reg12;
    reg92=reg125+reg92; reg9=reg53*reg123; reg13=reg52*reg99; reg35=reg119*reg121; reg67=reg17+reg67;
    reg97=reg93+reg97; reg17=reg44*reg71; reg117=reg40+reg117; reg126=reg77*reg126; reg40=reg89*reg121;
    reg107=reg57+reg107; reg57=reg62*reg110; reg82=reg118+reg82; reg63=reg49*reg78; reg112=reg120+reg112;
    reg93=reg89*reg71; reg108=reg49*reg12; reg114=reg53*reg23; reg118=reg48*reg116; reg131=reg130+reg131;
    reg120=reg72*reg96; reg129=reg10*reg129; reg74=reg80+reg74; reg69=reg19-reg69; reg19=reg50*reg5;
    reg80=reg52*reg23; reg125=reg48*reg12; reg130=reg89*reg111; T reg133=reg73*reg96; reg98=reg132-reg98;
    reg132=reg50*reg8; reg8=reg52*reg8; reg124=reg75+reg124; reg5=reg53*reg5; reg75=reg14*reg32;
    T reg134=reg27*reg18; T reg135=reg21*reg32; T reg136=reg22*reg46; reg68=reg95+reg68; reg29=reg28+reg29;
    reg15=reg16-reg15; reg16=reg85*reg123; reg28=reg55*reg99; reg71=reg119*reg71; reg95=reg85*reg78;
    T reg137=reg55*reg116; T reg138=(*f.m).alpha_2*reg43; T reg139=(*f.m).alpha_1*reg50; reg65=reg64*reg65; reg106=reg70*reg106;
    reg109=reg60+reg109; reg60=(*f.m).alpha_3*reg73; reg122=reg115+reg122; reg115=(*f.m).alpha_3*reg72; reg45=reg31+reg45;
    reg31=reg52*reg116; T reg140=reg25*reg110; reg100=reg113+reg100; reg113=reg53*reg78; T reg141=reg49*reg123;
    reg121=reg44*reg121; reg66=reg58+reg66; reg58=reg44*reg111; reg42=reg94+reg42; reg81=reg102-reg81;
    reg11=reg37+reg11; reg37=reg88*reg24; reg94=reg41*reg61; reg102=reg48*reg99; T reg142=reg50*reg81;
    reg115=reg122+reg115; reg31=reg113+reg31; reg113=reg43*reg98; reg8=reg5+reg8; reg7=reg50*reg7;
    reg93=reg112+reg93; reg102=reg141+reg102; reg138=reg139+reg138; reg5=(*f.m).alpha_3*reg84; reg35=reg67+reg35;
    reg111=reg119*reg111; reg106=reg65+reg106; reg126=reg117+reg126; reg120=reg131+reg120; reg65=reg101*reg94;
    reg67=(*f.m).alpha_1*reg64; reg112=reg72*reg11; reg117=(*f.m).alpha_2*reg70; reg122=reg10*reg2; reg125=reg80+reg125;
    reg80=reg50*reg98; reg131=reg70*reg2; reg12=reg70*reg12; reg96=reg103*reg96; reg23=reg64*reg23;
    reg4=reg4*reg22; reg14=reg14*reg27; reg70=reg43*reg81; reg118=reg63+reg118; reg19=reg69-reg19;
    reg61=reg29*reg61; reg63=reg68*reg83; reg110=reg77*reg110; reg69=reg73*reg11; reg108=reg114+reg108;
    reg58=reg42+reg58; reg13=reg9+reg13; reg9=reg104*reg94; reg133=reg97+reg133; reg121=reg66+reg121;
    reg140=reg100+reg140; reg17=reg45+reg17; reg42=reg44*reg37; reg57=reg107+reg57; reg45=reg89*reg37;
    reg137=reg95+reg137; reg71=reg124+reg71; reg129=reg74+reg129; reg66=reg105*reg98; reg28=reg16+reg28;
    reg16=reg105*reg81; reg40=reg82+reg40; reg132=reg15-reg132; reg135=reg136+reg135; reg60=reg109+reg60;
    reg130=reg92+reg130; reg75=reg134+reg75; reg15=reg27*reg46; reg74=reg26*reg32; reg82=reg49*reg19;
    reg92=reg48*reg132; reg142=reg31+reg142; reg31=reg51*reg122; reg95=reg51*reg131; reg37=reg119*reg37;
    reg11=reg103*reg11; reg12=reg23+reg12; reg110=reg129+reg110; reg96=reg106+reg96; reg94=reg10*reg94;
    reg80=reg13+reg80; reg13=reg10*reg131; reg16=reg137+reg16; reg23=reg10*reg122; reg5=reg138+reg5;
    reg117=reg67+reg117; reg103=(*f.m).alpha_3*reg103; reg67=(*f.m).alpha_1*reg4; reg97=(*f.m).alpha_2*reg14; reg8=reg7-reg8;
    reg7=reg40*reg71; reg100=reg121*reg71; reg106=reg93*reg35; reg107=reg17*reg35; reg109=reg58*reg60;
    reg114=reg121*reg115; reg14=reg14*reg0; reg124=reg130*reg60; reg129=reg40*reg115; reg134=reg41*reg0;
    reg136=reg53*reg19; reg137=reg52*reg132; reg111=reg126+reg111; reg113=reg102+reg113; reg102=reg39*reg131;
    reg101=reg101*reg61; reg112=reg125+reg112; reg125=reg62*reg63; reg65=reg120+reg65; reg45=reg57+reg45;
    reg104=reg104*reg61; reg69=reg108+reg69; reg57=reg25*reg63; reg9=reg133+reg9; reg42=reg140+reg42;
    reg108=reg75*reg24; reg83=reg135*reg83; reg70=reg118+reg70; reg118=reg39*reg122; reg27=reg26*reg27;
    reg66=reg28+reg66; reg22=reg21*reg22; reg18=reg32*reg18; reg74=reg15+reg74; reg50=reg50*reg8;
    reg137=reg136+reg137; reg15=reg29*reg1; reg21=reg27*reg1; reg26=reg56*reg14; reg95=reg80+reg95;
    reg2=reg64*reg2; reg57=reg9+reg57; reg24=reg74*reg24; reg118=reg70+reg118; reg9=reg86*reg134;
    reg46=reg32*reg46; reg28=reg85*reg19; reg32=reg55*reg132; reg13=reg66+reg13; reg64=reg41*reg14;
    reg23=reg16+reg23; reg16=reg41*reg134; reg103=reg117+reg103; reg97=reg67+reg97; reg18=(*f.m).alpha_3*reg18;
    reg66=(*f.m).alpha_1*reg22; reg27=(*f.m).alpha_2*reg27; reg106=reg7-reg106; reg107=reg100-reg107; reg7=reg121*reg93;
    reg67=reg17*reg40; reg114=reg109+reg114; reg42=reg42*reg5; reg129=reg124+reg129; reg45=reg45*reg5;
    reg70=reg111*reg60; reg80=reg35*reg115; reg31=reg142+reg31; reg100=reg56*reg134; reg92=reg82+reg92;
    reg43=reg43*reg8; reg61=reg10*reg61; reg11=reg12+reg11; reg63=reg77*reg63; reg94=reg96+reg94;
    reg37=reg110+reg37; reg102=reg113+reg102; reg62=reg62*reg83; reg101=reg112+reg101; reg12=reg86*reg14;
    reg82=reg44*reg108; reg104=reg69+reg104; reg25=reg25*reg83; reg69=reg89*reg108; reg125=reg65+reg125;
    reg43=reg92+reg43; reg67=reg7-reg67; reg7=reg90*reg21; reg42=reg114+reg42; reg65=reg51*reg2;
    reg92=reg17*reg103; reg96=reg91*reg15; reg50=reg137+reg50; reg100=reg31+reg100; reg45=reg129+reg45;
    reg31=reg93*reg103; reg80=reg70+reg80; reg37=reg37*reg5; reg64=reg13+reg64; reg13=reg29*reg21;
    reg105=reg105*reg8; reg32=reg28+reg32; reg16=reg23+reg16; reg23=reg29*reg15; reg28=reg90*reg15;
    reg18=reg97+reg18; reg12=reg102+reg12; reg9=reg118+reg9; reg27=reg66+reg27; reg46=(*f.m).alpha_3*reg46;
    reg66=reg58*reg106; reg70=reg39*reg2; reg97=reg130*reg107; reg26=reg95+reg26; reg95=reg91*reg21;
    reg83=reg77*reg83; reg61=reg11+reg61; reg108=reg119*reg108; reg63=reg94+reg63; reg89=reg89*reg24;
    reg62=reg101+reg62; reg69=reg125+reg69; reg44=reg44*reg24; reg25=reg104+reg25; reg82=reg57+reg82;
    reg0=reg4*reg0; reg37=reg80+reg37; reg70=reg43+reg70; reg86=reg86*reg0; reg4=reg71*reg103;
    reg108=reg63+reg108; reg11=1-var_inter[1]; reg89=reg62+reg89; reg69=reg69*reg18; reg31=reg45+reg31;
    reg44=reg25+reg44; reg25=1-var_inter[0]; reg7=reg12+reg7; reg92=reg42+reg92; reg82=reg82*reg18;
    reg12=reg10*reg2; reg28=reg9+reg28; reg105=reg32+reg105; reg1=reg22*reg1; reg65=reg50+reg65;
    reg56=reg56*reg0; reg46=reg27+reg46; reg97=reg66-reg97; reg9=reg111*reg67; reg22=reg130*reg71;
    reg27=reg93*reg111; reg95=reg26+reg95; reg26=reg58*reg71; reg83=reg61+reg83; reg24=reg119*reg24;
    reg13=reg64+reg13; reg96=reg100+reg96; reg23=reg16+reg23; reg16=reg17*reg111; reg32=reg58*reg35;
    reg42=reg40*reg111; reg43=elem.pos(1)[1]*var_inter[0]; reg45=reg25*elem.pos(0)[1]; reg89=reg89*reg46; reg69=reg31+reg69;
    reg31=reg25*elem.pos(0)[0]; reg50=elem.pos(1)[0]*var_inter[0]; reg57=reg11*elem.pos(1)[0]; reg61=reg11*elem.pos(0)[1]; reg62=reg11*elem.pos(1)[1];
    reg9=reg97+reg9; reg63=reg17*reg130; reg16=reg26-reg16; reg12=reg105+reg12; reg41=reg41*reg0;
    reg44=reg44*reg46; reg82=reg92+reg82; reg26=reg11*elem.pos(0)[0]; reg64=reg58*reg93; reg27=reg22-reg27;
    reg22=reg121*reg111; reg66=reg130*reg35; reg24=reg83+reg24; reg86=reg70+reg86; reg90=reg90*reg1;
    reg70=reg96*reg13; reg80=reg28*reg13; reg83=reg7*reg23; reg92=reg95*reg23; reg91=reg91*reg1;
    reg56=reg65+reg56; reg4=reg37+reg4; reg108=reg108*reg18; reg42=reg66-reg42; reg37=reg58*reg40;
    reg41=reg12+reg41; reg29=reg29*reg1; reg63=reg64-reg63; reg12=var_inter[1]*elem.pos(2)[0]; reg44=reg82+reg44;
    reg80=reg83-reg80; reg57=reg57-reg26; reg90=reg86+reg90; reg106=reg106/reg9; reg108=reg4+reg108;
    reg27=reg27/reg9; reg4=reg121*reg130; reg64=reg43+reg45; reg65=elem.pos(2)[1]*var_inter[0]; reg22=reg32-reg22;
    reg32=reg50+reg31; reg66=elem.pos(2)[0]*var_inter[0]; reg89=reg69+reg89; reg24=reg24*reg46; reg91=reg56+reg91;
    reg70=reg92-reg70; reg56=reg95*reg28; reg69=reg96*reg7; reg16=reg16/reg9; reg82=var_inter[1]*elem.pos(2)[1];
    reg62=reg62-reg61; reg107=reg107/reg9; reg42=reg42/reg9; reg22=reg22/reg9; reg107=reg107*reg89;
    reg27=reg27*reg44; reg16=reg16*reg89; reg83=elem.pos(3)[0]*reg25; reg12=reg57+reg12; reg57=elem.pos(3)[1]*reg25;
    reg86=var_inter[1]*elem.pos(3)[0]; reg92=var_inter[1]*elem.pos(3)[1]; reg69=reg56-reg69; reg82=reg62+reg82; reg66=reg66-reg32;
    reg65=reg65-reg64; reg24=reg108+reg24; reg29=reg41+reg29; reg41=reg80*reg91; reg56=reg90*reg70;
    reg4=reg37-reg4; reg63=reg63/reg9; reg67=reg67/reg9; reg106=reg106*reg44; reg37=reg29*reg69;
    reg56=reg41-reg56; reg82=reg82-reg92; reg12=reg12-reg86; reg41=reg90*reg13; reg67=reg67*reg24;
    reg83=reg66+reg83; reg89=reg22*reg89; reg22=reg7*reg29; reg107=reg106-reg107; reg44=reg42*reg44;
    reg27=reg16-reg27; reg9=reg4/reg9; reg63=reg63*reg24; reg4=reg95*reg29; reg65=reg57+reg65;
    reg16=reg91*reg13; reg42=reg91*reg23; reg57=reg91*reg7; reg4=reg16-reg4; reg16=reg28*reg29;
    reg22=reg41-reg22; reg41=reg95*reg90; reg62=reg96*reg29; reg66=reg90*reg23; reg37=reg56+reg37;
    reg89=reg44-reg89; reg24=reg9*reg24; reg63=reg27-reg63; reg9=reg82*reg83; reg27=reg12*reg65;
    reg107=reg67+reg107; reg44=1-(*f.m).resolution; reg56=(*f.m).resolution*reg115; reg9=reg27-reg9; reg16=reg66-reg16;
    reg22=reg22/reg37; reg107=reg44*reg107; reg62=reg42-reg62; reg4=reg4/reg37; reg27=reg91*reg28;
    reg41=reg57-reg41; reg42=reg96*reg90; reg89=reg24+reg89; reg24=(*f.m).resolution*reg60; reg63=reg44*reg63;
    reg63=reg56+reg63; reg82=reg82/reg9; reg83=reg83/reg9; reg12=reg12/reg9; reg65=reg65/reg9;
    reg56=(*f.m).resolution*reg103; reg57=(*f.m).resolution*reg22; reg66=(*f.m).resolution*reg4; reg89=reg44*reg89; reg107=reg24+reg107;
    reg62=reg62/reg37; reg70=reg70/reg37; reg42=reg27-reg42; reg35=reg35*reg44; reg16=reg16/reg37;
    reg41=reg41/reg37; reg111=reg111*reg44; reg80=reg80/reg37; reg58=reg58*reg44; reg24=(*f.m).resolution*reg80;
    reg27=(*f.m).resolution*reg41; reg69=reg69/reg37; reg37=reg42/reg37; reg42=(*f.m).resolution*reg70; reg71=reg71*reg44;
    reg67=(*f.m).resolution*reg62; reg66=reg35-reg66; reg35=(*f.m).resolution*reg16; reg111=reg57+reg111; reg63=(*f.m).deltaT*reg63;
    reg40=reg40*reg44; reg121=reg121*reg44; reg107=(*f.m).deltaT*reg107; reg89=reg56+reg89; reg130=reg130*reg44;
    reg56=reg25*reg82; reg57=reg11*reg83; reg94=reg25*reg12; reg97=var_inter[0]*reg82; reg100=var_inter[1]*reg83;
    reg101=var_inter[1]*reg65; reg102=reg11*reg65; reg104=var_inter[0]*reg12; reg105=reg66*reg63; reg71=reg27+reg71;
    reg27=reg111*reg107; reg106=(*f.m).resolution*reg37; reg40=reg67+reg40; reg35=reg130-reg35; reg67=reg56-reg102;
    reg108=reg102+reg97; reg109=reg104+reg57; reg110=reg57-reg94; reg112=(*f.m).resolution*reg69; reg42=reg121-reg42;
    reg113=reg100+reg94; reg58=reg24+reg58; reg89=(*f.m).deltaT*reg89; reg24=reg101+reg56; reg114=reg101-reg97;
    reg17=reg17*reg44; reg117=reg104-reg100; reg44=reg93*reg44; reg93=reg71*reg89; reg118=reg35*reg107;
    reg120=reg27+reg105; reg121=reg40*reg63; reg124=reg58*reg107; reg125=reg42*reg63; reg126=0.5*reg110;
    reg129=0.5*reg108; reg130=0.5*reg67; reg133=0.5*reg113; reg136=0.5*reg24; reg137=0.5*reg109;
    reg138=0.5*reg117; reg139=0.5*reg114; reg17=reg112+reg17; reg106=reg44-reg106; reg44=reg126*reg71;
    reg112=reg67*reg111; reg140=reg139*reg71; reg141=reg117*reg66; reg142=reg110*reg66; T reg143=reg130*reg71;
    T reg144=reg129*reg71; T reg145=reg109*reg66; T reg146=reg137*reg71; T reg147=reg138*reg71; T reg148=reg125+reg124;
    T reg149=reg17*reg89; T reg150=reg114*reg111; T reg151=reg108*reg111; T reg152=reg106*reg89; T reg153=reg118+reg121;
    T reg154=reg136*reg71; T reg155=reg120+reg93; T reg156=reg113*reg66; T reg157=reg133*reg71; T reg158=reg24*reg111;
    T reg159=reg110*reg42; T reg160=reg130*reg17; reg44=reg112+reg44; reg112=reg126*reg17; T reg161=reg67*reg58;
    reg157=reg157-reg158; T reg162=reg136*reg106; reg143=reg142+reg143; reg144=reg144-reg145; reg142=reg137*reg17;
    T reg163=reg113*reg42; T reg164=reg129*reg17; T reg165=reg109*reg42; reg151=reg151-reg146; T reg166=reg136*reg17;
    reg156=reg156-reg154; T reg167=reg130*reg106; T reg168=reg110*reg40; T reg169=reg67*reg35; T reg170=reg126*reg106;
    T reg171=reg108*reg58; T reg172=reg113*reg40; T reg173=2*reg155; T reg174=reg138*reg106; T reg175=reg117*reg40;
    T reg176=reg139*reg106; T reg177=reg117*reg42; T reg178=reg114*reg35; T reg179=reg153+reg152; T reg180=reg24*reg35;
    T reg181=reg139*reg17; T reg182=reg133*reg106; reg140=reg141+reg140; reg141=reg129*reg106; T reg183=reg109*reg40;
    reg147=reg150+reg147; reg150=reg133*reg17; T reg184=reg24*reg58; T reg185=reg108*reg35; T reg186=reg137*reg106;
    T reg187=reg138*reg17; T reg188=reg148+reg149; T reg189=reg114*reg58; reg182=reg182-reg180; T reg190=reg108*reg188;
    reg150=reg150-reg184; T reg191=reg109*reg179; reg176=reg175+reg176; reg171=reg171-reg142; reg187=reg189+reg187;
    reg157=2*reg157; reg147=2*reg147; reg175=reg137*reg173; reg172=reg172-reg162; reg163=reg163-reg166;
    reg151=2*reg151; reg189=reg129*reg173; reg44=2*reg44; T reg192=reg136*reg173; T reg193=reg113*reg179;
    reg185=reg185-reg186; reg160=reg159+reg160; reg112=reg161+reg112; reg159=reg24*reg188; reg140=2*reg140;
    reg167=reg168+reg167; reg161=reg133*reg173; reg181=reg177+reg181; reg170=reg169+reg170; reg141=reg141-reg183;
    reg144=2*reg144; reg174=reg178+reg174; reg156=2*reg156; reg143=2*reg143; reg164=reg164-reg165;
    reg168=reg139*reg151; reg169=reg129*reg44; reg177=reg109*reg170; reg178=reg137*reg157; T reg194=reg139*reg143;
    T reg195=reg117*reg167; T reg196=reg117*reg170; T reg197=reg117*reg185; T reg198=reg137*reg156; T reg199=reg108*reg163;
    T reg200=reg114*reg163; T reg201=reg108*reg150; T reg202=reg109*reg176; T reg203=reg113*reg185; T reg204=reg129*reg140;
    T reg205=reg109*reg174; T reg206=reg129*reg147; T reg207=reg109*reg182; T reg208=reg129*reg156; T reg209=reg109*reg172;
    T reg210=reg138*reg44; T reg211=reg114*reg112; T reg212=reg138*reg143; T reg213=reg114*reg160; T reg214=reg138*reg151;
    T reg215=reg114*reg171; T reg216=reg138*reg144; T reg217=reg109*reg141; T reg218=reg114*reg164; T reg219=reg138*reg147;
    T reg220=reg114*reg187; T reg221=reg138*reg140; T reg222=reg114*reg181; T reg223=reg138*reg157; T reg224=reg129*reg144;
    T reg225=reg109*reg185; T reg226=reg114*reg150; T reg227=reg138*reg156; T reg228=reg129*reg151; T reg229=reg109*reg167;
    T reg230=reg129*reg143; T reg231=reg136*reg44; T reg232=reg113*reg167; T reg233=reg136*reg143; T reg234=reg136*reg151;
    T reg235=reg113*reg141; T reg236=reg136*reg144; T reg237=reg113*reg174; T reg238=reg136*reg147; T reg239=reg113*reg176;
    T reg240=reg136*reg140; T reg241=reg113*reg182; T reg242=reg136*reg157; T reg243=reg113*reg172; T reg244=reg136*reg156;
    T reg245=reg67*reg188; T reg246=reg126*reg173; T reg247=reg130*reg173; T reg248=reg110*reg179; T reg249=reg175-reg190;
    T reg250=reg191-reg189; T reg251=reg138*reg173; T reg252=reg114*reg188; T reg253=reg117*reg179; T reg254=reg139*reg173;
    T reg255=reg159-reg161; T reg256=reg192-reg193; T reg257=reg117*reg141; T reg258=reg139*reg144; T reg259=reg117*reg174;
    T reg260=reg139*reg147; T reg261=reg117*reg176; T reg262=reg139*reg140; T reg263=reg117*reg182; T reg264=reg139*reg157;
    T reg265=reg117*reg172; T reg266=reg139*reg156; T reg267=reg133*reg44; T reg268=reg24*reg112; T reg269=reg133*reg143;
    T reg270=reg24*reg160; T reg271=reg133*reg151; T reg272=reg24*reg171; T reg273=reg133*reg144; T reg274=reg24*reg164;
    T reg275=reg133*reg147; T reg276=reg133*reg140; T reg277=reg24*reg181; T reg278=reg133*reg157; T reg279=reg24*reg150;
    T reg280=reg133*reg156; T reg281=reg24*reg163; T reg282=reg113*reg170; T reg283=reg67*reg164; T reg284=reg126*reg156;
    T reg285=reg67*reg187; T reg286=reg126*reg144; T reg287=reg137*reg144; T reg288=reg130*reg143; reg164=reg108*reg164;
    T reg289=reg126*reg140; reg176=reg110*reg176; reg156=reg130*reg156; T reg290=reg67*reg181; T reg291=reg139*reg44;
    T reg292=reg126*reg143; reg170=reg110*reg170; T reg293=reg67*reg112; T reg294=reg137*reg151; T reg295=reg108*reg171;
    T reg296=reg126*reg44; T reg297=reg108*reg160; reg167=reg110*reg167; reg174=reg110*reg174; T reg298=reg130*reg147;
    T reg299=reg130*reg140; T reg300=reg130*reg44; reg160=reg67*reg160; reg143=reg137*reg143; T reg301=reg130*reg157;
    reg185=reg110*reg185; reg144=reg130*reg144; T reg302=reg126*reg147; reg140=reg137*reg140; T reg303=reg126*reg157;
    reg150=reg67*reg150; reg157=reg129*reg157; reg181=reg108*reg181; T reg304=reg24*reg187; reg171=reg67*reg171;
    T reg305=reg130*reg151; reg141=reg110*reg141; reg151=reg126*reg151; reg172=reg110*reg172; reg182=reg110*reg182;
    reg112=reg108*reg112; reg147=reg137*reg147; reg187=reg108*reg187; reg44=reg137*reg44; reg163=reg67*reg163;
    T reg306=reg245+reg246; reg172=reg156+reg172; reg231=reg282-reg231; reg200=reg227+reg200; reg170=reg300+reg170;
    reg292=reg160+reg292; reg226=reg223+reg226; reg194=reg195+reg194; reg250=reg9*reg250; reg303=reg150+reg303;
    reg274=reg273-reg274; reg272=reg271-reg272; reg150=reg247+reg248; reg270=reg269-reg270; reg304=reg275-reg304;
    reg268=reg267-reg268; reg151=reg171+reg151; reg176=reg299+reg176; reg266=reg265+reg266; reg182=reg301+reg182;
    reg264=reg263+reg264; reg286=reg283+reg286; reg277=reg276-reg277; reg262=reg261+reg262; reg249=reg9*reg249;
    reg284=reg163+reg284; reg260=reg259+reg260; reg279=reg278-reg279; reg258=reg257+reg258; reg296=reg293+reg296;
    reg168=reg197+reg168; reg281=reg280-reg281; reg209=reg208-reg209; reg238=reg237-reg238; reg207=reg157-reg207;
    reg294=reg295-reg294; reg202=reg204-reg202; reg289=reg290+reg289; reg240=reg239-reg240; reg205=reg206-reg205;
    reg287=reg164-reg287; reg255=reg9*reg255; reg217=reg224-reg217; reg167=reg288+reg167; reg242=reg241-reg242;
    reg225=reg228-reg225; reg141=reg144+reg141; reg229=reg230-reg229; reg147=reg187-reg147; reg244=reg243-reg244;
    reg177=reg169-reg177; reg256=reg9*reg256; reg198=reg199-reg198; reg140=reg181-reg140; reg178=reg201-reg178;
    reg185=reg305+reg185; reg144=reg251+reg252; reg218=reg216+reg218; reg44=reg112-reg44; reg220=reg219+reg220;
    reg234=reg203-reg234; reg215=reg214+reg215; reg174=reg298+reg174; reg213=reg212+reg213; reg143=reg297-reg143;
    reg236=reg235-reg236; reg302=reg285+reg302; reg211=reg210+reg211; reg222=reg221+reg222; reg112=reg253+reg254;
    reg233=reg232-reg233; reg196=reg291+reg196; reg244=reg9*reg244; reg156=reg9*reg306; reg274=reg9*reg274;
    reg231=reg9*reg231; reg233=reg9*reg233; reg170=reg9*reg170; reg185=reg9*reg185; reg176=reg9*reg176;
    reg174=reg9*reg174; reg242=reg9*reg242; reg304=reg9*reg304; reg167=reg9*reg167; reg141=reg9*reg141;
    reg234=reg9*reg234; reg240=reg9*reg240; reg281=reg9*reg281; reg277=reg9*reg277; reg236=reg9*reg236;
    reg279=reg9*reg279; reg238=reg9*reg238; reg196=reg9*reg196; reg284=reg9*reg284; reg222=reg9*reg222;
    reg220=reg9*reg220; reg218=reg9*reg218; reg44=reg9*reg44; reg157=reg9*reg144; reg215=reg9*reg215;
    reg302=reg9*reg302; reg213=reg9*reg213; reg211=reg9*reg211; reg143=reg9*reg143; reg160=reg9*reg112;
    reg209=reg9*reg209; reg207=reg9*reg207; reg294=reg9*reg294; reg202=reg9*reg202; reg205=reg9*reg205;
    reg217=reg9*reg217; reg287=reg9*reg287; reg255=ponderation*reg255; reg225=reg9*reg225; reg289=reg9*reg289;
    reg229=reg9*reg229; reg147=reg9*reg147; reg177=reg9*reg177; reg198=reg9*reg198; reg256=ponderation*reg256;
    reg178=reg9*reg178; reg140=reg9*reg140; reg266=reg9*reg266; reg264=reg9*reg264; reg182=reg9*reg182;
    reg262=reg9*reg262; reg286=reg9*reg286; reg260=reg9*reg260; reg249=ponderation*reg249; reg258=reg9*reg258;
    reg151=reg9*reg151; reg296=reg9*reg296; reg268=reg9*reg268; reg168=reg9*reg168; reg163=reg9*reg150;
    reg270=reg9*reg270; reg272=reg9*reg272; reg303=reg9*reg303; reg194=reg9*reg194; reg292=reg9*reg292;
    reg200=reg9*reg200; reg250=ponderation*reg250; reg226=reg9*reg226; reg172=reg9*reg172; matrix(indices[0]+0,indices[2]+1)+=ponderation*reg289;
    matrix(indices[0]+1,indices[2]+1)+=ponderation*reg176; reg164=ponderation*reg163; sollicitation[indices[0]+1]+=reg164; sollicitation[indices[3]+1]+=-reg256; sollicitation[indices[3]+0]+=-reg255;
    matrix(indices[0]+1,indices[1]+1)+=ponderation*reg141; reg141=ponderation*reg160; sollicitation[indices[2]+1]+=reg141; sollicitation[indices[1]+0]+=-reg249; matrix(indices[0]+1,indices[2]+0)+=ponderation*reg174;
    matrix(indices[0]+0,indices[2]+0)+=ponderation*reg302; matrix(indices[0]+1,indices[3]+0)+=ponderation*reg182; matrix(indices[0]+0,indices[1]+1)+=ponderation*reg286; reg169=ponderation*reg157; sollicitation[indices[2]+0]+=reg169;
    reg171=ponderation*reg156; sollicitation[indices[0]+0]+=reg171; sollicitation[indices[1]+1]+=-reg250; matrix(indices[0]+0,indices[0]+0)+=ponderation*reg296; matrix(indices[2]+1,indices[0]+1)+=ponderation*reg194;
    matrix(indices[2]+0,indices[3]+1)+=ponderation*reg200; matrix(indices[0]+0,indices[0]+1)+=ponderation*reg292; matrix(indices[2]+0,indices[3]+0)+=ponderation*reg226; matrix(indices[2]+0,indices[2]+1)+=ponderation*reg222; matrix(indices[0]+1,indices[3]+1)+=ponderation*reg172;
    matrix(indices[2]+0,indices[2]+0)+=ponderation*reg220; matrix(indices[2]+0,indices[1]+1)+=ponderation*reg218; matrix(indices[2]+0,indices[1]+0)+=ponderation*reg215; matrix(indices[1]+0,indices[0]+0)+=ponderation*reg44; matrix(indices[2]+0,indices[0]+1)+=ponderation*reg213;
    matrix(indices[2]+0,indices[0]+0)+=ponderation*reg211; matrix(indices[1]+0,indices[0]+1)+=ponderation*reg143; matrix(indices[1]+1,indices[3]+1)+=ponderation*reg209; matrix(indices[1]+1,indices[3]+0)+=ponderation*reg207; matrix(indices[1]+1,indices[2]+1)+=ponderation*reg202;
    matrix(indices[1]+0,indices[1]+0)+=ponderation*reg294; matrix(indices[1]+1,indices[2]+0)+=ponderation*reg205; matrix(indices[1]+1,indices[1]+1)+=ponderation*reg217; matrix(indices[1]+0,indices[1]+1)+=ponderation*reg287; matrix(indices[1]+1,indices[1]+0)+=ponderation*reg225;
    matrix(indices[1]+1,indices[0]+1)+=ponderation*reg229; matrix(indices[1]+1,indices[0]+0)+=ponderation*reg177; matrix(indices[1]+0,indices[2]+0)+=ponderation*reg147; matrix(indices[1]+0,indices[3]+1)+=ponderation*reg198; matrix(indices[1]+0,indices[3]+0)+=ponderation*reg178;
    matrix(indices[1]+0,indices[2]+1)+=ponderation*reg140; matrix(indices[0]+1,indices[1]+0)+=ponderation*reg185; matrix(indices[3]+1,indices[3]+1)+=ponderation*reg244; matrix(indices[0]+1,indices[0]+1)+=ponderation*reg167; matrix(indices[3]+1,indices[3]+0)+=ponderation*reg242;
    matrix(indices[3]+1,indices[2]+1)+=ponderation*reg240; matrix(indices[3]+1,indices[2]+0)+=ponderation*reg238; matrix(indices[3]+1,indices[1]+1)+=ponderation*reg236; matrix(indices[2]+1,indices[0]+0)+=ponderation*reg196; matrix(indices[3]+1,indices[1]+0)+=ponderation*reg234;
    matrix(indices[3]+1,indices[0]+1)+=ponderation*reg233; matrix(indices[0]+1,indices[0]+0)+=ponderation*reg170; matrix(indices[3]+1,indices[0]+0)+=ponderation*reg231; matrix(indices[3]+0,indices[3]+1)+=ponderation*reg281; matrix(indices[0]+0,indices[3]+1)+=ponderation*reg284;
    matrix(indices[3]+0,indices[3]+0)+=ponderation*reg279; matrix(indices[3]+0,indices[2]+1)+=ponderation*reg277; matrix(indices[3]+0,indices[2]+0)+=ponderation*reg304; matrix(indices[0]+0,indices[3]+0)+=ponderation*reg303; matrix(indices[3]+0,indices[1]+1)+=ponderation*reg274;
    matrix(indices[3]+0,indices[1]+0)+=ponderation*reg272; matrix(indices[3]+0,indices[0]+1)+=ponderation*reg270; matrix(indices[3]+0,indices[0]+0)+=ponderation*reg268; matrix(indices[2]+1,indices[3]+1)+=ponderation*reg266; matrix(indices[0]+0,indices[1]+0)+=ponderation*reg151;
    matrix(indices[2]+1,indices[3]+0)+=ponderation*reg264; matrix(indices[2]+1,indices[2]+1)+=ponderation*reg262; matrix(indices[2]+1,indices[2]+0)+=ponderation*reg260; matrix(indices[2]+1,indices[1]+1)+=ponderation*reg258; matrix(indices[2]+1,indices[1]+0)+=ponderation*reg168;
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
    T reg3=reg1*reg0; reg2=1.0/reg2; T reg4=reg2*reg3; T reg5=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg6=1.0/(*f.m).elastic_modulus_3;
    T reg7=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1; T reg8=pow((*f.m).v1[0],2); T reg9=pow((*f.m).v1[1],2); T reg10=pow((*f.m).v2[1],2); T reg11=pow((*f.m).v2[0],2);
    T reg12=reg7*reg4; reg10=reg11+reg10; reg11=reg5*reg4; T reg13=reg6*reg4; T reg14=pow((*f.m).v2[2],2);
    reg9=reg8+reg9; reg8=pow((*f.m).v1[2],2); T reg15=1.0/(*f.m).elastic_modulus_2; T reg16=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; reg14=reg10+reg14;
    reg10=reg15*reg13; T reg17=reg5*reg11; T reg18=reg16*reg13; reg8=reg9+reg8; reg9=reg5*reg12;
    reg14=pow(reg14,0.5); T reg19=1.0/(*f.m).elastic_modulus_1; reg17=reg10-reg17; reg8=pow(reg8,0.5); reg9=reg18+reg9;
    reg10=reg16*reg11; T reg20=reg15*reg12; T reg21=reg16*reg9; T reg22=reg19*reg17; T reg23=reg10+reg20;
    T reg24=(*f.m).v2[1]/reg14; T reg25=(*f.m).v2[2]/reg14; T reg26=(*f.m).v1[2]/reg8; T reg27=(*f.m).v1[1]/reg8; reg14=(*f.m).v2[0]/reg14;
    T reg28=reg7*reg3; T reg29=reg16*reg4; T reg30=reg26*reg24; T reg31=reg27*reg25; T reg32=reg6*reg3;
    T reg33=reg7*reg12; reg3=reg5*reg3; reg13=reg19*reg13; reg8=(*f.m).v1[0]/reg8; T reg34=reg7*reg11;
    reg4=reg15*reg4; T reg35=reg2*reg0; T reg36=reg7*reg23; reg21=reg22-reg21; reg22=reg8*reg25;
    T reg37=reg26*reg14; T reg38=reg31-reg30; T reg39=2*reg8; T reg40=reg6*reg35; T reg41=2*reg14;
    T reg42=reg7*reg29; reg11=reg19*reg11; reg33=reg13-reg33; reg13=reg15*reg32; reg32=reg16*reg32;
    reg34=reg18+reg34; reg18=reg7*reg4; T reg43=reg5*reg35; reg12=reg16*reg12; reg36=reg21-reg36;
    reg35=reg7*reg35; reg21=reg5*reg3; T reg44=reg5*reg28; T reg45=reg2*reg1; T reg46=reg37-reg22;
    T reg47=reg7*reg45; reg28=reg15*reg28; T reg48=pow(reg8,2); T reg49=reg5*reg45; T reg50=reg8*reg24;
    T reg51=pow(reg27,2); T reg52=reg27*reg39; T reg53=reg15*reg40; reg40=reg16*reg40; T reg54=pow(reg14,2);
    T reg55=pow(reg24,2); T reg56=reg5*reg43; T reg57=reg24*reg41; T reg58=reg5*reg35; reg45=reg6*reg45;
    T reg59=2*reg38; reg44=reg32+reg44; reg21=reg13-reg21; reg17=reg17/reg36; reg34=reg34/reg36;
    reg9=reg9/reg36; reg18=reg10+reg18; reg33=reg33/reg36; reg4=reg19*reg4; reg42=reg11+reg42;
    reg13=reg27*reg14; reg12=reg11+reg12; reg29=reg16*reg29; reg3=reg16*reg3; reg11=reg55*reg9;
    reg32=reg51*reg17; T reg60=reg5*reg47; T reg61=reg5*reg49; T reg62=reg16*reg45; reg45=reg15*reg45;
    reg44=reg16*reg44; reg58=reg40+reg58; reg56=reg53-reg56; reg21=reg19*reg21; reg23=reg23/reg36;
    reg40=reg54*reg9; reg29=reg4-reg29; reg18=reg18/reg36; reg12=reg12/reg36; reg4=reg48*reg17;
    reg43=reg16*reg43; reg53=reg3+reg28; reg35=reg15*reg35; reg42=reg42/reg36; T reg63=pow(reg46,2);
    T reg64=pow(reg38,2); T reg65=reg57*reg9; T reg66=reg52*reg17; T reg67=reg46*reg59; T reg68=reg50-reg13;
    T reg69=reg57*reg33; T reg70=reg52*reg34; T reg71=reg48*reg34; T reg72=reg54*reg33; T reg73=pow(reg26,2);
    T reg74=reg51*reg34; T reg75=reg55*reg33; T reg76=pow(reg25,2); reg53=reg7*reg53; T reg77=reg63*reg23;
    T reg78=reg73*reg17; T reg79=reg76*reg9; T reg80=reg51*reg18; T reg81=reg54*reg42; T reg82=reg55*reg42;
    reg56=reg19*reg56; reg29=reg29/reg36; reg65=reg66+reg65; reg66=reg67*reg23; reg44=reg21-reg44;
    reg21=reg48*reg18; T reg83=reg67*reg12; reg69=reg70+reg69; reg72=reg71+reg72; reg70=reg64*reg12;
    reg75=reg74+reg75; reg71=reg63*reg12; reg74=reg73*reg34; T reg84=reg76*reg33; reg11=reg32+reg11;
    reg47=reg15*reg47; reg49=reg16*reg49; reg32=reg57*reg42; T reg85=reg64*reg23; reg4=reg40+reg4;
    reg60=reg62+reg60; reg40=reg52*reg18; reg61=reg45-reg61; reg45=pow(reg68,2); reg62=reg43+reg35;
    reg58=reg16*reg58; T reg86=reg67*reg29; reg70=reg72+reg70; reg32=reg40+reg32; reg85=reg4+reg85;
    reg82=reg80+reg82; reg4=reg49+reg47; reg40=reg8*reg27; reg72=reg14*reg24; reg80=reg25*reg41;
    T reg87=reg63*reg29; reg83=reg69+reg83; reg69=reg45*reg12; reg60=reg16*reg60; reg66=reg65+reg66;
    reg84=reg74+reg84; reg71=reg75+reg71; reg58=reg56-reg58; reg56=reg73*reg18; reg65=reg64*reg29;
    reg74=2*reg24; reg77=reg11+reg77; reg62=reg7*reg62; reg53=reg44-reg53; reg79=reg78+reg79;
    reg11=reg45*reg23; reg81=reg21+reg81; reg21=reg76*reg42; reg44=reg26*reg39; reg61=reg19*reg61;
    reg75=2*reg27; reg69=reg84+reg69; reg78=reg51*reg85; reg84=reg54*reg71; reg86=reg32+reg86;
    reg32=reg48*reg77; T reg88=reg55*reg70; T reg89=reg44*reg34; T reg90=reg80*reg33; reg87=reg82+reg87;
    reg82=reg45*reg29; reg21=reg56+reg21; reg65=reg81+reg65; reg56=reg55*reg83; reg81=reg46*reg38;
    T reg91=reg26*reg75; T reg92=2*reg26; T reg93=reg24*reg75; T reg94=reg14*reg39; T reg95=reg40*reg77;
    reg62=reg58-reg62; reg58=reg25*reg74; reg60=reg61-reg60; reg61=reg51*reg66; reg4=reg7*reg4;
    T reg96=reg8*reg46; T reg97=reg27*reg38; T reg98=reg55*reg71; reg59=reg68*reg59; reg50=reg13+reg50;
    reg13=2*reg46; T reg99=reg27*reg24; T reg100=reg8*reg14; T reg101=reg72*reg83; T reg102=reg40*reg66;
    T reg103=reg80*reg9; reg77=reg51*reg77; T reg104=reg48*reg85; T reg105=reg54*reg70; T reg106=reg44*reg17;
    reg83=reg54*reg83; reg66=reg48*reg66; reg11=reg79+reg11; reg71=reg72*reg71; reg53=reg53/reg36;
    reg79=reg44*reg18; reg9=reg58*reg9; reg17=reg91*reg17; T reg107=reg59*reg23; reg4=reg60-reg4;
    reg103=reg106+reg103; reg82=reg21+reg82; reg90=reg89+reg90; reg21=reg59*reg12; reg34=reg91*reg34;
    reg33=reg58*reg33; reg60=reg80*reg42; reg62=reg62/reg36; reg89=reg50*reg53; reg98=reg77+reg98;
    reg77=reg48*reg11; reg96=reg97+reg96; reg97=reg27*reg46; reg106=reg8*reg38; T reg108=reg26*reg25;
    reg13=reg68*reg13; reg71=reg95+reg71; reg95=reg81*reg87; reg83=reg66+reg83; reg66=reg64*reg86;
    reg101=reg102+reg101; reg102=reg81*reg86; reg105=reg104+reg105; reg104=reg64*reg65; T reg109=reg54*reg19;
    T reg110=reg55*reg16; T reg111=reg94*reg19; T reg112=reg93*reg16; T reg113=reg55*reg69; T reg114=reg51*reg11;
    T reg115=reg54*reg16; T reg116=reg55*reg15; T reg117=reg94*reg16; T reg118=reg93*reg15; T reg119=reg63*reg65;
    reg88=reg78+reg88; reg84=reg32+reg84; reg32=reg64*reg87; reg56=reg61+reg56; reg86=reg63*reg86;
    reg61=reg24*reg38; reg78=reg100*reg53; T reg120=reg14*reg46; T reg121=reg54*reg69; reg87=reg63*reg87;
    reg70=reg72*reg70; reg92=reg25*reg92; T reg122=reg99*reg53; reg85=reg40*reg85; T reg123=reg64*reg82;
    T reg124=reg59*reg29; reg110=reg109-reg110; reg109=reg76*reg7; T reg125=reg55*reg5; T reg126=reg94*reg78;
    reg104=reg105+reg104; reg75=reg46*reg75; reg113=reg114+reg113; reg105=reg94*reg89; reg66=reg83+reg66;
    reg83=reg93*reg89; reg114=reg92*reg5; reg117=reg118-reg117; reg119=reg88+reg119; reg115=reg116-reg115;
    reg18=reg91*reg18; reg88=reg76*reg5; reg12=reg13*reg12; reg33=reg34+reg33; reg34=reg93*reg78;
    reg121=reg77+reg121; reg21=reg90+reg21; reg42=reg58*reg42; reg86=reg56+reg86; reg112=reg111-reg112;
    reg56=reg92*reg7; reg77=reg54*reg7; reg65=reg81*reg65; reg36=reg4/reg36; reg4=reg14*reg38;
    reg90=reg24*reg46; reg111=reg108*reg53; reg116=reg106*reg62; reg118=reg26*reg68; reg120=reg61+reg120;
    reg61=reg97*reg62; reg70=reg85+reg70; reg95=reg71+reg95; reg71=reg50*reg122; reg85=reg96*reg62;
    reg11=reg40*reg11; T reg127=reg63*reg82; T reg128=reg94*reg7; T reg129=reg93*reg5; reg89=reg50*reg89;
    reg60=reg79+reg60; reg23=reg13*reg23; reg9=reg17+reg9; reg32=reg84+reg32; reg102=reg101+reg102;
    reg107=reg103+reg107; reg39=reg38*reg39; reg17=reg93*reg122; reg87=reg98+reg87; reg122=reg94*reg122;
    reg69=reg72*reg69; reg79=reg39*reg61; reg122=reg32+reg122; reg124=reg60+reg124; reg32=reg75*reg61;
    reg60=reg39*reg116; reg42=reg18+reg42; reg126=reg104+reg126; reg29=reg13*reg29; reg19=reg48*reg19;
    reg18=reg120*reg36; reg127=reg113+reg127; reg84=reg90*reg36; reg98=reg4*reg36; reg101=reg93*reg111;
    reg17=reg87+reg17; reg87=reg118*reg62; reg89=reg102+reg89; reg102=reg96*reg85; reg103=(*f.m).alpha_1*reg51;
    reg104=(*f.m).alpha_2*reg55; reg105=reg66+reg105; reg66=reg39*reg85; reg113=reg48*reg107; T reg130=reg51*reg16;
    T reg131=reg54*reg21; T reg132=reg75*reg116; reg109=reg110-reg109; reg56=reg112-reg56; reg34=reg119+reg34;
    reg15=reg51*reg15; reg88=reg115-reg88; reg114=reg117-reg114; reg110=reg76*reg6; reg125=reg77+reg125;
    reg92=reg92*reg6; reg129=reg128+reg129; reg77=(*f.m).alpha_1*reg48; reg12=reg33+reg12; reg23=reg9+reg23;
    reg16=reg48*reg16; reg83=reg86+reg83; reg85=reg75*reg85; reg9=reg51*reg107; reg33=reg55*reg21;
    reg123=reg121+reg123; reg86=reg94*reg111; reg74=reg46*reg74; reg41=reg38*reg41; reg112=reg25*reg68;
    reg65=reg70+reg65; reg78=reg50*reg78; reg70=reg8*reg68; reg115=reg26*reg38; reg37=reg22+reg37;
    reg71=reg95+reg71; reg61=reg96*reg61; reg69=reg11+reg69; reg82=reg81*reg82; reg11=(*f.m).alpha_2*reg54;
    reg132=reg34+reg132; reg22=reg74*reg98; reg34=reg54*reg12; reg102=reg89+reg102; reg89=reg120*reg18;
    reg107=reg40*reg107; reg21=reg72*reg21; reg95=reg55*reg88; reg117=reg54*reg109; reg119=reg51*reg114;
    reg121=reg48*reg56; reg128=reg73*reg7; reg130=reg19-reg130; reg19=reg73*reg5; reg16=reg15-reg16;
    reg7=reg48*reg7; reg5=reg51*reg5; reg125=reg110-reg125; reg129=reg92-reg129; reg15=reg51*reg88;
    reg92=reg48*reg109; reg32=reg17+reg32; reg17=reg74*reg84; reg101=reg127+reg101; reg110=reg75*reg87;
    reg85=reg83+reg85; reg83=reg74*reg18; reg33=reg9+reg33; reg9=reg63*reg124; reg127=reg55*reg114;
    T reg133=reg51*reg23; T reg134=reg55*reg12; T reg135=reg54*reg56; reg78=reg65+reg78; reg116=reg96*reg116;
    reg61=reg71+reg61; reg65=reg120*reg84; reg82=reg69+reg82; reg111=reg50*reg111; reg69=reg39*reg87;
    reg86=reg123+reg86; reg84=reg41*reg84; reg79=reg122+reg79; reg71=reg41*reg98; reg60=reg126+reg60;
    reg122=reg112*reg36; reg123=reg37*reg53; reg29=reg42+reg29; reg42=reg14*reg68; reg126=reg25*reg38;
    T reg136=reg27*reg68; T reg137=reg26*reg46; reg70=reg115+reg70; reg31=reg30+reg31; reg30=reg100*reg109;
    reg115=reg99*reg88; T reg138=reg100*reg56; T reg139=reg99*reg114; reg11=reg77+reg11; reg77=(*f.m).alpha_3*reg64;
    reg104=reg103+reg104; reg103=(*f.m).alpha_3*reg63; T reg140=(*f.m).alpha_1*reg73; T reg141=(*f.m).alpha_2*reg76; reg66=reg105+reg66;
    reg105=reg64*reg124; reg18=reg41*reg18; T reg142=reg48*reg23; reg131=reg113+reg131; reg87=reg96*reg87;
    reg111=reg82+reg111; reg82=reg108*reg129; reg139=reg138+reg139; reg65=reg61+reg65; reg6=reg73*reg6;
    reg5=reg7+reg5; reg7=reg108*reg125; reg115=reg30+reg115; reg30=reg76*reg129; reg98=reg120*reg98;
    reg116=reg78+reg116; reg127=reg135+reg127; reg136=reg137+reg136; reg61=reg94*reg123; reg42=reg126+reg42;
    reg78=reg25*reg46; reg113=reg63*reg29; reg134=reg133+reg134; reg126=(*f.m).alpha_1*reg40; reg133=(*f.m).alpha_2*reg72;
    reg135=(*f.m).alpha_3*reg45; reg119=reg121+reg119; reg121=reg73*reg129; reg128=reg130-reg128; reg141=reg140+reg141;
    reg95=reg117+reg95; reg117=reg76*reg125; reg12=reg72*reg12; reg23=reg40*reg23; reg19=reg16-reg19;
    reg124=reg81*reg124; reg21=reg107+reg21; reg103=reg104+reg103; reg89=reg102+reg89; reg34=reg142+reg34;
    reg16=reg64*reg29; reg77=reg11+reg77; reg15=reg92+reg15; reg105=reg131+reg105; reg11=reg74*reg122;
    reg110=reg101+reg110; reg53=reg31*reg53; reg92=reg73*reg125; reg101=reg50*reg2; reg102=reg70*reg62;
    reg72=reg72*reg2; reg83=reg85+reg83; reg18=reg66+reg18; reg17=reg32+reg17; reg14=reg14*reg25;
    reg9=reg33+reg9; reg32=reg93*reg123; reg33=reg24*reg68; reg71=reg60+reg71; reg8=reg8*reg26;
    reg69=reg86+reg69; reg22=reg132+reg22; reg84=reg79+reg84; reg60=reg41*reg122; reg66=(*f.m).alpha_3*reg81;
    reg133=reg126+reg133; reg79=reg57*reg72; reg85=reg83*reg65; reg86=reg18*reg65; reg123=reg50*reg123;
    reg124=reg21+reg124; reg21=reg52*reg72; reg104=reg57*reg101; reg11=reg110+reg11; reg92=reg15+reg92;
    reg62=reg136*reg62; reg12=reg23+reg12; reg135=reg141+reg135; reg29=reg81*reg29; reg60=reg69+reg60;
    reg15=reg84*reg89; reg23=reg17*reg89; reg117=reg95+reg117; reg69=reg42*reg36; reg121=reg119+reg121;
    reg81=reg55*reg19; reg95=reg54*reg128; reg107=reg52*reg101; reg110=(*f.m).alpha_2*reg14; reg119=(*f.m).alpha_1*reg8;
    reg25=reg24*reg25; reg24=reg50*reg72; reg7=reg115+reg7; reg26=reg27*reg26; reg38=reg68*reg38;
    reg30=reg127+reg30; reg98=reg116+reg98; reg32=reg9+reg32; reg9=reg51*reg19; reg27=reg48*reg128;
    reg115=reg39*reg102; reg116=reg75*reg102; reg61=reg105+reg61; reg105=reg22*reg77; reg126=reg17*reg103;
    reg93=reg93*reg53; reg113=reg134+reg113; reg33=reg78+reg33; reg16=reg34+reg16; reg94=reg94*reg53;
    reg34=reg71*reg77; reg78=reg50*reg101; reg127=reg84*reg103; reg122=reg120*reg122; reg87=reg111+reg87;
    reg5=reg6-reg5; reg82=reg139+reg82; reg14=reg14*reg1; reg6=reg37*reg1; reg127=reg34+reg127;
    reg34=reg25*reg0; reg11=reg11*reg135; reg126=reg105+reg126; reg105=reg31*reg0; reg86=reg15-reg86;
    reg15=reg98*reg77; reg2=reg40*reg2; reg40=reg84*reg83; reg111=reg18*reg17; reg9=reg27+reg9;
    reg73=reg73*reg5; reg66=reg133+reg66; reg85=reg23-reg85; reg25=(*f.m).alpha_2*reg25; reg23=reg65*reg103;
    reg60=reg60*reg135; reg110=reg119+reg110; reg38=(*f.m).alpha_3*reg38; reg27=(*f.m).alpha_1*reg26; reg75=reg75*reg62;
    reg93=reg113+reg93; reg113=reg74*reg69; reg24=reg7+reg24; reg7=reg37*reg14; reg116=reg32+reg116;
    reg46=reg68*reg46; reg122=reg87+reg122; reg78=reg82+reg78; reg32=reg37*reg6; reg68=reg99*reg19;
    reg82=reg100*reg128; reg123=reg124+reg123; reg102=reg96*reg102; reg87=reg80*reg14; reg79=reg117+reg79;
    reg21=reg92+reg21; reg92=reg44*reg14; reg107=reg121+reg107; reg117=reg44*reg6; reg94=reg16+reg94;
    reg39=reg39*reg62; reg81=reg95+reg81; reg76=reg76*reg5; reg16=reg41*reg69; reg115=reg61+reg115;
    reg36=reg33*reg36; reg104=reg30+reg104; reg30=reg80*reg6; reg53=reg50*reg53; reg29=reg12+reg29;
    reg68=reg82+reg68; reg12=reg58*reg105; reg30=reg104+reg30; reg23=reg15+reg23; reg122=reg122*reg135;
    reg15=reg83*reg66; reg11=reg126+reg11; reg92=reg21+reg92; reg21=reg91*reg34; reg117=reg107+reg117;
    reg61=reg91*reg105; reg38=reg110+reg38; reg76=reg81+reg76; reg81=reg57*reg2; reg25=reg27+reg25;
    reg46=(*f.m).alpha_3*reg46; reg27=reg71*reg85; reg82=reg22*reg86; reg95=reg31*reg105; reg32=reg78+reg32;
    reg111=reg40-reg111; reg40=reg31*reg34; reg7=reg24+reg7; reg108=reg108*reg5; reg60=reg127+reg60;
    reg24=reg18*reg66; reg87=reg79+reg87; reg78=reg58*reg34; reg75=reg93+reg75; reg113=reg116+reg113;
    reg41=reg41*reg36; reg39=reg94+reg39; reg16=reg115+reg16; reg74=reg74*reg36; reg102=reg123+reg102;
    reg69=reg120*reg69; reg1=reg8*reg1; reg53=reg29+reg53; reg73=reg9+reg73; reg8=reg52*reg2;
    reg62=reg96*reg62; reg69=reg102+reg69; reg74=reg75+reg74; reg62=reg53+reg62; reg78=reg87+reg78;
    reg36=reg120*reg36; reg9=reg18*reg98; reg41=reg39+reg41; reg29=reg71*reg89; reg12=reg30+reg12;
    reg30=reg83*reg98; reg39=reg22*reg89; reg53=reg98*reg111; reg82=reg27-reg82; reg46=reg25+reg46;
    reg95=reg32+reg95; reg40=reg7+reg40; reg7=reg50*reg2; reg108=reg68+reg108; reg0=reg26*reg0;
    reg122=reg23+reg122; reg23=1-var_inter[1]; reg25=reg89*reg66; reg21=reg92+reg21; reg113=reg113*reg38;
    reg15=reg11+reg15; reg11=1-var_inter[0]; reg16=reg16*reg38; reg80=reg80*reg1; reg8=reg73+reg8;
    reg81=reg76+reg81; reg44=reg44*reg1; reg24=reg60+reg24; reg61=reg117+reg61; reg26=reg11*elem.pos(0)[1];
    reg7=reg108+reg7; reg25=reg122+reg25; reg27=elem.pos(1)[1]*var_inter[0]; reg32=reg11*elem.pos(0)[0]; reg69=reg69*reg38;
    reg60=elem.pos(1)[0]*var_inter[0]; reg37=reg37*reg1; reg68=reg23*elem.pos(0)[1]; reg73=reg23*elem.pos(1)[1]; reg75=reg21*reg95;
    reg53=reg82+reg53; reg76=reg78*reg95; reg79=reg23*elem.pos(0)[0]; reg74=reg74*reg46; reg82=reg22*reg65;
    reg113=reg15+reg113; reg30=reg39-reg30; reg15=reg17*reg98; reg39=reg71*reg65; reg41=reg41*reg46;
    reg16=reg24+reg16; reg24=reg23*elem.pos(1)[0]; reg87=reg18*reg22; reg9=reg29-reg9; reg29=reg84*reg98;
    reg92=reg71*reg83; reg36=reg62+reg36; reg62=reg61*reg40; reg44=reg8+reg44; reg91=reg91*reg0;
    reg58=reg58*reg0; reg8=reg12*reg40; reg80=reg81+reg80; reg58=reg80+reg58; reg29=reg39-reg29;
    reg39=reg71*reg17; reg87=reg92-reg87; reg80=reg84*reg22; reg41=reg16+reg41; reg74=reg113+reg74;
    reg91=reg44+reg91; reg69=reg25+reg69; reg36=reg36*reg46; reg73=reg73-reg68; reg16=var_inter[1]*elem.pos(2)[1];
    reg25=reg61*reg78; reg44=reg21*reg12; reg62=reg75-reg62; reg75=reg60+reg32; reg81=var_inter[1]*elem.pos(2)[0];
    reg92=reg27+reg26; reg93=elem.pos(2)[1]*var_inter[0]; reg24=reg24-reg79; reg94=elem.pos(2)[0]*var_inter[0]; reg85=reg85/reg53;
    reg30=reg30/reg53; reg8=reg76-reg8; reg15=reg82-reg15; reg86=reg86/reg53; reg9=reg9/reg53;
    reg37=reg7+reg37; reg31=reg31*reg0; reg36=reg69+reg36; reg30=reg30*reg41; reg80=reg39-reg80;
    reg87=reg87/reg53; reg111=reg111/reg53; reg29=reg29/reg53; reg9=reg9*reg74; reg31=reg37+reg31;
    reg15=reg15/reg53; reg7=reg58*reg62; reg93=reg93-reg92; reg81=reg24+reg81; reg24=elem.pos(3)[1]*reg11;
    reg37=elem.pos(3)[0]*reg11; reg94=reg94-reg75; reg39=reg8*reg91; reg69=var_inter[1]*elem.pos(3)[0]; reg16=reg73+reg16;
    reg73=var_inter[1]*elem.pos(3)[1]; reg25=reg44-reg25; reg85=reg85*reg41; reg86=reg86*reg74; reg37=reg94+reg37;
    reg93=reg24+reg93; reg74=reg29*reg74; reg24=reg91*reg40; reg41=reg15*reg41; reg30=reg9-reg30;
    reg9=reg21*reg31; reg81=reg81-reg69; reg16=reg16-reg73; reg7=reg39-reg7; reg87=reg87*reg36;
    reg53=reg80/reg53; reg15=reg31*reg25; reg29=reg78*reg31; reg111=reg111*reg36; reg86=reg85-reg86;
    reg39=reg58*reg40; reg15=reg7+reg15; reg87=reg30-reg87; reg9=reg24-reg9; reg7=reg12*reg31;
    reg24=reg61*reg31; reg30=1-(*f.m).resolution; reg44=reg16*reg37; reg76=reg81*reg93; reg80=reg91*reg78;
    reg82=reg91*reg95; reg85=reg58*reg95; reg36=reg53*reg36; reg74=reg41-reg74; reg29=reg39-reg29;
    reg86=reg111+reg86; reg39=reg21*reg58; reg41=(*f.m).resolution*reg103; reg24=reg82-reg24; reg39=reg80-reg39;
    reg87=reg30*reg87; reg7=reg85-reg7; reg29=reg29/reg15; reg44=reg76-reg44; reg74=reg36+reg74;
    reg9=reg9/reg15; reg86=reg30*reg86; reg36=(*f.m).resolution*reg77; reg53=reg61*reg58; reg76=reg91*reg12;
    reg80=(*f.m).resolution*reg29; reg82=(*f.m).resolution*reg9; reg53=reg76-reg53; reg62=reg62/reg15; reg93=reg93/reg44;
    reg98=reg98*reg30; reg81=reg81/reg44; reg37=reg37/reg44; reg16=reg16/reg44; reg65=reg65*reg30;
    reg87=reg41+reg87; reg7=reg7/reg15; reg86=reg36+reg86; reg8=reg8/reg15; reg39=reg39/reg15;
    reg74=reg30*reg74; reg24=reg24/reg15; reg36=(*f.m).resolution*reg66; reg74=reg36+reg74; reg86=(*f.m).deltaT*reg86;
    reg82=reg65-reg82; reg84=reg84*reg30; reg22=reg22*reg30; reg36=var_inter[1]*reg93; reg41=var_inter[1]*reg37;
    reg71=reg71*reg30; reg65=(*f.m).resolution*reg39; reg76=(*f.m).resolution*reg24; reg85=(*f.m).resolution*reg7; reg87=(*f.m).deltaT*reg87;
    reg17=reg17*reg30; reg94=var_inter[0]*reg16; reg53=reg53/reg15; reg89=reg89*reg30; reg102=reg11*reg81;
    reg15=reg25/reg15; reg25=reg23*reg37; reg104=reg11*reg16; reg107=(*f.m).resolution*reg8; reg108=(*f.m).resolution*reg62;
    reg110=var_inter[0]*reg81; reg111=reg23*reg93; reg98=reg80+reg98; reg74=(*f.m).deltaT*reg74; reg18=reg18*reg30;
    reg80=reg41+reg102; reg113=reg36-reg94; reg115=(*f.m).resolution*reg53; reg116=reg110-reg41; reg30=reg83*reg30;
    reg83=reg36+reg104; reg117=reg25-reg102; reg119=(*f.m).resolution*reg15; reg121=reg111+reg94; reg122=reg104-reg111;
    reg71=reg107+reg71; reg108=reg84-reg108; reg85=reg22-reg85; reg22=reg98*reg86; reg84=reg82*reg87;
    reg107=reg110+reg25; reg17=reg76+reg17; reg89=reg65+reg89; reg18=reg119+reg18; reg115=reg30-reg115;
    reg30=reg85*reg86; reg65=reg108*reg87; reg76=reg71*reg86; reg119=reg22+reg84; reg123=reg89*reg74;
    reg124=0.5*reg107; reg126=0.5*reg122; reg127=0.5*reg117; reg130=0.5*reg80; reg131=0.5*reg83;
    reg132=0.5*reg116; reg133=0.5*reg113; reg134=0.5*reg121; reg137=reg17*reg87; reg138=reg132*reg89;
    reg139=reg113*reg98; reg140=reg121*reg98; reg141=reg124*reg89; reg142=reg126*reg89; T reg143=reg117*reg82;
    T reg144=reg107*reg82; T reg145=reg134*reg89; T reg146=reg127*reg89; T reg147=reg122*reg98; T reg148=reg30+reg137;
    T reg149=reg115*reg74; T reg150=reg65+reg76; T reg151=reg18*reg74; T reg152=reg131*reg89; T reg153=reg119+reg123;
    T reg154=reg80*reg82; T reg155=reg133*reg89; T reg156=reg116*reg82; T reg157=reg130*reg89; T reg158=reg83*reg98;
    T reg159=reg126*reg18; T reg160=reg130*reg115; T reg161=reg83*reg85; T reg162=reg117*reg108; T reg163=reg133*reg115;
    T reg164=reg116*reg17; reg146=reg147+reg146; reg147=reg127*reg18; T reg165=reg122*reg71; T reg166=reg132*reg115;
    T reg167=reg113*reg85; T reg168=reg117*reg17; T reg169=reg83*reg71; T reg170=reg134*reg115; T reg171=reg107*reg17;
    T reg172=reg150+reg151; T reg173=reg126*reg115; T reg174=reg124*reg115; T reg175=reg121*reg85; T reg176=reg130*reg18;
    T reg177=2*reg153; reg157=reg157-reg158; T reg178=reg131*reg18; T reg179=reg80*reg108; reg154=reg154-reg152;
    T reg180=reg148+reg149; T reg181=reg80*reg17; T reg182=reg132*reg18; T reg183=reg116*reg108; T reg184=reg131*reg115;
    T reg185=reg124*reg18; T reg186=reg107*reg108; reg142=reg143+reg142; reg143=reg113*reg71; T reg187=reg133*reg18;
    reg138=reg139+reg138; reg139=reg134*reg18; reg140=reg140-reg141; reg155=reg156+reg155; reg156=reg121*reg71;
    reg145=reg145-reg144; reg187=reg183+reg187; reg183=reg107*reg180; reg170=reg170-reg171; reg175=reg175-reg174;
    reg176=reg176-reg169; reg182=reg143+reg182; reg143=reg131*reg177; reg173=reg168+reg173; reg156=reg156-reg185;
    reg157=2*reg157; reg179=reg179-reg178; reg168=reg134*reg177; reg138=2*reg138; reg154=2*reg154;
    T reg188=reg121*reg172; T reg189=reg124*reg177; reg145=2*reg145; reg159=reg162+reg159; reg160=reg160-reg161;
    reg139=reg139-reg186; reg163=reg164+reg163; reg162=reg130*reg177; reg146=2*reg146; reg155=2*reg155;
    reg181=reg181-reg184; reg142=2*reg142; reg140=2*reg140; reg164=reg80*reg180; reg166=reg167+reg166;
    reg147=reg165+reg147; reg165=reg83*reg172; reg167=reg121*reg156; T reg190=reg121*reg139; T reg191=reg124*reg145;
    T reg192=reg124*reg138; T reg193=reg124*reg157; T reg194=reg121*reg187; T reg195=reg121*reg176; T reg196=reg121*reg179;
    T reg197=reg124*reg154; T reg198=reg134*reg145; T reg199=reg124*reg140; T reg200=reg107*reg170; T reg201=reg134*reg138;
    T reg202=reg124*reg155; T reg203=reg121*reg182; T reg204=reg113*reg179; T reg205=reg132*reg154; T reg206=reg133*reg177;
    T reg207=reg113*reg176; T reg208=reg165-reg162; T reg209=reg132*reg157; T reg210=reg143-reg164; T reg211=reg113*reg187;
    T reg212=reg132*reg155; T reg213=reg113*reg182; T reg214=reg132*reg138; T reg215=reg107*reg181; T reg216=reg134*reg154;
    T reg217=reg107*reg160; T reg218=reg107*reg163; T reg219=reg134*reg155; T reg220=reg107*reg166; T reg221=reg117*reg180;
    T reg222=reg189-reg188; T reg223=reg126*reg177; T reg224=reg127*reg177; T reg225=reg122*reg172; T reg226=reg131*reg154;
    T reg227=reg80*reg181; T reg228=reg83*reg179; T reg229=reg130*reg154; T reg230=reg83*reg176; T reg231=reg130*reg157;
    T reg232=reg183-reg168; T reg233=reg132*reg177; T reg234=reg113*reg172; T reg235=reg116*reg180; T reg236=reg133*reg154;
    T reg237=reg116*reg181; T reg238=reg133*reg157; T reg239=reg116*reg160; T reg240=reg133*reg155; T reg241=reg116*reg163;
    T reg242=reg127*reg142; T reg243=reg122*reg159; T reg244=reg127*reg145; T reg245=reg126*reg155; T reg246=reg122*reg182;
    T reg247=reg127*reg138; T reg248=reg117*reg173; T reg249=reg117*reg166; T reg250=reg126*reg142; T reg251=reg127*reg146;
    T reg252=reg126*reg154; T reg253=reg122*reg147; reg176=reg122*reg176; T reg254=reg126*reg138; T reg255=reg127*reg157;
    reg154=reg127*reg154; reg179=reg122*reg179; T reg256=reg122*reg187; T reg257=reg117*reg163; T reg258=reg127*reg155;
    T reg259=reg122*reg139; T reg260=reg134*reg157; T reg261=reg117*reg170; T reg262=reg126*reg140; T reg263=reg117*reg175;
    T reg264=reg126*reg157; reg181=reg117*reg181; T reg265=reg122*reg156; T reg266=reg117*reg160; T reg267=reg126*reg145;
    T reg268=reg127*reg140; T reg269=reg223+reg221; reg261=reg267+reg261; reg240=reg241+reg240; reg215=reg216-reg215;
    reg154=reg179+reg154; reg204=reg205+reg204; reg247=reg246+reg247; reg213=reg214+reg213; reg266=reg264+reg266;
    reg263=reg262+reg263; reg207=reg209+reg207; reg211=reg212+reg211; reg238=reg239+reg238; reg248=reg250+reg248;
    reg226=reg227-reg226; reg179=reg225+reg224; reg199=reg167-reg199; reg210=reg44*reg210; reg181=reg252+reg181;
    reg191=reg190-reg191; reg249=reg254+reg249; reg208=reg44*reg208; reg192=reg203-reg192; reg268=reg265+reg268;
    reg230=reg231-reg230; reg228=reg229-reg228; reg242=reg243+reg242; reg202=reg194-reg202; reg167=reg235+reg206;
    reg193=reg195-reg193; reg257=reg245+reg257; reg217=reg260-reg217; reg255=reg176+reg255; reg218=reg219-reg218;
    reg258=reg256+reg258; reg222=reg44*reg222; reg236=reg237+reg236; reg220=reg201-reg220; reg232=reg44*reg232;
    reg200=reg198-reg200; reg251=reg253+reg251; reg197=reg196-reg197; reg244=reg259+reg244; reg176=reg233+reg234;
    reg226=reg44*reg226; reg228=reg44*reg228; reg249=reg44*reg249; reg261=reg44*reg261; reg263=reg44*reg263;
    reg248=reg44*reg248; reg190=reg44*reg179; reg154=reg44*reg154; reg194=reg44*reg269; reg255=reg44*reg255;
    reg222=ponderation*reg222; reg232=ponderation*reg232; reg244=reg44*reg244; reg195=reg44*reg176; reg196=reg44*reg167;
    reg268=reg44*reg268; reg208=ponderation*reg208; reg210=ponderation*reg210; reg220=reg44*reg220; reg258=reg44*reg258;
    reg218=reg44*reg218; reg200=reg44*reg200; reg217=reg44*reg217; reg197=reg44*reg197; reg251=reg44*reg251;
    reg215=reg44*reg215; reg247=reg44*reg247; reg213=reg44*reg213; reg193=reg44*reg193; reg211=reg44*reg211;
    reg202=reg44*reg202; reg242=reg44*reg242; reg207=reg44*reg207; reg266=reg44*reg266; reg230=reg44*reg230;
    reg199=reg44*reg199; reg236=reg44*reg236; reg181=reg44*reg181; reg238=reg44*reg238; reg257=reg44*reg257;
    reg191=reg44*reg191; reg204=reg44*reg204; reg240=reg44*reg240; reg192=reg44*reg192; sollicitation[indices[3]+1]+=-reg210;
    matrix(indices[0]+0,indices[2]+1)+=ponderation*reg258; matrix(indices[1]+1,indices[1]+1)+=ponderation*reg200; matrix(indices[0]+1,indices[3]+1)+=ponderation*reg181; matrix(indices[1]+0,indices[2]+0)+=ponderation*reg192; matrix(indices[1]+0,indices[3]+1)+=ponderation*reg197;
    sollicitation[indices[1]+1]+=-reg232; matrix(indices[1]+0,indices[1]+0)+=ponderation*reg199; matrix(indices[0]+0,indices[1]+0)+=ponderation*reg268; reg181=ponderation*reg196; sollicitation[indices[2]+1]+=reg181;
    matrix(indices[0]+0,indices[0]+0)+=ponderation*reg251; matrix(indices[1]+0,indices[3]+0)+=ponderation*reg193; matrix(indices[0]+0,indices[0]+1)+=ponderation*reg242; reg192=ponderation*reg195; sollicitation[indices[2]+0]+=reg192;
    sollicitation[indices[3]+0]+=-reg208; matrix(indices[1]+0,indices[1]+1)+=ponderation*reg191; matrix(indices[1]+0,indices[2]+1)+=ponderation*reg202; matrix(indices[3]+0,indices[3]+0)+=ponderation*reg230; matrix(indices[0]+1,indices[2]+0)+=ponderation*reg249;
    matrix(indices[2]+1,indices[3]+1)+=ponderation*reg236; matrix(indices[3]+0,indices[3]+1)+=ponderation*reg228; matrix(indices[2]+1,indices[3]+0)+=ponderation*reg238; matrix(indices[3]+1,indices[3]+1)+=ponderation*reg226; matrix(indices[0]+1,indices[1]+1)+=ponderation*reg261;
    matrix(indices[0]+1,indices[2]+1)+=ponderation*reg257; matrix(indices[2]+1,indices[2]+1)+=ponderation*reg240; matrix(indices[2]+0,indices[3]+1)+=ponderation*reg204; matrix(indices[0]+1,indices[1]+0)+=ponderation*reg263; matrix(indices[2]+0,indices[3]+0)+=ponderation*reg207;
    matrix(indices[0]+1,indices[0]+1)+=ponderation*reg248; matrix(indices[0]+1,indices[3]+0)+=ponderation*reg266; matrix(indices[0]+0,indices[1]+1)+=ponderation*reg244; matrix(indices[2]+0,indices[2]+1)+=ponderation*reg211; reg191=ponderation*reg190;
    sollicitation[indices[0]+0]+=reg191; matrix(indices[2]+0,indices[2]+0)+=ponderation*reg213; matrix(indices[0]+0,indices[3]+1)+=ponderation*reg154; matrix(indices[1]+1,indices[3]+1)+=ponderation*reg215; matrix(indices[0]+0,indices[2]+0)+=ponderation*reg247;
    matrix(indices[1]+1,indices[3]+0)+=ponderation*reg217; reg154=ponderation*reg194; sollicitation[indices[0]+1]+=reg154; matrix(indices[0]+0,indices[3]+0)+=ponderation*reg255; matrix(indices[1]+1,indices[2]+1)+=ponderation*reg218;
    matrix(indices[1]+1,indices[2]+0)+=ponderation*reg220; sollicitation[indices[1]+0]+=-reg222;
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
    T reg0=pow((*f.m).v1[0],2); T reg1=pow((*f.m).v1[1],2); reg1=reg0+reg1; reg0=pow((*f.m).v1[2],2); T reg2=pow((*f.m).v2[1],2);
    T reg3=pow((*f.m).v2[0],2); reg0=reg1+reg0; reg1=pow((*f.m).v2[2],2); reg2=reg3+reg2; reg3=2*(*f.m).shear_modulus_13;
    T reg4=2*(*f.m).shear_modulus_23; reg4=1.0/reg4; reg1=reg2+reg1; reg2=2*(*f.m).shear_modulus_12; reg3=1.0/reg3;
    reg0=pow(reg0,0.5); T reg5=(*f.m).v1[1]/reg0; T reg6=reg3*reg4; reg1=pow(reg1,0.5); T reg7=(*f.m).v1[0]/reg0;
    reg2=1.0/reg2; reg0=(*f.m).v1[2]/reg0; T reg8=(*f.m).v2[0]/reg1; T reg9=(*f.m).v2[1]/reg1; T reg10=2*reg7;
    T reg11=2*reg5; T reg12=reg2*reg6; T reg13=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg14=1.0/(*f.m).elastic_modulus_3; T reg15=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg16=reg15*reg12; T reg17=reg13*reg12; T reg18=reg14*reg12; T reg19=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg20=1.0/(*f.m).elastic_modulus_2;
    T reg21=1.0/(*f.m).elastic_modulus_1; T reg22=2*reg0; T reg23=reg9*reg11; T reg24=reg8*reg10; T reg25=pow(reg9,2);
    T reg26=pow(reg8,2); reg1=(*f.m).v2[2]/reg1; T reg27=reg24*reg21; T reg28=reg23*reg20; T reg29=pow(reg1,2);
    T reg30=reg25*reg19; T reg31=reg26*reg19; T reg32=reg23*reg19; reg22=reg1*reg22; T reg33=reg13*reg16;
    T reg34=reg26*reg21; T reg35=reg20*reg18; T reg36=reg24*reg19; T reg37=reg13*reg17; T reg38=reg19*reg18;
    T reg39=reg25*reg20; T reg40=reg20*reg16; T reg41=reg19*reg17; reg36=reg28-reg36; reg28=reg23*reg13;
    T reg42=reg26*reg15; T reg43=pow(reg5,2); T reg44=pow(reg7,2); T reg45=reg22*reg13; reg31=reg39-reg31;
    reg39=reg29*reg13; reg32=reg27-reg32; reg37=reg35-reg37; reg27=reg24*reg15; reg35=reg22*reg15;
    reg33=reg38+reg33; reg30=reg34-reg30; reg34=reg29*reg15; T reg46=reg25*reg13; T reg47=reg44*reg21;
    T reg48=reg19*reg33; T reg49=reg5*reg8; T reg50=reg7*reg9; T reg51=reg21*reg37; T reg52=reg43*reg19;
    T reg53=pow(reg0,2); T reg54=reg29*reg14; T reg55=reg43*reg20; T reg56=reg7*reg8; reg34=reg30-reg34;
    reg46=reg42+reg46; reg30=reg5*reg9; reg42=reg41+reg40; reg28=reg27+reg28; reg22=reg22*reg14;
    reg35=reg32-reg35; reg27=reg44*reg19; reg45=reg36-reg45; reg39=reg31-reg39; reg31=reg26*reg35;
    reg32=reg0*reg1; reg36=reg49+reg50; T reg57=reg19*reg12; T reg58=reg15*reg16; reg18=reg21*reg18;
    T reg59=reg15*reg17; reg12=reg20*reg12; T reg60=reg53*reg13; T reg61=reg56*reg35; T reg62=reg30*reg45;
    T reg63=reg8*reg9; T reg64=reg15*reg42; T reg65=reg25*reg39; T reg66=reg26*reg34; T reg67=reg44*reg15;
    T reg68=reg43*reg39; reg28=reg22-reg28; reg22=reg15*reg6; T reg69=reg44*reg35; T reg70=reg43*reg45;
    reg46=reg54-reg46; reg48=reg51-reg48; reg52=reg47-reg52; reg47=reg53*reg15; reg51=reg43*reg13;
    reg54=reg14*reg6; T reg71=reg7*reg1; T reg72=reg0*reg8; reg6=reg13*reg6; T reg73=reg0*reg9;
    T reg74=reg5*reg1; T reg75=reg30*reg39; reg27=reg55-reg27; reg55=reg25*reg45; T reg76=2*reg8;
    T reg77=reg44*reg34; T reg78=reg2*reg4; T reg79=reg56*reg34; T reg80=reg32*reg46; T reg81=reg8*reg1;
    reg60=reg27-reg60; reg27=reg53*reg14; reg62=reg61+reg62; reg61=reg32*reg28; T reg82=reg74-reg73;
    T reg83=reg14*reg78; reg68=reg77+reg68; reg77=reg53*reg46; reg64=reg48-reg64; reg48=reg2*reg3;
    reg51=reg67+reg51; reg67=reg15*reg78; T reg84=reg53*reg28; T reg85=reg15*reg57; reg47=reg52-reg47;
    reg70=reg69+reg70; reg17=reg21*reg17; reg58=reg18-reg58; reg18=reg9*reg76; reg52=reg71+reg72;
    reg69=reg13*reg6; reg16=reg19*reg16; T reg86=reg19*reg54; reg59=reg38+reg59; reg55=reg31+reg55;
    reg31=reg15*reg12; reg54=reg20*reg54; reg38=reg36*reg2; T reg87=reg29*reg28; T reg88=reg63*reg2;
    reg78=reg13*reg78; reg75=reg79+reg75; reg79=reg13*reg22; T reg89=reg29*reg46; T reg90=reg5*reg10;
    reg65=reg66+reg65; reg51=reg27-reg51; reg69=reg54-reg69; reg27=reg13*reg67; reg54=reg15*reg48;
    reg22=reg20*reg22; reg66=reg13*reg48; reg6=reg19*reg6; T reg91=reg13*reg78; T reg92=reg19*reg83;
    reg83=reg20*reg83; reg48=reg14*reg48; reg57=reg19*reg57; reg79=reg86+reg79; reg33=reg33/reg64;
    reg16=reg17+reg16; reg31=reg41+reg31; reg58=reg58/reg64; reg12=reg21*reg12; reg85=reg17+reg85;
    reg14=2*reg82; reg17=2*reg9; reg86=reg1*reg76; T reg93=reg0*reg10; reg71=reg72-reg71;
    reg77=reg68+reg77; reg68=reg90*reg88; reg84=reg70+reg84; reg70=reg90*reg38; reg72=reg26*reg47;
    T reg94=reg25*reg60; reg89=reg65+reg89; reg65=reg18*reg88; reg87=reg55+reg87; reg55=reg18*reg38;
    reg80=reg75+reg80; reg75=reg36*reg88; reg61=reg62+reg61; reg62=reg36*reg38; reg74=reg73+reg74;
    reg59=reg59/reg64; reg73=reg44*reg47; T reg95=reg81*reg3; T reg96=reg52*reg3; T reg97=reg7*reg5;
    reg37=reg37/reg64; T reg98=reg9*reg1; T reg99=reg43*reg60; T reg100=reg20*reg48; T reg101=reg7*reg0;
    T reg102=reg44*reg59; T reg103=reg26*reg58; reg27=reg92+reg27; reg92=reg30*reg60; T reg104=reg56*reg47;
    T reg105=1-var_inter[1]; reg91=reg83-reg91; reg83=reg86*reg96; reg55=reg87+reg55; reg2=reg97*reg2;
    reg87=reg43*reg59; T reg106=1-var_inter[0]; T reg107=reg18*reg33; T reg108=reg90*reg37; T reg109=reg52*reg96;
    T reg110=reg25*reg33; T reg111=reg43*reg37; reg62=reg61+reg62; reg61=reg26*reg33; reg42=reg42/reg64;
    T reg112=reg13*reg54; T reg113=reg52*reg95; reg75=reg80+reg75; reg13=reg13*reg66; reg48=reg19*reg48;
    reg80=reg71*reg14; T reg114=pow(reg71,2); reg85=reg85/reg64; T reg115=pow(reg82,2); T reg116=reg1*reg17;
    reg16=reg16/reg64; T reg117=reg0*reg11; T reg118=reg53*reg51; reg99=reg73+reg99; reg57=reg12-reg57;
    reg12=reg44*reg37; reg68=reg77+reg68; reg73=reg93*reg95; reg77=reg18*reg58; reg78=reg19*reg78;
    T reg119=reg74*reg4; T reg120=reg25*reg58; T reg121=reg86*reg95; T reg122=reg6+reg22; reg65=reg89+reg65;
    reg79=reg19*reg79; reg89=reg29*reg51; reg69=reg21*reg69; reg94=reg72+reg94; reg72=reg98*reg4;
    reg31=reg31/reg64; T reg123=reg90*reg59; T reg124=reg93*reg96; reg70=reg84+reg70; reg67=reg20*reg67;
    reg84=reg25*reg85; reg109=reg62+reg109; reg62=reg74*reg119; T reg125=reg90*reg31; T reg126=reg18*reg85;
    T reg127=reg26*reg85; T reg128=reg44*reg31; T reg129=reg43*reg31; T reg130=reg90*reg2; reg118=reg99+reg118;
    reg73=reg68+reg73; reg68=reg117*reg72; reg124=reg70+reg124; reg70=reg117*reg119; reg89=reg94+reg89;
    reg94=reg18*reg2; reg121=reg65+reg121; reg65=reg116*reg72; reg3=reg101*reg3; reg83=reg55+reg83;
    reg55=reg116*reg119; reg92=reg104+reg92; reg99=reg32*reg51; reg113=reg75+reg113; reg75=reg74*reg72;
    reg79=reg69-reg79; reg122=reg15*reg122; reg69=reg105*elem.pos(1)[1]; reg91=reg21*reg91; reg104=reg105*elem.pos(0)[1];
    reg27=reg19*reg27; T reg131=reg78+reg67; reg13=reg100-reg13; reg112=reg48+reg112; reg48=elem.pos(1)[0]*var_inter[0];
    reg110=reg111+reg110; reg100=reg114*reg42; reg111=reg106*elem.pos(0)[0]; reg107=reg108+reg107; reg108=reg80*reg42;
    reg103=reg102+reg103; reg102=reg115*reg16; T reg132=reg106*elem.pos(0)[1]; T reg133=elem.pos(1)[1]*var_inter[0]; reg77=reg123+reg77;
    reg123=reg105*elem.pos(1)[0]; T reg134=reg105*elem.pos(0)[0]; T reg135=reg114*reg16; reg120=reg87+reg120; reg87=reg80*reg16;
    T reg136=reg5*reg0; reg57=reg57/reg64; reg12=reg61+reg12; reg61=reg115*reg42; reg66=reg19*reg66;
    reg54=reg20*reg54; reg55=reg83+reg55; reg62=reg109+reg62; reg20=reg48+reg111; reg123=reg123-reg134;
    reg99=reg92+reg99; reg83=reg36*reg2; reg92=var_inter[1]*elem.pos(2)[1]; reg109=var_inter[1]*elem.pos(2)[0]; reg69=reg69-reg104;
    reg75=reg113+reg75; reg70=reg124+reg70; reg94=reg89+reg94; reg89=reg86*reg3; reg68=reg73+reg68;
    reg73=elem.pos(2)[0]*var_inter[0]; reg65=reg121+reg65; reg113=reg93*reg3; reg130=reg118+reg130; reg4=reg136*reg4;
    reg118=elem.pos(2)[1]*var_inter[0]; reg122=reg79-reg122; reg27=reg91-reg27; reg131=reg15*reg131; reg13=reg21*reg13;
    reg112=reg19*reg112; reg19=reg66+reg54; reg61=reg12+reg61; reg100=reg110+reg100; reg108=reg107+reg108;
    reg12=reg133+reg132; reg87=reg77+reg87; reg127=reg128+reg127; reg135=reg120+reg135; reg102=reg103+reg102;
    reg21=reg114*reg57; reg84=reg129+reg84; reg77=reg115*reg57; reg126=reg125+reg126; reg79=reg80*reg57;
    reg91=reg52*reg3; reg92=reg69+reg92; reg69=reg7*reg71; reg103=reg5*reg82; reg107=var_inter[1]*elem.pos(3)[1];
    reg21=reg84+reg21; reg84=reg63*reg87; reg110=reg97*reg108; reg120=var_inter[1]*elem.pos(3)[0]; reg109=reg123+reg109;
    reg113=reg130+reg113; reg121=reg63*reg135; reg123=reg97*reg100; reg83=reg99+reg83; reg73=reg73-reg20;
    reg99=elem.pos(3)[0]*reg106; reg124=elem.pos(3)[1]*reg106; reg118=reg118-reg12; reg125=reg71*reg82; reg19=reg15*reg19;
    reg79=reg126+reg79; reg15=reg68*reg62; reg112=reg13-reg112; reg13=reg65*reg62; reg131=reg27-reg131;
    reg27=reg55*reg75; reg122=reg122/reg64; reg126=reg70*reg75; reg128=reg116*reg4; reg77=reg127+reg77;
    reg89=reg94+reg89; reg94=reg97*reg61; reg127=reg63*reg102; reg129=reg117*reg4; reg91=reg83+reg91;
    reg83=reg68*reg55; reg130=reg43*reg100; reg131=reg131/reg64; T reg137=reg26*reg87; T reg138=reg44*reg108;
    T reg139=reg25*reg135; T reg140=reg44*reg61; T reg141=reg26*reg102; reg99=reg73+reg99; reg73=reg9*reg82;
    reg87=reg25*reg87; T reg142=reg8*reg71; reg121=reg123+reg121; reg135=reg26*reg135; reg100=reg44*reg100;
    reg108=reg43*reg108; reg126=reg15-reg126; reg129=reg113+reg129; reg19=reg112-reg19; reg128=reg89+reg128;
    reg15=reg70*reg65; reg118=reg124+reg118; reg89=reg56*reg122; reg112=reg30*reg122; reg113=reg36*reg122;
    reg109=reg109-reg120; reg123=reg125*reg77; reg127=reg94+reg127; reg84=reg110+reg84; reg94=reg125*reg79;
    reg110=reg74*reg4; reg27=reg13-reg27; reg92=reg92-reg107; reg69=reg103+reg69; reg102=reg25*reg102;
    reg61=reg43*reg61; reg13=reg5*reg71; reg103=reg125*reg21; reg124=reg7*reg82; T reg143=reg115*reg77;
    T reg144=reg114*reg21; reg64=reg19/reg64; reg19=reg124*reg131; T reg145=reg13*reg131; reg139=reg130+reg139;
    reg130=reg27*reg129; T reg146=reg69*reg131; T reg147=reg128*reg126; reg141=reg140+reg141; reg103=reg121+reg103;
    reg121=reg36*reg112; reg140=reg109*reg118; T reg148=reg36*reg89; reg123=reg127+reg123; reg94=reg84+reg94;
    reg84=reg36*reg113; reg127=reg92*reg99; reg110=reg91+reg110; reg102=reg61+reg102; reg77=reg114*reg77;
    reg61=reg115*reg79; reg137=reg138+reg137; reg91=reg8*reg82; reg138=reg9*reg71; reg79=reg114*reg79;
    reg87=reg108+reg87; reg142=reg73+reg142; reg21=reg115*reg21; reg135=reg100+reg135; reg15=reg83-reg15;
    reg127=reg140-reg127; reg10=reg82*reg10; reg11=reg71*reg11; reg73=reg24*reg113; reg83=reg65*reg110;
    reg100=reg69*reg146; reg84=reg94+reg84; reg94=reg69*reg145; reg121=reg103+reg121; reg103=reg128*reg75;
    reg108=reg69*reg19; reg148=reg123+reg148; reg113=reg23*reg113; reg79=reg87+reg79; reg77=reg102+reg77;
    reg87=reg23*reg89; reg102=reg23*reg112; reg144=reg139+reg144; reg147=reg130-reg147; reg123=reg110*reg15;
    reg130=reg91*reg64; reg139=reg138*reg64; reg140=reg142*reg64; reg143=reg141+reg143; reg89=reg24*reg89;
    reg141=reg68*reg128; T reg149=reg129*reg65; reg21=reg135+reg21; reg61=reg137+reg61; reg135=reg129*reg75;
    reg137=reg68*reg110; reg112=reg24*reg112; T reg150=reg129*reg55; T reg151=reg70*reg128; reg118=reg118/reg127;
    reg109=reg109/reg127; reg137=reg135-reg137; reg99=reg99/reg127; reg92=reg92/reg127; reg135=reg128*reg62;
    reg141=reg149-reg141; reg149=reg70*reg110; reg123=reg147+reg123; reg147=reg129*reg62; T reg152=reg55*reg110;
    reg83=reg103-reg83; reg103=reg11*reg145; reg102=reg144+reg102; reg145=reg10*reg145; reg112=reg21+reg112;
    reg76=reg82*reg76; reg17=reg71*reg17; reg21=reg10*reg19; reg89=reg143+reg89; reg73=reg61+reg73;
    reg61=reg10*reg146; reg19=reg11*reg19; reg143=reg142*reg140; reg100=reg84+reg100; reg87=reg77+reg87;
    reg77=reg142*reg139; reg94=reg121+reg94; reg84=reg142*reg130; reg108=reg148+reg108; reg146=reg11*reg146;
    reg113=reg79+reg113; reg83=reg83/reg123; reg79=reg76*reg139; reg152=reg135-reg152; reg145=reg112+reg145;
    reg61=reg73+reg61; reg73=reg76*reg140; reg149=reg147-reg149; reg137=reg137/reg123; reg112=reg76*reg130;
    reg21=reg89+reg21; reg151=reg150-reg151; reg141=reg141/reg123; reg89=1-(*f.m).resolution; reg103=reg102+reg103;
    reg139=reg17*reg139; reg143=reg100+reg143; reg77=reg94+reg77; reg84=reg108+reg84; reg140=reg17*reg140;
    reg146=reg113+reg146; reg94=var_inter[1]*reg99; reg100=reg105*reg118; reg19=reg87+reg19; reg130=reg17*reg130;
    reg87=var_inter[1]*reg118; reg102=reg106*reg109; reg108=var_inter[0]*reg109; reg113=reg105*reg99; reg121=reg106*reg92;
    reg135=var_inter[0]*reg92; reg144=reg143*reg89; reg147=reg77*reg89; reg148=reg84*reg89; reg150=reg100+reg135;
    T reg153=(*f.m).resolution*reg141; T reg154=(*f.m).resolution*reg137; T reg155=reg113-reg102; reg139=reg103+reg139; reg103=(*f.m).resolution*reg83;
    T reg156=reg87+reg121; T reg157=reg108+reg113; T reg158=reg94+reg102; reg140=reg146+reg140; reg146=reg87-reg135;
    T reg159=reg108-reg94; reg152=reg152/reg123; reg73=reg61+reg73; reg27=reg27/reg123; reg130=reg19+reg130;
    reg79=reg145+reg79; reg126=reg126/reg123; reg149=reg149/reg123; reg15=reg15/reg123; reg19=reg121-reg100;
    reg123=reg151/reg123; reg112=reg21+reg112; reg21=0.5*reg155; reg61=0.5*reg150; reg145=(*f.m).resolution*reg149;
    reg151=0.5*reg159; T reg160=0.5*reg156; T reg161=0.5*reg158; T reg162=(*f.m).resolution*reg126; T reg163=0.5*reg19;
    T reg164=0.5*reg157; T reg165=0.5*reg146; T reg166=(*f.m).resolution*reg123; T reg167=reg140*reg89; T reg168=reg139*reg89;
    reg144=reg153+reg144; reg154=reg147-reg154; reg148=reg103+reg148; reg103=reg112*reg89; reg147=reg79*reg89;
    reg153=(*f.m).resolution*reg152; T reg169=(*f.m).resolution*reg15; T reg170=reg73*reg89; T reg171=reg130*reg89; T reg172=(*f.m).resolution*reg27;
    T reg173=reg151*reg144; T reg174=reg165*reg144; T reg175=reg146*reg148; T reg176=reg159*reg154; T reg177=reg161*reg144;
    T reg178=reg156*reg148; T reg179=reg21*reg144; T reg180=reg19*reg148; T reg181=reg61*reg144; T reg182=reg155*reg154;
    T reg183=reg158*reg154; T reg184=reg160*reg144; T reg185=reg163*reg144; reg166=reg167-reg166; reg168=reg145+reg168;
    reg153=reg171-reg153; reg145=reg150*reg148; reg167=reg164*reg144; reg171=reg157*reg154; reg103=reg172+reg103;
    reg162=reg147-reg162; reg170=reg169+reg170; reg147=reg150*reg103; reg169=reg160*reg166; reg172=reg156*reg153;
    T reg186=reg161*reg166; reg174=reg176+reg174; reg176=reg155*reg162; T reg187=reg163*reg170; T reg188=reg165*reg170;
    reg181=reg181-reg171; reg173=reg175+reg173; reg175=reg146*reg103; T reg189=reg151*reg170; T reg190=reg159*reg162;
    T reg191=reg156*reg103; T reg192=reg161*reg170; T reg193=reg164*reg170; T reg194=reg158*reg168; reg185=reg182+reg185;
    reg145=reg145-reg167; reg182=reg157*reg162; T reg195=reg61*reg170; T reg196=reg163*reg166; T reg197=reg150*reg153;
    T reg198=reg164*reg166; T reg199=reg155*reg168; T reg200=reg21*reg166; T reg201=reg19*reg153; T reg202=reg157*reg168;
    T reg203=reg61*reg166; reg183=reg183-reg184; T reg204=reg146*reg153; T reg205=reg151*reg166; T reg206=reg19*reg103;
    T reg207=reg21*reg170; T reg208=reg160*reg170; T reg209=reg158*reg162; T reg210=reg159*reg168; T reg211=reg165*reg166;
    reg179=reg180+reg179; reg177=reg177-reg178; reg185=2*reg185; reg187=reg176+reg187; reg207=reg206+reg207;
    reg174=2*reg174; reg179=2*reg179; reg194=reg194-reg169; reg183=2*reg183; reg181=2*reg181;
    reg186=reg186-reg172; reg192=reg192-reg191; reg197=reg197-reg198; reg211=reg210+reg211; reg209=reg209-reg208;
    reg189=reg175+reg189; reg200=reg201+reg200; reg173=2*reg173; reg203=reg203-reg202; reg147=reg147-reg193;
    reg145=2*reg145; reg177=2*reg177; reg205=reg204+reg205; reg196=reg199+reg196; reg195=reg195-reg182;
    reg188=reg190+reg188; reg175=reg159*reg211; reg176=reg150*reg192; reg180=reg165*reg173; reg190=reg165*reg174;
    reg199=reg159*reg205; reg201=reg161*reg183; reg204=reg165*reg181; reg206=reg159*reg203; reg210=reg156*reg209;
    T reg212=reg165*reg145; T reg213=reg159*reg197; T reg214=reg165*reg185; T reg215=reg159*reg196; T reg216=reg159*reg200;
    T reg217=reg146*reg209; T reg218=reg150*reg207; T reg219=reg164*reg179; T reg220=reg150*reg187; T reg221=reg164*reg185;
    T reg222=reg160*reg173; T reg223=reg165*reg183; T reg224=reg156*reg192; T reg225=reg161*reg179; T reg226=reg156*reg207;
    T reg227=reg150*reg147; T reg228=reg164*reg145; T reg229=reg158*reg205; T reg230=reg150*reg195; T reg231=reg161*reg177;
    T reg232=reg161*reg185; T reg233=reg156*reg187; T reg234=reg161*reg145; T reg235=reg164*reg181; T reg236=reg159*reg194;
    T reg237=reg150*reg189; T reg238=reg164*reg173; T reg239=reg150*reg188; T reg240=reg164*reg174; T reg241=reg165*reg177;
    T reg242=reg159*reg186; T reg243=reg160*reg185; T reg244=reg157*reg194; T reg245=reg61*reg183; T reg246=reg157*reg186;
    T reg247=reg158*reg197; T reg248=reg157*reg211; T reg249=reg61*reg174; T reg250=reg157*reg205; T reg251=reg61*reg173;
    T reg252=reg160*reg145; T reg253=reg157*reg203; T reg254=reg61*reg181; T reg255=reg157*reg197; T reg256=reg61*reg145;
    T reg257=reg164*reg177; T reg258=reg157*reg196; T reg259=reg61*reg185; T reg260=reg160*reg181; T reg261=reg150*reg209;
    T reg262=reg164*reg183; T reg263=reg158*reg203; T reg264=reg157*reg200; T reg265=reg61*reg179; T reg266=reg156*reg147;
    T reg267=reg161*reg181; T reg268=reg156*reg195; T reg269=reg161*reg173; T reg270=reg156*reg188; T reg271=reg161*reg174;
    T reg272=reg151*reg183; T reg273=reg146*reg192; T reg274=reg151*reg177; T reg275=reg158*reg200; T reg276=reg146*reg188;
    T reg277=reg151*reg174; T reg278=reg160*reg179; T reg279=reg146*reg189; T reg280=reg151*reg173; T reg281=reg146*reg195;
    T reg282=reg151*reg181; T reg283=reg146*reg147; T reg284=reg151*reg145; T reg285=reg146*reg187; T reg286=reg151*reg185;
    T reg287=reg158*reg196; T reg288=reg146*reg207; T reg289=reg151*reg179; T reg290=reg160*reg174; T reg291=reg21*reg174;
    reg188=reg19*reg188; T reg292=reg21*reg173; T reg293=reg19*reg189; T reg294=reg155*reg186; T reg295=reg163*reg177;
    reg186=reg158*reg186; T reg296=reg155*reg211; reg174=reg163*reg174; T reg297=reg160*reg177; reg205=reg155*reg205;
    reg173=reg163*reg173; reg203=reg155*reg203; T reg298=reg163*reg181; reg197=reg155*reg197; reg147=reg19*reg147;
    T reg299=reg21*reg145; reg195=reg19*reg195; reg181=reg21*reg181; T reg300=reg163*reg183; T reg301=reg165*reg179;
    reg192=reg19*reg192; T reg302=reg21*reg177; reg209=reg19*reg209; T reg303=reg21*reg183; reg183=reg160*reg183;
    T reg304=reg163*reg179; reg200=reg155*reg200; T reg305=reg158*reg194; T reg306=reg163*reg185; reg196=reg155*reg196;
    reg145=reg163*reg145; reg211=reg158*reg211; reg187=reg19*reg187; reg185=reg21*reg185; reg179=reg21*reg179;
    reg194=reg155*reg194; reg207=reg19*reg207; reg177=reg61*reg177; reg189=reg156*reg189; reg200=reg304+reg200;
    reg210=reg201-reg210; reg214=reg215+reg214; reg185=reg187+reg185; reg270=reg271-reg270; reg217=reg272+reg217;
    reg299=reg147+reg299; reg273=reg274+reg273; reg196=reg306+reg196; reg276=reg277+reg276; reg189=reg269-reg189;
    reg235=reg230-reg235; reg197=reg145+reg197; reg279=reg280+reg279; reg278=reg275-reg278; reg281=reg282+reg281;
    reg283=reg284+reg283; reg203=reg298+reg203; reg285=reg286+reg285; reg297=reg186-reg297; reg266=reg234-reg266;
    reg219=reg218-reg219; reg233=reg232-reg233; reg181=reg195+reg181; reg222=reg229-reg222; reg226=reg225-reg226;
    reg216=reg301+reg216; reg223=reg236+reg223; reg194=reg300+reg194; reg224=reg231-reg224; reg241=reg242+reg241;
    reg221=reg220-reg221; reg302=reg192+reg302; reg190=reg175+reg190; reg183=reg305-reg183; reg180=reg199+reg180;
    reg268=reg267-reg268; reg303=reg209+reg303; reg204=reg206+reg204; reg212=reg213+reg212; reg228=reg227-reg228;
    reg253=reg254-reg253; reg252=reg247-reg252; reg294=reg295+reg294; reg240=reg239-reg240; reg250=reg251-reg250;
    reg255=reg256-reg255; reg260=reg263-reg260; reg179=reg207+reg179; reg248=reg249-reg248; reg258=reg259-reg258;
    reg292=reg293+reg292; reg296=reg174+reg296; reg246=reg177-reg246; reg243=reg287-reg243; reg264=reg265-reg264;
    reg290=reg211-reg290; reg288=reg289+reg288; reg205=reg173+reg205; reg238=reg237-reg238; reg257=reg176-reg257;
    reg291=reg188+reg291; reg244=reg245-reg244; reg262=reg261-reg262; reg262=reg127*reg262; reg221=reg127*reg221;
    reg290=reg127*reg290; reg180=reg127*reg180; reg257=reg127*reg257; reg302=reg127*reg302; reg255=reg127*reg255;
    reg222=reg127*reg222; reg219=reg127*reg219; reg189=reg127*reg189; reg190=reg127*reg190; reg224=reg127*reg224;
    reg233=reg127*reg233; reg240=reg127*reg240; reg266=reg127*reg266; reg241=reg127*reg241; reg252=reg127*reg252;
    reg291=reg127*reg291; reg270=reg127*reg270; reg258=reg127*reg258; reg194=reg127*reg194; reg268=reg127*reg268;
    reg223=reg127*reg223; reg264=reg127*reg264; reg216=reg127*reg216; reg183=reg127*reg183; reg292=reg127*reg292;
    reg226=reg127*reg226; reg181=reg127*reg181; reg299=reg127*reg299; reg243=reg127*reg243; reg196=reg127*reg196;
    reg276=reg127*reg276; reg296=reg127*reg296; reg246=reg127*reg246; reg279=reg127*reg279; reg197=reg127*reg197;
    reg297=reg127*reg297; reg260=reg127*reg260; reg281=reg127*reg281; reg238=reg127*reg238; reg235=reg127*reg235;
    reg278=reg127*reg278; reg283=reg127*reg283; reg244=reg127*reg244; reg203=reg127*reg203; reg285=reg127*reg285;
    reg205=reg127*reg205; reg179=reg127*reg179; reg288=reg127*reg288; reg204=reg127*reg204; reg253=reg127*reg253;
    reg303=reg127*reg303; reg185=reg127*reg185; reg212=reg127*reg212; reg294=reg127*reg294; reg214=reg127*reg214;
    reg250=reg127*reg250; reg200=reg127*reg200; reg248=reg127*reg248; reg273=reg127*reg273; reg217=reg127*reg217;
    reg210=reg127*reg210; reg228=reg127*reg228; matrix(indices[3]+1,indices[1]+0)+=ponderation*reg252; matrix(indices[3]+1,indices[1]+1)+=ponderation*reg260; matrix(indices[3]+1,indices[3]+1)+=ponderation*reg183;
    matrix(indices[3]+1,indices[0]+1)+=ponderation*reg243; matrix(indices[3]+0,indices[3]+1)+=ponderation*reg210; matrix(indices[3]+0,indices[3]+0)+=ponderation*reg224; matrix(indices[3]+1,indices[2]+0)+=ponderation*reg222; matrix(indices[3]+1,indices[0]+0)+=ponderation*reg278;
    matrix(indices[3]+1,indices[3]+0)+=ponderation*reg297; matrix(indices[3]+0,indices[2]+1)+=ponderation*reg270; matrix(indices[3]+1,indices[2]+1)+=ponderation*reg290; matrix(indices[1]+1,indices[0]+0)+=ponderation*reg264; matrix(indices[1]+0,indices[3]+1)+=ponderation*reg262;
    matrix(indices[1]+0,indices[3]+0)+=ponderation*reg257; matrix(indices[1]+0,indices[2]+1)+=ponderation*reg240; matrix(indices[1]+0,indices[2]+0)+=ponderation*reg238; matrix(indices[1]+0,indices[1]+1)+=ponderation*reg235; matrix(indices[1]+0,indices[1]+0)+=ponderation*reg228;
    matrix(indices[1]+0,indices[0]+1)+=ponderation*reg221; matrix(indices[1]+0,indices[0]+0)+=ponderation*reg219; matrix(indices[0]+1,indices[3]+1)+=ponderation*reg194; matrix(indices[0]+0,indices[0]+1)+=ponderation*reg185; matrix(indices[0]+0,indices[0]+0)+=ponderation*reg179;
    matrix(indices[0]+0,indices[2]+1)+=ponderation*reg291; matrix(indices[0]+0,indices[2]+0)+=ponderation*reg292; matrix(indices[0]+0,indices[1]+1)+=ponderation*reg181; matrix(indices[0]+1,indices[3]+0)+=ponderation*reg294; matrix(indices[0]+1,indices[2]+1)+=ponderation*reg296;
    matrix(indices[0]+1,indices[2]+0)+=ponderation*reg205; matrix(indices[0]+1,indices[1]+1)+=ponderation*reg203; matrix(indices[0]+1,indices[1]+0)+=ponderation*reg197; matrix(indices[0]+1,indices[0]+1)+=ponderation*reg196; matrix(indices[0]+1,indices[0]+0)+=ponderation*reg200;
    matrix(indices[0]+0,indices[3]+1)+=ponderation*reg303; matrix(indices[0]+0,indices[3]+0)+=ponderation*reg302; matrix(indices[2]+1,indices[0]+0)+=ponderation*reg216; matrix(indices[0]+0,indices[1]+0)+=ponderation*reg299; matrix(indices[3]+0,indices[2]+0)+=ponderation*reg189;
    matrix(indices[3]+0,indices[1]+1)+=ponderation*reg268; matrix(indices[3]+0,indices[1]+0)+=ponderation*reg266; matrix(indices[3]+0,indices[0]+1)+=ponderation*reg233; matrix(indices[3]+0,indices[0]+0)+=ponderation*reg226; matrix(indices[2]+1,indices[3]+1)+=ponderation*reg223;
    matrix(indices[2]+1,indices[3]+0)+=ponderation*reg241; matrix(indices[2]+1,indices[2]+1)+=ponderation*reg190; matrix(indices[2]+1,indices[2]+0)+=ponderation*reg180; matrix(indices[2]+1,indices[1]+1)+=ponderation*reg204; matrix(indices[2]+1,indices[1]+0)+=ponderation*reg212;
    matrix(indices[2]+1,indices[0]+1)+=ponderation*reg214; matrix(indices[2]+0,indices[3]+1)+=ponderation*reg217; matrix(indices[2]+0,indices[3]+0)+=ponderation*reg273; matrix(indices[2]+0,indices[2]+1)+=ponderation*reg276; matrix(indices[2]+0,indices[2]+0)+=ponderation*reg279;
    matrix(indices[2]+0,indices[1]+1)+=ponderation*reg281; matrix(indices[2]+0,indices[1]+0)+=ponderation*reg283; matrix(indices[2]+0,indices[0]+1)+=ponderation*reg285; matrix(indices[2]+0,indices[0]+0)+=ponderation*reg288; matrix(indices[1]+1,indices[3]+1)+=ponderation*reg244;
    matrix(indices[1]+1,indices[3]+0)+=ponderation*reg246; matrix(indices[1]+1,indices[2]+1)+=ponderation*reg248; matrix(indices[1]+1,indices[2]+0)+=ponderation*reg250; matrix(indices[1]+1,indices[1]+1)+=ponderation*reg253; matrix(indices[1]+1,indices[1]+0)+=ponderation*reg255;
    matrix(indices[1]+1,indices[0]+1)+=ponderation*reg258;
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
    reg1=pow((*f.m).v1[2],2); T reg4=pow((*f.m).v2[2],2); reg2=reg3+reg2; reg3=2*(*f.m).shear_modulus_23; T reg5=2*(*f.m).shear_modulus_13;
    reg1=reg0+reg1; reg0=2*(*f.m).shear_modulus_12; reg5=1.0/reg5; reg4=reg2+reg4; reg1=pow(reg1,0.5);
    reg3=1.0/reg3; reg2=reg5*reg3; T reg6=(*f.m).v1[0]/reg1; T reg7=(*f.m).v1[1]/reg1; reg4=pow(reg4,0.5);
    reg0=1.0/reg0; T reg8=reg0*reg2; reg1=(*f.m).v1[2]/reg1; T reg9=(*f.m).v2[0]/reg4; T reg10=(*f.m).v2[1]/reg4;
    T reg11=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg12=2*reg7; T reg13=2*reg6; T reg14=1.0/(*f.m).elastic_modulus_3; T reg15=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg16=2*reg1; T reg17=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; T reg18=pow(reg9,2); T reg19=1.0/(*f.m).elastic_modulus_2; T reg20=pow(reg10,2);
    T reg21=reg15*reg8; T reg22=reg9*reg13; reg4=(*f.m).v2[2]/reg4; T reg23=1.0/(*f.m).elastic_modulus_1; T reg24=reg11*reg8;
    T reg25=reg14*reg8; T reg26=reg10*reg12; T reg27=reg17*reg25; T reg28=reg11*reg24; T reg29=reg19*reg25;
    T reg30=reg11*reg21; T reg31=reg26*reg19; reg16=reg4*reg16; T reg32=reg18*reg23; T reg33=reg18*reg17;
    T reg34=reg20*reg19; T reg35=reg20*reg17; T reg36=reg26*reg17; T reg37=pow(reg4,2); T reg38=reg22*reg17;
    T reg39=reg22*reg23; reg36=reg39-reg36; reg39=reg16*reg15; T reg40=reg37*reg15; reg35=reg32-reg35;
    reg38=reg31-reg38; reg31=reg18*reg15; reg32=reg20*reg11; T reg41=reg22*reg15; T reg42=reg26*reg11;
    T reg43=pow(reg7,2); T reg44=pow(reg6,2); T reg45=reg16*reg11; reg33=reg34-reg33; reg34=reg37*reg11;
    T reg46=reg19*reg21; T reg47=reg17*reg24; reg30=reg27+reg30; reg28=reg29-reg28; reg29=reg6*reg10;
    T reg48=reg7*reg9; T reg49=pow(reg1,2); T reg50=reg6*reg9; T reg51=reg7*reg10; T reg52=reg43*reg17;
    reg40=reg35-reg40; reg39=reg36-reg39; reg35=reg43*reg19; reg34=reg33-reg34; reg45=reg38-reg45;
    reg33=reg37*reg14; reg32=reg31+reg32; reg16=reg16*reg14; reg42=reg41+reg42; reg31=reg44*reg23;
    reg36=reg44*reg17; reg38=reg23*reg28; reg41=reg17*reg30; T reg53=reg47+reg46; T reg54=reg49*reg15;
    T reg55=reg48+reg29; T reg56=reg20*reg45; T reg57=reg18*reg39; reg52=reg31-reg52; reg41=reg38-reg41;
    reg31=reg49*reg11; reg38=reg9*reg10; T reg58=reg15*reg53; reg36=reg35-reg36; reg35=reg19*reg8;
    T reg59=reg44*reg15; T reg60=reg43*reg11; reg32=reg33-reg32; reg33=reg15*reg24; reg25=reg23*reg25;
    reg42=reg16-reg42; reg16=reg15*reg21; reg8=reg17*reg8; T reg61=reg44*reg40; T reg62=reg43*reg34;
    T reg63=reg44*reg39; T reg64=reg43*reg45; T reg65=reg0*reg3; T reg66=reg18*reg40; T reg67=reg20*reg34;
    T reg68=reg11*reg2; T reg69=reg14*reg2; reg2=reg15*reg2; T reg70=reg50*reg39; T reg71=reg51*reg34;
    T reg72=reg50*reg40; T reg73=reg51*reg45; T reg74=reg6*reg4; T reg75=reg1*reg9; T reg76=reg1*reg10;
    T reg77=reg7*reg4; T reg78=2*reg9; T reg79=reg1*reg4; T reg80=reg49*reg14; reg67=reg66+reg67;
    reg66=reg79*reg32; T reg81=reg37*reg32; T reg82=reg7*reg13; T reg83=reg77-reg76; reg71=reg72+reg71;
    reg31=reg36-reg31; reg36=reg74+reg75; reg54=reg52-reg54; reg52=reg14*reg65; reg73=reg70+reg73;
    reg70=reg11*reg2; reg72=reg0*reg5; reg58=reg41-reg58; reg41=reg11*reg68; T reg84=reg17*reg69;
    reg69=reg19*reg69; T reg85=reg79*reg42; T reg86=reg15*reg65; reg60=reg59+reg60; reg59=reg37*reg42;
    reg56=reg57+reg56; reg57=reg38*reg0; T reg87=reg55*reg0; reg16=reg25-reg16; reg24=reg23*reg24;
    reg33=reg27+reg33; reg62=reg61+reg62; reg25=reg49*reg32; reg27=reg10*reg78; reg21=reg17*reg21;
    reg65=reg11*reg65; reg61=reg9*reg4; T reg88=reg15*reg35; T reg89=reg49*reg42; reg64=reg63+reg64;
    reg63=reg15*reg8; reg68=reg17*reg68; reg63=reg24+reg63; reg16=reg16/reg58; reg21=reg24+reg21;
    reg8=reg17*reg8; reg35=reg23*reg35; reg24=reg11*reg72; T reg90=reg61*reg5; T reg91=reg36*reg5;
    T reg92=reg44*reg54; T reg93=reg43*reg31; reg25=reg62+reg25; reg62=reg82*reg57; reg89=reg64+reg89;
    reg64=reg82*reg87; T reg94=reg18*reg54; T reg95=reg20*reg31; T reg96=reg55*reg87; reg14=reg14*reg72;
    T reg97=reg11*reg86; T reg98=reg11*reg65; reg60=reg80-reg60; reg80=reg17*reg52; reg52=reg19*reg52;
    reg70=reg84+reg70; reg85=reg73+reg85; reg41=reg69-reg41; reg72=reg15*reg72; reg2=reg19*reg2;
    reg69=reg27*reg87; reg59=reg56+reg59; reg28=reg28/reg58; reg56=2*reg10; reg33=reg33/reg58;
    reg73=reg1*reg13; reg66=reg71+reg66; reg71=reg55*reg57; reg30=reg30/reg58; reg84=reg27*reg57;
    T reg99=reg4*reg78; reg77=reg76+reg77; reg76=reg6*reg7; reg88=reg47+reg88; reg74=reg75-reg74;
    reg75=2*reg83; T reg100=reg10*reg4; reg81=reg67+reg81; reg67=pow(reg83,2); reg97=reg80+reg97;
    reg84=reg81+reg84; reg80=pow(reg74,2); reg81=reg100*reg3; reg98=reg52-reg98; reg52=reg77*reg3;
    T reg101=reg4*reg56; T reg102=reg19*reg14; reg14=reg17*reg14; T reg103=reg44*reg33; T reg104=reg18*reg16;
    T reg105=reg27*reg30; T reg106=reg82*reg28; T reg107=reg6*reg1; reg64=reg89+reg64; reg89=reg73*reg90;
    T reg108=reg20*reg30; T reg109=reg43*reg28; reg62=reg25+reg62; reg25=reg73*reg91; T reg110=reg99*reg90;
    T reg111=reg43*reg33; T reg112=reg20*reg16; T reg113=reg11*reg72; T reg114=reg1*reg12; T reg115=reg49*reg60;
    reg93=reg92+reg93; reg11=reg11*reg24; reg96=reg85+reg96; reg65=reg17*reg65; reg69=reg59+reg69;
    reg8=reg35-reg8; reg21=reg21/reg58; reg71=reg66+reg71; reg35=reg36*reg90; reg63=reg63/reg58;
    reg59=1-var_inter[0]; reg66=reg99*reg91; reg88=reg88/reg58; reg85=reg44*reg28; reg92=1-var_inter[1];
    T reg116=reg74*reg75; reg0=reg76*reg0; T reg117=reg82*reg33; T reg118=reg27*reg16; T reg119=reg68+reg2;
    reg70=reg17*reg70; reg95=reg94+reg95; reg94=reg36*reg91; reg53=reg53/reg58; T reg120=reg18*reg30;
    reg41=reg23*reg41; T reg121=reg50*reg54; T reg122=reg51*reg31; T reg123=reg37*reg60; reg86=reg19*reg86;
    T reg124=reg114*reg52; T reg125=reg82*reg88; T reg126=reg27*reg63; T reg127=reg7*reg1; reg123=reg95+reg123;
    reg25=reg64+reg25; reg64=reg27*reg0; reg95=reg92*elem.pos(1)[0]; T reg128=reg92*elem.pos(0)[0]; T reg129=reg92*elem.pos(1)[1];
    T reg130=reg92*elem.pos(0)[1]; T reg131=elem.pos(1)[0]*var_inter[0]; T reg132=reg59*elem.pos(0)[0]; T reg133=reg77*reg52; T reg134=elem.pos(1)[1]*var_inter[0];
    T reg135=reg59*elem.pos(0)[1]; reg5=reg107*reg5; T reg136=reg101*reg81; reg115=reg93+reg115; reg93=reg82*reg0;
    reg110=reg84+reg110; reg66=reg69+reg66; reg89=reg62+reg89; reg62=reg114*reg81; reg94=reg96+reg94;
    reg112=reg111+reg112; reg69=reg67*reg21; reg104=reg103+reg104; reg84=reg116*reg53; reg105=reg106+reg105;
    reg96=reg80*reg53; reg108=reg109+reg108; reg113=reg14+reg113; reg11=reg102-reg11; reg14=reg65+reg86;
    reg97=reg17*reg97; reg98=reg23*reg98; reg102=reg101*reg52; reg119=reg15*reg119; reg70=reg41-reg70;
    reg72=reg19*reg72; reg24=reg17*reg24; reg122=reg121+reg122; reg19=reg67*reg53; reg85=reg120+reg85;
    reg41=reg79*reg60; reg35=reg71+reg35; reg8=reg8/reg58; reg71=reg77*reg81; reg103=reg20*reg63;
    reg106=reg43*reg88; reg109=reg18*reg63; reg111=reg44*reg88; reg120=reg116*reg21; reg118=reg117+reg118;
    reg117=reg80*reg21; reg14=reg15*reg14; reg3=reg127*reg3; reg136=reg110+reg136; reg110=reg99*reg5;
    reg97=reg98-reg97; reg119=reg70-reg119; reg70=reg131+reg132; reg98=reg67*reg8; reg133=reg94+reg133;
    reg102=reg66+reg102; reg66=var_inter[1]*elem.pos(2)[1]; reg129=reg129-reg130; reg94=reg116*reg8; reg126=reg125+reg126;
    reg121=var_inter[1]*elem.pos(2)[0]; reg95=reg95-reg128; reg69=reg104+reg69; reg41=reg122+reg41; reg104=reg55*reg0;
    reg122=reg80*reg8; reg103=reg106+reg103; reg71=reg35+reg71; reg35=elem.pos(2)[0]*var_inter[0]; reg84=reg105+reg84;
    reg62=reg89+reg62; reg96=reg108+reg96; reg120=reg118+reg120; reg117=reg112+reg117; reg19=reg85+reg19;
    reg85=reg24+reg72; reg124=reg25+reg124; reg25=elem.pos(2)[1]*var_inter[0]; reg89=reg134+reg135; reg105=reg73*reg5;
    reg93=reg115+reg93; reg113=reg17*reg113; reg64=reg123+reg64; reg11=reg23*reg11; reg109=reg111+reg109;
    reg85=reg15*reg85; reg15=reg7*reg83; reg98=reg109+reg98; reg17=reg102*reg71; reg122=reg103+reg122;
    reg23=reg62*reg133; reg103=reg76*reg96; reg106=reg38*reg117; reg110=reg64+reg110; reg104=reg41+reg104;
    reg41=reg36*reg5; reg64=reg136*reg133; reg108=reg76*reg84; reg109=reg38*reg120; reg14=reg97-reg14;
    reg25=reg25-reg89; reg97=reg74*reg83; reg35=reg35-reg70; reg111=elem.pos(3)[0]*reg59; reg112=reg124*reg71;
    reg115=elem.pos(3)[1]*reg59; reg105=reg93+reg105; reg94=reg126+reg94; reg66=reg129+reg66; reg93=reg6*reg74;
    reg118=reg114*reg3; reg119=reg119/reg58; reg113=reg11-reg113; reg11=var_inter[1]*elem.pos(3)[1]; reg123=var_inter[1]*elem.pos(3)[0];
    reg125=reg76*reg19; reg126=reg38*reg69; reg121=reg95+reg121; reg95=reg101*reg3; reg112=reg23-reg112;
    reg118=reg105+reg118; reg41=reg104+reg41; reg23=reg77*reg3; reg104=reg62*reg102; reg105=reg124*reg136;
    reg95=reg110+reg95; reg110=reg6*reg83; reg129=reg7*reg74; reg93=reg15+reg93; reg106=reg103+reg106;
    reg15=reg97*reg122; reg103=reg10*reg83; T reg137=reg9*reg74; T reg138=reg18*reg69; T reg139=reg44*reg84;
    T reg140=reg44*reg19; T reg141=reg43*reg96; T reg142=reg97*reg94; reg109=reg108+reg109; reg108=reg97*reg98;
    reg17=reg64-reg17; reg126=reg125+reg126; reg14=reg14/reg58; reg121=reg121-reg123; reg64=reg55*reg119;
    reg125=reg18*reg117; reg19=reg43*reg19; reg69=reg20*reg69; reg96=reg44*reg96; T reg143=reg51*reg119;
    reg117=reg20*reg117; T reg144=reg50*reg119; reg66=reg66-reg11; reg111=reg35+reg111; reg25=reg115+reg25;
    reg85=reg113-reg85; reg84=reg43*reg84; reg35=reg18*reg120; reg120=reg20*reg120; reg108=reg126+reg108;
    reg113=reg55*reg144; reg115=reg129*reg14; reg126=reg93*reg14; T reg145=reg67*reg122; T reg146=reg80*reg94;
    reg120=reg84+reg120; reg138=reg140+reg138; reg84=reg67*reg98; reg117=reg141+reg117; reg125=reg96+reg125;
    reg122=reg80*reg122; reg105=reg104-reg105; reg69=reg19+reg69; reg98=reg80*reg98; reg94=reg67*reg94;
    reg35=reg139+reg35; reg19=reg121*reg25; reg96=reg66*reg111; reg23=reg41+reg23; reg41=reg9*reg83;
    reg104=reg10*reg74; reg137=reg103+reg137; reg103=reg55*reg64; reg142=reg109+reg142; reg109=reg17*reg118;
    reg139=reg95*reg112; reg58=reg85/reg58; reg85=reg55*reg143; reg15=reg106+reg15; reg106=reg110*reg14;
    reg140=reg22*reg143; reg145=reg125+reg145; reg94=reg35+reg94; reg35=reg22*reg64; reg125=reg22*reg144;
    reg84=reg138+reg84; reg138=reg137*reg58; reg141=reg104*reg58; T reg147=reg41*reg58; reg12=reg74*reg12;
    reg13=reg83*reg13; reg96=reg19-reg96; reg19=reg62*reg95; T reg148=reg118*reg136; T reg149=reg62*reg23;
    T reg150=reg118*reg71; T reg151=reg136*reg23; T reg152=reg95*reg71; T reg153=reg23*reg105; reg139=reg109-reg139;
    reg122=reg117+reg122; reg143=reg26*reg143; reg109=reg93*reg115; reg85=reg15+reg85; reg15=reg93*reg106;
    reg113=reg108+reg113; reg103=reg142+reg103; reg108=reg93*reg126; reg144=reg26*reg144; reg98=reg69+reg98;
    reg146=reg120+reg146; reg64=reg26*reg64; reg149=reg150-reg149; reg69=reg118*reg102; reg15=reg113+reg15;
    reg113=reg137*reg147; reg117=reg12*reg126; reg64=reg146+reg64; reg109=reg85+reg109; reg85=reg137*reg141;
    reg120=reg124*reg23; reg56=reg74*reg56; reg153=reg139+reg153; reg78=reg83*reg78; reg108=reg103+reg108;
    reg103=reg137*reg138; reg151=reg152-reg151; reg139=reg95*reg133; reg142=reg102*reg23; reg146=reg118*reg133;
    reg121=reg121/reg96; reg25=reg25/reg96; reg150=reg13*reg115; reg140=reg145+reg140; reg145=reg12*reg106;
    reg144=reg98+reg144; reg143=reg122+reg143; reg66=reg66/reg96; reg35=reg94+reg35; reg19=reg148-reg19;
    reg106=reg13*reg106; reg125=reg84+reg125; reg126=reg13*reg126; reg115=reg12*reg115; reg111=reg111/reg96;
    reg84=reg124*reg95; reg151=reg151/reg153; reg94=reg92*reg25; reg142=reg139-reg142; reg98=var_inter[0]*reg121;
    reg122=var_inter[1]*reg25; reg139=var_inter[1]*reg111; reg148=1-(*f.m).resolution; reg152=var_inter[0]*reg66; reg120=reg146-reg120;
    reg146=reg59*reg121; reg149=reg149/reg153; reg19=reg19/reg153; reg84=reg69-reg84; reg126=reg35+reg126;
    reg35=reg78*reg138; reg69=reg78*reg141; reg150=reg140+reg150; reg145=reg144+reg145; reg140=reg56*reg147;
    reg147=reg78*reg147; reg106=reg125+reg106; reg115=reg143+reg115; reg141=reg56*reg141; reg117=reg64+reg117;
    reg64=reg92*reg111; reg125=reg59*reg66; reg103=reg108+reg103; reg85=reg109+reg85; reg113=reg15+reg113;
    reg138=reg56*reg138; reg105=reg105/reg153; reg35=reg126+reg35; reg69=reg150+reg69; reg84=reg84/reg153;
    reg15=reg98+reg64; reg142=reg142/reg153; reg108=reg139+reg146; reg17=reg17/reg153; reg109=reg122-reg152;
    reg112=reg112/reg153; reg126=reg98-reg139; reg140=reg145+reg140; reg143=(*f.m).resolution*reg19; reg144=reg122+reg125;
    reg145=reg64-reg146; reg138=reg117+reg138; reg153=reg120/reg153; reg117=reg103*reg148; reg120=reg85*reg148;
    reg150=reg113*reg148; T reg154=(*f.m).resolution*reg151; T reg155=(*f.m).resolution*reg149; reg141=reg115+reg141; reg115=reg125-reg94;
    T reg156=reg94+reg152; reg147=reg106+reg147; reg106=reg147*reg148; T reg157=(*f.m).resolution*reg84; T reg158=(*f.m).resolution*reg153;
    T reg159=(*f.m).resolution*reg142; T reg160=reg69*reg148; T reg161=reg35*reg148; T reg162=reg140*reg148; T reg163=(*f.m).resolution*reg105;
    T reg164=(*f.m).resolution*reg17; T reg165=(*f.m).resolution*reg112; T reg166=reg138*reg148; T reg167=reg141*reg148; T reg168=0.5*reg156;
    T reg169=0.5*reg109; T reg170=0.5*reg126; T reg171=0.5*reg144; T reg172=0.5*reg108; T reg173=0.5*reg115;
    reg150=reg154+reg150; reg155=reg120-reg155; reg117=reg143+reg117; reg120=0.5*reg15; reg143=0.5*reg145;
    reg154=reg15*reg155; T reg174=reg108*reg155; T reg175=reg120*reg117; reg159=reg162-reg159; reg162=reg168*reg117;
    T reg176=reg171*reg117; reg167=reg158+reg167; reg158=reg170*reg117; T reg177=reg109*reg150; reg157=reg166-reg157;
    reg166=reg156*reg150; T reg178=reg144*reg150; T reg179=reg172*reg117; reg161=reg163+reg161; reg165=reg160-reg165;
    reg106=reg164+reg106; reg160=reg145*reg155; reg163=reg126*reg155; reg164=reg173*reg117; T reg180=reg169*reg117;
    T reg181=reg115*reg150; T reg182=reg143*reg117; T reg183=reg172*reg161; T reg184=reg171*reg157; T reg185=reg109*reg106;
    T reg186=reg170*reg161; T reg187=reg143*reg161; T reg188=reg144*reg106; reg182=reg181+reg182; reg180=reg163+reg180;
    reg158=reg177+reg158; reg163=reg108*reg167; reg177=reg145*reg165; reg162=reg162-reg154; reg181=reg169*reg161;
    T reg189=reg173*reg161; T reg190=reg126*reg165; T reg191=reg15*reg165; reg164=reg160+reg164; reg160=reg156*reg106;
    T reg192=reg168*reg161; T reg193=reg156*reg159; T reg194=reg120*reg157; T reg195=reg173*reg157; T reg196=reg15*reg167;
    T reg197=reg168*reg157; T reg198=reg145*reg167; reg174=reg174-reg176; T reg199=reg109*reg159; T reg200=reg170*reg157;
    reg166=reg166-reg175; T reg201=reg171*reg161; T reg202=reg126*reg167; T reg203=reg169*reg157; T reg204=reg108*reg165;
    T reg205=reg144*reg159; T reg206=reg172*reg157; T reg207=reg115*reg106; T reg208=reg120*reg161; reg179=reg179-reg178;
    reg189=reg177+reg189; reg187=reg207+reg187; reg182=2*reg182; reg163=reg163-reg184; reg164=2*reg164;
    reg200=reg199+reg200; reg162=2*reg162; reg192=reg192-reg191; reg203=reg202+reg203; reg197=reg197-reg196;
    reg183=reg183-reg188; reg206=reg206-reg205; reg166=2*reg166; reg179=2*reg179; reg160=reg160-reg208;
    reg186=reg185+reg186; reg193=reg193-reg194; reg174=2*reg174; reg195=reg198+reg195; reg180=2*reg180;
    reg181=reg190+reg181; reg204=reg204-reg201; reg158=2*reg158; reg177=reg109*reg181; reg185=reg115*reg183;
    reg190=reg170*reg180; reg198=reg115*reg204; reg199=reg115*reg187; reg202=reg143*reg179; reg207=reg173*reg174;
    T reg209=reg15*reg206; T reg210=reg169*reg174; T reg211=reg143*reg162; T reg212=reg143*reg182; T reg213=reg143*reg174;
    T reg214=reg120*reg180; T reg215=reg156*reg183; T reg216=reg156*reg192; T reg217=reg168*reg174; T reg218=reg15*reg163;
    T reg219=reg120*reg162; T reg220=reg120*reg166; T reg221=reg156*reg160; T reg222=reg168*reg179; T reg223=reg172*reg179;
    T reg224=reg126*reg206; T reg225=reg145*reg163; T reg226=reg156*reg186; T reg227=reg115*reg160; T reg228=reg144*reg183;
    T reg229=reg170*reg158; T reg230=reg120*reg158; T reg231=reg143*reg164; T reg232=reg143*reg166; T reg233=reg109*reg186;
    T reg234=reg115*reg189; T reg235=reg115*reg192; T reg236=reg156*reg181; T reg237=reg168*reg162; T reg238=reg15*reg197;
    T reg239=reg173*reg180; T reg240=reg145*reg203; T reg241=reg169*reg180; T reg242=reg108*reg163; T reg243=reg173*reg179;
    T reg244=reg145*reg206; T reg245=reg15*reg203; T reg246=reg170*reg174; T reg247=reg171*reg174; T reg248=reg168*reg180;
    T reg249=reg115*reg186; T reg250=reg143*reg158; T reg251=reg168*reg158; T reg252=reg109*reg204; T reg253=reg15*reg200;
    T reg254=reg115*reg181; T reg255=reg143*reg180; T reg256=reg126*reg203; T reg257=reg120*reg179; T reg258=reg172*reg174;
    T reg259=reg173*reg164; T reg260=reg145*reg195; T reg261=reg169*reg179; T reg262=reg144*reg204; T reg263=reg173*reg166;
    T reg264=reg145*reg193; reg163=reg126*reg163; reg204=reg156*reg204; reg174=reg120*reg174; T reg265=reg173*reg162;
    T reg266=reg145*reg197; T reg267=reg170*reg179; reg183=reg109*reg183; T reg268=reg173*reg158; T reg269=reg145*reg200;
    reg174=reg204-reg174; reg214=reg236-reg214; reg238=reg237-reg238; reg230=reg226-reg230; reg209=reg222-reg209;
    reg219=reg216-reg219; reg257=reg215-reg257; reg245=reg248-reg245; reg253=reg251-reg253; reg261=reg224+reg261;
    reg210=reg163+reg210; reg255=reg254+reg255; reg252=reg246+reg252; reg250=reg249+reg250; reg247=reg242-reg247;
    reg244=reg243+reg244; reg240=reg239+reg240; reg183=reg267+reg183; reg269=reg268+reg269; reg266=reg265+reg266;
    reg262=reg258-reg262; reg241=reg256+reg241; reg264=reg263+reg264; reg260=reg259+reg260; reg177=reg190+reg177;
    reg220=reg221-reg220; reg218=reg217-reg218; reg225=reg207+reg225; reg231=reg234+reg231; reg232=reg227+reg232;
    reg228=reg223-reg228; reg233=reg229+reg233; reg212=reg199+reg212; reg211=reg235+reg211; reg202=reg185+reg202;
    reg213=reg198+reg213; reg233=reg96*reg233; reg245=reg96*reg245; reg183=reg96*reg183; reg177=reg96*reg177;
    reg218=reg96*reg218; reg209=reg96*reg209; reg241=reg96*reg241; reg252=reg96*reg252; reg255=reg96*reg255;
    reg247=reg96*reg247; reg250=reg96*reg250; reg244=reg96*reg244; reg240=reg96*reg240; reg269=reg96*reg269;
    reg262=reg96*reg262; reg266=reg96*reg266; reg264=reg96*reg264; reg260=reg96*reg260; reg213=reg96*reg213;
    reg202=reg96*reg202; reg211=reg96*reg211; reg212=reg96*reg212; reg228=reg96*reg228; reg232=reg96*reg232;
    reg253=reg96*reg253; reg261=reg96*reg261; reg238=reg96*reg238; reg174=reg96*reg174; reg257=reg96*reg257;
    reg214=reg96*reg214; reg230=reg96*reg230; reg210=reg96*reg210; reg219=reg96*reg219; reg220=reg96*reg220;
    reg225=reg96*reg225; reg231=reg96*reg231; matrix(indices[2]+1,indices[3]+1)+=ponderation*reg210; matrix(indices[3]+0,indices[3]+0)+=ponderation*reg228; matrix(indices[3]+0,indices[3]+1)+=ponderation*reg262;
    matrix(indices[2]+1,indices[3]+0)+=ponderation*reg261; matrix(indices[2]+1,indices[2]+1)+=ponderation*reg241; matrix(indices[3]+1,indices[3]+1)+=ponderation*reg247; matrix(indices[0]+0,indices[2]+1)+=ponderation*reg255; matrix(indices[0]+0,indices[2]+0)+=ponderation*reg250;
    matrix(indices[0]+0,indices[1]+1)+=ponderation*reg211; matrix(indices[0]+1,indices[3]+0)+=ponderation*reg244; matrix(indices[0]+1,indices[2]+1)+=ponderation*reg240; matrix(indices[0]+1,indices[2]+0)+=ponderation*reg269; matrix(indices[0]+1,indices[1]+1)+=ponderation*reg266;
    matrix(indices[0]+1,indices[1]+0)+=ponderation*reg264; matrix(indices[0]+1,indices[0]+1)+=ponderation*reg260; matrix(indices[0]+0,indices[3]+1)+=ponderation*reg213; matrix(indices[0]+0,indices[3]+0)+=ponderation*reg202; matrix(indices[0]+0,indices[0]+0)+=ponderation*reg212;
    matrix(indices[0]+0,indices[1]+0)+=ponderation*reg232; matrix(indices[0]+0,indices[0]+1)+=ponderation*reg231; matrix(indices[0]+1,indices[3]+1)+=ponderation*reg225; matrix(indices[1]+0,indices[1]+0)+=ponderation*reg220; matrix(indices[1]+0,indices[1]+1)+=ponderation*reg219;
    matrix(indices[1]+0,indices[2]+0)+=ponderation*reg230; matrix(indices[1]+0,indices[2]+1)+=ponderation*reg214; matrix(indices[1]+0,indices[3]+0)+=ponderation*reg257; matrix(indices[1]+0,indices[3]+1)+=ponderation*reg174; matrix(indices[1]+1,indices[1]+1)+=ponderation*reg238;
    matrix(indices[1]+1,indices[2]+0)+=ponderation*reg253; matrix(indices[1]+1,indices[2]+1)+=ponderation*reg245; matrix(indices[1]+1,indices[3]+0)+=ponderation*reg209; matrix(indices[1]+1,indices[3]+1)+=ponderation*reg218; matrix(indices[2]+0,indices[2]+0)+=ponderation*reg233;
    matrix(indices[2]+0,indices[2]+1)+=ponderation*reg177; matrix(indices[2]+0,indices[3]+0)+=ponderation*reg183; matrix(indices[2]+0,indices[3]+1)+=ponderation*reg252;
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
    T reg0=2*(*f.m).shear_modulus_23; T reg1=2*(*f.m).shear_modulus_13; reg0=1.0/reg0; reg1=1.0/reg1; T reg2=2*(*f.m).shear_modulus_12;
    reg2=1.0/reg2; T reg3=reg1*reg0; T reg4=pow((*f.m).v1[0],2); T reg5=pow((*f.m).v1[1],2); T reg6=reg2*reg3;
    T reg7=pow((*f.m).v2[0],2); T reg8=pow((*f.m).v2[1],2); T reg9=(*f.m).poisson_ratio_23/(*f.m).elastic_modulus_2; T reg10=1.0/(*f.m).elastic_modulus_3; T reg11=(*f.m).poisson_ratio_13/(*f.m).elastic_modulus_1;
    T reg12=pow((*f.m).v2[2],2); T reg13=reg10*reg6; T reg14=1.0/(*f.m).elastic_modulus_2; T reg15=(*f.m).poisson_ratio_12/(*f.m).elastic_modulus_1; reg8=reg7+reg8;
    reg7=reg9*reg6; reg5=reg4+reg5; reg4=pow((*f.m).v1[2],2); T reg16=reg11*reg6; T reg17=reg14*reg13;
    T reg18=reg9*reg7; T reg19=reg15*reg13; reg4=reg5+reg4; reg5=reg9*reg16; reg12=reg8+reg12;
    reg5=reg19+reg5; reg8=reg15*reg7; reg18=reg17-reg18; reg17=reg14*reg16; reg4=pow(reg4,0.5);
    T reg20=1.0/(*f.m).elastic_modulus_1; reg12=pow(reg12,0.5); T reg21=reg20*reg18; T reg22=(*f.m).v2[1]/reg12; T reg23=(*f.m).v2[2]/reg12;
    T reg24=reg15*reg5; T reg25=(*f.m).v1[1]/reg4; T reg26=(*f.m).v1[2]/reg4; T reg27=reg8+reg17; reg12=(*f.m).v2[0]/reg12;
    reg4=(*f.m).v1[0]/reg4; reg24=reg21-reg24; reg21=reg10*reg3; T reg28=reg11*reg3; reg3=reg9*reg3;
    T reg29=reg2*reg0; T reg30=reg15*reg6; T reg31=reg11*reg16; reg13=reg20*reg13; T reg32=reg25*reg23;
    T reg33=reg26*reg22; T reg34=reg11*reg7; T reg35=reg11*reg27; reg6=reg14*reg6; T reg36=reg11*reg29;
    reg35=reg24-reg35; reg24=reg2*reg1; T reg37=2*reg12; T reg38=reg4*reg23; T reg39=reg14*reg21;
    reg21=reg15*reg21; reg16=reg15*reg16; T reg40=reg9*reg29; T reg41=reg9*reg3; T reg42=reg9*reg28;
    T reg43=reg11*reg30; T reg44=2*reg4; reg7=reg20*reg7; reg31=reg13-reg31; reg13=reg11*reg6;
    reg29=reg10*reg29; T reg45=reg32-reg33; T reg46=reg26*reg12; reg34=reg19+reg34; reg19=reg9*reg36;
    T reg47=2*reg45; T reg48=reg9*reg40; T reg49=reg10*reg24; T reg50=reg22*reg37; T reg51=pow(reg22,2);
    T reg52=reg4*reg22; T reg53=pow(reg12,2); reg18=reg18/reg35; T reg54=reg25*reg12; T reg55=reg25*reg44;
    T reg56=pow(reg25,2); T reg57=pow(reg4,2); T reg58=reg46-reg38; reg30=reg15*reg30; reg16=reg7+reg16;
    T reg59=reg15*reg29; reg43=reg7+reg43; reg41=reg39-reg41; reg42=reg21+reg42; reg6=reg20*reg6;
    reg31=reg31/reg35; reg3=reg15*reg3; reg13=reg8+reg13; reg7=reg9*reg24; reg28=reg14*reg28;
    reg5=reg5/reg35; reg29=reg14*reg29; reg24=reg11*reg24; reg34=reg34/reg35; reg21=reg50*reg5;
    reg41=reg20*reg41; reg39=reg55*reg18; reg36=reg14*reg36; T reg60=reg58*reg47; reg48=reg29-reg48;
    reg42=reg15*reg42; reg29=reg3+reg28; reg19=reg59+reg19; reg59=reg51*reg5; T reg61=reg56*reg18;
    T reg62=reg14*reg49; reg49=reg15*reg49; T reg63=reg9*reg7; T reg64=reg9*reg24; reg13=reg13/reg35;
    T reg65=reg57*reg18; T reg66=reg50*reg31; T reg67=pow(reg26,2); T reg68=reg55*reg34; T reg69=reg52-reg54;
    reg43=reg43/reg35; reg16=reg16/reg35; reg30=reg6-reg30; reg6=pow(reg23,2); T reg70=reg53*reg5;
    T reg71=reg51*reg31; reg27=reg27/reg35; T reg72=reg56*reg34; T reg73=reg53*reg31; T reg74=pow(reg45,2);
    T reg75=pow(reg58,2); reg40=reg15*reg40; T reg76=reg57*reg34; reg63=reg62-reg63; reg62=reg40+reg36;
    reg19=reg15*reg19; reg29=reg11*reg29; reg30=reg30/reg35; reg42=reg41-reg42; reg48=reg20*reg48;
    reg24=reg14*reg24; reg7=reg15*reg7; reg41=reg74*reg27; reg65=reg70+reg65; reg70=pow(reg69,2);
    T reg77=reg60*reg27; reg73=reg76+reg73; reg21=reg39+reg21; reg39=reg74*reg16; reg71=reg72+reg71;
    reg72=reg75*reg16; reg76=reg67*reg34; T reg78=reg6*reg31; T reg79=reg6*reg5; T reg80=reg67*reg18;
    reg66=reg68+reg66; reg68=reg60*reg16; T reg81=reg75*reg27; T reg82=reg57*reg13; T reg83=reg53*reg43;
    reg59=reg61+reg59; reg61=reg56*reg13; T reg84=reg51*reg43; reg64=reg49+reg64; reg49=reg50*reg43;
    T reg85=reg55*reg13; reg29=reg42-reg29; reg42=2*reg22; T reg86=reg23*reg37; T reg87=2*reg25;
    T reg88=reg26*reg44; T reg89=reg60*reg30; reg77=reg21+reg77; reg49=reg85+reg49; reg39=reg73+reg39;
    reg21=reg6*reg43; reg72=reg71+reg72; reg71=reg67*reg13; reg73=reg75*reg30; reg78=reg76+reg78;
    reg84=reg61+reg84; reg61=reg70*reg16; reg68=reg66+reg68; reg66=reg74*reg30; reg83=reg82+reg83;
    reg62=reg11*reg62; reg76=reg12*reg22; reg82=reg4*reg25; reg81=reg59+reg81; reg19=reg48-reg19;
    reg79=reg80+reg79; reg41=reg65+reg41; reg48=reg70*reg27; reg59=reg7+reg24; reg63=reg20*reg63;
    reg64=reg15*reg64; reg65=reg58*reg45; reg29=reg29/reg35; reg80=reg53*reg68; reg85=reg57*reg77;
    reg66=reg83+reg66; reg83=reg53*reg72; T reg90=reg57*reg81; T reg91=2*reg26; reg62=reg19-reg62;
    reg19=reg53*reg39; T reg92=reg51*reg68; T reg93=reg22*reg87; reg73=reg84+reg73; reg84=reg56*reg77;
    reg21=reg71+reg21; reg71=reg70*reg30; T reg94=reg12*reg44; reg59=reg11*reg59; T reg95=reg57*reg41;
    reg64=reg63-reg64; reg89=reg49+reg89; reg49=reg86*reg5; reg63=reg88*reg18; T reg96=reg4*reg12;
    T reg97=reg51*reg39; T reg98=2*reg58; reg47=reg69*reg47; T reg99=reg25*reg22; T reg100=reg56*reg81;
    T reg101=reg51*reg72; reg52=reg54+reg52; reg54=reg26*reg87; T reg102=reg56*reg41; reg48=reg79+reg48;
    reg79=reg25*reg45; reg81=reg82*reg81; T reg103=reg4*reg58; reg61=reg78+reg61; reg77=reg82*reg77;
    reg68=reg76*reg68; reg78=reg23*reg42; reg72=reg76*reg72; T reg104=reg88*reg34; T reg105=reg86*reg31;
    T reg106=reg74*reg89; reg80=reg85+reg80; reg19=reg95+reg19; reg85=reg74*reg66; reg83=reg90+reg83;
    reg90=reg74*reg73; reg95=reg57*reg48; T reg107=reg53*reg61; reg105=reg104+reg105; reg104=reg47*reg16;
    reg34=reg54*reg34; reg31=reg78*reg31; reg62=reg62/reg35; reg59=reg64-reg59; reg71=reg21+reg71;
    reg5=reg78*reg5; reg18=reg54*reg18; reg21=reg47*reg27; reg49=reg63+reg49; reg63=reg88*reg13;
    reg64=reg86*reg43; T reg108=reg65*reg73; T reg109=reg96*reg29; T reg110=reg99*reg29; T reg111=reg52*reg29;
    reg103=reg79+reg103; reg79=reg94*reg20; T reg112=reg93*reg15; reg41=reg82*reg41; reg39=reg76*reg39;
    T reg113=reg25*reg58; T reg114=reg4*reg45; T reg115=reg75*reg89; reg92=reg84+reg92; reg84=reg22*reg45;
    T reg116=reg12*reg58; T reg117=reg26*reg23; T reg118=reg53*reg15; reg98=reg69*reg98; T reg119=reg51*reg15;
    reg97=reg102+reg97; reg72=reg81+reg72; reg101=reg100+reg101; reg81=reg51*reg61; reg100=reg56*reg48;
    reg73=reg75*reg73; reg102=reg93*reg14; reg68=reg77+reg68; reg89=reg65*reg89; reg91=reg23*reg91;
    reg77=reg94*reg15; T reg120=reg53*reg20; T reg121=reg51*reg14; T reg122=reg75*reg66; reg73=reg101+reg73;
    reg119=reg120-reg119; reg101=reg94*reg11; reg120=reg93*reg110; T reg123=reg6*reg11; T reg124=reg52*reg110;
    T reg125=reg75*reg71; reg81=reg100+reg81; reg61=reg76*reg61; reg77=reg102-reg77; reg100=reg91*reg9;
    reg48=reg82*reg48; reg118=reg121-reg118; reg102=reg53*reg11; reg16=reg98*reg16; reg31=reg34+reg31;
    reg34=reg6*reg9; reg104=reg105+reg104; reg105=reg51*reg9; reg85=reg19+reg85; reg112=reg79-reg112;
    reg19=reg91*reg11; reg79=reg94*reg109; reg121=reg22*reg58; reg64=reg63+reg64; reg116=reg84+reg116;
    reg63=reg47*reg30; reg84=reg113*reg62; reg107=reg95+reg107; reg35=reg59/reg35; reg59=reg74*reg71;
    reg95=reg52*reg111; reg89=reg68+reg89; reg44=reg45*reg44; reg13=reg54*reg13; reg87=reg58*reg87;
    reg43=reg78*reg43; reg68=reg114*reg62; T reg126=reg94*reg111; reg106=reg80+reg106; reg108=reg72+reg108;
    reg72=reg93*reg9; reg80=reg117*reg29; reg115=reg92+reg115; reg111=reg93*reg111; reg27=reg98*reg27;
    reg5=reg18+reg5; reg18=reg103*reg62; reg92=reg93*reg109; reg21=reg49+reg21; reg122=reg97+reg122;
    reg49=reg26*reg69; reg90=reg83+reg90; reg110=reg94*reg110; reg39=reg41+reg39; reg66=reg65*reg66;
    reg41=reg12*reg45; reg83=reg41*reg35; reg97=reg49*reg62; T reg127=reg121*reg35; reg61=reg48+reg61;
    reg71=reg65*reg71; reg72=reg101+reg72; reg48=reg51*reg104; reg101=reg56*reg21; reg124=reg108+reg124;
    reg108=reg103*reg84; reg126=reg106+reg126; reg106=reg94*reg80; reg59=reg107+reg59; reg107=reg44*reg18;
    T reg128=reg44*reg84; reg110=reg90+reg110; reg90=reg87*reg18; reg111=reg115+reg111; reg115=reg57*reg21;
    T reg129=reg53*reg104; T reg130=reg44*reg68; reg79=reg85+reg79; reg85=(*f.m).alpha_1*reg57; reg20=reg57*reg20;
    T reg131=reg116*reg35; reg84=reg87*reg84; reg19=reg112-reg19; reg123=reg119-reg123; reg120=reg73+reg120;
    reg73=reg56*reg15; reg112=(*f.m).alpha_2*reg53; reg42=reg58*reg42; reg119=reg87*reg68; reg27=reg5+reg27;
    reg92=reg122+reg92; reg46=reg38+reg46; reg5=reg26*reg45; reg38=reg4*reg69; reg122=reg23*reg69;
    reg66=reg39+reg66; reg37=reg45*reg37; reg95=reg89+reg95; reg109=reg52*reg109; reg18=reg103*reg18;
    reg91=reg91*reg10; reg30=reg98*reg30; reg43=reg13+reg43; reg63=reg64+reg63; reg15=reg57*reg15;
    reg13=(*f.m).alpha_1*reg56; reg39=(*f.m).alpha_2*reg51; reg105=reg102+reg105; reg64=reg6*reg10; reg89=reg93*reg80;
    reg125=reg81+reg125; reg100=reg77-reg100; reg34=reg118-reg34; reg16=reg31+reg16; reg14=reg56*reg14;
    reg31=reg57*reg27; reg109=reg66+reg109; reg66=reg53*reg16; reg48=reg101+reg48; reg68=reg103*reg68;
    reg77=reg75*reg63; reg81=reg42*reg131; reg90=reg111+reg90; reg107=reg126+reg107; reg101=reg37*reg131;
    reg119=reg92+reg119; reg92=reg42*reg83; reg84=reg120+reg84; reg102=reg42*reg127; reg129=reg115+reg129;
    reg89=reg125+reg89; reg111=reg87*reg97; reg115=reg51*reg16; reg118=reg74*reg63; reg120=reg56*reg27;
    reg125=reg26*reg58; reg38=reg5+reg38; reg32=reg33+reg32; reg131=reg116*reg131; reg5=(*f.m).alpha_1*reg67;
    reg33=(*f.m).alpha_2*reg6; reg21=reg82*reg21; reg104=reg76*reg104; reg126=reg56*reg100; T reg132=reg57*reg19;
    T reg133=reg56*reg34; T reg134=reg57*reg123; T reg135=reg51*reg100; T reg136=reg53*reg19; T reg137=reg67*reg11;
    reg73=reg20-reg73; reg20=reg67*reg9; reg15=reg14-reg15; reg11=reg57*reg11; reg9=reg56*reg9;
    reg72=reg91-reg72; reg14=reg99*reg34; reg91=reg99*reg100; reg105=reg64-reg105; reg64=reg96*reg19;
    T reg138=reg96*reg123; reg108=reg124+reg108; reg124=reg44*reg97; reg106=reg59+reg106; reg59=reg116*reg127;
    reg127=reg37*reg127; reg128=reg110+reg128; reg110=reg37*reg83; reg130=reg79+reg130; reg79=reg122*reg35;
    reg71=reg61+reg71; reg61=reg46*reg29; reg30=reg43+reg30; reg112=reg85+reg112; reg43=(*f.m).alpha_3*reg74;
    reg85=reg25*reg69; reg18=reg95+reg18; reg95=reg23*reg45; T reg139=reg12*reg69; T reg140=(*f.m).alpha_3*reg75;
    reg39=reg13+reg39; reg13=reg53*reg123; T reg141=reg51*reg34; reg80=reg52*reg80; T reg142=reg52*reg2;
    reg126=reg132+reg126; reg132=reg76*reg2; T reg143=reg93*reg61; T reg144=reg6*reg105; reg141=reg13+reg141;
    reg77=reg48+reg77; reg13=reg67*reg105; reg48=reg67*reg72; reg133=reg134+reg133; reg59=reg108+reg59;
    reg83=reg116*reg83; reg68=reg109+reg68; reg80=reg71+reg80; reg97=reg103*reg97; reg71=reg6*reg72;
    reg131=reg18+reg131; reg104=reg21+reg104; reg63=reg65*reg63; reg135=reg136+reg135; reg18=reg75*reg30;
    reg27=reg82*reg27; reg16=reg76*reg16; reg137=reg73-reg137; reg20=reg15-reg20; reg115=reg120+reg115;
    reg10=reg67*reg10; reg9=reg11+reg9; reg66=reg31+reg66; reg91=reg64+reg91; reg11=reg117*reg72;
    reg15=reg94*reg61; reg118=reg129+reg118; reg101=reg107+reg101; reg21=reg37*reg79; reg124=reg106+reg124;
    reg127=reg128+reg127; reg110=reg130+reg110; reg31=reg38*reg62; reg29=reg32*reg29; reg43=reg112+reg43;
    reg12=reg12*reg23; reg4=reg4*reg26; reg64=reg22*reg69; reg73=reg23*reg58; reg139=reg95+reg139;
    reg85=reg125+reg85; reg140=reg39+reg140; reg33=reg5+reg33; reg5=(*f.m).alpha_3*reg70; reg39=(*f.m).alpha_1*reg82;
    reg76=(*f.m).alpha_2*reg76; reg81=reg90+reg81; reg90=reg42*reg79; reg111=reg89+reg111; reg102=reg84+reg102;
    reg92=reg119+reg92; reg84=reg117*reg105; reg89=reg74*reg30; reg14=reg138+reg14; reg95=(*f.m).alpha_1*reg4;
    reg30=reg65*reg30; reg16=reg27+reg16; reg65=(*f.m).alpha_3*reg65; reg76=reg39+reg76; reg27=reg92*reg43;
    reg39=reg110*reg43; reg61=reg52*reg61; reg106=reg127*reg140; reg107=reg55*reg142; reg108=reg50*reg132;
    reg109=reg53*reg137; reg112=reg51*reg20; reg144=reg141+reg144; reg63=reg104+reg63; reg104=reg52*reg142;
    reg11=reg91+reg11; reg9=reg10-reg9; reg10=reg12*reg1; reg91=reg46*reg1; reg119=reg57*reg137;
    reg120=reg101*reg59; reg125=reg56*reg20; reg128=reg102*reg140; reg84=reg14+reg84; reg14=reg52*reg132;
    reg129=reg81*reg59; reg130=reg127*reg131; reg134=reg102*reg131; reg13=reg133+reg13; reg133=reg55*reg132;
    reg48=reg126+reg48; reg126=reg50*reg142; reg12=(*f.m).alpha_2*reg12; reg83=reg68+reg83; reg68=reg44*reg31;
    reg15=reg118+reg15; reg89=reg66+reg89; reg94=reg94*reg29; reg71=reg135+reg71; reg21=reg124+reg21;
    reg93=reg93*reg29; reg66=reg139*reg35; reg62=reg85*reg62; reg18=reg115+reg18; reg23=reg22*reg23;
    reg97=reg80+reg97; reg5=reg33+reg5; reg22=reg87*reg31; reg64=reg73+reg64; reg90=reg111+reg90;
    reg79=reg116*reg79; reg45=reg69*reg45; reg26=reg25*reg26; reg143=reg77+reg143; reg25=reg46*reg10;
    reg33=reg32*reg0; reg73=reg23*reg0; reg44=reg44*reg62; reg90=reg90*reg5; reg2=reg82*reg2;
    reg128=reg27+reg128; reg94=reg89+reg94; reg27=reg46*reg91; reg104=reg11+reg104; reg21=reg21*reg5;
    reg106=reg39+reg106; reg11=reg37*reg66; reg22=reg143+reg22; reg39=reg42*reg66; reg112=reg109+reg112;
    reg77=reg88*reg91; reg107=reg48+reg107; reg48=reg59*reg140; reg80=reg83*reg43; reg108=reg144+reg108;
    reg82=reg88*reg10; reg133=reg13+reg133; reg6=reg6*reg9; reg93=reg18+reg93; reg14=reg84+reg14;
    reg87=reg87*reg62; reg13=reg86*reg10; reg67=reg67*reg9; reg125=reg119+reg125; reg58=reg69*reg58;
    reg35=reg64*reg35; reg129=reg134-reg129; reg23=(*f.m).alpha_2*reg23; reg18=(*f.m).alpha_1*reg26; reg126=reg71+reg126;
    reg45=(*f.m).alpha_3*reg45; reg120=reg130-reg120; reg69=reg127*reg81; reg71=reg86*reg91; reg12=reg95+reg12;
    reg84=reg101*reg102; reg29=reg52*reg29; reg68=reg15+reg68; reg61=reg63+reg61; reg31=reg103*reg31;
    reg65=reg76+reg65; reg15=reg99*reg20; reg63=reg96*reg137; reg30=reg16+reg30; reg79=reg97+reg79;
    reg16=reg32*reg73; reg82=reg133+reg82; reg25=reg14+reg25; reg14=reg54*reg73; reg23=reg18+reg23;
    reg29=reg30+reg29; reg48=reg80+reg48; reg31=reg61+reg31; reg79=reg79*reg5; reg45=reg12+reg45;
    reg39=reg22+reg39; reg66=reg116*reg66; reg62=reg103*reg62; reg12=reg54*reg33; reg6=reg112+reg6;
    reg18=reg50*reg2; reg77=reg107+reg77; reg21=reg106+reg21; reg22=reg101*reg65; reg117=reg117*reg9;
    reg15=reg63+reg15; reg11=reg68+reg11; reg30=reg78*reg33; reg90=reg128+reg90; reg1=reg4*reg1;
    reg84=reg69-reg84; reg44=reg94+reg44; reg71=reg126+reg71; reg4=reg81*reg65; reg37=reg37*reg35;
    reg67=reg125+reg67; reg61=reg55*reg2; reg63=reg78*reg73; reg13=reg108+reg13; reg68=reg32*reg33;
    reg27=reg104+reg27; reg42=reg42*reg35; reg87=reg93+reg87; reg69=reg92*reg120; reg76=reg110*reg129;
    reg58=(*f.m).alpha_3*reg58; reg16=reg25+reg16; reg25=reg52*reg2; reg117=reg15+reg117; reg30=reg71+reg30;
    reg63=reg13+reg63; reg11=reg11*reg45; reg22=reg21+reg22; reg4=reg90+reg4; reg39=reg39*reg45;
    reg13=reg101*reg83; reg79=reg48+reg79; reg15=reg110*reg131; reg21=reg131*reg65; reg48=reg81*reg83;
    reg71=reg92*reg131; reg80=reg83*reg84; reg69=reg76-reg69; reg58=reg23+reg58; reg86=reg86*reg1;
    reg18=reg6+reg18; reg12=reg77+reg12; reg14=reg82+reg14; reg68=reg27+reg68; reg88=reg88*reg1;
    reg61=reg67+reg61; reg0=reg26*reg0; reg35=reg116*reg35; reg62=reg29+reg62; reg66=reg31+reg66;
    reg42=reg87+reg42; reg37=reg44+reg37; reg21=reg79+reg21; reg42=reg42*reg58; reg39=reg4+reg39;
    reg37=reg37*reg58; reg66=reg66*reg45; reg11=reg22+reg11; reg4=reg101*reg92; reg6=reg110*reg81;
    reg22=reg127*reg83; reg13=reg15-reg13; reg15=reg110*reg59; reg23=reg102*reg83; reg48=reg71-reg48;
    reg26=reg92*reg59; reg80=reg69+reg80; reg27=reg14*reg68; reg29=reg63*reg68; reg31=reg30*reg16;
    reg44=reg12*reg16; reg35=reg62+reg35; reg88=reg61+reg88; reg54=reg54*reg0; reg86=reg18+reg86;
    reg78=reg78*reg0; reg25=reg117+reg25; reg46=reg46*reg1; reg120=reg120/reg80; reg23=reg26-reg23;
    reg42=reg39+reg42; reg13=reg13/reg80; reg37=reg11+reg37; reg22=reg15-reg22; reg44=reg27-reg44;
    reg11=reg14*reg30; reg15=reg127*reg92; reg18=reg12*reg63; reg4=reg6-reg4; reg6=reg110*reg102;
    reg26=1-var_inter[0]; reg27=1-var_inter[1]; reg66=reg21+reg66; reg35=reg35*reg58; reg129=reg129/reg80;
    reg32=reg32*reg0; reg46=reg25+reg46; reg78=reg86+reg78; reg54=reg88+reg54; reg31=reg29-reg31;
    reg48=reg48/reg80; reg21=elem.pos(1)[0]*var_inter[0]; reg84=reg84/reg80; reg25=reg26*elem.pos(0)[0]; reg29=reg26*elem.pos(0)[1];
    reg22=reg22/reg80; reg39=elem.pos(1)[1]*var_inter[0]; reg23=reg23/reg80; reg4=reg4/reg80; reg15=reg6-reg15;
    reg6=reg27*elem.pos(0)[1]; reg61=reg27*elem.pos(1)[1]; reg62=reg27*elem.pos(0)[0]; reg67=reg27*elem.pos(1)[0]; reg35=reg66+reg35;
    reg129=reg129*reg37; reg120=reg120*reg42; reg48=reg48*reg37; reg13=reg13*reg42; reg66=reg78*reg44;
    reg18=reg11-reg18; reg32=reg46+reg32; reg11=reg31*reg54; reg4=reg4*reg35; reg48=reg13-reg48;
    reg61=reg61-reg6; reg13=elem.pos(2)[1]*var_inter[0]; reg46=reg39+reg29; reg69=reg54*reg16; reg37=reg23*reg37;
    reg23=var_inter[1]*elem.pos(2)[1]; reg71=reg63*reg32; reg42=reg22*reg42; reg22=reg21+reg25; reg80=reg15/reg80;
    reg15=reg78*reg16; reg66=reg11-reg66; reg11=reg32*reg18; reg67=reg67-reg62; reg76=var_inter[1]*elem.pos(2)[0];
    reg84=reg84*reg35; reg77=reg14*reg32; reg79=elem.pos(2)[0]*var_inter[0]; reg120=reg129-reg120; reg82=reg78*reg68;
    reg23=reg61+reg23; reg61=1-(*f.m).resolution; reg86=reg30*reg32; reg42=reg37-reg42; reg77=reg69-reg77;
    reg35=reg80*reg35; reg37=reg54*reg68; reg79=reg79-reg22; reg69=elem.pos(3)[0]*reg26; reg11=reg66+reg11;
    reg71=reg15-reg71; reg15=reg54*reg63; reg66=elem.pos(3)[1]*reg26; reg4=reg48-reg4; reg76=reg67+reg76;
    reg48=reg14*reg78; reg67=reg12*reg32; reg80=var_inter[1]*elem.pos(3)[0]; reg13=reg13-reg46; reg87=var_inter[1]*elem.pos(3)[1];
    reg120=reg84+reg120; reg84=(*f.m).resolution*reg140; reg76=reg76-reg80; reg88=(*f.m).resolution*reg43; reg42=reg35+reg42;
    reg67=reg37-reg67; reg77=reg77/reg11; reg35=reg54*reg30; reg37=reg12*reg78; reg4=reg61*reg4;
    reg48=reg15-reg48; reg13=reg66+reg13; reg69=reg79+reg69; reg71=reg71/reg11; reg86=reg82-reg86;
    reg120=reg61*reg120; reg23=reg23-reg87; reg15=(*f.m).resolution*reg71; reg66=(*f.m).resolution*reg65; reg83=reg83*reg61;
    reg59=reg59*reg61; reg42=reg61*reg42; reg120=reg88+reg120; reg79=(*f.m).resolution*reg77; reg4=reg84+reg4;
    reg44=reg44/reg11; reg67=reg67/reg11; reg37=reg35-reg37; reg86=reg86/reg11; reg35=reg76*reg13;
    reg48=reg48/reg11; reg82=reg23*reg69; reg31=reg31/reg11; reg82=reg35-reg82; reg127=reg127*reg61;
    reg92=reg92*reg61; reg102=reg102*reg61; reg131=reg131*reg61; reg83=reg15+reg83; reg79=reg59-reg79;
    reg37=reg37/reg11; reg11=reg18/reg11; reg42=reg66+reg42; reg120=(*f.m).deltaT*reg120; reg4=(*f.m).deltaT*reg4;
    reg15=(*f.m).resolution*reg86; reg18=(*f.m).resolution*reg67; reg35=(*f.m).resolution*reg44; reg59=(*f.m).resolution*reg31; reg110=reg110*reg61;
    reg66=(*f.m).resolution*reg48; reg131=reg66+reg131; reg42=(*f.m).deltaT*reg42; reg102=reg18+reg102; reg15=reg92-reg15;
    reg35=reg127-reg35; reg110=reg59+reg110; reg81=reg81*reg61; reg13=reg13/reg82; reg18=(*f.m).resolution*reg37;
    reg61=reg101*reg61; reg59=reg83*reg120; reg23=reg23/reg82; reg69=reg69/reg82; reg66=reg79*reg4;
    reg76=reg76/reg82; reg84=(*f.m).resolution*reg11; reg88=reg59+reg66; reg89=reg131*reg42; reg90=reg110*reg120;
    reg92=reg15*reg120; reg93=reg35*reg4; reg94=reg102*reg4; reg95=var_inter[1]*reg13; reg97=var_inter[1]*reg69;
    reg101=reg27*reg13; reg104=var_inter[0]*reg76; reg106=var_inter[0]*reg23; reg107=reg26*reg76; reg108=reg27*reg69;
    reg109=reg26*reg23; reg61=reg84+reg61; reg18=reg81-reg18; reg81=reg92+reg94; reg84=reg104+reg108;
    reg111=reg18*reg42; reg112=reg88+reg89; reg115=reg61*reg42; reg117=reg93+reg90; reg118=reg97+reg107;
    reg119=reg101+reg106; reg124=reg95+reg109; reg125=reg81+reg111; reg126=0.5*reg119; reg127=reg109-reg101;
    reg128=reg108-reg107; reg129=2*reg112; reg130=0.5*reg84; reg133=reg117+reg115; reg134=reg95-reg106;
    reg135=0.5*reg124; reg136=0.5*reg118; reg138=reg104-reg97; reg141=reg130*reg129; reg143=reg119*reg133;
    reg144=reg124*reg133; T reg145=reg136*reg129; T reg146=reg126*reg129; T reg147=reg84*reg125; T reg148=reg118*reg125;
    T reg149=reg135*reg129; T reg150=0.5*reg138; T reg151=0.5*reg128; T reg152=0.5*reg134; T reg153=0.5*reg127;
    T reg154=reg147-reg146; T reg155=reg138*reg125; T reg156=reg134*reg133; T reg157=reg152*reg129; T reg158=reg150*reg129;
    T reg159=reg141-reg143; T reg160=reg144-reg145; T reg161=reg128*reg125; T reg162=reg153*reg129; T reg163=reg149-reg148;
    T reg164=reg151*reg129; T reg165=reg127*reg133; T reg166=reg155+reg157; reg160=reg82*reg160; reg163=reg82*reg163;
    T reg167=reg158+reg156; reg154=reg82*reg154; T reg168=reg165+reg164; reg159=reg82*reg159; T reg169=reg162+reg161;
    reg154=ponderation*reg154; T reg170=reg82*reg167; reg163=ponderation*reg163; reg159=ponderation*reg159; T reg171=reg82*reg166;
    T reg172=reg82*reg168; T reg173=reg82*reg169; reg160=ponderation*reg160; sollicitation[indices[3]+1]+=-reg163; reg163=ponderation*reg172;
    sollicitation[indices[0]+0]+=reg163; sollicitation[indices[3]+0]+=-reg160; reg160=ponderation*reg173; sollicitation[indices[0]+1]+=reg160; T reg174=ponderation*reg171;
    sollicitation[indices[2]+1]+=reg174; sollicitation[indices[1]+0]+=-reg159; reg159=ponderation*reg170; sollicitation[indices[2]+0]+=reg159; sollicitation[indices[1]+1]+=-reg154;
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

#ifndef elasticity_orthotropy_stat_Qstat_read_material_to_mesh
#define elasticity_orthotropy_stat_Qstat_read_material_to_mesh
template<class TM, class T, bool wont_add_nz>
void read_material_to_mesh_(const XmlNode &n, Formulation<TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f){ 
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

    if(n.has_attribute("elastic_modulus_1"))  
        n.get_attribute("elastic_modulus_1", f.m->elastic_modulus_1 ); 
    else  
        std::cerr << "Warning using default value of elastic_modulus_1 : " << f.m->elastic_modulus_1 << std::endl; 

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

