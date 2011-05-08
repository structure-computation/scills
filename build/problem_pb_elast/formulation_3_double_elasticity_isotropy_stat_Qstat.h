
#include "formulation/formulation.h"
namespace LMT {
#ifndef ELASTICITY_ISOTROPY_STAT_QSTAT_3_NUM_0_CARAC_H
#define ELASTICITY_ISOTROPY_STAT_QSTAT_3_NUM_0_CARAC_H
#ifndef STRUCT_ELASTICITY_ISOTROPY_STAT_QSTAT
#define STRUCT_ELASTICITY_ISOTROPY_STAT_QSTAT
struct elasticity_isotropy_stat_Qstat {
  static const char *name() { return "elasticity_isotropy_stat_Qstat"; }
};
#endif // ELASTICITY_ISOTROPY_STAT_QSTAT

template<class P_T>
class CaracFormulation<elasticity_isotropy_stat_Qstat,3,P_T>  {
public:
  typedef P_T T;
  static const char *name() { return "elasticity_isotropy_stat_Qstat"; }
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
  
  static const unsigned nb_nodal_unknowns = 3;
  template<class TE,class TTs,class Tvec>
  inline static void set_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice) {
    node.dep[0]=vecs[0][indice+0]; node.dep[1]=vecs[0][indice+1]; node.dep[2]=vecs[0][indice+2];
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
    node.dep[1]=vecs[0][indice+1]; node.dep[0]=vecs[0][indice+0]; node.dep[2]=vecs[0][indice+2];
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_nodal_initial_conditions(const TE &node,const TTs &f,Tvec &vecs,unsigned indice) {
    vecs[0][indice+2]=node.dep[2]; vecs[1][indice+2]=node.dep[2]; vecs[2][indice+2]=node.dep[2]; vecs[3][indice+2]=node.dep[2]; vecs[4][indice+2]=node.dep[2];
    vecs[0][indice+1]=node.dep[1]; vecs[1][indice+1]=node.dep[1]; vecs[2][indice+1]=node.dep[1]; vecs[3][indice+1]=node.dep[1]; vecs[4][indice+1]=node.dep[1];
    vecs[0][indice+0]=node.dep[0]; vecs[1][indice+0]=node.dep[0]; vecs[2][indice+0]=node.dep[0]; vecs[3][indice+0]=node.dep[0]; vecs[4][indice+0]=node.dep[0];
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_nodal_error(const TE &node,const TTs &f,const Tvec &vecs,int indice) {
    T reg0=vecs[1][indice+1]-vecs[0][indice+1]; T reg1=vecs[1][indice+0]-vecs[0][indice+0]; T reg2=vecs[1][indice+2]-vecs[0][indice+2]; reg0=abs(reg0); reg1=abs(reg1);
    reg1=max(reg0,reg1); reg2=abs(reg2); return max(reg2,reg1);
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_nodal(const TE &node,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
    old_vec[indice+0]=vecs[1][indice+0]; old_vec[indice+2]=vecs[1][indice+2]; old_vec[indice+1]=vecs[1][indice+1];
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
#endif // ELASTICITY_ISOTROPY_STAT_QSTAT_3_NUM_0_CARAC_H
extern "C" void apply_on_elements_after_solve_0_elasticity_isotropy_stat_Qstat_Tetra_0( double * );
extern "C" void apply_on_elements_after_solve_1_elasticity_isotropy_stat_Qstat_Tetra_1( double * );
extern "C" void apply_on_elements_after_solve_2_elasticity_isotropy_stat_Qstat_Tetra_2( double * );
extern "C" void apply_on_elements_after_solve_3_elasticity_isotropy_stat_Qstat_Tetra_3( double * );
extern "C" void apply_on_elements_after_solve_4_elasticity_isotropy_stat_Qstat_Tetra_4( double * );
extern "C" void apply_on_elements_after_solve_5_elasticity_isotropy_stat_Qstat_Tetra_5( double * );
extern "C" void apply_on_elements_after_solve_6_elasticity_isotropy_stat_Qstat_Tetra_6( double * );
extern "C" void apply_on_elements_after_solve_7_elasticity_isotropy_stat_Qstat_Tetra_7( double * );
extern "C" void apply_on_elements_after_solve_8_elasticity_isotropy_stat_Qstat_Tetra_8( double * );
extern "C" void apply_on_elements_after_solve_9_elasticity_isotropy_stat_Qstat_Tetra_9( double * );
extern "C" void apply_on_elements_after_solve_10_elasticity_isotropy_stat_Qstat_Tetra_10( double * );
extern "C" void apply_on_elements_after_solve_11_elasticity_isotropy_stat_Qstat_Tetra_11( double * );
extern "C" void apply_on_elements_after_solve_12_elasticity_isotropy_stat_Qstat_Tetra_12( double * );
extern "C" void apply_on_elements_after_solve_13_elasticity_isotropy_stat_Qstat_Tetra_13( double * );
extern "C" void apply_on_elements_after_solve_14_elasticity_isotropy_stat_Qstat_Tetra_14( double * );
class Tetra;
template<unsigned A,class B,class C> class Node;
template<class A,class B,class C,class D,unsigned E> class Element;

// Carac for ...
template<class P_T_pos,class P_ND,class TED,unsigned nim,class TM,class T>
class CaracFormulationForElement<elasticity_isotropy_stat_Qstat,Element<Tetra,DefaultBehavior,Node<3,P_T_pos,P_ND>,TED,nim>,TM,T> {
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
    T reg0=1+(*f.m).poisson_ratio; T reg1=elem.pos(3)[2]-elem.pos(0)[2]; T reg2=elem.pos(3)[1]-elem.pos(0)[1]; T reg3=elem.pos(2)[2]-elem.pos(0)[2]; T reg4=elem.pos(2)[1]-elem.pos(0)[1];
    T reg5=elem.pos(1)[2]-elem.pos(0)[2]; T reg6=elem.pos(1)[1]-elem.pos(0)[1]; reg0=reg0/(*f.m).elastic_modulus; T reg7=reg5*reg2; T reg8=reg3*reg2;
    T reg9=reg6*reg1; T reg10=reg4*reg1; T reg11=elem.pos(1)[0]-elem.pos(0)[0]; reg8=reg10-reg8; reg7=reg9-reg7;
    reg9=reg6*reg3; reg10=reg5*reg4; T reg12=pow(reg0,2); T reg13=elem.pos(2)[0]-elem.pos(0)[0]; reg10=reg9-reg10;
    reg0=reg0*reg12; reg9=elem.pos(3)[0]-elem.pos(0)[0]; T reg14=reg13*reg7; T reg15=1.0/(*f.m).elastic_modulus; T reg16=reg11*reg8;
    T reg17=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg18=reg13*reg1; T reg19=reg3*reg9; T reg20=reg9*reg10; reg14=reg16-reg14;
    reg1=reg11*reg1; reg16=reg17*reg0; reg0=reg15*reg0; T reg21=reg5*reg9; reg5=reg5*reg13;
    reg20=reg14+reg20; reg14=reg13*reg2; T reg22=reg17*reg0; reg19=reg18-reg19; reg18=reg17*reg16;
    reg0=reg15*reg0; reg3=reg11*reg3; T reg23=reg4*reg9; reg2=reg11*reg2; reg21=reg1-reg21;
    reg9=reg6*reg9; reg1=PNODE(2).dep[0]-PNODE(0).dep[0]; T reg24=PNODE(1).dep[1]-PNODE(0).dep[1]; T reg25=PNODE(2).dep[1]-PNODE(0).dep[1]; T reg26=PNODE(1).dep[0]-PNODE(0).dep[0];
    reg13=reg6*reg13; reg5=reg3-reg5; reg4=reg11*reg4; reg8=reg8/reg20; reg19=reg19/reg20;
    reg23=reg14-reg23; reg7=reg7/reg20; reg21=reg21/reg20; reg9=reg2-reg9; reg16=reg15*reg16;
    reg22=reg18+reg22; reg0=reg0-reg18; reg23=reg23/reg20; reg2=reg21*reg25; reg3=reg19*reg24;
    reg6=PNODE(3).dep[1]-PNODE(0).dep[1]; reg11=reg15*reg0; reg14=reg7*reg1; T reg27=reg15*reg12; T reg28=PNODE(3).dep[0]-PNODE(0).dep[0];
    reg12=reg17*reg12; reg16=reg18+reg16; reg18=reg17*reg22; T reg29=PNODE(1).dep[2]-PNODE(0).dep[2]; T reg30=PNODE(2).dep[2]-PNODE(0).dep[2];
    T reg31=reg8*reg26; reg13=reg4-reg13; reg5=reg5/reg20; reg9=reg9/reg20; reg10=reg10/reg20;
    reg13=reg13/reg20; reg14=reg31-reg14; reg18=reg11-reg18; reg4=reg17*reg16; reg11=reg9*reg30;
    reg31=reg23*reg29; T reg32=PNODE(3).dep[2]-PNODE(0).dep[2]; reg3=reg2-reg3; reg2=reg5*reg6; T reg33=reg17*reg27;
    T reg34=reg17*reg12; reg27=reg15*reg27; T reg35=reg10*reg28; T reg36=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; T reg37=vectors[0][indices[1]+0]-vectors[0][indices[0]+0];
    T reg38=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg11=reg31-reg11; reg31=reg13*reg32; T reg39=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg2=reg3-reg2;
    reg3=reg7*reg25; T reg40=reg8*reg24; T reg41=reg21*reg1; T reg42=reg19*reg26; reg35=reg14+reg35;
    reg33=reg34+reg33; reg27=reg27-reg34; reg12=reg15*reg12; reg4=reg18-reg4; reg14=(*f.m).alpha*(*f.m).deltaT;
    reg3=reg40-reg3; reg18=reg10*reg6; reg35=reg35-reg14; reg40=reg34+reg12; T reg43=reg19*reg39;
    reg2=reg2-reg14; T reg44=reg8*reg29; T reg45=reg7*reg30; T reg46=vectors[0][indices[1]+2]-vectors[0][indices[0]+2]; reg11=reg31+reg11;
    reg31=reg8*reg37; T reg47=reg21*reg38; T reg48=vectors[0][indices[2]+2]-vectors[0][indices[0]+2]; T reg49=vectors[0][indices[3]+0]-vectors[0][indices[0]+0]; T reg50=reg7*reg36;
    reg26=reg23*reg26; reg33=reg17*reg33; reg27=reg15*reg27; reg15=vectors[0][indices[3]+1]-vectors[0][indices[0]+1]; T reg51=reg5*reg28;
    reg1=reg9*reg1; reg42=reg41-reg42; reg22=reg22/reg4; reg0=reg0/reg4; reg45=reg44-reg45;
    reg41=reg10*reg32; reg29=reg19*reg29; reg30=reg21*reg30; reg44=reg9*reg48; T reg52=vectors[0][indices[3]+2]-vectors[0][indices[0]+2];
    T reg53=reg23*reg46; reg50=reg31-reg50; reg31=reg0*reg35; reg16=reg16/reg4; T reg54=reg10*reg49;
    T reg55=reg5*reg15; reg33=reg27-reg33; reg43=reg47-reg43; reg18=reg3+reg18; reg11=reg11-reg14;
    reg40=reg17*reg40; reg1=reg26-reg1; reg24=reg23*reg24; reg25=reg9*reg25; reg3=reg22*reg35;
    reg28=reg13*reg28; reg51=reg42-reg51; reg26=reg22*reg2; reg27=reg0*reg2; reg42=reg11*reg16;
    reg47=reg16*reg2; reg44=reg53-reg44; reg27=reg3+reg27; reg54=reg50+reg54; elem.epsilon[0][0]=reg54;
    reg18=reg51+reg18; reg50=reg13*reg52; reg55=reg43-reg55; elem.epsilon[0][1]=reg55; reg32=reg5*reg32;
    reg29=reg30-reg29; reg41=reg45+reg41; reg26=reg31+reg26; reg25=reg24-reg25; reg6=reg13*reg6;
    reg40=reg33-reg40; reg1=reg28+reg1; reg44=reg50+reg44; elem.epsilon[0][2]=reg44; reg24=reg54+reg55;
    reg28=reg11*reg0; reg47=reg3+reg47; reg4=reg40/reg4; reg27=reg42+reg27; reg25=reg6+reg25;
    reg32=reg29-reg32; reg18=0.5*reg18; reg41=reg1+reg41; reg42=reg26+reg42; reg1=reg4*reg18;
    reg28=reg47+reg28; reg41=0.5*reg41; reg32=reg25+reg32; reg2=reg2*reg27; reg35=reg35*reg42;
    reg3=reg8*reg39; reg6=reg19*reg37; reg24=reg44+reg24; reg25=reg21*reg36; reg26=reg7*reg38;
    reg35=reg2+reg35; reg2=reg4*reg41; reg29=reg7*reg48; reg30=reg5*reg49; reg6=reg25-reg6;
    reg1=2*reg1; reg37=reg23*reg37; reg24=reg24/3; reg36=reg9*reg36; reg11=reg11*reg28;
    reg32=0.5*reg32; reg26=reg3-reg26; reg3=reg10*reg15; reg25=reg8*reg46; reg36=reg37-reg36;
    reg30=reg6-reg30; reg3=reg26+reg3; reg49=reg13*reg49; reg11=reg35+reg11; reg18=reg1*reg18;
    reg6=reg55-reg24; reg26=reg54-reg24; reg46=reg19*reg46; reg48=reg21*reg48; reg38=reg9*reg38;
    reg39=reg23*reg39; reg2=2*reg2; reg31=reg10*reg52; reg29=reg25-reg29; reg25=reg4*reg32;
    reg6=pow(reg6,2); reg36=reg49+reg36; reg3=reg30+reg3; reg38=reg39-reg38; reg15=reg13*reg15;
    reg31=reg29+reg31; reg46=reg48-reg46; reg52=reg5*reg52; reg24=reg44-reg24; reg26=pow(reg26,2);
    reg41=reg41*reg2; reg25=2*reg25; reg18=reg11+reg18; reg41=reg18+reg41; reg32=reg32*reg25;
    reg24=pow(reg24,2); reg11=0.5*reg3; elem.epsilon[0][3]=reg11; reg52=reg46-reg52; reg6=reg26+reg6;
    reg38=reg15+reg38; reg31=reg36+reg31; reg32=reg41+reg32; reg52=reg38+reg52; reg15=0.5*reg31;
    elem.epsilon[0][4]=reg15; reg3=reg3*reg11; reg24=reg6+reg24; reg32=reg20*reg32; reg3=reg24+reg3;
    reg6=0.5*reg52; elem.epsilon[0][5]=reg6; reg31=reg31*reg15; reg18=0.041666666666666664354*reg32; reg32=0.083333333333333328707*reg32;
    reg54=reg54-reg14; reg31=reg3+reg31; reg52=reg52*reg6; reg55=reg55-reg14; reg3=reg0*reg55;
    reg32=reg18+reg32; reg20=reg22*reg54; reg52=reg31+reg52; reg14=reg44-reg14; reg54=reg0*reg54;
    reg24=reg16*reg55; reg55=reg22*reg55; reg32=reg18+reg32; reg55=reg54+reg55; reg16=reg16*reg14;
    reg3=reg20+reg3; reg52=1.5*reg52; reg24=reg20+reg24; reg14=reg0*reg14; elem.sigma_von_mises=pow(reg52,0.5);
    elem.sigma[0][5]=reg4*reg6; elem.ener=reg32/2; elem.sigma[0][4]=reg4*reg15; elem.sigma[0][0]=reg55+reg16; elem.sigma[0][3]=reg4*reg11;
    elem.sigma[0][1]=reg16+reg3; elem.sigma[0][2]=reg24+reg14;
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=elem.pos(3)[2]-elem.pos(0)[2]; T reg3=elem.pos(1)[1]-elem.pos(0)[1];
    T reg4=elem.pos(1)[2]-elem.pos(0)[2]; T reg5=elem.pos(2)[1]-elem.pos(0)[1]; T reg6=elem.pos(2)[2]-elem.pos(0)[2]; T reg7=elem.pos(3)[1]-elem.pos(0)[1]; T reg8=reg4*reg7;
    T reg9=reg5*reg2; T reg10=reg3*reg2; T reg11=reg6*reg7; T reg12=1.0/(*f.m).elastic_modulus; T reg13=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    reg0=reg0*reg1; T reg14=reg12*reg1; T reg15=reg3*reg6; T reg16=reg13*reg0; reg0=reg12*reg0;
    reg8=reg10-reg8; reg11=reg9-reg11; reg9=reg4*reg5; reg1=reg13*reg1; reg10=elem.pos(1)[0]-elem.pos(0)[0];
    T reg17=elem.pos(2)[0]-elem.pos(0)[0]; T reg18=reg13*reg14; T reg19=reg13*reg1; reg14=reg12*reg14; T reg20=reg13*reg0;
    T reg21=reg13*reg16; reg9=reg15-reg9; reg0=reg12*reg0; reg15=reg17*reg8; T reg22=reg10*reg11;
    T reg23=elem.pos(3)[0]-elem.pos(0)[0]; reg14=reg14-reg19; reg1=reg12*reg1; reg18=reg19+reg18; reg16=reg12*reg16;
    reg20=reg21+reg20; reg15=reg22-reg15; reg22=reg23*reg9; reg0=reg0-reg21; T reg24=reg4*reg23;
    T reg25=reg17*reg2; T reg26=reg6*reg23; T reg27=reg17*reg7; reg7=reg10*reg7; reg2=reg10*reg2;
    T reg28=reg5*reg23; reg23=reg3*reg23; reg3=reg3*reg17; reg5=reg10*reg5; reg24=reg2-reg24;
    reg17=reg4*reg17; reg28=reg27-reg28; reg26=reg25-reg26; reg23=reg7-reg23; reg6=reg10*reg6;
    reg22=reg15+reg22; reg2=reg19+reg1; reg4=reg13*reg20; reg18=reg13*reg18; reg14=reg12*reg14;
    reg16=reg21+reg16; reg12=reg12*reg0; reg28=reg28/reg22; reg8=reg8/reg22; reg24=reg24/reg22;
    reg26=reg26/reg22; reg7=reg13*reg16; reg11=reg11/reg22; reg4=reg12-reg4; reg23=reg23/reg22;
    reg3=reg5-reg3; reg18=reg14-reg18; reg2=reg13*reg2; reg17=reg6-reg17; reg9=reg9/reg22;
    reg17=reg17/reg22; reg5=reg23-reg28; reg3=reg3/reg22; reg6=reg26-reg24; reg7=reg4-reg7;
    reg4=reg8-reg11; reg2=reg18-reg2; reg5=reg5-reg3; reg10=0.5*reg26; reg4=reg4-reg9;
    reg12=0.5*reg17; reg13=0.5*reg9; reg6=reg17+reg6; reg14=0.5*reg8; reg15=0.5*reg11;
    reg18=0.5*reg24; reg21=(*f.m).alpha*(*f.m).deltaT; reg2=reg2/reg7; reg16=reg16/reg7; reg20=reg20/reg7;
    reg7=reg0/reg7; reg0=0.5*reg23; reg25=0.5*reg6; reg27=reg2*reg13; T reg29=0.5*reg28;
    T reg30=reg20*reg21; T reg31=reg7*reg21; T reg32=0.5*reg4; T reg33=reg16*reg21; T reg34=0.5*reg5;
    T reg35=reg2*reg18; T reg36=reg2*reg12; T reg37=reg2*reg15; T reg38=0.5*reg3; T reg39=reg2*reg10;
    T reg40=reg2*reg14; T reg41=reg7*reg24; T reg42=reg7*reg26; T reg43=reg2*reg32; T reg44=reg33+reg30;
    T reg45=reg2*reg38; reg27=2*reg27; T reg46=reg31+reg30; T reg47=2*reg40; T reg48=2*reg39;
    T reg49=reg7*reg9; T reg50=reg7*reg17; reg35=2*reg35; T reg51=reg7*reg11; T reg52=reg2*reg0;
    reg37=2*reg37; T reg53=2*reg36; T reg54=reg7*reg8; T reg55=reg2*reg29; T reg56=1-var_inter[0];
    T reg57=reg7*reg3; T reg58=reg2*reg34; T reg59=reg7*reg23; T reg60=reg7*reg28; T reg61=reg2*reg25;
    T reg62=reg20*reg8; T reg63=reg12*reg35; T reg64=reg8*reg49; reg45=2*reg45; T reg65=reg9*reg54;
    T reg66=reg16*reg23; T reg67=reg15*reg47; T reg68=reg20*reg24; T reg69=reg26*reg41; T reg70=reg28*reg59;
    T reg71=reg10*reg35; T reg72=reg11*reg54; T reg73=reg18*reg48; T reg74=2*reg52; T reg75=reg8*reg51;
    T reg76=reg18*reg53; T reg77=reg7*reg6; T reg78=reg20*reg4; T reg79=reg23*reg57; T reg80=reg13*reg47;
    T reg81=reg23*reg60; T reg82=reg17*reg41; T reg83=reg16*reg3; T reg84=reg7*reg5; T reg85=reg16*reg26;
    T reg86=reg20*reg11; T reg87=reg14*reg27; T reg88=reg24*reg50; T reg89=reg16*reg24; T reg90=reg20*reg17;
    T reg91=reg20*reg9; T reg92=reg14*reg37; T reg93=reg24*reg42; T reg94=reg16*reg17; reg58=2*reg58;
    reg56=reg56-var_inter[1]; reg55=2*reg55; T reg95=reg16*reg5; T reg96=reg33+reg46; T reg97=reg20*reg26;
    reg61=2*reg61; T reg98=reg16*reg28; reg43=2*reg43; T reg99=reg7*reg4; T reg100=reg3*reg59;
    T reg101=reg31+reg44; T reg102=reg24*reg41; T reg103=reg14*reg35; T reg104=reg5*reg59; T reg105=reg8*reg68;
    T reg106=reg24*reg62; T reg107=reg93+reg92; T reg108=var_inter[2]*(*f.m).f_vol[1]; T reg109=reg17*reg96; T reg110=reg14*reg43;
    T reg111=var_inter[1]*(*f.m).f_vol[2]; T reg112=reg24*reg77; reg64=reg76+reg64; T reg113=reg0*reg45; T reg114=reg5*reg57;
    T reg115=reg0*reg27; T reg116=reg99*reg11; T reg117=reg10*reg61; T reg118=reg8*reg83; T reg119=reg99*reg4;
    T reg120=reg99*reg9; T reg121=reg6*reg50; T reg122=reg32*reg27; reg79=reg87+reg79; T reg123=reg23*reg91;
    T reg124=reg14*reg45; T reg125=reg23*reg59; T reg126=reg23*reg89; T reg127=reg18*reg74; T reg128=reg16*reg6;
    reg81=reg92+reg81; reg92=var_inter[1]*(*f.m).f_vol[0]; T reg129=reg8*reg96; T reg130=reg5*reg84; T reg131=reg23*reg86;
    T reg132=reg14*reg55; T reg133=reg23*reg84; T reg134=reg23*reg78; T reg135=var_inter[0]*(*f.m).f_vol[1]; T reg136=reg14*reg58;
    T reg137=reg5*reg60; reg87=reg88+reg87; T reg138=reg23*reg101; T reg139=reg5*reg62; T reg140=reg32*reg74;
    T reg141=reg0*reg35; T reg142=reg24*reg66; T reg143=reg14*reg47; T reg144=reg11*reg90; T reg145=reg10*reg27;
    T reg146=reg15*reg74; T reg147=reg28*reg62; T reg148=reg28*reg60; T reg149=reg10*reg55; T reg150=reg28*reg85;
    T reg151=reg28*reg84; T reg152=reg15*reg43; T reg153=reg26*reg77; T reg154=reg26*reg83; T reg155=reg29*reg53;
    T reg156=reg26*reg50; T reg157=reg15*reg27; T reg158=reg15*reg48; T reg159=reg26*reg86; T reg160=reg26*reg91;
    T reg161=reg15*reg53; T reg162=reg15*reg37; T reg163=reg26*reg42; reg69=reg67+reg69; T reg164=reg29*reg48;
    T reg165=reg26*reg98; T reg166=reg12*reg27; T reg167=reg9*reg90; T reg168=reg18*reg47; T reg169=reg11*reg51;
    T reg170=reg10*reg48; T reg171=reg8*reg54; T reg172=reg18*reg35; T reg173=reg0*reg37; T reg174=reg8*reg98;
    T reg175=reg11*reg97; T reg176=reg10*reg37; T reg177=reg0*reg55; reg75=reg73+reg75; T reg178=reg0*reg43;
    reg56=reg56-var_inter[2]; T reg179=reg72+reg71; T reg180=reg29*reg74; T reg181=reg8*reg95; reg99=reg99*reg8;
    T reg182=reg18*reg61; T reg183=reg28*reg57; T reg184=reg11*reg66; T reg185=reg12*reg61; T reg186=reg10*reg45;
    T reg187=reg28*reg94; T reg188=reg11*reg49; T reg189=reg10*reg53; T reg190=reg67+reg70; T reg191=reg20*reg6;
    T reg192=reg6*reg77; T reg193=reg38*reg74; T reg194=reg26*reg96; T reg195=reg32*reg43; T reg196=reg13*reg27;
    T reg197=reg13*reg37; T reg198=reg17*reg50; T reg199=reg3*reg94; T reg200=reg9*reg66; T reg201=reg17*reg86;
    T reg202=reg13*reg48; T reg203=reg38*reg47; reg84=reg3*reg84; T reg204=reg9*reg49; T reg205=reg12*reg53;
    T reg206=reg25*reg53; T reg207=reg34*reg47; T reg208=reg32*reg37; T reg209=reg6*reg42; reg57=reg3*reg57;
    T reg210=reg12*reg45; T reg211=reg4*reg66; reg49=reg4*reg49; reg77=reg17*reg77; T reg212=reg4*reg51;
    T reg213=reg38*reg53; T reg214=reg25*reg48; T reg215=reg17*reg83; T reg216=reg13*reg43; T reg217=reg13*reg74;
    reg82=reg80+reg82; T reg218=reg3*reg62; reg51=reg9*reg51; T reg219=reg12*reg48; T reg220=reg29*reg47;
    T reg221=reg17*reg98; reg60=reg3*reg60; T reg222=reg80+reg100; T reg223=reg9*reg97; T reg224=reg12*reg37;
    T reg225=reg38*reg48; T reg226=reg4*reg54; T reg227=reg13*reg53; T reg228=reg65+reg63; T reg229=reg3*reg85;
    T reg230=reg32*reg47; T reg231=reg25*reg35; reg41=reg6*reg41; T reg232=reg17*reg42; T reg233=reg25*reg61;
    T reg234=reg12*reg55; T reg235=reg17*reg91; T reg236=reg4*reg96; reg186=reg187+reg186; T reg237=reg8*reg97;
    T reg238=reg18*reg37; reg183=reg157+reg183; T reg239=reg75+reg177; reg57=reg196+reg57; T reg240=reg193+reg82;
    T reg241=reg38*reg35; reg235=reg227+reg235; reg99=reg182-reg99; T reg242=reg58*reg0; T reg243=reg18*reg43;
    reg178=reg181+reg178; T reg244=reg8*reg191; reg196=reg196+reg198; T reg245=reg17*reg66; T reg246=reg28*reg128;
    T reg247=reg58*reg15; T reg248=reg28*reg78; T reg249=reg155+reg154; T reg250=reg13*reg55; reg157=reg157+reg156;
    reg234=reg229+reg234; T reg251=reg3*reg91; reg160=reg161+reg160; reg63=reg63+reg222; T reg252=reg26*reg66;
    T reg253=reg29*reg35; T reg254=reg180+reg69; T reg255=reg12*reg74; T reg256=reg3*reg89; reg60=reg197+reg60;
    T reg257=reg3*reg101; T reg258=reg26*reg62; T reg259=reg15*reg35; T reg260=reg218+reg217; reg210=reg199+reg210;
    T reg261=reg15*reg45; T reg262=reg213+reg215; T reg263=reg28*reg91; reg71=reg71+reg190; T reg264=reg3*reg78;
    T reg265=reg13*reg58; T reg266=reg10*reg74; T reg267=reg28*reg89; T reg268=reg147+reg146; T reg269=reg3*reg128;
    reg148=reg162+reg148; T reg270=reg12*reg58; reg149=reg150+reg149; reg84=reg216+reg84; T reg271=reg15*reg55;
    T reg272=reg28*reg86; reg151=reg152+reg151; T reg273=reg58*reg10; T reg274=reg13*reg45; T reg275=reg3*reg86;
    T reg276=reg23*reg85; T reg277=reg194-reg135; T reg278=reg18*reg55; reg131=reg132+reg131; reg133=reg110+reg133;
    reg132=reg228+reg193; T reg279=reg23*reg128; T reg280=reg58*reg18; reg134=reg136+reg134; reg136=reg9*reg68;
    T reg281=reg12*reg47; T reg282=reg0*reg53; T reg283=reg24*reg83; T reg284=reg113+reg87; T reg285=reg129-reg92;
    T reg286=reg200+reg203; T reg287=reg24*reg96; T reg288=reg14*reg53; reg204=reg204+reg205; T reg289=reg24*reg91;
    T reg290=reg38*reg43; T reg291=reg9*reg95; T reg292=reg12*reg43; T reg293=reg9*reg191; T reg294=reg38*reg58;
    reg120=reg120-reg185; reg51=reg51+reg219; reg79=reg76+reg79; T reg295=reg38*reg55; T reg296=reg28*reg101;
    T reg297=reg23*reg94; T reg298=reg18*reg45; reg123=reg124+reg123; reg124=reg143+reg125; reg224=reg223+reg224;
    reg126=reg127+reg126; T reg299=reg9*reg98; T reg300=reg23*reg62; T reg301=reg14*reg74; reg81=reg73+reg81;
    T reg302=reg38*reg37; T reg303=reg38*reg61; T reg304=reg109-reg108; T reg305=reg17*reg95; T reg306=reg14*reg61;
    T reg307=reg24*reg78; reg115=reg118+reg115; reg201=reg202+reg201; T reg308=reg8*reg90; T reg309=reg18*reg27;
    reg113=reg64+reg113; T reg310=reg0*reg47; T reg311=reg8*reg66; reg105=reg168+reg105; reg197=reg197+reg232;
    T reg312=reg0*reg74; T reg313=reg172+reg171; T reg314=reg225+reg221; T reg315=reg13*reg35; reg173=reg174+reg173;
    T reg316=reg17*reg62; T reg317=reg6*reg96; reg141=reg142+reg141; reg102=reg102+reg143; T reg318=reg138-reg111;
    T reg319=reg11*reg96; T reg320=reg38*reg45; reg103=reg106+reg103; reg166=reg167+reg166; T reg321=reg38*reg27;
    T reg322=reg13*reg61; T reg323=reg17*reg78; T reg324=reg0*reg48; T reg325=reg24*reg98; reg177=reg177+reg107;
    T reg326=reg5*reg101; T reg327=reg9*reg96; T reg328=reg14*reg48; T reg329=reg24*reg86; T reg330=reg0*reg61;
    T reg331=reg24*reg95; reg110=reg112-reg110; reg77=reg216-reg77; reg169=reg169+reg170; reg112=reg4*reg68;
    reg216=reg29*reg43; T reg332=reg11*reg95; T reg333=reg25*reg47; T reg334=reg10*reg43; T reg335=reg11*reg191;
    T reg336=reg207+reg211; T reg337=reg58*reg29; reg116=reg116-reg117; T reg338=reg34*reg45; reg114=reg122+reg114;
    reg49=reg49-reg206; T reg339=reg4*reg90; T reg340=reg5*reg94; T reg341=reg25*reg45; T reg342=reg25*reg27;
    T reg343=reg32*reg45; T reg344=reg4*reg95; reg45=reg29*reg45; reg188=reg188+reg189; T reg345=reg34*reg55;
    reg212=reg212-reg214; T reg346=reg220+reg184; T reg347=reg4*reg97; T reg348=reg10*reg47; reg68=reg11*reg68;
    T reg349=reg25*reg37; T reg350=reg179+reg180; T reg351=reg34*reg37; T reg352=reg4*reg98; reg37=reg29*reg37;
    T reg353=reg11*reg98; reg176=reg175+reg176; T reg354=reg34*reg74; T reg355=reg231-reg226; T reg356=reg29*reg55;
    reg130=reg195+reg130; T reg357=reg208-reg209; reg128=reg5*reg128; T reg358=reg58*reg25; T reg359=reg34*reg48;
    reg98=reg6*reg98; T reg360=reg58*reg32; T reg361=reg5*reg78; T reg362=reg6*reg62; T reg363=reg6*reg83;
    T reg364=reg32*reg35; T reg365=reg34*reg53; reg41=reg41-reg230; reg122=reg122-reg121; reg35=reg34*reg35;
    T reg366=reg32*reg53; T reg367=reg6*reg66; T reg368=reg6*reg91; reg91=reg5*reg91; T reg369=reg34*reg27;
    T reg370=reg4*reg83; T reg371=reg230+reg104; T reg372=reg6*reg78; reg89=reg5*reg89; T reg373=reg25*reg74;
    T reg374=reg32*reg61; reg195=reg192+reg195; reg192=reg139+reg140; reg137=reg208+reg137; reg208=reg34*reg61;
    T reg375=reg6*reg95; T reg376=reg5*reg85; T reg377=reg25*reg55; T reg378=reg6*reg86; T reg379=reg32*reg48;
    reg55=reg32*reg55; reg86=reg5*reg86; T reg380=var_inter[2]*(*f.m).f_vol[2]; T reg381=var_inter[1]*(*f.m).f_vol[1]; reg162=reg162+reg163;
    T reg382=reg11*reg83; reg27=reg29*reg27; T reg383=reg25*reg43; reg191=reg4*reg191; T reg384=var_inter[2]*(*f.m).f_vol[0];
    reg159=reg158+reg159; T reg385=reg56*(*f.m).f_vol[0]; T reg386=reg56*(*f.m).f_vol[1]; T reg387=reg15*reg61; reg78=reg26*reg78;
    T reg388=reg56*(*f.m).f_vol[2]; reg95=reg26*reg95; reg61=reg29*reg61; reg119=reg119+reg233; T reg389=var_inter[0]*(*f.m).f_vol[0];
    T reg390=var_inter[0]*(*f.m).f_vol[2]; reg58=reg58*reg34; reg153=reg152-reg153; reg83=reg9*reg83; reg152=reg164+reg165;
    reg145=reg144+reg145; reg43=reg34*reg43; T reg391=reg22*reg152; reg102=reg312+reg102; reg368=reg368-reg366;
    reg374=reg372+reg374; reg197=reg295+reg197; reg372=reg22*reg192; T reg392=reg22*reg141; T reg393=reg381+reg287;
    reg157=reg45+reg157; reg289=reg289+reg288; reg377=reg377-reg376; T reg394=reg22*reg201; reg195=reg58+reg195;
    reg120=reg120+reg294; reg137=reg137-reg214; reg61=reg61-reg95; reg259=reg259+reg258; reg295=reg51+reg295;
    reg51=reg22*reg284; reg303=reg303-reg305; reg330=reg331-reg330; reg247=reg248+reg247; reg241=reg241+reg245;
    reg342=reg342-reg339; reg329=reg329+reg328; reg343=reg91+reg343; reg290=reg291+reg290; reg153=reg337+reg153;
    reg91=reg22*reg240; reg248=reg22*reg177; reg285=reg22*reg285; reg370=reg369+reg370; reg231=reg231-reg371;
    reg325=reg325+reg324; reg315=reg315+reg316; reg291=reg22*reg249; reg331=reg22*reg103; reg35=reg35-reg367;
    reg369=reg22*reg314; reg89=reg89-reg373; reg292=reg293-reg292; reg122=reg338+reg122; reg357=reg345+reg357;
    reg128=reg358+reg128; reg253=reg253+reg252; reg302=reg299+reg302; reg293=reg22*reg81; reg299=reg22*reg286;
    reg358=reg301+reg300; reg298=reg298+reg297; reg360=reg361+reg360; reg364=reg364-reg362; reg361=reg22*reg126;
    reg136=reg136+reg281; reg98=reg98-reg359; T reg395=reg22*reg254; reg172=reg172+reg124; T reg396=reg22*reg132;
    reg162=reg356+reg162; reg363=reg363-reg365; reg304=reg22*reg304; T reg397=reg22*reg123; reg283=reg283+reg282;
    reg77=reg294+reg77; reg41=reg41-reg354; reg318=reg22*reg318; reg375=reg208+reg375; reg208=reg22*reg134;
    reg323=reg322-reg323; reg55=reg86+reg55; reg279=reg280-reg279; reg86=reg22*reg160; reg280=reg22*reg224;
    reg294=reg22*reg159; reg133=reg182-reg133; reg321=reg83+reg321; reg378=reg378-reg379; reg130=reg233+reg130;
    reg83=reg22*reg131; reg182=reg22*reg166; reg233=reg22*reg79; reg278=reg278+reg276; reg204=reg204+reg320;
    reg322=reg380+reg257; T reg398=reg384+reg327; reg244=reg243-reg244; reg256=reg256+reg255; reg243=reg22*reg178;
    T reg399=reg22*reg260; reg37=reg353+reg37; reg353=reg22*reg149; reg352=reg351+reg352; reg351=reg22*reg239;
    reg60=reg219+reg60; T reg400=reg22*reg176; reg27=reg382+reg27; reg238=reg238+reg237; reg382=reg388+reg326;
    T reg401=reg22*reg234; reg356=reg169+reg356; reg169=reg22*reg173; reg355=reg355-reg354; T reg402=reg389+reg319;
    reg250=reg275+reg250; reg313=reg313+reg312; reg271=reg272+reg271; reg216=reg332+reg216; reg272=reg22*reg105;
    reg267=reg267+reg266; reg344=reg43+reg344; reg45=reg188+reg45; reg43=reg22*reg71; reg57=reg205+reg57;
    reg188=reg385+reg236; reg275=reg22*reg268; reg261=reg263+reg261; reg263=reg22*reg210; reg332=reg22*reg145;
    T reg403=reg22*reg346; T reg404=reg22*reg186; reg212=reg345+reg212; reg274=reg251+reg274; reg183=reg189+reg183;
    reg383=reg191+reg383; reg68=reg68+reg348; reg191=reg386+reg317; reg251=reg22*reg63; reg99=reg99-reg242;
    reg148=reg170+reg148; reg349=reg349-reg347; reg345=reg22*reg350; reg119=reg58+reg119; reg117=reg151-reg117;
    reg265=reg264+reg265; reg58=reg22*reg113; reg78=reg387-reg78; reg151=reg22*reg336; reg337=reg116+reg337;
    reg309=reg309+reg308; reg116=reg22*reg262; reg277=reg22*reg277; reg264=reg22*reg115; reg196=reg320+reg196;
    reg114=reg114-reg206; reg273=reg246-reg273; reg306=reg307-reg306; reg49=reg338+reg49; reg341=reg341-reg340;
    reg242=reg110-reg242; reg110=reg22*reg235; reg246=reg390+reg296; reg185=reg84-reg185; reg112=reg112-reg333;
    reg84=reg311+reg310; reg270=reg269-reg270; reg334=reg335-reg334; reg342=reg22*reg342; reg269=ponderation*reg91;
    reg185=reg22*reg185; reg307=reg22*reg188; reg302=reg22*reg302; reg320=reg22*reg393; reg364=reg22*reg364;
    reg335=ponderation*reg401; reg212=reg22*reg212; reg304=ponderation*reg304; reg274=reg22*reg274; reg270=reg22*reg270;
    reg338=ponderation*reg394; reg196=reg22*reg196; reg387=ponderation*reg251; T reg405=ponderation*reg396; T reg406=ponderation*reg263;
    T reg407=ponderation*reg280; reg49=reg22*reg49; T reg408=ponderation*reg110; reg41=reg22*reg41; reg295=reg22*reg295;
    T reg409=reg22*reg402; reg241=reg22*reg241; reg57=reg22*reg57; reg303=reg22*reg303; reg383=reg22*reg383;
    reg195=reg22*reg195; T reg410=reg22*reg322; reg344=reg22*reg344; T reg411=reg22*reg246; reg35=reg22*reg35;
    T reg412=ponderation*reg399; reg119=reg22*reg119; reg355=reg22*reg355; reg285=ponderation*reg285; reg318=ponderation*reg318;
    reg375=reg22*reg375; reg204=reg22*reg204; reg352=reg22*reg352; reg378=reg22*reg378; reg323=reg22*reg323;
    reg60=reg22*reg60; reg374=reg22*reg374; T reg413=ponderation*reg182; reg265=reg22*reg265; reg112=reg22*reg112;
    reg321=reg22*reg321; reg197=reg22*reg197; reg77=reg22*reg77; reg98=reg22*reg98; reg277=ponderation*reg277;
    T reg414=reg22*reg191; T reg415=ponderation*reg151; reg250=reg22*reg250; reg315=reg22*reg315; T reg416=reg22*reg382;
    reg136=reg22*reg136; reg349=reg22*reg349; T reg417=ponderation*reg116; reg370=reg22*reg370; reg256=reg22*reg256;
    T reg418=reg22*reg398; T reg419=ponderation*reg299; reg357=reg22*reg357; T reg420=ponderation*reg369; reg261=reg22*reg261;
    reg343=reg22*reg343; reg329=reg22*reg329; T reg421=ponderation*reg248; T reg422=ponderation*reg43; reg231=reg22*reg231;
    reg325=reg22*reg325; reg45=reg22*reg45; reg267=reg22*reg267; reg89=reg22*reg89; T reg423=ponderation*reg331;
    reg102=reg22*reg102; T reg424=ponderation*reg275; T reg425=ponderation*reg332; T reg426=ponderation*reg372; T reg427=ponderation*reg392;
    reg148=reg22*reg148; reg289=reg22*reg289; reg137=reg22*reg137; T reg428=ponderation*reg353; T reg429=ponderation*reg51;
    reg27=reg22*reg27; reg377=reg22*reg377; reg283=reg22*reg283; T reg430=ponderation*reg208; reg271=reg22*reg271;
    reg55=reg22*reg55; reg238=reg22*reg238; T reg431=ponderation*reg351; reg356=reg22*reg356; T reg432=ponderation*reg169;
    T reg433=ponderation*reg400; reg216=reg22*reg216; reg313=reg22*reg313; T reg434=ponderation*reg243; reg37=reg22*reg37;
    T reg435=ponderation*reg272; reg334=reg22*reg334; reg84=reg22*reg84; reg244=reg22*reg244; T reg436=ponderation*reg345;
    reg337=reg22*reg337; T reg437=ponderation*reg58; reg99=reg22*reg99; reg309=reg22*reg309; reg114=reg22*reg114;
    reg183=reg22*reg183; T reg438=ponderation*reg264; reg306=reg22*reg306; reg68=reg22*reg68; reg341=reg22*reg341;
    reg242=reg22*reg242; T reg439=ponderation*reg404; reg330=reg22*reg330; T reg440=ponderation*reg403; reg128=reg22*reg128;
    T reg441=ponderation*reg293; T reg442=ponderation*reg291; T reg443=reg22*reg358; reg157=reg22*reg157; reg360=reg22*reg360;
    T reg444=ponderation*reg361; reg61=reg22*reg61; reg172=reg22*reg172; T reg445=ponderation*reg86; reg363=reg22*reg363;
    T reg446=ponderation*reg397; T reg447=ponderation*reg294; reg298=reg22*reg298; reg253=reg22*reg253; reg122=reg22*reg122;
    T reg448=ponderation*reg233; T reg449=ponderation*reg395; reg162=reg22*reg162; reg368=reg22*reg368; reg120=reg22*reg120;
    reg292=reg22*reg292; reg259=reg22*reg259; T reg450=ponderation*reg391; reg290=reg22*reg290; reg279=reg22*reg279;
    reg117=reg22*reg117; reg78=reg22*reg78; reg133=reg22*reg133; reg273=reg22*reg273; reg130=reg22*reg130;
    T reg451=ponderation*reg83; reg247=reg22*reg247; reg153=reg22*reg153; reg278=reg22*reg278; T tmp_5_11=ponderation*reg183;
    T tmp_3_6=-reg436; T tmp_0_1=ponderation*reg383; T tmp_4_3=-reg447; T tmp_6_0=ponderation*reg99; T tmp_3_11=ponderation*reg27;
    reg27=ponderation*reg418; sollicitation[indices[3]+0]+=reg27; T tmp_5_1=ponderation*reg273; T tmp_11_7=ponderation*reg256; T tmp_4_7=-reg449;
    T tmp_6_1=ponderation*reg244; reg99=ponderation*reg416; sollicitation[indices[0]+2]+=reg99; T tmp_5_3=ponderation*reg271; T tmp_0_4=ponderation*reg349;
    T tmp_3_5=ponderation*reg37; T tmp_11_6=-reg412; T tmp_4_4=ponderation*reg162; sollicitation[indices[3]+1]+=-reg304; T tmp_5_6=-reg424;
    T tmp_6_2=-reg434; T tmp_4_6=ponderation*reg259; T tmp_3_4=-reg433; T tmp_5_2=ponderation*reg117; T tmp_11_5=ponderation*reg60;
    reg37=ponderation*reg409; sollicitation[indices[1]+0]+=reg37; reg60=ponderation*reg411; sollicitation[indices[1]+2]+=reg60; reg117=ponderation*reg410;
    sollicitation[indices[3]+2]+=reg117; T tmp_6_3=-reg431; T tmp_0_5=ponderation*reg352; T tmp_4_1=ponderation*reg153; T tmp_5_7=ponderation*reg267;
    sollicitation[indices[2]+0]+=-reg285; T tmp_3_10=-reg425; T tmp_11_11=ponderation*reg57; T tmp_4_10=ponderation*reg157; T tmp_5_5=ponderation*reg148;
    T tmp_3_9=ponderation*reg45; T tmp_5_8=-reg422; T tmp_0_2=ponderation*reg344; T tmp_3_8=-reg440; T tmp_11_10=-reg406;
    T tmp_5_0=ponderation*reg247; reg45=ponderation*reg414; sollicitation[indices[0]+1]+=reg45; reg57=ponderation*reg320; sollicitation[indices[2]+1]+=reg57;
    T tmp_5_9=ponderation*reg261; T tmp_4_2=ponderation*reg61; T tmp_4_0=ponderation*reg78; T tmp_4_9=-reg445; T tmp_5_10=-reg439;
    T tmp_5_4=-reg428; reg61=ponderation*reg307; sollicitation[indices[0]+0]+=reg61; T tmp_11_9=ponderation*reg274; T tmp_3_7=ponderation*reg68;
    T tmp_0_0=ponderation*reg119; T tmp_4_11=-reg442; T tmp_0_3=ponderation*reg212; sollicitation[indices[2]+2]+=-reg318; T tmp_4_8=ponderation*reg253;
    sollicitation[indices[1]+1]+=-reg277; T tmp_11_8=-reg387; T tmp_8_5=-reg441; T tmp_9_8=-reg419; T tmp_8_4=ponderation*reg278;
    T tmp_1_3=ponderation*reg378; T tmp_2_1=ponderation*reg128; T tmp_9_9=ponderation*reg204; T tmp_8_3=-reg451; T tmp_9_10=-reg413;
    T tmp_2_2=ponderation*reg130; T tmp_8_2=ponderation*reg133; T tmp_9_11=ponderation*reg321; T tmp_8_1=ponderation*reg279; T tmp_8_0=-reg430;
    T tmp_10_0=ponderation*reg323; T tmp_1_2=ponderation*reg375; T tmp_2_3=ponderation*reg55; T tmp_7_11=ponderation*reg283; T tmp_10_1=ponderation*reg77;
    T tmp_7_10=-reg429; T tmp_2_4=ponderation*reg377; T tmp_10_2=ponderation*reg303; T tmp_1_1=ponderation*reg195; T tmp_7_9=ponderation*reg289;
    T tmp_2_5=ponderation*reg137; T tmp_9_2=ponderation*reg290; T tmp_1_8=ponderation*reg35; T tmp_9_1=ponderation*reg292; T tmp_9_0=ponderation*reg120;
    T tmp_9_3=ponderation*reg295; T tmp_1_7=ponderation*reg41; T tmp_1_9=ponderation*reg368; T tmp_4_5=-reg450; T tmp_8_11=-reg448;
    T tmp_9_4=-reg407; T tmp_1_6=ponderation*reg364; T tmp_8_10=ponderation*reg298; T tmp_9_5=ponderation*reg302; T tmp_1_10=ponderation*reg122;
    T tmp_8_9=-reg446; T tmp_1_5=ponderation*reg98; T tmp_8_8=ponderation*reg172; T tmp_9_6=-reg405; T tmp_1_11=ponderation*reg363;
    T tmp_8_7=-reg444; T tmp_9_7=ponderation*reg136; T tmp_8_6=ponderation*reg443; T tmp_1_4=ponderation*reg357; T tmp_2_0=ponderation*reg360;
    T tmp_7_0=ponderation*reg306; T tmp_6_11=-reg438; T tmp_10_10=ponderation*reg196; T tmp_0_8=-reg415; T tmp_2_11=ponderation*reg114;
    T tmp_6_10=ponderation*reg309; T tmp_10_11=-reg417; T tmp_6_9=-reg437; T tmp_3_0=ponderation*reg337; T tmp_11_0=ponderation*reg265;
    T tmp_0_7=ponderation*reg112; T tmp_6_8=ponderation*reg84; T tmp_11_1=ponderation*reg270; T tmp_3_1=ponderation*reg334; T tmp_6_7=-reg435;
    T tmp_6_6=ponderation*reg313; T tmp_11_2=ponderation*reg185; T tmp_0_6=ponderation*reg355; T tmp_3_2=ponderation*reg216; T tmp_11_3=ponderation*reg250;
    T tmp_6_5=-reg432; T tmp_3_3=ponderation*reg356; T tmp_6_4=ponderation*reg238; T tmp_11_4=-reg335; T tmp_10_3=-reg338;
    T tmp_7_8=-reg427; T tmp_1_0=ponderation*reg374; T tmp_2_6=-reg426; T tmp_10_4=ponderation*reg197; T tmp_7_7=ponderation*reg102;
    T tmp_7_6=-reg423; T tmp_10_5=-reg420; T tmp_0_11=ponderation*reg370; T tmp_2_7=ponderation*reg89; T tmp_7_5=ponderation*reg325;
    T tmp_10_6=ponderation*reg315; T tmp_7_4=-reg421; T tmp_0_10=ponderation*reg342; T tmp_2_8=ponderation*reg231; T tmp_7_3=ponderation*reg329;
    T tmp_10_7=-reg269; T tmp_2_9=ponderation*reg343; T tmp_7_2=ponderation*reg330; T tmp_10_8=ponderation*reg241; T tmp_7_1=ponderation*reg242;
    T tmp_0_9=ponderation*reg49; T tmp_10_9=-reg408; T tmp_2_10=ponderation*reg341;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[0]+2) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_3;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_4;
    matrix(indices[0]+0,indices[1]+2) += tmp_0_5;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_7;
    matrix(indices[0]+0,indices[2]+2) += tmp_0_8;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_9;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_10;
    matrix(indices[0]+0,indices[3]+2) += tmp_0_11;
    matrix(indices[0]+1,indices[0]+0) += tmp_1_0;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[0]+2) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_3;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_4;
    matrix(indices[0]+1,indices[1]+2) += tmp_1_5;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_7;
    matrix(indices[0]+1,indices[2]+2) += tmp_1_8;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_9;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_10;
    matrix(indices[0]+1,indices[3]+2) += tmp_1_11;
    matrix(indices[0]+2,indices[0]+0) += tmp_2_0;
    matrix(indices[0]+2,indices[0]+1) += tmp_2_1;
    matrix(indices[0]+2,indices[0]+2) += tmp_2_2;
    matrix(indices[0]+2,indices[1]+0) += tmp_2_3;
    matrix(indices[0]+2,indices[1]+1) += tmp_2_4;
    matrix(indices[0]+2,indices[1]+2) += tmp_2_5;
    matrix(indices[0]+2,indices[2]+0) += tmp_2_6;
    matrix(indices[0]+2,indices[2]+1) += tmp_2_7;
    matrix(indices[0]+2,indices[2]+2) += tmp_2_8;
    matrix(indices[0]+2,indices[3]+0) += tmp_2_9;
    matrix(indices[0]+2,indices[3]+1) += tmp_2_10;
    matrix(indices[0]+2,indices[3]+2) += tmp_2_11;
    matrix(indices[1]+0,indices[0]+0) += tmp_3_0;
    matrix(indices[1]+0,indices[0]+1) += tmp_3_1;
    matrix(indices[1]+0,indices[0]+2) += tmp_3_2;
    matrix(indices[1]+0,indices[1]+0) += tmp_3_3;
    matrix(indices[1]+0,indices[1]+1) += tmp_3_4;
    matrix(indices[1]+0,indices[1]+2) += tmp_3_5;
    matrix(indices[1]+0,indices[2]+0) += tmp_3_6;
    matrix(indices[1]+0,indices[2]+1) += tmp_3_7;
    matrix(indices[1]+0,indices[2]+2) += tmp_3_8;
    matrix(indices[1]+0,indices[3]+0) += tmp_3_9;
    matrix(indices[1]+0,indices[3]+1) += tmp_3_10;
    matrix(indices[1]+0,indices[3]+2) += tmp_3_11;
    matrix(indices[1]+1,indices[0]+0) += tmp_4_0;
    matrix(indices[1]+1,indices[0]+1) += tmp_4_1;
    matrix(indices[1]+1,indices[0]+2) += tmp_4_2;
    matrix(indices[1]+1,indices[1]+0) += tmp_4_3;
    matrix(indices[1]+1,indices[1]+1) += tmp_4_4;
    matrix(indices[1]+1,indices[1]+2) += tmp_4_5;
    matrix(indices[1]+1,indices[2]+0) += tmp_4_6;
    matrix(indices[1]+1,indices[2]+1) += tmp_4_7;
    matrix(indices[1]+1,indices[2]+2) += tmp_4_8;
    matrix(indices[1]+1,indices[3]+0) += tmp_4_9;
    matrix(indices[1]+1,indices[3]+1) += tmp_4_10;
    matrix(indices[1]+1,indices[3]+2) += tmp_4_11;
    matrix(indices[1]+2,indices[0]+0) += tmp_5_0;
    matrix(indices[1]+2,indices[0]+1) += tmp_5_1;
    matrix(indices[1]+2,indices[0]+2) += tmp_5_2;
    matrix(indices[1]+2,indices[1]+0) += tmp_5_3;
    matrix(indices[1]+2,indices[1]+1) += tmp_5_4;
    matrix(indices[1]+2,indices[1]+2) += tmp_5_5;
    matrix(indices[1]+2,indices[2]+0) += tmp_5_6;
    matrix(indices[1]+2,indices[2]+1) += tmp_5_7;
    matrix(indices[1]+2,indices[2]+2) += tmp_5_8;
    matrix(indices[1]+2,indices[3]+0) += tmp_5_9;
    matrix(indices[1]+2,indices[3]+1) += tmp_5_10;
    matrix(indices[1]+2,indices[3]+2) += tmp_5_11;
    matrix(indices[2]+0,indices[0]+0) += tmp_6_0;
    matrix(indices[2]+0,indices[0]+1) += tmp_6_1;
    matrix(indices[2]+0,indices[0]+2) += tmp_6_2;
    matrix(indices[2]+0,indices[1]+0) += tmp_6_3;
    matrix(indices[2]+0,indices[1]+1) += tmp_6_4;
    matrix(indices[2]+0,indices[1]+2) += tmp_6_5;
    matrix(indices[2]+0,indices[2]+0) += tmp_6_6;
    matrix(indices[2]+0,indices[2]+1) += tmp_6_7;
    matrix(indices[2]+0,indices[2]+2) += tmp_6_8;
    matrix(indices[2]+0,indices[3]+0) += tmp_6_9;
    matrix(indices[2]+0,indices[3]+1) += tmp_6_10;
    matrix(indices[2]+0,indices[3]+2) += tmp_6_11;
    matrix(indices[2]+1,indices[0]+0) += tmp_7_0;
    matrix(indices[2]+1,indices[0]+1) += tmp_7_1;
    matrix(indices[2]+1,indices[0]+2) += tmp_7_2;
    matrix(indices[2]+1,indices[1]+0) += tmp_7_3;
    matrix(indices[2]+1,indices[1]+1) += tmp_7_4;
    matrix(indices[2]+1,indices[1]+2) += tmp_7_5;
    matrix(indices[2]+1,indices[2]+0) += tmp_7_6;
    matrix(indices[2]+1,indices[2]+1) += tmp_7_7;
    matrix(indices[2]+1,indices[2]+2) += tmp_7_8;
    matrix(indices[2]+1,indices[3]+0) += tmp_7_9;
    matrix(indices[2]+1,indices[3]+1) += tmp_7_10;
    matrix(indices[2]+1,indices[3]+2) += tmp_7_11;
    matrix(indices[2]+2,indices[0]+0) += tmp_8_0;
    matrix(indices[2]+2,indices[0]+1) += tmp_8_1;
    matrix(indices[2]+2,indices[0]+2) += tmp_8_2;
    matrix(indices[2]+2,indices[1]+0) += tmp_8_3;
    matrix(indices[2]+2,indices[1]+1) += tmp_8_4;
    matrix(indices[2]+2,indices[1]+2) += tmp_8_5;
    matrix(indices[2]+2,indices[2]+0) += tmp_8_6;
    matrix(indices[2]+2,indices[2]+1) += tmp_8_7;
    matrix(indices[2]+2,indices[2]+2) += tmp_8_8;
    matrix(indices[2]+2,indices[3]+0) += tmp_8_9;
    matrix(indices[2]+2,indices[3]+1) += tmp_8_10;
    matrix(indices[2]+2,indices[3]+2) += tmp_8_11;
    matrix(indices[3]+0,indices[0]+0) += tmp_9_0;
    matrix(indices[3]+0,indices[0]+1) += tmp_9_1;
    matrix(indices[3]+0,indices[0]+2) += tmp_9_2;
    matrix(indices[3]+0,indices[1]+0) += tmp_9_3;
    matrix(indices[3]+0,indices[1]+1) += tmp_9_4;
    matrix(indices[3]+0,indices[1]+2) += tmp_9_5;
    matrix(indices[3]+0,indices[2]+0) += tmp_9_6;
    matrix(indices[3]+0,indices[2]+1) += tmp_9_7;
    matrix(indices[3]+0,indices[2]+2) += tmp_9_8;
    matrix(indices[3]+0,indices[3]+0) += tmp_9_9;
    matrix(indices[3]+0,indices[3]+1) += tmp_9_10;
    matrix(indices[3]+0,indices[3]+2) += tmp_9_11;
    matrix(indices[3]+1,indices[0]+0) += tmp_10_0;
    matrix(indices[3]+1,indices[0]+1) += tmp_10_1;
    matrix(indices[3]+1,indices[0]+2) += tmp_10_2;
    matrix(indices[3]+1,indices[1]+0) += tmp_10_3;
    matrix(indices[3]+1,indices[1]+1) += tmp_10_4;
    matrix(indices[3]+1,indices[1]+2) += tmp_10_5;
    matrix(indices[3]+1,indices[2]+0) += tmp_10_6;
    matrix(indices[3]+1,indices[2]+1) += tmp_10_7;
    matrix(indices[3]+1,indices[2]+2) += tmp_10_8;
    matrix(indices[3]+1,indices[3]+0) += tmp_10_9;
    matrix(indices[3]+1,indices[3]+1) += tmp_10_10;
    matrix(indices[3]+1,indices[3]+2) += tmp_10_11;
    matrix(indices[3]+2,indices[0]+0) += tmp_11_0;
    matrix(indices[3]+2,indices[0]+1) += tmp_11_1;
    matrix(indices[3]+2,indices[0]+2) += tmp_11_2;
    matrix(indices[3]+2,indices[1]+0) += tmp_11_3;
    matrix(indices[3]+2,indices[1]+1) += tmp_11_4;
    matrix(indices[3]+2,indices[1]+2) += tmp_11_5;
    matrix(indices[3]+2,indices[2]+0) += tmp_11_6;
    matrix(indices[3]+2,indices[2]+1) += tmp_11_7;
    matrix(indices[3]+2,indices[2]+2) += tmp_11_8;
    matrix(indices[3]+2,indices[3]+0) += tmp_11_9;
    matrix(indices[3]+2,indices[3]+1) += tmp_11_10;
    matrix(indices[3]+2,indices[3]+2) += tmp_11_11;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=elem.pos(3)[1]-elem.pos(0)[1]; T reg2=elem.pos(3)[2]-elem.pos(0)[2]; T reg3=elem.pos(2)[2]-elem.pos(0)[2];
    T reg4=pow(reg0,2); T reg5=elem.pos(2)[1]-elem.pos(0)[1]; T reg6=elem.pos(1)[1]-elem.pos(0)[1]; T reg7=elem.pos(1)[2]-elem.pos(0)[2]; T reg8=reg3*reg1;
    T reg9=1.0/(*f.m).elastic_modulus; T reg10=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg11=reg6*reg2; T reg12=reg5*reg2; reg0=reg0*reg4;
    T reg13=reg7*reg1; T reg14=elem.pos(2)[0]-elem.pos(0)[0]; T reg15=reg10*reg0; T reg16=elem.pos(1)[0]-elem.pos(0)[0]; reg0=reg9*reg0;
    T reg17=reg10*reg4; T reg18=reg7*reg5; T reg19=reg6*reg3; reg13=reg11-reg13; reg4=reg9*reg4;
    reg8=reg12-reg8; reg11=reg9*reg4; reg12=reg10*reg17; reg4=reg10*reg4; T reg20=elem.pos(3)[0]-elem.pos(0)[0];
    T reg21=reg10*reg0; T reg22=reg10*reg15; reg0=reg9*reg0; reg18=reg19-reg18; reg19=reg14*reg13;
    T reg23=reg16*reg8; T reg24=reg14*reg2; T reg25=reg3*reg20; reg2=reg16*reg2; reg4=reg12+reg4;
    T reg26=reg20*reg18; reg19=reg23-reg19; reg11=reg11-reg12; reg17=reg9*reg17; reg0=reg0-reg22;
    reg15=reg9*reg15; reg23=reg7*reg20; reg21=reg22+reg21; reg25=reg24-reg25; reg24=reg14*reg1;
    reg7=reg7*reg14; reg3=reg16*reg3; reg26=reg19+reg26; reg19=reg6*reg20; reg23=reg2-reg23;
    reg2=reg10*reg21; T reg27=reg9*reg0; reg1=reg16*reg1; reg11=reg9*reg11; reg4=reg10*reg4;
    reg20=reg5*reg20; reg9=reg12+reg17; reg15=reg22+reg15; reg25=reg25/reg26; reg8=reg8/reg26;
    reg14=reg6*reg14; reg7=reg3-reg7; reg5=reg16*reg5; reg20=reg24-reg20; reg4=reg11-reg4;
    reg19=reg1-reg19; reg9=reg10*reg9; reg10=reg10*reg15; reg13=reg13/reg26; reg2=reg27-reg2;
    reg23=reg23/reg26; reg7=reg7/reg26; reg18=reg18/reg26; reg19=reg19/reg26; reg20=reg20/reg26;
    reg14=reg5-reg14; reg1=reg25-reg23; reg3=reg13-reg8; reg9=reg4-reg9; reg10=reg2-reg10;
    reg2=reg19-reg20; reg4=0.5*reg7; reg5=0.5*reg18; reg6=0.5*reg13; reg3=reg3-reg18;
    reg1=reg7+reg1; reg14=reg14/reg26; reg11=0.5*reg23; reg16=(*f.m).alpha*(*f.m).deltaT; reg9=reg9/reg10;
    reg15=reg15/reg10; reg21=reg21/reg10; reg10=reg0/reg10; reg2=reg2-reg14; reg0=0.5*reg3;
    reg22=0.5*reg25; reg24=reg9*reg5; reg27=0.5*reg1; T reg28=reg21*reg16; T reg29=0.5*reg8;
    T reg30=reg10*reg16; T reg31=0.5*reg19; T reg32=reg15*reg16; T reg33=reg9*reg4; T reg34=reg9*reg6;
    T reg35=0.5*reg14; T reg36=reg9*reg11; T reg37=2*reg34; reg24=2*reg24; T reg38=reg30+reg28;
    T reg39=reg10*reg19; T reg40=reg9*reg0; T reg41=reg9*reg22; T reg42=reg10*reg14; T reg43=reg9*reg29;
    T reg44=reg9*reg27; T reg45=2*reg33; reg36=2*reg36; T reg46=reg10*reg18; T reg47=1-var_inter[0];
    T reg48=reg10*reg23; T reg49=reg10*reg7; T reg50=reg9*reg31; T reg51=0.5*reg2; T reg52=0.5*reg20;
    T reg53=reg32+reg28; T reg54=reg9*reg35; T reg55=reg10*reg13; T reg56=reg15*reg20; T reg57=reg15*reg19;
    T reg58=reg9*reg52; T reg59=reg29*reg37; T reg60=reg21*reg23; T reg61=reg21*reg7; T reg62=reg10*reg8;
    T reg63=reg25*reg48; T reg64=reg32+reg38; reg54=2*reg54; T reg65=2*reg41; T reg66=reg20*reg39;
    T reg67=reg15*reg7; T reg68=2*reg50; T reg69=reg21*reg25; T reg70=reg8*reg55; T reg71=reg22*reg36;
    reg43=2*reg43; T reg72=reg30+reg53; T reg73=reg6*reg24; T reg74=reg23*reg49; T reg75=reg10*reg1;
    T reg76=reg21*reg18; T reg77=reg10*reg2; reg44=2*reg44; reg47=reg47-var_inter[1]; T reg78=reg10*reg25;
    T reg79=reg9*reg51; T reg80=reg15*reg14; T reg81=reg10*reg20; T reg82=reg13*reg46; reg40=2*reg40;
    T reg83=reg11*reg45; T reg84=reg10*reg3; T reg85=reg21*reg13; T reg86=reg19*reg42; T reg87=reg0*reg43;
    T reg88=reg1*reg78; T reg89=reg21*reg8; T reg90=reg22*reg43; T reg91=reg25*reg64; T reg92=reg70+reg71;
    T reg93=reg52*reg68; T reg94=reg8*reg57; T reg95=reg15*reg23; T reg96=reg0*reg68; T reg97=reg2*reg85;
    T reg98=reg2*reg81; T reg99=reg2*reg39; T reg100=reg13*reg64; T reg101=reg15*reg25; T reg102=reg2*reg42;
    T reg103=reg2*reg77; T reg104=reg8*reg62; T reg105=reg0*reg24; T reg106=reg22*reg65; T reg107=reg8*reg69;
    T reg108=reg1*reg49; T reg109=reg0*reg37; T reg110=reg1*reg48; T reg111=reg19*reg72; T reg112=reg11*reg36;
    T reg113=reg13*reg55; T reg114=reg11*reg37; T reg115=reg13*reg60; reg82=reg83+reg82; T reg116=reg31*reg54;
    T reg117=reg13*reg80; T reg118=reg31*reg24; T reg119=reg23*reg48; T reg120=reg6*reg37; T reg121=reg14*reg42;
    T reg122=reg23*reg57; T reg123=reg31*reg36; T reg124=reg74+reg73; T reg125=reg19*reg39; T reg126=reg7*reg80;
    T reg127=reg6*reg54; T reg128=reg19*reg76; T reg129=reg35*reg45; reg86=reg73+reg86; reg73=reg18*reg46;
    T reg130=reg4*reg45; T reg131=reg7*reg49; T reg132=reg5*reg24; T reg133=reg8*reg46; T reg134=reg22*reg45;
    T reg135=reg8*reg61; T reg136=reg22*reg24; T reg137=reg29*reg43; T reg138=reg25*reg78; T reg139=reg52*reg65;
    T reg140=reg25*reg56; T reg141=reg18*reg61; T reg142=reg4*reg24; reg63=reg59+reg63; T reg143=reg29*reg45;
    T reg144=reg25*reg76; T reg145=reg29*reg24; T reg146=reg25*reg49; T reg147=reg52*reg45; T reg148=reg25*reg80;
    T reg149=reg20*reg81; T reg150=reg20*reg85; T reg151=reg29*reg68; T reg152=reg59+reg66; T reg153=reg20*reg67;
    T reg154=reg22*reg54; reg42=reg20*reg42; T reg155=var_inter[0]*(*f.m).f_vol[1]; T reg156=reg27*reg65; T reg157=var_inter[2]*(*f.m).f_vol[1];
    T reg158=reg7*reg64; T reg159=reg0*reg40; T reg160=reg15*reg2; T reg161=var_inter[1]*(*f.m).f_vol[2]; reg46=reg3*reg46;
    T reg162=var_inter[1]*(*f.m).f_vol[0]; T reg163=reg27*reg36; T reg164=reg3*reg55; T reg165=reg27*reg45; reg79=2*reg79;
    T reg166=reg84*reg3; T reg167=reg3*reg57; T reg168=reg1*reg75; T reg169=reg52*reg37; T reg170=reg51*reg37;
    T reg171=reg21*reg1; T reg172=reg27*reg44; T reg173=reg3*reg62; reg58=2*reg58; reg47=reg47-var_inter[2];
    reg42=reg145+reg42; T reg174=reg8*reg56; T reg175=reg52*reg43; T reg176=reg22*reg37; T reg177=reg8*reg60;
    T reg178=reg31*reg68; T reg179=reg112+reg113; T reg180=reg3*reg160; T reg181=reg170+reg167; T reg182=reg92+reg93;
    T reg183=reg2*reg76; T reg184=reg0*reg54; T reg185=reg13*reg61; T reg186=reg27*reg54; T reg187=reg11*reg24;
    T reg188=reg2*reg67; T reg189=reg82+reg116; T reg190=reg3*reg171; reg46=reg46-reg165; reg102=reg105+reg102;
    T reg191=reg27*reg40; T reg192=reg31*reg37; reg104=reg104+reg106; T reg193=reg52*reg58; T reg194=reg13*reg57;
    T reg195=reg91-reg155; T reg196=reg51*reg54; reg90=reg107+reg90; reg115=reg114+reg115; T reg197=reg51*reg40;
    T reg198=reg150+reg151; T reg199=reg163-reg164; T reg200=reg139+reg140; T reg201=reg18*reg80; T reg202=reg29*reg36;
    T reg203=reg25*reg85; reg149=reg137+reg149; T reg204=reg14*reg72; T reg205=reg51*reg68; T reg206=reg1*reg64;
    T reg207=reg147+reg148; T reg208=reg93+reg63; T reg209=reg3*reg69; T reg210=reg27*reg43; T reg211=reg3*reg56;
    T reg212=reg52*reg36; T reg213=reg25*reg57; reg145=reg145+reg146; T reg214=reg51*reg43; reg144=reg143+reg144;
    T reg215=reg158-reg157; T reg216=reg169+reg94; reg154=reg153+reg154; T reg217=reg8*reg64; T reg218=reg51*reg58;
    reg133=reg133+reg134; T reg219=reg52*reg54; T reg220=reg29*reg54; T reg221=reg20*reg76; T reg222=reg27*reg37;
    T reg223=reg3*reg60; reg136=reg135+reg136; T reg224=reg3*reg64; reg71=reg71+reg152; T reg225=reg8*reg80;
    T reg226=reg52*reg24; reg173=reg173-reg156; T reg227=reg22*reg68; reg137=reg137+reg138; T reg228=reg20*reg95;
    T reg229=reg2*reg72; T reg230=reg23*reg64; reg105=reg105-reg108; T reg231=reg47*(*f.m).f_vol[0]; T reg232=reg47*(*f.m).f_vol[1];
    reg73=reg73+reg130; T reg233=reg51*reg45; T reg234=reg1*reg80; T reg235=reg120+reg125; T reg236=reg47*(*f.m).f_vol[2];
    T reg237=reg35*reg54; T reg238=reg1*reg56; T reg239=reg129+reg126; T reg240=reg51*reg65; T reg241=reg31*reg45;
    reg103=reg159+reg103; T reg242=reg23*reg80; reg142=reg141+reg142; T reg243=reg2*reg89; T reg244=reg18*reg64;
    T reg245=var_inter[2]*(*f.m).f_vol[2]; reg86=reg83+reg86; reg159=reg168+reg159; reg110=reg110-reg109; reg168=reg0*reg36;
    T reg246=var_inter[1]*(*f.m).f_vol[1]; T reg247=var_inter[2]*(*f.m).f_vol[0]; T reg248=reg51*reg36; T reg249=reg1*reg57; T reg250=reg1*reg85;
    T reg251=reg19*reg67; T reg252=reg11*reg54; T reg253=reg1*reg76; T reg254=reg0*reg45; T reg255=reg51*reg44;
    reg128=reg127+reg128; reg127=reg111-reg161; T reg256=reg23*reg76; T reg257=reg35*reg24; reg123=reg122+reg123;
    T reg258=reg79*reg51; T reg259=reg97+reg96; reg119=reg119+reg120; T reg260=reg27*reg68; T reg261=reg27*reg24;
    T reg262=reg2*reg95; reg166=reg166+reg172; reg121=reg132+reg121; T reg263=reg3*reg61; T reg264=reg0*reg65;
    T reg265=reg109+reg99; T reg266=reg1*reg89; reg118=reg117+reg118; T reg267=reg20*reg72; T reg268=reg0*reg58;
    reg116=reg116+reg124; T reg269=var_inter[0]*(*f.m).f_vol[0]; T reg270=reg87-reg88; T reg271=reg100-reg162; T reg272=reg27*reg58;
    T reg273=reg2*reg101; T reg274=var_inter[0]*(*f.m).f_vol[2]; reg98=reg87+reg98; reg132=reg132+reg131; reg87=reg6*reg45;
    reg24=reg51*reg24; reg80=reg3*reg80; T reg275=reg1*reg160; reg145=reg219+reg145; reg132=reg237+reg132;
    reg149=reg106+reg149; reg237=reg73+reg237; reg73=reg26*reg142; reg173=reg218+reg173; T reg276=reg26*reg207;
    reg257=reg201+reg257; reg210=reg210-reg209; reg201=reg232+reg206; reg166=reg258+reg166; reg121=reg130+reg121;
    reg187=reg187+reg185; T reg277=reg26*reg118; T reg278=reg26*reg189; reg119=reg178+reg119; T reg279=reg194+reg192;
    T reg280=reg26*reg123; reg191=reg190+reg191; reg256=reg256+reg87; reg190=reg26*reg115; reg179=reg179+reg178;
    T reg281=reg26*reg116; T reg282=reg26*reg239; reg180=reg197+reg180; reg42=reg134+reg42; reg242=reg242+reg241;
    reg197=reg245+reg204; T reg283=reg26*reg154; reg112=reg112+reg235; reg220=reg221+reg220; reg221=reg26*reg128;
    T reg284=reg26*reg71; reg252=reg252+reg251; T reg285=reg231+reg224; reg228=reg228+reg227; T reg286=reg26*reg198;
    T reg287=reg26*reg86; reg248=reg248-reg249; T reg288=reg26*reg136; reg253=reg253-reg254; reg219=reg133+reg219;
    reg223=reg223-reg222; reg133=reg246+reg230; reg105=reg196+reg105; T reg289=reg26*reg216; T reg290=reg269+reg217;
    reg234=reg234-reg233; reg177=reg177+reg176; T reg291=reg26*reg182; T reg292=reg26*reg181; T reg293=reg247+reg244;
    reg175=reg174+reg175; reg103=reg172+reg103; reg172=reg274+reg267; reg262=reg262-reg260; reg163=reg163-reg265;
    reg184=reg183+reg184; reg174=reg26*reg259; reg261=reg261-reg263; reg186=reg186-reg188; reg98=reg98-reg156;
    reg46=reg196+reg46; reg102=reg102-reg165; reg272=reg272-reg273; reg195=reg26*reg195; reg104=reg104+reg193;
    reg80=reg24+reg80; reg24=reg26*reg90; reg268=reg243+reg268; reg271=reg26*reg271; reg183=reg26*reg144;
    reg266=reg266-reg264; reg275=reg255+reg275; reg212=reg212+reg213; reg270=reg218+reg270; reg196=reg26*reg208;
    reg211=reg214+reg211; reg127=reg26*reg127; reg238=reg238-reg240; reg202=reg202+reg203; reg168=reg168-reg250;
    reg214=reg26*reg200; reg199=reg199-reg205; reg226=reg225+reg226; reg218=reg236+reg229; reg215=reg26*reg215;
    reg110=reg110-reg205; reg137=reg193+reg137; reg159=reg258+reg159; reg193=ponderation*reg221; reg252=reg26*reg252;
    reg272=reg26*reg272; reg225=reg26*reg293; reg243=reg26*reg133; reg255=ponderation*reg73; reg258=ponderation*reg280;
    reg98=reg26*reg98; reg261=reg26*reg261; reg275=reg26*reg275; reg127=ponderation*reg127; reg119=reg26*reg119;
    reg266=reg26*reg266; T reg294=ponderation*reg174; T reg295=reg26*reg172; reg257=reg26*reg257; reg248=reg26*reg248;
    reg262=reg26*reg262; reg105=reg26*reg105; reg242=reg26*reg242; reg271=ponderation*reg271; reg132=reg26*reg132;
    reg168=reg26*reg168; reg80=reg26*reg80; reg234=reg26*reg234; reg238=reg26*reg238; reg103=reg26*reg103;
    T reg296=ponderation*reg281; reg159=reg26*reg159; T reg297=ponderation*reg287; reg253=reg26*reg253; reg268=reg26*reg268;
    reg237=reg26*reg237; reg270=reg26*reg270; reg110=reg26*reg110; reg112=reg26*reg112; T reg298=ponderation*reg282;
    reg256=reg26*reg256; reg180=reg26*reg180; reg177=reg26*reg177; T reg299=ponderation*reg283; T reg300=ponderation*reg289;
    reg223=reg26*reg223; reg220=reg26*reg220; reg219=reg26*reg219; T reg301=ponderation*reg284; T reg302=ponderation*reg288;
    T reg303=reg26*reg218; reg226=reg26*reg226; reg228=reg26*reg228; reg199=reg26*reg199; T reg304=reg26*reg285;
    T reg305=ponderation*reg286; reg137=reg26*reg137; reg173=reg26*reg173; reg149=reg26*reg149; T reg306=ponderation*reg214;
    reg202=reg26*reg202; reg215=ponderation*reg215; T reg307=ponderation*reg276; reg211=reg26*reg211; T reg308=ponderation*reg196;
    reg145=reg26*reg145; reg212=reg26*reg212; reg210=reg26*reg210; T reg309=reg26*reg201; T reg310=ponderation*reg183;
    T reg311=ponderation*reg277; reg163=reg26*reg163; reg166=reg26*reg166; reg187=reg26*reg187; reg184=reg26*reg184;
    reg46=reg26*reg46; reg121=reg26*reg121; T reg312=ponderation*reg278; reg186=reg26*reg186; reg195=ponderation*reg195;
    reg279=reg26*reg279; reg102=reg26*reg102; reg104=reg26*reg104; T reg313=ponderation*reg190; T reg314=reg26*reg290;
    T reg315=ponderation*reg24; reg191=reg26*reg191; reg42=reg26*reg42; T reg316=reg26*reg197; T reg317=ponderation*reg292;
    reg179=reg26*reg179; T reg318=ponderation*reg291; reg175=reg26*reg175; sollicitation[indices[2]+2]+=-reg127; reg127=ponderation*reg225;
    sollicitation[indices[3]+0]+=reg127; sollicitation[indices[3]+1]+=-reg215; reg215=ponderation*reg309; sollicitation[indices[0]+1]+=reg215; sollicitation[indices[2]+0]+=-reg271;
    reg271=ponderation*reg314; sollicitation[indices[1]+0]+=reg271; sollicitation[indices[1]+1]+=-reg195; reg195=ponderation*reg295; sollicitation[indices[1]+2]+=reg195;
    T tmp_10_11=-reg298; reg298=ponderation*reg304; sollicitation[indices[0]+0]+=reg298; T reg319=ponderation*reg316; sollicitation[indices[3]+2]+=reg319;
    T reg320=ponderation*reg303; sollicitation[indices[0]+2]+=reg320; T tmp_11_11=ponderation*reg121; T tmp_10_10=ponderation*reg132; reg121=ponderation*reg243;
    sollicitation[indices[2]+1]+=reg121; T tmp_1_8=ponderation*reg248; T tmp_1_9=ponderation*reg253; T tmp_4_5=-reg306; T tmp_1_10=ponderation*reg105;
    T tmp_1_11=ponderation*reg234; T tmp_2_2=ponderation*reg103; T tmp_2_3=ponderation*reg268; T tmp_2_4=ponderation*reg272; T tmp_2_5=ponderation*reg98;
    T tmp_2_6=-reg294; T tmp_2_7=ponderation*reg262; T tmp_2_8=ponderation*reg163; T tmp_2_9=ponderation*reg184; T tmp_2_10=ponderation*reg186;
    T tmp_2_11=ponderation*reg102; T tmp_3_3=ponderation*reg104; T tmp_3_4=-reg315; T tmp_3_5=ponderation*reg175; T tmp_0_0=ponderation*reg166;
    T tmp_0_1=ponderation*reg191; T tmp_0_2=ponderation*reg180; T tmp_0_3=ponderation*reg173; T tmp_0_4=ponderation*reg210; T tmp_0_5=ponderation*reg211;
    T tmp_0_6=ponderation*reg199; T tmp_0_7=ponderation*reg223; T tmp_0_8=-reg317; T tmp_0_9=ponderation*reg46; T tmp_0_10=ponderation*reg261;
    T tmp_0_11=ponderation*reg80; T tmp_1_1=ponderation*reg159; T tmp_1_2=ponderation*reg275; T tmp_1_3=ponderation*reg266; T tmp_1_4=ponderation*reg270;
    T tmp_1_5=ponderation*reg238; T tmp_1_6=ponderation*reg168; T tmp_1_7=ponderation*reg110; T tmp_5_11=ponderation*reg42; T tmp_6_6=ponderation*reg179;
    T tmp_6_7=-reg313; T tmp_6_8=ponderation*reg279; T tmp_6_9=-reg312; T tmp_6_10=ponderation*reg187; T tmp_6_11=-reg311;
    T tmp_7_7=ponderation*reg119; T tmp_7_8=-reg258; T tmp_7_9=ponderation*reg256; T tmp_7_10=-reg296; T tmp_7_11=ponderation*reg242;
    T tmp_8_8=ponderation*reg112; T tmp_8_9=-reg193; T tmp_8_10=ponderation*reg252; T tmp_8_11=-reg297; T tmp_9_9=ponderation*reg237;
    T tmp_9_10=-reg255; T tmp_9_11=ponderation*reg257; T tmp_3_6=-reg318; T tmp_3_7=ponderation*reg177; T tmp_3_8=-reg300;
    T tmp_3_9=ponderation*reg219; T tmp_3_10=-reg302; T tmp_3_11=ponderation*reg226; T tmp_4_4=ponderation*reg137; T tmp_4_6=ponderation*reg202;
    T tmp_4_7=-reg308; T tmp_4_8=ponderation*reg212; T tmp_4_9=-reg310; T tmp_4_10=ponderation*reg145; T tmp_4_11=-reg307;
    T tmp_5_5=ponderation*reg149; T tmp_5_6=-reg305; T tmp_5_7=ponderation*reg228; T tmp_5_8=-reg301; T tmp_5_9=ponderation*reg220;
    T tmp_5_10=-reg299;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[0]+2) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_3;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_4;
    matrix(indices[0]+0,indices[1]+2) += tmp_0_5;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_7;
    matrix(indices[0]+0,indices[2]+2) += tmp_0_8;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_9;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_10;
    matrix(indices[0]+0,indices[3]+2) += tmp_0_11;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[0]+2) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_3;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_4;
    matrix(indices[0]+1,indices[1]+2) += tmp_1_5;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_7;
    matrix(indices[0]+1,indices[2]+2) += tmp_1_8;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_9;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_10;
    matrix(indices[0]+1,indices[3]+2) += tmp_1_11;
    matrix(indices[0]+2,indices[0]+2) += tmp_2_2;
    matrix(indices[0]+2,indices[1]+0) += tmp_2_3;
    matrix(indices[0]+2,indices[1]+1) += tmp_2_4;
    matrix(indices[0]+2,indices[1]+2) += tmp_2_5;
    matrix(indices[0]+2,indices[2]+0) += tmp_2_6;
    matrix(indices[0]+2,indices[2]+1) += tmp_2_7;
    matrix(indices[0]+2,indices[2]+2) += tmp_2_8;
    matrix(indices[0]+2,indices[3]+0) += tmp_2_9;
    matrix(indices[0]+2,indices[3]+1) += tmp_2_10;
    matrix(indices[0]+2,indices[3]+2) += tmp_2_11;
    matrix(indices[1]+0,indices[1]+0) += tmp_3_3;
    matrix(indices[1]+0,indices[1]+1) += tmp_3_4;
    matrix(indices[1]+0,indices[1]+2) += tmp_3_5;
    matrix(indices[1]+0,indices[2]+0) += tmp_3_6;
    matrix(indices[1]+0,indices[2]+1) += tmp_3_7;
    matrix(indices[1]+0,indices[2]+2) += tmp_3_8;
    matrix(indices[1]+0,indices[3]+0) += tmp_3_9;
    matrix(indices[1]+0,indices[3]+1) += tmp_3_10;
    matrix(indices[1]+0,indices[3]+2) += tmp_3_11;
    matrix(indices[1]+1,indices[1]+1) += tmp_4_4;
    matrix(indices[1]+1,indices[1]+2) += tmp_4_5;
    matrix(indices[1]+1,indices[2]+0) += tmp_4_6;
    matrix(indices[1]+1,indices[2]+1) += tmp_4_7;
    matrix(indices[1]+1,indices[2]+2) += tmp_4_8;
    matrix(indices[1]+1,indices[3]+0) += tmp_4_9;
    matrix(indices[1]+1,indices[3]+1) += tmp_4_10;
    matrix(indices[1]+1,indices[3]+2) += tmp_4_11;
    matrix(indices[1]+2,indices[1]+2) += tmp_5_5;
    matrix(indices[1]+2,indices[2]+0) += tmp_5_6;
    matrix(indices[1]+2,indices[2]+1) += tmp_5_7;
    matrix(indices[1]+2,indices[2]+2) += tmp_5_8;
    matrix(indices[1]+2,indices[3]+0) += tmp_5_9;
    matrix(indices[1]+2,indices[3]+1) += tmp_5_10;
    matrix(indices[1]+2,indices[3]+2) += tmp_5_11;
    matrix(indices[2]+0,indices[2]+0) += tmp_6_6;
    matrix(indices[2]+0,indices[2]+1) += tmp_6_7;
    matrix(indices[2]+0,indices[2]+2) += tmp_6_8;
    matrix(indices[2]+0,indices[3]+0) += tmp_6_9;
    matrix(indices[2]+0,indices[3]+1) += tmp_6_10;
    matrix(indices[2]+0,indices[3]+2) += tmp_6_11;
    matrix(indices[2]+1,indices[2]+1) += tmp_7_7;
    matrix(indices[2]+1,indices[2]+2) += tmp_7_8;
    matrix(indices[2]+1,indices[3]+0) += tmp_7_9;
    matrix(indices[2]+1,indices[3]+1) += tmp_7_10;
    matrix(indices[2]+1,indices[3]+2) += tmp_7_11;
    matrix(indices[2]+2,indices[2]+2) += tmp_8_8;
    matrix(indices[2]+2,indices[3]+0) += tmp_8_9;
    matrix(indices[2]+2,indices[3]+1) += tmp_8_10;
    matrix(indices[2]+2,indices[3]+2) += tmp_8_11;
    matrix(indices[3]+0,indices[3]+0) += tmp_9_9;
    matrix(indices[3]+0,indices[3]+1) += tmp_9_10;
    matrix(indices[3]+0,indices[3]+2) += tmp_9_11;
    matrix(indices[3]+1,indices[3]+1) += tmp_10_10;
    matrix(indices[3]+1,indices[3]+2) += tmp_10_11;
    matrix(indices[3]+2,indices[3]+2) += tmp_11_11;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=elem.pos(3)[2]-elem.pos(0)[2]; T reg3=elem.pos(3)[1]-elem.pos(0)[1];
    T reg4=elem.pos(2)[2]-elem.pos(0)[2]; T reg5=elem.pos(2)[1]-elem.pos(0)[1]; T reg6=elem.pos(1)[2]-elem.pos(0)[2]; T reg7=elem.pos(1)[1]-elem.pos(0)[1]; T reg8=reg5*reg2;
    reg0=reg0*reg1; T reg9=reg7*reg2; T reg10=reg4*reg3; T reg11=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg12=1.0/(*f.m).elastic_modulus;
    T reg13=reg6*reg3; reg10=reg8-reg10; reg13=reg9-reg13; reg8=elem.pos(2)[0]-elem.pos(0)[0]; reg9=reg7*reg4;
    T reg14=elem.pos(1)[0]-elem.pos(0)[0]; T reg15=reg6*reg5; T reg16=reg12*reg1; reg1=reg11*reg1; T reg17=reg11*reg0;
    reg0=reg12*reg0; T reg18=reg8*reg13; T reg19=reg14*reg10; T reg20=reg12*reg16; T reg21=reg11*reg1;
    T reg22=reg12*reg0; reg16=reg11*reg16; T reg23=elem.pos(3)[0]-elem.pos(0)[0]; T reg24=reg11*reg17; reg0=reg11*reg0;
    reg15=reg9-reg15; reg18=reg19-reg18; reg9=reg7*reg23; reg22=reg22-reg24; reg19=reg6*reg23;
    reg17=reg12*reg17; T reg25=reg5*reg23; T reg26=reg14*reg2; reg0=reg24+reg0; reg1=reg12*reg1;
    reg20=reg20-reg21; T reg27=reg8*reg3; T reg28=reg4*reg23; reg2=reg8*reg2; reg16=reg21+reg16;
    reg23=reg23*reg15; reg3=reg14*reg3; reg28=reg2-reg28; reg9=reg3-reg9; reg4=reg14*reg4;
    reg25=reg27-reg25; reg6=reg6*reg8; reg19=reg26-reg19; reg5=reg14*reg5; reg8=reg7*reg8;
    reg23=reg18+reg23; reg2=reg21+reg1; reg3=reg11*reg0; reg7=reg12*reg22; reg16=reg11*reg16;
    reg20=reg12*reg20; reg17=reg24+reg17; reg8=reg5-reg8; reg6=reg4-reg6; reg4=reg11*reg17;
    reg9=reg9/reg23; reg3=reg7-reg3; reg19=reg19/reg23; reg13=reg13/reg23; reg25=reg25/reg23;
    reg28=reg28/reg23; reg16=reg20-reg16; reg10=reg10/reg23; reg2=reg11*reg2; reg4=reg3-reg4;
    reg3=reg28-reg19; reg5=reg13-reg10; reg15=reg15/reg23; reg6=reg6/reg23; reg2=reg16-reg2;
    reg7=reg9-reg25; reg8=reg8/reg23; reg11=0.5*reg15; reg12=0.5*reg6; reg7=reg7-reg8;
    reg3=reg6+reg3; reg14=0.5*reg28; reg5=reg5-reg15; reg16=0.5*reg10; reg18=0.5*reg13;
    reg20=0.5*reg19; reg2=reg2/reg4; reg24=0.5*reg9; reg26=reg2*reg16; reg27=0.5*reg3;
    T reg29=reg2*reg11; T reg30=reg2*reg12; T reg31=reg2*reg14; T reg32=0.5*reg25; T reg33=0.5*reg5;
    T reg34=0.5*reg7; T reg35=reg2*reg18; T reg36=0.5*reg8; T reg37=reg2*reg20; reg22=reg22/reg4;
    T reg38=reg22*reg9; T reg39=reg2*reg32; T reg40=reg22*reg8; reg29=2*reg29; T reg41=reg22*reg10;
    T reg42=reg22*reg13; T reg43=reg2*reg24; T reg44=2*reg31; T reg45=2*reg30; T reg46=reg22*reg15;
    T reg47=reg2*reg36; T reg48=reg22*reg6; reg26=2*reg26; T reg49=2*reg35; reg37=2*reg37;
    T reg50=reg2*reg27; T reg51=reg22*reg28; T reg52=reg2*reg33; T reg53=reg22*reg25; T reg54=reg22*reg19;
    T reg55=reg2*reg34; reg0=reg0/reg4; reg4=reg17/reg4; reg17=reg0*reg19; T reg56=reg9*reg40;
    T reg57=reg9*reg53; T reg58=reg10*reg42; T reg59=reg0*reg15; T reg60=reg4*reg9; T reg61=reg20*reg44;
    T reg62=reg6*reg54; T reg63=reg11*reg49; T reg64=reg0*reg13; T reg65=reg13*reg41; T reg66=reg22*reg7;
    T reg67=reg22*reg3; T reg68=reg4*reg19; T reg69=reg19*reg51; T reg70=reg0*reg5; T reg71=reg13*reg46;
    T reg72=reg4*reg8; T reg73=reg4*reg6; T reg74=reg20*reg45; T reg75=reg18*reg26; T reg76=reg4*reg28;
    T reg77=reg0*reg10; T reg78=reg12*reg37; T reg79=reg15*reg42; T reg80=reg0*reg6; T reg81=reg19*reg48;
    T reg82=reg18*reg29; reg47=2*reg47; reg52=2*reg52; T reg83=reg28*reg54; reg39=2*reg39;
    reg50=2*reg50; T reg84=reg0*reg28; T reg85=reg4*reg25; T reg86=reg16*reg49; T reg87=reg25*reg38;
    T reg88=reg8*reg38; reg55=2*reg55; T reg89=reg22*reg5; T reg90=reg4*reg7; T reg91=2*reg43;
    T reg92=reg14*reg37; T reg93=reg18*reg49; T reg94=reg12*reg29; T reg95=reg15*reg80; T reg96=reg19*reg54;
    T reg97=reg19*reg60; T reg98=reg18*reg37; T reg99=reg19*reg64; T reg100=reg7*reg53; T reg101=reg69+reg75;
    reg83=reg86+reg83; T reg102=reg7*reg64; T reg103=reg33*reg91; T reg104=reg28*reg67; T reg105=reg16*reg52;
    T reg106=reg18*reg52; T reg107=reg19*reg67; T reg108=reg16*reg45; T reg109=reg9*reg59; T reg110=reg16*reg26;
    T reg111=reg18*reg47; T reg112=reg9*reg38; T reg113=reg9*reg68; T reg114=reg20*reg91; T reg115=reg3*reg48;
    T reg116=reg33*reg29; reg57=reg75+reg57; reg75=reg28*reg51; T reg117=reg32*reg44; T reg118=reg9*reg77;
    T reg119=reg18*reg39; T reg120=reg9*reg66; T reg121=reg9*reg70; T reg122=reg18*reg55; T reg123=reg4*reg3;
    T reg124=reg81+reg82; T reg125=reg7*reg66; T reg126=reg28*reg85; T reg127=reg28*reg77; T reg128=reg16*reg44;
    T reg129=reg24*reg37; T reg130=reg24*reg39; T reg131=reg25*reg76; reg65=reg61+reg65; T reg132=reg14*reg39;
    T reg133=reg12*reg50; T reg134=reg10*reg60; T reg135=reg25*reg53; T reg136=reg25*reg64; T reg137=reg24*reg52;
    T reg138=reg13*reg90; T reg139=reg10*reg41; T reg140=reg14*reg44; T reg141=reg89*reg13; T reg142=reg16*reg91;
    T reg143=reg10*reg84; T reg144=reg14*reg26; T reg145=reg20*reg50; T reg146=reg25*reg40; T reg147=reg14*reg47;
    T reg148=reg25*reg73; T reg149=reg86+reg87; T reg150=reg32*reg91; T reg151=reg58+reg92; T reg152=reg7*reg38;
    T reg153=reg24*reg29; T reg154=reg13*reg72; T reg155=reg24*reg47; reg71=reg74+reg71; T reg156=reg28*reg59;
    T reg157=reg16*reg29; T reg158=reg28*reg48; T reg159=reg32*reg45; T reg160=reg7*reg40; T reg161=reg13*reg17;
    T reg162=reg28*reg72; T reg163=reg14*reg29; T reg164=reg20*reg49; T reg165=reg10*reg80; T reg166=reg89*reg10;
    T reg167=reg14*reg50; T reg168=reg13*reg42; T reg169=reg20*reg37; T reg170=reg24*reg26; T reg171=reg25*reg66;
    T reg172=reg14*reg45; T reg173=reg13*reg85; T reg174=reg10*reg46; T reg175=reg36*reg45; T reg176=reg6*reg72;
    T reg177=reg32*reg49; reg66=reg8*reg66; T reg178=reg3*reg67; T reg179=reg33*reg52; T reg180=reg12*reg45;
    T reg181=reg15*reg46; T reg182=reg8*reg76; T reg183=reg12*reg39; T reg184=reg36*reg49; T reg185=reg15*reg60;
    T reg186=reg36*reg91; T reg187=reg79+reg78; T reg188=reg27*reg44; T reg189=reg5*reg41; reg53=reg8*reg53;
    T reg190=reg3*reg51; reg62=reg63+reg62; T reg191=reg34*reg49; T reg192=reg5*reg60; T reg193=reg11*reg45;
    reg46=reg5*reg46; T reg194=reg27*reg45; T reg195=reg6*reg59; T reg196=reg6*reg85; T reg197=reg36*reg44;
    T reg198=reg6*reg51; T reg199=reg11*reg29; T reg200=reg11*reg26; T reg201=reg6*reg77; T reg202=reg11*reg44;
    T reg203=reg6*reg48; T reg204=reg27*reg37; T reg205=reg5*reg42; reg67=reg6*reg67; T reg206=reg11*reg52;
    T reg207=reg63+reg88; T reg208=reg0*reg3; T reg209=reg8*reg73; T reg210=reg12*reg47; T reg211=reg89*reg15;
    reg41=reg15*reg41; reg54=reg3*reg54; T reg212=reg33*reg49; T reg213=reg12*reg44; T reg214=reg27*reg50;
    reg89=reg89*reg5; T reg215=reg11*reg91; reg56=reg82+reg56; reg82=reg33*reg26; T reg216=reg8*reg64;
    T reg217=reg15*reg84; T reg218=reg12*reg26; reg40=reg8*reg40; T reg219=reg199+reg203; T reg220=reg8*reg59;
    reg40=reg199+reg40; reg199=reg28*reg64; reg92=reg92+reg149; T reg221=reg25*reg59; T reg222=reg16*reg47;
    T reg223=reg16*reg37; T reg224=reg11*reg47; T reg225=reg15*reg72; reg147=reg148+reg147; T reg226=reg117+reg126;
    T reg227=reg36*reg37; reg146=reg157+reg146; T reg228=reg110+reg75; reg195=reg193+reg195; reg210=reg209+reg210;
    T reg229=reg13*reg208; T reg230=reg20*reg52; T reg231=reg55*reg24; reg141=reg145-reg141; T reg232=reg6*reg60;
    reg53=reg200+reg53; T reg233=reg159+reg162; T reg234=reg25*reg70; T reg235=reg55*reg16; reg157=reg157+reg158;
    T reg236=reg25*reg123; T reg237=reg55*reg14; reg183=reg182+reg183; reg171=reg105+reg171; T reg238=reg25*reg77;
    reg156=reg108+reg156; T reg239=reg216+reg215; T reg240=reg16*reg39; T reg241=reg8*reg68; T reg242=reg28*reg60;
    T reg243=reg32*reg37; T reg244=reg11*reg39; T reg245=reg8*reg77; reg132=reg131+reg132; reg66=reg206+reg66;
    T reg246=reg150+reg83; T reg247=reg12*reg55; reg135=reg110+reg135; reg110=reg8*reg123; T reg248=reg11*reg55;
    T reg249=reg136+reg142; T reg250=reg25*reg68; T reg251=reg14*reg91; T reg252=reg8*reg70; T reg253=reg12*reg91;
    T reg254=reg175+reg176; reg78=reg78+reg207; T reg255=reg55*reg20; reg121=reg122+reg121; reg41=reg41+reg213;
    reg122=reg36*reg39; T reg256=reg24*reg45; T reg257=reg19*reg72; T reg258=reg155+reg124; reg218=reg217+reg218;
    T reg259=reg18*reg45; T reg260=reg19*reg59; reg129=reg97+reg129; T reg261=reg15*reg85; T reg262=reg36*reg26;
    reg96=reg96+reg93; T reg263=reg187+reg186; reg98=reg99+reg98; T reg264=reg15*reg17; T reg265=reg12*reg49;
    T reg266=reg24*reg44; T reg267=reg20*reg47; reg109=reg111+reg109; reg111=reg9*reg73; T reg268=reg93+reg112;
    reg56=reg74+reg56; reg113=reg114+reg113; T reg269=reg9*reg64; T reg270=reg18*reg91; reg211=reg211-reg133;
    T reg271=reg36*reg55; reg57=reg61+reg57; T reg272=reg15*reg208; T reg273=reg12*reg52; T reg274=reg9*reg76;
    T reg275=reg20*reg39; reg118=reg119+reg118; reg119=reg15*reg90; T reg276=reg36*reg52; reg120=reg106+reg120;
    T reg277=reg9*reg123; reg67=reg206-reg67; reg206=reg36*reg50; T reg278=reg6*reg90; T reg279=reg24*reg49;
    T reg280=reg13*reg60; reg161=reg164+reg161; reg201=reg202+reg201; T reg281=reg24*reg91; T reg282=reg169+reg168;
    reg200=reg200+reg198; reg170=reg173+reg170; T reg283=reg13*reg84; T reg284=reg20*reg26; T reg285=reg65+reg130;
    T reg286=reg197+reg196; T reg287=reg11*reg37; T reg288=reg6*reg64; reg137=reg138+reg137; T reg289=reg186+reg62;
    T reg290=reg19*reg85; reg130=reg130+reg101; T reg291=reg185+reg184; T reg292=reg18*reg44; T reg293=reg19*reg77;
    T reg294=reg24*reg50; T reg295=reg19*reg90; reg181=reg181+reg180; reg106=reg107-reg106; reg107=reg36*reg47;
    reg94=reg95+reg94; T reg296=reg18*reg50; T reg297=reg36*reg29; T reg298=reg11*reg50; T reg299=reg19*reg70;
    reg153=reg154+reg153; T reg300=reg6*reg70; T reg301=reg13*reg80; T reg302=reg20*reg29; reg155=reg71+reg155;
    T reg303=reg204-reg205; T reg304=reg34*reg39; T reg305=reg14*reg52; reg54=reg54-reg212; T reg306=reg10*reg208;
    reg174=reg174+reg172; reg100=reg82+reg100; T reg307=reg34*reg91; T reg308=reg32*reg47; T reg309=reg33*reg37;
    reg208=reg5*reg208; T reg310=reg3*reg64; T reg311=reg32*reg26; T reg312=reg191+reg192; T reg313=reg102+reg103;
    T reg314=reg10*reg85; T reg315=reg7*reg70; T reg316=reg3*reg85; reg163=reg165+reg163; T reg317=reg34*reg44;
    reg139=reg139+reg140; T reg318=reg32*reg39; T reg319=reg7*reg77; T reg320=reg34*reg45; T reg321=reg151+reg150;
    T reg322=reg33*reg39; T reg323=reg116-reg115; reg85=reg5*reg85; T reg324=reg34*reg52; T reg325=reg34*reg26;
    T reg326=reg5*reg17; T reg327=reg32*reg52; T reg328=reg33*reg45; reg17=reg10*reg17; T reg329=reg14*reg49;
    T reg330=reg10*reg90; T reg331=reg3*reg59; reg39=reg27*reg39; T reg332=reg3*reg72; T reg333=reg3*reg60;
    reg52=reg27*reg52; T reg334=reg7*reg76; T reg335=reg27*reg49; T reg336=reg177+reg134; reg37=reg34*reg37;
    reg178=reg178+reg179; T reg337=reg55*reg27; reg104=reg105-reg104; reg105=reg212+reg152; T reg338=reg5*reg90;
    T reg339=reg55*reg34; T reg340=reg33*reg50; reg144=reg143+reg144; T reg341=reg5*reg84; T reg342=reg32*reg50;
    T reg343=reg28*reg90; T reg344=reg3*reg70; reg123=reg7*reg123; reg59=reg7*reg59; T reg345=reg5*reg72;
    T reg346=reg7*reg73; reg127=reg128+reg127; T reg347=reg5*reg80; T reg348=reg34*reg29; T reg349=reg27*reg47;
    T reg350=reg33*reg47; T reg351=reg27*reg29; reg82=reg82-reg190; T reg352=reg55*reg32; reg166=reg166-reg167;
    reg72=reg10*reg72; reg189=reg189-reg188; T reg353=reg27*reg91; T reg354=reg33*reg44; reg29=reg32*reg29;
    reg77=reg3*reg77; reg55=reg55*reg33; reg90=reg3*reg90; reg68=reg7*reg68; reg89=reg89+reg214;
    T reg355=reg16*reg50; reg70=reg28*reg70; reg125=reg179+reg125; reg50=reg34*reg50; reg47=reg34*reg47;
    reg160=reg116+reg160; reg46=reg46-reg194; reg26=reg27*reg26; reg125=reg214+reg125; reg260=reg260+reg259;
    reg116=reg23*reg254; reg85=reg325+reg85; reg179=reg23*reg129; reg214=reg23*reg258; reg325=reg23*reg155;
    reg302=reg302+reg301; reg204=reg204-reg105; T reg356=reg23*reg153; reg133=reg66-reg133; reg296=reg299-reg296;
    reg68=reg68-reg353; reg106=reg106-reg231; reg26=reg26-reg341; reg294=reg295-reg294; reg66=reg23*reg313;
    reg293=reg293+reg292; reg247=reg110-reg247; reg100=reg100-reg188; reg110=reg23*reg130; reg290=reg290+reg266;
    reg39=reg39-reg334; reg295=reg23*reg98; reg248=reg252+reg248; reg96=reg281+reg96; reg322=reg319+reg322;
    reg273=reg272-reg273; reg276=reg119+reg276; reg119=reg23*reg312; reg316=reg316-reg317; reg41=reg41+reg122;
    reg82=reg304+reg82; reg252=reg23*reg218; reg287=reg287+reg288; reg262=reg261+reg262; reg77=reg77-reg354;
    reg261=reg23*reg263; reg90=reg50+reg90; reg264=reg264+reg265; reg50=reg23*reg286; reg272=reg23*reg291;
    reg178=reg339+reg178; reg46=reg47+reg46; reg181=reg181+reg107; reg340=reg344+reg340; reg299=reg23*reg94;
    reg297=reg225+reg297; reg200=reg122+reg200; reg300=reg298-reg300; reg345=reg348+reg345; reg67=reg271+reg67;
    reg206=reg206-reg278; reg351=reg351-reg347; reg122=reg23*reg201; reg123=reg337+reg123; reg257=reg257+reg256;
    reg225=reg23*reg121; reg55=reg315+reg55; reg277=reg255-reg277; reg219=reg107+reg219; reg120=reg145-reg120;
    reg332=reg332-reg320; reg107=reg23*reg118; reg303=reg303-reg307; reg275=reg275+reg274; reg323=reg47+reg323;
    reg47=reg23*reg57; reg145=reg23*reg195; reg255=reg270+reg269; reg331=reg331-reg328; reg298=reg23*reg113;
    reg169=reg169+reg268; reg227=reg227+reg232; reg37=reg37-reg333; reg315=reg23*reg109; reg326=reg326-reg335;
    reg267=reg267+reg111; reg54=reg54-reg307; reg319=reg23*reg56; reg337=reg23*reg289; reg309=reg309-reg310;
    reg271=reg211+reg271; reg222=reg221+reg222; reg167=reg171-reg167; reg311=reg314+reg311; reg53=reg213+reg53;
    reg89=reg339+reg89; reg171=reg23*reg156; reg305=reg306-reg305; reg211=reg23*reg285; reg221=reg23*reg210;
    reg243=reg243+reg242; reg284=reg284+reg283; reg306=reg23*reg92; reg70=reg355-reg70; reg314=reg23*reg321;
    reg240=reg238+reg240; reg238=reg23*reg246; reg166=reg166+reg352; reg146=reg172+reg146; reg235=reg234+reg235;
    reg234=reg23*reg163; reg224=reg220+reg224; reg338=reg324+reg338; reg220=reg23*reg144; reg139=reg139+reg318;
    reg231=reg141-reg231; reg141=reg23*reg233; reg324=reg23*reg239; reg229=reg230-reg229; reg237=reg236-reg237;
    reg230=reg23*reg147; reg174=reg174+reg308; reg157=reg308+reg157; reg29=reg72+reg29; reg327=reg330+reg327;
    reg72=reg23*reg137; reg189=reg304+reg189; reg160=reg160-reg194; reg236=reg23*reg336; reg304=reg23*reg132;
    reg308=reg23*reg161; reg52=reg208+reg52; reg208=reg23*reg226; reg40=reg180+reg40; reg349=reg349-reg346;
    reg330=reg23*reg127; reg339=reg280+reg279; reg344=reg23*reg249; reg17=reg17+reg329; reg135=reg140+reg135;
    reg244=reg245+reg244; reg228=reg318+reg228; reg350=reg59+reg350; reg104=reg352+reg104; reg59=reg23*reg170;
    reg241=reg241+reg253; reg250=reg250+reg251; reg245=reg23*reg78; reg318=reg23*reg183; reg223=reg223+reg199;
    reg342=reg342-reg343; reg282=reg282+reg281; reg135=reg23*reg135; reg54=reg23*reg54; reg169=reg23*reg169;
    reg240=reg23*reg240; reg167=reg23*reg167; reg276=reg23*reg276; reg316=reg23*reg316; reg174=reg23*reg174;
    reg267=reg23*reg267; reg348=ponderation*reg245; reg37=reg23*reg37; reg352=ponderation*reg319; reg273=reg23*reg273;
    reg17=reg23*reg17; reg237=reg23*reg237; reg309=reg23*reg309; reg326=reg23*reg326; reg355=ponderation*reg236;
    T reg357=ponderation*reg315; reg271=reg23*reg271; T reg358=ponderation*reg304; T reg359=ponderation*reg337; reg70=reg23*reg70;
    T reg360=ponderation*reg221; T reg361=ponderation*reg238; reg340=reg23*reg340; reg181=reg23*reg181; reg104=reg23*reg104;
    T reg362=ponderation*reg299; reg200=reg23*reg200; reg223=reg23*reg223; reg297=reg23*reg297; reg46=reg23*reg46;
    reg342=reg23*reg342; reg300=reg23*reg300; reg345=reg23*reg345; T reg363=ponderation*reg208; reg67=reg23*reg67;
    reg40=reg23*reg40; reg228=reg23*reg228; reg206=reg23*reg206; T reg364=ponderation*reg122; reg351=reg23*reg351;
    T reg365=ponderation*reg330; reg235=reg23*reg235; reg224=reg23*reg224; reg41=reg23*reg41; reg287=reg23*reg287;
    reg82=reg23*reg82; T reg366=ponderation*reg141; T reg367=ponderation*reg234; T reg368=ponderation*reg252; T reg369=ponderation*reg119;
    reg77=reg23*reg77; reg89=reg23*reg89; reg262=reg23*reg262; reg157=reg23*reg157; reg90=reg23*reg90;
    reg29=reg23*reg29; T reg370=ponderation*reg261; T reg371=ponderation*reg50; T reg372=ponderation*reg171; reg264=reg23*reg264;
    reg243=reg23*reg243; reg178=reg23*reg178; T reg373=ponderation*reg272; reg322=reg23*reg322; T reg374=ponderation*reg72;
    reg96=reg23*reg96; reg327=reg23*reg327; reg189=reg23*reg189; reg53=reg23*reg53; reg296=reg23*reg296;
    T reg375=ponderation*reg179; T reg376=ponderation*reg308; T reg377=ponderation*reg116; reg125=reg23*reg125; reg229=reg23*reg229;
    reg260=reg23*reg260; reg338=reg23*reg338; reg231=reg23*reg231; reg123=reg23*reg123; T reg378=ponderation*reg214;
    reg85=reg23*reg85; reg139=reg23*reg139; reg146=reg23*reg146; reg247=reg23*reg247; T reg379=ponderation*reg66;
    T reg380=ponderation*reg59; reg293=reg23*reg293; reg294=reg23*reg294; reg166=reg23*reg166; reg26=reg23*reg26;
    reg100=reg23*reg100; T reg381=ponderation*reg318; reg160=reg23*reg160; T reg382=ponderation*reg110; reg106=reg23*reg106;
    reg284=reg23*reg284; reg282=reg23*reg282; reg290=reg23*reg290; reg248=reg23*reg248; reg39=reg23*reg39;
    T reg383=ponderation*reg211; T reg384=ponderation*reg295; reg305=reg23*reg305; reg68=reg23*reg68; reg302=reg23*reg302;
    T reg385=ponderation*reg107; reg339=reg23*reg339; reg244=reg23*reg244; reg275=reg23*reg275; T reg386=ponderation*reg145;
    reg323=reg23*reg323; reg250=reg23*reg250; T reg387=ponderation*reg47; T reg388=ponderation*reg314; reg303=reg23*reg303;
    T reg389=ponderation*reg325; reg241=reg23*reg241; T reg390=reg23*reg255; reg350=reg23*reg350; reg331=reg23*reg331;
    T reg391=ponderation*reg344; T reg392=ponderation*reg298; reg52=reg23*reg52; reg227=reg23*reg227; reg133=reg23*reg133;
    reg257=reg23*reg257; T reg393=ponderation*reg356; reg349=reg23*reg349; reg55=reg23*reg55; T reg394=ponderation*reg225;
    reg204=reg23*reg204; reg219=reg23*reg219; T reg395=ponderation*reg220; T reg396=ponderation*reg230; reg277=reg23*reg277;
    T reg397=ponderation*reg324; reg222=reg23*reg222; reg120=reg23*reg120; reg332=reg23*reg332; reg311=reg23*reg311;
    T reg398=ponderation*reg306; T tmp_11_2=ponderation*reg133; T tmp_0_3=ponderation*reg189; T tmp_11_3=ponderation*reg244; T tmp_0_9=ponderation*reg46;
    T tmp_10_3=-reg364; T tmp_11_1=ponderation*reg247; T tmp_11_11=ponderation*reg40; T tmp_10_8=ponderation*reg227; T tmp_11_7=ponderation*reg241;
    T tmp_0_6=ponderation*reg303; T tmp_10_7=-reg359; T tmp_10_9=-reg386; T tmp_0_1=ponderation*reg52; T tmp_11_8=-reg348;
    T tmp_11_6=-reg397; T tmp_0_7=ponderation*reg326; T tmp_0_0=ponderation*reg89; T tmp_10_10=ponderation*reg219; T tmp_0_5=ponderation*reg85;
    T tmp_10_6=ponderation*reg287; T tmp_11_9=ponderation*reg224; T tmp_10_11=-reg377; T tmp_0_8=-reg369; T tmp_11_5=ponderation*reg53;
    T tmp_10_5=-reg371; T tmp_0_2=ponderation*reg338; T tmp_11_0=ponderation*reg248; T tmp_0_4=ponderation*reg26; T tmp_11_4=-reg381;
    T tmp_11_10=-reg360; T tmp_10_4=ponderation*reg200; T tmp_3_4=-reg395; T tmp_5_10=-reg396; T tmp_5_11=ponderation*reg146;
    T tmp_3_3=ponderation*reg139; T tmp_6_0=ponderation*reg231; T tmp_6_1=ponderation*reg229; T tmp_3_2=ponderation*reg327; T tmp_6_2=-reg374;
    T tmp_3_1=ponderation*reg305; T tmp_6_3=-reg383; T tmp_6_4=ponderation*reg284; T tmp_3_0=ponderation*reg166; T tmp_6_5=-reg380;
    T tmp_2_11=ponderation*reg160; T tmp_6_6=ponderation*reg282; T tmp_2_10=ponderation*reg349; T tmp_6_7=-reg376; T tmp_6_8=ponderation*reg339;
    T tmp_2_9=ponderation*reg350; T tmp_6_9=-reg389; T tmp_6_10=ponderation*reg302; T tmp_2_8=ponderation*reg204; T tmp_6_11=-reg393;
    T tmp_2_7=ponderation*reg68; T tmp_7_0=ponderation*reg296; T tmp_7_1=ponderation*reg106; T tmp_2_6=-reg379; T tmp_7_2=ponderation*reg294;
    T tmp_4_3=-reg365; T tmp_4_4=ponderation*reg228; T tmp_4_2=ponderation*reg342; T tmp_4_1=ponderation*reg104; T tmp_4_6=ponderation*reg223;
    T tmp_4_0=ponderation*reg70; T tmp_4_7=-reg361; T tmp_4_8=ponderation*reg243; T tmp_3_11=ponderation*reg29; T tmp_4_9=-reg372;
    T tmp_4_10=ponderation*reg157; T tmp_3_10=-reg367; T tmp_4_11=-reg366; T tmp_5_0=ponderation*reg235; T tmp_3_9=ponderation*reg174;
    T tmp_5_1=ponderation*reg237; T tmp_5_2=ponderation*reg167; T tmp_3_8=-reg355; T tmp_5_3=ponderation*reg240; T tmp_5_4=-reg358;
    T tmp_3_7=ponderation*reg17; T tmp_5_5=ponderation*reg135; T tmp_5_6=-reg391; T tmp_3_6=-reg388; T tmp_5_7=ponderation*reg250;
    T tmp_3_5=ponderation*reg311; T tmp_5_8=-reg398; T tmp_5_9=ponderation*reg222; T tmp_8_9=-reg357; T tmp_1_7=ponderation*reg54;
    T tmp_8_10=ponderation*reg267; T tmp_8_11=-reg352; T tmp_1_6=ponderation*reg309; T tmp_9_0=ponderation*reg271; T tmp_9_1=ponderation*reg273;
    T tmp_1_5=ponderation*reg316; T tmp_9_2=ponderation*reg276; T tmp_1_4=ponderation*reg82; T tmp_9_3=ponderation*reg41; T tmp_9_4=-reg368;
    T tmp_1_3=ponderation*reg77; T tmp_9_5=ponderation*reg262; T tmp_1_2=ponderation*reg90; T tmp_9_6=-reg370; T tmp_9_7=ponderation*reg264;
    T tmp_1_1=ponderation*reg178; T tmp_9_8=-reg373; T tmp_1_0=ponderation*reg340; T tmp_9_9=ponderation*reg181; T tmp_9_10=-reg362;
    T tmp_9_11=ponderation*reg297; T tmp_0_11=ponderation*reg345; T tmp_10_0=ponderation*reg300; T tmp_10_1=ponderation*reg67; T tmp_0_10=ponderation*reg351;
    T tmp_10_2=ponderation*reg206; T tmp_7_3=ponderation*reg293; T tmp_2_5=ponderation*reg100; T tmp_7_4=-reg382; T tmp_2_4=ponderation*reg39;
    T tmp_7_5=ponderation*reg290; T tmp_7_6=-reg384; T tmp_2_3=ponderation*reg322; T tmp_7_7=ponderation*reg96; T tmp_7_8=-reg375;
    T tmp_2_2=ponderation*reg125; T tmp_7_9=ponderation*reg260; T tmp_2_1=ponderation*reg123; T tmp_7_10=-reg378; T tmp_7_11=ponderation*reg257;
    T tmp_2_0=ponderation*reg55; T tmp_8_0=-reg394; T tmp_8_1=ponderation*reg277; T tmp_1_11=ponderation*reg332; T tmp_8_2=ponderation*reg120;
    T tmp_8_3=-reg385; T tmp_1_10=ponderation*reg323; T tmp_8_4=ponderation*reg275; T tmp_8_5=-reg387; T tmp_4_5=-reg363;
    T tmp_1_9=ponderation*reg331; T tmp_8_6=ponderation*reg390; T tmp_8_7=-reg392; T tmp_1_8=ponderation*reg37; T tmp_8_8=ponderation*reg169;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[0]+2) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_3;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_4;
    matrix(indices[0]+0,indices[1]+2) += tmp_0_5;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_7;
    matrix(indices[0]+0,indices[2]+2) += tmp_0_8;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_9;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_10;
    matrix(indices[0]+0,indices[3]+2) += tmp_0_11;
    matrix(indices[0]+1,indices[0]+0) += tmp_1_0;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[0]+2) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_3;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_4;
    matrix(indices[0]+1,indices[1]+2) += tmp_1_5;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_7;
    matrix(indices[0]+1,indices[2]+2) += tmp_1_8;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_9;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_10;
    matrix(indices[0]+1,indices[3]+2) += tmp_1_11;
    matrix(indices[0]+2,indices[0]+0) += tmp_2_0;
    matrix(indices[0]+2,indices[0]+1) += tmp_2_1;
    matrix(indices[0]+2,indices[0]+2) += tmp_2_2;
    matrix(indices[0]+2,indices[1]+0) += tmp_2_3;
    matrix(indices[0]+2,indices[1]+1) += tmp_2_4;
    matrix(indices[0]+2,indices[1]+2) += tmp_2_5;
    matrix(indices[0]+2,indices[2]+0) += tmp_2_6;
    matrix(indices[0]+2,indices[2]+1) += tmp_2_7;
    matrix(indices[0]+2,indices[2]+2) += tmp_2_8;
    matrix(indices[0]+2,indices[3]+0) += tmp_2_9;
    matrix(indices[0]+2,indices[3]+1) += tmp_2_10;
    matrix(indices[0]+2,indices[3]+2) += tmp_2_11;
    matrix(indices[1]+0,indices[0]+0) += tmp_3_0;
    matrix(indices[1]+0,indices[0]+1) += tmp_3_1;
    matrix(indices[1]+0,indices[0]+2) += tmp_3_2;
    matrix(indices[1]+0,indices[1]+0) += tmp_3_3;
    matrix(indices[1]+0,indices[1]+1) += tmp_3_4;
    matrix(indices[1]+0,indices[1]+2) += tmp_3_5;
    matrix(indices[1]+0,indices[2]+0) += tmp_3_6;
    matrix(indices[1]+0,indices[2]+1) += tmp_3_7;
    matrix(indices[1]+0,indices[2]+2) += tmp_3_8;
    matrix(indices[1]+0,indices[3]+0) += tmp_3_9;
    matrix(indices[1]+0,indices[3]+1) += tmp_3_10;
    matrix(indices[1]+0,indices[3]+2) += tmp_3_11;
    matrix(indices[1]+1,indices[0]+0) += tmp_4_0;
    matrix(indices[1]+1,indices[0]+1) += tmp_4_1;
    matrix(indices[1]+1,indices[0]+2) += tmp_4_2;
    matrix(indices[1]+1,indices[1]+0) += tmp_4_3;
    matrix(indices[1]+1,indices[1]+1) += tmp_4_4;
    matrix(indices[1]+1,indices[1]+2) += tmp_4_5;
    matrix(indices[1]+1,indices[2]+0) += tmp_4_6;
    matrix(indices[1]+1,indices[2]+1) += tmp_4_7;
    matrix(indices[1]+1,indices[2]+2) += tmp_4_8;
    matrix(indices[1]+1,indices[3]+0) += tmp_4_9;
    matrix(indices[1]+1,indices[3]+1) += tmp_4_10;
    matrix(indices[1]+1,indices[3]+2) += tmp_4_11;
    matrix(indices[1]+2,indices[0]+0) += tmp_5_0;
    matrix(indices[1]+2,indices[0]+1) += tmp_5_1;
    matrix(indices[1]+2,indices[0]+2) += tmp_5_2;
    matrix(indices[1]+2,indices[1]+0) += tmp_5_3;
    matrix(indices[1]+2,indices[1]+1) += tmp_5_4;
    matrix(indices[1]+2,indices[1]+2) += tmp_5_5;
    matrix(indices[1]+2,indices[2]+0) += tmp_5_6;
    matrix(indices[1]+2,indices[2]+1) += tmp_5_7;
    matrix(indices[1]+2,indices[2]+2) += tmp_5_8;
    matrix(indices[1]+2,indices[3]+0) += tmp_5_9;
    matrix(indices[1]+2,indices[3]+1) += tmp_5_10;
    matrix(indices[1]+2,indices[3]+2) += tmp_5_11;
    matrix(indices[2]+0,indices[0]+0) += tmp_6_0;
    matrix(indices[2]+0,indices[0]+1) += tmp_6_1;
    matrix(indices[2]+0,indices[0]+2) += tmp_6_2;
    matrix(indices[2]+0,indices[1]+0) += tmp_6_3;
    matrix(indices[2]+0,indices[1]+1) += tmp_6_4;
    matrix(indices[2]+0,indices[1]+2) += tmp_6_5;
    matrix(indices[2]+0,indices[2]+0) += tmp_6_6;
    matrix(indices[2]+0,indices[2]+1) += tmp_6_7;
    matrix(indices[2]+0,indices[2]+2) += tmp_6_8;
    matrix(indices[2]+0,indices[3]+0) += tmp_6_9;
    matrix(indices[2]+0,indices[3]+1) += tmp_6_10;
    matrix(indices[2]+0,indices[3]+2) += tmp_6_11;
    matrix(indices[2]+1,indices[0]+0) += tmp_7_0;
    matrix(indices[2]+1,indices[0]+1) += tmp_7_1;
    matrix(indices[2]+1,indices[0]+2) += tmp_7_2;
    matrix(indices[2]+1,indices[1]+0) += tmp_7_3;
    matrix(indices[2]+1,indices[1]+1) += tmp_7_4;
    matrix(indices[2]+1,indices[1]+2) += tmp_7_5;
    matrix(indices[2]+1,indices[2]+0) += tmp_7_6;
    matrix(indices[2]+1,indices[2]+1) += tmp_7_7;
    matrix(indices[2]+1,indices[2]+2) += tmp_7_8;
    matrix(indices[2]+1,indices[3]+0) += tmp_7_9;
    matrix(indices[2]+1,indices[3]+1) += tmp_7_10;
    matrix(indices[2]+1,indices[3]+2) += tmp_7_11;
    matrix(indices[2]+2,indices[0]+0) += tmp_8_0;
    matrix(indices[2]+2,indices[0]+1) += tmp_8_1;
    matrix(indices[2]+2,indices[0]+2) += tmp_8_2;
    matrix(indices[2]+2,indices[1]+0) += tmp_8_3;
    matrix(indices[2]+2,indices[1]+1) += tmp_8_4;
    matrix(indices[2]+2,indices[1]+2) += tmp_8_5;
    matrix(indices[2]+2,indices[2]+0) += tmp_8_6;
    matrix(indices[2]+2,indices[2]+1) += tmp_8_7;
    matrix(indices[2]+2,indices[2]+2) += tmp_8_8;
    matrix(indices[2]+2,indices[3]+0) += tmp_8_9;
    matrix(indices[2]+2,indices[3]+1) += tmp_8_10;
    matrix(indices[2]+2,indices[3]+2) += tmp_8_11;
    matrix(indices[3]+0,indices[0]+0) += tmp_9_0;
    matrix(indices[3]+0,indices[0]+1) += tmp_9_1;
    matrix(indices[3]+0,indices[0]+2) += tmp_9_2;
    matrix(indices[3]+0,indices[1]+0) += tmp_9_3;
    matrix(indices[3]+0,indices[1]+1) += tmp_9_4;
    matrix(indices[3]+0,indices[1]+2) += tmp_9_5;
    matrix(indices[3]+0,indices[2]+0) += tmp_9_6;
    matrix(indices[3]+0,indices[2]+1) += tmp_9_7;
    matrix(indices[3]+0,indices[2]+2) += tmp_9_8;
    matrix(indices[3]+0,indices[3]+0) += tmp_9_9;
    matrix(indices[3]+0,indices[3]+1) += tmp_9_10;
    matrix(indices[3]+0,indices[3]+2) += tmp_9_11;
    matrix(indices[3]+1,indices[0]+0) += tmp_10_0;
    matrix(indices[3]+1,indices[0]+1) += tmp_10_1;
    matrix(indices[3]+1,indices[0]+2) += tmp_10_2;
    matrix(indices[3]+1,indices[1]+0) += tmp_10_3;
    matrix(indices[3]+1,indices[1]+1) += tmp_10_4;
    matrix(indices[3]+1,indices[1]+2) += tmp_10_5;
    matrix(indices[3]+1,indices[2]+0) += tmp_10_6;
    matrix(indices[3]+1,indices[2]+1) += tmp_10_7;
    matrix(indices[3]+1,indices[2]+2) += tmp_10_8;
    matrix(indices[3]+1,indices[3]+0) += tmp_10_9;
    matrix(indices[3]+1,indices[3]+1) += tmp_10_10;
    matrix(indices[3]+1,indices[3]+2) += tmp_10_11;
    matrix(indices[3]+2,indices[0]+0) += tmp_11_0;
    matrix(indices[3]+2,indices[0]+1) += tmp_11_1;
    matrix(indices[3]+2,indices[0]+2) += tmp_11_2;
    matrix(indices[3]+2,indices[1]+0) += tmp_11_3;
    matrix(indices[3]+2,indices[1]+1) += tmp_11_4;
    matrix(indices[3]+2,indices[1]+2) += tmp_11_5;
    matrix(indices[3]+2,indices[2]+0) += tmp_11_6;
    matrix(indices[3]+2,indices[2]+1) += tmp_11_7;
    matrix(indices[3]+2,indices[2]+2) += tmp_11_8;
    matrix(indices[3]+2,indices[3]+0) += tmp_11_9;
    matrix(indices[3]+2,indices[3]+1) += tmp_11_10;
    matrix(indices[3]+2,indices[3]+2) += tmp_11_11;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=elem.pos(2)[2]-elem.pos(0)[2]; T reg3=elem.pos(2)[1]-elem.pos(0)[1];
    T reg4=elem.pos(3)[1]-elem.pos(0)[1]; T reg5=elem.pos(3)[2]-elem.pos(0)[2]; T reg6=elem.pos(1)[2]-elem.pos(0)[2]; T reg7=elem.pos(1)[1]-elem.pos(0)[1]; T reg8=reg6*reg4;
    reg0=reg0*reg1; T reg9=reg2*reg4; T reg10=reg3*reg5; T reg11=reg7*reg5; T reg12=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg13=1.0/(*f.m).elastic_modulus; T reg14=elem.pos(2)[0]-elem.pos(0)[0]; T reg15=reg13*reg1; reg1=reg12*reg1; T reg16=elem.pos(1)[0]-elem.pos(0)[0];
    reg9=reg10-reg9; reg10=reg13*reg0; reg0=reg12*reg0; T reg17=reg6*reg3; reg8=reg11-reg8;
    reg11=reg7*reg2; reg17=reg11-reg17; reg11=elem.pos(3)[0]-elem.pos(0)[0]; T reg18=reg13*reg10; T reg19=reg12*reg15;
    T reg20=reg14*reg8; T reg21=reg12*reg0; T reg22=reg12*reg1; reg15=reg13*reg15; reg10=reg12*reg10;
    T reg23=reg16*reg9; T reg24=reg11*reg17; reg20=reg23-reg20; reg23=reg14*reg5; T reg25=reg2*reg11;
    reg5=reg16*reg5; T reg26=reg6*reg11; reg19=reg22+reg19; reg15=reg15-reg22; reg1=reg13*reg1;
    reg0=reg13*reg0; reg10=reg21+reg10; reg18=reg18-reg21; T reg27=reg7*reg11; reg26=reg5-reg26;
    reg5=reg16*reg4; reg2=reg16*reg2; reg11=reg3*reg11; reg25=reg23-reg25; reg4=reg14*reg4;
    reg6=reg6*reg14; reg23=reg12*reg10; T reg28=reg13*reg18; T reg29=reg22+reg1; reg24=reg20+reg24;
    reg19=reg12*reg19; reg0=reg21+reg0; reg15=reg13*reg15; reg26=reg26/reg24; reg27=reg5-reg27;
    reg29=reg12*reg29; reg8=reg8/reg24; reg11=reg4-reg11; reg19=reg15-reg19; reg25=reg25/reg24;
    reg9=reg9/reg24; reg12=reg12*reg0; reg14=reg7*reg14; reg23=reg28-reg23; reg3=reg16*reg3;
    reg6=reg2-reg6; reg6=reg6/reg24; reg2=reg25-reg26; reg14=reg3-reg14; reg17=reg17/reg24;
    reg3=reg8-reg9; reg27=reg27/reg24; reg11=reg11/reg24; reg12=reg23-reg12; reg29=reg19-reg29;
    reg3=reg3-reg17; reg29=reg29/reg12; reg2=reg6+reg2; reg4=reg27-reg11; reg5=0.5*reg6;
    reg7=0.5*reg17; reg13=0.5*reg8; reg15=0.5*reg26; reg14=reg14/reg24; reg16=0.5*reg9;
    reg19=0.5*reg27; reg20=reg29*reg15; reg21=0.5*reg2; reg23=reg29*reg13; reg28=reg29*reg5;
    T reg30=reg29*reg7; T reg31=0.5*reg25; T reg32=0.5*reg3; reg4=reg4-reg14; T reg33=0.5*reg14;
    reg18=reg18/reg12; T reg34=reg18*reg6; reg30=2*reg30; T reg35=2*reg28; T reg36=reg18*reg17;
    T reg37=reg18*reg26; T reg38=reg29*reg33; T reg39=reg29*reg32; T reg40=reg29*reg31; T reg41=reg18*reg27;
    T reg42=2*reg23; T reg43=reg18*reg14; reg20=2*reg20; reg10=reg10/reg12; T reg44=reg29*reg21;
    T reg45=reg18*reg8; reg12=reg0/reg12; reg0=reg29*reg16; T reg46=reg29*reg19; T reg47=0.5*reg4;
    T reg48=0.5*reg11; T reg49=reg10*reg26; T reg50=reg12*reg27; T reg51=2*reg46; reg38=2*reg38;
    T reg52=reg12*reg11; T reg53=reg13*reg30; T reg54=reg26*reg34; T reg55=reg8*reg36; T reg56=reg15*reg35;
    T reg57=reg11*reg41; T reg58=reg25*reg37; T reg59=reg16*reg42; T reg60=reg27*reg43; T reg61=reg31*reg20;
    T reg62=reg9*reg45; T reg63=reg12*reg6; T reg64=reg18*reg11; T reg65=reg18*reg4; T reg66=reg10*reg17;
    T reg67=reg10*reg8; T reg68=reg18*reg25; T reg69=reg18*reg2; T reg70=reg12*reg14; T reg71=reg10*reg6;
    T reg72=reg29*reg48; T reg73=reg18*reg9; T reg74=reg29*reg47; reg39=2*reg39; T reg75=reg18*reg3;
    reg0=2*reg0; T reg76=2*reg40; T reg77=reg10*reg25; reg44=2*reg44; T reg78=reg32*reg0;
    T reg79=reg16*reg35; T reg80=reg25*reg66; T reg81=reg8*reg49; T reg82=reg6*reg34; T reg83=reg15*reg42;
    T reg84=reg4*reg43; T reg85=reg7*reg30; T reg86=reg27*reg66; T reg87=reg16*reg30; T reg88=reg2*reg68;
    T reg89=reg9*reg73; T reg90=reg31*reg76; T reg91=reg10*reg9; T reg92=reg9*reg77; T reg93=reg31*reg0;
    reg74=2*reg74; reg55=reg56+reg55; T reg94=reg19*reg38; T reg95=reg32*reg39; T reg96=reg2*reg69;
    T reg97=reg48*reg42; reg60=reg53+reg60; T reg98=reg12*reg25; T reg99=reg4*reg64; T reg100=reg16*reg51;
    T reg101=reg11*reg67; T reg102=reg4*reg67; T reg103=reg32*reg51; T reg104=reg4*reg65; T reg105=reg59+reg57;
    T reg106=reg11*reg64; T reg107=reg27*reg41; T reg108=reg12*reg26; T reg109=reg32*reg30; T reg110=reg2*reg34;
    T reg111=reg11*reg63; T reg112=reg31*reg38; T reg113=reg10*reg2; T reg114=reg4*reg41; T reg115=reg25*reg70;
    T reg116=reg48*reg35; T reg117=reg11*reg43; T reg118=reg5*reg35; T reg119=reg13*reg38; T reg120=reg17*reg36;
    T reg121=reg32*reg42; T reg122=reg2*reg37; T reg123=reg15*reg20; T reg124=reg8*reg45; T reg125=reg25*reg34;
    T reg126=reg3*reg36; reg72=2*reg72; T reg127=reg9*reg50; T reg128=reg5*reg30; T reg129=reg26*reg37;
    T reg130=reg17*reg71; T reg131=reg3*reg50; T reg132=reg3*reg73; T reg133=reg47*reg42; T reg134=reg13*reg42;
    reg36=reg9*reg36; T reg135=reg31*reg35; T reg136=reg21*reg76; reg53=reg54+reg53; T reg137=reg26*reg50;
    T reg138=reg19*reg20; T reg139=reg21*reg20; T reg140=reg25*reg52; T reg141=reg9*reg71; T reg142=reg3*reg45;
    T reg143=reg31*reg30; T reg144=reg48*reg76; T reg145=reg25*reg68; reg43=reg14*reg43; T reg146=reg16*reg0;
    T reg147=reg62+reg61; T reg148=reg48*reg51; T reg149=reg33*reg35; T reg150=reg75*reg3; T reg151=reg21*reg44;
    T reg152=reg8*reg70; T reg153=reg6*reg70; T reg154=reg12*reg4; reg58=reg59+reg58; T reg155=reg21*reg35;
    T reg156=reg19*reg30; T reg157=reg116+reg115; T reg158=reg48*reg0; reg60=reg56+reg60; T reg159=reg4*reg98;
    T reg160=reg4*reg108; reg150=reg150+reg151; T reg161=reg21*reg51; T reg162=reg9*reg52; T reg163=reg33*reg38;
    reg143=reg141+reg143; T reg164=reg101+reg100; T reg165=reg3*reg113; reg106=reg146+reg106; T reg166=reg102+reg103;
    T reg167=reg15*reg38; reg99=reg78+reg99; T reg168=reg148+reg58; T reg169=reg27*reg63; T reg170=reg9*reg70;
    T reg171=reg48*reg30; reg146=reg146+reg145; T reg172=reg25*reg67; T reg173=reg31*reg42; T reg174=reg16*reg20;
    reg80=reg79+reg80; reg84=reg109+reg84; T reg175=reg9*reg49; T reg176=reg97+reg127; reg89=reg89+reg90;
    T reg177=reg48*reg72; reg86=reg119+reg86; reg119=reg4*reg63; T reg178=reg87+reg125; T reg179=reg25*reg50;
    T reg180=reg17*reg70; T reg181=reg21*reg38; T reg182=reg32*reg38; T reg183=reg4*reg66; reg93=reg92+reg93;
    reg36=reg36+reg135; T reg184=reg48*reg38; T reg185=reg48*reg20; T reg186=reg121+reg114; T reg187=reg147+reg148;
    T reg188=reg144+reg140; reg120=reg120+reg118; T reg189=reg3*reg71; T reg190=reg3*reg154; T reg191=reg21*reg30;
    T reg192=reg8*reg71; T reg193=reg47*reg30; T reg194=reg3*reg70; T reg195=reg15*reg30; T reg196=reg55+reg94;
    T reg197=reg47*reg39; reg96=reg96+reg95; T reg198=reg47*reg44; T reg199=reg2*reg154; T reg200=reg19*reg42;
    T reg201=reg2*reg91; T reg202=reg32*reg76; T reg203=reg8*reg50; T reg204=reg26*reg70; T reg205=reg85+reg82;
    reg78=reg78-reg88; reg81=reg83+reg81; T reg206=reg47*reg76; T reg207=reg47*reg0; T reg208=reg21*reg0;
    T reg209=reg3*reg52; T reg210=reg13*reg35; T reg211=reg26*reg66; T reg212=reg3*reg77; T reg213=reg47*reg51;
    reg43=reg85+reg43; reg138=reg137+reg138; reg85=reg139-reg142; reg132=reg132-reg136; T reg214=reg3*reg49;
    T reg215=reg21*reg42; reg129=reg129+reg134; T reg216=reg133+reg131; reg94=reg94+reg53; T reg217=reg149+reg153;
    T reg218=reg47*reg72; T reg219=reg47*reg38; reg126=reg126-reg155; reg156=reg152+reg156; T reg220=reg2*reg50;
    T reg221=reg32*reg72; T reg222=reg4*reg91; reg112=reg111+reg112; T reg223=reg21*reg39; T reg224=reg2*reg66;
    reg128=reg130+reg128; T reg225=reg32*reg35; reg30=reg33*reg30; T reg226=reg134+reg107; reg109=reg109-reg110;
    T reg227=reg74*reg47; reg104=reg95+reg104; reg95=reg16*reg38; T reg228=reg11*reg66; T reg229=reg47*reg35;
    reg70=reg2*reg70; reg61=reg61+reg105; T reg230=reg2*reg52; T reg231=reg19*reg35; T reg232=reg21*reg72;
    T reg233=reg2*reg67; T reg234=reg32*reg20; T reg235=reg19*reg51; T reg236=reg123+reg124; T reg237=reg11*reg108;
    reg122=reg122-reg121; T reg238=reg31*reg51; reg117=reg87+reg117; reg87=reg47*reg20; reg123=reg123+reg226;
    T reg239=reg24*reg164; T reg240=reg24*reg188; reg106=reg90+reg106; reg146=reg177+reg146; T reg241=reg24*reg138;
    reg237=reg237+reg238; reg211=reg211+reg210; reg129=reg235+reg129; T reg242=reg24*reg80; reg236=reg236+reg235;
    T reg243=reg24*reg81; reg185=reg185+reg179; reg204=reg204+reg231; T reg244=reg203+reg200; reg178=reg184+reg178;
    reg117=reg135+reg117; T reg245=reg24*reg86; T reg246=reg24*reg168; T reg247=reg24*reg196; T reg248=reg24*reg112;
    reg195=reg195+reg192; T reg249=reg24*reg94; reg95=reg228+reg95; reg228=reg24*reg157; reg174=reg174+reg172;
    T reg250=reg24*reg156; T reg251=reg24*reg61; reg126=reg219+reg126; reg84=reg84-reg155; reg190=reg197+reg190;
    reg181=reg181-reg119; reg191=reg191-reg189; reg194=reg193+reg194; reg182=reg183+reg182; reg139=reg139-reg186;
    reg160=reg160-reg161; reg96=reg227+reg96; reg183=reg24*reg166; reg120=reg120+reg163; reg199=reg198+reg199;
    reg99=reg99-reg136; reg205=reg163+reg205; reg201=reg201-reg202; reg232=reg232-reg159; reg221=reg222+reg221;
    reg78=reg218+reg78; reg104=reg151+reg104; reg151=reg24*reg128; reg70=reg70-reg229; reg230=reg230-reg206;
    reg109=reg219+reg109; reg234=reg234-reg233; reg223=reg165+reg223; reg224=reg224-reg225; reg30=reg180+reg30;
    reg122=reg122-reg213; reg87=reg87-reg220; reg208=reg208-reg212; reg171=reg170+reg171; reg209=reg207+reg209;
    reg163=reg24*reg143; reg167=reg167+reg169; reg184=reg36+reg184; reg132=reg218+reg132; reg36=reg24*reg176;
    reg85=reg85-reg213; reg175=reg175+reg173; reg214=reg214-reg215; reg165=reg24*reg187; reg150=reg227+reg150;
    reg43=reg118+reg43; reg177=reg89+reg177; reg89=reg24*reg60; reg170=reg24*reg93; reg180=reg24*reg216;
    reg193=reg24*reg217; reg158=reg162+reg158; reg122=reg24*reg122; reg117=reg24*reg117; reg162=ponderation*reg247;
    reg211=reg24*reg211; reg208=reg24*reg208; reg197=ponderation*reg250; reg234=reg24*reg234; reg126=reg24*reg126;
    reg223=reg24*reg223; reg190=reg24*reg190; reg236=reg24*reg236; reg230=reg24*reg230; reg209=reg24*reg209;
    reg198=ponderation*reg180; reg78=reg24*reg78; reg96=reg24*reg96; reg207=ponderation*reg249; reg205=reg24*reg205;
    reg129=reg24*reg129; reg244=reg24*reg244; reg218=ponderation*reg193; reg199=reg24*reg199; reg191=reg24*reg191;
    reg85=reg24*reg85; reg194=reg24*reg194; reg201=reg24*reg201; reg219=ponderation*reg243; reg214=reg24*reg214;
    reg43=reg24*reg43; reg195=reg24*reg195; reg132=reg24*reg132; reg222=ponderation*reg241; reg227=ponderation*reg228;
    reg160=reg24*reg160; reg139=reg24*reg139; reg178=reg24*reg178; reg182=reg24*reg182; reg181=reg24*reg181;
    T reg252=ponderation*reg242; reg150=reg24*reg150; reg84=reg24*reg84; T reg253=ponderation*reg89; reg185=reg24*reg185;
    reg177=reg24*reg177; T reg254=ponderation*reg170; T reg255=ponderation*reg246; reg158=reg24*reg158; T reg256=ponderation*reg165;
    T reg257=ponderation*reg245; reg174=reg24*reg174; reg175=reg24*reg175; T reg258=ponderation*reg36; reg167=reg24*reg167;
    T reg259=ponderation*reg240; reg184=reg24*reg184; reg146=reg24*reg146; T reg260=ponderation*reg163; reg171=reg24*reg171;
    reg30=reg24*reg30; reg204=reg24*reg204; T reg261=ponderation*reg248; reg87=reg24*reg87; reg224=reg24*reg224;
    reg95=reg24*reg95; reg109=reg24*reg109; T reg262=ponderation*reg151; T reg263=ponderation*reg251; reg70=reg24*reg70;
    reg104=reg24*reg104; reg237=reg24*reg237; reg221=reg24*reg221; T reg264=ponderation*reg239; reg232=reg24*reg232;
    reg120=reg24*reg120; reg99=reg24*reg99; T reg265=ponderation*reg183; reg106=reg24*reg106; reg123=reg24*reg123;
    T tmp_10_10=ponderation*reg205; T tmp_8_8=ponderation*reg123; T tmp_8_9=-reg257; T tmp_10_11=-reg218; T tmp_9_11=ponderation*reg30;
    T tmp_7_11=ponderation*reg204; T tmp_8_10=ponderation*reg167; T tmp_7_10=-reg207; T tmp_9_10=-reg262; T tmp_11_11=ponderation*reg43;
    T tmp_8_11=-reg253; T tmp_9_9=ponderation*reg120; T tmp_2_10=ponderation*reg181; T tmp_2_9=ponderation*reg182; T tmp_2_8=ponderation*reg139;
    T tmp_2_7=ponderation*reg160; T tmp_2_6=-reg265; T tmp_2_5=ponderation*reg99; T tmp_2_4=ponderation*reg232; T tmp_2_3=ponderation*reg221;
    T tmp_2_2=ponderation*reg104; T tmp_1_11=ponderation*reg70; T tmp_1_10=ponderation*reg109; T tmp_4_5=-reg259; T tmp_1_9=ponderation*reg224;
    T tmp_1_8=ponderation*reg87; T tmp_1_7=ponderation*reg122; T tmp_1_6=ponderation*reg234; T tmp_0_0=ponderation*reg150; T tmp_0_1=ponderation*reg223;
    T tmp_0_2=ponderation*reg190; T tmp_0_3=ponderation*reg132; T tmp_0_4=ponderation*reg208; T tmp_0_5=ponderation*reg209; T tmp_0_6=ponderation*reg85;
    T tmp_0_7=ponderation*reg214; T tmp_0_8=-reg198; T tmp_0_9=ponderation*reg126; T tmp_0_10=ponderation*reg191; T tmp_0_11=ponderation*reg194;
    T tmp_1_1=ponderation*reg96; T tmp_1_2=ponderation*reg199; T tmp_1_3=ponderation*reg201; T tmp_1_4=ponderation*reg78; T tmp_1_5=ponderation*reg230;
    T tmp_7_9=ponderation*reg211; T tmp_7_8=-reg222; T tmp_7_7=ponderation*reg129; T tmp_6_11=-reg197; T tmp_6_10=ponderation*reg195;
    T tmp_6_9=-reg162; T tmp_6_8=ponderation*reg244; T tmp_6_7=-reg219; T tmp_6_6=ponderation*reg236; T tmp_5_11=ponderation*reg117;
    T tmp_5_10=-reg261; T tmp_5_9=ponderation*reg95; T tmp_5_8=-reg263; T tmp_5_7=ponderation*reg237; T tmp_5_6=-reg264;
    T tmp_5_5=ponderation*reg106; T tmp_2_11=ponderation*reg84; T tmp_3_3=ponderation*reg177; T tmp_3_4=-reg254; T tmp_3_5=ponderation*reg158;
    T tmp_3_6=-reg256; T tmp_3_7=ponderation*reg175; T tmp_3_8=-reg258; T tmp_3_9=ponderation*reg184; T tmp_3_10=-reg260;
    T tmp_3_11=ponderation*reg171; T tmp_4_4=ponderation*reg146; T tmp_4_6=ponderation*reg174; T tmp_4_7=-reg255; T tmp_4_8=ponderation*reg185;
    T tmp_4_9=-reg252; T tmp_4_10=ponderation*reg178; T tmp_4_11=-reg227;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[0]+2) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_3;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_4;
    matrix(indices[0]+0,indices[1]+2) += tmp_0_5;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_7;
    matrix(indices[0]+0,indices[2]+2) += tmp_0_8;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_9;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_10;
    matrix(indices[0]+0,indices[3]+2) += tmp_0_11;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[0]+2) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_3;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_4;
    matrix(indices[0]+1,indices[1]+2) += tmp_1_5;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_7;
    matrix(indices[0]+1,indices[2]+2) += tmp_1_8;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_9;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_10;
    matrix(indices[0]+1,indices[3]+2) += tmp_1_11;
    matrix(indices[0]+2,indices[0]+2) += tmp_2_2;
    matrix(indices[0]+2,indices[1]+0) += tmp_2_3;
    matrix(indices[0]+2,indices[1]+1) += tmp_2_4;
    matrix(indices[0]+2,indices[1]+2) += tmp_2_5;
    matrix(indices[0]+2,indices[2]+0) += tmp_2_6;
    matrix(indices[0]+2,indices[2]+1) += tmp_2_7;
    matrix(indices[0]+2,indices[2]+2) += tmp_2_8;
    matrix(indices[0]+2,indices[3]+0) += tmp_2_9;
    matrix(indices[0]+2,indices[3]+1) += tmp_2_10;
    matrix(indices[0]+2,indices[3]+2) += tmp_2_11;
    matrix(indices[1]+0,indices[1]+0) += tmp_3_3;
    matrix(indices[1]+0,indices[1]+1) += tmp_3_4;
    matrix(indices[1]+0,indices[1]+2) += tmp_3_5;
    matrix(indices[1]+0,indices[2]+0) += tmp_3_6;
    matrix(indices[1]+0,indices[2]+1) += tmp_3_7;
    matrix(indices[1]+0,indices[2]+2) += tmp_3_8;
    matrix(indices[1]+0,indices[3]+0) += tmp_3_9;
    matrix(indices[1]+0,indices[3]+1) += tmp_3_10;
    matrix(indices[1]+0,indices[3]+2) += tmp_3_11;
    matrix(indices[1]+1,indices[1]+1) += tmp_4_4;
    matrix(indices[1]+1,indices[1]+2) += tmp_4_5;
    matrix(indices[1]+1,indices[2]+0) += tmp_4_6;
    matrix(indices[1]+1,indices[2]+1) += tmp_4_7;
    matrix(indices[1]+1,indices[2]+2) += tmp_4_8;
    matrix(indices[1]+1,indices[3]+0) += tmp_4_9;
    matrix(indices[1]+1,indices[3]+1) += tmp_4_10;
    matrix(indices[1]+1,indices[3]+2) += tmp_4_11;
    matrix(indices[1]+2,indices[1]+2) += tmp_5_5;
    matrix(indices[1]+2,indices[2]+0) += tmp_5_6;
    matrix(indices[1]+2,indices[2]+1) += tmp_5_7;
    matrix(indices[1]+2,indices[2]+2) += tmp_5_8;
    matrix(indices[1]+2,indices[3]+0) += tmp_5_9;
    matrix(indices[1]+2,indices[3]+1) += tmp_5_10;
    matrix(indices[1]+2,indices[3]+2) += tmp_5_11;
    matrix(indices[2]+0,indices[2]+0) += tmp_6_6;
    matrix(indices[2]+0,indices[2]+1) += tmp_6_7;
    matrix(indices[2]+0,indices[2]+2) += tmp_6_8;
    matrix(indices[2]+0,indices[3]+0) += tmp_6_9;
    matrix(indices[2]+0,indices[3]+1) += tmp_6_10;
    matrix(indices[2]+0,indices[3]+2) += tmp_6_11;
    matrix(indices[2]+1,indices[2]+1) += tmp_7_7;
    matrix(indices[2]+1,indices[2]+2) += tmp_7_8;
    matrix(indices[2]+1,indices[3]+0) += tmp_7_9;
    matrix(indices[2]+1,indices[3]+1) += tmp_7_10;
    matrix(indices[2]+1,indices[3]+2) += tmp_7_11;
    matrix(indices[2]+2,indices[2]+2) += tmp_8_8;
    matrix(indices[2]+2,indices[3]+0) += tmp_8_9;
    matrix(indices[2]+2,indices[3]+1) += tmp_8_10;
    matrix(indices[2]+2,indices[3]+2) += tmp_8_11;
    matrix(indices[3]+0,indices[3]+0) += tmp_9_9;
    matrix(indices[3]+0,indices[3]+1) += tmp_9_10;
    matrix(indices[3]+0,indices[3]+2) += tmp_9_11;
    matrix(indices[3]+1,indices[3]+1) += tmp_10_10;
    matrix(indices[3]+1,indices[3]+2) += tmp_10_11;
    matrix(indices[3]+2,indices[3]+2) += tmp_11_11;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim,unsigned symmetric_version>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); reg0=reg0*reg1; T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg3=1.0/(*f.m).elastic_modulus; T reg4=reg3*reg0; reg0=reg2*reg0; T reg5=reg3*reg4; T reg6=reg2*reg0;
    reg4=reg2*reg4; T reg7=elem.pos(3)[2]-elem.pos(0)[2]; T reg8=elem.pos(3)[1]-elem.pos(0)[1]; T reg9=elem.pos(2)[2]-elem.pos(0)[2]; T reg10=elem.pos(2)[1]-elem.pos(0)[1];
    T reg11=elem.pos(1)[2]-elem.pos(0)[2]; T reg12=elem.pos(1)[1]-elem.pos(0)[1]; reg0=reg3*reg0; reg4=reg6+reg4; reg5=reg5-reg6;
    T reg13=reg2*reg4; T reg14=reg11*reg8; reg0=reg6+reg0; reg6=reg3*reg5; T reg15=reg9*reg8;
    T reg16=reg10*reg7; T reg17=reg12*reg7; T reg18=elem.pos(2)[0]-elem.pos(0)[0]; T reg19=elem.pos(1)[0]-elem.pos(0)[0]; reg14=reg17-reg14;
    reg17=reg12*reg9; reg13=reg6-reg13; reg6=reg11*reg10; reg15=reg16-reg15; reg16=reg2*reg0;
    T reg20=reg19*reg15; T reg21=elem.pos(3)[0]-elem.pos(0)[0]; T reg22=reg18*reg14; reg16=reg13-reg16; reg6=reg17-reg6;
    reg13=reg9*reg21; reg17=reg11*reg21; T reg23=reg18*reg8; T reg24=reg18*reg7; T reg25=reg10*reg21;
    reg7=reg19*reg7; reg8=reg19*reg8; reg22=reg20-reg22; reg20=reg21*reg6; T reg26=(*f.m).alpha*(*f.m).deltaT;
    reg0=reg0/reg16; reg4=reg4/reg16; reg5=reg5/reg16; reg21=reg12*reg21; T reg27=reg4*reg26;
    T reg28=reg5*reg26; T reg29=reg0*reg26; reg12=reg12*reg18; reg17=reg7-reg17; reg25=reg23-reg25;
    reg21=reg8-reg21; reg9=reg19*reg9; reg18=reg11*reg18; reg10=reg19*reg10; reg13=reg24-reg13;
    reg20=reg22+reg20; reg7=reg28+reg27; reg8=reg29+reg27; reg12=reg10-reg12; reg10=1-var_inter[0];
    reg18=reg9-reg18; reg25=reg25/reg20; reg13=reg13/reg20; reg14=reg14/reg20; reg21=reg21/reg20;
    reg15=reg15/reg20; reg17=reg17/reg20; reg9=reg14-reg15; reg11=reg13-reg17; reg19=reg29+reg7;
    reg22=reg21-reg25; reg23=reg28+reg8; reg12=reg12/reg20; reg10=reg10-var_inter[1]; reg18=reg18/reg20;
    reg6=reg6/reg20; reg24=reg13*reg19; T reg30=reg21*reg23; T reg31=reg18*reg19; T reg32=reg14*reg19;
    T reg33=var_inter[1]*(*f.m).f_vol[0]; T reg34=var_inter[0]*(*f.m).f_vol[1]; reg11=reg18+reg11; reg9=reg9-reg6; reg10=reg10-var_inter[2];
    T reg35=var_inter[1]*(*f.m).f_vol[2]; reg22=reg22-reg12; T reg36=var_inter[2]*(*f.m).f_vol[1]; T reg37=reg25*reg23; T reg38=reg32-reg33;
    T reg39=reg24-reg34; T reg40=reg17*reg19; T reg41=reg15*reg19; T reg42=var_inter[2]*(*f.m).f_vol[2]; T reg43=reg30-reg35;
    T reg44=reg22*reg23; T reg45=var_inter[1]*(*f.m).f_vol[1]; T reg46=var_inter[2]*(*f.m).f_vol[0]; T reg47=var_inter[0]*(*f.m).f_vol[2]; T reg48=var_inter[0]*(*f.m).f_vol[0];
    T reg49=reg10*(*f.m).f_vol[2]; T reg50=reg10*(*f.m).f_vol[1]; T reg51=reg10*(*f.m).f_vol[0]; T reg52=reg9*reg19; T reg53=reg12*reg23;
    T reg54=reg31-reg36; T reg55=reg11*reg19; T reg56=reg6*reg19; reg54=reg20*reg54; T reg57=reg45+reg40;
    reg43=reg20*reg43; reg38=reg20*reg38; T reg58=reg42+reg53; T reg59=reg46+reg56; T reg60=reg47+reg37;
    T reg61=reg50+reg55; reg39=reg20*reg39; T reg62=reg51+reg52; T reg63=reg49+reg44; T reg64=reg48+reg41;
    T reg65=reg20*reg59; T reg66=reg20*reg61; T reg67=reg20*reg62; reg43=ponderation*reg43; T reg68=reg20*reg63;
    reg54=ponderation*reg54; T reg69=reg20*reg57; T reg70=reg20*reg64; reg38=ponderation*reg38; reg39=ponderation*reg39;
    T reg71=reg20*reg58; T reg72=reg20*reg60; sollicitation[indices[3]+1]+=-reg54; reg54=ponderation*reg67; sollicitation[indices[0]+0]+=reg54;
    T reg73=ponderation*reg71; sollicitation[indices[3]+2]+=reg73; T reg74=ponderation*reg65; sollicitation[indices[3]+0]+=reg74; T reg75=ponderation*reg66;
    sollicitation[indices[0]+1]+=reg75; sollicitation[indices[2]+2]+=-reg43; reg43=ponderation*reg68; sollicitation[indices[0]+2]+=reg43; T reg76=ponderation*reg69;
    sollicitation[indices[2]+1]+=reg76; sollicitation[indices[2]+0]+=-reg38; reg38=ponderation*reg70; sollicitation[indices[1]+0]+=reg38; T reg77=ponderation*reg72;
    sollicitation[indices[1]+2]+=reg77; sollicitation[indices[1]+1]+=-reg39;
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim,unsigned symmetric_version>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_residual( TD ponderation, const TD *var_inter,
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices ) { 
  #define PNODE(N) (*elem.node(N))
    T reg0=elem.pos(3)[1]-elem.pos(0)[1]; T reg1=elem.pos(3)[2]-elem.pos(0)[2]; T reg2=1+(*f.m).poisson_ratio; T reg3=elem.pos(2)[2]-elem.pos(0)[2]; T reg4=elem.pos(2)[1]-elem.pos(0)[1];
    T reg5=elem.pos(1)[2]-elem.pos(0)[2]; T reg6=elem.pos(1)[1]-elem.pos(0)[1]; T reg7=reg5*reg0; T reg8=reg4*reg1; T reg9=reg6*reg1;
    T reg10=reg3*reg0; reg2=reg2/(*f.m).elastic_modulus; reg7=reg9-reg7; reg9=reg6*reg3; T reg11=pow(reg2,2);
    reg10=reg8-reg10; reg8=reg5*reg4; T reg12=elem.pos(2)[0]-elem.pos(0)[0]; T reg13=elem.pos(1)[0]-elem.pos(0)[0]; T reg14=1.0/(*f.m).elastic_modulus;
    T reg15=reg13*reg10; T reg16=(*f.m).poisson_ratio/(*f.m).elastic_modulus; reg8=reg9-reg8; reg9=elem.pos(3)[0]-elem.pos(0)[0]; reg2=reg2*reg11;
    T reg17=reg12*reg7; T reg18=reg6*reg9; T reg19=reg16*reg2; reg2=reg14*reg2; reg17=reg15-reg17;
    reg15=reg9*reg8; T reg20=reg12*reg1; T reg21=reg3*reg9; T reg22=reg12*reg0; reg1=reg13*reg1;
    T reg23=reg4*reg9; reg9=reg5*reg9; reg0=reg13*reg0; reg6=reg6*reg12; reg4=reg13*reg4;
    reg12=reg5*reg12; reg3=reg13*reg3; reg23=reg22-reg23; reg18=reg0-reg18; reg9=reg1-reg9;
    reg21=reg20-reg21; reg15=reg17+reg15; reg0=reg16*reg2; reg1=reg16*reg19; reg5=reg16*reg11;
    reg2=reg14*reg2; reg11=reg14*reg11; reg10=reg10/reg15; reg7=reg7/reg15; reg21=reg21/reg15;
    reg2=reg2-reg1; reg13=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg17=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg23=reg23/reg15; reg6=reg4-reg6;
    reg12=reg3-reg12; reg3=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; reg4=reg16*reg11; reg20=vectors[0][indices[2]+2]-vectors[0][indices[0]+2]; reg22=reg16*reg5;
    reg11=reg14*reg11; T reg24=vectors[0][indices[1]+2]-vectors[0][indices[0]+2]; reg18=reg18/reg15; reg19=reg14*reg19; reg0=reg1+reg0;
    T reg25=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; reg9=reg9/reg15; reg12=reg12/reg15; T reg26=reg10*reg24; T reg27=reg18*reg25;
    T reg28=reg20*reg7; T reg29=reg23*reg3; T reg30=reg7*reg17; T reg31=reg10*reg13; T reg32=reg21*reg13;
    T reg33=reg9*reg17; T reg34=vectors[0][indices[3]+0]-vectors[0][indices[0]+0]; T reg35=reg7*reg25; T reg36=reg10*reg3; T reg37=vectors[0][indices[3]+2]-vectors[0][indices[0]+2];
    reg3=reg21*reg3; reg25=reg9*reg25; reg6=reg6/reg15; T reg38=vectors[0][indices[3]+1]-vectors[0][indices[0]+1]; T reg39=reg16*reg0;
    reg11=reg11-reg22; reg5=reg14*reg5; reg19=reg1+reg19; reg4=reg22+reg4; reg8=reg8/reg15;
    reg1=reg14*reg2; T reg40=reg12*reg38; T reg41=reg12*reg34; reg27=reg29-reg27; reg30=reg31-reg30;
    reg29=reg8*reg34; reg11=reg14*reg11; reg3=reg25-reg3; reg39=reg1-reg39; reg1=reg21*reg24;
    reg14=reg8*reg38; reg35=reg36-reg35; reg32=reg33-reg32; reg25=reg16*reg19; reg34=reg6*reg34;
    reg13=reg23*reg13; reg31=reg22+reg5; reg28=reg26-reg28; reg26=reg8*reg37; reg24=reg23*reg24;
    reg33=reg20*reg9; reg4=reg16*reg4; reg17=reg18*reg17; reg20=reg20*reg18; reg17=reg13-reg17;
    reg14=reg30+reg14; reg25=reg39-reg25; reg13=reg12*reg37; reg31=reg16*reg31; reg16=(*f.m).alpha*(*f.m).deltaT;
    reg38=reg6*reg38; reg29=reg35+reg29; reg1=reg33-reg1; reg4=reg11-reg4; reg37=reg6*reg37;
    reg27=reg34+reg27; reg26=reg28+reg26; reg41=reg3-reg41; reg40=reg32-reg40; reg20=reg24-reg20;
    reg29=reg29-reg16; reg40=reg40-reg16; reg20=reg37+reg20; reg14=reg41+reg14; reg17=reg38+reg17;
    reg13=reg1-reg13; reg2=reg2/reg25; reg0=reg0/reg25; reg19=reg19/reg25; reg31=reg4-reg31;
    reg26=reg27+reg26; reg13=reg17+reg13; reg1=reg2*reg29; reg3=reg0*reg40; reg26=0.5*reg26;
    reg25=reg31/reg25; reg20=reg20-reg16; reg29=reg0*reg29; reg4=reg2*reg40; reg11=reg21-reg9;
    reg40=reg19*reg40; reg17=reg7-reg10; reg14=0.5*reg14; reg11=reg12+reg11; reg24=reg18-reg23;
    reg26=reg25*reg26; reg3=reg1+reg3; reg17=reg17-reg8; reg14=reg25*reg14; reg13=0.5*reg13;
    reg1=reg2*reg20; reg20=reg19*reg20; reg4=reg29+reg4; reg40=reg29+reg40; reg24=reg24-reg6;
    reg27=0.5*reg7; reg28=1-var_inter[0]; reg1=reg40+reg1; reg29=0.5*reg21; reg13=reg25*reg13;
    reg14=2*reg14; reg26=2*reg26; reg30=0.5*reg10; reg4=reg20+reg4; reg31=0.5*reg9;
    reg32=0.5*reg12; reg3=reg20+reg3; reg20=0.5*reg8; reg33=0.5*reg11; reg34=0.5*reg17;
    reg35=reg20*reg26; reg36=reg12*reg4; reg37=reg9*reg4; reg38=reg34*reg26; reg39=reg1*reg24;
    reg40=reg14*reg30; reg41=reg34*reg14; T reg42=reg14*reg31; T reg43=reg6*reg1; T reg44=0.5*reg18;
    reg13=2*reg13; reg28=reg28-var_inter[1]; T reg45=reg10*reg3; T reg46=reg14*reg29; T reg47=0.5*reg23;
    T reg48=reg21*reg4; T reg49=reg3*reg17; T reg50=reg14*reg33; T reg51=reg27*reg26; T reg52=reg27*reg14;
    T reg53=reg8*reg3; T reg54=0.5*reg6; T reg55=reg7*reg3; T reg56=reg26*reg30; T reg57=0.5*reg24;
    T reg58=reg32*reg14; T reg59=reg20*reg14; T reg60=reg23*reg1; T reg61=reg18*reg1; T reg62=reg4*reg11;
    T reg63=reg13*reg31; T reg64=reg13*reg29; reg45=reg45-reg46; T reg65=reg13*reg32; T reg66=reg47*reg26;
    reg59=reg59-reg36; T reg67=reg13*reg47; reg38=reg39+reg38; reg53=reg53-reg58; reg37=reg37-reg52;
    reg39=reg54*reg26; reg28=reg28-var_inter[2]; reg56=reg60+reg56; reg50=reg49+reg50; reg49=reg13*reg57;
    reg41=reg62+reg41; reg60=reg51+reg61; reg40=reg40-reg48; reg62=reg13*reg33; T reg68=reg57*reg26;
    reg35=reg43+reg35; reg42=reg42-reg55; reg43=reg13*reg44; T reg69=reg26*reg44; T reg70=reg13*reg54;
    T reg71=var_inter[2]*(*f.m).f_vol[0]; reg39=reg53+reg39; reg53=var_inter[0]*(*f.m).f_vol[2]; reg63=reg63-reg60; T reg72=var_inter[1]*(*f.m).f_vol[0];
    T reg73=var_inter[2]*(*f.m).f_vol[2]; reg56=reg56-reg64; T reg74=var_inter[1]*(*f.m).f_vol[2]; T reg75=var_inter[1]*(*f.m).f_vol[1]; reg40=reg67+reg40;
    reg67=var_inter[0]*(*f.m).f_vol[0]; reg66=reg45+reg66; reg59=reg70+reg59; reg45=reg28*(*f.m).f_vol[2]; reg38=reg62+reg38;
    reg62=var_inter[2]*(*f.m).f_vol[1]; reg37=reg37-reg43; reg42=reg42-reg69; reg41=reg49+reg41; reg35=reg35-reg65;
    reg50=reg68+reg50; reg49=reg28*(*f.m).f_vol[0]; reg68=reg28*(*f.m).f_vol[1]; reg70=var_inter[0]*(*f.m).f_vol[1]; reg40=reg40-reg70;
    reg56=reg56-reg53; reg42=reg42-reg72; reg66=reg66-reg67; reg35=reg35-reg73; reg38=reg38-reg45;
    reg37=reg37-reg75; reg41=reg41-reg68; reg59=reg59-reg62; reg39=reg39-reg71; reg50=reg50-reg49;
    reg63=reg63-reg74; reg35=reg15*reg35; reg59=reg15*reg59; reg56=reg15*reg56; reg37=reg15*reg37;
    reg39=reg15*reg39; reg40=reg15*reg40; reg50=reg15*reg50; reg63=reg15*reg63; reg41=reg15*reg41;
    reg38=reg15*reg38; reg42=reg15*reg42; reg66=reg15*reg66; sollicitation[indices[2]+2]+=ponderation*reg63; sollicitation[indices[3]+2]+=ponderation*reg35;
    sollicitation[indices[2]+0]+=ponderation*reg42; sollicitation[indices[1]+2]+=ponderation*reg56; sollicitation[indices[3]+1]+=ponderation*reg59; sollicitation[indices[2]+1]+=ponderation*reg37; sollicitation[indices[3]+0]+=ponderation*reg39;
    sollicitation[indices[1]+1]+=ponderation*reg40; sollicitation[indices[0]+0]+=ponderation*reg50; sollicitation[indices[0]+1]+=ponderation*reg41; sollicitation[indices[0]+2]+=ponderation*reg38; sollicitation[indices[1]+0]+=ponderation*reg66;
  #undef PNODE
}
// 
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_false_true_true
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_false_true_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_true_true_true
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_true_true_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_false_true_false
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_false_true_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_true_true_false
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_true_true_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_symmetric_version_false_true
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_symmetric_version_false_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE ,unsigned symmetric_version>
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_symmetric_version_false_false
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_symmetric_version_false_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE ,unsigned symmetric_version>
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_RESIDUAL_elasticity_isotropy_stat_Qstat
#define ADD_NODAL_RESIDUAL_elasticity_isotropy_stat_Qstat
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE>
void add_nodal_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Tetra,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}

#ifndef elasticity_isotropy_stat_Qstat_read_material_to_mesh
#define elasticity_isotropy_stat_Qstat_read_material_to_mesh
template<class TM, class T, bool wont_add_nz>
void read_material_to_mesh_(const XmlNode &n, Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f){ 
    if(n.has_attribute("elastic_modulus"))  
        n.get_attribute("elastic_modulus", f.m->elastic_modulus ); 
    else  
        std::cerr << "Warning using default value of elastic_modulus : " << f.m->elastic_modulus << std::endl; 

    if(n.has_attribute("density"))  
        n.get_attribute("density", f.m->density ); 
    else  
        std::cerr << "Warning using default value of density : " << f.m->density << std::endl; 

    if(n.has_attribute("deltaT"))  
        n.get_attribute("deltaT", f.m->deltaT ); 
    else  
        std::cerr << "Warning using default value of deltaT : " << f.m->deltaT << std::endl; 

    if(n.has_attribute("poisson_ratio"))  
        n.get_attribute("poisson_ratio", f.m->poisson_ratio ); 
    else  
        std::cerr << "Warning using default value of poisson_ratio : " << f.m->poisson_ratio << std::endl; 

    if(n.has_attribute("alpha"))  
        n.get_attribute("alpha", f.m->alpha ); 
    else  
        std::cerr << "Warning using default value of alpha : " << f.m->alpha << std::endl; 

    if(n.has_attribute("resolution"))  
        n.get_attribute("resolution", f.m->resolution ); 
    else  
        std::cerr << "Warning using default value of resolution : " << f.m->resolution << std::endl; 

    if(n.has_attribute("f_vol"))  
        n.get_attribute("f_vol", f.m->f_vol ); 
    else  
        std::cerr << "Warning using default value of f_vol : " << f.m->f_vol << std::endl; 

  };
#endif // elasticity_isotropy_stat_Qstat_read_material_to_mesh
} // namespace LMT


#include "formulation/formulation.h"
namespace LMT {
#ifndef ELASTICITY_ISOTROPY_STAT_QSTAT_3_NUM_0_CARAC_H
#define ELASTICITY_ISOTROPY_STAT_QSTAT_3_NUM_0_CARAC_H
#ifndef STRUCT_ELASTICITY_ISOTROPY_STAT_QSTAT
#define STRUCT_ELASTICITY_ISOTROPY_STAT_QSTAT
struct elasticity_isotropy_stat_Qstat {
  static const char *name() { return "elasticity_isotropy_stat_Qstat"; }
};
#endif // ELASTICITY_ISOTROPY_STAT_QSTAT

template<class P_T>
class CaracFormulation<elasticity_isotropy_stat_Qstat,3,P_T>  {
public:
  typedef P_T T;
  static const char *name() { return "elasticity_isotropy_stat_Qstat"; }
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
  
  static const unsigned nb_nodal_unknowns = 3;
  template<class TE,class TTs,class Tvec>
  inline static void set_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice) {
    node.dep[2]=vecs[0][indice+2]; node.dep[0]=vecs[0][indice+0]; node.dep[1]=vecs[0][indice+1];
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
    node.dep[0]=vecs[0][indice+0]; node.dep[2]=vecs[0][indice+2]; node.dep[1]=vecs[0][indice+1];
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_nodal_initial_conditions(const TE &node,const TTs &f,Tvec &vecs,unsigned indice) {
    vecs[0][indice+2]=node.dep[2]; vecs[1][indice+2]=node.dep[2]; vecs[2][indice+2]=node.dep[2]; vecs[3][indice+2]=node.dep[2]; vecs[4][indice+2]=node.dep[2];
    vecs[0][indice+0]=node.dep[0]; vecs[1][indice+0]=node.dep[0]; vecs[2][indice+0]=node.dep[0]; vecs[3][indice+0]=node.dep[0]; vecs[4][indice+0]=node.dep[0];
    vecs[0][indice+1]=node.dep[1]; vecs[1][indice+1]=node.dep[1]; vecs[2][indice+1]=node.dep[1]; vecs[3][indice+1]=node.dep[1]; vecs[4][indice+1]=node.dep[1];
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_nodal_error(const TE &node,const TTs &f,const Tvec &vecs,int indice) {
    T reg0=vecs[1][indice+0]-vecs[0][indice+0]; T reg1=vecs[1][indice+1]-vecs[0][indice+1]; reg1=abs(reg1); T reg2=vecs[1][indice+2]-vecs[0][indice+2]; reg0=abs(reg0);
    reg0=max(reg1,reg0); reg2=abs(reg2); return max(reg2,reg0);
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_nodal(const TE &node,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
    old_vec[indice+2]=vecs[1][indice+2]; old_vec[indice+1]=vecs[1][indice+1]; old_vec[indice+0]=vecs[1][indice+0];
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
#endif // ELASTICITY_ISOTROPY_STAT_QSTAT_3_NUM_0_CARAC_H
extern "C" void apply_on_elements_after_solve_0_elasticity_isotropy_stat_Qstat_Wedge_0( double * );
extern "C" void apply_on_elements_after_solve_1_elasticity_isotropy_stat_Qstat_Wedge_1( double * );
extern "C" void apply_on_elements_after_solve_2_elasticity_isotropy_stat_Qstat_Wedge_2( double * );
extern "C" void apply_on_elements_after_solve_3_elasticity_isotropy_stat_Qstat_Wedge_3( double * );
extern "C" void apply_on_elements_after_solve_4_elasticity_isotropy_stat_Qstat_Wedge_4( double * );
extern "C" void apply_on_elements_after_solve_5_elasticity_isotropy_stat_Qstat_Wedge_5( double * );
extern "C" void apply_on_elements_after_solve_6_elasticity_isotropy_stat_Qstat_Wedge_6( double * );
extern "C" void apply_on_elements_after_solve_7_elasticity_isotropy_stat_Qstat_Wedge_7( double * );
extern "C" void apply_on_elements_after_solve_8_elasticity_isotropy_stat_Qstat_Wedge_8( double * );
extern "C" void apply_on_elements_after_solve_9_elasticity_isotropy_stat_Qstat_Wedge_9( double * );
extern "C" void apply_on_elements_after_solve_10_elasticity_isotropy_stat_Qstat_Wedge_10( double * );
extern "C" void apply_on_elements_after_solve_11_elasticity_isotropy_stat_Qstat_Wedge_11( double * );
extern "C" void apply_on_elements_after_solve_12_elasticity_isotropy_stat_Qstat_Wedge_12( double * );
extern "C" void apply_on_elements_after_solve_13_elasticity_isotropy_stat_Qstat_Wedge_13( double * );
extern "C" void apply_on_elements_after_solve_14_elasticity_isotropy_stat_Qstat_Wedge_14( double * );
class Wedge;
template<unsigned A,class B,class C> class Node;
template<class A,class B,class C,class D,unsigned E> class Element;

// Carac for ...
template<class P_T_pos,class P_ND,class TED,unsigned nim,class TM,class T>
class CaracFormulationForElement<elasticity_isotropy_stat_Qstat,Element<Wedge,DefaultBehavior,Node<3,P_T_pos,P_ND>,TED,nim>,TM,T> {
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
    T reg0=0.78867513459481286553*elem.pos(2)[2]; T reg1=0.78867513459481286553*elem.pos(0)[2]; T reg2=0.78867513459481286553*elem.pos(1)[2]; T reg3=0.78867513459481286553*elem.pos(2)[1]; T reg4=0.78867513459481286553*elem.pos(1)[1];
    T reg5=0.78867513459481286553*elem.pos(0)[1]; T reg6=0.5*elem.pos(0)[1]; T reg7=0.5*elem.pos(1)[1]; T reg8=0.5*elem.pos(1)[2]; T reg9=0.5*elem.pos(2)[2];
    T reg10=0.5*elem.pos(2)[1]; T reg11=0.5*elem.pos(0)[2]; T reg12=0.21132486540518713447*elem.pos(3)[2]; T reg13=0.21132486540518713447*elem.pos(3)[1]; T reg14=0.5*elem.pos(3)[1];
    T reg15=reg7+reg6; reg3=reg3-reg5; T reg16=0.5*elem.pos(4)[1]; T reg17=reg8+reg9; reg5=reg4-reg5;
    reg4=reg8+reg11; reg0=reg0-reg1; reg1=reg2-reg1; reg2=0.5*elem.pos(4)[2]; T reg18=reg7+reg10;
    T reg19=0.5*elem.pos(3)[2]; reg3=reg3-reg13; T reg20=0.21132486540518713447*elem.pos(5)[1]; T reg21=0.78867513459481286553*elem.pos(0)[0]; T reg22=0.78867513459481286553*elem.pos(1)[0];
    T reg23=0.21132486540518713447*elem.pos(2)[1]; T reg24=0.21132486540518713447*elem.pos(0)[1]; T reg25=0.21132486540518713447*elem.pos(2)[2]; T reg26=0.21132486540518713447*elem.pos(0)[2]; T reg27=0.21132486540518713447*elem.pos(1)[1];
    T reg28=0.21132486540518713447*elem.pos(1)[2]; reg4=reg19-reg4; T reg29=0.21132486540518713447*elem.pos(5)[2]; reg0=reg0-reg12; reg15=reg14-reg15;
    T reg30=0.78867513459481286553*elem.pos(2)[0]; reg13=reg5-reg13; reg5=0.21132486540518713447*elem.pos(4)[1]; reg12=reg1-reg12; reg1=0.21132486540518713447*elem.pos(4)[2];
    T reg31=reg6+reg10; T reg32=reg11+reg9; T reg33=reg16-reg18; T reg34=0.5*elem.pos(5)[1]; T reg35=reg2-reg17;
    T reg36=0.5*elem.pos(5)[2]; T reg37=0.78867513459481286553*elem.pos(3)[2]; reg25=reg25-reg26; T reg38=0.78867513459481286553*elem.pos(3)[1]; reg23=reg23-reg24;
    reg32=reg19-reg32; reg31=reg14-reg31; T reg39=0.5*elem.pos(2)[0]; T reg40=1+(*f.m).poisson_ratio; reg22=reg22-reg21;
    T reg41=0.21132486540518713447*elem.pos(3)[0]; reg3=reg20+reg3; reg4=reg4+reg2; reg0=reg29+reg0; reg15=reg15+reg16;
    reg20=0.5*elem.pos(0)[0]; reg21=reg30-reg21; reg29=0.5*elem.pos(1)[0]; reg5=reg13+reg5; reg1=reg12+reg1;
    reg26=reg28-reg26; reg33=reg34+reg33; reg24=reg27-reg24; reg35=reg36+reg35; reg21=reg21-reg41;
    reg26=reg26-reg37; reg31=reg34+reg31; reg12=0.78867513459481286553*elem.pos(4)[2]; reg41=reg22-reg41; reg13=0.21132486540518713447*elem.pos(4)[0];
    reg22=0.21132486540518713447*elem.pos(0)[0]; reg27=reg0*reg33; reg28=0.5*elem.pos(3)[0]; reg30=0.21132486540518713447*elem.pos(1)[0]; T reg42=0.21132486540518713447*elem.pos(5)[0];
    T reg43=reg29+reg20; T reg44=reg3*reg4; T reg45=reg0*reg15; T reg46=0.5*elem.pos(4)[0]; T reg47=reg39+reg29;
    T reg48=reg3*reg35; T reg49=reg4*reg5; reg40=reg40/(*f.m).elastic_modulus; T reg50=0.78867513459481286553*elem.pos(5)[2]; reg24=reg24-reg38;
    T reg51=reg15*reg1; reg38=reg23-reg38; reg23=0.21132486540518713447*elem.pos(2)[0]; T reg52=0.78867513459481286553*elem.pos(4)[1]; reg37=reg25-reg37;
    reg32=reg36+reg32; reg25=reg1*reg33; T reg53=reg5*reg35; T reg54=0.78867513459481286553*elem.pos(5)[1]; T reg55=reg3*reg1;
    T reg56=reg0*reg5; reg51=reg49-reg51; reg21=reg42+reg21; reg42=0.5*elem.pos(5)[0]; reg43=reg28-reg43;
    reg49=reg46-reg47; reg25=reg53-reg25; reg45=reg44-reg45; reg23=reg23-reg22; reg44=pow(reg40,2);
    reg27=reg48-reg27; reg52=reg24+reg52; reg12=reg26+reg12; reg38=reg54+reg38; reg13=reg41+reg13;
    reg24=reg0*reg31; reg26=reg3*reg32; reg41=reg5*reg32; reg48=reg1*reg31; reg53=reg39+reg20;
    reg22=reg30-reg22; reg30=0.78867513459481286553*elem.pos(3)[0]; reg37=reg50+reg37; reg50=0.78867513459481286553*PNODE(1).dep[1]; reg54=reg21*reg51;
    T reg57=0.78867513459481286553*PNODE(0).dep[1]; T reg58=0.78867513459481286553*PNODE(2).dep[1]; reg43=reg43+reg46; reg55=reg56-reg55; reg53=reg28-reg53;
    reg48=reg41-reg48; reg24=reg26-reg24; reg26=0.78867513459481286553*PNODE(1).dep[0]; reg41=0.78867513459481286553*PNODE(0).dep[0]; reg56=0.78867513459481286553*PNODE(2).dep[0];
    T reg59=reg21*reg25; T reg60=reg13*reg27; T reg61=reg4*reg38; reg40=reg40*reg44; T reg62=reg15*reg37;
    T reg63=(*f.m).poisson_ratio/(*f.m).elastic_modulus; reg49=reg49+reg42; T reg64=reg13*reg45; T reg65=0.78867513459481286553*elem.pos(5)[0]; reg23=reg23-reg30;
    T reg66=0.78867513459481286553*elem.pos(4)[0]; reg30=reg22-reg30; reg22=reg4*reg52; T reg67=reg15*reg12; T reg68=1.0/(*f.m).elastic_modulus;
    reg66=reg30+reg66; reg30=0.21132486540518713447*PNODE(3).dep[1]; T reg69=0.5*PNODE(2).dep[1]; reg58=reg58-reg57; T reg70=reg49*reg0;
    T reg71=reg21*reg35; T reg72=reg63*reg40; T reg73=reg49*reg1; T reg74=reg13*reg35; reg26=reg26-reg41;
    T reg75=0.21132486540518713447*PNODE(3).dep[0]; T reg76=0.5*PNODE(2).dep[0]; reg41=reg56-reg41; reg56=reg1*reg43; T reg77=reg13*reg4;
    T reg78=reg38*reg12; reg62=reg61-reg62; reg59=reg60-reg59; reg60=reg37*reg52; reg23=reg65+reg23;
    reg61=0.5*PNODE(0).dep[0]; reg65=0.5*PNODE(1).dep[0]; reg67=reg22-reg67; reg8=reg8-reg11; reg11=reg9-reg11;
    reg9=0.78867513459481286553*PNODE(0).dep[2]; reg22=0.78867513459481286553*PNODE(1).dep[2]; T reg79=reg43*reg55; T reg80=0.5*PNODE(0).dep[1]; T reg81=reg49*reg55;
    reg10=reg10-reg6; T reg82=reg35*reg38; reg54=reg64-reg54; reg64=reg33*reg37; T reg83=reg35*reg52;
    T reg84=reg33*reg12; T reg85=0.5*PNODE(1).dep[1]; reg6=reg7-reg6; reg7=reg0*reg43; T reg86=reg13*reg24;
    T reg87=reg4*reg21; reg53=reg42+reg53; reg57=reg50-reg57; reg40=reg68*reg40; reg50=reg21*reg48;
    T reg88=0.78867513459481286553*PNODE(2).dep[2]; T reg89=reg21*reg1; T reg90=reg23*reg67; T reg91=0.5*PNODE(3).dep[1]; T reg92=0.21132486540518713447*PNODE(2).dep[1];
    T reg93=reg32*reg38; T reg94=reg31*reg37; T reg95=0.5*PNODE(4).dep[1]; T reg96=reg80+reg85; reg78=reg60-reg78;
    reg60=0.21132486540518713447*PNODE(1).dep[0]; T reg97=0.21132486540518713447*PNODE(0).dep[0]; T reg98=0.21132486540518713447*PNODE(2).dep[0]; reg84=reg83-reg84; reg64=reg82-reg64;
    reg10=reg10-reg14; reg8=reg8-reg19; reg82=reg0*reg53; reg83=reg21*reg32; reg1=reg1*reg53;
    T reg99=reg13*reg32; T reg100=reg55*reg53; reg19=reg11-reg19; reg50=reg86-reg50; reg11=reg31*reg12;
    reg14=reg6-reg14; reg6=reg32*reg52; reg85=reg85+reg69; reg70=reg71-reg70; reg73=reg74-reg73;
    reg71=0.21132486540518713447*PNODE(1).dep[1]; reg74=reg65+reg76; reg86=0.5*PNODE(2).dep[2]; T reg101=reg49*reg5; T reg102=reg13*reg33;
    T reg103=reg49*reg3; T reg104=reg21*reg33; T reg105=reg66*reg62; T reg106=0.21132486540518713447*PNODE(0).dep[1]; reg26=reg26-reg75;
    T reg107=0.21132486540518713447*PNODE(5).dep[1]; reg56=reg77-reg56; reg77=reg13*reg15; T reg108=0.5*PNODE(1).dep[2]; T reg109=reg5*reg43;
    reg79=reg54+reg79; reg54=reg68*reg40; reg22=reg22-reg9; reg65=reg61+reg65; T reg110=0.21132486540518713447*PNODE(3).dep[2];
    T reg111=0.5*PNODE(3).dep[0]; T reg112=0.5*PNODE(4).dep[0]; reg75=reg41-reg75; reg9=reg88-reg9; reg40=reg63*reg40;
    reg41=reg63*reg72; reg88=0.21132486540518713447*PNODE(5).dep[0]; reg0=reg13*reg0; T reg113=reg3*reg43; T reg114=0.5*PNODE(0).dep[2];
    reg7=reg87-reg7; reg87=reg15*reg21; T reg115=0.21132486540518713447*PNODE(4).dep[1]; reg57=reg57-reg30; T reg116=0.21132486540518713447*PNODE(4).dep[0];
    reg30=reg58-reg30; reg59=reg81+reg59; reg100=reg50+reg100; reg65=reg111-reg65; reg50=reg13*reg31;
    reg58=reg21*reg31; reg81=reg3*reg53; T reg117=reg5*reg53; reg109=reg77-reg109; reg76=reg61+reg76;
    reg56=reg56/reg79; reg30=reg107+reg30; reg61=reg43*reg37; reg51=reg51/reg79; reg88=reg75+reg88;
    reg75=0.5*PNODE(5).dep[1]; reg77=reg108+reg114; reg70=reg70/reg59; reg107=reg95-reg85; reg73=reg73/reg59;
    T reg118=reg23*reg84; T reg119=reg66*reg64; reg10=reg34+reg10; reg34=0.5*PNODE(5).dep[0]; T reg120=reg112-reg74;
    reg26=reg116+reg26; reg19=reg36+reg19; reg29=reg29-reg20; reg25=reg25/reg59; reg27=reg27/reg59;
    reg20=reg39-reg20; reg36=0.5*PNODE(3).dep[2]; reg108=reg108+reg86; reg45=reg45/reg79; reg39=0.5*PNODE(4).dep[2];
    reg5=reg21*reg5; reg3=reg13*reg3; reg101=reg102-reg101; reg13=0.21132486540518713447*PNODE(5).dep[2]; reg9=reg9-reg110;
    reg103=reg104-reg103; reg21=reg43*reg78; reg102=0.78867513459481286553*PNODE(3).dep[1]; reg92=reg92-reg106; reg94=reg93-reg94;
    reg54=reg54-reg41; reg93=reg43*reg12; reg104=reg4*reg66; reg11=reg6-reg11; reg90=reg105-reg90;
    reg40=reg41+reg40; reg98=reg98-reg97; reg72=reg68*reg72; reg6=0.78867513459481286553*PNODE(3).dep[0]; reg97=reg60-reg97;
    reg8=reg2+reg8; reg110=reg22-reg110; reg14=reg16+reg14; reg2=0.21132486540518713447*PNODE(4).dep[2]; reg16=0.21132486540518713447*PNODE(2).dep[2];
    reg4=reg4*reg23; reg106=reg71-reg106; reg22=0.21132486540518713447*PNODE(1).dep[2]; reg60=0.21132486540518713447*PNODE(0).dep[2]; reg113=reg87-reg113;
    reg69=reg80+reg69; reg96=reg91-reg96; reg82=reg83-reg82; reg115=reg57+reg115; reg89=reg0-reg89;
    reg1=reg99-reg1; reg7=reg7/reg79; reg29=reg29-reg28; reg0=reg68*reg44; reg57=reg15*reg66;
    reg28=reg20-reg28; reg20=reg43*reg38; reg71=reg23*reg12; reg44=reg63*reg44; reg72=reg41+reg72;
    reg15=reg15*reg23; reg41=0.78867513459481286553*PNODE(3).dep[2]; reg80=reg33*reg19; reg107=reg75+reg107; reg83=reg49*reg78;
    reg118=reg119-reg118; reg22=reg22-reg60; reg87=reg35*reg10; reg21=reg90+reg21; reg76=reg111-reg76;
    reg48=reg48/reg100; reg24=reg24/reg100; reg61=reg4-reg61; reg86=reg114+reg86; reg117=reg50-reg117;
    reg1=reg1/reg100; reg110=reg2+reg110; reg81=reg58-reg81; reg82=reg82/reg100; reg69=reg91-reg69;
    reg2=0.78867513459481286553*PNODE(4).dep[1]; reg106=reg106-reg102; reg4=reg63*reg40; reg102=reg92-reg102; reg113=reg113/reg79;
    reg50=0.78867513459481286553*PNODE(5).dep[1]; reg93=reg104-reg93; reg58=0.78867513459481286553*PNODE(5).dep[0]; reg109=reg109/reg79; reg98=reg98-reg6;
    reg90=reg68*reg54; reg6=reg97-reg6; reg91=0.78867513459481286553*PNODE(4).dep[0]; reg60=reg16-reg60; reg13=reg9+reg13;
    reg9=reg33*reg8; reg5=reg3-reg5; reg3=reg35*reg14; reg43=reg43*reg52; reg16=reg66*reg37;
    reg120=reg120+reg34; reg92=reg55/reg59; reg97=reg88*reg25; reg99=reg26*reg27; reg104=0.5*PNODE(5).dep[2];
    reg105=reg39-reg108; reg101=reg101/reg59; reg111=reg26*reg45; reg103=reg103/reg59; reg114=reg88*reg51;
    reg116=reg55/reg79; reg65=reg112+reg65; reg112=reg56*reg30; reg119=reg35*reg66; T reg121=reg49*reg12;
    T reg122=reg35*reg23; T reg123=reg49*reg37; T reg124=reg23*reg11; T reg125=reg66*reg94; reg95=reg96+reg95;
    reg96=reg115*reg7; T reg126=reg89/reg79; T reg127=reg89/reg59; T reg128=reg115*reg70; T reg129=reg30*reg73;
    reg77=reg36-reg77; reg6=reg91+reg6; reg117=reg117/reg100; reg81=reg81/reg100; reg91=reg23*reg52;
    T reg130=reg66*reg38; T reg131=reg10*reg8; T reg132=0.78867513459481286553*PNODE(5).dep[2]; T reg133=reg110*reg113; T reg134=reg116*reg65;
    reg60=reg60-reg41; T reg135=reg49*reg38; reg114=reg111-reg114; reg111=reg33*reg23; reg83=reg118+reg83;
    reg118=reg110*reg103; T reg136=reg109*reg13; T reg137=reg19*reg14; T reg138=reg13*reg101; T reg139=reg5/reg79;
    T reg140=reg5/reg59; reg9=reg3-reg9; reg3=reg33*reg66; T reg141=reg115*reg82; T reg142=reg127*reg107;
    reg89=reg89/reg100; reg123=reg122-reg123; reg122=reg30*reg1; reg121=reg119-reg121; reg69=reg75+reg69;
    reg2=reg106+reg2; reg102=reg50+reg102; reg37=reg53*reg37; reg76=reg34+reg76; reg34=reg32*reg23;
    reg93=reg93/reg21; reg55=reg55/reg100; reg50=reg88*reg48; reg75=reg26*reg24; reg67=reg67/reg21;
    reg58=reg98+reg58; reg61=reg61/reg21; reg12=reg53*reg12; reg86=reg36-reg86; reg32=reg32*reg66;
    reg62=reg62/reg21; reg29=reg46+reg29; reg4=reg90-reg4; reg36=reg63*reg72; reg80=reg87-reg80;
    reg97=reg99-reg97; reg46=reg49*reg52; reg77=reg39+reg77; reg124=reg125-reg124; reg39=reg53*reg78;
    reg20=reg15-reg20; reg15=reg92*reg120; reg87=reg63*reg0; reg0=reg68*reg0; reg90=reg63*reg44;
    reg96=reg112-reg96; reg41=reg22-reg41; reg128=reg129-reg128; reg43=reg57-reg43; reg105=reg105+reg104;
    reg22=reg126*reg95; reg71=reg16-reg71; reg16=0.78867513459481286553*PNODE(4).dep[2]; reg28=reg42+reg28; reg36=reg4-reg36;
    reg4=reg89*reg69; reg22=reg96-reg22; reg42=reg78/reg21; reg123=reg123/reg83; reg57=reg58*reg67;
    reg39=reg124+reg39; reg96=reg93*reg102; reg23=reg31*reg23; reg41=reg16+reg41; reg43=reg43/reg21;
    reg66=reg31*reg66; reg52=reg53*reg52; reg16=reg88*reg56; reg31=reg26*reg7; reg98=reg28*reg9;
    reg99=reg45*reg115; reg132=reg60+reg132; reg60=(*f.m).deltaT*(*f.m).alpha; reg87=reg90+reg87; reg131=reg137-reg131;
    reg91=reg130-reg91; reg106=0.5*vectors[0][indices[2]+0]; reg112=reg51*reg30; reg20=reg20/reg21; reg0=reg0-reg90;
    reg119=reg6*reg62; reg44=reg68*reg44; reg38=reg53*reg38; reg12=reg32-reg12; reg50=reg75-reg50;
    reg84=reg84/reg83; reg86=reg104+reg86; reg5=reg5/reg100; reg32=reg13*reg117; reg64=reg64/reg83;
    reg53=reg110*reg81; reg134=reg114+reg134; reg135=reg111-reg135; reg75=reg139*reg77; reg136=reg133-reg136;
    reg138=reg118-reg138; reg104=reg2*reg61; reg111=0.5*vectors[0][indices[2]+1]; reg114=0.5*vectors[0][indices[0]+1]; reg118=reg140*reg105;
    reg124=reg71/reg21; reg46=reg3-reg46; reg3=reg80*reg29; reg125=0.5*vectors[0][indices[1]+1]; reg15=reg97+reg15;
    reg97=reg30*reg25; reg129=reg115*reg27; reg130=reg26*reg70; reg133=reg88*reg73; reg137=0.5*vectors[0][indices[0]+0];
    T reg143=0.5*vectors[0][indices[1]+0]; reg37=reg34-reg37; reg34=reg55*reg76; reg121=reg121/reg83; reg141=reg122-reg141;
    reg142=reg128-reg142; reg122=reg90+reg44; reg128=0.5*vectors[0][indices[2]+2]; reg38=reg23-reg38; reg23=reg125-reg114;
    T reg144=reg2*reg123; reg118=reg138+reg118; reg138=reg95*reg124; T reg145=reg49*reg8; T reg146=reg35*reg29;
    T reg147=0.5*vectors[0][indices[3]+1]; reg114=reg111-reg114; T reg148=reg71/reg83; reg104=reg96-reg104; reg52=reg66-reg52;
    reg66=reg143-reg137; reg96=reg102*reg121; reg130=reg133-reg130; reg133=reg120*reg127; T reg149=0.5*vectors[0][indices[3]+0];
    reg142=reg142-reg60; reg22=reg22-reg60; reg97=reg129-reg97; reg129=reg92*reg107; reg15=reg15-reg60;
    reg40=reg40/reg36; reg4=reg141-reg4; reg54=reg54/reg36; reg141=reg49*reg19; reg35=reg35*reg28;
    T reg150=reg26*reg103; T reg151=reg88*reg101; reg0=reg68*reg0; reg27=reg110*reg27; reg25=reg13*reg25;
    reg87=reg63*reg87; reg68=reg41*reg20; reg135=reg135/reg83; reg57=reg119-reg57; reg119=reg88*reg1;
    T reg152=0.5*vectors[0][indices[0]+2]; T reg153=reg49*reg131; T reg154=0.5*vectors[0][indices[1]+2]; T reg155=reg91/reg21; reg37=reg37/reg39;
    T reg156=reg5*reg86; T reg157=reg113*reg26; T reg158=reg6*reg64; reg112=reg99-reg112; reg99=reg116*reg95;
    reg137=reg106-reg137; T reg159=reg109*reg88; reg11=reg11/reg39; reg32=reg53-reg32; reg53=reg58*reg84;
    reg12=reg12/reg39; reg31=reg16-reg31; reg16=reg65*reg126; T reg160=reg65*reg42; reg94=reg94/reg39;
    T reg161=reg115*reg24; reg75=reg136+reg75; reg136=reg78/reg83; reg98=reg3-reg98; reg34=reg50+reg34;
    reg51=reg13*reg51; reg45=reg110*reg45; reg134=reg134-reg60; reg46=reg46/reg83; reg3=reg43*reg132;
    reg50=reg30*reg48; T reg162=reg26*reg82; T reg163=reg41*reg135; reg66=reg66-reg149; reg78=reg78/reg39;
    reg160=reg57+reg160; reg57=reg67*reg102; T reg164=reg120*reg136; reg52=reg52/reg39; T reg165=reg107*reg148;
    reg38=reg38/reg39; reg3=reg68-reg3; reg138=reg104-reg138; reg144=reg96-reg144; reg68=reg132*reg46;
    reg53=reg158-reg53; reg96=reg91/reg83; reg104=reg77*reg155; reg158=reg2*reg37; reg71=reg71/reg39;
    T reg166=0.5*vectors[0][indices[4]+0]; T reg167=reg102*reg12; T reg168=reg58*reg93; T reg169=reg6*reg61; T reg170=reg6*reg94;
    T reg171=reg58*reg11; T reg172=reg62*reg2; reg133=reg130-reg133; reg130=reg139*reg65; reg129=reg97+reg129;
    reg97=reg54*reg22; reg141=reg35-reg141; reg35=reg140*reg120; reg151=reg150-reg151; reg25=reg27-reg25;
    reg92=reg105*reg92; reg27=reg33*reg29; reg150=0.5*vectors[0][indices[4]+1]; reg23=reg23-reg147; reg103=reg115*reg103;
    reg101=reg30*reg101; reg73=reg13*reg73; reg70=reg110*reg70; reg99=reg112+reg99; reg16=reg31-reg16;
    reg147=reg114-reg147; reg31=0.5*vectors[0][indices[5]+1]; reg156=reg32+reg156; reg34=reg34-reg60; reg4=reg4-reg60;
    reg145=reg146-reg145; reg32=reg154-reg152; reg125=reg111+reg125; reg111=reg40*reg134; reg118=reg118-reg60;
    reg112=reg40*reg15; reg75=reg75-reg60; reg7=reg110*reg7; reg114=reg54*reg15; reg146=reg40*reg142;
    reg56=reg13*reg56; reg8=reg28*reg8; T reg173=reg54*reg142; reg109=reg109*reg30; reg113=reg113*reg115;
    T reg174=reg40*reg22; T reg175=0.5*vectors[0][indices[3]+2]; reg33=reg33*reg28; T reg176=reg54*reg134; reg152=reg128-reg152;
    T reg177=reg49*reg10; reg116=reg77*reg116; reg19=reg19*reg29; reg51=reg45-reg51; reg159=reg157-reg159;
    reg49=reg49*reg14; reg45=reg76*reg89; reg26=reg26*reg81; reg88=reg88*reg117; reg162=reg119-reg162;
    reg50=reg161-reg50; reg149=reg137-reg149; reg119=0.5*vectors[0][indices[5]+0]; reg137=reg55*reg69; reg153=reg98+reg153;
    reg106=reg143+reg106; reg122=reg63*reg122; reg87=reg0-reg87; reg48=reg13*reg48; reg24=reg110*reg24;
    reg72=reg72/reg36; reg32=reg32-reg175; reg49=reg27-reg49; reg81=reg115*reg81; reg80=reg80/reg153;
    reg174=reg176+reg174; reg104=reg3+reg104; reg175=reg152-reg175; reg128=reg154+reg128; reg116=reg51+reg116;
    reg0=reg132*reg52; reg165=reg144-reg165; reg55=reg86*reg55; reg159=reg130+reg159; reg48=reg24-reg48;
    reg67=reg132*reg67; reg62=reg41*reg62; reg3=reg43*reg58; reg24=0.5*vectors[0][indices[4]+2]; reg126=reg77*reg126;
    reg27=reg105*reg96; reg51=reg150-reg125; reg68=reg163-reg68; reg122=reg87-reg122; reg87=reg72*reg142;
    reg7=reg56-reg7; reg164=reg53+reg164; reg8=reg19-reg8; reg19=0.5*vectors[0][indices[5]+2]; reg53=reg41*reg38;
    reg146=reg114+reg146; reg56=reg72*reg118; reg14=reg28*reg14; reg29=reg10*reg29; reg109=reg113-reg109;
    reg82=reg110*reg82; reg1=reg13*reg1; reg158=reg167-reg158; reg173=reg112+reg173; reg10=reg69*reg71;
    reg13=reg72*reg75; reg117=reg30*reg117; reg139=reg139*reg95; reg149=reg149+reg119; reg9=reg9/reg153;
    reg28=reg58*reg121; reg140=reg140*reg107; reg137=reg50+reg137; reg30=reg40*reg4; reg171=reg170-reg171;
    reg50=reg54*reg34; reg101=reg103-reg101; reg99=reg16+reg99; reg45=reg162-reg45; reg138=reg138-reg60;
    reg70=reg73-reg70; reg127=reg105*reg127; reg16=reg72*reg22; reg160=reg160-reg60; reg147=reg31+reg147;
    reg145=reg145/reg153; reg73=reg102*reg84; reg156=reg156-reg60; reg98=reg2*reg64; reg103=reg40*reg34;
    reg110=reg20*reg6; reg177=reg33-reg177; reg91=reg91/reg39; reg66=reg166+reg66; reg129=reg133+reg129;
    reg88=reg26-reg88; reg141=reg141/reg153; reg26=reg54*reg4; reg97=reg111+reg97; reg33=reg95*reg42;
    reg57=reg172-reg57; reg166=reg166-reg106; reg113=reg65*reg124; reg169=reg168-reg169; reg114=reg5*reg76;
    reg151=reg35+reg151; reg92=reg25+reg92; reg25=reg6*reg123; reg35=reg76*reg78; reg150=reg23+reg150;
    reg126=reg7-reg126; reg7=reg6*reg135; reg25=reg28-reg25; reg23=reg58*reg46; reg28=reg86*reg91;
    reg115=reg120*reg148; reg73=reg98-reg73; reg164=reg164-reg60; reg97=reg13+reg97; reg109=reg139+reg109;
    reg99=0.5*reg99; reg49=reg49/reg153; reg98=reg107*reg136; reg32=reg24+reg32; reg116=reg159+reg116;
    reg35=reg171+reg35; reg0=reg53-reg0; reg165=reg165-reg60; reg177=reg177/reg153; reg129=0.5*reg129;
    reg33=reg57+reg33; reg113=reg169-reg113; reg92=reg151+reg92; reg53=reg54*reg138; reg57=reg40*reg138;
    reg130=reg54*reg160; reg101=reg140+reg101; reg127=reg70-reg127; reg70=reg145*reg147; reg133=reg72*reg4;
    reg30=reg50+reg30; reg50=reg72*reg156; reg166=reg119+reg166; reg26=reg103+reg26; reg131=reg131/reg153;
    reg119=reg149*reg9; reg139=reg40*reg160; reg137=reg45+reg137; reg88=reg114+reg88; reg55=reg48+reg55;
    reg104=reg104-reg60; reg45=reg66*reg80; reg5=reg5*reg69; reg117=reg81-reg117; reg82=reg1-reg82;
    reg89=reg86*reg89; reg36=reg122/reg36; reg64=reg41*reg64; reg51=reg31+reg51; reg27=reg68+reg27;
    reg87=reg112+reg87; reg1=reg54*reg118; reg8=reg8/reg153; reg84=reg132*reg84; reg14=reg29-reg14;
    reg146=reg146+reg56; reg13=reg174+reg13; reg10=reg158-reg10; reg29=reg58*reg12; reg31=reg6*reg37;
    reg48=reg2*reg94; reg68=reg102*reg11; reg61=reg41*reg61; reg93=reg132*reg93; reg43=reg43*reg102;
    reg20=reg20*reg2; reg24=reg24-reg128; reg81=reg54*reg75; reg173=reg56+reg173; reg16=reg111+reg16;
    reg56=reg150*reg141; reg42=reg77*reg42; reg67=reg62-reg67; reg3=reg110-reg3; reg65=reg65*reg155;
    reg175=reg175+reg19; reg24=reg19+reg24; reg119=reg45-reg119; reg19=reg131*reg166; reg56=reg70-reg56;
    reg11=reg132*reg11; reg94=reg41*reg94; reg58=reg58*reg52; reg6=reg6*reg38; reg45=reg69*reg78;
    reg68=reg48-reg68; reg48=reg76*reg71; reg31=reg29-reg31; reg10=reg10-reg60; reg29=reg49*reg175;
    reg153=reg14/reg153; reg14=reg8*reg51; reg62=reg32*reg177; reg35=reg35-reg60; reg28=reg0+reg28;
    reg97=reg22*reg97; reg0=reg54*reg165; reg22=reg40*reg165; reg70=reg54*reg164; reg110=reg36*reg99;
    reg111=reg40*reg164; reg27=reg27-reg60; reg116=0.5*reg116; reg89=reg82-reg89; reg126=reg109+reg126;
    reg57=reg130+reg57; reg1=reg87+reg1; reg146=reg15*reg146; reg173=reg142*reg173; reg15=reg36*reg129;
    reg124=reg77*reg124; reg117=reg5+reg117; reg33=reg113+reg33; reg55=reg88+reg55; reg3=reg65+reg3;
    reg42=reg67+reg42; reg137=0.5*reg137; reg53=reg139+reg53; reg26=reg50+reg26; reg50=reg30+reg50;
    reg5=reg54*reg156; reg133=reg103+reg133; reg127=reg101+reg127; reg155=reg95*reg155; reg92=0.5*reg92;
    reg30=reg72*reg104; reg43=reg20-reg43; reg61=reg93-reg61; reg136=reg105*reg136; reg84=reg64-reg84;
    reg46=reg102*reg46; reg135=reg2*reg135; reg23=reg7-reg23; reg120=reg120*reg96; reg7=reg72*reg138;
    reg121=reg132*reg121; reg81=reg16+reg81; reg123=reg41*reg123; reg13=reg134*reg13; reg115=reg25-reg115;
    reg98=reg73+reg98; reg127=0.5*reg127; reg136=reg84+reg136; reg173=reg146+reg173; reg52=reg102*reg52;
    reg7=reg139+reg7; reg124=reg61-reg124; reg48=reg31-reg48; reg78=reg86*reg78; reg38=reg2*reg38;
    reg43=reg155+reg43; reg96=reg107*reg96; reg45=reg68+reg45; reg19=reg119+reg19; elem.epsilon[0][0]=reg19;
    reg26=reg4*reg26; reg5=reg133+reg5; reg11=reg94-reg11; reg15=2*reg15; reg2=reg153*reg24;
    reg50=reg34*reg50; reg12=reg132*reg12; reg4=reg36*reg92; reg148=reg105*reg148; reg76=reg76*reg91;
    reg37=reg41*reg37; reg16=reg54*reg104; reg58=reg6-reg58; reg57=reg57+reg30; reg123=reg121-reg123;
    reg28=reg28-reg60; reg97=reg13+reg97; reg0=reg111+reg0; reg33=0.5*reg33; reg6=reg72*reg27;
    reg22=reg70+reg22; reg13=reg40*reg35; reg98=reg115+reg98; reg81=reg75*reg81; reg110=2*reg110;
    reg55=0.5*reg55; reg20=reg72*reg165; reg46=reg135-reg46; reg89=reg117+reg89; reg29=reg62-reg29;
    reg1=reg118*reg1; reg25=reg36*reg137; reg31=reg54*reg10; reg14=reg56-reg14; elem.epsilon[0][1]=reg14;
    reg53=reg30+reg53; reg30=reg40*reg10; reg34=reg54*reg35; reg23=reg120+reg23; reg126=0.5*reg126;
    reg42=reg3+reg42; reg3=reg36*reg116; reg41=reg36*reg55; reg5=reg156*reg5; reg56=reg36*reg33;
    reg42=0.5*reg42; reg25=2*reg25; reg71=reg86*reg71; reg26=reg50+reg26; reg37=reg12-reg37;
    reg2=reg29+reg2; elem.epsilon[0][2]=reg2; reg12=reg19+reg14; reg89=0.5*reg89; reg53=reg138*reg53;
    reg0=reg6+reg0; reg97=reg81+reg97; reg6=reg22+reg6; reg98=0.5*reg98; reg22=reg54*reg27;
    reg20=reg111+reg20; reg110=reg99*reg110; reg16=reg7+reg16; reg3=2*reg3; reg7=reg72*reg10;
    reg29=reg36*reg126; reg30=reg34+reg30; reg34=reg72*reg28; reg31=reg13+reg31; reg173=reg1+reg173;
    reg136=reg23+reg136; reg124=reg43+reg124; reg45=reg48+reg45; reg15=reg129*reg15; reg46=reg96+reg46;
    reg4=2*reg4; reg58=reg76+reg58; reg78=reg11+reg78; reg1=reg36*reg127; reg57=reg160*reg57;
    reg148=reg123-reg148; reg91=reg69*reg91; reg52=reg38-reg52; reg11=reg36*reg89; reg41=2*reg41;
    reg25=reg137*reg25; reg26=reg5+reg26; reg1=2*reg1; reg4=reg92*reg4; reg15=reg173+reg15;
    reg29=2*reg29; reg3=reg116*reg3; reg110=reg97+reg110; reg5=reg149*reg145; reg23=reg66*reg141;
    reg38=reg80*reg150; reg43=reg9*reg147; reg12=reg2+reg12; reg0=reg165*reg0; reg6=reg164*reg6;
    reg22=reg20+reg22; reg7=reg13+reg7; reg13=reg54*reg28; reg30=reg30+reg34; reg31=reg34+reg31;
    reg124=0.5*reg124; reg45=0.5*reg45; reg20=reg36*reg98; reg78=reg58+reg78; reg52=reg91+reg52;
    reg71=reg37-reg71; reg148=reg46+reg148; reg34=reg36*reg42; reg53=reg57+reg53; reg136=0.5*reg136;
    reg16=reg104*reg16; reg56=2*reg56; reg34=2*reg34; reg4=reg15+reg4; reg9=reg9*reg175;
    reg80=reg80*reg32; reg41=reg55*reg41; reg15=reg36*reg45; reg37=reg36*reg124; reg31=reg10*reg31;
    reg30=reg35*reg30; reg149=reg149*reg49; reg66=reg66*reg177; reg29=reg126*reg29; reg13=reg7+reg13;
    reg7=reg36*reg136; reg3=reg110+reg3; reg10=reg131*reg51; reg43=reg38-reg43; reg0=reg6+reg0;
    reg11=2*reg11; reg22=reg27*reg22; reg25=reg26+reg25; reg12=reg12/3; reg148=0.5*reg148;
    reg71=reg52+reg71; reg56=reg33*reg56; reg23=reg5-reg23; reg5=reg166*reg8; reg20=2*reg20;
    reg53=reg16+reg53; reg78=0.5*reg78; reg1=reg127*reg1; reg7=2*reg7; reg11=reg89*reg11;
    reg6=reg19-reg12; reg5=reg23-reg5; reg16=reg36*reg148; reg56=reg53+reg56; reg23=reg36*reg78;
    reg37=2*reg37; reg20=reg98*reg20; reg26=reg14-reg12; reg32=reg141*reg32; reg13=reg28*reg13;
    reg10=reg43+reg10; reg0=reg22+reg0; reg166=reg166*reg153; reg71=0.5*reg71; reg34=reg42*reg34;
    reg29=reg3+reg29; reg175=reg145*reg175; reg1=reg4+reg1; reg49=reg147*reg49; reg177=reg150*reg177;
    reg31=reg30+reg31; reg149=reg66-reg149; reg9=reg80-reg9; reg15=2*reg15; reg41=reg25+reg41;
    reg131=reg131*reg24; reg7=reg136*reg7; reg16=2*reg16; reg10=reg5+reg10; reg149=reg166+reg149;
    reg131=reg9+reg131; reg153=reg51*reg153; reg49=reg177-reg49; reg32=reg175-reg32; reg24=reg8*reg24;
    reg6=pow(reg6,2); reg26=pow(reg26,2); reg12=reg2-reg12; reg20=reg0+reg20; reg11=reg41+reg11;
    reg37=reg124*reg37; reg29=reg79*reg29; reg0=reg36*reg71; reg31=reg13+reg31; reg23=2*reg23;
    reg15=reg45*reg15; reg1=reg59*reg1; reg34=reg56+reg34; reg37=reg34+reg37; reg3=0.5*reg10;
    elem.epsilon[0][3]=reg3; reg29=0.083333333333333328707*reg29; reg15=reg31+reg15; reg1=0.083333333333333328707*reg1; reg131=reg149+reg131;
    reg23=reg78*reg23; reg49=reg153+reg49; reg0=2*reg0; reg24=reg32-reg24; reg11=reg100*reg11;
    reg26=reg6+reg26; reg12=pow(reg12,2); reg20=reg7+reg20; reg16=reg148*reg16; reg11=0.083333333333333328707*reg11;
    reg24=reg49+reg24; reg0=reg71*reg0; reg23=reg15+reg23; reg4=0.5*reg131; elem.epsilon[0][4]=reg4;
    reg12=reg26+reg12; reg1=reg29+reg1; reg10=reg10*reg3; reg37=reg21*reg37; reg16=reg20+reg16;
    reg11=reg1+reg11; reg131=reg131*reg4; reg10=reg12+reg10; reg0=reg23+reg0; reg1=0.5*reg24;
    elem.epsilon[0][5]=reg1; reg37=0.083333333333333328707*reg37; reg16=reg83*reg16; reg0=reg39*reg0; reg16=0.083333333333333328707*reg16;
    reg19=reg19-reg60; reg37=reg11+reg37; reg14=reg14-reg60; reg24=reg24*reg1; reg131=reg10+reg131;
    reg5=reg54*reg19; reg6=reg40*reg14; reg16=reg37+reg16; reg2=reg2-reg60; reg19=reg40*reg19;
    reg7=reg54*reg14; reg24=reg131+reg24; reg0=0.083333333333333328707*reg0; reg14=reg72*reg14; reg24=1.5*reg24;
    reg0=reg16+reg0; reg8=reg72*reg2; reg7=reg19+reg7; reg14=reg19+reg14; reg2=reg54*reg2;
    reg6=reg5+reg6; elem.sigma[0][5]=reg36*reg1; elem.sigma[0][0]=reg6+reg8; elem.sigma[0][4]=reg36*reg4; elem.sigma[0][1]=reg8+reg7;
    elem.sigma[0][2]=reg14+reg2; elem.ener=reg0/2; elem.sigma[0][3]=reg36*reg3; elem.sigma_von_mises=pow(reg24,0.5);
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1-var_inter[0]; reg0=reg0-var_inter[1]; T reg1=var_inter[0]*elem.pos(1)[2]; T reg2=reg0*elem.pos(0)[1]; T reg3=var_inter[0]*elem.pos(1)[1];
    T reg4=reg0*elem.pos(0)[2]; T reg5=1-var_inter[2]; T reg6=reg3+reg2; T reg7=var_inter[1]*elem.pos(2)[1]; T reg8=var_inter[1]*elem.pos(2)[2];
    T reg9=reg1+reg4; T reg10=elem.pos(3)[2]*reg0; T reg11=reg9+reg8; T reg12=reg5*elem.pos(1)[1]; T reg13=reg5*elem.pos(0)[1];
    T reg14=reg5*elem.pos(1)[2]; T reg15=reg5*elem.pos(0)[2]; T reg16=reg5*elem.pos(2)[1]; T reg17=reg5*elem.pos(2)[2]; T reg18=reg7+reg6;
    T reg19=elem.pos(3)[1]*reg0; reg17=reg17-reg15; T reg20=var_inter[0]*elem.pos(1)[0]; T reg21=elem.pos(0)[0]*reg0; reg12=reg12-reg13;
    T reg22=var_inter[2]*elem.pos(3)[1]; reg14=reg14-reg15; reg19=reg19-reg18; T reg23=var_inter[0]*elem.pos(4)[1]; T reg24=var_inter[2]*elem.pos(3)[2];
    reg16=reg16-reg13; reg10=reg10-reg11; T reg25=var_inter[0]*elem.pos(4)[2]; reg14=reg14-reg24; T reg26=var_inter[2]*elem.pos(4)[2];
    reg16=reg16-reg22; T reg27=var_inter[2]*elem.pos(4)[1]; reg12=reg12-reg22; T reg28=var_inter[2]*elem.pos(5)[2]; T reg29=reg5*elem.pos(2)[0];
    T reg30=var_inter[2]*elem.pos(5)[1]; reg25=reg10+reg25; reg10=var_inter[1]*elem.pos(5)[2]; reg23=reg19+reg23; reg19=var_inter[1]*elem.pos(5)[1];
    T reg31=1+(*f.m).poisson_ratio; T reg32=var_inter[1]*elem.pos(2)[0]; T reg33=reg20+reg21; T reg34=reg5*elem.pos(1)[0]; T reg35=reg5*elem.pos(0)[0];
    reg17=reg17-reg24; T reg36=reg33+reg32; T reg37=elem.pos(3)[0]*reg0; reg16=reg30+reg16; reg23=reg19+reg23;
    reg17=reg28+reg17; reg25=reg10+reg25; reg31=reg31/(*f.m).elastic_modulus; reg34=reg34-reg35; reg10=var_inter[2]*elem.pos(3)[0];
    reg27=reg12+reg27; reg26=reg14+reg26; reg29=reg29-reg35; reg12=var_inter[0]*elem.pos(4)[0]; reg37=reg37-reg36;
    reg14=pow(reg31,2); reg19=reg17*reg23; reg29=reg29-reg10; reg28=reg27*reg25; reg34=reg34-reg10;
    reg30=var_inter[2]*elem.pos(4)[0]; T reg38=reg16*reg25; T reg39=reg26*reg23; T reg40=var_inter[2]*elem.pos(5)[0]; T reg41=reg26*reg16;
    reg19=reg38-reg19; reg31=reg31*reg14; reg29=reg40+reg29; reg12=reg37+reg12; reg30=reg34+reg30;
    reg34=reg27*reg17; reg37=(*f.m).poisson_ratio/(*f.m).elastic_modulus; reg38=var_inter[1]*elem.pos(5)[0]; reg40=1.0/(*f.m).elastic_modulus; reg39=reg28-reg39;
    reg28=reg30*reg19; T reg42=reg40*reg14; reg41=reg34-reg41; reg34=reg37*reg31; reg31=reg40*reg31;
    T reg43=reg29*reg39; reg14=reg37*reg14; reg12=reg38+reg12; reg38=reg12*reg41; T reg44=reg16*reg12;
    T reg45=reg30*reg25; reg43=reg28-reg43; reg25=reg29*reg25; reg28=reg17*reg12; T reg46=reg29*reg23;
    T reg47=reg37*reg42; T reg48=reg37*reg14; reg42=reg40*reg42; T reg49=reg37*reg31; T reg50=reg37*reg34;
    reg23=reg30*reg23; T reg51=reg26*reg29; reg17=reg30*reg17; T reg52=reg27*reg12; reg31=reg40*reg31;
    reg12=reg26*reg12; reg49=reg50+reg49; reg31=reg31-reg50; reg28=reg25-reg28; reg44=reg46-reg44;
    reg34=reg40*reg34; reg14=reg40*reg14; reg42=reg42-reg48; reg47=reg48+reg47; reg38=reg43+reg38;
    reg29=reg27*reg29; reg51=reg17-reg51; reg16=reg30*reg16; reg12=reg45-reg12; reg52=reg23-reg52;
    reg42=reg40*reg42; reg47=reg37*reg47; reg34=reg50+reg34; reg17=reg48+reg14; reg40=reg40*reg31;
    reg12=reg12/reg38; reg52=reg52/reg38; reg41=reg41/reg38; reg28=reg28/reg38; reg23=reg37*reg49;
    reg51=reg51/reg38; reg19=reg19/reg38; reg39=reg39/reg38; reg44=reg44/reg38; reg29=reg16-reg29;
    reg47=reg42-reg47; reg16=var_inter[2]*reg52; reg17=reg37*reg17; reg25=reg5*reg19; reg26=reg5*reg39;
    reg27=var_inter[1]*reg41; reg30=var_inter[1]*reg51; reg42=reg5*reg12; reg43=reg5*reg28; reg45=reg5*reg44;
    reg46=reg5*reg52; reg50=var_inter[0]*reg51; T reg53=var_inter[0]*reg41; reg29=reg29/reg38; reg23=reg40-reg23;
    reg37=reg37*reg34; reg40=var_inter[2]*reg28; T reg54=var_inter[2]*reg12; T reg55=var_inter[2]*reg39; T reg56=var_inter[2]*reg44;
    T reg57=var_inter[2]*reg19; reg37=reg23-reg37; reg23=reg46-reg45; T reg58=reg0*reg29; T reg59=reg0*reg51;
    T reg60=reg43-reg42; T reg61=reg40-reg54; T reg62=reg57+reg53; T reg63=reg40+reg50; T reg64=reg0*reg41;
    T reg65=reg26-reg25; T reg66=reg16-reg56; T reg67=reg55-reg57; reg17=reg47-reg17; reg47=reg26+reg27;
    T reg68=var_inter[0]*reg29; T reg69=reg30+reg42; T reg70=var_inter[1]*reg29; T reg71=0.5*reg63; reg65=reg65-reg64;
    T reg72=reg45-reg68; T reg73=reg50-reg43; T reg74=reg68+reg56; reg67=reg67+reg64; reg60=reg60+reg59;
    reg66=reg66+reg58; reg61=reg61-reg59; T reg75=0.5*reg69; reg23=reg23-reg58; T reg76=reg70+reg46;
    T reg77=reg25-reg53; T reg78=0.5*reg47; reg17=reg17/reg37; T reg79=(*f.m).deltaT*(*f.m).alpha; T reg80=reg70-reg16;
    T reg81=reg54-reg30; T reg82=reg27-reg55; T reg83=0.5*reg62; reg34=reg34/reg37; reg31=reg31/reg37;
    reg37=reg49/reg37; reg49=reg17*reg78; T reg84=0.5*reg67; T reg85=0.5*reg81; T reg86=0.5*reg82;
    T reg87=0.5*reg74; T reg88=0.5*reg80; T reg89=0.5*reg60; T reg90=0.5*reg76; T reg91=reg17*reg83;
    T reg92=reg31*reg79; T reg93=0.5*reg23; T reg94=reg17*reg71; T reg95=reg17*reg75; T reg96=0.5*reg61;
    T reg97=0.5*reg77; T reg98=0.5*reg72; T reg99=0.5*reg73; T reg100=reg34*reg79; T reg101=0.5*reg66;
    T reg102=0.5*reg65; T reg103=reg37*reg79; T reg104=reg17*reg88; T reg105=reg31*reg76; T reg106=reg17*reg86;
    T reg107=2*reg94; T reg108=reg17*reg84; T reg109=reg17*reg89; T reg110=reg31*reg63; T reg111=reg17*reg101;
    T reg112=reg31*reg62; T reg113=reg17*reg87; T reg114=reg17*reg93; T reg115=reg17*reg96; T reg116=reg92+reg103;
    reg91=2*reg91; T reg117=reg17*reg85; T reg118=reg17*reg99; T reg119=reg31*reg74; T reg120=reg17*reg90;
    T reg121=reg17*reg98; T reg122=reg17*reg102; T reg123=reg31*reg47; reg95=2*reg95; T reg124=reg17*reg97;
    T reg125=2*reg49; T reg126=reg31*reg69; T reg127=reg103+reg100; reg108=2*reg108; T reg128=reg76*reg119;
    T reg129=reg71*reg95; T reg130=reg62*reg123; reg115=2*reg115; T reg131=reg34*reg74; T reg132=reg31*reg67;
    T reg133=reg34*reg66; T reg134=reg34*reg76; T reg135=reg34*reg72; reg111=2*reg111; T reg136=reg31*reg82;
    T reg137=2*reg120; reg114=2*reg114; reg117=2*reg117; reg109=2*reg109; T reg138=reg31*reg80;
    T reg139=reg37*reg63; T reg140=reg31*reg65; T reg141=reg37*reg69; T reg142=reg34*reg63; reg113=2*reg113;
    reg124=2*reg124; T reg143=reg31*reg60; T reg144=reg31*reg72; T reg145=reg31*reg61; reg121=2*reg121;
    T reg146=reg83*reg125; reg118=2*reg118; T reg147=reg31*reg77; T reg148=reg116+reg100; T reg149=reg126*reg63;
    T reg150=reg31*reg23; T reg151=reg34*reg23; T reg152=reg74*reg105; T reg153=reg92+reg127; reg122=2*reg122;
    T reg154=reg37*reg77; T reg155=reg37*reg67; T reg156=reg37*reg47; T reg157=reg37*reg82; T reg158=reg47*reg112;
    T reg159=reg75*reg107; T reg160=reg31*reg81; T reg161=reg31*reg73; reg104=2*reg104; reg106=2*reg106;
    T reg162=reg37*reg65; T reg163=var_inter[1]*reg5; T reg164=reg34*reg69; T reg165=var_inter[0]*var_inter[2]; T reg166=reg31*reg66;
    T reg167=reg37*reg62; T reg168=reg34*reg80; T reg169=reg69*reg110; T reg170=reg78*reg91; T reg171=reg66*reg166;
    T reg172=reg66*reg119; T reg173=reg99*reg109; T reg174=reg66*reg138; T reg175=reg74*reg166; T reg176=reg98*reg125;
    T reg177=reg62*reg140; T reg178=reg61*reg160; T reg179=reg60*reg160; T reg180=reg23*reg138; T reg181=reg74*reg144;
    T reg182=reg60*reg110; T reg183=reg156*reg74; T reg184=reg83*reg137; T reg185=reg77*reg134; T reg186=reg102*reg91;
    T reg187=reg84*reg106; T reg188=reg147*reg77; T reg189=reg146+reg152; T reg190=reg66*reg150; T reg191=reg77*reg123;
    T reg192=reg66*reg105; T reg193=reg84*reg137; T reg194=reg99*reg95; T reg195=reg66*reg144; T reg196=reg156*reg66;
    T reg197=reg34*reg73; reg149=reg146+reg149; T reg198=reg23*reg144; T reg199=reg63*reg161; T reg200=reg83*reg124;
    T reg201=reg156*reg23; T reg202=reg102*reg137; T reg203=reg63*reg143; T reg204=reg83*reg122; T reg205=reg23*reg105;
    T reg206=reg71*reg91; T reg207=reg62*reg139; T reg208=reg34*reg61; T reg209=reg71*reg107; T reg210=reg62*reg134;
    T reg211=reg62*reg112; T reg212=reg23*reg166; T reg213=reg71*reg115; T reg214=reg62*reg132; T reg215=reg87*reg125;
    T reg216=reg102*reg106; T reg217=reg74*reg150; T reg218=reg63*reg160; T reg219=reg83*reg106; T reg220=reg71*reg109;
    T reg221=reg34*reg81; T reg222=reg34*reg60; T reg223=reg87*reg107; T reg224=reg147*reg62; T reg225=reg118*reg71;
    T reg226=reg23*reg119; T reg227=reg63*reg131; T reg228=reg63*reg110; T reg229=reg83*reg91; T reg230=reg130+reg129;
    T reg231=reg63*reg167; T reg232=reg87*reg137; T reg233=reg83*reg107; T reg234=reg63*reg145; T reg235=reg83*reg108;
    T reg236=reg47*reg168; T reg237=reg90*reg106; T reg238=reg69*reg143; T reg239=reg78*reg122; T reg240=reg72*reg105;
    T reg241=reg69*reg161; T reg242=reg78*reg124; T reg243=reg97*reg137; T reg244=reg156*reg72; T reg245=reg156*reg69;
    T reg246=reg78*reg95; T reg247=reg72*reg144; T reg248=reg126*reg69; T reg249=reg78*reg125; T reg250=reg69*reg134;
    T reg251=reg90*reg95; T reg252=reg69*reg145; T reg253=reg78*reg108; T reg254=reg72*reg150; T reg255=reg169+reg170;
    T reg256=reg69*reg160; T reg257=reg78*reg106; T reg258=reg76*reg162; T reg259=reg78*reg114; T reg260=reg97*reg106;
    T reg261=reg73*reg160; T reg262=reg76*reg150; T reg263=reg151*reg47; T reg264=reg122*reg90; T reg265=reg118*reg75;
    T reg266=reg147*reg47; T reg267=reg47*reg140; T reg268=reg75*reg109; T reg269=reg72*reg138; T reg270=reg47*reg135;
    T reg271=reg124*reg90; T reg272=reg75*reg95; T reg273=reg47*reg123; T reg274=reg75*reg125; T reg275=reg47*reg141;
    T reg276=reg75*reg115; T reg277=reg72*reg119; T reg278=reg47*reg132; T reg279=reg47*reg133; T reg280=reg90*reg108;
    T reg281=reg71*reg117; T reg282=reg62*reg136; reg158=reg159+reg158; T reg283=reg90*reg113; T reg284=reg47*reg131;
    T reg285=reg90*reg91; T reg286=reg72*reg166; T reg287=reg75*reg117; T reg288=reg47*reg136; T reg289=reg124*reg97;
    T reg290=reg73*reg161; T reg291=reg67*reg134; T reg292=reg101*reg125; T reg293=reg67*reg132; T reg294=reg96*reg115;
    T reg295=reg122*reg97; T reg296=reg67*reg112; T reg297=reg96*reg107; T reg298=reg67*reg136; T reg299=reg96*reg117;
    T reg300=reg61*reg143; T reg301=reg122*reg84; T reg302=reg77*reg136; T reg303=reg99*reg117; T reg304=reg61*reg161;
    T reg305=reg124*reg84; T reg306=reg126*reg61; T reg307=reg84*reg125; T reg308=reg77*reg112; T reg309=reg99*reg107;
    T reg310=reg61*reg145; T reg311=reg84*reg108; T reg312=reg61*reg110; T reg313=reg84*reg91; T reg314=reg77*reg132;
    T reg315=reg99*reg115; T reg316=reg76*reg154; T reg317=reg78*reg121; T reg318=reg76*reg144; T reg319=reg75*reg137;
    T reg320=reg76*reg164; T reg321=reg97*reg91; T reg322=reg73*reg110; T reg323=reg76*reg105; T reg324=reg76*reg155;
    T reg325=reg78*reg111; T reg326=reg76*reg166; T reg327=reg76*reg167; T reg328=reg78*reg113; T reg329=reg97*reg108;
    T reg330=reg73*reg145; reg128=reg170+reg128; reg170=reg76*reg157; T reg331=reg78*reg104; T reg332=reg76*reg138;
    T reg333=reg67*reg140; T reg334=reg96*reg109; T reg335=reg97*reg125; T reg336=reg126*reg73; T reg337=reg147*reg67;
    T reg338=reg118*reg96; T reg339=reg67*reg123; T reg340=reg96*reg95; T reg341=reg86*reg91; T reg342=reg81*reg110;
    T reg343=reg37*reg81; T reg344=reg65*reg140; T reg345=reg89*reg109; T reg346=reg86*reg108; T reg347=reg81*reg145;
    T reg348=reg65*reg134; T reg349=reg93*reg125; T reg350=reg65*reg132; T reg351=reg89*reg115; T reg352=reg125*reg86;
    T reg353=reg126*reg81; T reg354=reg37*reg61; T reg355=reg124*reg86; T reg356=reg81*reg161; T reg357=reg60*reg143;
    T reg358=reg82*reg123; T reg359=reg85*reg95; T reg360=reg82*reg134; T reg361=reg88*reg125; reg132=reg82*reg132;
    T reg362=reg85*reg115; T reg363=reg82*reg112; T reg364=reg85*reg107; T reg365=reg89*reg117; T reg366=reg65*reg136;
    reg136=reg82*reg136; T reg367=reg85*reg117; T reg368=reg81*reg143; T reg369=reg89*reg107; reg112=reg65*reg112;
    T reg370=reg122*reg86; T reg371=reg80*reg138; T reg372=reg118*reg89; T reg373=reg147*reg65; T reg374=reg37*reg60;
    T reg375=reg47*reg148; T reg376=reg76*reg153; T reg377=reg163*(*f.m).f_vol[2]; T reg378=reg165*(*f.m).f_vol[1]; T reg379=reg63*reg148;
    T reg380=reg102*reg125; T reg381=reg118*reg99; reg126=reg126*reg60; reg161=reg60*reg161; T reg382=var_inter[2]*reg0;
    T reg383=var_inter[1]*var_inter[2]; T reg384=reg5*var_inter[0]; T reg385=reg5*reg0; reg160=reg81*reg160; T reg386=reg89*reg95;
    T reg387=reg65*reg123; T reg388=reg86*reg106; T reg389=reg23*reg150; reg150=reg80*reg150; T reg390=reg163*(*f.m).f_vol[0];
    reg144=reg80*reg144; T reg391=reg156*reg80; T reg392=reg137*reg86; reg145=reg60*reg145; T reg393=reg80*reg105;
    reg143=reg73*reg143; T reg394=reg77*reg140; reg166=reg80*reg166; T reg395=reg80*reg119; T reg396=reg37*reg73;
    T reg397=reg85*reg118; reg147=reg82*reg147; reg138=reg74*reg138; T reg398=reg102*reg108; reg119=reg74*reg119;
    reg140=reg82*reg140; T reg399=reg85*reg109; T reg400=reg124*reg102; T reg401=reg71*reg113; T reg402=reg74*reg142;
    T reg403=reg122*reg102; reg138=reg219+reg138; T reg404=reg121*reg71; reg317=reg316+reg317; reg316=reg121*reg75;
    T reg405=reg76*reg197; T reg406=reg74*reg197; reg318=reg242+reg318; T reg407=reg156*reg76; T reg408=reg78*reg137;
    reg144=reg355+reg144; T reg409=reg83*reg121; reg320=reg319+reg320; T reg410=reg85*reg121; T reg411=reg80*reg197;
    T reg412=reg121*reg86; T reg413=reg249+reg323; T reg414=reg74*reg154; T reg415=reg80*reg154; reg325=reg324+reg325;
    reg324=reg75*reg111; T reg416=reg76*reg208; T reg417=reg86*reg111; T reg418=reg80*reg155; T reg419=reg283+reg255;
    T reg420=reg69*reg131; T reg421=reg90*reg107; T reg422=reg69*reg157; T reg423=reg78*reg117; T reg424=reg352+reg393;
    T reg425=reg71*reg104; reg256=reg256-reg257; T reg426=reg69*reg168; T reg427=reg90*reg117; reg181=reg200+reg181;
    T reg428=reg81*reg167; reg259=reg258+reg259; reg258=reg75*reg114; T reg429=reg76*reg222; T reg430=reg85*reg137;
    T reg431=reg80*reg164; T reg432=reg391+reg392; reg262=reg239+reg262; T reg433=reg67*reg141; reg333=reg333+reg334;
    T reg434=reg101*reg114; reg140=reg140+reg399; T reg435=reg71*reg114; T reg436=reg122*reg96; T reg437=reg151*reg67;
    T reg438=reg122*reg101; T reg439=reg86*reg117; T reg440=reg81*reg157; reg337=reg337+reg338; T reg441=reg121*reg101;
    T reg442=reg396*reg67; T reg443=reg124*reg96; T reg444=reg67*reg135; T reg445=reg124*reg101; T reg446=reg88*reg107;
    T reg447=reg81*reg131; T reg448=reg101*reg137; T reg449=reg340-reg339; T reg450=reg341-reg342; T reg451=reg74*reg222;
    reg150=reg370+reg150; reg326=reg253+reg326; reg217=reg204+reg217; T reg452=reg85*reg114; reg328=reg327+reg328;
    reg327=reg75*reg113; T reg453=reg291+reg292; T reg454=reg76*reg142; T reg455=reg80*reg222; T reg456=reg86*reg114;
    T reg457=reg80*reg162; T reg458=reg88*reg117; reg128=reg159+reg128; T reg459=reg81*reg168; T reg460=reg86*reg107;
    T reg461=reg83*reg114; reg331=reg170+reg331; reg170=reg75*reg104; T reg462=reg76*reg221; reg160=reg160+reg388;
    reg332=reg257+reg332; reg257=reg96*reg125; T reg463=reg90*reg125; T reg464=reg74*reg155; T reg465=reg66*reg153;
    T reg466=reg61*reg148; reg278=reg276-reg278; T reg467=reg90*reg111; T reg468=reg75*reg108; T reg469=reg47*reg354;
    T reg470=reg67*reg148; reg129=reg129+reg189; reg280=reg279+reg280; T reg471=reg376-reg377; T reg472=reg69*reg148;
    T reg473=reg375-reg390; reg119=reg229+reg119; T reg474=reg72*reg153; reg283=reg158+reg283; T reg475=reg75*reg91;
    T reg476=reg47*reg139; T reg477=reg73*reg148; T reg478=reg77*reg148; reg285=reg284+reg285; T reg479=reg23*reg153;
    T reg480=reg83*reg113; reg264=reg263+reg264; T reg481=reg74*reg167; reg175=reg235+reg175; reg266=reg265-reg266;
    T reg482=reg80*reg153; T reg483=reg81*reg148; T reg484=reg121*reg90; T reg485=reg83*reg95; T reg486=reg156*reg63;
    reg401=reg402+reg401; T reg487=reg71*reg111; T reg488=reg74*reg208; T reg489=reg82*reg148; reg271=reg270+reg271;
    T reg490=reg74*reg153; T reg491=reg379-reg378; T reg492=reg272+reg273; T reg493=reg90*reg137; T reg494=reg83*reg111;
    T reg495=reg62*reg148; reg275=reg274+reg275; T reg496=reg47*reg134; T reg497=reg86*reg104; T reg498=reg118*reg90;
    T reg499=reg80*reg157; reg395=reg341+reg395; reg246=reg245+reg246; reg341=reg183+reg184; T reg500=reg85*reg113;
    T reg501=reg80*reg142; reg248=reg248+reg249; T reg502=reg86*reg113; T reg503=reg80*reg167; reg251=reg250+reg251;
    T reg504=reg69*reg155; T reg505=reg78*reg115; reg166=reg346+reg166; T reg506=reg74*reg221; reg253=reg252-reg253;
    reg252=reg69*reg133; T reg507=reg90*reg115; T reg508=reg69*reg167; T reg509=reg85*reg111; T reg510=reg78*reg107;
    T reg511=reg80*reg208; T reg512=reg60*reg148; reg288=reg287-reg288; T reg513=reg90*reg104; T reg514=reg75*reg106;
    T reg515=reg47*reg343; T reg516=reg65*reg148; T reg517=reg74*reg157; reg237=reg236+reg237; T reg518=reg69*reg162;
    T reg519=reg78*reg109; T reg520=reg83*reg104; reg371=reg388+reg371; reg388=reg71*reg137; reg239=reg238-reg239;
    reg238=reg151*reg69; T reg521=reg90*reg109; T reg522=reg69*reg154; T reg523=reg78*reg118; T reg524=reg85*reg104;
    T reg525=reg74*reg164; reg242=reg241-reg242; reg241=reg69*reg135; T reg526=reg80*reg221; reg177=reg177-reg220;
    T reg527=reg87*reg114; T reg528=reg374*reg62; T reg529=reg122*reg71; T reg530=reg88*reg104; T reg531=reg151*reg62;
    T reg532=reg122*reg87; reg136=reg136+reg367; T reg533=reg82*reg396; reg224=reg224-reg225; T reg534=reg121*reg87;
    T reg535=reg396*reg62; T reg536=reg124*reg71; T reg537=reg62*reg135; T reg538=reg87*reg115; T reg539=reg88*reg91;
    T reg540=reg124*reg87; T reg541=reg82*reg131; T reg542=reg85*reg91; T reg543=reg82*reg139; T reg544=reg63*reg133;
    reg234=reg235-reg234; reg235=reg230+reg232; T reg545=reg84*reg111; reg147=reg147+reg397; T reg546=reg66*reg208;
    T reg547=reg96*reg111; reg231=reg233+reg231; reg171=reg311+reg171; T reg548=reg66*reg167; T reg549=reg84*reg113;
    T reg550=reg66*reg142; T reg551=reg96*reg113; T reg552=reg88*reg121; reg172=reg313+reg172; T reg553=reg66*reg157;
    T reg554=reg84*reg104; T reg555=reg66*reg221; T reg556=reg96*reg104; T reg557=reg86*reg109; T reg558=reg81*reg162;
    reg174=reg187+reg174; T reg559=reg88*reg106; T reg560=reg82*reg168; T reg561=reg85*reg106; T reg562=reg82*reg343;
    T reg563=reg62*reg343; T reg564=reg71*reg106; T reg565=reg83*reg115; T reg566=reg360+reg361; T reg567=reg62*reg168;
    T reg568=reg87*reg106; T reg569=reg83*reg109; T reg570=reg87*reg95; T reg571=reg88*reg124; T reg572=reg85*reg125;
    reg203=reg204-reg203; reg204=reg151*reg63; T reg573=reg87*reg109; T reg574=reg83*reg118; T reg575=reg82*reg141;
    T reg576=reg88*reg137; T reg577=reg63*reg154; T reg578=reg359-reg358; T reg579=reg63*reg134; reg199=reg200-reg199;
    reg200=reg232+reg149; T reg580=reg63*reg135; T reg581=reg118*reg87; T reg582=reg62*reg141; T reg583=reg71*reg125;
    T reg584=reg88*reg113; reg363=reg363-reg364; T reg585=reg210+reg215; T reg586=reg85*reg124; T reg587=reg82*reg135;
    reg214=reg214-reg213; T reg588=reg87*reg111; T reg589=reg62*reg354; T reg590=reg88*reg108; T reg591=reg82*reg133;
    T reg592=reg71*reg108; T reg593=reg62*reg133; T reg594=reg87*reg108; T reg595=reg85*reg108; T reg596=reg82*reg354;
    T reg597=reg88*reg111; reg211=reg211+reg209; T reg598=reg87*reg113; reg132=reg132+reg362; T reg599=reg63*reg155;
    reg206=reg207+reg206; T reg600=reg62*reg131; T reg601=reg61*reg162; T reg602=reg84*reg109; reg218=reg219-reg218;
    reg219=reg88*reg114; T reg603=reg86*reg115; T reg604=reg82*reg374; reg300=reg300+reg301; T reg605=reg151*reg61;
    T reg606=reg101*reg109; T reg607=reg61*reg154; T reg608=reg118*reg84; T reg609=reg63*reg157; T reg610=reg81*reg155;
    reg304=reg304+reg305; T reg611=reg61*reg135; T reg612=reg118*reg101; T reg613=reg156*reg61; T reg614=reg84*reg95;
    T reg615=reg88*reg95; T reg616=reg81*reg134; reg306=reg306-reg307; reg353=reg353-reg352; T reg617=reg83*reg117;
    T reg618=reg88*reg115; reg293=reg293+reg294; T reg619=reg101*reg111; T reg620=reg67*reg354; T reg621=reg96*reg108;
    T reg622=reg67*reg133; T reg623=reg101*reg108; T reg624=reg81*reg133; T reg625=reg74*reg162; T reg626=reg87*reg117;
    reg296=reg296-reg297; T reg627=reg101*reg113; T reg628=reg67*reg139; T reg629=reg96*reg91; T reg630=reg67*reg131;
    T reg631=reg101*reg91; T reg632=reg63*reg168; reg346=reg347+reg346; reg298=reg298+reg299; reg347=reg101*reg104;
    T reg633=reg67*reg343; T reg634=reg96*reg106; T reg635=reg67*reg168; T reg636=reg101*reg106; reg187=reg178+reg187;
    reg178=reg61*reg168; T reg637=reg101*reg117; T reg638=reg66*reg162; T reg639=reg84*reg114; T reg640=reg66*reg222;
    T reg641=reg96*reg114; reg190=reg301+reg190; reg301=reg118*reg86; T reg642=reg81*reg154; T reg643=reg66*reg154;
    T reg644=reg121*reg84; T reg645=reg66*reg197; T reg646=reg121*reg96; T reg647=reg88*reg109; reg195=reg305+reg195;
    reg305=reg81*reg151; reg370=reg368+reg370; reg368=reg196+reg193; T reg648=reg66*reg164; T reg649=reg96*reg137;
    T reg650=reg307+reg192; T reg651=reg66*reg155; T reg652=reg61*reg134; T reg653=reg101*reg95; T reg654=reg61*reg155;
    T reg655=reg84*reg115; T reg656=reg227+reg223; T reg657=reg85*reg122; reg311=reg310+reg311; reg310=reg61*reg133;
    T reg658=reg101*reg115; T reg659=reg61*reg167; T reg660=reg84*reg107; T reg661=reg82*reg151; T reg662=reg95*reg86;
    reg313=reg313-reg312; T reg663=reg61*reg131; T reg664=reg156*reg81; T reg665=reg88*reg118; T reg666=reg101*reg107;
    T reg667=reg61*reg157; T reg668=reg84*reg117; T reg669=reg88*reg122; T reg670=reg81*reg135; reg355=reg356+reg355;
    reg229=reg229+reg228; reg356=reg151*reg77; T reg671=reg23*reg155; T reg672=reg65*reg133; reg179=reg179+reg216;
    T reg673=reg121*reg99; T reg674=reg72*reg197; T reg675=reg93*reg108; T reg676=reg102*reg111; T reg677=reg122*reg99;
    reg112=reg112-reg369; reg143=reg143+reg295; T reg678=reg93*reg113; T reg679=reg73*reg167; T reg680=reg97*reg107;
    T reg681=reg121*reg97; T reg682=reg72*reg154; T reg683=reg65*reg139; T reg684=reg89*reg91; T reg685=reg98*reg114;
    reg394=reg394+reg173; T reg686=reg65*reg131; T reg687=reg93*reg91; T reg688=reg102*reg117; T reg689=reg60*reg157;
    reg254=reg295+reg254; reg295=reg23*reg208; reg366=reg366+reg365; T reg690=reg77*reg168; T reg691=reg93*reg137;
    T reg692=reg98*reg106; T reg693=reg396*reg77; T reg694=reg65*reg141; T reg695=reg89*reg125; T reg696=reg124*reg99;
    T reg697=reg98*reg115; T reg698=reg99*reg137; T reg699=reg72*reg164; T reg700=reg383*(*f.m).f_vol[1]; T reg701=reg383*(*f.m).f_vol[2];
    T reg702=reg93*reg117; reg344=reg344+reg345; T reg703=reg93*reg114; T reg704=reg121*reg98; reg188=reg381+reg188;
    T reg705=reg348+reg349; T reg706=reg244+reg243; T reg707=reg73*reg162; T reg708=reg60*reg168; reg350=reg350+reg351;
    T reg709=reg93*reg111; reg247=reg289+reg247; T reg710=reg122*reg98; T reg711=reg97*reg109; T reg712=reg65*reg354;
    T reg713=reg89*reg108; T reg714=reg151*reg60; T reg715=reg93*reg109; T reg716=reg73*reg168; reg226=reg186+reg226;
    T reg717=reg60*reg154; T reg718=reg118*reg102; reg186=reg186-reg182; reg261=reg261+reg260; T reg719=reg89*reg113;
    reg289=reg290+reg289; reg290=reg23*reg142; T reg720=reg60*reg155; T reg721=reg102*reg115; T reg722=reg73*reg131;
    T reg723=reg102*reg113; T reg724=reg97*reg117; T reg725=reg73*reg157; reg145=reg145+reg398; T reg726=reg73*reg135;
    T reg727=reg60*reg133; T reg728=reg102*reg107; T reg729=reg93*reg115; T reg730=reg396*reg47; T reg731=reg124*reg75;
    T reg732=reg118*reg98; T reg733=reg98*reg107; T reg734=reg60*reg167; T reg735=reg93*reg104; T reg736=reg151*reg73;
    T reg737=reg98*reg109; T reg738=reg23*reg167; T reg739=reg99*reg114; T reg740=reg72*reg222; T reg741=reg97*reg95;
    T reg742=reg65*reg343; T reg743=reg89*reg106; T reg744=reg89*reg104; T reg745=reg23*reg221; reg168=reg65*reg168;
    T reg746=reg156*reg73; T reg747=reg93*reg106; T reg748=reg93*reg107; T reg749=reg97*reg114; T reg750=reg60*reg131;
    T reg751=reg60*reg162; T reg752=reg73*reg154; T reg753=reg72*reg162; reg109=reg102*reg109; T reg754=reg118*reg97;
    T reg755=reg102*reg104; reg357=reg357+reg403; T reg756=reg23*reg157; T reg757=reg321-reg322; reg117=reg98*reg117;
    T reg758=reg89*reg111; T reg759=reg99*reg108; reg221=reg72*reg221; reg131=reg77*reg131; T reg760=reg374*reg65;
    T reg761=reg98*reg91; T reg762=reg122*reg89; T reg763=reg73*reg134; T reg764=reg121*reg102; T reg765=reg97*reg104;
    reg157=reg72*reg157; T reg766=reg98*reg111; reg151=reg151*reg65; reg154=reg23*reg154; T reg767=reg122*reg93;
    reg314=reg315+reg314; reg373=reg373+reg372; T reg768=reg121*reg93; T reg769=reg201+reg202; reg277=reg321+reg277;
    reg396=reg396*reg65; reg321=reg185+reg176; T reg770=reg124*reg89; reg302=reg303+reg302; T reg771=reg98*reg104;
    T reg772=reg99*reg113; T reg773=reg72*reg142; T reg774=reg97*reg113; T reg775=reg382*(*f.m).f_vol[1]; T reg776=reg374*reg47;
    reg113=reg98*reg113; reg308=reg308-reg309; reg122=reg122*reg75; T reg777=reg63*reg162; reg121=reg121*reg89;
    T reg778=reg93*reg95; reg198=reg400+reg198; T reg779=reg99*reg91; T reg780=reg77*reg139; T reg781=reg90*reg114;
    reg267=reg268-reg267; reg197=reg23*reg197; T reg782=reg60*reg135; reg108=reg98*reg108; T reg783=reg165*(*f.m).f_vol[2];
    T reg784=reg383*(*f.m).f_vol[0]; T reg785=reg77*reg133; T reg786=reg163*(*f.m).f_vol[1]; T reg787=reg73*reg155; reg269=reg260+reg269;
    reg126=reg126-reg380; reg115=reg97*reg115; reg354=reg77*reg354; reg260=reg98*reg95; T reg788=reg99*reg104;
    reg95=reg102*reg95; T reg789=reg156*reg60; T reg790=reg374*reg67; reg180=reg216+reg180; reg216=reg194-reg191;
    T reg791=reg384*(*f.m).f_vol[0]; reg222=reg23*reg222; T reg792=reg99*reg111; reg208=reg72*reg208; T reg793=reg384*(*f.m).f_vol[1];
    T reg794=reg382*(*f.m).f_vol[0]; reg336=reg336-reg335; T reg795=reg384*(*f.m).f_vol[2]; T reg796=reg385*(*f.m).f_vol[0]; reg343=reg77*reg343;
    T reg797=reg385*(*f.m).f_vol[1]; T reg798=reg385*(*f.m).f_vol[2]; T reg799=reg89*reg137; reg111=reg97*reg111; reg155=reg72*reg155;
    T reg800=reg65*reg135; T reg801=reg124*reg98; reg124=reg124*reg93; reg133=reg73*reg133; reg135=reg77*reg135;
    reg403=reg389+reg403; reg162=reg23*reg162; reg389=reg386-reg387; T reg802=reg335+reg240; T reg803=reg99*reg125;
    reg330=reg330+reg329; T reg804=reg165*(*f.m).f_vol[0]; reg141=reg77*reg141; reg282=reg282-reg281; reg164=reg23*reg164;
    reg91=reg87*reg91; reg167=reg72*reg167; reg118=reg118*reg93; reg104=reg87*reg104; T reg805=reg60*reg134;
    T reg806=reg89*reg114; T reg807=reg382*(*f.m).f_vol[2]; reg374=reg374*reg77; T reg808=reg380+reg205; T reg809=reg98*reg137;
    reg114=reg102*reg114; reg400=reg161+reg400; reg106=reg99*reg106; reg286=reg329+reg286; reg212=reg398+reg212;
    reg161=reg38*reg769; reg777=reg569-reg777; reg164=reg164-reg799; reg199=reg534+reg199; reg214=reg214+reg588;
    reg211=reg211+reg598; reg676=reg671+reg676; reg568=reg567+reg568; reg351=reg212+reg351; reg203=reg527+reg203;
    reg592=reg589-reg592; reg573=reg573-reg204; reg198=reg372+reg198; reg577=reg574-reg577; reg758=reg295+reg758;
    reg594=reg593+reg594; reg564=reg563-reg564; reg386=reg386-reg808; reg212=reg38*reg206; reg188=reg188+reg704;
    reg648=reg648-reg649; reg295=reg38*reg368; reg693=reg696+reg693; reg195=reg338+reg195; reg646=reg645+reg646;
    reg801=reg135+reg801; reg644=reg643+reg644; reg190=reg334+reg190; reg641=reg640+reg641; reg216=reg216-reg809;
    reg639=reg638+reg639; reg637=reg178+reg637; reg141=reg141-reg803; reg187=reg347+reg187; reg668=reg667+reg668;
    reg135=reg38*reg321; reg663=reg663-reg666; reg313=reg627+reg313; reg314=reg314+reg766; reg659=reg659-reg660;
    reg658=reg310+reg658; reg354=reg759+reg354; reg311=reg619+reg311; reg655=reg654+reg655; reg108=reg785+reg108;
    reg653=reg653-reg652; reg306=reg306-reg448; reg308=reg308+reg113; reg178=reg38*reg585; reg730=reg731-reg730;
    reg582=reg582+reg583; reg723=reg738+reg723; reg310=reg38*reg235; reg719=reg719-reg290; reg540=reg537+reg540;
    reg536=reg535-reg536; reg226=reg226-reg369; reg534=reg224+reg534; reg755=reg756+reg755; reg532=reg531+reg532;
    reg529=reg528-reg529; reg527=reg177+reg527; reg744=reg745+reg744; reg365=reg180+reg365; reg174=reg299+reg174;
    reg556=reg555+reg556; reg554=reg553+reg554; reg394=reg394+reg685; reg172=reg172-reg297; reg551=reg551-reg550;
    reg677=reg374+reg677; reg549=reg548+reg549; reg171=reg294+reg171; reg710=reg356+reg710; reg547=reg546+reg547;
    reg545=reg651+reg545; reg340=reg340-reg650; reg694=reg694-reg695; reg450=reg584+reg450; reg447=reg447-reg446;
    reg389=reg389-reg691; reg439=reg440+reg439; reg345=reg403+reg345; reg160=reg530+reg160; reg458=reg459+reg458;
    reg124=reg800+reg124; reg456=reg457+reg456; reg452=reg455+reg452; reg150=reg399+reg150; reg412=reg415+reg412;
    reg410=reg411+reg410; reg144=reg397+reg144; reg177=reg38*reg432; reg431=reg431-reg430; reg359=reg359-reg424;
    reg557=reg558+reg557; reg370=reg219+reg370; reg675=reg672+reg675; reg647=reg305+reg647; reg301=reg642+reg301;
    reg713=reg712+reg713; reg355=reg552+reg355; reg665=reg670+reg665; reg350=reg350+reg709; reg662=reg662-reg664;
    reg353=reg353-reg576; reg180=reg38*reg705; reg615=reg615-reg616; reg603=reg610+reg603; reg344=reg344+reg703;
    reg346=reg597+reg346; reg618=reg624+reg618; reg428=reg428-reg460; reg224=reg791+reg478; reg294=reg793+reg477;
    reg299=reg795+reg474; reg473=reg38*reg473; reg95=reg95-reg789; reg305=reg786+reg472; reg471=reg38*reg471;
    reg329=reg794+reg470; reg334=reg775+reg466; reg338=reg807+reg465; reg356=reg804+reg495; reg400=reg768+reg400;
    reg491=reg38*reg491; reg372=reg783+reg490; reg374=reg784+reg489; reg397=reg700+reg483; reg398=reg701+reg482;
    reg417=reg418+reg417; reg509=reg511+reg509; reg118=reg782+reg118; reg166=reg362+reg166; reg502=reg503+reg502;
    reg500=reg500-reg501; reg126=reg126-reg691; reg395=reg395-reg364; reg770=reg396+reg770; reg497=reg499+reg497;
    reg524=reg526+reg524; reg768=reg373+reg768; reg371=reg367+reg371; reg767=reg151+reg767; reg151=reg796+reg516;
    reg362=reg797+reg512; reg762=reg760+reg762; reg367=reg798+reg479; reg461=reg625+reg461; reg179=reg735+reg179;
    reg435=reg451-reg435; reg220=reg217-reg220; reg688=reg689+reg688; reg409=reg414+reg409; reg404=reg406-reg404;
    reg225=reg181-reg225; reg750=reg750-reg748; reg181=reg38*reg341; reg525=reg525+reg388; reg186=reg678+reg186;
    reg217=reg38*reg129; reg734=reg734-reg728; reg494=reg464+reg494; reg487=reg488-reg487; reg213=reg175-reg213;
    reg729=reg727+reg729; reg581=reg581-reg580; reg121=reg197+reg121; reg175=reg38*reg200; reg570=reg579+reg570;
    reg764=reg154+reg764; reg599=reg565-reg599; reg234=reg588+reg234; reg538=reg538-reg544; reg167=reg774+reg167;
    reg154=reg38*reg231; reg806=reg222+reg806; reg229=reg598+reg229; reg197=reg38*reg656; reg162=reg114+reg162;
    reg609=reg617-reg609; reg218=reg104+reg218; reg702=reg708+reg702; reg626=reg626-reg632; reg578=reg578-reg576;
    reg575=reg575-reg572; reg747=reg168+reg747; reg114=reg38*reg566; reg743=reg742+reg743; reg597=reg132+reg597;
    reg595=reg596+reg595; reg590=reg591+reg590; reg735=reg366+reg735; reg584=reg363+reg584; reg542=reg542-reg543;
    reg687=reg686+reg687; reg539=reg541+reg539; reg684=reg684-reg683; reg530=reg136+reg530; reg561=reg562+reg561;
    reg559=reg560+reg559; reg678=reg112+reg678; reg480=reg481+reg480; reg112=reg38*reg401; reg145=reg709+reg145;
    reg119=reg209+reg119; reg520=reg517+reg520; reg721=reg720+reg721; reg425=reg506-reg425; reg281=reg138-reg281;
    reg778=reg778-reg805; reg219=reg140+reg219; reg718=reg717+reg718; reg657=reg604+reg657; reg669=reg661+reg669;
    reg715=reg714+reg715; reg552=reg147+reg552; reg586=reg533+reg586; reg357=reg703+reg357; reg571=reg587+reg571;
    reg109=reg751+reg109; reg515=reg514-reg515; reg629=reg629-reg628; reg111=reg155+reg111; reg132=reg38*reg237;
    reg336=reg336-reg809; reg333=reg333+reg434; reg627=reg296+reg627; reg711=reg707+reg711; reg519=reg518-reg519;
    reg429=reg258-reg429; reg194=reg194-reg802; reg239=reg239-reg781; reg261=reg771+reg261; reg136=reg38*reg259;
    reg699=reg699-reg698; reg623=reg622+reg623; reg521=reg238-reg521; reg143=reg685+reg143; reg722=reg722-reg733;
    reg634=reg633+reg634; reg462=reg170-reg462; reg343=reg106+reg343; reg600=reg91+reg600; reg91=reg38*reg283;
    reg405=reg316-reg405; reg347=reg298+reg347; reg475=reg475+reg476; reg332=reg287-reg332; reg286=reg315+reg286;
    reg106=reg38*reg285; reg138=reg38*reg317; reg631=reg630+reg631; reg724=reg725+reg724; reg692=reg690+reg692;
    reg792=reg208+reg792; reg288=reg288-reg513; reg262=reg268-reg262; reg337=reg337+reg441; reg248=reg493+reg248;
    reg140=reg38*reg251; reg433=reg433-reg257; reg420=reg420+reg421; reg754=reg752+reg754; reg681=reg682+reg681;
    reg505=reg504-reg505; reg147=reg38*reg419; reg732=reg726+reg732; reg739=reg740+reg739; reg253=reg253-reg467;
    reg449=reg449-reg448; reg443=reg442+reg443; reg289=reg704+reg289; reg254=reg173+reg254; reg507=reg252-reg507;
    reg508=reg508+reg510; reg445=reg444+reg445; reg436=reg790+reg436; reg523=reg522-reg523; reg427=reg426-reg427;
    reg621=reg620+reg621; reg155=reg38*reg706; reg242=reg242-reg484; reg117=reg716+reg117; reg741=reg741-reg746;
    reg498=reg241-reg498; reg438=reg437+reg438; reg513=reg256-reg513; reg619=reg293+reg619; reg737=reg736+reg737;
    reg247=reg381+reg247; reg168=reg38*reg246; reg423=reg422-reg423; reg673=reg674+reg673; reg170=reg38*reg453;
    reg749=reg753+reg749; reg788=reg221+reg788; reg608=reg607+reg608; reg173=reg38*reg271; reg679=reg679-reg680;
    reg606=reg605+reg606; reg761=reg131+reg761; reg765=reg157+reg765; reg492=reg492+reg493; reg272=reg272+reg413;
    reg115=reg787+reg115; reg131=reg38*reg275; reg157=reg38*reg128; reg300=reg434+reg300; reg277=reg277-reg309;
    reg208=reg496+reg463; reg221=reg38*reg320; reg757=reg113+reg757; reg614=reg614-reg613; reg776=reg122-reg776;
    reg113=reg38*reg328; reg122=reg38*reg264; reg326=reg276-reg326; reg697=reg133+reg697; reg612=reg611+reg612;
    reg781=reg267-reg781; reg484=reg266-reg484; reg330=reg766+reg330; reg304=reg441+reg304; reg416=reg324-reg416;
    reg779=reg779-reg780; reg269=reg303+reg269; reg485=reg485+reg486; reg327=reg327+reg454; reg133=reg38*reg325;
    reg467=reg278-reg467; reg772=reg772-reg773; reg104=reg282+reg104; reg771=reg302+reg771; reg636=reg635+reg636;
    reg222=reg38*reg280; reg602=reg601+reg602; reg238=reg38*reg331; reg260=reg260-reg763; reg469=reg468-reg469;
    reg241=reg407+reg408; reg318=reg265-reg318; reg117=reg38*reg117; reg480=reg38*reg480; reg405=reg38*reg405;
    reg252=ponderation*reg114; reg420=reg38*reg420; reg578=reg38*reg578; reg626=reg38*reg626; reg734=reg38*reg734;
    reg423=reg38*reg423; reg256=ponderation*reg217; reg697=reg38*reg697; reg722=reg38*reg722; reg747=reg38*reg747;
    reg575=reg38*reg575; reg461=reg38*reg461; reg584=reg38*reg584; reg162=reg38*reg162; reg609=reg38*reg609;
    reg508=reg38*reg508; reg735=reg38*reg735; reg590=reg38*reg590; reg739=reg38*reg739; reg326=reg38*reg326;
    reg702=reg38*reg702; reg487=reg38*reg487; reg729=reg38*reg729; reg318=reg38*reg318; reg595=reg38*reg595;
    reg218=reg38*reg218; reg213=reg38*reg213; reg494=reg38*reg494; reg258=ponderation*reg147; reg597=reg38*reg597;
    reg416=reg38*reg416; reg749=reg38*reg749; reg743=reg38*reg743; reg679=reg38*reg679; reg657=reg38*reg657;
    reg265=ponderation*reg136; reg404=reg38*reg404; reg119=reg38*reg119; reg219=reg38*reg219; reg757=reg38*reg757;
    reg750=reg38*reg750; reg718=reg38*reg718; reg721=reg38*reg721; reg225=reg38*reg225; reg429=reg38*reg429;
    reg778=reg38*reg778; reg281=reg38*reg281; reg266=ponderation*reg221; reg520=reg38*reg520; reg262=reg38*reg262;
    reg186=reg38*reg186; reg425=reg38*reg425; reg724=reg38*reg724; reg267=ponderation*reg181; reg179=reg38*reg179;
    reg435=reg38*reg435; reg109=reg38*reg109; reg571=reg38*reg571; reg268=ponderation*reg133; reg145=reg38*reg145;
    reg276=reg38*reg241; reg513=reg38*reg513; reg586=reg38*reg586; reg220=reg38*reg220; reg688=reg38*reg688;
    reg357=reg38*reg357; reg278=ponderation*reg112; reg552=reg38*reg552; reg409=reg38*reg409; reg272=reg38*reg272;
    reg427=reg38*reg427; reg715=reg38*reg715; reg669=reg38*reg669; reg282=ponderation*reg138; reg525=reg38*reg525;
    reg261=reg38*reg261; reg767=reg38*reg767; reg371=reg38*reg371; reg208=reg38*reg208; reg768=reg38*reg768;
    reg524=reg38*reg524; reg772=reg38*reg772; reg497=reg38*reg497; reg467=reg38*reg467; reg395=reg38*reg395;
    reg770=reg38*reg770; reg500=reg38*reg500; reg469=reg38*reg469; reg502=reg38*reg502; reg287=ponderation*reg222;
    reg95=reg38*reg95; reg166=reg38*reg166; reg118=reg38*reg118; reg509=reg38*reg509; reg600=reg38*reg600;
    reg126=reg38*reg126; reg417=reg38*reg417; reg293=ponderation*reg91; reg359=reg38*reg359; reg286=reg38*reg286;
    reg431=reg38*reg431; reg475=reg38*reg475; reg296=ponderation*reg177; reg298=ponderation*reg106; reg144=reg38*reg144;
    reg400=reg38*reg400; reg792=reg38*reg792; reg302=reg38*reg398; reg776=reg38*reg776; reg303=reg38*reg397;
    reg315=reg38*reg374; reg781=reg38*reg781; reg316=reg38*reg372; reg324=ponderation*reg122; reg491=ponderation*reg491;
    reg363=reg38*reg356; reg366=reg38*reg338; reg484=reg38*reg484; reg373=reg38*reg334; reg269=reg38*reg269;
    reg381=reg38*reg329; reg485=reg38*reg485; reg471=ponderation*reg471; reg396=reg38*reg305; reg788=reg38*reg788;
    reg473=ponderation*reg473; reg399=ponderation*reg173; reg403=reg38*reg299; reg406=reg38*reg294; reg765=reg38*reg765;
    reg411=reg38*reg224; reg414=reg38*reg367; reg492=reg38*reg492; reg762=reg38*reg762; reg415=reg38*reg362;
    reg418=reg38*reg151; reg422=ponderation*reg131; reg277=reg38*reg277; reg426=ponderation*reg180; reg353=reg38*reg353;
    reg498=reg38*reg498; reg247=reg38*reg247; reg662=reg38*reg662; reg350=reg38*reg350; reg665=reg38*reg665;
    reg434=ponderation*reg168; reg355=reg38*reg355; reg673=reg38*reg673; reg713=reg38*reg713; reg301=reg38*reg301;
    reg248=reg38*reg248; reg647=reg38*reg647; reg370=reg38*reg370; reg437=ponderation*reg140; reg675=reg38*reg675;
    reg557=reg38*reg557; reg681=reg38*reg681; reg559=reg38*reg559; reg505=reg38*reg505; reg678=reg38*reg678;
    reg561=reg38*reg561; reg530=reg38*reg530; reg253=reg38*reg253; reg684=reg38*reg684; reg254=reg38*reg254;
    reg539=reg38*reg539; reg507=reg38*reg507; reg542=reg38*reg542; reg687=reg38*reg687; reg410=reg38*reg410;
    reg288=reg38*reg288; reg412=reg38*reg412; reg111=reg38*reg111; reg150=reg38*reg150; reg515=reg38*reg515;
    reg452=reg38*reg452; reg456=reg38*reg456; reg440=ponderation*reg132; reg458=reg38*reg458; reg194=reg38*reg194;
    reg124=reg38*reg124; reg160=reg38*reg160; reg519=reg38*reg519; reg439=reg38*reg439; reg239=reg38*reg239;
    reg447=reg38*reg447; reg389=reg38*reg389; reg450=reg38*reg450; reg699=reg38*reg699; reg521=reg38*reg521;
    reg428=reg38*reg428; reg694=reg38*reg694; reg618=reg38*reg618; reg523=reg38*reg523; reg346=reg38*reg346;
    reg441=ponderation*reg155; reg603=reg38*reg603; reg242=reg38*reg242; reg344=reg38*reg344; reg615=reg38*reg615;
    reg438=reg38*reg438; reg676=reg38*reg676; reg216=reg38*reg216; reg636=reg38*reg636; reg641=reg38*reg641;
    reg634=reg38*reg634; reg594=reg38*reg594; reg732=reg38*reg732; reg758=reg38*reg758; reg190=reg38*reg190;
    reg592=reg38*reg592; reg337=reg38*reg337; reg801=reg38*reg801; reg644=reg38*reg644; reg214=reg38*reg214;
    reg351=reg38*reg351; reg443=reg38*reg443; reg347=reg38*reg347; reg646=reg38*reg646; reg442=ponderation*reg178;
    reg730=reg38*reg730; reg343=reg38*reg343; reg195=reg38*reg195; reg444=ponderation*reg161; reg663=reg38*reg663;
    reg203=reg38*reg203; reg336=reg38*reg336; reg451=ponderation*reg135; reg300=reg38*reg300; reg668=reg38*reg668;
    reg777=reg38*reg777; reg333=reg38*reg333; reg568=reg38*reg568; reg141=reg38*reg141; reg164=reg38*reg164;
    reg187=reg38*reg187; reg564=reg38*reg564; reg436=reg38*reg436; reg386=reg38*reg386; reg602=reg38*reg602;
    reg104=reg38*reg104; reg637=reg38*reg637; reg771=reg38*reg771; reg455=ponderation*reg212; reg741=reg38*reg741;
    reg639=reg38*reg639; reg211=reg38*reg211; reg710=reg38*reg710; reg529=reg38*reg529; reg457=ponderation*reg170;
    reg547=reg38*reg547; reg737=reg38*reg737; reg527=reg38*reg527; reg627=reg38*reg627; reg744=reg38*reg744;
    reg171=reg38*reg171; reg619=reg38*reg619; reg711=reg38*reg711; reg174=reg38*reg174; reg365=reg38*reg365;
    reg549=reg38*reg549; reg556=reg38*reg556; reg143=reg38*reg143; reg677=reg38*reg677; reg554=reg38*reg554;
    reg621=reg38*reg621; reg623=reg38*reg623; reg551=reg38*reg551; reg394=reg38*reg394; reg172=reg38*reg172;
    reg582=reg38*reg582; reg445=reg38*reg445; reg723=reg38*reg723; reg693=reg38*reg693; reg459=ponderation*reg310;
    reg289=reg38*reg289; reg631=reg38*reg631; reg719=reg38*reg719; reg464=ponderation*reg295; reg540=reg38*reg540;
    reg449=reg38*reg449; reg648=reg38*reg648; reg188=reg38*reg188; reg692=reg38*reg692; reg536=reg38*reg536;
    reg754=reg38*reg754; reg226=reg38*reg226; reg534=reg38*reg534; reg340=reg38*reg340; reg629=reg38*reg629;
    reg433=reg38*reg433; reg545=reg38*reg545; reg532=reg38*reg532; reg755=reg38*reg755; reg314=reg38*reg314;
    reg779=reg38*reg779; reg608=reg38*reg608; reg764=reg38*reg764; reg198=reg38*reg198; reg345=reg38*reg345;
    reg167=reg38*reg167; reg462=reg38*reg462; reg229=reg38*reg229; reg653=reg38*reg653; reg199=reg38*reg199;
    reg311=reg38*reg311; reg806=reg38*reg806; reg658=reg38*reg658; reg659=reg38*reg659; reg612=reg38*reg612;
    reg306=reg38*reg306; reg570=reg38*reg570; reg260=reg38*reg260; reg121=reg38*reg121; reg468=ponderation*reg154;
    reg327=reg38*reg327; reg581=reg38*reg581; reg538=reg38*reg538; reg655=reg38*reg655; reg614=reg38*reg614;
    reg115=reg38*reg115; reg481=ponderation*reg197; reg573=reg38*reg573; reg761=reg38*reg761; reg488=ponderation*reg175;
    reg308=reg38*reg308; reg332=reg38*reg332; reg330=reg38*reg330; reg499=ponderation*reg238; reg108=reg38*reg108;
    reg599=reg38*reg599; reg313=reg38*reg313; reg577=reg38*reg577; reg354=reg38*reg354; reg234=reg38*reg234;
    reg606=reg38*reg606; reg503=ponderation*reg113; reg304=reg38*reg304; reg504=ponderation*reg157; T tmp_11_8=ponderation*reg340;
    T tmp_17_2=ponderation*reg150; T tmp_6_16=ponderation*reg515; T tmp_1_5=ponderation*reg118; T tmp_3_17=ponderation*reg692; T tmp_17_1=ponderation*reg452;
    T tmp_6_2=-reg324; T tmp_5_9=ponderation*reg111; reg111=ponderation*reg363; sollicitation[indices[4]+0]+=reg111; T tmp_6_15=ponderation*reg288;
    T tmp_9_13=ponderation*reg629; T tmp_17_3=ponderation*reg412; T tmp_10_8=ponderation*reg653; T tmp_5_17=ponderation*reg269; T tmp_10_4=ponderation*reg304;
    T tmp_11_7=ponderation*reg648; reg118=ponderation*reg366; sollicitation[indices[3]+2]+=reg118; T tmp_17_4=ponderation*reg410; T tmp_3_3=ponderation*reg188;
    T tmp_11_6=-reg464; T tmp_6_14=-reg298; T tmp_17_5=ponderation*reg144; T tmp_6_3=ponderation*reg484; T tmp_7_0=ponderation*reg519;
    T tmp_3_1=ponderation*reg677; T tmp_16_15=ponderation*reg439; T tmp_4_0=ponderation*reg711; T tmp_6_0=ponderation*reg781; T tmp_5_7=ponderation*reg699;
    T tmp_6_1=ponderation*reg776; T tmp_11_12=ponderation*reg549; T tmp_16_14=ponderation*reg447; T tmp_3_12=ponderation*reg308; reg144=ponderation*reg303;
    sollicitation[indices[5]+1]+=reg144; T tmp_7_1=ponderation*reg239; T tmp_9_11=ponderation*reg623; T tmp_0_6=ponderation*reg389; T tmp_16_13=ponderation*reg450;
    T tmp_10_6=ponderation*reg614; T tmp_11_13=ponderation*reg551; T tmp_3_0=ponderation*reg394; reg150=ponderation*reg302; sollicitation[indices[5]+2]+=reg150;
    T tmp_11_9=ponderation*reg545; T tmp_0_4=ponderation*reg770; sollicitation[indices[4]+1]+=-reg491; T tmp_3_2=ponderation*reg710; T tmp_17_0=ponderation*reg456;
    T tmp_1_6=ponderation*reg95; T tmp_9_12=ponderation*reg627; T tmp_6_17=-reg440; T tmp_5_8=ponderation*reg194; T tmp_3_11=ponderation*reg108;
    T tmp_16_17=ponderation*reg458; T tmp_10_7=ponderation*reg306; T tmp_11_10=ponderation*reg547; reg95=ponderation*reg316; sollicitation[indices[4]+2]+=reg95;
    T tmp_0_5=ponderation*reg124; T tmp_16_16=ponderation*reg160; T tmp_10_5=ponderation*reg612; T tmp_11_11=ponderation*reg171; reg108=ponderation*reg315;
    sollicitation[indices[5]+0]+=reg108; T tmp_17_16=ponderation*reg524; T tmp_10_16=ponderation*reg187; T tmp_10_12=ponderation*reg659; T tmp_5_13=ponderation*reg772;
    T tmp_0_3=ponderation*reg768; reg124=ponderation*reg406; sollicitation[indices[1]+1]+=reg124; T tmp_17_15=ponderation*reg497; T tmp_0_0=ponderation*reg344;
    T tmp_10_17=ponderation*reg637; T tmp_6_5=-reg399; T tmp_3_6=ponderation*reg216; T tmp_17_14=ponderation*reg395; T tmp_9_17=ponderation*reg636;
    T tmp_6_9=ponderation*reg467; reg160=ponderation*reg403; sollicitation[indices[1]+2]+=reg160; T tmp_10_2=ponderation*reg606; T tmp_11_0=ponderation*reg639;
    T tmp_17_13=ponderation*reg500; T tmp_3_9=ponderation*reg314; T tmp_10_11=ponderation*reg658; T tmp_10_1=ponderation*reg300; T tmp_0_1=ponderation*reg762;
    reg171=ponderation*reg415; sollicitation[indices[0]+1]+=reg171; T tmp_6_6=ponderation*reg492; T tmp_3_8=-reg451; T tmp_5_14=ponderation*reg277;
    T tmp_10_14=ponderation*reg663; reg187=ponderation*reg418; sollicitation[indices[0]+0]+=reg187; reg188=ponderation*reg414; sollicitation[indices[0]+2]+=reg188;
    T tmp_6_7=-reg422; T tmp_10_13=ponderation*reg313; T tmp_10_15=ponderation*reg668; T tmp_0_2=ponderation*reg767; T tmp_3_14=ponderation*reg761;
    T tmp_3_7=ponderation*reg141; T tmp_17_17=ponderation*reg371; reg141=ponderation*reg411; sollicitation[indices[1]+0]+=reg141; T tmp_10_0=ponderation*reg602;
    T tmp_6_8=ponderation*reg208; T tmp_5_15=ponderation*reg765; T tmp_11_3=ponderation*reg644; T tmp_9_15=ponderation*reg347; sollicitation[indices[2]+2]+=-reg471;
    T tmp_10_10=ponderation*reg311; T tmp_17_8=ponderation*reg359; T tmp_11_4=ponderation*reg646; T tmp_6_12=-reg293; T tmp_5_11=ponderation*reg286;
    T tmp_3_4=ponderation*reg693; T tmp_17_7=ponderation*reg431; reg194=ponderation*reg381; sollicitation[indices[3]+0]+=reg194; T tmp_11_5=ponderation*reg195;
    T tmp_6_13=ponderation*reg475; T tmp_17_6=-reg296; T tmp_3_10=ponderation*reg354; T tmp_10_9=ponderation*reg655; T tmp_3_16=ponderation*reg343;
    reg195=ponderation*reg373; sollicitation[indices[3]+1]+=reg195; T tmp_5_10=ponderation*reg792; T tmp_9_14=ponderation*reg631; T tmp_6_10=ponderation*reg469;
    T tmp_17_12=ponderation*reg502; T tmp_11_1=ponderation*reg641; T tmp_3_15=ponderation*reg771; T tmp_9_16=ponderation*reg634; sollicitation[indices[2]+0]+=-reg473;
    T tmp_6_4=ponderation*reg730; T tmp_17_11=ponderation*reg166; T tmp_6_11=-reg287; T tmp_10_3=ponderation*reg608; T tmp_12_14=ponderation*reg600;
    reg166=ponderation*reg396; sollicitation[indices[2]+1]+=reg166; T tmp_11_2=ponderation*reg190; T tmp_17_10=ponderation*reg509; T tmp_3_5=ponderation*reg801;
    T tmp_1_4=ponderation*reg400; T tmp_1_7=ponderation*reg126; T tmp_5_16=ponderation*reg788; T tmp_13_6=ponderation*reg485; T tmp_3_13=ponderation*reg779;
    T tmp_17_9=ponderation*reg417; T tmp_4_15=ponderation*reg724; T tmp_14_16=ponderation*reg425; T tmp_2_6=-reg444; T tmp_13_1=ponderation*reg203;
    T tmp_14_15=ponderation*reg520; T tmp_8_2=ponderation*reg262; T tmp_4_7=ponderation*reg336; T tmp_1_9=ponderation*reg721; T tmp_14_14=ponderation*reg119;
    T tmp_8_17=ponderation*reg332; T tmp_13_2=ponderation*reg573; T tmp_14_13=-reg278; T tmp_2_5=ponderation*reg198; T tmp_8_3=-reg282;
    T tmp_1_10=ponderation*reg145; T tmp_4_14=ponderation*reg722; T tmp_14_12=ponderation*reg480; T tmp_13_3=ponderation*reg577; T tmp_8_16=ponderation*reg462;
    T tmp_14_11=ponderation*reg213; T tmp_8_4=ponderation*reg405; T tmp_13_4=ponderation*reg199; T tmp_1_11=ponderation*reg729; T tmp_14_10=ponderation*reg487;
    T tmp_2_4=ponderation*reg121; T tmp_7_16=ponderation*reg513; T tmp_15_4=ponderation*reg586; T tmp_9_1=ponderation*reg436; T tmp_12_13=-reg455;
    T tmp_1_1=ponderation*reg357; T tmp_15_3=ponderation*reg552; T tmp_2_8=ponderation*reg386; T tmp_12_15=ponderation*reg104; T tmp_7_17=ponderation*reg427;
    T tmp_4_16=ponderation*reg261; T tmp_15_2=ponderation*reg669; T tmp_1_2=ponderation*reg715; T tmp_12_16=ponderation*reg564; T tmp_15_1=ponderation*reg657;
    T tmp_2_7=ponderation*reg164; T tmp_8_0=-reg265; T tmp_4_6=ponderation*reg741; T tmp_15_0=ponderation*reg219; T tmp_12_17=ponderation*reg568;
    T tmp_9_0=ponderation*reg333; T tmp_1_3=ponderation*reg718; T tmp_14_17=ponderation*reg281; T tmp_13_0=ponderation*reg777; T tmp_8_1=ponderation*reg429;
    T tmp_1_8=ponderation*reg778; T tmp_8_8=ponderation*reg272; T tmp_5_12=ponderation*reg167; T tmp_4_9=ponderation*reg115; T tmp_14_2=ponderation*reg220;
    T tmp_13_11=ponderation*reg538; T tmp_1_15=ponderation*reg688; T tmp_14_1=ponderation*reg435; T tmp_8_13=ponderation*reg327; T tmp_8_9=-reg268;
    T tmp_4_11=ponderation*reg697; T tmp_14_0=ponderation*reg461; T tmp_13_12=-reg468; T tmp_1_16=ponderation*reg179; T tmp_2_1=ponderation*reg806;
    T tmp_13_17=ponderation*reg626; T tmp_8_12=-reg503; T tmp_8_10=ponderation*reg416; T tmp_13_16=ponderation*reg218; T tmp_13_13=ponderation*reg229;
    T tmp_2_0=ponderation*reg162; T tmp_1_17=ponderation*reg702; T tmp_13_15=ponderation*reg609; T tmp_8_11=ponderation*reg326; T tmp_4_10=ponderation*reg330;
    T tmp_13_14=-reg481; T tmp_13_5=ponderation*reg581; T tmp_14_9=ponderation*reg494; T tmp_8_5=ponderation*reg318; T tmp_4_13=ponderation*reg757;
    T tmp_8_15=-reg499; T tmp_14_8=-reg256; T tmp_1_12=ponderation*reg734; T tmp_13_7=-reg488; T tmp_8_6=ponderation*reg276;
    T tmp_14_7=ponderation*reg525; T tmp_2_3=ponderation*reg764; T tmp_4_8=ponderation*reg260; T tmp_1_13=ponderation*reg186; T tmp_14_6=-reg267;
    T tmp_13_8=ponderation*reg570; T tmp_8_14=-reg504; T tmp_8_7=-reg266; T tmp_14_5=ponderation*reg225; T tmp_4_12=ponderation*reg679;
    T tmp_13_9=ponderation*reg599; T tmp_1_14=ponderation*reg750; T tmp_14_4=ponderation*reg404; T tmp_2_2=ponderation*reg345; T tmp_13_10=ponderation*reg234;
    T tmp_14_3=ponderation*reg409; T tmp_9_8=-reg457; T tmp_16_5=ponderation*reg665; T tmp_0_9=ponderation*reg350; T tmp_7_6=-reg434;
    T tmp_16_4=ponderation*reg355; T tmp_12_1=ponderation*reg529; T tmp_5_4=ponderation*reg673; T tmp_2_15=ponderation*reg755; T tmp_16_3=ponderation*reg301;
    T tmp_12_2=ponderation*reg532; T tmp_4_2=ponderation*reg737; T tmp_0_10=ponderation*reg713; T tmp_9_7=ponderation*reg433; T tmp_16_2=ponderation*reg647;
    T tmp_7_7=ponderation*reg248; T tmp_5_3=ponderation*reg681; T tmp_2_14=ponderation*reg226; T tmp_16_1=ponderation*reg370; T tmp_12_3=ponderation*reg534;
    T tmp_0_11=ponderation*reg675; T tmp_7_8=-reg437; T tmp_16_0=ponderation*reg557; T tmp_12_4=ponderation*reg536; T tmp_9_6=ponderation*reg449;
    T tmp_15_17=ponderation*reg559; T tmp_16_12=ponderation*reg428; T tmp_11_14=ponderation*reg172; T tmp_7_2=ponderation*reg521; T tmp_0_7=ponderation*reg694;
    T tmp_16_11=ponderation*reg618; T tmp_9_10=ponderation*reg621; T tmp_7_3=ponderation*reg523; T tmp_16_10=ponderation*reg346; T tmp_11_15=ponderation*reg554;
    T tmp_5_6=-reg441; T tmp_4_1=ponderation*reg143; T tmp_16_9=ponderation*reg603; T tmp_11_16=ponderation*reg556; T tmp_2_17=ponderation*reg365;
    T tmp_7_4=ponderation*reg242; T tmp_16_8=ponderation*reg615; T tmp_9_9=ponderation*reg619; T tmp_5_5=ponderation*reg247; T tmp_11_17=ponderation*reg174;
    T tmp_16_7=ponderation*reg353; T tmp_2_16=ponderation*reg744; T tmp_0_8=-reg426; T tmp_7_5=ponderation*reg498; T tmp_16_6=ponderation*reg662;
    T tmp_12_0=ponderation*reg527; T tmp_0_15=ponderation*reg735; T tmp_2_11=ponderation*reg351; T tmp_15_10=ponderation*reg595; T tmp_12_9=ponderation*reg214;
    T tmp_7_13=-reg258; T tmp_15_9=ponderation*reg597; T tmp_9_3=ponderation*reg337; T tmp_5_0=ponderation*reg749; T tmp_2_10=ponderation*reg758;
    T tmp_0_16=ponderation*reg743; T tmp_15_8=-reg252; T tmp_12_10=ponderation*reg592; T tmp_7_14=ponderation*reg420; T tmp_15_7=ponderation*reg575;
    T tmp_12_11=ponderation*reg594; T tmp_4_5=ponderation*reg732; T tmp_0_17=ponderation*reg747; T tmp_9_2=ponderation*reg438; T tmp_15_6=ponderation*reg578;
    T tmp_2_9=ponderation*reg676; T tmp_7_15=ponderation*reg423; T tmp_4_17=ponderation*reg117; T tmp_1_0=ponderation*reg109; T tmp_15_5=ponderation*reg571;
    T tmp_12_12=ponderation*reg211; T tmp_2_13=ponderation*reg719; T tmp_0_12=ponderation*reg678; T tmp_15_16=ponderation*reg561; T tmp_12_5=ponderation*reg540;
    T tmp_7_9=ponderation*reg505; T tmp_5_2=ponderation*reg254; T tmp_15_15=ponderation*reg530; T tmp_4_3=ponderation*reg754; T tmp_2_12=ponderation*reg723;
    T tmp_7_10=ponderation*reg253; T tmp_0_13=ponderation*reg684; T tmp_15_14=ponderation*reg539; T tmp_12_6=-reg459; T tmp_9_5=ponderation*reg445;
    T tmp_15_13=ponderation*reg542; T tmp_7_11=ponderation*reg507; T tmp_12_7=ponderation*reg582; T tmp_0_14=ponderation*reg687; T tmp_15_12=ponderation*reg584;
    T tmp_5_1=ponderation*reg739; T tmp_4_4=ponderation*reg289; T tmp_9_4=ponderation*reg443; T tmp_7_12=ponderation*reg508; T tmp_15_11=ponderation*reg590;
    T tmp_12_8=-reg442;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[0]+2) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_3;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_4;
    matrix(indices[0]+0,indices[1]+2) += tmp_0_5;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_7;
    matrix(indices[0]+0,indices[2]+2) += tmp_0_8;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_9;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_10;
    matrix(indices[0]+0,indices[3]+2) += tmp_0_11;
    matrix(indices[0]+0,indices[4]+0) += tmp_0_12;
    matrix(indices[0]+0,indices[4]+1) += tmp_0_13;
    matrix(indices[0]+0,indices[4]+2) += tmp_0_14;
    matrix(indices[0]+0,indices[5]+0) += tmp_0_15;
    matrix(indices[0]+0,indices[5]+1) += tmp_0_16;
    matrix(indices[0]+0,indices[5]+2) += tmp_0_17;
    matrix(indices[0]+1,indices[0]+0) += tmp_1_0;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[0]+2) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_3;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_4;
    matrix(indices[0]+1,indices[1]+2) += tmp_1_5;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_7;
    matrix(indices[0]+1,indices[2]+2) += tmp_1_8;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_9;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_10;
    matrix(indices[0]+1,indices[3]+2) += tmp_1_11;
    matrix(indices[0]+1,indices[4]+0) += tmp_1_12;
    matrix(indices[0]+1,indices[4]+1) += tmp_1_13;
    matrix(indices[0]+1,indices[4]+2) += tmp_1_14;
    matrix(indices[0]+1,indices[5]+0) += tmp_1_15;
    matrix(indices[0]+1,indices[5]+1) += tmp_1_16;
    matrix(indices[0]+1,indices[5]+2) += tmp_1_17;
    matrix(indices[0]+2,indices[0]+0) += tmp_2_0;
    matrix(indices[0]+2,indices[0]+1) += tmp_2_1;
    matrix(indices[0]+2,indices[0]+2) += tmp_2_2;
    matrix(indices[0]+2,indices[1]+0) += tmp_2_3;
    matrix(indices[0]+2,indices[1]+1) += tmp_2_4;
    matrix(indices[0]+2,indices[1]+2) += tmp_2_5;
    matrix(indices[0]+2,indices[2]+0) += tmp_2_6;
    matrix(indices[0]+2,indices[2]+1) += tmp_2_7;
    matrix(indices[0]+2,indices[2]+2) += tmp_2_8;
    matrix(indices[0]+2,indices[3]+0) += tmp_2_9;
    matrix(indices[0]+2,indices[3]+1) += tmp_2_10;
    matrix(indices[0]+2,indices[3]+2) += tmp_2_11;
    matrix(indices[0]+2,indices[4]+0) += tmp_2_12;
    matrix(indices[0]+2,indices[4]+1) += tmp_2_13;
    matrix(indices[0]+2,indices[4]+2) += tmp_2_14;
    matrix(indices[0]+2,indices[5]+0) += tmp_2_15;
    matrix(indices[0]+2,indices[5]+1) += tmp_2_16;
    matrix(indices[0]+2,indices[5]+2) += tmp_2_17;
    matrix(indices[1]+0,indices[0]+0) += tmp_3_0;
    matrix(indices[1]+0,indices[0]+1) += tmp_3_1;
    matrix(indices[1]+0,indices[0]+2) += tmp_3_2;
    matrix(indices[1]+0,indices[1]+0) += tmp_3_3;
    matrix(indices[1]+0,indices[1]+1) += tmp_3_4;
    matrix(indices[1]+0,indices[1]+2) += tmp_3_5;
    matrix(indices[1]+0,indices[2]+0) += tmp_3_6;
    matrix(indices[1]+0,indices[2]+1) += tmp_3_7;
    matrix(indices[1]+0,indices[2]+2) += tmp_3_8;
    matrix(indices[1]+0,indices[3]+0) += tmp_3_9;
    matrix(indices[1]+0,indices[3]+1) += tmp_3_10;
    matrix(indices[1]+0,indices[3]+2) += tmp_3_11;
    matrix(indices[1]+0,indices[4]+0) += tmp_3_12;
    matrix(indices[1]+0,indices[4]+1) += tmp_3_13;
    matrix(indices[1]+0,indices[4]+2) += tmp_3_14;
    matrix(indices[1]+0,indices[5]+0) += tmp_3_15;
    matrix(indices[1]+0,indices[5]+1) += tmp_3_16;
    matrix(indices[1]+0,indices[5]+2) += tmp_3_17;
    matrix(indices[1]+1,indices[0]+0) += tmp_4_0;
    matrix(indices[1]+1,indices[0]+1) += tmp_4_1;
    matrix(indices[1]+1,indices[0]+2) += tmp_4_2;
    matrix(indices[1]+1,indices[1]+0) += tmp_4_3;
    matrix(indices[1]+1,indices[1]+1) += tmp_4_4;
    matrix(indices[1]+1,indices[1]+2) += tmp_4_5;
    matrix(indices[1]+1,indices[2]+0) += tmp_4_6;
    matrix(indices[1]+1,indices[2]+1) += tmp_4_7;
    matrix(indices[1]+1,indices[2]+2) += tmp_4_8;
    matrix(indices[1]+1,indices[3]+0) += tmp_4_9;
    matrix(indices[1]+1,indices[3]+1) += tmp_4_10;
    matrix(indices[1]+1,indices[3]+2) += tmp_4_11;
    matrix(indices[1]+1,indices[4]+0) += tmp_4_12;
    matrix(indices[1]+1,indices[4]+1) += tmp_4_13;
    matrix(indices[1]+1,indices[4]+2) += tmp_4_14;
    matrix(indices[1]+1,indices[5]+0) += tmp_4_15;
    matrix(indices[1]+1,indices[5]+1) += tmp_4_16;
    matrix(indices[1]+1,indices[5]+2) += tmp_4_17;
    matrix(indices[1]+2,indices[0]+0) += tmp_5_0;
    matrix(indices[1]+2,indices[0]+1) += tmp_5_1;
    matrix(indices[1]+2,indices[0]+2) += tmp_5_2;
    matrix(indices[1]+2,indices[1]+0) += tmp_5_3;
    matrix(indices[1]+2,indices[1]+1) += tmp_5_4;
    matrix(indices[1]+2,indices[1]+2) += tmp_5_5;
    matrix(indices[1]+2,indices[2]+0) += tmp_5_6;
    matrix(indices[1]+2,indices[2]+1) += tmp_5_7;
    matrix(indices[1]+2,indices[2]+2) += tmp_5_8;
    matrix(indices[1]+2,indices[3]+0) += tmp_5_9;
    matrix(indices[1]+2,indices[3]+1) += tmp_5_10;
    matrix(indices[1]+2,indices[3]+2) += tmp_5_11;
    matrix(indices[1]+2,indices[4]+0) += tmp_5_12;
    matrix(indices[1]+2,indices[4]+1) += tmp_5_13;
    matrix(indices[1]+2,indices[4]+2) += tmp_5_14;
    matrix(indices[1]+2,indices[5]+0) += tmp_5_15;
    matrix(indices[1]+2,indices[5]+1) += tmp_5_16;
    matrix(indices[1]+2,indices[5]+2) += tmp_5_17;
    matrix(indices[2]+0,indices[0]+0) += tmp_6_0;
    matrix(indices[2]+0,indices[0]+1) += tmp_6_1;
    matrix(indices[2]+0,indices[0]+2) += tmp_6_2;
    matrix(indices[2]+0,indices[1]+0) += tmp_6_3;
    matrix(indices[2]+0,indices[1]+1) += tmp_6_4;
    matrix(indices[2]+0,indices[1]+2) += tmp_6_5;
    matrix(indices[2]+0,indices[2]+0) += tmp_6_6;
    matrix(indices[2]+0,indices[2]+1) += tmp_6_7;
    matrix(indices[2]+0,indices[2]+2) += tmp_6_8;
    matrix(indices[2]+0,indices[3]+0) += tmp_6_9;
    matrix(indices[2]+0,indices[3]+1) += tmp_6_10;
    matrix(indices[2]+0,indices[3]+2) += tmp_6_11;
    matrix(indices[2]+0,indices[4]+0) += tmp_6_12;
    matrix(indices[2]+0,indices[4]+1) += tmp_6_13;
    matrix(indices[2]+0,indices[4]+2) += tmp_6_14;
    matrix(indices[2]+0,indices[5]+0) += tmp_6_15;
    matrix(indices[2]+0,indices[5]+1) += tmp_6_16;
    matrix(indices[2]+0,indices[5]+2) += tmp_6_17;
    matrix(indices[2]+1,indices[0]+0) += tmp_7_0;
    matrix(indices[2]+1,indices[0]+1) += tmp_7_1;
    matrix(indices[2]+1,indices[0]+2) += tmp_7_2;
    matrix(indices[2]+1,indices[1]+0) += tmp_7_3;
    matrix(indices[2]+1,indices[1]+1) += tmp_7_4;
    matrix(indices[2]+1,indices[1]+2) += tmp_7_5;
    matrix(indices[2]+1,indices[2]+0) += tmp_7_6;
    matrix(indices[2]+1,indices[2]+1) += tmp_7_7;
    matrix(indices[2]+1,indices[2]+2) += tmp_7_8;
    matrix(indices[2]+1,indices[3]+0) += tmp_7_9;
    matrix(indices[2]+1,indices[3]+1) += tmp_7_10;
    matrix(indices[2]+1,indices[3]+2) += tmp_7_11;
    matrix(indices[2]+1,indices[4]+0) += tmp_7_12;
    matrix(indices[2]+1,indices[4]+1) += tmp_7_13;
    matrix(indices[2]+1,indices[4]+2) += tmp_7_14;
    matrix(indices[2]+1,indices[5]+0) += tmp_7_15;
    matrix(indices[2]+1,indices[5]+1) += tmp_7_16;
    matrix(indices[2]+1,indices[5]+2) += tmp_7_17;
    matrix(indices[2]+2,indices[0]+0) += tmp_8_0;
    matrix(indices[2]+2,indices[0]+1) += tmp_8_1;
    matrix(indices[2]+2,indices[0]+2) += tmp_8_2;
    matrix(indices[2]+2,indices[1]+0) += tmp_8_3;
    matrix(indices[2]+2,indices[1]+1) += tmp_8_4;
    matrix(indices[2]+2,indices[1]+2) += tmp_8_5;
    matrix(indices[2]+2,indices[2]+0) += tmp_8_6;
    matrix(indices[2]+2,indices[2]+1) += tmp_8_7;
    matrix(indices[2]+2,indices[2]+2) += tmp_8_8;
    matrix(indices[2]+2,indices[3]+0) += tmp_8_9;
    matrix(indices[2]+2,indices[3]+1) += tmp_8_10;
    matrix(indices[2]+2,indices[3]+2) += tmp_8_11;
    matrix(indices[2]+2,indices[4]+0) += tmp_8_12;
    matrix(indices[2]+2,indices[4]+1) += tmp_8_13;
    matrix(indices[2]+2,indices[4]+2) += tmp_8_14;
    matrix(indices[2]+2,indices[5]+0) += tmp_8_15;
    matrix(indices[2]+2,indices[5]+1) += tmp_8_16;
    matrix(indices[2]+2,indices[5]+2) += tmp_8_17;
    matrix(indices[3]+0,indices[0]+0) += tmp_9_0;
    matrix(indices[3]+0,indices[0]+1) += tmp_9_1;
    matrix(indices[3]+0,indices[0]+2) += tmp_9_2;
    matrix(indices[3]+0,indices[1]+0) += tmp_9_3;
    matrix(indices[3]+0,indices[1]+1) += tmp_9_4;
    matrix(indices[3]+0,indices[1]+2) += tmp_9_5;
    matrix(indices[3]+0,indices[2]+0) += tmp_9_6;
    matrix(indices[3]+0,indices[2]+1) += tmp_9_7;
    matrix(indices[3]+0,indices[2]+2) += tmp_9_8;
    matrix(indices[3]+0,indices[3]+0) += tmp_9_9;
    matrix(indices[3]+0,indices[3]+1) += tmp_9_10;
    matrix(indices[3]+0,indices[3]+2) += tmp_9_11;
    matrix(indices[3]+0,indices[4]+0) += tmp_9_12;
    matrix(indices[3]+0,indices[4]+1) += tmp_9_13;
    matrix(indices[3]+0,indices[4]+2) += tmp_9_14;
    matrix(indices[3]+0,indices[5]+0) += tmp_9_15;
    matrix(indices[3]+0,indices[5]+1) += tmp_9_16;
    matrix(indices[3]+0,indices[5]+2) += tmp_9_17;
    matrix(indices[3]+1,indices[0]+0) += tmp_10_0;
    matrix(indices[3]+1,indices[0]+1) += tmp_10_1;
    matrix(indices[3]+1,indices[0]+2) += tmp_10_2;
    matrix(indices[3]+1,indices[1]+0) += tmp_10_3;
    matrix(indices[3]+1,indices[1]+1) += tmp_10_4;
    matrix(indices[3]+1,indices[1]+2) += tmp_10_5;
    matrix(indices[3]+1,indices[2]+0) += tmp_10_6;
    matrix(indices[3]+1,indices[2]+1) += tmp_10_7;
    matrix(indices[3]+1,indices[2]+2) += tmp_10_8;
    matrix(indices[3]+1,indices[3]+0) += tmp_10_9;
    matrix(indices[3]+1,indices[3]+1) += tmp_10_10;
    matrix(indices[3]+1,indices[3]+2) += tmp_10_11;
    matrix(indices[3]+1,indices[4]+0) += tmp_10_12;
    matrix(indices[3]+1,indices[4]+1) += tmp_10_13;
    matrix(indices[3]+1,indices[4]+2) += tmp_10_14;
    matrix(indices[3]+1,indices[5]+0) += tmp_10_15;
    matrix(indices[3]+1,indices[5]+1) += tmp_10_16;
    matrix(indices[3]+1,indices[5]+2) += tmp_10_17;
    matrix(indices[3]+2,indices[0]+0) += tmp_11_0;
    matrix(indices[3]+2,indices[0]+1) += tmp_11_1;
    matrix(indices[3]+2,indices[0]+2) += tmp_11_2;
    matrix(indices[3]+2,indices[1]+0) += tmp_11_3;
    matrix(indices[3]+2,indices[1]+1) += tmp_11_4;
    matrix(indices[3]+2,indices[1]+2) += tmp_11_5;
    matrix(indices[3]+2,indices[2]+0) += tmp_11_6;
    matrix(indices[3]+2,indices[2]+1) += tmp_11_7;
    matrix(indices[3]+2,indices[2]+2) += tmp_11_8;
    matrix(indices[3]+2,indices[3]+0) += tmp_11_9;
    matrix(indices[3]+2,indices[3]+1) += tmp_11_10;
    matrix(indices[3]+2,indices[3]+2) += tmp_11_11;
    matrix(indices[3]+2,indices[4]+0) += tmp_11_12;
    matrix(indices[3]+2,indices[4]+1) += tmp_11_13;
    matrix(indices[3]+2,indices[4]+2) += tmp_11_14;
    matrix(indices[3]+2,indices[5]+0) += tmp_11_15;
    matrix(indices[3]+2,indices[5]+1) += tmp_11_16;
    matrix(indices[3]+2,indices[5]+2) += tmp_11_17;
    matrix(indices[4]+0,indices[0]+0) += tmp_12_0;
    matrix(indices[4]+0,indices[0]+1) += tmp_12_1;
    matrix(indices[4]+0,indices[0]+2) += tmp_12_2;
    matrix(indices[4]+0,indices[1]+0) += tmp_12_3;
    matrix(indices[4]+0,indices[1]+1) += tmp_12_4;
    matrix(indices[4]+0,indices[1]+2) += tmp_12_5;
    matrix(indices[4]+0,indices[2]+0) += tmp_12_6;
    matrix(indices[4]+0,indices[2]+1) += tmp_12_7;
    matrix(indices[4]+0,indices[2]+2) += tmp_12_8;
    matrix(indices[4]+0,indices[3]+0) += tmp_12_9;
    matrix(indices[4]+0,indices[3]+1) += tmp_12_10;
    matrix(indices[4]+0,indices[3]+2) += tmp_12_11;
    matrix(indices[4]+0,indices[4]+0) += tmp_12_12;
    matrix(indices[4]+0,indices[4]+1) += tmp_12_13;
    matrix(indices[4]+0,indices[4]+2) += tmp_12_14;
    matrix(indices[4]+0,indices[5]+0) += tmp_12_15;
    matrix(indices[4]+0,indices[5]+1) += tmp_12_16;
    matrix(indices[4]+0,indices[5]+2) += tmp_12_17;
    matrix(indices[4]+1,indices[0]+0) += tmp_13_0;
    matrix(indices[4]+1,indices[0]+1) += tmp_13_1;
    matrix(indices[4]+1,indices[0]+2) += tmp_13_2;
    matrix(indices[4]+1,indices[1]+0) += tmp_13_3;
    matrix(indices[4]+1,indices[1]+1) += tmp_13_4;
    matrix(indices[4]+1,indices[1]+2) += tmp_13_5;
    matrix(indices[4]+1,indices[2]+0) += tmp_13_6;
    matrix(indices[4]+1,indices[2]+1) += tmp_13_7;
    matrix(indices[4]+1,indices[2]+2) += tmp_13_8;
    matrix(indices[4]+1,indices[3]+0) += tmp_13_9;
    matrix(indices[4]+1,indices[3]+1) += tmp_13_10;
    matrix(indices[4]+1,indices[3]+2) += tmp_13_11;
    matrix(indices[4]+1,indices[4]+0) += tmp_13_12;
    matrix(indices[4]+1,indices[4]+1) += tmp_13_13;
    matrix(indices[4]+1,indices[4]+2) += tmp_13_14;
    matrix(indices[4]+1,indices[5]+0) += tmp_13_15;
    matrix(indices[4]+1,indices[5]+1) += tmp_13_16;
    matrix(indices[4]+1,indices[5]+2) += tmp_13_17;
    matrix(indices[4]+2,indices[0]+0) += tmp_14_0;
    matrix(indices[4]+2,indices[0]+1) += tmp_14_1;
    matrix(indices[4]+2,indices[0]+2) += tmp_14_2;
    matrix(indices[4]+2,indices[1]+0) += tmp_14_3;
    matrix(indices[4]+2,indices[1]+1) += tmp_14_4;
    matrix(indices[4]+2,indices[1]+2) += tmp_14_5;
    matrix(indices[4]+2,indices[2]+0) += tmp_14_6;
    matrix(indices[4]+2,indices[2]+1) += tmp_14_7;
    matrix(indices[4]+2,indices[2]+2) += tmp_14_8;
    matrix(indices[4]+2,indices[3]+0) += tmp_14_9;
    matrix(indices[4]+2,indices[3]+1) += tmp_14_10;
    matrix(indices[4]+2,indices[3]+2) += tmp_14_11;
    matrix(indices[4]+2,indices[4]+0) += tmp_14_12;
    matrix(indices[4]+2,indices[4]+1) += tmp_14_13;
    matrix(indices[4]+2,indices[4]+2) += tmp_14_14;
    matrix(indices[4]+2,indices[5]+0) += tmp_14_15;
    matrix(indices[4]+2,indices[5]+1) += tmp_14_16;
    matrix(indices[4]+2,indices[5]+2) += tmp_14_17;
    matrix(indices[5]+0,indices[0]+0) += tmp_15_0;
    matrix(indices[5]+0,indices[0]+1) += tmp_15_1;
    matrix(indices[5]+0,indices[0]+2) += tmp_15_2;
    matrix(indices[5]+0,indices[1]+0) += tmp_15_3;
    matrix(indices[5]+0,indices[1]+1) += tmp_15_4;
    matrix(indices[5]+0,indices[1]+2) += tmp_15_5;
    matrix(indices[5]+0,indices[2]+0) += tmp_15_6;
    matrix(indices[5]+0,indices[2]+1) += tmp_15_7;
    matrix(indices[5]+0,indices[2]+2) += tmp_15_8;
    matrix(indices[5]+0,indices[3]+0) += tmp_15_9;
    matrix(indices[5]+0,indices[3]+1) += tmp_15_10;
    matrix(indices[5]+0,indices[3]+2) += tmp_15_11;
    matrix(indices[5]+0,indices[4]+0) += tmp_15_12;
    matrix(indices[5]+0,indices[4]+1) += tmp_15_13;
    matrix(indices[5]+0,indices[4]+2) += tmp_15_14;
    matrix(indices[5]+0,indices[5]+0) += tmp_15_15;
    matrix(indices[5]+0,indices[5]+1) += tmp_15_16;
    matrix(indices[5]+0,indices[5]+2) += tmp_15_17;
    matrix(indices[5]+1,indices[0]+0) += tmp_16_0;
    matrix(indices[5]+1,indices[0]+1) += tmp_16_1;
    matrix(indices[5]+1,indices[0]+2) += tmp_16_2;
    matrix(indices[5]+1,indices[1]+0) += tmp_16_3;
    matrix(indices[5]+1,indices[1]+1) += tmp_16_4;
    matrix(indices[5]+1,indices[1]+2) += tmp_16_5;
    matrix(indices[5]+1,indices[2]+0) += tmp_16_6;
    matrix(indices[5]+1,indices[2]+1) += tmp_16_7;
    matrix(indices[5]+1,indices[2]+2) += tmp_16_8;
    matrix(indices[5]+1,indices[3]+0) += tmp_16_9;
    matrix(indices[5]+1,indices[3]+1) += tmp_16_10;
    matrix(indices[5]+1,indices[3]+2) += tmp_16_11;
    matrix(indices[5]+1,indices[4]+0) += tmp_16_12;
    matrix(indices[5]+1,indices[4]+1) += tmp_16_13;
    matrix(indices[5]+1,indices[4]+2) += tmp_16_14;
    matrix(indices[5]+1,indices[5]+0) += tmp_16_15;
    matrix(indices[5]+1,indices[5]+1) += tmp_16_16;
    matrix(indices[5]+1,indices[5]+2) += tmp_16_17;
    matrix(indices[5]+2,indices[0]+0) += tmp_17_0;
    matrix(indices[5]+2,indices[0]+1) += tmp_17_1;
    matrix(indices[5]+2,indices[0]+2) += tmp_17_2;
    matrix(indices[5]+2,indices[1]+0) += tmp_17_3;
    matrix(indices[5]+2,indices[1]+1) += tmp_17_4;
    matrix(indices[5]+2,indices[1]+2) += tmp_17_5;
    matrix(indices[5]+2,indices[2]+0) += tmp_17_6;
    matrix(indices[5]+2,indices[2]+1) += tmp_17_7;
    matrix(indices[5]+2,indices[2]+2) += tmp_17_8;
    matrix(indices[5]+2,indices[3]+0) += tmp_17_9;
    matrix(indices[5]+2,indices[3]+1) += tmp_17_10;
    matrix(indices[5]+2,indices[3]+2) += tmp_17_11;
    matrix(indices[5]+2,indices[4]+0) += tmp_17_12;
    matrix(indices[5]+2,indices[4]+1) += tmp_17_13;
    matrix(indices[5]+2,indices[4]+2) += tmp_17_14;
    matrix(indices[5]+2,indices[5]+0) += tmp_17_15;
    matrix(indices[5]+2,indices[5]+1) += tmp_17_16;
    matrix(indices[5]+2,indices[5]+2) += tmp_17_17;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1-var_inter[0]; reg0=reg0-var_inter[1]; T reg1=var_inter[0]*elem.pos(1)[1]; T reg2=reg0*elem.pos(0)[1]; T reg3=var_inter[0]*elem.pos(1)[2];
    T reg4=reg0*elem.pos(0)[2]; T reg5=var_inter[1]*elem.pos(2)[1]; T reg6=var_inter[1]*elem.pos(2)[2]; T reg7=reg1+reg2; T reg8=reg3+reg4;
    T reg9=1-var_inter[2]; T reg10=reg9*elem.pos(0)[2]; T reg11=reg9*elem.pos(1)[2]; T reg12=elem.pos(3)[1]*reg0; T reg13=reg8+reg6;
    T reg14=elem.pos(3)[2]*reg0; T reg15=reg5+reg7; T reg16=reg9*elem.pos(2)[1]; T reg17=reg9*elem.pos(0)[1]; T reg18=reg9*elem.pos(1)[1];
    T reg19=reg9*elem.pos(2)[2]; T reg20=var_inter[0]*elem.pos(1)[0]; T reg21=elem.pos(0)[0]*reg0; reg18=reg18-reg17; T reg22=var_inter[2]*elem.pos(3)[1];
    reg12=reg12-reg15; T reg23=var_inter[0]*elem.pos(4)[1]; reg19=reg19-reg10; reg16=reg16-reg17; T reg24=var_inter[0]*elem.pos(4)[2];
    reg14=reg14-reg13; T reg25=var_inter[2]*elem.pos(3)[2]; reg11=reg11-reg10; T reg26=var_inter[2]*elem.pos(4)[1]; reg18=reg18-reg22;
    reg16=reg16-reg22; T reg27=var_inter[2]*elem.pos(4)[2]; reg11=reg11-reg25; T reg28=reg9*elem.pos(2)[0]; T reg29=var_inter[2]*elem.pos(5)[2];
    T reg30=var_inter[2]*elem.pos(5)[1]; T reg31=1+(*f.m).poisson_ratio; reg24=reg14+reg24; reg14=var_inter[1]*elem.pos(5)[2]; reg23=reg12+reg23;
    reg12=var_inter[1]*elem.pos(5)[1]; T reg32=var_inter[1]*elem.pos(2)[0]; T reg33=reg20+reg21; T reg34=reg9*elem.pos(1)[0]; T reg35=reg9*elem.pos(0)[0];
    reg19=reg19-reg25; T reg36=reg33+reg32; T reg37=elem.pos(3)[0]*reg0; reg23=reg12+reg23; reg19=reg29+reg19;
    reg16=reg30+reg16; reg24=reg14+reg24; reg31=reg31/(*f.m).elastic_modulus; reg28=reg28-reg35; reg27=reg11+reg27;
    reg26=reg18+reg26; reg11=var_inter[2]*elem.pos(3)[0]; reg34=reg34-reg35; reg12=pow(reg31,2); reg14=reg27*reg23;
    reg18=reg19*reg23; reg29=reg26*reg24; reg30=reg16*reg24; T reg38=var_inter[0]*elem.pos(4)[0]; reg37=reg37-reg36;
    T reg39=var_inter[2]*elem.pos(5)[0]; T reg40=var_inter[2]*elem.pos(4)[0]; reg28=reg28-reg11; reg34=reg34-reg11; T reg41=1.0/(*f.m).elastic_modulus;
    T reg42=(*f.m).poisson_ratio/(*f.m).elastic_modulus; reg28=reg39+reg28; reg39=reg27*reg16; reg31=reg31*reg12; T reg43=reg26*reg19;
    T reg44=var_inter[1]*elem.pos(5)[0]; reg14=reg29-reg14; reg38=reg37+reg38; reg40=reg34+reg40; reg18=reg30-reg18;
    reg29=reg40*reg18; reg30=reg28*reg14; reg39=reg43-reg39; reg34=reg42*reg12; reg12=reg41*reg12;
    reg38=reg44+reg38; reg37=reg42*reg31; reg31=reg41*reg31; reg43=reg40*reg23; reg44=reg38*reg39;
    T reg45=reg42*reg12; reg30=reg29-reg30; reg29=reg28*reg24; T reg46=reg41*reg31; T reg47=reg19*reg38;
    reg23=reg28*reg23; T reg48=reg42*reg37; reg31=reg42*reg31; reg24=reg40*reg24; T reg49=reg42*reg34;
    reg12=reg41*reg12; T reg50=reg16*reg38; T reg51=reg27*reg28; reg19=reg40*reg19; T reg52=reg26*reg38;
    reg38=reg27*reg38; reg51=reg19-reg51; reg16=reg40*reg16; reg37=reg41*reg37; reg31=reg48+reg31;
    reg34=reg41*reg34; reg46=reg46-reg48; reg28=reg26*reg28; reg45=reg49+reg45; reg52=reg43-reg52;
    reg38=reg24-reg38; reg50=reg23-reg50; reg44=reg30+reg44; reg12=reg12-reg49; reg47=reg29-reg47;
    reg50=reg50/reg44; reg37=reg48+reg37; reg28=reg16-reg28; reg12=reg41*reg12; reg51=reg51/reg44;
    reg18=reg18/reg44; reg41=reg41*reg46; reg38=reg38/reg44; reg16=reg49+reg34; reg47=reg47/reg44;
    reg19=reg42*reg31; reg39=reg39/reg44; reg45=reg42*reg45; reg52=reg52/reg44; reg14=reg14/reg44;
    reg16=reg42*reg16; reg45=reg12-reg45; reg12=reg9*reg14; reg23=reg9*reg18; reg24=var_inter[2]*reg52;
    reg26=var_inter[2]*reg50; reg27=var_inter[2]*reg47; reg29=var_inter[2]*reg38; reg30=var_inter[2]*reg14; reg40=var_inter[2]*reg18;
    reg43=reg9*reg38; reg48=reg9*reg47; reg42=reg42*reg37; reg19=reg41-reg19; reg41=var_inter[0]*reg51;
    T reg53=var_inter[0]*reg39; reg28=reg28/reg44; T reg54=var_inter[1]*reg51; T reg55=var_inter[0]*reg28; T reg56=var_inter[1]*reg39;
    T reg57=reg48-reg43; T reg58=reg0*reg28; T reg59=reg24-reg26; T reg60=reg9*reg50; T reg61=reg27-reg29;
    T reg62=reg0*reg51; T reg63=reg30-reg40; T reg64=var_inter[1]*reg28; reg16=reg45-reg16; reg45=reg12-reg23;
    T reg65=reg9*reg52; T reg66=reg0*reg39; T reg67=reg27+reg41; reg42=reg19-reg42; reg19=reg40+reg53;
    T reg68=reg41-reg48; T reg69=reg55+reg26; reg57=reg57+reg62; T reg70=0.5*reg67; reg63=reg63+reg66;
    reg61=reg61-reg62; reg59=reg59+reg58; T reg71=reg12+reg56; T reg72=reg54+reg43; T reg73=reg64+reg65;
    T reg74=reg65-reg60; T reg75=reg23-reg53; reg45=reg45-reg66; T reg76=reg64-reg24; reg16=reg16/reg42;
    T reg77=(*f.m).deltaT*(*f.m).alpha; reg37=reg37/reg42; reg46=reg46/reg42; reg42=reg31/reg42; reg31=reg29-reg54;
    T reg78=reg56-reg30; T reg79=0.5*reg19; T reg80=0.5*reg45; T reg81=reg60-reg55; T reg82=0.5*reg68;
    T reg83=0.5*reg71; T reg84=0.5*reg69; T reg85=0.5*reg78; T reg86=0.5*reg63; T reg87=reg46*reg77;
    T reg88=0.5*reg75; T reg89=0.5*reg76; T reg90=0.5*reg31; reg74=reg74-reg58; T reg91=0.5*reg61;
    T reg92=0.5*reg57; T reg93=0.5*reg59; T reg94=0.5*reg73; T reg95=0.5*reg72; T reg96=reg16*reg70;
    T reg97=reg16*reg79; T reg98=reg37*reg77; T reg99=reg42*reg77; T reg100=reg16*reg86; T reg101=reg16*reg82;
    T reg102=reg16*reg89; T reg103=0.5*reg74; T reg104=0.5*reg81; T reg105=reg16*reg90; T reg106=reg16*reg95;
    T reg107=reg99+reg98; T reg108=reg16*reg94; T reg109=reg46*reg67; T reg110=reg16*reg83; T reg111=reg46*reg69;
    T reg112=reg16*reg85; T reg113=reg87+reg99; T reg114=reg46*reg19; T reg115=reg16*reg93; T reg116=reg16*reg92;
    T reg117=reg16*reg84; T reg118=reg16*reg91; T reg119=reg16*reg88; reg97=2*reg97; T reg120=reg16*reg80;
    T reg121=2*reg96; reg117=2*reg117; T reg122=reg42*reg67; reg100=2*reg100; reg105=2*reg105;
    T reg123=reg37*reg59; T reg124=reg46*reg78; T reg125=var_inter[1]*reg9; T reg126=reg37*reg69; T reg127=var_inter[0]*var_inter[2];
    reg115=2*reg115; T reg128=reg37*reg76; reg118=2*reg118; T reg129=reg46*reg57; T reg130=reg46*reg63;
    T reg131=reg37*reg73; T reg132=reg73*reg111; T reg133=reg16*reg103; reg112=2*reg112; T reg134=reg46*reg68;
    reg116=2*reg116; T reg135=reg46*reg45; T reg136=reg42*reg63; T reg137=reg46*reg61; T reg138=2*reg110;
    T reg139=reg42*reg72; T reg140=reg42*reg19; T reg141=2*reg108; reg106=2*reg106; T reg142=reg83*reg97;
    T reg143=reg72*reg109; T reg144=reg46*reg71; T reg145=reg42*reg78; reg102=2*reg102; T reg146=reg46*reg31;
    T reg147=reg46*reg81; T reg148=reg113+reg98; T reg149=reg46*reg75; reg101=2*reg101; T reg150=reg46*reg72;
    T reg151=reg87+reg107; T reg152=reg42*reg71; T reg153=reg46*reg73; T reg154=reg16*reg104; T reg155=reg95*reg121;
    T reg156=reg46*reg74; reg120=2*reg120; T reg157=reg46*reg59; T reg158=reg71*reg114; T reg159=reg46*reg76;
    reg119=2*reg119; T reg160=reg120*reg80; T reg161=reg91*reg118; T reg162=reg63*reg130; T reg163=reg57*reg129;
    T reg164=reg71*reg144; T reg165=reg73*reg159; T reg166=reg83*reg102; T reg167=reg74*reg153; T reg168=reg73*reg145;
    T reg169=reg42*reg75; reg132=reg142+reg132; T reg170=reg95*reg138; T reg171=reg71*reg139; T reg172=reg119*reg80;
    T reg173=reg37*reg72; T reg174=reg83*reg117; T reg175=reg74*reg159; T reg176=reg73*reg140; T reg177=reg68*reg146;
    T reg178=reg88*reg112; T reg179=reg37*reg31; T reg180=reg81*reg147; T reg181=reg152*reg81; T reg182=reg88*reg141;
    T reg183=reg81*reg153; T reg184=reg86*reg100; T reg185=reg61*reg137; T reg186=reg74*reg111; T reg187=reg45*reg124;
    T reg188=reg81*reg157; T reg189=reg37*reg67; T reg190=reg19*reg124; T reg191=reg92*reg105; T reg192=reg91*reg105;
    T reg193=reg70*reg105; T reg194=reg63*reg124; T reg195=reg74*reg157; T reg196=reg91*reg121; T reg197=reg81*reg111;
    T reg198=reg63*reg114; T reg199=reg37*reg61; T reg200=reg81*reg159; T reg201=reg95*reg106; T reg202=reg57*reg146;
    T reg203=reg80*reg112; T reg204=reg72*reg131; T reg205=reg82*reg121; T reg206=reg75*reg114; T reg207=reg82*reg105;
    T reg208=reg83*reg138; T reg209=reg150*reg72; T reg210=reg94*reg97; T reg211=reg94*reg112; T reg212=reg95*reg105;
    T reg213=reg88*reg100; T reg214=reg68*reg137; T reg215=reg75*reg124; T reg216=reg75*reg130; T reg217=reg82*reg106;
    T reg218=reg75*reg144; T reg219=reg71*reg124; T reg220=reg88*reg138; T reg221=reg150*reg68; T reg222=reg119*reg88;
    T reg223=reg68*reg134; T reg224=reg71*reg128; T reg225=reg75*reg131; T reg226=reg82*reg118; T reg227=reg104*reg138;
    T reg228=reg149*reg75; T reg229=reg73*reg157; T reg230=reg88*reg97; T reg231=reg95*reg118; T reg232=reg80*reg100;
    T reg233=reg83*reg115; T reg234=reg73*reg136; T reg235=reg71*reg130; T reg236=reg80*reg141; T reg237=reg73*reg153;
    T reg238=reg152*reg74; T reg239=reg83*reg112; T reg240=reg72*reg146; T reg241=reg71*reg123; T reg242=reg74*reg147;
    T reg243=reg94*reg100; T reg244=reg57*reg109; T reg245=reg80*reg97; reg142=reg143+reg142; reg158=reg155+reg158;
    T reg246=reg37*reg68; T reg247=reg83*reg100; T reg248=reg94*reg117; T reg249=reg71*reg126; T reg250=reg72*reg137;
    T reg251=reg94*reg106; T reg252=reg68*reg109; T reg253=reg57*reg137; T reg254=reg69*reg111; T reg255=reg125*(*f.m).f_vol[0];
    T reg256=reg67*reg146; T reg257=reg79*reg112; T reg258=reg37*reg81; T reg259=reg84*reg121; T reg260=reg67*reg126;
    T reg261=reg67*reg109; T reg262=reg74*reg156; T reg263=reg79*reg97; T reg264=reg45*reg144; T reg265=reg92*reg106;
    T reg266=reg70*reg97; T reg267=reg19*reg122; T reg268=reg70*reg121; T reg269=reg19*reg114; T reg270=reg42*reg31;
    T reg271=reg45*reg135; T reg272=reg92*reg116; T reg273=reg59*reg159; reg133=2*reg133; T reg274=reg45*reg131;
    T reg275=reg103*reg138; T reg276=reg9*reg0; T reg277=reg9*var_inter[0]; T reg278=var_inter[1]*var_inter[2]; T reg279=var_inter[2]*reg0;
    T reg280=reg57*reg134; T reg281=reg150*reg57; T reg282=reg101*reg82; T reg283=reg80*reg138; T reg284=reg67*reg148;
    T reg285=reg127*(*f.m).f_vol[1]; T reg286=reg125*(*f.m).f_vol[2]; T reg287=reg73*reg151; T reg288=reg71*reg148; T reg289=reg42*reg57;
    T reg290=reg37*reg74; T reg291=reg149*reg45; T reg292=reg101*reg92; reg154=2*reg154; T reg293=reg76*reg159;
    T reg294=reg42*reg68; T reg295=reg85*reg112; T reg296=reg31*reg146; T reg297=reg90*reg105; reg124=reg78*reg124;
    reg159=reg69*reg159; reg146=reg61*reg146; T reg298=reg86*reg112; T reg299=reg42*reg61; T reg300=reg59*reg157;
    T reg301=reg92*reg118; T reg302=reg45*reg114; T reg303=reg92*reg121; T reg304=reg86*reg97; T reg305=reg45*reg130;
    T reg306=reg59*reg111; T reg307=reg61*reg109; T reg308=reg81*reg140; reg188=reg213+reg188; T reg309=reg57*reg148;
    T reg310=reg91*reg112; T reg311=reg104*reg117; T reg312=reg82*reg115; T reg313=reg81*reg199; T reg314=reg88*reg115;
    T reg315=reg81*reg136; T reg316=reg74*reg151; T reg317=reg75*reg148; T reg318=reg220+reg183; T reg319=reg68*reg148;
    T reg320=reg82*reg141; T reg321=reg81*reg173; T reg322=reg181+reg182; T reg323=reg81*reg151; T reg324=reg63*reg128;
    reg180=reg222+reg180; T reg325=reg93*reg121; T reg326=reg288-reg255; T reg327=reg104*reg105; T reg328=reg68*reg128;
    T reg329=reg86*reg121; reg296=reg296+reg295; T reg330=reg201+reg164; T reg331=reg31*reg128; T reg332=reg89*reg105;
    reg200=reg178+reg200; T reg333=reg86*reg105; T reg334=reg82*reg102; T reg335=reg81*reg179; T reg336=reg88*reg102;
    T reg337=reg61*reg145; T reg338=reg81*reg145; reg197=reg230+reg197; reg293=reg295+reg293; reg295=reg82*reg117;
    T reg339=reg81*reg189; reg194=reg194+reg192; T reg340=reg93*reg102; T reg341=reg84*reg102; reg190=reg190-reg193;
    T reg342=reg63*reg270; T reg343=reg45*reg148; reg206=reg206-reg205; T reg344=reg84*reg97; T reg345=reg284-reg285;
    reg221=reg221-reg220; T reg346=reg69*reg151; T reg347=reg78*reg148; T reg348=reg82*reg97; T reg349=reg88*reg106;
    T reg350=reg152*reg68; T reg351=reg101*reg104; T reg352=reg68*reg258; reg222=reg223+reg222; reg223=reg31*reg148;
    T reg353=reg75*reg122; T reg354=reg75*reg126; T reg355=reg61*reg123; T reg356=reg104*reg112; T reg357=reg75*reg128;
    T reg358=reg75*reg270; T reg359=reg82*reg112; T reg360=reg104*reg97; T reg361=reg61*reg140; T reg362=reg76*reg151;
    T reg363=reg104*reg102; reg215=reg207+reg215; T reg364=reg93*reg118; reg178=reg177+reg178; reg177=reg72*reg148;
    T reg365=reg93*reg112; T reg366=reg287-reg286; T reg367=reg63*reg148; T reg368=reg88*reg105; T reg369=reg68*reg145;
    T reg370=reg104*reg121; T reg371=reg61*reg126; T reg372=reg68*reg126; reg230=reg230-reg252; T reg373=reg61*reg148;
    T reg374=reg59*reg151; T reg375=reg88*reg121; T reg376=reg68*reg140; T reg377=reg104*reg118; T reg378=reg68*reg123;
    reg213=reg214+reg213; reg214=reg19*reg148; T reg379=reg304-reg307; T reg380=reg88*reg118; T reg381=reg68*reg136;
    T reg382=reg104*reg106; T reg383=reg68*reg131; reg185=reg185+reg184; T reg384=reg63*reg123; T reg385=reg208+reg237;
    T reg386=reg94*reg105; T reg387=reg72*reg128; reg240=reg240-reg239; reg266=reg267+reg266; T reg388=reg19*reg126;
    T reg389=reg19*reg270; T reg390=reg83*reg105; T reg391=reg72*reg145; T reg392=reg94*reg121; T reg393=reg72*reg126;
    T reg394=reg248+reg142; T reg395=reg70*reg112; T reg396=reg19*reg128; T reg397=reg93*reg100; T reg398=reg59*reg189;
    T reg399=reg84*reg112; T reg400=reg86*reg117; T reg401=reg83*reg121; T reg402=reg72*reg140; T reg403=reg94*reg118;
    T reg404=reg72*reg123; reg250=reg250-reg247; reg198=reg198-reg196; reg165=reg239+reg165; reg239=reg59*reg145;
    T reg405=reg73*reg179; T reg406=reg95*reg102; reg166=reg168+reg166; reg306=reg304+reg306; reg168=reg86*reg102;
    reg304=reg59*reg179; reg132=reg155+reg132; T reg407=reg91*reg102; T reg408=reg73*reg189; T reg409=reg95*reg117;
    T reg410=reg91*reg117; reg174=reg176+reg174; reg162=reg162+reg161; reg176=reg93*reg115; reg273=reg298+reg273;
    reg229=reg247+reg229; reg247=reg63*reg299; T reg411=reg73*reg199; T reg412=reg95*reg115; reg233=reg234+reg233;
    reg234=reg91*reg100; reg269=reg269+reg268; T reg413=reg84*reg117; T reg414=reg95*reg97; reg248=reg158+reg248;
    T reg415=reg69*reg145; T reg416=reg79*reg102; T reg417=reg69*reg179; T reg418=reg70*reg102; reg243=reg241+reg243;
    T reg419=reg93*reg105; reg159=reg257+reg159; T reg420=reg61*reg128; reg298=reg146+reg298; reg146=reg71*reg299;
    T reg421=reg95*reg100; T reg422=reg94*reg115; reg124=reg124+reg297; T reg423=reg89*reg102; reg235=reg231-reg235;
    T reg424=reg78*reg270; T reg425=reg90*reg112; T reg426=reg94*reg138; T reg427=reg78*reg128; T reg428=reg71*reg131;
    reg171=reg170+reg171; T reg429=reg89*reg112; T reg430=reg94*reg141; T reg431=reg59*reg140; T reg432=reg263+reg261;
    T reg433=reg83*reg118; T reg434=reg72*reg136; reg251=reg204+reg251; reg300=reg184+reg300; reg184=reg93*reg117;
    T reg435=reg260+reg259; T reg436=reg63*reg122; T reg437=reg79*reg105; reg209=reg209+reg208; T reg438=reg67*reg145;
    T reg439=reg91*reg97; reg211=reg224+reg211; T reg440=reg63*reg126; T reg441=reg93*reg97; T reg442=reg71*reg270;
    T reg443=reg95*reg112; reg256=reg257-reg256; reg257=reg94*reg102; reg219=reg212-reg219; T reg444=reg67*reg128;
    T reg445=reg84*reg105; reg210=reg249+reg210; reg254=reg263+reg254; reg263=reg71*reg122; T reg446=reg154*reg80;
    T reg447=reg74*reg169; T reg448=reg88*reg117; T reg449=reg290*reg45; T reg450=reg103*reg105; T reg451=reg57*reg128;
    T reg452=reg120*reg103; reg202=reg202+reg203; T reg453=reg80*reg105; T reg454=reg57*reg145; T reg455=reg103*reg121;
    reg291=reg291+reg292; T reg456=reg57*reg126; T reg457=reg154*reg103; T reg458=reg245-reg244; T reg459=reg80*reg121;
    T reg460=reg57*reg140; T reg461=reg103*reg118; T reg462=reg57*reg123; T reg463=reg294*reg45; reg253=reg253+reg232;
    T reg464=reg119*reg92; T reg465=reg80*reg118; T reg466=reg57*reg136; T reg467=reg101*reg80; T reg468=reg278*(*f.m).f_vol[0];
    T reg469=reg92*reg102; T reg470=reg74*reg179; T reg471=reg125*(*f.m).f_vol[1]; T reg472=reg80*reg102; T reg473=reg74*reg145;
    reg186=reg245+reg186; reg281=reg281-reg283; reg245=reg92*reg117; T reg474=reg74*reg189; T reg475=reg80*reg117;
    T reg476=reg80*reg106; T reg477=reg92*reg115; T reg478=reg74*reg199; T reg479=reg152*reg57; T reg480=reg80*reg115;
    T reg481=reg74*reg136; T reg482=reg92*reg141; T reg483=reg289*reg45; T reg484=reg74*reg173; T reg485=reg238+reg236;
    T reg486=reg120*reg92; reg242=reg172+reg242; T reg487=reg154*reg92; T reg488=reg74*reg246; reg195=reg232+reg195;
    reg232=reg92*reg100; T reg489=reg45*reg299; reg175=reg203+reg175; reg203=reg103*reg115; reg305=reg305+reg301;
    T reg490=reg277*(*f.m).f_vol[0]; T reg491=reg274+reg275; T reg492=reg103*reg133; reg271=reg271+reg272; T reg493=reg278*(*f.m).f_vol[2];
    T reg494=reg277*(*f.m).f_vol[1]; T reg495=reg278*(*f.m).f_vol[1]; T reg496=reg279*(*f.m).f_vol[0]; T reg497=reg92*reg138; T reg498=reg45*reg139;
    T reg499=reg277*(*f.m).f_vol[2]; T reg500=reg103*reg141; T reg501=reg276*(*f.m).f_vol[0]; T reg502=reg265-reg264; T reg503=reg276*(*f.m).f_vol[1];
    reg262=reg262+reg160; T reg504=reg119*reg103; T reg505=reg45*reg258; T reg506=reg276*(*f.m).f_vol[2]; T reg507=reg57*reg169;
    T reg508=reg279*(*f.m).f_vol[1]; T reg509=reg103*reg116; T reg510=reg290*reg57; T reg511=reg279*(*f.m).f_vol[2]; reg160=reg163+reg160;
    reg163=reg127*(*f.m).f_vol[0]; T reg512=reg103*reg112; reg128=reg45*reg128; reg112=reg92*reg112; reg270=reg45*reg270;
    T reg513=reg74*reg140; T reg514=reg103*reg102; reg187=reg187+reg191; T reg515=reg103*reg97; T reg516=reg45*reg126;
    T reg517=reg101*reg103; T reg518=reg92*reg97; T reg519=reg57*reg131; T reg520=reg45*reg122; T reg521=reg103*reg117;
    reg302=reg302-reg303; T reg522=reg283+reg167; reg172=reg280+reg172; reg280=reg103*reg100; T reg523=reg45*reg123;
    T reg524=reg217-reg218; T reg525=reg104*reg141; reg228=reg282+reg228; T reg526=reg103*reg106; reg216=reg226+reg216;
    T reg527=reg104*reg115; T reg528=reg119*reg82; T reg529=reg127*(*f.m).f_vol[2]; T reg530=reg294*reg75; T reg531=reg225+reg227;
    T reg532=reg82*reg100; T reg533=reg119*reg104; T reg534=reg75*reg299; T reg535=reg57*reg258; T reg536=reg75*reg139;
    T reg537=reg154*reg104; T reg538=reg82*reg138; T reg539=reg75*reg258; T reg540=reg104*reg100; T reg541=reg75*reg123;
    reg512=reg128+reg512; reg128=reg44*reg233; reg465=reg466+reg465; reg162=reg162+reg176; reg411=reg412-reg411;
    reg524=reg524-reg525; reg425=reg424+reg425; reg412=reg529+reg346; reg229=reg231-reg229; reg296=reg423+reg296;
    reg351=reg352+reg351; reg464=reg463+reg464; reg526=reg526-reg519; reg160=reg492+reg160; reg231=reg44*reg174;
    reg281=reg281-reg500; reg409=reg409+reg408; reg165=reg212-reg165; reg467=reg507+reg467; reg349=reg349-reg350;
    reg405=reg406-reg405; reg212=reg44*reg166; reg352=reg44*reg132; reg429=reg427+reg429; reg509=reg510+reg509;
    reg450=reg451+reg450; reg406=reg44*reg211; reg424=reg503+reg309; reg202=reg514+reg202; reg209=reg430+reg209;
    reg345=reg44*reg345; reg427=reg44*reg251; reg451=reg501+reg343; reg453=reg454+reg453; reg433=reg434-reg433;
    reg452=reg449+reg452; reg380=reg381+reg380; reg250=reg250-reg422; reg456=reg456-reg455; reg403=reg404-reg403;
    reg533=reg539+reg533; reg253=reg203+reg253; reg221=reg221-reg525; reg201=reg201+reg385; reg332=reg331+reg332;
    reg386=reg387-reg386; reg172=reg457+reg172; reg240=reg240-reg257; reg461=reg462+reg461; reg390=reg391-reg390;
    reg293=reg297+reg293; reg393=reg393+reg392; reg460=reg460-reg459; reg297=reg44*reg394; reg382=reg382-reg383;
    reg458=reg521+reg458; reg457=reg291+reg457; reg402=reg402+reg401; reg273=reg192+reg273; reg445=reg445-reg444;
    reg407=reg304+reg407; reg492=reg271+reg492; reg168=reg239+reg168; reg360=reg354+reg360; reg306=reg306-reg196;
    reg192=reg44*reg491; reg410=reg410-reg398; reg254=reg268+reg254; reg400=reg431+reg400; reg203=reg305+reg203;
    reg300=reg161+reg300; reg216=reg216+reg527; reg419=reg420+reg419; reg215=reg215+reg363; reg540=reg541+reg540;
    reg438=reg437-reg438; reg161=reg44*reg435; reg504=reg505+reg504; reg432=reg413+reg432; reg399=reg396+reg399;
    reg272=reg262+reg272; reg395=reg389-reg395; reg206=reg206+reg311; reg502=reg502-reg500; reg256=reg341+reg256;
    reg239=reg44*reg266; reg262=reg493+reg362; reg534=reg532+reg534; reg413=reg269+reg413; reg498=reg498-reg497;
    reg348=reg348-reg353; reg365=reg324+reg365; reg193=reg159-reg193; reg310=reg342+reg310; reg515=reg516+reg515;
    reg194=reg194+reg340; reg159=reg468+reg347; reg441=reg440+reg441; reg514=reg187+reg514; reg439=reg439-reg436;
    reg517=reg535+reg517; reg536=reg536-reg538; reg423=reg124+reg423; reg198=reg198+reg184; reg222=reg537+reg222;
    reg112=reg270+reg112; reg397=reg384+reg397; reg234=reg247+reg234; reg298=reg340+reg298; reg232=reg489+reg232;
    reg416=reg415+reg416; reg333=reg337+reg333; reg124=reg495+reg223; reg371=reg371-reg325; reg280=reg523+reg280;
    reg379=reg184+reg379; reg358=reg359+reg358; reg418=reg417-reg418; reg361=reg361-reg329; reg521=reg302+reg521;
    reg364=reg355+reg364; reg184=reg44*reg531; reg185=reg176+reg185; reg518=reg518-reg520; reg356=reg357+reg356;
    reg295=reg295-reg339; reg477=reg478+reg477; reg230=reg311+reg230; reg197=reg197-reg205; reg336=reg338+reg336;
    reg326=reg44*reg326; reg480=reg481+reg480; reg334=reg335+reg334; reg537=reg228+reg537; reg176=reg163+reg214;
    reg200=reg207+reg200; reg265=reg265-reg522; reg330=reg330+reg430; reg187=reg499+reg323; reg484=reg484-reg482;
    reg207=reg44*reg171; reg376=reg376-reg375; reg228=reg428+reg426; reg247=reg44*reg485; reg269=reg494+reg319;
    reg422=reg235-reg422; reg242=reg292+reg242; reg146=reg421-reg146; reg235=reg44*reg243; reg270=reg490+reg317;
    reg487=reg488+reg487; reg271=reg511+reg374; reg469=reg470+reg469; reg178=reg363+reg178; reg368=reg369+reg368;
    reg327=reg328+reg327; reg291=reg508+reg373; reg180=reg282+reg180; reg191=reg175+reg191; reg472=reg473+reg472;
    reg175=reg44*reg322; reg321=reg321-reg320; reg282=reg496+reg367; reg186=reg186-reg303; reg217=reg217-reg318;
    reg372=reg372-reg370; reg314=reg315+reg314; reg245=reg245-reg474; reg312=reg313+reg312; reg366=reg44*reg366;
    reg188=reg226+reg188; reg476=reg476-reg479; reg475=reg513+reg475; reg388=reg344+reg388; reg301=reg195+reg301;
    reg341=reg190+reg341; reg190=reg471+reg177; reg195=reg506+reg316; reg226=reg44*reg210; reg414=reg414+reg263;
    reg530=reg528+reg530; reg308=reg448+reg308; reg446=reg447+reg446; reg257=reg219-reg257; reg219=reg44*reg248;
    reg213=reg527+reg213; reg377=reg378+reg377; reg486=reg483+reg486; reg442=reg443-reg442; reg380=reg44*reg380;
    reg379=reg44*reg379; reg486=reg44*reg486; reg280=reg44*reg280; reg475=reg44*reg475; reg371=reg44*reg371;
    reg292=ponderation*reg161; reg188=reg44*reg188; reg215=reg44*reg215; reg302=ponderation*reg427; reg333=reg44*reg333;
    reg416=reg44*reg416; reg232=reg44*reg232; reg308=reg44*reg308; reg312=reg44*reg312; reg298=reg44*reg298;
    reg453=reg44*reg453; reg304=reg44*reg451; reg366=ponderation*reg366; reg419=reg44*reg419; reg314=reg44*reg314;
    reg203=reg44*reg203; reg254=reg44*reg254; reg300=reg44*reg300; reg305=reg44*reg159; reg456=reg44*reg456;
    reg310=reg44*reg310; reg193=reg44*reg193; reg356=reg44*reg356; reg452=reg44*reg452; reg230=reg44*reg230;
    reg365=reg44*reg365; reg295=reg44*reg295; reg518=reg44*reg518; reg272=reg44*reg272; reg477=reg44*reg477;
    reg185=reg44*reg185; reg433=reg44*reg433; reg311=ponderation*reg226; reg364=reg44*reg364; reg418=reg44*reg418;
    reg521=reg44*reg521; reg301=reg44*reg301; reg313=reg44*reg190; reg361=reg44*reg361; reg469=reg44*reg469;
    reg358=reg44*reg358; reg388=reg44*reg388; reg315=ponderation*reg184; reg257=reg44*reg257; reg324=reg44*reg271;
    reg498=reg44*reg498; reg191=reg44*reg191; reg438=reg44*reg438; reg413=reg44*reg413; reg328=ponderation*reg406;
    reg180=reg44*reg180; reg256=reg44*reg256; reg331=ponderation*reg239; reg472=reg44*reg472; reg502=reg44*reg502;
    reg341=reg44*reg341; reg206=reg44*reg206; reg335=reg44*reg291; reg327=reg44*reg327; reg395=reg44*reg395;
    reg534=reg44*reg534; reg337=reg44*reg262; reg530=reg44*reg530; reg450=reg44*reg450; reg399=reg44*reg399;
    reg368=reg44*reg368; reg504=reg44*reg504; reg178=reg44*reg178; reg432=reg44*reg432; reg245=reg44*reg245;
    reg216=reg44*reg216; reg442=reg44*reg442; reg400=reg44*reg400; reg338=reg44*reg124; reg172=reg44*reg172;
    reg217=reg44*reg217; reg410=reg44*reg410; reg360=reg44*reg360; reg340=ponderation*reg192; reg209=reg44*reg209;
    reg372=reg44*reg372; reg306=reg44*reg306; reg213=reg44*reg213; reg492=reg44*reg492; reg540=reg44*reg540;
    reg168=reg44*reg168; reg445=reg44*reg445; reg321=reg44*reg321; reg202=reg44*reg202; reg407=reg44*reg407;
    reg186=reg44*reg186; reg342=reg44*reg282; reg344=ponderation*reg175; reg273=reg44*reg273; reg348=reg44*reg348;
    reg354=reg44*reg424; reg355=ponderation*reg247; reg229=reg44*reg229; reg357=reg44*reg269; reg228=reg44*reg228;
    reg526=reg44*reg526; reg359=ponderation*reg219; reg363=ponderation*reg231; reg464=reg44*reg464; reg369=reg44*reg176;
    reg467=reg44*reg467; reg458=reg44*reg458; reg409=reg44*reg409; reg378=ponderation*reg207; reg349=reg44*reg349;
    reg376=reg44*reg376; reg381=ponderation*reg352; reg429=reg44*reg429; reg402=reg44*reg402; reg509=reg44*reg509;
    reg484=reg44*reg484; reg384=ponderation*reg212; reg524=reg44*reg524; reg387=reg44*reg412; reg330=reg44*reg330;
    reg405=reg44*reg405; reg393=reg44*reg393; reg533=reg44*reg533; reg390=reg44*reg390; reg457=reg44*reg457;
    reg461=reg44*reg461; reg389=ponderation*reg235; reg460=reg44*reg460; reg240=reg44*reg240; reg487=reg44*reg487;
    reg221=reg44*reg221; reg391=reg44*reg270; reg386=reg44*reg386; reg332=reg44*reg332; reg377=reg44*reg377;
    reg253=reg44*reg253; reg201=reg44*reg201; reg146=reg44*reg146; reg293=reg44*reg293; reg396=ponderation*reg297;
    reg242=reg44*reg242; reg404=ponderation*reg128; reg465=reg44*reg465; reg422=reg44*reg422; reg446=reg44*reg446;
    reg411=reg44*reg411; reg296=reg44*reg296; reg234=reg44*reg234; reg476=reg44*reg476; reg537=reg44*reg537;
    reg397=reg44*reg397; reg423=reg44*reg423; reg112=reg44*reg112; reg222=reg44*reg222; reg334=reg44*reg334;
    reg414=reg44*reg414; reg198=reg44*reg198; reg415=reg44*reg195; reg480=reg44*reg480; reg336=reg44*reg336;
    reg439=reg44*reg439; reg517=reg44*reg517; reg514=reg44*reg514; reg250=reg44*reg250; reg326=ponderation*reg326;
    reg441=reg44*reg441; reg536=reg44*reg536; reg281=reg44*reg281; reg515=reg44*reg515; reg197=reg44*reg197;
    reg194=reg44*reg194; reg403=reg44*reg403; reg351=reg44*reg351; reg417=reg44*reg187; reg200=reg44*reg200;
    reg345=ponderation*reg345; reg265=reg44*reg265; reg425=reg44*reg425; reg160=reg44*reg160; reg165=reg44*reg165;
    reg162=reg44*reg162; reg382=reg44*reg382; reg512=reg44*reg512; reg420=ponderation*reg391; sollicitation[indices[1]+0]+=reg420;
    reg421=ponderation*reg304; sollicitation[indices[0]+0]+=reg421; sollicitation[indices[2]+2]+=-reg366; reg366=ponderation*reg324; sollicitation[indices[3]+2]+=reg366;
    T tmp_16_17=ponderation*reg332; reg332=ponderation*reg335; sollicitation[indices[3]+1]+=reg332; T tmp_13_15=ponderation*reg438; T tmp_1_7=ponderation*reg281;
    reg281=ponderation*reg387; sollicitation[indices[4]+2]+=reg281; sollicitation[indices[2]+0]+=-reg326; T tmp_0_3=ponderation*reg457; reg326=ponderation*reg337;
    sollicitation[indices[5]+2]+=reg326; T tmp_15_16=ponderation*reg425; T tmp_14_14=ponderation*reg254; sollicitation[indices[4]+1]+=-reg345; T tmp_14_17=ponderation*reg193;
    T tmp_0_1=ponderation*reg486; T tmp_14_15=ponderation*reg416; reg193=ponderation*reg338; sollicitation[indices[5]+1]+=reg193; reg254=ponderation*reg357;
    sollicitation[indices[1]+1]+=reg254; T tmp_1_6=ponderation*reg476; T tmp_15_15=ponderation*reg423; reg345=ponderation*reg342; sollicitation[indices[3]+0]+=reg345;
    T tmp_0_0=ponderation*reg492; T tmp_13_17=ponderation*reg445; T tmp_15_17=ponderation*reg429; T tmp_17_17=ponderation*reg293; reg293=ponderation*reg313;
    sollicitation[indices[2]+1]+=reg293; reg416=ponderation*reg417; sollicitation[indices[1]+2]+=reg416; reg423=ponderation*reg415; sollicitation[indices[0]+2]+=reg423;
    T tmp_13_16=ponderation*reg256; T tmp_16_16=ponderation*reg296; T tmp_1_4=ponderation*reg172; reg172=ponderation*reg369; sollicitation[indices[4]+0]+=reg172;
    reg256=ponderation*reg354; sollicitation[indices[0]+1]+=reg256; T tmp_0_2=ponderation*reg452; T tmp_14_16=ponderation*reg418; reg296=ponderation*reg305;
    sollicitation[indices[5]+0]+=reg296; T tmp_1_5=ponderation*reg517; T tmp_2_6=-reg355; T tmp_6_7=-reg378; T tmp_6_6=ponderation*reg330;
    T tmp_2_7=ponderation*reg484; T tmp_2_8=ponderation*reg265; T tmp_5_17=ponderation*reg200; T tmp_5_16=ponderation*reg334; T tmp_5_15=ponderation*reg336;
    T tmp_2_9=ponderation*reg480; T tmp_5_14=ponderation*reg197; T tmp_5_13=ponderation*reg295; T tmp_2_10=ponderation*reg477; T tmp_2_11=ponderation*reg301;
    T tmp_12_14=ponderation*reg388; T tmp_5_11=ponderation*reg188; T tmp_2_12=ponderation*reg475; T tmp_5_10=ponderation*reg312; T tmp_5_9=ponderation*reg314;
    T tmp_1_15=ponderation*reg453; T tmp_7_7=ponderation*reg209; T tmp_1_16=ponderation*reg202; T tmp_6_17=-reg328; T tmp_6_16=ponderation*reg442;
    T tmp_1_17=ponderation*reg450; T tmp_6_15=ponderation*reg257; T tmp_5_12=ponderation*reg308; T tmp_6_14=-reg311; T tmp_2_2=ponderation*reg272;
    T tmp_6_13=ponderation*reg414; T tmp_6_12=-reg359; T tmp_2_3=ponderation*reg446; T tmp_6_11=-reg389; T tmp_2_4=ponderation*reg487;
    T tmp_6_10=ponderation*reg146; T tmp_6_9=ponderation*reg422; T tmp_2_5=ponderation*reg242; T tmp_6_8=ponderation*reg228; T tmp_4_9=ponderation*reg380;
    T tmp_4_8=ponderation*reg382; T tmp_3_5=ponderation*reg533; T tmp_4_7=ponderation*reg221; T tmp_4_6=ponderation*reg349; T tmp_3_6=ponderation*reg524;
    T tmp_4_5=ponderation*reg351; T tmp_4_4=ponderation*reg222; T tmp_3_7=ponderation*reg536; T tmp_3_17=ponderation*reg356; T tmp_3_16=ponderation*reg358;
    T tmp_3_8=-reg315; T tmp_3_15=ponderation*reg215; T tmp_3_14=ponderation*reg360; T tmp_3_9=ponderation*reg216; T tmp_3_13=ponderation*reg348;
    T tmp_3_12=ponderation*reg206; T tmp_3_10=ponderation*reg534; T tmp_3_11=ponderation*reg540; T tmp_2_13=ponderation*reg245; T tmp_5_8=ponderation*reg217;
    T tmp_5_7=ponderation*reg321; T tmp_2_14=ponderation*reg186; T tmp_5_6=-reg344; T tmp_5_5=ponderation*reg180; T tmp_2_15=ponderation*reg472;
    T tmp_4_17=ponderation*reg327; T tmp_4_16=ponderation*reg178; T tmp_2_16=ponderation*reg469; T tmp_4_15=ponderation*reg368; T tmp_4_14=ponderation*reg372;
    T tmp_2_17=ponderation*reg191; T tmp_4_13=ponderation*reg230; T tmp_4_12=ponderation*reg376; T tmp_3_3=ponderation*reg537; T tmp_4_11=ponderation*reg377;
    T tmp_4_10=ponderation*reg213; T tmp_3_4=ponderation*reg530; T tmp_0_9=ponderation*reg203; T tmp_10_17=ponderation*reg419; T tmp_10_16=ponderation*reg298;
    T tmp_0_10=ponderation*reg232; T tmp_10_15=ponderation*reg333; T tmp_10_14=ponderation*reg371; T tmp_10_13=ponderation*reg379; T tmp_0_11=ponderation*reg280;
    T tmp_10_12=ponderation*reg361; T tmp_10_11=ponderation*reg364; T tmp_0_12=ponderation*reg521; T tmp_10_10=ponderation*reg185; T tmp_9_17=ponderation*reg365;
    T tmp_0_13=ponderation*reg518; T tmp_9_16=ponderation*reg310; T tmp_9_15=ponderation*reg194; T tmp_0_14=ponderation*reg515; T tmp_9_14=ponderation*reg441;
    T tmp_0_4=ponderation*reg464; T tmp_13_14=-reg292; T tmp_13_13=ponderation*reg432; T tmp_0_5=ponderation*reg504; T tmp_12_17=ponderation*reg399;
    T tmp_12_16=ponderation*reg395; T tmp_12_15=ponderation*reg341; T tmp_12_13=-reg331; T tmp_0_6=ponderation*reg502; T tmp_12_12=ponderation*reg413;
    T tmp_0_7=ponderation*reg498; T tmp_11_17=ponderation*reg273; T tmp_11_16=ponderation*reg407; T tmp_11_15=ponderation*reg168; T tmp_11_14=ponderation*reg306;
    T tmp_11_13=ponderation*reg410; T tmp_0_8=-reg340; T tmp_11_12=ponderation*reg400; T tmp_11_11=ponderation*reg300; T tmp_8_10=ponderation*reg411;
    T tmp_8_9=-reg404; T tmp_1_9=ponderation*reg465; T tmp_8_8=ponderation*reg201; T tmp_1_10=ponderation*reg253; T tmp_7_17=ponderation*reg386;
    T tmp_7_16=ponderation*reg240; T tmp_1_11=ponderation*reg461; T tmp_7_15=ponderation*reg390; T tmp_7_14=ponderation*reg393; T tmp_7_13=-reg396;
    T tmp_1_12=ponderation*reg460; T tmp_7_12=ponderation*reg402; T tmp_1_13=ponderation*reg458; T tmp_7_11=ponderation*reg403; T tmp_7_10=ponderation*reg250;
    T tmp_1_14=ponderation*reg456; T tmp_7_9=ponderation*reg433; T tmp_7_8=-reg302; T tmp_9_13=ponderation*reg439; T tmp_0_15=ponderation*reg514;
    T tmp_9_12=ponderation*reg198; T tmp_9_11=ponderation*reg397; T tmp_0_16=ponderation*reg112; T tmp_9_10=ponderation*reg234; T tmp_9_9=ponderation*reg162;
    T tmp_0_17=ponderation*reg512; T tmp_8_17=ponderation*reg165; T tmp_1_1=ponderation*reg160; T tmp_8_16=ponderation*reg405; T tmp_8_15=-reg384;
    T tmp_1_2=ponderation*reg509; T tmp_8_14=-reg381; T tmp_8_13=ponderation*reg409; T tmp_1_3=ponderation*reg467; T tmp_8_12=-reg363;
    T tmp_1_8=ponderation*reg526; T tmp_8_11=ponderation*reg229;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[0]+2) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_3;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_4;
    matrix(indices[0]+0,indices[1]+2) += tmp_0_5;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_7;
    matrix(indices[0]+0,indices[2]+2) += tmp_0_8;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_9;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_10;
    matrix(indices[0]+0,indices[3]+2) += tmp_0_11;
    matrix(indices[0]+0,indices[4]+0) += tmp_0_12;
    matrix(indices[0]+0,indices[4]+1) += tmp_0_13;
    matrix(indices[0]+0,indices[4]+2) += tmp_0_14;
    matrix(indices[0]+0,indices[5]+0) += tmp_0_15;
    matrix(indices[0]+0,indices[5]+1) += tmp_0_16;
    matrix(indices[0]+0,indices[5]+2) += tmp_0_17;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[0]+2) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_3;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_4;
    matrix(indices[0]+1,indices[1]+2) += tmp_1_5;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_7;
    matrix(indices[0]+1,indices[2]+2) += tmp_1_8;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_9;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_10;
    matrix(indices[0]+1,indices[3]+2) += tmp_1_11;
    matrix(indices[0]+1,indices[4]+0) += tmp_1_12;
    matrix(indices[0]+1,indices[4]+1) += tmp_1_13;
    matrix(indices[0]+1,indices[4]+2) += tmp_1_14;
    matrix(indices[0]+1,indices[5]+0) += tmp_1_15;
    matrix(indices[0]+1,indices[5]+1) += tmp_1_16;
    matrix(indices[0]+1,indices[5]+2) += tmp_1_17;
    matrix(indices[0]+2,indices[0]+2) += tmp_2_2;
    matrix(indices[0]+2,indices[1]+0) += tmp_2_3;
    matrix(indices[0]+2,indices[1]+1) += tmp_2_4;
    matrix(indices[0]+2,indices[1]+2) += tmp_2_5;
    matrix(indices[0]+2,indices[2]+0) += tmp_2_6;
    matrix(indices[0]+2,indices[2]+1) += tmp_2_7;
    matrix(indices[0]+2,indices[2]+2) += tmp_2_8;
    matrix(indices[0]+2,indices[3]+0) += tmp_2_9;
    matrix(indices[0]+2,indices[3]+1) += tmp_2_10;
    matrix(indices[0]+2,indices[3]+2) += tmp_2_11;
    matrix(indices[0]+2,indices[4]+0) += tmp_2_12;
    matrix(indices[0]+2,indices[4]+1) += tmp_2_13;
    matrix(indices[0]+2,indices[4]+2) += tmp_2_14;
    matrix(indices[0]+2,indices[5]+0) += tmp_2_15;
    matrix(indices[0]+2,indices[5]+1) += tmp_2_16;
    matrix(indices[0]+2,indices[5]+2) += tmp_2_17;
    matrix(indices[1]+0,indices[1]+0) += tmp_3_3;
    matrix(indices[1]+0,indices[1]+1) += tmp_3_4;
    matrix(indices[1]+0,indices[1]+2) += tmp_3_5;
    matrix(indices[1]+0,indices[2]+0) += tmp_3_6;
    matrix(indices[1]+0,indices[2]+1) += tmp_3_7;
    matrix(indices[1]+0,indices[2]+2) += tmp_3_8;
    matrix(indices[1]+0,indices[3]+0) += tmp_3_9;
    matrix(indices[1]+0,indices[3]+1) += tmp_3_10;
    matrix(indices[1]+0,indices[3]+2) += tmp_3_11;
    matrix(indices[1]+0,indices[4]+0) += tmp_3_12;
    matrix(indices[1]+0,indices[4]+1) += tmp_3_13;
    matrix(indices[1]+0,indices[4]+2) += tmp_3_14;
    matrix(indices[1]+0,indices[5]+0) += tmp_3_15;
    matrix(indices[1]+0,indices[5]+1) += tmp_3_16;
    matrix(indices[1]+0,indices[5]+2) += tmp_3_17;
    matrix(indices[1]+1,indices[1]+1) += tmp_4_4;
    matrix(indices[1]+1,indices[1]+2) += tmp_4_5;
    matrix(indices[1]+1,indices[2]+0) += tmp_4_6;
    matrix(indices[1]+1,indices[2]+1) += tmp_4_7;
    matrix(indices[1]+1,indices[2]+2) += tmp_4_8;
    matrix(indices[1]+1,indices[3]+0) += tmp_4_9;
    matrix(indices[1]+1,indices[3]+1) += tmp_4_10;
    matrix(indices[1]+1,indices[3]+2) += tmp_4_11;
    matrix(indices[1]+1,indices[4]+0) += tmp_4_12;
    matrix(indices[1]+1,indices[4]+1) += tmp_4_13;
    matrix(indices[1]+1,indices[4]+2) += tmp_4_14;
    matrix(indices[1]+1,indices[5]+0) += tmp_4_15;
    matrix(indices[1]+1,indices[5]+1) += tmp_4_16;
    matrix(indices[1]+1,indices[5]+2) += tmp_4_17;
    matrix(indices[1]+2,indices[1]+2) += tmp_5_5;
    matrix(indices[1]+2,indices[2]+0) += tmp_5_6;
    matrix(indices[1]+2,indices[2]+1) += tmp_5_7;
    matrix(indices[1]+2,indices[2]+2) += tmp_5_8;
    matrix(indices[1]+2,indices[3]+0) += tmp_5_9;
    matrix(indices[1]+2,indices[3]+1) += tmp_5_10;
    matrix(indices[1]+2,indices[3]+2) += tmp_5_11;
    matrix(indices[1]+2,indices[4]+0) += tmp_5_12;
    matrix(indices[1]+2,indices[4]+1) += tmp_5_13;
    matrix(indices[1]+2,indices[4]+2) += tmp_5_14;
    matrix(indices[1]+2,indices[5]+0) += tmp_5_15;
    matrix(indices[1]+2,indices[5]+1) += tmp_5_16;
    matrix(indices[1]+2,indices[5]+2) += tmp_5_17;
    matrix(indices[2]+0,indices[2]+0) += tmp_6_6;
    matrix(indices[2]+0,indices[2]+1) += tmp_6_7;
    matrix(indices[2]+0,indices[2]+2) += tmp_6_8;
    matrix(indices[2]+0,indices[3]+0) += tmp_6_9;
    matrix(indices[2]+0,indices[3]+1) += tmp_6_10;
    matrix(indices[2]+0,indices[3]+2) += tmp_6_11;
    matrix(indices[2]+0,indices[4]+0) += tmp_6_12;
    matrix(indices[2]+0,indices[4]+1) += tmp_6_13;
    matrix(indices[2]+0,indices[4]+2) += tmp_6_14;
    matrix(indices[2]+0,indices[5]+0) += tmp_6_15;
    matrix(indices[2]+0,indices[5]+1) += tmp_6_16;
    matrix(indices[2]+0,indices[5]+2) += tmp_6_17;
    matrix(indices[2]+1,indices[2]+1) += tmp_7_7;
    matrix(indices[2]+1,indices[2]+2) += tmp_7_8;
    matrix(indices[2]+1,indices[3]+0) += tmp_7_9;
    matrix(indices[2]+1,indices[3]+1) += tmp_7_10;
    matrix(indices[2]+1,indices[3]+2) += tmp_7_11;
    matrix(indices[2]+1,indices[4]+0) += tmp_7_12;
    matrix(indices[2]+1,indices[4]+1) += tmp_7_13;
    matrix(indices[2]+1,indices[4]+2) += tmp_7_14;
    matrix(indices[2]+1,indices[5]+0) += tmp_7_15;
    matrix(indices[2]+1,indices[5]+1) += tmp_7_16;
    matrix(indices[2]+1,indices[5]+2) += tmp_7_17;
    matrix(indices[2]+2,indices[2]+2) += tmp_8_8;
    matrix(indices[2]+2,indices[3]+0) += tmp_8_9;
    matrix(indices[2]+2,indices[3]+1) += tmp_8_10;
    matrix(indices[2]+2,indices[3]+2) += tmp_8_11;
    matrix(indices[2]+2,indices[4]+0) += tmp_8_12;
    matrix(indices[2]+2,indices[4]+1) += tmp_8_13;
    matrix(indices[2]+2,indices[4]+2) += tmp_8_14;
    matrix(indices[2]+2,indices[5]+0) += tmp_8_15;
    matrix(indices[2]+2,indices[5]+1) += tmp_8_16;
    matrix(indices[2]+2,indices[5]+2) += tmp_8_17;
    matrix(indices[3]+0,indices[3]+0) += tmp_9_9;
    matrix(indices[3]+0,indices[3]+1) += tmp_9_10;
    matrix(indices[3]+0,indices[3]+2) += tmp_9_11;
    matrix(indices[3]+0,indices[4]+0) += tmp_9_12;
    matrix(indices[3]+0,indices[4]+1) += tmp_9_13;
    matrix(indices[3]+0,indices[4]+2) += tmp_9_14;
    matrix(indices[3]+0,indices[5]+0) += tmp_9_15;
    matrix(indices[3]+0,indices[5]+1) += tmp_9_16;
    matrix(indices[3]+0,indices[5]+2) += tmp_9_17;
    matrix(indices[3]+1,indices[3]+1) += tmp_10_10;
    matrix(indices[3]+1,indices[3]+2) += tmp_10_11;
    matrix(indices[3]+1,indices[4]+0) += tmp_10_12;
    matrix(indices[3]+1,indices[4]+1) += tmp_10_13;
    matrix(indices[3]+1,indices[4]+2) += tmp_10_14;
    matrix(indices[3]+1,indices[5]+0) += tmp_10_15;
    matrix(indices[3]+1,indices[5]+1) += tmp_10_16;
    matrix(indices[3]+1,indices[5]+2) += tmp_10_17;
    matrix(indices[3]+2,indices[3]+2) += tmp_11_11;
    matrix(indices[3]+2,indices[4]+0) += tmp_11_12;
    matrix(indices[3]+2,indices[4]+1) += tmp_11_13;
    matrix(indices[3]+2,indices[4]+2) += tmp_11_14;
    matrix(indices[3]+2,indices[5]+0) += tmp_11_15;
    matrix(indices[3]+2,indices[5]+1) += tmp_11_16;
    matrix(indices[3]+2,indices[5]+2) += tmp_11_17;
    matrix(indices[4]+0,indices[4]+0) += tmp_12_12;
    matrix(indices[4]+0,indices[4]+1) += tmp_12_13;
    matrix(indices[4]+0,indices[4]+2) += tmp_12_14;
    matrix(indices[4]+0,indices[5]+0) += tmp_12_15;
    matrix(indices[4]+0,indices[5]+1) += tmp_12_16;
    matrix(indices[4]+0,indices[5]+2) += tmp_12_17;
    matrix(indices[4]+1,indices[4]+1) += tmp_13_13;
    matrix(indices[4]+1,indices[4]+2) += tmp_13_14;
    matrix(indices[4]+1,indices[5]+0) += tmp_13_15;
    matrix(indices[4]+1,indices[5]+1) += tmp_13_16;
    matrix(indices[4]+1,indices[5]+2) += tmp_13_17;
    matrix(indices[4]+2,indices[4]+2) += tmp_14_14;
    matrix(indices[4]+2,indices[5]+0) += tmp_14_15;
    matrix(indices[4]+2,indices[5]+1) += tmp_14_16;
    matrix(indices[4]+2,indices[5]+2) += tmp_14_17;
    matrix(indices[5]+0,indices[5]+0) += tmp_15_15;
    matrix(indices[5]+0,indices[5]+1) += tmp_15_16;
    matrix(indices[5]+0,indices[5]+2) += tmp_15_17;
    matrix(indices[5]+1,indices[5]+1) += tmp_16_16;
    matrix(indices[5]+1,indices[5]+2) += tmp_16_17;
    matrix(indices[5]+2,indices[5]+2) += tmp_17_17;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1-var_inter[0]; reg0=reg0-var_inter[1]; T reg1=var_inter[0]*elem.pos(1)[1]; T reg2=reg0*elem.pos(0)[1]; T reg3=var_inter[0]*elem.pos(1)[2];
    T reg4=reg0*elem.pos(0)[2]; T reg5=reg3+reg4; T reg6=var_inter[1]*elem.pos(2)[2]; T reg7=reg1+reg2; T reg8=1-var_inter[2];
    T reg9=var_inter[1]*elem.pos(2)[1]; T reg10=reg9+reg7; T reg11=reg8*elem.pos(2)[2]; T reg12=reg8*elem.pos(1)[1]; T reg13=reg8*elem.pos(0)[1];
    T reg14=elem.pos(3)[1]*reg0; T reg15=reg8*elem.pos(2)[1]; T reg16=reg8*elem.pos(1)[2]; T reg17=reg8*elem.pos(0)[2]; T reg18=reg5+reg6;
    T reg19=elem.pos(3)[2]*reg0; reg12=reg12-reg13; T reg20=var_inter[2]*elem.pos(3)[1]; reg11=reg11-reg17; reg15=reg15-reg13;
    T reg21=var_inter[0]*elem.pos(1)[0]; T reg22=elem.pos(0)[0]*reg0; reg19=reg19-reg18; T reg23=var_inter[0]*elem.pos(4)[2]; reg14=reg14-reg10;
    T reg24=var_inter[0]*elem.pos(4)[1]; reg16=reg16-reg17; T reg25=var_inter[2]*elem.pos(3)[2]; T reg26=var_inter[2]*elem.pos(5)[1]; T reg27=reg8*elem.pos(2)[0];
    reg15=reg15-reg20; T reg28=var_inter[2]*elem.pos(5)[2]; reg11=reg11-reg25; T reg29=reg21+reg22; T reg30=var_inter[1]*elem.pos(2)[0];
    T reg31=var_inter[1]*elem.pos(5)[1]; T reg32=1+(*f.m).poisson_ratio; reg24=reg14+reg24; reg14=var_inter[1]*elem.pos(5)[2]; reg23=reg19+reg23;
    reg19=var_inter[2]*elem.pos(4)[2]; reg16=reg16-reg25; T reg33=var_inter[2]*elem.pos(4)[1]; reg12=reg12-reg20; T reg34=reg8*elem.pos(1)[0];
    T reg35=reg8*elem.pos(0)[0]; reg32=reg32/(*f.m).elastic_modulus; reg23=reg14+reg23; reg24=reg31+reg24; reg14=elem.pos(3)[0]*reg0;
    reg31=reg29+reg30; reg11=reg28+reg11; reg34=reg34-reg35; reg28=var_inter[2]*elem.pos(3)[0]; reg27=reg27-reg35;
    reg33=reg12+reg33; reg15=reg26+reg15; reg19=reg16+reg19; reg12=reg33*reg23; reg16=reg15*reg23;
    reg26=reg11*reg24; reg27=reg27-reg28; T reg36=reg19*reg24; T reg37=var_inter[2]*elem.pos(5)[0]; T reg38=pow(reg32,2);
    reg34=reg34-reg28; T reg39=var_inter[2]*elem.pos(4)[0]; reg14=reg14-reg31; T reg40=var_inter[0]*elem.pos(4)[0]; T reg41=var_inter[1]*elem.pos(5)[0];
    reg27=reg37+reg27; reg39=reg34+reg39; reg26=reg16-reg26; reg40=reg14+reg40; reg32=reg32*reg38;
    reg36=reg12-reg36; reg12=reg33*reg11; reg14=reg19*reg15; reg16=1.0/(*f.m).elastic_modulus; reg34=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    reg37=reg27*reg36; T reg42=reg39*reg26; reg14=reg12-reg14; reg12=reg34*reg38; T reg43=reg16*reg32;
    reg32=reg34*reg32; reg38=reg16*reg38; reg40=reg41+reg40; reg41=reg34*reg43; T reg44=reg34*reg32;
    reg43=reg16*reg43; T reg45=reg19*reg27; T reg46=reg39*reg24; T reg47=reg39*reg11; T reg48=reg15*reg40;
    T reg49=reg39*reg23; reg24=reg27*reg24; reg11=reg11*reg40; T reg50=reg16*reg38; reg23=reg27*reg23;
    T reg51=reg34*reg12; T reg52=reg40*reg14; reg38=reg34*reg38; reg37=reg42-reg37; reg42=reg33*reg40;
    reg40=reg19*reg40; reg38=reg51+reg38; reg43=reg43-reg44; reg40=reg49-reg40; reg27=reg33*reg27;
    reg41=reg44+reg41; reg32=reg16*reg32; reg45=reg47-reg45; reg15=reg39*reg15; reg52=reg37+reg52;
    reg12=reg16*reg12; reg50=reg50-reg51; reg48=reg24-reg48; reg11=reg23-reg11; reg42=reg46-reg42;
    reg32=reg44+reg32; reg38=reg34*reg38; reg50=reg16*reg50; reg19=reg51+reg12; reg40=reg40/reg52;
    reg42=reg42/reg52; reg14=reg14/reg52; reg45=reg45/reg52; reg27=reg15-reg27; reg26=reg26/reg52;
    reg11=reg11/reg52; reg48=reg48/reg52; reg36=reg36/reg52; reg15=reg34*reg41; reg16=reg16*reg43;
    reg23=reg8*reg40; reg24=reg8*reg11; reg33=reg8*reg48; reg37=var_inter[1]*reg45; reg39=reg8*reg42;
    reg44=var_inter[0]*reg14; reg46=var_inter[0]*reg45; reg47=var_inter[1]*reg14; reg49=reg8*reg36; T reg53=var_inter[2]*reg42;
    T reg54=var_inter[2]*reg48; T reg55=var_inter[2]*reg11; T reg56=var_inter[2]*reg40; T reg57=var_inter[2]*reg36; T reg58=var_inter[2]*reg26;
    reg27=reg27/reg52; reg38=reg50-reg38; reg15=reg16-reg15; reg16=reg34*reg32; reg50=reg8*reg26;
    reg19=reg34*reg19; reg34=reg49+reg47; T reg59=reg37+reg23; T reg60=var_inter[1]*reg27; T reg61=reg49-reg50;
    T reg62=reg0*reg14; T reg63=reg24-reg23; T reg64=reg0*reg45; T reg65=reg0*reg27; T reg66=reg39-reg33;
    reg19=reg38-reg19; reg38=reg53-reg54; T reg67=reg58+reg44; T reg68=var_inter[0]*reg27; T reg69=reg55+reg46;
    T reg70=reg57-reg58; reg16=reg15-reg16; reg15=reg55-reg56; T reg71=reg46-reg24; T reg72=reg33-reg68;
    reg63=reg63+reg64; reg70=reg70+reg62; reg15=reg15-reg64; reg61=reg61-reg62; T reg73=reg68+reg54;
    T reg74=0.5*reg69; T reg75=reg60+reg39; T reg76=0.5*reg59; reg38=reg38+reg65; T reg77=0.5*reg34;
    T reg78=reg60-reg53; T reg79=reg56-reg37; T reg80=reg47-reg57; T reg81=0.5*reg67; T reg82=reg50-reg44;
    reg66=reg66-reg65; reg19=reg19/reg16; T reg83=reg19*reg81; T reg84=reg19*reg74; T reg85=0.5*reg15;
    reg43=reg43/reg16; T reg86=0.5*reg73; T reg87=0.5*reg82; T reg88=0.5*reg72; T reg89=0.5*reg71;
    T reg90=0.5*reg70; T reg91=0.5*reg79; T reg92=0.5*reg38; T reg93=0.5*reg78; T reg94=0.5*reg63;
    T reg95=0.5*reg80; T reg96=reg19*reg77; T reg97=0.5*reg75; T reg98=0.5*reg66; T reg99=0.5*reg61;
    T reg100=reg19*reg76; T reg101=reg19*reg93; T reg102=reg19*reg91; reg83=2*reg83; T reg103=reg19*reg95;
    T reg104=reg19*reg86; T reg105=2*reg84; T reg106=reg43*reg34; T reg107=reg43*reg75; T reg108=reg43*reg67;
    T reg109=reg19*reg90; T reg110=reg19*reg92; reg100=2*reg100; T reg111=reg43*reg69; T reg112=reg19*reg97;
    T reg113=reg19*reg85; T reg114=2*reg96; T reg115=reg19*reg94; T reg116=reg43*reg73; T reg117=reg19*reg98;
    T reg118=reg19*reg88; reg32=reg32/reg16; reg16=reg41/reg16; reg41=reg19*reg89; T reg119=reg19*reg99;
    T reg120=reg43*reg59; T reg121=reg19*reg87; T reg122=reg67*reg106; T reg123=reg32*reg66; T reg124=reg74*reg100;
    T reg125=reg32*reg78; T reg126=reg75*reg116; reg110=2*reg110; T reg127=reg120*reg69; T reg128=reg43*reg82;
    T reg129=reg81*reg114; reg41=2*reg41; T reg130=reg43*reg78; T reg131=reg32*reg73; T reg132=reg43*reg72;
    reg113=2*reg113; reg119=2*reg119; T reg133=reg32*reg72; T reg134=reg34*reg108; T reg135=reg43*reg70;
    T reg136=reg43*reg63; T reg137=reg76*reg105; T reg138=2*reg112; T reg139=reg32*reg59; T reg140=reg32*reg75;
    T reg141=reg16*reg61; T reg142=reg43*reg80; T reg143=reg16*reg67; T reg144=reg43*reg61; reg115=2*reg115;
    reg101=2*reg101; reg104=2*reg104; T reg145=reg16*reg80; T reg146=reg43*reg38; T reg147=reg77*reg83;
    T reg148=reg43*reg15; T reg149=reg43*reg79; T reg150=reg59*reg111; T reg151=reg32*reg38; T reg152=reg73*reg107;
    T reg153=reg32*reg69; T reg154=reg16*reg59; reg118=2*reg118; T reg155=reg16*reg82; reg102=2*reg102;
    T reg156=reg43*reg66; reg117=2*reg117; T reg157=reg16*reg69; reg121=2*reg121; T reg158=reg43*reg71;
    T reg159=reg16*reg34; reg109=2*reg109; reg103=2*reg103; T reg160=reg16*reg70; T reg161=reg75*reg160;
    T reg162=reg75*reg107; T reg163=reg77*reg110; T reg164=reg75*reg146; T reg165=reg159*reg73; T reg166=reg59*reg149;
    T reg167=reg77*reg103; T reg168=reg128*reg82; T reg169=reg73*reg132; T reg170=reg75*reg141; T reg171=reg77*reg117;
    T reg172=reg75*reg156; T reg173=reg75*reg155; T reg174=reg77*reg118; T reg175=reg89*reg115; T reg176=reg66*reg130;
    T reg177=reg79*reg148; T reg178=reg75*reg132; T reg179=reg32*reg79; T reg180=reg95*reg109; T reg181=reg73*reg156;
    T reg182=reg76*reg138; T reg183=reg75*reg139; T reg184=reg66*reg116; T reg185=reg95*reg103; T reg186=reg81*reg83;
    T reg187=reg128*reg70; T reg188=reg41*reg85; T reg189=reg69*reg143; T reg190=reg81*reg105; T reg191=reg70*reg106;
    T reg192=reg85*reg100; T reg193=reg32*reg63; T reg194=reg70*reg140; T reg195=reg92*reg114; T reg196=reg70*reg135;
    T reg197=reg85*reg113; T reg198=reg99*reg103; T reg199=reg63*reg149; T reg200=reg69*reg148; T reg201=reg81*reg109;
    T reg202=reg78*reg156; T reg203=reg70*reg108; T reg204=reg85*reg105; T reg205=reg99*reg83; T reg206=reg63*reg111;
    T reg207=reg75*reg143; T reg208=reg77*reg104; T reg209=reg66*reg146; T reg210=reg79*reg111; T reg211=reg69*reg149;
    T reg212=reg32*reg15; T reg213=reg95*reg83; T reg214=reg81*reg103; reg126=reg147+reg126; T reg215=reg75*reg145;
    T reg216=reg77*reg101; T reg217=reg66*reg107; T reg218=reg75*reg130; T reg219=reg86*reg105; T reg220=reg70*reg144;
    T reg221=reg99*reg138; T reg222=reg159*reg66; T reg223=reg85*reg115; T reg224=reg69*reg131; T reg225=reg66*reg132;
    T reg226=reg69*reg111; T reg227=reg32*reg71; T reg228=reg79*reg149; T reg229=reg123*reg34; T reg230=reg119*reg97;
    T reg231=reg87*reg83; T reg232=reg71*reg111; T reg233=reg41*reg76; T reg234=reg128*reg34; T reg235=reg34*reg133;
    T reg236=reg121*reg97; T reg237=reg87*reg109; T reg238=reg71*reg148; T reg239=reg76*reg100; T reg240=reg34*reg106;
    T reg241=reg80*reg142; T reg242=reg91*reg102; T reg243=reg76*reg114; T reg244=reg34*reg154; T reg245=reg91*reg115;
    T reg246=reg80*reg144; T reg247=reg87*reg114; T reg248=reg120*reg71; T reg249=reg76*reg113; T reg250=reg34*reg135;
    T reg251=reg73*reg130; T reg252=reg93*reg114; T reg253=reg80*reg140; T reg254=reg72*reg107; T reg255=reg87*reg138;
    T reg256=reg159*reg72; T reg257=reg72*reg132; T reg258=reg80*reg135; T reg259=reg91*reg113; T reg260=reg72*reg146;
    T reg261=reg91*reg100; T reg262=reg80*reg106; T reg263=reg67*reg142; T reg264=reg72*reg156; T reg265=reg74*reg102;
    T reg266=reg72*reg116; T reg267=reg72*reg130; T reg268=reg76*reg115; T reg269=reg87*reg103; T reg270=reg71*reg149;
    T reg271=reg80*reg108; T reg272=reg91*reg105; T reg273=reg34*reg144; T reg274=reg91*reg41; T reg275=reg80*reg128;
    T reg276=reg73*reg146; T reg277=reg159*reg59; T reg278=reg77*reg100; T reg279=reg82*reg108; T reg280=reg89*reg105;
    T reg281=reg79*reg158; T reg282=reg120*reg59; T reg283=reg77*reg114; T reg284=reg121*reg95; T reg285=reg59*reg140;
    T reg286=reg97*reg100; T reg287=reg82*reg135; T reg288=reg89*reg113; T reg289=reg59*reg148; T reg290=reg77*reg109;
    T reg291=reg129+reg152; T reg292=reg88*reg114; T reg293=reg82*reg140; reg147=reg150+reg147; T reg294=reg82*reg106;
    T reg295=reg89*reg100; T reg296=reg81*reg138; T reg297=reg120*reg79; T reg298=reg114*reg95; T reg299=reg34*reg151;
    T reg300=reg97*reg109; T reg301=reg121*reg87; T reg302=reg71*reg158; reg134=reg137+reg134; T reg303=reg97*reg104;
    T reg304=reg34*reg131; T reg305=reg97*reg83; T reg306=reg76*reg102; T reg307=reg34*reg142; T reg308=reg119*reg87;
    T reg309=reg73*reg116; T reg310=reg79*reg136; T reg311=reg119*reg95; T reg312=reg34*reg125; T reg313=reg97*reg103;
    T reg314=reg74*reg104; T reg315=reg73*reg153; T reg316=reg59*reg136; T reg317=reg77*reg119; T reg318=reg82*reg142;
    T reg319=reg89*reg102; T reg320=reg59*reg158; T reg321=reg77*reg121; T reg322=reg38*reg130; T reg323=reg122+reg124;
    T reg324=reg86*reg138; T reg325=reg128*reg61; T reg326=reg90*reg114; T reg327=reg120*reg15; T reg328=reg41*reg94;
    T reg329=reg94*reg105; T reg330=reg16*reg71; T reg331=reg63*reg136; T reg332=reg119*reg99; T reg333=reg121*reg90;
    T reg334=reg15*reg158; T reg335=reg41*reg74; T reg336=reg61*reg108; reg128=reg128*reg67; T reg337=reg38*reg116;
    T reg338=reg138*reg95; T reg339=reg81*reg121; T reg340=reg159*reg78; T reg341=reg16*reg79; T reg342=reg67*reg140;
    T reg343=reg61*reg144; reg149=reg15*reg149; T reg344=reg61*reg106; T reg345=reg16*reg63; T reg346=reg66*reg156;
    T reg347=reg94*reg100; T reg348=reg61*reg142; reg120=reg120*reg63; T reg349=reg15*reg111; T reg350=reg94*reg102;
    T reg351=reg90*reg109; T reg352=reg63*reg158; T reg353=reg38*reg107; reg130=reg78*reg130; T reg354=reg41*reg89;
    T reg355=reg15*reg148; T reg356=reg99*reg114; T reg357=reg90*reg83; T reg358=reg90*reg138; T reg359=reg159*reg38;
    T reg360=reg74*reg115; T reg361=reg78*reg107; T reg362=reg81*reg119; T reg363=reg67*reg144; T reg364=reg69*reg136;
    T reg365=reg38*reg132; reg142=reg70*reg142; T reg366=reg85*reg102; T reg367=reg74*reg105; reg132=reg78*reg132;
    T reg368=reg94*reg113; T reg369=reg74*reg83; T reg370=reg67*reg135; T reg371=reg98*reg114; T reg372=reg86*reg114;
    T reg373=reg61*reg140; reg127=reg129+reg127; T reg374=reg16*reg15; T reg375=reg78*reg146; reg146=reg38*reg146;
    T reg376=reg67*reg157; T reg377=reg71*reg136; reg144=reg82*reg144; T reg378=reg99*reg109; T reg379=reg90*reg103;
    reg158=reg69*reg158; reg116=reg78*reg116; T reg380=reg94*reg115; T reg381=reg119*reg90; reg156=reg38*reg156;
    T reg382=reg121*reg99; reg108=reg67*reg108; T reg383=reg74*reg113; reg148=reg63*reg148; reg135=reg61*reg135;
    reg136=reg15*reg136; T reg384=reg74*reg101; T reg385=reg59*reg155; reg251=reg214+reg251; T reg386=reg86*reg117;
    T reg387=reg41*reg97; T reg388=reg76*reg109; T reg389=reg77*reg41; T reg390=reg34*reg374; reg276=reg201+reg276;
    reg156=reg381+reg156; reg300=reg299+reg300; T reg391=reg73*reg179; reg320=reg320-reg321; reg108=reg108+reg367;
    T reg392=reg74*reg110; T reg393=reg359+reg358; T reg394=reg59*reg133; reg314=reg315+reg314; reg313=reg312+reg313;
    T reg395=reg67*reg133; T reg396=reg59*reg141; T reg397=reg77*reg115; T reg398=reg34*reg341; T reg399=reg121*reg86;
    T reg400=reg118*reg90; T reg401=reg38*reg155; T reg402=reg123*reg67; reg316=reg316-reg317; T reg403=reg81*reg104;
    T reg404=reg76*reg103; T reg405=reg38*reg227; reg309=reg186+reg309; T reg406=reg123*reg59; T reg407=reg97*reg101;
    reg307=reg306-reg307; T reg408=reg118*reg85; T reg409=reg86*reg109; T reg410=reg119*reg74; reg305=reg304+reg305;
    T reg411=reg345*reg67; T reg412=reg67*reg151; T reg413=reg73*reg145; T reg414=reg81*reg101; T reg415=reg34*reg157;
    T reg416=reg76*reg83; T reg417=reg134+reg303; T reg418=reg97*reg115; T reg419=reg73*reg143; reg365=reg333+reg365;
    T reg420=reg86*reg101; T reg421=reg72*reg153; T reg422=reg89*reg104; T reg423=reg93*reg121; T reg424=reg38*reg143;
    reg266=reg231+reg266; T reg425=reg80*reg133; T reg426=reg72*reg145; T reg427=reg87*reg101; T reg428=reg91*reg121;
    T reg429=reg72*reg179; T reg430=reg89*reg101; T reg431=reg80*reg330; T reg432=reg90*reg104; reg267=reg269+reg267;
    T reg433=reg93*reg118; T reg434=reg38*reg153; T reg435=reg85*reg104; T reg436=reg74*reg114; T reg437=reg89*reg138;
    T reg438=reg38*reg160; T reg439=reg90*reg110; T reg440=reg247+reg254; T reg441=reg91*reg114; T reg442=reg72*reg160;
    T reg443=reg87*reg110; T reg444=reg38*reg212; T reg445=reg85*reg110; T reg446=reg80*reg154; T reg447=reg72*reg212;
    T reg448=reg89*reg110; T reg449=reg93*reg138; reg260=reg237+reg260; T reg450=reg261-reg262; T reg451=reg72*reg143;
    T reg452=reg86*reg83; reg146=reg351+reg146; T reg453=reg342+reg372; reg263=reg263-reg265; T reg454=reg80*reg345;
    reg236=reg235+reg236; reg322=reg379+reg322; reg363=reg363-reg360; T reg455=reg239+reg240; T reg456=reg97*reg138;
    T reg457=reg93*reg117; reg370=reg370-reg383; T reg458=reg86*reg110; reg246=reg246+reg245; T reg459=reg326+reg353;
    reg244=reg243+reg244; T reg460=reg85*reg138; T reg461=reg38*reg139; T reg462=reg34*reg140; T reg463=reg67*reg374;
    T reg464=reg97*reg114; T reg465=reg74*reg109; reg250=reg249-reg250; T reg466=reg97*reg110; reg337=reg357+reg337;
    reg273=reg268-reg273; T reg467=reg97*reg117; reg275=reg275+reg274; T reg468=reg67*reg154; T reg469=reg119*reg76;
    T reg470=reg345*reg34; T reg471=reg38*reg145; T reg472=reg90*reg101; reg230=reg229+reg230; T reg473=reg93*reg119;
    T reg474=reg323+reg324; T reg475=reg38*reg179; T reg476=reg85*reg101; reg234=reg233-reg234; T reg477=reg118*reg97;
    T reg478=reg80*reg123; T reg479=reg81*reg100; T reg480=reg159*reg69; T reg481=reg91*reg119; T reg482=reg41*reg92;
    reg218=reg167+reg218; T reg483=reg224+reg219; T reg484=reg15*reg133; T reg485=reg81*reg41; reg333=reg334+reg333;
    reg220=reg220+reg223; reg334=reg92*reg117; T reg486=reg69*reg155; T reg487=reg119*reg85; reg186=reg186+reg226;
    T reg488=reg123*reg70; T reg489=reg119*reg92; T reg490=reg41*reg90; T reg491=reg15*reg155; T reg492=reg92*reg115;
    T reg493=reg123*reg15; reg187=reg187+reg188; T reg494=reg118*reg92; T reg495=reg330*reg70; T reg496=reg90*reg113;
    T reg497=reg15*reg160; reg208=reg207+reg208; reg211=reg214-reg211; reg207=reg76*reg104; reg214=reg75*reg153;
    T reg498=reg121*reg74; T reg499=reg92*reg100; T reg500=reg15*reg140; reg364=reg362-reg364; reg327=reg327-reg326;
    reg126=reg137+reg126; T reg501=reg123*reg69; reg216=reg215+reg216; reg215=reg90*reg100; T reg502=reg159*reg15;
    T reg503=reg69*reg145; T reg504=reg81*reg102; T reg505=reg76*reg101; T reg506=reg75*reg179; T reg507=reg86*reg115;
    T reg508=reg92*reg110; T reg509=reg85*reg103; T reg510=reg70*reg374; T reg511=reg85*reg109; T reg512=reg70*reg341;
    T reg513=reg70*reg151; T reg514=reg92*reg109; T reg515=reg69*reg160; T reg516=reg81*reg113; T reg517=reg324+reg127;
    T reg518=reg92*reg101; reg203=reg203-reg204; T reg519=reg92*reg104; reg142=reg142+reg366; T reg520=reg70*reg157;
    T reg521=reg85*reg83; T reg522=reg86*reg100; T reg523=reg69*reg140; T reg524=reg70*reg131; T reg525=reg92*reg83;
    T reg526=reg121*reg85; reg189=reg190+reg189; reg158=reg339-reg158; T reg527=reg70*reg133; T reg528=reg121*reg92;
    reg381=reg136+reg381; reg136=reg192-reg191; T reg529=reg92*reg138; T reg530=reg86*reg113; T reg531=reg70*reg154;
    T reg532=reg85*reg114; T reg533=reg69*reg151; T reg534=reg90*reg115; T reg535=reg194+reg195; reg200=reg201-reg200;
    reg201=reg15*reg141; T reg536=reg69*reg133; reg196=reg196+reg197; T reg537=reg41*reg86; T reg538=reg92*reg103;
    T reg539=reg70*reg125; T reg540=reg97*reg113; T reg541=reg59*reg143; T reg542=reg77*reg105; T reg543=reg74*reg138;
    T reg544=reg73*reg139; T reg545=reg67*reg131; T reg546=reg90*reg102; T reg547=reg15*reg145; reg303=reg303+reg147;
    T reg548=reg59*reg131; T reg549=reg97*reg105; T reg550=reg165+reg296; T reg551=reg59*reg145; T reg552=reg77*reg102;
    reg128=reg128-reg335; T reg553=reg67*reg341; reg167=reg166-reg167; reg169=reg339+reg169; reg166=reg59*reg125;
    reg339=reg118*reg86; T reg554=reg86*reg104; T reg555=reg85*reg117; reg278=reg277+reg278; T reg556=reg73*reg212;
    T reg557=reg38*reg193; T reg558=reg119*reg86; reg282=reg282+reg283; T reg559=reg81*reg110; T reg560=reg90*reg117;
    T reg561=reg38*reg141; reg286=reg285+reg286; T reg562=reg73*reg160; T reg563=reg59*reg160; T reg564=reg77*reg113;
    reg124=reg124+reg291; reg289=reg289-reg290; T reg565=reg92*reg102; T reg566=reg15*reg125; reg369=reg376+reg369;
    reg379=reg149+reg379; reg149=reg59*reg151; T reg567=reg86*reg103; T reg568=reg92*reg113; T reg569=reg15*reg151;
    T reg570=reg159*reg75; T reg571=reg77*reg138; T reg572=reg74*reg117; T reg573=reg81*reg115; reg351=reg355+reg351;
    reg183=reg182+reg183; reg355=reg73*reg193; T reg574=reg330*reg67; T reg575=reg283+reg162; T reg576=reg81*reg117;
    reg163=reg161+reg163; reg161=reg73*reg141; T reg577=reg76*reg110; T reg578=reg75*reg212; T reg579=reg86*reg102;
    reg164=reg290+reg164; reg290=reg69*reg125; T reg580=reg97*reg102; T reg581=reg92*reg105; T reg582=reg15*reg131;
    T reg583=reg74*reg103; reg171=reg170+reg171; reg170=reg118*reg74; T reg584=reg76*reg117; T reg585=reg75*reg193;
    T reg586=reg73*reg227; reg357=reg357-reg349; reg172=reg317+reg172; reg317=reg81*reg118; reg174=reg173+reg174;
    reg173=reg73*reg155; T reg587=reg118*reg76; T reg588=reg67*reg125; T reg589=reg75*reg227; T reg590=reg90*reg105;
    reg181=reg362+reg181; reg362=reg15*reg143; reg178=reg321+reg178; reg321=reg89*reg114; T reg591=reg82*reg154;
    T reg592=reg61*reg125; T reg593=reg94*reg103; T reg594=reg293+reg292; T reg595=reg61*reg341; T reg596=reg93*reg41;
    T reg597=reg79*reg133; reg287=reg288+reg287; T reg598=reg66*reg143; T reg599=reg298+reg361; T reg600=reg88*reg110;
    T reg601=reg98*reg101; reg348=reg348+reg350; T reg602=reg89*reg109; T reg603=reg82*reg374; reg281=reg281+reg284;
    T reg604=reg82*reg151; T reg605=reg123*reg63; reg168=reg354+reg168; T reg606=reg118*reg88; T reg607=reg340+reg338;
    reg331=reg331+reg332; reg297=reg297-reg298; T reg608=reg121*reg89; T reg609=reg330*reg82; T reg610=reg82*reg133;
    T reg611=reg121*reg88; T reg612=reg99*reg115; T reg613=reg63*reg141; T reg614=reg78*reg139; T reg615=reg91*reg138;
    T reg616=reg295-reg294; T reg617=reg88*reg138; T reg618=reg100*reg95; T reg619=reg159*reg79; T reg620=reg98*reg103;
    reg318=reg319+reg318; T reg621=reg91*reg110; T reg622=reg88*reg101; T reg623=reg79*reg123; T reg624=reg98*reg109;
    T reg625=reg61*reg151; T reg626=reg89*reg103; T reg627=reg82*reg341; reg310=reg310+reg311; T reg628=reg82*reg125;
    T reg629=reg88*reg103; T reg630=reg94*reg109; T reg631=reg61*reg374; reg375=reg180+reg375; T reg632=reg71*reg141;
    T reg633=reg87*reg115; T reg634=reg98*reg110; reg377=reg377+reg308; reg135=reg135+reg368; T reg635=reg88*reg109;
    T reg636=reg98*reg83; T reg637=reg61*reg131; reg279=reg279-reg280; T reg638=reg88*reg104; T reg639=reg78*reg160;
    T reg640=reg95*reg110; T reg641=reg41*reg95; T reg642=reg94*reg83; T reg643=reg89*reg83; T reg644=reg61*reg157;
    T reg645=reg82*reg157; T reg646=reg79*reg155; T reg647=reg82*reg131; T reg648=reg88*reg83; T reg649=reg98*reg104;
    reg336=reg336-reg329; T reg650=reg78*reg212; T reg651=reg93*reg115; T reg652=reg63*reg125; T reg653=reg94*reg138;
    T reg654=reg93*reg105; T reg655=reg79*reg131; T reg656=reg78*reg193; reg199=reg199+reg198; T reg657=reg91*reg117;
    T reg658=reg66*reg160; T reg659=reg99*reg110; T reg660=reg213-reg210; reg212=reg66*reg212; T reg661=reg99*reg102;
    T reg662=reg94*reg110; T reg663=reg63*reg145; T reg664=reg98*reg105; T reg665=reg95*reg105; T reg666=reg63*reg131;
    T reg667=reg121*reg76; reg202=reg311+reg202; reg311=reg79*reg125; T reg668=reg87*reg104; reg228=reg228+reg185;
    T reg669=reg66*reg155; T reg670=reg118*reg99; T reg671=reg93*reg102; T reg672=reg66*reg227; T reg673=reg118*reg94;
    T reg674=reg94*reg117; T reg675=reg66*reg193; T reg676=reg78*reg141; reg225=reg382+reg225; T reg677=reg95*reg102;
    T reg678=reg79*reg145; T reg679=reg222+reg221; T reg680=reg66*reg141; T reg681=reg95*reg117; T reg682=reg66*reg139;
    T reg683=reg98*reg102; reg148=reg148+reg378; T reg684=reg66*reg179; T reg685=reg94*reg101; T reg686=reg91*reg118;
    T reg687=reg99*reg113; T reg688=reg63*reg160; T reg689=reg95*reg113; reg144=reg144+reg175; T reg690=reg88*reg117;
    T reg691=reg79*reg160; reg132=reg284+reg132; reg284=reg119*reg89; T reg692=reg41*reg99; T reg693=reg93*reg100;
    T reg694=reg63*reg155; T reg695=reg123*reg82; T reg696=reg119*reg88; T reg697=reg98*reg115; T reg698=reg79*reg140;
    T reg699=reg330*reg34; T reg700=reg79*reg143; T reg701=reg205-reg206; T reg702=reg99*reg104; T reg703=reg78*reg155;
    T reg704=reg66*reg153; T reg705=reg94*reg104; T reg706=reg99*reg105; T reg707=reg63*reg143; T reg708=reg118*reg95;
    T reg709=reg93*reg113; T reg710=reg79*reg151; reg184=reg205+reg184; reg205=reg98*reg113; T reg711=reg66*reg145;
    T reg712=reg63*reg151; T reg713=reg99*reg101; T reg714=reg78*reg227; reg180=reg177+reg180; reg193=reg72*reg193;
    reg177=reg89*reg117; T reg715=reg123*reg61; T reg716=reg356+reg217; T reg717=reg99*reg117; reg382=reg352+reg382;
    reg209=reg378+reg209; reg352=reg88*reg113; reg378=reg71*reg151; T reg718=reg80*reg131; T reg719=reg93*reg83;
    reg237=reg238+reg237; reg238=reg345*reg70; reg176=reg198+reg176; reg198=reg91*reg109; reg374=reg80*reg374;
    reg264=reg308+reg264; reg308=reg61*reg133; T reg720=reg119*reg94; T reg721=reg121*reg98; T reg722=reg345*reg61;
    reg113=reg87*reg113; reg160=reg71*reg160; T reg723=reg88*reg100; T reg724=reg72*reg155; T reg725=reg71*reg140;
    reg241=reg241+reg242; T reg726=reg93*reg101; T reg727=reg118*reg87; T reg728=reg159*reg63; reg330=reg330*reg61;
    T reg729=reg71*reg125; reg121=reg121*reg94; T reg730=reg95*reg101; T reg731=reg88*reg102; reg269=reg270+reg269;
    reg102=reg87*reg102; reg270=reg71*reg145; reg109=reg93*reg109; reg271=reg271-reg272; T reg732=reg118*reg98;
    reg151=reg80*reg151; T reg733=reg93*reg104; T reg734=reg72*reg141; reg325=reg325+reg328; T reg735=reg88*reg105;
    T reg736=reg87*reg117; reg145=reg78*reg145; reg131=reg71*reg131; reg179=reg78*reg179; reg116=reg213+reg116;
    reg213=reg41*reg98; reg101=reg91*reg101; reg231=reg231-reg232; T reg737=reg80*reg157; reg83=reg91*reg83;
    T reg738=reg63*reg140; T reg739=reg87*reg105; reg119=reg119*reg98; reg345=reg345*reg82; T reg740=reg71*reg143;
    reg118=reg118*reg89; reg130=reg185+reg130; reg120=reg120-reg356; reg185=reg78*reg153; reg154=reg61*reg154;
    T reg741=reg94*reg114; T reg742=reg63*reg133; reg133=reg71*reg133; reg257=reg301+reg257; reg301=reg302+reg301;
    reg125=reg80*reg125; reg302=reg93*reg103; T reg743=reg41*reg87; reg155=reg71*reg155; T reg744=reg95*reg104;
    reg343=reg343+reg380; reg117=reg98*reg117; T reg745=reg98*reg100; T reg746=reg256+reg255; reg143=reg78*reg143;
    T reg747=reg373+reg371; T reg748=reg88*reg115; T reg749=reg69*reg141; reg123=reg123*reg71; T reg750=reg253+reg252;
    reg141=reg79*reg141; reg115=reg95*reg115; reg139=reg72*reg139; T reg751=reg87*reg100; reg100=reg99*reg100;
    reg110=reg93*reg110; reg258=reg258+reg259; reg332=reg346+reg332; reg346=reg159*reg71; T reg752=reg347-reg344;
    reg248=reg248-reg247; reg341=reg80*reg341; T reg753=reg98*reg138; reg103=reg91*reg103; reg227=reg72*reg227;
    reg41=reg41*reg88; reg104=reg91*reg104; reg707=reg707-reg706; reg657=reg656+reg657; reg671=reg311+reg671;
    reg196=reg196+reg508; reg399=reg395+reg399; reg128=reg128+reg339; reg521=reg521-reg520; reg325=reg325+reg732;
    reg311=reg52*reg535; reg199=reg601+reg199; reg525=reg524+reg525; reg395=reg52*reg474; reg524=reg52*reg517;
    reg708=reg703+reg708; reg205=reg712+reg205; reg200=reg458+reg200; reg472=reg471+reg472; reg451=reg668+reg451;
    reg526=reg495+reg526; reg142=reg142+reg518; reg511=reg510+reg511; reg515=reg516-reg515; reg531=reg531-reg532;
    reg666=reg666-reg664; reg681=reg676+reg681; reg720=reg722+reg720; reg514=reg513+reg514; reg101=reg179+reg101;
    reg680=reg717+reg680; reg530=reg530-reg533; reg100=reg100-reg728; reg136=reg136-reg529; reg322=reg366+reg322;
    reg558=reg402+reg558; reg661=reg663+reg661; reg498=reg574-reg498; reg683=reg652+reg683; reg202=reg245+reg202;
    reg119=reg715+reg119; reg130=reg242+reg130; reg674=reg675+reg674; reg701=reg649+reg701; reg410=reg411-reg410;
    reg363=reg363+reg386; reg528=reg527+reg528; reg382=reg732+reg382; reg203=reg203+reg519; reg522=reg523+reg522;
    reg476=reg475+reg476; reg583=reg553-reg583; reg357=reg519+reg357; reg640=reg639+reg640; reg365=reg188+reg365;
    reg465=reg463-reg465; reg642=reg642-reg644; reg362=reg362-reg590; reg752=reg752-reg753; reg179=reg52*reg393;
    reg567=reg588+reg567; reg568=reg569+reg568; reg380=reg332+reg380; reg636=reg637+reg636; reg351=reg508+reg351;
    reg261=reg261-reg599; reg104=reg104-reg185; reg749=reg573-reg749; reg496=reg497+reg496; reg601=reg348+reg601;
    reg461=reg461-reg460; reg458=reg370+reg458; reg188=reg52*reg747; reg555=reg557+reg555; reg375=reg259+reg375;
    reg108=reg108+reg554; reg560=reg561+reg560; reg135=reg135+reg634; reg156=reg223+reg156; reg565=reg566+reg565;
    reg343=reg343+reg117; reg630=reg631+reg630; reg223=reg52*reg369; reg379=reg518+reg379; reg621=reg650+reg621;
    reg400=reg401+reg400; reg624=reg625+reg624; reg546=reg547+reg546; reg409=reg412+reg409; reg744=reg143+reg744;
    reg154=reg154-reg741; reg408=reg405+reg408; reg582=reg582-reg581; reg649=reg336+reg649; reg116=reg116-reg272;
    reg432=reg424+reg432; reg492=reg493+reg492; reg213=reg742+reg213; reg692=reg694+reg692; reg381=reg334+reg381;
    reg132=reg274+reg132; reg120=reg120-reg753; reg745=reg745-reg738; reg158=reg339+reg158; reg534=reg201+reg534;
    reg435=reg435-reg434; reg468=reg468+reg436; reg687=reg688+reg687; reg538=reg539+reg538; reg686=reg714+reg686;
    reg537=reg537-reg536; reg509=reg512+reg509; reg121=reg330+reg121; reg337=reg337-reg204; reg148=reg634+reg148;
    reg730=reg145+reg730; reg721=reg308+reg721; reg499=reg499-reg500; reg192=reg192-reg459; reg593=reg595+reg593;
    reg364=reg386+reg364; reg327=reg327-reg529; reg614=reg614-reg615; reg439=reg438+reg439; reg620=reg592+reg620;
    reg215=reg215-reg502; reg507=reg507-reg501; reg482=reg484+reg482; reg445=reg444+reg445; reg143=reg52*reg453;
    reg612=reg613+reg612; reg333=reg494+reg333; reg145=reg52*reg607; reg331=reg117+reg331; reg486=reg485-reg486;
    reg490=reg491+reg490; reg146=reg197+reg146; reg697=reg605+reg697; reg316=reg316-reg467; reg318=reg318+reg622;
    reg397=reg396-reg397; reg627=reg626+reg627; reg310=reg457+reg310; reg117=reg52*reg313; reg197=reg52*reg314;
    reg398=reg404-reg398; reg629=reg628+reg629; reg633=reg632+reg633; reg307=reg307-reg407; reg115=reg141+reg115;
    reg377=reg690+reg377; reg309=reg367+reg309; reg141=reg52*reg305; reg748=reg123+reg748; reg416=reg416+reg415;
    reg302=reg125+reg302; reg414=reg413+reg414; reg123=reg52*reg417; reg743=reg155+reg743; reg301=reg606+reg301;
    reg125=reg52*reg300; reg103=reg341+reg103; reg384=reg391-reg384; reg390=reg388-reg390; reg41=reg133+reg41;
    reg250=reg250-reg466; reg541=reg541+reg542; reg544=reg544+reg543; reg540=reg149-reg540; reg591=reg591-reg321;
    reg596=reg597+reg596; reg466=reg289-reg466; reg133=reg52*reg594; reg149=reg52*reg124; reg564=reg563-reg564;
    reg287=reg287+reg600; reg155=reg52*reg286; reg281=reg433+reg281; reg603=reg602+reg603; reg559=reg562+reg559;
    reg282=reg456+reg282; reg635=reg604+reg635; reg201=reg52*reg278; reg279=reg279+reg638; reg392=reg556-reg392;
    reg387=reg394-reg387; reg641=reg646+reg641; reg320=reg320-reg477; reg643=reg643-reg645; reg383=reg276-reg383;
    reg389=reg385-reg389; reg648=reg647+reg648; reg418=reg406-reg418; reg651=reg623+reg651; reg403=reg419+reg403;
    reg267=reg319+reg267; reg109=reg151+reg109; reg430=reg429+reg430; reg731=reg729+reg731; reg428=reg431+reg428;
    reg427=reg426+reg427; reg266=reg266-reg280; reg736=reg734+reg736; reg198=reg374+reg198; reg423=reg425+reg423;
    reg422=reg422-reg421; reg263=reg263+reg420; reg177=reg193+reg177; reg264=reg175+reg264; reg545=reg452+reg545;
    reg258=reg258+reg110; reg260=reg288+reg260; reg727=reg724+reg727; reg450=reg450-reg449; reg448=reg447+reg448;
    reg118=reg227+reg118; reg443=reg442+reg443; reg257=reg354+reg257; reg446=reg446-reg441; reg295=reg295-reg440;
    reg151=reg52*reg750; reg175=reg52*reg746; reg139=reg139-reg437; reg751=reg751-reg346; reg265=reg251-reg265;
    reg193=reg462+reg464; reg241=reg241+reg726; reg248=reg248-reg617; reg227=reg52*reg244; reg723=reg723-reg725;
    reg455=reg455+reg456; reg113=reg160+reg113; reg719=reg718+reg719; reg457=reg246+reg457; reg160=reg52*reg236;
    reg237=reg600+reg237; reg352=reg378+reg352; reg481=reg454+reg481; reg479=reg479+reg480; reg83=reg83-reg737;
    reg477=reg234-reg477; reg740=reg740-reg739; reg473=reg478+reg473; reg234=reg52*reg230; reg231=reg638+reg231;
    reg271=reg271+reg733; reg470=reg469-reg470; reg131=reg131-reg735; reg467=reg273-reg467; reg102=reg270+reg102;
    reg269=reg622+reg269; reg433=reg275+reg433; reg368=reg209+reg368; reg700=reg700-reg665; reg164=reg249-reg164;
    reg699=reg667-reg699; reg579=reg579-reg290; reg578=reg577-reg578; reg702=reg598+reg702; reg209=reg52*reg163;
    reg709=reg710+reg709; reg705=reg705-reg704; reg576=reg161+reg576; reg239=reg239+reg575; reg161=reg52*reg183;
    reg184=reg184-reg329; reg180=reg110+reg180; reg572=reg355-reg572; reg110=reg570+reg571; reg713=reg711+reg713;
    reg178=reg233-reg178; reg685=reg684+reg685; reg360=reg181-reg360; reg589=reg587-reg589; reg350=reg176+reg350;
    reg689=reg691+reg689; reg176=reg52*reg174; reg690=reg144+reg690; reg144=reg52*reg189; reg494=reg187+reg494;
    reg228=reg726+reg228; reg670=reg669+reg670; reg489=reg488+reg489; reg673=reg672+reg673; reg186=reg554+reg186;
    reg487=reg238+reg487; reg677=reg678+reg677; reg334=reg220+reg334; reg225=reg328+reg225; reg181=reg52*reg679;
    reg218=reg306-reg218; reg655=reg655-reg654; reg187=reg52*reg483; reg506=reg505-reg506; reg682=reg682-reg653;
    reg220=reg52*reg216; reg347=reg347-reg716; reg660=reg733+reg660; reg503=reg504-reg503; reg233=reg52*reg126;
    reg659=reg658+reg659; reg207=reg207+reg214; reg662=reg212+reg662; reg211=reg420+reg211; reg212=reg52*reg208;
    reg407=reg167-reg407; reg335=reg169-reg335; reg606=reg168+reg606; reg580=reg166-reg580; reg297=reg297-reg449;
    reg552=reg551-reg552; reg696=reg695+reg696; reg166=reg52*reg171; reg170=reg586-reg170; reg609=reg608+reg609;
    reg548=reg548+reg549; reg611=reg610+reg611; reg693=reg693-reg698; reg284=reg345+reg284; reg585=reg584-reg585;
    reg167=reg52*reg550; reg616=reg616-reg617; reg618=reg618-reg619; reg317=reg173+reg317; reg168=reg52*reg303;
    reg172=reg268-reg172; reg749=reg52*reg749; reg730=reg52*reg730; reg433=reg52*reg433; reg169=ponderation*reg395;
    reg281=reg52*reg281; reg109=reg52*reg109; reg681=reg52*reg681; reg364=reg52*reg364; reg660=reg52*reg660;
    reg335=reg52*reg335; reg130=reg52*reg130; reg173=ponderation*reg144; reg503=reg52*reg503; reg446=reg52*reg446;
    reg468=reg52*reg468; reg473=reg52*reg473; reg261=reg52*reg261; reg116=reg52*reg116; reg271=reg52*reg271;
    reg567=reg52*reg567; reg228=reg52*reg228; reg530=reg52*reg530; reg677=reg52*reg677; reg614=reg52*reg614;
    reg671=reg52*reg671; reg399=reg52*reg399; reg258=reg52*reg258; reg450=reg52*reg450; reg596=reg52*reg596;
    reg238=ponderation*reg167; reg186=reg52*reg186; reg423=reg52*reg423; reg101=reg52*reg101; reg618=reg52*reg618;
    reg544=reg52*reg544; reg655=reg52*reg655; reg198=reg52*reg198; reg498=reg52*reg498; reg507=reg52*reg507;
    reg200=reg52*reg200; reg428=reg52*reg428; reg242=ponderation*reg145; reg245=ponderation*reg151; reg246=ponderation*reg187;
    reg249=ponderation*reg149; reg384=reg52*reg384; reg251=ponderation*reg524; reg640=reg52*reg640; reg103=reg52*reg103;
    reg180=reg52*reg180; reg409=reg52*reg409; reg383=reg52*reg383; reg414=reg52*reg414; reg708=reg52*reg708;
    reg302=reg52*reg302; reg263=reg52*reg263; reg572=reg52*reg572; reg317=reg52*reg317; reg375=reg52*reg375;
    reg309=reg52*reg309; reg693=reg52*reg693; reg537=reg52*reg537; reg108=reg52*reg108; reg115=reg52*reg115;
    reg689=reg52*reg689; reg259=ponderation*reg223; reg360=reg52*reg360; reg268=ponderation*reg197; reg403=reg52*reg403;
    reg310=reg52*reg310; reg621=reg52*reg621; reg686=reg52*reg686; reg158=reg52*reg158; reg700=reg52*reg700;
    reg515=reg52*reg515; reg559=reg52*reg559; reg481=reg52*reg481; reg270=ponderation*reg143; reg657=reg52*reg657;
    reg83=reg52*reg83; reg297=reg52*reg297; reg641=reg52*reg641; reg211=reg52*reg211; reg457=reg52*reg457;
    reg104=reg52*reg104; reg170=reg52*reg170; reg486=reg52*reg486; reg719=reg52*reg719; reg583=reg52*reg583;
    reg458=reg52*reg458; reg522=reg52*reg522; reg579=reg52*reg579; reg202=reg52*reg202; reg709=reg52*reg709;
    reg392=reg52*reg392; reg265=reg52*reg265; reg241=reg52*reg241; reg132=reg52*reg132; reg651=reg52*reg651;
    reg465=reg52*reg465; reg576=reg52*reg576; reg744=reg52*reg744; reg207=reg52*reg207; reg422=reg52*reg422;
    reg273=ponderation*reg212; reg643=reg52*reg643; reg674=reg52*reg674; reg164=reg52*reg164; reg266=reg52*reg266;
    reg578=reg52*reg578; reg451=reg52*reg451; reg274=ponderation*reg209; reg279=reg52*reg279; reg239=reg52*reg239;
    reg380=reg52*reg380; reg427=reg52*reg427; reg275=ponderation*reg161; reg670=reg52*reg670; reg276=reg52*reg110;
    reg430=reg52*reg430; reg178=reg52*reg178; reg673=reg52*reg673; reg589=reg52*reg589; reg225=reg52*reg225;
    reg288=ponderation*reg176; reg635=reg52*reg635; reg172=reg52*reg172; reg267=reg52*reg267; reg289=ponderation*reg181;
    reg585=reg52*reg585; reg306=ponderation*reg166; reg148=reg52*reg148; reg196=reg52*reg196; reg295=reg52*reg295;
    reg308=ponderation*reg311; reg205=reg52*reg205; reg531=reg52*reg531; reg443=reg52*reg443; reg136=reg52*reg136;
    reg707=reg52*reg707; reg528=reg52*reg528; reg318=reg52*reg318; reg526=reg52*reg526; reg701=reg52*reg701;
    reg494=reg52*reg494; reg448=reg52*reg448; reg666=reg52*reg666; reg489=reg52*reg489; reg487=reg52*reg487;
    reg260=reg52*reg260; reg661=reg52*reg661; reg334=reg52*reg334; reg648=reg52*reg648; reg218=reg52*reg218;
    reg199=reg52*reg199; reg506=reg52*reg506; reg545=reg52*reg545; reg319=ponderation*reg220; reg683=reg52*reg683;
    reg328=ponderation*reg233; reg680=reg52*reg680; reg320=reg52*reg320; reg713=reg52*reg713; reg389=reg52*reg389;
    reg591=reg52*reg591; reg685=reg52*reg685; reg418=reg52*reg418; reg316=reg52*reg316; reg330=ponderation*reg160;
    reg350=reg52*reg350; reg397=reg52*reg397; reg332=ponderation*reg117; reg690=reg52*reg690; reg398=reg52*reg398;
    reg455=reg52*reg455; reg616=reg52*reg616; reg307=reg52*reg307; reg284=reg52*reg284; reg336=ponderation*reg141;
    reg339=ponderation*reg227; reg696=reg52*reg696; reg416=reg52*reg416; reg341=ponderation*reg123; reg611=reg52*reg611;
    reg606=reg52*reg606; reg345=ponderation*reg125; reg193=reg52*reg193; reg390=reg52*reg390; reg609=reg52*reg609;
    reg250=reg52*reg250; reg682=reg52*reg682; reg467=reg52*reg467; reg580=reg52*reg580; reg603=reg52*reg603;
    reg407=reg52*reg407; reg347=reg52*reg347; reg470=reg52*reg470; reg552=reg52*reg552; reg659=reg52*reg659;
    reg548=reg52*reg548; reg348=ponderation*reg168; reg662=reg52*reg662; reg287=reg52*reg287; reg541=reg52*reg541;
    reg354=ponderation*reg234; reg540=reg52*reg540; reg368=reg52*reg368; reg466=reg52*reg466; reg699=reg52*reg699;
    reg564=reg52*reg564; reg477=reg52*reg477; reg702=reg52*reg702; reg355=ponderation*reg155; reg282=reg52*reg282;
    reg705=reg52*reg705; reg366=ponderation*reg133; reg370=ponderation*reg201; reg479=reg52*reg479; reg387=reg52*reg387;
    reg184=reg52*reg184; reg461=reg52*reg461; reg41=reg52*reg41; reg131=reg52*reg131; reg374=ponderation*reg179;
    reg721=reg52*reg721; reg365=reg52*reg365; reg752=reg52*reg752; reg408=reg52*reg408; reg102=reg52*reg102;
    reg400=reg52*reg400; reg301=reg52*reg301; reg154=reg52*reg154; reg156=reg52*reg156; reg555=reg52*reg555;
    reg343=reg52*reg343; reg560=reg52*reg560; reg269=reg52*reg269; reg378=ponderation*reg188; reg565=reg52*reg565;
    reg379=reg52*reg379; reg135=reg52*reg135; reg731=reg52*reg731; reg546=reg52*reg546; reg743=reg52*reg743;
    reg630=reg52*reg630; reg582=reg52*reg582; reg357=reg52*reg357; reg624=reg52*reg624; reg736=reg52*reg736;
    reg723=reg52*reg723; reg128=reg52*reg128; reg113=reg52*reg113; reg558=reg52*reg558; reg237=reg52*reg237;
    reg410=reg52*reg410; reg248=reg52*reg248; reg363=reg52*reg363; reg720=reg52*reg720; reg322=reg52*reg322;
    reg352=reg52*reg352; reg476=reg52*reg476; reg119=reg52*reg119; reg472=reg52*reg472; reg751=reg52*reg751;
    reg740=reg52*reg740; reg337=reg52*reg337; reg325=reg52*reg325; reg435=reg52*reg435; reg121=reg52*reg121;
    reg432=reg52*reg432; reg100=reg52*reg100; reg146=reg52*reg146; reg445=reg52*reg445; reg213=reg52*reg213;
    reg231=reg52*reg231; reg439=reg52*reg439; reg120=reg52*reg120; reg192=reg52*reg192; reg382=reg52*reg382;
    reg593=reg52*reg593; reg333=reg52*reg333; reg490=reg52*reg490; reg620=reg52*reg620; reg118=reg52*reg118;
    reg492=reg52*reg492; reg633=reg52*reg633; reg381=reg52*reg381; reg612=reg52*reg612; reg534=reg52*reg534;
    reg331=reg52*reg331; reg538=reg52*reg538; reg257=reg52*reg257; reg509=reg52*reg509; reg697=reg52*reg697;
    reg142=reg52*reg142; reg385=ponderation*reg175; reg525=reg52*reg525; reg692=reg52*reg692; reg629=reg52*reg629;
    reg521=reg52*reg521; reg745=reg52*reg745; reg203=reg52*reg203; reg139=reg52*reg139; reg514=reg52*reg514;
    reg687=reg52*reg687; reg511=reg52*reg511; reg627=reg52*reg627; reg496=reg52*reg496; reg642=reg52*reg642;
    reg351=reg52*reg351; reg177=reg52*reg177; reg499=reg52*reg499; reg636=reg52*reg636; reg649=reg52*reg649;
    reg327=reg52*reg327; reg568=reg52*reg568; reg264=reg52*reg264; reg727=reg52*reg727; reg362=reg52*reg362;
    reg748=reg52*reg748; reg482=reg52*reg482; reg377=reg52*reg377; reg215=reg52*reg215; reg601=reg52*reg601;
    T tmp_3_5=ponderation*reg611; T tmp_4_2=ponderation*reg748; T tmp_4_6=ponderation*reg751; T tmp_3_16=ponderation*reg627; T tmp_3_15=ponderation*reg318;
    T tmp_3_11=ponderation*reg635; T tmp_3_10=ponderation*reg603; T tmp_16_1=ponderation*reg310; T tmp_4_7=ponderation*reg248; T tmp_16_6=ponderation*reg618;
    T tmp_16_3=ponderation*reg641; T tmp_16_4=ponderation*reg281; T tmp_4_4=ponderation*reg301; T tmp_4_0=ponderation*reg633; T tmp_16_2=ponderation*reg651;
    T tmp_3_8=-reg366; T tmp_3_13=ponderation*reg643; T tmp_15_16=ponderation*reg103; T tmp_3_9=ponderation*reg287; T tmp_3_7=ponderation*reg591;
    T tmp_16_0=ponderation*reg115; T tmp_4_5=ponderation*reg41; T tmp_3_12=ponderation*reg279; T tmp_16_5=ponderation*reg596; T tmp_3_14=ponderation*reg648;
    T tmp_15_15=ponderation*reg241; T tmp_4_1=ponderation*reg377; T tmp_15_17=ponderation*reg302; T tmp_3_6=ponderation*reg616; T tmp_3_17=ponderation*reg629;
    T tmp_4_3=ponderation*reg743; T tmp_0_10=ponderation*reg630; T tmp_0_11=ponderation*reg624; T tmp_17_9=ponderation*reg640; T tmp_0_12=ponderation*reg649;
    T tmp_0_13=ponderation*reg642; T tmp_17_8=ponderation*reg261; T tmp_0_14=ponderation*reg636; T tmp_0_15=ponderation*reg601; T tmp_17_7=ponderation*reg614;
    T tmp_0_16=ponderation*reg593; T tmp_0_17=ponderation*reg620; T tmp_17_6=-reg242; T tmp_1_0=ponderation*reg612; T tmp_1_1=ponderation*reg331;
    T tmp_17_5=ponderation*reg132; T tmp_1_2=ponderation*reg697; T tmp_1_3=ponderation*reg692; T tmp_17_4=ponderation*reg686; T tmp_1_8=ponderation*reg745;
    T tmp_1_9=ponderation*reg687; T tmp_17_3=ponderation*reg708; T tmp_1_10=ponderation*reg148; T tmp_17_17=ponderation*reg130; T tmp_1_6=ponderation*reg100;
    T tmp_1_5=ponderation*reg213; T tmp_1_4=ponderation*reg382; T tmp_0_0=ponderation*reg343; T tmp_17_16=ponderation*reg101; T tmp_0_1=ponderation*reg720;
    T tmp_0_2=ponderation*reg119; T tmp_17_15=ponderation*reg730; T tmp_0_3=ponderation*reg325; T tmp_17_14=ponderation*reg116; T tmp_1_7=ponderation*reg120;
    T tmp_17_13=ponderation*reg104; T tmp_0_4=ponderation*reg121; T tmp_0_5=ponderation*reg721; T tmp_17_12=ponderation*reg744; T tmp_0_6=ponderation*reg752;
    T tmp_0_7=ponderation*reg154; T tmp_17_11=ponderation*reg375; T tmp_0_8=-reg378; T tmp_0_9=ponderation*reg135; T tmp_17_10=ponderation*reg621;
    T tmp_16_13=ponderation*reg660; T tmp_2_8=ponderation*reg347; T tmp_2_9=ponderation*reg659; T tmp_16_12=ponderation*reg700; T tmp_2_10=ponderation*reg662;
    T tmp_2_11=ponderation*reg368; T tmp_16_11=ponderation*reg709; T tmp_2_12=ponderation*reg702; T tmp_2_13=ponderation*reg705; T tmp_16_10=ponderation*reg180;
    T tmp_2_14=ponderation*reg184; T tmp_2_15=ponderation*reg713; T tmp_16_9=ponderation*reg689; T tmp_2_16=ponderation*reg685; T tmp_2_17=ponderation*reg350;
    T tmp_16_8=ponderation*reg693; T tmp_3_0=ponderation*reg690; T tmp_3_1=ponderation*reg284; T tmp_16_7=ponderation*reg297; T tmp_3_2=ponderation*reg696;
    T tmp_3_3=ponderation*reg606; T tmp_3_4=ponderation*reg609; T tmp_1_11=ponderation*reg205; T tmp_17_2=ponderation*reg202; T tmp_1_12=ponderation*reg707;
    T tmp_1_13=ponderation*reg701; T tmp_17_1=ponderation*reg657; T tmp_1_14=ponderation*reg666; T tmp_1_15=ponderation*reg661; T tmp_17_0=ponderation*reg681;
    T tmp_1_16=ponderation*reg199; T tmp_1_17=ponderation*reg683; T tmp_16_17=ponderation*reg671; T tmp_2_0=ponderation*reg680; T tmp_2_1=ponderation*reg674;
    T tmp_16_16=ponderation*reg228; T tmp_5_12=ponderation*reg451; T tmp_2_2=ponderation*reg380; T tmp_2_3=ponderation*reg670; T tmp_16_15=ponderation*reg677;
    T tmp_2_4=ponderation*reg673; T tmp_2_5=ponderation*reg225; T tmp_16_14=ponderation*reg655; T tmp_2_6=-reg289; T tmp_2_7=ponderation*reg682;
    T tmp_9_5=ponderation*reg528; T tmp_9_6=ponderation*reg136; T tmp_13_10=ponderation*reg200; T tmp_9_7=ponderation*reg531; T tmp_9_8=-reg308;
    T tmp_9_9=ponderation*reg196; T tmp_13_9=ponderation*reg515; T tmp_9_10=ponderation*reg511; T tmp_9_11=ponderation*reg514; T tmp_13_8=ponderation*reg522;
    T tmp_9_12=ponderation*reg203; T tmp_9_13=ponderation*reg521; T tmp_13_7=-reg251; T tmp_9_14=ponderation*reg525; T tmp_9_15=ponderation*reg142;
    T tmp_13_5=ponderation*reg537; T tmp_9_16=ponderation*reg509; T tmp_9_17=ponderation*reg538; T tmp_13_4=ponderation*reg158; T tmp_10_0=ponderation*reg534;
    T tmp_10_1=ponderation*reg381; T tmp_13_3=ponderation*reg486; T tmp_10_2=ponderation*reg492; T tmp_10_3=ponderation*reg490; T tmp_13_2=ponderation*reg507;
    T tmp_8_5=ponderation*reg178; T tmp_8_6=ponderation*reg276; T tmp_14_0=ponderation*reg576; T tmp_8_7=-reg275; T tmp_8_8=ponderation*reg239;
    T tmp_13_17=ponderation*reg579; T tmp_8_9=-reg274; T tmp_8_10=ponderation*reg578; T tmp_13_16=ponderation*reg211; T tmp_8_11=ponderation*reg164;
    T tmp_8_12=-reg273; T tmp_8_13=ponderation*reg207; T tmp_13_15=ponderation*reg503; T tmp_8_14=-reg328; T tmp_8_15=-reg319;
    T tmp_13_14=-reg246; T tmp_8_16=ponderation*reg506; T tmp_8_17=ponderation*reg218; T tmp_13_13=ponderation*reg186; T tmp_9_0=ponderation*reg334;
    T tmp_9_1=ponderation*reg487; T tmp_9_2=ponderation*reg489; T tmp_13_12=-reg173; T tmp_9_3=ponderation*reg494; T tmp_9_4=ponderation*reg526;
    T tmp_13_11=ponderation*reg530; T tmp_11_4=ponderation*reg408; T tmp_12_10=ponderation*reg465; T tmp_11_5=ponderation*reg365; T tmp_11_6=-reg374;
    T tmp_12_9=ponderation*reg458; T tmp_11_7=ponderation*reg461; T tmp_11_8=ponderation*reg192; T tmp_11_9=ponderation*reg439; T tmp_12_8=-reg270;
    T tmp_11_10=ponderation*reg445; T tmp_11_11=ponderation*reg146; T tmp_12_7=ponderation*reg468; T tmp_11_12=ponderation*reg432; T tmp_11_13=ponderation*reg435;
    T tmp_12_6=-reg169; T tmp_11_14=ponderation*reg337; T tmp_11_15=ponderation*reg472; T tmp_11_16=ponderation*reg476; T tmp_11_17=ponderation*reg322;
    T tmp_12_5=ponderation*reg399; T tmp_12_0=ponderation*reg363; T tmp_12_1=ponderation*reg410; T tmp_12_4=ponderation*reg498; T tmp_12_2=ponderation*reg558;
    T tmp_12_3=ponderation*reg128; T tmp_10_4=ponderation*reg333; T tmp_10_5=ponderation*reg482; T tmp_13_1=ponderation*reg364; T tmp_10_6=ponderation*reg215;
    T tmp_10_7=ponderation*reg327; T tmp_10_8=ponderation*reg499; T tmp_13_0=ponderation*reg749; T tmp_10_9=ponderation*reg496; T tmp_12_17=ponderation*reg567;
    T tmp_10_10=ponderation*reg351; T tmp_10_11=ponderation*reg568; T tmp_12_16=ponderation*reg583; T tmp_10_12=ponderation*reg362; T tmp_10_13=ponderation*reg357;
    T tmp_12_15=ponderation*reg263; T tmp_10_14=ponderation*reg582; T tmp_12_13=-reg259; T tmp_10_15=ponderation*reg546; T tmp_10_16=ponderation*reg379;
    T tmp_10_17=ponderation*reg565; T tmp_12_12=ponderation*reg108; T tmp_11_0=ponderation*reg560; T tmp_11_1=ponderation*reg555; T tmp_11_2=ponderation*reg156;
    T tmp_12_11=ponderation*reg409; T tmp_11_3=ponderation*reg400; T tmp_5_8=ponderation*reg295; T tmp_15_6=ponderation*reg450; T tmp_5_9=ponderation*reg443;
    T tmp_5_10=ponderation*reg448; T tmp_5_11=ponderation*reg260; T tmp_15_5=ponderation*reg423; T tmp_12_14=ponderation*reg545; T tmp_5_13=ponderation*reg422;
    T tmp_15_4=ponderation*reg428; T tmp_5_14=ponderation*reg266; T tmp_5_15=ponderation*reg427; T tmp_15_3=ponderation*reg433; T tmp_5_16=ponderation*reg430;
    T tmp_5_17=ponderation*reg267; T tmp_6_0=ponderation*reg467; T tmp_15_2=ponderation*reg473; T tmp_6_1=ponderation*reg470; T tmp_6_2=-reg354;
    T tmp_15_1=ponderation*reg481; T tmp_6_3=ponderation*reg477; T tmp_13_6=ponderation*reg479; T tmp_6_4=ponderation*reg699; T tmp_15_0=ponderation*reg457;
    T tmp_6_5=-reg330; T tmp_6_6=ponderation*reg455; T tmp_4_8=ponderation*reg723; T tmp_15_14=ponderation*reg719; T tmp_4_9=ponderation*reg113;
    T tmp_4_10=ponderation*reg237; T tmp_15_13=ponderation*reg83; T tmp_4_11=ponderation*reg352; T tmp_4_12=ponderation*reg740; T tmp_15_12=ponderation*reg271;
    T tmp_4_13=ponderation*reg231; T tmp_4_14=ponderation*reg131; T tmp_4_15=ponderation*reg102; T tmp_15_11=ponderation*reg109; T tmp_4_16=ponderation*reg269;
    T tmp_4_17=ponderation*reg731; T tmp_15_10=ponderation*reg198; T tmp_5_0=ponderation*reg736; T tmp_5_1=ponderation*reg177; T tmp_15_9=ponderation*reg258;
    T tmp_5_2=ponderation*reg264; T tmp_5_3=ponderation*reg727; T tmp_5_4=ponderation*reg118; T tmp_15_8=-reg245; T tmp_5_5=ponderation*reg257;
    T tmp_5_6=-reg385; T tmp_15_7=ponderation*reg446; T tmp_5_7=ponderation*reg139; T tmp_7_6=-reg370; T tmp_7_7=ponderation*reg282;
    T tmp_14_8=-reg249; T tmp_7_8=-reg355; T tmp_7_9=ponderation*reg564; T tmp_7_10=ponderation*reg466; T tmp_14_7=ponderation*reg544;
    T tmp_7_11=ponderation*reg540; T tmp_7_12=ponderation*reg541; T tmp_14_6=-reg238; T tmp_7_13=-reg348; T tmp_7_14=ponderation*reg548;
    T tmp_14_5=ponderation*reg335; T tmp_7_15=ponderation*reg552; T tmp_7_16=ponderation*reg407; T tmp_14_4=ponderation*reg170; T tmp_7_17=ponderation*reg580;
    T tmp_8_0=-reg306; T tmp_14_3=ponderation*reg317; T tmp_8_1=ponderation*reg585; T tmp_8_2=ponderation*reg172; T tmp_14_2=ponderation*reg360;
    T tmp_8_3=-reg288; T tmp_8_4=ponderation*reg589; T tmp_14_1=ponderation*reg572; T tmp_6_7=-reg339; T tmp_14_17=ponderation*reg265;
    T tmp_6_8=ponderation*reg193; T tmp_6_9=ponderation*reg250; T tmp_14_16=ponderation*reg384; T tmp_6_10=ponderation*reg390; T tmp_6_11=-reg345;
    T tmp_14_15=ponderation*reg414; T tmp_6_12=-reg341; T tmp_6_13=ponderation*reg416; T tmp_14_14=ponderation*reg309; T tmp_6_14=-reg336;
    T tmp_6_15=ponderation*reg307; T tmp_14_13=-reg268; T tmp_6_16=ponderation*reg398; T tmp_6_17=-reg332; T tmp_14_12=ponderation*reg403;
    T tmp_7_0=ponderation*reg397; T tmp_7_1=ponderation*reg316; T tmp_14_11=ponderation*reg383; T tmp_7_2=ponderation*reg418; T tmp_7_3=ponderation*reg389;
    T tmp_14_10=ponderation*reg392; T tmp_7_4=ponderation*reg320; T tmp_7_5=ponderation*reg387; T tmp_14_9=ponderation*reg559;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[0]+2) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_3;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_4;
    matrix(indices[0]+0,indices[1]+2) += tmp_0_5;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_7;
    matrix(indices[0]+0,indices[2]+2) += tmp_0_8;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_9;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_10;
    matrix(indices[0]+0,indices[3]+2) += tmp_0_11;
    matrix(indices[0]+0,indices[4]+0) += tmp_0_12;
    matrix(indices[0]+0,indices[4]+1) += tmp_0_13;
    matrix(indices[0]+0,indices[4]+2) += tmp_0_14;
    matrix(indices[0]+0,indices[5]+0) += tmp_0_15;
    matrix(indices[0]+0,indices[5]+1) += tmp_0_16;
    matrix(indices[0]+0,indices[5]+2) += tmp_0_17;
    matrix(indices[0]+1,indices[0]+0) += tmp_1_0;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[0]+2) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_3;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_4;
    matrix(indices[0]+1,indices[1]+2) += tmp_1_5;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_7;
    matrix(indices[0]+1,indices[2]+2) += tmp_1_8;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_9;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_10;
    matrix(indices[0]+1,indices[3]+2) += tmp_1_11;
    matrix(indices[0]+1,indices[4]+0) += tmp_1_12;
    matrix(indices[0]+1,indices[4]+1) += tmp_1_13;
    matrix(indices[0]+1,indices[4]+2) += tmp_1_14;
    matrix(indices[0]+1,indices[5]+0) += tmp_1_15;
    matrix(indices[0]+1,indices[5]+1) += tmp_1_16;
    matrix(indices[0]+1,indices[5]+2) += tmp_1_17;
    matrix(indices[0]+2,indices[0]+0) += tmp_2_0;
    matrix(indices[0]+2,indices[0]+1) += tmp_2_1;
    matrix(indices[0]+2,indices[0]+2) += tmp_2_2;
    matrix(indices[0]+2,indices[1]+0) += tmp_2_3;
    matrix(indices[0]+2,indices[1]+1) += tmp_2_4;
    matrix(indices[0]+2,indices[1]+2) += tmp_2_5;
    matrix(indices[0]+2,indices[2]+0) += tmp_2_6;
    matrix(indices[0]+2,indices[2]+1) += tmp_2_7;
    matrix(indices[0]+2,indices[2]+2) += tmp_2_8;
    matrix(indices[0]+2,indices[3]+0) += tmp_2_9;
    matrix(indices[0]+2,indices[3]+1) += tmp_2_10;
    matrix(indices[0]+2,indices[3]+2) += tmp_2_11;
    matrix(indices[0]+2,indices[4]+0) += tmp_2_12;
    matrix(indices[0]+2,indices[4]+1) += tmp_2_13;
    matrix(indices[0]+2,indices[4]+2) += tmp_2_14;
    matrix(indices[0]+2,indices[5]+0) += tmp_2_15;
    matrix(indices[0]+2,indices[5]+1) += tmp_2_16;
    matrix(indices[0]+2,indices[5]+2) += tmp_2_17;
    matrix(indices[1]+0,indices[0]+0) += tmp_3_0;
    matrix(indices[1]+0,indices[0]+1) += tmp_3_1;
    matrix(indices[1]+0,indices[0]+2) += tmp_3_2;
    matrix(indices[1]+0,indices[1]+0) += tmp_3_3;
    matrix(indices[1]+0,indices[1]+1) += tmp_3_4;
    matrix(indices[1]+0,indices[1]+2) += tmp_3_5;
    matrix(indices[1]+0,indices[2]+0) += tmp_3_6;
    matrix(indices[1]+0,indices[2]+1) += tmp_3_7;
    matrix(indices[1]+0,indices[2]+2) += tmp_3_8;
    matrix(indices[1]+0,indices[3]+0) += tmp_3_9;
    matrix(indices[1]+0,indices[3]+1) += tmp_3_10;
    matrix(indices[1]+0,indices[3]+2) += tmp_3_11;
    matrix(indices[1]+0,indices[4]+0) += tmp_3_12;
    matrix(indices[1]+0,indices[4]+1) += tmp_3_13;
    matrix(indices[1]+0,indices[4]+2) += tmp_3_14;
    matrix(indices[1]+0,indices[5]+0) += tmp_3_15;
    matrix(indices[1]+0,indices[5]+1) += tmp_3_16;
    matrix(indices[1]+0,indices[5]+2) += tmp_3_17;
    matrix(indices[1]+1,indices[0]+0) += tmp_4_0;
    matrix(indices[1]+1,indices[0]+1) += tmp_4_1;
    matrix(indices[1]+1,indices[0]+2) += tmp_4_2;
    matrix(indices[1]+1,indices[1]+0) += tmp_4_3;
    matrix(indices[1]+1,indices[1]+1) += tmp_4_4;
    matrix(indices[1]+1,indices[1]+2) += tmp_4_5;
    matrix(indices[1]+1,indices[2]+0) += tmp_4_6;
    matrix(indices[1]+1,indices[2]+1) += tmp_4_7;
    matrix(indices[1]+1,indices[2]+2) += tmp_4_8;
    matrix(indices[1]+1,indices[3]+0) += tmp_4_9;
    matrix(indices[1]+1,indices[3]+1) += tmp_4_10;
    matrix(indices[1]+1,indices[3]+2) += tmp_4_11;
    matrix(indices[1]+1,indices[4]+0) += tmp_4_12;
    matrix(indices[1]+1,indices[4]+1) += tmp_4_13;
    matrix(indices[1]+1,indices[4]+2) += tmp_4_14;
    matrix(indices[1]+1,indices[5]+0) += tmp_4_15;
    matrix(indices[1]+1,indices[5]+1) += tmp_4_16;
    matrix(indices[1]+1,indices[5]+2) += tmp_4_17;
    matrix(indices[1]+2,indices[0]+0) += tmp_5_0;
    matrix(indices[1]+2,indices[0]+1) += tmp_5_1;
    matrix(indices[1]+2,indices[0]+2) += tmp_5_2;
    matrix(indices[1]+2,indices[1]+0) += tmp_5_3;
    matrix(indices[1]+2,indices[1]+1) += tmp_5_4;
    matrix(indices[1]+2,indices[1]+2) += tmp_5_5;
    matrix(indices[1]+2,indices[2]+0) += tmp_5_6;
    matrix(indices[1]+2,indices[2]+1) += tmp_5_7;
    matrix(indices[1]+2,indices[2]+2) += tmp_5_8;
    matrix(indices[1]+2,indices[3]+0) += tmp_5_9;
    matrix(indices[1]+2,indices[3]+1) += tmp_5_10;
    matrix(indices[1]+2,indices[3]+2) += tmp_5_11;
    matrix(indices[1]+2,indices[4]+0) += tmp_5_12;
    matrix(indices[1]+2,indices[4]+1) += tmp_5_13;
    matrix(indices[1]+2,indices[4]+2) += tmp_5_14;
    matrix(indices[1]+2,indices[5]+0) += tmp_5_15;
    matrix(indices[1]+2,indices[5]+1) += tmp_5_16;
    matrix(indices[1]+2,indices[5]+2) += tmp_5_17;
    matrix(indices[2]+0,indices[0]+0) += tmp_6_0;
    matrix(indices[2]+0,indices[0]+1) += tmp_6_1;
    matrix(indices[2]+0,indices[0]+2) += tmp_6_2;
    matrix(indices[2]+0,indices[1]+0) += tmp_6_3;
    matrix(indices[2]+0,indices[1]+1) += tmp_6_4;
    matrix(indices[2]+0,indices[1]+2) += tmp_6_5;
    matrix(indices[2]+0,indices[2]+0) += tmp_6_6;
    matrix(indices[2]+0,indices[2]+1) += tmp_6_7;
    matrix(indices[2]+0,indices[2]+2) += tmp_6_8;
    matrix(indices[2]+0,indices[3]+0) += tmp_6_9;
    matrix(indices[2]+0,indices[3]+1) += tmp_6_10;
    matrix(indices[2]+0,indices[3]+2) += tmp_6_11;
    matrix(indices[2]+0,indices[4]+0) += tmp_6_12;
    matrix(indices[2]+0,indices[4]+1) += tmp_6_13;
    matrix(indices[2]+0,indices[4]+2) += tmp_6_14;
    matrix(indices[2]+0,indices[5]+0) += tmp_6_15;
    matrix(indices[2]+0,indices[5]+1) += tmp_6_16;
    matrix(indices[2]+0,indices[5]+2) += tmp_6_17;
    matrix(indices[2]+1,indices[0]+0) += tmp_7_0;
    matrix(indices[2]+1,indices[0]+1) += tmp_7_1;
    matrix(indices[2]+1,indices[0]+2) += tmp_7_2;
    matrix(indices[2]+1,indices[1]+0) += tmp_7_3;
    matrix(indices[2]+1,indices[1]+1) += tmp_7_4;
    matrix(indices[2]+1,indices[1]+2) += tmp_7_5;
    matrix(indices[2]+1,indices[2]+0) += tmp_7_6;
    matrix(indices[2]+1,indices[2]+1) += tmp_7_7;
    matrix(indices[2]+1,indices[2]+2) += tmp_7_8;
    matrix(indices[2]+1,indices[3]+0) += tmp_7_9;
    matrix(indices[2]+1,indices[3]+1) += tmp_7_10;
    matrix(indices[2]+1,indices[3]+2) += tmp_7_11;
    matrix(indices[2]+1,indices[4]+0) += tmp_7_12;
    matrix(indices[2]+1,indices[4]+1) += tmp_7_13;
    matrix(indices[2]+1,indices[4]+2) += tmp_7_14;
    matrix(indices[2]+1,indices[5]+0) += tmp_7_15;
    matrix(indices[2]+1,indices[5]+1) += tmp_7_16;
    matrix(indices[2]+1,indices[5]+2) += tmp_7_17;
    matrix(indices[2]+2,indices[0]+0) += tmp_8_0;
    matrix(indices[2]+2,indices[0]+1) += tmp_8_1;
    matrix(indices[2]+2,indices[0]+2) += tmp_8_2;
    matrix(indices[2]+2,indices[1]+0) += tmp_8_3;
    matrix(indices[2]+2,indices[1]+1) += tmp_8_4;
    matrix(indices[2]+2,indices[1]+2) += tmp_8_5;
    matrix(indices[2]+2,indices[2]+0) += tmp_8_6;
    matrix(indices[2]+2,indices[2]+1) += tmp_8_7;
    matrix(indices[2]+2,indices[2]+2) += tmp_8_8;
    matrix(indices[2]+2,indices[3]+0) += tmp_8_9;
    matrix(indices[2]+2,indices[3]+1) += tmp_8_10;
    matrix(indices[2]+2,indices[3]+2) += tmp_8_11;
    matrix(indices[2]+2,indices[4]+0) += tmp_8_12;
    matrix(indices[2]+2,indices[4]+1) += tmp_8_13;
    matrix(indices[2]+2,indices[4]+2) += tmp_8_14;
    matrix(indices[2]+2,indices[5]+0) += tmp_8_15;
    matrix(indices[2]+2,indices[5]+1) += tmp_8_16;
    matrix(indices[2]+2,indices[5]+2) += tmp_8_17;
    matrix(indices[3]+0,indices[0]+0) += tmp_9_0;
    matrix(indices[3]+0,indices[0]+1) += tmp_9_1;
    matrix(indices[3]+0,indices[0]+2) += tmp_9_2;
    matrix(indices[3]+0,indices[1]+0) += tmp_9_3;
    matrix(indices[3]+0,indices[1]+1) += tmp_9_4;
    matrix(indices[3]+0,indices[1]+2) += tmp_9_5;
    matrix(indices[3]+0,indices[2]+0) += tmp_9_6;
    matrix(indices[3]+0,indices[2]+1) += tmp_9_7;
    matrix(indices[3]+0,indices[2]+2) += tmp_9_8;
    matrix(indices[3]+0,indices[3]+0) += tmp_9_9;
    matrix(indices[3]+0,indices[3]+1) += tmp_9_10;
    matrix(indices[3]+0,indices[3]+2) += tmp_9_11;
    matrix(indices[3]+0,indices[4]+0) += tmp_9_12;
    matrix(indices[3]+0,indices[4]+1) += tmp_9_13;
    matrix(indices[3]+0,indices[4]+2) += tmp_9_14;
    matrix(indices[3]+0,indices[5]+0) += tmp_9_15;
    matrix(indices[3]+0,indices[5]+1) += tmp_9_16;
    matrix(indices[3]+0,indices[5]+2) += tmp_9_17;
    matrix(indices[3]+1,indices[0]+0) += tmp_10_0;
    matrix(indices[3]+1,indices[0]+1) += tmp_10_1;
    matrix(indices[3]+1,indices[0]+2) += tmp_10_2;
    matrix(indices[3]+1,indices[1]+0) += tmp_10_3;
    matrix(indices[3]+1,indices[1]+1) += tmp_10_4;
    matrix(indices[3]+1,indices[1]+2) += tmp_10_5;
    matrix(indices[3]+1,indices[2]+0) += tmp_10_6;
    matrix(indices[3]+1,indices[2]+1) += tmp_10_7;
    matrix(indices[3]+1,indices[2]+2) += tmp_10_8;
    matrix(indices[3]+1,indices[3]+0) += tmp_10_9;
    matrix(indices[3]+1,indices[3]+1) += tmp_10_10;
    matrix(indices[3]+1,indices[3]+2) += tmp_10_11;
    matrix(indices[3]+1,indices[4]+0) += tmp_10_12;
    matrix(indices[3]+1,indices[4]+1) += tmp_10_13;
    matrix(indices[3]+1,indices[4]+2) += tmp_10_14;
    matrix(indices[3]+1,indices[5]+0) += tmp_10_15;
    matrix(indices[3]+1,indices[5]+1) += tmp_10_16;
    matrix(indices[3]+1,indices[5]+2) += tmp_10_17;
    matrix(indices[3]+2,indices[0]+0) += tmp_11_0;
    matrix(indices[3]+2,indices[0]+1) += tmp_11_1;
    matrix(indices[3]+2,indices[0]+2) += tmp_11_2;
    matrix(indices[3]+2,indices[1]+0) += tmp_11_3;
    matrix(indices[3]+2,indices[1]+1) += tmp_11_4;
    matrix(indices[3]+2,indices[1]+2) += tmp_11_5;
    matrix(indices[3]+2,indices[2]+0) += tmp_11_6;
    matrix(indices[3]+2,indices[2]+1) += tmp_11_7;
    matrix(indices[3]+2,indices[2]+2) += tmp_11_8;
    matrix(indices[3]+2,indices[3]+0) += tmp_11_9;
    matrix(indices[3]+2,indices[3]+1) += tmp_11_10;
    matrix(indices[3]+2,indices[3]+2) += tmp_11_11;
    matrix(indices[3]+2,indices[4]+0) += tmp_11_12;
    matrix(indices[3]+2,indices[4]+1) += tmp_11_13;
    matrix(indices[3]+2,indices[4]+2) += tmp_11_14;
    matrix(indices[3]+2,indices[5]+0) += tmp_11_15;
    matrix(indices[3]+2,indices[5]+1) += tmp_11_16;
    matrix(indices[3]+2,indices[5]+2) += tmp_11_17;
    matrix(indices[4]+0,indices[0]+0) += tmp_12_0;
    matrix(indices[4]+0,indices[0]+1) += tmp_12_1;
    matrix(indices[4]+0,indices[0]+2) += tmp_12_2;
    matrix(indices[4]+0,indices[1]+0) += tmp_12_3;
    matrix(indices[4]+0,indices[1]+1) += tmp_12_4;
    matrix(indices[4]+0,indices[1]+2) += tmp_12_5;
    matrix(indices[4]+0,indices[2]+0) += tmp_12_6;
    matrix(indices[4]+0,indices[2]+1) += tmp_12_7;
    matrix(indices[4]+0,indices[2]+2) += tmp_12_8;
    matrix(indices[4]+0,indices[3]+0) += tmp_12_9;
    matrix(indices[4]+0,indices[3]+1) += tmp_12_10;
    matrix(indices[4]+0,indices[3]+2) += tmp_12_11;
    matrix(indices[4]+0,indices[4]+0) += tmp_12_12;
    matrix(indices[4]+0,indices[4]+1) += tmp_12_13;
    matrix(indices[4]+0,indices[4]+2) += tmp_12_14;
    matrix(indices[4]+0,indices[5]+0) += tmp_12_15;
    matrix(indices[4]+0,indices[5]+1) += tmp_12_16;
    matrix(indices[4]+0,indices[5]+2) += tmp_12_17;
    matrix(indices[4]+1,indices[0]+0) += tmp_13_0;
    matrix(indices[4]+1,indices[0]+1) += tmp_13_1;
    matrix(indices[4]+1,indices[0]+2) += tmp_13_2;
    matrix(indices[4]+1,indices[1]+0) += tmp_13_3;
    matrix(indices[4]+1,indices[1]+1) += tmp_13_4;
    matrix(indices[4]+1,indices[1]+2) += tmp_13_5;
    matrix(indices[4]+1,indices[2]+0) += tmp_13_6;
    matrix(indices[4]+1,indices[2]+1) += tmp_13_7;
    matrix(indices[4]+1,indices[2]+2) += tmp_13_8;
    matrix(indices[4]+1,indices[3]+0) += tmp_13_9;
    matrix(indices[4]+1,indices[3]+1) += tmp_13_10;
    matrix(indices[4]+1,indices[3]+2) += tmp_13_11;
    matrix(indices[4]+1,indices[4]+0) += tmp_13_12;
    matrix(indices[4]+1,indices[4]+1) += tmp_13_13;
    matrix(indices[4]+1,indices[4]+2) += tmp_13_14;
    matrix(indices[4]+1,indices[5]+0) += tmp_13_15;
    matrix(indices[4]+1,indices[5]+1) += tmp_13_16;
    matrix(indices[4]+1,indices[5]+2) += tmp_13_17;
    matrix(indices[4]+2,indices[0]+0) += tmp_14_0;
    matrix(indices[4]+2,indices[0]+1) += tmp_14_1;
    matrix(indices[4]+2,indices[0]+2) += tmp_14_2;
    matrix(indices[4]+2,indices[1]+0) += tmp_14_3;
    matrix(indices[4]+2,indices[1]+1) += tmp_14_4;
    matrix(indices[4]+2,indices[1]+2) += tmp_14_5;
    matrix(indices[4]+2,indices[2]+0) += tmp_14_6;
    matrix(indices[4]+2,indices[2]+1) += tmp_14_7;
    matrix(indices[4]+2,indices[2]+2) += tmp_14_8;
    matrix(indices[4]+2,indices[3]+0) += tmp_14_9;
    matrix(indices[4]+2,indices[3]+1) += tmp_14_10;
    matrix(indices[4]+2,indices[3]+2) += tmp_14_11;
    matrix(indices[4]+2,indices[4]+0) += tmp_14_12;
    matrix(indices[4]+2,indices[4]+1) += tmp_14_13;
    matrix(indices[4]+2,indices[4]+2) += tmp_14_14;
    matrix(indices[4]+2,indices[5]+0) += tmp_14_15;
    matrix(indices[4]+2,indices[5]+1) += tmp_14_16;
    matrix(indices[4]+2,indices[5]+2) += tmp_14_17;
    matrix(indices[5]+0,indices[0]+0) += tmp_15_0;
    matrix(indices[5]+0,indices[0]+1) += tmp_15_1;
    matrix(indices[5]+0,indices[0]+2) += tmp_15_2;
    matrix(indices[5]+0,indices[1]+0) += tmp_15_3;
    matrix(indices[5]+0,indices[1]+1) += tmp_15_4;
    matrix(indices[5]+0,indices[1]+2) += tmp_15_5;
    matrix(indices[5]+0,indices[2]+0) += tmp_15_6;
    matrix(indices[5]+0,indices[2]+1) += tmp_15_7;
    matrix(indices[5]+0,indices[2]+2) += tmp_15_8;
    matrix(indices[5]+0,indices[3]+0) += tmp_15_9;
    matrix(indices[5]+0,indices[3]+1) += tmp_15_10;
    matrix(indices[5]+0,indices[3]+2) += tmp_15_11;
    matrix(indices[5]+0,indices[4]+0) += tmp_15_12;
    matrix(indices[5]+0,indices[4]+1) += tmp_15_13;
    matrix(indices[5]+0,indices[4]+2) += tmp_15_14;
    matrix(indices[5]+0,indices[5]+0) += tmp_15_15;
    matrix(indices[5]+0,indices[5]+1) += tmp_15_16;
    matrix(indices[5]+0,indices[5]+2) += tmp_15_17;
    matrix(indices[5]+1,indices[0]+0) += tmp_16_0;
    matrix(indices[5]+1,indices[0]+1) += tmp_16_1;
    matrix(indices[5]+1,indices[0]+2) += tmp_16_2;
    matrix(indices[5]+1,indices[1]+0) += tmp_16_3;
    matrix(indices[5]+1,indices[1]+1) += tmp_16_4;
    matrix(indices[5]+1,indices[1]+2) += tmp_16_5;
    matrix(indices[5]+1,indices[2]+0) += tmp_16_6;
    matrix(indices[5]+1,indices[2]+1) += tmp_16_7;
    matrix(indices[5]+1,indices[2]+2) += tmp_16_8;
    matrix(indices[5]+1,indices[3]+0) += tmp_16_9;
    matrix(indices[5]+1,indices[3]+1) += tmp_16_10;
    matrix(indices[5]+1,indices[3]+2) += tmp_16_11;
    matrix(indices[5]+1,indices[4]+0) += tmp_16_12;
    matrix(indices[5]+1,indices[4]+1) += tmp_16_13;
    matrix(indices[5]+1,indices[4]+2) += tmp_16_14;
    matrix(indices[5]+1,indices[5]+0) += tmp_16_15;
    matrix(indices[5]+1,indices[5]+1) += tmp_16_16;
    matrix(indices[5]+1,indices[5]+2) += tmp_16_17;
    matrix(indices[5]+2,indices[0]+0) += tmp_17_0;
    matrix(indices[5]+2,indices[0]+1) += tmp_17_1;
    matrix(indices[5]+2,indices[0]+2) += tmp_17_2;
    matrix(indices[5]+2,indices[1]+0) += tmp_17_3;
    matrix(indices[5]+2,indices[1]+1) += tmp_17_4;
    matrix(indices[5]+2,indices[1]+2) += tmp_17_5;
    matrix(indices[5]+2,indices[2]+0) += tmp_17_6;
    matrix(indices[5]+2,indices[2]+1) += tmp_17_7;
    matrix(indices[5]+2,indices[2]+2) += tmp_17_8;
    matrix(indices[5]+2,indices[3]+0) += tmp_17_9;
    matrix(indices[5]+2,indices[3]+1) += tmp_17_10;
    matrix(indices[5]+2,indices[3]+2) += tmp_17_11;
    matrix(indices[5]+2,indices[4]+0) += tmp_17_12;
    matrix(indices[5]+2,indices[4]+1) += tmp_17_13;
    matrix(indices[5]+2,indices[4]+2) += tmp_17_14;
    matrix(indices[5]+2,indices[5]+0) += tmp_17_15;
    matrix(indices[5]+2,indices[5]+1) += tmp_17_16;
    matrix(indices[5]+2,indices[5]+2) += tmp_17_17;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1-var_inter[0]; reg0=reg0-var_inter[1]; T reg1=reg0*elem.pos(0)[1]; T reg2=var_inter[0]*elem.pos(1)[1]; T reg3=var_inter[0]*elem.pos(1)[2];
    T reg4=reg0*elem.pos(0)[2]; T reg5=1-var_inter[2]; T reg6=reg2+reg1; T reg7=var_inter[1]*elem.pos(2)[1]; T reg8=reg3+reg4;
    T reg9=var_inter[1]*elem.pos(2)[2]; T reg10=reg7+reg6; T reg11=elem.pos(3)[1]*reg0; T reg12=reg5*elem.pos(2)[2]; T reg13=reg5*elem.pos(1)[1];
    T reg14=reg5*elem.pos(0)[1]; T reg15=reg5*elem.pos(1)[2]; T reg16=reg5*elem.pos(0)[2]; T reg17=reg5*elem.pos(2)[1]; T reg18=reg8+reg9;
    T reg19=elem.pos(3)[2]*reg0; reg12=reg12-reg16; reg15=reg15-reg16; reg17=reg17-reg14; T reg20=var_inter[2]*elem.pos(3)[1];
    reg13=reg13-reg14; T reg21=var_inter[2]*elem.pos(3)[2]; T reg22=elem.pos(0)[0]*reg0; T reg23=var_inter[0]*elem.pos(1)[0]; T reg24=var_inter[0]*elem.pos(4)[2];
    reg19=reg19-reg18; reg11=reg11-reg10; T reg25=var_inter[0]*elem.pos(4)[1]; T reg26=var_inter[2]*elem.pos(5)[2]; reg17=reg17-reg20;
    reg24=reg19+reg24; reg19=var_inter[2]*elem.pos(5)[1]; reg12=reg12-reg21; T reg27=reg5*elem.pos(2)[0]; T reg28=reg5*elem.pos(1)[0];
    T reg29=reg5*elem.pos(0)[0]; T reg30=var_inter[2]*elem.pos(4)[2]; reg15=reg15-reg21; T reg31=var_inter[1]*elem.pos(5)[1]; T reg32=reg23+reg22;
    T reg33=var_inter[2]*elem.pos(4)[1]; reg13=reg13-reg20; T reg34=var_inter[1]*elem.pos(5)[2]; T reg35=1+(*f.m).poisson_ratio; reg25=reg11+reg25;
    reg11=var_inter[1]*elem.pos(2)[0]; T reg36=elem.pos(3)[0]*reg0; reg12=reg26+reg12; reg25=reg31+reg25; reg24=reg34+reg24;
    reg26=reg32+reg11; reg17=reg19+reg17; reg27=reg27-reg29; reg30=reg15+reg30; reg33=reg13+reg33;
    reg35=reg35/(*f.m).elastic_modulus; reg28=reg28-reg29; reg13=var_inter[2]*elem.pos(3)[0]; reg15=reg12*reg25; reg19=reg33*reg24;
    reg31=reg17*reg24; reg34=pow(reg35,2); T reg37=reg30*reg25; reg27=reg27-reg13; T reg38=var_inter[2]*elem.pos(5)[0];
    T reg39=var_inter[2]*elem.pos(4)[0]; reg28=reg28-reg13; reg36=reg36-reg26; T reg40=var_inter[0]*elem.pos(4)[0]; reg27=reg38+reg27;
    reg38=var_inter[1]*elem.pos(5)[0]; reg39=reg28+reg39; reg35=reg35*reg34; reg28=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg41=1.0/(*f.m).elastic_modulus;
    T reg42=reg30*reg17; T reg43=reg33*reg12; reg37=reg19-reg37; reg40=reg36+reg40; reg15=reg31-reg15;
    reg42=reg43-reg42; reg19=reg41*reg35; reg35=reg28*reg35; reg31=reg28*reg34; reg34=reg41*reg34;
    reg36=reg27*reg37; reg43=reg39*reg15; reg40=reg38+reg40; reg38=reg28*reg34; T reg44=reg28*reg31;
    reg34=reg41*reg34; T reg45=reg28*reg19; T reg46=reg28*reg35; T reg47=reg30*reg40; reg19=reg41*reg19;
    T reg48=reg33*reg40; T reg49=reg39*reg12; T reg50=reg39*reg25; T reg51=reg17*reg40; T reg52=reg39*reg24;
    reg25=reg27*reg25; reg12=reg12*reg40; reg24=reg27*reg24; reg40=reg40*reg42; reg30=reg30*reg27;
    reg36=reg43-reg36; reg38=reg44+reg38; reg34=reg34-reg44; reg31=reg41*reg31; reg27=reg33*reg27;
    reg35=reg41*reg35; reg45=reg46+reg45; reg19=reg19-reg46; reg40=reg36+reg40; reg30=reg49-reg30;
    reg17=reg39*reg17; reg47=reg52-reg47; reg12=reg24-reg12; reg48=reg50-reg48; reg51=reg25-reg51;
    reg12=reg12/reg40; reg24=reg44+reg31; reg38=reg28*reg38; reg27=reg17-reg27; reg15=reg15/reg40;
    reg34=reg41*reg34; reg51=reg51/reg40; reg35=reg46+reg35; reg37=reg37/reg40; reg30=reg30/reg40;
    reg42=reg42/reg40; reg41=reg41*reg19; reg17=reg28*reg45; reg47=reg47/reg40; reg48=reg48/reg40;
    reg25=var_inter[2]*reg47; reg33=var_inter[2]*reg12; reg36=var_inter[2]*reg37; reg39=reg5*reg37; reg43=reg5*reg15;
    reg46=var_inter[2]*reg48; reg49=var_inter[2]*reg15; reg50=reg5*reg47; reg52=var_inter[2]*reg51; reg24=reg28*reg24;
    reg38=reg34-reg38; reg34=reg5*reg12; reg27=reg27/reg40; T reg53=var_inter[0]*reg42; T reg54=var_inter[0]*reg30;
    reg28=reg28*reg35; reg17=reg41-reg17; reg41=var_inter[0]*reg27; T reg55=reg36-reg49; T reg56=reg33+reg54;
    T reg57=reg33-reg25; T reg58=reg49+reg53; T reg59=reg46-reg52; T reg60=var_inter[1]*reg42; T reg61=var_inter[1]*reg30;
    T reg62=var_inter[1]*reg27; T reg63=reg39-reg43; T reg64=reg0*reg42; T reg65=reg34-reg50; T reg66=reg0*reg30;
    T reg67=reg0*reg27; T reg68=reg5*reg51; T reg69=reg5*reg48; reg24=reg38-reg24; reg28=reg17-reg28;
    reg55=reg55+reg64; reg17=0.5*reg58; reg38=reg60-reg36; T reg70=reg25-reg61; T reg71=reg54-reg34;
    T reg72=reg69-reg68; T reg73=reg62-reg46; T reg74=reg43-reg53; reg57=reg57-reg66; T reg75=reg41+reg52;
    T reg76=0.5*reg56; reg59=reg59+reg67; T reg77=reg39+reg60; T reg78=reg61+reg50; T reg79=reg62+reg69;
    reg65=reg65+reg66; reg24=reg24/reg28; reg63=reg63-reg64; T reg80=0.5*reg63; T reg81=0.5*reg55;
    T reg82=0.5*reg77; T reg83=0.5*reg59; reg72=reg72-reg67; T reg84=reg24*reg17; T reg85=0.5*reg38;
    T reg86=0.5*reg71; T reg87=0.5*reg70; T reg88=reg24*reg76; T reg89=reg68-reg41; T reg90=0.5*reg78;
    T reg91=0.5*reg75; T reg92=0.5*reg74; T reg93=0.5*reg73; T reg94=0.5*reg65; reg19=reg19/reg28;
    T reg95=0.5*reg79; T reg96=0.5*reg57; T reg97=reg24*reg93; T reg98=0.5*reg72; T reg99=reg24*reg80;
    reg84=2*reg84; T reg100=0.5*reg89; reg35=reg35/reg28; T reg101=reg24*reg87; reg28=reg45/reg28;
    reg45=reg19*reg75; T reg102=reg24*reg90; T reg103=reg24*reg95; T reg104=reg24*reg86; T reg105=reg24*reg91;
    T reg106=reg19*reg56; T reg107=2*reg88; T reg108=reg24*reg82; T reg109=reg24*reg81; T reg110=reg24*reg92;
    T reg111=reg19*reg58; T reg112=reg24*reg94; T reg113=reg24*reg85; T reg114=reg24*reg96; T reg115=reg24*reg83;
    T reg116=reg19*reg65; reg113=2*reg113; T reg117=reg19*reg71; T reg118=reg28*reg55; T reg119=reg35*reg73;
    T reg120=reg19*reg57; reg97=2*reg97; T reg121=reg19*reg38; T reg122=reg28*reg38; T reg123=reg19*reg70;
    T reg124=reg35*reg75; T reg125=reg19*reg72; T reg126=reg19*reg89; T reg127=reg28*reg56; T reg128=reg19*reg79;
    reg105=2*reg105; T reg129=reg19*reg59; T reg130=reg35*reg59; T reg131=reg19*reg73; reg109=2*reg109;
    reg115=2*reg115; reg114=2*reg114; T reg132=reg19*reg55; T reg133=reg35*reg79; T reg134=reg24*reg98;
    reg112=2*reg112; T reg135=reg19*reg63; T reg136=2*reg108; T reg137=reg28*reg78; T reg138=2*reg103;
    reg102=2*reg102; T reg139=reg90*reg107; T reg140=reg77*reg111; T reg141=reg19*reg77; reg101=2*reg101;
    T reg142=reg28*reg58; reg110=2*reg110; reg99=2*reg99; reg104=2*reg104; T reg143=reg24*reg100;
    T reg144=reg19*reg78; T reg145=reg79*reg45; T reg146=reg28*reg77; T reg147=reg19*reg74; T reg148=reg82*reg84;
    T reg149=reg78*reg106; T reg150=reg57*reg120; T reg151=reg81*reg109; T reg152=reg57*reg106; T reg153=reg81*reg84;
    T reg154=reg35*reg57; T reg155=reg72*reg128; T reg156=reg57*reg123; T reg157=reg81*reg113; T reg158=reg59*reg129;
    T reg159=reg35*reg78; T reg160=reg80*reg138; T reg161=reg72*reg129; T reg162=reg35*reg56; T reg163=reg96*reg101;
    T reg164=reg72*reg45; T reg165=reg35*reg70; T reg166=reg55*reg121; T reg167=reg72*reg131; T reg168=reg147*reg74;
    T reg169=reg86*reg102; T reg170=reg96*reg107; T reg171=reg74*reg141; T reg172=reg74*reg133; T reg173=reg100*reg136;
    T reg174=reg86*reg114; T reg175=reg73*reg131; T reg176=reg85*reg113; T reg177=reg65*reg116; T reg178=reg70*reg123;
    T reg179=reg87*reg101; T reg180=reg99*reg80; T reg181=reg28*reg74; T reg182=reg38*reg121; T reg183=reg110*reg80;
    T reg184=reg75*reg131; T reg185=reg75*reg45; T reg186=reg56*reg123; T reg187=reg17*reg113; T reg188=reg80*reg109;
    T reg189=reg91*reg107; T reg190=reg56*reg124; T reg191=reg56*reg106; T reg192=reg65*reg106; T reg193=reg80*reg84;
    T reg194=reg17*reg84; T reg195=reg76*reg84; T reg196=reg65*reg123; T reg197=reg58*reg127; T reg198=reg80*reg113;
    T reg199=reg76*reg107; T reg200=reg35*reg71; T reg201=reg58*reg111; T reg202=reg59*reg131; T reg203=reg72*reg126;
    T reg204=reg59*reg45; T reg205=reg146*reg72; T reg206=reg79*reg142; T reg207=reg79*reg129; T reg208=reg89*reg131;
    T reg209=reg82*reg115; T reg210=reg90*reg102; T reg211=reg77*reg141; T reg212=reg90*reg136; T reg213=reg77*reg137;
    T reg214=reg79*reg118; T reg215=reg79*reg128; T reg216=reg90*reg114; T reg217=reg77*reg132; T reg218=reg82*reg113;
    T reg219=reg78*reg123; T reg220=reg77*reg130; T reg221=reg95*reg109; T reg222=reg149+reg148; reg140=reg139+reg140;
    T reg223=reg95*reg105; T reg224=reg77*reg124; T reg225=reg95*reg84; T reg226=reg90*reg101; T reg227=reg77*reg121;
    T reg228=reg77*reg119; T reg229=reg82*reg109; T reg230=reg78*reg120; T reg231=reg95*reg102; T reg232=reg78*reg133;
    T reg233=reg82*reg136; T reg234=reg144*reg78; T reg235=reg95*reg113; T reg236=reg74*reg132; T reg237=reg86*reg107;
    T reg238=reg55*reg111; T reg239=reg96*reg114; T reg240=reg74*reg111; T reg241=reg86*reg101; T reg242=reg55*reg132;
    T reg243=reg74*reg121; reg131=reg79*reg131; T reg244=reg71*reg117; T reg245=reg110*reg92; T reg246=reg144*reg71;
    T reg247=reg92*reg136; T reg248=reg71*reg120; T reg249=reg92*reg109; T reg250=reg82*reg97; T reg251=reg71*reg106;
    T reg252=reg92*reg84; reg123=reg71*reg123; T reg253=reg92*reg113; T reg254=reg79*reg122; T reg255=reg89*reg126;
    T reg256=reg146*reg89; reg145=reg148+reg145; reg148=reg92*reg138; T reg257=reg82*reg105; T reg258=reg89*reg128;
    T reg259=reg89*reg129; T reg260=reg58*reg121; T reg261=reg76*reg101; T reg262=reg89*reg45; T reg263=reg65*reg120;
    T reg264=reg28*reg57; T reg265=reg63*reg111; T reg266=reg98*reg136; T reg267=reg63*reg133; T reg268=reg104*reg86;
    T reg269=reg147*reg63; T reg270=reg104*reg94; T reg271=reg80*reg136; reg143=2*reg143; T reg272=reg94*reg101;
    T reg273=reg28*reg65; reg134=2*reg134; reg121=reg63*reg121; T reg274=reg144*reg65; T reg275=reg28*reg71;
    T reg276=reg65*reg117; T reg277=reg94*reg112; T reg278=reg63*reg135; T reg279=reg28*reg70; T reg280=reg94*reg102;
    T reg281=reg35*reg72; T reg282=reg63*reg132; T reg283=reg94*reg114; T reg284=reg35*reg89; T reg285=reg94*reg107;
    T reg286=reg72*reg125; T reg287=reg63*reg141; T reg288=reg190+reg189; T reg289=reg91*reg84; T reg290=reg89*reg118;
    T reg291=reg89*reg142; T reg292=reg92*reg115; reg217=reg216-reg217; T reg293=reg76*reg113; reg259=reg249+reg259;
    reg274=reg274-reg271; T reg294=reg95*reg115; reg221=reg220+reg221; T reg295=reg86*reg115; T reg296=reg90*reg109;
    T reg297=reg77*reg264; T reg298=reg89*reg154; T reg299=reg65*reg284; T reg300=reg98*reg138; T reg301=reg91*reg113;
    reg213=reg212+reg213; T reg302=reg95*reg138; T reg303=reg58*reg119; T reg304=reg210+reg211; reg208=reg253+reg208;
    T reg305=reg63*reg137; T reg306=reg86*reg97; T reg307=reg89*reg165; T reg308=reg92*reg97; T reg309=reg89*reg122;
    T reg310=reg77*reg133; T reg311=reg95*reg136; T reg312=reg280-reg287; reg262=reg252+reg262; T reg313=reg94*reg136;
    T reg314=reg194+reg191; T reg315=reg86*reg105; T reg316=reg89*reg162; T reg317=reg91*reg97; reg260=reg260-reg261;
    T reg318=reg71*reg133; reg246=reg246-reg247; T reg319=reg267+reg266; T reg320=reg56*reg119; T reg321=reg92*reg102;
    T reg322=reg146*reg71; T reg323=reg91*reg101; T reg324=reg104*reg100; T reg325=reg71*reg284; reg244=reg244+reg245;
    T reg326=reg100*reg113; T reg327=reg74*reg119; T reg328=reg74*reg279; T reg329=reg86*reg113; reg185=reg194+reg185;
    reg194=reg100*reg97; reg243=reg241+reg243; T reg330=reg100*reg84; T reg331=reg74*reg124; T reg332=reg74*reg127;
    T reg333=reg86*reg84; T reg334=reg100*reg105; reg240=reg240-reg237; reg282=reg282+reg283; T reg335=reg98*reg115;
    T reg336=reg247+reg258; T reg337=reg86*reg138; T reg338=reg89*reg159; T reg339=reg256+reg148; reg255=reg245+reg255;
    reg278=reg278+reg277; reg245=reg17*reg101; T reg340=reg100*reg101; T reg341=reg71*reg119; T reg342=reg56*reg122;
    reg253=reg123+reg253; reg123=reg98*reg134; T reg343=reg92*reg101; T reg344=reg71*reg122; T reg345=reg100*reg107;
    T reg346=reg71*reg124; reg252=reg252-reg251; T reg347=reg92*reg107; T reg348=reg71*reg142; T reg349=reg100*reg114;
    T reg350=reg71*reg130; reg186=reg187-reg186; reg249=reg248+reg249; reg248=reg92*reg114; T reg351=reg71*reg118;
    T reg352=reg100*reg102; T reg353=reg99*reg94; T reg354=reg83*reg84; T reg355=reg55*reg124; reg204=reg153+reg204;
    T reg356=reg96*reg84; T reg357=reg55*reg127; T reg358=reg83*reg105; reg238=reg238-reg170; T reg359=reg83*reg109;
    T reg360=reg55*reg130; T reg361=reg96*reg109; T reg362=reg55*reg264; T reg363=reg59*reg122; T reg364=reg83*reg115;
    reg242=reg242+reg239; reg269=reg269+reg270; T reg365=reg81*reg97; reg131=reg218+reg131; T reg366=reg143*reg98;
    T reg367=reg79*reg165; T reg368=reg90*reg97; reg250=reg254+reg250; reg145=reg139+reg145; reg254=reg59*reg165;
    T reg369=reg96*reg97; T reg370=reg79*reg162; reg158=reg151+reg158; T reg371=reg83*reg101; T reg372=reg57*reg119;
    reg156=reg156+reg157; T reg373=reg59*reg142; T reg374=reg81*reg101; T reg375=reg57*reg122; T reg376=reg81*reg105;
    T reg377=reg83*reg107; T reg378=reg57*reg124; reg153=reg153-reg152; T reg379=reg281*reg63; T reg380=reg81*reg107;
    T reg381=reg57*reg142; T reg382=reg83*reg114; T reg383=reg57*reg130; T reg384=reg59*reg162; reg151=reg150+reg151;
    reg150=reg99*reg98; T reg385=reg96*reg105; T reg386=reg83*reg113; T reg387=reg55*reg119; T reg388=reg96*reg113;
    T reg389=reg55*reg279; T reg390=reg83*reg97; reg166=reg166+reg163; T reg391=reg65*reg133; T reg392=reg271+reg155;
    T reg393=reg82*reg114; T reg394=reg78*reg118; T reg395=reg146*reg65; reg231=reg232+reg231; reg276=reg276+reg183;
    reg161=reg188+reg161; T reg396=reg80*reg102; reg234=reg234+reg233; reg167=reg198+reg167; reg195=reg197+reg195;
    reg235=reg228+reg235; T reg397=reg77*reg279; T reg398=reg90*reg113; T reg399=reg95*reg97; reg227=reg226-reg227;
    T reg400=reg63*reg284; T reg401=reg58*reg124; reg225=reg224+reg225; T reg402=reg110*reg98; T reg403=reg77*reg127;
    T reg404=reg90*reg84; T reg405=reg140+reg223; T reg406=reg58*reg279; reg286=reg286+reg180; T reg407=reg90*reg105;
    reg257=reg206+reg257; reg206=reg275*reg63; T reg408=reg110*reg94; reg207=reg229+reg207; T reg409=reg79*reg154;
    T reg410=reg90*reg115; T reg411=reg273*reg63; reg209=reg214+reg209; reg202=reg157+reg202; reg157=reg233+reg215;
    reg214=reg95*reg101; T reg412=reg78*reg119; reg218=reg219-reg218; reg219=reg82*reg101; T reg413=reg78*reg122;
    T reg414=reg95*reg107; T reg415=reg78*reg124; reg223=reg223+reg222; reg201=reg201+reg199; T reg416=reg104*reg98;
    T reg417=reg91*reg105; T reg418=reg82*reg107; T reg419=reg78*reg142; T reg420=reg95*reg114; T reg421=reg78*reg130;
    reg229=reg230-reg229; reg230=reg94*reg97; T reg422=reg92*reg105; T reg423=reg98*reg101; T reg424=reg65*reg119;
    reg184=reg187+reg184; reg198=reg196+reg198; reg187=reg63*reg124; reg196=reg98*reg84; T reg425=reg98*reg102;
    reg168=reg268+reg168; T reg426=reg143*reg100; T reg427=reg80*reg101; T reg428=reg65*reg122; reg178=reg178+reg176;
    T reg429=reg110*reg86; T reg430=reg275*reg74; T reg431=reg98*reg107; T reg432=reg65*reg124; T reg433=reg193-reg192;
    T reg434=reg74*reg284; T reg435=reg110*reg100; T reg436=reg80*reg107; T reg437=reg65*reg142; T reg438=reg70*reg119;
    T reg439=reg98*reg114; T reg440=reg94*reg109; T reg441=reg63*reg264; T reg442=reg169-reg171; T reg443=reg65*reg130;
    T reg444=reg93*reg101; T reg445=reg100*reg138; T reg446=reg94*reg115; T reg447=reg72*reg154; T reg448=reg80*reg115;
    T reg449=reg72*reg118; T reg450=reg93*reg97; reg265=reg265-reg285; T reg451=reg98*reg105; T reg452=reg38*reg279;
    reg182=reg182+reg179; T reg453=reg80*reg105; T reg454=reg94*reg138; T reg455=reg72*reg159; T reg456=reg87*reg113;
    T reg457=reg72*reg162; T reg458=reg94*reg105; T reg459=reg205+reg160; reg203=reg183+reg203; reg183=reg98*reg109;
    T reg460=reg63*reg127; T reg461=reg94*reg84; reg164=reg193+reg164; reg193=reg143*reg94; T reg462=reg72*reg200;
    T reg463=reg38*reg119; T reg464=reg143*reg80; T reg465=reg72*reg122; T reg466=reg80*reg97; T reg467=reg72*reg181;
    T reg468=reg93*reg113; T reg469=reg63*reg130; T reg470=reg72*reg165; T reg471=reg65*reg118; T reg472=reg75*reg165;
    T reg473=reg172+reg173; T reg474=reg104*reg80; T reg475=reg65*reg181; T reg476=reg72*reg142; reg175=reg176+reg175;
    reg176=reg98*reg112; T reg477=reg281*reg65; reg236=reg174+reg236; T reg478=reg100*reg115; reg180=reg177+reg180;
    reg177=reg86*reg109; reg279=reg63*reg279; T reg479=reg94*reg113; T reg480=reg74*reg264; reg113=reg98*reg113;
    T reg481=reg17*reg97; reg119=reg63*reg119; T reg482=reg74*reg130; T reg483=reg100*reg109; T reg484=reg75*reg122;
    reg188=reg263+reg188; reg121=reg121+reg272; reg263=reg76*reg97; T reg485=reg98*reg97; T reg486=reg86*reg136;
    T reg487=reg74*reg137; T reg488=reg80*reg114; reg402=reg400+reg402; reg400=reg40*reg459; T reg489=reg40*reg225;
    reg378=reg378-reg377; reg407=reg407+reg370; reg178=reg450+reg178; reg456=reg452+reg456; reg227=reg227-reg399;
    reg203=reg270+reg203; reg444=reg438+reg444; reg397=reg398-reg397; reg376=reg373+reg376; reg270=reg40*reg195;
    reg439=reg443+reg439; reg373=reg40*reg235; reg153=reg358+reg153; reg193=reg462+reg193; reg488=reg471+reg488;
    reg234=reg302+reg234; reg180=reg123+reg180; reg381=reg381-reg380; reg461=reg461-reg460; reg175=reg179+reg175;
    reg179=reg40*reg231; reg188=reg335+reg188; reg398=reg40*reg213; reg353=reg411+reg353; reg411=reg310+reg311;
    reg448=reg449+reg448; reg158=reg239+reg158; reg361=reg362+reg361; reg217=reg217-reg294; reg280=reg280-reg392;
    reg297=reg296-reg297; reg293=reg406-reg293; reg371=reg372+reg371; reg277=reg286+reg277; reg242=reg242+reg364;
    reg239=reg40*reg221; reg156=reg390+reg156; reg455=reg455-reg454; reg433=reg451+reg433; reg359=reg360+reg359;
    reg113=reg119+reg113; reg374=reg375+reg374; reg119=reg40*reg405; reg451=reg265+reg451; reg479=reg279+reg479;
    reg404=reg404+reg403; reg437=reg437-reg436; reg399=reg218-reg399; reg474=reg475+reg474; reg385=reg385-reg384;
    reg214=reg412-reg214; reg388=reg389+reg388; reg202=reg163+reg202; reg210=reg210+reg157; reg427=reg428+reg427;
    reg390=reg166+reg390; reg274=reg274-reg300; reg163=reg40*reg250; reg166=reg40*reg209; reg425=reg425-reg391;
    reg196=reg187+reg196; reg354=reg355+reg354; reg409=reg410-reg409; reg408=reg206+reg408; reg121=reg121+reg485;
    reg269=reg269+reg366; reg207=reg216-reg207; reg432=reg432-reg431; reg356=reg356-reg357; reg187=reg40*reg145;
    reg464=reg467+reg464; reg131=reg226-reg131; reg393=reg394-reg393; reg369=reg254+reg369; reg206=reg40*reg257;
    reg382=reg383+reg382; reg365=reg363+reg365; reg294=reg229-reg294; reg468=reg463+reg468; reg420=reg421-reg420;
    reg176=reg477+reg176; reg201=reg201+reg417; reg291=reg422+reg291; reg151=reg364+reg151; reg367=reg368-reg367;
    reg358=reg238+reg358; reg419=reg419+reg418; reg423=reg424+reg423; reg216=reg40*reg223; reg150=reg379+reg150;
    reg204=reg204-reg170; reg415=reg415+reg414; reg198=reg485+reg198; reg219=reg413-reg219; reg386=reg387+reg386;
    reg416=reg299+reg416; reg323=reg323-reg320; reg164=reg164-reg285; reg401=reg289+reg401; reg263=reg472-reg263;
    reg244=reg426+reg244; reg314=reg417+reg314; reg435=reg434+reg435; reg255=reg268+reg255; reg305=reg305-reg313;
    reg272=reg167+reg272; reg349=reg350+reg349; reg183=reg469+reg183; reg260=reg260+reg317; reg326=reg327+reg326;
    reg458=reg458-reg457; reg276=reg366+reg276; reg348=reg348-reg347; reg315=reg315-reg316; reg340=reg341+reg340;
    reg328=reg329+reg328; reg335=reg282+reg335; reg352=reg352-reg318; reg338=reg338-reg337; reg230=reg470+reg230;
    reg167=reg40*reg288; reg169=reg169-reg336; reg246=reg246-reg445; reg218=reg40*reg319; reg487=reg487-reg486;
    reg321=reg321-reg322; reg442=reg442-reg445; reg292=reg290+reg292; reg248=reg351+reg248; reg466=reg465+reg466;
    reg295=reg298+reg295; reg226=reg40*reg339; reg324=reg325+reg324; reg261=reg184-reg261; reg186=reg317+reg186;
    reg184=reg40*reg473; reg249=reg478+reg249; reg259=reg174+reg259; reg330=reg331+reg330; reg480=reg177+reg480;
    reg342=reg245-reg342; reg283=reg161+reg283; reg208=reg241+reg208; reg301=reg303+reg301; reg343=reg344+reg343;
    reg333=reg333-reg332; reg123=reg278+reg123; reg450=reg182+reg450; reg304=reg304+reg302; reg252=reg334+reg252;
    reg346=reg346-reg345; reg334=reg240+reg334; reg430=reg429+reg430; reg483=reg482+reg483; reg312=reg312-reg300;
    reg446=reg447+reg446; reg262=reg262-reg237; reg306=reg307+reg306; reg396=reg396-reg395; reg478=reg236+reg478;
    reg185=reg199+reg185; reg440=reg441+reg440; reg243=reg243+reg194; reg253=reg194+reg253; reg308=reg309+reg308;
    reg453=reg476+reg453; reg426=reg168+reg426; reg481=reg484+reg481; reg188=reg40*reg188; reg367=reg40*reg367;
    reg186=reg40*reg186; reg433=reg40*reg433; reg361=reg40*reg361; reg161=ponderation*reg163; reg365=reg40*reg365;
    reg352=reg40*reg352; reg437=reg40*reg437; reg252=reg40*reg252; reg407=reg40*reg407; reg263=reg40*reg263;
    reg242=reg40*reg242; reg248=reg40*reg248; reg349=reg40*reg349; reg196=reg40*reg196; reg178=reg40*reg178;
    reg123=reg40*reg123; reg249=reg40*reg249; reg348=reg40*reg348; reg439=reg40*reg439; reg168=ponderation*reg187;
    reg131=reg40*reg131; reg435=reg40*reg435; reg269=reg40*reg269; reg335=reg40*reg335; reg151=reg40*reg151;
    reg176=reg40*reg176; reg328=reg40*reg328; reg382=reg40*reg382; reg243=reg40*reg243; reg381=reg40*reg381;
    reg376=reg40*reg376; reg180=reg40*reg180; reg478=reg40*reg478; reg153=reg40*reg153; reg330=reg40*reg330;
    reg378=reg40*reg378; reg185=reg40*reg185; reg175=reg40*reg175; reg113=reg40*reg113; reg333=reg40*reg333;
    reg374=reg40*reg374; reg353=reg40*reg353; reg480=reg40*reg480; reg156=reg40*reg156; reg158=reg40*reg158;
    reg334=reg40*reg334; reg479=reg40*reg479; reg371=reg40*reg371; reg483=reg40*reg483; reg442=reg40*reg442;
    reg359=reg40*reg359; reg204=reg40*reg204; reg246=reg40*reg246; reg488=reg40*reg488; reg358=reg40*reg358;
    reg321=reg40*reg321; reg487=reg40*reg487; reg356=reg40*reg356; reg444=reg40*reg444; reg425=reg40*reg425;
    reg324=reg40*reg324; reg354=reg40*reg354; reg150=reg40*reg150; reg121=reg40*reg121; reg390=reg40*reg390;
    reg385=reg40*reg385; reg474=reg40*reg474; reg244=reg40*reg244; reg388=reg40*reg388; reg174=ponderation*reg184;
    reg177=ponderation*reg218; reg323=reg40*reg323; reg386=reg40*reg386; reg326=reg40*reg326; reg481=reg40*reg481;
    reg283=reg40*reg283; reg193=reg40*reg193; reg234=reg40*reg234; reg314=reg40*reg314; reg259=reg40*reg259;
    reg164=reg40*reg164; reg208=reg40*reg208; reg297=reg40*reg297; reg182=ponderation*reg179; reg464=reg40*reg464;
    reg301=reg40*reg301; reg280=reg40*reg280; reg295=reg40*reg295; reg393=reg40*reg393; reg293=reg40*reg293;
    reg292=reg40*reg292; reg277=reg40*reg277; reg294=reg40*reg294; reg201=reg40*reg201; reg466=reg40*reg466;
    reg291=reg40*reg291; reg261=reg40*reg261; reg420=reg40*reg420; reg169=reg40*reg169; reg194=ponderation*reg119;
    reg308=reg40*reg308; reg453=reg40*reg453; reg402=reg40*reg402; reg404=reg40*reg404; reg260=reg40*reg260;
    reg262=reg40*reg262; reg229=ponderation*reg400; reg183=reg40*reg183; reg455=reg40*reg455; reg236=ponderation*reg489;
    reg315=reg40*reg315; reg451=reg40*reg451; reg306=reg40*reg306; reg456=reg40*reg456; reg203=reg40*reg203;
    reg227=reg40*reg227; reg238=ponderation*reg239; reg458=reg40*reg458; reg397=reg40*reg397; reg240=ponderation*reg270;
    reg276=reg40*reg276; reg305=reg40*reg305; reg241=ponderation*reg373; reg401=reg40*reg401; reg272=reg40*reg272;
    reg214=reg40*reg214; reg340=reg40*reg340; reg396=reg40*reg396; reg427=reg40*reg427; reg210=reg40*reg210;
    reg253=reg40*reg253; reg408=reg40*reg408; reg245=ponderation*reg166; reg411=reg40*reg411; reg446=reg40*reg446;
    reg426=reg40*reg426; reg448=reg40*reg448; reg450=reg40*reg450; reg343=reg40*reg343; reg409=reg40*reg409;
    reg254=ponderation*reg398; reg432=reg40*reg432; reg342=reg40*reg342; reg207=reg40*reg207; reg369=reg40*reg369;
    reg346=reg40*reg346; reg430=reg40*reg430; reg265=ponderation*reg206; reg338=reg40*reg338; reg419=reg40*reg419;
    reg217=reg40*reg217; reg461=reg40*reg461; reg468=reg40*reg468; reg274=reg40*reg274; reg423=reg40*reg423;
    reg268=ponderation*reg216; reg312=reg40*reg312; reg230=reg40*reg230; reg278=ponderation*reg167; reg279=ponderation*reg226;
    reg415=reg40*reg415; reg198=reg40*reg198; reg219=reg40*reg219; reg304=reg40*reg304; reg416=reg40*reg416;
    reg440=reg40*reg440; reg255=reg40*reg255; reg399=reg40*reg399; reg202=reg40*reg202; T tmp_15_15=ponderation*reg450;
    T tmp_14_14=ponderation*reg185; T tmp_11_11=ponderation*reg158; T tmp_12_17=ponderation*reg301; T tmp_12_16=ponderation*reg293; T tmp_17_17=ponderation*reg175;
    T tmp_13_15=ponderation*reg342; T tmp_11_16=ponderation*reg369; T tmp_14_16=ponderation*reg263; T tmp_11_15=ponderation*reg365; T tmp_1_5=ponderation*reg416;
    T tmp_16_16=ponderation*reg178; T tmp_13_14=-reg278; T tmp_11_17=ponderation*reg202; T tmp_15_17=ponderation*reg468; T tmp_0_0=ponderation*reg123;
    T tmp_11_14=ponderation*reg204; T tmp_13_16=ponderation*reg186; T tmp_14_15=ponderation*reg481; T tmp_12_12=ponderation*reg201; T tmp_14_17=ponderation*reg261;
    T tmp_11_13=ponderation*reg385; T tmp_16_17=ponderation*reg444; T tmp_13_13=ponderation*reg314; T tmp_15_16=ponderation*reg456; T tmp_13_17=ponderation*reg323;
    T tmp_12_13=-reg240; T tmp_1_6=ponderation*reg396; T tmp_11_12=ponderation*reg376; T tmp_1_4=ponderation*reg276; T tmp_12_15=ponderation*reg260;
    T tmp_3_15=ponderation*reg243; T tmp_3_14=ponderation*reg330; T tmp_3_13=ponderation*reg333; T tmp_3_12=ponderation*reg334; T tmp_3_11=ponderation*reg483;
    T tmp_3_10=ponderation*reg480; T tmp_3_9=ponderation*reg478; T tmp_3_8=-reg174; T tmp_0_9=ponderation*reg335; T tmp_3_7=ponderation*reg487;
    T tmp_3_6=ponderation*reg442; T tmp_3_5=ponderation*reg435; T tmp_3_4=ponderation*reg430; T tmp_3_3=ponderation*reg426; T tmp_2_17=ponderation*reg272;
    T tmp_0_10=ponderation*reg440; T tmp_2_16=ponderation*reg230; T tmp_5_5=ponderation*reg255; T tmp_4_17=ponderation*reg340; T tmp_4_16=ponderation*reg253;
    T tmp_4_15=ponderation*reg343; T tmp_4_14=ponderation*reg346; T tmp_4_13=ponderation*reg252; T tmp_4_12=ponderation*reg348; T tmp_4_11=ponderation*reg349;
    T tmp_4_10=ponderation*reg249; T tmp_4_9=ponderation*reg248; T tmp_4_8=ponderation*reg352; T tmp_4_7=ponderation*reg246; T tmp_4_6=ponderation*reg321;
    T tmp_4_5=ponderation*reg324; T tmp_4_4=ponderation*reg244; T tmp_0_8=-reg177; T tmp_3_17=ponderation*reg326; T tmp_3_16=ponderation*reg328;
    T tmp_1_16=ponderation*reg198; T tmp_0_13=ponderation*reg461; T tmp_1_15=ponderation*reg427; T tmp_1_14=ponderation*reg432; T tmp_1_13=ponderation*reg433;
    T tmp_0_14=ponderation*reg196; T tmp_1_12=ponderation*reg437; T tmp_1_11=ponderation*reg439; T tmp_1_10=ponderation*reg188; T tmp_1_9=ponderation*reg488;
    T tmp_1_8=ponderation*reg425; T tmp_0_15=ponderation*reg121; T tmp_1_3=ponderation*reg474; T tmp_1_2=ponderation*reg176; T tmp_1_1=ponderation*reg180;
    T tmp_0_17=ponderation*reg113; T tmp_0_16=ponderation*reg479; T tmp_2_15=ponderation*reg466; T tmp_2_14=ponderation*reg164; T tmp_2_13=ponderation*reg458;
    T tmp_2_12=ponderation*reg453; T tmp_2_11=ponderation*reg283; T tmp_0_11=ponderation*reg183; T tmp_2_10=ponderation*reg446; T tmp_2_9=ponderation*reg448;
    T tmp_2_8=ponderation*reg280; T tmp_2_7=ponderation*reg455; T tmp_2_6=-reg229; T tmp_0_12=ponderation*reg451; T tmp_2_5=ponderation*reg203;
    T tmp_2_4=ponderation*reg193; T tmp_2_3=ponderation*reg464; T tmp_2_2=ponderation*reg277; T tmp_5_12=ponderation*reg291; T tmp_1_17=ponderation*reg423;
    T tmp_9_9=ponderation*reg242; T tmp_8_17=ponderation*reg131; T tmp_8_16=ponderation*reg367; T tmp_8_15=-reg161; T tmp_8_14=-reg168;
    T tmp_0_3=ponderation*reg269; T tmp_8_13=ponderation*reg407; T tmp_8_12=-reg265; T tmp_8_11=ponderation*reg207; T tmp_8_10=ponderation*reg409;
    T tmp_8_9=-reg245; T tmp_8_8=ponderation*reg210; T tmp_7_17=ponderation*reg214; T tmp_7_16=ponderation*reg399; T tmp_7_15=ponderation*reg219;
    T tmp_7_14=ponderation*reg415; T tmp_7_13=-reg268; T tmp_10_17=ponderation*reg371; T tmp_10_16=ponderation*reg156; T tmp_0_1=ponderation*reg353;
    T tmp_10_15=ponderation*reg374; T tmp_10_14=ponderation*reg378; T tmp_10_13=ponderation*reg153; T tmp_10_12=ponderation*reg381; T tmp_10_11=ponderation*reg382;
    T tmp_10_10=ponderation*reg151; T tmp_9_17=ponderation*reg386; T tmp_9_16=ponderation*reg388; T tmp_9_15=ponderation*reg390; T tmp_0_2=ponderation*reg150;
    T tmp_9_14=ponderation*reg354; T tmp_9_13=ponderation*reg356; T tmp_9_12=ponderation*reg358; T tmp_9_11=ponderation*reg359; T tmp_9_10=ponderation*reg361;
    T tmp_6_8=ponderation*reg411; T tmp_6_7=-reg254; T tmp_6_6=ponderation*reg304; T tmp_0_6=ponderation*reg312; T tmp_5_17=ponderation*reg208;
    T tmp_5_16=ponderation*reg306; T tmp_5_15=ponderation*reg308; T tmp_5_14=ponderation*reg262; T tmp_5_13=ponderation*reg315; T tmp_12_14=ponderation*reg401;
    T tmp_5_11=ponderation*reg259; T tmp_0_7=ponderation*reg305; T tmp_5_10=ponderation*reg295; T tmp_5_9=ponderation*reg292; T tmp_5_8=ponderation*reg169;
    T tmp_5_7=ponderation*reg338; T tmp_5_6=-reg279; T tmp_1_7=ponderation*reg274; T tmp_7_12=ponderation*reg419; T tmp_7_11=ponderation*reg420;
    T tmp_7_10=ponderation*reg294; T tmp_7_9=ponderation*reg393; T tmp_7_8=-reg182; T tmp_7_7=ponderation*reg234; T tmp_6_17=-reg241;
    T tmp_0_4=ponderation*reg408; T tmp_6_16=ponderation*reg397; T tmp_6_15=ponderation*reg227; T tmp_6_14=-reg236; T tmp_6_13=ponderation*reg404;
    T tmp_6_12=-reg194; T tmp_0_5=ponderation*reg402; T tmp_6_11=-reg238; T tmp_6_10=ponderation*reg297; T tmp_6_9=ponderation*reg217;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[0]+2) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_3;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_4;
    matrix(indices[0]+0,indices[1]+2) += tmp_0_5;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_7;
    matrix(indices[0]+0,indices[2]+2) += tmp_0_8;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_9;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_10;
    matrix(indices[0]+0,indices[3]+2) += tmp_0_11;
    matrix(indices[0]+0,indices[4]+0) += tmp_0_12;
    matrix(indices[0]+0,indices[4]+1) += tmp_0_13;
    matrix(indices[0]+0,indices[4]+2) += tmp_0_14;
    matrix(indices[0]+0,indices[5]+0) += tmp_0_15;
    matrix(indices[0]+0,indices[5]+1) += tmp_0_16;
    matrix(indices[0]+0,indices[5]+2) += tmp_0_17;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[0]+2) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_3;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_4;
    matrix(indices[0]+1,indices[1]+2) += tmp_1_5;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_7;
    matrix(indices[0]+1,indices[2]+2) += tmp_1_8;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_9;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_10;
    matrix(indices[0]+1,indices[3]+2) += tmp_1_11;
    matrix(indices[0]+1,indices[4]+0) += tmp_1_12;
    matrix(indices[0]+1,indices[4]+1) += tmp_1_13;
    matrix(indices[0]+1,indices[4]+2) += tmp_1_14;
    matrix(indices[0]+1,indices[5]+0) += tmp_1_15;
    matrix(indices[0]+1,indices[5]+1) += tmp_1_16;
    matrix(indices[0]+1,indices[5]+2) += tmp_1_17;
    matrix(indices[0]+2,indices[0]+2) += tmp_2_2;
    matrix(indices[0]+2,indices[1]+0) += tmp_2_3;
    matrix(indices[0]+2,indices[1]+1) += tmp_2_4;
    matrix(indices[0]+2,indices[1]+2) += tmp_2_5;
    matrix(indices[0]+2,indices[2]+0) += tmp_2_6;
    matrix(indices[0]+2,indices[2]+1) += tmp_2_7;
    matrix(indices[0]+2,indices[2]+2) += tmp_2_8;
    matrix(indices[0]+2,indices[3]+0) += tmp_2_9;
    matrix(indices[0]+2,indices[3]+1) += tmp_2_10;
    matrix(indices[0]+2,indices[3]+2) += tmp_2_11;
    matrix(indices[0]+2,indices[4]+0) += tmp_2_12;
    matrix(indices[0]+2,indices[4]+1) += tmp_2_13;
    matrix(indices[0]+2,indices[4]+2) += tmp_2_14;
    matrix(indices[0]+2,indices[5]+0) += tmp_2_15;
    matrix(indices[0]+2,indices[5]+1) += tmp_2_16;
    matrix(indices[0]+2,indices[5]+2) += tmp_2_17;
    matrix(indices[1]+0,indices[1]+0) += tmp_3_3;
    matrix(indices[1]+0,indices[1]+1) += tmp_3_4;
    matrix(indices[1]+0,indices[1]+2) += tmp_3_5;
    matrix(indices[1]+0,indices[2]+0) += tmp_3_6;
    matrix(indices[1]+0,indices[2]+1) += tmp_3_7;
    matrix(indices[1]+0,indices[2]+2) += tmp_3_8;
    matrix(indices[1]+0,indices[3]+0) += tmp_3_9;
    matrix(indices[1]+0,indices[3]+1) += tmp_3_10;
    matrix(indices[1]+0,indices[3]+2) += tmp_3_11;
    matrix(indices[1]+0,indices[4]+0) += tmp_3_12;
    matrix(indices[1]+0,indices[4]+1) += tmp_3_13;
    matrix(indices[1]+0,indices[4]+2) += tmp_3_14;
    matrix(indices[1]+0,indices[5]+0) += tmp_3_15;
    matrix(indices[1]+0,indices[5]+1) += tmp_3_16;
    matrix(indices[1]+0,indices[5]+2) += tmp_3_17;
    matrix(indices[1]+1,indices[1]+1) += tmp_4_4;
    matrix(indices[1]+1,indices[1]+2) += tmp_4_5;
    matrix(indices[1]+1,indices[2]+0) += tmp_4_6;
    matrix(indices[1]+1,indices[2]+1) += tmp_4_7;
    matrix(indices[1]+1,indices[2]+2) += tmp_4_8;
    matrix(indices[1]+1,indices[3]+0) += tmp_4_9;
    matrix(indices[1]+1,indices[3]+1) += tmp_4_10;
    matrix(indices[1]+1,indices[3]+2) += tmp_4_11;
    matrix(indices[1]+1,indices[4]+0) += tmp_4_12;
    matrix(indices[1]+1,indices[4]+1) += tmp_4_13;
    matrix(indices[1]+1,indices[4]+2) += tmp_4_14;
    matrix(indices[1]+1,indices[5]+0) += tmp_4_15;
    matrix(indices[1]+1,indices[5]+1) += tmp_4_16;
    matrix(indices[1]+1,indices[5]+2) += tmp_4_17;
    matrix(indices[1]+2,indices[1]+2) += tmp_5_5;
    matrix(indices[1]+2,indices[2]+0) += tmp_5_6;
    matrix(indices[1]+2,indices[2]+1) += tmp_5_7;
    matrix(indices[1]+2,indices[2]+2) += tmp_5_8;
    matrix(indices[1]+2,indices[3]+0) += tmp_5_9;
    matrix(indices[1]+2,indices[3]+1) += tmp_5_10;
    matrix(indices[1]+2,indices[3]+2) += tmp_5_11;
    matrix(indices[1]+2,indices[4]+0) += tmp_5_12;
    matrix(indices[1]+2,indices[4]+1) += tmp_5_13;
    matrix(indices[1]+2,indices[4]+2) += tmp_5_14;
    matrix(indices[1]+2,indices[5]+0) += tmp_5_15;
    matrix(indices[1]+2,indices[5]+1) += tmp_5_16;
    matrix(indices[1]+2,indices[5]+2) += tmp_5_17;
    matrix(indices[2]+0,indices[2]+0) += tmp_6_6;
    matrix(indices[2]+0,indices[2]+1) += tmp_6_7;
    matrix(indices[2]+0,indices[2]+2) += tmp_6_8;
    matrix(indices[2]+0,indices[3]+0) += tmp_6_9;
    matrix(indices[2]+0,indices[3]+1) += tmp_6_10;
    matrix(indices[2]+0,indices[3]+2) += tmp_6_11;
    matrix(indices[2]+0,indices[4]+0) += tmp_6_12;
    matrix(indices[2]+0,indices[4]+1) += tmp_6_13;
    matrix(indices[2]+0,indices[4]+2) += tmp_6_14;
    matrix(indices[2]+0,indices[5]+0) += tmp_6_15;
    matrix(indices[2]+0,indices[5]+1) += tmp_6_16;
    matrix(indices[2]+0,indices[5]+2) += tmp_6_17;
    matrix(indices[2]+1,indices[2]+1) += tmp_7_7;
    matrix(indices[2]+1,indices[2]+2) += tmp_7_8;
    matrix(indices[2]+1,indices[3]+0) += tmp_7_9;
    matrix(indices[2]+1,indices[3]+1) += tmp_7_10;
    matrix(indices[2]+1,indices[3]+2) += tmp_7_11;
    matrix(indices[2]+1,indices[4]+0) += tmp_7_12;
    matrix(indices[2]+1,indices[4]+1) += tmp_7_13;
    matrix(indices[2]+1,indices[4]+2) += tmp_7_14;
    matrix(indices[2]+1,indices[5]+0) += tmp_7_15;
    matrix(indices[2]+1,indices[5]+1) += tmp_7_16;
    matrix(indices[2]+1,indices[5]+2) += tmp_7_17;
    matrix(indices[2]+2,indices[2]+2) += tmp_8_8;
    matrix(indices[2]+2,indices[3]+0) += tmp_8_9;
    matrix(indices[2]+2,indices[3]+1) += tmp_8_10;
    matrix(indices[2]+2,indices[3]+2) += tmp_8_11;
    matrix(indices[2]+2,indices[4]+0) += tmp_8_12;
    matrix(indices[2]+2,indices[4]+1) += tmp_8_13;
    matrix(indices[2]+2,indices[4]+2) += tmp_8_14;
    matrix(indices[2]+2,indices[5]+0) += tmp_8_15;
    matrix(indices[2]+2,indices[5]+1) += tmp_8_16;
    matrix(indices[2]+2,indices[5]+2) += tmp_8_17;
    matrix(indices[3]+0,indices[3]+0) += tmp_9_9;
    matrix(indices[3]+0,indices[3]+1) += tmp_9_10;
    matrix(indices[3]+0,indices[3]+2) += tmp_9_11;
    matrix(indices[3]+0,indices[4]+0) += tmp_9_12;
    matrix(indices[3]+0,indices[4]+1) += tmp_9_13;
    matrix(indices[3]+0,indices[4]+2) += tmp_9_14;
    matrix(indices[3]+0,indices[5]+0) += tmp_9_15;
    matrix(indices[3]+0,indices[5]+1) += tmp_9_16;
    matrix(indices[3]+0,indices[5]+2) += tmp_9_17;
    matrix(indices[3]+1,indices[3]+1) += tmp_10_10;
    matrix(indices[3]+1,indices[3]+2) += tmp_10_11;
    matrix(indices[3]+1,indices[4]+0) += tmp_10_12;
    matrix(indices[3]+1,indices[4]+1) += tmp_10_13;
    matrix(indices[3]+1,indices[4]+2) += tmp_10_14;
    matrix(indices[3]+1,indices[5]+0) += tmp_10_15;
    matrix(indices[3]+1,indices[5]+1) += tmp_10_16;
    matrix(indices[3]+1,indices[5]+2) += tmp_10_17;
    matrix(indices[3]+2,indices[3]+2) += tmp_11_11;
    matrix(indices[3]+2,indices[4]+0) += tmp_11_12;
    matrix(indices[3]+2,indices[4]+1) += tmp_11_13;
    matrix(indices[3]+2,indices[4]+2) += tmp_11_14;
    matrix(indices[3]+2,indices[5]+0) += tmp_11_15;
    matrix(indices[3]+2,indices[5]+1) += tmp_11_16;
    matrix(indices[3]+2,indices[5]+2) += tmp_11_17;
    matrix(indices[4]+0,indices[4]+0) += tmp_12_12;
    matrix(indices[4]+0,indices[4]+1) += tmp_12_13;
    matrix(indices[4]+0,indices[4]+2) += tmp_12_14;
    matrix(indices[4]+0,indices[5]+0) += tmp_12_15;
    matrix(indices[4]+0,indices[5]+1) += tmp_12_16;
    matrix(indices[4]+0,indices[5]+2) += tmp_12_17;
    matrix(indices[4]+1,indices[4]+1) += tmp_13_13;
    matrix(indices[4]+1,indices[4]+2) += tmp_13_14;
    matrix(indices[4]+1,indices[5]+0) += tmp_13_15;
    matrix(indices[4]+1,indices[5]+1) += tmp_13_16;
    matrix(indices[4]+1,indices[5]+2) += tmp_13_17;
    matrix(indices[4]+2,indices[4]+2) += tmp_14_14;
    matrix(indices[4]+2,indices[5]+0) += tmp_14_15;
    matrix(indices[4]+2,indices[5]+1) += tmp_14_16;
    matrix(indices[4]+2,indices[5]+2) += tmp_14_17;
    matrix(indices[5]+0,indices[5]+0) += tmp_15_15;
    matrix(indices[5]+0,indices[5]+1) += tmp_15_16;
    matrix(indices[5]+0,indices[5]+2) += tmp_15_17;
    matrix(indices[5]+1,indices[5]+1) += tmp_16_16;
    matrix(indices[5]+1,indices[5]+2) += tmp_16_17;
    matrix(indices[5]+2,indices[5]+2) += tmp_17_17;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim,unsigned symmetric_version>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1-var_inter[0]; reg0=reg0-var_inter[1]; T reg1=1+(*f.m).poisson_ratio; T reg2=reg0*elem.pos(0)[2]; T reg3=var_inter[0]*elem.pos(1)[2];
    T reg4=reg0*elem.pos(0)[1]; T reg5=var_inter[0]*elem.pos(1)[1]; T reg6=var_inter[1]*elem.pos(2)[2]; T reg7=reg3+reg2; T reg8=reg5+reg4;
    T reg9=var_inter[1]*elem.pos(2)[1]; T reg10=1-var_inter[2]; reg1=reg1/(*f.m).elastic_modulus; T reg11=reg7+reg6; T reg12=elem.pos(3)[2]*reg0;
    T reg13=reg10*elem.pos(1)[2]; T reg14=reg10*elem.pos(0)[2]; T reg15=pow(reg1,2); T reg16=reg10*elem.pos(0)[1]; T reg17=reg10*elem.pos(1)[1];
    T reg18=reg10*elem.pos(2)[1]; T reg19=reg9+reg8; T reg20=reg10*elem.pos(2)[2]; T reg21=elem.pos(3)[1]*reg0; reg17=reg17-reg16;
    T reg22=var_inter[2]*elem.pos(3)[1]; reg13=reg13-reg14; T reg23=var_inter[2]*elem.pos(3)[2]; T reg24=elem.pos(0)[0]*reg0; T reg25=var_inter[0]*elem.pos(1)[0];
    reg18=reg18-reg16; reg20=reg20-reg14; T reg26=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg27=1.0/(*f.m).elastic_modulus; T reg28=var_inter[0]*elem.pos(4)[2];
    reg12=reg12-reg11; reg1=reg1*reg15; reg21=reg21-reg19; T reg29=var_inter[0]*elem.pos(4)[1]; T reg30=var_inter[2]*elem.pos(5)[2];
    reg20=reg20-reg23; reg18=reg18-reg22; reg29=reg21+reg29; reg21=var_inter[2]*elem.pos(5)[1]; T reg31=reg10*elem.pos(1)[0];
    T reg32=reg10*elem.pos(2)[0]; T reg33=reg25+reg24; T reg34=var_inter[1]*elem.pos(2)[0]; T reg35=var_inter[2]*elem.pos(4)[2]; reg13=reg13-reg23;
    T reg36=var_inter[1]*elem.pos(5)[2]; T reg37=var_inter[1]*elem.pos(5)[1]; T reg38=reg26*reg1; reg1=reg27*reg1; reg17=reg17-reg22;
    T reg39=var_inter[2]*elem.pos(4)[1]; T reg40=reg10*elem.pos(0)[0]; reg28=reg12+reg28; reg28=reg36+reg28; reg12=reg33+reg34;
    reg20=reg30+reg20; reg30=reg26*reg1; reg36=reg26*reg38; reg29=reg37+reg29; reg1=reg27*reg1;
    reg37=elem.pos(3)[0]*reg0; reg31=reg31-reg40; T reg41=var_inter[2]*elem.pos(3)[0]; reg39=reg17+reg39; reg35=reg13+reg35;
    reg32=reg32-reg40; reg18=reg21+reg18; reg38=reg27*reg38; reg13=reg35*reg29; reg30=reg36+reg30;
    reg17=reg20*reg29; reg21=reg39*reg28; reg1=reg1-reg36; reg31=reg31-reg41; T reg42=var_inter[2]*elem.pos(4)[0];
    reg32=reg32-reg41; T reg43=reg18*reg28; T reg44=var_inter[0]*elem.pos(4)[0]; reg37=reg37-reg12; T reg45=var_inter[2]*elem.pos(5)[0];
    reg17=reg43-reg17; reg43=reg27*reg1; reg13=reg21-reg13; reg21=reg39*reg20; T reg46=reg35*reg18;
    T reg47=reg26*reg30; reg38=reg36+reg38; reg32=reg45+reg32; reg36=var_inter[1]*elem.pos(5)[0]; reg44=reg37+reg44;
    reg42=reg31+reg42; reg46=reg21-reg46; reg21=reg32*reg13; reg31=reg42*reg17; reg44=reg36+reg44;
    reg47=reg43-reg47; reg36=reg26*reg38; reg37=reg42*reg29; reg43=reg18*reg44; reg45=reg42*reg28;
    reg29=reg32*reg29; T reg48=reg20*reg44; T reg49=reg39*reg32; reg28=reg32*reg28; reg36=reg47-reg36;
    reg47=reg44*reg46; reg18=reg42*reg18; reg21=reg31-reg21; reg32=reg35*reg32; reg20=reg42*reg20;
    reg39=reg39*reg44; reg44=reg35*reg44; reg38=reg38/reg36; reg1=reg1/reg36; reg31=(*f.m).deltaT*(*f.m).alpha;
    reg30=reg30/reg36; reg49=reg18-reg49; reg43=reg29-reg43; reg32=reg20-reg32; reg48=reg28-reg48;
    reg39=reg37-reg39; reg47=reg21+reg47; reg44=reg45-reg44; reg49=reg49/reg47; reg32=reg32/reg47;
    reg46=reg46/reg47; reg39=reg39/reg47; reg44=reg44/reg47; reg13=reg13/reg47; reg18=reg1*reg31;
    reg20=reg30*reg31; reg17=reg17/reg47; reg48=reg48/reg47; reg43=reg43/reg47; reg21=reg38*reg31;
    reg28=reg10*reg48; reg29=reg10*reg44; reg35=var_inter[1]*reg49; reg37=var_inter[1]*reg46; reg42=reg10*reg13;
    reg45=var_inter[2]*reg39; T reg50=var_inter[2]*reg43; T reg51=var_inter[2]*reg17; T reg52=var_inter[2]*reg13; T reg53=var_inter[2]*reg48;
    T reg54=var_inter[2]*reg44; T reg55=var_inter[0]*reg32; T reg56=reg20+reg21; T reg57=reg18+reg20; T reg58=reg10*reg17;
    T reg59=reg10*reg39; T reg60=reg10*reg43; T reg61=reg42-reg58; T reg62=reg35+reg59; T reg63=reg0*reg46;
    T reg64=var_inter[1]*reg32; T reg65=reg42+reg37; T reg66=reg57+reg21; T reg67=reg45-reg50; T reg68=reg0*reg49;
    T reg69=reg18+reg56; T reg70=reg53-reg54; T reg71=reg59-reg60; T reg72=var_inter[0]*reg46; T reg73=reg28-reg29;
    T reg74=reg0*reg32; T reg75=reg52-reg51; T reg76=var_inter[0]*var_inter[2]; T reg77=var_inter[0]*reg49; T reg78=var_inter[1]*reg10;
    T reg79=reg53+reg55; T reg80=reg35-reg45; T reg81=reg54-reg64; reg73=reg73+reg74; T reg82=reg37-reg52;
    T reg83=reg76*(*f.m).f_vol[1]; T reg84=reg10*reg0; reg71=reg71-reg68; T reg85=reg65*reg66; T reg86=reg10*var_inter[0];
    T reg87=reg58-reg72; T reg88=var_inter[1]*var_inter[2]; T reg89=var_inter[2]*reg0; reg67=reg67+reg68; T reg90=reg60-reg77;
    T reg91=reg55-reg28; reg70=reg70-reg74; T reg92=reg77+reg50; T reg93=reg64+reg29; T reg94=reg51+reg72;
    T reg95=reg62*reg69; T reg96=reg78*(*f.m).f_vol[0]; T reg97=reg79*reg66; T reg98=reg78*(*f.m).f_vol[2]; reg75=reg75+reg63;
    reg61=reg61-reg63; T reg99=reg67*reg69; T reg100=reg71*reg69; T reg101=reg81*reg66; T reg102=reg87*reg66;
    T reg103=reg70*reg66; T reg104=reg80*reg69; T reg105=reg91*reg66; T reg106=reg82*reg66; T reg107=reg94*reg66;
    T reg108=reg75*reg66; T reg109=reg90*reg69; T reg110=reg92*reg69; T reg111=reg95-reg98; T reg112=reg85-reg96;
    T reg113=reg93*reg66; T reg114=reg97-reg83; T reg115=reg86*(*f.m).f_vol[1]; T reg116=reg86*(*f.m).f_vol[0]; T reg117=reg89*(*f.m).f_vol[0];
    T reg118=reg86*(*f.m).f_vol[2]; T reg119=reg76*(*f.m).f_vol[0]; T reg120=reg89*(*f.m).f_vol[2]; T reg121=reg89*(*f.m).f_vol[1]; T reg122=reg78*(*f.m).f_vol[1];
    T reg123=reg88*(*f.m).f_vol[0]; T reg124=reg76*(*f.m).f_vol[2]; T reg125=reg84*(*f.m).f_vol[0]; T reg126=reg84*(*f.m).f_vol[1]; T reg127=reg84*(*f.m).f_vol[2];
    T reg128=reg88*(*f.m).f_vol[1]; T reg129=reg88*(*f.m).f_vol[2]; T reg130=reg73*reg66; T reg131=reg61*reg66; T reg132=reg126+reg130;
    T reg133=reg120+reg99; T reg134=reg121+reg103; T reg135=reg117+reg108; T reg136=reg129+reg104; T reg137=reg119+reg107;
    reg111=reg47*reg111; T reg138=reg125+reg131; T reg139=reg127+reg100; T reg140=reg128+reg101; T reg141=reg122+reg113;
    reg114=reg47*reg114; reg112=reg47*reg112; T reg142=reg116+reg102; T reg143=reg124+reg110; T reg144=reg118+reg109;
    T reg145=reg123+reg106; T reg146=reg115+reg105; reg114=ponderation*reg114; T reg147=reg47*reg136; T reg148=reg47*reg137;
    T reg149=reg47*reg145; T reg150=reg47*reg140; T reg151=reg47*reg133; T reg152=reg47*reg143; T reg153=reg47*reg146;
    T reg154=reg47*reg142; T reg155=reg47*reg144; reg112=ponderation*reg112; T reg156=reg47*reg139; T reg157=reg47*reg141;
    reg111=ponderation*reg111; T reg158=reg47*reg132; T reg159=reg47*reg135; T reg160=reg47*reg138; T reg161=reg47*reg134;
    T reg162=ponderation*reg150; sollicitation[indices[5]+1]+=reg162; T reg163=ponderation*reg158; sollicitation[indices[0]+1]+=reg163; T reg164=ponderation*reg160;
    sollicitation[indices[0]+0]+=reg164; T reg165=ponderation*reg156; sollicitation[indices[0]+2]+=reg165; T reg166=ponderation*reg149; sollicitation[indices[5]+0]+=reg166;
    T reg167=ponderation*reg154; sollicitation[indices[1]+0]+=reg167; T reg168=ponderation*reg147; sollicitation[indices[5]+2]+=reg168; T reg169=ponderation*reg152;
    sollicitation[indices[4]+2]+=reg169; T reg170=ponderation*reg153; sollicitation[indices[1]+1]+=reg170; T reg171=ponderation*reg155; sollicitation[indices[1]+2]+=reg171;
    sollicitation[indices[4]+1]+=-reg114; sollicitation[indices[2]+0]+=-reg112; reg112=ponderation*reg157; sollicitation[indices[2]+1]+=reg112; reg114=ponderation*reg148;
    sollicitation[indices[4]+0]+=reg114; sollicitation[indices[2]+2]+=-reg111; reg111=ponderation*reg151; sollicitation[indices[3]+2]+=reg111; T reg172=ponderation*reg159;
    sollicitation[indices[3]+0]+=reg172; T reg173=ponderation*reg161; sollicitation[indices[3]+1]+=reg173;
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim,unsigned symmetric_version>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_residual( TD ponderation, const TD *var_inter,
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices ) { 
  #define PNODE(N) (*elem.node(N))
    T reg0=1-var_inter[0]; reg0=reg0-var_inter[1]; T reg1=var_inter[0]*elem.pos(1)[2]; T reg2=reg0*elem.pos(0)[2]; T reg3=reg0*elem.pos(0)[1];
    T reg4=var_inter[0]*elem.pos(1)[1]; T reg5=reg1+reg2; T reg6=var_inter[1]*elem.pos(2)[2]; T reg7=reg4+reg3; T reg8=var_inter[1]*elem.pos(2)[1];
    T reg9=1-var_inter[2]; T reg10=reg8+reg7; T reg11=elem.pos(3)[1]*reg0; T reg12=reg5+reg6; T reg13=reg9*elem.pos(1)[1];
    T reg14=reg9*elem.pos(0)[1]; T reg15=reg9*elem.pos(1)[2]; T reg16=reg9*elem.pos(0)[2]; T reg17=reg9*elem.pos(2)[1]; T reg18=reg9*elem.pos(2)[2];
    T reg19=elem.pos(3)[2]*reg0; T reg20=var_inter[2]*elem.pos(3)[1]; reg13=reg13-reg14; reg15=reg15-reg16; T reg21=var_inter[2]*elem.pos(3)[2];
    T reg22=var_inter[0]*elem.pos(1)[0]; T reg23=elem.pos(0)[0]*reg0; reg18=reg18-reg16; reg17=reg17-reg14; reg11=reg11-reg10;
    T reg24=var_inter[0]*elem.pos(4)[1]; reg19=reg19-reg12; T reg25=var_inter[0]*elem.pos(4)[2]; T reg26=var_inter[2]*elem.pos(4)[1]; reg13=reg13-reg20;
    reg25=reg19+reg25; reg19=var_inter[2]*elem.pos(5)[2]; T reg27=reg22+reg23; reg15=reg15-reg21; T reg28=var_inter[2]*elem.pos(4)[2];
    T reg29=reg9*elem.pos(0)[0]; T reg30=reg9*elem.pos(1)[0]; T reg31=var_inter[1]*elem.pos(2)[0]; T reg32=reg9*elem.pos(2)[0]; T reg33=var_inter[1]*elem.pos(5)[1];
    reg18=reg18-reg21; T reg34=var_inter[2]*elem.pos(5)[1]; T reg35=var_inter[1]*elem.pos(5)[2]; reg17=reg17-reg20; reg24=reg11+reg24;
    reg28=reg15+reg28; reg18=reg19+reg18; reg32=reg32-reg29; reg17=reg34+reg17; reg11=1+(*f.m).poisson_ratio;
    reg26=reg13+reg26; reg13=var_inter[2]*elem.pos(3)[0]; reg30=reg30-reg29; reg15=reg27+reg31; reg19=elem.pos(3)[0]*reg0;
    reg24=reg33+reg24; reg25=reg35+reg25; reg11=reg11/(*f.m).elastic_modulus; reg33=reg17*reg25; reg34=reg26*reg25;
    reg35=var_inter[0]*elem.pos(4)[0]; reg19=reg19-reg15; T reg36=reg18*reg24; T reg37=reg28*reg24; reg32=reg32-reg13;
    T reg38=var_inter[2]*elem.pos(5)[0]; T reg39=var_inter[2]*elem.pos(4)[0]; reg30=reg30-reg13; T reg40=var_inter[0]*vectors[0][indices[1]+2]; T reg41=reg0*vectors[0][indices[0]+2];
    reg36=reg33-reg36; reg37=reg34-reg37; reg33=reg26*reg18; reg34=reg28*reg17; T reg42=var_inter[0]*vectors[0][indices[1]+0];
    T reg43=reg0*vectors[0][indices[0]+0]; T reg44=pow(reg11,2); reg32=reg38+reg32; reg38=var_inter[0]*vectors[0][indices[1]+1]; T reg45=reg0*vectors[0][indices[0]+1];
    T reg46=var_inter[1]*elem.pos(5)[0]; reg39=reg30+reg39; reg35=reg19+reg35; reg35=reg46+reg35; reg19=reg32*reg37;
    reg30=reg9*vectors[0][indices[0]+1]; reg46=var_inter[1]*vectors[0][indices[2]+1]; T reg47=reg9*vectors[0][indices[1]+0]; reg34=reg33-reg34; reg33=reg9*vectors[0][indices[1]+1];
    T reg48=var_inter[1]*vectors[0][indices[2]+0]; reg42=reg43+reg42; reg38=reg45+reg38; reg43=reg9*vectors[0][indices[2]+0]; reg11=reg11*reg44;
    reg45=var_inter[1]*vectors[0][indices[2]+2]; T reg49=reg9*vectors[0][indices[2]+2]; T reg50=reg9*vectors[0][indices[0]+0]; reg41=reg40+reg41; reg40=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg51=reg9*vectors[0][indices[0]+2]; T reg52=reg9*vectors[0][indices[1]+2]; T reg53=1.0/(*f.m).elastic_modulus; T reg54=reg9*vectors[0][indices[2]+1]; T reg55=reg39*reg36;
    T reg56=reg40*reg11; T reg57=var_inter[2]*vectors[0][indices[3]+1]; reg42=reg48+reg42; reg11=reg53*reg11; reg33=reg33-reg30;
    reg30=reg54-reg30; reg48=reg0*vectors[0][indices[3]+1]; reg46=reg38+reg46; reg52=reg52-reg51; reg41=reg45+reg41;
    reg38=var_inter[2]*vectors[0][indices[3]+2]; reg45=reg0*vectors[0][indices[3]+2]; reg54=reg0*vectors[0][indices[3]+0]; reg51=reg49-reg51; reg49=reg32*reg25;
    reg47=reg47-reg50; reg50=reg43-reg50; reg19=reg55-reg19; reg43=reg35*reg34; reg55=reg26*reg35;
    T reg58=reg17*reg35; T reg59=reg39*reg24; reg25=reg39*reg25; reg24=reg32*reg24; T reg60=reg18*reg35;
    T reg61=var_inter[2]*vectors[0][indices[3]+0]; reg35=reg28*reg35; reg55=reg59-reg55; reg18=reg39*reg18; reg59=var_inter[2]*vectors[0][indices[5]+0];
    reg52=reg52-reg38; reg28=reg28*reg32; T reg62=var_inter[2]*vectors[0][indices[4]+2]; T reg63=var_inter[0]*vectors[0][indices[4]+0]; reg41=reg45-reg41;
    reg58=reg24-reg58; reg42=reg54-reg42; reg35=reg25-reg35; reg50=reg50-reg61; reg60=reg49-reg60;
    reg17=reg39*reg17; reg24=var_inter[2]*vectors[0][indices[5]+2]; reg25=reg53*reg11; reg32=reg26*reg32; reg26=reg40*reg56;
    reg11=reg40*reg11; reg39=var_inter[0]*vectors[0][indices[4]+2]; reg38=reg51-reg38; reg30=reg30-reg57; reg61=reg47-reg61;
    reg45=reg40*reg44; reg47=var_inter[2]*vectors[0][indices[5]+1]; reg44=reg53*reg44; reg57=reg33-reg57; reg43=reg19+reg43;
    reg46=reg48-reg46; reg19=var_inter[0]*vectors[0][indices[4]+1]; reg33=var_inter[2]*vectors[0][indices[4]+1]; reg48=var_inter[2]*vectors[0][indices[4]+0]; reg36=reg36/reg43;
    reg55=reg55/reg43; reg52=reg62+reg52; reg30=reg47+reg30; reg35=reg35/reg43; reg33=reg57+reg33;
    reg24=reg38+reg24; reg25=reg25-reg26; reg38=var_inter[1]*vectors[0][indices[5]+1]; reg32=reg17-reg32; reg11=reg26+reg11;
    reg56=reg53*reg56; reg17=reg53*reg44; reg47=reg40*reg45; reg44=reg40*reg44; reg19=reg46+reg19;
    reg48=reg61+reg48; reg59=reg50+reg59; reg46=var_inter[1]*vectors[0][indices[5]+0]; reg28=reg18-reg28; reg18=var_inter[1]*vectors[0][indices[5]+2];
    reg41=reg39+reg41; reg37=reg37/reg43; reg63=reg42+reg63; reg58=reg58/reg43; reg60=reg60/reg43;
    reg19=reg38+reg19; reg28=reg28/reg43; reg38=reg37*reg59; reg34=reg34/reg43; reg39=reg60*reg33;
    reg32=reg32/reg43; reg42=reg36*reg48; reg63=reg46+reg63; reg17=reg17-reg47; reg45=reg53*reg45;
    reg18=reg41+reg18; reg56=reg26+reg56; reg26=reg36*reg52; reg41=reg37*reg24; reg46=reg55*reg59;
    reg49=reg37*reg30; reg50=reg36*reg33; reg51=reg58*reg48; reg54=reg35*reg30; reg57=reg53*reg25;
    reg61=reg40*reg11; reg48=reg60*reg48; reg59=reg35*reg59; reg44=reg47+reg44; reg46=reg51-reg46;
    reg51=reg32*reg63; reg41=reg26-reg41; reg26=reg34*reg18; reg30=reg55*reg30; reg33=reg58*reg33;
    reg38=reg42-reg38; reg42=reg34*reg19; reg49=reg50-reg49; reg50=reg28*reg63; reg62=reg58*reg52;
    reg48=reg59-reg48; reg59=reg55*reg24; reg52=reg60*reg52; reg63=reg34*reg63; T reg64=reg47+reg45;
    reg44=reg40*reg44; reg17=reg53*reg17; reg39=reg54-reg39; reg24=reg35*reg24; reg53=reg28*reg19;
    reg61=reg57-reg61; reg54=reg40*reg56; reg54=reg61-reg54; reg44=reg17-reg44; reg17=(*f.m).deltaT*(*f.m).alpha;
    reg64=reg40*reg64; reg40=reg32*reg18; reg46=reg51+reg46; reg38=reg63+reg38; reg26=reg41+reg26;
    reg42=reg49+reg42; reg50=reg48-reg50; reg59=reg62-reg59; reg18=reg28*reg18; reg52=reg24-reg52;
    reg53=reg39-reg53; reg19=reg32*reg19; reg30=reg33-reg30; reg24=reg9*reg60; reg42=reg50+reg42;
    reg33=reg9*reg35; reg39=var_inter[2]*reg36; reg41=reg9*reg36; reg38=reg38-reg17; reg26=reg46+reg26;
    reg40=reg59+reg40; reg64=reg44-reg64; reg44=var_inter[2]*reg60; reg46=var_inter[2]*reg37; reg53=reg53-reg17;
    reg11=reg11/reg54; reg25=reg25/reg54; reg30=reg19+reg30; reg18=reg52-reg18; reg19=var_inter[2]*reg35;
    reg56=reg56/reg54; reg48=reg9*reg37; reg49=var_inter[0]*reg34; reg50=reg24-reg33; reg51=reg0*reg28;
    reg52=var_inter[0]*reg28; reg57=reg9*reg58; reg59=reg11*reg53; reg61=reg25*reg38; reg40=reg40-reg17;
    reg62=var_inter[2]*reg55; reg38=reg11*reg38; reg63=reg25*reg53; reg54=reg64/reg54; reg64=var_inter[2]*reg58;
    T reg65=reg9*reg55; reg18=reg30+reg18; reg30=var_inter[1]*reg34; T reg66=var_inter[1]*reg28; T reg67=reg46-reg39;
    reg26=0.5*reg26; reg53=reg56*reg53; T reg68=reg48-reg41; reg42=0.5*reg42; T reg69=reg0*reg34;
    T reg70=reg44-reg19; reg42=reg54*reg42; reg18=0.5*reg18; T reg71=reg41-reg49; reg59=reg61+reg59;
    reg67=reg67+reg69; reg26=reg54*reg26; reg70=reg70-reg51; reg61=reg48+reg30; T reg72=reg62-reg64;
    T reg73=reg66+reg33; T reg74=reg52-reg24; T reg75=reg65-reg57; T reg76=var_inter[1]*reg32; T reg77=reg0*reg32;
    reg50=reg50+reg51; T reg78=reg30-reg46; T reg79=reg44+reg52; T reg80=var_inter[0]*reg32; T reg81=reg25*reg40;
    reg53=reg53+reg38; reg40=reg56*reg40; reg68=reg68-reg69; reg63=reg38+reg63; reg38=reg19-reg66;
    T reg82=reg39+reg49; T reg83=0.5*reg74; T reg84=reg57-reg80; reg26=2*reg26; T reg85=0.5*reg71;
    T reg86=0.5*reg78; T reg87=0.5*reg38; T reg88=reg76-reg62; T reg89=0.5*reg73; T reg90=reg76+reg65;
    T reg91=0.5*reg61; reg18=reg54*reg18; T reg92=0.5*reg67; reg72=reg72+reg77; T reg93=0.5*reg70;
    T reg94=0.5*reg79; T reg95=reg80+reg64; T reg96=0.5*reg82; reg63=reg63+reg40; reg81=reg53+reg81;
    reg42=2*reg42; reg40=reg59+reg40; reg75=reg75-reg77; reg53=0.5*reg68; reg59=0.5*reg50;
    T reg97=reg40*reg71; T reg98=reg40*reg68; T reg99=reg26*reg85; T reg100=reg81*reg84; T reg101=0.5*reg90;
    T reg102=reg42*reg85; T reg103=reg42*reg59; T reg104=reg42*reg89; T reg105=reg40*reg61; reg18=2*reg18;
    T reg106=0.5*reg72; T reg107=reg42*reg93; T reg108=reg40*reg67; T reg109=0.5*reg75; T reg110=reg63*reg50;
    T reg111=reg42*reg53; T reg112=reg81*reg75; T reg113=reg26*reg53; T reg114=reg26*reg86; T reg115=reg92*reg26;
    T reg116=reg96*reg26; T reg117=reg95*reg81; T reg118=reg81*reg72; T reg119=reg78*reg40; T reg120=reg87*reg42;
    T reg121=reg96*reg42; T reg122=reg42*reg86; T reg123=reg82*reg40; T reg124=reg94*reg42; T reg125=reg79*reg63;
    T reg126=0.5*reg95; T reg127=0.5*reg88; T reg128=reg42*reg83; T reg129=0.5*reg84; T reg130=reg63*reg74;
    T reg131=reg38*reg63; T reg132=reg88*reg81; T reg133=reg92*reg42; T reg134=reg63*reg70; T reg135=reg81*reg90;
    T reg136=reg63*reg73; T reg137=reg91*reg42; T reg138=reg91*reg26; T reg139=reg26*reg101; T reg140=reg18*reg101;
    reg136=reg136-reg137; reg99=reg100+reg99; reg100=reg83*reg18; reg104=reg104-reg105; reg120=reg119+reg120;
    reg119=reg127*reg26; T reg141=reg9*reg0; T reg142=reg9*var_inter[0]; reg113=reg112+reg113; reg112=reg18*reg59;
    T reg143=reg18*reg109; reg111=reg110+reg111; reg110=var_inter[1]*var_inter[2]; T reg144=var_inter[1]*reg9; T reg145=var_inter[2]*reg0;
    T reg146=reg26*reg109; T reg147=var_inter[0]*var_inter[2]; reg103=reg98+reg103; reg121=reg121-reg125; reg98=reg126*reg18;
    T reg148=reg87*reg18; reg116=reg117+reg116; reg117=reg94*reg18; reg114=reg132+reg114; reg132=reg26*reg106;
    T reg149=reg93*reg18; reg115=reg118+reg115; reg118=reg127*reg18; reg122=reg131+reg122; reg107=reg108+reg107;
    reg123=reg123-reg124; reg108=reg106*reg18; reg131=reg126*reg26; reg133=reg134+reg133; reg97=reg128+reg97;
    reg128=reg26*reg129; reg102=reg130+reg102; reg130=reg129*reg18; reg134=reg135+reg138; T reg150=reg18*reg89;
    reg146=reg103+reg146; reg128=reg97+reg128; reg149=reg115+reg149; reg97=reg145*(*f.m).f_vol[2]; reg148=reg114+reg148;
    reg103=reg147*(*f.m).f_vol[1]; reg98=reg121+reg98; reg130=reg102+reg130; reg102=reg142*(*f.m).f_vol[1]; reg114=reg141*(*f.m).f_vol[1];
    reg115=reg141*(*f.m).f_vol[0]; reg121=reg141*(*f.m).f_vol[2]; T reg151=reg110*(*f.m).f_vol[2]; reg143=reg111+reg143; reg111=reg142*(*f.m).f_vol[0];
    reg131=reg123+reg131; reg123=reg147*(*f.m).f_vol[0]; reg112=reg113+reg112; reg100=reg99+reg100; reg99=reg142*(*f.m).f_vol[2];
    reg113=reg144*(*f.m).f_vol[2]; T reg152=reg110*(*f.m).f_vol[0]; reg150=reg150-reg134; reg119=reg120+reg119; reg120=reg110*(*f.m).f_vol[1];
    reg108=reg133+reg108; reg133=reg145*(*f.m).f_vol[1]; reg118=reg122+reg118; reg136=reg136-reg140; reg122=reg144*(*f.m).f_vol[1];
    reg104=reg104-reg139; T reg153=reg144*(*f.m).f_vol[0]; T reg154=reg145*(*f.m).f_vol[0]; T reg155=reg147*(*f.m).f_vol[2]; reg116=reg116-reg117;
    reg132=reg107+reg132; reg143=reg143-reg114; reg136=reg136-reg122; reg118=reg118-reg120; reg130=reg130-reg102;
    reg119=reg119-reg152; reg128=reg128-reg111; reg150=reg150-reg113; reg146=reg146-reg115; reg100=reg100-reg99;
    reg112=reg112-reg121; reg104=reg104-reg153; reg131=reg131-reg123; reg148=reg148-reg151; reg108=reg108-reg133;
    reg149=reg149-reg97; reg132=reg132-reg154; reg98=reg98-reg103; reg116=reg116-reg155; reg112=reg43*reg112;
    reg118=reg43*reg118; reg148=reg43*reg148; reg119=reg43*reg119; reg136=reg43*reg136; reg108=reg43*reg108;
    reg149=reg43*reg149; reg128=reg43*reg128; reg104=reg43*reg104; reg130=reg43*reg130; reg116=reg43*reg116;
    reg150=reg43*reg150; reg132=reg43*reg132; reg98=reg43*reg98; reg100=reg43*reg100; reg131=reg43*reg131;
    reg146=reg43*reg146; reg143=reg43*reg143; sollicitation[indices[4]+1]+=ponderation*reg98; sollicitation[indices[4]+0]+=ponderation*reg131; sollicitation[indices[4]+2]+=ponderation*reg116;
    sollicitation[indices[5]+1]+=ponderation*reg118; sollicitation[indices[2]+0]+=ponderation*reg104; sollicitation[indices[5]+2]+=ponderation*reg148; sollicitation[indices[2]+1]+=ponderation*reg136; sollicitation[indices[3]+1]+=ponderation*reg108;
    sollicitation[indices[3]+2]+=ponderation*reg149; sollicitation[indices[1]+0]+=ponderation*reg128; sollicitation[indices[5]+0]+=ponderation*reg119; sollicitation[indices[1]+1]+=ponderation*reg130; sollicitation[indices[2]+2]+=ponderation*reg150;
    sollicitation[indices[3]+0]+=ponderation*reg132; sollicitation[indices[0]+2]+=ponderation*reg112; sollicitation[indices[1]+2]+=ponderation*reg100; sollicitation[indices[0]+0]+=ponderation*reg146; sollicitation[indices[0]+1]+=ponderation*reg143;
  #undef PNODE
}
// 
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_false_true_true
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_false_true_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_true_true_true
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_true_true_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_false_true_false
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_false_true_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_true_true_false
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_true_true_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_symmetric_version_false_true
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_symmetric_version_false_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE ,unsigned symmetric_version>
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_symmetric_version_false_false
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_symmetric_version_false_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE ,unsigned symmetric_version>
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_RESIDUAL_elasticity_isotropy_stat_Qstat
#define ADD_NODAL_RESIDUAL_elasticity_isotropy_stat_Qstat
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE>
void add_nodal_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<4> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<4> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<4> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<4> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<4> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<4> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Wedge,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Triangle,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<4> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}

#ifndef elasticity_isotropy_stat_Qstat_read_material_to_mesh
#define elasticity_isotropy_stat_Qstat_read_material_to_mesh
template<class TM, class T, bool wont_add_nz>
void read_material_to_mesh_(const XmlNode &n, Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f){ 
    if(n.has_attribute("elastic_modulus"))  
        n.get_attribute("elastic_modulus", f.m->elastic_modulus ); 
    else  
        std::cerr << "Warning using default value of elastic_modulus : " << f.m->elastic_modulus << std::endl; 

    if(n.has_attribute("density"))  
        n.get_attribute("density", f.m->density ); 
    else  
        std::cerr << "Warning using default value of density : " << f.m->density << std::endl; 

    if(n.has_attribute("deltaT"))  
        n.get_attribute("deltaT", f.m->deltaT ); 
    else  
        std::cerr << "Warning using default value of deltaT : " << f.m->deltaT << std::endl; 

    if(n.has_attribute("poisson_ratio"))  
        n.get_attribute("poisson_ratio", f.m->poisson_ratio ); 
    else  
        std::cerr << "Warning using default value of poisson_ratio : " << f.m->poisson_ratio << std::endl; 

    if(n.has_attribute("alpha"))  
        n.get_attribute("alpha", f.m->alpha ); 
    else  
        std::cerr << "Warning using default value of alpha : " << f.m->alpha << std::endl; 

    if(n.has_attribute("resolution"))  
        n.get_attribute("resolution", f.m->resolution ); 
    else  
        std::cerr << "Warning using default value of resolution : " << f.m->resolution << std::endl; 

    if(n.has_attribute("f_vol"))  
        n.get_attribute("f_vol", f.m->f_vol ); 
    else  
        std::cerr << "Warning using default value of f_vol : " << f.m->f_vol << std::endl; 

  };
#endif // elasticity_isotropy_stat_Qstat_read_material_to_mesh
} // namespace LMT


#include "formulation/formulation.h"
namespace LMT {
#ifndef ELASTICITY_ISOTROPY_STAT_QSTAT_3_NUM_0_CARAC_H
#define ELASTICITY_ISOTROPY_STAT_QSTAT_3_NUM_0_CARAC_H
#ifndef STRUCT_ELASTICITY_ISOTROPY_STAT_QSTAT
#define STRUCT_ELASTICITY_ISOTROPY_STAT_QSTAT
struct elasticity_isotropy_stat_Qstat {
  static const char *name() { return "elasticity_isotropy_stat_Qstat"; }
};
#endif // ELASTICITY_ISOTROPY_STAT_QSTAT

template<class P_T>
class CaracFormulation<elasticity_isotropy_stat_Qstat,3,P_T>  {
public:
  typedef P_T T;
  static const char *name() { return "elasticity_isotropy_stat_Qstat"; }
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
  
  static const unsigned nb_nodal_unknowns = 3;
  template<class TE,class TTs,class Tvec>
  inline static void set_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice) {
    node.dep[0]=vecs[0][indice+0]; node.dep[2]=vecs[0][indice+2]; node.dep[1]=vecs[0][indice+1];
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
    node.dep[2]=vecs[0][indice+2]; node.dep[1]=vecs[0][indice+1]; node.dep[0]=vecs[0][indice+0];
  }
  template<class TE,class TTs,class Tvec>
  inline static void get_nodal_initial_conditions(const TE &node,const TTs &f,Tvec &vecs,unsigned indice) {
    vecs[0][indice+0]=node.dep[0]; vecs[1][indice+0]=node.dep[0]; vecs[2][indice+0]=node.dep[0]; vecs[3][indice+0]=node.dep[0]; vecs[4][indice+0]=node.dep[0];
    vecs[0][indice+1]=node.dep[1]; vecs[1][indice+1]=node.dep[1]; vecs[2][indice+1]=node.dep[1]; vecs[3][indice+1]=node.dep[1]; vecs[4][indice+1]=node.dep[1];
    vecs[0][indice+2]=node.dep[2]; vecs[1][indice+2]=node.dep[2]; vecs[2][indice+2]=node.dep[2]; vecs[3][indice+2]=node.dep[2]; vecs[4][indice+2]=node.dep[2];
  }
  template<class TE,class TTs,class Tvec>
  inline static T max_nodal_error(const TE &node,const TTs &f,const Tvec &vecs,int indice) {
    T reg0=vecs[1][indice+0]-vecs[0][indice+0]; T reg1=vecs[1][indice+1]-vecs[0][indice+1]; reg0=abs(reg0); reg1=abs(reg1); T reg2=vecs[1][indice+2]-vecs[0][indice+2];
    reg0=max(reg1,reg0); reg2=abs(reg2); return max(reg2,reg0);
  }
  template<class TE,class TTs,class Tvecs,class Tvec>
  inline static void set_old_vec_nodal(const TE &node,const TTs &f,const Tvecs &vecs,Tvec &old_vec,int indice) {
    old_vec[indice+1]=vecs[1][indice+1]; old_vec[indice+0]=vecs[1][indice+0]; old_vec[indice+2]=vecs[1][indice+2];
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
#endif // ELASTICITY_ISOTROPY_STAT_QSTAT_3_NUM_0_CARAC_H
extern "C" void apply_on_elements_after_solve_0_elasticity_isotropy_stat_Qstat_Hexa_0( double * );
extern "C" void apply_on_elements_after_solve_1_elasticity_isotropy_stat_Qstat_Hexa_1( double * );
extern "C" void apply_on_elements_after_solve_2_elasticity_isotropy_stat_Qstat_Hexa_2( double * );
extern "C" void apply_on_elements_after_solve_3_elasticity_isotropy_stat_Qstat_Hexa_3( double * );
extern "C" void apply_on_elements_after_solve_4_elasticity_isotropy_stat_Qstat_Hexa_4( double * );
extern "C" void apply_on_elements_after_solve_5_elasticity_isotropy_stat_Qstat_Hexa_5( double * );
extern "C" void apply_on_elements_after_solve_6_elasticity_isotropy_stat_Qstat_Hexa_6( double * );
extern "C" void apply_on_elements_after_solve_7_elasticity_isotropy_stat_Qstat_Hexa_7( double * );
extern "C" void apply_on_elements_after_solve_8_elasticity_isotropy_stat_Qstat_Hexa_8( double * );
extern "C" void apply_on_elements_after_solve_9_elasticity_isotropy_stat_Qstat_Hexa_9( double * );
extern "C" void apply_on_elements_after_solve_10_elasticity_isotropy_stat_Qstat_Hexa_10( double * );
extern "C" void apply_on_elements_after_solve_11_elasticity_isotropy_stat_Qstat_Hexa_11( double * );
extern "C" void apply_on_elements_after_solve_12_elasticity_isotropy_stat_Qstat_Hexa_12( double * );
extern "C" void apply_on_elements_after_solve_13_elasticity_isotropy_stat_Qstat_Hexa_13( double * );
extern "C" void apply_on_elements_after_solve_14_elasticity_isotropy_stat_Qstat_Hexa_14( double * );
class Hexa;
template<unsigned A,class B,class C> class Node;
template<class A,class B,class C,class D,unsigned E> class Element;

// Carac for ...
template<class P_T_pos,class P_ND,class TED,unsigned nim,class TM,class T>
class CaracFormulationForElement<elasticity_isotropy_stat_Qstat,Element<Hexa,DefaultBehavior,Node<3,P_T_pos,P_ND>,TED,nim>,TM,T> {
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
    T reg0=0.622008467928146233*elem.pos(1)[2]; T reg1=0.16666666666666668806*elem.pos(0)[2]; T reg2=0.62200846792814627674*elem.pos(1)[2]; T reg3=0.622008467928146233*elem.pos(1)[1]; T reg4=0.62200846792814627674*elem.pos(1)[1];
    T reg5=0.16666666666666664427*elem.pos(1)[1]; T reg6=0.16666666666666668806*elem.pos(0)[1]; T reg7=0.62200846792814627674*elem.pos(0)[1]; T reg8=0.62200846792814627674*elem.pos(0)[2]; T reg9=0.16666666666666664427*elem.pos(1)[2];
    reg3=reg3+reg6; reg0=reg0+reg1; reg9=reg9+reg8; T reg10=0.16666666666666663255*elem.pos(2)[1]; T reg11=0.622008467928146233*elem.pos(2)[1];
    T reg12=0.16666666666666663255*elem.pos(2)[2]; T reg13=0.16666666666666667632*elem.pos(1)[1]; T reg14=0.044658198738520458147*elem.pos(0)[1]; T reg15=0.044658198738520434687*elem.pos(2)[2]; T reg16=0.044658198738520458147*elem.pos(0)[2];
    T reg17=0.16666666666666667632*elem.pos(1)[2]; T reg18=0.16666666666666664427*elem.pos(2)[2]; reg5=reg5+reg7; T reg19=0.622008467928146233*elem.pos(2)[2]; reg8=reg2-reg8;
    reg2=0.16666666666666668806*elem.pos(1)[2]; T reg20=0.16666666666666668806*elem.pos(1)[1]; T reg21=0.044658198738520434687*elem.pos(2)[1]; T reg22=0.16666666666666664427*elem.pos(2)[1]; reg7=reg4-reg7;
    reg4=0.044658198738520446417*elem.pos(3)[2]; T reg23=0.16666666666666668806*elem.pos(3)[1]; reg2=reg2-reg1; T reg24=0.62200846792814627674*elem.pos(3)[2]; reg13=reg13+reg14;
    T reg25=reg18-reg9; T reg26=0.6220084679281461892*elem.pos(2)[1]; T reg27=reg19-reg0; reg15=reg9+reg15; reg20=reg20-reg6;
    reg9=0.044658198738520446417*elem.pos(1)[1]; T reg28=0.16666666666666664427*elem.pos(3)[1]; reg21=reg5+reg21; T reg29=0.16666666666666668806*elem.pos(0)[0]; T reg30=0.622008467928146233*elem.pos(1)[0];
    reg5=reg22-reg5; T reg31=0.62200846792814627674*elem.pos(3)[1]; reg8=reg18+reg8; reg18=0.62200846792814627674*elem.pos(0)[0]; T reg32=0.62200846792814627674*elem.pos(1)[0];
    T reg33=0.044658198738520446417*elem.pos(3)[1]; reg7=reg22+reg7; reg22=reg11-reg3; reg0=reg0+reg12; T reg34=0.16666666666666664427*elem.pos(3)[2];
    T reg35=0.6220084679281461892*elem.pos(2)[2]; reg17=reg17+reg16; T reg36=0.16666666666666664427*elem.pos(1)[0]; T reg37=0.044658198738520446417*elem.pos(1)[2]; T reg38=0.16666666666666668806*elem.pos(3)[2];
    reg3=reg3+reg10; reg9=reg6+reg9; reg31=reg5+reg31; reg15=reg34+reg15; reg5=0.16666666666666664427*elem.pos(4)[1];
    reg6=0.16666666666666664427*elem.pos(2)[0]; reg27=reg38+reg27; T reg39=0.622008467928146233*elem.pos(3)[2]; reg37=reg1+reg37; reg32=reg32-reg18;
    reg1=0.62200846792814627674*elem.pos(4)[2]; reg19=reg19+reg2; reg0=reg4+reg0; reg30=reg30+reg29; T reg40=0.622008467928146233*elem.pos(2)[0];
    reg34=reg8-reg34; reg8=0.16666666666666668806*elem.pos(1)[0]; reg3=reg33+reg3; reg7=reg7-reg28; reg36=reg18+reg36;
    reg18=0.16666666666666667632*elem.pos(3)[2]; reg35=reg17+reg35; T reg41=0.16666666666666668806*elem.pos(4)[1]; reg26=reg13+reg26; T reg42=0.16666666666666667632*elem.pos(3)[1];
    T reg43=0.044658198738520446417*elem.pos(4)[2]; T reg44=0.16666666666666668806*elem.pos(4)[2]; reg25=reg24+reg25; reg24=0.16666666666666664427*elem.pos(4)[2]; T reg45=0.622008467928146233*elem.pos(3)[1];
    reg11=reg11+reg20; reg28=reg21+reg28; reg21=0.62200846792814627674*elem.pos(4)[1]; reg22=reg22+reg23; T reg46=0.044658198738520446417*elem.pos(4)[1];
    T reg47=0.044658198738520434687*elem.pos(2)[0]; T reg48=reg40-reg30; reg27=reg27-reg43; T reg49=0.16666666666666668806*elem.pos(3)[0]; reg34=reg34-reg24;
    reg28=reg21-reg28; reg31=reg31-reg5; reg21=0.044658198738520434687*elem.pos(5)[1]; reg24=reg25-reg24; reg25=0.622008467928146233*elem.pos(5)[1];
    reg3=reg41-reg3; T reg50=0.044658198738520434687*elem.pos(5)[2]; T reg51=0.16666666666666664427*elem.pos(5)[2]; reg5=reg7-reg5; reg15=reg1-reg15;
    reg1=0.16666666666666663255*elem.pos(2)[0]; reg7=0.044658198738520446417*elem.pos(2)[1]; T reg52=reg6-reg36; T reg53=0.16666666666666663255*elem.pos(5)[2]; T reg54=0.16666666666666664427*elem.pos(5)[1];
    T reg55=0.62200846792814627674*elem.pos(3)[0]; reg22=reg22-reg46; T reg56=0.16666666666666663255*elem.pos(5)[1]; reg11=reg11-reg45; reg26=reg26+reg42;
    T reg57=0.044658198738520458147*elem.pos(4)[1]; reg35=reg18+reg35; reg19=reg19-reg39; T reg58=0.044658198738520458147*elem.pos(4)[2]; reg10=reg10+reg9;
    T reg59=0.044658198738520446417*elem.pos(2)[2]; T reg60=0.16666666666666667632*elem.pos(1)[0]; T reg61=0.044658198738520458147*elem.pos(0)[0]; reg12=reg12+reg37; T reg62=0.622008467928146233*elem.pos(5)[2];
    reg0=reg44-reg0; T reg63=0.16666666666666664427*elem.pos(3)[0]; reg6=reg32+reg6; reg8=reg8-reg29; reg26=reg57-reg26;
    reg32=0.16666666666666667632*elem.pos(5)[1]; reg3=reg3+reg25; reg57=0.16666666666666667632*elem.pos(5)[2]; reg28=reg28+reg54; reg40=reg40+reg8;
    reg34=reg51+reg34; reg35=reg58-reg35; reg5=reg54+reg5; reg54=0.044658198738520458147*elem.pos(1)[2]; reg52=reg55+reg52;
    reg9=reg7-reg9; reg7=reg20+reg7; reg30=reg30+reg1; reg20=0.16666666666666667632*elem.pos(2)[1]; reg55=0.16666666666666667632*elem.pos(2)[2];
    reg58=0.044658198738520446417*elem.pos(3)[0]; reg47=reg36+reg47; reg36=0.16666666666666663255*elem.pos(6)[2]; reg0=reg0+reg62; reg2=reg2+reg59;
    T reg64=0.16666666666666664427*elem.pos(4)[0]; reg6=reg6-reg63; reg21=reg31-reg21; reg31=0.044658198738520434687*elem.pos(6)[1]; T reg65=0.6220084679281461892*elem.pos(2)[0];
    reg60=reg60+reg61; reg12=reg39+reg12; reg39=0.044658198738520458147*elem.pos(1)[1]; T reg66=0.044658198738520446417*elem.pos(5)[2]; reg10=reg45+reg10;
    reg37=reg59-reg37; reg51=reg15+reg51; reg15=0.044658198738520434687*elem.pos(6)[2]; reg43=reg19-reg43; reg48=reg49+reg48;
    reg19=0.044658198738520446417*elem.pos(4)[0]; reg45=0.622008467928146233*elem.pos(3)[0]; reg27=reg27-reg53; reg50=reg24-reg50; reg46=reg11-reg46;
    reg11=0.16666666666666663255*elem.pos(6)[1]; reg24=0.044658198738520446417*elem.pos(5)[1]; reg59=0.044658198738520446417*elem.pos(1)[0]; reg22=reg22-reg56; reg4=reg2-reg4;
    reg3=reg11+reg3; reg2=0.16666666666666664427*PNODE(1).dep[0]; T reg67=0.62200846792814627674*PNODE(1).dep[0]; T reg68=0.62200846792814627674*PNODE(0).dep[0]; reg33=reg7-reg33;
    reg13=reg20-reg13; reg7=0.044658198738520458147*elem.pos(3)[1]; reg16=reg54-reg16; reg54=0.62200846792814627674*PNODE(0).dep[1]; T reg69=0.16666666666666664427*PNODE(1).dep[1];
    reg27=reg36+reg27; reg22=reg22+reg11; T reg70=0.044658198738520446417*elem.pos(7)[1]; T reg71=0.622008467928146233*elem.pos(4)[2]; reg37=reg38+reg37;
    reg9=reg23+reg9; reg23=0.622008467928146233*elem.pos(4)[1]; reg38=0.62200846792814627674*PNODE(1).dep[1]; reg0=reg0+reg36; T reg72=0.044658198738520446417*elem.pos(7)[2];
    reg21=reg21+reg31; reg14=reg39-reg14; reg39=0.16666666666666664427*elem.pos(5)[0]; reg6=reg6-reg64; reg65=reg60+reg65;
    reg12=reg44-reg12; reg44=0.16666666666666667632*elem.pos(3)[0]; reg43=reg43+reg66; reg10=reg41-reg10; reg41=0.622008467928146233*PNODE(1).dep[0];
    T reg73=0.16666666666666668806*PNODE(0).dep[0]; reg46=reg46+reg24; reg59=reg29+reg59; reg29=0.6220084679281461892*elem.pos(6)[1]; reg26=reg26+reg32;
    T reg74=0.044658198738520446417*elem.pos(2)[0]; T reg75=0.6220084679281461892*elem.pos(6)[2]; reg35=reg35+reg57; reg17=reg55-reg17; T reg76=0.044658198738520458147*elem.pos(3)[2];
    reg40=reg40-reg45; reg47=reg63+reg47; reg63=0.62200846792814627674*elem.pos(4)[0]; reg48=reg48-reg19; T reg77=0.16666666666666663255*elem.pos(5)[0];
    T reg78=0.044658198738520434687*elem.pos(7)[2]; reg34=reg15+reg34; T reg79=0.044658198738520434687*elem.pos(7)[1]; reg5=reg31+reg5; T reg80=0.16666666666666668806*elem.pos(4)[0];
    T reg81=0.044658198738520434687*elem.pos(5)[0]; reg64=reg52-reg64; reg30=reg58+reg30; reg28=reg31+reg28; reg31=0.16666666666666664427*elem.pos(7)[1];
    reg51=reg51+reg15; reg52=0.16666666666666664427*elem.pos(7)[2]; T reg82=0.622008467928146233*PNODE(1).dep[1]; T reg83=0.16666666666666668806*PNODE(0).dep[1]; reg50=reg15+reg50;
    reg15=0.16666666666666667632*elem.pos(4)[2]; reg17=reg76+reg17; reg76=0.16666666666666663255*elem.pos(6)[0]; reg48=reg48-reg77; reg16=reg55+reg16;
    reg65=reg44+reg65; reg3=reg70+reg3; reg82=reg83+reg82; reg55=0.16666666666666668806*PNODE(1).dep[0]; T reg84=0.16666666666666667632*elem.pos(4)[1];
    reg14=reg20+reg14; reg0=reg0+reg72; reg67=reg67-reg68; reg7=reg13+reg7; reg19=reg40-reg19;
    reg13=0.16666666666666664427*PNODE(2).dep[0]; reg4=reg4-reg71; reg20=0.044658198738520446417*elem.pos(5)[0]; reg40=0.044658198738520458147*elem.pos(4)[0]; T reg85=0.16666666666666663255*elem.pos(7)[1];
    reg46=reg11+reg46; T reg86=0.622008467928146233*PNODE(2).dep[0]; T reg87=0.16666666666666668806*PNODE(1).dep[1]; reg70=reg22+reg70; reg22=0.622008467928146233*elem.pos(5)[0];
    reg41=reg41+reg73; reg27=reg72+reg27; reg30=reg80-reg30; reg72=0.044658198738520458147*elem.pos(1)[0]; T reg88=0.16666666666666667632*elem.pos(7)[1];
    reg26=reg26+reg29; T reg89=0.622008467928146233*PNODE(2).dep[1]; reg43=reg36+reg43; T reg90=0.16666666666666667632*elem.pos(7)[2]; reg35=reg35+reg75;
    T reg91=0.16666666666666663255*elem.pos(7)[2]; T reg92=0.16666666666666668806*PNODE(0).dep[2]; T reg93=0.622008467928146233*PNODE(1).dep[2]; T reg94=0.16666666666666667632*elem.pos(2)[0]; reg1=reg1+reg59;
    reg10=reg24+reg10; reg12=reg66+reg12; reg6=reg6+reg39; reg24=0.044658198738520434687*elem.pos(6)[0]; reg33=reg33-reg23;
    reg21=reg21+reg31; reg51=reg51+reg52; reg50=reg52+reg50; reg28=reg31+reg28; reg81=reg64-reg81;
    reg79=reg5-reg79; reg78=reg34-reg78; reg47=reg63-reg47; reg2=reg68+reg2; reg5=0.16666666666666664427*PNODE(2).dep[1];
    reg69=reg54+reg69; reg23=reg9-reg23; reg54=reg38-reg54; reg71=reg37-reg71; reg59=reg74-reg59;
    reg9=0.16666666666666664427*PNODE(1).dep[2]; reg31=0.62200846792814627674*PNODE(0).dep[2]; reg34=1+(*f.m).poisson_ratio; reg37=0.62200846792814627674*PNODE(1).dep[2]; reg74=reg8+reg74;
    reg8=0.62200846792814627674*PNODE(3).dep[1]; reg17=reg17-reg15; reg38=reg27*reg3; reg47=reg39+reg47; reg39=0.6220084679281461892*elem.pos(5)[2];
    reg52=reg28*reg78; reg63=0.044658198738520458147*elem.pos(3)[0]; reg64=reg51*reg79; reg60=reg94-reg60; reg48=reg48+reg76;
    reg66=0.044658198738520446417*elem.pos(7)[0]; reg68=reg79*reg0; T reg95=reg78*reg3; T reg96=0.16666666666666667632*PNODE(1).dep[1]; T reg97=0.044658198738520458147*PNODE(0).dep[1];
    reg18=reg16-reg18; reg9=reg9+reg31; reg59=reg49+reg59; reg16=0.044658198738520434687*PNODE(2).dep[1]; reg53=reg71-reg53;
    reg49=0.16666666666666664427*PNODE(3).dep[1]; reg54=reg5+reg54; reg7=reg7-reg84; reg56=reg23-reg56; reg23=reg70*reg0;
    reg71=0.6220084679281461892*elem.pos(5)[1]; reg5=reg5-reg69; reg19=reg19+reg20; reg93=reg93+reg92; reg34=reg34/(*f.m).elastic_modulus;
    reg12=reg36+reg12; T reg98=0.622008467928146233*elem.pos(7)[2]; reg43=reg43-reg91; T reg99=0.622008467928146233*PNODE(2).dep[2]; reg10=reg11+reg10;
    T reg100=0.622008467928146233*elem.pos(7)[1]; reg46=reg46-reg85; T reg101=0.16666666666666668806*PNODE(3).dep[0]; T reg102=reg86-reg41; reg58=reg74-reg58;
    reg1=reg45+reg1; reg45=0.16666666666666663255*PNODE(2).dep[0]; reg61=reg72-reg61; reg26=reg26+reg88; reg72=0.622008467928146233*elem.pos(4)[0];
    reg74=0.16666666666666664427*elem.pos(7)[0]; reg81=reg24+reg81; T reg103=reg50*reg28; reg30=reg30+reg22; reg87=reg87-reg83;
    T reg104=0.16666666666666668806*PNODE(3).dep[1]; T reg105=reg21*reg51; T reg106=reg89-reg82; T reg107=0.16666666666666663255*PNODE(2).dep[1]; reg42=reg14-reg42;
    reg14=0.16666666666666667632*elem.pos(5)[0]; reg65=reg40-reg65; reg40=0.044658198738520434687*elem.pos(7)[0]; reg6=reg6+reg24; reg35=reg35+reg90;
    reg31=reg37-reg31; reg37=0.16666666666666664427*PNODE(2).dep[2]; reg4=reg62+reg4; reg33=reg25+reg33; reg67=reg13+reg67;
    reg25=0.16666666666666664427*PNODE(3).dep[0]; reg62=0.044658198738520458147*PNODE(0).dep[0]; T reg108=0.16666666666666667632*PNODE(1).dep[0]; T reg109=0.62200846792814627674*PNODE(3).dep[0]; reg13=reg13-reg2;
    T reg110=0.044658198738520434687*PNODE(2).dep[0]; reg55=reg55-reg73; T reg111=0.25*elem.pos(0)[2]; T reg112=0.25*elem.pos(1)[2]; T reg113=0.16666666666666668806*PNODE(1).dep[2];
    T reg114=0.25*elem.pos(0)[1]; T reg115=0.25*elem.pos(1)[1]; reg89=reg89+reg87; reg106=reg104+reg106; T reg116=0.044658198738520446417*PNODE(4).dep[1];
    reg13=reg109+reg13; reg109=0.044658198738520446417*PNODE(3).dep[1]; T reg117=0.044658198738520458147*PNODE(0).dep[2]; reg82=reg82+reg107; T reg118=0.6220084679281461892*elem.pos(6)[0];
    reg84=reg42-reg84; reg40=reg6-reg40; reg65=reg65+reg14; reg6=0.6220084679281461892*PNODE(2).dep[0]; reg108=reg108+reg62;
    reg42=0.16666666666666663255*PNODE(2).dep[2]; reg58=reg58-reg72; T reg119=0.16666666666666664427*PNODE(4).dep[0]; reg67=reg67-reg25; T reg120=0.25*elem.pos(2)[2];
    T reg121=reg21*reg78; reg38=reg23-reg38; reg23=reg50*reg79; reg47=reg24+reg47; reg52=reg64-reg52;
    reg86=reg86+reg55; reg24=0.622008467928146233*PNODE(3).dep[0]; reg48=reg48+reg66; reg95=reg68-reg95; reg81=reg81+reg74;
    reg110=reg2+reg110; reg2=0.16666666666666667632*PNODE(1).dep[2]; reg103=reg105-reg103; reg30=reg76+reg30; reg64=reg79*reg27;
    reg68=reg78*reg70; reg105=reg70*reg35; T reg122=0.62200846792814627674*PNODE(3).dep[2]; T reg123=reg37-reg9; T reg124=0.16666666666666667632*elem.pos(4)[0];
    reg60=reg63+reg60; reg63=0.16666666666666663255*elem.pos(7)[0]; reg19=reg76+reg19; reg39=reg17-reg39; reg4=reg36+reg4;
    reg17=0.044658198738520446417*PNODE(1).dep[0]; reg71=reg7-reg71; reg37=reg31+reg37; reg7=0.25*elem.pos(2)[1]; reg12=reg12+reg98;
    reg31=reg26*reg43; T reg125=reg115-reg114; T reg126=reg99-reg93; T reg127=0.16666666666666668806*PNODE(3).dep[2]; reg33=reg11+reg33;
    reg10=reg10+reg100; T reg128=reg112-reg111; T reg129=reg35*reg46; reg102=reg101+reg102; T reg130=0.044658198738520446417*PNODE(4).dep[0];
    reg1=reg80-reg1; reg80=0.044658198738520446417*PNODE(3).dep[0]; reg41=reg41+reg45; reg61=reg94+reg61; reg94=reg27*reg26;
    T reg131=0.044658198738520446417*PNODE(1).dep[1]; T reg132=pow(reg34,2); T reg133=0.16666666666666664427*PNODE(3).dep[2]; reg53=reg36+reg53; reg56=reg11+reg56;
    reg16=reg69+reg16; reg11=0.044658198738520434687*PNODE(2).dep[2]; reg36=0.16666666666666664427*PNODE(4).dep[1]; reg113=reg113-reg92; reg5=reg8+reg5;
    reg72=reg59-reg72; reg114=reg115+reg114; reg54=reg54-reg49; reg8=0.622008467928146233*PNODE(3).dep[1]; reg15=reg18-reg15;
    reg96=reg97+reg96; reg18=0.6220084679281461892*PNODE(2).dep[1]; reg111=reg112+reg111; reg125=reg7+reg125; reg77=reg72-reg77;
    reg59=reg21*reg12; reg126=reg126+reg127; reg69=0.044658198738520446417*PNODE(4).dep[2]; reg72=reg50*reg10; reg85=reg33-reg85;
    reg128=reg120+reg128; reg33=0.16666666666666664427*PNODE(5).dep[0]; reg112=reg46*reg12; reg115=reg43*reg10; T reg134=0.16666666666666668806*PNODE(4).dep[1];
    T reg135=0.25*elem.pos(0)[0]; reg99=reg99+reg113; reg82=reg109+reg82; T reg136=0.16666666666666667632*elem.pos(7)[0]; reg65=reg65+reg118;
    T reg137=0.62200846792814627674*PNODE(4).dep[1]; reg2=reg2+reg117; reg16=reg49+reg16; reg49=0.6220084679281461892*PNODE(2).dep[2]; reg6=reg6+reg108;
    reg93=reg93+reg42; T reg138=0.044658198738520446417*PNODE(3).dep[2]; reg58=reg22+reg58; reg53=reg98+reg53; reg67=reg67-reg119;
    reg22=0.25*elem.pos(1)[0]; reg31=reg129-reg31; reg98=0.6220084679281461892*elem.pos(5)[0]; reg60=reg60-reg124; reg11=reg9+reg11;
    reg9=reg7+reg114; reg91=reg4-reg91; reg19=reg19-reg63; reg39=reg75+reg39; reg4=0.16666666666666667632*PNODE(3).dep[1];
    reg17=reg73+reg17; reg71=reg29+reg71; reg73=0.16666666666666664427*PNODE(4).dep[2]; reg37=reg37-reg133; reg129=0.044658198738520446417*PNODE(1).dep[2];
    reg18=reg96+reg18; reg102=reg102-reg130; T reg139=0.16666666666666663255*PNODE(5).dep[0]; reg1=reg20+reg1; reg20=0.16666666666666668806*PNODE(4).dep[0];
    reg89=reg89-reg8; reg41=reg80+reg41; reg44=reg61-reg44; reg61=0.25*elem.pos(3)[1]; reg94=reg105-reg94;
    reg105=0.044658198738520446417*PNODE(2).dep[1]; reg114=reg7-reg114; reg34=reg34*reg132; reg131=reg83+reg131; reg15=reg57+reg15;
    reg7=(*f.m).poisson_ratio/(*f.m).elastic_modulus; reg122=reg123+reg122; reg57=1.0/(*f.m).elastic_modulus; reg5=reg5-reg36; reg83=0.16666666666666663255*PNODE(5).dep[1];
    reg47=reg74+reg47; reg74=0.044658198738520434687*PNODE(5).dep[1]; reg106=reg106-reg116; reg36=reg54-reg36; reg119=reg13-reg119;
    reg13=reg81*reg52; reg54=0.044658198738520434687*PNODE(5).dep[0]; reg123=0.16666666666666667632*PNODE(3).dep[0]; reg86=reg86-reg24; reg110=reg25+reg110;
    reg25=reg48*reg95; reg68=reg64-reg68; reg64=0.622008467928146233*PNODE(3).dep[2]; T reg140=0.25*elem.pos(3)[2]; T reg141=0.62200846792814627674*PNODE(4).dep[0];
    reg30=reg66+reg30; reg66=reg40*reg103; T reg142=reg120-reg111; T reg143=0.16666666666666664427*PNODE(5).dep[1]; reg84=reg32+reg84;
    reg32=reg40*reg38; T reg144=0.044658198738520446417*PNODE(2).dep[0]; reg120=reg111+reg120; reg121=reg23-reg121; reg23=reg70*reg43;
    reg111=reg27*reg46; reg56=reg100+reg56; reg25=reg32-reg25; reg32=reg50*reg46; reg100=reg21*reg43;
    T reg145=0.16666666666666667632*PNODE(2).dep[0]; T reg146=reg105-reg131; T reg147=reg78*reg47; T reg148=reg7*reg34; T reg149=reg51*reg81;
    reg34=reg57*reg34; T reg150=reg40*reg51; reg110=reg141-reg110; reg120=reg140+reg120; reg141=reg50*reg47;
    reg105=reg87+reg105; reg124=reg44-reg124; reg44=reg19*reg94; reg87=0.622008467928146233*elem.pos(7)[0]; reg1=reg76+reg1;
    reg41=reg20-reg41; reg129=reg92+reg129; reg114=reg114+reg61; reg18=reg4+reg18; reg92=0.044658198738520458147*PNODE(1).dep[1];
    reg39=reg90+reg39; reg90=0.044658198738520458147*PNODE(4).dep[1]; T reg151=reg51*reg85; reg133=reg11+reg133; reg11=0.044658198738520446417*PNODE(5).dep[0];
    reg13=reg66-reg13; reg66=reg47*reg121; reg74=reg5-reg74; reg98=reg60-reg98; reg5=reg28*reg91;
    reg45=reg45+reg17; reg37=reg37-reg73; reg116=reg89-reg116; reg60=0.6220084679281461892*elem.pos(7)[2]; reg15=reg75+reg15;
    reg71=reg88+reg71; reg75=0.044658198738520434687*PNODE(6).dep[1]; reg88=0.16666666666666667632*PNODE(2).dep[1]; reg17=reg144-reg17; reg58=reg76+reg58;
    reg73=reg122-reg73; reg130=reg86-reg130; reg86=0.044658198738520434687*PNODE(5).dep[2]; reg131=reg107+reg131; reg9=reg61+reg9;
    reg77=reg76+reg77; reg61=reg125-reg61; reg76=0.16666666666666663255*PNODE(5).dep[2]; reg126=reg126-reg69; reg89=0.044658198738520446417*PNODE(5).dep[1];
    reg107=0.044658198738520458147*PNODE(4).dep[0]; reg122=reg27*reg30; reg54=reg119-reg54; reg144=reg55+reg144; reg55=0.044658198738520458147*PNODE(1).dep[0];
    reg119=0.16666666666666667632*PNODE(3).dep[2]; reg125=reg0*reg48; reg49=reg2+reg49; reg16=reg137-reg16; reg137=0.16666666666666668806*PNODE(4).dep[2];
    T reg152=reg22-reg135; T reg153=reg48*reg31; T reg154=0.044658198738520434687*PNODE(6).dep[0]; reg67=reg33+reg67; T reg155=reg28*reg53;
    reg135=reg22+reg135; reg22=0.16666666666666663255*PNODE(6).dep[1]; reg36=reg143+reg36; T reg156=0.25*elem.pos(2)[0]; reg93=reg93+reg138;
    reg106=reg106-reg83; T reg157=0.622008467928146233*PNODE(5).dep[1]; T reg158=0.622008467928146233*PNODE(5).dep[0]; reg65=reg65+reg136; T reg159=0.25*elem.pos(4)[2];
    reg142=reg140+reg142; reg115=reg112-reg115; reg112=0.16666666666666663255*PNODE(6).dep[0]; reg102=reg102-reg139; reg99=reg99-reg64;
    T reg160=0.6220084679281461892*elem.pos(7)[1]; reg84=reg29+reg84; reg29=reg51*reg56; reg6=reg123+reg6; reg23=reg111-reg23;
    reg82=reg134-reg82; reg111=0.16666666666666664427*PNODE(5).dep[2]; T reg161=0.044658198738520446417*PNODE(2).dep[2]; T reg162=0.62200846792814627674*PNODE(4).dep[2]; T reg163=0.25*elem.pos(4)[1];
    T reg164=reg78*reg30; T reg165=reg40*reg0; reg140=reg128-reg140; reg128=reg30*reg68; reg72=reg59-reg72;
    reg59=reg43*reg65; reg80=reg144-reg80; reg144=reg27*reg65; reg63=reg58-reg63; reg18=reg90-reg18;
    reg58=reg48*reg35; reg113=reg113+reg161; reg90=reg19*reg35; reg42=reg42+reg129; reg69=reg99-reg69;
    reg99=reg65*reg23; reg130=reg11+reg130; reg109=reg105-reg109; reg6=reg107-reg6; reg105=reg81*reg115;
    reg1=reg1+reg87; reg107=reg19*reg72; T reg166=0.622008467928146233*PNODE(4).dep[1]; reg100=reg32-reg100; reg146=reg104+reg146;
    reg17=reg101+reg17; reg32=0.622008467928146233*PNODE(4).dep[0]; reg131=reg8+reg131; reg8=0.16666666666666667632*PNODE(5).dep[0]; reg49=reg49+reg119;
    reg45=reg24+reg45; reg24=0.044658198738520458147*PNODE(4).dep[2]; reg101=reg56*reg91; reg104=reg53*reg85; reg5=reg151-reg5;
    reg77=reg87+reg77; reg116=reg89+reg116; reg87=0.044658198738520446417*PNODE(5).dep[2]; reg153=reg44-reg153; reg155=reg29-reg155;
    reg129=reg161-reg129; reg29=0.16666666666666667632*PNODE(5).dep[1]; reg37=reg111+reg37; reg126=reg126-reg76; reg62=reg55-reg62;
    reg141=reg149-reg141; reg114=reg114-reg163; reg44=0.25*elem.pos(5)[1]; reg55=reg7*reg34; reg106=reg106+reg22;
    reg149=reg28*reg81; reg36=reg75+reg36; reg151=0.044658198738520434687*PNODE(7).dep[1]; reg102=reg102+reg112; reg161=0.044658198738520446417*PNODE(7).dep[0];
    reg120=reg159-reg120; T reg167=0.044658198738520446417*PNODE(7).dep[1]; reg124=reg14+reg124; reg67=reg67+reg154; reg41=reg158+reg41;
    reg54=reg154+reg54; reg14=reg7*reg148; T reg168=reg78*reg48; T reg169=reg40*reg27; reg82=reg157+reg82;
    T reg170=0.25*elem.pos(3)[0]; reg152=reg152+reg156; reg66=reg13+reg66; reg147=reg150-reg147; reg13=0.622008467928146233*PNODE(5).dep[2];
    reg93=reg137-reg93; reg150=reg156-reg135; T reg171=0.16666666666666663255*PNODE(6).dep[2]; T reg172=0.16666666666666664427*PNODE(7).dep[1]; reg74=reg74+reg75;
    T reg173=reg3*reg48; T reg174=reg70*reg30; reg122=reg125-reg122; reg125=0.044658198738520434687*PNODE(7).dep[0]; T reg175=reg40*reg3;
    T reg176=reg79*reg30; T reg177=reg21*reg47; reg96=reg88-reg96; T reg178=0.044658198738520458147*PNODE(3).dep[1]; T reg179=0.044658198738520434687*PNODE(6).dep[2];
    reg86=reg73-reg86; reg9=reg163-reg9; reg110=reg33+reg110; reg33=0.044658198738520458147*PNODE(1).dep[2]; reg163=reg61-reg163;
    reg61=reg3*reg91; reg73=reg0*reg85; reg140=reg140-reg159; reg98=reg118+reg98; reg133=reg162-reg133;
    reg162=reg3*reg39; T reg180=reg0*reg71; reg97=reg92-reg97; reg92=0.25*elem.pos(5)[2]; T reg181=reg40*reg50;
    T reg182=reg40*reg28; reg164=reg165-reg164; reg128=reg25+reg128; reg108=reg145-reg108; reg25=0.044658198738520458147*PNODE(3).dep[0];
    reg165=reg79*reg47; reg78=reg81*reg78; reg34=reg57*reg34; reg160=reg84-reg160; reg159=reg142-reg159;
    reg16=reg143+reg16; reg84=0.16666666666666664427*PNODE(7).dep[0]; reg142=0.16666666666666667632*PNODE(2).dep[2]; reg60=reg15-reg60; reg62=reg145+reg62;
    reg6=reg8+reg6; reg17=reg17-reg32; reg150=reg170+reg150; reg97=reg88+reg97; reg18=reg29+reg18;
    reg15=0.16666666666666667632*PNODE(5).dep[2]; reg117=reg33-reg117; reg33=0.6220084679281461892*PNODE(6).dep[1]; reg103=reg103/reg66; reg138=reg113-reg138;
    reg88=0.622008467928146233*PNODE(4).dep[2]; reg9=reg44+reg9; reg113=0.16666666666666663255*PNODE(7).dep[0]; reg130=reg112+reg130; reg133=reg111+reg133;
    reg110=reg154+reg110; reg111=reg81*reg79; reg143=reg40*reg21; reg145=reg63*reg155; reg154=0.25*elem.pos(6)[2];
    reg54=reg54+reg84; reg16=reg75+reg16; reg120=reg120+reg92; reg129=reg127+reg129; reg78=reg181-reg78;
    reg125=reg67-reg125; reg67=0.25*elem.pos(4)[0]; reg152=reg152-reg170; reg49=reg24-reg49; reg147=reg147/reg66;
    reg32=reg80-reg32; reg74=reg172+reg74; reg52=reg52/reg66; reg24=reg48*reg26; reg75=reg70*reg65;
    reg114=reg114-reg44; reg80=0.25*elem.pos(6)[1]; reg141=reg141/reg66; reg69=reg87+reg69; reg151=reg36-reg151;
    reg36=reg19*reg26; reg127=reg46*reg65; reg45=reg20-reg45; reg108=reg25+reg108; reg20=0.16666666666666667632*PNODE(4).dep[0];
    reg25=0.044658198738520434687*PNODE(7).dep[2]; reg41=reg112+reg41; reg37=reg179+reg37; reg181=reg50*reg1; reg102=reg102+reg161;
    reg95=reg95/reg128; T reg183=reg81*reg12; reg38=reg38/reg128; reg124=reg118+reg124; reg118=0.6220084679281461892*elem.pos(7)[0];
    T reg184=reg35*reg71; T reg185=reg26*reg39; reg34=reg34-reg14; reg126=reg171+reg126; reg162=reg180-reg162;
    reg101=reg104-reg101; reg98=reg136+reg98; reg61=reg73-reg61; reg42=reg64+reg42; reg64=reg85*reg39;
    reg73=reg91*reg71; reg165=reg182-reg165; reg86=reg179+reg86; reg96=reg178+reg96; reg104=0.16666666666666667632*PNODE(4).dep[1];
    reg136=reg77*reg5; reg131=reg134-reg131; reg2=reg142-reg2; reg134=0.044658198738520458147*PNODE(3).dep[2]; reg146=reg146-reg166;
    reg122=reg122/reg128; reg59=reg90-reg59; reg177=reg149-reg177; reg106=reg167+reg106; reg144=reg58-reg144;
    reg164=reg164/reg128; reg58=reg35*reg160; reg90=reg26*reg60; reg116=reg22+reg116; reg149=0.16666666666666663255*PNODE(7).dep[1];
    reg105=reg107-reg105; reg135=reg156+reg135; reg27=reg27*reg19; reg107=reg48*reg43; reg156=0.16666666666666664427*PNODE(7).dep[2];
    reg140=reg92+reg140; reg178=0.6220084679281461892*PNODE(6).dep[0]; reg163=reg44+reg163; reg176=reg175-reg176; reg44=reg43*reg1;
    reg175=reg19*reg12; reg174=reg173-reg174; reg92=reg159-reg92; reg159=0.044658198738520446417*PNODE(7).dep[2]; reg173=reg1*reg100;
    reg93=reg13+reg93; reg166=reg109-reg166; reg79=reg79*reg48; reg40=reg40*reg70; reg55=reg14+reg55;
    reg82=reg22+reg82; reg168=reg169-reg168; reg148=reg57*reg148; reg99=reg153+reg99; reg136=reg145-reg136;
    reg177=reg177/reg66; reg179=reg133+reg179; reg165=reg165/reg66; reg25=reg37-reg25; reg37=reg47*reg101;
    reg109=0.25*vectors[0][indices[0]+0]; reg133=0.25*vectors[0][indices[0]+1]; reg173=reg105+reg173; reg44=reg175-reg44; reg181=reg183-reg181;
    reg105=reg57*reg34; reg50=reg50*reg19; reg43=reg81*reg43; reg45=reg11+reg45; reg11=reg7*reg55;
    reg131=reg89+reg131; reg89=0.25*vectors[0][indices[1]+1]; reg145=reg46*reg1; reg153=reg19*reg10; reg169=0.25*vectors[0][indices[1]+0];
    reg42=reg137-reg42; reg137=reg21*reg1; reg175=reg81*reg10; reg111=reg143-reg111; reg9=reg80+reg9;
    reg16=reg172+reg16; reg78=reg78/reg66; reg143=reg141*reg151; reg172=reg147*reg74; reg152=reg152-reg67;
    reg180=0.25*elem.pos(5)[0]; reg150=reg150-reg67; reg163=reg80+reg163; reg140=reg154+reg140; reg182=reg57*reg132;
    reg135=reg170+reg135; reg132=reg7*reg132; reg148=reg14+reg148; reg14=reg151*reg122; reg59=reg59/reg99;
    reg170=reg164*reg106; reg144=reg144/reg99; reg116=reg116-reg149; reg90=reg58-reg90; reg107=reg27-reg107;
    reg27=reg39*reg160; reg58=0.16666666666666667632*PNODE(7).dep[1]; reg183=reg71*reg60; reg4=reg97-reg4; reg117=reg142+reg117;
    reg18=reg18+reg33; reg139=reg17-reg139; reg70=reg70*reg19; reg48=reg48*reg46; reg123=reg62-reg123;
    reg17=reg12*reg56; reg62=reg10*reg53; reg97=reg12*reg160; reg142=reg10*reg60; reg49=reg15+reg49;
    reg138=reg138-reg88; T reg186=reg63*reg162; T reg187=reg98*reg61; reg73=reg64-reg73; reg96=reg96-reg104;
    reg64=0.6220084679281461892*PNODE(5).dep[1]; reg134=reg2+reg134; reg2=0.16666666666666667632*PNODE(4).dep[2]; reg108=reg108-reg20; T reg188=0.6220084679281461892*PNODE(5).dep[0];
    reg41=reg161+reg41; reg68=reg68/reg128; reg161=reg95*reg102; T reg189=reg125*reg38; reg118=reg124-reg118;
    reg126=reg159+reg126; reg176=reg176/reg128; reg185=reg184-reg185; reg174=reg174/reg128; reg93=reg93+reg171;
    reg79=reg40-reg79; reg82=reg167+reg82; reg168=reg168/reg128; reg127=reg36-reg127; reg88=reg129-reg88;
    reg86=reg156+reg86; reg94=reg94/reg99; reg83=reg146-reg83; reg120=reg120+reg154; reg36=0.25*elem.pos(7)[2];
    reg110=reg84+reg110; reg121=reg121/reg66; reg130=reg130-reg113; reg31=reg31/reg99; reg40=reg51*reg77;
    reg84=reg47*reg53; reg124=reg52*reg54; reg129=reg103*reg125; reg6=reg6+reg178; reg146=0.16666666666666667632*PNODE(7).dep[0];
    reg92=reg154+reg92; reg166=reg157+reg166; reg154=0.25*elem.pos(7)[1]; reg69=reg171+reg69; reg80=reg114+reg80;
    reg75=reg24-reg75; reg24=reg47*reg91; reg114=0.16666666666666663255*PNODE(7).dep[2]; reg51=reg51*reg63; reg32=reg158+reg32;
    reg157=0.6220084679281461892*PNODE(6).dep[2]; reg137=reg175-reg137; reg188=reg108-reg188; reg108=0.622008467928146233*PNODE(7).dep[0]; reg37=reg136+reg37;
    reg84=reg40-reg84; reg187=reg186-reg187; reg45=reg112+reg45; reg40=0.25*vectors[0][indices[0]+2]; reg136=0.6220084679281461892*PNODE(5).dep[2];
    reg158=0.25*vectors[0][indices[1]+2]; reg134=reg134-reg2; reg83=reg22+reg83; reg24=reg51-reg24; reg51=reg30*reg73;
    reg167=reg30*reg39; reg175=reg0*reg98; reg0=reg0*reg63; reg64=reg96-reg64; reg92=reg36+reg92;
    reg96=reg30*reg91; reg184=reg133+reg89; reg145=reg153-reg145; reg133=reg89-reg133; reg72=reg72/reg173;
    reg115=reg115/reg173; reg142=reg97-reg142; reg89=reg53*reg160; reg97=reg56*reg60; reg62=reg17-reg62;
    reg20=reg123-reg20; reg152=reg152+reg180; reg17=0.25*elem.pos(6)[0]; reg32=reg112+reg32; reg80=reg80+reg154;
    reg150=reg150-reg180; reg139=reg112+reg139; reg76=reg88-reg76; reg119=reg117-reg119; reg104=reg4-reg104;
    reg183=reg27-reg183; reg4=reg47*reg85; reg166=reg22+reg166; reg27=reg63*reg53; reg88=reg77*reg91;
    reg112=reg118*reg185; reg117=reg28*reg77; reg123=0.25*vectors[0][indices[2]+1]; reg153=reg169+reg109; reg109=reg169-reg109;
    reg169=0.25*vectors[0][indices[2]+0]; reg120=reg120+reg36; reg135=reg67-reg135; reg9=reg154+reg9; reg67=reg98*reg90;
    reg36=reg140-reg36; reg138=reg13+reg138; reg28=reg28*reg63; reg154=reg163-reg154; reg6=reg6+reg146;
    reg13=reg57*reg182; reg140=reg7*reg132; reg182=reg7*reg182; reg14=reg170-reg14; reg163=reg168*reg82;
    reg79=reg79/reg128; reg159=reg93+reg159; reg93=reg25*reg174; reg170=reg176*reg126; reg161=reg189-reg161;
    reg186=reg68*reg41; reg107=reg107/reg99; reg44=reg44/reg173; reg181=reg181/reg173; reg124=reg129-reg124;
    reg129=reg7*reg148; reg11=reg105-reg11; reg105=reg121*reg110; reg189=reg165*reg86; T reg190=reg177*reg25;
    reg156=reg179+reg156; reg111=reg111/reg66; reg179=reg78*reg16; reg143=reg172-reg143; reg47=reg47*reg56;
    reg172=reg106*reg59; T reg191=reg144*reg116; reg19=reg21*reg19; reg127=reg127/reg99; reg48=reg70-reg48;
    reg18=reg58+reg18; reg46=reg81*reg46; reg21=reg94*reg130; reg70=reg102*reg31; reg131=reg22+reg131;
    reg22=0.622008467928146233*PNODE(7).dep[1]; reg43=reg50-reg43; reg69=reg69-reg114; reg42=reg87+reg42; reg50=0.16666666666666667632*PNODE(7).dep[2];
    reg23=reg23/reg99; reg49=reg49+reg157; reg75=reg75/reg99; reg81=reg126*reg127; reg49=reg49+reg50;
    reg87=reg75*reg69; reg76=reg171+reg76; T reg192=reg92*reg9; reg179=reg143-reg179; reg133=reg123+reg133;
    reg91=reg91*reg98; reg5=reg5/reg37; reg143=reg111*reg156; T reg193=reg63*reg39; reg113=reg32-reg113;
    reg167=reg175-reg167; reg32=reg164*reg102; reg51=reg187+reg51; reg96=reg0-reg96; reg0=reg151*reg38;
    reg175=reg125*reg122; reg187=0.25*vectors[0][indices[3]+1]; reg48=reg48/reg99; T reg194=reg123-reg184; reg155=reg155/reg37;
    reg64=reg33+reg64; T reg195=reg106*reg95; T reg196=reg77*reg85; T reg197=reg63*reg56; reg135=reg180+reg135;
    reg88=reg27-reg88; reg129=reg11-reg129; reg11=reg23*reg6; reg109=reg169+reg109; reg149=reg166-reg149;
    reg27=reg9*reg36; reg132=reg57*reg132; reg47=reg117-reg47; reg84=reg84/reg37; reg13=reg13-reg140;
    reg182=reg140+reg182; reg117=reg120*reg154; reg83=reg22+reg83; reg166=(*f.m).deltaT*(*f.m).alpha; reg180=0.25*vectors[0][indices[3]+0];
    T reg198=reg169-reg153; reg70=reg21-reg70; reg150=reg17+reg150; reg138=reg171+reg138; reg24=reg24/reg37;
    reg139=reg108+reg139; reg4=reg28-reg4; reg21=0.25*elem.pos(7)[0]; reg152=reg152+reg17; reg28=reg107*reg18;
    reg191=reg172-reg191; reg188=reg178+reg188; reg172=reg158+reg40; T reg199=reg54*reg115; T reg200=0.25*vectors[0][indices[2]+2];
    T reg201=reg147*reg54; T reg202=reg141*reg125; reg100=reg100/reg173; T reg203=reg151*reg103; T reg204=reg74*reg52;
    reg131=reg22+reg131; reg43=reg43/reg173; reg108=reg45+reg108; reg104=reg29+reg104; reg186=reg161+reg186;
    reg170=reg93-reg170; reg22=reg79*reg159; reg29=reg118*reg62; reg163=reg14-reg163; reg14=reg65*reg39;
    reg40=reg158-reg40; reg45=reg35*reg98; reg97=reg89-reg97; reg89=reg65*reg60; reg35=reg35*reg118;
    reg93=reg77*reg142; reg67=reg112-reg67; reg112=reg65*reg183; reg158=reg74*reg44; reg161=reg116*reg181;
    reg189=reg190-reg189; reg190=0.622008467928146233*PNODE(7).dep[2]; reg42=reg171+reg42; reg171=reg3*reg98; T reg205=reg30*reg71;
    reg137=reg137/reg173; reg3=reg3*reg63; reg20=reg8+reg20; reg30=reg30*reg85; reg46=reg19-reg46;
    reg145=reg145/reg173; reg136=reg134-reg136; reg8=reg130*reg72; reg2=reg119-reg2; reg105=reg124+reg105;
    reg19=reg80*reg120; reg13=reg57*reg13; reg57=0.6220084679281461892*PNODE(7).dep[0]; reg119=reg65*reg71; reg20=reg178+reg20;
    reg2=reg15+reg2; reg89=reg35-reg89; reg15=reg26*reg98; reg35=reg5*reg139; reg101=reg101/reg37;
    reg26=reg26*reg118; reg14=reg45-reg14; reg45=0.25*vectors[0][indices[4]+0]; reg135=reg17+reg135; reg103=reg25*reg103;
    reg17=reg92*reg154; reg124=reg80*reg36; reg134=reg177*reg125; reg178=reg165*reg54; reg65=reg65*reg160;
    reg34=reg34/reg129; reg104=reg33+reg104; reg33=reg155*reg113; T reg206=0.6220084679281461892*PNODE(7).dep[1]; reg109=reg109-reg180;
    T reg207=reg98*reg60; reg52=reg86*reg52; reg39=reg39*reg118; reg28=reg191-reg28; reg192=reg19-reg192;
    reg179=reg179-reg166; reg95=reg126*reg95; reg38=reg25*reg38; reg19=reg176*reg102; reg125=reg125*reg174;
    reg188=reg146+reg188; reg105=reg105-reg166; reg61=reg61/reg51; reg162=reg162/reg51; reg96=reg96/reg51;
    reg136=reg157+reg136; reg146=reg82*reg68; reg195=reg0-reg195; reg189=reg143+reg189; reg30=reg3-reg30;
    reg0=reg168*reg41; reg175=reg32-reg175; reg64=reg58+reg64; reg205=reg171-reg205; reg85=reg85*reg98;
    reg167=reg167/reg51; reg63=reg63*reg71; reg91=reg193-reg91; reg27=reg117-reg27; reg112=reg67+reg112;
    reg93=reg29-reg93; reg182=reg7*reg182; reg3=reg1*reg97; reg29=reg12*reg118; reg32=reg140+reg132;
    reg163=reg163-reg166; reg150=reg21+reg150; reg170=reg22+reg170; reg186=reg186-reg166; reg22=reg1*reg60;
    reg12=reg12*reg77; reg152=reg152-reg21; reg58=reg1*reg53; reg194=reg187+reg194; reg67=0.25*vectors[0][indices[4]+1];
    reg117=reg16*reg121; reg204=reg203-reg204; reg143=reg78*reg110; reg202=reg201-reg202; reg199=reg8-reg199;
    reg47=reg47/reg37; reg196=reg197-reg196; reg8=reg43*reg131; reg184=reg123+reg184; reg161=reg158-reg161;
    reg123=reg100*reg108; reg158=reg200-reg172; reg171=reg106*reg31; reg191=reg116*reg94; reg193=reg144*reg130;
    reg40=reg200+reg40; reg88=reg88/reg37; reg197=reg102*reg59; reg201=reg84*reg149; reg11=reg70+reg11;
    reg153=reg169+reg153; reg198=reg180+reg198; reg70=0.25*vectors[0][indices[3]+2]; reg55=reg55/reg129; reg169=reg24*reg83;
    reg203=reg48*reg49; reg42=reg42+reg190; T reg208=reg69*reg137; reg76=reg190+reg76; reg81=reg87-reg81;
    reg4=reg4/reg37; reg46=reg46/reg173; reg87=reg86*reg145; reg114=reg138-reg114; reg133=reg133-reg187;
    reg95=reg38-reg95; reg38=reg149*reg167; reg138=reg55*reg105; reg190=0.25*vectors[0][indices[4]+2]; reg68=reg159*reg68;
    T reg209=reg61*reg188; reg73=reg73/reg51; T reg210=reg34*reg179; T reg211=reg34*reg105; T reg212=reg55*reg179;
    reg200=reg172+reg200; reg102=reg102*reg127; reg172=reg75*reg130; reg143=reg202-reg143; reg117=reg204+reg117;
    reg202=0.25*vectors[0][indices[5]+1]; reg0=reg175-reg0; reg146=reg195+reg146; reg31=reg126*reg31; reg94=reg69*reg94;
    reg184=reg187+reg184; reg30=reg30/reg51; reg189=reg189-reg166; reg136=reg50+reg136; reg19=reg125-reg19;
    reg205=reg205/reg51; reg85=reg63-reg85; reg50=reg113*reg162; reg63=reg55*reg163; reg125=reg34*reg186;
    reg175=reg79*reg41; reg187=reg96*reg64; reg133=reg133-reg67; reg91=reg91/reg51; reg89=reg89/reg112;
    reg158=reg70+reg158; reg14=reg14/reg112; reg123=reg199+reg123; reg206=reg104-reg206; reg8=reg161-reg8;
    reg207=reg39-reg207; reg71=reg71*reg118; reg98=reg98*reg160; reg119=reg15-reg119; reg121=reg156*reg121;
    reg52=reg103-reg52; reg2=reg157+reg2; reg15=0.6220084679281461892*PNODE(7).dep[2]; reg87=reg208-reg87; reg65=reg26-reg65;
    reg185=reg185/reg112; reg26=reg111*reg110; reg178=reg134-reg178; reg39=reg46*reg42; reg57=reg20-reg57;
    reg90=reg90/reg112; reg194=reg194-reg67; reg20=reg1*reg160; reg103=reg10*reg118; reg1=reg1*reg56;
    reg10=reg10*reg77; reg60=reg77*reg60; reg53=reg53*reg118; reg58=reg12-reg58; reg22=reg29-reg22;
    reg12=reg34*reg163; reg29=reg55*reg186; reg170=reg170-reg166; reg104=reg74*reg115; reg134=reg116*reg72;
    reg157=reg130*reg181; reg161=reg54*reg44; reg3=reg93+reg3; reg40=reg40-reg70; reg141=reg141*reg25;
    reg147=reg147*reg86; reg165=reg74*reg165; reg177=reg151*reg177; reg35=reg33-reg35; reg33=reg47*reg114;
    reg93=reg110*reg101; reg195=reg18*reg23; reg109=reg109-reg45; reg171=reg191-reg171; reg196=reg196/reg37;
    reg191=reg107*reg6; reg193=reg197-reg193; reg197=reg192*reg152; reg148=reg148/reg129; reg153=reg180+reg153;
    reg180=0.25*vectors[0][indices[5]+0]; reg124=reg17-reg124; reg198=reg198-reg45; reg32=reg7*reg32; reg182=reg13-reg182;
    reg135=reg21+reg135; reg13=reg16*reg88; reg201=reg169-reg201; reg11=reg11-reg166; reg17=reg150*reg27;
    reg81=reg203+reg81; reg28=reg28-reg166; reg164=reg164*reg126; reg174=reg151*reg174; reg21=reg4*reg76;
    reg176=reg106*reg176; reg122=reg25*reg122; reg32=reg182-reg32; reg25=reg84*reg113; reg151=reg64*reg89;
    reg115=reg86*reg115; reg72=reg69*reg72; reg169=0.25*vectors[0][indices[6]+1]; reg194=reg194-reg202; reg182=reg55*reg28;
    reg111=reg16*reg111; reg13=reg201-reg13; reg158=reg158-reg190; reg70=reg200+reg70; reg199=reg24*reg139;
    reg200=reg148*reg179; reg123=reg123-reg166; reg212=reg211+reg212; reg102=reg172-reg102; reg172=reg135*reg124;
    reg210=reg138+reg210; reg198=reg198-reg180; reg201=reg14*reg206; reg203=reg41*reg73; reg209=reg50-reg209;
    reg50=reg48*reg6; reg204=reg120*reg152; reg208=reg148*reg170; reg211=reg131*reg100; reg104=reg134-reg104;
    reg134=reg36*reg135; T reg213=reg55*reg11; reg22=reg22/reg3; T reg214=reg43*reg108; reg130=reg130*reg137;
    reg54=reg54*reg145; reg157=reg161-reg157; reg58=reg58/reg3; reg17=reg197-reg17; reg161=reg34*reg28;
    reg79=reg82*reg79; reg60=reg53-reg60; reg53=reg185*reg57; reg118=reg56*reg118; reg160=reg77*reg160;
    reg40=reg40-reg190; reg117=reg143+reg117; reg1=reg10-reg1; reg78=reg78*reg156; reg68=reg95+reg68;
    reg10=reg34*reg11; reg20=reg103-reg20; reg141=reg147-reg141; reg62=reg62/reg3; reg56=reg83*reg5;
    reg142=reg142/reg3; reg165=reg177-reg165; reg77=reg149*reg155; reg121=reg52+reg121; reg52=reg148*reg189;
    reg184=reg67-reg184; reg195=reg171+reg195; reg109=reg180+reg109; reg67=reg114*reg205; reg146=reg0+reg146;
    reg0=0.25*vectors[0][indices[6]+0]; reg85=reg85/reg51; reg15=reg2-reg15; reg133=reg202+reg133; reg2=reg82*reg91;
    reg87=reg39+reg87; reg153=reg45-reg153; reg93=reg35+reg93; reg26=reg178+reg26; reg38=reg187-reg38;
    reg65=reg65/reg112; reg35=reg92*reg135; reg120=reg120*reg150; reg183=reg183/reg112; reg39=reg188*reg90;
    reg144=reg144*reg69; reg59=reg126*reg59; reg168=reg168*reg159; reg81=reg81-reg166; reg127=reg106*reg127;
    reg75=reg116*reg75; reg122=reg164-reg122; reg45=reg156*reg196; reg98=reg71-reg98; reg191=reg193-reg191;
    reg71=0.25*vectors[0][indices[5]+2]; reg207=reg207/reg112; reg95=reg30*reg136; reg31=reg94-reg31; reg63=reg125+reg63;
    reg21=reg33-reg21; reg8=reg8-reg166; reg176=reg174-reg176; reg175=reg19+reg175; reg23=reg49*reg23;
    reg119=reg119/reg112; reg19=reg148*reg163; reg12=reg29+reg12; reg13=reg13-reg166; reg211=reg104+reg211;
    reg33=reg18*reg207; reg94=0.25*vectors[0][indices[7]+0]; reg201=reg151-reg201; reg103=reg136*reg65; reg104=reg6*reg183;
    reg172=reg17+reg172; reg39=reg53-reg39; reg17=reg4*reg139; reg214=reg157-reg214; reg53=reg47*reg113;
    reg21=reg45+reg21; reg45=reg110*reg88; reg106=reg154*reg135; reg125=reg55*reg8; reg126=reg34*reg123;
    reg121=reg26+reg121; reg25=reg199-reg25; reg26=reg9*reg152; reg165=reg111+reg165; reg87=reg87-reg166;
    reg109=reg109+reg0; reg93=reg93-reg166; reg56=reg77-reg56; reg77=reg148*reg81; reg158=reg158-reg71;
    reg98=reg98/reg112; reg78=reg141-reg78; reg111=reg119*reg15; reg40=reg71+reg40; reg117=0.5*reg117;
    reg141=reg34*reg8; reg143=reg55*reg123; reg147=reg16*reg101; reg203=reg209+reg203; reg210=reg52+reg210;
    reg35=reg120-reg35; reg120=reg148*reg28; reg212=reg52+reg212; reg153=reg180+reg153; reg194=reg194+reg169;
    reg200=reg138+reg200; reg52=reg34*reg189; reg168=reg122-reg168; reg100=reg42*reg100; reg122=reg96*reg188;
    reg138=reg113*reg167; reg115=reg72-reg115; reg72=reg149*reg162; reg151=reg64*reg61; reg36=reg150*reg36;
    reg157=reg34*reg170; reg19=reg29+reg19; reg184=reg202+reg184; reg181=reg69*reg181; reg95=reg67-reg95;
    reg92=reg92*reg152; reg44=reg86*reg44; reg29=reg159*reg85; reg145=reg74*reg145; reg133=reg169+reg133;
    reg2=reg38-reg2; reg137=reg116*reg137; reg63=reg208+reg63; reg12=reg208+reg12; reg195=reg191+reg195;
    reg198=reg0+reg198; reg38=0.25*vectors[0][indices[6]+2]; reg155=reg114*reg155; reg67=reg83*reg22; reg68=reg175+reg68;
    reg5=reg76*reg5; reg69=reg206*reg58; reg176=reg79+reg176; reg161=reg213+reg161; reg60=reg60/reg3;
    reg50=reg102+reg50; reg160=reg118-reg160; reg23=reg31+reg23; reg9=reg9*reg150; reg134=reg204-reg134;
    reg129=reg32/reg129; reg146=0.5*reg146; reg70=reg190-reg70; reg107=reg107*reg49; reg182=reg10+reg182;
    reg10=reg46*reg108; reg54=reg130-reg54; reg144=reg59-reg144; reg97=reg97/reg3; reg127=reg75-reg127;
    reg31=reg139*reg142; reg32=reg57*reg62; reg48=reg18*reg48; reg20=reg20/reg3; reg59=0.25*vectors[0][indices[7]+1];
    reg1=reg1/reg3; reg135=reg80*reg135; reg153=reg0+reg153; reg158=reg38+reg158; reg109=reg109-reg94;
    reg133=reg133-reg59; reg35=reg35/reg172; reg134=reg134/reg172; reg194=reg59+reg194; reg135=reg9-reg135;
    reg0=0.25*vectors[0][indices[7]+2]; reg70=reg71+reg70; reg69=reg67-reg69; reg9=reg131*reg60; reg160=reg160/reg3;
    reg27=reg27/reg172; reg40=reg38+reg40; reg67=reg129*reg117; reg71=reg15*reg1; reg74=reg76*reg20;
    reg31=reg32-reg31; reg32=reg108*reg97; reg106=reg26-reg106; reg52=reg200+reg52; reg212=reg105*reg212;
    reg26=reg64*reg90; reg75=reg206*reg185; reg210=reg179*reg210; reg192=reg192/reg172; reg198=reg94+reg198;
    reg79=reg14*reg57; reg86=reg188*reg89; reg154=reg150*reg154; reg152=reg80*reg152; reg104=reg39+reg104;
    reg184=reg169+reg184; reg36=reg92-reg36; reg2=reg2-reg166; reg195=0.5*reg195; reg39=reg82*reg73;
    reg151=reg72-reg151; reg100=reg115+reg100; reg10=reg54+reg10; reg110=reg110*reg196; reg54=reg41*reg91;
    reg138=reg122-reg138; reg17=reg53-reg17; reg157=reg19+reg157; reg61=reg136*reg61; reg211=reg214+reg211;
    reg95=reg29+reg95; reg162=reg114*reg162; reg19=reg30*reg188; reg113=reg113*reg205; reg107=reg144-reg107;
    reg127=reg48+reg127; reg84=reg84*reg114; reg24=reg24*reg76; reg23=reg50+reg23; reg4=reg83*reg4;
    reg78=reg165+reg78; reg47=reg149*reg47; reg168=reg176+reg168; reg43=reg43*reg42; reg181=reg44-reg181;
    reg101=reg156*reg101; reg68=0.5*reg68; reg145=reg137-reg145; reg5=reg155-reg5; reg46=reg131*reg46;
    reg29=reg129*reg146; reg103=reg111-reg103; reg44=reg49*reg98; reg45=reg25-reg45; reg203=reg203-reg166;
    reg12=reg163*reg12; reg25=reg148*reg87; reg48=reg55*reg13; reg141=reg143+reg141; reg50=reg34*reg93;
    reg147=reg56+reg147; reg125=reg126+reg125; reg53=reg34*reg13; reg56=reg55*reg93; reg63=reg186*reg63;
    reg21=reg21-reg166; reg33=reg201-reg33; reg161=reg77+reg161; reg72=reg34*reg81; reg80=reg148*reg8;
    reg182=reg77+reg182; reg121=0.5*reg121; reg120=reg213+reg120; reg77=reg139*reg22; reg101=reg5+reg101;
    reg9=reg69-reg9; reg96=reg96*reg136; reg41=reg41*reg85; reg167=reg114*reg167; reg196=reg16*reg196;
    reg19=reg113-reg19; reg5=reg57*reg58; reg153=reg94+reg153; reg32=reg31+reg32; reg4=reg47-reg4;
    reg23=0.5*reg23; reg16=reg206*reg62; reg31=reg83*reg142; reg47=reg148*reg13; reg161=reg28*reg161;
    reg124=reg124/reg172; reg28=reg192*reg109; reg73=reg159*reg73; reg69=reg148*reg21; reg61=reg162-reg61;
    reg110=reg17+reg110; reg17=reg129*reg195; reg53=reg56+reg53; reg74=reg71-reg74; reg182=reg11*reg182;
    reg11=reg42*reg160; reg71=reg27*reg198; reg205=reg149*reg205; reg30=reg64*reg30; reg72=reg120+reg72;
    reg48=reg50+reg48; reg147=reg45+reg147; reg157=reg170*reg157; reg90=reg136*reg90; reg185=reg15*reg185;
    reg63=reg12+reg63; reg188=reg188*reg65; reg12=reg119*reg57; reg135=reg135/reg172; reg52=reg189*reg52;
    reg45=reg34*reg87; reg80=reg143+reg80; reg50=reg18*reg183; reg26=reg75-reg26; reg75=reg129*reg121;
    reg95=reg95-reg166; reg92=reg6*reg207; reg79=reg86-reg79; reg33=reg33-reg166; reg40=reg40-reg0;
    reg212=reg210+reg212; reg125=reg25+reg125; reg106=reg106/reg172; reg141=reg25+reg141; reg25=reg55*reg2;
    reg86=reg34*reg203; reg158=reg0+reg158; reg104=reg104-reg166; reg103=reg44+reg103; reg44=reg34*reg2;
    reg94=reg55*reg203; reg107=reg127+reg107; reg84=reg24-reg84; reg88=reg156*reg88; reg24=reg134*reg194;
    reg102=reg35*reg133; reg78=0.5*reg78; reg67=2*reg67; reg168=0.5*reg168; reg36=reg36/reg172;
    reg105=reg129*reg68; reg43=reg181-reg43; reg184=reg59+reg184; reg145=reg46+reg145; reg211=0.5*reg211;
    reg54=reg138-reg54; reg38=reg70+reg38; reg100=reg10+reg100; reg39=reg151+reg39; reg29=2*reg29;
    reg154=reg152-reg154; reg10=reg148*reg2; reg142=reg76*reg142; reg157=reg63+reg157; reg47=reg56+reg47;
    reg46=reg34*reg21; reg62=reg15*reg62; reg44=reg94+reg44; reg4=reg196+reg4; reg88=reg84-reg88;
    reg56=reg129*reg168; reg25=reg86+reg25; reg57=reg57*reg1; reg29=reg146*reg29; reg139=reg139*reg20;
    reg147=0.5*reg147; reg59=reg148*reg95; reg105=2*reg105; reg101=reg110+reg101; reg52=reg212+reg52;
    reg39=reg54+reg39; reg141=reg8*reg141; reg8=reg129*reg23; reg125=reg123*reg125; reg107=0.5*reg107;
    reg102=reg24-reg102; reg24=reg135*reg40; reg92=reg79-reg92; reg54=reg129*reg78; reg67=reg117*reg67;
    reg50=reg26+reg50; reg45=reg80+reg45; reg43=reg145+reg43; reg75=2*reg75; reg26=reg36*reg184;
    reg188=reg12-reg188; reg172=reg154/reg172; reg6=reg6*reg98; reg14=reg14*reg15; reg89=reg136*reg89;
    reg100=0.5*reg100; reg65=reg64*reg65; reg119=reg206*reg119; reg90=reg185-reg90; reg0=reg38+reg0;
    reg183=reg49*reg183; reg12=reg129*reg211; reg38=reg131*reg97; reg31=reg16-reg31; reg48=reg69+reg48;
    reg16=reg108*reg60; reg103=reg103-reg166; reg5=reg77-reg5; reg63=reg55*reg104; reg53=reg69+reg53;
    reg182=reg161+reg182; reg32=reg32-reg166; reg41=reg19+reg41; reg19=reg34*reg33; reg74=reg11+reg74;
    reg73=reg61+reg73; reg71=reg28-reg71; reg72=reg81*reg72; reg11=reg106*reg158; reg85=reg82*reg85;
    reg30=reg205-reg30; reg9=reg9-reg166; reg17=2*reg17; reg28=reg34*reg104; reg167=reg96-reg167;
    reg91=reg159*reg91; reg61=reg124*reg153; reg64=reg55*reg33; reg45=reg87*reg45; reg125=reg141+reg125;
    reg11=reg24-reg11; reg48=reg93*reg48; reg46=reg47+reg46; reg53=reg13*reg53; reg13=reg129*reg147;
    reg72=reg182+reg72; reg101=0.5*reg101; reg61=reg71+reg61; elem.epsilon[0][0]=reg61; reg17=reg195*reg17;
    reg8=2*reg8; reg24=reg129*reg107; reg88=reg4+reg88; reg56=2*reg56; reg43=0.5*reg43;
    reg105=reg68*reg105; reg29=reg157+reg29; reg26=reg102-reg26; elem.epsilon[0][1]=reg26; reg4=reg129*reg100;
    reg12=2*reg12; reg47=reg172*reg0; reg58=reg15*reg58; reg22=reg76*reg22; reg20=reg83*reg20;
    reg1=reg206*reg1; reg97=reg42*reg97; reg142=reg62-reg142; reg108=reg108*reg160; reg139=reg57-reg139;
    reg38=reg31+reg38; reg16=reg5-reg16; reg5=reg55*reg9; reg15=reg34*reg32; reg31=reg34*reg9;
    reg57=reg55*reg32; reg74=reg74-reg166; reg73=reg41+reg73; reg30=reg85+reg30; reg91=reg167-reg91;
    reg54=2*reg54; reg67=reg52+reg67; reg75=reg121*reg75; reg41=reg148*reg103; reg19=reg63+reg19;
    reg207=reg49*reg207; reg14=reg89-reg14; reg65=reg119-reg65; reg98=reg18*reg98; reg64=reg28+reg64;
    reg183=reg90+reg183; reg6=reg188+reg6; reg50=reg92+reg50; reg18=reg148*reg33; reg44=reg59+reg44;
    reg39=0.5*reg39; reg25=reg59+reg25; reg28=reg34*reg95; reg10=reg94+reg10; reg11=reg47+reg11;
    elem.epsilon[0][2]=reg11; reg13=2*reg13; reg108=reg139+reg108; reg56=reg168*reg56; reg64=reg41+reg64;
    reg12=reg211*reg12; reg19=reg41+reg19; reg46=reg21*reg46; reg73=0.5*reg73; reg38=reg16+reg38;
    reg48=reg53+reg48; reg16=reg61+reg26; reg65=reg98+reg65; reg28=reg10+reg28; reg4=2*reg4;
    reg10=reg148*reg9; reg21=reg148*reg74; reg207=reg14-reg207; reg5=reg15+reg5; reg105=reg29+reg105;
    reg24=2*reg24; reg31=reg57+reg31; reg18=reg63+reg18; reg14=reg34*reg103; reg88=0.5*reg88;
    reg8=reg23*reg8; reg60=reg42*reg60; reg58=reg22-reg58; reg15=reg129*reg43; reg75=reg67+reg75;
    reg20=reg1-reg20; reg44=reg2*reg44; reg50=0.5*reg50; reg160=reg131*reg160; reg1=reg129*reg39;
    reg17=reg72+reg17; reg45=reg125+reg45; reg91=reg30+reg91; reg54=reg78*reg54; reg2=reg129*reg101;
    reg25=reg203*reg25; reg97=reg142+reg97; reg183=reg6+reg183; reg15=2*reg15; reg6=reg198*reg134;
    reg8=reg17+reg8; reg207=reg65+reg207; reg91=0.5*reg91; reg54=reg75+reg54; reg19=reg33*reg19;
    reg17=reg109*reg35; reg1=2*reg1; reg22=reg192*reg133; reg23=reg129*reg73; reg29=reg129*reg88;
    reg60=reg58-reg60; reg14=reg18+reg14; reg20=reg160+reg20; reg18=reg129*reg50; reg30=reg27*reg194;
    reg2=2*reg2; reg64=reg104*reg64; reg97=reg108+reg97; reg183=0.5*reg183; reg56=reg105+reg56;
    reg24=reg107*reg24; reg13=reg147*reg13; reg12=reg45+reg12; reg46=reg48+reg46; reg25=reg44+reg25;
    reg38=0.5*reg38; reg33=reg34*reg74; reg10=reg57+reg10; reg4=reg100*reg4; reg5=reg21+reg5;
    reg16=reg11+reg16; reg28=reg95*reg28; reg31=reg21+reg31; reg1=reg39*reg1; reg30=reg22-reg30;
    reg21=reg124*reg184; reg15=reg43*reg15; reg22=reg153*reg36; reg17=reg6-reg17; reg28=reg25+reg28;
    reg207=0.5*reg207; reg4=reg12+reg4; reg6=reg129*reg183; reg18=2*reg18; reg14=reg103*reg14;
    reg64=reg19+reg64; reg13=reg46+reg13; reg56=reg128*reg56; reg23=2*reg23; reg97=0.5*reg97;
    reg27=reg27*reg158; reg60=reg20+reg60; reg192=reg192*reg40; reg2=reg101*reg2; reg12=reg129*reg91;
    reg31=reg9*reg31; reg16=reg16/3; reg198=reg198*reg106; reg109=reg109*reg135; reg5=reg32*reg5;
    reg33=reg10+reg33; reg29=2*reg29; reg54=reg66*reg54; reg24=reg8+reg24; reg8=reg129*reg38;
    reg124=reg124*reg0; reg8=2*reg8; reg6=2*reg6; reg60=0.5*reg60; reg9=reg129*reg97;
    reg135=reg133*reg135; reg18=reg50*reg18; reg2=reg13+reg2; reg106=reg194*reg106; reg56=0.125*reg56;
    reg14=reg64+reg14; reg40=reg35*reg40; reg29=reg88*reg29; reg158=reg134*reg158; reg21=reg30+reg21;
    reg24=reg99*reg24; reg22=reg17-reg22; reg15=reg4+reg15; reg54=0.125*reg54; reg198=reg109-reg198;
    reg153=reg153*reg172; reg12=2*reg12; reg23=reg73*reg23; reg27=reg192-reg27; reg33=reg74*reg33;
    reg5=reg31+reg5; reg4=reg61-reg16; reg1=reg28+reg1; reg10=reg129*reg207; reg13=reg26-reg16;
    reg4=pow(reg4,2); reg153=reg198+reg153; reg21=reg22+reg21; reg0=reg36*reg0; reg124=reg27+reg124;
    reg40=reg158-reg40; reg16=reg11-reg16; reg172=reg184*reg172; reg13=pow(reg13,2); reg106=reg135-reg106;
    reg17=reg129*reg60; reg56=reg54+reg56; reg29=reg2+reg29; reg9=2*reg9; reg8=reg38*reg8;
    reg33=reg5+reg33; reg23=reg1+reg23; reg12=reg91*reg12; reg24=0.125*reg24; reg15=reg173*reg15;
    reg10=2*reg10; reg6=reg183*reg6; reg18=reg14+reg18; reg0=reg40-reg0; reg10=reg207*reg10;
    reg12=reg23+reg12; reg8=reg33+reg8; reg124=reg153+reg124; reg1=0.5*reg21; elem.epsilon[0][3]=reg1;
    reg9=reg97*reg9; reg24=reg56+reg24; reg15=0.125*reg15; reg16=pow(reg16,2); reg13=reg4+reg13;
    reg29=reg37*reg29; reg17=2*reg17; reg106=reg172+reg106; reg6=reg18+reg6; reg10=reg6+reg10;
    reg15=reg24+reg15; reg17=reg60*reg17; reg2=0.5*reg124; elem.epsilon[0][4]=reg2; reg29=0.125*reg29;
    reg12=reg51*reg12; reg9=reg8+reg9; reg0=reg106+reg0; reg16=reg13+reg16; reg21=reg21*reg1;
    reg29=reg15+reg29; reg12=0.125*reg12; reg4=0.5*reg0; elem.epsilon[0][5]=reg4; reg124=reg124*reg2;
    reg21=reg16+reg21; reg10=reg112*reg10; reg17=reg9+reg17; reg26=reg26-reg166; reg61=reg61-reg166;
    reg0=reg0*reg4; reg12=reg29+reg12; reg124=reg21+reg124; reg10=0.125*reg10; reg17=reg3*reg17;
    reg17=0.125*reg17; reg3=reg55*reg61; reg5=reg55*reg26; reg61=reg34*reg61; reg11=reg11-reg166;
    reg6=reg34*reg26; reg10=reg12+reg10; reg0=reg124+reg0; reg26=reg148*reg26; reg17=reg10+reg17;
    reg6=reg3+reg6; reg0=1.5*reg0; reg26=reg3+reg26; reg3=reg34*reg11; reg11=reg148*reg11;
    reg5=reg61+reg5; elem.sigma_von_mises=pow(reg0,0.5); elem.ener=reg17/2; elem.sigma[0][5]=reg129*reg4; elem.sigma[0][4]=reg129*reg2;
    elem.sigma[0][3]=reg129*reg1; elem.sigma[0][2]=reg26+reg3; elem.sigma[0][1]=reg11+reg6; elem.sigma[0][0]=reg5+reg11;
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1-var_inter[1]; T reg1=1-var_inter[0]; T reg2=1-var_inter[2]; T reg3=reg2*reg1; T reg4=reg2*var_inter[0];
    T reg5=reg0*reg2; T reg6=reg0*reg1; T reg7=reg0*var_inter[0]; T reg8=elem.pos(0)[2]*reg3; T reg9=elem.pos(1)[2]*reg4;
    T reg10=reg2*var_inter[1]; T reg11=var_inter[0]*var_inter[1]; T reg12=elem.pos(0)[1]*reg3; T reg13=elem.pos(1)[1]*reg4; T reg14=elem.pos(1)[1]*reg5;
    T reg15=elem.pos(0)[1]*reg5; T reg16=elem.pos(1)[1]*reg7; T reg17=elem.pos(1)[2]*reg5; T reg18=elem.pos(0)[2]*reg5; T reg19=elem.pos(0)[1]*reg6;
    T reg20=reg7*elem.pos(1)[2]; T reg21=elem.pos(0)[2]*reg6; T reg22=var_inter[1]*reg1; T reg23=elem.pos(2)[2]*reg4; T reg24=reg16+reg19;
    T reg25=reg9+reg8; T reg26=elem.pos(2)[1]*reg4; T reg27=reg13+reg12; reg14=reg14-reg15; T reg28=elem.pos(2)[1]*reg10;
    reg17=reg17-reg18; T reg29=elem.pos(2)[2]*reg10; T reg30=reg11*elem.pos(2)[2]; T reg31=reg20+reg21; T reg32=reg11*elem.pos(2)[1];
    T reg33=elem.pos(3)[1]*reg10; T reg34=reg31+reg30; T reg35=elem.pos(1)[0]*reg5; T reg36=elem.pos(0)[0]*reg5; reg28=reg14+reg28;
    reg14=reg0*var_inter[2]; T reg37=var_inter[2]*reg1; T reg38=elem.pos(3)[1]*reg3; reg26=reg26-reg27; T reg39=elem.pos(3)[2]*reg3;
    T reg40=reg22*elem.pos(3)[1]; T reg41=reg24+reg32; reg23=reg23-reg25; reg29=reg17+reg29; reg17=elem.pos(0)[0]*reg3;
    T reg42=elem.pos(1)[0]*reg4; T reg43=elem.pos(3)[2]*reg10; T reg44=reg22*elem.pos(3)[2]; T reg45=elem.pos(4)[2]*reg14; reg29=reg29-reg43;
    T reg46=var_inter[2]*var_inter[0]; T reg47=reg42+reg17; T reg48=elem.pos(2)[0]*reg4; T reg49=elem.pos(4)[1]*reg14; reg28=reg28-reg33;
    T reg50=elem.pos(0)[0]*reg6; T reg51=elem.pos(1)[0]*reg7; T reg52=elem.pos(4)[2]*reg37; reg23=reg39+reg23; reg39=elem.pos(4)[2]*reg6;
    T reg53=reg44+reg34; reg35=reg35-reg36; T reg54=elem.pos(4)[1]*reg6; T reg55=elem.pos(2)[0]*reg10; T reg56=reg41+reg40;
    reg38=reg26+reg38; reg26=elem.pos(4)[1]*reg37; T reg57=reg7*elem.pos(5)[2]; T reg58=reg11*elem.pos(2)[0]; T reg59=reg51+reg50;
    T reg60=reg7*elem.pos(5)[1]; reg38=reg38-reg26; T reg61=elem.pos(3)[0]*reg3; reg39=reg39-reg53; T reg62=elem.pos(5)[2]*reg46;
    reg23=reg23-reg52; reg54=reg54-reg56; reg48=reg48-reg47; T reg63=elem.pos(5)[1]*reg46; T reg64=elem.pos(3)[0]*reg10;
    reg55=reg35+reg55; reg35=var_inter[2]*var_inter[1]; reg28=reg28-reg49; T reg65=elem.pos(5)[1]*reg14; T reg66=elem.pos(5)[2]*reg14;
    reg29=reg29-reg45; T reg67=reg11*elem.pos(6)[2]; reg57=reg39+reg57; reg39=elem.pos(6)[1]*reg46; reg38=reg38-reg63;
    reg55=reg55-reg64; T reg68=elem.pos(4)[0]*reg14; reg23=reg23-reg62; T reg69=elem.pos(6)[2]*reg46; T reg70=reg22*elem.pos(3)[0];
    T reg71=reg11*elem.pos(6)[1]; reg60=reg54+reg60; reg66=reg29+reg66; reg29=elem.pos(6)[1]*reg35; reg65=reg28+reg65;
    reg48=reg61+reg48; reg28=elem.pos(4)[0]*reg37; reg54=elem.pos(6)[2]*reg35; reg61=reg59+reg58; reg69=reg23+reg69;
    reg71=reg60+reg71; reg23=elem.pos(7)[2]*reg37; reg60=reg22*elem.pos(7)[1]; T reg72=elem.pos(4)[0]*reg6; T reg73=reg70+reg61;
    reg67=reg57+reg67; reg57=elem.pos(7)[1]*reg35; reg29=reg65+reg29; reg48=reg48-reg28; reg65=elem.pos(5)[0]*reg46;
    T reg74=reg22*elem.pos(7)[2]; T reg75=elem.pos(5)[0]*reg14; reg55=reg55-reg68; T reg76=elem.pos(7)[2]*reg35; reg54=reg66+reg54;
    reg39=reg38+reg39; reg38=elem.pos(7)[1]*reg37; reg38=reg39+reg38; reg29=reg29-reg57; reg39=1+(*f.m).poisson_ratio;
    reg66=reg7*elem.pos(5)[0]; reg72=reg72-reg73; reg48=reg48-reg65; T reg77=elem.pos(6)[0]*reg46; reg74=reg67+reg74;
    reg23=reg69+reg23; reg67=elem.pos(6)[0]*reg35; reg54=reg54-reg76; reg75=reg55+reg75; reg60=reg71+reg60;
    reg55=reg29*reg74; reg69=reg38*reg74; reg71=reg23*reg60; T reg78=reg54*reg60; reg39=reg39/(*f.m).elastic_modulus;
    reg67=reg75+reg67; reg75=reg11*elem.pos(6)[0]; T reg79=elem.pos(7)[0]*reg37; reg77=reg48+reg77; reg66=reg72+reg66;
    reg48=elem.pos(7)[0]*reg35; reg72=pow(reg39,2); reg79=reg77+reg79; reg71=reg69-reg71; reg67=reg67-reg48;
    reg78=reg55-reg78; reg55=reg29*reg23; reg69=reg54*reg38; reg75=reg66+reg75; reg66=reg22*elem.pos(7)[0];
    reg77=reg79*reg78; T reg80=reg67*reg71; T reg81=1.0/(*f.m).elastic_modulus; T reg82=(*f.m).poisson_ratio/(*f.m).elastic_modulus; reg69=reg55-reg69;
    reg66=reg75+reg66; reg39=reg39*reg72; reg77=reg80-reg77; reg55=reg54*reg66; reg75=reg66*reg69;
    reg80=reg82*reg39; reg39=reg81*reg39; T reg83=reg82*reg72; T reg84=reg79*reg74; T reg85=reg23*reg66;
    reg72=reg81*reg72; reg74=reg67*reg74; reg54=reg54*reg79; T reg86=reg79*reg60; reg85=reg84-reg85;
    reg23=reg67*reg23; reg84=reg38*reg66; reg60=reg67*reg60; reg75=reg77+reg75; reg77=reg82*reg39;
    T reg87=reg82*reg80; reg39=reg81*reg39; T reg88=reg81*reg72; T reg89=reg82*reg83; reg72=reg82*reg72;
    reg55=reg74-reg55; reg66=reg29*reg66; reg38=reg67*reg38; reg66=reg60-reg66; reg54=reg23-reg54;
    reg79=reg29*reg79; reg39=reg39-reg87; reg77=reg87+reg77; reg80=reg81*reg80; reg83=reg81*reg83;
    reg88=reg88-reg89; reg72=reg89+reg72; reg55=reg55/reg75; reg78=reg78/reg75; reg84=reg86-reg84;
    reg85=reg85/reg75; reg71=reg71/reg75; reg23=reg81*reg39; reg29=reg14*reg71; reg80=reg87+reg80;
    reg88=reg81*reg88; reg60=reg46*reg78; reg72=reg82*reg72; reg67=reg89+reg83; reg74=reg3*reg55;
    reg81=reg46*reg55; reg86=reg10*reg71; reg87=reg3*reg78; T reg90=reg10*reg85; reg66=reg66/reg75;
    reg69=reg69/reg75; reg54=reg54/reg75; reg84=reg84/reg75; reg79=reg38-reg79; reg38=reg82*reg77;
    T reg91=reg14*reg85; T reg92=reg37*reg55; T reg93=reg5*reg71; reg79=reg79/reg75; T reg94=reg14*reg84;
    T reg95=reg5*reg85; reg67=reg82*reg67; T reg96=reg10*reg84; reg72=reg88-reg72; reg88=reg37*reg78;
    T reg97=reg4*reg66; T reg98=reg37*reg66; T reg99=reg35*reg84; T reg100=reg7*reg69; T reg101=reg3*reg66;
    T reg102=reg5*reg84; T reg103=reg4*reg78; T reg104=reg60+reg29; T reg105=reg35*reg71; T reg106=reg35*reg85;
    T reg107=reg46*reg66; T reg108=reg7*reg54; T reg109=reg22*reg69; T reg110=reg87+reg86; T reg111=reg74+reg90;
    T reg112=reg22*reg54; reg38=reg23-reg38; reg82=reg82*reg80; reg23=reg4*reg55; T reg113=reg81+reg91;
    T reg114=reg11*reg79; T reg115=reg11*reg69; T reg116=reg86-reg103; T reg117=reg105-reg60; T reg118=reg97+reg102;
    T reg119=reg99-reg107; T reg120=reg87-reg93; T reg121=reg23-reg90; T reg122=reg101+reg96; T reg123=reg6*reg69;
    T reg124=reg22*reg79; T reg125=reg101-reg102; T reg126=reg6*reg79; T reg127=reg110+reg109; reg111=reg111+reg112;
    T reg128=reg7*reg79; T reg129=reg95+reg23; T reg130=reg107+reg94; reg82=reg38-reg82; reg38=reg91-reg92;
    reg104=reg104+reg100; T reg131=reg93+reg103; T reg132=reg92+reg106; T reg133=reg88-reg29; reg67=reg72-reg67;
    reg72=reg96-reg97; T reg134=reg11*reg54; T reg135=reg81-reg106; T reg136=reg6*reg54; T reg137=reg95-reg74;
    T reg138=reg98-reg94; T reg139=reg113+reg108; T reg140=reg98+reg99; T reg141=reg88+reg105; T reg142=reg124+reg122;
    reg133=reg133+reg123; reg38=reg38-reg136; reg116=reg116-reg115; reg138=reg138+reg126; T reg143=0.5*reg139;
    reg137=reg137+reg136; reg72=reg72-reg114; reg121=reg121+reg134; reg130=reg128+reg130; T reg144=reg124-reg140;
    T reg145=reg109-reg141; reg135=reg135-reg134; reg132=reg132-reg112; reg117=reg117+reg115; T reg146=(*f.m).deltaT*(*f.m).alpha;
    reg119=reg119+reg114; T reg147=0.5*reg104; reg67=reg67/reg82; T reg148=reg108-reg129; reg80=reg80/reg82;
    reg39=reg39/reg82; reg82=reg77/reg82; reg131=reg131-reg100; reg118=reg118-reg128; reg125=reg125-reg126;
    reg120=reg120-reg123; reg77=0.5*reg127; T reg149=0.5*reg111; T reg150=0.5*reg145; T reg151=0.5*reg117;
    T reg152=0.5*reg133; T reg153=0.5*reg144; T reg154=0.5*reg138; T reg155=0.5*reg132; T reg156=0.5*reg72;
    T reg157=0.5*reg116; T reg158=0.5*reg148; T reg159=0.5*reg130; T reg160=0.5*reg118; T reg161=0.5*reg38;
    T reg162=0.5*reg142; T reg163=0.5*reg125; T reg164=0.5*reg131; T reg165=reg80*reg146; T reg166=reg82*reg146;
    T reg167=0.5*reg137; T reg168=0.5*reg119; T reg169=reg39*reg146; T reg170=reg67*reg143; T reg171=reg67*reg147;
    T reg172=0.5*reg135; T reg173=0.5*reg121; T reg174=reg67*reg149; T reg175=reg67*reg77; T reg176=0.5*reg120;
    T reg177=reg67*reg154; T reg178=reg67*reg153; T reg179=reg67*reg155; T reg180=reg67*reg158; T reg181=reg67*reg151;
    T reg182=reg67*reg162; T reg183=reg67*reg172; T reg184=reg67*reg173; T reg185=reg67*reg161; T reg186=reg67*reg168;
    T reg187=reg67*reg152; T reg188=reg67*reg150; T reg189=reg67*reg157; T reg190=reg39*reg130; T reg191=reg39*reg111;
    T reg192=reg39*reg142; T reg193=2*reg175; T reg194=reg166+reg165; T reg195=reg67*reg167; T reg196=reg39*reg127;
    T reg197=reg39*reg139; T reg198=reg169+reg166; T reg199=reg67*reg176; T reg200=reg67*reg163; T reg201=reg67*reg156;
    reg174=2*reg174; T reg202=reg160*reg67; reg171=2*reg171; T reg203=reg67*reg159; T reg204=reg67*reg164;
    T reg205=reg39*reg104; T reg206=2*reg170; T reg207=reg118*reg39; T reg208=reg82*reg139; T reg209=reg82*reg111;
    reg187=2*reg187; T reg210=reg39*reg125; T reg211=reg39*reg133; reg185=2*reg185; T reg212=reg82*reg127;
    reg177=2*reg177; T reg213=reg77*reg171; T reg214=reg111*reg197; T reg215=reg80*reg142; T reg216=reg39*reg132;
    T reg217=reg82*reg117; reg188=2*reg188; T reg218=reg39*reg135; T reg219=reg39*reg145; reg179=2*reg179;
    T reg220=reg82*reg133; T reg221=reg139*reg191; reg178=2*reg178; T reg222=reg82*reg145; reg180=2*reg180;
    T reg223=reg131*reg39; reg204=2*reg204; T reg224=reg39*reg144; T reg225=reg118*reg80; reg202=2*reg202;
    T reg226=reg142*reg190; T reg227=reg39*reg119; reg200=2*reg200; T reg228=reg80*reg144; T reg229=reg82*reg120;
    T reg230=reg149*reg206; T reg231=reg127*reg205; T reg232=reg39*reg137; reg181=2*reg181; T reg233=reg80*reg139;
    T reg234=reg39*reg38; T reg235=reg39*reg117; reg183=2*reg183; T reg236=reg39*reg138; T reg237=reg82*reg104;
    reg186=2*reg186; T reg238=reg80*reg111; T reg239=reg130*reg192; T reg240=reg80*reg130; T reg241=reg39*reg72;
    T reg242=reg80*reg138; reg203=2*reg203; T reg243=2*reg182; T reg244=reg80*reg72; reg189=2*reg189;
    T reg245=reg39*reg116; reg184=2*reg184; reg201=2*reg201; T reg246=reg80*reg125; T reg247=reg143*reg174;
    T reg248=reg104*reg196; reg199=2*reg199; T reg249=reg198+reg165; T reg250=reg39*reg120; T reg251=reg169+reg194;
    T reg252=reg39*reg121; T reg253=reg82*reg116; T reg254=reg39*reg148; T reg255=reg131*reg82; reg195=2*reg195;
    T reg256=reg7*var_inter[2]; T reg257=reg22*reg2; T reg258=reg147*reg193; T reg259=reg80*reg119; T reg260=reg142*reg207;
    T reg261=reg250*reg104; T reg262=reg149*reg183; T reg263=reg202*reg77; T reg264=reg179*reg172; T reg265=reg219*reg117;
    T reg266=reg255*reg142; T reg267=reg176*reg243; T reg268=reg125*reg212; T reg269=reg142*reg210; T reg270=reg125*reg241;
    T reg271=reg127*reg235; T reg272=reg38*reg234; T reg273=reg200*reg77; T reg274=reg199*reg151; T reg275=reg135*reg232;
    T reg276=reg173*reg174; T reg277=reg138*reg224; T reg278=reg138*reg227; T reg279=reg80*reg121; T reg280=reg152*reg171;
    T reg281=reg38*reg197; T reg282=reg116*reg196; T reg283=reg138*reg190; T reg284=reg142*reg229; T reg285=reg151*reg204;
    T reg286=reg254*reg135; T reg287=reg188*reg77; T reg288=reg111*reg216; T reg289=reg152*reg181; T reg290=reg38*reg191;
    T reg291=reg142*reg220; T reg292=reg118*reg236; T reg293=reg125*reg190; T reg294=reg143*reg180; T reg295=reg172*reg174;
    T reg296=reg117*reg196; T reg297=reg142*reg192; T reg298=reg104*reg223; T reg299=reg142*reg238; T reg300=reg215*reg117;
    T reg301=reg168*reg193; T reg302=reg149*reg243; T reg303=reg185*reg172; T reg304=reg211*reg117; T reg305=reg118*reg210;
    T reg306=reg125*reg236; T reg307=reg152*reg187; T reg308=reg142*reg241; T reg309=reg172*reg206; T reg310=reg80*reg38;
    T reg311=reg117*reg205; T reg312=reg184*reg173; T reg313=reg245*reg116; T reg314=reg201*reg77; T reg315=reg253*reg142;
    T reg316=reg172*reg183; T reg317=reg195*reg143; T reg318=reg117*reg235; T reg319=reg125*reg192; T reg320=reg131*reg205;
    T reg321=reg151*reg181; T reg322=reg135*reg218; T reg323=reg152*reg243; T reg324=reg137*reg216; T reg325=reg138*reg212;
    T reg326=reg176*reg188; T reg327=reg187*reg77; T reg328=reg111*reg234; T reg329=reg151*reg188; T reg330=reg135*reg216;
    T reg331=reg173*reg206; T reg332=reg116*reg205; T reg333=reg162*reg174; T reg334=reg252*reg139; T reg335=reg111*reg215;
    T reg336=reg119*reg207; T reg337=reg164*reg171; T reg338=reg197*reg148; T reg339=reg138*reg241; T reg340=reg133*reg205;
    T reg341=reg138*reg207; T reg342=reg148*reg218; T reg343=reg164*reg181; T reg344=reg161*reg206; T reg345=reg77*reg193;
    T reg346=reg137*reg218; T reg347=reg138*reg210; T reg348=reg176*reg181; T reg349=reg111*reg191; T reg350=reg119*reg241;
    T reg351=reg148*reg216; T reg352=reg125*reg207; T reg353=reg206*reg158; T reg354=reg164*reg187; T reg355=reg234*reg148;
    T reg356=reg164*reg188; T reg357=reg189*reg151; T reg358=reg138*reg236; T reg359=reg252*reg135; T reg360=reg80*reg148;
    T reg361=reg156*reg193; T reg362=reg215*reg116; T reg363=reg151*reg193; T reg364=reg135*reg191; T reg365=reg38*reg218;
    T reg366=reg181*reg77; T reg367=reg189*reg147; T reg368=reg111*reg218; T reg369=reg185*reg173; T reg370=reg125*reg210;
    T reg371=reg138*reg192; T reg372=reg211*reg116; T reg373=reg151*reg187; T reg374=reg135*reg234; T reg375=reg80*reg137;
    T reg376=reg214+reg213; T reg377=reg151*reg171; T reg378=reg152*reg188; T reg379=reg135*reg197; T reg380=reg38*reg216;
    reg221=reg258+reg221; T reg381=reg133*reg215; T reg382=reg143*reg179; T reg383=reg104*reg219; T reg384=reg118*reg241;
    T reg385=reg143*reg183; T reg386=reg104*reg235; T reg387=reg130*reg210; T reg388=reg130*reg207; T reg389=reg133*reg223;
    T reg390=reg184*reg158; T reg391=reg131*reg245; T reg392=reg161*reg180; T reg393=reg118*reg227; T reg394=reg130*reg241;
    T reg395=reg139*reg197; T reg396=reg130*reg212; T reg397=reg243*reg147; T reg398=reg189*reg162; T reg399=reg258+reg239;
    T reg400=reg147*reg171; T reg401=reg199*reg152; T reg402=reg250*reg133; T reg403=reg139*reg237; T reg404=reg143*reg171;
    T reg405=reg195*reg161; T reg406=reg104*reg208; T reg407=reg180*reg158; T reg408=reg131*reg223; T reg409=reg130*reg236;
    T reg410=reg142*reg224; T reg411=reg130*reg233; T reg412=reg185*reg158; T reg413=reg147*reg204; T reg414=reg131*reg211;
    T reg415=reg139*reg240; T reg416=reg154*reg193; T reg417=reg160*reg193; T reg418=reg131*reg215; T reg419=reg159*reg206;
    T reg420=reg118*reg192; T reg421=reg147*reg181; T reg422=reg161*reg183; T reg423=reg133*reg235; T reg424=reg139*reg218;
    T reg425=reg139*reg232; T reg426=reg133*reg196; T reg427=reg161*reg174; T reg428=reg188*reg147; T reg429=reg139*reg216;
    T reg430=reg174*reg158; T reg431=reg254*reg139; T reg432=reg199*reg147; T reg433=reg131*reg196; T reg434=reg118*reg190;
    T reg435=reg133*reg245; T reg436=reg161*reg184; T reg437=reg131*reg219; T reg438=reg179*reg158; T reg439=reg161*reg179;
    T reg440=reg164*reg243; T reg441=reg118*reg212; T reg442=reg133*reg219; T reg443=reg118*reg207; T reg444=reg164*reg204;
    T reg445=reg147*reg206; T reg446=reg254*reg148; T reg447=reg159*reg243; T reg448=reg248+reg247; T reg449=reg161*reg185;
    T reg450=reg133*reg211; T reg451=reg152*reg189; T reg452=reg172*reg180; T reg453=reg117*reg223; T reg454=reg183*reg158;
    T reg455=reg125*reg227; T reg456=reg131*reg235; T reg457=reg77*reg203; T reg458=reg252*reg38; T reg459=reg139*reg234;
    T reg460=reg143*reg184; T reg461=reg104*reg245; T reg462=reg142*reg237; T reg463=reg164*reg189; T reg464=reg147*reg187;
    T reg465=reg80*reg135; T reg466=reg142*reg236; T reg467=reg184*reg172; T reg468=reg245*reg117; T reg469=reg180*reg173;
    T reg470=reg252*reg148; T reg471=reg223*reg116; T reg472=reg177*reg77; T reg473=reg152*reg193; T reg474=reg38*reg232;
    T reg475=reg149*reg174; T reg476=reg143*reg203; T reg477=reg130*reg190; T reg478=reg178*reg77; T reg479=reg118*reg224;
    T reg480=reg143*reg206; T reg481=reg142*reg222; T reg482=reg199*reg164; T reg483=reg232*reg148; T reg484=reg142*reg227;
    T reg485=reg104*reg205; T reg486=reg195*reg158; T reg487=reg131*reg250; T reg488=reg130*reg227; T reg489=reg195*reg173;
    T reg490=reg152*reg204; T reg491=reg254*reg38; T reg492=reg125*reg224; T reg493=reg143*reg185; T reg494=reg250*reg116;
    T reg495=reg130*reg224; T reg496=reg186*reg77; T reg497=reg80*reg132; T reg498=reg195*reg172; T reg499=reg104*reg211;
    T reg500=reg250*reg117; T reg501=reg142*reg217; reg226=reg213+reg226; reg213=reg159*reg193; T reg502=reg104*reg215;
    T reg503=reg132*reg252; T reg504=reg82*reg135; T reg505=reg127*reg196; T reg506=reg150*reg193; T reg507=reg120*reg235;
    T reg508=reg167*reg183; T reg509=reg132*reg191; T reg510=reg119*reg210; T reg511=reg127*reg244; T reg512=reg150*reg187;
    T reg513=reg132*reg234; T reg514=reg257*(*f.m).f_vol[2]; T reg515=reg157*reg193; T reg516=reg121*reg191; T reg517=reg150*reg171;
    T reg518=reg82*reg38; T reg519=reg132*reg197; T reg520=reg127*reg245; T reg521=reg149*reg184; T reg522=reg120*reg211;
    T reg523=reg167*reg185; T reg524=reg204*reg162; T reg525=reg150*reg181; T reg526=reg132*reg218; T reg527=reg127*reg225;
    T reg528=reg120*reg215; T reg529=reg256*(*f.m).f_vol[1]; T reg530=reg150*reg188; T reg531=reg176*reg187; T reg532=reg162*reg171;
    T reg533=reg127*reg240; T reg534=reg155*reg183; T reg535=reg145*reg235; T reg536=reg162*reg203; reg231=reg230+reg231;
    T reg537=reg137*reg232; T reg538=reg199*reg176; T reg539=reg157*reg204; T reg540=reg155*reg179; T reg541=reg145*reg219;
    T reg542=reg254*reg121; T reg543=reg187*reg162; T reg544=reg127*reg242; T reg545=reg257*(*f.m).f_vol[0]; T reg546=reg127*reg211;
    T reg547=reg150*reg199; T reg548=reg132*reg232; T reg549=reg149*reg185; T reg550=reg127*reg209; T reg551=reg150*reg204;
    T reg552=reg132*reg254; T reg553=reg82*reg148; T reg554=reg149*reg193; T reg555=reg157*reg189; T reg556=reg252*reg121;
    T reg557=reg167*reg180; T reg558=reg150*reg189; T reg559=reg144*reg227; T reg560=reg157*reg181; reg218=reg121*reg218;
    T reg561=reg144*reg224; T reg562=reg82*reg137; T reg563=reg250*reg120; T reg564=reg195*reg167; T reg565=reg72*reg192;
    T reg566=reg252*reg137; T reg567=reg176*reg189; T reg568=reg157*reg243; T reg569=reg72*reg212; T reg570=reg72*reg241;
    T reg571=reg127*reg249; T reg572=reg157*reg188; T reg573=reg142*reg251; T reg574=reg121*reg216; T reg575=reg72*reg207;
    T reg576=reg254*reg137; T reg577=reg176*reg204; T reg578=reg139*reg249; T reg579=reg72*reg210; T reg580=reg120*reg223;
    T reg581=reg2*reg11; T reg582=reg6*reg2; T reg583=reg7*reg2; T reg584=reg22*var_inter[2]; T reg585=var_inter[2]*reg11;
    T reg586=var_inter[2]*reg6; reg216=reg132*reg216; T reg587=reg127*reg223; T reg588=reg82*reg132; T reg589=reg149*reg180;
    T reg590=reg199*reg162; T reg591=reg157*reg187; T reg592=reg121*reg234; T reg593=reg120*reg219; T reg594=reg167*reg179;
    T reg595=reg127*reg246; reg210=reg144*reg210; T reg596=reg127*reg250; reg207=reg144*reg207; T reg597=reg195*reg149;
    T reg598=reg72*reg224; reg241=reg144*reg241; T reg599=reg144*reg212; T reg600=reg150*reg243; T reg601=reg82*reg121;
    T reg602=reg157*reg171; T reg603=reg144*reg192; T reg604=reg121*reg197; T reg605=reg72*reg227; T reg606=reg120*reg245;
    T reg607=reg167*reg184; T reg608=reg144*reg236; T reg609=reg72*reg190; T reg610=reg144*reg190; T reg611=reg72*reg236;
    T reg612=reg111*reg232; T reg613=reg179*reg173; T reg614=reg219*reg116; T reg615=reg155*reg180; reg223=reg145*reg223;
    T reg616=reg120*reg196; T reg617=reg167*reg174; T reg618=reg188*reg162; T reg619=reg127*reg228; T reg620=reg120*reg205;
    T reg621=reg167*reg206; T reg622=reg155*reg184; reg245=reg145*reg245; reg219=reg127*reg219; T reg623=reg149*reg179;
    T reg624=reg181*reg162; T reg625=reg155*reg174; T reg626=reg145*reg196; T reg627=reg127*reg259; T reg628=reg137*reg197;
    T reg629=reg176*reg171; T reg630=reg77*reg174; T reg631=reg191*reg148; T reg632=reg164*reg193; T reg633=reg119*reg212;
    T reg634=reg243*reg151; T reg635=reg111*reg212; reg236=reg119*reg236; T reg636=reg189*reg77; T reg637=reg183*reg173;
    reg235=reg235*reg116; reg252=reg252*reg111; reg190=reg119*reg190; reg191=reg137*reg191; T reg638=reg176*reg193;
    reg227=reg119*reg227; T reg639=reg204*reg77; reg254=reg254*reg111; reg224=reg119*reg224; T reg640=reg163*reg193;
    T reg641=reg155*reg195; reg250=reg145*reg250; T reg642=reg199*reg77; T reg643=reg153*reg193; T reg644=reg155*reg185;
    reg211=reg145*reg211; T reg645=reg199*reg157; T reg646=reg145*reg215; T reg647=reg155*reg206; T reg648=reg119*reg192;
    reg205=reg145*reg205; reg232=reg121*reg232; reg234=reg137*reg234; T reg649=reg152*reg180; T reg650=reg195*reg162;
    T reg651=reg255*reg38; T reg652=reg157*reg203; T reg653=reg72*reg237; reg611=reg591+reg611; reg254=reg254-reg639;
    T reg654=reg154*reg178; reg491=reg490+reg491; T reg655=reg72*reg310; T reg656=reg345+reg297; T reg657=reg111*reg242;
    T reg658=reg154*reg180; reg328=reg328-reg327; T reg659=reg177*reg173; T reg660=reg38*reg225; T reg661=reg157*reg177;
    T reg662=reg72*reg220; T reg663=reg243*reg172; T reg664=reg515+reg565; T reg665=reg152*reg184; T reg666=reg253*reg38;
    T reg667=reg111*reg225; T reg668=reg180*reg162; reg458=reg451+reg458; T reg669=reg154*reg184; T reg670=reg253*reg111;
    T reg671=reg184*reg77; T reg672=reg180*reg77; reg442=reg439+reg442; T reg673=reg161*reg188; T reg674=reg133*reg588;
    T reg675=reg149*reg201; T reg676=reg255*reg111; T reg677=reg536+reg376; T reg678=reg154*reg188; T reg679=reg133*reg228;
    T reg680=reg72*reg465; T reg681=reg111*reg240; T reg682=reg195*reg152; T reg683=reg38*reg229; T reg684=reg119*reg238;
    T reg685=reg186*reg173; T reg686=reg157*reg186; T reg687=reg77*reg206; T reg688=reg72*reg217; reg609=reg602+reg609;
    reg474=reg401+reg474; T reg689=reg111*reg237; reg605=reg560+reg605; T reg690=reg72*reg233; T reg691=reg72*reg222;
    T reg692=reg173*reg203; T reg693=reg162*reg206; T reg694=reg185*reg162; T reg695=reg195*reg154; T reg696=reg246*reg38;
    T reg697=reg38*reg240; reg575=reg539+reg575; T reg698=reg152*reg183; T reg699=reg38*reg217; T reg700=reg111*reg244;
    reg349=reg349+reg345; T reg701=reg184*reg162; T reg702=reg72*reg360; T reg703=reg202*reg173; T reg704=reg157*reg202;
    reg365=reg289+reg365; T reg705=reg154*reg183; T reg706=reg38*reg259; T reg707=reg152*reg179; T reg708=reg255*reg72;
    T reg709=reg38*reg222; reg579=reg645+reg579; reg532=reg533+reg532; T reg710=reg72*reg375; reg380=reg378+reg380;
    T reg711=reg243*reg77; T reg712=reg142*reg212; T reg713=reg173*reg200; T reg714=reg154*reg179; T reg715=reg38*reg228;
    T reg716=reg138*reg229; T reg717=reg152*reg200; T reg718=reg161*reg200; reg630=reg635+reg630; T reg719=reg187*reg168;
    reg308=reg636+reg308; T reg720=reg138*reg375; T reg721=reg38*reg244; T reg722=reg72*reg238; T reg723=reg243*reg173;
    T reg724=reg152*reg174; T reg725=reg38*reg212; T reg726=reg569+reg568; T reg727=reg185*reg77; reg290=reg290-reg473;
    reg536=reg231+reg536; reg299=reg302+reg299; T reg728=reg149*reg171; T reg729=reg154*reg174; T reg730=reg38*reg215;
    T reg731=reg152*reg185; T reg732=reg38*reg220; reg570=reg555+reg570; T reg733=reg127*reg208; T reg734=reg72*reg279;
    T reg735=reg142*reg279; T reg736=reg201*reg173; T reg737=reg111*reg220; reg333=reg335+reg333; reg272=reg307+reg272;
    T reg738=reg154*reg185; T reg739=reg38*reg242; T reg740=reg152*reg206; T reg741=reg38*reg237; reg636=reg252-reg636;
    reg252=reg157*reg201; T reg742=reg253*reg72; T reg743=reg280-reg281; T reg744=reg154*reg206; reg546=reg549-reg546;
    reg402=reg405+reg402; reg219=reg623-reg219; T reg745=reg119*reg375; T reg746=reg199*reg161; T reg747=reg562*reg133;
    T reg748=reg199*reg154; T reg749=reg172*reg200; T reg750=reg151*reg200; T reg751=reg246*reg133; T reg752=reg178*reg162;
    T reg753=reg177*reg162; T reg754=reg154*reg202; T reg755=reg149*reg187; T reg756=reg127*reg601; T reg757=reg149*reg189;
    reg457=reg462+reg457; reg389=reg392+reg389; reg462=reg201*reg162; T reg758=reg161*reg204; T reg759=reg133*reg553;
    T reg760=reg154*reg204; T reg761=reg133*reg225; T reg762=reg154*reg201; T reg763=reg127*reg518; T reg764=reg179*reg162;
    T reg765=reg111*reg228; reg435=reg436+reg435; T reg766=reg161*reg189; T reg767=reg133*reg601; T reg768=reg149*reg188;
    T reg769=reg162*reg193; T reg770=reg127*reg215; reg550=reg554+reg550; reg263=reg266+reg263; reg266=reg149*reg202;
    reg496=reg501+reg496; reg226=reg230+reg226; reg501=reg149*reg181; T reg771=reg127*reg504; T reg772=reg588*reg117;
    T reg773=reg149*reg186; T reg774=reg142*reg465; T reg775=reg142*reg360; reg269=reg642+reg269; T reg776=reg243*reg162;
    T reg777=reg151*reg177; reg484=reg366+reg484; T reg778=reg142*reg233; T reg779=reg149*reg203; T reg780=reg475+reg505;
    T reg781=reg142*reg375; reg478=reg481+reg478; reg481=reg149*reg200; T reg782=reg151*reg202; T reg783=reg149*reg178;
    T reg784=reg142*reg497; T reg785=reg255*reg119; reg510=reg510+reg274; reg410=reg287+reg410; T reg786=reg154*reg200;
    reg273=reg284+reg273; reg624=reg627+reg624; reg284=reg648+reg363; T reg787=reg133*reg518; T reg788=reg154*reg187;
    T reg789=reg183*reg77; reg642=reg612-reg642; reg472=reg291+reg472; reg291=reg133*reg242; reg612=reg154*reg203;
    T reg790=reg127*reg562; T reg791=reg149*reg199; reg340=reg340-reg344; T reg792=reg200*reg162; reg596=reg597-reg596;
    T reg793=reg161*reg171; T reg794=reg133*reg208; T reg795=reg154*reg171; reg543=reg544+reg543; T reg796=reg133*reg240;
    T reg797=reg154*reg186; T reg798=reg111*reg246; reg598=reg572+reg598; reg314=reg315+reg314; reg315=reg72*reg497;
    T reg799=reg178*reg173; reg423=reg422+reg423; T reg800=reg161*reg181; T reg801=reg111*reg217; T reg802=reg157*reg178;
    T reg803=reg133*reg504; T reg804=reg154*reg181; T reg805=reg133*reg259; T reg806=reg127*reg588; reg520=reg521-reg520;
    T reg807=reg154*reg189; reg287=reg288-reg287; reg288=reg133*reg244; T reg808=reg154*reg243; reg524=reg527+reg524;
    reg466=reg327+reg466; reg327=reg427-reg426; T reg809=reg179*reg77; reg618=reg619+reg618; T reg810=reg111*reg229;
    T reg811=reg111*reg222; T reg812=reg161*reg193; T reg813=reg133*reg209; reg260=reg639+reg260; reg639=reg127*reg553;
    T reg814=reg119*reg220; T reg815=reg149*reg204; T reg816=reg183*reg162; T reg817=reg202*reg162; reg587=reg589-reg587;
    T reg818=reg142*reg310; T reg819=reg416+reg381; T reg820=reg195*reg77; T reg821=reg154*reg177; T reg822=reg111*reg259;
    reg366=reg368-reg366; reg590=reg595+reg590; reg368=reg149*reg177; reg450=reg449+reg450; T reg823=reg161*reg187;
    T reg824=reg625-reg626; T reg825=reg153*reg189; T reg826=reg145*reg244; T reg827=reg145*reg601; T reg828=reg155*reg189;
    T reg829=reg153*reg201; reg245=reg622+reg245; T reg830=reg153*reg204; T reg831=reg145*reg225; T reg832=reg145*reg553;
    T reg833=reg155*reg204; T reg834=reg153*reg202; reg223=reg615+reg223; T reg835=reg153*reg199; T reg836=reg145*reg246;
    T reg837=reg145*reg562; T reg838=reg155*reg199; T reg839=reg153*reg200; reg250=reg641+reg250; reg224=reg329+reg224;
    T reg840=reg119*reg497; T reg841=reg178*reg172; T reg842=reg151*reg178; T reg843=reg119*reg222; reg227=reg321+reg227;
    T reg844=reg119*reg465; T reg845=reg145*reg588; T reg846=reg155*reg188; T reg847=reg153*reg178; reg541=reg540+reg541;
    T reg848=reg153*reg181; T reg849=reg145*reg259; T reg850=reg145*reg504; T reg851=reg155*reg181; T reg852=reg153*reg186;
    reg535=reg534+reg535; T reg853=reg153*reg171; T reg854=reg145*reg240; T reg855=reg145*reg208; T reg856=reg155*reg171;
    T reg857=reg153*reg203; reg205=reg205-reg647; T reg858=reg153*reg187; T reg859=reg145*reg242; T reg860=reg145*reg518;
    T reg861=reg155*reg187; T reg862=reg153*reg177; reg211=reg644+reg211; T reg863=reg646+reg643; T reg864=reg145*reg209;
    T reg865=reg155*reg193; T reg866=reg153*reg243; T reg867=reg259*reg135; T reg868=reg168*reg183; reg322=reg321+reg322;
    reg321=reg135*reg217; T reg869=reg151*reg183; T reg870=reg135*reg240; T reg871=reg168*reg206; T reg872=reg377-reg379;
    T reg873=reg135*reg237; T reg874=reg151*reg206; T reg875=reg135*reg242; T reg876=reg185*reg168; reg374=reg373+reg374;
    T reg877=reg135*reg220; T reg878=reg151*reg185; T reg879=reg215*reg135; T reg880=reg168*reg174; reg364=reg364-reg363;
    T reg881=reg135*reg212; T reg882=reg151*reg174; T reg883=reg244*reg135; T reg884=reg184*reg168; reg359=reg357+reg359;
    T reg885=reg253*reg135; T reg886=reg184*reg151; T reg887=reg135*reg225; T reg888=reg172*reg186; T reg889=reg151*reg186;
    T reg890=reg119*reg217; reg190=reg377+reg190; reg377=reg119*reg233; T reg891=reg172*reg203; T reg892=reg151*reg203;
    T reg893=reg119*reg237; reg236=reg373+reg236; reg373=reg119*reg310; T reg894=reg177*reg172; T reg895=reg633+reg634;
    reg350=reg357+reg350; reg357=reg119*reg279; T reg896=reg201*reg172; T reg897=reg201*reg151; T reg898=reg253*reg119;
    reg336=reg285+reg336; T reg899=reg119*reg360; T reg900=reg172*reg202; T reg901=reg119*reg229; T reg902=reg135*reg228;
    T reg903=reg179*reg168; reg330=reg329+reg330; reg329=reg135*reg222; T reg904=reg151*reg179; T reg905=reg155*reg178;
    T reg906=reg150*reg178; T reg907=reg144*reg222; reg559=reg525+reg559; T reg908=reg144*reg465; T reg909=reg155*reg186;
    T reg910=reg150*reg186; T reg911=reg144*reg217; reg610=reg517+reg610; T reg912=reg144*reg233; T reg913=reg155*reg203;
    T reg914=reg150*reg203; T reg915=reg144*reg237; reg608=reg512+reg608; T reg916=reg144*reg310; T reg917=reg155*reg177;
    T reg918=reg150*reg177; T reg919=reg144*reg220; T reg920=reg506+reg603; T reg921=reg144*reg238; T reg922=reg155*reg243;
    T reg923=reg599+reg600; reg241=reg558+reg241; T reg924=reg144*reg279; T reg925=reg155*reg201; T reg926=reg150*reg201;
    T reg927=reg144*reg251; T reg928=reg132*reg249; T reg929=reg145*reg249; T reg930=reg119*reg251; T reg931=reg135*reg249;
    T reg932=reg117*reg249; T reg933=reg130*reg251; T reg934=reg578-reg529; T reg935=reg104*reg249; T reg936=reg138*reg251;
    T reg937=reg38*reg249; T reg938=reg133*reg249; T reg939=reg573-reg514; T reg940=reg111*reg249; T reg941=reg571-reg545;
    T reg942=reg72*reg251; T reg943=reg121*reg249; T reg944=reg116*reg249; T reg945=reg118*reg251; T reg946=reg148*reg249;
    T reg947=reg131*reg249; T reg948=reg125*reg251; T reg949=reg137*reg249; T reg950=reg120*reg249; reg561=reg530+reg561;
    T reg951=reg144*reg497; T reg952=reg153*reg185; reg513=reg512+reg513; reg512=reg132*reg220; T reg953=reg150*reg185;
    T reg954=reg132*reg215; T reg955=reg153*reg174; reg509=reg509-reg506; T reg956=reg132*reg212; T reg957=reg150*reg174;
    T reg958=reg132*reg244; T reg959=reg153*reg184; reg503=reg558+reg503; reg558=reg132*reg253; T reg960=reg150*reg184;
    T reg961=reg132*reg225; T reg962=reg153*reg180; reg552=reg551+reg552; T reg963=reg132*reg255; T reg964=reg150*reg180;
    T reg965=reg132*reg246; T reg966=reg153*reg195; reg548=reg547+reg548; T reg967=reg132*reg229; T reg968=reg150*reg195;
    T reg969=reg153*reg188; T reg970=reg145*reg228; T reg971=reg144*reg253; reg207=reg551+reg207; reg551=reg144*reg360;
    T reg972=reg155*reg202; T reg973=reg150*reg202; T reg974=reg144*reg255; reg210=reg547+reg210; reg547=reg144*reg375;
    T reg975=reg155*reg200; T reg976=reg150*reg200; T reg977=reg144*reg229; T reg978=reg132*reg228; T reg979=reg153*reg179;
    reg216=reg530+reg216; reg530=reg132*reg222; T reg980=reg150*reg179; T reg981=reg132*reg259; T reg982=reg153*reg183;
    reg526=reg525+reg526; reg525=reg132*reg217; T reg983=reg150*reg183; T reg984=reg132*reg240; T reg985=reg153*reg206;
    reg517=reg517-reg519; T reg986=reg132*reg237; T reg987=reg150*reg206; T reg988=reg132*reg242; T reg989=reg104*reg228;
    T reg990=reg143*reg188; T reg991=reg104*reg588; T reg992=reg159*reg178; reg383=reg383-reg382; T reg993=reg159*reg181;
    T reg994=reg104*reg259; T reg995=reg143*reg181; T reg996=reg104*reg504; T reg997=reg159*reg186; reg386=reg386-reg385;
    T reg998=reg159*reg171; T reg999=reg104*reg240; reg404=reg406+reg404; T reg1000=reg159*reg203; reg485=reg485+reg480;
    T reg1001=reg159*reg187; T reg1002=reg104*reg242; T reg1003=reg143*reg187; T reg1004=reg104*reg518; T reg1005=reg159*reg177;
    reg499=reg499-reg493; T reg1006=reg502+reg213; T reg1007=reg143*reg193; T reg1008=reg104*reg209; T reg1009=reg448+reg447;
    T reg1010=reg139*reg242; T reg1011=reg159*reg185; reg459=reg464-reg459; T reg1012=reg139*reg220; T reg1013=reg147*reg185;
    T reg1014=reg139*reg215; T reg1015=reg159*reg174; T reg1016=reg447+reg221; T reg1017=reg139*reg212; T reg1018=reg147*reg174;
    T reg1019=reg139*reg244; T reg1020=reg159*reg184; reg334=reg367-reg334; T reg1021=reg253*reg139; T reg1022=reg184*reg147;
    T reg1023=reg139*reg225; T reg1024=reg159*reg180; reg431=reg413-reg431; T reg1025=reg255*reg139; T reg1026=reg147*reg180;
    T reg1027=reg246*reg139; T reg1028=reg195*reg159; reg425=reg432-reg425; T reg1029=reg139*reg229; T reg1030=reg195*reg147;
    T reg1031=reg159*reg188; T reg1032=reg138*reg217; reg283=reg280+reg283; reg280=reg138*reg233; T reg1033=reg161*reg203;
    T reg1034=reg152*reg203; T reg1035=reg138*reg237; reg358=reg307+reg358; reg307=reg138*reg310; T reg1036=reg161*reg177;
    T reg1037=reg152*reg177; T reg1038=reg138*reg220; T reg1039=reg473+reg371; T reg1040=reg138*reg238; T reg1041=reg161*reg243;
    T reg1042=reg325+reg323; reg339=reg451+reg339; reg451=reg138*reg279; T reg1043=reg161*reg201; T reg1044=reg152*reg201;
    T reg1045=reg253*reg138; reg341=reg490+reg341; reg490=reg138*reg360; T reg1046=reg161*reg202; T reg1047=reg152*reg202;
    T reg1048=reg255*reg138; reg347=reg401+reg347; reg401=reg159*reg189; T reg1049=reg104*reg244; T reg1050=reg143*reg189;
    T reg1051=reg104*reg601; T reg1052=reg159*reg201; reg461=reg461-reg460; T reg1053=reg159*reg204; T reg1054=reg104*reg225;
    T reg1055=reg143*reg204; T reg1056=reg104*reg553; T reg1057=reg159*reg202; reg298=reg298-reg294; T reg1058=reg199*reg159;
    T reg1059=reg246*reg104; T reg1060=reg199*reg143; T reg1061=reg562*reg104; T reg1062=reg159*reg200; reg261=reg261-reg317;
    reg277=reg378+reg277; reg378=reg138*reg497; T reg1063=reg161*reg178; T reg1064=reg152*reg178; T reg1065=reg138*reg222;
    reg278=reg289+reg278; reg289=reg138*reg465; T reg1066=reg161*reg186; T reg1067=reg152*reg186; T reg1068=reg117*reg242;
    T reg1069=reg518*reg117; T reg1070=reg187*reg172; T reg1071=reg177*reg168; reg304=reg303+reg304; T reg1072=reg300+reg301;
    T reg1073=reg117*reg209; T reg1074=reg172*reg193; T reg1075=reg243*reg168; T reg1076=reg295-reg296; T reg1077=reg189*reg168;
    T reg1078=reg244*reg117; T reg1079=reg601*reg117; T reg1080=reg189*reg172; T reg1081=reg201*reg168; reg468=reg467+reg468;
    T reg1082=reg168*reg204; T reg1083=reg117*reg225; T reg1084=reg117*reg553; T reg1085=reg172*reg204; T reg1086=reg168*reg202;
    reg453=reg452+reg453; T reg1087=reg199*reg168; T reg1088=reg246*reg117; T reg1089=reg562*reg117; T reg1090=reg199*reg172;
    T reg1091=reg168*reg180; reg286=reg285+reg286; reg285=reg255*reg135; T reg1092=reg151*reg180; T reg1093=reg246*reg135;
    T reg1094=reg195*reg168; reg275=reg274+reg275; reg274=reg135*reg229; T reg1095=reg195*reg151; T reg1096=reg188*reg168;
    T reg1097=reg117*reg228; T reg1098=reg188*reg172; T reg1099=reg178*reg168; reg265=reg264+reg265; T reg1100=reg168*reg181;
    T reg1101=reg259*reg117; T reg1102=reg117*reg504; T reg1103=reg172*reg181; T reg1104=reg168*reg186; reg318=reg316+reg318;
    T reg1105=reg168*reg171; T reg1106=reg117*reg240; T reg1107=reg117*reg208; T reg1108=reg172*reg171; T reg1109=reg168*reg203;
    reg311=reg311-reg309; T reg1110=reg130*reg279; T reg1111=reg201*reg147; T reg1112=reg253*reg130; reg388=reg413+reg388;
    reg413=reg143*reg202; T reg1113=reg130*reg360; T reg1114=reg147*reg202; T reg1115=reg255*reg130; reg387=reg432+reg387;
    reg432=reg143*reg200; T reg1116=reg130*reg375; T reg1117=reg147*reg200; T reg1118=reg130*reg229; T reg1119=reg139*reg228;
    T reg1120=reg159*reg179; reg429=reg428-reg429; T reg1121=reg139*reg222; T reg1122=reg179*reg147; T reg1123=reg139*reg259;
    T reg1124=reg159*reg183; reg424=reg421-reg424; T reg1125=reg139*reg217; T reg1126=reg147*reg183; T reg1127=reg419+reg415;
    T reg1128=reg400+reg395; reg403=reg445+reg403; T reg1129=reg168*reg200; reg500=reg498+reg500; reg495=reg428+reg495;
    reg428=reg143*reg178; T reg1130=reg130*reg497; T reg1131=reg178*reg147; T reg1132=reg130*reg222; reg488=reg421+reg488;
    reg421=reg143*reg186; T reg1133=reg130*reg465; T reg1134=reg147*reg186; T reg1135=reg130*reg217; reg477=reg400+reg477;
    reg476=reg411+reg476; reg400=reg147*reg203; T reg1136=reg130*reg237; reg409=reg464+reg409; reg464=reg143*reg177;
    T reg1137=reg130*reg310; T reg1138=reg147*reg177; T reg1139=reg130*reg220; reg247=reg247+reg399; T reg1140=reg143*reg243;
    T reg1141=reg130*reg238; T reg1142=reg396+reg397; reg394=reg367+reg394; reg367=reg143*reg201; reg614=reg613+reg614;
    T reg1143=reg195*reg176; T reg1144=reg148*reg217; T reg1145=reg164*reg183; T reg1146=reg246*reg116; T reg1147=reg131*reg588;
    T reg1148=reg188*reg173; T reg1149=reg120*reg259; T reg1150=reg199*reg156; T reg1151=reg167*reg200; T reg1152=reg588*reg116;
    T reg1153=reg163*reg181; T reg1154=reg137*reg229; T reg1155=reg125*reg375; T reg1156=reg418+reg417; T reg1157=reg156*reg188;
    T reg1158=reg116*reg228; T reg1159=reg118*reg279; T reg1160=reg585*(*f.m).f_vol[2]; T reg1161=reg240*reg148; T reg1162=reg131*reg504;
    T reg1163=reg176*reg186; T reg1164=reg585*(*f.m).f_vol[1]; T reg1165=reg160*reg206; T reg1166=reg259*reg148; T reg1167=reg163*reg179;
    T reg1168=reg137*reg228; T reg1169=reg160*reg183; T reg1170=reg181*reg173; reg588=reg120*reg588; reg471=reg469+reg471;
    T reg1171=reg504*reg116; T reg1172=reg164*reg201; T reg1173=reg120*reg228; T reg1174=reg167*reg188; T reg1175=reg156*reg202;
    reg414=reg414+reg412; T reg1176=reg156*reg181; T reg1177=reg160*reg199; T reg1178=reg163*reg174; T reg1179=reg259*reg116;
    reg342=reg343+reg342; T reg1180=reg125*reg229; T reg1181=reg160*reg203; reg593=reg594+reg593; T reg1182=reg163*reg178;
    T reg1183=reg156*reg178; T reg1184=reg176*reg200; T reg1185=reg201*reg158; T reg1186=reg164*reg206; T reg1187=reg193*reg158;
    T reg1188=reg157*reg180; T reg1189=reg255*reg121; T reg1190=reg176*reg178; T reg1191=reg131*reg209; T reg1192=reg156*reg200;
    reg606=reg607+reg606; reg542=reg539+reg542; reg539=reg137*reg212; T reg1193=reg195*reg163; T reg1194=reg246*reg137;
    T reg1195=reg167*reg202; T reg1196=reg163*reg201; T reg1197=reg583*(*f.m).f_vol[2]; T reg1198=reg242*reg148; T reg1199=reg160*reg185;
    T reg1200=reg131*reg208; T reg1201=reg156*reg180; T reg1202=reg225*reg121; T reg1203=reg176*reg174; T reg1204=reg583*(*f.m).f_vol[1];
    T reg1205=reg583*(*f.m).f_vol[0]; T reg1206=reg125*reg360; T reg1207=reg441+reg440; reg191=reg191-reg638; T reg1208=reg195*reg157;
    T reg1209=reg121*reg229; reg370=reg538+reg370; T reg1210=reg188*reg158; T reg1211=reg562*reg116; T reg1212=reg199*reg173;
    reg408=reg408+reg407; T reg1213=reg125*reg222; T reg1214=reg120*reg244; reg232=reg645+reg232; reg645=reg163*reg189;
    T reg1215=reg337-reg338; reg537=reg538+reg537; reg538=reg125*reg217; T reg1216=reg255*reg125; reg384=reg463+reg384;
    T reg1217=reg195*reg156; T reg1218=reg246*reg121; T reg1219=reg160*reg202; T reg1220=reg120*reg601; T reg1221=reg237*reg148;
    reg494=reg489+reg494; T reg1222=reg167*reg189; T reg1223=reg176*reg202; T reg1224=reg116*reg209; reg375=reg118*reg375;
    T reg1225=reg160*reg181; T reg1226=reg131*reg562; reg507=reg508+reg507; T reg1227=reg163*reg186; reg313=reg312+reg313;
    T reg1228=reg200*reg158; T reg1229=reg361+reg362; T reg1230=reg156*reg201; reg360=reg118*reg360; T reg1231=reg163*reg204;
    T reg1232=reg120*reg225; T reg1233=reg164*reg200; T reg1234=reg120*reg240; T reg1235=reg156*reg177; reg372=reg369+reg372;
    T reg1236=reg163*reg171; T reg1237=reg163*reg183; T reg1238=reg118*reg229; T reg1239=reg137*reg259; T reg1240=reg199*reg158;
    T reg1241=reg176*reg185; T reg1242=reg225*reg116; T reg1243=reg187*reg173; T reg1244=reg164*reg202; T reg1245=reg131*reg259;
    T reg1246=reg118*reg255; T reg1247=reg601*reg116; T reg1248=reg189*reg173; T reg1249=reg156*reg189; T reg1250=reg244*reg116;
    reg580=reg580+reg557; reg631=reg631-reg632; T reg1251=reg137*reg217; reg455=reg348+reg455; T reg1252=reg160*reg187;
    T reg1253=reg163*reg202; reg305=reg482+reg305; reg504=reg120*reg504; T reg1254=reg156*reg243; T reg1255=reg276-reg282;
    T reg1256=reg167*reg204; T reg1257=reg131*reg242; T reg1258=reg120*reg553; T reg1259=reg167*reg181; T reg1260=reg137*reg240;
    reg202=reg202*reg158; reg346=reg348+reg346; reg348=reg173*reg193; reg522=reg523+reg522; reg320=reg320-reg353;
    T reg1261=reg173*reg171; T reg1262=reg116*reg208; reg324=reg326+reg324; T reg1263=reg163*reg177; T reg1264=reg148*reg222;
    T reg1265=reg204*reg173; T reg1266=reg131*reg246; T reg1267=reg156*reg171; T reg1268=reg240*reg116; reg437=reg437+reg438;
    T reg1269=reg137*reg244; T reg1270=reg164*reg179; T reg1271=reg137*reg215; T reg1272=reg118*reg253; T reg1273=reg528+reg640;
    T reg1274=reg256*(*f.m).f_vol[2]; T reg1275=reg160*reg178; T reg1276=reg167*reg186; T reg1277=reg156*reg186; reg235=reg637+reg235;
    T reg1278=reg160*reg177; T reg1279=reg163*reg188; reg487=reg487+reg486; T reg1280=reg518*reg116; T reg1281=reg228*reg148;
    T reg1282=reg120*reg208; T reg1283=reg257*(*f.m).f_vol[1]; reg443=reg444+reg443; T reg1284=reg187*reg158; T reg1285=reg585*(*f.m).f_vol[0];
    T reg1286=reg125*reg465; T reg1287=reg156*reg187; T reg1288=reg116*reg242; T reg1289=reg160*reg179; T reg1290=reg156*reg204;
    T reg1291=reg160*reg200; T reg1292=reg176*reg179; T reg1293=reg131*reg518; reg518=reg120*reg518; T reg1294=reg167*reg187;
    T reg1295=reg137*reg222; T reg1296=reg156*reg203; reg332=reg332-reg331; T reg1297=reg584*(*f.m).f_vol[1]; T reg1298=reg584*(*f.m).f_vol[2];
    reg181=reg181*reg158; reg351=reg356+reg351; T reg1299=reg553*reg116; T reg1300=reg160*reg184; T reg1301=reg240*reg121;
    T reg1302=reg125*reg238; T reg1303=reg176*reg184; T reg1304=reg629-reg628; T reg1305=reg125*reg233; T reg1306=reg160*reg204;
    reg240=reg131*reg240; reg456=reg456+reg454; T reg1307=reg157*reg183; T reg1308=reg121*reg217; reg470=reg463+reg470;
    reg463=reg164*reg186; reg217=reg118*reg217; T reg1309=reg638+reg319; T reg1310=reg137*reg225; T reg1311=reg163*reg180;
    T reg1312=reg160*reg171; T reg1313=reg163*reg206; reg218=reg560+reg218; reg560=reg253*reg148; T reg1314=reg164*reg184;
    T reg1315=reg177*reg158; T reg1316=reg167*reg203; reg393=reg343+reg393; reg343=reg176*reg206; T reg1317=reg167*reg193;
    T reg1318=reg137*reg237; T reg1319=reg157*reg206; T reg1320=reg121*reg237; reg566=reg567+reg566; reg465=reg118*reg465;
    T reg1321=reg176*reg183; T reg1322=reg212*reg148; T reg1323=reg164*reg174; T reg1324=reg160*reg189; T reg1325=reg118*reg220;
    T reg1326=reg268+reg267; reg602=reg602-reg604; reg398=reg398+reg511; T reg1327=reg186*reg158; T reg1328=reg167*reg243;
    T reg1329=reg131*reg225; T reg1330=reg164*reg177; T reg1331=reg253*reg137; T reg1332=reg156*reg206; T reg1333=reg244*reg148;
    T reg1334=reg131*reg244; reg574=reg572+reg574; reg391=reg391+reg390; reg572=reg176*reg180; reg292=reg354+reg292;
    T reg1335=reg255*reg148; T reg1336=reg164*reg180; T reg1337=reg125*reg310; T reg1338=reg125*reg497; T reg1339=reg156*reg179;
    T reg1340=reg121*reg228; reg306=reg531+reg306; T reg1341=reg137*reg220; T reg1342=reg163*reg203; T reg1343=reg246*reg148;
    T reg1344=reg584*(*f.m).f_vol[0]; reg620=reg620-reg621; T reg1345=reg72*reg229; T reg1346=reg160*reg201; T reg1347=reg176*reg203;
    T reg1348=reg160*reg195; T reg1349=reg157*reg200; T reg1350=reg125*reg237; T reg1351=reg164*reg203; reg237=reg118*reg237;
    reg310=reg118*reg310; T reg1352=reg125*reg220; T reg1353=reg617-reg616; reg434=reg337+reg434; reg183=reg156*reg183;
    reg337=reg163*reg243; T reg1354=reg176*reg177; reg259=reg259*reg121; reg225=reg225*reg148; reg576=reg577+reg576;
    T reg1355=reg167*reg171; reg180=reg160*reg180; reg189=reg189*reg158; reg179=reg157*reg179; reg601=reg131*reg601;
    reg187=reg163*reg187; reg483=reg482+reg483; reg482=reg121*reg222; reg446=reg444+reg446; reg444=reg120*reg242;
    reg177=reg167*reg177; T reg1356=reg167*reg178; reg255=reg255*reg137; T reg1357=reg118*reg233; reg203=reg203*reg158;
    reg516=reg516-reg515; T reg1358=reg160*reg243; reg246=reg246*reg120; T reg1359=reg199*reg163; T reg1360=(*f.m).f_vol[0]*reg581;
    T reg1361=reg215*reg148; T reg1362=reg156*reg174; T reg1363=reg581*(*f.m).f_vol[1]; T reg1364=reg215*reg121; T reg1365=reg430-reg433;
    T reg1366=reg160*reg174; T reg1367=reg156*reg184; reg293=reg629+reg293; reg629=reg164*reg178; reg234=reg531+reg234;
    reg204=reg204*reg158; reg531=reg167*reg201; reg222=reg118*reg222; reg562=reg562*reg120; T reg1368=reg586*(*f.m).f_vol[1];
    reg188=reg160*reg188; T reg1369=reg164*reg185; reg174=reg157*reg174; T reg1370=reg121*reg212; T reg1371=reg586*(*f.m).f_vol[0];
    T reg1372=reg253*reg125; T reg1373=reg586*(*f.m).f_vol[2]; reg229=reg229*reg148; reg492=reg326+reg492; reg201=reg176*reg201;
    reg326=reg220*reg148; reg497=reg118*reg497; reg228=reg131*reg228; reg553=reg131*reg553; reg244=reg244*reg121;
    reg178=reg178*reg158; T reg1374=reg243*reg158; reg238=reg118*reg238; T reg1375=reg163*reg184; T reg1376=reg582*(*f.m).f_vol[2];
    T reg1377=reg632+reg420; reg563=reg564+reg563; T reg1378=reg186*reg162; reg195=reg195*reg164; T reg1379=reg121*reg242;
    T reg1380=reg156*reg185; T reg1381=(*f.m).f_vol[0]*reg582; reg200=reg163*reg200; reg209=reg120*reg209; reg171=reg171*reg158;
    T reg1382=reg581*(*f.m).f_vol[2]; reg355=reg354+reg355; reg354=reg582*(*f.m).f_vol[1]; reg479=reg356+reg479; reg253=reg253*reg121;
    reg270=reg567+reg270; reg356=reg157*reg185; reg220=reg121*reg220; reg352=reg577+reg352; reg556=reg555+reg556;
    reg271=reg262-reg271; reg555=reg256*(*f.m).f_vol[0]; reg185=reg163*reg185; reg279=reg125*reg279; reg242=reg137*reg242;
    reg199=reg199*reg167; reg186=reg160*reg186; reg592=reg591+reg592; reg184=reg157*reg184; reg384=reg390+reg384;
    reg360=reg202+reg360; reg450=reg821+reg450; reg1087=reg1088+reg1087; reg1244=reg1246+reg1244; reg1260=reg1260-reg1313;
    reg202=reg75*reg1207; reg195=reg229+reg195; reg840=reg841+reg840; reg458=reg762+reg458; reg292=reg412+reg292;
    reg442=reg654+reg442; reg425=reg1062+reg425; reg832=reg833+reg832; reg229=reg75*reg895; reg500=reg500+reg1129;
    reg373=reg894+reg373; reg227=reg316+reg227; reg620=reg1342+reg620; reg805=reg804+reg805; reg666=reg665+reg666;
    reg316=reg75*reg819; reg1024=reg1024-reg1023; reg842=reg843+reg842; reg835=reg836+reg835; reg1330=reg1325+reg1330;
    reg651=reg649+reg651; reg456=reg456+reg186; reg683=reg682+reg683; reg889=reg890+reg889; reg340=reg612+reg340;
    reg1190=reg1213+reg1190; reg837=reg838+reg837; reg209=reg209-reg1317; reg191=reg191-reg337; reg891=reg891-reg377;
    reg1028=reg1028-reg1027; reg796=reg795+reg796; reg1089=reg1090+reg1089; reg430=reg430-reg1377; reg793=reg793-reg794;
    reg696=reg695+reg696; reg1172=reg1272+reg1172; reg431=reg1057+reg431; reg474=reg786+reg474; reg250=reg250+reg839;
    reg190=reg190-reg309; reg223=reg223+reg834; reg803=reg800+reg803; reg1241=reg1341+reg1241; reg674=reg673+reg674;
    reg660=reg658+reg660; reg310=reg1315+reg310; reg236=reg303+reg236; reg1203=reg1203-reg539; reg787=reg823+reg787;
    reg181=reg1162+reg181; reg1025=reg1026-reg1025; reg491=reg754+reg491; reg224=reg264+reg224; reg291=reg788+reg291;
    reg844=reg888+reg844; reg679=reg678+reg679; reg423=reg797+reg423; reg443=reg407+reg443; reg892=reg893+reg892;
    reg238=reg238-reg1374; reg1353=reg1353-reg337; reg1159=reg1185+reg1159; reg1178=reg1178-reg1271; reg293=reg293-reg621;
    reg289=reg1066+reg289; reg1093=reg1094+reg1093; reg1077=reg1078+reg1077; reg278=reg422+reg278; reg275=reg1129+reg275;
    reg631=reg631-reg1358; reg279=reg531+reg279; reg1064=reg1065+reg1064; reg998=reg999+reg998; reg274=reg1095+reg274;
    reg271=reg271-reg1378; reg378=reg1063+reg378; reg277=reg439+reg277; reg1096=reg1097+reg1096; reg270=reg607+reg270;
    reg1316=reg1316-reg1305; reg1098=reg772+reg1098; reg1323=reg1323-reg1322; reg1062=reg261+reg1062; reg261=reg75*reg404;
    reg883=reg884+reg883; reg1037=reg1038+reg1037; reg359=reg1081+reg359; reg1206=reg1195+reg1206; reg307=reg1036+reg307;
    reg1079=reg1080+reg1079; reg355=reg1278+reg355; reg358=reg449+reg358; reg885=reg886+reg885; reg386=reg386+reg997;
    reg887=reg1091+reg887; reg1034=reg1035+reg1034; reg352=reg557+reg352; reg326=reg1369+reg326; reg1033=reg1033-reg280;
    reg286=reg1086+reg286; reg283=reg283-reg344; reg188=reg228+reg188; reg285=reg1092+reg285; reg1366=reg1366-reg1361;
    reg1067=reg1032+reg1067; reg201=reg1372+reg201; reg1108=reg1108-reg1107; reg225=reg180+reg225; reg401=reg1049+reg401;
    reg1354=reg1352+reg1354; reg311=reg311+reg1109; reg446=reg1219+reg446; reg180=reg75*reg1009; reg177=reg1337+reg177;
    reg1068=reg719+reg1068; reg1008=reg1008+reg1007; reg1347=reg1350+reg1347; reg1001=reg1002+reg1001; reg1335=reg1336+reg1335;
    reg228=reg75*reg1006; reg1069=reg1070+reg1069; reg483=reg1291+reg483; reg304=reg304+reg1071; reg1343=reg1348+reg1343;
    reg499=reg499+reg1005; reg306=reg523+reg306; reg1003=reg1004-reg1003; reg264=reg75*reg1072; reg265=reg265+reg1099;
    reg1060=reg1061-reg1060; reg303=reg75*reg1326; reg390=reg75*reg398; reg1058=reg1059+reg1058; reg1100=reg1101+reg1100;
    reg1333=reg1300+reg1333; reg1057=reg298+reg1057; reg1102=reg1103+reg1102; reg1302=reg1302-reg1328; reg470=reg1346+reg470;
    reg1055=reg1056-reg1055; reg1076=reg1076-reg1075; reg318=reg318+reg1104; reg1053=reg1054+reg1053; reg617=reg617-reg1309;
    reg485=reg485+reg1000; reg560=reg1314+reg560; reg461=reg461+reg1052; reg1105=reg1106+reg1105; reg1073=reg1073-reg1074;
    reg1050=reg1051-reg1050; reg1031=reg989+reg1031; reg902=reg903+reg902; reg1239=reg1237+reg1239; reg1281=reg1289+reg1281;
    reg743=reg612+reg743; reg1286=reg1276+reg1286; reg330=reg1099+reg330; reg697=reg697-reg744; reg1084=reg1085+reg1084;
    reg351=reg1275+reg351; reg699=reg698+reg699; reg329=reg904+reg329; reg1295=reg1292+reg1295; reg990=reg991-reg990;
    reg365=reg797+reg365; reg867=reg868+reg867; reg1264=reg1270+reg1264; reg706=reg705+reg706; reg322=reg1104+reg322;
    reg324=reg1182+reg324; reg709=reg707+reg709; reg721=reg669+reg721; reg455=reg508+reg455; reg350=reg467+reg350;
    reg724=reg724-reg725; reg305=reg486+reg305; reg290=reg290-reg808; reg357=reg896+reg357; reg1251=reg1321+reg1251;
    reg1029=reg1030-reg1029; reg729=reg729-reg730; reg897=reg898+reg897; reg375=reg1228+reg375; reg732=reg731+reg732;
    reg336=reg452+reg336; reg346=reg1227+reg346; reg1225=reg1245+reg1225; reg272=reg821+reg272; reg1086=reg453+reg1086;
    reg899=reg900+reg899; reg1233=reg1238+reg1233; reg739=reg738+reg739; reg741=reg741-reg740; reg1161=reg1161-reg1165;
    reg341=reg392+reg341; reg1155=reg1151+reg1155; reg993=reg994+reg993; reg1044=reg1045+reg1044; reg877=reg878+reg877;
    reg1215=reg1181+reg1215; reg451=reg1043+reg451; reg880=reg880-reg879; reg370=reg564+reg370; reg1210=reg1147+reg1210;
    reg339=reg436+reg339; reg364=reg364-reg1075; reg1221=reg1221-reg1186; reg298=reg75*reg1042; reg1081=reg468+reg1081;
    reg882=reg882-reg881; reg1040=reg1040-reg1041; reg1223=reg1216+reg1223; reg995=reg996-reg995; reg1198=reg1199+reg1198;
    reg427=reg427-reg1039; reg321=reg869+reg321; reg1166=reg1169+reg1166; reg380=reg654+reg380; reg1082=reg1083+reg1082;
    reg715=reg714+reg715; reg870=reg870-reg871; reg1168=reg1167+reg1168; reg383=reg383+reg992; reg342=reg186+reg342;
    reg717=reg716+reg717; reg1275=reg437+reg1275; reg872=reg1109+reg872; reg720=reg718+reg720; reg873=reg873-reg874;
    reg1144=reg1145+reg1144; reg347=reg405+reg347; reg1184=reg1180+reg1184; reg1047=reg1048+reg1047; reg875=reg876+reg875;
    reg1163=reg538+reg1163; reg490=reg1046+reg490; reg374=reg1071+reg374; reg186=reg770+reg769; reg608=reg644+reg608;
    reg1202=reg1201+reg1202; reg546=reg546-reg753; reg392=reg75*reg247; reg916=reg917+reg916; reg763=reg755-reg763;
    reg1124=reg1124-reg1123; reg918=reg919+reg918; reg542=reg1175+reg542; reg405=reg75*reg543; reg606=reg1196+reg606;
    reg625=reg625-reg920; reg1189=reg1188+reg1189; reg407=reg75*reg536; reg921=reg921-reg922; reg1218=reg1217+reg1218;
    reg728=reg728+reg733; reg412=reg75*reg923; reg422=reg75*reg532; reg1220=reg1222+reg1220; reg424=reg997+reg424;
    reg1191=reg1191-reg1187; reg516=reg516-reg1254; reg559=reg534+reg559; reg520=reg520-reg462; reg246=reg1359+reg246;
    reg1365=reg1365-reg1358; reg174=reg174-reg1370; reg756=reg757-reg756; reg908=reg909+reg908; reg901=reg750+reg901;
    reg1141=reg1141+reg1140; reg910=reg911+reg910; reg745=reg749+reg745; reg1121=reg1122-reg1121; reg610=reg610-reg647;
    reg244=reg1367+reg244; reg498=reg510+reg498; reg782=reg785+reg782; reg913=reg913-reg912; reg556=reg1230+reg556;
    reg780=reg780+reg776; reg914=reg915+reg914; reg253=reg184+reg253; reg184=reg75*reg550; reg1219=reg408+reg1219;
    reg210=reg641+reg210; reg614=reg1183+reg614; reg820=reg810-reg820; reg408=reg75*reg1127; reg547=reg975+reg547;
    reg642=reg642-reg792; reg1177=reg1266+reg1177; reg976=reg977+reg976; reg1179=reg1176+reg1179; reg650=reg798-reg650;
    reg593=reg1182+reg593; reg672=reg676-reg672; reg978=reg979+reg978; reg493=reg409-reg493; reg1171=reg1170+reg1171;
    reg254=reg254-reg817; reg216=reg847+reg216; reg588=reg1174+reg588; reg668=reg667-reg668; reg1128=reg1000+reg1128;
    reg530=reg980+reg530; reg235=reg1277+reg235; reg671=reg670-reg671; reg232=reg1192+reg232; reg684=reg684-reg663;
    reg241=reg622+reg241; reg295=reg295-reg284; reg924=reg925+reg924; reg1214=reg645+reg1214; reg1209=reg1208+reg1209;
    reg777=reg814+reg777; reg1138=reg1139+reg1138; reg771=reg501-reg771; reg926=reg971+reg926; reg207=reg615+reg207;
    reg1158=reg1157+reg1158; reg409=reg75*reg624; reg1125=reg1126-reg1125; reg551=reg972+reg551; reg219=reg219-reg752;
    reg436=reg75*reg1156; reg1152=reg1148+reg1152; reg806=reg768-reg806; reg973=reg974+reg973; reg1149=reg1153+reg1149;
    reg464=reg1137-reg464; reg437=reg75*reg618; reg482=reg179+reg482; reg570=reg312+reg570; reg934=reg75*reg934;
    reg179=reg555+reg935; reg312=reg75*reg726; reg259=reg183+reg259; reg722=reg722-reg723; reg183=reg1373+reg936;
    reg576=reg576+reg1253; reg276=reg276-reg664; reg439=reg1368+reg937; reg317=reg387-reg317; reg189=reg601+reg189;
    reg218=reg1277+reg218; reg661=reg662+reg661; reg387=reg1371+reg938; reg1306=reg1329+reg1306; reg939=reg75*reg939;
    reg655=reg659+reg655; reg367=reg1110-reg367; reg1308=reg1307+reg1308; reg611=reg369+reg611; reg369=reg1283+reg940;
    reg1346=reg391+reg1346; reg391=reg1298+reg927; reg1349=reg1345+reg1349; reg710=reg713+reg710; reg449=reg1297+reg928;
    reg579=reg489+reg579; reg413=reg1113-reg413; reg452=reg1344+reg929; reg1340=reg1339+reg1340; reg704=reg708+reg704;
    reg294=reg388-reg294; reg702=reg703+reg702; reg388=reg1160+reg930; reg453=reg1164+reg931; reg575=reg469+reg575;
    reg574=reg1183+reg574; reg252=reg742+reg252; reg467=reg1285+reg932; reg1111=reg1112+reg1111; reg468=reg1274+reg933;
    reg734=reg736+reg734; reg255=reg572+reg255; reg1114=reg1115+reg1114; reg469=reg1376+reg948; reg1196=reg566+reg1196;
    reg598=reg613+reg598; reg1120=reg1120-reg1119; reg486=reg354+reg949; reg1379=reg1380+reg1379; reg792=reg596-reg792;
    reg204=reg553+reg204; reg489=reg1381+reg950; reg592=reg1235+reg592; reg790=reg791-reg790; reg563=reg200+reg563;
    reg501=reg75*reg1142; reg508=reg75*reg590; reg561=reg540+reg561; reg220=reg356+reg220; reg817=reg587-reg817;
    reg562=reg199+reg562; reg951=reg905+reg951; reg1362=reg1362-reg1364; reg639=reg815-reg639; reg429=reg992+reg429;
    reg906=reg907+reg906; reg199=reg75*reg524; reg1310=reg1311+reg1310; reg652=reg653+reg652; reg941=reg75*reg941;
    reg432=reg1116-reg432; reg692=reg692-reg690; reg356=reg1382+reg942; reg1301=reg1301-reg1332; reg609=reg609-reg331;
    reg510=reg1363+reg943; reg686=reg688+reg686; reg523=reg1360+reg944; reg1331=reg1303+reg1331; reg602=reg1296+reg602;
    reg680=reg685+reg680; reg460=reg394-reg460; reg394=reg1197+reg945; reg1117=reg1118+reg1117; reg605=reg637+reg605;
    reg531=reg1204+reg946; reg1324=reg1334+reg1324; reg1320=reg1320-reg1319; reg802=reg691+reg802; reg534=reg1205+reg947;
    reg315=reg799+reg315; reg1146=reg1150+reg1146; reg466=reg549-reg466; reg1181=reg320+reg1181; reg848=reg849+reg848;
    reg320=reg75*reg457; reg1154=reg1143+reg1154; reg1211=reg1212+reg1211; reg779=reg779+reg778; reg850=reg851+reg850;
    reg385=reg488-reg385; reg488=reg75*reg226; reg535=reg535+reg852; reg537=reg200+reg537; reg494=reg1192+reg494;
    reg200=reg75*reg496; reg1018=reg1018+reg1017; reg853=reg854+reg853; reg774=reg773-reg774; reg479=reg438+reg479;
    reg484=reg262-reg484; reg856=reg856-reg855; reg1194=reg1193+reg1194; reg492=reg594+reg492; reg262=reg75*reg314;
    reg1134=reg1135+reg1134; reg965=reg966+reg965; reg735=reg675-reg735; reg548=reg839+reg548; reg1242=reg1290+reg1242;
    reg308=reg521-reg308; reg1232=reg1231+reg1232; reg1015=reg1015+reg1014; reg438=reg712+reg711; reg967=reg968+reg967;
    reg1299=reg1265+reg1299; reg521=reg75*reg299; reg969=reg970+reg969; reg475=reg475+reg656; reg1375=reg1269+reg1375;
    reg845=reg846+reg845; reg471=reg1175+reg471; reg538=reg75*reg472; reg421=reg1133-reg421; reg540=reg75*reg1016;
    reg818=reg368-reg818; reg847=reg541+reg847; reg1173=reg1279+reg1173; reg759=reg758+reg759; reg864=reg864-reg865;
    reg463=reg217+reg463; reg761=reg760+reg761; reg824=reg824-reg866; reg1342=reg1304+reg1342; reg435=reg762+reg435;
    reg825=reg826+reg825; reg434=reg434-reg353; reg767=reg766+reg767; reg1338=reg1356+reg1338; reg1021=reg1022-reg1021;
    reg827=reg828+reg827; reg288=reg807+reg288; reg1312=reg240+reg1312; reg203=reg203-reg1357; reg245=reg245+reg829;
    reg327=reg327-reg808; reg444=reg187+reg444; reg1351=reg237+reg1351; reg813=reg813-reg812; reg382=reg495-reg382;
    reg830=reg831+reg830; reg187=reg75*reg478; reg205=reg205+reg857; reg497=reg178+reg497; reg784=reg783-reg784;
    reg234=reg1263+reg234; reg858=reg859+reg858; reg410=reg623-reg410; reg1020=reg1020-reg1019; reg171=reg171-reg1200;
    reg860=reg861+reg860; reg629=reg222+reg629; reg402=reg786+reg402; reg1131=reg1132+reg1131; reg211=reg211+reg862;
    reg747=reg746+reg747; reg242=reg185+reg242; reg393=reg454+reg393; reg751=reg748+reg751; reg178=reg75*reg863;
    reg465=reg1327+reg465; reg389=reg754+reg389; reg1318=reg1318-reg343; reg334=reg1052+reg334; reg428=reg1130-reg428;
    reg988=reg952+reg988; reg689=reg689+reg687; reg1240=reg1226+reg1240; reg1280=reg1243+reg1280; reg513=reg862+reg513;
    reg185=reg75*reg677; reg1355=reg1355-reg1282; reg1011=reg1011-reg1010; reg372=reg1235+reg372; reg681=reg681+reg693;
    reg512=reg953+reg512; reg789=reg801-reg789; reg1284=reg1293+reg1284; reg955=reg955-reg954; reg1234=reg1236+reg1234;
    reg217=reg75*reg1229; reg1378=reg366-reg1378; reg222=reg75*reg476; reg509=reg509-reg866; reg816=reg822-reg816;
    reg1224=reg1224-reg348; reg809=reg811-reg809; reg957=reg957-reg956; reg981=reg982+reg981; reg462=reg636-reg462;
    reg1278=reg414+reg1278; reg1268=reg1267+reg1268; reg701=reg700-reg701; reg526=reg852+reg526; reg237=reg75*reg1273;
    reg400=reg1136+reg400; reg240=reg75*reg630; reg525=reg983+reg525; reg1261=reg1261-reg1262; reg349=reg776+reg349;
    reg984=reg984-reg985; reg522=reg1263+reg522; reg332=reg1296+reg332; reg366=reg75*reg333; reg368=reg75*reg403;
    reg727=reg737-reg727; reg517=reg857+reg517; reg1288=reg1287+reg1288; reg753=reg328-reg753; reg986=reg986-reg987;
    reg518=reg1294+reg518; reg694=reg657-reg694; reg503=reg829+reg503; reg328=reg75*reg273; reg504=reg1259+reg504;
    reg1250=reg1249+reg1250; reg781=reg481-reg781; reg558=reg960+reg558; reg269=reg597-reg269; reg580=reg1253+reg580;
    reg961=reg962+reg961; reg1247=reg1248+reg1247; reg414=reg75*reg263; reg1012=reg1013-reg1012; reg552=reg834+reg552;
    reg775=reg266-reg775; reg1258=reg1256+reg1258; reg1252=reg1257+reg1252; reg260=reg589-reg260; reg963=reg964+reg963;
    reg1291=reg487+reg1291; reg313=reg1230+reg313; reg507=reg1227+reg507; reg459=reg1005+reg459; reg752=reg287-reg752;
    reg477=reg480+reg477; reg958=reg959+reg958; reg1255=reg1255-reg1254; reg764=reg765-reg764; reg1079=reg75*reg1079;
    reg1347=reg75*reg1347; reg1086=reg75*reg1086; reg460=reg75*reg460; reg428=reg75*reg428; reg1084=reg75*reg1084;
    reg1219=reg75*reg1219; reg400=reg75*reg400; reg294=reg75*reg294; reg1076=reg75*reg1076; reg266=ponderation*reg222;
    reg287=ponderation*reg392; reg455=reg75*reg455; reg1087=reg75*reg1087; reg1141=reg75*reg1141; reg1134=reg75*reg1134;
    reg492=reg75*reg492; reg493=reg75*reg493; reg385=reg75*reg385; reg382=reg75*reg382; reg1081=reg75*reg1081;
    reg1177=reg75*reg1177; reg1138=reg75*reg1138; reg1077=reg75*reg1077; reg367=reg75*reg367; reg1089=reg75*reg1089;
    reg1286=reg75*reg1286; reg1131=reg75*reg1131; reg204=reg75*reg204; reg454=ponderation*reg501; reg1306=reg75*reg1306;
    reg1338=reg75*reg1338; reg293=reg75*reg293; reg421=reg75*reg421; reg477=reg75*reg477; reg1111=reg75*reg1111;
    reg1163=reg75*reg1163; reg1073=reg75*reg1073; reg464=reg75*reg464; reg500=reg75*reg500; reg1291=reg75*reg1291;
    reg1240=reg75*reg1240; reg1082=reg75*reg1082; reg1190=reg75*reg1190; reg1316=reg75*reg1316; reg1234=reg75*reg1234;
    reg955=reg75*reg955; reg512=reg75*reg512; reg1355=reg75*reg1355; reg513=reg75*reg513; reg988=reg75*reg988;
    reg986=reg75*reg986; reg517=reg75*reg517; reg522=reg75*reg522; reg984=reg75*reg984; reg525=reg75*reg525;
    reg481=ponderation*reg237; reg526=reg75*reg526; reg981=reg75*reg981; reg530=reg75*reg530; reg588=reg75*reg588;
    reg216=reg75*reg216; reg978=reg75*reg978; reg593=reg75*reg593; reg976=reg75*reg976; reg547=reg75*reg547;
    reg1149=reg75*reg1149; reg210=reg75*reg210; reg973=reg75*reg973; reg551=reg75*reg551; reg504=reg75*reg504;
    reg856=reg75*reg856; reg853=reg75*reg853; reg537=reg75*reg537; reg535=reg75*reg535; reg1154=reg75*reg1154;
    reg850=reg75*reg850; reg848=reg75*reg848; reg1173=reg75*reg1173; reg847=reg75*reg847; reg1375=reg75*reg1375;
    reg845=reg75*reg845; reg969=reg75*reg969; reg967=reg75*reg967; reg1232=reg75*reg1232; reg548=reg75*reg548;
    reg965=reg75*reg965; reg1258=reg75*reg1258; reg963=reg75*reg963; reg552=reg75*reg552; reg580=reg75*reg580;
    reg961=reg75*reg961; reg558=reg75*reg558; reg503=reg75*reg503; reg958=reg75*reg958; reg507=reg75*reg507;
    reg957=reg75*reg957; reg509=reg75*reg509; reg1196=reg75*reg1196; reg487=reg75*reg469; reg495=reg75*reg534;
    reg541=reg75*reg531; reg1331=reg75*reg1331; reg549=reg75*reg394; reg553=reg75*reg523; reg557=reg75*reg510;
    reg564=reg75*reg356; reg1310=reg75*reg1310; reg941=ponderation*reg941; reg566=reg75*reg369; reg939=ponderation*reg939;
    reg567=reg75*reg387; reg576=reg75*reg576; reg572=reg75*reg439; reg577=reg75*reg183; reg587=reg75*reg179;
    reg255=reg75*reg255; reg934=ponderation*reg934; reg589=reg75*reg468; reg591=reg75*reg467; reg594=reg75*reg453;
    reg596=reg75*reg388; reg597=reg75*reg452; reg601=reg75*reg449; reg607=reg75*reg391; reg207=reg75*reg207;
    reg926=reg75*reg926; reg1214=reg75*reg1214; reg924=reg75*reg924; reg241=reg75*reg241; reg1220=reg75*reg1220;
    reg612=ponderation*reg412; reg921=reg75*reg921; reg606=reg75*reg606; reg625=reg75*reg625; reg918=reg75*reg918;
    reg916=reg75*reg916; reg608=reg75*reg608; reg914=reg75*reg914; reg913=reg75*reg913; reg610=reg75*reg610;
    reg910=reg75*reg910; reg908=reg75*reg908; reg246=reg75*reg246; reg559=reg75*reg559; reg906=reg75*reg906;
    reg562=reg75*reg562; reg951=reg75*reg951; reg561=reg75*reg561; reg563=reg75*reg563; reg613=reg75*reg489;
    reg615=reg75*reg486; reg352=reg75*reg352; reg887=reg75*reg887; reg885=reg75*reg885; reg1206=reg75*reg1206;
    reg359=reg75*reg359; reg883=reg75*reg883; reg1223=reg75*reg1223; reg882=reg75*reg882; reg364=reg75*reg364;
    reg370=reg75*reg370; reg880=reg75*reg880; reg877=reg75*reg877; reg1155=reg75*reg1155; reg374=reg75*reg374;
    reg875=reg75*reg875; reg1184=reg75*reg1184; reg873=reg75*reg873; reg872=reg75*reg872; reg1168=reg75*reg1168;
    reg870=reg75*reg870; reg321=reg75*reg321; reg324=reg75*reg324; reg322=reg75*reg322; reg867=reg75*reg867;
    reg1295=reg75*reg1295; reg329=reg75*reg329; reg622=ponderation*reg264; reg306=reg75*reg306; reg304=reg75*reg304;
    reg1069=reg75*reg1069; reg177=reg75*reg177; reg1068=reg75*reg1068; reg1354=reg75*reg1354; reg311=reg75*reg311;
    reg1108=reg75*reg1108; reg1105=reg75*reg1105; reg617=reg75*reg617; reg318=reg75*reg318; reg1302=reg75*reg1302;
    reg1102=reg75*reg1102; reg1100=reg75*reg1100; reg623=ponderation*reg303; reg265=reg75*reg265; reg1098=reg75*reg1098;
    reg270=reg75*reg270; reg1096=reg75*reg1096; reg274=reg75*reg274; reg279=reg75*reg279; reg275=reg75*reg275;
    reg1093=reg75*reg1093; reg201=reg75*reg201; reg285=reg75*reg285; reg286=reg75*reg286; reg224=reg75*reg224;
    reg209=reg75*reg209; reg250=reg75*reg250; reg837=reg75*reg837; reg1353=reg75*reg1353; reg835=reg75*reg835;
    reg223=reg75*reg223; reg620=reg75*reg620; reg832=reg75*reg832; reg830=reg75*reg830; reg444=reg75*reg444;
    reg245=reg75*reg245; reg827=reg75*reg827; reg518=reg75*reg518; reg825=reg75*reg825; reg1342=reg75*reg1342;
    reg824=reg75*reg824; reg864=reg75*reg864; reg1318=reg75*reg1318; reg636=ponderation*reg178; reg242=reg75*reg242;
    reg211=reg75*reg211; reg860=reg75*reg860; reg234=reg75*reg234; reg858=reg75*reg858; reg205=reg75*reg205;
    reg1194=reg75*reg1194; reg330=reg75*reg330; reg1239=reg75*reg1239; reg902=reg75*reg902; reg782=reg75*reg782;
    reg899=reg75*reg899; reg346=reg75*reg346; reg336=reg75*reg336; reg897=reg75*reg897; reg1251=reg75*reg1251;
    reg357=reg75*reg357; reg350=reg75*reg350; reg1260=reg75*reg1260; reg637=ponderation*reg229; reg1241=reg75*reg1241;
    reg373=reg75*reg373; reg236=reg75*reg236; reg1178=reg75*reg1178; reg892=reg75*reg892; reg891=reg75*reg891;
    reg190=reg75*reg190; reg191=reg75*reg191; reg889=reg75*reg889; reg844=reg75*reg844; reg1203=reg75*reg1203;
    reg227=reg75*reg227; reg842=reg75*reg842; reg840=reg75*reg840; reg351=reg75*reg351; reg699=reg75*reg699;
    reg1179=reg75*reg1179; reg1264=reg75*reg1264; reg365=reg75*reg365; reg820=reg75*reg820; reg706=reg75*reg706;
    reg709=reg75*reg709; reg641=ponderation*reg437; reg614=reg75*reg614; reg1166=reg75*reg1166; reg380=reg75*reg380;
    reg806=reg75*reg806; reg715=reg75*reg715; reg342=reg75*reg342; reg717=reg75*reg717; reg219=reg75*reg219;
    reg1152=reg75*reg1152; reg720=reg75*reg720; reg1144=reg75*reg1144; reg347=reg75*reg347; reg644=ponderation*reg409;
    reg1047=reg75*reg1047; reg1158=reg75*reg1158; reg1161=reg75*reg1161; reg490=reg75*reg490; reg1244=reg75*reg1244;
    reg458=reg75*reg458; reg721=reg75*reg721; reg701=reg75*reg701; reg1268=reg75*reg1268; reg724=reg75*reg724;
    reg462=reg75*reg462; reg305=reg75*reg305; reg290=reg75*reg290; reg671=reg75*reg671; reg375=reg75*reg375;
    reg729=reg75*reg729; reg235=reg75*reg235; reg732=reg75*reg732; reg668=reg75*reg668; reg1233=reg75*reg1233;
    reg272=reg75*reg272; reg254=reg75*reg254; reg739=reg75*reg739; reg1171=reg75*reg1171; reg741=reg75*reg741;
    reg672=reg75*reg672; reg1281=reg75*reg1281; reg743=reg75*reg743; reg650=reg75*reg650; reg697=reg75*reg697;
    reg642=reg75*reg642; reg1033=reg75*reg1033; reg645=ponderation*reg405; reg283=reg75*reg283; reg763=reg75*reg763;
    reg1366=reg75*reg1366; reg1067=reg75*reg1067; reg542=reg75*reg542; reg289=reg75*reg289; reg546=reg75*reg546;
    reg1202=reg75*reg1202; reg631=reg75*reg631; reg278=reg75*reg278; reg186=reg75*reg186; reg1064=reg75*reg1064;
    reg271=reg75*reg271; reg378=reg75*reg378; reg649=ponderation*reg184; reg253=reg75*reg253; reg277=reg75*reg277;
    reg780=reg75*reg780; reg1323=reg75*reg1323; reg556=reg75*reg556; reg1062=reg75*reg1062; reg1060=reg75*reg1060;
    reg498=reg75*reg498; reg745=reg75*reg745; reg1333=reg75*reg1333; reg341=reg75*reg341; reg771=reg75*reg771;
    reg1044=reg75*reg1044; reg777=reg75*reg777; reg1209=reg75*reg1209; reg1215=reg75*reg1215; reg451=reg75*reg451;
    reg295=reg75*reg295; reg339=reg75*reg339; reg684=reg75*reg684; reg1221=reg75*reg1221; reg653=ponderation*reg298;
    reg232=reg75*reg232; reg1040=reg75*reg1040; reg654=ponderation*reg422; reg1198=reg75*reg1198; reg427=reg75*reg427;
    reg1037=reg75*reg1037; reg728=reg75*reg728; reg1218=reg75*reg1218; reg355=reg75*reg355; reg307=reg75*reg307;
    reg657=ponderation*reg407; reg358=reg75*reg358; reg1189=reg75*reg1189; reg326=reg75*reg326; reg1034=reg75*reg1034;
    reg402=reg75*reg402; reg260=reg75*reg260; reg313=reg75*reg313; reg393=reg75*reg393; reg747=reg75*reg747;
    reg751=reg75*reg751; reg775=reg75*reg775; reg658=ponderation*reg414; reg465=reg75*reg465; reg389=reg75*reg389;
    reg1247=reg75*reg1247; reg759=reg75*reg759; reg463=reg75*reg463; reg761=reg75*reg761; reg269=reg75*reg269;
    reg781=reg75*reg781; reg434=reg75*reg434; reg1250=reg75*reg1250; reg435=reg75*reg435; reg767=reg75*reg767;
    reg659=ponderation*reg328; reg203=reg75*reg203; reg288=reg75*reg288; reg764=reg75*reg764; reg1255=reg75*reg1255;
    reg327=reg75*reg327; reg818=reg75*reg818; reg1146=reg75*reg1146; reg466=reg75*reg466; reg662=ponderation*reg538;
    reg471=reg75*reg471; reg665=ponderation*reg320; reg1211=reg75*reg1211; reg779=reg75*reg779; reg475=reg75*reg475;
    reg667=ponderation*reg521; reg669=ponderation*reg488; reg494=reg75*reg494; reg670=ponderation*reg200; reg1299=reg75*reg1299;
    reg673=reg75*reg438; reg774=reg75*reg774; reg479=reg75*reg479; reg484=reg75*reg484; reg308=reg75*reg308;
    reg497=reg75*reg497; reg675=ponderation*reg187; reg735=reg75*reg735; reg1242=reg75*reg1242; reg784=reg75*reg784;
    reg676=ponderation*reg262; reg410=reg75*reg410; reg629=reg75*reg629; reg805=reg75*reg805; reg689=reg75*reg689;
    reg1280=reg75*reg1280; reg384=reg75*reg384; reg442=reg75*reg442; reg674=reg75*reg674; reg694=reg75*reg694;
    reg1159=reg75*reg1159; reg753=reg75*reg753; reg679=reg75*reg679; reg1288=reg75*reg1288; reg683=reg75*reg683;
    reg1172=reg75*reg1172; reg474=reg75*reg474; reg727=reg75*reg727; reg696=reg75*reg696; reg678=ponderation*reg366;
    reg332=reg75*reg332; reg651=reg75*reg651; reg443=reg75*reg443; reg491=reg75*reg491; reg349=reg75*reg349;
    reg682=ponderation*reg240; reg360=reg75*reg360; reg660=reg75*reg660; reg1261=reg75*reg1261; reg666=reg75*reg666;
    reg1351=reg75*reg1351; reg752=reg75*reg752; reg813=reg75*reg813; reg685=ponderation*reg316; reg809=reg75*reg809;
    reg292=reg75*reg292; reg1224=reg75*reg1224; reg450=reg75*reg450; reg816=reg75*reg816; reg310=reg75*reg310;
    reg787=reg75*reg787; reg291=reg75*reg291; reg1330=reg75*reg1330; reg1378=reg75*reg1378; reg688=ponderation*reg217;
    reg340=reg75*reg340; reg793=reg75*reg793; reg789=reg75*reg789; reg430=reg75*reg430; reg796=reg75*reg796;
    reg681=reg75*reg681; reg372=reg75*reg372; reg238=reg75*reg238; reg423=reg75*reg423; reg803=reg75*reg803;
    reg691=ponderation*reg185; reg695=ponderation*reg202; reg1125=reg75*reg1125; reg698=ponderation*reg408; reg188=reg75*reg188;
    reg1379=reg75*reg1379; reg386=reg75*reg386; reg700=ponderation*reg436; reg315=reg75*reg315; reg570=reg75*reg570;
    reg1210=reg75*reg1210; reg1128=reg75*reg1128; reg995=reg75*reg995; reg993=reg75*reg993; reg802=reg75*reg802;
    reg1278=reg75*reg1278; reg1320=reg75*reg1320; reg1275=reg75*reg1275; reg383=reg75*reg383; reg703=ponderation*reg368;
    reg259=reg75*reg259; reg605=reg75*reg605; reg220=reg75*reg220; reg1003=reg75*reg1003; reg575=reg75*reg575;
    reg483=reg75*reg483; reg1001=reg75*reg1001; reg1124=reg75*reg1124; reg790=reg75*reg790; reg195=reg75*reg195;
    reg252=reg75*reg252; reg592=reg75*reg592; reg792=reg75*reg792; reg424=reg75*reg424; reg485=reg75*reg485;
    reg1191=reg75*reg1191; reg705=ponderation*reg390; reg707=ponderation*reg261; reg482=reg75*reg482; reg598=reg75*reg598;
    reg734=reg75*reg734; reg998=reg75*reg998; reg1015=reg75*reg1015; reg692=reg75*reg692; reg456=reg75*reg456;
    reg431=reg75*reg431; reg1301=reg75*reg1301; reg218=reg75*reg218; reg661=reg75*reg661; reg708=ponderation*reg540;
    reg1181=reg75*reg1181; reg1024=reg75*reg1024; reg1312=reg75*reg1312; reg1021=reg75*reg1021; reg652=reg75*reg652;
    reg334=reg75*reg334; reg611=reg75*reg611; reg1308=reg75*reg1308; reg1018=reg75*reg1018; reg655=reg75*reg655;
    reg171=reg75*reg171; reg1020=reg75*reg1020; reg990=reg75*reg990; reg713=ponderation*reg312; reg1011=reg75*reg1011;
    reg1225=reg75*reg1225; reg1284=reg75*reg1284; reg1031=reg75*reg1031; reg680=reg75*reg680; reg722=reg75*reg722;
    reg1029=reg75*reg1029; reg459=reg75*reg459; reg602=reg75*reg602; reg686=reg75*reg686; reg181=reg75*reg181;
    reg425=reg75*reg425; reg276=reg75*reg276; reg1012=reg75*reg1012; reg1028=reg75*reg1028; reg609=reg75*reg609;
    reg1252=reg75*reg1252; reg1025=reg75*reg1025; reg520=reg75*reg520; reg579=reg75*reg579; reg432=reg75*reg432;
    reg225=reg75*reg225; reg1050=reg75*reg1050; reg516=reg75*reg516; reg401=reg75*reg401; reg714=ponderation*reg199;
    reg1340=reg75*reg1340; reg1117=reg75*reg1117; reg446=reg75*reg446; reg716=ponderation*reg180; reg704=reg75*reg704;
    reg1120=reg75*reg1120; reg1324=reg75*reg1324; reg1335=reg75*reg1335; reg1008=reg75*reg1008; reg1058=reg75*reg1058;
    reg244=reg75*reg244; reg413=reg75*reg413; reg1057=reg75*reg1057; reg1346=reg75*reg1346; reg1349=reg75*reg1349;
    reg470=reg75*reg470; reg1055=reg75*reg1055; reg901=reg75*reg901; reg710=reg75*reg710; reg1114=reg75*reg1114;
    reg1053=reg75*reg1053; reg756=reg75*reg756; reg560=reg75*reg560; reg174=reg75*reg174; reg461=reg75*reg461;
    reg317=reg75*reg317; reg189=reg75*reg189; reg429=reg75*reg429; reg817=reg75*reg817; reg1343=reg75*reg1343;
    reg499=reg75*reg499; reg718=ponderation*reg228; reg1121=reg75*reg1121; reg719=ponderation*reg508; reg1365=reg75*reg1365;
    reg702=reg75*reg702; reg574=reg75*reg574; reg639=reg75*reg639; reg1362=reg75*reg1362; T tmp_10_17=ponderation*reg681;
    T tmp_21_22=ponderation*reg845; reg681=ponderation*reg567; sollicitation[indices[4]+0]+=reg681; reg731=ponderation*reg601; sollicitation[indices[7]+1]+=reg731;
    T tmp_7_22=ponderation*reg574; T tmp_6_13=ponderation*reg1280; T tmp_8_1=ponderation*reg710; T tmp_8_11=ponderation*reg276; T tmp_22_11=ponderation*reg955;
    T tmp_10_18=ponderation*reg789; T tmp_6_3=ponderation*reg471; T tmp_11_10=-reg667; T tmp_6_11=-reg688; T tmp_8_6=ponderation*reg252;
    T tmp_7_19=ponderation*reg218; reg218=ponderation*reg572; sollicitation[indices[4]+1]+=reg218; T tmp_22_10=ponderation*reg509; T tmp_7_23=ponderation*reg1340;
    T tmp_22_4=ponderation*reg552; T tmp_6_2=ponderation*reg1146; T tmp_22_3=ponderation*reg963; T tmp_21_21=ponderation*reg847; T tmp_6_6=ponderation*reg313;
    T tmp_7_18=ponderation*reg1308; T tmp_22_13=ponderation*reg513; T tmp_10_15=ponderation*reg689; T tmp_8_12=ponderation*reg661; T tmp_6_12=ponderation*reg372;
    sollicitation[indices[3]+2]+=-reg939; T tmp_8_13=ponderation*reg655; T tmp_11_12=-reg662; reg252=ponderation*reg607; sollicitation[indices[7]+2]+=reg252;
    T tmp_10_16=-reg691; T tmp_1_4=ponderation*reg576; T tmp_0_16=ponderation*reg1355; T tmp_11_4=ponderation*reg775; reg276=ponderation*reg591;
    sollicitation[indices[6]+0]+=reg276; T tmp_8_0=ponderation*reg1349; T tmp_22_12=ponderation*reg512; T tmp_8_5=ponderation*reg575; T tmp_11_3=-reg658;
    T tmp_11_11=ponderation*reg475; T tmp_0_5=ponderation*reg1232; T tmp_7_20=ponderation*reg259; reg259=ponderation*reg587; sollicitation[indices[5]+0]+=reg259;
    T tmp_10_22=ponderation*reg752; T tmp_0_18=ponderation*reg507; T tmp_22_5=ponderation*reg961; reg313=ponderation*reg589; sollicitation[indices[5]+2]+=reg313;
    T tmp_7_21=ponderation*reg482; T tmp_11_8=ponderation*reg308; T tmp_1_3=ponderation*reg255; T tmp_22_7=ponderation*reg503; T tmp_11_6=-reg676;
    T tmp_6_5=ponderation*reg1242; T tmp_11_7=ponderation*reg735; T tmp_8_8=ponderation*reg570; reg255=ponderation*reg596; sollicitation[indices[6]+2]+=reg255;
    T tmp_8_3=ponderation*reg704; T tmp_10_23=ponderation*reg764; T tmp_6_8=ponderation*reg1250; T tmp_11_1=ponderation*reg781; T tmp_11_0=-reg659;
    T tmp_8_7=ponderation*reg734; T tmp_22_1=ponderation*reg548; T tmp_22_6=ponderation*reg558; sollicitation[indices[5]+1]+=-reg934; T tmp_10_19=ponderation*reg1378;
    T tmp_11_5=ponderation*reg260; T tmp_6_10=ponderation*reg1224; T tmp_8_10=ponderation*reg722; T tmp_0_17=ponderation*reg1234; T tmp_21_23=ponderation*reg969;
    T tmp_0_3=ponderation*reg580; T tmp_0_4=ponderation*reg1258; T tmp_22_2=ponderation*reg965; T tmp_8_9=-reg713; T tmp_11_9=ponderation*reg673;
    reg260=ponderation*reg594; sollicitation[indices[6]+1]+=reg260; T tmp_10_20=ponderation*reg816; T tmp_6_7=ponderation*reg1247; T tmp_22_9=ponderation*reg957;
    T tmp_11_2=ponderation*reg269; reg269=ponderation*reg577; sollicitation[indices[4]+2]+=reg269; reg308=ponderation*reg597; sollicitation[indices[7]+0]+=reg308;
    T tmp_8_4=ponderation*reg702; T tmp_6_4=ponderation*reg1299; T tmp_10_21=ponderation*reg809; T tmp_8_2=ponderation*reg579; T tmp_6_9=ponderation*reg1255;
    T tmp_22_8=ponderation*reg958; T tmp_22_0=ponderation*reg967; T tmp_23_23=ponderation*reg561; T tmp_23_7=ponderation*reg924; T tmp_7_0=ponderation*reg1209;
    T tmp_20_10=ponderation*reg684; T tmp_9_2=-reg719; T tmp_7_11=ponderation*reg1362; T tmp_23_8=ponderation*reg241; T tmp_0_7=ponderation*reg1220;
    T tmp_20_9=-reg637; T tmp_23_22=ponderation*reg951; T tmp_9_3=ponderation*reg817; T tmp_7_1=ponderation*reg232; T tmp_9_16=ponderation*reg728;
    T tmp_23_9=-reg612; T tmp_0_1=ponderation*reg562; T tmp_9_15=-reg657; T tmp_23_21=ponderation*reg906; T tmp_8_23=ponderation*reg598;
    T tmp_7_13=ponderation*reg592; T tmp_0_19=ponderation*reg504; T tmp_6_22=ponderation*reg1152; T tmp_23_4=ponderation*reg551; reg232=ponderation*reg613;
    sollicitation[indices[0]+0]+=reg232; T tmp_9_20=-reg644; T tmp_9_0=ponderation*reg792; T tmp_6_23=ponderation*reg1158; T tmp_23_5=ponderation*reg207;
    T tmp_9_19=ponderation*reg771; T tmp_0_8=ponderation*reg1214; T tmp_9_1=ponderation*reg790; T tmp_7_12=ponderation*reg220; T tmp_20_12=ponderation*reg777;
    T tmp_23_6=ponderation*reg926; T tmp_0_0=ponderation*reg563; T tmp_20_11=ponderation*reg295; T tmp_9_11=ponderation*reg186; T tmp_7_5=ponderation*reg1202;
    T tmp_23_14=ponderation*reg608; T tmp_9_10=-reg649; T tmp_0_2=ponderation*reg246; T tmp_9_7=ponderation*reg756; T tmp_23_18=ponderation*reg910;
    T tmp_9_9=ponderation*reg780; T tmp_23_15=ponderation*reg914; T tmp_7_6=ponderation*reg253; T tmp_7_7=ponderation*reg556; T tmp_20_0=ponderation*reg901;
    T tmp_23_16=ponderation*reg913; T tmp_20_2=ponderation*reg498; T tmp_7_8=ponderation*reg244; T tmp_20_1=ponderation*reg745; T tmp_23_17=ponderation*reg610;
    T tmp_9_4=ponderation*reg639; T tmp_23_10=ponderation*reg921; T tmp_0_6=ponderation*reg606; T tmp_7_2=ponderation*reg1218; T tmp_7_10=ponderation*reg516;
    T tmp_9_14=-reg645; T tmp_23_11=ponderation*reg625; T tmp_9_5=-reg714; T tmp_7_3=ponderation*reg1189; T tmp_9_13=ponderation*reg763;
    T tmp_23_20=ponderation*reg559; T tmp_23_12=ponderation*reg918; T tmp_7_4=ponderation*reg542; T tmp_7_9=ponderation*reg174; T tmp_9_12=ponderation*reg546;
    T tmp_23_13=ponderation*reg916; T tmp_9_6=ponderation*reg520; T tmp_23_19=ponderation*reg908; T tmp_22_17=ponderation*reg984; T tmp_10_9=-reg682;
    T tmp_8_16=ponderation*reg692; T tmp_22_18=ponderation*reg525; reg174=ponderation*reg557; sollicitation[indices[2]+1]+=reg174; T tmp_0_11=-reg481;
    T tmp_6_16=ponderation*reg1261; T tmp_10_8=ponderation*reg701; T tmp_8_17=ponderation*reg609; T tmp_7_16=ponderation*reg602; reg186=ponderation*reg553;
    sollicitation[indices[2]+0]+=reg186; T tmp_22_19=ponderation*reg526; T tmp_10_7=ponderation*reg462; T tmp_6_17=ponderation*reg1268; T tmp_8_18=ponderation*reg686;
    T tmp_10_6=ponderation*reg671; T tmp_22_20=ponderation*reg981; T tmp_10_14=ponderation*reg694; T tmp_22_14=ponderation*reg988; reg207=ponderation*reg566;
    sollicitation[indices[3]+1]+=reg207; T tmp_10_13=ponderation*reg753; T tmp_8_14=ponderation*reg611; T tmp_22_15=ponderation*reg986; sollicitation[indices[3]+0]+=-reg941;
    T tmp_6_14=ponderation*reg1288; T tmp_10_12=ponderation*reg727; T tmp_0_12=ponderation*reg522; T tmp_10_11=-reg678; T tmp_22_16=ponderation*reg517;
    T tmp_8_15=ponderation*reg652; T tmp_7_17=ponderation*reg1301; T tmp_1_5=ponderation*reg1310; T tmp_6_15=ponderation*reg332; T tmp_10_10=ponderation*reg349;
    reg220=ponderation*reg564; sollicitation[indices[2]+2]+=reg220; T tmp_10_1=ponderation*reg642; T tmp_8_21=ponderation*reg802; T tmp_23_0=ponderation*reg976;
    T tmp_6_20=ponderation*reg1179; reg241=ponderation*reg487; sollicitation[indices[0]+2]+=reg241; T tmp_10_0=ponderation*reg820; T tmp_23_1=ponderation*reg547;
    T tmp_8_22=ponderation*reg315; T tmp_0_20=ponderation*reg1149; T tmp_9_23=-reg641; T tmp_7_14=ponderation*reg1379; T tmp_23_2=ponderation*reg210;
    T tmp_1_7=ponderation*reg1196; T tmp_6_21=ponderation*reg614; T tmp_9_22=ponderation*reg806; reg210=ponderation*reg615; sollicitation[indices[0]+1]+=reg210;
    T tmp_9_21=ponderation*reg219; T tmp_23_3=ponderation*reg973; reg219=ponderation*reg549; sollicitation[indices[1]+2]+=reg219; T tmp_0_22=ponderation*reg588;
    T tmp_10_5=ponderation*reg668; T tmp_22_21=ponderation*reg530; T tmp_8_19=ponderation*reg680; T tmp_6_18=ponderation*reg235; T tmp_7_15=ponderation*reg1320;
    T tmp_10_4=ponderation*reg254; reg235=ponderation*reg541; sollicitation[indices[1]+1]+=reg235; T tmp_22_22=ponderation*reg216; T tmp_10_3=ponderation*reg672;
    T tmp_0_21=ponderation*reg593; T tmp_8_20=ponderation*reg605; T tmp_1_6=ponderation*reg1331; T tmp_6_19=ponderation*reg1171; T tmp_22_23=ponderation*reg978;
    T tmp_10_2=ponderation*reg650; reg216=ponderation*reg495; sollicitation[indices[1]+0]+=reg216; T tmp_4_3=ponderation*reg1335; T tmp_18_14=ponderation*reg1068;
    T tmp_15_9=-reg716; T tmp_2_12=ponderation*reg1354; T tmp_4_4=ponderation*reg446; T tmp_15_8=ponderation*reg401; T tmp_18_15=ponderation*reg311;
    T tmp_15_7=ponderation*reg1050; T tmp_18_16=ponderation*reg1108; T tmp_4_5=ponderation*reg225; T tmp_2_11=ponderation*reg617; T tmp_15_6=ponderation*reg461;
    T tmp_18_17=ponderation*reg1105; T tmp_15_5=ponderation*reg1053; T tmp_2_10=ponderation*reg1302; T tmp_4_6=ponderation*reg560; T tmp_15_4=ponderation*reg1055;
    T tmp_18_18=ponderation*reg318; T tmp_15_3=ponderation*reg1057; T tmp_18_19=ponderation*reg1102; T tmp_4_7=ponderation*reg470; T tmp_15_2=ponderation*reg1058;
    T tmp_15_17=ponderation*reg998; T tmp_3_23=ponderation*reg188; T tmp_18_8=ponderation*reg1077; T tmp_2_16=ponderation*reg1316; T tmp_15_16=-reg707;
    T tmp_15_15=ponderation*reg485; T tmp_18_9=ponderation*reg1076; T tmp_9_8=-reg705; T tmp_2_15=ponderation*reg1347; T tmp_4_0=ponderation*reg195;
    T tmp_18_10=ponderation*reg1073; T tmp_15_14=ponderation*reg1001; T tmp_15_13=ponderation*reg1003; T tmp_18_11=-reg622; T tmp_4_1=ponderation*reg483;
    T tmp_2_14=ponderation*reg306; T tmp_15_12=ponderation*reg499; T tmp_4_2=ponderation*reg1343; T tmp_18_12=ponderation*reg304; T tmp_15_11=-reg718;
    T tmp_2_13=ponderation*reg177; T tmp_15_10=ponderation*reg1008; T tmp_18_13=ponderation*reg1069; T tmp_14_18=ponderation*reg1067; T tmp_14_17=ponderation*reg283;
    T tmp_19_3=ponderation*reg285; T tmp_4_11=ponderation*reg1366; T tmp_2_5=ponderation*reg352; T tmp_14_16=ponderation*reg1033; T tmp_19_4=ponderation*reg286;
    T tmp_14_15=ponderation*reg1034; T tmp_4_12=ponderation*reg326; T tmp_19_5=ponderation*reg887; T tmp_14_14=ponderation*reg358; T tmp_2_4=ponderation*reg1206;
    T tmp_14_13=ponderation*reg307; T tmp_19_6=ponderation*reg885; T tmp_4_13=ponderation*reg355; T tmp_14_12=ponderation*reg1037; T tmp_19_7=ponderation*reg359;
    T tmp_14_11=ponderation*reg427; T tmp_2_3=ponderation*reg1223; T tmp_19_8=ponderation*reg883; T tmp_14_10=ponderation*reg1040; T tmp_4_14=ponderation*reg1198;
    T tmp_14_9=-reg653; T tmp_18_20=ponderation*reg1100; T tmp_2_9=-reg623; T tmp_4_8=ponderation*reg1333; T tmp_15_1=ponderation*reg1060;
    T tmp_15_0=ponderation*reg1062; T tmp_18_21=ponderation*reg265; T tmp_2_8=ponderation*reg270; T tmp_4_9=ponderation*reg1323; T tmp_18_22=ponderation*reg1098;
    T tmp_9_17=-reg654; T tmp_14_23=ponderation*reg277; T tmp_14_22=ponderation*reg378; T tmp_18_23=ponderation*reg1096; T tmp_2_7=ponderation*reg279;
    T tmp_14_21=ponderation*reg1064; T tmp_19_0=ponderation*reg274; T tmp_9_18=ponderation*reg271; T tmp_14_20=ponderation*reg278; T tmp_19_1=ponderation*reg275;
    T tmp_4_10=ponderation*reg631; T tmp_2_6=ponderation*reg201; T tmp_14_19=ponderation*reg289; T tmp_19_2=ponderation*reg1093; T tmp_17_11=-reg287;
    T tmp_3_10=ponderation*reg1191; T tmp_16_18=ponderation*reg1125; T tmp_17_12=ponderation*reg1138; T tmp_3_2=ponderation*reg1177; T tmp_16_17=-reg698;
    T tmp_17_13=ponderation*reg464; T tmp_3_11=-reg700; T tmp_16_16=ponderation*reg1128; T tmp_17_14=ponderation*reg493; T tmp_3_1=ponderation*reg1240;
    T tmp_16_15=-reg703; T tmp_17_15=ponderation*reg400; T tmp_3_12=ponderation*reg1278; T tmp_16_14=ponderation*reg1011; T tmp_16_13=ponderation*reg459;
    T tmp_17_16=-reg266; T tmp_3_0=ponderation*reg1291; T tmp_3_13=ponderation*reg1284; T tmp_16_12=ponderation*reg1012; T tmp_17_17=ponderation*reg477;
    T tmp_16_11=ponderation*reg1015; T tmp_17_18=ponderation*reg1134; T tmp_17_4=ponderation*reg413; T tmp_3_6=ponderation*reg1346; T tmp_17_5=ponderation*reg294;
    T tmp_17_3=ponderation*reg1114; T tmp_3_5=ponderation*reg1306; T tmp_17_2=ponderation*reg317; T tmp_17_6=ponderation*reg1111; T tmp_17_1=ponderation*reg432;
    T tmp_17_7=ponderation*reg367; T tmp_3_7=ponderation*reg189; T tmp_17_0=ponderation*reg1117; T tmp_16_23=ponderation*reg1120; T tmp_17_8=ponderation*reg460;
    T tmp_3_4=ponderation*reg204; T tmp_3_8=ponderation*reg1324; T tmp_16_22=ponderation*reg429; T tmp_17_9=-reg454; T tmp_16_21=ponderation*reg1121;
    T tmp_17_10=ponderation*reg1141; T tmp_3_3=ponderation*reg1219; T tmp_3_9=ponderation*reg1365; T tmp_16_20=ponderation*reg1124; T tmp_16_19=ponderation*reg424;
    T tmp_18_1=ponderation*reg1089; T tmp_2_20=ponderation*reg455; T tmp_16_1=ponderation*reg425; T tmp_18_2=ponderation*reg1087; T tmp_3_19=ponderation*reg181;
    T tmp_16_0=ponderation*reg1029; T tmp_2_19=ponderation*reg1286; T tmp_15_23=ponderation*reg1031; T tmp_18_3=ponderation*reg1086; T tmp_3_20=ponderation*reg1225;
    T tmp_15_22=ponderation*reg990; T tmp_18_4=ponderation*reg1084; T tmp_15_21=ponderation*reg383; T tmp_2_18=ponderation*reg1163; T tmp_18_5=ponderation*reg1082;
    T tmp_3_21=ponderation*reg1275; T tmp_15_20=ponderation*reg993; T tmp_15_19=ponderation*reg995; T tmp_18_6=ponderation*reg1081; T tmp_2_17=ponderation*reg293;
    T tmp_3_22=ponderation*reg1210; T tmp_15_18=ponderation*reg386; T tmp_18_7=ponderation*reg1079; T tmp_3_14=ponderation*reg1252; T tmp_16_10=-reg708;
    T tmp_17_19=ponderation*reg421; T tmp_2_23=ponderation*reg492; T tmp_3_15=ponderation*reg1181; T tmp_16_9=ponderation*reg1018; T tmp_17_20=ponderation*reg385;
    T tmp_16_8=ponderation*reg1020; T tmp_16_7=ponderation*reg334; T tmp_17_21=ponderation*reg1131; T tmp_2_22=ponderation*reg1338; T tmp_3_16=ponderation*reg171;
    T tmp_16_6=ponderation*reg1021; T tmp_17_22=ponderation*reg428; T tmp_16_5=ponderation*reg1024; T tmp_17_23=ponderation*reg382; T tmp_3_17=ponderation*reg1312;
    T tmp_2_21=ponderation*reg1190; T tmp_16_4=ponderation*reg431; T tmp_16_3=ponderation*reg1025; T tmp_18_0=ponderation*reg500; T tmp_3_18=ponderation*reg456;
    T tmp_16_2=ponderation*reg1028; T tmp_21_2=ponderation*reg835; T tmp_12_13=ponderation*reg787; T tmp_12_12=ponderation*reg450; T tmp_21_3=ponderation*reg223;
    T tmp_0_15=ponderation*reg620; T tmp_5_13=ponderation*reg310; T tmp_12_11=-reg685; T tmp_21_4=ponderation*reg832; T tmp_5_14=ponderation*reg292;
    T tmp_12_10=ponderation*reg813; T tmp_21_5=ponderation*reg830; T tmp_0_14=ponderation*reg444; T tmp_12_9=ponderation*reg327; T tmp_5_15=ponderation*reg1351;
    T tmp_12_8=ponderation*reg288; T tmp_21_6=ponderation*reg245; T tmp_5_16=ponderation*reg203; T tmp_12_7=ponderation*reg767; T tmp_21_7=ponderation*reg827;
    T tmp_0_13=ponderation*reg518; T tmp_12_6=ponderation*reg435; T tmp_21_8=ponderation*reg825; T tmp_1_16=ponderation*reg1342; T tmp_12_22=ponderation*reg674;
    T tmp_12_21=ponderation*reg442; T tmp_20_20=ponderation*reg227; T tmp_5_8=ponderation*reg384; T tmp_12_20=ponderation*reg805; T tmp_20_21=ponderation*reg842;
    T tmp_1_8=ponderation*reg1375; T tmp_5_9=-reg695; T tmp_12_19=ponderation*reg803; T tmp_20_22=ponderation*reg840; T tmp_12_18=ponderation*reg423;
    T tmp_20_23=ponderation*reg224; T tmp_0_10=ponderation*reg209; T tmp_5_10=ponderation*reg238; T tmp_12_17=ponderation*reg796; T tmp_12_16=ponderation*reg793;
    T tmp_21_0=ponderation*reg250; T tmp_5_11=ponderation*reg430; T tmp_12_15=ponderation*reg340; T tmp_21_1=ponderation*reg837; T tmp_0_9=ponderation*reg1353;
    T tmp_5_12=ponderation*reg1330; T tmp_12_14=ponderation*reg291; T tmp_1_12=ponderation*reg1241; T tmp_1_2=ponderation*reg1194; T tmp_21_15=ponderation*reg205;
    T tmp_5_22=ponderation*reg497; T tmp_11_20=ponderation*reg484; T tmp_21_16=ponderation*reg856; T tmp_11_19=ponderation*reg774; T tmp_5_23=ponderation*reg479;
    T tmp_11_18=-reg670; T tmp_21_17=ponderation*reg853; T tmp_1_1=ponderation*reg537; T tmp_11_17=-reg669; T tmp_6_0=ponderation*reg494;
    T tmp_21_18=ponderation*reg535; T tmp_1_0=ponderation*reg1154; T tmp_11_16=ponderation*reg779; T tmp_21_19=ponderation*reg850; T tmp_11_15=-reg665;
    T tmp_6_1=ponderation*reg1211; T tmp_11_14=ponderation*reg466; T tmp_21_20=ponderation*reg848; T tmp_0_23=ponderation*reg1173; T tmp_11_13=ponderation*reg818;
    T tmp_5_17=ponderation*reg434; T tmp_12_5=ponderation*reg761; T tmp_21_9=ponderation*reg824; T tmp_12_4=ponderation*reg759; T tmp_1_15=ponderation*reg1318;
    T tmp_5_18=ponderation*reg463; T tmp_21_10=ponderation*reg864; T tmp_12_3=ponderation*reg389; T tmp_5_19=ponderation*reg465; T tmp_12_2=ponderation*reg751;
    T tmp_21_11=-reg636; T tmp_1_14=ponderation*reg242; T tmp_12_1=ponderation*reg747; T tmp_5_20=ponderation*reg393; T tmp_21_12=ponderation*reg211;
    T tmp_12_0=ponderation*reg402; T tmp_11_23=ponderation*reg410; T tmp_21_13=ponderation*reg860; T tmp_1_13=ponderation*reg234; T tmp_5_21=ponderation*reg629;
    T tmp_11_22=ponderation*reg784; T tmp_21_14=ponderation*reg858; T tmp_11_21=-reg675; T tmp_14_0=ponderation*reg717; T tmp_19_16=ponderation*reg872;
    T tmp_13_23=ponderation*reg715; T tmp_4_19=ponderation*reg342; T tmp_19_17=ponderation*reg870; T tmp_13_22=ponderation*reg380; T tmp_13_21=ponderation*reg709;
    T tmp_19_18=ponderation*reg321; T tmp_1_22=ponderation*reg324; T tmp_4_20=ponderation*reg1166; T tmp_13_20=ponderation*reg706; T tmp_19_19=ponderation*reg322;
    T tmp_13_19=ponderation*reg365; T tmp_1_21=ponderation*reg1295; T tmp_19_20=ponderation*reg867; T tmp_4_21=ponderation*reg1264; T tmp_13_18=ponderation*reg699;
    T tmp_19_21=ponderation*reg329; T tmp_13_17=ponderation*reg697; T tmp_1_20=ponderation*reg1239; T tmp_4_22=ponderation*reg351; T tmp_13_16=ponderation*reg743;
    T tmp_19_22=ponderation*reg330; T tmp_19_9=ponderation*reg882; T tmp_2_2=ponderation*reg370; T tmp_14_8=ponderation*reg339; T tmp_19_10=ponderation*reg364;
    T tmp_4_15=ponderation*reg1221; T tmp_14_7=ponderation*reg451; T tmp_19_11=ponderation*reg880; T tmp_14_6=ponderation*reg1044; T tmp_2_1=ponderation*reg1155;
    T tmp_4_16=ponderation*reg1215; T tmp_19_12=ponderation*reg877; T tmp_14_5=ponderation*reg341; T tmp_14_4=ponderation*reg490; T tmp_19_13=ponderation*reg374;
    T tmp_4_17=ponderation*reg1161; T tmp_2_0=ponderation*reg1184; T tmp_14_3=ponderation*reg1047; T tmp_19_14=ponderation*reg875; T tmp_14_2=ponderation*reg347;
    T tmp_14_1=ponderation*reg720; T tmp_19_15=ponderation*reg873; T tmp_4_18=ponderation*reg1144; T tmp_1_23=ponderation*reg1168; T tmp_5_3=ponderation*reg1244;
    T tmp_13_6=ponderation*reg666; T tmp_13_5=ponderation*reg660; T tmp_20_13=ponderation*reg373; T tmp_5_4=ponderation*reg360; T tmp_13_4=ponderation*reg491;
    T tmp_20_14=ponderation*reg236; T tmp_1_11=ponderation*reg1178; T tmp_13_3=ponderation*reg651; T tmp_20_15=ponderation*reg892; T tmp_5_5=ponderation*reg443;
    T tmp_13_2=ponderation*reg696; T tmp_20_16=ponderation*reg891; T tmp_13_1=ponderation*reg474; T tmp_1_10=ponderation*reg191; T tmp_20_17=ponderation*reg190;
    T tmp_5_6=ponderation*reg1172; T tmp_13_0=ponderation*reg683; T tmp_12_23=ponderation*reg679; T tmp_20_18=ponderation*reg889; T tmp_1_9=ponderation*reg1203;
    T tmp_5_7=ponderation*reg1159; T tmp_20_19=ponderation*reg844; T tmp_13_15=ponderation*reg741; T tmp_19_23=ponderation*reg902; T tmp_4_23=ponderation*reg1281;
    T tmp_13_14=ponderation*reg739; T tmp_20_3=ponderation*reg782; T tmp_1_19=ponderation*reg346; T tmp_13_13=ponderation*reg272; T tmp_20_4=ponderation*reg899;
    T tmp_5_0=ponderation*reg1233; T tmp_13_12=ponderation*reg732; T tmp_20_5=ponderation*reg336; T tmp_13_11=ponderation*reg729; T tmp_1_18=ponderation*reg1251;
    T tmp_20_6=ponderation*reg897; T tmp_5_1=ponderation*reg375; T tmp_13_10=ponderation*reg290; T tmp_20_7=ponderation*reg357; T tmp_13_9=ponderation*reg724;
    T tmp_5_2=ponderation*reg305; T tmp_13_8=ponderation*reg721; T tmp_20_8=ponderation*reg350; T tmp_1_17=ponderation*reg1260; T tmp_13_7=ponderation*reg458;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[0]+2) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_3;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_4;
    matrix(indices[0]+0,indices[1]+2) += tmp_0_5;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_7;
    matrix(indices[0]+0,indices[2]+2) += tmp_0_8;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_9;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_10;
    matrix(indices[0]+0,indices[3]+2) += tmp_0_11;
    matrix(indices[0]+0,indices[4]+0) += tmp_0_12;
    matrix(indices[0]+0,indices[4]+1) += tmp_0_13;
    matrix(indices[0]+0,indices[4]+2) += tmp_0_14;
    matrix(indices[0]+0,indices[5]+0) += tmp_0_15;
    matrix(indices[0]+0,indices[5]+1) += tmp_0_16;
    matrix(indices[0]+0,indices[5]+2) += tmp_0_17;
    matrix(indices[0]+0,indices[6]+0) += tmp_0_18;
    matrix(indices[0]+0,indices[6]+1) += tmp_0_19;
    matrix(indices[0]+0,indices[6]+2) += tmp_0_20;
    matrix(indices[0]+0,indices[7]+0) += tmp_0_21;
    matrix(indices[0]+0,indices[7]+1) += tmp_0_22;
    matrix(indices[0]+0,indices[7]+2) += tmp_0_23;
    matrix(indices[0]+1,indices[0]+0) += tmp_1_0;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[0]+2) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_3;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_4;
    matrix(indices[0]+1,indices[1]+2) += tmp_1_5;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_7;
    matrix(indices[0]+1,indices[2]+2) += tmp_1_8;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_9;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_10;
    matrix(indices[0]+1,indices[3]+2) += tmp_1_11;
    matrix(indices[0]+1,indices[4]+0) += tmp_1_12;
    matrix(indices[0]+1,indices[4]+1) += tmp_1_13;
    matrix(indices[0]+1,indices[4]+2) += tmp_1_14;
    matrix(indices[0]+1,indices[5]+0) += tmp_1_15;
    matrix(indices[0]+1,indices[5]+1) += tmp_1_16;
    matrix(indices[0]+1,indices[5]+2) += tmp_1_17;
    matrix(indices[0]+1,indices[6]+0) += tmp_1_18;
    matrix(indices[0]+1,indices[6]+1) += tmp_1_19;
    matrix(indices[0]+1,indices[6]+2) += tmp_1_20;
    matrix(indices[0]+1,indices[7]+0) += tmp_1_21;
    matrix(indices[0]+1,indices[7]+1) += tmp_1_22;
    matrix(indices[0]+1,indices[7]+2) += tmp_1_23;
    matrix(indices[0]+2,indices[0]+0) += tmp_2_0;
    matrix(indices[0]+2,indices[0]+1) += tmp_2_1;
    matrix(indices[0]+2,indices[0]+2) += tmp_2_2;
    matrix(indices[0]+2,indices[1]+0) += tmp_2_3;
    matrix(indices[0]+2,indices[1]+1) += tmp_2_4;
    matrix(indices[0]+2,indices[1]+2) += tmp_2_5;
    matrix(indices[0]+2,indices[2]+0) += tmp_2_6;
    matrix(indices[0]+2,indices[2]+1) += tmp_2_7;
    matrix(indices[0]+2,indices[2]+2) += tmp_2_8;
    matrix(indices[0]+2,indices[3]+0) += tmp_2_9;
    matrix(indices[0]+2,indices[3]+1) += tmp_2_10;
    matrix(indices[0]+2,indices[3]+2) += tmp_2_11;
    matrix(indices[0]+2,indices[4]+0) += tmp_2_12;
    matrix(indices[0]+2,indices[4]+1) += tmp_2_13;
    matrix(indices[0]+2,indices[4]+2) += tmp_2_14;
    matrix(indices[0]+2,indices[5]+0) += tmp_2_15;
    matrix(indices[0]+2,indices[5]+1) += tmp_2_16;
    matrix(indices[0]+2,indices[5]+2) += tmp_2_17;
    matrix(indices[0]+2,indices[6]+0) += tmp_2_18;
    matrix(indices[0]+2,indices[6]+1) += tmp_2_19;
    matrix(indices[0]+2,indices[6]+2) += tmp_2_20;
    matrix(indices[0]+2,indices[7]+0) += tmp_2_21;
    matrix(indices[0]+2,indices[7]+1) += tmp_2_22;
    matrix(indices[0]+2,indices[7]+2) += tmp_2_23;
    matrix(indices[1]+0,indices[0]+0) += tmp_3_0;
    matrix(indices[1]+0,indices[0]+1) += tmp_3_1;
    matrix(indices[1]+0,indices[0]+2) += tmp_3_2;
    matrix(indices[1]+0,indices[1]+0) += tmp_3_3;
    matrix(indices[1]+0,indices[1]+1) += tmp_3_4;
    matrix(indices[1]+0,indices[1]+2) += tmp_3_5;
    matrix(indices[1]+0,indices[2]+0) += tmp_3_6;
    matrix(indices[1]+0,indices[2]+1) += tmp_3_7;
    matrix(indices[1]+0,indices[2]+2) += tmp_3_8;
    matrix(indices[1]+0,indices[3]+0) += tmp_3_9;
    matrix(indices[1]+0,indices[3]+1) += tmp_3_10;
    matrix(indices[1]+0,indices[3]+2) += tmp_3_11;
    matrix(indices[1]+0,indices[4]+0) += tmp_3_12;
    matrix(indices[1]+0,indices[4]+1) += tmp_3_13;
    matrix(indices[1]+0,indices[4]+2) += tmp_3_14;
    matrix(indices[1]+0,indices[5]+0) += tmp_3_15;
    matrix(indices[1]+0,indices[5]+1) += tmp_3_16;
    matrix(indices[1]+0,indices[5]+2) += tmp_3_17;
    matrix(indices[1]+0,indices[6]+0) += tmp_3_18;
    matrix(indices[1]+0,indices[6]+1) += tmp_3_19;
    matrix(indices[1]+0,indices[6]+2) += tmp_3_20;
    matrix(indices[1]+0,indices[7]+0) += tmp_3_21;
    matrix(indices[1]+0,indices[7]+1) += tmp_3_22;
    matrix(indices[1]+0,indices[7]+2) += tmp_3_23;
    matrix(indices[1]+1,indices[0]+0) += tmp_4_0;
    matrix(indices[1]+1,indices[0]+1) += tmp_4_1;
    matrix(indices[1]+1,indices[0]+2) += tmp_4_2;
    matrix(indices[1]+1,indices[1]+0) += tmp_4_3;
    matrix(indices[1]+1,indices[1]+1) += tmp_4_4;
    matrix(indices[1]+1,indices[1]+2) += tmp_4_5;
    matrix(indices[1]+1,indices[2]+0) += tmp_4_6;
    matrix(indices[1]+1,indices[2]+1) += tmp_4_7;
    matrix(indices[1]+1,indices[2]+2) += tmp_4_8;
    matrix(indices[1]+1,indices[3]+0) += tmp_4_9;
    matrix(indices[1]+1,indices[3]+1) += tmp_4_10;
    matrix(indices[1]+1,indices[3]+2) += tmp_4_11;
    matrix(indices[1]+1,indices[4]+0) += tmp_4_12;
    matrix(indices[1]+1,indices[4]+1) += tmp_4_13;
    matrix(indices[1]+1,indices[4]+2) += tmp_4_14;
    matrix(indices[1]+1,indices[5]+0) += tmp_4_15;
    matrix(indices[1]+1,indices[5]+1) += tmp_4_16;
    matrix(indices[1]+1,indices[5]+2) += tmp_4_17;
    matrix(indices[1]+1,indices[6]+0) += tmp_4_18;
    matrix(indices[1]+1,indices[6]+1) += tmp_4_19;
    matrix(indices[1]+1,indices[6]+2) += tmp_4_20;
    matrix(indices[1]+1,indices[7]+0) += tmp_4_21;
    matrix(indices[1]+1,indices[7]+1) += tmp_4_22;
    matrix(indices[1]+1,indices[7]+2) += tmp_4_23;
    matrix(indices[1]+2,indices[0]+0) += tmp_5_0;
    matrix(indices[1]+2,indices[0]+1) += tmp_5_1;
    matrix(indices[1]+2,indices[0]+2) += tmp_5_2;
    matrix(indices[1]+2,indices[1]+0) += tmp_5_3;
    matrix(indices[1]+2,indices[1]+1) += tmp_5_4;
    matrix(indices[1]+2,indices[1]+2) += tmp_5_5;
    matrix(indices[1]+2,indices[2]+0) += tmp_5_6;
    matrix(indices[1]+2,indices[2]+1) += tmp_5_7;
    matrix(indices[1]+2,indices[2]+2) += tmp_5_8;
    matrix(indices[1]+2,indices[3]+0) += tmp_5_9;
    matrix(indices[1]+2,indices[3]+1) += tmp_5_10;
    matrix(indices[1]+2,indices[3]+2) += tmp_5_11;
    matrix(indices[1]+2,indices[4]+0) += tmp_5_12;
    matrix(indices[1]+2,indices[4]+1) += tmp_5_13;
    matrix(indices[1]+2,indices[4]+2) += tmp_5_14;
    matrix(indices[1]+2,indices[5]+0) += tmp_5_15;
    matrix(indices[1]+2,indices[5]+1) += tmp_5_16;
    matrix(indices[1]+2,indices[5]+2) += tmp_5_17;
    matrix(indices[1]+2,indices[6]+0) += tmp_5_18;
    matrix(indices[1]+2,indices[6]+1) += tmp_5_19;
    matrix(indices[1]+2,indices[6]+2) += tmp_5_20;
    matrix(indices[1]+2,indices[7]+0) += tmp_5_21;
    matrix(indices[1]+2,indices[7]+1) += tmp_5_22;
    matrix(indices[1]+2,indices[7]+2) += tmp_5_23;
    matrix(indices[2]+0,indices[0]+0) += tmp_6_0;
    matrix(indices[2]+0,indices[0]+1) += tmp_6_1;
    matrix(indices[2]+0,indices[0]+2) += tmp_6_2;
    matrix(indices[2]+0,indices[1]+0) += tmp_6_3;
    matrix(indices[2]+0,indices[1]+1) += tmp_6_4;
    matrix(indices[2]+0,indices[1]+2) += tmp_6_5;
    matrix(indices[2]+0,indices[2]+0) += tmp_6_6;
    matrix(indices[2]+0,indices[2]+1) += tmp_6_7;
    matrix(indices[2]+0,indices[2]+2) += tmp_6_8;
    matrix(indices[2]+0,indices[3]+0) += tmp_6_9;
    matrix(indices[2]+0,indices[3]+1) += tmp_6_10;
    matrix(indices[2]+0,indices[3]+2) += tmp_6_11;
    matrix(indices[2]+0,indices[4]+0) += tmp_6_12;
    matrix(indices[2]+0,indices[4]+1) += tmp_6_13;
    matrix(indices[2]+0,indices[4]+2) += tmp_6_14;
    matrix(indices[2]+0,indices[5]+0) += tmp_6_15;
    matrix(indices[2]+0,indices[5]+1) += tmp_6_16;
    matrix(indices[2]+0,indices[5]+2) += tmp_6_17;
    matrix(indices[2]+0,indices[6]+0) += tmp_6_18;
    matrix(indices[2]+0,indices[6]+1) += tmp_6_19;
    matrix(indices[2]+0,indices[6]+2) += tmp_6_20;
    matrix(indices[2]+0,indices[7]+0) += tmp_6_21;
    matrix(indices[2]+0,indices[7]+1) += tmp_6_22;
    matrix(indices[2]+0,indices[7]+2) += tmp_6_23;
    matrix(indices[2]+1,indices[0]+0) += tmp_7_0;
    matrix(indices[2]+1,indices[0]+1) += tmp_7_1;
    matrix(indices[2]+1,indices[0]+2) += tmp_7_2;
    matrix(indices[2]+1,indices[1]+0) += tmp_7_3;
    matrix(indices[2]+1,indices[1]+1) += tmp_7_4;
    matrix(indices[2]+1,indices[1]+2) += tmp_7_5;
    matrix(indices[2]+1,indices[2]+0) += tmp_7_6;
    matrix(indices[2]+1,indices[2]+1) += tmp_7_7;
    matrix(indices[2]+1,indices[2]+2) += tmp_7_8;
    matrix(indices[2]+1,indices[3]+0) += tmp_7_9;
    matrix(indices[2]+1,indices[3]+1) += tmp_7_10;
    matrix(indices[2]+1,indices[3]+2) += tmp_7_11;
    matrix(indices[2]+1,indices[4]+0) += tmp_7_12;
    matrix(indices[2]+1,indices[4]+1) += tmp_7_13;
    matrix(indices[2]+1,indices[4]+2) += tmp_7_14;
    matrix(indices[2]+1,indices[5]+0) += tmp_7_15;
    matrix(indices[2]+1,indices[5]+1) += tmp_7_16;
    matrix(indices[2]+1,indices[5]+2) += tmp_7_17;
    matrix(indices[2]+1,indices[6]+0) += tmp_7_18;
    matrix(indices[2]+1,indices[6]+1) += tmp_7_19;
    matrix(indices[2]+1,indices[6]+2) += tmp_7_20;
    matrix(indices[2]+1,indices[7]+0) += tmp_7_21;
    matrix(indices[2]+1,indices[7]+1) += tmp_7_22;
    matrix(indices[2]+1,indices[7]+2) += tmp_7_23;
    matrix(indices[2]+2,indices[0]+0) += tmp_8_0;
    matrix(indices[2]+2,indices[0]+1) += tmp_8_1;
    matrix(indices[2]+2,indices[0]+2) += tmp_8_2;
    matrix(indices[2]+2,indices[1]+0) += tmp_8_3;
    matrix(indices[2]+2,indices[1]+1) += tmp_8_4;
    matrix(indices[2]+2,indices[1]+2) += tmp_8_5;
    matrix(indices[2]+2,indices[2]+0) += tmp_8_6;
    matrix(indices[2]+2,indices[2]+1) += tmp_8_7;
    matrix(indices[2]+2,indices[2]+2) += tmp_8_8;
    matrix(indices[2]+2,indices[3]+0) += tmp_8_9;
    matrix(indices[2]+2,indices[3]+1) += tmp_8_10;
    matrix(indices[2]+2,indices[3]+2) += tmp_8_11;
    matrix(indices[2]+2,indices[4]+0) += tmp_8_12;
    matrix(indices[2]+2,indices[4]+1) += tmp_8_13;
    matrix(indices[2]+2,indices[4]+2) += tmp_8_14;
    matrix(indices[2]+2,indices[5]+0) += tmp_8_15;
    matrix(indices[2]+2,indices[5]+1) += tmp_8_16;
    matrix(indices[2]+2,indices[5]+2) += tmp_8_17;
    matrix(indices[2]+2,indices[6]+0) += tmp_8_18;
    matrix(indices[2]+2,indices[6]+1) += tmp_8_19;
    matrix(indices[2]+2,indices[6]+2) += tmp_8_20;
    matrix(indices[2]+2,indices[7]+0) += tmp_8_21;
    matrix(indices[2]+2,indices[7]+1) += tmp_8_22;
    matrix(indices[2]+2,indices[7]+2) += tmp_8_23;
    matrix(indices[3]+0,indices[0]+0) += tmp_9_0;
    matrix(indices[3]+0,indices[0]+1) += tmp_9_1;
    matrix(indices[3]+0,indices[0]+2) += tmp_9_2;
    matrix(indices[3]+0,indices[1]+0) += tmp_9_3;
    matrix(indices[3]+0,indices[1]+1) += tmp_9_4;
    matrix(indices[3]+0,indices[1]+2) += tmp_9_5;
    matrix(indices[3]+0,indices[2]+0) += tmp_9_6;
    matrix(indices[3]+0,indices[2]+1) += tmp_9_7;
    matrix(indices[3]+0,indices[2]+2) += tmp_9_8;
    matrix(indices[3]+0,indices[3]+0) += tmp_9_9;
    matrix(indices[3]+0,indices[3]+1) += tmp_9_10;
    matrix(indices[3]+0,indices[3]+2) += tmp_9_11;
    matrix(indices[3]+0,indices[4]+0) += tmp_9_12;
    matrix(indices[3]+0,indices[4]+1) += tmp_9_13;
    matrix(indices[3]+0,indices[4]+2) += tmp_9_14;
    matrix(indices[3]+0,indices[5]+0) += tmp_9_15;
    matrix(indices[3]+0,indices[5]+1) += tmp_9_16;
    matrix(indices[3]+0,indices[5]+2) += tmp_9_17;
    matrix(indices[3]+0,indices[6]+0) += tmp_9_18;
    matrix(indices[3]+0,indices[6]+1) += tmp_9_19;
    matrix(indices[3]+0,indices[6]+2) += tmp_9_20;
    matrix(indices[3]+0,indices[7]+0) += tmp_9_21;
    matrix(indices[3]+0,indices[7]+1) += tmp_9_22;
    matrix(indices[3]+0,indices[7]+2) += tmp_9_23;
    matrix(indices[3]+1,indices[0]+0) += tmp_10_0;
    matrix(indices[3]+1,indices[0]+1) += tmp_10_1;
    matrix(indices[3]+1,indices[0]+2) += tmp_10_2;
    matrix(indices[3]+1,indices[1]+0) += tmp_10_3;
    matrix(indices[3]+1,indices[1]+1) += tmp_10_4;
    matrix(indices[3]+1,indices[1]+2) += tmp_10_5;
    matrix(indices[3]+1,indices[2]+0) += tmp_10_6;
    matrix(indices[3]+1,indices[2]+1) += tmp_10_7;
    matrix(indices[3]+1,indices[2]+2) += tmp_10_8;
    matrix(indices[3]+1,indices[3]+0) += tmp_10_9;
    matrix(indices[3]+1,indices[3]+1) += tmp_10_10;
    matrix(indices[3]+1,indices[3]+2) += tmp_10_11;
    matrix(indices[3]+1,indices[4]+0) += tmp_10_12;
    matrix(indices[3]+1,indices[4]+1) += tmp_10_13;
    matrix(indices[3]+1,indices[4]+2) += tmp_10_14;
    matrix(indices[3]+1,indices[5]+0) += tmp_10_15;
    matrix(indices[3]+1,indices[5]+1) += tmp_10_16;
    matrix(indices[3]+1,indices[5]+2) += tmp_10_17;
    matrix(indices[3]+1,indices[6]+0) += tmp_10_18;
    matrix(indices[3]+1,indices[6]+1) += tmp_10_19;
    matrix(indices[3]+1,indices[6]+2) += tmp_10_20;
    matrix(indices[3]+1,indices[7]+0) += tmp_10_21;
    matrix(indices[3]+1,indices[7]+1) += tmp_10_22;
    matrix(indices[3]+1,indices[7]+2) += tmp_10_23;
    matrix(indices[3]+2,indices[0]+0) += tmp_11_0;
    matrix(indices[3]+2,indices[0]+1) += tmp_11_1;
    matrix(indices[3]+2,indices[0]+2) += tmp_11_2;
    matrix(indices[3]+2,indices[1]+0) += tmp_11_3;
    matrix(indices[3]+2,indices[1]+1) += tmp_11_4;
    matrix(indices[3]+2,indices[1]+2) += tmp_11_5;
    matrix(indices[3]+2,indices[2]+0) += tmp_11_6;
    matrix(indices[3]+2,indices[2]+1) += tmp_11_7;
    matrix(indices[3]+2,indices[2]+2) += tmp_11_8;
    matrix(indices[3]+2,indices[3]+0) += tmp_11_9;
    matrix(indices[3]+2,indices[3]+1) += tmp_11_10;
    matrix(indices[3]+2,indices[3]+2) += tmp_11_11;
    matrix(indices[3]+2,indices[4]+0) += tmp_11_12;
    matrix(indices[3]+2,indices[4]+1) += tmp_11_13;
    matrix(indices[3]+2,indices[4]+2) += tmp_11_14;
    matrix(indices[3]+2,indices[5]+0) += tmp_11_15;
    matrix(indices[3]+2,indices[5]+1) += tmp_11_16;
    matrix(indices[3]+2,indices[5]+2) += tmp_11_17;
    matrix(indices[3]+2,indices[6]+0) += tmp_11_18;
    matrix(indices[3]+2,indices[6]+1) += tmp_11_19;
    matrix(indices[3]+2,indices[6]+2) += tmp_11_20;
    matrix(indices[3]+2,indices[7]+0) += tmp_11_21;
    matrix(indices[3]+2,indices[7]+1) += tmp_11_22;
    matrix(indices[3]+2,indices[7]+2) += tmp_11_23;
    matrix(indices[4]+0,indices[0]+0) += tmp_12_0;
    matrix(indices[4]+0,indices[0]+1) += tmp_12_1;
    matrix(indices[4]+0,indices[0]+2) += tmp_12_2;
    matrix(indices[4]+0,indices[1]+0) += tmp_12_3;
    matrix(indices[4]+0,indices[1]+1) += tmp_12_4;
    matrix(indices[4]+0,indices[1]+2) += tmp_12_5;
    matrix(indices[4]+0,indices[2]+0) += tmp_12_6;
    matrix(indices[4]+0,indices[2]+1) += tmp_12_7;
    matrix(indices[4]+0,indices[2]+2) += tmp_12_8;
    matrix(indices[4]+0,indices[3]+0) += tmp_12_9;
    matrix(indices[4]+0,indices[3]+1) += tmp_12_10;
    matrix(indices[4]+0,indices[3]+2) += tmp_12_11;
    matrix(indices[4]+0,indices[4]+0) += tmp_12_12;
    matrix(indices[4]+0,indices[4]+1) += tmp_12_13;
    matrix(indices[4]+0,indices[4]+2) += tmp_12_14;
    matrix(indices[4]+0,indices[5]+0) += tmp_12_15;
    matrix(indices[4]+0,indices[5]+1) += tmp_12_16;
    matrix(indices[4]+0,indices[5]+2) += tmp_12_17;
    matrix(indices[4]+0,indices[6]+0) += tmp_12_18;
    matrix(indices[4]+0,indices[6]+1) += tmp_12_19;
    matrix(indices[4]+0,indices[6]+2) += tmp_12_20;
    matrix(indices[4]+0,indices[7]+0) += tmp_12_21;
    matrix(indices[4]+0,indices[7]+1) += tmp_12_22;
    matrix(indices[4]+0,indices[7]+2) += tmp_12_23;
    matrix(indices[4]+1,indices[0]+0) += tmp_13_0;
    matrix(indices[4]+1,indices[0]+1) += tmp_13_1;
    matrix(indices[4]+1,indices[0]+2) += tmp_13_2;
    matrix(indices[4]+1,indices[1]+0) += tmp_13_3;
    matrix(indices[4]+1,indices[1]+1) += tmp_13_4;
    matrix(indices[4]+1,indices[1]+2) += tmp_13_5;
    matrix(indices[4]+1,indices[2]+0) += tmp_13_6;
    matrix(indices[4]+1,indices[2]+1) += tmp_13_7;
    matrix(indices[4]+1,indices[2]+2) += tmp_13_8;
    matrix(indices[4]+1,indices[3]+0) += tmp_13_9;
    matrix(indices[4]+1,indices[3]+1) += tmp_13_10;
    matrix(indices[4]+1,indices[3]+2) += tmp_13_11;
    matrix(indices[4]+1,indices[4]+0) += tmp_13_12;
    matrix(indices[4]+1,indices[4]+1) += tmp_13_13;
    matrix(indices[4]+1,indices[4]+2) += tmp_13_14;
    matrix(indices[4]+1,indices[5]+0) += tmp_13_15;
    matrix(indices[4]+1,indices[5]+1) += tmp_13_16;
    matrix(indices[4]+1,indices[5]+2) += tmp_13_17;
    matrix(indices[4]+1,indices[6]+0) += tmp_13_18;
    matrix(indices[4]+1,indices[6]+1) += tmp_13_19;
    matrix(indices[4]+1,indices[6]+2) += tmp_13_20;
    matrix(indices[4]+1,indices[7]+0) += tmp_13_21;
    matrix(indices[4]+1,indices[7]+1) += tmp_13_22;
    matrix(indices[4]+1,indices[7]+2) += tmp_13_23;
    matrix(indices[4]+2,indices[0]+0) += tmp_14_0;
    matrix(indices[4]+2,indices[0]+1) += tmp_14_1;
    matrix(indices[4]+2,indices[0]+2) += tmp_14_2;
    matrix(indices[4]+2,indices[1]+0) += tmp_14_3;
    matrix(indices[4]+2,indices[1]+1) += tmp_14_4;
    matrix(indices[4]+2,indices[1]+2) += tmp_14_5;
    matrix(indices[4]+2,indices[2]+0) += tmp_14_6;
    matrix(indices[4]+2,indices[2]+1) += tmp_14_7;
    matrix(indices[4]+2,indices[2]+2) += tmp_14_8;
    matrix(indices[4]+2,indices[3]+0) += tmp_14_9;
    matrix(indices[4]+2,indices[3]+1) += tmp_14_10;
    matrix(indices[4]+2,indices[3]+2) += tmp_14_11;
    matrix(indices[4]+2,indices[4]+0) += tmp_14_12;
    matrix(indices[4]+2,indices[4]+1) += tmp_14_13;
    matrix(indices[4]+2,indices[4]+2) += tmp_14_14;
    matrix(indices[4]+2,indices[5]+0) += tmp_14_15;
    matrix(indices[4]+2,indices[5]+1) += tmp_14_16;
    matrix(indices[4]+2,indices[5]+2) += tmp_14_17;
    matrix(indices[4]+2,indices[6]+0) += tmp_14_18;
    matrix(indices[4]+2,indices[6]+1) += tmp_14_19;
    matrix(indices[4]+2,indices[6]+2) += tmp_14_20;
    matrix(indices[4]+2,indices[7]+0) += tmp_14_21;
    matrix(indices[4]+2,indices[7]+1) += tmp_14_22;
    matrix(indices[4]+2,indices[7]+2) += tmp_14_23;
    matrix(indices[5]+0,indices[0]+0) += tmp_15_0;
    matrix(indices[5]+0,indices[0]+1) += tmp_15_1;
    matrix(indices[5]+0,indices[0]+2) += tmp_15_2;
    matrix(indices[5]+0,indices[1]+0) += tmp_15_3;
    matrix(indices[5]+0,indices[1]+1) += tmp_15_4;
    matrix(indices[5]+0,indices[1]+2) += tmp_15_5;
    matrix(indices[5]+0,indices[2]+0) += tmp_15_6;
    matrix(indices[5]+0,indices[2]+1) += tmp_15_7;
    matrix(indices[5]+0,indices[2]+2) += tmp_15_8;
    matrix(indices[5]+0,indices[3]+0) += tmp_15_9;
    matrix(indices[5]+0,indices[3]+1) += tmp_15_10;
    matrix(indices[5]+0,indices[3]+2) += tmp_15_11;
    matrix(indices[5]+0,indices[4]+0) += tmp_15_12;
    matrix(indices[5]+0,indices[4]+1) += tmp_15_13;
    matrix(indices[5]+0,indices[4]+2) += tmp_15_14;
    matrix(indices[5]+0,indices[5]+0) += tmp_15_15;
    matrix(indices[5]+0,indices[5]+1) += tmp_15_16;
    matrix(indices[5]+0,indices[5]+2) += tmp_15_17;
    matrix(indices[5]+0,indices[6]+0) += tmp_15_18;
    matrix(indices[5]+0,indices[6]+1) += tmp_15_19;
    matrix(indices[5]+0,indices[6]+2) += tmp_15_20;
    matrix(indices[5]+0,indices[7]+0) += tmp_15_21;
    matrix(indices[5]+0,indices[7]+1) += tmp_15_22;
    matrix(indices[5]+0,indices[7]+2) += tmp_15_23;
    matrix(indices[5]+1,indices[0]+0) += tmp_16_0;
    matrix(indices[5]+1,indices[0]+1) += tmp_16_1;
    matrix(indices[5]+1,indices[0]+2) += tmp_16_2;
    matrix(indices[5]+1,indices[1]+0) += tmp_16_3;
    matrix(indices[5]+1,indices[1]+1) += tmp_16_4;
    matrix(indices[5]+1,indices[1]+2) += tmp_16_5;
    matrix(indices[5]+1,indices[2]+0) += tmp_16_6;
    matrix(indices[5]+1,indices[2]+1) += tmp_16_7;
    matrix(indices[5]+1,indices[2]+2) += tmp_16_8;
    matrix(indices[5]+1,indices[3]+0) += tmp_16_9;
    matrix(indices[5]+1,indices[3]+1) += tmp_16_10;
    matrix(indices[5]+1,indices[3]+2) += tmp_16_11;
    matrix(indices[5]+1,indices[4]+0) += tmp_16_12;
    matrix(indices[5]+1,indices[4]+1) += tmp_16_13;
    matrix(indices[5]+1,indices[4]+2) += tmp_16_14;
    matrix(indices[5]+1,indices[5]+0) += tmp_16_15;
    matrix(indices[5]+1,indices[5]+1) += tmp_16_16;
    matrix(indices[5]+1,indices[5]+2) += tmp_16_17;
    matrix(indices[5]+1,indices[6]+0) += tmp_16_18;
    matrix(indices[5]+1,indices[6]+1) += tmp_16_19;
    matrix(indices[5]+1,indices[6]+2) += tmp_16_20;
    matrix(indices[5]+1,indices[7]+0) += tmp_16_21;
    matrix(indices[5]+1,indices[7]+1) += tmp_16_22;
    matrix(indices[5]+1,indices[7]+2) += tmp_16_23;
    matrix(indices[5]+2,indices[0]+0) += tmp_17_0;
    matrix(indices[5]+2,indices[0]+1) += tmp_17_1;
    matrix(indices[5]+2,indices[0]+2) += tmp_17_2;
    matrix(indices[5]+2,indices[1]+0) += tmp_17_3;
    matrix(indices[5]+2,indices[1]+1) += tmp_17_4;
    matrix(indices[5]+2,indices[1]+2) += tmp_17_5;
    matrix(indices[5]+2,indices[2]+0) += tmp_17_6;
    matrix(indices[5]+2,indices[2]+1) += tmp_17_7;
    matrix(indices[5]+2,indices[2]+2) += tmp_17_8;
    matrix(indices[5]+2,indices[3]+0) += tmp_17_9;
    matrix(indices[5]+2,indices[3]+1) += tmp_17_10;
    matrix(indices[5]+2,indices[3]+2) += tmp_17_11;
    matrix(indices[5]+2,indices[4]+0) += tmp_17_12;
    matrix(indices[5]+2,indices[4]+1) += tmp_17_13;
    matrix(indices[5]+2,indices[4]+2) += tmp_17_14;
    matrix(indices[5]+2,indices[5]+0) += tmp_17_15;
    matrix(indices[5]+2,indices[5]+1) += tmp_17_16;
    matrix(indices[5]+2,indices[5]+2) += tmp_17_17;
    matrix(indices[5]+2,indices[6]+0) += tmp_17_18;
    matrix(indices[5]+2,indices[6]+1) += tmp_17_19;
    matrix(indices[5]+2,indices[6]+2) += tmp_17_20;
    matrix(indices[5]+2,indices[7]+0) += tmp_17_21;
    matrix(indices[5]+2,indices[7]+1) += tmp_17_22;
    matrix(indices[5]+2,indices[7]+2) += tmp_17_23;
    matrix(indices[6]+0,indices[0]+0) += tmp_18_0;
    matrix(indices[6]+0,indices[0]+1) += tmp_18_1;
    matrix(indices[6]+0,indices[0]+2) += tmp_18_2;
    matrix(indices[6]+0,indices[1]+0) += tmp_18_3;
    matrix(indices[6]+0,indices[1]+1) += tmp_18_4;
    matrix(indices[6]+0,indices[1]+2) += tmp_18_5;
    matrix(indices[6]+0,indices[2]+0) += tmp_18_6;
    matrix(indices[6]+0,indices[2]+1) += tmp_18_7;
    matrix(indices[6]+0,indices[2]+2) += tmp_18_8;
    matrix(indices[6]+0,indices[3]+0) += tmp_18_9;
    matrix(indices[6]+0,indices[3]+1) += tmp_18_10;
    matrix(indices[6]+0,indices[3]+2) += tmp_18_11;
    matrix(indices[6]+0,indices[4]+0) += tmp_18_12;
    matrix(indices[6]+0,indices[4]+1) += tmp_18_13;
    matrix(indices[6]+0,indices[4]+2) += tmp_18_14;
    matrix(indices[6]+0,indices[5]+0) += tmp_18_15;
    matrix(indices[6]+0,indices[5]+1) += tmp_18_16;
    matrix(indices[6]+0,indices[5]+2) += tmp_18_17;
    matrix(indices[6]+0,indices[6]+0) += tmp_18_18;
    matrix(indices[6]+0,indices[6]+1) += tmp_18_19;
    matrix(indices[6]+0,indices[6]+2) += tmp_18_20;
    matrix(indices[6]+0,indices[7]+0) += tmp_18_21;
    matrix(indices[6]+0,indices[7]+1) += tmp_18_22;
    matrix(indices[6]+0,indices[7]+2) += tmp_18_23;
    matrix(indices[6]+1,indices[0]+0) += tmp_19_0;
    matrix(indices[6]+1,indices[0]+1) += tmp_19_1;
    matrix(indices[6]+1,indices[0]+2) += tmp_19_2;
    matrix(indices[6]+1,indices[1]+0) += tmp_19_3;
    matrix(indices[6]+1,indices[1]+1) += tmp_19_4;
    matrix(indices[6]+1,indices[1]+2) += tmp_19_5;
    matrix(indices[6]+1,indices[2]+0) += tmp_19_6;
    matrix(indices[6]+1,indices[2]+1) += tmp_19_7;
    matrix(indices[6]+1,indices[2]+2) += tmp_19_8;
    matrix(indices[6]+1,indices[3]+0) += tmp_19_9;
    matrix(indices[6]+1,indices[3]+1) += tmp_19_10;
    matrix(indices[6]+1,indices[3]+2) += tmp_19_11;
    matrix(indices[6]+1,indices[4]+0) += tmp_19_12;
    matrix(indices[6]+1,indices[4]+1) += tmp_19_13;
    matrix(indices[6]+1,indices[4]+2) += tmp_19_14;
    matrix(indices[6]+1,indices[5]+0) += tmp_19_15;
    matrix(indices[6]+1,indices[5]+1) += tmp_19_16;
    matrix(indices[6]+1,indices[5]+2) += tmp_19_17;
    matrix(indices[6]+1,indices[6]+0) += tmp_19_18;
    matrix(indices[6]+1,indices[6]+1) += tmp_19_19;
    matrix(indices[6]+1,indices[6]+2) += tmp_19_20;
    matrix(indices[6]+1,indices[7]+0) += tmp_19_21;
    matrix(indices[6]+1,indices[7]+1) += tmp_19_22;
    matrix(indices[6]+1,indices[7]+2) += tmp_19_23;
    matrix(indices[6]+2,indices[0]+0) += tmp_20_0;
    matrix(indices[6]+2,indices[0]+1) += tmp_20_1;
    matrix(indices[6]+2,indices[0]+2) += tmp_20_2;
    matrix(indices[6]+2,indices[1]+0) += tmp_20_3;
    matrix(indices[6]+2,indices[1]+1) += tmp_20_4;
    matrix(indices[6]+2,indices[1]+2) += tmp_20_5;
    matrix(indices[6]+2,indices[2]+0) += tmp_20_6;
    matrix(indices[6]+2,indices[2]+1) += tmp_20_7;
    matrix(indices[6]+2,indices[2]+2) += tmp_20_8;
    matrix(indices[6]+2,indices[3]+0) += tmp_20_9;
    matrix(indices[6]+2,indices[3]+1) += tmp_20_10;
    matrix(indices[6]+2,indices[3]+2) += tmp_20_11;
    matrix(indices[6]+2,indices[4]+0) += tmp_20_12;
    matrix(indices[6]+2,indices[4]+1) += tmp_20_13;
    matrix(indices[6]+2,indices[4]+2) += tmp_20_14;
    matrix(indices[6]+2,indices[5]+0) += tmp_20_15;
    matrix(indices[6]+2,indices[5]+1) += tmp_20_16;
    matrix(indices[6]+2,indices[5]+2) += tmp_20_17;
    matrix(indices[6]+2,indices[6]+0) += tmp_20_18;
    matrix(indices[6]+2,indices[6]+1) += tmp_20_19;
    matrix(indices[6]+2,indices[6]+2) += tmp_20_20;
    matrix(indices[6]+2,indices[7]+0) += tmp_20_21;
    matrix(indices[6]+2,indices[7]+1) += tmp_20_22;
    matrix(indices[6]+2,indices[7]+2) += tmp_20_23;
    matrix(indices[7]+0,indices[0]+0) += tmp_21_0;
    matrix(indices[7]+0,indices[0]+1) += tmp_21_1;
    matrix(indices[7]+0,indices[0]+2) += tmp_21_2;
    matrix(indices[7]+0,indices[1]+0) += tmp_21_3;
    matrix(indices[7]+0,indices[1]+1) += tmp_21_4;
    matrix(indices[7]+0,indices[1]+2) += tmp_21_5;
    matrix(indices[7]+0,indices[2]+0) += tmp_21_6;
    matrix(indices[7]+0,indices[2]+1) += tmp_21_7;
    matrix(indices[7]+0,indices[2]+2) += tmp_21_8;
    matrix(indices[7]+0,indices[3]+0) += tmp_21_9;
    matrix(indices[7]+0,indices[3]+1) += tmp_21_10;
    matrix(indices[7]+0,indices[3]+2) += tmp_21_11;
    matrix(indices[7]+0,indices[4]+0) += tmp_21_12;
    matrix(indices[7]+0,indices[4]+1) += tmp_21_13;
    matrix(indices[7]+0,indices[4]+2) += tmp_21_14;
    matrix(indices[7]+0,indices[5]+0) += tmp_21_15;
    matrix(indices[7]+0,indices[5]+1) += tmp_21_16;
    matrix(indices[7]+0,indices[5]+2) += tmp_21_17;
    matrix(indices[7]+0,indices[6]+0) += tmp_21_18;
    matrix(indices[7]+0,indices[6]+1) += tmp_21_19;
    matrix(indices[7]+0,indices[6]+2) += tmp_21_20;
    matrix(indices[7]+0,indices[7]+0) += tmp_21_21;
    matrix(indices[7]+0,indices[7]+1) += tmp_21_22;
    matrix(indices[7]+0,indices[7]+2) += tmp_21_23;
    matrix(indices[7]+1,indices[0]+0) += tmp_22_0;
    matrix(indices[7]+1,indices[0]+1) += tmp_22_1;
    matrix(indices[7]+1,indices[0]+2) += tmp_22_2;
    matrix(indices[7]+1,indices[1]+0) += tmp_22_3;
    matrix(indices[7]+1,indices[1]+1) += tmp_22_4;
    matrix(indices[7]+1,indices[1]+2) += tmp_22_5;
    matrix(indices[7]+1,indices[2]+0) += tmp_22_6;
    matrix(indices[7]+1,indices[2]+1) += tmp_22_7;
    matrix(indices[7]+1,indices[2]+2) += tmp_22_8;
    matrix(indices[7]+1,indices[3]+0) += tmp_22_9;
    matrix(indices[7]+1,indices[3]+1) += tmp_22_10;
    matrix(indices[7]+1,indices[3]+2) += tmp_22_11;
    matrix(indices[7]+1,indices[4]+0) += tmp_22_12;
    matrix(indices[7]+1,indices[4]+1) += tmp_22_13;
    matrix(indices[7]+1,indices[4]+2) += tmp_22_14;
    matrix(indices[7]+1,indices[5]+0) += tmp_22_15;
    matrix(indices[7]+1,indices[5]+1) += tmp_22_16;
    matrix(indices[7]+1,indices[5]+2) += tmp_22_17;
    matrix(indices[7]+1,indices[6]+0) += tmp_22_18;
    matrix(indices[7]+1,indices[6]+1) += tmp_22_19;
    matrix(indices[7]+1,indices[6]+2) += tmp_22_20;
    matrix(indices[7]+1,indices[7]+0) += tmp_22_21;
    matrix(indices[7]+1,indices[7]+1) += tmp_22_22;
    matrix(indices[7]+1,indices[7]+2) += tmp_22_23;
    matrix(indices[7]+2,indices[0]+0) += tmp_23_0;
    matrix(indices[7]+2,indices[0]+1) += tmp_23_1;
    matrix(indices[7]+2,indices[0]+2) += tmp_23_2;
    matrix(indices[7]+2,indices[1]+0) += tmp_23_3;
    matrix(indices[7]+2,indices[1]+1) += tmp_23_4;
    matrix(indices[7]+2,indices[1]+2) += tmp_23_5;
    matrix(indices[7]+2,indices[2]+0) += tmp_23_6;
    matrix(indices[7]+2,indices[2]+1) += tmp_23_7;
    matrix(indices[7]+2,indices[2]+2) += tmp_23_8;
    matrix(indices[7]+2,indices[3]+0) += tmp_23_9;
    matrix(indices[7]+2,indices[3]+1) += tmp_23_10;
    matrix(indices[7]+2,indices[3]+2) += tmp_23_11;
    matrix(indices[7]+2,indices[4]+0) += tmp_23_12;
    matrix(indices[7]+2,indices[4]+1) += tmp_23_13;
    matrix(indices[7]+2,indices[4]+2) += tmp_23_14;
    matrix(indices[7]+2,indices[5]+0) += tmp_23_15;
    matrix(indices[7]+2,indices[5]+1) += tmp_23_16;
    matrix(indices[7]+2,indices[5]+2) += tmp_23_17;
    matrix(indices[7]+2,indices[6]+0) += tmp_23_18;
    matrix(indices[7]+2,indices[6]+1) += tmp_23_19;
    matrix(indices[7]+2,indices[6]+2) += tmp_23_20;
    matrix(indices[7]+2,indices[7]+0) += tmp_23_21;
    matrix(indices[7]+2,indices[7]+1) += tmp_23_22;
    matrix(indices[7]+2,indices[7]+2) += tmp_23_23;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1-var_inter[1]; T reg1=1-var_inter[2]; T reg2=1-var_inter[0]; T reg3=reg0*reg2; T reg4=reg0*reg1;
    T reg5=reg0*var_inter[0]; T reg6=reg1*var_inter[0]; T reg7=reg1*reg2; T reg8=elem.pos(1)[1]*reg4; T reg9=elem.pos(0)[1]*reg4;
    T reg10=reg1*var_inter[1]; T reg11=elem.pos(1)[1]*reg5; T reg12=elem.pos(0)[1]*reg7; T reg13=elem.pos(0)[1]*reg3; T reg14=elem.pos(1)[2]*reg6;
    T reg15=elem.pos(1)[2]*reg4; T reg16=elem.pos(1)[1]*reg6; T reg17=elem.pos(0)[2]*reg4; T reg18=elem.pos(0)[2]*reg7; T reg19=var_inter[0]*var_inter[1];
    T reg20=elem.pos(0)[2]*reg3; T reg21=reg5*elem.pos(1)[2]; T reg22=reg19*elem.pos(2)[2]; T reg23=reg21+reg20; T reg24=elem.pos(2)[2]*reg6;
    T reg25=reg14+reg18; T reg26=reg11+reg13; T reg27=reg19*elem.pos(2)[1]; reg8=reg8-reg9; T reg28=elem.pos(2)[1]*reg10;
    T reg29=reg16+reg12; T reg30=elem.pos(2)[1]*reg6; T reg31=var_inter[1]*reg2; reg15=reg15-reg17; T reg32=elem.pos(2)[2]*reg10;
    reg32=reg15+reg32; reg15=elem.pos(3)[1]*reg7; reg28=reg8+reg28; reg30=reg30-reg29; reg8=elem.pos(3)[1]*reg10;
    T reg33=elem.pos(3)[2]*reg7; T reg34=reg23+reg22; T reg35=elem.pos(1)[0]*reg6; T reg36=elem.pos(0)[0]*reg7; T reg37=elem.pos(3)[2]*reg10;
    reg24=reg24-reg25; T reg38=reg0*var_inter[2]; T reg39=elem.pos(0)[0]*reg4; T reg40=var_inter[2]*reg2; T reg41=elem.pos(1)[0]*reg4;
    T reg42=reg26+reg27; T reg43=reg31*elem.pos(3)[2]; T reg44=reg31*elem.pos(3)[1]; T reg45=reg42+reg44; reg28=reg28-reg8;
    T reg46=elem.pos(4)[1]*reg38; reg15=reg30+reg15; reg30=elem.pos(4)[1]*reg40; T reg47=elem.pos(4)[2]*reg3; T reg48=elem.pos(4)[2]*reg38;
    T reg49=elem.pos(4)[1]*reg3; reg32=reg32-reg37; T reg50=var_inter[2]*var_inter[0]; reg41=reg41-reg39; T reg51=elem.pos(2)[0]*reg10;
    T reg52=elem.pos(2)[0]*reg6; T reg53=reg35+reg36; T reg54=elem.pos(4)[2]*reg40; reg24=reg33+reg24; reg33=elem.pos(1)[0]*reg5;
    T reg55=elem.pos(0)[0]*reg3; T reg56=reg43+reg34; T reg57=reg5*elem.pos(5)[1]; reg15=reg15-reg30; T reg58=elem.pos(5)[2]*reg38;
    reg32=reg32-reg48; reg52=reg52-reg53; T reg59=elem.pos(3)[0]*reg7; reg49=reg49-reg45; T reg60=elem.pos(3)[0]*reg10;
    reg51=reg41+reg51; reg41=elem.pos(5)[2]*reg50; reg24=reg24-reg54; T reg61=var_inter[2]*var_inter[1]; T reg62=reg33+reg55;
    T reg63=reg19*elem.pos(2)[0]; reg28=reg28-reg46; T reg64=elem.pos(5)[1]*reg38; T reg65=elem.pos(5)[1]*reg50; reg47=reg47-reg56;
    T reg66=reg5*elem.pos(5)[2]; T reg67=reg62+reg63; T reg68=elem.pos(6)[1]*reg50; reg15=reg15-reg65; reg57=reg49+reg57;
    reg49=reg31*elem.pos(3)[0]; reg24=reg24-reg41; T reg69=elem.pos(6)[2]*reg50; T reg70=reg19*elem.pos(6)[1]; reg66=reg47+reg66;
    reg47=reg19*elem.pos(6)[2]; T reg71=elem.pos(6)[1]*reg61; reg64=reg28+reg64; reg51=reg51-reg60; reg28=elem.pos(4)[0]*reg38;
    reg58=reg32+reg58; reg32=elem.pos(4)[0]*reg40; reg52=reg59+reg52; reg59=elem.pos(6)[2]*reg61; T reg72=elem.pos(5)[0]*reg38;
    reg51=reg51-reg28; reg69=reg24+reg69; reg24=elem.pos(7)[2]*reg40; T reg73=elem.pos(4)[0]*reg3; T reg74=reg49+reg67;
    T reg75=reg31*elem.pos(7)[2]; reg47=reg66+reg47; reg66=elem.pos(5)[0]*reg50; reg52=reg52-reg32; T reg76=reg31*elem.pos(7)[1];
    T reg77=elem.pos(7)[2]*reg61; reg59=reg58+reg59; reg70=reg57+reg70; reg57=elem.pos(7)[1]*reg61; reg71=reg64+reg71;
    reg68=reg15+reg68; reg15=elem.pos(7)[1]*reg40; reg58=1+(*f.m).poisson_ratio; reg76=reg70+reg76; reg64=reg5*elem.pos(5)[0];
    reg73=reg73-reg74; reg75=reg47+reg75; reg47=elem.pos(6)[0]*reg50; reg24=reg69+reg24; reg59=reg59-reg77;
    reg71=reg71-reg57; reg52=reg52-reg66; reg15=reg68+reg15; reg68=elem.pos(6)[0]*reg61; reg72=reg51+reg72;
    reg51=reg59*reg76; reg69=elem.pos(7)[0]*reg40; reg47=reg52+reg47; reg52=elem.pos(7)[0]*reg61; reg68=reg72+reg68;
    reg70=reg15*reg75; reg58=reg58/(*f.m).elastic_modulus; reg72=reg24*reg76; T reg78=reg19*elem.pos(6)[0]; reg64=reg73+reg64;
    reg73=reg71*reg75; reg72=reg70-reg72; reg51=reg73-reg51; reg70=reg71*reg24; reg73=reg59*reg15;
    T reg79=pow(reg58,2); reg69=reg47+reg69; reg78=reg64+reg78; reg68=reg68-reg52; reg47=reg31*elem.pos(7)[0];
    reg58=reg58*reg79; reg73=reg70-reg73; reg64=reg68*reg72; reg70=reg69*reg51; reg47=reg78+reg47;
    reg78=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg80=1.0/(*f.m).elastic_modulus; reg70=reg64-reg70; reg64=reg47*reg73; T reg81=reg69*reg75;
    T reg82=reg24*reg47; reg75=reg68*reg75; T reg83=reg59*reg47; T reg84=reg80*reg79; reg79=reg78*reg79;
    T reg85=reg80*reg58; reg58=reg78*reg58; T reg86=reg68*reg76; T reg87=reg78*reg85; T reg88=reg15*reg47;
    reg83=reg75-reg83; reg82=reg81-reg82; reg76=reg69*reg76; reg75=reg78*reg58; reg47=reg71*reg47;
    reg85=reg80*reg85; reg24=reg68*reg24; reg81=reg78*reg84; reg59=reg59*reg69; reg84=reg80*reg84;
    reg64=reg70+reg64; reg70=reg78*reg79; reg47=reg86-reg47; reg69=reg71*reg69; reg15=reg68*reg15;
    reg83=reg83/reg64; reg59=reg24-reg59; reg51=reg51/reg64; reg88=reg76-reg88; reg85=reg85-reg75;
    reg81=reg70+reg81; reg82=reg82/reg64; reg87=reg75+reg87; reg72=reg72/reg64; reg84=reg84-reg70;
    reg58=reg80*reg58; reg79=reg80*reg79; reg24=reg38*reg82; reg68=reg78*reg87; reg71=reg50*reg51;
    reg76=reg38*reg72; reg86=reg80*reg85; reg81=reg78*reg81; reg58=reg75+reg58; reg84=reg80*reg84;
    reg75=reg70+reg79; reg80=reg50*reg83; reg69=reg15-reg69; reg59=reg59/reg64; reg73=reg73/reg64;
    reg47=reg47/reg64; reg88=reg88/reg64; reg15=reg38*reg88; T reg89=reg7*reg83; T reg90=reg10*reg88;
    T reg91=reg40*reg51; T reg92=reg10*reg72; T reg93=reg61*reg88; T reg94=reg61*reg72; T reg95=reg61*reg82;
    T reg96=reg10*reg82; T reg97=reg7*reg47; T reg98=reg40*reg47; reg81=reg84-reg81; reg84=reg7*reg51;
    T reg99=reg4*reg72; reg75=reg78*reg75; T reg100=reg40*reg83; T reg101=reg50*reg47; T reg102=reg71+reg76;
    T reg103=reg80+reg24; T reg104=reg5*reg59; T reg105=reg6*reg51; T reg106=reg4*reg82; reg68=reg86-reg68;
    reg78=reg78*reg58; reg86=reg6*reg83; reg69=reg69/reg64; T reg107=reg5*reg73; T reg108=reg106+reg86;
    T reg109=reg86-reg96; T reg110=reg100+reg95; T reg111=reg31*reg69; T reg112=reg31*reg73; T reg113=reg84+reg92;
    T reg114=reg31*reg59; T reg115=reg89+reg96; T reg116=reg80-reg95; T reg117=reg94-reg71; T reg118=reg93-reg101;
    T reg119=reg19*reg59; T reg120=reg4*reg88; T reg121=reg3*reg69; T reg122=reg3*reg73; T reg123=reg84-reg99;
    T reg124=reg98-reg15; T reg125=reg5*reg69; T reg126=reg6*reg47; T reg127=reg103+reg104; T reg128=reg3*reg59;
    T reg129=reg106-reg89; reg102=reg102+reg107; T reg130=reg101+reg15; T reg131=reg24-reg100; T reg132=reg91-reg76;
    T reg133=reg19*reg73; reg78=reg68-reg78; reg68=reg97+reg90; T reg134=reg91+reg94; T reg135=reg98+reg93;
    T reg136=reg19*reg69; reg75=reg81-reg75; reg81=reg99+reg105; T reg137=reg92-reg105; reg110=reg110-reg114;
    reg81=reg81-reg107; reg124=reg124+reg121; T reg138=reg112-reg134; T reg139=reg111-reg135; reg123=reg123-reg122;
    T reg140=reg97-reg120; T reg141=0.5*reg127; T reg142=reg90-reg126; reg129=reg129+reg128; T reg143=reg126+reg120;
    T reg144=(*f.m).deltaT*(*f.m).alpha; reg130=reg125+reg130; reg75=reg75/reg78; reg87=reg87/reg78; reg131=reg131-reg128;
    reg85=reg85/reg78; reg132=reg132+reg122; reg78=reg58/reg78; reg115=reg115+reg114; reg137=reg137-reg133;
    reg58=0.5*reg102; T reg145=reg111+reg68; T reg146=reg113+reg112; T reg147=reg104-reg108; reg116=reg116-reg119;
    reg117=reg117+reg133; reg109=reg109+reg119; reg118=reg118+reg136; T reg148=0.5*reg131; T reg149=reg75*reg58;
    T reg150=0.5*reg118; T reg151=0.5*reg124; T reg152=0.5*reg132; T reg153=0.5*reg129; T reg154=0.5*reg123;
    T reg155=0.5*reg110; T reg156=0.5*reg130; reg142=reg142-reg136; T reg157=0.5*reg137; T reg158=0.5*reg116;
    T reg159=0.5*reg147; reg143=reg143-reg125; T reg160=0.5*reg109; T reg161=0.5*reg145; T reg162=0.5*reg146;
    T reg163=0.5*reg117; T reg164=reg78*reg144; T reg165=reg87*reg144; T reg166=reg85*reg144; T reg167=0.5*reg81;
    reg140=reg140-reg121; T reg168=0.5*reg115; T reg169=reg75*reg141; T reg170=0.5*reg138; T reg171=0.5*reg139;
    T reg172=reg75*reg167; T reg173=reg75*reg159; T reg174=reg85*reg102; T reg175=reg75*reg163; T reg176=reg75*reg170;
    T reg177=reg75*reg158; T reg178=reg75*reg150; T reg179=reg75*reg155; T reg180=reg75*reg171; T reg181=2*reg169;
    T reg182=reg75*reg168; T reg183=reg75*reg161; reg149=2*reg149; T reg184=reg75*reg157; T reg185=reg75*reg160;
    T reg186=reg75*reg151; T reg187=reg85*reg127; T reg188=reg85*reg130; T reg189=reg75*reg156; T reg190=0.5*reg140;
    T reg191=reg75*reg162; T reg192=reg75*reg152; T reg193=reg75*reg148; T reg194=reg165+reg164; T reg195=0.5*reg142;
    T reg196=reg75*reg154; T reg197=reg166+reg165; T reg198=0.5*reg143; T reg199=reg75*reg153; T reg200=reg87*reg127;
    T reg201=reg162*reg149; T reg202=reg78*reg118; T reg203=reg115*reg187; T reg204=reg145*reg188; T reg205=reg85*reg139;
    T reg206=reg85*reg142; T reg207=reg197+reg164; T reg208=2*reg183; T reg209=reg85*reg132; T reg210=reg87*reg138;
    T reg211=reg85*reg117; T reg212=reg85*reg115; T reg213=reg75*reg195; T reg214=reg78*reg145; reg175=2*reg175;
    T reg215=reg85*reg140; T reg216=reg87*reg115; T reg217=reg85*reg110; T reg218=2*reg191; reg173=2*reg173;
    T reg219=reg87*reg146; reg176=2*reg176; T reg220=reg81*reg85; reg186=2*reg186; T reg221=reg85*reg146;
    T reg222=reg85*reg138; T reg223=reg143*reg85; reg179=2*reg179; reg182=2*reg182; reg172=2*reg172;
    reg180=2*reg180; T reg224=reg85*reg123; reg193=2*reg193; T reg225=reg85*reg116; T reg226=reg85*reg124;
    T reg227=reg85*reg109; T reg228=reg85*reg118; T reg229=reg87*reg117; T reg230=reg87*reg102; T reg231=reg146*reg174;
    T reg232=reg168*reg181; T reg233=reg78*reg130; T reg234=reg78*reg139; T reg235=reg87*reg132; T reg236=reg198*reg75;
    T reg237=reg31*reg1; T reg238=reg85*reg131; T reg239=reg85*reg129; T reg240=reg5*var_inter[2]; reg189=2*reg189;
    reg199=2*reg199; reg177=2*reg177; T reg241=reg166+reg194; T reg242=reg85*reg147; reg184=2*reg184;
    reg192=2*reg192; T reg243=reg85*reg145; T reg244=reg78*reg124; T reg245=reg75*reg190; T reg246=reg85*reg137;
    reg178=2*reg178; reg185=2*reg185; reg196=2*reg196; T reg247=reg31*var_inter[2]; T reg248=reg157*reg192;
    T reg249=reg109*reg238; T reg250=reg132*reg174; T reg251=reg5*reg1; T reg252=reg148*reg181; T reg253=reg182*reg159;
    T reg254=reg3*reg1; T reg255=reg154*reg218; T reg256=reg81*reg221; T reg257=reg132*reg211; T reg258=reg109*reg187;
    T reg259=reg148*reg177; T reg260=reg129*reg212; T reg261=reg157*reg149; T reg262=reg227*reg109; T reg263=reg157*reg184;
    T reg264=reg152*reg175; T reg265=reg131*reg225; T reg266=reg129*reg187; T reg267=reg222*reg137; T reg268=reg154*reg149;
    T reg269=reg179*reg160; T reg270=reg81*reg174; T reg271=reg152*reg176; T reg272=reg131*reg217; T reg273=reg181*reg159;
    T reg274=reg124*reg226; T reg275=reg211*reg137; T reg276=reg124*reg188; T reg277=reg129*reg238; T reg278=reg154*reg192;
    T reg279=reg177*reg160; T reg280=reg124*reg228; T reg281=reg190*reg218; T reg282=reg81*reg214; T reg283=reg198*reg218;
    T reg284=reg148*reg179; T reg285=reg132*reg222; T reg286=reg109*reg212; T reg287=reg81*reg209; T reg288=reg123*reg221;
    T reg289=reg153*reg182; T reg290=reg1*reg19; T reg291=reg152*reg192; T reg292=reg131*reg238; T reg293=reg193*reg159;
    T reg294=reg157*reg218; T reg295=reg123*reg174; T reg296=reg153*reg181; T reg297=reg123*reg220; T reg298=reg152*reg149;
    T reg299=reg131*reg187; T reg300=reg168*reg179; T reg301=reg146*reg222; T reg302=reg140*reg243; T reg303=reg146*reg234;
    T reg304=reg176*reg161; T reg305=reg78*reg115; T reg306=reg115*reg212; T reg307=reg162*reg218; T reg308=reg78*reg110;
    T reg309=reg115*reg214; T reg310=reg161*reg182; T reg311=reg154*reg208; T reg312=reg140*reg219; T reg313=reg142*reg188;
    T reg314=reg115*reg238; T reg315=reg140*reg206; T reg316=reg192*reg162; T reg317=reg140*reg205; T reg318=reg78*reg109;
    T reg319=reg203+reg201; T reg320=reg142*reg226; T reg321=reg168*reg218; T reg322=reg146*reg216; T reg323=reg168*reg193;
    T reg324=reg146*reg209; T reg325=reg140*reg188; T reg326=reg78*reg116; T reg327=reg146*reg221; T reg328=reg142*reg205;
    T reg329=reg78*reg127; T reg330=reg146*reg244; T reg331=reg192*reg161; reg231=reg232+reg231; T reg332=reg161*reg189;
    T reg333=reg140*reg228; T reg334=reg140*reg226; T reg335=reg146*reg233; T reg336=reg161*reg149; T reg337=reg78*reg131;
    T reg338=reg146*reg202; T reg339=reg175*reg161; T reg340=reg142*reg228; T reg341=reg81*reg246; T reg342=reg185*reg159;
    reg204=reg201+reg204; reg201=reg109*reg217; T reg343=reg145*reg229; T reg344=reg178*reg162; T reg345=reg157*reg176;
    T reg346=reg145*reg228; T reg347=reg145*reg210; T reg348=reg129*reg225; T reg349=reg154*reg175; T reg350=reg180*reg162;
    T reg351=reg109*reg225; T reg352=reg145*reg205; T reg353=reg212*reg147; T reg354=reg167*reg218; T reg355=reg157*reg175;
    T reg356=reg148*reg193; T reg357=reg132*reg209; T reg358=var_inter[2]*reg3; T reg359=var_inter[2]*reg19; T reg360=reg81*reg220;
    T reg361=reg173*reg159; T reg362=reg140*reg223; T reg363=reg115*reg225; T reg364=reg175*reg162; T reg365=reg78*reg147;
    T reg366=reg142*reg243; T reg367=reg115*reg217; T reg368=reg176*reg162; T reg369=reg157*reg208; T reg370=reg145*reg243;
    T reg371=reg142*reg219; T reg372=reg140*reg215; T reg373=reg145*reg235; T reg374=reg186*reg162; T reg375=reg142*reg206;
    T reg376=reg145*reg226; T reg377=reg145*reg230; T reg378=reg162*reg189; T reg379=reg129*reg217; T reg380=reg154*reg176;
    T reg381=reg154*reg184; T reg382=reg227*reg129; T reg383=reg143*reg188; T reg384=reg176*reg58; T reg385=reg127*reg217;
    T reg386=reg87*reg131; T reg387=reg130*reg188; T reg388=reg167*reg184; T reg389=reg123*reg209; T reg390=reg153*reg193;
    T reg391=reg227*reg147; T reg392=reg143*reg226; T reg393=reg130*reg228; T reg394=reg123*reg214; T reg395=reg240*(*f.m).f_vol[1];
    T reg396=reg143*reg243; T reg397=reg130*reg205; T reg398=reg87*reg110; T reg399=reg158*reg177; T reg400=reg246*reg137;
    T reg401=reg145*reg241; T reg402=reg185*reg160; T reg403=reg58*reg149; T reg404=reg127*reg187; T reg405=reg123*reg211;
    T reg406=reg153*reg177; T reg407=reg156*reg181; T reg408=reg143*reg205; T reg409=reg127*reg233; T reg410=reg168*reg182;
    T reg411=reg167*reg172; T reg412=reg58*reg175; T reg413=reg127*reg225; T reg414=reg242*reg147; T reg415=reg237*(*f.m).f_vol[2];
    T reg416=reg146*reg207; T reg417=reg143*reg228; T reg418=reg87*reg137; T reg419=reg167*reg176; T reg420=reg118*reg228;
    T reg421=reg123*reg246; T reg422=reg153*reg185; reg213=2*reg213; T reg423=reg87*reg129; T reg424=reg139*reg205;
    T reg425=reg147*reg225; T reg426=reg167*reg192; T reg427=reg238*reg147; T reg428=reg167*reg175; T reg429=reg78*reg140;
    T reg430=reg110*reg217; T reg431=reg170*reg176; T reg432=reg187*reg147; T reg433=reg167*reg149; T reg434=reg138*reg222;
    T reg435=reg118*reg205; T reg436=reg155*reg179; T reg437=reg117*reg211; T reg438=reg167*reg208; T reg439=reg123*reg222;
    T reg440=reg153*reg179; T reg441=reg199*reg153; T reg442=reg143*reg219; T reg443=reg168*reg177; T reg444=reg146*reg211;
    T reg445=reg179*reg158; T reg446=reg222*reg117; T reg447=reg143*reg206; T reg448=reg224*reg123; T reg449=reg163*reg175;
    T reg450=reg116*reg225; T reg451=reg78*reg142; T reg452=reg143*reg223; T reg453=reg163*reg176; T reg454=reg87*reg109;
    T reg455=reg116*reg217; reg217=reg147*reg217; T reg456=reg177*reg159; reg205=reg124*reg205; T reg457=reg160*reg181;
    T reg458=reg160*reg182; T reg459=reg137*reg221; T reg460=reg87*reg147; T reg461=reg102*reg174; T reg462=reg179*reg159;
    T reg463=reg81*reg222; T reg464=reg141*reg177; T reg465=reg143*reg78; T reg466=reg102*reg211; T reg467=reg195*reg218;
    T reg468=reg237*(*f.m).f_vol[0]; T reg469=reg214*reg137; T reg470=reg81*reg87; reg245=2*reg245; T reg471=reg141*reg181;
    T reg472=reg193*reg160; T reg473=reg242*reg129; T reg474=reg154*reg172; T reg475=reg127*reg207; T reg476=reg102*reg200;
    T reg477=reg209*reg137; T reg478=reg141*reg149; reg222=reg102*reg222; T reg479=reg196*reg154; T reg480=reg129*reg239;
    T reg481=reg141*reg179; reg236=2*reg236; T reg482=reg137*reg174; T reg483=reg153*reg173; T reg484=reg81*reg211;
    T reg485=reg87*reg116; T reg486=reg213*reg159; T reg487=reg167*reg213; T reg488=reg143*reg418; reg482=reg482-reg457;
    T reg489=reg160*reg189; T reg490=reg142*reg329; T reg491=reg160*reg149; reg452=reg411+reg452; T reg492=reg195*reg192;
    T reg493=reg234*reg147; T reg494=reg137*reg200; T reg495=reg198*reg179; reg217=reg419+reg217; T reg496=reg410+reg327;
    T reg497=reg157*reg189; T reg498=reg386*reg137; T reg499=reg143*reg318; T reg500=reg192*reg160; T reg501=reg142*reg230;
    T reg502=reg208*reg161; reg320=reg248+reg320; T reg503=reg398*reg117; reg267=reg269+reg267; T reg504=reg176*reg160;
    T reg505=reg398*reg137; T reg506=reg142*reg337; reg447=reg388+reg447; T reg507=reg230*reg147; T reg508=reg180*reg160;
    T reg509=reg147*reg229; reg425=reg428+reg425; T reg510=reg198*reg177; T reg511=reg195*reg178; T reg512=reg202*reg147;
    T reg513=reg157*reg180; T reg514=reg142*reg308; T reg515=reg142*reg210; reg340=reg355+reg340; reg275=reg279+reg275;
    T reg516=reg233*reg137; T reg517=reg167*reg179; T reg518=reg175*reg160; T reg519=reg195*reg149; T reg520=reg485*reg137;
    T reg521=reg195*reg175; T reg522=reg167*reg177; T reg523=reg233*reg147; T reg524=reg142*reg326; T reg525=reg178*reg160;
    T reg526=reg157*reg178; T reg527=reg147*reg210; T reg528=reg195*reg189; T reg529=reg142*reg229; T reg530=reg202*reg137;
    T reg531=reg195*reg180; reg328=reg345+reg328; reg313=reg261+reg313; T reg532=reg198*reg181; T reg533=reg433-reg432;
    T reg534=reg137*reg244; T reg535=reg195*reg181; T reg536=reg233*reg109; T reg537=reg157*reg177; T reg538=reg195*reg208;
    T reg539=reg180*reg159; T reg540=reg167*reg180; T reg541=reg143*reg210; reg417=reg428+reg417; reg428=reg109*reg229;
    T reg542=reg458-reg459; T reg543=reg143*reg326; T reg544=reg178*reg159; T reg545=reg167*reg178; T reg546=reg143*reg229;
    reg351=reg355+reg351; reg355=reg195*reg177; reg383=reg433+reg383; reg433=reg160*reg218; T reg547=reg157*reg193;
    T reg548=reg109*reg235; reg400=reg402+reg400; T reg549=reg184*reg160; T reg550=reg195*reg213; reg249=reg248+reg249;
    reg248=reg195*reg193; T reg551=reg454*reg137; T reg552=reg195*reg184; reg408=reg419+reg408; reg419=reg109*reg244;
    T reg553=reg157*reg181; T reg554=reg109*reg230; T reg555=reg214*reg109; T reg556=reg195*reg182; reg286=reg286-reg294;
    T reg557=reg143*reg308; T reg558=reg451*reg137; reg261=reg261-reg258; T reg559=reg195*reg186; reg375=reg263+reg375;
    T reg560=reg167*reg186; T reg561=reg143*reg235; T reg562=reg354+reg396; T reg563=reg143*reg305; reg262=reg263+reg262;
    reg263=reg371+reg369; T reg564=reg208*reg160; T reg565=reg208*reg159; T reg566=reg142*reg305; T reg567=reg294+reg366;
    T reg568=reg442+reg438; T reg569=reg142*reg235; T reg570=reg157*reg186; T reg571=reg137*reg234; T reg572=reg195*reg176;
    reg477=reg472+reg477; T reg573=reg186*reg160; T reg574=reg202*reg109; T reg575=reg109*reg219; T reg576=reg157*reg179;
    T reg577=reg109*reg210; T reg578=reg143*reg329; T reg579=reg189*reg159; T reg580=reg167*reg189; T reg581=reg157*reg182;
    T reg582=reg137*reg216; T reg583=reg451*reg109; T reg584=reg143*reg230; reg201=reg345+reg201; reg392=reg426+reg392;
    reg345=reg143*reg337; T reg585=reg195*reg185; T reg586=reg467+reg469; T reg587=reg195*reg179; T reg588=reg186*reg159;
    T reg589=reg109*reg234; T reg590=reg102*reg485; T reg591=reg141*reg175; T reg592=reg102*reg202; T reg593=reg132*reg207;
    T reg594=reg156*reg175; reg222=reg222-reg481; T reg595=reg401-reg415; T reg596=reg156*reg180; T reg597=reg102*reg398;
    T reg598=reg141*reg176; T reg599=reg102*reg234; T reg600=reg156*reg176; T reg601=reg403+reg404; T reg602=reg115*reg207;
    T reg603=reg416-reg468; T reg604=reg407+reg409; T reg605=reg58*reg177; T reg606=reg127*reg229; reg413=reg412-reg413;
    T reg607=reg156*reg177; T reg608=reg127*reg202; T reg609=reg179*reg58; T reg610=reg127*reg210; T reg611=reg142*reg241;
    T reg612=reg109*reg207; T reg613=reg137*reg207; T reg614=reg143*reg241; T reg615=reg124*reg329; reg276=reg298+reg276;
    T reg616=reg124*reg229; T reg617=reg152*reg178; T reg618=reg148*reg178; T reg619=reg124*reg326; reg280=reg264+reg280;
    T reg620=reg124*reg210; T reg621=reg152*reg180; T reg622=reg148*reg180; T reg623=reg124*reg308; reg205=reg271+reg205;
    T reg624=reg118*reg241; T reg625=reg116*reg207; T reg626=reg117*reg207; reg461=reg461+reg471; T reg627=reg156*reg189;
    T reg628=reg130*reg241; T reg629=reg475-reg395; T reg630=reg102*reg207; T reg631=reg124*reg241; T reg632=reg131*reg207;
    reg478=reg476+reg478; T reg633=reg102*reg233; T reg634=reg156*reg149; reg466=reg466-reg464; T reg635=reg156*reg178;
    T reg636=reg176*reg158; T reg637=reg117*reg234; T reg638=reg176*reg150; reg450=reg449+reg450; T reg639=reg150*reg177;
    T reg640=reg202*reg116; T reg641=reg163*reg179; T reg642=reg116*reg210; reg455=reg453+reg455; T reg643=reg179*reg150;
    T reg644=reg116*reg234; reg420=reg449+reg420; reg424=reg431+reg424; reg449=reg118*reg210; T reg645=reg163*reg180;
    T reg646=reg180*reg158; T reg647=reg118*reg308; T reg648=reg110*reg234; T reg649=reg171*reg179; reg430=reg431+reg430;
    reg431=reg171*reg176; T reg650=reg138*reg234; T reg651=reg138*reg398; T reg652=reg155*reg176; T reg653=reg171*reg180;
    reg434=reg436+reg434; reg435=reg453+reg435; reg453=reg147*reg207; T reg654=reg81*reg207; T reg655=reg140*reg241;
    reg385=reg384-reg385; T reg656=reg156*reg179; T reg657=reg127*reg234; reg387=reg403+reg387; reg403=reg130*reg229;
    T reg658=reg129*reg207; T reg659=reg58*reg178; T reg660=reg130*reg326; T reg661=reg141*reg178; reg393=reg412+reg393;
    reg412=reg123*reg207; T reg662=reg130*reg210; T reg663=reg180*reg58; T reg664=reg130*reg308; T reg665=reg141*reg180;
    reg397=reg384+reg397; reg437=reg399+reg437; reg384=reg150*reg178; T reg666=reg158*reg175; T reg667=reg117*reg485;
    T reg668=reg202*reg117; T reg669=reg150*reg175; reg446=reg445+reg446; T reg670=reg180*reg150; T reg671=reg115*reg230;
    T reg672=reg162*reg181; T reg673=reg332+reg319; T reg674=reg115*reg233; T reg675=reg161*reg181; T reg676=reg115*reg229;
    T reg677=reg177*reg162; reg363=reg363-reg364; T reg678=reg115*reg202; T reg679=reg177*reg161; T reg680=reg115*reg210;
    T reg681=reg179*reg162; reg367=reg367-reg368; T reg682=reg115*reg234; T reg683=reg179*reg161; T reg684=reg307+reg370;
    reg374=reg373+reg374; reg373=reg168*reg186; T reg685=reg145*reg337; reg376=reg316+reg376; reg378=reg377+reg378;
    reg377=reg168*reg189; T reg686=reg145*reg329; reg204=reg232+reg204; reg344=reg343+reg344; reg343=reg168*reg178;
    T reg687=reg145*reg326; reg322=reg321+reg322; T reg688=reg146*reg214; T reg689=reg161*reg218; reg324=reg323-reg324;
    T reg690=reg186*reg161; T reg691=reg168*reg192; T reg692=reg146*reg386; reg331=reg330+reg331; reg332=reg231+reg332;
    T reg693=reg168*reg149; T reg694=reg146*reg200; reg336=reg335+reg336; T reg695=reg168*reg175; T reg696=reg146*reg485;
    reg339=reg338+reg339; reg301=reg300-reg301; T reg697=reg180*reg161; T reg698=reg168*reg176; T reg699=reg146*reg398;
    reg304=reg303+reg304; reg306=reg306+reg307; reg310=reg309+reg310; T reg700=reg115*reg235; T reg701=reg193*reg162;
    reg316=reg314-reg316; reg314=reg115*reg244; T reg702=reg193*reg161; T reg703=reg132*reg398; T reg704=reg151*reg176;
    T reg705=reg132*reg234; reg292=reg291+reg292; T reg706=reg151*reg193; T reg707=reg131*reg244; T reg708=reg152*reg181;
    T reg709=reg110*reg207; T reg710=reg131*reg230; reg298=reg298-reg299; T reg711=reg151*reg181; T reg712=reg131*reg233;
    T reg713=reg152*reg177; T reg714=reg131*reg229; reg265=reg264+reg265; reg264=reg151*reg177; T reg715=reg131*reg202;
    T reg716=reg138*reg207; T reg717=reg152*reg179; T reg718=reg131*reg210; reg272=reg271+reg272; reg271=reg151*reg179;
    T reg719=reg131*reg234; reg274=reg291+reg274; reg291=reg124*reg230; T reg720=reg152*reg189; T reg721=reg148*reg189;
    reg346=reg364+reg346; reg350=reg347+reg350; reg347=reg168*reg180; reg364=reg145*reg308; reg352=reg368+reg352;
    reg368=reg151*reg186; T reg722=reg139*reg241; reg357=reg356+reg357; T reg723=reg148*reg192; T reg724=reg132*reg386;
    T reg725=reg151*reg192; T reg726=reg132*reg244; T reg727=reg151*reg189; reg250=reg250-reg252; T reg728=reg148*reg149;
    T reg729=reg132*reg200; T reg730=reg151*reg149; T reg731=reg132*reg233; T reg732=reg151*reg178; reg257=reg259+reg257;
    T reg733=reg148*reg175; T reg734=reg132*reg485; T reg735=reg151*reg175; T reg736=reg132*reg202; T reg737=reg151*reg180;
    reg285=reg284+reg285; T reg738=reg148*reg176; T reg739=reg254*(*f.m).f_vol[1]; T reg740=reg81*reg485; T reg741=reg254*(*f.m).f_vol[2];
    T reg742=reg251*(*f.m).f_vol[0]; T reg743=reg429*reg129; T reg744=reg255+reg302; T reg745=reg198*reg178; reg484=reg484+reg456;
    T reg746=reg251*(*f.m).f_vol[1]; T reg747=reg251*(*f.m).f_vol[2]; T reg748=reg199*reg190; T reg749=reg198*reg149; T reg750=reg140*reg235;
    T reg751=reg190*reg213; reg421=reg422+reg421; T reg752=reg81*reg233; T reg753=reg129*reg451; T reg754=reg81*reg244;
    T reg755=reg198*reg192; T reg756=reg154*reg189; T reg757=reg190*reg184; T reg758=reg140*reg230; T reg759=reg123*reg234;
    reg270=reg270-reg273; reg334=reg278+reg334; T reg760=reg198*reg189; T reg761=reg123*reg454; T reg762=reg153*reg184;
    T reg763=reg153*reg186; reg480=reg479+reg480; T reg764=reg81*reg200; T reg765=reg149*reg159; T reg766=reg154*reg186;
    T reg767=reg418*reg140; reg414=reg411+reg414; reg411=reg154*reg213; T reg768=reg423*reg123; T reg769=reg129*reg230;
    T reg770=reg154*reg181; T reg771=reg153*reg213; T reg772=reg198*reg176; T reg773=reg81*reg234; T reg774=reg140*reg318;
    T reg775=reg196*reg190; T reg776=reg429*reg123; T reg777=reg176*reg159; T reg778=reg129*reg244; T reg779=reg81*reg398;
    reg315=reg381+reg315; T reg780=(*f.m).f_vol[0]*reg254; T reg781=reg140*reg305; T reg782=reg175*reg159; T reg783=reg153*reg208;
    T reg784=reg240*(*f.m).f_vol[0]; T reg785=reg290*(*f.m).f_vol[1]; reg277=reg278+reg277; reg278=reg81*reg202; T reg786=reg198*reg175;
    T reg787=(*f.m).f_vol[0]*reg290; T reg788=reg358*(*f.m).f_vol[2]; T reg789=reg312+reg311; T reg790=reg358*(*f.m).f_vol[1]; T reg791=reg358*(*f.m).f_vol[0];
    reg463=reg463+reg462; T reg792=reg198*reg180; T reg793=reg190*reg193; T reg794=reg123*reg485; T reg795=reg198*reg184;
    T reg796=reg140*reg210; T reg797=reg190*reg176; T reg798=reg81*reg451; T reg799=reg153*reg175; T reg800=reg154*reg180;
    T reg801=reg240*(*f.m).f_vol[2]; T reg802=reg394+reg281; T reg803=reg184*reg159; T reg804=reg81*reg454; reg405=reg406+reg405;
    T reg805=reg153*reg180; T reg806=reg140*reg308; T reg807=reg190*reg178; T reg808=reg190*reg186; T reg809=reg237*(*f.m).f_vol[1];
    T reg810=reg359*(*f.m).f_vol[0]; T reg811=reg123*reg200; T reg812=reg81*reg460; T reg813=reg172*reg159; T reg814=reg198*reg236;
    reg360=reg360+reg361; T reg815=reg123*reg386; T reg816=reg153*reg192; T reg817=reg190*reg149; T reg818=reg81*reg465;
    T reg819=reg198*reg172; T reg820=reg123*reg233; reg317=reg380+reg317; reg389=reg390+reg389; reg341=reg341+reg342;
    T reg821=reg198*reg213; T reg822=reg140*reg229; T reg823=reg190*reg172; T reg824=reg123*reg202; reg287=reg287+reg293;
    T reg825=reg190*reg175; T reg826=reg123*reg465; T reg827=reg198*reg186; reg325=reg268+reg325; T reg828=reg359*(*f.m).f_vol[2];
    T reg829=reg140*reg329; T reg830=reg359*(*f.m).f_vol[1]; T reg831=reg81*reg386; T reg832=reg192*reg159; T reg833=reg247*(*f.m).f_vol[1];
    T reg834=reg153*reg189; T reg835=reg247*(*f.m).f_vol[2]; T reg836=reg123*reg451; T reg837=reg253-reg256; T reg838=reg190*reg236;
    reg398=reg123*reg398; reg333=reg349+reg333; T reg839=reg198*reg208; reg297=reg297+reg483; reg176=reg153*reg176;
    reg439=reg440+reg439; T reg840=reg153*reg172; T reg841=reg81*reg216; T reg842=reg218*reg159; T reg843=reg140*reg326;
    T reg844=reg153*reg178; T reg845=reg190*reg180; T reg846=reg123*reg460; T reg847=reg154*reg178; T reg848=reg282+reg283;
    T reg849=reg190*reg179; reg234=reg129*reg234; T reg850=reg154*reg185; T reg851=reg167*reg193; T reg852=reg470*reg129;
    T reg853=reg129*reg219; T reg854=reg451*reg147; reg295=reg295-reg296; T reg855=reg190*reg189; T reg856=reg418*reg129;
    T reg857=reg235*reg147; T reg858=reg198*reg185; T reg859=reg129*reg229; T reg860=reg154*reg173; reg372=reg479+reg372;
    reg353=reg353-reg354; reg479=reg123*reg244; reg427=reg426+reg427; reg391=reg388+reg391; reg388=reg470*reg140;
    reg426=reg154*reg236; T reg861=reg140*reg337; T reg862=reg190*reg192; T reg863=reg153*reg218; T reg864=reg154*reg177;
    T reg865=reg154*reg179; T reg866=reg129*reg210; T reg867=reg129*reg202; T reg868=reg123*reg216; T reg869=reg178*reg161;
    reg473=reg474+reg473; T reg870=reg190*reg177; reg444=reg443-reg444; T reg871=reg289-reg288; T reg872=reg198*reg182;
    T reg873=reg190*reg173; T reg874=reg129*reg465; reg379=reg380+reg379; reg380=reg190*reg208; T reg875=reg153*reg149;
    T reg876=reg219*reg147; T reg877=reg167*reg182; T reg878=reg214*reg147; T reg879=reg190*reg185; reg348=reg349+reg348;
    reg349=reg154*reg182; T reg880=reg153*reg236; T reg881=reg140*reg365; T reg882=reg190*reg181; T reg883=reg198*reg193;
    T reg884=reg247*(*f.m).f_vol[0]; T reg885=reg190*reg245; reg268=reg268-reg266; T reg886=reg465*reg147; T reg887=reg154*reg193;
    reg362=reg474+reg362; reg474=reg244*reg147; T reg888=reg198*reg173; T reg889=reg290*(*f.m).f_vol[2]; T reg890=reg190*reg182;
    T reg891=reg129*reg214; reg448=reg441+reg448; T reg892=reg167*reg181; T reg893=reg196*reg153; reg382=reg381+reg382;
    reg381=reg129*reg233; T reg894=reg418*reg147; T reg895=reg167*reg185; T reg896=reg129*reg235; reg260=reg260-reg255;
    reg313=reg313-reg457; reg594=reg592+reg594; reg328=reg269+reg328; reg724=reg723+reg724; reg346=reg443-reg346;
    reg348=reg807+reg348; reg875=reg875-reg811; reg843=reg844+reg843; reg736=reg735+reg736; reg687=reg343-reg687;
    reg570=reg569+reg570; reg250=reg727+reg250; reg607=reg607-reg608; reg867=reg870+reg867; reg466=reg466+reg635;
    reg269=reg64*reg322; reg489=reg489-reg490; reg890=reg890-reg891; reg726=reg725+reg726; reg343=reg64*reg344;
    reg360=reg360+reg814; reg847=reg822+reg847; reg846=reg840+reg846; reg731=reg730+reg731; reg381=reg381-reg882;
    reg868=reg868-reg863; reg496=reg496+reg502; reg591=reg590-reg591; reg806=reg805+reg806; reg443=reg64*reg604;
    reg497=reg501+reg497; reg364=reg347-reg364; reg524=reg525+reg524; reg349=reg349-reg853; reg800=reg796+reg800;
    reg728=reg728-reg729; reg820=reg817+reg820; reg357=reg368+reg357; reg260=reg260-reg380; reg859=reg864+reg859;
    reg606=reg605-reg606; reg598=reg597-reg598; reg320=reg472+reg320; reg340=reg279+reg340; reg317=reg440+reg317;
    reg297=reg838+reg297; reg601=reg627+reg601; reg514=reg508+reg514; reg413=reg635+reg413; reg352=reg300-reg352;
    reg734=reg733+reg734; reg506=reg573+reg506; reg257=reg732+reg257; reg279=reg64*reg350; reg887=reg896+reg887;
    reg405=reg807+reg405; reg222=reg222+reg596; reg526=reg529+reg526; reg513=reg515+reg513; reg333=reg406+reg333;
    reg600=reg599+reg600; reg794=reg799+reg794; reg363=reg363-reg869; reg277=reg808+reg277; reg300=reg64*reg789;
    reg721=reg721-reg615; reg306=reg502+reg306; reg679=reg678-reg679; reg714=reg713+reg714; reg681=reg680-reg681;
    reg712=reg712-reg711; reg276=reg276-reg252; reg347=reg64*reg304; reg426=reg388+reg426; reg367=reg367-reg697;
    reg781=reg781-reg783; reg617=reg616+reg617; reg743=reg748+reg743; reg298=reg727+reg298; reg699=reg698-reg699;
    reg479=reg862+reg479; reg683=reg682-reg683; reg619=reg618+reg619; reg372=reg441+reg372; reg388=reg64*reg673;
    reg769=reg769-reg770; reg411=reg767+reg411; reg718=reg717+reg718; reg674=reg674+reg675; reg272=reg737+reg272;
    reg671=reg671+reg672; reg715=reg264+reg715; reg719=reg271+reg719; reg702=reg314-reg702; reg362=reg483+reg362;
    reg774=reg771+reg774; reg677=reg676-reg677; reg778=reg793+reg778; reg316=reg316-reg690; reg274=reg356+reg274;
    reg268=reg268+reg855; reg265=reg732+reg265; reg701=reg700-reg701; reg315=reg422+reg315; reg720=reg291+reg720;
    reg264=reg64*reg310; reg881=reg880+reg881; reg271=reg64*reg332; reg879=reg753+reg879; reg376=reg323-reg376;
    reg379=reg845+reg379; reg705=reg704+reg705; reg756=reg758+reg756; reg627=reg461+reg627; reg291=reg64*reg331;
    reg314=reg64*reg378; reg323=reg64*reg478; reg871=reg871-reg380; reg866=reg865+reg866; reg692=reg691-reg692;
    reg377=reg377+reg686; reg834=reg834-reg829; reg703=reg738+reg703; reg690=reg324-reg690; reg826=reg823+reg826;
    reg634=reg633+reg634; reg325=reg325-reg296; reg324=reg64*reg204; reg285=reg737+reg285; reg356=reg688+reg689;
    reg410=reg410+reg684; reg710=reg710-reg708; reg697=reg301-reg697; reg289=reg289-reg744; reg280=reg259+reg280;
    reg234=reg849+reg234; reg259=reg64*reg339; reg301=reg64*reg374; reg766=reg750+reg766; reg480=reg885+reg480;
    reg621=reg620+reg621; reg696=reg695-reg696; reg707=reg706+reg707; reg295=reg855+reg295; reg623=reg622+reg623;
    reg685=reg373-reg685; reg373=reg64*reg336; reg763=reg861+reg763; reg759=reg797+reg759; reg205=reg284+reg205;
    reg292=reg368+reg292; reg693=reg693+reg694; reg334=reg390+reg334; reg558=reg552+reg558; reg431=reg650+reg431;
    reg776=reg775+reg776; reg551=reg549+reg551; reg777=reg779+reg777; reg400=reg550+reg400; reg430=reg653+reg430;
    reg648=reg649+reg648; reg408=reg462+reg408; reg772=reg773+reg772; reg768=reg893+reg768; reg557=reg539+reg557;
    reg424=reg436+reg424; reg540=reg541+reg540; reg414=reg814+reg414; reg448=reg885+reg448; reg417=reg456+reg417;
    reg284=reg780+reg412; reg543=reg544+reg543; reg886=reg888+reg886; reg368=reg739+reg658; reg545=reg546+reg545;
    reg382=reg382+reg751; reg390=reg741+reg655; reg383=reg383-reg273; reg579=reg579-reg578; reg894=reg895+reg894;
    reg270=reg270+reg760; reg642=reg641+reg642; reg520=reg518+reg520; reg455=reg670+reg455; reg275=reg511+reg275;
    reg765=reg765-reg764; reg421=reg751+reg421; reg644=reg643+reg644; reg516=reg519+reg516; reg491=reg491-reg494;
    reg749=reg752+reg749; reg420=reg399+reg420; reg482=reg528+reg482; reg645=reg449+reg645; reg534=reg492+reg534;
    reg484=reg484+reg745; reg647=reg646+reg647; reg498=reg500+reg498; reg477=reg559+reg477; reg782=reg740+reg782;
    reg435=reg445+reg435; reg399=reg64*reg586; reg786=reg278+reg786; reg582=reg582-reg433; reg653=reg434+reg653;
    reg542=reg542-reg538; reg463=reg463+reg792; reg651=reg652+reg651; reg278=reg788+reg631; reg452=reg361+reg452;
    reg353=reg353-reg839; reg361=reg784+reg630; reg493=reg495+reg493; reg629=reg64*reg629; reg217=reg792+reg217;
    reg872=reg872-reg878; reg852=reg860+reg852; reg406=reg801+reg628; reg527=reg517+reg527; reg422=reg810+reg626;
    reg512=reg510+reg512; reg857=reg851+reg857; reg425=reg745+reg425; reg434=reg830+reg625; reg436=reg828+reg624;
    reg509=reg522+reg509; reg427=reg827+reg427; reg440=reg884+reg716; reg523=reg523-reg532; reg441=reg833+reg709;
    reg533=reg760+reg533; reg474=reg883+reg474; reg445=reg835+reg722; reg507=reg507-reg892; reg449=reg742+reg654;
    reg456=reg746+reg453; reg580=reg584+reg580; reg461=reg747+reg614; reg392=reg293+reg392; reg391=reg821+reg391;
    reg856=reg850+reg856; reg293=reg787+reg613; reg345=reg588+reg345; reg462=reg785+reg612; reg560=reg561+reg560;
    reg854=reg858+reg854; reg253=reg253-reg562; reg472=reg889+reg611; reg603=reg64*reg603; reg563=reg563-reg565;
    reg874=reg873+reg874; reg483=reg64*reg568; reg877=reg877-reg876; reg492=reg809+reg602; reg447=reg342+reg447;
    reg595=reg64*reg595; reg342=reg791+reg593; reg499=reg486+reg499; reg869=reg444-reg869; reg487=reg488+reg487;
    reg444=reg790+reg632; reg838=reg473+reg838; reg663=reg662+reg663; reg536=reg536-reg535; reg795=reg798+reg795;
    reg261=reg528+reg261; reg398=reg176+reg398; reg665=reg664-reg665; reg554=reg554-reg553; reg837=reg837-reg839;
    reg481=reg397-reg481; reg419=reg248+reg419; reg439=reg845+reg439; reg249=reg559+reg249; reg841=reg841-reg842;
    reg437=reg437+reg384; reg548=reg547+reg548; reg667=reg666+reg667; reg556=reg556-reg555; reg176=reg64*reg848;
    reg824=reg825+reg824; reg669=reg668+reg669; reg286=reg286-reg538; reg581=reg581-reg575; reg827=reg287+reg827;
    reg458=reg458-reg567; reg610=reg609-reg610; reg566=reg566-reg564; reg815=reg816+reg815; reg248=reg64*reg263;
    reg813=reg812+reg813; reg385=reg596+reg385; reg375=reg402+reg375; reg656=reg656-reg657; reg589=reg587+reg589;
    reg819=reg818+reg819; reg389=reg808+reg389; reg201=reg531+reg201; reg387=reg471+reg387; reg577=reg576+reg577;
    reg821=reg341+reg821; reg659=reg403+reg659; reg574=reg355+reg574; reg287=reg64*reg802; reg661=reg660-reg661;
    reg351=reg511+reg351; reg803=reg804+reg803; reg464=reg393-reg464; reg428=reg537+reg428; reg530=reg521+reg530;
    reg755=reg754+reg755; reg761=reg762+reg761; reg670=reg446+reg670; reg836=reg757+reg836; reg640=reg639+reg640;
    reg636=reg503+reg636; reg450=reg384+reg450; reg638=reg637+reg638; reg267=reg531+reg267; reg262=reg550+reg262;
    reg571=reg572+reg571; reg583=reg585+reg583; reg505=reg504+reg505; reg832=reg831+reg832; reg480=reg64*reg480;
    reg659=reg64*reg659; reg285=reg64*reg285; reg341=reg64*reg444; reg601=reg64*reg601; reg644=reg64*reg644;
    reg387=reg64*reg387; reg355=reg64*reg278; reg820=reg64*reg820; reg736=reg64*reg736; reg295=reg64*reg295;
    reg292=reg64*reg292; reg603=ponderation*reg603; reg420=reg64*reg420; reg705=reg64*reg705; reg600=reg64*reg600;
    reg661=reg64*reg661; reg384=reg64*reg492; reg405=reg64*reg405; reg595=ponderation*reg595; reg703=reg64*reg703;
    reg450=reg64*reg450; reg280=reg64*reg280; reg393=ponderation*reg287; reg871=reg64*reg871; reg397=reg64*reg342;
    reg838=reg64*reg838; reg731=reg64*reg731; reg402=reg64*reg434; reg385=reg64*reg385; reg205=reg64*reg205;
    reg728=reg64*reg728; reg761=reg64*reg761; reg403=reg64*reg436; reg413=reg64*reg413; reg642=reg64*reg642;
    reg446=reg64*reg440; reg260=reg64*reg260; reg623=reg64*reg623; reg610=reg64*reg610; reg759=reg64*reg759;
    reg250=reg64*reg250; reg473=reg64*reg441; reg607=reg64*reg607; reg486=reg64*reg445; reg488=reg64*reg361;
    reg852=reg64*reg852; reg868=reg64*reg868; reg734=reg64*reg734; reg495=ponderation*reg443; reg879=reg64*reg879;
    reg421=reg64*reg421; reg629=ponderation*reg629; reg656=reg64*reg656; reg257=reg64*reg257; reg621=reg64*reg621;
    reg500=reg64*reg406; reg640=reg64*reg640; reg389=reg64*reg389; reg606=reg64*reg606; reg501=reg64*reg422;
    reg455=reg64*reg455; reg349=reg64*reg349; reg875=reg64*reg875; reg826=reg64*reg826; reg648=reg64*reg648;
    reg721=reg64*reg721; reg824=reg64*reg824; reg667=reg64*reg667; reg424=reg64*reg424; reg503=ponderation*reg323;
    reg718=reg64*reg718; reg448=reg64*reg448; reg769=reg64*reg769; reg437=reg64*reg437; reg715=reg64*reg715;
    reg591=reg64*reg591; reg435=reg64*reg435; reg504=reg64*reg284; reg636=reg64*reg636; reg276=reg64*reg276;
    reg743=reg64*reg743; reg720=reg64*reg720; reg277=reg64*reg277; reg651=reg64*reg651; reg776=reg64*reg776;
    reg274=reg64*reg274; reg634=reg64*reg634; reg431=reg64*reg431; reg794=reg64*reg794; reg669=reg64*reg669;
    reg670=reg64*reg670; reg719=reg64*reg719; reg846=reg64*reg846; reg653=reg64*reg653; reg430=reg64*reg430;
    reg768=reg64*reg768; reg778=reg64*reg778; reg466=reg64*reg466; reg272=reg64*reg272; reg638=reg64*reg638;
    reg508=reg64*reg461; reg479=reg64*reg479; reg710=reg64*reg710; reg222=reg64*reg222; reg398=reg64*reg398;
    reg663=reg64*reg663; reg510=reg64*reg293; reg627=reg64*reg627; reg645=reg64*reg645; reg707=reg64*reg707;
    reg836=reg64*reg836; reg511=reg64*reg462; reg598=reg64*reg598; reg464=reg64*reg464; reg515=reg64*reg472;
    reg874=reg64*reg874; reg619=reg64*reg619; reg265=reg64*reg265; reg439=reg64*reg439; reg517=reg64*reg368;
    reg382=reg64*reg382; reg268=reg64*reg268; reg481=reg64*reg481; reg594=reg64*reg594; reg297=reg64*reg297;
    reg518=reg64*reg390; reg714=reg64*reg714; reg815=reg64*reg815; reg712=reg64*reg712; reg647=reg64*reg647;
    reg519=reg64*reg449; reg665=reg64*reg665; reg521=reg64*reg456; reg856=reg64*reg856; reg617=reg64*reg617;
    reg298=reg64*reg298; reg894=reg64*reg894; reg579=reg64*reg579; reg699=reg64*reg699; reg516=reg64*reg516;
    reg819=reg64*reg819; reg697=reg64*reg697; reg383=reg64*reg383; reg289=reg64*reg289; reg589=reg64*reg589;
    reg522=ponderation*reg259; reg766=reg64*reg766; reg886=reg64*reg886; reg696=reg64*reg696; reg545=reg64*reg545;
    reg375=reg64*reg375; reg491=reg64*reg491; reg525=ponderation*reg373; reg763=reg64*reg763; reg543=reg64*reg543;
    reg813=reg64*reg813; reg693=reg64*reg693; reg417=reg64*reg417; reg528=ponderation*reg248; reg529=ponderation*reg271;
    reg334=reg64*reg334; reg351=reg64*reg351; reg671=reg64*reg671; reg854=reg64*reg854; reg253=reg64*reg253;
    reg270=reg64*reg270; reg702=reg64*reg702; reg560=reg64*reg560; reg574=reg64*reg574; reg316=reg64*reg316;
    reg774=reg64*reg774; reg345=reg64*reg345; reg701=reg64*reg701; reg821=reg64*reg821; reg315=reg64*reg315;
    reg391=reg64*reg391; reg531=ponderation*reg264; reg392=reg64*reg392; reg577=reg64*reg577; reg306=reg64*reg306;
    reg275=reg64*reg275; reg537=ponderation*reg300; reg580=reg64*reg580; reg765=reg64*reg765; reg201=reg64*reg201;
    reg539=ponderation*reg347; reg781=reg64*reg781; reg328=reg64*reg328; reg843=reg64*reg843; reg506=reg64*reg506;
    reg514=reg64*reg514; reg558=reg64*reg558; reg498=reg64*reg498; reg513=reg64*reg513; reg333=reg64*reg333;
    reg340=reg64*reg340; reg463=reg64*reg463; reg542=reg64*reg542; reg317=reg64*reg317; reg524=reg64*reg524;
    reg320=reg64*reg320; reg526=reg64*reg526; reg786=reg64*reg786; reg800=reg64*reg800; reg582=reg64*reg582;
    reg477=reg64*reg477; reg313=reg64*reg313; reg541=ponderation*reg399; reg782=reg64*reg782; reg489=reg64*reg489;
    reg497=reg64*reg497; reg806=reg64*reg806; reg414=reg64*reg414; reg482=reg64*reg482; reg544=ponderation*reg291;
    reg756=reg64*reg756; reg540=reg64*reg540; reg566=reg64*reg566; reg749=reg64*reg749; reg692=reg64*reg692;
    reg557=reg64*reg557; reg690=reg64*reg690; reg834=reg64*reg834; reg772=reg64*reg772; reg408=reg64*reg408;
    reg356=reg64*reg356; reg325=reg64*reg325; reg360=reg64*reg360; reg458=reg64*reg458; reg534=reg64*reg534;
    reg546=ponderation*reg269; reg570=reg64*reg570; reg777=reg64*reg777; reg496=reg64*reg496; reg847=reg64*reg847;
    reg400=reg64*reg400; reg551=reg64*reg551; reg484=reg64*reg484; reg857=reg64*reg857; reg425=reg64*reg425;
    reg687=reg64*reg687; reg841=reg64*reg841; reg512=reg64*reg512; reg547=ponderation*reg343; reg867=reg64*reg867;
    reg549=ponderation*reg324; reg548=reg64*reg548; reg571=reg64*reg571; reg527=reg64*reg527; reg249=reg64*reg249;
    reg377=reg64*reg377; reg866=reg64*reg866; reg872=reg64*reg872; reg550=ponderation*reg314; reg217=reg64*reg217;
    reg832=reg64*reg832; reg376=reg64*reg376; reg379=reg64*reg379; reg837=reg64*reg837; reg493=reg64*reg493;
    reg419=reg64*reg419; reg685=reg64*reg685; reg505=reg64*reg505; reg581=reg64*reg581; reg726=reg64*reg726;
    reg890=reg64*reg890; reg507=reg64*reg507; reg583=reg64*reg583; reg724=reg64*reg724; reg474=reg64*reg474;
    reg357=reg64*reg357; reg887=reg64*reg887; reg381=reg64*reg381; reg533=reg64*reg533; reg827=reg64*reg827;
    reg552=ponderation*reg176; reg352=reg64*reg352; reg286=reg64*reg286; reg523=reg64*reg523; reg364=reg64*reg364;
    reg556=reg64*reg556; reg559=ponderation*reg279; reg859=reg64*reg859; reg427=reg64*reg427; reg509=reg64*reg509;
    reg346=reg64*reg346; reg348=reg64*reg348; reg262=reg64*reg262; reg367=reg64*reg367; reg499=reg64*reg499;
    reg681=reg64*reg681; reg426=reg64*reg426; reg261=reg64*reg261; reg755=reg64*reg755; reg536=reg64*reg536;
    reg679=reg64*reg679; reg447=reg64*reg447; reg363=reg64*reg363; reg881=reg64*reg881; reg877=reg64*reg877;
    reg561=ponderation*reg483; reg677=reg64*reg677; reg530=reg64*reg530; reg674=reg64*reg674; reg428=reg64*reg428;
    reg563=reg64*reg563; reg362=reg64*reg362; reg803=reg64*reg803; reg520=reg64*reg520; reg569=ponderation*reg388;
    reg411=reg64*reg411; reg572=ponderation*reg301; reg234=reg64*reg234; reg353=reg64*reg353; reg452=reg64*reg452;
    reg554=reg64*reg554; reg410=reg64*reg410; reg869=reg64*reg869; reg683=reg64*reg683; reg487=reg64*reg487;
    reg267=reg64*reg267; reg795=reg64*reg795; reg372=reg64*reg372; T tmp_20_22=ponderation*reg647; T tmp_3_13=ponderation*reg832;
    T tmp_6_16=ponderation*reg491; T tmp_0_6=ponderation*reg421; T tmp_3_15=ponderation*reg270; T tmp_6_21=ponderation*reg267; T tmp_6_22=ponderation*reg505;
    T tmp_6_12=ponderation*reg477; T tmp_6_14=ponderation*reg534; T tmp_6_15=ponderation*reg482; T tmp_18_21=ponderation*reg670; T tmp_6_19=ponderation*reg520;
    T tmp_19_21=ponderation*reg642; T tmp_6_23=ponderation*reg571; T tmp_7_8=ponderation*reg583; T tmp_3_17=ponderation*reg749; T tmp_0_8=ponderation*reg836;
    T tmp_0_7=ponderation*reg761; T tmp_19_22=ponderation*reg455; T tmp_20_20=ponderation*reg420; T tmp_18_23=ponderation*reg638; T tmp_6_17=ponderation*reg516;
    T tmp_20_21=ponderation*reg645; T tmp_18_22=ponderation*reg636; T tmp_7_7=ponderation*reg262; T tmp_3_14=ponderation*reg755; T tmp_19_20=ponderation*reg640;
    T tmp_6_13=ponderation*reg498; T tmp_19_19=ponderation*reg450; T tmp_3_16=ponderation*reg765; T tmp_6_18=ponderation*reg275; T tmp_6_20=ponderation*reg530;
    T tmp_3_12=ponderation*reg827; T tmp_19_23=ponderation*reg644; T tmp_3_18=ponderation*reg484; T tmp_4_23=ponderation*reg493; reg262=ponderation*reg355;
    sollicitation[indices[4]+2]+=reg262; T tmp_5_5=ponderation*reg452; reg267=ponderation*reg341; sollicitation[indices[4]+1]+=reg267; T tmp_9_18=ponderation*reg869;
    reg270=ponderation*reg397; sollicitation[indices[4]+0]+=reg270; T tmp_5_6=ponderation*reg487; T tmp_1_4=ponderation*reg838; sollicitation[indices[3]+2]+=-reg595;
    T tmp_5_7=ponderation*reg499; T tmp_9_17=-reg525; T tmp_5_8=ponderation*reg447; reg275=ponderation*reg384; sollicitation[indices[3]+1]+=reg275;
    T tmp_4_9=ponderation*reg877; sollicitation[indices[3]+0]+=-reg603; T tmp_5_9=-reg561; T tmp_5_10=ponderation*reg563; T tmp_1_5=ponderation*reg874;
    reg420=ponderation*reg515; sollicitation[indices[2]+2]+=reg420; T tmp_4_8=ponderation*reg854; reg421=ponderation*reg511; sollicitation[indices[2]+1]+=reg421;
    T tmp_5_11=ponderation*reg253; reg253=ponderation*reg486; sollicitation[indices[7]+2]+=reg253; T tmp_4_14=ponderation*reg474; T tmp_4_15=ponderation*reg507;
    reg447=ponderation*reg473; sollicitation[indices[7]+1]+=reg447; T tmp_4_16=ponderation*reg533; reg450=ponderation*reg446; sollicitation[indices[7]+0]+=reg450;
    T tmp_4_13=ponderation*reg427; T tmp_4_17=ponderation*reg523; reg427=ponderation*reg403; sollicitation[indices[6]+2]+=reg427; T tmp_4_18=ponderation*reg509;
    reg452=ponderation*reg402; sollicitation[indices[6]+1]+=reg452; T tmp_4_12=ponderation*reg857; reg455=ponderation*reg501; sollicitation[indices[6]+0]+=reg455;
    T tmp_4_19=ponderation*reg425; T tmp_4_20=ponderation*reg512; reg425=ponderation*reg500; sollicitation[indices[5]+2]+=reg425; T tmp_4_11=ponderation*reg872;
    T tmp_4_21=ponderation*reg527; sollicitation[indices[5]+1]+=-reg629; T tmp_4_22=ponderation*reg217; T tmp_1_3=ponderation*reg852; reg217=ponderation*reg488;
    sollicitation[indices[5]+0]+=reg217; T tmp_4_10=ponderation*reg353; T tmp_23_23=ponderation*reg424; T tmp_5_21=ponderation*reg540; T tmp_3_23=ponderation*reg772;
    T tmp_5_22=ponderation*reg557; T tmp_22_23=ponderation*reg648; T tmp_5_23=ponderation*reg408; T tmp_0_1=ponderation*reg768; T tmp_22_22=ponderation*reg430;
    T tmp_3_22=ponderation*reg777; T tmp_6_6=ponderation*reg400; T tmp_21_23=ponderation*reg431; T tmp_6_7=ponderation*reg551; T tmp_3_21=ponderation*reg463;
    T tmp_0_2=ponderation*reg776; T tmp_21_22=ponderation*reg651; T tmp_6_8=ponderation*reg558; T tmp_6_9=ponderation*reg542; T tmp_21_21=ponderation*reg653;
    T tmp_3_20=ponderation*reg786; T tmp_6_10=ponderation*reg582; T tmp_6_11=-reg541; T tmp_20_23=ponderation*reg435; T tmp_3_19=ponderation*reg782;
    T tmp_5_12=ponderation*reg560; reg353=ponderation*reg510; sollicitation[indices[2]+0]+=reg353; T tmp_4_7=ponderation*reg391; T tmp_5_13=ponderation*reg345;
    reg345=ponderation*reg508; sollicitation[indices[1]+2]+=reg345; T tmp_5_14=ponderation*reg392; reg391=ponderation*reg521; sollicitation[indices[1]+1]+=reg391;
    T tmp_4_6=ponderation*reg894; T tmp_1_6=ponderation*reg856; reg392=ponderation*reg519; sollicitation[indices[1]+0]+=reg392; T tmp_5_15=ponderation*reg580;
    reg400=ponderation*reg518; sollicitation[indices[0]+2]+=reg400; T tmp_5_16=ponderation*reg579; T tmp_5_17=ponderation*reg383; T tmp_4_5=ponderation*reg886;
    T tmp_1_7=ponderation*reg382; reg382=ponderation*reg517; sollicitation[indices[0]+1]+=reg382; T tmp_5_18=ponderation*reg545; reg383=ponderation*reg504;
    sollicitation[indices[0]+0]+=reg383; T tmp_5_19=ponderation*reg543; T tmp_4_4=ponderation*reg414; T tmp_5_20=ponderation*reg417; T tmp_0_0=ponderation*reg448;
    T tmp_14_14=ponderation*reg274; T tmp_10_12=ponderation*reg701; T tmp_2_7=ponderation*reg774; T tmp_13_23=ponderation*reg719; T tmp_10_13=ponderation*reg316;
    T tmp_10_14=ponderation*reg702; T tmp_1_14=ponderation*reg778; T tmp_13_22=ponderation*reg272; T tmp_2_6=ponderation*reg411; T tmp_10_15=ponderation*reg671;
    T tmp_13_21=ponderation*reg718; T tmp_10_16=-reg569; T tmp_2_5=ponderation*reg362; T tmp_1_15=ponderation*reg769; T tmp_13_20=ponderation*reg715;
    T tmp_10_17=ponderation*reg674; T tmp_10_18=ponderation*reg677; T tmp_13_19=ponderation*reg265; T tmp_2_4=ponderation*reg881; T tmp_1_16=ponderation*reg268;
    T tmp_10_19=ponderation*reg363; T tmp_13_18=ponderation*reg714; T tmp_10_20=ponderation*reg679; T tmp_2_3=ponderation*reg426; T tmp_14_22=ponderation*reg623;
    T tmp_0_23=ponderation*reg759; T tmp_2_12=ponderation*reg766; T tmp_14_21=ponderation*reg621; T tmp_9_19=ponderation*reg696; T tmp_14_20=ponderation*reg280;
    T tmp_9_20=-reg522; T tmp_2_11=ponderation*reg289; T tmp_1_1=ponderation*reg480; T tmp_14_19=ponderation*reg619; T tmp_9_21=ponderation*reg697;
    T tmp_2_10=ponderation*reg781; T tmp_14_18=ponderation*reg617; T tmp_9_22=ponderation*reg699; T tmp_1_2=ponderation*reg743; T tmp_14_17=ponderation*reg276;
    T tmp_9_23=-reg539; T tmp_2_9=-reg537; T tmp_1_12=ponderation*reg887; T tmp_14_16=ponderation*reg721; T tmp_10_10=ponderation*reg306;
    T tmp_2_8=ponderation*reg315; T tmp_14_15=ponderation*reg720; T tmp_10_11=-reg531; T tmp_1_13=ponderation*reg277; T tmp_1_20=ponderation*reg867;
    T tmp_11_17=-reg549; T tmp_11_18=-reg547; T tmp_12_20=ponderation*reg736; T tmp_1_19=ponderation*reg348; T tmp_11_19=ponderation*reg687;
    T tmp_0_10=ponderation*reg868; T tmp_12_19=ponderation*reg734; T tmp_11_20=ponderation*reg346; T tmp_1_8=ponderation*reg879; T tmp_12_18=ponderation*reg257;
    T tmp_1_18=ponderation*reg859; T tmp_11_21=-reg559; T tmp_11_22=ponderation*reg364; T tmp_12_17=ponderation*reg731; T tmp_1_17=ponderation*reg381;
    T tmp_1_9=ponderation*reg349; T tmp_11_23=ponderation*reg352; T tmp_12_16=ponderation*reg728; T tmp_12_12=ponderation*reg357; T tmp_12_15=ponderation*reg250;
    T tmp_12_13=ponderation*reg724; T tmp_1_10=ponderation*reg260; T tmp_1_11=ponderation*reg890; T tmp_12_14=ponderation*reg726; T tmp_0_13=ponderation*reg815;
    T tmp_13_17=ponderation*reg712; T tmp_10_21=ponderation*reg681; T tmp_13_16=ponderation*reg298; T tmp_10_22=ponderation*reg367; T tmp_2_2=ponderation*reg372;
    T tmp_0_14=ponderation*reg479; T tmp_13_15=ponderation*reg710; T tmp_10_23=ponderation*reg683; T tmp_11_11=ponderation*reg410; T tmp_1_23=ponderation*reg234;
    T tmp_13_14=ponderation*reg707; T tmp_11_12=-reg572; T tmp_13_13=ponderation*reg292; T tmp_11_13=ponderation*reg685; T tmp_1_22=ponderation*reg379;
    T tmp_0_15=ponderation*reg295; T tmp_12_23=ponderation*reg705; T tmp_11_14=ponderation*reg376; T tmp_1_21=ponderation*reg866; T tmp_11_15=-reg550;
    T tmp_12_22=ponderation*reg703; T tmp_11_16=ponderation*reg377; T tmp_0_9=ponderation*reg871; T tmp_12_21=ponderation*reg285; T tmp_7_18=ponderation*reg428;
    T tmp_17_19=ponderation*reg661; T tmp_7_19=ponderation*reg351; T tmp_3_6=ponderation*reg821; T tmp_17_18=ponderation*reg659; T tmp_7_20=ponderation*reg574;
    T tmp_0_11=-reg393; T tmp_17_17=ponderation*reg387; T tmp_7_21=ponderation*reg577; T tmp_3_5=ponderation*reg819; T tmp_7_22=ponderation*reg201;
    T tmp_16_23=ponderation*reg656; T tmp_7_23=ponderation*reg589; T tmp_0_12=ponderation*reg389; T tmp_16_22=ponderation*reg385; T tmp_8_8=ponderation*reg375;
    T tmp_3_4=ponderation*reg813; T tmp_16_21=ponderation*reg610; T tmp_8_9=-reg528; T tmp_8_10=ponderation*reg566; T tmp_3_3=ponderation*reg360;
    T tmp_16_20=ponderation*reg607; T tmp_8_11=ponderation*reg458; T tmp_16_19=ponderation*reg413; T tmp_0_19=ponderation*reg794; T tmp_7_9=ponderation*reg581;
    T tmp_3_11=-reg552; T tmp_18_20=ponderation*reg669; T tmp_7_10=ponderation*reg286; T tmp_18_19=ponderation*reg667; T tmp_7_11=ponderation*reg556;
    T tmp_3_10=ponderation*reg841; T tmp_0_20=ponderation*reg824; T tmp_18_18=ponderation*reg437; T tmp_7_12=ponderation*reg548; T tmp_7_13=ponderation*reg249;
    T tmp_3_9=ponderation*reg837; T tmp_0_21=ponderation*reg439; T tmp_17_23=ponderation*reg481; T tmp_7_14=ponderation*reg419; T tmp_17_22=ponderation*reg665;
    T tmp_7_15=ponderation*reg554; T tmp_3_8=ponderation*reg795; T tmp_17_21=ponderation*reg663; T tmp_7_16=ponderation*reg261; T tmp_0_22=ponderation*reg398;
    T tmp_7_17=ponderation*reg536; T tmp_17_20=ponderation*reg464; T tmp_3_7=ponderation*reg803; T tmp_15_20=ponderation*reg594; T tmp_0_3=ponderation*reg297;
    T tmp_8_23=ponderation*reg328; T tmp_15_19=ponderation*reg591; T tmp_2_18=ponderation*reg847; T tmp_9_9=ponderation*reg496; T tmp_15_18=ponderation*reg466;
    T tmp_9_10=-reg546; T tmp_2_17=ponderation*reg325; T tmp_0_4=ponderation*reg846; T tmp_15_17=ponderation*reg634; T tmp_9_11=ponderation*reg356;
    T tmp_2_16=ponderation*reg834; T tmp_15_16=-reg503; T tmp_9_12=ponderation*reg690; T tmp_0_5=ponderation*reg826; T tmp_9_13=ponderation*reg692;
    T tmp_2_15=ponderation*reg756; T tmp_15_15=ponderation*reg627; T tmp_9_14=-reg544; T tmp_2_14=ponderation*reg334; T tmp_14_23=ponderation*reg205;
    T tmp_9_15=-reg529; T tmp_9_16=ponderation*reg693; T tmp_2_13=ponderation*reg763; T tmp_8_12=ponderation*reg570; T tmp_2_23=ponderation*reg317;
    T tmp_0_16=ponderation*reg875; T tmp_8_13=ponderation*reg506; T tmp_16_18=ponderation*reg606; T tmp_8_14=ponderation*reg320; T tmp_2_22=ponderation*reg806;
    T tmp_16_17=-reg495; T tmp_8_15=ponderation*reg497; T tmp_0_17=ponderation*reg820; T tmp_16_16=ponderation*reg601; T tmp_8_16=ponderation*reg489;
    T tmp_8_17=ponderation*reg313; T tmp_2_21=ponderation*reg800; T tmp_15_23=ponderation*reg600; T tmp_8_18=ponderation*reg526; T tmp_0_18=ponderation*reg405;
    T tmp_15_22=ponderation*reg598; T tmp_8_19=ponderation*reg524; T tmp_2_20=ponderation*reg333; T tmp_8_20=ponderation*reg340; T tmp_15_21=ponderation*reg222;
    T tmp_8_21=ponderation*reg513; T tmp_2_19=ponderation*reg843; T tmp_8_22=ponderation*reg514;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[0]+2) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_3;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_4;
    matrix(indices[0]+0,indices[1]+2) += tmp_0_5;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_7;
    matrix(indices[0]+0,indices[2]+2) += tmp_0_8;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_9;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_10;
    matrix(indices[0]+0,indices[3]+2) += tmp_0_11;
    matrix(indices[0]+0,indices[4]+0) += tmp_0_12;
    matrix(indices[0]+0,indices[4]+1) += tmp_0_13;
    matrix(indices[0]+0,indices[4]+2) += tmp_0_14;
    matrix(indices[0]+0,indices[5]+0) += tmp_0_15;
    matrix(indices[0]+0,indices[5]+1) += tmp_0_16;
    matrix(indices[0]+0,indices[5]+2) += tmp_0_17;
    matrix(indices[0]+0,indices[6]+0) += tmp_0_18;
    matrix(indices[0]+0,indices[6]+1) += tmp_0_19;
    matrix(indices[0]+0,indices[6]+2) += tmp_0_20;
    matrix(indices[0]+0,indices[7]+0) += tmp_0_21;
    matrix(indices[0]+0,indices[7]+1) += tmp_0_22;
    matrix(indices[0]+0,indices[7]+2) += tmp_0_23;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[0]+2) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_3;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_4;
    matrix(indices[0]+1,indices[1]+2) += tmp_1_5;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_7;
    matrix(indices[0]+1,indices[2]+2) += tmp_1_8;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_9;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_10;
    matrix(indices[0]+1,indices[3]+2) += tmp_1_11;
    matrix(indices[0]+1,indices[4]+0) += tmp_1_12;
    matrix(indices[0]+1,indices[4]+1) += tmp_1_13;
    matrix(indices[0]+1,indices[4]+2) += tmp_1_14;
    matrix(indices[0]+1,indices[5]+0) += tmp_1_15;
    matrix(indices[0]+1,indices[5]+1) += tmp_1_16;
    matrix(indices[0]+1,indices[5]+2) += tmp_1_17;
    matrix(indices[0]+1,indices[6]+0) += tmp_1_18;
    matrix(indices[0]+1,indices[6]+1) += tmp_1_19;
    matrix(indices[0]+1,indices[6]+2) += tmp_1_20;
    matrix(indices[0]+1,indices[7]+0) += tmp_1_21;
    matrix(indices[0]+1,indices[7]+1) += tmp_1_22;
    matrix(indices[0]+1,indices[7]+2) += tmp_1_23;
    matrix(indices[0]+2,indices[0]+2) += tmp_2_2;
    matrix(indices[0]+2,indices[1]+0) += tmp_2_3;
    matrix(indices[0]+2,indices[1]+1) += tmp_2_4;
    matrix(indices[0]+2,indices[1]+2) += tmp_2_5;
    matrix(indices[0]+2,indices[2]+0) += tmp_2_6;
    matrix(indices[0]+2,indices[2]+1) += tmp_2_7;
    matrix(indices[0]+2,indices[2]+2) += tmp_2_8;
    matrix(indices[0]+2,indices[3]+0) += tmp_2_9;
    matrix(indices[0]+2,indices[3]+1) += tmp_2_10;
    matrix(indices[0]+2,indices[3]+2) += tmp_2_11;
    matrix(indices[0]+2,indices[4]+0) += tmp_2_12;
    matrix(indices[0]+2,indices[4]+1) += tmp_2_13;
    matrix(indices[0]+2,indices[4]+2) += tmp_2_14;
    matrix(indices[0]+2,indices[5]+0) += tmp_2_15;
    matrix(indices[0]+2,indices[5]+1) += tmp_2_16;
    matrix(indices[0]+2,indices[5]+2) += tmp_2_17;
    matrix(indices[0]+2,indices[6]+0) += tmp_2_18;
    matrix(indices[0]+2,indices[6]+1) += tmp_2_19;
    matrix(indices[0]+2,indices[6]+2) += tmp_2_20;
    matrix(indices[0]+2,indices[7]+0) += tmp_2_21;
    matrix(indices[0]+2,indices[7]+1) += tmp_2_22;
    matrix(indices[0]+2,indices[7]+2) += tmp_2_23;
    matrix(indices[1]+0,indices[1]+0) += tmp_3_3;
    matrix(indices[1]+0,indices[1]+1) += tmp_3_4;
    matrix(indices[1]+0,indices[1]+2) += tmp_3_5;
    matrix(indices[1]+0,indices[2]+0) += tmp_3_6;
    matrix(indices[1]+0,indices[2]+1) += tmp_3_7;
    matrix(indices[1]+0,indices[2]+2) += tmp_3_8;
    matrix(indices[1]+0,indices[3]+0) += tmp_3_9;
    matrix(indices[1]+0,indices[3]+1) += tmp_3_10;
    matrix(indices[1]+0,indices[3]+2) += tmp_3_11;
    matrix(indices[1]+0,indices[4]+0) += tmp_3_12;
    matrix(indices[1]+0,indices[4]+1) += tmp_3_13;
    matrix(indices[1]+0,indices[4]+2) += tmp_3_14;
    matrix(indices[1]+0,indices[5]+0) += tmp_3_15;
    matrix(indices[1]+0,indices[5]+1) += tmp_3_16;
    matrix(indices[1]+0,indices[5]+2) += tmp_3_17;
    matrix(indices[1]+0,indices[6]+0) += tmp_3_18;
    matrix(indices[1]+0,indices[6]+1) += tmp_3_19;
    matrix(indices[1]+0,indices[6]+2) += tmp_3_20;
    matrix(indices[1]+0,indices[7]+0) += tmp_3_21;
    matrix(indices[1]+0,indices[7]+1) += tmp_3_22;
    matrix(indices[1]+0,indices[7]+2) += tmp_3_23;
    matrix(indices[1]+1,indices[1]+1) += tmp_4_4;
    matrix(indices[1]+1,indices[1]+2) += tmp_4_5;
    matrix(indices[1]+1,indices[2]+0) += tmp_4_6;
    matrix(indices[1]+1,indices[2]+1) += tmp_4_7;
    matrix(indices[1]+1,indices[2]+2) += tmp_4_8;
    matrix(indices[1]+1,indices[3]+0) += tmp_4_9;
    matrix(indices[1]+1,indices[3]+1) += tmp_4_10;
    matrix(indices[1]+1,indices[3]+2) += tmp_4_11;
    matrix(indices[1]+1,indices[4]+0) += tmp_4_12;
    matrix(indices[1]+1,indices[4]+1) += tmp_4_13;
    matrix(indices[1]+1,indices[4]+2) += tmp_4_14;
    matrix(indices[1]+1,indices[5]+0) += tmp_4_15;
    matrix(indices[1]+1,indices[5]+1) += tmp_4_16;
    matrix(indices[1]+1,indices[5]+2) += tmp_4_17;
    matrix(indices[1]+1,indices[6]+0) += tmp_4_18;
    matrix(indices[1]+1,indices[6]+1) += tmp_4_19;
    matrix(indices[1]+1,indices[6]+2) += tmp_4_20;
    matrix(indices[1]+1,indices[7]+0) += tmp_4_21;
    matrix(indices[1]+1,indices[7]+1) += tmp_4_22;
    matrix(indices[1]+1,indices[7]+2) += tmp_4_23;
    matrix(indices[1]+2,indices[1]+2) += tmp_5_5;
    matrix(indices[1]+2,indices[2]+0) += tmp_5_6;
    matrix(indices[1]+2,indices[2]+1) += tmp_5_7;
    matrix(indices[1]+2,indices[2]+2) += tmp_5_8;
    matrix(indices[1]+2,indices[3]+0) += tmp_5_9;
    matrix(indices[1]+2,indices[3]+1) += tmp_5_10;
    matrix(indices[1]+2,indices[3]+2) += tmp_5_11;
    matrix(indices[1]+2,indices[4]+0) += tmp_5_12;
    matrix(indices[1]+2,indices[4]+1) += tmp_5_13;
    matrix(indices[1]+2,indices[4]+2) += tmp_5_14;
    matrix(indices[1]+2,indices[5]+0) += tmp_5_15;
    matrix(indices[1]+2,indices[5]+1) += tmp_5_16;
    matrix(indices[1]+2,indices[5]+2) += tmp_5_17;
    matrix(indices[1]+2,indices[6]+0) += tmp_5_18;
    matrix(indices[1]+2,indices[6]+1) += tmp_5_19;
    matrix(indices[1]+2,indices[6]+2) += tmp_5_20;
    matrix(indices[1]+2,indices[7]+0) += tmp_5_21;
    matrix(indices[1]+2,indices[7]+1) += tmp_5_22;
    matrix(indices[1]+2,indices[7]+2) += tmp_5_23;
    matrix(indices[2]+0,indices[2]+0) += tmp_6_6;
    matrix(indices[2]+0,indices[2]+1) += tmp_6_7;
    matrix(indices[2]+0,indices[2]+2) += tmp_6_8;
    matrix(indices[2]+0,indices[3]+0) += tmp_6_9;
    matrix(indices[2]+0,indices[3]+1) += tmp_6_10;
    matrix(indices[2]+0,indices[3]+2) += tmp_6_11;
    matrix(indices[2]+0,indices[4]+0) += tmp_6_12;
    matrix(indices[2]+0,indices[4]+1) += tmp_6_13;
    matrix(indices[2]+0,indices[4]+2) += tmp_6_14;
    matrix(indices[2]+0,indices[5]+0) += tmp_6_15;
    matrix(indices[2]+0,indices[5]+1) += tmp_6_16;
    matrix(indices[2]+0,indices[5]+2) += tmp_6_17;
    matrix(indices[2]+0,indices[6]+0) += tmp_6_18;
    matrix(indices[2]+0,indices[6]+1) += tmp_6_19;
    matrix(indices[2]+0,indices[6]+2) += tmp_6_20;
    matrix(indices[2]+0,indices[7]+0) += tmp_6_21;
    matrix(indices[2]+0,indices[7]+1) += tmp_6_22;
    matrix(indices[2]+0,indices[7]+2) += tmp_6_23;
    matrix(indices[2]+1,indices[2]+1) += tmp_7_7;
    matrix(indices[2]+1,indices[2]+2) += tmp_7_8;
    matrix(indices[2]+1,indices[3]+0) += tmp_7_9;
    matrix(indices[2]+1,indices[3]+1) += tmp_7_10;
    matrix(indices[2]+1,indices[3]+2) += tmp_7_11;
    matrix(indices[2]+1,indices[4]+0) += tmp_7_12;
    matrix(indices[2]+1,indices[4]+1) += tmp_7_13;
    matrix(indices[2]+1,indices[4]+2) += tmp_7_14;
    matrix(indices[2]+1,indices[5]+0) += tmp_7_15;
    matrix(indices[2]+1,indices[5]+1) += tmp_7_16;
    matrix(indices[2]+1,indices[5]+2) += tmp_7_17;
    matrix(indices[2]+1,indices[6]+0) += tmp_7_18;
    matrix(indices[2]+1,indices[6]+1) += tmp_7_19;
    matrix(indices[2]+1,indices[6]+2) += tmp_7_20;
    matrix(indices[2]+1,indices[7]+0) += tmp_7_21;
    matrix(indices[2]+1,indices[7]+1) += tmp_7_22;
    matrix(indices[2]+1,indices[7]+2) += tmp_7_23;
    matrix(indices[2]+2,indices[2]+2) += tmp_8_8;
    matrix(indices[2]+2,indices[3]+0) += tmp_8_9;
    matrix(indices[2]+2,indices[3]+1) += tmp_8_10;
    matrix(indices[2]+2,indices[3]+2) += tmp_8_11;
    matrix(indices[2]+2,indices[4]+0) += tmp_8_12;
    matrix(indices[2]+2,indices[4]+1) += tmp_8_13;
    matrix(indices[2]+2,indices[4]+2) += tmp_8_14;
    matrix(indices[2]+2,indices[5]+0) += tmp_8_15;
    matrix(indices[2]+2,indices[5]+1) += tmp_8_16;
    matrix(indices[2]+2,indices[5]+2) += tmp_8_17;
    matrix(indices[2]+2,indices[6]+0) += tmp_8_18;
    matrix(indices[2]+2,indices[6]+1) += tmp_8_19;
    matrix(indices[2]+2,indices[6]+2) += tmp_8_20;
    matrix(indices[2]+2,indices[7]+0) += tmp_8_21;
    matrix(indices[2]+2,indices[7]+1) += tmp_8_22;
    matrix(indices[2]+2,indices[7]+2) += tmp_8_23;
    matrix(indices[3]+0,indices[3]+0) += tmp_9_9;
    matrix(indices[3]+0,indices[3]+1) += tmp_9_10;
    matrix(indices[3]+0,indices[3]+2) += tmp_9_11;
    matrix(indices[3]+0,indices[4]+0) += tmp_9_12;
    matrix(indices[3]+0,indices[4]+1) += tmp_9_13;
    matrix(indices[3]+0,indices[4]+2) += tmp_9_14;
    matrix(indices[3]+0,indices[5]+0) += tmp_9_15;
    matrix(indices[3]+0,indices[5]+1) += tmp_9_16;
    matrix(indices[3]+0,indices[5]+2) += tmp_9_17;
    matrix(indices[3]+0,indices[6]+0) += tmp_9_18;
    matrix(indices[3]+0,indices[6]+1) += tmp_9_19;
    matrix(indices[3]+0,indices[6]+2) += tmp_9_20;
    matrix(indices[3]+0,indices[7]+0) += tmp_9_21;
    matrix(indices[3]+0,indices[7]+1) += tmp_9_22;
    matrix(indices[3]+0,indices[7]+2) += tmp_9_23;
    matrix(indices[3]+1,indices[3]+1) += tmp_10_10;
    matrix(indices[3]+1,indices[3]+2) += tmp_10_11;
    matrix(indices[3]+1,indices[4]+0) += tmp_10_12;
    matrix(indices[3]+1,indices[4]+1) += tmp_10_13;
    matrix(indices[3]+1,indices[4]+2) += tmp_10_14;
    matrix(indices[3]+1,indices[5]+0) += tmp_10_15;
    matrix(indices[3]+1,indices[5]+1) += tmp_10_16;
    matrix(indices[3]+1,indices[5]+2) += tmp_10_17;
    matrix(indices[3]+1,indices[6]+0) += tmp_10_18;
    matrix(indices[3]+1,indices[6]+1) += tmp_10_19;
    matrix(indices[3]+1,indices[6]+2) += tmp_10_20;
    matrix(indices[3]+1,indices[7]+0) += tmp_10_21;
    matrix(indices[3]+1,indices[7]+1) += tmp_10_22;
    matrix(indices[3]+1,indices[7]+2) += tmp_10_23;
    matrix(indices[3]+2,indices[3]+2) += tmp_11_11;
    matrix(indices[3]+2,indices[4]+0) += tmp_11_12;
    matrix(indices[3]+2,indices[4]+1) += tmp_11_13;
    matrix(indices[3]+2,indices[4]+2) += tmp_11_14;
    matrix(indices[3]+2,indices[5]+0) += tmp_11_15;
    matrix(indices[3]+2,indices[5]+1) += tmp_11_16;
    matrix(indices[3]+2,indices[5]+2) += tmp_11_17;
    matrix(indices[3]+2,indices[6]+0) += tmp_11_18;
    matrix(indices[3]+2,indices[6]+1) += tmp_11_19;
    matrix(indices[3]+2,indices[6]+2) += tmp_11_20;
    matrix(indices[3]+2,indices[7]+0) += tmp_11_21;
    matrix(indices[3]+2,indices[7]+1) += tmp_11_22;
    matrix(indices[3]+2,indices[7]+2) += tmp_11_23;
    matrix(indices[4]+0,indices[4]+0) += tmp_12_12;
    matrix(indices[4]+0,indices[4]+1) += tmp_12_13;
    matrix(indices[4]+0,indices[4]+2) += tmp_12_14;
    matrix(indices[4]+0,indices[5]+0) += tmp_12_15;
    matrix(indices[4]+0,indices[5]+1) += tmp_12_16;
    matrix(indices[4]+0,indices[5]+2) += tmp_12_17;
    matrix(indices[4]+0,indices[6]+0) += tmp_12_18;
    matrix(indices[4]+0,indices[6]+1) += tmp_12_19;
    matrix(indices[4]+0,indices[6]+2) += tmp_12_20;
    matrix(indices[4]+0,indices[7]+0) += tmp_12_21;
    matrix(indices[4]+0,indices[7]+1) += tmp_12_22;
    matrix(indices[4]+0,indices[7]+2) += tmp_12_23;
    matrix(indices[4]+1,indices[4]+1) += tmp_13_13;
    matrix(indices[4]+1,indices[4]+2) += tmp_13_14;
    matrix(indices[4]+1,indices[5]+0) += tmp_13_15;
    matrix(indices[4]+1,indices[5]+1) += tmp_13_16;
    matrix(indices[4]+1,indices[5]+2) += tmp_13_17;
    matrix(indices[4]+1,indices[6]+0) += tmp_13_18;
    matrix(indices[4]+1,indices[6]+1) += tmp_13_19;
    matrix(indices[4]+1,indices[6]+2) += tmp_13_20;
    matrix(indices[4]+1,indices[7]+0) += tmp_13_21;
    matrix(indices[4]+1,indices[7]+1) += tmp_13_22;
    matrix(indices[4]+1,indices[7]+2) += tmp_13_23;
    matrix(indices[4]+2,indices[4]+2) += tmp_14_14;
    matrix(indices[4]+2,indices[5]+0) += tmp_14_15;
    matrix(indices[4]+2,indices[5]+1) += tmp_14_16;
    matrix(indices[4]+2,indices[5]+2) += tmp_14_17;
    matrix(indices[4]+2,indices[6]+0) += tmp_14_18;
    matrix(indices[4]+2,indices[6]+1) += tmp_14_19;
    matrix(indices[4]+2,indices[6]+2) += tmp_14_20;
    matrix(indices[4]+2,indices[7]+0) += tmp_14_21;
    matrix(indices[4]+2,indices[7]+1) += tmp_14_22;
    matrix(indices[4]+2,indices[7]+2) += tmp_14_23;
    matrix(indices[5]+0,indices[5]+0) += tmp_15_15;
    matrix(indices[5]+0,indices[5]+1) += tmp_15_16;
    matrix(indices[5]+0,indices[5]+2) += tmp_15_17;
    matrix(indices[5]+0,indices[6]+0) += tmp_15_18;
    matrix(indices[5]+0,indices[6]+1) += tmp_15_19;
    matrix(indices[5]+0,indices[6]+2) += tmp_15_20;
    matrix(indices[5]+0,indices[7]+0) += tmp_15_21;
    matrix(indices[5]+0,indices[7]+1) += tmp_15_22;
    matrix(indices[5]+0,indices[7]+2) += tmp_15_23;
    matrix(indices[5]+1,indices[5]+1) += tmp_16_16;
    matrix(indices[5]+1,indices[5]+2) += tmp_16_17;
    matrix(indices[5]+1,indices[6]+0) += tmp_16_18;
    matrix(indices[5]+1,indices[6]+1) += tmp_16_19;
    matrix(indices[5]+1,indices[6]+2) += tmp_16_20;
    matrix(indices[5]+1,indices[7]+0) += tmp_16_21;
    matrix(indices[5]+1,indices[7]+1) += tmp_16_22;
    matrix(indices[5]+1,indices[7]+2) += tmp_16_23;
    matrix(indices[5]+2,indices[5]+2) += tmp_17_17;
    matrix(indices[5]+2,indices[6]+0) += tmp_17_18;
    matrix(indices[5]+2,indices[6]+1) += tmp_17_19;
    matrix(indices[5]+2,indices[6]+2) += tmp_17_20;
    matrix(indices[5]+2,indices[7]+0) += tmp_17_21;
    matrix(indices[5]+2,indices[7]+1) += tmp_17_22;
    matrix(indices[5]+2,indices[7]+2) += tmp_17_23;
    matrix(indices[6]+0,indices[6]+0) += tmp_18_18;
    matrix(indices[6]+0,indices[6]+1) += tmp_18_19;
    matrix(indices[6]+0,indices[6]+2) += tmp_18_20;
    matrix(indices[6]+0,indices[7]+0) += tmp_18_21;
    matrix(indices[6]+0,indices[7]+1) += tmp_18_22;
    matrix(indices[6]+0,indices[7]+2) += tmp_18_23;
    matrix(indices[6]+1,indices[6]+1) += tmp_19_19;
    matrix(indices[6]+1,indices[6]+2) += tmp_19_20;
    matrix(indices[6]+1,indices[7]+0) += tmp_19_21;
    matrix(indices[6]+1,indices[7]+1) += tmp_19_22;
    matrix(indices[6]+1,indices[7]+2) += tmp_19_23;
    matrix(indices[6]+2,indices[6]+2) += tmp_20_20;
    matrix(indices[6]+2,indices[7]+0) += tmp_20_21;
    matrix(indices[6]+2,indices[7]+1) += tmp_20_22;
    matrix(indices[6]+2,indices[7]+2) += tmp_20_23;
    matrix(indices[7]+0,indices[7]+0) += tmp_21_21;
    matrix(indices[7]+0,indices[7]+1) += tmp_21_22;
    matrix(indices[7]+0,indices[7]+2) += tmp_21_23;
    matrix(indices[7]+1,indices[7]+1) += tmp_22_22;
    matrix(indices[7]+1,indices[7]+2) += tmp_22_23;
    matrix(indices[7]+2,indices[7]+2) += tmp_23_23;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1-var_inter[2]; T reg1=1-var_inter[0]; T reg2=1-var_inter[1]; T reg3=reg2*reg1; T reg4=reg2*reg0;
    T reg5=reg2*var_inter[0]; T reg6=reg0*var_inter[0]; T reg7=reg0*reg1; T reg8=elem.pos(1)[1]*reg5; T reg9=elem.pos(0)[1]*reg3;
    T reg10=elem.pos(0)[2]*reg7; T reg11=reg0*var_inter[1]; T reg12=elem.pos(1)[2]*reg6; T reg13=elem.pos(0)[1]*reg7; T reg14=elem.pos(1)[1]*reg6;
    T reg15=elem.pos(1)[1]*reg4; T reg16=elem.pos(0)[1]*reg4; T reg17=elem.pos(1)[2]*reg4; T reg18=elem.pos(0)[2]*reg4; T reg19=var_inter[0]*var_inter[1];
    T reg20=elem.pos(0)[2]*reg3; T reg21=reg5*elem.pos(1)[2]; T reg22=reg21+reg20; reg15=reg15-reg16; T reg23=reg12+reg10;
    T reg24=elem.pos(2)[2]*reg6; T reg25=elem.pos(2)[1]*reg11; T reg26=reg19*elem.pos(2)[2]; T reg27=reg14+reg13; reg17=reg17-reg18;
    T reg28=reg8+reg9; T reg29=reg19*elem.pos(2)[1]; T reg30=var_inter[1]*reg1; T reg31=elem.pos(2)[1]*reg6; T reg32=elem.pos(2)[2]*reg11;
    reg31=reg31-reg27; T reg33=var_inter[2]*reg1; T reg34=reg22+reg26; reg25=reg15+reg25; reg15=elem.pos(3)[1]*reg11;
    T reg35=elem.pos(0)[0]*reg7; T reg36=elem.pos(1)[0]*reg6; reg32=reg17+reg32; reg17=elem.pos(3)[2]*reg11; T reg37=reg28+reg29;
    T reg38=reg30*elem.pos(3)[1]; reg24=reg24-reg23; T reg39=elem.pos(1)[0]*reg4; T reg40=reg30*elem.pos(3)[2]; T reg41=elem.pos(0)[0]*reg4;
    T reg42=elem.pos(3)[2]*reg7; T reg43=elem.pos(3)[1]*reg7; T reg44=reg2*var_inter[2]; T reg45=elem.pos(4)[1]*reg3; T reg46=reg40+reg34;
    T reg47=elem.pos(4)[2]*reg3; T reg48=reg37+reg38; reg32=reg32-reg17; T reg49=elem.pos(4)[2]*reg44; T reg50=elem.pos(0)[0]*reg3;
    T reg51=elem.pos(1)[0]*reg5; T reg52=elem.pos(4)[2]*reg33; reg24=reg42+reg24; reg42=reg36+reg35; T reg53=elem.pos(2)[0]*reg6;
    T reg54=var_inter[2]*var_inter[0]; T reg55=elem.pos(4)[1]*reg33; reg43=reg31+reg43; reg25=reg25-reg15; reg31=elem.pos(4)[1]*reg44;
    reg39=reg39-reg41; T reg56=elem.pos(2)[0]*reg11; T reg57=var_inter[2]*var_inter[1]; T reg58=elem.pos(5)[1]*reg54; reg43=reg43-reg55;
    T reg59=elem.pos(3)[0]*reg11; reg47=reg47-reg46; T reg60=reg5*elem.pos(5)[2]; T reg61=reg5*elem.pos(5)[1]; reg45=reg45-reg48;
    reg56=reg39+reg56; reg24=reg24-reg52; reg39=elem.pos(5)[2]*reg54; T reg62=reg19*elem.pos(2)[0]; reg32=reg32-reg49;
    T reg63=elem.pos(5)[2]*reg44; T reg64=reg51+reg50; T reg65=elem.pos(5)[1]*reg44; T reg66=elem.pos(3)[0]*reg7; reg25=reg25-reg31;
    reg53=reg53-reg42; T reg67=reg19*elem.pos(6)[1]; T reg68=reg30*elem.pos(3)[0]; reg61=reg45+reg61; reg45=elem.pos(6)[2]*reg54;
    T reg69=elem.pos(6)[1]*reg54; reg43=reg43-reg58; T reg70=reg64+reg62; reg24=reg24-reg39; reg53=reg66+reg53;
    reg66=reg19*elem.pos(6)[2]; T reg71=elem.pos(4)[0]*reg33; reg60=reg47+reg60; reg47=elem.pos(6)[2]*reg57; reg63=reg32+reg63;
    reg65=reg25+reg65; reg25=elem.pos(6)[1]*reg57; reg32=elem.pos(4)[0]*reg44; reg56=reg56-reg59; T reg72=reg68+reg70;
    T reg73=elem.pos(4)[0]*reg3; reg25=reg65+reg25; reg47=reg63+reg47; reg63=elem.pos(7)[1]*reg57; reg65=elem.pos(7)[2]*reg57;
    T reg74=elem.pos(7)[2]*reg33; reg45=reg24+reg45; reg24=reg30*elem.pos(7)[2]; reg66=reg60+reg66; reg53=reg53-reg71;
    reg60=elem.pos(5)[0]*reg54; reg67=reg61+reg67; reg61=reg30*elem.pos(7)[1]; T reg75=elem.pos(7)[1]*reg33; reg69=reg43+reg69;
    reg43=elem.pos(5)[0]*reg44; reg56=reg56-reg32; reg24=reg66+reg24; reg61=reg67+reg61; reg66=1+(*f.m).poisson_ratio;
    reg75=reg69+reg75; reg43=reg56+reg43; reg47=reg47-reg65; reg74=reg45+reg74; reg73=reg73-reg72;
    reg45=reg5*elem.pos(5)[0]; reg56=elem.pos(6)[0]*reg54; reg53=reg53-reg60; reg67=elem.pos(6)[0]*reg57; reg25=reg25-reg63;
    reg69=elem.pos(7)[0]*reg57; reg56=reg53+reg56; reg53=elem.pos(7)[0]*reg33; T reg76=reg75*reg24; reg67=reg43+reg67;
    reg43=reg25*reg24; T reg77=reg74*reg61; T reg78=reg47*reg61; reg45=reg73+reg45; reg73=reg19*elem.pos(6)[0];
    reg66=reg66/(*f.m).elastic_modulus; reg77=reg76-reg77; reg67=reg67-reg69; reg78=reg43-reg78; reg43=reg25*reg74;
    reg76=reg47*reg75; T reg79=pow(reg66,2); T reg80=reg30*elem.pos(7)[0]; reg53=reg56+reg53; reg73=reg45+reg73;
    reg45=reg67*reg77; reg66=reg66*reg79; reg56=reg53*reg78; T reg81=(*f.m).poisson_ratio/(*f.m).elastic_modulus; reg76=reg43-reg76;
    reg80=reg73+reg80; reg43=1.0/(*f.m).elastic_modulus; reg73=reg74*reg80; T reg82=reg53*reg24; T reg83=reg80*reg76;
    reg56=reg45-reg56; reg45=reg81*reg79; T reg84=reg43*reg66; reg24=reg67*reg24; reg66=reg81*reg66;
    T reg85=reg47*reg80; reg79=reg43*reg79; reg47=reg47*reg53; reg74=reg67*reg74; T reg86=reg25*reg80;
    reg85=reg24-reg85; reg24=reg81*reg45; T reg87=reg43*reg79; T reg88=reg81*reg84; T reg89=reg81*reg66;
    reg84=reg43*reg84; reg79=reg81*reg79; T reg90=reg67*reg61; reg83=reg56+reg83; reg61=reg53*reg61;
    reg73=reg82-reg73; reg80=reg75*reg80; reg66=reg43*reg66; reg88=reg89+reg88; reg45=reg43*reg45;
    reg84=reg84-reg89; reg87=reg87-reg24; reg79=reg24+reg79; reg77=reg77/reg83; reg73=reg73/reg83;
    reg85=reg85/reg83; reg53=reg25*reg53; reg47=reg74-reg47; reg75=reg67*reg75; reg80=reg61-reg80;
    reg86=reg90-reg86; reg78=reg78/reg83; reg25=reg54*reg85; reg66=reg89+reg66; reg87=reg43*reg87;
    reg56=reg54*reg78; reg79=reg81*reg79; reg61=reg44*reg73; reg67=reg24+reg45; reg74=reg7*reg85;
    reg82=reg11*reg77; reg89=reg44*reg77; reg90=reg7*reg78; T reg91=reg11*reg73; T reg92=reg81*reg88;
    reg53=reg75-reg53; reg86=reg86/reg83; reg76=reg76/reg83; reg43=reg43*reg84; reg47=reg47/reg83;
    reg80=reg80/reg83; reg75=reg30*reg76; T reg93=reg90+reg82; T reg94=reg54*reg86; T reg95=reg30*reg47;
    T reg96=reg74+reg91; T reg97=reg33*reg85; reg53=reg53/reg83; T reg98=reg7*reg86; T reg99=reg4*reg80;
    T reg100=reg33*reg78; T reg101=reg44*reg80; T reg102=reg5*reg76; T reg103=reg6*reg78; T reg104=reg4*reg77;
    T reg105=reg6*reg85; T reg106=reg5*reg47; reg92=reg43-reg92; reg43=reg25+reg61; T reg107=reg81*reg66;
    T reg108=reg6*reg86; T reg109=reg57*reg73; T reg110=reg57*reg80; reg79=reg87-reg79; reg67=reg81*reg67;
    reg81=reg56+reg89; reg87=reg4*reg73; T reg111=reg33*reg86; T reg112=reg57*reg77; T reg113=reg11*reg80;
    T reg114=reg3*reg47; T reg115=reg111-reg101; T reg116=reg87-reg74; T reg117=reg98+reg113; T reg118=reg100-reg89;
    T reg119=reg43+reg106; T reg120=reg113-reg108; T reg121=reg61-reg97; T reg122=reg94+reg101; reg96=reg96+reg95;
    reg81=reg81+reg102; T reg123=reg30*reg53; T reg124=reg93+reg75; reg107=reg92-reg107; reg92=reg25-reg109;
    T reg125=reg5*reg53; T reg126=reg112-reg56; T reg127=reg108+reg99; reg67=reg79-reg67; reg79=reg110-reg94;
    T reg128=reg100+reg112; T reg129=reg97+reg109; T reg130=reg111+reg110; T reg131=reg19*reg47; T reg132=reg82-reg103;
    T reg133=reg19*reg76; T reg134=reg87+reg105; T reg135=reg19*reg53; T reg136=reg105-reg91; T reg137=reg104+reg103;
    T reg138=reg98-reg99; T reg139=reg3*reg53; T reg140=reg90-reg104; T reg141=reg3*reg76; T reg142=reg123+reg117;
    reg92=reg92-reg131; reg116=reg116+reg114; reg115=reg115+reg139; reg126=reg126+reg133; T reg143=0.5*reg119;
    reg118=reg118+reg141; reg120=reg120-reg135; T reg144=0.5*reg81; reg140=reg140-reg141; reg79=reg79+reg135;
    reg121=reg121-reg114; reg138=reg138-reg139; reg136=reg136+reg131; reg122=reg125+reg122; reg132=reg132-reg133;
    T reg145=0.5*reg96; T reg146=reg106-reg134; T reg147=reg123-reg130; T reg148=reg75-reg128; reg129=reg129-reg95;
    reg67=reg67/reg107; reg137=reg137-reg102; T reg149=0.5*reg124; reg127=reg127-reg125; T reg150=0.5*reg122;
    T reg151=reg67*reg145; T reg152=0.5*reg121; T reg153=0.5*reg140; T reg154=0.5*reg127; T reg155=0.5*reg116;
    T reg156=0.5*reg146; T reg157=0.5*reg138; T reg158=0.5*reg137; T reg159=0.5*reg148; T reg160=0.5*reg136;
    T reg161=0.5*reg147; T reg162=reg67*reg149; T reg163=0.5*reg126; T reg164=0.5*reg79; T reg165=0.5*reg129;
    reg84=reg84/reg107; T reg166=0.5*reg92; T reg167=reg67*reg143; T reg168=0.5*reg118; T reg169=0.5*reg120;
    T reg170=0.5*reg132; T reg171=0.5*reg115; T reg172=0.5*reg142; T reg173=reg67*reg144; T reg174=reg67*reg160;
    T reg175=reg84*reg122; T reg176=reg67*reg150; T reg177=reg67*reg153; T reg178=reg67*reg169; T reg179=reg67*reg170;
    T reg180=reg67*reg158; T reg181=reg67*reg156; T reg182=reg67*reg163; T reg183=reg67*reg166; T reg184=reg67*reg172;
    T reg185=reg67*reg171; T reg186=reg67*reg152; T reg187=reg67*reg159; T reg188=reg67*reg164; T reg189=reg84*reg119;
    T reg190=reg154*reg67; T reg191=reg67*reg161; T reg192=reg67*reg165; reg66=reg66/reg107; T reg193=reg67*reg157;
    reg107=reg88/reg107; reg88=reg67*reg168; T reg194=2*reg162; T reg195=reg84*reg81; T reg196=reg67*reg155;
    T reg197=reg84*reg124; T reg198=reg84*reg96; reg151=2*reg151; T reg199=2*reg167; T reg200=reg84*reg142;
    reg173=2*reg173; T reg201=reg66*reg122; T reg202=reg119*reg198; T reg203=reg107*reg81; T reg204=reg84*reg121;
    T reg205=reg107*reg118; T reg206=reg107*reg124; T reg207=reg66*reg119; T reg208=reg84*reg116; T reg209=reg96*reg189;
    T reg210=reg149*reg173; T reg211=reg107*reg96; reg193=2*reg193; T reg212=reg107*reg140; T reg213=reg66*reg147;
    reg88=2*reg88; T reg214=reg84*reg118; reg186=2*reg186; T reg215=reg122*reg200; reg185=2*reg185;
    T reg216=reg107*reg119; T reg217=reg66*reg142; T reg218=reg84*reg115; T reg219=reg84*reg147; T reg220=reg84*reg79;
    T reg221=reg66*reg96; reg176=2*reg176; T reg222=reg142*reg175; T reg223=reg66*reg115; T reg224=reg84*reg120;
    T reg225=reg127*reg84; T reg226=reg127*reg66; T reg227=reg84*reg138; T reg228=reg81*reg197; T reg229=reg143*reg151;
    reg180=2*reg180; T reg230=reg84*reg129; T reg231=reg137*reg84; reg181=2*reg181; T reg232=reg107*reg148;
    reg190=2*reg190; reg182=2*reg182; T reg233=reg84*reg92; T reg234=reg84*reg126; reg183=2*reg183;
    T reg235=reg107*reg126; T reg236=reg145*reg199; T reg237=reg144*reg194; T reg238=reg124*reg195; reg188=2*reg188;
    reg179=2*reg179; T reg239=reg84*reg140; T reg240=reg66*reg79; T reg241=reg84*reg146; T reg242=reg107*reg132;
    reg177=2*reg177; T reg243=2*reg184; reg187=2*reg187; T reg244=reg84*reg148; T reg245=reg84*reg132;
    reg192=2*reg192; reg174=2*reg174; T reg246=reg66*reg120; reg178=2*reg178; T reg247=reg66*reg138;
    reg196=2*reg196; T reg248=reg137*reg107; T reg249=reg84*reg136; reg191=2*reg191; T reg250=reg145*reg174;
    T reg251=reg124*reg245; T reg252=reg124*reg246; T reg253=reg79*reg227; T reg254=reg180*reg172; T reg255=reg124*reg197;
    T reg256=reg145*reg194; T reg257=reg124*reg211; T reg258=reg145*reg186; T reg259=reg124*reg214; T reg260=reg129*reg230;
    T reg261=reg124*reg223; T reg262=reg88*reg172; T reg263=reg159*reg187; T reg264=reg129*reg233; T reg265=reg170*reg182;
    T reg266=reg136*reg233; T reg267=reg170*reg187; T reg268=reg136*reg230; T reg269=reg120*reg227; T reg270=reg120*reg225;
    T reg271=reg120*reg224; T reg272=reg120*reg206; T reg273=reg170*reg243; T reg274=reg120*reg200; T reg275=reg147*reg225;
    T reg276=reg120*reg218; T reg277=reg120*reg175; T reg278=reg120*reg220; T reg279=reg120*reg219; T reg280=reg196*reg145;
    T reg281=reg147*reg227; T reg282=reg124*reg239; T reg283=reg124*reg247; T reg284=reg177*reg172; T reg285=reg145*reg181;
    T reg286=reg124*reg231; T reg287=reg124*reg226; T reg288=reg249*reg96; T reg289=reg159*reg179; T reg290=reg179*reg149;
    T reg291=reg96*reg206; T reg292=reg149*reg151; T reg293=reg96*reg198; T reg294=reg149*reg194; T reg295=reg129*reg241;
    T reg296=reg96*reg217; T reg297=reg172*reg151; T reg298=reg96*reg204; T reg299=reg88*reg149; T reg300=reg159*reg180;
    T reg301=reg129*reg208; T reg302=reg209+reg210; T reg303=reg96*reg233; T reg304=reg182*reg149; T reg305=reg96*reg230;
    T reg306=reg159*reg177; T reg307=reg187*reg149; T reg308=reg142*reg212; T reg309=reg193*reg149; reg238=reg236+reg238;
    T reg310=reg172*reg176; T reg311=reg124*reg201; T reg312=reg159*reg182; T reg313=reg129*reg189; T reg314=reg159*reg173;
    T reg315=reg129*reg204; T reg316=reg159*reg88; T reg317=reg129*reg198; T reg318=reg159*reg194; T reg319=reg129*reg249;
    T reg320=reg172*reg173; T reg321=reg79*reg200; T reg322=reg124*reg240; T reg323=reg182*reg172; T reg324=reg145*reg192;
    T reg325=reg124*reg244; T reg326=reg124*reg213; T reg327=reg187*reg172; T reg328=reg96*reg208; T reg329=reg177*reg149;
    T reg330=reg241*reg96; T reg331=reg180*reg149; T reg332=reg199*reg156; T reg333=reg137*reg234; T reg334=reg183*reg156;
    T reg335=reg137*reg244; T reg336=reg147*reg200; T reg337=reg192*reg156; T reg338=reg179*reg172; T reg339=reg145*reg151;
    T reg340=reg177*reg158; T reg341=reg208*reg146; T reg342=reg158*reg180; T reg343=reg241*reg146; T reg344=reg158*reg179;
    T reg345=reg249*reg146; T reg346=reg145*reg183; T reg347=reg124*reg234; T reg348=reg158*reg88; T reg349=reg204*reg146;
    T reg350=reg158*reg173; T reg351=reg189*reg146; T reg352=reg158*reg182; T reg353=reg159*reg243; T reg354=reg147*reg219;
    T reg355=reg147*reg220; T reg356=reg138*reg175; T reg357=reg66*reg92; T reg358=reg138*reg220; T reg359=reg66*reg129;
    T reg360=reg147*reg175; T reg361=reg138*reg219; T reg362=reg137*reg239; T reg363=reg196*reg156; T reg364=reg147*reg218;
    T reg365=reg137*reg231; T reg366=reg181*reg156; T reg367=reg137*reg245; T reg368=reg174*reg156; T reg369=reg137*reg197;
    T reg370=reg151*reg156; T reg371=reg137*reg217; T reg372=reg154*reg194; T reg373=reg137*reg214; T reg374=reg186*reg156;
    T reg375=reg137*reg195; T reg376=reg217*reg132; T reg377=reg186*reg160; T reg378=reg214*reg132; T reg379=reg160*reg199;
    T reg380=reg132*reg195; T reg381=reg183*reg160; T reg382=reg234*reg132; T reg383=reg192*reg160; T reg384=reg244*reg132;
    T reg385=reg177*reg170; T reg386=reg147*reg224; T reg387=reg136*reg208; T reg388=reg170*reg180; T reg389=reg241*reg136;
    T reg390=reg170*reg179; T reg391=reg249*reg136; T reg392=reg170*reg194; T reg393=reg136*reg198; T reg394=reg170*reg88;
    T reg395=reg136*reg204; T reg396=reg170*reg173; T reg397=reg136*reg189; T reg398=reg146*reg233; T reg399=reg158*reg187;
    T reg400=reg146*reg230; T reg401=reg127*reg227; T reg402=reg127*reg225; T reg403=reg127*reg224; T reg404=reg127*reg206;
    T reg405=reg158*reg243; T reg406=reg127*reg200; T reg407=reg127*reg218; T reg408=reg127*reg175; T reg409=reg127*reg220;
    T reg410=reg127*reg219; T reg411=reg196*reg160; T reg412=reg239*reg132; T reg413=reg181*reg160; T reg414=reg147*reg206;
    T reg415=reg231*reg132; T reg416=reg174*reg160; T reg417=reg245*reg132; T reg418=reg160*reg151; T reg419=reg132*reg197;
    T reg420=reg169*reg194; T reg421=reg81*reg195; T reg422=reg143*reg199; T reg423=reg163*reg88; T reg424=reg81*reg216;
    T reg425=reg92*reg198; T reg426=reg163*reg194; T reg427=reg143*reg173; T reg428=reg249*reg92; T reg429=reg81*reg234;
    T reg430=reg143*reg183; T reg431=reg81*reg244; T reg432=reg143*reg192; T reg433=reg179*reg163; T reg434=reg177*reg144;
    T reg435=reg119*reg208; T reg436=reg144*reg180; T reg437=reg241*reg119; T reg438=reg179*reg144; T reg439=reg241*reg92;
    T reg440=reg249*reg119; T reg441=reg163*reg180; T reg442=reg92*reg208; T reg443=reg243*reg163; T reg444=reg79*reg206;
    T reg445=reg79*reg224; T reg446=reg79*reg225; T reg447=reg92*reg230; T reg448=reg163*reg187; T reg449=reg92*reg233;
    T reg450=reg196*reg143; T reg451=reg81*reg231; T reg452=reg143*reg181; T reg453=reg81*reg245; T reg454=reg163*reg182;
    T reg455=reg143*reg174; T reg456=reg92*reg189; T reg457=reg228+reg229; T reg458=reg150*reg243; T reg459=reg81*reg217;
    T reg460=reg150*reg194; T reg461=reg163*reg173; T reg462=reg92*reg204; T reg463=reg81*reg214; T reg464=reg143*reg186;
    T reg465=reg119*reg233; T reg466=reg187*reg144; T reg467=reg166*reg151; T reg468=reg245*reg126; T reg469=reg119*reg230;
    T reg470=reg122*reg227; T reg471=reg122*reg225; T reg472=reg122*reg224; T reg473=reg122*reg206; T reg474=reg243*reg144;
    T reg475=reg174*reg166; T reg476=reg126*reg231; T reg477=reg166*reg181; T reg478=reg239*reg126; T reg479=reg196*reg166;
    T reg480=reg237+reg215; T reg481=reg122*reg218; T reg482=reg122*reg207; T reg483=reg143*reg176; T reg484=reg122*reg175;
    T reg485=reg122*reg219; T reg486=reg122*reg220; T reg487=reg177*reg163; T reg488=reg244*reg126; T reg489=reg192*reg166;
    reg202=reg237+reg202; T reg490=reg144*reg88; T reg491=reg119*reg204; T reg492=reg144*reg199; T reg493=reg119*reg203;
    T reg494=reg144*reg173; T reg495=reg119*reg189; T reg496=reg150*reg199; T reg497=reg126*reg234; T reg498=reg166*reg183;
    T reg499=reg126*reg195; T reg500=reg166*reg199; T reg501=reg214*reg126; T reg502=reg186*reg166; T reg503=reg164*reg194;
    T reg504=reg217*reg126; T reg505=reg126*reg197; T reg506=reg119*reg201; T reg507=reg144*reg182; T reg508=reg161*reg194;
    T reg509=reg148*reg217; T reg510=reg149*reg176; T reg511=reg148*reg197; reg222=reg210+reg222; reg210=reg142*reg235;
    T reg512=reg165*reg151; T reg513=reg188*reg149; T reg514=reg142*reg220; T reg515=reg142*reg232; T reg516=reg191*reg149;
    T reg517=reg142*reg219; T reg518=reg196*reg152; T reg519=reg239*reg118; T reg520=reg152*reg181; T reg521=reg118*reg231;
    T reg522=reg152*reg174; T reg523=reg118*reg245; T reg524=reg152*reg151; T reg525=reg118*reg197; T reg526=reg171*reg194;
    T reg527=reg148*reg245; T reg528=reg142*reg227; T reg529=reg248*reg142; T reg530=reg190*reg149; T reg531=reg142*reg225;
    T reg532=reg148*reg244; T reg533=reg242*reg142; T reg534=reg165*reg192; T reg535=reg148*reg234; T reg536=reg178*reg149;
    T reg537=reg165*reg183; T reg538=reg142*reg224; T reg539=reg148*reg195; T reg540=reg165*reg199; T reg541=reg145*reg243;
    T reg542=reg142*reg221; T reg543=reg142*reg200; T reg544=reg142*reg205; T reg545=reg185*reg149; T reg546=reg148*reg214;
    T reg547=reg142*reg218; T reg548=reg165*reg186; T reg549=reg142*reg203; T reg550=reg121*reg233; T reg551=reg148*reg231;
    T reg552=reg168*reg187; T reg553=reg121*reg230; T reg554=reg115*reg227; T reg555=reg115*reg225; T reg556=reg115*reg224;
    T reg557=reg165*reg181; T reg558=reg115*reg206; T reg559=reg168*reg243; T reg560=reg115*reg200; T reg561=reg115*reg218;
    T reg562=reg115*reg175; T reg563=reg148*reg239; T reg564=reg165*reg196; T reg565=reg79*reg219; T reg566=reg79*reg220;
    reg175=reg79*reg175; reg220=reg115*reg220; reg219=reg115*reg219; T reg567=reg239*reg81; T reg568=reg79*reg218;
    T reg569=reg118*reg217; T reg570=reg152*reg186; T reg571=reg118*reg214; T reg572=reg152*reg199; T reg573=reg118*reg195;
    T reg574=reg152*reg183; T reg575=reg118*reg234; T reg576=reg152*reg192; T reg577=reg118*reg244; T reg578=reg177*reg168;
    T reg579=reg121*reg208; T reg580=reg165*reg174; T reg581=reg168*reg180; T reg582=reg241*reg121; T reg583=reg168*reg179;
    T reg584=reg249*reg121; T reg585=reg168*reg194; T reg586=reg121*reg198; T reg587=reg168*reg88; T reg588=reg121*reg204;
    T reg589=reg168*reg173; T reg590=reg121*reg189; T reg591=reg168*reg182; T reg592=reg155*reg151; reg244=reg140*reg244;
    T reg593=reg66*reg116; T reg594=reg155*reg192; T reg595=reg158*reg194; T reg596=reg198*reg146; reg204=reg116*reg204;
    T reg597=reg107*reg121; T reg598=reg138*reg200; T reg599=reg66*reg136; reg195=reg140*reg195; T reg600=reg153*reg182;
    T reg601=reg155*reg199; T reg602=reg155*reg174; reg241=reg241*reg116; reg233=reg116*reg233; reg218=reg138*reg218;
    reg245=reg140*reg245; T reg603=reg153*reg180; T reg604=reg196*reg155; T reg605=reg66*reg146; reg239=reg239*reg140;
    reg224=reg138*reg224; reg208=reg116*reg208; T reg606=reg157*reg194; T reg607=reg140*reg217; T reg608=reg177*reg153;
    T reg609=reg153*reg194; reg198=reg116*reg198; reg249=reg249*reg116; T reg610=reg107*reg129; T reg611=reg107*reg116;
    T reg612=reg155*reg186; reg214=reg140*reg214; T reg613=reg153*reg179; T reg614=reg66*reg121; T reg615=reg153*reg88;
    T reg616=reg140*reg197; reg225=reg138*reg225; T reg617=reg107*reg92; reg227=reg138*reg227; T reg618=reg116*reg189;
    T reg619=reg107*reg146; T reg620=reg107*reg136; reg231=reg140*reg231; reg234=reg140*reg234; T reg621=reg155*reg183;
    T reg622=reg153*reg187; T reg623=reg138*reg206; T reg624=reg155*reg181; T reg625=reg153*reg243; reg230=reg116*reg230;
    T reg626=reg153*reg173; T reg627=reg168*reg186; T reg628=reg121*reg217; T reg629=reg171*reg151; reg586=reg586-reg585;
    T reg630=reg153*reg199; T reg631=reg242*reg121; T reg632=reg526+reg569; T reg633=reg247*reg116; reg584=reg583+reg584;
    T reg634=reg171*reg174; T reg635=reg121*reg246; T reg636=reg121*reg206; T reg637=reg168*reg151; T reg638=reg171*reg179;
    T reg639=reg118*reg246; T reg640=reg168*reg183; T reg641=reg121*reg201; T reg642=reg171*reg199; T reg643=reg589-reg590;
    T reg644=reg121*reg203; T reg645=reg171*reg243; T reg646=reg168*reg199; T reg647=reg121*reg223; T reg648=reg524-reg525;
    T reg649=reg152*reg194; T reg650=reg171*reg186; reg588=reg587+reg588; T reg651=reg196*reg157; T reg652=reg121*reg205;
    T reg653=reg118*reg211; T reg654=reg116*reg203; T reg655=reg152*reg173; T reg656=reg121*reg212; T reg657=reg118*reg216;
    T reg658=reg196*reg168; T reg659=reg118*reg213; T reg660=reg171*reg173; T reg661=reg171*reg187; T reg662=reg118*reg610;
    T reg663=reg152*reg187; T reg664=reg118*reg201; T reg665=reg171*reg188; reg577=reg576+reg577; T reg666=reg157*reg186;
    reg575=reg574+reg575; T reg667=reg171*reg191; T reg668=reg118*reg240; T reg669=reg171*reg182; T reg670=reg118*reg617;
    T reg671=reg152*reg182; T reg672=reg168*reg174; T reg673=reg121*reg226; T reg674=reg171*reg181; T reg675=reg171*reg185;
    reg582=reg581+reg582; reg571=reg570+reg571; T reg676=reg248*reg121; T reg677=reg152*reg88; T reg678=reg168*reg181;
    T reg679=reg118*reg597; T reg680=reg247*reg121; T reg681=reg171*reg88; T reg682=reg196*reg171; T reg683=reg118*reg223;
    reg579=reg578+reg579; T reg684=reg171*reg176; T reg685=reg116*reg223; reg204=reg615+reg204; reg573=reg573-reg572;
    reg453=reg453-reg455; T reg686=reg157*reg180; T reg687=reg140*reg226; T reg688=reg150*reg180; T reg689=reg81*reg226;
    T reg690=reg143*reg180; T reg691=reg81*reg619; T reg692=reg150*reg190; reg451=reg451-reg452; T reg693=reg177*reg150;
    T reg694=reg247*reg81; T reg695=reg177*reg143; T reg696=reg611*reg81; T reg697=reg150*reg193; reg567=reg567-reg450;
    T reg698=reg116*reg246; reg219=reg552+reg219; T reg699=reg115*reg359; T reg700=reg152*reg191; T reg701=reg168*reg191;
    T reg702=reg115*reg232; reg220=reg591+reg220; T reg703=reg81*reg201; reg427=reg424+reg427; T reg704=reg157*reg190;
    T reg705=reg150*reg176; reg421=reg421+reg422; reg231=reg231+reg624; T reg706=reg150*reg88; T reg707=reg81*reg223;
    T reg708=reg143*reg88; T reg709=reg81*reg597; T reg710=reg150*reg185; reg463=reg463-reg464; T reg711=reg459+reg460;
    T reg712=reg155*reg180; T reg713=reg143*reg194; T reg714=reg81*reg211; T reg715=reg457+reg458; T reg716=reg140*reg619;
    T reg717=reg150*reg179; T reg718=reg81*reg246; T reg719=reg143*reg179; T reg720=reg81*reg620; T reg721=reg150*reg178;
    T reg722=reg152*reg178; T reg723=reg168*reg178; T reg724=reg242*reg115; reg555=reg581+reg555; reg581=reg115*reg605;
    T reg725=reg152*reg190; T reg726=reg168*reg190; T reg727=reg248*reg115; reg554=reg578+reg554; reg578=reg115*reg593;
    T reg728=reg152*reg193; T reg729=reg168*reg193; T reg730=reg115*reg212; T reg731=reg121*reg213; T reg732=reg171*reg192;
    reg553=reg552+reg553; reg552=reg121*reg232; T reg733=reg168*reg192; T reg734=reg121*reg240; T reg735=reg171*reg183;
    reg550=reg591+reg550; reg208=reg608+reg208; reg591=reg121*reg235; T reg736=reg115*reg357; T reg737=reg152*reg188;
    T reg738=reg168*reg188; T reg739=reg115*reg235; reg562=reg589+reg562; reg589=reg140*reg213; T reg740=reg115*reg207;
    T reg741=reg152*reg176; T reg742=reg168*reg176; T reg743=reg115*reg203; reg561=reg587+reg561; reg587=reg115*reg614;
    T reg744=reg152*reg185; T reg745=reg168*reg185; T reg746=reg115*reg205; T reg747=reg585+reg560; T reg748=reg115*reg221;
    T reg749=reg152*reg243; T reg750=reg558+reg559; T reg751=reg196*reg153; reg556=reg583+reg556; reg583=reg116*reg212;
    T reg752=reg115*reg599; T reg753=reg163*reg185; T reg754=reg79*reg205; T reg755=reg321+reg426; T reg756=reg79*reg221;
    T reg757=reg243*reg166; reg320=reg311+reg320; T reg758=reg116*reg201; T reg759=reg124*reg216; T reg760=reg145*reg173;
    T reg761=reg238+reg310; reg262=reg261+reg262; reg596=reg596-reg595; T reg762=reg124*reg597; T reg763=reg145*reg88;
    T reg764=reg185*reg172; reg259=reg258-reg259; T reg765=reg116*reg235; T reg766=reg172*reg194; T reg767=reg124*reg217;
    reg257=reg256+reg257; T reg768=reg610*reg126; T reg769=reg243*reg172; T reg770=reg242*reg96; T reg771=reg181*reg172;
    T reg772=reg96*reg226; reg330=reg330-reg331; reg198=reg198-reg609; T reg773=reg181*reg149; T reg774=reg248*reg96;
    T reg775=reg196*reg172; T reg776=reg96*reg247; reg328=reg328-reg329; T reg777=reg196*reg149; T reg778=reg96*reg212;
    reg327=reg326+reg327; T reg779=reg157*reg151; T reg780=reg116*reg217; T reg781=reg124*reg610; T reg782=reg145*reg187;
    T reg783=reg191*reg172; reg325=reg324-reg325; reg323=reg322+reg323; T reg784=reg153*reg186; T reg785=reg124*reg617;
    T reg786=reg145*reg182; reg282=reg280-reg282; T reg787=reg153*reg192; T reg788=reg116*reg232; reg279=reg267+reg279;
    T reg789=reg120*reg359; T reg790=reg191*reg160; T reg791=reg170*reg191; T reg792=reg120*reg232; reg278=reg265+reg278;
    T reg793=reg120*reg357; T reg794=reg188*reg160; T reg795=reg170*reg188; T reg796=reg120*reg235; reg277=reg396+reg277;
    T reg797=reg120*reg207; T reg798=reg160*reg176; T reg799=reg170*reg176; T reg800=reg120*reg203; reg276=reg394+reg276;
    T reg801=reg120*reg614; T reg802=reg185*reg160; T reg803=reg170*reg185; T reg804=reg120*reg205; T reg805=reg339+reg255;
    T reg806=reg163*reg190; T reg807=reg248*reg79; reg253=reg253+reg487; T reg808=reg79*reg593; T reg809=reg166*reg193;
    T reg810=reg163*reg193; reg233=reg600+reg233; T reg811=reg124*reg620; T reg812=reg145*reg179; T reg813=reg178*reg172;
    reg251=reg250-reg251; reg254=reg287+reg254; T reg814=reg157*reg183; T reg815=reg124*reg619; T reg816=reg145*reg180;
    T reg817=reg190*reg172; reg286=reg285-reg286; T reg818=reg116*reg240; reg284=reg283+reg284; T reg819=reg124*reg611;
    T reg820=reg145*reg177; T reg821=reg193*reg172; T reg822=reg142*reg357; T reg823=reg145*reg188; reg513=reg210+reg513;
    reg210=reg157*reg88; T reg824=reg140*reg223; reg222=reg236+reg222; T reg825=reg142*reg207; T reg826=reg145*reg176;
    reg510=reg549+reg510; reg549=reg157*reg176; reg547=reg299+reg547; T reg827=reg142*reg614; T reg828=reg145*reg185;
    reg545=reg544+reg545; reg195=reg195-reg601; reg544=reg294+reg543; reg542=reg541+reg542; T reg829=reg243*reg149;
    T reg830=reg142*reg206; reg538=reg290+reg538; T reg831=reg155*reg173; T reg832=reg142*reg599; T reg833=reg145*reg178;
    T reg834=reg118*reg620; T reg835=reg152*reg179; reg523=reg522+reg523; T reg836=reg171*reg178; T reg837=reg118*reg226;
    T reg838=reg171*reg180; T reg839=reg118*reg619; T reg840=reg152*reg180; reg521=reg520+reg521; T reg841=reg626-reg618;
    T reg842=reg171*reg190; T reg843=reg247*reg118; T reg844=reg177*reg171; T reg845=reg611*reg118; T reg846=reg177*reg152;
    reg519=reg518+reg519; T reg847=reg171*reg193; reg517=reg307+reg517; T reg848=reg142*reg359; T reg849=reg145*reg191;
    reg516=reg515+reg516; reg515=reg157*reg199; reg514=reg304+reg514; T reg850=reg96*reg235; T reg851=reg172*reg199;
    T reg852=reg96*reg201; reg310=reg310+reg302; T reg853=reg157*reg174; T reg854=reg149*reg199; T reg855=reg96*reg203;
    T reg856=reg186*reg172; T reg857=reg96*reg223; reg299=reg298-reg299; reg298=reg212*reg146; T reg858=reg186*reg149;
    T reg859=reg96*reg205; reg297=reg296+reg297; T reg860=reg153*reg151; reg293=reg293+reg294; T reg861=reg116*reg206;
    reg292=reg291+reg292; T reg862=reg88*reg164; T reg863=reg174*reg172; T reg864=reg96*reg246; reg290=reg288-reg290;
    reg288=reg174*reg149; reg536=reg533+reg536; reg533=reg157*reg243; reg531=reg331+reg531; reg331=reg142*reg605;
    T reg865=reg145*reg190; reg530=reg529+reg530; reg529=reg592-reg616; reg528=reg329+reg528; reg329=reg142*reg593;
    T reg866=reg145*reg193; reg309=reg308+reg309; reg308=reg153*reg183; T reg867=reg192*reg172; T reg868=reg96*reg213;
    reg307=reg305-reg307; reg305=reg155*reg194; T reg869=reg140*reg211; T reg870=reg192*reg149; T reg871=reg96*reg232;
    T reg872=reg183*reg172; T reg873=reg96*reg240; reg304=reg303-reg304; reg303=reg183*reg149; reg239=reg604+reg239;
    T reg874=reg148*reg211; T reg875=reg165*reg194; T reg876=reg161*reg243; T reg877=reg512-reg511; T reg878=reg161*reg179;
    T reg879=reg148*reg246; T reg880=reg148*reg620; T reg881=reg165*reg179; T reg882=reg161*reg178; reg527=reg580+reg527;
    T reg883=reg177*reg155; T reg884=reg161*reg180; T reg885=reg148*reg226; T reg886=reg148*reg619; T reg887=reg165*reg180;
    T reg888=reg161*reg190; reg551=reg557+reg551; T reg889=reg611*reg140; T reg890=reg161*reg177; T reg891=reg148*reg247;
    T reg892=reg148*reg611; T reg893=reg165*reg177; T reg894=reg161*reg191; reg532=reg534+reg532; reg249=reg613+reg249;
    T reg895=reg161*reg182; T reg896=reg148*reg240; T reg897=reg148*reg617; T reg898=reg165*reg182; T reg899=reg161*reg188;
    reg535=reg537+reg535; T reg900=reg161*reg173; T reg901=reg148*reg201; T reg902=reg148*reg216; T reg903=reg165*reg173;
    T reg904=reg161*reg176; reg539=reg539-reg540; T reg905=reg157*reg193; T reg906=reg161*reg88; T reg907=reg148*reg223;
    T reg908=reg148*reg597; T reg909=reg165*reg88; T reg910=reg161*reg185; reg546=reg548+reg546; T reg911=reg509+reg508;
    T reg912=reg157*reg178; reg445=reg433+reg445; T reg913=reg79*reg599; T reg914=reg178*reg166; T reg915=reg178*reg163;
    T reg916=reg242*reg79; reg446=reg441+reg446; T reg917=reg79*reg605; T reg918=reg166*reg190; T reg919=reg79*reg212;
    T reg920=reg92*reg213; T reg921=reg192*reg164; reg447=reg448+reg447; reg245=reg602+reg245; T reg922=reg92*reg232;
    T reg923=reg163*reg192; T reg924=reg240*reg92; T reg925=reg164*reg183; reg449=reg454+reg449; T reg926=reg92*reg235;
    T reg927=reg163*reg183; T reg928=reg92*reg201; T reg929=reg164*reg199; T reg930=reg161*reg193; reg563=reg564+reg563;
    reg565=reg448+reg565; reg448=reg177*reg157; T reg931=reg79*reg359; T reg932=reg191*reg166; T reg933=reg163*reg191;
    T reg934=reg79*reg232; reg566=reg454+reg566; reg454=reg247*reg140; T reg935=reg79*reg357; T reg936=reg166*reg188;
    T reg937=reg163*reg188; T reg938=reg79*reg235; reg175=reg461+reg175; T reg939=reg79*reg207; T reg940=reg166*reg176;
    T reg941=reg163*reg176; T reg942=reg79*reg203; reg568=reg423+reg568; T reg943=reg79*reg614; T reg944=reg185*reg166;
    T reg945=reg444+reg443; T reg946=reg153*reg181; T reg947=reg147*reg221; T reg948=reg165*reg243; T reg949=reg414+reg353;
    T reg950=reg248*reg116; reg386=reg289+reg386; T reg951=reg147*reg599; T reg952=reg165*reg178; T reg953=reg159*reg178;
    T reg954=reg147*reg242; reg275=reg300+reg275; T reg955=reg147*reg605; T reg956=reg165*reg190; T reg957=reg159*reg190;
    T reg958=reg147*reg248; reg281=reg306+reg281; T reg959=reg147*reg593; T reg960=reg165*reg193; T reg961=reg159*reg193;
    T reg962=reg147*reg212; T reg963=reg129*reg213; T reg964=reg161*reg192; reg260=reg263+reg260; reg354=reg263+reg354;
    reg263=reg147*reg359; T reg965=reg165*reg191; T reg966=reg159*reg191; T reg967=reg147*reg232; reg355=reg312+reg355;
    T reg968=reg147*reg357; T reg969=reg165*reg188; T reg970=reg159*reg188; T reg971=reg147*reg235; reg360=reg314+reg360;
    T reg972=reg116*reg205; T reg973=reg147*reg207; T reg974=reg165*reg176; T reg975=reg159*reg176; T reg976=reg147*reg203;
    reg364=reg316+reg364; T reg977=reg138*reg614; T reg978=reg147*reg614; T reg979=reg165*reg185; T reg980=reg159*reg185;
    T reg981=reg147*reg205; T reg982=reg318+reg336; T reg983=reg129*reg206; T reg984=reg159*reg151; T reg985=reg129*reg246;
    T reg986=reg161*reg174; reg319=reg289+reg319; reg289=reg153*reg174; T reg987=reg242*reg116; T reg988=reg129*reg242;
    T reg989=reg159*reg174; T reg990=reg129*reg226; T reg991=reg161*reg181; reg295=reg300+reg295; reg300=reg129*reg248;
    T reg992=reg159*reg181; T reg993=reg129*reg247; T reg994=reg161*reg196; reg301=reg306+reg301; reg306=reg129*reg212;
    T reg995=reg159*reg196; T reg996=reg161*reg187; T reg997=reg148*reg213; T reg998=reg148*reg610; T reg999=reg165*reg187;
    reg241=reg603+reg241; T reg1000=reg129*reg232; T reg1001=reg159*reg192; T reg1002=reg129*reg240; T reg1003=reg161*reg183;
    reg264=reg312+reg264; reg312=reg129*reg235; T reg1004=reg159*reg183; T reg1005=reg129*reg201; T reg1006=reg161*reg199;
    reg314=reg314-reg313; T reg1007=reg157*reg181; T reg1008=reg129*reg203; T reg1009=reg159*reg199; T reg1010=reg129*reg223;
    T reg1011=reg161*reg186; reg315=reg316+reg315; reg316=reg116*reg226; T reg1012=reg129*reg205; T reg1013=reg159*reg186;
    T reg1014=reg129*reg217; T reg1015=reg161*reg151; reg317=reg317-reg318; T reg1016=reg143*reg190; T reg1017=reg122*reg605;
    T reg1018=reg144*reg190; T reg1019=reg248*reg122; reg470=reg434+reg470; T reg1020=reg143*reg193; T reg1021=reg122*reg593;
    T reg1022=reg144*reg193; T reg1023=reg122*reg212; T reg1024=reg119*reg213; T reg1025=reg150*reg192; reg469=reg466-reg469;
    T reg1026=reg155*reg88; T reg1027=reg119*reg232; T reg1028=reg192*reg144; T reg1029=reg119*reg240; T reg1030=reg150*reg183;
    reg465=reg507-reg465; T reg1031=reg140*reg597; T reg1032=reg119*reg235; T reg1033=reg144*reg183; T reg1034=reg496+reg506;
    T reg1035=reg140*reg216; T reg1036=reg122*reg235; reg484=reg494+reg484; reg483=reg482+reg483; T reg1037=reg607+reg606;
    T reg1038=reg144*reg176; T reg1039=reg122*reg203; reg481=reg490+reg481; T reg1040=reg143*reg185; T reg1041=reg122*reg614;
    T reg1042=reg144*reg185; T reg1043=reg122*reg205; reg229=reg229+reg480; T reg1044=reg157*reg185; T reg1045=reg143*reg243;
    T reg1046=reg122*reg221; T reg1047=reg473+reg474; reg472=reg438+reg472; reg214=reg612+reg214; T reg1048=reg143*reg178;
    T reg1049=reg122*reg599; T reg1050=reg178*reg144; T reg1051=reg242*reg122; reg471=reg436+reg471; reg234=reg621+reg234;
    T reg1052=reg248*reg119; T reg1053=reg144*reg181; T reg1054=reg247*reg119; T reg1055=reg196*reg150; reg435=reg434-reg435;
    reg434=reg119*reg212; T reg1056=reg196*reg144; T reg1057=reg150*reg187; T reg1058=reg81*reg213; T reg1059=reg143*reg187;
    T reg1060=reg81*reg610; T reg1061=reg150*reg191; reg431=reg431-reg432; T reg1062=reg155*reg182; T reg1063=reg150*reg182;
    T reg1064=reg81*reg240; T reg1065=reg143*reg182; T reg1066=reg81*reg617; T reg1067=reg150*reg188; reg429=reg429-reg430;
    T reg1068=reg140*reg617; T reg1069=reg150*reg173; reg494=reg494+reg495; reg493=reg492+reg493; T reg1070=reg157*reg173;
    T reg1071=reg119*reg223; T reg1072=reg150*reg186; reg491=reg490-reg491; reg490=reg140*reg201; T reg1073=reg119*reg205;
    T reg1074=reg144*reg186; T reg1075=reg119*reg217; T reg1076=reg150*reg151; T reg1077=reg458+reg202; T reg1078=reg157*reg188;
    T reg1079=reg119*reg206; T reg1080=reg144*reg151; T reg1081=reg119*reg246; T reg1082=reg150*reg174; reg440=reg438-reg440;
    reg438=reg242*reg119; T reg1083=reg174*reg144; T reg1084=reg119*reg226; T reg1085=reg150*reg181; reg437=reg436-reg437;
    reg436=reg164*reg181; reg439=reg441+reg439; reg441=reg248*reg92; T reg1086=reg163*reg181; T reg1087=reg247*reg92;
    T reg1088=reg196*reg164; reg442=reg487+reg442; reg487=reg92*reg212; T reg1089=reg196*reg163; T reg1090=reg187*reg164;
    T reg1091=reg126*reg213; T reg1092=reg187*reg166; T reg1093=reg191*reg164; reg488=reg489+reg488; T reg1094=reg157*reg182;
    T reg1095=reg164*reg182; T reg1096=reg240*reg126; T reg1097=reg126*reg617; T reg1098=reg166*reg182; T reg1099=reg164*reg188;
    reg497=reg498+reg497; T reg1100=reg140*reg240; T reg1101=reg164*reg173; reg461=reg461-reg456; T reg1102=reg155*reg179;
    T reg1103=reg92*reg203; T reg1104=reg163*reg199; T reg1105=reg92*reg223; T reg1106=reg186*reg164; reg462=reg423+reg462;
    reg423=reg140*reg620; T reg1107=reg92*reg205; T reg1108=reg163*reg186; T reg1109=reg217*reg92; T reg1110=reg164*reg151;
    reg425=reg425-reg426; T reg1111=reg92*reg206; T reg1112=reg163*reg151; T reg1113=reg246*reg92; T reg1114=reg174*reg164;
    reg428=reg433+reg428; reg433=reg157*reg179; T reg1115=reg140*reg246; T reg1116=reg242*reg92; T reg1117=reg174*reg163;
    T reg1118=reg92*reg226; T reg1119=reg164*reg180; T reg1120=reg126*reg226; T reg1121=reg126*reg619; T reg1122=reg166*reg180;
    T reg1123=reg164*reg190; reg476=reg477+reg476; T reg1124=reg140*reg610; T reg1125=reg177*reg164; T reg1126=reg247*reg126;
    T reg1127=reg611*reg126; T reg1128=reg177*reg166; T reg1129=reg164*reg193; reg478=reg479+reg478; reg485=reg466+reg485;
    reg466=reg143*reg191; T reg1130=reg122*reg359; T reg1131=reg191*reg144; T reg1132=reg122*reg232; reg486=reg507+reg486;
    reg507=reg157*reg187; T reg1133=reg143*reg188; T reg1134=reg122*reg357; T reg1135=reg144*reg188; T reg1136=reg126*reg201;
    T reg1137=reg126*reg216; T reg1138=reg166*reg173; T reg1139=reg164*reg176; reg499=reg499-reg500; T reg1140=reg126*reg223;
    T reg1141=reg597*reg126; T reg1142=reg88*reg166; T reg1143=reg185*reg164; reg501=reg502+reg501; T reg1144=reg157*reg191;
    T reg1145=reg504+reg503; reg244=reg594+reg244; T reg1146=reg126*reg211; T reg1147=reg166*reg194; T reg1148=reg243*reg164;
    T reg1149=reg467-reg505; T reg1150=reg179*reg164; T reg1151=reg246*reg126; T reg1152=reg620*reg126; T reg1153=reg179*reg166;
    T reg1154=reg178*reg164; reg468=reg475+reg468; T reg1155=reg155*reg187; T reg1156=reg137*reg611; T reg1157=reg177*reg169;
    T reg1158=reg247*reg132; T reg1159=reg158*reg193; T reg1160=reg127*reg212; T reg1161=reg170*reg183; T reg1162=reg136*reg235;
    T reg1163=reg213*reg146; T reg1164=reg154*reg192; T reg1165=reg177*reg156; reg387=reg385+reg387; reg400=reg399+reg400;
    T reg1166=reg169*reg190; reg415=reg413+reg415; T reg1167=reg179*reg156; reg358=reg600+reg358; reg600=reg138*reg203;
    reg375=reg375-reg332; reg266=reg265+reg266; reg265=reg146*reg232; T reg1168=reg158*reg190; T reg1169=reg153*reg191;
    T reg1170=reg127*reg248; T reg1171=reg137*reg223; T reg1172=reg137*reg246; T reg1173=reg153*reg193; T reg1174=reg138*reg212;
    T reg1175=reg177*reg160; reg401=reg340+reg401; reg396=reg396-reg397; T reg1176=reg138*reg599; reg611=reg611*reg132;
    T reg1177=reg154*reg88; T reg1178=reg247*reg136; T reg1179=reg138*reg232; T reg1180=reg169*reg199; T reg1181=reg196*reg169;
    T reg1182=reg127*reg593; T reg1183=reg201*reg136; T reg1184=reg193*reg156; T reg1185=reg158*reg183; T reg1186=reg132*reg213;
    T reg1187=reg169*reg187; T reg1188=reg201*reg146; T reg1189=reg154*reg199; T reg1190=reg116*reg213; T reg1191=reg169*reg178;
    T reg1192=reg137*reg201; reg268=reg267+reg268; reg267=reg350-reg351; reg224=reg613+reg224; reg613=reg610*reg132;
    T reg1193=reg138*reg357; T reg1194=reg187*reg160; reg417=reg416+reg417; T reg1195=reg169*reg192; T reg1196=reg203*reg146;
    T reg1197=reg154*reg173; T reg1198=reg136*reg213; T reg1199=reg158*reg199; T reg1200=reg137*reg620; T reg1201=reg154*reg176;
    T reg1202=reg158*reg192; T reg1203=reg169*reg183; T reg1204=reg180*reg160; T reg1205=reg136*reg212; T reg1206=reg240*reg136;
    T reg1207=reg240*reg146; reg183=reg154*reg183; T reg1208=reg619*reg132; T reg1209=reg153*reg176; T reg1210=reg196*reg170;
    reg398=reg352+reg398; T reg1211=reg137*reg216; T reg1212=reg169*reg180; T reg1213=reg170*reg192; T reg1214=reg136*reg232;
    T reg1215=reg226*reg132; T reg1216=reg173*reg156; T reg1217=reg146*reg235; T reg1218=reg188*reg156; T reg1219=reg170*reg174;
    reg614=reg127*reg614; T reg1220=reg185*reg156; reg357=reg127*reg357; T reg1221=reg158*reg185; T reg1222=reg155*reg193;
    T reg1223=reg127*reg205; T reg1224=reg154*reg243; reg393=reg393-reg392; T reg1225=reg371+reg372; T reg1226=reg370-reg369;
    T reg1227=reg595+reg406; T reg1228=reg138*reg359; T reg1229=reg226*reg136; reg409=reg352+reg409; reg352=reg169*reg151;
    T reg1230=reg217*reg136; T reg1231=reg169*reg181; reg232=reg127*reg232; reg391=reg390+reg391; T reg1232=reg194*reg156;
    T reg1233=reg242*reg138; T reg1234=reg138*reg593; reg408=reg350+reg408; reg361=reg622+reg361; reg350=reg127*reg207;
    T reg1235=reg169*reg174; T reg1236=reg246*reg136; T reg1237=reg176*reg156; T reg1238=reg137*reg211; T reg1239=reg158*reg176;
    T reg1240=reg127*reg235; T reg1241=reg158*reg188; T reg1242=reg127*reg203; T reg1243=reg242*reg136; T reg1244=reg170*reg151;
    reg407=reg348+reg407; T reg1245=reg136*reg206; T reg1246=reg153*reg178; T reg1247=reg178*reg156; reg410=reg399+reg410;
    reg399=reg154*reg193; T reg1248=reg158*reg178; T reg1249=reg127*reg242; T reg1250=reg169*reg186; T reg1251=reg136*reg223;
    T reg1252=reg137*reg597; T reg1253=reg88*reg156; T reg1254=reg169*reg193; T reg1255=reg154*reg179; reg402=reg342+reg402;
    T reg1256=reg248*reg136; T reg1257=reg170*reg181; T reg1258=reg170*reg199; T reg1259=reg127*reg605; reg225=reg603+reg225;
    reg203=reg136*reg203; reg603=reg190*reg156; reg412=reg411+reg412; T reg1260=reg158*reg191; T reg1261=reg127*reg221;
    T reg1262=reg243*reg156; T reg1263=reg153*reg185; reg218=reg615+reg218; reg615=reg404+reg405; T reg1264=reg155*reg185;
    T reg1265=reg138*reg205; T reg1266=reg170*reg186; T reg1267=reg136*reg205; reg373=reg373+reg374; reg362=reg362+reg363;
    reg389=reg388+reg389; reg403=reg344+reg403; T reg1268=reg155*reg178; T reg1269=reg191*reg156; T reg1270=reg154*reg185;
    T reg1271=reg155*reg191; reg359=reg127*reg359; T reg1272=reg127*reg599; reg395=reg394+reg395; reg270=reg388+reg270;
    reg388=reg617*reg132; reg345=reg344+reg345; reg344=reg169*reg173; reg341=reg340+reg341; reg340=reg182*reg156;
    reg394=reg154*reg174; reg211=reg132*reg211; T reg1273=reg160*reg194; T reg1274=reg154*reg190; T reg1275=reg246*reg146;
    T reg1276=reg155*reg243; reg230=reg622+reg230; reg622=reg154*reg180; reg617=reg137*reg617; T reg1277=reg158*reg151;
    T reg1278=reg120*reg605; T reg1279=reg272+reg273; reg619=reg137*reg619; T reg1280=reg190*reg160; T reg1281=reg206*reg146;
    T reg1282=reg154*reg191; T reg1283=reg138*reg221; T reg1284=reg169*reg182; T reg1285=reg138*reg207; T reg1286=reg169*reg88;
    T reg1287=reg418-reg419; reg235=reg138*reg235; reg347=reg346-reg347; T reg1288=reg170*reg190; T reg1289=reg188*reg172;
    T reg1290=reg132*reg223; reg378=reg377+reg378; T reg1291=reg169*reg188; reg343=reg342+reg343; reg599=reg120*reg599;
    reg342=reg178*reg160; reg382=reg381+reg382; reg88=reg88*reg160; T reg1292=reg154*reg181; T reg1293=reg226*reg146;
    T reg1294=reg154*reg182; reg180=reg180*reg156; T reg1295=reg248*reg146; reg181=reg158*reg181; reg185=reg169*reg185;
    T reg1296=reg137*reg240; reg271=reg390+reg271; reg174=reg158*reg174; reg390=reg170*reg178; reg597=reg597*reg132;
    T reg1297=reg242*reg146; T reg1298=reg420+reg376; reg365=reg365+reg366; reg242=reg242*reg120; T reg1299=reg247*reg146;
    reg356=reg626+reg356; reg626=reg154*reg196; reg201=reg201*reg132; T reg1300=reg155*reg190; reg335=reg335+reg337;
    reg226=reg137*reg226; T reg1301=reg623+reg625; reg182=reg182*reg160; reg221=reg120*reg221; T reg1302=reg169*reg176;
    reg196=reg196*reg158; T reg1303=reg158*reg186; reg205=reg205*reg146; reg190=reg153*reg190; T reg1304=reg169*reg179;
    reg593=reg120*reg593; T reg1305=reg160*reg193; T reg1306=reg609+reg598; reg173=reg160*reg173; reg176=reg155*reg176;
    T reg1307=reg248*reg138; reg227=reg608+reg227; reg608=reg154*reg187; reg213=reg137*reg213; reg349=reg348+reg349;
    reg348=reg392+reg274; reg384=reg383+reg384; reg380=reg380-reg379; reg367=reg367+reg368; reg620=reg620*reg132;
    reg610=reg137*reg610; reg193=reg170*reg193; reg212=reg120*reg212; reg186=reg154*reg186; reg179=reg179*reg160;
    T reg1308=reg155*reg188; reg187=reg187*reg156; reg223=reg223*reg146; reg178=reg154*reg178; reg247=reg137*reg247;
    T reg1309=reg132*reg216; reg605=reg138*reg605; T reg1310=reg154*reg188; reg338=reg338+reg252; reg269=reg385+reg269;
    reg192=reg157*reg192; reg333=reg333+reg334; reg248=reg248*reg120; reg151=reg154*reg151; reg191=reg169*reg191;
    reg188=reg153*reg188; reg385=reg243*reg160; T reg1311=reg169*reg243; reg240=reg240*reg132; reg177=reg154*reg177;
    T reg1312=reg217*reg146; reg246=reg246*reg132; reg1012=reg1013+reg1012; reg998=reg999+reg998; reg717=reg718+reg717;
    reg357=reg1218+reg357; reg316=reg1007+reg316; reg1290=reg1286+reg1290; reg716=reg712+reg716; reg579=reg847+reg579;
    reg659=reg661+reg659; reg180=reg619+reg180; reg631=reg672+reg631; reg584=reg836+reg584; reg317=reg317-reg876;
    reg687=reg686+reg687; reg680=reg682+reg680; reg380=reg1302+reg380; reg693=reg694+reg693; reg1241=reg1240+reg1241;
    reg582=reg842+reg582; reg409=reg334+reg409; reg714=reg714+reg713; reg1015=reg1015-reg1014; reg688=reg689+reg688;
    reg453=reg453+reg721; reg532=reg532+reg894; reg656=reg658+reg656; reg334=reg83*reg711; reg463=reg463+reg710;
    reg895=reg896+reg895; reg690=reg691-reg690; reg597=reg88+reg597; reg719=reg720-reg719; reg1260=reg232+reg1260;
    reg897=reg898+reg897; reg676=reg678+reg676; reg1226=reg1226-reg1224; reg673=reg674+reg673; reg88=reg83*reg715;
    reg451=reg451+reg692; reg408=reg408-reg332; reg987=reg289+reg987; reg752=reg722+reg752; reg1208=reg1204+reg1208;
    reg246=reg1304+reg246; reg556=reg522+reg556; reg550=reg665+reg550; reg300=reg992+reg300; reg232=reg83*reg750;
    reg319=reg882+reg319; reg591=reg640+reg591; reg415=reg1166+reg415; reg748=reg748-reg749; reg524=reg524-reg747;
    reg745=reg746+reg745; reg641=reg641-reg642; reg1287=reg1287-reg1311; reg993=reg994+reg993; reg643=reg684+reg643;
    reg1167=reg1200+reg1167; reg587=reg744+reg587; reg729=reg730+reg729; reg417=reg1191+reg417; reg731=reg732+reg731;
    reg578=reg728+reg578; reg990=reg991+reg990; reg605=reg1300+reg605; reg554=reg518+reg554; reg367=reg367+reg178;
    reg553=reg667+reg553; reg620=reg179+reg620; reg988=reg989+reg988; reg726=reg727+reg726; reg581=reg725+reg581;
    reg583=reg751+reg583; reg1215=reg1212+reg1215; reg552=reg733+reg552; reg555=reg520+reg555; reg295=reg888+reg295;
    reg1264=reg977+reg1264; reg734=reg735+reg734; reg723=reg724+reg723; reg984=reg984-reg983; reg629=reg629-reg628;
    reg225=reg624+reg225; reg701=reg702+reg701; reg179=reg83*reg1298; reg586=reg586-reg645; reg1255=reg1172+reg1255;
    reg699=reg700+reg699; reg633=reg651+reg633; reg853=reg698+reg853; reg219=reg576+reg219; reg637=reg637-reg636;
    reg410=reg337+reg410; reg996=reg997+reg996; reg190=reg1307+reg190; reg567=reg567+reg697; reg635=reg634+reg635;
    reg695=reg696-reg695; reg359=reg1269+reg359; reg378=reg185+reg378; reg589=reg507+reg589; reg622=reg226+reg622;
    reg1158=reg1157+reg1158; reg561=reg570+reg561; reg742=reg743+reg742; reg208=reg905+reg208; reg301=reg930+reg301;
    reg644=reg644-reg646; reg741=reg741-reg740; reg985=reg986+reg985; reg611=reg1175+reg611; reg647=reg650+reg647;
    reg211=reg211-reg1273; reg588=reg675+reg588; reg562=reg562-reg572; reg738=reg739+reg738; reg412=reg1254+reg412;
    reg652=reg627+reg652; reg736=reg737+reg736; reg220=reg574+reg220; reg306=reg995+reg306; reg1146=reg1146-reg1147;
    reg930=reg563+reg930; reg592=reg592-reg1306; reg226=reg83*reg1145; reg205=reg1303+reg205; reg501=reg501+reg1143;
    reg151=reg151-reg1312; reg1141=reg1142+reg1141; reg1140=reg862+reg1140; reg565=reg489+reg565; reg499=reg499+reg1139;
    reg596=reg596-reg1224; reg1138=reg1138-reg1137; reg1101=reg1136+reg1101; reg931=reg932+reg931; reg1100=reg1094+reg1100;
    reg347=reg347-reg1289; reg497=reg497+reg1099; reg333=reg333+reg1310; reg1097=reg1098+reg1097; reg1277=reg1277-reg1281;
    reg1095=reg1096+reg1095; reg933=reg934+reg933; reg488=reg488+reg1093; reg1092=reg768+reg1092; reg1275=reg394+reg1275;
    reg1090=reg1091+reg1090; reg566=reg498+reg566; reg1217=reg1185+reg1217; reg1131=reg1132+reg1131; reg466=reg1130-reg466;
    reg1216=reg1216-reg1211; reg432=reg485-reg432; reg888=reg551+reg888; reg1124=reg1155+reg1124; reg1188=reg1188-reg1189;
    reg478=reg478+reg1129; reg1127=reg1128+reg1127; reg1125=reg1126+reg1125; reg889=reg883+reg889; reg267=reg1201+reg267;
    reg224=reg602+reg224; reg476=reg476+reg1123; reg1121=reg1122+reg1121; reg890=reg891+reg890; reg1196=reg1196-reg1199;
    reg1119=reg1120+reg1119; reg1197=reg1192+reg1197; reg468=reg468+reg1154; reg223=reg186+reg223; reg1152=reg1153+reg1152;
    reg892=reg893+reg892; reg1150=reg1151+reg1150; reg244=reg1144+reg244; reg349=reg1270+reg349; reg1149=reg1149-reg1148;
    reg462=reg1143+reg462; reg1105=reg1106+reg1105; reg1103=reg1103-reg1104; reg941=reg942+reg941; reg1299=reg626+reg1299;
    reg461=reg1139+reg461; reg928=reg928-reg929; reg341=reg399+reg341; reg926=reg927+reg926; reg568=reg502+reg568;
    reg449=reg1099+reg449; reg924=reg925+reg924; reg1283=reg1283-reg1276; reg922=reg923+reg922; reg245=reg912+reg245;
    reg943=reg944+reg943; reg447=reg1093+reg447; reg335=reg335+reg1282; reg186=reg83*reg338; reg920=reg921+reg920;
    reg289=reg83*reg945; reg917=reg918+reg917; reg608=reg213+reg608; reg446=reg477+reg446; reg915=reg916+reg915;
    reg913=reg914+reg913; reg187=reg610+reg187; reg445=reg475+reg445; reg487=reg1089+reg487; reg340=reg617+reg340;
    reg442=reg1129+reg442; reg345=reg178+reg345; reg1087=reg1088+reg1087; reg935=reg936+reg935; reg178=reg83*reg1301;
    reg441=reg1086+reg441; reg1115=reg433+reg1115; reg439=reg1123+reg439; reg1297=reg174+reg1297; reg1118=reg436+reg1118;
    reg1116=reg1117+reg1116; reg937=reg938+reg937; reg1293=reg1292+reg1293; reg428=reg1154+reg428; reg1113=reg1114+reg1113;
    reg1112=reg1112-reg1111; reg175=reg175-reg500; reg1294=reg1296+reg1294; reg343=reg1274+reg343; reg425=reg425-reg1148;
    reg454=reg448+reg454; reg1110=reg1110-reg1109; reg940=reg940-reg939; reg1107=reg1108+reg1107; reg423=reg1102+reg423;
    reg1295=reg181+reg1295; reg1055=reg1055-reg1054; reg1246=reg1233+reg1246; reg539=reg539+reg904; reg1052=reg1053-reg1052;
    reg234=reg1078+reg234; reg1261=reg1261-reg1262; reg437=reg692+reg437; reg1085=reg1085-reg1084; reg438=reg1083-reg438;
    reg174=reg83*reg615; reg440=reg721+reg440; reg906=reg907+reg906; reg1082=reg1082-reg1081; reg403=reg368+reg403;
    reg1080=reg1080+reg1079; reg181=reg83*reg1077; reg908=reg909+reg908; reg1076=reg1076+reg1075; reg1270=reg373+reg1270;
    reg1272=reg1247+reg1272; reg1073=reg1074-reg1073; reg490=reg1070+reg490; reg491=reg710+reg491; reg546=reg546+reg910;
    reg1248=reg1249+reg1248; reg1072=reg1072-reg1071; reg213=reg83*reg493; reg831=reg831-reg1035; reg708=reg709-reg708;
    reg1238=reg1238-reg1232; reg706=reg707+reg706; reg231=reg704+reg231; reg1237=reg1237-reg350; reg421=reg421+reg705;
    reg535=reg535+reg899; reg337=reg83*reg427; reg1239=reg1242+reg1239; reg1069=reg703+reg1069; reg1068=reg1062+reg1068;
    reg912=reg249+reg912; reg429=reg429+reg1067; reg407=reg374+reg407; reg1065=reg1066-reg1065; reg1063=reg1064+reg1063;
    reg900=reg901+reg900; reg614=reg1220+reg614; reg431=reg431+reg1061; reg1263=reg1265+reg1263; reg1059=reg1060-reg1059;
    reg1221=reg1223+reg1221; reg1057=reg1058+reg1057; reg903=reg903-reg902; reg434=reg1056-reg434; reg249=reg83*reg1225;
    reg435=reg697+reg435; reg370=reg370-reg1227; reg1050=reg1051+reg1050; reg878=reg879+reg878; reg1048=reg1049-reg1048;
    reg1163=reg1164+reg1163; reg455=reg472-reg455; reg368=reg83*reg1047; reg880=reg881+reg880; reg400=reg1282+reg400;
    reg1046=reg1046+reg1045; reg1176=reg1268+reg1176; reg373=reg83*reg229; reg1042=reg1043+reg1042; reg882=reg527+reg882;
    reg265=reg1202+reg265; reg1040=reg1041-reg1040; reg464=reg481-reg464; reg1201=reg375+reg1201; reg1038=reg1039+reg1038;
    reg374=reg83*reg1037; reg1207=reg183+reg1207; reg183=reg83*reg483; reg484=reg422+reg484; reg884=reg885+reg884;
    reg1135=reg1036+reg1135; reg398=reg1310+reg398; reg1133=reg1134-reg1133; reg430=reg486-reg430; reg886=reg887+reg886;
    reg494=reg705+reg494; reg1031=reg1026+reg1031; reg402=reg366+reg402; reg366=reg83*reg1034; reg375=reg83*reg911;
    reg1032=reg1033-reg1032; reg1253=reg1252+reg1253; reg1259=reg603+reg1259; reg465=reg1067+reg465; reg239=reg905+reg239;
    reg1030=reg1030-reg1029; reg1168=reg1170+reg1168; reg1027=reg1028-reg1027; reg874=reg874-reg875; reg469=reg1061+reg469;
    reg1025=reg1025-reg1024; reg401=reg363+reg401; reg1022=reg1023+reg1022; reg877=reg877-reg876; reg1020=reg1021-reg1020;
    reg450=reg470-reg450; reg1177=reg1171+reg1177; reg1182=reg1184+reg1182; reg1018=reg1019+reg1018; reg1016=reg1017-reg1016;
    reg214=reg1044+reg214; reg452=reg471-reg452; reg1159=reg1160+reg1159; reg363=reg83*reg262; reg963=reg964+reg963;
    reg827=reg828-reg827; reg1198=reg1195+reg1198; reg389=reg1166+reg389; reg762=reg763-reg762; reg547=reg258-reg547;
    reg1234=reg1222+reg1234; reg258=reg83*reg510; reg259=reg259-reg764; reg193=reg212+reg193; reg975=reg976+reg975;
    reg765=reg308+reg765; reg260=reg894+reg260; reg826=reg826+reg825; reg212=reg767+reg766; reg824=reg210+reg824;
    reg210=reg83*reg257; reg593=reg1305+reg593; reg784=reg972+reg784; reg978=reg979+reg978; reg538=reg250-reg538;
    reg756=reg756-reg757; reg1190=reg192+reg1190; reg1243=reg1219+reg1243; reg192=reg830+reg829; reg250=reg83*reg320;
    reg961=reg962+reg961; reg268=reg191+reg268; reg308=reg83*reg542; reg760=reg760+reg759; reg195=reg549+reg195;
    reg1209=reg600+reg1209; reg339=reg339+reg544; reg394=reg83*reg761; reg758=reg758-reg515; reg1229=reg1231+reg1229;
    reg433=reg83*reg545; reg364=reg548+reg364; reg387=reg1254+reg387; reg517=reg324-reg517; reg360=reg360-reg540;
    reg549=reg841+reg549; reg811=reg812-reg811; reg1165=reg1156+reg1165; reg519=reg847+reg519; reg251=reg251-reg813;
    reg1278=reg1280+reg1278; reg845=reg846+reg845; reg233=reg1078+reg233; reg264=reg899+reg264; reg1205=reg1210+reg1205;
    reg843=reg844+reg843; reg324=reg83*reg254; reg970=reg971+reg970; reg704=reg241+reg704; reg1186=reg1187+reg1186;
    reg521=reg842+reg521; reg815=reg816-reg815; reg241=reg83*reg222; reg399=reg362+reg399; reg1256=reg1257+reg1256;
    reg362=reg83*reg513; reg188=reg235+reg188; reg822=reg823-reg822; reg805=reg805+reg769; reg974=reg974-reg973;
    reg806=reg807+reg806; reg1000=reg1001+reg1000; reg1178=reg1181+reg1178; reg514=reg346-reg514; reg269=reg411+reg269;
    reg235=reg83*reg516; reg479=reg253+reg479; reg808=reg809+reg808; reg1288=reg248+reg1288; reg848=reg849-reg848;
    reg919=reg810+reg919; reg1002=reg1003+reg1002; reg1173=reg1174+reg1173; reg288=reg770-reg288; reg771=reg772-reg771;
    reg248=reg83*reg310; reg1183=reg1183-reg1180; reg253=reg83*reg949; reg330=reg330-reg817; reg1267=reg1266+reg1267;
    reg852=reg852+reg851; reg275=reg557+reg275; reg303=reg850-reg303; reg198=reg198-reg533; reg773=reg774-reg773;
    reg869=reg869-reg305; reg1289=reg304-reg1289; reg775=reg776-reg775; reg1162=reg1161+reg1162; reg358=reg621+reg358;
    reg352=reg352-reg1230; reg872=reg873-reg872; reg293=reg769+reg293; reg203=reg203-reg1258; reg860=reg860-reg861;
    reg951=reg952+reg951; reg304=reg83*reg297; reg1169=reg1179+reg1169; reg346=reg83*reg292; reg1251=reg1250+reg1251;
    reg858=reg859-reg858; reg386=reg580+reg386; reg196=reg298+reg196; reg764=reg299-reg764; reg863=reg864-reg863;
    reg813=reg290-reg813; reg856=reg857-reg856; reg396=reg1302+reg396; reg950=reg946+reg950; reg953=reg954+reg953;
    reg395=reg185+reg395; reg855=reg855+reg854; reg528=reg280-reg528; reg185=reg83*reg323; reg529=reg529-reg533;
    reg281=reg564+reg281; reg280=reg83*reg530; reg980=reg981+reg980; reg1236=reg1235+reg1236; reg331=reg865-reg331;
    reg785=reg786-reg785; reg218=reg612+reg218; reg531=reg285-reg531; reg753=reg754+reg753; reg1214=reg1213+reg1214;
    reg391=reg1191+reg391; reg285=reg83*reg536; reg467=reg467-reg755; reg959=reg960+reg959; reg832=reg833-reg832;
    reg1193=reg1308+reg1193; reg361=reg594+reg361; reg328=reg328-reg821; reg947=reg947-reg948; reg870=reg871-reg870;
    reg955=reg956+reg955; reg393=reg393-reg1311; reg307=reg307-reg783; reg777=reg778-reg777; reg1228=reg1271+reg1228;
    reg867=reg868-reg867; reg290=reg83*reg327; reg957=reg958+reg957; reg266=reg1291+reg266; reg298=reg83*reg309;
    reg512=reg512-reg982; reg781=reg782-reg781; reg1244=reg1244-reg1245; reg329=reg866-reg329; reg783=reg325-reg783;
    reg1206=reg1203+reg1206; reg779=reg779-reg780; reg798=reg798-reg797; reg279=reg383+reg279; reg653=reg653-reg649;
    reg177=reg247+reg177; reg799=reg800+reg799; reg355=reg537+reg355; reg648=reg648-reg645; reg1010=reg1011+reg1010;
    reg599=reg342+reg599; reg227=reg604+reg227; reg654=reg654-reg630; reg201=reg344+reg201; reg575=reg665+reg575;
    reg821=reg282-reg821; reg1005=reg1005-reg1006; reg670=reg671+reg670; reg639=reg638+reg639; reg384=reg191+reg384;
    reg834=reg835+reg834; reg390=reg242+reg390; reg819=reg820-reg819; reg276=reg377+reg276; reg221=reg221-reg385;
    reg668=reg669+reg668; reg685=reg666+reg685; reg793=reg794+reg793; reg230=reg1144+reg230; reg382=reg1291+reg382;
    reg966=reg967+reg966; reg683=reg681+reg683; reg278=reg381+reg278; reg573=reg684+reg573; reg679=reg677+reg679;
    reg795=reg796+reg795; reg388=reg182+reg388; reg1008=reg1008-reg1009; reg271=reg416+reg271; reg176=reg176-reg1285;
    reg791=reg792+reg791; reg571=reg675+reg571; reg277=reg277-reg379; reg655=reg655-reg657; reg182=reg83*reg1279;
    reg314=reg904+reg314; reg789=reg790+reg789; reg263=reg965+reg263; reg664=reg660+reg664; reg191=reg83*reg632;
    reg240=reg1284+reg240; reg788=reg787+reg788; reg356=reg356-reg601; reg837=reg838+reg837; reg801=reg802+reg801;
    reg803=reg804+reg803; reg354=reg534+reg354; reg315=reg910+reg315; reg523=reg836+reg523; reg418=reg418-reg348;
    reg270=reg413+reg270; reg173=reg173-reg1309; reg242=reg83*reg284; reg662=reg663+reg662; reg577=reg667+reg577;
    reg1274=reg365+reg1274; reg968=reg969+reg968; reg312=reg1004+reg312; reg818=reg814+reg818; reg817=reg286-reg817;
    reg613=reg1194+reg613; reg204=reg1044+reg204; reg839=reg840+reg839; reg277=reg83*reg277; reg781=reg83*reg781;
    reg512=reg83*reg512; reg1038=reg83*reg1038; reg464=reg83*reg464; reg882=reg83*reg882; reg798=reg83*reg798;
    reg777=reg83*reg777; reg1042=reg83*reg1042; reg187=reg83*reg187; reg265=reg83*reg265; reg1190=reg83*reg1190;
    reg203=reg83*reg203; reg196=reg83*reg196; reg247=ponderation*reg290; reg1040=reg83*reg1040; reg924=reg83*reg924;
    reg282=ponderation*reg374; reg874=reg83*reg874; reg1027=reg83*reg1027; reg454=reg83*reg454; reg266=reg83*reg266;
    reg928=reg83*reg928; reg418=reg83*reg418; reg430=reg83*reg430; reg785=reg83*reg785; reg966=reg83*reg966;
    reg886=reg83*reg886; reg980=reg83*reg980; reg1133=reg83*reg1133; reg465=reg83*reg465; reg1206=reg83*reg1206;
    reg1253=reg83*reg1253; reg926=reg83*reg926; reg1135=reg83*reg1135; reg568=reg83*reg568; reg795=reg83*reg795;
    reg445=reg83*reg445; reg398=reg83*reg398; reg286=ponderation*reg185; reg884=reg83*reg884; reg950=reg83*reg950;
    reg484=reg83*reg484; reg779=reg83*reg779; reg245=reg83*reg245; reg913=reg83*reg913; reg299=ponderation*reg183;
    reg1030=reg83*reg1030; reg783=reg83*reg783; reg1168=reg83*reg1168; reg1207=reg83*reg1207; reg449=reg83*reg449;
    reg1201=reg83*reg1201; reg325=ponderation*reg346; reg1177=reg83*reg1177; reg1022=reg83*reg1022; reg773=reg83*reg773;
    reg1050=reg83*reg1050; reg877=reg83*reg877; reg920=reg83*reg920; reg342=ponderation*reg253; reg335=reg83*reg335;
    reg878=reg83*reg878; reg452=reg83*reg452; reg396=reg83*reg396; reg1183=reg83*reg1183; reg1283=reg83*reg1283;
    reg806=reg83*reg806; reg1159=reg83*reg1159; reg344=ponderation*reg289; reg1016=reg83*reg1016; reg330=reg83*reg330;
    reg276=reg83*reg276; reg446=reg83*reg446; reg239=reg83*reg239; reg1018=reg83*reg1018; reg801=reg83*reg801;
    reg1020=reg83*reg1020; reg917=reg83*reg917; reg771=reg83*reg771; reg450=reg83*reg450; reg288=reg83*reg288;
    reg608=reg83*reg608; reg1182=reg83*reg1182; reg214=reg83*reg214; reg221=reg83*reg221; reg860=reg83*reg860;
    reg1176=reg83*reg1176; reg803=reg83*reg803; reg922=reg83*reg922; reg365=ponderation*reg373; reg943=reg83*reg943;
    reg947=reg83*reg947; reg592=reg83*reg592; reg198=reg83*reg198; reg386=reg83*reg386; reg1046=reg83*reg1046;
    reg377=ponderation*reg182; reg381=ponderation*reg186; reg469=reg83*reg469; reg400=reg83*reg400; reg328=reg83*reg328;
    reg354=reg83*reg354; reg358=reg83*reg358; reg1048=reg83*reg1048; reg813=reg83*reg813; reg1163=reg83*reg1163;
    reg775=reg83*reg775; reg1025=reg83*reg1025; reg401=reg83*reg401; reg455=reg83*reg455; reg863=reg83*reg863;
    reg915=reg83*reg915; reg1162=reg83*reg1162; reg880=reg83*reg880; reg799=reg83*reg799; reg383=ponderation*reg368;
    reg447=reg83*reg447; reg263=reg83*reg263; reg593=reg83*reg593; reg390=reg83*reg390; reg499=reg83*reg499;
    reg428=reg83*reg428; reg596=reg83*reg596; reg269=reg83*reg269; reg974=reg83*reg974; reg1138=reg83*reg1138;
    reg931=reg83*reg931; reg347=reg83*reg347; reg1293=reg83*reg1293; reg479=reg83*reg479; reg1101=reg83*reg1101;
    reg819=reg83*reg819; reg808=reg83*reg808; reg1116=reg83*reg1116; reg968=reg83*reg968; reg818=reg83*reg818;
    reg212=reg83*reg212; reg411=ponderation*reg226; reg413=ponderation*reg210; reg1112=reg83*reg1112; reg889=reg83*reg889;
    reg343=reg83*reg343; reg175=reg83*reg175; reg501=reg83*reg501; reg1113=reg83*reg1113; reg151=reg83*reg151;
    reg765=reg83*reg765; reg188=reg83*reg188; reg1141=reg83*reg1141; reg821=reg83*reg821; reg805=reg83*reg805;
    reg1140=reg83*reg1140; reg565=reg83*reg565; reg1209=reg83*reg1209; reg1100=reg83*reg1100; reg566=reg83*reg566;
    reg1278=reg83*reg1278; reg270=reg83*reg270; reg233=reg83*reg233; reg1115=reg83*reg1115; reg1090=reg83*reg1090;
    reg416=ponderation*reg324; reg340=reg83*reg340; reg441=reg83*reg441; reg356=reg83*reg356; reg487=reg83*reg487;
    reg817=reg83*reg817; reg1068=reg83*reg1068; reg345=reg83*reg345; reg970=reg83*reg970; reg442=reg83*reg442;
    reg1087=reg83*reg1087; reg935=reg83*reg935; reg815=reg83*reg815; reg919=reg83*reg919; reg497=reg83*reg497;
    reg333=reg83*reg333; reg1288=reg83*reg1288; reg1097=reg83*reg1097; reg937=reg83*reg937; reg933=reg83*reg933;
    reg1118=reg83*reg1118; reg1277=reg83*reg1277; reg360=reg83*reg360; reg436=ponderation*reg242; reg1095=reg83*reg1095;
    reg811=reg83*reg811; reg1297=reg83*reg1297; reg251=reg83*reg251; reg488=reg83*reg488; reg1275=reg83*reg1275;
    reg439=reg83*reg439; reg1092=reg83*reg1092; reg756=reg83*reg756; reg478=reg83*reg478; reg1299=reg83*reg1299;
    reg1216=reg83*reg1216; reg941=reg83*reg941; reg1105=reg83*reg1105; reg978=reg83*reg978; reg1127=reg83*reg1127;
    reg448=ponderation*reg250; reg267=reg83*reg267; reg1125=reg83*reg1125; reg784=reg83*reg784; reg791=reg83*reg791;
    reg462=reg83*reg462; reg760=reg83*reg760; reg230=reg83*reg230; reg268=reg83*reg268; reg476=reg83*reg476;
    reg1217=reg83*reg1217; reg753=reg83*reg753; reg1131=reg83*reg1131; reg341=reg83*reg341; reg793=reg83*reg793;
    reg461=reg83*reg461; reg1214=reg83*reg1214; reg1193=reg83*reg1193; reg466=reg83*reg466; reg467=reg83*reg467;
    reg888=reg83*reg888; reg788=reg83*reg788; reg432=reg83*reg432; reg271=reg83*reg271; reg278=reg83*reg278;
    reg1103=reg83*reg1103; reg1124=reg83*reg1124; reg176=reg83*reg176; reg1188=reg83*reg1188; reg1198=reg83*reg1198;
    reg1152=reg83*reg1152; reg762=reg83*reg762; reg599=reg83*reg599; reg425=reg83*reg425; reg1197=reg83*reg1197;
    reg1150=reg83*reg1150; reg244=reg83*reg244; reg259=reg83*reg259; reg349=reg83*reg349; reg279=reg83*reg279;
    reg423=reg83*reg423; reg1149=reg83*reg1149; reg193=reg83*reg193; reg224=reg83*reg224; reg1146=reg83*reg1146;
    reg930=reg83*reg930; reg975=reg83*reg975; reg205=reg83*reg205; reg890=reg83*reg890; reg1196=reg83*reg1196;
    reg1295=reg83*reg1295; reg1121=reg83*reg1121; reg470=ponderation*reg394; reg1294=reg83*reg1294; reg1107=reg83*reg1107;
    reg758=reg83*reg758; reg1119=reg83*reg1119; reg789=reg83*reg789; reg940=reg83*reg940; reg1110=reg83*reg1110;
    reg364=reg83*reg364; reg471=ponderation*reg363; reg223=reg83*reg223; reg468=reg83*reg468; reg355=reg83*reg355;
    reg892=reg83*reg892; reg472=ponderation*reg178; reg1205=reg83*reg1205; reg415=reg83*reg415; reg264=reg83*reg264;
    reg748=reg83*reg748; reg519=reg83*reg519; reg589=reg83*reg589; reg524=reg83*reg524; reg993=reg83*reg993;
    reg549=reg83*reg549; reg745=reg83*reg745; reg517=reg83*reg517; reg1158=reg83*reg1158; reg587=reg83*reg587;
    reg1167=reg83*reg1167; reg1002=reg83*reg1002; reg848=reg83*reg848; reg561=reg83*reg561; reg387=reg83*reg387;
    reg742=reg83*reg742; reg301=reg83*reg301; reg475=ponderation*reg235; reg611=reg83*reg611; reg613=reg83*reg613;
    reg312=reg83*reg312; reg726=reg83*reg726; reg367=reg83*reg367; reg1215=reg83*reg1215; reg581=reg83*reg581;
    reg295=reg83*reg295; reg839=reg83*reg839; reg704=reg83*reg704; reg555=reg83*reg555; reg521=reg83*reg521;
    reg1186=reg83*reg1186; reg723=reg83*reg723; reg1208=reg83*reg1208; reg1165=reg83*reg1165; reg752=reg83*reg752;
    reg987=reg83*reg987; reg843=reg83*reg843; reg556=reg83*reg556; reg300=reg83*reg300; reg845=reg83*reg845;
    reg477=ponderation*reg232; reg826=reg83*reg826; reg225=reg83*reg225; reg567=reg83*reg567; reg1255=reg83*reg1255;
    reg359=reg83*reg359; reg1234=reg83*reg1234; reg695=reg83*reg695; reg260=reg83*reg260; reg481=ponderation*reg258;
    reg693=reg83*reg693; reg998=reg83*reg998; reg687=reg83*reg687; reg547=reg83*reg547; reg451=reg83*reg451;
    reg1260=reg83*reg1260; reg389=reg83*reg389; reg827=reg83*reg827; reg690=reg83*reg690; reg963=reg83*reg963;
    reg688=reg83*reg688; reg532=reg83*reg532; reg409=reg83*reg409; reg741=reg83*reg741; reg1031=reg83*reg1031;
    reg514=reg83*reg514; reg562=reg83*reg562; reg738=reg83*reg738; reg412=reg83*reg412; reg1000=reg83*reg1000;
    reg822=reg83*reg822; reg736=reg83*reg736; reg306=reg83*reg306; reg853=reg83*reg853; reg1178=reg83*reg1178;
    reg220=reg83*reg220; reg399=reg83*reg399; reg485=ponderation*reg362; reg701=reg83*reg701; reg486=ponderation*reg241;
    reg699=reg83*reg699; reg1256=reg83*reg1256; reg410=reg83*reg410; reg219=reg83*reg219; reg996=reg83*reg996;
    reg824=reg83*reg824; reg631=reg83*reg631; reg317=reg83*reg317; reg633=reg83*reg633; reg378=reg83*reg378;
    reg664=reg83*reg664; reg584=reg83*reg584; reg1264=reg83*reg1264; reg655=reg83*reg655; reg635=reg83*reg635;
    reg180=reg83*reg180; reg190=reg83*reg190; reg227=reg83*reg227; reg637=reg83*reg637; reg316=reg83*reg316;
    reg573=reg83*reg573; reg586=reg83*reg586; reg489=ponderation*reg179; reg382=reg83*reg382; reg629=reg83*reg629;
    reg984=reg83*reg984; reg1008=reg83*reg1008; reg652=reg83*reg652; reg380=reg83*reg380; reg1274=reg83*reg1274;
    reg659=reg83*reg659; reg662=reg83*reg662; reg656=reg83*reg656; reg577=reg83*reg577; reg1012=reg83*reg1012;
    reg579=reg83*reg579; reg204=reg83*reg204; reg1290=reg83*reg1290; reg315=reg83*reg315; reg680=reg83*reg680;
    reg668=reg83*reg668; reg173=reg83*reg173; reg676=reg83*reg676; reg1015=reg83*reg1015; reg597=reg83*reg597;
    reg670=reg83*reg670; reg582=reg83*reg582; reg575=reg83*reg575; reg673=reg83*reg673; reg201=reg83*reg201;
    reg1010=reg83*reg1010; reg550=reg83*reg550; reg648=reg83*reg648; reg734=reg83*reg734; reg654=reg83*reg654;
    reg552=reg83*reg552; reg988=reg83*reg988; reg620=reg83*reg620; reg639=reg83*reg639; reg553=reg83*reg553;
    reg1005=reg83*reg1005; reg834=reg83*reg834; reg731=reg83*reg731; reg384=reg83*reg384; reg417=reg83*reg417;
    reg729=reg83*reg729; reg990=reg83*reg990; reg523=reg83*reg523; reg578=reg83*reg578; reg583=reg83*reg583;
    reg837=reg83*reg837; reg554=reg83*reg554; reg605=reg83*reg605; reg683=reg83*reg683; reg588=reg83*reg588;
    reg985=reg83*reg985; reg211=reg83*reg211; reg679=reg83*reg679; reg647=reg83*reg647; reg388=reg83*reg388;
    reg571=reg83*reg571; reg644=reg83*reg644; reg208=reg83*reg208; reg177=reg83*reg177; reg685=reg83*reg685;
    reg643=reg83*reg643; reg1287=reg83*reg1287; reg314=reg83*reg314; reg498=ponderation*reg191; reg641=reg83*reg641;
    reg319=reg83*reg319; reg240=reg83*reg240; reg591=reg83*reg591; reg653=reg83*reg653; reg622=reg83*reg622;
    reg246=reg83*reg246; reg1173=reg83*reg1173; reg1076=reg83*reg1076; reg1261=reg83*reg1261; reg900=reg83*reg900;
    reg529=reg83*reg529; reg855=reg83*reg855; reg1065=reg83*reg1065; reg1073=reg83*reg1073; reg528=reg83*reg528;
    reg281=reg83*reg281; reg912=reg83*reg912; reg546=reg83*reg546; reg437=reg83*reg437; reg429=reg83*reg429;
    reg1080=reg83*reg1080; reg407=reg83*reg407; reg352=reg83*reg352; reg331=reg83*reg331; reg852=reg83*reg852;
    reg1239=reg83*reg1239; reg403=reg83*reg403; reg502=ponderation*reg337; reg1072=reg83*reg1072; reg1238=reg83*reg1238;
    reg1236=reg83*reg1236; reg953=reg83*reg953; reg395=reg83*reg395; reg218=reg83*reg218; reg1069=reg83*reg1069;
    reg1248=reg83*reg1248; reg491=reg83*reg491; reg507=ponderation*reg280; reg856=reg83*reg856; reg518=ponderation*reg248;
    reg1057=reg83*reg1057; reg870=reg83*reg870; reg903=reg83*reg903; reg867=reg83*reg867; reg275=reg83*reg275;
    reg490=reg83*reg490; reg434=reg83*reg434; reg307=reg83*reg307; reg539=reg83*reg539; reg234=reg83*reg234;
    reg370=reg83*reg370; reg435=reg83*reg435; reg1267=reg83*reg1267; reg1055=reg83*reg1055; reg393=reg83*reg393;
    reg908=reg83*reg908; reg872=reg83*reg872; reg1063=reg83*reg1063; reg329=reg83*reg329; reg614=reg83*reg614;
    reg520=ponderation*reg249; reg1244=reg83*reg1244; reg1272=reg83*reg1272; reg1052=reg83*reg1052; reg955=reg83*reg955;
    reg431=reg83*reg431; reg1263=reg83*reg1263; reg522=ponderation*reg298; reg1221=reg83*reg1221; reg1059=reg83*reg1059;
    reg527=ponderation*reg181; reg1228=reg83*reg1228; reg957=reg83*reg957; reg494=reg83*reg494; reg714=reg83*reg714;
    reg534=reg83*reg192; reg1241=reg83*reg1241; reg961=reg83*reg961; reg402=reg83*reg402; reg895=reg83*reg895;
    reg537=ponderation*reg88; reg195=reg83*reg195; reg1251=reg83*reg1251; reg548=ponderation*reg366; reg551=ponderation*reg308;
    reg716=reg83*reg716; reg557=ponderation*reg375; reg1226=reg83*reg1226; reg438=reg83*reg438; reg563=ponderation*reg304;
    reg440=reg83*reg440; reg869=reg83*reg869; reg717=reg83*reg717; reg361=reg83*reg361; reg357=reg83*reg357;
    reg1032=reg83*reg1032; reg951=reg83*reg951; reg719=reg83*reg719; reg1259=reg83*reg1259; reg339=reg83*reg339;
    reg1246=reg83*reg1246; reg293=reg83*reg293; reg1229=reg83*reg1229; reg453=reg83*reg453; reg564=ponderation*reg433;
    reg463=reg83*reg463; reg832=reg83*reg832; reg708=reg83*reg708; reg570=ponderation*reg285; reg391=reg83*reg391;
    reg574=ponderation*reg213; reg706=reg83*reg706; reg1085=reg83*reg1085; reg764=reg83*reg764; reg959=reg83*reg959;
    reg1082=reg83*reg1082; reg1237=reg83*reg1237; reg1289=reg83*reg1289; reg1270=reg83*reg1270; reg531=reg83*reg531;
    reg421=reg83*reg421; reg535=reg83*reg535; reg576=ponderation*reg174; reg538=reg83*reg538; reg580=ponderation*reg334;
    reg858=reg83*reg858; reg1169=reg83*reg1169; reg231=reg83*reg231; reg906=reg83*reg906; reg897=reg83*reg897;
    reg408=reg83*reg408; reg1243=reg83*reg1243; reg303=reg83*reg303; reg831=reg83*reg831; T tmp_21_13=ponderation*reg908;
    T tmp_2_16=ponderation*reg176; T tmp_21_8=ponderation*reg878; T tmp_23_9=-reg342; T tmp_1_4=ponderation*reg704; T tmp_2_14=ponderation*reg218;
    T tmp_20_15=ponderation*reg941; T tmp_3_2=ponderation*reg177; T tmp_3_20=ponderation*reg1294; T tmp_22_14=ponderation*reg1010; T tmp_23_22=ponderation*reg263;
    T tmp_21_9=ponderation*reg877; T tmp_23_6=ponderation*reg953; T tmp_3_12=ponderation*reg1270; T tmp_3_13=ponderation*reg1253; T tmp_22_11=ponderation*reg1015;
    T tmp_22_13=ponderation*reg315; T tmp_3_21=ponderation*reg335; T tmp_3_3=ponderation*reg1274; T tmp_21_11=-reg557; T tmp_0_2=ponderation*reg454;
    T tmp_21_10=ponderation*reg874; T tmp_20_8=ponderation*reg445; T tmp_22_12=ponderation*reg1012; T tmp_23_7=ponderation*reg951; T tmp_23_23=ponderation*reg354;
    T tmp_22_8=ponderation*reg985; T tmp_3_4=ponderation*reg180; T tmp_2_10=ponderation*reg1283; T tmp_0_0=ponderation*reg239; T tmp_22_15=ponderation*reg1008;
    T tmp_20_14=ponderation*reg568; T tmp_23_21=ponderation*reg966; T tmp_23_5=ponderation*reg275; T tmp_22_9=ponderation*reg984; T tmp_1_5=ponderation*reg316;
    T tmp_2_20=ponderation*reg358; T tmp_2_21=ponderation*reg1169; T tmp_21_12=ponderation*reg546; T tmp_20_13=ponderation*reg943; T tmp_2_13=ponderation*reg1264;
    T tmp_22_10=ponderation*reg317; T tmp_23_8=ponderation*reg386; T tmp_0_1=ponderation*reg889; T tmp_22_21=ponderation*reg1000; T tmp_20_23=ponderation*reg565;
    T tmp_21_3=ponderation*reg888; T tmp_21_18=ponderation*reg535; T tmp_22_0=ponderation*reg306; T tmp_1_6=ponderation*reg987; T tmp_23_1=ponderation*reg959;
    T tmp_3_0=ponderation*reg399; T tmp_20_22=ponderation*reg931; T tmp_3_7=ponderation*reg1167; T tmp_23_12=ponderation*reg980; T tmp_23_16=ponderation*reg974;
    T tmp_22_20=ponderation*reg1002; T tmp_22_1=ponderation*reg301; T tmp_21_4=ponderation*reg886; T tmp_20_21=ponderation*reg933; T tmp_3_10=ponderation*reg1238;
    T tmp_2_19=ponderation*reg1193; T tmp_3_15=ponderation*reg1201; T tmp_22_2=ponderation*reg993; T tmp_21_1=ponderation*reg892; T tmp_21_21=ponderation*reg532;
    T tmp_22_23=ponderation*reg963; T tmp_2_12=ponderation*reg1263; T tmp_2_23=ponderation*reg361; T tmp_23_14=ponderation*reg364; T tmp_21_2=ponderation*reg890;
    T tmp_2_18=ponderation*reg188; T tmp_22_22=ponderation*reg260; T tmp_3_17=ponderation*reg1197; T tmp_21_20=ponderation*reg895; T tmp_23_13=ponderation*reg978;
    T tmp_21_22=ponderation*reg998; T tmp_21_0=ponderation*reg930; T tmp_3_8=ponderation*reg1255; T tmp_3_16=ponderation*reg1216; T tmp_3_9=ponderation*reg1226;
    T tmp_23_0=ponderation*reg961; T tmp_23_15=ponderation*reg975; T tmp_21_23=ponderation*reg996; T tmp_21_19=ponderation*reg897; T tmp_21_6=ponderation*reg882;
    T tmp_23_10=ponderation*reg947; T tmp_22_17=ponderation*reg1005; T tmp_23_19=ponderation*reg968; T tmp_22_5=ponderation*reg990; T tmp_2_11=ponderation*reg592;
    T tmp_20_17=ponderation*reg175; T tmp_3_11=-reg520; T tmp_2_15=ponderation*reg1209; T tmp_22_6=ponderation*reg988; T tmp_23_4=ponderation*reg955;
    T tmp_3_5=ponderation*reg622; T tmp_21_7=ponderation*reg880; T tmp_1_7=ponderation*reg912; T tmp_22_16=ponderation*reg314; T tmp_1_3=ponderation*reg950;
    T tmp_3_14=ponderation*reg1177; T tmp_20_16=ponderation*reg940; T tmp_21_14=ponderation*reg906; T tmp_22_7=ponderation*reg319; T tmp_23_20=ponderation*reg355;
    T tmp_3_18=ponderation*reg333; T tmp_23_17=ponderation*reg360; T tmp_20_20=ponderation*reg566; T tmp_22_19=ponderation*reg264; T tmp_21_17=ponderation*reg900;
    T tmp_21_5=ponderation*reg884; T tmp_23_2=ponderation*reg281; T tmp_22_3=ponderation*reg300; T tmp_2_17=ponderation*reg356; T tmp_23_11=ponderation*reg512;
    T tmp_2_22=ponderation*reg1228; T tmp_20_19=ponderation*reg935; T tmp_3_6=ponderation*reg367; T tmp_23_18=ponderation*reg970; T tmp_21_16=ponderation*reg903;
    T tmp_22_18=ponderation*reg312; T tmp_22_4=ponderation*reg295; T tmp_23_3=ponderation*reg957; T tmp_3_19=ponderation*reg340; T tmp_20_18=ponderation*reg937;
    T tmp_21_15=ponderation*reg539; T tmp_3_1=ponderation*reg1165; T tmp_11_23=ponderation*reg517; T tmp_2_1=ponderation*reg1234; T tmp_12_0=ponderation*reg519;
    T tmp_7_0=ponderation*reg1205; T tmp_12_1=ponderation*reg845; T tmp_12_2=ponderation*reg843; T tmp_6_23=ponderation*reg1186; T tmp_12_3=ponderation*reg521;
    T tmp_12_4=ponderation*reg839; T tmp_6_22=ponderation*reg613; T tmp_12_5=ponderation*reg837; T tmp_1_15=ponderation*reg654; T tmp_12_6=ponderation*reg523;
    T tmp_6_21=ponderation*reg384; T tmp_12_7=ponderation*reg834; T tmp_12_8=ponderation*reg639; T tmp_12_9=ponderation*reg648; T tmp_6_20=ponderation*reg240;
    T tmp_12_10=ponderation*reg653; T tmp_12_11=-reg498; T tmp_1_14=ponderation*reg685; T tmp_6_19=ponderation*reg388; T tmp_12_12=ponderation*reg571;
    T tmp_12_13=ponderation*reg679; T tmp_12_14=ponderation*reg683; T tmp_6_18=ponderation*reg382; T tmp_12_15=ponderation*reg573; T tmp_12_16=ponderation*reg655;
    T tmp_7_7=ponderation*reg391; T tmp_11_6=-reg570; T tmp_11_7=ponderation*reg832; T tmp_7_6=ponderation*reg1243; T tmp_11_8=ponderation*reg538;
    T tmp_11_9=ponderation*reg534; T tmp_0_15=ponderation*reg195; T tmp_11_10=-reg551; T tmp_7_5=ponderation*reg1229; T tmp_11_11=ponderation*reg339;
    T tmp_11_12=-reg564; T tmp_7_4=ponderation*reg389; T tmp_11_13=ponderation*reg827; T tmp_11_14=ponderation*reg547; T tmp_0_14=ponderation*reg824;
    T tmp_11_15=-reg481; T tmp_7_3=ponderation*reg1256; T tmp_11_16=ponderation*reg826; T tmp_11_17=-reg486; T tmp_11_18=-reg485;
    T tmp_7_2=ponderation*reg1178; T tmp_11_19=ponderation*reg822; T tmp_0_13=ponderation*reg1031; T tmp_11_20=ponderation*reg514; T tmp_11_21=-reg475;
    T tmp_7_1=ponderation*reg387; T tmp_11_22=ponderation*reg848; T tmp_1_16=ponderation*reg549; T tmp_13_10=ponderation*reg586; T tmp_13_11=ponderation*reg629;
    T tmp_2_3=ponderation*reg190; T tmp_13_12=ponderation*reg652; T tmp_1_1=ponderation*reg208; T tmp_6_10=ponderation*reg211; T tmp_13_13=ponderation*reg588;
    T tmp_13_14=ponderation*reg647; T tmp_13_15=ponderation*reg644; T tmp_6_9=ponderation*reg1287; T tmp_13_16=ponderation*reg643; T tmp_13_17=ponderation*reg641;
    T tmp_13_18=ponderation*reg591; T tmp_6_8=ponderation*reg246; T tmp_13_19=ponderation*reg550; T tmp_13_20=ponderation*reg734; T tmp_6_7=ponderation*reg620;
    T tmp_13_21=ponderation*reg552; T tmp_1_0=ponderation*reg583; T tmp_13_22=ponderation*reg553; T tmp_13_23=ponderation*reg731; T tmp_6_6=ponderation*reg417;
    T tmp_14_0=ponderation*reg729; T tmp_14_1=ponderation*reg578; T tmp_14_2=ponderation*reg554; T tmp_14_3=ponderation*reg726; T tmp_6_5=ponderation*reg1215;
    T tmp_14_4=ponderation*reg581; T tmp_2_2=ponderation*reg227; T tmp_12_17=ponderation*reg664; T tmp_6_17=ponderation*reg201; T tmp_1_13=ponderation*reg204;
    T tmp_12_18=ponderation*reg575; T tmp_12_19=ponderation*reg670; T tmp_6_16=ponderation*reg173; T tmp_12_20=ponderation*reg668; T tmp_12_21=ponderation*reg577;
    T tmp_6_15=ponderation*reg380; T tmp_12_22=ponderation*reg662; T tmp_12_23=ponderation*reg659; T tmp_13_0=ponderation*reg656; T tmp_6_14=ponderation*reg1290;
    T tmp_13_1=ponderation*reg579; T tmp_13_2=ponderation*reg680; T tmp_13_3=ponderation*reg676; T tmp_6_13=ponderation*reg597; T tmp_1_12=ponderation*reg784;
    T tmp_13_4=ponderation*reg582; T tmp_13_5=ponderation*reg673; T tmp_6_12=ponderation*reg378; T tmp_13_6=ponderation*reg631; T tmp_1_2=ponderation*reg633;
    T tmp_13_7=ponderation*reg584; T tmp_13_8=ponderation*reg635; T tmp_13_9=ponderation*reg637; T tmp_6_11=-reg489; T tmp_8_4=ponderation*reg1278;
    T tmp_9_5=-reg416; T tmp_9_6=ponderation*reg251; T tmp_9_7=ponderation*reg811; T tmp_8_3=ponderation*reg1288; T tmp_20_0=ponderation*reg919;
    T tmp_20_1=ponderation*reg808; T tmp_8_2=ponderation*reg269; T tmp_20_2=ponderation*reg479; T tmp_8_1=ponderation*reg593; T tmp_9_9=ponderation*reg805;
    T tmp_1_18=ponderation*reg765; T tmp_9_10=-reg413; T tmp_9_11=ponderation*reg212; T tmp_8_0=ponderation*reg193; T tmp_9_12=ponderation*reg259;
    T tmp_9_13=ponderation*reg762; T tmp_7_23=ponderation*reg1198; T tmp_9_14=-reg471; T tmp_1_17=ponderation*reg758; T tmp_9_15=-reg470;
    T tmp_7_22=ponderation*reg268; T tmp_9_16=ponderation*reg760; T tmp_20_9=-reg344; T tmp_20_10=ponderation*reg756; T tmp_7_21=ponderation*reg1214;
    T tmp_20_11=ponderation*reg467; T tmp_20_12=ponderation*reg753; T tmp_8_11=ponderation*reg418; T tmp_8_12=ponderation*reg803; T tmp_8_10=ponderation*reg221;
    T tmp_8_13=ponderation*reg801; T tmp_8_14=ponderation*reg276; T tmp_8_15=ponderation*reg799; T tmp_8_9=-reg377; T tmp_8_16=ponderation*reg798;
    T tmp_1_21=ponderation*reg788; T tmp_8_17=ponderation*reg277; T tmp_8_18=ponderation*reg795; T tmp_8_8=ponderation*reg271; T tmp_8_19=ponderation*reg793;
    T tmp_8_20=ponderation*reg278; T tmp_8_21=ponderation*reg791; T tmp_8_7=ponderation*reg599; T tmp_8_22=ponderation*reg789; T tmp_8_23=ponderation*reg279;
    T tmp_8_6=ponderation*reg390; T tmp_9_0=ponderation*reg821; T tmp_9_1=ponderation*reg819; T tmp_1_20=ponderation*reg818; T tmp_9_2=-reg436;
    T tmp_8_5=ponderation*reg270; T tmp_9_3=ponderation*reg817; T tmp_9_4=ponderation*reg815; T tmp_1_22=ponderation*reg230; T tmp_1_19=ponderation*reg233;
    T tmp_10_11=-reg563; T tmp_10_12=ponderation*reg858; T tmp_10_13=ponderation*reg764; T tmp_7_13=ponderation*reg395; T tmp_10_14=ponderation*reg856;
    T tmp_10_15=ponderation*reg855; T tmp_0_10=ponderation*reg869; T tmp_7_12=ponderation*reg1267; T tmp_10_16=-reg518; T tmp_10_17=ponderation*reg852;
    T tmp_10_18=ponderation*reg303; T tmp_7_11=ponderation*reg352; T tmp_10_19=ponderation*reg1289; T tmp_10_20=ponderation*reg872; T tmp_10_21=ponderation*reg870;
    T tmp_7_10=ponderation*reg393; T tmp_10_22=ponderation*reg307; T tmp_2_0=ponderation*reg1173; T tmp_10_23=ponderation*reg867; T tmp_0_9=ponderation*reg529;
    T tmp_7_9=ponderation*reg1244; T tmp_11_0=-reg522; T tmp_11_1=ponderation*reg329; T tmp_11_2=ponderation*reg528; T tmp_7_8=ponderation*reg1236;
    T tmp_11_3=-reg507; T tmp_11_4=ponderation*reg331; T tmp_11_5=ponderation*reg531; T tmp_9_19=ponderation*reg785; T tmp_7_20=ponderation*reg1206;
    T tmp_1_11=ponderation*reg779; T tmp_9_20=-reg286; T tmp_9_21=ponderation*reg783; T tmp_7_19=ponderation*reg266; T tmp_9_22=ponderation*reg781;
    T tmp_9_23=-reg247; T tmp_1_23=ponderation*reg1190; T tmp_10_0=ponderation*reg777; T tmp_1_10=ponderation*reg198; T tmp_7_18=ponderation*reg1162;
    T tmp_10_1=ponderation*reg328; T tmp_10_2=ponderation*reg775; T tmp_10_3=ponderation*reg773; T tmp_7_17=ponderation*reg1183; T tmp_10_4=ponderation*reg330;
    T tmp_10_5=ponderation*reg771; T tmp_7_16=ponderation*reg396; T tmp_10_6=ponderation*reg288; T tmp_10_7=ponderation*reg813; T tmp_10_8=ponderation*reg863;
    T tmp_1_9=ponderation*reg860; T tmp_7_15=ponderation*reg203; T tmp_10_9=-reg325; T tmp_10_10=ponderation*reg293; T tmp_1_8=ponderation*reg853;
    T tmp_7_14=ponderation*reg1251; T tmp_17_23=ponderation*reg432; T tmp_18_0=ponderation*reg478; T tmp_18_1=ponderation*reg1127; T tmp_4_16=ponderation*reg267;
    T tmp_18_2=ponderation*reg1125; T tmp_18_3=ponderation*reg476; T tmp_4_15=ponderation*reg1196; T tmp_18_4=ponderation*reg1121; T tmp_18_5=ponderation*reg1119;
    T tmp_0_21=ponderation*reg244; T tmp_4_14=ponderation*reg223; T tmp_18_6=ponderation*reg468; T tmp_18_7=ponderation*reg1152; T tmp_4_13=ponderation*reg349;
    T tmp_18_8=ponderation*reg1150; T tmp_18_9=ponderation*reg1149; T tmp_18_10=ponderation*reg1146; T tmp_2_8=ponderation*reg224; T tmp_4_12=ponderation*reg205;
    T tmp_18_11=-reg411; T tmp_4_11=ponderation*reg151; T tmp_18_12=ponderation*reg501; T tmp_18_13=ponderation*reg1141; T tmp_18_14=ponderation*reg1140;
    T tmp_0_20=ponderation*reg1100; T tmp_4_10=ponderation*reg596; T tmp_18_15=ponderation*reg499; T tmp_18_16=ponderation*reg1138; T tmp_17_5=ponderation*reg452;
    T tmp_17_6=ponderation*reg1050; T tmp_4_23=ponderation*reg1163; T tmp_17_7=ponderation*reg1048; T tmp_17_8=ponderation*reg455; T tmp_4_22=ponderation*reg400;
    T tmp_17_9=-reg383; T tmp_17_10=ponderation*reg1046; T tmp_0_11=-reg282; T tmp_17_11=-reg365; T tmp_4_21=ponderation*reg265;
    T tmp_17_12=ponderation*reg1042; T tmp_17_13=ponderation*reg1040; T tmp_17_14=ponderation*reg464; T tmp_4_20=ponderation*reg1207; T tmp_17_15=ponderation*reg1038;
    T tmp_17_16=-reg299; T tmp_17_17=ponderation*reg484; T tmp_4_19=ponderation*reg398; T tmp_17_18=ponderation*reg1135; T tmp_17_19=ponderation*reg1133;
    T tmp_2_7=ponderation*reg1176; T tmp_4_18=ponderation*reg1217; T tmp_17_20=ponderation*reg430; T tmp_17_21=ponderation*reg1131; T tmp_17_22=ponderation*reg466;
    T tmp_0_22=ponderation*reg1124; T tmp_4_17=ponderation*reg1188; T tmp_19_10=ponderation*reg425; T tmp_19_11=ponderation*reg1110; T tmp_2_9=-reg472;
    T tmp_4_3=ponderation*reg1295; T tmp_19_12=ponderation*reg1107; T tmp_19_13=ponderation*reg462; T tmp_19_14=ponderation*reg1105; T tmp_4_2=ponderation*reg1299;
    T tmp_19_15=ponderation*reg1103; T tmp_19_16=ponderation*reg461; T tmp_4_1=ponderation*reg341; T tmp_19_17=ponderation*reg928; T tmp_19_18=ponderation*reg926;
    T tmp_0_6=ponderation*reg245; T tmp_19_19=ponderation*reg449; T tmp_4_0=ponderation*reg196; T tmp_19_20=ponderation*reg924; T tmp_19_21=ponderation*reg922;
    T tmp_9_8=-reg381; T tmp_19_22=ponderation*reg447; T tmp_19_23=ponderation*reg920; T tmp_20_3=ponderation*reg806; T tmp_3_23=ponderation*reg608;
    T tmp_20_4=ponderation*reg917; T tmp_20_5=ponderation*reg446; T tmp_20_6=ponderation*reg915; T tmp_3_22=ponderation*reg187; T tmp_20_7=ponderation*reg913;
    T tmp_9_18=ponderation*reg347; T tmp_18_17=ponderation*reg1101; T tmp_18_18=ponderation*reg497; T tmp_9_17=-reg448; T tmp_4_9=ponderation*reg1277;
    T tmp_18_19=ponderation*reg1097; T tmp_18_20=ponderation*reg1095; T tmp_0_19=ponderation*reg1068; T tmp_4_8=ponderation*reg1275; T tmp_18_21=ponderation*reg488;
    T tmp_18_22=ponderation*reg1092; T tmp_18_23=ponderation*reg1090; T tmp_4_7=ponderation*reg345; T tmp_19_0=ponderation*reg487; T tmp_0_8=ponderation*reg1115;
    T tmp_19_1=ponderation*reg442; T tmp_19_2=ponderation*reg1087; T tmp_19_3=ponderation*reg441; T tmp_4_6=ponderation*reg1297; T tmp_19_4=ponderation*reg439;
    T tmp_19_5=ponderation*reg1118; T tmp_19_6=ponderation*reg1116; T tmp_4_5=ponderation*reg1293; T tmp_19_7=ponderation*reg428; T tmp_19_8=ponderation*reg1113;
    T tmp_4_4=ponderation*reg343; T tmp_19_9=ponderation*reg1112; T tmp_0_7=ponderation*reg423; T tmp_15_0=ponderation*reg567; T tmp_5_22=ponderation*reg359;
    T tmp_15_1=ponderation*reg695; T tmp_15_2=ponderation*reg693; T tmp_5_21=ponderation*reg1260; T tmp_15_3=ponderation*reg451; T tmp_15_4=ponderation*reg690;
    T tmp_15_5=ponderation*reg688; T tmp_5_20=ponderation*reg409; T tmp_15_6=ponderation*reg453; T tmp_2_5=ponderation*reg225; T tmp_5_19=ponderation*reg357;
    T tmp_15_7=ponderation*reg719; T tmp_15_8=ponderation*reg717; T tmp_0_4=ponderation*reg716; T tmp_15_9=-reg537; T tmp_5_18=ponderation*reg1241;
    T tmp_15_10=ponderation*reg714; T tmp_15_11=-reg580; T tmp_0_3=ponderation*reg231; T tmp_5_17=ponderation*reg408; T tmp_15_12=ponderation*reg463;
    T tmp_15_13=ponderation*reg708; T tmp_5_16=ponderation*reg1237; T tmp_15_14=ponderation*reg706; T tmp_15_15=ponderation*reg421; T tmp_5_15=ponderation*reg1239;
    T tmp_15_16=-reg502; T tmp_14_5=ponderation*reg555; T tmp_14_6=ponderation*reg723; T tmp_6_4=ponderation*reg1208; T tmp_14_7=ponderation*reg752;
    T tmp_14_8=ponderation*reg556; T tmp_0_23=ponderation*reg589; T tmp_6_3=ponderation*reg415; T tmp_14_9=-reg477; T tmp_14_10=ponderation*reg748;
    T tmp_14_11=ponderation*reg524; T tmp_2_4=ponderation*reg605; T tmp_14_12=ponderation*reg745; T tmp_6_2=ponderation*reg1158; T tmp_14_13=ponderation*reg587;
    T tmp_14_14=ponderation*reg561; T tmp_14_15=ponderation*reg742; T tmp_6_1=ponderation*reg611; T tmp_14_16=ponderation*reg741; T tmp_14_17=ponderation*reg562;
    T tmp_6_0=ponderation*reg412; T tmp_14_18=ponderation*reg738; T tmp_14_19=ponderation*reg736; T tmp_14_20=ponderation*reg220; T tmp_14_21=ponderation*reg701;
    T tmp_5_23=ponderation*reg410; T tmp_14_22=ponderation*reg699; T tmp_14_23=ponderation*reg219; T tmp_0_5=ponderation*reg687; T tmp_5_7=ponderation*reg1272;
    T tmp_16_11=ponderation*reg1076; T tmp_16_12=ponderation*reg1073; T tmp_5_6=ponderation*reg1248; T tmp_16_13=ponderation*reg491; T tmp_16_14=ponderation*reg1072;
    T tmp_0_16=ponderation*reg831; T tmp_16_15=-reg574; T tmp_5_5=ponderation*reg402; T tmp_16_16=ponderation*reg494; T tmp_16_17=-reg548;
    T tmp_5_4=ponderation*reg1259; T tmp_16_18=ponderation*reg1032; T tmp_16_19=ponderation*reg465; T tmp_5_3=ponderation*reg1168; T tmp_16_20=ponderation*reg1030;
    T tmp_16_21=ponderation*reg1027; T tmp_0_12=ponderation*reg214; T tmp_16_22=ponderation*reg469; T tmp_5_2=ponderation*reg401; T tmp_16_23=ponderation*reg1025;
    T tmp_17_0=ponderation*reg1022; T tmp_5_1=ponderation*reg1182; T tmp_17_1=ponderation*reg1020; T tmp_17_2=ponderation*reg450; T tmp_17_3=ponderation*reg1018;
    T tmp_5_0=ponderation*reg1159; T tmp_17_4=ponderation*reg1016; T tmp_15_17=ponderation*reg1069; T tmp_5_14=ponderation*reg407; T tmp_15_18=ponderation*reg429;
    T tmp_15_19=ponderation*reg1065; T tmp_5_13=ponderation*reg614; T tmp_15_20=ponderation*reg1063; T tmp_15_21=ponderation*reg431; T tmp_5_12=ponderation*reg1221;
    T tmp_15_22=ponderation*reg1059; T tmp_15_23=ponderation*reg1057; T tmp_5_11=ponderation*reg370; T tmp_16_0=ponderation*reg434; T tmp_0_18=ponderation*reg234;
    T tmp_16_1=ponderation*reg435; T tmp_16_2=ponderation*reg1055; T tmp_5_10=ponderation*reg1261; T tmp_16_3=ponderation*reg1052; T tmp_16_4=ponderation*reg437;
    T tmp_16_5=ponderation*reg1085; T tmp_5_9=-reg576; T tmp_16_6=ponderation*reg438; T tmp_16_7=ponderation*reg440; T tmp_2_6=ponderation*reg1246;
    T tmp_5_8=ponderation*reg403; T tmp_16_8=ponderation*reg1082; T tmp_16_9=ponderation*reg1080; T tmp_0_17=ponderation*reg490; T tmp_16_10=-reg527;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[0]+2) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_3;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_4;
    matrix(indices[0]+0,indices[1]+2) += tmp_0_5;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_7;
    matrix(indices[0]+0,indices[2]+2) += tmp_0_8;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_9;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_10;
    matrix(indices[0]+0,indices[3]+2) += tmp_0_11;
    matrix(indices[0]+0,indices[4]+0) += tmp_0_12;
    matrix(indices[0]+0,indices[4]+1) += tmp_0_13;
    matrix(indices[0]+0,indices[4]+2) += tmp_0_14;
    matrix(indices[0]+0,indices[5]+0) += tmp_0_15;
    matrix(indices[0]+0,indices[5]+1) += tmp_0_16;
    matrix(indices[0]+0,indices[5]+2) += tmp_0_17;
    matrix(indices[0]+0,indices[6]+0) += tmp_0_18;
    matrix(indices[0]+0,indices[6]+1) += tmp_0_19;
    matrix(indices[0]+0,indices[6]+2) += tmp_0_20;
    matrix(indices[0]+0,indices[7]+0) += tmp_0_21;
    matrix(indices[0]+0,indices[7]+1) += tmp_0_22;
    matrix(indices[0]+0,indices[7]+2) += tmp_0_23;
    matrix(indices[0]+1,indices[0]+0) += tmp_1_0;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[0]+2) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_3;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_4;
    matrix(indices[0]+1,indices[1]+2) += tmp_1_5;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_7;
    matrix(indices[0]+1,indices[2]+2) += tmp_1_8;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_9;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_10;
    matrix(indices[0]+1,indices[3]+2) += tmp_1_11;
    matrix(indices[0]+1,indices[4]+0) += tmp_1_12;
    matrix(indices[0]+1,indices[4]+1) += tmp_1_13;
    matrix(indices[0]+1,indices[4]+2) += tmp_1_14;
    matrix(indices[0]+1,indices[5]+0) += tmp_1_15;
    matrix(indices[0]+1,indices[5]+1) += tmp_1_16;
    matrix(indices[0]+1,indices[5]+2) += tmp_1_17;
    matrix(indices[0]+1,indices[6]+0) += tmp_1_18;
    matrix(indices[0]+1,indices[6]+1) += tmp_1_19;
    matrix(indices[0]+1,indices[6]+2) += tmp_1_20;
    matrix(indices[0]+1,indices[7]+0) += tmp_1_21;
    matrix(indices[0]+1,indices[7]+1) += tmp_1_22;
    matrix(indices[0]+1,indices[7]+2) += tmp_1_23;
    matrix(indices[0]+2,indices[0]+0) += tmp_2_0;
    matrix(indices[0]+2,indices[0]+1) += tmp_2_1;
    matrix(indices[0]+2,indices[0]+2) += tmp_2_2;
    matrix(indices[0]+2,indices[1]+0) += tmp_2_3;
    matrix(indices[0]+2,indices[1]+1) += tmp_2_4;
    matrix(indices[0]+2,indices[1]+2) += tmp_2_5;
    matrix(indices[0]+2,indices[2]+0) += tmp_2_6;
    matrix(indices[0]+2,indices[2]+1) += tmp_2_7;
    matrix(indices[0]+2,indices[2]+2) += tmp_2_8;
    matrix(indices[0]+2,indices[3]+0) += tmp_2_9;
    matrix(indices[0]+2,indices[3]+1) += tmp_2_10;
    matrix(indices[0]+2,indices[3]+2) += tmp_2_11;
    matrix(indices[0]+2,indices[4]+0) += tmp_2_12;
    matrix(indices[0]+2,indices[4]+1) += tmp_2_13;
    matrix(indices[0]+2,indices[4]+2) += tmp_2_14;
    matrix(indices[0]+2,indices[5]+0) += tmp_2_15;
    matrix(indices[0]+2,indices[5]+1) += tmp_2_16;
    matrix(indices[0]+2,indices[5]+2) += tmp_2_17;
    matrix(indices[0]+2,indices[6]+0) += tmp_2_18;
    matrix(indices[0]+2,indices[6]+1) += tmp_2_19;
    matrix(indices[0]+2,indices[6]+2) += tmp_2_20;
    matrix(indices[0]+2,indices[7]+0) += tmp_2_21;
    matrix(indices[0]+2,indices[7]+1) += tmp_2_22;
    matrix(indices[0]+2,indices[7]+2) += tmp_2_23;
    matrix(indices[1]+0,indices[0]+0) += tmp_3_0;
    matrix(indices[1]+0,indices[0]+1) += tmp_3_1;
    matrix(indices[1]+0,indices[0]+2) += tmp_3_2;
    matrix(indices[1]+0,indices[1]+0) += tmp_3_3;
    matrix(indices[1]+0,indices[1]+1) += tmp_3_4;
    matrix(indices[1]+0,indices[1]+2) += tmp_3_5;
    matrix(indices[1]+0,indices[2]+0) += tmp_3_6;
    matrix(indices[1]+0,indices[2]+1) += tmp_3_7;
    matrix(indices[1]+0,indices[2]+2) += tmp_3_8;
    matrix(indices[1]+0,indices[3]+0) += tmp_3_9;
    matrix(indices[1]+0,indices[3]+1) += tmp_3_10;
    matrix(indices[1]+0,indices[3]+2) += tmp_3_11;
    matrix(indices[1]+0,indices[4]+0) += tmp_3_12;
    matrix(indices[1]+0,indices[4]+1) += tmp_3_13;
    matrix(indices[1]+0,indices[4]+2) += tmp_3_14;
    matrix(indices[1]+0,indices[5]+0) += tmp_3_15;
    matrix(indices[1]+0,indices[5]+1) += tmp_3_16;
    matrix(indices[1]+0,indices[5]+2) += tmp_3_17;
    matrix(indices[1]+0,indices[6]+0) += tmp_3_18;
    matrix(indices[1]+0,indices[6]+1) += tmp_3_19;
    matrix(indices[1]+0,indices[6]+2) += tmp_3_20;
    matrix(indices[1]+0,indices[7]+0) += tmp_3_21;
    matrix(indices[1]+0,indices[7]+1) += tmp_3_22;
    matrix(indices[1]+0,indices[7]+2) += tmp_3_23;
    matrix(indices[1]+1,indices[0]+0) += tmp_4_0;
    matrix(indices[1]+1,indices[0]+1) += tmp_4_1;
    matrix(indices[1]+1,indices[0]+2) += tmp_4_2;
    matrix(indices[1]+1,indices[1]+0) += tmp_4_3;
    matrix(indices[1]+1,indices[1]+1) += tmp_4_4;
    matrix(indices[1]+1,indices[1]+2) += tmp_4_5;
    matrix(indices[1]+1,indices[2]+0) += tmp_4_6;
    matrix(indices[1]+1,indices[2]+1) += tmp_4_7;
    matrix(indices[1]+1,indices[2]+2) += tmp_4_8;
    matrix(indices[1]+1,indices[3]+0) += tmp_4_9;
    matrix(indices[1]+1,indices[3]+1) += tmp_4_10;
    matrix(indices[1]+1,indices[3]+2) += tmp_4_11;
    matrix(indices[1]+1,indices[4]+0) += tmp_4_12;
    matrix(indices[1]+1,indices[4]+1) += tmp_4_13;
    matrix(indices[1]+1,indices[4]+2) += tmp_4_14;
    matrix(indices[1]+1,indices[5]+0) += tmp_4_15;
    matrix(indices[1]+1,indices[5]+1) += tmp_4_16;
    matrix(indices[1]+1,indices[5]+2) += tmp_4_17;
    matrix(indices[1]+1,indices[6]+0) += tmp_4_18;
    matrix(indices[1]+1,indices[6]+1) += tmp_4_19;
    matrix(indices[1]+1,indices[6]+2) += tmp_4_20;
    matrix(indices[1]+1,indices[7]+0) += tmp_4_21;
    matrix(indices[1]+1,indices[7]+1) += tmp_4_22;
    matrix(indices[1]+1,indices[7]+2) += tmp_4_23;
    matrix(indices[1]+2,indices[0]+0) += tmp_5_0;
    matrix(indices[1]+2,indices[0]+1) += tmp_5_1;
    matrix(indices[1]+2,indices[0]+2) += tmp_5_2;
    matrix(indices[1]+2,indices[1]+0) += tmp_5_3;
    matrix(indices[1]+2,indices[1]+1) += tmp_5_4;
    matrix(indices[1]+2,indices[1]+2) += tmp_5_5;
    matrix(indices[1]+2,indices[2]+0) += tmp_5_6;
    matrix(indices[1]+2,indices[2]+1) += tmp_5_7;
    matrix(indices[1]+2,indices[2]+2) += tmp_5_8;
    matrix(indices[1]+2,indices[3]+0) += tmp_5_9;
    matrix(indices[1]+2,indices[3]+1) += tmp_5_10;
    matrix(indices[1]+2,indices[3]+2) += tmp_5_11;
    matrix(indices[1]+2,indices[4]+0) += tmp_5_12;
    matrix(indices[1]+2,indices[4]+1) += tmp_5_13;
    matrix(indices[1]+2,indices[4]+2) += tmp_5_14;
    matrix(indices[1]+2,indices[5]+0) += tmp_5_15;
    matrix(indices[1]+2,indices[5]+1) += tmp_5_16;
    matrix(indices[1]+2,indices[5]+2) += tmp_5_17;
    matrix(indices[1]+2,indices[6]+0) += tmp_5_18;
    matrix(indices[1]+2,indices[6]+1) += tmp_5_19;
    matrix(indices[1]+2,indices[6]+2) += tmp_5_20;
    matrix(indices[1]+2,indices[7]+0) += tmp_5_21;
    matrix(indices[1]+2,indices[7]+1) += tmp_5_22;
    matrix(indices[1]+2,indices[7]+2) += tmp_5_23;
    matrix(indices[2]+0,indices[0]+0) += tmp_6_0;
    matrix(indices[2]+0,indices[0]+1) += tmp_6_1;
    matrix(indices[2]+0,indices[0]+2) += tmp_6_2;
    matrix(indices[2]+0,indices[1]+0) += tmp_6_3;
    matrix(indices[2]+0,indices[1]+1) += tmp_6_4;
    matrix(indices[2]+0,indices[1]+2) += tmp_6_5;
    matrix(indices[2]+0,indices[2]+0) += tmp_6_6;
    matrix(indices[2]+0,indices[2]+1) += tmp_6_7;
    matrix(indices[2]+0,indices[2]+2) += tmp_6_8;
    matrix(indices[2]+0,indices[3]+0) += tmp_6_9;
    matrix(indices[2]+0,indices[3]+1) += tmp_6_10;
    matrix(indices[2]+0,indices[3]+2) += tmp_6_11;
    matrix(indices[2]+0,indices[4]+0) += tmp_6_12;
    matrix(indices[2]+0,indices[4]+1) += tmp_6_13;
    matrix(indices[2]+0,indices[4]+2) += tmp_6_14;
    matrix(indices[2]+0,indices[5]+0) += tmp_6_15;
    matrix(indices[2]+0,indices[5]+1) += tmp_6_16;
    matrix(indices[2]+0,indices[5]+2) += tmp_6_17;
    matrix(indices[2]+0,indices[6]+0) += tmp_6_18;
    matrix(indices[2]+0,indices[6]+1) += tmp_6_19;
    matrix(indices[2]+0,indices[6]+2) += tmp_6_20;
    matrix(indices[2]+0,indices[7]+0) += tmp_6_21;
    matrix(indices[2]+0,indices[7]+1) += tmp_6_22;
    matrix(indices[2]+0,indices[7]+2) += tmp_6_23;
    matrix(indices[2]+1,indices[0]+0) += tmp_7_0;
    matrix(indices[2]+1,indices[0]+1) += tmp_7_1;
    matrix(indices[2]+1,indices[0]+2) += tmp_7_2;
    matrix(indices[2]+1,indices[1]+0) += tmp_7_3;
    matrix(indices[2]+1,indices[1]+1) += tmp_7_4;
    matrix(indices[2]+1,indices[1]+2) += tmp_7_5;
    matrix(indices[2]+1,indices[2]+0) += tmp_7_6;
    matrix(indices[2]+1,indices[2]+1) += tmp_7_7;
    matrix(indices[2]+1,indices[2]+2) += tmp_7_8;
    matrix(indices[2]+1,indices[3]+0) += tmp_7_9;
    matrix(indices[2]+1,indices[3]+1) += tmp_7_10;
    matrix(indices[2]+1,indices[3]+2) += tmp_7_11;
    matrix(indices[2]+1,indices[4]+0) += tmp_7_12;
    matrix(indices[2]+1,indices[4]+1) += tmp_7_13;
    matrix(indices[2]+1,indices[4]+2) += tmp_7_14;
    matrix(indices[2]+1,indices[5]+0) += tmp_7_15;
    matrix(indices[2]+1,indices[5]+1) += tmp_7_16;
    matrix(indices[2]+1,indices[5]+2) += tmp_7_17;
    matrix(indices[2]+1,indices[6]+0) += tmp_7_18;
    matrix(indices[2]+1,indices[6]+1) += tmp_7_19;
    matrix(indices[2]+1,indices[6]+2) += tmp_7_20;
    matrix(indices[2]+1,indices[7]+0) += tmp_7_21;
    matrix(indices[2]+1,indices[7]+1) += tmp_7_22;
    matrix(indices[2]+1,indices[7]+2) += tmp_7_23;
    matrix(indices[2]+2,indices[0]+0) += tmp_8_0;
    matrix(indices[2]+2,indices[0]+1) += tmp_8_1;
    matrix(indices[2]+2,indices[0]+2) += tmp_8_2;
    matrix(indices[2]+2,indices[1]+0) += tmp_8_3;
    matrix(indices[2]+2,indices[1]+1) += tmp_8_4;
    matrix(indices[2]+2,indices[1]+2) += tmp_8_5;
    matrix(indices[2]+2,indices[2]+0) += tmp_8_6;
    matrix(indices[2]+2,indices[2]+1) += tmp_8_7;
    matrix(indices[2]+2,indices[2]+2) += tmp_8_8;
    matrix(indices[2]+2,indices[3]+0) += tmp_8_9;
    matrix(indices[2]+2,indices[3]+1) += tmp_8_10;
    matrix(indices[2]+2,indices[3]+2) += tmp_8_11;
    matrix(indices[2]+2,indices[4]+0) += tmp_8_12;
    matrix(indices[2]+2,indices[4]+1) += tmp_8_13;
    matrix(indices[2]+2,indices[4]+2) += tmp_8_14;
    matrix(indices[2]+2,indices[5]+0) += tmp_8_15;
    matrix(indices[2]+2,indices[5]+1) += tmp_8_16;
    matrix(indices[2]+2,indices[5]+2) += tmp_8_17;
    matrix(indices[2]+2,indices[6]+0) += tmp_8_18;
    matrix(indices[2]+2,indices[6]+1) += tmp_8_19;
    matrix(indices[2]+2,indices[6]+2) += tmp_8_20;
    matrix(indices[2]+2,indices[7]+0) += tmp_8_21;
    matrix(indices[2]+2,indices[7]+1) += tmp_8_22;
    matrix(indices[2]+2,indices[7]+2) += tmp_8_23;
    matrix(indices[3]+0,indices[0]+0) += tmp_9_0;
    matrix(indices[3]+0,indices[0]+1) += tmp_9_1;
    matrix(indices[3]+0,indices[0]+2) += tmp_9_2;
    matrix(indices[3]+0,indices[1]+0) += tmp_9_3;
    matrix(indices[3]+0,indices[1]+1) += tmp_9_4;
    matrix(indices[3]+0,indices[1]+2) += tmp_9_5;
    matrix(indices[3]+0,indices[2]+0) += tmp_9_6;
    matrix(indices[3]+0,indices[2]+1) += tmp_9_7;
    matrix(indices[3]+0,indices[2]+2) += tmp_9_8;
    matrix(indices[3]+0,indices[3]+0) += tmp_9_9;
    matrix(indices[3]+0,indices[3]+1) += tmp_9_10;
    matrix(indices[3]+0,indices[3]+2) += tmp_9_11;
    matrix(indices[3]+0,indices[4]+0) += tmp_9_12;
    matrix(indices[3]+0,indices[4]+1) += tmp_9_13;
    matrix(indices[3]+0,indices[4]+2) += tmp_9_14;
    matrix(indices[3]+0,indices[5]+0) += tmp_9_15;
    matrix(indices[3]+0,indices[5]+1) += tmp_9_16;
    matrix(indices[3]+0,indices[5]+2) += tmp_9_17;
    matrix(indices[3]+0,indices[6]+0) += tmp_9_18;
    matrix(indices[3]+0,indices[6]+1) += tmp_9_19;
    matrix(indices[3]+0,indices[6]+2) += tmp_9_20;
    matrix(indices[3]+0,indices[7]+0) += tmp_9_21;
    matrix(indices[3]+0,indices[7]+1) += tmp_9_22;
    matrix(indices[3]+0,indices[7]+2) += tmp_9_23;
    matrix(indices[3]+1,indices[0]+0) += tmp_10_0;
    matrix(indices[3]+1,indices[0]+1) += tmp_10_1;
    matrix(indices[3]+1,indices[0]+2) += tmp_10_2;
    matrix(indices[3]+1,indices[1]+0) += tmp_10_3;
    matrix(indices[3]+1,indices[1]+1) += tmp_10_4;
    matrix(indices[3]+1,indices[1]+2) += tmp_10_5;
    matrix(indices[3]+1,indices[2]+0) += tmp_10_6;
    matrix(indices[3]+1,indices[2]+1) += tmp_10_7;
    matrix(indices[3]+1,indices[2]+2) += tmp_10_8;
    matrix(indices[3]+1,indices[3]+0) += tmp_10_9;
    matrix(indices[3]+1,indices[3]+1) += tmp_10_10;
    matrix(indices[3]+1,indices[3]+2) += tmp_10_11;
    matrix(indices[3]+1,indices[4]+0) += tmp_10_12;
    matrix(indices[3]+1,indices[4]+1) += tmp_10_13;
    matrix(indices[3]+1,indices[4]+2) += tmp_10_14;
    matrix(indices[3]+1,indices[5]+0) += tmp_10_15;
    matrix(indices[3]+1,indices[5]+1) += tmp_10_16;
    matrix(indices[3]+1,indices[5]+2) += tmp_10_17;
    matrix(indices[3]+1,indices[6]+0) += tmp_10_18;
    matrix(indices[3]+1,indices[6]+1) += tmp_10_19;
    matrix(indices[3]+1,indices[6]+2) += tmp_10_20;
    matrix(indices[3]+1,indices[7]+0) += tmp_10_21;
    matrix(indices[3]+1,indices[7]+1) += tmp_10_22;
    matrix(indices[3]+1,indices[7]+2) += tmp_10_23;
    matrix(indices[3]+2,indices[0]+0) += tmp_11_0;
    matrix(indices[3]+2,indices[0]+1) += tmp_11_1;
    matrix(indices[3]+2,indices[0]+2) += tmp_11_2;
    matrix(indices[3]+2,indices[1]+0) += tmp_11_3;
    matrix(indices[3]+2,indices[1]+1) += tmp_11_4;
    matrix(indices[3]+2,indices[1]+2) += tmp_11_5;
    matrix(indices[3]+2,indices[2]+0) += tmp_11_6;
    matrix(indices[3]+2,indices[2]+1) += tmp_11_7;
    matrix(indices[3]+2,indices[2]+2) += tmp_11_8;
    matrix(indices[3]+2,indices[3]+0) += tmp_11_9;
    matrix(indices[3]+2,indices[3]+1) += tmp_11_10;
    matrix(indices[3]+2,indices[3]+2) += tmp_11_11;
    matrix(indices[3]+2,indices[4]+0) += tmp_11_12;
    matrix(indices[3]+2,indices[4]+1) += tmp_11_13;
    matrix(indices[3]+2,indices[4]+2) += tmp_11_14;
    matrix(indices[3]+2,indices[5]+0) += tmp_11_15;
    matrix(indices[3]+2,indices[5]+1) += tmp_11_16;
    matrix(indices[3]+2,indices[5]+2) += tmp_11_17;
    matrix(indices[3]+2,indices[6]+0) += tmp_11_18;
    matrix(indices[3]+2,indices[6]+1) += tmp_11_19;
    matrix(indices[3]+2,indices[6]+2) += tmp_11_20;
    matrix(indices[3]+2,indices[7]+0) += tmp_11_21;
    matrix(indices[3]+2,indices[7]+1) += tmp_11_22;
    matrix(indices[3]+2,indices[7]+2) += tmp_11_23;
    matrix(indices[4]+0,indices[0]+0) += tmp_12_0;
    matrix(indices[4]+0,indices[0]+1) += tmp_12_1;
    matrix(indices[4]+0,indices[0]+2) += tmp_12_2;
    matrix(indices[4]+0,indices[1]+0) += tmp_12_3;
    matrix(indices[4]+0,indices[1]+1) += tmp_12_4;
    matrix(indices[4]+0,indices[1]+2) += tmp_12_5;
    matrix(indices[4]+0,indices[2]+0) += tmp_12_6;
    matrix(indices[4]+0,indices[2]+1) += tmp_12_7;
    matrix(indices[4]+0,indices[2]+2) += tmp_12_8;
    matrix(indices[4]+0,indices[3]+0) += tmp_12_9;
    matrix(indices[4]+0,indices[3]+1) += tmp_12_10;
    matrix(indices[4]+0,indices[3]+2) += tmp_12_11;
    matrix(indices[4]+0,indices[4]+0) += tmp_12_12;
    matrix(indices[4]+0,indices[4]+1) += tmp_12_13;
    matrix(indices[4]+0,indices[4]+2) += tmp_12_14;
    matrix(indices[4]+0,indices[5]+0) += tmp_12_15;
    matrix(indices[4]+0,indices[5]+1) += tmp_12_16;
    matrix(indices[4]+0,indices[5]+2) += tmp_12_17;
    matrix(indices[4]+0,indices[6]+0) += tmp_12_18;
    matrix(indices[4]+0,indices[6]+1) += tmp_12_19;
    matrix(indices[4]+0,indices[6]+2) += tmp_12_20;
    matrix(indices[4]+0,indices[7]+0) += tmp_12_21;
    matrix(indices[4]+0,indices[7]+1) += tmp_12_22;
    matrix(indices[4]+0,indices[7]+2) += tmp_12_23;
    matrix(indices[4]+1,indices[0]+0) += tmp_13_0;
    matrix(indices[4]+1,indices[0]+1) += tmp_13_1;
    matrix(indices[4]+1,indices[0]+2) += tmp_13_2;
    matrix(indices[4]+1,indices[1]+0) += tmp_13_3;
    matrix(indices[4]+1,indices[1]+1) += tmp_13_4;
    matrix(indices[4]+1,indices[1]+2) += tmp_13_5;
    matrix(indices[4]+1,indices[2]+0) += tmp_13_6;
    matrix(indices[4]+1,indices[2]+1) += tmp_13_7;
    matrix(indices[4]+1,indices[2]+2) += tmp_13_8;
    matrix(indices[4]+1,indices[3]+0) += tmp_13_9;
    matrix(indices[4]+1,indices[3]+1) += tmp_13_10;
    matrix(indices[4]+1,indices[3]+2) += tmp_13_11;
    matrix(indices[4]+1,indices[4]+0) += tmp_13_12;
    matrix(indices[4]+1,indices[4]+1) += tmp_13_13;
    matrix(indices[4]+1,indices[4]+2) += tmp_13_14;
    matrix(indices[4]+1,indices[5]+0) += tmp_13_15;
    matrix(indices[4]+1,indices[5]+1) += tmp_13_16;
    matrix(indices[4]+1,indices[5]+2) += tmp_13_17;
    matrix(indices[4]+1,indices[6]+0) += tmp_13_18;
    matrix(indices[4]+1,indices[6]+1) += tmp_13_19;
    matrix(indices[4]+1,indices[6]+2) += tmp_13_20;
    matrix(indices[4]+1,indices[7]+0) += tmp_13_21;
    matrix(indices[4]+1,indices[7]+1) += tmp_13_22;
    matrix(indices[4]+1,indices[7]+2) += tmp_13_23;
    matrix(indices[4]+2,indices[0]+0) += tmp_14_0;
    matrix(indices[4]+2,indices[0]+1) += tmp_14_1;
    matrix(indices[4]+2,indices[0]+2) += tmp_14_2;
    matrix(indices[4]+2,indices[1]+0) += tmp_14_3;
    matrix(indices[4]+2,indices[1]+1) += tmp_14_4;
    matrix(indices[4]+2,indices[1]+2) += tmp_14_5;
    matrix(indices[4]+2,indices[2]+0) += tmp_14_6;
    matrix(indices[4]+2,indices[2]+1) += tmp_14_7;
    matrix(indices[4]+2,indices[2]+2) += tmp_14_8;
    matrix(indices[4]+2,indices[3]+0) += tmp_14_9;
    matrix(indices[4]+2,indices[3]+1) += tmp_14_10;
    matrix(indices[4]+2,indices[3]+2) += tmp_14_11;
    matrix(indices[4]+2,indices[4]+0) += tmp_14_12;
    matrix(indices[4]+2,indices[4]+1) += tmp_14_13;
    matrix(indices[4]+2,indices[4]+2) += tmp_14_14;
    matrix(indices[4]+2,indices[5]+0) += tmp_14_15;
    matrix(indices[4]+2,indices[5]+1) += tmp_14_16;
    matrix(indices[4]+2,indices[5]+2) += tmp_14_17;
    matrix(indices[4]+2,indices[6]+0) += tmp_14_18;
    matrix(indices[4]+2,indices[6]+1) += tmp_14_19;
    matrix(indices[4]+2,indices[6]+2) += tmp_14_20;
    matrix(indices[4]+2,indices[7]+0) += tmp_14_21;
    matrix(indices[4]+2,indices[7]+1) += tmp_14_22;
    matrix(indices[4]+2,indices[7]+2) += tmp_14_23;
    matrix(indices[5]+0,indices[0]+0) += tmp_15_0;
    matrix(indices[5]+0,indices[0]+1) += tmp_15_1;
    matrix(indices[5]+0,indices[0]+2) += tmp_15_2;
    matrix(indices[5]+0,indices[1]+0) += tmp_15_3;
    matrix(indices[5]+0,indices[1]+1) += tmp_15_4;
    matrix(indices[5]+0,indices[1]+2) += tmp_15_5;
    matrix(indices[5]+0,indices[2]+0) += tmp_15_6;
    matrix(indices[5]+0,indices[2]+1) += tmp_15_7;
    matrix(indices[5]+0,indices[2]+2) += tmp_15_8;
    matrix(indices[5]+0,indices[3]+0) += tmp_15_9;
    matrix(indices[5]+0,indices[3]+1) += tmp_15_10;
    matrix(indices[5]+0,indices[3]+2) += tmp_15_11;
    matrix(indices[5]+0,indices[4]+0) += tmp_15_12;
    matrix(indices[5]+0,indices[4]+1) += tmp_15_13;
    matrix(indices[5]+0,indices[4]+2) += tmp_15_14;
    matrix(indices[5]+0,indices[5]+0) += tmp_15_15;
    matrix(indices[5]+0,indices[5]+1) += tmp_15_16;
    matrix(indices[5]+0,indices[5]+2) += tmp_15_17;
    matrix(indices[5]+0,indices[6]+0) += tmp_15_18;
    matrix(indices[5]+0,indices[6]+1) += tmp_15_19;
    matrix(indices[5]+0,indices[6]+2) += tmp_15_20;
    matrix(indices[5]+0,indices[7]+0) += tmp_15_21;
    matrix(indices[5]+0,indices[7]+1) += tmp_15_22;
    matrix(indices[5]+0,indices[7]+2) += tmp_15_23;
    matrix(indices[5]+1,indices[0]+0) += tmp_16_0;
    matrix(indices[5]+1,indices[0]+1) += tmp_16_1;
    matrix(indices[5]+1,indices[0]+2) += tmp_16_2;
    matrix(indices[5]+1,indices[1]+0) += tmp_16_3;
    matrix(indices[5]+1,indices[1]+1) += tmp_16_4;
    matrix(indices[5]+1,indices[1]+2) += tmp_16_5;
    matrix(indices[5]+1,indices[2]+0) += tmp_16_6;
    matrix(indices[5]+1,indices[2]+1) += tmp_16_7;
    matrix(indices[5]+1,indices[2]+2) += tmp_16_8;
    matrix(indices[5]+1,indices[3]+0) += tmp_16_9;
    matrix(indices[5]+1,indices[3]+1) += tmp_16_10;
    matrix(indices[5]+1,indices[3]+2) += tmp_16_11;
    matrix(indices[5]+1,indices[4]+0) += tmp_16_12;
    matrix(indices[5]+1,indices[4]+1) += tmp_16_13;
    matrix(indices[5]+1,indices[4]+2) += tmp_16_14;
    matrix(indices[5]+1,indices[5]+0) += tmp_16_15;
    matrix(indices[5]+1,indices[5]+1) += tmp_16_16;
    matrix(indices[5]+1,indices[5]+2) += tmp_16_17;
    matrix(indices[5]+1,indices[6]+0) += tmp_16_18;
    matrix(indices[5]+1,indices[6]+1) += tmp_16_19;
    matrix(indices[5]+1,indices[6]+2) += tmp_16_20;
    matrix(indices[5]+1,indices[7]+0) += tmp_16_21;
    matrix(indices[5]+1,indices[7]+1) += tmp_16_22;
    matrix(indices[5]+1,indices[7]+2) += tmp_16_23;
    matrix(indices[5]+2,indices[0]+0) += tmp_17_0;
    matrix(indices[5]+2,indices[0]+1) += tmp_17_1;
    matrix(indices[5]+2,indices[0]+2) += tmp_17_2;
    matrix(indices[5]+2,indices[1]+0) += tmp_17_3;
    matrix(indices[5]+2,indices[1]+1) += tmp_17_4;
    matrix(indices[5]+2,indices[1]+2) += tmp_17_5;
    matrix(indices[5]+2,indices[2]+0) += tmp_17_6;
    matrix(indices[5]+2,indices[2]+1) += tmp_17_7;
    matrix(indices[5]+2,indices[2]+2) += tmp_17_8;
    matrix(indices[5]+2,indices[3]+0) += tmp_17_9;
    matrix(indices[5]+2,indices[3]+1) += tmp_17_10;
    matrix(indices[5]+2,indices[3]+2) += tmp_17_11;
    matrix(indices[5]+2,indices[4]+0) += tmp_17_12;
    matrix(indices[5]+2,indices[4]+1) += tmp_17_13;
    matrix(indices[5]+2,indices[4]+2) += tmp_17_14;
    matrix(indices[5]+2,indices[5]+0) += tmp_17_15;
    matrix(indices[5]+2,indices[5]+1) += tmp_17_16;
    matrix(indices[5]+2,indices[5]+2) += tmp_17_17;
    matrix(indices[5]+2,indices[6]+0) += tmp_17_18;
    matrix(indices[5]+2,indices[6]+1) += tmp_17_19;
    matrix(indices[5]+2,indices[6]+2) += tmp_17_20;
    matrix(indices[5]+2,indices[7]+0) += tmp_17_21;
    matrix(indices[5]+2,indices[7]+1) += tmp_17_22;
    matrix(indices[5]+2,indices[7]+2) += tmp_17_23;
    matrix(indices[6]+0,indices[0]+0) += tmp_18_0;
    matrix(indices[6]+0,indices[0]+1) += tmp_18_1;
    matrix(indices[6]+0,indices[0]+2) += tmp_18_2;
    matrix(indices[6]+0,indices[1]+0) += tmp_18_3;
    matrix(indices[6]+0,indices[1]+1) += tmp_18_4;
    matrix(indices[6]+0,indices[1]+2) += tmp_18_5;
    matrix(indices[6]+0,indices[2]+0) += tmp_18_6;
    matrix(indices[6]+0,indices[2]+1) += tmp_18_7;
    matrix(indices[6]+0,indices[2]+2) += tmp_18_8;
    matrix(indices[6]+0,indices[3]+0) += tmp_18_9;
    matrix(indices[6]+0,indices[3]+1) += tmp_18_10;
    matrix(indices[6]+0,indices[3]+2) += tmp_18_11;
    matrix(indices[6]+0,indices[4]+0) += tmp_18_12;
    matrix(indices[6]+0,indices[4]+1) += tmp_18_13;
    matrix(indices[6]+0,indices[4]+2) += tmp_18_14;
    matrix(indices[6]+0,indices[5]+0) += tmp_18_15;
    matrix(indices[6]+0,indices[5]+1) += tmp_18_16;
    matrix(indices[6]+0,indices[5]+2) += tmp_18_17;
    matrix(indices[6]+0,indices[6]+0) += tmp_18_18;
    matrix(indices[6]+0,indices[6]+1) += tmp_18_19;
    matrix(indices[6]+0,indices[6]+2) += tmp_18_20;
    matrix(indices[6]+0,indices[7]+0) += tmp_18_21;
    matrix(indices[6]+0,indices[7]+1) += tmp_18_22;
    matrix(indices[6]+0,indices[7]+2) += tmp_18_23;
    matrix(indices[6]+1,indices[0]+0) += tmp_19_0;
    matrix(indices[6]+1,indices[0]+1) += tmp_19_1;
    matrix(indices[6]+1,indices[0]+2) += tmp_19_2;
    matrix(indices[6]+1,indices[1]+0) += tmp_19_3;
    matrix(indices[6]+1,indices[1]+1) += tmp_19_4;
    matrix(indices[6]+1,indices[1]+2) += tmp_19_5;
    matrix(indices[6]+1,indices[2]+0) += tmp_19_6;
    matrix(indices[6]+1,indices[2]+1) += tmp_19_7;
    matrix(indices[6]+1,indices[2]+2) += tmp_19_8;
    matrix(indices[6]+1,indices[3]+0) += tmp_19_9;
    matrix(indices[6]+1,indices[3]+1) += tmp_19_10;
    matrix(indices[6]+1,indices[3]+2) += tmp_19_11;
    matrix(indices[6]+1,indices[4]+0) += tmp_19_12;
    matrix(indices[6]+1,indices[4]+1) += tmp_19_13;
    matrix(indices[6]+1,indices[4]+2) += tmp_19_14;
    matrix(indices[6]+1,indices[5]+0) += tmp_19_15;
    matrix(indices[6]+1,indices[5]+1) += tmp_19_16;
    matrix(indices[6]+1,indices[5]+2) += tmp_19_17;
    matrix(indices[6]+1,indices[6]+0) += tmp_19_18;
    matrix(indices[6]+1,indices[6]+1) += tmp_19_19;
    matrix(indices[6]+1,indices[6]+2) += tmp_19_20;
    matrix(indices[6]+1,indices[7]+0) += tmp_19_21;
    matrix(indices[6]+1,indices[7]+1) += tmp_19_22;
    matrix(indices[6]+1,indices[7]+2) += tmp_19_23;
    matrix(indices[6]+2,indices[0]+0) += tmp_20_0;
    matrix(indices[6]+2,indices[0]+1) += tmp_20_1;
    matrix(indices[6]+2,indices[0]+2) += tmp_20_2;
    matrix(indices[6]+2,indices[1]+0) += tmp_20_3;
    matrix(indices[6]+2,indices[1]+1) += tmp_20_4;
    matrix(indices[6]+2,indices[1]+2) += tmp_20_5;
    matrix(indices[6]+2,indices[2]+0) += tmp_20_6;
    matrix(indices[6]+2,indices[2]+1) += tmp_20_7;
    matrix(indices[6]+2,indices[2]+2) += tmp_20_8;
    matrix(indices[6]+2,indices[3]+0) += tmp_20_9;
    matrix(indices[6]+2,indices[3]+1) += tmp_20_10;
    matrix(indices[6]+2,indices[3]+2) += tmp_20_11;
    matrix(indices[6]+2,indices[4]+0) += tmp_20_12;
    matrix(indices[6]+2,indices[4]+1) += tmp_20_13;
    matrix(indices[6]+2,indices[4]+2) += tmp_20_14;
    matrix(indices[6]+2,indices[5]+0) += tmp_20_15;
    matrix(indices[6]+2,indices[5]+1) += tmp_20_16;
    matrix(indices[6]+2,indices[5]+2) += tmp_20_17;
    matrix(indices[6]+2,indices[6]+0) += tmp_20_18;
    matrix(indices[6]+2,indices[6]+1) += tmp_20_19;
    matrix(indices[6]+2,indices[6]+2) += tmp_20_20;
    matrix(indices[6]+2,indices[7]+0) += tmp_20_21;
    matrix(indices[6]+2,indices[7]+1) += tmp_20_22;
    matrix(indices[6]+2,indices[7]+2) += tmp_20_23;
    matrix(indices[7]+0,indices[0]+0) += tmp_21_0;
    matrix(indices[7]+0,indices[0]+1) += tmp_21_1;
    matrix(indices[7]+0,indices[0]+2) += tmp_21_2;
    matrix(indices[7]+0,indices[1]+0) += tmp_21_3;
    matrix(indices[7]+0,indices[1]+1) += tmp_21_4;
    matrix(indices[7]+0,indices[1]+2) += tmp_21_5;
    matrix(indices[7]+0,indices[2]+0) += tmp_21_6;
    matrix(indices[7]+0,indices[2]+1) += tmp_21_7;
    matrix(indices[7]+0,indices[2]+2) += tmp_21_8;
    matrix(indices[7]+0,indices[3]+0) += tmp_21_9;
    matrix(indices[7]+0,indices[3]+1) += tmp_21_10;
    matrix(indices[7]+0,indices[3]+2) += tmp_21_11;
    matrix(indices[7]+0,indices[4]+0) += tmp_21_12;
    matrix(indices[7]+0,indices[4]+1) += tmp_21_13;
    matrix(indices[7]+0,indices[4]+2) += tmp_21_14;
    matrix(indices[7]+0,indices[5]+0) += tmp_21_15;
    matrix(indices[7]+0,indices[5]+1) += tmp_21_16;
    matrix(indices[7]+0,indices[5]+2) += tmp_21_17;
    matrix(indices[7]+0,indices[6]+0) += tmp_21_18;
    matrix(indices[7]+0,indices[6]+1) += tmp_21_19;
    matrix(indices[7]+0,indices[6]+2) += tmp_21_20;
    matrix(indices[7]+0,indices[7]+0) += tmp_21_21;
    matrix(indices[7]+0,indices[7]+1) += tmp_21_22;
    matrix(indices[7]+0,indices[7]+2) += tmp_21_23;
    matrix(indices[7]+1,indices[0]+0) += tmp_22_0;
    matrix(indices[7]+1,indices[0]+1) += tmp_22_1;
    matrix(indices[7]+1,indices[0]+2) += tmp_22_2;
    matrix(indices[7]+1,indices[1]+0) += tmp_22_3;
    matrix(indices[7]+1,indices[1]+1) += tmp_22_4;
    matrix(indices[7]+1,indices[1]+2) += tmp_22_5;
    matrix(indices[7]+1,indices[2]+0) += tmp_22_6;
    matrix(indices[7]+1,indices[2]+1) += tmp_22_7;
    matrix(indices[7]+1,indices[2]+2) += tmp_22_8;
    matrix(indices[7]+1,indices[3]+0) += tmp_22_9;
    matrix(indices[7]+1,indices[3]+1) += tmp_22_10;
    matrix(indices[7]+1,indices[3]+2) += tmp_22_11;
    matrix(indices[7]+1,indices[4]+0) += tmp_22_12;
    matrix(indices[7]+1,indices[4]+1) += tmp_22_13;
    matrix(indices[7]+1,indices[4]+2) += tmp_22_14;
    matrix(indices[7]+1,indices[5]+0) += tmp_22_15;
    matrix(indices[7]+1,indices[5]+1) += tmp_22_16;
    matrix(indices[7]+1,indices[5]+2) += tmp_22_17;
    matrix(indices[7]+1,indices[6]+0) += tmp_22_18;
    matrix(indices[7]+1,indices[6]+1) += tmp_22_19;
    matrix(indices[7]+1,indices[6]+2) += tmp_22_20;
    matrix(indices[7]+1,indices[7]+0) += tmp_22_21;
    matrix(indices[7]+1,indices[7]+1) += tmp_22_22;
    matrix(indices[7]+1,indices[7]+2) += tmp_22_23;
    matrix(indices[7]+2,indices[0]+0) += tmp_23_0;
    matrix(indices[7]+2,indices[0]+1) += tmp_23_1;
    matrix(indices[7]+2,indices[0]+2) += tmp_23_2;
    matrix(indices[7]+2,indices[1]+0) += tmp_23_3;
    matrix(indices[7]+2,indices[1]+1) += tmp_23_4;
    matrix(indices[7]+2,indices[1]+2) += tmp_23_5;
    matrix(indices[7]+2,indices[2]+0) += tmp_23_6;
    matrix(indices[7]+2,indices[2]+1) += tmp_23_7;
    matrix(indices[7]+2,indices[2]+2) += tmp_23_8;
    matrix(indices[7]+2,indices[3]+0) += tmp_23_9;
    matrix(indices[7]+2,indices[3]+1) += tmp_23_10;
    matrix(indices[7]+2,indices[3]+2) += tmp_23_11;
    matrix(indices[7]+2,indices[4]+0) += tmp_23_12;
    matrix(indices[7]+2,indices[4]+1) += tmp_23_13;
    matrix(indices[7]+2,indices[4]+2) += tmp_23_14;
    matrix(indices[7]+2,indices[5]+0) += tmp_23_15;
    matrix(indices[7]+2,indices[5]+1) += tmp_23_16;
    matrix(indices[7]+2,indices[5]+2) += tmp_23_17;
    matrix(indices[7]+2,indices[6]+0) += tmp_23_18;
    matrix(indices[7]+2,indices[6]+1) += tmp_23_19;
    matrix(indices[7]+2,indices[6]+2) += tmp_23_20;
    matrix(indices[7]+2,indices[7]+0) += tmp_23_21;
    matrix(indices[7]+2,indices[7]+1) += tmp_23_22;
    matrix(indices[7]+2,indices[7]+2) += tmp_23_23;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1-var_inter[2]; T reg1=1-var_inter[0]; T reg2=1-var_inter[1]; T reg3=reg2*reg1; T reg4=reg2*reg0;
    T reg5=reg2*var_inter[0]; T reg6=reg0*var_inter[0]; T reg7=reg0*reg1; T reg8=elem.pos(1)[1]*reg5; T reg9=elem.pos(0)[1]*reg3;
    T reg10=elem.pos(0)[2]*reg7; T reg11=reg0*var_inter[1]; T reg12=elem.pos(1)[2]*reg6; T reg13=elem.pos(0)[1]*reg7; T reg14=elem.pos(1)[1]*reg6;
    T reg15=elem.pos(1)[1]*reg4; T reg16=elem.pos(0)[1]*reg4; T reg17=elem.pos(1)[2]*reg4; T reg18=elem.pos(0)[2]*reg4; T reg19=var_inter[0]*var_inter[1];
    T reg20=elem.pos(0)[2]*reg3; T reg21=reg5*elem.pos(1)[2]; T reg22=reg21+reg20; reg15=reg15-reg16; T reg23=reg12+reg10;
    T reg24=elem.pos(2)[2]*reg6; T reg25=elem.pos(2)[1]*reg11; T reg26=reg19*elem.pos(2)[2]; T reg27=reg14+reg13; reg17=reg17-reg18;
    T reg28=reg8+reg9; T reg29=reg19*elem.pos(2)[1]; T reg30=var_inter[1]*reg1; T reg31=elem.pos(2)[1]*reg6; T reg32=elem.pos(2)[2]*reg11;
    reg31=reg31-reg27; T reg33=var_inter[2]*reg1; T reg34=reg22+reg26; reg25=reg15+reg25; reg15=elem.pos(3)[1]*reg11;
    T reg35=elem.pos(0)[0]*reg7; T reg36=elem.pos(1)[0]*reg6; reg32=reg17+reg32; reg17=elem.pos(3)[2]*reg11; T reg37=reg28+reg29;
    T reg38=reg30*elem.pos(3)[1]; reg24=reg24-reg23; T reg39=elem.pos(1)[0]*reg4; T reg40=reg30*elem.pos(3)[2]; T reg41=elem.pos(0)[0]*reg4;
    T reg42=elem.pos(3)[2]*reg7; T reg43=elem.pos(3)[1]*reg7; T reg44=reg2*var_inter[2]; T reg45=elem.pos(4)[1]*reg3; T reg46=reg40+reg34;
    T reg47=elem.pos(4)[2]*reg3; T reg48=reg37+reg38; reg32=reg32-reg17; T reg49=elem.pos(4)[2]*reg44; T reg50=elem.pos(0)[0]*reg3;
    T reg51=elem.pos(1)[0]*reg5; T reg52=elem.pos(4)[2]*reg33; reg24=reg42+reg24; reg42=reg36+reg35; T reg53=elem.pos(2)[0]*reg6;
    T reg54=var_inter[2]*var_inter[0]; T reg55=elem.pos(4)[1]*reg33; reg43=reg31+reg43; reg25=reg25-reg15; reg31=elem.pos(4)[1]*reg44;
    reg39=reg39-reg41; T reg56=elem.pos(2)[0]*reg11; T reg57=var_inter[2]*var_inter[1]; T reg58=elem.pos(5)[1]*reg54; reg43=reg43-reg55;
    T reg59=elem.pos(3)[0]*reg11; reg47=reg47-reg46; T reg60=reg5*elem.pos(5)[2]; T reg61=reg5*elem.pos(5)[1]; reg45=reg45-reg48;
    reg56=reg39+reg56; reg24=reg24-reg52; reg39=elem.pos(5)[2]*reg54; T reg62=reg19*elem.pos(2)[0]; reg32=reg32-reg49;
    T reg63=elem.pos(5)[2]*reg44; T reg64=reg51+reg50; T reg65=elem.pos(5)[1]*reg44; T reg66=elem.pos(3)[0]*reg7; reg25=reg25-reg31;
    reg53=reg53-reg42; T reg67=reg19*elem.pos(6)[1]; T reg68=reg30*elem.pos(3)[0]; reg61=reg45+reg61; reg45=elem.pos(6)[2]*reg54;
    T reg69=elem.pos(6)[1]*reg54; reg43=reg43-reg58; T reg70=reg64+reg62; reg24=reg24-reg39; reg53=reg66+reg53;
    reg66=reg19*elem.pos(6)[2]; T reg71=elem.pos(4)[0]*reg33; reg60=reg47+reg60; reg47=elem.pos(6)[2]*reg57; reg63=reg32+reg63;
    reg65=reg25+reg65; reg25=elem.pos(6)[1]*reg57; reg32=elem.pos(4)[0]*reg44; reg56=reg56-reg59; T reg72=reg68+reg70;
    T reg73=elem.pos(4)[0]*reg3; reg25=reg65+reg25; reg47=reg63+reg47; reg63=elem.pos(7)[1]*reg57; reg65=elem.pos(7)[2]*reg57;
    T reg74=elem.pos(7)[2]*reg33; reg45=reg24+reg45; reg24=reg30*elem.pos(7)[2]; reg66=reg60+reg66; reg53=reg53-reg71;
    reg60=elem.pos(5)[0]*reg54; reg67=reg61+reg67; reg61=reg30*elem.pos(7)[1]; T reg75=elem.pos(7)[1]*reg33; reg69=reg43+reg69;
    reg43=elem.pos(5)[0]*reg44; reg56=reg56-reg32; reg24=reg66+reg24; reg61=reg67+reg61; reg66=1+(*f.m).poisson_ratio;
    reg75=reg69+reg75; reg43=reg56+reg43; reg47=reg47-reg65; reg74=reg45+reg74; reg73=reg73-reg72;
    reg45=reg5*elem.pos(5)[0]; reg56=elem.pos(6)[0]*reg54; reg53=reg53-reg60; reg67=elem.pos(6)[0]*reg57; reg25=reg25-reg63;
    reg69=elem.pos(7)[0]*reg57; reg56=reg53+reg56; reg53=elem.pos(7)[0]*reg33; T reg76=reg75*reg24; reg67=reg43+reg67;
    reg43=reg25*reg24; T reg77=reg74*reg61; T reg78=reg47*reg61; reg45=reg73+reg45; reg73=reg19*elem.pos(6)[0];
    reg66=reg66/(*f.m).elastic_modulus; reg77=reg76-reg77; reg67=reg67-reg69; reg78=reg43-reg78; reg43=reg25*reg74;
    reg76=reg47*reg75; T reg79=pow(reg66,2); T reg80=reg30*elem.pos(7)[0]; reg53=reg56+reg53; reg73=reg45+reg73;
    reg45=reg67*reg77; reg66=reg66*reg79; reg56=reg53*reg78; T reg81=(*f.m).poisson_ratio/(*f.m).elastic_modulus; reg76=reg43-reg76;
    reg80=reg73+reg80; reg43=1.0/(*f.m).elastic_modulus; reg73=reg74*reg80; T reg82=reg53*reg24; T reg83=reg80*reg76;
    reg56=reg45-reg56; reg45=reg81*reg79; T reg84=reg43*reg66; reg24=reg67*reg24; reg79=reg43*reg79;
    T reg85=reg47*reg80; reg66=reg81*reg66; reg47=reg47*reg53; reg74=reg67*reg74; T reg86=reg25*reg80;
    reg85=reg24-reg85; reg24=reg81*reg66; T reg87=reg81*reg84; T reg88=reg43*reg79; T reg89=reg81*reg45;
    reg79=reg81*reg79; reg84=reg43*reg84; T reg90=reg67*reg61; reg83=reg56+reg83; reg61=reg53*reg61;
    reg73=reg82-reg73; reg80=reg75*reg80; reg45=reg43*reg45; reg88=reg88-reg89; reg66=reg43*reg66;
    reg87=reg24+reg87; reg79=reg89+reg79; reg84=reg84-reg24; reg77=reg77/reg83; reg73=reg73/reg83;
    reg85=reg85/reg83; reg53=reg25*reg53; reg47=reg74-reg47; reg75=reg67*reg75; reg80=reg61-reg80;
    reg86=reg90-reg86; reg78=reg78/reg83; reg88=reg43*reg88; reg86=reg86/reg83; reg66=reg24+reg66;
    reg79=reg81*reg79; reg76=reg76/reg83; reg24=reg44*reg77; reg25=reg44*reg73; reg47=reg47/reg83;
    reg56=reg54*reg78; reg43=reg43*reg84; reg61=reg89+reg45; reg67=reg54*reg85; reg74=reg81*reg87;
    reg53=reg75-reg53; reg80=reg80/reg83; reg75=reg57*reg80; reg82=reg11*reg77; reg90=reg57*reg77;
    T reg91=reg11*reg80; T reg92=reg57*reg73; T reg93=reg11*reg73; T reg94=reg7*reg86; T reg95=reg7*reg78;
    T reg96=reg4*reg77; T reg97=reg44*reg80; T reg98=reg33*reg78; T reg99=reg33*reg85; T reg100=reg54*reg86;
    T reg101=reg56+reg24; T reg102=reg67+reg25; T reg103=reg33*reg86; reg53=reg53/reg83; T reg104=reg5*reg76;
    T reg105=reg81*reg66; reg74=reg43-reg74; reg43=reg6*reg85; T reg106=reg6*reg78; T reg107=reg5*reg47;
    T reg108=reg4*reg73; reg79=reg88-reg79; reg88=reg7*reg85; reg61=reg81*reg61; reg81=reg103+reg75;
    T reg109=reg6*reg86; T reg110=reg96+reg106; T reg111=reg3*reg47; T reg112=reg108-reg88; T reg113=reg5*reg53;
    reg101=reg101+reg104; T reg114=reg98-reg24; reg105=reg74-reg105; reg74=reg30*reg53; T reg115=reg19*reg53;
    T reg116=reg25-reg99; T reg117=reg108+reg43; T reg118=reg98+reg90; T reg119=reg88+reg93; T reg120=reg30*reg47;
    T reg121=reg30*reg76; T reg122=reg100+reg97; T reg123=reg95+reg82; T reg124=reg19*reg47; T reg125=reg103-reg97;
    T reg126=reg75-reg100; T reg127=reg90-reg56; T reg128=reg94+reg91; T reg129=reg67-reg92; T reg130=reg82-reg106;
    reg61=reg79-reg61; reg79=reg99+reg92; T reg131=reg43-reg93; T reg132=reg19*reg76; T reg133=reg102+reg107;
    T reg134=reg4*reg80; T reg135=reg3*reg53; T reg136=reg3*reg76; T reg137=reg95-reg96; T reg138=0.5*reg133;
    reg125=reg125+reg135; T reg139=reg91-reg109; reg122=reg113+reg122; reg116=reg116-reg111; reg130=reg130-reg132;
    reg114=reg114+reg136; T reg140=reg74+reg128; reg112=reg112+reg111; reg137=reg137-reg136; T reg141=reg107-reg117;
    T reg142=reg94-reg134; T reg143=0.5*reg101; reg129=reg129-reg124; reg127=reg127+reg132; reg126=reg126+reg115;
    reg131=reg131+reg124; reg110=reg110-reg104; T reg144=reg123+reg121; reg119=reg119+reg120; reg61=reg61/reg105;
    T reg145=reg109+reg134; reg79=reg79-reg120; T reg146=reg121-reg118; T reg147=reg74-reg81; T reg148=0.5*reg137;
    T reg149=0.5*reg147; T reg150=0.5*reg79; T reg151=0.5*reg146; T reg152=0.5*reg131; T reg153=0.5*reg116;
    reg142=reg142-reg135; T reg154=0.5*reg112; T reg155=0.5*reg140; T reg156=0.5*reg122; T reg157=0.5*reg141;
    T reg158=0.5*reg119; reg145=reg145-reg113; T reg159=0.5*reg110; T reg160=0.5*reg144; reg84=reg84/reg105;
    T reg161=reg61*reg138; T reg162=0.5*reg127; T reg163=reg61*reg143; T reg164=0.5*reg129; T reg165=0.5*reg126;
    T reg166=0.5*reg125; T reg167=0.5*reg114; T reg168=0.5*reg130; reg139=reg139-reg115; T reg169=reg84*reg122;
    T reg170=0.5*reg145; T reg171=reg61*reg151; T reg172=reg61*reg150; T reg173=reg61*reg149; T reg174=reg61*reg166;
    T reg175=reg61*reg153; T reg176=reg61*reg159; T reg177=reg61*reg167; T reg178=reg61*reg155; T reg179=reg61*reg160;
    T reg180=reg61*reg152; T reg181=reg61*reg157; T reg182=reg61*reg165; T reg183=reg61*reg162; T reg184=0.5*reg142;
    T reg185=reg61*reg168; T reg186=reg61*reg164; T reg187=reg84*reg133; T reg188=reg61*reg156; T reg189=2*reg161;
    T reg190=reg84*reg101; reg163=2*reg163; T reg191=reg61*reg158; T reg192=reg61*reg148; T reg193=0.5*reg139;
    reg87=reg87/reg105; T reg194=reg61*reg154; reg105=reg66/reg105; reg66=reg87*reg133; reg177=2*reg177;
    T reg195=reg61*reg184; T reg196=reg84*reg116; reg194=2*reg194; T reg197=reg140*reg169; T reg198=reg84*reg140;
    T reg199=reg84*reg114; reg175=2*reg175; T reg200=reg61*reg193; T reg201=reg158*reg189; reg174=2*reg174;
    reg188=2*reg188; T reg202=reg144*reg190; T reg203=reg87*reg114; T reg204=reg105*reg140; T reg205=reg84*reg141;
    T reg206=reg87*reg101; T reg207=reg84*reg125; reg171=2*reg171; T reg208=reg105*reg125; T reg209=reg105*reg147;
    reg176=2*reg176; T reg210=reg87*reg127; T reg211=reg87*reg144; T reg212=reg110*reg84; reg181=2*reg181;
    T reg213=reg84*reg112; T reg214=reg84*reg129; T reg215=reg87*reg146; reg183=2*reg183; T reg216=reg87*reg119;
    T reg217=reg84*reg79; T reg218=2*reg179; T reg219=reg84*reg127; reg186=2*reg186; T reg220=reg84*reg144;
    reg182=2*reg182; T reg221=reg84*reg142; T reg222=reg105*reg122; reg191=2*reg191; T reg223=reg145*reg84;
    T reg224=reg84*reg119; T reg225=reg84*reg139; T reg226=reg84*reg137; T reg227=reg84*reg147; T reg228=reg84*reg131;
    T reg229=reg119*reg187; T reg230=2*reg178; T reg231=reg160*reg163; reg192=2*reg192; T reg232=reg170*reg61;
    T reg233=reg105*reg126; T reg234=reg84*reg126; reg185=2*reg185; reg173=2*reg173; reg180=2*reg180;
    T reg235=reg84*reg146; T reg236=reg84*reg130; reg172=2*reg172; T reg237=reg191*reg157; T reg238=reg168*reg177;
    T reg239=reg187*reg141; T reg240=reg159*reg163; T reg241=reg159*reg218; T reg242=reg110*reg190; T reg243=reg224*reg141;
    T reg244=reg189*reg157; T reg245=reg168*reg163; T reg246=reg131*reg187; T reg247=reg112*reg224; T reg248=reg125*reg227;
    T reg249=reg196*reg141; T reg250=reg159*reg177; T reg251=reg110*reg204; T reg252=reg170*reg218; T reg253=reg159*reg171;
    T reg254=reg133*reg187; T reg255=reg143*reg163; T reg256=reg110*reg199; T reg257=reg175*reg157; T reg258=reg114*reg219;
    T reg259=reg172*reg164; T reg260=reg141*reg217; T reg261=reg131*reg196; T reg262=reg153*reg186; T reg263=reg125*reg234;
    T reg264=reg141*reg214; T reg265=reg156*reg189; T reg266=reg235*reg127; T reg267=reg159*reg183; T reg268=reg131*reg217;
    T reg269=reg153*reg189; T reg270=reg137*reg190; T reg271=reg172*reg157; T reg272=reg154*reg189; T reg273=reg228*reg141;
    T reg274=reg159*reg185; T reg275=reg101*reg219; T reg276=reg138*reg186; T reg277=reg162*reg183; T reg278=reg139*reg225;
    T reg279=reg158*reg191; T reg280=reg159*reg176; T reg281=reg129*reg214; T reg282=reg139*reg211; T reg283=reg168*reg230;
    T reg284=reg112*reg187; T reg285=reg205*reg141; T reg286=reg101*reg190; T reg287=reg138*reg189; T reg288=reg148*reg218;
    T reg289=reg110*reg219; T reg290=reg184*reg218; T reg291=reg168*reg183; T reg292=reg131*reg214; T reg293=reg186*reg157;
    T reg294=reg138*reg172; T reg295=reg101*reg235; T reg296=reg137*reg220; T reg297=reg144*reg219; T reg298=reg154*reg191;
    T reg299=reg158*reg186; T reg300=reg114*reg190; T reg301=reg168*reg171; T reg302=reg101*reg66; T reg303=reg138*reg163;
    T reg304=reg110*reg235; T reg305=reg142*reg234; T reg306=reg142*reg198; T reg307=reg167*reg171; T reg308=reg145*reg169;
    T reg309=reg116*reg217; T reg310=reg175*reg152; T reg311=reg199*reg130; T reg312=reg105*reg119; T reg313=reg105*reg79;
    T reg314=reg167*reg163; T reg315=reg148*reg230; T reg316=reg142*reg211; T reg317=reg122*reg234; T reg318=reg142*reg225;
    T reg319=reg142*reg227; T reg320=reg152*reg189; T reg321=reg145*reg207; T reg322=reg130*reg190; T reg323=reg142*reg169;
    T reg324=reg116*reg214; T reg325=reg167*reg183; T reg326=reg105*reg133; T reg327=reg145*reg227; T reg328=reg180*reg152;
    T reg329=reg236*reg130; T reg330=reg122*reg169; T reg331=reg142*reg207; T reg332=reg152*reg191; T reg333=reg130*reg220;
    T reg334=reg133*reg217; T reg335=reg105*reg116; T reg336=reg193*reg218; T reg337=reg204*reg130; T reg338=reg116*reg187;
    T reg339=reg105*reg129; T reg340=reg145*reg234; T reg341=reg171*reg143; T reg342=reg159*reg230; T reg343=reg110*reg236;
    T reg344=reg168*reg185; T reg345=reg228*reg131; T reg346=reg112*reg217; T reg347=reg180*reg157; T reg348=reg148*reg171;
    T reg349=reg114*reg235; T reg350=reg145*reg211; T reg351=reg145*reg225; T reg352=reg125*reg169; T reg353=reg153*reg172;
    T reg354=reg168*reg218; T reg355=reg131*reg224; T reg356=reg112*reg214; T reg357=reg145*reg223; T reg358=reg133*reg222;
    T reg359=reg148*reg183; T reg360=reg110*reg220; T reg361=reg105*reg131; T reg362=reg186*reg152; T reg363=reg219*reg130;
    T reg364=reg142*reg223; T reg365=reg110*reg212; T reg366=reg181*reg157; T reg367=reg116*reg196; T reg368=reg167*reg177;
    T reg369=reg105*reg141; T reg370=reg172*reg152; T reg371=reg122*reg227; T reg372=reg125*reg207; T reg373=reg235*reg130;
    T reg374=reg142*reg221; T reg375=reg145*reg198; T reg376=reg164*reg186; T reg377=reg133*reg214; T reg378=reg143*reg183;
    T reg379=reg127*reg219; T reg380=reg119*reg214; T reg381=reg148*reg176; T reg382=reg87*reg116; reg197=reg231+reg197;
    T reg383=reg144*reg208; T reg384=reg177*reg155; T reg385=reg137*reg199; T reg386=reg154*reg175; T reg387=reg126*reg227;
    reg202=reg201+reg202; T reg388=reg155*reg188; T reg389=reg137*reg204; T reg390=reg87*reg79; T reg391=reg144*reg222;
    T reg392=reg155*reg163; T reg393=reg150*reg172; T reg394=reg137*reg235; T reg395=reg154*reg172; T reg396=reg126*reg234;
    T reg397=reg87*reg141; T reg398=reg140*reg234; T reg399=reg139*reg227; T reg400=reg154*reg181; T reg401=reg144*reg220;
    reg232=2*reg232; T reg402=reg87*reg129; T reg403=reg158*reg218; T reg404=reg144*reg216; T reg405=reg137*reg219;
    T reg406=reg154*reg186; T reg407=reg182*reg160; T reg408=reg140*reg210; T reg409=reg158*reg175; T reg410=reg144*reg199;
    T reg411=reg110*reg87; T reg412=reg119*reg217; T reg413=reg183*reg160; T reg414=reg147*reg227; T reg415=reg87*reg130;
    reg231=reg229+reg231; T reg416=reg148*reg185; T reg417=reg228*reg112; T reg418=reg140*reg198; T reg419=reg194*reg154;
    T reg420=reg226*reg137; T reg421=reg177*reg160; T reg422=reg119*reg196; T reg423=reg79*reg217; T reg424=reg140*reg203;
    T reg425=reg119*reg204; T reg426=reg155*reg191; T reg427=reg105*reg142; T reg428=reg151*reg171; T reg429=reg174*reg160;
    T reg430=reg87*reg112; T reg431=reg146*reg235; T reg432=reg160*reg188; T reg433=reg144*reg233; T reg434=reg183*reg155;
    T reg435=reg140*reg206; T reg436=reg205*reg112; T reg437=reg158*reg172; reg235=reg144*reg235; T reg438=reg105*reg139;
    T reg439=reg140*reg207; T reg440=reg144*reg209; T reg441=reg87*reg131; T reg442=reg171*reg155; T reg443=reg119*reg224;
    T reg444=reg160*reg218; T reg445=reg137*reg236; T reg446=reg154*reg180; reg200=2*reg200; T reg447=reg171*reg160;
    reg217=reg129*reg217; T reg448=reg139*reg169; reg227=reg140*reg227; T reg449=reg162*reg171; T reg450=reg148*reg177;
    T reg451=reg112*reg196; T reg452=reg112*reg213; T reg453=reg192*reg148; T reg454=reg139*reg207; T reg455=reg153*reg175;
    T reg456=reg114*reg199; reg195=2*reg195; T reg457=reg137*reg212; T reg458=reg173*reg160; T reg459=reg139*reg198;
    T reg460=reg139*reg234; T reg461=reg148*reg163; T reg462=reg140*reg215; T reg463=reg145*reg105; T reg464=reg166*reg188;
    T reg465=reg463*reg141; T reg466=reg200*reg157; T reg467=reg145*reg361; reg454=reg238+reg454; reg351=reg274+reg351;
    T reg468=reg125*reg326; T reg469=reg153*reg188; T reg470=reg140*reg335; reg442=reg440+reg442; reg268=reg301+reg268;
    T reg471=reg170*reg181; reg275=reg275-reg276; T reg472=reg350+reg342; T reg473=reg230*reg157; T reg474=reg153*reg182;
    T reg475=reg186*reg160; reg412=reg412-reg447; reg260=reg253+reg260; T reg476=reg139*reg335; T reg477=reg170*reg172;
    T reg478=reg209*reg141; T reg479=reg167*reg182; T reg480=reg119*reg203; reg443=reg443+reg444; T reg481=reg166*reg174;
    T reg482=reg125*reg210; T reg483=reg230*reg152; reg357=reg280+reg357; T reg484=reg145*reg415; T reg485=reg159*reg200;
    reg352=reg314+reg352; T reg486=reg116*reg209; reg235=reg437-reg235; T reg487=reg175*reg160; T reg488=reg166*reg172;
    T reg489=reg145*reg206; T reg490=reg159*reg188; T reg491=reg188*reg157; reg227=reg447+reg227; reg447=reg145*reg326;
    reg309=reg307+reg309; T reg492=reg152*reg188; reg308=reg240+reg308; reg434=reg433+reg434; T reg493=reg145*reg210;
    T reg494=reg139*reg312; T reg495=reg159*reg182; T reg496=reg116*reg215; T reg497=reg145*reg312; T reg498=reg167*reg188;
    T reg499=reg156*reg182; T reg500=reg139*reg206; T reg501=reg125*reg206; T reg502=reg144*reg390; T reg503=reg241+reg375;
    T reg504=reg145*reg203; T reg505=reg159*reg174; reg372=reg368+reg372; T reg506=reg158*reg171; T reg507=reg174*reg157;
    reg439=reg421+reg439; T reg508=reg168*reg188; T reg509=reg145*reg335; T reg510=reg173*reg155; reg321=reg250+reg321;
    T reg511=reg114*reg208; T reg512=reg354+reg459; reg286=reg286+reg287; T reg513=reg170*reg191; T reg514=reg204*reg141;
    reg429=reg424+reg429; reg421=reg422-reg421; reg422=reg415*reg141; reg278=reg344+reg278; reg424=reg139*reg203;
    T reg515=reg119*reg208; T reg516=reg159*reg175; T reg517=reg203*reg141; T reg518=reg175*reg155; T reg519=reg444+reg418;
    reg249=reg250+reg249; reg248=reg307+reg248; reg250=reg156*reg163; reg307=reg101*reg222; T reg520=reg170*reg180;
    T reg521=reg438*reg141; reg426=reg425+reg426; reg273=reg274+reg273; reg274=reg114*reg382; T reg522=reg159*reg191;
    T reg523=reg211*reg141; reg303=reg302+reg303; T reg524=reg153*reg177; T reg525=reg158*reg174; T reg526=reg166*reg177;
    T reg527=reg282+reg283; reg297=reg299-reg297; T reg528=reg182*reg155; T reg529=reg156*reg188; T reg530=reg141*reg210;
    T reg531=reg159*reg180; T reg532=reg193*reg172; reg264=reg267+reg264; reg263=reg325+reg263; T reg533=reg388+reg231;
    T reg534=reg170*reg186; T reg535=reg172*reg155; T reg536=reg119*reg209; T reg537=reg119*reg222; T reg538=reg155*reg189;
    T reg539=reg233*reg141; T reg540=reg159*reg172; T reg541=reg119*reg210; T reg542=reg174*reg152; T reg543=reg141*reg215;
    T reg544=reg125*reg339; T reg545=reg168*reg174; T reg546=reg170*reg175; T reg547=reg208*reg141; T reg548=reg159*reg189;
    T reg549=reg206*reg141; T reg550=reg125*reg313; T reg551=reg153*reg173; T reg552=reg119*reg206; T reg553=reg131*reg209;
    reg456=reg455+reg456; T reg554=reg160*reg189; reg240=reg240-reg239; T reg555=reg170*reg189; T reg556=reg222*reg141;
    T reg557=reg167*reg173; T reg558=reg125*reg215; T reg559=reg159*reg186; T reg560=reg233*reg130; T reg561=reg114*reg209;
    T reg562=reg166*reg171; T reg563=reg193*reg173; T reg564=reg182*reg152; T reg565=reg139*reg339; reg373=reg370+reg373;
    T reg566=reg171*reg152; reg458=reg462+reg458; reg462=reg155*reg218; T reg567=reg390*reg130; T reg568=reg114*reg390;
    T reg569=reg153*reg171; T reg570=reg193*reg171; T reg571=reg130*reg209; T reg572=reg144*reg204; reg407=reg408+reg407;
    reg408=reg166*reg182; T reg573=reg152*reg163; T reg574=reg130*reg66; T reg575=reg166*reg175; T reg576=reg193*reg163;
    T reg577=reg222*reg130; reg367=reg368+reg367; reg368=reg193*reg182; T reg578=reg168*reg182; T reg579=reg119*reg233;
    reg363=reg362+reg363; T reg580=reg186*reg155; T reg581=reg172*reg160; T reg582=reg158*reg177; T reg583=reg174*reg155;
    T reg584=reg183*reg152; T reg585=reg402*reg130; reg410=reg409-reg410; T reg586=reg193*reg183; reg258=reg262+reg258;
    T reg587=reg193*reg191; T reg588=reg204*reg131; T reg589=reg168*reg175; T reg590=reg131*reg203; T reg591=reg140*reg339;
    T reg592=reg114*reg402; reg399=reg301+reg399; reg301=reg153*reg183; T reg593=reg131*reg208; T reg594=reg139*reg313;
    reg261=reg238+reg261; reg238=reg119*reg215; T reg595=reg193*reg175; T reg596=reg173*reg152; reg398=reg413+reg398;
    T reg597=reg168*reg173; T reg598=reg139*reg215; reg404=reg403+reg404; T reg599=reg131*reg206; reg345=reg344+reg345;
    reg349=reg353+reg349; reg344=reg158*reg182; T reg600=reg193*reg180; T reg601=reg390*reg127; T reg602=reg438*reg131;
    T reg603=reg168*reg189; T reg604=reg166*reg173; T reg605=reg168*reg191; T reg606=reg131*reg211; T reg607=reg230*reg155;
    T reg608=reg279+reg401; reg355=reg355-reg354; T reg609=reg114*reg233; T reg610=reg166*reg183; reg460=reg291+reg460;
    reg327=reg253+reg327; reg392=reg391+reg392; reg253=reg233*reg131; T reg611=reg193*reg200; T reg612=reg193*reg186;
    reg448=reg245+reg448; T reg613=reg153*reg163; reg329=reg328+reg329; T reg614=reg185*reg152; T reg615=reg441*reg130;
    T reg616=reg116*reg210; T reg617=reg158*reg188; T reg618=reg167*reg186; T reg619=reg193*reg185; T reg620=reg438*reg130;
    T reg621=reg144*reg66; T reg622=reg158*reg163; T reg623=reg193*reg230; T reg624=reg182*reg157; T reg625=reg139*reg326;
    T reg626=reg145*reg339; T reg627=reg167*reg172; reg432=reg435+reg432; reg300=reg300-reg269; reg340=reg267+reg340;
    reg267=reg144*reg402; reg435=reg116*reg233; T reg628=reg145*reg215; T reg629=reg159*reg173; T reg630=reg166*reg186;
    T reg631=reg173*reg157; T reg632=reg145*reg313; reg324=reg325+reg324; reg325=reg131*reg215; T reg633=reg168*reg172;
    T reg634=reg158*reg183; T reg635=reg166*reg163; T reg636=reg177*reg152; T reg637=reg382*reg130; T reg638=reg116*reg206;
    T reg639=reg193*reg177; T reg640=reg139*reg210; reg197=reg201+reg197; T reg641=reg144*reg382; T reg642=reg130*reg208;
    T reg643=reg167*reg189; T reg644=reg193*reg188; T reg645=reg114*reg222; T reg646=reg222*reg131; T reg647=reg193*reg189;
    reg413=reg380-reg413; reg322=reg322-reg320; reg380=reg116*reg208; reg245=reg245-reg246; T reg648=reg140*reg326;
    reg388=reg202+reg388; T reg649=reg332-reg333; T reg650=reg116*reg222; T reg651=reg166*reg189; T reg652=reg152*reg218;
    T reg653=reg130*reg216; T reg654=reg114*reg66; reg292=reg291+reg292; reg291=reg336+reg337; reg314=reg314-reg338;
    T reg655=reg140*reg313; T reg656=reg158*reg173; T reg657=reg193*reg174; reg384=reg383+reg384; T reg658=reg131*reg210;
    reg311=reg310+reg311; T reg659=reg168*reg186; T reg660=reg165*reg183; T reg661=reg233*reg127; reg243=reg243-reg241;
    T reg662=reg112*reg210; reg356=reg359+reg356; T reg663=reg127*reg402; T reg664=reg164*reg183; T reg665=reg184*reg186;
    T reg666=reg112*reg233; T reg667=reg148*reg172; T reg668=reg112*reg215; T reg669=reg165*reg182; reg379=reg376+reg379;
    reg346=reg348+reg346; T reg670=reg184*reg172; T reg671=reg112*reg209; reg374=reg453+reg374; reg371=reg341+reg371;
    T reg672=reg411*reg142; T reg673=reg148*reg232; T reg674=reg154*reg232; T reg675=reg184*reg188; reg270=reg270-reg272;
    T reg676=reg171*reg165; T reg677=reg154*reg163; T reg678=reg184*reg230; T reg679=reg298-reg296; T reg680=reg127*reg209;
    T reg681=reg148*reg186; T reg682=reg154*reg218; T reg683=reg137*reg216; T reg684=reg171*reg164; T reg685=reg184*reg180;
    T reg686=reg148*reg191; T reg687=reg112*reg211; T reg688=reg173*reg165; reg247=reg247-reg288; reg266=reg259+reg266;
    T reg689=reg184*reg191; T reg690=reg112*reg204; T reg691=reg148*reg175; T reg692=reg112*reg222; T reg693=reg122*reg210;
    T reg694=reg154*reg174; reg331=reg450+reg331; reg330=reg255+reg330; T reg695=reg142*reg206; T reg696=reg148*reg188;
    T reg697=reg154*reg188; T reg698=reg142*reg326; T reg699=reg133*reg209; T reg700=reg156*reg172; reg323=reg461+reg323;
    T reg701=reg142*reg210; T reg702=reg148*reg182; reg334=reg341-reg334; reg341=reg154*reg182; T reg703=reg142*reg339;
    reg305=reg359+reg305; reg359=reg133*reg215; T reg704=reg142*reg215; T reg705=reg148*reg173; T reg706=reg172*reg143;
    T reg707=reg142*reg369; T reg708=reg138*reg173; T reg709=reg122*reg313; reg364=reg381+reg364; T reg710=reg415*reg142;
    T reg711=reg148*reg200; T reg712=reg173*reg143; T reg713=reg122*reg215; T reg714=reg154*reg200; T reg715=reg142*reg361;
    reg318=reg416+reg318; reg317=reg378+reg317; T reg716=reg316+reg315; T reg717=reg154*reg230; T reg718=reg142*reg312;
    T reg719=reg138*reg182; T reg720=reg122*reg339; T reg721=reg288+reg306; T reg722=reg142*reg203; T reg723=reg148*reg174;
    T reg724=reg143*reg182; T reg725=reg184*reg200; reg445=reg446+reg445; T reg726=reg146*reg209; T reg727=reg154*reg185;
    T reg728=reg137*reg441; T reg729=reg146*reg390; T reg730=reg184*reg185; T reg731=reg137*reg438; T reg732=reg150*reg171;
    T reg733=reg149*reg173; T reg734=reg184*reg183; T reg735=reg137*reg233; reg431=reg393+reg431; T reg736=reg184*reg173;
    reg394=reg395+reg394; T reg737=reg154*reg171; T reg738=reg137*reg390; reg387=reg449+reg387; T reg739=reg184*reg171;
    T reg740=reg389+reg290; T reg741=reg184*reg174; T reg742=reg112*reg203; T reg743=reg142*reg335; T reg744=reg148*reg181;
    T reg745=reg411*reg112; reg436=reg381+reg436; reg414=reg428+reg414; reg381=reg184*reg181; T reg746=reg112*reg463;
    T reg747=reg148*reg180; T reg748=reg415*reg112; T reg749=reg79*reg209; T reg750=reg149*reg172; reg417=reg416+reg417;
    reg416=reg184*reg195; reg420=reg419+reg420; reg423=reg428+reg423; reg428=reg192*reg154; T reg751=reg430*reg137;
    T reg752=reg192*reg184; T reg753=reg427*reg137; T reg754=reg149*reg171; T reg755=reg112*reg438; T reg756=reg172*reg165;
    T reg757=reg137*reg209; reg452=reg453+reg452; reg217=reg449+reg217; reg449=reg194*reg184; reg453=reg427*reg112;
    reg451=reg450+reg451; reg450=reg129*reg215; T reg758=reg162*reg172; T reg759=reg184*reg175; T reg760=reg112*reg208;
    T reg761=reg148*reg189; T reg762=reg112*reg206; T reg763=reg233*reg129; T reg764=reg165*reg186; reg461=reg461-reg284;
    T reg765=reg184*reg189; reg281=reg277+reg281; T reg766=reg184*reg177; T reg767=reg137*reg208; reg385=reg386+reg385;
    T reg768=reg126*reg313; T reg769=reg154*reg177; T reg770=reg137*reg382; T reg771=reg137*reg66; T reg772=reg173*reg164;
    T reg773=reg184*reg163; T reg774=reg137*reg222; T reg775=reg162*reg173; T reg776=reg184*reg182; reg405=reg406+reg405;
    T reg777=reg126*reg215; T reg778=reg154*reg183; T reg779=reg137*reg402; T reg780=reg184*reg232; reg457=reg457+reg400;
    reg396=reg277+reg396; reg277=reg154*reg176; T reg781=reg137*reg397; T reg782=reg184*reg176; T reg783=reg137*reg463;
    T reg784=reg129*reg209; T reg785=reg170*reg182; reg289=reg289+reg293; reg304=reg304+reg271; T reg786=reg170*reg177;
    T reg787=reg170*reg173; T reg788=reg138*reg171; T reg789=reg170*reg232; T reg790=reg110*reg208; reg365=reg365+reg366;
    reg295=reg295-reg294; T reg791=reg101*reg209; T reg792=reg156*reg173; T reg793=reg110*reg390; T reg794=reg110*reg216;
    T reg795=reg171*reg157; T reg796=reg218*reg157; T reg797=reg156*reg171; T reg798=reg170*reg185; reg242=reg242-reg244;
    T reg799=reg183*reg157; T reg800=reg110*reg463; T reg801=reg170*reg188; T reg802=reg110*reg233; T reg803=reg176*reg157;
    T reg804=reg110*reg397; T reg805=reg170*reg176; T reg806=reg170*reg183; T reg807=reg265+reg358; T reg808=reg237-reg360;
    T reg809=reg156*reg183; T reg810=reg170*reg230; reg377=reg378-reg377; reg378=reg101*reg233; T reg811=reg110*reg402;
    T reg812=reg110*reg438; T reg813=reg110*reg222; T reg814=reg110*reg66; T reg815=reg251+reg252; T reg816=reg163*reg157;
    T reg817=reg143*reg186; reg171=reg170*reg171; T reg818=reg133*reg210; reg255=reg255+reg254; reg209=reg110*reg209;
    T reg819=reg110*reg441; T reg820=reg110*reg382; T reg821=reg185*reg157; T reg822=reg101*reg402; T reg823=reg133*reg233;
    T reg824=reg138*reg183; T reg825=reg156*reg186; reg285=reg280+reg285; reg280=reg177*reg157; reg343=reg343+reg347;
    reg390=reg101*reg390; reg256=reg256+reg257; T reg826=reg170*reg200; T reg827=reg154*reg173; T reg828=reg170*reg174;
    T reg829=reg142*reg313; T reg830=reg170*reg163; reg319=reg348+reg319; reg816=reg816-reg814; reg685=reg755+reg685;
    reg578=reg640+reg578; reg779=reg778+reg779; reg348=reg83*reg458; reg249=reg828+reg249; reg565=reg564+reg565;
    reg549=reg549-reg548; reg242=reg242+reg801; reg784=reg756+reg784; reg591=reg344-reg591; reg783=reg782+reg783;
    reg399=reg370+reg399; reg460=reg362+reg460; reg550=reg551+reg550; reg597=reg598+reg597; reg398=reg299-reg398;
    reg457=reg780+reg457; reg396=reg376+reg396; reg781=reg277+reg781; reg594=reg596+reg594; reg547=reg546+reg547;
    reg332=reg332-reg512; reg762=reg762-reg761; reg280=reg820+reg280; reg494=reg494-reg483; reg530=reg559+reg530;
    reg461=reg461+reg675; reg274=reg524+reg274; reg277=reg83*reg527; reg281=reg669+reg281; reg263=reg262+reg263;
    reg278=reg328+reg278; reg797=reg791+reg797; reg767=reg766+reg767; reg511=reg526+reg511; reg553=reg532+reg553;
    reg264=reg785+reg264; reg270=reg675+reg270; reg268=reg563+reg268; reg676=reg680+reg676; reg828=reg256+reg828;
    reg256=reg83*reg815; reg757=reg739+reg757; reg217=reg688+reg217; reg448=reg448-reg320; reg788=reg390-reg788;
    reg452=reg416+reg452; reg655=reg656-reg655; reg492=reg492-reg625; reg786=reg790+reg786; reg240=reg801+reg240;
    reg453=reg449+reg453; reg508=reg500+reg508; reg450=reg758+reg450; reg227=reg437-reg227; reg454=reg310+reg454;
    reg451=reg741+reg451; reg557=reg558+reg557; reg476=reg542+reg476; reg556=reg556-reg555; reg760=reg759+reg760;
    reg545=reg424+reg545; reg763=reg764+reg763; reg456=reg481+reg456; reg518=reg515-reg518; reg423=reg733+reg423;
    reg250=reg307+reg250; reg279=reg279+reg519; reg421=reg421-reg583; reg420=reg416+reg420; reg824=reg822-reg824;
    reg273=reg826+reg273; reg487=reg480-reg487; reg751=reg428+reg751; reg262=reg83*reg429; reg299=reg83*reg426;
    reg304=reg304+reg787; reg753=reg752+reg753; reg754=reg726+reg754; reg443=reg607+reg443; reg521=reg520+reg521;
    reg445=reg725+reg445; reg470=reg525-reg470; reg307=reg83*reg442; reg310=reg83*reg303; reg728=reg727+reg728;
    reg581=reg238-reg581; reg580=reg579-reg580; reg285=reg789+reg285; reg275=reg275+reg499; reg413=reg413-reg528;
    reg745=reg744+reg745; reg414=reg393+reg414; reg171=reg209+reg171; reg780=reg436+reg780; reg475=reg541-reg475;
    reg465=reg471+reg465; reg412=reg412-reg510; reg537=reg537+reg538; reg746=reg381+reg746; reg749=reg750+reg749;
    reg209=reg83*reg533; reg748=reg747+reg748; reg795=reg793+reg795; reg535=reg536-reg535; reg552=reg552+reg554;
    reg422=reg531+reg422; reg725=reg417+reg725; reg238=reg83*reg388; reg785=reg289+reg785; reg289=reg83*reg740;
    reg328=reg83*reg384; reg385=reg741+reg385; reg768=reg772+reg768; reg513=reg513-reg514; reg641=reg582-reg641;
    reg770=reg769+reg770; reg830=reg813+reg830; reg344=reg83*reg197; reg583=reg410-reg583; reg677=reg677-reg771;
    reg517=reg516+reg517; reg774=reg773+reg774; reg362=reg572+reg462; reg775=reg777+reg775; reg248=reg353+reg248;
    reg353=reg83*reg407; reg370=reg83*reg404; reg295=reg295+reg792; reg405=reg776+reg405; reg608=reg608+reg607;
    reg502=reg506-reg502; reg729=reg732+reg729; reg522=reg522-reg523; reg439=reg409-reg439; reg510=reg235-reg510;
    reg806=reg802+reg806; reg731=reg730+reg731; reg733=reg431+reg733; reg235=reg83*reg434; reg809=reg378+reg809;
    reg735=reg734+reg735; reg799=reg811+reg799; reg267=reg634-reg267; reg528=reg297-reg528; reg297=reg83*reg432;
    reg376=reg83*reg392; reg394=reg736+reg394; reg286=reg286+reg529; reg622=reg622+reg621; reg738=reg737+reg738;
    reg387=reg259+reg387; reg243=reg243-reg810; reg617=reg617+reg648; reg259=reg83*reg472; reg573=reg573-reg574;
    reg498=reg501+reg498; reg805=reg800+reg805; reg380=reg575+reg380; reg322=reg644+reg322; reg715=reg714+reg715;
    reg276=reg317-reg276; reg318=reg446+reg318; reg642=reg639+reg642; reg497=reg497-reg473; reg317=reg83*reg716;
    reg638=reg638-reg643; reg637=reg636+reg637; reg719=reg720-reg719; reg803=reg804+reg803; reg311=reg657+reg311;
    reg718=reg718-reg717; reg237=reg237-reg503; reg298=reg298-reg721; reg378=reg83*reg291; reg724=reg693+reg724;
    reg671=reg670+reg671; reg567=reg566+reg567; reg294=reg371-reg294; reg467=reg466+reg467; reg568=reg569+reg568;
    reg373=reg563+reg373; reg374=reg419+reg374; reg469=reg469-reg468; reg673=reg672+reg673; reg560=reg586+reg560;
    reg708=reg709-reg708; reg561=reg562+reg561; reg585=reg584+reg585; reg351=reg347+reg351; reg707=reg674+reg707;
    reg826=reg343+reg826; reg363=reg368+reg363; reg364=reg400+reg364; reg712=reg713+reg712; reg818=reg817-reg818;
    reg367=reg481+reg367; reg577=reg576+reg577; reg711=reg710+reg711; reg321=reg257+reg321; reg323=reg323-reg272;
    reg324=reg408+reg324; reg632=reg631+reg632; reg334=reg792+reg334; reg629=reg628+reg629; reg702=reg701+reg702;
    reg825=reg825-reg823; reg435=reg630+reg435; reg340=reg293+reg340; reg829=reg827+reg829; reg490=reg489+reg490;
    reg703=reg341+reg703; reg626=reg624+reg626; reg305=reg406+reg305; reg496=reg627+reg496; reg495=reg493+reg495;
    reg359=reg706-reg359; reg309=reg604+reg309; reg308=reg308-reg244; reg705=reg704+reg705; reg491=reg491-reg447;
    reg372=reg455+reg372; reg314=reg464+reg314; reg653=reg653-reg652; reg723=reg722+reg723; reg377=reg499+reg377;
    reg505=reg504+reg505; reg649=reg649-reg623; reg789=reg365+reg789; reg694=reg743+reg694; reg330=reg287+reg330;
    reg650=reg650-reg651; reg620=reg619+reg620; reg331=reg386+reg331; reg615=reg614+reg615; reg509=reg507+reg509;
    reg696=reg695+reg696; reg616=reg618+reg616; reg329=reg611+reg329; reg700=reg700-reg699; reg486=reg488+reg486;
    reg697=reg697-reg698; reg327=reg271+reg327; reg319=reg395+reg319; reg255=reg529+reg255; reg478=reg477+reg478;
    reg590=reg589+reg590; reg662=reg681+reg662; reg592=reg301+reg592; reg543=reg540+reg543; reg587=reg587-reg588;
    reg658=reg659+reg658; reg663=reg664+reg663; reg798=reg812+reg798; reg356=reg776+reg356; reg613=reg613-reg654;
    reg355=reg355-reg623; reg684=reg601+reg684; reg357=reg366+reg357; reg683=reg683-reg682; reg292=reg368+reg292;
    reg260=reg787+reg260; reg599=reg599-reg603; reg247=reg247-reg678; reg689=reg689-reg690; reg245=reg644+reg245;
    reg479=reg482+reg479; reg593=reg595+reg593; reg691=reg742+reg691; reg660=reg661+reg660; reg258=reg408+reg258;
    reg645=reg635+reg645; reg688=reg266+reg688; reg686=reg686-reg687; reg646=reg646-reg647; reg261=reg657+reg261;
    reg808=reg808-reg810; reg692=reg692-reg765; reg794=reg794-reg796; reg253=reg612+reg253; reg345=reg611+reg345;
    reg485=reg484+reg485; reg668=reg667+reg668; reg300=reg464+reg300; reg257=reg83*reg807; reg539=reg534+reg539;
    reg821=reg819+reg821; reg346=reg736+reg346; reg349=reg604+reg349; reg602=reg600+reg602; reg679=reg679-reg678;
    reg669=reg379+reg669; reg325=reg633+reg325; reg571=reg570+reg571; reg609=reg610+reg609; reg666=reg665+reg666;
    reg352=reg352-reg269; reg544=reg474+reg544; reg605=reg605-reg606; reg266=ponderation*reg262; reg300=reg83*reg300;
    reg700=reg83*reg700; reg275=reg83*reg275; reg616=reg83*reg616; reg581=reg83*reg581; reg496=reg83*reg496;
    reg645=reg83*reg645; reg825=reg83*reg825; reg255=reg83*reg255; reg754=reg83*reg754; reg824=reg83*reg824;
    reg456=reg83*reg456; reg377=reg83*reg377; reg511=reg83*reg511; reg749=reg83*reg749; reg412=reg83*reg412;
    reg613=reg83*reg613; reg324=reg83*reg324; reg435=reg83*reg435; reg676=reg83*reg676; reg334=reg83*reg334;
    reg535=reg83*reg535; reg263=reg83*reg263; reg797=reg83*reg797; reg250=reg83*reg250; reg281=reg83*reg281;
    reg423=reg83*reg423; reg688=reg83*reg688; reg359=reg83*reg359; reg309=reg83*reg309; reg274=reg83*reg274;
    reg544=reg83*reg544; reg684=reg83*reg684; reg414=reg83*reg414; reg279=reg83*reg279; reg486=reg83*reg486;
    reg217=reg83*reg217; reg498=reg83*reg498; reg768=reg83*reg768; reg617=reg83*reg617; reg367=reg83*reg367;
    reg271=ponderation*reg348; reg818=reg83*reg818; reg352=reg83*reg352; reg609=reg83*reg609; reg712=reg83*reg712;
    reg293=ponderation*reg344; reg669=reg83*reg669; reg775=reg83*reg775; reg248=reg83*reg248; reg561=reg83*reg561;
    reg398=reg83*reg398; reg784=reg83*reg784; reg708=reg83*reg708; reg469=reg83*reg469; reg301=ponderation*reg257;
    reg295=reg83*reg295; reg341=ponderation*reg353; reg396=reg83*reg396; reg550=reg83*reg550; reg568=reg83*reg568;
    reg591=reg83*reg591; reg294=reg83*reg294; reg349=reg83*reg349; reg650=reg83*reg650; reg330=reg83*reg330;
    reg479=reg83*reg479; reg343=ponderation*reg310; reg763=reg83*reg763; reg729=reg83*reg729; reg660=reg83*reg660;
    reg470=reg83*reg470; reg314=reg83*reg314; reg557=reg83*reg557; reg724=reg83*reg724; reg258=reg83*reg258;
    reg227=reg83*reg227; reg372=reg83*reg372; reg733=reg83*reg733; reg439=reg83*reg439; reg450=reg83*reg450;
    reg809=reg83*reg809; reg638=reg83*reg638; reg719=reg83*reg719; reg788=reg83*reg788; reg655=reg83*reg655;
    reg387=reg83*reg387; reg286=reg83*reg286; reg347=ponderation*reg297; reg592=reg83*reg592; reg380=reg83*reg380;
    reg276=reg83*reg276; reg663=reg83*reg663; reg662=reg83*reg662; reg649=reg83*reg649; reg356=reg83*reg356;
    reg620=reg83*reg620; reg615=reg83*reg615; reg666=reg83*reg666; reg329=reg83*reg329; reg668=reg83*reg668;
    reg327=reg83*reg327; reg632=reg83*reg632; reg346=reg83*reg346; reg629=reg83*reg629; reg340=reg83*reg340;
    reg671=reg83*reg671; reg626=reg83*reg626; reg374=reg83*reg374; reg495=reg83*reg495; reg308=reg83*reg308;
    reg673=reg83*reg673; reg491=reg83*reg491; reg490=reg83*reg490; reg707=reg83*reg707; reg321=reg83*reg321;
    reg509=reg83*reg509; reg364=reg83*reg364; reg505=reg83*reg505; reg711=reg83*reg711; reg237=reg83*reg237;
    reg497=reg83*reg497; reg715=reg83*reg715; reg365=ponderation*reg259; reg590=reg83*reg590; reg461=reg83*reg461;
    reg587=reg83*reg587; reg770=reg83*reg770; reg355=reg83*reg355; reg605=reg83*reg605; reg767=reg83*reg767;
    reg602=reg83*reg602; reg345=reg83*reg345; reg270=reg83*reg270; reg571=reg83*reg571; reg679=reg83*reg679;
    reg567=reg83*reg567; reg373=reg83*reg373; reg683=reg83*reg683; reg560=reg83*reg560; reg585=reg83*reg585;
    reg363=reg83*reg363; reg686=reg83*reg686; reg577=reg83*reg577; reg247=reg83*reg247; reg573=reg83*reg573;
    reg322=reg83*reg322; reg689=reg83*reg689; reg642=reg83*reg642; reg637=reg83*reg637; reg691=reg83*reg691;
    reg692=reg83*reg692; reg311=reg83*reg311; reg366=ponderation*reg378; reg653=reg83*reg653; reg522=reg83*reg522;
    reg305=reg83*reg305; reg521=reg83*reg521; reg705=reg83*reg705; reg273=reg83*reg273; reg422=reg83*reg422;
    reg829=reg83*reg829; reg465=reg83*reg465; reg285=reg83*reg285; reg319=reg83*reg319; reg171=reg83*reg171;
    reg795=reg83*reg795; reg789=reg83*reg789; reg304=reg83*reg304; reg803=reg83*reg803; reg806=reg83*reg806;
    reg799=reg83*reg799; reg805=reg83*reg805; reg785=reg83*reg785; reg830=reg83*reg830; reg826=reg83*reg826;
    reg816=reg83*reg816; reg242=reg83*reg242; reg821=reg83*reg821; reg786=reg83*reg786; reg798=reg83*reg798;
    reg280=reg83*reg280; reg828=reg83*reg828; reg808=reg83*reg808; reg368=ponderation*reg256; reg794=reg83*reg794;
    reg318=reg83*reg318; reg351=reg83*reg351; reg467=reg83*reg467; reg371=ponderation*reg317; reg485=reg83*reg485;
    reg357=reg83*reg357; reg718=reg83*reg718; reg478=reg83*reg478; reg260=reg83*reg260; reg298=reg83*reg298;
    reg543=reg83*reg543; reg723=reg83*reg723; reg539=reg83*reg539; reg264=reg83*reg264; reg694=reg83*reg694;
    reg530=reg83*reg530; reg331=reg83*reg331; reg556=reg83*reg556; reg240=reg83*reg240; reg696=reg83*reg696;
    reg549=reg83*reg549; reg547=reg83*reg547; reg697=reg83*reg697; reg249=reg83*reg249; reg517=reg83*reg517;
    reg323=reg83*reg323; reg513=reg83*reg513; reg702=reg83*reg702; reg243=reg83*reg243; reg528=reg83*reg528;
    reg703=reg83*reg703; reg502=reg83*reg502; reg565=reg83*reg565; reg677=reg83*reg677; reg753=reg83*reg753;
    reg578=reg83*reg578; reg379=ponderation*reg307; reg448=reg83*reg448; reg774=reg83*reg774; reg492=reg83*reg492;
    reg443=reg83*reg443; reg405=reg83*reg405; reg751=reg83*reg751; reg508=reg83*reg508; reg454=reg83*reg454;
    reg381=ponderation*reg299; reg457=reg83*reg457; reg420=reg83*reg420; reg476=reg83*reg476; reg487=reg83*reg487;
    reg545=reg83*reg545; reg781=reg83*reg781; reg332=reg83*reg332; reg725=reg83*reg725; reg421=reg83*reg421;
    reg386=ponderation*reg238; reg390=ponderation*reg328; reg779=reg83*reg779; reg731=reg83*reg731; reg641=reg83*reg641;
    reg583=reg83*reg583; reg622=reg83*reg622; reg735=reg83*reg735; reg362=reg83*reg362; reg728=reg83*reg728;
    reg394=reg83*reg394; reg393=ponderation*reg370; reg395=ponderation*reg376; reg267=reg83*reg267; reg608=reg83*reg608;
    reg738=reg83*reg738; reg399=reg83*reg399; reg400=ponderation*reg235; reg406=ponderation*reg289; reg594=reg83*reg594;
    reg445=reg83*reg445; reg510=reg83*reg510; reg597=reg83*reg597; reg385=reg83*reg385; reg460=reg83*reg460;
    reg552=reg83*reg552; reg268=reg83*reg268; reg452=reg83*reg452; reg325=reg83*reg325; reg760=reg83*reg760;
    reg745=reg83*reg745; reg245=reg83*reg245; reg408=ponderation*reg209; reg746=reg83*reg746; reg253=reg83*reg253;
    reg413=reg83*reg413; reg453=reg83*reg453; reg780=reg83*reg780; reg646=reg83*reg646; reg292=reg83*reg292;
    reg537=reg83*reg537; reg451=reg83*reg451; reg658=reg83*reg658; reg475=reg83*reg475; reg494=reg83*reg494;
    reg783=reg83*reg783; reg261=reg83*reg261; reg409=ponderation*reg277; reg518=reg83*reg518; reg685=reg83*reg685;
    reg278=reg83*reg278; reg762=reg83*reg762; reg593=reg83*reg593; reg599=reg83*reg599; reg748=reg83*reg748;
    reg553=reg83*reg553; reg580=reg83*reg580; reg757=reg83*reg757; T tmp_21_21=ponderation*reg733; T tmp_3_7=ponderation*reg821;
    T tmp_0_7=ponderation*reg728; T tmp_23_23=ponderation*reg414; T tmp_16_16=ponderation*reg255; T tmp_0_6=ponderation*reg445; T tmp_3_9=ponderation*reg808;
    T tmp_16_22=ponderation*reg334; T tmp_2_14=ponderation*reg331; T tmp_2_15=ponderation*reg696; T tmp_2_16=ponderation*reg697; T tmp_3_8=ponderation*reg798;
    T tmp_16_23=ponderation*reg700; T tmp_1_3=ponderation*reg745; T tmp_22_22=ponderation*reg423; T tmp_2_23=ponderation*reg319; T tmp_3_3=ponderation*reg789;
    T tmp_16_19=ponderation*reg377; T tmp_1_7=ponderation*reg725; T tmp_1_6=ponderation*reg748; T tmp_3_4=ponderation*reg803; T tmp_2_22=ponderation*reg829;
    T tmp_0_0=ponderation*reg420; T tmp_2_21=ponderation*reg705; T tmp_16_18=ponderation*reg818; T tmp_16_20=ponderation*reg825; T tmp_1_5=ponderation*reg746;
    T tmp_21_23=ponderation*reg754; T tmp_2_20=ponderation*reg305; T tmp_3_5=ponderation*reg805; T tmp_0_1=ponderation*reg751; T tmp_22_23=ponderation*reg749;
    T tmp_2_19=ponderation*reg703; T tmp_16_21=ponderation*reg359; T tmp_3_6=ponderation*reg826; T tmp_16_17=-reg301; T tmp_0_2=ponderation*reg753;
    T tmp_2_18=ponderation*reg702; T tmp_1_4=ponderation*reg780; T tmp_21_22=ponderation*reg729; T tmp_2_17=ponderation*reg323; T tmp_1_1=ponderation*reg452;
    T tmp_1_9=ponderation*reg686; T tmp_0_23=ponderation*reg757; T tmp_1_10=ponderation*reg247; T tmp_18_20=ponderation*reg660; T tmp_19_22=ponderation*reg217;
    T tmp_1_11=ponderation*reg689; T tmp_0_5=ponderation*reg783; T tmp_1_17=ponderation*reg692; T tmp_18_19=ponderation*reg663; T tmp_0_4=ponderation*reg781;
    T tmp_1_18=ponderation*reg662; T tmp_19_23=ponderation*reg784; T tmp_1_19=ponderation*reg356; T tmp_18_18=ponderation*reg669; T tmp_0_3=ponderation*reg457;
    T tmp_1_20=ponderation*reg666; T tmp_1_15=ponderation*reg762; T tmp_19_19=ponderation*reg281; T tmp_1_14=ponderation*reg760; T tmp_1_16=ponderation*reg461;
    T tmp_0_13=ponderation*reg770; T tmp_19_20=ponderation*reg763; T tmp_18_23=ponderation*reg676; T tmp_1_13=ponderation*reg451; T tmp_0_14=ponderation*reg767;
    T tmp_1_12=ponderation*reg691; T tmp_0_15=ponderation*reg270; T tmp_1_2=ponderation*reg453; T tmp_18_22=ponderation*reg684; T tmp_0_9=ponderation*reg679;
    T tmp_19_21=ponderation*reg450; T tmp_0_10=ponderation*reg683; T tmp_18_21=ponderation*reg688; T tmp_1_8=ponderation*reg685; T tmp_2_6=ponderation*reg711;
    T tmp_17_20=ponderation*reg276; T tmp_0_22=ponderation*reg738; T tmp_2_7=ponderation*reg715; T tmp_2_8=ponderation*reg318; T tmp_17_19=ponderation*reg719;
    T tmp_0_21=ponderation*reg394; T tmp_20_23=ponderation*reg387; T tmp_2_9=-reg371; T tmp_2_10=ponderation*reg718; T tmp_17_18=ponderation*reg724;
    T tmp_0_20=ponderation*reg735; T tmp_2_11=ponderation*reg298; T tmp_0_19=ponderation*reg779; T tmp_2_12=ponderation*reg723; T tmp_17_17=ponderation*reg330;
    T tmp_0_8=ponderation*reg731; T tmp_2_13=ponderation*reg694; T tmp_0_18=ponderation*reg405; T tmp_1_21=ponderation*reg668; T tmp_20_20=ponderation*reg396;
    T tmp_1_22=ponderation*reg346; T tmp_17_23=ponderation*reg294; T tmp_0_17=ponderation*reg774; T tmp_1_23=ponderation*reg671; T tmp_2_2=ponderation*reg374;
    T tmp_17_22=ponderation*reg708; T tmp_0_16=ponderation*reg677; T tmp_20_21=ponderation*reg775; T tmp_2_3=ponderation*reg673; T tmp_2_4=ponderation*reg707;
    T tmp_0_12=ponderation*reg385; T tmp_17_21=ponderation*reg712; T tmp_2_5=ponderation*reg364; T tmp_0_11=-reg406; T tmp_20_22=ponderation*reg768;
    T tmp_12_17=ponderation*reg645; T tmp_7_16=ponderation*reg245; T tmp_7_17=ponderation*reg646; T tmp_12_16=ponderation*reg613; T tmp_7_18=ponderation*reg658;
    T tmp_7_19=ponderation*reg292; T tmp_12_15=ponderation*reg300; T tmp_7_20=ponderation*reg253; T tmp_7_21=ponderation*reg325; T tmp_7_22=ponderation*reg268;
    T tmp_12_14=ponderation*reg511; T tmp_7_23=ponderation*reg553; T tmp_8_8=ponderation*reg278; T tmp_12_13=ponderation*reg274; T tmp_8_9=-reg409;
    T tmp_8_10=ponderation*reg494; T tmp_12_12=ponderation*reg456; T tmp_8_11=ponderation*reg332; T tmp_8_12=ponderation*reg545; T tmp_8_13=ponderation*reg476;
    T tmp_11_23=ponderation*reg227; T tmp_8_14=ponderation*reg454; T tmp_8_15=ponderation*reg508; T tmp_13_13=ponderation*reg367; T tmp_6_16=ponderation*reg573;
    T tmp_6_17=ponderation*reg577; T tmp_6_18=ponderation*reg363; T tmp_12_23=ponderation*reg561; T tmp_6_19=ponderation*reg585; T tmp_6_20=ponderation*reg560;
    T tmp_12_22=ponderation*reg568; T tmp_6_21=ponderation*reg373; T tmp_6_22=ponderation*reg567; T tmp_12_21=ponderation*reg349; T tmp_6_23=ponderation*reg571;
    T tmp_7_7=ponderation*reg345; T tmp_7_8=ponderation*reg602; T tmp_12_20=ponderation*reg609; T tmp_7_9=ponderation*reg605; T tmp_7_10=ponderation*reg355;
    T tmp_12_19=ponderation*reg592; T tmp_7_11=ponderation*reg587; T tmp_7_12=ponderation*reg590; T tmp_12_18=ponderation*reg258; T tmp_7_13=ponderation*reg261;
    T tmp_7_14=ponderation*reg593; T tmp_7_15=ponderation*reg599; T tmp_9_19=ponderation*reg267; T tmp_11_14=ponderation*reg439; T tmp_9_20=-reg400;
    T tmp_9_21=ponderation*reg510; T tmp_11_13=ponderation*reg470; T tmp_9_22=ponderation*reg502; T tmp_9_23=-reg379; T tmp_11_12=-reg266;
    T tmp_10_10=ponderation*reg443; T tmp_10_11=-reg381; T tmp_10_12=ponderation*reg487; T tmp_11_11=ponderation*reg279; T tmp_10_13=ponderation*reg421;
    T tmp_10_14=ponderation*reg518; T tmp_10_23=ponderation*reg535; T tmp_10_15=ponderation*reg552; T tmp_10_16=-reg408; T tmp_10_22=ponderation*reg412;
    T tmp_10_17=ponderation*reg537; T tmp_10_18=ponderation*reg475; T tmp_10_21=ponderation*reg581; T tmp_10_19=ponderation*reg413; T tmp_10_20=ponderation*reg580;
    T tmp_11_22=ponderation*reg655; T tmp_8_16=ponderation*reg492; T tmp_8_17=ponderation*reg448; T tmp_11_21=-reg271; T tmp_8_18=ponderation*reg578;
    T tmp_8_19=ponderation*reg565; T tmp_11_20=ponderation*reg398; T tmp_8_20=ponderation*reg460; T tmp_8_21=ponderation*reg597; T tmp_11_19=ponderation*reg591;
    T tmp_8_22=ponderation*reg594; T tmp_8_23=ponderation*reg399; T tmp_11_18=-reg341; T tmp_9_9=ponderation*reg608; T tmp_9_10=-reg393;
    T tmp_9_11=ponderation*reg362; T tmp_11_17=-reg293; T tmp_9_12=ponderation*reg583; T tmp_9_13=ponderation*reg641; T tmp_11_16=ponderation*reg617;
    T tmp_9_14=-reg390; T tmp_9_15=-reg386; T tmp_11_15=-reg347; T tmp_9_16=ponderation*reg622; T tmp_4_7=ponderation*reg273;
    T tmp_15_16=-reg343; T tmp_4_8=ponderation*reg521; T tmp_4_9=ponderation*reg522; T tmp_9_17=-reg395; T tmp_15_15=ponderation*reg286;
    T tmp_9_18=ponderation*reg528; T tmp_4_10=ponderation*reg243; T tmp_4_11=ponderation*reg513; T tmp_14_23=ponderation*reg248; T tmp_4_12=ponderation*reg517;
    T tmp_4_13=ponderation*reg249; T tmp_14_22=ponderation*reg550; T tmp_4_14=ponderation*reg547; T tmp_4_15=ponderation*reg549; T tmp_14_21=ponderation*reg557;
    T tmp_4_16=ponderation*reg240; T tmp_4_17=ponderation*reg556; T tmp_14_20=ponderation*reg263; T tmp_4_18=ponderation*reg530; T tmp_4_19=ponderation*reg264;
    T tmp_14_19=ponderation*reg544; T tmp_4_20=ponderation*reg539; T tmp_3_10=ponderation*reg794; T tmp_3_11=-reg368; T tmp_15_23=ponderation*reg797;
    T tmp_3_12=ponderation*reg828; T tmp_3_13=ponderation*reg280; T tmp_15_22=ponderation*reg788; T tmp_3_14=ponderation*reg786; T tmp_3_15=ponderation*reg242;
    T tmp_15_21=ponderation*reg295; T tmp_3_16=ponderation*reg816; T tmp_3_17=ponderation*reg830; T tmp_3_18=ponderation*reg785; T tmp_15_20=ponderation*reg809;
    T tmp_3_19=ponderation*reg799; T tmp_3_20=ponderation*reg806; T tmp_15_19=ponderation*reg824; T tmp_3_21=ponderation*reg304; T tmp_3_22=ponderation*reg795;
    T tmp_15_18=ponderation*reg275; T tmp_3_23=ponderation*reg171; T tmp_4_4=ponderation*reg285; T tmp_4_5=ponderation*reg465; T tmp_15_17=ponderation*reg250;
    T tmp_4_6=ponderation*reg422; T tmp_5_18=ponderation*reg495; T tmp_13_20=ponderation*reg435; T tmp_5_19=ponderation*reg626; T tmp_5_20=ponderation*reg340;
    T tmp_13_19=ponderation*reg324; T tmp_5_21=ponderation*reg629; T tmp_5_22=ponderation*reg632; T tmp_5_23=ponderation*reg327; T tmp_13_18=ponderation*reg616;
    T tmp_6_6=ponderation*reg329; T tmp_6_7=ponderation*reg615; T tmp_13_17=ponderation*reg650; T tmp_6_8=ponderation*reg620; T tmp_6_9=ponderation*reg649;
    T tmp_13_16=ponderation*reg314; T tmp_6_10=ponderation*reg653; T tmp_6_11=-reg366; T tmp_13_15=ponderation*reg638; T tmp_6_12=ponderation*reg311;
    T tmp_6_13=ponderation*reg637; T tmp_13_14=ponderation*reg380; T tmp_6_14=ponderation*reg642; T tmp_6_15=ponderation*reg322; T tmp_4_21=ponderation*reg543;
    T tmp_14_18=ponderation*reg479; T tmp_4_22=ponderation*reg260; T tmp_4_23=ponderation*reg478; T tmp_14_17=ponderation*reg352; T tmp_5_5=ponderation*reg357;
    T tmp_5_6=ponderation*reg485; T tmp_14_16=ponderation*reg469; T tmp_5_7=ponderation*reg467; T tmp_5_8=ponderation*reg351; T tmp_14_15=ponderation*reg498;
    T tmp_5_9=-reg365; T tmp_5_10=ponderation*reg497; T tmp_14_14=ponderation*reg372; T tmp_5_11=ponderation*reg237; T tmp_5_12=ponderation*reg505;
    T tmp_13_23=ponderation*reg486; T tmp_5_13=ponderation*reg509; T tmp_5_14=ponderation*reg321; T tmp_13_22=ponderation*reg309; T tmp_5_15=ponderation*reg490;
    T tmp_5_16=ponderation*reg491; T tmp_13_21=ponderation*reg496; T tmp_5_17=ponderation*reg308;
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    matrix(indices[0]+0,indices[0]+0) += tmp_0_0;
    matrix(indices[0]+0,indices[0]+1) += tmp_0_1;
    matrix(indices[0]+0,indices[0]+2) += tmp_0_2;
    matrix(indices[0]+0,indices[1]+0) += tmp_0_3;
    matrix(indices[0]+0,indices[1]+1) += tmp_0_4;
    matrix(indices[0]+0,indices[1]+2) += tmp_0_5;
    matrix(indices[0]+0,indices[2]+0) += tmp_0_6;
    matrix(indices[0]+0,indices[2]+1) += tmp_0_7;
    matrix(indices[0]+0,indices[2]+2) += tmp_0_8;
    matrix(indices[0]+0,indices[3]+0) += tmp_0_9;
    matrix(indices[0]+0,indices[3]+1) += tmp_0_10;
    matrix(indices[0]+0,indices[3]+2) += tmp_0_11;
    matrix(indices[0]+0,indices[4]+0) += tmp_0_12;
    matrix(indices[0]+0,indices[4]+1) += tmp_0_13;
    matrix(indices[0]+0,indices[4]+2) += tmp_0_14;
    matrix(indices[0]+0,indices[5]+0) += tmp_0_15;
    matrix(indices[0]+0,indices[5]+1) += tmp_0_16;
    matrix(indices[0]+0,indices[5]+2) += tmp_0_17;
    matrix(indices[0]+0,indices[6]+0) += tmp_0_18;
    matrix(indices[0]+0,indices[6]+1) += tmp_0_19;
    matrix(indices[0]+0,indices[6]+2) += tmp_0_20;
    matrix(indices[0]+0,indices[7]+0) += tmp_0_21;
    matrix(indices[0]+0,indices[7]+1) += tmp_0_22;
    matrix(indices[0]+0,indices[7]+2) += tmp_0_23;
    matrix(indices[0]+1,indices[0]+1) += tmp_1_1;
    matrix(indices[0]+1,indices[0]+2) += tmp_1_2;
    matrix(indices[0]+1,indices[1]+0) += tmp_1_3;
    matrix(indices[0]+1,indices[1]+1) += tmp_1_4;
    matrix(indices[0]+1,indices[1]+2) += tmp_1_5;
    matrix(indices[0]+1,indices[2]+0) += tmp_1_6;
    matrix(indices[0]+1,indices[2]+1) += tmp_1_7;
    matrix(indices[0]+1,indices[2]+2) += tmp_1_8;
    matrix(indices[0]+1,indices[3]+0) += tmp_1_9;
    matrix(indices[0]+1,indices[3]+1) += tmp_1_10;
    matrix(indices[0]+1,indices[3]+2) += tmp_1_11;
    matrix(indices[0]+1,indices[4]+0) += tmp_1_12;
    matrix(indices[0]+1,indices[4]+1) += tmp_1_13;
    matrix(indices[0]+1,indices[4]+2) += tmp_1_14;
    matrix(indices[0]+1,indices[5]+0) += tmp_1_15;
    matrix(indices[0]+1,indices[5]+1) += tmp_1_16;
    matrix(indices[0]+1,indices[5]+2) += tmp_1_17;
    matrix(indices[0]+1,indices[6]+0) += tmp_1_18;
    matrix(indices[0]+1,indices[6]+1) += tmp_1_19;
    matrix(indices[0]+1,indices[6]+2) += tmp_1_20;
    matrix(indices[0]+1,indices[7]+0) += tmp_1_21;
    matrix(indices[0]+1,indices[7]+1) += tmp_1_22;
    matrix(indices[0]+1,indices[7]+2) += tmp_1_23;
    matrix(indices[0]+2,indices[0]+2) += tmp_2_2;
    matrix(indices[0]+2,indices[1]+0) += tmp_2_3;
    matrix(indices[0]+2,indices[1]+1) += tmp_2_4;
    matrix(indices[0]+2,indices[1]+2) += tmp_2_5;
    matrix(indices[0]+2,indices[2]+0) += tmp_2_6;
    matrix(indices[0]+2,indices[2]+1) += tmp_2_7;
    matrix(indices[0]+2,indices[2]+2) += tmp_2_8;
    matrix(indices[0]+2,indices[3]+0) += tmp_2_9;
    matrix(indices[0]+2,indices[3]+1) += tmp_2_10;
    matrix(indices[0]+2,indices[3]+2) += tmp_2_11;
    matrix(indices[0]+2,indices[4]+0) += tmp_2_12;
    matrix(indices[0]+2,indices[4]+1) += tmp_2_13;
    matrix(indices[0]+2,indices[4]+2) += tmp_2_14;
    matrix(indices[0]+2,indices[5]+0) += tmp_2_15;
    matrix(indices[0]+2,indices[5]+1) += tmp_2_16;
    matrix(indices[0]+2,indices[5]+2) += tmp_2_17;
    matrix(indices[0]+2,indices[6]+0) += tmp_2_18;
    matrix(indices[0]+2,indices[6]+1) += tmp_2_19;
    matrix(indices[0]+2,indices[6]+2) += tmp_2_20;
    matrix(indices[0]+2,indices[7]+0) += tmp_2_21;
    matrix(indices[0]+2,indices[7]+1) += tmp_2_22;
    matrix(indices[0]+2,indices[7]+2) += tmp_2_23;
    matrix(indices[1]+0,indices[1]+0) += tmp_3_3;
    matrix(indices[1]+0,indices[1]+1) += tmp_3_4;
    matrix(indices[1]+0,indices[1]+2) += tmp_3_5;
    matrix(indices[1]+0,indices[2]+0) += tmp_3_6;
    matrix(indices[1]+0,indices[2]+1) += tmp_3_7;
    matrix(indices[1]+0,indices[2]+2) += tmp_3_8;
    matrix(indices[1]+0,indices[3]+0) += tmp_3_9;
    matrix(indices[1]+0,indices[3]+1) += tmp_3_10;
    matrix(indices[1]+0,indices[3]+2) += tmp_3_11;
    matrix(indices[1]+0,indices[4]+0) += tmp_3_12;
    matrix(indices[1]+0,indices[4]+1) += tmp_3_13;
    matrix(indices[1]+0,indices[4]+2) += tmp_3_14;
    matrix(indices[1]+0,indices[5]+0) += tmp_3_15;
    matrix(indices[1]+0,indices[5]+1) += tmp_3_16;
    matrix(indices[1]+0,indices[5]+2) += tmp_3_17;
    matrix(indices[1]+0,indices[6]+0) += tmp_3_18;
    matrix(indices[1]+0,indices[6]+1) += tmp_3_19;
    matrix(indices[1]+0,indices[6]+2) += tmp_3_20;
    matrix(indices[1]+0,indices[7]+0) += tmp_3_21;
    matrix(indices[1]+0,indices[7]+1) += tmp_3_22;
    matrix(indices[1]+0,indices[7]+2) += tmp_3_23;
    matrix(indices[1]+1,indices[1]+1) += tmp_4_4;
    matrix(indices[1]+1,indices[1]+2) += tmp_4_5;
    matrix(indices[1]+1,indices[2]+0) += tmp_4_6;
    matrix(indices[1]+1,indices[2]+1) += tmp_4_7;
    matrix(indices[1]+1,indices[2]+2) += tmp_4_8;
    matrix(indices[1]+1,indices[3]+0) += tmp_4_9;
    matrix(indices[1]+1,indices[3]+1) += tmp_4_10;
    matrix(indices[1]+1,indices[3]+2) += tmp_4_11;
    matrix(indices[1]+1,indices[4]+0) += tmp_4_12;
    matrix(indices[1]+1,indices[4]+1) += tmp_4_13;
    matrix(indices[1]+1,indices[4]+2) += tmp_4_14;
    matrix(indices[1]+1,indices[5]+0) += tmp_4_15;
    matrix(indices[1]+1,indices[5]+1) += tmp_4_16;
    matrix(indices[1]+1,indices[5]+2) += tmp_4_17;
    matrix(indices[1]+1,indices[6]+0) += tmp_4_18;
    matrix(indices[1]+1,indices[6]+1) += tmp_4_19;
    matrix(indices[1]+1,indices[6]+2) += tmp_4_20;
    matrix(indices[1]+1,indices[7]+0) += tmp_4_21;
    matrix(indices[1]+1,indices[7]+1) += tmp_4_22;
    matrix(indices[1]+1,indices[7]+2) += tmp_4_23;
    matrix(indices[1]+2,indices[1]+2) += tmp_5_5;
    matrix(indices[1]+2,indices[2]+0) += tmp_5_6;
    matrix(indices[1]+2,indices[2]+1) += tmp_5_7;
    matrix(indices[1]+2,indices[2]+2) += tmp_5_8;
    matrix(indices[1]+2,indices[3]+0) += tmp_5_9;
    matrix(indices[1]+2,indices[3]+1) += tmp_5_10;
    matrix(indices[1]+2,indices[3]+2) += tmp_5_11;
    matrix(indices[1]+2,indices[4]+0) += tmp_5_12;
    matrix(indices[1]+2,indices[4]+1) += tmp_5_13;
    matrix(indices[1]+2,indices[4]+2) += tmp_5_14;
    matrix(indices[1]+2,indices[5]+0) += tmp_5_15;
    matrix(indices[1]+2,indices[5]+1) += tmp_5_16;
    matrix(indices[1]+2,indices[5]+2) += tmp_5_17;
    matrix(indices[1]+2,indices[6]+0) += tmp_5_18;
    matrix(indices[1]+2,indices[6]+1) += tmp_5_19;
    matrix(indices[1]+2,indices[6]+2) += tmp_5_20;
    matrix(indices[1]+2,indices[7]+0) += tmp_5_21;
    matrix(indices[1]+2,indices[7]+1) += tmp_5_22;
    matrix(indices[1]+2,indices[7]+2) += tmp_5_23;
    matrix(indices[2]+0,indices[2]+0) += tmp_6_6;
    matrix(indices[2]+0,indices[2]+1) += tmp_6_7;
    matrix(indices[2]+0,indices[2]+2) += tmp_6_8;
    matrix(indices[2]+0,indices[3]+0) += tmp_6_9;
    matrix(indices[2]+0,indices[3]+1) += tmp_6_10;
    matrix(indices[2]+0,indices[3]+2) += tmp_6_11;
    matrix(indices[2]+0,indices[4]+0) += tmp_6_12;
    matrix(indices[2]+0,indices[4]+1) += tmp_6_13;
    matrix(indices[2]+0,indices[4]+2) += tmp_6_14;
    matrix(indices[2]+0,indices[5]+0) += tmp_6_15;
    matrix(indices[2]+0,indices[5]+1) += tmp_6_16;
    matrix(indices[2]+0,indices[5]+2) += tmp_6_17;
    matrix(indices[2]+0,indices[6]+0) += tmp_6_18;
    matrix(indices[2]+0,indices[6]+1) += tmp_6_19;
    matrix(indices[2]+0,indices[6]+2) += tmp_6_20;
    matrix(indices[2]+0,indices[7]+0) += tmp_6_21;
    matrix(indices[2]+0,indices[7]+1) += tmp_6_22;
    matrix(indices[2]+0,indices[7]+2) += tmp_6_23;
    matrix(indices[2]+1,indices[2]+1) += tmp_7_7;
    matrix(indices[2]+1,indices[2]+2) += tmp_7_8;
    matrix(indices[2]+1,indices[3]+0) += tmp_7_9;
    matrix(indices[2]+1,indices[3]+1) += tmp_7_10;
    matrix(indices[2]+1,indices[3]+2) += tmp_7_11;
    matrix(indices[2]+1,indices[4]+0) += tmp_7_12;
    matrix(indices[2]+1,indices[4]+1) += tmp_7_13;
    matrix(indices[2]+1,indices[4]+2) += tmp_7_14;
    matrix(indices[2]+1,indices[5]+0) += tmp_7_15;
    matrix(indices[2]+1,indices[5]+1) += tmp_7_16;
    matrix(indices[2]+1,indices[5]+2) += tmp_7_17;
    matrix(indices[2]+1,indices[6]+0) += tmp_7_18;
    matrix(indices[2]+1,indices[6]+1) += tmp_7_19;
    matrix(indices[2]+1,indices[6]+2) += tmp_7_20;
    matrix(indices[2]+1,indices[7]+0) += tmp_7_21;
    matrix(indices[2]+1,indices[7]+1) += tmp_7_22;
    matrix(indices[2]+1,indices[7]+2) += tmp_7_23;
    matrix(indices[2]+2,indices[2]+2) += tmp_8_8;
    matrix(indices[2]+2,indices[3]+0) += tmp_8_9;
    matrix(indices[2]+2,indices[3]+1) += tmp_8_10;
    matrix(indices[2]+2,indices[3]+2) += tmp_8_11;
    matrix(indices[2]+2,indices[4]+0) += tmp_8_12;
    matrix(indices[2]+2,indices[4]+1) += tmp_8_13;
    matrix(indices[2]+2,indices[4]+2) += tmp_8_14;
    matrix(indices[2]+2,indices[5]+0) += tmp_8_15;
    matrix(indices[2]+2,indices[5]+1) += tmp_8_16;
    matrix(indices[2]+2,indices[5]+2) += tmp_8_17;
    matrix(indices[2]+2,indices[6]+0) += tmp_8_18;
    matrix(indices[2]+2,indices[6]+1) += tmp_8_19;
    matrix(indices[2]+2,indices[6]+2) += tmp_8_20;
    matrix(indices[2]+2,indices[7]+0) += tmp_8_21;
    matrix(indices[2]+2,indices[7]+1) += tmp_8_22;
    matrix(indices[2]+2,indices[7]+2) += tmp_8_23;
    matrix(indices[3]+0,indices[3]+0) += tmp_9_9;
    matrix(indices[3]+0,indices[3]+1) += tmp_9_10;
    matrix(indices[3]+0,indices[3]+2) += tmp_9_11;
    matrix(indices[3]+0,indices[4]+0) += tmp_9_12;
    matrix(indices[3]+0,indices[4]+1) += tmp_9_13;
    matrix(indices[3]+0,indices[4]+2) += tmp_9_14;
    matrix(indices[3]+0,indices[5]+0) += tmp_9_15;
    matrix(indices[3]+0,indices[5]+1) += tmp_9_16;
    matrix(indices[3]+0,indices[5]+2) += tmp_9_17;
    matrix(indices[3]+0,indices[6]+0) += tmp_9_18;
    matrix(indices[3]+0,indices[6]+1) += tmp_9_19;
    matrix(indices[3]+0,indices[6]+2) += tmp_9_20;
    matrix(indices[3]+0,indices[7]+0) += tmp_9_21;
    matrix(indices[3]+0,indices[7]+1) += tmp_9_22;
    matrix(indices[3]+0,indices[7]+2) += tmp_9_23;
    matrix(indices[3]+1,indices[3]+1) += tmp_10_10;
    matrix(indices[3]+1,indices[3]+2) += tmp_10_11;
    matrix(indices[3]+1,indices[4]+0) += tmp_10_12;
    matrix(indices[3]+1,indices[4]+1) += tmp_10_13;
    matrix(indices[3]+1,indices[4]+2) += tmp_10_14;
    matrix(indices[3]+1,indices[5]+0) += tmp_10_15;
    matrix(indices[3]+1,indices[5]+1) += tmp_10_16;
    matrix(indices[3]+1,indices[5]+2) += tmp_10_17;
    matrix(indices[3]+1,indices[6]+0) += tmp_10_18;
    matrix(indices[3]+1,indices[6]+1) += tmp_10_19;
    matrix(indices[3]+1,indices[6]+2) += tmp_10_20;
    matrix(indices[3]+1,indices[7]+0) += tmp_10_21;
    matrix(indices[3]+1,indices[7]+1) += tmp_10_22;
    matrix(indices[3]+1,indices[7]+2) += tmp_10_23;
    matrix(indices[3]+2,indices[3]+2) += tmp_11_11;
    matrix(indices[3]+2,indices[4]+0) += tmp_11_12;
    matrix(indices[3]+2,indices[4]+1) += tmp_11_13;
    matrix(indices[3]+2,indices[4]+2) += tmp_11_14;
    matrix(indices[3]+2,indices[5]+0) += tmp_11_15;
    matrix(indices[3]+2,indices[5]+1) += tmp_11_16;
    matrix(indices[3]+2,indices[5]+2) += tmp_11_17;
    matrix(indices[3]+2,indices[6]+0) += tmp_11_18;
    matrix(indices[3]+2,indices[6]+1) += tmp_11_19;
    matrix(indices[3]+2,indices[6]+2) += tmp_11_20;
    matrix(indices[3]+2,indices[7]+0) += tmp_11_21;
    matrix(indices[3]+2,indices[7]+1) += tmp_11_22;
    matrix(indices[3]+2,indices[7]+2) += tmp_11_23;
    matrix(indices[4]+0,indices[4]+0) += tmp_12_12;
    matrix(indices[4]+0,indices[4]+1) += tmp_12_13;
    matrix(indices[4]+0,indices[4]+2) += tmp_12_14;
    matrix(indices[4]+0,indices[5]+0) += tmp_12_15;
    matrix(indices[4]+0,indices[5]+1) += tmp_12_16;
    matrix(indices[4]+0,indices[5]+2) += tmp_12_17;
    matrix(indices[4]+0,indices[6]+0) += tmp_12_18;
    matrix(indices[4]+0,indices[6]+1) += tmp_12_19;
    matrix(indices[4]+0,indices[6]+2) += tmp_12_20;
    matrix(indices[4]+0,indices[7]+0) += tmp_12_21;
    matrix(indices[4]+0,indices[7]+1) += tmp_12_22;
    matrix(indices[4]+0,indices[7]+2) += tmp_12_23;
    matrix(indices[4]+1,indices[4]+1) += tmp_13_13;
    matrix(indices[4]+1,indices[4]+2) += tmp_13_14;
    matrix(indices[4]+1,indices[5]+0) += tmp_13_15;
    matrix(indices[4]+1,indices[5]+1) += tmp_13_16;
    matrix(indices[4]+1,indices[5]+2) += tmp_13_17;
    matrix(indices[4]+1,indices[6]+0) += tmp_13_18;
    matrix(indices[4]+1,indices[6]+1) += tmp_13_19;
    matrix(indices[4]+1,indices[6]+2) += tmp_13_20;
    matrix(indices[4]+1,indices[7]+0) += tmp_13_21;
    matrix(indices[4]+1,indices[7]+1) += tmp_13_22;
    matrix(indices[4]+1,indices[7]+2) += tmp_13_23;
    matrix(indices[4]+2,indices[4]+2) += tmp_14_14;
    matrix(indices[4]+2,indices[5]+0) += tmp_14_15;
    matrix(indices[4]+2,indices[5]+1) += tmp_14_16;
    matrix(indices[4]+2,indices[5]+2) += tmp_14_17;
    matrix(indices[4]+2,indices[6]+0) += tmp_14_18;
    matrix(indices[4]+2,indices[6]+1) += tmp_14_19;
    matrix(indices[4]+2,indices[6]+2) += tmp_14_20;
    matrix(indices[4]+2,indices[7]+0) += tmp_14_21;
    matrix(indices[4]+2,indices[7]+1) += tmp_14_22;
    matrix(indices[4]+2,indices[7]+2) += tmp_14_23;
    matrix(indices[5]+0,indices[5]+0) += tmp_15_15;
    matrix(indices[5]+0,indices[5]+1) += tmp_15_16;
    matrix(indices[5]+0,indices[5]+2) += tmp_15_17;
    matrix(indices[5]+0,indices[6]+0) += tmp_15_18;
    matrix(indices[5]+0,indices[6]+1) += tmp_15_19;
    matrix(indices[5]+0,indices[6]+2) += tmp_15_20;
    matrix(indices[5]+0,indices[7]+0) += tmp_15_21;
    matrix(indices[5]+0,indices[7]+1) += tmp_15_22;
    matrix(indices[5]+0,indices[7]+2) += tmp_15_23;
    matrix(indices[5]+1,indices[5]+1) += tmp_16_16;
    matrix(indices[5]+1,indices[5]+2) += tmp_16_17;
    matrix(indices[5]+1,indices[6]+0) += tmp_16_18;
    matrix(indices[5]+1,indices[6]+1) += tmp_16_19;
    matrix(indices[5]+1,indices[6]+2) += tmp_16_20;
    matrix(indices[5]+1,indices[7]+0) += tmp_16_21;
    matrix(indices[5]+1,indices[7]+1) += tmp_16_22;
    matrix(indices[5]+1,indices[7]+2) += tmp_16_23;
    matrix(indices[5]+2,indices[5]+2) += tmp_17_17;
    matrix(indices[5]+2,indices[6]+0) += tmp_17_18;
    matrix(indices[5]+2,indices[6]+1) += tmp_17_19;
    matrix(indices[5]+2,indices[6]+2) += tmp_17_20;
    matrix(indices[5]+2,indices[7]+0) += tmp_17_21;
    matrix(indices[5]+2,indices[7]+1) += tmp_17_22;
    matrix(indices[5]+2,indices[7]+2) += tmp_17_23;
    matrix(indices[6]+0,indices[6]+0) += tmp_18_18;
    matrix(indices[6]+0,indices[6]+1) += tmp_18_19;
    matrix(indices[6]+0,indices[6]+2) += tmp_18_20;
    matrix(indices[6]+0,indices[7]+0) += tmp_18_21;
    matrix(indices[6]+0,indices[7]+1) += tmp_18_22;
    matrix(indices[6]+0,indices[7]+2) += tmp_18_23;
    matrix(indices[6]+1,indices[6]+1) += tmp_19_19;
    matrix(indices[6]+1,indices[6]+2) += tmp_19_20;
    matrix(indices[6]+1,indices[7]+0) += tmp_19_21;
    matrix(indices[6]+1,indices[7]+1) += tmp_19_22;
    matrix(indices[6]+1,indices[7]+2) += tmp_19_23;
    matrix(indices[6]+2,indices[6]+2) += tmp_20_20;
    matrix(indices[6]+2,indices[7]+0) += tmp_20_21;
    matrix(indices[6]+2,indices[7]+1) += tmp_20_22;
    matrix(indices[6]+2,indices[7]+2) += tmp_20_23;
    matrix(indices[7]+0,indices[7]+0) += tmp_21_21;
    matrix(indices[7]+0,indices[7]+1) += tmp_21_22;
    matrix(indices[7]+0,indices[7]+2) += tmp_21_23;
    matrix(indices[7]+1,indices[7]+1) += tmp_22_22;
    matrix(indices[7]+1,indices[7]+2) += tmp_22_23;
    matrix(indices[7]+2,indices[7]+2) += tmp_23_23;
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim,unsigned symmetric_version>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1-var_inter[2]; T reg1=1-var_inter[1]; T reg2=1-var_inter[0]; T reg3=reg1*reg2; T reg4=reg1*var_inter[0];
    T reg5=reg0*reg2; T reg6=reg0*var_inter[0]; T reg7=reg1*reg0; T reg8=elem.pos(1)[1]*reg4; T reg9=elem.pos(0)[1]*reg3;
    T reg10=elem.pos(0)[2]*reg5; T reg11=elem.pos(1)[2]*reg6; T reg12=reg0*var_inter[1]; T reg13=var_inter[0]*var_inter[1]; T reg14=elem.pos(0)[1]*reg5;
    T reg15=elem.pos(1)[1]*reg6; T reg16=elem.pos(0)[2]*reg7; T reg17=elem.pos(1)[2]*reg7; T reg18=elem.pos(1)[1]*reg7; T reg19=elem.pos(0)[1]*reg7;
    T reg20=elem.pos(0)[2]*reg3; T reg21=reg4*elem.pos(1)[2]; T reg22=reg13*elem.pos(2)[2]; T reg23=reg21+reg20; T reg24=reg11+reg10;
    T reg25=elem.pos(2)[2]*reg6; T reg26=var_inter[1]*reg2; T reg27=elem.pos(2)[1]*reg6; T reg28=reg15+reg14; T reg29=elem.pos(2)[1]*reg12;
    reg18=reg18-reg19; T reg30=reg8+reg9; T reg31=reg13*elem.pos(2)[1]; T reg32=elem.pos(2)[2]*reg12; reg17=reg17-reg16;
    T reg33=elem.pos(3)[2]*reg12; T reg34=elem.pos(1)[0]*reg7; T reg35=elem.pos(0)[0]*reg7; T reg36=var_inter[2]*reg2; T reg37=elem.pos(3)[1]*reg12;
    T reg38=elem.pos(0)[0]*reg5; T reg39=reg1*var_inter[2]; T reg40=elem.pos(1)[0]*reg6; reg29=reg18+reg29; reg32=reg17+reg32;
    reg27=reg27-reg28; reg17=elem.pos(3)[1]*reg5; reg18=reg23+reg22; T reg41=elem.pos(3)[2]*reg5; T reg42=reg26*elem.pos(3)[2];
    reg25=reg25-reg24; T reg43=reg26*elem.pos(3)[1]; T reg44=reg30+reg31; T reg45=elem.pos(4)[2]*reg3; T reg46=reg42+reg18;
    T reg47=elem.pos(4)[1]*reg3; T reg48=reg44+reg43; reg32=reg32-reg33; T reg49=elem.pos(4)[2]*reg39; T reg50=elem.pos(0)[0]*reg3;
    T reg51=elem.pos(1)[0]*reg4; T reg52=elem.pos(4)[2]*reg36; reg25=reg41+reg25; reg41=reg40+reg38; T reg53=elem.pos(2)[0]*reg6;
    T reg54=elem.pos(4)[1]*reg36; reg17=reg27+reg17; reg27=var_inter[2]*var_inter[0]; T reg55=1+(*f.m).poisson_ratio; T reg56=elem.pos(4)[1]*reg39;
    reg29=reg29-reg37; T reg57=elem.pos(2)[0]*reg12; reg34=reg34-reg35; T reg58=reg13*elem.pos(2)[0]; T reg59=elem.pos(3)[0]*reg5;
    T reg60=reg51+reg50; reg57=reg34+reg57; reg53=reg53-reg41; reg34=elem.pos(3)[0]*reg12; T reg61=elem.pos(5)[2]*reg27;
    reg25=reg25-reg52; T reg62=elem.pos(5)[1]*reg27; reg55=reg55/(*f.m).elastic_modulus; reg17=reg17-reg54; reg29=reg29-reg56;
    T reg63=elem.pos(5)[1]*reg39; T reg64=reg4*elem.pos(5)[2]; reg45=reg45-reg46; T reg65=reg4*elem.pos(5)[1]; reg47=reg47-reg48;
    T reg66=var_inter[2]*var_inter[1]; reg32=reg32-reg49; T reg67=elem.pos(5)[2]*reg39; T reg68=reg60+reg58; T reg69=reg13*elem.pos(6)[2];
    reg64=reg45+reg64; reg17=reg17-reg62; reg45=elem.pos(6)[1]*reg27; T reg70=reg13*elem.pos(6)[1]; reg65=reg47+reg65;
    reg25=reg25-reg61; reg47=elem.pos(6)[2]*reg27; T reg71=reg26*elem.pos(3)[0]; reg67=reg32+reg67; reg32=elem.pos(6)[2]*reg66;
    T reg72=elem.pos(6)[1]*reg66; reg63=reg29+reg63; reg29=elem.pos(4)[0]*reg39; reg57=reg57-reg34; reg53=reg59+reg53;
    reg59=elem.pos(4)[0]*reg36; T reg73=pow(reg55,2); T reg74=reg26*elem.pos(7)[2]; reg69=reg64+reg69; reg64=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg75=elem.pos(4)[0]*reg3; T reg76=reg71+reg68; reg72=reg63+reg72; reg63=elem.pos(5)[0]*reg39; reg57=reg57-reg29;
    T reg77=1.0/(*f.m).elastic_modulus; T reg78=elem.pos(7)[2]*reg66; reg45=reg17+reg45; reg32=reg67+reg32; reg17=elem.pos(7)[2]*reg36;
    reg47=reg25+reg47; reg25=elem.pos(7)[1]*reg36; reg53=reg53-reg59; reg67=elem.pos(5)[0]*reg27; reg55=reg55*reg73;
    T reg79=elem.pos(7)[1]*reg66; T reg80=reg26*elem.pos(7)[1]; reg70=reg65+reg70; reg72=reg72-reg79; reg63=reg57+reg63;
    reg80=reg70+reg80; reg57=elem.pos(6)[0]*reg66; reg65=reg4*elem.pos(5)[0]; reg75=reg75-reg76; reg32=reg32-reg78;
    reg17=reg47+reg17; reg53=reg53-reg67; reg47=elem.pos(6)[0]*reg27; reg74=reg69+reg74; reg69=reg64*reg55;
    reg55=reg77*reg55; reg25=reg45+reg25; reg45=reg32*reg80; reg70=reg17*reg80; T reg81=reg25*reg74;
    T reg82=reg72*reg74; T reg83=elem.pos(7)[0]*reg66; reg57=reg63+reg57; reg63=reg64*reg55; reg55=reg77*reg55;
    T reg84=elem.pos(7)[0]*reg36; reg47=reg53+reg47; reg53=reg64*reg69; reg65=reg75+reg65; reg75=reg13*elem.pos(6)[0];
    T reg85=reg32*reg25; T reg86=reg72*reg17; reg45=reg82-reg45; reg70=reg81-reg70; reg84=reg47+reg84;
    reg75=reg65+reg75; reg57=reg57-reg83; reg55=reg55-reg53; reg63=reg53+reg63; reg69=reg77*reg69;
    reg47=reg26*elem.pos(7)[0]; reg85=reg86-reg85; reg47=reg75+reg47; reg65=reg64*reg63; reg75=reg84*reg45;
    reg81=reg57*reg70; reg82=reg77*reg55; reg69=reg53+reg69; reg53=reg32*reg47; reg86=reg57*reg80;
    T reg87=reg72*reg47; T reg88=reg25*reg47; T reg89=reg57*reg74; reg80=reg84*reg80; T reg90=reg17*reg47;
    reg74=reg84*reg74; reg47=reg47*reg85; reg75=reg81-reg75; reg81=reg64*reg69; reg65=reg82-reg65;
    reg32=reg32*reg84; reg25=reg57*reg25; reg17=reg57*reg17; reg81=reg65-reg81; reg87=reg86-reg87;
    reg53=reg89-reg53; reg84=reg72*reg84; reg88=reg80-reg88; reg90=reg74-reg90; reg47=reg75+reg47;
    reg69=reg69/reg81; reg55=reg55/reg81; reg63=reg63/reg81; reg57=(*f.m).deltaT*(*f.m).alpha; reg84=reg25-reg84;
    reg53=reg53/reg47; reg45=reg45/reg47; reg87=reg87/reg47; reg88=reg88/reg47; reg32=reg17-reg32;
    reg90=reg90/reg47; reg70=reg70/reg47; reg17=reg5*reg87; reg25=reg27*reg53; reg84=reg84/reg47;
    reg65=reg12*reg88; reg85=reg85/reg47; reg72=reg5*reg45; reg74=reg55*reg57; reg75=reg63*reg57;
    reg80=reg69*reg57; reg82=reg12*reg70; reg86=reg39*reg90; reg32=reg32/reg47; reg89=reg26*reg84;
    T reg91=reg36*reg87; T reg92=reg66*reg90; T reg93=reg26*reg85; T reg94=reg72+reg82; T reg95=reg66*reg70;
    T reg96=reg12*reg90; T reg97=reg25+reg86; T reg98=reg7*reg88; T reg99=reg27*reg45; T reg100=reg7*reg70;
    T reg101=reg27*reg87; T reg102=reg39*reg88; T reg103=reg36*reg53; T reg104=reg36*reg45; T reg105=reg39*reg70;
    T reg106=reg75+reg80; T reg107=reg7*reg90; T reg108=reg6*reg87; T reg109=reg74+reg75; T reg110=reg6*reg45;
    T reg111=reg4*reg32; T reg112=reg6*reg53; T reg113=reg17+reg65; T reg114=reg5*reg53; T reg115=reg66*reg88;
    T reg116=reg74+reg106; T reg117=reg95-reg99; T reg118=reg109+reg80; T reg119=reg104-reg105; T reg120=reg115-reg101;
    T reg121=reg107-reg114; T reg122=reg25-reg92; T reg123=reg3*reg32; T reg124=reg4*reg85; T reg125=reg100+reg110;
    T reg126=reg13*reg32; T reg127=reg72-reg100; T reg128=reg3*reg85; T reg129=reg3*reg84; T reg130=reg82-reg110;
    T reg131=reg107+reg112; T reg132=reg13*reg85; T reg133=reg89+reg113; T reg134=reg17-reg98; T reg135=reg13*reg84;
    T reg136=reg112-reg96; T reg137=reg103+reg92; T reg138=reg91-reg102; T reg139=reg104+reg95; T reg140=reg91+reg115;
    T reg141=reg97+reg111; T reg142=reg65-reg108; T reg143=reg94+reg93; T reg144=reg4*var_inter[2]; T reg145=reg99+reg105;
    T reg146=reg4*reg84; T reg147=reg108+reg98; T reg148=reg101+reg102; T reg149=reg26*reg32; T reg150=reg114+reg96;
    T reg151=reg86-reg103; T reg152=reg26*reg0; reg130=reg130-reg132; reg137=reg137-reg149; T reg153=reg93-reg139;
    reg136=reg136+reg126; T reg154=reg89-reg140; T reg155=reg152*(*f.m).f_vol[2]; reg120=reg120+reg135; T reg156=var_inter[2]*reg3;
    T reg157=reg152*(*f.m).f_vol[0]; T reg158=var_inter[2]*reg13; T reg159=reg26*var_inter[2]; T reg160=reg4*reg0; T reg161=reg111-reg131;
    T reg162=reg3*reg0; T reg163=reg0*reg13; reg125=reg125-reg124; reg147=reg147-reg146; T reg164=reg143*reg118;
    T reg165=reg133*reg116; T reg166=reg141*reg118; reg119=reg119+reg128; reg134=reg134-reg129; reg127=reg127-reg128;
    reg117=reg117+reg132; reg122=reg122-reg126; reg150=reg150+reg149; reg148=reg146+reg148; T reg167=reg144*(*f.m).f_vol[1];
    reg151=reg151-reg123; reg138=reg138+reg129; reg145=reg145+reg124; reg121=reg121+reg123; reg142=reg142-reg135;
    T reg168=reg166-reg167; T reg169=reg119*reg118; T reg170=reg125*reg118; T reg171=reg137*reg118; T reg172=reg134*reg116;
    T reg173=(*f.m).f_vol[0]*reg162; T reg174=reg151*reg118; T reg175=reg154*reg116; T reg176=reg161*reg118; T reg177=reg153*reg118;
    T reg178=reg165-reg155; T reg179=reg147*reg116; T reg180=reg144*(*f.m).f_vol[0]; T reg181=reg148*reg116; T reg182=reg120*reg116;
    T reg183=reg130*reg118; T reg184=reg150*reg118; T reg185=reg136*reg118; T reg186=reg122*reg118; T reg187=reg164-reg157;
    T reg188=reg142*reg116; T reg189=reg117*reg118; T reg190=reg159*(*f.m).f_vol[0]; T reg191=reg163*(*f.m).f_vol[2]; T reg192=reg158*(*f.m).f_vol[0];
    T reg193=reg152*(*f.m).f_vol[1]; T reg194=reg144*(*f.m).f_vol[2]; T reg195=reg158*(*f.m).f_vol[2]; T reg196=reg158*(*f.m).f_vol[1]; T reg197=reg156*(*f.m).f_vol[0];
    T reg198=reg159*(*f.m).f_vol[1]; T reg199=reg159*(*f.m).f_vol[2]; T reg200=reg160*(*f.m).f_vol[2]; T reg201=reg160*(*f.m).f_vol[1]; T reg202=reg160*(*f.m).f_vol[0];
    T reg203=reg162*(*f.m).f_vol[2]; T reg204=reg156*(*f.m).f_vol[1]; T reg205=reg121*reg118; T reg206=reg156*(*f.m).f_vol[2]; T reg207=(*f.m).f_vol[0]*reg163;
    T reg208=reg163*(*f.m).f_vol[1]; T reg209=reg145*reg118; T reg210=reg138*reg116; T reg211=reg162*(*f.m).f_vol[1]; T reg212=reg127*reg118;
    T reg213=reg204+reg174; T reg214=reg180+reg209; reg187=reg47*reg187; T reg215=reg194+reg181; T reg216=reg193+reg184;
    T reg217=reg197+reg169; T reg218=reg206+reg210; reg168=reg47*reg168; reg178=reg47*reg178; T reg219=reg173+reg212;
    T reg220=reg199+reg175; T reg221=reg211+reg205; T reg222=reg203+reg172; T reg223=reg198+reg171; T reg224=reg202+reg170;
    T reg225=reg190+reg177; T reg226=reg201+reg176; T reg227=reg200+reg179; T reg228=reg195+reg182; T reg229=reg207+reg183;
    T reg230=reg196+reg186; T reg231=reg208+reg185; T reg232=reg191+reg188; T reg233=reg192+reg189; reg168=ponderation*reg168;
    T reg234=reg47*reg228; T reg235=reg47*reg225; T reg236=reg47*reg230; T reg237=reg47*reg215; T reg238=reg47*reg214;
    T reg239=reg47*reg223; T reg240=reg47*reg233; T reg241=reg47*reg220; T reg242=reg47*reg219; T reg243=reg47*reg221;
    T reg244=reg47*reg222; T reg245=reg47*reg224; T reg246=reg47*reg226; T reg247=reg47*reg227; T reg248=reg47*reg229;
    T reg249=reg47*reg231; T reg250=reg47*reg232; reg187=ponderation*reg187; T reg251=reg47*reg216; T reg252=reg47*reg218;
    reg178=ponderation*reg178; T reg253=reg47*reg217; T reg254=reg47*reg213; T reg255=ponderation*reg241; sollicitation[indices[7]+2]+=reg255;
    T reg256=ponderation*reg242; sollicitation[indices[0]+0]+=reg256; T reg257=ponderation*reg252; sollicitation[indices[4]+2]+=reg257; T reg258=ponderation*reg239;
    sollicitation[indices[7]+1]+=reg258; T reg259=ponderation*reg243; sollicitation[indices[0]+1]+=reg259; T reg260=ponderation*reg254; sollicitation[indices[4]+1]+=reg260;
    T reg261=ponderation*reg244; sollicitation[indices[0]+2]+=reg261; T reg262=ponderation*reg235; sollicitation[indices[7]+0]+=reg262; T reg263=ponderation*reg245;
    sollicitation[indices[1]+0]+=reg263; T reg264=ponderation*reg253; sollicitation[indices[4]+0]+=reg264; T reg265=ponderation*reg234; sollicitation[indices[6]+2]+=reg265;
    T reg266=ponderation*reg246; sollicitation[indices[1]+1]+=reg266; T reg267=ponderation*reg238; sollicitation[indices[5]+0]+=reg267; T reg268=ponderation*reg247;
    sollicitation[indices[1]+2]+=reg268; T reg269=ponderation*reg236; sollicitation[indices[6]+1]+=reg269; T reg270=ponderation*reg248; sollicitation[indices[2]+0]+=reg270;
    sollicitation[indices[3]+2]+=-reg178; reg178=ponderation*reg240; sollicitation[indices[6]+0]+=reg178; T reg271=ponderation*reg249; sollicitation[indices[2]+1]+=reg271;
    T reg272=ponderation*reg250; sollicitation[indices[2]+2]+=reg272; T reg273=ponderation*reg237; sollicitation[indices[5]+2]+=reg273; T reg274=ponderation*reg251;
    sollicitation[indices[3]+1]+=reg274; sollicitation[indices[3]+0]+=-reg187; sollicitation[indices[5]+1]+=-reg168;
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TMA,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim,unsigned symmetric_version>
void add_local_elem_matrix(TD ponderation,const TD *var_inter,
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TD,class T,class TM,bool wont_add_nz,class TVE,class TVEVE,class T_pos,class ND,class ED,unsigned nim>
void add_local_elem_residual( TD ponderation, const TD *var_inter,
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices ) { 
  #define PNODE(N) (*elem.node(N))
    T reg0=1-var_inter[2]; T reg1=1-var_inter[0]; T reg2=1-var_inter[1]; T reg3=reg0*reg1; T reg4=reg0*var_inter[0];
    T reg5=reg2*reg1; T reg6=reg2*var_inter[0]; T reg7=reg2*reg0; T reg8=elem.pos(1)[2]*reg7; T reg9=elem.pos(1)[1]*reg7;
    T reg10=elem.pos(0)[1]*reg5; T reg11=elem.pos(1)[1]*reg6; T reg12=elem.pos(0)[1]*reg7; T reg13=elem.pos(0)[2]*reg3; T reg14=elem.pos(1)[2]*reg4;
    T reg15=elem.pos(0)[1]*reg3; T reg16=elem.pos(1)[1]*reg4; T reg17=reg0*var_inter[1]; T reg18=elem.pos(0)[2]*reg7; T reg19=var_inter[0]*var_inter[1];
    T reg20=elem.pos(0)[2]*reg5; T reg21=reg6*elem.pos(1)[2]; T reg22=var_inter[1]*reg1; reg8=reg8-reg18; T reg23=elem.pos(2)[2]*reg17;
    T reg24=reg19*elem.pos(2)[2]; T reg25=reg21+reg20; T reg26=reg16+reg15; T reg27=elem.pos(2)[1]*reg4; T reg28=reg14+reg13;
    T reg29=elem.pos(2)[2]*reg4; T reg30=reg11+reg10; T reg31=elem.pos(2)[1]*reg17; reg9=reg9-reg12; T reg32=reg19*elem.pos(2)[1];
    T reg33=elem.pos(3)[2]*reg17; T reg34=elem.pos(0)[0]*reg3; T reg35=elem.pos(1)[0]*reg4; T reg36=reg30+reg32; T reg37=elem.pos(0)[0]*reg7;
    T reg38=elem.pos(1)[0]*reg7; reg23=reg8+reg23; reg8=reg25+reg24; T reg39=var_inter[2]*reg1; T reg40=reg22*elem.pos(3)[1];
    reg27=reg27-reg26; T reg41=elem.pos(3)[1]*reg3; T reg42=reg2*var_inter[2]; T reg43=elem.pos(3)[2]*reg3; T reg44=reg22*elem.pos(3)[2];
    T reg45=elem.pos(3)[1]*reg17; reg31=reg9+reg31; reg29=reg29-reg28; reg9=elem.pos(4)[1]*reg42; T reg46=elem.pos(4)[2]*reg5;
    T reg47=reg44+reg8; reg31=reg31-reg45; reg23=reg23-reg33; T reg48=elem.pos(4)[2]*reg42; T reg49=reg35+reg34;
    T reg50=elem.pos(2)[0]*reg4; T reg51=var_inter[2]*var_inter[0]; reg41=reg27+reg41; reg27=elem.pos(4)[1]*reg39; T reg52=elem.pos(4)[1]*reg5;
    T reg53=reg36+reg40; reg29=reg43+reg29; reg43=elem.pos(4)[2]*reg39; T reg54=elem.pos(1)[0]*reg6; T reg55=elem.pos(0)[0]*reg5;
    reg38=reg38-reg37; T reg56=elem.pos(2)[0]*reg17; T reg57=elem.pos(3)[0]*reg17; reg23=reg23-reg48; T reg58=elem.pos(5)[2]*reg42;
    T reg59=reg19*elem.pos(2)[0]; T reg60=reg54+reg55; T reg61=elem.pos(3)[0]*reg3; reg50=reg50-reg49; T reg62=elem.pos(5)[2]*reg51;
    reg29=reg29-reg43; reg46=reg46-reg47; reg41=reg41-reg27; T reg63=elem.pos(5)[1]*reg51; T reg64=reg6*elem.pos(5)[2];
    reg56=reg38+reg56; reg31=reg31-reg9; reg38=elem.pos(5)[1]*reg42; reg52=reg52-reg53; T reg65=reg6*elem.pos(5)[1];
    T reg66=var_inter[2]*var_inter[1]; reg65=reg52+reg65; reg52=reg19*elem.pos(6)[1]; T reg67=elem.pos(4)[0]*reg42; reg56=reg56-reg57;
    reg29=reg29-reg62; T reg68=elem.pos(6)[2]*reg51; T reg69=reg19*elem.pos(6)[2]; T reg70=elem.pos(4)[0]*reg39; reg64=reg46+reg64;
    reg50=reg61+reg50; reg46=reg22*elem.pos(3)[0]; reg38=reg31+reg38; reg31=elem.pos(6)[1]*reg66; reg61=reg60+reg59;
    T reg71=elem.pos(6)[2]*reg66; reg58=reg23+reg58; reg41=reg41-reg63; reg23=elem.pos(6)[1]*reg51; reg68=reg29+reg68;
    reg29=elem.pos(7)[2]*reg39; T reg72=reg46+reg61; reg52=reg65+reg52; reg65=elem.pos(4)[0]*reg5; T reg73=reg22*elem.pos(7)[1];
    T reg74=reg4*vectors[0][indices[1]+1]; T reg75=elem.pos(5)[0]*reg42; reg56=reg56-reg67; reg31=reg38+reg31; reg38=elem.pos(7)[1]*reg66;
    T reg76=reg7*vectors[0][indices[0]+0]; reg71=reg58+reg71; reg58=elem.pos(7)[2]*reg66; T reg77=reg4*vectors[0][indices[1]+0]; T reg78=reg3*vectors[0][indices[0]+1];
    T reg79=reg3*vectors[0][indices[0]+0]; reg50=reg50-reg70; T reg80=reg4*vectors[0][indices[1]+2]; T reg81=elem.pos(5)[0]*reg51; T reg82=reg3*vectors[0][indices[0]+2];
    T reg83=elem.pos(7)[1]*reg39; reg23=reg41+reg23; reg41=reg7*vectors[0][indices[1]+0]; T reg84=reg7*vectors[0][indices[1]+1]; T reg85=reg7*vectors[0][indices[0]+1];
    T reg86=reg7*vectors[0][indices[0]+2]; T reg87=reg7*vectors[0][indices[1]+2]; reg69=reg64+reg69; reg64=reg22*elem.pos(7)[2]; reg74=reg78+reg74;
    reg78=reg6*vectors[0][indices[1]+1]; T reg88=1+(*f.m).poisson_ratio; T reg89=reg17*vectors[0][indices[2]+0]; reg76=reg41-reg76; reg73=reg52+reg73;
    reg79=reg77+reg79; reg85=reg84-reg85; reg41=reg4*vectors[0][indices[2]+0]; reg64=reg69+reg64; reg82=reg80+reg82;
    reg52=reg5*vectors[0][indices[0]+0]; reg69=reg6*vectors[0][indices[1]+0]; reg77=reg17*vectors[0][indices[2]+1]; reg80=reg4*vectors[0][indices[2]+2]; reg84=reg6*vectors[0][indices[1]+2];
    reg83=reg23+reg83; reg86=reg87-reg86; reg31=reg31-reg38; reg75=reg56+reg75; reg23=elem.pos(6)[0]*reg66;
    reg29=reg68+reg29; reg65=reg65-reg72; reg56=reg6*elem.pos(5)[0]; reg71=reg71-reg58; reg50=reg50-reg81;
    reg68=reg5*vectors[0][indices[0]+1]; reg87=reg17*vectors[0][indices[2]+2]; T reg90=reg5*vectors[0][indices[0]+2]; T reg91=elem.pos(6)[0]*reg51; T reg92=reg4*vectors[0][indices[2]+1];
    reg74=reg92-reg74; reg92=reg3*vectors[0][indices[3]+2]; T reg93=reg19*vectors[0][indices[2]+0]; reg82=reg80-reg82; reg91=reg50+reg91;
    reg50=reg3*vectors[0][indices[3]+0]; reg76=reg89+reg76; reg80=reg17*vectors[0][indices[3]+0]; reg89=reg19*vectors[0][indices[2]+2]; reg90=reg84+reg90;
    reg84=reg19*vectors[0][indices[2]+1]; reg68=reg78+reg68; reg88=reg88/(*f.m).elastic_modulus; reg78=elem.pos(7)[0]*reg66; T reg94=reg3*vectors[0][indices[3]+1];
    reg23=reg75+reg23; reg75=reg17*vectors[0][indices[3]+2]; reg87=reg86+reg87; reg86=reg19*elem.pos(6)[0]; reg56=reg65+reg56;
    reg77=reg85+reg77; reg65=reg83*reg64; reg85=reg31*reg64; T reg95=reg29*reg73; reg79=reg41-reg79;
    reg41=reg71*reg73; T reg96=elem.pos(7)[0]*reg39; T reg97=reg17*vectors[0][indices[3]+1]; reg69=reg52+reg69; reg52=reg42*vectors[0][indices[4]+2];
    reg75=reg87-reg75; reg97=reg77-reg97; reg77=reg42*vectors[0][indices[4]+1]; reg87=reg22*vectors[0][indices[3]+0]; reg84=reg68+reg84;
    reg94=reg74+reg94; reg80=reg76-reg80; reg69=reg93+reg69; reg68=reg39*vectors[0][indices[4]+2]; reg92=reg82+reg92;
    reg74=reg42*vectors[0][indices[4]+0]; reg76=reg22*vectors[0][indices[3]+1]; reg82=reg39*vectors[0][indices[4]+0]; reg89=reg90+reg89; reg79=reg50+reg79;
    reg50=reg22*vectors[0][indices[3]+2]; reg41=reg85-reg41; reg85=reg31*reg29; reg95=reg65-reg95; reg65=reg71*reg83;
    reg90=pow(reg88,2); reg93=reg39*vectors[0][indices[4]+1]; T reg98=reg22*elem.pos(7)[0]; reg86=reg56+reg86; reg96=reg91+reg96;
    reg23=reg23-reg78; reg56=reg51*vectors[0][indices[5]+1]; reg68=reg92-reg68; reg91=reg96*reg41; reg92=reg42*vectors[0][indices[5]+0];
    T reg99=reg51*vectors[0][indices[5]+2]; T reg100=reg5*vectors[0][indices[4]+1]; T reg101=reg42*vectors[0][indices[5]+2]; reg88=reg88*reg90; reg82=reg79-reg82;
    reg79=reg23*reg95; reg98=reg86+reg98; reg86=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg102=1.0/(*f.m).elastic_modulus; reg50=reg89+reg50;
    reg89=reg42*vectors[0][indices[5]+1]; reg76=reg84+reg76; reg93=reg94-reg93; reg77=reg97-reg77; reg87=reg69+reg87;
    reg52=reg75-reg52; reg69=reg5*vectors[0][indices[4]+0]; reg75=reg5*vectors[0][indices[4]+2]; reg84=reg51*vectors[0][indices[5]+0]; reg65=reg85-reg65;
    reg74=reg80-reg74; reg80=reg23*reg64; reg85=reg83*reg98; reg94=reg23*reg73; reg84=reg82-reg84;
    reg50=reg75-reg50; reg75=reg71*reg98; reg56=reg93-reg56; reg82=reg51*vectors[0][indices[6]+1]; reg93=reg51*vectors[0][indices[6]+0];
    reg97=reg31*reg98; T reg103=reg66*vectors[0][indices[6]+0]; T reg104=reg66*vectors[0][indices[6]+1]; reg77=reg89+reg77; reg76=reg100-reg76;
    reg74=reg92+reg74; reg52=reg101+reg52; reg89=reg66*vectors[0][indices[6]+2]; reg73=reg96*reg73; reg92=reg29*reg98;
    reg64=reg96*reg64; reg98=reg98*reg65; reg100=reg6*vectors[0][indices[5]+2]; reg101=reg6*vectors[0][indices[5]+1]; reg91=reg79-reg91;
    reg87=reg69-reg87; reg69=reg86*reg88; reg79=reg6*vectors[0][indices[5]+0]; reg88=reg102*reg88; reg99=reg68-reg99;
    reg68=reg51*vectors[0][indices[6]+2]; T reg105=reg19*vectors[0][indices[6]+1]; reg76=reg101+reg76; reg101=reg102*reg90; reg92=reg64-reg92;
    reg64=reg66*vectors[0][indices[7]+1]; reg104=reg77+reg104; reg71=reg71*reg96; reg75=reg80-reg75; reg83=reg23*reg83;
    reg96=reg31*reg96; reg31=reg86*reg69; reg77=reg86*reg88; reg98=reg91+reg98; reg84=reg93+reg84;
    reg85=reg73-reg85; reg90=reg86*reg90; reg89=reg52+reg89; reg52=reg66*vectors[0][indices[7]+2]; reg73=reg19*vectors[0][indices[6]+0];
    reg79=reg87+reg79; reg80=reg19*vectors[0][indices[6]+2]; reg50=reg100+reg50; reg87=reg66*vectors[0][indices[7]+0]; reg29=reg23*reg29;
    reg97=reg94-reg97; reg23=reg39*vectors[0][indices[7]+1]; reg56=reg82+reg56; reg82=reg39*vectors[0][indices[7]+2]; reg68=reg99+reg68;
    reg88=reg102*reg88; reg74=reg103+reg74; reg91=reg39*vectors[0][indices[7]+0]; reg92=reg92/reg98; reg95=reg95/reg98;
    reg88=reg88-reg31; reg80=reg50+reg80; reg41=reg41/reg98; reg50=reg22*vectors[0][indices[7]+2]; reg93=reg86*reg101;
    reg85=reg85/reg98; reg94=reg86*reg90; reg101=reg102*reg101; reg77=reg31+reg77; reg69=reg102*reg69;
    reg87=reg74-reg87; reg74=reg22*vectors[0][indices[7]+0]; reg82=reg68+reg82; reg73=reg79+reg73; reg52=reg89-reg52;
    reg64=reg104-reg64; reg56=reg23+reg56; reg23=reg22*vectors[0][indices[7]+1]; reg71=reg29-reg71; reg96=reg83-reg96;
    reg91=reg84+reg91; reg75=reg75/reg98; reg105=reg76+reg105; reg97=reg97/reg98; reg29=reg41*reg91;
    reg101=reg101-reg94; reg68=reg75*reg56; reg90=reg102*reg90; reg93=reg94+reg93; reg96=reg96/reg98;
    reg105=reg23+reg105; reg23=reg86*reg77; reg74=reg73+reg74; reg65=reg65/reg98; reg73=reg102*reg88;
    reg76=reg75*reg91; reg79=reg92*reg64; reg83=reg92*reg87; reg84=reg95*reg64; reg89=reg41*reg56;
    reg99=reg95*reg87; reg100=reg41*reg82; reg103=reg95*reg52; reg71=reg71/reg98; reg69=reg31+reg69;
    reg50=reg80+reg50; reg91=reg97*reg91; reg87=reg85*reg87; reg31=reg85*reg52; reg80=reg65*reg74;
    reg23=reg73-reg23; reg73=reg86*reg69; reg104=reg75*reg82; reg52=reg92*reg52; reg82=reg97*reg82;
    T reg106=reg71*reg105; T reg107=reg65*reg105; reg89=reg84-reg89; reg79=reg68-reg79; reg68=reg71*reg74;
    reg83=reg76-reg83; reg56=reg97*reg56; reg64=reg85*reg64; reg76=reg65*reg50; reg100=reg103-reg100;
    reg101=reg102*reg101; reg93=reg86*reg93; reg74=reg96*reg74; reg84=reg94+reg90; reg91=reg87-reg91;
    reg29=reg99-reg29; reg87=reg96*reg50; reg82=reg31-reg82; reg106=reg79-reg106; reg93=reg101-reg93;
    reg29=reg80+reg29; reg107=reg89+reg107; reg31=(*f.m).deltaT*(*f.m).alpha; reg84=reg86*reg84; reg68=reg83-reg68;
    reg56=reg64-reg56; reg50=reg71*reg50; reg52=reg104-reg52; reg105=reg96*reg105; reg76=reg100+reg76;
    reg74=reg91+reg74; reg73=reg23-reg73; reg23=reg66*reg92; reg106=reg106-reg31; reg64=reg66*reg95;
    reg79=reg4*reg41; reg80=reg4*reg75; reg83=reg7*reg95; reg86=reg3*reg41; reg82=reg87+reg82;
    reg87=reg7*reg92; reg89=reg17*reg95; reg91=reg17*reg92; reg76=reg74+reg76; reg107=reg68+reg107;
    reg68=reg51*reg75; reg74=reg51*reg41; reg99=reg39*reg75; reg100=reg42*reg92; reg101=reg42*reg95;
    reg102=reg39*reg41; reg103=reg3*reg75; reg84=reg93-reg84; reg56=reg105+reg56; reg50=reg52-reg50;
    reg69=reg69/reg73; reg29=reg29-reg31; reg88=reg88/reg73; reg77=reg77/reg73; reg107=0.5*reg107;
    reg52=reg89-reg79; reg73=reg84/reg73; reg84=reg4*reg97; reg93=reg87-reg103; reg104=reg19*reg71;
    reg105=reg74+reg101; T reg108=reg39*reg97; T reg109=reg77*reg106; T reg110=reg51*reg97; T reg111=reg102+reg64;
    T reg112=reg22*reg65; T reg113=reg86+reg89; T reg114=reg100-reg99; T reg115=reg102-reg101; T reg116=reg99+reg23;
    T reg117=reg22*reg71; T reg118=reg103+reg91; T reg119=reg87+reg80; T reg120=reg88*reg29; T reg121=reg42*reg85;
    T reg122=reg80-reg91; T reg123=reg6*reg71; T reg124=reg83+reg79; T reg125=reg6*reg65; reg50=reg56+reg50;
    reg56=reg68-reg23; T reg126=reg69*reg106; T reg127=reg64-reg74; reg106=reg88*reg106; T reg128=reg7*reg85;
    T reg129=reg3*reg97; reg76=0.5*reg76; reg29=reg77*reg29; T reg130=reg19*reg65; T reg131=reg66*reg85;
    reg82=reg82-reg31; T reg132=reg5*reg65; T reg133=reg5*reg71; T reg134=reg86-reg83; T reg135=reg17*reg85;
    T reg136=reg68+reg100; T reg137=reg131-reg110; reg127=reg127+reg130; T reg138=reg123-reg119; reg118=reg118+reg117;
    reg56=reg56-reg104; reg115=reg115+reg132; reg114=reg114-reg133; T reg139=reg88*reg82; reg126=reg29+reg126;
    T reg140=reg113+reg112; T reg141=reg22*reg96; T reg142=reg110+reg121; reg106=reg29+reg106; reg109=reg120+reg109;
    reg105=reg105+reg125; reg29=reg108-reg121; reg82=reg69*reg82; reg122=reg122+reg104; reg52=reg52-reg130;
    reg120=reg129+reg135; T reg143=reg19*reg96; reg107=reg73*reg107; T reg144=reg6*reg96; T reg145=reg84+reg128;
    reg116=reg116-reg117; T reg146=reg112-reg111; reg124=reg124-reg125; reg93=reg93+reg133; T reg147=reg108+reg131;
    T reg148=reg129-reg128; T reg149=reg5*reg96; reg134=reg134-reg132; reg50=0.5*reg50; T reg150=reg136+reg123;
    reg76=reg73*reg76; T reg151=reg135-reg84; T reg152=0.5*reg127; reg139=reg126+reg139; reg137=reg137+reg143;
    reg142=reg144+reg142; reg126=0.5*reg105; T reg153=0.5*reg124; reg107=2*reg107; T reg154=0.5*reg116;
    reg76=2*reg76; T reg155=0.5*reg115; T reg156=0.5*reg52; T reg157=0.5*reg134; reg145=reg145-reg144;
    reg29=reg29+reg149; reg50=reg73*reg50; T reg158=0.5*reg146; T reg159=reg141+reg120; T reg160=0.5*reg56;
    T reg161=0.5*reg122; reg106=reg82+reg106; T reg162=0.5*reg93; reg151=reg151-reg143; T reg163=0.5*reg150;
    T reg164=0.5*reg114; T reg165=reg141-reg147; T reg166=0.5*reg138; reg82=reg109+reg82; reg148=reg148-reg149;
    reg109=0.5*reg118; T reg167=0.5*reg140; T reg168=reg165*reg139; T reg169=reg158*reg76; T reg170=reg107*reg161;
    T reg171=reg142*reg139; T reg172=reg126*reg76; T reg173=reg157*reg76; T reg174=reg116*reg106; T reg175=reg148*reg139;
    T reg176=reg158*reg107; T reg177=0.5*reg137; T reg178=reg124*reg82; T reg179=reg82*reg127; T reg180=reg82*reg52;
    T reg181=0.5*reg148; T reg182=reg162*reg107; T reg183=reg134*reg82; T reg184=reg107*reg160; T reg185=0.5*reg151;
    reg50=2*reg50; T reg186=reg109*reg107; T reg187=reg114*reg106; T reg188=reg164*reg107; T reg189=reg93*reg106;
    T reg190=reg157*reg107; T reg191=reg145*reg139; T reg192=reg153*reg76; T reg193=reg115*reg82; T reg194=0.5*reg165;
    T reg195=reg139*reg159; T reg196=reg155*reg76; T reg197=reg107*reg166; T reg198=reg29*reg139; T reg199=reg152*reg76;
    T reg200=reg137*reg139; T reg201=reg152*reg107; T reg202=reg156*reg76; T reg203=reg151*reg139; T reg204=reg106*reg56;
    T reg205=0.5*reg159; T reg206=reg107*reg167; T reg207=reg150*reg106; T reg208=reg126*reg107; T reg209=reg156*reg107;
    T reg210=reg106*reg122; T reg211=0.5*reg29; T reg212=reg163*reg107; T reg213=reg118*reg106; T reg214=0.5*reg145;
    T reg215=reg105*reg82; T reg216=reg154*reg107; T reg217=reg153*reg107; T reg218=reg106*reg138; T reg219=0.5*reg142;
    T reg220=reg155*reg107; T reg221=reg140*reg82; T reg222=reg146*reg82; T reg223=reg76*reg167; reg196=reg198+reg196;
    reg198=reg164*reg50; reg187=reg220+reg187; reg220=reg22*var_inter[2]; T reg224=reg6*reg0; T reg225=var_inter[2]*reg19;
    T reg226=reg5*reg0; T reg227=reg0*reg19; reg197=reg178+reg197; reg178=reg214*reg76; reg222=reg216+reg222;
    reg216=reg194*reg76; T reg228=reg194*reg50; T reg229=reg50*reg161; reg202=reg203+reg202; reg203=reg6*var_inter[2];
    T reg230=var_inter[2]*reg5; reg186=reg186-reg221; T reg231=reg22*reg0; T reg232=reg211*reg76; T reg233=reg50*reg205;
    T reg234=reg76*reg205; reg169=reg168+reg169; reg168=reg154*reg50; reg174=reg176+reg174; reg176=reg211*reg50;
    reg180=reg170+reg180; reg170=reg185*reg76; reg179=reg184+reg179; reg184=reg177*reg76; T reg235=reg219*reg50;
    T reg236=reg214*reg50; reg208=reg208-reg207; T reg237=reg50*reg177; reg173=reg175+reg173; reg213=reg213-reg206;
    reg204=reg201+reg204; reg193=reg188+reg193; reg175=reg50*reg160; reg218=reg217+reg218; reg188=reg50*reg162;
    reg183=reg182+reg183; reg182=reg181*reg76; reg210=reg209+reg210; reg201=reg163*reg50; reg209=reg185*reg50;
    reg217=reg109*reg50; T reg238=reg195+reg223; reg199=reg200+reg199; reg215=reg215-reg212; reg200=reg50*reg181;
    reg189=reg190+reg189; reg190=reg219*reg76; reg192=reg191+reg192; reg191=reg50*reg166; reg172=reg171+reg172;
    reg186=reg186-reg234; reg171=reg220*(*f.m).f_vol[0]; reg216=reg222+reg216; reg192=reg191+reg192; reg191=reg220*(*f.m).f_vol[2];
    reg222=reg231*(*f.m).f_vol[0]; T reg239=reg224*(*f.m).f_vol[2]; T reg240=reg203*(*f.m).f_vol[1]; reg208=reg235+reg208; reg204=reg237+reg204;
    reg235=reg230*(*f.m).f_vol[0]; reg193=reg232+reg193; reg210=reg209+reg210; reg209=reg227*(*f.m).f_vol[1]; reg232=reg225*(*f.m).f_vol[1];
    reg190=reg215+reg190; reg215=reg203*(*f.m).f_vol[0]; reg202=reg229+reg202; reg229=reg227*(*f.m).f_vol[2]; reg237=reg226*(*f.m).f_vol[1];
    T reg241=reg231*(*f.m).f_vol[2]; reg217=reg217-reg238; reg189=reg200+reg189; reg200=reg230*(*f.m).f_vol[2]; reg196=reg198+reg196;
    reg198=(*f.m).f_vol[0]*reg226; reg183=reg182+reg183; reg182=reg230*(*f.m).f_vol[1]; reg187=reg176+reg187; reg176=reg231*(*f.m).f_vol[1];
    reg213=reg213-reg233; reg173=reg188+reg173; reg188=reg226*(*f.m).f_vol[2]; reg169=reg168+reg169; reg172=reg172-reg201;
    reg168=reg203*(*f.m).f_vol[2]; T reg242=reg224*(*f.m).f_vol[1]; reg218=reg236+reg218; reg236=reg220*(*f.m).f_vol[1]; reg174=reg228+reg174;
    reg199=reg175+reg199; reg175=reg224*(*f.m).f_vol[0]; reg178=reg197+reg178; reg197=reg225*(*f.m).f_vol[2]; reg228=(*f.m).f_vol[0]*reg227;
    reg180=reg170+reg180; reg184=reg179+reg184; reg170=reg225*(*f.m).f_vol[0]; reg202=reg202-reg229; reg204=reg204-reg232;
    reg192=reg192-reg239; reg216=reg216-reg171; reg218=reg218-reg242; reg178=reg178-reg175; reg173=reg173-reg188;
    reg183=reg183-reg198; reg189=reg189-reg237; reg169=reg169-reg191; reg172=reg172-reg168; reg208=reg208-reg240;
    reg174=reg174-reg236; reg186=reg186-reg222; reg180=reg180-reg228; reg184=reg184-reg170; reg187=reg187-reg182;
    reg193=reg193-reg235; reg210=reg210-reg209; reg190=reg190-reg215; reg213=reg213-reg176; reg217=reg217-reg241;
    reg199=reg199-reg197; reg196=reg196-reg200; reg217=reg98*reg217; reg218=reg98*reg218; reg174=reg98*reg174;
    reg196=reg98*reg196; reg172=reg98*reg172; reg178=reg98*reg178; reg180=reg98*reg180; reg173=reg98*reg173;
    reg183=reg98*reg183; reg213=reg98*reg213; reg184=reg98*reg184; reg189=reg98*reg189; reg187=reg98*reg187;
    reg193=reg98*reg193; reg204=reg204*reg98; reg210=reg98*reg210; reg190=reg98*reg190; reg192=reg98*reg192;
    reg186=reg98*reg186; reg216=reg98*reg216; reg208=reg98*reg208; reg199=reg98*reg199; reg202=reg98*reg202;
    reg169=reg98*reg169; sollicitation[indices[6]+2]+=ponderation*reg199; sollicitation[indices[0]+0]+=ponderation*reg183; sollicitation[indices[3]+2]+=ponderation*reg217; sollicitation[indices[3]+1]+=ponderation*reg213;
    sollicitation[indices[5]+0]+=ponderation*reg190; sollicitation[indices[4]+2]+=ponderation*reg196; sollicitation[indices[6]+1]+=ponderation*reg204; sollicitation[indices[4]+1]+=ponderation*reg187; sollicitation[indices[2]+2]+=ponderation*reg202;
    sollicitation[indices[0]+1]+=ponderation*reg189; sollicitation[indices[4]+0]+=ponderation*reg193; sollicitation[indices[6]+0]+=ponderation*reg184; sollicitation[indices[2]+1]+=ponderation*reg210; sollicitation[indices[3]+0]+=ponderation*reg186;
    sollicitation[indices[1]+2]+=ponderation*reg192; sollicitation[indices[0]+2]+=ponderation*reg173; sollicitation[indices[2]+0]+=ponderation*reg180; sollicitation[indices[5]+2]+=ponderation*reg172; sollicitation[indices[1]+0]+=ponderation*reg178;
    sollicitation[indices[7]+0]+=ponderation*reg216; sollicitation[indices[7]+1]+=ponderation*reg174; sollicitation[indices[5]+1]+=ponderation*reg208; sollicitation[indices[7]+2]+=ponderation*reg169; sollicitation[indices[1]+1]+=ponderation*reg218;
  #undef PNODE
}
// 
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_false_true_true
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_false_true_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_true_true_true
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_true_true_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_false_true_false
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_false_true_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_true_true_false
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_true_true_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE >
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_symmetric_version_false_true
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_symmetric_version_false_true
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE ,unsigned symmetric_version>
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_symmetric_version_false_false
#define ADD_NODAL_MATRIX_elasticity_isotropy_stat_Qstat_symmetric_version_false_false
template<class TM,class T,bool wont_add_nz,class TMA,class TVE,class TVEVE ,unsigned symmetric_version>
void add_nodal_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ADD_NODAL_RESIDUAL_elasticity_isotropy_stat_Qstat
#define ADD_NODAL_RESIDUAL_elasticity_isotropy_stat_Qstat
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE>
void add_nodal_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<0> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<1> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<2> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<3> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<4> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<4> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<4> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<4> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<4> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<4> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<4> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<5> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<5> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<5> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<5> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<5> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TMA, class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2,unsigned symmetric_version>
void add_skin_elem_matrix(
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<5> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}
// 
template<class TM,class T,bool wont_add_nz,class TVE,class TVEVE, class T_pos, class ND,class ED, unsigned nim,class ED2,unsigned nim2>
void add_skin_elem_residual(
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Hexa,DefaultBehavior,Node<3,T_pos,ND>,ED,nim> &elem,
      const Element<Quad,DefaultBehavior,Node<3,T_pos,ND>,ED2,nim2> &skin_elem,
      const Number<5> &num_child,
      const unsigned *indices){
   #define PNODE(N) (*elem.node(N))
  #undef PNODE
}

#ifndef elasticity_isotropy_stat_Qstat_read_material_to_mesh
#define elasticity_isotropy_stat_Qstat_read_material_to_mesh
template<class TM, class T, bool wont_add_nz>
void read_material_to_mesh_(const XmlNode &n, Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f){ 
    if(n.has_attribute("elastic_modulus"))  
        n.get_attribute("elastic_modulus", f.m->elastic_modulus ); 
    else  
        std::cerr << "Warning using default value of elastic_modulus : " << f.m->elastic_modulus << std::endl; 

    if(n.has_attribute("density"))  
        n.get_attribute("density", f.m->density ); 
    else  
        std::cerr << "Warning using default value of density : " << f.m->density << std::endl; 

    if(n.has_attribute("deltaT"))  
        n.get_attribute("deltaT", f.m->deltaT ); 
    else  
        std::cerr << "Warning using default value of deltaT : " << f.m->deltaT << std::endl; 

    if(n.has_attribute("poisson_ratio"))  
        n.get_attribute("poisson_ratio", f.m->poisson_ratio ); 
    else  
        std::cerr << "Warning using default value of poisson_ratio : " << f.m->poisson_ratio << std::endl; 

    if(n.has_attribute("alpha"))  
        n.get_attribute("alpha", f.m->alpha ); 
    else  
        std::cerr << "Warning using default value of alpha : " << f.m->alpha << std::endl; 

    if(n.has_attribute("resolution"))  
        n.get_attribute("resolution", f.m->resolution ); 
    else  
        std::cerr << "Warning using default value of resolution : " << f.m->resolution << std::endl; 

    if(n.has_attribute("f_vol"))  
        n.get_attribute("f_vol", f.m->f_vol ); 
    else  
        std::cerr << "Warning using default value of f_vol : " << f.m->f_vol << std::endl; 

  };
#endif // elasticity_isotropy_stat_Qstat_read_material_to_mesh
} // namespace LMT

