
#include "formulation/formulation.h"
namespace LMT {
#ifndef ELASTICITY_ISOTROPY_STAT_QSTAT_2_NUM_0_CARAC_H
#define ELASTICITY_ISOTROPY_STAT_QSTAT_2_NUM_0_CARAC_H
#ifndef STRUCT_ELASTICITY_ISOTROPY_STAT_QSTAT
#define STRUCT_ELASTICITY_ISOTROPY_STAT_QSTAT
struct elasticity_isotropy_stat_Qstat {
  static const char *name() { return "elasticity_isotropy_stat_Qstat"; }
};
#endif // ELASTICITY_ISOTROPY_STAT_QSTAT

template<class P_T>
class CaracFormulation<elasticity_isotropy_stat_Qstat,2,P_T>  {
public:
  typedef P_T T;
  static const char *name() { return "elasticity_isotropy_stat_Qstat"; }
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
    T reg0=vecs[1][indice+0]-vecs[0][indice+0]; T reg1=vecs[1][indice+1]-vecs[0][indice+1]; reg0=abs(reg0); reg1=abs(reg1); return max(reg1,reg0);
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
#endif // ELASTICITY_ISOTROPY_STAT_QSTAT_2_NUM_0_CARAC_H
extern "C" void apply_on_elements_after_solve_0_elasticity_isotropy_stat_Qstat_Triangle_0( double * );
extern "C" void apply_on_elements_after_solve_1_elasticity_isotropy_stat_Qstat_Triangle_1( double * );
extern "C" void apply_on_elements_after_solve_2_elasticity_isotropy_stat_Qstat_Triangle_2( double * );
extern "C" void apply_on_elements_after_solve_3_elasticity_isotropy_stat_Qstat_Triangle_3( double * );
extern "C" void apply_on_elements_after_solve_4_elasticity_isotropy_stat_Qstat_Triangle_4( double * );
extern "C" void apply_on_elements_after_solve_5_elasticity_isotropy_stat_Qstat_Triangle_5( double * );
extern "C" void apply_on_elements_after_solve_6_elasticity_isotropy_stat_Qstat_Triangle_6( double * );
extern "C" void apply_on_elements_after_solve_7_elasticity_isotropy_stat_Qstat_Triangle_7( double * );
extern "C" void apply_on_elements_after_solve_8_elasticity_isotropy_stat_Qstat_Triangle_8( double * );
extern "C" void apply_on_elements_after_solve_9_elasticity_isotropy_stat_Qstat_Triangle_9( double * );
extern "C" void apply_on_elements_after_solve_10_elasticity_isotropy_stat_Qstat_Triangle_10( double * );
extern "C" void apply_on_elements_after_solve_11_elasticity_isotropy_stat_Qstat_Triangle_11( double * );
extern "C" void apply_on_elements_after_solve_12_elasticity_isotropy_stat_Qstat_Triangle_12( double * );
extern "C" void apply_on_elements_after_solve_13_elasticity_isotropy_stat_Qstat_Triangle_13( double * );
extern "C" void apply_on_elements_after_solve_14_elasticity_isotropy_stat_Qstat_Triangle_14( double * );
class Triangle;
template<unsigned A,class B,class C> class Node;
template<class A,class B,class C,class D,unsigned E> class Element;

// Carac for ...
template<class P_T_pos,class P_ND,class TED,unsigned nim,class TM,class T>
class CaracFormulationForElement<elasticity_isotropy_stat_Qstat,Element<Triangle,DefaultBehavior,Node<2,P_T_pos,P_ND>,TED,nim>,TM,T> {
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=reg1*reg0; T reg5=reg3*reg4; T reg6=reg3*reg1; reg1=reg1*reg2; reg4=reg4*reg2;
    T reg7=reg3*reg5; T reg8=reg3*reg4; T reg9=reg1*reg3; reg4=reg4*reg2; reg1=reg1*reg2;
    T reg10=reg3*reg6; reg6=reg6*reg2; reg1=reg1-reg10; reg4=reg4-reg7; reg8=reg8+reg7;
    reg5=reg5*reg2; reg9=reg9+reg10; T reg11=reg4*reg2; T reg12=reg6+reg10; reg1=reg1*reg2;
    T reg13=reg3*reg8; reg9=reg3*reg9; reg5=reg7+reg5; reg9=reg1-reg9; reg1=reg3*reg5;
    reg13=reg11-reg13; reg12=reg12*reg3; reg1=reg13-reg1; reg12=reg9-reg12; reg8=reg8/reg1;
    reg4=reg4/reg1; reg12=reg12/reg1; reg7=reg4*reg12; reg9=reg8*reg12; reg11=(*f.m).alpha*reg8;
    reg13=elem.pos(2)[0]-elem.pos(0)[0]; T reg14=elem.pos(1)[0]-elem.pos(0)[0]; T reg15=reg7*reg4; T reg16=reg8*reg9; T reg17=(*f.m).alpha*reg4;
    reg1=reg5/reg1; reg5=elem.pos(1)[1]-elem.pos(0)[1]; T reg18=elem.pos(2)[1]-elem.pos(0)[1]; reg16=reg15-reg16; reg1=(*f.m).alpha*reg1;
    reg15=reg14*reg18; reg17=reg11+reg17; reg11=reg13*reg5; reg1=reg17+reg1; reg9=reg9/reg16;
    reg11=reg15-reg11; reg16=reg7/reg16; reg13=reg13/reg11; reg5=reg5/reg11; reg14=reg14/reg11;
    reg7=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg18=reg18/reg11; reg15=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg17=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; T reg19=vectors[0][indices[2]+0]-vectors[0][indices[0]+0];
    T reg20=reg3*reg0; T reg21=reg2*reg0; reg9=reg1*reg9; reg16=reg1*reg16; reg1=reg17*reg18;
    T reg22=1-(*f.m).resolution; T reg23=reg5*reg19; reg9=reg16-reg9; reg16=reg15*reg13; T reg24=reg21*reg2;
    T reg25=reg7*reg14; T reg26=reg3*reg20; reg26=reg24-reg26; reg24=PNODE(1).dep[0]-PNODE(0).dep[0]; T reg27=PNODE(2).dep[0]-PNODE(0).dep[0];
    T reg28=pow(reg3,2); T reg29=PNODE(1).dep[1]-PNODE(0).dep[1]; T reg30=pow(reg2,2); reg9=reg22*reg9; reg23=reg1-reg23;
    elem.epsilon[0][0]=reg23; reg1=(*f.m).alpha*(*f.m).deltaT; reg16=reg25-reg16; elem.epsilon[0][1]=reg16; reg25=(*f.m).alpha*(*f.m).resolution;
    T reg31=PNODE(2).dep[1]-PNODE(0).dep[1]; T reg32=reg14*reg27; reg28=reg30-reg28; reg30=reg24*reg13; reg20=reg20/reg26;
    reg21=reg21/reg26; T reg33=reg29*reg13; T reg34=reg31*reg14; reg27=reg5*reg27; reg29=reg29*reg18;
    reg24=reg24*reg18; reg25=reg9+reg25; reg9=reg16-reg1; T reg35=reg23-reg1; reg31=reg31*reg5;
    T reg36=reg9*reg20; T reg37=reg35*reg21; reg30=reg32-reg30; reg32=(*f.m).resolution*reg21; T reg38=(*f.m).resolution*reg20;
    reg25=(*f.m).deltaT*reg25; reg33=reg34-reg33; reg27=reg24-reg27; reg31=reg29-reg31; reg35=reg35*reg20;
    reg9=reg9*reg21; reg26=reg28/reg26; reg4=reg22*reg4; reg8=reg22*reg8; reg15=reg15*reg18;
    reg9=reg35+reg9; reg7=reg7*reg5; reg12=reg22*reg12; reg22=reg33-reg25; reg31=reg30+reg31;
    reg36=reg37+reg36; reg17=reg17*reg13; reg19=reg14*reg19; reg24=reg27-reg25; reg38=reg8+reg38;
    reg32=reg4+reg32; reg4=(*f.m).resolution*reg26; reg8=reg32*reg22; reg7=reg15-reg7; reg15=reg38*reg22;
    reg28=reg24*reg38; reg29=reg36+reg9; reg12=reg4+reg12; reg31=0.5*reg31; reg17=reg19-reg17;
    reg4=reg24*reg32; reg29=reg29/3; reg8=reg28+reg8; reg19=reg12*reg31; reg15=reg4+reg15;
    reg7=reg17+reg7; reg7=0.5*reg7; elem.epsilon[0][2]=reg7; reg19=2*reg19; reg24=reg24*reg15;
    reg22=reg8*reg22; reg36=reg36-reg29; reg9=reg9-reg29; reg24=reg22+reg24; reg4=reg31*reg19;
    reg36=pow(reg36,2); reg17=reg7*reg26; reg9=pow(reg9,2); reg4=reg24+reg4; reg36=reg9+reg36;
    reg9=2*reg17; reg29=pow(reg29,2); reg4=reg11*reg4; reg36=reg29+reg36; reg9=reg17*reg9;
    reg23=reg23-reg25; reg25=reg16-reg25; reg11=0.33333333333333331483*reg4; reg9=reg36+reg9; reg4=0.16666666666666665741*reg4;
    reg16=reg25*reg38; reg25=reg25*reg32; reg38=reg38*reg23; reg23=reg32*reg23; reg9=1.5*reg9;
    reg11=reg4+reg11; elem.sigma[0][1]=reg25+reg38; elem.ener=reg11/2; elem.sigma_von_mises=pow(reg9,0.5); elem.sigma[0][0]=reg16+reg23;
    elem.sigma[0][2]=reg7*reg12;
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
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg1*reg0; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg3*reg2; T reg6=reg3*reg1; reg2=reg2*reg4; reg1=reg1*reg4;
    T reg7=reg3*reg2; reg2=reg2*reg4; T reg8=reg1*reg3; T reg9=reg3*reg6; reg1=reg1*reg4;
    T reg10=reg3*reg5; reg8=reg8+reg9; reg1=reg1-reg9; reg6=reg6*reg4; reg7=reg7+reg10;
    reg5=reg5*reg4; reg2=reg2-reg10; reg8=reg3*reg8; reg1=reg1*reg4; T reg11=reg3*reg7;
    T reg12=reg2*reg4; reg5=reg10+reg5; reg10=reg6+reg9; T reg13=reg3*reg5; reg11=reg12-reg11;
    reg8=reg1-reg8; reg10=reg10*reg3; reg13=reg11-reg13; reg10=reg8-reg10; reg7=reg7/reg13;
    reg10=reg10/reg13; reg2=reg2/reg13; reg1=reg2*reg10; reg8=reg7*reg10; reg11=reg1*reg2;
    reg12=reg7*reg8; T reg14=(*f.m).alpha*reg7; reg13=reg5/reg13; reg5=(*f.m).alpha*reg2; reg12=reg11-reg12;
    reg13=(*f.m).alpha*reg13; reg5=reg14+reg5; reg8=reg8/reg12; reg12=reg1/reg12; reg1=reg3*reg0;
    reg13=reg5+reg13; reg0=reg4*reg0; reg5=reg0*reg4; reg11=elem.pos(2)[0]-elem.pos(0)[0]; reg14=elem.pos(1)[0]-elem.pos(0)[0];
    T reg15=elem.pos(2)[1]-elem.pos(0)[1]; T reg16=elem.pos(1)[1]-elem.pos(0)[1]; reg4=pow(reg4,2); T reg17=pow(reg3,2); reg3=reg3*reg1;
    reg12=reg13*reg12; reg8=reg13*reg8; reg13=1-(*f.m).resolution; reg17=reg4-reg17; reg3=reg5-reg3;
    reg8=reg12-reg8; reg4=reg14*reg15; reg5=reg11*reg16; reg17=reg17/reg3; reg5=reg4-reg5;
    reg1=reg1/reg3; reg8=reg13*reg8; reg4=(*f.m).alpha*(*f.m).resolution; reg3=reg0/reg3; reg0=(*f.m).resolution*reg17;
    reg12=(*f.m).resolution*reg3; T reg18=(*f.m).resolution*reg1; reg10=reg13*reg10; reg7=reg13*reg7; reg4=reg8+reg4;
    reg2=reg13*reg2; reg15=reg15/reg5; reg14=reg14/reg5; reg16=reg16/reg5; reg11=reg11/reg5;
    reg8=reg11-reg14; reg10=reg0+reg10; reg0=0.5*reg16; reg13=0.5*reg15; reg12=reg2+reg12;
    reg18=reg7+reg18; reg2=0.5*reg11; reg7=0.5*reg14; reg4=(*f.m).deltaT*reg4; T reg19=reg16-reg15;
    T reg20=reg10*reg2; T reg21=0.5*reg19; T reg22=0.5*reg8; T reg23=reg0*reg10; T reg24=reg10*reg13;
    T reg25=reg10*reg7; T reg26=reg12*reg4; T reg27=reg18*reg4; T reg28=2*reg20; T reg29=reg15*reg12;
    reg25=2*reg25; T reg30=reg15*reg18; T reg31=reg26+reg27; T reg32=reg18*reg16; reg24=2*reg24;
    T reg33=reg10*reg21; T reg34=reg18*reg11; T reg35=reg12*reg16; T reg36=1-var_inter[0]; T reg37=reg12*reg11;
    T reg38=reg12*reg14; T reg39=reg22*reg10; T reg40=reg18*reg14; T reg41=2*reg23; T reg42=reg25*reg2;
    T reg43=reg35*reg15; reg36=reg36-var_inter[1]; T reg44=reg11*reg30; T reg45=reg28*reg13; T reg46=var_inter[1]*elem.f_vol_e[0];
    T reg47=reg38*reg11; T reg48=reg15*reg34; T reg49=reg2*reg24; T reg50=reg41*reg13; T reg51=reg8*reg12;
    T reg52=reg40*reg16; T reg53=reg41*reg7; T reg54=reg31*reg16; T reg55=reg0*reg24; T reg56=reg37*reg14;
    reg33=2*reg33; reg39=2*reg39; T reg57=reg8*reg18; T reg58=reg19*reg12; T reg59=reg31*reg11;
    T reg60=var_inter[0]*elem.f_vol_e[1]; T reg61=reg19*reg18; T reg62=reg28*reg7; T reg63=reg32*reg14; T reg64=reg0*reg25;
    T reg65=reg16*reg29; T reg66=reg32*reg11; T reg67=var_inter[0]*elem.f_vol_e[0]; T reg68=reg28*reg2; T reg69=reg15*reg29;
    T reg70=reg19*reg40; T reg71=reg36*elem.f_vol_e[0]; T reg72=reg22*reg41; T reg73=reg58*reg16; T reg74=reg15*reg31;
    T reg75=reg39*reg7; T reg76=reg13*reg25; reg55=reg56+reg55; T reg77=reg22*reg24; T reg78=reg19*reg34;
    T reg79=reg35*reg19; T reg80=reg21*reg24; T reg81=reg41*reg0; T reg82=reg38*reg14; T reg83=reg41*reg2;
    reg40=reg40*reg15; T reg84=reg22*reg39; T reg85=reg37*reg8; T reg86=reg19*reg58; reg42=reg43+reg42;
    T reg87=reg33*reg7; T reg88=reg57*reg16; reg38=reg38*reg8; T reg89=reg36*elem.f_vol_e[1]; T reg90=reg39*reg0;
    T reg91=reg61*reg14; T reg92=reg41*reg21; reg44=reg45+reg44; reg65=reg62+reg65; T reg93=reg51*reg11;
    T reg94=reg22*reg25; T reg95=reg28*reg21; T reg96=reg59-reg60; T reg97=reg8*reg30; T reg98=reg22*reg28;
    reg58=reg58*reg15; T reg99=var_inter[1]*elem.f_vol_e[1]; reg29=reg19*reg29; T reg100=reg39*reg13; T reg101=reg16*reg34;
    T reg102=reg61*reg11; T reg103=reg19*reg31; T reg104=reg33*reg13; reg52=reg53+reg52; T reg105=reg7*reg24;
    T reg106=reg39*reg21; T reg107=reg31*reg14; T reg108=reg37*reg11; T reg109=reg35*reg16; T reg110=reg57*reg15;
    reg49=reg48+reg49; reg24=reg13*reg24; reg64=reg63+reg64; reg61=reg61*reg8; T reg111=reg51*reg8;
    reg57=reg19*reg57; T reg112=reg22*reg33; T reg113=reg28*reg0; T reg114=reg54-reg46; T reg115=reg33*reg2;
    reg47=reg50+reg47; T reg116=reg32*reg8; T reg117=reg8*reg31; T reg118=reg25*reg21; reg25=reg25*reg7;
    reg39=reg39*reg2; reg51=reg51*reg14; T reg119=reg0*reg33; reg33=reg33*reg21; reg30=reg14*reg30;
    reg73=reg75-reg73; reg75=reg49*reg5; reg33=reg111+reg33; reg111=reg5*reg52; T reg120=reg99+reg107;
    reg102=reg100-reg102; reg100=reg5*reg65; reg97=reg97-reg95; reg77=reg77-reg78; reg106=reg61+reg106;
    reg38=reg38-reg92; reg24=reg24+reg108; reg93=reg104-reg93; reg105=reg105+reg101; reg61=reg5*reg44;
    reg29=reg29-reg98; reg30=reg30+reg113; reg90=reg91-reg90; reg119=reg51-reg119; reg96=reg96*reg5;
    reg118=reg118-reg116; reg88=reg87-reg88; reg51=reg47*reg5; reg87=reg71+reg103; reg91=reg5*reg42;
    reg40=reg40+reg83; reg104=reg64*reg5; reg112=reg57+reg112; reg39=reg58-reg39; reg82=reg82+reg81;
    reg94=reg94-reg79; reg86=reg84+reg86; reg57=reg55*reg5; reg80=reg80-reg85; reg76=reg76+reg66;
    reg70=reg70-reg72; reg69=reg69+reg68; reg115=reg110-reg115; reg58=reg74+reg67; reg114=reg114*reg5;
    reg25=reg25+reg109; reg84=reg117+reg89; reg110=ponderation*reg111; reg106=reg106*reg5; reg30=reg30*reg5;
    reg82=reg82*reg5; T reg121=reg84*reg5; reg90=reg5*reg90; T reg122=reg58*reg5; reg39=reg39*reg5;
    reg94=reg94*reg5; reg119=reg119*reg5; reg33=reg33*reg5; T reg123=ponderation*reg75; T reg124=ponderation*reg104;
    T reg125=ponderation*reg57; reg118=reg118*reg5; reg96=ponderation*reg96; reg112=reg112*reg5; reg88=reg5*reg88;
    T reg126=ponderation*reg91; reg73=reg5*reg73; T reg127=reg87*reg5; reg25=reg5*reg25; reg40=reg40*reg5;
    reg115=reg115*reg5; T reg128=reg5*reg120; T reg129=ponderation*reg51; T reg130=ponderation*reg100; reg24=reg24*reg5;
    reg105=reg105*reg5; reg97=reg5*reg97; reg114=ponderation*reg114; reg69=reg69*reg5; reg77=reg5*reg77;
    reg86=reg86*reg5; T reg131=ponderation*reg61; reg29=reg5*reg29; reg102=reg102*reg5; reg93=reg5*reg93;
    reg70=reg70*reg5; reg76=reg76*reg5; reg38=reg38*reg5; reg80=reg80*reg5; T reg132=ponderation*reg122;
    T vec_2=reg132; T tmp_4_2=-reg130; T tmp_0_5=ponderation*reg70; T tmp_1_3=ponderation*reg80; T tmp_5_3=-reg125;
    T tmp_4_0=ponderation*reg73; T tmp_3_5=-reg129; T tmp_0_1=ponderation*reg112; T tmp_0_2=ponderation*reg29; reg29=ponderation*reg127;
    T vec_0=reg29; T tmp_2_4=-reg126; T tmp_5_4=-reg124; T tmp_2_0=ponderation*reg39; T tmp_5_2=ponderation*reg30;
    T tmp_0_0=ponderation*reg86; T tmp_3_3=ponderation*reg24; T tmp_3_4=ponderation*reg76; T tmp_4_3=ponderation*reg105; T tmp_2_1=ponderation*reg115;
    T vec_4=-reg114; T vec_3=-reg96; reg24=ponderation*reg121; T vec_1=reg24; T tmp_4_4=ponderation*reg25;
    T tmp_5_5=ponderation*reg82; T tmp_3_0=ponderation*reg102; T tmp_4_1=ponderation*reg88; T tmp_1_4=ponderation*reg118; T tmp_2_3=-reg123;
    T tmp_5_1=ponderation*reg119; T tmp_5_0=ponderation*reg90; T tmp_4_5=-reg110; T tmp_0_4=reg94*ponderation; T tmp_3_2=-reg131;
    T tmp_3_1=ponderation*reg93; T tmp_0_3=ponderation*reg77; T tmp_1_2=ponderation*reg97; T tmp_2_2=ponderation*reg69; T tmp_1_5=ponderation*reg38;
    T tmp_1_1=ponderation*reg33; reg25=ponderation*reg128; T vec_5=reg25; T tmp_2_5=ponderation*reg40; T tmp_1_0=ponderation*reg106;
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg1*reg0; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg3*reg2; T reg6=reg3*reg1; reg2=reg2*reg4; reg1=reg1*reg4;
    T reg7=reg3*reg2; reg2=reg2*reg4; T reg8=reg1*reg3; T reg9=reg3*reg6; reg1=reg1*reg4;
    T reg10=reg3*reg5; reg8=reg8+reg9; reg1=reg1-reg9; reg6=reg6*reg4; reg7=reg7+reg10;
    reg5=reg5*reg4; reg2=reg2-reg10; reg8=reg3*reg8; reg1=reg1*reg4; T reg11=reg3*reg7;
    T reg12=reg2*reg4; reg5=reg10+reg5; reg10=reg6+reg9; T reg13=reg3*reg5; reg11=reg12-reg11;
    reg8=reg1-reg8; reg10=reg10*reg3; reg10=reg8-reg10; reg13=reg11-reg13; reg10=reg10/reg13;
    reg7=reg7/reg13; reg2=reg2/reg13; reg1=reg2*reg10; reg8=reg7*reg10; reg11=(*f.m).alpha*reg7;
    reg12=reg1*reg2; T reg14=reg7*reg8; reg13=reg5/reg13; reg5=(*f.m).alpha*reg2; reg14=reg12-reg14;
    reg13=(*f.m).alpha*reg13; reg5=reg11+reg5; reg1=reg1/reg14; reg11=reg3*reg0; reg14=reg8/reg14;
    reg13=reg5+reg13; reg0=reg4*reg0; reg5=elem.pos(2)[1]-elem.pos(0)[1]; reg8=elem.pos(1)[1]-elem.pos(0)[1]; reg12=pow(reg4,2);
    T reg15=pow(reg3,2); reg14=reg13*reg14; reg1=reg13*reg1; reg13=elem.pos(1)[0]-elem.pos(0)[0]; T reg16=elem.pos(2)[0]-elem.pos(0)[0];
    reg4=reg0*reg4; reg3=reg3*reg11; T reg17=1-(*f.m).resolution; reg3=reg4-reg3; reg4=reg16*reg8;
    T reg18=reg13*reg5; reg15=reg12-reg15; reg14=reg1-reg14; reg14=reg17*reg14; reg15=reg15/reg3;
    reg4=reg18-reg4; reg11=reg11/reg3; reg3=reg0/reg3; reg0=(*f.m).alpha*(*f.m).resolution; reg5=reg5/reg4;
    reg0=reg14+reg0; reg2=reg17*reg2; reg1=(*f.m).resolution*reg3; reg12=(*f.m).resolution*reg15; reg10=reg17*reg10;
    reg7=reg17*reg7; reg14=(*f.m).resolution*reg11; reg16=reg16/reg4; reg8=reg8/reg4; reg13=reg13/reg4;
    reg10=reg12+reg10; reg12=0.5*reg8; reg17=reg8-reg5; reg18=0.5*reg13; reg0=(*f.m).deltaT*reg0;
    T reg19=reg16-reg13; T reg20=0.5*reg5; reg14=reg7+reg14; reg1=reg2+reg1; reg2=0.5*reg19;
    reg7=reg10*reg18; T reg21=reg14*reg0; T reg22=reg1*reg0; T reg23=reg10*reg20; T reg24=0.5*reg16;
    T reg25=reg12*reg10; T reg26=0.5*reg17; T reg27=reg2*reg10; T reg28=reg22+reg21; T reg29=reg14*reg16;
    T reg30=reg14*reg13; T reg31=2*reg25; T reg32=reg10*reg26; reg7=2*reg7; T reg33=1-var_inter[0];
    T reg34=reg1*reg8; T reg35=reg10*reg24; T reg36=reg1*reg13; reg23=2*reg23; T reg37=reg28*reg8;
    T reg38=reg30*reg8; T reg39=reg7*reg24; T reg40=reg5*reg1; T reg41=reg31*reg20; T reg42=reg36*reg16;
    T reg43=reg31*reg18; reg32=2*reg32; T reg44=var_inter[0]*elem.f_vol_e[1]; T reg45=reg19*reg14; T reg46=reg5*reg14;
    reg33=reg33-var_inter[1]; T reg47=2*reg35; T reg48=reg17*reg1; T reg49=var_inter[1]*elem.f_vol_e[0]; reg27=2*reg27;
    T reg50=reg5*reg29; T reg51=reg24*reg23; T reg52=reg19*reg1; T reg53=reg1*reg16; T reg54=reg14*reg8;
    T reg55=reg34*reg5; T reg56=reg28*reg16; T reg57=reg53*reg19; T reg58=reg33*elem.f_vol_e[1]; T reg59=reg31*reg12;
    T reg60=reg47*reg26; T reg61=reg19*reg46; T reg62=reg17*reg28; T reg63=reg2*reg7; T reg64=reg26*reg23;
    T reg65=reg33*elem.f_vol_e[0]; reg38=reg43+reg38; T reg66=reg36*reg13; T reg67=reg31*reg24; T reg68=reg54*reg16;
    reg39=reg55+reg39; T reg69=reg30*reg5; T reg70=reg20*reg7; T reg71=reg31*reg26; T reg72=reg53*reg16;
    T reg73=reg20*reg23; T reg74=reg2*reg31; T reg75=reg7*reg26; T reg76=reg54*reg19; reg30=reg17*reg30;
    T reg77=reg56-reg44; T reg78=reg5*reg40; T reg79=reg47*reg24; T reg80=var_inter[0]*elem.f_vol_e[0]; T reg81=reg2*reg27;
    T reg82=reg17*reg45; reg36=reg36*reg19; T reg83=reg17*reg48; T reg84=reg2*reg23; reg42=reg41+reg42;
    T reg85=reg7*reg18; T reg86=reg19*reg28; reg51=reg50+reg51; T reg87=reg52*reg19; T reg88=reg2*reg32;
    T reg89=reg2*reg47; T reg90=reg32*reg26; T reg91=reg5*reg28; T reg92=reg17*reg29; T reg93=var_inter[1]*elem.f_vol_e[1];
    T reg94=reg37-reg49; T reg95=reg34*reg8; T reg96=reg34*reg17; T reg97=reg17*reg40; T reg98=reg28*reg13;
    reg84=reg84-reg92; T reg99=reg65+reg62; reg61=reg61-reg60; reg78=reg78+reg79; T reg100=reg86+reg58;
    reg63=reg63-reg96; T reg101=reg91+reg80; T reg102=reg93+reg98; reg36=reg36-reg71; T reg103=reg4*reg39;
    reg83=reg81+reg83; reg69=reg69+reg67; reg66=reg66+reg59; reg94=reg94*reg4; reg97=reg97-reg89;
    reg90=reg87+reg90; reg70=reg70+reg68; reg81=reg51*reg4; reg73=reg73+reg72; reg75=reg75-reg76;
    reg64=reg64-reg57; reg77=reg77*reg4; reg30=reg30-reg74; reg88=reg82+reg88; reg85=reg85+reg95;
    reg82=reg42*reg4; reg87=reg4*reg38; T reg104=ponderation*reg81; T reg105=ponderation*reg87; T reg106=ponderation*reg82;
    reg97=reg4*reg97; reg63=reg63*reg4; reg36=reg36*reg4; reg69=reg69*reg4; T reg107=reg100*reg4;
    reg78=reg78*reg4; T reg108=reg99*reg4; reg88=reg88*reg4; reg30=reg30*reg4; reg75=reg75*reg4;
    reg73=reg73*reg4; reg90=reg90*reg4; reg77=ponderation*reg77; reg94=ponderation*reg94; reg83=reg83*reg4;
    T reg109=reg101*reg4; T reg110=ponderation*reg103; reg66=reg66*reg4; T reg111=reg4*reg102; reg85=reg4*reg85;
    reg64=reg64*reg4; reg61=reg4*reg61; reg70=reg70*reg4; reg84=reg4*reg84; T tmp_2_5=ponderation*reg69;
    T tmp_1_2=ponderation*reg61; reg61=ponderation*reg107; T vec_1=reg61; T tmp_0_3=ponderation*reg84; T tmp_1_1=ponderation*reg90;
    T tmp_1_5=ponderation*reg36; reg36=ponderation*reg108; T vec_0=reg36; reg69=ponderation*reg111; T vec_5=reg69;
    T tmp_0_2=ponderation*reg97; T tmp_0_1=ponderation*reg88; T tmp_3_5=-reg106; T tmp_0_4=reg63*ponderation; T tmp_0_5=ponderation*reg30;
    T tmp_1_4=ponderation*reg75; T tmp_4_5=-reg105; T tmp_2_2=ponderation*reg78; T tmp_5_5=ponderation*reg66; T tmp_4_4=ponderation*reg85;
    T vec_3=-reg77; T vec_4=-reg94; T tmp_0_0=ponderation*reg83; reg30=ponderation*reg109; T vec_2=reg30;
    T tmp_1_3=ponderation*reg64; T tmp_3_4=ponderation*reg70; T tmp_2_4=-reg110; T tmp_3_3=ponderation*reg73; T tmp_2_3=-reg104;
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=reg1*reg0;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg4*reg3; T reg6=reg4*reg1; reg3=reg3*reg2; reg1=reg1*reg2;
    T reg7=reg1*reg2; T reg8=reg3*reg2; T reg9=reg4*reg6; reg1=reg1*reg4; reg3=reg4*reg3;
    T reg10=reg4*reg5; reg8=reg8-reg10; reg1=reg1+reg9; reg3=reg3+reg10; reg7=reg7-reg9;
    reg6=reg6*reg2; reg5=reg5*reg2; T reg11=reg6+reg9; T reg12=reg8*reg2; reg5=reg10+reg5;
    reg1=reg4*reg1; reg10=reg4*reg0; T reg13=reg4*reg3; reg0=reg2*reg0; reg7=reg7*reg2;
    T reg14=reg4*reg5; T reg15=pow(reg2,2); reg1=reg7-reg1; reg7=reg4*reg10; reg2=reg0*reg2;
    T reg16=elem.pos(2)[0]-elem.pos(0)[0]; T reg17=elem.pos(1)[0]-elem.pos(0)[0]; T reg18=pow(reg4,2); T reg19=elem.pos(1)[1]-elem.pos(0)[1]; T reg20=elem.pos(2)[1]-elem.pos(0)[1];
    reg13=reg12-reg13; reg4=reg11*reg4; reg11=reg16*reg19; reg12=reg17*reg20; reg7=reg2-reg7;
    reg18=reg15-reg18; reg4=reg1-reg4; reg14=reg13-reg14; reg11=reg12-reg11; reg18=reg18/reg7;
    reg1=1-(*f.m).resolution; reg4=reg4/reg14; reg10=reg10/reg7; reg20=reg20/reg11; reg8=reg8/reg14;
    reg3=reg3/reg14; reg2=(*f.m).resolution*reg18; reg7=reg0/reg7; reg17=reg17/reg11; reg19=reg19/reg11;
    reg16=reg16/reg11; reg0=reg1*reg4; reg12=reg16-reg17; reg13=0.5*reg20; reg0=reg2+reg0;
    reg2=reg1*reg3; reg15=0.5*reg16; T reg21=0.5*reg17; T reg22=(*f.m).resolution*reg10; T reg23=(*f.m).resolution*reg7;
    T reg24=0.5*reg19; T reg25=reg19-reg20; T reg26=reg1*reg8; T reg27=0.5*reg12; T reg28=reg0*reg15;
    reg23=reg26+reg23; reg22=reg2+reg22; reg2=reg24*reg0; reg26=reg0*reg21; T reg29=0.5*reg25;
    T reg30=reg0*reg13; T reg31=reg23*reg16; T reg32=reg20*reg22; T reg33=reg0*reg29; reg26=2*reg26;
    T reg34=reg20*reg23; T reg35=reg22*reg19; reg30=2*reg30; T reg36=reg22*reg16; T reg37=2*reg28;
    T reg38=reg23*reg19; T reg39=reg23*reg17; T reg40=2*reg2; T reg41=reg27*reg0; T reg42=reg22*reg17;
    T reg43=reg38*reg20; T reg44=reg31*reg17; T reg45=reg24*reg30; T reg46=reg26*reg15; T reg47=reg20*reg36;
    T reg48=reg15*reg30; T reg49=reg25*reg23; reg41=2*reg41; T reg50=reg12*reg23; T reg51=reg40*reg21;
    T reg52=reg35*reg17; T reg53=reg42*reg19; T reg54=reg39*reg16; reg33=2*reg33; T reg55=reg37*reg13;
    T reg56=reg16*reg32; T reg57=reg12*reg22; T reg58=reg24*reg26; T reg59=reg40*reg13; T reg60=reg25*reg22;
    T reg61=reg19*reg34; T reg62=reg37*reg21; T reg63=reg33*reg21; T reg64=reg42*reg20; T reg65=reg33*reg15;
    T reg66=reg40*reg15; reg46=reg43+reg46; T reg67=reg57*reg19; T reg68=reg57*reg20; T reg69=reg39*reg17;
    T reg70=reg40*reg24; T reg71=reg27*reg33; T reg72=reg29*reg30; T reg73=reg27*reg26; T reg74=reg50*reg16;
    T reg75=reg33*reg13; T reg76=reg37*reg29; T reg77=reg12*reg32; reg56=reg55+reg56; T reg78=reg27*reg40;
    reg42=reg25*reg42; T reg79=reg20*reg34; reg53=reg51+reg53; T reg80=reg60*reg17; T reg81=reg41*reg24;
    T reg82=reg37*reg15; T reg83=reg27*reg41; T reg84=reg25*reg49; T reg85=reg26*reg21; T reg86=reg41*reg15;
    T reg87=reg35*reg16; T reg88=reg13*reg26; T reg89=reg60*reg16; reg45=reg44+reg45; T reg90=reg41*reg13;
    T reg91=reg31*reg12; reg60=reg60*reg12; T reg92=reg41*reg29; reg48=reg47+reg48; T reg93=reg49*reg20;
    T reg94=reg25*reg36; T reg95=reg31*reg16; reg39=reg39*reg12; T reg96=reg13*reg30; reg41=reg41*reg21;
    reg49=reg49*reg19; T reg97=reg27*reg30; T reg98=reg35*reg12; T reg99=reg37*reg24; reg32=reg17*reg32;
    reg26=reg26*reg29; T reg100=reg40*reg29; T reg101=reg24*reg33; T reg102=reg50*reg17; reg34=reg25*reg34;
    T reg103=reg38*reg25; T reg104=reg38*reg19; reg54=reg59+reg54; reg57=reg25*reg57; reg30=reg21*reg30;
    reg33=reg33*reg29; reg58=reg52+reg58; reg50=reg50*reg12; T reg105=reg27*reg37; T reg106=reg19*reg36;
    reg61=reg62+reg61; reg79=reg79+reg82; reg26=reg26-reg98; reg81=reg80-reg81; reg85=reg85+reg104;
    reg30=reg30+reg106; reg71=reg57+reg71; reg57=reg11*reg61; reg74=reg75-reg74; reg97=reg97-reg94;
    reg75=reg11*reg53; reg77=reg77-reg76; reg65=reg68-reg65; reg92=reg60+reg92; reg42=reg42-reg78;
    reg73=reg73-reg103; reg84=reg83+reg84; reg60=reg54*reg11; reg34=reg34-reg105; reg33=reg50+reg33;
    reg72=reg72-reg91; reg67=reg63-reg67; reg50=reg11*reg46; reg63=reg48*reg11; reg64=reg64+reg66;
    reg68=reg58*reg11; reg69=reg69+reg70; reg96=reg96+reg95; reg39=reg39-reg100; reg101=reg102-reg101;
    reg49=reg41-reg49; reg89=reg90-reg89; reg88=reg88+reg87; reg32=reg32+reg99; reg41=reg45*reg11;
    reg80=reg11*reg56; reg86=reg93-reg86; reg39=reg39*reg11; reg73=reg73*reg11; reg34=reg11*reg34;
    reg83=ponderation*reg60; reg89=reg89*reg11; reg81=reg11*reg81; reg67=reg11*reg67; reg90=ponderation*reg63;
    reg74=reg11*reg74; reg49=reg11*reg49; reg93=ponderation*reg75; reg64=reg64*reg11; reg30=reg30*reg11;
    reg33=reg33*reg11; reg84=reg84*reg11; reg102=ponderation*reg41; reg96=reg96*reg11; reg72=reg72*reg11;
    reg86=reg86*reg11; T reg107=ponderation*reg50; T reg108=ponderation*reg68; reg69=reg69*reg11; T reg109=ponderation*reg80;
    reg32=reg32*reg11; reg101=reg101*reg11; reg88=reg88*reg11; reg71=reg71*reg11; reg79=reg79*reg11;
    reg97=reg11*reg97; reg85=reg11*reg85; reg26=reg26*reg11; T reg110=ponderation*reg57; reg42=reg42*reg11;
    reg65=reg65*reg11; reg92=reg92*reg11; reg77=reg11*reg77; T tmp_0_1=ponderation*reg71; T tmp_1_0=ponderation*reg92;
    T tmp_2_2=ponderation*reg79; T tmp_1_5=ponderation*reg39; T tmp_0_5=ponderation*reg42; T tmp_3_5=-reg83; T tmp_4_2=-reg110;
    T tmp_3_2=-reg109; T tmp_5_3=-reg102; T tmp_5_1=ponderation*reg101; T tmp_4_0=ponderation*reg49; T tmp_0_2=ponderation*reg34;
    T tmp_1_4=ponderation*reg26; T tmp_3_0=ponderation*reg89; T tmp_0_3=ponderation*reg97; T tmp_3_1=ponderation*reg74; T tmp_1_2=ponderation*reg77;
    T tmp_1_1=ponderation*reg33; T tmp_2_1=ponderation*reg65; T tmp_5_5=ponderation*reg69; T tmp_4_5=-reg93; T tmp_4_4=ponderation*reg85;
    T tmp_3_4=ponderation*reg88; T tmp_5_0=ponderation*reg81; T tmp_3_3=ponderation*reg96; T tmp_5_2=ponderation*reg32; T tmp_0_4=reg73*ponderation;
    T tmp_2_3=-reg90; T tmp_5_4=-reg108; T tmp_2_4=-reg107; T tmp_4_1=ponderation*reg67; T tmp_1_3=ponderation*reg72;
    T tmp_2_0=ponderation*reg86; T tmp_0_0=ponderation*reg84; T tmp_2_5=ponderation*reg64; T tmp_4_3=ponderation*reg30;
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=reg1*reg0;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg4*reg1; T reg6=reg4*reg3; reg1=reg1*reg2; reg3=reg3*reg2;
    T reg7=reg4*reg3; reg3=reg3*reg2; T reg8=reg4*reg5; T reg9=reg1*reg4; T reg10=reg4*reg6;
    reg1=reg1*reg2; reg1=reg1-reg8; reg3=reg3-reg10; reg6=reg6*reg2; reg5=reg5*reg2;
    reg9=reg9+reg8; reg7=reg7+reg10; T reg11=reg3*reg2; T reg12=reg2*reg0; T reg13=reg4*reg7;
    reg9=reg4*reg9; reg6=reg10+reg6; reg0=reg4*reg0; reg10=reg5+reg8; reg1=reg1*reg2;
    T reg14=reg4*reg6; reg9=reg1-reg9; reg1=pow(reg2,2); T reg15=reg4*reg0; reg2=reg12*reg2;
    T reg16=elem.pos(2)[0]-elem.pos(0)[0]; T reg17=elem.pos(1)[0]-elem.pos(0)[0]; T reg18=pow(reg4,2); T reg19=elem.pos(2)[1]-elem.pos(0)[1]; T reg20=elem.pos(1)[1]-elem.pos(0)[1];
    reg4=reg10*reg4; reg13=reg11-reg13; reg18=reg1-reg18; reg4=reg9-reg4; reg15=reg2-reg15;
    reg1=reg16*reg20; reg14=reg13-reg14; reg2=reg17*reg19; reg9=1-(*f.m).resolution; reg18=reg18/reg15;
    reg1=reg2-reg1; reg4=reg4/reg14; reg0=reg0/reg15; reg17=reg17/reg1; reg15=reg12/reg15;
    reg20=reg20/reg1; reg16=reg16/reg1; reg2=reg9*reg4; reg7=reg7/reg14; reg3=reg3/reg14;
    reg19=reg19/reg1; reg10=(*f.m).resolution*reg18; reg11=reg16-reg17; reg12=0.5*reg17; reg2=reg10+reg2;
    reg10=reg20-reg19; reg13=(*f.m).resolution*reg0; T reg21=reg9*reg3; T reg22=reg9*reg7; T reg23=(*f.m).resolution*reg15;
    T reg24=0.5*reg20; T reg25=0.5*reg19; T reg26=0.5*reg10; T reg27=reg2*reg25; reg13=reg22+reg13;
    reg23=reg21+reg23; reg21=0.5*reg11; reg22=reg2*reg12; T reg28=0.5*reg16; T reg29=reg24*reg2;
    T reg30=reg2*reg26; reg22=2*reg22; T reg31=reg23*reg17; T reg32=reg2*reg28; reg27=2*reg27;
    T reg33=reg13*reg16; T reg34=reg21*reg2; T reg35=reg13*reg17; T reg36=2*reg29; T reg37=reg23*reg20;
    T reg38=reg19*reg13; T reg39=reg23*reg16; T reg40=reg11*reg13; reg30=2*reg30; T reg41=reg36*reg12;
    T reg42=reg35*reg20; T reg43=reg13*reg20; T reg44=reg36*reg25; T reg45=reg31*reg16; T reg46=reg19*reg23;
    T reg47=reg10*reg23; reg34=2*reg34; T reg48=reg22*reg28; T reg49=reg37*reg19; T reg50=reg19*reg33;
    T reg51=2*reg32; T reg52=reg28*reg27; T reg53=reg11*reg23; T reg54=reg10*reg40; T reg55=reg25*reg22;
    T reg56=reg43*reg16; T reg57=reg21*reg27; reg45=reg44+reg45; T reg58=reg39*reg11; T reg59=reg10*reg46;
    T reg60=reg36*reg24; T reg61=reg22*reg12; T reg62=reg10*reg47; T reg63=reg21*reg34; T reg64=reg10*reg33;
    T reg65=reg51*reg28; T reg66=reg19*reg46; T reg67=reg10*reg35; reg42=reg41+reg42; T reg68=reg21*reg36;
    T reg69=reg37*reg10; T reg70=reg53*reg11; T reg71=reg30*reg26; reg52=reg50+reg52; T reg72=reg37*reg20;
    T reg73=reg21*reg22; reg48=reg49+reg48; T reg74=reg21*reg30; T reg75=reg39*reg16; T reg76=reg25*reg27;
    T reg77=reg21*reg51; T reg78=reg22*reg26; T reg79=reg43*reg11; T reg80=reg26*reg27; T reg81=reg31*reg11;
    T reg82=reg51*reg26; T reg83=reg11*reg38; reg35=reg35*reg19; T reg84=reg36*reg28; T reg85=reg36*reg26;
    reg31=reg31*reg17; reg57=reg57-reg64; reg83=reg83-reg82; reg73=reg73-reg69; T reg86=reg1*reg42;
    reg59=reg59-reg77; reg80=reg80-reg58; T reg87=reg45*reg1; T reg88=reg52*reg1; T reg89=reg1*reg48;
    reg71=reg70+reg71; reg35=reg35+reg84; reg76=reg76+reg75; reg31=reg31+reg60; reg81=reg81-reg85;
    reg66=reg66+reg65; reg55=reg55+reg56; reg74=reg54+reg74; reg78=reg78-reg79; reg67=reg67-reg68;
    reg61=reg61+reg72; reg62=reg63+reg62; reg71=reg71*reg1; reg54=ponderation*reg86; reg63=ponderation*reg88;
    reg70=ponderation*reg87; reg76=reg76*reg1; reg78=reg78*reg1; reg67=reg67*reg1; reg59=reg1*reg59;
    reg62=reg62*reg1; reg74=reg74*reg1; reg80=reg80*reg1; reg81=reg81*reg1; reg57=reg1*reg57;
    reg61=reg1*reg61; reg83=reg1*reg83; T reg90=ponderation*reg89; reg35=reg35*reg1; reg66=reg66*reg1;
    reg73=reg73*reg1; reg31=reg31*reg1; reg55=reg55*reg1; T tmp_1_5=ponderation*reg81; T tmp_2_5=ponderation*reg35;
    T tmp_0_1=ponderation*reg74; T tmp_3_5=-reg70; T tmp_0_2=ponderation*reg59; T tmp_0_5=ponderation*reg67; T tmp_1_4=ponderation*reg78;
    T tmp_2_2=ponderation*reg66; T tmp_2_3=-reg63; T tmp_0_4=reg73*ponderation; T tmp_4_5=-reg54; T tmp_0_0=ponderation*reg62;
    T tmp_1_3=ponderation*reg80; T tmp_0_3=ponderation*reg57; T tmp_1_2=ponderation*reg83; T tmp_2_4=-reg90; T tmp_1_1=ponderation*reg71;
    T tmp_3_3=ponderation*reg76; T tmp_3_4=ponderation*reg55; T tmp_5_5=ponderation*reg31; T tmp_4_4=ponderation*reg61;
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=reg1*reg0;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg2*reg3; T reg6=reg2*reg1; reg3=reg3*reg4; reg1=reg1*reg4;
    T reg7=reg2*reg3; T reg8=reg2*reg5; T reg9=reg1*reg4; T reg10=reg2*reg6; reg3=reg3*reg4;
    reg1=reg1*reg2; reg5=reg5*reg4; reg3=reg3-reg8; reg1=reg1+reg10; reg7=reg7+reg8;
    reg9=reg9-reg10; reg6=reg6*reg4; reg9=reg9*reg4; reg5=reg8+reg5; reg8=reg3*reg4;
    reg1=reg2*reg1; T reg11=reg6+reg10; T reg12=reg2*reg7; reg1=reg9-reg1; reg9=reg2*reg5;
    reg11=reg11*reg2; reg12=reg8-reg12; reg11=reg1-reg11; reg9=reg12-reg9; reg7=reg7/reg9;
    reg3=reg3/reg9; reg11=reg11/reg9; reg1=reg7*reg11; reg8=reg3*reg11; reg9=reg5/reg9;
    reg5=(*f.m).alpha*reg3; reg12=reg7*reg1; T reg13=reg8*reg3; T reg14=(*f.m).alpha*reg7; reg9=(*f.m).alpha*reg9;
    reg12=reg13-reg12; reg5=reg14+reg5; reg9=reg5+reg9; reg5=reg2*reg0; reg8=reg8/reg12;
    reg12=reg1/reg12; reg0=reg4*reg0; reg12=reg9*reg12; reg8=reg9*reg8; reg1=reg0*reg4;
    reg9=reg2*reg5; reg9=reg1-reg9; reg1=1-(*f.m).resolution; reg12=reg8-reg12; reg12=reg1*reg12;
    reg5=reg5/reg9; reg8=(*f.m).alpha*(*f.m).resolution; reg0=reg0/reg9; reg13=elem.pos(2)[1]-elem.pos(0)[1]; reg14=elem.pos(1)[1]-elem.pos(0)[1];
    T reg15=(*f.m).resolution*reg0; T reg16=(*f.m).resolution*reg5; T reg17=elem.pos(1)[0]-elem.pos(0)[0]; T reg18=elem.pos(2)[0]-elem.pos(0)[0]; reg8=reg12+reg8;
    reg3=reg1*reg3; reg7=reg1*reg7; reg16=reg7+reg16; reg15=reg3+reg15; reg8=(*f.m).deltaT*reg8;
    reg3=reg18*reg14; reg7=reg17*reg13; reg12=reg15*reg8; T reg19=reg16*reg8; reg3=reg7-reg3;
    reg7=reg12+reg19; T reg20=1-var_inter[0]; reg18=reg18/reg3; reg13=reg13/reg3; reg14=reg14/reg3;
    reg17=reg17/reg3; T reg21=var_inter[0]*elem.f_vol_e[1]; T reg22=reg14-reg13; T reg23=reg7*reg18; reg20=reg20-var_inter[1];
    T reg24=var_inter[1]*elem.f_vol_e[0]; T reg25=reg7*reg14; T reg26=reg18-reg17; T reg27=reg20*elem.f_vol_e[0]; T reg28=reg23-reg21;
    T reg29=var_inter[1]*elem.f_vol_e[1]; T reg30=reg13*reg7; T reg31=var_inter[0]*elem.f_vol_e[0]; T reg32=reg26*reg7; T reg33=reg20*elem.f_vol_e[1];
    T reg34=reg7*reg17; T reg35=reg22*reg7; T reg36=reg25-reg24; T reg37=reg29+reg34; reg36=reg36*reg3;
    T reg38=reg27+reg35; T reg39=reg30+reg31; T reg40=reg32+reg33; reg28=reg28*reg3; T reg41=reg3*reg37;
    T reg42=reg39*reg3; T reg43=reg38*reg3; reg28=ponderation*reg28; T reg44=reg40*reg3; reg36=ponderation*reg36;
    T reg45=ponderation*reg42; T vec_2=reg45; T reg46=ponderation*reg44; T vec_1=reg46; T vec_3=-reg28;
    reg28=ponderation*reg43; T vec_0=reg28; T vec_4=-reg36; reg36=ponderation*reg41; T vec_5=reg36;
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices ) { 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=reg1*reg0; T reg5=reg2*reg4; T reg6=reg2*reg1; reg4=reg4*reg3; reg1=reg1*reg3;
    T reg7=reg2*reg5; T reg8=reg1*reg3; T reg9=reg2*reg6; T reg10=reg2*reg4; reg4=reg4*reg3;
    reg1=reg1*reg2; reg5=reg5*reg3; reg6=reg6*reg3; reg4=reg4-reg7; reg1=reg1+reg9;
    reg8=reg8-reg9; reg10=reg10+reg7; reg8=reg8*reg3; T reg11=reg2*reg10; T reg12=reg6+reg9;
    reg1=reg2*reg1; T reg13=reg4*reg3; reg5=reg7+reg5; reg1=reg8-reg1; reg7=reg2*reg5;
    reg11=reg13-reg11; reg12=reg12*reg2; reg12=reg1-reg12; reg7=reg11-reg7; reg12=reg12/reg7;
    reg4=reg4/reg7; reg10=reg10/reg7; reg1=reg10*reg12; reg8=reg4*reg12; reg7=reg5/reg7;
    reg5=(*f.m).alpha*reg10; reg11=reg10*reg1; reg13=(*f.m).alpha*reg4; T reg14=reg8*reg4; reg11=reg14-reg11;
    reg14=elem.pos(1)[1]-elem.pos(0)[1]; T reg15=elem.pos(1)[0]-elem.pos(0)[0]; T reg16=elem.pos(2)[0]-elem.pos(0)[0]; reg13=reg5+reg13; reg5=elem.pos(2)[1]-elem.pos(0)[1];
    reg7=(*f.m).alpha*reg7; reg1=reg1/reg11; reg11=reg8/reg11; reg8=reg16*reg14; T reg17=reg15*reg5;
    reg7=reg13+reg7; reg8=reg17-reg8; reg13=reg2*reg0; reg1=reg7*reg1; reg11=reg7*reg11;
    reg0=reg3*reg0; reg5=reg5/reg8; reg1=reg11-reg1; reg7=pow(reg2,2); reg11=vectors[0][indices[2]+0]-vectors[0][indices[0]+0];
    reg2=reg2*reg13; reg17=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; T reg18=1-(*f.m).resolution; T reg19=pow(reg3,2); reg15=reg15/reg8;
    reg14=reg14/reg8; reg16=reg16/reg8; T reg20=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg3=reg0*reg3; T reg21=vectors[0][indices[1]+1]-vectors[0][indices[0]+1];
    reg1=reg18*reg1; reg2=reg3-reg2; reg3=reg20*reg14; T reg22=reg5*reg21; T reg23=reg16*reg17;
    T reg24=reg11*reg15; T reg25=(*f.m).alpha*(*f.m).resolution; reg7=reg19-reg7; reg25=reg1+reg25; reg23=reg24-reg23;
    reg17=reg5*reg17; reg11=reg11*reg14; reg7=reg7/reg2; reg0=reg0/reg2; reg3=reg22-reg3;
    reg21=reg21*reg16; reg20=reg20*reg15; reg2=reg13/reg2; reg1=(*f.m).resolution*reg0; reg13=(*f.m).resolution*reg2;
    reg25=(*f.m).deltaT*reg25; reg19=(*f.m).resolution*reg7; reg4=reg18*reg4; reg21=reg20-reg21; reg23=reg3+reg23;
    reg11=reg17-reg11; reg10=reg18*reg10; reg12=reg18*reg12; reg13=reg10+reg13; reg1=reg4+reg1;
    reg23=0.5*reg23; reg11=reg11-reg25; reg21=reg21-reg25; reg12=reg19+reg12; reg3=reg16-reg15;
    reg23=reg12*reg23; reg4=reg14-reg5; reg10=reg21*reg1; reg17=reg11*reg13; reg11=reg11*reg1;
    reg21=reg21*reg13; reg18=0.5*reg15; reg23=2*reg23; reg19=0.5*reg16; reg10=reg17+reg10;
    reg17=1-var_inter[0]; reg21=reg11+reg21; reg11=0.5*reg4; reg20=0.5*reg14; reg22=0.5*reg5;
    reg24=0.5*reg3; T reg26=reg23*reg19; T reg27=reg24*reg23; T reg28=reg20*reg23; T reg29=reg23*reg18;
    T reg30=reg3*reg10; T reg31=reg5*reg21; T reg32=reg14*reg21; T reg33=reg23*reg22; T reg34=reg16*reg10;
    T reg35=reg4*reg21; T reg36=reg23*reg11; reg17=reg17-var_inter[1]; T reg37=reg15*reg10; T reg38=var_inter[1]*elem.f_vol_e[1];
    reg37=reg37-reg28; T reg39=var_inter[0]*elem.f_vol_e[1]; T reg40=reg17*elem.f_vol_e[0]; reg35=reg27+reg35; reg27=var_inter[1]*elem.f_vol_e[0];
    reg29=reg29-reg32; reg33=reg33-reg34; T reg41=reg17*elem.f_vol_e[1]; reg36=reg30+reg36; reg31=reg31-reg26;
    reg30=var_inter[0]*elem.f_vol_e[0]; reg29=reg29-reg27; reg35=reg35-reg40; reg31=reg31-reg30; reg37=reg37-reg38;
    reg33=reg33-reg39; reg36=reg36-reg41; reg37=reg37*reg8; reg33=reg8*reg33; reg35=reg8*reg35;
    reg36=reg8*reg36; reg29=reg29*reg8; reg31=reg31*reg8; T vec_2=ponderation*reg31; T vec_4=ponderation*reg29;
    T vec_0=ponderation*reg35; T vec_5=ponderation*reg37; T vec_3=ponderation*reg33; T vec_1=ponderation*reg36;
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
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
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
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
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
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
#ifndef ELASTICITY_ISOTROPY_STAT_QSTAT_2_NUM_0_CARAC_H
#define ELASTICITY_ISOTROPY_STAT_QSTAT_2_NUM_0_CARAC_H
#ifndef STRUCT_ELASTICITY_ISOTROPY_STAT_QSTAT
#define STRUCT_ELASTICITY_ISOTROPY_STAT_QSTAT
struct elasticity_isotropy_stat_Qstat {
  static const char *name() { return "elasticity_isotropy_stat_Qstat"; }
};
#endif // ELASTICITY_ISOTROPY_STAT_QSTAT

template<class P_T>
class CaracFormulation<elasticity_isotropy_stat_Qstat,2,P_T>  {
public:
  typedef P_T T;
  static const char *name() { return "elasticity_isotropy_stat_Qstat"; }
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
    T reg0=vecs[1][indice+0]-vecs[0][indice+0]; T reg1=vecs[1][indice+1]-vecs[0][indice+1]; reg0=abs(reg0); reg1=abs(reg1); return max(reg1,reg0);
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
#endif // ELASTICITY_ISOTROPY_STAT_QSTAT_2_NUM_0_CARAC_H
extern "C" void apply_on_elements_after_solve_0_elasticity_isotropy_stat_Qstat_Quad_0( double * );
extern "C" void apply_on_elements_after_solve_1_elasticity_isotropy_stat_Qstat_Quad_1( double * );
extern "C" void apply_on_elements_after_solve_2_elasticity_isotropy_stat_Qstat_Quad_2( double * );
extern "C" void apply_on_elements_after_solve_3_elasticity_isotropy_stat_Qstat_Quad_3( double * );
extern "C" void apply_on_elements_after_solve_4_elasticity_isotropy_stat_Qstat_Quad_4( double * );
extern "C" void apply_on_elements_after_solve_5_elasticity_isotropy_stat_Qstat_Quad_5( double * );
extern "C" void apply_on_elements_after_solve_6_elasticity_isotropy_stat_Qstat_Quad_6( double * );
extern "C" void apply_on_elements_after_solve_7_elasticity_isotropy_stat_Qstat_Quad_7( double * );
extern "C" void apply_on_elements_after_solve_8_elasticity_isotropy_stat_Qstat_Quad_8( double * );
extern "C" void apply_on_elements_after_solve_9_elasticity_isotropy_stat_Qstat_Quad_9( double * );
extern "C" void apply_on_elements_after_solve_10_elasticity_isotropy_stat_Qstat_Quad_10( double * );
extern "C" void apply_on_elements_after_solve_11_elasticity_isotropy_stat_Qstat_Quad_11( double * );
extern "C" void apply_on_elements_after_solve_12_elasticity_isotropy_stat_Qstat_Quad_12( double * );
extern "C" void apply_on_elements_after_solve_13_elasticity_isotropy_stat_Qstat_Quad_13( double * );
extern "C" void apply_on_elements_after_solve_14_elasticity_isotropy_stat_Qstat_Quad_14( double * );
class Quad;
template<unsigned A,class B,class C> class Node;
template<class A,class B,class C,class D,unsigned E> class Element;

// Carac for ...
template<class P_T_pos,class P_ND,class TED,unsigned nim,class TM,class T>
class CaracFormulationForElement<elasticity_isotropy_stat_Qstat,Element<Quad,DefaultBehavior,Node<2,P_T_pos,P_ND>,TED,nim>,TM,T> {
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg2*reg4; reg4=reg3*reg4; T reg6=reg3*reg1; reg1=reg2*reg1;
    T reg7=reg3*reg4; T reg8=reg3*reg5; T reg9=reg2*reg1; reg5=reg2*reg5; T reg10=reg3*reg6;
    reg1=reg3*reg1; reg5=reg5-reg7; reg8=reg7+reg8; reg4=reg2*reg4; reg6=reg2*reg6;
    reg9=reg9-reg10; reg1=reg10+reg1; reg1=reg3*reg1; T reg11=reg3*reg8; T reg12=reg10+reg6;
    reg9=reg2*reg9; reg4=reg7+reg4; reg7=reg2*reg5; reg11=reg7-reg11; reg7=reg3*reg4;
    reg12=reg3*reg12; reg1=reg9-reg1; reg12=reg1-reg12; reg7=reg11-reg7; reg8=reg8/reg7;
    reg5=reg5/reg7; reg12=reg12/reg7; reg1=0.78867513459481286553*elem.pos(0)[0]; reg9=0.5*elem.pos(1)[0]; reg11=0.5*elem.pos(1)[1];
    T reg13=0.5*elem.pos(0)[1]; T reg14=0.5*elem.pos(0)[0]; T reg15=reg5*reg12; T reg16=reg8*reg12; T reg17=0.78867513459481286553*elem.pos(1)[1];
    T reg18=0.78867513459481286553*elem.pos(0)[1]; T reg19=0.21132486540518713447*elem.pos(0)[1]; T reg20=0.21132486540518713447*elem.pos(1)[0]; T reg21=0.21132486540518713447*elem.pos(0)[0]; T reg22=0.78867513459481286553*elem.pos(1)[0];
    T reg23=0.21132486540518713447*elem.pos(1)[1]; T reg24=(*f.m).alpha*reg8; T reg25=0.78867513459481286553*elem.pos(2)[0]; T reg26=(*f.m).alpha*reg5; T reg27=reg9-reg14;
    T reg28=reg8*reg16; T reg29=reg20-reg21; reg20=reg20+reg1; T reg30=0.21132486540518713447*elem.pos(2)[0]; T reg31=0.5*elem.pos(2)[0];
    T reg32=reg23+reg18; T reg33=0.21132486540518713447*elem.pos(2)[1]; T reg34=reg5*reg15; T reg35=reg11+reg13; T reg36=0.5*elem.pos(2)[1];
    T reg37=0.78867513459481286553*elem.pos(2)[1]; reg21=reg22+reg21; reg7=reg4/reg7; reg23=reg23-reg19; reg13=reg11-reg13;
    reg19=reg17+reg19; reg9=reg14+reg9; reg4=0.78867513459481286553*elem.pos(3)[0]; reg32=reg33-reg32; reg11=0.5*elem.pos(3)[1];
    reg35=reg36-reg35; reg23=reg23+reg37; reg21=reg25-reg21; reg14=0.78867513459481286553*elem.pos(3)[1]; reg29=reg25+reg29;
    reg25=0.21132486540518713447*elem.pos(3)[0]; T reg38=0.5*elem.pos(3)[0]; reg7=(*f.m).alpha*reg7; reg24=reg26+reg24; reg1=reg22-reg1;
    reg22=0.21132486540518713447*elem.pos(3)[1]; reg19=reg37-reg19; reg27=reg31+reg27; reg9=reg31-reg9; reg28=reg34-reg28;
    reg18=reg17-reg18; reg20=reg30-reg20; reg13=reg36+reg13; reg17=0.5*vectors[0][indices[0]+1]; reg26=0.5*vectors[0][indices[1]+1];
    reg21=reg21+reg25; reg31=0.5*vectors[0][indices[1]+0]; reg35=reg35+reg11; reg11=reg13-reg11; reg9=reg38+reg9;
    reg23=reg23-reg14; reg38=reg27-reg38; reg30=reg1+reg30; reg1=0.21132486540518713447*PNODE(0).dep[0]; reg13=0.21132486540518713447*PNODE(1).dep[0];
    reg33=reg18+reg33; reg18=0.21132486540518713447*PNODE(1).dep[1]; reg27=0.78867513459481286553*PNODE(0).dep[0]; reg19=reg19+reg22; reg7=reg24+reg7;
    reg16=reg16/reg28; reg24=0.78867513459481286553*PNODE(0).dep[1]; reg29=reg29-reg4; reg34=0.78867513459481286553*PNODE(1).dep[0]; reg36=0.21132486540518713447*PNODE(0).dep[1];
    reg37=0.78867513459481286553*PNODE(1).dep[1]; reg28=reg15/reg28; reg32=reg14+reg32; reg14=0.5*vectors[0][indices[0]+0]; reg20=reg4+reg20;
    reg4=reg31-reg14; reg15=0.5*vectors[0][indices[2]+0]; T reg39=reg36+reg37; reg36=reg18-reg36; T reg40=0.78867513459481286553*PNODE(2).dep[1];
    T reg41=0.5*vectors[0][indices[2]+1]; T reg42=reg21*reg23; reg14=reg31+reg14; reg31=0.78867513459481286553*PNODE(2).dep[0]; T reg43=reg13-reg1;
    reg25=reg30-reg25; reg22=reg33-reg22; reg16=reg16*reg7; reg7=reg28*reg7; reg13=reg13+reg27;
    reg28=0.21132486540518713447*PNODE(2).dep[0]; reg18=reg18+reg24; reg30=0.21132486540518713447*PNODE(2).dep[1]; reg33=reg23*reg20; T reg44=reg29*reg32;
    T reg45=reg17+reg26; T reg46=reg9*reg11; T reg47=reg35*reg38; reg17=reg26-reg17; reg26=reg29*reg19;
    reg1=reg1+reg34; T reg48=reg3*reg0; T reg49=reg2*reg0; reg27=reg34-reg27; reg34=1-(*f.m).resolution;
    T reg50=0.21132486540518713447*PNODE(3).dep[1]; T reg51=reg2*reg49; reg42=reg26-reg42; reg24=reg37-reg24; reg26=reg21*reg22;
    reg37=reg3*reg48; T reg52=reg25*reg19; reg39=reg40-reg39; T reg53=0.78867513459481286553*PNODE(3).dep[0]; reg43=reg43+reg31;
    reg13=reg28-reg13; reg18=reg30-reg18; reg33=reg44-reg33; reg44=0.21132486540518713447*PNODE(3).dep[0]; reg1=reg31-reg1;
    reg17=reg41+reg17; reg31=0.78867513459481286553*PNODE(3).dep[1]; reg46=reg47-reg46; reg45=reg41-reg45; reg41=0.5*vectors[0][indices[3]+1];
    reg36=reg40+reg36; reg14=reg15-reg14; reg15=reg4+reg15; reg16=reg7-reg16; reg4=0.5*vectors[0][indices[3]+0];
    reg7=pow(reg2,2); reg37=reg51-reg37; reg40=pow(reg3,2); reg47=reg22*reg20; reg51=reg25*reg32;
    T reg54=reg29/reg33; T reg55=reg21/reg42; reg28=reg27+reg28; reg18=reg31+reg18; reg27=reg20/reg33;
    T reg56=reg32/reg33; T reg57=reg23/reg33; T reg58=(*f.m).resolution*(*f.m).alpha; reg13=reg53+reg13; reg30=reg24+reg30;
    reg26=reg52-reg26; reg24=reg19/reg42; reg16=reg16*reg34; reg15=reg15-reg4; reg35=reg35/reg46;
    reg31=reg36-reg31; reg11=reg11/reg46; reg45=reg41+reg45; reg9=reg9/reg46; reg41=reg17-reg41;
    reg46=reg38/reg46; reg53=reg43-reg53; reg23=reg23/reg42; reg39=reg50+reg39; reg29=reg29/reg42;
    reg1=reg1+reg44; reg14=reg4+reg14; reg4=reg18*reg57; reg17=reg31*reg56; reg36=reg54*reg13;
    reg38=reg53*reg27; reg43=reg25/reg26; reg21=reg21/reg26; reg57=reg57*reg13; reg52=reg11*reg14;
    T reg59=reg35*reg15; reg50=reg30-reg50; reg56=reg53*reg56; reg27=reg31*reg27; reg19=reg19/reg26;
    reg54=reg54*reg18; reg30=reg46*reg45; T reg60=reg9*reg41; reg44=reg28-reg44; reg28=reg24*reg53;
    reg53=reg55*reg53; T reg61=reg29*reg1; T reg62=reg22/reg26; T reg63=reg39*reg23; reg24=reg31*reg24;
    reg23=reg23*reg1; reg31=reg55*reg31; reg29=reg29*reg39; reg40=reg7-reg40; reg47=reg51-reg47;
    reg49=reg49/reg37; reg48=reg48/reg37; reg16=reg58+reg16; reg7=(*f.m).resolution*reg48; reg51=reg19*reg44;
    reg55=reg62*reg1; reg58=reg21*reg44; reg1=reg43*reg1; reg62=reg39*reg62; reg37=reg40/reg37;
    reg4=reg17-reg4; reg60=reg30-reg60; elem.epsilon[0][1]=reg60; reg17=(*f.m).deltaT*(*f.m).alpha; reg19=reg50*reg19;
    reg16=(*f.m).deltaT*reg16; reg57=reg56-reg57; reg31=reg29-reg31; reg53=reg61-reg53; reg20=reg20/reg47;
    reg52=reg59-reg52; elem.epsilon[0][0]=reg52; reg63=reg24-reg63; reg25=reg25/reg47; reg43=reg39*reg43;
    reg38=reg36-reg38; reg23=reg28-reg23; reg5=reg5*reg34; reg21=reg21*reg50; reg32=reg32/reg47;
    reg27=reg54-reg27; reg22=reg22/reg47; reg8=reg8*reg34; reg24=(*f.m).resolution*reg49; reg27=reg27-reg16;
    reg28=reg52-reg17; reg29=reg60-reg17; reg57=reg57-reg16; reg62=reg19-reg62; reg58=reg1-reg58;
    reg1=reg22*reg13; reg19=reg32*reg44; reg23=reg23-reg16; reg30=reg20*reg50; reg53=reg63+reg53;
    reg36=reg25*reg18; reg4=reg38+reg4; reg38=(*f.m).resolution*reg37; reg21=reg43-reg21; reg34=reg12*reg34;
    reg55=reg51-reg55; reg13=reg25*reg13; reg44=reg20*reg44; reg18=reg22*reg18; reg50=reg32*reg50;
    reg8=reg7+reg8; reg5=reg24+reg5; reg31=reg31-reg16; reg58=reg62+reg58; reg7=reg31*reg8;
    reg55=reg55-reg16; reg12=reg5*reg23; reg18=reg50-reg18; reg20=reg8*reg27; reg44=reg13-reg44;
    reg30=reg36-reg30; reg1=reg19-reg1; reg13=reg5*reg57; reg19=reg8*reg57; reg21=reg21-reg16;
    reg22=reg5*reg27; reg38=reg34+reg38; reg24=reg28*reg49; reg4=0.5*reg4; reg25=reg29*reg48;
    reg28=reg28*reg48; reg29=reg29*reg49; reg53=0.5*reg53; reg32=reg8*reg23; reg34=reg31*reg5;
    reg58=0.5*reg58; reg9=reg15*reg9; reg7=reg12+reg7; reg1=reg1-reg16; reg22=reg19+reg22;
    reg34=reg32+reg34; reg29=reg28+reg29; reg20=reg13+reg20; reg25=reg24+reg25; reg12=reg38*reg53;
    reg13=reg8*reg21; reg44=reg18+reg44; reg15=reg38*reg4; reg18=reg5*reg55; reg30=reg30-reg16;
    reg19=reg5*reg21; reg24=reg8*reg55; reg41=reg35*reg41; reg45=reg11*reg45; reg46=reg14*reg46;
    reg34=reg31*reg34; reg15=2*reg15; reg12=2*reg12; reg7=reg23*reg7; reg22=reg27*reg22;
    reg44=0.5*reg44; reg11=reg25+reg29; reg20=reg57*reg20; reg18=reg13+reg18; reg13=reg30*reg8;
    reg14=reg1*reg5; reg19=reg24+reg19; reg45=reg41-reg45; reg9=reg46-reg9; reg23=reg58*reg38;
    reg24=reg30*reg5; reg27=reg1*reg8; reg55=reg18*reg55; reg20=reg22+reg20; reg23=2*reg23;
    reg18=reg44*reg38; reg13=reg14+reg13; reg12=reg53*reg12; reg7=reg34+reg7; reg24=reg27+reg24;
    reg19=reg21*reg19; reg11=reg11/3; reg15=reg4*reg15; reg9=reg45+reg9; reg19=reg55+reg19;
    reg24=reg30*reg24; reg13=reg1*reg13; reg20=reg15+reg20; reg9=0.5*reg9; elem.epsilon[0][2]=reg9;
    reg23=reg58*reg23; reg18=2*reg18; reg25=reg25-reg11; reg29=reg29-reg11; reg12=reg7+reg12;
    reg1=reg37*reg9; reg13=reg24+reg13; reg25=pow(reg25,2); reg29=pow(reg29,2); reg18=reg44*reg18;
    reg33=reg20*reg33; reg23=reg19+reg23; reg12=reg42*reg12; reg4=2*reg1; reg11=pow(reg11,2);
    reg29=reg25+reg29; reg13=reg18+reg13; reg26=reg23*reg26; reg33=0.25*reg33; reg12=0.25*reg12;
    reg4=reg1*reg4; reg11=reg29+reg11; reg12=reg33+reg12; reg47=reg13*reg47; reg26=0.25*reg26;
    reg12=reg26+reg12; reg4=reg11+reg4; reg60=reg60-reg16; reg52=reg52-reg16; reg47=0.25*reg47;
    reg1=reg5*reg52; reg4=1.5*reg4; reg7=reg8*reg60; reg52=reg8*reg52; reg60=reg5*reg60;
    reg12=reg47+reg12; elem.sigma_von_mises=pow(reg4,0.5); elem.sigma[0][0]=reg1+reg7; elem.ener=reg12/2; elem.sigma[0][1]=reg52+reg60;
    elem.sigma[0][2]=reg38*reg9;
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
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg2*reg4; T reg6=reg3*reg1; reg4=reg3*reg4; reg1=reg2*reg1;
    T reg7=reg3*reg1; T reg8=reg2*reg5; T reg9=reg3*reg4; reg5=reg3*reg5; reg1=reg2*reg1;
    T reg10=reg3*reg6; reg8=reg8-reg9; reg5=reg9+reg5; reg4=reg2*reg4; reg6=reg2*reg6;
    reg1=reg1-reg10; reg7=reg10+reg7; T reg11=reg2*reg8; reg4=reg9+reg4; reg1=reg2*reg1;
    reg9=reg3*reg5; reg7=reg3*reg7; T reg12=reg10+reg6; reg12=reg3*reg12; T reg13=reg3*reg4;
    reg9=reg11-reg9; reg7=reg1-reg7; reg12=reg7-reg12; reg1=1-var_inter[1]; reg7=1-var_inter[0];
    reg13=reg9-reg13; reg8=reg8/reg13; reg12=reg12/reg13; reg5=reg5/reg13; reg9=var_inter[0]*elem.pos(1)[0];
    reg11=var_inter[0]*elem.pos(1)[1]; T reg14=reg7*elem.pos(0)[1]; T reg15=reg7*elem.pos(0)[0]; T reg16=reg1*elem.pos(1)[0]; T reg17=reg1*elem.pos(0)[0];
    T reg18=reg1*elem.pos(1)[1]; T reg19=reg1*elem.pos(0)[1]; reg16=reg16-reg17; T reg20=var_inter[1]*elem.pos(2)[1]; T reg21=var_inter[0]*elem.pos(2)[0];
    T reg22=reg8*reg12; T reg23=var_inter[1]*elem.pos(2)[0]; T reg24=reg5*reg12; T reg25=var_inter[0]*elem.pos(2)[1]; reg18=reg18-reg19;
    T reg26=reg15+reg9; T reg27=reg11+reg14; reg25=reg25-reg27; T reg28=reg7*elem.pos(3)[0]; reg21=reg21-reg26;
    T reg29=var_inter[1]*elem.pos(3)[1]; T reg30=reg7*elem.pos(3)[1]; reg20=reg18+reg20; reg16=reg23+reg16; reg18=var_inter[1]*elem.pos(3)[0];
    reg13=reg4/reg13; reg4=(*f.m).alpha*reg5; reg23=(*f.m).alpha*reg8; T reg31=reg5*reg24; T reg32=reg8*reg22;
    reg20=reg20-reg29; reg13=(*f.m).alpha*reg13; reg4=reg23+reg4; reg28=reg21+reg28; reg16=reg16-reg18;
    reg30=reg25+reg30; reg31=reg32-reg31; reg21=reg28*reg20; reg23=reg30*reg16; reg22=reg22/reg31;
    reg31=reg24/reg31; reg24=reg2*reg0; reg0=reg3*reg0; reg13=reg4+reg13; reg21=reg23-reg21;
    reg4=pow(reg3,2); reg23=pow(reg2,2); reg3=reg3*reg0; reg2=reg2*reg24; reg22=reg22*reg13;
    reg13=reg31*reg13; reg16=reg16/reg21; reg30=reg30/reg21; reg28=reg28/reg21; reg20=reg20/reg21;
    reg4=reg23-reg4; reg3=reg2-reg3; reg2=1-(*f.m).resolution; reg13=reg22-reg13; reg22=var_inter[0]*reg16;
    reg23=var_inter[1]*reg30; reg25=var_inter[0]*reg20; reg31=var_inter[1]*reg28; reg32=reg1*reg28; T reg33=reg7*reg16;
    T reg34=reg7*reg20; reg24=reg24/reg3; reg0=reg0/reg3; reg3=reg4/reg3; reg4=(*f.m).resolution*(*f.m).alpha;
    reg13=reg13*reg2; T reg35=reg1*reg30; T reg36=(*f.m).resolution*reg0; reg8=reg8*reg2; T reg37=(*f.m).resolution*reg24;
    T reg38=reg33+reg31; T reg39=reg34+reg23; reg5=reg5*reg2; reg2=reg12*reg2; reg12=reg32+reg22;
    T reg40=reg35+reg25; T reg41=(*f.m).resolution*reg3; reg13=reg4+reg13; reg4=0.5*reg38; T reg42=reg32-reg33;
    reg5=reg36+reg5; reg8=reg37+reg8; reg41=reg2+reg41; reg2=reg34-reg35; reg36=reg22-reg31;
    reg37=reg23-reg25; T reg43=0.5*reg39; reg13=(*f.m).deltaT*reg13; T reg44=0.5*reg40; T reg45=0.5*reg12;
    T reg46=reg5*reg13; T reg47=reg8*reg13; T reg48=reg41*reg44; T reg49=reg41*reg4; T reg50=reg41*reg43;
    T reg51=0.5*reg42; T reg52=0.5*reg2; T reg53=0.5*reg37; T reg54=reg41*reg45; T reg55=0.5*reg36;
    T reg56=reg8*reg39; reg49=2*reg49; T reg57=reg5*reg40; T reg58=reg8*reg12; T reg59=reg5*reg38;
    T reg60=reg8*reg38; T reg61=reg41*reg53; T reg62=reg47+reg46; T reg63=reg41*reg55; T reg64=reg41*reg51;
    T reg65=reg41*reg52; T reg66=reg5*reg39; T reg67=reg8*reg40; T reg68=2*reg54; T reg69=reg7*var_inter[1];
    reg48=2*reg48; T reg70=reg1*var_inter[0]; T reg71=reg5*reg12; T reg72=2*reg50; T reg73=reg5*reg37;
    T reg74=reg5*reg36; reg63=2*reg63; T reg75=reg12*reg60; T reg76=reg70*elem.f_vol_e[1]; T reg77=reg8*reg42;
    reg61=2*reg61; T reg78=reg44*reg72; T reg79=reg8*reg36; T reg80=reg8*reg37; T reg81=reg44*reg68;
    T reg82=reg5*reg2; T reg83=reg45*reg48; T reg84=reg12*reg57; T reg85=reg40*reg71; T reg86=reg39*reg59;
    T reg87=reg4*reg72; T reg88=reg40*reg56; reg65=2*reg65; T reg89=reg5*reg42; reg64=2*reg64;
    T reg90=reg38*reg58; T reg91=reg48*reg43; T reg92=reg12*reg62; T reg93=reg8*reg2; T reg94=reg69*elem.f_vol_e[0];
    T reg95=reg1*reg7; T reg96=var_inter[0]*var_inter[1]; T reg97=reg39*reg62; T reg98=reg39*reg67; T reg99=reg4*reg68;
    T reg100=reg45*reg49; T reg101=reg38*reg66; T reg102=reg49*reg43; T reg103=reg42*reg60; T reg104=reg45*reg64;
    T reg105=reg2*reg80; T reg106=reg40*reg89; T reg107=reg4*reg63; T reg108=reg45*reg65; T reg109=reg39*reg71;
    T reg110=reg51*reg48; T reg111=reg40*reg67; T reg112=reg45*reg68; T reg113=reg4*reg48; T reg114=reg40*reg93;
    T reg115=reg52*reg72; T reg116=reg36*reg73; T reg117=reg53*reg48; T reg118=reg36*reg58; T reg119=reg53*reg68;
    T reg120=reg36*reg57; T reg121=reg53*reg65; T reg122=reg36*reg77; T reg123=reg53*reg64; T reg124=reg53*reg63;
    T reg125=reg36*reg82; T reg126=reg55*reg72; T reg127=reg37*reg59; T reg128=reg55*reg49; T reg129=reg37*reg56;
    T reg130=reg55*reg61; T reg131=reg39*reg56; T reg132=reg36*reg79; T reg133=reg4*reg49; T reg134=reg40*reg74;
    T reg135=reg45*reg63; T reg136=reg53*reg61; T reg137=reg40*reg80; reg83=reg85+reg83; T reg138=reg51*reg64;
    T reg139=reg2*reg89; T reg140=reg51*reg65; T reg141=reg2*reg67; T reg142=reg51*reg68; T reg143=reg12*reg73;
    T reg144=reg44*reg63; T reg145=reg12*reg58; T reg146=reg44*reg48; reg84=reg81+reg84; T reg147=reg38*reg57;
    T reg148=reg12*reg77; T reg149=reg2*reg71; T reg150=reg2*reg62; T reg151=reg42*reg62; T reg152=reg40*reg62;
    T reg153=reg92-reg76; T reg154=reg37*reg62; T reg155=reg36*reg62; reg100=reg88+reg100; T reg156=reg40*reg59;
    T reg157=reg45*reg72; T reg158=reg44*reg64; T reg159=reg12*reg82; T reg160=reg44*reg65; T reg161=reg51*reg63;
    T reg162=reg45*reg61; T reg163=reg2*reg74; T reg164=reg51*reg61; T reg165=reg2*reg56; T reg166=reg51*reg49;
    reg59=reg2*reg59; T reg167=reg51*reg72; T reg168=reg12*reg79; T reg169=reg44*reg49; T reg170=reg12*reg66;
    reg75=reg78+reg75; T reg171=reg37*reg93; T reg172=reg44*reg61; T reg173=reg55*reg64; T reg174=reg37*reg89;
    T reg175=reg55*reg65; reg67=reg37*reg67; T reg176=reg55*reg68; T reg177=reg37*reg71; T reg178=reg55*reg48;
    T reg179=reg37*reg80; T reg180=reg55*reg63; T reg181=reg37*reg74; T reg182=reg69*elem.f_vol_e[1]; T reg183=reg2*reg93;
    T reg184=reg38*reg60; T reg185=reg72*reg43; T reg186=reg4*reg64; T reg187=reg52*reg65; T reg188=reg42*reg77;
    T reg189=reg4*reg61; reg86=reg86+reg87; T reg190=reg97-reg94; T reg191=reg53*reg72; T reg192=reg43*reg68;
    T reg193=reg52*reg64; reg60=reg36*reg60; reg91=reg90+reg91; T reg194=reg42*reg82; reg82=reg38*reg82;
    reg64=reg43*reg64; T reg195=reg96*elem.f_vol_e[1]; T reg196=reg53*reg49; T reg197=reg96*elem.f_vol_e[0]; reg77=reg38*reg77;
    T reg198=reg36*reg66; T reg199=reg43*reg65; T reg200=reg70*elem.f_vol_e[0]; T reg201=reg38*reg62; reg74=reg39*reg74;
    T reg202=reg95*elem.f_vol_e[1]; T reg203=reg95*elem.f_vol_e[0]; T reg204=reg42*reg79; T reg205=reg52*reg61; reg80=reg39*reg80;
    T reg206=reg52*reg63; reg98=reg99+reg98; T reg207=reg38*reg73; reg73=reg42*reg73; reg63=reg63*reg43;
    reg48=reg52*reg48; reg79=reg38*reg79; reg61=reg61*reg43; reg89=reg39*reg89; reg65=reg4*reg65;
    reg93=reg39*reg93; reg57=reg42*reg57; T reg208=reg52*reg68; reg102=reg101+reg102; reg49=reg52*reg49;
    T reg209=reg42*reg66; T reg210=reg42*reg58; reg59=reg59-reg167; reg48=reg48-reg210; reg161=reg105+reg161;
    reg164=reg163+reg164; reg169=reg169+reg170; reg57=reg57-reg208; reg168=reg172-reg168; reg193=reg194+reg193;
    reg206=reg73+reg206; reg110=reg110-reg149; reg166=reg166-reg165; reg187=reg188+reg187; reg73=reg21*reg75;
    reg105=reg202+reg151; reg163=reg200+reg152; reg172=reg102*reg21; reg153=reg21*reg153; reg188=reg197+reg154;
    reg194=reg195+reg155; reg61=reg79-reg61; reg79=reg21*reg100; reg156=reg156+reg157; reg159=reg158-reg159;
    reg63=reg207-reg63; reg148=reg160-reg148; reg158=reg21*reg84; reg146=reg146+reg145; reg80=reg107-reg80;
    reg143=reg144-reg143; reg173=reg171+reg173; reg107=reg86*reg21; reg175=reg174+reg175; reg199=reg77-reg199;
    reg67=reg67-reg176; reg178=reg178-reg177; reg180=reg179+reg180; reg64=reg82-reg64; reg77=reg91*reg21;
    reg147=reg192+reg147; reg138=reg183+reg138; reg74=reg189-reg74; reg140=reg139+reg140; reg133=reg133+reg131;
    reg141=reg141-reg142; reg184=reg184+reg185; reg82=reg203+reg150; reg111=reg111+reg112; reg121=reg122+reg121;
    reg93=reg186-reg93; reg122=reg21*reg83; reg123=reg125+reg123; reg190=reg21*reg190; reg108=reg106-reg108;
    reg120=reg120-reg119; reg106=reg201+reg182; reg60=reg60-reg191; reg135=reg137-reg135; reg104=reg114-reg104;
    reg49=reg49-reg209; reg117=reg117-reg118; reg127=reg127-reg126; reg103=reg103-reg115; reg136=reg132+reg136;
    reg181=reg130+reg181; reg128=reg128-reg129; reg113=reg113+reg109; reg196=reg196-reg198; reg114=reg21*reg98;
    reg162=reg134-reg162; reg205=reg204+reg205; reg124=reg116+reg124; reg89=reg65-reg89; reg178=reg21*reg178;
    reg60=reg21*reg60; reg181=reg21*reg181; reg65=ponderation*reg122; reg64=reg64*reg21; reg141=reg21*reg141;
    reg184=reg184*reg21; reg133=reg133*reg21; reg140=reg21*reg140; reg180=reg21*reg180; reg74=reg74*reg21;
    reg162=reg21*reg162; reg196=reg21*reg196; reg147=reg147*reg21; reg199=reg199*reg21; reg135=reg21*reg135;
    reg138=reg21*reg138; reg143=reg21*reg143; reg116=ponderation*reg114; reg146=reg21*reg146; reg113=reg21*reg113;
    reg80=reg80*reg21; reg125=ponderation*reg158; reg130=ponderation*reg77; reg49=reg21*reg49; reg148=reg21*reg148;
    reg89=reg21*reg89; reg159=reg21*reg159; reg63=reg63*reg21; reg156=reg21*reg156; reg103=reg21*reg103;
    reg132=ponderation*reg79; reg104=reg21*reg104; reg61=reg61*reg21; reg134=reg21*reg194; reg137=reg21*reg188;
    reg93=reg21*reg93; reg153=ponderation*reg153; reg108=reg21*reg108; reg139=ponderation*reg172; reg144=reg21*reg163;
    reg160=reg21*reg105; reg171=reg21*reg82; reg111=reg21*reg111; reg164=reg21*reg164; reg190=ponderation*reg190;
    reg169=reg21*reg169; reg124=reg21*reg124; reg123=reg21*reg123; reg57=reg21*reg57; reg168=reg21*reg168;
    reg121=reg21*reg121; reg174=ponderation*reg73; reg48=reg21*reg48; reg127=reg21*reg127; reg173=reg21*reg173;
    reg117=reg21*reg117; reg179=ponderation*reg107; reg193=reg21*reg193; reg161=reg21*reg161; reg187=reg21*reg187;
    reg110=reg21*reg110; reg120=reg21*reg120; reg67=reg21*reg67; reg128=reg21*reg128; reg206=reg21*reg206;
    reg136=reg21*reg136; reg205=reg21*reg205; reg59=reg21*reg59; reg175=reg21*reg175; reg166=reg21*reg166;
    reg183=reg21*reg106; T tmp_2_6=-reg132; T tmp_2_0=ponderation*reg104; T tmp_1_1=ponderation*reg187; T tmp_1_7=ponderation*reg103;
    reg103=ponderation*reg134; T vec_5=reg103; T tmp_0_6=ponderation*reg166; T tmp_7_5=ponderation*reg61; reg61=ponderation*reg137;
    T vec_4=reg61; T tmp_0_5=ponderation*reg164; T tmp_6_2=-reg116; T tmp_3_4=ponderation*reg143; T tmp_0_3=ponderation*reg110;
    T tmp_1_4=ponderation*reg206; T tmp_3_3=ponderation*reg146; T tmp_1_5=ponderation*reg205; T tmp_6_4=ponderation*reg80; T tmp_3_2=-reg125;
    T tmp_6_3=ponderation*reg113; T tmp_0_4=ponderation*reg161; T tmp_7_3=-reg130; T tmp_1_3=ponderation*reg48; T tmp_3_1=ponderation*reg148;
    T tmp_1_6=ponderation*reg49; T tmp_6_1=ponderation*reg89; T tmp_5_4=ponderation*reg124; T tmp_3_0=ponderation*reg159; T tmp_1_2=ponderation*reg57;
    T tmp_5_3=ponderation*reg117; T tmp_7_4=ponderation*reg63; T tmp_2_7=ponderation*reg156; T tmp_6_6=ponderation*reg133; T tmp_0_1=ponderation*reg140;
    T tmp_3_7=-reg174; T tmp_6_7=-reg179; T tmp_7_1=ponderation*reg199; T tmp_4_0=ponderation*reg173; T tmp_0_0=ponderation*reg138;
    T tmp_2_4=ponderation*reg135; T tmp_5_6=ponderation*reg196; T tmp_7_2=ponderation*reg147; reg48=ponderation*reg183; T vec_7=reg48;
    T tmp_6_5=ponderation*reg74; T tmp_4_1=ponderation*reg175; T tmp_2_5=ponderation*reg162; T tmp_4_6=ponderation*reg128; T tmp_4_4=ponderation*reg180;
    T tmp_4_2=ponderation*reg67; T tmp_7_0=ponderation*reg64; T tmp_4_3=ponderation*reg178; T tmp_5_5=ponderation*reg136; T tmp_4_5=ponderation*reg181;
    T tmp_6_0=ponderation*reg93; T tmp_5_2=ponderation*reg120; T vec_3=-reg153; T tmp_1_0=ponderation*reg193; T tmp_0_7=ponderation*reg59;
    reg49=ponderation*reg144; T vec_2=reg49; T tmp_2_1=ponderation*reg108; T tmp_7_6=-reg139; reg57=ponderation*reg160;
    T vec_1=reg57; T tmp_5_1=ponderation*reg121; reg59=ponderation*reg171; T vec_0=reg59; T tmp_3_5=ponderation*reg168;
    T tmp_2_2=ponderation*reg111; T tmp_5_7=ponderation*reg60; T vec_6=-reg190; T tmp_7_7=ponderation*reg184; T tmp_0_2=ponderation*reg141;
    T tmp_5_0=ponderation*reg123; T tmp_3_6=ponderation*reg169; T tmp_2_3=-reg65; T tmp_4_7=ponderation*reg127;
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg4*reg1; T reg6=reg3*reg2; reg2=reg4*reg2; reg1=reg3*reg1;
    T reg7=reg3*reg1; T reg8=reg3*reg5; T reg9=reg4*reg2; T reg10=reg3*reg6; reg2=reg3*reg2;
    reg5=reg4*reg5; reg1=reg4*reg1; reg6=reg4*reg6; reg8=reg7+reg8; reg2=reg10+reg2;
    reg5=reg5-reg7; reg9=reg9-reg10; reg5=reg4*reg5; T reg11=reg3*reg2; T reg12=reg7+reg1;
    T reg13=reg4*reg9; reg6=reg10+reg6; reg8=reg3*reg8; reg10=reg3*reg6; reg12=reg3*reg12;
    reg11=reg13-reg11; reg8=reg5-reg8; reg12=reg8-reg12; reg5=1-var_inter[0]; reg10=reg11-reg10;
    reg8=1-var_inter[1]; reg11=reg8*elem.pos(0)[0]; reg13=reg8*elem.pos(1)[0]; T reg14=reg5*elem.pos(0)[1]; T reg15=var_inter[0]*elem.pos(1)[1];
    T reg16=reg5*elem.pos(0)[0]; T reg17=reg8*elem.pos(1)[1]; T reg18=reg8*elem.pos(0)[1]; T reg19=var_inter[0]*elem.pos(1)[0]; reg2=reg2/reg10;
    reg9=reg9/reg10; reg12=reg12/reg10; T reg20=reg9*reg12; T reg21=var_inter[1]*elem.pos(2)[0]; reg13=reg13-reg11;
    reg17=reg17-reg18; T reg22=var_inter[1]*elem.pos(2)[1]; T reg23=var_inter[0]*elem.pos(2)[0]; T reg24=var_inter[0]*elem.pos(2)[1]; T reg25=reg15+reg14;
    T reg26=reg2*reg12; T reg27=reg16+reg19; T reg28=var_inter[1]*elem.pos(3)[1]; reg22=reg17+reg22; reg17=reg9*reg20;
    T reg29=reg2*reg26; T reg30=(*f.m).alpha*reg9; T reg31=var_inter[1]*elem.pos(3)[0]; reg13=reg21+reg13; reg21=(*f.m).alpha*reg2;
    reg10=reg6/reg10; reg23=reg23-reg27; reg6=reg5*elem.pos(3)[0]; T reg32=reg5*elem.pos(3)[1]; reg24=reg24-reg25;
    reg21=reg30+reg21; reg13=reg13-reg31; reg10=(*f.m).alpha*reg10; reg29=reg17-reg29; reg32=reg24+reg32;
    reg22=reg22-reg28; reg6=reg23+reg6; reg17=reg3*reg0; reg26=reg26/reg29; reg29=reg20/reg29;
    reg0=reg4*reg0; reg10=reg21+reg10; reg20=reg32*reg13; reg21=reg6*reg22; reg29=reg29*reg10;
    reg23=pow(reg4,2); reg21=reg20-reg21; reg20=reg3*reg17; reg3=pow(reg3,2); reg4=reg4*reg0;
    reg10=reg26*reg10; reg24=1-(*f.m).resolution; reg20=reg4-reg20; reg3=reg23-reg3; reg22=reg22/reg21;
    reg6=reg6/reg21; reg13=reg13/reg21; reg32=reg32/reg21; reg10=reg29-reg10; reg3=reg3/reg20;
    reg4=var_inter[1]*reg6; reg17=reg17/reg20; reg23=var_inter[1]*reg32; reg20=reg0/reg20; reg0=reg5*reg22;
    reg26=reg5*reg13; reg29=(*f.m).resolution*(*f.m).alpha; reg10=reg10*reg24; reg30=var_inter[0]*reg22; T reg33=reg8*reg32;
    T reg34=reg8*reg6; T reg35=var_inter[0]*reg13; T reg36=(*f.m).resolution*reg3; reg9=reg9*reg24; T reg37=reg33+reg30;
    reg12=reg12*reg24; reg24=reg2*reg24; reg2=(*f.m).resolution*reg20; T reg38=(*f.m).resolution*reg17; reg10=reg29+reg10;
    reg29=reg0+reg23; T reg39=reg26+reg4; reg10=(*f.m).deltaT*reg10; reg36=reg12+reg36; reg12=reg34-reg26;
    reg24=reg38+reg24; reg9=reg2+reg9; reg2=reg0-reg33; reg38=reg34+reg35; T reg40=0.5*reg37;
    T reg41=0.5*reg29; T reg42=reg23-reg30; T reg43=0.5*reg39; T reg44=reg35-reg4; T reg45=0.5*reg42;
    T reg46=reg9*reg10; T reg47=0.5*reg2; T reg48=reg36*reg41; T reg49=0.5*reg12; T reg50=0.5*reg44;
    T reg51=reg24*reg10; T reg52=reg36*reg43; T reg53=0.5*reg38; T reg54=reg36*reg40; T reg55=reg36*reg47;
    T reg56=reg36*reg49; T reg57=reg24*reg38; T reg58=reg46+reg51; T reg59=reg36*reg53; reg54=2*reg54;
    T reg60=reg9*reg39; T reg61=reg36*reg50; T reg62=reg36*reg45; T reg63=reg9*reg29; reg52=2*reg52;
    T reg64=2*reg48; T reg65=reg8*var_inter[0]; T reg66=reg24*reg39; T reg67=reg5*var_inter[1]; T reg68=reg24*reg12;
    reg56=2*reg56; T reg69=reg37*reg63; T reg70=reg53*reg54; T reg71=reg24*reg42; T reg72=reg9*reg2;
    T reg73=reg37*reg57; T reg74=reg67*elem.f_vol_e[0]; T reg75=reg9*reg44; T reg76=reg38*reg60; T reg77=reg8*reg5;
    T reg78=reg40*reg64; T reg79=reg9*reg42; reg61=2*reg61; T reg80=reg24*reg29; T reg81=reg24*reg44;
    T reg82=reg29*reg58; reg62=2*reg62; T reg83=reg65*elem.f_vol_e[1]; T reg84=reg43*reg64; T reg85=reg53*reg52;
    T reg86=reg9*reg12; T reg87=var_inter[0]*var_inter[1]; T reg88=2*reg59; T reg89=reg29*reg66; T reg90=reg24*reg37;
    T reg91=reg38*reg58; T reg92=reg9*reg37; reg55=2*reg55; T reg93=reg9*reg38; T reg94=reg49*reg61;
    T reg95=reg53*reg64; T reg96=reg37*reg66; T reg97=reg53*reg62; T reg98=reg2*reg79; reg85=reg69+reg85;
    T reg99=reg2*reg81; T reg100=reg49*reg54; T reg101=reg44*reg58; T reg102=reg40*reg54; T reg103=reg45*reg62;
    T reg104=reg44*reg75; T reg105=reg50*reg64; T reg106=reg42*reg66; T reg107=reg50*reg52; T reg108=reg42*reg63;
    T reg109=reg50*reg62; T reg110=reg42*reg58; T reg111=reg91-reg83; T reg112=reg29*reg63; T reg113=reg43*reg52;
    T reg114=reg37*reg81; T reg115=reg53*reg61; T reg116=reg49*reg56; T reg117=reg2*reg68; T reg118=reg2*reg72;
    T reg119=reg49*reg55; T reg120=reg67*elem.f_vol_e[1]; T reg121=reg2*reg92; T reg122=reg42*reg81; T reg123=reg50*reg61;
    T reg124=reg42*reg79; reg76=reg78+reg76; T reg125=reg49*reg88; T reg126=reg2*reg57; T reg127=reg38*reg80;
    T reg128=reg40*reg52; T reg129=reg38*reg75; T reg130=reg2*reg58; T reg131=reg49*reg64; T reg132=reg12*reg58;
    T reg133=reg37*reg58; reg66=reg2*reg66; T reg134=reg49*reg52; T reg135=reg40*reg62; T reg136=reg38*reg71;
    T reg137=reg2*reg63; T reg138=reg40*reg61; T reg139=reg49*reg62; T reg140=reg38*reg93; T reg141=reg65*elem.f_vol_e[0];
    T reg142=reg87*elem.f_vol_e[0]; T reg143=reg87*elem.f_vol_e[1]; T reg144=reg12*reg86; T reg145=reg47*reg55; T reg146=reg47*reg64;
    T reg147=reg12*reg60; T reg148=reg12*reg90; T reg149=reg47*reg88; T reg150=reg12*reg93; T reg151=reg47*reg54;
    T reg152=reg12*reg71; T reg153=reg47*reg61; T reg154=reg45*reg64; T reg155=reg44*reg60; T reg156=reg12*reg75;
    T reg157=reg47*reg62; T reg158=reg45*reg52; T reg159=reg44*reg80; T reg160=reg12*reg80; T reg161=reg47*reg52;
    T reg162=reg82-reg74; T reg163=reg39*reg58; reg89=reg89+reg84; reg70=reg73+reg70; T reg164=reg64*reg41;
    reg60=reg39*reg60; T reg165=reg77*elem.f_vol_e[0]; T reg166=reg37*reg79; T reg167=reg37*reg92; T reg168=reg77*elem.f_vol_e[1];
    T reg169=reg53*reg88; reg123=reg124+reg123; reg96=reg96+reg95; reg60=reg60+reg164; reg102=reg102+reg140;
    reg124=reg21*reg76; reg157=reg156+reg157; reg128=reg128+reg127; reg136=reg138-reg136; reg129=reg135-reg129;
    reg161=reg161-reg160; reg162=reg21*reg162; reg66=reg66-reg131; reg135=reg21*reg85; reg153=reg152+reg153;
    reg113=reg113+reg112; reg151=reg151-reg150; reg138=reg143+reg101; reg116=reg118+reg116; reg118=reg142+reg110;
    reg148=reg148-reg149; reg152=reg89*reg21; reg119=reg117+reg119; reg111=reg21*reg111; reg117=reg141+reg133;
    reg145=reg144+reg145; reg144=reg168+reg132; reg121=reg121-reg125; reg156=reg165+reg130; reg107=reg107-reg108;
    reg167=reg167+reg169; reg158=reg158-reg159; reg139=reg99+reg139; reg106=reg106-reg105; reg122=reg109+reg122;
    reg94=reg98+reg94; reg103=reg104+reg103; reg155=reg155-reg154; reg97=reg114-reg97; reg98=reg163+reg120;
    reg147=reg147-reg146; reg99=reg21*reg70; reg134=reg134-reg137; reg100=reg100-reg126; reg115=reg166-reg115;
    reg103=reg21*reg103; reg121=reg21*reg121; reg145=reg21*reg145; reg147=reg21*reg147; reg106=reg21*reg106;
    reg104=reg21*reg156; reg107=reg21*reg107; reg109=reg21*reg144; reg119=reg21*reg119; reg100=reg21*reg100;
    reg114=reg21*reg117; reg166=ponderation*reg152; reg167=reg21*reg167; reg111=ponderation*reg111; reg122=reg21*reg122;
    T reg170=reg21*reg118; T reg171=reg21*reg138; reg113=reg113*reg21; T reg172=ponderation*reg135; reg97=reg21*reg97;
    reg96=reg21*reg96; reg60=reg60*reg21; T reg173=ponderation*reg99; reg102=reg21*reg102; reg115=reg21*reg115;
    reg136=reg21*reg136; reg162=ponderation*reg162; reg134=reg21*reg134; reg66=reg21*reg66; T reg174=reg21*reg98;
    reg129=reg21*reg129; reg128=reg21*reg128; reg157=reg21*reg157; reg139=reg21*reg139; T reg175=ponderation*reg124;
    reg158=reg21*reg158; reg123=reg21*reg123; reg153=reg21*reg153; reg94=reg21*reg94; reg148=reg21*reg148;
    reg151=reg21*reg151; reg155=reg21*reg155; reg116=reg21*reg116; reg161=reg21*reg161; T tmp_0_1=ponderation*reg119;
    T tmp_0_5=ponderation*reg139; T tmp_5_5=ponderation*reg103; T tmp_2_2=ponderation*reg167; T tmp_2_5=ponderation*reg97; reg97=ponderation*reg174;
    T vec_7=reg97; T tmp_2_6=-reg172; T tmp_3_5=ponderation*reg129; T tmp_7_7=ponderation*reg60; T tmp_1_2=ponderation*reg148;
    T tmp_2_7=ponderation*reg96; reg60=ponderation*reg171; T vec_5=reg60; T tmp_0_7=ponderation*reg66; T tmp_0_3=ponderation*reg100;
    T tmp_1_6=ponderation*reg161; T tmp_3_3=ponderation*reg102; T tmp_2_3=-reg173; T tmp_2_4=ponderation*reg115; T tmp_0_6=ponderation*reg134;
    T tmp_3_4=ponderation*reg136; T vec_6=-reg162; T tmp_0_2=ponderation*reg121; T tmp_1_3=ponderation*reg151; T tmp_1_7=ponderation*reg147;
    T tmp_4_7=ponderation*reg106; T tmp_0_4=ponderation*reg94; T tmp_1_4=ponderation*reg153; reg66=ponderation*reg104; T vec_0=reg66;
    T tmp_4_6=ponderation*reg107; T tmp_4_4=ponderation*reg123; reg94=ponderation*reg109; T vec_1=reg94; T tmp_5_6=ponderation*reg158;
    T tmp_5_7=ponderation*reg155; T tmp_3_7=-reg175; reg96=ponderation*reg114; T vec_2=reg96; T tmp_6_7=-reg166;
    T tmp_1_1=ponderation*reg145; T tmp_0_0=ponderation*reg116; T vec_3=-reg111; T tmp_4_5=ponderation*reg122; T tmp_1_5=ponderation*reg157;
    T tmp_3_6=ponderation*reg128; reg100=ponderation*reg170; T vec_4=reg100; T tmp_6_6=ponderation*reg113;
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<false> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1-var_inter[0]; T reg3=1-var_inter[1];
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg6=var_inter[0]*elem.pos(1)[0]; T reg7=reg0*reg1; T reg8=reg3*elem.pos(0)[1];
    T reg9=reg3*elem.pos(1)[1]; T reg10=reg3*elem.pos(0)[0]; T reg11=reg3*elem.pos(1)[0]; T reg12=reg2*elem.pos(0)[1]; T reg13=var_inter[0]*elem.pos(1)[1];
    T reg14=reg2*elem.pos(0)[0]; T reg15=var_inter[1]*elem.pos(2)[1]; reg9=reg9-reg8; T reg16=var_inter[0]*elem.pos(2)[0]; T reg17=reg4*reg7;
    reg11=reg11-reg10; T reg18=var_inter[1]*elem.pos(2)[0]; reg7=reg5*reg7; T reg19=var_inter[0]*elem.pos(2)[1]; T reg20=reg13+reg12;
    T reg21=reg14+reg6; T reg22=reg4*reg1; reg1=reg5*reg1; T reg23=reg5*reg22; reg15=reg9+reg15;
    reg9=var_inter[1]*elem.pos(3)[1]; T reg24=reg4*reg17; reg16=reg16-reg21; T reg25=var_inter[1]*elem.pos(3)[0]; reg11=reg18+reg11;
    reg18=reg5*reg1; reg22=reg4*reg22; T reg26=reg2*elem.pos(3)[1]; T reg27=reg2*elem.pos(3)[0]; reg19=reg19-reg20;
    T reg28=reg5*reg7; reg17=reg5*reg17; reg11=reg11-reg25; reg26=reg19+reg26; reg27=reg16+reg27;
    reg17=reg28+reg17; reg24=reg24-reg28; reg7=reg4*reg7; reg23=reg18+reg23; reg1=reg4*reg1;
    reg22=reg22-reg18; reg15=reg15-reg9; reg16=reg5*reg17; reg19=reg18+reg1; reg22=reg4*reg22;
    T reg29=reg5*reg0; reg7=reg28+reg7; reg28=reg27*reg15; T reg30=reg26*reg11; reg0=reg4*reg0;
    reg23=reg5*reg23; T reg31=reg4*reg24; T reg32=pow(reg5,2); reg16=reg31-reg16; reg31=reg5*reg7;
    reg23=reg22-reg23; reg28=reg30-reg28; reg22=reg4*reg0; reg4=pow(reg4,2); reg19=reg5*reg19;
    reg5=reg5*reg29; reg26=reg26/reg28; reg31=reg16-reg31; reg19=reg23-reg19; reg11=reg11/reg28;
    reg32=reg4-reg32; reg5=reg22-reg5; reg15=reg15/reg28; reg27=reg27/reg28; reg4=var_inter[0]*reg11;
    reg16=1-(*f.m).resolution; reg19=reg19/reg31; reg22=reg2*reg15; reg23=var_inter[1]*reg27; reg30=reg2*reg11;
    T reg33=var_inter[1]*reg26; T reg34=reg3*reg27; reg32=reg32/reg5; T reg35=var_inter[0]*reg15; T reg36=reg3*reg26;
    T reg37=reg34+reg4; reg17=reg17/reg31; T reg38=reg22+reg33; T reg39=reg30+reg23; reg24=reg24/reg31;
    T reg40=(*f.m).resolution*reg32; reg29=reg29/reg5; T reg41=reg36+reg35; reg5=reg0/reg5; reg0=reg19*reg16;
    T reg42=reg34-reg30; T reg43=0.5*reg37; T reg44=0.5*reg41; T reg45=reg17*reg16; T reg46=reg33-reg35;
    reg40=reg0+reg40; reg0=reg4-reg23; T reg47=(*f.m).resolution*reg5; T reg48=reg22-reg36; T reg49=(*f.m).resolution*reg29;
    T reg50=0.5*reg39; T reg51=0.5*reg38; T reg52=reg24*reg16; T reg53=0.5*reg42; T reg54=0.5*reg48;
    T reg55=reg40*reg43; T reg56=0.5*reg0; T reg57=0.5*reg46; T reg58=reg40*reg51; T reg59=reg40*reg50;
    T reg60=reg40*reg44; reg45=reg49+reg45; reg52=reg47+reg52; reg47=reg45*reg38; reg49=reg52*reg37;
    reg60=2*reg60; T reg61=reg40*reg57; T reg62=reg45*reg41; T reg63=reg52*reg39; T reg64=reg52*reg38;
    reg59=2*reg59; T reg65=reg40*reg56; T reg66=reg45*reg39; T reg67=2*reg58; T reg68=reg40*reg53;
    T reg69=reg40*reg54; T reg70=reg52*reg41; T reg71=2*reg55; T reg72=reg45*reg37; T reg73=reg52*reg46;
    reg65=2*reg65; T reg74=reg45*reg46; T reg75=reg45*reg0; reg61=2*reg61; T reg76=reg44*reg71;
    T reg77=reg52*reg42; T reg78=reg52*reg0; T reg79=reg43*reg60; T reg80=reg41*reg72; T reg81=reg52*reg48;
    T reg82=reg50*reg67; reg68=2*reg68; T reg83=reg45*reg42; T reg84=reg38*reg66; T reg85=reg39*reg49;
    reg69=2*reg69; T reg86=reg60*reg51; T reg87=reg41*reg64; T reg88=reg59*reg51; T reg89=reg39*reg47;
    T reg90=reg43*reg59; T reg91=reg45*reg48; T reg92=reg50*reg71; T reg93=reg38*reg70; T reg94=reg44*reg67;
    T reg95=reg37*reg63; T reg96=reg37*reg62; T reg97=reg50*reg65; T reg98=reg42*reg63; T reg99=reg0*reg78;
    T reg100=reg54*reg67; T reg101=reg43*reg71; T reg102=reg53*reg60; T reg103=reg41*reg81; T reg104=reg57*reg61;
    T reg105=reg43*reg68; T reg106=reg41*reg83; T reg107=reg41*reg70; T reg108=reg43*reg69; T reg109=reg57*reg65;
    T reg110=reg0*reg74; T reg111=reg57*reg60; T reg112=reg0*reg49; T reg113=reg57*reg71; T reg114=reg0*reg62;
    T reg115=reg57*reg69; T reg116=reg0*reg77; T reg117=reg57*reg68; T reg118=reg0*reg91; T reg119=reg56*reg67;
    T reg120=reg46*reg66; T reg121=reg56*reg59; T reg122=reg46*reg64; T reg123=reg56*reg61; T reg124=reg38*reg64;
    T reg125=reg50*reg59; T reg126=reg41*reg75; T reg127=reg43*reg65; T reg128=reg41*reg73; reg79=reg80+reg79;
    T reg129=reg53*reg71; T reg130=reg48*reg70; T reg131=reg53*reg69; T reg132=reg48*reg83; T reg133=reg53*reg68;
    T reg134=reg48*reg81; T reg135=reg46*reg75; T reg136=reg56*reg65; T reg137=reg46*reg73; T reg138=reg56*reg60;
    T reg139=reg46*reg72; T reg140=reg56*reg71; reg70=reg46*reg70; T reg141=reg56*reg69; T reg142=reg46*reg83;
    T reg143=reg56*reg68; T reg144=reg46*reg81; reg95=reg94+reg95; T reg145=reg37*reg47; T reg146=reg44*reg59;
    T reg147=reg37*reg78; T reg148=reg53*reg67; T reg149=reg48*reg73; T reg150=reg53*reg65; T reg151=reg43*reg61;
    T reg152=reg48*reg75; T reg153=reg53*reg61; T reg154=reg48*reg64; T reg155=reg53*reg59; T reg156=reg48*reg66;
    T reg157=reg44*reg61; T reg158=reg37*reg74; T reg159=reg44*reg65; T reg160=reg37*reg49; T reg161=reg44*reg60;
    reg96=reg76+reg96; T reg162=reg39*reg62; T reg163=reg37*reg77; T reg164=reg44*reg69; T reg165=reg37*reg91;
    T reg166=reg44*reg68; T reg167=reg43*reg67; reg66=reg41*reg66; reg90=reg87+reg90; T reg168=reg48*reg72;
    T reg169=reg50*reg69; T reg170=reg67*reg51; T reg171=reg42*reg47; reg62=reg42*reg62; T reg172=reg50*reg61;
    reg84=reg84+reg82; T reg173=reg54*reg59; T reg174=reg54*reg69; reg81=reg38*reg81; T reg175=reg50*reg68;
    T reg176=reg51*reg71; T reg177=reg42*reg77; reg86=reg85+reg86; T reg178=reg54*reg68; T reg179=reg0*reg47;
    T reg180=reg39*reg91; reg91=reg42*reg91; T reg181=reg57*reg67; T reg182=reg0*reg63; reg68=reg51*reg68;
    reg77=reg39*reg77; reg59=reg57*reg59; reg69=reg51*reg69; reg75=reg38*reg75; T reg183=reg38*reg72;
    reg73=reg38*reg73; T reg184=reg50*reg60; T reg185=reg42*reg78; T reg186=reg54*reg61; T reg187=reg54*reg65;
    T reg188=reg39*reg74; reg74=reg42*reg74; reg65=reg65*reg51; reg93=reg92+reg93; reg60=reg54*reg60;
    reg83=reg38*reg83; reg63=reg39*reg63; T reg189=reg54*reg71; reg88=reg89+reg88; T reg190=reg42*reg49;
    reg61=reg61*reg51; reg78=reg39*reg78; reg102=reg102-reg168; reg178=reg91+reg178; reg174=reg177+reg174;
    reg153=reg152+reg153; reg150=reg149+reg150; reg155=reg155-reg154; reg62=reg62-reg189; reg60=reg60-reg190;
    reg156=reg156-reg148; reg91=reg84*reg28; reg187=reg74+reg187; reg104=reg99+reg104; reg158=reg159-reg158;
    reg73=reg97-reg73; reg161=reg161+reg160; reg74=reg28*reg96; reg163=reg164-reg163; reg65=reg188-reg65;
    reg165=reg166-reg165; reg66=reg66+reg167; reg97=reg28*reg90; reg61=reg78-reg61; reg130=reg130-reg129;
    reg78=reg88*reg28; reg131=reg132+reg131; reg63=reg63+reg170; reg133=reg134+reg133; reg125=reg125+reg124;
    reg75=reg172-reg75; reg136=reg137+reg136; reg162=reg176+reg162; reg138=reg138-reg139; reg70=reg70-reg140;
    reg99=reg86*reg28; reg141=reg142+reg141; reg143=reg144+reg143; reg68=reg180-reg68; reg132=reg28*reg95;
    reg146=reg146+reg145; reg69=reg77-reg69; reg147=reg157-reg147; reg105=reg103-reg105; reg173=reg173-reg171;
    reg115=reg116+reg115; reg182=reg182-reg181; reg117=reg118+reg117; reg127=reg128-reg127; reg151=reg126-reg151;
    reg114=reg114-reg113; reg135=reg123+reg135; reg83=reg169-reg83; reg120=reg120-reg119; reg77=reg28*reg79;
    reg81=reg175-reg81; reg184=reg184+reg183; reg107=reg107+reg101; reg121=reg121-reg122; reg109=reg110+reg109;
    reg98=reg98-reg100; reg111=reg111-reg112; reg108=reg106-reg108; reg103=reg28*reg93; reg59=reg59-reg179;
    reg186=reg185+reg186; reg135=reg28*reg135; reg143=reg28*reg143; reg151=reg28*reg151; reg66=reg28*reg66;
    reg106=ponderation*reg132; reg65=reg65*reg28; reg141=reg28*reg141; reg68=reg68*reg28; reg165=reg28*reg165;
    reg182=reg28*reg182; reg110=ponderation*reg78; reg131=reg28*reg131; reg107=reg28*reg107; reg108=reg28*reg108;
    reg63=reg63*reg28; reg130=reg28*reg130; reg81=reg28*reg81; reg133=reg28*reg133; reg125=reg125*reg28;
    reg116=ponderation*reg77; reg75=reg75*reg28; reg69=reg69*reg28; reg61=reg61*reg28; reg136=reg28*reg136;
    reg83=reg28*reg83; reg127=reg28*reg127; reg138=reg28*reg138; reg118=ponderation*reg97; reg162=reg162*reg28;
    reg105=reg28*reg105; reg70=reg28*reg70; reg123=ponderation*reg99; reg155=reg28*reg155; reg117=reg28*reg117;
    reg102=reg28*reg102; reg184=reg28*reg184; reg126=ponderation*reg74; reg174=reg28*reg174; reg153=reg28*reg153;
    reg115=reg28*reg115; reg60=reg28*reg60; reg73=reg73*reg28; reg186=reg28*reg186; reg173=reg28*reg173;
    reg62=reg28*reg62; reg161=reg28*reg161; reg114=reg28*reg114; reg150=reg28*reg150; reg98=reg28*reg98;
    reg146=reg28*reg146; reg59=reg28*reg59; reg121=reg28*reg121; reg147=reg28*reg147; reg128=ponderation*reg91;
    reg163=reg28*reg163; reg187=reg28*reg187; reg156=reg28*reg156; reg109=reg28*reg109; reg120=reg28*reg120;
    reg104=reg28*reg104; reg134=ponderation*reg103; reg158=reg28*reg158; reg111=reg28*reg111; reg178=reg28*reg178;
    T tmp_6_3=ponderation*reg184; T tmp_0_2=ponderation*reg130; T tmp_3_4=ponderation*reg158; T tmp_7_6=-reg110; T tmp_3_3=ponderation*reg161;
    T tmp_7_5=ponderation*reg61; T tmp_1_6=ponderation*reg173; T tmp_6_4=ponderation*reg73; T tmp_2_0=ponderation*reg105; T tmp_3_2=-reg126;
    T tmp_2_6=-reg118; T tmp_6_2=-reg134; T tmp_6_1=ponderation*reg83; T tmp_7_3=-reg123; T tmp_3_1=ponderation*reg163;
    T tmp_2_7=ponderation*reg66; T tmp_7_4=ponderation*reg65; T tmp_1_7=ponderation*reg98; T tmp_3_0=ponderation*reg165; T tmp_3_6=ponderation*reg146;
    T tmp_3_5=ponderation*reg147; T tmp_4_6=ponderation*reg121; T tmp_0_7=ponderation*reg156; T tmp_6_7=-reg128; T tmp_4_7=ponderation*reg120;
    T tmp_5_5=ponderation*reg104; T tmp_0_6=ponderation*reg155; T tmp_1_0=ponderation*reg178; T tmp_5_0=ponderation*reg117; T tmp_0_5=ponderation*reg153;
    T tmp_1_1=ponderation*reg174; T tmp_5_1=ponderation*reg115; T tmp_1_5=ponderation*reg186; T tmp_1_2=ponderation*reg62; T tmp_0_4=ponderation*reg150;
    T tmp_5_2=ponderation*reg114; T tmp_1_3=ponderation*reg60; T tmp_0_3=ponderation*reg102; T tmp_5_3=ponderation*reg111; T tmp_1_4=ponderation*reg187;
    T tmp_5_4=ponderation*reg109; T tmp_2_1=ponderation*reg108; T tmp_0_1=ponderation*reg131; T tmp_2_2=ponderation*reg107; T tmp_7_7=ponderation*reg63;
    T tmp_6_0=ponderation*reg81; T tmp_0_0=ponderation*reg133; T tmp_6_6=ponderation*reg125; T tmp_7_1=ponderation*reg69; T tmp_4_4=ponderation*reg136;
    T tmp_2_3=-reg116; T tmp_4_3=ponderation*reg138; T tmp_2_4=ponderation*reg127; T tmp_7_2=ponderation*reg162; T tmp_4_2=ponderation*reg70;
    T tmp_5_7=ponderation*reg182; T tmp_4_1=ponderation*reg141; T tmp_6_5=ponderation*reg75; T tmp_4_0=ponderation*reg143; T tmp_2_5=ponderation*reg151;
    T tmp_3_7=-reg106; T tmp_4_5=ponderation*reg135; T tmp_7_0=ponderation*reg68; T tmp_5_6=ponderation*reg59;
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<true> &matrix_is_sym,
      const Number<true> &assemble_mat,
      const Number<false> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=1-var_inter[1]; T reg2=pow(reg0,2); T reg3=1-var_inter[0];
    T reg4=reg1*elem.pos(0)[1]; T reg5=reg1*elem.pos(1)[1]; T reg6=reg1*elem.pos(0)[0]; T reg7=reg1*elem.pos(1)[0]; T reg8=var_inter[0]*elem.pos(1)[0];
    T reg9=reg3*elem.pos(0)[1]; T reg10=var_inter[0]*elem.pos(1)[1]; T reg11=reg3*elem.pos(0)[0]; T reg12=reg0*reg2; T reg13=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg14=1.0/(*f.m).elastic_modulus; T reg15=var_inter[0]*elem.pos(2)[0]; T reg16=var_inter[0]*elem.pos(2)[1]; T reg17=reg10+reg9; T reg18=var_inter[1]*elem.pos(2)[1];
    T reg19=reg11+reg8; reg5=reg5-reg4; T reg20=reg13*reg12; T reg21=var_inter[1]*elem.pos(2)[0]; reg7=reg7-reg6;
    T reg22=reg13*reg2; reg12=reg14*reg12; reg2=reg14*reg2; T reg23=reg13*reg2; T reg24=reg3*elem.pos(3)[1];
    reg16=reg16-reg17; T reg25=var_inter[1]*elem.pos(3)[1]; reg18=reg5+reg18; reg5=reg14*reg12; T reg26=reg3*elem.pos(3)[0];
    reg15=reg15-reg19; T reg27=var_inter[1]*elem.pos(3)[0]; T reg28=reg13*reg22; reg2=reg14*reg2; T reg29=reg13*reg20;
    reg12=reg13*reg12; reg7=reg21+reg7; reg7=reg7-reg27; reg18=reg18-reg25; reg24=reg16+reg24;
    reg26=reg15+reg26; reg23=reg28+reg23; reg5=reg5-reg29; reg12=reg29+reg12; reg20=reg14*reg20;
    reg2=reg2-reg28; reg22=reg14*reg22; reg15=reg14*reg0; reg0=reg13*reg0; reg16=reg28+reg22;
    reg21=reg14*reg5; reg20=reg29+reg20; reg2=reg14*reg2; reg23=reg13*reg23; reg29=reg13*reg12;
    T reg30=reg26*reg18; T reg31=reg24*reg7; T reg32=pow(reg13,2); reg30=reg31-reg30; reg31=pow(reg14,2);
    T reg33=reg13*reg0; reg23=reg2-reg23; reg29=reg21-reg29; reg16=reg13*reg16; reg14=reg14*reg15;
    reg13=reg13*reg20; reg13=reg29-reg13; reg16=reg23-reg16; reg33=reg14-reg33; reg32=reg31-reg32;
    reg24=reg24/reg30; reg7=reg7/reg30; reg26=reg26/reg30; reg18=reg18/reg30; reg2=reg1*reg24;
    reg32=reg32/reg33; reg16=reg16/reg13; reg14=var_inter[1]*reg24; reg21=1-(*f.m).resolution; reg23=var_inter[1]*reg26;
    reg29=reg3*reg18; reg31=reg3*reg7; T reg34=var_inter[0]*reg18; T reg35=reg16*reg21; T reg36=(*f.m).resolution*reg32;
    T reg37=reg1*reg26; reg0=reg0/reg33; T reg38=var_inter[0]*reg7; T reg39=reg2+reg34; reg33=reg15/reg33;
    reg5=reg5/reg13; reg15=reg31+reg23; T reg40=reg29+reg14; reg12=reg12/reg13; T reg41=0.5*reg15;
    T reg42=reg37-reg31; T reg43=reg29-reg2; T reg44=reg38-reg23; T reg45=reg5*reg21; T reg46=0.5*reg40;
    T reg47=reg14-reg34; T reg48=reg12*reg21; T reg49=(*f.m).resolution*reg33; T reg50=(*f.m).resolution*reg0; T reg51=0.5*reg39;
    T reg52=reg37+reg38; reg36=reg35+reg36; reg35=0.5*reg44; T reg53=0.5*reg52; T reg54=0.5*reg47;
    T reg55=reg36*reg46; T reg56=reg36*reg41; T reg57=reg36*reg51; T reg58=0.5*reg43; T reg59=0.5*reg42;
    reg45=reg49+reg45; reg48=reg50+reg48; reg49=reg45*reg40; reg50=reg36*reg54; reg56=2*reg56;
    T reg60=reg48*reg15; T reg61=2*reg55; T reg62=reg45*reg15; T reg63=reg36*reg35; reg57=2*reg57;
    T reg64=reg36*reg59; T reg65=reg36*reg58; T reg66=reg36*reg53; T reg67=reg48*reg52; reg50=2*reg50;
    T reg68=reg39*reg49; T reg69=2*reg66; T reg70=reg40*reg60; T reg71=reg48*reg40; T reg72=reg48*reg44;
    reg63=2*reg63; T reg73=reg45*reg44; T reg74=reg51*reg61; T reg75=reg52*reg62; T reg76=reg45*reg47;
    T reg77=reg48*reg47; T reg78=reg53*reg56; T reg79=reg45*reg42; T reg80=reg53*reg57; T reg81=reg45*reg39;
    T reg82=reg39*reg67; reg65=2*reg65; T reg83=reg48*reg39; T reg84=reg48*reg42; T reg85=reg45*reg43;
    T reg86=reg41*reg61; T reg87=reg45*reg52; reg64=2*reg64; T reg88=reg35*reg56; T reg89=reg35*reg61;
    T reg90=reg47*reg49; T reg91=reg59*reg63; T reg92=reg43*reg76; T reg93=reg59*reg57; T reg94=reg44*reg73;
    T reg95=reg54*reg50; T reg96=reg53*reg50; T reg97=reg35*reg50; T reg98=reg43*reg72; T reg99=reg47*reg60;
    T reg100=reg40*reg49; T reg101=reg41*reg56; T reg102=reg39*reg72; T reg103=reg52*reg77; T reg104=reg51*reg63;
    T reg105=reg52*reg87; T reg106=reg51*reg57; T reg107=reg53*reg61; T reg108=reg39*reg60; reg78=reg68+reg78;
    T reg109=reg43*reg67; T reg110=reg43*reg81; T reg111=reg59*reg65; T reg112=reg43*reg84; T reg113=reg59*reg64;
    T reg114=reg43*reg85; T reg115=reg59*reg69; T reg116=reg47*reg72; T reg117=reg35*reg63; T reg118=reg47*reg76;
    reg75=reg74+reg75; T reg119=reg52*reg71; T reg120=reg51*reg56; T reg121=reg52*reg73; T reg122=reg51*reg50;
    T reg123=reg59*reg61; reg60=reg43*reg60; T reg124=reg59*reg56; T reg125=reg43*reg49; T reg126=reg59*reg50;
    T reg127=reg53*reg69; T reg128=reg42*reg73; T reg129=reg58*reg57; T reg130=reg42*reg87; T reg131=reg58*reg50;
    T reg132=reg39*reg81; T reg133=reg58*reg69; T reg134=reg42*reg83; T reg135=reg42*reg62; reg80=reg82+reg80;
    T reg136=reg42*reg77; T reg137=reg58*reg65; T reg138=reg42*reg79; T reg139=reg42*reg71; T reg140=reg58*reg56;
    T reg141=reg58*reg63; T reg142=reg54*reg61; reg70=reg70+reg86; T reg143=reg44*reg62; T reg144=reg39*reg76;
    T reg145=reg58*reg61; T reg146=reg44*reg71; T reg147=reg61*reg46; T reg148=reg53*reg63; reg62=reg15*reg62;
    T reg149=reg54*reg56; reg120=reg120+reg119; reg121=reg122-reg121; reg131=reg128+reg131; reg60=reg60-reg123;
    reg140=reg140-reg139; reg124=reg124-reg125; reg149=reg149-reg146; reg103=reg104-reg103; reg106=reg106+reg105;
    reg62=reg62+reg147; reg108=reg108+reg107; reg104=reg30*reg78; reg101=reg101+reg100; reg122=reg70*reg30;
    reg110=reg110-reg115; reg137=reg138+reg137; reg111=reg112+reg111; reg134=reg134-reg133; reg113=reg114+reg113;
    reg129=reg129-reg130; reg117=reg118+reg117; reg112=reg30*reg75; reg141=reg136+reg141; reg91=reg92+reg91;
    reg88=reg88-reg90; reg135=reg135-reg145; reg92=reg30*reg80; reg99=reg99-reg89; reg116=reg97+reg116;
    reg143=reg143-reg142; reg132=reg132+reg127; reg95=reg94+reg95; reg126=reg98+reg126; reg96=reg102-reg96;
    reg148=reg144-reg148; reg93=reg93-reg109; reg134=reg30*reg134; reg95=reg30*reg95; reg132=reg30*reg132;
    reg113=reg30*reg113; reg99=reg30*reg99; reg137=reg30*reg137; reg111=reg30*reg111; reg129=reg30*reg129;
    reg88=reg30*reg88; reg94=ponderation*reg122; reg110=reg30*reg110; reg97=ponderation*reg92; reg101=reg101*reg30;
    reg116=reg30*reg116; reg98=ponderation*reg104; reg108=reg30*reg108; reg62=reg62*reg30; reg96=reg30*reg96;
    reg106=reg30*reg106; reg148=reg30*reg148; reg103=reg30*reg103; reg60=reg30*reg60; reg140=reg30*reg140;
    reg131=reg30*reg131; reg91=reg30*reg91; reg121=reg30*reg121; reg143=reg30*reg143; reg120=reg30*reg120;
    reg124=reg30*reg124; reg141=reg30*reg141; reg117=reg30*reg117; reg149=reg30*reg149; reg93=reg30*reg93;
    reg102=ponderation*reg112; reg135=reg30*reg135; reg126=reg30*reg126; T tmp_4_5=ponderation*reg116; T tmp_2_4=ponderation*reg148;
    T tmp_3_4=ponderation*reg103; T tmp_1_5=ponderation*reg131; T tmp_2_6=-reg98; T tmp_3_3=ponderation*reg106; T tmp_5_7=ponderation*reg143;
    T tmp_7_7=ponderation*reg62; T tmp_2_3=-reg97; T tmp_5_6=ponderation*reg149; T tmp_2_7=ponderation*reg108; T tmp_0_5=ponderation*reg126;
    T tmp_5_5=ponderation*reg95; T tmp_2_5=ponderation*reg96; T tmp_0_6=ponderation*reg124; T tmp_4_4=ponderation*reg117; T tmp_1_3=ponderation*reg129;
    T tmp_1_7=ponderation*reg135; T tmp_3_7=-reg102; T tmp_1_2=ponderation*reg134; T tmp_0_3=ponderation*reg93; T tmp_0_0=ponderation*reg113;
    T tmp_4_7=ponderation*reg99; T tmp_1_4=ponderation*reg141; T tmp_3_6=ponderation*reg120; T tmp_1_1=ponderation*reg137; T tmp_2_2=ponderation*reg132;
    T tmp_1_6=ponderation*reg140; T tmp_0_1=ponderation*reg111; T tmp_4_6=ponderation*reg88; T tmp_3_5=ponderation*reg121; T tmp_6_7=-reg94;
    T tmp_0_4=ponderation*reg91; T tmp_0_2=ponderation*reg110; T tmp_0_7=ponderation*reg60; T tmp_6_6=ponderation*reg101;
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TMA &matrix,
      TVE &sollicitation,
      TVEVE &vectors,
      const Number<symmetric_version> &matrix_is_sym,
      const Number<false> &assemble_mat,
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg4*reg1; reg1=reg3*reg1; T reg6=reg3*reg2; reg2=reg4*reg2;
    T reg7=reg3*reg6; T reg8=reg3*reg2; T reg9=reg4*reg5; T reg10=reg3*reg1; reg2=reg4*reg2;
    reg5=reg3*reg5; reg2=reg2-reg7; reg8=reg7+reg8; reg6=reg4*reg6; reg5=reg10+reg5;
    reg1=reg4*reg1; reg9=reg9-reg10; T reg11=reg4*reg2; reg6=reg7+reg6; reg9=reg4*reg9;
    reg5=reg3*reg5; reg7=reg3*reg8; T reg12=reg10+reg1; T reg13=reg3*reg6; reg7=reg11-reg7;
    reg5=reg9-reg5; reg12=reg3*reg12; reg12=reg5-reg12; reg13=reg7-reg13; reg12=reg12/reg13;
    reg8=reg8/reg13; reg2=reg2/reg13; reg5=reg8*reg12; reg7=reg2*reg12; reg9=(*f.m).alpha*reg8;
    reg11=1-var_inter[1]; T reg14=(*f.m).alpha*reg2; T reg15=reg8*reg5; T reg16=reg2*reg7; reg13=reg6/reg13;
    reg6=1-var_inter[0]; T reg17=reg6*elem.pos(0)[0]; reg13=(*f.m).alpha*reg13; reg9=reg14+reg9; reg14=var_inter[0]*elem.pos(1)[1];
    T reg18=reg6*elem.pos(0)[1]; T reg19=var_inter[0]*elem.pos(1)[0]; reg15=reg16-reg15; reg16=reg11*elem.pos(0)[1]; T reg20=reg11*elem.pos(1)[1];
    T reg21=reg11*elem.pos(0)[0]; T reg22=reg11*elem.pos(1)[0]; T reg23=reg4*reg0; T reg24=reg17+reg19; reg0=reg3*reg0;
    T reg25=reg14+reg18; T reg26=var_inter[0]*elem.pos(2)[1]; reg13=reg9+reg13; reg9=var_inter[0]*elem.pos(2)[0]; reg7=reg7/reg15;
    reg15=reg5/reg15; reg5=var_inter[1]*elem.pos(2)[0]; T reg27=var_inter[1]*elem.pos(2)[1]; reg20=reg20-reg16; reg22=reg22-reg21;
    T reg28=reg4*reg23; T reg29=reg3*reg0; reg26=reg26-reg25; T reg30=reg6*elem.pos(3)[1]; reg15=reg15*reg13;
    reg13=reg7*reg13; reg22=reg5+reg22; reg5=var_inter[1]*elem.pos(3)[0]; reg27=reg20+reg27; reg7=var_inter[1]*elem.pos(3)[1];
    reg9=reg9-reg24; reg20=reg6*elem.pos(3)[0]; reg15=reg13-reg15; reg27=reg27-reg7; reg22=reg22-reg5;
    reg29=reg28-reg29; reg13=1-(*f.m).resolution; reg30=reg26+reg30; reg20=reg9+reg20; reg0=reg0/reg29;
    reg23=reg23/reg29; reg9=reg30*reg22; reg26=reg20*reg27; reg15=reg15*reg13; reg28=(*f.m).resolution*(*f.m).alpha;
    reg2=reg2*reg13; reg8=reg8*reg13; T reg31=(*f.m).resolution*reg0; T reg32=(*f.m).resolution*reg23; reg26=reg9-reg26;
    reg15=reg28+reg15; reg30=reg30/reg26; reg22=reg22/reg26; reg2=reg32+reg2; reg8=reg31+reg8;
    reg15=(*f.m).deltaT*reg15; reg20=reg20/reg26; reg27=reg27/reg26; reg9=var_inter[1]*reg30; reg28=var_inter[0]*reg22;
    reg31=reg11*reg20; reg32=reg6*reg27; T reg33=reg8*reg15; T reg34=reg2*reg15; T reg35=reg11*reg30;
    T reg36=reg11*var_inter[0]; T reg37=reg34+reg33; T reg38=reg6*var_inter[1]; T reg39=reg6*reg22; T reg40=var_inter[1]*reg20;
    T reg41=reg31+reg28; T reg42=reg32+reg9; T reg43=var_inter[0]*reg27; T reg44=reg9-reg43; T reg45=reg28-reg40;
    T reg46=reg39+reg40; T reg47=reg42*reg37; T reg48=reg41*reg37; T reg49=reg38*elem.f_vol_e[0]; T reg50=reg36*elem.f_vol_e[1];
    T reg51=reg35+reg43; T reg52=reg11*reg6; T reg53=var_inter[0]*var_inter[1]; T reg54=reg31-reg39; T reg55=reg32-reg35;
    T reg56=reg54*reg37; T reg57=reg51*reg37; T reg58=reg55*reg37; T reg59=reg38*elem.f_vol_e[1]; T reg60=reg48-reg50;
    T reg61=reg46*reg37; T reg62=reg44*reg37; T reg63=reg45*reg37; T reg64=reg47-reg49; T reg65=reg52*elem.f_vol_e[0];
    T reg66=reg52*elem.f_vol_e[1]; T reg67=reg53*elem.f_vol_e[1]; T reg68=reg53*elem.f_vol_e[0]; T reg69=reg36*elem.f_vol_e[0]; T reg70=reg68+reg62;
    reg60=reg26*reg60; T reg71=reg67+reg63; T reg72=reg69+reg57; T reg73=reg66+reg56; reg64=reg26*reg64;
    T reg74=reg65+reg58; T reg75=reg61+reg59; T reg76=reg26*reg71; T reg77=reg26*reg70; reg64=ponderation*reg64;
    reg60=ponderation*reg60; T reg78=reg26*reg75; T reg79=reg26*reg72; T reg80=reg26*reg74; T reg81=reg26*reg73;
    T reg82=ponderation*reg78; T vec_7=reg82; T reg83=ponderation*reg79; T vec_2=reg83; T vec_3=-reg60;
    T vec_6=-reg64; reg60=ponderation*reg81; T vec_1=reg60; reg64=ponderation*reg77; T vec_4=reg64;
    T reg84=ponderation*reg76; T vec_5=reg84; T reg85=ponderation*reg80; T vec_0=reg85;
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
      TVE &sollicitation,
      TVEVE &vectors,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices ) { 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg3*reg4; T reg6=reg3*reg1; reg1=reg2*reg1; reg4=reg2*reg4;
    T reg7=reg2*reg4; T reg8=reg3*reg1; reg4=reg3*reg4; T reg9=reg3*reg6; reg1=reg2*reg1;
    T reg10=reg3*reg5; reg1=reg1-reg9; reg8=reg9+reg8; reg6=reg2*reg6; reg5=reg2*reg5;
    reg4=reg10+reg4; reg7=reg7-reg10; reg1=reg2*reg1; T reg11=reg3*reg4; T reg12=reg2*reg7;
    reg5=reg10+reg5; reg10=reg9+reg6; reg8=reg3*reg8; reg11=reg12-reg11; reg12=reg3*reg5;
    reg10=reg3*reg10; reg8=reg1-reg8; reg10=reg8-reg10; reg1=1-var_inter[0]; reg12=reg11-reg12;
    reg8=1-var_inter[1]; reg4=reg4/reg12; reg11=reg1*elem.pos(0)[1]; T reg13=var_inter[0]*elem.pos(1)[1]; T reg14=var_inter[0]*elem.pos(1)[0];
    T reg15=reg1*elem.pos(0)[0]; T reg16=reg8*elem.pos(0)[1]; T reg17=reg8*elem.pos(1)[1]; reg7=reg7/reg12; reg10=reg10/reg12;
    T reg18=reg8*elem.pos(0)[0]; T reg19=reg8*elem.pos(1)[0]; reg17=reg17-reg16; T reg20=var_inter[1]*elem.pos(2)[1]; T reg21=var_inter[0]*elem.pos(2)[0];
    reg19=reg19-reg18; T reg22=var_inter[1]*elem.pos(2)[0]; T reg23=var_inter[0]*elem.pos(2)[1]; T reg24=reg13+reg11; T reg25=reg15+reg14;
    T reg26=reg7*reg10; T reg27=reg4*reg10; T reg28=var_inter[1]*elem.pos(3)[1]; reg20=reg17+reg20; reg17=(*f.m).alpha*reg4;
    T reg29=(*f.m).alpha*reg7; T reg30=var_inter[1]*elem.pos(3)[0]; reg23=reg23-reg24; T reg31=reg1*elem.pos(3)[1]; T reg32=reg1*elem.pos(3)[0];
    reg21=reg21-reg25; reg12=reg5/reg12; reg5=reg7*reg26; T reg33=reg4*reg27; reg19=reg22+reg19;
    reg12=(*f.m).alpha*reg12; reg22=reg1*vectors[0][indices[0]+1]; reg17=reg29+reg17; reg29=var_inter[0]*vectors[0][indices[1]+1]; T reg34=reg8*vectors[0][indices[1]+1];
    T reg35=reg8*vectors[0][indices[0]+1]; reg33=reg5-reg33; reg5=reg8*vectors[0][indices[0]+0]; T reg36=reg1*vectors[0][indices[0]+0]; T reg37=var_inter[0]*vectors[0][indices[1]+0];
    T reg38=reg8*vectors[0][indices[1]+0]; reg32=reg21+reg32; reg31=reg23+reg31; reg19=reg19-reg30; reg20=reg20-reg28;
    reg26=reg26/reg33; reg35=reg34-reg35; reg21=var_inter[1]*vectors[0][indices[2]+1]; reg33=reg27/reg33; reg29=reg22+reg29;
    reg22=var_inter[1]*vectors[0][indices[2]+0]; reg23=var_inter[0]*vectors[0][indices[2]+0]; reg12=reg17+reg12; reg37=reg36+reg37; reg17=var_inter[0]*vectors[0][indices[2]+1];
    reg5=reg38-reg5; reg27=reg31*reg19; reg34=reg32*reg20; reg5=reg22+reg5; reg22=var_inter[1]*vectors[0][indices[3]+0];
    reg29=reg17-reg29; reg26=reg26*reg12; reg17=reg1*vectors[0][indices[3]+1]; reg35=reg21+reg35; reg12=reg33*reg12;
    reg21=var_inter[1]*vectors[0][indices[3]+1]; reg34=reg27-reg34; reg27=reg2*reg0; reg33=reg1*vectors[0][indices[3]+0]; reg37=reg23-reg37;
    reg0=reg3*reg0; reg20=reg20/reg34; reg32=reg32/reg34; reg33=reg37+reg33; reg19=reg19/reg34;
    reg31=reg31/reg34; reg22=reg5-reg22; reg5=1-(*f.m).resolution; reg23=reg2*reg27; reg12=reg26-reg12;
    reg26=pow(reg3,2); reg21=reg35-reg21; reg29=reg17+reg29; reg3=reg3*reg0; reg2=pow(reg2,2);
    reg17=reg22*reg32; reg35=reg19*reg33; reg26=reg2-reg26; reg3=reg23-reg3; reg12=reg12*reg5;
    reg2=(*f.m).resolution*(*f.m).alpha; reg23=reg20*reg29; reg36=reg31*reg21; reg33=reg20*reg33; reg21=reg32*reg21;
    reg29=reg19*reg29; reg23=reg36-reg23; reg17=reg35-reg17; reg12=reg2+reg12; reg26=reg26/reg3;
    reg0=reg0/reg3; reg3=reg27/reg3; reg22=reg22*reg31; reg10=reg10*reg5; reg7=reg7*reg5;
    reg2=(*f.m).resolution*reg26; reg5=reg4*reg5; reg33=reg22-reg33; reg4=(*f.m).resolution*reg3; reg22=(*f.m).resolution*reg0;
    reg12=(*f.m).deltaT*reg12; reg21=reg29-reg21; reg17=reg23+reg17; reg23=var_inter[1]*reg32; reg27=reg1*reg20;
    reg29=var_inter[1]*reg31; reg33=reg33-reg12; reg17=0.5*reg17; reg35=reg1*reg19; reg36=reg8*reg31;
    reg21=reg21-reg12; reg2=reg10+reg2; reg5=reg22+reg5; reg10=var_inter[0]*reg20; reg7=reg4+reg7;
    reg4=var_inter[0]*reg19; reg22=reg8*reg32; reg37=reg7*reg33; reg38=reg35+reg23; T reg39=reg27+reg29;
    T reg40=reg27-reg36; T reg41=reg5*reg21; reg17=reg2*reg17; T reg42=reg4-reg23; reg33=reg5*reg33;
    reg21=reg7*reg21; T reg43=reg36+reg10; T reg44=reg22+reg4; T reg45=reg29-reg10; T reg46=reg22-reg35;
    T reg47=0.5*reg46; T reg48=0.5*reg42; reg21=reg33+reg21; reg41=reg37+reg41; reg17=2*reg17;
    reg33=0.5*reg44; reg37=0.5*reg40; T reg49=0.5*reg45; T reg50=0.5*reg38; T reg51=0.5*reg39;
    T reg52=0.5*reg43; T reg53=reg33*reg17; T reg54=reg1*var_inter[1]; T reg55=reg43*reg41; T reg56=reg38*reg21;
    T reg57=reg17*reg51; T reg58=var_inter[0]*var_inter[1]; T reg59=reg8*var_inter[0]; T reg60=reg37*reg17; T reg61=reg46*reg21;
    T reg62=reg47*reg17; T reg63=reg45*reg41; T reg64=reg8*reg1; T reg65=reg48*reg17; T reg66=reg50*reg17;
    T reg67=reg42*reg21; T reg68=reg49*reg17; T reg69=reg40*reg41; T reg70=reg52*reg17; T reg71=reg39*reg41;
    T reg72=reg44*reg21; T reg73=reg59*elem.f_vol_e[0]; reg65=reg63+reg65; reg62=reg69+reg62; reg63=reg64*elem.f_vol_e[0];
    reg69=reg58*elem.f_vol_e[0]; T reg74=reg54*elem.f_vol_e[1]; reg56=reg56-reg57; reg55=reg55-reg53; T reg75=reg59*elem.f_vol_e[1];
    reg60=reg61+reg60; reg61=reg64*elem.f_vol_e[1]; T reg76=reg54*elem.f_vol_e[0]; reg68=reg67+reg68; reg67=reg58*elem.f_vol_e[1];
    reg70=reg70-reg72; reg66=reg66-reg71; reg56=reg56-reg74; reg66=reg66-reg76; reg62=reg62-reg63;
    reg68=reg68-reg67; reg60=reg60-reg61; reg65=reg65-reg69; reg55=reg55-reg73; reg70=reg70-reg75;
    reg56=reg34*reg56; reg62=reg34*reg62; reg60=reg34*reg60; reg66=reg34*reg66; reg55=reg34*reg55;
    reg70=reg34*reg70; reg65=reg34*reg65; reg68=reg34*reg68; T vec_2=ponderation*reg55; T vec_3=ponderation*reg70;
    T vec_6=ponderation*reg66; T vec_1=ponderation*reg60; T vec_4=ponderation*reg65; T vec_0=ponderation*reg62; T vec_5=ponderation*reg68;
    T vec_7=ponderation*reg56;
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
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
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
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
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
    pthread_mutex_lock( &( f.mutex_assemble_matrix ) );
    pthread_mutex_unlock( &( f.mutex_assemble_matrix ) );
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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
      const Formulation<TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,T,wont_add_nz> &f,
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

