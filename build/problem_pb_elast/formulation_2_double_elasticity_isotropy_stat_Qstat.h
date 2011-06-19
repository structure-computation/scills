
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
    vecs[0][indice+1]=node.dep[1]; vecs[1][indice+1]=node.dep[1]; vecs[2][indice+1]=node.dep[1]; vecs[3][indice+1]=node.dep[1]; vecs[4][indice+1]=node.dep[1];
    vecs[0][indice+0]=node.dep[0]; vecs[1][indice+0]=node.dep[0]; vecs[2][indice+0]=node.dep[0]; vecs[3][indice+0]=node.dep[0]; vecs[4][indice+0]=node.dep[0];
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg3*reg2; T reg6=reg1*reg4; reg1=reg3*reg1; reg2=reg2*reg4;
    T reg7=reg3*reg2; reg2=reg2*reg4; T reg8=reg3*reg5; T reg9=reg6*reg4; reg6=reg3*reg6;
    T reg10=reg3*reg1; reg6=reg6+reg10; reg7=reg8+reg7; reg2=reg2-reg8; reg1=reg1*reg4;
    reg5=reg5*reg4; reg9=reg9-reg10; T reg11=reg3*reg7; T reg12=reg1+reg10; reg9=reg9*reg4;
    reg6=reg3*reg6; reg5=reg8+reg5; reg8=reg2*reg4; T reg13=reg3*reg5; reg12=reg12*reg3;
    reg6=reg9-reg6; reg11=reg8-reg11; reg13=reg11-reg13; reg12=reg6-reg12; reg2=reg2/reg13;
    reg7=reg7/reg13; reg12=reg12/reg13; reg6=reg2*reg12; reg8=reg12*reg7; reg9=(*f.m).alpha*reg7;
    reg11=elem.pos(2)[0]-elem.pos(0)[0]; T reg14=elem.pos(2)[1]-elem.pos(0)[1]; T reg15=reg8*reg7; reg13=reg5/reg13; reg5=(*f.m).alpha*reg2;
    T reg16=reg2*reg6; T reg17=elem.pos(1)[0]-elem.pos(0)[0]; T reg18=elem.pos(1)[1]-elem.pos(0)[1]; T reg19=reg18*reg11; reg13=(*f.m).alpha*reg13;
    reg5=reg9+reg5; reg15=reg16-reg15; reg9=reg14*reg17; reg8=reg8/reg15; reg13=reg5+reg13;
    reg19=reg9-reg19; reg15=reg6/reg15; reg8=reg8*reg13; reg18=reg18/reg19; reg17=reg17/reg19;
    reg5=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; reg15=reg13*reg15; reg6=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg14=reg14/reg19; reg11=reg11/reg19;
    reg9=reg0*reg4; reg13=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg16=reg3*reg0; T reg20=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; T reg21=1-(*f.m).resolution;
    T reg22=reg6*reg17; T reg23=reg16*reg3; reg8=reg15-reg8; reg15=reg14*reg5; T reg24=reg18*reg20;
    T reg25=reg13*reg11; T reg26=reg9*reg4; T reg27=PNODE(2).dep[1]-PNODE(0).dep[1]; reg25=reg22-reg25; elem.epsilon[0][1]=reg25;
    reg22=PNODE(2).dep[0]-PNODE(0).dep[0]; T reg28=pow(reg4,2); T reg29=(*f.m).alpha*(*f.m).resolution; reg8=reg8*reg21; reg23=reg26-reg23;
    reg26=(*f.m).alpha*(*f.m).deltaT; T reg30=pow(reg3,2); T reg31=PNODE(1).dep[1]-PNODE(0).dep[1]; T reg32=PNODE(1).dep[0]-PNODE(0).dep[0]; reg24=reg15-reg24;
    elem.epsilon[0][0]=reg24; reg15=reg22*reg17; reg22=reg22*reg18; reg16=reg16/reg23; reg8=reg29+reg8;
    reg29=reg18*reg27; T reg33=reg32*reg11; T reg34=reg11*reg31; reg30=reg28-reg30; reg28=reg25-reg26;
    T reg35=reg24-reg26; reg27=reg17*reg27; reg32=reg32*reg14; reg9=reg9/reg23; reg31=reg14*reg31;
    T reg36=reg28*reg16; T reg37=reg9*reg35; reg22=reg32-reg22; reg8=(*f.m).deltaT*reg8; reg7=reg21*reg7;
    reg28=reg9*reg28; reg35=reg35*reg16; reg34=reg27-reg34; reg29=reg31-reg29; reg2=reg2*reg21;
    reg23=reg30/reg23; reg27=(*f.m).resolution*reg9; reg30=(*f.m).resolution*reg16; reg33=reg15-reg33; reg15=reg34-reg8;
    reg33=reg29+reg33; reg13=reg13*reg14; reg6=reg6*reg18; reg21=reg12*reg21; reg28=reg35+reg28;
    reg12=reg22-reg8; reg36=reg37+reg36; reg20=reg17*reg20; reg5=reg11*reg5; reg29=(*f.m).resolution*reg23;
    reg30=reg7+reg30; reg27=reg2+reg27; reg2=reg12*reg30; reg21=reg29+reg21; reg33=0.5*reg33;
    reg6=reg13-reg6; reg7=reg36+reg28; reg13=reg12*reg27; reg5=reg20-reg5; reg20=reg15*reg30;
    reg29=reg15*reg27; reg7=reg7/3; reg20=reg13+reg20; reg13=reg33*reg21; reg5=reg6+reg5;
    reg29=reg2+reg29; reg13=2*reg13; reg28=reg28-reg7; reg5=0.5*reg5; elem.epsilon[0][2]=reg5;
    reg15=reg29*reg15; reg12=reg20*reg12; reg36=reg36-reg7; reg2=reg5*reg23; reg6=reg33*reg13;
    reg15=reg12+reg15; reg28=pow(reg28,2); reg36=pow(reg36,2); reg12=2*reg2; reg6=reg15+reg6;
    reg7=pow(reg7,2); reg28=reg36+reg28; reg7=reg28+reg7; reg2=reg12*reg2; reg6=reg19*reg6;
    reg12=0.16666666666666665741*reg6; reg24=reg24-reg8; reg7=reg2+reg7; reg6=0.33333333333333331483*reg6; reg8=reg25-reg8;
    reg2=reg27*reg8; reg15=reg24*reg30; reg8=reg30*reg8; reg27=reg24*reg27; reg7=1.5*reg7;
    reg12=reg6+reg12; elem.sigma[0][0]=reg27+reg8; elem.sigma[0][2]=reg5*reg21; elem.sigma_von_mises=pow(reg7,0.5); elem.sigma[0][1]=reg15+reg2;
    elem.ener=reg12/2;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg2*reg4; T reg6=reg1*reg3; reg1=reg2*reg1; reg4=reg4*reg3;
    T reg7=reg4*reg3; T reg8=reg2*reg5; T reg9=reg6*reg3; reg4=reg2*reg4; reg6=reg2*reg6;
    T reg10=reg2*reg1; reg6=reg6+reg10; reg5=reg5*reg3; reg1=reg1*reg3; reg9=reg9-reg10;
    reg7=reg7-reg8; reg4=reg8+reg4; reg5=reg8+reg5; reg8=reg2*reg4; T reg11=reg7*reg3;
    T reg12=reg1+reg10; reg6=reg2*reg6; reg9=reg9*reg3; reg6=reg9-reg6; reg9=reg2*reg5;
    reg12=reg12*reg2; reg8=reg11-reg8; reg12=reg6-reg12; reg9=reg8-reg9; reg7=reg7/reg9;
    reg12=reg12/reg9; reg4=reg4/reg9; reg6=reg7*reg12; reg8=reg12*reg4; reg11=(*f.m).alpha*reg7;
    T reg13=reg7*reg6; T reg14=reg8*reg4; reg9=reg5/reg9; reg5=(*f.m).alpha*reg4; reg9=(*f.m).alpha*reg9;
    reg14=reg13-reg14; reg11=reg5+reg11; reg6=reg6/reg14; reg5=reg0*reg3; reg9=reg11+reg9;
    reg0=reg2*reg0; reg14=reg8/reg14; reg6=reg9*reg6; reg8=elem.pos(2)[1]-elem.pos(0)[1]; reg11=elem.pos(2)[0]-elem.pos(0)[0];
    reg13=elem.pos(1)[0]-elem.pos(0)[0]; T reg15=reg0*reg2; reg2=pow(reg2,2); reg9=reg14*reg9; reg14=elem.pos(1)[1]-elem.pos(0)[1];
    T reg16=reg5*reg3; reg3=pow(reg3,2); T reg17=reg14*reg11; reg2=reg3-reg2; reg3=1-(*f.m).resolution;
    reg15=reg16-reg15; reg9=reg6-reg9; reg6=reg8*reg13; reg0=reg0/reg15; reg17=reg6-reg17;
    reg9=reg9*reg3; reg2=reg2/reg15; reg15=reg5/reg15; reg5=(*f.m).alpha*(*f.m).resolution; reg8=reg8/reg17;
    reg12=reg12*reg3; reg11=reg11/reg17; reg14=reg14/reg17; reg4=reg3*reg4; reg6=(*f.m).resolution*reg15;
    reg9=reg5+reg9; reg5=(*f.m).resolution*reg2; reg3=reg7*reg3; reg7=(*f.m).resolution*reg0; reg13=reg13/reg17;
    reg7=reg4+reg7; reg12=reg5+reg12; reg4=0.5*reg8; reg5=0.5*reg14; reg6=reg3+reg6;
    reg3=reg11-reg13; reg16=0.5*reg13; reg9=(*f.m).deltaT*reg9; T reg18=0.5*reg11; T reg19=reg14-reg8;
    T reg20=reg18*reg12; T reg21=reg6*reg9; T reg22=reg7*reg9; T reg23=reg16*reg12; T reg24=reg4*reg12;
    T reg25=0.5*reg3; T reg26=0.5*reg19; T reg27=reg5*reg12; T reg28=reg14*reg6; reg23=2*reg23;
    T reg29=reg13*reg6; T reg30=reg11*reg6; T reg31=reg8*reg6; T reg32=reg11*reg7; T reg33=2*reg20;
    T reg34=reg13*reg7; T reg35=reg21+reg22; T reg36=2*reg27; reg24=2*reg24; T reg37=1-var_inter[0];
    T reg38=reg14*reg7; T reg39=reg26*reg12; T reg40=reg8*reg7; T reg41=reg12*reg25; reg39=2*reg39;
    T reg42=reg35*reg11; T reg43=reg19*reg6; T reg44=reg5*reg24; T reg45=reg5*reg23; T reg46=reg38*reg13;
    T reg47=reg3*reg6; T reg48=reg14*reg31; T reg49=reg16*reg33; T reg50=reg33*reg4; T reg51=reg11*reg40;
    T reg52=reg36*reg16; reg37=reg37-var_inter[1]; T reg53=reg34*reg14; T reg54=reg11*reg29; T reg55=reg36*reg4;
    T reg56=reg13*reg30; T reg57=var_inter[0]*elem.f_vol_e[1]; T reg58=reg35*reg14; T reg59=reg8*reg32; T reg60=reg18*reg24;
    T reg61=var_inter[1]*elem.f_vol_e[0]; T reg62=reg19*reg7; reg41=2*reg41; T reg63=reg3*reg7; T reg64=reg18*reg23;
    T reg65=reg8*reg28; T reg66=reg47*reg13; T reg67=reg5*reg39; T reg68=reg18*reg41; T reg69=reg3*reg62;
    T reg70=var_inter[0]*elem.f_vol_e[0]; T reg71=reg33*reg25; T reg72=reg26*reg41; T reg73=reg36*reg25; T reg74=reg3*reg30;
    T reg75=reg3*reg40; T reg76=reg19*reg31; T reg77=reg26*reg24; reg53=reg52+reg53; T reg78=reg38*reg11;
    T reg79=reg18*reg33; reg31=reg8*reg31; T reg80=reg18*reg39; T reg81=reg8*reg63; reg60=reg59+reg60;
    T reg82=reg37*elem.f_vol_e[0]; T reg83=reg41*reg25; T reg84=reg28*reg14; T reg85=reg32*reg14; T reg86=reg16*reg24;
    T reg87=reg19*reg63; T reg88=reg39*reg25; T reg89=reg19*reg35; T reg90=reg13*reg62; T reg91=reg35*reg8;
    T reg92=reg3*reg35; reg62=reg11*reg62; T reg93=reg4*reg41; T reg94=reg36*reg18; T reg95=reg8*reg34;
    reg48=reg49+reg48; reg63=reg14*reg63; T reg96=reg16*reg39; T reg97=reg5*reg41; T reg98=reg11*reg47;
    T reg99=reg4*reg39; T reg100=reg19*reg32; T reg101=reg19*reg28; T reg102=reg24*reg25; T reg103=reg37*elem.f_vol_e[1];
    T reg104=reg5*reg36; T reg105=reg13*reg29; reg45=reg46+reg45; T reg106=reg35*reg13; T reg107=reg58-reg61;
    T reg108=reg5*reg33; reg40=reg13*reg40; T reg109=reg14*reg43; reg41=reg16*reg41; reg54=reg54+reg55;
    T reg110=reg3*reg38; T reg111=reg26*reg23; reg34=reg19*reg34; T reg112=reg23*reg25; reg64=reg65+reg64;
    T reg113=reg42-reg57; reg44=reg44+reg56; T reg114=reg23*reg4; T reg115=reg11*reg30; reg24=reg24*reg4;
    reg47=reg3*reg47; reg39=reg26*reg39; reg29=reg3*reg29; T reg116=reg26*reg36; T reg117=reg26*reg33;
    reg51=reg50+reg51; T reg118=reg19*reg43; reg43=reg8*reg43; reg23=reg16*reg23; T reg119=var_inter[1]*elem.f_vol_e[1];
    reg86=reg86+reg85; T reg120=reg17*reg64; reg105=reg105+reg104; reg47=reg39+reg47; reg88=reg87+reg88;
    reg39=reg17*reg60; reg87=reg17*reg45; reg112=reg112-reg101; reg68=reg43-reg68; reg31=reg31+reg79;
    reg29=reg29-reg116; reg80=reg81-reg80; reg77=reg77-reg74; reg43=reg106+reg119; reg109=reg41-reg109;
    reg111=reg111-reg110; reg107=reg17*reg107; reg41=reg17*reg54; reg114=reg114+reg78; reg40=reg40+reg108;
    reg81=reg17*reg51; reg23=reg23+reg84; reg67=reg66-reg67; reg66=reg17*reg44; reg62=reg93-reg62;
    reg34=reg34-reg73; reg97=reg90-reg97; reg90=reg103+reg92; reg93=reg70+reg91; reg95=reg95+reg94;
    reg76=reg76-reg71; T reg121=reg82+reg89; T reg122=reg17*reg53; reg24=reg24+reg115; reg75=reg75-reg117;
    reg102=reg102-reg100; reg113=reg17*reg113; reg63=reg96-reg63; reg96=reg17*reg48; reg98=reg99-reg98;
    reg69=reg72+reg69; reg83=reg118+reg83; reg24=reg17*reg24; reg34=reg17*reg34; reg72=ponderation*reg81;
    reg99=ponderation*reg41; reg111=reg17*reg111; reg29=reg17*reg29; reg69=reg17*reg69; reg114=reg17*reg114;
    reg118=ponderation*reg66; reg97=reg17*reg97; reg77=reg17*reg77; T reg123=ponderation*reg122; reg31=reg17*reg31;
    reg112=reg17*reg112; T reg124=reg17*reg93; T reg125=ponderation*reg39; reg113=ponderation*reg113; reg83=reg17*reg83;
    T reg126=ponderation*reg120; reg75=reg17*reg75; reg63=reg17*reg63; reg98=reg17*reg98; T reg127=reg121*reg17;
    reg95=reg17*reg95; reg62=reg17*reg62; reg76=reg17*reg76; T reg128=reg17*reg90; reg102=reg17*reg102;
    T reg129=ponderation*reg96; reg105=reg17*reg105; reg86=reg17*reg86; T reg130=ponderation*reg87; reg47=reg17*reg47;
    reg68=reg17*reg68; T reg131=reg17*reg43; reg80=reg17*reg80; reg109=reg17*reg109; reg88=reg17*reg88;
    reg107=ponderation*reg107; reg23=reg17*reg23; reg40=reg17*reg40; reg67=reg17*reg67; T tmp_0_3=ponderation*reg102;
    T tmp_4_0=ponderation*reg109; reg102=ponderation*reg128; sollicitation[indices[0]+1]+=reg102; T tmp_4_5=-reg123; T tmp_5_3=-reg118;
    T tmp_5_1=ponderation*reg67; reg67=ponderation*reg127; sollicitation[indices[0]+0]+=reg67; T tmp_5_0=ponderation*reg97; T tmp_0_2=ponderation*reg76;
    T tmp_2_5=ponderation*reg95; T tmp_1_5=ponderation*reg29; T tmp_3_3=ponderation*reg24; T tmp_3_0=ponderation*reg62; T tmp_0_5=ponderation*reg34;
    T tmp_3_1=ponderation*reg98; T tmp_2_4=-reg126; T tmp_4_4=ponderation*reg23; T tmp_3_5=-reg99; T tmp_4_1=ponderation*reg63;
    T tmp_2_1=ponderation*reg80; T tmp_1_3=ponderation*reg77; T tmp_1_1=ponderation*reg47; reg23=ponderation*reg131; sollicitation[indices[2]+1]+=reg23;
    T tmp_2_0=ponderation*reg68; T tmp_2_2=ponderation*reg31; T tmp_0_4=ponderation*reg112; sollicitation[indices[2]+0]+=-reg107; T tmp_3_4=ponderation*reg114;
    T tmp_5_4=-reg130; T tmp_1_2=ponderation*reg75; T tmp_3_2=-reg72; T tmp_4_3=ponderation*reg86; T tmp_5_2=ponderation*reg40;
    T tmp_2_3=-reg125; T tmp_5_5=ponderation*reg105; T tmp_4_2=-reg129; T tmp_0_0=ponderation*reg83; sollicitation[indices[1]+1]+=-reg113;
    T tmp_0_1=ponderation*reg88; T tmp_1_0=ponderation*reg69; T tmp_1_4=ponderation*reg111; reg24=ponderation*reg124; sollicitation[indices[1]+0]+=reg24;
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
      const Number<true> &assemble_vec,
      const Element<Triangle,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg2*reg4; T reg6=reg1*reg3; reg1=reg2*reg1; reg4=reg4*reg3;
    T reg7=reg4*reg3; T reg8=reg2*reg5; T reg9=reg6*reg3; reg4=reg2*reg4; reg6=reg2*reg6;
    T reg10=reg2*reg1; reg6=reg6+reg10; reg5=reg5*reg3; reg1=reg1*reg3; reg9=reg9-reg10;
    reg7=reg7-reg8; reg4=reg8+reg4; reg5=reg8+reg5; reg8=reg2*reg4; T reg11=reg7*reg3;
    T reg12=reg1+reg10; reg6=reg2*reg6; reg9=reg9*reg3; reg6=reg9-reg6; reg9=reg2*reg5;
    reg12=reg12*reg2; reg8=reg11-reg8; reg12=reg6-reg12; reg9=reg8-reg9; reg7=reg7/reg9;
    reg12=reg12/reg9; reg4=reg4/reg9; reg6=reg7*reg12; reg8=reg12*reg4; reg11=(*f.m).alpha*reg7;
    T reg13=reg7*reg6; T reg14=reg8*reg4; reg9=reg5/reg9; reg5=(*f.m).alpha*reg4; reg11=reg5+reg11;
    reg9=(*f.m).alpha*reg9; reg14=reg13-reg14; reg5=reg0*reg3; reg9=reg11+reg9; reg6=reg6/reg14;
    reg0=reg2*reg0; reg14=reg8/reg14; reg8=pow(reg2,2); reg6=reg9*reg6; reg11=elem.pos(2)[0]-elem.pos(0)[0];
    reg13=elem.pos(2)[1]-elem.pos(0)[1]; T reg15=elem.pos(1)[1]-elem.pos(0)[1]; T reg16=pow(reg3,2); reg9=reg14*reg9; reg3=reg5*reg3;
    reg2=reg0*reg2; reg14=elem.pos(1)[0]-elem.pos(0)[0]; reg8=reg16-reg8; reg2=reg3-reg2; reg3=1-(*f.m).resolution;
    reg9=reg6-reg9; reg6=reg13*reg14; reg16=reg15*reg11; reg0=reg0/reg2; reg16=reg6-reg16;
    reg5=reg5/reg2; reg6=(*f.m).alpha*(*f.m).resolution; reg9=reg9*reg3; reg2=reg8/reg2; reg12=reg12*reg3;
    reg8=(*f.m).resolution*reg5; reg15=reg15/reg16; reg13=reg13/reg16; reg14=reg14/reg16; reg4=reg3*reg4;
    T reg17=(*f.m).resolution*reg2; reg11=reg11/reg16; reg3=reg7*reg3; reg7=(*f.m).resolution*reg0; reg9=reg6+reg9;
    reg12=reg17+reg12; reg6=reg11-reg14; reg7=reg4+reg7; reg8=reg3+reg8; reg3=0.5*reg13;
    reg4=0.5*reg14; reg17=0.5*reg15; T reg18=reg15-reg13; reg9=(*f.m).deltaT*reg9; T reg19=reg8*reg9;
    T reg20=reg7*reg9; T reg21=0.5*reg11; T reg22=reg4*reg12; T reg23=reg3*reg12; T reg24=0.5*reg6;
    T reg25=reg17*reg12; T reg26=0.5*reg18; T reg27=reg19+reg20; T reg28=reg26*reg12; T reg29=reg21*reg12;
    T reg30=reg14*reg7; T reg31=reg14*reg8; T reg32=2*reg25; T reg33=reg12*reg24; reg23=2*reg23;
    T reg34=reg11*reg7; reg22=2*reg22; T reg35=1-var_inter[0]; T reg36=reg15*reg8; T reg37=reg18*reg8;
    T reg38=reg30*reg15; T reg39=reg32*reg4; T reg40=reg13*reg7; reg35=reg35-var_inter[1]; T reg41=reg11*reg8;
    T reg42=var_inter[1]*elem.f_vol_e[0]; reg33=2*reg33; T reg43=reg13*reg8; T reg44=2*reg29; T reg45=reg6*reg7;
    T reg46=reg21*reg22; T reg47=reg13*reg36; reg28=2*reg28; T reg48=reg32*reg3; T reg49=reg27*reg11;
    T reg50=reg27*reg15; T reg51=reg11*reg31; T reg52=reg13*reg34; T reg53=reg21*reg23; T reg54=reg15*reg7;
    T reg55=var_inter[0]*elem.f_vol_e[1]; T reg56=reg6*reg8; T reg57=reg18*reg43; T reg58=reg54*reg11; T reg59=reg44*reg24;
    T reg60=var_inter[0]*elem.f_vol_e[0]; T reg61=reg18*reg45; T reg62=reg13*reg30; T reg63=reg32*reg21; reg46=reg47+reg46;
    T reg64=reg18*reg34; T reg65=reg6*reg27; reg51=reg51+reg48; T reg66=reg21*reg44; T reg67=reg13*reg43;
    T reg68=reg27*reg13; T reg69=reg28*reg24; T reg70=reg18*reg27; T reg71=reg49-reg55; T reg72=reg17*reg32;
    T reg73=reg36*reg15; T reg74=reg14*reg31; T reg75=reg32*reg24; T reg76=reg22*reg3; T reg77=reg18*reg36;
    T reg78=reg26*reg44; T reg79=reg26*reg22; T reg80=reg50-reg42; T reg81=var_inter[1]*elem.f_vol_e[1]; reg38=reg39+reg38;
    T reg82=reg18*reg37; T reg83=reg23*reg24; T reg84=reg35*elem.f_vol_e[0]; T reg85=reg27*reg14; T reg86=reg35*elem.f_vol_e[1];
    T reg87=reg26*reg32; reg31=reg6*reg31; T reg88=reg6*reg54; T reg89=reg6*reg56; T reg90=reg4*reg22;
    T reg91=reg23*reg3; T reg92=reg22*reg24; T reg93=reg33*reg24; T reg94=reg11*reg41; reg53=reg52+reg53;
    T reg95=reg6*reg41; T reg96=reg26*reg23; T reg97=reg26*reg28; reg30=reg18*reg30; T reg98=reg6*reg40;
    T reg99=reg16*reg38; T reg100=reg16*reg53; reg67=reg67+reg66; reg76=reg76+reg58; reg93=reg82+reg93;
    reg90=reg90+reg73; reg30=reg30-reg75; reg82=reg60+reg68; reg71=reg16*reg71; T reg101=reg16*reg46;
    reg91=reg91+reg94; reg92=reg92-reg77; reg96=reg96-reg95; reg79=reg79-reg88; T reg102=reg16*reg51;
    reg80=reg16*reg80; T reg103=reg85+reg81; reg74=reg74+reg72; reg83=reg83-reg64; reg57=reg57-reg59;
    reg89=reg97+reg89; reg31=reg31-reg87; reg62=reg62+reg63; reg69=reg61+reg69; reg61=reg86+reg65;
    reg97=reg84+reg70; reg98=reg98-reg78; T reg104=reg97*reg16; reg80=ponderation*reg80; reg67=reg16*reg67;
    T reg105=ponderation*reg102; reg79=reg16*reg79; T reg106=ponderation*reg100; reg62=reg16*reg62; reg96=reg16*reg96;
    reg92=reg16*reg92; reg91=reg16*reg91; reg93=reg16*reg93; reg30=reg16*reg30; T reg107=ponderation*reg101;
    T reg108=reg16*reg61; reg71=ponderation*reg71; T reg109=ponderation*reg99; T reg110=reg16*reg82; reg57=reg16*reg57;
    reg89=reg16*reg89; reg83=reg16*reg83; reg69=reg16*reg69; reg98=reg16*reg98; reg74=reg16*reg74;
    reg76=reg16*reg76; reg31=reg16*reg31; reg90=reg16*reg90; T reg111=reg16*reg103; T reg112=ponderation*reg110;
    sollicitation[indices[1]+0]+=reg112; T tmp_4_5=-reg109; T tmp_4_4=ponderation*reg90; T tmp_3_3=ponderation*reg91; reg90=ponderation*reg108;
    sollicitation[indices[0]+1]+=reg90; T tmp_0_5=ponderation*reg30; T tmp_1_1=ponderation*reg89; sollicitation[indices[1]+1]+=-reg71; T tmp_2_3=-reg106;
    T tmp_1_5=ponderation*reg31; T tmp_1_2=ponderation*reg98; T tmp_0_2=ponderation*reg57; T tmp_0_0=ponderation*reg93; reg30=ponderation*reg104;
    sollicitation[indices[0]+0]+=reg30; T tmp_0_3=ponderation*reg83; T tmp_0_4=ponderation*reg92; T tmp_2_5=ponderation*reg62; T tmp_1_3=ponderation*reg96;
    T tmp_5_5=ponderation*reg74; T tmp_1_4=ponderation*reg79; T tmp_3_4=ponderation*reg76; T tmp_2_2=ponderation*reg67; T tmp_0_1=ponderation*reg69;
    T tmp_3_5=-reg105; T tmp_2_4=-reg107; reg31=ponderation*reg111; sollicitation[indices[2]+1]+=reg31; sollicitation[indices[2]+0]+=-reg80;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg2*reg4; T reg6=reg1*reg3; reg1=reg2*reg1; reg4=reg4*reg3;
    T reg7=reg2*reg1; T reg8=reg4*reg3; T reg9=reg2*reg6; T reg10=reg2*reg5; reg6=reg6*reg3;
    reg4=reg2*reg4; reg9=reg9+reg7; reg6=reg6-reg7; reg5=reg5*reg3; reg1=reg1*reg3;
    reg8=reg8-reg10; reg4=reg10+reg4; T reg11=reg2*reg4; reg5=reg10+reg5; reg6=reg6*reg3;
    reg10=reg1+reg7; T reg12=reg8*reg3; T reg13=reg0*reg3; reg9=reg2*reg9; reg0=reg2*reg0;
    T reg14=reg13*reg3; reg10=reg10*reg2; T reg15=elem.pos(2)[0]-elem.pos(0)[0]; T reg16=elem.pos(2)[1]-elem.pos(0)[1]; T reg17=reg2*reg5;
    reg9=reg6-reg9; reg6=pow(reg2,2); T reg18=elem.pos(1)[0]-elem.pos(0)[0]; reg11=reg12-reg11; reg2=reg0*reg2;
    reg3=pow(reg3,2); reg12=elem.pos(1)[1]-elem.pos(0)[1]; reg2=reg14-reg2; reg14=reg16*reg18; reg17=reg11-reg17;
    reg10=reg9-reg10; reg9=reg12*reg15; reg6=reg3-reg6; reg10=reg10/reg17; reg9=reg14-reg9;
    reg3=1-(*f.m).resolution; reg6=reg6/reg2; reg11=(*f.m).resolution*reg6; reg8=reg8/reg17; reg13=reg13/reg2;
    reg14=reg10*reg3; reg16=reg16/reg9; reg15=reg15/reg9; reg18=reg18/reg9; reg4=reg4/reg17;
    reg2=reg0/reg2; reg12=reg12/reg9; reg0=(*f.m).resolution*reg2; T reg19=0.5*reg16; T reg20=reg3*reg4;
    T reg21=reg15-reg18; T reg22=reg8*reg3; T reg23=0.5*reg18; T reg24=0.5*reg15; T reg25=(*f.m).resolution*reg13;
    T reg26=reg12-reg16; T reg27=0.5*reg12; reg14=reg11+reg14; reg11=0.5*reg21; T reg28=0.5*reg26;
    reg25=reg22+reg25; reg0=reg20+reg0; reg20=reg23*reg14; reg22=reg27*reg14; T reg29=reg19*reg14;
    T reg30=reg24*reg14; reg20=2*reg20; T reg31=reg16*reg25; T reg32=reg12*reg25; T reg33=reg15*reg0;
    reg29=2*reg29; T reg34=reg18*reg25; T reg35=reg15*reg25; T reg36=reg14*reg11; T reg37=2*reg30;
    T reg38=reg16*reg0; T reg39=reg18*reg0; T reg40=reg12*reg0; T reg41=reg28*reg14; T reg42=2*reg22;
    T reg43=reg39*reg12; T reg44=reg42*reg19; T reg45=reg26*reg25; reg41=2*reg41; T reg46=reg40*reg18;
    T reg47=reg16*reg32; T reg48=reg37*reg19; T reg49=reg15*reg38; T reg50=reg15*reg34; T reg51=reg27*reg20;
    T reg52=reg12*reg31; T reg53=reg21*reg25; T reg54=reg23*reg37; T reg55=reg42*reg23; T reg56=reg27*reg29;
    T reg57=reg24*reg20; T reg58=reg21*reg0; reg36=2*reg36; T reg59=reg24*reg29; T reg60=reg18*reg35;
    T reg61=reg26*reg0; T reg62=reg16*reg33; reg59=reg62+reg59; T reg63=reg28*reg29; T reg64=reg12*reg58;
    T reg65=reg21*reg35; T reg66=reg23*reg41; T reg67=reg26*reg58; T reg68=reg37*reg11; T reg69=reg26*reg31;
    T reg70=reg40*reg15; T reg71=reg36*reg11; reg31=reg16*reg31; T reg72=reg24*reg41; reg58=reg16*reg58;
    T reg73=reg28*reg36; T reg74=reg24*reg37; T reg75=reg42*reg11; T reg76=reg27*reg41; T reg77=reg32*reg12;
    T reg78=reg53*reg18; T reg79=reg33*reg12; T reg80=reg27*reg36; T reg81=reg18*reg61; T reg82=reg23*reg29;
    reg43=reg55+reg43; T reg83=reg41*reg11; T reg84=reg15*reg61; T reg85=reg19*reg36; T reg86=reg42*reg24;
    T reg87=reg16*reg39; reg61=reg21*reg61; reg52=reg54+reg52; T reg88=reg15*reg35; T reg89=reg21*reg40;
    reg51=reg46+reg51; T reg90=reg21*reg53; T reg91=reg28*reg41; reg49=reg48+reg49; T reg92=reg28*reg20;
    T reg93=reg20*reg19; T reg94=reg18*reg34; reg56=reg56+reg60; reg53=reg15*reg53; T reg95=reg27*reg42;
    reg41=reg19*reg41; reg39=reg26*reg39; T reg96=reg26*reg33; reg57=reg47+reg57; T reg97=reg29*reg11;
    T reg98=reg26*reg32; T reg99=reg20*reg11; T reg100=reg21*reg38; T reg101=reg23*reg36; T reg102=reg12*reg45;
    reg36=reg24*reg36; T reg103=reg16*reg45; reg20=reg23*reg20; reg50=reg50+reg44; reg38=reg18*reg38;
    reg29=reg29*reg19; T reg104=reg27*reg37; T reg105=reg28*reg37; reg45=reg26*reg45; T reg106=reg28*reg42;
    reg34=reg21*reg34; reg38=reg38+reg104; reg99=reg99-reg98; reg97=reg97-reg96; reg92=reg92-reg89;
    T reg107=reg9*reg59; reg36=reg103-reg36; reg102=reg101-reg102; reg101=reg9*reg51; reg76=reg78-reg76;
    reg20=reg20+reg77; reg94=reg94+reg95; reg78=reg9*reg50; reg63=reg63-reg65; reg72=reg58-reg72;
    reg93=reg93+reg70; reg31=reg31+reg74; reg34=reg34-reg106; reg39=reg39-reg75; reg90=reg91+reg90;
    reg64=reg66-reg64; reg29=reg29+reg88; reg80=reg81-reg80; reg83=reg67+reg83; reg87=reg87+reg86;
    reg58=reg9*reg43; reg53=reg41-reg53; reg84=reg85-reg84; reg41=reg9*reg56; reg66=reg9*reg49;
    reg61=reg73+reg61; reg67=reg9*reg52; reg71=reg45+reg71; reg45=reg9*reg57; reg100=reg100-reg105;
    reg69=reg69-reg68; reg82=reg82+reg79; reg61=reg9*reg61; reg73=ponderation*reg66; reg39=reg9*reg39;
    reg93=reg9*reg93; reg81=ponderation*reg78; reg92=reg9*reg92; reg29=reg9*reg29; reg80=reg9*reg80;
    reg63=reg9*reg63; reg85=ponderation*reg58; reg31=reg9*reg31; reg91=ponderation*reg41; reg71=reg9*reg71;
    reg99=reg9*reg99; reg103=ponderation*reg107; T reg108=ponderation*reg45; reg84=reg9*reg84; reg87=reg9*reg87;
    reg53=reg9*reg53; reg69=reg9*reg69; T reg109=ponderation*reg67; reg82=reg9*reg82; reg97=reg9*reg97;
    reg94=reg9*reg94; reg64=reg9*reg64; reg36=reg9*reg36; reg90=reg9*reg90; reg72=reg9*reg72;
    reg100=reg9*reg100; reg102=reg9*reg102; T reg110=ponderation*reg101; reg38=reg9*reg38; reg83=reg9*reg83;
    reg34=reg9*reg34; reg20=reg9*reg20; reg76=reg9*reg76; T tmp_1_2=ponderation*reg100; T tmp_4_5=-reg85;
    T tmp_3_1=ponderation*reg53; T tmp_3_3=ponderation*reg29; T tmp_1_5=ponderation*reg34; T tmp_2_4=-reg108; T tmp_1_0=ponderation*reg61;
    T tmp_5_0=ponderation*reg80; T tmp_4_1=ponderation*reg64; T tmp_0_5=ponderation*reg39; T tmp_4_0=ponderation*reg102; T tmp_0_1=ponderation*reg83;
    T tmp_3_2=-reg73; T tmp_1_1=ponderation*reg90; T tmp_3_5=-reg81; T tmp_4_4=ponderation*reg20; T tmp_1_4=ponderation*reg92;
    T tmp_3_4=ponderation*reg93; T tmp_5_1=ponderation*reg76; T tmp_5_2=ponderation*reg38; T tmp_1_3=ponderation*reg63; T tmp_5_4=-reg110;
    T tmp_2_1=ponderation*reg72; T tmp_2_0=ponderation*reg36; T tmp_2_2=ponderation*reg31; T tmp_0_4=ponderation*reg99; T tmp_5_5=ponderation*reg94;
    T tmp_0_3=ponderation*reg97; T tmp_4_3=ponderation*reg82; T tmp_4_2=-reg109; T tmp_0_0=ponderation*reg71; T tmp_0_2=ponderation*reg69;
    T tmp_2_3=-reg103; T tmp_5_3=-reg91; T tmp_3_0=ponderation*reg84; T tmp_2_5=ponderation*reg87;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg4*reg2; T reg6=reg1*reg3; reg1=reg4*reg1; reg2=reg2*reg3;
    T reg7=reg6*reg3; T reg8=reg4*reg5; T reg9=reg2*reg3; reg6=reg4*reg6; T reg10=reg4*reg1;
    reg2=reg4*reg2; reg5=reg5*reg3; reg1=reg1*reg3; reg6=reg6+reg10; reg7=reg7-reg10;
    reg9=reg9-reg8; reg2=reg8+reg2; T reg11=reg9*reg3; reg5=reg8+reg5; reg8=reg4*reg2;
    T reg12=reg4*reg0; T reg13=reg1+reg10; reg0=reg0*reg3; reg6=reg4*reg6; reg7=reg7*reg3;
    T reg14=pow(reg3,2); T reg15=elem.pos(1)[1]-elem.pos(0)[1]; T reg16=pow(reg4,2); T reg17=reg4*reg5; reg13=reg13*reg4;
    T reg18=elem.pos(1)[0]-elem.pos(0)[0]; T reg19=elem.pos(2)[0]-elem.pos(0)[0]; T reg20=elem.pos(2)[1]-elem.pos(0)[1]; reg8=reg11-reg8; reg6=reg7-reg6;
    reg3=reg0*reg3; reg4=reg12*reg4; reg7=reg20*reg18; reg11=reg15*reg19; reg13=reg6-reg13;
    reg16=reg14-reg16; reg17=reg8-reg17; reg4=reg3-reg4; reg11=reg7-reg11; reg3=1-(*f.m).resolution;
    reg16=reg16/reg4; reg13=reg13/reg17; reg2=reg2/reg17; reg20=reg20/reg11; reg9=reg9/reg17;
    reg12=reg12/reg4; reg4=reg0/reg4; reg0=reg13*reg3; reg19=reg19/reg11; reg6=(*f.m).resolution*reg16;
    reg18=reg18/reg11; reg15=reg15/reg11; reg7=0.5*reg20; reg8=0.5*reg15; reg14=reg9*reg3;
    T reg21=0.5*reg18; T reg22=reg3*reg2; T reg23=(*f.m).resolution*reg4; T reg24=reg15-reg20; T reg25=reg19-reg18;
    reg0=reg6+reg0; reg6=(*f.m).resolution*reg12; T reg26=reg21*reg0; T reg27=reg7*reg0; T reg28=0.5*reg19;
    T reg29=0.5*reg24; reg23=reg14+reg23; reg6=reg22+reg6; reg14=reg8*reg0; reg22=0.5*reg25;
    T reg30=reg18*reg23; T reg31=reg28*reg0; T reg32=reg0*reg22; T reg33=reg29*reg0; reg27=2*reg27;
    T reg34=reg19*reg6; T reg35=reg15*reg23; T reg36=reg18*reg6; T reg37=2*reg14; reg26=2*reg26;
    T reg38=reg20*reg6; T reg39=reg19*reg23; reg33=2*reg33; T reg40=reg20*reg23; T reg41=2*reg31;
    T reg42=reg25*reg6; reg32=2*reg32; T reg43=reg19*reg30; T reg44=reg15*reg6; T reg45=reg28*reg27;
    T reg46=reg20*reg34; T reg47=reg37*reg21; T reg48=reg36*reg15; T reg49=reg20*reg35; T reg50=reg28*reg26;
    T reg51=reg37*reg7; T reg52=reg25*reg23; T reg53=reg24*reg23; T reg54=reg37*reg28; T reg55=reg32*reg22;
    T reg56=reg24*reg53; T reg57=reg8*reg37; T reg58=reg18*reg30; reg50=reg49+reg50; T reg59=reg35*reg15;
    reg43=reg43+reg51; T reg60=reg20*reg40; T reg61=reg28*reg41; reg45=reg46+reg45; T reg62=reg25*reg44;
    T reg63=reg29*reg26; T reg64=reg37*reg22; T reg65=reg29*reg41; T reg66=reg24*reg36; T reg67=reg44*reg19;
    T reg68=reg29*reg33; T reg69=reg25*reg52; reg48=reg47+reg48; T reg70=reg26*reg7; T reg71=reg26*reg22;
    T reg72=reg24*reg34; T reg73=reg21*reg26; T reg74=reg27*reg7; T reg75=reg25*reg38; T reg76=reg19*reg39;
    T reg77=reg29*reg27; T reg78=reg25*reg39; T reg79=reg33*reg22; reg30=reg25*reg30; T reg80=reg24*reg42;
    T reg81=reg29*reg37; T reg82=reg24*reg40; reg36=reg20*reg36; T reg83=reg24*reg35; T reg84=reg41*reg22;
    T reg85=reg27*reg22; reg71=reg71-reg83; reg36=reg36+reg54; reg55=reg56+reg55; reg56=reg11*reg50;
    reg75=reg75-reg65; reg69=reg68+reg69; reg60=reg60+reg61; reg68=reg11*reg48; reg73=reg73+reg59;
    reg82=reg82-reg84; reg58=reg58+reg57; reg66=reg66-reg64; reg77=reg77-reg78; reg30=reg30-reg81;
    reg63=reg63-reg62; reg79=reg80+reg79; reg74=reg74+reg76; reg85=reg85-reg72; reg80=reg11*reg45;
    reg70=reg70+reg67; T reg86=reg11*reg43; reg60=reg11*reg60; reg66=reg11*reg66; T reg87=ponderation*reg80;
    reg74=reg11*reg74; reg55=reg11*reg55; T reg88=ponderation*reg68; reg79=reg11*reg79; reg85=reg11*reg85;
    reg36=reg11*reg36; reg58=reg11*reg58; reg75=reg11*reg75; T reg89=ponderation*reg86; reg63=reg11*reg63;
    reg30=reg11*reg30; reg77=reg11*reg77; reg82=reg11*reg82; reg73=reg11*reg73; reg69=reg11*reg69;
    reg71=reg11*reg71; T reg90=ponderation*reg56; reg70=reg11*reg70; T tmp_0_2=ponderation*reg82; T tmp_4_5=-reg88;
    T tmp_1_1=ponderation*reg69; T tmp_1_5=ponderation*reg30; T tmp_1_2=ponderation*reg75; T tmp_0_3=ponderation*reg85; T tmp_0_5=ponderation*reg66;
    T tmp_0_1=ponderation*reg79; T tmp_5_5=ponderation*reg58; T tmp_0_0=ponderation*reg55; T tmp_2_4=-reg90; T tmp_2_5=ponderation*reg36;
    T tmp_3_5=-reg89; T tmp_3_3=ponderation*reg74; T tmp_1_4=ponderation*reg63; T tmp_1_3=ponderation*reg77; T tmp_2_2=ponderation*reg60;
    T tmp_4_4=ponderation*reg73; T tmp_0_4=ponderation*reg71; T tmp_2_3=-reg87; T tmp_3_4=ponderation*reg70;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg1*reg4; T reg6=reg3*reg2; reg2=reg2*reg4; reg1=reg3*reg1;
    T reg7=reg3*reg1; T reg8=reg5*reg4; T reg9=reg3*reg2; reg5=reg3*reg5; T reg10=reg3*reg6;
    reg2=reg2*reg4; reg9=reg10+reg9; reg2=reg2-reg10; reg5=reg5+reg7; reg1=reg1*reg4;
    reg6=reg6*reg4; reg8=reg8-reg7; reg6=reg10+reg6; reg10=reg1+reg7; T reg11=reg2*reg4;
    reg8=reg8*reg4; reg5=reg3*reg5; T reg12=reg3*reg9; reg5=reg8-reg5; reg8=reg3*reg6;
    reg10=reg10*reg3; reg12=reg11-reg12; reg8=reg12-reg8; reg10=reg5-reg10; reg10=reg10/reg8;
    reg9=reg9/reg8; reg2=reg2/reg8; reg5=reg2*reg10; reg11=reg10*reg9; reg12=reg11*reg9;
    T reg13=(*f.m).alpha*reg2; reg8=reg6/reg8; reg6=reg2*reg5; T reg14=(*f.m).alpha*reg9; reg13=reg14+reg13;
    reg8=(*f.m).alpha*reg8; reg12=reg6-reg12; reg6=reg3*reg0; reg11=reg11/reg12; reg8=reg13+reg8;
    reg0=reg0*reg4; reg12=reg5/reg12; reg12=reg8*reg12; reg5=reg6*reg3; reg13=reg0*reg4;
    reg8=reg11*reg8; reg11=1-(*f.m).resolution; reg5=reg13-reg5; reg8=reg12-reg8; reg6=reg6/reg5;
    reg0=reg0/reg5; reg8=reg8*reg11; reg12=(*f.m).alpha*(*f.m).resolution; reg13=(*f.m).resolution*reg0; reg14=elem.pos(1)[1]-elem.pos(0)[1];
    T reg15=elem.pos(1)[0]-elem.pos(0)[0]; reg2=reg2*reg11; T reg16=(*f.m).resolution*reg6; T reg17=elem.pos(2)[1]-elem.pos(0)[1]; T reg18=elem.pos(2)[0]-elem.pos(0)[0];
    reg8=reg12+reg8; reg9=reg11*reg9; reg12=reg14*reg18; reg13=reg2+reg13; reg16=reg9+reg16;
    reg8=(*f.m).deltaT*reg8; reg2=reg17*reg15; reg9=reg13*reg8; T reg19=reg16*reg8; reg12=reg2-reg12;
    reg18=reg18/reg12; reg17=reg17/reg12; reg2=1-var_inter[0]; reg15=reg15/reg12; reg14=reg14/reg12;
    T reg20=reg9+reg19; T reg21=reg14-reg17; reg2=reg2-var_inter[1]; T reg22=reg20*reg18; T reg23=reg20*reg14;
    T reg24=var_inter[1]*elem.f_vol_e[0]; T reg25=var_inter[0]*elem.f_vol_e[1]; T reg26=reg18-reg15; T reg27=reg20*reg15; T reg28=reg2*elem.f_vol_e[0];
    T reg29=reg23-reg24; T reg30=reg22-reg25; T reg31=reg2*elem.f_vol_e[1]; T reg32=var_inter[1]*elem.f_vol_e[1]; T reg33=var_inter[0]*elem.f_vol_e[0];
    T reg34=reg21*reg20; T reg35=reg26*reg20; T reg36=reg20*reg17; reg30=reg12*reg30; T reg37=reg33+reg36;
    T reg38=reg27+reg32; T reg39=reg28+reg34; reg29=reg12*reg29; T reg40=reg31+reg35; reg30=ponderation*reg30;
    reg29=ponderation*reg29; T reg41=reg12*reg38; T reg42=reg12*reg37; T reg43=reg39*reg12; T reg44=reg12*reg40;
    sollicitation[indices[1]+1]+=-reg30; reg30=ponderation*reg42; sollicitation[indices[1]+0]+=reg30; sollicitation[indices[2]+0]+=-reg29; reg29=ponderation*reg41;
    sollicitation[indices[2]+1]+=reg29; T reg45=ponderation*reg43; sollicitation[indices[0]+0]+=reg45; T reg46=ponderation*reg44; sollicitation[indices[0]+1]+=reg46;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=reg0*reg1;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg3*reg4; T reg6=reg2*reg1; reg3=reg2*reg3; reg1=reg1*reg4;
    T reg7=reg2*reg5; T reg8=reg2*reg6; T reg9=reg2*reg1; reg5=reg5*reg4; T reg10=reg2*reg3;
    reg1=reg1*reg4; reg9=reg9+reg8; reg7=reg10+reg7; reg5=reg5-reg10; reg1=reg1-reg8;
    reg6=reg6*reg4; reg3=reg3*reg4; T reg11=reg2*reg7; reg3=reg10+reg3; reg1=reg1*reg4;
    reg10=reg6+reg8; T reg12=reg5*reg4; reg9=reg2*reg9; reg11=reg12-reg11; reg9=reg1-reg9;
    reg1=reg2*reg3; reg10=reg10*reg2; reg10=reg9-reg10; reg1=reg11-reg1; reg7=reg7/reg1;
    reg5=reg5/reg1; reg10=reg10/reg1; reg9=reg5*reg10; reg11=reg10*reg7; reg12=(*f.m).alpha*reg5;
    T reg13=reg5*reg9; T reg14=reg11*reg7; reg1=reg3/reg1; reg3=(*f.m).alpha*reg7; T reg15=elem.pos(1)[0]-elem.pos(0)[0];
    T reg16=elem.pos(2)[1]-elem.pos(0)[1]; T reg17=elem.pos(2)[0]-elem.pos(0)[0]; reg12=reg3+reg12; reg1=(*f.m).alpha*reg1; reg3=elem.pos(1)[1]-elem.pos(0)[1];
    reg14=reg13-reg14; reg13=reg16*reg15; reg9=reg9/reg14; reg1=reg12+reg1; reg12=reg3*reg17;
    reg14=reg11/reg14; reg14=reg14*reg1; reg11=reg2*reg0; reg9=reg1*reg9; reg0=reg0*reg4;
    reg12=reg13-reg12; reg1=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; reg13=1-(*f.m).resolution; reg15=reg15/reg12; reg3=reg3/reg12;
    T reg18=reg0*reg4; T reg19=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; T reg20=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg14=reg9-reg14; reg17=reg17/reg12;
    reg4=pow(reg4,2); reg9=reg11*reg2; reg16=reg16/reg12; reg2=pow(reg2,2); T reg21=vectors[0][indices[1]+1]-vectors[0][indices[0]+1];
    T reg22=reg21*reg16; T reg23=reg20*reg3; reg2=reg4-reg2; reg4=reg19*reg15; reg9=reg18-reg9;
    reg18=reg17*reg1; reg14=reg14*reg13; T reg24=(*f.m).alpha*(*f.m).resolution; reg21=reg21*reg17; reg20=reg20*reg15;
    reg19=reg19*reg3; reg23=reg22-reg23; reg18=reg4-reg18; reg1=reg16*reg1; reg2=reg2/reg9;
    reg0=reg0/reg9; reg14=reg24+reg14; reg9=reg11/reg9; reg4=(*f.m).resolution*reg0; reg5=reg5*reg13;
    reg19=reg1-reg19; reg18=reg23+reg18; reg7=reg13*reg7; reg14=(*f.m).deltaT*reg14; reg21=reg20-reg21;
    reg13=reg10*reg13; reg1=(*f.m).resolution*reg2; reg10=(*f.m).resolution*reg9; reg21=reg21-reg14; reg13=reg1+reg13;
    reg10=reg7+reg10; reg4=reg5+reg4; reg18=0.5*reg18; reg19=reg19-reg14; reg18=reg18*reg13;
    reg1=reg17-reg15; reg5=reg3-reg16; reg7=reg19*reg10; reg11=reg21*reg10; reg21=reg21*reg4;
    reg19=reg19*reg4; reg20=1-var_inter[0]; reg22=0.5*reg17; reg23=0.5*reg15; reg18=2*reg18;
    reg24=0.5*reg5; reg11=reg19+reg11; reg21=reg7+reg21; reg7=0.5*reg3; reg19=0.5*reg16;
    T reg25=0.5*reg1; T reg26=reg17*reg21; T reg27=reg5*reg11; T reg28=reg24*reg18; T reg29=reg18*reg7;
    T reg30=reg18*reg23; T reg31=reg18*reg25; T reg32=reg1*reg21; T reg33=reg16*reg11; T reg34=reg18*reg19;
    T reg35=reg18*reg22; reg20=reg20-var_inter[1]; T reg36=reg11*reg3; T reg37=reg21*reg15; reg33=reg33-reg35;
    T reg38=var_inter[0]*elem.f_vol_e[0]; reg31=reg27+reg31; reg27=reg20*elem.f_vol_e[0]; reg34=reg34-reg26; T reg39=var_inter[1]*elem.f_vol_e[1];
    reg37=reg37-reg29; reg32=reg28+reg32; reg28=reg20*elem.f_vol_e[1]; T reg40=var_inter[1]*elem.f_vol_e[0]; reg30=reg30-reg36;
    T reg41=var_inter[0]*elem.f_vol_e[1]; reg34=reg34-reg41; reg32=reg32-reg28; reg33=reg33-reg38; reg30=reg30-reg40;
    reg31=reg31-reg27; reg37=reg37-reg39; reg31=reg31*reg12; reg37=reg12*reg37; reg34=reg12*reg34;
    reg32=reg12*reg32; reg33=reg12*reg33; reg30=reg30*reg12; sollicitation[indices[0]+0]+=reg31*ponderation; sollicitation[indices[2]+1]+=ponderation*reg37;
    sollicitation[indices[0]+1]+=ponderation*reg32; sollicitation[indices[1]+0]+=ponderation*reg33; sollicitation[indices[2]+0]+=ponderation*reg30; sollicitation[indices[1]+1]+=ponderation*reg34;
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
    T reg0=vecs[1][indice+1]-vecs[0][indice+1]; T reg1=vecs[1][indice+0]-vecs[0][indice+0]; reg1=abs(reg1); reg0=abs(reg0); return max(reg0,reg1);
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
    T reg4=reg0*reg1; T reg5=reg2*reg4; T reg6=reg3*reg1; reg1=reg2*reg1; reg4=reg3*reg4;
    T reg7=reg3*reg1; T reg8=reg2*reg5; T reg9=reg3*reg4; reg5=reg3*reg5; reg1=reg2*reg1;
    T reg10=reg3*reg6; reg8=reg8-reg9; reg5=reg9+reg5; reg4=reg2*reg4; reg7=reg10+reg7;
    reg6=reg2*reg6; reg1=reg1-reg10; reg7=reg3*reg7; T reg11=reg10+reg6; T reg12=reg3*reg5;
    T reg13=reg2*reg8; reg4=reg9+reg4; reg1=reg2*reg1; reg7=reg1-reg7; reg1=reg3*reg4;
    reg12=reg13-reg12; reg11=reg3*reg11; reg1=reg12-reg1; reg11=reg7-reg11; reg8=reg8/reg1;
    reg5=reg5/reg1; reg11=reg11/reg1; reg7=0.5*elem.pos(0)[0]; reg9=0.21132486540518713447*elem.pos(0)[0]; reg12=0.21132486540518713447*elem.pos(1)[0];
    reg13=0.78867513459481286553*elem.pos(0)[0]; T reg14=0.5*elem.pos(1)[0]; T reg15=0.78867513459481286553*elem.pos(1)[1]; T reg16=0.21132486540518713447*elem.pos(0)[1]; T reg17=0.78867513459481286553*elem.pos(0)[1];
    T reg18=reg8*reg11; T reg19=reg5*reg11; T reg20=0.21132486540518713447*elem.pos(1)[1]; T reg21=0.78867513459481286553*elem.pos(1)[0]; T reg22=0.5*elem.pos(1)[1];
    T reg23=0.5*elem.pos(0)[1]; T reg24=reg20+reg17; reg1=reg4/reg1; reg4=reg22+reg23; T reg25=0.21132486540518713447*elem.pos(2)[1];
    T reg26=0.5*elem.pos(2)[1]; T reg27=0.21132486540518713447*elem.pos(2)[0]; T reg28=reg7+reg14; T reg29=reg5*(*f.m).alpha; T reg30=reg9+reg21;
    reg23=reg22-reg23; reg22=reg8*(*f.m).alpha; T reg31=reg5*reg19; reg20=reg20-reg16; T reg32=0.78867513459481286553*elem.pos(2)[0];
    T reg33=reg8*reg18; T reg34=0.5*elem.pos(2)[0]; reg7=reg14-reg7; reg9=reg12-reg9; reg12=reg12+reg13;
    reg16=reg15+reg16; reg14=0.78867513459481286553*elem.pos(2)[1]; reg31=reg33-reg31; reg29=reg22+reg29; reg1=reg1*(*f.m).alpha;
    reg28=reg34-reg28; reg13=reg21-reg13; reg23=reg26+reg23; reg17=reg15-reg17; reg15=0.5*elem.pos(3)[0];
    reg4=reg26-reg4; reg34=reg7+reg34; reg7=0.5*elem.pos(3)[1]; reg21=0.78867513459481286553*elem.pos(3)[0]; reg9=reg9+reg32;
    reg22=0.21132486540518713447*elem.pos(3)[0]; reg30=reg32-reg30; reg12=reg27-reg12; reg16=reg14-reg16; reg26=0.78867513459481286553*elem.pos(3)[1];
    reg20=reg14+reg20; reg24=reg25-reg24; reg14=0.21132486540518713447*elem.pos(3)[1]; reg9=reg9-reg21; reg32=0.21132486540518713447*PNODE(1).dep[1];
    reg27=reg13+reg27; reg13=0.5*vectors[0][indices[0]+0]; reg24=reg26+reg24; reg23=reg23-reg7; reg16=reg16+reg14;
    reg18=reg18/reg31; reg33=0.5*vectors[0][indices[0]+1]; T reg35=0.5*vectors[0][indices[1]+1]; reg28=reg15+reg28; T reg36=0.5*vectors[0][indices[1]+0];
    reg1=reg29+reg1; reg30=reg30+reg22; reg26=reg20-reg26; reg31=reg19/reg31; reg19=0.78867513459481286553*PNODE(0).dep[1];
    reg12=reg21+reg12; reg20=0.21132486540518713447*PNODE(0).dep[1]; reg21=0.78867513459481286553*PNODE(1).dep[0]; reg29=0.78867513459481286553*PNODE(1).dep[1]; T reg37=0.78867513459481286553*PNODE(0).dep[0];
    reg7=reg4+reg7; reg4=0.21132486540518713447*PNODE(0).dep[0]; T reg38=0.21132486540518713447*PNODE(1).dep[0]; reg15=reg34-reg15; reg25=reg17+reg25;
    reg17=reg36-reg13; reg34=0.5*vectors[0][indices[2]+0]; T reg39=0.5*vectors[0][indices[2]+1]; T reg40=reg33+reg35; reg33=reg35-reg33;
    reg13=reg36+reg13; reg35=reg7*reg15; reg36=reg9*reg16; T reg41=reg26*reg30; T reg42=reg9*reg24;
    T reg43=reg20+reg29; T reg44=0.78867513459481286553*PNODE(2).dep[1]; T reg45=reg26*reg12; T reg46=reg3*reg0; T reg47=reg2*reg0;
    T reg48=reg19+reg32; T reg49=0.21132486540518713447*PNODE(2).dep[1]; T reg50=reg4+reg21; reg18=reg18*reg1; reg1=reg31*reg1;
    reg31=reg23*reg28; reg20=reg32-reg20; reg22=reg27-reg22; reg14=reg25-reg14; reg4=reg38-reg4;
    reg25=0.78867513459481286553*PNODE(2).dep[0]; reg38=reg37+reg38; reg27=0.21132486540518713447*PNODE(2).dep[0]; reg48=reg49-reg48; reg32=reg30*reg14;
    reg17=reg17+reg34; T reg51=0.5*vectors[0][indices[3]+0]; reg13=reg34-reg13; reg45=reg42-reg45; reg34=0.5*vectors[0][indices[3]+1];
    reg42=reg16*reg22; reg19=reg29-reg19; reg37=reg21-reg37; reg38=reg27-reg38; reg31=reg35-reg31;
    reg4=reg4+reg25; reg21=0.78867513459481286553*PNODE(3).dep[0]; reg29=0.21132486540518713447*PNODE(3).dep[0]; reg35=0.78867513459481286553*PNODE(3).dep[1]; reg20=reg44+reg20;
    reg50=reg25-reg50; reg1=reg18-reg1; reg18=1-(*f.m).resolution; reg25=reg3*reg46; reg40=reg39-reg40;
    T reg52=0.21132486540518713447*PNODE(3).dep[1]; reg43=reg44-reg43; reg33=reg39+reg33; reg41=reg36-reg41; reg36=reg2*reg47;
    reg49=reg19+reg49; reg38=reg38+reg21; reg21=reg4-reg21; reg43=reg52+reg43; reg4=reg9/reg41;
    reg27=reg37+reg27; reg19=reg24/reg45; reg37=reg12/reg45; reg9=reg9/reg45; reg39=pow(reg3,2);
    reg50=reg29+reg50; reg25=reg36-reg25; reg36=pow(reg2,2); reg48=reg48+reg35; reg44=reg26/reg41;
    reg40=reg34+reg40; reg15=reg15/reg31; reg23=reg23/reg31; reg13=reg51+reg13; reg34=reg33-reg34;
    reg51=reg17-reg51; reg28=reg28/reg31; reg17=reg14*reg12; reg33=reg22*reg24; reg31=reg7/reg31;
    reg7=(*f.m).resolution*(*f.m).alpha; reg1=reg18*reg1; reg26=reg26/reg45; reg35=reg20-reg35; reg20=reg30/reg41;
    reg32=reg42-reg32; reg42=reg16/reg41; T reg53=reg43*reg44; T reg54=reg35*reg42; T reg55=reg20*reg21;
    T reg56=reg38*reg26; T reg57=reg40*reg15; T reg58=reg13*reg23; T reg59=reg4*reg50; reg44=reg50*reg44;
    reg20=reg35*reg20; T reg60=reg31*reg51; reg42=reg42*reg21; T reg61=reg34*reg28; T reg62=reg19*reg21;
    reg39=reg36-reg39; reg36=reg48*reg9; reg1=reg7+reg1; reg7=reg37*reg35; reg47=reg47/reg25;
    T reg63=reg22/reg32; reg46=reg46/reg25; reg52=reg49-reg52; reg30=reg30/reg32; reg4=reg43*reg4;
    reg49=reg14/reg32; reg16=reg16/reg32; reg9=reg9*reg38; reg29=reg27-reg29; reg21=reg37*reg21;
    reg26=reg48*reg26; reg17=reg33-reg17; reg35=reg19*reg35; reg5=reg5*reg18; reg58=reg60-reg58;
    elem.epsilon[0][0]=reg58; reg8=reg8*reg18; reg19=reg63*reg50; reg27=reg29*reg30; reg33=(*f.m).resolution*reg46;
    reg37=(*f.m).resolution*reg47; reg1=(*f.m).deltaT*reg1; reg20=reg4-reg20; reg14=reg14/reg17; reg24=reg24/reg17;
    reg61=reg57-reg61; elem.epsilon[0][1]=reg61; reg4=reg52*reg16; reg57=reg43*reg49; reg53=reg54-reg53;
    reg22=reg22/reg17; reg7=reg36-reg7; reg12=reg12/reg17; reg30=reg52*reg30; reg63=reg43*reg63;
    reg56=reg62-reg56; reg55=reg59-reg55; reg44=reg42-reg44; reg21=reg9-reg21; reg9=(*f.m).deltaT*(*f.m).alpha;
    reg26=reg35-reg26; reg50=reg49*reg50; reg16=reg29*reg16; reg25=reg39/reg25; reg37=reg8+reg37;
    reg57=reg4-reg57; reg4=reg22*reg48; reg20=reg20-reg1; reg8=reg12*reg52; reg18=reg11*reg18;
    reg5=reg33+reg5; reg48=reg14*reg48; reg50=reg16-reg50; reg52=reg24*reg52; reg22=reg22*reg38;
    reg12=reg29*reg12; reg11=reg61-reg9; reg16=(*f.m).resolution*reg25; reg26=reg21+reg26; reg21=reg58-reg9;
    reg27=reg19-reg27; reg56=reg56-reg1; reg30=reg63-reg30; reg38=reg14*reg38; reg24=reg29*reg24;
    reg44=reg44-reg1; reg7=reg7-reg1; reg53=reg55+reg53; reg30=reg30-reg1; reg8=reg4-reg8;
    reg4=reg56*reg5; reg14=reg37*reg44; reg26=0.5*reg26; reg19=reg5*reg44; reg29=reg47*reg21;
    reg48=reg52-reg48; reg33=reg37*reg20; reg35=reg47*reg11; reg53=0.5*reg53; reg57=reg27+reg57;
    reg38=reg24-reg38; reg18=reg16+reg18; reg16=reg5*reg20; reg21=reg46*reg21; reg24=reg7*reg37;
    reg50=reg50-reg1; reg27=reg56*reg37; reg36=reg7*reg5; reg12=reg22-reg12; reg11=reg46*reg11;
    reg22=reg50*reg5; reg39=reg30*reg37; reg57=0.5*reg57; reg42=reg50*reg37; reg43=reg30*reg5;
    reg33=reg19+reg33; reg19=reg26*reg18; reg36=reg27+reg36; reg11=reg29+reg11; reg4=reg24+reg4;
    reg38=reg38-reg1; reg24=reg18*reg53; reg35=reg21+reg35; reg15=reg13*reg15; reg28=reg51*reg28;
    reg34=reg31*reg34; reg40=reg23*reg40; reg16=reg14+reg16; reg48=reg12+reg48; reg8=reg8-reg1;
    reg12=reg11+reg35; reg16=reg44*reg16; reg39=reg22+reg39; reg33=reg20*reg33; reg19=2*reg19;
    reg56=reg36*reg56; reg7=reg4*reg7; reg4=reg38*reg5; reg28=reg15-reg28; reg40=reg34-reg40;
    reg13=reg8*reg37; reg24=2*reg24; reg14=reg38*reg37; reg15=reg8*reg5; reg48=0.5*reg48;
    reg43=reg42+reg43; reg20=reg57*reg18; reg19=reg26*reg19; reg56=reg7+reg56; reg30=reg39*reg30;
    reg43=reg50*reg43; reg20=2*reg20; reg12=reg12/3; reg40=reg28+reg40; reg13=reg4+reg13;
    reg53=reg24*reg53; reg4=reg48*reg18; reg15=reg14+reg15; reg16=reg33+reg16; reg4=2*reg4;
    reg40=0.5*reg40; elem.epsilon[0][2]=reg40; reg16=reg53+reg16; reg43=reg30+reg43; reg15=reg38*reg15;
    reg20=reg57*reg20; reg11=reg11-reg12; reg35=reg35-reg12; reg19=reg56+reg19; reg13=reg8*reg13;
    reg11=pow(reg11,2); reg16=reg41*reg16; reg35=pow(reg35,2); reg4=reg48*reg4; reg45=reg19*reg45;
    reg15=reg13+reg15; reg20=reg43+reg20; reg7=reg40*reg25; reg45=0.25*reg45; reg16=0.25*reg16;
    reg35=reg11+reg35; reg12=pow(reg12,2); reg8=2*reg7; reg15=reg4+reg15; reg32=reg20*reg32;
    reg12=reg35+reg12; reg16=reg45+reg16; reg8=reg7*reg8; reg15=reg17*reg15; reg32=0.25*reg32;
    reg8=reg12+reg8; reg16=reg32+reg16; reg15=0.25*reg15; reg61=reg61-reg1; reg58=reg58-reg1;
    reg8=1.5*reg8; reg15=reg16+reg15; reg4=reg61*reg37; reg7=reg58*reg5; reg61=reg61*reg5;
    reg58=reg58*reg37; elem.sigma[0][2]=reg40*reg18; elem.sigma_von_mises=pow(reg8,0.5); elem.ener=reg15/2; elem.sigma[0][1]=reg7+reg4;
    elem.sigma[0][0]=reg58+reg61;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg4*reg2; reg2=reg3*reg2; T reg6=reg3*reg1; reg1=reg4*reg1;
    T reg7=reg4*reg5; T reg8=reg4*reg1; T reg9=reg3*reg6; T reg10=reg3*reg2; reg1=reg3*reg1;
    reg5=reg3*reg5; reg6=reg4*reg6; reg8=reg8-reg9; reg1=reg9+reg1; reg5=reg10+reg5;
    reg2=reg4*reg2; reg7=reg7-reg10; reg1=reg3*reg1; reg2=reg10+reg2; reg10=reg4*reg7;
    reg8=reg4*reg8; T reg11=reg9+reg6; T reg12=reg3*reg5; reg11=reg3*reg11; reg1=reg8-reg1;
    reg8=reg3*reg2; reg12=reg10-reg12; reg8=reg12-reg8; reg10=1-var_inter[1]; reg11=reg1-reg11;
    reg1=1-var_inter[0]; reg7=reg7/reg8; reg5=reg5/reg8; reg12=elem.pos(1)[0]*var_inter[0]; T reg13=reg1*elem.pos(0)[0];
    T reg14=reg10*elem.pos(1)[0]; T reg15=reg10*elem.pos(0)[0]; T reg16=elem.pos(1)[1]*var_inter[0]; T reg17=reg1*elem.pos(0)[1]; reg11=reg11/reg8;
    T reg18=reg10*elem.pos(0)[1]; T reg19=reg10*elem.pos(1)[1]; reg19=reg19-reg18; T reg20=elem.pos(2)[1]*var_inter[1]; T reg21=reg13+reg12;
    T reg22=elem.pos(2)[0]*var_inter[0]; T reg23=reg5*reg11; T reg24=reg7*reg11; T reg25=reg16+reg17; T reg26=elem.pos(2)[1]*var_inter[0];
    T reg27=elem.pos(2)[0]*var_inter[1]; reg14=reg14-reg15; reg19=reg20+reg19; reg20=elem.pos(3)[1]*var_inter[1]; T reg28=elem.pos(3)[1]*reg1;
    reg8=reg2/reg8; reg2=reg5*reg23; T reg29=elem.pos(3)[0]*reg1; reg22=reg22-reg21; T reg30=reg7*(*f.m).alpha;
    T reg31=elem.pos(3)[0]*var_inter[1]; reg26=reg26-reg25; reg27=reg14+reg27; reg14=reg5*(*f.m).alpha; T reg32=reg7*reg24;
    reg2=reg32-reg2; reg29=reg22+reg29; reg28=reg26+reg28; reg14=reg30+reg14; reg19=reg19-reg20;
    reg27=reg27-reg31; reg8=reg8*(*f.m).alpha; reg22=reg3*reg0; reg0=reg4*reg0; reg26=reg19*reg29;
    reg30=reg27*reg28; reg23=reg23/reg2; reg2=reg24/reg2; reg8=reg14+reg8; reg26=reg30-reg26;
    reg14=pow(reg4,2); reg4=reg4*reg0; reg24=pow(reg3,2); reg3=reg3*reg22; reg23=reg23*reg8;
    reg8=reg2*reg8; reg19=reg19/reg26; reg2=1-(*f.m).resolution; reg29=reg29/reg26; reg24=reg14-reg24;
    reg28=reg28/reg26; reg27=reg27/reg26; reg23=reg8-reg23; reg3=reg4-reg3; reg4=reg10*reg28;
    reg8=var_inter[1]*reg28; reg14=var_inter[1]*reg29; reg22=reg22/reg3; reg30=reg1*reg27; reg32=var_inter[0]*reg19;
    T reg33=reg1*reg19; T reg34=reg10*reg29; T reg35=var_inter[0]*reg27; reg0=reg0/reg3; reg23=reg2*reg23;
    T reg36=(*f.m).resolution*(*f.m).alpha; reg3=reg24/reg3; reg24=(*f.m).resolution*reg3; T reg37=reg4+reg32; T reg38=(*f.m).resolution*reg22;
    T reg39=reg34+reg35; reg5=reg5*reg2; reg11=reg11*reg2; reg2=reg7*reg2; reg7=(*f.m).resolution*reg0;
    T reg40=reg14+reg30; T reg41=reg33+reg8; reg23=reg36+reg23; reg23=(*f.m).deltaT*reg23; reg36=reg34-reg30;
    T reg42=reg8-reg32; T reg43=reg35-reg14; reg7=reg2+reg7; reg2=0.5*reg37; T reg44=reg33-reg4;
    T reg45=0.5*reg41; T reg46=0.5*reg39; reg5=reg38+reg5; reg38=0.5*reg40; reg11=reg24+reg11;
    reg24=reg11*reg38; T reg47=0.5*reg44; T reg48=reg11*reg45; T reg49=0.5*reg36; T reg50=0.5*reg42;
    T reg51=0.5*reg43; T reg52=reg11*reg2; T reg53=reg5*reg23; T reg54=reg7*reg23; T reg55=reg11*reg46;
    T reg56=reg7*reg39; T reg57=reg5*reg41; T reg58=reg7*reg40; T reg59=reg1*var_inter[1]; T reg60=reg7*reg41;
    reg24=2*reg24; T reg61=reg5*reg40; T reg62=2*reg48; T reg63=reg5*reg37; T reg64=reg11*reg49;
    T reg65=reg11*reg50; T reg66=reg11*reg47; T reg67=reg10*var_inter[0]; T reg68=reg7*reg37; T reg69=reg53+reg54;
    T reg70=2*reg55; T reg71=reg5*reg39; reg52=2*reg52; T reg72=reg11*reg51; T reg73=reg7*reg43;
    reg65=2*reg65; T reg74=reg5*reg42; T reg75=reg68*reg41; T reg76=reg70*reg38; T reg77=reg5*reg44;
    T reg78=reg7*reg36; T reg79=reg70*reg2; T reg80=reg24*reg46; T reg81=reg60*reg37; T reg82=reg52*reg46;
    T reg83=reg71*reg37; T reg84=reg59*elem.f_vol_e[0]; T reg85=reg61*reg41; T reg86=reg69*reg41; T reg87=reg7*reg44;
    T reg88=reg62*reg38; T reg89=reg63*reg39; reg64=2*reg64; T reg90=reg5*reg36; T reg91=reg5*reg43;
    T reg92=reg56*reg40; T reg93=reg52*reg45; reg66=2*reg66; T reg94=reg10*reg1; T reg95=var_inter[0]*var_inter[1];
    T reg96=reg69*reg39; reg72=2*reg72; T reg97=reg58*reg39; T reg98=reg7*reg42; T reg99=reg62*reg2;
    T reg100=reg67*elem.f_vol_e[1]; T reg101=reg57*reg40; T reg102=reg24*reg45; T reg103=reg95*elem.f_vol_e[0]; T reg104=reg69*reg37;
    T reg105=reg69*reg36; T reg106=reg69*reg44; T reg107=reg62*reg45; T reg108=reg65*reg50; T reg109=reg57*reg43;
    T reg110=reg24*reg50; T reg111=reg58*reg43; T reg112=reg62*reg50; T reg113=reg64*reg38; T reg114=reg87*reg41;
    T reg115=reg65*reg46; T reg116=reg91*reg44; T reg117=reg65*reg49; T reg118=reg91*reg41; T reg119=reg60*reg44;
    T reg120=reg67*elem.f_vol_e[0]; T reg121=reg94*elem.f_vol_e[1]; T reg122=reg24*reg49; T reg123=reg61*reg44; T reg124=reg94*elem.f_vol_e[0];
    T reg125=reg59*elem.f_vol_e[1]; T reg126=reg62*reg49; T reg127=reg77*reg36; T reg128=reg64*reg47; T reg129=reg95*elem.f_vol_e[1];
    T reg130=reg78*reg36; T reg131=reg66*reg47; T reg132=reg58*reg36; T reg133=reg62*reg47; T reg134=reg87*reg37;
    T reg135=reg64*reg46; T reg136=reg90*reg37; T reg137=reg66*reg46; T reg138=reg87*reg42; T reg139=reg64*reg51;
    T reg140=reg90*reg42; T reg141=reg66*reg51; T reg142=reg68*reg42; T reg143=reg70*reg51; T reg144=reg71*reg42;
    T reg145=reg52*reg51; T reg146=reg98*reg42; T reg147=reg72*reg51; T reg148=reg61*reg42; T reg149=reg24*reg51;
    T reg150=reg91*reg42; T reg151=reg65*reg51; T reg152=reg60*reg42; reg58=reg58*reg40; reg102=reg101+reg102;
    T reg153=reg65*reg45; T reg154=reg73*reg40; T reg155=reg72*reg45; T reg156=reg74*reg40; reg93=reg92+reg93;
    T reg157=reg63*reg36; T reg158=reg70*reg47; T reg159=reg56*reg36; T reg160=reg52*reg47; T reg161=reg60*reg41;
    T reg162=reg65*reg38; T reg163=reg74*reg36; T reg164=reg72*reg47; T reg165=reg64*reg50; T reg166=reg73*reg36;
    T reg167=reg65*reg47; T reg168=reg77*reg43; T reg169=reg57*reg36; T reg170=reg24*reg47; T reg171=reg62*reg51;
    T reg172=reg72*reg46; T reg173=reg98*reg44; reg91=reg91*reg37; T reg174=reg24*reg38; T reg175=reg52*reg49;
    reg80=reg81+reg80; T reg176=reg71*reg44; reg61=reg61*reg37; T reg177=reg62*reg46; T reg178=reg70*reg49;
    T reg179=reg68*reg44; T reg180=reg64*reg2; T reg181=reg66*reg49; T reg182=reg90*reg44; T reg183=reg77*reg39;
    T reg184=reg64*reg49; reg87=reg87*reg44; T reg185=reg69*reg40; reg64=reg64*reg45; reg77=reg77*reg40;
    T reg186=reg78*reg40; T reg187=reg98*reg41; T reg188=reg72*reg38; T reg189=reg66*reg45; T reg190=reg71*reg41;
    T reg191=reg52*reg38; T reg192=reg70*reg45; reg75=reg76+reg75; reg90=reg90*reg41; T reg193=reg66*reg38;
    reg68=reg68*reg37; T reg194=reg70*reg46; reg85=reg88+reg85; reg82=reg83+reg82; T reg195=reg72*reg49;
    reg98=reg98*reg37; reg97=reg99+reg97; T reg196=reg57*reg39; reg24=reg24*reg2; T reg197=reg56*reg43;
    T reg198=reg52*reg50; T reg199=reg73*reg39; reg65=reg65*reg2; T reg200=reg74*reg43; reg74=reg74*reg39;
    T reg201=reg72*reg2; reg72=reg72*reg50; T reg202=reg56*reg39; reg52=reg52*reg2; reg89=reg89+reg79;
    T reg203=reg63*reg40; T reg204=reg78*reg39; T reg205=reg66*reg2; reg73=reg73*reg43; reg78=reg78*reg43;
    reg66=reg66*reg50; T reg206=reg86-reg84; T reg207=reg69*reg43; T reg208=reg69*reg42; reg63=reg63*reg43;
    T reg209=reg70*reg50; T reg210=reg96-reg100; reg90=reg193-reg90; reg151=reg150+reg151; reg147=reg146+reg147;
    reg184=reg87+reg184; reg87=reg75*reg26; reg146=reg207+reg129; reg145=reg145-reg144; reg142=reg142-reg143;
    reg150=reg89*reg26; reg58=reg58+reg107; reg141=reg140+reg141; reg191=reg191+reg190; reg140=reg185+reg125;
    reg139=reg138+reg139; reg204=reg205-reg204; reg138=reg97*reg26; reg206=reg206*reg26; reg187=reg188-reg187;
    reg135=reg134-reg135; reg132=reg132-reg133; reg118=reg162-reg118; reg175=reg175-reg176; reg24=reg24+reg196;
    reg134=reg106+reg124; reg153=reg154-reg153; reg155=reg156-reg155; reg154=reg105+reg121; reg199=reg65-reg199;
    reg181=reg182+reg181; reg210=reg210*reg26; reg65=reg26*reg93; reg195=reg173+reg195; reg74=reg201-reg74;
    reg174=reg174+reg161; reg156=reg104+reg120; reg162=reg85*reg26; reg179=reg179-reg178; reg173=reg208+reg103;
    reg165=reg168+reg165; reg148=reg148-reg171; reg52=reg52+reg202; reg149=reg149-reg152; reg168=reg26*reg102;
    reg164=reg163+reg164; reg203=reg203+reg192; reg160=reg160-reg159; reg137=reg136-reg137; reg157=reg157-reg158;
    reg68=reg68+reg194; reg136=reg82*reg26; reg172=reg98-reg172; reg115=reg91-reg115; reg91=reg80*reg26;
    reg61=reg61+reg177; reg131=reg130+reg131; reg128=reg127+reg128; reg183=reg180-reg183; reg123=reg123-reg126;
    reg66=reg78+reg66; reg122=reg122-reg119; reg63=reg63-reg209; reg117=reg116+reg117; reg198=reg198-reg197;
    reg114=reg113-reg114; reg111=reg111-reg112; reg72=reg200+reg72; reg110=reg110-reg109; reg108=reg73+reg108;
    reg64=reg77-reg64; reg167=reg166+reg167; reg189=reg186-reg189; reg170=reg170-reg169; reg135=reg26*reg135;
    reg206=ponderation*reg206; reg123=reg123*reg26; reg73=reg140*reg26; reg165=reg26*reg165; reg184=reg184*reg26;
    reg128=reg128*reg26; reg77=reg26*reg156; reg78=reg26*reg154; reg61=reg61*reg26; reg181=reg181*reg26;
    reg98=reg26*reg134; reg131=reg131*reg26; reg167=reg26*reg167; reg58=reg26*reg58; reg179=reg179*reg26;
    reg113=ponderation*reg168; reg116=ponderation*reg91; reg118=reg26*reg118; reg115=reg115*reg26; reg153=reg26*reg153;
    reg175=reg175*reg26; reg183=reg183*reg26; reg108=reg108*reg26; reg204=reg204*reg26; reg72=reg72*reg26;
    reg110=reg110*reg26; reg127=ponderation*reg150; reg170=reg26*reg170; reg52=reg52*reg26; reg111=reg111*reg26;
    reg74=reg74*reg26; reg132=reg26*reg132; reg198=reg198*reg26; reg199=reg199*reg26; reg117=reg117*reg26;
    reg24=reg24*reg26; reg63=reg63*reg26; reg130=ponderation*reg138; reg187=reg187*reg26; reg210=ponderation*reg210;
    reg64=reg64*reg26; reg122=reg122*reg26; reg163=reg173*reg26; reg66=reg66*reg26; reg166=reg146*reg26;
    reg203=reg203*reg26; reg174=reg26*reg174; reg68=reg68*reg26; reg191=reg191*reg26; reg180=ponderation*reg162;
    reg114=reg114*reg26; reg148=reg26*reg148; reg157=reg26*reg157; reg149=reg26*reg149; reg137=reg137*reg26;
    reg189=reg189*reg26; reg151=reg26*reg151; reg141=reg26*reg141; reg90=reg90*reg26; reg147=reg26*reg147;
    reg160=reg26*reg160; reg145=reg26*reg145; reg182=ponderation*reg87; reg142=reg26*reg142; reg164=reg26*reg164;
    reg172=reg172*reg26; reg139=reg26*reg139; reg186=ponderation*reg136; reg195=reg195*reg26; reg188=ponderation*reg65;
    reg155=reg26*reg155; T tmp_5_3=ponderation*reg198; T tmp_3_4=ponderation*reg74; T tmp_7_5=ponderation*reg153; T tmp_4_4=ponderation*reg147;
    T tmp_1_7=ponderation*reg132; T tmp_6_1=ponderation*reg90; T tmp_3_5=ponderation*reg199; T tmp_5_2=ponderation*reg63; T tmp_6_4=ponderation*reg187;
    T tmp_2_3=-reg186; T tmp_3_6=ponderation*reg24; T tmp_0_5=ponderation*reg117; T tmp_2_5=ponderation*reg115; T tmp_4_5=ponderation*reg151;
    T tmp_7_3=-reg188; T tmp_7_4=ponderation*reg155; T tmp_1_4=ponderation*reg164; T tmp_3_0=ponderation*reg183; T tmp_0_3=ponderation*reg175;
    T tmp_5_5=ponderation*reg108; T tmp_5_4=ponderation*reg72; T tmp_4_2=ponderation*reg142; T tmp_3_1=ponderation*reg204; T tmp_6_2=-reg182;
    T tmp_1_6=ponderation*reg170; T tmp_7_1=ponderation*reg189; T tmp_3_2=-reg127; T tmp_5_6=ponderation*reg110; T tmp_2_4=ponderation*reg172;
    T tmp_4_3=ponderation*reg145; T tmp_4_1=ponderation*reg141; T tmp_3_3=ponderation*reg52; T tmp_1_3=ponderation*reg160; T tmp_5_7=ponderation*reg111;
    T tmp_7_0=ponderation*reg64; T tmp_0_2=ponderation*reg179; T tmp_0_0=ponderation*reg184; T tmp_1_5=ponderation*reg167; T tmp_4_0=ponderation*reg139;
    reg24=ponderation*reg77; sollicitation[indices[1]+0]+=reg24; T tmp_2_6=-reg116; T tmp_1_0=ponderation*reg128; T tmp_2_7=ponderation*reg61;
    reg52=ponderation*reg78; sollicitation[indices[0]+1]+=reg52; T tmp_2_2=ponderation*reg68; T tmp_6_6=ponderation*reg174; reg61=ponderation*reg98;
    sollicitation[indices[0]+0]+=reg61; T tmp_7_7=ponderation*reg58; T tmp_3_7=-reg130; T tmp_0_1=ponderation*reg181; T tmp_7_2=ponderation*reg203;
    T tmp_1_1=ponderation*reg131; T tmp_0_4=ponderation*reg195; T tmp_5_1=ponderation*reg66; T tmp_6_5=ponderation*reg118; sollicitation[indices[1]+1]+=-reg210;
    T tmp_2_1=ponderation*reg137; T tmp_6_3=ponderation*reg191; T tmp_4_6=ponderation*reg149; reg58=ponderation*reg163; sollicitation[indices[2]+0]+=reg58;
    T tmp_0_6=ponderation*reg122; T tmp_1_2=ponderation*reg157; T tmp_6_0=ponderation*reg114; reg63=ponderation*reg166; sollicitation[indices[2]+1]+=reg63;
    T tmp_5_0=ponderation*reg165; sollicitation[indices[3]+0]+=-reg206; T tmp_2_0=ponderation*reg135; T tmp_4_7=ponderation*reg148; T tmp_7_6=-reg113;
    reg64=ponderation*reg73; sollicitation[indices[3]+1]+=reg64; T tmp_0_7=ponderation*reg123; T tmp_6_7=-reg180;
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
      const Number<true> &assemble_vec,
      const Element<Quad,DefaultBehavior,Node<2,T_pos,ND>,ED,nim> &elem,
      const unsigned *indices){ 
  #define PNODE(N) (*elem.node(N))
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg4*reg2; reg2=reg3*reg2; T reg6=reg3*reg1; reg1=reg4*reg1;
    T reg7=reg4*reg5; T reg8=reg3*reg2; reg5=reg3*reg5; T reg9=reg4*reg1; T reg10=reg3*reg6;
    reg1=reg3*reg1; reg9=reg9-reg10; reg1=reg10+reg1; reg6=reg4*reg6; reg7=reg7-reg8;
    reg5=reg8+reg5; reg2=reg4*reg2; T reg11=reg10+reg6; T reg12=reg4*reg7; reg2=reg8+reg2;
    reg8=reg3*reg5; reg1=reg3*reg1; reg9=reg4*reg9; reg11=reg3*reg11; reg1=reg9-reg1;
    reg9=reg3*reg2; reg8=reg12-reg8; reg12=1-var_inter[0]; reg11=reg1-reg11; reg9=reg8-reg9;
    reg1=1-var_inter[1]; reg8=reg1*elem.pos(0)[0]; T reg13=reg1*elem.pos(0)[1]; T reg14=reg1*elem.pos(1)[0]; T reg15=reg1*elem.pos(1)[1];
    reg5=reg5/reg9; reg7=reg7/reg9; T reg16=reg12*elem.pos(0)[1]; T reg17=elem.pos(1)[1]*var_inter[0]; T reg18=elem.pos(1)[0]*var_inter[0];
    T reg19=reg12*elem.pos(0)[0]; reg11=reg11/reg9; T reg20=reg17+reg16; T reg21=elem.pos(2)[0]*var_inter[1]; T reg22=elem.pos(2)[1]*var_inter[0];
    T reg23=reg7*reg11; T reg24=reg5*reg11; T reg25=elem.pos(2)[0]*var_inter[0]; T reg26=reg19+reg18; reg15=reg15-reg13;
    reg14=reg14-reg8; T reg27=elem.pos(2)[1]*var_inter[1]; reg9=reg2/reg9; reg2=elem.pos(3)[1]*reg12; reg22=reg22-reg20;
    T reg28=elem.pos(3)[0]*reg12; reg25=reg25-reg26; reg15=reg27+reg15; reg27=elem.pos(3)[1]*var_inter[1]; T reg29=reg7*(*f.m).alpha;
    T reg30=reg5*(*f.m).alpha; T reg31=elem.pos(3)[0]*var_inter[1]; T reg32=reg5*reg24; reg21=reg14+reg21; reg14=reg7*reg23;
    reg2=reg22+reg2; reg9=reg9*(*f.m).alpha; reg30=reg29+reg30; reg28=reg25+reg28; reg15=reg15-reg27;
    reg32=reg14-reg32; reg21=reg21-reg31; reg23=reg23/reg32; reg32=reg24/reg32; reg14=reg21*reg2;
    reg22=reg15*reg28; reg9=reg30+reg9; reg24=reg3*reg0; reg0=reg4*reg0; reg32=reg32*reg9;
    reg25=pow(reg3,2); reg22=reg14-reg22; reg3=reg3*reg24; reg14=reg4*reg0; reg9=reg23*reg9;
    reg4=pow(reg4,2); reg32=reg9-reg32; reg3=reg14-reg3; reg25=reg4-reg25; reg4=1-(*f.m).resolution;
    reg15=reg15/reg22; reg2=reg2/reg22; reg21=reg21/reg22; reg28=reg28/reg22; reg9=reg1*reg2;
    reg14=reg12*reg15; reg23=reg12*reg21; reg29=var_inter[0]*reg15; reg30=(*f.m).resolution*(*f.m).alpha; reg32=reg4*reg32;
    reg24=reg24/reg3; reg0=reg0/reg3; T reg33=var_inter[1]*reg2; reg3=reg25/reg3; reg25=var_inter[1]*reg28;
    T reg34=reg9+reg29; T reg35=reg1*reg28; T reg36=var_inter[0]*reg21; reg5=reg5*reg4; T reg37=(*f.m).resolution*reg24;
    reg32=reg30+reg32; reg30=(*f.m).resolution*reg3; reg11=reg11*reg4; T reg38=reg14+reg33; reg4=reg7*reg4;
    reg7=reg25+reg23; T reg39=(*f.m).resolution*reg0; T reg40=reg36-reg25; T reg41=reg33-reg29; reg32=(*f.m).deltaT*reg32;
    reg11=reg30+reg11; reg5=reg37+reg5; reg30=0.5*reg7; reg37=reg35+reg36; reg39=reg4+reg39;
    reg4=reg35-reg23; T reg42=reg14-reg9; T reg43=0.5*reg38; T reg44=0.5*reg34; T reg45=0.5*reg37;
    T reg46=reg39*reg32; T reg47=reg5*reg32; T reg48=0.5*reg40; T reg49=0.5*reg41; T reg50=reg11*reg44;
    T reg51=reg11*reg43; T reg52=0.5*reg42; T reg53=0.5*reg4; T reg54=reg11*reg30; T reg55=reg11*reg49;
    T reg56=reg39*reg38; reg54=2*reg54; T reg57=reg5*reg7; T reg58=2*reg51; T reg59=reg39*reg7;
    T reg60=reg47+reg46; T reg61=reg11*reg52; T reg62=reg11*reg45; T reg63=reg11*reg53; T reg64=reg5*reg37;
    reg50=2*reg50; T reg65=reg11*reg48; T reg66=reg12*var_inter[1]; T reg67=reg1*var_inter[0]; reg55=2*reg55;
    T reg68=reg39*reg41; T reg69=reg39*reg4; T reg70=reg5*reg34; reg65=2*reg65; T reg71=2*reg62;
    T reg72=reg54*reg45; T reg73=reg39*reg34; reg61=2*reg61; T reg74=reg5*reg4; T reg75=reg39*reg37;
    T reg76=reg5*reg41; T reg77=reg39*reg40; T reg78=reg5*reg38; reg63=2*reg63; T reg79=reg60*reg37;
    T reg80=reg39*reg42; T reg81=reg60*reg38; T reg82=reg59*reg37; T reg83=reg58*reg44; T reg84=reg67*elem.f_vol_e[1];
    T reg85=reg66*elem.f_vol_e[0]; T reg86=reg1*reg12; T reg87=var_inter[0]*var_inter[1]; T reg88=reg5*reg40; T reg89=reg58*reg30;
    T reg90=reg57*reg38; T reg91=reg64*reg34; T reg92=reg50*reg45; T reg93=reg56*reg34; T reg94=reg58*reg43;
    T reg95=reg59*reg7; T reg96=reg56*reg38; T reg97=reg63*reg53; T reg98=reg58*reg48; T reg99=reg57*reg41;
    T reg100=reg54*reg48; T reg101=reg56*reg41; T reg102=reg55*reg48; T reg103=reg88*reg41; T reg104=reg65*reg48;
    T reg105=reg68*reg41; T reg106=reg58*reg52; T reg107=reg59*reg4; T reg108=reg54*reg52; T reg109=reg78*reg4;
    T reg110=reg55*reg52; T reg111=reg77*reg4; T reg112=reg65*reg52; T reg113=reg76*reg4; T reg114=reg50*reg52;
    T reg115=reg75*reg4; T reg116=reg71*reg52; T reg117=reg70*reg4; reg92=reg91+reg92; T reg118=reg74*reg42;
    T reg119=reg61*reg53; T reg120=reg73*reg42; T reg121=reg68*reg34; T reg122=reg67*elem.f_vol_e[0]; reg72=reg93+reg72;
    T reg123=reg86*elem.f_vol_e[1]; T reg124=reg50*reg44; T reg125=reg86*elem.f_vol_e[0]; T reg126=reg75*reg37; T reg127=reg65*reg44;
    T reg128=reg66*elem.f_vol_e[1]; T reg129=reg76*reg37; reg90=reg89+reg90; T reg130=reg87*elem.f_vol_e[1]; T reg131=reg55*reg44;
    T reg132=reg87*elem.f_vol_e[0]; T reg133=reg77*reg37; T reg134=reg54*reg44; T reg135=reg78*reg37; reg82=reg83+reg82;
    T reg136=reg79-reg84; T reg137=reg60*reg41; T reg138=reg60*reg40; T reg139=reg81-reg85; T reg140=reg60*reg7;
    T reg141=reg65*reg53; T reg142=reg73*reg34; T reg143=reg60*reg34; T reg144=reg60*reg4; T reg145=reg60*reg42;
    T reg146=reg80*reg42; T reg147=reg71*reg45; T reg148=reg61*reg52; T reg149=reg69*reg4; T reg150=reg88*reg34;
    T reg151=reg58*reg53; T reg152=reg57*reg42; T reg153=reg54*reg53; T reg154=reg54*reg30; T reg155=reg68*reg42;
    T reg156=reg56*reg42; T reg157=reg55*reg53; T reg158=reg88*reg42; T reg159=reg55*reg45; T reg160=reg58*reg49;
    reg59=reg59*reg40; T reg161=reg54*reg49; T reg162=reg78*reg40; T reg163=reg55*reg49; T reg164=reg77*reg40;
    T reg165=reg58*reg45; reg57=reg57*reg34; T reg166=reg65*reg45; T reg167=reg71*reg53; T reg168=reg64*reg42;
    T reg169=reg50*reg53; reg117=reg117-reg116; reg104=reg105+reg104; reg134=reg134+reg135; reg169=reg169-reg168;
    reg102=reg103+reg102; reg59=reg59-reg160; reg103=reg82*reg22; reg100=reg100-reg101; reg136=reg136*reg22;
    reg105=reg137+reg132; reg99=reg99-reg98; T reg170=reg138+reg130; reg153=reg153-reg156; reg114=reg114-reg115;
    reg154=reg154+reg96; reg141=reg155+reg141; reg155=reg92*reg22; reg139=reg139*reg22; reg95=reg95+reg94;
    T reg171=reg140+reg128; reg142=reg142+reg147; reg166=reg121-reg166; reg159=reg150-reg159; reg121=reg145+reg125;
    reg150=reg143+reg122; T reg172=reg144+reg123; reg157=reg158+reg157; reg110=reg111+reg110; reg163=reg164+reg163;
    reg148=reg149+reg148; reg129=reg127-reg129; reg124=reg124+reg126; reg120=reg120-reg167; reg108=reg108-reg109;
    reg161=reg161-reg162; reg112=reg113+reg112; reg97=reg146+reg97; reg111=reg72*reg22; reg133=reg131-reg133;
    reg107=reg107-reg106; reg119=reg118+reg119; reg57=reg57+reg165; reg113=reg90*reg22; reg152=reg152-reg151;
    reg159=reg159*reg22; reg118=reg171*reg22; reg127=reg22*reg150; reg57=reg57*reg22; reg139=ponderation*reg139;
    reg124=reg124*reg22; reg141=reg141*reg22; reg131=reg170*reg22; reg163=reg163*reg22; reg129=reg129*reg22;
    reg119=reg119*reg22; reg146=ponderation*reg155; reg149=reg105*reg22; reg136=ponderation*reg136; reg161=reg161*reg22;
    reg59=reg59*reg22; reg158=ponderation*reg103; reg133=reg133*reg22; reg164=ponderation*reg113; reg134=reg134*reg22;
    reg120=reg120*reg22; reg152=reg152*reg22; reg107=reg22*reg107; reg104=reg22*reg104; reg108=reg22*reg108;
    reg117=reg22*reg117; reg169=reg169*reg22; reg102=reg22*reg102; reg166=reg166*reg22; reg142=reg142*reg22;
    reg100=reg22*reg100; reg110=reg22*reg110; reg153=reg153*reg22; reg99=reg22*reg99; reg154=reg22*reg154;
    reg114=reg22*reg114; T reg173=reg22*reg172; reg157=reg157*reg22; T reg174=reg22*reg121; T reg175=ponderation*reg111;
    reg148=reg148*reg22; reg112=reg22*reg112; reg97=reg97*reg22; reg95=reg22*reg95; T tmp_0_7=ponderation*reg152;
    T tmp_2_3=-reg146; T tmp_5_6=ponderation*reg161; T tmp_1_6=ponderation*reg108; T tmp_2_2=ponderation*reg142; T tmp_2_4=ponderation*reg166;
    T tmp_0_4=ponderation*reg141; T tmp_1_3=ponderation*reg114; T tmp_3_4=ponderation*reg129; T tmp_1_4=ponderation*reg112; T tmp_1_2=ponderation*reg117;
    T tmp_2_6=-reg175; T tmp_1_1=ponderation*reg148; T tmp_0_1=ponderation*reg119; T tmp_5_5=ponderation*reg163; T tmp_1_5=ponderation*reg110;
    T tmp_3_3=ponderation*reg124; T tmp_2_7=ponderation*reg57; reg57=ponderation*reg173; sollicitation[indices[0]+1]+=reg57; reg108=ponderation*reg174;
    sollicitation[indices[0]+0]+=reg108; reg110=ponderation*reg127; sollicitation[indices[1]+0]+=reg110; T tmp_0_2=ponderation*reg120; reg112=ponderation*reg118;
    sollicitation[indices[3]+1]+=reg112; T tmp_0_5=ponderation*reg157; T tmp_7_7=ponderation*reg95; T tmp_2_5=ponderation*reg159; sollicitation[indices[3]+0]+=-reg139;
    T tmp_0_0=ponderation*reg97; T tmp_6_6=ponderation*reg154; reg95=ponderation*reg131; sollicitation[indices[2]+1]+=reg95; T tmp_4_7=ponderation*reg99;
    reg97=ponderation*reg149; sollicitation[indices[2]+0]+=reg97; T tmp_0_6=ponderation*reg153; T tmp_5_7=ponderation*reg59; sollicitation[indices[1]+1]+=-reg136;
    T tmp_4_6=ponderation*reg100; T tmp_0_3=ponderation*reg169; T tmp_4_5=ponderation*reg102; T tmp_6_7=-reg164; T tmp_3_6=ponderation*reg134;
    T tmp_4_4=ponderation*reg104; T tmp_3_7=-reg158; T tmp_1_7=ponderation*reg107; T tmp_3_5=ponderation*reg133;
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
    T reg4=reg0*reg1; T reg5=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg6=1.0/(*f.m).elastic_modulus; T reg7=reg3*elem.pos(0)[0]; T reg8=reg3*elem.pos(1)[0];
    T reg9=reg3*elem.pos(1)[1]; T reg10=reg3*elem.pos(0)[1]; T reg11=reg2*elem.pos(0)[0]; T reg12=elem.pos(1)[0]*var_inter[0]; T reg13=reg2*elem.pos(0)[1];
    T reg14=elem.pos(1)[1]*var_inter[0]; T reg15=reg11+reg12; T reg16=elem.pos(2)[0]*var_inter[0]; T reg17=reg6*reg4; reg4=reg5*reg4;
    T reg18=elem.pos(2)[1]*var_inter[0]; T reg19=reg5*reg1; reg1=reg6*reg1; reg8=reg8-reg7; T reg20=elem.pos(2)[0]*var_inter[1];
    T reg21=elem.pos(2)[1]*var_inter[1]; T reg22=reg14+reg13; reg9=reg9-reg10; reg18=reg18-reg22; T reg23=elem.pos(3)[0]*var_inter[1];
    reg20=reg8+reg20; reg8=elem.pos(3)[1]*reg2; reg9=reg21+reg9; reg21=reg6*reg17; T reg24=elem.pos(3)[1]*var_inter[1];
    T reg25=reg5*reg4; T reg26=reg5*reg1; reg17=reg5*reg17; T reg27=elem.pos(3)[0]*reg2; reg16=reg16-reg15;
    T reg28=reg5*reg19; reg1=reg6*reg1; reg21=reg21-reg25; reg4=reg6*reg4; reg8=reg18+reg8;
    reg17=reg25+reg17; reg20=reg20-reg23; reg26=reg28+reg26; reg1=reg1-reg28; reg19=reg6*reg19;
    reg27=reg16+reg27; reg9=reg9-reg24; reg16=reg6*reg0; reg0=reg5*reg0; reg18=reg28+reg19;
    reg26=reg5*reg26; T reg29=reg9*reg27; T reg30=reg20*reg8; reg4=reg25+reg4; reg1=reg6*reg1;
    reg25=reg5*reg17; T reg31=reg6*reg21; T reg32=pow(reg5,2); T reg33=pow(reg6,2); T reg34=reg5*reg0;
    reg26=reg1-reg26; reg29=reg30-reg29; reg1=reg5*reg4; reg25=reg31-reg25; reg18=reg5*reg18;
    reg6=reg6*reg16; reg20=reg20/reg29; reg8=reg8/reg29; reg27=reg27/reg29; reg9=reg9/reg29;
    reg18=reg26-reg18; reg32=reg33-reg32; reg1=reg25-reg1; reg34=reg6-reg34; reg32=reg32/reg34;
    reg5=var_inter[1]*reg8; reg18=reg18/reg1; reg6=reg3*reg8; reg25=reg2*reg20; reg26=reg3*reg27;
    reg30=var_inter[0]*reg20; reg31=var_inter[1]*reg27; reg33=var_inter[0]*reg9; T reg35=1-(*f.m).resolution; T reg36=reg2*reg9;
    T reg37=(*f.m).resolution*reg32; T reg38=reg18*reg35; T reg39=reg26+reg30; reg16=reg16/reg34; reg34=reg0/reg34;
    reg0=reg6+reg33; T reg40=reg36+reg5; T reg41=reg31+reg25; reg17=reg17/reg1; reg21=reg21/reg1;
    T reg42=reg26-reg25; T reg43=0.5*reg39; T reg44=0.5*reg0; T reg45=reg36-reg6; T reg46=reg17*reg35;
    T reg47=(*f.m).resolution*reg34; T reg48=reg5-reg33; T reg49=reg30-reg31; T reg50=reg21*reg35; T reg51=(*f.m).resolution*reg16;
    reg38=reg37+reg38; reg37=0.5*reg41; T reg52=0.5*reg40; T reg53=reg38*reg37; T reg54=0.5*reg45;
    T reg55=0.5*reg42; reg46=reg47+reg46; reg47=reg38*reg52; T reg56=reg38*reg43; T reg57=0.5*reg48;
    T reg58=reg38*reg44; reg51=reg50+reg51; reg50=0.5*reg49; T reg59=reg46*reg0; reg53=2*reg53;
    T reg60=reg51*reg40; T reg61=reg38*reg54; T reg62=reg51*reg41; T reg63=reg51*reg39; T reg64=reg46*reg41;
    T reg65=reg38*reg57; T reg66=reg46*reg40; T reg67=reg51*reg0; reg58=2*reg58; T reg68=2*reg56;
    T reg69=reg38*reg55; T reg70=2*reg47; T reg71=reg38*reg50; T reg72=reg46*reg39; T reg73=reg68*reg44;
    T reg74=reg67*reg40; T reg75=reg51*reg42; T reg76=reg68*reg37; T reg77=reg51*reg49; reg71=2*reg71;
    T reg78=reg46*reg48; T reg79=reg66*reg41; T reg80=reg53*reg52; T reg81=reg58*reg43; reg69=2*reg69;
    T reg82=reg72*reg0; T reg83=reg51*reg45; T reg84=reg62*reg39; T reg85=reg70*reg44; T reg86=reg46*reg49;
    T reg87=reg63*reg41; T reg88=reg58*reg52; T reg89=reg59*reg39; reg65=2*reg65; reg61=2*reg61;
    T reg90=reg46*reg42; T reg91=reg70*reg37; T reg92=reg64*reg40; T reg93=reg53*reg43; T reg94=reg60*reg0;
    T reg95=reg51*reg48; T reg96=reg46*reg45; T reg97=reg58*reg57; T reg98=reg78*reg49; T reg99=reg62*reg49;
    T reg100=reg53*reg57; T reg101=reg70*reg57; T reg102=reg69*reg37; T reg103=reg83*reg40; T reg104=reg65*reg43;
    T reg105=reg86*reg45; T reg106=reg66*reg49; T reg107=reg65*reg57; T reg108=reg71*reg57; T reg109=reg86*reg40;
    T reg110=reg65*reg55; T reg111=reg77*reg49; T reg112=reg60*reg45; T reg113=reg53*reg55; T reg114=reg61*reg54;
    T reg115=reg64*reg45; T reg116=reg75*reg42; T reg117=reg69*reg54; T reg118=reg96*reg42; T reg119=reg70*reg55;
    T reg120=reg61*reg50; T reg121=reg67*reg48; T reg122=reg68*reg50; T reg123=reg72*reg48; T reg124=reg65*reg52;
    T reg125=reg77*reg41; T reg126=reg71*reg52; T reg127=reg78*reg41; reg88=reg87+reg88; T reg128=reg60*reg40;
    T reg129=reg65*reg37; T reg130=reg69*reg57; T reg131=reg96*reg49; T reg132=reg70*reg50; T reg133=reg64*reg48;
    T reg134=reg53*reg50; T reg135=reg60*reg48; T reg136=reg65*reg50; T reg137=reg86*reg48; T reg138=reg71*reg50;
    T reg139=reg95*reg48; T reg140=reg58*reg50; T reg141=reg70*reg52; T reg142=reg62*reg41; reg80=reg79+reg80;
    T reg143=reg59*reg42; T reg144=reg68*reg54; T reg145=reg63*reg42; T reg146=reg58*reg54; T reg147=reg78*reg42;
    T reg148=reg71*reg54; T reg149=reg77*reg42; T reg150=reg65*reg54; T reg151=reg66*reg42; T reg152=reg53*reg54;
    reg62=reg62*reg42; T reg153=reg70*reg54; T reg154=reg83*reg0; T reg155=reg69*reg43; T reg156=reg90*reg0;
    T reg157=reg61*reg43; T reg158=reg83*reg48; T reg159=reg69*reg50; T reg160=reg90*reg48; T reg161=reg90*reg45;
    T reg162=reg61*reg55; T reg163=reg70*reg43; reg64=reg64*reg0; T reg164=reg67*reg45; T reg165=reg68*reg55;
    reg93=reg94+reg93; T reg166=reg72*reg45; T reg167=reg58*reg55; T reg168=reg53*reg37; T reg169=reg95*reg45;
    T reg170=reg71*reg55; reg86=reg86*reg0; reg92=reg91+reg92; T reg171=reg71*reg43; T reg172=reg95*reg0;
    T reg173=reg61*reg37; reg81=reg82+reg81; reg90=reg90*reg40; T reg174=reg68*reg43; reg67=reg67*reg0;
    reg74=reg76+reg74; T reg175=reg68*reg52; T reg176=reg58*reg37; T reg177=reg72*reg40; T reg178=reg61*reg52;
    T reg179=reg75*reg41; T reg180=reg71*reg37; T reg181=reg69*reg52; reg95=reg95*reg40; T reg182=reg96*reg41;
    T reg183=reg63*reg49; T reg184=reg61*reg44; T reg185=reg75*reg39; T reg186=reg59*reg41; reg89=reg89+reg73;
    reg58=reg58*reg44; T reg187=reg63*reg39; T reg188=reg68*reg57; reg71=reg71*reg44; reg78=reg78*reg39;
    reg59=reg59*reg49; reg65=reg65*reg44; reg77=reg77*reg39; reg53=reg53*reg44; T reg189=reg66*reg39;
    reg84=reg85+reg84; T reg190=reg69*reg44; reg96=reg96*reg39; reg69=reg69*reg55; reg83=reg83*reg45;
    reg75=reg75*reg49; reg61=reg61*reg57; reg58=reg58+reg187; T reg191=reg29*reg80; reg62=reg62-reg153;
    reg133=reg133-reg132; T reg192=reg74*reg29; reg152=reg152-reg151; reg126=reg127-reg126; reg127=reg89*reg29;
    T reg193=reg84*reg29; reg150=reg149+reg150; reg124=reg125-reg124; reg185=reg184-reg185; reg176=reg176+reg177;
    reg148=reg147+reg148; reg69=reg83+reg69; reg109=reg129-reg109; reg95=reg180-reg95; reg146=reg146-reg145;
    reg130=reg131+reg130; reg134=reg134-reg135; reg53=reg53+reg189; reg136=reg137+reg136; reg164=reg164-reg165;
    reg167=reg167-reg166; reg77=reg65-reg77; reg138=reg139+reg138; reg170=reg169+reg170; reg140=reg140-reg123;
    reg168=reg168+reg128; reg121=reg121-reg122; reg65=reg92*reg29; reg142=reg142+reg141; reg78=reg71-reg78;
    reg120=reg160+reg120; reg159=reg158+reg159; reg71=reg29*reg88; reg162=reg161+reg162; reg90=reg173-reg90;
    reg155=reg154-reg155; reg108=reg98+reg108; reg157=reg156-reg157; reg67=reg67+reg174; reg59=reg59-reg188;
    reg107=reg111+reg107; reg83=reg81*reg29; reg171=reg172-reg171; reg100=reg100-reg106; reg114=reg116+reg114;
    reg117=reg118+reg117; reg61=reg75+reg61; reg104=reg86-reg104; reg99=reg99-reg101; reg115=reg115-reg119;
    reg75=reg93*reg29; reg103=reg102-reg103; reg113=reg113-reg112; reg64=reg64+reg163; reg96=reg190-reg96;
    reg110=reg105+reg110; reg97=reg97-reg183; reg143=reg143-reg144; reg181=reg182-reg181; reg178=reg179-reg178;
    reg186=reg186+reg175; reg86=ponderation*reg75; reg143=reg29*reg143; reg164=reg164*reg29; reg99=reg99*reg29;
    reg152=reg29*reg152; reg133=reg29*reg133; reg115=reg115*reg29; reg53=reg53*reg29; reg104=reg104*reg29;
    reg59=reg59*reg29; reg134=reg29*reg134; reg181=reg181*reg29; reg167=reg167*reg29; reg148=reg29*reg148;
    reg61=reg61*reg29; reg136=reg29*reg136; reg97=reg97*reg29; reg117=reg117*reg29; reg110=reg110*reg29;
    reg69=reg69*reg29; reg109=reg29*reg109; reg124=reg29*reg124; reg95=reg95*reg29; reg130=reg29*reg130;
    reg98=ponderation*reg191; reg126=reg29*reg126; reg64=reg64*reg29; reg96=reg96*reg29; reg146=reg29*reg146;
    reg102=ponderation*reg71; reg113=reg113*reg29; reg162=reg162*reg29; reg105=ponderation*reg193; reg186=reg186*reg29;
    reg142=reg29*reg142; reg168=reg29*reg168; reg120=reg29*reg120; reg107=reg107*reg29; reg103=reg103*reg29;
    reg159=reg29*reg159; reg111=ponderation*reg83; reg58=reg58*reg29; reg116=ponderation*reg192; reg150=reg29*reg150;
    reg185=reg185*reg29; reg67=reg67*reg29; reg155=reg29*reg155; reg178=reg178*reg29; reg90=reg90*reg29;
    reg118=ponderation*reg127; reg62=reg29*reg62; reg108=reg108*reg29; reg157=reg157*reg29; reg140=reg29*reg140;
    reg176=reg176*reg29; reg114=reg114*reg29; reg100=reg100*reg29; reg170=reg170*reg29; reg121=reg29*reg121;
    reg171=reg171*reg29; reg78=reg78*reg29; reg138=reg29*reg138; reg125=ponderation*reg65; reg77=reg77*reg29;
    T tmp_3_2=-reg118; T tmp_6_5=ponderation*reg109; T tmp_5_6=ponderation*reg100; T tmp_3_5=ponderation*reg77; T tmp_5_4=ponderation*reg108;
    T tmp_3_1=ponderation*reg185; T tmp_5_3=ponderation*reg97; T tmp_7_6=-reg98; T tmp_3_4=ponderation*reg78; T tmp_5_0=ponderation*reg130;
    T tmp_3_0=ponderation*reg96; T tmp_3_3=ponderation*reg58; T tmp_5_2=ponderation*reg59; T tmp_7_7=ponderation*reg142; T tmp_3_6=ponderation*reg53;
    T tmp_5_5=ponderation*reg107; T tmp_5_7=ponderation*reg99; T tmp_5_1=ponderation*reg61; T tmp_4_2=ponderation*reg121; T tmp_1_1=ponderation*reg114;
    T tmp_4_1=ponderation*reg120; T tmp_2_3=-reg111; T tmp_6_7=-reg125; T tmp_4_0=ponderation*reg159; T tmp_6_0=ponderation*reg103;
    T tmp_3_7=-reg105; T tmp_2_2=ponderation*reg67; T tmp_2_0=ponderation*reg155; T tmp_6_1=ponderation*reg90; T tmp_1_7=ponderation*reg62;
    T tmp_2_1=ponderation*reg157; T tmp_1_6=ponderation*reg152; T tmp_6_2=-reg116; T tmp_7_1=ponderation*reg178; T tmp_1_5=ponderation*reg150;
    T tmp_7_0=ponderation*reg181; T tmp_6_3=ponderation*reg176; T tmp_1_4=ponderation*reg148; T tmp_1_3=ponderation*reg146; T tmp_1_2=ponderation*reg143;
    T tmp_6_4=ponderation*reg95; T tmp_7_5=ponderation*reg124; T tmp_0_5=ponderation*reg110; T tmp_0_0=ponderation*reg69; T tmp_2_7=ponderation*reg64;
    T tmp_7_4=ponderation*reg126; T tmp_7_3=-reg102; T tmp_0_1=ponderation*reg162; T tmp_0_6=ponderation*reg113; T tmp_7_2=ponderation*reg186;
    T tmp_6_6=ponderation*reg168; T tmp_2_6=-reg86; T tmp_2_5=ponderation*reg104; T tmp_0_2=ponderation*reg164; T tmp_4_7=ponderation*reg133;
    T tmp_0_7=ponderation*reg115; T tmp_4_6=ponderation*reg134; T tmp_0_3=ponderation*reg167; T tmp_4_5=ponderation*reg136; T tmp_1_0=ponderation*reg117;
    T tmp_4_4=ponderation*reg138; T tmp_2_4=ponderation*reg171; T tmp_4_3=ponderation*reg140; T tmp_0_4=ponderation*reg170;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1-var_inter[1]; T reg3=1-var_inter[0];
    T reg4=reg3*elem.pos(0)[1]; T reg5=elem.pos(1)[1]*var_inter[0]; T reg6=reg0*reg1; T reg7=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg8=1.0/(*f.m).elastic_modulus;
    T reg9=elem.pos(1)[0]*var_inter[0]; T reg10=reg3*elem.pos(0)[0]; T reg11=reg2*elem.pos(0)[1]; T reg12=reg2*elem.pos(1)[1]; T reg13=reg2*elem.pos(1)[0];
    T reg14=reg2*elem.pos(0)[0]; T reg15=reg5+reg4; T reg16=elem.pos(2)[1]*var_inter[0]; T reg17=reg7*reg6; reg6=reg8*reg6;
    T reg18=reg8*reg1; T reg19=elem.pos(2)[0]*var_inter[0]; T reg20=reg10+reg9; reg1=reg7*reg1; T reg21=elem.pos(2)[0]*var_inter[1];
    T reg22=elem.pos(2)[1]*var_inter[1]; reg13=reg13-reg14; reg12=reg12-reg11; T reg23=elem.pos(3)[0]*reg3; T reg24=reg8*reg18;
    T reg25=reg7*reg1; reg19=reg19-reg20; reg18=reg7*reg18; reg12=reg22+reg12; reg22=elem.pos(3)[1]*var_inter[1];
    T reg26=reg8*reg6; T reg27=elem.pos(3)[1]*reg3; T reg28=elem.pos(3)[0]*var_inter[1]; T reg29=reg7*reg17; reg6=reg7*reg6;
    reg16=reg16-reg15; reg21=reg13+reg21; reg12=reg12-reg22; reg27=reg16+reg27; reg24=reg24-reg25;
    reg18=reg25+reg18; reg21=reg21-reg28; reg23=reg19+reg23; reg1=reg8*reg1; reg17=reg8*reg17;
    reg26=reg26-reg29; reg6=reg29+reg6; reg13=reg7*reg6; reg16=reg8*reg26; reg24=reg8*reg24;
    reg18=reg7*reg18; reg19=reg21*reg27; T reg30=reg12*reg23; T reg31=reg25+reg1; T reg32=reg7*reg0;
    reg0=reg8*reg0; reg17=reg29+reg17; reg29=reg8*reg0; reg31=reg7*reg31; reg18=reg24-reg18;
    reg24=reg7*reg32; T reg33=pow(reg7,2); reg8=pow(reg8,2); reg30=reg19-reg30; reg13=reg16-reg13;
    reg7=reg7*reg17; reg7=reg13-reg7; reg12=reg12/reg30; reg23=reg23/reg30; reg31=reg18-reg31;
    reg21=reg21/reg30; reg27=reg27/reg30; reg24=reg29-reg24; reg33=reg8-reg33; reg33=reg33/reg24;
    reg8=1-(*f.m).resolution; reg31=reg31/reg7; reg13=var_inter[0]*reg12; reg16=reg3*reg12; reg18=reg2*reg27;
    reg19=reg3*reg21; reg29=var_inter[1]*reg27; T reg34=var_inter[1]*reg23; reg0=reg0/reg24; T reg35=(*f.m).resolution*reg33;
    T reg36=reg31*reg8; reg26=reg26/reg7; T reg37=reg16+reg29; reg6=reg6/reg7; T reg38=var_inter[0]*reg21;
    T reg39=reg18+reg13; T reg40=reg34+reg19; T reg41=reg2*reg23; reg24=reg32/reg24; reg32=reg16-reg18;
    T reg42=reg29-reg13; T reg43=0.5*reg37; reg36=reg35+reg36; reg35=(*f.m).resolution*reg24; T reg44=reg6*reg8;
    T reg45=(*f.m).resolution*reg0; T reg46=reg38-reg34; T reg47=0.5*reg39; T reg48=0.5*reg40; T reg49=reg41-reg19;
    T reg50=reg41+reg38; T reg51=reg26*reg8; T reg52=0.5*reg46; T reg53=0.5*reg42; T reg54=0.5*reg50;
    T reg55=reg36*reg48; reg44=reg35+reg44; reg45=reg51+reg45; reg35=reg36*reg43; reg51=0.5*reg49;
    T reg56=0.5*reg32; T reg57=reg36*reg47; T reg58=reg36*reg54; T reg59=reg36*reg56; T reg60=reg45*reg40;
    T reg61=reg44*reg50; T reg62=reg36*reg51; T reg63=reg36*reg53; T reg64=reg45*reg37; reg55=2*reg55;
    T reg65=reg36*reg52; T reg66=reg44*reg40; T reg67=2*reg35; reg57=2*reg57; T reg68=reg55*reg54;
    T reg69=reg64*reg39; T reg70=reg67*reg47; T reg71=reg60*reg50; T reg72=reg44*reg37; T reg73=reg44*reg39;
    T reg74=2*reg58; T reg75=reg45*reg49; T reg76=reg45*reg32; T reg77=reg45*reg39; reg63=2*reg63;
    T reg78=reg45*reg42; T reg79=reg44*reg42; T reg80=reg67*reg48; T reg81=reg66*reg37; reg62=2*reg62;
    reg65=2*reg65; T reg82=reg45*reg50; T reg83=reg45*reg46; T reg84=reg44*reg49; reg59=2*reg59;
    T reg85=reg61*reg39; T reg86=reg57*reg54; T reg87=reg44*reg46; T reg88=reg64*reg32; T reg89=reg82*reg49;
    T reg90=reg55*reg51; T reg91=reg87*reg42; T reg92=reg63*reg52; T reg93=reg64*reg42; T reg94=reg65*reg52;
    T reg95=reg57*reg56; T reg96=reg66*reg32; T reg97=reg79*reg49; T reg98=reg67*reg51; T reg99=reg83*reg49;
    T reg100=reg63*reg56; T reg101=reg59*reg51; T reg102=reg77*reg32; T reg103=reg74*reg51; T reg104=reg84*reg32;
    T reg105=reg65*reg56; T reg106=reg62*reg51; T reg107=reg72*reg49; T reg108=reg55*reg56; T reg109=reg60*reg49;
    T reg110=reg73*reg49; T reg111=reg59*reg56; T reg112=reg75*reg49; T reg113=reg74*reg56; T reg114=reg67*reg56;
    T reg115=reg78*reg42; T reg116=reg66*reg39; reg71=reg70+reg71; reg68=reg69+reg68; T reg117=reg55*reg48;
    T reg118=reg87*reg39; T reg119=reg72*reg50; T reg120=reg55*reg47; T reg121=reg65*reg54; T reg122=reg78*reg39;
    reg86=reg85+reg86; T reg123=reg83*reg50; T reg124=reg63*reg47; T reg125=reg74*reg54; T reg126=reg77*reg39;
    T reg127=reg79*reg50; T reg128=reg65*reg47; reg81=reg80+reg81; T reg129=reg61*reg32; T reg130=reg57*reg51;
    T reg131=reg82*reg50; T reg132=reg78*reg32; T reg133=reg65*reg51; T reg134=reg57*reg47; T reg135=reg60*reg46;
    T reg136=reg67*reg52; T reg137=reg67*reg53; reg66=reg66*reg42; T reg138=reg76*reg32; T reg139=reg55*reg53;
    T reg140=reg55*reg52; T reg141=reg72*reg46; T reg142=reg64*reg37; T reg143=reg63*reg54; reg60=reg60*reg40;
    T reg144=reg63*reg53; T reg145=reg83*reg46; T reg146=reg67*reg43; T reg147=reg87*reg32; T reg148=reg63*reg51;
    T reg149=reg67*reg54; reg92=reg91+reg92; reg134=reg134+reg131; reg95=reg95-reg89; reg140=reg140-reg93;
    reg105=reg97+reg105; reg94=reg115+reg94; reg127=reg128-reg127; reg66=reg66-reg136; reg123=reg124-reg123;
    reg100=reg99+reg100; reg117=reg117+reg142; reg109=reg109-reg114; reg120=reg120+reg119; reg91=reg71*reg30;
    reg108=reg108-reg107; reg60=reg60+reg146; reg90=reg90-reg88; reg143=reg118-reg143; reg96=reg96-reg98;
    reg126=reg126+reg125; reg148=reg147+reg148; reg97=reg86*reg30; reg111=reg112+reg111; reg99=reg68*reg30;
    reg135=reg135-reg137; reg112=reg81*reg30; reg130=reg130-reg129; reg116=reg116+reg149; reg102=reg102-reg103;
    reg139=reg139-reg141; reg101=reg104+reg101; reg133=reg132+reg133; reg121=reg122-reg121; reg144=reg145+reg144;
    reg110=reg110-reg113; reg106=reg138+reg106; reg60=reg30*reg60; reg116=reg116*reg30; reg148=reg148*reg30;
    reg92=reg30*reg92; reg104=ponderation*reg91; reg121=reg121*reg30; reg135=reg135*reg30; reg140=reg30*reg140;
    reg115=ponderation*reg97; reg117=reg30*reg117; reg144=reg144*reg30; reg123=reg123*reg30; reg139=reg139*reg30;
    reg66=reg30*reg66; reg143=reg143*reg30; reg95=reg30*reg95; reg133=reg133*reg30; reg110=reg30*reg110;
    reg106=reg106*reg30; reg105=reg30*reg105; reg120=reg120*reg30; reg101=reg101*reg30; reg130=reg130*reg30;
    reg100=reg30*reg100; reg134=reg134*reg30; reg102=reg102*reg30; reg108=reg30*reg108; reg111=reg111*reg30;
    reg118=ponderation*reg112; reg96=reg96*reg30; reg109=reg30*reg109; reg122=ponderation*reg99; reg127=reg127*reg30;
    reg90=reg90*reg30; reg94=reg30*reg94; reg126=reg126*reg30; T tmp_2_4=ponderation*reg121; T tmp_3_5=ponderation*reg123;
    T tmp_2_3=-reg115; T tmp_2_2=ponderation*reg126; T tmp_3_4=ponderation*reg127; T tmp_6_7=-reg118; T tmp_3_3=ponderation*reg134;
    T tmp_0_3=ponderation*reg130; T tmp_0_4=ponderation*reg133; T tmp_1_2=ponderation*reg110; T tmp_1_3=ponderation*reg95; T tmp_1_4=ponderation*reg105;
    T tmp_0_0=ponderation*reg106; T tmp_0_1=ponderation*reg101; T tmp_1_5=ponderation*reg100; T tmp_0_2=ponderation*reg102; T tmp_1_1=ponderation*reg111;
    T tmp_1_6=ponderation*reg108; T tmp_0_7=ponderation*reg96; T tmp_1_7=ponderation*reg109; T tmp_3_7=-reg104; T tmp_0_6=ponderation*reg90;
    T tmp_4_4=ponderation*reg94; T tmp_0_5=ponderation*reg148; T tmp_4_5=ponderation*reg92; T tmp_5_7=ponderation*reg135; T tmp_4_6=ponderation*reg140;
    T tmp_5_6=ponderation*reg139; T tmp_4_7=ponderation*reg66; T tmp_5_5=ponderation*reg144; T tmp_6_6=ponderation*reg117; T tmp_2_7=ponderation*reg116;
    T tmp_7_7=ponderation*reg60; T tmp_2_6=-reg122; T tmp_2_5=ponderation*reg143; T tmp_3_6=ponderation*reg120;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg3*reg4; T reg6=reg2*reg1; reg4=reg2*reg4; reg1=reg3*reg1;
    T reg7=reg2*reg4; T reg8=reg3*reg5; reg4=reg3*reg4; T reg9=reg3*reg6; T reg10=reg3*reg1;
    reg6=reg2*reg6; reg1=reg2*reg1; reg7=reg7-reg8; reg9=reg10+reg9; reg4=reg8+reg4;
    reg6=reg6-reg10; reg5=reg2*reg5; T reg11=reg3*reg4; reg9=reg3*reg9; T reg12=reg10+reg1;
    reg6=reg2*reg6; T reg13=reg2*reg7; reg5=reg8+reg5; reg12=reg3*reg12; reg9=reg6-reg9;
    reg6=reg3*reg5; reg11=reg13-reg11; reg6=reg11-reg6; reg12=reg9-reg12; reg4=reg4/reg6;
    reg12=reg12/reg6; reg7=reg7/reg6; reg8=reg4*reg12; reg9=reg7*reg12; reg11=reg4*reg8;
    reg6=reg5/reg6; reg5=1-var_inter[1]; reg13=reg7*reg9; T reg14=1-var_inter[0]; T reg15=reg7*(*f.m).alpha;
    T reg16=reg4*(*f.m).alpha; T reg17=reg5*elem.pos(0)[0]; T reg18=elem.pos(1)[0]*var_inter[0]; reg16=reg15+reg16; reg15=elem.pos(1)[1]*var_inter[0];
    T reg19=reg14*elem.pos(0)[1]; reg6=reg6*(*f.m).alpha; T reg20=reg5*elem.pos(1)[0]; T reg21=reg5*elem.pos(1)[1]; T reg22=reg5*elem.pos(0)[1];
    T reg23=reg14*elem.pos(0)[0]; reg11=reg13-reg11; reg9=reg9/reg11; reg13=elem.pos(2)[0]*var_inter[1]; T reg24=reg3*reg0;
    reg20=reg20-reg17; reg0=reg2*reg0; reg11=reg8/reg11; reg8=elem.pos(2)[1]*var_inter[0]; T reg25=reg15+reg19;
    T reg26=elem.pos(2)[1]*var_inter[1]; reg21=reg21-reg22; reg6=reg16+reg6; reg16=reg23+reg18; T reg27=elem.pos(2)[0]*var_inter[0];
    reg13=reg20+reg13; reg20=elem.pos(3)[0]*var_inter[1]; reg21=reg26+reg21; reg26=elem.pos(3)[1]*var_inter[1]; reg27=reg27-reg16;
    T reg28=reg3*reg24; T reg29=reg2*reg0; T reg30=elem.pos(3)[0]*reg14; reg9=reg9*reg6; reg6=reg11*reg6;
    reg11=elem.pos(3)[1]*reg14; reg8=reg8-reg25; reg21=reg21-reg26; reg11=reg8+reg11; reg6=reg9-reg6;
    reg8=1-(*f.m).resolution; reg13=reg13-reg20; reg30=reg27+reg30; reg28=reg29-reg28; reg9=reg13*reg11;
    reg27=reg21*reg30; reg29=(*f.m).resolution*(*f.m).alpha; reg6=reg8*reg6; reg0=reg0/reg28; reg24=reg24/reg28;
    reg27=reg9-reg27; reg6=reg29+reg6; reg9=(*f.m).resolution*reg0; reg7=reg7*reg8; reg4=reg4*reg8;
    reg29=(*f.m).resolution*reg24; reg21=reg21/reg27; reg6=(*f.m).deltaT*reg6; reg4=reg29+reg4; reg9=reg7+reg9;
    reg11=reg11/reg27; reg13=reg13/reg27; reg30=reg30/reg27; reg7=reg9*reg6; reg29=reg4*reg6;
    T reg31=reg14*reg21; T reg32=var_inter[1]*reg11; T reg33=reg5*reg30; T reg34=var_inter[0]*reg13; T reg35=var_inter[1]*reg30;
    T reg36=reg31+reg32; T reg37=reg14*reg13; T reg38=var_inter[0]*reg21; T reg39=reg33+reg34; T reg40=reg5*var_inter[0];
    T reg41=reg5*reg11; T reg42=reg29+reg7; T reg43=reg14*var_inter[1]; T reg44=reg34-reg35; T reg45=reg32-reg38;
    T reg46=reg35+reg37; T reg47=reg43*elem.f_vol_e[0]; T reg48=reg40*elem.f_vol_e[1]; T reg49=reg31-reg41; T reg50=reg42*reg39;
    T reg51=reg33-reg37; T reg52=reg42*reg36; T reg53=reg5*reg14; T reg54=reg41+reg38; T reg55=var_inter[0]*var_inter[1];
    T reg56=reg42*reg45; T reg57=reg50-reg48; T reg58=reg42*reg44; T reg59=reg55*elem.f_vol_e[0]; T reg60=reg55*elem.f_vol_e[1];
    T reg61=reg43*elem.f_vol_e[1]; T reg62=reg53*elem.f_vol_e[0]; T reg63=reg53*elem.f_vol_e[1]; T reg64=reg40*elem.f_vol_e[0]; T reg65=reg42*reg51;
    T reg66=reg42*reg54; T reg67=reg42*reg49; T reg68=reg52-reg47; T reg69=reg42*reg46; T reg70=reg65+reg63;
    T reg71=reg66+reg64; T reg72=reg67+reg62; T reg73=reg69+reg61; T reg74=reg58+reg60; T reg75=reg56+reg59;
    reg57=reg57*reg27; reg68=reg68*reg27; reg57=ponderation*reg57; T reg76=reg73*reg27; T reg77=reg74*reg27;
    T reg78=reg27*reg70; T reg79=reg27*reg72; reg68=ponderation*reg68; T reg80=reg75*reg27; T reg81=reg27*reg71;
    T reg82=ponderation*reg78; sollicitation[indices[0]+1]+=reg82; T reg83=ponderation*reg80; sollicitation[indices[2]+0]+=reg83; T reg84=ponderation*reg81;
    sollicitation[indices[1]+0]+=reg84; T reg85=ponderation*reg79; sollicitation[indices[0]+0]+=reg85; T reg86=ponderation*reg77; sollicitation[indices[2]+1]+=reg86;
    sollicitation[indices[1]+1]+=-reg57; reg57=ponderation*reg76; sollicitation[indices[3]+1]+=reg57; sollicitation[indices[3]+0]+=-reg68;
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
    T reg4=reg0*reg1; T reg5=reg2*reg4; reg4=reg3*reg4; T reg6=reg2*reg1; reg1=reg3*reg1;
    T reg7=reg3*reg5; T reg8=reg3*reg4; reg5=reg2*reg5; T reg9=reg3*reg6; T reg10=reg3*reg1;
    reg6=reg2*reg6; reg5=reg5-reg8; reg7=reg8+reg7; reg4=reg2*reg4; reg9=reg10+reg9;
    reg6=reg6-reg10; reg1=reg2*reg1; T reg11=reg10+reg1; T reg12=reg3*reg7; reg6=reg2*reg6;
    reg4=reg8+reg4; reg8=reg2*reg5; reg9=reg3*reg9; reg11=reg3*reg11; reg9=reg6-reg9;
    reg6=reg3*reg4; reg12=reg8-reg12; reg11=reg9-reg11; reg8=1-var_inter[0]; reg6=reg12-reg6;
    reg9=1-var_inter[1]; reg12=reg9*elem.pos(1)[1]; reg11=reg11/reg6; T reg13=reg9*elem.pos(0)[1]; T reg14=reg8*elem.pos(0)[1];
    T reg15=elem.pos(1)[1]*var_inter[0]; T reg16=elem.pos(1)[0]*var_inter[0]; T reg17=reg8*elem.pos(0)[0]; T reg18=reg9*elem.pos(1)[0]; reg7=reg7/reg6;
    T reg19=reg9*elem.pos(0)[0]; reg5=reg5/reg6; reg12=reg12-reg13; T reg20=elem.pos(2)[1]*var_inter[0]; reg18=reg18-reg19;
    T reg21=reg15+reg14; T reg22=elem.pos(2)[0]*var_inter[1]; T reg23=reg5*reg11; T reg24=reg7*reg11; T reg25=elem.pos(2)[1]*var_inter[1];
    T reg26=elem.pos(2)[0]*var_inter[0]; T reg27=reg17+reg16; T reg28=elem.pos(3)[1]*var_inter[1]; T reg29=elem.pos(3)[0]*var_inter[1]; reg22=reg18+reg22;
    reg12=reg25+reg12; reg18=reg7*(*f.m).alpha; reg25=reg5*(*f.m).alpha; reg6=reg4/reg6; reg4=reg7*reg24;
    T reg30=reg5*reg23; T reg31=elem.pos(3)[1]*reg8; reg20=reg20-reg21; T reg32=elem.pos(3)[0]*reg8; reg26=reg26-reg27;
    T reg33=var_inter[0]*vectors[0][indices[1]+1]; reg22=reg22-reg29; T reg34=reg8*vectors[0][indices[0]+1]; T reg35=reg8*vectors[0][indices[0]+0]; T reg36=reg9*vectors[0][indices[0]+0];
    T reg37=var_inter[0]*vectors[0][indices[1]+0]; T reg38=reg9*vectors[0][indices[1]+0]; reg6=reg6*(*f.m).alpha; reg12=reg12-reg28; reg18=reg25+reg18;
    reg25=reg9*vectors[0][indices[1]+1]; T reg39=reg9*vectors[0][indices[0]+1]; reg32=reg26+reg32; reg31=reg20+reg31; reg4=reg30-reg4;
    reg37=reg35+reg37; reg20=var_inter[0]*vectors[0][indices[2]+0]; reg26=var_inter[1]*vectors[0][indices[2]+0]; reg30=reg22*reg31; reg33=reg34+reg33;
    reg34=reg12*reg32; reg23=reg23/reg4; reg35=var_inter[0]*vectors[0][indices[2]+1]; reg6=reg18+reg6; reg4=reg24/reg4;
    reg36=reg38-reg36; reg18=var_inter[1]*vectors[0][indices[2]+1]; reg39=reg25-reg39; reg33=reg35-reg33; reg34=reg30-reg34;
    reg26=reg36+reg26; reg24=var_inter[1]*vectors[0][indices[3]+0]; reg25=var_inter[1]*vectors[0][indices[3]+1]; reg18=reg39+reg18; reg30=reg8*vectors[0][indices[3]+1];
    reg35=reg3*reg0; reg0=reg2*reg0; reg36=reg8*vectors[0][indices[3]+0]; reg23=reg23*reg6; reg6=reg4*reg6;
    reg37=reg20-reg37; reg12=reg12/reg34; reg32=reg32/reg34; reg4=pow(reg3,2); reg22=reg22/reg34;
    reg33=reg30+reg33; reg31=reg31/reg34; reg24=reg26-reg24; reg25=reg18-reg25; reg36=reg37+reg36;
    reg6=reg23-reg6; reg18=1-(*f.m).resolution; reg20=reg2*reg0; reg2=pow(reg2,2); reg3=reg3*reg35;
    reg4=reg2-reg4; reg6=reg18*reg6; reg2=(*f.m).resolution*(*f.m).alpha; reg23=reg36*reg22; reg26=reg32*reg24;
    reg30=reg25*reg31; reg3=reg20-reg3; reg20=reg12*reg33; reg0=reg0/reg3; reg33=reg22*reg33;
    reg24=reg31*reg24; reg36=reg36*reg12; reg20=reg30-reg20; reg25=reg25*reg32; reg6=reg2+reg6;
    reg26=reg23-reg26; reg4=reg4/reg3; reg3=reg35/reg3; reg2=(*f.m).resolution*reg0; reg36=reg24-reg36;
    reg5=reg5*reg18; reg6=(*f.m).deltaT*reg6; reg11=reg11*reg18; reg18=reg7*reg18; reg7=(*f.m).resolution*reg3;
    reg20=reg26+reg20; reg23=(*f.m).resolution*reg4; reg25=reg33-reg25; reg24=var_inter[1]*reg32; reg26=var_inter[0]*reg12;
    reg30=var_inter[1]*reg31; reg33=var_inter[0]*reg22; reg11=reg23+reg11; reg18=reg7+reg18; reg7=reg9*reg31;
    reg23=reg9*reg32; reg2=reg5+reg2; reg25=reg25-reg6; reg5=reg8*reg22; reg20=0.5*reg20;
    reg36=reg36-reg6; reg35=reg8*reg12; reg37=reg23+reg33; reg38=reg7+reg26; reg39=reg36*reg2;
    T reg40=reg35+reg30; T reg41=reg35-reg7; T reg42=reg24+reg5; T reg43=reg23-reg5; T reg44=reg25*reg18;
    reg36=reg36*reg18; reg20=reg20*reg11; reg25=reg25*reg2; T reg45=reg33-reg24; T reg46=reg30-reg26;
    reg39=reg44+reg39; reg44=0.5*reg42; T reg47=0.5*reg43; T reg48=0.5*reg40; T reg49=0.5*reg41;
    reg25=reg36+reg25; reg36=0.5*reg45; T reg50=0.5*reg46; T reg51=0.5*reg38; reg20=2*reg20;
    T reg52=0.5*reg37; T reg53=reg9*var_inter[0]; T reg54=reg39*reg41; T reg55=var_inter[0]*var_inter[1]; T reg56=reg9*reg8;
    T reg57=reg8*var_inter[1]; T reg58=reg20*reg47; T reg59=reg25*reg43; T reg60=reg20*reg49; T reg61=reg39*reg38;
    T reg62=reg25*reg37; T reg63=reg20*reg52; T reg64=reg20*reg50; T reg65=reg25*reg45; T reg66=reg20*reg36;
    T reg67=reg39*reg46; T reg68=reg39*reg40; T reg69=reg20*reg44; T reg70=reg25*reg42; T reg71=reg20*reg51;
    T reg72=reg20*reg48; reg64=reg65+reg64; reg65=reg57*elem.f_vol_e[1]; T reg73=reg55*elem.f_vol_e[0]; T reg74=reg56*elem.f_vol_e[1];
    reg60=reg59+reg60; reg59=reg55*elem.f_vol_e[1]; reg61=reg61-reg63; T reg75=reg53*elem.f_vol_e[0]; reg66=reg67+reg66;
    reg70=reg70-reg72; reg67=reg56*elem.f_vol_e[0]; reg58=reg54+reg58; reg69=reg69-reg68; reg54=reg57*elem.f_vol_e[0];
    T reg76=reg53*elem.f_vol_e[1]; reg71=reg71-reg62; reg61=reg61-reg75; reg66=reg66-reg73; reg71=reg71-reg76;
    reg70=reg70-reg65; reg69=reg69-reg54; reg58=reg58-reg67; reg64=reg64-reg59; reg60=reg60-reg74;
    reg60=reg34*reg60; reg64=reg34*reg64; reg61=reg34*reg61; reg58=reg34*reg58; reg66=reg34*reg66;
    reg70=reg34*reg70; reg71=reg34*reg71; reg69=reg34*reg69; sollicitation[indices[0]+1]+=ponderation*reg60; sollicitation[indices[3]+1]+=ponderation*reg70;
    sollicitation[indices[3]+0]+=ponderation*reg69; sollicitation[indices[2]+0]+=ponderation*reg66; sollicitation[indices[1]+1]+=ponderation*reg71; sollicitation[indices[2]+1]+=ponderation*reg64; sollicitation[indices[1]+0]+=ponderation*reg61;
    sollicitation[indices[0]+0]+=ponderation*reg58;
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

