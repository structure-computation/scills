
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg4*reg1; T reg6=reg3*reg2; reg2=reg4*reg2; reg1=reg3*reg1;
    T reg7=reg4*reg2; T reg8=reg4*reg5; reg2=reg3*reg2; T reg9=reg3*reg6; T reg10=reg3*reg1;
    reg5=reg3*reg5; reg2=reg2+reg9; reg7=reg7-reg9; reg8=reg8-reg10; reg5=reg10+reg5;
    reg6=reg4*reg6; reg1=reg4*reg1; reg6=reg9+reg6; reg9=reg10+reg1; T reg11=reg3*reg2;
    T reg12=reg4*reg7; reg5=reg5*reg3; reg8=reg4*reg8; reg9=reg9*reg3; T reg13=elem.pos(1)[0]-elem.pos(0)[0];
    reg11=reg12-reg11; reg12=reg3*reg6; T reg14=elem.pos(1)[1]-elem.pos(0)[1]; T reg15=elem.pos(2)[0]-elem.pos(0)[0]; T reg16=elem.pos(2)[1]-elem.pos(0)[1];
    reg5=reg8-reg5; reg8=reg13*reg16; T reg17=reg14*reg15; reg9=reg5-reg9; reg12=reg11-reg12;
    reg7=reg7/reg12; reg17=reg8-reg17; reg2=reg2/reg12; reg9=reg9/reg12; reg5=reg7*reg9;
    reg8=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg11=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; T reg18=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; T reg19=reg3*reg0; T reg20=reg2*reg9;
    reg16=reg16/reg17; T reg21=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg15=reg15/reg17; T reg22=reg4*reg0; reg13=reg13/reg17;
    reg14=reg14/reg17; T reg23=reg18*reg14; T reg24=reg8*reg13; T reg25=reg15*reg21; T reg26=reg2*reg20;
    T reg27=reg3*reg19; T reg28=reg16*reg11; reg12=reg6/reg12; reg6=reg5*reg7; T reg29=(*f.m).alpha*reg2;
    T reg30=(*f.m).alpha*reg7; T reg31=reg22*reg4; reg25=reg24-reg25; elem.epsilon[0][1]=reg25; reg23=reg28-reg23;
    elem.epsilon[0][0]=reg23; reg24=(*f.m).deltaT*(*f.m).alpha; reg29=reg30+reg29; reg12=(*f.m).alpha*reg12; reg26=reg6-reg26;
    reg27=reg31-reg27; reg6=reg25-reg24; reg28=reg23-reg24; reg5=reg5/reg26; reg22=reg22/reg27;
    reg12=reg29+reg12; reg26=reg20/reg26; reg19=reg19/reg27; reg26=reg26*reg12; reg20=reg22*reg28;
    reg12=reg5*reg12; reg5=reg6*reg22; reg28=reg19*reg28; reg6=reg6*reg19; reg26=reg12-reg26;
    reg20=reg6+reg20; reg6=1-(*f.m).resolution; reg5=reg28+reg5; reg12=(*f.m).resolution*(*f.m).alpha; reg28=pow(reg3,2);
    reg29=reg20*reg4; reg30=PNODE(2).dep[1]-PNODE(0).dep[1]; reg20=reg20*reg3; reg31=pow(reg4,2); T reg32=reg5*reg4;
    reg5=reg5*reg3; reg26=reg6*reg26; T reg33=PNODE(1).dep[0]-PNODE(0).dep[0]; T reg34=PNODE(1).dep[1]-PNODE(0).dep[1]; T reg35=PNODE(2).dep[0]-PNODE(0).dep[0];
    reg28=reg31-reg28; reg32=reg32-reg20; reg31=reg30*reg13; T reg36=reg16*reg33; T reg37=reg15*reg34;
    reg30=reg30*reg14; T reg38=reg35*reg14; reg33=reg15*reg33; reg35=reg35*reg13; reg34=reg16*reg34;
    reg29=reg29-reg5; reg26=reg12+reg26; reg11=reg15*reg11; reg18=reg18*reg13; reg29=reg29+reg24;
    reg21=reg16*reg21; reg32=reg32+reg24; reg2=reg2*reg6; reg20=reg5+reg20; reg7=reg7*reg6;
    reg8=reg8*reg14; reg5=(*f.m).resolution*reg22; reg12=(*f.m).resolution*reg19; reg27=reg28/reg27; reg37=reg31-reg37;
    reg30=reg34-reg30; reg38=reg36-reg38; reg33=reg35-reg33; reg26=(*f.m).deltaT*reg26; reg6=reg9*reg6;
    reg30=reg33+reg30; reg7=reg5+reg7; reg8=reg21-reg8; reg11=reg18-reg11; reg5=(*f.m).resolution*reg27;
    reg2=reg12+reg2; reg9=reg38-reg26; reg12=reg29+reg32; reg20=reg24-reg20; reg18=reg37-reg26;
    reg30=0.5*reg30; reg6=reg5+reg6; reg8=reg11+reg8; reg5=reg9*reg7; reg11=reg2*reg18;
    reg21=reg2*reg9; reg12=reg12+reg20; reg28=reg18*reg7; reg11=reg5+reg11; reg12=reg12/3;
    reg28=reg21+reg28; reg8=0.5*reg8; elem.epsilon[0][2]=reg8; reg5=reg6*reg30; reg29=reg29-reg12;
    reg32=reg32-reg12; reg5=2*reg5; reg9=reg9*reg11; reg21=reg27*reg8; reg18=reg18*reg28;
    reg31=reg5*reg30; reg9=reg18+reg9; reg12=reg20-reg12; reg32=pow(reg32,2); reg29=pow(reg29,2);
    reg21=reg21*reg0; reg31=reg9+reg31; reg32=reg29+reg32; reg12=pow(reg12,2); reg9=2*reg21;
    reg12=reg32+reg12; reg17=reg31*reg17; reg9=reg21*reg9; reg23=reg23-reg26; reg26=reg25-reg26;
    reg9=reg12+reg9; reg12=0.16666666666666665741*reg17; reg17=0.33333333333333331483*reg17; reg9=1.5*reg9; reg18=reg26*reg7;
    reg20=reg2*reg23; reg26=reg2*reg26; reg7=reg23*reg7; reg17=reg12+reg17; elem.sigma_von_mises=pow(reg9,0.5);
    elem.sigma[0][2]=reg6*reg8; elem.sigma[0][1]=reg20+reg18; elem.sigma[0][0]=reg7+reg26; elem.ener=reg17/2;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg2*reg4; reg4=reg3*reg4; T reg6=reg3*reg1; reg1=reg2*reg1;
    T reg7=reg2*reg5; reg5=reg3*reg5; T reg8=reg3*reg4; T reg9=reg2*reg1; T reg10=reg3*reg6;
    reg1=reg3*reg1; reg1=reg10+reg1; reg7=reg7-reg8; reg5=reg5+reg8; reg4=reg2*reg4;
    reg6=reg2*reg6; reg9=reg9-reg10; reg1=reg1*reg3; T reg11=reg3*reg5; T reg12=reg2*reg7;
    reg4=reg8+reg4; reg9=reg2*reg9; reg8=reg10+reg6; reg1=reg9-reg1; reg8=reg8*reg3;
    reg9=reg3*reg4; reg11=reg12-reg11; reg8=reg1-reg8; reg9=reg11-reg9; reg8=reg8/reg9;
    reg7=reg7/reg9; reg5=reg5/reg9; reg1=reg7*reg8; reg11=reg5*reg8; reg12=(*f.m).alpha*reg7;
    T reg13=(*f.m).alpha*reg5; T reg14=reg1*reg7; reg9=reg4/reg9; reg4=reg5*reg11; reg13=reg12+reg13;
    reg9=(*f.m).alpha*reg9; reg4=reg14-reg4; reg1=reg1/reg4; reg12=reg3*reg0; reg14=reg2*reg0;
    reg9=reg13+reg9; reg4=reg11/reg4; reg11=pow(reg2,2); reg13=reg3*reg12; T reg15=reg14*reg2;
    reg4=reg4*reg9; reg9=reg1*reg9; reg1=pow(reg3,2); T reg16=elem.pos(1)[0]-elem.pos(0)[0]; T reg17=elem.pos(1)[1]-elem.pos(0)[1];
    T reg18=elem.pos(2)[0]-elem.pos(0)[0]; T reg19=elem.pos(2)[1]-elem.pos(0)[1]; reg13=reg15-reg13; reg4=reg9-reg4; reg9=reg16*reg19;
    reg15=1-(*f.m).resolution; T reg20=reg17*reg18; reg1=reg11-reg1; reg20=reg9-reg20; reg9=(*f.m).resolution*(*f.m).alpha;
    reg4=reg15*reg4; reg1=reg1/reg13; reg12=reg12/reg13; reg13=reg14/reg13; reg11=(*f.m).resolution*reg1;
    reg8=reg8*reg15; reg5=reg5*reg15; reg15=reg7*reg15; reg7=(*f.m).resolution*reg13; reg19=reg19/reg20;
    reg14=(*f.m).resolution*reg12; reg18=reg18/reg20; reg17=reg17/reg20; reg16=reg16/reg20; reg4=reg9+reg4;
    reg5=reg14+reg5; reg15=reg7+reg15; reg8=reg11+reg8; reg4=(*f.m).deltaT*reg4; reg7=0.5*reg17;
    reg9=reg18-reg16; reg11=0.5*reg18; reg14=0.5*reg19; T reg21=0.5*reg16; T reg22=reg17-reg19;
    T reg23=reg8*reg7; T reg24=reg5*reg4; T reg25=reg8*reg11; T reg26=reg4*reg15; T reg27=reg8*reg14;
    T reg28=reg8*reg21; T reg29=0.5*reg22; T reg30=0.5*reg9; T reg31=reg26+reg24; T reg32=reg16*reg15;
    T reg33=reg5*reg16; T reg34=2*reg23; T reg35=reg19*reg5; T reg36=reg18*reg5; reg27=2*reg27;
    T reg37=reg17*reg15; reg28=2*reg28; T reg38=2*reg25; T reg39=1-var_inter[0]; T reg40=reg19*reg15;
    T reg41=reg8*reg29; T reg42=reg8*reg30; T reg43=reg18*reg15; T reg44=reg5*reg17; T reg45=reg9*reg15;
    T reg46=reg18*reg31; T reg47=reg5*reg22; T reg48=var_inter[0]*elem.f_vol_e[1]; T reg49=var_inter[1]*elem.f_vol_e[0]; T reg50=reg17*reg31;
    T reg51=reg22*reg15; T reg52=reg11*reg28; reg41=2*reg41; T reg53=reg19*reg37; T reg54=reg5*reg9;
    reg42=2*reg42; T reg55=reg43*reg16; T reg56=reg7*reg27; T reg57=reg44*reg16; T reg58=reg7*reg28;
    T reg59=reg11*reg27; T reg60=reg19*reg36; T reg61=reg18*reg35; T reg62=reg14*reg38; reg39=reg39-var_inter[1];
    T reg63=reg21*reg38; T reg64=reg40*reg17; T reg65=reg34*reg21; T reg66=reg18*reg32; T reg67=reg34*reg14;
    T reg68=reg33*reg17; reg59=reg60+reg59; T reg69=reg21*reg41; T reg70=reg44*reg9; T reg71=reg28*reg29;
    T reg72=reg37*reg17; T reg73=reg11*reg38; T reg74=reg21*reg28; T reg75=reg19*reg40; T reg76=reg51*reg17;
    T reg77=reg19*reg51; T reg78=reg21*reg42; T reg79=reg11*reg42; T reg80=var_inter[0]*elem.f_vol_e[0]; T reg81=reg11*reg41;
    T reg82=reg19*reg54; T reg83=reg18*reg47; T reg84=reg14*reg42; T reg85=reg28*reg30; T reg86=reg21*reg27;
    T reg87=reg39*elem.f_vol_e[1]; T reg88=reg14*reg41; T reg89=reg41*reg29; T reg90=reg45*reg9; T reg91=reg18*reg45;
    T reg92=reg35*reg9; T reg93=reg38*reg29; T reg94=reg42*reg29; T reg95=reg47*reg9; T reg96=reg39*elem.f_vol_e[0];
    reg64=reg63+reg64; T reg97=reg36*reg17; T reg98=reg34*reg30; T reg99=reg33*reg22; T reg100=reg43*reg9;
    reg52=reg53+reg52; reg61=reg62+reg61; T reg101=reg27*reg29; T reg102=reg54*reg17; reg28=reg14*reg28;
    reg56=reg55+reg56; T reg103=reg50-reg49; reg68=reg65+reg68; T reg104=reg16*reg31; T reg105=reg19*reg31;
    reg54=reg54*reg22; T reg106=reg27*reg30; T reg107=reg7*reg38; reg35=reg35*reg16; T reg108=reg41*reg30;
    reg40=reg40*reg22; reg41=reg7*reg41; reg45=reg45*reg16; T reg109=reg46-reg48; T reg110=reg36*reg22;
    T reg111=reg7*reg42; reg47=reg47*reg16; T reg112=var_inter[1]*elem.f_vol_e[1]; T reg113=reg34*reg29; T reg114=reg32*reg9;
    reg66=reg66+reg67; T reg115=reg22*reg31; T reg116=reg34*reg11; T reg117=reg34*reg7; reg32=reg32*reg16;
    reg51=reg51*reg22; reg33=reg19*reg33; T reg118=reg37*reg22; reg58=reg57+reg58; reg27=reg14*reg27;
    T reg119=reg18*reg43; T reg120=reg38*reg30; T reg121=reg9*reg31; reg42=reg42*reg30; T reg122=reg18*reg44;
    reg92=reg92-reg93; reg40=reg40-reg120; T reg123=reg104+reg112; reg108=reg54+reg108; reg86=reg86+reg97;
    reg106=reg106-reg110; reg74=reg72+reg74; reg101=reg101-reg100; reg85=reg85-reg118; reg42=reg51+reg42;
    reg83=reg84-reg83; reg33=reg33+reg116; reg79=reg77-reg79; reg51=reg66*reg20; reg114=reg114-reg113;
    reg71=reg71-reg70; reg54=reg96+reg115; reg77=reg20*reg58; reg84=reg87+reg121; T reg124=reg20*reg56;
    reg27=reg27+reg119; T reg125=reg80+reg105; reg35=reg35+reg107; reg41=reg45-reg41; reg109=reg20*reg109;
    reg111=reg47-reg111; reg45=reg20*reg68; reg28=reg28+reg122; reg103=reg20*reg103; reg47=reg20*reg52;
    reg102=reg69-reg102; reg99=reg99-reg98; reg91=reg88-reg91; reg69=reg20*reg59; reg94=reg95+reg94;
    reg88=reg20*reg61; reg95=reg64*reg20; reg76=reg78-reg76; reg75=reg75+reg73; reg81=reg82-reg81;
    reg89=reg90+reg89; reg32=reg32+reg117; reg78=ponderation*reg95; reg109=ponderation*reg109; reg103=ponderation*reg103;
    reg102=reg102*reg20; reg82=reg20*reg125; reg74=reg74*reg20; reg90=reg20*reg84; reg106=reg20*reg106;
    reg86=reg86*reg20; reg85=reg20*reg85; T reg126=reg20*reg54; reg76=reg76*reg20; reg32=reg20*reg32;
    T reg127=ponderation*reg51; T reg128=ponderation*reg77; T reg129=ponderation*reg69; reg91=reg20*reg91; reg99=reg20*reg99;
    T reg130=ponderation*reg124; reg94=reg20*reg94; reg27=reg20*reg27; T reg131=ponderation*reg88; reg89=reg20*reg89;
    reg92=reg20*reg92; reg35=reg20*reg35; reg71=reg20*reg71; reg101=reg20*reg101; reg75=reg20*reg75;
    reg41=reg20*reg41; reg40=reg20*reg40; reg79=reg20*reg79; reg108=reg20*reg108; reg28=reg20*reg28;
    T reg132=reg20*reg123; reg114=reg20*reg114; T reg133=ponderation*reg47; T reg134=ponderation*reg45; reg42=reg20*reg42;
    reg81=reg20*reg81; reg33=reg20*reg33; reg111=reg20*reg111; reg83=reg20*reg83; T tmp_5_4=-reg128;
    T tmp_3_2=-reg131; T tmp_1_4=ponderation*reg71; T tmp_5_5=ponderation*reg32; reg32=ponderation*reg126; sollicitation[indices[0]+0]+=reg32;
    T tmp_5_3=-reg130; reg71=ponderation*reg90; sollicitation[indices[0]+1]+=reg71; T tmp_5_2=ponderation*reg35; T tmp_5_1=ponderation*reg41;
    reg35=ponderation*reg82; sollicitation[indices[1]+0]+=reg35; T tmp_5_0=ponderation*reg111; T tmp_4_5=-reg134; sollicitation[indices[1]+1]+=-reg109;
    T tmp_3_3=ponderation*reg27; T tmp_3_4=ponderation*reg28; sollicitation[indices[2]+0]+=-reg103; T tmp_0_2=ponderation*reg40; T tmp_0_1=ponderation*reg108;
    T tmp_4_4=ponderation*reg74; reg27=ponderation*reg132; sollicitation[indices[2]+1]+=reg27; T tmp_0_0=ponderation*reg42; T tmp_0_3=ponderation*reg106;
    T tmp_2_4=-reg133; T tmp_4_3=ponderation*reg86; T tmp_2_5=ponderation*reg33; T tmp_1_3=ponderation*reg101; T tmp_0_4=ponderation*reg85;
    T tmp_3_0=ponderation*reg83; T tmp_1_2=ponderation*reg92; T tmp_1_1=ponderation*reg89; T tmp_4_2=-reg78; T tmp_1_0=ponderation*reg94;
    T tmp_0_5=ponderation*reg99; T tmp_4_1=ponderation*reg102; T tmp_3_1=ponderation*reg91; T tmp_2_3=-reg129; T tmp_4_0=ponderation*reg76;
    T tmp_2_2=ponderation*reg75; T tmp_2_1=ponderation*reg81; T tmp_3_5=-reg127; T tmp_2_0=ponderation*reg79; T tmp_1_5=ponderation*reg114;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg3*reg1; reg1=reg2*reg1; T reg6=reg2*reg4; reg4=reg3*reg4;
    T reg7=reg3*reg6; T reg8=reg3*reg1; reg6=reg2*reg6; T reg9=reg3*reg4; reg1=reg2*reg1;
    T reg10=reg3*reg5; reg1=reg1-reg10; reg5=reg2*reg5; reg4=reg2*reg4; reg7=reg7+reg9;
    reg8=reg10+reg8; reg6=reg6-reg9; T reg11=reg10+reg5; T reg12=reg2*reg6; reg4=reg9+reg4;
    reg9=reg3*reg7; reg8=reg8*reg3; reg1=reg2*reg1; reg9=reg12-reg9; reg12=reg3*reg4;
    reg8=reg1-reg8; reg11=reg11*reg3; reg11=reg8-reg11; reg12=reg9-reg12; reg7=reg7/reg12;
    reg11=reg11/reg12; reg6=reg6/reg12; reg1=reg7*reg11; reg8=reg6*reg11; reg9=reg8*reg6;
    reg12=reg4/reg12; reg4=(*f.m).alpha*reg6; T reg13=(*f.m).alpha*reg7; T reg14=reg7*reg1; reg12=(*f.m).alpha*reg12;
    reg13=reg4+reg13; reg14=reg9-reg14; reg4=reg2*reg0; reg12=reg13+reg12; reg9=reg3*reg0;
    reg8=reg8/reg14; reg14=reg1/reg14; reg1=elem.pos(1)[1]-elem.pos(0)[1]; reg13=elem.pos(2)[0]-elem.pos(0)[0]; T reg15=elem.pos(1)[0]-elem.pos(0)[0];
    T reg16=pow(reg2,2); T reg17=elem.pos(2)[1]-elem.pos(0)[1]; reg14=reg14*reg12; T reg18=reg4*reg2; T reg19=pow(reg3,2);
    reg12=reg8*reg12; reg8=reg3*reg9; reg19=reg16-reg19; reg8=reg18-reg8; reg16=reg1*reg13;
    reg18=reg15*reg17; T reg20=1-(*f.m).resolution; reg14=reg12-reg14; reg4=reg4/reg8; reg9=reg9/reg8;
    reg8=reg19/reg8; reg14=reg20*reg14; reg12=(*f.m).resolution*(*f.m).alpha; reg16=reg18-reg16; reg6=reg6*reg20;
    reg18=(*f.m).resolution*reg4; reg19=(*f.m).resolution*reg9; reg1=reg1/reg16; reg15=reg15/reg16; reg7=reg7*reg20;
    reg20=reg11*reg20; reg13=reg13/reg16; reg14=reg12+reg14; reg17=reg17/reg16; reg11=(*f.m).resolution*reg8;
    reg12=reg13-reg15; T reg21=0.5*reg17; T reg22=0.5*reg15; T reg23=0.5*reg1; T reg24=reg1-reg17;
    reg7=reg19+reg7; reg20=reg11+reg20; reg14=(*f.m).deltaT*reg14; reg6=reg18+reg6; reg11=reg20*reg22;
    reg18=0.5*reg12; reg19=0.5*reg13; T reg25=reg7*reg14; T reg26=reg20*reg23; T reg27=reg20*reg21;
    T reg28=reg14*reg6; T reg29=0.5*reg24; T reg30=reg15*reg6; reg27=2*reg27; T reg31=reg13*reg7;
    T reg32=reg20*reg19; T reg33=reg20*reg29; T reg34=reg7*reg15; T reg35=reg20*reg18; reg11=2*reg11;
    T reg36=2*reg26; T reg37=reg1*reg6; T reg38=reg28+reg25; T reg39=1-var_inter[0]; reg33=2*reg33;
    T reg40=reg1*reg38; T reg41=reg7*reg12; T reg42=reg17*reg37; T reg43=reg19*reg11; reg35=2*reg35;
    T reg44=reg13*reg38; T reg45=reg24*reg6; T reg46=reg12*reg6; T reg47=var_inter[1]*elem.f_vol_e[0]; T reg48=reg34*reg1;
    T reg49=var_inter[0]*elem.f_vol_e[1]; T reg50=reg36*reg22; T reg51=reg17*reg7; T reg52=reg7*reg1; T reg53=reg13*reg6;
    reg39=reg39-var_inter[1]; T reg54=reg36*reg21; T reg55=2*reg32; T reg56=reg13*reg30; T reg57=reg17*reg6;
    T reg58=reg17*reg31; T reg59=reg19*reg27; T reg60=reg11*reg18; T reg61=reg52*reg12; T reg62=reg27*reg29;
    T reg63=reg11*reg29; T reg64=reg39*elem.f_vol_e[0]; T reg65=reg37*reg24; T reg66=var_inter[0]*elem.f_vol_e[0]; T reg67=reg39*elem.f_vol_e[1];
    T reg68=reg27*reg18; T reg69=reg40-reg47; T reg70=reg31*reg24; T reg71=reg24*reg38; T reg72=var_inter[1]*elem.f_vol_e[1];
    T reg73=reg55*reg18; T reg74=reg45*reg24; T reg75=reg12*reg38; T reg76=reg17*reg38; T reg77=reg44-reg49;
    T reg78=reg35*reg18; T reg79=reg57*reg24; T reg80=reg33*reg18; T reg81=reg41*reg24; T reg82=reg22*reg11;
    T reg83=reg17*reg57; T reg84=reg19*reg55; T reg85=reg15*reg38; T reg86=reg36*reg29; T reg87=reg37*reg1;
    reg59=reg58+reg59; T reg88=reg30*reg12; T reg89=reg17*reg34; T reg90=reg36*reg23; reg43=reg42+reg43;
    reg30=reg30*reg15; reg34=reg34*reg24; T reg91=reg36*reg18; T reg92=reg36*reg19; T reg93=reg53*reg12;
    T reg94=reg21*reg11; T reg95=reg13*reg52; reg56=reg56+reg54; T reg96=reg13*reg53; T reg97=reg55*reg29;
    T reg98=reg51*reg12; reg48=reg50+reg48; T reg99=reg21*reg27; T reg100=reg33*reg29; T reg101=reg46*reg12;
    reg63=reg63-reg61; T reg102=reg85+reg72; reg68=reg68-reg70; T reg103=reg16*reg48; reg99=reg99+reg96;
    reg88=reg88-reg86; reg30=reg30+reg90; reg94=reg94+reg95; reg79=reg79-reg73; reg60=reg60-reg65;
    reg89=reg89+reg92; reg62=reg62-reg93; reg69=reg16*reg69; reg98=reg98-reg97; reg77=reg16*reg77;
    reg100=reg101+reg100; reg101=reg56*reg16; reg34=reg34-reg91; T reg104=reg16*reg43; reg78=reg74+reg78;
    reg74=reg66+reg76; reg83=reg83+reg84; T reg105=reg67+reg75; T reg106=reg16*reg59; reg82=reg87+reg82;
    reg80=reg81+reg80; reg81=reg64+reg71; reg99=reg16*reg99; reg69=ponderation*reg69; reg82=reg82*reg16;
    reg30=reg16*reg30; reg77=ponderation*reg77; T reg107=reg16*reg81; reg89=reg16*reg89; T reg108=ponderation*reg101;
    T reg109=reg16*reg105; T reg110=reg16*reg74; T reg111=reg16*reg102; T reg112=ponderation*reg104; T reg113=ponderation*reg106;
    reg80=reg16*reg80; reg83=reg16*reg83; reg79=reg16*reg79; reg78=reg16*reg78; reg34=reg16*reg34;
    reg88=reg16*reg88; reg68=reg16*reg68; reg100=reg16*reg100; reg63=reg16*reg63; reg60=reg16*reg60;
    reg98=reg16*reg98; reg94=reg16*reg94; reg62=reg16*reg62; T reg114=ponderation*reg103; sollicitation[indices[2]+0]+=-reg69;
    T tmp_2_3=-reg113; reg69=ponderation*reg109; sollicitation[indices[0]+1]+=reg69; T tmp_3_5=-reg108; T tmp_1_2=ponderation*reg98;
    T tmp_0_0=ponderation*reg78; reg78=ponderation*reg110; sollicitation[indices[1]+0]+=reg78; T tmp_1_3=ponderation*reg62; T tmp_0_5=ponderation*reg34;
    T tmp_1_1=ponderation*reg100; sollicitation[indices[1]+1]+=-reg77; T tmp_2_2=ponderation*reg83; reg34=ponderation*reg107; sollicitation[indices[0]+0]+=reg34;
    T tmp_0_1=ponderation*reg80; T tmp_4_4=ponderation*reg82; T tmp_0_2=ponderation*reg79; T tmp_5_5=ponderation*reg30; T tmp_1_5=ponderation*reg88;
    T tmp_0_3=ponderation*reg68; reg30=ponderation*reg111; sollicitation[indices[2]+1]+=reg30; T tmp_1_4=ponderation*reg63; T tmp_2_4=-reg112;
    T tmp_0_4=ponderation*reg60; T tmp_4_5=-reg114; T tmp_2_5=ponderation*reg89; T tmp_3_4=ponderation*reg94; T tmp_3_3=ponderation*reg99;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg4*reg2; reg2=reg3*reg2; T reg6=reg3*reg1; reg1=reg4*reg1;
    T reg7=reg3*reg5; reg5=reg4*reg5; T reg8=reg3*reg2; T reg9=reg4*reg1; T reg10=reg3*reg6;
    reg1=reg3*reg1; reg9=reg9-reg10; reg6=reg4*reg6; reg1=reg10+reg1; reg2=reg4*reg2;
    reg5=reg5-reg8; reg7=reg7+reg8; T reg11=reg3*reg0; T reg12=reg3*reg7; reg1=reg1*reg3;
    T reg13=reg4*reg5; reg2=reg8+reg2; reg8=reg4*reg0; T reg14=reg10+reg6; reg9=reg4*reg9;
    T reg15=elem.pos(2)[1]-elem.pos(0)[1]; T reg16=elem.pos(2)[0]-elem.pos(0)[0]; T reg17=elem.pos(1)[0]-elem.pos(0)[0]; T reg18=elem.pos(1)[1]-elem.pos(0)[1]; T reg19=reg3*reg11;
    T reg20=pow(reg4,2); T reg21=pow(reg3,2); T reg22=reg3*reg2; reg12=reg13-reg12; reg14=reg14*reg3;
    reg13=reg8*reg4; reg1=reg9-reg1; reg22=reg12-reg22; reg9=reg18*reg16; reg12=reg17*reg15;
    reg19=reg13-reg19; reg14=reg1-reg14; reg21=reg20-reg21; reg9=reg12-reg9; reg14=reg14/reg22;
    reg21=reg21/reg19; reg1=1-(*f.m).resolution; reg17=reg17/reg9; reg18=reg18/reg9; reg7=reg7/reg22;
    reg5=reg5/reg22; reg11=reg11/reg19; reg12=reg14*reg1; reg19=reg8/reg19; reg8=(*f.m).resolution*reg21;
    reg16=reg16/reg9; reg15=reg15/reg9; reg13=0.5*reg16; reg12=reg8+reg12; reg8=0.5*reg18;
    reg20=(*f.m).resolution*reg11; T reg23=0.5*reg17; T reg24=0.5*reg15; T reg25=reg16-reg17; T reg26=reg7*reg1;
    T reg27=(*f.m).resolution*reg19; T reg28=reg5*reg1; T reg29=reg18-reg15; T reg30=reg12*reg24; T reg31=reg12*reg23;
    reg28=reg27+reg28; reg27=reg12*reg13; reg26=reg20+reg26; reg20=reg12*reg8; T reg32=0.5*reg29;
    T reg33=0.5*reg25; T reg34=reg26*reg18; T reg35=reg16*reg26; T reg36=reg12*reg33; reg30=2*reg30;
    T reg37=reg17*reg28; T reg38=2*reg20; T reg39=reg18*reg28; reg31=2*reg31; T reg40=reg12*reg32;
    T reg41=reg15*reg28; T reg42=2*reg27; T reg43=reg26*reg17; T reg44=reg15*reg26; T reg45=reg16*reg28;
    T reg46=reg13*reg31; T reg47=reg26*reg25; reg40=2*reg40; T reg48=reg38*reg24; reg36=2*reg36;
    T reg49=reg29*reg28; T reg50=reg15*reg39; T reg51=reg13*reg30; T reg52=reg15*reg35; T reg53=reg26*reg29;
    T reg54=reg25*reg28; T reg55=reg8*reg31; T reg56=reg34*reg17; T reg57=reg8*reg30; T reg58=reg45*reg17;
    T reg59=reg41*reg18; T reg60=reg23*reg42; T reg61=reg43*reg18; T reg62=reg24*reg42; T reg63=reg16*reg44;
    T reg64=reg16*reg37; T reg65=reg38*reg23; T reg66=reg35*reg18; T reg67=reg42*reg32; T reg68=reg37*reg17;
    T reg69=reg38*reg8; reg55=reg56+reg55; T reg70=reg31*reg33; T reg71=reg45*reg25; T reg72=reg44*reg25;
    T reg73=reg30*reg32; T reg74=reg23*reg30; reg59=reg60+reg59; reg37=reg37*reg25; T reg75=reg38*reg32;
    T reg76=reg40*reg32; T reg77=reg54*reg25; reg46=reg50+reg46; reg64=reg64+reg48; reg51=reg52+reg51;
    T reg78=reg23*reg36; T reg79=reg49*reg18; T reg80=reg43*reg29; T reg81=reg23*reg40; T reg82=reg47*reg18;
    T reg83=reg13*reg42; T reg84=reg15*reg41; T reg85=reg38*reg33; T reg86=reg13*reg40; T reg87=reg15*reg47;
    T reg88=reg53*reg25; T reg89=reg36*reg32; T reg90=reg13*reg36; T reg91=reg15*reg49; T reg92=reg23*reg31;
    reg41=reg41*reg29; reg47=reg47*reg29; reg63=reg62+reg63; T reg93=reg42*reg33; T reg94=reg16*reg54;
    T reg95=reg16*reg34; T reg96=reg24*reg40; T reg97=reg24*reg31; T reg98=reg35*reg29; reg61=reg65+reg61;
    T reg99=reg36*reg33; T reg100=reg16*reg53; T reg101=reg24*reg36; reg53=reg53*reg17; reg36=reg8*reg36;
    reg49=reg49*reg29; reg57=reg58+reg57; T reg102=reg39*reg18; T reg103=reg16*reg45; T reg104=reg39*reg29;
    T reg105=reg34*reg25; T reg106=reg40*reg33; reg31=reg31*reg32; T reg107=reg8*reg42; reg44=reg44*reg17;
    reg43=reg15*reg43; T reg108=reg24*reg30; reg40=reg8*reg40; reg54=reg54*reg17; T reg109=reg38*reg13;
    reg30=reg30*reg33; T reg110=reg9*reg46; reg106=reg47+reg106; reg68=reg68+reg69; reg99=reg49+reg99;
    reg80=reg80-reg85; reg41=reg41-reg93; reg89=reg88+reg89; reg76=reg77+reg76; reg30=reg30-reg98;
    reg72=reg72-reg67; reg70=reg70-reg104; reg73=reg73-reg71; reg74=reg74+reg66; reg40=reg54-reg40;
    reg36=reg53-reg36; reg44=reg44+reg107; reg47=reg9*reg61; reg49=reg59*reg9; reg53=reg9*reg57;
    reg54=reg9*reg55; reg97=reg97+reg95; reg82=reg81-reg82; reg31=reg31-reg105; reg108=reg108+reg103;
    reg92=reg102+reg92; reg37=reg37-reg75; reg79=reg78-reg79; reg43=reg43+reg109; reg77=reg9*reg51;
    reg100=reg101-reg100; reg84=reg84+reg83; reg78=reg64*reg9; reg86=reg87-reg86; reg81=reg9*reg63;
    reg94=reg96-reg94; reg90=reg91-reg90; reg43=reg9*reg43; reg36=reg9*reg36; reg94=reg9*reg94;
    reg30=reg9*reg30; reg87=ponderation*reg47; reg88=ponderation*reg81; reg70=reg9*reg70; reg92=reg92*reg9;
    reg99=reg9*reg99; reg74=reg74*reg9; reg100=reg9*reg100; reg97=reg9*reg97; reg41=reg9*reg41;
    reg108=reg9*reg108; reg106=reg9*reg106; reg91=ponderation*reg110; reg96=ponderation*reg77; reg68=reg9*reg68;
    reg84=reg9*reg84; reg86=reg9*reg86; reg80=reg9*reg80; reg101=ponderation*reg78; reg90=reg9*reg90;
    reg89=reg9*reg89; reg37=reg9*reg37; reg79=reg79*reg9; reg76=reg9*reg76; reg31=reg9*reg31;
    reg40=reg9*reg40; reg44=reg9*reg44; T reg111=ponderation*reg49; reg73=reg9*reg73; T reg112=ponderation*reg53;
    reg82=reg82*reg9; reg72=reg9*reg72; T reg113=ponderation*reg54; T tmp_3_5=-reg101; T tmp_4_4=ponderation*reg92;
    T tmp_4_1=ponderation*reg82; T tmp_4_3=ponderation*reg74; T tmp_0_3=ponderation*reg30; T tmp_4_0=ponderation*reg79; T tmp_0_4=ponderation*reg70;
    T tmp_4_2=-reg111; T tmp_5_5=ponderation*reg68; T tmp_2_3=-reg96; T tmp_2_2=ponderation*reg84; T tmp_2_1=ponderation*reg86;
    T tmp_0_5=ponderation*reg80; T tmp_2_0=ponderation*reg90; T tmp_1_0=ponderation*reg89; T tmp_1_5=ponderation*reg37; T tmp_1_4=ponderation*reg31;
    T tmp_1_1=ponderation*reg76; T tmp_5_4=-reg113; T tmp_1_2=ponderation*reg72; T tmp_5_3=-reg112; T tmp_1_3=ponderation*reg73;
    T tmp_5_2=ponderation*reg44; T tmp_5_1=ponderation*reg40; T tmp_5_0=ponderation*reg36; T tmp_4_5=-reg87; T tmp_3_4=ponderation*reg97;
    T tmp_0_0=ponderation*reg99; T tmp_3_3=ponderation*reg108; T tmp_3_2=-reg88; T tmp_0_1=ponderation*reg106; T tmp_3_1=ponderation*reg94;
    T tmp_3_0=ponderation*reg100; T tmp_0_2=ponderation*reg41; T tmp_2_5=ponderation*reg43; T tmp_2_4=-reg91;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg3*reg1; reg1=reg4*reg1; T reg6=reg3*reg2; reg2=reg4*reg2;
    T reg7=reg4*reg2; reg2=reg3*reg2; T reg8=reg3*reg6; T reg9=reg4*reg1; reg1=reg3*reg1;
    T reg10=reg3*reg5; reg5=reg4*reg5; reg7=reg7-reg8; reg9=reg9-reg10; reg2=reg2+reg8;
    reg1=reg10+reg1; reg6=reg4*reg6; T reg11=reg4*reg0; T reg12=reg3*reg2; T reg13=reg10+reg5;
    T reg14=reg4*reg7; reg9=reg4*reg9; reg6=reg8+reg6; reg1=reg1*reg3; reg8=reg3*reg0;
    T reg15=reg3*reg6; T reg16=elem.pos(2)[1]-elem.pos(0)[1]; T reg17=reg11*reg4; T reg18=elem.pos(2)[0]-elem.pos(0)[0]; T reg19=elem.pos(1)[1]-elem.pos(0)[1];
    T reg20=elem.pos(1)[0]-elem.pos(0)[0]; reg12=reg14-reg12; reg13=reg13*reg3; reg1=reg9-reg1; reg9=reg3*reg8;
    reg14=pow(reg3,2); T reg21=pow(reg4,2); reg15=reg12-reg15; reg12=reg19*reg18; T reg22=reg20*reg16;
    reg9=reg17-reg9; reg14=reg21-reg14; reg13=reg1-reg13; reg12=reg22-reg12; reg1=1-(*f.m).resolution;
    reg13=reg13/reg15; reg14=reg14/reg9; reg19=reg19/reg12; reg17=reg13*reg1; reg11=reg11/reg9;
    reg21=(*f.m).resolution*reg14; reg20=reg20/reg12; reg9=reg8/reg9; reg7=reg7/reg15; reg2=reg2/reg15;
    reg16=reg16/reg12; reg18=reg18/reg12; reg8=reg2*reg1; reg22=reg7*reg1; T reg23=(*f.m).resolution*reg11;
    T reg24=(*f.m).resolution*reg9; T reg25=reg18-reg20; T reg26=0.5*reg16; T reg27=0.5*reg20; T reg28=0.5*reg19;
    T reg29=reg19-reg16; reg17=reg21+reg17; reg21=0.5*reg18; T reg30=0.5*reg25; T reg31=0.5*reg29;
    T reg32=reg17*reg28; T reg33=reg17*reg26; T reg34=reg17*reg27; reg8=reg24+reg8; reg22=reg23+reg22;
    reg34=2*reg34; reg23=2*reg32; reg24=reg19*reg22; T reg35=reg20*reg22; reg33=2*reg33;
    T reg36=reg8*reg20; T reg37=reg18*reg8; T reg38=reg17*reg31; T reg39=reg17*reg30; T reg40=reg17*reg21;
    T reg41=reg23*reg27; reg39=2*reg39; T reg42=reg29*reg22; T reg43=reg8*reg19; T reg44=reg16*reg37;
    T reg45=reg21*reg33; T reg46=reg18*reg22; T reg47=reg16*reg24; T reg48=reg21*reg34; T reg49=reg16*reg8;
    T reg50=reg25*reg22; T reg51=reg8*reg25; reg38=2*reg38; T reg52=reg23*reg26; T reg53=reg36*reg19;
    T reg54=reg16*reg22; T reg55=2*reg40; T reg56=reg18*reg35; T reg57=reg34*reg30; T reg58=reg50*reg25;
    T reg59=reg38*reg31; T reg60=reg24*reg29; T reg61=reg33*reg30; T reg62=reg49*reg25; T reg63=reg55*reg31;
    T reg64=reg37*reg29; T reg65=reg55*reg30; T reg66=reg46*reg25; T reg67=reg33*reg31; T reg68=reg23*reg21;
    T reg69=reg54*reg29; T reg70=reg43*reg25; T reg71=reg34*reg31; T reg72=reg16*reg36; T reg73=reg38*reg30;
    T reg74=reg42*reg29; T reg75=reg51*reg29; T reg76=reg39*reg30; T reg77=reg27*reg34; T reg78=reg16*reg54;
    T reg79=reg23*reg31; T reg80=reg35*reg25; T reg81=reg21*reg55; T reg82=reg23*reg28; T reg83=reg18*reg46;
    reg45=reg44+reg45; reg35=reg35*reg20; reg53=reg41+reg53; T reg84=reg26*reg33; reg56=reg56+reg52;
    T reg85=reg18*reg43; reg48=reg47+reg48; T reg86=reg24*reg19; reg36=reg36*reg29; T reg87=reg26*reg34;
    T reg88=reg23*reg30; reg80=reg80-reg79; reg87=reg87+reg85; reg76=reg74+reg76; reg77=reg86+reg77;
    reg71=reg71-reg70; reg57=reg57-reg60; reg61=reg61-reg64; reg73=reg75+reg73; reg74=reg12*reg53;
    reg75=reg56*reg12; reg69=reg69-reg65; reg72=reg72+reg68; T reg89=reg12*reg45; reg62=reg62-reg63;
    reg67=reg67-reg66; T reg90=reg12*reg48; reg78=reg78+reg81; reg59=reg58+reg59; reg35=reg35+reg82;
    reg36=reg36-reg88; reg84=reg84+reg83; reg59=reg12*reg59; reg61=reg12*reg61; reg72=reg12*reg72;
    reg58=ponderation*reg75; reg87=reg12*reg87; reg36=reg12*reg36; reg62=reg12*reg62; reg57=reg12*reg57;
    reg69=reg12*reg69; T reg91=ponderation*reg89; T reg92=ponderation*reg74; reg77=reg77*reg12; reg73=reg12*reg73;
    reg67=reg12*reg67; reg78=reg12*reg78; reg71=reg12*reg71; reg84=reg12*reg84; T reg93=ponderation*reg90;
    reg76=reg12*reg76; reg35=reg12*reg35; reg80=reg12*reg80; T tmp_0_4=ponderation*reg57; T tmp_5_5=ponderation*reg35;
    T tmp_3_5=-reg58; T tmp_4_4=ponderation*reg77; T tmp_2_4=-reg93; T tmp_2_5=ponderation*reg72; T tmp_3_3=ponderation*reg84;
    T tmp_3_4=ponderation*reg87; T tmp_4_5=-reg92; T tmp_1_4=ponderation*reg71; T tmp_1_5=ponderation*reg80; T tmp_2_2=ponderation*reg78;
    T tmp_2_3=-reg91; T tmp_0_5=ponderation*reg36; T tmp_1_1=ponderation*reg59; T tmp_1_2=ponderation*reg62; T tmp_1_3=ponderation*reg67;
    T tmp_0_0=ponderation*reg76; T tmp_0_1=ponderation*reg73; T tmp_0_2=ponderation*reg69; T tmp_0_3=ponderation*reg61;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg2*reg1; reg1=reg3*reg1; T reg6=reg3*reg4; reg4=reg2*reg4;
    T reg7=reg3*reg4; T reg8=reg2*reg5; reg5=reg3*reg5; T reg9=reg3*reg6; T reg10=reg3*reg1;
    reg4=reg2*reg4; reg8=reg8-reg10; reg1=reg2*reg1; reg5=reg10+reg5; reg4=reg4-reg9;
    reg7=reg7+reg9; reg6=reg2*reg6; T reg11=reg2*reg4; T reg12=reg3*reg7; reg8=reg2*reg8;
    reg6=reg9+reg6; reg5=reg5*reg3; reg9=reg10+reg1; T reg13=reg3*reg6; reg12=reg11-reg12;
    reg5=reg8-reg5; reg9=reg9*reg3; reg13=reg12-reg13; reg9=reg5-reg9; reg7=reg7/reg13;
    reg4=reg4/reg13; reg9=reg9/reg13; reg5=reg4*reg9; reg8=reg7*reg9; reg11=reg5*reg4;
    reg13=reg6/reg13; reg6=reg7*reg8; reg12=(*f.m).alpha*reg7; T reg14=(*f.m).alpha*reg4; reg13=(*f.m).alpha*reg13;
    reg6=reg11-reg6; reg12=reg14+reg12; reg11=reg3*reg0; reg5=reg5/reg6; reg14=reg2*reg0;
    reg13=reg12+reg13; reg6=reg8/reg6; reg8=reg3*reg11; reg12=reg14*reg2; reg5=reg5*reg13;
    reg13=reg6*reg13; reg8=reg12-reg8; reg6=1-(*f.m).resolution; reg13=reg5-reg13; reg14=reg14/reg8;
    reg5=(*f.m).resolution*(*f.m).alpha; reg11=reg11/reg8; reg13=reg6*reg13; reg12=elem.pos(2)[0]-elem.pos(0)[0]; T reg15=elem.pos(1)[1]-elem.pos(0)[1];
    T reg16=elem.pos(1)[0]-elem.pos(0)[0]; T reg17=elem.pos(2)[1]-elem.pos(0)[1]; reg13=reg5+reg13; reg5=(*f.m).resolution*reg14; T reg18=(*f.m).resolution*reg11;
    reg4=reg4*reg6; reg7=reg7*reg6; T reg19=reg16*reg17; reg13=(*f.m).deltaT*reg13; T reg20=reg15*reg12;
    reg7=reg18+reg7; reg4=reg5+reg4; reg20=reg19-reg20; reg5=reg7*reg13; reg18=reg13*reg4;
    reg16=reg16/reg20; reg15=reg15/reg20; reg19=reg18+reg5; T reg21=1-var_inter[0]; reg12=reg12/reg20;
    reg17=reg17/reg20; T reg22=var_inter[0]*elem.f_vol_e[1]; T reg23=var_inter[1]*elem.f_vol_e[0]; T reg24=reg12*reg19; T reg25=reg15-reg17;
    T reg26=reg15*reg19; T reg27=reg12-reg16; reg21=reg21-var_inter[1]; T reg28=reg24-reg22; T reg29=reg25*reg19;
    T reg30=reg27*reg19; T reg31=reg17*reg19; T reg32=reg16*reg19; T reg33=reg26-reg23; T reg34=var_inter[0]*elem.f_vol_e[0];
    T reg35=var_inter[1]*elem.f_vol_e[1]; T reg36=reg21*elem.f_vol_e[0]; T reg37=reg21*elem.f_vol_e[1]; reg28=reg20*reg28; reg33=reg20*reg33;
    T reg38=reg34+reg31; T reg39=reg37+reg30; T reg40=reg36+reg29; T reg41=reg32+reg35; T reg42=reg20*reg39;
    T reg43=reg20*reg38; T reg44=reg20*reg40; reg28=ponderation*reg28; T reg45=reg20*reg41; reg33=ponderation*reg33;
    T reg46=ponderation*reg42; sollicitation[indices[0]+1]+=reg46; T reg47=ponderation*reg44; sollicitation[indices[0]+0]+=reg47; T reg48=ponderation*reg43;
    sollicitation[indices[1]+0]+=reg48; sollicitation[indices[1]+1]+=-reg28; reg28=ponderation*reg45; sollicitation[indices[2]+1]+=reg28; sollicitation[indices[2]+0]+=-reg33;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg2*reg4; reg4=reg3*reg4; T reg6=reg2*reg1; reg1=reg3*reg1;
    T reg7=reg3*reg4; T reg8=reg3*reg5; T reg9=reg2*reg6; T reg10=reg3*reg1; reg5=reg2*reg5;
    reg6=reg3*reg6; reg6=reg10+reg6; reg9=reg9-reg10; reg5=reg5-reg7; reg8=reg8+reg7;
    reg4=reg2*reg4; reg1=reg2*reg1; reg9=reg2*reg9; T reg11=reg3*reg8; T reg12=reg2*reg5;
    reg4=reg7+reg4; reg6=reg6*reg3; reg7=reg10+reg1; reg6=reg9-reg6; reg9=reg3*reg4;
    reg7=reg7*reg3; reg11=reg12-reg11; reg9=reg11-reg9; reg7=reg6-reg7; reg7=reg7/reg9;
    reg5=reg5/reg9; reg8=reg8/reg9; reg6=reg5*reg7; reg11=reg8*reg7; reg12=reg8*reg11;
    T reg13=(*f.m).alpha*reg8; T reg14=(*f.m).alpha*reg5; reg9=reg4/reg9; reg4=reg6*reg5; reg9=(*f.m).alpha*reg9;
    reg13=reg14+reg13; reg12=reg4-reg12; reg4=elem.pos(1)[0]-elem.pos(0)[0]; reg14=elem.pos(1)[1]-elem.pos(0)[1]; T reg15=elem.pos(2)[1]-elem.pos(0)[1];
    T reg16=elem.pos(2)[0]-elem.pos(0)[0]; reg11=reg11/reg12; reg12=reg6/reg12; reg9=reg13+reg9; reg6=reg14*reg16;
    reg13=reg4*reg15; T reg17=reg3*reg0; reg11=reg11*reg9; reg9=reg12*reg9; reg6=reg13-reg6;
    reg12=reg2*reg0; reg11=reg9-reg11; reg9=1-(*f.m).resolution; reg13=reg3*reg17; T reg18=pow(reg2,2);
    T reg19=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; T reg20=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg16=reg16/reg6; T reg21=reg12*reg2; T reg22=vectors[0][indices[2]+0]-vectors[0][indices[0]+0];
    reg15=reg15/reg6; T reg23=pow(reg3,2); reg14=reg14/reg6; reg4=reg4/reg6; T reg24=vectors[0][indices[1]+0]-vectors[0][indices[0]+0];
    reg13=reg21-reg13; reg21=reg20*reg14; T reg25=reg15*reg19; T reg26=reg16*reg24; T reg27=reg22*reg4;
    reg11=reg9*reg11; T reg28=(*f.m).resolution*(*f.m).alpha; reg23=reg18-reg23; reg17=reg17/reg13; reg11=reg28+reg11;
    reg24=reg15*reg24; reg22=reg22*reg14; reg20=reg20*reg4; reg19=reg16*reg19; reg12=reg12/reg13;
    reg13=reg23/reg13; reg26=reg27-reg26; reg21=reg25-reg21; reg11=(*f.m).deltaT*reg11; reg7=reg7*reg9;
    reg21=reg26+reg21; reg8=reg8*reg9; reg22=reg24-reg22; reg18=(*f.m).resolution*reg13; reg23=(*f.m).resolution*reg17;
    reg9=reg5*reg9; reg19=reg20-reg19; reg5=(*f.m).resolution*reg12; reg8=reg23+reg8; reg7=reg18+reg7;
    reg9=reg5+reg9; reg21=0.5*reg21; reg22=reg22-reg11; reg19=reg19-reg11; reg5=reg19*reg9;
    reg21=reg7*reg21; reg18=reg16-reg4; reg20=reg14-reg15; reg23=reg8*reg22; reg22=reg22*reg9;
    reg19=reg8*reg19; reg24=0.5*reg20; reg25=0.5*reg18; reg26=0.5*reg16; reg27=0.5*reg14;
    reg19=reg22+reg19; reg5=reg23+reg5; reg22=0.5*reg4; reg23=1-var_inter[0]; reg28=0.5*reg15;
    reg21=2*reg21; T reg29=reg5*reg18; T reg30=reg21*reg24; T reg31=reg21*reg26; T reg32=reg21*reg25;
    T reg33=reg16*reg5; T reg34=reg19*reg20; T reg35=reg21*reg28; T reg36=reg15*reg19; T reg37=reg21*reg27;
    T reg38=reg5*reg4; T reg39=reg19*reg14; T reg40=reg21*reg22; reg23=reg23-var_inter[1]; reg36=reg36-reg31;
    T reg41=reg23*elem.f_vol_e[1]; reg30=reg29+reg30; reg40=reg40-reg39; reg29=var_inter[1]*elem.f_vol_e[0]; reg38=reg38-reg37;
    T reg42=var_inter[1]*elem.f_vol_e[1]; T reg43=var_inter[0]*elem.f_vol_e[1]; reg35=reg35-reg33; reg32=reg34+reg32; reg34=reg23*elem.f_vol_e[0];
    T reg44=var_inter[0]*elem.f_vol_e[0]; reg35=reg35-reg43; reg38=reg38-reg42; reg30=reg30-reg41; reg32=reg32-reg34;
    reg40=reg40-reg29; reg36=reg36-reg44; reg36=reg6*reg36; reg35=reg6*reg35; reg40=reg6*reg40;
    reg38=reg6*reg38; reg30=reg6*reg30; reg32=reg6*reg32; sollicitation[indices[1]+1]+=ponderation*reg35; sollicitation[indices[2]+0]+=ponderation*reg40;
    sollicitation[indices[0]+1]+=ponderation*reg30; sollicitation[indices[2]+1]+=ponderation*reg38; sollicitation[indices[1]+0]+=ponderation*reg36; sollicitation[indices[0]+0]+=ponderation*reg32;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg4*reg2; reg2=reg3*reg2; T reg6=reg4*reg1; reg1=reg3*reg1;
    T reg7=reg3*reg5; T reg8=reg3*reg1; T reg9=reg3*reg6; T reg10=reg3*reg2; reg5=reg4*reg5;
    reg6=reg4*reg6; reg5=reg5-reg10; reg7=reg10+reg7; reg1=reg4*reg1; reg6=reg6-reg8;
    reg2=reg4*reg2; reg9=reg8+reg9; T reg11=0.5*elem.pos(0)[0]; T reg12=reg8+reg1; reg9=reg3*reg9;
    reg6=reg4*reg6; reg2=reg10+reg2; reg10=0.5*elem.pos(1)[0]; T reg13=0.5*elem.pos(1)[1]; T reg14=0.5*elem.pos(0)[1];
    T reg15=reg3*reg7; T reg16=reg4*reg5; T reg17=reg13-reg14; T reg18=reg10+reg11; T reg19=0.5*elem.pos(2)[0];
    reg11=reg10-reg11; reg10=reg3*reg2; reg15=reg16-reg15; reg13=reg14+reg13; reg9=reg6-reg9;
    reg12=reg3*reg12; reg6=0.5*elem.pos(2)[1]; reg13=reg6-reg13; reg11=reg11+reg19; reg14=0.5*elem.pos(3)[0];
    reg16=0.5*elem.pos(3)[1]; reg17=reg6+reg17; reg18=reg19-reg18; reg10=reg15-reg10; reg12=reg9-reg12;
    reg7=reg7/reg10; reg13=reg16+reg13; reg6=0.5*vectors[0][indices[0]+1]; reg9=0.5*vectors[0][indices[1]+1]; reg12=reg12/reg10;
    reg11=reg11-reg14; reg5=reg5/reg10; reg18=reg14+reg18; reg16=reg17-reg16; reg14=0.5*vectors[0][indices[0]+0];
    reg15=0.5*vectors[0][indices[1]+0]; reg17=0.21132486540518713447*elem.pos(1)[1]; reg19=0.78867513459481286553*elem.pos(0)[0]; T reg20=0.78867513459481286553*elem.pos(1)[0]; T reg21=reg5*reg12;
    T reg22=reg7*reg12; T reg23=0.78867513459481286553*elem.pos(1)[1]; T reg24=0.21132486540518713447*elem.pos(0)[1]; T reg25=0.21132486540518713447*elem.pos(0)[0]; T reg26=0.21132486540518713447*elem.pos(1)[0];
    T reg27=0.78867513459481286553*elem.pos(0)[1]; T reg28=reg9-reg6; T reg29=reg13*reg11; reg6=reg9+reg6; reg9=reg16*reg18;
    T reg30=0.5*vectors[0][indices[2]+1]; T reg31=0.5*vectors[0][indices[2]+0]; T reg32=reg15-reg14; reg14=reg15+reg14; reg28=reg30+reg28;
    reg6=reg30-reg6; reg15=0.5*vectors[0][indices[3]+1]; reg14=reg31-reg14; reg30=0.5*vectors[0][indices[3]+0]; reg10=reg2/reg10;
    reg2=(*f.m).alpha*reg7; T reg33=(*f.m).alpha*reg5; reg32=reg31+reg32; reg31=reg7*reg22; T reg34=reg5*reg21;
    reg9=reg29-reg9; reg29=reg26+reg19; T reg35=0.21132486540518713447*elem.pos(2)[0]; T reg36=reg17+reg27; T reg37=0.21132486540518713447*elem.pos(2)[1];
    T reg38=reg24+reg23; reg24=reg17-reg24; reg26=reg26-reg25; reg17=0.78867513459481286553*elem.pos(2)[0]; T reg39=0.78867513459481286553*elem.pos(2)[1];
    reg25=reg25+reg20; reg14=reg30+reg14; reg16=reg16/reg9; reg25=reg17-reg25; reg30=reg32-reg30;
    reg27=reg23-reg27; reg23=reg4*reg0; reg19=reg20-reg19; reg13=reg13/reg9; reg29=reg35-reg29;
    reg24=reg39+reg24; reg20=0.78867513459481286553*elem.pos(3)[1]; reg36=reg37-reg36; reg18=reg18/reg9; reg9=reg11/reg9;
    reg6=reg15+reg6; reg17=reg26+reg17; reg11=0.78867513459481286553*elem.pos(3)[0]; reg26=0.21132486540518713447*elem.pos(3)[0]; reg15=reg28-reg15;
    reg10=(*f.m).alpha*reg10; reg2=reg33+reg2; reg28=0.21132486540518713447*elem.pos(3)[1]; reg31=reg34-reg31; reg32=reg3*reg0;
    reg38=reg39-reg38; reg25=reg25+reg26; reg33=0.78867513459481286553*PNODE(1).dep[1]; reg34=0.21132486540518713447*PNODE(0).dep[1]; reg19=reg35+reg19;
    reg35=reg13*reg30; reg27=reg37+reg27; reg37=reg16*reg14; reg39=0.21132486540518713447*PNODE(1).dep[1]; T reg40=0.21132486540518713447*PNODE(1).dep[0];
    T reg41=0.21132486540518713447*PNODE(0).dep[0]; T reg42=reg6*reg9; reg38=reg28+reg38; reg21=reg21/reg31; T reg43=0.78867513459481286553*PNODE(0).dep[1];
    reg31=reg22/reg31; reg10=reg2+reg10; reg2=reg3*reg32; reg17=reg17-reg11; reg29=reg11+reg29;
    reg36=reg20+reg36; reg20=reg24-reg20; reg11=0.78867513459481286553*PNODE(0).dep[0]; reg22=0.78867513459481286553*PNODE(1).dep[0]; reg24=reg4*reg23;
    T reg44=reg15*reg18; reg31=reg31*reg10; reg10=reg21*reg10; reg21=reg40+reg11; T reg45=0.78867513459481286553*PNODE(2).dep[0];
    T reg46=(*f.m).deltaT*(*f.m).alpha; T reg47=reg39-reg34; reg37=reg35-reg37; elem.epsilon[0][0]=reg37; reg40=reg40-reg41;
    reg35=0.21132486540518713447*PNODE(2).dep[0]; reg28=reg27-reg28; reg39=reg39+reg43; reg27=reg17*reg38; T reg48=0.21132486540518713447*PNODE(2).dep[1];
    reg2=reg24-reg2; reg26=reg19-reg26; reg19=reg20*reg29; reg34=reg33+reg34; reg24=reg17*reg36;
    T reg49=0.78867513459481286553*PNODE(2).dep[1]; reg44=reg42-reg44; elem.epsilon[0][1]=reg44; reg42=reg20*reg25; reg41=reg41+reg22;
    T reg50=0.21132486540518713447*PNODE(3).dep[0]; T reg51=reg37-reg46; T reg52=reg44-reg46; reg41=reg45-reg41; reg23=reg23/reg2;
    reg32=reg32/reg2; T reg53=1-(*f.m).resolution; reg31=reg10-reg31; reg10=0.78867513459481286553*PNODE(3).dep[0]; reg42=reg27-reg42;
    reg27=0.21132486540518713447*PNODE(3).dep[1]; reg34=reg49-reg34; reg47=reg49+reg47; reg49=0.78867513459481286553*PNODE(3).dep[1]; reg40=reg45+reg40;
    reg21=reg35-reg21; reg45=reg38*reg26; reg43=reg33-reg43; reg19=reg24-reg19; reg11=reg22-reg11;
    reg22=reg25*reg28; reg39=reg48-reg39; reg24=reg25/reg42; reg33=reg20/reg42; T reg54=pow(reg4,2);
    reg35=reg11+reg35; reg11=pow(reg3,2); reg39=reg49+reg39; T reg55=reg17/reg19; T reg56=reg29/reg19;
    T reg57=reg36/reg19; reg20=reg20/reg19; T reg58=reg38/reg42; reg40=reg40-reg10; T reg59=(*f.m).resolution*(*f.m).alpha;
    reg31=reg53*reg31; reg21=reg10+reg21; reg10=reg36*reg26; T reg60=reg29*reg28; reg41=reg41+reg50;
    reg43=reg48+reg43; reg48=reg51*reg23; reg22=reg45-reg22; reg45=reg52*reg23; reg17=reg17/reg42;
    reg52=reg52*reg32; reg34=reg27+reg34; reg51=reg51*reg32; reg49=reg47-reg49; reg47=reg40*reg56;
    T reg61=reg55*reg21; T reg62=reg39*reg20; T reg63=reg33*reg41; T reg64=reg49*reg57; reg60=reg10-reg60;
    reg10=reg49*reg24; reg57=reg40*reg57; reg25=reg25/reg22; T reg65=reg26/reg22; reg20=reg20*reg21;
    reg48=reg52+reg48; reg31=reg59+reg31; reg50=reg35-reg50; reg33=reg34*reg33; reg45=reg51+reg45;
    reg35=reg34*reg17; reg11=reg54-reg11; reg51=reg28/reg22; reg38=reg38/reg22; reg52=reg58*reg40;
    reg58=reg49*reg58; reg55=reg39*reg55; reg27=reg43-reg27; reg17=reg17*reg41; reg56=reg49*reg56;
    reg40=reg24*reg40; reg24=reg34*reg65; reg20=reg57-reg20; reg43=reg27*reg25; reg49=reg38*reg50;
    reg25=reg50*reg25; reg65=reg41*reg65; reg34=reg51*reg34; reg38=reg38*reg27; reg10=reg35-reg10;
    reg41=reg51*reg41; reg62=reg64-reg62; reg47=reg61-reg47; reg35=reg4*reg48; reg63=reg52-reg63;
    reg51=reg3*reg45; reg33=reg58-reg33; reg45=reg4*reg45; reg48=reg3*reg48; reg40=reg17-reg40;
    reg2=reg11/reg2; reg11=(*f.m).resolution*reg32; reg56=reg55-reg56; reg7=reg7*reg53; reg28=reg28/reg60;
    reg5=reg5*reg53; reg26=reg26/reg60; reg36=reg36/reg60; reg29=reg29/reg60; reg17=(*f.m).resolution*reg23;
    reg31=(*f.m).deltaT*reg31; reg63=reg63-reg31; reg52=reg26*reg39; reg35=reg35-reg51; reg54=reg29*reg27;
    reg20=reg20-reg31; reg45=reg45-reg48; reg43=reg24-reg43; reg7=reg11+reg7; reg40=reg33+reg40;
    reg29=reg50*reg29; reg26=reg26*reg21; reg39=reg28*reg39; reg25=reg65-reg25; reg27=reg36*reg27;
    reg11=(*f.m).resolution*reg2; reg41=reg49-reg41; reg47=reg62+reg47; reg17=reg5+reg17; reg34=reg38-reg34;
    reg53=reg12*reg53; reg56=reg56-reg31; reg21=reg28*reg21; reg36=reg50*reg36; reg10=reg10-reg31;
    reg9=reg14*reg9; reg47=0.5*reg47; reg5=reg7*reg56; reg12=reg17*reg20; reg14=reg17*reg56;
    reg24=reg7*reg20; reg40=0.5*reg40; reg28=reg7*reg10; reg33=reg17*reg63; reg38=reg17*reg10;
    reg49=reg7*reg63; reg53=reg11+reg53; reg18=reg30*reg18; reg35=reg35+reg46; reg45=reg45+reg46;
    reg25=reg34+reg25; reg15=reg13*reg15; reg43=reg43-reg31; reg54=reg52-reg54; reg41=reg41-reg31;
    reg29=reg26-reg29; reg6=reg16*reg6; reg21=reg36-reg21; reg48=reg51+reg48; reg39=reg27-reg39;
    reg54=reg54-reg31; reg11=reg35+reg45; reg13=reg53*reg40; reg5=reg12+reg5; reg25=0.5*reg25;
    reg6=reg15-reg6; reg48=reg46-reg48; reg14=reg24+reg14; reg12=reg43*reg7; reg29=reg39+reg29;
    reg28=reg33+reg28; reg15=reg41*reg17; reg38=reg49+reg38; reg16=reg43*reg17; reg24=reg41*reg7;
    reg21=reg21-reg31; reg26=reg53*reg47; reg18=reg9-reg18; reg29=0.5*reg29; reg9=reg25*reg53;
    reg27=reg21*reg7; reg30=reg54*reg7; reg33=reg21*reg17; reg34=reg54*reg17; reg18=reg6+reg18;
    reg26=2*reg26; reg16=reg24+reg16; reg38=reg10*reg38; reg28=reg63*reg28; reg5=reg20*reg5;
    reg13=2*reg13; reg12=reg15+reg12; reg14=reg56*reg14; reg11=reg11+reg48; reg5=reg14+reg5;
    reg16=reg43*reg16; reg9=2*reg9; reg12=reg41*reg12; reg26=reg47*reg26; reg11=reg11/3;
    reg13=reg40*reg13; reg28=reg38+reg28; reg18=0.5*reg18; elem.epsilon[0][2]=reg18; reg6=reg29*reg53;
    reg34=reg27+reg34; reg30=reg33+reg30; reg13=reg28+reg13; reg34=reg54*reg34; reg35=reg35-reg11;
    reg12=reg16+reg12; reg26=reg5+reg26; reg9=reg25*reg9; reg45=reg45-reg11; reg5=reg18*reg2;
    reg6=2*reg6; reg30=reg21*reg30; reg26=reg19*reg26; reg5=reg0*reg5; reg35=pow(reg35,2);
    reg11=reg48-reg11; reg45=pow(reg45,2); reg30=reg34+reg30; reg9=reg12+reg9; reg6=reg29*reg6;
    reg13=reg42*reg13; reg45=reg35+reg45; reg11=pow(reg11,2); reg13=0.25*reg13; reg22=reg9*reg22;
    reg9=2*reg5; reg6=reg30+reg6; reg26=0.25*reg26; reg11=reg45+reg11; reg22=0.25*reg22;
    reg6=reg60*reg6; reg26=reg13+reg26; reg9=reg5*reg9; reg37=reg37-reg31; reg26=reg22+reg26;
    reg9=reg11+reg9; reg6=0.25*reg6; reg44=reg44-reg31; reg9=1.5*reg9; reg5=reg37*reg17;
    reg37=reg37*reg7; reg10=reg44*reg17; reg6=reg26+reg6; reg44=reg44*reg7; elem.ener=reg6/2;
    elem.sigma_von_mises=pow(reg9,0.5); elem.sigma[0][2]=reg18*reg53; elem.sigma[0][1]=reg37+reg10; elem.sigma[0][0]=reg5+reg44;
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
    T reg4=reg0*reg1; T reg5=reg3*reg1; T reg6=reg3*reg4; reg4=reg2*reg4; reg1=reg2*reg1;
    T reg7=reg2*reg4; T reg8=reg3*reg1; T reg9=reg3*reg6; reg4=reg3*reg4; T reg10=reg3*reg5;
    reg1=reg2*reg1; reg5=reg2*reg5; reg8=reg10+reg8; reg1=reg1-reg10; reg7=reg7-reg9;
    reg4=reg9+reg4; reg6=reg2*reg6; T reg11=reg10+reg5; reg6=reg9+reg6; reg8=reg3*reg8;
    reg1=reg2*reg1; reg9=reg2*reg7; T reg12=reg3*reg4; reg11=reg3*reg11; reg8=reg1-reg8;
    reg12=reg9-reg12; reg1=reg3*reg6; reg11=reg8-reg11; reg1=reg12-reg1; reg8=1-var_inter[1];
    reg9=1-var_inter[0]; reg7=reg7/reg1; reg12=reg8*elem.pos(0)[0]; T reg13=reg8*elem.pos(1)[0]; reg11=reg11/reg1;
    T reg14=reg9*elem.pos(0)[0]; T reg15=reg9*elem.pos(0)[1]; T reg16=reg8*elem.pos(0)[1]; reg4=reg4/reg1; T reg17=reg8*elem.pos(1)[1];
    T reg18=elem.pos(1)[0]*var_inter[0]; T reg19=elem.pos(1)[1]*var_inter[0]; T reg20=reg18+reg14; reg17=reg17-reg16; T reg21=elem.pos(2)[0]*var_inter[0];
    T reg22=reg7*reg11; T reg23=reg4*reg11; T reg24=elem.pos(2)[1]*var_inter[0]; T reg25=reg15+reg19; T reg26=var_inter[1]*elem.pos(2)[1];
    T reg27=var_inter[1]*elem.pos(2)[0]; reg13=reg13-reg12; T reg28=var_inter[1]*elem.pos(3)[1]; reg1=reg6/reg1; reg6=elem.pos(3)[0]*reg9;
    reg24=reg24-reg25; T reg29=var_inter[1]*elem.pos(3)[0]; reg27=reg13+reg27; reg13=reg7*reg22; reg21=reg21-reg20;
    T reg30=reg4*reg23; T reg31=elem.pos(3)[1]*reg9; T reg32=(*f.m).alpha*reg7; reg26=reg17+reg26; reg17=(*f.m).alpha*reg4;
    reg26=reg26-reg28; reg30=reg13-reg30; reg1=(*f.m).alpha*reg1; reg24=reg31+reg24; reg27=reg27-reg29;
    reg17=reg32+reg17; reg6=reg21+reg6; reg13=reg26*reg6; reg21=reg27*reg24; reg1=reg17+reg1;
    reg23=reg23/reg30; reg30=reg22/reg30; reg17=reg2*reg0; reg22=reg3*reg0; reg31=reg3*reg22;
    reg32=pow(reg2,2); reg13=reg21-reg13; reg21=pow(reg3,2); T reg33=reg2*reg17; reg23=reg23*reg1;
    reg1=reg30*reg1; reg27=reg27/reg13; reg6=reg6/reg13; reg26=reg26/reg13; reg24=reg24/reg13;
    reg23=reg1-reg23; reg1=1-(*f.m).resolution; reg31=reg33-reg31; reg21=reg32-reg21; reg30=var_inter[1]*reg24;
    reg17=reg17/reg31; reg22=reg22/reg31; reg31=reg21/reg31; reg21=var_inter[0]*reg27; reg32=var_inter[0]*reg26;
    reg33=var_inter[1]*reg6; T reg34=reg8*reg6; T reg35=reg9*reg27; T reg36=reg9*reg26; T reg37=reg8*reg24;
    reg23=reg1*reg23; T reg38=(*f.m).resolution*(*f.m).alpha; reg7=reg7*reg1; reg11=reg11*reg1; T reg39=reg34+reg21;
    reg1=reg4*reg1; reg4=(*f.m).resolution*reg17; T reg40=(*f.m).resolution*reg22; T reg41=(*f.m).resolution*reg31; T reg42=reg37+reg32;
    T reg43=reg35+reg33; reg23=reg38+reg23; reg38=reg36+reg30; T reg44=reg30-reg32; T reg45=0.5*reg38;
    T reg46=0.5*reg43; T reg47=reg21-reg33; T reg48=0.5*reg42; T reg49=0.5*reg39; T reg50=reg34-reg35;
    T reg51=reg36-reg37; reg23=(*f.m).deltaT*reg23; reg11=reg41+reg11; reg1=reg40+reg1; reg4=reg7+reg4;
    reg7=reg4*reg23; reg40=reg1*reg23; reg41=reg11*reg45; T reg52=reg11*reg46; T reg53=reg11*reg48;
    T reg54=reg11*reg49; T reg55=0.5*reg50; T reg56=0.5*reg51; T reg57=0.5*reg47; T reg58=0.5*reg44;
    T reg59=reg4*reg43; T reg60=reg1*reg38; T reg61=reg4*reg39; T reg62=reg1*reg42; T reg63=2*reg41;
    T reg64=reg8*var_inter[0]; T reg65=2*reg54; T reg66=reg1*reg43; reg52=2*reg52; T reg67=var_inter[1]*reg9;
    T reg68=reg4*reg38; T reg69=reg1*reg39; reg53=2*reg53; T reg70=reg11*reg58; T reg71=reg11*reg57;
    T reg72=reg7+reg40; T reg73=reg11*reg55; T reg74=reg11*reg56; T reg75=reg4*reg42; T reg76=reg62*reg39;
    reg71=2*reg71; T reg77=reg1*reg47; T reg78=reg52*reg49; T reg79=reg68*reg42; T reg80=reg61*reg43;
    reg70=2*reg70; T reg81=reg53*reg45; T reg82=reg64*elem.f_vol_e[1]; T reg83=reg67*elem.f_vol_e[0]; T reg84=reg4*reg51;
    T reg85=reg53*reg49; reg73=2*reg73; T reg86=reg1*reg50; T reg87=reg69*reg42; T reg88=reg1*reg51;
    T reg89=reg4*reg50; T reg90=reg65*reg48; reg74=2*reg74; T reg91=reg63*reg48; T reg92=reg1*reg44;
    T reg93=reg4*reg47; T reg94=reg59*reg39; T reg95=reg8*reg9; T reg96=var_inter[1]*var_inter[0]; T reg97=reg65*reg46;
    T reg98=reg38*reg75; T reg99=reg63*reg46; T reg100=reg66*reg38; T reg101=reg60*reg43; T reg102=reg52*reg45;
    T reg103=reg72*reg39; T reg104=reg72*reg38; T reg105=reg4*reg44; T reg106=reg57*reg74; T reg107=reg44*reg75;
    T reg108=reg65*reg57; T reg109=reg69*reg44; T reg110=reg53*reg57; T reg111=reg105*reg44; T reg112=reg71*reg57;
    T reg113=reg77*reg38; T reg114=reg62*reg43; T reg115=reg70*reg46; T reg116=reg103-reg82; T reg117=reg72*reg42;
    T reg118=reg72*reg50; T reg119=reg72*reg51; T reg120=reg63*reg45; T reg121=reg59*reg43; T reg122=reg52*reg46;
    reg102=reg101+reg102; T reg123=reg70*reg45; T reg124=reg93*reg43; T reg125=reg71*reg45; reg100=reg99+reg100;
    T reg126=reg71*reg55; T reg127=reg88*reg43; T reg128=reg77*reg42; T reg129=reg68*reg38; reg78=reg79+reg78;
    T reg130=reg66*reg42; T reg131=reg63*reg49; T reg132=reg48*reg73; T reg133=reg88*reg39; T reg134=reg48*reg74;
    T reg135=reg89*reg39; T reg136=reg65*reg45; T reg137=reg72*reg44; reg76=reg90+reg76; T reg138=reg53*reg48;
    T reg139=reg61*reg39; T reg140=reg71*reg48; T reg141=reg92*reg39; T reg142=reg70*reg48; T reg143=reg93*reg39;
    T reg144=reg52*reg48; T reg145=reg60*reg39; reg94=reg91+reg94; T reg146=reg44*reg84; T reg147=reg57*reg73;
    T reg148=reg44*reg86; T reg149=reg65*reg58; T reg150=reg62*reg47; T reg151=reg95*elem.f_vol_e[0]; T reg152=reg95*elem.f_vol_e[1];
    T reg153=reg58*reg74; T reg154=reg89*reg47; T reg155=reg64*elem.f_vol_e[0]; T reg156=reg58*reg73; T reg157=reg88*reg47;
    T reg158=reg63*reg57; T reg159=reg66*reg44; T reg160=reg96*elem.f_vol_e[0]; T reg161=reg96*elem.f_vol_e[1]; T reg162=reg52*reg57;
    T reg163=reg68*reg44; T reg164=reg67*elem.f_vol_e[1]; T reg165=reg51*reg84; T reg166=reg70*reg57; T reg167=reg55*reg73;
    T reg168=reg77*reg44; T reg169=reg51*reg86; T reg170=reg55*reg74; T reg171=reg105*reg51; T reg172=reg51*reg75;
    T reg173=reg45*reg74; T reg174=reg89*reg43; T reg175=reg45*reg73; T reg176=reg105*reg38; T reg177=reg71*reg46;
    T reg178=reg69*reg38; T reg179=reg53*reg46; reg98=reg97+reg98; T reg180=reg38*reg86; T reg181=reg46*reg74;
    T reg182=reg38*reg84; T reg183=reg46*reg73; T reg184=reg63*reg58; T reg185=reg59*reg47; T reg186=reg52*reg58;
    T reg187=reg60*reg47; T reg188=reg70*reg58; T reg189=reg93*reg47; T reg190=reg71*reg58; T reg191=reg92*reg47;
    T reg192=reg53*reg58; T reg193=reg61*reg47; reg81=reg80+reg81; T reg194=reg92*reg43; reg85=reg87+reg85;
    T reg195=reg72*reg43; reg88=reg88*reg50; T reg196=reg70*reg56; T reg197=reg72*reg47; T reg198=reg56*reg73;
    reg93=reg93*reg50; reg89=reg89*reg50; reg84=reg42*reg84; T reg199=reg56*reg74; T reg200=reg65*reg49;
    reg73=reg49*reg73; reg75=reg42*reg75; reg62=reg62*reg50; T reg201=reg65*reg56; T reg202=reg71*reg56;
    reg92=reg92*reg50; T reg203=reg70*reg49; T reg204=reg104-reg83; reg74=reg49*reg74; reg86=reg42*reg86;
    T reg205=reg61*reg50; T reg206=reg53*reg56; reg71=reg71*reg49; T reg207=reg52*reg56; reg59=reg59*reg50;
    T reg208=reg68*reg51; reg52=reg52*reg55; reg105=reg105*reg42; reg70=reg70*reg55; T reg209=reg63*reg56;
    T reg210=reg65*reg55; T reg211=reg69*reg51; T reg212=reg60*reg50; reg66=reg66*reg51; reg53=reg53*reg55;
    T reg213=reg63*reg55; reg77=reg77*reg51; reg113=reg115-reg113; reg207=reg207-reg212; reg173=reg174-reg173;
    reg116=reg116*reg13; reg175=reg127-reg175; reg176=reg177-reg176; reg126=reg171+reg126; reg202=reg92+reg202;
    reg92=reg100*reg13; reg125=reg194-reg125; reg53=reg53-reg211; reg115=reg117+reg155; reg123=reg124-reg123;
    reg124=reg102*reg13; reg127=reg118+reg152; reg171=reg119+reg151; reg196=reg93+reg196; reg121=reg121+reg120;
    reg170=reg169+reg170; reg70=reg77+reg70; reg167=reg165+reg167; reg166=reg168+reg166; reg172=reg172-reg210;
    reg162=reg162-reg163; reg52=reg52-reg208; reg159=reg159-reg158; reg77=reg195+reg164; reg156=reg157+reg156;
    reg66=reg66-reg213; reg153=reg154+reg153; reg150=reg150-reg149; reg93=reg13*reg81; reg198=reg88+reg198;
    reg192=reg192-reg193; reg190=reg191+reg190; reg199=reg89+reg199; reg188=reg189+reg188; reg186=reg186-reg187;
    reg204=reg204*reg13; reg185=reg185-reg184; reg182=reg183-reg182; reg62=reg62-reg201; reg180=reg181-reg180;
    reg88=reg98*reg13; reg179=reg179+reg178; reg206=reg206-reg205; reg89=reg85*reg13; reg154=reg137+reg160;
    reg74=reg86-reg74; reg147=reg146+reg147; reg133=reg132-reg133; reg73=reg84-reg73; reg106=reg148+reg106;
    reg84=reg78*reg13; reg141=reg140-reg141; reg143=reg142-reg143; reg138=reg138+reg139; reg107=reg107-reg108;
    reg75=reg75+reg200; reg86=reg94*reg13; reg203=reg128-reg203; reg144=reg144+reg145; reg122=reg122+reg129;
    reg135=reg134-reg135; reg114=reg136+reg114; reg128=reg76*reg13; reg130=reg130+reg131; reg59=reg59-reg209;
    reg112=reg111+reg112; reg111=reg197+reg161; reg71=reg105-reg71; reg110=reg110-reg109; reg176=reg176*reg13;
    reg116=ponderation*reg116; reg52=reg52*reg13; reg62=reg62*reg13; reg130=reg130*reg13; reg180=reg180*reg13;
    reg105=reg77*reg13; reg132=ponderation*reg89; reg53=reg53*reg13; reg143=reg143*reg13; reg134=ponderation*reg88;
    reg179=reg179*reg13; reg144=reg144*reg13; reg159=reg159*reg13; reg140=ponderation*reg93; reg75=reg75*reg13;
    reg150=reg150*reg13; reg198=reg198*reg13; reg114=reg114*reg13; reg135=reg135*reg13; reg192=reg192*reg13;
    reg142=ponderation*reg128; reg153=reg153*reg13; reg190=reg190*reg13; reg138=reg138*reg13; reg188=reg188*reg13;
    reg199=reg199*reg13; reg66=reg66*reg13; reg186=reg186*reg13; reg133=reg133*reg13; reg204=ponderation*reg204;
    reg185=reg185*reg13; reg156=reg156*reg13; reg141=reg141*reg13; reg74=reg74*reg13; reg182=reg182*reg13;
    reg147=reg147*reg13; reg202=reg202*reg13; reg125=reg125*reg13; reg146=reg127*reg13; reg148=reg154*reg13;
    reg123=reg123*reg13; reg167=reg13*reg167; reg126=reg126*reg13; reg170=reg13*reg170; reg157=reg111*reg13;
    reg106=reg106*reg13; reg70=reg70*reg13; reg165=ponderation*reg124; reg168=reg171*reg13; reg107=reg107*reg13;
    reg203=reg203*reg13; reg59=reg59*reg13; reg196=reg196*reg13; reg121=reg121*reg13; reg113=reg113*reg13;
    reg162=reg162*reg13; reg112=reg112*reg13; reg206=reg206*reg13; reg207=reg207*reg13; reg175=reg175*reg13;
    reg71=reg71*reg13; reg110=reg110*reg13; reg169=ponderation*reg84; reg174=ponderation*reg86; reg173=reg173*reg13;
    reg177=reg115*reg13; reg172=reg172*reg13; reg122=reg122*reg13; reg166=reg166*reg13; reg73=reg73*reg13;
    reg181=ponderation*reg92; T tmp_3_1=ponderation*reg135; T tmp_5_1=ponderation*reg153; T tmp_0_1=ponderation*reg170; T tmp_2_4=ponderation*reg71;
    T tmp_0_7=ponderation*reg66; T tmp_4_4=ponderation*reg112; reg66=ponderation*reg148; sollicitation[indices[2]+0]+=reg66; T tmp_2_5=ponderation*reg203;
    T tmp_3_0=ponderation*reg133; T tmp_0_0=ponderation*reg167; T tmp_5_0=ponderation*reg156; T tmp_0_5=ponderation*reg70; T tmp_0_2=ponderation*reg172;
    T tmp_4_5=ponderation*reg166; T tmp_4_7=ponderation*reg159; T tmp_2_7=ponderation*reg130; T tmp_2_6=-reg169; T tmp_2_3=-reg132;
    T tmp_0_6=ponderation*reg52; T tmp_4_6=ponderation*reg162; reg52=ponderation*reg105; sollicitation[indices[3]+1]+=reg52; T tmp_6_5=ponderation*reg113;
    T tmp_2_0=ponderation*reg73; T tmp_7_0=ponderation*reg175; T tmp_1_3=ponderation*reg206; T tmp_7_1=ponderation*reg173; T tmp_3_7=-reg174;
    T tmp_6_6=ponderation*reg122; T tmp_6_7=-reg181; T tmp_4_0=ponderation*reg147; T tmp_7_4=ponderation*reg125; T tmp_1_4=ponderation*reg202;
    reg70=ponderation*reg157; sollicitation[indices[2]+1]+=reg70; T tmp_7_5=ponderation*reg123; T tmp_0_4=ponderation*reg126; T tmp_4_1=ponderation*reg106;
    T tmp_1_7=ponderation*reg59; T tmp_7_6=-reg165; T tmp_7_7=ponderation*reg121; T tmp_1_5=ponderation*reg196; reg59=ponderation*reg168;
    sollicitation[indices[0]+0]+=reg59; T tmp_4_2=ponderation*reg107; reg71=ponderation*reg146; sollicitation[indices[0]+1]+=reg71; reg73=ponderation*reg177;
    sollicitation[indices[1]+0]+=reg73; T tmp_4_3=ponderation*reg110; T tmp_1_6=ponderation*reg207; T tmp_7_3=-reg140; T tmp_2_2=ponderation*reg75;
    T tmp_5_2=ponderation*reg150; T tmp_7_2=ponderation*reg114; T tmp_1_0=ponderation*reg198; T tmp_5_3=ponderation*reg192; sollicitation[indices[3]+0]+=-reg204;
    T tmp_3_2=-reg142; T tmp_5_4=ponderation*reg190; T tmp_5_5=ponderation*reg188; T tmp_1_1=ponderation*reg199; T tmp_3_3=ponderation*reg138;
    T tmp_5_6=ponderation*reg186; T tmp_2_1=ponderation*reg74; T tmp_5_7=ponderation*reg185; T tmp_3_4=ponderation*reg141; T tmp_6_0=ponderation*reg182;
    T tmp_6_1=ponderation*reg180; T tmp_0_3=ponderation*reg53; T tmp_1_2=ponderation*reg62; T tmp_6_2=-reg134; T tmp_3_5=ponderation*reg143;
    T tmp_6_3=ponderation*reg179; T tmp_3_6=ponderation*reg144; sollicitation[indices[1]+1]+=-reg116; T tmp_6_4=ponderation*reg176;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg2*reg1; T reg6=reg3*reg4; reg4=reg2*reg4; reg1=reg3*reg1;
    T reg7=reg3*reg1; T reg8=reg2*reg4; T reg9=reg2*reg5; T reg10=reg3*reg6; reg1=reg2*reg1;
    reg6=reg2*reg6; reg5=reg3*reg5; reg7=reg7-reg9; reg1=reg9+reg1; reg4=reg3*reg4;
    reg10=reg10-reg8; reg6=reg8+reg6; T reg11=reg2*reg6; reg7=reg3*reg7; reg1=reg2*reg1;
    reg4=reg8+reg4; reg8=reg3*reg10; T reg12=reg9+reg5; T reg13=reg2*reg4; reg11=reg8-reg11;
    reg1=reg7-reg1; reg12=reg2*reg12; reg13=reg11-reg13; reg7=1-var_inter[0]; reg12=reg1-reg12;
    reg1=1-var_inter[1]; reg10=reg10/reg13; reg8=elem.pos(1)[0]*var_inter[0]; reg12=reg12/reg13; reg11=reg7*elem.pos(0)[0];
    reg6=reg6/reg13; T reg14=reg1*elem.pos(0)[1]; T reg15=reg1*elem.pos(1)[1]; T reg16=reg7*elem.pos(0)[1]; T reg17=reg1*elem.pos(1)[0];
    T reg18=reg1*elem.pos(0)[0]; T reg19=elem.pos(1)[1]*var_inter[0]; T reg20=reg6*reg12; T reg21=reg10*reg12; T reg22=var_inter[1]*elem.pos(2)[0];
    reg17=reg17-reg18; T reg23=var_inter[1]*elem.pos(2)[1]; T reg24=reg16+reg19; reg15=reg15-reg14; T reg25=elem.pos(2)[0]*var_inter[0];
    T reg26=elem.pos(2)[1]*var_inter[0]; T reg27=reg8+reg11; T reg28=var_inter[1]*elem.pos(3)[0]; T reg29=(*f.m).alpha*reg6; reg22=reg17+reg22;
    reg23=reg15+reg23; reg15=(*f.m).alpha*reg10; reg17=elem.pos(3)[1]*reg7; T reg30=reg6*reg20; reg26=reg26-reg24;
    T reg31=var_inter[1]*elem.pos(3)[1]; T reg32=reg10*reg21; reg13=reg4/reg13; reg25=reg25-reg27; reg4=elem.pos(3)[0]*reg7;
    reg30=reg32-reg30; reg29=reg15+reg29; reg26=reg17+reg26; reg22=reg22-reg28; reg13=(*f.m).alpha*reg13;
    reg23=reg23-reg31; reg4=reg25+reg4; reg15=reg2*reg0; reg17=reg3*reg0; reg25=reg22*reg26;
    reg32=reg23*reg4; reg21=reg21/reg30; reg13=reg29+reg13; reg30=reg20/reg30; reg30=reg30*reg13;
    reg13=reg21*reg13; reg20=pow(reg2,2); reg21=pow(reg3,2); reg29=reg2*reg15; reg32=reg25-reg32;
    reg25=reg3*reg17; reg22=reg22/reg32; reg26=reg26/reg32; reg4=reg4/reg32; reg23=reg23/reg32;
    reg29=reg25-reg29; reg20=reg21-reg20; reg21=1-(*f.m).resolution; reg30=reg13-reg30; reg13=var_inter[1]*reg4;
    reg25=var_inter[1]*reg26; T reg33=var_inter[0]*reg23; T reg34=reg7*reg22; T reg35=reg7*reg23; T reg36=reg1*reg26;
    reg30=reg21*reg30; T reg37=(*f.m).resolution*(*f.m).alpha; reg17=reg17/reg29; reg15=reg15/reg29; reg29=reg20/reg29;
    reg20=reg34+reg13; T reg38=reg35+reg25; T reg39=var_inter[0]*reg22; T reg40=reg36+reg33; T reg41=reg1*reg4;
    T reg42=(*f.m).resolution*reg29; T reg43=(*f.m).resolution*reg17; T reg44=(*f.m).resolution*reg15; reg6=reg6*reg21; reg30=reg37+reg30;
    reg12=reg12*reg21; reg21=reg10*reg21; reg10=0.5*reg38; reg37=0.5*reg20; reg43=reg21+reg43;
    reg21=reg39-reg13; reg6=reg44+reg6; reg44=reg25-reg33; reg12=reg42+reg12; reg42=0.5*reg40;
    T reg45=reg41+reg39; reg30=(*f.m).deltaT*reg30; T reg46=reg41-reg34; T reg47=reg35-reg36; T reg48=reg12*reg10;
    T reg49=reg12*reg37; T reg50=reg12*reg42; T reg51=reg6*reg30; T reg52=reg43*reg30; T reg53=0.5*reg46;
    T reg54=0.5*reg47; T reg55=0.5*reg45; T reg56=0.5*reg21; T reg57=0.5*reg44; T reg58=reg43*reg20;
    T reg59=reg12*reg55; T reg60=var_inter[1]*reg7; T reg61=reg1*var_inter[0]; T reg62=reg12*reg54; T reg63=reg12*reg53;
    T reg64=reg52+reg51; T reg65=reg6*reg45; T reg66=2*reg48; T reg67=reg6*reg20; reg49=2*reg49;
    T reg68=reg43*reg38; T reg69=reg12*reg57; T reg70=reg12*reg56; reg50=2*reg50; T reg71=reg43*reg47;
    T reg72=reg60*elem.f_vol_e[0]; reg63=2*reg63; T reg73=reg6*reg46; T reg74=reg61*elem.f_vol_e[1]; reg62=2*reg62;
    T reg75=reg43*reg40; T reg76=reg58*reg45; T reg77=reg66*reg42; T reg78=reg49*reg55; T reg79=reg68*reg40;
    T reg80=reg50*reg55; T reg81=reg65*reg40; T reg82=reg1*reg7; T reg83=var_inter[1]*var_inter[0]; T reg84=reg66*reg37;
    T reg85=reg67*reg38; T reg86=reg64*reg45; T reg87=reg64*reg38; T reg88=2*reg59; T reg89=reg43*reg44;
    reg70=2*reg70; T reg90=reg6*reg21; reg69=2*reg69; T reg91=reg43*reg46; T reg92=reg6*reg40;
    T reg93=reg43*reg45; T reg94=reg6*reg44; T reg95=reg43*reg21; T reg96=reg6*reg38; T reg97=reg49*reg53;
    T reg98=reg82*elem.f_vol_e[0]; T reg99=reg67*reg47; T reg100=reg66*reg53; T reg101=reg82*elem.f_vol_e[1]; T reg102=reg49*reg57;
    T reg103=reg96*reg21; T reg104=reg91*reg46; T reg105=reg54*reg62; T reg106=reg40*reg75; T reg107=reg92*reg46;
    T reg108=reg88*reg54; T reg109=reg61*elem.f_vol_e[0]; T reg110=reg69*reg57; T reg111=reg88*reg55; T reg112=reg93*reg46;
    T reg113=reg50*reg54; T reg114=reg95*reg21; T reg115=reg94*reg46; T reg116=reg70*reg54; T reg117=reg66*reg56;
    T reg118=reg70*reg56; T reg119=reg89*reg44; T reg120=reg67*reg44; T reg121=reg83*elem.f_vol_e[0]; T reg122=reg64*reg40;
    T reg123=reg86-reg74; T reg124=reg64*reg46; T reg125=reg64*reg44; T reg126=reg64*reg47; T reg127=reg64*reg21;
    T reg128=reg66*reg10; T reg129=reg87-reg72; T reg130=reg58*reg20; T reg131=reg49*reg37; T reg132=reg64*reg20;
    reg85=reg84+reg85; T reg133=reg88*reg53; T reg134=reg65*reg47; T reg135=reg50*reg53; T reg136=reg49*reg54;
    T reg137=reg58*reg46; T reg138=reg89*reg47; T reg139=reg66*reg54; T reg140=reg70*reg53; T reg141=reg69*reg55;
    T reg142=reg90*reg47; T reg143=reg66*reg57; T reg144=reg69*reg53; reg58=reg58*reg21; T reg145=reg68*reg47;
    T reg146=reg53*reg63; T reg147=reg90*reg40; T reg148=reg50*reg42; T reg149=reg93*reg45; T reg150=reg68*reg38;
    T reg151=reg69*reg56; T reg152=reg53*reg62; T reg153=reg70*reg42; T reg154=reg94*reg45; T reg155=reg66*reg55;
    T reg156=reg70*reg55; T reg157=reg47*reg71; T reg158=reg89*reg40; T reg159=reg69*reg42; T reg160=reg95*reg45;
    reg67=reg67*reg40; reg78=reg79+reg78; reg76=reg77+reg76; T reg161=reg83*elem.f_vol_e[1]; T reg162=reg47*reg75;
    T reg163=reg95*reg46; T reg164=reg69*reg54; T reg165=reg49*reg56; T reg166=reg47*reg73; reg80=reg81+reg80;
    T reg167=reg96*reg45; T reg168=reg90*reg44; T reg169=reg49*reg42; T reg170=reg96*reg46; T reg171=reg68*reg44;
    T reg172=reg60*elem.f_vol_e[1]; reg131=reg131+reg150; reg136=reg136-reg170; reg135=reg135-reg134; reg162=reg162-reg133;
    reg141=reg147-reg141; reg147=reg85*reg32; T reg173=reg132+reg172; reg152=reg166+reg152; reg156=reg158-reg156;
    reg129=reg129*reg32; reg130=reg130+reg128; reg158=reg127+reg161; reg166=reg126+reg98; T reg174=reg125+reg121;
    T reg175=reg80*reg32; reg123=reg123*reg32; reg164=reg163+reg164; reg163=reg124+reg101; T reg176=reg122+reg109;
    T reg177=reg76*reg32; reg165=reg165-reg171; reg118=reg119+reg118; reg120=reg120-reg117; reg116=reg115+reg116;
    reg169=reg169+reg167; reg106=reg106+reg111; reg113=reg113-reg112; reg160=reg159-reg160; reg107=reg107-reg108;
    reg154=reg153-reg154; reg110=reg114+reg110; reg105=reg104+reg105; reg151=reg168+reg151; reg148=reg148+reg149;
    reg140=reg138+reg140; reg58=reg58-reg143; reg137=reg137-reg139; reg104=reg78*reg32; reg144=reg142+reg144;
    reg99=reg99-reg100; reg102=reg102-reg103; reg97=reg97-reg145; reg146=reg157+reg146; reg67=reg67+reg155;
    reg102=reg102*reg32; reg152=reg32*reg152; reg130=reg130*reg32; reg120=reg120*reg32; reg146=reg32*reg146;
    reg110=reg110*reg32; reg114=reg166*reg32; reg131=reg131*reg32; reg115=ponderation*reg147; reg151=reg151*reg32;
    reg58=reg58*reg32; reg119=reg163*reg32; reg165=reg165*reg32; reg140=reg140*reg32; reg137=reg137*reg32;
    reg144=reg144*reg32; reg97=reg97*reg32; reg99=reg99*reg32; reg105=reg105*reg32; reg107=reg107*reg32;
    reg113=reg113*reg32; reg116=reg116*reg32; reg118=reg118*reg32; reg138=ponderation*reg177; reg106=reg106*reg32;
    reg169=reg169*reg32; reg160=reg160*reg32; reg154=reg154*reg32; reg148=reg148*reg32; reg67=reg67*reg32;
    reg142=ponderation*reg104; reg141=reg141*reg32; reg156=reg156*reg32; reg153=ponderation*reg175; reg164=reg164*reg32;
    reg157=reg173*reg32; reg159=reg174*reg32; reg123=ponderation*reg123; reg162=reg162*reg32; reg129=ponderation*reg129;
    reg136=reg136*reg32; reg135=reg135*reg32; reg168=reg176*reg32; T reg178=reg158*reg32; T tmp_4_5=ponderation*reg151;
    reg151=ponderation*reg159; sollicitation[indices[2]+0]+=reg151; T tmp_3_5=ponderation*reg160; T tmp_3_4=ponderation*reg154; T tmp_7_7=ponderation*reg130;
    T tmp_3_6=ponderation*reg169; T tmp_3_3=ponderation*reg148; sollicitation[indices[1]+1]+=-reg123; T tmp_2_7=ponderation*reg67; reg67=ponderation*reg114;
    sollicitation[indices[0]+0]+=reg67; T tmp_2_6=-reg142; T tmp_0_0=ponderation*reg146; T tmp_2_5=ponderation*reg141; T tmp_4_4=ponderation*reg118;
    T tmp_2_4=ponderation*reg156; reg118=ponderation*reg168; sollicitation[indices[1]+0]+=reg118; T tmp_2_3=-reg153; T tmp_0_1=ponderation*reg152;
    reg123=ponderation*reg119; sollicitation[indices[0]+1]+=reg123; T tmp_2_2=ponderation*reg106; T tmp_0_3=ponderation*reg135; reg106=ponderation*reg157;
    sollicitation[indices[3]+1]+=reg106; T tmp_5_7=ponderation*reg58; T tmp_6_6=ponderation*reg131; T tmp_0_4=ponderation*reg140; T tmp_0_2=ponderation*reg162;
    T tmp_5_6=ponderation*reg102; T tmp_0_5=ponderation*reg144; T tmp_6_7=-reg115; T tmp_0_6=ponderation*reg97; T tmp_0_7=ponderation*reg99;
    T tmp_1_7=ponderation*reg137; T tmp_5_5=ponderation*reg110; T tmp_1_5=ponderation*reg164; T tmp_1_1=ponderation*reg105; T tmp_1_2=ponderation*reg107;
    sollicitation[indices[3]+0]+=-reg129; T tmp_4_7=ponderation*reg120; T tmp_1_3=ponderation*reg113; T tmp_1_6=ponderation*reg136; T tmp_3_7=-reg138;
    reg58=ponderation*reg178; sollicitation[indices[2]+1]+=reg58; T tmp_4_6=ponderation*reg165; T tmp_1_4=ponderation*reg116;
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
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=1.0/(*f.m).elastic_modulus; T reg6=reg0*reg1; T reg7=elem.pos(1)[1]*var_inter[0]; T reg8=reg2*elem.pos(0)[1];
    T reg9=reg2*elem.pos(0)[0]; T reg10=elem.pos(1)[0]*var_inter[0]; T reg11=reg3*elem.pos(1)[1]; T reg12=reg3*elem.pos(0)[1]; T reg13=reg3*elem.pos(1)[0];
    T reg14=reg3*elem.pos(0)[0]; T reg15=reg5*reg1; reg1=reg4*reg1; T reg16=elem.pos(2)[1]*var_inter[0]; T reg17=reg8+reg7;
    reg13=reg13-reg14; T reg18=var_inter[1]*elem.pos(2)[0]; T reg19=reg5*reg6; reg6=reg4*reg6; reg11=reg11-reg12;
    T reg20=var_inter[1]*elem.pos(2)[1]; T reg21=elem.pos(2)[0]*var_inter[0]; T reg22=reg10+reg9; T reg23=reg4*reg15; T reg24=reg4*reg1;
    reg15=reg5*reg15; T reg25=reg4*reg19; reg20=reg11+reg20; reg11=var_inter[1]*elem.pos(3)[1]; reg16=reg16-reg17;
    reg21=reg21-reg22; T reg26=elem.pos(3)[0]*reg2; T reg27=elem.pos(3)[1]*reg2; T reg28=var_inter[1]*elem.pos(3)[0]; reg18=reg13+reg18;
    reg13=reg4*reg6; reg19=reg5*reg19; reg26=reg21+reg26; reg20=reg20-reg11; reg16=reg27+reg16;
    reg19=reg19-reg13; reg23=reg24+reg23; reg18=reg18-reg28; reg25=reg13+reg25; reg6=reg5*reg6;
    reg15=reg15-reg24; reg1=reg5*reg1; reg21=reg5*reg19; reg27=reg4*reg25; T reg29=reg18*reg16;
    T reg30=reg20*reg26; reg6=reg13+reg6; reg15=reg5*reg15; reg13=reg24+reg1; T reg31=reg4*reg0;
    reg23=reg4*reg23; T reg32=reg5*reg0; T reg33=reg4*reg31; T reg34=pow(reg5,2); T reg35=pow(reg4,2);
    reg23=reg15-reg23; reg27=reg21-reg27; reg15=reg4*reg6; reg21=reg5*reg32; reg13=reg4*reg13;
    reg30=reg29-reg30; reg13=reg23-reg13; reg16=reg16/reg30; reg18=reg18/reg30; reg33=reg21-reg33;
    reg15=reg27-reg15; reg35=reg34-reg35; reg26=reg26/reg30; reg20=reg20/reg30; reg21=reg2*reg18;
    reg23=reg3*reg26; reg27=reg2*reg20; reg29=reg3*reg16; reg35=reg35/reg33; reg34=var_inter[0]*reg20;
    T reg36=var_inter[0]*reg18; T reg37=1-(*f.m).resolution; reg13=reg13/reg15; T reg38=var_inter[1]*reg16; T reg39=var_inter[1]*reg26;
    T reg40=reg29+reg34; T reg41=reg21+reg39; T reg42=reg23+reg36; T reg43=reg27+reg38; reg32=reg32/reg33;
    reg33=reg31/reg33; reg19=reg19/reg15; reg25=reg25/reg15; reg31=(*f.m).resolution*reg35; T reg44=reg13*reg37;
    T reg45=0.5*reg43; T reg46=0.5*reg41; T reg47=reg36-reg39; T reg48=reg38-reg34; T reg49=(*f.m).resolution*reg33;
    T reg50=reg25*reg37; T reg51=0.5*reg40; T reg52=0.5*reg42; T reg53=reg19*reg37; T reg54=(*f.m).resolution*reg32;
    T reg55=reg27-reg29; T reg56=reg23-reg21; reg44=reg31+reg44; reg31=reg44*reg51; T reg57=reg44*reg46;
    T reg58=reg44*reg45; T reg59=reg44*reg52; T reg60=0.5*reg56; T reg61=0.5*reg55; reg50=reg49+reg50;
    reg49=0.5*reg48; T reg62=0.5*reg47; reg54=reg53+reg54; reg53=reg54*reg42; T reg63=reg50*reg41;
    T reg64=reg44*reg60; T reg65=2*reg59; T reg66=reg44*reg61; T reg67=2*reg58; T reg68=reg50*reg43;
    T reg69=reg50*reg42; T reg70=reg44*reg49; T reg71=reg50*reg40; T reg72=reg44*reg62; T reg73=reg54*reg40;
    T reg74=reg54*reg41; reg57=2*reg57; T reg75=reg54*reg43; reg31=2*reg31; T reg76=reg43*reg73;
    T reg77=reg50*reg48; reg70=2*reg70; T reg78=reg31*reg52; T reg79=reg69*reg40; T reg80=reg50*reg47;
    T reg81=reg65*reg46; T reg82=reg71*reg42; T reg83=reg63*reg43; T reg84=reg53*reg41; T reg85=reg67*reg46;
    T reg86=reg50*reg55; T reg87=reg65*reg51; T reg88=reg31*reg45; T reg89=reg54*reg56; T reg90=reg57*reg52;
    T reg91=reg75*reg40; reg66=2*reg66; T reg92=reg54*reg55; T reg93=reg54*reg47; reg72=2*reg72;
    T reg94=reg68*reg41; T reg95=reg54*reg48; T reg96=reg57*reg45; reg64=2*reg64; T reg97=reg50*reg56;
    T reg98=reg74*reg42; T reg99=reg67*reg51; T reg100=reg67*reg49; T reg101=reg75*reg43; T reg102=reg52*reg64;
    T reg103=reg72*reg46; T reg104=reg40*reg92; T reg105=reg74*reg47; T reg106=reg95*reg43; T reg107=reg45*reg64;
    T reg108=reg57*reg49; T reg109=reg46*reg64; T reg110=reg80*reg40; T reg111=reg40*reg97; T reg112=reg69*reg43;
    T reg113=reg31*reg46; T reg114=reg72*reg52; T reg115=reg52*reg66; T reg116=reg95*reg40; reg76=reg81+reg76;
    T reg117=reg43*reg92; T reg118=reg46*reg66; T reg119=reg40*reg73; reg78=reg79+reg78; T reg120=reg65*reg52;
    T reg121=reg43*reg97; T reg122=reg53*reg42; T reg123=reg72*reg62; T reg124=reg95*reg48; T reg125=reg72*reg51;
    T reg126=reg77*reg42; T reg127=reg31*reg62; T reg128=reg69*reg48; T reg129=reg70*reg51; T reg130=reg93*reg42;
    T reg131=reg65*reg62; T reg132=reg48*reg73; T reg133=reg57*reg51; T reg134=reg68*reg42; T reg135=reg62*reg66;
    T reg136=reg48*reg97; reg98=reg99+reg98; T reg137=reg62*reg64; T reg138=reg48*reg92; reg90=reg91+reg90;
    T reg139=reg89*reg41; T reg140=reg45*reg66; T reg141=reg63*reg40; T reg142=reg67*reg52; T reg143=reg86*reg41;
    T reg144=reg51*reg64; T reg145=reg86*reg42; reg83=reg85+reg83; T reg146=reg51*reg66; T reg147=reg89*reg42;
    T reg148=reg72*reg45; T reg149=reg93*reg41; T reg150=reg65*reg45; T reg151=reg70*reg46; reg82=reg87+reg82;
    T reg152=reg71*reg41; T reg153=reg80*reg43; T reg154=reg31*reg51; reg88=reg84+reg88; reg97=reg55*reg97;
    T reg155=reg61*reg66; T reg156=reg89*reg56; T reg157=reg77*reg41; T reg158=reg49*reg66; reg89=reg89*reg47;
    T reg159=reg61*reg64; T reg160=reg86*reg56; T reg161=reg69*reg55; T reg162=reg60*reg64; reg64=reg49*reg64;
    reg86=reg86*reg47; T reg163=reg67*reg60; T reg164=reg70*reg45; T reg165=reg63*reg55; T reg166=reg31*reg60;
    reg92=reg55*reg92; T reg167=reg67*reg62; reg63=reg63*reg48; T reg168=reg57*reg60; reg95=reg95*reg55;
    T reg169=reg75*reg55; T reg170=reg57*reg62; T reg171=reg75*reg48; T reg172=reg72*reg60; T reg173=reg80*reg48;
    T reg174=reg70*reg60; reg80=reg80*reg55; T reg175=reg70*reg52; T reg176=reg70*reg62; T reg177=reg67*reg61;
    T reg178=reg74*reg56; T reg179=reg68*reg47; reg96=reg94+reg96; T reg180=reg70*reg49; T reg181=reg57*reg61;
    T reg182=reg68*reg56; T reg183=reg93*reg47; reg73=reg55*reg73; reg57=reg57*reg46; T reg184=reg72*reg49;
    reg70=reg70*reg61; reg93=reg93*reg56; T reg185=reg77*reg47; reg74=reg74*reg41; T reg186=reg71*reg56;
    T reg187=reg65*reg61; T reg188=reg65*reg60; reg66=reg60*reg66; reg71=reg71*reg47; T reg189=reg53*reg56;
    T reg190=reg31*reg61; T reg191=reg65*reg49; T reg192=reg67*reg45; T reg193=reg53*reg47; reg77=reg77*reg56;
    reg72=reg72*reg61; reg31=reg31*reg49; reg127=reg127-reg128; reg135=reg136+reg135; reg132=reg132-reg131;
    reg137=reg138+reg137; reg162=reg92+reg162; reg66=reg97+reg66; reg108=reg108-reg179; reg180=reg183+reg180;
    reg105=reg105-reg100; reg184=reg185+reg184; reg117=reg109-reg117; reg121=reg118-reg121; reg31=reg31-reg193;
    reg71=reg71-reg191; reg92=reg76*reg30; reg113=reg113+reg112; reg158=reg89+reg158; reg106=reg103-reg106;
    reg89=reg30*reg88; reg107=reg143-reg107; reg140=reg139-reg140; reg64=reg86+reg64; reg86=reg83*reg30;
    reg63=reg63-reg167; reg148=reg157-reg148; reg170=reg170-reg171; reg176=reg173+reg176; reg153=reg151-reg153;
    reg123=reg124+reg123; reg119=reg119+reg120; reg154=reg154+reg122; reg70=reg93+reg70; reg73=reg73-reg188;
    reg93=reg82*reg30; reg152=reg150+reg152; reg181=reg181-reg182; reg97=reg90*reg30; reg147=reg146-reg147;
    reg115=reg111-reg115; reg165=reg165-reg163; reg145=reg144-reg145; reg178=reg178-reg177; reg172=reg95+reg172;
    reg168=reg168-reg169; reg166=reg166-reg161; reg102=reg104-reg102; reg141=reg141+reg142; reg155=reg156+reg155;
    reg95=reg98*reg30; reg186=reg186-reg187; reg103=reg96*reg30; reg114=reg116-reg114; reg133=reg133+reg134;
    reg104=reg78*reg30; reg164=reg149-reg164; reg190=reg190-reg189; reg174=reg80+reg174; reg130=reg129-reg130;
    reg57=reg57+reg101; reg175=reg110-reg175; reg74=reg74+reg192; reg159=reg160+reg159; reg126=reg125-reg126;
    reg72=reg77+reg72; reg117=reg117*reg30; reg115=reg115*reg30; reg170=reg170*reg30; reg77=ponderation*reg104;
    reg119=reg119*reg30; reg121=reg121*reg30; reg168=reg168*reg30; reg174=reg174*reg30; reg155=reg155*reg30;
    reg80=ponderation*reg89; reg152=reg152*reg30; reg186=reg186*reg30; reg158=reg158*reg30; reg190=reg190*reg30;
    reg71=reg71*reg30; reg159=reg159*reg30; reg72=reg72*reg30; reg31=reg31*reg30; reg64=reg64*reg30;
    reg70=reg70*reg30; reg184=reg184*reg30; reg165=reg165*reg30; reg181=reg181*reg30; reg180=reg180*reg30;
    reg178=reg178*reg30; reg108=reg108*reg30; reg63=reg63*reg30; reg102=reg102*reg30; reg105=reg105*reg30;
    reg164=reg164*reg30; reg109=ponderation*reg86; reg145=reg145*reg30; reg162=reg30*reg162; reg147=reg147*reg30;
    reg148=reg148*reg30; reg73=reg73*reg30; reg110=ponderation*reg93; reg123=reg123*reg30; reg154=reg154*reg30;
    reg74=reg74*reg30; reg127=reg127*reg30; reg126=reg126*reg30; reg132=reg132*reg30; reg130=reg130*reg30;
    reg66=reg30*reg66; reg135=reg135*reg30; reg133=reg133*reg30; reg111=ponderation*reg103; reg137=reg137*reg30;
    reg116=ponderation*reg95; reg106=reg106*reg30; reg175=reg175*reg30; reg153=reg153*reg30; reg113=reg113*reg30;
    reg172=reg172*reg30; reg107=reg107*reg30; reg118=ponderation*reg97; reg176=reg176*reg30; reg114=reg114*reg30;
    reg140=reg140*reg30; reg166=reg166*reg30; reg141=reg141*reg30; reg124=ponderation*reg92; reg57=reg57*reg30;
    T tmp_0_5=ponderation*reg174; T tmp_0_1=ponderation*reg66; T tmp_7_3=-reg80; T tmp_1_0=ponderation*reg159; T tmp_4_5=ponderation*reg176;
    T tmp_7_5=ponderation*reg164; T tmp_7_6=-reg111; T tmp_0_4=ponderation*reg172; T tmp_4_6=ponderation*reg170; T tmp_7_7=ponderation*reg74;
    T tmp_5_0=ponderation*reg64; T tmp_4_4=ponderation*reg123; T tmp_0_7=ponderation*reg165; T tmp_0_0=ponderation*reg162; T tmp_0_3=ponderation*reg166;
    T tmp_4_7=ponderation*reg63; T tmp_0_2=ponderation*reg73; T tmp_0_6=ponderation*reg168; T tmp_6_3=ponderation*reg113; T tmp_6_4=ponderation*reg106;
    T tmp_6_5=ponderation*reg153; T tmp_2_5=ponderation*reg175; T tmp_7_0=ponderation*reg107; T tmp_2_6=-reg118; T tmp_7_1=ponderation*reg140;
    T tmp_6_6=ponderation*reg57; T tmp_2_7=ponderation*reg141; T tmp_6_7=-reg109; T tmp_3_0=ponderation*reg145; T tmp_7_4=ponderation*reg148;
    T tmp_3_1=ponderation*reg147; T tmp_3_2=-reg110; T tmp_3_3=ponderation*reg154; T tmp_4_3=ponderation*reg127; T tmp_3_4=ponderation*reg126;
    T tmp_4_2=ponderation*reg132; T tmp_3_5=ponderation*reg130; T tmp_4_1=ponderation*reg135; T tmp_3_6=ponderation*reg133; T tmp_4_0=ponderation*reg137;
    T tmp_3_7=-reg116; T tmp_1_1=ponderation*reg155; T tmp_7_2=ponderation*reg152; T tmp_5_1=ponderation*reg158; T tmp_1_2=ponderation*reg186;
    T tmp_5_2=ponderation*reg71; T tmp_1_3=ponderation*reg190; T tmp_5_3=ponderation*reg31; T tmp_1_4=ponderation*reg72; T tmp_5_4=ponderation*reg184;
    T tmp_1_5=ponderation*reg70; T tmp_5_5=ponderation*reg180; T tmp_1_6=ponderation*reg181; T tmp_5_6=ponderation*reg108; T tmp_1_7=ponderation*reg178;
    T tmp_5_7=ponderation*reg105; T tmp_2_0=ponderation*reg102; T tmp_6_0=ponderation*reg117; T tmp_2_1=ponderation*reg115; T tmp_6_1=ponderation*reg121;
    T tmp_2_2=ponderation*reg119; T tmp_2_3=-reg77; T tmp_6_2=-reg124; T tmp_2_4=ponderation*reg114;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1-var_inter[0]; T reg3=1-var_inter[1];
    T reg4=elem.pos(1)[0]*var_inter[0]; T reg5=reg2*elem.pos(0)[0]; T reg6=1.0/(*f.m).elastic_modulus; T reg7=elem.pos(1)[1]*var_inter[0]; T reg8=reg3*elem.pos(1)[1];
    T reg9=reg2*elem.pos(0)[1]; T reg10=reg3*elem.pos(0)[1]; T reg11=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg12=reg0*reg1; T reg13=reg3*elem.pos(0)[0];
    T reg14=reg3*elem.pos(1)[0]; T reg15=reg6*reg12; reg12=reg11*reg12; T reg16=elem.pos(2)[1]*var_inter[0]; T reg17=reg9+reg7;
    T reg18=reg11*reg1; reg1=reg6*reg1; reg14=reg14-reg13; T reg19=var_inter[1]*elem.pos(2)[0]; T reg20=elem.pos(2)[0]*var_inter[0];
    T reg21=reg4+reg5; reg8=reg8-reg10; T reg22=var_inter[1]*elem.pos(2)[1]; T reg23=reg6*reg15; reg19=reg14+reg19;
    reg14=var_inter[1]*elem.pos(3)[0]; T reg24=reg6*reg1; T reg25=var_inter[1]*elem.pos(3)[1]; reg22=reg8+reg22; reg8=reg11*reg18;
    T reg26=reg11*reg12; reg1=reg11*reg1; reg15=reg11*reg15; reg16=reg16-reg17; reg20=reg20-reg21;
    T reg27=elem.pos(3)[0]*reg2; T reg28=elem.pos(3)[1]*reg2; reg15=reg26+reg15; reg12=reg6*reg12; reg18=reg6*reg18;
    reg24=reg24-reg8; reg1=reg8+reg1; reg27=reg20+reg27; reg23=reg23-reg26; reg19=reg19-reg14;
    reg22=reg22-reg25; reg16=reg28+reg16; reg20=reg6*reg23; reg12=reg26+reg12; reg26=reg11*reg15;
    reg24=reg6*reg24; reg1=reg11*reg1; reg28=reg6*reg0; T reg29=reg22*reg27; T reg30=reg19*reg16;
    T reg31=reg11*reg0; T reg32=reg8+reg18; reg1=reg24-reg1; reg29=reg30-reg29; reg26=reg20-reg26;
    reg20=reg11*reg12; reg24=pow(reg11,2); reg30=pow(reg6,2); T reg33=reg11*reg31; reg32=reg11*reg32;
    T reg34=reg6*reg28; reg22=reg22/reg29; reg27=reg27/reg29; reg19=reg19/reg29; reg20=reg26-reg20;
    reg16=reg16/reg29; reg24=reg30-reg24; reg33=reg34-reg33; reg32=reg1-reg32; reg1=var_inter[1]*reg16;
    reg26=reg3*reg16; reg30=var_inter[0]*reg22; reg34=var_inter[1]*reg27; T reg35=reg2*reg19; T reg36=reg2*reg22;
    T reg37=1-(*f.m).resolution; reg32=reg32/reg20; reg24=reg24/reg33; T reg38=reg26+reg30; T reg39=var_inter[0]*reg19;
    T reg40=reg3*reg27; reg28=reg28/reg33; reg33=reg31/reg33; reg31=(*f.m).resolution*reg24; reg23=reg23/reg20;
    reg15=reg15/reg20; T reg41=reg36+reg1; T reg42=reg35+reg34; T reg43=reg32*reg37; T reg44=reg1-reg30;
    T reg45=0.5*reg41; T reg46=reg39-reg34; T reg47=0.5*reg38; T reg48=reg40+reg39; T reg49=0.5*reg42;
    T reg50=(*f.m).resolution*reg33; T reg51=reg40-reg35; T reg52=reg36-reg26; T reg53=reg15*reg37; T reg54=reg23*reg37;
    reg43=reg31+reg43; reg31=(*f.m).resolution*reg28; T reg55=reg43*reg47; T reg56=0.5*reg44; T reg57=0.5*reg46;
    T reg58=reg43*reg49; reg31=reg54+reg31; reg53=reg50+reg53; reg50=reg43*reg45; reg54=0.5*reg48;
    T reg59=0.5*reg51; T reg60=0.5*reg52; T reg61=reg43*reg59; T reg62=reg43*reg60; T reg63=reg43*reg54;
    T reg64=reg31*reg42; reg55=2*reg55; T reg65=reg43*reg56; T reg66=2*reg50; reg58=2*reg58;
    T reg67=reg43*reg57; T reg68=reg53*reg42; T reg69=reg31*reg41; T reg70=reg53*reg48; T reg71=reg53*reg51;
    T reg72=reg70*reg38; T reg73=reg55*reg54; T reg74=reg53*reg44; reg61=2*reg61; T reg75=reg31*reg44;
    T reg76=reg68*reg41; T reg77=reg31*reg52; T reg78=reg69*reg38; T reg79=reg58*reg54; T reg80=reg64*reg48;
    reg67=2*reg67; T reg81=reg66*reg47; reg65=2*reg65; T reg82=reg53*reg46; T reg83=reg53*reg41;
    T reg84=reg31*reg38; T reg85=reg31*reg51; T reg86=reg53*reg38; T reg87=reg31*reg46; T reg88=reg66*reg49;
    reg62=2*reg62; T reg89=reg31*reg48; T reg90=2*reg63; T reg91=reg82*reg52; T reg92=reg87*reg48;
    T reg93=reg66*reg57; T reg94=reg66*reg59; T reg95=reg58*reg47; T reg96=reg83*reg48; T reg97=reg65*reg59;
    T reg98=reg87*reg46; T reg99=reg65*reg47; T reg100=reg68*reg52; T reg101=reg68*reg44; T reg102=reg85*reg51;
    T reg103=reg89*reg51; T reg104=reg55*reg60; T reg105=reg58*reg59; T reg106=reg58*reg56; T reg107=reg83*reg46;
    T reg108=reg74*reg51; T reg109=reg67*reg60; T reg110=reg69*reg52; T reg111=reg90*reg60; T reg112=reg86*reg51;
    T reg113=reg67*reg57; T reg114=reg75*reg44; T reg115=reg64*reg46; T reg116=reg65*reg56; T reg117=reg66*reg56;
    T reg118=reg60*reg62; reg80=reg81+reg80; T reg119=reg58*reg49; T reg120=reg67*reg54; T reg121=reg75*reg38;
    T reg122=reg59*reg61; T reg123=reg64*reg42; reg73=reg72+reg73; T reg124=reg55*reg59; T reg125=reg90*reg54;
    T reg126=reg38*reg84; T reg127=reg52*reg71; T reg128=reg70*reg52; T reg129=reg66*reg60; reg64=reg64*reg51;
    T reg130=reg59*reg62; T reg131=reg66*reg45; T reg132=reg58*reg60; T reg133=reg83*reg51; T reg134=reg90*reg59;
    T reg135=reg65*reg60; T reg136=reg87*reg51; T reg137=reg52*reg84; T reg138=reg65*reg54; T reg139=reg74*reg48;
    T reg140=reg67*reg47; reg76=reg88+reg76; T reg141=reg58*reg57; T reg142=reg89*reg48; T reg143=reg55*reg47;
    T reg144=reg69*reg44; T reg145=reg67*reg59; T reg146=reg66*reg54; reg68=reg68*reg38; T reg147=reg65*reg57;
    reg79=reg78+reg79; T reg148=reg52*reg77; T reg149=reg75*reg52; T reg150=reg69*reg41; T reg151=reg82*reg38;
    T reg152=reg82*reg44; reg105=reg105-reg110; reg145=reg149+reg145; reg97=reg91+reg97; reg124=reg124-reg128;
    reg123=reg123+reg131; reg91=reg76*reg29; reg137=reg137-reg134; reg139=reg140-reg139; reg141=reg141-reg144;
    reg143=reg143+reg142; reg68=reg68+reg146; reg147=reg152+reg147; reg140=reg79*reg29; reg138=reg151-reg138;
    reg119=reg119+reg150; reg120=reg121-reg120; reg122=reg148+reg122; reg121=reg73*reg29; reg126=reg126+reg125;
    reg64=reg64-reg129; reg109=reg108+reg109; reg130=reg127+reg130; reg132=reg132-reg133; reg135=reg136+reg135;
    reg92=reg99-reg92; reg101=reg101-reg93; reg95=reg95+reg96; reg100=reg100-reg94; reg99=reg80*reg29;
    reg116=reg98+reg116; reg113=reg114+reg113; reg115=reg115-reg117; reg104=reg104-reg103; reg106=reg106-reg107;
    reg118=reg102+reg118; reg112=reg112-reg111; reg123=reg123*reg29; reg115=reg115*reg29; reg130=reg29*reg130;
    reg119=reg119*reg29; reg106=reg106*reg29; reg122=reg29*reg122; reg116=reg116*reg29; reg98=ponderation*reg91;
    reg101=reg101*reg29; reg147=reg147*reg29; reg141=reg141*reg29; reg143=reg143*reg29; reg68=reg68*reg29;
    reg139=reg139*reg29; reg102=ponderation*reg140; reg92=reg92*reg29; reg95=reg95*reg29; reg138=reg138*reg29;
    reg108=ponderation*reg99; reg120=reg120*reg29; reg113=reg113*reg29; reg114=ponderation*reg121; reg104=reg104*reg29;
    reg112=reg112*reg29; reg126=reg126*reg29; reg118=reg118*reg29; reg100=reg100*reg29; reg109=reg109*reg29;
    reg137=reg137*reg29; reg64=reg64*reg29; reg105=reg105*reg29; reg124=reg124*reg29; reg132=reg132*reg29;
    reg135=reg135*reg29; reg97=reg97*reg29; reg145=reg145*reg29; T tmp_4_5=ponderation*reg147; T tmp_1_4=ponderation*reg109;
    T tmp_2_6=-reg102; T tmp_2_5=ponderation*reg138; T tmp_0_1=ponderation*reg130; T tmp_2_4=ponderation*reg120; T tmp_1_5=ponderation*reg135;
    T tmp_2_3=-reg114; T tmp_1_6=ponderation*reg132; T tmp_2_2=ponderation*reg126; T tmp_0_0=ponderation*reg122; T tmp_4_4=ponderation*reg113;
    T tmp_1_7=ponderation*reg64; T tmp_7_7=ponderation*reg123; T tmp_0_2=ponderation*reg137; T tmp_0_3=ponderation*reg124; T tmp_0_4=ponderation*reg145;
    T tmp_6_7=-reg98; T tmp_0_5=ponderation*reg97; T tmp_6_6=ponderation*reg119; T tmp_0_6=ponderation*reg105; T tmp_5_7=ponderation*reg115;
    T tmp_0_7=ponderation*reg100; T tmp_1_1=ponderation*reg118; T tmp_5_6=ponderation*reg106; T tmp_1_2=ponderation*reg112; T tmp_1_3=ponderation*reg104;
    T tmp_5_5=ponderation*reg116; T tmp_3_7=-reg108; T tmp_3_6=ponderation*reg95; T tmp_4_7=ponderation*reg101; T tmp_3_5=ponderation*reg92;
    T tmp_3_4=ponderation*reg139; T tmp_4_6=ponderation*reg141; T tmp_3_3=ponderation*reg143; T tmp_2_7=ponderation*reg68;
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
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg3*reg1; reg1=reg4*reg1; T reg6=reg3*reg2; reg2=reg4*reg2;
    T reg7=reg3*reg2; T reg8=reg4*reg1; T reg9=reg3*reg5; reg1=reg3*reg1; T reg10=reg3*reg6;
    reg2=reg4*reg2; reg2=reg2-reg10; reg7=reg10+reg7; reg6=reg4*reg6; reg5=reg4*reg5;
    reg8=reg8-reg9; reg1=reg9+reg1; T reg11=reg9+reg5; T reg12=reg4*reg2; reg1=reg3*reg1;
    T reg13=reg3*reg7; reg8=reg4*reg8; reg6=reg10+reg6; reg1=reg8-reg1; reg8=reg3*reg6;
    reg11=reg3*reg11; reg13=reg12-reg13; reg11=reg1-reg11; reg8=reg13-reg8; reg7=reg7/reg8;
    reg11=reg11/reg8; reg2=reg2/reg8; reg1=reg2*reg11; reg10=reg7*reg11; reg8=reg6/reg8;
    reg6=(*f.m).alpha*reg2; reg12=reg2*reg1; reg13=1-var_inter[1]; T reg14=1-var_inter[0]; T reg15=(*f.m).alpha*reg7;
    T reg16=reg7*reg10; reg8=(*f.m).alpha*reg8; reg15=reg6+reg15; reg6=reg13*elem.pos(1)[1]; T reg17=reg13*elem.pos(0)[1];
    T reg18=reg13*elem.pos(1)[0]; T reg19=elem.pos(1)[1]*var_inter[0]; T reg20=reg14*elem.pos(0)[1]; T reg21=elem.pos(1)[0]*var_inter[0]; T reg22=reg14*elem.pos(0)[0];
    reg16=reg12-reg16; reg12=reg13*elem.pos(0)[0]; T reg23=reg20+reg19; T reg24=elem.pos(2)[1]*var_inter[0]; reg6=reg6-reg17;
    T reg25=var_inter[1]*elem.pos(2)[1]; T reg26=reg21+reg22; T reg27=elem.pos(2)[0]*var_inter[0]; T reg28=var_inter[1]*elem.pos(2)[0]; reg18=reg18-reg12;
    reg8=reg15+reg8; reg15=reg3*reg0; T reg29=reg4*reg0; reg1=reg1/reg16; reg16=reg10/reg16;
    reg24=reg24-reg23; reg10=reg4*reg29; T reg30=reg3*reg15; T reg31=elem.pos(3)[1]*reg14; T reg32=elem.pos(3)[0]*reg14;
    reg27=reg27-reg26; reg1=reg1*reg8; reg25=reg6+reg25; reg6=var_inter[1]*elem.pos(3)[1]; reg28=reg18+reg28;
    reg18=var_inter[1]*elem.pos(3)[0]; reg8=reg16*reg8; reg16=1-(*f.m).resolution; reg8=reg1-reg8; reg30=reg10-reg30;
    reg24=reg31+reg24; reg28=reg28-reg18; reg25=reg25-reg6; reg32=reg27+reg32; reg29=reg29/reg30;
    reg8=reg16*reg8; reg1=(*f.m).resolution*(*f.m).alpha; reg15=reg15/reg30; reg10=reg25*reg32; reg27=reg28*reg24;
    reg8=reg1+reg8; reg10=reg27-reg10; reg1=(*f.m).resolution*reg29; reg27=(*f.m).resolution*reg15; reg2=reg2*reg16;
    reg7=reg7*reg16; reg25=reg25/reg10; reg1=reg2+reg1; reg7=reg27+reg7; reg32=reg32/reg10;
    reg8=(*f.m).deltaT*reg8; reg24=reg24/reg10; reg28=reg28/reg10; reg2=reg1*reg8; reg27=reg7*reg8;
    reg31=var_inter[1]*reg24; T reg33=var_inter[0]*reg28; T reg34=reg13*reg32; T reg35=reg14*reg25; T reg36=reg13*var_inter[0];
    T reg37=var_inter[0]*reg25; T reg38=reg34+reg33; T reg39=var_inter[1]*reg14; T reg40=reg35+reg31; T reg41=var_inter[1]*reg32;
    T reg42=reg2+reg27; T reg43=reg14*reg28; T reg44=reg13*reg24; T reg45=reg43+reg41; T reg46=reg39*elem.f_vol_e[0];
    T reg47=reg33-reg41; T reg48=reg31-reg37; T reg49=reg36*elem.f_vol_e[1]; T reg50=reg44+reg37; T reg51=reg34-reg43;
    T reg52=reg35-reg44; T reg53=reg42*reg38; T reg54=reg13*reg14; T reg55=var_inter[1]*var_inter[0]; T reg56=reg42*reg40;
    T reg57=reg53-reg49; T reg58=reg36*elem.f_vol_e[0]; T reg59=reg42*reg51; T reg60=reg42*reg48; T reg61=reg42*reg52;
    T reg62=reg55*elem.f_vol_e[0]; T reg63=reg42*reg47; T reg64=reg42*reg50; T reg65=reg55*elem.f_vol_e[1]; T reg66=reg54*elem.f_vol_e[0];
    T reg67=reg56-reg46; T reg68=reg54*elem.f_vol_e[1]; T reg69=reg42*reg45; T reg70=reg39*elem.f_vol_e[1]; T reg71=reg64+reg58;
    T reg72=reg59+reg68; reg57=reg57*reg10; T reg73=reg61+reg66; T reg74=reg60+reg62; T reg75=reg63+reg65;
    reg67=reg67*reg10; T reg76=reg69+reg70; T reg77=reg73*reg10; T reg78=reg72*reg10; T reg79=reg71*reg10;
    reg57=ponderation*reg57; T reg80=reg74*reg10; T reg81=reg75*reg10; reg67=ponderation*reg67; T reg82=reg76*reg10;
    T reg83=ponderation*reg78; sollicitation[indices[0]+1]+=reg83; T reg84=ponderation*reg79; sollicitation[indices[1]+0]+=reg84; T reg85=ponderation*reg77;
    sollicitation[indices[0]+0]+=reg85; sollicitation[indices[1]+1]+=-reg57; reg57=ponderation*reg82; sollicitation[indices[3]+1]+=reg57; T reg86=ponderation*reg80;
    sollicitation[indices[2]+0]+=reg86; T reg87=ponderation*reg81; sollicitation[indices[2]+1]+=reg87; sollicitation[indices[3]+0]+=-reg67;
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
    T reg4=reg0*reg1; T reg5=reg2*reg4; reg4=reg3*reg4; T reg6=reg3*reg1; reg1=reg2*reg1;
    T reg7=reg2*reg5; T reg8=reg3*reg4; T reg9=reg3*reg1; reg1=reg2*reg1; T reg10=reg3*reg6;
    reg5=reg3*reg5; reg6=reg2*reg6; reg7=reg7-reg8; reg4=reg2*reg4; reg1=reg1-reg10;
    reg9=reg10+reg9; reg5=reg8+reg5; reg1=reg2*reg1; reg4=reg8+reg4; reg9=reg3*reg9;
    reg8=reg10+reg6; T reg11=reg2*reg7; T reg12=reg3*reg5; reg8=reg3*reg8; reg9=reg1-reg9;
    reg1=reg3*reg4; reg12=reg11-reg12; reg11=1-var_inter[0]; T reg13=1-var_inter[1]; reg8=reg9-reg8;
    reg1=reg12-reg1; reg9=elem.pos(1)[0]*var_inter[0]; reg12=reg11*elem.pos(0)[0]; T reg14=elem.pos(1)[1]*var_inter[0]; T reg15=reg11*elem.pos(0)[1];
    T reg16=reg13*elem.pos(1)[1]; T reg17=reg13*elem.pos(0)[1]; T reg18=reg13*elem.pos(0)[0]; T reg19=reg13*elem.pos(1)[0]; reg8=reg8/reg1;
    reg5=reg5/reg1; reg7=reg7/reg1; T reg20=elem.pos(2)[1]*var_inter[0]; T reg21=reg15+reg14; T reg22=reg7*reg8;
    T reg23=elem.pos(2)[0]*var_inter[0]; T reg24=reg9+reg12; T reg25=var_inter[1]*elem.pos(2)[1]; reg16=reg16-reg17; T reg26=reg5*reg8;
    T reg27=var_inter[1]*elem.pos(2)[0]; reg19=reg19-reg18; T reg28=reg5*reg26; T reg29=reg7*reg22; T reg30=(*f.m).alpha*reg7;
    reg1=reg4/reg1; reg4=(*f.m).alpha*reg5; reg27=reg19+reg27; reg19=var_inter[1]*elem.pos(3)[0]; T reg31=var_inter[1]*elem.pos(3)[1];
    reg20=reg20-reg21; reg25=reg16+reg25; reg23=reg23-reg24; reg16=elem.pos(3)[1]*reg11; T reg32=elem.pos(3)[0]*reg11;
    T reg33=reg13*vectors[0][indices[0]+1]; T reg34=reg13*vectors[0][indices[1]+1]; reg27=reg27-reg19; reg28=reg29-reg28; reg25=reg25-reg31;
    reg29=reg11*vectors[0][indices[0]+0]; reg32=reg23+reg32; reg20=reg16+reg20; reg16=reg13*vectors[0][indices[0]+0]; reg23=reg13*vectors[0][indices[1]+0];
    T reg35=var_inter[0]*vectors[0][indices[1]+1]; T reg36=reg11*vectors[0][indices[0]+1]; reg4=reg30+reg4; reg30=var_inter[0]*vectors[0][indices[1]+0]; reg1=(*f.m).alpha*reg1;
    T reg37=var_inter[0]*vectors[0][indices[2]+0]; T reg38=var_inter[0]*vectors[0][indices[2]+1]; T reg39=var_inter[1]*vectors[0][indices[2]+0]; reg1=reg4+reg1; reg4=reg27*reg20;
    T reg40=reg25*reg32; reg36=reg35+reg36; reg16=reg23-reg16; reg26=reg26/reg28; reg29=reg30+reg29;
    reg28=reg22/reg28; reg33=reg34-reg33; reg22=var_inter[1]*vectors[0][indices[2]+1]; reg28=reg28*reg1; reg36=reg38-reg36;
    reg1=reg26*reg1; reg23=reg11*vectors[0][indices[3]+1]; reg26=reg2*reg0; reg30=var_inter[1]*vectors[0][indices[3]+1]; reg34=var_inter[1]*vectors[0][indices[3]+0];
    reg22=reg33+reg22; reg16=reg39+reg16; reg33=reg3*reg0; reg29=reg37-reg29; reg40=reg4-reg40;
    reg4=reg11*vectors[0][indices[3]+0]; reg32=reg32/reg40; reg27=reg27/reg40; reg20=reg20/reg40; reg34=reg16-reg34;
    reg16=1-(*f.m).resolution; reg30=reg22-reg30; reg1=reg28-reg1; reg25=reg25/reg40; reg36=reg23+reg36;
    reg4=reg29+reg4; reg22=reg2*reg26; reg23=reg3*reg33; reg28=pow(reg3,2); reg29=pow(reg2,2);
    reg28=reg29-reg28; reg23=reg22-reg23; reg22=reg34*reg32; reg29=reg20*reg30; reg35=reg25*reg36;
    reg37=reg27*reg4; reg38=(*f.m).resolution*(*f.m).alpha; reg1=reg16*reg1; reg22=reg37-reg22; reg35=reg29-reg35;
    reg30=reg32*reg30; reg36=reg27*reg36; reg4=reg25*reg4; reg1=reg38+reg1; reg26=reg26/reg23;
    reg34=reg34*reg20; reg33=reg33/reg23; reg23=reg28/reg23; reg4=reg34-reg4; reg28=(*f.m).resolution*reg23;
    reg29=(*f.m).resolution*reg33; reg5=reg5*reg16; reg8=reg8*reg16; reg16=reg7*reg16; reg7=(*f.m).resolution*reg26;
    reg22=reg35+reg22; reg30=reg36-reg30; reg1=(*f.m).deltaT*reg1; reg34=var_inter[1]*reg20; reg35=var_inter[1]*reg32;
    reg22=0.5*reg22; reg36=var_inter[0]*reg25; reg37=reg11*reg27; reg38=var_inter[0]*reg27; reg39=reg13*reg32;
    reg30=reg30-reg1; reg7=reg16+reg7; reg16=reg11*reg25; reg4=reg4-reg1; reg5=reg29+reg5;
    reg8=reg28+reg8; reg28=reg13*reg20; reg29=reg16+reg34; reg22=reg8*reg22; T reg41=reg16-reg28;
    T reg42=reg39-reg37; T reg43=reg34-reg36; T reg44=reg7*reg30; T reg45=reg5*reg4; reg30=reg5*reg30;
    reg4=reg7*reg4; T reg46=reg38-reg35; T reg47=reg37+reg35; T reg48=reg28+reg36; T reg49=reg39+reg38;
    T reg50=0.5*reg46; T reg51=0.5*reg41; T reg52=0.5*reg48; T reg53=0.5*reg49; T reg54=0.5*reg29;
    T reg55=0.5*reg43; T reg56=0.5*reg42; reg22=2*reg22; reg44=reg45+reg44; reg30=reg4+reg30;
    reg4=0.5*reg47; reg45=reg22*reg54; T reg57=reg44*reg47; T reg58=reg22*reg53; T reg59=reg22*reg4;
    T reg60=reg44*reg49; T reg61=reg30*reg29; T reg62=reg22*reg52; T reg63=reg30*reg43; T reg64=reg22*reg55;
    T reg65=reg22*reg50; T reg66=reg44*reg46; T reg67=reg44*reg42; T reg68=reg22*reg51; T reg69=reg22*reg56;
    T reg70=reg30*reg41; T reg71=var_inter[1]*var_inter[0]; T reg72=reg13*var_inter[0]; T reg73=var_inter[1]*reg11; T reg74=reg13*reg11;
    T reg75=reg30*reg48; reg59=reg59-reg61; T reg76=reg73*elem.f_vol_e[0]; T reg77=reg71*elem.f_vol_e[1]; reg64=reg66+reg64;
    reg57=reg57-reg45; reg69=reg70+reg69; reg66=reg73*elem.f_vol_e[1]; reg70=reg74*elem.f_vol_e[0]; T reg78=reg71*elem.f_vol_e[0];
    reg65=reg63+reg65; reg75=reg75-reg58; reg63=reg72*elem.f_vol_e[0]; reg68=reg67+reg68; reg67=reg74*elem.f_vol_e[1];
    reg62=reg62-reg60; T reg79=reg72*elem.f_vol_e[1]; reg57=reg57-reg66; reg75=reg75-reg63; reg59=reg59-reg76;
    reg68=reg68-reg67; reg62=reg62-reg79; reg69=reg69-reg70; reg64=reg64-reg77; reg65=reg65-reg78;
    reg57=reg40*reg57; reg69=reg40*reg69; reg65=reg40*reg65; reg75=reg40*reg75; reg59=reg40*reg59;
    reg68=reg40*reg68; reg62=reg40*reg62; reg64=reg40*reg64; sollicitation[indices[0]+1]+=ponderation*reg68; sollicitation[indices[2]+0]+=ponderation*reg65;
    sollicitation[indices[3]+0]+=ponderation*reg59; sollicitation[indices[1]+0]+=ponderation*reg75; sollicitation[indices[0]+0]+=ponderation*reg69; sollicitation[indices[1]+1]+=ponderation*reg62; sollicitation[indices[2]+1]+=ponderation*reg64;
    sollicitation[indices[3]+1]+=ponderation*reg57;
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

