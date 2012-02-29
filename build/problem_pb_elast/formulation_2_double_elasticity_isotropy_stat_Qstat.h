
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
    node.dep[0]=vecs[0][indice+0]; node.dep[1]=vecs[0][indice+1];
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
    node.dep[0]=vecs[0][indice+0]; node.dep[1]=vecs[0][indice+1];
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg1*reg0; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg4*reg1; T reg6=reg2*reg3; reg1=reg1*reg3; reg2=reg4*reg2;
    T reg7=reg1*reg3; T reg8=reg5*reg3; T reg9=reg2*reg3; reg5=reg4*reg5; T reg10=reg6*reg3;
    reg2=reg4*reg2; reg8=reg7+reg8; reg9=reg9+reg10; reg1=reg4*reg1; reg2=reg2-reg10;
    reg6=reg4*reg6; reg5=reg5-reg7; reg5=reg4*reg5; reg10=reg6+reg10; reg6=reg9*reg3;
    T reg11=reg7+reg1; T reg12=reg4*reg2; reg8=reg3*reg8; T reg13=elem.pos(1)[1]-elem.pos(0)[1]; T reg14=elem.pos(1)[0]-elem.pos(0)[0];
    T reg15=elem.pos(2)[1]-elem.pos(0)[1]; T reg16=reg10*reg3; reg8=reg5-reg8; reg6=reg12-reg6; reg11=reg3*reg11;
    reg5=elem.pos(2)[0]-elem.pos(0)[0]; reg12=reg5*reg13; reg16=reg6-reg16; reg11=reg8-reg11; reg6=reg15*reg14;
    reg9=reg9/reg16; reg2=reg2/reg16; reg11=reg11/reg16; reg12=reg6-reg12; reg6=vectors[0][indices[2]+1]-vectors[0][indices[0]+1];
    reg5=reg5/reg12; reg8=reg3*reg0; T reg17=reg9*reg11; T reg18=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; T reg19=vectors[0][indices[2]+0]-vectors[0][indices[0]+0];
    T reg20=reg2*reg11; T reg21=reg4*reg0; T reg22=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg14=reg14/reg12; reg13=reg13/reg12;
    reg15=reg15/reg12; T reg23=reg9*reg17; T reg24=reg22*reg5; T reg25=reg2*reg20; T reg26=reg15*reg18;
    T reg27=reg4*reg21; T reg28=reg6*reg14; T reg29=reg19*reg13; T reg30=(*f.m).alpha*reg9; T reg31=reg8*reg3;
    T reg32=reg2*(*f.m).alpha; reg16=reg10/reg16; reg31=reg27-reg31; reg16=(*f.m).alpha*reg16; reg29=reg26-reg29;
    elem.epsilon[0][0]=reg29; reg30=reg32+reg30; reg24=reg28-reg24; elem.epsilon[0][1]=reg24; reg10=(*f.m).deltaT*(*f.m).alpha;
    reg23=reg25-reg23; reg8=reg8/reg31; reg21=reg21/reg31; reg17=reg17/reg23; reg23=reg20/reg23;
    reg30=reg16+reg30; reg16=reg29-reg10; reg20=reg24-reg10; reg25=reg16*reg21; reg26=reg20*reg8;
    reg23=reg30*reg23; reg17=reg30*reg17; reg16=reg16*reg8; reg20=reg20*reg21; reg26=reg25+reg26;
    reg20=reg16+reg20; reg16=1-(*f.m).resolution; reg17=reg23-reg17; reg23=reg4*reg20; reg25=reg4*reg26;
    reg27=pow(reg4,2); reg20=reg20*reg3; reg26=reg26*reg3; reg28=PNODE(1).dep[1]-PNODE(0).dep[1]; reg30=PNODE(2).dep[1]-PNODE(0).dep[1];
    reg32=pow(reg3,2); T reg33=PNODE(1).dep[0]-PNODE(0).dep[0]; T reg34=(*f.m).alpha*(*f.m).resolution; reg17=reg16*reg17; T reg35=PNODE(2).dep[0]-PNODE(0).dep[0];
    T reg36=reg33*reg15; T reg37=reg30*reg13; reg25=reg25-reg20; reg23=reg23-reg26; reg30=reg30*reg14;
    T reg38=reg28*reg5; reg32=reg27-reg32; reg27=reg35*reg14; reg28=reg28*reg15; reg33=reg33*reg5;
    reg17=reg34+reg17; reg35=reg35*reg13; reg35=reg36-reg35; reg38=reg30-reg38; reg26=reg20+reg26;
    reg37=reg28-reg37; reg31=reg32/reg31; reg17=(*f.m).deltaT*reg17; reg19=reg19*reg14; reg18=reg5*reg18;
    reg22=reg22*reg15; reg6=reg6*reg13; reg9=reg16*reg9; reg33=reg27-reg33; reg2=reg2*reg16;
    reg20=(*f.m).resolution*reg21; reg27=(*f.m).resolution*reg8; reg23=reg10+reg23; reg25=reg25+reg10; reg18=reg19-reg18;
    reg19=reg25+reg23; reg28=(*f.m).resolution*reg31; reg6=reg22-reg6; reg9=reg27+reg9; reg20=reg2+reg20;
    reg2=reg35-reg17; reg22=reg38-reg17; reg37=reg33+reg37; reg26=reg10-reg26; reg11=reg16*reg11;
    reg19=reg26+reg19; reg28=reg11+reg28; reg11=reg20*reg2; reg16=reg20*reg22; reg27=reg9*reg2;
    reg30=reg9*reg22; reg37=0.5*reg37; reg18=reg6+reg18; reg30=reg11+reg30; reg19=reg19/3;
    reg16=reg27+reg16; reg6=reg28*reg37; reg18=0.5*reg18; elem.epsilon[0][2]=reg18; reg6=2*reg6;
    reg2=reg2*reg30; reg23=reg23-reg19; reg11=reg18*reg31; reg25=reg25-reg19; reg22=reg22*reg16;
    reg11=reg0*reg11; reg23=pow(reg23,2); reg27=reg6*reg37; reg2=reg22+reg2; reg25=pow(reg25,2);
    reg19=reg26-reg19; reg2=reg27+reg2; reg22=2*reg11; reg19=pow(reg19,2); reg23=reg25+reg23;
    reg22=reg11*reg22; reg23=reg19+reg23; reg12=reg2*reg12; reg29=reg29-reg17; reg23=reg22+reg23;
    reg17=reg24-reg17; reg2=0.16666666666666665741*reg12; reg12=0.33333333333333331483*reg12; reg11=reg17*reg9; reg17=reg17*reg20;
    reg23=1.5*reg23; reg12=reg2+reg12; reg20=reg20*reg29; reg29=reg9*reg29; elem.sigma[0][2]=reg28*reg18;
    elem.sigma[0][0]=reg11+reg20; elem.sigma[0][1]=reg29+reg17; elem.ener=reg12/2; elem.sigma_von_mises=pow(reg23,0.5);
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg1*reg0; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg1*reg4; reg1=reg3*reg1; T reg6=reg3*reg2; reg2=reg2*reg4;
    T reg7=reg3*reg1; T reg8=reg5*reg4; T reg9=reg2*reg4; T reg10=reg3*reg6; reg6=reg6*reg4;
    reg1=reg1*reg4; reg5=reg3*reg5; reg10=reg10-reg9; reg7=reg7-reg8; reg2=reg3*reg2;
    reg1=reg8+reg1; reg6=reg6+reg9; T reg11=reg3*reg10; reg9=reg2+reg9; reg2=reg6*reg4;
    T reg12=reg8+reg5; reg7=reg3*reg7; reg1=reg4*reg1; reg12=reg4*reg12; T reg13=reg9*reg4;
    reg2=reg11-reg2; reg1=reg7-reg1; reg13=reg2-reg13; reg12=reg1-reg12; reg6=reg6/reg13;
    reg10=reg10/reg13; reg12=reg12/reg13; reg1=reg6*reg12; reg2=reg10*reg12; reg7=reg10*(*f.m).alpha;
    reg11=(*f.m).alpha*reg6; T reg14=reg6*reg1; reg13=reg9/reg13; reg9=reg10*reg2; reg13=(*f.m).alpha*reg13;
    reg11=reg7+reg11; reg14=reg9-reg14; reg1=reg1/reg14; reg7=reg3*reg0; reg11=reg13+reg11;
    reg9=reg4*reg0; reg14=reg2/reg14; reg2=pow(reg3,2); reg13=elem.pos(1)[0]-elem.pos(0)[0]; T reg15=elem.pos(1)[1]-elem.pos(0)[1];
    reg1=reg11*reg1; reg14=reg11*reg14; reg11=reg3*reg7; T reg16=pow(reg4,2); T reg17=elem.pos(2)[0]-elem.pos(0)[0];
    T reg18=reg9*reg4; T reg19=elem.pos(2)[1]-elem.pos(0)[1]; reg16=reg2-reg16; reg2=1-(*f.m).resolution; reg1=reg14-reg1;
    reg18=reg11-reg18; reg11=reg19*reg13; reg14=reg17*reg15; reg1=reg2*reg1; reg16=reg16/reg18;
    T reg20=(*f.m).alpha*(*f.m).resolution; reg14=reg11-reg14; reg9=reg9/reg18; reg18=reg7/reg18; reg1=reg20+reg1;
    reg6=reg2*reg6; reg12=reg2*reg12; reg2=reg10*reg2; reg17=reg17/reg14; reg19=reg19/reg14;
    reg13=reg13/reg14; reg7=(*f.m).resolution*reg9; reg10=(*f.m).resolution*reg18; reg15=reg15/reg14; reg11=(*f.m).resolution*reg16;
    reg10=reg2+reg10; reg6=reg7+reg6; reg2=0.5*reg19; reg11=reg12+reg11; reg7=0.5*reg15;
    reg12=0.5*reg13; reg20=0.5*reg17; T reg21=reg17-reg13; reg1=(*f.m).deltaT*reg1; T reg22=reg15-reg19;
    T reg23=reg11*reg20; T reg24=reg6*reg1; T reg25=0.5*reg21; T reg26=reg10*reg1; T reg27=0.5*reg22;
    T reg28=reg11*reg2; T reg29=reg11*reg12; T reg30=reg11*reg7; T reg31=reg6*reg17; reg28=2*reg28;
    T reg32=reg6*reg15; T reg33=reg10*reg19; T reg34=reg11*reg27; T reg35=2*reg30; reg29=2*reg29;
    T reg36=reg10*reg13; T reg37=reg6*reg13; T reg38=reg6*reg19; T reg39=reg10*reg15; T reg40=1-var_inter[0];
    T reg41=reg11*reg25; T reg42=2*reg23; T reg43=reg10*reg17; T reg44=reg24+reg26; T reg45=reg29*reg20;
    T reg46=reg44*reg15; T reg47=reg22*reg6; T reg48=reg6*reg21; reg34=2*reg34; T reg49=var_inter[1]*(*f.m).f_vol[0];
    T reg50=reg22*reg10; T reg51=reg19*reg39; T reg52=var_inter[0]*(*f.m).f_vol[1]; T reg53=reg37*reg15; T reg54=reg7*reg29;
    T reg55=reg35*reg12; T reg56=reg36*reg17; T reg57=reg10*reg21; T reg58=reg32*reg13; T reg59=reg2*reg35;
    T reg60=reg44*reg17; T reg61=reg28*reg20; T reg62=reg42*reg2; T reg63=reg31*reg19; reg40=reg40-var_inter[1];
    T reg64=reg7*reg28; T reg65=reg43*reg13; T reg66=reg33*reg15; T reg67=reg42*reg12; reg41=2*reg41;
    T reg68=reg38*reg17; T reg69=reg44*reg13; T reg70=reg41*reg27; T reg71=reg40*(*f.m).f_vol[1]; T reg72=reg37*reg19;
    T reg73=reg42*reg27; T reg74=reg21*reg38; T reg75=reg47*reg21; T reg76=reg48*reg15; T reg77=reg34*reg12;
    T reg78=reg50*reg15; T reg79=reg41*reg12; T reg80=reg22*reg39; T reg81=reg35*reg7; T reg82=reg36*reg13;
    T reg83=reg42*reg25; T reg84=reg22*reg31; reg64=reg65+reg64; T reg85=reg25*reg28; T reg86=reg57*reg17;
    T reg87=reg42*reg7; T reg88=reg2*reg34; T reg89=reg43*reg21; T reg90=reg35*reg27; T reg91=reg46-reg49;
    T reg92=reg44*reg19; T reg93=reg12*reg29; T reg94=reg50*reg19; reg68=reg62+reg68; T reg95=reg44*reg21;
    T reg96=reg25*reg29; reg53=reg55+reg53; T reg97=reg34*reg27; reg45=reg51+reg45; T reg98=reg25*reg34;
    T reg99=reg22*reg48; T reg100=reg39*reg15; T reg101=reg57*reg21; T reg102=reg25*reg35; reg37=reg22*reg37;
    reg48=reg48*reg19; T reg103=reg34*reg20; T reg104=reg44*reg22; T reg105=reg27*reg29; T reg106=reg41*reg2;
    reg38=reg38*reg13; T reg107=reg41*reg7; T reg108=reg22*reg33; reg56=reg56+reg59; reg34=reg34*reg7;
    reg57=reg57*reg13; T reg109=reg12*reg28; T reg110=reg31*reg15; T reg111=reg47*reg13; reg33=reg33*reg19;
    reg54=reg54+reg58; reg47=reg47*reg17; T reg112=reg42*reg20; reg50=reg22*reg50; T reg113=var_inter[0]*(*f.m).f_vol[0];
    T reg114=reg32*reg17; T reg115=reg27*reg28; reg29=reg2*reg29; reg61=reg63+reg61; T reg116=reg41*reg25;
    T reg117=reg32*reg21; reg41=reg41*reg20; reg66=reg67+reg66; T reg118=reg43*reg17; T reg119=reg60-reg52;
    T reg120=reg35*reg20; T reg121=reg40*(*f.m).f_vol[0]; reg28=reg2*reg28; T reg122=var_inter[1]*(*f.m).f_vol[1]; reg36=reg36*reg21;
    reg108=reg108-reg83; reg116=reg50+reg116; reg50=reg95+reg71; reg33=reg33+reg112; T reg123=reg92+reg113;
    reg119=reg119*reg14; reg36=reg36-reg90; T reg124=reg53*reg14; reg115=reg115-reg89; T reg125=reg56*reg14;
    reg91=reg91*reg14; T reg126=reg64*reg14; reg74=reg74-reg73; reg41=reg94-reg41; reg94=reg104+reg121;
    reg98=reg99+reg98; reg103=reg48-reg103; reg38=reg38+reg87; reg86=reg88-reg86; reg82=reg82+reg81;
    reg72=reg72+reg120; reg48=reg14*reg61; reg76=reg77-reg76; reg78=reg79-reg78; reg37=reg37-reg102;
    reg28=reg28+reg118; reg97=reg101+reg97; reg85=reg85-reg84; reg29=reg29+reg114; reg77=reg122+reg69;
    reg107=reg111-reg107; reg70=reg75+reg70; reg105=reg105-reg117; reg75=reg66*reg14; reg34=reg57-reg34;
    reg57=reg68*reg14; reg79=reg45*reg14; reg109=reg109+reg110; reg47=reg106-reg47; reg96=reg96-reg80;
    reg88=reg54*reg14; reg93=reg93+reg100; reg47=reg47*reg14; reg99=reg123*reg14; reg37=reg37*reg14;
    reg96=reg96*reg14; reg101=ponderation*reg124; reg85=reg14*reg85; reg86=reg86*reg14; reg36=reg36*reg14;
    reg29=reg29*reg14; reg109=reg14*reg109; reg106=ponderation*reg75; reg70=reg70*reg14; reg93=reg93*reg14;
    reg33=reg33*reg14; reg111=ponderation*reg57; T reg127=ponderation*reg88; reg107=reg107*reg14; reg116=reg116*reg14;
    T reg128=ponderation*reg79; T reg129=ponderation*reg48; reg34=reg34*reg14; T reg130=reg77*reg14; reg105=reg105*reg14;
    reg97=reg97*reg14; reg78=reg78*reg14; reg72=reg72*reg14; reg28=reg28*reg14; reg82=reg82*reg14;
    reg76=reg76*reg14; reg115=reg115*reg14; T reg131=ponderation*reg126; reg108=reg108*reg14; reg74=reg74*reg14;
    reg41=reg41*reg14; reg38=reg38*reg14; reg91=ponderation*reg91; reg98=reg98*reg14; T reg132=reg94*reg14;
    reg103=reg103*reg14; T reg133=reg50*reg14; T reg134=ponderation*reg125; reg119=ponderation*reg119; T tmp_5_0=ponderation*reg107;
    T tmp_2_0=ponderation*reg41; T tmp_5_1=ponderation*reg34; T tmp_5_4=-reg127; T tmp_4_2=-reg106; T tmp_1_4=ponderation*reg105;
    T tmp_0_1=ponderation*reg98; T tmp_2_5=ponderation*reg72; T tmp_3_4=ponderation*reg29; T tmp_1_3=ponderation*reg115; reg29=reg99*ponderation;
    T vec_2=reg29; T tmp_4_5=-reg101; T tmp_3_3=ponderation*reg28; T tmp_2_1=ponderation*reg103; T tmp_0_5=ponderation*reg37;
    T tmp_5_3=-reg131; T tmp_3_1=ponderation*reg86; T tmp_0_2=ponderation*reg108; T tmp_5_2=ponderation*reg38; T tmp_1_2=ponderation*reg74;
    T vec_3=-reg119; T tmp_1_5=ponderation*reg36; reg28=ponderation*reg133; T vec_1=reg28; T tmp_3_0=ponderation*reg47;
    reg34=ponderation*reg132; T vec_0=reg34; T tmp_4_4=ponderation*reg93; T tmp_0_4=ponderation*reg96; T tmp_3_2=-reg111;
    T tmp_4_3=ponderation*reg109; T tmp_0_0=ponderation*reg116; T tmp_2_4=-reg128; T tmp_1_0=ponderation*reg70; T tmp_2_3=-reg129;
    reg36=ponderation*reg130; T vec_5=reg36; T vec_4=-reg91; T tmp_1_1=ponderation*reg97; T tmp_4_0=ponderation*reg78;
    T tmp_0_3=ponderation*reg85; T tmp_3_5=-reg134; T tmp_2_2=ponderation*reg33; T tmp_5_5=ponderation*reg82; T tmp_4_1=ponderation*reg76;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg1*reg0; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg1*reg4; reg1=reg3*reg1; T reg6=reg3*reg2; reg2=reg2*reg4;
    T reg7=reg3*reg1; T reg8=reg5*reg4; T reg9=reg2*reg4; T reg10=reg3*reg6; reg6=reg6*reg4;
    reg1=reg1*reg4; reg5=reg3*reg5; reg10=reg10-reg9; reg7=reg7-reg8; reg2=reg3*reg2;
    reg1=reg8+reg1; reg6=reg6+reg9; T reg11=reg3*reg10; reg9=reg2+reg9; reg2=reg6*reg4;
    T reg12=reg8+reg5; reg7=reg3*reg7; reg1=reg4*reg1; reg12=reg4*reg12; T reg13=reg9*reg4;
    reg2=reg11-reg2; reg1=reg7-reg1; reg13=reg2-reg13; reg12=reg1-reg12; reg10=reg10/reg13;
    reg12=reg12/reg13; reg6=reg6/reg13; reg1=reg10*reg12; reg2=reg6*reg12; reg7=reg6*reg2;
    reg11=(*f.m).alpha*reg6; T reg14=reg10*(*f.m).alpha; reg13=reg9/reg13; reg9=reg10*reg1; reg13=(*f.m).alpha*reg13;
    reg11=reg14+reg11; reg7=reg9-reg7; reg1=reg1/reg7; reg7=reg2/reg7; reg2=reg3*reg0;
    reg9=reg4*reg0; reg11=reg13+reg11; reg13=pow(reg4,2); reg14=elem.pos(2)[0]-elem.pos(0)[0]; T reg15=reg3*reg2;
    T reg16=reg9*reg4; reg1=reg11*reg1; T reg17=elem.pos(1)[1]-elem.pos(0)[1]; T reg18=elem.pos(1)[0]-elem.pos(0)[0]; reg7=reg11*reg7;
    reg11=elem.pos(2)[1]-elem.pos(0)[1]; T reg19=pow(reg3,2); T reg20=reg11*reg18; T reg21=reg14*reg17; T reg22=1-(*f.m).resolution;
    reg7=reg1-reg7; reg13=reg19-reg13; reg16=reg15-reg16; reg7=reg22*reg7; reg13=reg13/reg16;
    reg2=reg2/reg16; reg16=reg9/reg16; reg21=reg20-reg21; reg1=(*f.m).alpha*(*f.m).resolution; reg17=reg17/reg21;
    reg18=reg18/reg21; reg14=reg14/reg21; reg10=reg10*reg22; reg9=(*f.m).resolution*reg13; reg7=reg1+reg7;
    reg1=(*f.m).resolution*reg16; reg15=(*f.m).resolution*reg2; reg11=reg11/reg21; reg12=reg22*reg12; reg6=reg22*reg6;
    reg9=reg12+reg9; reg15=reg10+reg15; reg6=reg1+reg6; reg1=0.5*reg11; reg10=0.5*reg17;
    reg7=(*f.m).deltaT*reg7; reg12=0.5*reg18; reg19=reg17-reg11; reg20=reg14-reg18; reg22=reg6*reg7;
    T reg23=reg9*reg10; T reg24=reg9*reg1; T reg25=0.5*reg14; T reg26=0.5*reg20; T reg27=0.5*reg19;
    T reg28=reg9*reg12; T reg29=reg15*reg7; T reg30=reg15*reg18; T reg31=reg15*reg17; T reg32=reg22+reg29;
    reg24=2*reg24; T reg33=reg6*reg14; reg28=2*reg28; T reg34=reg6*reg18; T reg35=reg9*reg26;
    T reg36=reg9*reg27; T reg37=1-var_inter[0]; T reg38=2*reg23; T reg39=reg9*reg25; T reg40=reg38*reg12;
    T reg41=reg34*reg17; reg36=2*reg36; T reg42=var_inter[1]*(*f.m).f_vol[0]; T reg43=reg6*reg11; T reg44=reg32*reg17;
    reg37=reg37-var_inter[1]; T reg45=reg6*reg20; T reg46=var_inter[0]*(*f.m).f_vol[1]; T reg47=reg15*reg14; T reg48=reg24*reg25;
    T reg49=reg19*reg15; T reg50=2*reg39; reg35=2*reg35; T reg51=reg1*reg38; T reg52=reg15*reg20;
    T reg53=reg33*reg11; T reg54=reg30*reg14; T reg55=reg15*reg11; T reg56=reg28*reg25; T reg57=reg32*reg14;
    T reg58=reg11*reg31; T reg59=reg6*reg17; reg41=reg40+reg41; T reg60=reg31*reg17; T reg61=reg34*reg11;
    T reg62=reg38*reg25; T reg63=reg19*reg49; T reg64=reg35*reg26; T reg65=reg26*reg36; T reg66=reg27*reg24;
    T reg67=reg19*reg45; T reg68=reg32*reg20; T reg69=var_inter[0]*(*f.m).f_vol[0]; T reg70=reg32*reg11; T reg71=reg26*reg28;
    reg56=reg58+reg56; T reg72=reg44-reg42; T reg73=reg47*reg20; T reg74=reg50*reg26; T reg75=reg36*reg27;
    T reg76=reg52*reg20; T reg77=reg57-reg46; T reg78=reg19*reg33; T reg79=reg26*reg24; T reg80=reg32*reg19;
    T reg81=reg38*reg27; T reg82=reg12*reg28; T reg83=reg1*reg24; T reg84=reg47*reg14; T reg85=reg19*reg31;
    T reg86=var_inter[1]*(*f.m).f_vol[1]; T reg87=reg1*reg28; T reg88=reg59*reg14; T reg89=reg30*reg20; T reg90=reg37*(*f.m).f_vol[0];
    reg30=reg30*reg18; T reg91=reg38*reg10; T reg92=reg27*reg28; reg48=reg53+reg48; T reg93=reg20*reg43;
    T reg94=reg50*reg27; T reg95=reg50*reg25; T reg96=reg59*reg20; T reg97=reg37*(*f.m).f_vol[1]; T reg98=reg26*reg38;
    reg34=reg19*reg34; T reg99=reg19*reg55; reg54=reg54+reg51; T reg100=reg32*reg18; T reg101=reg55*reg11;
    reg89=reg89-reg81; reg66=reg66-reg73; reg83=reg83+reg84; T reg102=reg80+reg90; T reg103=reg41*reg21;
    T reg104=reg86+reg100; reg87=reg87+reg88; reg64=reg63+reg64; reg65=reg67+reg65; reg63=reg68+reg97;
    reg34=reg34-reg98; reg82=reg82+reg60; reg67=reg21*reg48; reg71=reg71-reg85; reg101=reg101+reg95;
    T reg105=reg70+reg69; reg77=reg77*reg21; T reg106=reg54*reg21; T reg107=reg56*reg21; reg75=reg76+reg75;
    reg72=reg72*reg21; reg61=reg61+reg62; reg92=reg92-reg96; reg93=reg93-reg94; reg99=reg99-reg74;
    reg79=reg79-reg78; reg30=reg30+reg91; reg76=reg105*reg21; reg79=reg21*reg79; reg71=reg71*reg21;
    reg66=reg66*reg21; reg65=reg65*reg21; T reg108=ponderation*reg106; reg101=reg101*reg21; reg92=reg92*reg21;
    reg61=reg61*reg21; reg83=reg83*reg21; T reg109=ponderation*reg67; reg30=reg30*reg21; T reg110=ponderation*reg107;
    reg75=reg75*reg21; reg93=reg93*reg21; T reg111=reg104*reg21; reg99=reg99*reg21; T reg112=ponderation*reg103;
    reg89=reg89*reg21; reg72=ponderation*reg72; reg77=ponderation*reg77; T reg113=reg63*reg21; reg87=reg87*reg21;
    T reg114=reg102*reg21; reg34=reg34*reg21; reg82=reg82*reg21; reg64=reg64*reg21; T tmp_3_4=ponderation*reg87;
    T tmp_2_4=-reg110; T tmp_1_2=ponderation*reg93; T tmp_2_3=-reg109; T tmp_4_5=-reg112; reg87=reg76*ponderation;
    T vec_2=reg87; T tmp_0_2=ponderation*reg99; T tmp_2_5=ponderation*reg61; T tmp_0_5=ponderation*reg34; T tmp_1_5=ponderation*reg89;
    T vec_3=-reg77; reg34=ponderation*reg113; T vec_1=reg34; T tmp_0_4=ponderation*reg71; reg61=ponderation*reg114;
    T vec_0=reg61; T tmp_4_4=ponderation*reg82; T tmp_0_0=ponderation*reg64; reg64=ponderation*reg111; T vec_5=reg64;
    T vec_4=-reg72; T tmp_1_1=ponderation*reg75; T tmp_3_5=-reg108; T tmp_5_5=ponderation*reg30; T tmp_0_3=ponderation*reg79;
    T tmp_3_3=ponderation*reg83; T tmp_2_2=ponderation*reg101; T tmp_1_4=ponderation*reg92; T tmp_0_1=ponderation*reg65; T tmp_1_3=ponderation*reg66;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=reg1*reg0; T reg5=reg3*reg1; T reg6=reg4*reg2; reg1=reg1*reg2; reg4=reg3*reg4;
    T reg7=reg1*reg2; T reg8=reg5*reg2; T reg9=reg4*reg2; T reg10=reg6*reg2; reg5=reg3*reg5;
    reg4=reg3*reg4; reg9=reg9+reg10; reg1=reg3*reg1; reg4=reg4-reg10; reg5=reg5-reg7;
    reg6=reg3*reg6; reg8=reg7+reg8; reg5=reg3*reg5; reg8=reg2*reg8; T reg11=reg3*reg4;
    T reg12=reg3*reg0; reg10=reg6+reg10; reg6=reg2*reg0; T reg13=reg7+reg1; T reg14=reg9*reg2;
    T reg15=pow(reg3,2); T reg16=reg10*reg2; T reg17=elem.pos(1)[1]-elem.pos(0)[1]; reg8=reg5-reg8; reg5=pow(reg2,2);
    T reg18=reg3*reg12; T reg19=elem.pos(1)[0]-elem.pos(0)[0]; reg14=reg11-reg14; reg13=reg2*reg13; reg11=reg6*reg2;
    T reg20=elem.pos(2)[1]-elem.pos(0)[1]; T reg21=elem.pos(2)[0]-elem.pos(0)[0]; reg13=reg8-reg13; reg11=reg18-reg11; reg8=reg21*reg17;
    reg18=reg20*reg19; reg16=reg14-reg16; reg5=reg15-reg5; reg8=reg18-reg8; reg13=reg13/reg16;
    reg14=1-(*f.m).resolution; reg5=reg5/reg11; reg9=reg9/reg16; reg15=reg14*reg13; reg12=reg12/reg11;
    reg4=reg4/reg16; reg20=reg20/reg8; reg21=reg21/reg8; reg18=(*f.m).resolution*reg5; reg11=reg6/reg11;
    reg19=reg19/reg8; reg17=reg17/reg8; reg6=reg17-reg20; T reg22=0.5*reg21; T reg23=reg4*reg14;
    T reg24=0.5*reg19; T reg25=0.5*reg17; T reg26=reg21-reg19; T reg27=0.5*reg20; reg18=reg15+reg18;
    reg15=reg14*reg9; T reg28=(*f.m).resolution*reg11; T reg29=(*f.m).resolution*reg12; T reg30=reg18*reg22; reg29=reg23+reg29;
    reg15=reg28+reg15; reg23=0.5*reg26; reg28=reg18*reg24; T reg31=reg18*reg27; T reg32=reg18*reg25;
    T reg33=0.5*reg6; T reg34=reg29*reg20; T reg35=reg29*reg19; T reg36=reg15*reg20; T reg37=reg15*reg17;
    T reg38=2*reg32; T reg39=reg18*reg33; T reg40=reg18*reg23; T reg41=reg29*reg17; T reg42=reg29*reg21;
    T reg43=2*reg30; T reg44=reg15*reg19; reg28=2*reg28; T reg45=reg15*reg21; reg31=2*reg31;
    T reg46=reg43*reg27; T reg47=reg31*reg22; T reg48=reg45*reg20; T reg49=reg25*reg28; T reg50=reg43*reg24;
    T reg51=reg34*reg17; T reg52=reg35*reg21; T reg53=reg15*reg26; T reg54=reg37*reg19; reg39=2*reg39;
    T reg55=reg38*reg24; T reg56=reg44*reg17; T reg57=reg6*reg15; T reg58=reg42*reg19; T reg59=reg36*reg21;
    T reg60=reg25*reg31; T reg61=reg6*reg29; reg40=2*reg40; T reg62=reg29*reg26; T reg63=reg27*reg38;
    T reg64=reg20*reg41; T reg65=reg28*reg22; T reg66=reg40*reg22; T reg67=reg43*reg25; T reg68=reg40*reg23;
    reg59=reg46+reg59; reg56=reg55+reg56; T reg69=reg44*reg20; T reg70=reg57*reg26; T reg71=reg23*reg39;
    T reg72=reg33*reg31; T reg73=reg6*reg53; T reg74=reg23*reg38; T reg75=reg35*reg19; reg44=reg6*reg44;
    T reg76=reg40*reg33; T reg77=reg38*reg25; T reg78=reg43*reg33; T reg79=reg26*reg36; T reg80=reg37*reg26;
    T reg81=reg53*reg17; T reg82=reg39*reg24; T reg83=reg61*reg17; T reg84=reg40*reg24; T reg85=reg27*reg39;
    T reg86=reg62*reg21; T reg87=reg45*reg17; T reg88=reg24*reg31; reg47=reg48+reg47; T reg89=reg39*reg33;
    T reg90=reg6*reg41; T reg91=reg62*reg26; T reg92=reg57*reg19; T reg93=reg40*reg25; T reg94=reg43*reg23;
    T reg95=reg41*reg17; T reg96=reg23*reg28; T reg97=reg38*reg22; T reg98=reg61*reg20; T reg99=reg24*reg28;
    T reg100=reg38*reg33; T reg101=reg23*reg31; T reg102=reg6*reg45; reg40=reg40*reg27; reg57=reg57*reg21;
    reg60=reg58+reg60; T reg103=reg42*reg26; reg61=reg6*reg61; T reg104=reg43*reg22; T reg105=reg39*reg22;
    T reg106=reg34*reg20; reg53=reg53*reg20; reg49=reg49+reg54; reg62=reg62*reg19; reg39=reg39*reg25;
    reg36=reg36*reg19; T reg107=reg33*reg28; reg52=reg52+reg63; reg35=reg35*reg26; reg31=reg27*reg31;
    T reg108=reg42*reg21; reg34=reg6*reg34; reg28=reg27*reg28; T reg109=reg37*reg21; reg51=reg50+reg51;
    reg65=reg64+reg65; T reg110=reg59*reg8; reg68=reg61+reg68; reg107=reg107-reg80; reg96=reg96-reg90;
    reg76=reg70+reg76; reg36=reg36+reg67; reg75=reg75+reg77; reg61=reg51*reg8; reg70=reg56*reg8;
    T reg111=reg49*reg8; reg44=reg44-reg74; reg39=reg62-reg39; reg57=reg40-reg57; reg40=reg60*reg8;
    reg35=reg35-reg100; reg28=reg28+reg109; reg62=reg52*reg8; T reg112=reg8*reg47; reg69=reg69+reg97;
    reg86=reg85-reg86; reg106=reg106+reg104; reg31=reg31+reg108; reg105=reg53-reg105; reg34=reg34-reg94;
    reg72=reg72-reg103; reg81=reg82-reg81; reg89=reg91+reg89; reg53=reg65*reg8; reg93=reg92-reg93;
    reg79=reg79-reg78; reg88=reg88+reg87; reg101=reg101-reg102; reg71=reg73+reg71; reg99=reg99+reg95;
    reg83=reg84-reg83; reg66=reg98-reg66; reg73=ponderation*reg111; reg79=reg79*reg8; reg28=reg28*reg8;
    reg101=reg8*reg101; reg35=reg35*reg8; reg106=reg106*reg8; reg105=reg105*reg8; reg34=reg34*reg8;
    reg75=reg75*reg8; reg76=reg76*reg8; reg96=reg96*reg8; reg82=ponderation*reg53; reg88=reg8*reg88;
    reg71=reg71*reg8; reg89=reg89*reg8; reg31=reg31*reg8; reg72=reg72*reg8; reg81=reg81*reg8;
    reg68=reg68*reg8; reg86=reg86*reg8; reg83=reg83*reg8; reg99=reg99*reg8; reg84=ponderation*reg112;
    reg107=reg107*reg8; reg85=ponderation*reg110; reg57=reg57*reg8; reg91=ponderation*reg62; reg69=reg69*reg8;
    reg93=reg93*reg8; reg92=ponderation*reg70; reg98=ponderation*reg40; reg39=reg39*reg8; reg44=reg44*reg8;
    reg66=reg66*reg8; reg36=reg36*reg8; T reg113=ponderation*reg61; T tmp_4_2=-reg113; T tmp_0_3=ponderation*reg101;
    T tmp_3_1=ponderation*reg86; T tmp_3_0=ponderation*reg57; T tmp_2_1=ponderation*reg105; T tmp_0_1=ponderation*reg71; T tmp_4_5=-reg92;
    T tmp_2_2=ponderation*reg106; T tmp_2_3=-reg84; T tmp_5_5=ponderation*reg75; T tmp_3_5=-reg91; T tmp_4_1=ponderation*reg81;
    T tmp_4_0=ponderation*reg83; T tmp_1_1=ponderation*reg89; T tmp_3_3=ponderation*reg31; T tmp_1_3=ponderation*reg72; T tmp_0_0=ponderation*reg68;
    T tmp_2_4=-reg82; T tmp_3_2=-reg85; T tmp_2_0=ponderation*reg66; T tmp_4_4=ponderation*reg99; T tmp_1_4=ponderation*reg107;
    T tmp_0_5=ponderation*reg44; T tmp_5_2=ponderation*reg36; T tmp_5_1=ponderation*reg39; T tmp_5_3=-reg98; T tmp_5_0=ponderation*reg93;
    T tmp_1_2=ponderation*reg79; T tmp_2_5=ponderation*reg69; T tmp_5_4=-reg73; T tmp_3_4=ponderation*reg28; T tmp_1_5=ponderation*reg35;
    T tmp_1_0=ponderation*reg76; T tmp_0_4=ponderation*reg96; T tmp_4_3=ponderation*reg88; T tmp_0_2=ponderation*reg34;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=reg1*reg0; T reg5=reg3*reg1; T reg6=reg4*reg2; reg1=reg1*reg2; reg4=reg3*reg4;
    T reg7=reg1*reg2; T reg8=reg5*reg2; T reg9=reg4*reg2; T reg10=reg6*reg2; reg5=reg3*reg5;
    reg4=reg3*reg4; reg9=reg9+reg10; reg1=reg3*reg1; reg4=reg4-reg10; reg5=reg5-reg7;
    reg6=reg3*reg6; reg8=reg7+reg8; reg5=reg3*reg5; reg8=reg2*reg8; T reg11=reg3*reg4;
    T reg12=reg3*reg0; reg10=reg6+reg10; reg6=reg2*reg0; T reg13=reg7+reg1; T reg14=reg9*reg2;
    T reg15=pow(reg3,2); T reg16=reg10*reg2; T reg17=elem.pos(1)[1]-elem.pos(0)[1]; reg8=reg5-reg8; reg5=pow(reg2,2);
    T reg18=reg3*reg12; T reg19=elem.pos(1)[0]-elem.pos(0)[0]; reg14=reg11-reg14; reg13=reg2*reg13; reg11=reg6*reg2;
    T reg20=elem.pos(2)[1]-elem.pos(0)[1]; T reg21=elem.pos(2)[0]-elem.pos(0)[0]; reg5=reg15-reg5; reg15=reg20*reg19; T reg22=reg21*reg17;
    reg13=reg8-reg13; reg16=reg14-reg16; reg11=reg18-reg11; reg8=1-(*f.m).resolution; reg5=reg5/reg11;
    reg13=reg13/reg16; reg22=reg15-reg22; reg14=reg8*reg13; reg20=reg20/reg22; reg21=reg21/reg22;
    reg12=reg12/reg11; reg9=reg9/reg16; reg11=reg6/reg11; reg19=reg19/reg22; reg17=reg17/reg22;
    reg6=(*f.m).resolution*reg5; reg4=reg4/reg16; reg15=0.5*reg19; reg18=0.5*reg17; T reg23=(*f.m).resolution*reg11;
    T reg24=reg21-reg19; T reg25=0.5*reg20; T reg26=reg4*reg8; reg6=reg14+reg6; reg14=(*f.m).resolution*reg12;
    T reg27=reg8*reg9; T reg28=reg17-reg20; T reg29=0.5*reg28; T reg30=0.5*reg21; reg27=reg23+reg27;
    reg14=reg26+reg14; reg23=reg6*reg25; reg26=0.5*reg24; T reg31=reg6*reg15; T reg32=reg6*reg18;
    T reg33=reg14*reg19; T reg34=reg14*reg17; reg31=2*reg31; T reg35=reg6*reg26; T reg36=reg27*reg21;
    T reg37=2*reg32; reg23=2*reg23; T reg38=reg6*reg30; T reg39=reg6*reg29; T reg40=reg27*reg19;
    T reg41=reg31*reg30; T reg42=reg27*reg24; T reg43=reg37*reg15; T reg44=reg27*reg20; T reg45=reg20*reg34;
    T reg46=reg40*reg17; reg39=2*reg39; T reg47=reg23*reg30; T reg48=reg36*reg20; T reg49=reg14*reg21;
    T reg50=reg25*reg37; T reg51=reg14*reg24; T reg52=reg27*reg17; T reg53=2*reg38; T reg54=reg28*reg14;
    T reg55=reg33*reg21; reg35=2*reg35; T reg56=reg14*reg20; T reg57=reg26*reg31; T reg58=reg28*reg40;
    T reg59=reg15*reg31; T reg60=reg37*reg29; T reg61=reg26*reg37; T reg62=reg28*reg42; T reg63=reg34*reg17;
    T reg64=reg26*reg39; T reg65=reg26*reg23; T reg66=reg28*reg36; reg47=reg48+reg47; reg46=reg43+reg46;
    T reg67=reg51*reg24; T reg68=reg28*reg34; T reg69=reg39*reg29; reg55=reg55+reg50; T reg70=reg28*reg56;
    T reg71=reg56*reg20; T reg72=reg53*reg30; T reg73=reg52*reg24; T reg74=reg29*reg23; T reg75=reg29*reg31;
    T reg76=reg35*reg26; T reg77=reg33*reg24; T reg78=reg25*reg23; T reg79=reg49*reg21; T reg80=reg28*reg54;
    T reg81=reg37*reg30; T reg82=reg25*reg31; T reg83=reg53*reg26; T reg84=reg53*reg29; T reg85=reg24*reg44;
    T reg86=reg52*reg21; reg33=reg33*reg19; T reg87=reg37*reg18; T reg88=reg49*reg24; reg41=reg45+reg41;
    reg40=reg40*reg20; reg64=reg62+reg64; reg77=reg77-reg60; reg74=reg74-reg88; reg62=reg41*reg22;
    reg57=reg57-reg68; reg76=reg80+reg76; reg80=reg46*reg22; T reg89=reg22*reg47; reg85=reg85-reg84;
    reg70=reg70-reg83; reg40=reg40+reg81; reg82=reg82+reg86; reg65=reg65-reg66; reg78=reg78+reg79;
    reg75=reg75-reg73; reg69=reg67+reg69; reg59=reg59+reg63; reg33=reg33+reg87; reg58=reg58-reg61;
    reg71=reg71+reg72; reg67=reg55*reg22; T reg90=ponderation*reg80; reg65=reg22*reg65; reg70=reg70*reg22;
    reg59=reg59*reg22; reg33=reg33*reg22; T reg91=ponderation*reg89; reg85=reg85*reg22; reg78=reg78*reg22;
    reg40=reg40*reg22; reg82=reg82*reg22; T reg92=ponderation*reg67; reg64=reg64*reg22; reg58=reg58*reg22;
    reg77=reg77*reg22; reg69=reg69*reg22; reg74=reg74*reg22; reg71=reg71*reg22; reg57=reg57*reg22;
    reg76=reg76*reg22; reg75=reg75*reg22; T reg93=ponderation*reg62; T tmp_3_5=-reg92; T tmp_1_1=ponderation*reg69;
    T tmp_0_0=ponderation*reg76; T tmp_4_5=-reg90; T tmp_2_3=-reg91; T tmp_1_2=ponderation*reg85; T tmp_2_4=-reg93;
    T tmp_5_5=ponderation*reg33; T tmp_0_2=ponderation*reg70; T tmp_2_2=ponderation*reg71; T tmp_0_3=ponderation*reg65; T tmp_4_4=ponderation*reg59;
    T tmp_3_3=ponderation*reg78; T tmp_2_5=ponderation*reg40; T tmp_3_4=ponderation*reg82; T tmp_1_4=ponderation*reg75; T tmp_0_4=ponderation*reg57;
    T tmp_0_5=ponderation*reg58; T tmp_1_3=ponderation*reg74; T tmp_1_5=ponderation*reg77; T tmp_0_1=ponderation*reg64;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=reg1*reg0;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg2*reg1; reg1=reg1*reg4; T reg6=reg2*reg3; reg3=reg3*reg4;
    T reg7=reg2*reg5; T reg8=reg1*reg4; T reg9=reg3*reg4; T reg10=reg2*reg6; reg6=reg6*reg4;
    reg5=reg5*reg4; reg1=reg2*reg1; reg10=reg10-reg9; reg7=reg7-reg8; reg3=reg2*reg3;
    reg5=reg8+reg5; reg6=reg6+reg9; T reg11=reg2*reg10; reg9=reg3+reg9; reg3=reg8+reg1;
    T reg12=reg6*reg4; reg7=reg2*reg7; reg5=reg4*reg5; reg3=reg4*reg3; T reg13=reg9*reg4;
    reg12=reg11-reg12; reg5=reg7-reg5; reg3=reg5-reg3; reg13=reg12-reg13; reg3=reg3/reg13;
    reg10=reg10/reg13; reg6=reg6/reg13; reg5=reg6*reg3; reg7=reg10*reg3; reg13=reg9/reg13;
    reg9=reg10*reg7; reg11=(*f.m).alpha*reg6; reg12=reg6*reg5; T reg14=reg10*(*f.m).alpha; reg11=reg14+reg11;
    reg12=reg9-reg12; reg13=(*f.m).alpha*reg13; reg5=reg5/reg12; reg9=reg4*reg0; reg12=reg7/reg12;
    reg11=reg13+reg11; reg7=reg2*reg0; reg13=reg9*reg4; reg5=reg11*reg5; reg12=reg11*reg12;
    reg11=reg2*reg7; reg14=1-(*f.m).resolution; reg5=reg12-reg5; reg13=reg11-reg13; reg5=reg14*reg5;
    reg7=reg7/reg13; reg11=(*f.m).alpha*(*f.m).resolution; reg9=reg9/reg13; reg12=elem.pos(2)[1]-elem.pos(0)[1]; T reg15=elem.pos(1)[0]-elem.pos(0)[0];
    T reg16=elem.pos(1)[1]-elem.pos(0)[1]; T reg17=(*f.m).resolution*reg7; T reg18=elem.pos(2)[0]-elem.pos(0)[0]; reg6=reg14*reg6; reg5=reg11+reg5;
    reg10=reg10*reg14; reg11=(*f.m).resolution*reg9; T reg19=reg18*reg16; reg17=reg10+reg17; reg6=reg11+reg6;
    reg10=reg12*reg15; reg5=(*f.m).deltaT*reg5; reg11=reg17*reg5; reg19=reg10-reg19; reg10=reg6*reg5;
    T reg20=1-var_inter[0]; reg16=reg16/reg19; reg15=reg15/reg19; reg12=reg12/reg19; T reg21=reg10+reg11;
    reg18=reg18/reg19; reg20=reg20-var_inter[1]; T reg22=reg16-reg12; T reg23=var_inter[0]*(*f.m).f_vol[1]; T reg24=reg21*reg16;
    T reg25=reg18-reg15; T reg26=reg21*reg18; T reg27=var_inter[1]*(*f.m).f_vol[0]; T reg28=reg21*reg22; T reg29=var_inter[1]*(*f.m).f_vol[1];
    T reg30=reg20*(*f.m).f_vol[0]; T reg31=reg21*reg15; T reg32=reg26-reg23; T reg33=reg20*(*f.m).f_vol[1]; T reg34=reg21*reg25;
    T reg35=reg21*reg12; T reg36=reg24-reg27; T reg37=var_inter[0]*(*f.m).f_vol[0]; T reg38=reg29+reg31; T reg39=reg34+reg33;
    T reg40=reg35+reg37; reg32=reg32*reg19; reg36=reg36*reg19; T reg41=reg28+reg30; T reg42=reg40*reg19;
    reg36=ponderation*reg36; reg32=ponderation*reg32; T reg43=reg39*reg19; T reg44=reg41*reg19; T reg45=reg38*reg19;
    T vec_3=-reg32; reg32=ponderation*reg43; T vec_1=reg32; T reg46=ponderation*reg44; T vec_0=reg46;
    T reg47=ponderation*reg45; T vec_5=reg47; T reg48=reg42*ponderation; T vec_2=reg48; T vec_4=-reg36;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=reg1*reg0;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg1*reg2; T reg6=reg3*reg2; reg3=reg4*reg3; reg1=reg4*reg1;
    T reg7=reg1*reg2; T reg8=reg5*reg2; reg1=reg4*reg1; T reg9=reg4*reg3; T reg10=reg6*reg2;
    reg3=reg3*reg2; reg5=reg4*reg5; reg6=reg4*reg6; reg9=reg9-reg10; reg3=reg3+reg10;
    reg7=reg8+reg7; reg1=reg1-reg8; reg10=reg6+reg10; reg6=reg3*reg2; reg7=reg2*reg7;
    T reg11=reg8+reg5; reg1=reg4*reg1; T reg12=reg4*reg9; reg11=reg2*reg11; reg6=reg12-reg6;
    reg12=reg10*reg2; reg7=reg1-reg7; reg11=reg7-reg11; reg12=reg6-reg12; reg11=reg11/reg12;
    reg9=reg9/reg12; reg3=reg3/reg12; reg1=reg3*reg11; reg6=reg9*reg11; reg12=reg10/reg12;
    reg7=(*f.m).alpha*reg3; reg10=reg9*(*f.m).alpha; T reg13=reg9*reg6; T reg14=reg3*reg1; reg14=reg13-reg14;
    reg13=elem.pos(1)[0]-elem.pos(0)[0]; T reg15=elem.pos(1)[1]-elem.pos(0)[1]; reg7=reg10+reg7; reg12=(*f.m).alpha*reg12; reg10=elem.pos(2)[0]-elem.pos(0)[0];
    T reg16=elem.pos(2)[1]-elem.pos(0)[1]; T reg17=reg10*reg15; T reg18=reg16*reg13; reg7=reg12+reg7; reg1=reg1/reg14;
    reg14=reg6/reg14; reg1=reg7*reg1; reg17=reg18-reg17; reg6=reg4*reg0; reg12=reg2*reg0;
    reg14=reg7*reg14; reg7=reg4*reg6; reg18=pow(reg2,2); T reg19=pow(reg4,2); T reg20=vectors[0][indices[1]+0]-vectors[0][indices[0]+0];
    reg16=reg16/reg17; reg10=reg10/reg17; T reg21=1-(*f.m).resolution; T reg22=reg12*reg2; reg1=reg14-reg1;
    reg15=reg15/reg17; reg13=reg13/reg17; reg14=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; T reg23=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; T reg24=vectors[0][indices[1]+1]-vectors[0][indices[0]+1];
    T reg25=reg24*reg16; T reg26=reg20*reg10; T reg27=reg23*reg13; T reg28=(*f.m).alpha*(*f.m).resolution; reg1=reg21*reg1;
    T reg29=reg14*reg15; reg22=reg7-reg22; reg18=reg19-reg18; reg29=reg25-reg29; reg26=reg27-reg26;
    reg18=reg18/reg22; reg12=reg12/reg22; reg22=reg6/reg22; reg1=reg28+reg1; reg23=reg23*reg15;
    reg14=reg14*reg13; reg24=reg24*reg10; reg20=reg20*reg16; reg11=reg21*reg11; reg3=reg21*reg3;
    reg6=(*f.m).resolution*reg22; reg7=(*f.m).resolution*reg12; reg24=reg14-reg24; reg29=reg26+reg29; reg21=reg9*reg21;
    reg9=(*f.m).resolution*reg18; reg23=reg20-reg23; reg1=(*f.m).deltaT*reg1; reg23=reg23-reg1; reg9=reg11+reg9;
    reg24=reg24-reg1; reg6=reg21+reg6; reg3=reg7+reg3; reg29=0.5*reg29; reg7=reg24*reg6;
    reg11=reg15-reg16; reg14=reg23*reg3; reg23=reg23*reg6; reg19=reg10-reg13; reg29=reg9*reg29;
    reg24=reg24*reg3; reg7=reg14+reg7; reg14=0.5*reg15; reg20=0.5*reg16; reg21=0.5*reg10;
    reg24=reg23+reg24; reg29=2*reg29; reg23=0.5*reg19; reg25=0.5*reg11; reg26=1-var_inter[0];
    reg27=0.5*reg13; reg28=reg11*reg24; T reg30=reg29*reg21; T reg31=reg24*reg16; T reg32=reg23*reg29;
    reg26=reg26-var_inter[1]; T reg33=reg7*reg13; T reg34=reg7*reg10; T reg35=reg7*reg19; T reg36=reg29*reg25;
    T reg37=reg24*reg15; T reg38=reg20*reg29; T reg39=reg29*reg27; T reg40=reg29*reg14; T reg41=reg26*(*f.m).f_vol[0];
    reg32=reg28+reg32; reg28=var_inter[0]*(*f.m).f_vol[0]; reg31=reg31-reg30; reg36=reg35+reg36; reg35=var_inter[1]*(*f.m).f_vol[1];
    reg33=reg33-reg40; T reg42=reg26*(*f.m).f_vol[1]; reg39=reg39-reg37; T reg43=var_inter[1]*(*f.m).f_vol[0]; T reg44=var_inter[0]*(*f.m).f_vol[1];
    reg38=reg38-reg34; reg32=reg32-reg41; reg33=reg33-reg35; reg39=reg39-reg43; reg36=reg36-reg42;
    reg38=reg38-reg44; reg31=reg31-reg28; reg36=reg36*reg17; reg31=reg17*reg31; reg39=reg39*reg17;
    reg33=reg33*reg17; reg32=reg32*reg17; reg38=reg38*reg17; T vec_3=ponderation*reg38; T vec_0=ponderation*reg32;
    T vec_2=ponderation*reg31; T vec_5=ponderation*reg33; T vec_1=ponderation*reg36; T vec_4=ponderation*reg39;
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
    node.dep[0]=vecs[0][indice+0]; node.dep[1]=vecs[0][indice+1];
  }
  template<class TE,class TTs,class Tvec>
  inline static void set_nodal_unknowns(TE &node,const TTs &f,const Tvec &vecs,unsigned indice,T partial_ts) {
    node.dep[0]=vecs[0][indice+0]; node.dep[1]=vecs[0][indice+1];
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=reg0*reg1;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg2*reg1; reg1=reg4*reg1; T reg6=reg4*reg3; reg3=reg2*reg3;
    T reg7=reg4*reg5; T reg8=reg4*reg6; T reg9=reg4*reg3; T reg10=reg4*reg1; reg5=reg2*reg5;
    reg3=reg2*reg3; reg1=reg2*reg1; reg5=reg5-reg10; reg6=reg2*reg6; reg7=reg10+reg7;
    reg3=reg3-reg8; reg9=reg8+reg9; T reg11=0.5*elem.pos(1)[0]; T reg12=reg4*reg9; T reg13=0.5*elem.pos(0)[1];
    reg6=reg8+reg6; reg8=0.5*elem.pos(0)[0]; T reg14=0.5*elem.pos(1)[1]; reg5=reg2*reg5; reg7=reg4*reg7;
    T reg15=reg2*reg3; T reg16=reg10+reg1; T reg17=reg13+reg14; T reg18=0.5*elem.pos(2)[1]; T reg19=reg11-reg8;
    T reg20=0.5*elem.pos(2)[0]; reg8=reg11+reg8; reg13=reg14-reg13; reg12=reg15-reg12; reg11=reg4*reg6;
    reg16=reg4*reg16; reg7=reg5-reg7; reg8=reg20-reg8; reg13=reg13+reg18; reg5=0.5*elem.pos(3)[0];
    reg11=reg12-reg11; reg17=reg18-reg17; reg12=0.5*elem.pos(3)[1]; reg20=reg19+reg20; reg16=reg7-reg16;
    reg8=reg5+reg8; reg13=reg13-reg12; reg5=reg20-reg5; reg7=0.5*vectors[0][indices[1]+0]; reg14=0.5*vectors[0][indices[0]+0];
    reg3=reg3/reg11; reg9=reg9/reg11; reg15=0.5*vectors[0][indices[0]+1]; reg18=0.5*vectors[0][indices[1]+1]; reg16=reg16/reg11;
    reg12=reg17+reg12; reg17=0.78867513459481286553*elem.pos(0)[0]; reg19=reg5*reg12; reg20=0.21132486540518713447*elem.pos(0)[1]; T reg21=0.78867513459481286553*elem.pos(1)[1];
    T reg22=0.21132486540518713447*elem.pos(0)[0]; T reg23=0.21132486540518713447*elem.pos(1)[0]; T reg24=reg18-reg15; T reg25=reg13*reg8; T reg26=0.78867513459481286553*elem.pos(1)[0];
    T reg27=0.78867513459481286553*elem.pos(0)[1]; reg18=reg15+reg18; reg15=0.5*vectors[0][indices[2]+1]; T reg28=0.5*vectors[0][indices[2]+0]; T reg29=reg7-reg14;
    T reg30=reg3*reg16; T reg31=0.21132486540518713447*elem.pos(1)[1]; T reg32=reg9*reg16; reg14=reg7+reg14; reg24=reg15+reg24;
    reg7=reg3*reg30; T reg33=reg9*reg32; T reg34=(*f.m).alpha*reg3; T reg35=(*f.m).alpha*reg9; T reg36=reg23+reg17;
    T reg37=reg31-reg20; T reg38=0.78867513459481286553*elem.pos(2)[0]; reg23=reg23-reg22; T reg39=0.78867513459481286553*elem.pos(2)[1]; reg22=reg22+reg26;
    reg20=reg20+reg21; reg29=reg28+reg29; T reg40=0.5*vectors[0][indices[3]+0]; T reg41=0.21132486540518713447*elem.pos(2)[1]; reg31=reg27+reg31;
    reg25=reg19-reg25; reg11=reg6/reg11; reg14=reg28-reg14; reg6=0.5*vectors[0][indices[3]+1]; reg19=0.21132486540518713447*elem.pos(2)[0];
    reg18=reg15-reg18; reg22=reg38-reg22; reg15=0.78867513459481286553*elem.pos(3)[1]; reg12=reg12/reg25; reg8=reg8/reg25;
    reg38=reg23+reg38; reg23=0.21132486540518713447*elem.pos(3)[0]; reg37=reg39+reg37; reg28=reg2*reg0; T reg42=0.78867513459481286553*elem.pos(3)[0];
    T reg43=reg4*reg0; reg36=reg19-reg36; reg11=(*f.m).alpha*reg11; reg35=reg34+reg35; reg29=reg29-reg40;
    reg13=reg13/reg25; reg33=reg7-reg33; reg14=reg40+reg14; reg24=reg24-reg6; reg18=reg6+reg18;
    reg25=reg5/reg25; reg5=0.21132486540518713447*elem.pos(3)[1]; reg20=reg39-reg20; reg27=reg21-reg27; reg31=reg41-reg31;
    reg17=reg26-reg17; reg6=0.78867513459481286553*PNODE(0).dep[1]; reg7=0.21132486540518713447*PNODE(1).dep[1]; reg30=reg30/reg33; reg33=reg32/reg33;
    reg11=reg35+reg11; reg21=0.78867513459481286553*PNODE(1).dep[0]; reg26=reg8*reg24; reg32=0.21132486540518713447*PNODE(0).dep[0]; reg34=0.21132486540518713447*PNODE(1).dep[0];
    reg17=reg19+reg17; reg22=reg22+reg23; reg37=reg37-reg15; reg19=reg4*reg43; reg27=reg41+reg27;
    reg35=reg2*reg28; reg39=reg18*reg25; reg40=reg13*reg14; reg36=reg42+reg36; reg41=reg12*reg29;
    reg15=reg31+reg15; reg31=0.78867513459481286553*PNODE(0).dep[0]; reg42=reg38-reg42; reg38=0.78867513459481286553*PNODE(1).dep[1]; reg20=reg20+reg5;
    T reg44=0.21132486540518713447*PNODE(0).dep[1]; T reg45=(*f.m).deltaT*(*f.m).alpha; reg19=reg35-reg19; reg35=reg20*reg42; T reg46=reg37*reg22;
    reg23=reg17-reg23; reg17=reg15*reg42; T reg47=reg36*reg37; reg5=reg27-reg5; reg27=0.78867513459481286553*PNODE(2).dep[0];
    reg40=reg41-reg40; elem.epsilon[0][0]=reg40; reg41=0.78867513459481286553*PNODE(2).dep[1]; T reg48=0.21132486540518713447*PNODE(2).dep[1]; T reg49=reg7+reg6;
    reg7=reg7-reg44; T reg50=reg31+reg34; reg30=reg30*reg11; T reg51=reg32+reg21; T reg52=0.21132486540518713447*PNODE(2).dep[0];
    reg11=reg33*reg11; reg26=reg39-reg26; elem.epsilon[0][1]=reg26; reg32=reg34-reg32; reg44=reg44+reg38;
    reg33=1-(*f.m).resolution; reg11=reg30-reg11; reg30=reg26-reg45; reg34=reg40-reg45; reg39=reg23*reg20;
    reg7=reg7+reg41; T reg53=0.78867513459481286553*PNODE(3).dep[1]; T reg54=reg5*reg22; reg47=reg17-reg47; reg44=reg41-reg44;
    reg50=reg52-reg50; reg17=0.21132486540518713447*PNODE(3).dep[1]; reg49=reg48-reg49; reg41=0.21132486540518713447*PNODE(3).dep[0]; reg51=reg27-reg51;
    T reg55=0.78867513459481286553*PNODE(3).dep[0]; reg32=reg27+reg32; reg46=reg35-reg46; reg31=reg21-reg31; reg6=reg38-reg6;
    reg28=reg28/reg19; reg43=reg43/reg19; reg49=reg53+reg49; reg21=reg15*reg23; reg52=reg31+reg52;
    reg53=reg7-reg53; reg7=reg22/reg46; reg50=reg50+reg55; reg27=reg34*reg28; reg54=reg39-reg54;
    reg31=reg30*reg43; reg35=pow(reg2,2); reg38=reg37/reg47; reg34=reg34*reg43; reg30=reg30*reg28;
    reg39=reg15/reg47; reg48=reg6+reg48; reg6=reg5*reg36; T reg56=reg20/reg46; reg11=reg33*reg11;
    T reg57=pow(reg4,2); reg55=reg32-reg55; reg32=(*f.m).resolution*(*f.m).alpha; reg37=reg37/reg46; reg51=reg51+reg41;
    T reg58=reg36/reg47; T reg59=reg42/reg46; reg44=reg17+reg44; reg42=reg42/reg47; T reg60=reg44*reg59;
    T reg61=reg53*reg39; T reg62=reg38*reg49; T reg63=reg37*reg44; T reg64=reg53*reg56; T reg65=reg7*reg55;
    reg59=reg51*reg59; reg7=reg53*reg7; reg57=reg35-reg57; reg37=reg37*reg51; reg56=reg56*reg55;
    reg11=reg32+reg11; reg6=reg21-reg6; reg22=reg22/reg54; reg17=reg48-reg17; reg21=reg23/reg54;
    reg53=reg53*reg58; reg32=reg49*reg42; reg30=reg34+reg30; reg34=reg5/reg54; reg58=reg58*reg55;
    reg42=reg50*reg42; reg31=reg27+reg31; reg41=reg52-reg41; reg20=reg20/reg54; reg38=reg38*reg50;
    reg55=reg39*reg55; reg65=reg59-reg65; reg27=(*f.m).resolution*reg28; reg53=reg32-reg53; reg58=reg42-reg58;
    reg32=reg34*reg51; reg35=(*f.m).resolution*reg43; reg39=reg4*reg30; reg42=reg20*reg17; reg3=reg3*reg33;
    reg48=reg41*reg22; reg9=reg9*reg33; reg5=reg5/reg6; reg52=reg2*reg31; reg7=reg60-reg7;
    reg23=reg23/reg6; reg19=reg57/reg19; reg62=reg61-reg62; reg11=(*f.m).deltaT*reg11; reg15=reg15/reg6;
    reg51=reg21*reg51; reg20=reg20*reg41; reg38=reg55-reg38; reg22=reg17*reg22; reg36=reg36/reg6;
    reg37=reg56-reg37; reg34=reg34*reg44; reg44=reg21*reg44; reg63=reg64-reg63; reg31=reg4*reg31;
    reg30=reg2*reg30; reg32=reg20-reg32; reg20=(*f.m).resolution*reg19; reg21=reg5*reg49; reg55=reg15*reg17;
    reg17=reg36*reg17; reg36=reg36*reg41; reg49=reg23*reg49; reg23=reg23*reg50; reg58=reg62+reg58;
    reg34=reg42-reg34; reg41=reg15*reg41; reg9=reg35+reg9; reg37=reg37-reg11; reg27=reg3+reg27;
    reg63=reg65+reg63; reg30=reg30-reg31; reg52=reg52-reg39; reg53=reg53-reg11; reg48=reg51-reg48;
    reg38=reg38-reg11; reg50=reg5*reg50; reg22=reg44-reg22; reg7=reg7-reg11; reg16=reg33*reg16;
    reg32=reg32-reg11; reg34=reg48+reg34; reg16=reg20+reg16; reg29=reg8*reg29; reg50=reg41-reg50;
    reg25=reg14*reg25; reg3=reg7*reg27; reg5=reg9*reg37; reg8=reg53*reg9; reg14=reg7*reg9;
    reg15=reg38*reg9; reg20=reg53*reg27; reg52=reg45+reg52; reg63=0.5*reg63; reg33=reg27*reg37;
    reg18=reg13*reg18; reg24=reg12*reg24; reg30=reg45+reg30; reg17=reg49-reg17; reg36=reg23-reg36;
    reg31=reg39+reg31; reg21=reg55-reg21; reg12=reg38*reg27; reg58=0.5*reg58; reg22=reg22-reg11;
    reg13=reg52+reg30; reg23=reg63*reg16; reg31=reg45-reg31; reg35=reg58*reg16; reg8=reg12+reg8;
    reg34=0.5*reg34; reg12=reg32*reg9; reg39=reg22*reg27; reg41=reg32*reg27; reg21=reg36+reg21;
    reg29=reg25-reg29; reg25=reg22*reg9; reg17=reg17-reg11; reg18=reg24-reg18; reg20=reg15+reg20;
    reg50=reg50-reg11; reg33=reg14+reg33; reg3=reg5+reg3; reg20=reg53*reg20; reg41=reg25+reg41;
    reg39=reg12+reg39; reg8=reg38*reg8; reg5=reg34*reg16; reg35=2*reg35; reg13=reg13+reg31;
    reg23=2*reg23; reg21=0.5*reg21; reg3=reg7*reg3; reg37=reg33*reg37; reg18=reg29+reg18;
    reg7=reg50*reg27; reg12=reg17*reg27; reg14=reg17*reg9; reg15=reg50*reg9; reg22=reg39*reg22;
    reg24=reg21*reg16; reg32=reg41*reg32; reg12=reg15+reg12; reg14=reg7+reg14; reg20=reg8+reg20;
    reg18=0.5*reg18; elem.epsilon[0][2]=reg18; reg3=reg37+reg3; reg23=reg63*reg23; reg35=reg58*reg35;
    reg13=reg13/3; reg5=2*reg5; reg7=reg18*reg19; reg5=reg34*reg5; reg23=reg3+reg23;
    reg20=reg35+reg20; reg52=reg52-reg13; reg22=reg32+reg22; reg14=reg50*reg14; reg24=2*reg24;
    reg12=reg17*reg12; reg30=reg30-reg13; reg47=reg20*reg47; reg46=reg23*reg46; reg5=reg22+reg5;
    reg30=pow(reg30,2); reg12=reg14+reg12; reg13=reg31-reg13; reg24=reg21*reg24; reg7=reg0*reg7;
    reg52=pow(reg52,2); reg47=0.25*reg47; reg3=2*reg7; reg46=0.25*reg46; reg13=pow(reg13,2);
    reg54=reg5*reg54; reg24=reg12+reg24; reg30=reg52+reg30; reg3=reg7*reg3; reg54=0.25*reg54;
    reg24=reg6*reg24; reg47=reg46+reg47; reg13=reg30+reg13; reg47=reg54+reg47; reg40=reg40-reg11;
    reg26=reg26-reg11; reg24=0.25*reg24; reg3=reg13+reg3; reg5=reg40*reg27; reg6=reg26*reg9;
    reg40=reg40*reg9; reg26=reg26*reg27; reg24=reg47+reg24; reg3=1.5*reg3; elem.sigma[0][2]=reg18*reg16;
    elem.sigma[0][0]=reg6+reg5; elem.sigma[0][1]=reg40+reg26; elem.sigma_von_mises=pow(reg3,0.5); elem.ener=reg24/2;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg3*reg2; reg2=reg4*reg2; T reg6=reg3*reg1; reg1=reg4*reg1;
    T reg7=reg3*reg5; T reg8=reg4*reg2; reg5=reg4*reg5; T reg9=reg4*reg6; T reg10=reg4*reg1;
    reg6=reg3*reg6; reg6=reg6-reg10; reg1=reg3*reg1; reg9=reg10+reg9; reg5=reg8+reg5;
    reg7=reg7-reg8; reg2=reg3*reg2; T reg11=reg10+reg1; reg9=reg4*reg9; reg6=reg3*reg6;
    reg2=reg8+reg2; reg8=reg4*reg5; T reg12=reg3*reg7; reg8=reg12-reg8; reg12=reg4*reg2;
    reg11=reg4*reg11; reg9=reg6-reg9; reg6=1-var_inter[1]; T reg13=1-var_inter[0]; reg12=reg8-reg12;
    reg11=reg9-reg11; reg8=elem.pos(1)[0]*var_inter[0]; reg9=reg13*elem.pos(0)[0]; reg11=reg11/reg12; reg7=reg7/reg12;
    T reg14=reg6*elem.pos(1)[1]; reg5=reg5/reg12; T reg15=reg6*elem.pos(0)[1]; T reg16=reg6*elem.pos(1)[0]; T reg17=reg6*elem.pos(0)[0];
    T reg18=elem.pos(0)[1]*reg13; T reg19=elem.pos(1)[1]*var_inter[0]; reg14=reg14-reg15; T reg20=var_inter[1]*elem.pos(2)[1]; T reg21=elem.pos(2)[1]*var_inter[0];
    T reg22=reg18+reg19; reg16=reg16-reg17; T reg23=var_inter[1]*elem.pos(2)[0]; T reg24=reg7*reg11; T reg25=reg9+reg8;
    T reg26=reg5*reg11; T reg27=elem.pos(2)[0]*var_inter[0]; T reg28=(*f.m).alpha*reg7; T reg29=reg5*reg26; T reg30=var_inter[1]*elem.pos(3)[1];
    T reg31=(*f.m).alpha*reg5; reg20=reg14+reg20; reg23=reg16+reg23; reg14=var_inter[1]*elem.pos(3)[0]; reg27=reg27-reg25;
    reg12=reg2/reg12; reg2=reg13*elem.pos(3)[0]; reg16=reg13*elem.pos(3)[1]; reg21=reg21-reg22; T reg32=reg7*reg24;
    reg16=reg21+reg16; reg2=reg27+reg2; reg31=reg28+reg31; reg20=reg20-reg30; reg12=(*f.m).alpha*reg12;
    reg23=reg23-reg14; reg29=reg32-reg29; reg26=reg26/reg29; reg12=reg31+reg12; reg21=reg4*reg0;
    reg29=reg24/reg29; reg24=reg3*reg0; reg27=reg20*reg2; reg28=reg23*reg16; reg29=reg29*reg12;
    reg12=reg26*reg12; reg26=pow(reg3,2); reg31=reg3*reg24; reg32=pow(reg4,2); reg27=reg28-reg27;
    reg28=reg4*reg21; T reg33=1-(*f.m).resolution; reg12=reg29-reg12; reg28=reg31-reg28; reg16=reg16/reg27;
    reg23=reg23/reg27; reg2=reg2/reg27; reg20=reg20/reg27; reg32=reg26-reg32; reg26=reg6*reg16;
    reg29=var_inter[1]*reg16; reg31=reg13*reg20; reg12=reg33*reg12; T reg34=(*f.m).resolution*(*f.m).alpha; T reg35=reg6*reg2;
    T reg36=reg13*reg23; T reg37=var_inter[0]*reg23; T reg38=var_inter[1]*reg2; reg32=reg32/reg28; T reg39=var_inter[0]*reg20;
    reg24=reg24/reg28; reg28=reg21/reg28; reg21=(*f.m).resolution*reg32; T reg40=reg26+reg39; T reg41=reg35+reg37;
    T reg42=reg36+reg38; reg12=reg34+reg12; reg34=(*f.m).resolution*reg28; reg5=reg5*reg33; reg11=reg33*reg11;
    T reg43=(*f.m).resolution*reg24; T reg44=reg31+reg29; reg33=reg7*reg33; reg7=0.5*reg44; T reg45=0.5*reg41;
    T reg46=reg29-reg39; T reg47=0.5*reg40; T reg48=reg37-reg38; reg43=reg33+reg43; reg5=reg34+reg5;
    reg11=reg21+reg11; reg12=(*f.m).deltaT*reg12; reg21=0.5*reg42; reg33=reg35-reg36; reg34=reg31-reg26;
    T reg49=0.5*reg33; T reg50=0.5*reg46; T reg51=reg7*reg11; T reg52=0.5*reg34; T reg53=reg47*reg11;
    T reg54=reg45*reg11; T reg55=reg21*reg11; T reg56=0.5*reg48; T reg57=reg43*reg12; T reg58=reg5*reg12;
    T reg59=2*reg54; reg55=2*reg55; T reg60=reg44*reg5; T reg61=reg44*reg43; T reg62=reg49*reg11;
    T reg63=2*reg51; T reg64=var_inter[1]*reg13; T reg65=reg6*var_inter[0]; T reg66=reg52*reg11; T reg67=reg50*reg11;
    T reg68=reg56*reg11; T reg69=reg42*reg5; T reg70=reg40*reg5; T reg71=reg57+reg58; T reg72=reg41*reg5;
    reg53=2*reg53; T reg73=reg42*reg43; T reg74=reg41*reg43; T reg75=reg40*reg43; T reg76=reg34*reg43;
    reg66=2*reg66; T reg77=reg33*reg5; T reg78=reg64*(*f.m).f_vol[0]; T reg79=reg21*reg59; T reg80=reg44*reg75;
    T reg81=reg21*reg63; T reg82=reg65*(*f.m).f_vol[1]; reg62=2*reg62; T reg83=reg44*reg69; T reg84=reg41*reg71;
    T reg85=reg46*reg5; T reg86=reg6*reg13; T reg87=reg7*reg55; T reg88=reg60*reg42; T reg89=reg59*reg47;
    T reg90=var_inter[1]*var_inter[0]; T reg91=reg33*reg43; T reg92=reg42*reg74; T reg93=reg34*reg5; T reg94=reg7*reg53;
    T reg95=reg40*reg72; T reg96=reg45*reg53; T reg97=reg48*reg43; T reg98=reg41*reg70; T reg99=reg48*reg5;
    T reg100=reg44*reg71; reg67=2*reg67; T reg101=reg73*reg41; T reg102=reg47*reg63; T reg103=reg46*reg43;
    T reg104=reg40*reg61; reg68=2*reg68; T reg105=reg45*reg55; T reg106=reg62*reg47; T reg107=reg73*reg42;
    reg87=reg88+reg87; T reg108=reg41*reg93; T reg109=reg7*reg67; reg83=reg81+reg83; T reg110=reg66*reg47;
    T reg111=reg97*reg42; reg98=reg89+reg98; T reg112=reg41*reg91; T reg113=reg42*reg70; T reg114=reg56*reg55;
    T reg115=reg46*reg61; T reg116=reg56*reg63; T reg117=reg46*reg69; T reg118=reg62*reg50; T reg119=reg48*reg93;
    T reg120=reg66*reg50; T reg121=reg48*reg91; T reg122=reg59*reg50; T reg123=reg48*reg70; T reg124=reg21*reg53;
    T reg125=reg44*reg72; T reg126=reg21*reg68; T reg127=reg44*reg103; T reg128=reg21*reg67; T reg129=reg21*reg55;
    T reg130=reg44*reg61; T reg131=reg33*reg71; T reg132=reg34*reg71; T reg133=reg7*reg63; reg70=reg33*reg70;
    T reg134=reg52*reg59; T reg135=reg33*reg91; T reg136=reg52*reg66; T reg137=reg33*reg93; T reg138=reg52*reg62;
    T reg139=reg34*reg69; T reg140=reg49*reg63; T reg141=reg34*reg61; T reg142=reg49*reg55; T reg143=reg34*reg99;
    T reg144=reg49*reg67; T reg145=reg45*reg67; T reg146=reg34*reg103; T reg147=reg49*reg68; T reg148=reg34*reg72;
    T reg149=reg49*reg53; T reg150=reg34*reg75; T reg151=reg49*reg59; T reg152=reg45*reg63; reg69=reg40*reg69;
    reg105=reg104+reg105; T reg153=reg44*reg99; T reg154=reg40*reg99; T reg155=reg45*reg68; T reg156=reg40*reg103;
    T reg157=reg7*reg68; T reg158=reg42*reg85; reg94=reg92+reg94; T reg159=reg7*reg59; T reg160=reg7*reg66;
    reg91=reg42*reg91; T reg161=reg7*reg62; reg93=reg42*reg93; T reg162=reg33*reg85; T reg163=reg52*reg68;
    T reg164=reg33*reg74; T reg165=reg52*reg53; T reg166=reg21*reg66; T reg167=reg44*reg77; reg101=reg102+reg101;
    reg80=reg79+reg80; T reg168=reg86*(*f.m).f_vol[0]; T reg169=reg90*(*f.m).f_vol[1]; T reg170=reg64*(*f.m).f_vol[1]; T reg171=reg49*reg62;
    T reg172=reg34*reg76; T reg173=reg49*reg66; T reg174=reg34*reg77; T reg175=reg60*reg41; T reg176=reg47*reg55;
    T reg177=reg100-reg78; T reg178=reg52*reg55; T reg179=reg86*(*f.m).f_vol[1]; T reg180=reg65*(*f.m).f_vol[0]; T reg181=reg60*reg33;
    T reg182=reg90*(*f.m).f_vol[0]; T reg183=reg97*reg41; T reg184=reg50*reg68; T reg185=reg48*reg85; T reg186=reg48*reg74;
    T reg187=reg50*reg53; T reg188=reg46*reg75; T reg189=reg50*reg67; T reg190=reg97*reg48; T reg191=reg59*reg56;
    reg55=reg50*reg55; T reg192=reg60*reg48; T reg193=reg77*reg46; T reg194=reg66*reg56; T reg195=reg52*reg67;
    T reg196=reg50*reg63; T reg197=reg73*reg48; reg97=reg97*reg33; T reg198=reg76*reg46; T reg199=reg21*reg62;
    T reg200=reg44*reg76; T reg201=reg62*reg56; T reg202=reg42*reg71; T reg203=reg59*reg45; T reg204=reg47*reg53;
    T reg205=reg41*reg74; T reg206=reg40*reg71; reg62=reg62*reg45; reg73=reg73*reg33; T reg207=reg52*reg63;
    T reg208=reg46*reg71; reg96=reg95+reg96; T reg209=reg84-reg82; T reg210=reg56*reg67; T reg211=reg46*reg72;
    reg75=reg40*reg75; T reg212=reg47*reg68; reg85=reg41*reg85; reg76=reg76*reg40; reg53=reg56*reg53;
    reg68=reg56*reg68; reg99=reg46*reg99; reg77=reg77*reg40; T reg213=reg48*reg71; reg103=reg46*reg103;
    reg66=reg66*reg45; reg67=reg47*reg67; reg70=reg70-reg134; reg75=reg75+reg203; reg198=reg201+reg198;
    reg139=reg139-reg140; reg62=reg76-reg62; reg193=reg194+reg193; reg157=reg158-reg157; reg137=reg138+reg137;
    reg135=reg136+reg135; reg66=reg77-reg66; reg187=reg187-reg186; reg146=reg147+reg146; reg178=reg178-reg181;
    reg183=reg67-reg183; reg160=reg91-reg160; reg161=reg93-reg161; reg85=reg212-reg85; reg149=reg149-reg148;
    reg162=reg163+reg162; reg143=reg144+reg143; reg176=reg176+reg175; reg113=reg159+reg113; reg97=reg195+reg97;
    reg73=reg73-reg207; reg204=reg204+reg205; reg150=reg150-reg151; reg67=reg101*reg27; reg142=reg142-reg141;
    reg165=reg165-reg164; reg76=reg27*reg94; reg77=reg96*reg27; reg177=reg177*reg27; reg69=reg69+reg152;
    reg174=reg173+reg174; reg103=reg68+reg103; reg68=reg83*reg27; reg108=reg106-reg108; reg119=reg118+reg119;
    reg129=reg129+reg130; reg188=reg188-reg191; reg91=reg213+reg169; reg112=reg110-reg112; reg107=reg107+reg133;
    reg121=reg120+reg121; reg127=reg126-reg127; reg93=reg27*reg98; reg106=reg208+reg182; reg124=reg124+reg125;
    reg53=reg53-reg211; reg109=reg111-reg109; reg110=reg80*reg27; reg111=reg27*reg87; reg123=reg123-reg122;
    reg209=reg209*reg27; reg153=reg128-reg153; reg114=reg114-reg115; reg190=reg189+reg190; reg55=reg55-reg192;
    reg99=reg210+reg99; reg145=reg154-reg145; reg197=reg197-reg196; reg155=reg156-reg155; reg118=reg202+reg170;
    reg200=reg199-reg200; reg185=reg184+reg185; reg120=reg179+reg131; reg126=reg27*reg105; reg167=reg166-reg167;
    reg128=reg168+reg132; reg117=reg117-reg116; reg172=reg171+reg172; reg136=reg206+reg180; reg165=reg27*reg165;
    reg138=reg27*reg120; reg144=ponderation*reg111; reg73=reg73*reg27; reg135=reg27*reg135; reg62=reg62*reg27;
    reg107=reg27*reg107; reg70=reg27*reg70; reg147=reg27*reg128; reg154=reg136*reg27; reg156=ponderation*reg77;
    reg162=reg27*reg162; reg155=reg27*reg155; reg157=reg27*reg157; reg158=reg118*reg27; reg153=reg27*reg153;
    reg145=reg27*reg145; reg163=ponderation*reg76; reg166=ponderation*reg126; reg177=ponderation*reg177; reg97=reg97*reg27;
    reg113=reg27*reg113; reg69=reg27*reg69; reg160=reg27*reg160; reg108=reg27*reg108; reg171=reg91*reg27;
    reg112=reg27*reg112; reg161=reg27*reg161; reg178=reg178*reg27; reg173=reg106*reg27; reg109=reg27*reg109;
    reg209=ponderation*reg209; reg146=reg27*reg146; reg183=reg183*reg27; reg53=reg53*reg27; reg123=reg123*reg27;
    reg149=reg27*reg149; reg176=reg176*reg27; reg184=ponderation*reg110; reg150=reg27*reg150; reg189=ponderation*reg67;
    reg124=reg124*reg27; reg198=reg198*reg27; reg193=reg193*reg27; reg127=reg127*reg27; reg187=reg187*reg27;
    reg188=reg188*reg27; reg185=reg185*reg27; reg129=reg129*reg27; reg190=reg190*reg27; reg55=reg55*reg27;
    reg197=reg197*reg27; reg194=ponderation*reg68; reg200=reg200*reg27; reg174=reg174*reg27; reg167=reg167*reg27;
    reg172=reg172*reg27; reg99=reg99*reg27; reg137=reg27*reg137; reg66=reg66*reg27; reg114=reg114*reg27;
    reg139=reg27*reg139; reg75=reg75*reg27; reg195=ponderation*reg93; reg117=reg117*reg27; reg142=reg27*reg142;
    reg103=reg103*reg27; reg204=reg204*reg27; reg121=reg121*reg27; reg85=reg85*reg27; reg143=reg27*reg143;
    reg119=reg119*reg27; T tmp_4_3=ponderation*reg53; T tmp_7_6=-reg144; T tmp_2_7=ponderation*reg69; T tmp_0_0=ponderation*reg172;
    T tmp_4_5=ponderation*reg99; T tmp_5_1=ponderation*reg121; reg53=ponderation*reg171; T vec_5=reg53; reg69=ponderation*reg147;
    T vec_0=reg69; T tmp_0_1=ponderation*reg174; T tmp_4_6=ponderation*reg114; T tmp_6_6=ponderation*reg129; T tmp_3_0=ponderation*reg108;
    reg99=ponderation*reg154; T vec_2=reg99; T tmp_6_5=ponderation*reg153; T tmp_4_4=ponderation*reg103; reg103=ponderation*reg138;
    T vec_1=reg103; T tmp_6_4=ponderation*reg127; T tmp_3_1=ponderation*reg112; T vec_3=-reg209; T tmp_5_0=ponderation*reg119;
    reg108=ponderation*reg173; T vec_4=reg108; T tmp_4_7=ponderation*reg117; T tmp_4_2=ponderation*reg188; T tmp_6_3=ponderation*reg124;
    T tmp_7_5=ponderation*reg109; T tmp_7_7=ponderation*reg107; T tmp_5_2=ponderation*reg123; T tmp_6_2=-reg184; T tmp_3_6=ponderation*reg176;
    T tmp_7_1=ponderation*reg160; T tmp_0_3=ponderation*reg149; T tmp_3_5=ponderation*reg183; T tmp_7_0=ponderation*reg161; T tmp_0_4=ponderation*reg146;
    T tmp_3_4=ponderation*reg85; T tmp_1_6=ponderation*reg178; T tmp_6_7=-reg194; T tmp_3_3=ponderation*reg204; T tmp_0_5=ponderation*reg143;
    T tmp_1_3=ponderation*reg165; T tmp_3_2=-reg195; T tmp_0_6=ponderation*reg142; T tmp_1_7=ponderation*reg73; T tmp_1_2=ponderation*reg70;
    T tmp_2_2=ponderation*reg75; T tmp_0_7=ponderation*reg139; T tmp_2_1=ponderation*reg66; T tmp_2_0=ponderation*reg62; T tmp_1_1=ponderation*reg135;
    T tmp_1_0=ponderation*reg137; T tmp_6_1=ponderation*reg167; T tmp_2_6=-reg166; T tmp_6_0=ponderation*reg200; T vec_6=-reg177;
    T tmp_5_7=ponderation*reg197; T tmp_2_5=ponderation*reg145; T tmp_5_6=ponderation*reg55; T tmp_5_5=ponderation*reg190; reg55=ponderation*reg158;
    T vec_7=reg55; T tmp_2_4=ponderation*reg155; T tmp_5_4=ponderation*reg185; T tmp_5_3=ponderation*reg187; T tmp_1_4=ponderation*reg162;
    T tmp_2_3=-reg156; T tmp_7_4=ponderation*reg157; T tmp_4_1=ponderation*reg193; T tmp_4_0=ponderation*reg198; T tmp_7_3=-reg163;
    T tmp_3_7=-reg189; T tmp_0_2=ponderation*reg150; T tmp_1_5=ponderation*reg97; T tmp_7_2=ponderation*reg113;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg3*reg4; reg4=reg2*reg4; T reg6=reg2*reg1; reg1=reg3*reg1;
    T reg7=reg3*reg5; T reg8=reg2*reg4; T reg9=reg3*reg6; T reg10=reg3*reg1; reg6=reg2*reg6;
    reg4=reg3*reg4; reg1=reg2*reg1; reg6=reg6-reg10; reg9=reg10+reg9; reg8=reg8-reg7;
    reg4=reg7+reg4; reg5=reg2*reg5; reg5=reg7+reg5; reg7=reg10+reg1; T reg11=reg3*reg4;
    reg9=reg3*reg9; reg6=reg2*reg6; T reg12=reg2*reg8; T reg13=reg3*reg5; reg7=reg3*reg7;
    reg11=reg12-reg11; reg9=reg6-reg9; reg6=1-var_inter[1]; reg7=reg9-reg7; reg9=1-var_inter[0];
    reg13=reg11-reg13; reg11=reg6*elem.pos(0)[1]; reg12=reg6*elem.pos(1)[1]; T reg14=reg9*elem.pos(0)[0]; T reg15=elem.pos(1)[0]*var_inter[0];
    T reg16=elem.pos(0)[1]*reg9; T reg17=elem.pos(1)[1]*var_inter[0]; reg7=reg7/reg13; reg4=reg4/reg13; T reg18=reg6*elem.pos(1)[0];
    T reg19=reg6*elem.pos(0)[0]; reg8=reg8/reg13; reg12=reg12-reg11; T reg20=var_inter[1]*elem.pos(2)[1]; T reg21=reg8*reg7;
    T reg22=reg4*reg7; T reg23=elem.pos(2)[1]*var_inter[0]; T reg24=reg14+reg15; T reg25=elem.pos(2)[0]*var_inter[0]; T reg26=reg16+reg17;
    T reg27=var_inter[1]*elem.pos(2)[0]; reg18=reg18-reg19; reg27=reg18+reg27; reg18=reg9*elem.pos(3)[1]; T reg28=var_inter[1]*elem.pos(3)[0];
    T reg29=reg8*reg21; reg20=reg12+reg20; reg12=reg4*reg22; reg23=reg23-reg26; T reg30=var_inter[1]*elem.pos(3)[1];
    T reg31=(*f.m).alpha*reg8; reg13=reg5/reg13; reg5=(*f.m).alpha*reg4; T reg32=reg9*elem.pos(3)[0]; reg25=reg25-reg24;
    reg20=reg20-reg30; reg32=reg25+reg32; reg18=reg23+reg18; reg27=reg27-reg28; reg12=reg29-reg12;
    reg5=reg31+reg5; reg13=(*f.m).alpha*reg13; reg23=reg2*reg0; reg25=reg20*reg32; reg29=reg27*reg18;
    reg31=reg3*reg0; reg13=reg5+reg13; reg21=reg21/reg12; reg12=reg22/reg12; reg5=pow(reg3,2);
    reg25=reg29-reg25; reg22=reg2*reg23; reg29=pow(reg2,2); T reg33=reg3*reg31; reg12=reg12*reg13;
    reg13=reg21*reg13; reg21=1-(*f.m).resolution; reg12=reg13-reg12; reg27=reg27/reg25; reg20=reg20/reg25;
    reg32=reg32/reg25; reg18=reg18/reg25; reg5=reg29-reg5; reg33=reg22-reg33; reg13=reg9*reg20;
    reg22=reg9*reg27; reg29=var_inter[0]*reg20; T reg34=reg6*reg18; T reg35=var_inter[1]*reg32; reg31=reg31/reg33;
    reg5=reg5/reg33; T reg36=var_inter[1]*reg18; reg33=reg23/reg33; reg23=(*f.m).resolution*(*f.m).alpha; reg12=reg21*reg12;
    T reg37=reg22+reg35; T reg38=reg6*reg32; reg12=reg23+reg12; reg23=(*f.m).resolution*reg33; reg8=reg8*reg21;
    reg7=reg21*reg7; reg21=reg4*reg21; reg4=(*f.m).resolution*reg31; T reg39=(*f.m).resolution*reg5; T reg40=reg34+reg29;
    T reg41=reg13+reg36; T reg42=var_inter[0]*reg27; T reg43=reg38-reg22; T reg44=reg36-reg29; T reg45=reg13-reg34;
    T reg46=reg42-reg35; T reg47=0.5*reg40; T reg48=reg38+reg42; reg12=(*f.m).deltaT*reg12; T reg49=0.5*reg41;
    reg23=reg8+reg23; reg8=0.5*reg37; reg21=reg4+reg21; reg7=reg39+reg7; reg4=reg23*reg12;
    reg39=reg47*reg7; T reg50=reg21*reg12; T reg51=reg49*reg7; T reg52=reg8*reg7; T reg53=0.5*reg48;
    T reg54=0.5*reg44; T reg55=0.5*reg46; T reg56=0.5*reg45; T reg57=0.5*reg43; T reg58=reg56*reg7;
    T reg59=2*reg51; T reg60=var_inter[1]*reg9; T reg61=reg37*reg21; T reg62=reg57*reg7; reg39=2*reg39;
    T reg63=reg53*reg7; T reg64=reg48*reg21; T reg65=reg6*var_inter[0]; T reg66=reg4+reg50; T reg67=reg55*reg7;
    reg52=2*reg52; T reg68=reg54*reg7; T reg69=reg37*reg23; T reg70=reg41*reg23; T reg71=reg8*reg59;
    reg58=2*reg58; T reg72=reg43*reg21; T reg73=reg41*reg61; T reg74=reg53*reg39; T reg75=reg40*reg64;
    T reg76=reg44*reg21; T reg77=reg48*reg66; T reg78=reg40*reg70; T reg79=2*reg63; T reg80=reg43*reg23;
    T reg81=reg53*reg52; T reg82=reg41*reg66; T reg83=reg65*(*f.m).f_vol[1]; T reg84=reg41*reg21; T reg85=reg48*reg23;
    T reg86=reg40*reg21; T reg87=reg46*reg23; T reg88=reg40*reg23; T reg89=reg6*reg9; T reg90=reg44*reg23;
    T reg91=reg45*reg23; T reg92=reg69*reg48; reg62=2*reg62; T reg93=reg46*reg21; T reg94=reg60*(*f.m).f_vol[0];
    T reg95=var_inter[1]*var_inter[0]; reg68=2*reg68; T reg96=reg47*reg59; reg67=2*reg67; T reg97=reg56*reg67;
    T reg98=reg69*reg43; T reg99=reg57*reg68; T reg100=reg53*reg68; T reg101=reg56*reg59; T reg102=reg43*reg85;
    T reg103=reg84*reg43; T reg104=reg40*reg88; T reg105=reg79*reg53; T reg106=reg43*reg76; T reg107=reg45*reg93;
    T reg108=reg45*reg70; reg74=reg75+reg74; T reg109=reg47*reg39; T reg110=reg48*reg85; T reg111=reg40*reg90;
    reg73=reg71+reg73; T reg112=reg57*reg52; T reg113=reg57*reg79; T reg114=reg45*reg88; T reg115=reg57*reg59;
    T reg116=reg40*reg66; T reg117=reg45*reg61; T reg118=reg57*reg39; T reg119=reg77-reg83; T reg120=reg44*reg66;
    T reg121=reg56*reg58; T reg122=reg46*reg66; T reg123=reg43*reg80; T reg124=reg82-reg94; T reg125=reg45*reg64;
    T reg126=reg56*reg79; T reg127=reg37*reg66; T reg128=reg43*reg86; T reg129=reg56*reg68; T reg130=reg57*reg67;
    T reg131=reg87*reg43; T reg132=reg56*reg39; T reg133=reg56*reg52; T reg134=reg45*reg90; T reg135=reg44*reg90;
    T reg136=reg55*reg67; T reg137=reg47*reg67; T reg138=reg95*(*f.m).f_vol[0]; T reg139=reg48*reg76; T reg140=reg40*reg61;
    T reg141=reg47*reg68; T reg142=reg87*reg48; T reg143=reg65*(*f.m).f_vol[0]; T reg144=reg53*reg59; T reg145=reg47*reg52;
    T reg146=reg84*reg48; T reg147=reg89*(*f.m).f_vol[1]; T reg148=reg69*reg37; T reg149=reg45*reg72; T reg150=reg57*reg58;
    T reg151=reg49*reg59; T reg152=reg45*reg91; T reg153=reg57*reg62; T reg154=reg60*(*f.m).f_vol[1]; T reg155=reg95*(*f.m).f_vol[1];
    reg92=reg96+reg92; T reg156=reg54*reg68; T reg157=reg87*reg46; T reg158=reg45*reg66; T reg159=reg54*reg52;
    T reg160=reg84*reg46; T reg161=reg43*reg66; T reg162=reg54*reg59; reg69=reg69*reg46; T reg163=reg89*(*f.m).f_vol[0];
    T reg164=reg8*reg52; T reg165=reg55*reg52; T reg166=reg44*reg70; T reg167=reg44*reg93; T reg168=reg40*reg93;
    T reg169=reg41*reg70; reg61=reg44*reg61; T reg170=reg55*reg68; reg81=reg78+reg81; T reg171=reg55*reg59;
    T reg172=reg53*reg67; reg114=reg114-reg113; reg152=reg153+reg152; reg61=reg61-reg171; reg148=reg148+reg151;
    reg145=reg145+reg146; reg153=reg74*reg25; reg100=reg168-reg100; reg118=reg118-reg125; reg168=reg92*reg25;
    reg109=reg109+reg110; reg157=reg156+reg157; reg106=reg97+reg106; reg97=reg163+reg158; reg134=reg130+reg134;
    reg112=reg112-reg108; reg159=reg159-reg160; reg130=reg73*reg25; reg164=reg164+reg169; reg107=reg99+reg107;
    reg172=reg111-reg172; reg69=reg69-reg162; reg99=reg147+reg161; reg135=reg136+reg135; reg111=reg25*reg81;
    reg131=reg129+reg131; reg133=reg133-reg103; reg128=reg128-reg126; reg139=reg137-reg139; reg129=reg127+reg154;
    reg132=reg132-reg102; reg124=reg124*reg25; reg167=reg170+reg167; reg123=reg121+reg123; reg142=reg141-reg142;
    reg121=reg122+reg155; reg98=reg98-reg101; reg104=reg104+reg105; reg136=reg116+reg143; reg149=reg150+reg149;
    reg165=reg165-reg166; reg119=reg119*reg25; reg117=reg117-reg115; reg137=reg120+reg138; reg140=reg140+reg144;
    reg123=reg25*reg123; reg112=reg25*reg112; reg117=reg25*reg117; reg141=ponderation*reg153; reg148=reg25*reg148;
    reg150=ponderation*reg111; reg156=reg25*reg97; reg140=reg25*reg140; reg128=reg25*reg128; reg100=reg25*reg100;
    reg172=reg25*reg172; reg132=reg25*reg132; reg170=reg25*reg99; reg104=reg104*reg25; reg142=reg142*reg25;
    reg98=reg98*reg25; reg109=reg109*reg25; reg133=reg133*reg25; reg139=reg139*reg25; reg131=reg131*reg25;
    reg135=reg135*reg25; reg149=reg149*reg25; reg106=reg25*reg106; T reg173=reg129*reg25; reg152=reg152*reg25;
    reg145=reg145*reg25; reg124=ponderation*reg124; T reg174=ponderation*reg130; T reg175=reg121*reg25; reg165=reg165*reg25;
    reg107=reg25*reg107; reg69=reg69*reg25; reg134=reg25*reg134; reg61=reg61*reg25; reg118=reg25*reg118;
    reg159=reg159*reg25; reg114=reg25*reg114; T reg176=reg136*reg25; reg157=reg157*reg25; reg119=ponderation*reg119;
    reg164=reg164*reg25; T reg177=ponderation*reg168; T reg178=reg137*reg25; reg167=reg167*reg25; T tmp_4_5=ponderation*reg167;
    T tmp_4_4=ponderation*reg135; T tmp_2_5=ponderation*reg100; T tmp_2_6=-reg150; T tmp_3_4=ponderation*reg139; T tmp_3_5=ponderation*reg142;
    T tmp_0_1=ponderation*reg149; T tmp_2_7=ponderation*reg140; T tmp_0_0=ponderation*reg152; T tmp_7_7=ponderation*reg148; T tmp_3_6=ponderation*reg145;
    T tmp_3_7=-reg177; T tmp_5_5=ponderation*reg157; reg100=ponderation*reg156; T vec_0=reg100; T tmp_5_6=ponderation*reg159;
    T tmp_5_7=ponderation*reg69; T tmp_0_5=ponderation*reg107; T tmp_0_4=ponderation*reg134; T tmp_0_6=ponderation*reg112; T tmp_0_3=ponderation*reg118;
    T tmp_0_2=ponderation*reg114; reg69=ponderation*reg170; T vec_1=reg69; reg107=ponderation*reg176; T vec_2=reg107;
    T vec_3=-reg119; T tmp_0_7=ponderation*reg117; reg112=ponderation*reg178; T vec_4=reg112; reg114=ponderation*reg175;
    T vec_5=reg114; T tmp_1_1=ponderation*reg123; T vec_6=-reg124; reg117=ponderation*reg173; T vec_7=reg117;
    T tmp_1_4=ponderation*reg106; T tmp_1_2=ponderation*reg128; T tmp_1_5=ponderation*reg131; T tmp_1_6=ponderation*reg133; T tmp_1_3=ponderation*reg132;
    T tmp_1_7=ponderation*reg98; T tmp_2_2=ponderation*reg104; T tmp_6_7=-reg174; T tmp_3_3=ponderation*reg109; T tmp_2_3=-reg141;
    T tmp_6_6=ponderation*reg164; T tmp_4_7=ponderation*reg61; T tmp_2_4=ponderation*reg172; T tmp_4_6=ponderation*reg165;
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
    T reg4=elem.pos(1)[0]*var_inter[0]; T reg5=reg2*elem.pos(0)[0]; T reg6=reg3*elem.pos(1)[1]; T reg7=elem.pos(0)[1]*reg2; T reg8=elem.pos(1)[1]*var_inter[0];
    T reg9=reg0*reg1; T reg10=reg3*elem.pos(0)[1]; T reg11=reg3*elem.pos(1)[0]; T reg12=reg3*elem.pos(0)[0]; T reg13=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg14=1.0/(*f.m).elastic_modulus; T reg15=reg14*reg9; reg9=reg13*reg9; T reg16=elem.pos(2)[1]*var_inter[0]; T reg17=reg7+reg8;
    T reg18=elem.pos(2)[0]*var_inter[0]; T reg19=reg5+reg4; T reg20=reg14*reg1; T reg21=var_inter[1]*elem.pos(2)[1]; reg6=reg6-reg10;
    reg1=reg13*reg1; reg11=reg11-reg12; T reg22=var_inter[1]*elem.pos(2)[0]; reg18=reg18-reg19; T reg23=reg2*elem.pos(3)[0];
    reg16=reg16-reg17; T reg24=reg2*elem.pos(3)[1]; T reg25=var_inter[1]*elem.pos(3)[1]; reg22=reg11+reg22; reg11=reg13*reg20;
    reg21=reg6+reg21; reg6=reg13*reg1; T reg26=reg13*reg15; T reg27=reg13*reg9; reg15=reg14*reg15;
    reg20=reg14*reg20; T reg28=var_inter[1]*elem.pos(3)[0]; reg22=reg22-reg28; reg21=reg21-reg25; reg23=reg18+reg23;
    reg24=reg16+reg24; reg26=reg27+reg26; reg15=reg15-reg27; reg9=reg14*reg9; reg1=reg14*reg1;
    reg11=reg6+reg11; reg20=reg20-reg6; reg16=reg13*reg26; reg9=reg27+reg9; reg18=reg21*reg23;
    reg27=reg14*reg15; T reg29=reg22*reg24; reg20=reg14*reg20; reg11=reg13*reg11; T reg30=reg6+reg1;
    T reg31=reg13*reg0; T reg32=reg14*reg0; T reg33=pow(reg13,2); reg11=reg20-reg11; reg20=reg13*reg31;
    T reg34=pow(reg14,2); T reg35=reg13*reg9; reg18=reg29-reg18; reg29=reg14*reg32; reg30=reg13*reg30;
    reg16=reg27-reg16; reg33=reg34-reg33; reg21=reg21/reg18; reg23=reg23/reg18; reg30=reg11-reg30;
    reg35=reg16-reg35; reg22=reg22/reg18; reg20=reg29-reg20; reg24=reg24/reg18; reg11=1-(*f.m).resolution;
    reg30=reg30/reg35; reg16=var_inter[1]*reg24; reg27=var_inter[0]*reg22; reg29=var_inter[1]*reg23; reg34=reg3*reg23;
    T reg36=var_inter[0]*reg21; reg33=reg33/reg20; T reg37=reg2*reg22; T reg38=reg3*reg24; T reg39=reg2*reg21;
    reg32=reg32/reg20; reg20=reg31/reg20; reg31=(*f.m).resolution*reg33; T reg40=reg39+reg16; T reg41=reg37+reg29;
    T reg42=reg38+reg36; T reg43=reg34+reg27; reg26=reg26/reg35; reg15=reg15/reg35; T reg44=reg11*reg30;
    T reg45=reg16-reg36; T reg46=(*f.m).resolution*reg20; T reg47=reg15*reg11; T reg48=reg26*reg11; T reg49=reg27-reg29;
    T reg50=0.5*reg42; T reg51=0.5*reg43; T reg52=0.5*reg41; T reg53=0.5*reg40; T reg54=(*f.m).resolution*reg32;
    T reg55=reg34-reg37; reg44=reg31+reg44; reg31=reg39-reg38; T reg56=reg53*reg44; reg54=reg47+reg54;
    reg48=reg46+reg48; reg46=reg52*reg44; reg47=0.5*reg31; T reg57=0.5*reg49; T reg58=reg51*reg44;
    T reg59=reg50*reg44; T reg60=0.5*reg55; T reg61=0.5*reg45; T reg62=reg60*reg44; T reg63=reg41*reg54;
    T reg64=2*reg58; T reg65=reg47*reg44; T reg66=reg40*reg54; T reg67=2*reg56; T reg68=reg41*reg48;
    reg46=2*reg46; T reg69=reg42*reg48; T reg70=reg61*reg44; T reg71=reg43*reg54; T reg72=reg57*reg44;
    T reg73=reg43*reg48; reg59=2*reg59; T reg74=reg42*reg54; T reg75=reg40*reg48; reg65=2*reg65;
    T reg76=reg55*reg48; T reg77=reg52*reg67; T reg78=reg40*reg68; T reg79=reg40*reg74; T reg80=reg31*reg54;
    reg62=2*reg62; T reg81=reg42*reg66; T reg82=reg51*reg46; T reg83=reg53*reg59; T reg84=reg41*reg71;
    T reg85=reg43*reg69; T reg86=reg42*reg73; T reg87=reg51*reg59; T reg88=reg45*reg48; T reg89=reg75*reg41;
    T reg90=reg49*reg54; T reg91=reg53*reg46; T reg92=reg64*reg50; T reg93=reg55*reg54; T reg94=reg31*reg48;
    T reg95=reg50*reg67; T reg96=reg52*reg64; reg72=2*reg72; T reg97=reg45*reg54; reg70=2*reg70;
    T reg98=reg49*reg48; T reg99=reg63*reg43; T reg100=reg57*reg67; reg82=reg81+reg82; T reg101=reg42*reg68;
    T reg102=reg51*reg67; T reg103=reg62*reg50; T reg104=reg43*reg94; T reg105=reg65*reg50; T reg106=reg45*reg68;
    T reg107=reg62*reg61; T reg108=reg43*reg93; T reg109=reg41*reg69; T reg110=reg49*reg94; reg85=reg92+reg85;
    T reg111=reg90*reg41; T reg112=reg40*reg73; T reg113=reg65*reg61; T reg114=reg53*reg70; T reg115=reg52*reg72;
    T reg116=reg49*reg93; T reg117=reg64*reg61; T reg118=reg40*reg97; reg91=reg89+reg91; reg78=reg77+reg78;
    T reg119=reg49*reg69; T reg120=reg63*reg41; T reg121=reg52*reg70; T reg122=reg53*reg67; T reg123=reg40*reg66;
    T reg124=reg52*reg59; T reg125=reg52*reg46; T reg126=reg55*reg93; T reg127=reg47*reg65; T reg128=reg55*reg94;
    T reg129=reg47*reg62; reg68=reg31*reg68; T reg130=reg60*reg67; T reg131=reg31*reg66; T reg132=reg60*reg46;
    T reg133=reg31*reg98; T reg134=reg60*reg70; T reg135=reg51*reg70; T reg136=reg31*reg97; T reg137=reg60*reg72;
    T reg138=reg31*reg73; T reg139=reg60*reg59; T reg140=reg31*reg74; T reg141=reg60*reg64; T reg142=reg40*reg98;
    T reg143=reg42*reg98; T reg144=reg51*reg72; T reg145=reg42*reg97; T reg146=reg53*reg72; T reg147=reg41*reg88;
    reg83=reg84+reg83; T reg148=reg53*reg64; T reg149=reg53*reg65; reg93=reg41*reg93; T reg150=reg53*reg62;
    reg94=reg41*reg94; T reg151=reg55*reg88; T reg152=reg47*reg72; T reg153=reg55*reg71; T reg154=reg47*reg59;
    reg69=reg55*reg69; T reg155=reg47*reg64; T reg156=reg50*reg72; reg97=reg45*reg97; T reg157=reg57*reg72;
    T reg158=reg43*reg88; T reg159=reg47*reg70; T reg160=reg90*reg55; T reg161=reg76*reg45; T reg162=reg47*reg46;
    T reg163=reg52*reg62; T reg164=reg50*reg70; T reg165=reg90*reg43; T reg166=reg63*reg55; T reg167=reg45*reg73;
    T reg168=reg57*reg59; T reg169=reg40*reg80; T reg170=reg65*reg57; T reg171=reg50*reg46; T reg172=reg75*reg43;
    T reg173=reg52*reg65; reg99=reg95+reg99; T reg174=reg31*reg76; T reg175=reg60*reg65; T reg176=reg40*reg76;
    T reg177=reg47*reg67; T reg178=reg31*reg80; T reg179=reg60*reg62; reg79=reg96+reg79; T reg180=reg62*reg57;
    T reg181=reg80*reg45; reg72=reg61*reg72; reg88=reg49*reg88; reg76=reg76*reg42; T reg182=reg45*reg66;
    T reg183=reg57*reg46; reg65=reg65*reg51; reg62=reg62*reg51; T reg184=reg61*reg70; reg90=reg90*reg49;
    T reg185=reg42*reg74; T reg186=reg49*reg71; reg98=reg45*reg98; T reg187=reg64*reg51; T reg188=reg61*reg59;
    reg63=reg63*reg49; T reg189=reg61*reg67; T reg190=reg75*reg55; T reg191=reg43*reg71; reg59=reg50*reg59;
    reg80=reg80*reg42; T reg192=reg64*reg57; reg87=reg86+reg87; T reg193=reg75*reg49; reg46=reg61*reg46;
    reg70=reg57*reg70; reg74=reg45*reg74; reg104=reg103-reg104; reg139=reg139-reg138; reg188=reg188-reg186;
    reg162=reg162-reg190; reg140=reg140-reg141; reg136=reg137+reg136; reg181=reg180+reg181; reg161=reg170+reg161;
    reg144=reg145-reg144; reg146=reg147-reg146; reg62=reg80-reg62; reg80=reg18*reg83; reg65=reg76-reg65;
    reg142=reg121-reg142; reg109=reg148+reg109; reg149=reg93-reg149; reg185=reg185+reg187; reg135=reg143-reg135;
    reg150=reg94-reg150; reg151=reg152+reg151; reg76=reg87*reg18; reg154=reg154-reg153; reg69=reg69-reg155;
    reg59=reg59+reg191; reg93=reg18*reg82; reg126=reg127+reg126; reg158=reg156-reg158; reg128=reg129+reg128;
    reg166=reg166-reg177; reg68=reg68-reg130; reg165=reg164-reg165; reg101=reg101+reg102; reg171=reg171+reg172;
    reg132=reg132-reg131; reg94=reg99*reg18; reg133=reg134+reg133; reg176=reg173-reg176; reg103=reg18*reg85;
    reg178=reg179+reg178; reg98=reg70+reg98; reg116=reg113+reg116; reg174=reg175+reg174; reg119=reg119-reg117;
    reg160=reg159+reg160; reg74=reg74-reg192; reg70=reg79*reg18; reg113=reg78*reg18; reg114=reg111-reg114;
    reg111=reg18*reg91; reg125=reg125+reg123; reg168=reg168-reg167; reg97=reg157+reg97; reg124=reg124+reg112;
    reg118=reg115-reg118; reg108=reg105-reg108; reg88=reg72+reg88; reg183=reg183-reg182; reg90=reg184+reg90;
    reg46=reg46-reg193; reg106=reg106-reg100; reg120=reg120+reg122; reg63=reg63-reg189; reg110=reg107+reg110;
    reg169=reg163-reg169; reg124=reg124*reg18; reg59=reg59*reg18; reg62=reg62*reg18; reg109=reg18*reg109;
    reg69=reg18*reg69; reg166=reg166*reg18; reg97=reg97*reg18; reg110=reg110*reg18; reg144=reg18*reg144;
    reg72=ponderation*reg76; reg105=ponderation*reg70; reg154=reg18*reg154; reg98=reg98*reg18; reg146=reg18*reg146;
    reg183=reg183*reg18; reg107=ponderation*reg103; reg115=ponderation*reg111; reg142=reg18*reg142; reg119=reg119*reg18;
    reg106=reg106*reg18; reg135=reg18*reg135; reg185=reg185*reg18; reg150=reg18*reg150; reg121=ponderation*reg80;
    reg151=reg18*reg151; reg116=reg116*reg18; reg120=reg18*reg120; reg65=reg65*reg18; reg149=reg18*reg149;
    reg171=reg171*reg18; reg160=reg160*reg18; reg127=ponderation*reg94; reg178=reg178*reg18; reg133=reg18*reg133;
    reg176=reg176*reg18; reg181=reg181*reg18; reg104=reg18*reg104; reg136=reg18*reg136; reg169=reg169*reg18;
    reg63=reg63*reg18; reg161=reg161*reg18; reg162=reg162*reg18; reg139=reg18*reg139; reg46=reg46*reg18;
    reg90=reg90*reg18; reg108=reg18*reg108; reg88=reg88*reg18; reg140=reg18*reg140; reg188=reg188*reg18;
    reg129=ponderation*reg93; reg126=reg18*reg126; reg118=reg118*reg18; reg158=reg158*reg18; reg128=reg18*reg128;
    reg168=reg168*reg18; reg125=reg125*reg18; reg114=reg18*reg114; reg68=reg18*reg68; reg174=reg174*reg18;
    reg132=reg18*reg132; reg101=reg18*reg101; reg74=reg74*reg18; reg134=ponderation*reg113; reg165=reg165*reg18;
    T tmp_1_7=ponderation*reg166; T tmp_1_6=ponderation*reg162; T tmp_3_0=ponderation*reg104; T tmp_1_5=ponderation*reg160; T tmp_2_6=-reg129;
    T tmp_2_5=ponderation*reg135; T tmp_7_7=ponderation*reg120; T tmp_3_1=ponderation*reg108; T tmp_7_5=ponderation*reg114; T tmp_2_4=ponderation*reg144;
    T tmp_2_7=ponderation*reg101; T tmp_7_6=-reg115; T tmp_1_4=ponderation*reg151; T tmp_0_3=ponderation*reg139; T tmp_0_2=ponderation*reg140;
    T tmp_5_3=ponderation*reg188; T tmp_5_4=ponderation*reg88; T tmp_5_5=ponderation*reg90; T tmp_5_6=ponderation*reg46; T tmp_5_7=ponderation*reg63;
    T tmp_6_0=ponderation*reg169; T tmp_6_1=ponderation*reg176; T tmp_0_0=ponderation*reg178; T tmp_0_1=ponderation*reg174; T tmp_4_2=ponderation*reg74;
    T tmp_6_6=ponderation*reg125; T tmp_6_5=ponderation*reg142; T tmp_4_3=ponderation*reg168; T tmp_6_4=ponderation*reg118; T tmp_6_3=ponderation*reg124;
    T tmp_6_2=-reg105; T tmp_5_2=ponderation*reg119; T tmp_4_4=ponderation*reg97; T tmp_5_1=ponderation*reg116; T tmp_5_0=ponderation*reg110;
    T tmp_4_5=ponderation*reg98; T tmp_4_7=ponderation*reg106; T tmp_4_6=ponderation*reg183; T tmp_2_0=ponderation*reg62; T tmp_2_3=-reg72;
    T tmp_7_4=ponderation*reg146; T tmp_7_3=-reg121; T tmp_2_1=ponderation*reg65; T tmp_7_2=ponderation*reg109; T tmp_7_1=ponderation*reg149;
    T tmp_7_0=ponderation*reg150; T tmp_2_2=ponderation*reg185; T tmp_6_7=-reg134; T tmp_1_3=ponderation*reg154; T tmp_3_2=-reg107;
    T tmp_1_2=ponderation*reg69; T tmp_3_3=ponderation*reg59; T tmp_1_1=ponderation*reg126; T tmp_1_0=ponderation*reg128; T tmp_3_4=ponderation*reg158;
    T tmp_0_7=ponderation*reg68; T tmp_3_5=ponderation*reg165; T tmp_0_6=ponderation*reg132; T tmp_3_6=ponderation*reg171; T tmp_0_5=ponderation*reg133;
    T tmp_3_7=-reg127; T tmp_4_0=ponderation*reg181; T tmp_0_4=ponderation*reg136; T tmp_4_1=ponderation*reg161;
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
    T reg4=elem.pos(1)[0]*var_inter[0]; T reg5=reg3*elem.pos(0)[0]; T reg6=elem.pos(0)[1]*reg3; T reg7=elem.pos(1)[1]*var_inter[0]; T reg8=reg1*elem.pos(1)[1];
    T reg9=reg1*elem.pos(0)[1]; T reg10=reg0*reg2; T reg11=reg1*elem.pos(1)[0]; T reg12=reg1*elem.pos(0)[0]; T reg13=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg14=1.0/(*f.m).elastic_modulus; T reg15=reg13*reg2; reg8=reg8-reg9; T reg16=var_inter[1]*elem.pos(2)[1]; T reg17=reg5+reg4;
    T reg18=elem.pos(2)[0]*var_inter[0]; T reg19=reg6+reg7; T reg20=elem.pos(2)[1]*var_inter[0]; T reg21=reg13*reg10; reg10=reg14*reg10;
    reg11=reg11-reg12; T reg22=var_inter[1]*elem.pos(2)[0]; reg2=reg14*reg2; T reg23=reg14*reg10; T reg24=var_inter[1]*elem.pos(3)[1];
    reg16=reg8+reg16; reg22=reg11+reg22; reg8=var_inter[1]*elem.pos(3)[0]; reg18=reg18-reg17; reg11=reg3*elem.pos(3)[0];
    T reg25=reg14*reg2; T reg26=reg13*reg15; T reg27=reg13*reg21; reg10=reg13*reg10; reg20=reg20-reg19;
    T reg28=reg3*elem.pos(3)[1]; reg2=reg13*reg2; reg11=reg18+reg11; reg2=reg26+reg2; reg28=reg20+reg28;
    reg23=reg23-reg27; reg10=reg27+reg10; reg16=reg16-reg24; reg25=reg25-reg26; reg15=reg14*reg15;
    reg22=reg22-reg8; reg21=reg14*reg21; reg21=reg27+reg21; reg18=reg22*reg28; reg20=reg14*reg0;
    reg27=reg13*reg0; T reg29=reg13*reg10; T reg30=reg26+reg15; T reg31=reg16*reg11; reg25=reg14*reg25;
    T reg32=reg14*reg23; reg2=reg13*reg2; T reg33=reg13*reg27; T reg34=pow(reg13,2); T reg35=reg13*reg21;
    reg29=reg32-reg29; reg32=pow(reg14,2); reg31=reg18-reg31; reg30=reg13*reg30; reg18=reg14*reg20;
    reg2=reg25-reg2; reg34=reg32-reg34; reg16=reg16/reg31; reg11=reg11/reg31; reg28=reg28/reg31;
    reg30=reg2-reg30; reg33=reg18-reg33; reg22=reg22/reg31; reg35=reg29-reg35; reg30=reg30/reg35;
    reg2=var_inter[1]*reg28; reg18=var_inter[1]*reg11; reg34=reg34/reg33; reg25=1-(*f.m).resolution; reg29=reg3*reg16;
    reg32=reg1*reg28; T reg36=var_inter[0]*reg16; T reg37=reg3*reg22; T reg38=reg37+reg18; T reg39=reg29+reg2;
    T reg40=reg1*reg11; reg23=reg23/reg35; reg27=reg27/reg33; T reg41=reg32+reg36; reg33=reg20/reg33;
    reg20=var_inter[0]*reg22; T reg42=(*f.m).resolution*reg34; reg10=reg10/reg35; T reg43=reg25*reg30; T reg44=(*f.m).resolution*reg33;
    T reg45=reg10*reg25; T reg46=reg29-reg32; T reg47=reg23*reg25; T reg48=0.5*reg39; T reg49=0.5*reg38;
    T reg50=reg40-reg37; T reg51=(*f.m).resolution*reg27; T reg52=reg40+reg20; T reg53=0.5*reg41; T reg54=reg20-reg18;
    T reg55=reg2-reg36; reg43=reg42+reg43; reg44=reg47+reg44; reg42=reg49*reg43; reg45=reg51+reg45;
    reg47=0.5*reg46; reg51=0.5*reg50; T reg56=reg48*reg43; T reg57=reg53*reg43; T reg58=0.5*reg52;
    T reg59=0.5*reg54; T reg60=0.5*reg55; T reg61=reg47*reg43; T reg62=reg59*reg43; T reg63=reg51*reg43;
    T reg64=reg52*reg45; reg57=2*reg57; T reg65=reg38*reg44; T reg66=reg60*reg43; reg42=2*reg42;
    T reg67=reg39*reg44; T reg68=reg58*reg43; T reg69=2*reg56; T reg70=reg38*reg45; reg61=2*reg61;
    T reg71=reg41*reg44; T reg72=reg50*reg45; T reg73=2*reg68; T reg74=reg39*reg45; T reg75=reg54*reg44;
    T reg76=reg39*reg70; T reg77=reg55*reg45; T reg78=reg49*reg69; T reg79=reg58*reg42; T reg80=reg52*reg44;
    T reg81=reg41*reg45; T reg82=reg58*reg57; T reg83=reg41*reg64; T reg84=reg50*reg44; T reg85=reg41*reg67;
    T reg86=reg65*reg52; T reg87=reg53*reg69; T reg88=reg54*reg45; reg66=2*reg66; T reg89=reg55*reg44;
    reg62=2*reg62; reg63=2*reg63; T reg90=reg46*reg44; T reg91=reg73*reg58; T reg92=reg74*reg50;
    T reg93=reg41*reg71; reg82=reg83+reg82; reg76=reg78+reg76; T reg94=reg65*reg50; T reg95=reg47*reg42;
    T reg96=reg47*reg69; T reg97=reg50*reg77; T reg98=reg47*reg62; T reg99=reg50*reg80; T reg100=reg47*reg57;
    T reg101=reg50*reg81; T reg102=reg47*reg73; T reg103=reg50*reg84; T reg104=reg47*reg61; T reg105=reg46*reg70;
    T reg106=reg51*reg69; T reg107=reg41*reg89; T reg108=reg46*reg67; T reg109=reg51*reg42; T reg110=reg58*reg62;
    T reg111=reg46*reg88; T reg112=reg51*reg66; T reg113=reg58*reg66; T reg114=reg46*reg89; T reg115=reg51*reg62;
    T reg116=reg46*reg64; T reg117=reg51*reg57; T reg118=reg46*reg71; T reg119=reg51*reg73; T reg120=reg41*reg88;
    T reg121=reg47*reg66; T reg122=reg75*reg50; T reg123=reg59*reg42; T reg124=reg75*reg52; T reg125=reg53*reg66;
    T reg126=reg48*reg69; T reg127=reg41*reg70; T reg128=reg46*reg90; T reg129=reg51*reg63; T reg130=reg55*reg88;
    T reg131=reg59*reg66; T reg132=reg52*reg80; T reg133=reg60*reg42; T reg134=reg52*reg77; T reg135=reg74*reg54;
    T reg136=reg53*reg62; T reg137=reg53*reg57; T reg138=reg65*reg54; T reg139=reg59*reg62; T reg140=reg55*reg89;
    T reg141=reg60*reg66; T reg142=reg75*reg54; reg86=reg87+reg86; reg65=reg65*reg38; T reg143=reg39*reg67;
    T reg144=reg49*reg42; reg79=reg85+reg79; T reg145=reg58*reg69; T reg146=reg74*reg52; T reg147=reg53*reg42;
    T reg148=reg51*reg61; reg70=reg55*reg70; T reg149=reg59*reg69; T reg150=reg46*reg72; T reg151=reg60*reg69;
    T reg152=reg55*reg67; reg128=reg129+reg128; reg65=reg65+reg126; reg117=reg117-reg116; reg105=reg105-reg106;
    reg127=reg127+reg145; reg134=reg136-reg134; reg138=reg138-reg151; reg114=reg115+reg114; reg110=reg107-reg110;
    reg124=reg125-reg124; reg111=reg112+reg111; reg107=reg86*reg31; reg147=reg147+reg146; reg109=reg109-reg108;
    reg112=reg31*reg79; reg100=reg100-reg99; reg123=reg123-reg152; reg142=reg141+reg142; reg70=reg70-reg149;
    reg93=reg93+reg91; reg97=reg98+reg97; reg94=reg94-reg96; reg130=reg131+reg130; reg101=reg101-reg102;
    reg95=reg95-reg92; reg133=reg133-reg135; reg140=reg139+reg140; reg113=reg120-reg113; reg122=reg121+reg122;
    reg144=reg144+reg143; reg137=reg137+reg132; reg98=reg82*reg31; reg103=reg104+reg103; reg150=reg148+reg150;
    reg118=reg118-reg119; reg104=reg76*reg31; reg109=reg31*reg109; reg100=reg31*reg100; reg115=ponderation*reg107;
    reg142=reg142*reg31; reg65=reg31*reg65; reg101=reg31*reg101; reg138=reg138*reg31; reg120=ponderation*reg98;
    reg105=reg31*reg105; reg103=reg31*reg103; reg133=reg133*reg31; reg144=reg144*reg31; reg121=ponderation*reg104;
    reg93=reg93*reg31; reg70=reg70*reg31; reg123=reg123*reg31; reg125=ponderation*reg112; reg130=reg130*reg31;
    reg94=reg94*reg31; reg113=reg31*reg113; reg95=reg95*reg31; reg140=reg140*reg31; reg122=reg122*reg31;
    reg97=reg31*reg97; reg150=reg150*reg31; reg118=reg31*reg118; reg128=reg128*reg31; reg147=reg147*reg31;
    reg111=reg31*reg111; reg137=reg137*reg31; reg117=reg31*reg117; reg110=reg31*reg110; reg124=reg124*reg31;
    reg134=reg134*reg31; reg114=reg31*reg114; reg127=reg31*reg127; T tmp_7_7=ponderation*reg65; T tmp_2_7=ponderation*reg127;
    T tmp_2_5=ponderation*reg113; T tmp_6_7=-reg121; T tmp_2_6=-reg125; T tmp_2_3=-reg120; T tmp_2_4=ponderation*reg110;
    T tmp_3_7=-reg115; T tmp_5_5=ponderation*reg142; T tmp_5_6=ponderation*reg133; T tmp_5_7=ponderation*reg138; T tmp_3_6=ponderation*reg147;
    T tmp_3_5=ponderation*reg124; T tmp_3_4=ponderation*reg134; T tmp_3_3=ponderation*reg137; T tmp_0_0=ponderation*reg128; T tmp_0_1=ponderation*reg150;
    T tmp_4_4=ponderation*reg140; T tmp_4_5=ponderation*reg130; T tmp_4_6=ponderation*reg123; T tmp_4_7=ponderation*reg70; T tmp_6_6=ponderation*reg144;
    T tmp_2_2=ponderation*reg93; T tmp_1_7=ponderation*reg94; T tmp_1_6=ponderation*reg95; T tmp_1_5=ponderation*reg122; T tmp_1_4=ponderation*reg97;
    T tmp_0_2=ponderation*reg118; T tmp_0_3=ponderation*reg117; T tmp_0_4=ponderation*reg114; T tmp_0_5=ponderation*reg111; T tmp_0_6=ponderation*reg109;
    T tmp_0_7=ponderation*reg105; T tmp_1_1=ponderation*reg103; T tmp_1_2=ponderation*reg101; T tmp_1_3=ponderation*reg100;
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
    T reg7=reg3*reg2; T reg8=reg3*reg6; reg2=reg4*reg2; T reg9=reg4*reg5; T reg10=reg3*reg1;
    reg5=reg3*reg5; reg7=reg8+reg7; reg2=reg2-reg8; reg6=reg4*reg6; reg5=reg10+reg5;
    reg9=reg9-reg10; reg1=reg4*reg1; T reg11=reg4*reg2; T reg12=reg3*reg7; reg6=reg8+reg6;
    reg9=reg4*reg9; reg5=reg3*reg5; reg8=reg10+reg1; reg5=reg9-reg5; reg8=reg3*reg8;
    reg9=reg3*reg6; reg12=reg11-reg12; reg8=reg5-reg8; reg9=reg12-reg9; reg8=reg8/reg9;
    reg2=reg2/reg9; reg7=reg7/reg9; reg5=reg2*reg8; reg11=reg7*reg8; reg9=reg6/reg9;
    reg6=1-var_inter[1]; reg12=reg2*reg5; T reg13=reg7*reg11; T reg14=1-var_inter[0]; T reg15=(*f.m).alpha*reg2;
    T reg16=(*f.m).alpha*reg7; T reg17=elem.pos(1)[1]*var_inter[0]; T reg18=elem.pos(0)[1]*reg14; T reg19=elem.pos(1)[0]*var_inter[0]; T reg20=reg14*elem.pos(0)[0];
    reg9=(*f.m).alpha*reg9; reg16=reg15+reg16; reg15=reg6*elem.pos(1)[1]; T reg21=reg6*elem.pos(0)[1]; reg13=reg12-reg13;
    reg12=reg6*elem.pos(0)[0]; T reg22=reg6*elem.pos(1)[0]; reg11=reg11/reg13; reg9=reg16+reg9; reg16=reg3*reg0;
    reg13=reg5/reg13; reg5=reg4*reg0; reg15=reg15-reg21; T reg23=var_inter[1]*elem.pos(2)[1]; T reg24=var_inter[1]*elem.pos(2)[0];
    T reg25=reg20+reg19; T reg26=elem.pos(2)[0]*var_inter[0]; reg22=reg22-reg12; T reg27=reg18+reg17; T reg28=elem.pos(2)[1]*var_inter[0];
    reg28=reg28-reg27; T reg29=reg14*elem.pos(3)[1]; reg23=reg15+reg23; reg15=var_inter[1]*elem.pos(3)[1]; T reg30=var_inter[1]*elem.pos(3)[0];
    reg24=reg22+reg24; reg22=reg3*reg16; reg26=reg26-reg25; T reg31=reg14*elem.pos(3)[0]; reg13=reg13*reg9;
    reg9=reg11*reg9; reg11=reg4*reg5; reg22=reg11-reg22; reg11=1-(*f.m).resolution; reg9=reg13-reg9;
    reg29=reg28+reg29; reg31=reg26+reg31; reg24=reg24-reg30; reg23=reg23-reg15; reg5=reg5/reg22;
    reg13=reg24*reg29; reg16=reg16/reg22; reg9=reg11*reg9; reg26=(*f.m).resolution*(*f.m).alpha; reg28=reg23*reg31;
    reg9=reg26+reg9; reg26=(*f.m).resolution*reg5; reg2=reg2*reg11; reg7=reg7*reg11; T reg32=(*f.m).resolution*reg16;
    reg28=reg13-reg28; reg9=(*f.m).deltaT*reg9; reg23=reg23/reg28; reg31=reg31/reg28; reg7=reg32+reg7;
    reg24=reg24/reg28; reg29=reg29/reg28; reg26=reg2+reg26; reg2=reg26*reg9; reg13=var_inter[0]*reg24;
    reg32=reg6*reg31; T reg33=var_inter[1]*reg29; T reg34=reg14*reg23; T reg35=reg7*reg9; T reg36=reg6*reg29;
    T reg37=reg2+reg35; T reg38=var_inter[1]*reg14; T reg39=reg34+reg33; T reg40=reg6*var_inter[0]; T reg41=var_inter[1]*reg31;
    T reg42=reg14*reg24; T reg43=var_inter[0]*reg23; T reg44=reg32+reg13; T reg45=var_inter[1]*var_inter[0]; T reg46=reg44*reg37;
    T reg47=reg32-reg42; T reg48=reg42+reg41; T reg49=reg6*reg14; T reg50=reg34-reg36; T reg51=reg38*(*f.m).f_vol[0];
    T reg52=reg36+reg43; T reg53=reg33-reg43; T reg54=reg40*(*f.m).f_vol[1]; T reg55=reg13-reg41; T reg56=reg39*reg37;
    T reg57=reg52*reg37; T reg58=reg48*reg37; T reg59=reg47*reg37; T reg60=reg46-reg54; T reg61=reg49*(*f.m).f_vol[1];
    T reg62=reg45*(*f.m).f_vol[0]; T reg63=reg56-reg51; T reg64=reg53*reg37; T reg65=reg40*(*f.m).f_vol[0]; T reg66=reg45*(*f.m).f_vol[1];
    T reg67=reg38*(*f.m).f_vol[1]; T reg68=reg50*reg37; T reg69=reg55*reg37; T reg70=reg49*(*f.m).f_vol[0]; T reg71=reg57+reg65;
    reg60=reg60*reg28; reg63=reg63*reg28; T reg72=reg64+reg62; T reg73=reg70+reg68; T reg74=reg61+reg59;
    T reg75=reg58+reg67; T reg76=reg69+reg66; T reg77=reg28*reg74; T reg78=reg28*reg73; T reg79=reg71*reg28;
    reg60=ponderation*reg60; T reg80=reg72*reg28; T reg81=reg76*reg28; reg63=ponderation*reg63; T reg82=reg75*reg28;
    T reg83=ponderation*reg77; T vec_1=reg83; T reg84=ponderation*reg79; T vec_2=reg84; T vec_3=-reg60;
    reg60=ponderation*reg80; T vec_4=reg60; T reg85=ponderation*reg81; T vec_5=reg85; T reg86=ponderation*reg78;
    T vec_0=reg86; T vec_6=-reg63; reg63=ponderation*reg82; T vec_7=reg63;
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
    T reg4=reg0*reg1; T reg5=reg3*reg1; T reg6=reg2*reg4; reg4=reg3*reg4; reg1=reg2*reg1;
    T reg7=reg3*reg1; T reg8=reg3*reg5; T reg9=reg3*reg6; T reg10=reg3*reg4; reg1=reg2*reg1;
    reg6=reg2*reg6; reg5=reg2*reg5; reg1=reg1-reg8; reg7=reg8+reg7; reg4=reg2*reg4;
    reg6=reg6-reg10; reg9=reg10+reg9; reg7=reg3*reg7; reg1=reg2*reg1; T reg11=reg3*reg9;
    T reg12=reg2*reg6; T reg13=reg8+reg5; reg4=reg10+reg4; reg13=reg3*reg13; reg10=reg3*reg4;
    reg11=reg12-reg11; reg7=reg1-reg7; reg1=1-var_inter[1]; reg13=reg7-reg13; reg7=1-var_inter[0];
    reg10=reg11-reg10; reg9=reg9/reg10; reg6=reg6/reg10; reg13=reg13/reg10; reg11=elem.pos(1)[1]*var_inter[0];
    reg12=elem.pos(0)[1]*reg7; T reg14=elem.pos(1)[0]*var_inter[0]; T reg15=reg7*elem.pos(0)[0]; T reg16=reg1*elem.pos(0)[0]; T reg17=reg1*elem.pos(1)[0];
    T reg18=reg1*elem.pos(1)[1]; T reg19=reg1*elem.pos(0)[1]; reg17=reg17-reg16; T reg20=var_inter[1]*elem.pos(2)[0]; reg18=reg18-reg19;
    T reg21=var_inter[1]*elem.pos(2)[1]; T reg22=reg15+reg14; T reg23=elem.pos(2)[0]*var_inter[0]; T reg24=reg12+reg11; T reg25=elem.pos(2)[1]*var_inter[0];
    T reg26=reg9*reg13; T reg27=reg6*reg13; reg20=reg17+reg20; reg17=var_inter[1]*elem.pos(3)[0]; T reg28=reg6*reg27;
    reg21=reg18+reg21; reg18=var_inter[1]*elem.pos(3)[1]; T reg29=(*f.m).alpha*reg9; reg23=reg23-reg22; T reg30=(*f.m).alpha*reg6;
    T reg31=reg7*elem.pos(3)[0]; T reg32=reg9*reg26; reg10=reg4/reg10; reg4=reg7*elem.pos(3)[1]; reg25=reg25-reg24;
    T reg33=var_inter[0]*vectors[0][indices[1]+0]; T reg34=reg7*vectors[0][indices[0]+0]; T reg35=reg1*vectors[0][indices[0]+1]; reg32=reg28-reg32; reg28=reg1*vectors[0][indices[1]+1];
    T reg36=var_inter[0]*vectors[0][indices[1]+1]; T reg37=reg7*vectors[0][indices[0]+1]; T reg38=reg1*vectors[0][indices[0]+0]; reg29=reg30+reg29; reg10=(*f.m).alpha*reg10;
    reg4=reg25+reg4; reg31=reg23+reg31; reg21=reg21-reg18; reg20=reg20-reg17; reg23=reg1*vectors[0][indices[1]+0];
    reg25=reg20*reg4; reg33=reg34+reg33; reg30=var_inter[1]*vectors[0][indices[2]+1]; reg34=var_inter[0]*vectors[0][indices[2]+0]; reg27=reg27/reg32;
    reg32=reg26/reg32; reg36=reg37+reg36; reg26=var_inter[0]*vectors[0][indices[2]+1]; reg38=reg23-reg38; reg23=reg21*reg31;
    reg35=reg28-reg35; reg10=reg29+reg10; reg28=var_inter[1]*vectors[0][indices[2]+0]; reg23=reg25-reg23; reg25=var_inter[1]*vectors[0][indices[3]+0];
    reg38=reg28+reg38; reg32=reg32*reg10; reg10=reg27*reg10; reg27=var_inter[1]*vectors[0][indices[3]+1]; reg35=reg30+reg35;
    reg28=reg3*reg0; reg33=reg34-reg33; reg29=reg2*reg0; reg30=reg7*vectors[0][indices[3]+0]; reg34=reg7*vectors[0][indices[3]+1];
    reg36=reg26-reg36; reg25=reg38-reg25; reg27=reg35-reg27; reg4=reg4/reg23; reg31=reg31/reg23;
    reg21=reg21/reg23; reg20=reg20/reg23; reg30=reg33+reg30; reg26=pow(reg2,2); reg33=reg2*reg29;
    reg32=reg10-reg32; reg10=reg3*reg28; reg35=1-(*f.m).resolution; reg36=reg34+reg36; reg34=pow(reg3,2);
    reg37=reg25*reg31; reg38=reg20*reg30; T reg39=reg21*reg36; reg32=reg35*reg32; T reg40=(*f.m).resolution*(*f.m).alpha;
    T reg41=reg4*reg27; reg10=reg33-reg10; reg34=reg26-reg34; reg36=reg20*reg36; reg29=reg29/reg10;
    reg32=reg40+reg32; reg27=reg31*reg27; reg34=reg34/reg10; reg37=reg38-reg37; reg25=reg25*reg4;
    reg10=reg28/reg10; reg39=reg41-reg39; reg30=reg30*reg21; reg26=(*f.m).resolution*reg34; reg28=(*f.m).resolution*reg10;
    reg37=reg39+reg37; reg9=reg9*reg35; reg30=reg25-reg30; reg13=reg35*reg13; reg27=reg36-reg27;
    reg32=(*f.m).deltaT*reg32; reg35=reg6*reg35; reg6=(*f.m).resolution*reg29; reg25=reg1*reg31; reg33=var_inter[0]*reg21;
    reg36=var_inter[0]*reg20; reg37=0.5*reg37; reg13=reg26+reg13; reg30=reg30-reg32; reg9=reg28+reg9;
    reg6=reg35+reg6; reg26=reg7*reg21; reg28=reg7*reg20; reg35=var_inter[1]*reg31; reg38=reg1*reg4;
    reg27=reg27-reg32; reg39=var_inter[1]*reg4; reg40=reg26-reg38; reg41=reg25-reg28; T reg42=reg30*reg6;
    T reg43=reg28+reg35; T reg44=reg25+reg36; T reg45=reg26+reg39; T reg46=reg39-reg33; T reg47=reg38+reg33;
    reg37=reg37*reg13; T reg48=reg36-reg35; T reg49=reg6*reg27; reg30=reg30*reg9; reg27=reg9*reg27;
    T reg50=0.5*reg43; T reg51=0.5*reg45; T reg52=0.5*reg40; reg49=reg30+reg49; reg30=0.5*reg46;
    reg27=reg42+reg27; reg42=0.5*reg48; T reg53=0.5*reg44; T reg54=0.5*reg41; reg37=2*reg37;
    T reg55=0.5*reg47; T reg56=reg45*reg27; T reg57=reg37*reg50; T reg58=reg1*var_inter[0]; T reg59=reg41*reg49;
    T reg60=var_inter[1]*var_inter[0]; T reg61=reg37*reg52; T reg62=reg40*reg27; T reg63=var_inter[1]*reg7; T reg64=reg1*reg7;
    T reg65=reg37*reg54; T reg66=reg37*reg53; T reg67=reg47*reg27; T reg68=reg37*reg55; T reg69=reg44*reg49;
    T reg70=reg43*reg49; T reg71=reg37*reg51; T reg72=reg37*reg42; T reg73=reg48*reg49; T reg74=reg37*reg30;
    T reg75=reg46*reg27; T reg76=reg58*(*f.m).f_vol[1]; reg68=reg68-reg69; reg73=reg74+reg73; reg74=reg60*(*f.m).f_vol[0];
    T reg77=reg58*(*f.m).f_vol[0]; reg67=reg67-reg66; reg75=reg72+reg75; reg72=reg64*(*f.m).f_vol[1]; reg59=reg61+reg59;
    reg57=reg57-reg56; reg61=reg63*(*f.m).f_vol[0]; T reg78=reg60*(*f.m).f_vol[1]; reg70=reg70-reg71; T reg79=reg64*(*f.m).f_vol[0];
    reg62=reg65+reg62; reg65=reg63*(*f.m).f_vol[1]; reg57=reg57-reg61; reg59=reg59-reg72; reg62=reg62-reg79;
    reg75=reg75-reg74; reg67=reg67-reg77; reg70=reg70-reg65; reg73=reg73-reg78; reg68=reg68-reg76;
    reg70=reg70*reg23; reg75=reg75*reg23; reg57=reg57*reg23; reg73=reg73*reg23; reg68=reg68*reg23;
    reg67=reg67*reg23; reg62=reg62*reg23; reg59=reg59*reg23; T vec_3=ponderation*reg68; T vec_0=ponderation*reg62;
    T vec_5=ponderation*reg73; T vec_6=ponderation*reg57; T vec_2=ponderation*reg67; T vec_4=ponderation*reg75; T vec_7=ponderation*reg70;
    T vec_1=ponderation*reg59;
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

