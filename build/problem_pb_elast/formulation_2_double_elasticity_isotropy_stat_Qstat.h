
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
    T reg4=reg0*reg1; T reg5=reg2*reg1; reg1=reg3*reg1; T reg6=reg2*reg4; reg4=reg3*reg4;
    T reg7=reg3*reg5; T reg8=reg3*reg1; reg5=reg2*reg5; T reg9=reg2*reg6; reg6=reg3*reg6;
    T reg10=reg3*reg4; reg7=reg8+reg7; reg5=reg5-reg8; reg9=reg9-reg10; reg6=reg10+reg6;
    reg4=reg2*reg4; reg1=reg2*reg1; T reg11=reg2*reg9; reg5=reg2*reg5; reg4=reg10+reg4;
    reg10=reg3*reg6; reg7=reg3*reg7; T reg12=reg8+reg1; T reg13=elem.pos(2)[1]-elem.pos(0)[1]; T reg14=elem.pos(2)[0]-elem.pos(0)[0];
    reg10=reg11-reg10; reg11=reg3*reg4; reg7=reg5-reg7; reg12=reg3*reg12; reg5=elem.pos(1)[1]-elem.pos(0)[1];
    T reg15=elem.pos(1)[0]-elem.pos(0)[0]; T reg16=reg15*reg13; T reg17=reg5*reg14; reg12=reg7-reg12; reg11=reg10-reg11;
    reg6=reg6/reg11; reg12=reg12/reg11; reg17=reg16-reg17; reg9=reg9/reg11; reg7=vectors[0][indices[2]+1]-vectors[0][indices[0]+1];
    reg10=reg9*reg12; reg16=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; T reg18=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; T reg19=reg2*reg0; T reg20=vectors[0][indices[1]+0]-vectors[0][indices[0]+0];
    reg14=reg14/reg17; T reg21=reg3*reg0; T reg22=reg6*reg12; reg5=reg5/reg17; reg15=reg15/reg17;
    reg13=reg13/reg17; T reg23=reg15*reg7; T reg24=reg14*reg16; T reg25=reg2*reg19; T reg26=(*f.m).alpha*reg9;
    T reg27=reg6*reg22; T reg28=(*f.m).alpha*reg6; T reg29=reg9*reg10; T reg30=reg5*reg18; T reg31=reg13*reg20;
    reg11=reg4/reg11; reg4=reg3*reg21; reg4=reg25-reg4; reg27=reg29-reg27; reg24=reg23-reg24;
    elem.epsilon[0][1]=reg24; reg30=reg31-reg30; elem.epsilon[0][0]=reg30; reg28=reg26+reg28; reg11=(*f.m).alpha*reg11;
    reg23=(*f.m).deltaT*(*f.m).alpha; reg22=reg22/reg27; reg25=reg24-reg23; reg26=reg30-reg23; reg21=reg21/reg4;
    reg11=reg28+reg11; reg27=reg10/reg27; reg19=reg19/reg4; reg10=reg19*reg25; reg27=reg27*reg11;
    reg11=reg22*reg11; reg22=reg21*reg26; reg25=reg21*reg25; reg26=reg19*reg26; reg11=reg27-reg11;
    reg27=1-(*f.m).resolution; reg25=reg26+reg25; reg22=reg10+reg22; reg10=reg2*reg25; reg11=reg27*reg11;
    reg25=reg3*reg25; reg26=reg3*reg22; reg28=(*f.m).alpha*(*f.m).resolution; reg22=reg2*reg22; reg29=PNODE(2).dep[0]-PNODE(0).dep[0];
    reg31=pow(reg3,2); T reg32=PNODE(1).dep[0]-PNODE(0).dep[0]; T reg33=PNODE(2).dep[1]-PNODE(0).dep[1]; T reg34=pow(reg2,2); T reg35=PNODE(1).dep[1]-PNODE(0).dep[1];
    reg11=reg28+reg11; reg28=reg33*reg15; T reg36=reg35*reg14; reg33=reg33*reg5; reg22=reg22-reg25;
    reg35=reg35*reg13; reg10=reg10-reg26; T reg37=reg29*reg15; reg31=reg34-reg31; reg34=reg32*reg14;
    reg29=reg29*reg5; reg32=reg13*reg32; reg34=reg37-reg34; reg29=reg32-reg29; reg7=reg5*reg7;
    reg16=reg13*reg16; reg20=reg20*reg14; reg18=reg15*reg18; reg33=reg35-reg33; reg4=reg31/reg4;
    reg36=reg28-reg36; reg11=(*f.m).deltaT*reg11; reg28=(*f.m).resolution*reg19; reg9=reg9*reg27; reg31=(*f.m).resolution*reg21;
    reg10=reg23+reg10; reg22=reg23+reg22; reg25=reg26+reg25; reg6=reg6*reg27; reg28=reg9+reg28;
    reg9=reg29-reg11; reg7=reg16-reg7; reg6=reg31+reg6; reg33=reg34+reg33; reg16=reg10+reg22;
    reg12=reg27*reg12; reg25=reg23-reg25; reg26=(*f.m).resolution*reg4; reg27=reg36-reg11; reg20=reg18-reg20;
    reg7=reg20+reg7; reg33=0.5*reg33; reg16=reg16+reg25; reg18=reg27*reg28; reg20=reg9*reg6;
    reg31=reg27*reg6; reg32=reg9*reg28; reg12=reg26+reg12; reg18=reg20+reg18; reg20=reg33*reg12;
    reg16=reg16/3; reg31=reg32+reg31; reg7=0.5*reg7; elem.epsilon[0][2]=reg7; reg26=reg7*reg4;
    reg9=reg9*reg31; reg10=reg10-reg16; reg22=reg22-reg16; reg20=2*reg20; reg27=reg27*reg18;
    reg32=reg20*reg33; reg10=pow(reg10,2); reg22=pow(reg22,2); reg27=reg9+reg27; reg16=reg25-reg16;
    reg26=reg0*reg26; reg32=reg27+reg32; reg16=pow(reg16,2); reg22=reg10+reg22; reg9=2*reg26;
    reg16=reg22+reg16; reg9=reg26*reg9; reg32=reg17*reg32; reg24=reg24-reg11; reg9=reg16+reg9;
    reg10=0.33333333333333331483*reg32; reg32=0.16666666666666665741*reg32; reg11=reg30-reg11; reg16=reg28*reg24; reg17=reg11*reg6;
    reg24=reg6*reg24; reg32=reg10+reg32; reg28=reg11*reg28; reg9=1.5*reg9; elem.sigma_von_mises=pow(reg9,0.5);
    elem.ener=reg32/2; elem.sigma[0][2]=reg12*reg7; elem.sigma[0][0]=reg28+reg24; elem.sigma[0][1]=reg17+reg16;
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
    T reg4=reg0*reg1; T reg5=reg3*reg1; T reg6=reg3*reg4; reg4=reg2*reg4; reg1=reg2*reg1;
    T reg7=reg2*reg1; T reg8=reg3*reg5; T reg9=reg3*reg4; T reg10=reg3*reg6; reg4=reg2*reg4;
    reg1=reg3*reg1; reg6=reg2*reg6; reg9=reg10+reg9; reg4=reg4-reg10; reg5=reg2*reg5;
    reg7=reg7-reg8; reg1=reg8+reg1; reg1=reg3*reg1; reg7=reg2*reg7; T reg11=reg2*reg4;
    T reg12=reg8+reg5; reg6=reg10+reg6; reg10=reg3*reg9; reg10=reg11-reg10; reg12=reg3*reg12;
    reg11=reg3*reg6; reg1=reg7-reg1; reg11=reg10-reg11; reg12=reg1-reg12; reg4=reg4/reg11;
    reg9=reg9/reg11; reg12=reg12/reg11; reg1=reg4*reg12; reg7=reg9*reg12; reg10=reg9*reg7;
    reg11=reg6/reg11; reg6=(*f.m).alpha*reg4; T reg13=(*f.m).alpha*reg9; T reg14=reg4*reg1; reg10=reg14-reg10;
    reg13=reg6+reg13; reg11=(*f.m).alpha*reg11; reg1=reg1/reg10; reg10=reg7/reg10; reg11=reg13+reg11;
    reg6=reg3*reg0; reg7=reg2*reg0; reg10=reg10*reg11; reg11=reg1*reg11; reg1=pow(reg3,2);
    reg13=elem.pos(2)[1]-elem.pos(0)[1]; reg14=elem.pos(2)[0]-elem.pos(0)[0]; T reg15=elem.pos(1)[1]-elem.pos(0)[1]; T reg16=elem.pos(1)[0]-elem.pos(0)[0]; T reg17=reg2*reg7;
    T reg18=pow(reg2,2); T reg19=reg3*reg6; T reg20=1-(*f.m).resolution; T reg21=reg15*reg14; T reg22=reg16*reg13;
    reg10=reg11-reg10; reg1=reg18-reg1; reg19=reg17-reg19; reg11=(*f.m).alpha*(*f.m).resolution; reg1=reg1/reg19;
    reg21=reg22-reg21; reg10=reg20*reg10; reg6=reg6/reg19; reg19=reg7/reg19; reg7=(*f.m).resolution*reg1;
    reg17=(*f.m).resolution*reg6; reg15=reg15/reg21; reg16=reg16/reg21; reg13=reg13/reg21; reg12=reg20*reg12;
    reg14=reg14/reg21; reg9=reg9*reg20; reg18=(*f.m).resolution*reg19; reg20=reg4*reg20; reg10=reg11+reg10;
    reg12=reg7+reg12; reg10=(*f.m).deltaT*reg10; reg9=reg17+reg9; reg18=reg20+reg18; reg4=0.5*reg13;
    reg7=0.5*reg14; reg11=0.5*reg15; reg17=reg14-reg16; reg20=reg15-reg13; reg22=0.5*reg16;
    T reg23=reg7*reg12; T reg24=reg11*reg12; T reg25=0.5*reg20; T reg26=reg22*reg12; T reg27=reg4*reg12;
    T reg28=0.5*reg17; T reg29=reg10*reg9; T reg30=reg10*reg18; reg26=2*reg26; T reg31=reg9*reg14;
    T reg32=reg15*reg18; T reg33=reg16*reg9; reg27=2*reg27; T reg34=2*reg24; T reg35=reg16*reg18;
    T reg36=reg15*reg9; T reg37=reg18*reg14; T reg38=2*reg23; T reg39=reg28*reg12; T reg40=reg13*reg9;
    T reg41=reg13*reg18; T reg42=reg30+reg29; T reg43=reg25*reg12; T reg44=1-var_inter[0]; reg44=reg44-var_inter[1];
    T reg45=reg4*reg34; T reg46=reg35*reg14; T reg47=reg40*reg14; T reg48=reg4*reg38; T reg49=reg17*reg18;
    T reg50=reg26*reg11; T reg51=reg42*reg14; T reg52=reg15*reg42; T reg53=reg16*reg37; T reg54=reg27*reg11;
    T reg55=reg16*reg36; T reg56=var_inter[0]*(*f.m).f_vol[1]; T reg57=reg20*reg9; T reg58=reg31*reg13; T reg59=reg7*reg27;
    T reg60=reg32*reg13; T reg61=reg7*reg26; T reg62=var_inter[1]*(*f.m).f_vol[0]; T reg63=reg20*reg18; reg39=2*reg39;
    T reg64=reg17*reg9; reg43=2*reg43; T reg65=reg22*reg38; T reg66=reg41*reg15; T reg67=reg22*reg34;
    T reg68=reg33*reg15; reg68=reg67+reg68; T reg69=reg25*reg38; T reg70=reg17*reg49; T reg71=reg57*reg17;
    T reg72=reg25*reg43; T reg73=reg64*reg13; T reg74=reg7*reg43; T reg75=reg41*reg13; T reg76=reg7*reg38;
    reg59=reg58+reg59; T reg77=reg57*reg16; reg61=reg60+reg61; T reg78=reg33*reg13; T reg79=reg7*reg34;
    T reg80=reg20*reg42; T reg81=reg34*reg11; T reg82=reg16*reg35; reg50=reg55+reg50; T reg83=reg4*reg39;
    reg57=reg57*reg14; reg54=reg53+reg54; T reg84=reg11*reg38; T reg85=reg16*reg40; T reg86=reg11*reg43;
    T reg87=reg16*reg49; T reg88=reg16*reg42; T reg89=reg52-reg62; T reg90=reg20*reg31; T reg91=reg27*reg28;
    T reg92=reg63*reg15; T reg93=reg22*reg39; reg46=reg45+reg46; T reg94=reg20*reg32; T reg95=reg36*reg14;
    T reg96=reg4*reg26; T reg97=reg37*reg14; T reg98=reg4*reg27; reg47=reg48+reg47; T reg99=reg26*reg28;
    reg49=reg49*reg14; T reg100=reg4*reg43; reg33=reg20*reg33; T reg101=reg34*reg28; T reg102=reg7*reg39;
    T reg103=reg63*reg13; T reg104=reg22*reg26; reg35=reg17*reg35; T reg105=reg25*reg34; T reg106=reg17*reg36;
    reg26=reg25*reg26; T reg107=reg17*reg37; T reg108=reg25*reg39; T reg109=reg25*reg27; reg40=reg17*reg40;
    T reg110=reg64*reg15; T reg111=reg13*reg42; T reg112=reg44*(*f.m).f_vol[1]; T reg113=reg44*(*f.m).f_vol[0]; reg63=reg20*reg63;
    T reg114=reg22*reg43; reg27=reg22*reg27; T reg115=reg31*reg15; T reg116=reg28*reg39; T reg117=var_inter[0]*(*f.m).f_vol[0];
    reg39=reg11*reg39; T reg118=reg51-reg56; T reg119=var_inter[1]*(*f.m).f_vol[1]; reg66=reg65+reg66; T reg120=reg28*reg38;
    reg41=reg20*reg41; reg64=reg20*reg64; T reg121=reg32*reg15; reg43=reg28*reg43; T reg122=reg17*reg42;
    T reg123=reg59*reg21; reg33=reg33-reg101; reg102=reg103-reg102; reg103=reg113+reg80; T reg124=reg68*reg21;
    T reg125=reg61*reg21; reg49=reg100-reg49; reg110=reg114-reg110; reg75=reg75+reg76; reg39=reg77-reg39;
    reg74=reg73-reg74; reg35=reg35-reg105; reg116=reg63+reg116; reg41=reg41-reg120; reg71=reg108+reg71;
    reg26=reg26-reg106; reg70=reg72+reg70; reg109=reg109-reg107; reg40=reg40-reg69; reg43=reg64+reg43;
    reg104=reg121+reg104; reg89=reg21*reg89; reg91=reg91-reg90; reg63=reg119+reg88; reg92=reg93-reg92;
    reg118=reg21*reg118; reg86=reg87-reg86; reg27=reg27+reg115; reg85=reg85+reg84; reg64=reg21*reg46;
    reg57=reg83-reg57; reg72=reg21*reg54; reg96=reg96+reg95; reg73=reg117+reg111; reg77=reg21*reg50;
    reg98=reg98+reg97; reg99=reg99-reg94; reg82=reg82+reg81; reg83=reg112+reg122; reg87=reg21*reg47;
    reg93=reg66*reg21; reg78=reg78+reg79; reg49=reg21*reg49; reg43=reg43*reg21; reg41=reg41*reg21;
    reg104=reg104*reg21; reg109=reg21*reg109; reg91=reg91*reg21; reg92=reg21*reg92; reg99=reg99*reg21;
    reg110=reg110*reg21; reg26=reg21*reg26; reg27=reg27*reg21; reg33=reg33*reg21; reg100=ponderation*reg64;
    reg108=ponderation*reg87; reg35=reg21*reg35; reg114=ponderation*reg93; reg96=reg21*reg96; reg98=reg21*reg98;
    T reg126=reg21*reg103; T reg127=ponderation*reg125; reg57=reg57*reg21; reg85=reg21*reg85; T reg128=ponderation*reg123;
    T reg129=ponderation*reg124; T reg130=ponderation*reg72; reg75=reg75*reg21; reg86=reg21*reg86; reg74=reg74*reg21;
    reg39=reg39*reg21; reg118=ponderation*reg118; reg102=reg21*reg102; reg78=reg78*reg21; T reg131=reg21*reg73;
    reg71=reg71*reg21; reg116=reg116*reg21; T reg132=reg21*reg63; reg70=reg21*reg70; reg82=reg21*reg82;
    T reg133=reg21*reg83; reg89=ponderation*reg89; reg40=reg21*reg40; T reg134=ponderation*reg77; T tmp_5_3=-reg130;
    T tmp_4_2=-reg114; T tmp_3_3=ponderation*reg98; T tmp_3_4=ponderation*reg96; reg96=ponderation*reg131; sollicitation[indices[1]+0]+=reg96;
    T tmp_0_3=ponderation*reg91; T tmp_5_2=ponderation*reg85; T tmp_3_5=-reg100; T tmp_5_1=ponderation*reg86; reg85=ponderation*reg132;
    sollicitation[indices[2]+1]+=reg85; T tmp_4_0=ponderation*reg92; sollicitation[indices[1]+1]+=-reg118; T tmp_4_3=ponderation*reg27; sollicitation[indices[2]+0]+=-reg89;
    T tmp_0_2=ponderation*reg41; T tmp_3_0=ponderation*reg57; T tmp_1_2=ponderation*reg40; T tmp_0_1=ponderation*reg43; T tmp_5_0=ponderation*reg39;
    T tmp_1_1=ponderation*reg70; T tmp_1_3=ponderation*reg109; T tmp_0_5=ponderation*reg33; T tmp_0_0=ponderation*reg116; T tmp_1_4=ponderation*reg26;
    T tmp_1_0=ponderation*reg71; T tmp_2_0=ponderation*reg102; T tmp_1_5=ponderation*reg35; T tmp_2_1=ponderation*reg74; T tmp_2_2=ponderation*reg75;
    T tmp_2_3=-reg128; T tmp_4_1=ponderation*reg110; T tmp_0_4=ponderation*reg99; T tmp_4_5=-reg129; T tmp_4_4=ponderation*reg104;
    T tmp_2_4=-reg127; T tmp_3_1=ponderation*reg49; reg26=ponderation*reg126; sollicitation[indices[0]+0]+=reg26; T tmp_2_5=ponderation*reg78;
    T tmp_5_5=ponderation*reg82; T tmp_3_2=-reg108; reg27=ponderation*reg133; sollicitation[indices[0]+1]+=reg27; T tmp_5_4=-reg134;
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
    T reg4=reg0*reg1; T reg5=reg2*reg4; reg4=reg3*reg4; T reg6=reg2*reg1; reg1=reg3*reg1;
    T reg7=reg2*reg5; T reg8=reg3*reg6; T reg9=reg3*reg1; reg6=reg2*reg6; T reg10=reg3*reg4;
    reg5=reg3*reg5; reg8=reg9+reg8; reg7=reg7-reg10; reg5=reg10+reg5; reg4=reg2*reg4;
    reg1=reg2*reg1; reg6=reg6-reg9; T reg11=reg3*reg5; T reg12=reg9+reg1; reg4=reg10+reg4;
    reg10=reg2*reg7; reg8=reg3*reg8; reg6=reg2*reg6; reg12=reg3*reg12; reg8=reg6-reg8;
    reg6=reg3*reg4; reg11=reg10-reg11; reg6=reg11-reg6; reg12=reg8-reg12; reg12=reg12/reg6;
    reg7=reg7/reg6; reg5=reg5/reg6; reg8=reg7*reg12; reg10=reg5*reg12; reg11=(*f.m).alpha*reg5;
    T reg13=(*f.m).alpha*reg7; T reg14=reg5*reg10; reg6=reg4/reg6; reg4=reg7*reg8; reg14=reg4-reg14;
    reg11=reg13+reg11; reg6=(*f.m).alpha*reg6; reg10=reg10/reg14; reg14=reg8/reg14; reg4=reg3*reg0;
    reg8=reg2*reg0; reg6=reg11+reg6; reg10=reg10*reg6; reg6=reg14*reg6; reg11=elem.pos(2)[1]-elem.pos(0)[1];
    reg13=elem.pos(2)[0]-elem.pos(0)[0]; reg14=elem.pos(1)[1]-elem.pos(0)[1]; T reg15=elem.pos(1)[0]-elem.pos(0)[0]; T reg16=pow(reg2,2); T reg17=pow(reg3,2);
    T reg18=reg3*reg4; T reg19=reg2*reg8; T reg20=1-(*f.m).resolution; reg10=reg6-reg10; reg18=reg19-reg18;
    reg6=reg14*reg13; reg19=reg15*reg11; reg17=reg16-reg17; reg8=reg8/reg18; reg6=reg19-reg6;
    reg4=reg4/reg18; reg10=reg20*reg10; reg18=reg17/reg18; reg16=(*f.m).alpha*(*f.m).resolution; reg17=(*f.m).resolution*reg18;
    reg19=(*f.m).resolution*reg4; reg10=reg16+reg10; reg12=reg20*reg12; reg5=reg5*reg20; reg14=reg14/reg6;
    reg13=reg13/reg6; reg11=reg11/reg6; reg20=reg7*reg20; reg7=(*f.m).resolution*reg8; reg15=reg15/reg6;
    reg16=0.5*reg11; T reg21=0.5*reg15; reg10=(*f.m).deltaT*reg10; T reg22=0.5*reg14; T reg23=reg14-reg11;
    reg12=reg17+reg12; reg5=reg19+reg5; reg7=reg20+reg7; reg17=reg13-reg15; reg19=reg21*reg12;
    reg20=0.5*reg23; T reg24=0.5*reg17; T reg25=0.5*reg13; T reg26=reg16*reg12; T reg27=reg10*reg5;
    T reg28=reg10*reg7; T reg29=reg22*reg12; T reg30=reg28+reg27; T reg31=reg15*reg7; reg26=2*reg26;
    T reg32=reg5*reg13; T reg33=1-var_inter[0]; T reg34=reg25*reg12; T reg35=reg14*reg7; T reg36=reg24*reg12;
    reg19=2*reg19; T reg37=2*reg29; T reg38=reg20*reg12; T reg39=reg15*reg5; T reg40=reg14*reg5;
    T reg41=reg11*reg7; T reg42=reg7*reg13; reg38=2*reg38; T reg43=reg17*reg5; T reg44=reg11*reg5;
    T reg45=reg17*reg7; T reg46=reg32*reg11; T reg47=reg25*reg26; T reg48=reg35*reg11; T reg49=reg25*reg19;
    T reg50=reg30*reg13; T reg51=var_inter[1]*(*f.m).f_vol[0]; T reg52=reg14*reg30; T reg53=reg23*reg7; reg36=2*reg36;
    T reg54=2*reg34; T reg55=reg16*reg37; T reg56=reg31*reg13; reg33=reg33-var_inter[1]; T reg57=var_inter[0]*(*f.m).f_vol[1];
    T reg58=reg39*reg14; T reg59=reg21*reg37; reg47=reg46+reg47; reg49=reg48+reg49; reg56=reg55+reg56;
    T reg60=reg39*reg11; T reg61=reg25*reg37; T reg62=var_inter[1]*(*f.m).f_vol[1]; T reg63=reg19*reg24; T reg64=reg23*reg53;
    T reg65=reg23*reg35; T reg66=reg23*reg30; T reg67=reg24*reg36; T reg68=reg37*reg22; T reg69=reg15*reg31;
    T reg70=reg17*reg40; T reg71=reg20*reg37; reg31=reg17*reg31; T reg72=reg20*reg19; T reg73=reg17*reg42;
    T reg74=reg21*reg19; T reg75=reg20*reg26; T reg76=reg37*reg24; reg39=reg23*reg39; T reg77=reg17*reg44;
    T reg78=reg16*reg26; T reg79=reg42*reg13; T reg80=reg20*reg54; T reg81=reg17*reg45; T reg82=reg20*reg38;
    T reg83=reg41*reg11; T reg84=reg25*reg54; T reg85=reg16*reg19; T reg86=reg40*reg13; T reg87=reg11*reg30;
    T reg88=reg17*reg30; T reg89=var_inter[0]*(*f.m).f_vol[0]; T reg90=reg23*reg32; T reg91=reg35*reg14; T reg92=reg50-reg57;
    T reg93=reg23*reg43; T reg94=reg24*reg38; T reg95=reg26*reg24; reg58=reg59+reg58; T reg96=reg52-reg51;
    T reg97=reg33*(*f.m).f_vol[1]; T reg98=reg23*reg41; T reg99=reg15*reg30; T reg100=reg24*reg54; T reg101=reg33*(*f.m).f_vol[0];
    T reg102=reg6*reg56; reg83=reg83+reg84; T reg103=reg97+reg88; reg31=reg31-reg71; reg39=reg39-reg76;
    reg78=reg78+reg79; T reg104=reg101+reg66; reg81=reg82+reg81; reg98=reg98-reg100; reg95=reg95-reg90;
    reg72=reg72-reg70; reg77=reg77-reg80; reg94=reg93+reg94; reg85=reg85+reg86; reg75=reg75-reg73;
    reg96=reg6*reg96; reg74=reg91+reg74; reg82=reg47*reg6; reg93=reg89+reg87; reg69=reg69+reg68;
    T reg105=reg58*reg6; T reg106=reg62+reg99; T reg107=reg49*reg6; reg60=reg60+reg61; reg92=reg6*reg92;
    reg67=reg64+reg67; reg63=reg63-reg65; reg75=reg6*reg75; reg96=ponderation*reg96; reg98=reg98*reg6;
    reg72=reg6*reg72; reg74=reg74*reg6; reg64=ponderation*reg105; reg69=reg6*reg69; T reg108=reg6*reg106;
    reg31=reg6*reg31; reg78=reg6*reg78; reg95=reg95*reg6; reg94=reg94*reg6; reg92=ponderation*reg92;
    reg60=reg60*reg6; reg85=reg6*reg85; reg77=reg6*reg77; reg67=reg67*reg6; T reg109=ponderation*reg107;
    reg81=reg6*reg81; reg63=reg63*reg6; T reg110=reg6*reg104; reg39=reg39*reg6; T reg111=reg6*reg93;
    T reg112=ponderation*reg82; T reg113=ponderation*reg102; reg83=reg83*reg6; T reg114=reg6*reg103; T tmp_3_3=ponderation*reg78;
    T tmp_3_4=ponderation*reg85; T tmp_0_2=ponderation*reg98; T tmp_4_4=ponderation*reg74; T tmp_3_5=-reg113; T tmp_4_5=-reg64;
    reg64=ponderation*reg108; sollicitation[indices[2]+1]+=reg64; T tmp_5_5=ponderation*reg69; sollicitation[indices[2]+0]+=-reg96; sollicitation[indices[1]+1]+=-reg92;
    T tmp_2_5=ponderation*reg60; T tmp_2_4=-reg109; reg60=ponderation*reg111; sollicitation[indices[1]+0]+=reg60; T tmp_2_3=-reg112;
    reg69=ponderation*reg114; sollicitation[indices[0]+1]+=reg69; T tmp_2_2=ponderation*reg83; reg74=ponderation*reg110; sollicitation[indices[0]+0]+=reg74;
    T tmp_0_5=ponderation*reg39; T tmp_1_1=ponderation*reg81; T tmp_0_0=ponderation*reg67; T tmp_1_2=ponderation*reg77; T tmp_1_3=ponderation*reg75;
    T tmp_0_1=ponderation*reg94; T tmp_1_4=ponderation*reg72; T tmp_1_5=ponderation*reg31; T tmp_0_3=ponderation*reg95; T tmp_0_4=ponderation*reg63;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg2*reg4; T reg6=reg2*reg1; reg1=reg3*reg1; reg4=reg3*reg4;
    T reg7=reg2*reg5; T reg8=reg3*reg6; T reg9=reg3*reg1; reg6=reg2*reg6; T reg10=reg3*reg4;
    reg5=reg3*reg5; reg8=reg9+reg8; reg7=reg7-reg10; reg5=reg10+reg5; reg4=reg2*reg4;
    reg1=reg2*reg1; reg6=reg6-reg9; T reg11=reg2*reg7; reg4=reg10+reg4; reg6=reg2*reg6;
    reg8=reg3*reg8; reg10=reg3*reg5; T reg12=reg9+reg1; T reg13=reg3*reg0; T reg14=reg2*reg0;
    T reg15=pow(reg2,2); T reg16=reg3*reg13; reg10=reg11-reg10; reg11=reg3*reg4; T reg17=pow(reg3,2);
    reg8=reg6-reg8; reg6=reg2*reg14; reg12=reg3*reg12; T reg18=elem.pos(2)[1]-elem.pos(0)[1]; T reg19=elem.pos(2)[0]-elem.pos(0)[0];
    T reg20=elem.pos(1)[0]-elem.pos(0)[0]; T reg21=elem.pos(1)[1]-elem.pos(0)[1]; reg17=reg15-reg17; reg11=reg10-reg11; reg16=reg6-reg16;
    reg6=reg21*reg19; reg12=reg8-reg12; reg8=reg20*reg18; reg12=reg12/reg11; reg17=reg17/reg16;
    reg10=1-(*f.m).resolution; reg6=reg8-reg6; reg8=(*f.m).resolution*reg17; reg14=reg14/reg16; reg7=reg7/reg11;
    reg15=reg10*reg12; reg5=reg5/reg11; reg19=reg19/reg6; reg16=reg13/reg16; reg21=reg21/reg6;
    reg18=reg18/reg6; reg20=reg20/reg6; reg13=reg19-reg20; T reg22=0.5*reg21; T reg23=reg21-reg18;
    T reg24=0.5*reg19; T reg25=0.5*reg18; T reg26=0.5*reg20; reg15=reg8+reg15; reg8=reg5*reg10;
    T reg27=(*f.m).resolution*reg14; T reg28=reg7*reg10; T reg29=(*f.m).resolution*reg16; reg27=reg28+reg27; reg8=reg29+reg8;
    reg28=0.5*reg13; reg29=reg25*reg15; T reg30=reg24*reg15; T reg31=reg26*reg15; T reg32=reg22*reg15;
    T reg33=0.5*reg23; T reg34=reg18*reg27; T reg35=2*reg30; T reg36=reg20*reg27; T reg37=reg21*reg8;
    T reg38=reg27*reg19; T reg39=reg18*reg8; reg31=2*reg31; T reg40=reg21*reg27; reg29=2*reg29;
    T reg41=reg8*reg19; T reg42=reg20*reg8; T reg43=2*reg32; T reg44=reg33*reg15; T reg45=reg28*reg15;
    T reg46=reg36*reg19; T reg47=reg25*reg43; T reg48=reg31*reg22; T reg49=reg20*reg37; T reg50=reg23*reg8;
    T reg51=reg13*reg8; T reg52=reg29*reg22; T reg53=reg20*reg38; reg44=2*reg44; T reg54=reg26*reg43;
    T reg55=reg42*reg21; T reg56=reg39*reg19; T reg57=reg25*reg35; T reg58=reg34*reg21; T reg59=reg26*reg35;
    T reg60=reg13*reg27; reg45=2*reg45; T reg61=reg24*reg29; T reg62=reg41*reg18; T reg63=reg23*reg27;
    T reg64=reg24*reg31; T reg65=reg40*reg18; T reg66=reg20*reg60; T reg67=reg33*reg31; T reg68=reg23*reg63;
    T reg69=reg13*reg37; T reg70=reg42*reg18; T reg71=reg24*reg43; T reg72=reg26*reg31; T reg73=reg13*reg36;
    T reg74=reg23*reg41; T reg75=reg33*reg43; T reg76=reg50*reg19; T reg77=reg25*reg45; reg61=reg62+reg61;
    reg64=reg65+reg64; T reg78=reg24*reg35; T reg79=reg34*reg18; T reg80=reg43*reg22; reg36=reg20*reg36;
    T reg81=reg24*reg44; T reg82=reg51*reg18; reg48=reg49+reg48; T reg83=reg13*reg60; T reg84=reg33*reg35;
    reg52=reg53+reg52; T reg85=reg33*reg44; T reg86=reg13*reg39; T reg87=reg22*reg35; reg39=reg20*reg39;
    T reg88=reg33*reg29; T reg89=reg13*reg38; T reg90=reg22*reg44; reg55=reg54+reg55; T reg91=reg25*reg29;
    T reg92=reg38*reg19; T reg93=reg40*reg21; T reg94=reg25*reg31; T reg95=reg37*reg19; T reg96=reg41*reg21;
    T reg97=reg26*reg29; reg58=reg59+reg58; reg46=reg47+reg46; T reg98=reg26*reg45; T reg99=reg51*reg21;
    T reg100=reg26*reg44; reg31=reg31*reg28; T reg101=reg63*reg21; T reg102=reg43*reg28; T reg103=reg28*reg35;
    reg42=reg23*reg42; reg34=reg23*reg34; reg29=reg29*reg28; reg63=reg63*reg18; T reg104=reg24*reg45;
    T reg105=reg23*reg40; T reg106=reg28*reg45; T reg107=reg50*reg13; T reg108=reg33*reg45; reg51=reg23*reg51;
    T reg109=reg28*reg44; reg56=reg57+reg56; reg50=reg50*reg20; reg45=reg22*reg45; reg44=reg25*reg44;
    reg60=reg60*reg19; reg83=reg85+reg83; reg85=reg6*reg56; reg79=reg79+reg78; reg73=reg73-reg75;
    reg42=reg42-reg102; reg101=reg98-reg101; reg104=reg63-reg104; reg86=reg86-reg84; reg81=reg82-reg81;
    reg88=reg88-reg89; reg107=reg108+reg107; reg63=reg6*reg46; reg60=reg44-reg60; reg67=reg67-reg69;
    reg94=reg94+reg95; reg91=reg91+reg92; reg44=reg58*reg6; reg70=reg70+reg71; reg90=reg66-reg90;
    reg76=reg77-reg76; reg39=reg39+reg87; reg29=reg29-reg74; reg99=reg100-reg99; reg66=reg6*reg52;
    reg31=reg31-reg105; reg97=reg97+reg96; reg106=reg68+reg106; reg68=reg6*reg48; reg77=reg61*reg6;
    reg82=reg55*reg6; reg72=reg93+reg72; reg109=reg51+reg109; reg51=reg64*reg6; reg34=reg34-reg103;
    reg45=reg50-reg45; reg36=reg36+reg80; reg50=ponderation*reg85; reg72=reg72*reg6; reg67=reg6*reg67;
    reg98=ponderation*reg44; reg70=reg70*reg6; reg91=reg6*reg91; reg76=reg76*reg6; reg73=reg6*reg73;
    reg94=reg6*reg94; reg100=ponderation*reg82; reg42=reg42*reg6; reg97=reg97*reg6; reg106=reg106*reg6;
    reg45=reg45*reg6; reg29=reg29*reg6; reg60=reg6*reg60; reg34=reg34*reg6; reg108=ponderation*reg77;
    reg109=reg109*reg6; reg79=reg79*reg6; reg81=reg81*reg6; reg36=reg6*reg36; T reg110=ponderation*reg68;
    reg101=reg6*reg101; reg104=reg6*reg104; reg31=reg31*reg6; T reg111=ponderation*reg66; reg83=reg6*reg83;
    reg39=reg6*reg39; reg86=reg6*reg86; reg107=reg107*reg6; T reg112=ponderation*reg63; T reg113=ponderation*reg51;
    reg99=reg99*reg6; reg88=reg6*reg88; reg90=reg6*reg90; T tmp_0_1=ponderation*reg109; T tmp_4_1=ponderation*reg99;
    T tmp_4_3=ponderation*reg97; T tmp_0_4=ponderation*reg31; T tmp_3_4=ponderation*reg94; T tmp_4_0=ponderation*reg101; T tmp_3_3=ponderation*reg91;
    T tmp_3_5=-reg112; T tmp_4_2=-reg98; T tmp_0_2=ponderation*reg34; T tmp_2_2=ponderation*reg79; T tmp_2_3=-reg108;
    T tmp_5_5=ponderation*reg36; T tmp_2_1=ponderation*reg81; T tmp_5_4=-reg110; T tmp_2_0=ponderation*reg104; T tmp_5_3=-reg111;
    T tmp_1_1=ponderation*reg83; T tmp_5_2=ponderation*reg39; T tmp_1_2=ponderation*reg86; T tmp_5_1=ponderation*reg90; T tmp_1_0=ponderation*reg107;
    T tmp_1_3=ponderation*reg88; T tmp_2_4=-reg113; T tmp_1_4=ponderation*reg67; T tmp_2_5=ponderation*reg70; T tmp_1_5=ponderation*reg73;
    T tmp_3_0=ponderation*reg76; T tmp_0_3=ponderation*reg29; T tmp_0_0=ponderation*reg106; T tmp_0_5=ponderation*reg42; T tmp_5_0=ponderation*reg45;
    T tmp_3_1=ponderation*reg60; T tmp_4_5=-reg100; T tmp_4_4=ponderation*reg72; T tmp_3_2=-reg50;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg2*reg4; T reg6=reg2*reg1; reg1=reg3*reg1; reg4=reg3*reg4;
    T reg7=reg2*reg5; T reg8=reg3*reg6; T reg9=reg3*reg1; reg6=reg2*reg6; T reg10=reg3*reg4;
    reg5=reg3*reg5; reg8=reg9+reg8; reg7=reg7-reg10; reg5=reg10+reg5; reg4=reg2*reg4;
    reg1=reg2*reg1; reg6=reg6-reg9; T reg11=reg2*reg7; reg4=reg10+reg4; reg6=reg2*reg6;
    reg8=reg3*reg8; reg10=reg3*reg5; T reg12=reg9+reg1; T reg13=reg3*reg0; T reg14=reg2*reg0;
    T reg15=pow(reg3,2); T reg16=reg3*reg4; reg10=reg11-reg10; reg11=reg3*reg13; T reg17=pow(reg2,2);
    reg8=reg6-reg8; reg6=reg2*reg14; reg12=reg3*reg12; T reg18=elem.pos(1)[0]-elem.pos(0)[0]; T reg19=elem.pos(1)[1]-elem.pos(0)[1];
    T reg20=elem.pos(2)[0]-elem.pos(0)[0]; T reg21=elem.pos(2)[1]-elem.pos(0)[1]; reg12=reg8-reg12; reg11=reg6-reg11; reg16=reg10-reg16;
    reg15=reg17-reg15; reg6=reg19*reg20; reg8=reg18*reg21; reg12=reg12/reg16; reg10=1-(*f.m).resolution;
    reg6=reg8-reg6; reg15=reg15/reg11; reg8=reg10*reg12; reg20=reg20/reg6; reg17=(*f.m).resolution*reg15;
    reg14=reg14/reg11; reg7=reg7/reg16; reg5=reg5/reg16; reg11=reg13/reg11; reg19=reg19/reg6;
    reg18=reg18/reg6; reg21=reg21/reg6; reg13=(*f.m).resolution*reg11; T reg22=reg7*reg10; T reg23=(*f.m).resolution*reg14;
    T reg24=reg5*reg10; T reg25=reg19-reg21; reg8=reg17+reg8; reg17=0.5*reg18; T reg26=0.5*reg21;
    T reg27=reg20-reg18; T reg28=0.5*reg19; T reg29=0.5*reg25; T reg30=reg28*reg8; T reg31=reg17*reg8;
    T reg32=reg26*reg8; T reg33=0.5*reg27; T reg34=0.5*reg20; reg24=reg13+reg24; reg23=reg22+reg23;
    reg32=2*reg32; reg13=reg19*reg23; reg31=2*reg31; reg22=reg24*reg20; T reg35=reg18*reg23;
    T reg36=reg34*reg8; T reg37=reg29*reg8; T reg38=reg33*reg8; T reg39=reg18*reg24; T reg40=2*reg30;
    T reg41=reg22*reg21; T reg42=reg13*reg21; T reg43=2*reg36; T reg44=reg21*reg23; T reg45=reg25*reg23;
    T reg46=reg27*reg23; reg38=2*reg38; T reg47=reg35*reg20; T reg48=reg34*reg31; T reg49=reg26*reg40;
    T reg50=reg34*reg32; T reg51=reg21*reg24; T reg52=reg27*reg24; T reg53=reg23*reg20; T reg54=reg39*reg19;
    T reg55=reg17*reg40; T reg56=reg19*reg24; reg37=2*reg37; T reg57=reg33*reg38; T reg58=reg25*reg52;
    T reg59=reg25*reg45; T reg60=reg40*reg33; T reg61=reg34*reg40; T reg62=reg17*reg31; T reg63=reg39*reg21;
    reg39=reg25*reg39; T reg64=reg34*reg43; T reg65=reg44*reg21; T reg66=reg29*reg37; T reg67=reg27*reg46;
    T reg68=reg27*reg35; T reg69=reg29*reg43; T reg70=reg27*reg51; T reg71=reg29*reg40; T reg72=reg29*reg32;
    T reg73=reg27*reg53; T reg74=reg29*reg31; T reg75=reg27*reg56; reg48=reg42+reg48; T reg76=reg13*reg19;
    T reg77=reg32*reg33; T reg78=reg53*reg20; T reg79=reg25*reg44; T reg80=reg25*reg22; T reg81=reg31*reg33;
    T reg82=reg26*reg32; reg47=reg49+reg47; T reg83=reg26*reg31; T reg84=reg56*reg20; reg50=reg41+reg50;
    T reg85=reg33*reg43; T reg86=reg25*reg13; reg54=reg55+reg54; T reg87=reg33*reg37; reg35=reg18*reg35;
    T reg88=reg40*reg28; reg81=reg81-reg86; reg68=reg68-reg71; reg67=reg66+reg67; reg70=reg70-reg69;
    reg62=reg76+reg62; reg66=reg48*reg6; T reg89=reg54*reg6; reg72=reg72-reg73; reg74=reg74-reg75;
    reg83=reg83+reg84; reg79=reg79-reg85; reg82=reg82+reg78; reg65=reg65+reg64; reg35=reg35+reg88;
    reg77=reg77-reg80; reg63=reg63+reg61; T reg90=reg6*reg47; reg39=reg39-reg60; T reg91=reg50*reg6;
    reg57=reg59+reg57; reg87=reg58+reg87; reg83=reg6*reg83; reg58=ponderation*reg90; reg74=reg6*reg74;
    reg39=reg39*reg6; reg68=reg6*reg68; reg82=reg6*reg82; reg59=ponderation*reg89; reg62=reg62*reg6;
    reg81=reg81*reg6; reg77=reg77*reg6; T reg92=ponderation*reg91; reg35=reg6*reg35; reg87=reg87*reg6;
    reg79=reg79*reg6; reg57=reg57*reg6; reg65=reg65*reg6; reg67=reg6*reg67; reg70=reg6*reg70;
    reg63=reg63*reg6; reg72=reg6*reg72; T reg93=ponderation*reg66; T tmp_3_5=-reg58; T tmp_0_0=ponderation*reg57;
    T tmp_0_1=ponderation*reg87; T tmp_2_5=ponderation*reg63; T tmp_0_4=ponderation*reg81; T tmp_2_4=-reg93; T tmp_3_4=ponderation*reg83;
    T tmp_4_4=ponderation*reg62; T tmp_4_5=-reg59; T tmp_3_3=ponderation*reg82; T tmp_0_3=ponderation*reg77; T tmp_2_3=-reg92;
    T tmp_5_5=ponderation*reg35; T tmp_0_2=ponderation*reg79; T tmp_0_5=ponderation*reg39; T tmp_2_2=ponderation*reg65; T tmp_1_5=ponderation*reg68;
    T tmp_1_1=ponderation*reg67; T tmp_1_2=ponderation*reg70; T tmp_1_4=ponderation*reg74; T tmp_1_3=ponderation*reg72;
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
    T reg4=reg0*reg1; T reg5=reg3*reg4; reg4=reg2*reg4; T reg6=reg2*reg1; reg1=reg3*reg1;
    T reg7=reg3*reg4; T reg8=reg3*reg5; reg4=reg2*reg4; T reg9=reg3*reg6; T reg10=reg3*reg1;
    reg6=reg2*reg6; reg6=reg6-reg10; reg9=reg10+reg9; reg1=reg2*reg1; reg5=reg2*reg5;
    reg7=reg8+reg7; reg4=reg4-reg8; T reg11=reg2*reg4; reg5=reg8+reg5; reg6=reg2*reg6;
    reg9=reg3*reg9; reg8=reg3*reg7; T reg12=reg10+reg1; T reg13=reg3*reg5; reg9=reg6-reg9;
    reg8=reg11-reg8; reg12=reg3*reg12; reg13=reg8-reg13; reg12=reg9-reg12; reg4=reg4/reg13;
    reg12=reg12/reg13; reg7=reg7/reg13; reg6=reg4*reg12; reg8=reg7*reg12; reg13=reg5/reg13;
    reg5=reg7*reg8; reg9=reg4*reg6; reg11=(*f.m).alpha*reg7; T reg14=(*f.m).alpha*reg4; reg11=reg14+reg11;
    reg13=(*f.m).alpha*reg13; reg5=reg9-reg5; reg6=reg6/reg5; reg5=reg8/reg5; reg8=reg2*reg0;
    reg9=reg3*reg0; reg13=reg11+reg13; reg11=reg2*reg8; reg5=reg5*reg13; reg13=reg6*reg13;
    reg6=reg3*reg9; reg6=reg11-reg6; reg11=1-(*f.m).resolution; reg5=reg13-reg5; reg13=(*f.m).alpha*(*f.m).resolution;
    reg5=reg11*reg5; reg8=reg8/reg6; reg9=reg9/reg6; reg14=elem.pos(1)[1]-elem.pos(0)[1]; reg4=reg4*reg11;
    T reg15=elem.pos(1)[0]-elem.pos(0)[0]; T reg16=(*f.m).resolution*reg8; reg7=reg7*reg11; T reg17=(*f.m).resolution*reg9; T reg18=elem.pos(2)[0]-elem.pos(0)[0];
    T reg19=elem.pos(2)[1]-elem.pos(0)[1]; reg5=reg13+reg5; reg13=reg14*reg18; T reg20=reg15*reg19; reg5=(*f.m).deltaT*reg5;
    reg7=reg17+reg7; reg16=reg4+reg16; reg13=reg20-reg13; reg4=reg5*reg7; reg17=reg5*reg16;
    reg14=reg14/reg13; reg18=reg18/reg13; reg20=reg17+reg4; reg19=reg19/reg13; T reg21=1-var_inter[0];
    reg15=reg15/reg13; reg21=reg21-var_inter[1]; T reg22=reg14-reg19; T reg23=var_inter[0]*(*f.m).f_vol[1]; T reg24=reg18-reg15;
    T reg25=var_inter[1]*(*f.m).f_vol[0]; T reg26=reg14*reg20; T reg27=reg20*reg18; T reg28=reg21*(*f.m).f_vol[0]; T reg29=reg21*(*f.m).f_vol[1];
    T reg30=var_inter[0]*(*f.m).f_vol[0]; T reg31=var_inter[1]*(*f.m).f_vol[1]; T reg32=reg24*reg20; T reg33=reg19*reg20; T reg34=reg22*reg20;
    T reg35=reg15*reg20; T reg36=reg26-reg25; T reg37=reg27-reg23; T reg38=reg30+reg33; T reg39=reg31+reg35;
    T reg40=reg29+reg32; reg37=reg13*reg37; reg36=reg13*reg36; T reg41=reg28+reg34; T reg42=reg13*reg39;
    reg36=ponderation*reg36; T reg43=reg13*reg41; T reg44=reg13*reg40; reg37=ponderation*reg37; T reg45=reg13*reg38;
    T reg46=ponderation*reg45; sollicitation[indices[1]+0]+=reg46; T reg47=ponderation*reg43; sollicitation[indices[0]+0]+=reg47; sollicitation[indices[2]+0]+=-reg36;
    sollicitation[indices[1]+1]+=-reg37; reg36=ponderation*reg44; sollicitation[indices[0]+1]+=reg36; reg37=ponderation*reg42; sollicitation[indices[2]+1]+=reg37;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg3*reg2; reg2=reg4*reg2; T reg6=reg4*reg1; reg1=reg3*reg1;
    T reg7=reg3*reg5; T reg8=reg4*reg2; reg2=reg3*reg2; T reg9=reg3*reg6; T reg10=reg3*reg1;
    reg6=reg4*reg6; reg1=reg4*reg1; reg9=reg10+reg9; reg6=reg6-reg10; reg5=reg4*reg5;
    reg2=reg7+reg2; reg8=reg8-reg7; reg5=reg7+reg5; reg7=reg10+reg1; T reg11=reg3*reg2;
    T reg12=reg4*reg8; reg6=reg4*reg6; reg9=reg3*reg9; reg11=reg12-reg11; reg12=reg3*reg5;
    reg9=reg6-reg9; reg7=reg3*reg7; reg7=reg9-reg7; reg12=reg11-reg12; reg7=reg7/reg12;
    reg8=reg8/reg12; reg2=reg2/reg12; reg6=reg2*reg7; reg9=reg8*reg7; reg11=(*f.m).alpha*reg2;
    T reg13=(*f.m).alpha*reg8; T reg14=reg2*reg6; T reg15=reg8*reg9; reg12=reg5/reg12; reg5=elem.pos(2)[1]-elem.pos(0)[1];
    reg14=reg15-reg14; reg12=(*f.m).alpha*reg12; reg11=reg13+reg11; reg13=elem.pos(1)[0]-elem.pos(0)[0]; reg15=elem.pos(1)[1]-elem.pos(0)[1];
    T reg16=elem.pos(2)[0]-elem.pos(0)[0]; T reg17=reg15*reg16; reg9=reg9/reg14; reg14=reg6/reg14; reg12=reg11+reg12;
    reg6=reg13*reg5; reg11=reg3*reg0; T reg18=reg4*reg0; reg17=reg6-reg17; reg9=reg9*reg12;
    reg12=reg14*reg12; reg15=reg15/reg17; reg6=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; reg14=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg13=reg13/reg17;
    reg5=reg5/reg17; T reg19=pow(reg3,2); T reg20=reg4*reg18; T reg21=1-(*f.m).resolution; reg12=reg9-reg12;
    reg9=pow(reg4,2); T reg22=reg3*reg11; reg16=reg16/reg17; T reg23=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; T reg24=vectors[0][indices[2]+0]-vectors[0][indices[0]+0];
    reg12=reg21*reg12; reg19=reg9-reg19; reg22=reg20-reg22; reg9=(*f.m).alpha*(*f.m).resolution; reg20=reg6*reg15;
    T reg25=reg24*reg13; T reg26=reg23*reg16; T reg27=reg14*reg5; reg11=reg11/reg22; reg20=reg27-reg20;
    reg26=reg25-reg26; reg12=reg9+reg12; reg19=reg19/reg22; reg22=reg18/reg22; reg14=reg14*reg16;
    reg6=reg6*reg13; reg23=reg23*reg5; reg24=reg24*reg15; reg12=(*f.m).deltaT*reg12; reg8=reg8*reg21;
    reg9=(*f.m).resolution*reg22; reg18=(*f.m).resolution*reg11; reg25=(*f.m).resolution*reg19; reg24=reg23-reg24; reg2=reg2*reg21;
    reg14=reg6-reg14; reg7=reg21*reg7; reg20=reg26+reg20; reg24=reg24-reg12; reg7=reg25+reg7;
    reg14=reg14-reg12; reg20=0.5*reg20; reg9=reg8+reg9; reg2=reg18+reg2; reg6=reg14*reg9;
    reg8=reg24*reg2; reg14=reg14*reg2; reg24=reg24*reg9; reg18=reg16-reg13; reg21=reg15-reg5;
    reg20=reg20*reg7; reg23=0.5*reg5; reg25=0.5*reg18; reg26=0.5*reg16; reg27=0.5*reg15;
    reg20=2*reg20; T reg28=0.5*reg21; T reg29=1-var_inter[0]; reg14=reg24+reg14; reg6=reg8+reg6;
    reg8=0.5*reg13; reg24=reg6*reg16; T reg30=reg23*reg20; T reg31=reg14*reg5; T reg32=reg26*reg20;
    T reg33=reg18*reg6; T reg34=reg28*reg20; T reg35=reg21*reg14; reg29=reg29-var_inter[1]; T reg36=reg6*reg13;
    T reg37=reg20*reg27; T reg38=reg20*reg25; T reg39=reg14*reg15; T reg40=reg8*reg20; T reg41=var_inter[1]*(*f.m).f_vol[1];
    reg36=reg36-reg37; reg38=reg35+reg38; reg35=reg29*(*f.m).f_vol[0]; reg33=reg34+reg33; reg34=reg29*(*f.m).f_vol[1];
    T reg42=var_inter[1]*(*f.m).f_vol[0]; reg30=reg30-reg24; T reg43=var_inter[0]*(*f.m).f_vol[1]; reg40=reg40-reg39; reg31=reg31-reg32;
    T reg44=var_inter[0]*(*f.m).f_vol[0]; reg33=reg33-reg34; reg31=reg31-reg44; reg36=reg36-reg41; reg30=reg30-reg43;
    reg40=reg40-reg42; reg38=reg38-reg35; reg36=reg36*reg17; reg38=reg38*reg17; reg30=reg30*reg17;
    reg31=reg31*reg17; reg33=reg33*reg17; reg40=reg40*reg17; sollicitation[indices[0]+0]+=ponderation*reg38; sollicitation[indices[1]+0]+=ponderation*reg31;
    sollicitation[indices[1]+1]+=ponderation*reg30; sollicitation[indices[2]+0]+=ponderation*reg40; sollicitation[indices[0]+1]+=ponderation*reg33; sollicitation[indices[2]+1]+=ponderation*reg36;
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
    vecs[0][indice+1]=node.dep[1]; vecs[1][indice+1]=node.dep[1]; vecs[2][indice+1]=node.dep[1]; vecs[3][indice+1]=node.dep[1]; vecs[4][indice+1]=node.dep[1];
    vecs[0][indice+0]=node.dep[0]; vecs[1][indice+0]=node.dep[0]; vecs[2][indice+0]=node.dep[0]; vecs[3][indice+0]=node.dep[0]; vecs[4][indice+0]=node.dep[0];
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg2*reg4; reg4=reg3*reg4; T reg6=reg2*reg1; reg1=reg3*reg1;
    T reg7=reg3*reg4; T reg8=reg2*reg5; reg4=reg2*reg4; T reg9=reg3*reg1; T reg10=reg2*reg6;
    reg1=reg2*reg1; reg1=reg10+reg1; reg5=reg3*reg5; reg6=reg3*reg6; reg4=reg8+reg4;
    reg7=reg7-reg8; reg9=reg9-reg10; reg9=reg3*reg9; reg5=reg8+reg5; reg8=0.5*elem.pos(1)[0];
    reg1=reg2*reg1; T reg11=reg3*reg7; T reg12=0.5*elem.pos(1)[1]; T reg13=0.5*elem.pos(0)[0]; T reg14=0.5*elem.pos(0)[1];
    T reg15=reg2*reg4; T reg16=reg10+reg6; T reg17=0.5*elem.pos(2)[0]; T reg18=reg12-reg14; reg1=reg9-reg1;
    reg16=reg2*reg16; reg9=reg8-reg13; reg15=reg11-reg15; reg11=reg2*reg5; reg12=reg14+reg12;
    reg14=0.5*elem.pos(2)[1]; reg8=reg13+reg8; reg18=reg14+reg18; reg8=reg17-reg8; reg13=0.5*elem.pos(3)[0];
    reg17=reg9+reg17; reg12=reg14-reg12; reg9=0.5*elem.pos(3)[1]; reg16=reg1-reg16; reg11=reg15-reg11;
    reg1=0.5*vectors[0][indices[0]+1]; reg14=0.5*vectors[0][indices[1]+1]; reg4=reg4/reg11; reg12=reg9+reg12; reg8=reg13+reg8;
    reg15=0.5*vectors[0][indices[1]+0]; reg9=reg18-reg9; reg18=0.5*vectors[0][indices[0]+0]; reg13=reg17-reg13; reg16=reg16/reg11;
    reg7=reg7/reg11; reg17=reg14+reg1; T reg19=reg15-reg18; T reg20=reg12*reg13; reg1=reg14-reg1;
    reg18=reg15+reg18; reg14=reg7*reg16; reg15=0.5*vectors[0][indices[2]+1]; T reg21=reg4*reg16; T reg22=0.5*vectors[0][indices[2]+0];
    T reg23=reg9*reg8; T reg24=0.21132486540518713447*elem.pos(1)[1]; T reg25=0.78867513459481286553*elem.pos(1)[1]; T reg26=0.21132486540518713447*elem.pos(0)[1]; T reg27=0.78867513459481286553*elem.pos(1)[0];
    T reg28=0.78867513459481286553*elem.pos(0)[1]; T reg29=0.21132486540518713447*elem.pos(0)[0]; T reg30=0.21132486540518713447*elem.pos(1)[0]; T reg31=0.78867513459481286553*elem.pos(0)[0]; reg19=reg19+reg22;
    T reg32=0.5*vectors[0][indices[3]+0]; reg11=reg5/reg11; reg5=0.21132486540518713447*elem.pos(2)[0]; T reg33=reg24+reg28; reg23=reg20-reg23;
    reg20=0.5*vectors[0][indices[3]+1]; T reg34=0.21132486540518713447*elem.pos(2)[1]; T reg35=reg29+reg27; reg17=reg15-reg17; reg24=reg24-reg26;
    reg18=reg22-reg18; reg22=reg30+reg31; reg26=reg26+reg25; T reg36=0.78867513459481286553*elem.pos(2)[1]; T reg37=reg7*reg14;
    reg1=reg15+reg1; reg15=reg4*reg21; reg29=reg30-reg29; reg30=reg4*(*f.m).alpha; T reg38=reg7*(*f.m).alpha;
    T reg39=0.78867513459481286553*elem.pos(2)[0]; reg8=reg8/reg23; reg17=reg17+reg20; reg15=reg37-reg15; reg30=reg38+reg30;
    reg11=reg11*(*f.m).alpha; reg22=reg5-reg22; reg31=reg27-reg31; reg28=reg25-reg28; reg20=reg1-reg20;
    reg33=reg34-reg33; reg13=reg13/reg23; reg29=reg29+reg39; reg1=0.21132486540518713447*elem.pos(3)[0]; reg35=reg39-reg35;
    reg25=0.78867513459481286553*elem.pos(3)[0]; reg27=0.21132486540518713447*elem.pos(3)[1]; reg37=0.78867513459481286553*elem.pos(3)[1]; reg24=reg36+reg24; reg26=reg36-reg26;
    reg36=reg2*reg0; reg38=reg3*reg0; reg12=reg12/reg23; reg19=reg19-reg32; reg23=reg9/reg23;
    reg18=reg32+reg18; reg9=reg2*reg36; reg14=reg14/reg15; reg35=reg35+reg1; reg32=0.78867513459481286553*PNODE(0).dep[0];
    reg24=reg24-reg37; reg39=0.78867513459481286553*PNODE(0).dep[1]; T reg40=0.78867513459481286553*PNODE(1).dep[0]; T reg41=0.78867513459481286553*PNODE(1).dep[1]; reg26=reg27+reg26;
    reg31=reg5+reg31; reg5=0.21132486540518713447*PNODE(0).dep[1]; reg28=reg34+reg28; reg29=reg29-reg25; reg33=reg37+reg33;
    reg34=0.21132486540518713447*PNODE(0).dep[0]; reg37=reg23*reg18; T reg42=reg12*reg19; T reg43=reg17*reg13; T reg44=0.21132486540518713447*PNODE(1).dep[0];
    T reg45=reg3*reg38; T reg46=reg20*reg8; reg15=reg21/reg15; reg21=0.21132486540518713447*PNODE(1).dep[1]; reg11=reg30+reg11;
    reg22=reg25+reg22; reg25=reg29*reg33; reg30=(*f.m).deltaT*(*f.m).alpha; T reg47=reg41+reg5; T reg48=reg24*reg22;
    T reg49=reg21+reg39; T reg50=0.21132486540518713447*PNODE(2).dep[0]; T reg51=0.21132486540518713447*PNODE(2).dep[1]; reg15=reg15*reg11; reg11=reg14*reg11;
    reg9=reg45-reg9; reg14=reg44+reg32; reg45=0.78867513459481286553*PNODE(2).dep[1]; reg44=reg44-reg34; reg37=reg42-reg37;
    elem.epsilon[0][0]=reg37; reg42=0.78867513459481286553*PNODE(2).dep[0]; reg46=reg43-reg46; elem.epsilon[0][1]=reg46; reg5=reg21-reg5;
    reg21=reg24*reg35; reg34=reg34+reg40; reg43=reg29*reg26; reg1=reg31-reg1; reg27=reg28-reg27;
    reg15=reg11-reg15; reg11=0.21132486540518713447*PNODE(3).dep[1]; reg44=reg44+reg42; reg47=reg45-reg47; reg28=1-(*f.m).resolution;
    reg31=0.78867513459481286553*PNODE(3).dep[1]; reg5=reg45+reg5; reg48=reg25-reg48; reg25=reg26*reg1; reg38=reg38/reg9;
    reg21=reg43-reg21; reg32=reg40-reg32; reg49=reg51-reg49; reg40=0.78867513459481286553*PNODE(3).dep[0]; reg43=reg37-reg30;
    reg34=reg42-reg34; reg42=reg46-reg30; reg14=reg50-reg14; reg45=0.21132486540518713447*PNODE(3).dep[0]; T reg52=reg35*reg27;
    reg36=reg36/reg9; reg39=reg41-reg39; reg34=reg45+reg34; reg47=reg47+reg11; reg41=reg26/reg21;
    T reg53=reg35/reg21; reg5=reg5-reg31; T reg54=reg24/reg21; reg15=reg28*reg15; T reg55=reg29/reg21;
    reg44=reg44-reg40; reg49=reg31+reg49; reg52=reg25-reg52; reg25=reg38*reg42; reg31=reg43*reg36;
    reg42=reg42*reg36; reg43=reg38*reg43; reg39=reg51+reg39; reg32=reg50+reg32; reg14=reg40+reg14;
    reg40=pow(reg3,2); reg50=pow(reg2,2); reg51=(*f.m).alpha*(*f.m).resolution; reg29=reg29/reg48; T reg56=reg22*reg27;
    T reg57=reg22/reg48; T reg58=reg33*reg1; T reg59=reg33/reg48; reg24=reg24/reg48; T reg60=reg5*reg53;
    reg26=reg26/reg52; reg35=reg35/reg52; reg11=reg39-reg11; reg45=reg32-reg45; reg32=reg41*reg44;
    reg39=reg27/reg52; reg42=reg43+reg42; reg56=reg58-reg56; reg25=reg31+reg25; reg31=reg1/reg52;
    reg43=reg47*reg55; reg58=reg24*reg14; T reg61=reg5*reg59; reg50=reg40-reg50; reg24=reg49*reg24;
    reg40=reg5*reg57; T reg62=reg49*reg29; reg53=reg53*reg44; reg55=reg55*reg34; reg15=reg51+reg15;
    reg51=reg47*reg54; reg41=reg5*reg41; reg29=reg29*reg14; reg59=reg44*reg59; reg54=reg54*reg34;
    reg57=reg44*reg57; reg5=reg34*reg31; reg44=reg47*reg39; T reg63=reg35*reg45; reg40=reg62-reg40;
    reg1=reg1/reg56; reg22=reg22/reg56; reg33=reg33/reg56; reg62=reg11*reg26; reg39=reg34*reg39;
    reg26=reg26*reg45; reg58=reg59-reg58; reg35=reg11*reg35; reg24=reg61-reg24; reg53=reg55-reg53;
    reg51=reg41-reg51; reg54=reg32-reg54; reg31=reg47*reg31; reg57=reg29-reg57; reg15=(*f.m).deltaT*reg15;
    reg29=reg38*(*f.m).resolution; reg7=reg7*reg28; reg27=reg27/reg56; reg32=reg2*reg42; reg34=reg3*reg25;
    reg4=reg4*reg28; reg41=(*f.m).resolution*reg36; reg25=reg2*reg25; reg42=reg3*reg42; reg9=reg50/reg9;
    reg60=reg43-reg60; reg42=reg42-reg25; reg34=reg34-reg32; reg35=reg31-reg35; reg31=reg11*reg33;
    reg53=reg51+reg53; reg4=reg41+reg4; reg29=reg7+reg29; reg7=reg27*reg49; reg41=reg14*reg1;
    reg28=reg16*reg28; reg16=reg45*reg22; reg40=reg40-reg15; reg43=(*f.m).resolution*reg9; reg57=reg24+reg57;
    reg58=reg58-reg15; reg60=reg60-reg15; reg54=reg54-reg15; reg33=reg45*reg33; reg22=reg11*reg22;
    reg1=reg49*reg1; reg14=reg27*reg14; reg39=reg26-reg39; reg44=reg62-reg44; reg63=reg5-reg63;
    reg5=reg29*reg60; reg11=reg4*reg54; reg32=reg25+reg32; reg7=reg31-reg7; reg63=reg44+reg63;
    reg16=reg41-reg16; reg57=0.5*reg57; reg22=reg1-reg22; reg1=reg4*reg40; reg24=reg29*reg58;
    reg25=reg4*reg58; reg26=reg29*reg40; reg20=reg12*reg20; reg35=reg35-reg15; reg14=reg33-reg14;
    reg42=reg30+reg42; reg8=reg19*reg8; reg39=reg39-reg15; reg53=0.5*reg53; reg12=reg4*reg60;
    reg34=reg30+reg34; reg17=reg23*reg17; reg19=reg29*reg54; reg13=reg18*reg13; reg28=reg43+reg28;
    reg26=reg25+reg26; reg1=reg24+reg1; reg22=reg22-reg15; reg17=reg20-reg17; reg18=reg28*reg57;
    reg32=reg30-reg32; reg8=reg13-reg8; reg63=0.5*reg63; reg13=reg4*reg39; reg20=reg4*reg35;
    reg23=reg42+reg34; reg24=reg29*reg35; reg25=reg29*reg39; reg5=reg11+reg5; reg16=reg7+reg16;
    reg7=reg28*reg53; reg12=reg19+reg12; reg14=reg14-reg15; reg11=reg4*reg22; reg19=reg14*reg29;
    reg8=reg17+reg8; reg7=2*reg7; reg17=reg29*reg22; reg27=reg14*reg4; reg18=2*reg18;
    reg16=0.5*reg16; reg23=reg23+reg32; reg31=reg28*reg63; reg12=reg54*reg12; reg5=reg60*reg5;
    reg24=reg13+reg24; reg1=reg58*reg1; reg20=reg25+reg20; reg26=reg40*reg26; reg31=2*reg31;
    reg12=reg5+reg12; reg17=reg27+reg17; reg24=reg35*reg24; reg7=reg53*reg7; reg8=0.5*reg8;
    elem.epsilon[0][2]=reg8; reg20=reg39*reg20; reg11=reg19+reg11; reg1=reg26+reg1; reg23=reg23/3;
    reg5=reg16*reg28; reg18=reg57*reg18; reg13=reg8*reg9; reg20=reg24+reg20; reg34=reg34-reg23;
    reg42=reg42-reg23; reg5=2*reg5; reg18=reg1+reg18; reg7=reg12+reg7; reg11=reg14*reg11;
    reg31=reg63*reg31; reg22=reg17*reg22; reg18=reg48*reg18; reg23=reg32-reg23; reg11=reg22+reg11;
    reg42=pow(reg42,2); reg7=reg21*reg7; reg34=pow(reg34,2); reg31=reg20+reg31; reg5=reg16*reg5;
    reg13=reg0*reg13; reg31=reg52*reg31; reg5=reg11+reg5; reg1=2*reg13; reg34=reg42+reg34;
    reg23=pow(reg23,2); reg7=0.25*reg7; reg18=0.25*reg18; reg18=reg7+reg18; reg31=0.25*reg31;
    reg56=reg5*reg56; reg1=reg13*reg1; reg23=reg34+reg23; reg1=reg23+reg1; reg37=reg37-reg15;
    reg56=0.25*reg56; reg46=reg46-reg15; reg31=reg18+reg31; reg31=reg56+reg31; reg5=reg37*reg29;
    reg7=reg46*reg29; reg1=1.5*reg1; reg46=reg46*reg4; reg37=reg37*reg4; elem.ener=reg31/2;
    elem.sigma[0][2]=reg8*reg28; elem.sigma[0][1]=reg37+reg7; elem.sigma_von_mises=pow(reg1,0.5); elem.sigma[0][0]=reg5+reg46;
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
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg3*reg2; reg2=reg4*reg2; T reg6=reg3*reg1; reg1=reg4*reg1;
    T reg7=reg3*reg5; T reg8=reg3*reg6; T reg9=reg4*reg1; T reg10=reg3*reg2; reg2=reg4*reg2;
    reg1=reg3*reg1; reg2=reg2-reg7; reg1=reg8+reg1; reg10=reg7+reg10; reg9=reg9-reg8;
    reg6=reg4*reg6; reg5=reg4*reg5; reg1=reg3*reg1; T reg11=reg8+reg6; reg5=reg7+reg5;
    reg9=reg4*reg9; reg7=reg3*reg10; T reg12=reg4*reg2; T reg13=reg3*reg5; reg7=reg12-reg7;
    reg11=reg3*reg11; reg1=reg9-reg1; reg9=1-var_inter[1]; reg11=reg1-reg11; reg1=1-var_inter[0];
    reg13=reg7-reg13; reg7=reg9*elem.pos(0)[1]; reg11=reg11/reg13; reg12=reg9*elem.pos(0)[0]; T reg14=reg9*elem.pos(1)[0];
    T reg15=reg1*elem.pos(0)[1]; T reg16=reg9*elem.pos(1)[1]; T reg17=reg1*elem.pos(0)[0]; reg2=reg2/reg13; T reg18=elem.pos(1)[0]*var_inter[0];
    reg10=reg10/reg13; T reg19=elem.pos(1)[1]*var_inter[0]; T reg20=reg15+reg19; T reg21=reg10*reg11; reg16=reg16-reg7;
    T reg22=elem.pos(2)[1]*var_inter[1]; T reg23=reg17+reg18; T reg24=elem.pos(2)[0]*var_inter[0]; T reg25=elem.pos(2)[1]*var_inter[0]; T reg26=elem.pos(2)[0]*var_inter[1];
    reg14=reg14-reg12; T reg27=reg2*reg11; T reg28=elem.pos(3)[1]*var_inter[1]; T reg29=elem.pos(3)[0]*var_inter[1]; reg26=reg14+reg26;
    reg22=reg16+reg22; reg13=reg5/reg13; reg5=reg10*reg21; reg14=reg10*(*f.m).alpha; reg16=reg2*(*f.m).alpha;
    reg24=reg24-reg23; T reg30=elem.pos(3)[1]*reg1; reg25=reg25-reg20; T reg31=elem.pos(3)[0]*reg1; T reg32=reg2*reg27;
    reg31=reg24+reg31; reg26=reg26-reg29; reg22=reg22-reg28; reg25=reg30+reg25; reg13=reg13*(*f.m).alpha;
    reg14=reg16+reg14; reg5=reg32-reg5; reg16=reg26*reg25; reg24=reg22*reg31; reg13=reg14+reg13;
    reg21=reg21/reg5; reg5=reg27/reg5; reg14=reg4*reg0; reg27=reg3*reg0; reg30=reg4*reg14;
    reg32=pow(reg4,2); T reg33=reg3*reg27; reg24=reg16-reg24; reg21=reg21*reg13; reg13=reg5*reg13;
    reg5=pow(reg3,2); reg21=reg13-reg21; reg26=reg26/reg24; reg33=reg30-reg33; reg25=reg25/reg24;
    reg31=reg31/reg24; reg22=reg22/reg24; reg13=1-(*f.m).resolution; reg5=reg32-reg5; reg16=reg9*reg31;
    reg5=reg5/reg33; reg30=var_inter[0]*reg22; reg32=var_inter[1]*reg31; T reg34=(*f.m).alpha*(*f.m).resolution; reg21=reg13*reg21;
    reg27=reg27/reg33; T reg35=var_inter[1]*reg25; T reg36=var_inter[0]*reg26; T reg37=reg1*reg26; T reg38=reg9*reg25;
    reg33=reg14/reg33; reg14=reg1*reg22; T reg39=(*f.m).resolution*reg5; T reg40=reg35+reg14; T reg41=reg32+reg37;
    reg21=reg34+reg21; reg34=(*f.m).resolution*reg27; reg10=reg10*reg13; reg11=reg11*reg13; T reg42=reg33*(*f.m).resolution;
    reg13=reg2*reg13; reg2=reg30+reg38; T reg43=reg36+reg16; T reg44=reg36-reg32; T reg45=reg35-reg30;
    reg11=reg39+reg11; reg10=reg34+reg10; reg21=(*f.m).deltaT*reg21; reg42=reg13+reg42; reg13=0.5*reg2;
    reg34=0.5*reg43; reg39=0.5*reg41; T reg46=reg14-reg38; T reg47=0.5*reg40; T reg48=reg16-reg37;
    T reg49=reg42*reg21; T reg50=reg34*reg11; T reg51=reg10*reg21; T reg52=0.5*reg44; T reg53=0.5*reg45;
    T reg54=0.5*reg48; T reg55=reg47*reg11; T reg56=0.5*reg46; T reg57=reg13*reg11; T reg58=reg39*reg11;
    T reg59=reg40*reg42; T reg60=reg1*var_inter[1]; T reg61=reg40*reg10; T reg62=2*reg55; reg58=2*reg58;
    T reg63=reg49+reg51; T reg64=reg11*reg56; T reg65=reg2*reg42; T reg66=reg41*reg10; T reg67=reg53*reg11;
    T reg68=reg11*reg54; T reg69=2*reg50; T reg70=reg41*reg42; T reg71=reg2*reg10; T reg72=reg43*reg42;
    T reg73=reg43*reg10; T reg74=reg52*reg11; T reg75=reg9*var_inter[0]; reg57=2*reg57; T reg76=reg2*reg59;
    T reg77=reg13*reg69; T reg78=reg40*reg66; T reg79=reg39*reg69; T reg80=reg40*reg65; T reg81=reg41*reg61;
    T reg82=reg42*reg48; T reg83=reg43*reg71; T reg84=reg42*reg46; T reg85=var_inter[0]*var_inter[1]; T reg86=reg9*reg1;
    T reg87=reg34*reg58; T reg88=reg45*reg42; T reg89=reg75*(*f.m).f_vol[1]; reg68=2*reg68; T reg90=reg10*reg48;
    T reg91=reg45*reg10; T reg92=reg60*(*f.m).f_vol[0]; T reg93=reg63*reg43; T reg94=reg39*reg62; T reg95=reg47*reg57;
    T reg96=reg44*reg10; reg67=2*reg67; T reg97=reg10*reg46; reg64=2*reg64; T reg98=reg13*reg62;
    T reg99=reg44*reg42; T reg100=reg63*reg40; T reg101=reg2*reg73; T reg102=reg41*reg72; T reg103=reg34*reg57;
    T reg104=reg70*reg43; reg74=2*reg74; T reg105=reg47*reg58; T reg106=reg47*reg74; T reg107=reg41*reg91;
    reg95=reg102+reg95; T reg108=reg86*(*f.m).f_vol[1]; T reg109=reg75*(*f.m).f_vol[0]; T reg110=reg85*(*f.m).f_vol[0]; T reg111=reg85*(*f.m).f_vol[1];
    T reg112=reg60*(*f.m).f_vol[1]; T reg113=reg84*reg46; T reg114=reg68*reg54; T reg115=reg64*reg56; T reg116=reg90*reg46;
    T reg117=reg64*reg54; T reg118=reg82*reg48; T reg119=reg65*reg46; T reg120=reg69*reg54; T reg121=reg73*reg46;
    T reg122=reg57*reg54; T reg123=reg68*reg56; T reg124=reg97*reg48; T reg125=reg88*reg46; T reg126=reg74*reg54;
    T reg127=reg34*reg67; T reg128=reg96*reg46; T reg129=reg67*reg54; T reg130=reg62*reg54; T reg131=reg59*reg46;
    T reg132=reg58*reg54; T reg133=reg66*reg46; T reg134=reg58*reg56; T reg135=reg61*reg48; T reg136=reg41*reg99;
    T reg137=reg47*reg67; reg105=reg81+reg105; T reg138=reg67*reg56; T reg139=reg70*reg41; T reg140=reg47*reg62;
    T reg141=reg99*reg48; T reg142=reg63*reg46; T reg143=reg63*reg48; T reg144=reg74*reg56; T reg145=reg91*reg48;
    T reg146=reg63*reg2; T reg147=reg93-reg89; T reg148=reg57*reg56; T reg149=reg63*reg45; T reg150=reg72*reg48;
    T reg151=reg63*reg44; T reg152=reg100-reg92; T reg153=reg69*reg56; T reg154=reg71*reg48; T reg155=reg63*reg41;
    reg78=reg94+reg78; T reg156=reg40*reg96; T reg157=reg39*reg67; T reg158=reg41*reg71; reg71=reg44*reg71;
    T reg159=reg53*reg69; T reg160=reg34*reg62; T reg161=reg13*reg57; T reg162=reg52*reg69; T reg163=reg2*reg66;
    T reg164=reg44*reg72; T reg165=reg53*reg57; T reg166=reg45*reg65; reg87=reg76+reg87; T reg167=reg44*reg91;
    T reg168=reg43*reg72; T reg169=reg52*reg64; T reg170=reg53*reg74; T reg171=reg45*reg90; T reg172=reg40*reg59;
    T reg173=reg44*reg99; T reg174=reg2*reg96; T reg175=reg53*reg67; T reg176=reg44*reg61; T reg177=reg45*reg88;
    T reg178=reg52*reg74; reg83=reg83+reg77; reg96=reg45*reg96; T reg179=reg52*reg67; T reg180=reg86*(*f.m).f_vol[0];
    T reg181=reg45*reg59; T reg182=reg52*reg58; T reg183=reg43*reg82; reg66=reg45*reg66; T reg184=reg52*reg62;
    T reg185=reg13*reg64; T reg186=reg44*reg97; T reg187=reg53*reg68; T reg188=reg52*reg57; T reg189=reg43*reg97;
    T reg190=reg44*reg82; T reg191=reg53*reg64; T reg192=reg13*reg68; T reg193=reg45*reg73; reg80=reg79+reg80;
    T reg194=reg43*reg61; T reg195=reg40*reg88; reg57=reg39*reg57; T reg196=reg40*reg73; T reg197=reg13*reg58;
    T reg198=reg39*reg58; T reg199=reg34*reg68; T reg200=reg2*reg84; T reg201=reg13*reg74; reg99=reg43*reg99;
    reg97=reg41*reg97; T reg202=reg47*reg68; reg67=reg13*reg67; reg82=reg41*reg82; T reg203=reg62*reg56;
    T reg204=reg47*reg64; T reg205=reg70*reg48; T reg206=reg47*reg69; reg91=reg43*reg91; T reg207=reg34*reg74;
    reg58=reg53*reg58; reg88=reg2*reg88; T reg208=reg52*reg68; reg103=reg101+reg103; T reg209=reg45*reg84;
    reg70=reg70*reg44; T reg210=reg53*reg62; T reg211=reg34*reg69; reg68=reg39*reg68; reg65=reg2*reg65;
    reg84=reg40*reg84; reg104=reg98+reg104; T reg212=reg39*reg64; reg64=reg34*reg64; T reg213=reg2*reg90;
    reg90=reg40*reg90; reg74=reg39*reg74; reg115=reg118+reg115; reg118=reg83*reg24; reg114=reg113+reg114;
    reg161=reg161+reg168; reg205=reg205-reg203; reg132=reg132-reg131; reg195=reg74-reg195; reg199=reg200-reg199;
    reg133=reg133-reg130; reg129=reg128+reg129; reg64=reg213-reg64; reg65=reg65+reg211; reg126=reg125+reg126;
    reg74=reg103*reg24; reg207=reg88-reg207; reg122=reg122-reg121; reg198=reg172+reg198; reg127=reg174-reg127;
    reg88=reg87*reg24; reg123=reg124+reg123; reg119=reg119-reg120; reg163=reg163+reg160; reg189=reg192-reg189;
    reg117=reg116+reg117; reg183=reg185-reg183; reg147=reg147*reg24; reg170=reg167+reg170; reg113=reg146+reg109;
    reg175=reg173+reg175; reg116=reg143+reg108; reg144=reg145+reg144; reg124=reg142+reg180; reg58=reg58-reg176;
    reg70=reg70-reg210; reg156=reg157-reg156; reg139=reg139+reg140; reg84=reg68-reg84; reg138=reg141+reg138;
    reg68=reg105*reg24; reg90=reg212-reg90; reg125=reg80*reg24; reg137=reg136-reg137; reg57=reg57+reg196;
    reg128=reg78*reg24; reg134=reg134-reg135; reg106=reg107-reg106; reg202=reg97-reg202; reg97=reg95*reg24;
    reg204=reg82-reg204; reg158=reg206+reg158; reg91=reg201-reg91; reg99=reg67-reg99; reg197=reg197+reg194;
    reg67=reg104*reg24; reg208=reg209+reg208; reg169=reg171+reg169; reg166=reg166-reg162; reg188=reg188-reg193;
    reg178=reg177+reg178; reg179=reg96+reg179; reg182=reg182-reg181; reg66=reg66-reg184; reg187=reg186+reg187;
    reg148=reg148-reg150; reg165=reg165-reg164; reg82=reg149+reg110; reg71=reg71-reg159; reg96=reg151+reg111;
    reg154=reg154-reg153; reg152=reg152*reg24; reg191=reg190+reg191; reg107=reg155+reg112; reg138=reg138*reg24;
    reg154=reg154*reg24; reg144=reg144*reg24; reg123=reg123*reg24; reg115=reg115*reg24; reg134=reg134*reg24;
    reg195=reg24*reg195; reg133=reg133*reg24; reg148=reg148*reg24; reg136=ponderation*reg67; reg57=reg57*reg24;
    reg141=ponderation*reg128; reg197=reg197*reg24; reg202=reg202*reg24; reg204=reg204*reg24; reg99=reg99*reg24;
    reg158=reg158*reg24; reg145=ponderation*reg97; reg91=reg91*reg24; reg198=reg198*reg24; reg161=reg161*reg24;
    reg106=reg106*reg24; reg137=reg137*reg24; reg157=ponderation*reg118; reg167=ponderation*reg68; reg183=reg183*reg24;
    reg139=reg139*reg24; reg182=reg182*reg24; reg179=reg179*reg24; reg66=reg66*reg24; reg187=reg187*reg24;
    reg178=reg178*reg24; reg191=reg191*reg24; reg71=reg71*reg24; reg188=reg188*reg24; reg165=reg165*reg24;
    reg170=reg170*reg24; reg166=reg166*reg24; reg175=reg175*reg24; reg58=reg58*reg24; reg169=reg169*reg24;
    reg70=reg70*reg24; reg84=reg84*reg24; reg208=reg208*reg24; reg90=reg90*reg24; reg171=ponderation*reg125;
    reg119=reg119*reg24; reg173=ponderation*reg74; reg174=reg107*reg24; reg156=reg156*reg24; reg122=reg122*reg24;
    reg152=ponderation*reg152; reg207=reg207*reg24; reg177=reg96*reg24; reg114=reg114*reg24; reg126=reg126*reg24;
    reg127=reg127*reg24; reg185=reg82*reg24; reg64=reg64*reg24; reg199=reg199*reg24; reg147=ponderation*reg147;
    reg65=reg65*reg24; reg186=ponderation*reg88; reg129=reg129*reg24; reg190=reg113*reg24; reg117=reg117*reg24;
    reg163=reg163*reg24; reg192=reg116*reg24; reg205=reg205*reg24; reg132=reg132*reg24; reg200=reg124*reg24;
    reg189=reg189*reg24; T tmp_4_1=ponderation*reg169; T tmp_6_3=ponderation*reg57; T tmp_4_2=ponderation*reg166; T tmp_2_0=ponderation*reg199;
    T tmp_4_3=ponderation*reg188; T tmp_1_7=ponderation*reg205; T tmp_4_4=ponderation*reg178; T tmp_6_4=ponderation*reg195; T tmp_2_1=ponderation*reg64;
    T tmp_4_0=ponderation*reg208; T tmp_2_2=ponderation*reg65; T tmp_3_7=-reg136; T tmp_2_3=-reg173; T tmp_3_6=ponderation*reg197;
    T tmp_3_5=ponderation*reg99; T tmp_2_4=ponderation*reg207; T tmp_3_4=ponderation*reg91; T tmp_2_5=ponderation*reg127; T tmp_3_3=ponderation*reg161;
    T tmp_2_6=-reg186; T tmp_3_2=-reg157; T tmp_3_1=ponderation*reg183; T tmp_2_7=ponderation*reg163; T tmp_3_0=ponderation*reg189;
    T tmp_7_5=ponderation*reg137; T tmp_7_6=-reg167; T tmp_7_7=ponderation*reg139; T tmp_1_4=ponderation*reg144; reg57=ponderation*reg200;
    sollicitation[indices[0]+0]+=reg57; reg64=ponderation*reg192; sollicitation[indices[0]+1]+=reg64; T tmp_1_3=ponderation*reg148; reg65=ponderation*reg190;
    sollicitation[indices[1]+0]+=reg65; sollicitation[indices[1]+1]+=-reg147; reg91=ponderation*reg185; sollicitation[indices[2]+0]+=reg91; T tmp_1_2=ponderation*reg154;
    reg99=ponderation*reg177; sollicitation[indices[2]+1]+=reg99; sollicitation[indices[3]+0]+=-reg152; reg127=ponderation*reg174; sollicitation[indices[3]+1]+=reg127;
    T tmp_1_1=ponderation*reg115; T tmp_0_0=ponderation*reg114; T tmp_0_1=ponderation*reg117; T tmp_1_0=ponderation*reg123; T tmp_0_2=ponderation*reg119;
    T tmp_0_3=ponderation*reg122; T tmp_0_7=ponderation*reg133; T tmp_0_4=ponderation*reg126; T tmp_0_5=ponderation*reg129; T tmp_0_6=ponderation*reg132;
    T tmp_4_5=ponderation*reg179; T tmp_4_6=ponderation*reg182; T tmp_4_7=ponderation*reg66; T tmp_5_0=ponderation*reg187; T tmp_5_1=ponderation*reg191;
    T tmp_5_2=ponderation*reg71; T tmp_5_3=ponderation*reg165; T tmp_5_4=ponderation*reg170; T tmp_5_5=ponderation*reg175; T tmp_5_6=ponderation*reg58;
    T tmp_5_7=ponderation*reg70; T tmp_6_0=ponderation*reg84; T tmp_6_1=ponderation*reg90; T tmp_6_2=-reg171; T tmp_6_5=ponderation*reg156;
    T tmp_6_7=-reg141; T tmp_7_0=ponderation*reg202; T tmp_7_1=ponderation*reg204; T tmp_7_2=ponderation*reg158; T tmp_7_3=-reg145;
    T tmp_1_6=ponderation*reg134; T tmp_6_6=ponderation*reg198; T tmp_7_4=ponderation*reg106; T tmp_1_5=ponderation*reg138;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg2*reg1; reg1=reg3*reg1; T reg6=reg2*reg4; reg4=reg3*reg4;
    T reg7=reg3*reg1; T reg8=reg2*reg5; reg5=reg3*reg5; T reg9=reg3*reg6; T reg10=reg3*reg4;
    reg6=reg2*reg6; reg4=reg2*reg4; reg8=reg8-reg7; reg9=reg10+reg9; reg6=reg6-reg10;
    reg5=reg7+reg5; reg1=reg2*reg1; T reg11=reg2*reg6; T reg12=reg3*reg9; reg5=reg3*reg5;
    reg4=reg10+reg4; reg10=reg7+reg1; reg8=reg2*reg8; reg10=reg3*reg10; T reg13=reg3*reg4;
    reg12=reg11-reg12; reg5=reg8-reg5; reg13=reg12-reg13; reg8=1-var_inter[1]; reg10=reg5-reg10;
    reg5=1-var_inter[0]; reg10=reg10/reg13; reg11=reg8*elem.pos(0)[0]; reg12=reg8*elem.pos(1)[0]; T reg14=reg8*elem.pos(0)[1];
    T reg15=reg8*elem.pos(1)[1]; reg6=reg6/reg13; T reg16=reg5*elem.pos(0)[0]; T reg17=elem.pos(1)[0]*var_inter[0]; reg9=reg9/reg13;
    T reg18=reg5*elem.pos(0)[1]; T reg19=elem.pos(1)[1]*var_inter[0]; T reg20=reg6*reg10; reg12=reg12-reg11; T reg21=reg9*reg10;
    T reg22=elem.pos(2)[1]*var_inter[0]; T reg23=reg18+reg19; T reg24=elem.pos(2)[0]*var_inter[1]; T reg25=elem.pos(2)[0]*var_inter[0]; T reg26=reg16+reg17;
    T reg27=elem.pos(2)[1]*var_inter[1]; reg15=reg15-reg14; reg24=reg12+reg24; reg12=elem.pos(3)[0]*var_inter[1]; T reg28=elem.pos(3)[1]*var_inter[1];
    reg13=reg4/reg13; reg27=reg15+reg27; reg25=reg25-reg26; reg4=elem.pos(3)[0]*reg5; reg15=elem.pos(3)[1]*reg5;
    reg22=reg22-reg23; T reg29=reg6*reg20; T reg30=reg9*reg21; T reg31=reg6*(*f.m).alpha; T reg32=reg9*(*f.m).alpha;
    reg22=reg15+reg22; reg4=reg25+reg4; reg27=reg27-reg28; reg30=reg29-reg30; reg24=reg24-reg12;
    reg32=reg31+reg32; reg13=reg13*(*f.m).alpha; reg15=reg27*reg4; reg13=reg32+reg13; reg25=reg24*reg22;
    reg20=reg20/reg30; reg29=reg2*reg0; reg31=reg3*reg0; reg30=reg21/reg30; reg30=reg30*reg13;
    reg13=reg20*reg13; reg20=pow(reg2,2); reg15=reg25-reg15; reg21=pow(reg3,2); reg25=reg3*reg31;
    reg32=reg2*reg29; reg27=reg27/reg15; reg4=reg4/reg15; reg24=reg24/reg15; reg22=reg22/reg15;
    reg21=reg20-reg21; reg25=reg32-reg25; reg20=1-(*f.m).resolution; reg30=reg13-reg30; reg31=reg31/reg25;
    reg21=reg21/reg25; reg25=reg29/reg25; reg13=reg5*reg24; reg29=reg5*reg27; reg32=reg8*reg22;
    T reg33=var_inter[1]*reg22; reg30=reg20*reg30; T reg34=var_inter[1]*reg4; T reg35=(*f.m).alpha*(*f.m).resolution; T reg36=var_inter[0]*reg27;
    T reg37=reg8*reg4; reg30=reg35+reg30; reg35=reg25*(*f.m).resolution; reg6=reg6*reg20; reg10=reg10*reg20;
    reg20=reg9*reg20; reg9=(*f.m).resolution*reg31; T reg38=reg36+reg32; T reg39=(*f.m).resolution*reg21; T reg40=var_inter[0]*reg24;
    T reg41=reg33+reg29; T reg42=reg34+reg13; T reg43=reg40-reg34; T reg44=reg33-reg36; reg35=reg6+reg35;
    reg20=reg9+reg20; reg10=reg39+reg10; reg6=0.5*reg41; reg9=0.5*reg42; reg30=(*f.m).deltaT*reg30;
    reg39=reg29-reg32; T reg45=0.5*reg38; T reg46=reg37-reg13; T reg47=reg40+reg37; T reg48=0.5*reg47;
    T reg49=reg9*reg10; T reg50=reg6*reg10; T reg51=0.5*reg43; T reg52=reg20*reg30; T reg53=reg35*reg30;
    T reg54=reg45*reg10; T reg55=0.5*reg39; T reg56=0.5*reg46; T reg57=0.5*reg44; T reg58=reg53+reg52;
    T reg59=reg51*reg10; T reg60=reg57*reg10; T reg61=reg41*reg35; reg54=2*reg54; reg49=2*reg49;
    T reg62=reg47*reg20; T reg63=reg42*reg20; T reg64=2*reg50; T reg65=reg8*var_inter[0]; T reg66=reg48*reg10;
    T reg67=reg10*reg55; T reg68=reg5*var_inter[1]; T reg69=reg10*reg56; T reg70=reg42*reg35; T reg71=reg20*reg46;
    reg60=2*reg60; T reg72=reg58*reg41; T reg73=reg43*reg20; reg59=2*reg59; reg67=2*reg67;
    reg69=2*reg69; T reg74=reg38*reg35; T reg75=reg41*reg63; T reg76=reg44*reg35; T reg77=reg35*reg39;
    T reg78=reg68*(*f.m).f_vol[0]; T reg79=2*reg66; T reg80=reg65*(*f.m).f_vol[1]; T reg81=var_inter[0]*var_inter[1]; T reg82=reg8*reg5;
    T reg83=reg38*reg62; T reg84=reg48*reg54; T reg85=reg38*reg61; T reg86=reg48*reg49; T reg87=reg70*reg47;
    T reg88=reg45*reg64; T reg89=reg9*reg64; T reg90=reg41*reg20; T reg91=reg43*reg35; T reg92=reg44*reg20;
    T reg93=reg58*reg47; T reg94=reg35*reg46; T reg95=reg38*reg20; T reg96=reg47*reg35; reg86=reg85+reg86;
    T reg97=reg54*reg55; T reg98=reg49*reg56; T reg99=reg67*reg55; T reg100=reg63*reg39; T reg101=reg82*(*f.m).f_vol[0];
    T reg102=reg41*reg61; T reg103=reg38*reg73; T reg104=reg48*reg59; T reg105=reg38*reg76; reg84=reg83+reg84;
    T reg106=reg94*reg46; T reg107=reg74*reg39; T reg108=reg96*reg46; T reg109=reg48*reg79; T reg110=reg38*reg74;
    T reg111=reg64*reg55; T reg112=reg67*reg56; T reg113=reg70*reg46; T reg114=reg64*reg56; T reg115=reg71*reg39;
    T reg116=reg73*reg39; T reg117=reg48*reg60; T reg118=reg60*reg56; T reg119=reg59*reg56; T reg120=reg92*reg46;
    T reg121=reg59*reg55; T reg122=reg76*reg39; T reg123=reg79*reg55; T reg124=reg91*reg46; T reg125=reg54*reg56;
    T reg126=reg60*reg55; T reg127=reg61*reg39; T reg128=reg62*reg39; T reg129=reg90*reg46; T reg130=reg95*reg46;
    T reg131=reg49*reg55; T reg132=reg47*reg92; T reg133=reg45*reg59; T reg134=reg47*reg96; T reg135=reg79*reg56;
    T reg136=reg45*reg54; T reg137=reg48*reg64; T reg138=reg38*reg63; T reg139=reg70*reg43; reg75=reg89+reg75;
    T reg140=reg57*reg64; reg87=reg88+reg87; T reg141=reg58*reg42; T reg142=reg47*reg90; T reg143=reg45*reg49;
    T reg144=reg72-reg78; T reg145=reg47*reg91; T reg146=reg9*reg49; T reg147=reg58*reg43; T reg148=reg57*reg49;
    reg70=reg70*reg42; T reg149=reg6*reg64; T reg150=reg43*reg91; T reg151=reg43*reg90; T reg152=reg58*reg44;
    T reg153=reg57*reg60; T reg154=reg93-reg80; T reg155=reg58*reg39; T reg156=reg58*reg46; T reg157=reg58*reg38;
    reg63=reg44*reg63; T reg158=reg51*reg49; T reg159=reg69*reg56; T reg160=reg44*reg61; T reg161=reg77*reg39;
    T reg162=reg68*(*f.m).f_vol[1]; T reg163=reg81*(*f.m).f_vol[1]; T reg164=reg81*(*f.m).f_vol[0]; T reg165=reg65*(*f.m).f_vol[0]; T reg166=reg51*reg60;
    T reg167=reg44*reg73; T reg168=reg82*(*f.m).f_vol[1]; T reg169=reg51*reg64; T reg170=reg45*reg60; T reg171=reg51*reg59;
    T reg172=reg44*reg76; reg63=reg63-reg169; reg130=reg130-reg123; reg110=reg110+reg109; T reg173=reg156+reg168;
    reg158=reg158-reg160; reg97=reg97-reg108; T reg174=reg155+reg101; reg121=reg120+reg121; reg153=reg150+reg153;
    reg120=reg87*reg15; reg70=reg70+reg149; reg126=reg124+reg126; reg117=reg103-reg117; reg166=reg167+reg166;
    reg103=reg86*reg15; reg124=reg84*reg15; reg143=reg143+reg142; reg171=reg172+reg171; reg131=reg131-reg129;
    reg150=reg75*reg15; reg136=reg136+reg134; reg104=reg105-reg104; reg145=reg170-reg145; reg146=reg102+reg146;
    reg138=reg138+reg137; reg144=reg144*reg15; reg118=reg116+reg118; reg139=reg139-reg140; reg119=reg122+reg119;
    reg105=reg141+reg162; reg116=reg147+reg163; reg122=reg152+reg164; reg98=reg98-reg127; reg113=reg113-reg111;
    reg125=reg125-reg128; reg154=reg154*reg15; reg100=reg100-reg114; reg132=reg133-reg132; reg148=reg148-reg151;
    reg107=reg107-reg135; reg133=reg157+reg165; reg99=reg106+reg99; reg112=reg115+reg112; reg159=reg161+reg159;
    reg145=reg145*reg15; reg144=ponderation*reg144; reg158=reg158*reg15; reg139=reg139*reg15; reg159=reg159*reg15;
    reg119=reg119*reg15; reg143=reg143*reg15; reg104=reg104*reg15; reg106=ponderation*reg103; reg115=reg105*reg15;
    reg112=reg112*reg15; reg161=ponderation*reg120; reg125=reg125*reg15; reg166=reg166*reg15; reg148=reg148*reg15;
    reg117=reg117*reg15; reg107=reg107*reg15; reg132=reg132*reg15; reg171=reg171*reg15; reg153=reg153*reg15;
    reg167=reg173*reg15; reg130=reg130*reg15; reg99=reg99*reg15; reg110=reg110*reg15; reg170=reg133*reg15;
    reg97=reg97*reg15; reg172=reg174*reg15; reg113=reg113*reg15; reg100=reg100*reg15; reg121=reg121*reg15;
    reg70=reg70*reg15; reg154=ponderation*reg154; reg126=reg126*reg15; reg138=reg138*reg15; reg63=reg63*reg15;
    T reg175=reg116*reg15; reg136=reg136*reg15; reg118=reg118*reg15; T reg176=ponderation*reg150; reg146=reg146*reg15;
    T reg177=ponderation*reg124; reg131=reg131*reg15; T reg178=reg122*reg15; reg98=reg98*reg15; T tmp_2_3=-reg177;
    T tmp_1_7=ponderation*reg113; T tmp_2_2=ponderation*reg110; T tmp_4_7=ponderation*reg63; T tmp_5_6=ponderation*reg148; T tmp_2_4=ponderation*reg104;
    T tmp_5_7=ponderation*reg139; T tmp_5_5=ponderation*reg153; T tmp_1_2=ponderation*reg130; reg63=ponderation*reg167; sollicitation[indices[0]+1]+=reg63;
    T tmp_1_1=ponderation*reg99; reg99=ponderation*reg170; sollicitation[indices[1]+0]+=reg99; T tmp_0_7=ponderation*reg100; sollicitation[indices[1]+1]+=-reg154;
    T tmp_0_6=ponderation*reg98; reg98=ponderation*reg178; sollicitation[indices[2]+0]+=reg98; T tmp_0_5=ponderation*reg118; reg100=ponderation*reg175;
    sollicitation[indices[2]+1]+=reg100; sollicitation[indices[3]+0]+=-reg144; T tmp_0_4=ponderation*reg119; reg104=ponderation*reg115; sollicitation[indices[3]+1]+=reg104;
    T tmp_0_3=ponderation*reg125; T tmp_3_4=ponderation*reg132; T tmp_0_2=ponderation*reg107; T tmp_0_1=ponderation*reg112; T tmp_0_0=ponderation*reg159;
    T tmp_4_6=ponderation*reg158; T tmp_4_5=ponderation*reg166; T tmp_4_4=ponderation*reg171; T tmp_3_7=-reg161; T tmp_2_5=ponderation*reg117;
    T tmp_3_6=ponderation*reg143; T tmp_2_6=-reg106; T tmp_3_5=ponderation*reg145; T tmp_2_7=ponderation*reg138; T tmp_3_3=ponderation*reg136;
    T tmp_6_7=-reg176; T tmp_1_6=ponderation*reg131; T tmp_6_6=ponderation*reg146; T tmp_1_5=ponderation*reg126; T tmp_7_7=ponderation*reg70;
    T tmp_1_4=ponderation*reg121; T tmp_1_3=ponderation*reg97; reg70=ponderation*reg172; sollicitation[indices[0]+0]+=reg70;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=1-var_inter[1]; T reg2=pow(reg0,2); T reg3=1-var_inter[0];
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg0*reg2; T reg6=1.0/(*f.m).elastic_modulus; T reg7=reg1*elem.pos(0)[0]; T reg8=reg1*elem.pos(1)[0];
    T reg9=reg1*elem.pos(0)[1]; T reg10=reg1*elem.pos(1)[1]; T reg11=reg3*elem.pos(0)[0]; T reg12=elem.pos(1)[0]*var_inter[0]; T reg13=reg3*elem.pos(0)[1];
    T reg14=elem.pos(1)[1]*var_inter[0]; reg8=reg8-reg7; T reg15=elem.pos(2)[0]*var_inter[1]; T reg16=reg6*reg5; reg5=reg4*reg5;
    reg10=reg10-reg9; T reg17=elem.pos(2)[1]*var_inter[1]; T reg18=reg4*reg2; T reg19=elem.pos(2)[1]*var_inter[0]; T reg20=reg11+reg12;
    T reg21=elem.pos(2)[0]*var_inter[0]; reg2=reg6*reg2; T reg22=reg13+reg14; T reg23=reg4*reg5; T reg24=reg4*reg16;
    reg16=reg6*reg16; T reg25=reg6*reg2; T reg26=reg4*reg18; reg2=reg4*reg2; T reg27=elem.pos(3)[1]*reg3;
    reg15=reg8+reg15; reg8=elem.pos(3)[0]*var_inter[1]; T reg28=elem.pos(3)[1]*var_inter[1]; T reg29=elem.pos(3)[0]*reg3; reg21=reg21-reg20;
    reg19=reg19-reg22; reg17=reg10+reg17; reg18=reg6*reg18; reg25=reg25-reg26; reg29=reg21+reg29;
    reg2=reg26+reg2; reg5=reg6*reg5; reg24=reg23+reg24; reg16=reg16-reg23; reg17=reg17-reg28;
    reg15=reg15-reg8; reg19=reg27+reg19; reg10=reg17*reg29; reg21=reg4*reg24; reg2=reg4*reg2;
    reg27=reg26+reg18; reg25=reg6*reg25; T reg30=reg4*reg0; T reg31=reg15*reg19; T reg32=reg6*reg0;
    reg5=reg23+reg5; reg23=reg6*reg16; T reg33=reg4*reg30; reg10=reg31-reg10; reg31=pow(reg4,2);
    reg27=reg4*reg27; T reg34=reg6*reg32; reg2=reg25-reg2; reg25=pow(reg6,2); reg21=reg23-reg21;
    reg23=reg4*reg5; reg29=reg29/reg10; reg23=reg21-reg23; reg17=reg17/reg10; reg19=reg19/reg10;
    reg27=reg2-reg27; reg15=reg15/reg10; reg31=reg25-reg31; reg33=reg34-reg33; reg2=reg1*reg29;
    reg21=reg3*reg15; reg25=reg3*reg17; reg34=reg1*reg19; T reg35=1-(*f.m).resolution; reg31=reg31/reg33;
    T reg36=var_inter[1]*reg29; T reg37=var_inter[1]*reg19; T reg38=var_inter[0]*reg15; T reg39=var_inter[0]*reg17; reg27=reg27/reg23;
    T reg40=reg37+reg25; T reg41=reg36+reg21; reg30=reg30/reg33; reg24=reg24/reg23; T reg42=reg38+reg2;
    T reg43=(*f.m).resolution*reg31; T reg44=reg27*reg35; reg16=reg16/reg23; T reg45=reg39+reg34; reg33=reg32/reg33;
    reg32=reg33*(*f.m).resolution; T reg46=reg16*reg35; reg44=reg43+reg44; reg43=reg38-reg36; T reg47=reg37-reg39;
    T reg48=reg25-reg34; T reg49=0.5*reg45; T reg50=0.5*reg41; T reg51=0.5*reg42; T reg52=reg2-reg21;
    T reg53=(*f.m).resolution*reg30; T reg54=reg24*reg35; T reg55=0.5*reg40; T reg56=reg50*reg44; T reg57=reg55*reg44;
    T reg58=0.5*reg52; T reg59=reg49*reg44; T reg60=0.5*reg43; T reg61=0.5*reg48; reg54=reg53+reg54;
    reg53=0.5*reg47; reg32=reg46+reg32; reg46=reg51*reg44; T reg62=reg60*reg44; T reg63=reg41*reg32;
    reg59=2*reg59; T reg64=reg41*reg54; T reg65=2*reg57; reg56=2*reg56; T reg66=reg42*reg32;
    T reg67=reg42*reg54; T reg68=2*reg46; T reg69=reg45*reg32; T reg70=reg40*reg32; T reg71=reg53*reg44;
    T reg72=reg44*reg58; T reg73=reg44*reg61; T reg74=reg40*reg54; T reg75=reg45*reg54; T reg76=reg51*reg56;
    reg72=2*reg72; T reg77=reg49*reg65; T reg78=reg63*reg42; T reg79=reg54*reg48; T reg80=reg40*reg69;
    T reg81=reg50*reg68; T reg82=reg32*reg52; T reg83=reg49*reg68; reg71=2*reg71; T reg84=reg43*reg54;
    T reg85=reg42*reg75; reg62=2*reg62; T reg86=reg47*reg32; T reg87=reg45*reg67; T reg88=reg45*reg70;
    T reg89=reg51*reg59; T reg90=reg41*reg66; T reg91=reg55*reg59; T reg92=reg50*reg65; T reg93=reg41*reg74;
    T reg94=reg55*reg56; reg73=2*reg73; T reg95=reg40*reg64; T reg96=reg54*reg52; T reg97=reg43*reg32;
    T reg98=reg47*reg54; T reg99=reg32*reg48; reg89=reg87+reg89; T reg100=reg67*reg48; reg91=reg90+reg91;
    T reg101=reg60*reg56; T reg102=reg43*reg82; T reg103=reg53*reg73; T reg104=reg47*reg70; T reg105=reg50*reg62;
    T reg106=reg55*reg68; T reg107=reg43*reg75; T reg108=reg59*reg58; T reg109=reg55*reg73; T reg110=reg41*reg82;
    T reg111=reg47*reg67; T reg112=reg53*reg68; T reg113=reg86*reg48; T reg114=reg55*reg72; T reg115=reg41*reg79;
    T reg116=reg62*reg61; T reg117=reg98*reg52; T reg118=reg62*reg58; T reg119=reg72*reg58; T reg120=reg60*reg59;
    T reg121=reg63*reg52; reg85=reg85+reg83; T reg122=reg43*reg66; T reg123=reg51*reg71; T reg124=reg51*reg68;
    T reg125=reg74*reg52; T reg126=reg41*reg75; T reg127=reg50*reg71; T reg128=reg40*reg84; T reg129=reg96*reg48;
    reg95=reg92+reg95; T reg130=reg47*reg64; T reg131=reg51*reg62; T reg132=reg60*reg65; T reg133=reg56*reg61;
    T reg134=reg55*reg65; T reg135=reg63*reg41; T reg136=reg45*reg86; T reg137=reg73*reg58; T reg138=reg42*reg82;
    reg94=reg93+reg94; T reg139=reg99*reg48; T reg140=reg71*reg61; T reg141=reg45*reg84; T reg142=reg43*reg79;
    T reg143=reg69*reg48; T reg144=reg55*reg71; T reg145=reg41*reg97; T reg146=reg47*reg69; T reg147=reg53*reg72;
    T reg148=reg97*reg52; T reg149=reg68*reg58; T reg150=reg49*reg73; T reg151=reg60*reg68; T reg152=reg55*reg62;
    T reg153=reg41*reg98; T reg154=reg60*reg73; T reg155=reg47*reg86; T reg156=reg43*reg97; T reg157=reg50*reg56;
    T reg158=reg65*reg58; T reg159=reg60*reg62; T reg160=reg45*reg64; T reg161=reg40*reg67; T reg162=reg50*reg59;
    T reg163=reg53*reg71; T reg164=reg79*reg52; T reg165=reg72*reg61; reg80=reg81+reg80; T reg166=reg47*reg84;
    reg86=reg40*reg86; T reg167=reg49*reg72; T reg168=reg51*reg65; T reg169=reg43*reg74; T reg170=reg40*reg96;
    T reg171=reg50*reg73; T reg172=reg51*reg72; reg82=reg82*reg52; T reg173=reg73*reg61; T reg174=reg40*reg99;
    T reg175=reg50*reg72; T reg176=reg45*reg99; T reg177=reg53*reg56; T reg178=reg68*reg61; T reg179=reg53*reg65;
    reg63=reg63*reg43; reg75=reg75*reg52; T reg180=reg49*reg59; T reg181=reg42*reg66; reg84=reg84*reg48;
    T reg182=reg53*reg59; reg69=reg45*reg69; T reg183=reg49*reg62; T reg184=reg42*reg98; T reg185=reg71*reg58;
    T reg186=reg49*reg71; reg97=reg42*reg97; T reg187=reg40*reg70; T reg188=reg65*reg61; T reg189=reg70*reg48;
    T reg190=reg49*reg56; T reg191=reg42*reg74; reg98=reg43*reg98; reg62=reg53*reg62; reg56=reg56*reg58;
    reg71=reg60*reg71; reg59=reg59*reg61; reg78=reg77+reg78; reg79=reg42*reg79; reg76=reg88+reg76;
    reg73=reg51*reg73; reg99=reg47*reg99; reg72=reg60*reg72; T reg192=reg66*reg52; reg64=reg64*reg48;
    T reg193=reg45*reg96; reg96=reg47*reg96; reg119=reg139+reg119; reg128=reg127-reg128; reg127=reg76*reg10;
    reg133=reg133-reg125; reg140=reg148+reg140; reg116=reg117+reg116; reg59=reg59-reg192; reg75=reg75-reg178;
    reg173=reg82+reg173; reg123=reg141-reg123; reg165=reg164+reg165; reg64=reg64-reg158; reg157=reg187+reg157;
    reg56=reg56-reg189; reg185=reg84+reg185; reg118=reg113+reg118; reg108=reg108-reg100; reg143=reg143-reg149;
    reg137=reg129+reg137; reg82=reg95*reg10; reg69=reg69+reg124; reg84=reg85*reg10; reg180=reg180+reg181;
    reg184=reg183-reg184; reg97=reg186-reg97; reg190=reg190+reg191; reg73=reg193-reg73; reg113=reg78*reg10;
    reg72=reg99+reg72; reg154=reg96+reg154; reg162=reg162+reg161; reg96=reg80*reg10; reg172=reg176-reg172;
    reg170=reg171-reg170; reg174=reg175-reg174; reg146=reg146-reg151; reg120=reg120-reg111; reg86=reg105-reg86;
    reg159=reg155+reg159; reg71=reg166+reg71; reg101=reg101-reg104; reg130=reg130-reg132; reg147=reg142+reg147;
    reg103=reg102+reg103; reg107=reg107-reg112; reg182=reg182-reg122; reg62=reg98+reg62; reg121=reg121-reg188;
    reg163=reg156+reg163; reg177=reg177-reg169; reg63=reg63-reg179; reg126=reg106+reg126; reg131=reg136-reg131;
    reg135=reg135+reg134; reg98=reg89*reg10; reg99=reg94*reg10; reg102=reg91*reg10; reg109=reg110-reg109;
    reg79=reg167-reg79; reg138=reg150-reg138; reg144=reg145-reg144; reg152=reg153-reg152; reg160=reg160+reg168;
    reg114=reg115-reg114; reg107=reg107*reg10; reg135=reg135*reg10; reg103=reg103*reg10; reg182=reg182*reg10;
    reg59=reg59*reg10; reg105=ponderation*reg99; reg62=reg62*reg10; reg137=reg137*reg10; reg163=reg163*reg10;
    reg75=reg75*reg10; reg121=reg121*reg10; reg177=reg177*reg10; reg144=reg144*reg10; reg146=reg146*reg10;
    reg86=reg10*reg86; reg160=reg160*reg10; reg120=reg120*reg10; reg131=reg131*reg10; reg110=ponderation*reg127;
    reg159=reg159*reg10; reg79=reg79*reg10; reg133=reg133*reg10; reg71=reg71*reg10; reg101=reg101*reg10;
    reg119=reg119*reg10; reg140=reg140*reg10; reg130=reg130*reg10; reg128=reg128*reg10; reg147=reg147*reg10;
    reg116=reg116*reg10; reg115=ponderation*reg113; reg56=reg56*reg10; reg117=ponderation*reg102; reg190=reg190*reg10;
    reg126=reg126*reg10; reg73=reg73*reg10; reg97=reg97*reg10; reg185=reg185*reg10; reg184=reg184*reg10;
    reg109=reg109*reg10; reg180=reg180*reg10; reg129=ponderation*reg84; reg108=reg108*reg10; reg118=reg118*reg10;
    reg138=reg138*reg10; reg114=reg114*reg10; reg136=ponderation*reg82; reg69=reg69*reg10; reg123=reg123*reg10;
    reg173=reg173*reg10; reg63=reg63*reg10; reg174=reg174*reg10; reg139=ponderation*reg98; reg170=reg170*reg10;
    reg165=reg165*reg10; reg152=reg152*reg10; reg141=ponderation*reg96; reg172=reg172*reg10; reg157=reg157*reg10;
    reg64=reg64*reg10; reg162=reg162*reg10; reg154=reg154*reg10; reg143=reg143*reg10; reg72=reg72*reg10;
    T tmp_6_3=ponderation*reg162; T tmp_2_1=ponderation*reg73; T tmp_2_4=ponderation*reg131; T tmp_2_5=ponderation*reg123; T tmp_2_0=ponderation*reg172;
    T tmp_2_2=ponderation*reg69; T tmp_2_3=-reg139; T tmp_1_7=ponderation*reg121; T tmp_6_4=ponderation*reg86; T tmp_3_1=ponderation*reg138;
    T tmp_3_2=-reg129; T tmp_3_3=ponderation*reg180; T tmp_3_4=ponderation*reg184; T tmp_3_5=ponderation*reg97; T tmp_3_6=ponderation*reg190;
    T tmp_3_7=-reg115; T tmp_4_0=ponderation*reg72; T tmp_6_5=ponderation*reg128; T tmp_6_2=-reg141; T tmp_6_1=ponderation*reg170;
    T tmp_6_0=ponderation*reg174; T tmp_5_7=ponderation*reg63; T tmp_5_6=ponderation*reg177; T tmp_5_5=ponderation*reg163; T tmp_5_4=ponderation*reg62;
    T tmp_5_3=ponderation*reg182; T tmp_5_2=ponderation*reg107; T tmp_5_1=ponderation*reg103; T tmp_5_0=ponderation*reg147; T tmp_4_7=ponderation*reg130;
    T tmp_4_6=ponderation*reg101; T tmp_4_5=ponderation*reg71; T tmp_4_4=ponderation*reg159; T tmp_4_3=ponderation*reg120; T tmp_4_2=ponderation*reg146;
    T tmp_4_1=ponderation*reg154; T tmp_2_6=-reg110; T tmp_1_6=ponderation*reg133; T tmp_1_5=ponderation*reg140; T tmp_1_4=ponderation*reg116;
    T tmp_1_3=ponderation*reg59; T tmp_1_2=ponderation*reg75; T tmp_1_1=ponderation*reg173; T tmp_1_0=ponderation*reg165; T tmp_0_7=ponderation*reg64;
    T tmp_0_6=ponderation*reg56; T tmp_0_5=ponderation*reg185; T tmp_0_4=ponderation*reg118; T tmp_0_3=ponderation*reg108; T tmp_0_2=ponderation*reg143;
    T tmp_0_1=ponderation*reg137; T tmp_0_0=ponderation*reg119; T tmp_2_7=ponderation*reg160; T tmp_3_0=ponderation*reg79; T tmp_7_7=ponderation*reg135;
    T tmp_7_6=-reg105; T tmp_7_5=ponderation*reg144; T tmp_7_4=ponderation*reg152; T tmp_6_6=ponderation*reg157; T tmp_7_3=-reg117;
    T tmp_7_2=ponderation*reg126; T tmp_7_1=ponderation*reg109; T tmp_7_0=ponderation*reg114; T tmp_6_7=-reg136;
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
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg0*reg2; T reg6=1.0/(*f.m).elastic_modulus; T reg7=reg1*elem.pos(0)[0]; T reg8=reg1*elem.pos(1)[0];
    T reg9=reg1*elem.pos(0)[1]; T reg10=reg1*elem.pos(1)[1]; T reg11=reg3*elem.pos(0)[0]; T reg12=elem.pos(1)[0]*var_inter[0]; T reg13=reg3*elem.pos(0)[1];
    T reg14=elem.pos(1)[1]*var_inter[0]; T reg15=reg6*reg5; reg5=reg4*reg5; T reg16=reg4*reg2; reg2=reg6*reg2;
    reg8=reg8-reg7; T reg17=elem.pos(2)[0]*var_inter[1]; reg10=reg10-reg9; T reg18=reg13+reg14; T reg19=elem.pos(2)[0]*var_inter[0];
    T reg20=reg11+reg12; T reg21=elem.pos(2)[1]*var_inter[0]; T reg22=elem.pos(2)[1]*var_inter[1]; T reg23=reg6*reg15; T reg24=elem.pos(3)[1]*reg3;
    T reg25=reg4*reg5; reg15=reg4*reg15; T reg26=reg6*reg2; T reg27=reg4*reg16; T reg28=elem.pos(3)[0]*reg3;
    reg2=reg4*reg2; reg19=reg19-reg20; reg22=reg10+reg22; reg10=elem.pos(3)[1]*var_inter[1]; T reg29=elem.pos(3)[0]*var_inter[1];
    reg21=reg21-reg18; reg17=reg8+reg17; reg16=reg6*reg16; reg26=reg26-reg27; reg28=reg19+reg28;
    reg2=reg27+reg2; reg5=reg6*reg5; reg15=reg25+reg15; reg23=reg23-reg25; reg22=reg22-reg10;
    reg17=reg17-reg29; reg21=reg24+reg21; reg5=reg25+reg5; reg8=reg6*reg23; reg26=reg6*reg26;
    reg19=reg22*reg28; reg24=reg17*reg21; reg2=reg4*reg2; reg25=reg4*reg15; T reg30=reg6*reg0;
    T reg31=reg27+reg16; T reg32=reg4*reg0; T reg33=pow(reg6,2); reg19=reg24-reg19; reg24=reg4*reg32;
    T reg34=pow(reg4,2); T reg35=reg4*reg5; reg25=reg8-reg25; reg31=reg4*reg31; reg8=reg6*reg30;
    reg2=reg26-reg2; reg28=reg28/reg19; reg22=reg22/reg19; reg35=reg25-reg35; reg31=reg2-reg31;
    reg17=reg17/reg19; reg24=reg8-reg24; reg34=reg33-reg34; reg21=reg21/reg19; reg2=var_inter[1]*reg28;
    reg8=reg3*reg17; reg25=var_inter[0]*reg22; reg26=reg3*reg22; reg33=reg1*reg21; reg34=reg34/reg24;
    T reg36=var_inter[1]*reg21; T reg37=1-(*f.m).resolution; reg31=reg31/reg35; reg23=reg23/reg35; reg15=reg15/reg35;
    reg30=reg30/reg24; T reg38=reg1*reg28; T reg39=reg31*reg37; T reg40=reg36+reg26; T reg41=reg2+reg8;
    reg24=reg32/reg24; reg32=reg25+reg33; T reg42=var_inter[0]*reg17; T reg43=(*f.m).resolution*reg34; T reg44=0.5*reg41;
    T reg45=reg15*reg37; T reg46=reg26-reg33; T reg47=reg36-reg25; T reg48=reg42-reg2; T reg49=0.5*reg40;
    T reg50=reg38-reg8; reg39=reg43+reg39; reg43=(*f.m).resolution*reg24; T reg51=0.5*reg32; T reg52=reg30*(*f.m).resolution;
    T reg53=reg42+reg38; T reg54=reg23*reg37; reg45=reg43+reg45; reg43=0.5*reg48; T reg55=reg51*reg39;
    T reg56=0.5*reg47; T reg57=reg44*reg39; reg52=reg54+reg52; reg54=reg49*reg39; T reg58=0.5*reg50;
    T reg59=0.5*reg53; T reg60=0.5*reg46; T reg61=reg39*reg58; T reg62=reg39*reg60; T reg63=reg59*reg39;
    T reg64=reg53*reg45; reg55=2*reg55; T reg65=reg43*reg39; T reg66=reg56*reg39; T reg67=reg40*reg52;
    reg57=2*reg57; T reg68=reg41*reg45; T reg69=2*reg54; T reg70=reg41*reg52; T reg71=reg40*reg68;
    T reg72=reg52*reg46; reg61=2*reg61; T reg73=reg45*reg50; reg62=2*reg62; T reg74=reg32*reg52;
    T reg75=2*reg63; T reg76=reg47*reg52; reg65=2*reg65; T reg77=reg48*reg45; reg66=2*reg66;
    T reg78=reg52*reg50; T reg79=reg32*reg45; T reg80=reg53*reg52; T reg81=reg47*reg45; T reg82=reg48*reg52;
    T reg83=reg40*reg45; T reg84=reg70*reg53; T reg85=reg51*reg69; T reg86=reg59*reg55; T reg87=reg32*reg64;
    T reg88=reg32*reg67; T reg89=reg59*reg57; T reg90=reg44*reg69; T reg91=reg70*reg48; T reg92=reg56*reg69;
    T reg93=reg74*reg46; T reg94=reg80*reg50; T reg95=reg32*reg74; T reg96=reg75*reg58; T reg97=reg59*reg69;
    T reg98=reg70*reg50; T reg99=reg47*reg76; T reg100=reg43*reg65; T reg101=reg64*reg46; reg86=reg87+reg86;
    T reg102=reg55*reg58; T reg103=reg75*reg60; T reg104=reg67*reg46; T reg105=reg47*reg77; T reg106=reg57*reg58;
    T reg107=reg66*reg58; T reg108=reg62*reg60; T reg109=reg59*reg75; T reg110=reg44*reg57; T reg111=reg56*reg66;
    T reg112=reg77*reg46; T reg113=reg59*reg66; T reg114=reg43*reg57; T reg115=reg47*reg67; T reg116=reg32*reg68;
    T reg117=reg48*reg82; T reg118=reg65*reg58; T reg119=reg78*reg50; T reg120=reg79*reg50; T reg121=reg76*reg46;
    T reg122=reg43*reg66; T reg123=reg40*reg67; T reg124=reg68*reg46; T reg125=reg53*reg81; T reg126=reg51*reg65;
    T reg127=reg56*reg57; T reg128=reg57*reg60; reg68=reg47*reg68; T reg129=reg53*reg80; T reg130=reg82*reg50;
    T reg131=reg66*reg60; reg89=reg88+reg89; reg70=reg70*reg41; T reg132=reg51*reg55; reg71=reg90+reg71;
    T reg133=reg49*reg69; T reg134=reg48*reg83; T reg135=reg83*reg50; T reg136=reg69*reg58; reg84=reg85+reg84;
    T reg137=reg55*reg60; T reg138=reg62*reg58; T reg139=reg32*reg76; T reg140=reg59*reg65; T reg141=reg73*reg46;
    T reg142=reg69*reg60; T reg143=reg53*reg83; T reg144=reg51*reg57; T reg145=reg61*reg58; T reg146=reg43*reg69;
    T reg147=reg81*reg50; T reg148=reg72*reg46; T reg149=reg65*reg60; T reg150=reg53*reg82; T reg151=reg51*reg66;
    T reg152=reg32*reg77; reg98=reg98-reg142; reg111=reg117+reg111; reg108=reg119+reg108; reg120=reg120-reg103;
    reg68=reg68-reg146; reg137=reg137-reg94; reg116=reg116+reg97; reg149=reg147+reg149; reg128=reg128-reg135;
    reg131=reg130+reg131; reg70=reg70+reg133; reg117=reg89*reg19; reg119=reg71*reg19; reg132=reg132+reg129;
    reg113=reg152-reg113; reg110=reg123+reg110; reg125=reg126-reg125; reg150=reg151-reg150; reg140=reg139-reg140;
    reg145=reg148+reg145; reg144=reg144+reg143; reg138=reg141+reg138; reg126=reg84*reg19; reg93=reg93-reg96;
    reg124=reg124-reg136; reg127=reg127-reg134; reg106=reg106-reg104; reg91=reg91-reg92; reg95=reg95+reg109;
    reg107=reg112+reg107; reg114=reg114-reg115; reg118=reg121+reg118; reg122=reg105+reg122; reg102=reg102-reg101;
    reg100=reg99+reg100; reg99=reg86*reg19; reg95=reg95*reg19; reg140=reg140*reg19; reg128=reg128*reg19;
    reg98=reg98*reg19; reg113=reg113*reg19; reg105=ponderation*reg117; reg112=ponderation*reg99; reg107=reg107*reg19;
    reg118=reg118*reg19; reg102=reg102*reg19; reg93=reg93*reg19; reg122=reg122*reg19; reg138=reg138*reg19;
    reg145=reg145*reg19; reg116=reg116*reg19; reg100=reg100*reg19; reg70=reg70*reg19; reg121=ponderation*reg126;
    reg110=reg110*reg19; reg130=ponderation*reg119; reg132=reg132*reg19; reg144=reg144*reg19; reg125=reg125*reg19;
    reg150=reg150*reg19; reg114=reg114*reg19; reg131=reg131*reg19; reg149=reg149*reg19; reg68=reg68*reg19;
    reg137=reg137*reg19; reg111=reg111*reg19; reg120=reg120*reg19; reg108=reg108*reg19; reg91=reg91*reg19;
    reg106=reg106*reg19; reg124=reg124*reg19; reg127=reg127*reg19; T tmp_5_5=ponderation*reg111; T tmp_1_7=ponderation*reg98;
    T tmp_3_7=-reg121; T tmp_2_2=ponderation*reg95; T tmp_5_6=ponderation*reg127; T tmp_4_7=ponderation*reg68; T tmp_2_3=-reg112;
    T tmp_3_6=ponderation*reg144; T tmp_4_5=ponderation*reg122; T tmp_4_4=ponderation*reg100; T tmp_4_6=ponderation*reg114; T tmp_5_7=ponderation*reg91;
    T tmp_3_5=ponderation*reg150; T tmp_1_5=ponderation*reg131; T tmp_1_4=ponderation*reg149; T tmp_1_3=ponderation*reg137; T tmp_1_6=ponderation*reg128;
    T tmp_1_2=ponderation*reg120; T tmp_1_1=ponderation*reg108; T tmp_2_6=-reg105; T tmp_0_7=ponderation*reg124; T tmp_0_6=ponderation*reg106;
    T tmp_0_5=ponderation*reg107; T tmp_0_4=ponderation*reg118; T tmp_0_3=ponderation*reg102; T tmp_0_2=ponderation*reg93; T tmp_0_1=ponderation*reg138;
    T tmp_0_0=ponderation*reg145; T tmp_2_7=ponderation*reg116; T tmp_7_7=ponderation*reg70; T tmp_6_6=ponderation*reg110; T tmp_6_7=-reg130;
    T tmp_2_5=ponderation*reg113; T tmp_3_3=ponderation*reg132; T tmp_3_4=ponderation*reg125; T tmp_2_4=ponderation*reg140;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg3*reg1; reg1=reg2*reg1; T reg6=reg3*reg4; reg4=reg2*reg4;
    T reg7=reg3*reg6; T reg8=reg2*reg4; reg6=reg2*reg6; T reg9=reg3*reg5; T reg10=reg2*reg1;
    reg5=reg2*reg5; reg7=reg7-reg8; reg6=reg8+reg6; reg4=reg3*reg4; reg1=reg3*reg1;
    reg5=reg10+reg5; reg9=reg9-reg10; reg4=reg8+reg4; reg9=reg3*reg9; reg5=reg2*reg5;
    reg8=reg3*reg7; T reg11=reg10+reg1; T reg12=reg2*reg6; reg11=reg2*reg11; reg5=reg9-reg5;
    reg9=reg2*reg4; reg12=reg8-reg12; reg11=reg5-reg11; reg9=reg12-reg9; reg11=reg11/reg9;
    reg6=reg6/reg9; reg7=reg7/reg9; reg5=reg7*reg11; reg8=reg6*reg11; reg12=reg6*(*f.m).alpha;
    T reg13=reg7*(*f.m).alpha; T reg14=reg6*reg8; T reg15=1-var_inter[1]; reg9=reg4/reg9; reg4=1-var_inter[0];
    T reg16=reg7*reg5; reg12=reg13+reg12; reg9=reg9*(*f.m).alpha; reg13=reg4*elem.pos(0)[1]; T reg17=elem.pos(1)[1]*var_inter[0];
    T reg18=elem.pos(1)[0]*var_inter[0]; T reg19=reg4*elem.pos(0)[0]; T reg20=reg15*elem.pos(0)[0]; reg14=reg16-reg14; reg16=reg15*elem.pos(1)[1];
    T reg21=reg15*elem.pos(0)[1]; T reg22=reg15*elem.pos(1)[0]; T reg23=elem.pos(2)[1]*var_inter[1]; reg16=reg16-reg21; T reg24=elem.pos(2)[0]*var_inter[1];
    reg22=reg22-reg20; T reg25=reg19+reg18; T reg26=elem.pos(2)[0]*var_inter[0]; reg9=reg12+reg9; reg12=reg13+reg17;
    T reg27=elem.pos(2)[1]*var_inter[0]; T reg28=reg3*reg0; reg8=reg8/reg14; T reg29=reg2*reg0; reg14=reg5/reg14;
    reg14=reg14*reg9; reg9=reg8*reg9; reg5=reg3*reg28; reg8=reg2*reg29; reg27=reg27-reg12;
    T reg30=elem.pos(3)[1]*reg4; T reg31=elem.pos(3)[0]*reg4; reg26=reg26-reg25; reg23=reg16+reg23; reg16=elem.pos(3)[1]*var_inter[1];
    T reg32=elem.pos(3)[0]*var_inter[1]; reg24=reg22+reg24; reg23=reg23-reg16; reg31=reg26+reg31; reg22=1-(*f.m).resolution;
    reg9=reg14-reg9; reg24=reg24-reg32; reg27=reg30+reg27; reg8=reg5-reg8; reg29=reg29/reg8;
    reg5=reg23*reg31; reg14=reg24*reg27; reg26=(*f.m).alpha*(*f.m).resolution; reg9=reg22*reg9; reg28=reg28/reg8;
    reg5=reg14-reg5; reg14=(*f.m).resolution*reg29; reg6=reg6*reg22; reg9=reg26+reg9; reg7=reg7*reg22;
    reg26=reg28*(*f.m).resolution; reg9=(*f.m).deltaT*reg9; reg6=reg14+reg6; reg27=reg27/reg5; reg26=reg7+reg26;
    reg24=reg24/reg5; reg31=reg31/reg5; reg23=reg23/reg5; reg7=reg15*reg31; reg14=reg4*reg23;
    reg30=reg6*reg9; T reg33=reg26*reg9; T reg34=var_inter[1]*reg27; T reg35=var_inter[0]*reg24; T reg36=reg4*reg24;
    T reg37=reg15*reg27; T reg38=reg33+reg30; T reg39=var_inter[0]*reg23; T reg40=reg35+reg7; T reg41=reg4*var_inter[1];
    T reg42=reg15*var_inter[0]; T reg43=var_inter[1]*reg31; T reg44=reg34+reg14; T reg45=reg7-reg36; T reg46=reg41*(*f.m).f_vol[0];
    T reg47=reg14-reg37; T reg48=reg43+reg36; T reg49=reg38*reg44; T reg50=reg38*reg40; T reg51=var_inter[0]*var_inter[1];
    T reg52=reg15*reg4; T reg53=reg35-reg43; T reg54=reg34-reg39; T reg55=reg42*(*f.m).f_vol[1]; T reg56=reg39+reg37;
    T reg57=reg50-reg55; T reg58=reg38*reg54; T reg59=reg38*reg53; T reg60=reg38*reg56; T reg61=reg38*reg45;
    T reg62=reg38*reg47; T reg63=reg41*(*f.m).f_vol[1]; T reg64=reg51*(*f.m).f_vol[1]; T reg65=reg51*(*f.m).f_vol[0]; T reg66=reg42*(*f.m).f_vol[0];
    T reg67=reg52*(*f.m).f_vol[1]; T reg68=reg38*reg48; T reg69=reg52*(*f.m).f_vol[0]; T reg70=reg49-reg46; T reg71=reg68+reg63;
    T reg72=reg62+reg69; T reg73=reg61+reg67; reg70=reg70*reg5; T reg74=reg59+reg64; T reg75=reg58+reg65;
    reg57=reg57*reg5; T reg76=reg60+reg66; T reg77=reg74*reg5; T reg78=reg75*reg5; T reg79=reg76*reg5;
    T reg80=reg71*reg5; reg57=ponderation*reg57; T reg81=reg72*reg5; reg70=ponderation*reg70; T reg82=reg73*reg5;
    sollicitation[indices[3]+0]+=-reg70; reg70=ponderation*reg82; sollicitation[indices[0]+1]+=reg70; sollicitation[indices[1]+1]+=-reg57; reg57=ponderation*reg80;
    sollicitation[indices[3]+1]+=reg57; T reg83=ponderation*reg77; sollicitation[indices[2]+1]+=reg83; T reg84=ponderation*reg81; sollicitation[indices[0]+0]+=reg84;
    T reg85=ponderation*reg78; sollicitation[indices[2]+0]+=reg85; T reg86=ponderation*reg79; sollicitation[indices[1]+0]+=reg86;
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
    T reg4=reg0*reg1; T reg5=reg3*reg1; T reg6=reg3*reg4; reg4=reg2*reg4; reg1=reg2*reg1;
    T reg7=reg3*reg4; T reg8=reg3*reg6; reg4=reg2*reg4; T reg9=reg2*reg1; T reg10=reg3*reg5;
    reg1=reg3*reg1; reg6=reg2*reg6; reg7=reg8+reg7; reg4=reg4-reg8; reg5=reg2*reg5;
    reg9=reg9-reg10; reg1=reg10+reg1; reg1=reg3*reg1; reg9=reg2*reg9; T reg11=reg2*reg4;
    reg6=reg8+reg6; reg8=reg3*reg7; T reg12=reg10+reg5; reg8=reg11-reg8; reg11=reg3*reg6;
    reg12=reg3*reg12; reg1=reg9-reg1; reg12=reg1-reg12; reg11=reg8-reg11; reg1=1-var_inter[0];
    reg8=1-var_inter[1]; reg9=reg8*elem.pos(0)[1]; T reg13=elem.pos(1)[1]*var_inter[0]; reg4=reg4/reg11; reg12=reg12/reg11;
    reg7=reg7/reg11; T reg14=reg1*elem.pos(0)[1]; T reg15=reg8*elem.pos(1)[1]; T reg16=reg8*elem.pos(1)[0]; T reg17=reg8*elem.pos(0)[0];
    T reg18=reg1*elem.pos(0)[0]; T reg19=elem.pos(1)[0]*var_inter[0]; T reg20=elem.pos(2)[0]*var_inter[1]; reg16=reg16-reg17; T reg21=reg7*reg12;
    T reg22=reg4*reg12; T reg23=reg18+reg19; T reg24=elem.pos(2)[1]*var_inter[1]; reg15=reg15-reg9; T reg25=elem.pos(2)[0]*var_inter[0];
    T reg26=elem.pos(2)[1]*var_inter[0]; T reg27=reg14+reg13; T reg28=reg4*reg22; T reg29=reg7*reg21; T reg30=reg4*(*f.m).alpha;
    T reg31=reg7*(*f.m).alpha; reg24=reg15+reg24; reg26=reg26-reg27; reg11=reg6/reg11; reg6=elem.pos(3)[1]*reg1;
    reg15=elem.pos(3)[0]*reg1; T reg32=elem.pos(3)[1]*var_inter[1]; reg25=reg25-reg23; reg20=reg16+reg20; reg16=elem.pos(3)[0]*var_inter[1];
    T reg33=reg8*vectors[0][indices[0]+0]; reg24=reg24-reg32; reg26=reg6+reg26; reg29=reg28-reg29; reg15=reg25+reg15;
    reg31=reg30+reg31; reg11=reg11*(*f.m).alpha; reg6=var_inter[0]*vectors[0][indices[1]+1]; reg25=reg1*vectors[0][indices[0]+1]; reg28=reg8*vectors[0][indices[1]+1];
    reg30=reg8*vectors[0][indices[0]+1]; reg20=reg20-reg16; T reg34=reg1*vectors[0][indices[0]+0]; T reg35=var_inter[0]*vectors[0][indices[1]+0]; T reg36=reg8*vectors[0][indices[1]+0];
    T reg37=var_inter[1]*vectors[0][indices[2]+1]; reg22=reg22/reg29; reg35=reg34+reg35; reg29=reg21/reg29; reg30=reg28-reg30;
    reg21=reg20*reg26; reg25=reg6+reg25; reg6=reg24*reg15; reg28=var_inter[0]*vectors[0][indices[2]+1]; reg11=reg31+reg11;
    reg33=reg36-reg33; reg31=var_inter[1]*vectors[0][indices[2]+0]; reg34=var_inter[0]*vectors[0][indices[2]+0]; reg36=reg3*reg0; T reg38=reg1*vectors[0][indices[3]+0];
    T reg39=reg2*reg0; T reg40=reg1*vectors[0][indices[3]+1]; reg25=reg28-reg25; reg22=reg22*reg11; reg11=reg29*reg11;
    reg31=reg33+reg31; reg28=var_inter[1]*vectors[0][indices[3]+0]; reg30=reg37+reg30; reg29=var_inter[1]*vectors[0][indices[3]+1]; reg35=reg34-reg35;
    reg6=reg21-reg6; reg26=reg26/reg6; reg40=reg25+reg40; reg21=1-(*f.m).resolution; reg20=reg20/reg6;
    reg25=pow(reg3,2); reg11=reg22-reg11; reg22=pow(reg2,2); reg33=reg3*reg36; reg34=reg2*reg39;
    reg29=reg30-reg29; reg38=reg35+reg38; reg15=reg15/reg6; reg28=reg31-reg28; reg24=reg24/reg6;
    reg11=reg21*reg11; reg30=(*f.m).alpha*(*f.m).resolution; reg33=reg34-reg33; reg25=reg22-reg25; reg22=reg40*reg24;
    reg31=reg28*reg15; reg34=reg38*reg20; reg35=reg29*reg26; reg22=reg35-reg22; reg31=reg34-reg31;
    reg11=reg30+reg11; reg38=reg38*reg24; reg28=reg28*reg26; reg39=reg39/reg33; reg40=reg40*reg20;
    reg29=reg29*reg15; reg25=reg25/reg33; reg33=reg36/reg33; reg11=(*f.m).deltaT*reg11; reg29=reg40-reg29;
    reg31=reg22+reg31; reg38=reg28-reg38; reg22=reg39*(*f.m).resolution; reg4=reg4*reg21; reg28=(*f.m).resolution*reg25;
    reg12=reg12*reg21; reg21=reg7*reg21; reg7=(*f.m).resolution*reg33; reg38=reg38-reg11; reg30=var_inter[0]*reg24;
    reg29=reg29-reg11; reg34=var_inter[0]*reg20; reg31=0.5*reg31; reg35=reg1*reg24; reg36=reg8*reg15;
    reg37=reg1*reg20; reg40=var_inter[1]*reg15; reg22=reg4+reg22; reg21=reg7+reg21; reg12=reg28+reg12;
    reg4=var_inter[1]*reg26; reg7=reg8*reg26; reg28=reg4-reg30; reg31=reg12*reg31; T reg41=reg34-reg40;
    T reg42=reg22*reg29; T reg43=reg21*reg38; reg29=reg21*reg29; T reg44=reg4+reg35; T reg45=reg40+reg37;
    T reg46=reg30+reg7; T reg47=reg34+reg36; reg38=reg22*reg38; T reg48=reg35-reg7; T reg49=reg36-reg37;
    reg31=2*reg31; T reg50=0.5*reg47; T reg51=0.5*reg41; reg42=reg43+reg42; reg43=0.5*reg49;
    reg29=reg38+reg29; reg38=0.5*reg28; T reg52=0.5*reg46; T reg53=0.5*reg45; T reg54=0.5*reg44;
    T reg55=0.5*reg48; T reg56=reg31*reg43; T reg57=reg1*var_inter[1]; T reg58=reg8*reg1; T reg59=var_inter[0]*var_inter[1];
    T reg60=reg42*reg49; T reg61=reg31*reg55; T reg62=reg29*reg48; T reg63=reg8*var_inter[0]; T reg64=reg28*reg29;
    T reg65=reg51*reg31; T reg66=reg41*reg42; T reg67=reg38*reg31; T reg68=reg44*reg29; T reg69=reg53*reg31;
    T reg70=reg45*reg42; T reg71=reg54*reg31; T reg72=reg52*reg31; T reg73=reg47*reg42; T reg74=reg46*reg29;
    T reg75=reg50*reg31; T reg76=reg63*(*f.m).f_vol[1]; reg65=reg64+reg65; reg64=reg59*(*f.m).f_vol[0]; reg74=reg74-reg75;
    T reg77=reg63*(*f.m).f_vol[0]; reg67=reg66+reg67; reg66=reg59*(*f.m).f_vol[1]; reg69=reg69-reg68; T reg78=reg57*(*f.m).f_vol[0];
    reg56=reg62+reg56; reg70=reg70-reg71; reg61=reg60+reg61; reg72=reg72-reg73; reg60=reg58*(*f.m).f_vol[1];
    reg62=reg58*(*f.m).f_vol[0]; T reg79=reg57*(*f.m).f_vol[1]; reg56=reg56-reg62; reg69=reg69-reg78; reg72=reg72-reg76;
    reg70=reg70-reg79; reg74=reg74-reg77; reg65=reg65-reg64; reg67=reg67-reg66; reg61=reg61-reg60;
    reg72=reg72*reg6; reg56=reg56*reg6; reg65=reg65*reg6; reg61=reg61*reg6; reg67=reg67*reg6;
    reg74=reg74*reg6; reg69=reg69*reg6; reg70=reg70*reg6; sollicitation[indices[0]+1]+=ponderation*reg61; sollicitation[indices[2]+0]+=ponderation*reg65;
    sollicitation[indices[2]+1]+=ponderation*reg67; sollicitation[indices[3]+1]+=ponderation*reg70; sollicitation[indices[0]+0]+=ponderation*reg56; sollicitation[indices[1]+1]+=ponderation*reg72; sollicitation[indices[3]+0]+=ponderation*reg69;
    sollicitation[indices[1]+0]+=ponderation*reg74;
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

