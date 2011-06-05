
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg1*reg3; reg1=reg1*reg2; T reg6=reg4*reg3; reg4=reg4*reg2;
    T reg7=reg2*reg1; T reg8=reg3*reg5; T reg9=reg6*reg2; reg5=reg2*reg5; reg6=reg6*reg3;
    T reg10=reg4*reg2; reg9=reg9+reg10; reg6=reg6-reg10; reg5=reg7+reg5; reg8=reg8-reg7;
    reg4=reg4*reg3; reg1=reg3*reg1; reg5=reg2*reg5; T reg11=reg3*reg6; T reg12=reg7+reg1;
    reg10=reg4+reg10; reg4=reg2*reg9; reg8=reg3*reg8; T reg13=reg2*reg10; T reg14=elem.pos(2)[0]-elem.pos(0)[0];
    T reg15=elem.pos(2)[1]-elem.pos(0)[1]; T reg16=elem.pos(1)[1]-elem.pos(0)[1]; reg4=reg11-reg4; reg11=elem.pos(1)[0]-elem.pos(0)[0]; reg5=reg8-reg5;
    reg12=reg2*reg12; reg8=reg14*reg16; reg12=reg5-reg12; reg5=reg15*reg11; reg13=reg4-reg13;
    reg8=reg5-reg8; reg12=reg12/reg13; reg9=reg9/reg13; reg6=reg6/reg13; reg4=reg12*reg9;
    reg5=reg0*reg3; T reg17=reg12*reg6; T reg18=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; T reg19=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; T reg20=vectors[0][indices[1]+0]-vectors[0][indices[0]+0];
    reg16=reg16/reg8; reg14=reg14/reg8; T reg21=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; reg11=reg11/reg8; reg15=reg15/reg8;
    T reg22=reg0*reg2; T reg23=(*f.m).alpha*reg9; T reg24=reg20*reg15; T reg25=reg19*reg11; T reg26=reg16*reg21;
    T reg27=reg14*reg18; T reg28=reg17*reg6; T reg29=reg4*reg9; reg13=reg10/reg13; reg10=reg2*reg22;
    T reg30=reg3*reg5; T reg31=(*f.m).alpha*reg6; reg27=reg25-reg27; elem.epsilon[0][1]=reg27; reg13=(*f.m).alpha*reg13;
    reg31=reg23+reg31; reg10=reg30-reg10; reg26=reg24-reg26; elem.epsilon[0][0]=reg26; reg29=reg28-reg29;
    reg23=(*f.m).deltaT*(*f.m).alpha; reg22=reg22/reg10; reg13=reg31+reg13; reg4=reg4/reg29; reg24=reg27-reg23;
    reg25=reg26-reg23; reg5=reg5/reg10; reg29=reg17/reg29; reg29=reg13*reg29; reg4=reg13*reg4;
    reg13=reg24*reg5; reg17=reg25*reg5; reg24=reg22*reg24; reg25=reg22*reg25; reg13=reg25+reg13;
    reg24=reg17+reg24; reg17=1-(*f.m).resolution; reg4=reg29-reg4; reg25=(*f.m).resolution*(*f.m).alpha; reg28=reg2*reg13;
    reg29=reg2*reg24; reg13=reg3*reg13; reg30=PNODE(2).dep[0]-PNODE(0).dep[0]; reg31=PNODE(1).dep[0]-PNODE(0).dep[0]; T reg32=pow(reg2,2);
    T reg33=PNODE(2).dep[1]-PNODE(0).dep[1]; T reg34=PNODE(1).dep[1]-PNODE(0).dep[1]; T reg35=pow(reg3,2); reg24=reg3*reg24; reg4=reg4*reg17;
    reg13=reg13-reg29; reg24=reg24-reg28; T reg36=reg33*reg11; T reg37=reg11*reg30; T reg38=reg15*reg31;
    T reg39=reg34*reg15; reg30=reg16*reg30; reg33=reg33*reg16; reg34=reg34*reg14; reg32=reg35-reg32;
    reg31=reg14*reg31; reg25=reg4+reg25; reg30=reg38-reg30; reg31=reg37-reg31; reg25=(*f.m).deltaT*reg25;
    reg13=reg23+reg13; reg18=reg18*reg15; reg24=reg23+reg24; reg29=reg28+reg29; reg10=reg32/reg10;
    reg20=reg14*reg20; reg21=reg21*reg11; reg33=reg39-reg33; reg34=reg36-reg34; reg19=reg16*reg19;
    reg9=reg17*reg9; reg6=reg17*reg6; reg4=(*f.m).resolution*reg5; reg28=(*f.m).resolution*reg22; reg19=reg18-reg19;
    reg18=reg34-reg25; reg32=reg30-reg25; reg35=reg13+reg24; reg20=reg21-reg20; reg21=(*f.m).resolution*reg10;
    reg29=reg23-reg29; reg12=reg17*reg12; reg33=reg31+reg33; reg6=reg4+reg6; reg28=reg9+reg28;
    reg4=reg28*reg32; reg9=reg6*reg18; reg17=reg28*reg18; reg33=0.5*reg33; reg12=reg21+reg12;
    reg21=reg6*reg32; reg35=reg35+reg29; reg19=reg20+reg19; reg21=reg17+reg21; reg35=reg35/3;
    reg19=0.5*reg19; elem.epsilon[0][2]=reg19; reg17=reg33*reg12; reg9=reg4+reg9; reg18=reg9*reg18;
    reg32=reg21*reg32; reg17=2*reg17; reg24=reg24-reg35; reg4=reg10*reg19; reg13=reg13-reg35;
    reg24=pow(reg24,2); reg4=reg0*reg4; reg18=reg32+reg18; reg20=reg33*reg17; reg35=reg29-reg35;
    reg13=pow(reg13,2); reg35=pow(reg35,2); reg18=reg20+reg18; reg13=reg24+reg13; reg20=2*reg4;
    reg35=reg13+reg35; reg20=reg4*reg20; reg18=reg8*reg18; reg4=0.33333333333333331483*reg18; reg27=reg27-reg25;
    reg18=0.16666666666666665741*reg18; reg20=reg35+reg20; reg25=reg26-reg25; reg8=reg27*reg28; reg4=reg18+reg4;
    reg13=reg25*reg6; reg28=reg25*reg28; reg6=reg27*reg6; reg20=1.5*reg20; elem.ener=reg4/2;
    elem.sigma[0][0]=reg13+reg8; elem.sigma[0][1]=reg28+reg6; elem.sigma_von_mises=pow(reg20,0.5); elem.sigma[0][2]=reg19*reg12;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg1*reg4; T reg6=reg2*reg4; reg1=reg1*reg3; reg2=reg2*reg3;
    T reg7=reg3*reg1; reg1=reg4*reg1; T reg8=reg4*reg5; T reg9=reg2*reg3; reg2=reg2*reg4;
    T reg10=reg6*reg4; reg6=reg6*reg3; reg7=reg7-reg8; reg1=reg8+reg1; reg5=reg3*reg5;
    reg9=reg9-reg10; reg2=reg2+reg10; T reg11=reg4*reg2; reg1=reg4*reg1; T reg12=reg3*reg9;
    T reg13=reg8+reg5; reg7=reg3*reg7; reg10=reg6+reg10; reg6=reg4*reg10; reg1=reg7-reg1;
    reg11=reg12-reg11; reg13=reg4*reg13; reg6=reg11-reg6; reg13=reg1-reg13; reg9=reg9/reg6;
    reg2=reg2/reg6; reg13=reg13/reg6; reg1=reg13*reg2; reg7=reg13*reg9; reg11=reg7*reg9;
    reg12=reg1*reg2; T reg14=(*f.m).alpha*reg9; reg6=reg10/reg6; reg10=(*f.m).alpha*reg2; reg12=reg11-reg12;
    reg6=(*f.m).alpha*reg6; reg14=reg10+reg14; reg7=reg7/reg12; reg10=reg0*reg3; reg12=reg1/reg12;
    reg1=reg0*reg4; reg6=reg14+reg6; reg12=reg6*reg12; reg7=reg6*reg7; reg6=elem.pos(2)[1]-elem.pos(0)[1];
    reg11=elem.pos(2)[0]-elem.pos(0)[0]; reg14=pow(reg4,2); T reg15=reg3*reg10; T reg16=elem.pos(1)[0]-elem.pos(0)[0]; T reg17=reg4*reg1;
    T reg18=elem.pos(1)[1]-elem.pos(0)[1]; T reg19=pow(reg3,2); T reg20=1-(*f.m).resolution; reg12=reg7-reg12; reg7=reg6*reg16;
    T reg21=reg11*reg18; reg17=reg15-reg17; reg14=reg19-reg14; reg14=reg14/reg17; reg21=reg7-reg21;
    reg12=reg12*reg20; reg1=reg1/reg17; reg17=reg10/reg17; reg7=(*f.m).resolution*(*f.m).alpha; reg6=reg6/reg21;
    reg11=reg11/reg21; reg18=reg18/reg21; reg10=(*f.m).resolution*reg14; reg15=(*f.m).resolution*reg1; reg13=reg20*reg13;
    reg19=(*f.m).resolution*reg17; reg9=reg20*reg9; reg7=reg12+reg7; reg2=reg20*reg2; reg16=reg16/reg21;
    reg12=0.5*reg18; reg20=0.5*reg16; T reg22=reg18-reg6; T reg23=reg11-reg16; T reg24=0.5*reg6;
    T reg25=0.5*reg11; reg7=(*f.m).deltaT*reg7; reg9=reg19+reg9; reg13=reg10+reg13; reg15=reg2+reg15;
    reg2=reg15*reg7; reg10=0.5*reg23; reg19=reg9*reg7; T reg26=0.5*reg22; T reg27=reg24*reg13;
    T reg28=reg12*reg13; T reg29=reg20*reg13; T reg30=reg25*reg13; T reg31=1-var_inter[0]; T reg32=reg2+reg19;
    T reg33=reg15*reg6; T reg34=reg18*reg15; reg27=2*reg27; T reg35=reg9*reg16; T reg36=2*reg30;
    T reg37=reg9*reg6; T reg38=2*reg28; T reg39=reg15*reg16; T reg40=reg11*reg9; T reg41=reg11*reg15;
    reg29=2*reg29; T reg42=reg18*reg9; T reg43=reg10*reg13; T reg44=reg26*reg13; T reg45=reg25*reg27;
    T reg46=reg32*reg11; T reg47=reg41*reg6; T reg48=reg22*reg15; T reg49=reg36*reg24; T reg50=reg11*reg33;
    T reg51=reg42*reg6; T reg52=reg22*reg9; reg43=2*reg43; T reg53=reg40*reg16; T reg54=reg12*reg27;
    T reg55=reg25*reg29; T reg56=reg34*reg16; T reg57=reg12*reg29; T reg58=reg23*reg9; T reg59=reg11*reg35;
    T reg60=var_inter[0]*(*f.m).f_vol[1]; T reg61=reg24*reg38; T reg62=var_inter[1]*(*f.m).f_vol[0]; T reg63=reg18*reg39; T reg64=reg23*reg15;
    T reg65=reg20*reg38; T reg66=reg37*reg18; T reg67=reg36*reg20; reg31=reg31-var_inter[1]; T reg68=reg32*reg18;
    reg44=2*reg44; T reg69=reg20*reg43; T reg70=reg18*reg52; T reg71=var_inter[1]*(*f.m).f_vol[1]; T reg72=reg39*reg6;
    T reg73=reg25*reg38; T reg74=reg44*reg20; T reg75=reg64*reg18; T reg76=reg24*reg43; T reg77=reg11*reg48;
    T reg78=reg20*reg29; T reg79=reg46-reg60; T reg80=reg18*reg42; T reg81=reg52*reg6; T reg82=reg18*reg41;
    T reg83=reg26*reg36; T reg84=reg37*reg22; T reg85=reg10*reg36; T reg86=reg36*reg25; reg45=reg47+reg45;
    T reg87=reg23*reg40; T reg88=reg26*reg27; reg37=reg37*reg6; T reg89=reg32*reg16; T reg90=reg10*reg27;
    T reg91=reg23*reg34; T reg92=reg26*reg29; T reg93=reg23*reg35; T reg94=reg22*reg41; T reg95=reg26*reg38;
    T reg96=reg68-reg62; T reg97=reg10*reg29; reg55=reg51+reg55; reg52=reg22*reg52; reg39=reg22*reg39;
    T reg98=reg10*reg43; T reg99=reg24*reg27; T reg100=reg40*reg11; T reg101=reg64*reg6; T reg102=reg48*reg16;
    T reg103=var_inter[0]*(*f.m).f_vol[0]; T reg104=reg12*reg43; reg29=reg24*reg29; reg54=reg53+reg54; T reg105=reg34*reg11;
    T reg106=reg44*reg25; reg59=reg61+reg59; T reg107=reg32*reg23; T reg108=reg36*reg12; T reg109=reg33*reg16;
    T reg110=reg32*reg6; T reg111=reg22*reg42; T reg112=reg44*reg12; T reg113=reg58*reg16; reg64=reg64*reg22;
    reg33=reg23*reg33; reg27=reg20*reg27; T reg114=reg25*reg43; reg66=reg67+reg66; T reg115=reg26*reg44;
    T reg116=reg23*reg58; reg63=reg65+reg63; T reg117=reg31*(*f.m).f_vol[1]; T reg118=reg31*(*f.m).f_vol[0]; T reg119=reg44*reg24;
    reg43=reg26*reg43; reg48=reg23*reg48; reg58=reg11*reg58; T reg120=reg10*reg38; T reg121=reg12*reg38;
    reg35=reg35*reg16; reg44=reg10*reg44; T reg122=reg32*reg22; reg50=reg49+reg50; reg57=reg56+reg57;
    reg97=reg97-reg111; reg39=reg39-reg120; T reg123=reg103+reg110; T reg124=reg117+reg107; reg29=reg29+reg105;
    reg99=reg99+reg100; T reg125=reg21*reg50; T reg126=reg122+reg118; reg58=reg119-reg58; reg119=reg66*reg21;
    reg27=reg27+reg82; reg77=reg76-reg77; reg72=reg72+reg73; reg76=reg21*reg55; reg79=reg79*reg21;
    reg90=reg90-reg94; reg96=reg96*reg21; T reg127=reg21*reg45; T reg128=reg89+reg71; reg37=reg37+reg86;
    reg70=reg69-reg70; reg52=reg98+reg52; reg75=reg74-reg75; reg78=reg80+reg78; reg69=reg21*reg59;
    reg33=reg33-reg83; reg92=reg92-reg91; reg115=reg116+reg115; reg43=reg48+reg43; reg48=reg21*reg63;
    reg35=reg35+reg121; reg88=reg88-reg87; reg74=reg21*reg57; reg84=reg84-reg85; reg64=reg44+reg64;
    reg106=reg101-reg106; reg112=reg113-reg112; reg109=reg109+reg108; reg114=reg81-reg114; reg93=reg93-reg95;
    reg44=reg21*reg54; reg104=reg102-reg104; reg104=reg21*reg104; reg81=reg21*reg123; reg70=reg70*reg21;
    reg72=reg21*reg72; reg99=reg21*reg99; reg75=reg75*reg21; reg98=ponderation*reg44; reg112=reg21*reg112;
    reg78=reg21*reg78; reg27=reg27*reg21; reg101=reg21*reg124; reg102=ponderation*reg119; reg113=ponderation*reg74;
    reg115=reg21*reg115; reg116=ponderation*reg125; reg43=reg21*reg43; T reg129=reg126*reg21; reg77=reg21*reg77;
    reg109=reg21*reg109; reg35=reg21*reg35; T reg130=ponderation*reg48; reg29=reg21*reg29; reg58=reg21*reg58;
    reg37=reg21*reg37; T reg131=reg128*reg21; reg33=reg21*reg33; T reg132=ponderation*reg127; reg64=reg64*reg21;
    reg88=reg88*reg21; reg84=reg84*reg21; reg39=reg21*reg39; reg96=ponderation*reg96; reg92=reg92*reg21;
    reg93=reg21*reg93; reg90=reg21*reg90; T reg133=ponderation*reg76; reg52=reg21*reg52; reg106=reg21*reg106;
    reg79=ponderation*reg79; reg97=reg21*reg97; reg114=reg21*reg114; T reg134=ponderation*reg69; T tmp_3_2=-reg116;
    T tmp_5_0=ponderation*reg104; reg104=ponderation*reg81; sollicitation[indices[1]+0]+=reg104; T tmp_4_5=-reg130; T tmp_5_3=-reg98;
    T tmp_0_1=ponderation*reg64; T tmp_5_1=ponderation*reg112; T tmp_1_5=ponderation*reg93; T tmp_5_2=ponderation*reg109; T tmp_3_3=ponderation*reg99;
    T tmp_2_1=ponderation*reg106; reg64=ponderation*reg131; sollicitation[indices[2]+1]+=reg64; T tmp_0_5=ponderation*reg39; T tmp_0_4=ponderation*reg97;
    reg39=ponderation*reg101; sollicitation[indices[0]+1]+=reg39; T tmp_3_4=ponderation*reg29; T tmp_2_0=ponderation*reg114; T tmp_5_5=ponderation*reg35;
    T tmp_2_4=-reg133; T tmp_3_5=-reg134; T tmp_4_0=ponderation*reg70; T tmp_2_5=ponderation*reg72; sollicitation[indices[1]+1]+=-reg79;
    T tmp_0_0=ponderation*reg52; T tmp_4_1=ponderation*reg75; T tmp_0_3=ponderation*reg90; T tmp_4_3=ponderation*reg27; T tmp_4_4=ponderation*reg78;
    T tmp_1_4=ponderation*reg92; T tmp_1_1=ponderation*reg115; T tmp_4_2=-reg102; T tmp_1_3=ponderation*reg88; T tmp_1_0=ponderation*reg43;
    sollicitation[indices[2]+0]+=-reg96; reg27=ponderation*reg129; sollicitation[indices[0]+0]+=reg27; T tmp_0_2=ponderation*reg84; T tmp_3_0=ponderation*reg77;
    T tmp_2_3=-reg132; T tmp_1_2=ponderation*reg33; T tmp_3_1=ponderation*reg58; T tmp_5_4=-reg113; T tmp_2_2=ponderation*reg37;
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
    T reg4=reg0*reg1; T reg5=reg4*reg3; reg4=reg4*reg2; T reg6=reg1*reg2; reg1=reg1*reg3;
    T reg7=reg3*reg1; T reg8=reg2*reg6; T reg9=reg4*reg2; T reg10=reg5*reg3; reg6=reg3*reg6;
    reg4=reg4*reg3; reg5=reg5*reg2; reg1=reg2*reg1; reg4=reg4+reg10; reg8=reg8-reg7;
    reg6=reg7+reg6; reg9=reg9-reg10; T reg11=reg3*reg4; reg10=reg5+reg10; reg5=reg2*reg9;
    T reg12=reg7+reg1; reg8=reg2*reg8; reg6=reg3*reg6; reg11=reg5-reg11; reg5=reg3*reg10;
    reg6=reg8-reg6; reg12=reg3*reg12; reg5=reg11-reg5; reg12=reg6-reg12; reg4=reg4/reg5;
    reg12=reg12/reg5; reg9=reg9/reg5; reg6=reg12*reg4; reg8=reg12*reg9; reg11=reg8*reg9;
    T reg13=(*f.m).alpha*reg9; reg5=reg10/reg5; reg10=reg6*reg4; T reg14=(*f.m).alpha*reg4; reg5=(*f.m).alpha*reg5;
    reg13=reg14+reg13; reg10=reg11-reg10; reg11=reg0*reg3; reg6=reg6/reg10; reg14=reg0*reg2;
    reg5=reg13+reg5; reg10=reg8/reg10; reg8=pow(reg2,2); reg13=reg3*reg11; T reg15=elem.pos(2)[1]-elem.pos(0)[1];
    T reg16=elem.pos(2)[0]-elem.pos(0)[0]; T reg17=elem.pos(1)[1]-elem.pos(0)[1]; T reg18=reg2*reg14; T reg19=pow(reg3,2); T reg20=elem.pos(1)[0]-elem.pos(0)[0];
    reg10=reg5*reg10; reg6=reg5*reg6; reg5=reg15*reg20; T reg21=1-(*f.m).resolution; reg13=reg18-reg13;
    reg6=reg10-reg6; reg19=reg8-reg19; reg8=reg16*reg17; reg19=reg19/reg13; reg14=reg14/reg13;
    reg8=reg5-reg8; reg6=reg6*reg21; reg5=(*f.m).resolution*(*f.m).alpha; reg13=reg11/reg13; reg12=reg21*reg12;
    reg10=(*f.m).resolution*reg13; reg11=(*f.m).resolution*reg14; reg9=reg21*reg9; reg4=reg21*reg4; reg15=reg15/reg8;
    reg20=reg20/reg8; reg5=reg6+reg5; reg6=(*f.m).resolution*reg19; reg17=reg17/reg8; reg16=reg16/reg8;
    reg18=reg16-reg20; reg21=0.5*reg15; T reg22=reg17-reg15; T reg23=0.5*reg20; reg12=reg6+reg12;
    reg10=reg4+reg10; reg9=reg11+reg9; reg5=(*f.m).deltaT*reg5; reg4=0.5*reg17; reg6=reg9*reg5;
    reg11=reg10*reg5; T reg24=0.5*reg18; T reg25=0.5*reg22; T reg26=reg4*reg12; T reg27=reg23*reg12;
    T reg28=0.5*reg16; T reg29=reg21*reg12; T reg30=reg9*reg20; T reg31=2*reg26; T reg32=1-var_inter[0];
    T reg33=reg11+reg6; T reg34=reg28*reg12; T reg35=reg10*reg20; T reg36=reg17*reg9; reg27=2*reg27;
    T reg37=reg16*reg10; T reg38=reg24*reg12; T reg39=reg25*reg12; reg29=2*reg29; T reg40=reg33*reg16;
    T reg41=var_inter[1]*(*f.m).f_vol[0]; T reg42=var_inter[0]*(*f.m).f_vol[1]; T reg43=reg10*reg15; T reg44=reg16*reg30; T reg45=reg21*reg31;
    T reg46=reg37*reg15; T reg47=reg28*reg27; T reg48=reg28*reg29; reg38=2*reg38; T reg49=reg22*reg9;
    T reg50=reg36*reg15; T reg51=reg18*reg10; reg39=2*reg39; T reg52=reg23*reg31; T reg53=reg17*reg35;
    T reg54=2*reg34; reg32=reg32-var_inter[1]; T reg55=reg9*reg15; T reg56=reg18*reg9; T reg57=reg33*reg17;
    T reg58=reg17*reg10; T reg59=reg16*reg9; T reg60=reg22*reg36; T reg61=reg55*reg22; T reg62=reg24*reg54;
    T reg63=reg23*reg27; T reg64=reg24*reg31; reg44=reg45+reg44; T reg65=reg22*reg35; T reg66=reg33*reg18;
    T reg67=reg55*reg15; T reg68=reg25*reg31; T reg69=reg18*reg30; T reg70=reg54*reg28; reg48=reg46+reg48;
    T reg71=reg21*reg29; T reg72=reg33*reg22; T reg73=reg59*reg16; T reg74=reg32*(*f.m).f_vol[0]; T reg75=reg32*(*f.m).f_vol[1];
    T reg76=reg17*reg36; T reg77=reg40-reg42; T reg78=reg57-reg41; T reg79=reg28*reg31; reg35=reg35*reg15;
    reg47=reg50+reg47; T reg80=reg21*reg27; T reg81=reg33*reg20; T reg82=reg24*reg27; T reg83=reg58*reg16;
    T reg84=reg22*reg37; T reg85=reg24*reg29; T reg86=reg24*reg39; T reg87=reg51*reg22; T reg88=var_inter[0]*(*f.m).f_vol[0];
    T reg89=reg33*reg15; reg30=reg30*reg20; T reg90=reg4*reg31; reg53=reg52+reg53; T reg91=reg18*reg56;
    T reg92=reg25*reg27; T reg93=reg18*reg58; T reg94=reg25*reg39; T reg95=var_inter[1]*(*f.m).f_vol[1]; T reg96=reg24*reg38;
    T reg97=reg25*reg54; T reg98=reg25*reg29; T reg99=reg18*reg59; T reg100=reg22*reg49; T reg101=reg18*reg43;
    reg85=reg85-reg84; reg80=reg80+reg83; T reg102=reg8*reg53; T reg103=reg81+reg95; reg100=reg96+reg100;
    reg96=reg8*reg47; reg35=reg35+reg79; reg78=reg78*reg8; reg30=reg30+reg90; reg71=reg71+reg73;
    reg77=reg77*reg8; reg94=reg91+reg94; reg101=reg101-reg97; reg91=reg72+reg74; T reg104=reg75+reg66;
    reg61=reg61-reg62; reg69=reg69-reg68; reg98=reg98-reg99; reg82=reg82-reg60; reg92=reg92-reg93;
    T reg105=reg8*reg44; reg65=reg65-reg64; reg63=reg76+reg63; T reg106=reg88+reg89; reg67=reg67+reg70;
    T reg107=reg8*reg48; reg87=reg86+reg87; reg71=reg8*reg71; reg78=ponderation*reg78; reg92=reg92*reg8;
    reg80=reg8*reg80; reg82=reg8*reg82; reg61=reg61*reg8; reg77=ponderation*reg77; reg86=reg8*reg104;
    reg98=reg98*reg8; reg100=reg8*reg100; reg94=reg8*reg94; reg69=reg8*reg69; T reg108=reg91*reg8;
    T reg109=ponderation*reg107; reg101=reg8*reg101; reg67=reg8*reg67; reg85=reg8*reg85; T reg110=reg8*reg106;
    T reg111=reg103*reg8; reg87=reg87*reg8; reg30=reg8*reg30; T reg112=ponderation*reg102; reg65=reg8*reg65;
    reg63=reg8*reg63; T reg113=ponderation*reg96; T reg114=ponderation*reg105; reg35=reg8*reg35; T reg115=ponderation*reg110;
    sollicitation[indices[1]+0]+=reg115; T tmp_3_3=ponderation*reg71; T tmp_3_4=ponderation*reg80; T tmp_5_5=ponderation*reg30; reg30=ponderation*reg86;
    sollicitation[indices[0]+1]+=reg30; T tmp_1_1=ponderation*reg94; T tmp_4_5=-reg112; T tmp_4_4=ponderation*reg63; T tmp_0_0=ponderation*reg100;
    T tmp_3_5=-reg114; T tmp_1_4=ponderation*reg92; T tmp_1_3=ponderation*reg98; T tmp_1_2=ponderation*reg101; T tmp_1_5=ponderation*reg69;
    T tmp_0_4=ponderation*reg82; T tmp_0_1=ponderation*reg87; T tmp_0_5=ponderation*reg65; T tmp_2_2=ponderation*reg67; T tmp_2_3=-reg109;
    T tmp_0_2=ponderation*reg61; T tmp_0_3=ponderation*reg85; reg61=ponderation*reg111; sollicitation[indices[2]+1]+=reg61; sollicitation[indices[2]+0]+=-reg78;
    T tmp_2_4=-reg113; T tmp_2_5=ponderation*reg35; sollicitation[indices[1]+1]+=-reg77; reg35=ponderation*reg108; sollicitation[indices[0]+0]+=reg35;
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
    T reg4=reg0*reg1; T reg5=reg4*reg2; T reg6=reg1*reg3; reg4=reg4*reg3; reg1=reg1*reg2;
    T reg7=reg5*reg3; reg5=reg5*reg2; T reg8=reg3*reg1; T reg9=reg4*reg3; T reg10=reg3*reg6;
    reg1=reg2*reg1; reg5=reg5-reg9; reg7=reg7+reg9; reg6=reg2*reg6; reg8=reg10+reg8;
    reg1=reg1-reg10; reg4=reg4*reg2; reg1=reg2*reg1; reg9=reg4+reg9; reg4=reg0*reg2;
    T reg11=reg3*reg7; T reg12=reg10+reg6; T reg13=reg0*reg3; reg8=reg3*reg8; T reg14=reg2*reg5;
    T reg15=pow(reg3,2); T reg16=elem.pos(1)[1]-elem.pos(0)[1]; reg11=reg14-reg11; reg14=reg2*reg4; T reg17=reg3*reg9;
    T reg18=reg3*reg13; T reg19=elem.pos(1)[0]-elem.pos(0)[0]; T reg20=elem.pos(2)[0]-elem.pos(0)[0]; T reg21=elem.pos(2)[1]-elem.pos(0)[1]; T reg22=pow(reg2,2);
    reg12=reg3*reg12; reg8=reg1-reg8; reg12=reg8-reg12; reg17=reg11-reg17; reg1=reg20*reg16;
    reg18=reg14-reg18; reg8=reg21*reg19; reg15=reg22-reg15; reg12=reg12/reg17; reg11=1-(*f.m).resolution;
    reg15=reg15/reg18; reg1=reg8-reg1; reg8=reg11*reg12; reg14=(*f.m).resolution*reg15; reg13=reg13/reg18;
    reg18=reg4/reg18; reg16=reg16/reg1; reg20=reg20/reg1; reg21=reg21/reg1; reg7=reg7/reg17;
    reg5=reg5/reg17; reg19=reg19/reg1; reg4=reg11*reg7; reg22=reg11*reg5; T reg23=reg20-reg19;
    T reg24=(*f.m).resolution*reg18; T reg25=0.5*reg20; T reg26=0.5*reg21; T reg27=reg16-reg21; reg8=reg14+reg8;
    reg14=(*f.m).resolution*reg13; T reg28=0.5*reg16; T reg29=0.5*reg19; T reg30=reg25*reg8; T reg31=0.5*reg23;
    T reg32=0.5*reg27; reg22=reg24+reg22; reg14=reg4+reg14; reg4=reg26*reg8; reg24=reg28*reg8;
    T reg33=reg29*reg8; reg4=2*reg4; T reg34=reg20*reg14; T reg35=reg31*reg8; reg33=2*reg33;
    T reg36=reg22*reg19; T reg37=reg16*reg14; T reg38=reg16*reg22; T reg39=reg22*reg21; T reg40=2*reg24;
    T reg41=reg14*reg19; T reg42=2*reg30; T reg43=reg20*reg22; T reg44=reg32*reg8; T reg45=reg14*reg21;
    T reg46=reg43*reg19; T reg47=reg27*reg14; T reg48=reg28*reg4; T reg49=reg37*reg19; T reg50=reg28*reg33;
    T reg51=reg26*reg40; reg44=2*reg44; T reg52=reg23*reg14; T reg53=reg23*reg22; T reg54=reg27*reg22;
    reg35=2*reg35; T reg55=reg20*reg36; T reg56=reg25*reg33; T reg57=reg29*reg40; T reg58=reg16*reg41;
    T reg59=reg20*reg45; T reg60=reg39*reg16; T reg61=reg38*reg21; T reg62=reg42*reg29; T reg63=reg42*reg26;
    T reg64=reg25*reg4; T reg65=reg34*reg21; reg55=reg51+reg55; reg56=reg61+reg56; T reg66=reg44*reg26;
    T reg67=reg20*reg53; T reg68=reg31*reg44; T reg69=reg43*reg20; T reg70=reg52*reg27; T reg71=reg26*reg4;
    T reg72=reg41*reg21; reg59=reg63+reg59; reg60=reg62+reg60; T reg73=reg27*reg54; T reg74=reg25*reg40;
    T reg75=reg31*reg35; T reg76=reg26*reg35; T reg77=reg20*reg47; T reg78=reg39*reg21; reg41=reg27*reg41;
    T reg79=reg42*reg25; T reg80=reg31*reg40; T reg81=reg27*reg38; T reg82=reg29*reg4; T reg83=reg37*reg20;
    reg64=reg65+reg64; T reg84=reg26*reg33; T reg85=reg44*reg25; T reg86=reg52*reg21; T reg87=reg31*reg4;
    T reg88=reg25*reg35; T reg89=reg54*reg21; T reg90=reg27*reg34; T reg91=reg31*reg33; T reg92=reg29*reg33;
    T reg93=reg32*reg40; T reg94=reg23*reg36; T reg95=reg16*reg34; T reg96=reg16*reg38; T reg97=reg29*reg35;
    reg54=reg16*reg54; T reg98=reg44*reg29; T reg99=reg28*reg40; reg36=reg36*reg19; reg52=reg52*reg16;
    reg50=reg49+reg50; reg58=reg57+reg58; T reg100=reg47*reg19; T reg101=reg28*reg35; reg48=reg46+reg48;
    reg39=reg39*reg27; T reg102=reg31*reg42; T reg103=reg53*reg19; T reg104=reg44*reg28; T reg105=reg42*reg28;
    T reg106=reg45*reg19; reg47=reg23*reg47; reg35=reg32*reg35; reg33=reg32*reg33; reg4=reg32*reg4;
    T reg107=reg32*reg42; T reg108=reg23*reg43; reg45=reg23*reg45; reg44=reg32*reg44; T reg109=reg23*reg37;
    reg53=reg23*reg53; reg71=reg71+reg69; reg101=reg100-reg101; reg100=reg1*reg64; reg39=reg39-reg102;
    reg84=reg84+reg83; reg78=reg78+reg79; reg4=reg4-reg108; reg41=reg41-reg80; reg104=reg103-reg104;
    reg33=reg33-reg109; reg87=reg87-reg90; reg103=reg1*reg58; reg92=reg96+reg92; T reg110=reg1*reg56;
    reg72=reg72+reg74; reg52=reg98-reg52; reg77=reg76-reg77; reg76=reg1*reg59; reg98=reg60*reg1;
    reg67=reg66-reg67; reg66=reg1*reg55; reg54=reg97-reg54; reg45=reg45-reg107; reg73=reg75+reg73;
    reg44=reg53+reg44; reg70=reg68+reg70; reg35=reg47+reg35; reg36=reg36+reg99; reg82=reg82+reg95;
    reg47=reg1*reg50; reg94=reg94-reg93; reg88=reg89-reg88; reg53=reg1*reg48; reg85=reg86-reg85;
    reg91=reg91-reg81; reg106=reg106+reg105; reg94=reg1*reg94; reg72=reg1*reg72; reg52=reg52*reg1;
    reg68=ponderation*reg47; reg75=ponderation*reg98; reg91=reg1*reg91; reg54=reg54*reg1; reg82=reg82*reg1;
    reg36=reg1*reg36; reg70=reg70*reg1; reg77=reg1*reg77; reg86=ponderation*reg66; reg35=reg1*reg35;
    reg67=reg1*reg67; reg89=ponderation*reg76; reg33=reg33*reg1; reg41=reg1*reg41; reg71=reg1*reg71;
    reg73=reg1*reg73; reg44=reg1*reg44; reg84=reg1*reg84; reg4=reg4*reg1; reg104=reg1*reg104;
    reg106=reg1*reg106; reg85=reg1*reg85; reg97=ponderation*reg100; reg39=reg39*reg1; reg101=reg1*reg101;
    reg88=reg1*reg88; reg87=reg1*reg87; T reg111=ponderation*reg103; reg78=reg1*reg78; T reg112=ponderation*reg53;
    reg45=reg1*reg45; reg92=reg1*reg92; T reg113=ponderation*reg110; T tmp_0_5=ponderation*reg41; T tmp_1_4=ponderation*reg33;
    T tmp_1_0=ponderation*reg35; T tmp_2_5=ponderation*reg72; T tmp_3_2=-reg89; T tmp_0_2=ponderation*reg39; T tmp_0_0=ponderation*reg73;
    T tmp_5_0=ponderation*reg101; T tmp_3_3=ponderation*reg71; T tmp_2_3=-reg97; T tmp_5_1=ponderation*reg104; T tmp_2_4=-reg113;
    T tmp_4_1=ponderation*reg52; T tmp_3_4=ponderation*reg84; T tmp_1_3=ponderation*reg4; T tmp_5_2=ponderation*reg106; T tmp_1_1=ponderation*reg44;
    T tmp_4_4=ponderation*reg92; T tmp_5_4=-reg68; T tmp_4_0=ponderation*reg54; T tmp_4_2=-reg75; T tmp_0_1=ponderation*reg70;
    T tmp_5_5=ponderation*reg36; T tmp_4_3=ponderation*reg82; T tmp_0_3=ponderation*reg87; T tmp_2_1=ponderation*reg85; T tmp_3_5=-reg86;
    T tmp_5_3=-reg112; T tmp_1_5=ponderation*reg94; T tmp_3_0=ponderation*reg77; T tmp_4_5=-reg111; T tmp_0_4=ponderation*reg91;
    T tmp_2_0=ponderation*reg88; T tmp_2_2=ponderation*reg78; T tmp_3_1=ponderation*reg67; T tmp_1_2=ponderation*reg45;
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
    T reg4=reg0*reg1; T reg5=reg1*reg3; T reg6=reg4*reg3; reg4=reg4*reg2; reg1=reg1*reg2;
    T reg7=reg3*reg1; T reg8=reg6*reg3; reg1=reg2*reg1; T reg9=reg3*reg5; T reg10=reg4*reg2;
    reg4=reg4*reg3; reg7=reg9+reg7; reg1=reg1-reg9; reg5=reg2*reg5; reg10=reg10-reg8;
    reg4=reg4+reg8; reg6=reg6*reg2; reg1=reg2*reg1; reg8=reg6+reg8; reg6=reg0*reg2;
    T reg11=reg9+reg5; T reg12=reg2*reg10; T reg13=reg3*reg4; T reg14=reg0*reg3; reg7=reg3*reg7;
    T reg15=elem.pos(2)[1]-elem.pos(0)[1]; T reg16=elem.pos(2)[0]-elem.pos(0)[0]; T reg17=pow(reg2,2); T reg18=elem.pos(1)[0]-elem.pos(0)[0]; T reg19=reg3*reg14;
    reg7=reg1-reg7; reg1=reg3*reg8; T reg20=reg2*reg6; reg13=reg12-reg13; reg11=reg3*reg11;
    reg12=elem.pos(1)[1]-elem.pos(0)[1]; T reg21=pow(reg3,2); reg11=reg7-reg11; reg1=reg13-reg1; reg7=reg16*reg12;
    reg19=reg20-reg19; reg13=reg15*reg18; reg21=reg17-reg21; reg17=1-(*f.m).resolution; reg21=reg21/reg19;
    reg7=reg13-reg7; reg11=reg11/reg1; reg6=reg6/reg19; reg19=reg14/reg19; reg4=reg4/reg1;
    reg15=reg15/reg7; reg10=reg10/reg1; reg18=reg18/reg7; reg13=reg17*reg11; reg12=reg12/reg7;
    reg16=reg16/reg7; reg14=(*f.m).resolution*reg21; reg20=reg16-reg18; T reg22=0.5*reg15; T reg23=0.5*reg18;
    T reg24=0.5*reg12; reg13=reg14+reg13; reg14=(*f.m).resolution*reg19; T reg25=reg17*reg4; T reg26=reg17*reg10;
    T reg27=(*f.m).resolution*reg6; T reg28=reg12-reg15; T reg29=0.5*reg28; reg14=reg25+reg14; reg26=reg27+reg26;
    reg25=reg22*reg13; reg27=reg23*reg13; T reg30=reg24*reg13; T reg31=0.5*reg16; T reg32=0.5*reg20;
    reg27=2*reg27; T reg33=reg31*reg13; T reg34=reg16*reg14; T reg35=2*reg30; reg25=2*reg25;
    T reg36=reg26*reg18; T reg37=reg32*reg13; T reg38=reg12*reg26; T reg39=reg14*reg18; T reg40=reg29*reg13;
    T reg41=2*reg33; reg37=2*reg37; T reg42=reg28*reg26; T reg43=reg26*reg15; T reg44=reg31*reg25;
    T reg45=reg34*reg15; T reg46=reg20*reg14; T reg47=reg16*reg36; T reg48=reg22*reg35; T reg49=reg12*reg14;
    T reg50=reg23*reg35; T reg51=reg12*reg39; reg40=2*reg40; T reg52=reg16*reg26; T reg53=reg31*reg27;
    T reg54=reg38*reg15; T reg55=reg14*reg15; T reg56=reg20*reg26; T reg57=reg23*reg27; reg47=reg48+reg47;
    T reg58=reg49*reg16; T reg59=reg28*reg38; T reg60=reg22*reg27; T reg61=reg52*reg16; T reg62=reg22*reg25;
    T reg63=reg32*reg25; T reg64=reg32*reg35; T reg65=reg28*reg34; T reg66=reg28*reg39; T reg67=reg12*reg38;
    T reg68=reg32*reg27; T reg69=reg32*reg40; T reg70=reg31*reg35; reg44=reg45+reg44; T reg71=reg43*reg15;
    T reg72=reg41*reg31; reg39=reg39*reg15; reg53=reg54+reg53; T reg73=reg36*reg18; T reg74=reg24*reg35;
    T reg75=reg20*reg52; T reg76=reg29*reg25; T reg77=reg20*reg56; T reg78=reg28*reg42; T reg79=reg43*reg28;
    T reg80=reg32*reg41; reg51=reg50+reg51; T reg81=reg32*reg37; T reg82=reg29*reg40; T reg83=reg20*reg49;
    T reg84=reg46*reg28; reg36=reg20*reg36; T reg85=reg29*reg35; T reg86=reg20*reg55; T reg87=reg29*reg27;
    T reg88=reg29*reg41; reg63=reg63-reg65; reg82=reg77+reg82; reg73=reg73+reg74; reg77=reg7*reg53;
    reg39=reg39+reg70; reg76=reg76-reg75; T reg89=reg7*reg51; reg57=reg67+reg57; reg62=reg62+reg61;
    T reg90=reg7*reg47; reg60=reg60+reg58; reg87=reg87-reg83; reg68=reg68-reg59; reg86=reg86-reg88;
    reg79=reg79-reg80; reg84=reg69+reg84; reg36=reg36-reg85; reg66=reg66-reg64; reg69=reg7*reg44;
    reg78=reg81+reg78; reg71=reg71+reg72; reg73=reg7*reg73; reg81=ponderation*reg90; reg39=reg7*reg39;
    reg57=reg7*reg57; reg66=reg7*reg66; reg76=reg76*reg7; reg62=reg7*reg62; T reg91=ponderation*reg89;
    reg68=reg7*reg68; reg60=reg7*reg60; T reg92=ponderation*reg77; reg78=reg7*reg78; reg86=reg7*reg86;
    reg71=reg7*reg71; reg84=reg84*reg7; reg36=reg7*reg36; reg63=reg7*reg63; T reg93=ponderation*reg69;
    reg82=reg7*reg82; reg87=reg87*reg7; reg79=reg79*reg7; T tmp_0_4=ponderation*reg68; T tmp_3_5=-reg81;
    T tmp_3_3=ponderation*reg62; T tmp_1_5=ponderation*reg36; T tmp_3_4=ponderation*reg60; T tmp_1_4=ponderation*reg87; T tmp_1_3=ponderation*reg76;
    T tmp_5_5=ponderation*reg73; T tmp_4_4=ponderation*reg57; T tmp_0_1=ponderation*reg84; T tmp_4_5=-reg91; T tmp_2_5=ponderation*reg39;
    T tmp_0_5=ponderation*reg66; T tmp_0_0=ponderation*reg78; T tmp_1_2=ponderation*reg86; T tmp_2_4=-reg92; T tmp_2_2=ponderation*reg71;
    T tmp_0_3=ponderation*reg63; T tmp_2_3=-reg93; T tmp_1_1=ponderation*reg82; T tmp_0_2=ponderation*reg79;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg2*reg3; reg2=reg2*reg4; T reg6=reg1*reg3; reg1=reg1*reg4;
    T reg7=reg5*reg4; reg5=reg5*reg3; T reg8=reg3*reg6; T reg9=reg4*reg1; reg6=reg4*reg6;
    T reg10=reg2*reg4; reg5=reg5-reg10; reg2=reg2*reg3; reg1=reg3*reg1; reg7=reg7+reg10;
    reg8=reg8-reg9; reg6=reg9+reg6; reg8=reg3*reg8; reg10=reg2+reg10; reg2=reg3*reg5;
    T reg11=reg9+reg1; reg6=reg4*reg6; T reg12=reg4*reg7; reg11=reg4*reg11; T reg13=reg4*reg10;
    reg6=reg8-reg6; reg12=reg2-reg12; reg11=reg6-reg11; reg13=reg12-reg13; reg11=reg11/reg13;
    reg7=reg7/reg13; reg5=reg5/reg13; reg2=reg11*reg5; reg6=reg11*reg7; reg8=(*f.m).alpha*reg7;
    reg12=reg2*reg5; T reg14=reg6*reg7; reg13=reg10/reg13; reg10=(*f.m).alpha*reg5; reg13=(*f.m).alpha*reg13;
    reg10=reg8+reg10; reg14=reg12-reg14; reg8=reg0*reg4; reg2=reg2/reg14; reg12=reg0*reg3;
    reg14=reg6/reg14; reg13=reg10+reg13; reg6=reg4*reg8; reg10=reg3*reg12; reg2=reg13*reg2;
    reg14=reg13*reg14; reg6=reg10-reg6; reg14=reg2-reg14; reg2=1-(*f.m).resolution; reg14=reg14*reg2;
    reg10=(*f.m).resolution*(*f.m).alpha; reg8=reg8/reg6; reg12=reg12/reg6; reg13=elem.pos(1)[0]-elem.pos(0)[0]; reg10=reg14+reg10;
    reg14=elem.pos(2)[1]-elem.pos(0)[1]; T reg15=elem.pos(2)[0]-elem.pos(0)[0]; reg7=reg2*reg7; reg5=reg2*reg5; T reg16=(*f.m).resolution*reg12;
    T reg17=elem.pos(1)[1]-elem.pos(0)[1]; T reg18=(*f.m).resolution*reg8; T reg19=reg14*reg13; reg18=reg7+reg18; reg5=reg16+reg5;
    reg7=reg15*reg17; reg10=(*f.m).deltaT*reg10; reg7=reg19-reg7; reg16=reg5*reg10; reg19=reg18*reg10;
    reg14=reg14/reg7; T reg20=1-var_inter[0]; reg15=reg15/reg7; reg17=reg17/reg7; reg13=reg13/reg7;
    T reg21=reg19+reg16; T reg22=reg21*reg15; T reg23=var_inter[1]*(*f.m).f_vol[0]; T reg24=var_inter[0]*(*f.m).f_vol[1]; reg20=reg20-var_inter[1];
    T reg25=reg15-reg13; T reg26=reg17-reg14; T reg27=reg21*reg17; T reg28=reg20*(*f.m).f_vol[1]; T reg29=reg21*reg14;
    T reg30=reg21*reg25; T reg31=reg20*(*f.m).f_vol[0]; T reg32=reg22-reg24; T reg33=reg27-reg23; T reg34=reg21*reg13;
    T reg35=var_inter[1]*(*f.m).f_vol[1]; T reg36=var_inter[0]*(*f.m).f_vol[0]; T reg37=reg21*reg26; T reg38=reg36+reg29; T reg39=reg28+reg30;
    T reg40=reg37+reg31; reg32=reg32*reg7; reg33=reg33*reg7; T reg41=reg34+reg35; reg32=ponderation*reg32;
    T reg42=reg40*reg7; T reg43=reg41*reg7; reg33=ponderation*reg33; T reg44=reg7*reg39; T reg45=reg7*reg38;
    sollicitation[indices[1]+1]+=-reg32; reg32=ponderation*reg43; sollicitation[indices[2]+1]+=reg32; T reg46=ponderation*reg42; sollicitation[indices[0]+0]+=reg46;
    T reg47=ponderation*reg44; sollicitation[indices[0]+1]+=reg47; sollicitation[indices[2]+0]+=-reg33; reg33=ponderation*reg45; sollicitation[indices[1]+0]+=reg33;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=reg0*reg1; T reg5=reg1*reg3; reg1=reg1*reg2; T reg6=reg4*reg2; reg4=reg4*reg3;
    T reg7=reg6*reg2; T reg8=reg2*reg1; T reg9=reg3*reg5; reg5=reg2*reg5; T reg10=reg4*reg3;
    reg4=reg4*reg2; reg6=reg6*reg3; reg10=reg10-reg7; reg1=reg3*reg1; reg5=reg8+reg5;
    reg9=reg9-reg8; reg4=reg4+reg7; reg5=reg2*reg5; reg7=reg6+reg7; reg9=reg3*reg9;
    reg6=reg8+reg1; T reg11=reg2*reg4; T reg12=reg3*reg10; reg11=reg12-reg11; reg12=reg2*reg7;
    reg6=reg2*reg6; reg5=reg9-reg5; reg12=reg11-reg12; reg6=reg5-reg6; reg6=reg6/reg12;
    reg4=reg4/reg12; reg10=reg10/reg12; reg5=reg6*reg4; reg9=reg6*reg10; reg11=reg5*reg4;
    T reg13=reg9*reg10; T reg14=(*f.m).alpha*reg10; T reg15=(*f.m).alpha*reg4; reg12=reg7/reg12; reg11=reg13-reg11;
    reg7=elem.pos(2)[0]-elem.pos(0)[0]; reg14=reg15+reg14; reg13=elem.pos(2)[1]-elem.pos(0)[1]; reg15=elem.pos(1)[0]-elem.pos(0)[0]; reg12=(*f.m).alpha*reg12;
    T reg16=elem.pos(1)[1]-elem.pos(0)[1]; T reg17=reg7*reg16; T reg18=reg13*reg15; reg9=reg9/reg11; reg11=reg5/reg11;
    reg12=reg14+reg12; reg11=reg12*reg11; reg9=reg12*reg9; reg5=reg0*reg3; reg17=reg18-reg17;
    reg12=reg0*reg2; reg14=vectors[0][indices[2]+0]-vectors[0][indices[0]+0]; reg13=reg13/reg17; reg18=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; reg15=reg15/reg17;
    T reg19=pow(reg3,2); T reg20=reg2*reg12; T reg21=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; T reg22=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; reg7=reg7/reg17;
    reg16=reg16/reg17; T reg23=reg3*reg5; T reg24=pow(reg2,2); reg11=reg9-reg11; reg9=1-(*f.m).resolution;
    T reg25=(*f.m).resolution*(*f.m).alpha; T reg26=reg22*reg13; T reg27=reg7*reg18; T reg28=reg14*reg15; T reg29=reg16*reg21;
    reg24=reg19-reg24; reg20=reg23-reg20; reg11=reg11*reg9; reg22=reg22*reg7; reg21=reg21*reg15;
    reg24=reg24/reg20; reg5=reg5/reg20; reg25=reg11+reg25; reg18=reg18*reg13; reg14=reg16*reg14;
    reg29=reg26-reg29; reg20=reg12/reg20; reg27=reg28-reg27; reg6=reg9*reg6; reg29=reg27+reg29;
    reg11=(*f.m).resolution*reg20; reg12=(*f.m).resolution*reg24; reg14=reg18-reg14; reg25=(*f.m).deltaT*reg25; reg4=reg9*reg4;
    reg10=reg9*reg10; reg9=(*f.m).resolution*reg5; reg22=reg21-reg22; reg29=0.5*reg29; reg14=reg14-reg25;
    reg6=reg12+reg6; reg22=reg22-reg25; reg11=reg4+reg11; reg10=reg9+reg10; reg4=reg22*reg10;
    reg9=reg11*reg14; reg14=reg10*reg14; reg12=reg16-reg13; reg18=reg7-reg15; reg22=reg22*reg11;
    reg29=reg6*reg29; reg29=2*reg29; reg14=reg22+reg14; reg4=reg9+reg4; reg9=0.5*reg7;
    reg19=0.5*reg13; reg21=0.5*reg12; reg22=0.5*reg15; reg23=0.5*reg18; reg26=0.5*reg16;
    reg27=1-var_inter[0]; reg28=reg9*reg29; T reg30=reg14*reg13; T reg31=reg7*reg4; T reg32=reg19*reg29;
    T reg33=reg12*reg14; T reg34=reg21*reg29; reg27=reg27-var_inter[1]; T reg35=reg16*reg14; T reg36=reg22*reg29;
    T reg37=reg18*reg4; T reg38=reg23*reg29; T reg39=reg4*reg15; T reg40=reg26*reg29; reg34=reg37+reg34;
    reg37=reg27*(*f.m).f_vol[1]; T reg41=reg27*(*f.m).f_vol[0]; reg33=reg38+reg33; reg30=reg30-reg28; reg39=reg39-reg40;
    reg38=var_inter[1]*(*f.m).f_vol[1]; T reg42=var_inter[1]*(*f.m).f_vol[0]; reg36=reg36-reg35; T reg43=var_inter[0]*(*f.m).f_vol[1]; reg32=reg32-reg31;
    T reg44=var_inter[0]*(*f.m).f_vol[0]; reg34=reg34-reg37; reg32=reg32-reg43; reg39=reg39-reg38; reg33=reg33-reg41;
    reg36=reg36-reg42; reg30=reg30-reg44; reg34=reg34*reg17; reg32=reg17*reg32; reg39=reg17*reg39;
    reg33=reg33*reg17; reg30=reg17*reg30; reg36=reg17*reg36; sollicitation[indices[2]+1]+=ponderation*reg39; sollicitation[indices[2]+0]+=ponderation*reg36;
    sollicitation[indices[0]+0]+=ponderation*reg33; sollicitation[indices[1]+1]+=ponderation*reg32; sollicitation[indices[0]+1]+=ponderation*reg34; sollicitation[indices[1]+0]+=ponderation*reg30;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=(*f.m).poisson_ratio/(*f.m).elastic_modulus;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg4*reg2; reg2=reg3*reg2; T reg6=reg3*reg1; reg1=reg4*reg1;
    T reg7=reg4*reg5; reg5=reg3*reg5; T reg8=reg1*reg4; T reg9=reg6*reg3; T reg10=reg3*reg2;
    reg1=reg1*reg3; reg1=reg9+reg1; reg8=reg8-reg9; reg6=reg6*reg4; reg2=reg4*reg2;
    reg5=reg10+reg5; reg7=reg7-reg10; reg1=reg1*reg3; T reg11=0.5*elem.pos(1)[0]; T reg12=reg9+reg6;
    T reg13=reg7*reg4; reg8=reg8*reg4; T reg14=reg5*reg3; reg2=reg10+reg2; reg10=0.5*elem.pos(0)[1];
    T reg15=0.5*elem.pos(1)[1]; T reg16=0.5*elem.pos(0)[0]; T reg17=reg15-reg10; reg12=reg12*reg3; reg1=reg8-reg1;
    reg8=reg11+reg16; reg14=reg13-reg14; reg16=reg11-reg16; reg11=0.5*elem.pos(2)[0]; reg13=reg2*reg3;
    reg15=reg10+reg15; reg10=0.5*elem.pos(2)[1]; reg13=reg14-reg13; reg8=reg11-reg8; reg17=reg10+reg17;
    reg12=reg1-reg12; reg15=reg10-reg15; reg1=0.5*elem.pos(3)[1]; reg10=0.5*elem.pos(3)[0]; reg16=reg11+reg16;
    reg15=reg15+reg1; reg5=reg5/reg13; reg7=reg7/reg13; reg11=0.5*vectors[0][indices[1]+1]; reg14=0.5*vectors[0][indices[0]+1];
    T reg18=0.5*vectors[0][indices[0]+0]; T reg19=0.5*vectors[0][indices[1]+0]; reg8=reg10+reg8; reg10=reg16-reg10; reg12=reg12/reg13;
    reg1=reg17-reg1; reg16=reg7*reg12; reg17=reg5*reg12; T reg20=0.21132486540518713447*elem.pos(0)[1]; T reg21=0.78867513459481286553*elem.pos(1)[1];
    T reg22=0.21132486540518713447*elem.pos(1)[0]; T reg23=0.21132486540518713447*elem.pos(0)[0]; T reg24=0.21132486540518713447*elem.pos(1)[1]; T reg25=0.78867513459481286553*elem.pos(1)[0]; T reg26=reg15*reg10;
    T reg27=reg1*reg8; T reg28=reg19-reg18; T reg29=0.5*vectors[0][indices[2]+0]; T reg30=0.78867513459481286553*elem.pos(0)[0]; reg18=reg19+reg18;
    reg19=0.78867513459481286553*elem.pos(0)[1]; T reg31=reg11-reg14; T reg32=0.5*vectors[0][indices[2]+1]; reg11=reg14+reg11; reg28=reg28+reg29;
    reg14=0.5*vectors[0][indices[3]+0]; reg18=reg29-reg18; reg13=reg2/reg13; reg11=reg32-reg11; reg2=0.5*vectors[0][indices[3]+1];
    reg29=reg23+reg25; reg31=reg32+reg31; reg32=reg24-reg20; reg23=reg22-reg23; T reg33=0.78867513459481286553*elem.pos(2)[1];
    T reg34=0.78867513459481286553*elem.pos(2)[0]; reg20=reg20+reg21; T reg35=reg7*reg16; reg24=reg19+reg24; T reg36=reg7*(*f.m).alpha;
    T reg37=reg5*reg17; T reg38=reg5*(*f.m).alpha; T reg39=0.21132486540518713447*elem.pos(2)[1]; reg22=reg30+reg22; reg27=reg26-reg27;
    reg26=0.21132486540518713447*elem.pos(2)[0]; reg10=reg10/reg27; reg8=reg8/reg27; T reg40=0.21132486540518713447*elem.pos(3)[0]; reg29=reg34-reg29;
    T reg41=0.78867513459481286553*elem.pos(3)[1]; reg32=reg33+reg32; reg15=reg15/reg27; reg31=reg31-reg2; reg19=reg21-reg19;
    reg21=0.78867513459481286553*elem.pos(3)[0]; reg23=reg34+reg23; reg30=reg25-reg30; reg25=0.21132486540518713447*elem.pos(3)[1]; reg20=reg33-reg20;
    reg33=reg4*reg0; reg34=reg3*reg0; reg22=reg26-reg22; reg28=reg28-reg14; reg27=reg1/reg27;
    reg13=reg13*(*f.m).alpha; reg38=reg36+reg38; reg2=reg11+reg2; reg24=reg39-reg24; reg37=reg35-reg37;
    reg18=reg14+reg18; reg1=0.21132486540518713447*PNODE(1).dep[1]; reg11=0.78867513459481286553*PNODE(1).dep[1]; reg14=0.21132486540518713447*PNODE(0).dep[1]; reg20=reg20+reg25;
    reg35=0.78867513459481286553*PNODE(1).dep[0]; reg16=reg16/reg37; reg36=reg8*reg31; T reg42=reg27*reg18; T reg43=reg2*reg10;
    reg29=reg29+reg40; T reg44=reg33*reg4; reg23=reg23-reg21; reg37=reg17/reg37; reg13=reg38+reg13;
    reg17=0.21132486540518713447*PNODE(0).dep[0]; reg38=reg28*reg15; T reg45=0.21132486540518713447*PNODE(1).dep[0]; reg32=reg32-reg41; T reg46=0.78867513459481286553*PNODE(0).dep[0];
    reg41=reg24+reg41; reg21=reg22+reg21; reg19=reg39+reg19; reg22=reg34*reg3; reg30=reg26+reg30;
    reg24=0.78867513459481286553*PNODE(0).dep[1]; reg26=0.78867513459481286553*PNODE(2).dep[1]; reg16=reg16*reg13; reg13=reg37*reg13; reg37=0.21132486540518713447*PNODE(2).dep[1];
    reg39=reg21*reg32; T reg47=reg14+reg11; T reg48=reg41*reg23; T reg49=reg45-reg17; reg45=reg46+reg45;
    reg42=reg38-reg42; elem.epsilon[0][0]=reg42; reg38=reg24+reg1; reg14=reg1-reg14; reg36=reg43-reg36;
    elem.epsilon[0][1]=reg36; reg1=(*f.m).deltaT*(*f.m).alpha; reg40=reg30-reg40; reg30=reg20*reg23; reg25=reg19-reg25;
    reg19=0.78867513459481286553*PNODE(2).dep[0]; reg22=reg44-reg22; reg17=reg17+reg35; reg43=0.21132486540518713447*PNODE(2).dep[0]; reg44=reg32*reg29;
    reg33=reg33/reg22; reg45=reg43-reg45; reg49=reg49+reg19; reg17=reg19-reg17; reg19=0.78867513459481286553*PNODE(3).dep[0];
    T reg50=0.21132486540518713447*PNODE(3).dep[0]; T reg51=0.78867513459481286553*PNODE(3).dep[1]; reg14=reg26+reg14; T reg52=reg20*reg40; T reg53=reg29*reg25;
    reg46=reg35-reg46; reg35=1-(*f.m).resolution; reg44=reg30-reg44; reg30=0.21132486540518713447*PNODE(3).dep[1]; reg24=reg11-reg24;
    reg47=reg26-reg47; reg11=reg42-reg1; reg26=reg36-reg1; reg39=reg48-reg39; reg13=reg16-reg13;
    reg38=reg37-reg38; reg34=reg34/reg22; reg47=reg47+reg30; reg16=reg23/reg44; reg48=reg29/reg44;
    T reg54=reg41/reg39; reg17=reg50+reg17; T reg55=(*f.m).resolution*(*f.m).alpha; T reg56=reg32/reg44; T reg57=pow(reg3,2);
    reg53=reg52-reg53; reg46=reg43+reg46; reg43=reg26*reg33; reg52=reg20/reg44; reg24=reg37+reg24;
    reg37=reg11*reg33; T reg58=reg21/reg39; T reg59=reg21*reg25; T reg60=reg41*reg40; reg23=reg23/reg39;
    reg38=reg38+reg51; reg26=reg26*reg34; T reg61=pow(reg4,2); reg11=reg11*reg34; reg45=reg45+reg19;
    reg32=reg32/reg39; reg51=reg14-reg51; reg19=reg49-reg19; reg13=reg35*reg13; reg43=reg11+reg43;
    reg57=reg61-reg57; reg11=reg19*reg48; reg14=reg17*reg16; reg49=reg54*reg19; reg61=reg32*reg45;
    reg13=reg55+reg13; reg54=reg54*reg51; reg55=reg38*reg23; T reg62=reg58*reg51; reg30=reg24-reg30;
    reg29=reg29/reg53; reg24=reg40/reg53; T reg63=reg25/reg53; reg50=reg46-reg50; reg20=reg20/reg53;
    reg58=reg58*reg19; reg23=reg45*reg23; reg32=reg32*reg38; reg46=reg56*reg47; reg56=reg56*reg17;
    reg19=reg52*reg19; reg52=reg52*reg51; reg16=reg47*reg16; reg51=reg48*reg51; reg59=reg60-reg59;
    reg26=reg37+reg26; reg56=reg19-reg56; reg13=(*f.m).deltaT*reg13; reg62=reg55-reg62; reg19=reg26*reg4;
    reg37=reg29*reg30; reg48=reg20*reg30; reg55=reg47*reg63; reg32=reg54-reg32; reg54=reg17*reg24;
    reg29=reg50*reg29; reg58=reg23-reg58; reg61=reg49-reg61; reg41=reg41/reg59; reg25=reg25/reg59;
    reg24=reg47*reg24; reg21=reg21/reg59; reg20=reg20*reg50; reg40=reg40/reg59; reg63=reg17*reg63;
    reg17=(*f.m).resolution*reg33; reg46=reg52-reg46; reg11=reg14-reg11; reg26=reg26*reg3; reg7=reg7*reg35;
    reg22=reg57/reg22; reg51=reg16-reg51; reg14=reg4*reg43; reg43=reg3*reg43; reg5=reg5*reg35;
    reg16=(*f.m).resolution*reg34; reg14=reg14-reg26; reg61=reg61-reg13; reg62=reg62-reg13; reg23=reg41*reg50;
    reg50=reg21*reg50; reg47=reg45*reg40; reg49=reg38*reg25; reg41=reg41*reg30; reg25=reg45*reg25;
    reg30=reg21*reg30; reg40=reg38*reg40; reg11=reg46+reg11; reg63=reg20-reg63; reg56=reg56-reg13;
    reg17=reg7+reg17; reg5=reg16+reg5; reg51=reg51-reg13; reg35=reg12*reg35; reg37=reg24-reg37;
    reg19=reg19-reg43; reg29=reg54-reg29; reg58=reg32+reg58; reg7=(*f.m).resolution*reg22; reg55=reg48-reg55;
    reg37=reg37-reg13; reg49=reg41-reg49; reg19=reg1+reg19; reg12=reg61*reg17; reg26=reg43+reg26;
    reg8=reg28*reg8; reg25=reg23-reg25; reg11=0.5*reg11; reg14=reg1+reg14; reg10=reg18*reg10;
    reg15=reg31*reg15; reg2=reg27*reg2; reg50=reg47-reg50; reg29=reg55+reg29; reg58=0.5*reg58;
    reg16=reg62*reg17; reg18=reg61*reg5; reg35=reg7+reg35; reg7=reg56*reg17; reg20=reg51*reg5;
    reg21=reg56*reg5; reg23=reg51*reg17; reg24=reg62*reg5; reg63=reg63-reg13; reg30=reg40-reg30;
    reg29=0.5*reg29; reg27=reg35*reg58; reg28=reg37*reg17; reg31=reg5*reg63; reg32=reg37*reg5;
    reg16=reg18+reg16; reg18=reg17*reg63; reg8=reg10-reg8; reg24=reg12+reg24; reg2=reg15-reg2;
    reg26=reg1-reg26; reg25=reg25-reg13; reg30=reg30-reg13; reg20=reg7+reg20; reg7=reg19+reg14;
    reg50=reg49+reg50; reg23=reg21+reg23; reg10=reg11*reg35; reg23=reg51*reg23; reg28=reg31+reg28;
    reg7=reg7+reg26; reg18=reg32+reg18; reg16=reg62*reg16; reg20=reg56*reg20; reg24=reg61*reg24;
    reg8=reg2+reg8; reg10=2*reg10; reg27=2*reg27; reg2=reg25*reg17; reg12=reg30*reg5;
    reg50=0.5*reg50; reg15=reg30*reg17; reg21=reg25*reg5; reg31=reg29*reg35; reg10=reg11*reg10;
    reg8=0.5*reg8; elem.epsilon[0][2]=reg8; reg11=reg35*reg50; reg15=reg21+reg15; reg12=reg2+reg12;
    reg63=reg18*reg63; reg7=reg7/3; reg31=2*reg31; reg28=reg37*reg28; reg16=reg24+reg16;
    reg23=reg20+reg23; reg27=reg58*reg27; reg15=reg30*reg15; reg12=reg25*reg12; reg19=reg19-reg7;
    reg14=reg14-reg7; reg2=reg8*reg22; reg11=2*reg11; reg23=reg10+reg23; reg27=reg16+reg27;
    reg31=reg29*reg31; reg28=reg63+reg28; reg15=reg12+reg15; reg27=reg39*reg27; reg31=reg28+reg31;
    reg19=pow(reg19,2); reg14=pow(reg14,2); reg11=reg50*reg11; reg7=reg26-reg7; reg2=reg0*reg2;
    reg44=reg23*reg44; reg14=reg19+reg14; reg7=pow(reg7,2); reg10=2*reg2; reg53=reg31*reg53;
    reg44=0.25*reg44; reg27=0.25*reg27; reg11=reg15+reg11; reg53=0.25*reg53; reg27=reg44+reg27;
    reg10=reg2*reg10; reg11=reg59*reg11; reg7=reg14+reg7; reg11=0.25*reg11; reg10=reg7+reg10;
    reg42=reg42-reg13; reg36=reg36-reg13; reg27=reg53+reg27; reg11=reg27+reg11; reg2=reg42*reg17;
    reg42=reg42*reg5; reg10=1.5*reg10; reg7=reg36*reg5; reg36=reg36*reg17; elem.sigma[0][2]=reg8*reg35;
    elem.sigma_von_mises=pow(reg10,0.5); elem.sigma[0][1]=reg42+reg36; elem.sigma[0][0]=reg2+reg7; elem.ener=reg11/2;
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
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg4*reg2; T reg6=reg3*reg1; reg2=reg3*reg2; reg1=reg4*reg1;
    T reg7=reg4*reg5; reg5=reg3*reg5; T reg8=reg3*reg2; T reg9=reg1*reg4; T reg10=reg6*reg3;
    reg1=reg1*reg3; reg7=reg7-reg8; reg5=reg8+reg5; reg2=reg4*reg2; reg1=reg10+reg1;
    reg6=reg6*reg4; reg9=reg9-reg10; T reg11=reg10+reg6; reg2=reg8+reg2; reg1=reg1*reg3;
    reg9=reg9*reg4; reg8=reg5*reg3; T reg12=reg7*reg4; reg1=reg9-reg1; reg11=reg11*reg3;
    reg8=reg12-reg8; reg9=reg2*reg3; reg9=reg8-reg9; reg11=reg1-reg11; reg1=1-var_inter[1];
    reg8=1-var_inter[0]; reg12=reg1*elem.pos(1)[0]; T reg13=reg1*elem.pos(0)[0]; T reg14=reg1*elem.pos(0)[1]; T reg15=reg1*elem.pos(1)[1];
    reg7=reg7/reg9; T reg16=elem.pos(0)[0]*reg8; T reg17=elem.pos(1)[0]*var_inter[0]; T reg18=reg8*elem.pos(0)[1]; T reg19=elem.pos(1)[1]*var_inter[0];
    reg11=reg11/reg9; reg5=reg5/reg9; T reg20=reg7*reg11; T reg21=elem.pos(2)[0]*var_inter[1]; T reg22=reg5*reg11;
    reg12=reg12-reg13; T reg23=elem.pos(2)[1]*var_inter[1]; reg15=reg15-reg14; T reg24=elem.pos(2)[0]*var_inter[0]; T reg25=reg16+reg17;
    T reg26=reg18+reg19; T reg27=elem.pos(2)[1]*var_inter[0]; T reg28=elem.pos(3)[1]*reg8; T reg29=reg5*reg22; T reg30=reg7*reg20;
    reg27=reg27-reg26; T reg31=reg7*(*f.m).alpha; T reg32=reg5*(*f.m).alpha; T reg33=elem.pos(3)[0]*reg8; reg24=reg24-reg25;
    T reg34=elem.pos(3)[1]*var_inter[1]; reg9=reg2/reg9; reg15=reg23+reg15; reg2=elem.pos(3)[0]*var_inter[1]; reg12=reg21+reg12;
    reg29=reg30-reg29; reg28=reg27+reg28; reg33=reg24+reg33; reg15=reg15-reg34; reg12=reg12-reg2;
    reg9=reg9*(*f.m).alpha; reg32=reg31+reg32; reg21=reg12*reg28; reg23=reg15*reg33; reg9=reg32+reg9;
    reg20=reg20/reg29; reg29=reg22/reg29; reg22=reg4*reg0; reg24=reg3*reg0; reg27=pow(reg3,2);
    reg30=reg24*reg3; reg31=pow(reg4,2); reg23=reg21-reg23; reg29=reg29*reg9; reg9=reg20*reg9;
    reg20=reg22*reg4; reg21=1-(*f.m).resolution; reg33=reg33/reg23; reg30=reg20-reg30; reg29=reg9-reg29;
    reg28=reg28/reg23; reg15=reg15/reg23; reg12=reg12/reg23; reg27=reg31-reg27; reg9=reg1*reg28;
    reg24=reg24/reg30; reg27=reg27/reg30; reg29=reg21*reg29; reg20=(*f.m).resolution*(*f.m).alpha; reg31=var_inter[0]*reg15;
    reg30=reg22/reg30; reg22=var_inter[1]*reg33; reg32=var_inter[1]*reg28; T reg35=var_inter[0]*reg12; T reg36=reg1*reg33;
    T reg37=reg8*reg12; T reg38=reg8*reg15; T reg39=reg38+reg32; T reg40=reg37+reg22; T reg41=reg36+reg35;
    T reg42=reg31+reg9; reg29=reg20+reg29; reg20=(*f.m).resolution*reg30; T reg43=(*f.m).resolution*reg27; T reg44=(*f.m).resolution*reg24;
    reg5=reg5*reg21; reg11=reg11*reg21; reg21=reg7*reg21; reg7=reg35-reg22; T reg45=reg32-reg31;
    T reg46=0.5*reg40; T reg47=0.5*reg42; T reg48=0.5*reg41; reg20=reg21+reg20; reg21=0.5*reg39;
    reg11=reg43+reg11; reg29=(*f.m).deltaT*reg29; reg43=reg38-reg9; reg5=reg44+reg5; reg44=reg36-reg37;
    T reg49=0.5*reg45; T reg50=0.5*reg7; T reg51=reg21*reg11; T reg52=reg46*reg11; T reg53=0.5*reg44;
    T reg54=reg5*reg29; T reg55=0.5*reg43; T reg56=reg20*reg29; T reg57=reg11*reg48; T reg58=reg11*reg47;
    T reg59=reg40*reg5; T reg60=reg54+reg56; T reg61=reg5*reg42; T reg62=reg11*reg50; T reg63=reg11*reg49;
    reg58=2*reg58; T reg64=reg20*reg41; T reg65=reg5*reg41; reg52=2*reg52; T reg66=2*reg51;
    T reg67=2*reg57; T reg68=reg5*reg39; T reg69=reg20*reg39; T reg70=reg40*reg20; T reg71=reg53*reg11;
    T reg72=reg11*reg55; T reg73=reg20*reg42; T reg74=reg1*var_inter[0]; T reg75=reg8*var_inter[1]; T reg76=reg52*reg48;
    T reg77=reg74*(*f.m).f_vol[1]; T reg78=reg5*reg45; reg63=2*reg63; T reg79=reg73*reg39; T reg80=reg46*reg67;
    T reg81=reg20*reg7; reg72=2*reg72; T reg82=reg21*reg52; T reg83=reg40*reg68; T reg84=reg75*(*f.m).f_vol[0];
    T reg85=reg44*reg5; T reg86=reg40*reg64; T reg87=reg43*reg5; T reg88=reg70*reg41; T reg89=reg66*reg47;
    T reg90=reg65*reg42; T reg91=reg21*reg58; T reg92=reg44*reg20; T reg93=reg61*reg41; reg71=2*reg71;
    T reg94=reg67*reg47; T reg95=reg60*reg39; T reg96=reg43*reg20; T reg97=reg58*reg48; T reg98=reg60*reg41;
    reg62=2*reg62; T reg99=reg5*reg7; T reg100=reg59*reg39; T reg101=var_inter[0]*var_inter[1]; T reg102=reg1*reg8;
    T reg103=reg46*reg66; T reg104=reg20*reg45; T reg105=reg69*reg42; T reg106=reg102*(*f.m).f_vol[0]; T reg107=reg104*reg43;
    T reg108=reg78*reg41; T reg109=reg81*reg41; T reg110=reg81*reg7; T reg111=reg40*reg81; T reg112=reg63*reg49;
    T reg113=reg62*reg47; T reg114=reg21*reg63; T reg115=reg62*reg49; T reg116=reg78*reg7; T reg117=reg62*reg53;
    T reg118=reg68*reg7; T reg119=reg40*reg61; T reg120=reg52*reg49; T reg121=reg63*reg47; T reg122=reg43*reg60;
    T reg123=reg63*reg48; T reg124=reg70*reg7; T reg125=reg66*reg49; reg82=reg83+reg82; T reg126=reg21*reg66;
    T reg127=reg40*reg70; T reg128=reg44*reg60; T reg129=reg46*reg71; T reg130=reg65*reg39; T reg131=reg60*reg45;
    T reg132=reg63*reg53; T reg133=reg46*reg62; T reg134=reg104*reg39; T reg135=reg69*reg43; T reg136=reg40*reg87;
    T reg137=reg21*reg71; T reg138=reg98-reg77; T reg139=reg40*reg92; T reg140=reg21*reg72; T reg141=reg46*reg52;
    T reg142=reg99*reg39; T reg143=reg52*reg53; T reg144=reg46*reg63; T reg145=reg66*reg53; T reg146=reg60*reg42;
    T reg147=reg59*reg43; reg93=reg94+reg93; T reg148=reg58*reg47; T reg149=reg21*reg62; T reg150=reg40*reg78;
    T reg151=reg72*reg55; T reg152=reg92*reg44; reg91=reg86+reg91; T reg153=reg71*reg55; T reg154=reg87*reg44;
    T reg155=reg96*reg39; T reg156=reg21*reg67; T reg157=reg99*reg43; T reg158=reg46*reg72; T reg159=reg85*reg39;
    T reg160=reg64*reg41; T reg161=reg60*reg7; reg79=reg80+reg79; T reg162=reg46*reg58; reg76=reg105+reg76;
    reg78=reg44*reg78; T reg163=reg62*reg55; T reg164=reg71*reg50; T reg165=reg74*(*f.m).f_vol[0]; T reg166=reg96*reg45;
    reg81=reg44*reg81; T reg167=reg63*reg55; T reg168=reg59*reg42; T reg169=reg66*reg48; T reg170=reg44*reg68;
    T reg171=reg52*reg55; reg88=reg89+reg88; T reg172=reg102*(*f.m).f_vol[1]; reg70=reg44*reg70; T reg173=reg66*reg55;
    T reg174=reg75*(*f.m).f_vol[1]; T reg175=reg40*reg60; T reg176=reg96*reg42; T reg177=reg71*reg48; T reg178=reg95-reg84;
    T reg179=reg99*reg42; T reg180=reg62*reg48; T reg181=reg65*reg45; T reg182=reg58*reg50; T reg183=reg104*reg42;
    reg97=reg90+reg97; reg104=reg104*reg45; reg62=reg62*reg50; T reg184=reg67*reg50; T reg185=reg73*reg45;
    reg99=reg99*reg45; reg63=reg63*reg50; T reg186=reg101*(*f.m).f_vol[1]; T reg187=reg101*(*f.m).f_vol[0]; reg100=reg100+reg103;
    T reg188=reg44*reg61; T reg189=reg67*reg55; T reg190=reg72*reg50; T reg191=reg85*reg45; T reg192=reg44*reg64;
    T reg193=reg58*reg55; T reg194=reg71*reg47; reg59=reg59*reg45; T reg195=reg66*reg50; T reg196=reg72*reg53;
    T reg197=reg87*reg7; T reg198=reg71*reg49; reg87=reg87*reg41; T reg199=reg73*reg43; T reg200=reg68*reg41;
    T reg201=reg92*reg7; T reg202=reg72*reg49; T reg203=reg52*reg47; T reg204=reg67*reg53; reg61=reg61*reg7;
    T reg205=reg67*reg49; T reg206=reg65*reg43; T reg207=reg72*reg47; T reg208=reg64*reg7; T reg209=reg58*reg49;
    reg92=reg92*reg41; reg58=reg58*reg53; reg96=reg96*reg43; reg73=reg73*reg42; T reg210=reg67*reg48;
    reg72=reg72*reg48; reg71=reg71*reg53; T reg211=reg69*reg45; T reg212=reg85*reg42; reg52=reg52*reg50;
    reg85=reg85*reg43; T reg213=reg69*reg39; reg149=reg150-reg149; reg143=reg143-reg135; reg178=reg178*reg23;
    reg141=reg213+reg141; reg150=reg97*reg23; reg58=reg58-reg206; T reg214=reg93*reg23; reg87=reg194-reg87;
    reg180=reg183-reg180; reg142=reg144-reg142; reg144=reg91*reg23; reg199=reg199-reg204; reg117=reg107+reg117;
    reg107=reg174+reg175; reg168=reg168+reg169; reg132=reg157+reg132; reg157=reg100*reg23; reg123=reg179-reg123;
    reg71=reg96+reg71; reg119=reg156+reg119; reg96=reg82*reg23; reg114=reg111-reg114; reg196=reg85+reg196;
    reg92=reg207-reg92; reg85=reg76*reg23; reg148=reg148+reg160; reg203=reg203+reg200; reg59=reg59-reg195;
    reg198=reg197+reg198; reg202=reg201+reg202; reg109=reg121-reg109; reg61=reg61-reg205; reg209=reg209-reg208;
    reg115=reg116+reg115; reg112=reg110+reg112; reg108=reg113-reg108; reg120=reg120-reg118; reg124=reg124-reg125;
    reg110=reg106+reg122; reg111=reg161+reg186; reg127=reg127+reg126; reg151=reg152+reg151; reg113=reg131+reg187;
    reg185=reg185-reg184; reg182=reg182-reg181; reg190=reg191+reg190; reg62=reg104+reg62; reg63=reg99+reg63;
    reg188=reg188-reg189; reg193=reg193-reg192; reg163=reg78+reg163; reg164=reg166+reg164; reg167=reg81+reg167;
    reg171=reg171-reg170; reg70=reg70-reg173; reg78=reg88*reg23; reg177=reg176-reg177; reg72=reg212-reg72;
    reg73=reg73+reg210; reg52=reg52-reg211; reg81=reg128+reg172; reg147=reg147-reg145; reg140=reg139-reg140;
    reg137=reg136-reg137; reg134=reg133-reg134; reg162=reg162+reg130; reg99=reg146+reg165; reg104=reg79*reg23;
    reg159=reg158-reg159; reg155=reg129-reg155; reg138=reg138*reg23; reg153=reg154+reg153; reg141=reg23*reg141;
    reg171=reg171*reg23; reg116=ponderation*reg85; reg121=reg107*reg23; reg70=reg70*reg23; reg129=ponderation*reg78;
    reg178=ponderation*reg178; reg134=reg134*reg23; reg177=reg177*reg23; reg132=reg132*reg23; reg72=reg72*reg23;
    reg168=reg168*reg23; reg133=reg113*reg23; reg73=reg73*reg23; reg71=reg71*reg23; reg203=reg203*reg23;
    reg142=reg142*reg23; reg52=reg52*reg23; reg162=reg162*reg23; reg136=ponderation*reg144; reg180=reg180*reg23;
    reg185=reg185*reg23; reg143=reg143*reg23; reg190=reg190*reg23; reg139=ponderation*reg150; reg182=reg182*reg23;
    reg147=reg147*reg23; reg62=reg62*reg23; reg123=reg123*reg23; reg152=reg81*reg23; reg63=reg63*reg23;
    reg140=reg140*reg23; reg188=reg188*reg23; reg153=reg153*reg23; reg193=reg193*reg23; reg164=reg164*reg23;
    reg163=reg163*reg23; reg137=reg137*reg23; reg151=reg151*reg23; reg167=reg167*reg23; reg92=reg92*reg23;
    reg58=reg58*reg23; reg127=reg127*reg23; reg115=reg115*reg23; reg108=reg108*reg23; reg112=reg112*reg23;
    reg155=reg155*reg23; reg117=reg117*reg23; reg138=ponderation*reg138; reg120=reg120*reg23; reg154=ponderation*reg157;
    reg148=reg148*reg23; reg119=reg119*reg23; reg124=reg124*reg23; reg158=ponderation*reg214; reg166=ponderation*reg96;
    reg176=reg110*reg23; reg114=reg114*reg23; reg179=reg111*reg23; reg149=reg149*reg23; reg59=reg59*reg23;
    reg196=reg196*reg23; reg87=reg87*reg23; reg198=reg198*reg23; reg183=ponderation*reg104; reg202=reg202*reg23;
    reg109=reg109*reg23; reg199=reg199*reg23; reg209=reg209*reg23; reg159=reg159*reg23; reg61=reg61*reg23;
    reg191=reg99*reg23; T tmp_2_6=-reg116; T tmp_3_5=ponderation*reg109; T tmp_4_1=ponderation*reg190; reg109=ponderation*reg179;
    sollicitation[indices[2]+1]+=reg109; T tmp_3_2=-reg158; T tmp_6_5=ponderation*reg142; reg116=ponderation*reg191; sollicitation[indices[1]+0]+=reg116;
    sollicitation[indices[1]+1]+=-reg138; reg138=ponderation*reg133; sollicitation[indices[2]+0]+=reg138; T tmp_2_7=ponderation*reg168; T tmp_3_4=ponderation*reg108;
    T tmp_3_0=ponderation*reg87; T tmp_3_1=ponderation*reg92; T tmp_2_5=ponderation*reg123; T tmp_3_7=-reg129; T tmp_3_3=ponderation*reg148;
    T tmp_4_0=ponderation*reg164; reg87=ponderation*reg152; sollicitation[indices[0]+1]+=reg87; T tmp_3_6=ponderation*reg203; T tmp_6_6=ponderation*reg141;
    T tmp_0_3=ponderation*reg58; T tmp_5_4=ponderation*reg115; T tmp_5_5=ponderation*reg112; T tmp_0_4=ponderation*reg117; T tmp_5_6=ponderation*reg120;
    T tmp_7_6=-reg166; T tmp_5_7=ponderation*reg124; reg58=ponderation*reg176; sollicitation[indices[0]+0]+=reg58; T tmp_7_5=ponderation*reg114;
    T tmp_7_7=ponderation*reg127; T tmp_7_4=ponderation*reg149; T tmp_1_1=ponderation*reg151; T tmp_7_3=-reg136; T tmp_1_0=ponderation*reg153;
    T tmp_0_7=ponderation*reg147; T tmp_7_2=ponderation*reg119; T tmp_6_0=ponderation*reg155; T tmp_6_1=ponderation*reg159; T tmp_6_2=-reg183;
    T tmp_6_3=ponderation*reg162; T tmp_0_5=ponderation*reg132; T tmp_6_4=ponderation*reg134; T tmp_6_7=-reg154; T tmp_7_0=ponderation*reg137;
    T tmp_7_1=ponderation*reg140; T tmp_0_6=ponderation*reg143; T tmp_2_4=ponderation*reg180; T tmp_4_2=ponderation*reg185; T tmp_2_3=-reg139;
    T tmp_4_3=ponderation*reg182; T tmp_4_4=ponderation*reg62; T tmp_4_5=ponderation*reg63; T tmp_1_2=ponderation*reg188; T tmp_1_3=ponderation*reg193;
    T tmp_1_4=ponderation*reg163; T tmp_1_5=ponderation*reg167; T tmp_1_6=ponderation*reg171; reg62=ponderation*reg121; sollicitation[indices[3]+1]+=reg62;
    T tmp_1_7=ponderation*reg70; sollicitation[indices[3]+0]+=-reg178; T tmp_2_0=ponderation*reg177; T tmp_2_1=ponderation*reg72; T tmp_2_2=ponderation*reg73;
    T tmp_0_0=ponderation*reg71; T tmp_4_6=ponderation*reg52; T tmp_4_7=ponderation*reg59; T tmp_0_1=ponderation*reg196; T tmp_5_0=ponderation*reg198;
    T tmp_5_1=ponderation*reg202; T tmp_0_2=ponderation*reg199; T tmp_5_2=ponderation*reg61; T tmp_5_3=ponderation*reg209;
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
    T reg4=reg0*reg1; T reg5=reg3*reg4; T reg6=reg2*reg1; reg1=reg3*reg1; reg4=reg2*reg4;
    T reg7=reg1*reg3; T reg8=reg2*reg4; T reg9=reg6*reg2; reg6=reg6*reg3; reg4=reg3*reg4;
    T reg10=reg3*reg5; reg6=reg7+reg6; reg8=reg8-reg10; reg4=reg10+reg4; reg5=reg2*reg5;
    reg1=reg1*reg2; reg9=reg9-reg7; reg9=reg9*reg2; reg5=reg10+reg5; reg6=reg6*reg3;
    reg10=reg4*reg3; T reg11=reg7+reg1; T reg12=reg8*reg2; reg11=reg11*reg3; reg6=reg9-reg6;
    reg10=reg12-reg10; reg9=reg5*reg3; reg9=reg10-reg9; reg11=reg6-reg11; reg6=1-var_inter[1];
    reg10=1-var_inter[0]; reg11=reg11/reg9; reg12=reg6*elem.pos(0)[0]; T reg13=reg6*elem.pos(1)[0]; T reg14=reg6*elem.pos(0)[1];
    T reg15=reg6*elem.pos(1)[1]; T reg16=elem.pos(0)[0]*reg10; T reg17=elem.pos(1)[0]*var_inter[0]; reg4=reg4/reg9; T reg18=elem.pos(1)[1]*var_inter[0];
    reg8=reg8/reg9; T reg19=reg10*elem.pos(0)[1]; T reg20=reg4*reg11; T reg21=reg8*reg11; T reg22=reg19+reg18;
    T reg23=elem.pos(2)[0]*var_inter[1]; reg13=reg13-reg12; T reg24=elem.pos(2)[1]*var_inter[0]; T reg25=elem.pos(2)[1]*var_inter[1]; T reg26=reg16+reg17;
    reg15=reg15-reg14; T reg27=elem.pos(2)[0]*var_inter[0]; reg9=reg5/reg9; reg5=reg8*reg21; T reg28=reg4*reg20;
    T reg29=reg8*(*f.m).alpha; T reg30=reg4*(*f.m).alpha; T reg31=elem.pos(3)[1]*reg10; reg24=reg24-reg22; T reg32=elem.pos(3)[0]*reg10;
    reg27=reg27-reg26; T reg33=elem.pos(3)[1]*var_inter[1]; reg15=reg25+reg15; reg25=elem.pos(3)[0]*var_inter[1]; reg13=reg23+reg13;
    reg31=reg24+reg31; reg32=reg27+reg32; reg15=reg15-reg33; reg13=reg13-reg25; reg30=reg29+reg30;
    reg9=reg9*(*f.m).alpha; reg28=reg5-reg28; reg5=reg2*reg0; reg23=reg3*reg0; reg24=reg13*reg31;
    reg9=reg30+reg9; reg20=reg20/reg28; reg28=reg21/reg28; reg21=reg15*reg32; reg21=reg24-reg21;
    reg24=pow(reg2,2); reg28=reg28*reg9; reg9=reg20*reg9; reg20=pow(reg3,2); reg27=reg23*reg3;
    reg29=reg5*reg2; reg15=reg15/reg21; reg32=reg32/reg21; reg27=reg29-reg27; reg29=1-(*f.m).resolution;
    reg13=reg13/reg21; reg31=reg31/reg21; reg9=reg28-reg9; reg20=reg24-reg20; reg24=var_inter[1]*reg31;
    reg28=var_inter[0]*reg15; reg30=(*f.m).resolution*(*f.m).alpha; T reg34=var_inter[1]*reg32; reg23=reg23/reg27; reg5=reg5/reg27;
    reg9=reg29*reg9; T reg35=reg6*reg31; reg27=reg20/reg27; reg20=reg10*reg13; T reg36=reg10*reg15;
    T reg37=reg28+reg35; T reg38=var_inter[0]*reg13; reg9=reg30+reg9; reg30=reg20+reg34; T reg39=reg36+reg24;
    T reg40=(*f.m).resolution*reg5; reg8=reg8*reg29; reg11=reg11*reg29; reg29=reg4*reg29; reg4=(*f.m).resolution*reg23;
    T reg41=(*f.m).resolution*reg27; T reg42=reg6*reg32; T reg43=reg38-reg34; T reg44=reg42+reg38; reg11=reg41+reg11;
    reg29=reg4+reg29; reg4=reg42-reg20; reg40=reg8+reg40; reg8=0.5*reg39; reg41=reg24-reg28;
    T reg45=reg36-reg35; T reg46=0.5*reg30; reg9=(*f.m).deltaT*reg9; T reg47=0.5*reg37; T reg48=reg29*reg9;
    T reg49=0.5*reg45; T reg50=0.5*reg44; T reg51=0.5*reg43; T reg52=0.5*reg41; T reg53=0.5*reg4;
    T reg54=reg46*reg11; T reg55=reg40*reg9; T reg56=reg11*reg47; T reg57=reg8*reg11; T reg58=reg53*reg11;
    reg56=2*reg56; T reg59=reg11*reg49; T reg60=reg11*reg50; T reg61=reg29*reg44; T reg62=reg30*reg40;
    T reg63=reg11*reg51; T reg64=reg11*reg52; T reg65=reg40*reg39; reg54=2*reg54; T reg66=reg30*reg29;
    T reg67=2*reg57; T reg68=reg48+reg55; T reg69=reg10*var_inter[1]; T reg70=reg6*var_inter[0]; T reg71=reg40*reg41;
    T reg72=reg69*(*f.m).f_vol[0]; reg63=2*reg63; T reg73=reg29*reg43; T reg74=reg70*(*f.m).f_vol[1]; T reg75=reg62*reg44;
    T reg76=reg61*reg37; T reg77=reg67*reg47; T reg78=reg56*reg50; reg64=2*reg64; T reg79=reg66*reg39;
    T reg80=reg65*reg37; T reg81=reg54*reg50; T reg82=reg68*reg39; T reg83=reg4*reg40; T reg84=reg6*reg10;
    T reg85=var_inter[0]*var_inter[1]; T reg86=reg46*reg67; T reg87=reg68*reg44; reg59=2*reg59; T reg88=reg29*reg41;
    T reg89=reg40*reg37; T reg90=reg40*reg43; T reg91=2*reg60; T reg92=reg4*reg29; T reg93=reg29*reg37;
    T reg94=reg45*reg40; T reg95=reg40*reg44; T reg96=reg29*reg39; reg58=2*reg58; T reg97=reg65*reg39;
    T reg98=reg54*reg53; T reg99=reg94*reg45; reg81=reg80+reg81; T reg100=reg66*reg45; T reg101=reg66*reg37;
    T reg102=reg69*(*f.m).f_vol[1]; T reg103=reg67*reg50; T reg104=reg67*reg53; T reg105=reg92*reg45; T reg106=reg46*reg54;
    T reg107=reg83*reg4; T reg108=reg59*reg49; T reg109=reg58*reg53; T reg110=reg56*reg47; T reg111=reg95*reg44;
    T reg112=reg84*(*f.m).f_vol[0]; T reg113=reg84*(*f.m).f_vol[1]; T reg114=reg56*reg53; T reg115=reg30*reg68; T reg116=reg71*reg45;
    T reg117=reg70*(*f.m).f_vol[0]; T reg118=reg61*reg45; T reg119=reg63*reg53; T reg120=reg8*reg67; T reg121=reg30*reg62;
    T reg122=reg64*reg50; T reg123=reg91*reg53; T reg124=reg73*reg45; T reg125=reg82-reg72; T reg126=reg85*(*f.m).f_vol[0];
    T reg127=reg85*(*f.m).f_vol[1]; T reg128=reg89*reg45; reg78=reg76+reg78; T reg129=reg64*reg53; T reg130=reg71*reg37;
    T reg131=reg63*reg50; T reg132=reg65*reg45; T reg133=reg73*reg37; T reg134=reg59*reg53; reg79=reg79+reg86;
    T reg135=reg67*reg52; T reg136=reg62*reg43; T reg137=reg96*reg44; T reg138=reg91*reg50; reg75=reg77+reg75;
    T reg139=reg63*reg49; T reg140=reg71*reg41; T reg141=reg54*reg52; T reg142=reg63*reg51; T reg143=reg96*reg43;
    T reg144=reg4*reg88; T reg145=reg73*reg41; T reg146=reg64*reg51; T reg147=reg64*reg52; T reg148=reg90*reg43;
    T reg149=reg65*reg41; T reg150=reg54*reg51; T reg151=reg4*reg93; T reg152=reg91*reg49; T reg153=reg67*reg51;
    reg66=reg66*reg41; T reg154=reg4*reg95; T reg155=reg56*reg49; T reg156=reg45*reg68; T reg157=reg87-reg74;
    T reg158=reg68*reg41; reg62=reg4*reg62; T reg159=reg67*reg49; T reg160=reg4*reg68; T reg161=reg54*reg49;
    T reg162=reg68*reg43; T reg163=reg4*reg96; T reg164=reg68*reg37; T reg165=reg64*reg49; T reg166=reg63*reg47;
    T reg167=reg88*reg44; T reg168=reg4*reg90; T reg169=reg64*reg47; T reg170=reg90*reg44; T reg171=reg89*reg37;
    T reg172=reg54*reg47; reg165=reg168+reg165; reg128=reg128-reg123; reg109=reg99+reg109; reg134=reg105+reg134;
    reg139=reg144+reg139; reg155=reg155-reg154; reg161=reg161-reg163; reg114=reg114-reg118; reg106=reg97+reg106;
    reg157=reg157*reg21; reg99=reg158+reg126; reg110=reg110+reg111; reg101=reg101+reg103; reg105=reg162+reg127;
    reg144=reg81*reg21; reg122=reg133-reg122; reg167=reg166-reg167; reg131=reg130-reg131; reg130=reg78*reg21;
    reg170=reg169-reg170; reg172=reg172+reg137; reg133=reg75*reg21; reg142=reg140+reg142; reg146=reg145+reg146;
    reg140=reg102+reg115; reg125=reg125*reg21; reg151=reg151-reg152; reg147=reg148+reg147; reg66=reg66-reg153;
    reg129=reg124+reg129; reg124=reg160+reg113; reg145=reg164+reg117; reg150=reg150-reg149; reg98=reg98-reg132;
    reg121=reg121+reg120; reg171=reg171+reg138; reg141=reg141-reg143; reg148=reg79*reg21; reg100=reg100-reg104;
    reg108=reg107+reg108; reg136=reg136-reg135; reg107=reg112+reg156; reg62=reg62-reg159; reg119=reg116+reg119;
    reg66=reg66*reg21; reg116=reg140*reg21; reg101=reg101*reg21; reg146=reg146*reg21; reg108=reg108*reg21;
    reg131=reg131*reg21; reg147=reg147*reg21; reg142=reg142*reg21; reg166=reg105*reg21; reg136=reg136*reg21;
    reg167=reg167*reg21; reg98=reg98*reg21; reg168=ponderation*reg133; reg169=ponderation*reg144; reg122=reg122*reg21;
    T reg173=ponderation*reg130; T reg174=reg145*reg21; reg172=reg172*reg21; reg141=reg141*reg21; reg100=reg100*reg21;
    reg110=reg110*reg21; reg170=reg170*reg21; T reg175=reg107*reg21; reg119=reg119*reg21; reg62=reg62*reg21;
    T reg176=ponderation*reg148; reg128=reg128*reg21; reg106=reg21*reg106; reg139=reg139*reg21; reg157=ponderation*reg157;
    reg134=reg134*reg21; reg171=reg171*reg21; reg121=reg121*reg21; reg165=reg165*reg21; reg155=reg155*reg21;
    reg109=reg109*reg21; T reg177=reg99*reg21; reg150=reg150*reg21; reg151=reg151*reg21; reg129=reg129*reg21;
    reg161=reg161*reg21; reg114=reg114*reg21; reg125=ponderation*reg125; T reg178=reg124*reg21; T tmp_1_1=ponderation*reg108;
    T tmp_6_6=ponderation*reg106; reg106=ponderation*reg174; sollicitation[indices[1]+0]+=reg106; T tmp_2_5=ponderation*reg122; T tmp_3_3=ponderation*reg110;
    T tmp_7_7=ponderation*reg121; reg108=ponderation*reg166; sollicitation[indices[2]+1]+=reg108; sollicitation[indices[1]+1]+=-reg157; T tmp_2_6=-reg169;
    reg110=ponderation*reg175; sollicitation[indices[0]+0]+=reg110; reg121=ponderation*reg178; sollicitation[indices[0]+1]+=reg121; reg122=ponderation*reg177;
    sollicitation[indices[2]+0]+=reg122; T tmp_2_7=ponderation*reg101; T tmp_0_3=ponderation*reg114; T tmp_1_6=ponderation*reg161; T tmp_1_5=ponderation*reg165;
    T tmp_0_2=ponderation*reg128; T tmp_0_4=ponderation*reg119; T tmp_1_4=ponderation*reg139; T tmp_0_1=ponderation*reg134; T tmp_1_7=ponderation*reg62;
    T tmp_1_3=ponderation*reg155; T tmp_0_0=ponderation*reg109; T tmp_2_2=ponderation*reg171; T tmp_4_6=ponderation*reg150; T tmp_1_2=ponderation*reg151;
    T tmp_0_5=ponderation*reg129; sollicitation[indices[3]+0]+=-reg125; T tmp_4_5=ponderation*reg146; reg62=ponderation*reg116; sollicitation[indices[3]+1]+=reg62;
    T tmp_4_7=ponderation*reg66; T tmp_4_4=ponderation*reg142; T tmp_5_5=ponderation*reg147; T tmp_3_7=-reg168; T tmp_0_6=ponderation*reg98;
    T tmp_3_6=ponderation*reg172; T tmp_5_6=ponderation*reg141; T tmp_3_5=ponderation*reg170; T tmp_6_7=-reg176; T tmp_0_7=ponderation*reg100;
    T tmp_2_3=-reg173; T tmp_3_4=ponderation*reg167; T tmp_2_4=ponderation*reg131; T tmp_5_7=ponderation*reg136;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1-var_inter[1]; T reg3=1-var_inter[0];
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=1.0/(*f.m).elastic_modulus; T reg6=reg0*reg1; T reg7=elem.pos(1)[1]*var_inter[0]; T reg8=reg3*elem.pos(0)[1];
    T reg9=elem.pos(1)[0]*var_inter[0]; T reg10=elem.pos(0)[0]*reg3; T reg11=reg2*elem.pos(1)[1]; T reg12=reg2*elem.pos(0)[1]; T reg13=reg2*elem.pos(1)[0];
    T reg14=reg2*elem.pos(0)[0]; reg11=reg11-reg12; T reg15=elem.pos(2)[0]*var_inter[0]; T reg16=reg5*reg1; T reg17=reg10+reg9;
    T reg18=elem.pos(2)[1]*var_inter[0]; T reg19=elem.pos(2)[1]*var_inter[1]; reg13=reg13-reg14; reg1=reg4*reg1; T reg20=reg8+reg7;
    T reg21=elem.pos(2)[0]*var_inter[1]; T reg22=reg5*reg6; reg6=reg4*reg6; T reg23=reg16*reg4; T reg24=reg1*reg4;
    reg16=reg16*reg5; reg13=reg21+reg13; reg21=elem.pos(3)[0]*var_inter[1]; reg11=reg19+reg11; reg19=elem.pos(3)[1]*var_inter[1];
    T reg25=reg4*reg6; T reg26=elem.pos(3)[1]*reg3; reg18=reg18-reg20; reg15=reg15-reg17; T reg27=elem.pos(3)[0]*reg3;
    T reg28=reg4*reg22; reg22=reg5*reg22; reg26=reg18+reg26; reg6=reg5*reg6; reg28=reg25+reg28;
    reg23=reg24+reg23; reg13=reg13-reg21; reg16=reg16-reg24; reg22=reg22-reg25; reg27=reg15+reg27;
    reg11=reg11-reg19; reg1=reg1*reg5; reg16=reg16*reg5; reg23=reg23*reg4; reg15=reg24+reg1;
    reg18=reg13*reg26; T reg29=reg4*reg0; T reg30=reg28*reg4; T reg31=reg5*reg0; reg6=reg25+reg6;
    reg25=reg11*reg27; T reg32=reg22*reg5; T reg33=pow(reg4,2); T reg34=pow(reg5,2); reg15=reg15*reg4;
    T reg35=reg31*reg5; reg25=reg18-reg25; reg23=reg16-reg23; reg30=reg32-reg30; reg16=reg6*reg4;
    reg18=reg29*reg4; reg16=reg30-reg16; reg15=reg23-reg15; reg18=reg35-reg18; reg33=reg34-reg33;
    reg11=reg11/reg25; reg27=reg27/reg25; reg13=reg13/reg25; reg26=reg26/reg25; reg15=reg15/reg16;
    reg33=reg33/reg18; reg23=1-(*f.m).resolution; reg30=reg2*reg26; reg32=var_inter[1]*reg27; reg34=var_inter[1]*reg26;
    reg35=var_inter[0]*reg13; T reg36=var_inter[0]*reg11; T reg37=reg3*reg11; T reg38=reg3*reg13; T reg39=reg2*reg27;
    T reg40=reg37+reg34; T reg41=reg39+reg35; T reg42=reg36+reg30; reg31=reg31/reg18; reg18=reg29/reg18;
    reg22=reg22/reg16; reg28=reg28/reg16; reg29=(*f.m).resolution*reg33; T reg43=reg15*reg23; T reg44=reg38+reg32;
    T reg45=(*f.m).resolution*reg31; T reg46=reg35-reg32; T reg47=reg34-reg36; T reg48=0.5*reg42; T reg49=0.5*reg41;
    T reg50=reg39-reg38; T reg51=(*f.m).resolution*reg18; T reg52=reg28*reg23; T reg53=reg22*reg23; T reg54=0.5*reg40;
    T reg55=reg37-reg30; reg43=reg29+reg43; reg29=0.5*reg44; T reg56=reg43*reg49; reg45=reg53+reg45;
    reg53=reg43*reg48; T reg57=0.5*reg55; reg52=reg51+reg52; reg51=reg29*reg43; T reg58=reg54*reg43;
    T reg59=0.5*reg46; T reg60=0.5*reg47; T reg61=0.5*reg50; T reg62=reg45*reg40; reg51=2*reg51;
    T reg63=reg44*reg45; T reg64=reg52*reg40; T reg65=reg43*reg60; T reg66=reg43*reg59; reg53=2*reg53;
    T reg67=reg45*reg41; T reg68=reg52*reg41; T reg69=2*reg56; T reg70=reg52*reg42; T reg71=reg45*reg42;
    T reg72=reg43*reg57; T reg73=reg61*reg43; T reg74=2*reg58; T reg75=reg44*reg52; reg66=2*reg66;
    T reg76=reg71*reg40; T reg77=reg52*reg46; T reg78=reg45*reg46; T reg79=reg55*reg52; reg65=2*reg65;
    reg73=2*reg73; T reg80=reg54*reg51; T reg81=reg44*reg64; T reg82=reg68*reg42; T reg83=reg50*reg52;
    T reg84=reg54*reg53; T reg85=reg44*reg67; T reg86=reg69*reg48; T reg87=reg70*reg41; T reg88=reg55*reg45;
    T reg89=reg53*reg49; T reg90=reg29*reg74; reg72=2*reg72; T reg91=reg63*reg41; T reg92=reg74*reg48;
    T reg93=reg50*reg45; T reg94=reg51*reg49; T reg95=reg62*reg42; T reg96=reg75*reg40; T reg97=reg52*reg47;
    T reg98=reg29*reg69; T reg99=reg45*reg47; T reg100=reg93*reg50; T reg101=reg64*reg46; T reg102=reg51*reg60;
    T reg103=reg79*reg50; T reg104=reg63*reg46; T reg105=reg75*reg55; T reg106=reg74*reg60; T reg107=reg29*reg73;
    T reg108=reg51*reg61; reg76=reg98+reg76; T reg109=reg29*reg53; T reg110=reg88*reg40; T reg111=reg54*reg72;
    T reg112=reg74*reg61; T reg113=reg44*reg93; T reg114=reg54*reg73; T reg115=reg29*reg72; T reg116=reg72*reg57;
    T reg117=reg44*reg79; T reg118=reg73*reg57; T reg119=reg99*reg40; T reg120=reg29*reg66; T reg121=reg68*reg40;
    T reg122=reg83*reg40; T reg123=reg73*reg61; reg89=reg82+reg89; T reg124=reg99*reg42; T reg125=reg66*reg49;
    T reg126=reg77*reg42; reg96=reg96+reg90; reg94=reg95+reg94; T reg127=reg75*reg42; T reg128=reg74*reg49;
    T reg129=reg73*reg48; T reg130=reg79*reg41; T reg131=reg72*reg48; T reg132=reg93*reg41; T reg133=reg44*reg70;
    reg87=reg86+reg87; T reg134=reg53*reg48; T reg135=reg67*reg41; T reg136=reg88*reg55; T reg137=reg29*reg51;
    T reg138=reg77*reg40; T reg139=reg29*reg65; T reg140=reg62*reg55; T reg141=reg54*reg69; reg84=reg85+reg84;
    T reg142=reg44*reg97; T reg143=reg54*reg66; T reg144=reg44*reg78; T reg145=reg54*reg65; reg80=reg81+reg80;
    T reg146=reg44*reg63; T reg147=reg54*reg74; T reg148=reg65*reg61; T reg149=reg77*reg55; T reg150=reg65*reg49;
    T reg151=reg66*reg61; T reg152=reg99*reg55; T reg153=reg53*reg61; T reg154=reg68*reg55; T reg155=reg69*reg61;
    T reg156=reg71*reg55; T reg157=reg72*reg61; T reg158=reg83*reg55; T reg159=reg62*reg40; T reg160=reg74*reg59;
    T reg161=reg50*reg67; T reg162=reg69*reg57; T reg163=reg50*reg70; T reg164=reg65*reg59; reg77=reg77*reg47;
    reg79=reg79*reg46; T reg165=reg73*reg60; T reg166=reg66*reg59; reg99=reg99*reg47; T reg167=reg53*reg59;
    T reg168=reg68*reg47; T reg169=reg69*reg59; reg93=reg93*reg46; T reg170=reg72*reg60; T reg171=reg71*reg47;
    T reg172=reg72*reg59; T reg173=reg73*reg49; T reg174=reg88*reg42; T reg175=reg83*reg42; reg72=reg72*reg49;
    T reg176=reg74*reg57; reg71=reg71*reg42; T reg177=reg69*reg49; reg63=reg50*reg63; T reg178=reg51*reg57;
    T reg179=reg50*reg64; T reg180=reg62*reg47; T reg181=reg51*reg59; T reg182=reg65*reg57; T reg183=reg50*reg78;
    T reg184=reg66*reg57; T reg185=reg50*reg97; T reg186=reg53*reg57; reg75=reg75*reg47; T reg187=reg65*reg48;
    T reg188=reg78*reg41; reg53=reg53*reg60; T reg189=reg64*reg41; reg91=reg92+reg91; reg51=reg51*reg48;
    T reg190=reg67*reg46; reg78=reg78*reg46; T reg191=reg97*reg41; T reg192=reg66*reg48; reg88=reg88*reg47;
    reg97=reg97*reg46; reg73=reg73*reg59; reg66=reg66*reg60; reg65=reg65*reg60; T reg193=reg69*reg60;
    reg70=reg70*reg46; reg83=reg83*reg47; reg148=reg149+reg148; reg149=reg96*reg25; T reg194=reg87*reg25;
    T reg195=reg80*reg25; reg182=reg183+reg182; reg188=reg187-reg188; reg134=reg134+reg135; reg178=reg178-reg179;
    reg145=reg144-reg145; reg143=reg142-reg143; reg191=reg192-reg191; reg142=reg84*reg25; reg63=reg63-reg176;
    reg133=reg141+reg133; reg137=reg159+reg137; reg146=reg146+reg147; reg172=reg83+reg172; reg83=reg89*reg25;
    reg125=reg124-reg125; reg171=reg171-reg169; reg73=reg88+reg73; reg167=reg167-reg168; reg123=reg136+reg123;
    reg150=reg126-reg150; reg166=reg99+reg166; reg88=reg91*reg25; reg99=reg94*reg25; reg157=reg158+reg157;
    reg164=reg77+reg164; reg127=reg127+reg128; reg156=reg156-reg155; reg163=reg163-reg162; reg51=reg51+reg189;
    reg130=reg129-reg130; reg153=reg153-reg154; reg186=reg186-reg161; reg151=reg152+reg151; reg184=reg185+reg184;
    reg132=reg131-reg132; reg70=reg70-reg193; reg181=reg181-reg180; reg114=reg117-reg114; reg118=reg103+reg118;
    reg105=reg105-reg112; reg119=reg120-reg119; reg170=reg93+reg170; reg116=reg100+reg116; reg108=reg108-reg140;
    reg111=reg113-reg111; reg53=reg53-reg190; reg165=reg79+reg165; reg66=reg97+reg66; reg109=reg109+reg121;
    reg110=reg107-reg110; reg71=reg71+reg177; reg173=reg174-reg173; reg104=reg104-reg106; reg122=reg115-reg122;
    reg102=reg102-reg101; reg72=reg175-reg72; reg138=reg139-reg138; reg65=reg78+reg65; reg77=reg76*reg25;
    reg75=reg75-reg160; reg78=ponderation*reg83; reg105=reg105*reg25; reg164=reg164*reg25; reg157=reg157*reg25;
    reg170=reg170*reg25; reg172=reg172*reg25; reg122=reg122*reg25; reg108=reg108*reg25; reg171=reg171*reg25;
    reg165=reg165*reg25; reg166=reg166*reg25; reg167=reg167*reg25; reg123=reg123*reg25; reg110=reg110*reg25;
    reg65=reg65*reg25; reg137=reg25*reg137; reg138=reg138*reg25; reg146=reg146*reg25; reg134=reg134*reg25;
    reg102=reg102*reg25; reg191=reg191*reg25; reg104=reg104*reg25; reg79=ponderation*reg194; reg66=reg66*reg25;
    reg188=reg188*reg25; reg132=reg132*reg25; reg53=reg53*reg25; reg130=reg130*reg25; reg51=reg51*reg25;
    reg127=reg127*reg25; reg116=reg116*reg25; reg93=ponderation*reg99; reg70=reg70*reg25; reg97=ponderation*reg88;
    reg150=reg150*reg25; reg118=reg118*reg25; reg125=reg125*reg25; reg73=reg73*reg25; reg186=reg186*reg25;
    reg181=reg181*reg25; reg151=reg151*reg25; reg119=reg119*reg25; reg184=reg184*reg25; reg148=reg148*reg25;
    reg100=ponderation*reg149; reg114=reg114*reg25; reg182=reg182*reg25; reg103=ponderation*reg195; reg111=reg111*reg25;
    reg145=reg145*reg25; reg178=reg178*reg25; reg143=reg143*reg25; reg71=reg71*reg25; reg107=ponderation*reg142;
    reg72=reg72*reg25; reg63=reg63*reg25; reg133=reg133*reg25; reg173=reg173*reg25; reg156=reg156*reg25;
    reg163=reg163*reg25; reg113=ponderation*reg77; reg75=reg75*reg25; reg109=reg109*reg25; reg153=reg153*reg25;
    T tmp_3_2=-reg79; T tmp_0_2=ponderation*reg156; T tmp_5_7=ponderation*reg104; T tmp_7_4=ponderation*reg143; T tmp_1_6=ponderation*reg178;
    T tmp_3_5=ponderation*reg188; T tmp_3_1=ponderation*reg132; T tmp_4_4=ponderation*reg166; T tmp_7_1=ponderation*reg111; T tmp_5_3=ponderation*reg53;
    T tmp_3_0=ponderation*reg130; T tmp_7_5=ponderation*reg145; T tmp_6_0=ponderation*reg110; T tmp_4_5=ponderation*reg164; T tmp_5_5=ponderation*reg65;
    T tmp_6_6=ponderation*reg137; T tmp_2_0=ponderation*reg173; T tmp_6_5=ponderation*reg138; T tmp_4_7=ponderation*reg75; T tmp_7_2=ponderation*reg133;
    T tmp_1_7=ponderation*reg63; T tmp_7_7=ponderation*reg146; T tmp_6_1=ponderation*reg122; T tmp_3_3=ponderation*reg134; T tmp_0_1=ponderation*reg157;
    T tmp_5_6=ponderation*reg102; T tmp_7_3=-reg107; T tmp_3_4=ponderation*reg191; T tmp_5_4=ponderation*reg66; T tmp_2_1=ponderation*reg72;
    T tmp_6_7=-reg100; T tmp_2_4=ponderation*reg125; T tmp_6_4=ponderation*reg119; T tmp_0_7=ponderation*reg105; T tmp_1_4=ponderation*reg184;
    T tmp_4_0=ponderation*reg73; T tmp_2_3=-reg78; T tmp_0_4=ponderation*reg151; T tmp_5_1=ponderation*reg170; T tmp_0_3=ponderation*reg153;
    T tmp_5_0=ponderation*reg165; T tmp_4_1=ponderation*reg172; T tmp_4_6=ponderation*reg181; T tmp_0_6=ponderation*reg108; T tmp_1_3=ponderation*reg186;
    T tmp_6_3=ponderation*reg109; T tmp_4_2=ponderation*reg171; T tmp_1_1=ponderation*reg116; T tmp_2_2=ponderation*reg71; T tmp_2_7=ponderation*reg127;
    T tmp_7_0=ponderation*reg114; T tmp_3_6=ponderation*reg51; T tmp_7_6=-reg103; T tmp_2_6=-reg93; T tmp_1_5=ponderation*reg182;
    T tmp_5_2=ponderation*reg70; T tmp_1_2=ponderation*reg163; T tmp_0_0=ponderation*reg123; T tmp_3_7=-reg97; T tmp_2_5=ponderation*reg150;
    T tmp_4_3=ponderation*reg167; T tmp_1_0=ponderation*reg118; T tmp_0_5=ponderation*reg148; T tmp_6_2=-reg113;
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
    T reg4=reg1*elem.pos(1)[0]; T reg5=elem.pos(0)[0]*reg3; T reg6=elem.pos(1)[0]*var_inter[0]; T reg7=reg1*elem.pos(1)[1]; T reg8=reg3*elem.pos(0)[1];
    T reg9=elem.pos(1)[1]*var_inter[0]; T reg10=1.0/(*f.m).elastic_modulus; T reg11=reg1*elem.pos(0)[0]; T reg12=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg13=reg1*elem.pos(0)[1];
    T reg14=reg0*reg2; T reg15=reg8+reg9; T reg16=elem.pos(2)[1]*var_inter[1]; reg7=reg7-reg13; T reg17=elem.pos(2)[1]*var_inter[0];
    reg4=reg4-reg11; T reg18=elem.pos(2)[0]*var_inter[1]; T reg19=reg10*reg2; T reg20=elem.pos(2)[0]*var_inter[0]; T reg21=reg10*reg14;
    reg2=reg12*reg2; reg14=reg12*reg14; T reg22=reg5+reg6; reg7=reg16+reg7; reg16=elem.pos(3)[1]*var_inter[1];
    T reg23=elem.pos(3)[0]*var_inter[1]; reg4=reg18+reg4; reg18=reg19*reg12; T reg24=reg2*reg12; reg19=reg19*reg10;
    T reg25=reg12*reg14; reg20=reg20-reg22; T reg26=elem.pos(3)[0]*reg3; T reg27=reg12*reg21; reg21=reg10*reg21;
    reg17=reg17-reg15; T reg28=elem.pos(3)[1]*reg3; reg14=reg10*reg14; reg27=reg25+reg27; reg2=reg2*reg10;
    reg28=reg17+reg28; reg4=reg4-reg23; reg19=reg19-reg24; reg21=reg21-reg25; reg26=reg20+reg26;
    reg7=reg7-reg16; reg18=reg24+reg18; reg17=reg10*reg0; reg20=reg12*reg0; T reg29=reg24+reg2;
    T reg30=reg27*reg12; reg18=reg18*reg12; reg19=reg19*reg10; reg14=reg25+reg14; reg25=reg7*reg26;
    T reg31=reg4*reg28; T reg32=reg21*reg10; T reg33=pow(reg10,2); reg30=reg32-reg30; reg32=reg14*reg12;
    T reg34=reg20*reg12; reg18=reg19-reg18; reg19=reg17*reg10; reg29=reg29*reg12; reg25=reg31-reg25;
    reg31=pow(reg12,2); reg7=reg7/reg25; reg26=reg26/reg25; reg4=reg4/reg25; reg28=reg28/reg25;
    reg32=reg30-reg32; reg31=reg33-reg31; reg34=reg19-reg34; reg29=reg18-reg29; reg18=1-(*f.m).resolution;
    reg19=var_inter[1]*reg26; reg31=reg31/reg34; reg30=var_inter[0]*reg7; reg33=var_inter[1]*reg28; T reg35=reg1*reg28;
    reg29=reg29/reg32; T reg36=reg3*reg4; T reg37=reg3*reg7; T reg38=var_inter[0]*reg4; reg27=reg27/reg32;
    T reg39=reg30+reg35; reg21=reg21/reg32; T reg40=reg1*reg26; T reg41=reg36+reg19; T reg42=reg29*reg18;
    T reg43=reg37+reg33; reg17=reg17/reg34; T reg44=(*f.m).resolution*reg31; reg34=reg20/reg34; reg42=reg44+reg42;
    reg20=(*f.m).resolution*reg34; reg44=reg27*reg18; T reg45=0.5*reg41; T reg46=reg21*reg18; T reg47=(*f.m).resolution*reg17;
    T reg48=reg40+reg38; T reg49=reg40-reg36; T reg50=reg37-reg35; T reg51=0.5*reg43; T reg52=0.5*reg39;
    T reg53=reg38-reg19; T reg54=reg33-reg30; T reg55=0.5*reg49; T reg56=reg42*reg52; T reg57=0.5*reg54;
    reg47=reg46+reg47; reg46=0.5*reg53; reg44=reg20+reg44; reg20=reg45*reg42; T reg58=reg51*reg42;
    T reg59=0.5*reg50; T reg60=0.5*reg48; T reg61=reg42*reg46; reg56=2*reg56; T reg62=reg42*reg57;
    T reg63=reg41*reg47; T reg64=reg44*reg48; T reg65=reg42*reg60; T reg66=reg47*reg43; T reg67=reg42*reg59;
    reg20=2*reg20; T reg68=reg41*reg44; T reg69=2*reg58; T reg70=reg55*reg42; T reg71=reg20*reg60;
    T reg72=reg64*reg39; reg62=2*reg62; T reg73=reg49*reg47; T reg74=reg44*reg53; T reg75=reg44*reg43;
    T reg76=reg68*reg43; T reg77=reg45*reg69; T reg78=reg56*reg60; reg70=2*reg70; T reg79=reg69*reg52;
    T reg80=reg63*reg48; T reg81=reg49*reg44; reg67=2*reg67; T reg82=reg44*reg39; T reg83=reg47*reg39;
    T reg84=2*reg65; T reg85=reg47*reg48; T reg86=reg44*reg54; reg61=2*reg61; T reg87=reg47*reg54;
    T reg88=reg50*reg47; T reg89=reg47*reg53; T reg90=reg66*reg39; T reg91=reg68*reg50; T reg92=reg74*reg50;
    T reg93=reg69*reg55; T reg94=reg70*reg55; T reg95=reg73*reg49; T reg96=reg67*reg59; T reg97=reg62*reg60;
    T reg98=reg61*reg55; reg78=reg72+reg78; T reg99=reg87*reg39; T reg100=reg61*reg60; T reg101=reg74*reg39;
    T reg102=reg66*reg43; T reg103=reg20*reg55; reg76=reg76+reg77; reg71=reg90+reg71; T reg104=reg81*reg50;
    T reg105=reg66*reg50; T reg106=reg87*reg50; T reg107=reg67*reg55; T reg108=reg62*reg55; T reg109=reg51*reg69;
    T reg110=reg83*reg50; T reg111=reg41*reg63; T reg112=reg84*reg55; T reg113=reg64*reg50; T reg114=reg56*reg55;
    T reg115=reg87*reg54; T reg116=reg62*reg57; T reg117=reg89*reg53; T reg118=reg61*reg46; T reg119=reg61*reg59;
    T reg120=reg45*reg20; T reg121=reg74*reg54; T reg122=reg62*reg46; T reg123=reg49*reg86; T reg124=reg69*reg46;
    T reg125=reg68*reg54; T reg126=reg56*reg52; T reg127=reg56*reg59; T reg128=reg20*reg46; T reg129=reg66*reg54;
    T reg130=reg85*reg48; T reg131=reg49*reg82; T reg132=reg84*reg59; T reg133=reg84*reg60; T reg134=reg83*reg39;
    T reg135=reg69*reg59; T reg136=reg49*reg85; reg68=reg68*reg39; T reg137=reg62*reg59; T reg138=reg61*reg52;
    T reg139=reg69*reg60; T reg140=reg86*reg48; T reg141=reg49*reg89; T reg142=reg69*reg57; T reg143=reg63*reg53;
    T reg144=reg49*reg75; T reg145=reg62*reg52; T reg146=reg89*reg48; T reg147=reg20*reg59; T reg148=reg20*reg52;
    T reg149=reg75*reg48; T reg150=reg20*reg57; T reg151=reg75*reg53; reg80=reg79+reg80; T reg152=reg88*reg50;
    reg63=reg49*reg63; reg127=reg127-reg136; reg119=reg123+reg119; reg114=reg114-reg113; reg137=reg141+reg137;
    reg110=reg110-reg112; reg131=reg131-reg132; reg97=reg101-reg97; reg107=reg104+reg107; reg122=reg121+reg122;
    reg118=reg115+reg118; reg94=reg152+reg94; reg120=reg102+reg120; reg101=reg80*reg25; reg148=reg148+reg149;
    reg146=reg145-reg146; reg104=reg78*reg25; reg140=reg138-reg140; reg100=reg99-reg100; reg116=reg117+reg116;
    reg99=reg71*reg25; reg103=reg103-reg105; reg91=reg91-reg93; reg125=reg125-reg124; reg128=reg128-reg129;
    reg150=reg150-reg151; reg134=reg134+reg133; reg126=reg126+reg130; reg96=reg95+reg96; reg63=reg63-reg135;
    reg95=reg76*reg25; reg98=reg106+reg98; reg68=reg68+reg139; reg108=reg92+reg108; reg143=reg143-reg142;
    reg111=reg111+reg109; reg147=reg147-reg144; reg116=reg116*reg25; reg111=reg111*reg25; reg68=reg68*reg25;
    reg92=ponderation*reg101; reg100=reg100*reg25; reg91=reg91*reg25; reg120=reg25*reg120; reg148=reg148*reg25;
    reg106=ponderation*reg99; reg126=reg126*reg25; reg143=reg143*reg25; reg115=ponderation*reg104; reg146=reg146*reg25;
    reg96=reg96*reg25; reg150=reg150*reg25; reg140=reg140*reg25; reg98=reg98*reg25; reg137=reg137*reg25;
    reg119=reg119*reg25; reg114=reg114*reg25; reg147=reg147*reg25; reg108=reg108*reg25; reg127=reg127*reg25;
    reg117=ponderation*reg95; reg97=reg97*reg25; reg110=reg110*reg25; reg63=reg63*reg25; reg131=reg131*reg25;
    reg134=reg134*reg25; reg122=reg122*reg25; reg107=reg107*reg25; reg128=reg128*reg25; reg125=reg125*reg25;
    reg118=reg118*reg25; reg103=reg103*reg25; reg94=reg94*reg25; T tmp_6_6=ponderation*reg120; T tmp_2_6=-reg106;
    T tmp_2_7=ponderation*reg68; T tmp_2_5=ponderation*reg97; T tmp_3_3=ponderation*reg126; T tmp_3_4=ponderation*reg140; T tmp_3_5=ponderation*reg146;
    T tmp_3_6=ponderation*reg148; T tmp_3_7=-reg92; T tmp_4_4=ponderation*reg118; T tmp_4_5=ponderation*reg122; T tmp_1_2=ponderation*reg131;
    T tmp_1_3=ponderation*reg127; T tmp_1_4=ponderation*reg119; T tmp_1_5=ponderation*reg137; T tmp_1_6=ponderation*reg147; T tmp_1_7=ponderation*reg63;
    T tmp_2_2=ponderation*reg134; T tmp_4_6=ponderation*reg128; T tmp_4_7=ponderation*reg125; T tmp_5_5=ponderation*reg116; T tmp_5_6=ponderation*reg150;
    T tmp_5_7=ponderation*reg143; T tmp_7_7=ponderation*reg111; T tmp_1_1=ponderation*reg96; T tmp_0_7=ponderation*reg91; T tmp_0_6=ponderation*reg103;
    T tmp_6_7=-reg117; T tmp_0_5=ponderation*reg108; T tmp_0_4=ponderation*reg98; T tmp_0_3=ponderation*reg114; T tmp_0_2=ponderation*reg110;
    T tmp_0_1=ponderation*reg107; T tmp_0_0=ponderation*reg94; T tmp_2_3=-reg115; T tmp_2_4=ponderation*reg100;
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
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg4*reg2; T reg6=reg4*reg1; reg2=reg3*reg2; reg1=reg3*reg1;
    T reg7=reg3*reg2; T reg8=reg3*reg5; T reg9=reg6*reg3; reg5=reg4*reg5; T reg10=reg1*reg3;
    reg6=reg6*reg4; reg2=reg4*reg2; reg8=reg7+reg8; reg5=reg5-reg7; reg9=reg10+reg9;
    reg1=reg1*reg4; reg6=reg6-reg10; T reg11=reg8*reg3; reg2=reg7+reg2; reg6=reg6*reg4;
    reg7=reg10+reg1; T reg12=reg5*reg4; reg9=reg9*reg3; reg7=reg7*reg3; reg11=reg12-reg11;
    reg12=reg2*reg3; reg9=reg6-reg9; reg12=reg11-reg12; reg7=reg9-reg7; reg8=reg8/reg12;
    reg5=reg5/reg12; reg7=reg7/reg12; reg6=reg5*reg7; reg9=reg8*reg7; reg12=reg2/reg12;
    reg2=reg5*reg6; reg11=reg8*reg9; T reg13=1-var_inter[0]; T reg14=reg5*(*f.m).alpha; T reg15=reg8*(*f.m).alpha;
    T reg16=1-var_inter[1]; T reg17=elem.pos(1)[1]*var_inter[0]; T reg18=reg13*elem.pos(0)[1]; T reg19=elem.pos(1)[0]*var_inter[0]; T reg20=elem.pos(0)[0]*reg13;
    reg12=reg12*(*f.m).alpha; reg15=reg14+reg15; reg11=reg2-reg11; reg2=reg16*elem.pos(0)[0]; reg14=reg16*elem.pos(1)[0];
    T reg21=reg16*elem.pos(0)[1]; T reg22=reg16*elem.pos(1)[1]; T reg23=elem.pos(2)[1]*var_inter[0]; T reg24=elem.pos(2)[1]*var_inter[1]; reg14=reg14-reg2;
    T reg25=reg20+reg19; T reg26=reg18+reg17; T reg27=elem.pos(2)[0]*var_inter[1]; T reg28=reg3*reg0; T reg29=reg4*reg0;
    reg12=reg15+reg12; reg6=reg6/reg11; reg11=reg9/reg11; reg9=elem.pos(2)[0]*var_inter[0]; reg22=reg22-reg21;
    reg11=reg11*reg12; reg12=reg6*reg12; reg6=reg28*reg3; reg14=reg27+reg14; reg15=elem.pos(3)[0]*var_inter[1];
    reg27=elem.pos(3)[1]*reg13; reg23=reg23-reg26; reg22=reg24+reg22; reg24=elem.pos(3)[1]*var_inter[1]; T reg30=elem.pos(3)[0]*reg13;
    reg9=reg9-reg25; T reg31=reg29*reg4; reg11=reg12-reg11; reg12=1-(*f.m).resolution; reg30=reg9+reg30;
    reg27=reg23+reg27; reg22=reg22-reg24; reg14=reg14-reg15; reg6=reg31-reg6; reg9=reg14*reg27;
    reg28=reg28/reg6; reg11=reg12*reg11; reg23=reg22*reg30; reg31=(*f.m).resolution*(*f.m).alpha; reg29=reg29/reg6;
    reg11=reg31+reg11; reg31=(*f.m).resolution*reg29; reg5=reg5*reg12; reg8=reg8*reg12; reg23=reg9-reg23;
    reg9=(*f.m).resolution*reg28; reg31=reg5+reg31; reg8=reg9+reg8; reg22=reg22/reg23; reg11=(*f.m).deltaT*reg11;
    reg30=reg30/reg23; reg14=reg14/reg23; reg27=reg27/reg23; reg5=reg16*reg30; reg9=reg8*reg11;
    T reg32=reg31*reg11; T reg33=var_inter[1]*reg27; T reg34=var_inter[0]*reg14; T reg35=reg13*reg22; T reg36=reg35+reg33;
    T reg37=reg13*reg14; T reg38=var_inter[1]*reg30; T reg39=reg13*var_inter[1]; T reg40=reg5+reg34; T reg41=var_inter[0]*reg22;
    T reg42=reg16*var_inter[0]; T reg43=reg16*reg27; T reg44=reg9+reg32; T reg45=reg16*reg13; T reg46=var_inter[0]*var_inter[1];
    T reg47=reg44*reg40; T reg48=reg41+reg43; T reg49=reg42*(*f.m).f_vol[1]; T reg50=reg37+reg38; T reg51=reg39*(*f.m).f_vol[0];
    T reg52=reg33-reg41; T reg53=reg34-reg38; T reg54=reg35-reg43; T reg55=reg5-reg37; T reg56=reg44*reg36;
    T reg57=reg46*(*f.m).f_vol[0]; T reg58=reg46*(*f.m).f_vol[1]; T reg59=reg44*reg52; T reg60=reg42*(*f.m).f_vol[0]; T reg61=reg47-reg49;
    T reg62=reg44*reg53; T reg63=reg45*(*f.m).f_vol[1]; T reg64=reg50*reg44; T reg65=reg45*(*f.m).f_vol[0]; T reg66=reg54*reg44;
    T reg67=reg55*reg44; T reg68=reg39*(*f.m).f_vol[1]; T reg69=reg44*reg48; T reg70=reg56-reg51; T reg71=reg62+reg58;
    reg70=reg70*reg23; T reg72=reg59+reg57; T reg73=reg65+reg66; T reg74=reg67+reg63; reg61=reg61*reg23;
    T reg75=reg68+reg64; T reg76=reg69+reg60; reg61=ponderation*reg61; reg70=ponderation*reg70; T reg77=reg73*reg23;
    T reg78=reg76*reg23; T reg79=reg72*reg23; T reg80=reg75*reg23; T reg81=reg71*reg23; T reg82=reg74*reg23;
    T reg83=ponderation*reg78; sollicitation[indices[1]+0]+=reg83; sollicitation[indices[3]+0]+=-reg70; sollicitation[indices[1]+1]+=-reg61; reg61=ponderation*reg77;
    sollicitation[indices[0]+0]+=reg61; reg70=ponderation*reg82; sollicitation[indices[0]+1]+=reg70; T reg84=ponderation*reg79; sollicitation[indices[2]+0]+=reg84;
    T reg85=ponderation*reg80; sollicitation[indices[3]+1]+=reg85; T reg86=ponderation*reg81; sollicitation[indices[2]+1]+=reg86;
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
    T reg4=reg0*reg1; T reg5=reg2*reg1; reg1=reg3*reg1; T reg6=reg3*reg4; reg4=reg2*reg4;
    T reg7=reg5*reg2; T reg8=reg1*reg3; T reg9=reg3*reg4; reg5=reg5*reg3; T reg10=reg3*reg6;
    reg4=reg2*reg4; reg5=reg8+reg5; reg4=reg4-reg10; reg9=reg10+reg9; reg6=reg2*reg6;
    reg1=reg1*reg2; reg7=reg7-reg8; T reg11=reg9*reg3; reg7=reg7*reg2; T reg12=reg8+reg1;
    reg6=reg10+reg6; reg5=reg5*reg3; reg10=reg4*reg2; reg12=reg12*reg3; reg5=reg7-reg5;
    reg7=reg6*reg3; reg11=reg10-reg11; reg7=reg11-reg7; reg10=1-var_inter[1]; reg12=reg5-reg12;
    reg5=1-var_inter[0]; reg9=reg9/reg7; reg4=reg4/reg7; reg11=reg10*elem.pos(0)[1]; T reg13=reg10*elem.pos(1)[1];
    T reg14=elem.pos(0)[0]*reg5; T reg15=elem.pos(1)[0]*var_inter[0]; T reg16=reg5*elem.pos(0)[1]; reg12=reg12/reg7; T reg17=reg10*elem.pos(1)[0];
    T reg18=reg10*elem.pos(0)[0]; T reg19=elem.pos(1)[1]*var_inter[0]; reg17=reg17-reg18; T reg20=reg16+reg19; T reg21=elem.pos(2)[1]*var_inter[0];
    T reg22=elem.pos(2)[1]*var_inter[1]; T reg23=elem.pos(2)[0]*var_inter[1]; T reg24=reg4*reg12; reg13=reg13-reg11; T reg25=elem.pos(2)[0]*var_inter[0];
    T reg26=reg9*reg12; T reg27=reg14+reg15; reg7=reg6/reg7; reg6=reg4*(*f.m).alpha; reg17=reg23+reg17;
    reg23=reg4*reg24; T reg28=reg9*(*f.m).alpha; T reg29=reg9*reg26; T reg30=elem.pos(3)[1]*reg5; reg21=reg21-reg20;
    T reg31=elem.pos(3)[0]*reg5; reg25=reg25-reg27; T reg32=elem.pos(3)[1]*var_inter[1]; reg13=reg22+reg13; reg22=elem.pos(3)[0]*var_inter[1];
    T reg33=reg10*vectors[0][indices[0]+1]; reg29=reg23-reg29; reg23=reg10*vectors[0][indices[1]+1]; T reg34=reg5*vectors[0][indices[0]+0]; reg28=reg6+reg28;
    reg6=reg10*vectors[0][indices[0]+0]; reg7=reg7*(*f.m).alpha; T reg35=var_inter[0]*vectors[0][indices[1]+0]; T reg36=var_inter[0]*vectors[0][indices[1]+1]; T reg37=reg5*vectors[0][indices[0]+1];
    reg17=reg17-reg22; T reg38=reg10*vectors[0][indices[1]+0]; reg13=reg13-reg32; reg31=reg25+reg31; reg30=reg21+reg30;
    reg26=reg26/reg29; reg34=reg35+reg34; reg33=reg23-reg33; reg29=reg24/reg29; reg21=var_inter[1]*vectors[0][indices[2]+1];
    reg23=reg17*reg30; reg24=var_inter[0]*vectors[0][indices[2]+0]; reg7=reg28+reg7; reg36=reg37+reg36; reg25=var_inter[1]*vectors[0][indices[2]+0];
    reg6=reg38-reg6; reg28=reg13*reg31; reg35=var_inter[0]*vectors[0][indices[2]+1]; reg34=reg24-reg34; reg26=reg26*reg7;
    reg24=var_inter[1]*vectors[0][indices[3]+1]; reg33=reg21+reg33; reg6=reg25+reg6; reg21=var_inter[1]*vectors[0][indices[3]+0]; reg25=reg5*vectors[0][indices[3]+0];
    reg37=reg5*vectors[0][indices[3]+1]; reg36=reg35-reg36; reg28=reg23-reg28; reg7=reg29*reg7; reg23=reg2*reg0;
    reg29=reg3*reg0; reg35=pow(reg3,2); reg38=reg23*reg2; reg21=reg6-reg21; reg6=pow(reg2,2);
    T reg39=1-(*f.m).resolution; reg24=reg33-reg24; reg33=reg29*reg3; reg37=reg36+reg37; reg34=reg25+reg34;
    reg13=reg13/reg28; reg31=reg31/reg28; reg26=reg7-reg26; reg17=reg17/reg28; reg30=reg30/reg28;
    reg7=reg24*reg30; reg25=reg13*reg37; reg36=reg17*reg34; T reg40=reg21*reg31; reg33=reg38-reg33;
    reg38=(*f.m).resolution*(*f.m).alpha; reg26=reg39*reg26; reg35=reg6-reg35; reg34=reg13*reg34; reg37=reg17*reg37;
    reg26=reg38+reg26; reg24=reg24*reg31; reg29=reg29/reg33; reg25=reg7-reg25; reg40=reg36-reg40;
    reg23=reg23/reg33; reg21=reg21*reg30; reg33=reg35/reg33; reg24=reg37-reg24; reg4=reg4*reg39;
    reg6=(*f.m).resolution*reg33; reg7=(*f.m).resolution*reg23; reg34=reg21-reg34; reg40=reg25+reg40; reg26=(*f.m).deltaT*reg26;
    reg21=(*f.m).resolution*reg29; reg9=reg9*reg39; reg39=reg12*reg39; reg24=reg24-reg26; reg34=reg34-reg26;
    reg39=reg6+reg39; reg9=reg21+reg9; reg7=reg4+reg7; reg4=var_inter[0]*reg13; reg6=var_inter[0]*reg17;
    reg12=var_inter[1]*reg31; reg21=var_inter[1]*reg30; reg25=reg5*reg13; reg40=0.5*reg40; reg35=reg10*reg30;
    reg36=reg10*reg31; reg37=reg5*reg17; reg38=reg4+reg35; T reg41=reg25-reg35; reg40=reg40*reg39;
    T reg42=reg24*reg9; T reg43=reg36+reg6; T reg44=reg21-reg4; reg24=reg24*reg7; T reg45=reg7*reg34;
    T reg46=reg25+reg21; T reg47=reg37+reg12; T reg48=reg36-reg37; T reg49=reg6-reg12; reg34=reg9*reg34;
    T reg50=0.5*reg46; T reg51=0.5*reg43; T reg52=0.5*reg44; T reg53=0.5*reg38; T reg54=0.5*reg49;
    reg40=2*reg40; reg24=reg34+reg24; reg34=0.5*reg41; reg45=reg42+reg45; reg42=0.5*reg47;
    T reg55=0.5*reg48; T reg56=reg10*reg5; T reg57=var_inter[0]*var_inter[1]; T reg58=reg45*reg44; T reg59=reg45*reg46;
    T reg60=reg5*var_inter[1]; T reg61=reg40*reg54; T reg62=reg24*reg49; T reg63=reg40*reg52; T reg64=reg42*reg40;
    T reg65=reg47*reg24; T reg66=reg10*var_inter[0]; T reg67=reg40*reg51; T reg68=reg45*reg38; T reg69=reg40*reg55;
    T reg70=reg24*reg43; T reg71=reg40*reg34; T reg72=reg24*reg48; T reg73=reg40*reg53; T reg74=reg45*reg41;
    T reg75=reg50*reg40; reg64=reg64-reg59; T reg76=reg60*(*f.m).f_vol[0]; reg71=reg72+reg71; reg72=reg56*(*f.m).f_vol[1];
    reg65=reg65-reg75; reg69=reg74+reg69; reg74=reg60*(*f.m).f_vol[1]; reg68=reg68-reg67; T reg77=reg66*(*f.m).f_vol[0];
    reg73=reg73-reg70; T reg78=reg56*(*f.m).f_vol[0]; reg61=reg58+reg61; reg58=reg57*(*f.m).f_vol[0]; T reg79=reg57*(*f.m).f_vol[1];
    reg63=reg62+reg63; reg62=reg66*(*f.m).f_vol[1]; reg64=reg64-reg76; reg61=reg61-reg58; reg73=reg73-reg62;
    reg69=reg69-reg78; reg71=reg71-reg72; reg65=reg65-reg74; reg63=reg63-reg79; reg68=reg68-reg77;
    reg73=reg73*reg28; reg65=reg65*reg28; reg61=reg61*reg28; reg69=reg69*reg28; reg64=reg64*reg28;
    reg71=reg71*reg28; reg68=reg68*reg28; reg63=reg63*reg28; sollicitation[indices[3]+1]+=ponderation*reg65; sollicitation[indices[0]+0]+=ponderation*reg69;
    sollicitation[indices[0]+1]+=ponderation*reg71; sollicitation[indices[3]+0]+=ponderation*reg64; sollicitation[indices[2]+1]+=ponderation*reg63; sollicitation[indices[1]+0]+=ponderation*reg68; sollicitation[indices[2]+0]+=ponderation*reg61;
    sollicitation[indices[1]+1]+=ponderation*reg73;
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

