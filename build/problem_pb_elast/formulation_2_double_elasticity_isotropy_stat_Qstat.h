
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=reg1*reg0;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg3*reg2; reg3=reg4*reg3; T reg6=reg1*reg2; reg1=reg4*reg1;
    T reg7=reg4*reg3; T reg8=reg4*reg5; T reg9=reg4*reg6; T reg10=reg4*reg1; reg6=reg2*reg6;
    reg5=reg5*reg2; reg8=reg7+reg8; reg3=reg3*reg2; reg1=reg2*reg1; reg6=reg6-reg10;
    reg9=reg10+reg9; reg5=reg5-reg7; reg6=reg2*reg6; reg9=reg4*reg9; reg3=reg7+reg3;
    reg7=reg5*reg2; T reg11=reg10+reg1; T reg12=reg4*reg8; T reg13=reg4*reg3; T reg14=elem.pos(2)[1]-elem.pos(0)[1];
    T reg15=elem.pos(1)[0]-elem.pos(0)[0]; T reg16=elem.pos(1)[1]-elem.pos(0)[1]; T reg17=elem.pos(2)[0]-elem.pos(0)[0]; reg12=reg7-reg12; reg11=reg4*reg11;
    reg9=reg6-reg9; reg11=reg9-reg11; reg13=reg12-reg13; reg6=reg15*reg14; reg7=reg16*reg17;
    reg8=reg8/reg13; reg7=reg6-reg7; reg11=reg11/reg13; reg5=reg5/reg13; reg14=reg14/reg7;
    reg15=reg15/reg7; reg6=reg2*reg0; reg17=reg17/reg7; reg9=reg8*reg11; reg12=vectors[0][indices[2]+1]-vectors[0][indices[0]+1];
    T reg18=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; T reg19=reg5*reg11; reg16=reg16/reg7; T reg20=reg4*reg0; T reg21=vectors[0][indices[2]+0]-vectors[0][indices[0]+0];
    T reg22=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; T reg23=reg9*reg8; T reg24=(*f.m).alpha*reg5; T reg25=(*f.m).alpha*reg8; T reg26=reg15*reg12;
    T reg27=reg17*reg18; T reg28=reg16*reg21; T reg29=reg22*reg14; reg13=reg3/reg13; reg3=reg4*reg20;
    T reg30=reg5*reg19; T reg31=reg6*reg2; reg28=reg29-reg28; elem.epsilon[0][0]=reg28; reg23=reg30-reg23;
    reg25=reg24+reg25; reg13=(*f.m).alpha*reg13; reg24=(*f.m).alpha*(*f.m).deltaT; reg27=reg26-reg27; elem.epsilon[0][1]=reg27;
    reg3=reg31-reg3; reg13=reg25+reg13; reg20=reg20/reg3; reg25=reg27-reg24; reg9=reg9/reg23;
    reg26=reg28-reg24; reg23=reg19/reg23; reg6=reg6/reg3; reg19=reg25*reg20; reg23=reg23*reg13;
    reg13=reg9*reg13; reg9=reg26*reg6; reg26=reg26*reg20; reg25=reg25*reg6; reg29=1-(*f.m).resolution;
    reg13=reg23-reg13; reg19=reg9+reg19; reg25=reg26+reg25; reg9=reg2*reg25; reg23=PNODE(1).dep[1]-PNODE(0).dep[1];
    reg26=reg4*reg19; reg30=PNODE(2).dep[1]-PNODE(0).dep[1]; reg25=reg4*reg25; reg19=reg2*reg19; reg13=reg29*reg13;
    reg31=(*f.m).resolution*(*f.m).alpha; T reg32=PNODE(1).dep[0]-PNODE(0).dep[0]; T reg33=pow(reg4,2); T reg34=PNODE(2).dep[0]-PNODE(0).dep[0]; T reg35=pow(reg2,2);
    reg33=reg35-reg33; reg19=reg19-reg25; reg9=reg9-reg26; reg13=reg31+reg13; reg31=reg32*reg14;
    reg35=reg16*reg34; reg32=reg17*reg32; T reg36=reg17*reg23; T reg37=reg15*reg30; reg30=reg16*reg30;
    reg23=reg23*reg14; reg34=reg15*reg34; reg13=(*f.m).deltaT*reg13; reg30=reg23-reg30; reg23=(*f.m).resolution*reg20;
    reg26=reg25+reg26; reg9=reg24+reg9; reg5=reg5*reg29; reg19=reg24+reg19; reg32=reg34-reg32;
    reg36=reg37-reg36; reg35=reg31-reg35; reg25=(*f.m).resolution*reg6; reg8=reg29*reg8; reg3=reg33/reg3;
    reg21=reg15*reg21; reg22=reg17*reg22; reg12=reg16*reg12; reg18=reg18*reg14; reg31=(*f.m).resolution*reg3;
    reg30=reg32+reg30; reg8=reg23+reg8; reg12=reg18-reg12; reg22=reg21-reg22; reg25=reg5+reg25;
    reg11=reg29*reg11; reg5=reg19+reg9; reg26=reg24-reg26; reg18=reg35-reg13; reg21=reg36-reg13;
    reg11=reg31+reg11; reg5=reg5+reg26; reg23=reg8*reg21; reg29=reg25*reg18; reg12=reg22+reg12;
    reg22=reg25*reg21; reg30=0.5*reg30; reg31=reg8*reg18; reg32=reg30*reg11; reg22=reg31+reg22;
    reg12=0.5*reg12; elem.epsilon[0][2]=reg12; reg5=reg5/3; reg23=reg29+reg23; reg21=reg21*reg22;
    reg32=2*reg32; reg19=reg19-reg5; reg9=reg9-reg5; reg18=reg18*reg23; reg29=reg12*reg3;
    reg19=pow(reg19,2); reg9=pow(reg9,2); reg5=reg26-reg5; reg29=reg0*reg29; reg26=reg30*reg32;
    reg21=reg18+reg21; reg18=2*reg29; reg9=reg19+reg9; reg5=pow(reg5,2); reg26=reg21+reg26;
    reg18=reg29*reg18; reg5=reg9+reg5; reg26=reg7*reg26; reg7=0.16666666666666665741*reg26; reg26=0.33333333333333331483*reg26;
    reg27=reg27-reg13; reg13=reg28-reg13; reg18=reg5+reg18; reg5=reg13*reg25; reg26=reg7+reg26;
    reg7=reg27*reg8; reg8=reg13*reg8; reg25=reg27*reg25; reg18=1.5*reg18; elem.sigma[0][0]=reg5+reg7;
    elem.sigma[0][1]=reg8+reg25; elem.sigma[0][2]=reg12*reg11; elem.sigma_von_mises=pow(reg18,0.5); elem.ener=reg26/2;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=reg1*reg0;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg1*reg2; reg1=reg4*reg1; T reg6=reg4*reg3; reg3=reg3*reg2;
    T reg7=reg4*reg1; T reg8=reg2*reg5; T reg9=reg4*reg3; T reg10=reg4*reg6; reg3=reg3*reg2;
    reg5=reg4*reg5; reg5=reg7+reg5; reg3=reg3-reg10; reg8=reg8-reg7; reg9=reg10+reg9;
    reg1=reg2*reg1; reg6=reg6*reg2; T reg11=reg4*reg9; reg5=reg4*reg5; reg8=reg2*reg8;
    T reg12=reg3*reg2; T reg13=reg7+reg1; reg6=reg10+reg6; reg13=reg4*reg13; reg5=reg8-reg5;
    reg11=reg12-reg11; reg8=reg4*reg6; reg8=reg11-reg8; reg13=reg5-reg13; reg3=reg3/reg8;
    reg9=reg9/reg8; reg13=reg13/reg8; reg5=reg3*reg13; reg10=reg9*reg13; reg11=reg3*reg5;
    reg12=reg10*reg9; T reg14=(*f.m).alpha*reg3; T reg15=(*f.m).alpha*reg9; reg8=reg6/reg8; reg12=reg11-reg12;
    reg15=reg14+reg15; reg8=(*f.m).alpha*reg8; reg6=reg2*reg0; reg5=reg5/reg12; reg12=reg10/reg12;
    reg8=reg15+reg8; reg10=reg4*reg0; reg11=elem.pos(2)[1]-elem.pos(0)[1]; reg14=elem.pos(2)[0]-elem.pos(0)[0]; reg15=elem.pos(1)[1]-elem.pos(0)[1];
    T reg16=elem.pos(1)[0]-elem.pos(0)[0]; T reg17=pow(reg4,2); T reg18=reg4*reg10; reg5=reg5*reg8; reg8=reg12*reg8;
    reg12=reg6*reg2; T reg19=pow(reg2,2); T reg20=reg16*reg11; reg18=reg12-reg18; reg12=1-(*f.m).resolution;
    reg8=reg5-reg8; reg5=reg15*reg14; reg17=reg19-reg17; reg8=reg12*reg8; reg19=(*f.m).resolution*(*f.m).alpha;
    reg5=reg20-reg5; reg10=reg10/reg18; reg17=reg17/reg18; reg18=reg6/reg18; reg11=reg11/reg5;
    reg8=reg19+reg8; reg13=reg12*reg13; reg9=reg12*reg9; reg6=(*f.m).resolution*reg18; reg12=reg3*reg12;
    reg16=reg16/reg5; reg14=reg14/reg5; reg15=reg15/reg5; reg3=(*f.m).resolution*reg10; reg19=(*f.m).resolution*reg17;
    reg13=reg19+reg13; reg9=reg3+reg9; reg3=0.5*reg14; reg8=(*f.m).deltaT*reg8; reg19=reg14-reg16;
    reg20=reg15-reg11; reg6=reg12+reg6; reg12=0.5*reg11; T reg21=0.5*reg16; T reg22=0.5*reg15;
    T reg23=reg6*reg8; T reg24=reg9*reg8; T reg25=0.5*reg20; T reg26=0.5*reg19; T reg27=reg3*reg13;
    T reg28=reg22*reg13; T reg29=reg12*reg13; T reg30=reg21*reg13; T reg31=reg23+reg24; T reg32=reg26*reg13;
    T reg33=reg25*reg13; T reg34=reg11*reg6; T reg35=2*reg28; T reg36=reg16*reg9; reg30=2*reg30;
    T reg37=reg15*reg6; reg29=2*reg29; T reg38=reg14*reg9; T reg39=reg11*reg9; T reg40=reg14*reg6;
    T reg41=reg15*reg9; T reg42=reg16*reg6; T reg43=2*reg27; T reg44=1-var_inter[0]; T reg45=reg14*reg31;
    T reg46=reg20*reg9; reg44=reg44-var_inter[1]; T reg47=reg6*reg19; T reg48=reg29*reg22; T reg49=reg16*reg41;
    T reg50=reg30*reg22; T reg51=reg16*reg40; T reg52=reg36*reg15; T reg53=reg35*reg21; T reg54=reg15*reg31;
    T reg55=reg38*reg11; T reg56=reg29*reg3; T reg57=reg37*reg11; T reg58=reg30*reg3; T reg59=reg43*reg12;
    T reg60=reg14*reg39; T reg61=var_inter[0]*(*f.m).f_vol[1]; T reg62=var_inter[1]*(*f.m).f_vol[0]; T reg63=reg6*reg20; reg32=2*reg32;
    T reg64=reg19*reg9; reg33=2*reg33; T reg65=reg35*reg12; T reg66=reg14*reg42; T reg67=reg15*reg34;
    T reg68=reg43*reg21; T reg69=reg35*reg3; T reg70=reg36*reg11; T reg71=reg37*reg20; reg58=reg57+reg58;
    reg67=reg68+reg67; reg56=reg55+reg56; T reg72=reg31*reg19; T reg73=reg11*reg64; T reg74=reg33*reg3;
    T reg75=reg29*reg26; T reg76=reg43*reg3; T reg77=reg11*reg34; T reg78=reg16*reg31; T reg79=reg43*reg26;
    T reg80=reg38*reg20; T reg81=reg54-reg62; T reg82=reg14*reg41; T reg83=reg14*reg40; T reg84=reg45-reg61;
    T reg85=reg29*reg12; T reg86=reg31*reg11; T reg87=reg31*reg20; T reg88=reg35*reg26; T reg89=reg42*reg19;
    reg34=reg34*reg20; T reg90=reg21*reg32; T reg91=reg35*reg25; T reg92=reg33*reg26; T reg93=reg15*reg63;
    T reg94=reg64*reg20; reg36=reg36*reg20; T reg95=reg32*reg26; reg66=reg65+reg66; T reg96=reg21*reg33;
    T reg97=reg63*reg20; reg64=reg15*reg64; T reg98=var_inter[1]*(*f.m).f_vol[1]; T reg99=reg30*reg21; reg63=reg11*reg63;
    T reg100=reg30*reg26; T reg101=reg47*reg14; T reg102=reg12*reg33; T reg103=reg32*reg3; T reg104=reg46*reg14;
    T reg105=reg12*reg32; T reg106=reg46*reg16; T reg107=reg44*(*f.m).f_vol[1]; T reg108=reg30*reg25; T reg109=reg22*reg32;
    reg30=reg30*reg12; reg50=reg49+reg50; T reg110=reg39*reg19; T reg111=var_inter[0]*(*f.m).f_vol[0]; reg42=reg16*reg42;
    T reg112=reg35*reg22; T reg113=reg41*reg19; T reg114=reg43*reg25; T reg115=reg29*reg21; T reg116=reg47*reg16;
    T reg117=reg22*reg33; reg60=reg59+reg60; T reg118=reg43*reg22; reg39=reg16*reg39; T reg119=reg40*reg19;
    T reg120=reg44*(*f.m).f_vol[0]; reg29=reg29*reg25; reg33=reg33*reg25; reg47=reg47*reg19; T reg121=reg37*reg15;
    T reg122=reg38*reg15; reg48=reg48+reg51; reg52=reg53+reg52; reg32=reg32*reg25; reg46=reg46*reg19;
    reg101=reg102-reg101; reg99=reg121+reg99; reg117=reg116-reg117; reg29=reg29-reg119; reg102=reg72+reg107;
    reg100=reg100-reg71; reg92=reg94+reg92; reg85=reg85+reg83; reg64=reg96-reg64; reg36=reg36-reg88;
    reg115=reg115+reg122; reg109=reg106-reg109; reg93=reg90-reg93; reg90=reg5*reg66; reg89=reg89-reg91;
    reg94=reg5*reg52; reg95=reg97+reg95; reg108=reg108-reg113; reg33=reg47+reg33; reg47=reg5*reg67;
    reg81=reg5*reg81; reg96=reg5*reg48; reg84=reg5*reg84; reg97=reg78+reg98; reg30=reg30+reg82;
    reg34=reg34-reg79; reg106=reg87+reg120; reg74=reg73-reg74; reg77=reg77+reg76; reg73=reg5*reg50;
    reg116=reg5*reg56; reg75=reg75-reg80; T reg123=reg86+reg111; reg42=reg42+reg112; T reg124=reg5*reg58;
    reg32=reg46+reg32; reg46=reg5*reg60; reg70=reg70+reg69; reg110=reg110-reg114; reg103=reg63-reg103;
    reg39=reg39+reg118; reg104=reg105-reg104; reg34=reg34*reg5; reg63=ponderation*reg90; reg29=reg5*reg29;
    reg33=reg5*reg33; reg85=reg5*reg85; reg30=reg5*reg30; reg75=reg75*reg5; reg100=reg100*reg5;
    reg108=reg5*reg108; reg110=reg5*reg110; reg81=ponderation*reg81; reg105=reg5*reg106; reg84=ponderation*reg84;
    T reg125=reg5*reg123; reg32=reg32*reg5; T reg126=ponderation*reg96; T reg127=ponderation*reg47; T reg128=ponderation*reg46;
    T reg129=ponderation*reg73; reg42=reg5*reg42; reg39=reg5*reg39; reg117=reg5*reg117; reg109=reg5*reg109;
    T reg130=ponderation*reg94; T reg131=reg5*reg102; reg99=reg5*reg99; reg115=reg5*reg115; reg92=reg5*reg92;
    reg93=reg5*reg93; reg95=reg5*reg95; reg89=reg5*reg89; reg36=reg36*reg5; reg64=reg5*reg64;
    reg101=reg5*reg101; reg104=reg5*reg104; reg70=reg5*reg70; reg103=reg5*reg103; T reg132=ponderation*reg124;
    T reg133=ponderation*reg116; reg77=reg5*reg77; reg74=reg5*reg74; T reg134=reg5*reg97; T tmp_1_4=ponderation*reg108;
    T tmp_3_2=-reg128; T tmp_1_5=ponderation*reg89; T tmp_2_0=ponderation*reg103; reg89=ponderation*reg131; sollicitation[indices[0]+1]+=reg89;
    reg103=ponderation*reg105; sollicitation[indices[0]+0]+=reg103; T tmp_2_2=ponderation*reg77; T tmp_2_1=ponderation*reg74; reg74=ponderation*reg134;
    sollicitation[indices[2]+1]+=reg74; sollicitation[indices[2]+0]+=-reg81; sollicitation[indices[1]+1]+=-reg84; reg77=ponderation*reg125; sollicitation[indices[1]+0]+=reg77;
    T tmp_1_0=ponderation*reg32; T tmp_5_5=ponderation*reg42; T tmp_5_3=-reg126; T tmp_4_2=-reg127; T tmp_5_4=-reg129;
    T tmp_5_2=ponderation*reg39; T tmp_5_1=ponderation*reg117; T tmp_5_0=ponderation*reg109; T tmp_4_5=-reg130; T tmp_4_4=ponderation*reg99;
    T tmp_4_3=ponderation*reg115; T tmp_1_3=ponderation*reg29; T tmp_3_3=ponderation*reg85; T tmp_1_2=ponderation*reg110; T tmp_1_1=ponderation*reg33;
    T tmp_0_2=ponderation*reg34; T tmp_0_3=ponderation*reg75; T tmp_3_4=ponderation*reg30; T tmp_0_4=ponderation*reg100; T tmp_3_5=-reg63;
    T tmp_0_1=ponderation*reg92; T tmp_0_0=ponderation*reg95; T tmp_4_0=ponderation*reg93; T tmp_0_5=ponderation*reg36; T tmp_4_1=ponderation*reg64;
    T tmp_3_1=ponderation*reg101; T tmp_3_0=ponderation*reg104; T tmp_2_5=ponderation*reg70; T tmp_2_4=-reg132; T tmp_2_3=-reg133;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=reg1*reg0;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg1*reg2; reg1=reg4*reg1; T reg6=reg4*reg3; reg3=reg3*reg2;
    T reg7=reg4*reg1; T reg8=reg2*reg5; reg5=reg4*reg5; T reg9=reg4*reg3; T reg10=reg4*reg6;
    reg3=reg3*reg2; reg6=reg6*reg2; reg9=reg10+reg9; reg3=reg3-reg10; reg5=reg7+reg5;
    reg8=reg8-reg7; reg1=reg2*reg1; T reg11=reg4*reg9; T reg12=reg7+reg1; T reg13=reg3*reg2;
    reg5=reg4*reg5; reg8=reg2*reg8; reg6=reg10+reg6; reg5=reg8-reg5; reg12=reg4*reg12;
    reg11=reg13-reg11; reg8=reg4*reg6; reg12=reg5-reg12; reg8=reg11-reg8; reg9=reg9/reg8;
    reg3=reg3/reg8; reg12=reg12/reg8; reg5=reg9*reg12; reg10=reg3*reg12; reg8=reg6/reg8;
    reg6=(*f.m).alpha*reg9; reg11=(*f.m).alpha*reg3; reg13=reg5*reg9; T reg14=reg3*reg10; reg8=(*f.m).alpha*reg8;
    reg6=reg11+reg6; reg13=reg14-reg13; reg11=reg4*reg0; reg10=reg10/reg13; reg13=reg5/reg13;
    reg5=reg2*reg0; reg8=reg6+reg8; reg6=pow(reg2,2); reg14=reg5*reg2; T reg15=reg4*reg11;
    T reg16=elem.pos(1)[1]-elem.pos(0)[1]; T reg17=elem.pos(2)[0]-elem.pos(0)[0]; T reg18=elem.pos(2)[1]-elem.pos(0)[1]; T reg19=pow(reg4,2); T reg20=elem.pos(1)[0]-elem.pos(0)[0];
    reg10=reg10*reg8; reg8=reg13*reg8; reg15=reg14-reg15; reg13=1-(*f.m).resolution; reg14=reg16*reg17;
    reg19=reg6-reg19; reg6=reg20*reg18; reg8=reg10-reg8; reg8=reg13*reg8; reg5=reg5/reg15;
    reg19=reg19/reg15; reg10=(*f.m).resolution*(*f.m).alpha; reg15=reg11/reg15; reg14=reg6-reg14; reg6=(*f.m).resolution*reg15;
    reg11=(*f.m).resolution*reg19; T reg21=(*f.m).resolution*reg5; reg20=reg20/reg14; reg18=reg18/reg14; reg17=reg17/reg14;
    reg16=reg16/reg14; reg9=reg13*reg9; reg3=reg3*reg13; reg12=reg13*reg12; reg8=reg10+reg8;
    reg21=reg3+reg21; reg8=(*f.m).deltaT*reg8; reg9=reg6+reg9; reg3=0.5*reg16; reg12=reg11+reg12;
    reg6=0.5*reg20; reg10=0.5*reg18; reg11=reg16-reg18; reg13=reg17-reg20; T reg22=reg6*reg12;
    T reg23=0.5*reg17; T reg24=reg21*reg8; T reg25=reg9*reg8; T reg26=0.5*reg11; T reg27=0.5*reg13;
    T reg28=reg10*reg12; T reg29=reg3*reg12; T reg30=reg16*reg21; T reg31=reg27*reg12; T reg32=reg20*reg21;
    T reg33=reg20*reg9; T reg34=1-var_inter[0]; reg28=2*reg28; T reg35=reg24+reg25; reg22=2*reg22;
    T reg36=reg26*reg12; T reg37=reg17*reg9; T reg38=reg23*reg12; T reg39=2*reg29; T reg40=reg18*reg21;
    T reg41=reg13*reg9; reg36=2*reg36; reg31=2*reg31; T reg42=reg21*reg11; T reg43=reg17*reg35;
    T reg44=reg16*reg35; T reg45=var_inter[1]*(*f.m).f_vol[0]; T reg46=var_inter[0]*(*f.m).f_vol[1]; reg34=reg34-var_inter[1]; T reg47=reg22*reg23;
    T reg48=reg30*reg18; T reg49=reg37*reg18; T reg50=reg28*reg23; T reg51=2*reg38; T reg52=reg18*reg9;
    T reg53=reg16*reg9; T reg54=reg21*reg13; T reg55=reg39*reg10; T reg56=reg17*reg32; T reg57=reg39*reg6;
    T reg58=reg17*reg21; T reg59=reg33*reg16; T reg60=reg28*reg10; T reg61=reg40*reg11; T reg62=reg17*reg58;
    T reg63=reg37*reg11; T reg64=reg39*reg23; T reg65=reg33*reg18; T reg66=reg22*reg6; T reg67=reg51*reg27;
    reg47=reg48+reg47; T reg68=reg39*reg26; T reg69=reg22*reg10; T reg70=reg35*reg13; T reg71=reg34*(*f.m).f_vol[0];
    T reg72=reg35*reg11; T reg73=reg34*(*f.m).f_vol[1]; T reg74=var_inter[1]*(*f.m).f_vol[1]; T reg75=var_inter[0]*(*f.m).f_vol[0]; T reg76=reg42*reg11;
    T reg77=reg39*reg3; T reg78=reg20*reg32; T reg79=reg31*reg27; T reg80=reg54*reg13; T reg81=reg36*reg27;
    T reg82=reg41*reg11; T reg83=reg39*reg27; reg33=reg33*reg11; T reg84=reg28*reg27; T reg85=reg52*reg13;
    T reg86=reg30*reg11; T reg87=reg51*reg26; T reg88=reg35*reg18; T reg89=reg36*reg26; reg59=reg57+reg59;
    T reg90=reg43-reg46; T reg91=reg58*reg13; T reg92=reg44-reg45; T reg93=reg28*reg26; T reg94=reg30*reg16;
    T reg95=reg17*reg53; reg32=reg32*reg13; reg50=reg49+reg50; reg56=reg55+reg56; T reg96=reg22*reg26;
    T reg97=reg51*reg23; T reg98=reg18*reg40; T reg99=reg53*reg13; T reg100=reg22*reg27; T reg101=reg20*reg35;
    reg85=reg85-reg87; reg61=reg61-reg67; reg96=reg96-reg99; reg81=reg82+reg81; reg93=reg93-reg91;
    reg32=reg32-reg68; reg89=reg80+reg89; reg84=reg84-reg63; reg80=reg88+reg75; reg82=reg14*reg59;
    reg90=reg14*reg90; reg66=reg94+reg66; reg92=reg14*reg92; T reg102=reg14*reg56; T reg103=reg101+reg74;
    reg98=reg98+reg97; reg69=reg69+reg95; T reg104=reg14*reg50; reg60=reg60+reg62; T reg105=reg14*reg47;
    T reg106=reg70+reg73; reg100=reg100-reg86; T reg107=reg72+reg71; reg33=reg33-reg83; reg65=reg65+reg64;
    reg79=reg76+reg79; reg78=reg78+reg77; reg76=ponderation*reg82; reg85=reg14*reg85; reg100=reg100*reg14;
    reg66=reg14*reg66; reg33=reg33*reg14; reg93=reg14*reg93; T reg108=ponderation*reg102; reg78=reg14*reg78;
    reg96=reg14*reg96; reg69=reg14*reg69; T reg109=reg14*reg107; reg60=reg14*reg60; reg32=reg14*reg32;
    T reg110=reg14*reg106; reg84=reg84*reg14; reg79=reg14*reg79; T reg111=ponderation*reg104; reg89=reg14*reg89;
    reg81=reg14*reg81; T reg112=reg14*reg80; T reg113=ponderation*reg105; reg90=ponderation*reg90; reg98=reg14*reg98;
    reg65=reg14*reg65; reg61=reg61*reg14; reg92=ponderation*reg92; T reg114=reg14*reg103; T tmp_2_3=-reg111;
    reg111=ponderation*reg110; sollicitation[indices[0]+1]+=reg111; T reg115=ponderation*reg109; sollicitation[indices[0]+0]+=reg115; T tmp_0_2=ponderation*reg61;
    T tmp_2_4=-reg113; T tmp_5_5=ponderation*reg78; T tmp_2_5=ponderation*reg65; T tmp_0_5=ponderation*reg33; T tmp_0_1=ponderation*reg81;
    T tmp_0_4=ponderation*reg100; T tmp_0_0=ponderation*reg79; T tmp_1_5=ponderation*reg32; T tmp_2_2=ponderation*reg98; T tmp_3_3=ponderation*reg60;
    reg32=ponderation*reg114; sollicitation[indices[2]+1]+=reg32; T tmp_1_4=ponderation*reg96; T tmp_3_4=ponderation*reg69; sollicitation[indices[2]+0]+=-reg92;
    T tmp_3_5=-reg108; T tmp_1_3=ponderation*reg93; sollicitation[indices[1]+1]+=-reg90; reg33=ponderation*reg112; sollicitation[indices[1]+0]+=reg33;
    T tmp_4_4=ponderation*reg66; T tmp_1_2=ponderation*reg85; T tmp_1_1=ponderation*reg89; T tmp_4_5=-reg76; T tmp_0_3=ponderation*reg84;
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
    T reg4=reg1*reg0; T reg5=reg3*reg1; reg1=reg1*reg2; T reg6=reg4*reg2; reg4=reg3*reg4;
    T reg7=reg2*reg1; reg1=reg3*reg1; T reg8=reg3*reg5; T reg9=reg6*reg2; T reg10=reg3*reg4;
    reg6=reg3*reg6; reg7=reg7-reg8; reg5=reg2*reg5; reg9=reg9-reg10; reg1=reg8+reg1;
    reg6=reg10+reg6; reg4=reg4*reg2; reg1=reg3*reg1; T reg11=reg3*reg0; reg7=reg2*reg7;
    T reg12=reg2*reg0; reg4=reg10+reg4; reg10=reg3*reg6; T reg13=reg9*reg2; T reg14=reg8+reg5;
    T reg15=elem.pos(2)[1]-elem.pos(0)[1]; T reg16=elem.pos(2)[0]-elem.pos(0)[0]; T reg17=elem.pos(1)[1]-elem.pos(0)[1]; T reg18=elem.pos(1)[0]-elem.pos(0)[0]; T reg19=pow(reg3,2);
    T reg20=reg3*reg11; T reg21=pow(reg2,2); T reg22=reg12*reg2; reg14=reg3*reg14; reg1=reg7-reg1;
    reg10=reg13-reg10; reg7=reg3*reg4; reg14=reg1-reg14; reg19=reg21-reg19; reg20=reg22-reg20;
    reg7=reg10-reg7; reg1=reg18*reg15; reg10=reg17*reg16; reg13=1-(*f.m).resolution; reg19=reg19/reg20;
    reg14=reg14/reg7; reg10=reg1-reg10; reg11=reg11/reg20; reg20=reg12/reg20; reg6=reg6/reg7;
    reg1=reg13*reg14; reg15=reg15/reg10; reg17=reg17/reg10; reg16=reg16/reg10; reg18=reg18/reg10;
    reg12=(*f.m).resolution*reg19; reg9=reg9/reg7; reg21=reg9*reg13; reg22=(*f.m).resolution*reg11; T reg23=0.5*reg17;
    reg1=reg12+reg1; reg12=0.5*reg18; T reg24=reg17-reg15; T reg25=0.5*reg15; T reg26=reg16-reg18;
    T reg27=0.5*reg16; T reg28=reg13*reg6; T reg29=(*f.m).resolution*reg20; T reg30=reg23*reg1; reg28=reg22+reg28;
    reg22=0.5*reg24; T reg31=0.5*reg26; T reg32=reg12*reg1; reg29=reg21+reg29; reg21=reg27*reg1;
    T reg33=reg25*reg1; T reg34=reg18*reg28; T reg35=reg22*reg1; T reg36=reg18*reg29; T reg37=reg16*reg28;
    T reg38=reg17*reg28; T reg39=reg15*reg28; T reg40=reg15*reg29; T reg41=reg31*reg1; reg33=2*reg33;
    T reg42=reg16*reg29; reg32=2*reg32; T reg43=2*reg30; T reg44=reg17*reg29; T reg45=2*reg21;
    T reg46=reg18*reg38; T reg47=reg32*reg23; T reg48=reg17*reg40; T reg49=reg45*reg12; T reg50=reg18*reg42;
    T reg51=reg24*reg28; T reg52=reg33*reg23; reg35=2*reg35; T reg53=reg29*reg26; T reg54=reg26*reg28;
    T reg55=reg45*reg25; T reg56=reg37*reg15; T reg57=reg16*reg36; T reg58=reg43*reg25; T reg59=reg33*reg27;
    T reg60=reg16*reg39; T reg61=reg34*reg17; reg41=2*reg41; T reg62=reg44*reg15; T reg63=reg32*reg27;
    T reg64=reg43*reg12; T reg65=reg29*reg24; T reg66=reg23*reg41; T reg67=reg44*reg24; T reg68=reg34*reg24;
    T reg69=reg53*reg18; T reg70=reg23*reg35; T reg71=reg45*reg23; T reg72=reg18*reg39; T reg73=reg32*reg31;
    T reg74=reg25*reg35; T reg75=reg53*reg16; T reg76=reg51*reg16; T reg77=reg25*reg41; T reg78=reg43*reg31;
    T reg79=reg43*reg27; reg34=reg34*reg15; T reg80=reg65*reg24; reg63=reg62+reg63; T reg81=reg41*reg31;
    reg59=reg56+reg59; T reg82=reg54*reg24; T reg83=reg45*reg27; T reg84=reg15*reg40; T reg85=reg35*reg31;
    T reg86=reg51*reg26; T reg87=reg41*reg22; reg40=reg40*reg24; T reg88=reg12*reg35; T reg89=reg38*reg26;
    T reg90=reg32*reg22; T reg91=reg17*reg65; T reg92=reg12*reg41; T reg93=reg36*reg26; T reg94=reg43*reg22;
    reg57=reg58+reg57; T reg95=reg32*reg12; reg65=reg15*reg65; T reg96=reg16*reg38; reg32=reg32*reg25;
    reg41=reg41*reg27; T reg97=reg15*reg54; T reg98=reg16*reg42; T reg99=reg33*reg25; T reg100=reg35*reg27;
    reg60=reg55+reg60; reg47=reg46+reg47; T reg101=reg45*reg31; reg36=reg18*reg36; T reg102=reg43*reg23;
    T reg103=reg44*reg17; reg35=reg35*reg22; T reg104=reg37*reg24; T reg105=reg37*reg17; T reg106=reg33*reg12;
    reg39=reg39*reg26; T reg107=reg45*reg22; reg61=reg64+reg61; reg53=reg53*reg26; T reg108=reg33*reg31;
    reg52=reg52+reg50; T reg109=reg42*reg26; reg48=reg49+reg48; reg33=reg33*reg22; reg51=reg51*reg18;
    reg54=reg17*reg54; reg108=reg108-reg104; reg75=reg74-reg75; reg68=reg68-reg78; reg40=reg40-reg101;
    reg100=reg97-reg100; reg73=reg73-reg67; reg35=reg53+reg35; reg41=reg65-reg41; reg81=reg80+reg81;
    reg93=reg93-reg94; reg39=reg39-reg107; reg90=reg90-reg89; reg33=reg33-reg109; reg85=reg82+reg85;
    reg66=reg51-reg66; reg32=reg32+reg96; reg70=reg69-reg70; reg51=reg10*reg61; reg72=reg72+reg71;
    reg95=reg103+reg95; reg53=reg10*reg57; reg87=reg86+reg87; reg99=reg99+reg98; reg106=reg106+reg105;
    reg84=reg84+reg83; reg65=reg10*reg60; reg69=reg10*reg47; reg74=reg10*reg59; reg54=reg88-reg54;
    reg76=reg77-reg76; reg36=reg36+reg102; reg91=reg92-reg91; reg34=reg34+reg79; reg77=reg10*reg48;
    reg80=reg10*reg63; reg82=reg10*reg52; reg91=reg10*reg91; reg33=reg10*reg33; reg32=reg10*reg32;
    reg108=reg108*reg10; reg86=ponderation*reg51; reg54=reg10*reg54; reg88=ponderation*reg77; reg93=reg10*reg93;
    reg95=reg10*reg95; reg90=reg10*reg90; reg92=ponderation*reg53; reg40=reg40*reg10; reg106=reg10*reg106;
    reg97=ponderation*reg82; reg35=reg10*reg35; reg39=reg10*reg39; reg76=reg10*reg76; reg75=reg10*reg75;
    reg36=reg10*reg36; reg68=reg68*reg10; reg34=reg10*reg34; T reg110=ponderation*reg80; T reg111=ponderation*reg69;
    T reg112=ponderation*reg74; reg81=reg10*reg81; reg84=reg10*reg84; T reg113=ponderation*reg65; reg100=reg10*reg100;
    reg85=reg10*reg85; reg87=reg87*reg10; reg66=reg10*reg66; reg70=reg10*reg70; reg41=reg10*reg41;
    reg72=reg10*reg72; reg99=reg10*reg99; reg73=reg73*reg10; T tmp_5_5=ponderation*reg36; T tmp_1_5=ponderation*reg93;
    T tmp_3_3=ponderation*reg99; T tmp_5_4=-reg111; T tmp_1_4=ponderation*reg90; T tmp_2_0=ponderation*reg41; T tmp_3_5=-reg92;
    T tmp_3_4=ponderation*reg32; T tmp_3_2=-reg113; T tmp_3_0=ponderation*reg76; T tmp_3_1=ponderation*reg75; T tmp_2_5=ponderation*reg34;
    T tmp_2_4=-reg110; T tmp_2_3=-reg112; T tmp_0_0=ponderation*reg81; T tmp_2_2=ponderation*reg84; T tmp_2_1=ponderation*reg100;
    T tmp_0_5=ponderation*reg68; T tmp_0_1=ponderation*reg85; T tmp_1_0=ponderation*reg87; T tmp_5_2=ponderation*reg72; T tmp_5_1=ponderation*reg70;
    T tmp_0_4=ponderation*reg73; T tmp_5_0=ponderation*reg66; T tmp_4_5=-reg86; T tmp_0_3=ponderation*reg108; T tmp_4_4=ponderation*reg95;
    T tmp_4_3=ponderation*reg106; T tmp_0_2=ponderation*reg40; T tmp_1_1=ponderation*reg35; T tmp_4_2=-reg88; T tmp_1_2=ponderation*reg39;
    T tmp_5_3=-reg97; T tmp_4_1=ponderation*reg54; T tmp_1_3=ponderation*reg33; T tmp_4_0=ponderation*reg91;
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
    T reg4=reg1*reg0; T reg5=reg3*reg1; reg1=reg1*reg2; T reg6=reg4*reg2; reg4=reg3*reg4;
    T reg7=reg2*reg1; reg1=reg3*reg1; T reg8=reg3*reg5; T reg9=reg6*reg2; T reg10=reg3*reg4;
    reg6=reg3*reg6; reg7=reg7-reg8; reg5=reg2*reg5; reg9=reg9-reg10; reg1=reg8+reg1;
    reg6=reg10+reg6; reg4=reg4*reg2; reg1=reg3*reg1; T reg11=reg3*reg0; reg7=reg2*reg7;
    T reg12=reg2*reg0; reg4=reg10+reg4; reg10=reg8+reg5; T reg13=reg3*reg6; T reg14=reg9*reg2;
    T reg15=reg12*reg2; T reg16=pow(reg2,2); T reg17=reg3*reg11; T reg18=pow(reg3,2); T reg19=elem.pos(2)[1]-elem.pos(0)[1];
    T reg20=elem.pos(2)[0]-elem.pos(0)[0]; T reg21=elem.pos(1)[1]-elem.pos(0)[1]; T reg22=elem.pos(1)[0]-elem.pos(0)[0]; reg10=reg3*reg10; T reg23=reg3*reg4;
    reg1=reg7-reg1; reg13=reg14-reg13; reg23=reg13-reg23; reg18=reg16-reg18; reg17=reg15-reg17;
    reg7=reg21*reg20; reg13=reg22*reg19; reg10=reg1-reg10; reg7=reg13-reg7; reg18=reg18/reg17;
    reg1=1-(*f.m).resolution; reg10=reg10/reg23; reg12=reg12/reg17; reg19=reg19/reg7; reg17=reg11/reg17;
    reg6=reg6/reg23; reg11=reg1*reg10; reg22=reg22/reg7; reg20=reg20/reg7; reg21=reg21/reg7;
    reg13=(*f.m).resolution*reg18; reg9=reg9/reg23; reg11=reg13+reg11; reg13=reg9*reg1; reg14=reg20-reg22;
    reg15=reg21-reg19; reg16=(*f.m).resolution*reg12; T reg24=reg1*reg6; T reg25=(*f.m).resolution*reg17; T reg26=0.5*reg19;
    T reg27=0.5*reg22; T reg28=0.5*reg21; reg16=reg13+reg16; reg24=reg25+reg24; reg13=0.5*reg20;
    reg25=reg28*reg11; T reg29=reg26*reg11; T reg30=0.5*reg15; T reg31=0.5*reg14; T reg32=reg27*reg11;
    T reg33=reg20*reg24; T reg34=reg22*reg16; reg29=2*reg29; T reg35=2*reg25; T reg36=reg31*reg11;
    T reg37=reg22*reg24; T reg38=reg13*reg11; reg32=2*reg32; T reg39=reg30*reg11; T reg40=reg21*reg16;
    T reg41=reg16*reg15; reg39=2*reg39; T reg42=reg14*reg24; T reg43=reg19*reg16; reg36=2*reg36;
    T reg44=reg32*reg13; T reg45=reg40*reg19; T reg46=reg19*reg24; T reg47=reg29*reg13; T reg48=reg33*reg19;
    T reg49=reg20*reg16; T reg50=reg16*reg14; T reg51=reg35*reg26; T reg52=reg20*reg34; T reg53=2*reg38;
    T reg54=reg35*reg27; T reg55=reg21*reg24; T reg56=reg37*reg21; T reg57=reg40*reg21; T reg58=reg29*reg31;
    T reg59=reg35*reg31; T reg60=reg34*reg14; T reg61=reg35*reg30; T reg62=reg41*reg15; T reg63=reg29*reg26;
    T reg64=reg20*reg49; T reg65=reg36*reg31; T reg66=reg40*reg15; T reg67=reg32*reg27; T reg68=reg53*reg31;
    T reg69=reg42*reg15; T reg70=reg32*reg26; T reg71=reg20*reg55; T reg72=reg39*reg31; T reg73=reg32*reg31;
    reg52=reg51+reg52; T reg74=reg43*reg15; T reg75=reg37*reg15; T reg76=reg19*reg43; T reg77=reg53*reg13;
    T reg78=reg29*reg30; T reg79=reg53*reg30; T reg80=reg46*reg14; T reg81=reg55*reg14; reg47=reg48+reg47;
    reg56=reg54+reg56; T reg82=reg49*reg14; T reg83=reg35*reg28; T reg84=reg39*reg30; reg44=reg45+reg44;
    T reg85=reg50*reg14; T reg86=reg32*reg30; T reg87=reg33*reg15; reg37=reg37*reg19; reg34=reg22*reg34;
    T reg88=reg35*reg13; reg73=reg73-reg66; reg78=reg78-reg82; reg80=reg80-reg79; reg60=reg60-reg61;
    reg84=reg85+reg84; reg74=reg74-reg68; reg58=reg58-reg87; reg86=reg86-reg81; reg85=reg7*reg56;
    T reg89=reg7*reg44; reg37=reg37+reg88; reg34=reg34+reg83; reg67=reg57+reg67; reg75=reg75-reg59;
    reg63=reg63+reg64; reg65=reg62+reg65; reg76=reg76+reg77; reg72=reg69+reg72; reg62=reg7*reg47;
    reg69=reg7*reg52; reg70=reg70+reg71; reg80=reg7*reg80; reg75=reg75*reg7; reg76=reg7*reg76;
    reg78=reg7*reg78; T reg90=ponderation*reg85; reg67=reg7*reg67; reg86=reg7*reg86; reg34=reg7*reg34;
    reg60=reg7*reg60; reg63=reg7*reg63; T reg91=ponderation*reg69; reg70=reg7*reg70; reg58=reg58*reg7;
    reg73=reg73*reg7; reg84=reg7*reg84; T reg92=ponderation*reg62; reg72=reg7*reg72; reg37=reg7*reg37;
    reg74=reg74*reg7; reg65=reg7*reg65; T reg93=ponderation*reg89; T tmp_5_5=ponderation*reg34; T tmp_0_3=ponderation*reg58;
    T tmp_1_5=ponderation*reg60; T tmp_0_0=ponderation*reg65; T tmp_3_3=ponderation*reg63; T tmp_3_5=-reg91; T tmp_0_1=ponderation*reg72;
    T tmp_3_4=ponderation*reg70; T tmp_0_4=ponderation*reg73; T tmp_2_5=ponderation*reg37; T tmp_1_4=ponderation*reg86; T tmp_2_4=-reg93;
    T tmp_4_4=ponderation*reg67; T tmp_0_2=ponderation*reg74; T tmp_1_3=ponderation*reg78; T tmp_4_5=-reg90; T tmp_2_3=-reg92;
    T tmp_0_5=ponderation*reg75; T tmp_1_2=ponderation*reg80; T tmp_1_1=ponderation*reg84; T tmp_2_2=ponderation*reg76;
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
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg4*reg1; reg1=reg1*reg2; T reg6=reg4*reg3; reg3=reg3*reg2;
    T reg7=reg2*reg1; T reg8=reg4*reg3; T reg9=reg4*reg6; reg3=reg3*reg2; reg1=reg4*reg1;
    T reg10=reg4*reg5; reg6=reg6*reg2; reg5=reg2*reg5; reg8=reg9+reg8; reg3=reg3-reg9;
    reg7=reg7-reg10; reg1=reg10+reg1; reg1=reg4*reg1; T reg11=reg3*reg2; reg7=reg2*reg7;
    T reg12=reg10+reg5; T reg13=reg4*reg8; reg6=reg9+reg6; reg12=reg4*reg12; reg13=reg11-reg13;
    reg9=reg4*reg6; reg1=reg7-reg1; reg12=reg1-reg12; reg9=reg13-reg9; reg12=reg12/reg9;
    reg8=reg8/reg9; reg3=reg3/reg9; reg1=reg3*reg12; reg7=reg8*reg12; reg11=reg3*reg1;
    reg13=reg7*reg8; T reg14=(*f.m).alpha*reg3; T reg15=(*f.m).alpha*reg8; reg9=reg6/reg9; reg13=reg11-reg13;
    reg15=reg14+reg15; reg9=(*f.m).alpha*reg9; reg7=reg7/reg13; reg13=reg1/reg13; reg9=reg15+reg9;
    reg1=reg4*reg0; reg6=reg2*reg0; reg13=reg13*reg9; reg9=reg7*reg9; reg7=reg6*reg2;
    reg11=reg4*reg1; reg9=reg13-reg9; reg11=reg7-reg11; reg7=1-(*f.m).resolution; reg9=reg7*reg9;
    reg6=reg6/reg11; reg1=reg1/reg11; reg13=(*f.m).resolution*(*f.m).alpha; reg8=reg7*reg8; reg14=(*f.m).resolution*reg6;
    reg3=reg3*reg7; reg15=(*f.m).resolution*reg1; T reg16=elem.pos(1)[0]-elem.pos(0)[0]; T reg17=elem.pos(1)[1]-elem.pos(0)[1]; reg9=reg13+reg9;
    reg13=elem.pos(2)[1]-elem.pos(0)[1]; T reg18=elem.pos(2)[0]-elem.pos(0)[0]; reg8=reg15+reg8; reg9=(*f.m).deltaT*reg9; reg15=reg16*reg13;
    T reg19=reg17*reg18; reg14=reg3+reg14; reg19=reg15-reg19; reg3=reg14*reg9; reg15=reg8*reg9;
    reg13=reg13/reg19; T reg20=1-var_inter[0]; reg16=reg16/reg19; reg18=reg18/reg19; reg17=reg17/reg19;
    T reg21=reg3+reg15; T reg22=var_inter[1]*(*f.m).f_vol[0]; T reg23=reg17-reg13; T reg24=reg18-reg16; T reg25=var_inter[0]*(*f.m).f_vol[1];
    T reg26=reg18*reg21; T reg27=reg17*reg21; reg20=reg20-var_inter[1]; T reg28=reg21*reg24; T reg29=reg20*(*f.m).f_vol[0];
    T reg30=reg21*reg13; T reg31=reg20*(*f.m).f_vol[1]; T reg32=reg27-reg22; T reg33=var_inter[0]*(*f.m).f_vol[0]; T reg34=reg16*reg21;
    T reg35=reg26-reg25; T reg36=reg21*reg23; T reg37=var_inter[1]*(*f.m).f_vol[1]; reg32=reg19*reg32; T reg38=reg30+reg33;
    T reg39=reg34+reg37; T reg40=reg36+reg29; T reg41=reg28+reg31; reg35=reg19*reg35; T reg42=reg19*reg41;
    reg35=ponderation*reg35; T reg43=reg19*reg40; T reg44=reg19*reg38; reg32=ponderation*reg32; T reg45=reg19*reg39;
    sollicitation[indices[1]+1]+=-reg35; reg35=ponderation*reg43; sollicitation[indices[0]+0]+=reg35; sollicitation[indices[2]+0]+=-reg32; reg32=ponderation*reg42;
    sollicitation[indices[0]+1]+=reg32; T reg46=ponderation*reg45; sollicitation[indices[2]+1]+=reg46; T reg47=ponderation*reg44; sollicitation[indices[1]+0]+=reg47;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=1.0/(*f.m).elastic_modulus; T reg3=reg1*reg0;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg1*reg2; reg1=reg4*reg1; T reg6=reg4*reg3; reg3=reg3*reg2;
    T reg7=reg4*reg1; T reg8=reg4*reg5; T reg9=reg4*reg3; T reg10=reg4*reg6; reg3=reg3*reg2;
    reg5=reg2*reg5; reg1=reg2*reg1; reg5=reg5-reg7; reg3=reg3-reg10; reg9=reg10+reg9;
    reg8=reg7+reg8; reg6=reg6*reg2; T reg11=reg3*reg2; T reg12=reg7+reg1; reg8=reg4*reg8;
    reg5=reg2*reg5; T reg13=reg4*reg9; reg6=reg10+reg6; reg8=reg5-reg8; reg5=reg4*reg6;
    reg12=reg4*reg12; reg13=reg11-reg13; reg12=reg8-reg12; reg5=reg13-reg5; reg3=reg3/reg5;
    reg9=reg9/reg5; reg12=reg12/reg5; reg8=reg9*reg12; reg10=reg3*reg12; reg11=(*f.m).alpha*reg3;
    reg5=reg6/reg5; reg6=(*f.m).alpha*reg9; reg13=reg3*reg10; T reg14=reg8*reg9; reg6=reg11+reg6;
    reg5=(*f.m).alpha*reg5; reg14=reg13-reg14; reg11=elem.pos(1)[0]-elem.pos(0)[0]; reg13=elem.pos(1)[1]-elem.pos(0)[1]; T reg15=elem.pos(2)[0]-elem.pos(0)[0];
    T reg16=elem.pos(2)[1]-elem.pos(0)[1]; reg10=reg10/reg14; reg14=reg8/reg14; reg5=reg6+reg5; reg6=reg11*reg16;
    reg8=reg13*reg15; T reg17=reg4*reg0; reg10=reg10*reg5; reg5=reg14*reg5; reg8=reg6-reg8;
    reg6=reg2*reg0; reg14=1-(*f.m).resolution; reg5=reg10-reg5; reg10=vectors[0][indices[1]+0]-vectors[0][indices[0]+0]; T reg18=vectors[0][indices[2]+0]-vectors[0][indices[0]+0];
    T reg19=vectors[0][indices[1]+1]-vectors[0][indices[0]+1]; T reg20=vectors[0][indices[2]+1]-vectors[0][indices[0]+1]; T reg21=pow(reg4,2); T reg22=reg4*reg17; T reg23=pow(reg2,2);
    T reg24=reg6*reg2; reg13=reg13/reg8; reg15=reg15/reg8; reg11=reg11/reg8; reg16=reg16/reg8;
    T reg25=reg19*reg16; T reg26=reg15*reg10; T reg27=reg11*reg18; T reg28=reg13*reg20; T reg29=(*f.m).resolution*(*f.m).alpha;
    reg22=reg24-reg22; reg21=reg23-reg21; reg5=reg14*reg5; reg19=reg15*reg19; reg20=reg11*reg20;
    reg10=reg16*reg10; reg18=reg13*reg18; reg5=reg29+reg5; reg6=reg6/reg22; reg21=reg21/reg22;
    reg22=reg17/reg22; reg28=reg25-reg28; reg26=reg27-reg26; reg18=reg10-reg18; reg5=(*f.m).deltaT*reg5;
    reg19=reg20-reg19; reg28=reg26+reg28; reg12=reg14*reg12; reg9=reg14*reg9; reg10=(*f.m).resolution*reg6;
    reg17=(*f.m).resolution*reg21; reg14=reg3*reg14; reg3=(*f.m).resolution*reg22; reg12=reg17+reg12; reg9=reg3+reg9;
    reg28=0.5*reg28; reg10=reg14+reg10; reg19=reg19-reg5; reg18=reg18-reg5; reg3=reg19*reg10;
    reg14=reg18*reg9; reg17=reg13-reg16; reg18=reg18*reg10; reg20=reg15-reg11; reg28=reg28*reg12;
    reg19=reg19*reg9; reg28=2*reg28; reg23=1-var_inter[0]; reg3=reg14+reg3; reg14=0.5*reg13;
    reg19=reg18+reg19; reg18=0.5*reg15; reg24=0.5*reg16; reg25=0.5*reg17; reg26=0.5*reg11;
    reg27=0.5*reg20; reg29=reg19*reg16; T reg30=reg28*reg25; T reg31=reg28*reg18; T reg32=reg15*reg3;
    T reg33=reg3*reg20; T reg34=reg19*reg17; T reg35=reg28*reg27; reg23=reg23-var_inter[1]; T reg36=reg28*reg24;
    T reg37=reg11*reg3; T reg38=reg28*reg14; T reg39=reg13*reg19; T reg40=reg28*reg26; T reg41=var_inter[1]*(*f.m).f_vol[0];
    reg37=reg37-reg38; T reg42=var_inter[0]*(*f.m).f_vol[0]; reg29=reg29-reg31; T reg43=var_inter[1]*(*f.m).f_vol[1]; T reg44=reg23*(*f.m).f_vol[1];
    reg30=reg33+reg30; reg40=reg40-reg39; reg33=var_inter[0]*(*f.m).f_vol[1]; reg36=reg36-reg32; T reg45=reg23*(*f.m).f_vol[0];
    reg35=reg34+reg35; reg36=reg36-reg33; reg29=reg29-reg42; reg37=reg37-reg43; reg30=reg30-reg44;
    reg35=reg35-reg45; reg40=reg40-reg41; reg40=reg8*reg40; reg37=reg8*reg37; reg30=reg8*reg30;
    reg36=reg8*reg36; reg29=reg8*reg29; reg35=reg8*reg35; sollicitation[indices[2]+0]+=ponderation*reg40; sollicitation[indices[1]+0]+=ponderation*reg29;
    sollicitation[indices[1]+1]+=ponderation*reg36; sollicitation[indices[0]+1]+=ponderation*reg30; sollicitation[indices[2]+1]+=ponderation*reg37; sollicitation[indices[0]+0]+=ponderation*reg35;
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
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg3*reg1; reg1=reg4*reg1; T reg6=reg3*reg2; reg2=reg4*reg2;
    T reg7=reg4*reg1; T reg8=reg3*reg5; T reg9=reg3*reg2; T reg10=reg3*reg6; reg1=reg3*reg1;
    reg2=reg4*reg2; reg6=reg4*reg6; reg9=reg10+reg9; reg2=reg2-reg10; reg5=reg4*reg5;
    reg7=reg7-reg8; reg1=reg8+reg1; T reg11=0.5*elem.pos(1)[1]; T reg12=0.5*elem.pos(0)[1]; reg6=reg10+reg6;
    reg10=0.5*elem.pos(1)[0]; T reg13=0.5*elem.pos(0)[0]; reg7=reg4*reg7; T reg14=reg8+reg5; reg1=reg3*reg1;
    T reg15=reg4*reg2; T reg16=reg3*reg9; T reg17=reg10-reg13; reg10=reg13+reg10; reg13=reg11-reg12;
    T reg18=0.5*elem.pos(2)[0]; reg11=reg12+reg11; reg12=0.5*elem.pos(2)[1]; T reg19=reg3*reg6; reg1=reg7-reg1;
    reg14=reg3*reg14; reg16=reg15-reg16; reg10=reg18-reg10; reg13=reg12+reg13; reg14=reg1-reg14;
    reg1=0.5*elem.pos(3)[0]; reg17=reg18+reg17; reg7=0.5*elem.pos(3)[1]; reg19=reg16-reg19; reg11=reg12-reg11;
    reg14=reg14/reg19; reg12=0.5*vectors[0][indices[0]+1]; reg2=reg2/reg19; reg9=reg9/reg19; reg15=0.5*vectors[0][indices[1]+1];
    reg16=0.5*vectors[0][indices[0]+0]; reg18=0.5*vectors[0][indices[1]+0]; reg10=reg1+reg10; reg11=reg11+reg7; reg7=reg13-reg7;
    reg1=reg17-reg1; reg13=reg9*reg14; reg17=reg2*reg14; T reg20=0.21132486540518713447*elem.pos(1)[1]; T reg21=reg15-reg12;
    T reg22=0.21132486540518713447*elem.pos(0)[0]; T reg23=0.21132486540518713447*elem.pos(1)[0]; T reg24=0.21132486540518713447*elem.pos(0)[1]; T reg25=0.78867513459481286553*elem.pos(1)[1]; T reg26=reg11*reg1;
    T reg27=reg7*reg10; T reg28=0.5*vectors[0][indices[2]+0]; T reg29=reg18-reg16; reg16=reg18+reg16; reg18=0.78867513459481286553*elem.pos(1)[0];
    T reg30=0.78867513459481286553*elem.pos(0)[0]; T reg31=0.5*vectors[0][indices[2]+1]; reg15=reg12+reg15; reg12=0.78867513459481286553*elem.pos(0)[1]; T reg32=reg9*(*f.m).alpha;
    T reg33=reg2*(*f.m).alpha; T reg34=reg9*reg13; T reg35=reg2*reg17; T reg36=0.78867513459481286553*elem.pos(2)[1]; T reg37=reg30+reg23;
    T reg38=0.21132486540518713447*elem.pos(2)[0]; T reg39=reg12+reg20; reg19=reg6/reg19; reg6=0.21132486540518713447*elem.pos(2)[1]; reg27=reg26-reg27;
    reg29=reg28+reg29; reg26=0.5*vectors[0][indices[3]+0]; reg16=reg28-reg16; reg28=reg22+reg18; T reg40=0.5*vectors[0][indices[3]+1];
    reg15=reg31-reg15; reg21=reg31+reg21; reg20=reg20-reg24; reg24=reg24+reg25; reg22=reg23-reg22;
    reg23=0.78867513459481286553*elem.pos(2)[0]; reg30=reg18-reg30; reg10=reg10/reg27; reg15=reg40+reg15; reg12=reg25-reg12;
    reg22=reg23+reg22; reg1=reg1/reg27; reg16=reg26+reg16; reg7=reg7/reg27; reg28=reg23-reg28;
    reg26=reg29-reg26; reg18=0.21132486540518713447*elem.pos(3)[1]; reg24=reg36-reg24; reg27=reg11/reg27; reg11=reg3*reg0;
    reg23=reg4*reg0; reg25=0.21132486540518713447*elem.pos(3)[0]; reg29=0.78867513459481286553*elem.pos(3)[0]; reg40=reg21-reg40; reg19=(*f.m).alpha*reg19;
    reg32=reg33+reg32; reg21=0.78867513459481286553*elem.pos(3)[1]; reg39=reg6-reg39; reg34=reg35-reg34; reg37=reg38-reg37;
    reg20=reg36+reg20; reg39=reg39+reg21; reg28=reg28+reg25; reg31=0.21132486540518713447*PNODE(1).dep[0]; reg33=0.21132486540518713447*PNODE(0).dep[0];
    reg37=reg37+reg29; reg12=reg6+reg12; reg6=0.78867513459481286553*PNODE(0).dep[1]; reg35=0.78867513459481286553*PNODE(1).dep[0]; reg30=reg38+reg30;
    reg36=0.78867513459481286553*PNODE(0).dep[0]; reg24=reg24+reg18; reg29=reg22-reg29; reg19=reg32+reg19; reg13=reg13/reg34;
    reg34=reg17/reg34; reg17=reg10*reg40; reg22=reg3*reg11; reg32=reg1*reg15; reg38=0.21132486540518713447*PNODE(1).dep[1];
    T reg41=reg27*reg26; T reg42=0.78867513459481286553*PNODE(1).dep[1]; T reg43=0.21132486540518713447*PNODE(0).dep[1]; T reg44=reg4*reg23; reg21=reg20-reg21;
    reg20=reg7*reg16; T reg45=0.21132486540518713447*PNODE(2).dep[0]; T reg46=reg38+reg6; T reg47=reg31+reg36; T reg48=0.21132486540518713447*PNODE(2).dep[1];
    T reg49=reg37*reg21; reg18=reg12-reg18; reg12=(*f.m).deltaT*(*f.m).alpha; reg22=reg44-reg22; reg34=reg34*reg19;
    reg19=reg13*reg19; reg13=reg39*reg29; reg44=reg24*reg29; reg25=reg30-reg25; reg30=reg21*reg28;
    T reg50=0.78867513459481286553*PNODE(2).dep[1]; T reg51=reg43+reg42; T reg52=reg33+reg35; reg17=reg32-reg17; elem.epsilon[0][1]=reg17;
    reg20=reg41-reg20; elem.epsilon[0][0]=reg20; reg32=0.78867513459481286553*PNODE(2).dep[0]; reg33=reg31-reg33; reg43=reg38-reg43;
    reg31=1-(*f.m).resolution; reg49=reg13-reg49; reg13=0.78867513459481286553*PNODE(3).dep[0]; reg52=reg32-reg52; reg46=reg48-reg46;
    reg23=reg23/reg22; reg11=reg11/reg22; reg19=reg34-reg19; reg34=0.21132486540518713447*PNODE(3).dep[0]; reg47=reg45-reg47;
    reg36=reg35-reg36; reg6=reg42-reg6; reg35=reg25*reg24; reg33=reg32+reg33; reg32=reg18*reg28;
    reg38=reg20-reg12; reg30=reg44-reg30; reg41=0.21132486540518713447*PNODE(3).dep[1]; reg43=reg43+reg50; reg42=reg17-reg12;
    reg44=0.78867513459481286553*PNODE(3).dep[1]; reg51=reg50-reg51; reg43=reg43-reg44; reg33=reg33-reg13; reg50=reg29/reg30;
    T reg53=reg24/reg30; reg32=reg35-reg32; reg35=reg28/reg30; T reg54=(*f.m).alpha*(*f.m).resolution; reg19=reg31*reg19;
    reg6=reg48+reg6; reg48=reg42*reg23; T reg55=reg38*reg11; reg36=reg45+reg36; reg29=reg29/reg49;
    reg51=reg51+reg41; reg52=reg34+reg52; reg42=reg42*reg11; reg46=reg44+reg46; reg38=reg38*reg23;
    reg44=reg37/reg49; reg45=pow(reg4,2); T reg56=pow(reg3,2); T reg57=reg39/reg49; T reg58=reg21/reg49;
    T reg59=reg39*reg25; T reg60=reg37*reg18; reg21=reg21/reg30; reg47=reg13+reg47; reg56=reg45-reg56;
    reg28=reg28/reg32; reg13=reg25/reg32; reg41=reg6-reg41; reg24=reg24/reg32; reg6=reg51*reg50;
    reg60=reg59-reg60; reg34=reg36-reg34; reg36=reg35*reg43; reg45=reg18/reg32; reg59=reg53*reg33;
    T reg61=reg33*reg57; T reg62=reg58*reg47; T reg63=reg43*reg44; reg48=reg55+reg48; reg55=reg29*reg46;
    T reg64=reg21*reg52; reg42=reg38+reg42; reg19=reg54+reg19; reg29=reg29*reg47; reg44=reg33*reg44;
    reg50=reg52*reg50; reg57=reg43*reg57; reg58=reg46*reg58; reg33=reg35*reg33; reg21=reg51*reg21;
    reg53=reg43*reg53; reg9=reg9*reg31; reg35=(*f.m).resolution*reg11; reg38=reg52*reg13; reg43=reg28*reg34;
    reg54=reg41*reg24; reg36=reg6-reg36; reg6=reg51*reg45; reg22=reg56/reg22; reg63=reg55-reg63;
    reg33=reg50-reg33; reg21=reg53-reg21; reg64=reg59-reg64; reg18=reg18/reg60; reg39=reg39/reg60;
    reg25=reg25/reg60; reg37=reg37/reg60; reg50=reg4*reg48; reg53=reg4*reg42; reg28=reg28*reg41;
    reg62=reg61-reg62; reg55=(*f.m).resolution*reg23; reg24=reg24*reg34; reg45=reg52*reg45; reg48=reg3*reg48;
    reg19=(*f.m).deltaT*reg19; reg42=reg3*reg42; reg2=reg2*reg31; reg58=reg57-reg58; reg44=reg29-reg44;
    reg13=reg51*reg13; reg29=reg34*reg39; reg58=reg44+reg58; reg44=reg41*reg37; reg51=reg47*reg18;
    reg21=reg33+reg21; reg9=reg35+reg9; reg53=reg53-reg48; reg36=reg36-reg19; reg18=reg46*reg18;
    reg39=reg41*reg39; reg47=reg47*reg25; reg37=reg34*reg37; reg31=reg14*reg31; reg55=reg2+reg55;
    reg2=(*f.m).resolution*reg22; reg50=reg50-reg42; reg45=reg24-reg45; reg43=reg38-reg43; reg6=reg54-reg6;
    reg63=reg63-reg19; reg62=reg62-reg19; reg28=reg13-reg28; reg64=reg64-reg19; reg25=reg46*reg25;
    reg15=reg7*reg15; reg40=reg27*reg40; reg53=reg53+reg12; reg10=reg26*reg10; reg1=reg16*reg1;
    reg18=reg39-reg18; reg42=reg48+reg42; reg31=reg2+reg31; reg50=reg50+reg12; reg6=reg43+reg6;
    reg45=reg45-reg19; reg28=reg28-reg19; reg58=0.5*reg58; reg2=reg63*reg9; reg7=reg62*reg55;
    reg44=reg25-reg44; reg51=reg29-reg51; reg13=reg63*reg55; reg14=reg62*reg9; reg16=reg64*reg9;
    reg24=reg36*reg55; reg37=reg47-reg37; reg21=0.5*reg21; reg25=reg36*reg9; reg26=reg64*reg55;
    reg6=0.5*reg6; reg42=reg12-reg42; reg27=reg53+reg50; reg29=reg28*reg9; reg33=reg45*reg55;
    reg34=reg28*reg55; reg35=reg45*reg9; reg38=reg58*reg31; reg39=reg21*reg31; reg2=reg7+reg2;
    reg15=reg40-reg15; reg24=reg16+reg24; reg25=reg26+reg25; reg51=reg51-reg19; reg10=reg1-reg10;
    reg44=reg44-reg19; reg13=reg14+reg13; reg18=reg37+reg18; reg24=reg36*reg24; reg25=reg64*reg25;
    reg38=2*reg38; reg34=reg35+reg34; reg13=reg63*reg13; reg2=reg62*reg2; reg29=reg33+reg29;
    reg39=2*reg39; reg27=reg27+reg42; reg15=reg10+reg15; reg18=0.5*reg18; reg1=reg51*reg9;
    reg7=reg6*reg31; reg10=reg44*reg55; reg14=reg51*reg55; reg16=reg44*reg9; reg38=reg58*reg38;
    reg2=reg13+reg2; reg10=reg1+reg10; reg15=0.5*reg15; elem.epsilon[0][2]=reg15; reg34=reg28*reg34;
    reg1=reg18*reg31; reg25=reg24+reg25; reg7=2*reg7; reg29=reg45*reg29; reg39=reg21*reg39;
    reg27=reg27/3; reg16=reg14+reg16; reg1=2*reg1; reg39=reg25+reg39; reg10=reg44*reg10;
    reg16=reg51*reg16; reg29=reg34+reg29; reg38=reg2+reg38; reg53=reg53-reg27; reg2=reg15*reg22;
    reg50=reg50-reg27; reg7=reg6*reg7; reg1=reg18*reg1; reg27=reg42-reg27; reg2=reg0*reg2;
    reg38=reg49*reg38; reg30=reg39*reg30; reg50=pow(reg50,2); reg7=reg29+reg7; reg53=pow(reg53,2);
    reg16=reg10+reg16; reg30=0.25*reg30; reg50=reg53+reg50; reg1=reg16+reg1; reg6=2*reg2;
    reg7=reg32*reg7; reg38=0.25*reg38; reg27=pow(reg27,2); reg6=reg2*reg6; reg7=0.25*reg7;
    reg1=reg60*reg1; reg38=reg30+reg38; reg27=reg50+reg27; reg6=reg27+reg6; reg17=reg17-reg19;
    reg20=reg20-reg19; reg7=reg38+reg7; reg1=0.25*reg1; reg1=reg7+reg1; reg2=reg17*reg55;
    reg7=reg20*reg9; reg17=reg17*reg9; reg20=reg20*reg55; reg6=1.5*reg6; elem.sigma[0][2]=reg15*reg31;
    elem.sigma_von_mises=pow(reg6,0.5); elem.sigma[0][0]=reg20+reg17; elem.sigma[0][1]=reg7+reg2; elem.ener=reg1/2;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg3=reg0*reg1;
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg2*reg1; T reg6=reg4*reg3; reg3=reg2*reg3; reg1=reg4*reg1;
    T reg7=reg4*reg6; T reg8=reg2*reg3; T reg9=reg2*reg1; reg6=reg2*reg6; T reg10=reg2*reg5;
    reg1=reg4*reg1; reg3=reg4*reg3; reg1=reg1-reg10; reg6=reg8+reg6; reg5=reg4*reg5;
    reg9=reg10+reg9; reg7=reg7-reg8; reg1=reg4*reg1; reg9=reg2*reg9; T reg11=reg10+reg5;
    T reg12=reg4*reg7; reg3=reg8+reg3; reg8=reg2*reg6; reg11=reg2*reg11; reg8=reg12-reg8;
    reg9=reg1-reg9; reg1=reg2*reg3; reg11=reg9-reg11; reg9=1-var_inter[1]; reg12=1-var_inter[0];
    reg1=reg8-reg1; reg6=reg6/reg1; reg7=reg7/reg1; reg8=elem.pos(1)[1]*var_inter[0]; T reg13=reg12*elem.pos(0)[1];
    T reg14=reg9*elem.pos(0)[0]; T reg15=reg9*elem.pos(1)[0]; T reg16=reg12*elem.pos(0)[0]; T reg17=elem.pos(1)[0]*var_inter[0]; reg11=reg11/reg1;
    T reg18=reg9*elem.pos(0)[1]; T reg19=reg9*elem.pos(1)[1]; reg15=reg15-reg14; T reg20=var_inter[1]*elem.pos(2)[0]; reg19=reg19-reg18;
    T reg21=var_inter[1]*elem.pos(2)[1]; T reg22=reg6*reg11; T reg23=reg7*reg11; T reg24=elem.pos(2)[1]*var_inter[0]; T reg25=elem.pos(2)[0]*var_inter[0];
    T reg26=reg17+reg16; T reg27=reg13+reg8; reg24=reg24-reg27; T reg28=elem.pos(3)[1]*reg12; reg1=reg3/reg1;
    reg3=reg7*(*f.m).alpha; T reg29=reg6*(*f.m).alpha; reg25=reg25-reg26; T reg30=elem.pos(3)[0]*reg12; T reg31=var_inter[1]*elem.pos(3)[1];
    T reg32=reg7*reg23; T reg33=reg6*reg22; reg21=reg19+reg21; reg19=var_inter[1]*elem.pos(3)[0]; reg15=reg20+reg15;
    reg33=reg32-reg33; reg28=reg24+reg28; reg15=reg15-reg19; reg21=reg21-reg31; reg30=reg25+reg30;
    reg29=reg3+reg29; reg1=(*f.m).alpha*reg1; reg23=reg23/reg33; reg3=reg4*reg0; reg20=reg2*reg0;
    reg33=reg22/reg33; reg1=reg29+reg1; reg22=reg21*reg30; reg24=reg28*reg15; reg25=reg4*reg3;
    reg22=reg24-reg22; reg24=reg2*reg20; reg29=pow(reg4,2); reg32=pow(reg2,2); reg23=reg23*reg1;
    reg1=reg33*reg1; reg15=reg15/reg22; reg32=reg29-reg32; reg29=1-(*f.m).resolution; reg21=reg21/reg22;
    reg1=reg23-reg1; reg24=reg25-reg24; reg30=reg30/reg22; reg28=reg28/reg22; reg20=reg20/reg24;
    reg23=var_inter[0]*reg21; reg25=var_inter[0]*reg15; reg33=(*f.m).alpha*(*f.m).resolution; reg1=reg29*reg1; reg32=reg32/reg24;
    T reg34=var_inter[1]*reg28; T reg35=reg9*reg28; reg24=reg3/reg24; reg3=reg12*reg21; T reg36=reg12*reg15;
    T reg37=reg9*reg30; T reg38=var_inter[1]*reg30; T reg39=reg3+reg34; reg11=reg11*reg29; reg6=reg6*reg29;
    T reg40=(*f.m).resolution*reg32; T reg41=(*f.m).resolution*reg20; T reg42=reg36+reg38; reg29=reg7*reg29; reg7=(*f.m).resolution*reg24;
    T reg43=reg35+reg23; reg1=reg33+reg1; reg33=reg37+reg25; T reg44=0.5*reg42; T reg45=0.5*reg39;
    T reg46=0.5*reg33; T reg47=reg37-reg36; T reg48=reg3-reg35; reg1=(*f.m).deltaT*reg1; T reg49=reg25-reg38;
    reg11=reg40+reg11; reg6=reg41+reg6; reg40=0.5*reg43; reg7=reg29+reg7; reg29=reg34-reg23;
    reg41=0.5*reg47; T reg50=reg6*reg1; T reg51=reg7*reg1; T reg52=reg11*reg40; T reg53=reg11*reg44;
    T reg54=reg11*reg46; T reg55=reg11*reg45; T reg56=0.5*reg29; T reg57=0.5*reg49; T reg58=0.5*reg48;
    T reg59=reg6*reg42; T reg60=reg11*reg56; T reg61=reg51+reg50; T reg62=reg7*reg33; T reg63=reg7*reg39;
    T reg64=reg6*reg43; T reg65=reg9*var_inter[0]; T reg66=2*reg55; T reg67=reg7*reg42; reg53=2*reg53;
    T reg68=reg6*reg39; T reg69=var_inter[1]*reg12; T reg70=reg11*reg41; T reg71=reg11*reg57; reg52=2*reg52;
    T reg72=reg6*reg33; T reg73=2*reg54; T reg74=reg7*reg43; T reg75=reg11*reg58; T reg76=reg73*reg44;
    T reg77=reg64*reg33; T reg78=reg7*reg47; T reg79=reg53*reg46; T reg80=reg63*reg43; T reg81=reg73*reg40;
    T reg82=reg66*reg44; T reg83=reg61*reg33; T reg84=reg74*reg39; T reg85=reg62*reg42; T reg86=reg52*reg45;
    T reg87=reg66*reg40; T reg88=reg67*reg33; T reg89=reg6*reg29; T reg90=reg72*reg43; T reg91=reg53*reg45;
    T reg92=reg68*reg42; T reg93=reg52*reg46; T reg94=reg7*reg49; reg75=2*reg75; T reg95=reg9*reg12;
    T reg96=reg6*reg48; T reg97=reg59*reg39; T reg98=var_inter[1]*var_inter[0]; T reg99=reg6*reg49; T reg100=reg61*reg39;
    reg60=2*reg60; T reg101=reg65*(*f.m).f_vol[1]; reg70=2*reg70; T reg102=reg7*reg48; T reg103=reg6*reg47;
    T reg104=reg69*(*f.m).f_vol[0]; reg71=2*reg71; T reg105=reg7*reg29; T reg106=reg66*reg57; T reg107=reg75*reg57;
    T reg108=reg96*reg49; T reg109=reg74*reg29; T reg110=reg59*reg29; T reg111=reg73*reg57; T reg112=reg72*reg29;
    T reg113=reg52*reg57; T reg114=reg105*reg29; T reg115=reg71*reg57; T reg116=reg99*reg29; T reg117=reg60*reg57;
    T reg118=reg53*reg57; T reg119=reg63*reg29; T reg120=reg70*reg40; T reg121=reg96*reg33; T reg122=reg75*reg40;
    T reg123=reg78*reg33; T reg124=reg64*reg42; reg77=reg81+reg77; T reg125=reg52*reg40; T reg126=reg62*reg33;
    T reg127=reg71*reg40; T reg128=reg89*reg33; T reg129=reg60*reg40; T reg130=reg94*reg33; T reg131=reg53*reg40;
    T reg132=reg68*reg33; reg88=reg87+reg88; T reg133=reg29*reg102; T reg134=reg70*reg57; T reg135=reg103*reg29;
    reg91=reg92+reg91; T reg136=reg67*reg42; T reg137=reg66*reg45; T reg138=reg95*(*f.m).f_vol[0]; T reg139=reg65*(*f.m).f_vol[0];
    T reg140=reg61*reg48; T reg141=reg61*reg47; T reg142=reg61*reg43; T reg143=reg83-reg101; T reg144=reg61*reg29;
    T reg145=reg61*reg49; T reg146=reg52*reg44; T reg147=reg72*reg39; T reg148=reg71*reg44; T reg149=reg105*reg39;
    T reg150=reg103*reg39; T reg151=reg48*reg102; T reg152=reg98*(*f.m).f_vol[0]; T reg153=reg98*(*f.m).f_vol[1]; T reg154=reg69*(*f.m).f_vol[1];
    T reg155=reg70*reg56; T reg156=reg78*reg49; T reg157=reg75*reg56; T reg158=reg64*reg49; T reg159=reg73*reg56;
    T reg160=reg62*reg49; T reg161=reg52*reg56; T reg162=reg95*(*f.m).f_vol[1]; T reg163=reg89*reg49; T reg164=reg71*reg56;
    T reg165=reg94*reg49; T reg166=reg60*reg56; T reg167=reg68*reg49; T reg168=reg53*reg56; T reg169=reg67*reg49;
    T reg170=reg66*reg56; T reg171=reg89*reg42; T reg172=reg71*reg45; T reg173=reg94*reg42; T reg174=reg60*reg45;
    T reg175=reg60*reg41; T reg176=reg99*reg48; reg93=reg90+reg93; T reg177=reg60*reg46; reg89=reg89*reg47;
    T reg178=reg71*reg41; T reg179=reg105*reg48; T reg180=reg71*reg58; T reg181=reg52*reg41; T reg182=reg72*reg48;
    T reg183=reg73*reg41; reg105=reg105*reg43; reg71=reg71*reg46; T reg184=reg74*reg48; T reg185=reg43*reg102;
    T reg186=reg75*reg41; T reg187=reg103*reg48; reg103=reg103*reg43; T reg188=reg75*reg46; T reg189=reg73*reg58;
    reg64=reg64*reg47; T reg190=reg75*reg58; T reg191=reg62*reg47; T reg192=reg78*reg47; reg52=reg52*reg58;
    reg74=reg74*reg43; T reg193=reg70*reg58; T reg194=reg96*reg47; T reg195=reg73*reg46; T reg196=reg66*reg41;
    T reg197=reg59*reg48; T reg198=reg53*reg41; T reg199=reg63*reg48; T reg200=reg70*reg46; T reg201=reg100-reg104;
    T reg202=reg66*reg46; reg59=reg59*reg43; T reg203=reg53*reg58; T reg204=reg61*reg42; reg67=reg67*reg47;
    T reg205=reg75*reg44; T reg206=reg66*reg58; T reg207=reg70*reg44; reg75=reg75*reg45; T reg208=reg73*reg45;
    reg79=reg80+reg79; T reg209=reg68*reg47; reg86=reg85+reg86; T reg210=reg60*reg44; reg53=reg53*reg44;
    T reg211=reg70*reg41; reg102=reg39*reg102; reg78=reg78*reg42; T reg212=reg99*reg43; reg99=reg99*reg39;
    reg94=reg94*reg47; reg70=reg70*reg45; reg96=reg96*reg42; T reg213=reg63*reg39; reg60=reg60*reg58;
    reg97=reg82+reg97; reg84=reg84+reg76; reg164=reg163+reg164; reg203=reg203-reg209; reg174=reg173-reg174;
    reg64=reg64-reg189; reg163=reg141+reg162; reg173=reg140+reg138; reg172=reg171-reg172; reg60=reg94+reg60;
    reg168=reg168-reg167; reg169=reg169-reg170; reg52=reg52-reg191; reg180=reg89+reg180; reg136=reg136+reg137;
    reg166=reg165+reg166; reg89=reg22*reg91; reg94=reg145+reg153; reg201=reg22*reg201; reg165=reg204+reg154;
    reg102=reg207-reg102; reg75=reg78-reg75; reg124=reg208+reg124; reg78=reg22*reg86; reg53=reg53+reg213;
    reg171=reg22*reg84; reg207=reg22*reg97; reg70=reg96-reg70; reg150=reg205-reg150; reg151=reg211+reg151;
    reg186=reg187+reg186; reg184=reg184-reg183; reg181=reg181-reg182; reg178=reg179+reg178; reg175=reg176+reg175;
    reg198=reg198-reg199; reg197=reg197-reg196; reg193=reg194+reg193; reg149=reg148-reg149; reg146=reg146+reg147;
    reg190=reg192+reg190; reg96=reg144+reg152; reg143=reg22*reg143; reg148=reg142+reg139; reg134=reg133+reg134;
    reg157=reg156+reg157; reg99=reg210-reg99; reg133=reg22*reg93; reg107=reg135+reg107; reg155=reg108+reg155;
    reg200=reg185-reg200; reg109=reg109-reg111; reg123=reg122-reg123; reg74=reg74+reg195; reg113=reg113-reg112;
    reg110=reg110-reg106; reg125=reg125+reg126; reg177=reg212-reg177; reg108=reg22*reg79; reg115=reg114+reg115;
    reg118=reg118-reg119; reg117=reg116+reg117; reg114=reg22*reg77; reg188=reg103-reg188; reg71=reg105-reg71;
    reg121=reg120-reg121; reg161=reg161-reg160; reg158=reg158-reg159; reg103=reg22*reg88; reg131=reg131+reg132;
    reg130=reg129-reg130; reg67=reg67-reg206; reg59=reg59+reg202; reg128=reg127-reg128; reg113=reg22*reg113;
    reg125=reg22*reg125; reg193=reg22*reg193; reg131=reg22*reg131; reg149=reg22*reg149; reg130=reg22*reg130;
    reg105=reg22*reg96; reg71=reg22*reg71; reg116=ponderation*reg78; reg177=reg22*reg177; reg190=reg22*reg190;
    reg146=reg22*reg146; reg151=reg22*reg151; reg115=reg22*reg115; reg188=reg22*reg188; reg181=reg22*reg181;
    reg120=ponderation*reg133; reg134=reg22*reg134; reg122=ponderation*reg207; reg178=reg22*reg178; reg128=reg22*reg128;
    reg184=reg22*reg184; reg107=reg22*reg107; reg99=reg22*reg99; reg175=reg22*reg175; reg53=reg22*reg53;
    reg74=reg22*reg74; reg127=ponderation*reg103; reg198=reg22*reg198; reg109=reg22*reg109; reg129=ponderation*reg171;
    reg197=reg22*reg197; reg70=reg22*reg70; reg186=reg22*reg186; reg135=ponderation*reg89; reg155=reg22*reg155;
    reg150=reg22*reg150; reg156=ponderation*reg108; reg174=reg22*reg174; reg180=reg22*reg180; reg176=reg22*reg165;
    reg172=reg22*reg172; reg157=reg22*reg157; reg121=reg22*reg121; reg102=reg22*reg102; reg169=reg22*reg169;
    reg67=reg22*reg67; reg60=reg22*reg60; reg201=ponderation*reg201; reg168=reg22*reg168; reg158=reg22*reg158;
    reg166=reg22*reg166; reg161=reg22*reg161; reg179=reg22*reg94; reg164=reg22*reg164; reg59=reg22*reg59;
    reg203=reg22*reg203; reg185=reg22*reg163; reg64=reg22*reg64; reg118=reg22*reg118; reg187=reg22*reg148;
    reg192=reg22*reg173; reg75=reg22*reg75; reg200=reg22*reg200; reg110=reg22*reg110; reg124=reg22*reg124;
    reg123=reg22*reg123; reg136=reg22*reg136; reg117=reg22*reg117; reg52=reg22*reg52; reg143=ponderation*reg143;
    reg194=ponderation*reg114; reg205=ponderation*reg179; sollicitation[indices[2]+1]+=reg205; T tmp_7_0=ponderation*reg70; T tmp_3_2=-reg194;
    T tmp_7_2=ponderation*reg124; T tmp_2_7=ponderation*reg59; T tmp_6_5=ponderation*reg99; T tmp_2_5=ponderation*reg177; T tmp_3_4=ponderation*reg128;
    sollicitation[indices[3]+0]+=-reg201; T tmp_7_1=ponderation*reg75; T tmp_6_7=-reg122; T tmp_3_0=ponderation*reg121; T tmp_3_1=ponderation*reg123;
    reg59=ponderation*reg176; sollicitation[indices[3]+1]+=reg59; T tmp_2_6=-reg156; T tmp_3_3=ponderation*reg125; T tmp_6_6=ponderation*reg53;
    T tmp_6_0=ponderation*reg102; T tmp_6_1=ponderation*reg150; T tmp_6_2=-reg129; T tmp_4_5=ponderation*reg117; reg53=ponderation*reg187;
    sollicitation[indices[1]+0]+=reg53; reg70=ponderation*reg185; sollicitation[indices[0]+1]+=reg70; T tmp_4_6=ponderation*reg118; reg75=ponderation*reg192;
    sollicitation[indices[0]+0]+=reg75; T tmp_2_0=ponderation*reg200; T tmp_1_2=ponderation*reg64; T tmp_7_7=ponderation*reg136; T tmp_4_7=ponderation*reg110;
    T tmp_1_3=ponderation*reg52; T tmp_7_6=-reg135; T tmp_5_0=ponderation*reg155; T tmp_7_5=ponderation*reg174; T tmp_1_7=ponderation*reg67;
    T tmp_1_4=ponderation*reg180; T tmp_7_4=ponderation*reg172; T tmp_5_1=ponderation*reg157; T tmp_7_3=-reg116; T tmp_5_7=ponderation*reg169;
    T tmp_5_6=ponderation*reg168; T tmp_1_5=ponderation*reg60; T tmp_5_5=ponderation*reg166; T tmp_5_2=ponderation*reg158; T tmp_5_4=ponderation*reg164;
    T tmp_5_3=ponderation*reg161; T tmp_1_6=ponderation*reg203; T tmp_3_5=ponderation*reg130; T tmp_2_4=ponderation*reg71; T tmp_0_0=ponderation*reg151;
    T tmp_3_6=ponderation*reg131; T tmp_0_1=ponderation*reg186; T tmp_0_2=ponderation*reg184; T tmp_3_7=-reg127; T tmp_2_3=-reg120;
    T tmp_0_3=ponderation*reg181; T tmp_4_0=ponderation*reg134; T tmp_0_4=ponderation*reg178; T tmp_4_1=ponderation*reg107; T tmp_0_5=ponderation*reg175;
    T tmp_2_2=ponderation*reg74; T tmp_0_6=ponderation*reg198; T tmp_4_2=ponderation*reg109; T tmp_0_7=ponderation*reg197; T tmp_4_3=ponderation*reg113;
    T tmp_6_4=ponderation*reg149; T tmp_2_1=ponderation*reg188; T tmp_1_0=ponderation*reg193; T tmp_6_3=ponderation*reg146; T tmp_4_4=ponderation*reg115;
    reg52=ponderation*reg105; sollicitation[indices[2]+0]+=reg52; T tmp_1_1=ponderation*reg190; sollicitation[indices[1]+1]+=-reg143;
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
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg4*reg1; reg1=reg3*reg1; T reg6=reg4*reg2; reg2=reg3*reg2;
    T reg7=reg3*reg1; T reg8=reg3*reg5; T reg9=reg4*reg6; T reg10=reg3*reg2; reg6=reg3*reg6;
    reg5=reg4*reg5; reg9=reg9-reg10; reg1=reg4*reg1; reg6=reg10+reg6; reg2=reg4*reg2;
    reg5=reg5-reg7; reg8=reg7+reg8; reg5=reg4*reg5; T reg11=reg4*reg9; T reg12=reg7+reg1;
    T reg13=reg3*reg6; reg2=reg10+reg2; reg8=reg3*reg8; reg12=reg3*reg12; reg8=reg5-reg8;
    reg5=reg3*reg2; reg13=reg11-reg13; reg10=1-var_inter[1]; reg5=reg13-reg5; reg12=reg8-reg12;
    reg8=1-var_inter[0]; reg11=reg10*elem.pos(1)[0]; reg12=reg12/reg5; reg13=reg10*elem.pos(0)[0]; T reg14=reg10*elem.pos(0)[1];
    T reg15=reg10*elem.pos(1)[1]; T reg16=elem.pos(1)[1]*var_inter[0]; T reg17=reg8*elem.pos(0)[1]; reg9=reg9/reg5; reg6=reg6/reg5;
    T reg18=reg8*elem.pos(0)[0]; T reg19=elem.pos(1)[0]*var_inter[0]; T reg20=elem.pos(2)[1]*var_inter[0]; T reg21=reg17+reg16; T reg22=reg6*reg12;
    T reg23=reg9*reg12; reg15=reg15-reg14; T reg24=var_inter[1]*elem.pos(2)[1]; T reg25=elem.pos(2)[0]*var_inter[0]; T reg26=reg19+reg18;
    reg11=reg11-reg13; T reg27=var_inter[1]*elem.pos(2)[0]; T reg28=reg9*reg23; reg24=reg15+reg24; reg15=var_inter[1]*elem.pos(3)[1];
    T reg29=reg6*reg22; reg25=reg25-reg26; T reg30=elem.pos(3)[0]*reg8; T reg31=var_inter[1]*elem.pos(3)[0]; reg11=reg27+reg11;
    reg5=reg2/reg5; reg2=reg6*(*f.m).alpha; reg27=elem.pos(3)[1]*reg8; T reg32=reg9*(*f.m).alpha; reg20=reg20-reg21;
    reg24=reg24-reg15; reg29=reg28-reg29; reg5=(*f.m).alpha*reg5; reg27=reg20+reg27; reg2=reg32+reg2;
    reg30=reg25+reg30; reg11=reg11-reg31; reg5=reg2+reg5; reg2=reg27*reg11; reg23=reg23/reg29;
    reg20=reg3*reg0; reg25=reg24*reg30; reg29=reg22/reg29; reg22=reg4*reg0; reg23=reg23*reg5;
    reg5=reg29*reg5; reg28=pow(reg3,2); reg25=reg2-reg25; reg2=reg4*reg22; reg29=pow(reg4,2);
    reg32=reg3*reg20; reg11=reg11/reg25; reg30=reg30/reg25; reg27=reg27/reg25; reg24=reg24/reg25;
    T reg33=1-(*f.m).resolution; reg5=reg23-reg5; reg32=reg2-reg32; reg28=reg29-reg28; reg5=reg33*reg5;
    reg20=reg20/reg32; reg2=(*f.m).alpha*(*f.m).resolution; reg23=reg10*reg27; reg29=reg8*reg24; T reg34=var_inter[0]*reg24;
    reg28=reg28/reg32; T reg35=reg8*reg11; reg32=reg22/reg32; reg22=var_inter[1]*reg30; T reg36=var_inter[1]*reg27;
    reg9=reg9*reg33; T reg37=reg29+reg36; T reg38=reg35+reg22; T reg39=(*f.m).resolution*reg32; reg5=reg2+reg5;
    reg12=reg12*reg33; reg33=reg6*reg33; reg2=(*f.m).resolution*reg20; reg6=reg10*reg30; T reg40=(*f.m).resolution*reg28;
    T reg41=reg23+reg34; T reg42=var_inter[0]*reg11; reg39=reg9+reg39; reg9=reg6-reg35; T reg43=reg42-reg22;
    T reg44=0.5*reg37; T reg45=reg29-reg23; reg5=(*f.m).deltaT*reg5; T reg46=reg6+reg42; reg12=reg40+reg12;
    reg40=reg36-reg34; reg33=reg2+reg33; reg2=0.5*reg41; T reg47=0.5*reg38; T reg48=0.5*reg46;
    T reg49=reg12*reg2; T reg50=reg12*reg47; T reg51=0.5*reg40; T reg52=0.5*reg43; T reg53=0.5*reg45;
    T reg54=reg39*reg5; T reg55=0.5*reg9; T reg56=reg12*reg44; T reg57=reg33*reg5; T reg58=reg12*reg53;
    T reg59=reg10*var_inter[0]; T reg60=var_inter[1]*reg8; T reg61=reg39*reg38; T reg62=reg12*reg55; T reg63=reg54+reg57;
    T reg64=2*reg56; T reg65=reg12*reg48; T reg66=reg33*reg38; T reg67=reg33*reg46; reg50=2*reg50;
    T reg68=reg39*reg37; reg49=2*reg49; T reg69=reg12*reg51; T reg70=reg12*reg52; T reg71=reg63*reg37;
    T reg72=reg60*(*f.m).f_vol[0]; T reg73=reg33*reg37; T reg74=reg39*reg45; T reg75=var_inter[1]*var_inter[0]; T reg76=reg10*reg8;
    T reg77=reg63*reg46; T reg78=reg67*reg41; T reg79=reg49*reg48; T reg80=reg68*reg41; T reg81=reg50*reg48;
    T reg82=reg61*reg46; T reg83=reg64*reg2; reg58=2*reg58; T reg84=reg39*reg9; T reg85=reg39*reg41;
    T reg86=reg33*reg9; T reg87=2*reg65; reg70=2*reg70; T reg88=reg33*reg43; reg69=2*reg69;
    T reg89=reg33*reg41; reg62=2*reg62; T reg90=reg39*reg40; T reg91=reg39*reg46; T reg92=reg66*reg37;
    T reg93=reg64*reg47; T reg94=reg59*(*f.m).f_vol[1]; T reg95=reg33*reg40; T reg96=reg39*reg43; T reg97=reg66*reg40;
    T reg98=reg50*reg52; T reg99=reg64*reg52; T reg100=reg68*reg40; T reg101=reg96*reg43; T reg102=reg69*reg52;
    T reg103=reg88*reg40; T reg104=reg70*reg52; T reg105=reg90*reg40; reg82=reg83+reg82; T reg106=reg73*reg46;
    T reg107=reg50*reg2; T reg108=reg96*reg46; T reg109=reg69*reg2; T reg110=reg95*reg46; T reg111=reg75*(*f.m).f_vol[1];
    T reg112=reg75*(*f.m).f_vol[0]; T reg113=reg60*(*f.m).f_vol[1]; T reg114=reg59*(*f.m).f_vol[0]; T reg115=reg76*(*f.m).f_vol[1]; T reg116=reg45*reg74;
    T reg117=reg63*reg43; T reg118=reg63*reg40; T reg119=reg77-reg94; T reg120=reg63*reg41; T reg121=reg63*reg9;
    T reg122=reg63*reg45; T reg123=reg76*(*f.m).f_vol[0]; T reg124=reg64*reg44; T reg125=reg61*reg38; T reg126=reg64*reg51;
    T reg127=reg61*reg43; T reg128=reg50*reg51; T reg129=reg73*reg43; T reg130=reg69*reg51; T reg131=reg69*reg48;
    T reg132=reg88*reg45; reg81=reg80+reg81; T reg133=reg69*reg55; T reg134=reg68*reg45; T reg135=reg50*reg55;
    T reg136=reg66*reg45; T reg137=reg64*reg55; T reg138=reg84*reg9; T reg139=reg88*reg41; T reg140=reg58*reg53;
    T reg141=reg89*reg9; T reg142=reg87*reg53; T reg143=reg91*reg9; T reg144=reg49*reg53; T reg145=reg95*reg9;
    T reg146=reg70*reg48; T reg147=reg70*reg53; T reg148=reg96*reg9; T reg149=reg90*reg41; T reg150=reg69*reg53;
    T reg151=reg73*reg9; T reg152=reg50*reg53; reg61=reg61*reg9; reg79=reg78+reg79; T reg153=reg64*reg53;
    T reg154=reg85*reg41; T reg155=reg87*reg48; T reg156=reg50*reg47; T reg157=reg68*reg37; T reg158=reg91*reg46;
    T reg159=reg49*reg2; reg92=reg93+reg92; T reg160=reg63*reg38; T reg161=reg62*reg55; T reg162=reg71-reg72;
    T reg163=reg86*reg45; T reg164=reg58*reg55; T reg165=reg85*reg45; T reg166=reg70*reg55; T reg167=reg90*reg45;
    T reg168=reg49*reg55; T reg169=reg67*reg45; T reg170=reg70*reg2; reg66=reg66*reg41; T reg171=reg64*reg48;
    T reg172=reg87*reg55; T reg173=reg117+reg111; reg162=reg25*reg162; reg119=reg25*reg119; reg152=reg152-reg151;
    reg61=reg61-reg153; T reg174=reg120+reg114; T reg175=reg118+reg112; T reg176=reg160+reg113; reg150=reg148+reg150;
    reg154=reg154+reg155; reg156=reg156+reg157; reg147=reg145+reg147; reg144=reg144-reg143; reg145=reg25*reg92;
    reg141=reg141-reg142; reg116=reg161+reg116; reg140=reg138+reg140; reg164=reg163+reg164; reg136=reg136-reg137;
    reg135=reg135-reg134; reg165=reg165-reg172; reg133=reg132+reg133; reg168=reg168-reg169; reg166=reg167+reg166;
    reg107=reg107+reg106; reg128=reg128-reg129; reg130=reg101+reg130; reg131=reg139-reg131; reg146=reg149-reg146;
    reg127=reg127-reg126; reg108=reg109-reg108; reg97=reg97-reg99; reg125=reg125+reg124; reg159=reg159+reg158;
    reg98=reg98-reg100; reg101=reg25*reg81; reg109=reg25*reg79; reg132=reg122+reg123; reg104=reg105+reg104;
    reg110=reg170-reg110; reg105=reg121+reg115; reg138=reg25*reg82; reg66=reg66+reg171; reg102=reg103+reg102;
    reg107=reg25*reg107; reg140=reg25*reg140; reg165=reg25*reg165; reg130=reg25*reg130; reg103=ponderation*reg138;
    reg66=reg25*reg66; reg139=ponderation*reg145; reg166=reg25*reg166; reg104=reg25*reg104; reg136=reg25*reg136;
    reg97=reg25*reg97; reg164=reg25*reg164; reg168=reg25*reg168; reg135=reg25*reg135; reg102=reg25*reg102;
    reg131=reg25*reg131; reg148=ponderation*reg101; reg116=reg25*reg116; reg133=reg25*reg133; reg98=reg25*reg98;
    reg149=reg25*reg132; reg161=reg25*reg175; reg110=reg25*reg110; reg61=reg25*reg61; reg162=ponderation*reg162;
    reg152=reg25*reg152; reg163=ponderation*reg109; reg159=reg25*reg159; reg119=ponderation*reg119; reg125=reg25*reg125;
    reg167=reg25*reg173; reg150=reg25*reg150; reg170=reg25*reg176; reg154=reg25*reg154; reg147=reg25*reg147;
    T reg177=reg25*reg174; reg127=reg25*reg127; reg108=reg25*reg108; reg144=reg25*reg144; T reg178=reg25*reg105;
    reg128=reg25*reg128; reg156=reg25*reg156; reg141=reg25*reg141; reg146=reg25*reg146; T tmp_0_1=ponderation*reg164;
    reg164=ponderation*reg167; sollicitation[indices[2]+1]+=reg164; T tmp_3_3=ponderation*reg159; T tmp_3_4=ponderation*reg110; sollicitation[indices[3]+0]+=-reg162;
    T tmp_0_0=ponderation*reg116; T tmp_2_7=ponderation*reg66; T tmp_3_6=ponderation*reg107; reg66=ponderation*reg170; sollicitation[indices[3]+1]+=reg66;
    T tmp_6_7=-reg139; T tmp_3_5=ponderation*reg108; T tmp_6_6=ponderation*reg156; reg107=ponderation*reg178; sollicitation[indices[0]+1]+=reg107;
    reg108=ponderation*reg177; sollicitation[indices[1]+0]+=reg108; reg110=ponderation*reg149; sollicitation[indices[0]+0]+=reg110; sollicitation[indices[1]+1]+=-reg119;
    T tmp_1_7=ponderation*reg61; reg61=ponderation*reg161; sollicitation[indices[2]+0]+=reg61; T tmp_2_2=ponderation*reg154; T tmp_1_6=ponderation*reg152;
    T tmp_7_7=ponderation*reg125; T tmp_1_5=ponderation*reg150; T tmp_2_3=-reg163; T tmp_5_7=ponderation*reg127; T tmp_1_4=ponderation*reg147;
    T tmp_1_3=ponderation*reg144; T tmp_5_6=ponderation*reg128; T tmp_1_2=ponderation*reg141; T tmp_2_4=ponderation*reg146; T tmp_5_5=ponderation*reg130;
    T tmp_1_1=ponderation*reg140; T tmp_0_7=ponderation*reg136; T tmp_4_7=ponderation*reg97; T tmp_0_6=ponderation*reg135; T tmp_2_5=ponderation*reg131;
    T tmp_4_6=ponderation*reg98; T tmp_0_5=ponderation*reg133; T tmp_4_5=ponderation*reg102; T tmp_0_4=ponderation*reg166; T tmp_0_3=ponderation*reg168;
    T tmp_2_6=-reg148; T tmp_4_4=ponderation*reg104; T tmp_0_2=ponderation*reg165; T tmp_3_7=-reg103;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=1-var_inter[0]; T reg2=1-var_inter[1]; T reg3=pow(reg0,2);
    T reg4=reg1*elem.pos(0)[1]; T reg5=reg2*elem.pos(0)[0]; T reg6=reg2*elem.pos(1)[0]; T reg7=elem.pos(1)[1]*var_inter[0]; T reg8=1.0/(*f.m).elastic_modulus;
    T reg9=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg10=reg0*reg3; T reg11=reg1*elem.pos(0)[0]; T reg12=elem.pos(1)[0]*var_inter[0]; T reg13=reg2*elem.pos(0)[1];
    T reg14=reg2*elem.pos(1)[1]; T reg15=reg12+reg11; T reg16=elem.pos(2)[0]*var_inter[0]; T reg17=reg8*reg10; reg10=reg9*reg10;
    T reg18=reg9*reg3; T reg19=elem.pos(2)[1]*var_inter[0]; reg3=reg8*reg3; T reg20=reg4+reg7; T reg21=var_inter[1]*elem.pos(2)[1];
    reg14=reg14-reg13; reg6=reg6-reg5; T reg22=var_inter[1]*elem.pos(2)[0]; T reg23=reg9*reg3; T reg24=reg9*reg18;
    reg3=reg8*reg3; reg19=reg19-reg20; T reg25=elem.pos(3)[1]*reg1; T reg26=reg9*reg17; T reg27=reg9*reg10;
    reg17=reg8*reg17; T reg28=var_inter[1]*elem.pos(3)[1]; reg21=reg14+reg21; reg16=reg16-reg15; reg6=reg22+reg6;
    reg14=elem.pos(3)[0]*reg1; reg22=var_inter[1]*elem.pos(3)[0]; reg26=reg27+reg26; reg17=reg17-reg27; reg14=reg16+reg14;
    reg10=reg8*reg10; reg25=reg19+reg25; reg6=reg6-reg22; reg21=reg21-reg28; reg18=reg8*reg18;
    reg23=reg24+reg23; reg3=reg3-reg24; reg3=reg8*reg3; reg10=reg27+reg10; reg23=reg9*reg23;
    reg16=reg24+reg18; reg19=reg8*reg17; reg27=reg25*reg6; T reg29=reg9*reg0; T reg30=reg21*reg14;
    T reg31=reg9*reg26; T reg32=reg8*reg0; T reg33=reg8*reg32; reg23=reg3-reg23; reg16=reg9*reg16;
    reg3=pow(reg9,2); T reg34=reg9*reg29; T reg35=pow(reg8,2); reg31=reg19-reg31; reg19=reg9*reg10;
    reg30=reg27-reg30; reg3=reg35-reg3; reg19=reg31-reg19; reg21=reg21/reg30; reg16=reg23-reg16;
    reg25=reg25/reg30; reg34=reg33-reg34; reg6=reg6/reg30; reg14=reg14/reg30; reg3=reg3/reg34;
    reg23=1-(*f.m).resolution; reg27=reg2*reg25; reg31=reg1*reg21; reg33=reg1*reg6; reg35=reg2*reg14;
    T reg36=var_inter[0]*reg21; T reg37=var_inter[0]*reg6; T reg38=var_inter[1]*reg25; reg16=reg16/reg19; T reg39=var_inter[1]*reg14;
    reg32=reg32/reg34; reg34=reg29/reg34; reg29=(*f.m).resolution*reg3; T reg40=reg16*reg23; reg17=reg17/reg19;
    reg26=reg26/reg19; T reg41=reg33+reg39; T reg42=reg31+reg38; T reg43=reg35+reg37; T reg44=reg27+reg36;
    T reg45=0.5*reg42; T reg46=0.5*reg41; T reg47=reg37-reg39; T reg48=reg38-reg36; T reg49=0.5*reg44;
    T reg50=0.5*reg43; T reg51=reg35-reg33; T reg52=reg26*reg23; T reg53=(*f.m).resolution*reg34; T reg54=reg17*reg23;
    T reg55=(*f.m).resolution*reg32; T reg56=reg31-reg27; reg40=reg29+reg40; reg29=reg40*reg46; T reg57=reg40*reg45;
    reg55=reg54+reg55; reg54=0.5*reg48; T reg58=reg40*reg49; T reg59=0.5*reg47; reg52=reg53+reg52;
    reg53=reg40*reg50; T reg60=0.5*reg51; T reg61=0.5*reg56; T reg62=reg55*reg43; reg29=2*reg29;
    T reg63=reg55*reg42; T reg64=reg52*reg41; T reg65=reg40*reg54; T reg66=2*reg57; T reg67=reg40*reg61;
    T reg68=reg40*reg59; T reg69=reg52*reg42; T reg70=reg40*reg60; T reg71=2*reg53; T reg72=reg52*reg43;
    T reg73=reg55*reg44; reg58=2*reg58; T reg74=reg52*reg44; T reg75=reg55*reg41; T reg76=reg52*reg48;
    T reg77=reg52*reg56; T reg78=reg66*reg49; T reg79=reg55*reg47; T reg80=reg74*reg43; T reg81=reg69*reg41;
    T reg82=reg55*reg51; T reg83=reg29*reg45; T reg84=reg29*reg50; T reg85=reg63*reg44; T reg86=reg75*reg43;
    T reg87=reg71*reg49; T reg88=reg64*reg42; T reg89=reg66*reg46; reg70=2*reg70; T reg90=reg73*reg42;
    T reg91=reg52*reg51; reg67=2*reg67; T reg92=reg58*reg45; T reg93=reg62*reg41; T reg94=reg55*reg48;
    reg68=2*reg68; T reg95=reg55*reg56; T reg96=reg52*reg47; T reg97=reg72*reg44; T reg98=reg58*reg50;
    T reg99=reg71*reg46; reg65=2*reg65; T reg100=reg72*reg48; T reg101=reg71*reg59; T reg102=reg73*reg48;
    T reg103=reg67*reg59; T reg104=reg91*reg48; T reg105=reg96*reg42; T reg106=reg70*reg59; T reg107=reg48*reg95;
    reg86=reg78+reg86; reg84=reg85+reg84; T reg108=reg96*reg44; T reg109=reg69*reg43; T reg110=reg77*reg43;
    T reg111=reg70*reg49; T reg112=reg67*reg49; T reg113=reg82*reg43; T reg114=reg74*reg41; T reg115=reg66*reg50;
    reg80=reg87+reg80; T reg116=reg64*reg44; T reg117=reg58*reg49; T reg118=reg62*reg43; T reg119=reg68*reg49;
    T reg120=reg76*reg43; T reg121=reg65*reg49; T reg122=reg94*reg44; T reg123=reg68*reg50; T reg124=reg79*reg43;
    T reg125=reg29*reg49; T reg126=reg79*reg47; T reg127=reg65*reg54; T reg128=reg69*reg47; T reg129=reg29*reg54;
    T reg130=reg75*reg47; T reg131=reg66*reg54; T reg132=reg76*reg41; T reg133=reg68*reg45; T reg134=reg79*reg41;
    T reg135=reg65*reg45; reg83=reg81+reg83; T reg136=reg75*reg41; T reg137=reg66*reg45; T reg138=reg58*reg46;
    T reg139=reg72*reg42; T reg140=reg68*reg46; T reg141=reg94*reg42; T reg142=reg91*reg42; T reg143=reg56*reg95;
    T reg144=reg58*reg59; T reg145=reg94*reg48; T reg146=reg68*reg59; T reg147=reg96*reg48; T reg148=reg65*reg59;
    T reg149=reg63*reg48; T reg150=reg29*reg59; T reg151=reg64*reg48; T reg152=reg66*reg59; T reg153=reg77*reg47;
    T reg154=reg70*reg54; T reg155=reg82*reg47; T reg156=reg67*reg54; T reg157=reg74*reg47; T reg158=reg71*reg54;
    T reg159=reg62*reg47; T reg160=reg58*reg54; T reg161=reg76*reg47; T reg162=reg68*reg54; reg79=reg79*reg51;
    T reg163=reg58*reg60; T reg164=reg72*reg56; T reg165=reg65*reg61; T reg166=reg71*reg60; T reg167=reg73*reg56;
    T reg168=reg69*reg51; T reg169=reg29*reg61; T reg170=reg67*reg60; T reg171=reg91*reg56; T reg172=reg70*reg60;
    reg75=reg75*reg51; T reg173=reg66*reg61; T reg174=reg42*reg95; T reg175=reg82*reg41; T reg176=reg70*reg45;
    T reg177=reg77*reg41; reg95=reg44*reg95; T reg178=reg70*reg50; reg77=reg77*reg51; T reg179=reg70*reg61;
    reg82=reg82*reg51; T reg180=reg66*reg60; T reg181=reg67*reg61; reg64=reg64*reg56; reg74=reg74*reg51;
    T reg182=reg71*reg61; T reg183=reg29*reg60; T reg184=reg63*reg56; T reg185=reg62*reg51; reg58=reg58*reg61;
    T reg186=reg65*reg60; reg96=reg96*reg56; T reg187=reg65*reg50; reg76=reg76*reg51; T reg188=reg68*reg61;
    reg68=reg68*reg60; reg94=reg94*reg56; reg29=reg29*reg46; T reg189=reg63*reg42; reg70=reg70*reg46;
    reg65=reg65*reg46; reg98=reg97+reg98; T reg190=reg67*reg50; T reg191=reg67*reg45; reg90=reg90+reg99;
    reg91=reg91*reg44; reg92=reg93+reg92; reg73=reg73*reg44; reg88=reg89+reg88; T reg192=reg71*reg45;
    reg67=reg67*reg46; T reg193=reg71*reg50; reg187=reg108-reg187; reg68=reg94+reg68; reg154=reg153+reg154;
    reg114=reg192+reg114; reg138=reg138+reg139; reg94=reg30*reg84; reg148=reg147+reg148; reg151=reg151-reg152;
    reg191=reg175-reg191; reg58=reg58-reg185; reg73=reg73+reg193; reg186=reg96+reg186; reg141=reg140-reg141;
    reg150=reg150-reg149; reg130=reg130-reg131; reg129=reg129-reg128; reg179=reg77+reg179; reg133=reg132-reg133;
    reg127=reg126+reg127; reg123=reg122-reg123; reg64=reg64-reg180; reg174=reg70-reg174; reg162=reg161+reg162;
    reg181=reg82+reg181; reg135=reg134-reg135; reg160=reg160-reg159; reg70=reg30*reg98; reg77=reg30*reg83;
    reg157=reg157-reg158; reg105=reg65-reg105; reg183=reg183-reg184; reg156=reg155+reg156; reg74=reg74-reg182;
    reg136=reg136+reg137; reg170=reg171+reg170; reg106=reg107+reg106; reg110=reg111-reg110; reg190=reg91-reg190;
    reg65=reg30*reg86; reg169=reg169-reg168; reg143=reg172+reg143; reg142=reg67-reg142; reg125=reg125+reg109;
    reg113=reg112-reg113; reg29=reg29+reg189; reg124=reg121-reg124; reg75=reg75-reg173; reg176=reg177-reg176;
    reg67=reg30*reg90; reg120=reg119-reg120; reg82=reg30*reg80; reg91=reg30*reg88; reg117=reg117+reg118;
    reg178=reg95-reg178; reg163=reg163-reg164; reg144=reg144-reg100; reg116=reg116+reg115; reg102=reg102-reg101;
    reg95=reg30*reg92; reg165=reg79+reg165; reg167=reg167-reg166; reg188=reg76+reg188; reg103=reg104+reg103;
    reg146=reg145+reg146; reg76=ponderation*reg91; reg174=reg30*reg174; reg64=reg30*reg64; reg130=reg30*reg130;
    reg114=reg30*reg114; reg68=reg30*reg68; reg79=ponderation*reg67; reg133=reg30*reg133; reg163=reg30*reg163;
    reg29=reg30*reg29; reg176=reg30*reg176; reg142=reg30*reg142; reg191=reg30*reg191; reg141=reg30*reg141;
    reg135=reg30*reg135; reg183=reg30*reg183; reg96=ponderation*reg77; reg167=reg30*reg167; reg143=reg30*reg143;
    reg186=reg30*reg186; reg138=reg30*reg138; reg104=ponderation*reg95; reg136=reg30*reg136; reg170=reg30*reg170;
    reg156=reg30*reg156; reg74=reg30*reg74; reg116=reg30*reg116; reg120=reg30*reg120; reg154=reg30*reg154;
    reg190=reg30*reg190; reg165=reg30*reg165; reg110=reg30*reg110; reg151=reg30*reg151; reg103=reg30*reg103;
    reg150=reg30*reg150; reg58=reg30*reg58; reg113=reg30*reg113; reg148=reg30*reg148; reg178=reg30*reg178;
    reg102=reg30*reg102; reg146=reg30*reg146; reg188=reg30*reg188; reg107=ponderation*reg82; reg117=reg30*reg117;
    reg144=reg30*reg144; reg125=reg30*reg125; reg129=reg30*reg129; reg108=ponderation*reg70; reg179=reg30*reg179;
    reg123=reg30*reg123; reg124=reg30*reg124; reg127=reg30*reg127; reg105=reg30*reg105; reg169=reg30*reg169;
    reg187=reg30*reg187; reg162=reg30*reg162; reg160=reg30*reg160; reg111=ponderation*reg65; reg181=reg30*reg181;
    reg75=reg30*reg75; reg73=reg30*reg73; reg112=ponderation*reg94; reg157=reg30*reg157; reg106=reg30*reg106;
    T tmp_7_0=ponderation*reg176; T tmp_6_7=-reg76; T tmp_3_3=ponderation*reg117; T tmp_1_7=ponderation*reg75; T tmp_3_5=ponderation*reg124;
    T tmp_3_4=ponderation*reg120; T tmp_2_3=-reg108; T tmp_2_4=ponderation*reg123; T tmp_6_1=ponderation*reg142; T tmp_6_0=ponderation*reg174;
    T tmp_2_5=ponderation*reg187; T tmp_2_2=ponderation*reg73; T tmp_7_1=ponderation*reg191; T tmp_2_6=-reg112; T tmp_7_2=ponderation*reg114;
    T tmp_2_7=ponderation*reg116; T tmp_2_1=ponderation*reg190; T tmp_6_5=ponderation*reg105; T tmp_3_0=ponderation*reg110; T tmp_3_1=ponderation*reg113;
    T tmp_2_0=ponderation*reg178; T tmp_6_2=-reg79; T tmp_6_6=ponderation*reg29; T tmp_3_2=-reg107; T tmp_4_7=ponderation*reg151;
    T tmp_1_2=ponderation*reg74; T tmp_0_5=ponderation*reg186; T tmp_7_7=ponderation*reg136; T tmp_5_0=ponderation*reg154; T tmp_5_1=ponderation*reg156;
    T tmp_7_6=-reg96; T tmp_5_2=ponderation*reg157; T tmp_1_1=ponderation*reg181; T tmp_0_6=ponderation*reg183; T tmp_7_5=ponderation*reg135;
    T tmp_5_3=ponderation*reg160; T tmp_5_4=ponderation*reg162; T tmp_7_4=ponderation*reg133; T tmp_1_0=ponderation*reg179; T tmp_5_5=ponderation*reg127;
    T tmp_0_7=ponderation*reg64; T tmp_7_3=-reg104; T tmp_5_7=ponderation*reg130; T tmp_5_6=ponderation*reg129; T tmp_3_6=ponderation*reg125;
    T tmp_1_6=ponderation*reg169; T tmp_0_0=ponderation*reg143; T tmp_3_7=-reg111; T tmp_4_0=ponderation*reg106; T tmp_0_1=ponderation*reg170;
    T tmp_1_5=ponderation*reg165; T tmp_4_1=ponderation*reg103; T tmp_0_2=ponderation*reg167; T tmp_4_2=ponderation*reg102; T tmp_1_4=ponderation*reg188;
    T tmp_0_3=ponderation*reg163; T tmp_4_3=ponderation*reg144; T tmp_4_4=ponderation*reg146; T tmp_0_4=ponderation*reg68; T tmp_4_5=ponderation*reg148;
    T tmp_1_3=ponderation*reg58; T tmp_6_4=ponderation*reg141; T tmp_4_6=ponderation*reg150; T tmp_6_3=ponderation*reg138;
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
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg2*elem.pos(0)[0]; T reg6=elem.pos(1)[0]*var_inter[0]; T reg7=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg8=elem.pos(1)[1]*var_inter[0];
    T reg9=reg3*elem.pos(1)[1]; T reg10=reg3*elem.pos(0)[1]; T reg11=reg2*elem.pos(0)[1]; T reg12=reg0*reg1; T reg13=reg3*elem.pos(0)[0];
    T reg14=reg3*elem.pos(1)[0]; T reg15=reg11+reg8; T reg16=reg6+reg5; T reg17=elem.pos(2)[0]*var_inter[0]; T reg18=reg4*reg1;
    T reg19=var_inter[1]*elem.pos(2)[1]; reg9=reg9-reg10; reg1=reg7*reg1; T reg20=elem.pos(2)[1]*var_inter[0]; T reg21=reg4*reg12;
    reg14=reg14-reg13; reg12=reg7*reg12; T reg22=var_inter[1]*elem.pos(2)[0]; T reg23=reg4*reg18; T reg24=reg7*reg1;
    T reg25=reg7*reg21; reg18=reg7*reg18; T reg26=reg7*reg12; reg21=reg4*reg21; T reg27=elem.pos(3)[0]*reg2;
    reg17=reg17-reg16; T reg28=var_inter[1]*elem.pos(3)[1]; reg19=reg9+reg19; reg9=elem.pos(3)[1]*reg2; T reg29=var_inter[1]*elem.pos(3)[0];
    reg14=reg22+reg14; reg20=reg20-reg15; reg12=reg4*reg12; reg9=reg20+reg9; reg25=reg26+reg25;
    reg21=reg21-reg26; reg1=reg4*reg1; reg23=reg23-reg24; reg14=reg14-reg29; reg18=reg24+reg18;
    reg19=reg19-reg28; reg27=reg17+reg27; reg17=reg24+reg1; reg18=reg7*reg18; reg23=reg4*reg23;
    reg20=reg7*reg0; reg22=reg9*reg14; reg12=reg26+reg12; reg26=reg19*reg27; T reg30=reg4*reg0;
    T reg31=reg4*reg21; T reg32=reg7*reg25; T reg33=reg7*reg20; T reg34=pow(reg4,2); reg32=reg31-reg32;
    reg31=pow(reg7,2); reg26=reg22-reg26; reg17=reg7*reg17; reg22=reg4*reg30; T reg35=reg7*reg12;
    reg18=reg23-reg18; reg17=reg18-reg17; reg33=reg22-reg33; reg31=reg34-reg31; reg19=reg19/reg26;
    reg27=reg27/reg26; reg14=reg14/reg26; reg35=reg32-reg35; reg9=reg9/reg26; reg18=1-(*f.m).resolution;
    reg22=reg3*reg9; reg31=reg31/reg33; reg23=reg2*reg19; reg32=reg2*reg14; reg17=reg17/reg35;
    reg34=var_inter[1]*reg9; T reg36=var_inter[1]*reg27; T reg37=var_inter[0]*reg19; T reg38=reg17*reg18; T reg39=(*f.m).resolution*reg31;
    T reg40=var_inter[0]*reg14; T reg41=reg22+reg37; reg20=reg20/reg33; T reg42=reg23+reg34; T reg43=reg32+reg36;
    T reg44=reg3*reg27; reg21=reg21/reg35; reg25=reg25/reg35; reg33=reg30/reg33; reg30=(*f.m).resolution*reg20;
    T reg45=reg25*reg18; T reg46=0.5*reg43; T reg47=0.5*reg42; T reg48=reg40-reg36; T reg49=reg34-reg37;
    T reg50=0.5*reg41; T reg51=reg44+reg40; T reg52=reg44-reg32; T reg53=reg23-reg22; reg38=reg39+reg38;
    reg39=reg21*reg18; T reg54=(*f.m).resolution*reg33; T reg55=0.5*reg49; T reg56=0.5*reg48; reg54=reg39+reg54;
    reg39=0.5*reg51; T reg57=0.5*reg53; T reg58=reg38*reg47; T reg59=0.5*reg52; reg45=reg30+reg45;
    reg30=reg38*reg50; T reg60=reg38*reg46; reg30=2*reg30; T reg61=reg38*reg59; T reg62=reg38*reg56;
    T reg63=reg54*reg42; T reg64=reg45*reg51; T reg65=reg38*reg57; T reg66=reg54*reg43; T reg67=2*reg58;
    reg60=2*reg60; T reg68=reg38*reg55; T reg69=reg45*reg43; T reg70=reg38*reg39; T reg71=reg66*reg51;
    T reg72=reg45*reg42; T reg73=reg54*reg48; T reg74=reg60*reg39; T reg75=reg45*reg49; T reg76=reg54*reg53;
    T reg77=reg67*reg50; T reg78=reg64*reg41; T reg79=reg63*reg41; T reg80=reg30*reg39; reg62=2*reg62;
    T reg81=reg45*reg48; T reg82=reg54*reg49; reg68=2*reg68; T reg83=2*reg70; T reg84=reg54*reg41;
    reg65=2*reg65; T reg85=reg54*reg52; T reg86=reg45*reg52; T reg87=reg45*reg41; T reg88=reg67*reg46;
    T reg89=reg69*reg42; reg61=2*reg61; T reg90=reg54*reg51; T reg91=reg72*reg51; T reg92=reg60*reg50;
    T reg93=reg73*reg51; T reg94=reg68*reg50; T reg95=reg75*reg51; T reg96=reg62*reg50; T reg97=reg90*reg51;
    T reg98=reg30*reg50; T reg99=reg67*reg39; T reg100=reg69*reg41; reg74=reg79+reg74; T reg101=reg81*reg41;
    T reg102=reg66*reg48; T reg103=reg67*reg55; T reg104=reg60*reg55; T reg105=reg72*reg48; T reg106=reg68*reg55;
    T reg107=reg73*reg48; T reg108=reg67*reg56; T reg109=reg66*reg43; T reg110=reg69*reg49; T reg111=reg60*reg56;
    T reg112=reg63*reg49; T reg113=reg68*reg56; T reg114=reg81*reg49; T reg115=reg67*reg47; T reg116=reg62*reg56;
    T reg117=reg82*reg49; reg71=reg77+reg71; T reg118=reg53*reg76; T reg119=reg82*reg53; reg80=reg78+reg80;
    T reg120=reg62*reg59; T reg121=reg67*reg57; T reg122=reg68*reg39; T reg123=reg81*reg53; T reg124=reg68*reg59;
    T reg125=reg63*reg53; T reg126=reg60*reg59; T reg127=reg62*reg57; T reg128=reg83*reg39; T reg129=reg84*reg41;
    reg69=reg69*reg53; T reg130=reg67*reg59; T reg131=reg83*reg57; T reg132=reg85*reg52; T reg133=reg65*reg57;
    T reg134=reg75*reg52; T reg135=reg87*reg52; T reg136=reg30*reg57; T reg137=reg60*reg46; T reg138=reg63*reg42;
    T reg139=reg72*reg52; reg89=reg88+reg89; T reg140=reg90*reg52; T reg141=reg62*reg39; T reg142=reg61*reg59;
    T reg143=reg68*reg57; T reg144=reg86*reg53; T reg145=reg82*reg41; T reg146=reg65*reg59; T reg147=reg73*reg52;
    T reg148=reg84*reg53; T reg149=reg83*reg59; T reg150=reg60*reg57; T reg151=reg64*reg53; T reg152=reg30*reg59;
    reg66=reg66*reg52; reg136=reg136-reg140; reg109=reg109+reg115; reg127=reg134+reg127; reg102=reg102-reg103;
    reg104=reg104-reg105; reg106=reg107+reg106; reg143=reg147+reg143; reg135=reg135-reg131; reg133=reg132+reg133;
    reg69=reg69-reg130; reg126=reg126-reg125; reg124=reg123+reg124; reg120=reg119+reg120; reg152=reg152-reg151;
    reg148=reg148-reg149; reg146=reg144+reg146; reg118=reg142+reg118; reg107=reg26*reg89; reg137=reg137+reg138;
    reg119=reg26*reg71; reg66=reg66-reg121; reg116=reg117+reg116; reg100=reg100+reg99; reg117=reg26*reg74;
    reg92=reg92+reg91; reg113=reg114+reg113; reg141=reg145-reg141; reg93=reg94-reg93; reg129=reg129+reg128;
    reg110=reg110-reg108; reg98=reg98+reg97; reg122=reg101-reg122; reg94=reg26*reg80; reg95=reg96-reg95;
    reg150=reg150-reg139; reg111=reg111-reg112; reg133=reg26*reg133; reg98=reg26*reg98; reg66=reg26*reg66;
    reg152=reg26*reg152; reg95=reg26*reg95; reg69=reg26*reg69; reg146=reg26*reg146; reg96=ponderation*reg94;
    reg92=reg26*reg92; reg129=reg26*reg129; reg126=reg26*reg126; reg100=reg26*reg100; reg93=reg26*reg93;
    reg120=reg26*reg120; reg148=reg26*reg148; reg124=reg26*reg124; reg110=reg26*reg110; reg143=reg26*reg143;
    reg122=reg26*reg122; reg106=reg26*reg106; reg111=reg26*reg111; reg104=reg26*reg104; reg137=reg26*reg137;
    reg127=reg26*reg127; reg141=reg26*reg141; reg102=reg26*reg102; reg113=reg26*reg113; reg150=reg26*reg150;
    reg101=ponderation*reg119; reg135=reg26*reg135; reg118=reg26*reg118; reg116=reg26*reg116; reg114=ponderation*reg117;
    reg136=reg26*reg136; reg123=ponderation*reg107; reg109=reg26*reg109; T tmp_2_7=ponderation*reg100; T tmp_6_6=ponderation*reg137;
    T tmp_6_7=-reg123; T tmp_2_3=-reg96; T tmp_2_5=ponderation*reg122; T tmp_0_1=ponderation*reg146; T tmp_0_2=ponderation*reg148;
    T tmp_0_0=ponderation*reg118; T tmp_2_6=-reg114; T tmp_2_4=ponderation*reg141; T tmp_4_7=ponderation*reg110; T tmp_1_5=ponderation*reg143;
    T tmp_4_6=ponderation*reg111; T tmp_5_5=ponderation*reg106; T tmp_5_6=ponderation*reg104; T tmp_1_4=ponderation*reg127; T tmp_4_5=ponderation*reg113;
    T tmp_5_7=ponderation*reg102; T tmp_7_7=ponderation*reg109; T tmp_1_6=ponderation*reg150; T tmp_1_3=ponderation*reg136; T tmp_4_4=ponderation*reg116;
    T tmp_1_2=ponderation*reg135; T tmp_3_7=-reg101; T tmp_1_1=ponderation*reg133; T tmp_1_7=ponderation*reg66; T tmp_0_7=ponderation*reg69;
    T tmp_3_6=ponderation*reg92; T tmp_0_6=ponderation*reg126; T tmp_3_5=ponderation*reg93; T tmp_0_5=ponderation*reg124; T tmp_3_4=ponderation*reg95;
    T tmp_0_4=ponderation*reg120; T tmp_2_2=ponderation*reg129; T tmp_3_3=ponderation*reg98; T tmp_0_3=ponderation*reg152;
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
    T reg4=1.0/(*f.m).elastic_modulus; T reg5=reg3*reg2; reg2=reg4*reg2; T reg6=reg3*reg1; reg1=reg4*reg1;
    T reg7=reg3*reg2; T reg8=reg3*reg5; reg2=reg4*reg2; T reg9=reg4*reg1; T reg10=reg3*reg6;
    reg1=reg3*reg1; reg9=reg9-reg10; reg6=reg4*reg6; reg7=reg8+reg7; reg1=reg10+reg1;
    reg2=reg2-reg8; reg5=reg4*reg5; reg1=reg3*reg1; reg9=reg4*reg9; reg5=reg8+reg5;
    reg8=reg3*reg7; T reg11=reg4*reg2; T reg12=reg10+reg6; T reg13=reg3*reg5; reg8=reg11-reg8;
    reg1=reg9-reg1; reg12=reg3*reg12; reg12=reg1-reg12; reg13=reg8-reg13; reg12=reg12/reg13;
    reg2=reg2/reg13; reg7=reg7/reg13; reg1=reg7*reg12; reg8=reg2*reg12; reg9=reg7*(*f.m).alpha;
    reg11=reg2*(*f.m).alpha; T reg14=1-var_inter[1]; T reg15=1-var_inter[0]; T reg16=reg7*reg1; reg13=reg5/reg13;
    reg5=reg2*reg8; T reg17=reg15*elem.pos(0)[0]; T reg18=elem.pos(1)[0]*var_inter[0]; reg16=reg5-reg16; reg9=reg11+reg9;
    reg5=reg15*elem.pos(0)[1]; reg13=(*f.m).alpha*reg13; reg11=elem.pos(1)[1]*var_inter[0]; T reg19=reg14*elem.pos(1)[0]; T reg20=reg14*elem.pos(0)[0];
    T reg21=reg14*elem.pos(1)[1]; T reg22=reg14*elem.pos(0)[1]; reg13=reg9+reg13; reg9=reg3*reg0; T reg23=elem.pos(2)[1]*var_inter[0];
    T reg24=reg5+reg11; T reg25=elem.pos(2)[0]*var_inter[0]; reg1=reg1/reg16; T reg26=reg18+reg17; T reg27=reg4*reg0;
    T reg28=var_inter[1]*elem.pos(2)[1]; reg21=reg21-reg22; reg16=reg8/reg16; reg8=var_inter[1]*elem.pos(2)[0]; reg19=reg19-reg20;
    reg1=reg1*reg13; reg13=reg16*reg13; reg16=reg4*reg27; T reg29=reg3*reg9; reg23=reg23-reg24;
    T reg30=elem.pos(3)[0]*reg15; reg25=reg25-reg26; T reg31=var_inter[1]*elem.pos(3)[1]; T reg32=elem.pos(3)[1]*reg15; reg28=reg21+reg28;
    reg21=var_inter[1]*elem.pos(3)[0]; reg19=reg8+reg19; reg30=reg25+reg30; reg28=reg28-reg31; reg29=reg16-reg29;
    reg19=reg19-reg21; reg1=reg13-reg1; reg32=reg23+reg32; reg8=1-(*f.m).resolution; reg9=reg9/reg29;
    reg13=reg32*reg19; reg16=(*f.m).alpha*(*f.m).resolution; reg1=reg8*reg1; reg27=reg27/reg29; reg23=reg28*reg30;
    reg1=reg16+reg1; reg23=reg13-reg23; reg13=(*f.m).resolution*reg9; reg7=reg7*reg8; reg16=(*f.m).resolution*reg27;
    reg2=reg2*reg8; reg32=reg32/reg23; reg16=reg2+reg16; reg28=reg28/reg23; reg30=reg30/reg23;
    reg7=reg13+reg7; reg1=(*f.m).deltaT*reg1; reg19=reg19/reg23; reg2=var_inter[1]*reg32; reg13=reg14*reg30;
    reg25=var_inter[0]*reg19; T reg33=reg15*reg28; T reg34=reg7*reg1; T reg35=reg16*reg1; T reg36=reg35+reg34;
    T reg37=reg33+reg2; T reg38=reg13+reg25; T reg39=var_inter[0]*reg28; T reg40=reg15*reg19; T reg41=var_inter[1]*reg30;
    T reg42=var_inter[1]*reg15; T reg43=reg14*var_inter[0]; T reg44=reg14*reg32; T reg45=reg44+reg39; T reg46=reg13-reg40;
    T reg47=reg33-reg44; T reg48=reg36*reg38; T reg49=reg43*(*f.m).f_vol[1]; T reg50=var_inter[1]*var_inter[0]; T reg51=reg14*reg15;
    T reg52=reg2-reg39; T reg53=reg25-reg41; T reg54=reg36*reg37; T reg55=reg42*(*f.m).f_vol[0]; T reg56=reg40+reg41;
    T reg57=reg43*(*f.m).f_vol[0]; T reg58=reg42*(*f.m).f_vol[1]; T reg59=reg51*(*f.m).f_vol[0]; T reg60=reg51*(*f.m).f_vol[1]; T reg61=reg50*(*f.m).f_vol[0];
    T reg62=reg50*(*f.m).f_vol[1]; T reg63=reg54-reg55; T reg64=reg36*reg56; T reg65=reg36*reg47; T reg66=reg36*reg46;
    T reg67=reg36*reg45; T reg68=reg48-reg49; T reg69=reg36*reg53; T reg70=reg36*reg52; T reg71=reg69+reg62;
    reg63=reg23*reg63; T reg72=reg64+reg58; T reg73=reg65+reg59; T reg74=reg70+reg61; T reg75=reg66+reg60;
    reg68=reg23*reg68; T reg76=reg67+reg57; reg68=ponderation*reg68; T reg77=reg23*reg71; reg63=ponderation*reg63;
    T reg78=reg23*reg76; T reg79=reg23*reg72; T reg80=reg23*reg74; T reg81=reg23*reg73; T reg82=reg23*reg75;
    T reg83=ponderation*reg79; sollicitation[indices[3]+1]+=reg83; T reg84=ponderation*reg80; sollicitation[indices[2]+0]+=reg84; sollicitation[indices[3]+0]+=-reg63;
    reg63=ponderation*reg78; sollicitation[indices[1]+0]+=reg63; sollicitation[indices[1]+1]+=-reg68; reg68=ponderation*reg77; sollicitation[indices[2]+1]+=reg68;
    T reg85=ponderation*reg81; sollicitation[indices[0]+0]+=reg85; T reg86=ponderation*reg82; sollicitation[indices[0]+1]+=reg86;
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
    T reg0=1+(*f.m).poisson_ratio; reg0=reg0/(*f.m).elastic_modulus; T reg1=pow(reg0,2); T reg2=reg0*reg1; T reg3=1.0/(*f.m).elastic_modulus;
    T reg4=(*f.m).poisson_ratio/(*f.m).elastic_modulus; T reg5=reg4*reg1; T reg6=reg3*reg2; reg1=reg3*reg1; reg2=reg4*reg2;
    T reg7=reg3*reg6; T reg8=reg4*reg1; T reg9=reg4*reg2; reg6=reg4*reg6; T reg10=reg4*reg5;
    reg1=reg3*reg1; reg5=reg3*reg5; reg2=reg3*reg2; reg1=reg1-reg10; reg8=reg10+reg8;
    reg6=reg9+reg6; reg7=reg7-reg9; T reg11=reg10+reg5; reg2=reg9+reg2; reg8=reg4*reg8;
    reg9=reg3*reg7; T reg12=reg4*reg6; reg1=reg3*reg1; T reg13=reg4*reg2; reg12=reg9-reg12;
    reg8=reg1-reg8; reg11=reg4*reg11; reg1=1-var_inter[0]; reg9=1-var_inter[1]; reg13=reg12-reg13;
    reg11=reg8-reg11; reg7=reg7/reg13; reg8=elem.pos(1)[1]*var_inter[0]; reg12=reg1*elem.pos(0)[1]; reg6=reg6/reg13;
    T reg14=reg1*elem.pos(0)[0]; T reg15=elem.pos(1)[0]*var_inter[0]; T reg16=reg9*elem.pos(1)[1]; T reg17=reg9*elem.pos(0)[1]; T reg18=reg9*elem.pos(0)[0];
    T reg19=reg9*elem.pos(1)[0]; reg11=reg11/reg13; T reg20=elem.pos(2)[1]*var_inter[0]; T reg21=reg12+reg8; T reg22=reg6*reg11;
    T reg23=reg7*reg11; T reg24=var_inter[1]*elem.pos(2)[0]; reg19=reg19-reg18; reg16=reg16-reg17; T reg25=var_inter[1]*elem.pos(2)[1];
    T reg26=elem.pos(2)[0]*var_inter[0]; T reg27=reg15+reg14; T reg28=reg7*reg23; T reg29=reg6*reg22; T reg30=reg7*(*f.m).alpha;
    T reg31=reg6*(*f.m).alpha; reg13=reg2/reg13; reg2=elem.pos(3)[0]*reg1; reg26=reg26-reg27; T reg32=var_inter[1]*elem.pos(3)[1];
    reg25=reg16+reg25; reg16=var_inter[1]*elem.pos(3)[0]; reg19=reg24+reg19; reg20=reg20-reg21; reg24=elem.pos(3)[1]*reg1;
    reg29=reg28-reg29; reg19=reg19-reg16; reg25=reg25-reg32; reg24=reg20+reg24; reg2=reg26+reg2;
    reg31=reg30+reg31; reg13=(*f.m).alpha*reg13; reg20=reg9*vectors[0][indices[0]+0]; reg26=reg9*vectors[0][indices[0]+1]; reg28=reg1*vectors[0][indices[0]+0];
    reg30=var_inter[0]*vectors[0][indices[1]+0]; T reg33=reg9*vectors[0][indices[1]+1]; T reg34=reg1*vectors[0][indices[0]+1]; T reg35=var_inter[0]*vectors[0][indices[1]+1]; T reg36=reg9*vectors[0][indices[1]+0];
    T reg37=var_inter[0]*vectors[0][indices[2]+1]; reg35=reg34+reg35; reg23=reg23/reg29; reg29=reg22/reg29; reg22=var_inter[1]*vectors[0][indices[2]+1];
    reg34=var_inter[0]*vectors[0][indices[2]+0]; reg26=reg33-reg26; reg13=reg31+reg13; reg31=reg24*reg19; reg33=var_inter[1]*vectors[0][indices[2]+0];
    T reg38=reg25*reg2; reg28=reg30+reg28; reg20=reg36-reg20; reg22=reg26+reg22; reg26=var_inter[1]*vectors[0][indices[3]+1];
    reg38=reg31-reg38; reg28=reg34-reg28; reg23=reg23*reg13; reg13=reg29*reg13; reg35=reg37-reg35;
    reg29=reg1*vectors[0][indices[3]+0]; reg30=var_inter[1]*vectors[0][indices[3]+0]; reg31=reg3*reg0; reg33=reg20+reg33; reg20=reg4*reg0;
    reg34=reg1*vectors[0][indices[3]+1]; reg35=reg34+reg35; reg13=reg23-reg13; reg30=reg33-reg30; reg23=1-(*f.m).resolution;
    reg29=reg28+reg29; reg25=reg25/reg38; reg28=reg3*reg31; reg24=reg24/reg38; reg26=reg22-reg26;
    reg2=reg2/reg38; reg22=reg4*reg20; reg33=pow(reg3,2); reg19=reg19/reg38; reg34=pow(reg4,2);
    reg22=reg28-reg22; reg13=reg23*reg13; reg28=(*f.m).alpha*(*f.m).resolution; reg34=reg33-reg34; reg33=reg25*reg35;
    reg36=reg19*reg29; reg37=reg2*reg30; T reg39=reg24*reg26; reg31=reg31/reg22; reg26=reg2*reg26;
    reg29=reg25*reg29; reg37=reg36-reg37; reg33=reg39-reg33; reg20=reg20/reg22; reg30=reg24*reg30;
    reg22=reg34/reg22; reg13=reg28+reg13; reg35=reg19*reg35; reg33=reg37+reg33; reg26=reg35-reg26;
    reg29=reg30-reg29; reg13=(*f.m).deltaT*reg13; reg28=(*f.m).resolution*reg22; reg30=(*f.m).resolution*reg20; reg34=(*f.m).resolution*reg31;
    reg7=reg7*reg23; reg6=reg6*reg23; reg23=reg11*reg23; reg11=var_inter[1]*reg2; reg35=var_inter[1]*reg24;
    reg36=var_inter[0]*reg19; reg37=var_inter[0]*reg25; reg39=reg9*reg2; T reg40=reg1*reg19; reg33=0.5*reg33;
    T reg41=reg1*reg25; reg26=reg26-reg13; reg34=reg7+reg34; reg6=reg30+reg6; reg23=reg28+reg23;
    reg7=reg9*reg24; reg29=reg29-reg13; reg28=reg36-reg11; reg30=reg35-reg37; T reg42=reg39+reg36;
    T reg43=reg41+reg35; T reg44=reg7+reg37; T reg45=reg39-reg40; T reg46=reg41-reg7; T reg47=reg40+reg11;
    T reg48=reg34*reg29; reg33=reg23*reg33; T reg49=reg6*reg26; reg29=reg6*reg29; reg26=reg34*reg26;
    T reg50=0.5*reg47; T reg51=0.5*reg30; T reg52=0.5*reg28; reg49=reg48+reg49; reg48=0.5*reg44;
    T reg53=0.5*reg43; T reg54=0.5*reg42; reg26=reg29+reg26; reg29=0.5*reg46; reg33=2*reg33;
    T reg55=0.5*reg45; T reg56=reg49*reg43; T reg57=reg49*reg46; T reg58=reg33*reg50; T reg59=reg33*reg51;
    T reg60=reg26*reg28; T reg61=reg33*reg52; T reg62=reg9*reg1; T reg63=reg9*var_inter[0]; T reg64=reg49*reg30;
    T reg65=reg33*reg55; T reg66=reg33*reg48; T reg67=reg26*reg42; T reg68=reg26*reg47; T reg69=reg26*reg45;
    T reg70=reg33*reg29; T reg71=reg33*reg53; T reg72=var_inter[1]*var_inter[0]; T reg73=var_inter[1]*reg1; T reg74=reg49*reg44;
    T reg75=reg33*reg54; T reg76=reg72*(*f.m).f_vol[1]; reg59=reg60+reg59; reg68=reg68-reg71; reg60=reg63*(*f.m).f_vol[0];
    reg58=reg58-reg56; T reg77=reg73*(*f.m).f_vol[0]; T reg78=reg72*(*f.m).f_vol[0]; reg61=reg64+reg61; reg74=reg74-reg75;
    reg64=reg62*(*f.m).f_vol[0]; reg66=reg66-reg67; T reg79=reg73*(*f.m).f_vol[1]; T reg80=reg62*(*f.m).f_vol[1]; reg65=reg57+reg65;
    reg70=reg69+reg70; reg57=reg63*(*f.m).f_vol[1]; reg70=reg70-reg80; reg59=reg59-reg76; reg65=reg65-reg64;
    reg74=reg74-reg60; reg61=reg61-reg78; reg68=reg68-reg79; reg58=reg58-reg77; reg66=reg66-reg57;
    reg68=reg38*reg68; reg58=reg38*reg58; reg59=reg38*reg59; reg61=reg38*reg61; reg66=reg38*reg66;
    reg74=reg38*reg74; reg70=reg38*reg70; reg65=reg38*reg65; sollicitation[indices[2]+1]+=ponderation*reg59; sollicitation[indices[2]+0]+=ponderation*reg61;
    sollicitation[indices[1]+1]+=ponderation*reg66; sollicitation[indices[3]+0]+=ponderation*reg58; sollicitation[indices[1]+0]+=ponderation*reg74; sollicitation[indices[0]+1]+=ponderation*reg70; sollicitation[indices[0]+0]+=ponderation*reg65;
    sollicitation[indices[3]+1]+=ponderation*reg68;
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

