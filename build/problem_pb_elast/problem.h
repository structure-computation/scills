#ifndef PROBLEM_pb_elast_H
#define PROBLEM_pb_elast_H
#ifndef has_formulation_elasticity_isotropy_stat_Qstat
#define has_formulation_elasticity_isotropy_stat_Qstat
#endif
#ifndef has_formulation_elasticity_orthotropy_stat_Qstat
#define has_formulation_elasticity_orthotropy_stat_Qstat
#endif
#include "mesh_carac.h"
#include "formulation/problem_ancestor.h"
namespace LMT {

template<class T,unsigned dim> class Problem_pb_elast;

class Problem_pb_elast_double_2 : public ProblemAncestor<double> {
public:
    typedef Mesh<Mesh_carac_pb_elast<double,2> > TM;
    Problem_pb_elast_double_2() {}
    Problem_pb_elast_double_2( TM &m, bool use_tim_davis=false ) {
        if ( use_tim_davis ) {
            formulation_elasticity_isotropy_stat_Qstat = new_formulation_elasticity_isotropy_stat_Qstat( Number<true>(), m );
            formulation_elasticity_orthotropy_stat_Qstat = new_formulation_elasticity_orthotropy_stat_Qstat( Number<true>(), m );
        } else {
            formulation_elasticity_isotropy_stat_Qstat = new_formulation_elasticity_isotropy_stat_Qstat( Number<false>(), m );
            formulation_elasticity_orthotropy_stat_Qstat = new_formulation_elasticity_orthotropy_stat_Qstat( Number<false>(), m );
        }
    }
    virtual unsigned nb_formulations() const { return 2; }
    virtual FormulationAncestor<double> *formulation_nb(unsigned i) {
        switch(i) {
          case 0: return formulation_elasticity_isotropy_stat_Qstat;
          case 1: return formulation_elasticity_orthotropy_stat_Qstat;
          default: return NULL;
        }
    }
    static FormulationAncestor<double> *new_formulation_elasticity_isotropy_stat_Qstat( Number<false>, TM &m );
    static FormulationAncestor<double> *new_formulation_elasticity_isotropy_stat_Qstat( Number<true >, TM &m );
    static FormulationAncestor<double> *new_formulation_elasticity_orthotropy_stat_Qstat( Number<false>, TM &m );
    static FormulationAncestor<double> *new_formulation_elasticity_orthotropy_stat_Qstat( Number<true >, TM &m );
    FormulationAncestor<double> *formulation_elasticity_isotropy_stat_Qstat;
    FormulationAncestor<double> *formulation_elasticity_orthotropy_stat_Qstat;
};
template<> class Problem_pb_elast<double,2> : public Problem_pb_elast_double_2 {
public:
    Problem_pb_elast<double,2>(TM &m,bool use_tim_davis=false):Problem_pb_elast_double_2(m,use_tim_davis) {}
};

class Problem_pb_elast_double_3 : public ProblemAncestor<double> {
public:
    typedef Mesh<Mesh_carac_pb_elast<double,3> > TM;
    Problem_pb_elast_double_3() {}
    Problem_pb_elast_double_3( TM &m, bool use_tim_davis=false ) {
        if ( use_tim_davis ) {
            formulation_elasticity_isotropy_stat_Qstat = new_formulation_elasticity_isotropy_stat_Qstat( Number<true>(), m );
            formulation_elasticity_orthotropy_stat_Qstat = new_formulation_elasticity_orthotropy_stat_Qstat( Number<true>(), m );
        } else {
            formulation_elasticity_isotropy_stat_Qstat = new_formulation_elasticity_isotropy_stat_Qstat( Number<false>(), m );
            formulation_elasticity_orthotropy_stat_Qstat = new_formulation_elasticity_orthotropy_stat_Qstat( Number<false>(), m );
        }
    }
    virtual unsigned nb_formulations() const { return 2; }
    virtual FormulationAncestor<double> *formulation_nb(unsigned i) {
        switch(i) {
          case 0: return formulation_elasticity_isotropy_stat_Qstat;
          case 1: return formulation_elasticity_orthotropy_stat_Qstat;
          default: return NULL;
        }
    }
    static FormulationAncestor<double> *new_formulation_elasticity_isotropy_stat_Qstat( Number<false>, TM &m );
    static FormulationAncestor<double> *new_formulation_elasticity_isotropy_stat_Qstat( Number<true >, TM &m );
    static FormulationAncestor<double> *new_formulation_elasticity_orthotropy_stat_Qstat( Number<false>, TM &m );
    static FormulationAncestor<double> *new_formulation_elasticity_orthotropy_stat_Qstat( Number<true >, TM &m );
    FormulationAncestor<double> *formulation_elasticity_isotropy_stat_Qstat;
    FormulationAncestor<double> *formulation_elasticity_orthotropy_stat_Qstat;
};
template<> class Problem_pb_elast<double,3> : public Problem_pb_elast_double_3 {
public:
    Problem_pb_elast<double,3>(TM &m,bool use_tim_davis=false):Problem_pb_elast_double_3(m,use_tim_davis) {}
};

} // namespace LMT
#endif // PROBLEM_pb_elast_H
