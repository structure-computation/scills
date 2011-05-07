#include "problem.h"
#include "formulation_2_double_elasticity_isotropy_stat_Qstat.h"
namespace LMT {
FormulationAncestor<Problem_pb_elast_double_2::T> *Problem_pb_elast_double_2::new_formulation_elasticity_isotropy_stat_Qstat( Number<false>, Problem_pb_elast_double_2::TM &m ) { return new Formulation<Problem_pb_elast_double_2::TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,double,false>(m); }
FormulationAncestor<Problem_pb_elast_double_2::T> *Problem_pb_elast_double_2::new_formulation_elasticity_isotropy_stat_Qstat( Number<true >, Problem_pb_elast_double_2::TM &m ) { return new Formulation<Problem_pb_elast_double_2::TM,elasticity_isotropy_stat_Qstat,DefaultBehavior,double,true >(m); }
}
