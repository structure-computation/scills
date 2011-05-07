#include "problem.h"
#include "formulation_3_double_elasticity_orthotropy_stat_Qstat.h"
namespace LMT {
FormulationAncestor<Problem_pb_elast_type0_3::T> *Problem_pb_elast_type0_3::new_formulation_elasticity_orthotropy_stat_Qstat( Number<false>, Problem_pb_elast_type0_3::TM &m ) { return new Formulation<Problem_pb_elast_type0_3::TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,double,false>(m); }
FormulationAncestor<Problem_pb_elast_type0_3::T> *Problem_pb_elast_type0_3::new_formulation_elasticity_orthotropy_stat_Qstat( Number<true >, Problem_pb_elast_type0_3::TM &m ) { return new Formulation<Problem_pb_elast_type0_3::TM,elasticity_orthotropy_stat_Qstat,DefaultBehavior,double,true >(m); }
}
