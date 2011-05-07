#ifndef ASSIGN_DEP_H
#define ASSIGN_DEP_H

// struct mise_a_zero_old_quantities{
// 	template<class TE> void operator()(TE &e) const{
// 	e.sigma_old.set(0.);
// 	e.epsilon_old.set(0.);
// 	}
// };
// 
// struct MAZ_old_quantities_SST{
// 	template<class SST> void operator()(SST &S) const{
//    	apply(S.mesh.elem_list,mise_a_zero_old_quantities());
//     	}
// };
// 
// 
// struct update_old_quantities{
// 	template<class TE> void operator()(TE &e) const{
// 	e.sigma_old=e.sigma;
// 	e.epsilon_old=e.epsilon;
// 	}
// };

#include "assign_material.h"

struct assign_dep_contrainte{
template<class SST> void operator()(SST &S) const{
assign_material_on_element(S);
S.f->set_mesh(S.mesh.m);
S.f->get_result()=S.q;
S.f->update_variables();
S.f->call_after_solve();
//S.f.vectors[0]=S.q;
//S.f.update_variables();
//S.f.call_after_solve();
}
};


template<class SST> void assign_dep_cont_slave(SST &S,Vec<typename SST::T> &q){
   S.mesh.load();
   assign_material_on_element(S);
   S.f->set_mesh(S.mesh.m);
   S.f->get_result()=q;
   S.f->update_variables();
   S.f->call_after_solve();
}

struct calcul_sigma_epsilon_MAJ{
	template<class SST,class MULTI> void operator()(SST &S,MULTI &process) const {
	//assignation du champ q calcule au champ dep du maillage de la SST, calcul des contraintes-deformation 
	assign_dep_cont_slave(S,S.t[process.intm].q);
	//mise a jour des quantites sigma_old et epsilon_old (visco uniquement)
	//apply(S.mesh.elem_list,update_old_quantities());  
  }
};

#endif //ASSIGN_DEP_H

