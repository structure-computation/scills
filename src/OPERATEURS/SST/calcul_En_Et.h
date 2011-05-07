#ifndef CALCUL_EN_ET_H
#define CALCUL_EN_ET_H

using namespace LMT;
using namespace std;

//Attention fichier generer dans : generation_auto.py
template<class SST> void calcul_En_Et(SST &S, typename SST::T &En, typename SST::T &Et){
//formulation isotrope ou viscoelastique 
if (S.f->get_name()=="elasticity_isotropy_stat_Qstat" or S.f->get_name()=="elasticity_viscosity_Qstat"){
	 En = S.mesh.elastic_modulus;
	 Et = S.mesh.elastic_modulus;
}
}
#endif //CALCUL_EN_ET_H 
