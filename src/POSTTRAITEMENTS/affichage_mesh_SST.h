#ifndef AFFICH_MESH_SST_H
#define AFFICH_MESH_SST_H

using namespace LMT;

#include "mesh/calculate_measure.h"
#include "mesh/remove_doubles.h"
#include "containers/evaluate_nb_cycles.h"

//----------------------------------------------------------------------------
//************************
// eclatement des SST
//***********************

/// modification du champ permettant une visu avec eclatement des SST
struct modif_qtrans {
    template<class TN>
    void operator()(TN &node, Point &trans) const {
        node.qtrans=node.pos-trans;
    }
};

struct eclat_SST_struct {
    template<class SST>
    void operator()(SST &S, Scalar ecl) const {
        /// calcul cdg du maillage
        S.mesh.load();
        Point G=barycenter_constant_rho(*S.mesh.m);
        /// translation du maillage
        Point trans=G*ecl;
        /// modification du champ qtrans
        apply(S.mesh->node_list,modif_qtrans(),trans);
    }
};


template<class TV1>
void eclat_SST(TV1 &S,Scalar ecl) {
    apply(S,eclat_SST_struct(),ecl);
};


struct Projection_num_num_proc_on_skin{
  template<class TE, class TMS, class TM> void operator()(TE &e, const TMS &mskin, const TM &m) const{
    const typename TM::EA *ea = mskin.get_parents_of(e)[0];
    typedef typename TM::TElemList::template SubType<0>::T TypeParent;
    const TypeParent &parent = static_cast<const TypeParent &>( *ea );
    e.numsst_skin = parent.numsst;
    e.num_proc_skin = parent.num_proc;
    e.typmat_skin = parent.typmat;
    e.id_group_skin = parent.id_group;
  }
};

//----------------------------------------------------------------------------
// *****************************
// affichage des SST
//********************************
/** \ingroup  Post_Traitement
\brief Affichage du maillage des sous-structures
 
On crée un fichier (tmp/paraview0.vtu) dans lequel il est possible de visualiser le maillage avec éclaté (qtrans) , numéro des sous-structures (numsst) et type de matériau (typmat)
*/
template<class TV1>
void affich_SST(TV1 &S,Process &process) {

    Sc2String typemail=process.affichage->type_affichage;    
    Sc2String nom_generique = process.affichage->repertoire_save +"Geometry/sst";

    int tmp=system(("mkdir -p "+process.affichage->repertoire_save+"Geometry").c_str());
    int tmp2=system(("mkdir -p "+process.affichage->repertoire_save+"Geometry/sst").c_str());
    //ecriture fichier paraview generique 
//     ostringstream sp;
//     sp<<"./tmp/paraview_"<<process.parallelisation->rank<<"_";
//     Sc2String strp(sp.str());
    
    //eclate des ssts
    double ecl=1.0;
    eclat_SST(S,ecl);

///ecriture des maillages pour chaque sst dans un repertoire séparé
        for(unsigned i=0;i<S.size();++i) {
		DisplayParaview dp;
		SstMesh::TM meshglob;
		meshglob.append(*S[i].mesh.m);
		S[i].mesh.unload();
                
                Sc2String strp;
                strp << nom_generique << "/sst_id_" << S[i].id;
                
		if (typemail=="Sinterieur") {
		  dp.add_mesh(meshglob,strp,Vec<Sc2String>("qtrans","typmat","numsst","num_proc","id_group"));
		} else if (typemail=="Sbord") {
		  meshglob.sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
		  meshglob.sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
		  meshglob.update_skin();
		  apply(meshglob.skin.elem_list,Projection_num_num_proc_on_skin(),meshglob.skin,meshglob);
		  dp.add_mesh(meshglob.skin,strp,Vec<Sc2String>("qtrans","typmat_skin","numsst_skin","num_proc_skin","id_group_skin"));
		}
// 		///modification du nom et deplacement du fichier generique
// 		ostringstream ss;
// 		ss<<nom_generique << "/proc_"<< process.parallelisation->rank << "_sst_id_"<<S[i].id<<".vtu";
// 		Sc2String namefile(ss.str());
// 		int tmp2=system(("mv "+strp+"0.vtu "+namefile).c_str());  

	}
};



#endif //AFFICH_MESH_SST_H
