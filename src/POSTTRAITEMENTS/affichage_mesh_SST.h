#ifndef AFFICH_MESH_SST_H
#define AFFICH_MESH_SST_H

using namespace LMT;
using namespace std;

#include "mesh/calculate_measure.h"
#include "mesh/remove_doubles.h"
#include "containers/evaluate_nb_cycles.h"

//----------------------------------------------------------------------------
//************************
// eclatement des SST
//***********************

// modification du champ permettant une visu avec eclatement des SST
struct modif_qtrans {
    template<class TN>
    void operator()(TN &node, typename TN::Pvec &trans) const {
        node.qtrans=node.pos-trans;
    }
};

struct eclat_SST_struct {
    template<class SST>
    void operator()(SST &S, TYPEREEL ecl) const {
        typedef Vec<TYPEREEL,DIM> TV;
        // calcul cdg du maillage
        S.mesh.load();
        TV G=barycenter_constant_rho(*S.mesh.m);
        // translation du maillage
        TV trans=G*ecl;
        // modification du champ qtrans
        apply(S.mesh->node_list,modif_qtrans(),trans);
    }
};


template<class TV1>
void eclat_SST(TV1 &S,TYPEREEL ecl) {
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
void affich_SST(TV1 &S,Param &process) {

    string typemail=process.affichage->type_affichage;    
    string nom_generique = process.affichage->repertoire_save +"results/Geometry_sst";

    int tmp=system(("mkdir -p "+process.affichage->repertoire_save+"results").c_str());

    //ecriture fichier paraview generique 
    ostringstream sp;
    sp<<"./tmp/paraview_"<<process.rank<<"_";
    string strp(sp.str());

    //eclate des ssts
    double ecl=1.0;
    eclat_SST(S,ecl);

    //lecture des maillages
    typename TV1::template SubType<0>::T::TMESH::TM meshglob;
    for(unsigned i=0;i<S.size();++i){
        meshglob.append(*S[i].mesh.m);
        S[i].mesh.unload();
    }
    //ecriture du maillage
    DisplayParaview dp;
    if (typemail=="Sinterieur") {
       dp.add_mesh(meshglob,strp,Vec<string>("qtrans","typmat","numsst","num_proc"));
    } else if (typemail=="Sbord") {
        meshglob.sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
        meshglob.sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
        meshglob.update_skin();
        apply(meshglob.skin.elem_list,Projection_num_num_proc_on_skin(),meshglob.skin,meshglob);
        dp.add_mesh(meshglob.skin,strp,Vec<string>("qtrans","typmat_skin","numsst_skin","num_proc_skin"));
    }
    

    if(process.affichage->save=="display") dp.exec();
    //modification du nom et deplacement du fichier generique
    ostringstream ss;
    ss<<nom_generique << "_proc_"<<S[0].num_proc<<".vtu";
    string namefile(ss.str());
    int tmp2=system(("mv "+strp+"0.vtu "+namefile).c_str());


};



#endif //AFFICH_MESH_SST_H
