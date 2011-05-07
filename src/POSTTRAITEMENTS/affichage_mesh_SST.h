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
    void operator()(SST &S, typename SST::T ecl) const {
        typedef Vec<typename SST::T,SST::dim> TV;
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
    string save=process.affichage->save;
    string nom_generique = process.affichage->repertoire_save +typemail;

    system(("mkdir -p "+process.affichage->repertoire_save).c_str());

    ostringstream ss;
    ss<<nom_generique << "_"<<process.rank<< "_";
    string namefile(ss.str());


    double ecl=1.0;
    eclat_SST(S,ecl);
    DisplayParaview dp;
    /*    if (typemail=="Sinterieur") {
            typename TV1::template SubType<0>
            ::T::TMESH meshglob;
            for(unsigned i=0;i<S.size();++i)
                meshglob.append(S[i].mesh);
            dp.add_mesh(meshglob,namefile,Vec<string>("qtrans","typmat","numsst","num_proc"));
        } else if (typemail=="Sbord") {
            typename TV1::template SubType<0>
            ::T::TMESHedge meshglob;
            for(unsigned i=0;i<S.size();++i) {
                for(unsigned j=0;j<S[i].edge.size();++j)
                    meshglob.append(S[i].edge[j].mesh);
            }
            dp.add_mesh(meshglob,namefile,Vec<string>("qtrans","typmat","numsst","num_proc"));
        }*/
    typename TV1::template SubType<0>::T::TMESH::TM meshglob;
    for(unsigned i=0;i<S.size();++i){
        meshglob.append(*S[i].mesh.m);
        S[i].mesh.unload();
    }
    if (typemail=="Sinterieur") {
       dp.add_mesh(meshglob,namefile,Vec<string>("qtrans","typmat","numsst","num_proc"));
    } else if (typemail=="Sbord") {
/*        TicToc tic1;
        tic1.start();
        cout << process.rank << " nbnode avant : " << meshglob.node_list.size() << endl;
        meshglob.sub_mesh(Number<1>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
        meshglob.sub_mesh(Number<2>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
        remove_doubles(meshglob,1e-8, true);
        cout << process.rank << " skin : " << meshglob.skin.node_list.size() << endl;
        cout << process.rank << " Remove double : " ;
        tic1.stop();
        tic1.start();*/
        meshglob.sub_mesh(Number<1>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
        meshglob.sub_mesh(Number<2>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
        meshglob.update_skin();
/*        cout << process.rank << " Update skin : " ;
        tic1.stop();
        cout << process.rank << " nbnode apres : " << meshglob.node_list.size() << endl;*/
        apply(meshglob.skin.elem_list,Projection_num_num_proc_on_skin(),meshglob.skin,meshglob);
        dp.add_mesh(meshglob.skin,namefile,Vec<string>("qtrans","typmat_skin","numsst_skin","num_proc_skin"));
    }

    if (save=="display" and process.size==1 and process.affichage->command_file=="")//si j'ai un fichier de commande, je suppose que je ne veux pas visualiser les affichages tout de suite...bon pourquoi pas...
        dp.exec();


};



#endif //AFFICH_MESH_SST_H
