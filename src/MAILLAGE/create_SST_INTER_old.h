#include<fstream>
#include <sstream>
#include "mesh/read_avs.h"
#include "mesh/read_geof.h"
#include "mesh/write_avs.h"
#include "mesh/read_mshadv.h"
#include "containers/algo.h"
#include <map>
// #include "definition_PARAM_AFFICHAGE.h"
// #include "affichage_mesh_SST.h"

//#include "mesh/ordering.h"
//#include "containers/evaluate_nb_cycles.h"

#include "mpi.h"

extern "C" {
    // #include "metis.h"
    void METIS_PartGraphRecursive(int *, long long int *, long long int *, long long int *, long long int *, int *, int *, int *, int *, int *, long long int *);
}

using namespace LMT;
using namespace std;

/** \defgroup Maillage_geometrie_sst Géométrie et maillage des Sous-structures
\ingroup Maillage_geometrie
*/
/** \defgroup Maillage_geometrie_inter Géométrie et maillage des Interfaces
\ingroup Maillage_geometrie
*/


/** \ingroup Maillage_geometrie_sst
\brief Création du nombre de Sst et affectation du numero du materiau par l'intermédiaire du  fichier de qualification donné dans le xml
 
Le fichier de qualification (STRUCTURE::nom_fichier_qualification_materiaux) contient un numéro pour chaque sous-structure. Ce numéro correspond à l'identificateur des matériaux donnés dans le fichier xml.
 
A cette etape, on spécifie la taille du vecteur de Sous-structures par l'intermédiaire du STRUCTURE::nb_maillage renseigné dans le xml et on assigne le numéro du matériau pour chaque sous-structure.
*/
template<class TV1>
void create_SST_typmat(STRUCTURE &structure,TV1 &S,Param &process) {
    //recherche du nom du fichier de qualification pour affectation d'un materiau
    string namequalif=structure.nom_fichier_qualification_materiaux;
    string namein=structure.repertoire_des_maillages;
    namein.append(namequalif);
    if (process.rank == 0)
        cout << "\t Fichier de qualification : " << namein << endl;
    std::ifstream f(namein.c_str());
    Vec< unsigned > num_materiau;
    //num_materiau.resize(structure.nb_maillages);//a priori >> en dessous fait un push_back !
    f>>num_materiau;
    //affectation du numero aux Ssts
    S.resize(structure.nb_maillages);
    for(unsigned i=0;i<S.size();i++) {
        S[i].typmat=num_materiau[i];
        S[i].num = i;
    }
}


/**\ingroup Maillage_geometrie_sst
\brief Lecture et assignation d'un maillage par sst
 
Le nom générique renseigné dans le fichier xml est décliné selon le nombre de Ssts du problème (ex : nom0, nom1...) et on ajoute l'extension donnée dans le fichier xml ".avs" par defaut. Le nom générique doit nécessairement ere indicé au départ par 0.
On crée ensuite une boite contenant tout le maillage (repérée par ses 2 points extrémité) (create_box_mesh()) et on applique à chaque élément le type de matériau ainsi que le numéro de la sst (cf. apply_mat_elem) (pour le post traitement).
 
*/
///Fonction generique
template<class TE, class TM, class TR>
void add_new_elem(TE &e, TM &m, TR &rep_nodes) {}
///Fonctions specialisees
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Tetra,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Tetra(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Wedge,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Wedge(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Hexa,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Hexa(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Triangle,TNB,TN,TD,NET> &e, TM &m, TR&rep_nodes) {
    m.add_element(Triangle(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Quad,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Quad(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Bar,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Bar(),DefaultBehavior(),rep_nodes.ptr() );
}

// template <class TV, class T>
// unsigned find_with_index(TV &V, T val) {
//     for(unsigned i=0;i<V.size();++i)
//       if ( (unsigned)V[i] == val )
//         return i;
//     return false;
// }

struct repart_elem_decoup {
  template<class TE, class TV1, class TT, class TM>
      void operator()(TE &e,TV1 &S, TT &temp, TM &M1) const {
        int &i = *temp.d;

        typedef typename TM::TNode TNode;
        Vec<long long int> &mrepart=*temp.a;
        //Vec<Vec<int> > nodemap=*temp.c;
        Vec<map<int, TNode * > > &nodemap= *temp.c;

        Vec<TNode *> repnodes;
        //on prend les noeuds de l element en cours et on a la map
        for(unsigned ii=0;ii<e.nb_nodes;ii++) {
          repnodes.resize(e.nb_nodes);
          typename map<int, TNode * >::const_iterator iter = nodemap[mrepart[i]].find( e.node(ii)->number_in_original_mesh() );
          //si c est pas dans la liste on ajoute le noeud
          if (iter == nodemap[mrepart[i]].end())
            nodemap[mrepart[i]][e.node(ii)->number_in_original_mesh()]=S[mrepart[i]].mesh->add_node(M1.node_list[e.node(ii)->number_in_original_mesh()].pos);
          
          repnodes[ii]=nodemap[mrepart[i]][e.node(ii)->number_in_original_mesh()];
        }
        //maintenant on ajoute l'element au maillage
        add_new_elem(e,*S[mrepart[i]].mesh.m, repnodes) ;
        i++;
      }
};

template<class TN>
struct passapply {
  Vec<map<int,TN * > > *c;
  Vec<long long int> *a;
//   Vec<Vec<int> > *c;
  int *d;
};


template<class TV1>
void create_maillage_SST(STRUCTURE &structure,TV1 &S,Param &process) {
    if (structure.nom_des_maillages != "decoup" and structure.nom_des_maillages != "decoup_auto") {
        for(unsigned i=0;i<S.size();++i) {
            // lecture du nom du maillage volumique de SST
            ostringstream ss;
            ss << structure.repertoire_des_maillages << structure.nom_des_maillages << i << structure.extension;
            string namein(ss.str());
            if (process.rank == 0)
                cout <<"\t " <<  namein << endl;

            //ouverture fichier correspondant
            std::ifstream f(namein.c_str());
            if (structure.extension==".avs") {
                S[i].mesh.name=namein;
                if (process.rank == 0) S[i].mesh.load();

            } else if (structure.extension==".geof") {
                S[i].mesh.name=namein;
                if (process.rank == 0) S[i].mesh.load();

            } else{
                if (process.rank == 0)
                    cout << "Extension non reconnue : modifier create_maillage_SST" << endl;
                assert(0);
            }

            //creation de la boite englobant le maillage de la Sst (utile pour la suite)
            if (process.rank == 0) S[i].box=create_box_mesh(*S[i].mesh.m);
            if (process.size > 1)
            for(unsigned j=0 ;j<S[i].box.size() ;j++ ){
                MPI_Bcast(S[i].box[j].ptr(),S[i].box[j].size(),MPI_DOUBLE, 0,MPI_COMM_WORLD);
            }
            if (process.rank == 0) S[i].mesh.unload();

        }
    } else if (structure.nom_des_maillages == "decoup") {//ben faut faire le decoupage du maillage propose
        //lecture du maillage propose
        typename TV1::template SubType<0>::T::TMESH::TM M1;
        ostringstream ss;
        ss << structure.repertoire_des_maillages << structure.nom_des_maillages << "0" << structure.extension;
        string namein(ss.str());
        if (process.rank == 0)
            cout <<"\t " <<  namein << endl;

        //ouverture fichier correspondant
        std::ifstream f(namein.c_str());
        if (structure.extension==".avs") {
            //lecture du fichier avs
            read_avs(M1,f);
        } else if (structure.extension==".mshadv") {
            //lecture du fichier msh
          read_mshadv(M1,f);
        } else if (structure.extension==".geof") {
            //lecture du fichier msh
            read_geof(M1,f);
        } else{
            if (process.rank == 0)
                cout << "Extension non reconnue : modifier create_maillage_SST" << endl;
            assert(0);
        }
        cout << "Fin lecture maillage" << endl;
        //creation de la structure qui va bien pour metis
        M1.update_elem_neighbours();
        cout << "Fin update neighbourgs" << endl;

        Vec<long long int> melem,melemvoisin,melemvoisinok;
        melem.resize(M1.elem_list.size()+1);
        melemvoisin.resize(0);
        melemvoisin.reserve(100*M1.elem_list.size());//pourquoi 15 ? parce que ! 100 c est mieux ..
//         for(unsigned i=0;i<M1.elem_list.size();++i) {
//           melem[i]=melemvoisin.size();
//           for( SimpleConstIterator<typename TV1::template SubType<0>::T::TMESH::TM::EA *> iter=M1.get_elem_neighbours( (const typename TV1::template SubType<0>::T::TMESH::TM::EA *) M1.elem_list[i] ); iter; ++iter ){
//             melemvoisin.push_back((*iter)->number);}
//           
//         }
        
        std::cout << "Probleme de compilation : revoir LMTpp sur zone commentee dans le code " << endl;
        assert(0);
        
        M1.clear_elem_parents();
        M1.clear_elem_children();
        
        melem[M1.elem_list.size()]=melemvoisin.size();
        melemvoisinok.resize(melemvoisin.size());
        melemvoisinok=melemvoisin;
        melemvoisin.free();
        cout << "Fin generation de la liste des voisins" << endl;

        Vec<long long int> mrepart;
        Vec<int> mopts;
        mopts.resize(5);
        mopts.set(0);
        int nbcut,nbelem=M1.elem_list.size(),wgtflag=0,npart=structure.nb_maillages,numflag=0;
        mrepart.resize(nbelem);
        cout << nbelem << " " << melem.size() << " " << melemvoisinok.size() << endl;
        //decoupage du maillage
        METIS_PartGraphRecursive(&nbelem,melem.ptr(),melemvoisinok.ptr(),NULL,NULL,&wgtflag,&numflag,&npart,mopts.ptr(),&nbcut,mrepart.ptr());
        cout << "Fin METIS : " << endl;

        //posttraitement du decoupage
//         Vec<Vec<int> > repartelem;
//         repartelem.resize(npart);
        //Vec<Vec<int> > nodemap;
        Vec<map<int, typename TV1::template SubType<0>::T::TMESH::TM::TNode *> > nodemap ;
        nodemap.resize(npart);
        
        passapply<typename TV1::template SubType<0>::T::TMESH::TM::TNode> temp;
        temp.a=&mrepart;
        temp.c=&nodemap;
        int iter=0;
        temp.d=&iter;
        
        for(unsigned i=0 ;i<S.size() ;i++ ){
#ifdef PRINT_ALLOC
            total_allocated[ typeid(typename TV1::template SubType<0>::T::TMESH::TM).name() ] += sizeof(typename TV1::template SubType<0>::T::TMESH::TM);
#endif
            S[i].mesh.m=new typename TV1::template SubType<0>::T::TMESH::TM;
        }

        apply(M1.elem_list,repart_elem_decoup(),S,temp,M1);

        for(unsigned i=0;i<S.size();++i){
          ostringstream ss;
          ss << structure.repertoire_des_maillages << "decoup_auto_en_" << structure.nb_maillages << "_" << i << ".avs";
          string namein(ss.str());
          write_avs(*S[i].mesh.m,namein,Vec<string>(),Ascii());
        }

        cout << "Maillage sauvegarde" << endl;

        for(unsigned i=0;i<S.size();++i){
            S[i].box=create_box_mesh(*S[i].mesh.m);
            S[i].mesh.unload();
        }
        M1.free();
    } else if (structure.nom_des_maillages == "decoup_auto") {//ben faut faire le decoupage du maillage propose
        //lecture du maillage propose
        typename TV1::template SubType<0>::T::TMESH::TM M1;
        ostringstream ss;
        ss << structure.repertoire_des_maillages << structure.nom_des_maillages << "0" << structure.extension;
        string namein(ss.str());
        if (process.rank == 0)
            cout <<"\t " <<  namein << endl;

        //ouverture fichier correspondant
        std::ifstream f(namein.c_str());
        if (structure.extension==".avs") {
            //lecture du fichier avs
            read_avs(M1,f);
        } else if (structure.extension==".mshadv") {
            //lecture du fichier msh
            read_mshadv(M1,f);
        } else if (structure.extension==".geof") {
            //lecture du fichier msh
            read_geof(M1,f);
        }else{
            if (process.rank == 0)
                cout << "Extension non reconnue : modifier create_maillage_SST" << endl;
            assert(0);
        }
        cout << "Fin lecture maillage" << endl;
        typedef typename TV1::template SubType<0>::T::TMESH::TM::Pvec Pvec;
        Pvec mini;
        Pvec maxi;
        get_min_max(M1.node_list,ExtractDM<pos_DM>(),mini,maxi);
        
        Vec<long long int> mrepart;
        int nbelem=M1.elem_list.size(),npart=(int)(pow(structure.nb_maillages*1.0,1.0/3.0));
        mrepart.resize(nbelem);
        mrepart.set(0);
        
//         cout << npart << endl;
        
        for(unsigned i=0;i<M1.elem_list.size();i++) {
            Pvec G=center(*M1.elem_list[i]);
            for( unsigned k=0; k<mini.size();k++ ){
                mrepart[i]+=(int)((G[k]-mini[k])/(maxi[k]-mini[k])*npart)*(int)(pow(npart*1.0,k*1.0));
            }
        }
//         cout << mrepart << endl;

        Vec<map<int, typename TV1::template SubType<0>::T::TMESH::TM::TNode *> > nodemap ;
        nodemap.resize(structure.nb_maillages);
        
        passapply<typename TV1::template SubType<0>::T::TMESH::TM::TNode> temp;
        temp.a=&mrepart;
        temp.c=&nodemap;
        int iter=0;
        temp.d=&iter;
        
        for(unsigned i=0 ;i<S.size() ;i++ ){
#ifdef PRINT_ALLOC
        total_allocated[ typeid(typename TV1::template SubType<0>::T::TMESH::TM).name() ] += sizeof(typename TV1::template SubType<0>::T::TMESH::TM);
#endif
            S[i].mesh.m=new typename TV1::template SubType<0>::T::TMESH::TM;
        }

        apply(M1.elem_list,repart_elem_decoup(),S,temp,M1);
        M1.free();

        for(unsigned i=0;i<S.size();++i){
            ostringstream ss;
            if (S[i].mesh->elem_list.size() == 0) {
                erase_elem_nb2(S,i);
                i--;
            }else{
               ss << structure.repertoire_des_maillages << "decoup_auto2_en_" << structure.nb_maillages << "_" << i << ".avs";
               string namein(ss.str());
               write_avs(*S[i].mesh.m,namein,Vec<string>(),Ascii());
            }
        }
        structure.nb_maillages=S.size();
        cout << "Maillage sauvergarde PUTAIN" << endl;

        int nbelemnew=0;
        for(unsigned i=0;i<S.size();++i){
            S[i].box=create_box_mesh(*S[i].mesh.m);
            S[i].mesh.unload();
            nbelemnew+=S[i].mesh->elem_list.size();
        }
        cout << nbelem << " et " << nbelemnew << endl;
    }
}


/** \ingroup Maillage_geometrie_inter
\brief Création des interfaces comprises entre les Sst. Par défaut elles sont considérées comme parfaites.
 
Le maillage de l'interface est obtenu par intersection des maillages de peau des Ssts voisines. Après une intersection entre les boites, on sélectionne les sous-structures voisines de la sous-structure sélectionnée.
Enfin, on effectue l'intersection entre les maillages de peau des sous-structures voisines.
 
*/
template<class TV1, class TV2>
void create_perfect_interfaces(TV1 &S, TV2 &Inter, Param process) {
    if (process.size == 1) {
    for(unsigned i=0;i<S.size();i++) {
        //1ere etape : recherche des sst dont les boites s'intersectent
        Vec<unsigned> vois;
        for(unsigned j=i+1;j<S.size();j++)
            if (intersection_box(S[i].box,S[j].box)==1)
                vois.push_back(j);

        //update_elem_neighbours
//         cout << S[i].mesh.name << endl;
        S[i].mesh->update_skin();

        //2eme etape : recherche des elements communs de bord pour chacune des Sst voisines et la Sst selectionnee
        for(unsigned j=0;j<vois.size();j++) {
            unsigned num=vois[j];
            //update_elem_neighbours
//             cout << "\t" << S[num].mesh.name << endl;
            S[num].mesh->update_skin();
            //creation du nouveau maillage
            typename TV2::template SubType<0>::T::TMESH meshnew;
            intersection_meshes(S[i].mesh->skin,S[num].mesh->skin,meshnew);

            if(meshnew.elem_list.size()>0) {
                //creation d'une interface et du cote correspondant sur les Ssts
                typename TV2::template SubType<0>::T internew;
                internew.side.resize(2);
                //for(unsigned q=0;q<2;++q) {
#ifdef PRINT_ALLOC
                total_allocated[ typeid(typename TV2::template SubType<0>::T::TMESH).name() ] += sizeof(typename TV2::template SubType<0>::T::TMESH);
#endif
                internew.side[0].mesh=new typename TV2::template SubType<0>::T::TMESH;
                internew.side[0].mesh->append(meshnew);
                internew.side[0].mesh->sub_mesh(Number<1>()).elem_list.change_hash_size( *internew.side[0].mesh,1);
                internew.side[0].mesh->sub_mesh(Number<2>()).elem_list.change_hash_size( *internew.side[0].mesh,1);
                internew.side[0].mesh->sub_mesh(Number<0>()).elem_list.change_hash_size( *internew.side[0].mesh,1);
                internew.side[0].mesh->elem_list.change_hash_size( *internew.side[0].mesh,1);
                internew.side[1].mesh=internew.side[0].mesh;
                //}
                //ajout des numeros des Sst voisines et cotes correspondants
                typename TV1::template SubType<0>::T::Edge edge;
                edge.internum=Inter.size();
                edge.datanum=0;
                S[i].edge.push_back(edge);
                edge.datanum=1;
                S[num].edge.push_back(edge);
                S[i].vois.push_back(num);
                S[num].vois.push_back(i);
                internew.vois=Vec<unsigned>(i,S[i].edge.size()-1,num,S[num].edge.size()-1);
                internew.side[0].vois=Vec<unsigned>(i,S[i].edge.size()-1);
                internew.side[1].vois=Vec<unsigned>(num,S[num].edge.size()-1);
                //modification des proprietes de l'interface
                internew.type="Int";
                internew.comp="Parfait";
                internew.num=meshnew.elem_list.size();
                Inter.push_back(internew);
            }
            S[num].mesh.unload();
        }
        S[i].mesh.unload();

    }
    } else {
        //1ere etape : creation du vecteur de truc à tester
        Vec<int> totest;
        if (process.rank == 0){
            totest.reserve(S.size()*5);
            for(unsigned i=0;i<S.size();i++) {
                for(unsigned j=i+1;j<S.size();j++)
                    if (intersection_box(S[i].box,S[j].box)==1){
                    totest.push_back(i);
                    totest.push_back(j);
                    }
            }
        }
        //2eme etape : envoyer a tout le monde le vecteur d element a tester
        unsigned size=totest.size();
        MPI_Bcast(&size,1,MPI_INT, 0 ,MPI_COMM_WORLD);
        totest.resize(size);
        MPI_Bcast(totest.ptr(),totest.size(),MPI_INT, 0 ,MPI_COMM_WORLD);

        //3eme etape : chacun test un bout et stocke le resultat
        Vec<unsigned> resu_inter;
        resu_inter.reserve(3*totest.size()/2/process.size);
//         cout << process.rank << " : " << totest << endl;
        for( unsigned i=0+(unsigned)(size/2/process.size)*process.rank;i<(unsigned)(size/2/process.size)*(1+process.rank);i++ ){
//             cout << process.rank << " : " << size << " " << process.size << " " << (unsigned)(size/2/process.size) << " " << 2*i << " " << 2*i+1 << endl;
            S[totest[2*i]].mesh->update_skin();S[totest[2*i+1]].mesh->update_skin();
            typename TV2::template SubType<0>::T::TMESH meshnew;
            intersection_meshes(S[totest[2*i]].mesh->skin,S[totest[2*i+1]].mesh->skin,meshnew);
            if (meshnew.elem_list.size() > 0 ) {
                resu_inter.push_back(totest[2*i]);
                resu_inter.push_back(totest[2*i+1]);
                resu_inter.push_back(meshnew.elem_list.size());
            }
            if ((2*i+3)<totest.size())
                if (totest[2*i+1]!=totest[2*i+3])
                    S[totest[2*i]].mesh.unload();
            if ((i+1)==(unsigned)(size/2/process.size)*(1+process.rank))
                S[totest[2*i]].mesh.unload();
            S[totest[2*i+1]].mesh.unload();
        }
        if (process.rank == process.size-1) {
            for( unsigned i=(unsigned)(size/2/process.size)*(1+ process.rank);i<size/2;i++){
//                 cout << process.rank << " pouet : " << 2*i << " " << 2*i+1 << endl;
                S[totest[2*i]].mesh->update_skin();S[totest[2*i+1]].mesh->update_skin();
                typename TV2::template SubType<0>::T::TMESH meshnew;
                intersection_meshes(S[totest[2*i]].mesh->skin,S[totest[2*i+1]].mesh->skin,meshnew);
                if (meshnew.elem_list.size() > 0 ) {
                    resu_inter.push_back(totest[2*i]);
                    resu_inter.push_back(totest[2*i+1]);
                    resu_inter.push_back(meshnew.elem_list.size());
                }
                if ((2*i+3)<totest.size())
                    if (totest[2*i+1]!=totest[2*i+3])
                        S[totest[2*i]].mesh.unload();
                if ((i+1)==(unsigned)(size/2/process.size)*(1+process.rank))
                    S[totest[2*i]].mesh.unload();
                S[totest[2*i+1]].mesh.unload();
            }
        }
        //4eme etape : on rappatrie les resultats sur tout les pros
        size=0;
        int size2=resu_inter.size();
//         cout << process.rank << " : " << resu_inter << endl;
        MPI_Allreduce(&size2,&size,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        Vec<unsigned> resu_inter_complet;
        resu_inter_complet.resize(size);
//         cout << process.rank << " : " << "size : " << size << endl;
        
        Vec<int> rcount,displs1;
        rcount.resize(process.size);
        displs1.resize(process.size);
        displs1.set(0);
        MPI_Allgather(&size2,1,MPI_INT,rcount.ptr(),1,MPI_INT,MPI_COMM_WORLD);
//         cout << process.rank << " : " << rcount << endl;
        for( unsigned i=1;i<displs1.size() ;i++ ){
            displs1[i]=displs1[i-1]+rcount[i-1];
        }
//         cout << process.rank << " : " << displs1 << endl;
        MPI_Allgatherv(resu_inter.ptr(),resu_inter.size(),MPI_INT,resu_inter_complet.ptr(), rcount.ptr(),displs1.ptr(),MPI_INT,MPI_COMM_WORLD);
        resu_inter.free();
//         cout << process.rank << " : " << resu_inter_complet << endl;
        
        //penser a mettre le S.mesh.elem_list_size à jour a partir du rank 0
        for (unsigned i=0;i<S.size();i++){
            MPI_Bcast(&S[i].mesh.node_list_size,1,MPI_INT,0,MPI_COMM_WORLD);
            MPI_Bcast(&S[i].mesh.elem_list_size,1,MPI_INT,0,MPI_COMM_WORLD);
        }

        //5eme etape modification du vecteur d interfaces, creation des edges + on met la taille de l interface dans num
        Inter.resize(resu_inter_complet.size()/3);
        for( unsigned i=0;i<Inter.size() ;i++ ){
            //creation d'une interface et du cote correspondant sur les Ssts
            Inter[i].side.resize(2);
//             for(unsigned q=0;q<2;++q)
//                 internew.side[q].mesh=meshnew;
            //ajout des numeros des Sst voisines et cotes correspondants
            typename TV1::template SubType<0>::T::Edge edge;
            edge.internum=i;
            edge.datanum=0;
            S[resu_inter_complet[3*i]].edge.push_back(edge);
            edge.datanum=1;
            S[resu_inter_complet[3*i+1]].edge.push_back(edge);

            S[resu_inter_complet[3*i]].vois.push_back(resu_inter_complet[3*i+1]);
            S[resu_inter_complet[3*i+1]].vois.push_back(resu_inter_complet[3*i]);
            
            Inter[i].vois=Vec<unsigned>(resu_inter_complet[3*i],S[resu_inter_complet[3*i]].edge.size()-1,resu_inter_complet[3*i+1],S[resu_inter_complet[3*i+1]].edge.size()-1);
            Inter[i].side[0].vois=Vec<unsigned>(resu_inter_complet[3*i],S[resu_inter_complet[3*i]].edge.size()-1);
            Inter[i].side[1].vois=Vec<unsigned>(resu_inter_complet[3*i+1],S[resu_inter_complet[3*i+1]].edge.size()-1);
            Inter[i].type="Int";
            Inter[i].comp="Parfait";
            Inter[i].num=resu_inter_complet[3*i+2];
        }
    }
}

/** \ingroup Maillage_geometrie_inter
\brief Création des interfaces de jeu physique.
 
On crée les nouvelles interfaces à partir des noms et numéros donnés dans le fichier de données. 
*/
template<class TV1, class TV2>
void create_gap_interfaces(STRUCTURE &structure, TV1 &S, TV2 &Inter,Param &process) {
    for(unsigned i=0;i<structure.inter_jeu.size();i++) {
        typename TV2::template SubType<0>::T internew;
        internew.side.resize(2);

        for( unsigned j=0;j<2;j++) {
            //creation du nom de fichier
            ostringstream ss;
            ss << structure.repertoire_des_maillages << structure.nom_maillages_jeu << structure.inter_jeu[i][j] <<"-"<<structure.inter_jeu[i][0] <<"-"<< structure.inter_jeu[i][1] << structure.extension;
            string namein(ss.str());
            if (process.rank == 0)
                cout <<"\t " <<  namein << endl;

            //ouverture fichier correspondant
            std::ifstream f(namein.c_str());
            if (structure.extension==".avs") {
                //lecture du fichier avs
#ifdef PRINT_ALLOC
                total_allocated[ typeid(typename TV2::template SubType<0>::T::TMESH).name() ] += sizeof(typename TV2::template SubType<0>::T::TMESH);
#endif
                internew.side[j].mesh=new typename TV2::template SubType<0>::T::TMESH;
                read_avs(*internew.side[j].mesh,f);
                //remise à 1 de la table de hashage pour prendre moins de place...
                internew.side[j].mesh->sub_mesh(Number<1>()).elem_list.change_hash_size( *internew.side[j].mesh,1);
                internew.side[j].mesh->sub_mesh(Number<2>()).elem_list.change_hash_size( *internew.side[j].mesh,1);
                internew.side[j].mesh->sub_mesh(Number<0>()).elem_list.change_hash_size( *internew.side[j].mesh,1);
                internew.side[j].mesh->elem_list.change_hash_size( *internew.side[j].mesh,1);

            } else if (structure.extension==".geof") {
                //lecture du fichier avs
#ifdef PRINT_ALLOC
                total_allocated[ typeid(typename TV2::template SubType<0>::T::TMESH).name() ] += sizeof(typename TV2::template SubType<0>::T::TMESH);
#endif
                internew.side[j].mesh=new typename TV2::template SubType<0>::T::TMESH;
                read_geof(*internew.side[j].mesh,f);
                //remise à 1 de la table de hashage pour prendre moins de place...
                internew.side[j].mesh->sub_mesh(Number<1>()).elem_list.change_hash_size( *internew.side[j].mesh,1);
                internew.side[j].mesh->sub_mesh(Number<2>()).elem_list.change_hash_size( *internew.side[j].mesh,1);
                internew.side[j].mesh->sub_mesh(Number<0>()).elem_list.change_hash_size( *internew.side[j].mesh,1);
                internew.side[j].mesh->elem_list.change_hash_size( *internew.side[j].mesh,1);

            } else {
                if (process.rank == 0)
                    cout << "Extension non reconnue : modifier create_maillage_SST" << endl;
                assert(0);
            }
            //ajout des numeros des Sst voisines et cotes correspondants
            typename TV1::template SubType<0>
            ::T::Edge edge;
            edge.internum=Inter.size();
            edge.datanum=j;
            S[structure.inter_jeu[i][j]].edge.push_back(edge);
        }
        internew.vois=Vec<unsigned>(structure.inter_jeu[i][0],S[structure.inter_jeu[i][0]].edge.size()-1,structure.inter_jeu[i][1],S[structure.inter_jeu[i][1]].edge.size()-1);
        internew.side[0].vois=Vec<unsigned>(structure.inter_jeu[i][0],S[structure.inter_jeu[i][0]].edge.size()-1);
        internew.side[1].vois=Vec<unsigned>(structure.inter_jeu[i][1],S[structure.inter_jeu[i][1]].edge.size()-1);
        internew.type="Int";
        internew.comp="Contact_jeu_physique";
        internew.num=internew.side[0].mesh->elem_list.size();
        Inter.push_back(internew);
        S[structure.inter_jeu[i][0]].vois.push_back(structure.inter_jeu[i][1]);
        S[structure.inter_jeu[i][1]].vois.push_back(structure.inter_jeu[i][0]);
    }
}

template <class TV1>
unsigned nb_zero(TV1 V,double eps){
    unsigned zob=0;
    for(unsigned i=0;i<V.size();i++)
        if (abs(V[i])< eps)
            zob++;
    return zob;
        }
/** \ingroup Maillage_geometrie_inter
\brief Création des interfaces contenues dans une condition aux limites donnée
 
On recherche les Ssts ayant une intersection avec chacune des conditions aux limites (intersection de boites). On extrait alors du maillage de peau de chaque Ssts sélectionnée le maillage contenu dans une condition aux limites donnée. 
*/
template<class TV1, class TV2, class TV5>
void create_interfaces_CL(TV1 &S, TV2 &Inter, TV5 &CL, Param &process) {
    typedef Vec<TYPEREEL, TV1::template SubType<0>::T::dim> Pvec;

    if (process.size ==1){
    for(unsigned l=0;l<CL.size();l++) {
        //recherche des ssts ayant des pts de boite communs.
        Vec<unsigned> numS, numS1;
        if (CL[l].sst_num.size() == 0) {
            for(unsigned i=0;i<S.size();i++)
                if(intersection_box(S[i].box,CL[l].box)==1)
                    numS.push_back(i);
        } else {
            for(unsigned i=0;i<CL[l].sst_num.size();i++){
                if(intersection_box(S[CL[l].sst_num[i]].box,CL[l].box)==1)
                    numS.push_back(CL[l].sst_num[i]);
            }
        }

        if (CL[l].comp=="periodique")
            for(unsigned i=0;i<S.size();i++)
                if(intersection_box(S[i].box,CL[l].box1)==1)
                    numS1.push_back(i);
        //extraction du maillage de peau de ces ssts contenu dans la boite
        for(unsigned k=0;k<numS.size();k++) {
            unsigned i=numS[k];
            S[i].mesh->update_skin();
            //extraction et creation du maillage
            typename TV2::template SubType<0>::T::TMESH meshnew;
            intersection_mesh_box(S[i].mesh->skin,CL[l].box,meshnew);
            
            if(meshnew.elem_list.size()>0) {
                //creation d'une interface de bord
                typename TV2::template SubType<0>::T internew;
#ifdef PRINT_ALLOC
                total_allocated[ typeid(typename TV2::template SubType<0>::T::TMESH).name() ] += sizeof(typename TV2::template SubType<0>::T::TMESH);
#endif
                internew.side.resize(1);
                internew.side[0].mesh=new typename TV2::template SubType<0>::T::TMESH;
                internew.side[0].mesh->append(meshnew);
                internew.side[0].mesh->sub_mesh(Number<1>()).elem_list.change_hash_size( *internew.side[0].mesh,1);
                internew.side[0].mesh->sub_mesh(Number<2>()).elem_list.change_hash_size( *internew.side[0].mesh,1);
                internew.side[0].mesh->sub_mesh(Number<0>()).elem_list.change_hash_size( *internew.side[0].mesh,1);
                internew.side[0].mesh->elem_list.change_hash_size( *internew.side[0].mesh,1);


                if (CL[l].comp=="periodique"){//on recherche l'interface qui va bien sur les autres SST
                    internew.side.resize(2);
                    for(unsigned kk=0;kk<numS1.size();kk++){
                        typename TV2::template SubType<0>::T::TMESH meshnew1;
                        S[numS1[kk]].mesh->update_skin();
                        intersection_mesh_box(S[numS1[kk]].mesh->skin,CL[l].box1,meshnew1);
                        S[numS1[kk]].mesh.unload();
                        if(meshnew1.elem_list.size()>0) {//tester si meshnew1 et meshnew sont bien generee l un de l autre par translation a priori suivant un des axes classiques...
                            Pvec G1 = barycenter_constant_rho(meshnew);
                            Pvec G2 = barycenter_constant_rho(meshnew1);
                            Pvec G1G2 = G2-G1;
                             if (nb_zero(G1G2,1e-8) == G1.size()-1 ) {//ouais y'a bien une translation entre les deux maillages
#ifdef PRINT_ALLOC
                                total_allocated[ typeid(typename TV2::template SubType<0>::T::TMESH).name() ] += sizeof(typename TV2::template SubType<0>::T::TMESH);
#endif
                                internew.side[1].mesh=new typename TV2::template SubType<0>::T::TMESH;
                                internew.side[1].mesh->append(meshnew1);
                                internew.side[1].mesh->sub_mesh(Number<1>()).elem_list.change_hash_size( *internew.side[0].mesh,1);
                                internew.side[1].mesh->sub_mesh(Number<2>()).elem_list.change_hash_size( *internew.side[0].mesh,1);
                                internew.side[1].mesh->sub_mesh(Number<0>()).elem_list.change_hash_size( *internew.side[0].mesh,1);
                                internew.side[1].mesh->elem_list.change_hash_size( *internew.side[0].mesh,1);
                                
                                //ajout des numeros de la Sst voisine et cote correspondant
                                typename TV1::template SubType<0>::T::Edge edge;
                                edge.internum=Inter.size();
                                edge.datanum=0;
                                int numedge1=S[i].edge.size();
                                S[i].edge.push_back(edge);
                                S[i].vois.push_back(numS1[kk]);
                                edge.datanum=1;
                                int numedge2= S[numS1[kk]].edge.size();
                                S[numS1[kk]].edge.push_back(edge);
                                S[numS1[kk]].vois.push_back(i);
                                internew.vois=Vec<unsigned>(i,numedge1,kk,numedge2);
                                internew.side[0].vois=Vec<unsigned>(i,numedge1);
                                internew.side[1].vois=Vec<unsigned>(numS1[kk],numedge2);
                                 //modification des proprietes de l'interface
                                internew.type="Ext";
                                internew.comp=CL[l].comp;
                                internew.refCL=l;
                                internew.num=meshnew.elem_list.size();
                                Inter.push_back(internew);

                            }
                        }
                    }
                    if (internew.side[1].mesh==NULL) {cout << "On ne trouve pas de maillage pour la periodicite" << endl;assert(0);}
                } else {
                    //ajout des numeros de la Sst voisine et cote correspondant
                    typename TV1::template SubType<0>::T::Edge edge;
                    edge.internum=Inter.size();
                    edge.datanum=0;
                    S[i].edge.push_back(edge);
                    S[i].vois.push_back(-1);
                    internew.vois=Vec<unsigned>(i,S[i].edge.size()-1);
                    internew.side[0].vois=Vec<unsigned>(i,S[i].edge.size()-1);
                    //modification des proprietes de l'interface
                    internew.type="Ext";
                    internew.comp=CL[l].comp;
                    internew.refCL=l;
                    internew.num=meshnew.elem_list.size();
                    Inter.push_back(internew);
                }
            }
            S[i].mesh.unload();
        }
    }
    } else {
        //1ere etape : creation du vecteur de truc à tester
        Vec<int> totest;
        if (process.rank == 0){
            totest.reserve(S.size()*5);
            for(unsigned l=0;l<CL.size();l++) {
                for(unsigned i=0;i<S.size();i++) {
                    if(intersection_box(S[i].box,CL[l].box)==1) {
                        totest.push_back(l);
                        totest.push_back(i);
                    }
                }
            }
        }
        //2eme etape : envoyer a tout le monde le vecteur d element a tester
        unsigned size=totest.size();
        MPI_Bcast(&size,1,MPI_INT, 0 ,MPI_COMM_WORLD);
        totest.resize(size);
        MPI_Bcast(totest.ptr(),totest.size(),MPI_INT, 0 ,MPI_COMM_WORLD);

        //3eme etape : chacun test un bout et stocke le resultat
        Vec<unsigned> resu_inter;
        resu_inter.reserve(3*totest.size()/2/process.size);
//         cout << process.rank << " : " << totest << endl;
        for( unsigned i=0+(unsigned)(size/2/process.size)*process.rank;i<(unsigned)(size/2/process.size)*(1+process.rank);i++ ){
//             cout << process.rank << " : " << size << " " << process.size << " " << (unsigned)(size/2/process.size) << " " << 2*i << " " << 2*i+1 << endl;
            S[totest[2*i+1]].mesh->update_skin();
            typename TV2::template SubType<0>::T::TMESH meshnew;
            intersection_mesh_box(S[totest[2*i+1]].mesh->skin,CL[totest[2*i]].box,meshnew);

            if (meshnew.elem_list.size() > 0 ) {
                if (CL[totest[2*i]].comp=="periodique"){//on recherche l'interface qui va bien sur les autres SST
                    Vec<unsigned> numS1;
                    for(unsigned ii=0;ii<S.size();ii++)
                        if(intersection_box(S[ii].box,CL[totest[2*i]].box1)==1)
                            numS1.push_back(ii);

                    unsigned sizebefore=resu_inter.size();
                    for(unsigned kk=0;kk<numS1.size();kk++){
                        typename TV2::template SubType<0>::T::TMESH meshnew1;
                        S[numS1[kk]].mesh->update_skin();
                        intersection_mesh_box(S[numS1[kk]].mesh->skin,CL[totest[2*i]].box1,meshnew1);
                        S[numS1[kk]].mesh.unload();
                        if(meshnew1.elem_list.size()>0) {//tester si meshnew1 et meshnew sont bien generee l un de l autre par translation a priori suivant un des axes classiques...
                            Pvec G1 = barycenter_constant_rho(meshnew);
                            Pvec G2 = barycenter_constant_rho(meshnew1);
                            Pvec G1G2 = G2-G1;
                            if (nb_zero(G1G2,1e-8) == G1.size()-1 ) {//ouais y'a bien une translation entre les deux maillages
                                resu_inter.push_back(totest[2*i]);
                                resu_inter.push_back(totest[2*i+1]);
                                resu_inter.push_back(numS1[kk]);
//                                     resu_inter.push_back(meshnew.elem_list.size());
                            }
                        }
                    }
                    if (sizebefore==resu_inter.size()) {cout << "On ne trouve pas de maillage pour la periodicite" << endl;assert(0);}
                } else {
                    resu_inter.push_back(totest[2*i]);
                    resu_inter.push_back(totest[2*i+1]);
                    resu_inter.push_back(meshnew.elem_list.size());
                }
            }
            S[totest[2*i+1]].mesh.unload();
        }
        if (process.rank == process.size-1) {
            for( unsigned i=(unsigned)(size/2/process.size)*(1+ process.rank);i<size/2;i++){
//                 cout << process.rank << " pouet : " << 2*i << " " << 2*i+1 << endl;
                S[totest[2*i+1]].mesh->update_skin();
                typename TV2::template SubType<0>::T::TMESH meshnew;
                intersection_mesh_box(S[totest[2*i+1]].mesh->skin,CL[totest[2*i]].box,meshnew);

                if (meshnew.elem_list.size() > 0 ) {
                    if (CL[totest[2*i]].comp=="periodique"){//on recherche l'interface qui va bien sur les autres SST
                        Vec<unsigned> numS1;
                        for(unsigned ii=0;ii<S.size();ii++)
                            if(intersection_box(S[ii].box,CL[totest[2*i]].box1)==1)
                                numS1.push_back(ii);

                        unsigned sizebefore=resu_inter.size();
                        for(unsigned kk=0;kk<numS1.size();kk++){
                            typename TV2::template SubType<0>::T::TMESH meshnew1;
                            S[numS1[kk]].mesh->update_skin();
                            intersection_mesh_box(S[numS1[kk]].mesh->skin,CL[totest[2*i]].box1,meshnew1);
                            S[numS1[kk]].mesh.unload();
                            if(meshnew1.elem_list.size()>0) {//tester si meshnew1 et meshnew sont bien generee l un de l autre par translation a priori suivant un des axes classiques...
                                Pvec G1 = barycenter_constant_rho(meshnew);
                                Pvec G2 = barycenter_constant_rho(meshnew1);
                                Pvec G1G2 = G2-G1;
                                if (nb_zero(G1G2,1e-8) == G1.size()-1 ) {//ouais y'a bien une translation entre les deux maillages
                                    resu_inter.push_back(totest[2*i]);
                                    resu_inter.push_back(totest[2*i+1]);
                                    resu_inter.push_back(numS1[kk]);
//                                     resu_inter.push_back(meshnew.elem_list.size());
                                }
                            }
                        }
                        if (sizebefore==resu_inter.size()) {cout << "On ne trouve pas de maillage pour la periodicite" << endl;assert(0);}
                    } else {
                        resu_inter.push_back(totest[2*i]);
                        resu_inter.push_back(totest[2*i+1]);
                        resu_inter.push_back(meshnew.elem_list.size());
                    }
                }
                S[totest[2*i+1]].mesh.unload();
            }
        }
        //4eme etape : on rappatrie les resultats sur tout les pros
        size=0;
        int size2=resu_inter.size();
//         cout << process.rank << " : " << resu_inter << endl;
        MPI_Allreduce(&size2,&size,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        Vec<unsigned> resu_inter_complet;
        resu_inter_complet.resize(size);
//         cout << process.rank << " : " << "size : " << size << endl;
        
        Vec<int> rcount,displs1;
        rcount.resize(process.size);
        displs1.resize(process.size);
        displs1.set(0);
        MPI_Allgather(&size2,1,MPI_INT,rcount.ptr(),1,MPI_INT,MPI_COMM_WORLD);
//         cout << process.rank << " : " << rcount << endl;
        for( unsigned i=1;i<displs1.size() ;i++ ){
            displs1[i]=displs1[i-1]+rcount[i-1];
        }
//         cout << process.rank << " : " << displs1 << endl;
        MPI_Allgatherv(resu_inter.ptr(),resu_inter.size(),MPI_INT,resu_inter_complet.ptr(), rcount.ptr(),displs1.ptr(),MPI_INT,MPI_COMM_WORLD);
        resu_inter.free();
//         cout << process.rank << " : " << resu_inter_complet << endl;

        //5eme etape modification du vecteur d interfaces, creation des edges + on met la taille de l interface dans num
        unsigned sizebefore=Inter.size();
        Inter.resize(Inter.size()+resu_inter_complet.size()/3);
        for(unsigned i=0;i<Inter.size()-sizebefore ;i++){
            if (CL[resu_inter_complet[3*i]].comp != "periodique"){
                //creation d'une interface de bord
                Inter[i+sizebefore].side.resize(1);
                //ajout des numeros des Sst voisines et cotes correspondants
                typename TV1::template SubType<0>::T::Edge edge;
                edge.internum=i+sizebefore;
                edge.datanum=0;
                S[resu_inter_complet[3*i+1]].edge.push_back(edge);
                S[resu_inter_complet[3*i+1]].vois.push_back(-1);
    
                Inter[i+sizebefore].vois=Vec<unsigned>(resu_inter_complet[3*i+1],S[resu_inter_complet[3*i+1]].edge.size()-1);
                Inter[i+sizebefore].side[0].vois=Vec<unsigned>(resu_inter_complet[3*i+1],S[resu_inter_complet[3*i+1]].edge.size()-1);
                //modification des proprietes de l'interface
                Inter[i+sizebefore].type="Ext";
                Inter[i+sizebefore].comp=CL[resu_inter_complet[3*i]].comp;
                Inter[i+sizebefore].refCL=resu_inter_complet[3*i];
                Inter[i+sizebefore].num=resu_inter_complet[3*i+2];
            } else {
                //creation d'une interface de bord
                Inter[i+sizebefore].side.resize(2);
                //ajout des numeros des Sst voisines et cotes correspondants
                typename TV1::template SubType<0>::T::Edge edge;
                edge.internum=i+sizebefore;
                edge.datanum=0;
                int numedge1=S[resu_inter_complet[3*i+1]].edge.size();
                S[resu_inter_complet[3*i+1]].edge.push_back(edge);
                S[resu_inter_complet[3*i+1]].vois.push_back(resu_inter_complet[3*i+2]);
                
                edge.datanum=1;
                int numedge2=S[resu_inter_complet[3*i+2]].edge.size();
                S[resu_inter_complet[3*i+2]].edge.push_back(edge);
                S[resu_inter_complet[3*i+2]].vois.push_back(resu_inter_complet[3*i+1]);
                Inter[i+sizebefore].vois=Vec<unsigned>(resu_inter_complet[3*i+1],numedge1,resu_inter_complet[3*i+2],numedge2);
                Inter[i+sizebefore].side[0].vois=Vec<unsigned>(resu_inter_complet[3*i+1],numedge1);
                Inter[i+sizebefore].side[1].vois=Vec<unsigned>(resu_inter_complet[3*i+2],numedge2);
                 //modification des proprietes de l'interface
                Inter[i+sizebefore].type="Ext";
                Inter[i+sizebefore].comp=CL[resu_inter_complet[3*i]].comp;
                Inter[i+sizebefore].refCL=resu_inter_complet[3*i];
                Inter[i+sizebefore].num=1;

            }
        }
    }
}


/** \ingroup  Maillage_geometrie_sst
\brief Recherche de l'appartenance d'un élément à une boite définie par ses deux noeuds extrémité.
 
On recherche ici si le centre de gravité de l'élément est contenu dans la boite. et on construit la liste d'élément vérifiant ce critère.
*/
struct extract_mesh_box {
    template<class TE,class TV, class TElist>
    void operator() (TE &e, TV &box, TElist &new_elem) const {
        if (pt_in_box(center(e),box)==1)
            new_elem.push_back(&e);
    }
};

/** \ingroup   Maillage_geometrie_sst
\brief Construction d'une base a partir d'une boite (ligne en 2d, plan en 3d) contenant la normale, puis un ou deux vecteurs orthogonaux décrivant la boite.
 
Cette procédure fonctionne pour tous les cas possibles de segments en 2d. En 3d, il faudrit modifier la routine pour prendre en compte des normales autres que x, y et z... A faire
*/
template<class T>
Vec<Vec<T,2>,2> create_base_box(Vec<Vec<T,2>,2> &box) {
    typedef Vec<T,2> Pvec;
    Pvec t = (box[1]-box[0]);
    t /=norm_2(t);
    Pvec n(-t[1],t[0]);
    Vec<Pvec,2> base;
    base[0]=t;
    base[1]=n;
    return base;
}
template<class T>
Vec<Vec<T,3>,3> create_base_box(Vec<Vec<T,3>,2> &box) {
    typedef Vec<T,3> Pvec;
    T eps=1e-6;
    Pvec v = box[1]-box[0];
    Vec<Pvec,3> base;
    if(abs(v[0])<=eps) {//normale x
        base[0]=Pvec(0.,1.,0.);
        base[1]=Pvec(0.,0.,1.);
        base[2]=Pvec(1.,0.,0.); //normale
    } else if(abs(v[1])<=eps) {//normale y
        base[0]=Pvec(0.,0.,1.);
        base[1]=Pvec(1.,0.,0.);
        base[2]=Pvec(0.,1.,0.); //normale
    } else if(abs(v[2])<=eps) {//normale z
        base[0]=Pvec(1.,0.,0.);
        base[1]=Pvec(0.,1.,0.);
        base[2]=Pvec(0.,0.,1.); //normale
    } else {
        cout << "Plan dont la normale est differente de x, y ou z : Il est necessaire de specifier la normale dans le xml et d'adapter cette procedure en fonction" << endl;
        assert(0);
    }
    return base;
}

template<class T>
void erase_elem_nb2(Vec<T> &v, unsigned i) {
    v[i]=v.back();
    v.pop_back();
} ///< erase element number i. This procedure maintain the order

/** \ingroup  Maillage_geometrie_sst
\brief Fonction permettant de déterminer le maillage d'intersection de deux maillages donnés, très rapidement.
 
Pour chaque maillage, on construit une boite donnée par ses deux points extrémités (create_box_mesh()) et contenant tout le maillage. On détermine ensuite les éléments du maillage suivant contenus dans cette boite. On obtient ainsi une sélection d'éléments potentiellement identiques. Dans de nombreux cas la sélection est nulle (si les maillages n'ont qu'une arête commune par exemple), ce qui accélère très fortement la recherche de l'intersection.
 
Dans un second temps, on regarde pour chaque élément sélectionné quel est l'élément correspondant dans le second maillage (position des centre de gravité identique).
 
A partir des éléments communs obtenus on construit un nouveau maillage de Bar, Triangle ou Quad (create_mesh_from_list_elem()).
*/
template<class TM, class TM2>
void intersection_meshes(TM &m1,TM &m2, TM2 &m3, const double eps=1e-6) {

    //creation de la boite du maillage m2
    Vec<typename TM::Pvec,2> box2 = create_box_mesh(m2);
    //creation de la boite du maillage m1
    Vec<typename TM::Pvec,2> box1 = create_box_mesh(m1);

    Vec<unsigned> elem_list1;
    elem_list1.reserve(m1.elem_list.size());
    for(unsigned i=0;i<m1.elem_list.size();i++)
        if (pt_in_box(center(*m1.elem_list[i]),box2)==1)
            elem_list1.push_back(i);

    Vec<unsigned > elem_list2;
    elem_list2.reserve(m2.elem_list.size());
    for(unsigned i=0;i<m2.elem_list.size();i++)
        if (pt_in_box(center(*m2.elem_list[i]),box1)==1)
            elem_list2.push_back(i);

    Vec<unsigned> elem_list;
    if(elem_list1.size()>=elem_list2.size()) {
        elem_list.reserve(elem_list2.size());
        for (unsigned i=0;i<elem_list2.size();i++) {
            typename TM::Pvec g2 = center(*m2.elem_list[elem_list2[i]]);
            for( unsigned j=0;j<elem_list1.size() ;j++ ) {
                if ( (length(center(*m1.elem_list[elem_list1[j]])-g2) <= eps ) ) {
                    elem_list.push_back(elem_list2[i]);
                    erase_elem_nb2(elem_list1,j);
                    break;
                }
            }
        }
        create_mesh_from_list_elem(m2,elem_list,m3);
    } else {
        elem_list.reserve(elem_list1.size());
        for (unsigned i=0;i<elem_list1.size();i++) {
            typename TM::Pvec g1 = center(*m1.elem_list[elem_list1[i]]);
            for( unsigned j=0;j<elem_list2.size() ;j++ ) {
                if ( (length(center(*m2.elem_list[elem_list2[j]])-g1) <= eps ) ) {
                    elem_list.push_back(elem_list1[i]);
                    erase_elem_nb2(elem_list2,j);
                    break;
                }
            }
        }
        create_mesh_from_list_elem(m1,elem_list,m3);
    }

}


/**\ingroup Maillage_geometrie_inter
\brief Intersection d'un maillage de peau avec une boite definie par deux points extrémité : création d'un nouveau maillage a partir de l'intersection 
*/
template<class TM1, class TV, class TM2>
void intersection_mesh_box(TM1 &m, Vec<TV,2> &box, TM2 &m2) {
    //creation de la base associee a la boite2 :
    Vec<typename TM1::Pvec,TM1::dim> base = create_base_box(box);

    //recherche des elements contenus dans la boite
    Vec<unsigned> num_elem;
    num_elem.reserve(m.elem_list.size());
    for(unsigned i=0;i<m.elem_list.size();i++) {
        typename TM1::Pvec G;
        G = center(*m.elem_list[i]);
        if (pt_in_box(G,box,base)==1)
            num_elem.push_back(i);
    }
    create_mesh_from_list_elem(m,num_elem,m2);
}

/**\ingroup Maillage_geometrie_inter
\brief Procédure pour la creation d'un nouveau maillage à partir du maillage original et de la liste des éléments sélectionnés pour créer le nouveau maillage
*/
template<class TM1, class TM2, class TV>
void create_mesh_from_list_elem(TM1 &m,TV &num_elem,TM2 &m2) {
    // construction des noeuds du nouveau maillage commun
    std::map<unsigned,unsigned> corr_nodes;
    for(unsigned i=0;i<num_elem.size();++i) {
        for(unsigned n=0;n<m.elem_list[num_elem[i]]->nb_nodes_virtual();++n) {
            unsigned num_node = m.elem_list[num_elem[i]]->node_virtual(n)->number_in_original_mesh();
            if ( corr_nodes.find( num_node ) == corr_nodes.end() ) {
                corr_nodes[ num_node ] = m2.node_list.size();
                m2.add_node( m.elem_list[num_elem[i]]->node_virtual(n)->pos );
            }
        }
    }

    for(unsigned i=0;i<num_elem.size();++i) {

        if (m.elem_list[num_elem[i]]->nb_nodes_virtual()==2) {
            //bar
            m2.add_element( Bar(), DefaultBehavior(),
                            corr_nodes[m.elem_list[num_elem[i]]->node_virtual(0)->number_in_original_mesh()],
                            corr_nodes[m.elem_list[num_elem[i]]->node_virtual(1)->number_in_original_mesh()]
                          );
        } else if(m.elem_list[num_elem[i]]->nb_nodes_virtual()==3) {
            //triangle
            m2.add_element( Triangle(), DefaultBehavior(),
                            corr_nodes[m.elem_list[num_elem[i]]->node_virtual(0)->number_in_original_mesh()],
                            corr_nodes[m.elem_list[num_elem[i]]->node_virtual(1)->number_in_original_mesh()],
                            corr_nodes[m.elem_list[num_elem[i]]->node_virtual(2)->number_in_original_mesh()]
                          );
        } else if(m.elem_list[num_elem[i]]->nb_nodes_virtual()==4) {
            //quad
            m2.add_element( Quad(), DefaultBehavior(),
                            corr_nodes[m.elem_list[num_elem[i]]->node_virtual(0)->number_in_original_mesh()],
                            corr_nodes[m.elem_list[num_elem[i]]->node_virtual(1)->number_in_original_mesh()],
                            corr_nodes[m.elem_list[num_elem[i]]->node_virtual(2)->number_in_original_mesh()],
                            corr_nodes[m.elem_list[num_elem[i]]->node_virtual(3)->number_in_original_mesh()]
                          );

        } else {
            cout << "reperage automatique des interfaces non implemente pour ce type d'element" << endl;
            assert(0);

        }
    }
}



//intersection de maillages pour determination du maillage de peau en commun : creation d'un nouveau maillage a partir de l'intersection
template<class TM, class TM1, class TM2>
void intersection_mesh_INTER(TM &m1, TM1 &m2, TM2 &m3) {
    typename IntersectionCarac<typename TM::TElemList,typename TM1::TSkin::TElemList>::T inter = intersection_ptr( m1.elem_list, m2.skin.elem_list, EquElem() );
    // construction des noeuds du nouveau maillage commun
    std::map<unsigned,unsigned> corr_nodes;
    for(unsigned i=0;i<inter.size();++i) {
        for(unsigned n=0;n<inter[i].first->nb_nodes_virtual();++n) {
            unsigned num_node = inter[i].first->node_virtual(n)->number_in_original_mesh();
            if ( corr_nodes.find( num_node ) == corr_nodes.end() ) {
                corr_nodes[ num_node ] = m3.node_list.size();
                m3.add_node( inter[i].first->node_virtual(n)->pos );
            }
        }
    }

    for(unsigned i=0;i<inter.size();++i) {
        //std::cout << center(*inter[i].first) << " -- " << center(*inter[i].second) << std::endl;

        if (inter[i].first->nb_nodes_virtual()==2) {
            //bar
            m3.add_element( Bar(), DefaultBehavior(),
                            corr_nodes[inter[i].first->node_virtual(0)->number_in_original_mesh()],
                            corr_nodes[inter[i].first->node_virtual(1)->number_in_original_mesh()]
                          );

        } else if(inter[i].first->nb_nodes_virtual()==3) {
            //triangle
            m3.add_element( Triangle(), DefaultBehavior(),
                            corr_nodes[inter[i].first->node_virtual(0)->number_in_original_mesh()],
                            corr_nodes[inter[i].first->node_virtual(1)->number_in_original_mesh()],
                            corr_nodes[inter[i].first->node_virtual(2)->number_in_original_mesh()]
                          );

        } else if(inter[i].first->nb_nodes_virtual()==4) {
            //quad
            m3.add_element( Quad(), DefaultBehavior(),
                            corr_nodes[inter[i].first->node_virtual(0)->number_in_original_mesh()],
                            corr_nodes[inter[i].first->node_virtual(1)->number_in_original_mesh()],
                            corr_nodes[inter[i].first->node_virtual(2)->number_in_original_mesh()],
                            corr_nodes[inter[i].first->node_virtual(3)->number_in_original_mesh()]
                          );

        } else {
            cout << "reperage automatique des interfaces non implemente pour ce type d'element" << endl;
            assert(0);

        }
    }

}


struct make_interface_inter {
    template <class TV1, class TV2>
            void operator()(TV2 &SubI,TV1 &S) const {
                if (SubI.comp=="Parfait"){
                    if (SubI.side[0].mesh == NULL){
                        unsigned numi=SubI.vois[0];
                        unsigned numj=SubI.vois[2];
                   
                        S[numi].mesh->update_skin();
                        S[numj].mesh->update_skin();
#ifdef PRINT_ALLOC
                        total_allocated[ typeid(typename TV2::TMESH).name() ] += sizeof(typename TV2::TMESH);
#endif
                        SubI.side[0].mesh=new typename TV2::TMESH;
                        intersection_meshes(S[numi].mesh->skin,S[numj].mesh->skin,*SubI.side[0].mesh);
                        SubI.side[0].mesh->sub_mesh(Number<1>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
                        SubI.side[0].mesh->sub_mesh(Number<2>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
                        SubI.side[0].mesh->sub_mesh(Number<0>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
                        SubI.side[0].mesh->elem_list.change_hash_size( *SubI.side[0].mesh,1);

                        SubI.side[1].mesh=SubI.side[0].mesh;

                        S[numi].mesh.unload();
                        S[numj].mesh.unload();
                    }
                }
            }
};

struct make_interface_CL {
    template <class TV1, class TV2, class TV3>
            void operator()(TV2 &SubI,TV1 &S, TV3 &CL) const {
                if (SubI.type=="Ext"){
                    if (SubI.side[0].mesh==NULL){
                        S[SubI.vois[0]].mesh->update_skin();
#ifdef PRINT_ALLOC
                        total_allocated[ typeid(typename TV2::TMESH).name() ] += sizeof(typename TV2::TMESH);
#endif
                        SubI.side[0].mesh=new typename TV2::TMESH;
                        intersection_mesh_box(S[SubI.vois[0]].mesh->skin,CL[SubI.refCL].box,*SubI.side[0].mesh);
                        SubI.side[0].mesh->sub_mesh(Number<1>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
                        SubI.side[0].mesh->sub_mesh(Number<2>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
                        SubI.side[0].mesh->sub_mesh(Number<0>()).elem_list.change_hash_size( *SubI.side[0].mesh,1);
                        SubI.side[0].mesh->elem_list.change_hash_size( *SubI.side[0].mesh,1);
                        S[SubI.vois[0]].mesh.unload();
                    }
                }
                if (SubI.comp == "periodique"){
                    if (SubI.side[1].mesh==NULL){
                        S[SubI.vois[2]].mesh->update_skin();
#ifdef PRINT_ALLOC
                        total_allocated[ typeid(typename TV2::TMESH).name() ] += sizeof(typename TV2::TMESH);
#endif
                        SubI.side[1].mesh=new typename TV2::TMESH;
                        intersection_mesh_box(S[SubI.vois[2]].mesh->skin,CL[SubI.refCL].box1,*SubI.side[1].mesh);
                        SubI.side[1].mesh->sub_mesh(Number<1>()).elem_list.change_hash_size( *SubI.side[1].mesh,1);
                        SubI.side[1].mesh->sub_mesh(Number<2>()).elem_list.change_hash_size( *SubI.side[1].mesh,1);
                        SubI.side[1].mesh->sub_mesh(Number<0>()).elem_list.change_hash_size( *SubI.side[1].mesh,1);
                        SubI.side[1].mesh->elem_list.change_hash_size( *SubI.side[0].mesh,1);
                        S[SubI.vois[2]].mesh.unload();
                    }
                }
            }
};

struct mesh_unload {
    template <class TV1>
    void operator()(TV1 &S) const {
        S.mesh.unload();
    }
};

struct erase_useless {
    template <class TV1>
    void operator()(TV1 &Inter) const {
        if (Inter.side.size() == 0) cout << "side nul !" << endl;
        else if (Inter.side[0].mesh->elem_list.size() == 0) Inter.side.free();
    }
};
//******************************************************************************************************
// Programme permettant la construction des SST et Interfaces aveclecture de maillage issu de castem
//******************************************************************************************************
/** \ingroup Maillage_geometrie
\brief Création des Sous-structures et des Interfaces (maillages).
 
On construit ici le vecteur des sous-structures et des interfaces à partir des informations du fichier xml.
 
- Assignation du type de matériau aux sous-structures create_SST_typmat().
- Lecture du maillage des sous-structures create_maillage_SST().
- creation des interfaces parfaites create_perfect_interfaces().
- creation des interfaces avec jeu physique create_gap_interfaces().
- creation des interfaces de conditions aux limites. create_interfaces_CL().
*/
                     
#ifndef INFO_TIME
#define INFO_TIME
#endif 
#include "containers/evaluate_nb_cycles.h"

template<class TV1, class TV2, class TV3, class TV4,class TV5>
void create_SST_INTER(STRUCTURE &structure,TV1 &S,TV2 &Inter, TV5 &CL, Param &process, TV3 &Stot,TV3 &SubS,TV4 &SubI) {
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    TicToc tic1;
    if (process.rank==0) tic1.start();
#endif
    if (process.rank == 0)
        cout << "Creation de la geometrie des SST" <<endl;
    create_SST_typmat(structure,S,process);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) cout << "create_SST_typmat : " ;
    if (process.rank==0) tic1.stop();;
    if (process.rank==0) tic1.start();
#endif
    if (process.rank == 0)
        cout << "\t Lecture maillages des SST" <<endl;
    create_maillage_SST(structure,S,process);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) cout << "Creation maillage : " ;
    if (process.rank==0) tic1.stop();;
    if (process.rank==0) tic1.start();
#endif

    
//     affich_SST(S,process);
    //mise a jour des elements et noeuds de peau
//     if (process.rank == 0)
//         cout << "\t Update skin" << endl;
//     for(unsigned i=0;i<S.size();i++) {
//         S[i].mesh.update_skin();
//     }

    if (process.rank == 0)
        cout << "Creation des Interfaces" <<endl;
    if (process.rank == 0)
        cout << "\tParfaites" << endl;
    create_perfect_interfaces(S,Inter,process);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) cout << "Creation interfaces parfaites : " ;
    if (process.rank==0) tic1.stop();;
    if (process.rank==0) tic1.start();
#endif

    if (process.rank == 0)
        cout << "\tJeux physiques" << endl;
    create_gap_interfaces(structure,S,Inter,process);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) cout << "Creation interfaces jeu physique : " ;
    if (process.rank==0) tic1.stop();;
    if (process.rank==0) tic1.start();
#endif


    //create_other_interfaces(structure,S,Inter);
    if (process.rank == 0)
        cout << "\tConditions aux limites" << endl;
    create_interfaces_CL(S,Inter,CL, process);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) cout << "Creation interfaces conditions limites : " ;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) tic1.start();
#endif
    
    //make subs et subi
    //assignation_mpi...
    mpi_repartition(S, Inter,process,Stot,SubS,SubI);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) cout << "MPI Repartition : " ;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) tic1.start();
#endif

    //create mesh of interface perfect (only necessary)
    if (process.size > 1) {
        apply(SubI,make_interface_inter(),S);
        apply(SubI,make_interface_CL(),S, CL);
        apply(S,mesh_unload());
        //         apply(Inter,erase_useless());//fait planter la suite car meme si y'a pas besoin des maillages, on a besoin des voisins et chose comme ca...
    }

    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) cout << "\t Assignation des numeros aux interfaces " << endl;
    //assignation du numero de l'interface
    for(unsigned i=0;i<Inter.size();i++){
        Inter[i].num=i;
    }

};

