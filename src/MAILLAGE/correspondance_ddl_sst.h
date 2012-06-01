#ifndef CALC_SST_CORRESP
#define CALC_SST_CORRESP


using namespace LMT;


//********************************************
// calcul des correspondance ddlbord -> ddl
//********************************************
/**\ingroup Operateurs_sst
 *  \brief Fonction permettant de retrouver la correspondance entre les ddl de bord des SST et les ddl de la SST
 
 En faisant l'intersection du maillage de bord correspondant au maillage d'interface du cote correspondant avec le maillage de la sous-structure, on determine la correspondance entre les noeuds. On verifie si l'intersection trouvee est de la taille du maillage de bord. Si tel n'est pas le cas, une erreur est envoyee (arrive souvent quand on oublie de donner les fichiers de maillages). 
 
 Connaissant les numeros des noeuds dans le maillage de la sous-structure, on determine les ddls correspondant que l'on stocke dans Sst::Edge::repddledge.
 */
template<class SST>
void Calc_SST_Correspddl(SST &S, Process &process) {
    for(unsigned j=0;j<S.edge.size();++j) {
        S.edge[j].repddledge.resize(S.edge[j].mesh->node_list.size()*DIM);
        typename IntersectionCarac<typename SST::TMESHedge::TNodeList, typename SST::TMESH::TM::TNodeList >::T inter=
            intersection_ptr(S.edge[j].mesh->node_list,S.mesh->node_list,Function<DistBetweenNodes>() < 0.0001);

        Vec<unsigned> repnode;
        repnode.resize(S.edge[j].mesh->node_list.size());

        if (inter.size() != repnode.size()) {
            std::cout<< "Calc_SST_Correspddl - attention noeud non repere - probleme de correspondance - cote " << j << endl;
            std::cout << "meshini - size node " <<S.mesh.node_list_size << " numero " << S.num << endl;
            std::cout << "intersize "<< inter.size() << " edge size "  << S.edge[j].mesh->node_list.size()<< " internum " << S.edge[j].internum << endl;
            assert(0);
        }

        for(unsigned i=0;i<inter.size();++i)
            repnode[inter[i].first->number_in_original_mesh()]=inter[i].second->number_in_original_mesh();
        for(unsigned i=0;i<repnode.size();++i)
            S.edge[j].repddledge[range(i*DIM,(i+1)*DIM)]=range(repnode[i]*DIM,(repnode[i]+1)*DIM);
    }
}

#endif //CALC_SST_CORRESP