#ifndef ASSIGNATION_MPI
#define ASSIGNATION_MPI
//
// C++ Interface: assignation_mpi
//
// Description:
//
//
// Author: Alain CAIGNOT <alain.caignot@lmt.ens-cachan.fr>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "../DEFINITIONS/ParallelisationData.h"
#include "../MAILLAGE/calculate_measure_G_SST.h"
#include "mpi_lmt_functions.h"
#include <fstream>
#include "crout.h"
#include "Process.h"
#include "../GEOMETRY/GeometryUser.h"
#include "../../LMT/include/containers/algo.h"

extern "C" {
    // #include "metis.h"
//     void METIS_PartGraphRecursive(int *, long long int *, long long int *, long long int *, long long int *, int *, int *, int *, int *, int *, long long int *);
    void METIS_PartGraphRecursive(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
}

using namespace LMT;


extern Crout crout;

template <class TP,class T1, class T2>
void  definition_mpi_param(TP &process,T1 &argc, T2 &argv) {
    // D�marrage de MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &process.parallelisation->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process.parallelisation->size);
}

template <class TS,class TI, class TP, class T1, class T2>
void mpi_repartition(TS &S, TI &Inter,TP &process,T1 &Stot, T1 &SubS,T2 &SubI, GeometryUser &geometry_user) {
    ///Lecture du fichier de repartion et creation des vecteurs de pointeurs de SST
    if (process.parallelisation->is_multi_cpu()) {//il s agit de creer un vecteur de vecteur qui donnera le numero des SST a parcourir par processeur
        Sc2String namein=process.structure->repertoire_des_maillages;
        Sc2String namempi="mpi_repartition";
        namein.append(namempi);
        std::ifstream f(namein.data());
        if (f) {//lecture du fichier mpi_repartition s il existe
            process.parallelisation->repartition_sst.resize(process.parallelisation->size);
            for( int i=1 ;i<process.parallelisation->size ; i++) {
                Sc2String str;
                getline(f,str);
                std::istringstream s(str);
                s >> process.parallelisation->repartition_sst[i];
            }
        } else {//repartition automatique avec les routines METIS
            Vec<int> msst,minter,mwsst, mwinter;
            msst.resize(S.size()+1);
            mwsst.resize(S.size());
            minter.resize(0);
            minter.reserve(Inter.size());
            mwinter.resize(0);
            mwinter.reserve(Inter.size());
            for( unsigned i=0 ;i<S.size() ; i++) {
                mwsst[i]=S[i].mesh.node_list_size;
                msst[i]=minter.size();
                for (unsigned j=0;j<S[i].vois.size();j++) {
                    if (S[i].vois[j] != -1) {
                        minter.push_back(S[i].vois[j]);
//                         mwinter.push_back(Inter[S[i].edge[j].internum].side[0].nodeeq.size());
                        mwinter.push_back(Inter[S[i].edge[j].internum].num);
                    }
                }
            }
            msst[S.size()]=minter.size();

            Vec<int> mrepart;
            Vec<int> mopts;
            mopts.resize(5);
            mopts.set(0);
            int nbcut,nbsst=S.size(),wgtflag=3,npart=process.parallelisation->size-1,numflag=0;
            mrepart.resize(nbsst);

            std::cout << "Proc " << process.parallelisation->rank << " : " << "msst = " << msst << std::endl;
            std::cout << "Proc " << process.parallelisation->rank << " : " << "mwsst = "  << mwsst << std::endl;
            std::cout << "Proc " << process.parallelisation->rank << " : " << "minter = "  << minter << std::endl;
            std::cout << "Proc " << process.parallelisation->rank << " : " << "mwinter = "  << mwinter << std::endl;

            METIS_PartGraphRecursive(&nbsst,msst.ptr(),minter.ptr(),mwsst.ptr(),mwinter.ptr(),&wgtflag,&numflag,&npart,mopts.ptr(),&nbcut,mrepart.ptr());
            std::cout << "Proc " << process.parallelisation->rank << " : " << "mrepart = " << mrepart << std::endl;
            process.parallelisation->repartition_sst.resize(process.parallelisation->size);
            for( unsigned i=0 ;i<mrepart.size() ; i++) {
                process.parallelisation->repartition_sst[mrepart[i]+1].push_back(i);
//                 process.parallelisation->repartition_sst[mrepart[i]].push_back(i);
            }
        }
    } else {// si on est pas en MPI on devra parcourir toutes les SST
        process.parallelisation->repartition_sst.resize(1);
        process.parallelisation->repartition_sst[0] = range(S.size());
    }
    ///On compare si on a bien toutes les SST une seule fois en faisant une norme 2
    Vec<int> temp;
    for( unsigned i=0;i<process.parallelisation->repartition_sst.size() ;i++ ) {
        temp.append(process.parallelisation->repartition_sst[i]);
    }
    sort(temp);

    if (std::abs(norm_2(1.0*temp)-norm_2(1.0*range(S.size())))> 0.000001 ) {
        std::cout << "Toutes les SST ne sont pas presentes ou le sont en double : verifier le nombre de processeur ou le fichier mpi_repartition" << std::endl;
        assert(0);
    }
    if ((unsigned)process.parallelisation->size-1 > S.size()) {
        std::cout << "Il faut avoir au moins autant de proc que de SST a repartir" << std::endl;
        assert(0);
    }



    ///On assigne le num�ro des pro contenant les SST
    for( unsigned i=0;i<process.parallelisation->repartition_sst.size() ;i++ )
        for( unsigned j=0;j<process.parallelisation->repartition_sst[i].size() ; j++)
            S[process.parallelisation->repartition_sst[i][j]].num_proc = i;

    ///On cree le vecteur de pointeur de SST
    for( unsigned j=0;j<process.parallelisation->repartition_sst[process.parallelisation->rank].size() ; j++) {
        SubS.push_back(&S[process.parallelisation->repartition_sst[process.parallelisation->rank][j]]);
    }
    for( unsigned i=0;i<S.size() ;i++ ) {
        Stot.push_back(&S[i]);
        S[i].mesh.typmat=S[i].typmat;
        S[i].mesh.numsst=S[i].num;
        S[i].mesh.num_proc=S[i].num_proc;
        S[i].mesh.id_group=S[i].id_group;
        // assignation du materiau aux elements du maillage ainsi que le numero de la sst (pour affichage)
//         apply(S[i].mesh->elem_list,apply_mat_elem(),S[i].typmat,S[i].num,S[i].num_proc);
    }
    if (not process.parallelisation->is_multi_cpu())
        SubS=Stot;


    ///Creation du vecteur de la structure d'interface a envoyer sur les pro et sur le master
        //std::cout << process.parallelisation->rank << " " << process.parallelisation->repartition_sst[process.parallelisation->rank] << std::endl;
    for(unsigned i=0;i<process.parallelisation->repartition_sst[process.parallelisation->rank].size();i++) {
        int numsst1 = process.parallelisation->repartition_sst[process.parallelisation->rank][i];
        for( unsigned j=0;j<S[numsst1].edge.size() ;j++ ) {
            unsigned side=S[numsst1].edge[j].datanum;
            unsigned num=S[numsst1].edge[j].internum;
            int numsst2;
            if (Inter[num].type != "Ext" or Inter[num].comp == "periodique") {
                if ( Inter[num].vois[0] == numsst1)
                    numsst2 = Inter[num].vois[2];
                else
                    numsst2 = Inter[num].vois[0];
                //sur quel processeur est numsst2 ?
                for( unsigned k1=0;k1 <process.parallelisation->repartition_sst.size()  ;k1++ ) {
                    if (find(process.parallelisation->repartition_sst[k1],LMT::_1==numsst2)) {
                        if (k1 == (unsigned)process.parallelisation->rank)
                            break;
                        else {
                            ParallelisationData::INTERTOEXCHANGE temp;
                            temp.num = num;
                            temp.sidetosend = side;
                            temp.to = k1;
                            process.parallelisation->intertoexchange.push_back(temp);
                            //std::cout << process.parallelisation->rank << " doit envoyer l interface "<< num << " cote " << side << " a " << k1 << std::endl;
                        }
                    }
                }
            }
            if (side==0) {
                ParallelisationData::INTERTOEXCHANGE temp;
                temp.num = num;
                temp.sidetosend = side;
                temp.to = 0;
                process.parallelisation->intertoexchangeformaster.push_back(temp);
            }
            if (side==1) {
                ParallelisationData::INTERTOEXCHANGE temp;
                temp.num = num;
                temp.sidetosend = side;
                temp.to = 0;
                process.parallelisation->intertoexchangeformaster.push_back(temp);
            }
        }
    }
    ///Creation du vecteur a envoyer classe par pro
    for( unsigned i=0;i<process.parallelisation->intertoexchange.size() ;i++ ) {
        int found=0;
        for( unsigned j=0;j<process.parallelisation->intertoexchangebypro.size() ;j++ ) {
            if (process.parallelisation->intertoexchange[i].to==process.parallelisation->intertoexchangebypro[j].to) {
                process.parallelisation->intertoexchangebypro[j].inter.push_back(process.parallelisation->intertoexchange[i]);
                found+=1;
                break;
            }
        }
        if (found==0) {
            ParallelisationData::INTERTOEXCHANGEBYPRO temp;
            temp.inter.push_back(process.parallelisation->intertoexchange[i]);
            temp.to = process.parallelisation->intertoexchange[i].to;
            process.parallelisation->intertoexchangebypro.push_back(temp);
        }
    }
    ///Le vecteur par pro peut avoir les interfaces positionn�es dans un ordre diff�rent, il faut donc les classer par ordre croissant
    for( unsigned i=0; i<process.parallelisation->intertoexchangebypro.size();i++ ) {
        Vec<int> internum;
        for(unsigned j=0 ;j< process.parallelisation->intertoexchangebypro[i].inter.size();j++ ) {
            internum.push_back(process.parallelisation->intertoexchangebypro[i].inter[j].num);
        }
        Vec<ParallelisationData::INTERTOEXCHANGE> temp;
        temp.resize(internum.size());
        temp=process.parallelisation->intertoexchangebypro[i].inter[sort_with_index(internum)];
        process.parallelisation->intertoexchangebypro[i].inter=temp;
    }
    ///Affichage de qui envoie quoi � qui
    //     for( unsigned i=0; i<process.parallelisation->intertoexchangebypro.size();i++ ) {
    //         for(unsigned j=0 ;j< process.parallelisation->intertoexchangebypro[i].inter.size();j++ ) {
    //             std::cout << process.parallelisation->rank << " envoie sur le pro " << process.parallelisation->intertoexchangebypro[i].to << " l interface " << process.parallelisation->intertoexchangebypro[i].inter[j].num << endl;
    //         }
    //     }

    //     std::cout << "Taille intertoexchangebypro " << process.parallelisation->rank << " " << process.parallelisation->intertoexchangebypro.size() << std::endl;
    //     std::cout << "Taille intertoexchange " << process.parallelisation->rank << " " << process.parallelisation->intertoexchange.size() << std::endl;
    //     std::cout << "Taille intertoexchangeformaster " << process.parallelisation->rank << " " << process.parallelisation->intertoexchangeformaster.size() << std::endl;

    ///Creation du sous vecteur d interface
    Vec<int> interput;
    for( unsigned i=0;i<SubS.size() ;i++ ) {
        for( unsigned j=0;j<SubS[i].edge.size() ;j++ ) {
            if (!find(interput,LMT::_1==SubS[i].edge[j].internum)) {
                SubI.push_back(&Inter[SubS[i].edge[j].internum]);
                interput.push_back(SubS[i].edge[j].internum);
            }
        }
    }
//     if (process.parallelisation->is_local_cpu())///utile seulement si je dois assigner qq vecteurs resultats pour les interfaces -- devrait etre obsolete
//         for( unsigned i=0;i<Inter.size() ;i++ )
//             SubI.push_back(&Inter[i]);

    for(unsigned i=0;i<process.parallelisation->repartition_sst[process.parallelisation->rank].size();i++)
        process.parallelisation->listsst.push_back(process.parallelisation->repartition_sst[process.parallelisation->rank][i]);

    for(unsigned i=0;i<process.parallelisation->repartition_sst[process.parallelisation->rank].size();i++)
        for( unsigned j=0;j<S[process.parallelisation->repartition_sst[process.parallelisation->rank][i]].edge.size() ;j++ )
            process.parallelisation->listinter.push_back(S[process.parallelisation->repartition_sst[process.parallelisation->rank][i]].edge[j].internum);



/*    std::cout << "Vecteur d interface " << process.parallelisation->rank << " " << SubI.data << std::endl;
    for( unsigned i=0;i<Inter.size() ;i++ )
        std::cout<< &Inter[i] ;
    std::cout << std::endl;*/
    ///Affichage de la repartion des SST et Interface en terme de nombre de noeuds par pro et � echanger
    if (not process.parallelisation->is_local_cpu()) {
        for( unsigned k=0;k<(unsigned) Inter.size() ;k++ ) {
            if (Inter[k].type=="Int") {
                int numsst1= Inter[k].vois[0];
                int numsst2= Inter[k].vois[2];
                for(unsigned i=0;i<(unsigned) process.parallelisation->size;i++) {
                    if ((find(process.parallelisation->repartition_sst[i],LMT::_1==numsst1)) and (!find(process.parallelisation->repartition_sst[i],LMT::_1==numsst2))) {
                        process.parallelisation->listinter.push_back(k);
                        break;
                    }
                    if ((!find(process.parallelisation->repartition_sst[i],LMT::_1==numsst1)) and (find(process.parallelisation->repartition_sst[i],LMT::_1==numsst2))) {
                        process.parallelisation->listinter.push_back(k);
                        break;
                    }
                }
            }
        }
        std::cout << "Liste des interfaces qui vont etre envoyees  : " << process.parallelisation->listinter << std::endl;
    }
    Vec<int> nodeparpro,nodeinterparpro,nodeinterparprotoexchange;
    if (process.parallelisation->is_master_cpu()) {
        nodeparpro.resize(process.parallelisation->size);
        nodeparpro.set(0);
        nodeinterparpro.resize(process.parallelisation->size);
        nodeinterparpro.set(0);
        nodeinterparprotoexchange.resize(process.parallelisation->size);
        nodeinterparprotoexchange.set(0);
        for (unsigned i=0;i<process.parallelisation->repartition_sst.size();i++) {
            for( unsigned j=0;j< process.parallelisation->repartition_sst[i].size();j++ ) {
                nodeparpro[i]+=S[process.parallelisation->repartition_sst[i][j]].mesh.node_list_size;
                unsigned numsst=process.parallelisation->repartition_sst[i][j];
                for( unsigned k=0;k<S[numsst].edge.size() ;k++ ) {
                    unsigned internum=S[numsst].edge[k].internum;
//                     unsigned datanum=S[numsst].edge[k].datanum;
//                     nodeinterparpro[i]+=Inter[internum].side[datanum].nodeeq.size();
                    nodeinterparpro[i]+=Inter[internum].num;
                    if (find(process.parallelisation->listinter,LMT::_1==internum))
                        nodeinterparprotoexchange[i]+=Inter[internum].num;
//                     nodeinterparprotoexchange[i]+=Inter[internum].side[datanum].nodeeq.size();
                }
            }
        }
        std::cout << "Nb noeud par processeur                      : " << nodeparpro << std::endl;
        std::cout << "Nb noeud interface equivalent par processeur : " << nodeinterparpro << std::endl;
        std::cout << "Nb noeud interface a envoye   par processeur : " << nodeinterparprotoexchange << std::endl;
    }
    if (not process.parallelisation->is_local_cpu()) {
        process.parallelisation->listinter.resize(0);
    }


    //repartition des groups d'elements en fonction de la repartition des sst pour mpi
    geometry_user.repartition_mpi_group_elements.resize(process.parallelisation->repartition_sst.size());
    for(unsigned i_proc=0;i_proc<geometry_user.repartition_mpi_group_elements.size();i_proc++){
        geometry_user.repartition_mpi_group_elements[i_proc].resize(process.parallelisation->repartition_sst[i_proc].size());
        for(unsigned i_sst=0;i_sst<geometry_user.repartition_mpi_group_elements[i_proc].size();i_sst++){
            int id_sst=process.parallelisation->repartition_sst[i_proc][i_sst];
            geometry_user.repartition_mpi_group_elements[i_proc][i_sst]=id_sst;
            geometry_user.find_group_elements(id_sst)->processor_rank=i_proc;
        }
    }
    

    //    ////Tri du vecteur display_fields de facon a le mettre dans le meme ordre que les donnees par noeud et par element pour que la creation des PVTU se fasse dans le bon ordre
    //     if (process.affichage->type_affichage=="Sinterieur") {
    //     typedef typename TS::template SubType<0>::T::TMESH TM;
    //     const char *names[TM::TNode::nb_params+(TM::TNode::nb_params==0)];
    //     DM::get_names<typename TM::TNode>( names );
    //     Vec<Sc2String> display_fields_temp,display_fields=process.affichage->display_fields;
    //     for(unsigned i=0;i<TM::TNode::nb_params;++i)
    //       if ( std::find(display_fields.begin(),display_fields.end(),Sc2String(names[i]))!=display_fields.end())
    //         display_fields_temp.push_back(Sc2String(names[i]));
    //
    //     Data_vtk_extract_elem<true> dve;
    //     S[0].mesh.elem_list.apply_static(dve);
    //     apply( S[0].mesh.elem_list, dve, display_fields );
    //     for(typename Data_vtk_extract_elem<true>::Map::const_iterator iter=dve.mapd.begin();iter!=dve.mapd.end();++iter)
    //       if ( std::find(display_fields.begin(),display_fields.end(),iter->first)!=display_fields.end())
    //         display_fields_temp.push_back(iter->first);
    //
    //     process.affichage->display_fields=display_fields_temp;
    //     }
    //     if (process.affichage->type_affichage=="Sbord") {
    //       typedef typename TS::template SubType<0>::T::TMESH::TSkin TM;
    //       const char *names[TM::TNode::nb_params+(TM::TNode::nb_params==0)];
    //       DM::get_names<typename TM::TNode>( names );
    //       Vec<Sc2String> display_fields_temp,display_fields=process.affichage->display_fields;
    //       for(unsigned i=0;i<TM::TNode::nb_params;++i)
    //         if ( std::find(display_fields.begin(),display_fields.end(),Sc2String(names[i]))!=display_fields.end())
    //           display_fields_temp.push_back(Sc2String(names[i]));
    //
    //       Data_vtk_extract_elem<true> dve;
    //       S[0].mesh.skin.elem_list.apply_static(dve);
    //       apply( S[0].mesh.skin.elem_list, dve, display_fields );
    //       for(typename Data_vtk_extract_elem<true>::Map::const_iterator iter=dve.mapd.begin();iter!=dve.mapd.end();++iter)
    //         if ( std::find(display_fields.begin(),display_fields.end(),iter->first)!=display_fields.end())
    //           display_fields_temp.push_back(iter->first);
    //
    //       process.affichage->display_fields=display_fields_temp;
    //     }
    //
    //    ////Creation du vecteur display_fields de facon a le mettre dans le meme ordre que les donnees par noeud et par element pour que la creation des PVTU se fasse dans le bon ordre pour les interfaces
    //     typedef typename TI::template SubType<0>::T::TMESH TM1;
    //     const char *names1[TM1::TNode::nb_params+(TM1::TNode::nb_params==0)];
    //     DM::get_names<typename TM1::TNode>( names1 );
    //     Vec<Sc2String> display_fields_temp;
    //     unsigned nb_comp[TM1::TNode::nb_params+(TM1::TNode::nb_params==0)];
    //     DM::get_nb_comp<typename TM1::TNode>( nb_comp );
    //     for(unsigned i=0;i<TM1::TNode::nb_params;++i)
    //       if ( names1[i]!="pos" and nb_comp[i] )
    //         display_fields_temp.push_back(Sc2String(names1[i]));
    //
    //     Data_vtk_extract_elem<true> dve1;
    //     Inter[0].side[0].mesh.elem_list.apply_static(dve1);
    //     Vec<Sc2String> display_fields1("all");
    //     apply( Inter[0].side[0].mesh.elem_list, dve1, display_fields1 );
    //     for(typename Data_vtk_extract_elem<true>::Map::const_iterator iter=dve1.mapd.begin();iter!=dve1.mapd.end();++iter)
    //       if ( iter->second.nb_comp )
    //         display_fields_temp.push_back(iter->first);
    //
    //     process.affichage->display_fields_inter=display_fields_temp;
    //     //std::cout << "DEBUG vecteur display_fields_inter : " << process.affichage->display_fields_inter << std::endl;

}


template<class TV1,class TV2,class TV3,class TV4,class TV5>
void memory_free(TV1 &S,TV2 &Inter,TV3 &CL,TV4 &SstMat,TV5 &InterMat,Process &process) {
    std::cout << process.parallelisation->rank << " : " << process.parallelisation->listsst.size() << std::endl;
    int nbedge=0,nbvois=0;
    for(unsigned i=0;i<S.size();i++){
        if(!find(process.parallelisation->listsst,LMT::_1==i)) {
            if(not process.parallelisation->is_master_cpu()) {
                nbedge+=S[i].edge.size();
                nbvois+=S[i].vois.size();
/*                std::cout << process.parallelisation->rank << " : Edge size " << S[i].edge.size() << std::endl;
                std::cout << process.parallelisation->rank << " : Vois size " << S[i].vois.size() << std::*endl;*/
                S[i].free();
            }
        }
    }
//    std::cout << process.parallelisation->rank << " : nbedge " << nbedge << " et nbvois " << nbvois << std::endl;
    
#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.parallelisation->rank)+" : Vidage SST : ").c_str(),1);
#endif
    
    for(unsigned i=0;i<S.size();i++){
        S[i].mesh.unload();
    }
    
#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.parallelisation->rank)+" : Vidage SST maillage : ").c_str(),1);
#endif
      
//    std::cout << process.parallelisation->rank << " : " << process.parallelisation->listinter.size() << std::endl;
    for(unsigned i=0;i<Inter.size();i++){
        if(!find(process.parallelisation->listinter,LMT::_1==i)){
            if(not process.parallelisation->is_master_cpu()) {
                Inter[i].free();
                //std::cout << process.parallelisation->rank << " : On efface l interface " << i << std::endl;
            }
        }
    }
    
#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.parallelisation->rank)+" : Vidage Inter : ").c_str(),1);
#endif
    if (process.parallelisation->is_local_cpu()){
        for(unsigned i=0;i<S.size();i++){
            S[i].free();
        }
    }
    
#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.parallelisation->rank)+" : Vidage SST proc 0 : ").c_str(),1);
#endif
    
    if (process.parallelisation->is_local_cpu()){
        for(unsigned i=0;i<Inter.size();i++){
            Inter[i].free();
        }
    }

#ifdef PRINT_ALLOC
    disp_alloc((to_string(process.parallelisation->rank)+" : Vidage Inter proc 0 : ").c_str(),1);
#endif
    
    /*
    for(unsigned i=0;i<CL.size();i++)
        CL[i].free();
    for(unsigned i=0;i<SstMat.size();i++)
        SstMat[i].free();
    for(unsigned i=0;i<InterMat.size();i++)
        InterMat[i].free();
    */
}



#endif //ENDIF ASSIGNATION_MPI
