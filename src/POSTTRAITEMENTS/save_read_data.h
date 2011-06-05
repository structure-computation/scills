#include "mise_a_jour_quantites.h"

struct Projection_elements_0_skin{
  template<class TE, class TMS, class TM> void operator()(TE &e, const TMS &mskin, const TM &m) const{
    const typename TM::EA *ea = mskin.get_parents_of(e)[0];
    typedef typename TM::TElemList::template SubType<0>::T TypeParent;
    const TypeParent &parent = static_cast<const TypeParent &>( *ea );
    e.numsst_skin = parent.numsst;
    e.num_proc_skin = parent.num_proc;
    e.typmat_skin = parent.typmat;
  }
};

struct Projection_sigma_epsilon_on_skin_sst{
   template<class TE, class TMS, class TM> void operator()(TE &e, const TMS &mskin, const TM &m) const{
      const typename TM::EA *ea = mskin.get_parents_of(e)[0];
      typedef typename TM::TElemList::template SubType<0>::T TypeParent;
      const TypeParent &parent = static_cast<const TypeParent &>( *ea );
      e.sigma_skin = parent.sigma[0];
      e.sigma_mises_skin = parent.sigma_mises[0];
      e.epsilon_skin = parent.epsilon[0];
      e.numsst_skin = parent.numsst;
      e.num_proc_skin = parent.num_proc;
      e.typmat_skin = parent.typmat;

//#warning A modifier pour les comportements orthotrope
      #ifdef FORMUORTHO
        e.sigma_local_skin = parent.sigma_local[0];
      #endif
   }
};

template<class TSST>
void extract_nb_previous_nodes(TSST &S, BasicVec<int> &nb_previous_nodes) {
    for(unsigned i=0;i<S.size();i++) {
        S[i].mesh->update_skin();
        nb_previous_nodes.push_back(S[i].mesh->skin.node_list.size()+nb_previous_nodes[i]);
    }
}

template<class TSST>
void calcul_fields_on_sst(TSST &S, Param &process) {
     //assignation des deplacements a partir du deplacement au piquet de temps imic + calcul des champs a partir de ce deplacement
    if(process.nom_calcul=="incr")
        assign_dep_cont_slave(S,S.t[1].q); 
    else if(process.nom_calcul=="latin"){
        std::cout << "calcul_fields_on_sst non defini pour strategie latin pure " << std::endl; assert(0);
        //assign_dep_cont_slave(S,S.t[1].q);
    }
    else{std::cout << "Type de calcul non reconnu dans save_geometry_sst " << std::endl;assert(0);}
    S.mesh->update_skin();
    apply(S.mesh->skin.elem_list,Projection_elements_0_skin(),S.mesh->skin, *S.mesh.m);
    apply(S.mesh->skin.elem_list,Projection_sigma_epsilon_on_skin_sst(),S.mesh->skin,*S.mesh.m);
}

struct Extract_connectivities_on_element{
   template<class TE, class TS> void operator()(TE &e, TS &S, int nb_nodes_tot, map<int,int> &correspondance_number_in_original_mesh_local_mesh) const{
      for(unsigned i=0;i<e.nb_nodes;i++){
          S.mesh_connectivities[i][e.number]=correspondance_number_in_original_mesh_local_mesh[e.node(i)->number_in_original_mesh()]+nb_nodes_tot;
      }
      S.num_processor[e.number]=e.num_proc_skin;
      S.num_group[e.number]=e.numsst_skin;
      S.material[e.number]=e.typmat_skin;

   }
};

struct Extract_data_on_element{
   template<class TE, class TS> void operator()(TE &e, TS &S ) const{
        int nb_comp=DIM*(DIM+1)/2;
        for(unsigned i_comp=0;i_comp<nb_comp;i_comp++){
            S.sigma[i_comp][e.number]=e.sigma_skin[i_comp];
            S.epsilon[i_comp][e.number]=e.epsilon_skin[i_comp];
            S.sigma_mises[e.number]=e.sigma_mises_skin;
        }
   }
};


template<class TSST>
void create_hdf_geometry_data(TSST &S, Param &process, int nb_previous_nodes ) {
    //sauvegarde des coordonnées des noeuds
    int nb_nodes=S.mesh->skin.node_list.size();
    for(int d=0;d<DIM;d++){
        S.nodes[d].resize(nb_nodes);
        for(unsigned i=0;i<nb_nodes;i++){
            S.nodes[d][i]=S.mesh->skin.node_list[i].pos[d];
        }
    }
    
    //creation map pour connaitre les noeuds locaux
    map<int,int> correspondance_number_in_original_mesh_local_mesh;
    for(unsigned i=0;i<nb_nodes;i++){
          correspondance_number_in_original_mesh_local_mesh[ S.mesh->skin.node_list[i].number_in_original_mesh()]=i;
    }
    
    //sauvegarde des connectivites du maillage de peau
    //definition des tailles
    S.mesh_connectivities.resize(S.nb_nodes_by_element);
    for(int ne=0;ne<S.nb_nodes_by_element;ne++)
        S.mesh_connectivities[ne].resize(S.mesh->skin.elem_list.size());
    
    apply(S.mesh->skin.elem_list,Projection_elements_0_skin(),S.mesh->skin, *S.mesh.m);
    S.material.resize(S.mesh->skin.elem_list.size());
    S.num_processor.resize(S.mesh->skin.elem_list.size());
    S.num_group.resize(S.mesh->skin.elem_list.size());
    //extraction
    apply(S.mesh->skin.elem_list,Extract_connectivities_on_element(), S, nb_previous_nodes, correspondance_number_in_original_mesh_local_mesh);
    
}

template<class TSST>
void create_hdf_fields_data(TSST &S, Param &process ) {
    //sauvegarde des deplacements des noeuds
    int nb_nodes=S.mesh->skin.node_list.size();
    for(int d=0;d<DIM;d++){
        S.dep_nodes[d].resize(nb_nodes);
        for(unsigned i=0;i<nb_nodes;i++){
            S.dep_nodes[d][i]=S.mesh->skin.node_list[i].dep[d];
        }
    }
    
    //sauvegarde des donnees sur les elements
    //definition des tailles
    int nb_comp=DIM*(DIM+1)/2;
    for(unsigned i_comp=0;i_comp<nb_comp;i_comp++){
        S.sigma[i_comp].resize(S.mesh->skin.elem_list.size());
        S.epsilon[i_comp].resize(S.mesh->skin.elem_list.size());
        S.sigma_mises.resize(S.mesh->skin.elem_list.size());
    }
    //extraction
    apply(S.mesh->skin.elem_list,Extract_data_on_element(), S);
}


///sauvegarde de la geometrie utilise pour l'affichage des champs
template<class TSST>
void save_elements_hdf(TSST &S, Param &process, Hdf &hdf_file) {

    String name_geometry; name_geometry << "/Level_0/Geometry";
    
    String name_list ;
    name_list<< name_geometry << "/elements_1/list_" << S.num ;
    for (unsigned i_connect=0;i_connect<S.nb_nodes_by_element;i_connect++) {
        String name_connect;
        name_connect << name_list << "/c"<<i_connect;
        S.mesh_connectivities[i_connect].write_to( hdf_file, name_connect );
    }
    String name_field = name_list + "/num_proc";
    S.num_processor.write_to(hdf_file,name_field);
    name_field = name_list + "/material";
    S.material.write_to(hdf_file,name_field);        
    name_field = name_list + "/num_group";
    S.num_group.write_to(hdf_file,name_field);
    
    String type_elements;
    int pattern_id=0;
#if DIM == 2
    type_elements="Bar";
    pattern_id=0;
#else
    type_elements="Triangle";
    pattern_id=1;
#endif
    hdf_file.add_tag(name_list,"base",type_elements.c_str());
    
}


///sauvegarde de la geometrie utilise pour l'affichage des champs
template<class TSST>
void save_fields_hdf(TSST &S, Param &process, Hdf &hdf ) {

    String name_fields; name_fields << "/Level_0/Fields/pt_"<< process.temps->pt_cur;
    String name_sigma, name_epsilon ;
    name_sigma<< name_fields << "/sigma/list_" << S.num ;
    name_epsilon<< name_fields << "/epsilon/list_" << S.num ;
#if DIM==2
    BasicVec<String> tensor_fields= BasicVec<String>("/xx","/yy","/xy");
#else
    BasicVec<String> tensor_fields= BasicVec<String>("/xx","/yy","/zz","/xy","/xz","/yz");        
#endif
    for(unsigned i_comp=0;i_comp<tensor_fields.size();i_comp++){
        String name_sigma_field, name_epsilon_field;
        name_sigma_field=name_sigma+tensor_fields[i_comp];
        name_epsilon_field=name_epsilon+tensor_fields[i_comp];
        S.sigma[i_comp].write_to(hdf,name_sigma_field.c_str());
        S.epsilon[i_comp].write_to(hdf,name_epsilon_field.c_str());
    }
    String name_sigma_mises;
    name_sigma_mises <<name_fields<<"/sigma_mises/list_" << S.num ;
    S.sigma_mises.write_to(hdf,name_sigma_mises.c_str());
}

#include "utils_2.h"

template<class TSST>
void write_hdf_geometry(TSST &SubS, Param &process ) {
    //chaque processeur calcul stocke les noeuds de ces sst
    BasicVec<int> nb_previous_nodes;
    nb_previous_nodes.push_back(0);
    extract_nb_previous_nodes(SubS, nb_previous_nodes);   
    //creation des donnees hdf et sauvegarde d'un fichier pour chaque processeur
    for(unsigned i=0;i<SubS.size();i++) create_hdf_geometry_data(SubS[i],process,nb_previous_nodes[i]);
    String name_hdf ; name_hdf << process.affichage->name_hdf <<"_"<< process.rank<<".h5";
    if(FileExists(name_hdf.c_str())){ String command = "rm -rf "+name_hdf; int syst_rm=system(command.c_str());}
    Hdf hdf_file( name_hdf.c_str() );
    //ecriture des noeuds (concatenation en attendant la possibilite d'ecrire à la suite en hdf)
    process.affichage->name_geometry="/Level_0/Geometry";
    BasicVec<BasicVec<TYPE>,DIM> nodes;
    BasicVec<String> name_direction("x","y","z");
    for(unsigned d=0;d<DIM;d++) {
        nodes[d].resize(nb_previous_nodes[SubS.size()]);
        for(unsigned i=0;i<SubS.size();i++) 
            for(unsigned j=0;j<SubS[i].nodes[d].size();j++)
                nodes[d][j+nb_previous_nodes[i]]=SubS[i].nodes[d][j];
        String name_dim;  name_dim << process.affichage->name_geometry << "/local_nodes/" << name_direction[d];
        nodes[d].write_to(hdf_file,name_dim);
    }
    //ecriture des elements
    
    for(unsigned i=0;i<SubS.size();i++) {
        SubS[i].mesh->update_skin();
        save_elements_hdf(SubS[i], process, hdf_file);
    }
/*    
    String name_num_proc;
    for(unsigned i=0;i<SubS.size();i++) {name_num_proc <<process.affichage->name_geometry<<"/num_proc/list_" << SubS[i].num ; SubS[i].num_processor.write_to(hdf_file,name_num_proc.c_str());}*/

}

template<class TSST>
void write_hdf_fields(TSST &SubS, Param &process ) {
    //chaque processeur calcul stocke les noeuds de ces sst
    BasicVec<int> nb_previous_nodes;
    nb_previous_nodes.push_back(0);
    extract_nb_previous_nodes(SubS, nb_previous_nodes);   
    //ouverture d'un fichier pour chaque processeur
    String name_hdf ; name_hdf << process.affichage->name_hdf <<"_"<< process.rank<<".h5";
    Hdf hdf_file( name_hdf.c_str() );
    //ecriture des noeuds (concatenation en attendant la possibilite d'ecrire à la suite en hdf)
    process.affichage->name_geometry="/Level_0/Geometry";
    process.affichage->name_fields="/Level_0/Fields";
    //calcul des champs sur le maillage a partir de la solution et écriture des champs hdf
    for(unsigned i=0;i<SubS.size();i++){
        calcul_fields_on_sst(SubS[i],process);
        create_hdf_fields_data(SubS[i],process);
    }
    
    //concatenation des noeuds et ecriture dans le hdf
    BasicVec<BasicVec<TYPE>, DIM> dep_nodes;
    String name_displacements; name_displacements<< process.affichage->name_fields <<"/pt_"<<process.temps->pt_cur <<"/displacements";
    BasicVec<String> displacements_coor= BasicVec<String>("/x","/y","/z");
    for(unsigned d=0;d<DIM;d++) {
        dep_nodes[d].resize(nb_previous_nodes[SubS.size()]);
        for(unsigned i=0;i<SubS.size();i++) 
            for(unsigned j=0;j<SubS[i].nodes[d].size();j++)
                dep_nodes[d][j+nb_previous_nodes[i]]=SubS[i].dep_nodes[d][j];
        String name_displacement_coor; name_displacement_coor=name_displacements+displacements_coor[d];
        dep_nodes[d].write_to(hdf_file,name_displacement_coor.c_str());
    }

    //ecriture des champs par elements dans le hdf
    for(unsigned i=0;i<SubS.size();i++) {
        save_fields_hdf(SubS[i],process, hdf_file );
    }
    String name_fields ;
    name_fields<< "/Level_0/Fields/pt_"<< process.temps->pt_cur ;
    int i_step=process.temps->step_cur;
    int i_pt=process.temps->time_step[i_step].pt_cur;
    TYPE val_time=process.temps->time_step[i_step].t_ini+(i_pt+1)*process.temps->time_step[i_step].dt ;
    hdf_file.write_tag(name_fields,"time",val_time);

}


/** \ingroup  PostTraitements
\brief Sauvegarde des quantites désirées pour les interfaces
 
L'utilisateur choisi les quantites qu'il souhaite sauvegarder dans un fichier nommé save_inter. On boucle sur les interfaces et  on indique tout d'abord le nom puis on sauvegarde pour chaque coté de l'interface sélectionnée tous les pas de temps.
*/
template<class TV2, class TV1>
void save_data_inter(TV2 &Inter,TV1 &S, Param &process, Vec<string> &fields_to_save) {
    typedef TYPEREEL T;
    
    //nom du fichier de sauvegarde
    system(("mkdir -p "+process.affichage->repertoire_save).c_str());

    ifstream is;
    std::ostringstream ss1;
    ss1<<process.affichage->repertoire_save<<"save_inter_"<<process.rank;
    string name(ss1.str());
    ofstream os( name.c_str() );

    os << process.size << std::endl;
    if (process.size==1 or process.rank>0)
      for( unsigned i=0;i<S.size() ;i++ ){
      for( unsigned kk=0;kk<S[i].edge.size() ;kk++ ){
        unsigned q=S[i].edge[kk].internum;
        unsigned j=S[i].edge[kk].datanum;
//    for(unsigned q=0;q<Inter.size();q++) {
        os<< "Inter " << Inter[q].num <<std::endl;
//        for(unsigned j=0;j<Inter[q].side.size();j++) {
            os << "Side " << j << std::endl;
            for(unsigned pt=1;pt<=process.temps->nbpastemps;pt++) {
                os << "Pastemps " << pt<<std::endl;
                for(unsigned k=0;k<fields_to_save.size();k++) {
                    if(fields_to_save[k]=="Fchap") {
                        os << "Fchap" <<std::endl;
                        if(process.nom_calcul=="latin") os.write( (char *)Inter[q].side[j].t[pt].Fchap.ptr(), sizeof(T)*Inter[q].side[j].t[pt].Fchap.size() );
                        if(process.nom_calcul=="incr") os.write( (char *)Inter[q].side[j].t_post[pt].Fchap.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].Fchap.size() );
                        //if (pt == 2) std::cout << process.rank << " inter " << q << " pas de temps" << pt << " " << norm_2(Inter[q].side[j].t[pt].Fchap) << " " << Inter[q].side[j].t[pt].Fchap.size() << std::endl;
                    } else if(fields_to_save[k]=="F") {
                        os << "F" <<std::endl;
                        if(process.nom_calcul=="latin") os.write( (char *)Inter[q].side[j].t[pt].F.ptr(), sizeof(T)*Inter[q].side[j].t[pt].F.size() );
                        if(process.nom_calcul=="incr") os.write( (char *)Inter[q].side[j].t_post[pt].F.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].F.size() );
                    } else if(fields_to_save[k]=="Wchap") {
                        os << "Wchap" <<std::endl;
                        if(process.nom_calcul=="latin") os.write( (char *)Inter[q].side[j].t[pt].Wchap.ptr(), sizeof(T)*Inter[q].side[j].t[pt].Wchap.size() );
                        if(process.nom_calcul=="incr") os.write( (char *)Inter[q].side[j].t_post[pt].Wchap.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].Wchap.size() );
                    } else if(fields_to_save[k]=="Wpchap") {
                        os << "Wpchap" <<std::endl;
                        if(process.nom_calcul=="latin") os.write( (char *)Inter[q].side[j].t[pt].Wpchap.ptr(), sizeof(T)*Inter[q].side[j].t[pt].Wpchap.size() );
                        if(process.nom_calcul=="incr") os.write( (char *)Inter[q].side[j].t_post[pt].Wpchap.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].Wpchap.size() );
                    } else if(fields_to_save[k]=="W") {
                        os << "W" <<std::endl;
                        if(process.nom_calcul=="latin") os.write( (char *)Inter[q].side[j].t[pt].W.ptr(), sizeof(T)*Inter[q].side[j].t[pt].W.size() );
                        if(process.nom_calcul=="incr") os.write( (char *)Inter[q].side[j].t_post[pt].W.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].W.size() );
                    } else if(fields_to_save[k]=="Wp") {
                        os << "Wp" <<std::endl;
                        if(process.nom_calcul=="latin") os.write( (char *)Inter[q].side[j].t[pt].Wp.ptr(), sizeof(T)*Inter[q].side[j].t[pt].Wp.size() );
                        if(process.nom_calcul=="incr") os.write( (char *)Inter[q].side[j].t_post[pt].Wp.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].Wp.size() );
                    } else if(fields_to_save[k]=="WtildeM") {
                        os << "WtildeM" <<std::endl;
                        if(process.nom_calcul=="latin") os.write( (char *)Inter[q].side[j].t[pt].WtildeM.ptr(), sizeof(T)*Inter[q].side[j].t[pt].WtildeM.size() );
                        if(process.nom_calcul=="incr") os.write( (char *)Inter[q].side[j].t_post[pt].WtildeM.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].WtildeM.size() );
                    } else {
                        std::cout << "Erreur dans le choix des champs a sauvegarder " << std::endl;
                        assert(0);
                    }
                }
                os << "finchamps"<<std::endl;
            }
            os << "finPastemps" <<std::endl;
         os << "finSide" <<std::endl;
         }
    }
    os << "fin"<<std::endl;
    os.close();
}

/** \ingroup  PostTraitements
\brief Lecture des quantites désirées pour les interfaces
 
L'utilisateur choisi les quantites qu'il souhaite sauvegarder dans un fichier nommé save_inter. On boucle sur les interfaces et  on indique tout d'abord le nom puis on sauvegarde pour chaque coté de l'interface sélectionnée tous les pas de temps.
*/
template<class TV2>
void read_data_inter(TV2 &Inter, Param &process) {
    //typedef typename TV2::template SubType<0>::T T;
    typedef TYPEREEL T;
    
    
    //nom du fichier de sauvegarde
    ifstream is;
    std::ostringstream ss1;
    ss1<<process.affichage->repertoire_save<<"save_inter_"<<process.rank;
    string name(ss1.str());
    is.open( name.c_str() );


    unsigned toutlire=1;
    string str;
    if (is ) {
      getline(is,str);
      if (process.size != 1 and process.size != atoi(str.c_str()) ) {//ca va pas aller
        std::cout << "Pour lire les donnees soit on le fait monoprocesseur soit on doit avoir le meme nombre de pro que le calcul initial" << std::endl;
        assert(0);
      } else if (process.size != atoi(str.c_str()))//le master doit lire tout les fichiers
        toutlire=atoi(str.c_str());
    } else {
        std::cout << "Aucun fichier de sauvegarde a lire" << std::endl;
        assert(0);
    }
    is.close();
    
    
    for( unsigned i=0;i<toutlire ;i++ ){
      std::ostringstream ss2;
      if (toutlire==1) ss2<<process.affichage->repertoire_save<<"save_inter_"<<process.rank;
      else ss2<<process.affichage->repertoire_save<<"save_inter_"<<i;
      string name2(ss2.str());
      is.open( name2.c_str() );
      getline(is,str);

    while(is) {
        //lecture numero de l'interface
        getline(is,str);

        if(str=="fin") {
            break;
        } else {
            istringstream s(str);
            string nom; //nom = inter
            s >> nom;
            unsigned q; //numero de l'interface
            s >> q ;

            while(true) {
                //lecture numero du cote
                string str2;
                getline(is,str2);
                //std::cout << str2 << std::endl;
                if(str2=="finSide") {
                    break;
                } else {
                    istringstream s2(str2);
                    string nom2; //nom = Side
                    s2 >> nom2;
                    unsigned j; //numero du cote
                    s2 >> j ;

                    while(true) {
                        //lecture du pas de temps
                        string str3;
                        getline(is,str3);
                        if(str3=="finPastemps") {
                            break;
                        } else {
                            istringstream s3(str3);
                            string nom3; //nom = Pastemps
                            s3 >> nom3;
                            unsigned pt; //numero du pas de temps
                            s3 >> pt ;

                            while(true) {
                                //lecture nom du champ
                                string str4;
                                getline(is,str4);
                                if(str4=="Fchap") {
                                  if(process.nom_calcul=="latin") is.read( (char *)Inter[q].side[j].t[pt].Fchap.ptr(), sizeof(T)*Inter[q].side[j].t[pt].Fchap.size() );
                                  if(process.nom_calcul=="incr") is.read( (char *)Inter[q].side[j].t_post[pt].Fchap.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].Fchap.size() );
                                } else if(str4=="F") {
                                  if(process.nom_calcul=="latin") is.read( (char *)Inter[q].side[j].t[pt].F.ptr(), sizeof(T)*Inter[q].side[j].t[pt].F.size() );
                                  if(process.nom_calcul=="incr") is.read( (char *)Inter[q].side[j].t_post[pt].F.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].F.size() );
                                  //std::cout << Inter[q].side[j].t_post[pt].F << std::endl;
                                } else if(str4=="Wchap") {
                                  if(process.nom_calcul=="latin") is.read( (char *)Inter[q].side[j].t[pt].Wchap.ptr(), sizeof(T)*Inter[q].side[j].t[pt].Wchap.size() );
                                  if(process.nom_calcul=="incr") is.read( (char *)Inter[q].side[j].t_post[pt].Wchap.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].Wchap.size() );
                                } else if(str4=="Wpchap") {
                                  if(process.nom_calcul=="latin") is.read( (char *)Inter[q].side[j].t[pt].Wpchap.ptr(), sizeof(T)*Inter[q].side[j].t[pt].Wpchap.size() );
                                  if(process.nom_calcul=="incr") is.read( (char *)Inter[q].side[j].t_post[pt].Wpchap.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].Wpchap.size() );
                                } else if(str4=="W") {
                                  if(process.nom_calcul=="latin") is.read( (char *)Inter[q].side[j].t[pt].W.ptr(), sizeof(T)*Inter[q].side[j].t[pt].W.size() );
                                  if(process.nom_calcul=="incr") is.read( (char *)Inter[q].side[j].t_post[pt].W.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].W.size() );
                                } else if(str4=="Wp") {
                                  if(process.nom_calcul=="latin") is.read( (char *)Inter[q].side[j].t[pt].Wp.ptr(), sizeof(T)*Inter[q].side[j].t[pt].Wp.size() );
                                  if(process.nom_calcul=="incr") is.read( (char *)Inter[q].side[j].t_post[pt].Wp.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].Wp.size() );
                                } else if(str4=="WtildeM") {
                                  if(process.nom_calcul=="latin") is.read( (char *)Inter[q].side[j].t[pt].WtildeM.ptr(), sizeof(T)*Inter[q].side[j].t[pt].WtildeM.size() );
                                  if(process.nom_calcul=="incr") is.read( (char *)Inter[q].side[j].t_post[pt].WtildeM.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].WtildeM.size() );
                                } else if(str4=="finchamps") {
                                    break;
                                }  else {
                                    std::cout << "Erreur dans le choix des champs a lire " << std::endl;
                                    assert(0);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    is.close();
    }
}


/** \ingroup  PostTraitements
\brief Sauvegarde des quantites désirées pour les sous-structures
 
*/
template<class TV1>
void save_data_sst(TV1 &S, Param &process, Vec<string> &fields_to_save) {
    typedef TYPEREEL T;

    //nom du fichier de sauvegarde
    system(("mkdir -p "+process.affichage->repertoire_save).c_str());

    ifstream is;
    std::ostringstream ss1;
    ss1<<process.affichage->repertoire_save<<"save_sst_"<<process.rank;
    string name(ss1.str());
    ofstream os( name.c_str() );

    os << process.size << std::endl;
    if (process.latin->save_depl_SST ==1)
    for(unsigned q=0;q<S.size();q++) {
        os<< "S " << S[q].num <<std::endl;
        for(unsigned pt=1;pt<=process.temps->nbpastemps;pt++) {
            os << "Pastemps " << pt <<std::endl;
            for(unsigned k=0;k<fields_to_save.size();k++) {
                if(fields_to_save[k]=="q") {
                    os << "q" << std::endl;
                    if(process.nom_calcul=="latin") os.write( (char *)S[q].t[pt].q.ptr(), sizeof(T)*S[q].t[pt].q.size() );
                    if(process.nom_calcul=="incr")  os.write( (char *)S[q].t_post[pt].q.ptr(), sizeof(T)*S[q].t_post[pt].q.size() );
                    //std::cout << process.rank << " sst " << q << " pas de temps" << pt << " " << norm_2(S[q].t[pt].q) << " " << S[q].t[pt].q.size() << std::endl;
                } else {
                    std::cout << "Erreur dans le choix des champs a sauvegarder " << std::endl;
                    assert(0);
                }
            }
            os << "finchamps"<<std::endl;
        }
        os << "finPastemps" <<std::endl;
    }
    os << "fin"<<std::endl;
    os.close();
}

/** \ingroup  PostTraitements
\brief Lecture des quantites désirées pour les interfaces
 */
template<class TV1>
void read_data_sst(TV1 &S, Param &process) {
    typedef TYPEREEL T;
    //nom du fichier de sauvegarde
    ifstream is;
    std::ostringstream ss1;
    ss1<<process.affichage->repertoire_save<<"save_sst_"<<process.rank;
    string name(ss1.str());
    is.open( name.c_str() );


    unsigned toutlire=1;
    string str;
    if (is) {
      getline(is,str);
      if (process.size != 1 and process.size != atoi(str.c_str()) ) {//ca va pas aller
        std::cout << "Pour lire les donnees soit on le fait monoprocesseur soit on doit avoir le meme nombre de pro que le calcul initial" << std::endl;
        assert(0);
      } else if (process.size != atoi(str.c_str()))//le master doit lire tout les fichiers
        toutlire=atoi(str.c_str());
    } else {
      std::cout << "Aucun fichier de sauvegarde a lire" << std::endl;
      assert(0);
    }
    is.close();
    
    
    for( unsigned i=0;i<toutlire ;i++ ){
      std::ostringstream ss2;
      if (toutlire==1) ss2<<process.affichage->repertoire_save<<"save_sst_"<<process.rank;
      else ss2<<process.affichage->repertoire_save<<"save_sst_"<<i;
      string name2(ss2.str());
      is.open( name2.c_str() );
      getline(is,str);

      while(true) {
        //lecture numero de la sst
        string str;
        getline(is,str);

        if(str=="fin") {
          is.close();
          break;
        } else {
            istringstream s(str);
            string nom; //nom = Sst
            s >> nom;
            unsigned q; //numero de la sst
            s >> q ;

            while(true) {
                //lecture du pas de temps
                string str3;
                getline(is,str3);
                if(str3=="finPastemps") {
                    break;
                } else {
                    istringstream s3(str3);
                    string nom3; //nom = Pastemps
                    s3 >> nom3;
                    unsigned pt; //numero du pas de temps
                    s3 >> pt ;

                    while(true) {
                        //lecture nom du champ
                        string str4;
                        getline(is,str4);
                        //std::cout << str4 << std::endl;
                        if(str4=="q") {
                          if(process.nom_calcul=="latin") is.read( (char *)S[q].t[pt].q.ptr(), sizeof(T)*S[q].t[pt].q.size() );
                          if(process.nom_calcul=="incr") is.read( (char *)S[q].t_post[pt].q.ptr(), sizeof(T)*S[q].t_post[pt].q.size() );
                         } else if(str4=="finchamps") {
                            break;
                        } else {
                            std::cout << "Erreur dans le choix des champs a lire " << std::endl;
                            assert(0);
                        }
                    }
                }
            }
        }
    }
    }
}
