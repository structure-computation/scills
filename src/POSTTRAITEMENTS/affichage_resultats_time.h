#ifndef AFFICHAGE_RESULTATS_TIME_H
#define AFFICHAGE_RESULTATS_TIME_H

#include "Process.h"
#include "../DEFINITIONS/Sst.h"
#include "../DEFINITIONS/Interface.h"
#include "../DEFINITIONS/MultiResolutionData.h"
#include "../COMPUTE/DataUser.h"

#include "affichage_mesh_SST.h"
#include "displayparaview2.h"
#include "../ITERATIONS/manipulate_quantities.h"
#include "../UTILS/Sc2String.h"

#include "../../LMT/include/containers/vec.h"
#include "../../LMT/include/containers/vecpointedvalues.h"
#include "../../LMT/include/mesh/remove_doubles.h"

using namespace LMT;

struct Projection_sigma_epsilon_on_skin{
    template<class TE, class TMS, class TM> void operator()(TE &e, const TMS &mskin, const TM &m) const{
	const typename TM::EA *ea = mskin.get_parents_of(e)[0];
	typedef typename TM::TElemList::template SubType<0>::T TypeParent;
	const TypeParent &parent = static_cast<const TypeParent &>( *ea );
	e.sigma_skin = parent.sigma[0];
	e.epsilon_skin = parent.epsilon[0];
	e.numsst_skin = parent.numsst;
	e.num_proc_skin = parent.num_proc;
	e.typmat_skin = parent.typmat;
        e.id_group_skin = parent.id_group;
	
	//e.epsilon_skin[0] =parent.Yde[0];
	//e.epsilon_skin[1] =parent.Yde[1];

//#warning A modifier pour les comportements orthotrope
      #ifdef FORMUORTHO
	e.sigma_local_skin = parent.sigma_local[0];
      #endif
  }
};
#include "../MAILLAGE/mesh_data_accessors_sst.h"
///procedure permettant de sortir les fichiers paraview pour le volume et la peau des sst en une seule passe.
void write_paraview_results(PointedSubstructures &S,Process &process, DataUser &data_user) {
    std::cout << "Sauvegarde des resultats sur les sous-structures, pas de temps" << process.temps->pt_cur << std::endl;
    ///preparation des noms et des repertoires pour ecriture des resultats
    Sc2String name_multiresolution="";
    if(process.multiresolution->nb_calculs>1)
	name_multiresolution<<"_resolution_"<<process.multiresolution->m;
    Vec<Sc2String,2> directory_names=Vec<Sc2String>(process.affichage->repertoire_save +"results/sst_bulk",process.affichage->repertoire_save +"results/sst_skin"); 
    Vec<Sc2String,2> generic_names;
    for(int i=0;i<2;i++) {
	int tmp=system(("mkdir -p "+directory_names[i]).c_str()); //creation des repertoires
	generic_names[i] = directory_names[i] + "/" + process.affichage->name_data + name_multiresolution ; //nom generique du fichier vtu
    }
    
    ///eclatement de chaque sous-structure
    double ecl=1;
    eclat_SST(S,ecl);
    
    //for(unsigned imic=1;imic<process.temps->nbpastemps+1;imic++){     /// TMP, test sauvegarde a la fin de chaque pas de temps
//         ///creation du maillage global
//         SstMesh::TM meshglob;
//         for(unsigned i=0;i<S.size();++i){
//             //rebuild_state(S[i],S[i].t_post[imic], process);     /// TMP, test sauvegarde a la fin de chaque pas de temps
//             rebuild_state(S[i],S[i].t_post[process.temps->pt_cur], process);
//             meshglob.append(*S[i].mesh.m);
//             S[i].mesh.unload();
//         }
//         
//         ///extraction des quantites pour la peau
//         meshglob.sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
//         meshglob.sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
//         //remove_doubles(meshglob,1e-8, true);
//         meshglob.update_skin();
//         apply(meshglob.skin.elem_list,Projection_sigma_epsilon_on_skin(),meshglob.skin,meshglob);
//         
//         ///ecriture des fichiers vtu
//         for(unsigned i=0;i<2;i++){
//             DisplayParaview dp;
//             
//             ///nom du fichier paraview genere (volume ou peau)
//             ostringstream ss;
//             //ss<<generic_names[i] << "_proc_"<<S[0].num_proc<<"_time_"<<imic<<".vtu";      /// TMP, test sauvegarde a la fin de chaque pas de temps
//             ss<<generic_names[i] << "_proc_"<<S[0].num_proc<<"_time_"<<process.temps->pt_cur<<".vtu";   /// TMP, test sauvegarde a la fin de chaque pas de temps
//             Sc2String namefile(ss.str());
//             
//             ///ecriture fichier paraview generique de tous les champs (volume ou peau)
//             ostringstream sp;
//             sp<<"./tmp/paraview_"<<process.parallelisation->rank<<"_";
//             Sc2String strp(sp.str());
//             
//             if(process.parallelisation->is_multi_cpu()) process.affichage->save="save";
//             if(i==0) dp.add_mesh(meshglob,namefile.c_str(),process.affichage->display_fields_sst_bulk);
//             else dp.add_mesh(meshglob.skin,namefile.c_str(),process.affichage->display_fields_sst_skin);
//             
//             if(process.affichage->save=="display") dp.exec();
//             
//             ///modification du nom et deplacement du fichier generique
//             
//             //int tmp2=system(("mv "+strp+"0.vtu "+namefile).c_str());
//         }
//         
//     //}     /// TMP, test sauvegarde a la fin de chaque pas de temps
    
    
	///ecriture des maillages pour chaque sst dans un repertoire s�par�
	for(unsigned i=0;i<S.size();++i) {
	    SstMesh::TM meshglob;
	    rebuild_state(S[i],S[i].t_post[process.temps->pt_cur], process);
	    meshglob.append(*S[i].mesh.m);
	    S[i].mesh.unload();
	    ///extraction des quantites pour la peau
	    meshglob.sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
	    meshglob.sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
	    meshglob.update_skin();
	    apply(meshglob.skin.elem_list,Projection_sigma_epsilon_on_skin(),meshglob.skin,meshglob);
	    for(unsigned j=0;j<2;j++){
	      DisplayParaview dp;
	      
	      ///nom du fichier paraview genere (volume ou peau)
	      ostringstream ss;
	      //ss<<generic_names[i] << "_proc_"<<S[0].num_proc<<"_time_"<<imic<<".vtu";      /// TMP, test sauvegarde a la fin de chaque pas de temps
	      ss << generic_names[j] << "_sst_id_" << S[i].id << "_time_"<<process.temps->pt_cur;   /// TMP, test sauvegarde a la fin de chaque pas de temps
	      Sc2String namefile(ss.str());
	      
	      if(process.parallelisation->is_multi_cpu()) process.affichage->save="save";
	      if(j==0) dp.add_mesh(meshglob,namefile.c_str(),process.affichage->display_fields_sst_bulk);
	      else dp.add_mesh(meshglob.skin,namefile.c_str(),process.affichage->display_fields_sst_skin);
	      
	      if(process.affichage->save=="display") dp.exec();
	    }
	

	}
    
    
    
}



/**\ingroup Post_traitement 
\brief Procedure permettant d'afficher les champs apr�s calcul sur les sous-structures pour un calcul it�ratif (statique ou quasistatique)

Pour un calcul quasistatique en latin, on a stock� dans le vecteur t la solution a convergence pour chaque pas de temps, c'est donc celle ci qui est affich�e et pour laquelle on calcule contrainte d�formation et autre quantit�.
*/
template<class TV1> void affich_SST_resultat_latin(TV1 &S,Process &process, DataUser &data_user) {
    Sc2String name_multiresolution="";
    if(process.multiresolution->nb_calculs>1)
	name_multiresolution<<"resolution_"<<process.multiresolution->m<<"_";
  
    Sc2String typemail=process.affichage->type_affichage;
  
    if(process.affichage->type_affichage== "Sinterieur"){
	process.affichage->display_fields.resize(10);
	process.affichage->display_fields[0]= "dep";
	process.affichage->display_fields[1]= "qtrans";
	process.affichage->display_fields[2]= "sigma";
	process.affichage->display_fields[3]= "epsilon";
	process.affichage->display_fields[4]= "ener";
	process.affichage->display_fields[5]= "sigma_mises";
	process.affichage->display_fields[6]= "numsst";
	process.affichage->display_fields[7]= "f_vol_e";
	process.affichage->display_fields[8]= "num_proc";
        process.affichage->display_fields[9]= "id_group";
    }
    else if(process.affichage->type_affichage== "Sbord"){
	process.affichage->display_fields.resize(9);
	process.affichage->display_fields[0]= "dep";
	process.affichage->display_fields[1]= "qtrans";
	process.affichage->display_fields[2]= "sigma_skin";
	process.affichage->display_fields[3]= "epsilon_skin";
	process.affichage->display_fields[5]= "sigma_mises_skin";
	process.affichage->display_fields[6]= "numsst_skin";
	process.affichage->display_fields[7]= "num_proc_skin";
        process.affichage->display_fields[8]= "id_group_skin";
    }
  
  
  Sc2String name_directory;   
  if(typemail=="Sinterieur") name_directory=process.affichage->repertoire_save +"results/sst_bulk";
  else name_directory=process.affichage->repertoire_save +"results/sst_skin";
  int tmp=system(("mkdir -p "+name_directory).c_str());
  Sc2String nom_generique = name_directory+"/"+name_multiresolution.c_str()+ process.affichage->name_data;
  Sc2String save=process.affichage->save;
    
  //int tmp=system(("mkdir -p "+process.affichage->repertoire_save).c_str());
  
  //eclatement de chaque sous-structure
  double ecl=1;
  eclat_SST(S,ecl);
  if (process.parallelisation->is_local_cpu()) std::cout << "Champs affiches pour les resultats : " << process.affichage->display_fields <<std::endl;
  
/*   for(unsigned i=0;i<S.size();++i){
      //mise a zeros des quantites sigma_old et epsilon_old (visco uniquement)
      //apply(S[i].mesh.elem_list,mise_a_zero_old_quantities());
  }*/



  if (typemail=="Sinterieur") {
      for(unsigned imic=1;imic<process.temps->nbpastemps+1;imic++){
	DisplayParaview dp;
	typename TV1::template SubType<0>::T::TMESH::TM meshglob;
	for(unsigned i=0;i<S.size();++i){
	    if(process.nom_calcul=="incr")
		rebuild_state(S[i],S[i].t_post[imic].q, process);
	    else if(process.nom_calcul=="latin")
		rebuild_state(S[i],S[i].t[imic].q, process);
	    else{std::cout << "Type de calcul non reconnu dans affich_SST_resultat " << std::endl;assert(0);}
	    meshglob.append(*S[i].mesh.m);
	    S[i].mesh.unload();
	}
	    
	ostringstream ss;
	if (not process.parallelisation->is_multi_cpu()) ss<<nom_generique << "_time_"<<imic<<".vtu";
	else ss << nom_generique << "_" << S[0].num_proc << "_time_" << imic;
	Sc2String namefile(ss.str());
	
	
	//ecriture fichier paraview
	if(process.parallelisation->is_multi_cpu()) process.affichage->save="save";
	dp.add_mesh(meshglob,namefile.c_str(),process.affichage->display_fields);
	if(process.affichage->save=="display") dp.exec();

	
      }   
  }
  else if (typemail=="Sbord"){
      //sortie des valeurs de sigma et epsilon et deplacement sur la peau uniquement
      for(unsigned imic=1;imic<process.temps->nbpastemps+1;imic++){
	DisplayParaview dp;
	typename TV1::template SubType<0>::T::TMESH::TM meshglob;
	for(unsigned i=0;i<S.size();++i)
	{
	    if(process.nom_calcul=="incr")
		rebuild_state(S[i],S[i].t_post[imic].q, process); 
	    else if(process.nom_calcul=="latin")
		rebuild_state(S[i],S[i].t[imic].q, process);
	    else{std::cout << "Type de calcul non reconnu dans affich_SST_resultat " << std::endl;assert(0);}
//             S[i].mesh->update_skin();
	    meshglob.append(*S[i].mesh.m);
	    S[i].mesh.unload();
	}
	meshglob.sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
	meshglob.sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
	meshglob.update_skin();
	apply(meshglob.skin.elem_list,Projection_sigma_epsilon_on_skin(),meshglob.skin,meshglob);
	
	ostringstream ss;
	if (not process.parallelisation->is_multi_cpu()) ss<<nom_generique << "_time_"<<imic<<".vtu";
	else ss<<nom_generique << "_"<<S[0].num_proc<<"_time_"<<imic<<".vtu";
	Sc2String namefile(ss.str());
	
	
	//ecriture fichier paraview
	if(process.parallelisation->is_multi_cpu()) process.affichage->save="save";
	dp.add_mesh(meshglob.skin,namefile.c_str(),process.affichage->display_fields);
	if(process.affichage->save=="display") dp.exec();
	
      }
  }
  

/*   //Tri du vecteur display_fields de facon a le mettre dans le meme ordre que les donnees par noeud et par element pour que la creation des PVTU se fasse dans le bon ordre
  typedef typename TV1::template SubType<0>::T::TMESH TM;
  const char *names[TM::TNode::nb_params+(TM::TNode::nb_params==0)];
  DM::get_names<typename TM::TNode>( names );
  Vec<Sc2String> display_fields_temp,display_fields=process.affichage->display_fields;
  for(unsigned i=0;i<TM::TNode::nb_params;++i)
    if ( std::find(display_fields.begin(),display_fields.end(),Sc2String(names[i]))!=display_fields.end())
      display_fields_temp.push_back(Sc2String(names[i]));
  
    Data_vtk_extract_elem<true> dve;
    S[0].mesh.elem_list.apply_static(dve);
    apply( S[0].mesh.elem_list, dve, display_fields );
    for(typename Data_vtk_extract_elem<true>::Map::const_iterator iter=dve.mapd.begin();iter!=dve.mapd.end();++iter)
      if ( std::find(display_fields.begin(),display_fields.end(),iter->first)!=display_fields.end())
	display_fields_temp.push_back(iter->first);

    process.affichage->display_fields=display_fields_temp;*/
    
};


#include "../MAILLAGE/mesh_data_accessors_inter.h"


/** \ingroup  Post_Traitement 
* \brief Assignation des efforts et d�placements d'interface
* 
* Assignation des grandeurs de l'interface dans son maillage (pour un schema d'integration latin)
* Recupere les efforts (F), les deplacements (W) et l'endommagement (d) dans les resultats de l'interface 'Inter' (sur le bord indice 'data'),
* calcul les sauts en deplacement normal (dWn) et tangentiel (dWt) ainsi que la dissipation associee (dissipation) si utile.
* 'pt' et 'dt' sont l'indice du pas de temps (pour acceder a 't_post') et sa duree (pour le calcul des sauts et de la dissipation) 
*/
template<class INTER,class TV1> void assignation_INTER_F_W_latin(INTER &Inter,TV1 &S,unsigned data=0,unsigned pt=1,double dt=1.0){
    Point normale = Inter.G-S[Inter.vois[data*2]].G;
    Scalar sign = +1.0;
    if (dot(Inter.side[data].neq[range(0,(int)DIM)],normale)<=0.00001){
	sign = -1.0;
    }
    
    /// Assignation des deplacements et contraintes
    //std::cout << "\tAssignation de F : " << std::endl;
    Vector F = sign*Inter.side[data].t[pt].Fchap;
    upload_F(Inter.side[data],F);
    //std::cout << "\tAssignation de W : " << std::endl;
    Vector W = sign*Inter.side[data].t[pt].Wpchap;
    upload_W(Inter.side[data],W);
    
    /// Assignation du type d'element
    Vec<int,3> data_inter; //data interface to pass to visu : type, num, group_id
    data_inter[0] = Inter.get_type_elem();
    data_inter[1] = Inter.num; //
    data_inter[2] = 0; //group_id
    if (data_inter[0]==0) {
      data_inter[2]=Inter.edge_id;
    }
    else{
      data_inter[2]=Inter.id_link;
    }
    int numelem1=0;    
    apply(Inter.side[data].mesh->elem_list,apply_type_elem_interface(),data_inter,numelem1); 
    numelem1=0;
    if(Inter.comp == "Contact_jeu_physique" or Inter.comp == Interface::comp_periodique) {
	apply(Inter.side[1-data].mesh->elem_list,apply_type_elem_interface(),data_inter,numelem1);
    }

    /// Si les deux cotes de l'interface peuvent se decoller
    if(Inter.comp == "Contact" or
      Inter.comp == "Contact_jeu" or
      Inter.comp == "Contact_jeu_physique" or
      Inter.comp == Interface::comp_contact_parfait or
      Inter.comp == Interface::comp_contact_elastique or
      Inter.comp == Interface::comp_cassable_parfait or
      Inter.comp == Interface::comp_cassable_elastique or
      Inter.comp == Interface::comp_cohesive) {
	/// Calcul des sauts en deplacement normal et tangentiel
	Vector Un,Ut,tmp1,tmp0;
	Un.resize(Inter.side[0].nodeeq.size());Un.set(0.);
	Ut.resize(Inter.side[0].nodeeq.size());Ut.set(0.);
	tmp0.resize(Inter.side[0].nodeeq.size());tmp0.set(0.);
	tmp1.resize(Inter.side[0].nodeeq.size());tmp1.set(0.);
	tmp0=Inter.side[0].Pn(Inter.side[0].t[pt].Wchap);
	tmp1=Inter.side[1].Pn(Inter.side[1].t[pt].Wchap);
	Un=tmp0[Inter.side[0].ddlcorresp]-tmp1[Inter.side[1].ddlcorresp];
	tmp0=Inter.side[0].Pt(Inter.side[0].t[pt].Wchap);
	tmp1=Inter.side[1].Pt(Inter.side[1].t[pt].Wchap);
	Ut=tmp0[Inter.side[0].ddlcorresp]-tmp1[Inter.side[1].ddlcorresp];
	Point a0;a0.set(0.);
	for( unsigned k=0;k< Un.size()/DIM; k++){
	    if (norm_2(Un[range(k*DIM,(k+1)*DIM)]) > 1e-8 or norm_2(Ut[range(k*DIM,(k+1)*DIM)]) < 1e-8) {
		Ut[range(k*DIM,(k+1)*DIM)]=a0;
	    }
	}
	//std::cout << "\tAssignation de dWn : " << std::endl;
	upload_dWn(Inter.side[data],Un);
	//std::cout << "\tAssignation de dWt : " << std::endl;
	upload_dWt(Inter.side[data],Ut);
	
	/// Calcul de la dissipation par element puis assignation
	Vector dissi_inter;
	dissi_inter.resize(Inter.side[0].t[0].Fchap.size());
	dissi_inter.set(0.);
	/*for(unsigned j=0 ;j<pt ;j++ ) {
	    dissi_inter+=dt/2.0*(
			(Inter.side[0].Pt(Inter.side[0].t[j+1].Fchap)+Inter.side[0].Pt(Inter.side[0].t[j].Fchap))*Inter.side[0].Pt(Inter.side[0].t[j+1].Wpchap)+
			(Inter.side[1].Pt(Inter.side[1].t[j+1].Fchap)+Inter.side[1].Pt(Inter.side[1].t[j].Fchap))*Inter.side[1].Pt(Inter.side[1].t[j+1].Wpchap));
	}*/
	dissi_inter= Inter.side[0].Pt(Inter.side[0].t[pt].Fchap)*(Inter.side[0].M*(Inter.side[0].Pt(Inter.side[0].t[pt].Wpchap)-Inter.side[1].Pt(Inter.side[1].t[pt].Wpchap)));
	//std::cout << "\tAssignation de dissipation : " << std::endl;
	upload_dissipation(Inter.side[data],dissi_inter);
  }
}


/** \ingroup  Post_Traitement 
* \brief Assignation des efforts et d�placements d'interface
* 
* Assignation des grandeurs de l'interface dans son maillage (pour un schema d'integration incremental)
* Recupere les efforts (F), les deplacements (W) et l'endommagement (d) dans les resultats de l'interface 'Inter' (sur le bord indice 'data'),
* calcul les sauts en deplacement normal (dWn) et tangentiel (dWt) ainsi que la dissipation associee (dissipation) si utile.
* 'pt' et 'dt' sont l'indice du pas de temps (pour acceder a 't_post') et sa duree (pour le calcul des sauts et de la dissipation) 
*/
template<class INTER,class TV1> void assignation_INTER_F_W_incr(INTER &Inter,TV1 &S,unsigned data=0,unsigned pt=1,double dt=1.0){
    Point normale = Inter.G-S[Inter.vois[data*2]].G;
    Scalar sign = +1.0;
    if (dot(Inter.side[data].neq[range(0,(int)DIM)],normale)<=0.00001){
	sign = -1.0;
    }
    
    /// Assignation des deplacements et contraintes  
    //std::cout << "\tAssignation de F : " << std::endl;
    Vector F;
    const unsigned i_max = Inter.side[data].t_post[pt].Fchap.size();
    F.resize(i_max);
    for(unsigned i = 0; i < i_max; i++){
	F[i] = sign*Inter.side[data].t_post[pt].Fchap[i];
    }
    upload_F(Inter.side[data],F);
    //std::cout << "\tAssignation de W : " << std::endl;
    Vector W = sign*Inter.side[data].t_post[pt].Wpchap;
    upload_W(Inter.side[data],W);
    
    /// Assignation du type d'element
    Vec<int,3> data_inter; //data interface to pass to visu : type, num, group_id
    data_inter[0] = Inter.get_type_elem();
    data_inter[1] = Inter.num; //
    data_inter[2] = 0; //group_id
    if (data_inter[0]==0) {
      data_inter[2]=Inter.edge_id;
    }
    else{
      data_inter[2]=Inter.id_link;
    }
    int numelem1=0;
    apply(Inter.side[data].mesh->elem_list,apply_type_elem_interface(),data_inter,numelem1); 
    numelem1=0;
    if(Inter.comp == "Contact_jeu_physique" or Inter.comp == Interface::comp_periodique) {
	apply(Inter.side[1-data].mesh->elem_list,apply_type_elem_interface(),data_inter,numelem1);
    }
    
//     apply(Inter.side[data].mesh->elem_list,apply_type_elem_interface(),type,Inter.num,numelem1);
//     numelem1=0;
//     if(Inter.comp == "Contact_jeu_physique" or Inter.comp == Interface::comp_periodique) {
//         apply(Inter.side[1-data].mesh->elem_list,apply_type_elem_interface(),type,Inter.num,numelem1);
//     }
    
    /// Si les deux cotes de l'interface peuvent se decoller
    if(Inter.comp == "Contact" or
      Inter.comp == "Contact_jeu" or
      Inter.comp == "Contact_jeu_physique" or
      Inter.comp == Interface::comp_contact_parfait or
      Inter.comp == Interface::comp_contact_elastique or
      Inter.comp == Interface::comp_cassable_parfait or
      Inter.comp == Interface::comp_cassable_elastique or
      Inter.comp == Interface::comp_cohesive) {
      
	
	/// Calcul des sauts en deplacement normal et tangentiel
	Vector Un,Ut,tmp1,tmp0;
	Un.resize(Inter.side[0].nodeeq.size());Un.set(0.);
	Ut.resize(Inter.side[0].nodeeq.size());Ut.set(0.);
	tmp0.resize(Inter.side[0].nodeeq.size());tmp0.set(0.);
	tmp1.resize(Inter.side[0].nodeeq.size());tmp1.set(0.);
	
	tmp0=Inter.side[0].Pn(Inter.side[0].t_post[pt].Wchap);
	tmp1=Inter.side[1].Pn(Inter.side[1].t_post[pt].Wchap);
	Un=tmp0[Inter.side[0].ddlcorresp]-tmp1[Inter.side[1].ddlcorresp];
	tmp0=Inter.side[0].Pt(Inter.side[0].t_post[pt].Wchap);
	tmp1=Inter.side[1].Pt(Inter.side[1].t_post[pt].Wchap);
	Ut=tmp0[Inter.side[0].ddlcorresp]-tmp1[Inter.side[1].ddlcorresp];
	
	Point a0;a0.set(0.);
	
	for( unsigned k=0;k< Un.size()/DIM; k++){
	    if (norm_2(Un[range(k*DIM,(k+1)*DIM)]) > 1e-8 or norm_2(Ut[range(k*DIM,(k+1)*DIM)]) < 1e-8) {
		Ut[range(k*DIM,(k+1)*DIM)]=a0;
	    }
	}
	PRINT(Inter.id);
	
	//std::cout << "\tAssignation de dWn : " << std::endl;
	upload_dWn(Inter.side[data],Un);
	//std::cout << "\tAssignation de dWt : " << std::endl;
	upload_dWt(Inter.side[data],Ut);
	
	/// Calcul de la dissipation par element puis assignation
	/*  TMP, test sauvegarde a la fin de chaque pas de temps
	Vector dissi_inter;
	dissi_inter.resize(Inter.side[0].t_post[0].Fchap.size());
	dissi_inter.set(0.);
	for(unsigned j=0 ;j<pt ;j++ ) {
	    dissi_inter+=dt/2.0*(
		(Inter.side[0].Pt(Inter.side[0].t_post[j+1].Fchap)+Inter.side[0].Pt(Inter.side[0].t_post[j].Fchap))*Inter.side[0].Pt(Inter.side[0].t_post[j+1].Wpchap)+
		(Inter.side[1].Pt(Inter.side[1].t_post[j+1].Fchap)+Inter.side[1].Pt(Inter.side[1].t_post[j].Fchap))*Inter.side[1].Pt(Inter.side[1].t_post[j+1].Wpchap));
	}
	//std::cout << "\tAssignation de dissipation : " << std::endl;
	upload_dissipation(Inter.side[data],dissi_inter);
	//*/
    }
}


template<class TV2> int find_inter_in_subi(TV2 &Inter, int &num){
    for(unsigned q=0;q<Inter.size();q++)
	if (Inter[q].num == num)
	    return q;

    return -1;
}


template<class TV2,class TV1> void affich_inter_data_time(TV2 &Inter, TV1 &S, Process &process){
//     unsigned pt=process.affichage->pt;
//     unsigned data = process.affichage->side;    /// Choix du cote a sauvegarder
//     DisplayParaview dp;
//     InterfaceMesh meshglob;
//     if(find(process.affichage->num_inter_select,LMT::_1==-1)==1){
// 	for(unsigned q=0;q<Inter.size();q++){
// 	    unsigned side = (Inter[q].type == Interface::type_ext)? 0 : data;
// 	    if (S[Inter[q].vois[data*2]].num_proc==process.parallelisation->rank){
// 		/// Chargement des resultats sur le cote 'side' de l'interface
// 		assignation_INTER_F_W_incr(Inter[q],S,side,pt,process.temps->dt);
// 		meshglob.append(*Inter[q].side[side].mesh);
// 		/// Si necessaire, inclusion du deuxieme cote de l'interface
// 		if ( Inter[q].comp=="Contact_jeu_physique" or Inter[q].comp == Interface::comp_periodique) {
// 		    assignation_INTER_F_W_incr(Inter[q],S,1-side,pt,process.temps->dt);
// 		    meshglob.append(*Inter[q].side[1-side].mesh);
// 		}
// 	    }
// 	}
//     } else {
// 	for(unsigned q=0;q<process.affichage->num_inter_select.size();q++){
// 	    int testinter = find_inter_in_subi(Inter, process.affichage->num_inter_select[q]);
// 	    if (testinter != -1) {
// 		unsigned qs = testinter;
// 		unsigned side = (Inter[qs].type == Interface::type_ext)? 0 : data;
// 		if (S[Inter[qs].vois[data*2]].num_proc==process.parallelisation->rank){
// 		    assignation_INTER_F_W_incr(Inter[qs],S,side,pt,process.temps->dt);
// 		    meshglob.append(*Inter[qs].side[side].mesh);
// 		    if ( Inter[qs].comp=="Contact_jeu_physique" or Inter[qs].comp == Interface::comp_periodique) {
// 			assignation_INTER_F_W_incr(Inter[qs],S,1-side,pt,process.temps->dt);
// 			meshglob.append(*Inter[qs].side[1-side].mesh);
// 		    }
// 		}
// 	    }
// 	}
//     }
// 
//     ostringstream sp;
//     sp<<"./tmp/paraview_"<<process.parallelisation->rank<<"_";
//     Sc2String strp(sp.str());
//     
//     ///preparation des noms et des repertoires pour ecriture des resultats
//     Sc2String save_directory=process.affichage->repertoire_save+"results/inter/";
//     Sc2String base_filename= save_directory;
//     if(process.multiresolution->nb_calculs>1)
// 	base_filename<<"resolution_"<<process.multiresolution->m<<"_";
//     base_filename << "proc_" << process.parallelisation->rank << "_time_";
//     Sc2String namefile = base_filename;
//     namefile << process.temps->pt_cur << ".vtu";
//     
//     dp.add_mesh(meshglob,namefile.c_str(),process.affichage->display_fields_inter);
//     //dp.add_mesh(meshglob,strp.c_str(),Vec<Sc2String>("num","type","qtrans","F","W","dWn","dWt","d","dissipation"));
//     if(process.affichage->save=="display")
// 	dp.exec();
    
    unsigned pt=process.affichage->pt;
    unsigned data = process.affichage->side;    /// Choix du cote a sauvegarder

 
    
    if(find(process.affichage->num_inter_select,LMT::_1==-1)==1){
	for(unsigned q=0;q<Inter.size();q++){
	    unsigned side = (Inter[q].type == Interface::type_ext)? 0 : data;
	    if (S[Inter[q].vois[data*2]].num_proc==process.parallelisation->rank){
	        InterfaceMesh meshglob;
		/// Chargement des resultats sur le cote 'side' de l'interface
		assignation_INTER_F_W_incr(Inter[q],S,side,pt,process.temps->dt);
		meshglob.append(*Inter[q].side[side].mesh);
		/// Si necessaire, inclusion du deuxieme cote de l'interface
		if ( Inter[q].comp=="Contact_jeu_physique" or Inter[q].comp == Interface::comp_periodique) {
		    assignation_INTER_F_W_incr(Inter[q],S,1-side,pt,process.temps->dt);
		    meshglob.append(*Inter[q].side[1-side].mesh);
		}
		DisplayParaview dp;
		ostringstream sp;
		sp<<"./tmp/paraview_"<<process.parallelisation->rank<<"_";
		Sc2String strp(sp.str());
		
		///preparation des noms et des repertoires pour ecriture des resultats
		Sc2String save_directory=process.affichage->repertoire_save+"results/inter/";
		Sc2String base_filename= save_directory;
		base_filename << "result";
		if(process.multiresolution->nb_calculs>1)
		    base_filename <<"_resolution_" << process.multiresolution->m ;
		base_filename << "_inter_id_" << Inter[q].id <<"_time_";
		Sc2String namefile = base_filename;
		namefile << process.temps->pt_cur;	
		dp.add_mesh(meshglob,namefile.c_str(),process.affichage->display_fields_inter);
		//int tmp2=system(("mv "+strp+"0.vtu "+namefile).c_str()); 
	    }
	}
    } else {
	for(unsigned q=0;q<process.affichage->num_inter_select.size();q++){
	    int testinter = find_inter_in_subi(Inter, process.affichage->num_inter_select[q]);
	    if (testinter != -1) {
		unsigned qs = testinter;
		unsigned side = (Inter[qs].type == Interface::type_ext)? 0 : data;
		if (S[Inter[qs].vois[data*2]].num_proc==process.parallelisation->rank){
		    InterfaceMesh meshglob;
		    assignation_INTER_F_W_incr(Inter[qs],S,side,pt,process.temps->dt);
		    meshglob.append(*Inter[qs].side[side].mesh);
		    if ( Inter[qs].comp=="Contact_jeu_physique" or Inter[qs].comp == Interface::comp_periodique) {
			assignation_INTER_F_W_incr(Inter[qs],S,1-side,pt,process.temps->dt);
			meshglob.append(*Inter[qs].side[1-side].mesh);
		    }
		    DisplayParaview dp;
		    ostringstream sp;
		    sp<<"./tmp/paraview_"<<process.parallelisation->rank<<"_";
		    Sc2String strp(sp.str());
		    
		    ///preparation des noms et des repertoires pour ecriture des resultats
		    Sc2String save_directory=process.affichage->repertoire_save+"results/inter/";
		    Sc2String base_filename= save_directory;
		    base_filename << "result";
		    if(process.multiresolution->nb_calculs>1)
			base_filename<<"_resolution_"<<process.multiresolution->m;
		    base_filename << "_inter_id_" << Inter[q].id << "_time_";
		    Sc2String namefile = base_filename;
		    namefile << process.temps->pt_cur;
		    dp.add_mesh(meshglob,namefile.c_str(),process.affichage->display_fields_inter);
		    //int tmp2=system(("mv "+strp+"0.vtu "+namefile).c_str());
		}
	    }
	}
    }

 
}

/**\ingroup Post_traitement 
\brief Procedure permettant d'afficher les champs apr�s calcul sur les interfaces pour chaque piquet de temps
*/
template<class TV2,class TV1> void affich_INTER_resultat(TV2 &Inter,TV1 &S,Process &process) {
    std::cout << "Sauvegarde des resultats sur les interfaces, pas de temps" << process.temps->pt_cur << std::endl;
    Sc2String save_directory=process.affichage->repertoire_save+"results/inter/";

    int tmp=system(("mkdir -p "+save_directory).c_str());

    ///preparation des noms et des repertoires pour ecriture des resultats
    Sc2String base_filename= save_directory;
    if(process.multiresolution->nb_calculs>1)
	base_filename<<"resolution_"<<process.multiresolution->m<<"_";
    base_filename << "proc_" << process.parallelisation->rank << "_time_";
    
    if(process.parallelisation->is_multi_cpu()) process.affichage->save="save";

    //for(unsigned i=1;i<process.temps->nbpastemps+1;i++){  /// TMP, test sauvegarde a la fin de chaque pas de temps
	//process.temps->pt=i;      /// TMP, test sauvegarde a la fin de chaque pas de temps
	//process.affichage->pt=i;      /// TMP, test sauvegarde a la fin de chaque pas de temps
	process.affichage->pt=process.temps->pt_cur;        /// TMP, test sauvegarde a la fin de chaque pas de temps
	affich_inter_data_time(Inter, S, process);
//         ostringstream sp;
//         sp<<"./tmp/paraview_"<<process.parallelisation->rank<<"_";
//         Sc2String strp(sp.str());
//         Sc2String namefile = base_filename;
//         //namefile << i << ".vtu";    /// TMP, test sauvegarde a la fin de chaque pas de temps
//         namefile << process.temps->pt_cur << ".vtu";
//         int tmp2=system(("mv "+strp+"0.vtu "+namefile).c_str());
    //}     /// TMP, test sauvegarde a la fin de chaque pas de temps

};

#endif //AFFICHAGE_RESULTATS_TIME_H

