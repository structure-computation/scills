#include "mesh/remove_doubles.h"

using namespace LMT;
using namespace std;

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
        e.numsst_skin = parent.numsst;
        e.num_proc_skin = parent.num_proc;
        e.typmat_skin = parent.typmat;
      //e.epsilon_skin[0] =parent.Yde[0];
      //e.epsilon_skin[1] =parent.Yde[1];

//#warning A modifier pour les comportements orthotrope
      #ifdef FORMUORTHO
        e.sigma_local_skin = parent.sigma_local[0];
      #endif
   }
};

//procedure permettant de sortir les fichiers paraview pour le volume et la peau des sst en une seule passe.
template<class TV1> void write_paraview_results(TV1 &S,Param &process, DataUser &data_user) {
    
    //preparation des noms et des repertoires pour ecriture des resultats
    String name_multiresolution="";
    if(data_user.options.Multiresolution_on==1)
        name_multiresolution<<"resolution_"<<data_user.options.Multiresolution_current_resolution<<"_";
    Vec<string,2> directory_names=Vec<string>(process.affichage->repertoire_save +"results/sst_bulk",process.affichage->repertoire_save +"results/sst_skin"); 
    Vec<string,2> generic_names;
    for(int i=0;i<2;i++) {
        int tmp=system(("mkdir -p "+directory_names[i]).c_str()); //creation des repertoires
        generic_names[i]=directory_names[i]+"/"+name_multiresolution.c_str()+ process.affichage->name_data; //nom generique du fichier vtu
    }
   
   //eclatement de chaque sous-structure
   double ecl=1;
   eclat_SST(S,ecl);
   
    for(unsigned imic=1;imic<process.temps->nbpastemps+1;imic++){
        //creation du maillage global
        typename TV1::template SubType<0>::T::TMESH::TM meshglob;
        for(unsigned i=0;i<S.size();++i){
            if(process.nom_calcul=="incr")
                assign_dep_cont_slave(S[i],S[i].t_post[imic].q, data_user);
            else if(process.nom_calcul=="latin")
                assign_dep_cont_slave(S[i],S[i].t[imic].q, data_user);
            else{std::cout << "Type de calcul non reconnu dans affich_SST_resultat " << std::endl;assert(0);}
            meshglob.append(*S[i].mesh.m);
            S[i].mesh.unload();
        }
     
        //extraction des quantites pour la peau
        meshglob.sub_mesh(LMT::Number<1>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
        meshglob.sub_mesh(LMT::Number<2>()).elem_list.change_hash_size( meshglob, meshglob.elem_list.size() /2 +1);
    /*         remove_doubles(meshglob,1e-8, true);*/
        meshglob.update_skin();
        apply(meshglob.skin.elem_list,Projection_sigma_epsilon_on_skin(),meshglob.skin,meshglob);

        //ecriture des fichiers vtu
        for(unsigned i=0;i<2;i++){
            DisplayParaview dp;
            //ecriture fichier paraview generique de tous les champs (volume ou peau)
            ostringstream sp;
            sp<<"./tmp/paraview_"<<process.rank<<"_";
            string strp(sp.str());
            if(process.size > 1) process.affichage->save="save";
            if(i==0) dp.add_mesh(meshglob,strp.c_str(),process.affichage->display_fields_sst_bulk);
            else dp.add_mesh(meshglob.skin,strp.c_str(),process.affichage->display_fields_sst_skin);
            
            if(process.affichage->save=="display") dp.exec();
            //modification du nom et deplacement du fichier generique
            ostringstream ss;
            ss<<generic_names[i] << "_proc_"<<S[0].num_proc<<"_time_"<<imic<<".vtu";
            string namefile(ss.str());
            int tmp2=system(("mv "+strp+"0.vtu "+namefile).c_str());
        }
    } 
   
}



/**\ingroup Post_traitement 
 \brief Procedure permettant d'afficher les champs après calcul sur les sous-structures pour un calcul itératif (statique ou quasistatique)
 
 Pour un calcul quasistatique en latin, on a stocké dans le vecteur t la solution a convergence pour chaque pas de temps, c'est donc celle ci qui est affichée et pour laquelle on calcule contrainte déformation et autre quantité.
 */
template<class TV1> void affich_SST_resultat_latin(TV1 &S,Param &process, DataUser &data_user) {

   String name_multiresolution="";
   if(data_user.options.Multiresolution_on==1)
       name_multiresolution<<"resolution_"<<data_user.options.Multiresolution_current_resolution<<"_";
   
   string typemail=process.affichage->type_affichage;
   
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
    }
    else if(process.affichage->type_affichage== "Sbord"){
        process.affichage->display_fields.resize(8);
        process.affichage->display_fields[0]= "dep";
        process.affichage->display_fields[1]= "qtrans";
        process.affichage->display_fields[2]= "sigma_skin";
        process.affichage->display_fields[3]= "epsilon_skin";
        process.affichage->display_fields[5]= "sigma_mises_skin";
        process.affichage->display_fields[6]= "numsst_skin";
        process.affichage->display_fields[7]= "num_proc_skin";
    }
   
   
   string name_directory;   
   if(typemail=="Sinterieur") name_directory=process.affichage->repertoire_save +"results/sst_bulk";
   else name_directory=process.affichage->repertoire_save +"results/sst_skin";
   int tmp=system(("mkdir -p "+name_directory).c_str());
   string nom_generique = name_directory+"/"+name_multiresolution.c_str()+ process.affichage->name_data;
   string save=process.affichage->save;
     
   //int tmp=system(("mkdir -p "+process.affichage->repertoire_save).c_str());
   
   //eclatement de chaque sous-structure
   double ecl=1;
   eclat_SST(S,ecl);
   if (process.rank == 0 or process.rank == 1) std::cout << "Champs affiches pour les resultats : " << process.affichage->display_fields <<std::endl;
   
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
               assign_dep_cont_slave(S[i],S[i].t_post[imic].q, data_user);
            else if(process.nom_calcul=="latin")
               assign_dep_cont_slave(S[i],S[i].t[imic].q, data_user);
            else{std::cout << "Type de calcul non reconnu dans affich_SST_resultat " << std::endl;assert(0);}
            meshglob.append(*S[i].mesh.m);
            S[i].mesh.unload();
         }
             
         ostringstream ss;
         if (process.size == 1) ss<<nom_generique << "_time_"<<imic<<".vtu";
         else ss<<nom_generique << "_proc_"<<S[0].num_proc<<"_time_"<<imic<<".vtu";
         string namefile(ss.str());
         
         ostringstream sp;
         sp<<"./tmp/paraview_"<<process.rank<<"_";
         string strp(sp.str());
         //ecriture fichier paraview
         if(process.size > 1) process.affichage->save="save";
         dp.add_mesh(meshglob,strp.c_str(),process.affichage->display_fields);
         if(process.affichage->save=="display") dp.exec();
         int tmp2=system(("mv "+strp+"0.vtu "+namefile).c_str());
         
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
               assign_dep_cont_slave(S[i],S[i].t_post[imic].q, data_user); 
            else if(process.nom_calcul=="latin")
               assign_dep_cont_slave(S[i],S[i].t[imic].q, data_user);
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
         if (process.size == 1) ss<<nom_generique << "_time_"<<imic<<".vtu";
         else ss<<nom_generique << "_proc_"<<S[0].num_proc<<"_time_"<<imic<<".vtu";
         string namefile(ss.str());
         
         ostringstream sp;
         sp<<"./tmp/paraview_"<<process.rank<<"_";
         string strp(sp.str());
         //ecriture fichier paraview
         if(process.size > 1) process.affichage->save="save";
         dp.add_mesh(meshglob.skin,strp.c_str(),process.affichage->display_fields);
         if(process.affichage->save=="display") dp.exec();
         int tmp2=system(("mv "+strp+"0.vtu "+namefile).c_str());
         
      }
   }
   

/*   //Tri du vecteur display_fields de facon a le mettre dans le meme ordre que les donnees par noeud et par element pour que la creation des PVTU se fasse dans le bon ordre
   typedef typename TV1::template SubType<0>::T::TMESH TM;
   const char *names[TM::TNode::nb_params+(TM::TNode::nb_params==0)];
   DM::get_names<typename TM::TNode>( names );
   Vec<string> display_fields_temp,display_fields=process.affichage->display_fields;
   for(unsigned i=0;i<TM::TNode::nb_params;++i)
     if ( std::find(display_fields.begin(),display_fields.end(),std::string(names[i]))!=display_fields.end())
       display_fields_temp.push_back(std::string(names[i]));
   
     Data_vtk_extract_elem<true> dve;
     S[0].mesh.elem_list.apply_static(dve);
     apply( S[0].mesh.elem_list, dve, display_fields );
     for(typename Data_vtk_extract_elem<true>::Map::const_iterator iter=dve.mapd.begin();iter!=dve.mapd.end();++iter)
       if ( std::find(display_fields.begin(),display_fields.end(),iter->first)!=display_fields.end())
         display_fields_temp.push_back(iter->first);

     process.affichage->display_fields=display_fields_temp;*/
     
};



template<unsigned dim>
struct assign_effort{
   template<class TE, class TV> void operator()(TE &e, unsigned &compt, TV &F) const{
   Vec<unsigned> rep=range(compt*dim,(compt+1)*dim);
   e.F=F[rep];
   compt+=1;
   }
};

template<unsigned dim>
struct assign_deplacement{
   template<class TE, class TV> void operator()(TE &e, unsigned &compt, TV &W) const{
   Vec<unsigned> rep=range(compt*dim,(compt+1)*dim);
   e.W=W[rep];
   compt+=1;
   }
};

template<unsigned dim>
struct assign_saut_normal{
   template<class TE, class TV> void operator()(TE &e, unsigned &compt, TV &W) const{
   Vec<unsigned> rep=range(compt*dim,(compt+1)*dim);
   e.dWn=W[rep];
   compt+=1;
   }
};

template<unsigned dim>
struct assign_saut_tangent{
   template<class TE, class TV> void operator()(TE &e, unsigned &compt, TV &W) const{
   Vec<unsigned> rep=range(compt*dim,(compt+1)*dim);
   e.dWt=W[rep];
   compt+=1;
   }
};

template<unsigned dim>
struct assign_dissipation{
   template<class TE, class TV> void operator()(TE &e, unsigned &compt, TV &W) const{
   Vec<unsigned> rep=range(compt*dim,(compt+1)*dim);
   e.dissipation=W[rep];
   compt+=1;
   }
};

/** \ingroup  Post_Traitement 
\brief Assignation des efforts et déplacements d'interface

On affiche ici les valeurs de l'effort et du deplacement en chaque élément de l'interface pour le coté considéré (par défaut 0). 
Cette procédure s'applique pour une interface seulement.
*/
template<class INTER,class TV1> void assignation_INTER_F_W_latin(INTER &Inter,TV1 &S,unsigned data=0,unsigned pt=1,double dt=1.0){

   Vec<TYPEREEL, DIM> normale = Inter.G-S[Inter.vois[data*2]].G;
   double sign=1.;
   if (dot(Inter.side[data].neq[range(0,(int)DIM)],normale)<=0.00001)   sign=-1.;
   //assignation des deplacements et contraintes      
   unsigned numelem=0;
   apply(Inter.side[data].mesh->elem_list,assign_effort<DIM>(),numelem,sign*Inter.side[data].t[pt].Fchap);
   numelem=0;
   apply(Inter.side[data].mesh->elem_list,assign_deplacement<DIM>(),numelem,sign*Inter.side[data].t[pt].Wpchap);

   int type=0;
   if (Inter.type=="Ext" and Inter.comp=="depl"){type=0;}
   else if (Inter.type=="Ext" and Inter.comp=="effort"){type=1;}
   else if (Inter.type=="Ext" and ( Inter.comp=="sym" )){type=2;}
   else if (Inter.type=="Ext" and ( Inter.comp=="depl_normal")){type=3;}
   else if (Inter.type=="Int" and Inter.comp=="Parfait"){type=4;}
   else if (Inter.type=="Int" and (Inter.comp=="Contact" or Inter.comp=="Contact_jeu" or Inter.comp=="Contact_jeu_physique" or Inter.comp=="Contact_ep") ){type=5;}
   else if (Inter.type=="Int" and Inter.comp=="Jeu_impose"){type=6;}
   else if (Inter.type=="Ext" and Inter.comp=="periodique"){type=7;}
   else {type=8;}
   int numelem1=0;
   apply(Inter.side[data].mesh->elem_list,apply_type_elem_interface(),type,Inter.num,numelem1);
   numelem1=0;
   if ( Inter.comp=="Contact_jeu_physique" or Inter.comp=="periodique") apply(Inter.side[1-data].mesh->elem_list,apply_type_elem_interface(),type,Inter.num,numelem1);

   if (Inter.comp=="Contact" or Inter.comp=="Contact_jeu" or Inter.comp=="Contact_jeu_physique" or Inter.comp=="Contact_ep") {
   Vec<double> Un,Ut,tmp1,tmp0;
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
   Vec<double,DIM> a0;a0.set(0.);
   for( unsigned k=0;k< Un.size()/DIM; k++){
      if (norm_2(Un[range(k*DIM,(k+1)*DIM)]) > 1e-8 or norm_2(Ut[range(k*DIM,(k+1)*DIM)]) < 1e-8) {
        Ut[range(k*DIM,(k+1)*DIM)]=a0;
      }
   }
   numelem=0;
   apply(Inter.side[data].mesh->elem_list,assign_saut_normal<DIM>(),numelem,Un);
   numelem=0;
   apply(Inter.side[data].mesh->elem_list,assign_saut_tangent<DIM>(),numelem,Ut);
   
   
   //calcul de la dissipation par element puis assignation
   Vec<double> dissi_inter;
   dissi_inter.resize(Inter.side[0].t[0].Fchap.size());
   dissi_inter.set(0.);
/*    for(unsigned j=0 ;j<pt ;j++ ) {
        dissi_inter+=dt/2.0*(
                    (Inter.side[0].Pt(Inter.side[0].t[j+1].Fchap)+Inter.side[0].Pt(Inter.side[0].t[j].Fchap))*Inter.side[0].Pt(Inter.side[0].t[j+1].Wpchap)+
                    (Inter.side[1].Pt(Inter.side[1].t[j+1].Fchap)+Inter.side[1].Pt(Inter.side[1].t[j].Fchap))*Inter.side[1].Pt(Inter.side[1].t[j+1].Wpchap));
    }*/
   dissi_inter= Inter.side[0].Pt(Inter.side[0].t[pt].Fchap)*(Inter.side[0].M*(Inter.side[0].Pt(Inter.side[0].t[pt].Wpchap)-Inter.side[1].Pt(Inter.side[1].t[pt].Wpchap)));
   numelem=0;
    apply(Inter.side[data].mesh->elem_list,assign_dissipation<DIM>(),numelem,dissi_inter);

   }

}
template<class INTER,class TV1> void assignation_INTER_F_W_incr(INTER &Inter,TV1 &S,unsigned data=0,unsigned pt=1,double dt=1.0){

   Vec<TYPEREEL, DIM> normale = Inter.G-S[Inter.vois[data*2]].G;
   double sign=1.;
   if (dot(Inter.side[data].neq[range(0,(int)DIM)],normale)<=0.00001)   sign=-1.;
   //assignation des deplacements et contraintes  
   //PRINT("affichage des champs par interfaces");
   
   unsigned numelem=0;
   apply(Inter.side[data].mesh->elem_list,assign_effort<DIM>(),numelem,sign*Inter.side[data].t_post[pt].Fchap);
   numelem=0;
   apply(Inter.side[data].mesh->elem_list,assign_deplacement<DIM>(),numelem,sign*Inter.side[data].t_post[pt].Wpchap);

   int type=0;
   if (Inter.type=="Ext" and Inter.comp=="depl"){type=0;}
   else if (Inter.type=="Ext" and Inter.comp=="effort"){type=1;}
   else if (Inter.type=="Ext" and ( Inter.comp=="sym" )){type=2;}
   else if (Inter.type=="Ext" and ( Inter.comp=="depl_normal")){type=3;}
   else if (Inter.type=="Int" and Inter.comp=="Parfait"){type=4;}
   else if (Inter.type=="Int" and (Inter.comp=="Contact" or Inter.comp=="Contact_jeu" or Inter.comp=="Contact_jeu_physique" or Inter.comp=="Contact_ep") ){type=5;}
   else if (Inter.type=="Int" and Inter.comp=="Jeu_impose"){type=6;}
   else if (Inter.type=="Ext" and Inter.comp=="periodique"){type=7;}
   else {type=8;}
   int numelem1=0;
   apply(Inter.side[data].mesh->elem_list,apply_type_elem_interface(),type,Inter.num,numelem1);
   numelem1=0;
   if ( Inter.comp=="Contact_jeu_physique" or Inter.comp=="periodique") apply(Inter.side[1-data].mesh->elem_list,apply_type_elem_interface(),type,Inter.num,numelem1);

   if (Inter.comp=="Contact" or Inter.comp=="Contact_jeu" or Inter.comp=="Contact_jeu_physique" or Inter.comp=="Contact_ep") {
   Vec<double> Un,Ut,tmp1,tmp0;
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
   Vec<double,DIM> a0;a0.set(0.);
   for( unsigned k=0;k< Un.size()/DIM; k++){
      if (norm_2(Un[range(k*DIM,(k+1)*DIM)]) > 1e-8 or norm_2(Ut[range(k*DIM,(k+1)*DIM)]) < 1e-8) {
        Ut[range(k*DIM,(k+1)*DIM)]=a0;
      }
   }

   numelem=0;
   apply(Inter.side[data].mesh->elem_list,assign_saut_normal<DIM>(),numelem,Un);
   numelem=0;
   apply(Inter.side[data].mesh->elem_list,assign_saut_tangent<DIM>(),numelem,Ut);
   
   //calcul de la dissipation par element puis assignation
   Vec<double> dissi_inter;
   dissi_inter.resize(Inter.side[0].t_post[0].Fchap.size());
   dissi_inter.set(0.);
   for(unsigned j=0 ;j<pt ;j++ ) {
       dissi_inter+=dt/2.0*(
               (Inter.side[0].Pt(Inter.side[0].t_post[j+1].Fchap)+Inter.side[0].Pt(Inter.side[0].t_post[j].Fchap))*Inter.side[0].Pt(Inter.side[0].t_post[j+1].Wpchap)+
               (Inter.side[1].Pt(Inter.side[1].t_post[j+1].Fchap)+Inter.side[1].Pt(Inter.side[1].t_post[j].Fchap))*Inter.side[1].Pt(Inter.side[1].t_post[j+1].Wpchap));
   }
   numelem=0;
   apply(Inter.side[data].mesh->elem_list,assign_dissipation<DIM>(),numelem,dissi_inter);
   }
  
}

template<class TV2> int find_inter_in_subi(TV2 &Inter, int &num){
    for(unsigned q=0;q<Inter.size();q++)
        if (Inter[q].num == num)
            return q;

    return -1;
}

template<class TV2,class TV1> void affich_inter_data_time(TV2 &Inter, TV1 &S, Param &process){

      unsigned pt=process.affichage->pt;
      unsigned data = process.affichage->side;// choix du cote
      DisplayParaview dp;
      typename TV2::template SubType<0>::T::TMESH meshglob;
      if(find(process.affichage->num_inter_select,LMT::_1==-1)==1){
         for(unsigned q=0;q<Inter.size();q++){
            unsigned side;
            if(Inter[q].type=="Ext") side=0;
            else side=data;
            if (S[Inter[q].vois[data*2]].num_proc==process.rank){
            if(process.nom_calcul=="incr")
                assignation_INTER_F_W_incr(Inter[q],S,side,pt,process.temps->dt);
            else if(process.nom_calcul=="latin")
                assignation_INTER_F_W_latin(Inter[q],S,side,pt,process.temps->dt);
               
            meshglob.append(*Inter[q].side[side].mesh);
            if ( Inter[q].comp=="Contact_jeu_physique" or Inter[q].comp=="periodique") {
                if(process.nom_calcul=="incr")
                    assignation_INTER_F_W_incr(Inter[q],S,1-side,pt,process.temps->dt);
                else if(process.nom_calcul=="latin")
                    assignation_INTER_F_W_latin(Inter[q],S,1-side,pt,process.temps->dt); 
                meshglob.append(*Inter[q].side[1-side].mesh);
            }
            }
         }
      }
      else {
         for(unsigned q=0;q<process.affichage->num_inter_select.size();q++){
            int testinter = find_inter_in_subi(Inter, process.affichage->num_inter_select[q]);
            if (testinter != -1) {
            unsigned qs = testinter;
            std::cout << qs << " " << Inter[qs].num << " " << process.affichage->num_inter_select[q] << std::endl;
            unsigned side;
            if(Inter[qs].type=="Ext") side=0;
            else side=data;
            if (S[Inter[qs].vois[data*2]].num_proc==process.rank){
              if(process.nom_calcul=="incr")
                  assignation_INTER_F_W_incr(Inter[qs],S,side,pt,process.temps->dt);
            else if(process.nom_calcul=="latin")
                assignation_INTER_F_W_latin(Inter[qs],S,side,pt,process.temps->dt); 
            meshglob.append(*Inter[qs].side[side].mesh);
            if ( Inter[qs].comp=="Contact_jeu_physique" or Inter[qs].comp=="periodique") {
                if(process.nom_calcul=="incr")
                    assignation_INTER_F_W_incr(Inter[qs],S,1-side,pt,process.temps->dt);
                else if(process.nom_calcul=="latin")
                    assignation_INTER_F_W_latin(Inter[qs],S,1-side,pt,process.temps->dt); 
                meshglob.append(*Inter[qs].side[1-side].mesh);
            }
            }
            }
         }
      }

      ostringstream sp;
      sp<<"./tmp/paraview_"<<process.rank<<"_";
      string strp(sp.str());
      //dp.add_mesh(meshglob,strp.c_str(),process.affichage->display_fields);
      dp.add_mesh(meshglob,strp.c_str(),Vec<string>("num","type","qtrans","F","W"));
      if(process.affichage->save=="display") dp.exec();

}

/**\ingroup Post_traitement 
 \brief Procedure permettant d'afficher les champs après calcul sur les interfaces pour chaque piquet de temps
 */
template<class TV2,class TV1> void affich_INTER_resultat(TV2 &Inter,TV1 &S,Param &process) {
   
   string save_directory=process.affichage->repertoire_save+"results/inter/";
   
   int tmp=system(("mkdir -p "+save_directory).c_str());
   
   if(process.size > 1) process.affichage->save="save";

   for(unsigned i=1;i<process.temps->nbpastemps+1;i++){
      process.temps->pt=i;
      process.affichage->pt=i;
      affich_inter_data_time(Inter, S, process);
      string nom_generique = save_directory + process.affichage->name_data;
      ostringstream sp;
      sp<<"./tmp/paraview_"<<process.rank<<"_";
      string strp(sp.str());
      ostringstream ss;
      /*if (process.size == 1) ss<<nom_generique << "_time_"<<i<<".vtu";
      else */
      ss<<nom_generique << "_proc_"<<process.rank<<"_time_"<<i<<".vtu";
      string namefile(ss.str());
      int tmp2=system(("mv "+strp+"0.vtu "+namefile).c_str());
   }

};
