#include "DataUser.h"

using namespace Metil;
using namespace LMT;
using namespace std;


inline void read_data_process(Param &process, DataUser &data_user) {
    process.sousint = false;
    process.type_sousint = "h";
    process.rbm.bloq = false;
    process.rbm.mvts_bloques.resize(3);
    process.rbm.mvts_bloques[0]= "Ty";
    process.rbm.mvts_bloques[1]= "Tx";
    process.rbm.mvts_bloques[2]= "Rz";
    process.nb_threads = 1;
    process.save_data = true;
    process.read_data = false;
    process.reprise_calcul = 0;
    process.properties->deltaT = 0;
    
    process.multiscale->multiechelle = data_user.options.multiechelle;
    process.multiscale->type_base_macro = 3;
    process.multiscale->opti_multi = 0;
    process.multiscale->erreur_macro = 0.000001;
    
    process.latin->nbitermax = data_user.options.LATIN_nb_iter_max;
    process.latin->facteur_relaxation = 0.8;
    process.latin->critere_erreur = data_user.options.LATIN_crit_error;
    process.latin->critere_erreur_diss = 0;
    process.latin->type_error = "ddr";
    if (process.latin->type_error=="dissipation") process.latin->critere_erreur_diss = 0;
    
    process.latin->ktype = "scalaire_auto_CL";
    process.latin->kfact = 1;
    process.latin->copydirection = 0; 
    process.latin->list_error= 1;
    
    process.affichage->interactivite= 0;
    process.affichage->affich_resultat= 1;
    process.affichage->type_affichage= "Sinterieur";
    process.affichage->display_error= 0;
    process.affichage->affich_mesh= 1;
    process.affichage->save= "save";
    if(process.affichage->type_affichage== "Sinterieur"){
        process.affichage->display_fields.resize(6);
        process.affichage->display_fields[0]= "dep";
        process.affichage->display_fields[1]= "qtrans";
        process.affichage->display_fields[2]= "sigma";
        process.affichage->display_fields[3]= "epsilon";
        process.affichage->display_fields[4]= "ener";
        process.affichage->display_fields[5]= "sigma_mises";
    }
    else if(process.affichage->type_affichage== "Sbord"){
        process.affichage->display_fields.resize(6);
        process.affichage->display_fields[0]= "dep";
        process.affichage->display_fields[1]= "qtrans";
        process.affichage->display_fields[2]= "sigma_skin";
        process.affichage->display_fields[3]= "epsilon_skin";
        process.affichage->display_fields[5]= "sigma_mises_skin";
    }
    
    process.affichage->repertoire_save= data_user.calcul_path + "/";
    process.affichage->name_data= "result";
    process.affichage->command_file= "No";
    
//     std::cout<< "data_user.options.Temp_statique = " << data_user.options.Temp_statique << std::endl;
    
    if(data_user.options.Temp_statique == "statique"){
        process.temps->type_de_calcul= "stat";
    }else if(data_user.options.Temp_statique == "quasistatique"){
        process.temps->type_de_calcul= "Qstat";
    }
    
    if (process.temps->type_de_calcul=="stat") {
        if (process.rank==0) std::cout << "************************" << std::endl;
        if (process.rank==0) std::cout << "     STATIQUE     " << std::endl;
        if (process.rank==0) std::cout << "************************" << std::endl;
        process.temps->nbpastemps=1;
        process.temps->dt=1;
        process.nom_calcul="incr";
        process.temps->nb_step = 1;
        process.temps->time_step.resize(process.temps->nb_step);
        process.temps->nbpastemps=process.temps->nb_step;
        for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
            process.temps->time_step[i_step].dt=1;
            process.temps->time_step[i_step].t_ini=0;
            process.temps->time_step[i_step].t_fin=0;
            process.temps->time_step[i_step].nb_time_step=1;
        }
    
    }else if(process.temps->type_de_calcul=="Qstat") {
        if (process.rank==0) std::cout << "************************" << std::endl;
        if (process.rank==0) std::cout << "     QUASISTATIQUE      " << std::endl;
        if (process.rank==0) std::cout << "************************" << std::endl;
        process.nom_calcul="incr";
        process.temps->nb_step = data_user.time_step.size();
        process.temps->time_step.resize(process.temps->nb_step);
        process.temps->nbpastemps= 0;
        for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
            process.temps->time_step[i_step].dt = data_user.time_step[i_step].dt;
            process.temps->time_step[i_step].t_ini = data_user.time_step[i_step].ti;
            process.temps->time_step[i_step].t_fin = data_user.time_step[i_step].tf;
            process.temps->time_step[i_step].nb_time_step = data_user.time_step[i_step].nb_time_step;
            process.temps->nbpastemps += process.temps->time_step[i_step].nb_time_step;
        }
    }
};

/*// lecture des parametres de calcul
inline void read_data_process(Param &process, const XmlNode &n, DataUser &data_user) {
   process.sousint=0;
   process.type_sousint="h";
   process.rbm.bloq=0;
   process.rbm.mvts_bloques="";
   process.nb_threads=1;
   process.save_data=1;
   process.read_data=0;
   process.reprise_calcul=0;
   process.properties->deltaT=0;
   process.multiscale->type_base_macro=0;
   process.multiscale->multiechelle=data_user.multiechelle;
   process.multiscale->opti_multi=0;
   process.multiscale->erreur_macro=0;
   
   
//    nm.get_attribute("nb_macro_identique",process.multiscale->nbmacro_identique);   
//    if(process.multiscale->nbmacro_identique==0){
//       XmlNode nmacro = n.get_element("fcts_macro");
//       unsigned nbfcts= nmacro.nb_elements("inter");
//       
//       for(unsigned i=0;i<nbfcts;i++){
//          Vec<unsigned> num_inter;
//          XmlNode ninter = nmacro.get_element("inter",i);
//          unsigned nbmacro;
//          ninter.get_attribute("nbmacro",nbmacro);
//          ninter.get_attribute("num",num_inter);
//          process.multiscale->inter_correspondantes.append(num_inter);
//          Vec<unsigned> vecnbmacro;
//          vecnbmacro.resize(num_inter.size());
//          vecnbmacro.set(nbmacro);
//          process.multiscale->nb_fcts_macro_par_inter.append(vecnbmacro);         
//       }
//    }
//    
   
   process.latin->nbitermax=data_user.options.LATIN_nb_iter_max;
   process.latin->facteur_relaxation=0.8;
   process.latin->critere_erreur=data_user.options.LATIN_crit_error;
   process.latin->type_error="ddr";
   if (process.latin->type_error=="dissipation") process.latin->critere_erreur_diss = 0;
   process.latin->ktype
   process.latin->kfact
   process.latin->copydirection
   
   process.latin->list_error
   process.affichage->interactivite
   process.affichage->affich_resultat
   process.affichage->type_affichage
   process.affichage->display_error
   process.affichage->affich_mesh
   process.affichage->save
   process.affichage->display_fields
   process.affichage->repertoire_save
   process.affichage->name_data
   process.affichage->command_file
   
   process.temps->type_de_calcul
   
   
   nm.get_attribute("nbitermax",process.latin->nbitermax);
   nm.get_attribute("facteur_relaxation",process.latin->facteur_relaxation);
   nm.get_attribute("critere_erreur",process.latin->critere_erreur);
   nm.get_attribute("critere_erreur_diss",process.latin->critere_erreur_diss);
   nm.get_attribute("type_erreur",process.latin->type_error);
   if (process.latin->type_error=="dissipation") process.latin->critere_erreur_diss = 0;
   XmlNode nt =n.get_element("direction_recherche");
   nt.get_attribute("ktype",process.latin->ktype);
   nt.get_attribute("kfact",process.latin->kfact);
   nt.get_attribute("copydirection",process.latin->copydirection);

   XmlNode na =n.get_element("parametres_affichage");
   na.get_attribute("list_error",process.latin->list_error);   
   na.get_attribute("interactivite",process.affichage->interactivite);
   na.get_attribute("affich_resultat",process.affichage->affich_resultat);
   na.get_attribute("type_affichage",process.affichage->type_affichage);
   na.get_attribute("display_error",process.affichage->display_error);
   na.get_attribute("affich_mesh",process.affichage->affich_mesh);
   //na.get_attribute("affich_depl_pt",process.affichage->affich_depl_pt);
   //na.get_attribute("coor_point",process.affichage->coor_point);
   na.get_attribute("save_or_display",process.affichage->save);
   na.get_attribute("display_fields",process.affichage->display_fields);
   na.get_attribute("repertoire_save",process.affichage->repertoire_save);
   na.get_attribute("name_data",process.affichage->name_data);
   na.get_attribute("command_file",process.affichage->command_file);

   XmlNode ntp =n.get_element("parametres_temporels");

   XmlNode ntp_prop =ntp.get_element("proprietes");
   ntp_prop.get_attribute("type_de_calcul",process.temps->type_de_calcul);
   
   if (process.rank==0) cout << "************************" << endl;
   if (process.temps->type_de_calcul=="stat") {
      if (process.rank==0) cout << "     STATIQUE     " << endl;
      if (process.rank==0) cout << "************************" << endl;
      if (process.rank==0) cout << " Rq : 1 seul pas de temps automatiquement, dt=1 par defaut " << endl;
      process.temps->dt=1;
      process.nom_calcul="incr";
      ntp_prop.get_attribute("nb_step",process.temps->nb_step);
      process.temps->time_step.resize(process.temps->nb_step);
      process.temps->nbpastemps=process.temps->nb_step;
      for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
          process.temps->time_step[i_step].dt=1;
          process.temps->time_step[i_step].t_ini=0;
          process.temps->time_step[i_step].t_fin=0;
          process.temps->time_step[i_step].nb_time_step=1;
      }
      if (process.rank==0) cout << " Rq : Attention la valeur de la fonction spatiale sera tout de meme modulee par la fonction temporelle" << endl;
   } else if(process.temps->type_de_calcul=="Qstat") {
      if (process.rank==0) cout << "     QUASISTATIQUE          " << endl;
      if (process.rank==0) cout << "************************" << endl;
      ntp_prop.get_attribute("nb_step",process.temps->nb_step);
      process.temps->time_step.resize(process.temps->nb_step);
      process.temps->nbpastemps= 0;
      for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
          XmlNode n_step =ntp.get_element("step",i_step);
          n_step.get_attribute("dt",process.temps->time_step[i_step].dt);
          n_step.get_attribute("t_ini",process.temps->time_step[i_step].t_ini);
          n_step.get_attribute("t_fin",process.temps->time_step[i_step].t_fin);
          n_step.get_attribute("nb_time_step",process.temps->time_step[i_step].nb_time_step);
          process.temps->nbpastemps+=process.temps->time_step[i_step].nb_time_step;
      }
      
//       ntp.get_attribute("nbpastemps",process.temps->nbpastemps);
//       ntp.get_attribute("pasdetemps",process.temps->dt);
      nm.get_attribute("save_depl_SST",process.latin->save_depl_SST);      
      nm.get_attribute("nom_calcul",process.nom_calcul);
      //ntp.get_attribute("theta",process.temps->theta);
   } else {
      std::cout << "Type de calcul non defini " << std::endl;
      assert(0);
   }
   
cout << "sortie " << endl;
   system(("mkdir -p "+process.affichage->repertoire_save).c_str());


};*/

