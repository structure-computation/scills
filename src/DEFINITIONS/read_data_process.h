using namespace LMT;
using namespace std;

// lecture des parametres de calcul
inline void read_data_process(Param &process, const XmlNode &n, DataUser &data_user) {
//    XmlNode nm = n.get_element("parametres");
//    nm.get_attribute("sous_integration",process.sousint);
//    nm.get_attribute("type_sous_integration",process.type_sousint);
//    nm.get_attribute("blocage_modes_rigides",process.rbm.bloq);
//    nm.get_attribute("mvts_bloques",process.rbm.mvts_bloques);   
//    nm.get_attribute("nb_threads",process.nb_threads);
//    nm.get_attribute("save_data",process.save_data);
//    nm.get_attribute("read_data",process.read_data);
//    nm.get_attribute("reprise_calcul",process.reprise_calcul);
//    nm.get_attribute("deltaT",process.properties->deltaT);
//       
//    nm.get_attribute("type_base_macro",process.multiscale->type_base_macro);
//    nm.get_attribute("multiechelle",process.multiscale->multiechelle);
//    nm.get_attribute("opti_multi",process.multiscale->opti_multi);
//    nm.get_attribute("erreur_macro",process.multiscale->erreur_macro);

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
   
   process.latin->nbitermax
   process.latin->facteur_relaxation
   process.latin->critere_erreur
   process.latin->type_error
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
   ntp.get_attribute("type_de_calcul",process.temps->type_de_calcul);

   if (process.rank==0) cout << "************************" << endl;
   if (process.temps->type_de_calcul=="stat") {
      if (process.rank==0) cout << "     STATIQUE     " << endl;
      if (process.rank==0) cout << "************************" << endl;
      if (process.rank==0) cout << " Rq : 1 seul pas de temps automatiquement, dt=1 par defaut " << endl;
      process.temps->nbpastemps=1;
      process.temps->dt=1;
      process.nom_calcul="incr";
      if (process.rank==0) cout << " Rq : Attention la valeur de la fonction spatiale sera tout de meme modulee par la fonction temporelle" << endl;
   } else if(process.temps->type_de_calcul=="Qstat") {
      if (process.rank==0) cout << "     QUASISTATIQUE          " << endl;
      if (process.rank==0) cout << "************************" << endl;
      ntp.get_attribute("nbpastemps",process.temps->nbpastemps);
      ntp.get_attribute("pasdetemps",process.temps->dt);
      nm.get_attribute("save_depl_SST",process.latin->save_depl_SST);      
      nm.get_attribute("nom_calcul",process.nom_calcul);
      //ntp.get_attribute("theta",process.temps->theta);
   } else {
      cout << "Type de calcul non defini " << endl;
      assert(0);
   }
   

   system(("mkdir -p "+process.affichage->repertoire_save).c_str());


};

