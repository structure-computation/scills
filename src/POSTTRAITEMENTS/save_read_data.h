#include "../ITERATIONS/manipulate_quantities.h"


/** \ingroup  PostTraitements
\brief Sauvegarde des quantites desirees pour les interfaces
 
L'utilisateur choisi les quantites qu'il souhaite sauvegarder dans un fichier nomme save_inter. On boucle sur les interfaces et  on indique tout d'abord le nom puis on sauvegarde pour chaque cote de l'interface selectionnee tous les pas de temps.
*/
template<class TV2, class TV1>
void save_data_inter(TV2 &Inter,TV1 &S, Process &process, Vec<Sc2String> &fields_to_save) {
    typedef TYPEREEL T;
    
    //nom du fichier de sauvegarde
    system(("mkdir -p "+process.affichage->repertoire_save).c_str());

    ifstream is;
    std::ostringstream ss1;
    ss1<<process.affichage->repertoire_save<<"save_inter_"<<process.parallelisation->rank;
    Sc2String name(ss1.str());
    ofstream os( name.c_str() );

    os << process.parallelisation->size << std::endl;
    if (process.parallelisation->is_local_cpu())
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
                        //if (pt == 2) std::cout << process.parallelisation->rank << " inter " << q << " pas de temps" << pt << " " << norm_2(Inter[q].side[j].t[pt].Fchap) << " " << Inter[q].side[j].t[pt].Fchap.size() << std::endl;
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
\brief Lecture des quantites desirees pour les interfaces
 
L'utilisateur choisi les quantites qu'il souhaite sauvegarder dans un fichier nomme save_inter. On boucle sur les interfaces et  on indique tout d'abord le nom puis on sauvegarde pour chaque cote de l'interface selectionnee tous les pas de temps.
*/
template<class TV2>
void read_data_inter(TV2 &Inter, Process &process) {
    //typedef typename TV2::template SubType<0>::T T;
    typedef TYPEREEL T;
    
    
    //nom du fichier de sauvegarde
    ifstream is;
    std::ostringstream ss1;
    ss1<<process.affichage->repertoire_save<<"save_inter_"<<process.parallelisation->rank;
    Sc2String name(ss1.str());
    is.open( name.c_str() );


    unsigned toutlire=1;
    Sc2String str;
    if (is ) {
      getline(is,str);
      if ((not process.parallelisation->is_multi_cpu()) and process.parallelisation->size != atoi(str.c_str()) ) {//ca va pas aller
        std::cout << "Pour lire les donnees soit on le fait monoprocesseur soit on doit avoir le meme nombre de pro que le calcul initial" << std::endl;
        assert(0);
      } else if (process.parallelisation->size != atoi(str.c_str()))//le master doit lire tout les fichiers
        toutlire=atoi(str.c_str());
    } else {
        std::cout << "Aucun fichier de sauvegarde a lire" << std::endl;
        assert(0);
    }
    is.close();
    
    
    for( unsigned i=0;i<toutlire ;i++ ){
      std::ostringstream ss2;
      if (toutlire==1) ss2<<process.affichage->repertoire_save<<"save_inter_"<<process.parallelisation->rank;
      else ss2<<process.affichage->repertoire_save<<"save_inter_"<<i;
      Sc2String name2(ss2.str());
      is.open( name2.c_str() );
      getline(is,str);

    while(is) {
        //lecture numero de l'interface
        getline(is,str);

        if(str=="fin") {
            break;
        } else {
            istringstream s(str);
            Sc2String nom; //nom = inter
            s >> nom;
            unsigned q; //numero de l'interface
            s >> q ;

            while(true) {
                //lecture numero du cote
                Sc2String str2;
                getline(is,str2);
                //std::cout << str2 << std::endl;
                if(str2=="finSide") {
                    break;
                } else {
                    istringstream s2(str2);
                    Sc2String nom2; //nom = Side
                    s2 >> nom2;
                    unsigned j; //numero du cote
                    s2 >> j ;

                    while(true) {
                        //lecture du pas de temps
                        Sc2String str3;
                        getline(is,str3);
                        if(str3=="finPastemps") {
                            break;
                        } else {
                            istringstream s3(str3);
                            Sc2String nom3; //nom = Pastemps
                            s3 >> nom3;
                            unsigned pt; //numero du pas de temps
                            s3 >> pt ;

                            while(true) {
                                //lecture nom du champ
                                Sc2String str4;
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
\brief Sauvegarde des quantites desirees pour les sous-structures
 
*/
template<class TV1>
void save_data_sst(TV1 &S, Process &process, Vec<Sc2String> &fields_to_save) {
    typedef TYPEREEL T;

    //nom du fichier de sauvegarde
    system(("mkdir -p "+process.affichage->repertoire_save).c_str());

    ifstream is;
    std::ostringstream ss1;
    ss1<<process.affichage->repertoire_save<<"save_sst_"<<process.parallelisation->rank;
    Sc2String name(ss1.str());
    ofstream os( name.c_str() );

    os << process.parallelisation->size << std::endl;
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
                    //std::cout << process.parallelisation->rank << " sst " << q << " pas de temps" << pt << " " << norm_2(S[q].t[pt].q) << " " << S[q].t[pt].q.size() << std::endl;
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
\brief Lecture des quantites desirees pour les interfaces
 */
template<class TV1>
void read_data_sst(TV1 &S, Process &process) {
    typedef TYPEREEL T;
    //nom du fichier de sauvegarde
    ifstream is;
    std::ostringstream ss1;
    ss1<<process.affichage->repertoire_save<<"save_sst_"<<process.parallelisation->rank;
    Sc2String name(ss1.str());
    is.open( name.c_str() );


    unsigned toutlire=1;
    Sc2String str;
    if (is) {
      getline(is,str);
      if (not(process.parallelisation->is_multi_cpu()) and process.parallelisation->size != atoi(str.c_str()) ) {//ca va pas aller
        std::cout << "Pour lire les donnees soit on le fait monoprocesseur soit on doit avoir le meme nombre de pro que le calcul initial" << std::endl;
        assert(0);
      } else if (process.parallelisation->size != atoi(str.c_str()))//le master doit lire tout les fichiers
        toutlire=atoi(str.c_str());
    } else {
      std::cout << "Aucun fichier de sauvegarde a lire" << std::endl;
      assert(0);
    }
    is.close();
    
    
    for( unsigned i=0;i<toutlire ;i++ ){
      std::ostringstream ss2;
      if (toutlire==1) ss2<<process.affichage->repertoire_save<<"save_sst_"<<process.parallelisation->rank;
      else ss2<<process.affichage->repertoire_save<<"save_sst_"<<i;
      Sc2String name2(ss2.str());
      is.open( name2.c_str() );
      getline(is,str);

      while(true) {
        //lecture numero de la sst
        Sc2String str;
        getline(is,str);

        if(str=="fin") {
          is.close();
          break;
        } else {
            istringstream s(str);
            Sc2String nom; //nom = Sst
            s >> nom;
            unsigned q; //numero de la sst
            s >> q ;

            while(true) {
                //lecture du pas de temps
                Sc2String str3;
                getline(is,str3);
                if(str3=="finPastemps") {
                    break;
                } else {
                    istringstream s3(str3);
                    Sc2String nom3; //nom = Pastemps
                    s3 >> nom3;
                    unsigned pt; //numero du pas de temps
                    s3 >> pt ;

                    while(true) {
                        //lecture nom du champ
                        Sc2String str4;
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
