
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

    os << process.size << endl;
    if (process.size==1 or process.rank>0)
      for( unsigned i=0;i<S.size() ;i++ ){
      for( unsigned kk=0;kk<S[i].edge.size() ;kk++ ){
        unsigned q=S[i].edge[kk].internum;
        unsigned j=S[i].edge[kk].datanum;
//    for(unsigned q=0;q<Inter.size();q++) {
        os<< "Inter " << Inter[q].num <<endl;
//        for(unsigned j=0;j<Inter[q].side.size();j++) {
            os << "Side " << j << endl;
            for(unsigned pt=1;pt<=process.temps->nbpastemps;pt++) {
                os << "Pastemps " << pt<<endl;
                for(unsigned k=0;k<fields_to_save.size();k++) {
                    if(fields_to_save[k]=="Fchap") {
                        os << "Fchap" <<endl;
                        if(process.nom_calcul=="latin") os.write( (char *)Inter[q].side[j].t[pt].Fchap.ptr(), sizeof(T)*Inter[q].side[j].t[pt].Fchap.size() );
                        if(process.nom_calcul=="incr") os.write( (char *)Inter[q].side[j].t_post[pt].Fchap.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].Fchap.size() );
                        //if (pt == 2) cout << process.rank << " inter " << q << " pas de temps" << pt << " " << norm_2(Inter[q].side[j].t[pt].Fchap) << " " << Inter[q].side[j].t[pt].Fchap.size() << endl;
                    } else if(fields_to_save[k]=="F") {
                        os << "F" <<endl;
                        if(process.nom_calcul=="latin") os.write( (char *)Inter[q].side[j].t[pt].F.ptr(), sizeof(T)*Inter[q].side[j].t[pt].F.size() );
                        if(process.nom_calcul=="incr") os.write( (char *)Inter[q].side[j].t_post[pt].F.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].F.size() );
                    } else if(fields_to_save[k]=="Wchap") {
                        os << "Wchap" <<endl;
                        if(process.nom_calcul=="latin") os.write( (char *)Inter[q].side[j].t[pt].Wchap.ptr(), sizeof(T)*Inter[q].side[j].t[pt].Wchap.size() );
                        if(process.nom_calcul=="incr") os.write( (char *)Inter[q].side[j].t_post[pt].Wchap.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].Wchap.size() );
                    } else if(fields_to_save[k]=="Wpchap") {
                        os << "Wpchap" <<endl;
                        if(process.nom_calcul=="latin") os.write( (char *)Inter[q].side[j].t[pt].Wpchap.ptr(), sizeof(T)*Inter[q].side[j].t[pt].Wpchap.size() );
                        if(process.nom_calcul=="incr") os.write( (char *)Inter[q].side[j].t_post[pt].Wpchap.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].Wpchap.size() );
                    } else if(fields_to_save[k]=="W") {
                        os << "W" <<endl;
                        if(process.nom_calcul=="latin") os.write( (char *)Inter[q].side[j].t[pt].W.ptr(), sizeof(T)*Inter[q].side[j].t[pt].W.size() );
                        if(process.nom_calcul=="incr") os.write( (char *)Inter[q].side[j].t_post[pt].W.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].W.size() );
                    } else if(fields_to_save[k]=="Wp") {
                        os << "Wp" <<endl;
                        if(process.nom_calcul=="latin") os.write( (char *)Inter[q].side[j].t[pt].Wp.ptr(), sizeof(T)*Inter[q].side[j].t[pt].Wp.size() );
                        if(process.nom_calcul=="incr") os.write( (char *)Inter[q].side[j].t_post[pt].Wp.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].Wp.size() );
                    } else if(fields_to_save[k]=="WtildeM") {
                        os << "WtildeM" <<endl;
                        if(process.nom_calcul=="latin") os.write( (char *)Inter[q].side[j].t[pt].WtildeM.ptr(), sizeof(T)*Inter[q].side[j].t[pt].WtildeM.size() );
                        if(process.nom_calcul=="incr") os.write( (char *)Inter[q].side[j].t_post[pt].WtildeM.ptr(), sizeof(T)*Inter[q].side[j].t_post[pt].WtildeM.size() );
                    } else {
                        cout << "Erreur dans le choix des champs a sauvegarder " << endl;
                        assert(0);
                    }
                }
                os << "finchamps"<<endl;
            }
            os << "finPastemps" <<endl;
         os << "finSide" <<endl;
         }
    }
    os << "fin"<<endl;
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
        cout << "Pour lire les donnees soit on le fait monoprocesseur soit on doit avoir le meme nombre de pro que le calcul initial" << endl;
        assert(0);
      } else if (process.size != atoi(str.c_str()))//le master doit lire tout les fichiers
        toutlire=atoi(str.c_str());
    } else {
        cout << "Aucun fichier de sauvegarde a lire" << endl;
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
                //cout << str2 << endl;
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
                                  //cout << Inter[q].side[j].t_post[pt].F << endl;
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
                                    cout << "Erreur dans le choix des champs a lire " << endl;
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

    os << process.size << endl;
    if (process.latin->save_depl_SST ==1)
    for(unsigned q=0;q<S.size();q++) {
        os<< "S " << S[q].num <<endl;
        for(unsigned pt=1;pt<=process.temps->nbpastemps;pt++) {
            os << "Pastemps " << pt <<endl;
            for(unsigned k=0;k<fields_to_save.size();k++) {
                if(fields_to_save[k]=="q") {
                    os << "q" << endl;
                    if(process.nom_calcul=="latin") os.write( (char *)S[q].t[pt].q.ptr(), sizeof(T)*S[q].t[pt].q.size() );
                    if(process.nom_calcul=="incr")  os.write( (char *)S[q].t_post[pt].q.ptr(), sizeof(T)*S[q].t_post[pt].q.size() );
                    //cout << process.rank << " sst " << q << " pas de temps" << pt << " " << norm_2(S[q].t[pt].q) << " " << S[q].t[pt].q.size() << endl;
                } else {
                    cout << "Erreur dans le choix des champs a sauvegarder " << endl;
                    assert(0);
                }
            }
            os << "finchamps"<<endl;
        }
        os << "finPastemps" <<endl;
    }
    os << "fin"<<endl;
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
        cout << "Pour lire les donnees soit on le fait monoprocesseur soit on doit avoir le meme nombre de pro que le calcul initial" << endl;
        assert(0);
      } else if (process.size != atoi(str.c_str()))//le master doit lire tout les fichiers
        toutlire=atoi(str.c_str());
    } else {
      cout << "Aucun fichier de sauvegarde a lire" << endl;
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
                        //cout << str4 << endl;
                        if(str4=="q") {
                          if(process.nom_calcul=="latin") is.read( (char *)S[q].t[pt].q.ptr(), sizeof(T)*S[q].t[pt].q.size() );
                          if(process.nom_calcul=="incr") is.read( (char *)S[q].t_post[pt].q.ptr(), sizeof(T)*S[q].t_post[pt].q.size() );
                         } else if(str4=="finchamps") {
                            break;
                        } else {
                            cout << "Erreur dans le choix des champs a lire " << endl;
                            assert(0);
                        }
                    }
                }
            }
        }
    }
    }
}
