
using namespace LMT;
using namespace std;

inline void parametre_inconnu(Vec<string> &v,Param &process) {
    if (process.rank==0)
        cout << "Commande Inconnue. Taper help pour avoir les commandes disponibles." << endl;
}

template <class T1,class T2>
int findalain(Mat<T1> &m, T2 &s) {
    for(unsigned i=0 ;i<m.nb_rows() ;i++ ) {
        if (m(i,0)==s)
            return i;
    }
    return -1;
}

template <class T2>
int findalain2(Vec<Vec<string> > &m, T2 &s) {
    for(unsigned i=0 ;i<m.size() ;i++ ) {
        if (m[i][0]==s)
            return i;
    }
    return -1;
}

inline void help(Vec<string> &v, Param &process) {
    //    Mat<string> m(10,3);
    //    m(0,0)="help";m(0,1)="Affiche les commandes disponibles";m(0,2)="help";
    //    m(1,0)="evol";m(1,1)="Trace l'evolution d'un point";m(1,2)="evol coor1 coor2 (coor3)";
    //    m(2,0)="trac_sst";m(2,1)="Trace la deformee + champs donnes";m(2,2)="trac_sst champ1 champ2 champ3 ... ou trac_sst ";
    //    m(3,0)="trac_inter";m(3,1)="Trace les quantites sur les interfaces ";m(3,2)="trac_inter cote pas_de_temps num1 num2 num3 ou -1 pour afficher toutes les interfaces";
    //    m(4,0)="mesh_sst";m(4,1)="Trace des maillages des ssts ";m(4,2)="mesh_sst pour afficher toutes les ssts";
    //    m(5,0)="mesh_inter";m(5,1)="Trace des maillages des interfaces ";m(5,2)="mesh_inter cote pour afficher toutes les interfaces";
    //    m(6,0)="trac_sst_temps";m(6,1)="Trace de la deformee + champs au cours des pas de temps";m(6,2)="trac_sst_temps champ1 champ2 champ3 ... ou trac_sst_temps ";
    //    m(7,0)="trac_inter_temps";m(7,1)="Trace de la deformee + champs au cours des pas de temps";m(7,2)="trac_inter_temps cote  num1 num2 num3 ou -1 pour afficher toutes les interfaces ou trac_inter_temps seul pour le cote 0 et toutes les interfaces";
    //    m(8,0)="trac_error";m(8,1)="Trace de l'erreur latin au cours des iterations ";m(8,2)="trac_error";
    //    m(9,0)="trac_ener";m(9,1)="Trace de la dissipation ou energie imposee au cours des pas de temps à partir des quantités chapeaux ou quantités n";m(9,2)="trac_ener (0)ener_dissip/(1)ener_impo (0)Qchap/(1)Qn";

    Vec<Vec<string> > m;
    Vec<string> help_cur;
    help_cur.resize(3);
    help_cur[0]="help";
    help_cur[1]="Affiche les commandes disponibles";
    help_cur[2]="help";
    m.push_back(help_cur);
    help_cur[0]="evol";
    help_cur[1]="Trace l'evolution d'un point";
    help_cur[2]="evol coor1 coor2 (coor3)";
    m.push_back(help_cur);
    help_cur[0]="evol_inter";
    help_cur[1]="Trace l'evolution des champs d'un élément d'une interface donnée";
    help_cur[2]="evol_inter num_inter num_elem (0)Fchap/(1)Wchap/(2)Wpchap/(3)All";
    m.push_back(help_cur);
    help_cur[0]="mesh_sst";
    help_cur[1]="Trace des maillages des ssts";
    help_cur[2]="mesh_sst pour afficher toutes les ssts";
    m.push_back(help_cur);
    help_cur[0]="mesh_inter";
    help_cur[1]="Trace des maillages des interfaces";
    help_cur[2]="mesh_inter cote pour afficher toutes les interfaces";
    m.push_back(help_cur);
    help_cur[0]="trac_sst_temps";
    help_cur[1]="Trace de la deformee + champs au cours des pas de temps";
    help_cur[2]="trac_sst_temps champ1 champ2 champ3 ... ou trac_sst_temps";
    m.push_back(help_cur);
    help_cur[0]="trac_inter_temps";
    help_cur[1]="Trace de la deformee + champs au cours des pas de temps";
    help_cur[2]="trac_inter_temps cote  num1 num2 num3 ou -1 pour afficher toutes les interfaces ou trac_inter_temps seul pour le cote 0 et toutes les interfaces";
    m.push_back(help_cur);
    help_cur[0]="trac_error";
    help_cur[1]="Trace de l'erreur latin au cours des iterations";
    help_cur[2]="trac_error";
    m.push_back(help_cur);
    help_cur[0]="trac_ener";
    help_cur[1]="Trace de la dissipation ou energie imposee au cours des pas de temps a partir des quantites chapeaux ou quantites n";
    help_cur[2]="trac_ener (0)ener_dissip/(1)ener_impo/(2)Ft2/(3)Fn/(4)Un/(5)Ut   (0)Qchap/(1)Qn";
    m.push_back(help_cur);
    help_cur[0]="trac_sst";
    help_cur[1]="Trace la deformee + champs donnes";
    help_cur[2]="trac_sst champ1 champ2 champ3 ... ou trac_sst";
    m.push_back(help_cur);
    help_cur[0]="trac_inter";
    help_cur[1]="Trace les quantites sur les interfaces";
    help_cur[2]="trac_inter cote pas_de_temps num1 num2 num3 ou -1 pour afficher toutes les interfaces";
    m.push_back(help_cur);
    help_cur[0]="modif_param";
    help_cur[1]="Modification des parametres de calcul";
    help_cur[2]="modif_param + (critere_erreur + valeur et/ou nbitermax + valeur)";
    m.push_back(help_cur);
    help_cur[0]="calcul";
    help_cur[1]="Lancement d'un nouveau calcul iteratif en latin (non implemente en MPI)";
    help_cur[2]="calcul +(latin ou incr (non implemente))";
    m.push_back(help_cur);

    /*   if (v.size() > 1) {
          if (findalain2(m,v[1]) > 0 )
             cout << m(findalain2(m,v[1]),2) << endl;
          else {
             cout << "Liste des commandes disponibles :" << endl;
             for(unsigned i=0 ;i<m.nb_rows() ;i++ )
                cout  << setw(15) << m(i,0) << "   "  << m(i,1) << endl;
          }
       } else {
          cout << "Liste des commandes disponibles :" << endl;
          for(unsigned i=0 ;i<m.nb_rows() ;i++ )
             cout  << setw(15) << m(i,0) << "   "  << m(i,1) << endl;
       }*/
    if (v.size() > 1) {
        if (findalain2(m,v[1]) > 0 )
            if (process.rank==0) {
                cout << m[findalain2(m,v[1])][2] << endl;
            } else {
                if (process.rank==0)
                    cout << "Liste des commandes disponibles :" << endl;
                for(unsigned i=0 ;i<m.size() ;i++ )
                    if (process.rank==0)
                        cout  << setw(18) << m[i][0] << "   "  << m[i][1] << endl;
            }
    } else {
        if (process.rank==0)
            cout << "Liste des commandes disponibles :" << endl;
        for(unsigned i=0 ;i<m.size() ;i++ )
            if (process.rank==0)
                cout  << setw(18) << m[i][0] << "   "  << m[i][1] << endl;
    }
}

/** \ingroup  Post_Traitements
 \brief trace de l'evolution d'un point
 */
inline void evol(Vec<string> &v,Param &process) {
    //modification de la coordonnee du point
    if (v.size()!=process.dim+1) {
        if (process.rank==0)
            cout << "Mauvais choix de point" << endl;
        process.affichage->coor_point.resize(process.dim,0.);
        if (process.rank==0)
            cout << "On choisit par défaut : " << process.affichage->coor_point << endl;
    } else {
        process.affichage->coor_point.resize(process.dim);
        for(unsigned i=1;i<v.size();i++)
            process.affichage->coor_point[i-1]=atof(v[i].c_str()) ;
    }
    process.affichage->affich_depl_pt=1;
}

/** \ingroup  Post_Traitements
\brief trace de l'evolution d'un point
 */
inline void evol_inter(Vec<string> &v,Param &process) {
//modification de la coordonnee du point
    if (v.size()!=4) {
        if (process.rank==0)
            cout << "Mauvaise utilisation : help evol_inter" << endl;
        if (process.rank==0)
            cout << "On choisit par défaut : l'interface 0, élement 0 et l'affichage de toutes les composantes." << process.affichage->coor_point << endl;
        process.affichage->param_ener.set(0);
        process.affichage->param_ener[2]=3;
    } else {
        for(unsigned i=1;i<v.size();i++)
            process.affichage->param_ener[i-1]=atoi(v[i].c_str()) ;
    }
}

 /** \ingroup  Post_Traitements
\brief Affichage des Résultats pour les sous-structures
*/
inline void trac_sst(Vec<string> &v,Param &process) {
    if (process.dim == 2 ) process.affichage->type_affichage="Sinterieur";
    else process.affichage->type_affichage="Sbord";
    if (v.size()>1) {
        Vec<string> choixchamps;
        for(unsigned i=1;i<v.size();i++)
            choixchamps.push_back(v[i]);
        process.affichage->display_fields=choixchamps;
    }
    process.affichage->affich_resultat=1;
    process.affichage->save="display";
}

/** \ingroup  Post_Traitements
\brief Affichage des Résultats pour les sous-structures pour tous les pas de temps
*/
inline void trac_sst_temps(Vec<string> &v,Param &process) {
    if (process.dim == 2 ) process.affichage->type_affichage="Sinterieur";
    else process.affichage->type_affichage="Sbord";
    if (v.size()>1) {
        Vec<string> choixchamps;
        for(unsigned i=1;i<v.size();i++)
            choixchamps.push_back(v[i]);
        process.affichage->display_fields=choixchamps;
    }
    process.affichage->affich_resultat=1;
    process.affichage->save="save";
}

/** \ingroup  Post_Traitements
\brief Affichage des Résultats pour les interfaces
*/
inline void trac_inter(Vec<string> &v,Param &process) {

    process.affichage->affich_inter_data=1;
    process.affichage->type_affichage="Inter";
    process.affichage->save="display";

    if(v.size()<4) {
        if (process.rank==0)
            cout << "Mauvais nombre d'arguments" << endl;
        if (process.rank==0)
            cout << "Par défaut on trace le cote 0, le pas de temps 1 et toutes les interfaces" << endl;
        process.affichage->num_inter_select=Vec<int>(-1);
        process.affichage->side=0;
        process.affichage->pt=1;
    } else {
        string side = v[1];
        std::istringstream i1( v[1] );
        i1 >> process.affichage->side;
        istringstream i2( v[2] );
        i2 >> process.affichage->pt;
        Vec<int> Num_Inter;
        for(unsigned i=3;i<v.size();i++) {
            istringstream ii( v[i] );
            unsigned num;
            ii >> num;
            Num_Inter.push_back(num);
        }
        process.affichage->num_inter_select=Num_Inter;
    }
}

/** \ingroup  Post_Traitements
\brief Affichage des Résultats pour les interfaces pour tous les pas de temps
*/
inline void trac_inter_temps(Vec<string> &v,Param &process) {
    process.affichage->type_affichage="Inter";
    process.affichage->affich_resultat=1;
    process.affichage->save="save";
    if(v.size()<3) {
        if (process.rank==0)
            cout << "Mauvais nombre d'arguments" << endl;
        if (process.rank==0)
            cout << "Par défaut on trace le cote 0 et toutes les interfaces" << endl;
        process.affichage->num_inter_select=Vec<int>(-1);
        process.affichage->side=0;
        process.affichage->pt=1;
    } else {
        string side = v[1];
        std::istringstream i1( v[1] );
        i1 >> process.affichage->side;
        Vec<int> Num_Inter;
        for(unsigned i=2;i<v.size();i++) {
            istringstream ii( v[i] );
            unsigned num;
            ii >> num;
            Num_Inter.push_back(num);
        }
        process.affichage->num_inter_select=Num_Inter;
    }
}

/** \ingroup  Post_Traitements
\brief Affichage des maillages des sous-structures
*/
inline void trac_mesh_sst(Vec<string> &v,Param &process) {
    if (v.size()>1) {
        Vec<string> choixchamps;
        for(unsigned i=1;i<v.size();i++)
            choixchamps.push_back(v[i]);
        process.affichage->display_fields=choixchamps;
    }
    process.affichage->type_affichage="Sbord";
    process.affichage->affich_mesh=1;
    process.affichage->save="display";
}

/** \ingroup  Post_Traitements
\brief Affichage des maillages des interfaces
*/
inline void trac_mesh_inter(Vec<string> &v,Param &process) {
    process.affichage->affich_mesh=1;
    process.affichage->type_affichage="Inter";
    process.affichage->save="display";

    if(v.size()<2) {
        if (process.rank==0)
            cout << "Mauvais nombre d'arguments" << endl;
        if (process.rank==0)
            cout << "Par défaut on trace le cote 0" << endl;
        process.affichage->side=0;
    } else {
        istringstream i1(v[1]);
        i1 >> process.affichage->side;
    }
}

/** \ingroup  Post_Traitements
\brief Affichage de l'erreur latin
*/
inline void trac_error(Vec<string> &v,Param &process) {
    process.affichage->display_error=1;
}

/** \ingroup  Post_Traitements
 \brief Trace de l'energie de dissipation ou imposée pour les quantités chapeua ou quantités n au cours du temps
 */
inline void trac_ener(Vec<string> &v,Param &process) {
    if (v.size()!=3) {
        if (process.rank==0)
            cout << "Parametres par defaut : dissip - Qchap" << endl;
        process.affichage->param_ener.set(0);
    } else {
        Vec<int,2> num;
        for(unsigned i=1;i<3;i++) {
            istringstream ii( v[i] );
            ii >> num[i-1];
        }
        if(num[0]==0 or num[0]==1 or num[0]==2 or num[0]==3 or num[0]==4 or num[0]==5 or num[0]==6) {
            process.affichage->param_ener[0]=num[0];
        } else {
            if (process.rank==0)
                cout << "Deuxieme argument non pris en compte -> par defaut ener_dissip" << endl;
            process.affichage->param_ener[0]=0;
        }

        if(num[1]==0 or num[1]==1) {
            process.affichage->param_ener[1]=num[1];
        } else {
            if (process.rank==0)
                cout << "Troisieme argument non pris en compte -> par defaut Qchap" << endl;
            process.affichage->param_ener[1]=0;
        }
    }
}

inline void trac_cl(Vec<string> &v,Param &process) {
    //rien a faire
}

template <class TV1, class TI>
void affichage_cl(TV1 &S, TI &Inter,Param &process) {
    for(unsigned i=0;i<S.size();i++)
        for( unsigned j=0;j<S[i].edge.size() ;j++ ) {
            unsigned internum=S[i].edge[j].internum;
            if (Inter[internum].comp=="effort") {
                cout << "Interface " << internum << " de type " << Inter[internum].comp << endl;
                for( unsigned pt=0;pt<Inter[internum].side[0].t.size() ;pt++ ) {
                    cout << "Pas de temps " << setw(3) << pt << " : " ;
                    if (process.dim == 2) {
                        cout << setw(14) << mean(Inter[internum].side[0].t[pt].F[range(0,2,(int)Inter[internum].side[0].t[pt].F.size())]) << " ";
                        cout << setw(14) << mean(Inter[internum].side[0].t[pt].F[range(1,2,(int)Inter[internum].side[0].t[pt].F.size())]) << endl;
                    }
                    if (process.dim == 3) {
                        cout << setw(14) << mean(Inter[internum].side[0].t[pt].F[range(0,3,(int)Inter[internum].side[0].t[pt].F.size())]) << " ";
                        cout << setw(14) << mean(Inter[internum].side[0].t[pt].F[range(1,3,(int)Inter[internum].side[0].t[pt].F.size())]) << " ";
                        cout << setw(14) << mean(Inter[internum].side[0].t[pt].F[range(2,3,(int)Inter[internum].side[0].t[pt].F.size())]) << endl;
                    }
                }
            }
            if (Inter[internum].comp=="depl" or Inter[internum].comp=="depl_normal") {
                cout << "Interface " << internum << " de type " << Inter[internum].comp << endl;
                for( unsigned pt=0;pt<Inter[internum].side[0].t.size() ;pt++ ) {
                    cout << "Pas de temps " << setw(3) << pt << " : " ;
                    if (process.dim == 2) {
                        cout << setw(14) << mean(Inter[internum].side[0].t[pt].Wp[range(0,2,(int)Inter[internum].side[0].t[pt].Wp.size())]) << " ";
                        cout << setw(14) << mean(Inter[internum].side[0].t[pt].Wp[range(1,2,(int)Inter[internum].side[0].t[pt].Wp.size())]) << endl;
                    }
                    if (process.dim == 3) {
                        cout << setw(14) << mean(Inter[internum].side[0].t[pt].Wp[range(0,3,(int)Inter[internum].side[0].t[pt].Wp.size())]) << " ";
                        cout << setw(14) << mean(Inter[internum].side[0].t[pt].Wp[range(1,3,(int)Inter[internum].side[0].t[pt].Wp.size())]) << " ";
                        cout << setw(14) << mean(Inter[internum].side[0].t[pt].Wp[range(2,3,(int)Inter[internum].side[0].t[pt].Wp.size())]) << endl;
                    }
                }
            }
            if (Inter[internum].comp=="Jeu_impose") {
                cout << "Interface " << internum << " de type " << Inter[internum].comp << endl;
                for( unsigned pt=0;pt<Inter[internum].side[0].t.size() ;pt++ ) {
                    cout << "Pas de temps " << setw(3) << pt << " : " ;
                    if (process.dim == 2) {
                      cout << setw(14) << mean(Inter[internum].side[1].t[pt].Wp[range(0,2,(int)Inter[internum].side[0].t[pt].Wp.size())]-Inter[internum].side[0].t[pt].Wp[range(0,2,(int)Inter[internum].side[0].t[pt].Wp.size())]) << " ";
                      cout << setw(14) << mean(Inter[internum].side[1].t[pt].Wp[range(1,2,(int)Inter[internum].side[0].t[pt].Wp.size())]-Inter[internum].side[0].t[pt].Wp[range(1,2,(int)Inter[internum].side[0].t[pt].Wp.size())]) << endl;
                    }
                    if (process.dim == 3) {
                      cout << setw(14) << mean(Inter[internum].side[1].t[pt].Wp[range(0,3,(int)Inter[internum].side[0].t[pt].Wp.size())]-Inter[internum].side[0].t[pt].Wp[range(0,3,(int)Inter[internum].side[0].t[pt].Wp.size())]) << " ";
                      cout << setw(14) << mean(Inter[internum].side[1].t[pt].Wp[range(1,3,(int)Inter[internum].side[0].t[pt].Wp.size())]-Inter[internum].side[0].t[pt].Wp[range(1,3,(int)Inter[internum].side[0].t[pt].Wp.size())]) << " ";
                      cout << setw(14) << mean(Inter[internum].side[1].t[pt].Wp[range(2,3,(int)Inter[internum].side[0].t[pt].Wp.size())]-Inter[internum].side[0].t[pt].Wp[range(2,3,(int)Inter[internum].side[0].t[pt].Wp.size())]) << endl;
                    }
                }
            }
        }
}



/** \ingroup  Post_Traitements
 \brief Modification des parametres de calcul
 
 Possibilite de modifier l'erreur, le nombre d'iterations
 */
inline void modif_param(Vec<string> &v,Param &process) {

    if (v.size()==1) {
        if (process.rank==0)
            cout << "On garde les parametres actuels : " << endl;
        if (process.rank==0)
            cout << "Nombre d'iterations max : " << process.latin->nbitermax << endl;
        if (process.rank==0)
            cout << "Critere d'erreur : " << process.latin->critere_erreur << endl;
    } else if(v.size()==3) {
        if(v[1]=="nbitermax") {
            if (process.rank==0)
                cout << "Modification des parametres " << endl;
            if (process.rank==0)
                cout << "Nombre d'iterations max actuel :" <<  process.latin->nbitermax<< endl;
            process.latin->nbitermax=atoi(v[2].c_str());
            if (process.rank==0)
                cout << "Nouveau Nombre d'iterations max :" << process.latin->nbitermax << endl;
            process.latin->error.resize(process.latin->nbitermax+1,0.);
        } else if (v[1]=="critere_erreur") {
            if (process.rank==0)
                cout << "Modification des parametres " << endl;
            if (process.rank==0)
                cout << "Critere d'erreur actuel : " <<  process.latin->critere_erreur<< endl;
            process.latin->critere_erreur=atof(v[2].c_str());
            if (process.rank==0)
                cout << "Nouveau Critere d'erreur :" << process.latin->critere_erreur << endl;
        } else {
            if (process.rank==0)
                cout << "Mauvais choix de modification a faire" << endl;
        }
    } else if(v.size()==5) {
        if(v[1]=="nbitermax") {
            if (process.rank==0)
                cout << "Modification des parametres " << endl;
            if (process.rank==0)
                cout << "Nombre d'iterations max actuel :" <<  process.latin->nbitermax<< endl;
            process.latin->nbitermax=atoi(v[2].c_str());
            if (process.rank==0)
                cout << "Nouveau Nombre d'iterations max :" << process.latin->nbitermax << endl;
            process.latin->error.resize(process.latin->nbitermax+1,0.);
        } else if (v[1]=="critere_erreur") {
            if (process.rank==0)
                cout << "Modification des parametres " << endl;
            if (process.rank==0)
                cout << "Critere d'erreur actuel : " <<  process.latin->critere_erreur<< endl;
            process.latin->critere_erreur=atof(v[2].c_str());
            if (process.rank==0)
                cout << "Nouveau Critere d'erreur :" << process.latin->critere_erreur << endl;
        } else {
            if (process.rank==0)
                cout << "Mauvais choix de modification a faire" << endl;
        }
        if(v[3]=="nbitermax") {
            if (process.rank==0)
                cout << "Modification des parametres " << endl;
            if (process.rank==0)
                cout << "Nombre d'iterations max actuel :" <<  process.latin->nbitermax<< endl;
            process.latin->nbitermax=atoi(v[4].c_str());
            if (process.rank==0)
                cout << "Nouveau Nombre d'iterations max :" << process.latin->nbitermax << endl;
            process.latin->error.resize(process.latin->nbitermax+1,0.);
        } else if (v[3]=="critere_erreur") {
            if (process.rank==0)
                cout << "Modification des parametres " << endl;
            if (process.rank==0)
                cout << "Critere d'erreur actuel : " <<  process.latin->critere_erreur<< endl;
            process.latin->critere_erreur=atof(v[4].c_str());
            if (process.rank==0)
                cout << "Nouveau Critere d'erreur :" << process.latin->critere_erreur << endl;
        } else {
            if (process.rank==0)
                cout << "Mauvais choix de modification a faire" << endl;
        }
    } else {
        if (process.rank==0)
            cout << "Mauvais choix de modification a faire" << endl;
    }
}

/** \ingroup  Post_Traitements
 \brief Modification des parametres de calcul
 
 Possibilite de modifier l'erreur, le nombre d'iterations
 */
inline void modif_param_inter(Vec<string> &v,Param &process) {
    process.reprise_calcul=2;
}

/** \ingroup  Post_Traitements
 \brief Lancement d'un nouveau calcul multiechelle
 
 Possibilite de modifier l'erreur, le nombre d'iterations
 */
inline void calcul(Vec<string> &v,Param &process) {
    /*   if(v.size()!=2){cout << "Nouveau calcul en latin" << endl;}
       else{
          if(v[1]=="latin"){process.nom_calcul="latin";}
          else if(v[1]=="incr"){cout << "Non implemente pour un calcul incremental : on lance un calcul latin" << endl; process.nom_calcul="latin";}
          else{cout << "Type de calcul non implemente : on lance un calcul latin" << endl; process.nom_calcul="latin";}
       }*/

    //process.latin->alloc_quantites=0;
}
